!  SVN:$Id: ice_dyn_evp.F90 1228 2017-05-23 21:33:34Z tcraig $
!=======================================================================
!
! Elastic-viscous-plastic sea ice dynamics model
! Computes ice velocity and deformation
!
! See:
!
! Hunke, E. C., and J. K. Dukowicz (1997). An elastic-viscous-plastic model
! for sea ice dynamics. {\em J. Phys. Oceanogr.}, {\bf 27}, 1849--1867.
!
! Hunke, E. C. (2001).  Viscous-Plastic Sea Ice Dynamics with the EVP Model:
! Linearization Issues. {\em Journal of Computational Physics}, {\bf 170},
! 18--38.
!
! Hunke, E. C., and J. K. Dukowicz (2002).  The Elastic-Viscous-Plastic
! Sea Ice Dynamics Model in General Orthogonal Curvilinear Coordinates
! on a Sphere---Incorporation of Metric Terms. {\em Monthly Weather Review},
! {\bf 130}, 1848--1865.
!
! Hunke, E. C., and J. K. Dukowicz (2003).  The sea ice momentum
! equation in the free drift regime.  Los Alamos Tech. Rep. LA-UR-03-2219.
!
! Bouillon, S., T. Fichefet, V. Legat and G. Madec (submitted 2013).  The 
! revised elastic-viscous-plastic method.  Ocean Modelling.
!
! author: Elizabeth C. Hunke, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb (LANL)
! 2004: Block structure added by William Lipscomb
! 2005: Removed boundary calls for stress arrays (WHL)
! 2006: Streamlined for efficiency by Elizabeth Hunke
!       Converted to free source form (F90)

      module ice_dyn_vp

      use ice_kinds_mod
      use ice_blocks, only: nx_block, ny_block
      use ice_constants, only: field_loc_center, field_loc_NEcorner, &
          field_type_scalar, field_type_vector
      use ice_constants, only: c0, p027, p055, p111, p166, &
          p222, p25, p333, p5, c1
      use ice_domain, only: nblocks, distrb_info
      use ice_domain_size, only: max_blocks
      use ice_dyn_shared, only: dyn_prep1, dyn_prep2, dyn_finish, &
          ecci, cosw, sinw, fcor_blk, uvel_init,  &
          vvel_init, basal_stress_coeff, basalstress, Ktens
      use ice_fileunits, only: nu_diag
      use ice_flux, only: fm
      use ice_global_reductions, only: global_sum, global_sums
      use ice_grid, only: dxt, dyt, dxhy, dyhx, cxp, cyp, cxm, cym, uarear
      use ice_exit, only: abort_ice
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_ice_strength, icepack_query_parameters

      implicit none
      private
      public :: imp_solver, init_vp

      ! namelist parameters

      integer (kind=int_kind), public :: &
         maxits_nonlin  , & ! max nb of iteration for nonlinear solver
         precond        , & ! preconditioner for fgmres: 1: identity, 2: diagonal 3: pgmres + diag
         im_fgmres      , & ! size of fgmres Krylov subspace
         im_pgmres      , & ! size of pgmres Krylov subspace
         maxits_fgmres  , & ! max nb of iteration for fgmres
         maxits_pgmres  , & ! max nb of iteration for pgmres
         monitor_fgmres , & ! print fgmres residual norm
         monitor_pgmres , & ! print pgmres residual norm
         algo_nonlin    , & ! nonlinear algorithm: 1: Picard iteration, 2: Anderson acceleration (andacc)
         fpfunc_andacc  , & ! fixed point function for Anderson acceleration: 1: g(x) = FMGRES(A(x),b(x)), 2: g(x) = x - A(x)x + b(x)
         im_andacc      , & ! size of Anderson minimization matrix (number of saved previous residuals)
         start_andacc       ! acceleration delay factor (acceleration starts at this iteration)

      logical (kind=log_kind), public :: &
         monitor_nonlin , & ! print nonlinear residual norm
         use_mean_vrel      ! use mean of previous 2 iterates to compute vrel

      real (kind=dbl_kind), public :: &
         gammaNL        , & ! nonlinear stopping criterion: gammaNL*res(k=0)
         gamma          , & ! fgmres stopping criterion: gamma*res(k)
         epsprecond     , & ! pgmres stopping criterion: epsprecond*res(k)
         damping_andacc , & ! damping factor for Anderson acceleration
         reltol_andacc      ! relative tolerance for Anderson acceleration

      ! mmodule variables

      integer (kind=int_kind), allocatable :: &
         icellt(:)    , & ! no. of cells where icetmask = 1
         icellu(:)        ! no. of cells where iceumask = 1

      integer (kind=int_kind), allocatable :: &
         indxti(:,:)  , & ! compressed index in i-direction
         indxtj(:,:)  , & ! compressed index in j-direction
         indxui(:,:)  , & ! compressed index in i-direction
         indxuj(:,:)      ! compressed index in j-direction

!=======================================================================

      contains

!=======================================================================

! Initialize parameters and variables needed for the vp dynamics
! author: Philippe Blain, ECCC

      subroutine init_vp (dt)
      
      use ice_blocks, only: get_block, block
      use ice_boundary, only: ice_HaloUpdate
      use ice_constants, only: c1, &
          field_loc_center, field_type_scalar
      use ice_domain, only: blocks_ice, halo_info
      use ice_dyn_shared, only: init_evp
      use ice_grid, only: tarea, tinyarea
      
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step
      
      ! local variables
      
      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      type (block) :: &
         this_block           ! block information for current block
         
      real (kind=dbl_kind) :: &
         puny_vp = 2e-09_dbl_kind      ! special puny value for computing tinyarea
      
      ! Initialize variables shared with evp
      call init_evp(dt)
      
      ! Initialize module variables
      allocate(icellt(max_blocks), icellu(max_blocks))
      allocate(indxti(nx_block*ny_block, max_blocks), &
               indxtj(nx_block*ny_block, max_blocks), &
               indxui(nx_block*ny_block, max_blocks), &
               indxuj(nx_block*ny_block, max_blocks))
      
      ! Redefine tinyarea using a different puny value
      
      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            tinyarea(i,j,iblk) = puny_vp*tarea(i,j,iblk)
         enddo
         enddo
      enddo                     ! iblk
      !$OMP END PARALLEL DO
      
      call ice_HaloUpdate (tinyarea,           halo_info, &
                           field_loc_center,   field_type_scalar, &
                           fillValue=c1)
      
      end subroutine init_vp
!=======================================================================

! Viscous-plastic dynamics driver
!
#ifdef CICE_IN_NEMO
! Wind stress is set during this routine from the values supplied
! via NEMO (unless calc_strair is true).  These values are supplied 
! rotated on u grid and multiplied by aice.  strairxT = 0 in this 
! case so operations in dyn_prep1 are pointless but carried out to 
! minimise code changes.
#endif
!
! author: JF Lemieux, A. Qaddouri and F. Dupont ECCC

      subroutine imp_solver (dt)

      use ice_arrays_column, only: Cdn_ocn
      use ice_boundary, only: ice_halo, ice_HaloMask, ice_HaloUpdate, &
          ice_HaloDestroy, ice_HaloUpdate_stress
      use ice_blocks, only: block, get_block, nx_block, ny_block
      use ice_domain, only: blocks_ice, halo_info, maskhalo_dyn
      use ice_domain_size, only: max_blocks, ncat
      use ice_flux, only: rdg_conv, rdg_shear, strairxT, strairyT, &
          strairx, strairy, uocn, vocn, ss_tltx, ss_tlty, iceumask, fm, &
          strtltx, strtlty, strocnx, strocny, strintx, strinty, taubx, tauby, &
          strocnxT, strocnyT, strax, stray, &
          Tbu, hwater, &
          stressp_1, stressp_2, stressp_3, stressp_4, &
          stressm_1, stressm_2, stressm_3, stressm_4, &
          stress12_1, stress12_2, stress12_3, stress12_4
      use ice_grid, only: tmask, umask, dxt, dyt, dxhy, dyhx, cxp, cyp, cxm, cym, &
          tarear, to_ugrid, t2ugrid_vector, u2tgrid_vector, &
          grid_type
      use ice_state, only: aice, vice, vsno, uvel, vvel, divu, shear, &
          aice_init, aice0, aicen, vicen, strength
      use ice_timers, only: timer_dynamics, timer_bound, &
          ice_timer_start, ice_timer_stop

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) :: &
         ntot           , & ! size of problem for fgmres (for given cpu)
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, ij

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         tmass    , & ! total mass of ice and snow (kg/m^2)
         waterx   , & ! for ocean stress calculation, x (m/s)
         watery   , & ! for ocean stress calculation, y (m/s)
         forcex   , & ! work array: combined atm stress and ocn tilt, x
         forcey   , & ! work array: combined atm stress and ocn tilt, y
         bxfix    , & ! part of bx that is constant during Picard 
         byfix    , & ! part of by that is constant during Picard
         fpresx   , & ! x fixed point residual vector, fx = uvel - uprev_k
         fpresy   , & ! y fixed point residual vector, fy = vvel - vprev_k
         aiu      , & ! ice fraction on u-grid
         umass    , & ! total mass of ice and snow (u grid)
         umassdti     ! mass of U-cell/dte (kg/m^2 s)
         
      real (kind=dbl_kind), allocatable :: fld2(:,:,:,:)

      logical (kind=log_kind) :: calc_strair

      integer (kind=int_kind), dimension (nx_block,ny_block,max_blocks) :: &
         icetmask, &  ! ice extent mask (T-cell)
         halomask     ! generic halo mask

      type (ice_halo) :: &
         halo_info_mask !  ghost cell update info for masked halo

      type (block) :: &
         this_block           ! block information for current block
      
      real (kind=dbl_kind), allocatable :: &
         bvec(:)     , & ! right-hand-side vector
         sol(:)      , & ! solution vector
         diagvec(:)      ! diagonal vector
      
      character(len=*), parameter :: subname = '(imp_solver)'
      
      call ice_timer_start(timer_dynamics) ! dynamics

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      allocate(fld2(nx_block,ny_block,2,max_blocks))
      
       ! This call is needed only if dt changes during runtime.
!      call set_evp_parameters (dt)

      !-----------------------------------------------------------------
      ! boundary updates
      ! commented out because the ghost cells are freshly 
      ! updated after cleanup_itd
      !-----------------------------------------------------------------

!      call ice_timer_start(timer_bound)
!      call ice_HaloUpdate (aice,              halo_info, &
!                           field_loc_center,  field_type_scalar)
!      call ice_HaloUpdate (vice,              halo_info, &
!                           field_loc_center,  field_type_scalar)
!      call ice_HaloUpdate (vsno,              halo_info, &
!                           field_loc_center,  field_type_scalar)
!      call ice_timer_stop(timer_bound)

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks

         do j = 1, ny_block 
         do i = 1, nx_block 
            rdg_conv (i,j,iblk) = c0 
            rdg_shear(i,j,iblk) = c0 
            divu (i,j,iblk) = c0 
            shear(i,j,iblk) = c0 
         enddo
         enddo

      !-----------------------------------------------------------------
      ! preparation for dynamics 
      !-----------------------------------------------------------------

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         call dyn_prep1 (nx_block,           ny_block,           & 
                         ilo, ihi,           jlo, jhi,           &
                         aice    (:,:,iblk), vice    (:,:,iblk), & 
                         vsno    (:,:,iblk), tmask   (:,:,iblk), & 
                         strairxT(:,:,iblk), strairyT(:,:,iblk), & 
                         strairx (:,:,iblk), strairy (:,:,iblk), & 
                         tmass   (:,:,iblk), icetmask(:,:,iblk))

      enddo                     ! iblk
      !$OMP END PARALLEL DO

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (icetmask,          halo_info, &
                           field_loc_center,  field_type_scalar)
      call ice_timer_stop(timer_bound)

      !-----------------------------------------------------------------
      ! convert fields from T to U grid
      !-----------------------------------------------------------------

      call to_ugrid(tmass,umass)
      call to_ugrid(aice_init, aiu)

      !----------------------------------------------------------------
      ! Set wind stress to values supplied via NEMO or other forcing
      ! This wind stress is rotated on u grid and multiplied by aice
      !----------------------------------------------------------------
      call icepack_query_parameters(calc_strair_out=calc_strair)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (.not. calc_strair) then       
         strairx(:,:,:) = strax(:,:,:)
         strairy(:,:,:) = stray(:,:,:)
      else
         call t2ugrid_vector(strairx)
         call t2ugrid_vector(strairy)
      endif      

! tcraig, tcx, threading here leads to some non-reproducbile results and failures in icepack_ice_strength
! need to do more debugging
      !$TCXOMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks

      !-----------------------------------------------------------------
      ! more preparation for dynamics
      !-----------------------------------------------------------------

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         call dyn_prep2 (nx_block,             ny_block,             & 
                         ilo, ihi,             jlo, jhi,             &
                         icellt(iblk),         icellu(iblk),         & 
                         indxti      (:,iblk), indxtj      (:,iblk), & 
                         indxui      (:,iblk), indxuj      (:,iblk), & 
                         aiu       (:,:,iblk), umass     (:,:,iblk), & 
                         umassdti  (:,:,iblk), fcor_blk  (:,:,iblk), & 
                         umask     (:,:,iblk),                       & 
                         uocn      (:,:,iblk), vocn      (:,:,iblk), & 
                         strairx   (:,:,iblk), strairy   (:,:,iblk), & 
                         ss_tltx   (:,:,iblk), ss_tlty   (:,:,iblk), &  
                         icetmask  (:,:,iblk), iceumask  (:,:,iblk), & 
                         fm        (:,:,iblk), dt,                   & 
                         strtltx   (:,:,iblk), strtlty   (:,:,iblk), & 
                         strocnx   (:,:,iblk), strocny   (:,:,iblk), & 
                         strintx   (:,:,iblk), strinty   (:,:,iblk), & 
                         taubx     (:,:,iblk), tauby     (:,:,iblk), & 
                         waterx    (:,:,iblk), watery    (:,:,iblk), & 
                         forcex    (:,:,iblk), forcey    (:,:,iblk), & 
                         stressp_1 (:,:,iblk), stressp_2 (:,:,iblk), & 
                         stressp_3 (:,:,iblk), stressp_4 (:,:,iblk), & 
                         stressm_1 (:,:,iblk), stressm_2 (:,:,iblk), & 
                         stressm_3 (:,:,iblk), stressm_4 (:,:,iblk), & 
                         stress12_1(:,:,iblk), stress12_2(:,:,iblk), & 
                         stress12_3(:,:,iblk), stress12_4(:,:,iblk), & 
                         uvel_init (:,:,iblk), vvel_init (:,:,iblk), &
                         uvel      (:,:,iblk), vvel      (:,:,iblk), &
                         Tbu       (:,:,iblk))

         call calc_bfix (nx_block            , ny_block,             & 
                         icellu(iblk)        ,                       & 
                         indxui      (:,iblk), indxuj      (:,iblk), & 
                         umassdti  (:,:,iblk),                       & 
                         forcex    (:,:,iblk), forcey    (:,:,iblk), & 
                         uvel_init (:,:,iblk), vvel_init (:,:,iblk), &
                         bxfix     (:,:,iblk), byfix     (:,:,iblk))                         
                         
      !-----------------------------------------------------------------
      ! ice strength
      !-----------------------------------------------------------------

         strength(:,:,iblk) = c0  ! initialize
         do ij = 1, icellt(iblk)
            i = indxti(ij, iblk)
            j = indxtj(ij, iblk)
            call icepack_ice_strength (ncat,                 &
                                      aice    (i,j,  iblk), & 
                                      vice    (i,j,  iblk), & 
                                      aice0   (i,j,  iblk), & 
                                      aicen   (i,j,:,iblk), &  
                                      vicen   (i,j,:,iblk), & 
                                      strength(i,j,  iblk) )
         enddo  ! ij

         ! load velocity into array for boundary updates JFL move?
         fld2(:,:,1,iblk) = uvel(:,:,iblk)
         fld2(:,:,2,iblk) = vvel(:,:,iblk)

      enddo  ! iblk
      !$TCXOMP END PARALLEL DO

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (strength,           halo_info, &
                           field_loc_center,   field_type_scalar)
      ! velocities may have changed in dyn_prep2 ! JFL prends en compte la grille spherique qui se referme sur elle meme...
      call ice_HaloUpdate (fld2,               halo_info, & 
                           field_loc_NEcorner, field_type_vector)
      call ice_timer_stop(timer_bound)

      ! unload
      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks
         uvel(:,:,iblk) = fld2(:,:,1,iblk)
         vvel(:,:,iblk) = fld2(:,:,2,iblk)
      enddo
      !$OMP END PARALLEL DO

      if (maskhalo_dyn) then
         call ice_timer_start(timer_bound)
         halomask = 0
         where (iceumask) halomask = 1
         call ice_HaloUpdate (halomask,          halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_timer_stop(timer_bound)
         call ice_HaloMask(halo_info_mask, halo_info, halomask)
      endif

      !-----------------------------------------------------------------
      ! basal stress coefficients (landfast ice)
      !-----------------------------------------------------------------
      
      if (basalstress) then
         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks
            call basal_stress_coeff (nx_block,         ny_block,       &
                                     icellu  (iblk),                   &
                                     indxui(:,iblk),   indxuj(:,iblk), &
                                     vice(:,:,iblk),   aice(:,:,iblk), &
                                     hwater(:,:,iblk), Tbu(:,:,iblk))
         enddo
         !$OMP END PARALLEL DO
      endif
      
      !-----------------------------------------------------------------
      ! calc size of problem (ntot) and allocate arrays and vectors
      !-----------------------------------------------------------------
      
      ntot=0
      do iblk = 1, nblocks
        ntot = ntot + icellu(iblk)      
      enddo
      ntot = 2*ntot ! times 2 because of u and v
      
      allocate(bvec(ntot), sol(ntot), diagvec(ntot))
      
      !-----------------------------------------------------------------
      ! Start of nonlinear iteration
      !-----------------------------------------------------------------
      if (algo_nonlin == 1) then
         call picard_solver (icellt,   icellu,  &
                             indxti,   indxtj,  &
                             indxui,   indxuj,  &
                             fld2,              &
                             aiu,      ntot,    &
                             waterx,   watery,  & 
                             bxfix,    byfix,   &
                             umassdti, bvec,    & 
                             sol,      diagvec, &
                             fpresx,   fpresy,  &
                             halo_info_mask)
      elseif (algo_nonlin == 2) then
         call anderson_solver (icellt,   icellu,  &
                               indxti,   indxtj,  &
                               indxui,   indxuj,  &
                               fld2,              &
                               aiu,      ntot,    &
                               waterx,   watery,  & 
                               bxfix,    byfix,   &
                               umassdti, bvec,    & 
                               sol,      diagvec, &
                               fpresx,   fpresy,  &
                               halo_info_mask)
      endif
      !-----------------------------------------------------------------
      ! End of nonlinear iteration
      !-----------------------------------------------------------------

      deallocate(bvec, sol, diagvec)
      
      deallocate(fld2)
      if (maskhalo_dyn) call ice_HaloDestroy(halo_info_mask)
      
      !-----------------------------------------------------------------
      ! Compute deformations
      !-----------------------------------------------------------------
      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks
         call deformations (nx_block,             ny_block,             &
                            icellt(iblk),                               &
                            indxti      (:,iblk), indxtj      (:,iblk), &
                            uvel      (:,:,iblk), vvel      (:,:,iblk), &
                            dxt       (:,:,iblk), dyt       (:,:,iblk), &
                            cxp       (:,:,iblk), cyp       (:,:,iblk), &
                            cxm       (:,:,iblk), cym       (:,:,iblk), &
                            tarear    (:,:,iblk),                       &
                            shear     (:,:,iblk), divu      (:,:,iblk), &
                            rdg_conv  (:,:,iblk), rdg_shear (:,:,iblk))
      enddo
      !$OMP END PARALLEL DO
      
      ! phb: here we do halo updates for stresses (stressp_i, stressm_i, stress12_i, i=1..4),
      !      but stresses have not been updated ! (should be done in deformations ?)
      ! Force symmetry across the tripole seam
      if (trim(grid_type) == 'tripole') then
      if (maskhalo_dyn) then
         !-------------------------------------------------------
         ! set halomask to zero because ice_HaloMask always keeps
         ! local copies AND tripole zipper communication
         !-------------------------------------------------------
         halomask = 0
         call ice_HaloMask(halo_info_mask, halo_info, halomask)

         call ice_HaloUpdate_stress(stressp_1, stressp_3, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_3, stressp_1, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_2, stressp_4, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_4, stressp_2, halo_info_mask, &
                              field_loc_center,  field_type_scalar)

         call ice_HaloUpdate_stress(stressm_1, stressm_3, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_3, stressm_1, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_2, stressm_4, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_4, stressm_2, halo_info_mask, &
                              field_loc_center,  field_type_scalar)

         call ice_HaloUpdate_stress(stress12_1, stress12_3, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_3, stress12_1, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_2, stress12_4, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_4, stress12_2, halo_info_mask, &
                              field_loc_center,  field_type_scalar)

         call ice_HaloDestroy(halo_info_mask)
      else
         call ice_HaloUpdate_stress(stressp_1, stressp_3, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_3, stressp_1, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_2, stressp_4, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_4, stressp_2, halo_info, &
                              field_loc_center,  field_type_scalar)

         call ice_HaloUpdate_stress(stressm_1, stressm_3, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_3, stressm_1, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_2, stressm_4, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_4, stressm_2, halo_info, &
                              field_loc_center,  field_type_scalar)

         call ice_HaloUpdate_stress(stress12_1, stress12_3, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_3, stress12_1, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_2, stress12_4, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_4, stress12_2, halo_info, &
                              field_loc_center,  field_type_scalar)
      endif   ! maskhalo
      endif   ! tripole

      !-----------------------------------------------------------------
      ! ice-ocean stress
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks

         call dyn_finish                               & 
              (nx_block,           ny_block,           & 
               icellu      (iblk), Cdn_ocn (:,:,iblk), & 
               indxui    (:,iblk), indxuj    (:,iblk), & 
               uvel    (:,:,iblk), vvel    (:,:,iblk), & 
               uocn    (:,:,iblk), vocn    (:,:,iblk), & 
               aiu     (:,:,iblk), fm      (:,:,iblk), & 
               strintx (:,:,iblk), strinty (:,:,iblk), &
               strairx (:,:,iblk), strairy (:,:,iblk), &
               strocnx (:,:,iblk), strocny (:,:,iblk), & 
               strocnxT(:,:,iblk), strocnyT(:,:,iblk))

      enddo
      !$OMP END PARALLEL DO

      call u2tgrid_vector(strocnxT)    ! shift
      call u2tgrid_vector(strocnyT)

      call ice_timer_stop(timer_dynamics)    ! dynamics

      end subroutine imp_solver

!=======================================================================

! Solve nonlinear equation using Picard iterative solver
!
! author: JF Lemieux, A. Qaddouri and F. Dupont ECCC

      subroutine picard_solver (icellt,   icellu,  &
                                indxti,   indxtj,  &
                                indxui,   indxuj,  &
                                fld2,              &
                                aiu,      ntot,    &
                                waterx,   watery,  & 
                                bxfix,    byfix,   &
                                umassdti, bvec,    & 
                                sol,      diagvec, &
                                fpresx,   fpresy,  &
                                halo_info_mask)

      use ice_arrays_column, only: Cdn_ocn
      use ice_blocks, only: nx_block, ny_block
      use ice_boundary, only: ice_halo, ice_HaloUpdate
      use ice_domain, only: halo_info, maskhalo_dyn
      use ice_domain_size, only: max_blocks
      use ice_flux, only: uocn, vocn, fm, Tbu
      use ice_grid, only: dxt, dyt, dxhy, dyhx, cxp, cyp, cxm, cym, &
          tarear, uarear, tinyarea
      use ice_state, only: uvel, vvel, strength
      use ice_timers, only: timer_bound, ice_timer_start, ice_timer_stop

      integer (kind=int_kind), intent(in) :: & 
         ntot         ! size of problem for fgmres (for given cpu)

      integer (kind=int_kind), dimension(max_blocks), intent(in) :: & 
         icellt   , & ! no. of cells where icetmask = 1
         icellu       ! no. of cells where iceumask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block, max_blocks), intent(in) :: &
         indxti   , & ! compressed index in i-direction
         indxtj   , & ! compressed index in j-direction
         indxui   , & ! compressed index in i-direction
         indxuj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block,2,max_blocks), intent(inout) :: &
         fld2        ! work array for boundary updates

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), intent(in) :: &
         aiu      , & ! ice fraction on u-grid
         waterx   , & ! for ocean stress calculation, x (m/s)
         watery   , & ! for ocean stress calculation, y (m/s)
         bxfix    , & ! part of bx that is constant during Picard
         byfix    , & ! part of by that is constant during Picard
         umassdti     ! mass of U-cell/dte (kg/m^2 s)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), intent(inout) :: &
         fpresx   , & ! x fixed point residual vector, fx = uvel - uprev_k
         fpresy       ! y fixed point residual vector, fy = vvel - vprev_k

      real (kind=dbl_kind), dimension (ntot), intent(inout) :: &
         bvec     , & ! RHS vector for FGMRES
         sol      , & ! solution vector for FGMRES
         diagvec      ! diagonal of matrix A for preconditioners

      type (ice_halo) :: &
         halo_info_mask !  ghost cell update info for masked halo

      ! local variables

      integer (kind=int_kind) :: & 
         kOL            , & ! outer loop iteration
         iblk           , & ! block index
         icode          , & ! code for fgmres solver
         its            , & ! iteration nb for fgmres
         fgmres_its     , & ! final nb of fgmres iterations
         ierr               ! code for pgmres preconditioner !phb: needed?

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         uprev_k  , & ! uvel at previous Picard iteration
         vprev_k  , & ! vvel at previous Picard iteration
         vrel     , & ! coeff for tauw 
         Cb       , & ! seabed stress coeff
         bx       , & ! b vector
         by       , & ! b vector
         Diagu    , & ! Diagonal (u component) of the matrix A
         Diagv    , & ! Diagonal (v component) of the matrix A
         Au       , & ! matvec, Fx = Au - bx
         Av       , & ! matvec, Fy = Av - by
         Fx       , & ! x residual vector, Fx = Au - bx
         Fy           ! y residual vector, Fy = Av - by

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks,4):: &
         zetaD      ! zetaD = 2zeta (viscous coeff)

      real (kind=dbl_kind), dimension(nx_block,ny_block,8):: &
         stPrtmp,   & ! doit etre (nx_block,ny_block,max_blocks,8)???? PAs besoin des 2? reuse?
         Dstrtmp

      real (kind=dbl_kind), dimension (max_blocks) :: &
         L2norm       ! to compute l^2 norm of grid function

      real (kind=dbl_kind), allocatable :: & 
         vv(:,:), ww(:,:)  ! work arrays for FGMRES

      real (kind=dbl_kind), allocatable :: &
         wk11(:), wk22(:)  ! work vectors for FGMRES

      real (kind=dbl_kind) :: & 
         conv     , & ! ratio of current residual and initial residual for FGMRES !phb: needed for fgmres2
         tol      , & ! tolerance for nonlinear convergence: gammaNL * initial residual norm
         nlres_norm  , & ! norm of current nonlinear residual : F(x) = A(x)x -b(x)
         res_norm     ! residual norm for FGMRES

      character(len=*), parameter :: subname = '(picard_solver)'

      ! ! Allocate space for FGMRES work arrays
      ! allocate(wk11(ntot), wk22(ntot))
      ! allocate(vv(ntot,im_fgmres+1), ww(ntot,im_fgmres))
      ! 
      ! ! Start iterations
      ! do kOL = 1,maxits_nonlin        ! outer loop 
      ! 
      ! !-----------------------------------------------------------------
      ! ! Calc zetaD, Pr, Cb and vrel = f(uprev_k, vprev_k)
      ! !-----------------------------------------------------------------
      ! 
      ! !$OMP PARALLEL DO PRIVATE(iblk)
      !    do iblk = 1, nblocks
      ! 
      !       uprev_k(:,:,iblk) = uvel(:,:,iblk)
      !       vprev_k(:,:,iblk) = vvel(:,:,iblk)
      ! 
      !       call calc_zeta_Pr (nx_block           , ny_block,           &
      !                          icellt(iblk),                            & 
      !                          indxti   (:,iblk)  , indxtj(:,iblk),     & 
      !                          uprev_k  (:,:,iblk), vprev_k (:,:,iblk), & 
      !                          dxt      (:,:,iblk), dyt   (:,:,iblk),   & 
      !                          dxhy     (:,:,iblk), dyhx  (:,:,iblk),   & 
      !                          cxp      (:,:,iblk), cyp   (:,:,iblk),   & 
      !                          cxm      (:,:,iblk), cym   (:,:,iblk),   & 
      !                          tinyarea (:,:,iblk),                     & 
      !                          strength (:,:,iblk), zetaD (:,:,iblk,:) ,&
      !                          stPrtmp  (:,:,:) )                      
      ! 
      !       call calc_vrel_Cb (nx_block           , ny_block,           &
      !                          icellu       (iblk), Cdn_ocn (:,:,iblk), & 
      !                          indxui     (:,iblk), indxuj    (:,iblk), &
      !                          aiu      (:,:,iblk), Tbu     (:,:,iblk), &
      !                          uocn     (:,:,iblk), vocn    (:,:,iblk), &     
      !                          uprev_k  (:,:,iblk), vprev_k (:,:,iblk), & 
      !                          vrel     (:,:,iblk), Cb      (:,:,iblk))
      ! 
      ! !     prepare b vector (RHS)                                                
      !       call calc_bvec (nx_block           , ny_block,           &
      !                       icellu       (iblk),                     & 
      !                       indxui     (:,iblk), indxuj    (:,iblk), &
      !                       stPrtmp  (:,:,:)   , Cdn_ocn (:,:,iblk), &
      !                       aiu      (:,:,iblk), uarear  (:,:,iblk), & 
      !                       uocn     (:,:,iblk), vocn    (:,:,iblk), &     
      !                       waterx   (:,:,iblk), watery  (:,:,iblk), & 
      !                       uprev_k  (:,:,iblk), vprev_k (:,:,iblk), & 
      !                       bxfix    (:,:,iblk), byfix   (:,:,iblk), &
      !                       bx       (:,:,iblk), by      (:,:,iblk), &
      !                       vrel     (:,:,iblk))
      ! 
      ! !     prepare precond matrix
      !      if (precond .gt. 1) then
      ! 
      !      call formDiag_step1  (nx_block           , ny_block,       & ! D term due to rheology
      !                            icellu       (iblk),                 &
      !                            indxui     (:,iblk), indxuj(:,iblk), &
      !                            dxt      (:,:,iblk), dyt (:,:,iblk), & 
      !                            dxhy     (:,:,iblk), dyhx(:,:,iblk), & 
      !                            cxp      (:,:,iblk), cyp (:,:,iblk), & 
      !                            cxm      (:,:,iblk), cym (:,:,iblk), & 
      !                            zetaD (:,:,iblk,:) , Dstrtmp (:,:,:) )
      ! 
      !      call formDiag_step2 (nx_block           , ny_block,           &
      !                           icellu       (iblk),                     & 
      !                           indxui     (:,iblk), indxuj    (:,iblk), &
      !                           Dstrtmp  (:,:,:)   , vrel    (:,:,iblk), &
      !                           umassdti (:,:,iblk),                     & 
      !                           uarear   (:,:,iblk), Cb      (:,:,iblk), & 
      !                           Diagu    (:,:,iblk), Diagv   (:,:,iblk))         
      ! 
      !     endif                     
      ! 
      !    enddo
      !    !$OMP END PARALLEL DO                            
      ! 
      !    ! Compute nonlinear residual norm
      !    !$OMP PARALLEL DO PRIVATE(iblk)
      !    do iblk = 1, nblocks
      !       call matvec (nx_block             , ny_block,            &
      !                    icellu   (iblk)      , icellt   (iblk)    , & 
      !                    indxui   (:,iblk)    , indxuj   (:,iblk)  , &
      !                    indxti   (:,iblk)    , indxtj   (:,iblk)  , &
      !                    dxt      (:,:,iblk)  , dyt      (:,:,iblk), & 
      !                    dxhy     (:,:,iblk)  , dyhx     (:,:,iblk), & 
      !                    cxp      (:,:,iblk)  , cyp      (:,:,iblk), & 
      !                    cxm      (:,:,iblk)  , cym      (:,:,iblk), & 
      !                    uvel     (:,:,iblk)  , vvel     (:,:,iblk), &      
      !                    vrel     (:,:,iblk)  , Cb       (:,:,iblk), &  
      !                    zetaD    (:,:,iblk,:),                      &
      !                    umassdti (:,:,iblk)  , fm       (:,:,iblk), & 
      !                    uarear   (:,:,iblk)  ,                      & 
      !                    Au       (:,:,iblk)  , Av       (:,:,iblk))
      !       call residual_vec (nx_block           , ny_block,           &
      !                          icellu       (iblk),                     & 
      !                          indxui     (:,iblk), indxuj    (:,iblk), &
      !                          bx       (:,:,iblk), by      (:,:,iblk), &
      !                          Au       (:,:,iblk), Av      (:,:,iblk), &
      !                          Fx       (:,:,iblk), Fy      (:,:,iblk), &
      !                          L2norm(iblk))
      !    enddo
      !    !$OMP END PARALLEL DO
      !    nlres_norm = sqrt(sum(L2norm))
      !    if (monitor_nonlin) then
      !       write(nu_diag, '(a,i4,a,d26.16)') "monitor_nonlin: iter_nonlin= ", kOL, &
      !                                         " nonlin_res_L2norm= ", nlres_norm
      !    endif
      !    ! Compute relative tolerance at first iteration
      !    if (kOL == 1) then
      !       tol = gammaNL*nlres_norm
      !    endif
      !    ! Check for nonlinear convergence
      !    if (nlres_norm < tol) then
      !       exit
      !    endif
      ! 
      ! !-----------------------------------------------------------------------
      ! !     prep F G M R E S 
      ! !-----------------------------------------------------------------------                             
      ! 
      ! icode  = 0
      ! !      its    = 0 
      ! 
      !    ! form b vector from matrices (nblocks matrices)      
      !    call arrays_to_vec (nx_block, ny_block, nblocks,    &
      !                        max_blocks, icellu (:), ntot,   & 
      !                        indxui      (:,:), indxuj(:,:), &
      !                        bx        (:,:,:), by  (:,:,:), &
      !                        bvec(:))
      !    ! form sol vector for fgmres (sol is iniguess at the beginning)        
      !    call arrays_to_vec (nx_block, ny_block, nblocks,      &
      !                        max_blocks, icellu (:), ntot,   &  
      !                        indxui    (:,:), indxuj(:,:),     &
      !                        uprev_k (:,:,:), vprev_k (:,:,:), &
      !                        sol(:))
      ! 
      !    ! form matrix diagonal as a vector from Diagu and Diagv arrays      
      !    call arrays_to_vec (nx_block, ny_block, nblocks,    &
      !                        max_blocks, icellu (:), ntot,   & 
      !                        indxui      (:,:), indxuj(:,:), &
      !                        Diagu     (:,:,:), Diagv(:,:,:),&
      !                        diagvec(:))                             
      ! 
      ! !-----------------------------------------------------------------------
      ! !     F G M R E S   L O O P
      ! !-----------------------------------------------------------------------
      ! 1    continue
      ! !-----------------------------------------------------------------------
      ! 
      ! call fgmres (ntot,im_fgmres,bvec,sol,its,vv,ww,wk11,wk22, &
      !              gamma, maxits_fgmres, monitor_fgmres,   &
      !              icode,fgmres_its, res_norm)
      ! 
      ! if (icode == 1) then
      ! 
      !    if (precond .eq. 1) then
      ! 
      !      wk22(:)=wk11(:) ! precond=identity
      ! 
      !    elseif (precond .eq. 2) then ! use diagonal of A for precond step
      ! 
      !      call precond_diag (ntot,            & 
      !                         diagvec (:),     &
      !                         wk11 (:), wk22 (:) )
      ! 
      !    elseif (precond .eq. 3) then
      ! 
      !     call pgmres (nx_block,    ny_block,    nblocks       , &
      !                  max_blocks         , icellu   (:)       , & 
      !                  indxui   (:,:)     , indxuj   (:,:)     , &
      !                  icellt   (:)                            , & 
      !                  indxti   (:,:)     , indxtj   (:,:)     , &
      !                  dxt      (:,:,:)   , dyt      (:,:,:)   , & 
      !                  dxhy     (:,:,:)   , dyhx     (:,:,:)   , & 
      !                  cxp      (:,:,:)   , cyp      (:,:,:)   , & 
      !                  cxm      (:,:,:)   , cym      (:,:,:)   , & 
      !                  vrel     (:,:,:)   , Cb       (:,:,:)   , &  
      !                  zetaD    (:,:,:,:) ,                      &
      !                  umassdti (:,:,:)   , fm       (:,:,:)   , & 
      !                  uarear   (:,:,:)   , diagvec(:)         , &
      !                  wk22     (:)       , wk11(:)            , &
      !                  ntot               , im_pgmres          , &
      !                  epsprecond         , maxits_pgmres      , &
      !                  monitor_pgmres     , ierr )         
      !    endif ! precond
      ! 
      !    goto 1
      ! 
      ! elseif (icode >= 2) then
      ! 
      !    call vec_to_arrays (nx_block, ny_block, nblocks,      &
      !                        max_blocks, icellu (:), ntot,     & 
      !                        indxui    (:,:), indxuj(:,:),     &
      !                        wk11 (:),                         &
      !                        uvel (:,:,:), vvel (:,:,:))    
      ! 
      !    ! JFL halo update could be in subroutine...                    
      !    !$OMP PARALLEL DO PRIVATE(iblk) 
      !    do iblk = 1, nblocks                             
      !       fld2(:,:,1,iblk) = uvel(:,:,iblk)
      !       fld2(:,:,2,iblk) = vvel(:,:,iblk)            
      !    enddo
      !    !$OMP END PARALLEL DO                           
      ! 
      !    call ice_HaloUpdate (fld2,               halo_info, & 
      !                         field_loc_NEcorner, field_type_vector)
      ! 
      !    !$OMP PARALLEL DO PRIVATE(iblk)
      !    do iblk = 1, nblocks
      !       uvel(:,:,iblk) = fld2(:,:,1,iblk)
      !       vvel(:,:,iblk) = fld2(:,:,2,iblk)
      !    enddo
      !    !$OMP END PARALLEL DO                             
      ! 
      !    !$OMP PARALLEL DO PRIVATE(iblk)
      !    do iblk = 1, nblocks                                  
      ! 
      !     call matvec (nx_block             , ny_block,            &
      !                  icellu   (iblk)      , icellt   (iblk)    , & 
      !                  indxui   (:,iblk)    , indxuj   (:,iblk)  , &
      !                  indxti   (:,iblk)    , indxtj   (:,iblk)  , &
      !                  dxt      (:,:,iblk)  , dyt      (:,:,iblk), & 
      !                  dxhy     (:,:,iblk)  , dyhx     (:,:,iblk), & 
      !                  cxp      (:,:,iblk)  , cyp      (:,:,iblk), & 
      !                  cxm      (:,:,iblk)  , cym      (:,:,iblk), &
      !                  uvel     (:,:,iblk)  , vvel     (:,:,iblk), &      
      !                  vrel     (:,:,iblk)  , Cb       (:,:,iblk), &  
      !                  zetaD    (:,:,iblk,:),                      &
      !                  umassdti (:,:,iblk)  , fm       (:,:,iblk), & 
      !                  uarear   (:,:,iblk)  ,                      & 
      !                  Au       (:,:,iblk)  , Av       (:,:,iblk))                         
      ! 
      !    enddo
      !    !$OMP END PARALLEL DO 
      ! 
      !    ! form wk2 from Au and Av arrays        
      !    call arrays_to_vec (nx_block, ny_block, nblocks,      &
      !                        max_blocks, icellu (:), ntot,     & 
      !                        indxui    (:,:), indxuj(:,:),     &
      !                        Au      (:,:,:), Av    (:,:,:),   &
      !                        wk22(:))    
      ! 
      !       goto 1
      ! 
      ! endif ! icode
      ! 
      ! !-----------------------------------------------------------------------
      ! !     Put vector sol in uvel and vvel arrays
      ! !-----------------------------------------------------------------------
      ! 
      !    call vec_to_arrays (nx_block, ny_block, nblocks,      &
      !                        max_blocks, icellu (:), ntot,     & 
      !                        indxui    (:,:), indxuj(:,:),     &
      !                        sol (:),                          &
      !                        uvel (:,:,:), vvel (:,:,:))    
      ! 
      !    !$OMP PARALLEL DO PRIVATE(iblk)
      ! !         do iblk = 1, nblocks
      ! !              uvel(:,:,iblk) = (c1-krelax)*uprev_k(:,:,iblk) + krelax*uvel(:,:,iblk)
      ! !              vvel(:,:,iblk) = (c1-krelax)*vprev_k(:,:,iblk) + krelax*vvel(:,:,iblk)
      ! !         enddo
      !    !$OMP END PARALLEL DO  
      ! 
      !    !$OMP PARALLEL DO PRIVATE(iblk)
      !    do iblk = 1, nblocks                             
      ! 
      !       ! load velocity into array for boundary updates
      !       fld2(:,:,1,iblk) = uvel(:,:,iblk)
      !       fld2(:,:,2,iblk) = vvel(:,:,iblk)            
      ! 
      !    enddo
      !    !$OMP END PARALLEL DO                           
      ! 
      !    call ice_timer_start(timer_bound)
      !    if (maskhalo_dyn) then
      !       call ice_HaloUpdate (fld2,               halo_info_mask, &
      !                            field_loc_NEcorner, field_type_vector)
      !    else
      !       call ice_HaloUpdate (fld2,               halo_info, &
      !                            field_loc_NEcorner, field_type_vector)
      !    endif
      !    call ice_timer_stop(timer_bound)
      ! 
      !    ! unload
      !    !$OMP PARALLEL DO PRIVATE(iblk)
      !    do iblk = 1, nblocks
      !       uvel(:,:,iblk) = fld2(:,:,1,iblk)
      !       vvel(:,:,iblk) = fld2(:,:,2,iblk)
      !    enddo
      !    !$OMP END PARALLEL DO
      ! 
      !    ! Compute fixed point residual norm
      !    !$OMP PARALLEL DO PRIVATE(iblk)
      !    do iblk = 1, nblocks
      !       fpresx(:,:,iblk) = uvel(:,:,iblk) - uprev_k(:,:,iblk)
      !       fpresy(:,:,iblk) = vvel(:,:,iblk) - vprev_k(:,:,iblk)
      !       call calc_L2norm_squared (nx_block        , ny_block,         &
      !                                 icellu    (iblk),                   & 
      !                                 indxui  (:,iblk), indxuj  (:,iblk), &
      !                                 fpresx(:,:,iblk), fpresy(:,:,iblk), &
      !                                 L2norm    (iblk))
      !    enddo
      !    !$OMP END PARALLEL DO
      !    if (monitor_nonlin) then
      !       write(nu_diag, '(a,i4,a,d26.16)') "monitor_nonlin: iter_nonlin= ", kOL, &
      !                                         " fixed_point_res_L2norm= ", sqrt(sum(L2norm))
      !    endif
      ! 
      ! enddo                     ! outer loop
      ! 
      ! ! deallocate FGMRES work arrays
      ! deallocate(wk11, wk22, vv, ww)

      end subroutine picard_solver

!=======================================================================

! Solve nonlinear equation using fixed point iteration, accelerated with 
! Anderson acceleration
!
! author: P. Blain ECCC

      subroutine anderson_solver (icellt,   icellu,  &
                                  indxti,   indxtj,  &
                                  indxui,   indxuj,  &
                                  fld2,              &
                                  aiu,      ntot,    &
                                  waterx,   watery,  & 
                                  bxfix,    byfix,   &
                                  umassdti, bvec,    & 
                                  sol,      diagvec, &
                                  fpresx,   fpresy,  &
                                  halo_info_mask)

      use ice_arrays_column, only: Cdn_ocn
      use ice_blocks, only: nx_block, ny_block
      use ice_boundary, only: ice_halo, ice_HaloUpdate
      use ice_constants, only: c1
      use ice_domain, only: halo_info, maskhalo_dyn
      use ice_domain_size, only: max_blocks
      use ice_flux, only:   uocn, vocn, fm, Tbu
      use ice_grid, only: dxt, dyt, dxhy, dyhx, cxp, cyp, cxm, cym, &
          uarear, tinyarea
      use ice_state, only: uvel, vvel, strength
      use ice_timers, only: timer_bound, ice_timer_start, ice_timer_stop

      integer (kind=int_kind), intent(in) :: & 
         ntot         ! size of problem for fgmres (for given cpu)

      integer (kind=int_kind), dimension(max_blocks), intent(in) :: & 
         icellt   , & ! no. of cells where icetmask = 1
         icellu       ! no. of cells where iceumask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block, max_blocks), intent(in) :: &
         indxti   , & ! compressed index in i-direction
         indxtj   , & ! compressed index in j-direction
         indxui   , & ! compressed index in i-direction
         indxuj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block,2,max_blocks), intent(inout) :: &
         fld2        ! work array for boundary updates

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), intent(in) :: &
         aiu      , & ! ice fraction on u-grid
         waterx   , & ! for ocean stress calculation, x (m/s)
         watery   , & ! for ocean stress calculation, y (m/s)
         bxfix    , & ! part of bx that is constant during Picard
         byfix    , & ! part of by that is constant during Picard
         umassdti     ! mass of U-cell/dte (kg/m^2 s)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), intent(inout) :: &
         fpresx   , & ! x fixed point residual vector, fx = uvel - uprev_k
         fpresy       ! y fixed point residual vector, fy = vvel - vprev_k

      real (kind=dbl_kind), dimension (ntot), intent(inout) :: &
         bvec     , & ! RHS vector for FGMRES
         sol      , & ! current approximate solution
         diagvec      ! diagonal of matrix A for preconditioners

      type (ice_halo) :: &
         halo_info_mask !  ghost cell update info for masked halo

      ! local variables

      integer (kind=int_kind) :: & 
         it_nl      , & ! nonlinear loop iteration index
         res_num    , & ! current number of stored residuals
         j          , & ! iteration index for QR update
         iblk       , & ! block index
         nbiter         ! number of FGMRES iterations performed

      integer (kind=int_kind), parameter :: &
         inc = 1        ! increment value for BLAS calls

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         uprev_k  , & ! uvel at previous Picard iteration
         vprev_k  , & ! vvel at previous Picard iteration
         ulin     , & ! uvel to linearize vrel
         vlin     , & ! vvel to linearize vrel
         vrel     , & ! coeff for tauw 
         Cb       , & ! seabed stress coeff
         bx       , & ! b vector
         by       , & ! b vector
         Diagu    , & ! Diagonal (u component) of the matrix A
         Diagv    , & ! Diagonal (v component) of the matrix A
         Au       , & ! matvec, Fx = Au - bx
         Av       , & ! matvec, Fy = Av - by
         Fx       , & ! x residual vector, Fx = Au - bx
         Fy       , & ! y residual vector, Fy = Av - by
         solx     , & ! solution of FGMRES (x components)
         soly         ! solution of FGMRES (y components)

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks,4):: &
         zetaD        ! zetaD = 2zeta (viscous coeff)

      real (kind=dbl_kind), dimension(nx_block,ny_block,8):: &
         stPrtmp,   & ! doit etre (nx_block,ny_block,max_blocks,8)???? PAs besoin des 2? reuse?
         Dstrtmp

      real (kind=dbl_kind), dimension (max_blocks) :: &
         L2norm       ! to compute l^2 norm of grid function

      real (kind=dbl_kind), dimension (ntot) :: &
         res        , & ! current residual
         res_old    , & ! previous residual
         res_diff   , & ! difference between current and previous residuals
         fpfunc     , & ! current value of fixed point function
         fpfunc_old , & ! previous value of fixed point function
         Fvec       , & ! (Fx,Fy) (nonlinear residual) as vector
         tmp            ! temporary vector for BLAS calls
      
      real (kind=dbl_kind), dimension(ntot,im_andacc) :: &
         Q        , & ! Q factor for QR factorization of F (residuals) matrix
         G_diff       ! Matrix containing the differences of g(x) (fixed point function) evaluations
      
      real (kind=dbl_kind), dimension(im_andacc,im_andacc) :: &
         R            ! R factor for QR factorization of F (residuals) matrix
      
      real (kind=dbl_kind), dimension(im_andacc) :: &
         rhs_tri  , & ! right hand side vector for matrix-vector product
         coeffs       ! coeffs used to combine previous solutions

      real (kind=dbl_kind) :: & 
         tol         , & ! tolerance for fixed point convergence: reltol_andacc * (initial fixed point residual norm)
         tol_nl      , & ! tolerance for nonlinear convergence: gammaNL * (initial nonlinear residual norm)
         fpres_norm  , & ! norm of current fixed point residual : f(x) = g(x) - x
         prog_norm   , & ! norm of difference between current and previous solution
         nlres_norm  , & ! norm of current nonlinear residual : F(x) = A(x)x -b(x)
         ddot, dnrm2 , & ! BLAS functions
         conv            ! needed for FGMRES !phb keep ?

      character(len=*), parameter :: subname = '(anderson_solver)'
      
      ! Initialization
      res_num = 0
      
      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks
         uprev_k(:,:,iblk) = uvel(:,:,iblk)
         vprev_k(:,:,iblk) = vvel(:,:,iblk)
      enddo
      !$OMP END PARALLEL DO
      
      ! Start iterations
      do it_nl = 0, maxits_nonlin        ! nonlinear iteration loop 
         ! Compute quantities needed for computing PDE residual and g(x) (fixed point map)
         !-----------------------------------------------------------------
         ! Calc zetaD, Pr, Cb and vrel = f(uprev_k, vprev_k)
         !-----------------------------------------------------------------
         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks
            
            if (use_mean_vrel) then
               ulin(:,:,iblk) = p5*uprev_k(:,:,iblk) + p5*uvel(:,:,iblk)
               vlin(:,:,iblk) = p5*vprev_k(:,:,iblk) + p5*vvel(:,:,iblk)
            else
               ulin(:,:,iblk) = uvel(:,:,iblk)
               vlin(:,:,iblk) = vvel(:,:,iblk)
            endif
            uprev_k(:,:,iblk) = uvel(:,:,iblk)
            vprev_k(:,:,iblk) = vvel(:,:,iblk)
            
            call calc_zeta_Pr (nx_block           , ny_block,           &
                               icellt(iblk),                            & 
                               indxti   (:,iblk)  , indxtj(:,iblk),     & 
                               uprev_k  (:,:,iblk), vprev_k (:,:,iblk), & 
                               dxt      (:,:,iblk), dyt   (:,:,iblk),   & 
                               dxhy     (:,:,iblk), dyhx  (:,:,iblk),   & 
                               cxp      (:,:,iblk), cyp   (:,:,iblk),   & 
                               cxm      (:,:,iblk), cym   (:,:,iblk),   & 
                               tinyarea (:,:,iblk),                     & 
                               strength (:,:,iblk), zetaD (:,:,iblk,:) ,&
                               stPrtmp  (:,:,:) )                      
            
            call calc_vrel_Cb (nx_block           , ny_block,           &
                               icellu       (iblk), Cdn_ocn (:,:,iblk), & 
                               indxui     (:,iblk), indxuj    (:,iblk), &
                               aiu      (:,:,iblk), Tbu     (:,:,iblk), &
                               uocn     (:,:,iblk), vocn    (:,:,iblk), &     
                               ulin     (:,:,iblk), vlin    (:,:,iblk), &
                               vrel     (:,:,iblk), Cb      (:,:,iblk))
            
            ! prepare b vector (RHS)
            call calc_bvec (nx_block           , ny_block,           &
                            icellu       (iblk),                     & 
                            indxui     (:,iblk), indxuj    (:,iblk), &
                            stPrtmp  (:,:,:)   , Cdn_ocn (:,:,iblk), &
                            aiu      (:,:,iblk), uarear  (:,:,iblk), & 
                            uocn     (:,:,iblk), vocn    (:,:,iblk), &     
                            waterx   (:,:,iblk), watery  (:,:,iblk), & 
                            ulin     (:,:,iblk), vlin    (:,:,iblk), &
                            bxfix    (:,:,iblk), byfix   (:,:,iblk), &
                            bx       (:,:,iblk), by      (:,:,iblk), &
                            vrel     (:,:,iblk))
            
            ! Compute nonlinear residual norm (PDE residual)
            call matvec (nx_block             , ny_block,            &
                         icellu   (iblk)      , icellt   (iblk)    , &
                         indxui   (:,iblk)    , indxuj   (:,iblk)  , &
                         indxti   (:,iblk)    , indxtj   (:,iblk)  , &
                         dxt      (:,:,iblk)  , dyt      (:,:,iblk), &
                         dxhy     (:,:,iblk)  , dyhx     (:,:,iblk), &
                         cxp      (:,:,iblk)  , cyp      (:,:,iblk), &
                         cxm      (:,:,iblk)  , cym      (:,:,iblk), &
                         uprev_k  (:,:,iblk)  , vprev_k  (:,:,iblk), &
                         vrel     (:,:,iblk)  , Cb       (:,:,iblk), &
                         zetaD    (:,:,iblk,:),                      &
                         umassdti (:,:,iblk)  , fm       (:,:,iblk), &
                         uarear   (:,:,iblk)  ,                      &
                         Au       (:,:,iblk)  , Av       (:,:,iblk))
            call residual_vec (nx_block           , ny_block,           &
                               icellu       (iblk),                     & 
                               indxui     (:,iblk), indxuj    (:,iblk), &
                               bx       (:,:,iblk), by      (:,:,iblk), &
                               Au       (:,:,iblk), Av      (:,:,iblk), &
                               Fx       (:,:,iblk), Fy      (:,:,iblk), &
                               L2norm(iblk))
         enddo
         !$OMP END PARALLEL DO
         nlres_norm = sqrt(sum(L2norm))  ! phb: change after parallelization
         if (monitor_nonlin) then
            write(nu_diag, '(a,i4,a,d26.16)') "monitor_nonlin: iter_nonlin= ", it_nl, &
                                              " nonlin_res_L2norm= ", nlres_norm
         endif
         ! Compute relative tolerance at first iteration
         if (it_nl == 0) then
            tol_nl = gammaNL*nlres_norm
         endif
         
         ! Check for nonlinear convergence
         if (nlres_norm < tol_nl) then
            exit
         endif
         
         ! Put initial guess for FGMRES in solx,soly
         solx = uprev_k
         soly = vprev_k
         
         ! Compute fixed point map g(x)
         if (fpfunc_andacc == 1) then
            ! g_1(x) = FGMRES(A(x), b(x))
            
            ! Prepare precond matrix
            if (precond .gt. 1) then
               !$OMP PARALLEL DO PRIVATE(iblk)
               do iblk = 1, nblocks
                  call formDiag_step1 (nx_block           , ny_block,       & ! D term due to rheology
                                       icellu       (iblk),                 &
                                       indxui     (:,iblk), indxuj(:,iblk), &
                                       dxt      (:,:,iblk), dyt (:,:,iblk), & 
                                       dxhy     (:,:,iblk), dyhx(:,:,iblk), & 
                                       cxp      (:,:,iblk), cyp (:,:,iblk), & 
                                       cxm      (:,:,iblk), cym (:,:,iblk), & 
                                       zetaD (:,:,iblk,:) , Dstrtmp (:,:,:) )
                  call formDiag_step2 (nx_block           , ny_block,           &
                                       icellu       (iblk),                     & 
                                       indxui     (:,iblk), indxuj    (:,iblk), &
                                       Dstrtmp  (:,:,:)   , vrel    (:,:,iblk), &
                                       umassdti (:,:,iblk),                     & 
                                       uarear   (:,:,iblk), Cb      (:,:,iblk), & 
                                       Diagu    (:,:,iblk), Diagv   (:,:,iblk))
               enddo
               !$OMP END PARALLEL DO
            endif
            
            ! FGMRES linear solver
            call fgmres (zetaD,               &
                         Cb,         vrel,    &
                         umassdti,            &
                         solx,       soly,    &
                         bx,         by,      &
                         Diagu,      Diagv,   &
                         gamma, im_fgmres, &
                         maxits_fgmres, nbiter, conv)
            ! Put FGMRES solution solx,soly in fpfunc vector
            call arrays_to_vec (nx_block, ny_block, nblocks,      &
                                max_blocks, icellu (:), ntot,     &
                                indxui    (:,:), indxuj(:,:),     &
                                solx    (:,:,:), soly  (:,:,:),   &
                                fpfunc(:))
            ! FGMRES linear solver (solution is in fpfunc)
            ! fpfunc = sol
            ! call fgmres_solver (ntot,   bvec,     &
            !                     fpfunc, diagvec,  &
            !                     icellt, icellu,   &
            !                     indxti, indxtj,   &
            !                     indxui, indxuj,   &
            !                     zetaD,            &
            !                     Cb,     vrel,     &
            !                     aiu,    umassdti, &
            !                     fld2)

         elseif (fpfunc_andacc == 2) then
            ! g_2(x) = x - A(x)x + b(x) = x - F(x)
         endif

         ! Compute residual
         res = fpfunc - sol
         fpres_norm = dnrm2(size(res), res, inc)
         if (monitor_nonlin) then
            ! commented code is to compare fixed_point_res_L2norm BFB with progress_res_L2norm
            ! (should be BFB if Picard iteration is used)
            ! call vec_to_arrays (nx_block, ny_block, nblocks,      &
            !                     max_blocks, icellu (:), ntot,     & 
            !                     indxui    (:,:), indxuj(:,:),     &
            !                     res (:),                          &
            !                     fpresx (:,:,:), fpresy (:,:,:))
            ! !$OMP PARALLEL DO PRIVATE(iblk)
            ! do iblk = 1, nblocks
            ! call calc_L2norm_squared (nx_block        , ny_block,         &
            !                           icellu    (iblk),                   & 
            !                           indxui  (:,iblk), indxuj  (:,iblk), &
            !                           fpresx(:,:,iblk), fpresy(:,:,iblk), &
            !                           L2norm    (iblk))
            ! enddo
            ! !$OMP END PARALLEL DO
            ! write(nu_diag, '(a,i4,a,d26.16)') "monitor_nonlin: iter_nonlin= ", it_nl, &
            !                                   " fixed_point_res_L2norm= ", sqrt(sum(L2norm))
            write(nu_diag, '(a,i4,a,d26.16)') "monitor_nonlin: iter_nonlin= ", it_nl, &
                                              " fixed_point_res_L2norm= ", fpres_norm
         endif
         
         ! Not used for now (only nonlinear residual is checked)
         ! ! Store initial residual norm
         ! if (it_nl == 0) then
         !    tol = reltol_andacc*fpres_norm
         ! endif
         ! 
         ! ! Check residual 
         ! if (fpres_norm < tol) then
         !    exit
         ! endif

         if (im_andacc == 0 .or. it_nl < start_andacc) then 
            ! Simple fixed point (Picard) iteration in this case
            sol = fpfunc
         else
            ! Begin Anderson acceleration
            if (it_nl > start_andacc) then
               ! Update residual difference vector
               res_diff = res - res_old
               ! Update fixed point function difference matrix
               if (res_num < im_andacc) then
                  ! Add column
                  G_diff(:,res_num+1) = fpfunc - fpfunc_old
               else
                  ! Delete first column and add column
                  G_diff(:,1:res_num-1) = G_diff(:,2:res_num)
                  G_diff(:,res_num) = fpfunc - fpfunc_old
               endif
               res_num = res_num + 1
            endif
            res_old = res
            fpfunc_old = fpfunc
            if (res_num == 0) then
               sol = fpfunc
            else
               if (res_num == 1) then
                  ! Initialize QR factorization
                  R(1,1) = dnrm2(size(res_diff), res_diff, inc)
                  Q(:,1) = res_diff/R(1,1)
               else
                  if (res_num > im_andacc) then
                     ! Update factorization since 1st column was deleted
                     call qr_delete(Q,R)
                     res_num = res_num - 1
                  endif
                  ! Update QR factorization for new column
                  do j = 1, res_num - 1
                     R(j,res_num) = ddot(ntot, Q(:,j), inc, res_diff, inc)
                     res_diff = res_diff - R(j,res_num) * Q(:,j)
                  enddo
                  R(res_num, res_num) = dnrm2(size(res_diff) ,res_diff, inc)
                  Q(:,res_num) = res_diff / R(res_num, res_num)
               endif
               ! phb: here, drop more columns to improve conditioning
               ! if (droptol) then
               
               ! endif
               ! Solve least square problem for coefficients
               ! 1. Compute rhs_tri = Q^T * res
               call dgemv ('t', size(Q,1), res_num, c1, Q(:,1:res_num), size(Q,1), res, inc, c0, rhs_tri, inc)
               ! 2. Solve R*coeffs = rhs_tri, puts result in rhs_tri
               call dtrsv ('u', 'n', 'n', res_num, R(1:res_num,1:res_num), res_num, rhs_tri, inc)
               coeffs = rhs_tri
               ! Update approximate solution: x = fpfunc - G_diff*coeffs, puts result in fpfunc
               call dgemv ('n', size(G_diff,1), res_num, -c1, G_diff(:,1:res_num), size(G_diff,1), coeffs, inc, c1, fpfunc, inc)
               sol = fpfunc
               ! Apply damping
               if (damping_andacc > 0 .and. damping_andacc /= 1) then
                  ! x = x - (1-beta) (res - Q*R*coeffs)
                  
                  ! tmp = R*coeffs
                  call dgemv ('n', res_num, res_num, c1, R(1:res_num,1:res_num), res_num, coeffs, inc, c0, tmp, inc)
                  ! res = res - Q*tmp
                  call dgemv ('n', size(Q,1), res_num, -c1, Q(:,1:res_num), size(Q,1), tmp, inc, c1, res, inc)
                  ! x = x - (1-beta)*res
                  sol = sol - (1-damping_andacc)*res
               endif
            endif
         endif
         
         !-----------------------------------------------------------------------
         !     Put vector sol in uvel and vvel arrays
         !-----------------------------------------------------------------------
         call vec_to_arrays (nx_block, ny_block, nblocks,      &
                             max_blocks, icellu (:), ntot,     & 
                             indxui    (:,:), indxuj(:,:),     &
                             sol (:),                          &
                             uvel (:,:,:), vvel (:,:,:))
         ! Load velocity into array for boundary updates
         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks
            fld2(:,:,1,iblk) = uvel(:,:,iblk)
            fld2(:,:,2,iblk) = vvel(:,:,iblk)
         enddo
         !$OMP END PARALLEL DO

         call ice_timer_start(timer_bound)
         if (maskhalo_dyn) then
            call ice_HaloUpdate (fld2,               halo_info_mask, &
                                 field_loc_NEcorner, field_type_vector)
         else
            call ice_HaloUpdate (fld2,               halo_info, &
                                 field_loc_NEcorner, field_type_vector)
         endif
         call ice_timer_stop(timer_bound)

         ! Unload
         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks
            uvel(:,:,iblk) = fld2(:,:,1,iblk)
            vvel(:,:,iblk) = fld2(:,:,2,iblk)
         enddo
         !$OMP END PARALLEL DO
         
         ! Compute fixed point residual norm
         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks
            fpresx(:,:,iblk) = uvel(:,:,iblk) - uprev_k(:,:,iblk)
            fpresy(:,:,iblk) = vvel(:,:,iblk) - vprev_k(:,:,iblk)
            call calc_L2norm_squared (nx_block        , ny_block,         &
                                      icellu    (iblk),                   & 
                                      indxui  (:,iblk), indxuj  (:,iblk), &
                                      fpresx(:,:,iblk), fpresy(:,:,iblk), &
                                      L2norm    (iblk))
         enddo
         !$OMP END PARALLEL DO
         if (monitor_nonlin) then
            write(nu_diag, '(a,i4,a,d26.16)') "monitor_nonlin: iter_nonlin= ", it_nl, &
                                              " progress_res_L2norm= ", sqrt(sum(L2norm))
         endif
         
      enddo ! nonlinear iteration loop
      
      end subroutine anderson_solver

!=======================================================================

! Driver for the FGMRES linear solver

      subroutine fgmres_solver (ntot,   bvec,     &
                                sol,    diagvec,  &
                                icellt, icellu,   &
                                indxti, indxtj,   &
                                indxui, indxuj,   &
                                zetaD,            &
                                Cb,     vrel,     &
                                aiu,    umassdti, & 
                                fld2)

      use ice_blocks, only: nx_block, ny_block
      use ice_boundary, only: ice_HaloUpdate
      use ice_domain, only: halo_info
      use ice_domain_size, only: max_blocks
      use ice_flux, only:  fm
      use ice_grid, only: dxt, dyt, dxhy, dyhx, cxp, cyp, cxm, cym, &
          tarear, uarear, tinyarea
      use ice_state, only: uvel, vvel
      
      integer (kind=int_kind), intent(in) :: & 
         ntot         ! size of problem for fgmres (for given cpu)
      
      real (kind=dbl_kind), dimension (ntot), intent(in) :: &
         bvec     , & ! RHS vector for FGMRES
         diagvec      ! diagonal of matrix A for preconditioners
      
      real (kind=dbl_kind), dimension (ntot), intent(inout) :: &
         sol          ! solution vector for FGMRES
      
      integer (kind=int_kind), dimension(max_blocks), intent(in) :: & 
         icellt   , & ! no. of cells where icetmask = 1
         icellu       ! no. of cells where iceumask = 1
      
      integer (kind=int_kind), dimension (nx_block*ny_block, max_blocks), intent(in) :: &
         indxti   , & ! compressed index in i-direction
         indxtj   , & ! compressed index in j-direction
         indxui   , & ! compressed index in i-direction
         indxuj       ! compressed index in j-direction
      
      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), intent(in) :: &
         vrel     , & ! coeff for tauw 
         Cb       , & ! seabed stress coeff
         aiu      , & ! ice fraction on u-grid
         umassdti     ! mass of U-cell/dte (kg/m^2 s)
      
      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks,4), intent(in) :: &
         zetaD      ! zetaD = 2zeta (viscous coeff)
      
      real (kind=dbl_kind), dimension (nx_block,ny_block,2,max_blocks), intent(inout) :: &
         fld2        ! work array for boundary updates
      
      ! local variables

      integer (kind=int_kind) :: &
         iblk           , & ! block index
         icode          , & ! code for fgmres solver
         its            , & ! iteration nb for fgmres
         fgmres_its     , & ! final nb of fgmres iterations
         ierr               ! code for pgmres preconditioner !phb: needed?
      
      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         Au       , & ! matvec, Fx = Au - bx
         Av           ! matvec, Fy = Av - by
      
      real (kind=dbl_kind), allocatable :: & 
         vv(:,:), ww(:,:)  ! work arrays for FGMRES

      real (kind=dbl_kind), allocatable :: &
         wk11(:), wk22(:)  ! work vectors for FGMRES
         
      real (kind=dbl_kind) :: &
         res_norm          ! residual norm for FGMRES
      
      character(len=*), parameter :: subname = '(fgmres_solver)'

      ! ! Allocate space for FGMRES work arrays
      ! allocate(wk11(ntot), wk22(ntot))
      ! allocate(vv(ntot,im_fgmres+1), ww(ntot,im_fgmres))
      ! 
      ! !-----------------------------------------------------------------------
      ! !     prep F G M R E S 
      ! !-----------------------------------------------------------------------
      ! 
      ! icode  = 0
      ! 
      ! !-----------------------------------------------------------------------
      ! !     F G M R E S   L O O P
      ! !-----------------------------------------------------------------------
      ! 1    continue
      ! !-----------------------------------------------------------------------
      ! 
      ! 
      ! call fgmres (ntot,im_fgmres,bvec,sol,its,vv,ww,wk11,wk22, &
      !              gamma, maxits_fgmres,monitor_fgmres, &
      !              icode, fgmres_its, res_norm)
      ! 
      ! if (icode == 1) then
      ! 
      !    if (precond .eq. 1) then
      ! 
      !      wk22(:)=wk11(:) ! precond=identity
      ! 
      !    elseif (precond .eq. 2) then ! use diagonal of A for precond step
      ! 
      !      call precond_diag (ntot,            & 
      !                         diagvec (:),     &
      !                         wk11 (:), wk22 (:) )
      ! 
      !    elseif (precond .eq. 3) then
      ! 
      !     call pgmres (nx_block,    ny_block,    nblocks       , &
      !                  max_blocks         , icellu   (:)       , & 
      !                  indxui   (:,:)     , indxuj   (:,:)     , &
      !                  icellt   (:)       ,                      & 
      !                  indxti   (:,:)     , indxtj   (:,:)     , &
      !                  dxt      (:,:,:)   , dyt      (:,:,:)   , & 
      !                  dxhy     (:,:,:)   , dyhx     (:,:,:)   , & 
      !                  cxp      (:,:,:)   , cyp      (:,:,:)   , & 
      !                  cxm      (:,:,:)   , cym      (:,:,:)   , &
      !                  vrel     (:,:,:)   , Cb       (:,:,:)   , &  
      !                  zetaD    (:,:,:,:) ,                      &
      !                  umassdti (:,:,:)   , fm       (:,:,:)   , & 
      !                  uarear   (:,:,:)   , diagvec(:)         , &
      !                  wk22     (:)       , wk11(:)            , &
      !                  ntot               , im_pgmres          , &
      !                  epsprecond         , maxits_pgmres      , &
      !                  monitor_pgmres     , ierr )         
      !    endif ! precond
      ! 
      !    goto 1
      ! 
      ! elseif (icode >= 2) then
      ! 
      !    call vec_to_arrays (nx_block, ny_block, nblocks,      &
      !                        max_blocks, icellu (:), ntot,     & 
      !                        indxui    (:,:), indxuj(:,:),     &
      !                        wk11 (:),                         &
      !                        uvel (:,:,:), vvel (:,:,:))    
      ! 
      !    ! JFL halo update could be in subroutine...                    
      !    !$OMP PARALLEL DO PRIVATE(iblk) 
      !    do iblk = 1, nblocks                             
      !       fld2(:,:,1,iblk) = uvel(:,:,iblk)
      !       fld2(:,:,2,iblk) = vvel(:,:,iblk)            
      !    enddo
      !    !$OMP END PARALLEL DO                           
      ! 
      !    call ice_HaloUpdate (fld2,               halo_info, & 
      !                         field_loc_NEcorner, field_type_vector)
      ! 
      !    !$OMP PARALLEL DO PRIVATE(iblk)
      !    do iblk = 1, nblocks
      !       uvel(:,:,iblk) = fld2(:,:,1,iblk)
      !       vvel(:,:,iblk) = fld2(:,:,2,iblk)
      !    enddo
      !    !$OMP END PARALLEL DO                             
      ! 
      !    !$OMP PARALLEL DO PRIVATE(iblk)
      !    do iblk = 1, nblocks                                  
      ! 
      !     call matvec (nx_block             , ny_block,            &
      !                  icellu   (iblk)      , icellt   (iblk)    , &
      !                  indxui   (:,iblk)    , indxuj   (:,iblk)  , &
      !                  indxti   (:,iblk)    , indxtj   (:,iblk)  , &
      !                  dxt      (:,:,iblk)  , dyt      (:,:,iblk), & 
      !                  dxhy     (:,:,iblk)  , dyhx     (:,:,iblk), & 
      !                  cxp      (:,:,iblk)  , cyp      (:,:,iblk), & 
      !                  cxm      (:,:,iblk)  , cym      (:,:,iblk), &
      !                  uvel     (:,:,iblk)  , vvel     (:,:,iblk), &      
      !                  vrel     (:,:,iblk)  , Cb       (:,:,iblk), &  
      !                  zetaD    (:,:,iblk,:),                      &
      !                  umassdti (:,:,iblk)  , fm       (:,:,iblk), & 
      !                  uarear   (:,:,iblk)  ,                      & 
      !                  Au       (:,:,iblk)  , Av       (:,:,iblk))                         
      ! 
      !    enddo
      !    !$OMP END PARALLEL DO 
      ! 
      !    ! form wk2 from Au and Av arrays        
      !    call arrays_to_vec (nx_block, ny_block, nblocks,      &
      !                        max_blocks, icellu (:), ntot,     & 
      !                        indxui    (:,:), indxuj(:,:),     &
      !                        Au      (:,:,:), Av    (:,:,:),   &
      !                        wk22(:))    
      ! 
      !       goto 1
      ! 
      ! endif ! icode
      ! 
      ! deallocate(wk11, wk22, vv, ww)
         
      end subroutine fgmres_solver

!=======================================================================

! Computes the viscous coefficients (in fact zetaD=2*zeta) and dPr/dx. 

      subroutine calc_zeta_Pr  (nx_block,   ny_block,   & 
                                icellt,                 & 
                                indxti,     indxtj,     & 
                                uvel,       vvel,       & 
                                dxt,        dyt,        & 
                                dxhy,       dyhx,       & 
                                cxp,        cyp,        & 
                                cxm,        cym,        & 
                                tinyarea,               & 
                                strength,   zetaD,      &
                                stPr)

      use ice_dyn_shared, only: strain_rates

      integer (kind=int_kind), intent(in) :: & 
         nx_block, ny_block, & ! block dimensions
         icellt                ! no. of cells where icetmask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), & 
         intent(in) :: &
         indxti   , & ! compressed index in i-direction
         indxtj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         strength , & ! ice strength (N/m)
         uvel     , & ! x-component of velocity (m/s)
         vvel     , & ! y-component of velocity (m/s)
         dxt      , & ! width of T-cell through the middle (m)
         dyt      , & ! height of T-cell through the middle (m)
         dxhy     , & ! 0.5*(HTE - HTE)
         dyhx     , & ! 0.5*(HTN - HTN)
         cyp      , & ! 1.5*HTE - 0.5*HTE
         cxp      , & ! 1.5*HTN - 0.5*HTN
         cym      , & ! 0.5*HTE - 1.5*HTE
         cxm      , & ! 0.5*HTN - 1.5*HTN
         tinyarea     ! puny*tarea
         
      real (kind=dbl_kind), dimension(nx_block,ny_block,4), & 
         intent(out) :: &
         zetaD          ! 2*zeta
         
      real (kind=dbl_kind), dimension(nx_block,ny_block,8), & 
         intent(out) :: &
         stPr          ! stress Pr combinations

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &
        divune, divunw, divuse, divusw                , & ! divergence
        tensionne, tensionnw, tensionse, tensionsw    , & ! tension
        shearne, shearnw, shearse, shearsw            , & ! shearing
        Deltane, Deltanw, Deltase, Deltasw            , & ! Delt
        ssigpn, ssigps, ssigpe, ssigpw, ssigp1, ssigp2, &
        csigpne, csigpnw, csigpsw, csigpse            , &
        stressp_1, stressp_2, stressp_3, stressp_4    , &
        strp_tmp
        
      logical :: capping ! of the viscous coeff  

      character(len=*), parameter :: subname = '(calc_zeta_Pr)'

      capping = .false.
      
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icellt
         i = indxti(ij)
         j = indxtj(ij)

      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------
         call strain_rates (nx_block,   ny_block,   &
                            i,          j,          &
                            uvel,       vvel,       &
                            dxt,        dyt,        &
                            cxp,        cyp,        &
                            cxm,        cym,        &
                            divune,     divunw,     &
                            divuse,     divusw,     &
                            tensionne,  tensionnw,  &
                            tensionse,  tensionsw,  &
                            shearne,    shearnw,    &
                            shearse,    shearsw,    &
                            Deltane,    Deltanw,    &
                            Deltase,    Deltasw     )

         if (capping) then
         
          zetaD(i,j,1) = strength(i,j)/max(Deltane,tinyarea(i,j))
          zetaD(i,j,2) = strength(i,j)/max(Deltanw,tinyarea(i,j))
          zetaD(i,j,3) = strength(i,j)/max(Deltasw,tinyarea(i,j))
          zetaD(i,j,4) = strength(i,j)/max(Deltase,tinyarea(i,j))
          
         else

          zetaD(i,j,1) = strength(i,j)/(Deltane + tinyarea(i,j))
          zetaD(i,j,2) = strength(i,j)/(Deltanw + tinyarea(i,j))
          zetaD(i,j,3) = strength(i,j)/(Deltasw + tinyarea(i,j))
          zetaD(i,j,4) = strength(i,j)/(Deltase + tinyarea(i,j))
          
         
         endif
         
      !-----------------------------------------------------------------
      ! the stresses                            ! kg/s^2
      ! (1) northeast, (2) northwest, (3) southwest, (4) southeast
      !-----------------------------------------------------------------

         stressp_1 = -zetaD(i,j,1)*(Deltane*(c1-Ktens))
         stressp_2 = -zetaD(i,j,2)*(Deltanw*(c1-Ktens))
         stressp_3 = -zetaD(i,j,3)*(Deltasw*(c1-Ktens))
         stressp_4 = -zetaD(i,j,4)*(Deltase*(c1-Ktens))
         
      !-----------------------------------------------------------------
      ! combinations of the Pr related stresses for the momentum equation ! kg/s^2
      !-----------------------------------------------------------------

         ssigpn  = stressp_1 + stressp_2
         ssigps  = stressp_3 + stressp_4
         ssigpe  = stressp_1 + stressp_4
         ssigpw  = stressp_2 + stressp_3
         ssigp1  =(stressp_1 + stressp_3)*p055
         ssigp2  =(stressp_2 + stressp_4)*p055

         csigpne = p111*stressp_1 + ssigp2 + p027*stressp_3
         csigpnw = p111*stressp_2 + ssigp1 + p027*stressp_4
         csigpsw = p111*stressp_3 + ssigp2 + p027*stressp_1
         csigpse = p111*stressp_4 + ssigp1 + p027*stressp_2
         
      !-----------------------------------------------------------------
      ! for dF/dx (u momentum)
      !-----------------------------------------------------------------
         strp_tmp  = p25*dyt(i,j)*(p333*ssigpn  + p166*ssigps)

         ! northeast (i,j)
         stPr(i,j,1) = -strp_tmp &
              + dxhy(i,j)*(-csigpne)

         ! northwest (i+1,j)
         stPr(i,j,2) = strp_tmp  &
              + dxhy(i,j)*(-csigpnw)

         strp_tmp  = p25*dyt(i,j)*(p333*ssigps  + p166*ssigpn)

         ! southeast (i,j+1)
         stPr(i,j,3) = -strp_tmp &
              + dxhy(i,j)*(-csigpse)

         ! southwest (i+1,j+1)
         stPr(i,j,4) = strp_tmp  &
              + dxhy(i,j)*(-csigpsw)

      !-----------------------------------------------------------------
      ! for dF/dy (v momentum)
      !-----------------------------------------------------------------
         strp_tmp  = p25*dxt(i,j)*(p333*ssigpe  + p166*ssigpw)

         ! northeast (i,j)
         stPr(i,j,5) = -strp_tmp &
              - dyhx(i,j)*(csigpne)

         ! southeast (i,j+1)
         stPr(i,j,6) = strp_tmp  &
              - dyhx(i,j)*(csigpse)

         strp_tmp  = p25*dxt(i,j)*(p333*ssigpw  + p166*ssigpe)

         ! northwest (i+1,j)
         stPr(i,j,7) = -strp_tmp &
              - dyhx(i,j)*(csigpnw)

         ! southwest (i+1,j+1)
         stPr(i,j,8) = strp_tmp  &
              - dyhx(i,j)*(csigpsw)

      enddo                     ! ij

      end subroutine calc_zeta_Pr      
      
!=======================================================================

! Computes VP stress without the rep. pressure Pr (included in b vector)

      subroutine stress_prime_vpOLD (nx_block,   ny_block,   & 
                                  icellt,                 & 
                                  indxti,     indxtj,     & 
                                  uvel,       vvel,       & 
                                  dxt,        dyt,        & 
                                  dxhy,       dyhx,       & 
                                  cxp,        cyp,        & 
                                  cxm,        cym,        & 
                                  zetaD,                  & 
                                  str )

      use ice_dyn_shared, only: strain_rates

      integer (kind=int_kind), intent(in) :: & 
         nx_block, ny_block, & ! block dimensions
         icellt                ! no. of cells where icetmask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), & 
         intent(in) :: &
         indxti   , & ! compressed index in i-direction
         indxtj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         uvel     , & ! x-component of velocity (m/s)
         vvel     , & ! y-component of velocity (m/s)
         dxt      , & ! width of T-cell through the middle (m)
         dyt      , & ! height of T-cell through the middle (m)
         dxhy     , & ! 0.5*(HTE - HTE)
         dyhx     , & ! 0.5*(HTN - HTN)
         cyp      , & ! 1.5*HTE - 0.5*HTE
         cxp      , & ! 1.5*HTN - 0.5*HTN
         cym      , & ! 0.5*HTE - 1.5*HTE
         cxm          ! 0.5*HTN - 1.5*HTN
         
      real (kind=dbl_kind), dimension(nx_block,ny_block,4), & 
         intent(in) :: &
         zetaD          ! 2*zeta   

      real (kind=dbl_kind), dimension(nx_block,ny_block,8), & 
         intent(out) :: &
         str          ! stress combinations

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &
        divune, divunw, divuse, divusw            , & ! divergence
        tensionne, tensionnw, tensionse, tensionsw, & ! tension
        shearne, shearnw, shearse, shearsw        , & ! shearing
        Deltane, Deltanw, Deltase, Deltasw        , & ! Delt
        ssigpn, ssigps, ssigpe, ssigpw            , &
        ssigmn, ssigms, ssigme, ssigmw            , &
        ssig12n, ssig12s, ssig12e, ssig12w        , &
        ssigp1, ssigp2, ssigm1, ssigm2, ssig121, ssig122, &
        csigpne, csigpnw, csigpse, csigpsw        , &
        csigmne, csigmnw, csigmse, csigmsw        , &
        csig12ne, csig12nw, csig12se, csig12sw    , &
        str12ew, str12we, str12ns, str12sn        , &
        strp_tmp, strm_tmp
        
      real (kind=dbl_kind) :: &
         stressp_1, stressp_2, stressp_3, stressp_4 , & ! sigma11+sigma22 (without Pr)
         stressm_1, stressm_2, stressm_3, stressm_4 , & ! sigma11-sigma22
         stress12_1,stress12_2,stress12_3,stress12_4    ! sigma12        

      character(len=*), parameter :: subname = '(stress_prime_vpOLD)'

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      str(:,:,:) = c0

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu

      do ij = 1, icellt
         i = indxti(ij)
         j = indxtj(ij)

      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------
         call strain_rates (nx_block,   ny_block,   &
                            i,          j,          &
                            uvel,       vvel,       &
                            dxt,        dyt,        &
                            cxp,        cyp,        &
                            cxm,        cym,        &
                            divune,     divunw,     &
                            divuse,     divusw,     &
                            tensionne,  tensionnw,  &
                            tensionse,  tensionsw,  &
                            shearne,    shearnw,    &
                            shearse,    shearsw,    &
                            Deltane,    Deltanw,    &
                            Deltase,    Deltasw     )

      !-----------------------------------------------------------------
      ! the stresses                            ! kg/s^2
      ! (1) northeast, (2) northwest, (3) southwest, (4) southeast
      ! JFL commented part of stressp is for the rep pressure Pr
      !-----------------------------------------------------------------

         stressp_1 = zetaD(i,j,1)*(divune*(c1+Ktens))! - Deltane*(c1-Ktens))
         stressp_2 = zetaD(i,j,2)*(divunw*(c1+Ktens))! - Deltanw*(c1-Ktens))
         stressp_3 = zetaD(i,j,3)*(divusw*(c1+Ktens))! - Deltasw*(c1-Ktens))
         stressp_4 = zetaD(i,j,4)*(divuse*(c1+Ktens))! - Deltase*(c1-Ktens))
         
         stressm_1 = zetaD(i,j,1)*tensionne*(c1+Ktens)*ecci
         stressm_2 = zetaD(i,j,2)*tensionnw*(c1+Ktens)*ecci
         stressm_3 = zetaD(i,j,3)*tensionsw*(c1+Ktens)*ecci
         stressm_4 = zetaD(i,j,4)*tensionse*(c1+Ktens)*ecci
         
         stress12_1 = zetaD(i,j,1)*shearne*p5*(c1+Ktens)*ecci
         stress12_2 = zetaD(i,j,2)*shearnw*p5*(c1+Ktens)*ecci
         stress12_3 = zetaD(i,j,3)*shearsw*p5*(c1+Ktens)*ecci
         stress12_4 = zetaD(i,j,4)*shearse*p5*(c1+Ktens)*ecci

      !-----------------------------------------------------------------
      ! combinations of the stresses for the momentum equation ! kg/s^2
      !-----------------------------------------------------------------

         ssigpn  = stressp_1 + stressp_2
         ssigps  = stressp_3 + stressp_4
         ssigpe  = stressp_1 + stressp_4
         ssigpw  = stressp_2 + stressp_3
         ssigp1  =(stressp_1 + stressp_3)*p055
         ssigp2  =(stressp_2 + stressp_4)*p055

         ssigmn  = stressm_1 + stressm_2
         ssigms  = stressm_3 + stressm_4
         ssigme  = stressm_1 + stressm_4
         ssigmw  = stressm_2 + stressm_3
         ssigm1  =(stressm_1 + stressm_3)*p055
         ssigm2  =(stressm_2 + stressm_4)*p055

         ssig12n = stress12_1 + stress12_2
         ssig12s = stress12_3 + stress12_4
         ssig12e = stress12_1 + stress12_4
         ssig12w = stress12_2 + stress12_3
         ssig121 =(stress12_1 + stress12_3)*p111
         ssig122 =(stress12_2 + stress12_4)*p111

         csigpne = p111*stressp_1 + ssigp2 + p027*stressp_3
         csigpnw = p111*stressp_2 + ssigp1 + p027*stressp_4
         csigpsw = p111*stressp_3 + ssigp2 + p027*stressp_1
         csigpse = p111*stressp_4 + ssigp1 + p027*stressp_2
         
         csigmne = p111*stressm_1 + ssigm2 + p027*stressm_3
         csigmnw = p111*stressm_2 + ssigm1 + p027*stressm_4
         csigmsw = p111*stressm_3 + ssigm2 + p027*stressm_1
         csigmse = p111*stressm_4 + ssigm1 + p027*stressm_2
         
         csig12ne = p222*stress12_1 + ssig122 &
                  + p055*stress12_3
         csig12nw = p222*stress12_2 + ssig121 &
                  + p055*stress12_4
         csig12sw = p222*stress12_3 + ssig122 &
                  + p055*stress12_1
         csig12se = p222*stress12_4 + ssig121 &
                  + p055*stress12_2

         str12ew = p5*dxt(i,j)*(p333*ssig12e + p166*ssig12w)
         str12we = p5*dxt(i,j)*(p333*ssig12w + p166*ssig12e)
         str12ns = p5*dyt(i,j)*(p333*ssig12n + p166*ssig12s)
         str12sn = p5*dyt(i,j)*(p333*ssig12s + p166*ssig12n)

      !-----------------------------------------------------------------
      ! for dF/dx (u momentum)
      !-----------------------------------------------------------------
         strp_tmp  = p25*dyt(i,j)*(p333*ssigpn  + p166*ssigps)
         strm_tmp  = p25*dyt(i,j)*(p333*ssigmn  + p166*ssigms)

         ! northeast (i,j)
         str(i,j,1) = -strp_tmp - strm_tmp - str12ew &
              + dxhy(i,j)*(-csigpne + csigmne) + dyhx(i,j)*csig12ne

         ! northwest (i+1,j)
         str(i,j,2) = strp_tmp + strm_tmp - str12we &
              + dxhy(i,j)*(-csigpnw + csigmnw) + dyhx(i,j)*csig12nw

         strp_tmp  = p25*dyt(i,j)*(p333*ssigps  + p166*ssigpn)
         strm_tmp  = p25*dyt(i,j)*(p333*ssigms  + p166*ssigmn)

         ! southeast (i,j+1)
         str(i,j,3) = -strp_tmp - strm_tmp + str12ew &
              + dxhy(i,j)*(-csigpse + csigmse) + dyhx(i,j)*csig12se

         ! southwest (i+1,j+1)
         str(i,j,4) = strp_tmp + strm_tmp + str12we &
              + dxhy(i,j)*(-csigpsw + csigmsw) + dyhx(i,j)*csig12sw

      !-----------------------------------------------------------------
      ! for dF/dy (v momentum)
      !-----------------------------------------------------------------
         strp_tmp  = p25*dxt(i,j)*(p333*ssigpe  + p166*ssigpw)
         strm_tmp  = p25*dxt(i,j)*(p333*ssigme  + p166*ssigmw)

         ! northeast (i,j)
         str(i,j,5) = -strp_tmp + strm_tmp - str12ns &
              - dyhx(i,j)*(csigpne + csigmne) + dxhy(i,j)*csig12ne

         ! southeast (i,j+1)
         str(i,j,6) = strp_tmp - strm_tmp - str12sn &
              - dyhx(i,j)*(csigpse + csigmse) + dxhy(i,j)*csig12se

         strp_tmp  = p25*dxt(i,j)*(p333*ssigpw  + p166*ssigpe)
         strm_tmp  = p25*dxt(i,j)*(p333*ssigmw  + p166*ssigme)

         ! northwest (i+1,j)
         str(i,j,7) = -strp_tmp + strm_tmp + str12ns &
              - dyhx(i,j)*(csigpnw + csigmnw) + dxhy(i,j)*csig12nw

         ! southwest (i+1,j+1)
         str(i,j,8) = strp_tmp - strm_tmp + str12sn &
              - dyhx(i,j)*(csigpsw + csigmsw) + dxhy(i,j)*csig12sw

      enddo                     ! ij

      end subroutine stress_prime_vpOLD      
      
      
!=======================================================================

! Computes the VP stress (as diagnostic)

      subroutine stress_vp (nx_block,   ny_block,   & 
                            icellt,                 & 
                            indxti,     indxtj,     & 
                            uvel,       vvel,       & 
                            dxt,        dyt,        & 
                            dxhy,       dyhx,       & 
                            cxp,        cyp,        & 
                            cxm,        cym,        & 
                            zetaD,                  & 
                            stressp_1,  stressp_2,  & 
                            stressp_3,  stressp_4,  & 
                            stressm_1,  stressm_2,  & 
                            stressm_3,  stressm_4,  & 
                            stress12_1, stress12_2, & 
                            stress12_3, stress12_4, & 
                            str )

      use ice_dyn_shared, only: strain_rates

      integer (kind=int_kind), intent(in) :: & 
         nx_block, ny_block, & ! block dimensions
         icellt                ! no. of cells where icetmask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), & 
         intent(in) :: &
         indxti   , & ! compressed index in i-direction
         indxtj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         uvel     , & ! x-component of velocity (m/s)
         vvel     , & ! y-component of velocity (m/s)
         dxt      , & ! width of T-cell through the middle (m)
         dyt      , & ! height of T-cell through the middle (m)
         dxhy     , & ! 0.5*(HTE - HTE)
         dyhx     , & ! 0.5*(HTN - HTN)
         cyp      , & ! 1.5*HTE - 0.5*HTE
         cxp      , & ! 1.5*HTN - 0.5*HTN
         cym      , & ! 0.5*HTE - 1.5*HTE
         cxm          ! 0.5*HTN - 1.5*HTN
         
      real (kind=dbl_kind), dimension(nx_block,ny_block,4), & 
         intent(in) :: &
         zetaD          ! 2*zeta   

      real (kind=dbl_kind), dimension (nx_block,ny_block), & 
         intent(inout) :: &
         stressp_1, stressp_2, stressp_3, stressp_4 , & ! sigma11+sigma22
         stressm_1, stressm_2, stressm_3, stressm_4 , & ! sigma11-sigma22
         stress12_1,stress12_2,stress12_3,stress12_4    ! sigma12

      real (kind=dbl_kind), dimension(nx_block,ny_block,8), & 
         intent(out) :: &
         str          ! stress combinations

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &
        divune, divunw, divuse, divusw            , & ! divergence
        tensionne, tensionnw, tensionse, tensionsw, & ! tension
        shearne, shearnw, shearse, shearsw        , & ! shearing
        Deltane, Deltanw, Deltase, Deltasw        , & ! Delt
        ssigpn, ssigps, ssigpe, ssigpw            , &
        ssigmn, ssigms, ssigme, ssigmw            , &
        ssig12n, ssig12s, ssig12e, ssig12w        , &
        ssigp1, ssigp2, ssigm1, ssigm2, ssig121, ssig122, &
        csigpne, csigpnw, csigpse, csigpsw        , &
        csigmne, csigmnw, csigmse, csigmsw        , &
        csig12ne, csig12nw, csig12se, csig12sw    , &
        str12ew, str12we, str12ns, str12sn        , &
        strp_tmp, strm_tmp

        character(len=*), parameter :: subname = '(stress_vp)'

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      str(:,:,:) = c0

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu

      do ij = 1, icellt
         i = indxti(ij)
         j = indxtj(ij)

      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------
         call strain_rates (nx_block,   ny_block,   &
                            i,          j,          &
                            uvel,       vvel,       &
                            dxt,        dyt,        &
                            cxp,        cyp,        &
                            cxm,        cym,        &
                            divune,     divunw,     &
                            divuse,     divusw,     &
                            tensionne,  tensionnw,  &
                            tensionse,  tensionsw,  &
                            shearne,    shearnw,    &
                            shearse,    shearsw,    &
                            Deltane,    Deltanw,    &
                            Deltase,    Deltasw     )

      !-----------------------------------------------------------------
      ! the stresses                            ! kg/s^2
      ! (1) northeast, (2) northwest, (3) southwest, (4) southeast
      !-----------------------------------------------------------------

         stressp_1(i,j) = zetaD(i,j,1)*(divune*(c1+Ktens) - Deltane*(c1-Ktens))
         stressp_2(i,j) = zetaD(i,j,2)*(divunw*(c1+Ktens) - Deltanw*(c1-Ktens))
         stressp_3(i,j) = zetaD(i,j,3)*(divusw*(c1+Ktens) - Deltasw*(c1-Ktens))
         stressp_4(i,j) = zetaD(i,j,4)*(divuse*(c1+Ktens) - Deltase*(c1-Ktens))
         
         stressm_1(i,j) = zetaD(i,j,1)*tensionne*(c1+Ktens)*ecci
         stressm_2(i,j) = zetaD(i,j,2)*tensionnw*(c1+Ktens)*ecci
         stressm_3(i,j) = zetaD(i,j,3)*tensionsw*(c1+Ktens)*ecci
         stressm_4(i,j) = zetaD(i,j,4)*tensionse*(c1+Ktens)*ecci
         
         stress12_1(i,j) = zetaD(i,j,1)*shearne*p5*(c1+Ktens)*ecci
         stress12_2(i,j) = zetaD(i,j,2)*shearnw*p5*(c1+Ktens)*ecci
         stress12_3(i,j) = zetaD(i,j,3)*shearsw*p5*(c1+Ktens)*ecci
         stress12_4(i,j) = zetaD(i,j,4)*shearse*p5*(c1+Ktens)*ecci

      !-----------------------------------------------------------------
      ! combinations of the stresses for the momentum equation ! kg/s^2
      !-----------------------------------------------------------------

         ssigpn  = stressp_1(i,j) + stressp_2(i,j)
         ssigps  = stressp_3(i,j) + stressp_4(i,j)
         ssigpe  = stressp_1(i,j) + stressp_4(i,j)
         ssigpw  = stressp_2(i,j) + stressp_3(i,j)
         ssigp1  =(stressp_1(i,j) + stressp_3(i,j))*p055
         ssigp2  =(stressp_2(i,j) + stressp_4(i,j))*p055

         ssigmn  = stressm_1(i,j) + stressm_2(i,j)
         ssigms  = stressm_3(i,j) + stressm_4(i,j)
         ssigme  = stressm_1(i,j) + stressm_4(i,j)
         ssigmw  = stressm_2(i,j) + stressm_3(i,j)
         ssigm1  =(stressm_1(i,j) + stressm_3(i,j))*p055
         ssigm2  =(stressm_2(i,j) + stressm_4(i,j))*p055

         ssig12n = stress12_1(i,j) + stress12_2(i,j)
         ssig12s = stress12_3(i,j) + stress12_4(i,j)
         ssig12e = stress12_1(i,j) + stress12_4(i,j)
         ssig12w = stress12_2(i,j) + stress12_3(i,j)
         ssig121 =(stress12_1(i,j) + stress12_3(i,j))*p111
         ssig122 =(stress12_2(i,j) + stress12_4(i,j))*p111

         csigpne = p111*stressp_1(i,j) + ssigp2 + p027*stressp_3(i,j)
         csigpnw = p111*stressp_2(i,j) + ssigp1 + p027*stressp_4(i,j)
         csigpsw = p111*stressp_3(i,j) + ssigp2 + p027*stressp_1(i,j)
         csigpse = p111*stressp_4(i,j) + ssigp1 + p027*stressp_2(i,j)
         
         csigmne = p111*stressm_1(i,j) + ssigm2 + p027*stressm_3(i,j)
         csigmnw = p111*stressm_2(i,j) + ssigm1 + p027*stressm_4(i,j)
         csigmsw = p111*stressm_3(i,j) + ssigm2 + p027*stressm_1(i,j)
         csigmse = p111*stressm_4(i,j) + ssigm1 + p027*stressm_2(i,j)
         
         csig12ne = p222*stress12_1(i,j) + ssig122 &
                  + p055*stress12_3(i,j)
         csig12nw = p222*stress12_2(i,j) + ssig121 &
                  + p055*stress12_4(i,j)
         csig12sw = p222*stress12_3(i,j) + ssig122 &
                  + p055*stress12_1(i,j)
         csig12se = p222*stress12_4(i,j) + ssig121 &
                  + p055*stress12_2(i,j)

         str12ew = p5*dxt(i,j)*(p333*ssig12e + p166*ssig12w)
         str12we = p5*dxt(i,j)*(p333*ssig12w + p166*ssig12e)
         str12ns = p5*dyt(i,j)*(p333*ssig12n + p166*ssig12s)
         str12sn = p5*dyt(i,j)*(p333*ssig12s + p166*ssig12n)

      !-----------------------------------------------------------------
      ! for dF/dx (u momentum)
      !-----------------------------------------------------------------
         strp_tmp  = p25*dyt(i,j)*(p333*ssigpn  + p166*ssigps)
         strm_tmp  = p25*dyt(i,j)*(p333*ssigmn  + p166*ssigms)

         ! northeast (i,j)
         str(i,j,1) = -strp_tmp - strm_tmp - str12ew &
              + dxhy(i,j)*(-csigpne + csigmne) + dyhx(i,j)*csig12ne

         ! northwest (i+1,j)
         str(i,j,2) = strp_tmp + strm_tmp - str12we &
              + dxhy(i,j)*(-csigpnw + csigmnw) + dyhx(i,j)*csig12nw

         strp_tmp  = p25*dyt(i,j)*(p333*ssigps  + p166*ssigpn)
         strm_tmp  = p25*dyt(i,j)*(p333*ssigms  + p166*ssigmn)

         ! southeast (i,j+1)
         str(i,j,3) = -strp_tmp - strm_tmp + str12ew &
              + dxhy(i,j)*(-csigpse + csigmse) + dyhx(i,j)*csig12se

         ! southwest (i+1,j+1)
         str(i,j,4) = strp_tmp + strm_tmp + str12we &
              + dxhy(i,j)*(-csigpsw + csigmsw) + dyhx(i,j)*csig12sw

      !-----------------------------------------------------------------
      ! for dF/dy (v momentum)
      !-----------------------------------------------------------------
         strp_tmp  = p25*dxt(i,j)*(p333*ssigpe  + p166*ssigpw)
         strm_tmp  = p25*dxt(i,j)*(p333*ssigme  + p166*ssigmw)

         ! northeast (i,j)
         str(i,j,5) = -strp_tmp + strm_tmp - str12ns &
              - dyhx(i,j)*(csigpne + csigmne) + dxhy(i,j)*csig12ne

         ! southeast (i,j+1)
         str(i,j,6) = strp_tmp - strm_tmp - str12sn &
              - dyhx(i,j)*(csigpse + csigmse) + dxhy(i,j)*csig12se

         strp_tmp  = p25*dxt(i,j)*(p333*ssigpw  + p166*ssigpe)
         strm_tmp  = p25*dxt(i,j)*(p333*ssigmw  + p166*ssigme)

         ! northwest (i+1,j)
         str(i,j,7) = -strp_tmp + strm_tmp + str12ns &
              - dyhx(i,j)*(csigpnw + csigmnw) + dxhy(i,j)*csig12nw

         ! southwest (i+1,j+1)
         str(i,j,8) = strp_tmp - strm_tmp + str12sn &
              - dyhx(i,j)*(csigpsw + csigmsw) + dxhy(i,j)*csig12sw

      enddo                     ! ij

      end subroutine stress_vp
      
!=======================================================================

! calc deformations for mechanical redistribution

      subroutine deformations (nx_block,   ny_block,   & 
                               icellt,                 & 
                               indxti,     indxtj,     & 
                               uvel,       vvel,       & 
                               dxt,        dyt,        & 
                               cxp,        cyp,        & 
                               cxm,        cym,        & 
                               tarear,                 & 
                               shear,      divu,       & 
                               rdg_conv,   rdg_shear )

      integer (kind=int_kind), intent(in) :: & 
         nx_block, ny_block, & ! block dimensions
         icellt                ! no. of cells where icetmask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), & 
         intent(in) :: &
         indxti   , & ! compressed index in i-direction
         indxtj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         uvel     , & ! x-component of velocity (m/s)
         vvel     , & ! y-component of velocity (m/s)
         dxt      , & ! width of T-cell through the middle (m)
         dyt      , & ! height of T-cell through the middle (m)
         cyp      , & ! 1.5*HTE - 0.5*HTE
         cxp      , & ! 1.5*HTN - 0.5*HTN
         cym      , & ! 0.5*HTE - 1.5*HTE
         cxm      , & ! 0.5*HTN - 1.5*HTN
         tarear       ! 1/tarea
         
      real (kind=dbl_kind), dimension (nx_block,ny_block), & 
         intent(inout) :: &
         shear    , & ! strain rate II component (1/s)
         divu     , & ! strain rate I component, velocity divergence (1/s)
         rdg_conv , & ! convergence term for ridging (1/s)
         rdg_shear    ! shear term for ridging (1/s)

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &
        divune, divunw, divuse, divusw            , & ! divergence
        tensionne, tensionnw, tensionse, tensionsw, & ! tension
        shearne, shearnw, shearse, shearsw        , & ! shearing
        Deltane, Deltanw, Deltase, Deltasw        , & ! Delt
        tmp
        
      character(len=*), parameter :: subname = '(deformations)'
      
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu

      do ij = 1, icellt
         i = indxti(ij)
         j = indxtj(ij)

      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------
         ! divergence  =  e_11 + e_22
         divune    = cyp(i,j)*uvel(i  ,j  ) - dyt(i,j)*uvel(i-1,j  ) &
                   + cxp(i,j)*vvel(i  ,j  ) - dxt(i,j)*vvel(i  ,j-1)
         divunw    = cym(i,j)*uvel(i-1,j  ) + dyt(i,j)*uvel(i  ,j  ) &
                   + cxp(i,j)*vvel(i-1,j  ) - dxt(i,j)*vvel(i-1,j-1)
         divusw    = cym(i,j)*uvel(i-1,j-1) + dyt(i,j)*uvel(i  ,j-1) &
                   + cxm(i,j)*vvel(i-1,j-1) + dxt(i,j)*vvel(i-1,j  )
         divuse    = cyp(i,j)*uvel(i  ,j-1) - dyt(i,j)*uvel(i-1,j-1) &
                   + cxm(i,j)*vvel(i  ,j-1) + dxt(i,j)*vvel(i  ,j  )

         ! tension strain rate  =  e_11 - e_22
         tensionne = -cym(i,j)*uvel(i  ,j  ) - dyt(i,j)*uvel(i-1,j  ) &
                   +  cxm(i,j)*vvel(i  ,j  ) + dxt(i,j)*vvel(i  ,j-1)
         tensionnw = -cyp(i,j)*uvel(i-1,j  ) + dyt(i,j)*uvel(i  ,j  ) &
                   +  cxm(i,j)*vvel(i-1,j  ) + dxt(i,j)*vvel(i-1,j-1)
         tensionsw = -cyp(i,j)*uvel(i-1,j-1) + dyt(i,j)*uvel(i  ,j-1) &
                   +  cxp(i,j)*vvel(i-1,j-1) - dxt(i,j)*vvel(i-1,j  )
         tensionse = -cym(i,j)*uvel(i  ,j-1) - dyt(i,j)*uvel(i-1,j-1) &
                   +  cxp(i,j)*vvel(i  ,j-1) - dxt(i,j)*vvel(i  ,j  )

         ! shearing strain rate  =  e_12
         shearne = -cym(i,j)*vvel(i  ,j  ) - dyt(i,j)*vvel(i-1,j  ) &
                 -  cxm(i,j)*uvel(i  ,j  ) - dxt(i,j)*uvel(i  ,j-1)
         shearnw = -cyp(i,j)*vvel(i-1,j  ) + dyt(i,j)*vvel(i  ,j  ) &
                 -  cxm(i,j)*uvel(i-1,j  ) - dxt(i,j)*uvel(i-1,j-1)
         shearsw = -cyp(i,j)*vvel(i-1,j-1) + dyt(i,j)*vvel(i  ,j-1) &
                 -  cxp(i,j)*uvel(i-1,j-1) + dxt(i,j)*uvel(i-1,j  )
         shearse = -cym(i,j)*vvel(i  ,j-1) - dyt(i,j)*vvel(i-1,j-1) &
                 -  cxp(i,j)*uvel(i  ,j-1) + dxt(i,j)*uvel(i  ,j  )
         
         ! Delta (in the denominator of zeta, eta)
         Deltane = sqrt(divune**2 + ecci*(tensionne**2 + shearne**2))
         Deltanw = sqrt(divunw**2 + ecci*(tensionnw**2 + shearnw**2))
         Deltasw = sqrt(divusw**2 + ecci*(tensionsw**2 + shearsw**2))
         Deltase = sqrt(divuse**2 + ecci*(tensionse**2 + shearse**2))

      !-----------------------------------------------------------------
      ! deformations for mechanical redistribution
      !-----------------------------------------------------------------

         divu(i,j) = p25*(divune + divunw + divuse + divusw) * tarear(i,j)
            tmp = p25*(Deltane + Deltanw + Deltase + Deltasw)   * tarear(i,j)
            rdg_conv(i,j)  = -min(divu(i,j),c0)
            rdg_shear(i,j) = p5*(tmp-abs(divu(i,j))) 

            ! diagnostic only
            ! shear = sqrt(tension**2 + shearing**2)
         shear(i,j) = p25*tarear(i,j)*sqrt( &
                   (tensionne + tensionnw + tensionse + tensionsw)**2 &
                +  (shearne +   shearnw +   shearse +   shearsw)**2)

      enddo                     ! ij

      end subroutine deformations            
      
!=======================================================================

      subroutine calc_vrel_Cb (nx_block,   ny_block, &
                               icellu,     Cw,       &
                               indxui,     indxuj,   &
                               aiu,        Tbu,      &
                               uocn,       vocn,     &
                               uvel,       vvel,     &
                               vrel,       Cb)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellu                ! total count when iceumask is true

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxui  , & ! compressed index in i-direction
         indxuj      ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         Tbu,      & ! coefficient for basal stress (N/m^2)
         aiu     , & ! ice fraction on u-grid
         uocn    , & ! ocean current, x-direction (m/s)
         vocn        ! ocean current, y-direction (m/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         uvel    , & ! x-component of velocity (m/s)
         vvel        ! y-component of velocity (m/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         vrel      , & ! coeff for tauw ! jfl
         Cb            ! seabed stress coeff ! jfl
         
      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         Cw                   ! ocean-ice neutral drag coefficient

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &
         utp, vtp          , & ! utp = uvel, vtp = vvel !jfl needed?
         rhow                  !

      real (kind=dbl_kind) :: &
         u0 = 5e-5_dbl_kind    ! residual velocity for basal stress (m/s)
         
      character(len=*), parameter :: subname = '(calc_vrel_Cb)'
      
      !-----------------------------------------------------------------
      ! integrate the momentum equation
      !-----------------------------------------------------------------

      call icepack_query_parameters(rhow_out=rhow)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message="subname", &
         file=__FILE__, line=__LINE__)

      do ij =1, icellu
         i = indxui(ij)
         j = indxuj(ij)

         utp = uvel(i,j)
         vtp = vvel(i,j)

         ! (magnitude of relative ocean current)*rhow*drag*aice
         vrel(i,j) = aiu(i,j)*rhow*Cw(i,j)*sqrt((uocn(i,j) - utp)**2 + &
                                           (vocn(i,j) - vtp)**2)  ! m/s
      
         Cb(i,j)  = Tbu(i,j) / (sqrt(utp**2 + vtp**2) + u0) ! for basal stress

      enddo                     ! ij

      end subroutine calc_vrel_Cb      
      
!=======================================================================

      subroutine matvecOLD (nx_block,   ny_block, &
                         icellu,               &
                         indxui,     indxuj,   &
                         str,                  &
                         vrel,                 &
                         umassdti,   fm,       &
                         uarear,     Cb,       &
                         uvel,       vvel,     &
                         Au,         Av)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellu                ! total count when iceumask is true

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxui  , & ! compressed index in i-direction
         indxuj      ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         vrel,     & ! coefficient for tauw
         Cb,       & ! coefficient for basal stress
         umassdti, & ! mass of U-cell/dt (kg/m^2 s)
         fm      , & ! Coriolis param. * mass in U-cell (kg/s)
         uarear      ! 1/uarea

      real (kind=dbl_kind), dimension(nx_block,ny_block,8), &
         intent(in) :: &
         str

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         uvel    , & ! x-component of velocity (m/s)
         vvel        ! y-component of velocity (m/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         Au      , & ! matvec, Fx = Au - bx (N/m^2)! jfl
         Av          ! matvec, Fy = Av - by (N/m^2)! jfl    

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &
         utp, vtp          , & ! utp = uvel, vtp = vvel
         ccaimp,ccb        , & ! intermediate variables
         strintx, strinty

      character(len=*), parameter :: subname = '(matvecOLD)'

      !-----------------------------------------------------------------
      ! integrate the momentum equation
      !-----------------------------------------------------------------

      do ij =1, icellu
         i = indxui(ij)
         j = indxuj(ij)

         utp = uvel(i,j)
         vtp = vvel(i,j)

         ccaimp = umassdti(i,j) + vrel(i,j) * cosw + Cb(i,j) ! kg/m^2 s
               
         ccb = fm(i,j) + sign(c1,fm(i,j)) * vrel(i,j) * sinw ! kg/m^2 s

         ! divergence of the internal stress tensor
         strintx = uarear(i,j)* &
             (str(i,j,1) + str(i+1,j,2) + str(i,j+1,3) + str(i+1,j+1,4))
         strinty = uarear(i,j)* &
             (str(i,j,5) + str(i,j+1,6) + str(i+1,j,7) + str(i+1,j+1,8))

         Au(i,j) = ccaimp*utp - ccb*vtp - strintx
         Av(i,j) = ccaimp*vtp + ccb*utp - strinty
         
      enddo                     ! ij

      end subroutine matvecOLD
      
!=======================================================================

      subroutine matvec (nx_block,   ny_block, &
                         icellu,     icellt ,  &
                         indxui,     indxuj,   &
                         indxti,     indxtj,   &
                         dxt,        dyt,      & 
                         dxhy,       dyhx,     & 
                         cxp,        cyp,      & 
                         cxm,        cym,      & 
                         uvel,       vvel,     &
                         vrel,       Cb,       &
                         zetaD,                & 
                         umassdti,   fm,       &
                         uarear,               &
                         Au,         Av)

      use ice_dyn_shared, only: strain_rates

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellu,             & ! total count when iceumask is true
         icellt                ! no. of cells where icetmask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxui  , & ! compressed index in i-direction
         indxuj  , & ! compressed index in j-direction
         indxti  , & ! compressed index in i-direction
         indxtj      ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         dxt      , & ! width of T-cell through the middle (m)
         dyt      , & ! height of T-cell through the middle (m)
         dxhy     , & ! 0.5*(HTE - HTE)
         dyhx     , & ! 0.5*(HTN - HTN)
         cyp      , & ! 1.5*HTE - 0.5*HTE
         cxp      , & ! 1.5*HTN - 0.5*HTN
         cym      , & ! 0.5*HTE - 1.5*HTE
         cxm          ! 0.5*HTN - 1.5*HTN
         
      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         uvel    , & ! x-component of velocity (m/s)
         vvel    , & ! y-component of velocity (m/s)
         vrel    , & ! coefficient for tauw
         Cb      , & ! coefficient for basal stress
         umassdti, & ! mass of U-cell/dt (kg/m^2 s)
         fm      , & ! Coriolis param. * mass in U-cell (kg/s)
         uarear      ! 1/uarea

      real (kind=dbl_kind), dimension(nx_block,ny_block,4), & 
         intent(in) :: &
         zetaD          ! 2*zeta   
         
      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(out) :: &
         Au      , & ! matvec, Fx = Au - bx (N/m^2)! jfl
         Av          ! matvec, Fy = Av - by (N/m^2)! jfl    

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind), dimension(nx_block,ny_block,8):: &
         str    
         
      real (kind=dbl_kind) :: &
         utp, vtp          , & ! utp = uvel, vtp = vvel
         ccaimp,ccb        , & ! intermediate variables
         strintx, strinty

      real (kind=dbl_kind) :: &
        divune, divunw, divuse, divusw            , & ! divergence
        tensionne, tensionnw, tensionse, tensionsw, & ! tension
        shearne, shearnw, shearse, shearsw        , & ! shearing
        Deltane, Deltanw, Deltase, Deltasw        , & ! Delt
        ssigpn, ssigps, ssigpe, ssigpw            , &
        ssigmn, ssigms, ssigme, ssigmw            , &
        ssig12n, ssig12s, ssig12e, ssig12w        , &
        ssigp1, ssigp2, ssigm1, ssigm2, ssig121, ssig122, &
        csigpne, csigpnw, csigpse, csigpsw        , &
        csigmne, csigmnw, csigmse, csigmsw        , &
        csig12ne, csig12nw, csig12se, csig12sw    , &
        str12ew, str12we, str12ns, str12sn        , &
        strp_tmp, strm_tmp
        
      real (kind=dbl_kind) :: &
         stressp_1, stressp_2, stressp_3, stressp_4 , & ! sigma11+sigma22 (without Pr)
         stressm_1, stressm_2, stressm_3, stressm_4 , & ! sigma11-sigma22
         stress12_1,stress12_2,stress12_3,stress12_4    ! sigma12        

      character(len=*), parameter :: subname = '(matvec)'

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      str(:,:,:) = c0

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu

      do ij = 1, icellt
         i = indxti(ij)
         j = indxtj(ij)

      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------
         call strain_rates (nx_block,   ny_block,   &
                            i,          j,          &
                            uvel,       vvel,       &
                            dxt,        dyt,        &
                            cxp,        cyp,        &
                            cxm,        cym,        &
                            divune,     divunw,     &
                            divuse,     divusw,     &
                            tensionne,  tensionnw,  &
                            tensionse,  tensionsw,  &
                            shearne,    shearnw,    &
                            shearse,    shearsw,    &
                            Deltane,    Deltanw,    &
                            Deltase,    Deltasw     )

      !-----------------------------------------------------------------
      ! the stresses                            ! kg/s^2
      ! (1) northeast, (2) northwest, (3) southwest, (4) southeast
      ! JFL commented part of stressp is for the rep pressure Pr
      !-----------------------------------------------------------------

         stressp_1 = zetaD(i,j,1)*(divune*(c1+Ktens))! - Deltane*(c1-Ktens))
         stressp_2 = zetaD(i,j,2)*(divunw*(c1+Ktens))! - Deltanw*(c1-Ktens))
         stressp_3 = zetaD(i,j,3)*(divusw*(c1+Ktens))! - Deltasw*(c1-Ktens))
         stressp_4 = zetaD(i,j,4)*(divuse*(c1+Ktens))! - Deltase*(c1-Ktens))
         
         stressm_1 = zetaD(i,j,1)*tensionne*(c1+Ktens)*ecci
         stressm_2 = zetaD(i,j,2)*tensionnw*(c1+Ktens)*ecci
         stressm_3 = zetaD(i,j,3)*tensionsw*(c1+Ktens)*ecci
         stressm_4 = zetaD(i,j,4)*tensionse*(c1+Ktens)*ecci
         
         stress12_1 = zetaD(i,j,1)*shearne*p5*(c1+Ktens)*ecci
         stress12_2 = zetaD(i,j,2)*shearnw*p5*(c1+Ktens)*ecci
         stress12_3 = zetaD(i,j,3)*shearsw*p5*(c1+Ktens)*ecci
         stress12_4 = zetaD(i,j,4)*shearse*p5*(c1+Ktens)*ecci

      !-----------------------------------------------------------------
      ! combinations of the stresses for the momentum equation ! kg/s^2
      !-----------------------------------------------------------------

         ssigpn  = stressp_1 + stressp_2
         ssigps  = stressp_3 + stressp_4
         ssigpe  = stressp_1 + stressp_4
         ssigpw  = stressp_2 + stressp_3
         ssigp1  =(stressp_1 + stressp_3)*p055
         ssigp2  =(stressp_2 + stressp_4)*p055

         ssigmn  = stressm_1 + stressm_2
         ssigms  = stressm_3 + stressm_4
         ssigme  = stressm_1 + stressm_4
         ssigmw  = stressm_2 + stressm_3
         ssigm1  =(stressm_1 + stressm_3)*p055
         ssigm2  =(stressm_2 + stressm_4)*p055

         ssig12n = stress12_1 + stress12_2
         ssig12s = stress12_3 + stress12_4
         ssig12e = stress12_1 + stress12_4
         ssig12w = stress12_2 + stress12_3
         ssig121 =(stress12_1 + stress12_3)*p111
         ssig122 =(stress12_2 + stress12_4)*p111

         csigpne = p111*stressp_1 + ssigp2 + p027*stressp_3
         csigpnw = p111*stressp_2 + ssigp1 + p027*stressp_4
         csigpsw = p111*stressp_3 + ssigp2 + p027*stressp_1
         csigpse = p111*stressp_4 + ssigp1 + p027*stressp_2
         
         csigmne = p111*stressm_1 + ssigm2 + p027*stressm_3
         csigmnw = p111*stressm_2 + ssigm1 + p027*stressm_4
         csigmsw = p111*stressm_3 + ssigm2 + p027*stressm_1
         csigmse = p111*stressm_4 + ssigm1 + p027*stressm_2
         
         csig12ne = p222*stress12_1 + ssig122 &
                  + p055*stress12_3
         csig12nw = p222*stress12_2 + ssig121 &
                  + p055*stress12_4
         csig12sw = p222*stress12_3 + ssig122 &
                  + p055*stress12_1
         csig12se = p222*stress12_4 + ssig121 &
                  + p055*stress12_2

         str12ew = p5*dxt(i,j)*(p333*ssig12e + p166*ssig12w)
         str12we = p5*dxt(i,j)*(p333*ssig12w + p166*ssig12e)
         str12ns = p5*dyt(i,j)*(p333*ssig12n + p166*ssig12s)
         str12sn = p5*dyt(i,j)*(p333*ssig12s + p166*ssig12n)

      !-----------------------------------------------------------------
      ! for dF/dx (u momentum)
      !-----------------------------------------------------------------
         strp_tmp  = p25*dyt(i,j)*(p333*ssigpn  + p166*ssigps)
         strm_tmp  = p25*dyt(i,j)*(p333*ssigmn  + p166*ssigms)

         ! northeast (i,j)
         str(i,j,1) = -strp_tmp - strm_tmp - str12ew &
              + dxhy(i,j)*(-csigpne + csigmne) + dyhx(i,j)*csig12ne

         ! northwest (i+1,j)
         str(i,j,2) = strp_tmp + strm_tmp - str12we &
              + dxhy(i,j)*(-csigpnw + csigmnw) + dyhx(i,j)*csig12nw

         strp_tmp  = p25*dyt(i,j)*(p333*ssigps  + p166*ssigpn)
         strm_tmp  = p25*dyt(i,j)*(p333*ssigms  + p166*ssigmn)

         ! southeast (i,j+1)
         str(i,j,3) = -strp_tmp - strm_tmp + str12ew &
              + dxhy(i,j)*(-csigpse + csigmse) + dyhx(i,j)*csig12se

         ! southwest (i+1,j+1)
         str(i,j,4) = strp_tmp + strm_tmp + str12we &
              + dxhy(i,j)*(-csigpsw + csigmsw) + dyhx(i,j)*csig12sw

      !-----------------------------------------------------------------
      ! for dF/dy (v momentum)
      !-----------------------------------------------------------------
         strp_tmp  = p25*dxt(i,j)*(p333*ssigpe  + p166*ssigpw)
         strm_tmp  = p25*dxt(i,j)*(p333*ssigme  + p166*ssigmw)

         ! northeast (i,j)
         str(i,j,5) = -strp_tmp + strm_tmp - str12ns &
              - dyhx(i,j)*(csigpne + csigmne) + dxhy(i,j)*csig12ne

         ! southeast (i,j+1)
         str(i,j,6) = strp_tmp - strm_tmp - str12sn &
              - dyhx(i,j)*(csigpse + csigmse) + dxhy(i,j)*csig12se

         strp_tmp  = p25*dxt(i,j)*(p333*ssigpw  + p166*ssigpe)
         strm_tmp  = p25*dxt(i,j)*(p333*ssigmw  + p166*ssigme)

         ! northwest (i+1,j)
         str(i,j,7) = -strp_tmp + strm_tmp + str12ns &
              - dyhx(i,j)*(csigpnw + csigmnw) + dxhy(i,j)*csig12nw

         ! southwest (i+1,j+1)
         str(i,j,8) = strp_tmp - strm_tmp + str12sn &
              - dyhx(i,j)*(csigpsw + csigmsw) + dxhy(i,j)*csig12sw
  
      enddo ! ij - icellt
         
      !-----------------------------------------------------------------
      ! Form Au and Av
      !-----------------------------------------------------------------

      do ij =1, icellu
         i = indxui(ij)
         j = indxuj(ij)

         utp = uvel(i,j)
         vtp = vvel(i,j)

         ccaimp = umassdti(i,j) + vrel(i,j) * cosw + Cb(i,j) ! kg/m^2 s
               
         ccb = fm(i,j) + sign(c1,fm(i,j)) * vrel(i,j) * sinw ! kg/m^2 s

         ! divergence of the internal stress tensor
         strintx = uarear(i,j)* &
             (str(i,j,1) + str(i+1,j,2) + str(i,j+1,3) + str(i+1,j+1,4))
         strinty = uarear(i,j)* &
             (str(i,j,5) + str(i,j+1,6) + str(i+1,j,7) + str(i+1,j+1,8))

         Au(i,j) = ccaimp*utp - ccb*vtp - strintx
         Av(i,j) = ccaimp*vtp + ccb*utp - strinty
         
      enddo ! ij - icellu               

      end subroutine matvec

!=======================================================================      
      
     subroutine calc_bfix  (nx_block,   ny_block,   & 
                            icellu,                 & 
                            indxui,     indxuj,     & 
                            umassdti,               & 
                            forcex,     forcey,     &
                            uvel_init,  vvel_init,  &
                            bxfix,      byfix)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellu                ! no. of cells where iceumask = 1
         
      integer (kind=int_kind), dimension (nx_block*ny_block), & 
         intent(in) :: &
         indxui   , & ! compressed index in i-direction
         indxuj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), & 
         intent(in) :: &
         uvel_init,& ! x-component of velocity (m/s), beginning of time step
         vvel_init,& ! y-component of velocity (m/s), beginning of time step
         umassdti, & ! mass of U-cell/dt (kg/m^2 s)
         forcex  , & ! work array: combined atm stress and ocn tilt, x
         forcey      ! work array: combined atm stress and ocn tilt, y

      real (kind=dbl_kind), dimension (nx_block,ny_block), & 
         intent(out) :: &
         bxfix   , & ! bx = taux + bxfix !jfl
         byfix       ! by = tauy + byfix !jfl

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      character(len=*), parameter :: subname = '(calc_bfix)'

      !-----------------------------------------------------------------
      ! Define variables for momentum equation
      !-----------------------------------------------------------------

      do ij = 1, icellu
         i = indxui(ij)
         j = indxuj(ij)

         bxfix(i,j) = umassdti(i,j)*uvel_init(i,j) + forcex(i,j)
         byfix(i,j) = umassdti(i,j)*vvel_init(i,j) + forcey(i,j)
         
      enddo

      end subroutine calc_bfix            
      
!=======================================================================

      subroutine calc_bvec (nx_block,   ny_block, &
                       icellu,               &
                       indxui,     indxuj,   &
                       stPr,       Cw,       &
                       aiu,        uarear,   &
                       uocn,       vocn,     &
                       waterx,     watery,   &
                       uvel,       vvel,     &
                       bxfix,      byfix,    &
                       bx,         by,       &
                       vrel)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellu                ! total count when iceumask is true

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxui  , & ! compressed index in i-direction
         indxuj      ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         uvel    , & ! x-component of velocity (m/s)
         vvel    , & ! y-component of velocity (m/s)
         Cw      , & ! ocean-ice neutral drag coefficient
         aiu     , & ! ice fraction on u-grid
         uarear  , & ! 1/uarea
         waterx  , & ! for ocean stress calculation, x (m/s)
         watery  , & ! for ocean stress calculation, y (m/s)
         bxfix   , & ! bx = taux + bxfix !jfl
         byfix   , & ! by = tauy + byfix !jfl
         uocn    , & ! ocean current, x-direction (m/s)
         vocn    , & ! ocean current, y-direction (m/s)
         vrel        ! relative ice-ocean velocity
         
      real (kind=dbl_kind), dimension(nx_block,ny_block,8), &
         intent(in) :: &
         stPr
         
      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         bx      , & ! b vector, bx = taux + bxfix (N/m^2) !jfl
         by          ! b vector, by = tauy + byfix (N/m^2) !jfl
         
      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &
         utp, vtp          , & ! utp = uvel, vtp = vvel !jfl needed?
         taux, tauy        , & ! part of ocean stress term
         strintx, strinty  , & ! divergence of the internal stress tensor (only Pr part)
         rhow                  !
         
      character(len=*), parameter :: subname = '(calc_bvec)'
      
      !-----------------------------------------------------------------
      ! calc b vector
      !-----------------------------------------------------------------
      
      call icepack_query_parameters(rhow_out=rhow)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message="subname", &
         file=__FILE__, line=__LINE__)

      do ij =1, icellu
         i = indxui(ij)
         j = indxuj(ij)

         utp = uvel(i,j)
         vtp = vvel(i,j)

         ! ice/ocean stress
         taux = vrel(i,j)*waterx(i,j) ! NOTE this is not the entire
         tauy = vrel(i,j)*watery(i,j) ! ocn stress term
         
         ! divergence of the internal stress tensor (only Pr part, i.e. dPr/dx)
         strintx = uarear(i,j)* &
             (stPr(i,j,1) + stPr(i+1,j,2) + stPr(i,j+1,3) + stPr(i+1,j+1,4))
         strinty = uarear(i,j)* &
             (stPr(i,j,5) + stPr(i,j+1,6) + stPr(i+1,j,7) + stPr(i+1,j+1,8))

         bx(i,j) = bxfix(i,j) + taux + strintx
         by(i,j) = byfix(i,j) + tauy + strinty
         
      enddo                     ! ij

      end subroutine calc_bvec
      
      !=======================================================================

      subroutine residual_vec (nx_block,   ny_block, &
                               icellu,               &
                               indxui,     indxuj,   &
                               bx,         by,       &
                               Au,         Av,       &
                               Fx,         Fy,       &
                               sum_squared)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellu                ! total count when iceumask is true

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxui  , & ! compressed index in i-direction
         indxuj      ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         bx       , & ! b vector, bx = taux + bxfix (N/m^2) !jfl
         by       , & ! b vector, by = tauy + byfix (N/m^2) !jfl
         Au       , & ! matvec, Fx = Au - bx (N/m^2) ! jfl
         Av           ! matvec, Fy = Av - by (N/m^2) ! jfl

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         Fx      , & ! x residual vector, Fx = Au - bx (N/m^2)
         Fy          ! y residual vector, Fy = Av - by (N/m^2)
         
      real (kind=dbl_kind), intent(out), optional :: &
         sum_squared ! sum of squared residual vector components
      
      integer (kind=int_kind) :: &
         i, j, ij

      character(len=*), parameter :: subname = '(residual_vec)'

      !-----------------------------------------------------------------
      ! compute residual and sum its squared components
      !-----------------------------------------------------------------

      if (present(sum_squared)) then
         sum_squared = c0
      endif
         
      do ij =1, icellu
         i = indxui(ij)
         j = indxuj(ij)

         Fx(i,j) = Au(i,j) - bx(i,j)
         Fy(i,j) = Av(i,j) - by(i,j)
         if (present(sum_squared)) then
            sum_squared = sum_squared + Fx(i,j)**2 + Fy(i,j)**2
         endif
      enddo                     ! ij

      end subroutine residual_vec
      
!=======================================================================

      subroutine formDiag_step1  (nx_block,   ny_block,   & 
                                  icellu,                 & 
                                  indxui,     indxuj,     & 
                                  dxt,        dyt,        & 
                                  dxhy,       dyhx,       & 
                                  cxp,        cyp,        & 
                                  cxm,        cym,        & 
                                  zetaD,      Dstr )

      integer (kind=int_kind), intent(in) :: & 
         nx_block, ny_block, & ! block dimensions
         icellu                ! no. of cells where icetmask = 1 JFL

      integer (kind=int_kind), dimension (nx_block*ny_block), & 
         intent(in) :: &
         indxui   , & ! compressed index in i-direction JFL
         indxuj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         dxt      , & ! width of T-cell through the middle (m)
         dyt      , & ! height of T-cell through the middle (m)
         dxhy     , & ! 0.5*(HTE - HTE)
         dyhx     , & ! 0.5*(HTN - HTN)
         cyp      , & ! 1.5*HTE - 0.5*HTE
         cxp      , & ! 1.5*HTN - 0.5*HTN
         cym      , & ! 0.5*HTE - 1.5*HTE
         cxm          ! 0.5*HTN - 1.5*HTN
         
      real (kind=dbl_kind), dimension(nx_block,ny_block,4), & 
         intent(in) :: &
         zetaD          ! 2*zeta      
         
      real (kind=dbl_kind), dimension(nx_block,ny_block,8), & 
         intent(inout) :: &
         Dstr          ! intermediate calc for diagonal components of matrix A associated 
                       ! with rheology term         

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij, iu, ju, di, dj, cc

      real (kind=dbl_kind) :: &
        divune, divunw, divuse, divusw            , & ! divergence
        tensionne, tensionnw, tensionse, tensionsw, & ! tension
        shearne, shearnw, shearse, shearsw        , & ! shearing
        uij,ui1j,uij1,ui1j1,vij,vi1j,vij1,vi1j1   , & ! c0 or c1
        stressp_1, stressp_2, stressp_3, stressp_4, &
        stressm_1, stressm_2, stressm_3, stressm_4, &
        stress12_1,stress12_2,stress12_3,stress12_4,&
        ssigpn, ssigps, ssigpe, ssigpw            , &
        ssigmn, ssigms, ssigme, ssigmw            , &
        ssig12n, ssig12s, ssig12e, ssig12w        , &
        ssigp1, ssigp2, ssigm1, ssigm2, ssig121, ssig122, &
        csigpne, csigpnw, csigpse, csigpsw        , &
        csigmne, csigmnw, csigmse, csigmsw        , &
        csig12ne, csig12nw, csig12se, csig12sw    , &
        str12ew, str12we, str12ns, str12sn        , &
        strp_tmp, strm_tmp

      character(len=*), parameter :: subname = '(formDiag_step1)'

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu

      Dstr(:,:,:) = c0 ! BE careful: Dstr contains 4 terms for u and 4 terms for v. These 8 
                       ! come from the surrounding T cells but are all refrerenced to the i,j (u point)  
       
 !      Dstr(i,j,1) corresponds to str(i,j,1)
 !      Dstr(i,j,2) corresponds to str(i+1,j,2)
 !      Dstr(i,j,3) corresponds to str(i,j+1,3) 
 !      Dstr(i,j,4) corresponds to str(i+1,j+1,4))
 !      Dstr(i,j,5) corresponds to str(i,j,5) 
 !      Dstr(i,j,6) corresponds to str(i,j+1,6)
 !      Dstr(i,j,7) corresponds to str(i+1,j,7) 
 !      Dstr(i,j,8) corresponds to str(i+1,j+1,8))
             
      do cc=1, 8 ! 4 for u and 4 for v
      
       if (cc .eq. 1) then     ! u comp, T cell i,j
        uij   = c1
        ui1j  = c0
        uij1  = c0
        ui1j1 = c0
        vij   = c0
        vi1j  = c0
        vij1  = c0
        vi1j1 = c0
        di    = 0
        dj    = 0
       elseif (cc .eq. 2) then ! u comp, T cell i+1,j
        uij   = c0
        ui1j  = c1
        uij1  = c0
        ui1j1 = c0
        vij   = c0
        vi1j  = c0
        vij1  = c0
        vi1j1 = c0
        di    = 1
        dj    = 0
       elseif (cc .eq. 3) then ! u comp, T cell i,j+1
        uij   = c0
        ui1j  = c0
        uij1  = c1
        ui1j1 = c0
        vij   = c0
        vi1j  = c0
        vij1  = c0
        vi1j1 = c0
        di    = 0
        dj    = 1
       elseif (cc .eq. 4) then ! u comp, T cell i+1,j+1
        uij   = c0
        ui1j  = c0
        uij1  = c0
        ui1j1 = c1
        vij   = c0
        vi1j  = c0
        vij1  = c0
        vi1j1 = c0
        di    = 1
        dj    = 1
       elseif (cc .eq. 5) then ! v comp, T cell i,j
        uij   = c0
        ui1j  = c0
        uij1  = c0
        ui1j1 = c0
        vij   = c1
        vi1j  = c0
        vij1  = c0
        vi1j1 = c0
        di    = 0
        dj    = 0
       elseif (cc .eq. 6) then ! v comp, T cell i,j+1
        uij   = c0
        ui1j  = c0
        uij1  = c0
        ui1j1 = c0
        vij   = c0
        vi1j  = c0
        vij1  = c1
        vi1j1 = c0
        di    = 0
        dj    = 1
       elseif (cc .eq. 7) then ! v comp, T cell i+1,j
        uij   = c0
        ui1j  = c0
        uij1  = c0
        ui1j1 = c0
        vij   = c0
        vi1j  = c1
        vij1  = c0
        vi1j1 = c0
        di    = 1
        dj    = 0
       elseif (cc .eq. 8) then ! v comp, T cell i+1,j+1
        uij   = c0
        ui1j  = c0
        uij1  = c0
        ui1j1 = c0
        vij   = c0
        vi1j  = c0
        vij1  = c0
        vi1j1 = c1
        di    = 1
        dj    = 1
       endif 
       
      do ij = 1, icellu
      
         iu = indxui(ij)
         ju = indxuj(ij)
         i=iu+di
         j=ju+dj
          
      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------
         ! divergence  =  e_11 + e_22
         divune    = cyp(i,j)*uij - dyt(i,j)*ui1j &
                   + cxp(i,j)*vij - dxt(i,j)*vij1
         divunw    = cym(i,j)*ui1j + dyt(i,j)*uij &
                   + cxp(i,j)*vi1j - dxt(i,j)*vi1j1
         divusw    = cym(i,j)*ui1j1 + dyt(i,j)*uij1 &
                   + cxm(i,j)*vi1j1 + dxt(i,j)*vi1j
         divuse    = cyp(i,j)*uij1 - dyt(i,j)*ui1j1 &
                   + cxm(i,j)*vij1 + dxt(i,j)*vij

         ! tension strain rate  =  e_11 - e_22
         tensionne = -cym(i,j)*uij - dyt(i,j)*ui1j &
                   +  cxm(i,j)*vij + dxt(i,j)*vij1
         tensionnw = -cyp(i,j)*ui1j + dyt(i,j)*uij &
                   +  cxm(i,j)*vi1j + dxt(i,j)*vi1j1
         tensionsw = -cyp(i,j)*ui1j1 + dyt(i,j)*uij1 &
                   +  cxp(i,j)*vi1j1 - dxt(i,j)*vi1j
         tensionse = -cym(i,j)*uij1  - dyt(i,j)*ui1j1 &
                   +  cxp(i,j)*vij1 - dxt(i,j)*vij

         ! shearing strain rate  =  e_12
         shearne = -cym(i,j)*vij - dyt(i,j)*vi1j &
                 -  cxm(i,j)*uij - dxt(i,j)*uij1
         shearnw = -cyp(i,j)*vi1j + dyt(i,j)*vij &
                 -  cxm(i,j)*ui1j - dxt(i,j)*ui1j1
         shearsw = -cyp(i,j)*vi1j1 + dyt(i,j)*vij1 &
                 -  cxp(i,j)*ui1j1 + dxt(i,j)*ui1j
         shearse = -cym(i,j)*vij1 - dyt(i,j)*vi1j1 &
                 -  cxp(i,j)*uij1 + dxt(i,j)*uij
         
      !-----------------------------------------------------------------
      ! the stresses                            ! kg/s^2
      ! (1) northeast, (2) northwest, (3) southwest, (4) southeast
      !-----------------------------------------------------------------
         
         stressp_1 = zetaD(i,j,1)*divune*(c1+Ktens)
         stressp_2 = zetaD(i,j,2)*divunw*(c1+Ktens)
         stressp_3 = zetaD(i,j,3)*divusw*(c1+Ktens)
         stressp_4 = zetaD(i,j,4)*divuse*(c1+Ktens)
         
         stressm_1 = zetaD(i,j,1)*tensionne*(c1+Ktens)*ecci
         stressm_2 = zetaD(i,j,2)*tensionnw*(c1+Ktens)*ecci
         stressm_3 = zetaD(i,j,3)*tensionsw*(c1+Ktens)*ecci
         stressm_4 = zetaD(i,j,4)*tensionse*(c1+Ktens)*ecci
                          
         stress12_1 = zetaD(i,j,1)*shearne*p5*(c1+Ktens)*ecci
         stress12_2 = zetaD(i,j,2)*shearnw*p5*(c1+Ktens)*ecci
         stress12_3 = zetaD(i,j,3)*shearsw*p5*(c1+Ktens)*ecci
         stress12_4 = zetaD(i,j,4)*shearse*p5*(c1+Ktens)*ecci

      !-----------------------------------------------------------------
      ! combinations of the stresses for the momentum equation ! kg/s^2
      !-----------------------------------------------------------------

         ssigpn  = stressp_1 + stressp_2
         ssigps  = stressp_3 + stressp_4
         ssigpe  = stressp_1 + stressp_4
         ssigpw  = stressp_2 + stressp_3
         ssigp1  =(stressp_1 + stressp_3)*p055
         ssigp2  =(stressp_2 + stressp_4)*p055

         ssigmn  = stressm_1 + stressm_2
         ssigms  = stressm_3 + stressm_4
         ssigme  = stressm_1 + stressm_4
         ssigmw  = stressm_2 + stressm_3
         ssigm1  =(stressm_1 + stressm_3)*p055
         ssigm2  =(stressm_2 + stressm_4)*p055

         ssig12n = stress12_1 + stress12_2
         ssig12s = stress12_3 + stress12_4
         ssig12e = stress12_1 + stress12_4
         ssig12w = stress12_2 + stress12_3
         ssig121 =(stress12_1 + stress12_3)*p111
         ssig122 =(stress12_2 + stress12_4)*p111

         csigpne = p111*stressp_1 + ssigp2 + p027*stressp_3
         csigpnw = p111*stressp_2 + ssigp1 + p027*stressp_4
         csigpsw = p111*stressp_3 + ssigp2 + p027*stressp_1
         csigpse = p111*stressp_4 + ssigp1 + p027*stressp_2
         
         csigmne = p111*stressm_1 + ssigm2 + p027*stressm_3
         csigmnw = p111*stressm_2 + ssigm1 + p027*stressm_4
         csigmsw = p111*stressm_3 + ssigm2 + p027*stressm_1
         csigmse = p111*stressm_4 + ssigm1 + p027*stressm_2
         
         csig12ne = p222*stress12_1 + ssig122 &
                  + p055*stress12_3
         csig12nw = p222*stress12_2 + ssig121 &
                  + p055*stress12_4
         csig12sw = p222*stress12_3 + ssig122 &
                  + p055*stress12_1
         csig12se = p222*stress12_4 + ssig121 &
                  + p055*stress12_2

         str12ew = p5*dxt(i,j)*(p333*ssig12e + p166*ssig12w)
         str12we = p5*dxt(i,j)*(p333*ssig12w + p166*ssig12e)
         str12ns = p5*dyt(i,j)*(p333*ssig12n + p166*ssig12s)
         str12sn = p5*dyt(i,j)*(p333*ssig12s + p166*ssig12n)

      !-----------------------------------------------------------------
      ! for dF/dx (u momentum)
      !-----------------------------------------------------------------
      
         if (cc .eq. 1) then ! T cell i,j
         
          strp_tmp  = p25*dyt(i,j)*(p333*ssigpn  + p166*ssigps)
          strm_tmp  = p25*dyt(i,j)*(p333*ssigmn  + p166*ssigms)

         ! northeast (i,j)
         Dstr(iu,ju,1) = -strp_tmp - strm_tmp - str12ew &
              + dxhy(i,j)*(-csigpne + csigmne) + dyhx(i,j)*csig12ne
              
         elseif (cc .eq. 2) then ! T cell i+1,j     
          
          strp_tmp  = p25*dyt(i,j)*(p333*ssigpn  + p166*ssigps)
          strm_tmp  = p25*dyt(i,j)*(p333*ssigmn  + p166*ssigms)
          
         ! northwest (i+1,j)
         Dstr(iu,ju,2) = strp_tmp + strm_tmp - str12we &
              + dxhy(i,j)*(-csigpnw + csigmnw) + dyhx(i,j)*csig12nw

         elseif (cc .eq. 3) then ! T cell i,j+1        
              
         strp_tmp  = p25*dyt(i,j)*(p333*ssigps  + p166*ssigpn)
         strm_tmp  = p25*dyt(i,j)*(p333*ssigms  + p166*ssigmn)

         ! southeast (i,j+1)
         Dstr(iu,ju,3) = -strp_tmp - strm_tmp + str12ew &
              + dxhy(i,j)*(-csigpse + csigmse) + dyhx(i,j)*csig12se

         elseif (cc .eq. 4) then ! T cell i+1,j+1     
              
         strp_tmp  = p25*dyt(i,j)*(p333*ssigps  + p166*ssigpn)
         strm_tmp  = p25*dyt(i,j)*(p333*ssigms  + p166*ssigmn)     
              
         ! southwest (i+1,j+1)
         Dstr(iu,ju,4) = strp_tmp + strm_tmp + str12we &
              + dxhy(i,j)*(-csigpsw + csigmsw) + dyhx(i,j)*csig12sw

      !-----------------------------------------------------------------
      ! for dF/dy (v momentum)
      !-----------------------------------------------------------------
         
         elseif (cc .eq. 5) then ! T cell i,j 
         
         strp_tmp  = p25*dxt(i,j)*(p333*ssigpe  + p166*ssigpw)
         strm_tmp  = p25*dxt(i,j)*(p333*ssigme  + p166*ssigmw)

         ! northeast (i,j)
         Dstr(iu,ju,5) = -strp_tmp + strm_tmp - str12ns &
              - dyhx(i,j)*(csigpne + csigmne) + dxhy(i,j)*csig12ne

         elseif (cc .eq. 6) then ! T cell i,j+1      
              
         strp_tmp  = p25*dxt(i,j)*(p333*ssigpe  + p166*ssigpw)
         strm_tmp  = p25*dxt(i,j)*(p333*ssigme  + p166*ssigmw)              
              
         ! southeast (i,j+1)
         Dstr(iu,ju,6) = strp_tmp - strm_tmp - str12sn &
              - dyhx(i,j)*(csigpse + csigmse) + dxhy(i,j)*csig12se

         elseif (cc .eq. 7) then ! T cell i,j+1     
              
         strp_tmp  = p25*dxt(i,j)*(p333*ssigpw  + p166*ssigpe)
         strm_tmp  = p25*dxt(i,j)*(p333*ssigmw  + p166*ssigme)

         ! northwest (i+1,j)
         Dstr(iu,ju,7) = -strp_tmp + strm_tmp + str12ns &
              - dyhx(i,j)*(csigpnw + csigmnw) + dxhy(i,j)*csig12nw

         elseif (cc .eq. 8) then ! T cell i+1,j+1        
              
         strp_tmp  = p25*dxt(i,j)*(p333*ssigpw  + p166*ssigpe)
         strm_tmp  = p25*dxt(i,j)*(p333*ssigmw  + p166*ssigme)     
              
         ! southwest (i+1,j+1)
         Dstr(iu,ju,8) = strp_tmp - strm_tmp + str12sn &
              - dyhx(i,j)*(csigpsw + csigmsw) + dxhy(i,j)*csig12sw
              
         endif     

      enddo                     ! ij
      
      enddo                     ! cc

      end subroutine formDiag_step1      
      
      
!=======================================================================

      subroutine formDiag_step2 (nx_block,   ny_block, &
                                 icellu,               &
                                 indxui,     indxuj,   &
                                 Dstr,       vrel,     &
                                 umassdti,             &
                                 uarear,     Cb,       &
                                 Diagu,      Diagv )

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellu                ! total count when iceumask is true

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxui  , & ! compressed index in i-direction
         indxuj      ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         vrel,     & ! coefficient for tauw
         Cb,       & ! coefficient for basal stress
         umassdti, & ! mass of U-cell/dt (kg/m^2 s)
         uarear      ! 1/uarea

      real (kind=dbl_kind), dimension(nx_block,ny_block,8), &
         intent(in) :: &
         Dstr

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         Diagu   , & ! matvec, Fx = Au - bx (N/m^2)
         Diagv       ! matvec, Fy = Av - by (N/m^2)

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &
         ccaimp             , & ! intermediate variables
         strintx, strinty

      character(len=*), parameter :: subname = '(formDiag_step2)'

      !-----------------------------------------------------------------
      ! integrate the momentum equation
      !-----------------------------------------------------------------

      strintx=c0
      strinty=c0
         
! BE careful: Dstr contains 4 terms for u and 4 terms for v. These 8 
! come from the surrounding T cells but are all refrerenced to the i,j (u point)  
       
 !      Dstr(i,j,1) corresponds to str(i,j,1)
 !      Dstr(i,j,2) corresponds to str(i+1,j,2)
 !      Dstr(i,j,3) corresponds to str(i,j+1,3) 
 !      Dstr(i,j,4) corresponds to str(i+1,j+1,4))
 !      Dstr(i,j,5) corresponds to str(i,j,5) 
 !      Dstr(i,j,6) corresponds to str(i,j+1,6)
 !      Dstr(i,j,7) corresponds to str(i+1,j,7) 
 !      Dstr(i,j,8) corresponds to str(i+1,j+1,8))
         
      do ij =1, icellu
         i = indxui(ij)
         j = indxuj(ij)

         ccaimp = umassdti(i,j) + vrel(i,j) * cosw + Cb(i,j) ! kg/m^2 s
         
         strintx = uarear(i,j)* &
             (Dstr(i,j,1) + Dstr(i,j,2) + Dstr(i,j,3) + Dstr(i,j,4))
         strinty = uarear(i,j)* &
             (Dstr(i,j,5) + Dstr(i,j,6) + Dstr(i,j,7) + Dstr(i,j,8))

         Diagu(i,j) = ccaimp - strintx
         Diagv(i,j) = ccaimp - strinty
         
      enddo                     ! ij

      end subroutine formDiag_step2     
      
!=======================================================================
      
     subroutine precond_diag   (ntot,    &
                                diagvec, &
                                wk1, wk2)

      integer (kind=int_kind), intent(in) :: &
         ntot                  ! size of problem for fgmres

      real (kind=dbl_kind), dimension (ntot), intent(in) :: &
         diagvec, wk1
         
      real (kind=dbl_kind), dimension (ntot), intent(out) :: &
         wk2 

      ! local variables

      integer (kind=int_kind) :: &
         i         

      character(len=*), parameter :: subname = '(precond_diag)'

      !-----------------------------------------------------------------
      ! form vector (converts from max_blocks arrays to single vector
      !-----------------------------------------------------------------

      wk2(:)=c0
      
      do i=1, ntot

	wk2(i) = wk1(i)/diagvec(i)
      
      enddo! i

      end subroutine precond_diag
      
!=======================================================================      

      subroutine calc_L2norm_squared (nx_block,   ny_block, &
                                      icellu,               &
                                      indxui,     indxuj,   &
                                      tpu,        tpv,      &
                                      L2norm)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellu                ! total count when iceumask is true

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxui  , & ! compressed index in i-direction
         indxuj      ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         tpu     , & ! x-component of vector grid function
         tpv         ! y-component of vector grid function

      real (kind=dbl_kind), intent(out) :: &
         L2norm      ! squared l^2 norm of vector grid function (tpu,tpv)

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      character(len=*), parameter :: subname = '(calc_L2norm_squared)'

      !-----------------------------------------------------------------
      ! compute squared l^2 norm of vector grid function (tpu,tpv)
      !-----------------------------------------------------------------

     L2norm = c0
      
      do ij =1, icellu
         i = indxui(ij)
         j = indxuj(ij)
         
         L2norm = L2norm + tpu(i,j)**2 + tpv(i,j)**2
         
      enddo ! ij
      
      end subroutine calc_L2norm_squared
      
            !=======================================================================

      subroutine arrays_to_vec (nx_block, ny_block, nblocks, &
                                max_blocks, icellu,   ntot,  &
                                indxui,   indxuj,            &
                                tpu,      tpv,               &
                                outvec)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         nblocks,            & ! nb of blocks
         max_blocks,         & ! max nb of blocks
         ntot                  ! size of problem for fgmres

      integer (kind=int_kind), dimension (max_blocks), intent(in) :: &
         icellu          
         
      integer (kind=int_kind), dimension (nx_block*ny_block, max_blocks), &
         intent(in) :: &
         indxui  , & ! compressed index in i-direction
         indxuj      ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block, max_blocks), intent(in) :: &
         tpu     , & ! x-component of vector
         tpv         ! y-component of vector         
         
      real (kind=dbl_kind), dimension (ntot), intent(inout) :: &
         outvec

      ! local variables

      integer (kind=int_kind) :: &
         i, j, iblk, tot, ij
         

      character(len=*), parameter :: subname = '(arrays_to_vec)'

      !-----------------------------------------------------------------
      ! form vector (converts from max_blocks arrays to single vector
      !-----------------------------------------------------------------

      outvec(:)=c0
      tot=0
      
      do iblk=1, nblocks
       do ij =1, icellu(iblk)
          i = indxui(ij,iblk)
          j = indxuj(ij,iblk)
          tot=tot+1
          outvec(tot)=tpu(i,j,iblk)
          tot=tot+1
          outvec(tot)=tpv(i,j,iblk)
       enddo
      enddo! ij

      end subroutine arrays_to_vec
      
            !=======================================================================

      subroutine vec_to_arrays (nx_block, ny_block, nblocks, &
                                max_blocks, icellu,   ntot,  &
                                indxui,   indxuj,            &
                                invec,                       &
                                tpu,      tpv)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         nblocks,            & ! nb of blocks
         max_blocks,         & ! max nb of blocks
         ntot                  ! size of problem for fgmres

      integer (kind=int_kind), dimension (max_blocks), intent(in) :: &
         icellu          
         
      integer (kind=int_kind), dimension (nx_block*ny_block, max_blocks), &
         intent(in) :: &
         indxui  , & ! compressed index in i-direction
         indxuj      ! compressed index in j-direction
         
      real (kind=dbl_kind), dimension (ntot), intent(in) :: &
         invec         

      real (kind=dbl_kind), dimension (nx_block,ny_block, max_blocks), intent(inout) :: &
         tpu     , & ! x-component of vector
         tpv         ! y-component of vector         
         
      ! local variables

      integer (kind=int_kind) :: &
         i, j, iblk, tot, ij
         

      character(len=*), parameter :: subname = '(vec_to_arrays)'

      !-----------------------------------------------------------------
      ! form arrays (converts from vector to the max_blocks arrays
      !-----------------------------------------------------------------

      tpu(:,:,:)=c0
      tpv(:,:,:)=c0
      tot=0
      
      do iblk=1, nblocks
       do ij =1, icellu(iblk)
          i = indxui(ij,iblk)
          j = indxuj(ij,iblk)
          tot=tot+1
          tpu(i,j,iblk)=invec(tot)
          tot=tot+1
          tpv(i,j,iblk)=invec(tot)
       enddo
      enddo! ij

      end subroutine vec_to_arrays
      
!      JFL ROUTINE POUR CALC STRESS OCN POUR COUPLAGE
      
!=======================================================================

! Update Q and R factor after deletion of the 1st column of G_diff
!
! author: P. Blain ECCC
      subroutine qr_delete(Q, R)
      
      real (kind=dbl_kind), intent(inout) :: &
         Q(:,:),  & ! Q factor
         R(:,:)     ! R factor
      
      ! local variables
      
      integer (kind=int_kind) :: &
         i, j, k, & ! loop indices
         m, n       ! size of Q matrix
         
      real (kind=dbl_kind) :: &
         temp, c, s
               
      character(len=*), parameter :: subname = '(qr_delete)'
      
      n = size(Q,1)
      m = size(Q,2)
      do i = 1, m-1
         temp = sqrt(R(i,i+1)**2 + R(i+1,i+1)**2)
         c = R(i,i+1)/temp
         s = R(i+1,i+1)/temp
         R(i,i+1) = temp
         R(i+1,i+1) = 0
         if (i < m-1) then
            do j = i+2, m
               temp = c*R(i,j) + s*R(i+1,j)
               R(i+1,j) = -s*R(i,j) + c*R(i+1,j)
               R(i,j) = temp
            enddo
         endif
         do k = 1, n
            temp = c*Q(k,i) + s*Q(k,i+1);
            Q(k,i+1) = -s*Q(k,i) + c*Q(k,i+1);
            Q(k,i) = temp
         enddo
      enddo
      R(:,1:m-1) = R(:,2:m)
      
      end subroutine qr_delete

!=======================================================================

! FGMRES: Flexible generalized minimum residual method (with restarts). 
! Solves A x = b using GMRES with a varying (right) preconditioner
!
! authors: Stéphane Gaudreault, Abdessamad Qaddouri, Philippe Blain, ECCC

      subroutine fgmres (zetaD,               &
                         Cb,         vrel,    &
                         umassdti,            &
                         solx,       soly,    &
                         bx,         by,      &
                         diagx,      diagy,   &
                         tolerance, maxinner, &
                         maxouter, nbiter, conv)

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks,4), intent(in) :: &
         zetaD   ! zetaD = 2*zeta (viscous coefficient)

      real (kind=dbl_kind), dimension(nx_block, ny_block, max_blocks), intent(in) :: &
         vrel  , & ! coefficient for tauw 
         Cb    , & ! seabed stress coefficient
         umassdti  ! mass of U-cell/dte (kg/m^2 s)

      real (kind=dbl_kind), dimension(nx_block, ny_block, max_blocks), intent(inout) :: &
         solx     , & ! Initial guess on input, approximate solution on output (x components)
         soly         ! Initial guess on input, approximate solution on output (y components)

      real (kind=dbl_kind), dimension(nx_block, ny_block, max_blocks), intent(in) :: &
         bx       , & ! Right hand side of the linear system (x components)
         by           ! Right hand side of the linear system (y components)

      real (kind=dbl_kind), dimension(nx_block, ny_block, max_blocks), intent(in) :: &
         diagx    , & ! Diagonal of the system matrix (x components)
         diagy        ! Diagonal of the system matrix (y components)

      real (kind=dbl_kind), intent(in) :: &
         tolerance   ! Tolerance to achieve. The algorithm terminates when the relative
                     ! residual is below tolerance

      integer (kind=int_kind), intent(in) :: &
         maxinner    ! Restart the method every maxinner inner iterations

      integer (kind=int_kind), intent(in) :: &
         maxouter    ! Maximum number of outer iterations
                     ! Iteration will stop after maxinner*maxouter steps
                     ! even if the specified tolerance has not been achieved

      integer (kind=int_kind), intent(out) :: &
         nbiter      ! Total number of iteration performed

      real (kind=dbl_kind), intent(out) :: &
         conv        ! !phb DESCRIBE IF WE KEEP

      ! local variables

      integer (kind=int_kind) :: &
         iblk    , & ! block index
         ij      , & ! index for indx[t|u][i|j]
         i, j        ! grid indices
         
      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         workspace_x , & ! work vector (x components)
         workspace_y , & ! work vector (y components)
         Fx          , & ! residual vector (x components), Fx = Au - bx (N/m^2)
         Fy              ! residual vector (y components), Fy = Av - by (N/m^2)

      real (kind=dbl_kind), dimension (max_blocks) :: &
         norm_squared   ! array to accumulate squared norm of grid function over blocks

      real (kind=dbl_kind), dimension(nx_block, ny_block, max_blocks, maxinner+1) :: &
         arnoldi_basis_x , & ! arnoldi basis (x components) !phb == vv
         arnoldi_basis_y     ! arnoldi basis (y components)

      real (kind=dbl_kind), dimension(nx_block, ny_block, max_blocks, maxinner) :: &
         wwx             , & ! !phb FIND BETTER NAME (x components)
         wwy                 ! !phb FIND BETTER NAME (y components)

      real (kind=dbl_kind) :: &
         norm_residual   , & ! current L^2 norm of residual vector
         inverse_norm    , & ! inverse of the norm of a vector
         nu, t               ! local temporary values

      integer (kind=int_kind) :: &
         initer          , & ! inner (Arnoldi) loop counter
         outiter         , & ! outer (restarts) loop counter
         nextit          , & ! nextit == initer+1
         it, k, ii, jj       ! reusable loop counters

      real (kind=dbl_kind), dimension(maxinner+1) :: &
         rot_cos         , & ! cosine elements of Givens rotations 
         rot_sin         , & ! sine elements of Givens rotations
         rhs_hess            ! right hand side vector of the Hessenberg (least squares) system

      real (kind=dbl_kind), dimension(maxinner+1, maxinner) :: &
         hessenberg        ! system matrix of the Hessenberg (least squares) system

      integer (kind=int_kind) :: &
         precond_type ! type of preconditioner

      real (kind=dbl_kind) :: relative_tolerance, r0 !phb DESCRIBE if we keep

      real (kind=dbl_kind) :: &
         local_dot         ! local value to accumulate dot product computations
         
      real (kind=dbl_kind), dimension(maxinner) :: &
         dotprod_local     ! local array to accumulate several dot product computations

      character(len=*), parameter :: subname = '(fgmres)'

      ! Here we go !

      outiter = 0
      nbiter = 0
      
      conv = c1
      
      precond_type = precond
      
      ! Residual of the initial iterate
      
      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks
         call matvec (nx_block               , ny_block,              &
                      icellu         (iblk)  , icellt         (iblk), &
                      indxui       (:,iblk)  , indxuj       (:,iblk), &
                      indxti       (:,iblk)  , indxtj       (:,iblk), &
                      dxt        (:,:,iblk)  , dyt        (:,:,iblk), &
                      dxhy       (:,:,iblk)  , dyhx       (:,:,iblk), &
                      cxp        (:,:,iblk)  , cyp        (:,:,iblk), &
                      cxm        (:,:,iblk)  , cym        (:,:,iblk), &
                      solx       (:,:,iblk)  , soly       (:,:,iblk), &
                      vrel       (:,:,iblk)  , Cb         (:,:,iblk), &
                      zetaD      (:,:,iblk,:),                        &
                      umassdti   (:,:,iblk)  , fm         (:,:,iblk), &
                      uarear     (:,:,iblk)  ,                        &
                      workspace_x(:,:,iblk)  , workspace_y(:,:,iblk))
         call residual_vec (nx_block             , ny_block             , &
                            icellu         (iblk),                        & 
                            indxui       (:,iblk), indxuj    (:,iblk)   , &
                            bx         (:,:,iblk), by      (:,:,iblk)   , &
                            workspace_x(:,:,iblk), workspace_y(:,:,iblk), &
                            arnoldi_basis_x (:,:,iblk, 1),                &
                            arnoldi_basis_y (:,:,iblk, 1))
      enddo
      !$OMP END PARALLEL DO
      
      ! Start outer (restarts) loop
      do
         ! Compute norm of initial residual
         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks
            call calc_L2norm_squared(nx_block,        ny_block,        &
                                     icellu   (iblk),                  &
                                     indxui (:,iblk), indxuj(:, iblk), &
                                     arnoldi_basis_x(:,:,iblk, 1),     &
                                     arnoldi_basis_y(:,:,iblk, 1),     &
                                     norm_squared(iblk))

         enddo
         !$OMP END PARALLEL DO 
         norm_residual = sqrt(global_sum(sum(norm_squared), distrb_info))
         
         ! Current guess is a good enough solution
         if (norm_residual < tolerance) then
            return
         end if
         
         ! Normalize the first Arnoldi vector
         inverse_norm = c1 / norm_residual
         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks
            do ij =1, icellu(iblk)
               i = indxui(ij, iblk)
               j = indxuj(ij, iblk)

               arnoldi_basis_x(i, j, iblk, 1) = arnoldi_basis_x(i, j, iblk, 1) * inverse_norm
               arnoldi_basis_y(i, j, iblk, 1) = arnoldi_basis_y(i, j, iblk, 1) * inverse_norm
            enddo ! ij
         enddo
         !$OMP END PARALLEL DO
         
         if (outiter == 0) then
            relative_tolerance = tolerance * norm_residual
            r0 = norm_residual
         end if
         
         conv = norm_residual / r0
         
         ! Initialize 1-st term of RHS of Hessenberg system
         rhs_hess(1)  = norm_residual
         rhs_hess(2:) = c0
         
         initer = 0
      
         ! Start of inner (Arnoldi) loop
         do
            
            nbiter = nbiter + 1
            initer = initer + 1
            nextit = initer + 1
            
            ! precondition the current Arnoldi vector
            call precondition(zetaD,                &
                              Cb,         vrel,     &
                              umassdti,             &
                              arnoldi_basis_x(:,:,:,initer), &
                              arnoldi_basis_y(:,:,:,initer), &
                              workspace_x , workspace_y    , &
                              precond_type, diagx, diagy)
            ! !phb DESCRIBE ww
            wwx(:,:,:,initer) = workspace_x
            wwy(:,:,:,initer) = workspace_y
            
            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks
               call matvec (nx_block               , ny_block,              &
                            icellu         (iblk)  , icellt         (iblk), &
                            indxui       (:,iblk)  , indxuj       (:,iblk), &
                            indxti       (:,iblk)  , indxtj       (:,iblk), &
                            dxt        (:,:,iblk)  , dyt        (:,:,iblk), &
                            dxhy       (:,:,iblk)  , dyhx       (:,:,iblk), &
                            cxp        (:,:,iblk)  , cyp        (:,:,iblk), &
                            cxm        (:,:,iblk)  , cym        (:,:,iblk), &
                            workspace_x(:,:,iblk)  , workspace_y(:,:,iblk), &
                            vrel       (:,:,iblk)  , Cb         (:,:,iblk), &
                            zetaD      (:,:,iblk,:),                        &
                            umassdti   (:,:,iblk)  , fm         (:,:,iblk), &
                            uarear     (:,:,iblk)  ,                        &
                            arnoldi_basis_x(:,:,iblk,nextit),               &
                            arnoldi_basis_y(:,:,iblk,nextit))
            enddo
            !$OMP END PARALLEL DO
            
            ! Classical Gram-Schmidt orthogonalisation process
            ! First loop of Gram-Schmidt (compute coefficients)
            dotprod_local = c0
            do it=1,initer
               local_dot = c0
               
               !$OMP PARALLEL DO PRIVATE(iblk)
               do iblk = 1, nblocks
                  do ij =1, icellu(iblk)
                     i = indxui(ij, iblk)
                     j = indxuj(ij, iblk)
                     
                     local_dot = local_dot + (arnoldi_basis_x(i, j, iblk, it) * arnoldi_basis_x(i, j, iblk, nextit)) + &
                                             (arnoldi_basis_y(i, j, iblk, it) * arnoldi_basis_y(i, j, iblk, nextit))
                  enddo ! ij
               enddo
               !$OMP END PARALLEL DO
               
                dotprod_local(it) = local_dot
            end do
            
            hessenberg(1:initer,initer) = global_sums(dotprod_local(1:initer), distrb_info)
            
            ! Second loop of Gram-Schmidt (orthonormalize)
            do it = 1, initer
               !$OMP PARALLEL DO PRIVATE(iblk)
               do iblk = 1, nblocks
                  do ij =1, icellu(iblk)
                     i = indxui(ij, iblk)
                     j = indxuj(ij, iblk)
                     
                     arnoldi_basis_x(i, j, iblk, nextit) = arnoldi_basis_x(i, j, iblk, nextit) &
                                                           - hessenberg(it,initer) * arnoldi_basis_x(i, j, iblk, it)
                     arnoldi_basis_y(i, j, iblk, nextit) = arnoldi_basis_y(i, j, iblk, nextit) &
                                                           - hessenberg(it,initer) * arnoldi_basis_y(i, j, iblk, it)
                  enddo ! ij
               enddo
               !$OMP END PARALLEL DO 
            end do
            
            ! Compute norm of new Arnoldi vector and update Hessenberg matrix
            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks
               call calc_L2norm_squared(nx_block      ,  ny_block        , &
                                        icellu   (iblk),                   &
                                        indxui (:,iblk), indxuj(:, iblk) , &
                                        arnoldi_basis_x(:,:,iblk, nextit), &
                                        arnoldi_basis_y(:,:,iblk, nextit), &
                                        norm_squared(iblk))
            enddo
            !$OMP END PARALLEL DO 
            hessenberg(nextit,initer) = sqrt(global_sum(sum(norm_squared), distrb_info))
            
            ! Watch out for happy breakdown
            if (.not. almost_zero( hessenberg(nextit,initer) ) ) then
               ! Normalize next Arnoldi vector
               inverse_norm = c1 / hessenberg(nextit,initer)
               !$OMP PARALLEL DO PRIVATE(iblk)
               do iblk = 1, nblocks
                  do ij =1, icellu(iblk)
                     i = indxui(ij, iblk)
                     j = indxuj(ij, iblk)
                     
                     arnoldi_basis_x(i, j, iblk, nextit) = arnoldi_basis_x(i, j, iblk, nextit)*inverse_norm
                     arnoldi_basis_y(i, j, iblk, nextit) = arnoldi_basis_y(i, j, iblk, nextit)*inverse_norm
                  enddo ! ij
               enddo
               !$OMP END PARALLEL DO 
            end if
            
            ! Apply previous Givens rotation to the last column of the Hessenberg matrix
            if (initer > 1) then
               do k = 2, initer
                  t = hessenberg(k-1, initer)
                  hessenberg(k-1, initer) =  rot_cos(k-1)*t + rot_sin(k-1)*hessenberg(k, initer)
                  hessenberg(k,   initer) = -rot_sin(k-1)*t + rot_cos(k-1)*hessenberg(k, initer)
               end do
            end if
            
            ! Compute new Givens rotation
            nu = sqrt(hessenberg(initer,initer)**2 + hessenberg(nextit,initer)**2)
            if (.not. almost_zero(nu)) then
               rot_cos(initer) = hessenberg(initer,initer) / nu
               rot_sin(initer) = hessenberg(nextit,initer) / nu
               
               rhs_hess(nextit) = -rot_sin(initer) * rhs_hess(initer)
               rhs_hess(initer) =  rot_cos(initer) * rhs_hess(initer)
               
               hessenberg(initer,initer) = rot_cos(initer) * hessenberg(initer,initer) + rot_sin(initer) * hessenberg(nextit,initer)
            end if
            
            ! Check for convergence
            norm_residual = abs(rhs_hess(nextit))
            conv = norm_residual / r0
             if ((initer >= maxinner) .or. (norm_residual <= relative_tolerance)) then
               exit
            endif
            
         end do ! end of inner (Arnoldi) loop
      
         ! At this point either the maximum number of inner iterations
         ! was reached or the absolute residual is below the scaled tolerance.

         ! Solve the (now upper triangular) system "hessenberg * sol_hess = rhs_hess"
         ! (sol_hess is stored in rhs_hess)
         rhs_hess(initer) = rhs_hess(initer) / hessenberg(initer,initer)
         do ii = 2, initer
            k  = initer - ii + 1
            t  = rhs_hess(k)
            do j = k + 1, initer
               t = t - hessenberg(k,j) * rhs_hess(j)
            end do
            rhs_hess(k) = t / hessenberg(k,k)
         end do
         
         ! Form linear combination to get new solution iterate
         do it = 1, initer
            t = rhs_hess(it)
            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks
               do ij =1, icellu(iblk)
                  i = indxui(ij, iblk)
                  j = indxuj(ij, iblk)

                  solx(i, j, iblk) = solx(i, j, iblk) + t * wwx(i, j, iblk, it)
                  soly(i, j, iblk) = soly(i, j, iblk) + t * wwy(i, j, iblk, it)
               enddo ! ij
            enddo
            !$OMP END PARALLEL DO
         end do
         
         ! Increment outer loop counter and check for convergence
         outiter = outiter + 1
         if (norm_residual <= relative_tolerance .or. outiter > maxouter) then
            return
         end if

         ! Solution is not convergent : compute residual vector and continue.
         
         ! The residual vector is computed here using (see Saad p. 177) :
         ! \begin{equation}
         !    r =  V_{m+1} * Q_m^T * (\gamma_{m+1} * e_{m+1})
         ! \end{equation}
         ! where : 
         ! $r$ is the residual
         ! $V_{m+1}$ is a matrix whose columns are the Arnoldi vectors from 1 to nextit (m+1)
         ! $Q_m$ is the product of the Givens rotation : Q_m = G_m G_{m-1} ... G_1
         ! $gamma_{m+1}$ is the last element of rhs_hess
         ! $e_{m+1})$ is the unit vector (0, 0, ..., 1)^T \in \reals^{m+1}
         
         ! Apply the Givens rotation in reverse order to g := \gamma_{m+1} * e_{m+1}, 
         ! store the result in rhs_hess
         do it = 1, initer
            jj = nextit - it + 1
            rhs_hess(jj-1) = -rot_sin(jj-1) * rhs_hess(jj) ! + rot_cos(jj-1) * g(jj-1) (== 0)
            rhs_hess(jj)   =  rot_cos(jj-1) * rhs_hess(jj) ! + rot_sin(jj-1) * g(jj-1) (== 0)
         end do
         
         ! Compute the residual by multiplying V_{m+1} and rhs_hess
         workspace_x = c0
         workspace_y = c0
         do it = 1, nextit
            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks
               do ij =1, icellu(iblk)
                  i = indxui(ij, iblk)
                  j = indxuj(ij, iblk)

                  workspace_x(i, j, iblk) = rhs_hess(it) * arnoldi_basis_x(i, j, iblk, it)
                  workspace_y(i, j, iblk) = rhs_hess(it) * arnoldi_basis_y(i, j, iblk, it)
               enddo ! ij
            enddo
            !$OMP END PARALLEL DO
            arnoldi_basis_x(:,:,:,1) = workspace_x
            arnoldi_basis_y(:,:,:,1) = workspace_y
         end do
      end do ! end of outer (restarts) loop
      
      return
      end subroutine fgmres

!=======================================================================

! PGMRES: Right-preconditioned generalized minimum residual method (with restarts). 
! Solves A x = b using GMRES with a right preconditioner
!
! authors: Stéphane Gaudreault, Abdessamad Qaddouri, Philippe Blain, ECCC

      subroutine pgmres (zetaD,                &
                         Cb,         vrel,     &
                         umassdti,             &
                         solx,       soly,     &
                         bx,         by,       &
                         diagx,      diagy,    &
                         tolerance,  maxinner, &
                         maxouter, nbiter, conv)

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks,4), intent(in) :: &
         zetaD   ! zetaD = 2*zeta (viscous coefficient)

      real (kind=dbl_kind), dimension(nx_block, ny_block, max_blocks), intent(in) :: &
         vrel  , & ! coefficient for tauw 
         Cb    , & ! seabed stress coefficient
         umassdti  ! mass of U-cell/dte (kg/m^2 s)

      real (kind=dbl_kind), dimension(nx_block, ny_block, max_blocks), intent(inout) :: &
         solx     , & ! Initial guess on input, approximate solution on output (x components)
         soly         ! Initial guess on input, approximate solution on output (y components)

      real (kind=dbl_kind), dimension(nx_block, ny_block, max_blocks), intent(in) :: &
         bx       , & ! Right hand side of the linear system (x components)
         by           ! Right hand side of the linear system (y components)

      real (kind=dbl_kind), dimension(nx_block, ny_block, max_blocks), intent(in) :: &
         diagx    , & ! Diagonal of the system matrix (x components)
         diagy        ! Diagonal of the system matrix (y components)

      real (kind=dbl_kind), intent(in) :: &
         tolerance   ! Tolerance to achieve. The algorithm terminates when the relative
                     ! residual is below tolerance

      integer (kind=int_kind), intent(in) :: &
         maxinner    ! Restart the method every maxinner inner iterations

      integer (kind=int_kind), intent(in) :: &
         maxouter    ! Maximum number of outer iterations
                     ! Iteration will stop after maxinner*maxouter steps
                     ! even if the specified tolerance has not been achieved

      integer (kind=int_kind), intent(out) :: &
         nbiter      ! Total number of iteration performed

      real (kind=dbl_kind), intent(out) :: &
         conv        ! !phb DESCRIBE IF WE KEEP

      ! local variables

      integer (kind=int_kind) :: &
         iblk    , & ! block index
         ij      , & ! index for indx[t|u][i|j]
         i, j        ! grid indices
         
      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         workspace_x , & ! work vector (x components)
         workspace_y , & ! work vector (y components)
         Fx          , & ! residual vector (x components), Fx = Au - bx (N/m^2)
         Fy              ! residual vector (y components), Fy = Av - by (N/m^2)

      real (kind=dbl_kind), dimension (max_blocks) :: &
         norm_squared   ! array to accumulate squared norm of grid function over blocks

      real (kind=dbl_kind), dimension(nx_block, ny_block, max_blocks, maxinner+1) :: &
         arnoldi_basis_x , & ! arnoldi basis (x components) !phb == vv
         arnoldi_basis_y     ! arnoldi basis (y components)

      real (kind=dbl_kind) :: &
         norm_residual   , & ! current L^2 norm of residual vector
         inverse_norm    , & ! inverse of the norm of a vector
         nu, t               ! local temporary values

      integer (kind=int_kind) :: &
         initer          , & ! inner (Arnoldi) loop counter
         outiter         , & ! outer (restarts) loop counter
         nextit          , & ! nextit == initer+1
         it, k, ii, jj       ! reusable loop counters

      real (kind=dbl_kind), dimension(maxinner+1) :: &
         rot_cos         , & ! cosine elements of Givens rotations 
         rot_sin         , & ! sine elements of Givens rotations
         rhs_hess            ! right hand side vector of the Hessenberg (least squares) system

      real (kind=dbl_kind), dimension(maxinner+1, maxinner) :: &
         hessenberg        ! system matrix of the Hessenberg (least squares) system

      integer (kind=int_kind) :: &
         precond_type ! type of preconditioner

      real (kind=dbl_kind) :: relative_tolerance, r0 !phb DESCRIBE if we keep

      real (kind=dbl_kind) :: &
         local_dot         ! local value to accumulate dot product computations
         
      real (kind=dbl_kind), dimension(maxinner) :: &
         dotprod_local     ! local array to accumulate several dot product computations

      character(len=*), parameter :: subname = '(pgmres)'
      
      ! Here we go !

      outiter = 0
      nbiter = 0
      
      conv = c1
      
      precond_type = 2 ! Jacobi preconditioner
      
      ! Residual of the initial iterate
      
      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks
         call matvec (nx_block               , ny_block,              &
                      icellu         (iblk)  , icellt         (iblk), &
                      indxui       (:,iblk)  , indxuj       (:,iblk), &
                      indxti       (:,iblk)  , indxtj       (:,iblk), &
                      dxt        (:,:,iblk)  , dyt        (:,:,iblk), &
                      dxhy       (:,:,iblk)  , dyhx       (:,:,iblk), &
                      cxp        (:,:,iblk)  , cyp        (:,:,iblk), &
                      cxm        (:,:,iblk)  , cym        (:,:,iblk), &
                      solx       (:,:,iblk)  , soly       (:,:,iblk), &
                      vrel       (:,:,iblk)  , Cb         (:,:,iblk), &
                      zetaD      (:,:,iblk,:),                        &
                      umassdti   (:,:,iblk)  , fm         (:,:,iblk), &
                      uarear     (:,:,iblk)  ,                        &
                      workspace_x(:,:,iblk)  , workspace_y(:,:,iblk))
         call residual_vec (nx_block             , ny_block             , &
                            icellu         (iblk),                        & 
                            indxui       (:,iblk), indxuj    (:,iblk)   , &
                            bx         (:,:,iblk), by      (:,:,iblk)   , &
                            workspace_x(:,:,iblk), workspace_y(:,:,iblk), &
                            arnoldi_basis_x (:,:,iblk, 1),                &
                            arnoldi_basis_y (:,:,iblk, 1))
      enddo
      !$OMP END PARALLEL DO
      
      ! Start outer (restarts) loop
      do
         ! Compute norm of initial residual
         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks
            call calc_L2norm_squared(nx_block,        ny_block,        &
                                     icellu   (iblk),                  &
                                     indxui (:,iblk), indxuj(:, iblk), &
                                     arnoldi_basis_x(:,:,iblk, 1),     &
                                     arnoldi_basis_y(:,:,iblk, 1),     &
                                     norm_squared(iblk))

         enddo
         !$OMP END PARALLEL DO 
         norm_residual = sqrt(global_sum(sum(norm_squared), distrb_info))
         
         ! Current guess is a good enough solution
         if (norm_residual < tolerance) then
            return
         end if
         
         ! Normalize the first Arnoldi vector
         inverse_norm = c1 / norm_residual
         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks
            do ij =1, icellu(iblk)
               i = indxui(ij, iblk)
               j = indxuj(ij, iblk)

               arnoldi_basis_x(i, j, iblk, 1) = arnoldi_basis_x(i, j, iblk, 1) * inverse_norm
               arnoldi_basis_y(i, j, iblk, 1) = arnoldi_basis_y(i, j, iblk, 1) * inverse_norm
            enddo ! ij
         enddo
         !$OMP END PARALLEL DO
         
         if (outiter == 0) then
            relative_tolerance = tolerance * norm_residual
            r0 = norm_residual
         end if
         
         conv = norm_residual / r0
         
         ! Initialize 1-st term of RHS of Hessenberg system
         rhs_hess(1)  = norm_residual
         rhs_hess(2:) = c0
         
         initer = 0
      
         ! Start of inner (Arnoldi) loop
         do
            
            nbiter = nbiter + 1
            initer = initer + 1
            nextit = initer + 1
            
            ! precondition the current Arnoldi vector
            call precondition(zetaD,                &
                              Cb,         vrel,     &
                              umassdti,             &
                              arnoldi_basis_x(:,:,:,initer), &
                              arnoldi_basis_y(:,:,:,initer), &
                              workspace_x , workspace_y    , &
                              precond_type, diagx, diagy)
            
            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks
               call matvec (nx_block               , ny_block,              &
                            icellu         (iblk)  , icellt         (iblk), &
                            indxui       (:,iblk)  , indxuj       (:,iblk), &
                            indxti       (:,iblk)  , indxtj       (:,iblk), &
                            dxt        (:,:,iblk)  , dyt        (:,:,iblk), &
                            dxhy       (:,:,iblk)  , dyhx       (:,:,iblk), &
                            cxp        (:,:,iblk)  , cyp        (:,:,iblk), &
                            cxm        (:,:,iblk)  , cym        (:,:,iblk), &
                            workspace_x(:,:,iblk)  , workspace_y(:,:,iblk), &
                            vrel       (:,:,iblk)  , Cb         (:,:,iblk), &
                            zetaD      (:,:,iblk,:),                        &
                            umassdti   (:,:,iblk)  , fm         (:,:,iblk), &
                            uarear     (:,:,iblk)  ,                        &
                            arnoldi_basis_x(:,:,iblk,nextit),               &
                            arnoldi_basis_y(:,:,iblk,nextit))
            enddo
            !$OMP END PARALLEL DO
            
            ! Classical Gram-Schmidt orthogonalisation process
            ! First loop of Gram-Schmidt (compute coefficients)
            dotprod_local = c0
            do it=1,initer
               local_dot = c0
               
               !$OMP PARALLEL DO PRIVATE(iblk)
               do iblk = 1, nblocks
                  do ij =1, icellu(iblk)
                     i = indxui(ij, iblk)
                     j = indxuj(ij, iblk)
                     
                     local_dot = local_dot + (arnoldi_basis_x(i, j, iblk, it) * arnoldi_basis_x(i, j, iblk, nextit)) + &
                                             (arnoldi_basis_y(i, j, iblk, it) * arnoldi_basis_y(i, j, iblk, nextit))
                  enddo ! ij
               enddo
               !$OMP END PARALLEL DO
               
                dotprod_local(it) = local_dot
            end do
            
            hessenberg(1:initer,initer) = global_sums(dotprod_local(1:initer), distrb_info)
            
            ! Second loop of Gram-Schmidt (orthonormalize)
            do it = 1, initer
               !$OMP PARALLEL DO PRIVATE(iblk)
               do iblk = 1, nblocks
                  do ij =1, icellu(iblk)
                     i = indxui(ij, iblk)
                     j = indxuj(ij, iblk)
                     
                     arnoldi_basis_x(i, j, iblk, nextit) = arnoldi_basis_x(i, j, iblk, nextit) &
                                                           - hessenberg(it,initer) * arnoldi_basis_x(i, j, iblk, it)
                     arnoldi_basis_y(i, j, iblk, nextit) = arnoldi_basis_y(i, j, iblk, nextit) &
                                                           - hessenberg(it,initer) * arnoldi_basis_y(i, j, iblk, it)
                  enddo ! ij
               enddo
               !$OMP END PARALLEL DO 
            end do
            
            ! Compute norm of new Arnoldi vector and update Hessenberg matrix
            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks
               call calc_L2norm_squared(nx_block      ,  ny_block        , &
                                        icellu   (iblk),                   &
                                        indxui (:,iblk), indxuj(:, iblk) , &
                                        arnoldi_basis_x(:,:,iblk, nextit), &
                                        arnoldi_basis_y(:,:,iblk, nextit), &
                                        norm_squared(iblk))
            enddo
            !$OMP END PARALLEL DO 
            hessenberg(nextit,initer) = sqrt(global_sum(sum(norm_squared), distrb_info))
            
            ! Watch out for happy breakdown
            if (.not. almost_zero( hessenberg(nextit,initer) ) ) then
               ! Normalize next Arnoldi vector
               inverse_norm = c1 / hessenberg(nextit,initer)
               !$OMP PARALLEL DO PRIVATE(iblk)
               do iblk = 1, nblocks
                  do ij =1, icellu(iblk)
                     i = indxui(ij, iblk)
                     j = indxuj(ij, iblk)
                     
                     arnoldi_basis_x(i, j, iblk, nextit) = arnoldi_basis_x(i, j, iblk, nextit)*inverse_norm
                     arnoldi_basis_y(i, j, iblk, nextit) = arnoldi_basis_y(i, j, iblk, nextit)*inverse_norm
                  enddo ! ij
               enddo
               !$OMP END PARALLEL DO 
            end if
            
            ! Apply previous Givens rotation to the last column of the Hessenberg matrix
            if (initer > 1) then
               do k = 2, initer
                  t = hessenberg(k-1, initer)
                  hessenberg(k-1, initer) =  rot_cos(k-1)*t + rot_sin(k-1)*hessenberg(k, initer)
                  hessenberg(k,   initer) = -rot_sin(k-1)*t + rot_cos(k-1)*hessenberg(k, initer)
               end do
            end if
            
            ! Compute new Givens rotation
            nu = sqrt(hessenberg(initer,initer)**2 + hessenberg(nextit,initer)**2)
            if (.not. almost_zero(nu)) then
               rot_cos(initer) = hessenberg(initer,initer) / nu
               rot_sin(initer) = hessenberg(nextit,initer) / nu
               
               rhs_hess(nextit) = -rot_sin(initer) * rhs_hess(initer)
               rhs_hess(initer) =  rot_cos(initer) * rhs_hess(initer)
               
               hessenberg(initer,initer) = rot_cos(initer) * hessenberg(initer,initer) + rot_sin(initer) * hessenberg(nextit,initer)
            end if
            
            ! Check for convergence
            norm_residual = abs(rhs_hess(nextit))
            conv = norm_residual / r0
             if ((initer >= maxinner) .or. (norm_residual <= relative_tolerance)) then
               exit
            endif
            
         end do ! end of inner (Arnoldi) loop
      
         ! At this point either the maximum number of inner iterations
         ! was reached or the absolute residual is below the scaled tolerance.

         ! Solve the (now upper triangular) system "hessenberg * sol_hess = rhs_hess"
         ! (sol_hess is stored in rhs_hess)
         rhs_hess(initer) = rhs_hess(initer) / hessenberg(initer,initer)
         do ii = 2, initer
            k  = initer - ii + 1
            t  = rhs_hess(k)
            do j = k + 1, initer
               t = t - hessenberg(k,j) * rhs_hess(j)
            end do
            rhs_hess(k) = t / hessenberg(k,k)
         end do
         
         ! Form linear combination to get new solution iterate
         workspace_x = c0
         workspace_y = c0
         do it = 1, initer
            t = rhs_hess(it)
            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks
               do ij =1, icellu(iblk)
                  i = indxui(ij, iblk)
                  j = indxuj(ij, iblk)

                  workspace_x(i, j, iblk) = workspace_x(i, j, iblk) + t * arnoldi_basis_x(i, j, iblk, it)
                  workspace_y(i, j, iblk) = workspace_y(i, j, iblk) + t * arnoldi_basis_y(i, j, iblk, it)
               enddo ! ij
            enddo
            !$OMP END PARALLEL DO
         end do
         
         ! Call preconditioner
         call precondition(zetaD,                &
                           Cb,         vrel,     &
                           umassdti,             &
                           workspace_x , workspace_y, &
                           workspace_x , workspace_y, &
                           precond_type, diagx, diagy)
         
         solx = solx + workspace_x
         soly = soly + workspace_y
         
         ! Increment outer loop counter and check for convergence
         outiter = outiter + 1
         if (norm_residual <= relative_tolerance .or. outiter > maxouter) then
            return
         end if

         ! Solution is not convergent : compute residual vector and continue.
         
         ! The residual vector is computed here using (see Saad p. 177) :
         ! \begin{equation}
         !    r =  V_{m+1} * Q_m^T * (\gamma_{m+1} * e_{m+1})
         ! \end{equation}
         ! where : 
         ! $r$ is the residual
         ! $V_{m+1}$ is a matrix whose columns are the Arnoldi vectors from 1 to nextit (m+1)
         ! $Q_m$ is the product of the Givens rotation : Q_m = G_m G_{m-1} ... G_1
         ! $gamma_{m+1}$ is the last element of rhs_hess
         ! $e_{m+1})$ is the unit vector (0, 0, ..., 1)^T \in \reals^{m+1}
         
         ! Apply the Givens rotation in reverse order to g := \gamma_{m+1} * e_{m+1}, 
         ! store the result in rhs_hess
         do it = 1, initer
            jj = nextit - it + 1
            rhs_hess(jj-1) = -rot_sin(jj-1) * rhs_hess(jj) ! + rot_cos(jj-1) * g(jj-1) (== 0)
            rhs_hess(jj)   =  rot_cos(jj-1) * rhs_hess(jj) ! + rot_sin(jj-1) * g(jj-1) (== 0)
         end do
         
         ! Compute the residual by multiplying V_{m+1} and rhs_hess
         workspace_x = c0
         workspace_y = c0
         do it = 1, nextit
            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks
               do ij =1, icellu(iblk)
                  i = indxui(ij, iblk)
                  j = indxuj(ij, iblk)

                  workspace_x(i, j, iblk) = rhs_hess(it) * arnoldi_basis_x(i, j, iblk, it)
                  workspace_y(i, j, iblk) = rhs_hess(it) * arnoldi_basis_y(i, j, iblk, it)
               enddo ! ij
            enddo
            !$OMP END PARALLEL DO
            arnoldi_basis_x(:,:,:,1) = workspace_x
            arnoldi_basis_y(:,:,:,1) = workspace_y
         end do
      end do ! end of outer (restarts) loop
      
      end subroutine pgmres

!=======================================================================

! Generic routine to precondition a vector
!
! authors: Philippe Blain, ECCC

      subroutine precondition(zetaD,                &
                              Cb,         vrel,     &
                              umassdti,             &
                              vx,         vy,       &
                              wx,         wy,       &
                              precond_type,         &
                              diagx,      diagy)

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks,4), intent(in) :: &
         zetaD   ! zetaD = 2*zeta (viscous coefficient)

      real (kind=dbl_kind), dimension(nx_block, ny_block, max_blocks), intent(in) :: &
         vrel  , & ! coefficient for tauw 
         Cb    , & ! seabed stress coefficient
         umassdti  ! mass of U-cell/dte (kg/m^2 s)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), intent(in) :: &
         vx       , & ! input vector (x components)
         vy           ! input vector (y components)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), intent(out) :: &
         wx       , & ! preconditionned vector (x components)
         wy           ! preconditionned vector (y components)

      integer (kind=int_kind), intent(in) :: &
         precond_type ! type of preconditioner

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), intent(in) :: &
         diagx    , & ! diagonal of the system matrix (x components)
         diagy        ! diagonal of the system matrix (y components)

      ! local variables

      integer (kind=int_kind) :: &
         iblk    , & ! block index
         ij      , & ! compressed index
         i, j        ! grid indices

      real (kind=dbl_kind) :: &
         tolerance   ! Tolerance for pgmres

      integer (kind=int_kind) :: &
         maxinner    ! Restart parameter for pgmres

      integer (kind=int_kind) :: &
         maxouter    ! Maximum number of outer iterations for pgmres

      integer (kind=int_kind) :: &
         nbiter      ! Total number of iteration pgmres performed

      real (kind=dbl_kind) :: &
         conv        ! !phb DESCRIBE IF WE KEEP for pgmres

      character(len=*), parameter :: subname = '(precondition)'

      if     (precond_type == 1) then ! identity (no preconditioner)
         wx = vx
         wy = vy
      elseif (precond_type == 2) then ! Jacobi preconditioner (diagonal)
         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks
            do ij =1, icellu(iblk)
               i = indxui(ij, iblk)
               j = indxuj(ij, iblk)

               wx(i,j,iblk) = vx(i,j,iblk) / diagx(i,j,iblk)
               wy(i,j,iblk) = vy(i,j,iblk) / diagy(i,j,iblk)
            enddo ! ij
         enddo
         !$OMP END PARALLEL DO 
      elseif (precond_type == 3) then ! PGMRES (Jacobi-preconditioned GMRES)
         ! Initialize preconditioned vector to 0 !phb try with wx = vx or vx/diagx
         wx = c0
         wy = c0
         tolerance = epsprecond
         maxinner = im_pgmres
         maxouter = maxits_pgmres
         call pgmres (zetaD,                &
                      Cb,         vrel,     &
                      umassdti,             &
                      wx,         wy,       &
                      vx,         vy,       &
                      diagx,      diagy,    &
                      tolerance,  maxinner, &
                      maxouter, nbiter, conv)
      else
         
      endif
      end subroutine precondition

!=======================================================================

! Check if value A is close to zero, up to machine precision
!
!author
!     Stéphane Gaudreault, ECCC -- June 2014
!
!revision
!     v4-80 - Gaudreault S.         - gfortran compatibility
!     2019  - Philippe Blain, ECCC  - converted to CICE standards

      logical function almost_zero(A) result(retval)

      real (kind=dbl_kind), intent(in) :: A

      ! local variables

      character(len=*), parameter :: subname = '(almost_zero)'

      integer (kind=int8_kind) :: aBit
      integer (kind=int8_kind), parameter :: two_complement = int(Z'80000000', kind=int8_kind)
      aBit = 0
      aBit = transfer(A, aBit)
      if (aBit < 0) then
         aBit = two_complement - aBit
      end if
      ! lexicographic order test with a tolerance of 1 adjacent float
      retval = (abs(aBit) <= 1)
      
      end function almost_zero

!=======================================================================

      end module ice_dyn_vp

!=======================================================================
