!=======================================================================
!
! Viscous-plastic sea ice dynamics model
! Computes ice velocity and deformation
!
! See:
!
! Lemieux, J.‐F., Tremblay, B., Thomas, S., Sedláček, J., and Mysak, L. A. (2008),
! Using the preconditioned Generalized Minimum RESidual (GMRES) method to solve
! the sea‐ice momentum equation, J. Geophys. Res., 113, C10004, doi:10.1029/2007JC004680.
!
! Hibler, W. D., and Ackley, S. F. (1983), Numerical simulation of the Weddell Sea pack ice,
! J. Geophys. Res., 88( C5), 2873– 2887, doi:10.1029/JC088iC05p02873.
!
! Y. Saad. A Flexible Inner-Outer Preconditioned GMRES Algorithm. SIAM J. Sci. Comput.,
! 14(2):461–469, 1993. URL: https://doi.org/10.1137/0914028, doi:10.1137/0914028.
!
! C. T. Kelley, Iterative Methods for Linear and Nonlinear Equations, SIAM, 1995.
! (https://www.siam.org/books/textbooks/fr16_book.pdf)
!
! Y. Saad, Iterative Methods for Sparse Linear Systems. SIAM, 2003.
! (http://www-users.cs.umn.edu/~saad/IterMethBook_2ndEd.pdf)
!
! Walker, H. F., & Ni, P. (2011). Anderson Acceleration for Fixed-Point Iterations.
! SIAM Journal on Numerical Analysis, 49(4), 1715–1735. https://doi.org/10.1137/10078356X
!
! Fang, H., & Saad, Y. (2009). Two classes of multisecant methods for nonlinear acceleration.
! Numerical Linear Algebra with Applications, 16(3), 197–221. https://doi.org/10.1002/nla.617
!
! Birken, P. (2015) Termination criteria for inexact fixed‐point schemes.
! Numer. Linear Algebra Appl., 22: 702– 716. doi: 10.1002/nla.1982.
!
! authors: JF Lemieux, ECCC, Philppe Blain, ECCC
!

      module ice_dyn_vp

      use ice_kinds_mod
      use ice_blocks, only: nx_block, ny_block
      use ice_boundary, only: ice_halo
      use ice_communicate, only: my_task, master_task, get_num_procs
      use ice_constants, only: field_loc_center, field_loc_NEcorner, &
          field_type_scalar, field_type_vector
      use ice_constants, only: c0, p027, p055, p111, p166, &
          p222, p25, p333, p5, c1
      use ice_domain, only: nblocks, distrb_info
      use ice_domain_size, only: max_blocks
      use ice_dyn_shared, only: dyn_prep1, dyn_prep2, dyn_finish, &
          ecci, cosw, sinw, fcor_blk, uvel_init,  &
          vvel_init, basal_stress_coeff, basalstress, Ktens, &
          stack_velocity_field,  unstack_velocity_field
      use ice_fileunits, only: nu_diag
      use ice_flux, only: fm
      use ice_global_reductions, only: global_sum, global_allreduce_sum
      use ice_grid, only: dxt, dyt, dxhy, dyhx, cxp, cyp, cxm, cym, uarear
      use ice_exit, only: abort_ice
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_ice_strength, icepack_query_parameters

      implicit none
      private
      public :: implicit_solver, init_vp

      ! namelist parameters

      integer (kind=int_kind), public :: &
         maxits_nonlin  , & ! max nb of iteration for nonlinear solver
         dim_fgmres     , & ! size of fgmres Krylov subspace
         dim_pgmres     , & ! size of pgmres Krylov subspace
         maxits_fgmres  , & ! max nb of iteration for fgmres
         maxits_pgmres  , & ! max nb of iteration for pgmres
         fpfunc_andacc  , & ! fixed point function for Anderson acceleration: 1: g(x) = FMGRES(A(x),b(x)), 2: g(x) = x - A(x)x + b(x)
         dim_andacc     , & ! size of Anderson minimization matrix (number of saved previous residuals)
         start_andacc       ! acceleration delay factor (acceleration starts at this iteration)

      logical (kind=log_kind), public :: &
         monitor_nonlin , & ! print nonlinear residual norm
         monitor_fgmres , & ! print fgmres residual norm
         monitor_pgmres , & ! print pgmres residual norm
         use_mean_vrel      ! use mean of previous 2 iterates to compute vrel (see Hibler and Ackley 1983)

      real (kind=dbl_kind), public :: &
         reltol_nonlin  , & ! nonlinear stopping criterion: reltol_nonlin*res(k=0)
         reltol_fgmres  , & ! fgmres stopping criterion: reltol_fgmres*res(k)
         reltol_pgmres  , & ! pgmres stopping criterion: reltol_pgmres*res(k)
         damping_andacc , & ! damping factor for Anderson acceleration
         reltol_andacc      ! relative tolerance for Anderson acceleration

      character (len=char_len), public :: &
         precond        , & ! preconditioner for fgmres: 'ident' (identity), 'diag' (diagonal), 'pgmres' (Jacobi-preconditioned GMRES)
         algo_nonlin    , & ! nonlinear algorithm: 'picard' (Picard iteration), 'anderson' (Anderson acceleration)
         ortho_type         ! type of orthogonalization for FGMRES ('cgs' or 'mgs')

      ! module variables

      integer (kind=int_kind), allocatable :: &
         icellt(:)    , & ! no. of cells where icetmask = 1
         icellu(:)        ! no. of cells where iceumask = 1

      integer (kind=int_kind), allocatable :: &
         indxti(:,:)  , & ! compressed index in i-direction
         indxtj(:,:)  , & ! compressed index in j-direction
         indxui(:,:)  , & ! compressed index in i-direction
         indxuj(:,:)      ! compressed index in j-direction

      real (kind=dbl_kind), allocatable :: & 
         fld2(:,:,:,:)    ! work array for boundary updates

!=======================================================================

      contains

!=======================================================================

! Initialize parameters and variables needed for the vp dynamics
! author: Philippe Blain, ECCC

      subroutine init_vp

      use ice_blocks, only: get_block, block
      use ice_boundary, only: ice_HaloUpdate
      use ice_constants, only: c1, &
          field_loc_center, field_type_scalar
      use ice_domain, only: blocks_ice, halo_info
      use ice_grid, only: tarea, tinyarea

      ! local variables

      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      type (block) :: &
         this_block           ! block information for current block

      real (kind=dbl_kind) :: &
         min_strain_rate = 2e-09_dbl_kind      ! used for recomputing tinyarea
      
      ! Initialize module variables
      allocate(icellt(max_blocks), icellu(max_blocks))
      allocate(indxti(nx_block*ny_block, max_blocks), &
               indxtj(nx_block*ny_block, max_blocks), &
               indxui(nx_block*ny_block, max_blocks), &
               indxuj(nx_block*ny_block, max_blocks))
      allocate(fld2(nx_block,ny_block,2,max_blocks))
      
      ! Redefine tinyarea using min_strain_rate
      
      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            tinyarea(i,j,iblk) = min_strain_rate*tarea(i,j,iblk)
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

      subroutine implicit_solver (dt)

      use ice_arrays_column, only: Cdn_ocn
      use ice_boundary, only: ice_HaloMask, ice_HaloUpdate, &
          ice_HaloDestroy, ice_HaloUpdate_stress
      use ice_blocks, only: block, get_block, nx_block, ny_block
      use ice_domain, only: blocks_ice, halo_info, maskhalo_dyn
      use ice_domain_size, only: max_blocks, ncat
      use ice_dyn_shared, only: deformations
      use ice_flux, only: rdg_conv, rdg_shear, strairxT, strairyT, &
          strairx, strairy, uocn, vocn, ss_tltx, ss_tlty, iceumask, fm, &
          strtltx, strtlty, strocnx, strocny, strintx, strinty, taubx, tauby, &
          strocnxT, strocnyT, strax, stray, &
          Tbu, hwater, &
          stressp_1, stressp_2, stressp_3, stressp_4, &
          stressm_1, stressm_2, stressm_3, stressm_4, &
          stress12_1, stress12_2, stress12_3, stress12_4
      use ice_grid, only: tmask, umask, dxt, dyt, cxp, cyp, cxm, cym, &
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
         ntot           , & ! size of problem for Anderson
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
         Cb       , & ! seabed stress coefficient
         fpresx   , & ! fixed point residual vector, x components: fx = uvel - uprev_k
         fpresy   , & ! fixed point residual vector, y components: fy = vvel - vprev_k
         aiu      , & ! ice fraction on u-grid
         umass    , & ! total mass of ice and snow (u grid)
         umassdti     ! mass of U-cell/dte (kg/m^2 s)

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks,4):: &
         zetaD        ! zetaD = 2zeta (viscous coeff)

      logical (kind=log_kind) :: calc_strair

      integer (kind=int_kind), dimension (nx_block,ny_block,max_blocks) :: &
         icetmask, &  ! ice extent mask (T-cell)
         halomask     ! generic halo mask

      type (ice_halo) :: &
         halo_info_mask !  ghost cell update info for masked halo

      type (block) :: &
         this_block           ! block information for current block

      real (kind=dbl_kind), allocatable :: &
         sol(:)          ! solution vector

      character(len=*), parameter :: subname = '(implicit_solver)'

      call ice_timer_start(timer_dynamics) ! dynamics

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------
      
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

         call calc_bfix (nx_block            , ny_block            , &
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
                                       strength(i,j,  iblk))
         enddo  ! ij

      enddo  ! iblk
      !$TCXOMP END PARALLEL DO

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (strength,           halo_info, &
                           field_loc_center,   field_type_scalar)
      ! velocities may have changed in dyn_prep2
      call stack_velocity_field(uvel, vvel, fld2)
      call ice_HaloUpdate (fld2,               halo_info, &
                           field_loc_NEcorner, field_type_vector)
      call unstack_velocity_field(fld2, uvel, vvel)
      call ice_timer_stop(timer_bound)

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
            call basal_stress_coeff (nx_block        , ny_block      , &
                                     icellu    (iblk),                 &
                                     indxui  (:,iblk), indxuj(:,iblk), &
                                     vice  (:,:,iblk), aice(:,:,iblk), &
                                     hwater(:,:,iblk), Tbu (:,:,iblk))
         enddo
         !$OMP END PARALLEL DO
      endif
      
      !-----------------------------------------------------------------
      ! calc size of problem (ntot) and allocate solution vector
      !-----------------------------------------------------------------
      
      ntot = 0
      do iblk = 1, nblocks
         ntot = ntot + icellu(iblk)
      enddo
      ntot = 2 * ntot ! times 2 because of u and v
      
      allocate(sol(ntot))
      
      !-----------------------------------------------------------------
      ! Start of nonlinear iteration
      !-----------------------------------------------------------------
      call anderson_solver (icellt  , icellu, &
                            indxti  , indxtj, &
                            indxui  , indxuj, &
                            aiu     , ntot  , &
                            waterx  , watery, &
                            bxfix   , byfix , &
                            umassdti, sol   , &
                            fpresx  , fpresy, &
                            zetaD   , Cb    , &
                            halo_info_mask)
      !-----------------------------------------------------------------
      ! End of nonlinear iteration
      !-----------------------------------------------------------------

      deallocate(sol)
      
      if (maskhalo_dyn) call ice_HaloDestroy(halo_info_mask)

      !-----------------------------------------------------------------
      ! Compute stresses
      !-----------------------------------------------------------------
      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks
         call stress_vp (nx_block            , ny_block            , &
                         icellt(iblk)        ,                       &
                         indxti      (:,iblk), indxtj      (:,iblk), &
                         uvel      (:,:,iblk), vvel      (:,:,iblk), &
                         dxt       (:,:,iblk), dyt       (:,:,iblk), &
                         cxp       (:,:,iblk), cyp       (:,:,iblk), &
                         cxm       (:,:,iblk), cym       (:,:,iblk), &
                         zetaD   (:,:,iblk,:),                       &
                         stressp_1 (:,:,iblk), stressp_2 (:,:,iblk), &
                         stressp_3 (:,:,iblk), stressp_4 (:,:,iblk), &
                         stressm_1 (:,:,iblk), stressm_2 (:,:,iblk), &
                         stressm_3 (:,:,iblk), stressm_4 (:,:,iblk), &
                         stress12_1(:,:,iblk), stress12_2(:,:,iblk), &
                         stress12_3(:,:,iblk), stress12_4(:,:,iblk))
      enddo ! iblk
      !$OMP END PARALLEL DO

      !-----------------------------------------------------------------
      ! Compute deformations
      !-----------------------------------------------------------------
      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks
         call deformations (nx_block            , ny_block            , &
                            icellt(iblk)        ,                       &
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
      
      !-----------------------------------------------------------------
      ! Compute seabed stress (diagnostic)
      !-----------------------------------------------------------------
      if (basalstress) then
         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks
            call calc_seabed_stress (nx_block            ,  ny_block           , &
                                     icellu(iblk)        ,                       &
                                     indxui      (:,iblk), indxuj      (:,iblk), &
                                     uvel      (:,:,iblk), vvel      (:,:,iblk), &
                                     Cb        (:,:,iblk),                       &
                                     taubx     (:,:,iblk), tauby     (:,:,iblk))
         enddo
         !$OMP END PARALLEL DO
      endif
      
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

      end subroutine implicit_solver

!=======================================================================

! Solve the nonlinear equation F(u,v) = 0, where
! F(u,v) := A(u,v) * (u,v) - b(u,v)
! using Anderson acceleration (accelerated fixed point (Picard) iteration)
!
! author: JF Lemieux, A. Qaddouri, F. Dupont and P. Blain ECCC
!
! Anderson algorithm adadpted from:
! H. F. Walker, “Anderson Acceleration: Algorithms and Implementations”
!   [Online]. Available: https://users.wpi.edu/~walker/Papers/anderson_accn_algs_imps.pdf

      subroutine anderson_solver (icellt  , icellu, &
                                  indxti  , indxtj, &
                                  indxui  , indxuj, &
                                  aiu     , ntot  , &
                                  waterx  , watery, &
                                  bxfix   , byfix , &
                                  umassdti, sol   , &
                                  fpresx  , fpresy, &
                                  zetaD   , Cb    , &
                                  halo_info_mask)

      use ice_arrays_column, only: Cdn_ocn
      use ice_blocks, only: nx_block, ny_block
      use ice_boundary, only: ice_HaloUpdate
      use ice_constants, only: c1
      use ice_domain, only: maskhalo_dyn, halo_info
      use ice_domain_size, only: max_blocks
      use ice_flux, only:   uocn, vocn, fm, Tbu
      use ice_grid, only: dxt, dyt, dxhy, dyhx, cxp, cyp, cxm, cym, &
          uarear, tinyarea
      use ice_state, only: uvel, vvel, strength
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bound

      integer (kind=int_kind), intent(in) :: &
         ntot         ! size of problem for Anderson

      integer (kind=int_kind), dimension(max_blocks), intent(in) :: &
         icellt   , & ! no. of cells where icetmask = 1
         icellu       ! no. of cells where iceumask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block, max_blocks), intent(in) :: &
         indxti   , & ! compressed index in i-direction
         indxtj   , & ! compressed index in j-direction
         indxui   , & ! compressed index in i-direction
         indxuj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), intent(in) :: &
         aiu      , & ! ice fraction on u-grid
         waterx   , & ! for ocean stress calculation, x (m/s)
         watery   , & ! for ocean stress calculation, y (m/s)
         bxfix    , & ! part of bx that is constant during Picard
         byfix    , & ! part of by that is constant during Picard
         umassdti     ! mass of U-cell/dte (kg/m^2 s)

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks,4), intent(out) :: &
         zetaD        ! zetaD = 2zeta (viscous coeff)

      type (ice_halo), intent(in) :: &
         halo_info_mask !  ghost cell update info for masked halo

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), intent(inout) :: &
         fpresx   , & ! fixed point residual vector, x components: fx = uvel - uprev_k
         fpresy   , & ! fixed point residual vector, y components: fy = vvel - vprev_k
         Cb           ! seabed stress coefficient

      real (kind=dbl_kind), dimension (ntot), intent(inout) :: &
         sol          ! current approximate solution

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
         bx       , & ! b vector
         by       , & ! b vector
         diagx    , & ! Diagonal (x component) of the matrix A
         diagy    , & ! Diagonal (y component) of the matrix A
         Au       , & ! matvec, Fx = bx - Au
         Av       , & ! matvec, Fy = by - Av
         Fx       , & ! x residual vector, Fx = bx - Au
         Fy       , & ! y residual vector, Fy = by - Av
         solx     , & ! solution of FGMRES (x components)
         soly         ! solution of FGMRES (y components)

      real (kind=dbl_kind), dimension(nx_block,ny_block,8):: &
         stress_Pr,   & ! x,y-derivatives of the replacement pressure
         diag_rheo      ! contributions of the rhelogy term to the diagonal

      real (kind=dbl_kind), dimension (max_blocks) :: &
         L2norm       ! array used to compute l^2 norm of grid function

      real (kind=dbl_kind), dimension (ntot) :: &
         res        , & ! current residual
         res_old    , & ! previous residual
         res_diff   , & ! difference between current and previous residuals
         fpfunc     , & ! current value of fixed point function
         fpfunc_old , & ! previous value of fixed point function
         tmp            ! temporary vector for BLAS calls

      real (kind=dbl_kind), dimension(ntot,dim_andacc) :: &
         Q        , & ! Q factor for QR factorization of F (residuals) matrix
         G_diff       ! Matrix containing the differences of g(x) (fixed point function) evaluations

      real (kind=dbl_kind), dimension(dim_andacc,dim_andacc) :: &
         R            ! R factor for QR factorization of F (residuals) matrix

      real (kind=dbl_kind), dimension(dim_andacc) :: &
         rhs_tri  , & ! right hand side vector for matrix-vector product
         coeffs       ! coeffs used to combine previous solutions

      real (kind=dbl_kind) :: &
         ! tol         , & ! tolerance for fixed point convergence: reltol_andacc * (initial fixed point residual norm)  [unused for now]
         tol_nl      , & ! tolerance for nonlinear convergence: reltol_nonlin * (initial nonlinear residual norm)
         fpres_norm  , & ! norm of current fixed point residual : f(x) = g(x) - x
         prog_norm   , & ! norm of difference between current and previous solution
         nlres_norm      ! norm of current nonlinear residual : F(x) = A(x)x -b(x)

#ifdef USE_LAPACK
      real (kind=dbl_kind) :: &
         ddot, dnrm2     ! external BLAS functions
#endif

      character(len=*), parameter :: subname = '(anderson_solver)'

      ! Initialization
      res_num = 0
      L2norm  = c0
      
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
         ! Calc zetaD, dPr/dx, dPr/dy, Cb and vrel = f(uprev_k, vprev_k)
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
            
            call calc_zeta_dPr (nx_block           , ny_block          , &
                                icellt       (iblk),                     &
                                indxti     (:,iblk), indxtj    (:,iblk), &
                                uprev_k  (:,:,iblk), vprev_k (:,:,iblk), &
                                dxt      (:,:,iblk), dyt     (:,:,iblk), &
                                dxhy     (:,:,iblk), dyhx    (:,:,iblk), &
                                cxp      (:,:,iblk), cyp     (:,:,iblk), &
                                cxm      (:,:,iblk), cym     (:,:,iblk), &
                                tinyarea (:,:,iblk),                     &
                                strength (:,:,iblk), zetaD (:,:,iblk,:), &
                                stress_Pr  (:,:,:))
            
            call calc_vrel_Cb (nx_block           , ny_block          , &
                               icellu       (iblk), Cdn_ocn (:,:,iblk), &
                               indxui     (:,iblk), indxuj    (:,iblk), &
                               aiu      (:,:,iblk), Tbu     (:,:,iblk), &
                               uocn     (:,:,iblk), vocn    (:,:,iblk), &
                               ulin     (:,:,iblk), vlin    (:,:,iblk), &
                               vrel     (:,:,iblk), Cb      (:,:,iblk))
            
            ! prepare b vector (RHS)
            call calc_bvec (nx_block           , ny_block          , &
                            icellu       (iblk),                     &
                            indxui     (:,iblk), indxuj    (:,iblk), &
                            stress_Pr   (:,:,:), uarear  (:,:,iblk), &
                            waterx   (:,:,iblk), watery  (:,:,iblk), &
                            bxfix    (:,:,iblk), byfix   (:,:,iblk), &
                            bx       (:,:,iblk), by      (:,:,iblk), &
                            vrel     (:,:,iblk))
            
            ! Compute nonlinear residual norm (PDE residual)
            call matvec (nx_block             , ny_block           , &
                         icellu       (iblk)  , icellt       (iblk), &
                         indxui     (:,iblk)  , indxuj     (:,iblk), &
                         indxti     (:,iblk)  , indxtj     (:,iblk), &
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
            call residual_vec (nx_block           , ny_block          , &
                               icellu       (iblk),                     &
                               indxui     (:,iblk), indxuj    (:,iblk), &
                               bx       (:,:,iblk), by      (:,:,iblk), &
                               Au       (:,:,iblk), Av      (:,:,iblk), &
                               Fx       (:,:,iblk), Fy      (:,:,iblk), &
                               L2norm       (iblk))
         enddo
         !$OMP END PARALLEL DO
         nlres_norm = sqrt(global_sum(sum(L2norm), distrb_info))
         if (my_task == master_task .and. monitor_nonlin) then
            write(nu_diag, '(a,i4,a,d26.16)') "monitor_nonlin: iter_nonlin= ", it_nl, &
                                              " nonlin_res_L2norm= ", nlres_norm
         endif
         ! Compute relative tolerance at first iteration
         if (it_nl == 0) then
            tol_nl = reltol_nonlin*nlres_norm
         endif
         
         ! Check for nonlinear convergence
         if (nlres_norm < tol_nl) then
            exit
         endif
         
         ! Put initial guess for FGMRES in solx,soly and sol (needed for anderson)
         solx = uprev_k
         soly = vprev_k
         call arrays_to_vec (nx_block       , ny_block       , &
                             nblocks        , max_blocks     , &
                             icellu      (:), ntot           , &
                             indxui    (:,:), indxuj    (:,:), &
                             uprev_k (:,:,:), vprev_k (:,:,:), &
                             sol         (:))
         
         ! Compute fixed point map g(x)
         if (fpfunc_andacc == 1) then
            ! g_1(x) = FGMRES(A(x), b(x))
            
            ! Prepare diagonal for preconditioner
            if (precond == 'diag' .or. precond == 'pgmres') then
               !$OMP PARALLEL DO PRIVATE(iblk)
               do iblk = 1, nblocks
                  ! first compute diagonal contributions due to rheology term
                  call formDiag_step1 (nx_block             , ny_block      , &
                                       icellu       (iblk)  ,                 &
                                       indxui     (:,iblk)  , indxuj(:,iblk), &
                                       dxt      (:,:,iblk)  , dyt (:,:,iblk), &
                                       dxhy     (:,:,iblk)  , dyhx(:,:,iblk), &
                                       cxp      (:,:,iblk)  , cyp (:,:,iblk), &
                                       cxm      (:,:,iblk)  , cym (:,:,iblk), &
                                       zetaD    (:,:,iblk,:), diag_rheo(:,:,:))
                  ! second compute the full diagonal
                  call formDiag_step2 (nx_block           , ny_block          , &
                                       icellu       (iblk),                     &
                                       indxui     (:,iblk), indxuj    (:,iblk), &
                                       diag_rheo   (:,:,:), vrel    (:,:,iblk), &
                                       umassdti (:,:,iblk),                     &
                                       uarear   (:,:,iblk), Cb      (:,:,iblk), &
                                       diagx    (:,:,iblk), diagy   (:,:,iblk))
               enddo
               !$OMP END PARALLEL DO
            endif
            
            ! FGMRES linear solver
            call fgmres (zetaD         ,             &
                         Cb            , vrel      , &
                         umassdti      ,             &
                         halo_info_mask,             &
                         bx            , by        , &
                         diagx         , diagy     , &
                         reltol_fgmres , dim_fgmres, &
                         maxits_fgmres ,             &
                         solx          , soly      , &
                         nbiter)
            ! Put FGMRES solution solx,soly in fpfunc vector (needed for Anderson)
            call arrays_to_vec (nx_block       , ny_block     , &
                                nblocks        , max_blocks   , &
                                icellu      (:), ntot         , &
                                indxui    (:,:), indxuj  (:,:), &
                                solx    (:,:,:), soly  (:,:,:), &
                                fpfunc      (:))
         elseif (fpfunc_andacc == 2) then
            ! g_2(x) = x - A(x)x + b(x) = x - F(x)
            call abort_ice(error_message=subname // " Fixed point function g_2(x) not yet implemented (fpfunc_andacc = 2)" , &
               file=__FILE__, line=__LINE__)
         endif

         ! Compute fixed point residual f(x) = g(x) - x
         res = fpfunc - sol
#ifdef USE_LAPACK
         fpres_norm = global_sum(dnrm2(size(res), res, inc)**2, distrb_info)
#else
         call vec_to_arrays (nx_block      , ny_block    , &
                             nblocks       , max_blocks  , &
                             icellu     (:), ntot        , &
                             indxui   (:,:), indxuj(:,:) , &
                             res        (:),               &
                             fpresx (:,:,:), fpresy (:,:,:))
         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks
            call calc_L2norm_squared (nx_block        , ny_block        , &
                                      icellu    (iblk),                   &
                                      indxui  (:,iblk), indxuj  (:,iblk), &
                                      fpresx(:,:,iblk), fpresy(:,:,iblk), &
                                      L2norm    (iblk))
         enddo
         !$OMP END PARALLEL DO
         fpres_norm = sqrt(global_sum(sum(L2norm), distrb_info))
#endif
         if (my_task == master_task .and. monitor_nonlin) then
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

         if (dim_andacc == 0 .or. it_nl < start_andacc) then
            ! Simple fixed point (Picard) iteration in this case
            sol = fpfunc
         else
#ifdef USE_LAPACK
            ! Begin Anderson acceleration
            if (get_num_procs() > 1) then
               ! Anderson solver is not yet parallelized; abort
               if (my_task == master_task) then
                  call abort_ice(error_message=subname // " Anderson solver (algo_nonlin = 'anderson') is not yet parallelized, and nprocs > 1 " , &
                     file=__FILE__, line=__LINE__)
               endif
            endif
            if (it_nl > start_andacc) then
               ! Update residual difference vector
               res_diff = res - res_old
               ! Update fixed point function difference matrix
               if (res_num < dim_andacc) then
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
                  if (res_num > dim_andacc) then
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
               ! TODO: here, drop more columns to improve conditioning
               ! if (droptol) then
               
               ! endif
               ! Solve least square problem for coefficients
               ! 1. Compute rhs_tri = Q^T * res
               call dgemv ('t', size(Q,1), res_num, c1, Q(:,1:res_num), size(Q,1), res, inc, c0, rhs_tri, inc)
               ! 2. Solve R*coeffs = rhs_tri, put result in rhs_tri
               call dtrsv ('u', 'n', 'n', res_num, R(1:res_num,1:res_num), res_num, rhs_tri, inc)
               coeffs = rhs_tri
               ! Update approximate solution: x = fpfunc - G_diff*coeffs, put result in fpfunc
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
#else
            ! Anderson solver is not usable without LAPACK; abort
            call abort_ice(error_message=subname // " CICE was not compiled with LAPACK, and Anderson solver was chosen (algo_nonlin = 'anderson')" , &
               file=__FILE__, line=__LINE__)
#endif
         endif
         
         !-----------------------------------------------------------------------
         !     Put vector sol in uvel and vvel arrays
         !-----------------------------------------------------------------------
         call vec_to_arrays (nx_block    , ny_block    , &
                             nblocks     , max_blocks  , &
                             icellu   (:), ntot        , &
                             indxui (:,:), indxuj (:,:), &
                             sol      (:),               &
                             uvel (:,:,:), vvel (:,:,:))
         
         ! Do halo update so that halo cells contain up to date info for advection
         call stack_velocity_field(uvel, vvel, fld2)
         call ice_timer_start(timer_bound)
         if (maskhalo_dyn) then
            call ice_HaloUpdate (fld2,               halo_info_mask, &
                                 field_loc_NEcorner, field_type_vector)
         else
            call ice_HaloUpdate (fld2,               halo_info, &
                                 field_loc_NEcorner, field_type_vector)
         endif
         call ice_timer_stop(timer_bound)
         call unstack_velocity_field(fld2, uvel, vvel)
         
         ! Compute "progress" residual norm
         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks
            fpresx(:,:,iblk) = uvel(:,:,iblk) - uprev_k(:,:,iblk)
            fpresy(:,:,iblk) = vvel(:,:,iblk) - vprev_k(:,:,iblk)
            call calc_L2norm_squared (nx_block        , ny_block        , &
                                      icellu    (iblk),                   &
                                      indxui  (:,iblk), indxuj  (:,iblk), &
                                      fpresx(:,:,iblk), fpresy(:,:,iblk), &
                                      L2norm    (iblk))
         enddo
         !$OMP END PARALLEL DO
         prog_norm = sqrt(global_sum(sum(L2norm), distrb_info))
         if (my_task == master_task .and. monitor_nonlin) then
            write(nu_diag, '(a,i4,a,d26.16)') "monitor_nonlin: iter_nonlin= ", it_nl, &
                                              " progress_res_L2norm= ", prog_norm
         endif
         
      enddo ! nonlinear iteration loop
      
      end subroutine anderson_solver

!=======================================================================

! Computes the viscous coefficients (in fact zetaD=2*zeta) and dPr/dx, dPr/dy

      subroutine calc_zeta_dPr (nx_block, ny_block, &
                                icellt  ,           &
                                indxti  , indxtj  , &
                                uvel    , vvel    , &
                                dxt     , dyt     , &
                                dxhy    , dyhx    , &
                                cxp     , cyp     , &
                                cxm     , cym     , &
                                tinyarea,           &
                                strength, zetaD   , &
                                stPr)

      use ice_dyn_shared, only: strain_rates

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellt                ! no. of cells where icetmask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
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
         tinyarea     ! min_strain_rate*tarea

      real (kind=dbl_kind), dimension(nx_block,ny_block,4), intent(out) :: &
         zetaD          ! 2*zeta

      real (kind=dbl_kind), dimension(nx_block,ny_block,8), intent(out) :: &
         stPr          ! stress combinations from replacement pressure

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

      character(len=*), parameter :: subname = '(calc_zeta_dPr)'

      ! Initialize

      capping = .false.
      
      ! Initialize stPr and zetaD to zero (for cells where icetmask is false)
      stPr  = c0
      zetaD = c0

      do ij = 1, icellt
         i = indxti(ij)
         j = indxtj(ij)

      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------
         call strain_rates (nx_block , ny_block , &
                            i        , j        , &
                            uvel     , vvel     , &
                            dxt      , dyt      , &
                            cxp      , cyp      , &
                            cxm      , cym      , &
                            divune   , divunw   , &
                            divuse   , divusw   , &
                            tensionne, tensionnw, &
                            tensionse, tensionsw, &
                            shearne  , shearnw  , &
                            shearse  , shearsw  , &
                            Deltane  , Deltanw  , &
                            Deltase  , Deltasw)

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

      end subroutine calc_zeta_dPr

!=======================================================================

! Computes the VP stresses (as diagnostic)

      subroutine stress_vp (nx_block  , ny_block  , &
                            icellt    ,             &
                            indxti    , indxtj    , &
                            uvel      , vvel      , &
                            dxt       , dyt       , &
                            cxp       , cyp       , &
                            cxm       , cym       , &
                            zetaD     ,             &
                            stressp_1 , stressp_2 , &
                            stressp_3 , stressp_4 , &
                            stressm_1 , stressm_2 , &
                            stressm_3 , stressm_4 , &
                            stress12_1, stress12_2, &
                            stress12_3, stress12_4)

      use ice_dyn_shared, only: strain_rates

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellt                ! no. of cells where icetmask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
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
         cxm          ! 0.5*HTN - 1.5*HTN

      real (kind=dbl_kind), dimension(nx_block,ny_block,4), intent(in) :: &
         zetaD          ! 2*zeta

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         stressp_1, stressp_2, stressp_3, stressp_4 , & ! sigma11+sigma22
         stressm_1, stressm_2, stressm_3, stressm_4 , & ! sigma11-sigma22
         stress12_1,stress12_2,stress12_3,stress12_4    ! sigma12

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &
        divune, divunw, divuse, divusw            , & ! divergence
        tensionne, tensionnw, tensionse, tensionsw, & ! tension
        shearne, shearnw, shearse, shearsw        , & ! shearing
        Deltane, Deltanw, Deltase, Deltasw            ! Delt

      character(len=*), parameter :: subname = '(stress_vp)'

      do ij = 1, icellt
         i = indxti(ij)
         j = indxtj(ij)

      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------
         call strain_rates (nx_block , ny_block , &
                            i        , j        , &
                            uvel     , vvel     , &
                            dxt      , dyt      , &
                            cxp      , cyp      , &
                            cxm      , cym      , &
                            divune   , divunw   , &
                            divuse   , divusw   , &
                            tensionne, tensionnw, &
                            tensionse, tensionsw, &
                            shearne  , shearnw  , &
                            shearse  , shearsw  , &
                            Deltane  , Deltanw  , &
                            Deltase  , Deltasw)

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

      enddo                     ! ij

      end subroutine stress_vp

!=======================================================================

! Compute vrel and seabed stress coefficients

      subroutine calc_vrel_Cb (nx_block, ny_block, &
                               icellu  , Cw      , &
                               indxui  , indxuj  , &
                               aiu     , Tbu     , &
                               uocn    , vocn    , &
                               uvel    , vvel    , &
                               vrel    , Cb)

      use ice_dyn_shared, only: u0 ! residual velocity for basal stress (m/s)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellu                ! total count when iceumask is true

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxui  , & ! compressed index in i-direction
         indxuj      ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         Tbu,      & ! coefficient for basal stress (N/m^2)
         aiu     , & ! ice fraction on u-grid
         uocn    , & ! ocean current, x-direction (m/s)
         vocn    , & ! ocean current, y-direction (m/s)
         Cw          ! ocean-ice neutral drag coefficient

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         uvel    , & ! x-component of velocity (m/s)
         vvel        ! y-component of velocity (m/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         vrel      , & ! coeff for tauw
         Cb            ! seabed stress coeff

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &
         rhow                  !

      character(len=*), parameter :: subname = '(calc_vrel_Cb)'

      call icepack_query_parameters(rhow_out=rhow)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      do ij = 1, icellu
         i = indxui(ij)
         j = indxuj(ij)

         ! (magnitude of relative ocean current)*rhow*drag*aice
         vrel(i,j) = aiu(i,j)*rhow*Cw(i,j)*sqrt((uocn(i,j) - uvel(i,j))**2 + &
                                                (vocn(i,j) - vvel(i,j))**2)  ! m/s
      
         Cb(i,j)  = Tbu(i,j) / (sqrt(uvel(i,j)**2 + vvel(i,j)**2) + u0) ! for seabed stress
      enddo                     ! ij

      end subroutine calc_vrel_Cb

!=======================================================================

! Compute seabed stress (diagnostic)

      subroutine calc_seabed_stress (nx_block, ny_block, &
                                     icellu  ,           &
                                     indxui  , indxuj  , &
                                     uvel    , vvel    , &
                                     Cb      ,           &
                                     taubx   , tauby)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellu                ! total count when iceumask is true

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxui  , & ! compressed index in i-direction
         indxuj      ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         uvel    , & ! x-component of velocity (m/s)
         vvel    , & ! y-component of velocity (m/s)
         Cb          ! seabed stress coefficient

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out) :: &
         taubx   , & ! seabed stress, x-direction (N/m^2)
         tauby       ! seabed stress, y-direction (N/m^2)

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      character(len=*), parameter :: subname = '(calc_seabed_stress)'

      do ij = 1, icellu
         i = indxui(ij)
         j = indxuj(ij)
         
         taubx(i,j) = -uvel(i,j)*Cb(i,j)
         tauby(i,j) = -vvel(i,j)*Cb(i,j)
      enddo                     ! ij

      end subroutine calc_seabed_stress

!=======================================================================

! Computes the matrix vector product A(u,v) * (u,v)
! Au = A(u,v)_[x] * uvel (x components of  A(u,v) * (u,v))
! Av = A(u,v)_[y] * vvel (y components of  A(u,v) * (u,v))

      subroutine matvec (nx_block, ny_block, &
                         icellu  , icellt  , &
                         indxui  , indxuj  , &
                         indxti  , indxtj  , &
                         dxt     , dyt     , &
                         dxhy    , dyhx    , &
                         cxp     , cyp     , &
                         cxm     , cym     , &
                         uvel    , vvel    , &
                         vrel    , Cb      , &
                         zetaD   ,           &
                         umassdti, fm      , &
                         uarear  ,           &
                         Au      , Av)

      use ice_dyn_shared, only: strain_rates

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellu,             & ! total count when iceumask is true
         icellt                ! no. of cells where icetmask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
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

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         uvel    , & ! x-component of velocity (m/s)
         vvel    , & ! y-component of velocity (m/s)
         vrel    , & ! coefficient for tauw
         Cb      , & ! coefficient for basal stress
         umassdti, & ! mass of U-cell/dt (kg/m^2 s)
         fm      , & ! Coriolis param. * mass in U-cell (kg/s)
         uarear      ! 1/uarea

      real (kind=dbl_kind), dimension(nx_block,ny_block,4), intent(in) :: &
         zetaD          ! 2*zeta

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         Au      , & ! matvec, Fx = bx - Au (N/m^2)
         Av          ! matvec, Fy = by - Av (N/m^2)

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind), dimension(nx_block,ny_block,8):: &
         str

      real (kind=dbl_kind) :: &
         ccaimp,ccb        , & ! intermediate variables
         strintx, strinty      ! divergence of the internal stress tensor

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

      do ij = 1, icellt
         i = indxti(ij)
         j = indxtj(ij)

      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------
         call strain_rates (nx_block , ny_block , &
                            i        , j        , &
                            uvel     , vvel     , &
                            dxt      , dyt      , &
                            cxp      , cyp      , &
                            cxm      , cym      , &
                            divune   , divunw   , &
                            divuse   , divusw   , &
                            tensionne, tensionnw, &
                            tensionse, tensionsw, &
                            shearne  , shearnw  , &
                            shearse  , shearsw  , &
                            Deltane  , Deltanw  , &
                            Deltase  , Deltasw)

      !-----------------------------------------------------------------
      ! the stresses                            ! kg/s^2
      ! (1) northeast, (2) northwest, (3) southwest, (4) southeast
      ! NOTE: commented part of stressp is from the replacement pressure Pr
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

      do ij = 1, icellu
         i = indxui(ij)
         j = indxuj(ij)

         ccaimp = umassdti(i,j) + vrel(i,j) * cosw + Cb(i,j) ! kg/m^2 s
         
         ccb = fm(i,j) + sign(c1,fm(i,j)) * vrel(i,j) * sinw ! kg/m^2 s

         ! divergence of the internal stress tensor
         strintx = uarear(i,j)* &
             (str(i,j,1) + str(i+1,j,2) + str(i,j+1,3) + str(i+1,j+1,4))
         strinty = uarear(i,j)* &
             (str(i,j,5) + str(i,j+1,6) + str(i+1,j,7) + str(i+1,j+1,8))

         Au(i,j) = ccaimp*uvel(i,j) - ccb*vvel(i,j) - strintx
         Av(i,j) = ccaimp*vvel(i,j) + ccb*uvel(i,j) - strinty
      enddo ! ij - icellu

      end subroutine matvec

!=======================================================================

! Compute the constant component of b(u,v) i.e. the part of b(u,v) that
! does not depend on (u,v) and thus do not change during the nonlinear iteration

     subroutine calc_bfix  (nx_block , ny_block , &
                            icellu   ,            &
                            indxui   , indxuj   , &
                            umassdti ,            &
                            forcex   , forcey   , &
                            uvel_init, vvel_init, &
                            bxfix    , byfix)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellu                ! no. of cells where iceumask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxui   , & ! compressed index in i-direction
         indxuj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         uvel_init,& ! x-component of velocity (m/s), beginning of time step
         vvel_init,& ! y-component of velocity (m/s), beginning of time step
         umassdti, & ! mass of U-cell/dt (kg/m^2 s)
         forcex  , & ! work array: combined atm stress and ocn tilt, x
         forcey      ! work array: combined atm stress and ocn tilt, y

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out) :: &
         bxfix   , & ! bx = taux + bxfix
         byfix       ! by = tauy + byfix

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      character(len=*), parameter :: subname = '(calc_bfix)'

      do ij = 1, icellu
         i = indxui(ij)
         j = indxuj(ij)

         bxfix(i,j) = umassdti(i,j)*uvel_init(i,j) + forcex(i,j)
         byfix(i,j) = umassdti(i,j)*vvel_init(i,j) + forcey(i,j)
      enddo

      end subroutine calc_bfix

!=======================================================================

! Compute the vector b(u,v), i.e. the part of the nonlinear function F(u,v)
! that cannot be written as A(u,v)*(u,v), where A(u,v) is a matrix with entries
! depending on (u,v)

      subroutine calc_bvec (nx_block, ny_block, &
                            icellu  ,           &
                            indxui  , indxuj  , &
                            stPr    , uarear  , &
                            waterx  , watery  , &
                            bxfix   , byfix   , &
                            bx      , by      , &
                            vrel)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellu                ! total count when iceumask is true

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxui  , & ! compressed index in i-direction
         indxuj      ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         uarear  , & ! 1/uarea
         waterx  , & ! for ocean stress calculation, x (m/s)
         watery  , & ! for ocean stress calculation, y (m/s)
         bxfix   , & ! bx = taux + bxfix
         byfix   , & ! by = tauy + byfix
         vrel        ! relative ice-ocean velocity

      real (kind=dbl_kind), dimension(nx_block,ny_block,8), intent(in) :: &
         stPr

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out) :: &
         bx      , & ! b vector, bx = taux + bxfix (N/m^2)
         by          ! b vector, by = tauy + byfix (N/m^2)

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &
         taux, tauy        , & ! part of ocean stress term
         strintx, strinty  , & ! divergence of the internal stress tensor (only Pr contributions)
         rhow                  !

      character(len=*), parameter :: subname = '(calc_bvec)'
      
      !-----------------------------------------------------------------
      ! calc b vector
      !-----------------------------------------------------------------
      
      call icepack_query_parameters(rhow_out=rhow)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      do ij = 1, icellu
         i = indxui(ij)
         j = indxuj(ij)

         ! ice/ocean stress
         taux = vrel(i,j)*waterx(i,j) ! NOTE this is not the entire
         tauy = vrel(i,j)*watery(i,j) ! ocn stress term
         
         ! divergence of the internal stress tensor (only Pr part, i.e. dPr/dx, dPr/dy)
         strintx = uarear(i,j)* &
             (stPr(i,j,1) + stPr(i+1,j,2) + stPr(i,j+1,3) + stPr(i+1,j+1,4))
         strinty = uarear(i,j)* &
             (stPr(i,j,5) + stPr(i,j+1,6) + stPr(i+1,j,7) + stPr(i+1,j+1,8))

         bx(i,j) = bxfix(i,j) + taux + strintx
         by(i,j) = byfix(i,j) + tauy + strinty
      enddo                     ! ij

      end subroutine calc_bvec

!=======================================================================

! Compute the non linear residual F(u,v) = b(u,v) - A(u,v) * (u,v),
! with Au, Av precomputed as
! Au = A(u,v)_[x] * uvel (x components of  A(u,v) * (u,v))
! Av = A(u,v)_[y] * vvel (y components of  A(u,v) * (u,v))

      subroutine residual_vec (nx_block   , ny_block, &
                               icellu     ,           &
                               indxui     , indxuj  , &
                               bx         , by      , &
                               Au         , Av      , &
                               Fx         , Fy      , &
                               sum_squared)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellu                ! total count when iceumask is true

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxui  , & ! compressed index in i-direction
         indxuj      ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         bx       , & ! b vector, bx = taux + bxfix (N/m^2)
         by       , & ! b vector, by = tauy + byfix (N/m^2)
         Au       , & ! matvec, Fx = bx - Au (N/m^2)
         Av           ! matvec, Fy = by - Av (N/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         Fx      , & ! x residual vector, Fx = bx - Au (N/m^2)
         Fy          ! y residual vector, Fy = by - Av (N/m^2)

      real (kind=dbl_kind), intent(out), optional :: &
         sum_squared ! sum of squared residual vector components

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      character(len=*), parameter :: subname = '(residual_vec)'

      !-----------------------------------------------------------------
      ! compute residual and sum its squared components
      !-----------------------------------------------------------------

      if (present(sum_squared)) then
         sum_squared = c0
      endif
      
      do ij = 1, icellu
         i = indxui(ij)
         j = indxuj(ij)

         Fx(i,j) = bx(i,j) - Au(i,j)
         Fy(i,j) = by(i,j) - Av(i,j)
         if (present(sum_squared)) then
            sum_squared = sum_squared + Fx(i,j)**2 + Fy(i,j)**2
         endif
      enddo                     ! ij

      end subroutine residual_vec

!=======================================================================

! Form the diagonal of the matrix A(u,v) (first part of the computation)
! Part 1: compute the contributions to the diagonal from the rheology term

      subroutine formDiag_step1  (nx_block, ny_block, &
                                  icellu  ,           &
                                  indxui  , indxuj  , &
                                  dxt     , dyt     , &
                                  dxhy    , dyhx    , &
                                  cxp     , cyp     , &
                                  cxm     , cym     , &
                                  zetaD   , Drheo)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellu                ! no. of cells where icetmask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxui   , & ! compressed index in i-direction
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

      real (kind=dbl_kind), dimension(nx_block,ny_block,4), intent(in) :: &
         zetaD          ! 2*zeta

      real (kind=dbl_kind), dimension(nx_block,ny_block,8), intent(out) :: &
         Drheo          ! intermediate value for diagonal components of matrix A associated
                       ! with rheology term

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij, iu, ju, di, dj, cc

      real (kind=dbl_kind) :: &
        divune, divunw, divuse, divusw            , & ! divergence
        tensionne, tensionnw, tensionse, tensionsw, & ! tension
        shearne, shearnw, shearse, shearsw        , & ! shearing
        uij,ui1j,uij1,ui1j1,vij,vi1j,vij1,vi1j1   , & ! == c0 or c1
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

      Drheo(:,:,:) = c0
      
      ! Be careful: Drheo contains 4 terms for u and 4 terms for v.
      ! These 8 terms come from the surrounding T cells but are all
      ! refrerenced to the i,j (u point) :
      
      ! Drheo(i,j,1) corresponds to str(i,j,1)
      ! Drheo(i,j,2) corresponds to str(i+1,j,2)
      ! Drheo(i,j,3) corresponds to str(i,j+1,3)
      ! Drheo(i,j,4) corresponds to str(i+1,j+1,4))
      ! Drheo(i,j,5) corresponds to str(i,j,5)
      ! Drheo(i,j,6) corresponds to str(i,j+1,6)
      ! Drheo(i,j,7) corresponds to str(i+1,j,7)
      ! Drheo(i,j,8) corresponds to str(i+1,j+1,8))
      
      do cc = 1, 8 ! 4 for u and 4 for v
      
         if (cc == 1) then     ! u comp, T cell i,j
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
         elseif (cc == 2) then ! u comp, T cell i+1,j
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
         elseif (cc == 3) then ! u comp, T cell i,j+1
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
         elseif (cc == 4) then ! u comp, T cell i+1,j+1
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
         elseif (cc == 5) then ! v comp, T cell i,j
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
         elseif (cc == 6) then ! v comp, T cell i,j+1
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
         elseif (cc == 7) then ! v comp, T cell i+1,j
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
         elseif (cc == 8) then ! v comp, T cell i+1,j+1
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
            i  = iu + di
            j  = ju + dj
             
         !-----------------------------------------------------------------
         ! strain rates
         ! NOTE these are actually strain rates * area  (m^2/s)
         !-----------------------------------------------------------------
            ! divergence  =  e_11 + e_22
            divune    = cyp(i,j)*uij   - dyt(i,j)*ui1j  &
                      + cxp(i,j)*vij   - dxt(i,j)*vij1
            divunw    = cym(i,j)*ui1j  + dyt(i,j)*uij   &
                      + cxp(i,j)*vi1j  - dxt(i,j)*vi1j1
            divusw    = cym(i,j)*ui1j1 + dyt(i,j)*uij1  &
                      + cxm(i,j)*vi1j1 + dxt(i,j)*vi1j
            divuse    = cyp(i,j)*uij1  - dyt(i,j)*ui1j1 &
                      + cxm(i,j)*vij1  + dxt(i,j)*vij

            ! tension strain rate  =  e_11 - e_22
            tensionne = -cym(i,j)*uij   - dyt(i,j)*ui1j  &
                      +  cxm(i,j)*vij   + dxt(i,j)*vij1
            tensionnw = -cyp(i,j)*ui1j  + dyt(i,j)*uij   &
                      +  cxm(i,j)*vi1j  + dxt(i,j)*vi1j1
            tensionsw = -cyp(i,j)*ui1j1 + dyt(i,j)*uij1  &
                      +  cxp(i,j)*vi1j1 - dxt(i,j)*vi1j
            tensionse = -cym(i,j)*uij1  - dyt(i,j)*ui1j1 &
                      +  cxp(i,j)*vij1  - dxt(i,j)*vij

            ! shearing strain rate  =  2*e_12
            shearne = -cym(i,j)*vij   - dyt(i,j)*vi1j  &
                    -  cxm(i,j)*uij   - dxt(i,j)*uij1
            shearnw = -cyp(i,j)*vi1j  + dyt(i,j)*vij   &
                    -  cxm(i,j)*ui1j  - dxt(i,j)*ui1j1
            shearsw = -cyp(i,j)*vi1j1 + dyt(i,j)*vij1  &
                    -  cxp(i,j)*ui1j1 + dxt(i,j)*ui1j
            shearse = -cym(i,j)*vij1  - dyt(i,j)*vi1j1 &
                    -  cxp(i,j)*uij1  + dxt(i,j)*uij
            
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
         
            if (cc == 1) then ! T cell i,j
            
               strp_tmp  = p25*dyt(i,j)*(p333*ssigpn  + p166*ssigps)
               strm_tmp  = p25*dyt(i,j)*(p333*ssigmn  + p166*ssigms)

               ! northeast (i,j)
               Drheo(iu,ju,1) = -strp_tmp - strm_tmp - str12ew &
                  + dxhy(i,j)*(-csigpne + csigmne) + dyhx(i,j)*csig12ne
               
            elseif (cc == 2) then ! T cell i+1,j
               
               strp_tmp  = p25*dyt(i,j)*(p333*ssigpn  + p166*ssigps)
               strm_tmp  = p25*dyt(i,j)*(p333*ssigmn  + p166*ssigms)
               
               ! northwest (i+1,j)
               Drheo(iu,ju,2) = strp_tmp + strm_tmp - str12we &
                  + dxhy(i,j)*(-csigpnw + csigmnw) + dyhx(i,j)*csig12nw

            elseif (cc == 3) then ! T cell i,j+1
               
               strp_tmp  = p25*dyt(i,j)*(p333*ssigps  + p166*ssigpn)
               strm_tmp  = p25*dyt(i,j)*(p333*ssigms  + p166*ssigmn)

               ! southeast (i,j+1)
               Drheo(iu,ju,3) = -strp_tmp - strm_tmp + str12ew &
                  + dxhy(i,j)*(-csigpse + csigmse) + dyhx(i,j)*csig12se

            elseif (cc == 4) then ! T cell i+1,j+1
                 
               strp_tmp  = p25*dyt(i,j)*(p333*ssigps  + p166*ssigpn)
               strm_tmp  = p25*dyt(i,j)*(p333*ssigms  + p166*ssigmn)
               
               ! southwest (i+1,j+1)
               Drheo(iu,ju,4) = strp_tmp + strm_tmp + str12we &
                  + dxhy(i,j)*(-csigpsw + csigmsw) + dyhx(i,j)*csig12sw

         !-----------------------------------------------------------------
         ! for dF/dy (v momentum)
         !-----------------------------------------------------------------
            
            elseif (cc == 5) then ! T cell i,j
               
               strp_tmp  = p25*dxt(i,j)*(p333*ssigpe  + p166*ssigpw)
               strm_tmp  = p25*dxt(i,j)*(p333*ssigme  + p166*ssigmw)

               ! northeast (i,j)
               Drheo(iu,ju,5) = -strp_tmp + strm_tmp - str12ns &
                  - dyhx(i,j)*(csigpne + csigmne) + dxhy(i,j)*csig12ne

            elseif (cc == 6) then ! T cell i,j+1
               
               strp_tmp  = p25*dxt(i,j)*(p333*ssigpe  + p166*ssigpw)
               strm_tmp  = p25*dxt(i,j)*(p333*ssigme  + p166*ssigmw)
               
               ! southeast (i,j+1)
               Drheo(iu,ju,6) = strp_tmp - strm_tmp - str12sn &
                  - dyhx(i,j)*(csigpse + csigmse) + dxhy(i,j)*csig12se

            elseif (cc == 7) then ! T cell i,j+1
               
               strp_tmp  = p25*dxt(i,j)*(p333*ssigpw  + p166*ssigpe)
               strm_tmp  = p25*dxt(i,j)*(p333*ssigmw  + p166*ssigme)

               ! northwest (i+1,j)
               Drheo(iu,ju,7) = -strp_tmp + strm_tmp + str12ns &
                  - dyhx(i,j)*(csigpnw + csigmnw) + dxhy(i,j)*csig12nw

            elseif (cc == 8) then ! T cell i+1,j+1
               
               strp_tmp  = p25*dxt(i,j)*(p333*ssigpw  + p166*ssigpe)
               strm_tmp  = p25*dxt(i,j)*(p333*ssigmw  + p166*ssigme)
               
               ! southwest (i+1,j+1)
               Drheo(iu,ju,8) = strp_tmp - strm_tmp + str12sn &
                  - dyhx(i,j)*(csigpsw + csigmsw) + dxhy(i,j)*csig12sw
               
            endif

         enddo                     ! ij

      enddo                     ! cc

      end subroutine formDiag_step1

!=======================================================================

! Form the diagonal of the matrix A(u,v) (second part of the computation)
! Part 2: compute diagonal

      subroutine formDiag_step2 (nx_block, ny_block, &
                                 icellu  ,           &
                                 indxui  , indxuj  , &
                                 Drheo   , vrel    , &
                                 umassdti,           &
                                 uarear  , Cb      , &
                                 diagx   , diagy)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellu                ! total count when iceumask is true

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxui  , & ! compressed index in i-direction
         indxuj      ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         vrel,     & ! coefficient for tauw
         Cb,       & ! coefficient for basal stress
         umassdti, & ! mass of U-cell/dt (kg/m^2 s)
         uarear      ! 1/uarea

      real (kind=dbl_kind), dimension(nx_block,ny_block,8), intent(in) :: &
         Drheo

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out) :: &
         diagx   , & ! Diagonal (x component) of the matrix A
         diagy       ! Diagonal (y component) of the matrix A

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &
         ccaimp             , & ! intermediate variables
         strintx, strinty       ! diagonal contributions to the divergence

      character(len=*), parameter :: subname = '(formDiag_step2)'

      !-----------------------------------------------------------------
      ! integrate the momentum equation
      !-----------------------------------------------------------------

      strintx = c0
      strinty = c0
      
      ! Be careful: Drheo contains 4 terms for u and 4 terms for v.
      ! These 8 terms come from the surrounding T cells but are all
      ! refrerenced to the i,j (u point) :

      ! Drheo(i,j,1) corresponds to str(i,j,1)
      ! Drheo(i,j,2) corresponds to str(i+1,j,2)
      ! Drheo(i,j,3) corresponds to str(i,j+1,3)
      ! Drheo(i,j,4) corresponds to str(i+1,j+1,4))
      ! Drheo(i,j,5) corresponds to str(i,j,5)
      ! Drheo(i,j,6) corresponds to str(i,j+1,6)
      ! Drheo(i,j,7) corresponds to str(i+1,j,7)
      ! Drheo(i,j,8) corresponds to str(i+1,j+1,8))
      
      do ij = 1, icellu
         i = indxui(ij)
         j = indxuj(ij)

         ccaimp = umassdti(i,j) + vrel(i,j) * cosw + Cb(i,j) ! kg/m^2 s
         
         strintx = uarear(i,j)* &
             (Drheo(i,j,1) + Drheo(i,j,2) + Drheo(i,j,3) + Drheo(i,j,4))
         strinty = uarear(i,j)* &
             (Drheo(i,j,5) + Drheo(i,j,6) + Drheo(i,j,7) + Drheo(i,j,8))

         diagx(i,j) = ccaimp - strintx
         diagy(i,j) = ccaimp - strinty
      enddo                     ! ij

      end subroutine formDiag_step2

!=======================================================================

! Compute squared l^2 norm of a grid function (tpu,tpv)

      subroutine calc_L2norm_squared (nx_block, ny_block, &
                                      icellu  ,           &
                                      indxui  , indxuj  , &
                                      tpu     , tpv     , &
                                      L2norm)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellu                ! total count when iceumask is true

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
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
      
      do ij = 1, icellu
         i = indxui(ij)
         j = indxuj(ij)
         
         L2norm = L2norm + tpu(i,j)**2 + tpv(i,j)**2
      enddo ! ij
      
      end subroutine calc_L2norm_squared

!=======================================================================

! Convert a grid function (tpu,tpv) to a one dimensional vector

      subroutine arrays_to_vec (nx_block, ny_block  , &
                                nblocks , max_blocks, &
                                icellu  , ntot      , &
                                indxui  , indxuj    , &
                                tpu     , tpv       , &
                                outvec)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         nblocks,            & ! nb of blocks
         max_blocks,         & ! max nb of blocks
         ntot                  ! size of problem for Anderson

      integer (kind=int_kind), dimension (max_blocks), intent(in) :: &
         icellu

      integer (kind=int_kind), dimension (nx_block*ny_block, max_blocks), intent(in) :: &
         indxui  , & ! compressed index in i-direction
         indxuj      ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block, max_blocks), intent(in) :: &
         tpu     , & ! x-component of vector
         tpv         ! y-component of vector

      real (kind=dbl_kind), dimension (ntot), intent(out) :: &
         outvec      ! output 1D vector

      ! local variables

      integer (kind=int_kind) :: &
         i, j, iblk, tot, ij

      character(len=*), parameter :: subname = '(arrays_to_vec)'

      !-----------------------------------------------------------------
      ! form vector (converts from max_blocks arrays to single vector)
      !-----------------------------------------------------------------

      outvec(:) = c0
      tot = 0
      
      do iblk = 1, nblocks
         do ij = 1, icellu(iblk)
            i = indxui(ij, iblk)
            j = indxuj(ij, iblk)
            tot = tot + 1
            outvec(tot) = tpu(i, j, iblk)
            tot = tot + 1
            outvec(tot) = tpv(i, j, iblk)
         enddo
      enddo ! ij

      end subroutine arrays_to_vec

!=======================================================================

! Convert one dimensional vector to a grid function (tpu,tpv)

      subroutine vec_to_arrays (nx_block, ny_block  , &
                                nblocks , max_blocks, &
                                icellu  , ntot      , &
                                indxui  , indxuj    , &
                                invec   ,             &
                                tpu     , tpv)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         nblocks,            & ! nb of blocks
         max_blocks,         & ! max nb of blocks
         ntot                  ! size of problem for Anderson

      integer (kind=int_kind), dimension (max_blocks), intent(in) :: &
         icellu

      integer (kind=int_kind), dimension (nx_block*ny_block, max_blocks), intent(in) :: &
         indxui  , & ! compressed index in i-direction
         indxuj      ! compressed index in j-direction
         
      real (kind=dbl_kind), dimension (ntot), intent(in) :: &
         invec       ! input 1D vector

      real (kind=dbl_kind), dimension (nx_block,ny_block, max_blocks), intent(out) :: &
         tpu     , & ! x-component of vector
         tpv         ! y-component of vector

      ! local variables

      integer (kind=int_kind) :: &
         i, j, iblk, tot, ij

      character(len=*), parameter :: subname = '(vec_to_arrays)'

      !-----------------------------------------------------------------
      ! form arrays (converts from vector to the max_blocks arrays)
      !-----------------------------------------------------------------

      tpu(:,:,:) = c0
      tpv(:,:,:) = c0
      tot = 0
      
      do iblk = 1, nblocks
         do ij = 1, icellu(iblk)
            i = indxui(ij, iblk)
            j = indxuj(ij, iblk)
            tot = tot + 1
            tpu(i, j, iblk) = invec(tot)
            tot = tot + 1
            tpv(i, j, iblk) = invec(tot)
         enddo
      enddo! ij

      end subroutine vec_to_arrays

!=======================================================================

! Update Q and R factors after deletion of the 1st column of G_diff
!
! author: P. Blain ECCC
!
! adapted from :
! H. F. Walker, “Anderson Acceleration: Algorithms and Implementations”
!   [Online]. Available: https://users.wpi.edu/~walker/Papers/anderson_accn_algs_imps.pdf

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
      
      n = size(Q, 1)
      m = size(Q, 2)
      do i = 1, m-1
         temp = sqrt(R(i, i+1)**2 + R(i+1, i+1)**2)
         c = R(i  , i+1) / temp
         s = R(i+1, i+1) / temp
         R(i  , i+1) = temp
         R(i+1, i+1) = 0
         if (i < m-1) then
            do j = i+2, m
               temp      =  c*R(i, j) + s*R(i+1, j)
               R(i+1, j) = -s*R(i, j) + c*R(i+1, j)
               R(i  , j) = temp
            enddo
         endif
         do k = 1, n
            temp      =  c*Q(k, i) + s*Q(k, i+1);
            Q(k, i+1) = -s*Q(k, i) + c*Q(k, i+1);
            Q(k, i)   = temp
         enddo
      enddo
      R(:, 1:m-1) = R(:, 2:m)
      
      end subroutine qr_delete

!=======================================================================

! FGMRES: Flexible generalized minimum residual method (with restarts).
! Solves the linear system A x = b using GMRES with a varying (right) preconditioner
!
! authors: Stéphane Gaudreault, Abdessamad Qaddouri, Philippe Blain, ECCC

      subroutine fgmres (zetaD    ,           &
                         Cb       , vrel    , &
                         umassdti ,           &
                         halo_info_mask     , &
                         bx       , by      , &
                         diagx    , diagy   , &
                         tolerance, maxinner, &
                         maxouter ,           &
                         solx     , soly    , &
                         nbiter)

      use ice_boundary, only: ice_HaloUpdate
      use ice_domain, only: maskhalo_dyn, halo_info
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bound

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks,4), intent(in) :: &
         zetaD   ! zetaD = 2*zeta (viscous coefficient)

      real (kind=dbl_kind), dimension(nx_block, ny_block, max_blocks), intent(in) :: &
         vrel  , & ! coefficient for tauw
         Cb    , & ! seabed stress coefficient
         umassdti  ! mass of U-cell/dte (kg/m^2 s)

      type (ice_halo), intent(in) :: &
         halo_info_mask !  ghost cell update info for masked halo

      real (kind=dbl_kind), dimension(nx_block, ny_block, max_blocks), intent(in) :: &
         bx       , & ! Right hand side of the linear system (x components)
         by       , & ! Right hand side of the linear system (y components)
         diagx    , & ! Diagonal of the system matrix (x components)
         diagy        ! Diagonal of the system matrix (y components)

      real (kind=dbl_kind), intent(in) :: &
         tolerance   ! Tolerance to achieve. The algorithm terminates when the relative
                     ! residual is below tolerance

      integer (kind=int_kind), intent(in) :: &
         maxinner, & ! Restart the method every maxinner inner (Arnoldi) iterations
         maxouter    ! Maximum number of outer (restarts) iterations
                     ! Iteration will stop after maxinner*maxouter Arnoldi steps
                     ! even if the specified tolerance has not been achieved

      real (kind=dbl_kind), dimension(nx_block, ny_block, max_blocks), intent(inout) :: &
         solx     , & ! Initial guess on input, approximate solution on output (x components)
         soly         ! Initial guess on input, approximate solution on output (y components)

      integer (kind=int_kind), intent(out) :: &
         nbiter      ! Total number of Arnoldi iterations performed

      ! local variables

      integer (kind=int_kind) :: &
         iblk    , & ! block index
         ij      , & ! index for indx[t|u][i|j]
         i, j        ! grid indices

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         workspace_x , & ! work vector (x components)
         workspace_y     ! work vector (y components)

      real (kind=dbl_kind), dimension (max_blocks) :: &
         norm_squared   ! array to accumulate squared norm of grid function over blocks

      real (kind=dbl_kind), dimension(nx_block, ny_block, max_blocks, maxinner+1) :: &
         arnoldi_basis_x , & ! Arnoldi basis (x components)
         arnoldi_basis_y     ! Arnoldi basis (y components)

      real (kind=dbl_kind), dimension(nx_block, ny_block, max_blocks, maxinner) :: &
         orig_basis_x , & ! original basis (x components)
         orig_basis_y     ! original basis (y components)

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

      character (len=char_len) :: &
         precond_type ! type of preconditioner

      real (kind=dbl_kind) :: &
         relative_tolerance  ! relative_tolerance, i.e. tolerance*norm(initial residual)

      character(len=*), parameter :: subname = '(fgmres)'

      ! Here we go !

      ! Initialize
      outiter = 0
      nbiter = 0
      
      norm_squared = c0
      precond_type = precond
      
      ! Cells with no ice should be zero-initialized
      workspace_x = c0
      workspace_y = c0
      arnoldi_basis_x = c0
      arnoldi_basis_y = c0
      
      ! Residual of the initial iterate
      
      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks
         call matvec (nx_block               , ny_block             , &
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
                            indxui       (:,iblk), indxuj       (:,iblk), &
                            bx         (:,:,iblk), by         (:,:,iblk), &
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
            call calc_L2norm_squared(nx_block       , ny_block      , &
                                     icellu   (iblk),                 &
                                     indxui (:,iblk), indxuj(:,iblk), &
                                     arnoldi_basis_x(:,:,iblk, 1)   , &
                                     arnoldi_basis_y(:,:,iblk, 1)   , &
                                     norm_squared(iblk))

         enddo
         !$OMP END PARALLEL DO
         norm_residual = sqrt(global_sum(sum(norm_squared), distrb_info))
         
         if (my_task == master_task .and. monitor_fgmres) then
            write(nu_diag, '(a,i4,a,d26.16)') "monitor_fgmres: iter_fgmres= ", nbiter, &
                                              " fgmres_L2norm= ", norm_residual
         endif
         
         ! Current guess is a good enough solution TODO: reactivate and test this
         ! if (norm_residual < tolerance) then
         !    return
         ! end if
         
         ! Normalize the first Arnoldi vector
         inverse_norm = c1 / norm_residual
         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks
            do ij = 1, icellu(iblk)
               i = indxui(ij, iblk)
               j = indxuj(ij, iblk)

               arnoldi_basis_x(i, j, iblk, 1) = arnoldi_basis_x(i, j, iblk, 1) * inverse_norm
               arnoldi_basis_y(i, j, iblk, 1) = arnoldi_basis_y(i, j, iblk, 1) * inverse_norm
            enddo ! ij
         enddo
         !$OMP END PARALLEL DO
         
         if (outiter == 0) then
            relative_tolerance = tolerance * norm_residual
         end if
         
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
            call precondition(zetaD       ,                  &
                              Cb          , vrel           , &
                              umassdti    ,                  &
                              arnoldi_basis_x(:,:,:,initer), &
                              arnoldi_basis_y(:,:,:,initer), &
                              diagx       , diagy          , &
                              precond_type,                  &
                              workspace_x , workspace_y)
            orig_basis_x(:,:,:,initer) = workspace_x
            orig_basis_y(:,:,:,initer) = workspace_y
            
            ! Update workspace with boundary values
            call stack_velocity_field(workspace_x, workspace_y, fld2)
            call ice_timer_start(timer_bound)
            if (maskhalo_dyn) then
               call ice_HaloUpdate (fld2,               halo_info_mask, &
                                    field_loc_NEcorner, field_type_vector)
            else
               call ice_HaloUpdate (fld2,               halo_info, &
                                    field_loc_NEcorner, field_type_vector)
            endif
            call ice_timer_stop(timer_bound)
            call unstack_velocity_field(fld2, workspace_x, workspace_y)

            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks
               call matvec (nx_block               , ny_block             , &
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
            
            ! Orthogonalize the new vector
            call orthogonalize(ortho_type     , initer         , &
                               nextit         , maxinner       , &
                               arnoldi_basis_x, arnoldi_basis_y, &
                               hessenberg)
            
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
                  do ij = 1, icellu(iblk)
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
            
            ! Compute and apply new Givens rotation
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
            
            if (my_task == master_task .and. monitor_fgmres) then
               write(nu_diag, '(a,i4,a,d26.16)') "monitor_fgmres: iter_fgmres= ", nbiter, &
                                                 " fgmres_L2norm= ", norm_residual
            endif
            
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
               do ij = 1, icellu(iblk)
                  i = indxui(ij, iblk)
                  j = indxuj(ij, iblk)

                  solx(i, j, iblk) = solx(i, j, iblk) + t * orig_basis_x(i, j, iblk, it)
                  soly(i, j, iblk) = soly(i, j, iblk) + t * orig_basis_y(i, j, iblk, it)
               enddo ! ij
            enddo
            !$OMP END PARALLEL DO
         end do
         
         ! Increment outer loop counter and check for convergence
         outiter = outiter + 1
         if (norm_residual <= relative_tolerance .or. outiter >= maxouter) then
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
               do ij = 1, icellu(iblk)
                  i = indxui(ij, iblk)
                  j = indxuj(ij, iblk)

                  workspace_x(i, j, iblk) = workspace_x(i, j, iblk) + rhs_hess(it) * arnoldi_basis_x(i, j, iblk, it)
                  workspace_y(i, j, iblk) = workspace_x(i, j, iblk) + rhs_hess(it) * arnoldi_basis_y(i, j, iblk, it)
               enddo ! ij
            enddo
            !$OMP END PARALLEL DO
            arnoldi_basis_x(:,:,:,1) = workspace_x
            arnoldi_basis_y(:,:,:,1) = workspace_y
         end do
      end do ! end of outer (restarts) loop
      
      end subroutine fgmres

!=======================================================================

! PGMRES: Right-preconditioned generalized minimum residual method (with restarts).
! Solves the linear A x = b using GMRES with a right preconditioner
!
! authors: Stéphane Gaudreault, Abdessamad Qaddouri, Philippe Blain, ECCC

      subroutine pgmres (zetaD    ,           &
                         Cb       , vrel    , &
                         umassdti ,           &
                         bx       , by      , &
                         diagx    , diagy   , &
                         tolerance, maxinner, &
                         maxouter ,           &
                         solx     , soly    , &
                         nbiter)

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks,4), intent(in) :: &
         zetaD   ! zetaD = 2*zeta (viscous coefficient)

      real (kind=dbl_kind), dimension(nx_block, ny_block, max_blocks), intent(in) :: &
         vrel  , & ! coefficient for tauw
         Cb    , & ! seabed stress coefficient
         umassdti  ! mass of U-cell/dte (kg/m^2 s)

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
         maxinner, & ! Restart the method every maxinner inner (Arnoldi) iterations
         maxouter    ! Maximum number of outer (restarts) iterations
                     ! Iteration will stop after maxinner*maxouter Arnoldi steps
                     ! even if the specified tolerance has not been achieved

      real (kind=dbl_kind), dimension(nx_block, ny_block, max_blocks), intent(inout) :: &
         solx     , & ! Initial guess on input, approximate solution on output (x components)
         soly         ! Initial guess on input, approximate solution on output (y components)

      integer (kind=int_kind), intent(out) :: &
         nbiter      ! Total number of Arnoldi iterations performed

      ! local variables

      integer (kind=int_kind) :: &
         iblk    , & ! block index
         ij      , & ! index for indx[t|u][i|j]
         i, j        ! grid indices

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         workspace_x , & ! work vector (x components)
         workspace_y     ! work vector (y components)

      real (kind=dbl_kind), dimension (max_blocks) :: &
         norm_squared   ! array to accumulate squared norm of grid function over blocks

      real (kind=dbl_kind), dimension(nx_block, ny_block, max_blocks, maxinner+1) :: &
         arnoldi_basis_x , & ! Arnoldi basis (x components)
         arnoldi_basis_y     ! Arnoldi basis (y components)

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
         hessenberg          ! system matrix of the Hessenberg (least squares) system

      character(len=char_len) :: &
         precond_type    , & ! type of preconditioner
         ortho_type          ! type of orthogonalization

      real (kind=dbl_kind) :: &
         relative_tolerance  ! relative_tolerance, i.e. tolerance*norm(initial residual)

      character(len=*), parameter :: subname = '(pgmres)'
      
      ! Here we go !

      ! Initialize
      outiter = 0
      nbiter = 0
      
      norm_squared = c0
      precond_type = 'diag' ! Jacobi preconditioner
      ortho_type = 'cgs' ! classical gram-schmidt TODO: try with MGS
      
      ! Cells with no ice should be zero-initialized
      workspace_x = c0
      workspace_y = c0
      arnoldi_basis_x = c0
      arnoldi_basis_y = c0
      
      ! Residual of the initial iterate
      
      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks
         call matvec (nx_block               , ny_block             , &
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
                            indxui       (:,iblk), indxuj       (:,iblk), &
                            bx         (:,:,iblk), by         (:,:,iblk), &
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
            call calc_L2norm_squared(nx_block       , ny_block       , &
                                     icellu   (iblk),                  &
                                     indxui (:,iblk), indxuj(:, iblk), &
                                     arnoldi_basis_x(:,:,iblk, 1),     &
                                     arnoldi_basis_y(:,:,iblk, 1),     &
                                     norm_squared(iblk))

         enddo
         !$OMP END PARALLEL DO
         norm_residual = sqrt(global_sum(sum(norm_squared), distrb_info))
         
         if (my_task == master_task .and. monitor_pgmres) then
            write(nu_diag, '(a,i4,a,d26.16)') "monitor_pgmres: iter_pgmres= ", nbiter, &
                                              " pgmres_L2norm= ", norm_residual
         endif
         
         ! Current guess is a good enough solution
         ! if (norm_residual < tolerance) then
         !    return
         ! end if
         
         ! Normalize the first Arnoldi vector
         inverse_norm = c1 / norm_residual
         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks
            do ij = 1, icellu(iblk)
               i = indxui(ij, iblk)
               j = indxuj(ij, iblk)

               arnoldi_basis_x(i, j, iblk, 1) = arnoldi_basis_x(i, j, iblk, 1) * inverse_norm
               arnoldi_basis_y(i, j, iblk, 1) = arnoldi_basis_y(i, j, iblk, 1) * inverse_norm
            enddo ! ij
         enddo
         !$OMP END PARALLEL DO
         
         if (outiter == 0) then
            relative_tolerance = tolerance * norm_residual
         end if
         
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
            call precondition(zetaD       ,                  &
                              Cb          , vrel           , &
                              umassdti    ,                  &
                              arnoldi_basis_x(:,:,:,initer), &
                              arnoldi_basis_y(:,:,:,initer), &
                              diagx       , diagy          , &
                              precond_type,                  &
                              workspace_x , workspace_y)
            
            ! NOTE: halo updates for (workspace_x, workspace_y)
            ! are skipped here for efficiency since this is just a preconditioner
            
            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks
               call matvec (nx_block               , ny_block             , &
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
            
            ! Orthogonalize the new vector
            call orthogonalize(ortho_type     , initer         , &
                               nextit         , maxinner       , &
                               arnoldi_basis_x, arnoldi_basis_y, &
                               hessenberg)
            
            ! Compute norm of new Arnoldi vector and update Hessenberg matrix
            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks
               call calc_L2norm_squared(nx_block       , ny_block        , &
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
                  do ij = 1, icellu(iblk)
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
            
            ! Compute and apply new Givens rotation
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
            
            if (my_task == master_task .and. monitor_pgmres) then
               write(nu_diag, '(a,i4,a,d26.16)') "monitor_pgmres: iter_pgmres= ", nbiter, &
                                                 " pgmres_L2norm= ", norm_residual
            endif
            
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
               do ij = 1, icellu(iblk)
                  i = indxui(ij, iblk)
                  j = indxuj(ij, iblk)

                  workspace_x(i, j, iblk) = workspace_x(i, j, iblk) + t * arnoldi_basis_x(i, j, iblk, it)
                  workspace_y(i, j, iblk) = workspace_y(i, j, iblk) + t * arnoldi_basis_y(i, j, iblk, it)
               enddo ! ij
            enddo
            !$OMP END PARALLEL DO
         end do
         
         ! Call preconditioner
         call precondition(zetaD       ,              &
                           Cb          , vrel       , &
                           umassdti    ,              &
                           workspace_x , workspace_y, &
                           diagx       , diagy      , &
                           precond_type,              &
                           workspace_x , workspace_y)
         
         solx = solx + workspace_x
         soly = soly + workspace_y
         
         ! Increment outer loop counter and check for convergence
         outiter = outiter + 1
         if (norm_residual <= relative_tolerance .or. outiter >= maxouter) then
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
               do ij = 1, icellu(iblk)
                  i = indxui(ij, iblk)
                  j = indxuj(ij, iblk)

                  workspace_x(i, j, iblk) = workspace_x(i, j, iblk) + rhs_hess(it) * arnoldi_basis_x(i, j, iblk, it)
                  workspace_y(i, j, iblk) = workspace_x(i, j, iblk) + rhs_hess(it) * arnoldi_basis_y(i, j, iblk, it)
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

      subroutine precondition(zetaD       ,        &
                              Cb          , vrel , &
                              umassdti    ,        &
                              vx          , vy   , &
                              diagx       , diagy, &
                              precond_type,        &
                              wx          , wy)

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks,4), intent(in) :: &
         zetaD   ! zetaD = 2*zeta (viscous coefficient)

      real (kind=dbl_kind), dimension(nx_block, ny_block, max_blocks), intent(in) :: &
         vrel  , & ! coefficient for tauw
         Cb    , & ! seabed stress coefficient
         umassdti  ! mass of U-cell/dte (kg/m^2 s)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), intent(in) :: &
         vx       , & ! input vector (x components)
         vy           ! input vector (y components)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), intent(in) :: &
         diagx    , & ! diagonal of the system matrix (x components)
         diagy        ! diagonal of the system matrix (y components)

      character (len=char_len), intent(in) :: &
         precond_type ! type of preconditioner

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), intent(inout) :: &
         wx       , & ! preconditionned vector (x components)
         wy           ! preconditionned vector (y components)

      ! local variables

      integer (kind=int_kind) :: &
         iblk    , & ! block index
         ij      , & ! compressed index
         i, j        ! grid indices

      real (kind=dbl_kind) :: &
         tolerance   ! Tolerance for PGMRES

      integer (kind=int_kind) :: &
         maxinner    ! Restart parameter for PGMRES

      integer (kind=int_kind) :: &
         maxouter    ! Maximum number of outer iterations for PGMRES

      integer (kind=int_kind) :: &
         nbiter      ! Total number of iteration PGMRES performed

      character(len=*), parameter :: subname = '(precondition)'

      if     (precond_type == 'ident') then ! identity (no preconditioner)
         wx = vx
         wy = vy
      elseif (precond_type == 'diag') then ! Jacobi preconditioner (diagonal)
         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks
            do ij = 1, icellu(iblk)
               i = indxui(ij, iblk)
               j = indxuj(ij, iblk)

               wx(i,j,iblk) = vx(i,j,iblk) / diagx(i,j,iblk)
               wy(i,j,iblk) = vy(i,j,iblk) / diagy(i,j,iblk)
            enddo ! ij
         enddo
         !$OMP END PARALLEL DO
      elseif (precond_type == 'pgmres') then ! PGMRES (Jacobi-preconditioned GMRES)
         ! Initialize preconditioned vector to 0 ! TODO: try with wx = vx or vx/diagx
         wx = c0
         wy = c0
         tolerance = reltol_pgmres
         maxinner = dim_pgmres
         maxouter = maxits_pgmres
         call pgmres (zetaD,               &
                      Cb       , vrel    , &
                      umassdti ,           &
                      vx       , vy      , &
                      diagx    , diagy   , &
                      tolerance, maxinner, &
                      maxouter ,           &
                      wx       , wy      , &
                      nbiter)
      else
         call abort_ice(error_message='wrong preconditioner in ' // subname, &
            file=__FILE__, line=__LINE__)
      endif
      end subroutine precondition

!=======================================================================

! Generic routine to orthogonalize a vector (arnoldi_basis_[xy](:, :, :, nextit))
! against a set of vectors (arnoldi_basis_[xy](:, :, :, 1:initer))
!
! authors: Philippe Blain, ECCC

      subroutine orthogonalize(ortho_type     , initer         , &
                               nextit         , maxinner       , &
                               arnoldi_basis_x, arnoldi_basis_y, &
                               hessenberg)

      character(len=*), intent(in) :: &
         ortho_type ! type of orthogonalization

      integer (kind=int_kind), intent(in) :: &
         initer  , & ! inner (Arnoldi) loop counter
         nextit  , & ! nextit == initer+1
         maxinner    ! Restart the method every maxinner inner iterations

      real (kind=dbl_kind), dimension(nx_block, ny_block, max_blocks, maxinner+1), intent(inout) :: &
         arnoldi_basis_x , & ! arnoldi basis (x components)
         arnoldi_basis_y     ! arnoldi basis (y components)

      real (kind=dbl_kind), dimension(maxinner+1, maxinner), intent(inout) :: &
         hessenberg        ! system matrix of the Hessenberg (least squares) system

      ! local variables

      integer (kind=int_kind) :: &
         it      , & ! reusable loop counter
         iblk    , & ! block index
         ij      , & ! compressed index
         i, j        ! grid indices

      real (kind=dbl_kind), dimension (max_blocks) :: &
         local_dot      ! local array value to accumulate dot product of grid function over blocks

      real (kind=dbl_kind), dimension(maxinner) :: &
         dotprod_local  ! local array to accumulate several dot product computations

      character(len=*), parameter :: subname = '(orthogonalize)'

      if (trim(ortho_type) == 'cgs') then ! Classical Gram-Schmidt
         ! Classical Gram-Schmidt orthogonalisation process
         ! First loop of Gram-Schmidt (compute coefficients)
         dotprod_local = c0
         do it = 1, initer
            local_dot = c0
            
            !$OMP PARALLEL DO PRIVATE(iblk, ij, i, j)
            do iblk = 1, nblocks
               do ij = 1, icellu(iblk)
                  i = indxui(ij, iblk)
                  j = indxuj(ij, iblk)
                  
                  local_dot(iblk) = local_dot(iblk) + &
                                    (arnoldi_basis_x(i, j, iblk, it) * arnoldi_basis_x(i, j, iblk, nextit)) + &
                                    (arnoldi_basis_y(i, j, iblk, it) * arnoldi_basis_y(i, j, iblk, nextit))
               enddo ! ij
            enddo
            !$OMP END PARALLEL DO
            
            dotprod_local(it) = sum(local_dot)
         end do

         hessenberg(1:initer, initer) = global_allreduce_sum(dotprod_local(1:initer), distrb_info)

         ! Second loop of Gram-Schmidt (orthonormalize)
         do it = 1, initer
            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks
               do ij = 1, icellu(iblk)
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
      elseif (trim(ortho_type) == 'mgs') then ! Modified Gram-Schmidt
         ! Modified Gram-Schmidt orthogonalisation process
         do it = 1, initer
            local_dot = c0
            
            !$OMP PARALLEL DO PRIVATE(iblk, ij, i, j)
            do iblk = 1, nblocks
               do ij = 1, icellu(iblk)
                  i = indxui(ij, iblk)
                  j = indxuj(ij, iblk)
                  
                  local_dot(iblk) = local_dot(iblk) + &
                                    (arnoldi_basis_x(i, j, iblk, it) * arnoldi_basis_x(i, j, iblk, nextit)) + &
                                    (arnoldi_basis_y(i, j, iblk, it) * arnoldi_basis_y(i, j, iblk, nextit))
               enddo ! ij
            enddo
            !$OMP END PARALLEL DO
            
            hessenberg(it,initer) = global_sum(sum(local_dot), distrb_info)
            
            !$OMP PARALLEL DO PRIVATE(iblk, ij, i, j)
            do iblk = 1, nblocks
               do ij = 1, icellu(iblk)
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
      else
         call abort_ice(error_message='wrong orthonalization in ' // subname, &
         file=__FILE__, line=__LINE__)
      endif
      
   end subroutine orthogonalize

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
