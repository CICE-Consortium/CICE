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
      use ice_constants, only: field_loc_center, field_loc_NEcorner, &
          field_type_scalar, field_type_vector
      use ice_constants, only: c0, c4, p027, p055, p111, p166, &
          p2, p222, p25, p333, p5, c1
      use ice_dyn_shared, only: evp_prep1, evp_prep2, evp_finish, &
          yield_curve, ecci, cosw, sinw, fcor_blk, uvel_init,  &
          vvel_init, basal_stress_coeff, basalstress, Ktens
      use ice_fileunits, only: nu_diag
      use ice_exit, only: abort_ice
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_ice_strength, icepack_query_parameters
#ifdef CICE_IN_NEMO
      use icepack_intfc, only: calc_strair
#endif

      implicit none
      private
      public :: imp_solver

!=======================================================================

      contains

!=======================================================================

! Viscous-plastic dynamics driver
!
#ifdef CICE_IN_NEMO
! Wind stress is set during this routine from the values supplied
! via NEMO (unless calc_strair is true).  These values are supplied 
! rotated on u grid and multiplied by aice.  strairxT = 0 in this 
! case so operations in evp_prep1 are pointless but carried out to 
! minimise code changes.
#endif
!
! author: JF Lemieux, F. Dupont and A. Qaddouri, ECCC

      subroutine imp_solver (dt)

      use ice_arrays_column, only: Cdn_ocn
      use ice_boundary, only: ice_halo, ice_HaloMask, ice_HaloUpdate, &
          ice_HaloDestroy, ice_HaloUpdate_stress
      use ice_blocks, only: block, get_block, nx_block, ny_block
      use ice_domain, only: nblocks, blocks_ice, halo_info, maskhalo_dyn
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
          tarear, uarear, tinyarea, to_ugrid, t2ugrid_vector, u2tgrid_vector, &
          grid_type
      use ice_state, only: aice, vice, vsno, uvel, vvel, divu, shear, &
          aice_init, aice0, aicen, vicen, strength
      use ice_timers, only: timer_dynamics, timer_bound, &
          ice_timer_start, ice_timer_stop

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) :: & 
         kOL            , & ! outer loop iteration
         kmax           , & ! jfl put in namelist
         ntot           , & ! size of problem for fgmres (for given cpu)
         icode          , & ! for fgmres
         iconvNL        , & ! code for NL convergence criterion
         iout           , & ! for printing fgmres info
         its            , & ! iteration nb for fgmres
         ischmi         , & ! Quesse ca!?!?! jfl
         maxits         , & ! max nb of iteration for fgmres
         fgmres_its     , & ! final nb of fgmres_its
         im_fgmres      , & ! for size of Krylov subspace
         precond        , & ! 1: identity, 2: diagonal
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, ij

      integer (kind=int_kind), dimension(max_blocks) :: & 
         icellt   , & ! no. of cells where icetmask = 1
         icellu       ! no. of cells where iceumask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block, max_blocks) :: &
         indxti   , & ! compressed index in i-direction
         indxtj   , & ! compressed index in j-direction
         indxui   , & ! compressed index in i-direction
         indxuj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         tmass    , & ! total mass of ice and snow (kg/m^2)
         waterx   , & ! for ocean stress calculation, x (m/s)
         watery   , & ! for ocean stress calculation, y (m/s)
         forcex   , & ! work array: combined atm stress and ocn tilt, x
         forcey   , & ! work array: combined atm stress and ocn tilt, y
         bxfix    , & ! part of bx that is constant during Picard 
         byfix    , & ! part of by that is constant during Picard 
         bx       , & ! b vector
         by       , & ! b vector
         Au       , & ! matvec, Fx = Au - bx
         Av       , & ! matvec, Fy = Av - by
         Diagu    , & ! Diagonal (u component) of the matrix A
         Diagv    , & ! Diagonal (v component) of the matrix A
         Fx       , & ! x residual vector, Fx = Au - bx 
         Fy       , & ! y residual vector, Fy = Av - by 
         uprev_k  , & ! uvel at previous Picard iteration
         vprev_k  , & ! vvel at previous Picard iteration
         ulin     , & ! uvel to linearize vrel
         vlin     , & ! vvel to linearize vrel
         vrel     , & ! coeff for tauw 
         Cb       , & ! seabed stress coeff
         aiu      , & ! ice fraction on u-grid
         umass    , & ! total mass of ice and snow (u grid)
         umassdti     ! mass of U-cell/dte (kg/m^2 s)
         
      real (kind=dbl_kind), allocatable :: fld2(:,:,:,:)
      
      real (kind=dbl_kind), allocatable :: bvec(:), sol(:), diagvec(:), wk11(:), wk22(:)
      real (kind=dbl_kind), allocatable :: vv(:,:), ww(:,:)
      
      real (kind=dbl_kind), dimension (max_blocks) :: L2norm
      real (kind=dbl_kind) :: conv, gamma, gammaNL, tolNL, krelax

      real (kind=dbl_kind), dimension(nx_block,ny_block,8):: &
         strtmp,    & ! stress combinations for momentum equation !JFL CHECK PAS SUR QUE OK
         stPrtmp,   & ! doit etre (nx_block,ny_block,max_blocks,8)???? PAs besoin des 3? reuse?
         Dstrtmp
         
      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks,4):: &
         zetaD      ! zetaD = 2zeta (viscous coeff)

      integer (kind=int_kind), dimension (nx_block,ny_block,max_blocks) :: &
         icetmask, &  ! ice extent mask (T-cell)
         halomask     ! generic halo mask

      type (ice_halo) :: &
         halo_info_mask !  ghost cell update info for masked halo

      type (block) :: &
         this_block           ! block information for current block
      
      call ice_timer_start(timer_dynamics) ! dynamics

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      allocate(fld2(nx_block,ny_block,2,max_blocks))
      
      !-----------------------------------------------------------------
      ! Define a few things for FGMRES and Picard solver
      !-----------------------------------------------------------------
      
      im_fgmres = 50 
      maxits = 50    
      kmax=2
      gammaNL=1e-2_dbl_kind !linear stopping criterion: gamma*(res_ini)
      gamma=1e-6_dbl_kind   !nonlinear stopping criterion:
      iconvNL=0 ! equals 1 when NL convergence is reached
      krelax=c1
      precond=2 ! 1: identity, 2: diagonal

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
      ! preparation for dynamics JFL change names of evp_prep1 and 2
      !-----------------------------------------------------------------

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         call evp_prep1 (nx_block,           ny_block,           & 
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

#ifdef CICE_IN_NEMO
      !----------------------------------------------------------------
      ! Set wind stress to values supplied via NEMO
      ! This wind stress is rotated on u grid and multiplied by aice
      !----------------------------------------------------------------
      if (.not. calc_strair) then       
         strairx(:,:,:) = strax(:,:,:)
         strairy(:,:,:) = stray(:,:,:)
      else
#endif
      call t2ugrid_vector(strairx)
      call t2ugrid_vector(strairy)
#ifdef CICE_IN_NEMO
      endif      
#endif

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

         call evp_prep2 (nx_block,             ny_block,             & 
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

      !-----------------------------------------------------------------
      ! calc size of problem (ntot) and allocate arrays and vectors
      !-----------------------------------------------------------------
      
      ntot=0
      do iblk = 1, nblocks
        ntot = ntot + icellu(iblk)      
      enddo
      ntot = 2*ntot ! times 2 because of u and v
      
      allocate(bvec(ntot), sol(ntot), diagvec(ntot), wk11(ntot), wk22(ntot))
      allocate(vv(ntot,im_fgmres+1), ww(ntot,im_fgmres))
      
      !-----------------------------------------------------------------
      
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message="subname", &
         file=__FILE__, line=__LINE__)

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (strength,           halo_info, &
                           field_loc_center,   field_type_scalar)
      ! velocities may have changed in evp_prep2
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
       do iblk = 1, nblocks
         call basal_stress_coeff (nx_block,         ny_block,       &
                                  icellu  (iblk),                   &
                                  indxui(:,iblk),   indxuj(:,iblk), &
                                  vice(:,:,iblk),   aice(:,:,iblk), &
                                  hwater(:,:,iblk), Tbu(:,:,iblk))
       enddo                           
      endif
      
      do kOL = 1,kmax        ! outer loop 
      
      !-----------------------------------------------------------------
      ! Calc zetaD, vrel, Cb and vrel = f(uprev_k, vprev_k)
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk,strtmp)
         do iblk = 1, nblocks

	    if (kOL .eq. 1) then
	      ulin(:,:,iblk) = uvel(:,:,iblk)
	      vlin(:,:,iblk) = vvel(:,:,iblk)
	    else
	      ulin(:,:,iblk) = p5*uprev_k(:,:,iblk) + p5*uvel(:,:,iblk)
	      vlin(:,:,iblk) = p5*vprev_k(:,:,iblk) + p5*vvel(:,:,iblk)
	    endif
         
            uprev_k(:,:,iblk) = uvel(:,:,iblk)
            vprev_k(:,:,iblk) = vvel(:,:,iblk)
            
            call calc_zeta_Pr (nx_block           , ny_block,           &
                               kOL                , icellt(iblk),       & 
                               indxti   (:,iblk)  , indxtj(:,iblk),     & 
                               uprev_k  (:,:,iblk), vprev_k (:,:,iblk), & 
                               dxt      (:,:,iblk), dyt   (:,:,iblk),   & 
                               dxhy     (:,:,iblk), dyhx  (:,:,iblk),   & 
                               cxp      (:,:,iblk), cyp   (:,:,iblk),   & 
                               cxm      (:,:,iblk), cym   (:,:,iblk),   & 
                               tarear   (:,:,iblk), tinyarea (:,:,iblk),& 
                               strength (:,:,iblk), zetaD (:,:,iblk,:) ,&
                               stPrtmp  (:,:,:) )                      
            
            call calc_vrel_Cb (nx_block           , ny_block,           &
                               icellu       (iblk), Cdn_ocn (:,:,iblk), & 
                               indxui     (:,iblk), indxuj    (:,iblk), &
                               kOL                ,                     &
                               aiu      (:,:,iblk), Tbu     (:,:,iblk), &
                               uocn     (:,:,iblk), vocn    (:,:,iblk), &     
                               ulin     (:,:,iblk), vlin    (:,:,iblk), & 
                               vrel     (:,:,iblk), Cb      (:,:,iblk))

!     prepare b vector (RHS)                                                
            call calc_bvec (nx_block           , ny_block,           &
                            icellu       (iblk),                     & 
                            indxui     (:,iblk), indxuj    (:,iblk), &
                            kOL                , Cdn_ocn (:,:,iblk), &
                            aiu      (:,:,iblk), uarear  (:,:,iblk), & 
                            uocn     (:,:,iblk), vocn    (:,:,iblk), &     
                            waterx   (:,:,iblk), watery  (:,:,iblk), & 
                            ulin     (:,:,iblk), vlin    (:,:,iblk), & 
                            bxfix    (:,:,iblk), byfix   (:,:,iblk), &
                            bx       (:,:,iblk), by      (:,:,iblk), &
                            stPrtmp  (:,:,:))

!     prepare precond matrix                                                        
           call formDiag_step1 (nx_block           , ny_block,       &
                                icellt       (iblk), 1             , & ! for u comp
                                indxti     (:,iblk), indxtj(:,iblk), &
                                dxt      (:,:,iblk), dyt (:,:,iblk), & 
                                dxhy     (:,:,iblk), dyhx(:,:,iblk), & 
                                cxp      (:,:,iblk), cyp (:,:,iblk), & 
                                cxm      (:,:,iblk), cym (:,:,iblk), & 
                                zetaD (:,:,iblk,:) , Dstrtmp (:,:,:) )
                                
           call formDiag_step1 (nx_block           , ny_block,       &
                                icellt       (iblk), 2             , & ! for v comp
                                indxti     (:,iblk), indxtj(:,iblk), &
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

!-----------------------------------------------------------------------
!     prep F G M R E S 
!-----------------------------------------------------------------------                             
         
      icode  = 0
      iout   = 2 !0: nothing printed, 1: 1st ite only, 2: all iterations
!      its    = 0 
      ischmi = 0 
         
         ! form b vector from matrices (nblocks matrices)      
         call arrays_to_vec (nx_block, ny_block, nblocks,    &
                             max_blocks, icellu (:), ntot,   & 
                             indxui      (:,:), indxuj(:,:), &
                             bx        (:,:,:), by  (:,:,:), &
                             bvec(:))
         ! form sol vector for fgmres (sol is iniguess at the beginning)        
         call arrays_to_vec (nx_block, ny_block, nblocks,      &
                             max_blocks, icellu (:), ntot,   &  
                             indxui    (:,:), indxuj(:,:),     &
                             uprev_k (:,:,:), vprev_k (:,:,:), &
                             sol(:))
                             
         ! form matrix diagonal as a vector from Diagu and Diagv arrays      
         call arrays_to_vec (nx_block, ny_block, nblocks,    &
                             max_blocks, icellu (:), ntot,   & 
                             indxui      (:,:), indxuj(:,:), &
                             Diagu     (:,:,:), Diagv(:,:,:),&
                             diagvec(:))                             

!-----------------------------------------------------------------------
!     F G M R E S   L O O P
!-----------------------------------------------------------------------
 1    continue
!-----------------------------------------------------------------------

      !call fgmres2( ntot,im_fgmres,bvec,sol,ischmi,vv,ww,wk11,wk22, &
      !                     sol_eps, maxits,its,conv,icode )
                           
      call fgmres (ntot,im_fgmres,bvec,sol,its,vv,ww,wk11,wk22, &
                   gamma, gammaNL, tolNL, maxits,iout,icode,iconvNL,fgmres_its,kOL)                     

      if (iconvNL .eq. 1) exit             
                   
      if (icode == 1) then

!         if (sol2D_precond_S == 'JACOBI')   then
!               call pre_jacobi2D ( wk22,wk11,Prec_xevec_8,niloc,njloc,&
!                                   F_nk,Prec_ai_8,Prec_bi_8,Prec_ci_8 )
!         else
!            call dcopy (nloc, wk11, 1, wk22, 1) ! precond=identity
!         endif

         if (precond .eq. 1) then

           wk22(:)=wk11(:) ! precond=identity
           
         elseif (precond .eq. 2) then ! use diagonal of A for precond step
          
           call precond_diag (ntot,            & 
                              diagvec (:),     &
                              wk11 (:), wk22 (:) )
         endif
         
         goto 1

      else

         if (icode >= 2) then

!            if (Lun_debug_L.and.print_conv_L) write(lun_out, 199) conv,its
!            call sol_matvec ( wk22, wk11, Minx, Maxx, Miny, Maxy, &
!                           nil,njl, F_nk, minx1,maxx1,minx2,maxx2 )

         call vec_to_arrays (nx_block, ny_block, nblocks,      &
                             max_blocks, icellu (:), ntot,     & 
                             indxui    (:,:), indxuj(:,:),     &
                             wk11 (:),                         &
                             uvel (:,:,:), vvel (:,:,:))    

         !$OMP PARALLEL DO PRIVATE(iblk,strtmp)
         do iblk = 1, nblocks                                  
         
            call stress_vp (nx_block,             ny_block,             & 
                            kOL,                  icellt(iblk),         & 
                            indxti      (:,iblk), indxtj      (:,iblk), & 
                            uvel      (:,:,iblk), vvel      (:,:,iblk), &     
                            dxt       (:,:,iblk), dyt       (:,:,iblk), & 
                            dxhy      (:,:,iblk), dyhx      (:,:,iblk), & 
                            cxp       (:,:,iblk), cyp       (:,:,iblk), & 
                            cxm       (:,:,iblk), cym       (:,:,iblk), & 
                            tarear    (:,:,iblk), tinyarea  (:,:,iblk), & 
                            zetaD     (:,:,iblk,:),                     &
                            shear     (:,:,iblk), divu      (:,:,iblk), & 
                            rdg_conv  (:,:,iblk), rdg_shear (:,:,iblk), & 
                            strtmp    (:,:,:))                                                             
                               
            call matvec (nx_block           , ny_block,           &
                         icellu       (iblk),                     & 
                         indxui     (:,iblk), indxuj    (:,iblk), &
                         kOL                ,                     &
                         aiu      (:,:,iblk), strtmp  (:,:,:),    &
                         vrel     (:,:,iblk),                     &
                         umassdti (:,:,iblk), fm      (:,:,iblk), & 
                         uarear   (:,:,iblk), Cb      (:,:,iblk), & 
                         uvel     (:,:,iblk), vvel    (:,:,iblk), &
                         Au       (:,:,iblk), Av      (:,:,iblk))
                         
         enddo
         !$OMP END PARALLEL DO 
        
         ! form wk2 from Au and Av arrays        
         call arrays_to_vec (nx_block, ny_block, nblocks,      &
                             max_blocks, icellu (:), ntot,     & 
                             indxui    (:,:), indxuj(:,:),     &
                             Au      (:,:,:), Av    (:,:,:),   &
                             wk22(:))    

            goto 1

         endif

      endif

! 199  format (3x,'Iterative FGMRES solver convergence criteria: ',1pe14.7,' at iteration', i3)

!     deallocate (wk11,wk22,rhs1,sol1,vv_8,ww_8)         
         
!            call calc_L2norm (nx_block           , ny_block,          &
!                             icellu       (iblk),                     & 
!                             indxui     (:,iblk), indxuj    (:,iblk), &
!                             uvel     (:,:,iblk), vvel    (:,:,iblk))                                                     
           
!            call residual_vec (nx_block           , ny_block,           &
!                               icellu       (iblk),                     & 
!                               indxui     (:,iblk), indxuj    (:,iblk), &
!                               bx       (:,:,iblk), by      (:,:,iblk), &
!                               Au       (:,:,iblk), Av      (:,:,iblk), &
!                               Fx       (:,:,iblk), Fy      (:,:,iblk), &
!                               L2norm(iblk))
  
!            call precondD  (nx_block,             ny_block,             & 
!                            kOL                 , icellt(iblk),         & 
!                            indxti      (:,iblk), indxtj      (:,iblk), & 
!                            dxt       (:,:,iblk), dyt       (:,:,iblk), & 
!                            dxhy      (:,:,iblk), dyhx      (:,:,iblk), & 
!                            cxp       (:,:,iblk), cyp       (:,:,iblk), & 
!                            cxm       (:,:,iblk), cym       (:,:,iblk), & 
!                            uarear    (:,:,iblk),                       &
!                            vrel      (:,:,iblk), Cb        (:,:,iblk), &
!                            umassdti  (:,:,iblk), zetaD   (:,:,iblk,:), &
!                            Diagu     (:,:,iblk), Diagv   (:,:,iblk))

!-----------------------------------------------------------------------
!     Put vector sol in uvel and vvel arrays
!-----------------------------------------------------------------------

         call vec_to_arrays (nx_block, ny_block, nblocks,      &
                             max_blocks, icellu (:), ntot,     & 
                             indxui    (:,:), indxuj(:,:),     &
                             sol (:),                          &
                             uvel (:,:,:), vvel (:,:,:))    

         !$OMP PARALLEL DO PRIVATE(iblk,strtmp)
         do iblk = 1, nblocks
              uvel(:,:,iblk) = (c1-krelax)*uprev_k(:,:,iblk) + krelax*uvel(:,:,iblk)
              vvel(:,:,iblk) = (c1-krelax)*vprev_k(:,:,iblk) + krelax*vvel(:,:,iblk)
         enddo
         !$OMP END PARALLEL DO  
                             
         !$OMP PARALLEL DO PRIVATE(iblk,strtmp)
         do iblk = 1, nblocks                             
                            
            ! load velocity into array for boundary updates
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

         ! unload
         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks
            uvel(:,:,iblk) = fld2(:,:,1,iblk)
            vvel(:,:,iblk) = fld2(:,:,2,iblk)
         enddo
         !$OMP END PARALLEL DO
         
      enddo                     ! outer loop

      deallocate(bvec, sol, diagvec, wk11, wk22, vv, ww)
      deallocate(fld2)
      if (maskhalo_dyn) call ice_HaloDestroy(halo_info_mask)

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

         call evp_finish                               & 
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

! Computes the viscous coefficients (in fact zetaD=2*zeta) and dPr/dx. 

      subroutine calc_zeta_Pr  (nx_block,   ny_block,   & 
                                kOL,        icellt,     & 
                                indxti,     indxtj,     & 
                                uvel,       vvel,       & 
                                dxt,        dyt,        & 
                                dxhy,       dyhx,       & 
                                cxp,        cyp,        & 
                                cxm,        cym,        & 
                                tarear,     tinyarea,   & 
                                strength,   zetaD,      &
                                stPr)

      integer (kind=int_kind), intent(in) :: & 
         nx_block, ny_block, & ! block dimensions
         kOL               , & ! subcycling step
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
         tarear   , & ! 1/tarea
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

      subroutine stress_vp (nx_block,   ny_block,   & 
                            kOL,        icellt,     & 
                            indxti,     indxtj,     & 
                            uvel,       vvel,       & 
                            dxt,        dyt,        & 
                            dxhy,       dyhx,       & 
                            cxp,        cyp,        & 
                            cxm,        cym,        & 
                            tarear,     tinyarea,   & 
                            zetaD,                  & 
                            shear,      divu,       & 
                            rdg_conv,   rdg_shear,  & 
                            str )

      integer (kind=int_kind), intent(in) :: & 
         nx_block, ny_block, & ! block dimensions
         kOL               , & ! subcycling step
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
         cxm      , & ! 0.5*HTN - 1.5*HTN
         tarear   , & ! 1/tarea
         tinyarea     ! puny*tarea
         
      real (kind=dbl_kind), dimension(nx_block,ny_block,4), & 
         intent(in) :: &
         zetaD          ! 2*zeta   

      real (kind=dbl_kind), dimension (nx_block,ny_block), & 
         intent(inout) :: &
         shear    , & ! strain rate II component (1/s)
         divu     , & ! strain rate I component, velocity divergence (1/s)
         rdg_conv , & ! convergence term for ridging (1/s)
         rdg_shear    ! shear term for ridging (1/s)

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
        puny                                      , & ! puny
        ssigpn, ssigps, ssigpe, ssigpw            , &
        ssigmn, ssigms, ssigme, ssigmw            , &
        ssig12n, ssig12s, ssig12e, ssig12w        , &
        ssigp1, ssigp2, ssigm1, ssigm2, ssig121, ssig122, &
        csigpne, csigpnw, csigpse, csigpsw        , &
        csigmne, csigmnw, csigmse, csigmsw        , &
        csig12ne, csig12nw, csig12se, csig12sw    , &
        str12ew, str12we, str12ns, str12sn        , &
        strp_tmp, strm_tmp, tmp
        
      real (kind=dbl_kind) :: &
         stressp_1, stressp_2, stressp_3, stressp_4 , & ! sigma11+sigma22 (without Pr)
         stressm_1, stressm_2, stressm_3, stressm_4 , & ! sigma11-sigma22
         stress12_1,stress12_2,stress12_3,stress12_4    ! sigma12        

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
      ! on last subcycle, save quantities for mechanical redistribution
      !-----------------------------------------------------------------
         if (kOL == 100) then ! jfl MODIF
            divu(i,j) = p25*(divune + divunw + divuse + divusw) * tarear(i,j)
            tmp = p25*(Deltane + Deltanw + Deltase + Deltasw)   * tarear(i,j)
            rdg_conv(i,j)  = -min(divu(i,j),c0)
            rdg_shear(i,j) = p5*(tmp-abs(divu(i,j))) 

            ! diagnostic only
            ! shear = sqrt(tension**2 + shearing**2)
            shear(i,j) = p25*tarear(i,j)*sqrt( &
                 (tensionne + tensionnw + tensionse + tensionsw)**2 &
                +  (shearne +   shearnw +   shearse +   shearsw)**2)

         endif

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

      end subroutine stress_vp      
      
      
!=======================================================================

! Computes the VP stress (as diagnostic)

      subroutine Diagstress_vp (nx_block,   ny_block,   & 
                            kOL,        icellt,     & 
                            indxti,     indxtj,     & 
                            uvel,       vvel,       & 
                            dxt,        dyt,        & 
                            dxhy,       dyhx,       & 
                            cxp,        cyp,        & 
                            cxm,        cym,        & 
                            tarear,     tinyarea,   & 
                            zetaD,                  & 
                            stressp_1,  stressp_2,  & 
                            stressp_3,  stressp_4,  & 
                            stressm_1,  stressm_2,  & 
                            stressm_3,  stressm_4,  & 
                            stress12_1, stress12_2, & 
                            stress12_3, stress12_4, & 
                            shear,      divu,       & 
                            rdg_conv,   rdg_shear,  & 
                            str )

      integer (kind=int_kind), intent(in) :: & 
         nx_block, ny_block, & ! block dimensions
         kOL               , & ! subcycling step
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
         cxm      , & ! 0.5*HTN - 1.5*HTN
         tarear   , & ! 1/tarea
         tinyarea     ! puny*tarea
         
      real (kind=dbl_kind), dimension(nx_block,ny_block,4), & 
         intent(in) :: &
         zetaD          ! 2*zeta   

      real (kind=dbl_kind), dimension (nx_block,ny_block), & 
         intent(inout) :: &
         stressp_1, stressp_2, stressp_3, stressp_4 , & ! sigma11+sigma22
         stressm_1, stressm_2, stressm_3, stressm_4 , & ! sigma11-sigma22
         stress12_1,stress12_2,stress12_3,stress12_4    ! sigma12

      real (kind=dbl_kind), dimension (nx_block,ny_block), & 
         intent(inout) :: &
         shear    , & ! strain rate II component (1/s)
         divu     , & ! strain rate I component, velocity divergence (1/s)
         rdg_conv , & ! convergence term for ridging (1/s)
         rdg_shear    ! shear term for ridging (1/s)

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
        puny                                      , & ! puny
        ssigpn, ssigps, ssigpe, ssigpw            , &
        ssigmn, ssigms, ssigme, ssigmw            , &
        ssig12n, ssig12s, ssig12e, ssig12w        , &
        ssigp1, ssigp2, ssigm1, ssigm2, ssig121, ssig122, &
        csigpne, csigpnw, csigpse, csigpsw        , &
        csigmne, csigmnw, csigmse, csigmsw        , &
        csig12ne, csig12nw, csig12se, csig12sw    , &
        str12ew, str12we, str12ns, str12sn        , &
        strp_tmp, strm_tmp, tmp

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
      ! on last subcycle, save quantities for mechanical redistribution
      !-----------------------------------------------------------------
         if (kOL == 100) then ! jfl MODIF
            divu(i,j) = p25*(divune + divunw + divuse + divusw) * tarear(i,j)
            tmp = p25*(Deltane + Deltanw + Deltase + Deltasw)   * tarear(i,j)
            rdg_conv(i,j)  = -min(divu(i,j),c0)
            rdg_shear(i,j) = p5*(tmp-abs(divu(i,j))) 

            ! diagnostic only
            ! shear = sqrt(tension**2 + shearing**2)
            shear(i,j) = p25*tarear(i,j)*sqrt( &
                 (tensionne + tensionnw + tensionse + tensionsw)**2 &
                +  (shearne +   shearnw +   shearse +   shearsw)**2)

         endif

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

      end subroutine Diagstress_vp
      
!=======================================================================

      subroutine calc_vrel_Cb (nx_block,   ny_block, &
                               icellu,     Cw,       &
                               indxui,     indxuj,   &
                               kOL,                  &
                               aiu,        Tbu,      &
                               uocn,       vocn,     &
                               uvel,       vvel,     &
                               vrel,       Cb)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellu,             & ! total count when iceumask is true
         kOL                   ! outer loop iteration

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

      subroutine matvec (nx_block,   ny_block, &
                         icellu,               &
                         indxui,     indxuj,   &
                         kOL,                  &
                         aiu,        str,      &
                         vrel,                 &
                         umassdti,   fm,       &
                         uarear,     Cb,       &
                         uvel,       vvel,     &
                         Au,         Av)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellu,             & ! total count when iceumask is true
         kOL                   ! outer loop iteration

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxui  , & ! compressed index in i-direction
         indxuj      ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         vrel,     & ! coefficient for tauw
         Cb,       & ! coefficient for basal stress
         aiu     , & ! ice fraction on u-grid
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

      !-----------------------------------------------------------------
      ! integrate the momentum equation
      !-----------------------------------------------------------------

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message="subname", &
         file=__FILE__, line=__LINE__)

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
         
      ! calculate basal stress component for outputs ! jfl move this
!         if (ksub == ndte) then ! on last subcycling iteration
!          if ( basalstress ) then
!           taubx(i,j) = -uvel(i,j)*Tbu(i,j) / (sqrt(uold**2 + vold**2) + u0)
!           tauby(i,j) = -vvel(i,j)*Tbu(i,j) / (sqrt(uold**2 + vold**2) + u0)
!          endif
!         endif

      enddo                     ! ij

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
                       kOL,        Cw,       &
                       aiu,        uarear,   &
                       uocn,       vocn,     &
                       waterx,     watery,   &
                       uvel,       vvel,     &
                       bxfix,      byfix,    &
                       bx,         by,       &
                       stPr)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellu,             & ! total count when iceumask is true
         kOL                   ! outer loop iteration

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
         vocn        ! ocean current, y-direction (m/s)
         
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
         vrel              , & ! relative ice-ocean velocity
         utp, vtp          , & ! utp = uvel, vtp = vvel !jfl needed?
         taux, tauy        , & ! part of ocean stress term
         strintx, strinty  , & ! divergence of the internal stress tensor (only Pr part)
         rhow                  !
         
      !-----------------------------------------------------------------
      ! calc b vector
      !-----------------------------------------------------------------

      !JFL vrel could be sent here (already calc before...
      
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
         vrel = aiu(i,j)*rhow*Cw(i,j)*sqrt((uocn(i,j) - utp)**2 + &
                                           (vocn(i,j) - vtp)**2)  ! m/s
         ! ice/ocean stress
         taux = vrel*waterx(i,j) ! NOTE this is not the entire
         tauy = vrel*watery(i,j) ! ocn stress term
         
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
                               L2normtp)

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
         
      real (kind=dbl_kind), intent(inout) :: &
         L2normtp    ! (L2norm)^2
      
      integer (kind=int_kind) :: &
         i, j, ij

      !-----------------------------------------------------------------
      ! calc b vector
      !-----------------------------------------------------------------

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message="subname", &
         file=__FILE__, line=__LINE__)

      L2normtp=c0
         
      do ij =1, icellu
         i = indxui(ij)
         j = indxuj(ij)

         Fx(i,j) = Au(i,j) - bx(i,j)
         Fy(i,j) = Av(i,j) - by(i,j)
         L2normtp = L2normtp + Fx(i,j)**2 + Fy(i,j)**2
!         Fres(2*ij-1) = Au(i,j) - bx(i,j)
!         Fres(2*ij)   = Av(i,j) - by(i,j)
         
      enddo                     ! ij
      
!      do ij = 1, ntot
      
!	Ftp(ij) = Aw(ij) - bvec(ij)
      
!      enddo
      
!       L2norm = sqrt(DOT_PRODUCT(Fres,Fres))
!       print *, 'ici L2norm', sqrt(L2normtp)

      end subroutine residual_vec
      
!=======================================================================

      subroutine formDiag_step1  (nx_block,   ny_block,   & 
                                  icellt,     velcode,    & 
                                  indxti,     indxtj,     & 
                                  dxt,        dyt,        & 
                                  dxhy,       dyhx,       & 
                                  cxp,        cyp,        & 
                                  cxm,        cym,        & 
                                  zetaD,      Dstr )

      integer (kind=int_kind), intent(in) :: & 
         nx_block, ny_block, & ! block dimensions
         icellt,             & ! no. of cells where icetmask = 1
         velcode               ! 1: u comp, 2: v comp

      integer (kind=int_kind), dimension (nx_block*ny_block), & 
         intent(in) :: &
         indxti   , & ! compressed index in i-direction
         indxtj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         dxt      , & ! width of T-cell through the middle (m)
         dyt      , & ! height of T-cell through the middle (m)
         dxhy     , & ! 0.5*(HTE - HTE)
         dyhx     , & ! 0.5*(HTN - HTN)
         cyp      , & ! 1.5*HTE - 0.5*HTE
         cxp      , & ! 1.5*HTN - 0.5*HTN
         cym      , & ! 0.5*HTE - 1.5*HTE
         cxm          ! 0.5*HTN - 1.5*HTN

      real (kind=dbl_kind), dimension(nx_block,ny_block,8), & 
         intent(out) :: &
         Dstr          ! stress combinations
         
      real (kind=dbl_kind), dimension(nx_block,ny_block,4), & 
         intent(in) :: &
         zetaD          ! 2*zeta      

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &
        divune, divunw, divuse, divusw            , & ! divergence
        tensionne, tensionnw, tensionse, tensionsw, & ! tension
        shearne, shearnw, shearse, shearsw        , & ! shearing
        Deltane, Deltanw, Deltase, Deltasw        , & ! Delt
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
        strp_tmp, strm_tmp, tmp

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu

      if (velcode .eq. 1) then
       
       uij   = c1
       ui1j  = c0
       uij1  = c0
       ui1j1 = c0
       
       vij   = c0
       vi1j  = c0
       vij1  = c0
       vi1j1 = c0
       
       Dstr(:,:,:) = c0
       
      elseif (velcode .eq. 2) then
       
       uij   = c0
       ui1j  = c0
       uij1  = c0
       ui1j1 = c0
       
       vij   = c1
       vi1j  = c0
       vij1  = c0
       vi1j1 = c0
       
      endif
      

      do ij = 1, icellt
         i = indxti(ij)
         j = indxtj(ij)

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

          if (velcode .eq. 1) then
         
      !-----------------------------------------------------------------
      ! for dF/dx (u momentum)
      !-----------------------------------------------------------------
         strp_tmp  = p25*dyt(i,j)*(p333*ssigpn  + p166*ssigps)
         strm_tmp  = p25*dyt(i,j)*(p333*ssigmn  + p166*ssigms)

         ! northeast (i,j)
         Dstr(i,j,1) = -strp_tmp - strm_tmp - str12ew &
              + dxhy(i,j)*(-csigpne + csigmne) + dyhx(i,j)*csig12ne

         ! northwest (i+1,j)
         Dstr(i,j,2) = strp_tmp + strm_tmp - str12we &
              + dxhy(i,j)*(-csigpnw + csigmnw) + dyhx(i,j)*csig12nw

         strp_tmp  = p25*dyt(i,j)*(p333*ssigps  + p166*ssigpn)
         strm_tmp  = p25*dyt(i,j)*(p333*ssigms  + p166*ssigmn)

         ! southeast (i,j+1)
         Dstr(i,j,3) = -strp_tmp - strm_tmp + str12ew &
              + dxhy(i,j)*(-csigpse + csigmse) + dyhx(i,j)*csig12se

         ! southwest (i+1,j+1)
         Dstr(i,j,4) = strp_tmp + strm_tmp + str12we &
              + dxhy(i,j)*(-csigpsw + csigmsw) + dyhx(i,j)*csig12sw

        elseif (velcode .eq. 2) then      
              
      !-----------------------------------------------------------------
      ! for dF/dy (v momentum)
      !-----------------------------------------------------------------
         strp_tmp  = p25*dxt(i,j)*(p333*ssigpe  + p166*ssigpw)
         strm_tmp  = p25*dxt(i,j)*(p333*ssigme  + p166*ssigmw)

         ! northeast (i,j)
         Dstr(i,j,5) = -strp_tmp + strm_tmp - str12ns &
              - dyhx(i,j)*(csigpne + csigmne) + dxhy(i,j)*csig12ne

         ! southeast (i,j+1)
         Dstr(i,j,6) = strp_tmp - strm_tmp - str12sn &
              - dyhx(i,j)*(csigpse + csigmse) + dxhy(i,j)*csig12se

         strp_tmp  = p25*dxt(i,j)*(p333*ssigpw  + p166*ssigpe)
         strm_tmp  = p25*dxt(i,j)*(p333*ssigmw  + p166*ssigme)

         ! northwest (i+1,j)
         Dstr(i,j,7) = -strp_tmp + strm_tmp + str12ns &
              - dyhx(i,j)*(csigpnw + csigmnw) + dxhy(i,j)*csig12nw

         ! southwest (i+1,j+1)
         Dstr(i,j,8) = strp_tmp - strm_tmp + str12sn &
              - dyhx(i,j)*(csigpsw + csigmsw) + dxhy(i,j)*csig12sw

         endif     
              
      enddo                     ! ij

      end subroutine formDiag_step1      
      
!=======================================================================

! Calc diagonal term related to rheology for precond 

      subroutine OLDformDiag_step1  (nx_block,   ny_block,   & 
                                  icellt,     & 
                                  indxti,     indxtj,     & 
                                  dxt,        dyt,        & 
                                  dxhy,       dyhx,       & 
                                  cxp,        cyp,        & 
                                  cxm,        cym,        & 
                                  zetaD,      Dstr )

      integer (kind=int_kind), intent(in) :: & 
         nx_block, ny_block, & ! block dimensions
         icellt                ! no. of cells where icetmask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), & 
         intent(in) :: &
         indxti   , & ! compressed index in i-direction
         indxtj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         dxt      , & ! width of T-cell through the middle (m)
         dyt      , & ! height of T-cell through the middle (m)
         dxhy     , & ! 0.5*(HTE - HTE)
         dyhx     , & ! 0.5*(HTN - HTN)
         cyp      , & ! 1.5*HTE - 0.5*HTE
         cxp      , & ! 1.5*HTN - 0.5*HTN
         cym      , & ! 0.5*HTE - 1.5*HTE
         cxm          ! 0.5*HTN - 1.5*HTN

      real (kind=dbl_kind), dimension(nx_block,ny_block,8), & 
         intent(out) :: &
         Dstr          ! stress combinations
         
      real (kind=dbl_kind), dimension(nx_block,ny_block,4), & 
         intent(in) :: &
         zetaD          ! 2*zeta      

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &
        divuneDu, divunwDu, divuseDu, divuswDu            , & ! divergence
        divuneDv, divunwDv, divuseDv, divuswDv            , & ! divergence
        tensionneDu, tensionnwDu, tensionseDu, tensionswDu, & ! tension
        tensionneDv, tensionnwDv, tensionseDv, tensionswDv, & ! tension
        shearneDu, shearnwDu, shearseDu, shearswDu        , & ! shearing
        shearneDv, shearnwDv, shearseDv, shearswDv        , & ! shearing
        stressp_1u, stressp_2u, stressp_3u, stressp_4u    , & 
        stressp_1v, stressp_2v, stressp_3v, stressp_4v    , & 
        stressm_1u, stressm_2u, stressm_3u, stressm_4u    , & 
        stressm_1v, stressm_2v, stressm_3v, stressm_4v    , & 
        stress12_1u, stress12_2u, stress12_3u, stress12_4u, & 
        stress12_1v, stress12_2v, stress12_3v, stress12_4v, & 
        ssigpnu, ssigpsu, ssigpeu, ssigpwu                , &
        ssigpnv, ssigpsv, ssigpev, ssigpwv                , &
        ssigmnu, ssigmsu, ssigmeu, ssigmwu                , &
        ssigmnv, ssigmsv, ssigmev, ssigmwv                , &
        ssig12nu, ssig12su, ssig12eu, ssig12wu            , &
        ssig12nv, ssig12sv, ssig12ev, ssig12wv            , &
        ssigp1u, ssigp2u, ssigm1u, ssigm2u, ssig121u, ssig122u, &
        ssigp1v, ssigp2v, ssigm1v, ssigm2v, ssig121v, ssig122v, &
        csigpneu, csigpnwu, csigpseu, csigpswu            , &
        csigpnev, csigpnwv, csigpsev, csigpswv            , &
        csigmneu, csigmnwu, csigmseu, csigmswu            , &
        csigmnev, csigmnwv, csigmsev, csigmswv            , &
        csig12neu, csig12nwu, csig12seu, csig12swu        , &
        csig12nev, csig12nwv, csig12sev, csig12swv        , &
        str12ewu, str12weu, str12nsu, str12snu            , &
        str12ewv, str12wev, str12nsv, str12snv            , &
        strp_tmpu, strm_tmpu, strp_tmpv, strm_tmpv        , &
        str1, str2, str3, str4, str5, str6, str7, str8  
        
      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      Dstr(:,:,:) = c0

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icellt
         i = indxti(ij)
         j = indxtj(ij)

      !-----------------------------------------------------------------
      ! JFL watchout if currently on LHS or RHS
      !-----------------------------------------------------------------
         ! divergence  =  e_11 + e_22
         divuneDu  = cyp(i,j)
         divuneDv  = cxp(i,j)
         divunwDu  = dyt(i,j)
         divuseDv  = dxt(i,j)

         ! tension strain rate  =  e_11 - e_22
         tensionneDu = -cym(i,j)
         tensionneDv = cxm(i,j)
         tensionnwDu = dyt(i,j)
         tensionseDv = - dxt(i,j)

         ! shearing strain rate  =  e_12
         shearneDu = -cxm(i,j)
         shearneDv = -cym(i,j)
         shearnwDv = dyt(i,j)
         shearseDu = dxt(i,j)
         
      !-----------------------------------------------------------------
      ! the stresses                            ! kg/s^2
      ! (1) northeast, (2) northwest, (3) southwest, (4) southeast
      !-----------------------------------------------------------------

! IMProve: delete stress coeff not needed instead of setting them to 0.
!          no need for divuneDu...just plug them directly in eqs below. 
      
         stressp_1u = zetaD(i,j,1)*divuneDu*(c1+Ktens)
         stressp_1v = zetaD(i,j,1)*divuneDv*(c1+Ktens)
         stressp_2u = zetaD(i,j,2)*divunwDu*(c1+Ktens)
         stressp_2v = c0
         stressp_3u = c0
         stressp_3v = c0
         stressp_4u = c0
         stressp_4v = zetaD(i,j,4)*divuseDv*(c1+Ktens)
         
         stressm_1u = zetaD(i,j,1)*tensionneDu*(c1+Ktens)*ecci
         stressm_1v = zetaD(i,j,1)*tensionneDv*(c1+Ktens)*ecci
         stressm_2u = zetaD(i,j,2)*tensionnwDu*(c1+Ktens)*ecci
         stressm_2v = c0
         stressm_3u = c0
         stressm_3v = c0
         stressm_4u = c0
         stressm_4v = zetaD(i,j,4)*tensionseDv*(c1+Ktens)*ecci
         
         stress12_1u = zetaD(i,j,1)*shearneDu*p5*(c1+Ktens)*ecci
         stress12_1v = zetaD(i,j,1)*shearneDv*p5*(c1+Ktens)*ecci
         stress12_2u = c0
         stress12_2v = zetaD(i,j,2)*shearnwDv*p5*(c1+Ktens)*ecci
         stress12_3u = c0
         stress12_3v = c0
         stress12_4u = zetaD(i,j,4)*shearseDu*p5*(c1+Ktens)*ecci
         stress12_4v = c0

      !-----------------------------------------------------------------
      ! combinations of the stresses for the momentum equation ! kg/s^2
      !-----------------------------------------------------------------

         ssigpnu  = stressp_1u + stressp_2u
         ssigpnv  = stressp_1v + stressp_2v
         ssigpsu  = stressp_3u + stressp_4u
         ssigpsv  = stressp_3v + stressp_4v
         ssigpeu  = stressp_1u + stressp_4v
         ssigpev  = stressp_1v + stressp_4v
         ssigpwu  = stressp_2u + stressp_3u
         ssigpwv  = stressp_2v + stressp_3v
         ssigp1u  =(stressp_1u + stressp_3u)*p055
         ssigp1v  =(stressp_1v + stressp_3v)*p055
         ssigp2u  =(stressp_2u + stressp_4u)*p055
         ssigp2v  =(stressp_2v + stressp_4v)*p055

         ssigmnu  = stressm_1u + stressm_2u
         ssigmnv  = stressm_1v + stressm_2v
         ssigmsu  = stressm_3u + stressm_4u
         ssigmsv  = stressm_3v + stressm_4v
         ssigmeu  = stressm_1u + stressm_4u
         ssigmev  = stressm_1v + stressm_4v
         ssigmwu  = stressm_2u + stressm_3u
         ssigmwv  = stressm_2v + stressm_3v
         ssigm1u  =(stressm_1u + stressm_3u)*p055
         ssigm1v  =(stressm_1v + stressm_3v)*p055
         ssigm2u  =(stressm_2u + stressm_4u)*p055
         ssigm2v  =(stressm_2v + stressm_4v)*p055

         ssig12nu = stress12_1u + stress12_2u
         ssig12nv = stress12_1v + stress12_2v
         ssig12su = stress12_3u + stress12_4u
         ssig12sv = stress12_3v + stress12_4v
         ssig12eu = stress12_1u + stress12_4u
         ssig12ev = stress12_1v + stress12_4v
         ssig12wu = stress12_2u + stress12_3u
         ssig12wv = stress12_2v + stress12_3v
         ssig121u =(stress12_1u + stress12_3u)*p111
         ssig121v =(stress12_1v + stress12_3v)*p111
         ssig122u =(stress12_2u + stress12_4u)*p111
         ssig122v =(stress12_2v + stress12_4v)*p111

         csigpneu = p111*stressp_1u + ssigp2u + p027*stressp_3u
         csigpnev = p111*stressp_1v + ssigp2v + p027*stressp_3v
         csigpnwu = p111*stressp_2u + ssigp1u + p027*stressp_4u
         csigpnwv = p111*stressp_2v + ssigp1v + p027*stressp_4v
         csigpswu = p111*stressp_3u + ssigp2u + p027*stressp_1u
         csigpswv = p111*stressp_3v + ssigp2v + p027*stressp_1v
         csigpseu = p111*stressp_4u + ssigp1u + p027*stressp_2u
         csigpsev = p111*stressp_4v + ssigp1v + p027*stressp_2v
         
         csigmneu = p111*stressm_1u + ssigm2u + p027*stressm_3u
         csigmnev = p111*stressm_1v + ssigm2v + p027*stressm_3v
         csigmnwu = p111*stressm_2u + ssigm1u + p027*stressm_4u
         csigmnwv = p111*stressm_2v + ssigm1v + p027*stressm_4v
         csigmswu = p111*stressm_3u + ssigm2u + p027*stressm_1u
         csigmswv = p111*stressm_3v + ssigm2v + p027*stressm_1v
         csigmseu = p111*stressm_4u + ssigm1u + p027*stressm_2u
         csigmsev = p111*stressm_4v + ssigm1v + p027*stressm_2v
         
         csig12neu = p222*stress12_1u + ssig122u &
                  + p055*stress12_3u
         csig12nev = p222*stress12_1v + ssig122v &
                  + p055*stress12_3v             
         csig12nwu = p222*stress12_2u + ssig121u &
                  + p055*stress12_4u
         csig12nwv = p222*stress12_2v + ssig121v &
                  + p055*stress12_4v
         csig12swu = p222*stress12_3u + ssig122u &
                  + p055*stress12_1u
         csig12swv = p222*stress12_3v + ssig122v &
                  + p055*stress12_1v
         csig12seu = p222*stress12_4u + ssig121u &
                  + p055*stress12_2u
         csig12sev = p222*stress12_4v + ssig121v &
                  + p055*stress12_2v         

         str12ewu = p5*dxt(i,j)*(p333*ssig12eu + p166*ssig12wu)
         str12ewv = p5*dxt(i,j)*(p333*ssig12ev + p166*ssig12wv)
         str12weu = p5*dxt(i,j)*(p333*ssig12wu + p166*ssig12eu)
         str12wev = p5*dxt(i,j)*(p333*ssig12wv + p166*ssig12ev)
         str12nsu = p5*dyt(i,j)*(p333*ssig12nu + p166*ssig12su)
         str12nsv = p5*dyt(i,j)*(p333*ssig12nv + p166*ssig12sv)
         str12snu = p5*dyt(i,j)*(p333*ssig12su + p166*ssig12nu)
         str12snv = p5*dyt(i,j)*(p333*ssig12sv + p166*ssig12nv)

      !-----------------------------------------------------------------
      ! for dF/dx (u momentum)
      !-----------------------------------------------------------------
         strp_tmpu  = p25*dyt(i,j)*(p333*ssigpnu  + p166*ssigpsu)
         strm_tmpu  = p25*dyt(i,j)*(p333*ssigmnu  + p166*ssigmsu)

         ! northeast (i,j)
         Dstr(i,j,1) = -strp_tmpu - strm_tmpu - str12ewu &
              + dxhy(i,j)*(-csigpneu + csigmneu) + dyhx(i,j)*csig12neu

         ! northwest (i+1,j)
         Dstr(i,j,2) = strp_tmpu + strm_tmpu - str12weu &
              + dxhy(i,j)*(-csigpnwu + csigmnwu) + dyhx(i,j)*csig12nwu

         strp_tmpu  = p25*dyt(i,j)*(p333*ssigpsu  + p166*ssigpnu)
         strm_tmpu  = p25*dyt(i,j)*(p333*ssigmsu  + p166*ssigmnu)

         ! southeast (i,j+1)
         Dstr(i,j,3) = -strp_tmpu - strm_tmpu + str12ewu &
              + dxhy(i,j)*(-csigpseu + csigmseu) + dyhx(i,j)*csig12seu

         ! southwest (i+1,j+1)
         Dstr(i,j,4) = strp_tmpu + strm_tmpu + str12weu &
              + dxhy(i,j)*(-csigpswu + csigmswu) + dyhx(i,j)*csig12swu

      !-----------------------------------------------------------------
      ! for dF/dy (v momentum)
      !-----------------------------------------------------------------
         strp_tmpv  = p25*dxt(i,j)*(p333*ssigpev  + p166*ssigpwv)
         strm_tmpv  = p25*dxt(i,j)*(p333*ssigmev  + p166*ssigmwv)

         ! northeast (i,j)
         Dstr(i,j,5) = -strp_tmpv + strm_tmpv - str12nsv &
              - dyhx(i,j)*(csigpnev + csigmnev) + dxhy(i,j)*csig12nev

         ! southeast (i,j+1)
         Dstr(i,j,6) = strp_tmpv - strm_tmpv - str12snv &
              - dyhx(i,j)*(csigpsev + csigmsev) + dxhy(i,j)*csig12sev

         strp_tmpv  = p25*dxt(i,j)*(p333*ssigpwv  + p166*ssigpev)
         strm_tmpv  = p25*dxt(i,j)*(p333*ssigmwv  + p166*ssigmev)

         ! northwest (i+1,j)
         Dstr(i,j,7) = -strp_tmpv + strm_tmpv + str12nsv &
              - dyhx(i,j)*(csigpnwv + csigmnwv) + dxhy(i,j)*csig12nwv

         ! southwest (i+1,j+1)
         Dstr(i,j,8) = strp_tmpv - strm_tmpv + str12snv &
              - dyhx(i,j)*(csigpswv + csigmswv) + dxhy(i,j)*csig12swv

      enddo                     ! ij

      end subroutine OLDformDiag_step1
      
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

      !-----------------------------------------------------------------
      ! integrate the momentum equation
      !-----------------------------------------------------------------

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message="subname", &
         file=__FILE__, line=__LINE__)

      do ij =1, icellu
         i = indxui(ij)
         j = indxuj(ij)

         ccaimp = umassdti(i,j) + vrel(i,j) * cosw + Cb(i,j) ! kg/m^2 s
               
         ! divergence of the internal stress tensor
         strintx = uarear(i,j)* &
             (Dstr(i,j,1) + Dstr(i+1,j,2) + Dstr(i,j+1,3) + Dstr(i+1,j+1,4))
         strinty = uarear(i,j)* &
             (Dstr(i,j,5) + Dstr(i,j+1,6) + Dstr(i+1,j,7) + Dstr(i+1,j+1,8))

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

      !-----------------------------------------------------------------
      ! form vector (converts from max_blocks arrays to single vector
      !-----------------------------------------------------------------

      wk2(:)=c0
      
      do i=1, ntot

	wk2(i) = wk1(i)/diagvec(i)
      
      enddo! i

      end subroutine precond_diag
      
!=======================================================================      

      subroutine calc_L2norm (nx_block,   ny_block, &
                              icellu,               &
                              indxui,     indxuj,   &
                              tpu,        tpv )

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellu                ! total count when iceumask is true

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxui  , & ! compressed index in i-direction
         indxuj      ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         tpu     , & ! x-component of vector
         tpv         ! y-component of vector         

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij
         
      real (kind=dbl_kind) :: L2norm

      !-----------------------------------------------------------------
      ! form vector
      !-----------------------------------------------------------------

     L2norm = c0
      
      do ij =1, icellu
         i = indxui(ij)
         j = indxuj(ij)
         
         L2norm = L2norm + tpu(i,j)**2
         L2norm = L2norm + tpv(i,j)**2
         
      enddo                     ! ij

      L2norm = sqrt(L2norm)
      
      end subroutine calc_L2norm
      
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

      end module ice_dyn_vp

!=======================================================================
