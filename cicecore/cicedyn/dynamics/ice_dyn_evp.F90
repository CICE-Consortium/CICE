!=======================================================================
!
! Elastic-viscous-plastic sea ice dynamics model
! Computes ice velocity and deformation
!
! See:
!
! Hunke, E. C., and J. K. Dukowicz (1997). An elastic-viscous-plastic model
! for sea ice dynamics. J. Phys. Oceanogr., 27, 1849-1867.
!
! Hunke, E. C. (2001).  Viscous-Plastic Sea Ice Dynamics with the EVP Model:
! Linearization Issues. J. Comput. Phys., 170, 18-38.
!
! Hunke, E. C., and J. K. Dukowicz (2002).  The Elastic-Viscous-Plastic
! Sea Ice Dynamics Model in General Orthogonal Curvilinear Coordinates
! on a Sphere - Incorporation of Metric Terms. Mon. Weather Rev.,
! 130, 1848-1865.
!
! Hunke, E. C., and J. K. Dukowicz (2003).  The sea ice momentum
! equation in the free drift regime.  Los Alamos Tech. Rep. LA-UR-03-2219.
!
! Hibler, W. D. (1979). A dynamic thermodynamic sea ice model. J. Phys.
! Oceanogr., 9, 817-846.
!
! Bouillon, S., T. Fichefet, V. Legat and G. Madec (2013).  The
! elastic-viscous-plastic method revisited.  Ocean Model., 71, 2-12.
!
! author: Elizabeth C. Hunke, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb (LANL)
! 2004: Block structure added by William Lipscomb
! 2005: Removed boundary calls for stress arrays (WHL)
! 2006: Streamlined for efficiency by Elizabeth Hunke
!       Converted to free source form (F90)

      module ice_dyn_evp

      use ice_kinds_mod
      use ice_communicate, only: my_task, master_task
      use ice_constants, only: field_loc_center, field_loc_NEcorner, &
          field_loc_Nface, field_loc_Eface, &
          field_type_scalar, field_type_vector
      use ice_constants, only: c0, p027, p055, p111, p166, &
          p222, p25, p333, p5, c1
      use ice_dyn_shared, only: stepu, stepuv_CD, stepu_C, stepv_C, &
          dyn_prep1, dyn_prep2, dyn_finish, &
          ndte, yield_curve, ecci, denom1, arlx1i, fcor_blk, fcorE_blk, fcorN_blk, &
          uvel_init, vvel_init, uvelE_init, vvelE_init, uvelN_init, vvelN_init, &
          seabed_stress_factor_LKD, seabed_stress_factor_prob, seabed_stress_method, &
          seabed_stress, Ktens, revp
      use ice_fileunits, only: nu_diag
      use ice_exit, only: abort_ice
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_ice_strength, icepack_query_parameters

      implicit none
      private
! all c or cd
      real (kind=dbl_kind), allocatable :: &
         uocnN    (:,:,:) , & ! i ocean current (m/s)
         vocnN    (:,:,:) , & ! j ocean current (m/s)
         ss_tltxN (:,:,:) , & ! sea surface slope, x-direction (m/m)
         ss_tltyN (:,:,:) , & ! sea surface slope, y-direction (m/m)
         cdn_ocnN (:,:,:) , & ! ocn drag coefficient
         waterxN  (:,:,:) , & ! for ocean stress calculation, x (m/s)
         wateryN  (:,:,:) , & ! for ocean stress calculation, y (m/s)
         forcexN  (:,:,:) , & ! work array: combined atm stress and ocn tilt, x
         forceyN  (:,:,:) , & ! work array: combined atm stress and ocn tilt, y
         aiN      (:,:,:) , & ! ice fraction on N-grid
         nmass    (:,:,:) , & ! total mass of ice and snow (N grid)
         nmassdti (:,:,:)     ! mass of N-cell/dte (kg/m^2 s)
! all c or d
      real (kind=dbl_kind), allocatable :: &
         uocnE    (:,:,:) , & ! i ocean current (m/s)
         vocnE    (:,:,:) , & ! j ocean current (m/s)
         ss_tltxE (:,:,:) , & ! sea surface slope, x-direction (m/m)
         ss_tltyE (:,:,:) , & ! sea surface slope, y-direction (m/m)
         cdn_ocnE (:,:,:) , & ! ocn drag coefficient
         waterxE  (:,:,:) , & ! for ocean stress calculation, x (m/s)
         wateryE  (:,:,:) , & ! for ocean stress calculation, y (m/s)
         forcexE  (:,:,:) , & ! work array: combined atm stress and ocn tilt, x
         forceyE  (:,:,:) , & ! work array: combined atm stress and ocn tilt, y
         aiE      (:,:,:) , & ! ice fraction on E-grid
         emass    (:,:,:) , & ! total mass of ice and snow (E grid)
         emassdti (:,:,:)     ! mass of E-cell/dte (kg/m^2 s)

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         ratiodxN    , & ! - dxN(i+1,j)   / dxN(i,j)
         ratiodyE    , & ! - dyE(i  ,j+1) / dyE(i,j)
         ratiodxNr   , & !   1 / ratiodxN
         ratiodyEr       !   1 / ratiodyE

      real (kind=dbl_kind), allocatable :: &
         strengthU(:,:,:) , & ! strength averaged to U points
         divergU  (:,:,:) , & ! div array on U points, differentiate from divu
         tensionU (:,:,:) , & ! tension array on U points
         shearU   (:,:,:) , & ! shear array on U points
         deltaU   (:,:,:) , & ! delta array on U points
         zetax2T  (:,:,:) , & ! zetax2 = 2*zeta (bulk viscosity)
         zetax2U  (:,:,:) , & ! zetax2T averaged to U points
         etax2T   (:,:,:) , & ! etax2  = 2*eta  (shear viscosity)
         etax2U   (:,:,:)     ! etax2T averaged to U points

      real (kind=dbl_kind), allocatable :: &
         uocnU    (:,:,:) , & ! i ocean current (m/s)
         vocnU    (:,:,:) , & ! j ocean current (m/s)
         ss_tltxU (:,:,:) , & ! sea surface slope, x-direction (m/m)
         ss_tltyU (:,:,:) , & ! sea surface slope, y-direction (m/m)
         cdn_ocnU (:,:,:) , & ! ocn drag coefficient
         tmass    (:,:,:) , & ! total mass of ice and snow (kg/m^2)
         waterxU  (:,:,:) , & ! for ocean stress calculation, x (m/s)
         wateryU  (:,:,:) , & ! for ocean stress calculation, y (m/s)
         forcexU  (:,:,:) , & ! work array: combined atm stress and ocn tilt, x
         forceyU  (:,:,:) , & ! work array: combined atm stress and ocn tilt, y
         umass    (:,:,:) , & ! total mass of ice and snow (u grid)
         umassdti (:,:,:)     ! mass of U-cell/dte (kg/m^2 s)

      public :: evp, init_evp

!=======================================================================

      contains

!=======================================================================
! Elastic-viscous-plastic dynamics driver
!
      subroutine init_evp
      use ice_blocks, only: get_block, nx_block, ny_block, nghost, block
      use ice_domain_size, only: max_blocks
      use ice_domain, only: nblocks, blocks_ice
      use ice_grid, only: grid_ice, dyT, dxT, uarear, tmask, G_HTE, G_HTN, dxN, dyE
      use ice_calendar, only: dt_dyn
      use ice_dyn_shared, only: init_dyn_shared, evp_algorithm
      use ice_dyn_evp1d, only: dyn_evp1d_init

!allocate c and cd grid var. Follow structucre of eap
      integer (int_kind) :: ierr

      character(len=*), parameter :: subname = '(init_evp)'

     type (block) :: &
         this_block   ! block information for current block

      integer (kind=int_kind) :: &
         i, j, iblk            , & ! block index
         ilo,ihi,jlo,jhi           ! beginning and end of physical domain


      call init_dyn_shared(dt_dyn)

      if (evp_algorithm == "shared_mem_1d" ) then
         call dyn_evp1d_init
      endif

      allocate( uocnU    (nx_block,ny_block,max_blocks), & ! i ocean current (m/s)
                vocnU    (nx_block,ny_block,max_blocks), & ! j ocean current (m/s)
                ss_tltxU (nx_block,ny_block,max_blocks), & ! sea surface slope, x-direction (m/m)
                ss_tltyU (nx_block,ny_block,max_blocks), & ! sea surface slope, y-direction (m/m)
                cdn_ocnU (nx_block,ny_block,max_blocks), & ! ocn drag coefficient
                tmass    (nx_block,ny_block,max_blocks), & ! total mass of ice and snow (kg/m^2)
                waterxU  (nx_block,ny_block,max_blocks), & ! for ocean stress calculation, x (m/s)
                wateryU  (nx_block,ny_block,max_blocks), & ! for ocean stress calculation, y (m/s)
                forcexU  (nx_block,ny_block,max_blocks), & ! work array: combined atm stress and ocn tilt, x
                forceyU  (nx_block,ny_block,max_blocks), & ! work array: combined atm stress and ocn tilt, y
                umass    (nx_block,ny_block,max_blocks), & ! total mass of ice and snow (u grid)
                umassdti (nx_block,ny_block,max_blocks), & ! mass of U-cell/dte (kg/m^2 s)
                stat=ierr)
      if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory B-Grid evp')


      if (grid_ice == 'CD' .or. grid_ice == 'C') then

         allocate( strengthU(nx_block,ny_block,max_blocks), &
                   divergU  (nx_block,ny_block,max_blocks), &
                   tensionU (nx_block,ny_block,max_blocks), &
                   shearU   (nx_block,ny_block,max_blocks), &
                   deltaU   (nx_block,ny_block,max_blocks), &
                   zetax2T  (nx_block,ny_block,max_blocks), &
                   zetax2U  (nx_block,ny_block,max_blocks), &
                   etax2T   (nx_block,ny_block,max_blocks), &
                   etax2U   (nx_block,ny_block,max_blocks), &
                   stat=ierr)
         if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory U evp')

         allocate( uocnN    (nx_block,ny_block,max_blocks), &
                   vocnN    (nx_block,ny_block,max_blocks), &
                   ss_tltxN (nx_block,ny_block,max_blocks), &
                   ss_tltyN (nx_block,ny_block,max_blocks), &
                   cdn_ocnN (nx_block,ny_block,max_blocks), &
                   waterxN  (nx_block,ny_block,max_blocks), &
                   wateryN  (nx_block,ny_block,max_blocks), &
                   forcexN  (nx_block,ny_block,max_blocks), &
                   forceyN  (nx_block,ny_block,max_blocks), &
                   aiN      (nx_block,ny_block,max_blocks), &
                   nmass    (nx_block,ny_block,max_blocks), &
                   nmassdti (nx_block,ny_block,max_blocks), &
                   stat=ierr)
         if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory N evp')

         allocate( uocnE    (nx_block,ny_block,max_blocks), &
                   vocnE    (nx_block,ny_block,max_blocks), &
                   ss_tltxE (nx_block,ny_block,max_blocks), &
                   ss_tltyE (nx_block,ny_block,max_blocks), &
                   cdn_ocnE (nx_block,ny_block,max_blocks), &
                   waterxE  (nx_block,ny_block,max_blocks), &
                   wateryE  (nx_block,ny_block,max_blocks), &
                   forcexE  (nx_block,ny_block,max_blocks), &
                   forceyE  (nx_block,ny_block,max_blocks), &
                   aiE      (nx_block,ny_block,max_blocks), &
                   emass    (nx_block,ny_block,max_blocks), &
                   emassdti (nx_block,ny_block,max_blocks), &
                   stat=ierr)
         if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory E evp')

         allocate( ratiodxN (nx_block,ny_block,max_blocks), &
                   ratiodyE (nx_block,ny_block,max_blocks), &
                   ratiodxNr(nx_block,ny_block,max_blocks), &
                   ratiodyEr(nx_block,ny_block,max_blocks), &
                   stat=ierr)
         if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory ratio')

         !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
         do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            do j = jlo, jhi
            do i = ilo, ihi
               ratiodxN (i,j,iblk) = - dxN(i+1,j  ,iblk) / dxN(i,j,iblk)
               ratiodyE (i,j,iblk) = - dyE(i  ,j+1,iblk) / dyE(i,j,iblk)
               ratiodxNr(i,j,iblk) =   c1 / ratiodxN(i,j,iblk)
               ratiodyEr(i,j,iblk) =   c1 / ratiodyE(i,j,iblk)
            enddo
            enddo
         enddo                     ! iblk
         !$OMP END PARALLEL DO

      endif

      end subroutine init_evp

!=======================================================================
#ifdef CICE_IN_NEMO
! Wind stress is set during this routine from the values supplied
! via NEMO (unless calc_strair is true).  These values are supplied
! rotated on u grid and multiplied by aice.  strairxT = 0 in this
! case so operations in dyn_prep1 are pointless but carried out to
! minimise code changes.
#endif
!
! author: Elizabeth C. Hunke, LANL

      subroutine evp (dt)

      use ice_arrays_column, only: Cdn_ocn
      use ice_boundary, only: ice_halo, ice_HaloMask, ice_HaloUpdate, &
          ice_HaloDestroy, ice_HaloUpdate_stress
      use ice_blocks, only: block, get_block, nx_block, ny_block, nghost
      use ice_domain, only: nblocks, blocks_ice, halo_info, maskhalo_dyn, &
          ns_boundary_type
      use ice_domain_size, only: max_blocks, ncat
      use ice_flux, only: rdg_conv, rdg_shear, strairxT, strairyT, &
          strairxU, strairyU, uocn, vocn, ss_tltx, ss_tlty, fmU, &
          strtltxU, strtltyU, strocnxU, strocnyU, strintxU, strintyU, taubxU, taubyU, &
          strax, stray, &
          TbU, hwater, &
          strairxN, strairyN, fmN, &
          strtltxN, strtltyN, strocnxN, strocnyN, strintxN, strintyN, taubxN, taubyN, &
          TbN, &
          strairxE, strairyE, fmE, &
          strtltxE, strtltyE, strocnxE, strocnyE, strintxE, strintyE, taubxE, taubyE, &
          TbE, &
          stressp_1, stressp_2, stressp_3, stressp_4, &
          stressm_1, stressm_2, stressm_3, stressm_4, &
          stress12_1, stress12_2, stress12_3, stress12_4, &
          stresspT, stressmT, stress12T, &
          stresspU, stressmU, stress12U
      use ice_grid, only: tmask, umask, umaskCD, nmask, emask, uvm, epm, npm, &
          dxE, dxN, dxT, dxU, dyE, dyN, dyT, dyU, &
          tarear, uarear, earear, narear, grid_average_X2Y, uarea, &
          grid_ice, grid_atm_dynu, grid_atm_dynv, grid_ocn_dynu, grid_ocn_dynv
      use ice_state, only: aice, aiU, vice, vsno, uvel, vvel, uvelN, vvelN, &
          uvelE, vvelE, divu, shear, vort, &
          aice_init, aice0, aicen, vicen, strength
      use ice_timers, only: timer_dynamics, timer_bound, &
          ice_timer_start, ice_timer_stop, timer_evp
      use ice_dyn_shared, only: evp_algorithm, stack_fields, unstack_fields, &
          DminTarea, visc_method, deformations, deformationsC_T, deformationsCD_T, &
          strain_rates_U, dxhy, dyhx, cxp, cyp, cxm, cym, &
          iceTmask, iceUmask, iceEmask, iceNmask, &
          dyn_haloUpdate, fld2, fld3, fld4
      use ice_dyn_evp1d, only: dyn_evp1d_run

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) :: &
         ksub           , & ! subcycle step
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, ij           ! local indices

      integer (kind=int_kind), dimension(max_blocks) :: &
         icellT   , & ! no. of cells where iceTmask = .true.
         icellN   , & ! no. of cells where iceNmask = .true.
         icellE   , & ! no. of cells where iceEmask = .true.
         icellU       ! no. of cells where iceUmask = .true.

      integer (kind=int_kind), dimension (nx_block*ny_block, max_blocks) :: &
         indxTi   , & ! compressed index in i-direction
         indxTj   , & ! compressed index in j-direction
         indxEi   , & ! compressed index in i-direction
         indxEj   , & ! compressed index in j-direction
         indxNi   , & ! compressed index in i-direction
         indxNj   , & ! compressed index in j-direction
         indxUi   , & ! compressed index in i-direction
         indxUj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension(nx_block,ny_block,8):: &
         strtmp       ! stress combinations for momentum equation

      logical (kind=log_kind) :: &
         calc_strair  ! calculate air/ice stress

      integer (kind=int_kind), dimension (nx_block,ny_block,max_blocks) :: &
         halomask     ! generic halo mask

      type (ice_halo) :: &
         halo_info_mask !  ghost cell update info for masked halo

      type (block) :: &
         this_block   ! block information for current block

      character(len=*), parameter :: subname = '(evp)'

      call ice_timer_start(timer_dynamics) ! dynamics

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      if (grid_ice == 'CD' .or. grid_ice == 'C') then

         strengthU(:,:,:) = c0
         divergU  (:,:,:) = c0
         tensionU (:,:,:) = c0
         shearU   (:,:,:) = c0
         deltaU   (:,:,:) = c0
         zetax2T  (:,:,:) = c0
         zetax2U  (:,:,:) = c0
         etax2T   (:,:,:) = c0
         etax2U   (:,:,:) = c0

      endif

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

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block) SCHEDULE(runtime)
      do iblk = 1, nblocks

         do j = 1, ny_block
         do i = 1, nx_block
            rdg_conv (i,j,iblk) = c0
            rdg_shear(i,j,iblk) = c0
            divu (i,j,iblk) = c0
            shear(i,j,iblk) = c0
            vort (i,j,iblk) = c0
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
                         tmass   (:,:,iblk), iceTmask(:,:,iblk))

      enddo                     ! iblk
      !$OMP END PARALLEL DO

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (iceTmask,          halo_info, &
                           field_loc_center,  field_type_scalar)
      call ice_timer_stop(timer_bound)

      !-----------------------------------------------------------------
      ! convert fields from T to U grid
      !-----------------------------------------------------------------

      call stack_fields(tmass, aice_init, cdn_ocn, fld3)
      call ice_HaloUpdate (fld3,             halo_info, &
                           field_loc_center, field_type_scalar)
      call stack_fields(uocn, vocn, ss_tltx, ss_tlty, fld4)
      call ice_HaloUpdate (fld4,             halo_info, &
                           field_loc_center, field_type_vector)
      call unstack_fields(fld3, tmass, aice_init, cdn_ocn)
      call unstack_fields(fld4, uocn, vocn, ss_tltx, ss_tlty)

      call grid_average_X2Y('S', tmass    , 'T'          , umass   , 'U')
      call grid_average_X2Y('S', aice_init, 'T'          , aiU     , 'U')
      call grid_average_X2Y('S', cdn_ocn  , 'T'          , cdn_ocnU, 'U')
      call grid_average_X2Y('S', uocn     , grid_ocn_dynu, uocnU   , 'U')
      call grid_average_X2Y('S', vocn     , grid_ocn_dynv, vocnU   , 'U')
      call grid_average_X2Y('S', ss_tltx  , grid_ocn_dynu, ss_tltxU, 'U')
      call grid_average_X2Y('S', ss_tlty  , grid_ocn_dynv, ss_tltyU, 'U')

      if (grid_ice == 'CD' .or. grid_ice == 'C') then
         call grid_average_X2Y('S', tmass    , 'T'          , emass   , 'E')
         call grid_average_X2Y('S', aice_init, 'T'          , aie     , 'E')
         call grid_average_X2Y('S', cdn_ocn  , 'T'          , cdn_ocnE, 'E')
         call grid_average_X2Y('S', uocn     , grid_ocn_dynu, uocnE   , 'E')
         call grid_average_X2Y('S', vocn     , grid_ocn_dynv, vocnE   , 'E')
         call grid_average_X2Y('S', ss_tltx  , grid_ocn_dynu, ss_tltxE, 'E')
         call grid_average_X2Y('S', ss_tlty  , grid_ocn_dynv, ss_tltyE, 'E')
         call grid_average_X2Y('S', tmass    , 'T'          , nmass   , 'N')
         call grid_average_X2Y('S', aice_init, 'T'          , ain     , 'N')
         call grid_average_X2Y('S', cdn_ocn  , 'T'          , cdn_ocnN, 'N')
         call grid_average_X2Y('S', uocn     , grid_ocn_dynu, uocnN   , 'N')
         call grid_average_X2Y('S', vocn     , grid_ocn_dynv, vocnN   , 'N')
         call grid_average_X2Y('S', ss_tltx  , grid_ocn_dynu, ss_tltxN, 'N')
         call grid_average_X2Y('S', ss_tlty  , grid_ocn_dynv, ss_tltyN, 'N')
      endif
      !----------------------------------------------------------------
      ! Set wind stress to values supplied via NEMO or other forcing
      ! Map T to U, N, E as needed
      ! This wind stress is rotated on u grid and multiplied by aice in coupler
      !----------------------------------------------------------------
      call icepack_query_parameters(calc_strair_out=calc_strair)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (.not. calc_strair) then
         call grid_average_X2Y('F', strax, grid_atm_dynu, strairxU, 'U')
         call grid_average_X2Y('F', stray, grid_atm_dynv, strairyU, 'U')
      else
         call ice_HaloUpdate (strairxT,         halo_info, &
                              field_loc_center, field_type_vector)
         call ice_HaloUpdate (strairyT,         halo_info, &
                              field_loc_center, field_type_vector)
         call grid_average_X2Y('F', strairxT, 'T', strairxU, 'U')
         call grid_average_X2Y('F', strairyT, 'T', strairyU, 'U')
      endif

      if (grid_ice == 'CD' .or. grid_ice == 'C') then
         if (.not. calc_strair) then
            call grid_average_X2Y('F', strax   , grid_atm_dynu, strairxN, 'N')
            call grid_average_X2Y('F', stray   , grid_atm_dynv, strairyN, 'N')
            call grid_average_X2Y('F', strax   , grid_atm_dynu, strairxE, 'E')
            call grid_average_X2Y('F', stray   , grid_atm_dynv, strairyE, 'E')
         else
            call grid_average_X2Y('F', strairxT, 'T'          , strairxN, 'N')
            call grid_average_X2Y('F', strairyT, 'T'          , strairyN, 'N')
            call grid_average_X2Y('F', strairxT, 'T'          , strairxE, 'E')
            call grid_average_X2Y('F', strairyT, 'T'          , strairyE, 'E')
         endif
      endif

      if (trim(grid_ice) == 'B') then
         !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block,ij,i,j) SCHEDULE(runtime)
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
                            icellT        (iblk), icellU        (iblk), &
                            indxTi      (:,iblk), indxTj      (:,iblk), &
                            indxUi      (:,iblk), indxUj      (:,iblk), &
                            aiU       (:,:,iblk), umass     (:,:,iblk), &
                            umassdti  (:,:,iblk), fcor_blk  (:,:,iblk), &
                            umask     (:,:,iblk),                       &
                            uocnU     (:,:,iblk), vocnU     (:,:,iblk), &
                            strairxU  (:,:,iblk), strairyU  (:,:,iblk), &
                            ss_tltxU  (:,:,iblk), ss_tltyU  (:,:,iblk), &
                            iceTmask  (:,:,iblk), iceUmask  (:,:,iblk), &
                            fmU       (:,:,iblk), dt                  , &
                            strtltxU  (:,:,iblk), strtltyU  (:,:,iblk), &
                            strocnxU  (:,:,iblk), strocnyU  (:,:,iblk), &
                            strintxU  (:,:,iblk), strintyU  (:,:,iblk), &
                            taubxU    (:,:,iblk), taubyU    (:,:,iblk), &
                            waterxU   (:,:,iblk), wateryU   (:,:,iblk), &
                            forcexU   (:,:,iblk), forceyU   (:,:,iblk), &
                            stressp_1 (:,:,iblk), stressp_2 (:,:,iblk), &
                            stressp_3 (:,:,iblk), stressp_4 (:,:,iblk), &
                            stressm_1 (:,:,iblk), stressm_2 (:,:,iblk), &
                            stressm_3 (:,:,iblk), stressm_4 (:,:,iblk), &
                            stress12_1(:,:,iblk), stress12_2(:,:,iblk), &
                            stress12_3(:,:,iblk), stress12_4(:,:,iblk), &
                            uvel_init (:,:,iblk), vvel_init (:,:,iblk), &
                            uvel      (:,:,iblk), vvel      (:,:,iblk), &
                            TbU       (:,:,iblk))

            !-----------------------------------------------------------------
            ! ice strength
            !-----------------------------------------------------------------

            strength(:,:,iblk) = c0  ! initialize
            do ij = 1, icellT(iblk)
               i = indxTi(ij, iblk)
               j = indxTj(ij, iblk)
               call icepack_ice_strength(aice     = aice    (i,j,  iblk), &
                                         vice     = vice    (i,j,  iblk), &
                                         aice0    = aice0   (i,j,  iblk), &
                                         aicen    = aicen   (i,j,:,iblk), &
                                         vicen    = vicen   (i,j,:,iblk), &
                                         strength = strength(i,j,  iblk))
            enddo  ! ij

         enddo  ! iblk
         !$OMP END PARALLEL DO
      elseif (trim(grid_ice) == 'CD' .or. grid_ice == 'C') then
         !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block,ij,i,j) SCHEDULE(runtime)
         do iblk = 1, nblocks

            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            call dyn_prep2 (nx_block,             ny_block,             &
                            ilo, ihi,             jlo, jhi,             &
                            icellT        (iblk), icellU        (iblk), &
                            indxTi      (:,iblk), indxTj      (:,iblk), &
                            indxUi      (:,iblk), indxUj      (:,iblk), &
                            aiU       (:,:,iblk), umass     (:,:,iblk), &
                            umassdti  (:,:,iblk), fcor_blk  (:,:,iblk), &
                            umaskCD   (:,:,iblk),                       &
                            uocnU     (:,:,iblk), vocnU     (:,:,iblk), &
                            strairxU  (:,:,iblk), strairyU  (:,:,iblk), &
                            ss_tltxU  (:,:,iblk), ss_tltyU  (:,:,iblk), &
                            iceTmask  (:,:,iblk), iceUmask  (:,:,iblk), &
                            fmU       (:,:,iblk), dt,                   &
                            strtltxU  (:,:,iblk), strtltyU  (:,:,iblk), &
                            strocnxU  (:,:,iblk), strocnyU  (:,:,iblk), &
                            strintxU  (:,:,iblk), strintyU  (:,:,iblk), &
                            taubxU    (:,:,iblk), taubyU    (:,:,iblk), &
                            waterxU   (:,:,iblk), wateryU   (:,:,iblk), &
                            forcexU   (:,:,iblk), forceyU   (:,:,iblk), &
                            stressp_1 (:,:,iblk), stressp_2 (:,:,iblk), &
                            stressp_3 (:,:,iblk), stressp_4 (:,:,iblk), &
                            stressm_1 (:,:,iblk), stressm_2 (:,:,iblk), &
                            stressm_3 (:,:,iblk), stressm_4 (:,:,iblk), &
                            stress12_1(:,:,iblk), stress12_2(:,:,iblk), &
                            stress12_3(:,:,iblk), stress12_4(:,:,iblk), &
                            uvel_init (:,:,iblk), vvel_init (:,:,iblk), &
                            uvel      (:,:,iblk), vvel      (:,:,iblk), &
                            TbU       (:,:,iblk))

            !-----------------------------------------------------------------
            ! ice strength
            !-----------------------------------------------------------------

            strength(:,:,iblk) = c0  ! initialize
            do ij = 1, icellT(iblk)
               i = indxTi(ij, iblk)
               j = indxTj(ij, iblk)
               call icepack_ice_strength(aice     = aice    (i,j,  iblk), &
                                         vice     = vice    (i,j,  iblk), &
                                         aice0    = aice0   (i,j,  iblk), &
                                         aicen    = aicen   (i,j,:,iblk), &
                                         vicen    = vicen   (i,j,:,iblk), &
                                         strength = strength(i,j,  iblk) )
            enddo  ! ij


            !-----------------------------------------------------------------
            ! more preparation for dynamics on N grid
            !-----------------------------------------------------------------

            call dyn_prep2 (nx_block,             ny_block,             &
                            ilo, ihi,             jlo, jhi,             &
                            icellT        (iblk), icellN        (iblk), &
                            indxTi      (:,iblk), indxTj      (:,iblk), &
                            indxNi      (:,iblk), indxNj      (:,iblk), &
                            aiN       (:,:,iblk), nmass     (:,:,iblk), &
                            nmassdti  (:,:,iblk), fcorN_blk (:,:,iblk), &
                            nmask     (:,:,iblk),                       &
                            uocnN     (:,:,iblk), vocnN     (:,:,iblk), &
                            strairxN  (:,:,iblk), strairyN  (:,:,iblk), &
                            ss_tltxN  (:,:,iblk), ss_tltyN  (:,:,iblk), &
                            iceTmask  (:,:,iblk), iceNmask  (:,:,iblk), &
                            fmN       (:,:,iblk), dt                  , &
                            strtltxN  (:,:,iblk), strtltyN  (:,:,iblk), &
                            strocnxN  (:,:,iblk), strocnyN  (:,:,iblk), &
                            strintxN  (:,:,iblk), strintyN  (:,:,iblk), &
                            taubxN    (:,:,iblk), taubyN    (:,:,iblk), &
                            waterxN   (:,:,iblk), wateryN   (:,:,iblk), &
                            forcexN   (:,:,iblk), forceyN   (:,:,iblk), &
                            stressp_1 (:,:,iblk), stressp_2 (:,:,iblk), &
                            stressp_3 (:,:,iblk), stressp_4 (:,:,iblk), &
                            stressm_1 (:,:,iblk), stressm_2 (:,:,iblk), &
                            stressm_3 (:,:,iblk), stressm_4 (:,:,iblk), &
                            stress12_1(:,:,iblk), stress12_2(:,:,iblk), &
                            stress12_3(:,:,iblk), stress12_4(:,:,iblk), &
                            uvelN_init(:,:,iblk), vvelN_init(:,:,iblk), &
                            uvelN     (:,:,iblk), vvelN     (:,:,iblk), &
                            TbN       (:,:,iblk))

            !-----------------------------------------------------------------
            ! more preparation for dynamics on E grid
            !-----------------------------------------------------------------

            call dyn_prep2 (nx_block,             ny_block,             &
                            ilo, ihi,             jlo, jhi,             &
                            icellT        (iblk), icellE        (iblk), &
                            indxTi      (:,iblk), indxTj      (:,iblk), &
                            indxEi      (:,iblk), indxEj      (:,iblk), &
                            aiE       (:,:,iblk), emass     (:,:,iblk), &
                            emassdti  (:,:,iblk), fcorE_blk (:,:,iblk), &
                            emask     (:,:,iblk),                       &
                            uocnE     (:,:,iblk), vocnE     (:,:,iblk), &
                            strairxE  (:,:,iblk), strairyE  (:,:,iblk), &
                            ss_tltxE  (:,:,iblk), ss_tltyE  (:,:,iblk), &
                            iceTmask  (:,:,iblk), iceEmask  (:,:,iblk), &
                            fmE       (:,:,iblk), dt                  , &
                            strtltxE  (:,:,iblk), strtltyE  (:,:,iblk), &
                            strocnxE  (:,:,iblk), strocnyE  (:,:,iblk), &
                            strintxE  (:,:,iblk), strintyE  (:,:,iblk), &
                            taubxE    (:,:,iblk), taubyE    (:,:,iblk), &
                            waterxE   (:,:,iblk), wateryE   (:,:,iblk), &
                            forcexE   (:,:,iblk), forceyE   (:,:,iblk), &
                            stressp_1 (:,:,iblk), stressp_2 (:,:,iblk), &
                            stressp_3 (:,:,iblk), stressp_4 (:,:,iblk), &
                            stressm_1 (:,:,iblk), stressm_2 (:,:,iblk), &
                            stressm_3 (:,:,iblk), stressm_4 (:,:,iblk), &
                            stress12_1(:,:,iblk), stress12_2(:,:,iblk), &
                            stress12_3(:,:,iblk), stress12_4(:,:,iblk), &
                            uvelE_init(:,:,iblk), vvelE_init(:,:,iblk), &
                            uvelE     (:,:,iblk), vvelE     (:,:,iblk), &
                            TbE       (:,:,iblk))


            do i=1,nx_block
            do j=1,ny_block
               if (.not.iceUmask(i,j,iblk)) then
                  stresspU (i,j,iblk) = c0
                  stressmU (i,j,iblk) = c0
                  stress12U(i,j,iblk) = c0
               endif
               if (.not.iceTmask(i,j,iblk)) then
                  stresspT (i,j,iblk) = c0
                  stressmT (i,j,iblk) = c0
                  stress12T(i,j,iblk) = c0
               endif
            enddo
            enddo
         enddo  ! iblk
         !$OMP END PARALLEL DO

      endif ! grid_ice

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (grid_ice == 'CD' .or. grid_ice == 'C') then

         call ice_timer_start(timer_bound)
         call ice_HaloUpdate (uvelE,           halo_info, &
                              field_loc_Eface, field_type_vector)
         call ice_HaloUpdate (vvelN,           halo_info, &
                              field_loc_Nface, field_type_vector)
         call ice_timer_stop(timer_bound)

         if (grid_ice == 'C') then
            call grid_average_X2Y('A', uvelE, 'E', uvelN, 'N')
            call grid_average_X2Y('A', vvelN, 'N', vvelE, 'E')
            uvelN(:,:,:) = uvelN(:,:,:)*npm(:,:,:)
            vvelE(:,:,:) = vvelE(:,:,:)*epm(:,:,:)
         endif

         call ice_timer_start(timer_bound)
         call ice_HaloUpdate (uvelN,           halo_info, &
                              field_loc_Nface, field_type_vector)
         call ice_HaloUpdate (vvelE,           halo_info, &
                              field_loc_Eface, field_type_vector)
         call ice_timer_stop(timer_bound)

         call grid_average_X2Y('A', uvelE, 'E', uvel, 'U')
         call grid_average_X2Y('A', vvelN, 'N', vvel, 'U')
         uvel(:,:,:) = uvel(:,:,:)*uvm(:,:,:)
         vvel(:,:,:) = vvel(:,:,:)*uvm(:,:,:)
      endif

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (strength,         halo_info, &
                           field_loc_center, field_type_scalar)

      ! velocities may have changed in dyn_prep2
      call stack_fields(uvel, vvel, fld2)
      call ice_HaloUpdate (fld2,               halo_info, &
                           field_loc_NEcorner, field_type_vector)
      call unstack_fields(fld2, uvel, vvel)
      call ice_timer_stop(timer_bound)

      if (maskhalo_dyn) then
         halomask = 0
         if (grid_ice == 'B') then
            where (iceUmask) halomask = 1
         elseif (grid_ice == 'C' .or. grid_ice == 'CD') then
            !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block,i,j) SCHEDULE(runtime)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo,jhi
               do i = ilo,ihi
                  if (iceTmask(i  ,j  ,iblk) .or. &
                      iceTmask(i-1,j  ,iblk) .or. &
                      iceTmask(i+1,j  ,iblk) .or. &
                      iceTmask(i  ,j-1,iblk) .or. &
                      iceTmask(i  ,j+1,iblk)) then
                     halomask(i,j,iblk) = 1
                  endif
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO
         endif
         call ice_timer_start(timer_bound)
         call ice_HaloUpdate (halomask,         halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_timer_stop(timer_bound)
         call ice_HaloMask(halo_info_mask, halo_info, halomask)
      endif

      !-----------------------------------------------------------------
      ! seabed stress factor TbU (TbU is part of Cb coefficient)
      !-----------------------------------------------------------------

      if (seabed_stress) then

         if (grid_ice == "B") then

            if ( seabed_stress_method == 'LKD' ) then
               !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
               do iblk = 1, nblocks
                  call seabed_stress_factor_LKD (nx_block        , ny_block      , &
                                                 icellU    (iblk),                 &
                                                 indxUi  (:,iblk), indxUj(:,iblk), &
                                                 vice  (:,:,iblk), aice(:,:,iblk), &
                                                 hwater(:,:,iblk), TbU (:,:,iblk))
               enddo
               !$OMP END PARALLEL DO

            elseif ( seabed_stress_method == 'probabilistic' ) then
               !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
               do iblk = 1, nblocks
                  call seabed_stress_factor_prob (nx_block    , ny_block      ,                 &
                                                  icellT(iblk), indxTi(:,iblk), indxTj(:,iblk), &
                                                  icellU(iblk), indxUi(:,iblk), indxUj(:,iblk), &
                                                  aicen(:,:,:,iblk), vicen(:,:,:,iblk)        , &
                                                  hwater (:,:,iblk), TbU    (:,:,iblk))
               enddo
               !$OMP END PARALLEL DO
            endif

         elseif (grid_ice == "C" .or. grid_ice == "CD") then

            if ( seabed_stress_method == 'LKD' ) then
               !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
               do iblk = 1, nblocks
                  call seabed_stress_factor_LKD (nx_block        , ny_block,       &
                                                 icellE    (iblk),                 &
                                                 indxEi  (:,iblk), indxEj(:,iblk), &
                                                 vice  (:,:,iblk), aice(:,:,iblk), &
                                                 hwater(:,:,iblk), TbE (:,:,iblk), &
                                                 grid_location='E')
                  call seabed_stress_factor_LKD (nx_block        , ny_block,       &
                                                 icellN    (iblk),                 &
                                                 indxNi  (:,iblk), indxNj(:,iblk), &
                                                 vice  (:,:,iblk), aice(:,:,iblk), &
                                                 hwater(:,:,iblk), TbN (:,:,iblk), &
                                                 grid_location='N')
               enddo
               !$OMP END PARALLEL DO

            elseif ( seabed_stress_method == 'probabilistic' ) then
               !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
               do iblk = 1, nblocks
                  call seabed_stress_factor_prob (nx_block    , ny_block                      , &
                                                  icellT(iblk), indxTi(:,iblk), indxTj(:,iblk), &
                                                  icellU(iblk), indxUi(:,iblk), indxUj(:,iblk), &
                                                  aicen(:,:,:,iblk), vicen(:,:,:,iblk)        , &
                                                  hwater (:,:,iblk), TbU    (:,:,iblk)        , &
                                                  TbE    (:,:,iblk), TbN    (:,:,iblk)        , &
                                                  icellE(iblk), indxEi(:,iblk), indxEj(:,iblk), &
                                                  icellN(iblk), indxNi(:,iblk), indxNj(:,iblk)  )
               enddo
               !$OMP END PARALLEL DO
            endif

         endif

      endif

      call ice_timer_start(timer_evp)

      if (grid_ice == "B") then

         if (evp_algorithm == "shared_mem_1d" ) then

            call dyn_evp1d_run(stressp_1 , stressp_2, stressp_3  , stressp_4 , &
                               stressm_1 , stressm_2 , stressm_3 , stressm_4 , &
                               stress12_1, stress12_2, stress12_3, stress12_4, &
                               strength  ,                                     &
                               cdn_ocnU  , aiu       , uocnU     , vocnU     , &
                               waterxU   , wateryU   , forcexU   , forceyU   , &
                               umassdti  , fmU       , strintxU  , strintyU  , &
                               Tbu       , taubxU    , taubyU    , uvel      , &
                               vvel      , icetmask  , iceUmask)

         else ! evp_algorithm == standard_2d (Standard CICE)
            do ksub = 1,ndte        ! subcycling

               !$OMP PARALLEL DO PRIVATE(iblk,strtmp) SCHEDULE(runtime)
               do iblk = 1, nblocks

                  !-----------------------------------------------------------------
                  ! stress tensor equation, total surface stress
                  !-----------------------------------------------------------------
                  call stress (nx_block            , ny_block            , &
                                                     icellT        (iblk), &
                               indxTi      (:,iblk), indxTj      (:,iblk), &
                               uvel      (:,:,iblk), vvel      (:,:,iblk), &
                               dxT       (:,:,iblk), dyT       (:,:,iblk), &
                               dxhy      (:,:,iblk), dyhx      (:,:,iblk), &
                               cxp       (:,:,iblk), cyp       (:,:,iblk), &
                               cxm       (:,:,iblk), cym       (:,:,iblk), &
                                                     DminTarea (:,:,iblk), &
                               strength  (:,:,iblk),                       &
                               stressp_1 (:,:,iblk), stressp_2 (:,:,iblk), &
                               stressp_3 (:,:,iblk), stressp_4 (:,:,iblk), &
                               stressm_1 (:,:,iblk), stressm_2 (:,:,iblk), &
                               stressm_3 (:,:,iblk), stressm_4 (:,:,iblk), &
                               stress12_1(:,:,iblk), stress12_2(:,:,iblk), &
                               stress12_3(:,:,iblk), stress12_4(:,:,iblk), &
                               strtmp    (:,:,:))

                  !-----------------------------------------------------------------
                  ! momentum equation
                  !-----------------------------------------------------------------
                  call stepu (nx_block           , ny_block          , &
                              icellU       (iblk), Cdn_ocnU(:,:,iblk), &
                              indxUi     (:,iblk), indxUj    (:,iblk), &
                              aiU      (:,:,iblk), strtmp  (:,:,:),    &
                              uocnU    (:,:,iblk), vocnU   (:,:,iblk), &
                              waterxU  (:,:,iblk), wateryU (:,:,iblk), &
                              forcexU  (:,:,iblk), forceyU (:,:,iblk), &
                              umassdti (:,:,iblk), fmU     (:,:,iblk), &
                              uarear   (:,:,iblk),                     &
                              strintxU (:,:,iblk), strintyU(:,:,iblk), &
                              taubxU   (:,:,iblk), taubyU  (:,:,iblk), &
                              uvel_init(:,:,iblk), vvel_init(:,:,iblk),&
                              uvel     (:,:,iblk), vvel    (:,:,iblk), &
                              TbU      (:,:,iblk))

               enddo  ! iblk
               !$OMP END PARALLEL DO

               ! U fields at NE corner
               ! calls ice_haloUpdate, controls bundles and masks
               call dyn_haloUpdate (halo_info,          halo_info_mask,    &
                                    field_loc_NEcorner, field_type_vector, &
                                    uvel, vvel)

            enddo  ! sub cycling
         endif ! evp algorithm

         !-----------------------------------------------------------------
         ! save quantities for mechanical redistribution
         !-----------------------------------------------------------------

         !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
         do iblk = 1, nblocks
            call deformations (nx_block          , ny_block           , &
                               icellT      (iblk),                      &
                               indxTi    (:,iblk), indxTj     (:,iblk), &
                               uvel    (:,:,iblk), vvel     (:,:,iblk), &
                               dxT     (:,:,iblk), dyT      (:,:,iblk), &
                               dxU     (:,:,iblk), dyU      (:,:,iblk), &
                               cxp     (:,:,iblk), cyp      (:,:,iblk), &
                               cxm     (:,:,iblk), cym      (:,:,iblk), &
                               tarear  (:,:,iblk), vort     (:,:,iblk), &
                               shear   (:,:,iblk), divu     (:,:,iblk), &
                               rdg_conv(:,:,iblk), rdg_shear(:,:,iblk) )
         enddo
         !$OMP END PARALLEL DO

      elseif (grid_ice == "C") then

         do ksub = 1,ndte        ! subcycling

            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks

            !-----------------------------------------------------------------
            ! strain rates at U point
            ! NOTE these are actually strain rates * area  (m^2/s)
            !-----------------------------------------------------------------
               call strain_rates_U (nx_block          , ny_block           , &
                                    icellU      (iblk),                      &
                                    indxUi    (:,iblk), indxUj     (:,iblk), &
                                    uvelE   (:,:,iblk), vvelE    (:,:,iblk), &
                                    uvelN   (:,:,iblk), vvelN    (:,:,iblk), &
                                    uvel    (:,:,iblk), vvel     (:,:,iblk), &
                                    dxE     (:,:,iblk), dyN      (:,:,iblk), &
                                    dxU     (:,:,iblk), dyU      (:,:,iblk), &
                                    ratiodxN(:,:,iblk), ratiodxNr(:,:,iblk), &
                                    ratiodyE(:,:,iblk), ratiodyEr(:,:,iblk), &
                                    epm     (:,:,iblk), npm      (:,:,iblk), &
                                    divergU (:,:,iblk), tensionU (:,:,iblk), &
                                    shearU  (:,:,iblk), deltaU   (:,:,iblk)  )

            enddo  ! iblk
            !$OMP END PARALLEL DO

            ! calls ice_haloUpdate, controls bundles and masks
            call dyn_haloUpdate (halo_info,          halo_info_mask,    &
                                 field_loc_NEcorner, field_type_scalar, &
                                 shearU)

            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks
               call stressC_T (nx_block           , ny_block            , &
                                                    icellT        (iblk), &
                               indxTi      (:,iblk), indxTj      (:,iblk), &
                               uvelE     (:,:,iblk), vvelE     (:,:,iblk), &
                               uvelN     (:,:,iblk), vvelN     (:,:,iblk), &
                               dxN       (:,:,iblk), dyE       (:,:,iblk), &
                               dxT       (:,:,iblk), dyT       (:,:,iblk), &
                               uarea     (:,:,iblk), DminTarea (:,:,iblk), &
                               strength  (:,:,iblk), shearU    (:,:,iblk), &
                               zetax2T   (:,:,iblk), etax2T    (:,:,iblk), &
                               stresspT  (:,:,iblk), stressmT  (:,:,iblk), &
                               stress12T (:,:,iblk))

            enddo
            !$OMP END PARALLEL DO

            ! calls ice_haloUpdate, controls bundles and masks
            call dyn_haloUpdate (halo_info,        halo_info_mask,    &
                                 field_loc_center, field_type_scalar, &
                                 zetax2T, etax2T, stresspT, stressmT)

            if (visc_method == 'avg_strength') then
               call grid_average_X2Y('S', strength, 'T', strengthU, 'U')
            elseif (visc_method == 'avg_zeta') then
               call grid_average_X2Y('S', etax2T  , 'T', etax2U   , 'U')
            endif

            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks
               call stressC_U (nx_block           , ny_block            , &
                                                    icellU        (iblk), &
                              indxUi      (:,iblk), indxUj      (:,iblk), &
                              uarea     (:,:,iblk),                       &
                              etax2U    (:,:,iblk), deltaU    (:,:,iblk), &
                              strengthU (:,:,iblk), shearU    (:,:,iblk), &
                              stress12U (:,:,iblk))
            enddo
            !$OMP END PARALLEL DO

            ! calls ice_haloUpdate, controls bundles and masks
            call dyn_haloUpdate (halo_info         , halo_info_mask,    &
                                 field_loc_NEcorner, field_type_scalar, &
                                 stress12U)

            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks

               call div_stress_Ex (nx_block            , ny_block            , &
                                                         icellE        (iblk), &
                                   indxEi      (:,iblk), indxEj      (:,iblk), &
                                   dxE       (:,:,iblk), dyE       (:,:,iblk), &
                                   dxU       (:,:,iblk), dyT       (:,:,iblk), &
                                   earear    (:,:,iblk)                      , &
                                   stresspT  (:,:,iblk), stressmT  (:,:,iblk), &
                                   stress12U (:,:,iblk), strintxE  (:,:,iblk)  )

               call div_stress_Ny (nx_block            , ny_block            , &
                                                         icellN        (iblk), &
                                   indxNi      (:,iblk), indxNj      (:,iblk), &
                                   dxN       (:,:,iblk), dyN       (:,:,iblk), &
                                   dxT       (:,:,iblk), dyU       (:,:,iblk), &
                                   narear    (:,:,iblk)                      , &
                                   stresspT  (:,:,iblk), stressmT  (:,:,iblk), &
                                   stress12U (:,:,iblk), strintyN  (:,:,iblk)  )

            enddo
            !$OMP END PARALLEL DO

            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks

                call stepu_C (nx_block            , ny_block            , & ! u, E point
                              icellE        (iblk), Cdn_ocnE  (:,:,iblk), &
                              indxEi      (:,iblk), indxEj      (:,iblk), &
                                                    aiE       (:,:,iblk), &
                              uocnE     (:,:,iblk), vocnE     (:,:,iblk), &
                              waterxE   (:,:,iblk), forcexE   (:,:,iblk), &
                              emassdti  (:,:,iblk), fmE       (:,:,iblk), &
                              strintxE  (:,:,iblk), taubxE    (:,:,iblk), &
                              uvelE_init(:,:,iblk),                       &
                              uvelE     (:,:,iblk), vvelE     (:,:,iblk), &
                                 TbE       (:,:,iblk))

                call stepv_C (nx_block,             ny_block,             & ! v, N point
                              icellN        (iblk), Cdn_ocnN  (:,:,iblk), &
                              indxNi      (:,iblk), indxNj      (:,iblk), &
                                                    aiN       (:,:,iblk), &
                              uocnN     (:,:,iblk), vocnN     (:,:,iblk), &
                              wateryN   (:,:,iblk), forceyN   (:,:,iblk), &
                              nmassdti  (:,:,iblk), fmN       (:,:,iblk), &
                              strintyN  (:,:,iblk), taubyN    (:,:,iblk), &
                              vvelN_init(:,:,iblk),                       &
                              uvelN     (:,:,iblk), vvelN     (:,:,iblk), &
                              TbN       (:,:,iblk))
            enddo
            !$OMP END PARALLEL DO

            ! calls ice_haloUpdate, controls bundles and masks
            call dyn_haloUpdate (halo_info,       halo_info_mask,    &
                                 field_loc_Eface, field_type_vector, &
                                 uvelE)
            call dyn_haloUpdate (halo_info,       halo_info_mask,    &
                                 field_loc_Nface, field_type_vector, &
                                 vvelN)

            call grid_average_X2Y('A', uvelE, 'E', uvelN, 'N')
            call grid_average_X2Y('A', vvelN, 'N', vvelE, 'E')
            uvelN(:,:,:) = uvelN(:,:,:)*npm(:,:,:)
            vvelE(:,:,:) = vvelE(:,:,:)*epm(:,:,:)

            ! calls ice_haloUpdate, controls bundles and masks
            call dyn_haloUpdate (halo_info,       halo_info_mask,    &
                                 field_loc_Nface, field_type_vector, &
                                 uvelN)
            call dyn_haloUpdate (halo_info,       halo_info_mask,    &
                                 field_loc_Eface, field_type_vector, &
                                 vvelE)

            call grid_average_X2Y('A', uvelE, 'E', uvel, 'U')
            call grid_average_X2Y('A', vvelN, 'N', vvel, 'U')

            uvel(:,:,:) = uvel(:,:,:)*uvm(:,:,:)
            vvel(:,:,:) = vvel(:,:,:)*uvm(:,:,:)
            ! U fields at NE corner
            ! calls ice_haloUpdate, controls bundles and masks
            call dyn_haloUpdate (halo_info,          halo_info_mask,    &
                                 field_loc_NEcorner, field_type_vector, &
                                 uvel, vvel)

         enddo                     ! subcycling

         !-----------------------------------------------------------------
         ! save quantities for mechanical redistribution
         !-----------------------------------------------------------------

         !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
         do iblk = 1, nblocks
            call deformationsC_T (nx_block          , ny_block           , &
                                  icellT      (iblk),                      &
                                  indxTi    (:,iblk), indxTj     (:,iblk), &
                                  uvelE   (:,:,iblk), vvelE    (:,:,iblk), &
                                  uvelN   (:,:,iblk), vvelN    (:,:,iblk), &
                                  dxN     (:,:,iblk), dyE      (:,:,iblk), &
                                  dxT     (:,:,iblk), dyT      (:,:,iblk), &
                                  tarear  (:,:,iblk), uarea    (:,:,iblk), &
                                  shearU  (:,:,iblk), vort     (:,:,iblk), &
                                  shear   (:,:,iblk), divu     (:,:,iblk), &
                                  rdg_conv(:,:,iblk), rdg_shear(:,:,iblk))
         enddo
         !$OMP END PARALLEL DO

      elseif (grid_ice == "CD") then

         do ksub = 1,ndte        ! subcycling

            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks
               call stressCD_T (nx_block            , ny_block            , &
                                                      icellT        (iblk), &
                                indxTi      (:,iblk), indxTj      (:,iblk), &
                                uvelE     (:,:,iblk), vvelE     (:,:,iblk), &
                                uvelN     (:,:,iblk), vvelN     (:,:,iblk), &
                                dxN       (:,:,iblk), dyE       (:,:,iblk), &
                                dxT       (:,:,iblk), dyT       (:,:,iblk), &
                                                      DminTarea (:,:,iblk), &
                                strength  (:,:,iblk),                       &
                                zetax2T   (:,:,iblk), etax2T    (:,:,iblk), &
                                stresspT  (:,:,iblk), stressmT  (:,:,iblk), &
                                stress12T (:,:,iblk) )

            enddo
            !$OMP END PARALLEL DO

            ! calls ice_haloUpdate, controls bundles and masks
            call dyn_haloUpdate (halo_info,        halo_info_mask,    &
                                 field_loc_center, field_type_scalar, &
                                 zetax2T, etax2T)

            if (visc_method == 'avg_strength') then
               call grid_average_X2Y('S', strength, 'T', strengthU, 'U')
            elseif (visc_method == 'avg_zeta') then
               call grid_average_X2Y('S', zetax2T , 'T', zetax2U  , 'U')
               call grid_average_X2Y('S', etax2T  , 'T', etax2U   , 'U')
            endif

            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks
            !-----------------------------------------------------------------
            ! strain rates at U point
            ! NOTE these are actually strain rates * area  (m^2/s)
            !-----------------------------------------------------------------
               call strain_rates_U (nx_block           , ny_block           , &
                                                         icellU       (iblk), &
                                    indxUi     (:,iblk), indxUj     (:,iblk), &
                                    uvelE    (:,:,iblk), vvelE    (:,:,iblk), &
                                    uvelN    (:,:,iblk), vvelN    (:,:,iblk), &
                                    uvel     (:,:,iblk), vvel     (:,:,iblk), &
                                    dxE      (:,:,iblk), dyN      (:,:,iblk), &
                                    dxU      (:,:,iblk), dyU      (:,:,iblk), &
                                    ratiodxN (:,:,iblk), ratiodxNr(:,:,iblk), &
                                    ratiodyE (:,:,iblk), ratiodyEr(:,:,iblk), &
                                    epm      (:,:,iblk), npm      (:,:,iblk), &
                                    divergU  (:,:,iblk), tensionU (:,:,iblk), &
                                    shearU   (:,:,iblk), DeltaU   (:,:,iblk)  )

               call stressCD_U     (nx_block           , ny_block           , &
                                                         icellU       (iblk), &
                                    indxUi     (:,iblk), indxUj     (:,iblk), &
                                    uarea    (:,:,iblk),                      &
                                    zetax2U  (:,:,iblk), etax2U   (:,:,iblk), &
                                    strengthU(:,:,iblk),                      &
                                    divergU  (:,:,iblk), tensionU (:,:,iblk), &
                                    shearU   (:,:,iblk), DeltaU   (:,:,iblk), &
                                    stresspU (:,:,iblk), stressmU (:,:,iblk), &
                                    stress12U(:,:,iblk))
            enddo
            !$OMP END PARALLEL DO

            ! calls ice_haloUpdate, controls bundles and masks
            call dyn_haloUpdate (halo_info,         halo_info_mask,    &
                                 field_loc_center,  field_type_scalar, &
                                 stresspT, stressmT, stress12T)
            call dyn_haloUpdate (halo_info,         halo_info_mask,    &
                                 field_loc_NEcorner,field_type_scalar, &
                                 stresspU, stressmU, stress12U)

            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks

               call div_stress_Ex (nx_block            , ny_block            , &
                                                         icellE        (iblk), &
                                   indxEi      (:,iblk), indxEj      (:,iblk), &
                                   dxE       (:,:,iblk), dyE       (:,:,iblk), &
                                   dxU       (:,:,iblk), dyT       (:,:,iblk), &
                                   earear    (:,:,iblk)                      , &
                                   stresspT  (:,:,iblk), stressmT  (:,:,iblk), &
                                   stress12U (:,:,iblk), strintxE  (:,:,iblk)  )

               call div_stress_Ey (nx_block            , ny_block            , &
                                                         icellE        (iblk), &
                                   indxEi      (:,iblk), indxEj      (:,iblk), &
                                   dxE       (:,:,iblk), dyE       (:,:,iblk), &
                                   dxU       (:,:,iblk), dyT       (:,:,iblk), &
                                   earear    (:,:,iblk)                      , &
                                   stresspU  (:,:,iblk), stressmU  (:,:,iblk), &
                                   stress12T (:,:,iblk), strintyE  (:,:,iblk)  )

               call div_stress_Nx (nx_block            , ny_block            , &
                                                         icellN        (iblk), &
                                   indxNi      (:,iblk), indxNj      (:,iblk), &
                                   dxN       (:,:,iblk), dyN       (:,:,iblk), &
                                   dxT       (:,:,iblk), dyU       (:,:,iblk), &
                                   narear    (:,:,iblk)                      , &
                                   stresspU  (:,:,iblk), stressmU  (:,:,iblk), &
                                   stress12T (:,:,iblk), strintxN  (:,:,iblk)  )

               call div_stress_Ny (nx_block            , ny_block            , &
                                                         icellN        (iblk), &
                                   indxNi      (:,iblk), indxNj      (:,iblk), &
                                   dxN       (:,:,iblk), dyN       (:,:,iblk), &
                                   dxT       (:,:,iblk), dyU       (:,:,iblk), &
                                   narear    (:,:,iblk)                      , &
                                   stresspT  (:,:,iblk), stressmT  (:,:,iblk), &
                                   stress12U (:,:,iblk), strintyN  (:,:,iblk)  )

            enddo
            !$OMP END PARALLEL DO

            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks

               call stepuv_CD (nx_block            , ny_block            , & ! E point
                               icellE        (iblk), Cdn_ocnE  (:,:,iblk), &
                               indxEi      (:,iblk), indxEj      (:,iblk), &
                                                     aiE       (:,:,iblk), &
                               uocnE     (:,:,iblk), vocnE     (:,:,iblk), &
                               waterxE   (:,:,iblk), wateryE   (:,:,iblk), &
                               forcexE   (:,:,iblk), forceyE   (:,:,iblk), &
                               emassdti  (:,:,iblk), fmE       (:,:,iblk), &
                               strintxE  (:,:,iblk), strintyE  (:,:,iblk), &
                               taubxE    (:,:,iblk), taubyE    (:,:,iblk), &
                               uvelE_init(:,:,iblk), vvelE_init(:,:,iblk), &
                               uvelE     (:,:,iblk), vvelE     (:,:,iblk), &
                               TbE       (:,:,iblk))

               call stepuv_CD (nx_block            , ny_block            , & ! N point
                               icellN        (iblk), Cdn_ocnN  (:,:,iblk), &
                               indxNi      (:,iblk), indxNj      (:,iblk), &
                                                     aiN       (:,:,iblk), &
                               uocnN     (:,:,iblk), vocnN     (:,:,iblk), &
                               waterxN   (:,:,iblk), wateryN   (:,:,iblk), &
                               forcexN   (:,:,iblk), forceyN   (:,:,iblk), &
                               nmassdti  (:,:,iblk), fmN       (:,:,iblk), &
                               strintxN  (:,:,iblk), strintyN  (:,:,iblk), &
                               taubxN    (:,:,iblk), taubyN    (:,:,iblk), &
                               uvelN_init(:,:,iblk), vvelN_init(:,:,iblk), &
                               uvelN     (:,:,iblk), vvelN     (:,:,iblk), &
                               TbN       (:,:,iblk))
            enddo
            !$OMP END PARALLEL DO

            ! calls ice_haloUpdate, controls bundles and masks
            call dyn_haloUpdate (halo_info,       halo_info_mask,    &
                                 field_loc_Eface, field_type_vector, &
                                 uvelE, vvelE)
            call dyn_haloUpdate (halo_info,       halo_info_mask,    &
                                 field_loc_Nface, field_type_vector, &
                                 uvelN, vvelN)

            call grid_average_X2Y('A', uvelE, 'E', uvel, 'U')
            call grid_average_X2Y('A', vvelN, 'N', vvel, 'U')

            uvel(:,:,:) = uvel(:,:,:)*uvm(:,:,:)
            vvel(:,:,:) = vvel(:,:,:)*uvm(:,:,:)
               ! U fields at NE corner
               ! calls ice_haloUpdate, controls bundles and masks
            call dyn_haloUpdate (halo_info,          halo_info_mask,    &
                                 field_loc_NEcorner, field_type_vector, &
                                 uvel, vvel)

         enddo                     ! subcycling

         !-----------------------------------------------------------------
         ! save quantities for mechanical redistribution
         !-----------------------------------------------------------------

         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks
            call deformationsCD_T (nx_block          , ny_block           , &
                                   icellT      (iblk),                      &
                                   indxTi    (:,iblk), indxTj     (:,iblk), &
                                   uvelE   (:,:,iblk), vvelE    (:,:,iblk), &
                                   uvelN   (:,:,iblk), vvelN    (:,:,iblk), &
                                   dxN     (:,:,iblk), dyE      (:,:,iblk), &
                                   dxT     (:,:,iblk), dyT      (:,:,iblk), &
                                   tarear  (:,:,iblk), vort     (:,:,iblk), &
                                   shear   (:,:,iblk), divu     (:,:,iblk), &
                                   rdg_conv(:,:,iblk), rdg_shear(:,:,iblk))
         enddo
         !$OMP END PARALLEL DO
      endif   ! grid_ice

      call ice_timer_stop(timer_evp)

      if (maskhalo_dyn) then
         call ice_HaloDestroy(halo_info_mask)
      endif

      ! Force symmetry across the tripole seam
      if (trim(ns_boundary_type) == 'tripole' .or. &
          trim(ns_boundary_type) == 'tripoleT') then
         ! TODO: C/CD-grid
      if (maskhalo_dyn) then
         !-------------------------------------------------------
         ! set halomask to zero because ice_HaloMask always keeps
         ! local copies AND tripole zipper communication
         !-------------------------------------------------------
         halomask = 0
         call ice_HaloMask(halo_info_mask, halo_info, halomask)

         call ice_HaloUpdate_stress(stressp_1 , stressp_3 , halo_info_mask, &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_3 , stressp_1 , halo_info_mask, &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_2 , stressp_4 , halo_info_mask, &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_4 , stressp_2 , halo_info_mask, &
                                    field_loc_center,  field_type_scalar)

         call ice_HaloUpdate_stress(stressm_1 , stressm_3 , halo_info_mask, &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_3 , stressm_1 , halo_info_mask, &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_2 , stressm_4 , halo_info_mask, &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_4 , stressm_2 , halo_info_mask, &
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
         call ice_HaloUpdate_stress(stressp_1 , stressp_3 , halo_info,  &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_3 , stressp_1 , halo_info,  &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_2 , stressp_4 , halo_info,  &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_4 , stressp_2 , halo_info,  &
                                    field_loc_center,  field_type_scalar)

         call ice_HaloUpdate_stress(stressm_1 , stressm_3 , halo_info,  &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_3 , stressm_1 , halo_info,  &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_2 , stressm_4 , halo_info,  &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_4 , stressm_2 , halo_info,  &
                                    field_loc_center,  field_type_scalar)

         call ice_HaloUpdate_stress(stress12_1, stress12_3, halo_info,  &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_3, stress12_1, halo_info,  &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_2, stress12_4, halo_info,  &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_4, stress12_2, halo_info,  &
                                    field_loc_center,  field_type_scalar)
      endif   ! maskhalo
      endif   ! tripole

      !-----------------------------------------------------------------
      ! ice-ocean stress
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
      do iblk = 1, nblocks
         call dyn_finish                               &
              (nx_block          , ny_block          , &
               icellU      (iblk), Cdn_ocnU(:,:,iblk), &
               indxUi    (:,iblk), indxUj    (:,iblk), &
               uvel    (:,:,iblk), vvel    (:,:,iblk), &
               uocnU   (:,:,iblk), vocnU   (:,:,iblk), &
               aiU     (:,:,iblk), fmU     (:,:,iblk), &
               strocnxU(:,:,iblk), strocnyU(:,:,iblk))
      enddo
      !$OMP END PARALLEL DO

      if (grid_ice == 'CD' .or. grid_ice == 'C') then

         !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
         do iblk = 1, nblocks

            call dyn_finish                               &
                 (nx_block          , ny_block          , &
                  icellN      (iblk), Cdn_ocnN(:,:,iblk), &
                  indxNi    (:,iblk), indxNj    (:,iblk), &
                  uvelN   (:,:,iblk), vvelN   (:,:,iblk), &
                  uocnN   (:,:,iblk), vocnN   (:,:,iblk), &
                  aiN     (:,:,iblk), fmN     (:,:,iblk), &
                  strocnxN(:,:,iblk), strocnyN(:,:,iblk))

            call dyn_finish                               &
                 (nx_block          , ny_block          , &
                  icellE      (iblk), Cdn_ocnE(:,:,iblk), &
                  indxEi    (:,iblk), indxEj    (:,iblk), &
                  uvelE   (:,:,iblk), vvelE   (:,:,iblk), &
                  uocnE   (:,:,iblk), vocnE   (:,:,iblk), &
                  aiE     (:,:,iblk), fmE     (:,:,iblk), &
                  strocnxE(:,:,iblk), strocnyE(:,:,iblk))

         enddo
         !$OMP END PARALLEL DO

      endif

      if (grid_ice == 'CD' .or. grid_ice == 'C') then
         call ice_HaloUpdate (strintxE,        halo_info, &
                              field_loc_Eface, field_type_vector)
         call ice_HaloUpdate (strintyN,        halo_info, &
                              field_loc_Nface, field_type_vector)
         call ice_timer_stop(timer_bound)
         call grid_average_X2Y('S', strintxE, 'E', strintxU, 'U')    ! diagnostic
         call grid_average_X2Y('S', strintyN, 'N', strintyU, 'U')    ! diagnostic
      endif

      call ice_timer_stop(timer_dynamics)    ! dynamics

      end subroutine evp

!=======================================================================
! Computes the rates of strain and internal stress components for
! each of the four corners on each T-grid cell.
! Computes stress terms for the momentum equation
!
! author: Elizabeth C. Hunke, LANL

      subroutine stress (nx_block,   ny_block,   &
                         icellT,                 &
                         indxTi,     indxTj,     &
                         uvel,       vvel,       &
                         dxT,        dyT,        &
                         dxhy,       dyhx,       &
                         cxp,        cyp,        &
                         cxm,        cym,        &
                                     DminTarea,  &
                         strength,               &
                         stressp_1,  stressp_2,  &
                         stressp_3,  stressp_4,  &
                         stressm_1,  stressm_2,  &
                         stressm_3,  stressm_4,  &
                         stress12_1, stress12_2, &
                         stress12_3, stress12_4, &
                         str )

      use ice_dyn_shared, only: strain_rates, visc_replpress

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellT                ! no. of cells where iceTmask = .true.

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxTi   , & ! compressed index in i-direction
         indxTj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         strength , & ! ice strength (N/m)
         uvel     , & ! x-component of velocity (m/s)
         vvel     , & ! y-component of velocity (m/s)
         dxT      , & ! width of T-cell through the middle (m)
         dyT      , & ! height of T-cell through the middle (m)
         dxhy     , & ! 0.5*(HTE - HTW)
         dyhx     , & ! 0.5*(HTN - HTS)
         cyp      , & ! 1.5*HTE - 0.5*HTW
         cxp      , & ! 1.5*HTN - 0.5*HTS
         cym      , & ! 0.5*HTE - 1.5*HTW
         cxm      , & ! 0.5*HTN - 1.5*HTS
         DminTarea    ! deltaminEVP*tarea

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         stressp_1, stressp_2, stressp_3, stressp_4 , & ! sigma11+sigma22
         stressm_1, stressm_2, stressm_3, stressm_4 , & ! sigma11-sigma22
         stress12_1,stress12_2,stress12_3,stress12_4    ! sigma12

      real (kind=dbl_kind), dimension(nx_block,ny_block,8), intent(out) :: &
         str          ! stress combinations

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &
        divune, divunw, divuse, divusw            , & ! divergence
        tensionne, tensionnw, tensionse, tensionsw, & ! tension
        shearne, shearnw, shearse, shearsw        , & ! shearing
        Deltane, Deltanw, Deltase, Deltasw        , & ! Delt
        zetax2ne, zetax2nw, zetax2se, zetax2sw    , & ! 2 x zeta (bulk visc)
        etax2ne, etax2nw, etax2se, etax2sw        , & ! 2 x eta (shear visc)
        rep_prsne, rep_prsnw, rep_prsse, rep_prssw, & ! replacement pressure
!        puny                                      , & ! puny
        ssigpn, ssigps, ssigpe, ssigpw            , &
        ssigmn, ssigms, ssigme, ssigmw            , &
        ssig12n, ssig12s, ssig12e, ssig12w        , &
        ssigp1, ssigp2, ssigm1, ssigm2, ssig121, ssig122, &
        csigpne, csigpnw, csigpse, csigpsw        , &
        csigmne, csigmnw, csigmse, csigmsw        , &
        csig12ne, csig12nw, csig12se, csig12sw    , &
        str12ew, str12we, str12ns, str12sn        , &
        strp_tmp, strm_tmp

      character(len=*), parameter :: subname = '(stress)'

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      str(:,:,:) = c0

      do ij = 1, icellT
         i = indxTi(ij)
         j = indxTj(ij)

         !-----------------------------------------------------------------
         ! strain rates
         ! NOTE these are actually strain rates * area  (m^2/s)
         !-----------------------------------------------------------------

         call strain_rates (nx_block,   ny_block,   &
                            i,          j,          &
                            uvel,       vvel,       &
                            dxT,        dyT,        &
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
         ! viscosities and replacement pressure
         !-----------------------------------------------------------------

         call visc_replpress (strength(i,j), DminTarea(i,j), Deltane, &
                              zetax2ne, etax2ne, rep_prsne)

         call visc_replpress (strength(i,j), DminTarea(i,j), Deltanw, &
                              zetax2nw, etax2nw, rep_prsnw)

         call visc_replpress (strength(i,j), DminTarea(i,j), Deltasw, &
                              zetax2sw, etax2sw, rep_prssw)

         call visc_replpress (strength(i,j), DminTarea(i,j), Deltase, &
                              zetax2se, etax2se, rep_prsse)

         !-----------------------------------------------------------------
         ! the stresses                            ! kg/s^2
         ! (1) northeast, (2) northwest, (3) southwest, (4) southeast
         !-----------------------------------------------------------------

         ! NOTE: for comp. efficiency 2 x zeta and 2 x eta are used in the code
         stressp_1 (i,j) = (stressp_1 (i,j)*(c1-arlx1i*revp) &
                           + arlx1i*(zetax2ne*divune - rep_prsne)) * denom1
         stressp_2 (i,j) = (stressp_2 (i,j)*(c1-arlx1i*revp) &
                           + arlx1i*(zetax2nw*divunw - rep_prsnw)) * denom1
         stressp_3 (i,j) = (stressp_3 (i,j)*(c1-arlx1i*revp) &
                           + arlx1i*(zetax2sw*divusw - rep_prssw)) * denom1
         stressp_4 (i,j) = (stressp_4 (i,j)*(c1-arlx1i*revp) &
                           + arlx1i*(zetax2se*divuse - rep_prsse)) * denom1

         stressm_1 (i,j) = (stressm_1 (i,j)*(c1-arlx1i*revp) &
                           + arlx1i*etax2ne*tensionne) * denom1
         stressm_2 (i,j) = (stressm_2 (i,j)*(c1-arlx1i*revp) &
                           + arlx1i*etax2nw*tensionnw) * denom1
         stressm_3 (i,j) = (stressm_3 (i,j)*(c1-arlx1i*revp) &
                           + arlx1i*etax2sw*tensionsw) * denom1
         stressm_4 (i,j) = (stressm_4 (i,j)*(c1-arlx1i*revp) &
                           + arlx1i*etax2se*tensionse) * denom1

         stress12_1(i,j) = (stress12_1(i,j)*(c1-arlx1i*revp) &
                           + arlx1i*p5*etax2ne*shearne) * denom1
         stress12_2(i,j) = (stress12_2(i,j)*(c1-arlx1i*revp) &
                           + arlx1i*p5*etax2nw*shearnw) * denom1
         stress12_3(i,j) = (stress12_3(i,j)*(c1-arlx1i*revp) &
                           + arlx1i*p5*etax2sw*shearsw) * denom1
         stress12_4(i,j) = (stress12_4(i,j)*(c1-arlx1i*revp) &
                           + arlx1i*p5*etax2se*shearse) * denom1

         !-----------------------------------------------------------------
         ! Eliminate underflows.
         ! The following code is commented out because it is relatively
         ! expensive and most compilers include a flag that accomplishes
         ! the same thing more efficiently.  This code is cheaper than
         ! handling underflows if the compiler lacks a flag; uncomment
         ! it in that case.  The compiler flag is often described with the
         ! phrase "flush to zero".
         !-----------------------------------------------------------------

!         call icepack_query_parameters(puny_out=puny)
!         call icepack_warnings_flush(nu_diag)
!         if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
!            file=__FILE__, line=__LINE__)

!         stressp_1(i,j) = sign(max(abs(stressp_1(i,j)),puny),stressp_1(i,j))
!         stressp_2(i,j) = sign(max(abs(stressp_2(i,j)),puny),stressp_2(i,j))
!         stressp_3(i,j) = sign(max(abs(stressp_3(i,j)),puny),stressp_3(i,j))
!         stressp_4(i,j) = sign(max(abs(stressp_4(i,j)),puny),stressp_4(i,j))

!         stressm_1(i,j) = sign(max(abs(stressm_1(i,j)),puny),stressm_1(i,j))
!         stressm_2(i,j) = sign(max(abs(stressm_2(i,j)),puny),stressm_2(i,j))
!         stressm_3(i,j) = sign(max(abs(stressm_3(i,j)),puny),stressm_3(i,j))
!         stressm_4(i,j) = sign(max(abs(stressm_4(i,j)),puny),stressm_4(i,j))

!         stress12_1(i,j) = sign(max(abs(stress12_1(i,j)),puny),stress12_1(i,j))
!         stress12_2(i,j) = sign(max(abs(stress12_2(i,j)),puny),stress12_2(i,j))
!         stress12_3(i,j) = sign(max(abs(stress12_3(i,j)),puny),stress12_3(i,j))
!         stress12_4(i,j) = sign(max(abs(stress12_4(i,j)),puny),stress12_4(i,j))

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

         str12ew = p5*dxT(i,j)*(p333*ssig12e + p166*ssig12w)
         str12we = p5*dxT(i,j)*(p333*ssig12w + p166*ssig12e)
         str12ns = p5*dyT(i,j)*(p333*ssig12n + p166*ssig12s)
         str12sn = p5*dyT(i,j)*(p333*ssig12s + p166*ssig12n)

         !-----------------------------------------------------------------
         ! for dF/dx (u momentum)
         !-----------------------------------------------------------------
         strp_tmp  = p25*dyT(i,j)*(p333*ssigpn  + p166*ssigps)
         strm_tmp  = p25*dyT(i,j)*(p333*ssigmn  + p166*ssigms)

         ! northeast (i,j)
         str(i,j,1) = -strp_tmp - strm_tmp - str12ew &
              + dxhy(i,j)*(-csigpne + csigmne) + dyhx(i,j)*csig12ne

         ! northwest (i+1,j)
         str(i,j,2) = strp_tmp + strm_tmp - str12we &
              + dxhy(i,j)*(-csigpnw + csigmnw) + dyhx(i,j)*csig12nw

         strp_tmp  = p25*dyT(i,j)*(p333*ssigps  + p166*ssigpn)
         strm_tmp  = p25*dyT(i,j)*(p333*ssigms  + p166*ssigmn)

         ! southeast (i,j+1)
         str(i,j,3) = -strp_tmp - strm_tmp + str12ew &
              + dxhy(i,j)*(-csigpse + csigmse) + dyhx(i,j)*csig12se

         ! southwest (i+1,j+1)
         str(i,j,4) = strp_tmp + strm_tmp + str12we &
              + dxhy(i,j)*(-csigpsw + csigmsw) + dyhx(i,j)*csig12sw

         !-----------------------------------------------------------------
         ! for dF/dy (v momentum)
         !-----------------------------------------------------------------
         strp_tmp  = p25*dxT(i,j)*(p333*ssigpe  + p166*ssigpw)
         strm_tmp  = p25*dxT(i,j)*(p333*ssigme  + p166*ssigmw)

         ! northeast (i,j)
         str(i,j,5) = -strp_tmp + strm_tmp - str12ns &
              - dyhx(i,j)*(csigpne + csigmne) + dxhy(i,j)*csig12ne

         ! southeast (i,j+1)
         str(i,j,6) = strp_tmp - strm_tmp - str12sn &
              - dyhx(i,j)*(csigpse + csigmse) + dxhy(i,j)*csig12se

         strp_tmp  = p25*dxT(i,j)*(p333*ssigpw  + p166*ssigpe)
         strm_tmp  = p25*dxT(i,j)*(p333*ssigmw  + p166*ssigme)

         ! northwest (i+1,j)
         str(i,j,7) = -strp_tmp + strm_tmp + str12ns &
              - dyhx(i,j)*(csigpnw + csigmnw) + dxhy(i,j)*csig12nw

         ! southwest (i+1,j+1)
         str(i,j,8) = strp_tmp - strm_tmp + str12sn &
              - dyhx(i,j)*(csigpsw + csigmsw) + dxhy(i,j)*csig12sw

      enddo                     ! ij

      end subroutine stress

!=======================================================================
! Computes the strain rates and internal stress components for C grid
!
! author: JF Lemieux, ECCC
! updated: D. Bailey, NCAR
! Nov 2021
!
! Bouillon, S., T. Fichefet, V. Legat and G. Madec (2013). The
! elastic-viscous-plastic method revisited. Ocean Model., 71, 2-12.
!
! Kimmritz, M., S. Danilov and M. Losch (2016). The adaptive EVP method
! for solving the sea ice momentum equation. Ocean Model., 101, 59-67.

      subroutine stressC_T  (nx_block,    ny_block  , &
                                          icellT    , &
                             indxTi     , indxTj    , &
                             uvelE      , vvelE     , &
                             uvelN      , vvelN     , &
                             dxN        , dyE       , &
                             dxT        , dyT       , &
                             uarea      , DminTarea , &
                             strength   , shearU    , &
                             zetax2T    , etax2T    , &
                             stresspT   , stressmT  , &
                             stress12T)

      use ice_dyn_shared, only: strain_rates_T, &
                                visc_replpress, e_factor

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellT                ! no. of cells where iceTmask = .true.

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxTi   , & ! compressed index in i-direction
         indxTj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         uvelE    , & ! x-component of velocity (m/s) at the E point
         vvelE    , & ! y-component of velocity (m/s) at the E point
         uvelN    , & ! x-component of velocity (m/s) at the N point
         vvelN    , & ! y-component of velocity (m/s) at the N point
         dxN      , & ! width of N-cell through the middle (m)
         dyE      , & ! height of E-cell through the middle (m)
         dxT      , & ! width of T-cell through the middle (m)
         dyT      , & ! height of T-cell through the middle (m)
         strength , & ! ice strength (N/m)
         shearU   , & ! shearU local for this routine
         uarea    , & ! area of u cell
         DminTarea    ! deltaminEVP*tarea

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         zetax2T  , & ! zetax2 = 2*zeta (bulk viscosity)
         etax2T   , & ! etax2  = 2*eta  (shear viscosity)
         stresspT , & ! sigma11+sigma22
         stressmT , & ! sigma11-sigma22
         stress12T    ! sigma12

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
        divT      , & ! divergence at T point
        tensionT      ! tension at T point

      real (kind=dbl_kind) :: &
        shearTsqr , & ! strain rates squared at T point
        shearT    , & ! strain rate at T point
        DeltaT    , & ! delt at T point
        uareaavgr , & ! 1 / uarea avg
        rep_prsT      ! replacement pressure at T point

      character(len=*), parameter :: subname = '(stressC_T)'

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      call strain_rates_T (nx_block   , ny_block     , &
                           icellT     ,                &
                           indxTi(:)  , indxTj  (:)  , &
                           uvelE (:,:), vvelE   (:,:), &
                           uvelN (:,:), vvelN   (:,:), &
                           dxN   (:,:), dyE     (:,:), &
                           dxT   (:,:), dyT     (:,:), &
                           divT  (:,:), tensionT(:,:))

      do ij = 1, icellT
         i = indxTi(ij)
         j = indxTj(ij)

         !-----------------------------------------------------------------
         ! Square of shear strain rate at T obtained from interpolation of
         ! U point values (Bouillon et al., 2013, Kimmritz et al., 2016
         !-----------------------------------------------------------------

         uareaavgr = c1/(uarea(i,j)+uarea(i,j-1)+uarea(i-1,j-1)+uarea(i-1,j))

         shearTsqr = (shearU(i  ,j  )**2 * uarea(i  ,j  )  &
                    + shearU(i  ,j-1)**2 * uarea(i  ,j-1)  &
                    + shearU(i-1,j-1)**2 * uarea(i-1,j-1)  &
                    + shearU(i-1,j  )**2 * uarea(i-1,j  )) &
                    * uareaavgr

         shearT    = (shearU(i  ,j  )    * uarea(i  ,j  )  &
                    + shearU(i  ,j-1)    * uarea(i  ,j-1)  &
                    + shearU(i-1,j-1)    * uarea(i-1,j-1)  &
                    + shearU(i-1,j  )    * uarea(i-1,j  )) &
                    * uareaavgr

         DeltaT = sqrt(divT(i,j)**2 + e_factor*(tensionT(i,j)**2 + shearTsqr))

         !-----------------------------------------------------------------
         ! viscosities and replacement pressure at T point
         !-----------------------------------------------------------------

         call visc_replpress (strength(i,j), DminTarea(i,j), DeltaT, &
                              zetax2T (i,j), etax2T(i,j), rep_prsT)

         !-----------------------------------------------------------------
         ! the stresses                            ! kg/s^2
         !-----------------------------------------------------------------

         ! NOTE: for comp. efficiency 2 x zeta and 2 x eta are used in the code

         stresspT(i,j)  = (stresspT (i,j)*(c1-arlx1i*revp) &
                          + arlx1i*(zetax2T(i,j)*divT(i,j) - rep_prsT)) * denom1

         stressmT(i,j)  = (stressmT (i,j)*(c1-arlx1i*revp) &
                          + arlx1i*etax2T(i,j)*tensionT(i,j)) * denom1

         stress12T(i,j) = (stress12T(i,j)*(c1-arlx1i*revp) &
                          + arlx1i*p5*etax2T(i,j)*shearT    ) * denom1

      enddo                     ! ij

      end subroutine stressC_T

!=======================================================================
!
! Computes the strain rates and internal stress components for U points
!
! author: JF Lemieux, ECCC
! Nov 2021
!
! Bouillon, S., T. Fichefet, V. Legat and G. Madec (2013). The
! elastic-viscous-plastic method revisited. Ocean Model., 71, 2-12.
!
! Kimmritz, M., S. Danilov and M. Losch (2016). The adaptive EVP method
! for solving the sea ice momentum equation. Ocean Model., 101, 59-67.

      subroutine stressC_U  (nx_block , ny_block  ,&
                                        icellU    ,&
                             indxUi   , indxUj    ,&
                             uarea    ,            &
                             etax2U   , deltaU    ,&
                             strengthU, shearU    ,&
                             stress12U)

      use ice_dyn_shared, only: visc_replpress, &
                                visc_method, deltaminEVP

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellU                ! no. of cells where iceUmask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxUi   , & ! compressed index in i-direction
         indxUj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         uarea    , & ! area of U point
         etax2U   , & ! 2*eta at the U point
         shearU   , & ! shearU array
         deltaU   , & ! deltaU array
         strengthU    ! ice strength at the U point

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         stress12U    ! sigma12

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &
         lzetax2U , & ! bulk viscosity at U point
         letax2U  , & ! shear viscosity at U point
         lrep_prsU, & ! replacement pressure at U point
         DminUarea    ! Dmin on U

      character(len=*), parameter :: subname = '(stressC_U)'

      !-----------------------------------------------------------------
      ! viscosities and replacement pressure at U point
      ! avg_zeta: Bouillon et al. 2013, C1 method of Kimmritz et al. 2016
      ! avg_strength: C2 method of Kimmritz et al. 2016
      ! if outside do and stress12U equation repeated in each loop for performance
      !-----------------------------------------------------------------

      if (visc_method == 'avg_zeta') then
         do ij = 1, icellU
            i = indxUi(ij)
            j = indxUj(ij)
            stress12U(i,j) = (stress12U(i,j)*(c1-arlx1i*revp) &
                             + arlx1i*p5*etax2U(i,j)*shearU(i,j)) * denom1
         enddo

      elseif (visc_method == 'avg_strength') then
         do ij = 1, icellU
            i = indxUi(ij)
            j = indxUj(ij)
            DminUarea = deltaminEVP*uarea(i,j)
            ! only need etax2U here, but other terms are calculated with etax2U
            ! minimal extra calculations here even though it seems like there is
            call visc_replpress (strengthU(i,j), DminUarea, deltaU(i,j), &
                                 lzetax2U      , letax2U  , lrep_prsU)
            stress12U(i,j) = (stress12U(i,j)*(c1-arlx1i*revp) &
                             + arlx1i*p5*letax2U*shearU(i,j)) * denom1
         enddo

      endif

      end subroutine stressC_U

!=======================================================================
! Computes the strain rates and internal stress components for T points
!
! author: JF Lemieux, ECCC
! Nov 2021

      subroutine stressCD_T (nx_block,   ny_block , &
                                         icellT   , &
                             indxTi  ,   indxTj   , &
                             uvelE   ,   vvelE    , &
                             uvelN   ,   vvelN    , &
                             dxN     ,   dyE      , &
                             dxT     ,   dyT      , &
                                         DminTarea, &
                             strength,              &
                             zetax2T ,   etax2T   , &
                             stresspT,   stressmT , &
                             stress12T)

      use ice_dyn_shared, only: strain_rates_T, &
                                visc_replpress

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellT                ! no. of cells where iceTmask = .true.

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxTi   , & ! compressed index in i-direction
         indxTj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         uvelE    , & ! x-component of velocity (m/s) at the E point
         vvelE    , & ! y-component of velocity (m/s) at the N point
         uvelN    , & ! x-component of velocity (m/s) at the E point
         vvelN    , & ! y-component of velocity (m/s) at the N point
         dxN      , & ! width of N-cell through the middle (m)
         dyE      , & ! height of E-cell through the middle (m)
         dxT      , & ! width of T-cell through the middle (m)
         dyT      , & ! height of T-cell through the middle (m)
         strength , & ! ice strength (N/m)
         DminTarea    ! deltaminEVP*tarea

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         zetax2T  , & ! zetax2 = 2*zeta (bulk viscosity)
         etax2T   , & ! etax2  = 2*eta  (shear viscosity)
         stresspT , & ! sigma11+sigma22
         stressmT , & ! sigma11-sigma22
         stress12T    ! sigma12

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
        divT      , & ! divergence at T point
        tensionT  , & ! tension at T point
        shearT    , & ! shear at T point
        DeltaT        ! delt at T point

      real (kind=dbl_kind) :: &
        rep_prsT      ! replacement pressure at T point

      character(len=*), parameter :: subname = '(stressCD_T)'

      !-----------------------------------------------------------------
      ! strain rates at T point
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------

      call strain_rates_T (nx_block   , ny_block     , &
                           icellT     ,                &
                           indxTi(:)  , indxTj  (:)  , &
                           uvelE (:,:), vvelE   (:,:), &
                           uvelN (:,:), vvelN   (:,:), &
                           dxN   (:,:), dyE     (:,:), &
                           dxT   (:,:), dyT     (:,:), &
                           divT  (:,:), tensionT(:,:), &
                           shearT(:,:), DeltaT  (:,:))

      do ij = 1, icellT
         i = indxTi(ij)
         j = indxTj(ij)

         !-----------------------------------------------------------------
         ! viscosities and replacement pressure at T point
         !-----------------------------------------------------------------

         call visc_replpress (strength(i,j), DminTarea(i,j), DeltaT(i,j), &
                              zetax2T (i,j), etax2T(i,j), rep_prsT)

         !-----------------------------------------------------------------
         ! the stresses                            ! kg/s^2
         !-----------------------------------------------------------------

         ! NOTE: for comp. efficiency 2 x zeta and 2 x eta are used in the code

         stresspT(i,j)  = (stresspT (i,j)*(c1-arlx1i*revp) &
                          + arlx1i*(zetax2T(i,j)*divT(i,j) - rep_prsT)) * denom1

         stressmT(i,j)  = (stressmT (i,j)*(c1-arlx1i*revp) &
                          + arlx1i*etax2T(i,j)*tensionT(i,j)) * denom1

         stress12T(i,j) = (stress12T(i,j)*(c1-arlx1i*revp) &
                          + arlx1i*p5*etax2T(i,j)*shearT(i,j)) * denom1

      enddo                     ! ij

      end subroutine stressCD_T

!=======================================================================
! Computes the strain rates and internal stress components for U points
!
! author: JF Lemieux, ECCC
! Nov 2021

      subroutine stressCD_U (nx_block,   ny_block, &
                                         icellU  , &
                             indxUi    , indxUj  , &
                             uarea     ,           &
                             zetax2U   , etax2U  , &
                             strengthU ,           &
                             divergU   , tensionU, &
                             shearU    , deltaU  , &
                             stresspU  , stressmU, &
                             stress12U)

      use ice_dyn_shared, only: visc_replpress, &
                                visc_method, deltaminEVP

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellU                ! no. of cells where iceUmask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxUi   , & ! compressed index in i-direction
         indxUj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         uarea    , & ! area of U-cell (m^2)
         zetax2U  , & ! 2*zeta at U point
         etax2U   , & ! 2*eta at U point
         strengthU, & ! ice strength at U point
         divergU  , & ! div at U point
         tensionU , & ! tension at U point
         shearU   , & ! shear at U point
         deltaU       ! delt at U point

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         stresspU , & ! sigma11+sigma22
         stressmU , & ! sigma11-sigma22
         stress12U    ! sigma12

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &
        lzetax2U  , & ! bulk viscosity at U point
        letax2U   , & ! shear viscosity at U point
        lrep_prsU , & ! replacement pressure at U point
        DminUarea     ! Dmin on U

      character(len=*), parameter :: subname = '(stressCD_U)'

      do ij = 1, icellU
         i = indxUi(ij)
         j = indxUj(ij)

         !-----------------------------------------------------------------
         ! viscosities and replacement pressure at U point
         ! avg_zeta: Bouillon et al. 2013, C1 method of Kimmritz et al. 2016
         ! avg_strength: C2 method of Kimmritz et al. 2016
         !-----------------------------------------------------------------

         if (visc_method == 'avg_zeta') then
            lzetax2U  = zetax2U(i,j)
            letax2U   = etax2U(i,j)
            lrep_prsU = (c1-Ktens)/(c1+Ktens)*lzetax2U*deltaU(i,j)

         elseif (visc_method == 'avg_strength') then
            DminUarea = deltaminEVP*uarea(i,j)
            ! only need etax2U here, but other terms are calculated with etax2U
            ! minimal extra calculations here even though it seems like there is
            call visc_replpress (strengthU(i,j), DminUarea, deltaU(i,j), &
                                 lzetax2U      , letax2U  , lrep_prsU )
         endif

         !-----------------------------------------------------------------
         ! the stresses                            ! kg/s^2
         !-----------------------------------------------------------------

         ! NOTE: for comp. efficiency 2 x zeta and 2 x eta are used in the code

         stresspU(i,j)  = (stresspU (i,j)*(c1-arlx1i*revp) &
                          + arlx1i*(lzetax2U*divergU(i,j) - lrep_prsU)) * denom1

         stressmU(i,j)  = (stressmU (i,j)*(c1-arlx1i*revp) &
                          + arlx1i*letax2U*tensionU(i,j)) * denom1

         stress12U(i,j) = (stress12U(i,j)*(c1-arlx1i*revp) &
                          + arlx1i*p5*letax2U*shearU(i,j)) * denom1

      enddo                     ! ij

      end subroutine stressCD_U

!=======================================================================
! Computes divergence of stress tensor at the E or N point for the mom equation
!
! author: JF Lemieux, ECCC
! Nov 2021
!
! Hunke, E. C., and J. K. Dukowicz (2002).  The Elastic-Viscous-Plastic
! Sea Ice Dynamics Model in General Orthogonal Curvilinear Coordinates
! on a Sphere - Incorporation of Metric Terms. Mon. Weather Rev.,
! 130, 1848-1865.
!
! Bouillon, S., M. Morales Maqueda, V. Legat and T. Fichefet (2009). An
! elastic-viscous-plastic sea ice model formulated on Arakawa B and C grids.
! Ocean Model., 27, 174-184.

      subroutine div_stress_Ex(nx_block, ny_block, &
                                         icell   , &
                               indxi   , indxj   , &
                               dxE     , dyE     , &
                               dxU     , dyT     , &
                               arear   ,           &
                               stressp , stressm , &
                               stress12,           &
                               strintx )


      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icell                 ! no. of cells where epm (or npm) = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxi   , & ! compressed index in i-direction
         indxj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         dxE     , & ! width of E or N-cell through the middle (m)
         dyE     , & ! height of E or N-cell through the middle (m)
         dxU     , & ! width of T or U-cell through the middle (m)
         dyT     , & ! height of T or U-cell through the middle (m)
         arear       ! earear or narear

      real (kind=dbl_kind), optional, dimension (nx_block,ny_block), intent(in) :: &
         stressp , & ! stressp  (U or T) used for strintx calculation
         stressm , & ! stressm  (U or T) used for strintx calculation
         stress12    ! stress12 (U or T) used for strintx calculation

      real (kind=dbl_kind), optional, dimension (nx_block,ny_block), intent(out) :: &
         strintx     ! div of stress tensor for u component

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      character(len=*), parameter :: subname = '(div_stress_Ex)'

      do ij = 1, icell
         i = indxi(ij)
         j = indxj(ij)
         strintx(i,j) = arear(i,j) * &
              ( p5 * dyE(i,j)  * ( stressp(i+1,j  )  - stressp (i  ,j  ) ) &
              + (p5/ dyE(i,j)) * ( (dyT(i+1,j  )**2) * stressm (i+1,j  )   &
                                  -(dyT(i  ,j  )**2) * stressm (i  ,j  ) ) &
              + (c1/ dxE(i,j)) * ( (dxU(i  ,j  )**2) * stress12(i  ,j  )   &
                                  -(dxU(i  ,j-1)**2) * stress12(i  ,j-1) ) )
      enddo

      end subroutine div_stress_Ex

!=======================================================================
      subroutine div_stress_Ey(nx_block, ny_block, &
                                         icell   , &
                               indxi   , indxj   , &
                               dxE     , dyE     , &
                               dxU     , dyT     , &
                               arear   ,           &
                               stressp , stressm , &
                               stress12,           &
                               strinty )

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icell                 ! no. of cells where epm (or npm) = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxi   , & ! compressed index in i-direction
         indxj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         dxE     , & ! width of E or N-cell through the middle (m)
         dyE     , & ! height of E or N-cell through the middle (m)
         dxU     , & ! width of T or U-cell through the middle (m)
         dyT     , & ! height of T or U-cell through the middle (m)
         arear         ! earear or narear

      real (kind=dbl_kind), optional, dimension (nx_block,ny_block), intent(in) :: &
         stressp , & ! stressp  (U or T) used for strinty calculation
         stressm , & ! stressm  (U or T) used for strinty calculation
         stress12    ! stress12 (U or T) used for strinty calculation

      real (kind=dbl_kind), optional, dimension (nx_block,ny_block), intent(out) :: &
         strinty     ! div of stress tensor for v component

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      character(len=*), parameter :: subname = '(div_stress_Ey)'

      do ij = 1, icell
         i = indxi(ij)
         j = indxj(ij)
         strinty(i,j) = arear(i,j) * &
              ( p5 * dxE(i,j)  * ( stressp(i  ,j  )  - stressp (i  ,j-1) ) &
              - (p5/ dxE(i,j)) * ( (dxU(i  ,j  )**2) * stressm (i  ,j  )   &
                                  -(dxU(i  ,j-1)**2) * stressm (i  ,j-1) ) &
              + (c1/ dyE(i,j)) * ( (dyT(i+1,j  )**2) * stress12(i+1,j  )   &
                                  -(dyT(i  ,j  )**2) * stress12(i  ,j  ) ) )
      enddo

      end subroutine div_stress_Ey

!=======================================================================
      subroutine div_stress_Nx(nx_block, ny_block, &
                                         icell   , &
                               indxi   , indxj   , &
                               dxN     , dyN     , &
                               dxT     , dyU     , &
                               arear   ,           &
                               stressp , stressm , &
                               stress12,           &
                               strintx )

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icell                 ! no. of cells where epm (or npm) = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxi   , & ! compressed index in i-direction
         indxj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         dxN     , & ! width of E or N-cell through the middle (m)
         dyN     , & ! height of E or N-cell through the middle (m)
         dxT     , & ! width of T or U-cell through the middle (m)
         dyU     , & ! height of T or U-cell through the middle (m)
         arear       ! earear or narear

      real (kind=dbl_kind), optional, dimension (nx_block,ny_block), intent(in) :: &
         stressp , & ! stressp  (U or T) used for strintx calculation
         stressm , & ! stressm  (U or T) used for strintx calculation
         stress12    ! stress12 (U or T) used for strintx calculation

      real (kind=dbl_kind), optional, dimension (nx_block,ny_block), intent(out) :: &
         strintx     ! div of stress tensor for u component

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      character(len=*), parameter :: subname = '(div_stress_Nx)'

      do ij = 1, icell
         i = indxi(ij)
         j = indxj(ij)
         strintx(i,j) = arear(i,j) * &
              ( p5 * dyN(i,j)  * ( stressp(i  ,j  )  - stressp (i-1,j  ) ) &
              + (p5/ dyN(i,j)) * ( (dyU(i  ,j  )**2) * stressm (i  ,j  )   &
                                  -(dyU(i-1,j  )**2) * stressm (i-1,j  ) ) &
              + (c1/ dxN(i,j)) * ( (dxT(i  ,j+1)**2) * stress12(i  ,j+1)   &
                                  -(dxT(i  ,j  )**2) * stress12(i  ,j  ) ) )
      enddo

      end subroutine div_stress_Nx

!=======================================================================
      subroutine div_stress_Ny(nx_block, ny_block, &
                                         icell   , &
                               indxi   , indxj   , &
                               dxN     , dyN     , &
                               dxT     , dyU     , &
                               arear   ,           &
                               stressp , stressm , &
                               stress12,           &
                               strinty )

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icell                 ! no. of cells where epm (or npm) = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxi   , & ! compressed index in i-direction
         indxj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         dxN     , & ! width of E or N-cell through the middle (m)
         dyN     , & ! height of E or N-cell through the middle (m)
         dxT     , & ! width of T or U-cell through the middle (m)
         dyU     , & ! height of T or U-cell through the middle (m)
         arear       ! earear or narear

      real (kind=dbl_kind), optional, dimension (nx_block,ny_block), intent(in) :: &
         stressp , & ! stressp  (U or T) used for strinty calculation
         stressm , & ! stressm  (U or T) used for strinty calculation
         stress12    ! stress12 (U or T) used for strinty calculation

      real (kind=dbl_kind), optional, dimension (nx_block,ny_block), intent(out) :: &
         strinty     ! div of stress tensor for v component

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      character(len=*), parameter :: subname = '(div_stress_Ny)'

      do ij = 1, icell
         i = indxi(ij)
         j = indxj(ij)
         strinty(i,j) = arear(i,j) * &
              ( p5 * dxN(i,j)  * ( stressp(i  ,j+1)  - stressp (i  ,j  ) ) &
              - (p5/ dxN(i,j)) * ( (dxT(i  ,j+1)**2) * stressm (i  ,j+1)   &
                                  -(dxT(i  ,j  )**2) * stressm (i  ,j  ) ) &
              + (c1/ dyN(i,j)) * ( (dyU(i  ,j  )**2) * stress12(i  ,j  )   &
                                  -(dyU(i-1,j  )**2) * stress12(i-1,j  ) ) )
      enddo

      end subroutine div_stress_Ny

!=======================================================================

      end module ice_dyn_evp

!=======================================================================
