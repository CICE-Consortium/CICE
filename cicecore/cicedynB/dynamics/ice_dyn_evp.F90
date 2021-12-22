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
      use ice_dyn_shared, only: stepu, step_vel, dyn_prep1, dyn_prep2, dyn_finish, &
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
      public :: evp

!=======================================================================

      contains

!=======================================================================

! Elastic-viscous-plastic dynamics driver
!
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
      use ice_domain, only: nblocks, blocks_ice, halo_info, maskhalo_dyn
      use ice_domain_size, only: max_blocks, ncat, nx_global, ny_global
      use ice_flux, only: rdg_conv, rdg_shear, strairxT, strairyT, &
          strairx, strairy, uocn, vocn, ss_tltx, ss_tlty, iceumask, fm, &
          strtltx, strtlty, strocnx, strocny, strintx, strinty, taubx, tauby, &
          strocnxT, strocnyT, strax, stray, &
          Tbu, hwater, &
          strairxN, strairyN, icenmask, fmN, &
          strtltxN, strtltyN, strocnxN, strocnyN, strintxN, strintyN, taubxN, taubyN, &
          TbN, &
          strairxE, strairyE, iceemask, fmE, &
          strtltxE, strtltyE, strocnxE, strocnyE, strintxE, strintyE, taubxE, taubyE, &
          TbE, &
          stressp_1, stressp_2, stressp_3, stressp_4, &
          stressm_1, stressm_2, stressm_3, stressm_4, &
          stress12_1, stress12_2, stress12_3, stress12_4, &
          stresspT, stressmT, stress12T, &
          stresspU, stressmU, stress12U
      use ice_grid, only: hm, tmask, umask, umaskCD, nmask, emask, uvm, epm, npm, &
          dxe, dxn, dxt, dxu, dye, dyn, dyt, dyu, &
          ratiodxN, ratiodxNr, ratiodyE, ratiodyEr, & 
          dxhy, dyhx, cxp, cyp, cxm, cym, &
          tarear, uarear, earear, narear, tinyarea, grid_average_X2Y, tarea, &
          grid_type, grid_ice, &
          grid_atm_dynu, grid_atm_dynv, grid_ocn_dynu, grid_ocn_dynv
      use ice_state, only: aice, vice, vsno, uvel, vvel, uvelN, vvelN, &
          uvelE, vvelE, divu, shear, &
          aice_init, aice0, aicen, vicen, strength
      use ice_timers, only: timer_dynamics, timer_bound, &
          ice_timer_start, ice_timer_stop, timer_evp_1d, timer_evp_2d
      use ice_dyn_evp_1d, only: ice_dyn_evp_1d_copyin, ice_dyn_evp_1d_kernel, &
          ice_dyn_evp_1d_copyout
      use ice_dyn_shared, only: evp_algorithm, stack_velocity_field, unstack_velocity_field

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) :: & 
         ksub           , & ! subcycle step
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, ij

      integer (kind=int_kind), dimension(max_blocks) :: & 
         icellt   , & ! no. of cells where icetmask = 1
         icelln   , & ! no. of cells where icenmask = 1
         icelle   , & ! no. of cells where iceemask = 1
         icellu       ! no. of cells where iceumask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block, max_blocks) :: &
         indxti   , & ! compressed index in i-direction
         indxtj   , & ! compressed index in j-direction
         indxei   , & ! compressed index in i-direction
         indxej   , & ! compressed index in j-direction
         indxni   , & ! compressed index in i-direction
         indxnj   , & ! compressed index in j-direction
         indxui   , & ! compressed index in i-direction
         indxuj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         uocnU    , & ! i ocean current (m/s)
         vocnU    , & ! j ocean current (m/s)
         ss_tltxU , & ! sea surface slope, x-direction (m/m)
         ss_tltyU , & ! sea surface slope, y-direction (m/m)
         tmass    , & ! total mass of ice and snow (kg/m^2)
         waterx   , & ! for ocean stress calculation, x (m/s)
         watery   , & ! for ocean stress calculation, y (m/s)
         forcex   , & ! work array: combined atm stress and ocn tilt, x
         forcey   , & ! work array: combined atm stress and ocn tilt, y
         aiu      , & ! ice fraction on u-grid
         umass    , & ! total mass of ice and snow (u grid)
         umassdti     ! mass of U-cell/dte (kg/m^2 s)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         uocnN    , & ! i ocean current (m/s)
         vocnN    , & ! j ocean current (m/s)
         ss_tltxN , & ! sea surface slope, x-direction (m/m)
         ss_tltyN , & ! sea surface slope, y-direction (m/m)
         waterxN  , & ! for ocean stress calculation, x (m/s)
         wateryN  , & ! for ocean stress calculation, y (m/s)
         forcexN  , & ! work array: combined atm stress and ocn tilt, x
         forceyN  , & ! work array: combined atm stress and ocn tilt, y
         aiN      , & ! ice fraction on N-grid
         nmass    , & ! total mass of ice and snow (N grid)
         nmassdti     ! mass of N-cell/dte (kg/m^2 s)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         uocnE    , & ! i ocean current (m/s)
         vocnE    , & ! j ocean current (m/s)
         ss_tltxE , & ! sea surface slope, x-direction (m/m)
         ss_tltyE , & ! sea surface slope, y-direction (m/m)
         waterxE  , & ! for ocean stress calculation, x (m/s)
         wateryE  , & ! for ocean stress calculation, y (m/s)
         forcexE  , & ! work array: combined atm stress and ocn tilt, x
         forceyE  , & ! work array: combined atm stress and ocn tilt, y
         aiE      , & ! ice fraction on E-grid
         emass    , & ! total mass of ice and snow (E grid)
         emassdti     ! mass of E-cell/dte (kg/m^2 s)

      real (kind=dbl_kind), allocatable :: fld2(:,:,:,:)

      real (kind=dbl_kind), allocatable :: &
         zetax2T(:,:,:), & ! zetax2 = 2*zeta (bulk viscous coeff)
         etax2T(:,:,:)     ! etax2  = 2*eta  (shear viscous coeff)
      
      real (kind=dbl_kind), dimension(nx_block,ny_block,8):: &
         strtmp       ! stress combinations for momentum equation

      logical (kind=log_kind) :: calc_strair

      integer (kind=int_kind), dimension (nx_block,ny_block,max_blocks) :: &
         icetmask, &  ! ice extent mask (T-cell)
         halomask     ! generic halo mask

      type (ice_halo) :: &
         halo_info_mask !  ghost cell update info for masked halo

      type (block) :: &
         this_block           ! block information for current block

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1, &      ! temporary
         work2         ! temporary

      logical (kind=log_kind), save :: first_time = .true.
      
      character(len=*), parameter :: subname = '(evp)'

      call ice_timer_start(timer_dynamics) ! dynamics

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      allocate(fld2(nx_block,ny_block,2,max_blocks))

      if (grid_ice == 'CD') then

         allocate(zetax2T(nx_block,ny_block,max_blocks))
         allocate(etax2T(nx_block,ny_block,max_blocks))
         zetax2T(:,:,:) = c0
         etax2T (:,:,:) = c0
         
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

      call grid_average_X2Y('F',tmass,'T',umass,'U')
      call grid_average_X2Y('F',aice_init,'T',aiu,'U')
      call grid_average_X2Y('S',uocn,grid_ocn_dynu,uocnU,'U')
      call grid_average_X2Y('S',vocn,grid_ocn_dynv,vocnU,'U')
      call grid_average_X2Y('S',ss_tltx,grid_ocn_dynu,ss_tltxU,'U')
      call grid_average_X2Y('S',ss_tlty,grid_ocn_dynv,ss_tltyU,'U')

      if (grid_ice == 'CD') then
         call grid_average_X2Y('F',tmass,'T',emass,'E')
         call grid_average_X2Y('F',aice_init,'T', aie,'E')
         call grid_average_X2Y('F',tmass,'T',nmass,'N')
         call grid_average_X2Y('F',aice_init,'T', ain,'N')
         call grid_average_X2Y('S',uocn,grid_ocn_dynu,uocnN,'N')
         call grid_average_X2Y('S',vocn,grid_ocn_dynv,vocnN,'N')
         call grid_average_X2Y('S',uocn,grid_ocn_dynu,uocnE,'E')
         call grid_average_X2Y('S',vocn,grid_ocn_dynv,vocnE,'E')
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
         call grid_average_X2Y('F', strax, grid_atm_dynu, strairx, 'U')
         call grid_average_X2Y('F', stray, grid_atm_dynv, strairy, 'U')
      else
         call ice_HaloUpdate (strairxT,         halo_info, &
                              field_loc_center, field_type_vector)
         call ice_HaloUpdate (strairyT,         halo_info, &
                              field_loc_center, field_type_vector)
         call grid_average_X2Y('F',strairxT,'T',strairx,'U')
         call grid_average_X2Y('F',strairyT,'T',strairy,'U')
      endif

      if (grid_ice == 'CD') then
         if (.not. calc_strair) then
            call grid_average_X2Y('F', strax, grid_atm_dynu, strairxN, 'N')
            call grid_average_X2Y('F', stray, grid_atm_dynv, strairyN, 'N')
            call grid_average_X2Y('F', strax, grid_atm_dynu, strairxE, 'E')
            call grid_average_X2Y('F', stray, grid_atm_dynv, strairyE, 'E')
         else
            call grid_average_X2Y('F',strairxT,'T',strairxN,'N')
            call grid_average_X2Y('F',strairyT,'T',strairyN,'N')
            call grid_average_X2Y('F',strairxT,'T',strairxE,'E')
            call grid_average_X2Y('F',strairyT,'T',strairyE,'E')
         endif      
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
         
         if (trim(grid_ice) == 'B') then
            call dyn_prep2 (nx_block,             ny_block,             & 
                            ilo, ihi,             jlo, jhi,             &
                            icellt(iblk),         icellu(iblk),         & 
                            indxti      (:,iblk), indxtj      (:,iblk), & 
                            indxui      (:,iblk), indxuj      (:,iblk), & 
                            aiu       (:,:,iblk), umass     (:,:,iblk), & 
                            umassdti  (:,:,iblk), fcor_blk  (:,:,iblk), & 
                            umask     (:,:,iblk),                       & 
                            uocnU     (:,:,iblk), vocnU     (:,:,iblk), &
                            strairx   (:,:,iblk), strairy   (:,:,iblk), & 
                            ss_tltxU  (:,:,iblk), ss_tltyU  (:,:,iblk), &  
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

         elseif (trim(grid_ice) == 'CD') then
            call dyn_prep2 (nx_block,             ny_block,             &
                            ilo, ihi,             jlo, jhi,             &
                            icellt(iblk),         icellu(iblk),         &
                            indxti      (:,iblk), indxtj      (:,iblk), &
                            indxui      (:,iblk), indxuj      (:,iblk), &
                            aiu       (:,:,iblk), umass     (:,:,iblk), &
                            umassdti  (:,:,iblk), fcor_blk  (:,:,iblk), &
                            umaskCD   (:,:,iblk),                       &
                            uocnU     (:,:,iblk), vocnU     (:,:,iblk), &
                            strairx   (:,:,iblk), strairy   (:,:,iblk), &
                            ss_tltxU  (:,:,iblk), ss_tltyU  (:,:,iblk), &
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
         endif

      !-----------------------------------------------------------------
      ! ice strength
      !-----------------------------------------------------------------

         strength(:,:,iblk) = c0  ! initialize
         do ij = 1, icellt(iblk)
            i = indxti(ij, iblk)
            j = indxtj(ij, iblk)
            call icepack_ice_strength(ncat     = ncat,                 &
                                      aice     = aice    (i,j,  iblk), & 
                                      vice     = vice    (i,j,  iblk), & 
                                      aice0    = aice0   (i,j,  iblk), & 
                                      aicen    = aicen   (i,j,:,iblk), &  
                                      vicen    = vicen   (i,j,:,iblk), & 
                                      strength = strength(i,j,  iblk) )
         enddo  ! ij

      enddo  ! iblk
      !$TCXOMP END PARALLEL DO

      if (grid_ice == 'CD') then

      !$TCXOMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks

      !-----------------------------------------------------------------
      ! more preparation for dynamics on N grid
      !-----------------------------------------------------------------

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         call dyn_prep2 (nx_block,             ny_block,             & 
                         ilo, ihi,             jlo, jhi,             &
                         icellt(iblk),         icelln(iblk),         & 
                         indxti      (:,iblk), indxtj      (:,iblk), & 
                         indxni      (:,iblk), indxnj      (:,iblk), & 
                         aiN       (:,:,iblk), nmass     (:,:,iblk), & 
                         nmassdti  (:,:,iblk), fcorN_blk  (:,:,iblk),& 
                         nmask     (:,:,iblk),                       &
                         uocnN     (:,:,iblk), vocnN     (:,:,iblk), & 
                         strairxN  (:,:,iblk), strairyN  (:,:,iblk), & 
                         ss_tltxN  (:,:,iblk), ss_tltyN  (:,:,iblk), &  
                         icetmask  (:,:,iblk), icenmask  (:,:,iblk), & 
                         fmN       (:,:,iblk), dt,                   & 
                         strtltxN  (:,:,iblk), strtltyN  (:,:,iblk), & 
                         strocnxN  (:,:,iblk), strocnyN  (:,:,iblk), & 
                         strintxN  (:,:,iblk), strintyN  (:,:,iblk), & 
                         taubxN    (:,:,iblk), taubyN    (:,:,iblk), & 
                         waterxN   (:,:,iblk), wateryN   (:,:,iblk), & 
                         forcexN   (:,:,iblk), forceyN   (:,:,iblk), & 
                         stresspT  (:,:,iblk), stressp_2 (:,:,iblk), & 
                         stressp_3 (:,:,iblk), stressp_4 (:,:,iblk), & 
                         stressmT  (:,:,iblk), stressm_2 (:,:,iblk), & 
                         stressm_3 (:,:,iblk), stressm_4 (:,:,iblk), & 
                         stress12T (:,:,iblk), stress12_2(:,:,iblk), & 
                         stress12_3(:,:,iblk), stress12_4(:,:,iblk), & 
                         uvelN_init (:,:,iblk), vvelN_init (:,:,iblk), &
                         uvelN      (:,:,iblk), vvelN      (:,:,iblk), &
                         TbN       (:,:,iblk))

      !-----------------------------------------------------------------
      ! more preparation for dynamics on E grid
      !-----------------------------------------------------------------

         call dyn_prep2 (nx_block,             ny_block,             & 
                         ilo, ihi,             jlo, jhi,             &
                         icellt(iblk),         icelle(iblk),         & 
                         indxti      (:,iblk), indxtj      (:,iblk), & 
                         indxei      (:,iblk), indxej      (:,iblk), & 
                         aiE       (:,:,iblk), emass     (:,:,iblk), & 
                         emassdti  (:,:,iblk), fcorE_blk  (:,:,iblk),& 
                         emask     (:,:,iblk),                       & 
                         uocnE     (:,:,iblk), vocnE     (:,:,iblk), & 
                         strairxE  (:,:,iblk), strairyE  (:,:,iblk), & 
                         ss_tltxE  (:,:,iblk), ss_tltyE  (:,:,iblk), &  
                         icetmask  (:,:,iblk), iceemask  (:,:,iblk), & 
                         fmE       (:,:,iblk), dt,                   & 
                         strtltxE  (:,:,iblk), strtltyE  (:,:,iblk), & 
                         strocnxE  (:,:,iblk), strocnyE  (:,:,iblk), & 
                         strintxE  (:,:,iblk), strintyE  (:,:,iblk), & 
                         taubxE    (:,:,iblk), taubyE    (:,:,iblk), & 
                         waterxE   (:,:,iblk), wateryE   (:,:,iblk), & 
                         forcexE   (:,:,iblk), forceyE   (:,:,iblk), & 
                         stresspU  (:,:,iblk), stressp_2 (:,:,iblk), & 
                         stressp_3 (:,:,iblk), stressp_4 (:,:,iblk), & 
                         stressmU  (:,:,iblk), stressm_2 (:,:,iblk), & 
                         stressm_3 (:,:,iblk), stressm_4 (:,:,iblk), & 
                         stress12U (:,:,iblk), stress12_2(:,:,iblk), & 
                         stress12_3(:,:,iblk), stress12_4(:,:,iblk), & 
                         uvelE_init (:,:,iblk), vvelE_init (:,:,iblk), &
                         uvelE      (:,:,iblk), vvelE      (:,:,iblk), &
                         TbE       (:,:,iblk))

      enddo  ! iblk
      !$TCXOMP END PARALLEL DO
      
      endif ! grid_ice

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (grid_ice == 'CD') then

         call ice_timer_start(timer_bound)
         ! velocities may have changed in dyn_prep2
         call stack_velocity_field(uvelN, vvelN, fld2)
         call ice_HaloUpdate (fld2,               halo_info, &
                              field_loc_Nface, field_type_vector)
         call unstack_velocity_field(fld2, uvelN, vvelN)
         ! velocities may have changed in dyn_prep2
         call stack_velocity_field(uvelE, vvelE, fld2)
         call ice_HaloUpdate (fld2,               halo_info, &
                              field_loc_Eface, field_type_vector)
         call unstack_velocity_field(fld2, uvelE, vvelE)
         call ice_timer_stop(timer_bound)

         call grid_average_X2Y('S',uvelE,'E',uvel,'U')
         call grid_average_X2Y('S',vvelN,'N',vvel,'U')
      endif

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
      ! seabed stress factor Tbu (Tbu is part of Cb coefficient)
      !-----------------------------------------------------------------
      
      if (seabed_stress) then

         ! tcraig, causes abort with pgi compiler on cheyenne
         !$TCXOMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks

            select case (trim(grid_ice))
            case('B')

               if ( seabed_stress_method == 'LKD' ) then

                  call seabed_stress_factor_LKD (nx_block,         ny_block,       &
                                                 icellu  (iblk),                   &
                                                 indxui(:,iblk),   indxuj(:,iblk), &
                                                 vice(:,:,iblk),   aice(:,:,iblk), &
                                                 hwater(:,:,iblk), Tbu(:,:,iblk))

               elseif ( seabed_stress_method == 'probabilistic' ) then

                  call seabed_stress_factor_prob (nx_block,         ny_block,                   &
                                                  icellt(iblk), indxti(:,iblk), indxtj(:,iblk), &
                                                  icellu(iblk), indxui(:,iblk), indxuj(:,iblk), &
                                                  aicen(:,:,:,iblk), vicen(:,:,:,iblk),         &
                                                  hwater(:,:,iblk), Tbu(:,:,iblk))
               endif

         case('CD')

            if ( seabed_stress_method == 'LKD' ) then

               call seabed_stress_factor_LKD (nx_block,         ny_block,       &
                                              icelle  (iblk),                   &
                                              indxei(:,iblk),   indxej(:,iblk), &
                                              vice(:,:,iblk),   aice(:,:,iblk), &
                                              hwater(:,:,iblk), TbE(:,:,iblk))
               call seabed_stress_factor_LKD (nx_block,         ny_block,       &
                                              icelln  (iblk),                   &
                                              indxni(:,iblk),   indxnj(:,iblk), &
                                              vice(:,:,iblk),   aice(:,:,iblk), &
                                              hwater(:,:,iblk), TbN(:,:,iblk))

            elseif ( seabed_stress_method == 'probabilistic' ) then

               call seabed_stress_factor_prob (nx_block,         ny_block,                   &
                                               icellt(iblk), indxti(:,iblk), indxtj(:,iblk), &
                                               icellu(iblk), indxui(:,iblk), indxuj(:,iblk), &
                                               aicen(:,:,:,iblk), vicen(:,:,:,iblk),         &
                                               hwater(:,:,iblk), Tbu(:,:,iblk),              &
                                               TbE(:,:,iblk),    TbN(:,:,iblk),              &
                                               icelle(iblk), indxei(:,iblk), indxej(:,iblk), &
                                               icelln(iblk), indxni(:,iblk), indxnj(:,iblk)  )
            endif

         end select
         
         enddo
       !$TCXOMP END PARALLEL DO
      endif

      call ice_timer_start(timer_evp_2d)

      if (evp_algorithm == "shared_mem_1d" ) then

         if (first_time .and. my_task == master_task) then
            write(nu_diag,'(3a)') subname,' Entering evp_algorithm version ',evp_algorithm
            first_time = .false.
         endif
         if (trim(grid_type) == 'tripole') then
            call abort_ice(trim(subname)//' &
               & Kernel not tested on tripole grid. Set evp_algorithm=standard_2d')
         endif

         call ice_dyn_evp_1d_copyin(                                                &
            nx_block,ny_block,nblocks,nx_global+2*nghost,ny_global+2*nghost, &
            icetmask, iceumask,                                           &
            cdn_ocn,aiu,uocnU,vocnU,forcex,forcey,Tbu,        &
            umassdti,fm,uarear,tarear,strintx,strinty,uvel_init,vvel_init,&
            strength,uvel,vvel,dxt,dyt,                                   &
            stressp_1 ,stressp_2, stressp_3, stressp_4,                   &
            stressm_1 ,stressm_2, stressm_3, stressm_4,                   &
            stress12_1,stress12_2,stress12_3,stress12_4                   )
         call ice_timer_start(timer_evp_1d)
         call ice_dyn_evp_1d_kernel()
         call ice_timer_stop(timer_evp_1d)
         call ice_dyn_evp_1d_copyout(                                      &
            nx_block,ny_block,nblocks,nx_global+2*nghost,ny_global+2*nghost,&
!strocn            uvel,vvel, strocnx,strocny, strintx,strinty,                  &
            uvel,vvel, strintx,strinty,                                   &
            stressp_1, stressp_2, stressp_3, stressp_4,                   &
            stressm_1, stressm_2, stressm_3, stressm_4,                   &
            stress12_1,stress12_2,stress12_3,stress12_4,                  &
            divu,rdg_conv,rdg_shear,shear,taubx,tauby                     )

      else ! evp_algorithm == standard_2d (Standard CICE)

         do ksub = 1,ndte        ! subcycling

         !-----------------------------------------------------------------
         ! stress tensor equation, total surface stress
         !-----------------------------------------------------------------

            !$TCXOMP PARALLEL DO PRIVATE(iblk,strtmp)
            do iblk = 1, nblocks

               select case (grid_ice)
               case('B')
                  call stress (nx_block,             ny_block,             &
                               ksub,                 icellt(iblk),         &
                               indxti      (:,iblk), indxtj      (:,iblk), &
                               uvel      (:,:,iblk), vvel      (:,:,iblk), &
                               dxt       (:,:,iblk), dyt       (:,:,iblk), &
                               dxhy      (:,:,iblk), dyhx      (:,:,iblk), &
                               cxp       (:,:,iblk), cyp       (:,:,iblk), &
                               cxm       (:,:,iblk), cym       (:,:,iblk), &
                               tarear    (:,:,iblk), tinyarea  (:,:,iblk), &
                               strength  (:,:,iblk),                       &
                               stressp_1 (:,:,iblk), stressp_2 (:,:,iblk), &
                               stressp_3 (:,:,iblk), stressp_4 (:,:,iblk), &
                               stressm_1 (:,:,iblk), stressm_2 (:,:,iblk), &
                               stressm_3 (:,:,iblk), stressm_4 (:,:,iblk), &
                               stress12_1(:,:,iblk), stress12_2(:,:,iblk), &
                               stress12_3(:,:,iblk), stress12_4(:,:,iblk), &
                               shear     (:,:,iblk), divu      (:,:,iblk), &
                               rdg_conv  (:,:,iblk), rdg_shear (:,:,iblk), &
                               strtmp    (:,:,:) )

            !-----------------------------------------------------------------
            ! momentum equation
            !-----------------------------------------------------------------
                  call stepu (nx_block,            ny_block,           &
                              icellu       (iblk), Cdn_ocn (:,:,iblk), &
                              indxui     (:,iblk), indxuj    (:,iblk), &
                              ksub,                                    &
                              aiu      (:,:,iblk), strtmp  (:,:,:),    &
                              uocnU    (:,:,iblk), vocnU   (:,:,iblk), &
                              waterx   (:,:,iblk), watery  (:,:,iblk), &
                              forcex   (:,:,iblk), forcey  (:,:,iblk), &
                              umassdti (:,:,iblk), fm      (:,:,iblk), &
                              uarear   (:,:,iblk),                     &
                              strintx  (:,:,iblk), strinty (:,:,iblk), &
                              taubx    (:,:,iblk), tauby   (:,:,iblk), &
                              uvel_init(:,:,iblk), vvel_init(:,:,iblk),&
                              uvel     (:,:,iblk), vvel    (:,:,iblk), &
                              Tbu      (:,:,iblk))

               case('CD')

                  call stress_T (nx_block,             ny_block,             &
                                 ksub,                 icellt(iblk),         &
                                 indxti      (:,iblk), indxtj      (:,iblk), &
                                 uvelE     (:,:,iblk), vvelE     (:,:,iblk), &
                                 uvelN     (:,:,iblk), vvelN     (:,:,iblk), &
                                 dxN       (:,:,iblk), dyE       (:,:,iblk), &
                                 dxT       (:,:,iblk), dyT       (:,:,iblk), &
                                 tarear    (:,:,iblk), tinyarea  (:,:,iblk), &
                                 strength  (:,:,iblk),                       &
                                 zetax2T   (:,:,iblk), etax2T    (:,:,iblk), &
                                 stresspT  (:,:,iblk), stressmT  (:,:,iblk), &
                                 stress12T (:,:,iblk),                       &
                                 shear     (:,:,iblk), divu      (:,:,iblk), &
                                 rdg_conv  (:,:,iblk), rdg_shear (:,:,iblk)  )

                  call stress_U (nx_block,             ny_block,             &
                                 ksub,                 icellu(iblk),         &
                                 indxui      (:,iblk), indxuj      (:,iblk), &
                                 uvelE     (:,:,iblk), vvelE     (:,:,iblk), &
                                 uvelN     (:,:,iblk), vvelN     (:,:,iblk), &
                                 uvel      (:,:,iblk), vvel      (:,:,iblk), &
                                 dxE       (:,:,iblk), dyN       (:,:,iblk), &
                                 dxU       (:,:,iblk), dyU       (:,:,iblk), &
                                 tarea     (:,:,iblk),                       &
                                 ratiodxN  (:,:,iblk), ratiodxNr (:,:,iblk), &
                                 ratiodyE  (:,:,iblk), ratiodyEr (:,:,iblk), &
                                 epm       (:,:,iblk), npm       (:,:,iblk), &
                                 hm        (:,:,iblk), uvm       (:,:,iblk), &
                                 zetax2T   (:,:,iblk), etax2T    (:,:,iblk), &
                                 stresspU  (:,:,iblk), stressmU  (:,:,iblk), &
                                 stress12U (:,:,iblk))                   

                  call div_stress (nx_block,             ny_block,             & ! E point
                                   ksub,                 icelle(iblk),         &
                                   indxei      (:,iblk), indxej      (:,iblk), &
                                   dxE       (:,:,iblk), dyE       (:,:,iblk), &
                                   dxU       (:,:,iblk), dyT       (:,:,iblk), &
                                   earear    (:,:,iblk),                       &
                                   stresspT  (:,:,iblk), stressmT  (:,:,iblk), &
                                   stress12U (:,:,iblk),                       &
                                   stresspU  (:,:,iblk), stressmU  (:,:,iblk), &
                                   stress12T (:,:,iblk),                       &
                                   strintxE  (:,:,iblk), strintyE  (:,:,iblk), &
                                   'E')

                   call div_stress (nx_block,             ny_block,            & ! N point
                                   ksub,                 icelln(iblk),         &
                                   indxni      (:,iblk), indxnj      (:,iblk), &
                                   dxN       (:,:,iblk), dyN       (:,:,iblk), &
                                   dxT       (:,:,iblk), dyU       (:,:,iblk), &
                                   narear    (:,:,iblk),                       &
                                   stresspU  (:,:,iblk), stressmU  (:,:,iblk), &
                                   stress12T (:,:,iblk),                       &
                                   stresspT  (:,:,iblk), stressmT  (:,:,iblk), &
                                   stress12U (:,:,iblk),                       &
                                   strintxN  (:,:,iblk), strintyN  (:,:,iblk), &
                                   'N')
                  
                  call step_vel (nx_block,             ny_block,             & ! E point
                                 icelle        (iblk), Cdn_ocn   (:,:,iblk), &
                                 indxei      (:,iblk), indxej      (:,iblk), &
                                 ksub,                 aiE       (:,:,iblk), &
                                 uocnE     (:,:,iblk), vocnE     (:,:,iblk), &
                                 waterxE   (:,:,iblk), wateryE   (:,:,iblk), &
                                 forcexE   (:,:,iblk), forceyE   (:,:,iblk), &
                                 emassdti  (:,:,iblk), fmE       (:,:,iblk), &
                                 strintxE  (:,:,iblk), strintyE  (:,:,iblk), &
                                 taubxE    (:,:,iblk), taubyE    (:,:,iblk), &
                                 uvelE_init(:,:,iblk), vvelE_init(:,:,iblk), &
                                 uvelE     (:,:,iblk), vvelE     (:,:,iblk), &
                                 TbE       (:,:,iblk))

                  call step_vel (nx_block,             ny_block,             & ! N point
                                 icelln        (iblk), Cdn_ocn   (:,:,iblk), &
                                 indxni      (:,iblk), indxnj      (:,iblk), &
                                 ksub,                 aiN       (:,:,iblk), &
                                 uocnN     (:,:,iblk), vocnN     (:,:,iblk), &
                                 waterxN   (:,:,iblk), wateryN   (:,:,iblk), &
                                 forcexN   (:,:,iblk), forceyN   (:,:,iblk), &
                                 nmassdti  (:,:,iblk), fmN       (:,:,iblk), &
                                 strintxN  (:,:,iblk), strintyN  (:,:,iblk), &
                                 taubxN    (:,:,iblk), taubyN    (:,:,iblk), &
                                 uvelN_init(:,:,iblk), vvelN_init(:,:,iblk), &
                                 uvelN     (:,:,iblk), vvelN     (:,:,iblk), &
                                 TbN       (:,:,iblk))

               end select

            enddo
            !$TCXOMP END PARALLEL DO

            if (grid_ice == 'CD') then

               call ice_timer_start(timer_bound)
               call stack_velocity_field(uvelN, vvelN, fld2)
               call ice_HaloUpdate (fld2,               halo_info, &
                                    field_loc_Nface, field_type_vector)
               call unstack_velocity_field(fld2, uvelN, vvelN)
               call stack_velocity_field(uvelE, vvelE, fld2)
               call ice_HaloUpdate (fld2,               halo_info, &
                                    field_loc_Eface, field_type_vector)
               call unstack_velocity_field(fld2, uvelE, vvelE)
               call ice_timer_stop(timer_bound)

               call grid_average_X2Y('S',uvelE,'E',uvel,'U')
               call grid_average_X2Y('S',vvelN,'N',vvel,'U')

            endif

            call ice_timer_start(timer_bound)
            call stack_velocity_field(uvel, vvel, fld2)
            if (maskhalo_dyn) then
               call ice_HaloUpdate (fld2,               halo_info_mask, &
                                    field_loc_NEcorner, field_type_vector)
            else
               call ice_HaloUpdate (fld2,               halo_info, &
                                    field_loc_NEcorner, field_type_vector)
            endif
            call unstack_velocity_field(fld2, uvel, vvel)
            call ice_timer_stop(timer_bound)
         
         enddo                     ! subcycling
      endif  ! evp_algorithm

      call ice_timer_stop(timer_evp_2d)

      deallocate(fld2)
      if (grid_ice == 'CD') then
         deallocate(zetax2T, etax2T)
      endif
      
      if (maskhalo_dyn) call ice_HaloDestroy(halo_info_mask)

      ! Force symmetry across the tripole seam
      if (trim(grid_type) == 'tripole') then
         ! TODO: CD-grid
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
               uocnU   (:,:,iblk), vocnU   (:,:,iblk), &
               aiu     (:,:,iblk), fm      (:,:,iblk), & 
               strintx (:,:,iblk), strinty (:,:,iblk), &
               strairx (:,:,iblk), strairy (:,:,iblk), &
               strocnx (:,:,iblk), strocny (:,:,iblk))

      enddo
      !$OMP END PARALLEL DO

      if (grid_ice == 'CD') then

         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks
   
            call dyn_finish                               & 
                 (nx_block,           ny_block,           & 
                  icelln      (iblk), Cdn_ocn (:,:,iblk), & 
                  indxni    (:,iblk), indxnj    (:,iblk), & 
                  uvelN   (:,:,iblk), vvelN   (:,:,iblk), & 
                  uocnN   (:,:,iblk), vocnN   (:,:,iblk), & 
                  aiN     (:,:,iblk), fmN     (:,:,iblk), & 
                  strintxN(:,:,iblk), strintyN(:,:,iblk), &
                  strairxN(:,:,iblk), strairyN(:,:,iblk), &
                  strocnxN(:,:,iblk), strocnyN(:,:,iblk))

            call dyn_finish                               & 
                 (nx_block,           ny_block,           & 
                  icelle      (iblk), Cdn_ocn (:,:,iblk), & 
                  indxei    (:,iblk), indxej    (:,iblk), & 
                  uvelE   (:,:,iblk), vvelE   (:,:,iblk), & 
                  uocnE   (:,:,iblk), vocnE   (:,:,iblk), & 
                  aiE     (:,:,iblk), fmE     (:,:,iblk), & 
                  strintxE(:,:,iblk), strintyE(:,:,iblk), &
                  strairxE(:,:,iblk), strairyE(:,:,iblk), &
                  strocnxE(:,:,iblk), strocnyE(:,:,iblk))

         enddo
         !$OMP END PARALLEL DO

      endif

      ! strocn computed on U, N, E as needed. Map strocn U divided by aiu to T
      ! TODO: This should be done elsewhere as part of generalization?
      ! TODO: Rename strocn[x,y]T since it's different than strocn[x,y][U,N,E]
      ! conservation requires aiu be divided before averaging
      work1 = c0
      work2 = c0
      !$OMP PARALLEL DO PRIVATE(iblk,ij,i,j)
      do iblk = 1, nblocks
      do ij = 1, icellu(iblk)
         i = indxui(ij,iblk)
         j = indxuj(ij,iblk)
         work1(i,j,iblk) = strocnx(i,j,iblk)/aiu(i,j,iblk)
         work2(i,j,iblk) = strocny(i,j,iblk)/aiu(i,j,iblk)
      enddo
      enddo
      !$OMP END PARALLEL DO
      call ice_HaloUpdate (work1,              halo_info, &
                           field_loc_NEcorner, field_type_vector)
      call ice_HaloUpdate (work2,              halo_info, &
                           field_loc_NEcorner, field_type_vector)
      call grid_average_X2Y('F',work1,'U',strocnxT,'T')    ! shift
      call grid_average_X2Y('F',work2,'U',strocnyT,'T')

      call ice_timer_stop(timer_dynamics)    ! dynamics

      end subroutine evp

!=======================================================================

! Computes the rates of strain and internal stress components for
! each of the four corners on each T-grid cell.
! Computes stress terms for the momentum equation
!
! author: Elizabeth C. Hunke, LANL

      subroutine stress (nx_block,   ny_block,   & 
                         ksub,       icellt,     & 
                         indxti,     indxtj,     & 
                         uvel,       vvel,       & 
                         dxt,        dyt,        & 
                         dxhy,       dyhx,       & 
                         cxp,        cyp,        & 
                         cxm,        cym,        & 
                         tarear,     tinyarea,   & 
                         strength,               & 
                         stressp_1,  stressp_2,  & 
                         stressp_3,  stressp_4,  & 
                         stressm_1,  stressm_2,  & 
                         stressm_3,  stressm_4,  & 
                         stress12_1, stress12_2, & 
                         stress12_3, stress12_4, & 
                         shear,      divu,       & 
                         rdg_conv,   rdg_shear,  & 
                         str )

      use ice_dyn_shared, only: strain_rates, deformations, viscous_coeffs_and_rep_pressure_T
        
      integer (kind=int_kind), intent(in) :: & 
         nx_block, ny_block, & ! block dimensions
         ksub              , & ! subcycling step
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
         dxhy     , & ! 0.5*(HTE - HTW)
         dyhx     , & ! 0.5*(HTN - HTS)
         cyp      , & ! 1.5*HTE - 0.5*HTW
         cxp      , & ! 1.5*HTN - 0.5*HTS
         cym      , & ! 0.5*HTE - 1.5*HTW
         cxm      , & ! 0.5*HTN - 1.5*HTS
         tarear   , & ! 1/tarea
         tinyarea     ! puny*tarea

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         stressp_1, stressp_2, stressp_3, stressp_4 , & ! sigma11+sigma22
         stressm_1, stressm_2, stressm_3, stressm_4 , & ! sigma11-sigma22
         stress12_1,stress12_2,stress12_3,stress12_4    ! sigma12

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         shear    , & ! strain rate II component (1/s)
         divu     , & ! strain rate I component, velocity divergence (1/s)
         rdg_conv , & ! convergence term for ridging (1/s)
         rdg_shear    ! shear term for ridging (1/s)

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
        zetax2ne, zetax2nw, zetax2se, zetax2sw    , & ! 2 x zeta (visc coeff) 
        etax2ne, etax2nw, etax2se, etax2sw        , & ! 2 x eta (visc coeff)
        rep_prsne, rep_prsnw, rep_prsse, rep_prssw, & ! replacement pressure
!       puny                                      , & ! puny
        ssigpn, ssigps, ssigpe, ssigpw            , &
        ssigmn, ssigms, ssigme, ssigmw            , &
        ssig12n, ssig12s, ssig12e, ssig12w        , &
        ssigp1, ssigp2, ssigm1, ssigm2, ssig121, ssig122, &
        csigpne, csigpnw, csigpse, csigpsw        , &
        csigmne, csigmnw, csigmse, csigmsw        , &
        csig12ne, csig12nw, csig12se, csig12sw    , &
        str12ew, str12we, str12ns, str12sn        , &
        strp_tmp, strm_tmp, tmp

      real(kind=dbl_kind),parameter :: capping = c1 ! of the viscous coef

      character(len=*), parameter :: subname = '(stress)'

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
      ! viscous coefficients and replacement pressure
      !-----------------------------------------------------------------

         call viscous_coeffs_and_rep_pressure_T (strength(i,j), tinyarea(i,j),&
                                                 Deltane,       zetax2ne,     &
                                                 etax2ne,       rep_prsne,    &
                                                 capping)
 
         call viscous_coeffs_and_rep_pressure_T (strength(i,j), tinyarea(i,j),&
                                                 Deltanw,       zetax2nw,     &
                                                 etax2nw,       rep_prsnw,    &
                                                 capping)

         call viscous_coeffs_and_rep_pressure_T (strength(i,j), tinyarea(i,j),&
                                                 Deltasw,       zetax2sw,     &
                                                 etax2sw,       rep_prssw,    &
                                                 capping)

         call viscous_coeffs_and_rep_pressure_T (strength(i,j), tinyarea(i,j),&
                                                 Deltase,       zetax2se,     &
                                                 etax2se,       rep_prsse,    &
                                                 capping)

         
      !-----------------------------------------------------------------
      ! the stresses                            ! kg/s^2
      ! (1) northeast, (2) northwest, (3) southwest, (4) southeast
      !-----------------------------------------------------------------

      ! NOTE: for comp. efficiency 2 x zeta and 2 x eta are used in the code 
         
         stressp_1(i,j) = (stressp_1(i,j)*(c1-arlx1i*revp) + &
                           arlx1i*(zetax2ne*divune - rep_prsne)) * denom1
         stressp_2(i,j) = (stressp_2(i,j)*(c1-arlx1i*revp) + &
                           arlx1i*(zetax2nw*divunw - rep_prsnw)) * denom1
         stressp_3(i,j) = (stressp_3(i,j)*(c1-arlx1i*revp) + &
                           arlx1i*(zetax2sw*divusw - rep_prssw)) * denom1
         stressp_4(i,j) = (stressp_4(i,j)*(c1-arlx1i*revp) + &
                           arlx1i*(zetax2se*divuse - rep_prsse)) * denom1

         stressm_1(i,j) = (stressm_1(i,j)*(c1-arlx1i*revp) + &
                           arlx1i*etax2ne*tensionne) * denom1
         stressm_2(i,j) = (stressm_2(i,j)*(c1-arlx1i*revp) + &
                           arlx1i*etax2nw*tensionnw) * denom1
         stressm_3(i,j) = (stressm_3(i,j)*(c1-arlx1i*revp) + &
                           arlx1i*etax2sw*tensionsw) * denom1
         stressm_4(i,j) = (stressm_4(i,j)*(c1-arlx1i*revp) + &
                           arlx1i*etax2se*tensionse) * denom1

         stress12_1(i,j) = (stress12_1(i,j)*(c1-arlx1i*revp) + &
                            arlx1i*p5*etax2ne*shearne) * denom1
         stress12_2(i,j) = (stress12_2(i,j)*(c1-arlx1i*revp) + &
                            arlx1i*p5*etax2nw*shearnw) * denom1
         stress12_3(i,j) = (stress12_3(i,j)*(c1-arlx1i*revp) + &
                            arlx1i*p5*etax2sw*shearsw) * denom1
         stress12_4(i,j) = (stress12_4(i,j)*(c1-arlx1i*revp) + &
                            arlx1i*p5*etax2se*shearse) * denom1

      !-----------------------------------------------------------------
      ! Eliminate underflows.
      ! The following code is commented out because it is relatively 
      ! expensive and most compilers include a flag that accomplishes
      ! the same thing more efficiently.  This code is cheaper than
      ! handling underflows if the compiler lacks a flag; uncomment
      ! it in that case.  The compiler flag is often described with the 
      ! phrase "flush to zero".
      !-----------------------------------------------------------------

!      call icepack_query_parameters(puny_out=puny)
!      call icepack_warnings_flush(nu_diag)
!      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
!         file=__FILE__, line=__LINE__)

!      stressp_1(i,j) = sign(max(abs(stressp_1(i,j)),puny),stressp_1(i,j))
!      stressp_2(i,j) = sign(max(abs(stressp_2(i,j)),puny),stressp_2(i,j))
!      stressp_3(i,j) = sign(max(abs(stressp_3(i,j)),puny),stressp_3(i,j))
!      stressp_4(i,j) = sign(max(abs(stressp_4(i,j)),puny),stressp_4(i,j))

!      stressm_1(i,j) = sign(max(abs(stressm_1(i,j)),puny),stressm_1(i,j))
!      stressm_2(i,j) = sign(max(abs(stressm_2(i,j)),puny),stressm_2(i,j))
!      stressm_3(i,j) = sign(max(abs(stressm_3(i,j)),puny),stressm_3(i,j))
!      stressm_4(i,j) = sign(max(abs(stressm_4(i,j)),puny),stressm_4(i,j))

!      stress12_1(i,j) = sign(max(abs(stress12_1(i,j)),puny),stress12_1(i,j))
!      stress12_2(i,j) = sign(max(abs(stress12_2(i,j)),puny),stress12_2(i,j))
!      stress12_3(i,j) = sign(max(abs(stress12_3(i,j)),puny),stress12_3(i,j))
!      stress12_4(i,j) = sign(max(abs(stress12_4(i,j)),puny),stress12_4(i,j))

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

      !-----------------------------------------------------------------
      ! on last subcycle, save quantities for mechanical redistribution
      !-----------------------------------------------------------------
      if (ksub == ndte) then
         call deformations (nx_block  , ny_block  , &
                            icellt    ,             &
                            indxti    , indxtj    , &
                            uvel      , vvel      , &
                            dxt       , dyt       , &
                            cxp       , cyp       , &
                            cxm       , cym       , &
                            tarear    ,             &
                            shear     , divu      , &
                            rdg_conv  , rdg_shear )

      endif

      end subroutine stress

!=======================================================================

! Computes the strain rates and internal stress components for T points
      
! author: JF Lemieux, ECCC
! Nov 2021      

      subroutine stress_T   (nx_block,   ny_block,   & 
                             ksub,       icellt,     & 
                             indxti,     indxtj,     &
                             uvelE,      vvelE,      &
                             uvelN,      vvelN,      &
                             dxN,        dyE,        &
                             dxT,        dyT,        &
                             tarear,     tinyarea,   & 
                             strength,               &
                             zetax2T,    etax2T,     &
                             stresspT,   stressmT,   & 
                             stress12T,              & 
                             shear,      divu,       & 
                             rdg_conv,   rdg_shear   )

        use ice_dyn_shared, only: strain_rates_T, deformations_T, &
                                  viscous_coeffs_and_rep_pressure_T
        
      integer (kind=int_kind), intent(in) :: & 
         nx_block, ny_block, & ! block dimensions
         ksub              , & ! subcycling step
         icellt                ! no. of cells where icetmask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxti   , & ! compressed index in i-direction
         indxtj       ! compressed index in j-direction


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
         tarear   , & ! 1/tarea
         tinyarea     ! puny*tarea

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         zetax2T  , & ! zetax2 = 2*zeta (bulk viscous coeff)
         etax2T   , & ! etax2  = 2*eta  (shear viscous coeff)
         stresspT , & ! sigma11+sigma22
         stressmT , & ! sigma11-sigma22
         stress12T    ! sigma12

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         shear    , & ! strain rate II component (1/s)
         divu     , & ! strain rate I component, velocity divergence (1/s)
         rdg_conv , & ! convergence term for ridging (1/s)
         rdg_shear    ! shear term for ridging (1/s)

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &
        divT, tensionT, shearT, DeltaT, & ! strain rates at T point
        rep_prsT                          ! replacement pressure at T point

      real(kind=dbl_kind), parameter :: capping = c1 ! of the viscous coef

      character(len=*), parameter :: subname = '(stress_T)'

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------


      do ij = 1, icellt
         i = indxti(ij)
         j = indxtj(ij)

      !-----------------------------------------------------------------
      ! strain rates at T point
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------

         call strain_rates_T (nx_block,   ny_block,   &
                              i,          j,          &
                              uvelE,      vvelE,      &
                              uvelN,      vvelN,      &
                              dxN,        dyE,        &
                              dxT,        dyT,        &
                              divT,       tensionT,   &
                              shearT,     DeltaT      )

      !-----------------------------------------------------------------
      ! viscous coefficients and replacement pressure at T point
      !-----------------------------------------------------------------

         call viscous_coeffs_and_rep_pressure_T (strength(i,j),           &
                                                 tinyarea(i,j),           &
                                                 DeltaT,                  &
                                                 zetax2T(i,j),etax2T(i,j),&
                                                 rep_prsT, capping        )
         
      !-----------------------------------------------------------------
      ! the stresses                            ! kg/s^2
      !-----------------------------------------------------------------

      ! NOTE: for comp. efficiency 2 x zeta and 2 x eta are used in the code 
         
         stresspT(i,j)  = (stresspT(i,j)*(c1-arlx1i*revp) + &
                           arlx1i*(zetax2T(i,j)*divT - rep_prsT)) * denom1

         stressmT(i,j)  = (stressmT(i,j)*(c1-arlx1i*revp) + &
                           arlx1i*etax2T(i,j)*tensionT) * denom1

         stress12T(i,j) = (stress12T(i,j)*(c1-arlx1i*revp) + &
                            arlx1i*p5*etax2T(i,j)*shearT) * denom1

      enddo                     ! ij

      !-----------------------------------------------------------------
      ! on last subcycle, save quantities for mechanical redistribution
      !-----------------------------------------------------------------
      if (ksub == ndte) then

         call deformations_T (nx_block  , ny_block  , &
                              icellt    ,             &
                              indxti    , indxtj    , &
                              uvelE,      vvelE,      &
                              uvelN,      vvelN,      &
                              dxN,        dyE,        &
                              dxT,        dyT,        &
                              tarear    ,             &
                              shear     , divu      , &
                              rdg_conv  , rdg_shear   )

      endif

    end subroutine stress_T

!=======================================================================

! Computes the strain rates and internal stress components for U points
      
! author: JF Lemieux, ECCC
! Nov 2021      

      subroutine stress_U   (nx_block,   ny_block,  & 
                             ksub,       icellu,    &  
                             indxui,     indxuj,    &
                             uvelE,      vvelE,     &
                             uvelN,      vvelN,     &
                             uvelU,      vvelU,     &
                             dxE,        dyN,       &
                             dxU,        dyU,       &
                             tarea,                 &
                             ratiodxN,   ratiodxNr, &
                             ratiodyE,   ratiodyEr, &
                             epm,  npm, hm, uvm,    &
                             zetax2T,    etax2T,    &
                             stresspU,   stressmU,  & 
                             stress12U            )

      use ice_dyn_shared, only: strain_rates_U, &
                                viscous_coeffs_and_rep_pressure_T2U
        
      integer (kind=int_kind), intent(in) :: & 
         nx_block, ny_block, & ! block dimensions
         ksub              , & ! subcycling step
         icellu                ! no. of cells where iceumask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxui   , & ! compressed index in i-direction
         indxuj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         uvelE    , & ! x-component of velocity (m/s) at the E point
         vvelE    , & ! y-component of velocity (m/s) at the E point
         uvelN    , & ! x-component of velocity (m/s) at the N point
         vvelN    , & ! y-component of velocity (m/s) at the N point
         uvelU    , & ! x-component of velocity (m/s) at the U point
         vvelU    , & ! y-component of velocity (m/s) at the U point
         dxE      , & ! width  of E-cell through the middle (m)
         dyN      , & ! height of N-cell through the middle (m)
         dxU      , & ! width  of U-cell through the middle (m)
         dyU      , & ! height of U-cell through the middle (m)
         tarea    , & ! area of T-cell (m^2)
         ratiodxN , & ! -dxN(i+1,j)/dxN(i,j) factor for BCs across coastline
         ratiodxNr, & ! -dxN(i,j)/dxN(i+1,j) factor for BCs across coastline
         ratiodyE , & ! -dyE(i,j+1)/dyE(i,j) factor for BCs across coastline
         ratiodyEr, & ! -dyE(i,j)/dyE(i,j+1) factor for BCs across coastline
         epm      , & ! E-cell mask
         npm      , & ! N-cell mask
         hm       , & ! T-cell mask
         uvm      , & ! U-cell mask
         zetax2T  , & ! 2*zeta at the T point
         etax2T       ! 2*eta at the T point
      
      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         stresspU , & ! sigma11+sigma22
         stressmU , & ! sigma11-sigma22
         stress12U    ! sigma12

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &
        divU, tensionU, shearU, DeltaU, & ! strain rates at U point
        zetax2U, etax2U, rep_prsU         ! replacement pressure at U point

      character(len=*), parameter :: subname = '(stress_U)'

      do ij = 1, icellu
         i = indxui(ij)
         j = indxuj(ij)

      !-----------------------------------------------------------------
      ! strain rates at T point
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------

         call strain_rates_U (nx_block,   ny_block,  &
                              i,          j,         &
                              uvelE,      vvelE,     &
                              uvelN,      vvelN,     &
                              uvelU,      vvelU,     &
                              dxE,        dyN,       &
                              dxU,        dyU,       &
                              ratiodxN,   ratiodxNr, &
                              ratiodyE,   ratiodyEr, &
                              epm,  npm,  uvm,       &
                              divU,       tensionU,  &
                              shearU,     DeltaU     )
         
      !-----------------------------------------------------------------
      ! viscous coefficients and replacement pressure at T point
      !-----------------------------------------------------------------

         call viscous_coeffs_and_rep_pressure_T2U (zetax2T(i  ,j  ), zetax2T(i  ,j+1), &
                                                   zetax2T(i+1,j+1), zetax2T(i+1,j  ), &
                                                   etax2T (i  ,j  ), etax2T (i  ,j+1), &
                                                   etax2T (i+1,j+1), etax2T (i+1,j  ), &
                                                   hm     (i  ,j  ), hm     (i  ,j+1), &
                                                   hm     (i+1,j+1), hm     (i+1,j  ), &
                                                   tarea  (i  ,j  ), tarea  (i  ,j+1), &
                                                   tarea  (i+1,j+1), tarea  (i+1,j  ), &
                                                   DeltaU,zetax2U, etax2U, rep_prsU)

      !-----------------------------------------------------------------
      ! the stresses                            ! kg/s^2
      !-----------------------------------------------------------------

      ! NOTE: for comp. efficiency 2 x zeta and 2 x eta are used in the code 
         
         stresspU(i,j)  = (stresspU(i,j)*(c1-arlx1i*revp) + &
                           arlx1i*(zetax2U*divU - rep_prsU)) * denom1

         stressmU(i,j)  = (stressmU(i,j)*(c1-arlx1i*revp) + &
                           arlx1i*etax2U*tensionU) * denom1

         stress12U(i,j) = (stress12U(i,j)*(c1-arlx1i*revp) + &
                            arlx1i*p5*etax2U*shearU) * denom1

      enddo                     ! ij

    end subroutine stress_U
    
!=======================================================================

! Computes divergence of stress tensor at the E or N point for the mom equation
      
! author: JF Lemieux, ECCC
! Nov 2021      

      subroutine div_stress  (nx_block,   ny_block,   & 
                              ksub,       icell,      & 
                              indxi,     indxj,       &
                              dxE_N,   dyE_N,         &
                              dxT_U,   dyT_U,         &
                              arear,                  &
                              stresspF1,   stressmF1, & 
                              stress12F1,             &
                              stresspF2,   stressmF2, &
                              stress12F2,             &
                              F1, F2,                 &
                              grid_location)

        
      integer (kind=int_kind), intent(in) :: & 
         nx_block, ny_block, & ! block dimensions
         ksub              , & ! subcycling step
         icell                 ! no. of cells where epm (or npm) = 1 

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxi   , & ! compressed index in i-direction
         indxj       ! compressed index in j-direction


      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         dxE_N , & ! width of E or N-cell through the middle (m)
         dyE_N , & ! height of E or N-cell through the middle (m)
         dxT_U , & ! width of T or U-cell through the middle (m)
         dyT_U , & ! height of T or U-cell through the middle (m)
         arear , & ! earear or narear
         stresspF1  , & ! stressp  (U or T) used for F1 calculation
         stressmF1  , & ! stressm  (U or T) used for F1 calculation 
         stress12F1 , & ! stress12 (U or T) used for F1 calculation 
         stresspF2  , & ! stressp  (U or T) used for F2 calculation 
         stressmF2  , & ! stressm  (U or T) used for F2 calculation 
         stress12F2     ! stress12 (U or T) used for F2 calculation 

      character(len=*), intent(in) :: &
         grid_location ! E (East) or N (North) ! TO BE IMPROVED!!!!
      
      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out) :: &
         F1      , & ! div of stress tensor for u component
         F2          ! div of stress tensor for v component

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      character(len=*), parameter :: subname = '(div_stress)'

!!! Instead of having the if statements below we could define for example
!    i+ci, j+cj where ci, cj would change with grid_position
      
      do ij = 1, icell
         i = indxi(ij)
         j = indxj(ij)

      !-----------------------------------------------------------------
      ! F1,F2 : div of stress tensor for u,v components
      !-----------------------------------------------------------------

         select case (trim(grid_location))
         case('E')
            
            F1(i,j) = arear(i,j) * &
                 ( p5 * dyE_N(i,j) * ( stresspF1(i+1,j)-stresspF1(i,j) )     &
                 + (p5/dyE_N(i,j)) * ( (dyT_U(i+1,j)**2) * stressmF1(i+1,j)  &
                                      -(dyT_U(i,j)**2)*stressmF1(i,j) )      &
                 + (c1/dxE_N(i,j)) * ( (dxT_U(i,j)**2) * stress12F1(i,j)     &
                                      -(dxT_U(i,j-1)**2)*stress12F1(i,j-1) ) )

            F2(i,j) = arear(i,j) * &
                 ( p5 * dxE_N(i,j) * ( stresspF2(i,j)-stresspF2(i,j-1) )     &
                 - (p5/dxE_N(i,j)) * ( (dxT_U(i,j)**2) * stressmF2(i,j)      &
                                      -(dxT_U(i,j-1)**2)*stressmF2(i,j-1) )  &
                 + (c1/dyE_N(i,j)) * ( (dyT_U(i+1,j)**2) * stress12F2(i+1,j) &
                                      -(dyT_U(i,j)**2)*stress12F2(i,j) ) )

         case('N')

            F1(i,j) = arear(i,j) * &
                 ( p5 * dyE_N(i,j) * ( stresspF1(i,j)-stresspF1(i-1,j) )     &
                 + (p5/dyE_N(i,j)) * ( (dyT_U(i,j)**2) * stressmF1(i,j)      &
                                      -(dyT_U(i-1,j)**2)*stressmF1(i-1,j) )  &
                 + (c1/dxE_N(i,j)) * ( (dxT_U(i,j+1)**2) * stress12F1(i,j+1) &
                                      -(dxT_U(i,j)**2)*stress12F1(i,j) ) )

            F2(i,j) = arear(i,j) * &
                 ( p5 * dxE_N(i,j) * ( stresspF2(i,j+1)-stresspF2(i,j) )     &
                 - (p5/dxE_N(i,j)) * ( (dxT_U(i,j+1)**2) * stressmF2(i,j+1)  &
                                      -(dxT_U(i,j)**2)*stressmF2(i,j) )      &
                 + (c1/dyE_N(i,j)) * ( (dyT_U(i,j)**2) * stress12F2(i,j)     &
                                      -(dyT_U(i-1,j)**2)*stress12F2(i-1,j) ) )
         case default
            call abort_ice(subname // ' unknown grid_location: ' // grid_location)
         end select

      enddo                     ! ij

    end subroutine div_stress
    
!=======================================================================

      end module ice_dyn_evp

!=======================================================================
