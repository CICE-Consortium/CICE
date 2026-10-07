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
          dynOff, dynCnt, dynMap, dynNbr, dynMaxCnt, &
          uvel_slot, vvel_slot, aiU_slot, fmU_slot, TbU_slot, &
          umassdti_slot, waterxU_slot, wateryU_slot, &
          forcexU_slot, forceyU_slot, &
          uocnU_slot, vocnU_slot, cdn_ocnU_slot, uarear_slot, dynIJ, &
          pack_to_slots, unpack_from_slots, &
          nbr_W, nbr_S, nbr_SW, nbr_E, nbr_N, nbr_NE, &
          skipT, skipU, skipN, skipE, &
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
         waterxN  (:,:)   , & ! for ocean stress calculation, x (m/s)
         wateryN  (:,:)   , & ! for ocean stress calculation, y (m/s)
         forcexN  (:,:)   , & ! work array: combined atm stress and ocn tilt, x
         forceyN  (:,:)   , & ! work array: combined atm stress and ocn tilt, y
         rheofactN(:,:,:) , & ! mult. factor = 1, set to 0 if aiN <= rheo_area_min
         nmassdti (:,:)       ! mass of N-cell/dte (kg/m^2 s)
! all c or d
      real (kind=dbl_kind), allocatable :: &
         uocnE    (:,:,:) , & ! i ocean current (m/s)
         vocnE    (:,:,:) , & ! j ocean current (m/s)
         ss_tltxE (:,:,:) , & ! sea surface slope, x-direction (m/m)
         ss_tltyE (:,:,:) , & ! sea surface slope, y-direction (m/m)
         cdn_ocnE (:,:,:) , & ! ocn drag coefficient
         waterxE  (:,:)   , & ! for ocean stress calculation, x (m/s)
         wateryE  (:,:)   , & ! for ocean stress calculation, y (m/s)
         forcexE  (:,:)   , & ! work array: combined atm stress and ocn tilt, x
         forceyE  (:,:)   , & ! work array: combined atm stress and ocn tilt, y
         rheofactE(:,:,:) , & ! mult. factor = 1, set to 0 if aiE <= rheo_area_min
         emassdti (:,:)       ! mass of E-cell/dte (kg/m^2 s)

      ! the i and j dimensions are flattened: element (i,j,iblk) of the block
      ! layout lives at ((j-1)*nx_block + i, iblk).
      ! TARGET is required: these are the targets of the rank-3 aliases below
      real (kind=dbl_kind), allocatable, target :: &
         nmass    (:,:)   , & ! total mass of ice and snow (N grid)
         emass    (:,:)   , & ! total mass of ice and snow (E grid)
         aiN      (:,:)   , & ! ice fraction on N-grid
         aiE      (:,:)       ! ice fraction on E-grid

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         ratiodxN    , & ! - dxN(i+1,j)   / dxN(i,j)
         ratiodyE    , & ! - dyE(i  ,j+1) / dyE(i,j)
         ratiodxNr   , & !   1 / ratiodxN
         ratiodyEr       !   1 / ratiodyE

      real (kind=dbl_kind), allocatable :: &
         strengthU(:,:,:) , & ! strength averaged to U points
         divergU  (:,:) , & ! div array on U points, differentiate from divu
         tensionU (:,:) , & ! tension array on U points
         shearU   (:,:) , & ! shear array on U points
         deltaU   (:,:) , & ! delta array on U points
         zetax2T  (:,:) , & ! zetax2 = 2*zeta (bulk viscosity)
         zetax2U  (:,:) , & ! zetax2T averaged to U points
         etax2T   (:,:) , & ! etax2  = 2*eta  (shear viscosity)
         etax2U   (:,:)     ! etax2T averaged to U points

      real (kind=dbl_kind), allocatable :: &
         uocnU    (:,:,:) , & ! i ocean current (m/s)
         vocnU    (:,:,:) , & ! j ocean current (m/s)
         ss_tltxU (:,:,:) , & ! sea surface slope, x-direction (m/m)
         ss_tltyU (:,:,:) , & ! sea surface slope, y-direction (m/m)
         cdn_ocnU (:,:,:) , & ! ocn drag coefficient
         waterxU  (:,:)   , & ! for ocean stress calculation, x (m/s)
         wateryU  (:,:)   , & ! for ocean stress calculation, y (m/s)
         forcexU  (:,:)   , & ! work array: combined atm stress and ocn tilt, x
         forceyU  (:,:)       ! work array: combined atm stress and ocn tilt, y

      ! the i and j dimensions are flattened: element (i,j,iblk) of the block
      ! layout lives at ((j-1)*nx_block + i, iblk)
      real (kind=dbl_kind), allocatable, target :: &
         umass    (:,:)   , & ! total mass of ice and snow (u grid)
         umassdti (:,:)       ! mass of U-cell/dte (kg/m^2 s)

      ! the i and j dimensions are flattened: element (i,j,iblk) of the block
      ! layout lives at ((j-1)*nx_block + i, iblk)
      real (kind=dbl_kind), allocatable, target :: &
         tmass    (:,:)       ! total mass of ice and snow (kg/m^2)

      public :: evp, init_evp

!=======================================================================

      contains

!=======================================================================
! Elastic-viscous-plastic dynamics driver
!
      subroutine init_evp
      use ice_blocks, only: get_block, nx_block, ny_block, block
      use ice_domain_size, only: max_blocks
      use ice_domain, only: nblocks, blocks_ice
      use ice_grid, only: grid_ice, dxN, dyE
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
         ijb                   , & ! 1D index of (i,j) in the flattened block layout
         ilo,ihi,jlo,jhi           ! beginning and end of physical domain


      call init_dyn_shared(dt_dyn)

      ! Guarded on grid_ice as well as evp_algorithm, to match the run-time
      ! call below.  dyn_evp1d_run sits inside "if (grid_ice == B)", but this
      ! init did not, so a C or CD run with shared_mem_1d used to build the
      ! whole 1d machinery -- the static gathers, convert_2d_1d_init, the halo
      ! pairing and some forty persistent arrays of length navel on the master
      ! task -- for a solver that is then never called.  The 1d solver has
      ! never been implemented for C/CD; until it is, do not pay for it.
      if (evp_algorithm == "shared_mem_1d" .and. grid_ice == "B") then
         call dyn_evp1d_init
      endif

      allocate( uocnU    (nx_block,ny_block,max_blocks), & ! i ocean current (m/s)
                vocnU    (nx_block,ny_block,max_blocks), & ! j ocean current (m/s)
                ss_tltxU (nx_block,ny_block,max_blocks), & ! sea surface slope, x-direction (m/m)
                ss_tltyU (nx_block,ny_block,max_blocks), & ! sea surface slope, y-direction (m/m)
                cdn_ocnU (nx_block,ny_block,max_blocks), & ! ocn drag coefficient
                tmass    (nx_block*ny_block,max_blocks), & ! total mass of ice and snow (kg/m^2)
                waterxU  (nx_block*ny_block,max_blocks), & ! for ocean stress calculation, x (m/s)
                wateryU  (nx_block*ny_block,max_blocks), & ! for ocean stress calculation, y (m/s)
                forcexU  (nx_block*ny_block,max_blocks), & ! work array: combined atm stress and ocn tilt, x
                forceyU  (nx_block*ny_block,max_blocks), & ! work array: combined atm stress and ocn tilt, y
                umass    (nx_block*ny_block,max_blocks), & ! total mass of ice and snow (u grid)
                umassdti (nx_block*ny_block,max_blocks), & ! mass of U-cell/dte (kg/m^2 s)
                stat=ierr)
      if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory B-Grid evp')

      if (grid_ice == 'CD' .or. grid_ice == 'C') then

         allocate( strengthU(nx_block,ny_block,max_blocks), &
                   divergU  (nx_block*ny_block,max_blocks), &
                   tensionU (nx_block*ny_block,max_blocks), &
                   shearU   (nx_block*ny_block,max_blocks), &
                   deltaU   (nx_block*ny_block,max_blocks), &
                   zetax2T  (nx_block*ny_block,max_blocks), &
                   zetax2U  (nx_block*ny_block,max_blocks), &
                   etax2T   (nx_block*ny_block,max_blocks), &
                   etax2U   (nx_block*ny_block,max_blocks), &
                   stat=ierr)
         if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory U evp')

         allocate( uocnN    (nx_block,ny_block,max_blocks), &
                   vocnN    (nx_block,ny_block,max_blocks), &
                   ss_tltxN (nx_block,ny_block,max_blocks), &
                   ss_tltyN (nx_block,ny_block,max_blocks), &
                   cdn_ocnN (nx_block,ny_block,max_blocks), &
                   waterxN  (nx_block*ny_block,max_blocks), &
                   wateryN  (nx_block*ny_block,max_blocks), &
                   forcexN  (nx_block*ny_block,max_blocks), &
                   forceyN  (nx_block*ny_block,max_blocks), &
                   aiN      (nx_block*ny_block,max_blocks), &
                   nmass    (nx_block*ny_block,max_blocks), &
                   nmassdti (nx_block*ny_block,max_blocks), &
                   rheofactN(nx_block,ny_block,max_blocks), &
                   stat=ierr)
         if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory N evp')

         allocate( uocnE    (nx_block,ny_block,max_blocks), &
                   vocnE    (nx_block,ny_block,max_blocks), &
                   ss_tltxE (nx_block,ny_block,max_blocks), &
                   ss_tltyE (nx_block,ny_block,max_blocks), &
                   cdn_ocnE (nx_block,ny_block,max_blocks), &
                   waterxE  (nx_block*ny_block,max_blocks), &
                   wateryE  (nx_block*ny_block,max_blocks), &
                   forcexE  (nx_block*ny_block,max_blocks), &
                   forceyE  (nx_block*ny_block,max_blocks), &
                   aiE      (nx_block*ny_block,max_blocks), &
                   emass    (nx_block*ny_block,max_blocks), &
                   emassdti (nx_block*ny_block,max_blocks), &
                   rheofactE(nx_block,ny_block,max_blocks), &
                   stat=ierr)
         if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory E evp')

         allocate( ratiodxN (nx_block,ny_block,max_blocks), &
                   ratiodyE (nx_block,ny_block,max_blocks), &
                   ratiodxNr(nx_block,ny_block,max_blocks), &
                   ratiodyEr(nx_block,ny_block,max_blocks), &
                   stat=ierr)
         if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory ratio')

         !$OMP PARALLEL DO PRIVATE(iblk,i,j,ijb,ilo,ihi,jlo,jhi,this_block)
         do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            do j = jlo, jhi
            do i = ilo, ihi
               ijb = (j-1)*nx_block + i
               ratiodxN (i,j,iblk) = - dxN(ijb+1       ,iblk) / dxN(ijb,iblk)
               ratiodyE (i,j,iblk) = - dyE(ijb+nx_block,iblk) / dyE(ijb,iblk)
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
          ice_HaloUpdate2DR8_compressed, ice_HaloUpdate3DR8_compressed, &
          ice_HaloDestroy, ice_HaloUpdate_stress
      use ice_blocks, only: block, get_block, nx_block, ny_block, ns_boundary_type
      use ice_domain, only: nblocks, blocks_ice, halo_info, maskhalo_dyn
      use ice_domain_size, only: max_blocks
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
      use ice_restoring, only: ice_restoring_halo, restore_ice
      use ice_state, only: aice, aiU, vice, vsno, uvel, vvel, uvelN, vvelN, &
          uvelE, vvelE, divu, shear, vort, &
          aice_init, aice0, aicen, vicen, strength
      use ice_timers, only: timer_dynamics, timer_bound, &
          ice_timer_start, ice_timer_stop, timer_evp
      use ice_dyn_shared, only: evp_algorithm, stack_fields, unstack_fields, &
          DminTarea, visc_method, deformations, deformationsC_T, deformationsCD_T, &
          strain_rates_U, dxhy, dyhx, cxp, cyp, cxm, cym, &
          iceTmask, iceUmask, iceEmask, iceNmask, &
          dyn_haloUpdate, fld2, fld3, fld4, fld2c
      ! the compressed layout, and the C/CD slot arrays the subcycle runs on
      use ice_dyn_shared, only: dynMap, dynOff, dynCnt, dynNbr
      use ice_dyn_shared, only: uvelE_slot, vvelE_slot, uvelN_slot, vvelN_slot, &
          shearU_slot, stresspT_slot, stressmT_slot, stress12T_slot, &
          stress12U_slot, stresspU_slot, stressmU_slot, &
          uvel_slot, vvel_slot, &
          strintxU_slot, strintyU_slot, taubxU_slot, taubyU_slot, &
          strintxE_slot, strintyE_slot, strintxN_slot, strintyN_slot, &
          taubxE_slot, taubyE_slot, taubxN_slot, taubyN_slot, &
          aiE_slot, cdn_ocnE_slot, uocnE_slot, vocnE_slot, &
          waterxE_slot, wateryE_slot, forcexE_slot, forceyE_slot, &
          emassdti_slot, fmE_slot, TbE_slot, &
          aiN_slot, cdn_ocnN_slot, uocnN_slot, vocnN_slot, &
          waterxN_slot, wateryN_slot, forcexN_slot, forceyN_slot, &
          nmassdti_slot, fmN_slot, TbN_slot
      use ice_dyn_shared, only: uarea_slot, dxN_slot, dyE_slot, &
          dxU_slot, dyT_slot, dxT_slot, dyU_slot, &
          dxE_slot, dyN_slot, epm_slot, npm_slot, uvm_slot, &
          tarea_slot, hm_slot, earea_slot, narea_slot, &
          uareaavg_slot, wsumN_slot, wsumE_slot, &
          wsumNW_slot, wsumSE_slot, wsumNE_slot, &
          zetax2T_slot, etax2T_slot, zetax2U_slot, etax2U_slot, &
          strength_slot, strengthU_slot
      use ice_dyn_shared, only: grid_average_X2YS_compressed, &
          grid_average_X2YA_compressed, fld3c, fld4c
      use ice_grid, only: pack_to_slots, unpack_from_slots
      use ice_boundary, only: ice_HaloUpdate2DR8_compressed, &
          ice_HaloUpdate3DR8_compressed
      use ice_domain, only: maskhalo_dyn
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

      ! i and j flattened: element (i,j,n) at ((j-1)*nx_block + i, n)
      ! compressed, indexed by block-local slot 1..cnt; the LAST element is
      ! the "no such neighbour" dead slot the stepu stencil falls back on
      real (kind=dbl_kind), dimension(dynMaxCnt+1,8):: &
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
         divergU  (:,:) = c0
         tensionU (:,:) = c0
         shearU   (:,:) = c0
         deltaU   (:,:) = c0
         zetax2T  (:,:) = c0
         zetax2U  (:,:) = c0
         etax2T   (:,:) = c0
         etax2U   (:,:) = c0

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
                         vsno    (:,:,iblk), tmask   (:,  iblk), &
                         tmass   (:,iblk),   iceTmask(:,iblk))

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
         call grid_average_X2Y('S', aice_init, 'T'          , aiE     , 'E')
         call grid_average_X2Y('S', cdn_ocn  , 'T'          , cdn_ocnE, 'E')
         call grid_average_X2Y('S', uocn     , grid_ocn_dynu, uocnE   , 'E')
         call grid_average_X2Y('S', vocn     , grid_ocn_dynv, vocnE   , 'E')
         call grid_average_X2Y('S', ss_tltx  , grid_ocn_dynu, ss_tltxE, 'E')
         call grid_average_X2Y('S', ss_tlty  , grid_ocn_dynv, ss_tltyE, 'E')
         call grid_average_X2Y('S', tmass    , 'T'          , nmass   , 'N')
         call grid_average_X2Y('S', aice_init, 'T'          , aiN     , 'N')
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
                            aiU       (:,iblk),   umass     (:,iblk),   &
                            umassdti  (:,iblk),   fcor_blk  (:,iblk)  , &
                            umask     (:,iblk)  ,                       &
                            uocnU     (:,:,iblk), vocnU     (:,:,iblk), &
                            strairxU  (:,:,iblk), strairyU  (:,:,iblk), &
                            ss_tltxU  (:,:,iblk), ss_tltyU  (:,:,iblk), &
                            iceTmask  (:,iblk),   iceUmask  (:,iblk),   &
                            fmU       (:,iblk), dt                  , &
                            strtltxU  (:,:,iblk), strtltyU  (:,:,iblk), &
                            strocnxU  (:,:,iblk), strocnyU  (:,:,iblk), &
                            strintxU  (:,:,iblk), strintyU  (:,:,iblk), &
                            taubxU    (:,iblk)  , taubyU    (:,iblk)  , &
                            waterxU   (:,iblk)  , wateryU   (:,iblk)  , &
                            forcexU   (:,iblk)  , forceyU   (:,iblk)  , &
                            stressp_1 (:,iblk), stressp_2 (:,iblk), &
                            stressp_3 (:,iblk), stressp_4 (:,iblk), &
                            stressm_1 (:,iblk), stressm_2 (:,iblk), &
                            stressm_3 (:,iblk), stressm_4 (:,iblk), &
                            stress12_1(:,iblk), stress12_2(:,iblk), &
                            stress12_3(:,iblk), stress12_4(:,iblk), &
                            uvel_init (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),   vvel_init (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),   &
                            uvel      (:,:,iblk), vvel      (:,:,iblk), &
                           dynMap(:,iblk), skipU(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)), &
                            TbU       (:,iblk),                         &
                            Xmassdti_slot = umassdti_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                            waterx_slot   = waterxU_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                            watery_slot   = wateryU_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                            forcex_slot   = forcexU_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                            forcey_slot   = forceyU_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))

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
                            aiU       (:,iblk),   umass     (:,iblk),   &
                            umassdti  (:,iblk),   fcor_blk  (:,iblk)  , &
                            umaskCD   (:,iblk)  ,                       &
                            uocnU     (:,:,iblk), vocnU     (:,:,iblk), &
                            strairxU  (:,:,iblk), strairyU  (:,:,iblk), &
                            ss_tltxU  (:,:,iblk), ss_tltyU  (:,:,iblk), &
                            iceTmask  (:,iblk),   iceUmask  (:,iblk),   &
                            fmU       (:,iblk), dt,                   &
                            strtltxU  (:,:,iblk), strtltyU  (:,:,iblk), &
                            strocnxU  (:,:,iblk), strocnyU  (:,:,iblk), &
                            strintxU  (:,:,iblk), strintyU  (:,:,iblk), &
                            taubxU    (:,iblk)  , taubyU    (:,iblk)  , &
                            waterxU   (:,iblk)  , wateryU   (:,iblk)  , &
                            forcexU   (:,iblk)  , forceyU   (:,iblk)  , &
                            stressp_1 (:,iblk), stressp_2 (:,iblk), &
                            stressp_3 (:,iblk), stressp_4 (:,iblk), &
                            stressm_1 (:,iblk), stressm_2 (:,iblk), &
                            stressm_3 (:,iblk), stressm_4 (:,iblk), &
                            stress12_1(:,iblk), stress12_2(:,iblk), &
                            stress12_3(:,iblk), stress12_4(:,iblk), &
                            uvel_init (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),   vvel_init (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),   &
                            uvel      (:,:,iblk), vvel      (:,:,iblk), &
                           dynMap(:,iblk), skipU(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)), &
                            TbU       (:,iblk))

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
                            aiN       (:,iblk),   nmass     (:,iblk),   &
                            nmassdti  (:,iblk),   fcorN_blk (:,iblk)  , &
                            nmask     (:,iblk)  ,                       &
                            uocnN     (:,:,iblk), vocnN     (:,:,iblk), &
                            strairxN  (:,:,iblk), strairyN  (:,:,iblk), &
                            ss_tltxN  (:,:,iblk), ss_tltyN  (:,:,iblk), &
                            iceTmask  (:,iblk),   iceNmask  (:,iblk),   &
                            fmN       (:,iblk), dt                  , &
                            strtltxN  (:,:,iblk), strtltyN  (:,:,iblk), &
                            strocnxN  (:,:,iblk), strocnyN  (:,:,iblk), &
                            strintxN  (:,:,iblk), strintyN  (:,:,iblk), &
                            taubxN    (:,iblk)  , taubyN    (:,iblk)  , &
                            waterxN   (:,iblk)  , wateryN   (:,iblk)  , &
                            forcexN   (:,iblk)  , forceyN   (:,iblk)  , &
                            stressp_1 (:,iblk), stressp_2 (:,iblk), &
                            stressp_3 (:,iblk), stressp_4 (:,iblk), &
                            stressm_1 (:,iblk), stressm_2 (:,iblk), &
                            stressm_3 (:,iblk), stressm_4 (:,iblk), &
                            stress12_1(:,iblk), stress12_2(:,iblk), &
                            stress12_3(:,iblk), stress12_4(:,iblk), &
                            uvelN_init(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),   vvelN_init(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),   &
                            uvelN     (:,:,iblk), vvelN     (:,:,iblk), &
                           dynMap(:,iblk), skipN(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)), &
                            TbN       (:,iblk), rheofactN (:,:,iblk))

            !-----------------------------------------------------------------
            ! more preparation for dynamics on E grid
            !-----------------------------------------------------------------

            call dyn_prep2 (nx_block,             ny_block,             &
                            ilo, ihi,             jlo, jhi,             &
                            icellT        (iblk), icellE        (iblk), &
                            indxTi      (:,iblk), indxTj      (:,iblk), &
                            indxEi      (:,iblk), indxEj      (:,iblk), &
                            aiE       (:,iblk),   emass     (:,iblk),   &
                            emassdti  (:,iblk),   fcorE_blk (:,iblk)  , &
                            emask     (:,iblk)  ,                       &
                            uocnE     (:,:,iblk), vocnE     (:,:,iblk), &
                            strairxE  (:,:,iblk), strairyE  (:,:,iblk), &
                            ss_tltxE  (:,:,iblk), ss_tltyE  (:,:,iblk), &
                            iceTmask  (:,iblk),   iceEmask  (:,iblk),   &
                            fmE       (:,iblk), dt                  , &
                            strtltxE  (:,:,iblk), strtltyE  (:,:,iblk), &
                            strocnxE  (:,:,iblk), strocnyE  (:,:,iblk), &
                            strintxE  (:,:,iblk), strintyE  (:,:,iblk), &
                            taubxE    (:,iblk)  , taubyE    (:,iblk)  , &
                            waterxE   (:,iblk)  , wateryE   (:,iblk)  , &
                            forcexE   (:,iblk)  , forceyE   (:,iblk)  , &
                            stressp_1 (:,iblk), stressp_2 (:,iblk), &
                            stressp_3 (:,iblk), stressp_4 (:,iblk), &
                            stressm_1 (:,iblk), stressm_2 (:,iblk), &
                            stressm_3 (:,iblk), stressm_4 (:,iblk), &
                            stress12_1(:,iblk), stress12_2(:,iblk), &
                            stress12_3(:,iblk), stress12_4(:,iblk), &
                            uvelE_init(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),   vvelE_init(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),   &
                            uvelE     (:,:,iblk), vvelE     (:,:,iblk), &
                           dynMap(:,iblk), skipE(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)), &
                            TbE       (:,iblk), rheofactE (:,:,iblk))

            do i=1,nx_block
            do j=1,ny_block
               if (.not.iceUmask((j-1)*nx_block+i,iblk)) then
                  stresspU ((j-1)*nx_block+i,iblk) = c0
                  stressmU ((j-1)*nx_block+i,iblk) = c0
                  stress12U((j-1)*nx_block+i,iblk) = c0
               endif
               if (.not.iceTmask((j-1)*nx_block+i,iblk)) then
                  stresspT ((j-1)*nx_block+i,iblk) = c0
                  stressmT ((j-1)*nx_block+i,iblk) = c0
                  stress12T((j-1)*nx_block+i,iblk) = c0
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
            uvelN(:,:,:) = uvelN(:,:,:)*reshape(npm, shape(uvelN))
            vvelE(:,:,:) = vvelE(:,:,:)*reshape(epm, shape(vvelE))
         endif

         call ice_timer_start(timer_bound)
         call ice_HaloUpdate (uvelN,           halo_info, &
                              field_loc_Nface, field_type_vector)
         call ice_HaloUpdate (vvelE,           halo_info, &
                              field_loc_Eface, field_type_vector)
         call ice_timer_stop(timer_bound)

         call grid_average_X2Y('A', uvelE, 'E', uvel, 'U')
         call grid_average_X2Y('A', vvelN, 'N', vvel, 'U')
         uvel(:,:,:) = uvel(:,:,:)*reshape(uvm, shape(uvel))
         vvel(:,:,:) = vvel(:,:,:)*reshape(uvm, shape(vvel))
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
            do iblk = 1, nblocks
               do j = 1, ny_block
               do i = 1, nx_block
                  if (iceUmask((j-1)*nx_block+i,iblk)) halomask(i,j,iblk) = 1
               enddo
               enddo
            enddo
         elseif (grid_ice == 'C' .or. grid_ice == 'CD') then
            !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block,i,j,ij) SCHEDULE(runtime)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo,jhi
               do i = ilo,ihi
                  ij = (j-1)*nx_block + i
                  if (iceTmask(ij         ,iblk) .or. &
                      iceTmask(ij-1       ,iblk) .or. &
                      iceTmask(ij+1       ,iblk) .or. &
                      iceTmask(ij-nx_block,iblk) .or. &
                      iceTmask(ij+nx_block,iblk)) then
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
                                                 hwater(:,iblk), TbU (:,iblk))
               enddo
               !$OMP END PARALLEL DO

            elseif ( seabed_stress_method == 'probabilistic' ) then
               !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
               do iblk = 1, nblocks
                  call seabed_stress_factor_prob (nx_block    , ny_block      ,                 &
                                                  icellT(iblk), indxTi(:,iblk), indxTj(:,iblk), &
                                                  icellU(iblk), indxUi(:,iblk), indxUj(:,iblk), &
                                                  aicen(:,:,:,iblk), vicen(:,:,:,iblk)        , &
                                                  hwater (:,iblk), TbU    (:,iblk))
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
                                                 hwater(:,iblk), TbE (:,iblk), &
                                                 grid_location='E')
                  call seabed_stress_factor_LKD (nx_block        , ny_block,       &
                                                 icellN    (iblk),                 &
                                                 indxNi  (:,iblk), indxNj(:,iblk), &
                                                 vice  (:,:,iblk), aice(:,:,iblk), &
                                                 hwater(:,iblk), TbN (:,iblk), &
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
                                                  hwater (:,iblk), TbU    (:,iblk)        , &
                                                  TbE    (:,iblk), TbN    (:,iblk)        , &
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

            call dyn_evp1d_run(dt        ,                                     &
                               stressp_1 , stressp_2, stressp_3  , stressp_4 , &
                               stressm_1 , stressm_2 , stressm_3 , stressm_4 , &
                               stress12_1, stress12_2, stress12_3, stress12_4, &
                               strength  ,                                     &
                               cdn_ocnU  , aiU       , uocnU     , vocnU     , &
                               waterxU   , wateryU   , forcexU   , forceyU   , &
                               umass     , strintxU  , strintyU  , &
                               Tbu       , taubxU    , taubyU    , uvel      , &
                               vvel      , iceTmask  , iceUmask)

         else ! evp_algorithm == standard_2d (Standard CICE)

            ! Pack the velocities into the compressed layout once, subcycle on
            ! them there, and unpack once at the end.  uvel/vvel themselves
            ! stay in the block layout for restart, history, transport and the
            ! couplers; only the ndte subcycles run compressed.
            !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
            do iblk = 1, nblocks
               call pack_to_slots(uvel(:,:,iblk), dynMap(:,iblk), uvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
               call pack_to_slots(vvel(:,:,iblk), dynMap(:,iblk), vvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
               ! constant through the subcycling, so packed once and never
               ! unpacked; the block-layout originals stay for the consumers
               ! outside the dynamics
               call pack_to_slots(aiU (:,iblk), dynMap(:,iblk), aiU_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
               call pack_to_slots(fmU (:,iblk), dynMap(:,iblk), fmU_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
               call pack_to_slots(TbU (:,iblk), dynMap(:,iblk), TbU_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
               ! umassdti, waterxU, wateryU, forcexU and forceyU are not packed
               ! here: dyn_prep2 already filled their compressed copies in the
               ! same loop that computed them, which costs one store instead of
               ! a second traversal of the block.
               call pack_to_slots(uocnU(:,:,iblk), dynMap(:,iblk), uocnU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
               call pack_to_slots(vocnU(:,:,iblk), dynMap(:,iblk), vocnU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
               call pack_to_slots(Cdn_ocnU(:,:,iblk), dynMap(:,iblk), cdn_ocnU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            enddo
            !$OMP END PARALLEL DO

            do ksub = 1,ndte        ! subcycling

               !$OMP PARALLEL DO PRIVATE(iblk,strtmp) SCHEDULE(runtime)
               do iblk = 1, nblocks

                  !-----------------------------------------------------------------
                  ! stress tensor equation, total surface stress
                  !-----------------------------------------------------------------
                  call stress (nx_block            , ny_block            , &
                                                     icellT        (iblk), &
                               indxTi      (:,iblk), indxTj      (:,iblk), &
                               uvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), vvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               dxT       (:, iblk), dyT       (:, iblk), &
                               dxhy      (:,iblk), dyhx      (:,iblk), &
                               cxp       (:,iblk), cyp       (:,iblk), &
                               cxm       (:,iblk), cym       (:,iblk), &
                                                     DminTarea (:,iblk)  , &
                               strength  (:,:,iblk),                       &
                               stressp_1 (:,iblk), stressp_2 (:,iblk), &
                               stressp_3 (:,iblk), stressp_4 (:,iblk), &
                               stressm_1 (:,iblk), stressm_2 (:,iblk), &
                               stressm_3 (:,iblk), stressm_4 (:,iblk), &
                               stress12_1(:,iblk), stress12_2(:,iblk), &
                               stress12_3(:,iblk), stress12_4(:,iblk), &
                               dynMap  (:,iblk),                      &
                               dynNbr(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk),:), &
                               strtmp    (:,:))

                  !-----------------------------------------------------------------
                  ! momentum equation
                  !-----------------------------------------------------------------
                  call stepu (&
                              dynCnt       (iblk),                     &
                              cdn_ocnU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),  &
                              skipU(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)), &
                              aiU_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), strtmp  (:,:),      &
                              uocnU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), vocnU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              waterxU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), wateryU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              forcexU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), forceyU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              umassdti_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), fmU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              uarear_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),    &
                              strintxU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              strintyU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              taubxU_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              taubyU_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              uvel_init(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),   vvel_init(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),  &
                              uvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), vvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              dynNbr(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk),:), &
                              TbU_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))

               enddo  ! iblk
               !$OMP END PARALLEL DO

               ! U fields at NE corner.  Exchanged in the compressed layout,
               ! so the velocities stay compressed across all ndte subcycles,
               ! and bundled, so the pair costs one message rather than two --
               ! the compressed counterpart of what stack_fields does for the
               ! block layout.
               call ice_timer_start(timer_bound)
               call stack_fields(uvel_slot, vvel_slot, fld2c)
               if (maskhalo_dyn) then
                  call ice_HaloUpdate3DR8_compressed(fld2c, halo_info_mask, &
                       field_loc_NEcorner, field_type_vector, dynMap, dynOff)
               else
                  call ice_HaloUpdate3DR8_compressed(fld2c, halo_info, &
                       field_loc_NEcorner, field_type_vector, dynMap, dynOff)
               endif
               call unstack_fields(fld2c, uvel_slot, vvel_slot)
               call ice_timer_stop(timer_bound)

               call ice_restoring_halo(setfld='velocity', &
                                       uvel_s=uvel_slot, vvel_s=vvel_slot)

            enddo  ! sub cycling

            !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
            do iblk = 1, nblocks
               call unpack_from_slots(uvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), uvel(:,:,iblk))
               call unpack_from_slots(vvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), vvel(:,:,iblk))
               ! strintxU/strintyU and taubxU/taubyU are outputs only; stepu
               ! writes them to slots every pass and this is where the block
               ! arrays history reads get their one and only write.
               call unpack_from_slots(strintxU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), strintxU(:,:,iblk))
               call unpack_from_slots(strintyU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), strintyU(:,:,iblk))
               call unpack_from_slots(taubxU_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), taubxU  (:,iblk))
               call unpack_from_slots(taubyU_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), taubyU  (:,iblk))
            enddo
            !$OMP END PARALLEL DO
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
                               dxT     (:, iblk), dyT      (:, iblk), &
                               dxU     (:, iblk), dyU      (:, iblk), &
                               cxp     (:,iblk), cyp      (:,iblk), &
                               cxm     (:,iblk), cym      (:,iblk), &
                               tarear  (:,  iblk), vort     (:,:,iblk), &
                               shear   (:,:,iblk), divu     (:,:,iblk), &
                               rdg_conv(:,:,iblk), rdg_shear(:,:,iblk) )
         enddo
         !$OMP END PARALLEL DO

      elseif (grid_ice == "C") then

         ! ------------------------------------------------------------------
         ! Pack once, subcycle compressed, unpack once -- the shape the EAP
         ! conversion uses.  The block arrays stay authoritative outside the
         ! subcycle, which is what restart, history and the couplers read.
         !
         ! shearU is not packed: strain_rates_U is the only thing
         ! that writes it and it writes straight into the slots.
         ! ------------------------------------------------------------------
         !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
         do iblk = 1, nblocks
            call pack_to_slots(uvelE(:,:,iblk), dynMap(:,iblk), uvelE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(vvelE(:,:,iblk), dynMap(:,iblk), vvelE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(uvelN(:,:,iblk), dynMap(:,iblk), uvelN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(vvelN(:,:,iblk), dynMap(:,iblk), vvelN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            ! ⚠ uvel and vvel too.  strain_rates_U reads them on the FIRST
            ! subcycle pass -- as uvelU/vvelU -- before the grid_average
            ! shifts write them later in the same pass.  Without this the
            ! first pass used whatever uvel_slot held from the last timestep,
            ! which is the previous values in a continuous run and zeros
            ! after a restart: the restart test could not reproduce.
            call pack_to_slots(uvel (:,:,iblk), dynMap(:,iblk), uvel_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(vvel (:,:,iblk), dynMap(:,iblk), vvel_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(stresspT (:,iblk), dynMap(:,iblk), stresspT_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(stressmT (:,iblk), dynMap(:,iblk), stressmT_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(stress12T(:,iblk), dynMap(:,iblk), stress12T_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(stress12U(:,iblk), dynMap(:,iblk), stress12U_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            ! constant through the subcycling, so packed once
            call pack_to_slots(strength(:,:,iblk), dynMap(:,iblk), strength_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            ! The step kernels read these every pass but nothing changes
            ! them during the subcycle, so one pack each is enough.
            call pack_to_slots(aiE(:,iblk), dynMap(:,iblk), aiE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(Cdn_ocnE(:,:,iblk), dynMap(:,iblk), cdn_ocne_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(uocnE(:,:,iblk), dynMap(:,iblk), uocnE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(vocnE(:,:,iblk), dynMap(:,iblk), vocnE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(waterxE(:,iblk), dynMap(:,iblk), waterxE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(wateryE(:,iblk), dynMap(:,iblk), wateryE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(forcexE(:,iblk), dynMap(:,iblk), forcexE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(forceyE(:,iblk), dynMap(:,iblk), forceyE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(emassdti(:,iblk), dynMap(:,iblk), emassdti_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(fmE(:,iblk), dynMap(:,iblk), fmE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(TbE(:,iblk), dynMap(:,iblk), TbE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(aiN(:,iblk), dynMap(:,iblk), aiN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(Cdn_ocnN(:,:,iblk), dynMap(:,iblk), cdn_ocnn_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(uocnN(:,:,iblk), dynMap(:,iblk), uocnN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(vocnN(:,:,iblk), dynMap(:,iblk), vocnN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(waterxN(:,iblk), dynMap(:,iblk), waterxN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(wateryN(:,iblk), dynMap(:,iblk), wateryN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(forcexN(:,iblk), dynMap(:,iblk), forcexN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(forceyN(:,iblk), dynMap(:,iblk), forceyN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(nmassdti(:,iblk), dynMap(:,iblk), nmassdti_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(fmN(:,iblk), dynMap(:,iblk), fmN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(TbN(:,iblk), dynMap(:,iblk), TbN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         enddo
         !$OMP END PARALLEL DO

         do ksub = 1,ndte        ! subcycling

            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks

            !-----------------------------------------------------------------
            ! strain rates at U point
            ! NOTE these are actually strain rates * area  (m^2/s)
            !-----------------------------------------------------------------
               call strain_rates_U (nx_block, ny_block      , &
                                    icellU      (iblk),                      &
                                    indxUi    (:,iblk), indxUj     (:,iblk), &
                                    dynNbr(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:), &
                                    dynMap    (:,iblk),                      &
                                    uvelE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                    vvelE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                    uvelN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                    vvelN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                    uvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                    vvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                    dxE_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                    dyN_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                    dxU     (:, iblk), dyU      (:, iblk), &
                                    ratiodxN(:,:,iblk), ratiodxNr(:,:,iblk), &
                                    ratiodyE(:,:,iblk), ratiodyEr(:,:,iblk), &
                                    epm_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                    npm_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                    divergU (:,iblk), tensionU (:,iblk), &
                                    shearU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                    deltaU   (:,iblk)  )

            enddo  ! iblk
            !$OMP END PARALLEL DO

            call ice_timer_start(timer_bound)
            if (maskhalo_dyn) then
               call ice_HaloUpdate2DR8_compressed(shearU_slot, halo_info_mask, &
                    field_loc_NEcorner, field_type_scalar, dynMap, dynOff)
            else
               call ice_HaloUpdate2DR8_compressed(shearU_slot, halo_info, &
                    field_loc_NEcorner, field_type_scalar, dynMap, dynOff)
            endif
            call ice_timer_stop(timer_bound)

            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks
               call stressC_T (nx_block           , ny_block            , &
                                                    icellT        (iblk), &
                               indxTi      (:,iblk), indxTj      (:,iblk), &
                               dynNbr(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:), &
                               dynMap      (:,iblk),                       &
                               uvelE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               vvelN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               dxN_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               dyE_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               dxT       (:, iblk), dyT       (:, iblk), &
                               uarea_slot   (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               uareaavg_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               DminTarea (:,iblk)  , &
                               strength_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               shearU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               zetax2T_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               etax2T_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               stresspT_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               stressmT_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               stress12T_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))

            enddo
            !$OMP END PARALLEL DO

            ! All four go compressed in ONE exchange.  zetax2T and etax2T
            ! have no block form any more -- nothing outside the dynamics
            ! reads them -- so the block halo that used to carry them is gone.
            call ice_timer_start(timer_bound)
            call stack_fields(zetax2T_slot, etax2T_slot, &
                              stresspT_slot, stressmT_slot, fld4c)
            if (maskhalo_dyn) then
               call ice_HaloUpdate3DR8_compressed(fld4c, halo_info_mask, &
                    field_loc_center, field_type_scalar, dynMap, dynOff)
            else
               call ice_HaloUpdate3DR8_compressed(fld4c, halo_info, &
                    field_loc_center, field_type_scalar, dynMap, dynOff)
            endif
            call unstack_fields(fld4c, zetax2T_slot, etax2T_slot, &
                                stresspT_slot, stressmT_slot)
            call ice_timer_stop(timer_bound)

            ! The T-to-U shift on the compressed layout.  'NE' is the compass
            ! direction of the T2US stencil; hm is a real 0/1 mask there, not
            ! a logical, and tarea is its weight.
            if (visc_method == 'avg_strength') then
               call grid_average_X2YS_compressed('NE', strength_slot, &
                    tarea_slot, hm_slot, wsumNE_slot, strengthU_slot)
            elseif (visc_method == 'avg_zeta') then
               call grid_average_X2YS_compressed('NE', etax2T_slot, &
                    tarea_slot, hm_slot, wsumNE_slot, etax2U_slot)
            endif

            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks
               call stressC_U (nx_block           , ny_block            , &
                                                    icellU        (iblk), &
                              indxUi      (:,iblk), indxUj      (:,iblk), &
                              dynMap      (:,iblk),                       &
                              uarea_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              etax2U_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              deltaU    (:,iblk), &
                              strengthU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              shearU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              stress12U_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            enddo
            !$OMP END PARALLEL DO

            call ice_timer_start(timer_bound)
            if (maskhalo_dyn) then
               call ice_HaloUpdate2DR8_compressed(stress12U_slot, halo_info_mask, &
                    field_loc_NEcorner, field_type_scalar, dynMap, dynOff)
            else
               call ice_HaloUpdate2DR8_compressed(stress12U_slot, halo_info, &
                    field_loc_NEcorner, field_type_scalar, dynMap, dynOff)
            endif
            call ice_timer_stop(timer_bound)

            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks

               call div_stress_Ex (nx_block            , ny_block            , &
                                                         icellE        (iblk), &
                                   indxEi      (:,iblk), indxEj      (:,iblk), &
                                   dynNbr(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:), &
                                   dynMap      (:,iblk),                       &
                                   dxE       (:, iblk), dyE       (:, iblk), &
                                   dxU_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                   dyT_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                   earear    (:,  iblk), rheofactE (:,:,iblk), &
                                   stresspT_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                   stressmT_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                   stress12U_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                   strintxE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1) )

               call div_stress_Ny (nx_block            , ny_block            , &
                                                         icellN        (iblk), &
                                   indxNi      (:,iblk), indxNj      (:,iblk), &
                                   dynNbr(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:), &
                                   dynMap      (:,iblk),                       &
                                   dxN       (:, iblk), dyN       (:, iblk), &
                                   dxT_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                   dyU_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                   narear    (:,  iblk), rheofactN (:,:,iblk), &
                                   stresspT_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                   stressmT_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                   stress12U_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                   strintyN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1) )

            enddo
            !$OMP END PARALLEL DO

            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks

                call stepu_C (dynCnt(iblk), skipE(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),        & ! u, E point
                              cdn_ocnE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              aiE_slot     (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              uocnE_slot   (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), vocnE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              waterxE_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), forcexE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              emassdti_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), fmE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              strintxE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              taubxE_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              uvelE_init(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),                         &
                              uvelE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              vvelE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              TbE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))

                call stepv_C (dynCnt(iblk), skipN(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),        & ! v, N point
                              cdn_ocnN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              aiN_slot     (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              uocnN_slot   (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), vocnN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              wateryN_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), forceyN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              nmassdti_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), fmN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              strintyN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              taubyN_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              vvelN_init(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),                         &
                              uvelN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              vvelN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              TbN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            enddo
            !$OMP END PARALLEL DO

            ! Shifts, masks and halo exchanges all on the compressed layout,
            ! so the velocities never leave slot form inside the subcycle.
            ! Directions are each stencil's compass sense: E2NA is NW, N2EA is
            ! SE, E2UA is N, N2UA is E; the weight is the source cell's area.
            call ice_timer_start(timer_bound)
            if (maskhalo_dyn) then
               call ice_HaloUpdate2DR8_compressed(uvelE_slot, halo_info_mask, &
                    field_loc_Eface, field_type_vector, dynMap, dynOff)
            else
               call ice_HaloUpdate2DR8_compressed(uvelE_slot, halo_info, &
                    field_loc_Eface, field_type_vector, dynMap, dynOff)
            endif
            if (maskhalo_dyn) then
               call ice_HaloUpdate2DR8_compressed(vvelN_slot, halo_info_mask, &
                    field_loc_Nface, field_type_vector, dynMap, dynOff)
            else
               call ice_HaloUpdate2DR8_compressed(vvelN_slot, halo_info, &
                    field_loc_Nface, field_type_vector, dynMap, dynOff)
            endif
            call ice_timer_stop(timer_bound)

            call grid_average_X2YA_compressed('NW', uvelE_slot, earea_slot, wsumNW_slot, uvelN_slot)
            call grid_average_X2YA_compressed('SE', vvelN_slot, narea_slot, wsumSE_slot, vvelE_slot)
            uvelN_slot(:) = uvelN_slot(:)*npm_slot(:)
            vvelE_slot(:) = vvelE_slot(:)*epm_slot(:)

            call ice_timer_start(timer_bound)
            if (maskhalo_dyn) then
               call ice_HaloUpdate2DR8_compressed(uvelN_slot, halo_info_mask, &
                    field_loc_Nface, field_type_vector, dynMap, dynOff)
            else
               call ice_HaloUpdate2DR8_compressed(uvelN_slot, halo_info, &
                    field_loc_Nface, field_type_vector, dynMap, dynOff)
            endif
            if (maskhalo_dyn) then
               call ice_HaloUpdate2DR8_compressed(vvelE_slot, halo_info_mask, &
                    field_loc_Eface, field_type_vector, dynMap, dynOff)
            else
               call ice_HaloUpdate2DR8_compressed(vvelE_slot, halo_info, &
                    field_loc_Eface, field_type_vector, dynMap, dynOff)
            endif
            call ice_timer_stop(timer_bound)

            call grid_average_X2YA_compressed('N', uvelE_slot, earea_slot, wsumN_slot, uvel_slot)
            call grid_average_X2YA_compressed('E', vvelN_slot, narea_slot, wsumE_slot, vvel_slot)

            uvel_slot(:) = uvel_slot(:)*uvm_slot(:)
            vvel_slot(:) = vvel_slot(:)*uvm_slot(:)

            ! U fields at NE corner, bundled into one compressed exchange
            call ice_timer_start(timer_bound)
            call stack_fields(uvel_slot, vvel_slot, fld2c)
            if (maskhalo_dyn) then
               call ice_HaloUpdate3DR8_compressed(fld2c, halo_info_mask, &
                    field_loc_NEcorner, field_type_vector, dynMap, dynOff)
            else
               call ice_HaloUpdate3DR8_compressed(fld2c, halo_info, &
                    field_loc_NEcorner, field_type_vector, dynMap, dynOff)
            endif
            call unstack_fields(fld2c, uvel_slot, vvel_slot)
            call ice_timer_stop(timer_bound)

            call ice_restoring_halo(setfld='velocity', &
                                    uvel_s=uvel_slot, vvel_s=vvel_slot)

            ! ⚠ stepu_C/stepv_C and the four grid_average_X2Y shifts above

         enddo                     ! subcycling

         ! Unpack what the subcycle carried compressed.  deformationsC_T below
         ! reads the block arrays, so this has to come first.
         !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
         do iblk = 1, nblocks
            call unpack_from_slots(stresspT_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), stresspT (:,iblk))
            call unpack_from_slots(stressmT_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), stressmT (:,iblk))
            call unpack_from_slots(stress12T_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), stress12T(:,iblk))
            call unpack_from_slots(stress12U_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), stress12U(:,iblk))
            call unpack_from_slots(shearU_slot   (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), shearU   (:,iblk))
            call unpack_from_slots(uvelE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), uvelE(:,:,iblk))
            call unpack_from_slots(vvelE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), vvelE(:,:,iblk))
            call unpack_from_slots(uvelN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), uvelN(:,:,iblk))
            call unpack_from_slots(vvelN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), vvelN(:,:,iblk))
            call unpack_from_slots(uvel_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), uvel (:,:,iblk))
            call unpack_from_slots(vvel_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), vvel (:,:,iblk))
            ! strint* and taub* are outputs only: the kernels write slots
            ! every pass and this is their one and only block write.
            call unpack_from_slots(strintxE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), strintxE(:,:,iblk))
            call unpack_from_slots(strintyN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), strintyN(:,:,iblk))
            call unpack_from_slots(taubxE_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), taubxE  (:,iblk))
            call unpack_from_slots(taubyN_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), taubyN  (:,iblk))
         enddo
         !$OMP END PARALLEL DO

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
                                  dxN     (:, iblk), dyE      (:, iblk), &
                                  dxT     (:, iblk), dyT      (:, iblk), &
                                  tarear  (:,  iblk), uarea    (:,  iblk), &
                                  shearU  (:,iblk), vort     (:,:,iblk), &
                                  shear   (:,:,iblk), divu     (:,:,iblk), &
                                  rdg_conv(:,:,iblk), rdg_shear(:,:,iblk))
         enddo
         !$OMP END PARALLEL DO

      elseif (grid_ice == "CD") then

         ! Pack once, subcycle compressed, unpack once -- see the C branch.
         !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
         do iblk = 1, nblocks
            call pack_to_slots(uvelE(:,:,iblk), dynMap(:,iblk), uvelE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(vvelE(:,:,iblk), dynMap(:,iblk), vvelE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(uvelN(:,:,iblk), dynMap(:,iblk), uvelN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(vvelN(:,:,iblk), dynMap(:,iblk), vvelN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            ! ⚠ uvel and vvel too.  strain_rates_U reads them on the FIRST
            ! subcycle pass -- as uvelU/vvelU -- before the grid_average
            ! shifts write them later in the same pass.  Without this the
            ! first pass used whatever uvel_slot held from the last timestep,
            ! which is the previous values in a continuous run and zeros
            ! after a restart: the restart test could not reproduce.
            call pack_to_slots(uvel (:,:,iblk), dynMap(:,iblk), uvel_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(vvel (:,:,iblk), dynMap(:,iblk), vvel_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(stresspT (:,iblk), dynMap(:,iblk), stresspT_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(stressmT (:,iblk), dynMap(:,iblk), stressmT_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(stress12T(:,iblk), dynMap(:,iblk), stress12T_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(stresspU (:,iblk), dynMap(:,iblk), stresspU_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(stressmU (:,iblk), dynMap(:,iblk), stressmU_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(stress12U(:,iblk), dynMap(:,iblk), stress12U_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            ! constant through the subcycling, so packed once
            call pack_to_slots(strength(:,:,iblk), dynMap(:,iblk), strength_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            ! The step kernels read these every pass but nothing changes
            ! them during the subcycle, so one pack each is enough.
            call pack_to_slots(aiE(:,iblk), dynMap(:,iblk), aiE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(Cdn_ocnE(:,:,iblk), dynMap(:,iblk), cdn_ocne_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(uocnE(:,:,iblk), dynMap(:,iblk), uocnE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(vocnE(:,:,iblk), dynMap(:,iblk), vocnE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(waterxE(:,iblk), dynMap(:,iblk), waterxE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(wateryE(:,iblk), dynMap(:,iblk), wateryE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(forcexE(:,iblk), dynMap(:,iblk), forcexE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(forceyE(:,iblk), dynMap(:,iblk), forceyE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(emassdti(:,iblk), dynMap(:,iblk), emassdti_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(fmE(:,iblk), dynMap(:,iblk), fmE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(TbE(:,iblk), dynMap(:,iblk), TbE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(aiN(:,iblk), dynMap(:,iblk), aiN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(Cdn_ocnN(:,:,iblk), dynMap(:,iblk), cdn_ocnn_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(uocnN(:,:,iblk), dynMap(:,iblk), uocnN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(vocnN(:,:,iblk), dynMap(:,iblk), vocnN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(waterxN(:,iblk), dynMap(:,iblk), waterxN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(wateryN(:,iblk), dynMap(:,iblk), wateryN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(forcexN(:,iblk), dynMap(:,iblk), forcexN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(forceyN(:,iblk), dynMap(:,iblk), forceyN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(nmassdti(:,iblk), dynMap(:,iblk), nmassdti_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(fmN(:,iblk), dynMap(:,iblk), fmN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(TbN(:,iblk), dynMap(:,iblk), TbN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         enddo
         !$OMP END PARALLEL DO

         do ksub = 1,ndte        ! subcycling

            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks
               call stressCD_T (nx_block            , ny_block            , &
                                                      icellT        (iblk), &
                                indxTi      (:,iblk), indxTj      (:,iblk), &
                                dynNbr(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:), &
                                dynMap      (:,iblk),                       &
                                uvelE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                vvelE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                uvelN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                vvelN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                dxN_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                dyE_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                dxT       (:, iblk), dyT       (:, iblk), &
                                                      DminTarea (:,iblk)  , &
                                strength_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                zetax2T_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                etax2T_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                stresspT_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                stressmT_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                stress12T_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1) )

            enddo
            !$OMP END PARALLEL DO

            ! calls ice_haloUpdate, controls bundles and masks
            call ice_timer_start(timer_bound)
            call stack_fields(zetax2T_slot, etax2T_slot, fld2c)
            if (maskhalo_dyn) then
               call ice_HaloUpdate3DR8_compressed(fld2c, halo_info_mask, &
                    field_loc_center, field_type_scalar, dynMap, dynOff)
            else
               call ice_HaloUpdate3DR8_compressed(fld2c, halo_info, &
                    field_loc_center, field_type_scalar, dynMap, dynOff)
            endif
            call unstack_fields(fld2c, zetax2T_slot, etax2T_slot)
            call ice_timer_stop(timer_bound)

            ! T-to-U on the compressed layout; 'NE' is the T2US stencil's
            ! compass direction, tarea its weight and hm its (real 0/1) mask.
            if (visc_method == 'avg_strength') then
               call grid_average_X2YS_compressed('NE', strength_slot, &
                    tarea_slot, hm_slot, wsumNE_slot, strengthU_slot)
            elseif (visc_method == 'avg_zeta') then
               call grid_average_X2YS_compressed('NE', zetax2T_slot, &
                    tarea_slot, hm_slot, wsumNE_slot, zetax2U_slot)
               call grid_average_X2YS_compressed('NE', etax2T_slot, &
                    tarea_slot, hm_slot, wsumNE_slot, etax2U_slot)
            endif

            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks
            !-----------------------------------------------------------------
            ! strain rates at U point
            ! NOTE these are actually strain rates * area  (m^2/s)
            !-----------------------------------------------------------------
               call strain_rates_U (nx_block, ny_block      , &
                                                         icellU       (iblk), &
                                    indxUi     (:,iblk), indxUj     (:,iblk), &
                                    dynNbr(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:), &
                                    dynMap     (:,iblk),                      &
                                    uvelE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                    vvelE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                    uvelN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                    vvelN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                    uvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                    vvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                    dxE_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                    dyN_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                    dxU      (:, iblk), dyU      (:, iblk), &
                                    ratiodxN (:,:,iblk), ratiodxNr(:,:,iblk), &
                                    ratiodyE (:,:,iblk), ratiodyEr(:,:,iblk), &
                                    epm_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                    npm_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                    divergU  (:,iblk), tensionU (:,iblk), &
                                    shearU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                    DeltaU   (:,iblk)  )

               call stressCD_U     (nx_block           , ny_block           , &
                                                         icellU       (iblk), &
                                    indxUi     (:,iblk), indxUj     (:,iblk), &
                                    dynMap     (:,iblk),                      &
                                    uarea_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                    zetax2U_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                    etax2U_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                    strengthU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                    divergU  (:,iblk), tensionU (:,iblk), &
                                    shearU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                    DeltaU   (:,iblk), &
                                    stresspU_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                    stressmU_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                    stress12U_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            enddo
            !$OMP END PARALLEL DO

            ! calls ice_haloUpdate, controls bundles and masks
            ! Six stresses, two compressed exchanges -- one per field
            ! location, which is why they cannot be merged into one.
            call ice_timer_start(timer_bound)
            call stack_fields(stresspT_slot, stressmT_slot, &
                              stress12T_slot, fld3c)
            if (maskhalo_dyn) then
               call ice_HaloUpdate3DR8_compressed(fld3c, halo_info_mask, &
                    field_loc_center, field_type_scalar, dynMap, dynOff)
            else
               call ice_HaloUpdate3DR8_compressed(fld3c, halo_info, &
                    field_loc_center, field_type_scalar, dynMap, dynOff)
            endif
            call unstack_fields(fld3c, stresspT_slot, stressmT_slot, &
                                stress12T_slot)

            call stack_fields(stresspU_slot, stressmU_slot, &
                              stress12U_slot, fld3c)
            if (maskhalo_dyn) then
               call ice_HaloUpdate3DR8_compressed(fld3c, halo_info_mask, &
                    field_loc_NEcorner, field_type_scalar, dynMap, dynOff)
            else
               call ice_HaloUpdate3DR8_compressed(fld3c, halo_info, &
                    field_loc_NEcorner, field_type_scalar, dynMap, dynOff)
            endif
            call unstack_fields(fld3c, stresspU_slot, stressmU_slot, &
                                stress12U_slot)
            call ice_timer_stop(timer_bound)

            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks

               call div_stress_Ex (nx_block            , ny_block            , &
                                                         icellE        (iblk), &
                                   indxEi      (:,iblk), indxEj      (:,iblk), &
                                   dynNbr(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:), &
                                   dynMap      (:,iblk),                       &
                                   dxE       (:, iblk), dyE       (:, iblk), &
                                   dxU_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                   dyT_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                   earear    (:,  iblk), rheofactE (:,:,iblk), &
                                   stresspT_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                   stressmT_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                   stress12U_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                   strintxE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1) )

               call div_stress_Ey (nx_block            , ny_block            , &
                                                         icellE        (iblk), &
                                   indxEi      (:,iblk), indxEj      (:,iblk), &
                                   dynNbr(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:), &
                                   dynMap      (:,iblk),                       &
                                   dxE       (:, iblk), dyE       (:, iblk), &
                                   dxU_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                   dyT_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                   earear    (:,  iblk), rheofactE (:,:,iblk), &
                                   stresspU_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                   stressmU_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                   stress12T_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                   strintyE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1) )

               call div_stress_Nx (nx_block            , ny_block            , &
                                                         icellN        (iblk), &
                                   indxNi      (:,iblk), indxNj      (:,iblk), &
                                   dynNbr(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:), &
                                   dynMap      (:,iblk),                       &
                                   dxN       (:, iblk), dyN       (:, iblk), &
                                   dxT_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                   dyU_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                   narear    (:,  iblk), rheofactN (:,:,iblk), &
                                   stresspU_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                   stressmU_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                   stress12T_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                   strintxN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1) )

               call div_stress_Ny (nx_block            , ny_block            , &
                                                         icellN        (iblk), &
                                   indxNi      (:,iblk), indxNj      (:,iblk), &
                                   dynNbr(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:), &
                                   dynMap      (:,iblk),                       &
                                   dxN       (:, iblk), dyN       (:, iblk), &
                                   dxT_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                   dyU_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                   narear    (:,  iblk), rheofactN (:,:,iblk), &
                                   stresspT_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                   stressmT_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                   stress12U_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                   strintyN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1) )

            enddo
            !$OMP END PARALLEL DO

            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks

               call stepuv_CD (dynCnt(iblk), skipE(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),       & ! E point
                               cdn_ocnE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               aiE_slot     (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               uocnE_slot   (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), vocnE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               waterxE_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), wateryE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               forcexE_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), forceyE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               emassdti_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), fmE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               strintxE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               strintyE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               taubxE_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               taubyE_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               uvelE_init(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),   vvelE_init(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),   &
                               uvelE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               vvelE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               TbE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))

               call stepuv_CD (dynCnt(iblk), skipN(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),       & ! N point
                               cdn_ocnN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               aiN_slot     (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               uocnN_slot   (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), vocnN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               waterxN_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), wateryN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               forcexN_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), forceyN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               nmassdti_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), fmN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               strintxN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               strintyN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               taubxN_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               taubyN_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               uvelN_init(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),   vvelN_init(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),   &
                               uvelN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               vvelN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               TbN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            enddo
            !$OMP END PARALLEL DO

            ! CD evolves all four velocities independently, so unlike the C
            ! branch there is no E2NA/N2EA pair here -- only the two shifts
            ! that build the U-point velocities.  Each face pair goes in one
            ! compressed bundle.
            call ice_timer_start(timer_bound)
            call stack_fields(uvelE_slot, vvelE_slot, fld2c)
            if (maskhalo_dyn) then
               call ice_HaloUpdate3DR8_compressed(fld2c, halo_info_mask, &
                    field_loc_Eface, field_type_vector, dynMap, dynOff)
            else
               call ice_HaloUpdate3DR8_compressed(fld2c, halo_info, &
                    field_loc_Eface, field_type_vector, dynMap, dynOff)
            endif
            call unstack_fields(fld2c, uvelE_slot, vvelE_slot)
            call stack_fields(uvelN_slot, vvelN_slot, fld2c)
            if (maskhalo_dyn) then
               call ice_HaloUpdate3DR8_compressed(fld2c, halo_info_mask, &
                    field_loc_Nface, field_type_vector, dynMap, dynOff)
            else
               call ice_HaloUpdate3DR8_compressed(fld2c, halo_info, &
                    field_loc_Nface, field_type_vector, dynMap, dynOff)
            endif
            call unstack_fields(fld2c, uvelN_slot, vvelN_slot)
            call ice_timer_stop(timer_bound)

            call grid_average_X2YA_compressed('N', uvelE_slot, earea_slot, wsumN_slot, uvel_slot)
            call grid_average_X2YA_compressed('E', vvelN_slot, narea_slot, wsumE_slot, vvel_slot)

            uvel_slot(:) = uvel_slot(:)*uvm_slot(:)
            vvel_slot(:) = vvel_slot(:)*uvm_slot(:)

            ! U fields at NE corner, bundled into one compressed exchange
            call ice_timer_start(timer_bound)
            call stack_fields(uvel_slot, vvel_slot, fld2c)
            if (maskhalo_dyn) then
               call ice_HaloUpdate3DR8_compressed(fld2c, halo_info_mask, &
                    field_loc_NEcorner, field_type_vector, dynMap, dynOff)
            else
               call ice_HaloUpdate3DR8_compressed(fld2c, halo_info, &
                    field_loc_NEcorner, field_type_vector, dynMap, dynOff)
            endif
            call unstack_fields(fld2c, uvel_slot, vvel_slot)
            call ice_timer_stop(timer_bound)

            call ice_restoring_halo(setfld='velocity', &
                                    uvel_s=uvel_slot, vvel_s=vvel_slot)


         enddo                     ! subcycling

         ! Unpack what the subcycle carried compressed; deformationsCD_T below
         ! reads the block arrays.
         !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
         do iblk = 1, nblocks
            call unpack_from_slots(stresspT_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), stresspT (:,iblk))
            call unpack_from_slots(stressmT_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), stressmT (:,iblk))
            call unpack_from_slots(stress12T_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), stress12T(:,iblk))
            call unpack_from_slots(stresspU_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), stresspU (:,iblk))
            call unpack_from_slots(stressmU_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), stressmU (:,iblk))
            call unpack_from_slots(stress12U_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), stress12U(:,iblk))
            call unpack_from_slots(shearU_slot   (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), shearU   (:,iblk))
            call unpack_from_slots(uvelE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), uvelE(:,:,iblk))
            call unpack_from_slots(vvelE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), vvelE(:,:,iblk))
            call unpack_from_slots(uvelN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), uvelN(:,:,iblk))
            call unpack_from_slots(vvelN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), vvelN(:,:,iblk))
            call unpack_from_slots(uvel_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), uvel (:,:,iblk))
            call unpack_from_slots(vvel_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), vvel (:,:,iblk))
            ! strint* and taub* are outputs only: the kernels write slots
            ! every pass and this is their one and only block write.
            call unpack_from_slots(strintxE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), strintxE(:,:,iblk))
            call unpack_from_slots(strintyE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), strintyE(:,:,iblk))
            call unpack_from_slots(strintxN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), strintxN(:,:,iblk))
            call unpack_from_slots(strintyN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), strintyN(:,:,iblk))
            call unpack_from_slots(taubxE_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), taubxE  (:,iblk))
            call unpack_from_slots(taubyE_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), taubyE  (:,iblk))
            call unpack_from_slots(taubxN_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), taubxN  (:,iblk))
            call unpack_from_slots(taubyN_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), taubyN  (:,iblk))
         enddo
         !$OMP END PARALLEL DO

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
                                   dxN     (:, iblk), dyE      (:, iblk), &
                                   dxT     (:, iblk), dyT      (:, iblk), &
                                   tarear  (:,  iblk), vort     (:,:,iblk), &
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
               aiU     (:,iblk),   fmU     (:,iblk), &
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
                  aiN     (:,iblk),   fmN     (:,iblk), &
                  strocnxN(:,:,iblk), strocnyN(:,:,iblk))

            call dyn_finish                               &
                 (nx_block          , ny_block          , &
                  icellE      (iblk), Cdn_ocnE(:,:,iblk), &
                  indxEi    (:,iblk), indxEj    (:,iblk), &
                  uvelE   (:,:,iblk), vvelE   (:,:,iblk), &
                  uocnE   (:,:,iblk), vocnE   (:,:,iblk), &
                  aiE     (:,iblk),   fmE     (:,iblk), &
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
                         Xmap,       nbr,        &
                         str )

      use ice_dyn_shared, only: strain_rates, visc_replpress

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellT                ! no. of cells where iceTmask = .true.

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxTi   , & ! compressed index in i-direction
         indxTj       ! compressed index in j-direction

      ! uvel and vvel are COMPRESSED: this block's slice, 1-based, so cells
      ! are 1..cnt and the last element is the block's dead slot
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         uvel     , & ! x-component of velocity (m/s)
         vvel         ! y-component of velocity (m/s)

      ! nbr(is,d) is the block-local slot one step in direction d from is
      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nbr          ! neighbour slots

      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(in) :: &
         strength , & ! ice strength (N/m)
         dxT      , & ! width of T-cell through the middle (m)
         dyT      , & ! height of T-cell through the middle (m)
         dxhy     , & ! 0.5*(HTE - HTW)
         dyhx     , & ! 0.5*(HTN - HTS)
         cyp      , & ! 1.5*HTE - 0.5*HTW
         cxp      , & ! 1.5*HTN - 0.5*HTS
         cym      , & ! 0.5*HTE - 1.5*HTW
         cxm          ! 0.5*HTN - 1.5*HTS

      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(inout) :: &
         stressp_1, stressp_2, stressp_3, stressp_4 , & ! sigma11+sigma22
         stressm_1, stressm_2, stressm_3, stressm_4 , & ! sigma11-sigma22
         stress12_1,stress12_2,stress12_3,stress12_4    ! sigma12

      ! str is COMPRESSED and indexed by BLOCK-LOCAL slot 1..cnt; the LAST
      ! element is the "no such neighbour" dead slot stepu falls back on
      real (kind=dbl_kind), dimension(:,:), intent(out), contiguous :: &
         str          ! stress combinations

      ! Xmap holds each cell's block-local slot, see ice_dyn_shared
      integer (kind=int_kind), dimension(nx_block*ny_block), intent(in) :: &
         Xmap

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij       , & ! ij is the compressed cell counter
         ijb            , & ! 1D index of (i,j) in the block
         is             , & ! block-local compressed slot of (i,j)
         isW, isS, isSW     ! its W, S and SW neighbour slots

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

      ! DminTarea is assumed shape, 1D over the block in Fortran storage order,
      ! i.e. element (i,j) is stored at (j-1)*nx_block+i
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         DminTarea    ! deltaminEVP*tarea

      character(len=*), parameter :: subname = '(stress)'

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      str(:,:) = c0

      do ij = 1, icellT
         i = indxTi(ij)
         j = indxTj(ij)
         ijb = (j-1)*nx_block + i
          is  = Xmap(ijb)

          ! each neighbour slot was fetched twice, once for uvel and once for
          ! vvel; read them once instead
          isW  = nbr(is,nbr_W )
          isS  = nbr(is,nbr_S )
          isSW = nbr(is,nbr_SW)

         !-----------------------------------------------------------------
         ! strain rates
         ! NOTE these are actually strain rates * area  (m^2/s)
         !-----------------------------------------------------------------

         call strain_rates (nx_block,   ny_block,   &
                            ijb,                    &
                            uvel(is),   uvel(isW),  &
                            uvel(isS),  uvel(isSW), &
                            vvel(is),   vvel(isW),  &
                            vvel(isS),  vvel(isSW), &
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

         call visc_replpress (strength(ijb), DminTarea(ijb), Deltane, &
                              zetax2ne, etax2ne, rep_prsne)

         call visc_replpress (strength(ijb), DminTarea(ijb), Deltanw, &
                              zetax2nw, etax2nw, rep_prsnw)

         call visc_replpress (strength(ijb), DminTarea(ijb), Deltasw, &
                              zetax2sw, etax2sw, rep_prssw)

         call visc_replpress (strength(ijb), DminTarea(ijb), Deltase, &
                              zetax2se, etax2se, rep_prsse)

         !-----------------------------------------------------------------
         ! the stresses                            ! kg/s^2
         ! (1) northeast, (2) northwest, (3) southwest, (4) southeast
         !-----------------------------------------------------------------

         ! NOTE: for comp. efficiency 2 x zeta and 2 x eta are used in the code
         stressp_1 (ijb) = (stressp_1 (ijb)*(c1-arlx1i*revp) &
                           + arlx1i*(zetax2ne*divune - rep_prsne)) * denom1
         stressp_2 (ijb) = (stressp_2 (ijb)*(c1-arlx1i*revp) &
                           + arlx1i*(zetax2nw*divunw - rep_prsnw)) * denom1
         stressp_3 (ijb) = (stressp_3 (ijb)*(c1-arlx1i*revp) &
                           + arlx1i*(zetax2sw*divusw - rep_prssw)) * denom1
         stressp_4 (ijb) = (stressp_4 (ijb)*(c1-arlx1i*revp) &
                           + arlx1i*(zetax2se*divuse - rep_prsse)) * denom1

         stressm_1 (ijb) = (stressm_1 (ijb)*(c1-arlx1i*revp) &
                           + arlx1i*etax2ne*tensionne) * denom1
         stressm_2 (ijb) = (stressm_2 (ijb)*(c1-arlx1i*revp) &
                           + arlx1i*etax2nw*tensionnw) * denom1
         stressm_3 (ijb) = (stressm_3 (ijb)*(c1-arlx1i*revp) &
                           + arlx1i*etax2sw*tensionsw) * denom1
         stressm_4 (ijb) = (stressm_4 (ijb)*(c1-arlx1i*revp) &
                           + arlx1i*etax2se*tensionse) * denom1

         stress12_1(ijb) = (stress12_1(ijb)*(c1-arlx1i*revp) &
                           + arlx1i*p5*etax2ne*shearne) * denom1
         stress12_2(ijb) = (stress12_2(ijb)*(c1-arlx1i*revp) &
                           + arlx1i*p5*etax2nw*shearnw) * denom1
         stress12_3(ijb) = (stress12_3(ijb)*(c1-arlx1i*revp) &
                           + arlx1i*p5*etax2sw*shearsw) * denom1
         stress12_4(ijb) = (stress12_4(ijb)*(c1-arlx1i*revp) &
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

!         stressp_1(ijb) = sign(max(abs(stressp_1(ijb)),puny),stressp_1(ijb))
!         stressp_2(ijb) = sign(max(abs(stressp_2(ijb)),puny),stressp_2(ijb))
!         stressp_3(ijb) = sign(max(abs(stressp_3(ijb)),puny),stressp_3(ijb))
!         stressp_4(ijb) = sign(max(abs(stressp_4(ijb)),puny),stressp_4(ijb))

!         stressm_1(ijb) = sign(max(abs(stressm_1(ijb)),puny),stressm_1(ijb))
!         stressm_2(ijb) = sign(max(abs(stressm_2(ijb)),puny),stressm_2(ijb))
!         stressm_3(ijb) = sign(max(abs(stressm_3(ijb)),puny),stressm_3(ijb))
!         stressm_4(ijb) = sign(max(abs(stressm_4(ijb)),puny),stressm_4(ijb))

!         stress12_1(ijb) = sign(max(abs(stress12_1(ijb)),puny),stress12_1(ijb))
!         stress12_2(ijb) = sign(max(abs(stress12_2(ijb)),puny),stress12_2(ijb))
!         stress12_3(ijb) = sign(max(abs(stress12_3(ijb)),puny),stress12_3(ijb))
!         stress12_4(ijb) = sign(max(abs(stress12_4(ijb)),puny),stress12_4(ijb))

         !-----------------------------------------------------------------
         ! combinations of the stresses for the momentum equation ! kg/s^2
         !-----------------------------------------------------------------

         ssigpn  = stressp_1(ijb) + stressp_2(ijb)
         ssigps  = stressp_3(ijb) + stressp_4(ijb)
         ssigpe  = stressp_1(ijb) + stressp_4(ijb)
         ssigpw  = stressp_2(ijb) + stressp_3(ijb)
         ssigp1  =(stressp_1(ijb) + stressp_3(ijb))*p055
         ssigp2  =(stressp_2(ijb) + stressp_4(ijb))*p055

         ssigmn  = stressm_1(ijb) + stressm_2(ijb)
         ssigms  = stressm_3(ijb) + stressm_4(ijb)
         ssigme  = stressm_1(ijb) + stressm_4(ijb)
         ssigmw  = stressm_2(ijb) + stressm_3(ijb)
         ssigm1  =(stressm_1(ijb) + stressm_3(ijb))*p055
         ssigm2  =(stressm_2(ijb) + stressm_4(ijb))*p055

         ssig12n = stress12_1(ijb) + stress12_2(ijb)
         ssig12s = stress12_3(ijb) + stress12_4(ijb)
         ssig12e = stress12_1(ijb) + stress12_4(ijb)
         ssig12w = stress12_2(ijb) + stress12_3(ijb)
         ssig121 =(stress12_1(ijb) + stress12_3(ijb))*p111
         ssig122 =(stress12_2(ijb) + stress12_4(ijb))*p111

         csigpne = p111*stressp_1(ijb) + ssigp2 + p027*stressp_3(ijb)
         csigpnw = p111*stressp_2(ijb) + ssigp1 + p027*stressp_4(ijb)
         csigpsw = p111*stressp_3(ijb) + ssigp2 + p027*stressp_1(ijb)
         csigpse = p111*stressp_4(ijb) + ssigp1 + p027*stressp_2(ijb)

         csigmne = p111*stressm_1(ijb) + ssigm2 + p027*stressm_3(ijb)
         csigmnw = p111*stressm_2(ijb) + ssigm1 + p027*stressm_4(ijb)
         csigmsw = p111*stressm_3(ijb) + ssigm2 + p027*stressm_1(ijb)
         csigmse = p111*stressm_4(ijb) + ssigm1 + p027*stressm_2(ijb)

         csig12ne = p222*stress12_1(ijb) + ssig122 &
                  + p055*stress12_3(ijb)
         csig12nw = p222*stress12_2(ijb) + ssig121 &
                  + p055*stress12_4(ijb)
         csig12sw = p222*stress12_3(ijb) + ssig122 &
                  + p055*stress12_1(ijb)
         csig12se = p222*stress12_4(ijb) + ssig121 &
                  + p055*stress12_2(ijb)

         str12ew = p5*dxT(ijb)*(p333*ssig12e + p166*ssig12w)
         str12we = p5*dxT(ijb)*(p333*ssig12w + p166*ssig12e)
         str12ns = p5*dyT(ijb)*(p333*ssig12n + p166*ssig12s)
         str12sn = p5*dyT(ijb)*(p333*ssig12s + p166*ssig12n)

         !-----------------------------------------------------------------
         ! for dF/dx (u momentum)
         !-----------------------------------------------------------------
         strp_tmp  = p25*dyT(ijb)*(p333*ssigpn  + p166*ssigps)
         strm_tmp  = p25*dyT(ijb)*(p333*ssigmn  + p166*ssigms)

         ! northeast (i,j)
         str(is,1) = -strp_tmp - strm_tmp - str12ew &
              + dxhy(ijb)*(-csigpne + csigmne) + dyhx(ijb)*csig12ne

         ! northwest (i+1,j)
         str(is,2) = strp_tmp + strm_tmp - str12we &
              + dxhy(ijb)*(-csigpnw + csigmnw) + dyhx(ijb)*csig12nw

         strp_tmp  = p25*dyT(ijb)*(p333*ssigps  + p166*ssigpn)
         strm_tmp  = p25*dyT(ijb)*(p333*ssigms  + p166*ssigmn)

         ! southeast (i,j+1)
         str(is,3) = -strp_tmp - strm_tmp + str12ew &
              + dxhy(ijb)*(-csigpse + csigmse) + dyhx(ijb)*csig12se

         ! southwest (i+1,j+1)
         str(is,4) = strp_tmp + strm_tmp + str12we &
              + dxhy(ijb)*(-csigpsw + csigmsw) + dyhx(ijb)*csig12sw

         !-----------------------------------------------------------------
         ! for dF/dy (v momentum)
         !-----------------------------------------------------------------
         strp_tmp  = p25*dxT(ijb)*(p333*ssigpe  + p166*ssigpw)
         strm_tmp  = p25*dxT(ijb)*(p333*ssigme  + p166*ssigmw)

         ! northeast (i,j)
         str(is,5) = -strp_tmp + strm_tmp - str12ns &
              - dyhx(ijb)*(csigpne + csigmne) + dxhy(ijb)*csig12ne

         ! southeast (i,j+1)
         str(is,6) = strp_tmp - strm_tmp - str12sn &
              - dyhx(ijb)*(csigpse + csigmse) + dxhy(ijb)*csig12se

         strp_tmp  = p25*dxT(ijb)*(p333*ssigpw  + p166*ssigpe)
         strm_tmp  = p25*dxT(ijb)*(p333*ssigmw  + p166*ssigme)

         ! northwest (i+1,j)
         str(is,7) = -strp_tmp + strm_tmp + str12ns &
              - dyhx(ijb)*(csigpnw + csigmnw) + dxhy(ijb)*csig12nw

         ! southwest (i+1,j+1)
         str(is,8) = strp_tmp - strm_tmp + str12sn &
              - dyhx(ijb)*(csigpsw + csigmsw) + dxhy(ijb)*csig12sw

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
                             nbr        , Xmap      , &
                             uvelE_s    , vvelN_s   , &
                             dxN_s      , dyE_s     , &
                             dxT        , dyT       , &
                             uarea_s    , uareaavg_s, &
                             DminTarea  ,             &
                             strength_s , shearU_s  , &
                             zetax2T_s  , etax2T_s  , &
                             stresspT_s , stressmT_s, &
                             stress12T_s)

      use ice_dyn_shared, only: strain_rates_Tdt_compressed, &
                                visc_replpress, e_factor, &
                                nbr_W, nbr_S, nbr_SW

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellT                ! no. of cells where iceTmask = .true.

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxTi   , & ! compressed index in i-direction
         indxTj       ! compressed index in j-direction

      ! nbr(is,d) is the block-local slot one step in direction d from slot is
      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nbr          ! neighbour slots

      ! Xmap is assumed shape, 1D over the block in Fortran storage order,
      ! i.e. element (i,j) is stored at (j-1)*nx_block+i
      integer (kind=int_kind), dimension (:), intent(in) :: &
         Xmap         ! static compressed slot of each cell

      ! Compressed, block-local slots 1..cnt with the dead slot last.  shearU
      ! and uarea are read at W, S and SW here; the velocities and their
      ! metrics at W and S inside strain_rates_Tdt_compressed.
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         uvelE_s  , & ! x-component of velocity (m/s) at the E point
         vvelN_s  , & ! y-component of velocity (m/s) at the N point
         dxN_s    , & ! width of N-cell through the middle (m)
         dyE_s    , & ! height of E-cell through the middle (m)
         uarea_s  , & ! area of U-cell (m^2)
         uareaavg_s,& ! precomputed sum of the 4 U areas round this T cell
         shearU_s     ! shearU local for this routine

      ! Pointwise, so block-indexed.
      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(in) :: &
         dxT      , & ! width of T-cell through the middle (m)
         dyT          ! height of T-cell through the middle (m)

      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         strength_s   ! ice strength (N/m), compressed

      ! Written here, read one cell away by div_stress_*, so compressed.
      real (kind=dbl_kind), dimension (:), intent(inout), contiguous :: &
         stresspT_s , & ! sigma11+sigma22
         stressmT_s , & ! sigma11-sigma22
         stress12T_s    ! sigma12

      ! Compressed because grid_average_X2YS_compressed shifts etax2T from T
      ! to U.  Nothing outside the dynamics reads zetax2T or etax2T -- checked
      ! across the tree -- so they have no block form at all.
      real (kind=dbl_kind), dimension (:), intent(inout), contiguous :: &
         zetax2T_s , & ! zetax2 = 2*zeta (bulk viscosity)
         etax2T_s      ! etax2  = 2*eta  (shear viscosity)

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij       , & ! ij is the compressed cell counter
         ijb            , & ! 1D index of (i,j) in the block
         is, isW, isS, isSW ! own, W, S and SW slots

      real (kind=dbl_kind), dimension (nx_block*ny_block) :: &
        divT      , & ! divergence at T point
        tensionT      ! tension at T point

      real (kind=dbl_kind) :: &
        shearTsqr , & ! strain rates squared at T point
        shearT    , & ! strain rate at T point
        DeltaT    , & ! delt at T point
        uareaavgr , & ! 1 / uarea avg
        rep_prsT      ! replacement pressure at T point

      ! DminTarea is assumed shape, 1D over the block in Fortran storage order,
      ! i.e. element (i,j) is stored at (j-1)*nx_block+i
      real (kind=dbl_kind), dimension (:), intent(in) :: &
         DminTarea    ! deltaminEVP*tarea

      character(len=*), parameter :: subname = '(stressC_T)'

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      call strain_rates_Tdt_compressed (nx_block , ny_block , &
                           icellT     ,                &
                           indxTi(:)  , indxTj  (:)  , &
                           nbr        , Xmap         , &
                           uvelE_s    , vvelN_s      , &
                           dxN_s      , dyE_s        , &
                           dxT   (:)  , dyT     (:)  , &
                           divT  (:)  , tensionT(:)  )

      do ij = 1, icellT
         i = indxTi(ij)
         j = indxTj(ij)
         ijb = (j-1)*nx_block + i

         !-----------------------------------------------------------------
         ! Square of shear strain rate at T obtained from interpolation of
         ! U point values (Bouillon et al., 2013, Kimmritz et al., 2016
         !-----------------------------------------------------------------

         is   = Xmap(ijb)
         isS  = nbr(is,nbr_S )
         isSW = nbr(is,nbr_SW)
         isW  = nbr(is,nbr_W )

         ! The denominator is precomputed in BLOCK space -- see ice_grid.
         ! Summing it here would contribute 0 for an unmapped neighbour, where
         ! the block path uses that cell's real area; uarea is geometry, not
         ! state, and is nonzero even over land.  Numerator order is self, S,
         ! SW, W and must stay that way.
         uareaavgr = c1/uareaavg_s(is)

         shearTsqr = (shearU_s(is  )**2 * uarea_s(is  )  &
                    + shearU_s(isS )**2 * uarea_s(isS )  &
                    + shearU_s(isSW)**2 * uarea_s(isSW)  &
                    + shearU_s(isW )**2 * uarea_s(isW )) &
                    * uareaavgr

         shearT    = (shearU_s(is  )    * uarea_s(is  )  &
                    + shearU_s(isS )    * uarea_s(isS )  &
                    + shearU_s(isSW)    * uarea_s(isSW)  &
                    + shearU_s(isW )    * uarea_s(isW )) &
                    * uareaavgr

         DeltaT = sqrt(divT(ijb)**2 + e_factor*(tensionT(ijb)**2 + shearTsqr))

         !-----------------------------------------------------------------
         ! viscosities and replacement pressure at T point
         !-----------------------------------------------------------------

         call visc_replpress (strength_s(is), DminTarea(ijb), DeltaT, &
                              zetax2T_s (is), etax2T_s(is), rep_prsT)

         !-----------------------------------------------------------------
         ! the stresses                            ! kg/s^2
         !-----------------------------------------------------------------

         ! NOTE: for comp. efficiency 2 x zeta and 2 x eta are used in the code

         stresspT_s(is)  = (stresspT_s (is)*(c1-arlx1i*revp) &
                          + arlx1i*(zetax2T_s(is)*divT(ijb) - rep_prsT)) * denom1

         stressmT_s(is)  = (stressmT_s (is)*(c1-arlx1i*revp) &
                          + arlx1i*etax2T_s(is)*tensionT(ijb)) * denom1

         stress12T_s(is) = (stress12T_s(is)*(c1-arlx1i*revp) &
                          + arlx1i*p5*etax2T_s(is)*shearT    ) * denom1

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
                             Xmap     ,            &
                             uarea_s  ,            &
                             etax2U_s , deltaU    ,&
                             strengthU_s, shearU_s,&
                             stress12U_s)

      use ice_dyn_shared, only: visc_replpress, &
                                visc_method, deltaminEVP

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellU                ! no. of cells where iceUmask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxUi   , & ! compressed index in i-direction
         indxUj       ! compressed index in j-direction

      ! Everything here is POINTWISE.  The compressed arrays appear only
      ! because other kernels read them at a neighbour -- shearU by stressC_T
      ! at W/S/SW, stress12U by div_stress_* -- so they live in slots and this
      ! routine has to address them that way.
      ! Xmap is assumed shape, 1D over the block in Fortran storage order,
      ! i.e. element (i,j) is stored at (j-1)*nx_block+i
      integer (kind=int_kind), dimension (:), intent(in) :: &
         Xmap         ! static compressed slot of each cell

      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(in) :: &
         deltaU       ! deltaU array

      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         uarea_s    , & ! area of U point
         shearU_s   , & ! shearU array
         etax2U_s   , & ! 2*eta at the U point
         strengthU_s    ! ice strength at the U point

      real (kind=dbl_kind), dimension (:), intent(inout), contiguous :: &
         stress12U_s  ! sigma12

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij       , & ! ij is the compressed cell counter
         ijb            , & ! 1D index of (i,j) in the block
         is                 ! compressed slot of (i,j)

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
            ijb = (j-1)*nx_block + i
            is  = Xmap(ijb)
            stress12U_s(is) = (stress12U_s(is)*(c1-arlx1i*revp) &
                             + arlx1i*p5*etax2U_s(is)*shearU_s(is)) * denom1
         enddo

      elseif (visc_method == 'avg_strength') then
         do ij = 1, icellU
            i = indxUi(ij)
            j = indxUj(ij)
            ijb = (j-1)*nx_block + i
            is  = Xmap(ijb)
            DminUarea = deltaminEVP*uarea_s(is)
            ! only need etax2U here, but other terms are calculated with etax2U
            ! minimal extra calculations here even though it seems like there is
            call visc_replpress (strengthU_s(is), DminUarea, deltaU(ijb), &
                                 lzetax2U      , letax2U  , lrep_prsU)
            stress12U_s(is) = (stress12U_s(is)*(c1-arlx1i*revp) &
                             + arlx1i*p5*letax2U*shearU_s(is)) * denom1
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
                             nbr     ,   Xmap     , &
                             uvelE_s ,   vvelE_s  , &
                             uvelN_s ,   vvelN_s  , &
                             dxN_s   ,   dyE_s    , &
                             dxT     ,   dyT      , &
                                         DminTarea, &
                             strength_s,            &
                             zetax2T_s,  etax2T_s , &
                             stresspT_s, stressmT_s, &
                             stress12T_s)

      use ice_dyn_shared, only: strain_rates_Tdtsd_compressed, &
                                visc_replpress

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellT                ! no. of cells where iceTmask = .true.

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxTi   , & ! compressed index in i-direction
         indxTj       ! compressed index in j-direction

      ! nbr(is,d) is the block-local slot one step in direction d from slot is
      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nbr          ! neighbour slots

      ! Xmap is assumed shape, 1D over the block in Fortran storage order,
      ! i.e. element (i,j) is stored at (j-1)*nx_block+i
      integer (kind=int_kind), dimension (:), intent(in) :: &
         Xmap         ! static compressed slot of each cell

      ! Compressed: strain_rates_Tdtsd_compressed reaches W and S on these.
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         uvelE_s  , & ! x-component of velocity (m/s) at the E point
         vvelE_s  , & ! y-component of velocity (m/s) at the E point
         uvelN_s  , & ! x-component of velocity (m/s) at the N point
         vvelN_s  , & ! y-component of velocity (m/s) at the N point
         dxN_s    , & ! width of N-cell through the middle (m)
         dyE_s        ! height of E-cell through the middle (m)

      ! Pointwise, so block-indexed.
      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(in) :: &
         dxT      , & ! width of T-cell through the middle (m)
         dyT          ! height of T-cell through the middle (m)

      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         strength_s   ! ice strength (N/m), compressed

      ! Compressed because grid_average_X2YS_compressed shifts etax2T from T
      ! to U.  Nothing outside the dynamics reads zetax2T or etax2T -- checked
      ! across the tree -- so they have no block form at all.
      real (kind=dbl_kind), dimension (:), intent(inout), contiguous :: &
         zetax2T_s , & ! zetax2 = 2*zeta (bulk viscosity)
         etax2T_s      ! etax2  = 2*eta  (shear viscosity)

      real (kind=dbl_kind), dimension (:), intent(inout), contiguous :: &
         stresspT_s , & ! sigma11+sigma22
         stressmT_s , & ! sigma11-sigma22
         stress12T_s    ! sigma12

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij       , & ! ij is the compressed cell counter
         ijb            , & ! 1D index of (i,j) in the block
         is                 ! compressed slot of (i,j)

      real (kind=dbl_kind), dimension (nx_block*ny_block) :: &
        divT      , & ! divergence at T point
        tensionT  , & ! tension at T point
        shearT    , & ! shear at T point
        DeltaT        ! delt at T point

      real (kind=dbl_kind) :: &
        rep_prsT      ! replacement pressure at T point

      ! DminTarea is assumed shape, 1D over the block in Fortran storage order,
      ! i.e. element (i,j) is stored at (j-1)*nx_block+i
      real (kind=dbl_kind), dimension (:), intent(in) :: &
         DminTarea    ! deltaminEVP*tarea

      character(len=*), parameter :: subname = '(stressCD_T)'

      !-----------------------------------------------------------------
      ! strain rates at T point
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------

      call strain_rates_Tdtsd_compressed (nx_block, ny_block , &
                           icellT     ,                &
                           indxTi(:)  , indxTj  (:)  , &
                           nbr        , Xmap         , &
                           uvelE_s    , vvelE_s      , &
                           uvelN_s    , vvelN_s      , &
                           dxN_s      , dyE_s        , &
                           dxT   (:)  , dyT     (:)  , &
                           divT  (:)  , tensionT(:)  , &
                           shearT(:)  , DeltaT  (:)  )

      do ij = 1, icellT
         i = indxTi(ij)
         j = indxTj(ij)
         ijb = (j-1)*nx_block + i
         is  = Xmap(ijb)

         !-----------------------------------------------------------------
         ! viscosities and replacement pressure at T point
         !-----------------------------------------------------------------

         call visc_replpress (strength_s(is), DminTarea(ijb), DeltaT(ijb), &
                              zetax2T_s (is), etax2T_s(is), rep_prsT)

         !-----------------------------------------------------------------
         ! the stresses                            ! kg/s^2
         !-----------------------------------------------------------------

         ! NOTE: for comp. efficiency 2 x zeta and 2 x eta are used in the code

         stresspT_s(is)  = (stresspT_s (is)*(c1-arlx1i*revp) &
                          + arlx1i*(zetax2T_s(is)*divT(ijb) - rep_prsT)) * denom1

         stressmT_s(is)  = (stressmT_s (is)*(c1-arlx1i*revp) &
                          + arlx1i*etax2T_s(is)*tensionT(ijb)) * denom1

         stress12T_s(is) = (stress12T_s(is)*(c1-arlx1i*revp) &
                          + arlx1i*p5*etax2T_s(is)*shearT(ijb)) * denom1

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
                             Xmap      ,           &
                             uarea_s   ,           &
                             zetax2U_s , etax2U_s, &
                             strengthU_s,          &
                             divergU   , tensionU, &
                             shearU_s  , deltaU  , &
                             stresspU_s, stressmU_s, &
                             stress12U_s)

      use ice_dyn_shared, only: visc_replpress, &
                                visc_method, deltaminEVP

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellU                ! no. of cells where iceUmask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxUi   , & ! compressed index in i-direction
         indxUj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         zetax2U_s  , & ! 2*zeta at U point, compressed
         etax2U_s   , & ! 2*eta at U point, compressed
         strengthU_s    ! ice strength at U point, compressed

      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(in) :: &
         divergU  , & ! div at U point
         tensionU , & ! tension at U point
         deltaU       ! delt at U point

      ! Xmap is assumed shape, 1D over the block in Fortran storage order,
      ! i.e. element (i,j) is stored at (j-1)*nx_block+i
      integer (kind=int_kind), dimension (:), intent(in) :: &
         Xmap         ! static compressed slot of each cell

      ! Pointwise here; compressed because div_stress_* read them a cell away.
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         uarea_s  , & ! area of U-cell (m^2)
         shearU_s     ! shear at U point

      real (kind=dbl_kind), dimension (:), intent(inout), contiguous :: &
         stresspU_s , & ! sigma11+sigma22
         stressmU_s , & ! sigma11-sigma22
         stress12U_s    ! sigma12

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij       , & ! ij is the compressed cell counter
         ijb            , & ! 1D index of (i,j) in the block
         is                 ! compressed slot of (i,j)

      real (kind=dbl_kind) :: &
        lzetax2U  , & ! bulk viscosity at U point
        letax2U   , & ! shear viscosity at U point
        lrep_prsU , & ! replacement pressure at U point
        DminUarea     ! Dmin on U

      character(len=*), parameter :: subname = '(stressCD_U)'

      do ij = 1, icellU
         i = indxUi(ij)
         j = indxUj(ij)
         ijb = (j-1)*nx_block + i
         is  = Xmap(ijb)

         !-----------------------------------------------------------------
         ! viscosities and replacement pressure at U point
         ! avg_zeta: Bouillon et al. 2013, C1 method of Kimmritz et al. 2016
         ! avg_strength: C2 method of Kimmritz et al. 2016
         !-----------------------------------------------------------------

         if (visc_method == 'avg_zeta') then
            lzetax2U  = zetax2U_s(is)
            letax2U   = etax2U_s(is)
            lrep_prsU = (c1-Ktens)/(c1+Ktens)*lzetax2U*deltaU(ijb)

         elseif (visc_method == 'avg_strength') then
            DminUarea = deltaminEVP*uarea_s(is)
            ! only need etax2U here, but other terms are calculated with etax2U
            ! minimal extra calculations here even though it seems like there is
            call visc_replpress (strengthU_s(is), DminUarea, deltaU(ijb), &
                                 lzetax2U      , letax2U  , lrep_prsU )
         endif

         !-----------------------------------------------------------------
         ! the stresses                            ! kg/s^2
         !-----------------------------------------------------------------

         ! NOTE: for comp. efficiency 2 x zeta and 2 x eta are used in the code

         stresspU_s(is)  = (stresspU_s (is)*(c1-arlx1i*revp) &
                          + arlx1i*(lzetax2U*divergU(ijb) - lrep_prsU)) * denom1

         stressmU_s(is)  = (stressmU_s (is)*(c1-arlx1i*revp) &
                          + arlx1i*letax2U*tensionU(ijb)) * denom1

         stress12U_s(is) = (stress12U_s(is)*(c1-arlx1i*revp) &
                          + arlx1i*p5*letax2U*shearU_s(is)) * denom1

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

      subroutine div_stress_Ex(nx_block, ny_block , &
                                         icell    , &
                               indxi   , indxj    , &
                               nbr     , Xmap     , &
                               dxE     , dyE      , &
                               dxU_s   , dyT_s    , &
                               arear   , rheofactE, &
                               stressp_s, stressm_s, &
                               stress12_s,          &
                               strintx_s )

      use ice_dyn_shared, only: nbr_E, nbr_S

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icell                 ! no. of cells where epm (or npm) = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxi   , & ! compressed index in i-direction
         indxj       ! compressed index in j-direction

      ! nbr(is,d) is the block-local slot one step in direction d from slot is
      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nbr          ! neighbour slots

      ! Every block array below is 1D over the block in Fortran storage order,
      ! i.e. element (i,j) lives at (j-1)*nx_block+i.  This routine used to
      ! declare its metrics and strint as (nx_block,ny_block) while callers
      ! passed flat slices, so it worked only by sequence association;
      ! indexing everything by ijb removes that implicit reshape.
      integer (kind=int_kind), dimension (:), intent(in) :: &
         Xmap         ! static compressed slot of each cell

      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(in) :: &
         dxE     , & ! width of E-cell through the middle (m)
         dyE     , & ! height of E-cell through the middle (m)
         rheofactE   ! mult. factor = 1, set to 0 if aiE <= rheo_area_min

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         arear       ! earear or narear

      ! Read at a neighbour -- dyT to the E, dxU to the S -- so compressed.
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         dxU_s   , & ! width of U-cell through the middle (m), compressed
         dyT_s       ! height of T-cell through the middle (m), compressed

      real (kind=dbl_kind), optional, dimension (:), intent(in), contiguous :: &
         stressp_s , & ! stressp  (U or T) used for strintx calculation
         stressm_s , & ! stressm  (U or T) used for strintx calculation
         stress12_s    ! stress12 (U or T) used for strintx calculation

      ! Output only: the step kernel reads it in the same pass and nothing
      ! else does, so it stays compressed.  The block array that history
      ! reads is written once, after the subcycle.
      real (kind=dbl_kind), optional, dimension (:), intent(out), contiguous :: &
         strintx_s   ! div of stress tensor for u component

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij, ijb, is, isA, isB

      character(len=*), parameter :: subname = '(div_stress_Ex)'

      do ij = 1, icell
         i = indxi(ij)
         j = indxj(ij)
         ijb = (j-1)*nx_block + i
         is  = Xmap(ijb)
         isA = nbr(is,nbr_E)   ! (i+1,j  )
         isB = nbr(is,nbr_S)   ! (i  ,j-1)
         strintx_s(is) = rheofactE(ijb) * arear(ijb) * &
              ( p5 * dyE(ijb)  * ( stressp_s(isA)  - stressp_s (is ) ) &
              + (p5/ dyE(ijb)) * ( (dyT_s(isA)**2) * stressm_s (isA)   &
                                  -(dyT_s(is )**2) * stressm_s (is ) ) &
              + (c1/ dxE(ijb)) * ( (dxU_s(is )**2) * stress12_s(is )   &
                                  -(dxU_s(isB)**2) * stress12_s(isB) ) )
      enddo

      end subroutine div_stress_Ex

!=======================================================================

      subroutine div_stress_Ey(nx_block, ny_block , &
                                         icell    , &
                               indxi   , indxj    , &
                               nbr     , Xmap     , &
                               dxE     , dyE      , &
                               dxU_s   , dyT_s    , &
                               arear   , rheofactE, &
                               stressp_s, stressm_s, &
                               stress12_s,          &
                               strinty_s )

      use ice_dyn_shared, only: nbr_S, nbr_E

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icell                 ! no. of cells where epm (or npm) = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxi   , & ! compressed index in i-direction
         indxj       ! compressed index in j-direction

      ! nbr(is,d) is the block-local slot one step in direction d from slot is
      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nbr          ! neighbour slots

      ! Every block array below is 1D over the block in Fortran storage order,
      ! i.e. element (i,j) lives at (j-1)*nx_block+i.  This routine used to
      ! declare its metrics and strint as (nx_block,ny_block) while callers
      ! passed flat slices, so it worked only by sequence association;
      ! indexing everything by ijb removes that implicit reshape.
      integer (kind=int_kind), dimension (:), intent(in) :: &
         Xmap         ! static compressed slot of each cell

      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(in) :: &
         dxE     , & ! width of E-cell through the middle (m)
         dyE     , & ! height of E-cell through the middle (m)
         rheofactE   ! mult. factor = 1, set to 0 if aiE <= rheo_area_min

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         arear       ! earear or narear

      ! Read at a neighbour -- dxU to the S, dyT to the E -- so compressed.
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         dxU_s   , & ! width of U-cell through the middle (m), compressed
         dyT_s       ! height of T-cell through the middle (m), compressed

      real (kind=dbl_kind), optional, dimension (:), intent(in), contiguous :: &
         stressp_s , & ! stressp  (U or T) used for strinty calculation
         stressm_s , & ! stressm  (U or T) used for strinty calculation
         stress12_s    ! stress12 (U or T) used for strinty calculation

      ! Output only: the step kernel reads it in the same pass and nothing
      ! else does, so it stays compressed.  The block array that history
      ! reads is written once, after the subcycle.
      real (kind=dbl_kind), optional, dimension (:), intent(out), contiguous :: &
         strinty_s   ! div of stress tensor for v component

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij, ijb, is, isA, isB

      character(len=*), parameter :: subname = '(div_stress_Ey)'

      do ij = 1, icell
         i = indxi(ij)
         j = indxj(ij)
         ijb = (j-1)*nx_block + i
         is  = Xmap(ijb)
         isA = nbr(is,nbr_S)   ! (i  ,j-1)
         isB = nbr(is,nbr_E)   ! (i+1,j  )
         strinty_s(is) = rheofactE(ijb) * arear(ijb) * &
              ( p5 * dxE(ijb)  * ( stressp_s(is )  - stressp_s (isA) ) &
              - (p5/ dxE(ijb)) * ( (dxU_s(is )**2) * stressm_s (is )   &
                                  -(dxU_s(isA)**2) * stressm_s (isA) ) &
              + (c1/ dyE(ijb)) * ( (dyT_s(isB)**2) * stress12_s(isB)   &
                                  -(dyT_s(is )**2) * stress12_s(is ) ) )
      enddo

      end subroutine div_stress_Ey

!=======================================================================

      subroutine div_stress_Nx(nx_block, ny_block , &
                                         icell    , &
                               indxi   , indxj    , &
                               nbr     , Xmap     , &
                               dxN     , dyN      , &
                               dxT_s   , dyU_s    , &
                               arear   , rheofactN, &
                               stressp_s, stressm_s, &
                               stress12_s,          &
                               strintx_s )

      use ice_dyn_shared, only: nbr_W, nbr_N

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icell                 ! no. of cells where epm (or npm) = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxi   , & ! compressed index in i-direction
         indxj       ! compressed index in j-direction

      ! nbr(is,d) is the block-local slot one step in direction d from slot is
      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nbr          ! neighbour slots

      ! Every block array below is 1D over the block in Fortran storage order,
      ! i.e. element (i,j) lives at (j-1)*nx_block+i.  This routine used to
      ! declare its metrics and strint as (nx_block,ny_block) while callers
      ! passed flat slices, so it worked only by sequence association;
      ! indexing everything by ijb removes that implicit reshape.
      integer (kind=int_kind), dimension (:), intent(in) :: &
         Xmap         ! static compressed slot of each cell

      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(in) :: &
         dxN     , & ! width of N-cell through the middle (m)
         dyN     , & ! height of N-cell through the middle (m)
         rheofactN   ! mult. factor = 1, set to 0 if aiN <= rheo_area_min

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         arear       ! earear or narear

      ! Read at a neighbour -- dyU to the W, dxT to the N -- so compressed.
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         dxT_s   , & ! width of T-cell through the middle (m), compressed
         dyU_s       ! height of U-cell through the middle (m), compressed

      real (kind=dbl_kind), optional, dimension (:), intent(in), contiguous :: &
         stressp_s , & ! stressp  (U or T) used for strintx calculation
         stressm_s , & ! stressm  (U or T) used for strintx calculation
         stress12_s    ! stress12 (U or T) used for strintx calculation

      ! Output only: the step kernel reads it in the same pass and nothing
      ! else does, so it stays compressed.  The block array that history
      ! reads is written once, after the subcycle.
      real (kind=dbl_kind), optional, dimension (:), intent(out), contiguous :: &
         strintx_s   ! div of stress tensor for u component

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij, ijb, is, isA, isB

      character(len=*), parameter :: subname = '(div_stress_Nx)'

      do ij = 1, icell
         i = indxi(ij)
         j = indxj(ij)
         ijb = (j-1)*nx_block + i
         is  = Xmap(ijb)
         isA = nbr(is,nbr_W)   ! (i-1,j  )
         isB = nbr(is,nbr_N)   ! (i  ,j+1)
         strintx_s(is) = rheofactN(ijb) * arear(ijb) * &
              ( p5 * dyN(ijb)  * ( stressp_s(is )  - stressp_s (isA) ) &
              + (p5/ dyN(ijb)) * ( (dyU_s(is )**2) * stressm_s (is )   &
                                  -(dyU_s(isA)**2) * stressm_s (isA) ) &
              + (c1/ dxN(ijb)) * ( (dxT_s(isB)**2) * stress12_s(isB)   &
                                  -(dxT_s(is )**2) * stress12_s(is ) ) )
      enddo

      end subroutine div_stress_Nx

!=======================================================================

      subroutine div_stress_Ny(nx_block, ny_block , &
                                         icell    , &
                               indxi   , indxj    , &
                               nbr     , Xmap     , &
                               dxN     , dyN      , &
                               dxT_s   , dyU_s    , &
                               arear   , rheofactN, &
                               stressp_s, stressm_s, &
                               stress12_s,          &
                               strinty_s )

      use ice_dyn_shared, only: nbr_N, nbr_W

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icell                 ! no. of cells where epm (or npm) = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxi   , & ! compressed index in i-direction
         indxj       ! compressed index in j-direction

      ! nbr(is,d) is the block-local slot one step in direction d from slot is
      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nbr          ! neighbour slots

      ! Every block array below is 1D over the block in Fortran storage order,
      ! i.e. element (i,j) lives at (j-1)*nx_block+i.  This routine used to
      ! declare its metrics and strint as (nx_block,ny_block) while callers
      ! passed flat slices, so it worked only by sequence association;
      ! indexing everything by ijb removes that implicit reshape.
      integer (kind=int_kind), dimension (:), intent(in) :: &
         Xmap         ! static compressed slot of each cell

      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(in) :: &
         dxN     , & ! width of N-cell through the middle (m)
         dyN     , & ! height of N-cell through the middle (m)
         rheofactN   ! mult. factor = 1, set to 0 if aiN <= rheo_area_min

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         arear       ! earear or narear

      ! Read at a neighbour -- dxT to the N, dyU to the W -- so compressed.
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         dxT_s   , & ! width of T-cell through the middle (m), compressed
         dyU_s       ! height of U-cell through the middle (m), compressed

      real (kind=dbl_kind), optional, dimension (:), intent(in), contiguous :: &
         stressp_s , & ! stressp  (U or T) used for strinty calculation
         stressm_s , & ! stressm  (U or T) used for strinty calculation
         stress12_s    ! stress12 (U or T) used for strinty calculation

      ! Output only: the step kernel reads it in the same pass and nothing
      ! else does, so it stays compressed.  The block array that history
      ! reads is written once, after the subcycle.
      real (kind=dbl_kind), optional, dimension (:), intent(out), contiguous :: &
         strinty_s   ! div of stress tensor for v component

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij, ijb, is, isA, isB

      character(len=*), parameter :: subname = '(div_stress_Ny)'

      do ij = 1, icell
         i = indxi(ij)
         j = indxj(ij)
         ijb = (j-1)*nx_block + i
         is  = Xmap(ijb)
         isA = nbr(is,nbr_N)   ! (i  ,j+1)
         isB = nbr(is,nbr_W)   ! (i-1,j  )
         strinty_s(is) = rheofactN(ijb) * arear(ijb) * &
              ( p5 * dxN(ijb)  * ( stressp_s(isA)  - stressp_s (is ) ) &
              - (p5/ dxN(ijb)) * ( (dxT_s(isA)**2) * stressm_s (isA)   &
                                  -(dxT_s(is )**2) * stressm_s (is ) ) &
              + (c1/ dyN(ijb)) * ( (dyU_s(is )**2) * stress12_s(is )   &
                                  -(dyU_s(isB)**2) * stress12_s(isB) ) )
      enddo

      end subroutine div_stress_Ny

!=======================================================================

      end module ice_dyn_evp

!=======================================================================
