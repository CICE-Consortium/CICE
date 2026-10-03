!=======================================================================
!
! Elastic-anisotropic sea ice dynamics model
! Computes ice velocity and deformation
!
! See:
!
! Wilchinsky, A.V. and D.L. Feltham (2006). Modelling the rheology of
! sea ice as a collection of diamond-shaped floes.
! Journal of Non-Newtonian Fluid Mechanics, 138(1), 22-32.
!
! Tsamados, M., D.L. Feltham, and A.V. Wilchinsky (2013). Impact on new
! anisotropic rheology on simulations of Arctic sea ice. JGR, 118, 91-107.
!
! authors: Michel Tsamados, CPOM
!          David Schroeder, CPOM

      module ice_dyn_eap

      use ice_kinds_mod
      use ice_blocks, only: nx_block, ny_block
      use ice_communicate, only: my_task, master_task
      use ice_domain_size, only: max_blocks, ncat
      use ice_constants, only: c0, c1, c2, c3, c4, c12, p1, p2, p5, &
          p001, p027, p055, p111, p166, p222, p25, p333
      use ice_fileunits, only: nu_diag, nu_dump_eap, nu_restart_eap
      use ice_exit, only: abort_ice
      use ice_flux, only: rdg_shear
!      use ice_timers, only:  &
!          ice_timer_start, ice_timer_stop, &
!          timer_tmp1, timer_tmp2, timer_tmp3, timer_tmp4, &
!          timer_tmp5, timer_tmp6, timer_tmp7, timer_tmp8, timer_tmp9

      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_parameters
      use icepack_intfc, only: icepack_ice_strength

      implicit none
      private
      public :: eap, init_eap, write_restart_eap, read_restart_eap

      ! Look-up table needed for calculating structure tensor
      integer (int_kind), parameter :: &
         nx_yield =  41, &
         ny_yield =  41, &
         na_yield =  21

      real (kind=dbl_kind), dimension (nx_yield,ny_yield,na_yield) :: &
         s11r, s12r, s22r, s11s, s12s, s22s

      real (kind=dbl_kind), dimension (:,:,:), allocatable :: &
         a11_1, a11_2, a11_3, a11_4,       & ! components of
         a12_1, a12_2, a12_3, a12_4          ! structure tensor

      ! Compressed copies of the state the EAP subcycle carries, so the whole
      ! subcycle can run on slots the way the B-grid one does: stress_eap and
      ! stepa both read and write these every subcycle, and the block arrays
      ! are only needed either side of the loop, for restart and history.
      !
      ! All three remaining EAP kernels -- stepa, update_stress_rdg and
      ! calc_ffrac -- are POINTWISE, so none of this needs a neighbour
      ! direction; the gain is contiguous iteration, as e151ee5 got for stepu,
      ! not removing an indirection.
      !
      ! The averaged a11/a12 are deliberately NOT here.  stress_eap reads
      ! a11_1..4 and a12_1..4, never the averages, so those two are outputs
      ! rather than subcycle state and stepa can write them through dynIJ(is),
      ! the way stepu writes taubxU.
      !
      ! Copy, never relayout -- the fff0dcd rule.
      real (kind=dbl_kind), dimension (:), allocatable :: &
         a11_1_slot, a11_2_slot, a11_3_slot, a11_4_slot, & ! structure tensor,
         a12_1_slot, a12_2_slot, a12_3_slot, a12_4_slot    ! compressed

      real (kind=dbl_kind), dimension (:), allocatable :: &
         stressp_1_slot , stressp_2_slot , stressp_3_slot , stressp_4_slot , &
         stressm_1_slot , stressm_2_slot , stressm_3_slot , stressm_4_slot , &
         stress12_1_slot, stress12_2_slot, stress12_3_slot, stress12_4_slot

      ! history
      real (kind=dbl_kind), dimension(:,:,:), allocatable, public :: &
         e11          , & ! components of strain rate tensor (1/s)
         e12          , &
         e22          , &
         yieldstress11, & ! components of yield stress tensor (kg/s^2)
         yieldstress12, &
         yieldstress22, &
         s11          , & ! components of stress tensor (kg/s^2)
         s12          , &
         s22          , &
         a11          , & ! components of structure tensor ()
         a12

      ! private for reuse, set in init_eap

      real (kind=dbl_kind) :: &
         puny, pi, pi2, piq, pih

      real (kind=dbl_kind), parameter :: &
         kfriction = 0.45_dbl_kind

      real (kind=dbl_kind), save :: &
         invdx, invdy, invda, invsin


!=======================================================================

      contains

!=======================================================================
! Elastic-anisotropic-plastic dynamics driver
! based on subroutine evp

      subroutine eap (dt)

#ifdef CICE_IN_NEMO
! Wind stress is set during this routine from the values supplied
! via NEMO (unless calc_strair is true).  These values are supplied
! rotated on u grid and multiplied by aice.  strairxT = 0 in this
! case so operations in dyn_prep1 are pointless but carried out to
! minimise code changes.
#endif

      use ice_arrays_column, only: Cdn_ocn
      use ice_boundary, only: ice_halo, ice_HaloMask, ice_HaloUpdate, &
          ice_HaloUpdate2DR8_compressed, ice_HaloUpdate3DR8_compressed, &
          ice_HaloDestroy
      use ice_blocks, only: block, get_block
      use ice_constants, only: field_loc_center, field_loc_NEcorner, &
          field_type_scalar, field_type_vector
      use ice_domain, only: nblocks, blocks_ice, halo_info, maskhalo_dyn
      use ice_dyn_shared, only: fcor_blk, ndte, dtei, &
          denom1, uvel_init, vvel_init, dynOff, dynCnt, dynMap, dynNbr, &
          nbr_W, nbr_S, nbr_SW, &
          dynMaxCnt, skipU, arlx1i, &
          uvel_slot, vvel_slot, aiU_slot, fmU_slot, TbU_slot, &
          umassdti_slot, waterxU_slot, wateryU_slot, &
          forcexU_slot, forceyU_slot, &
          uocnU_slot, vocnU_slot, cdn_ocnU_slot, uarear_slot, dynIJ, &
          pack_to_slots, unpack_from_slots, &
          dyn_prep1, dyn_prep2, stepu, dyn_finish, &
          seabed_stress_factor_LKD, seabed_stress_factor_prob, &
          seabed_stress_method, seabed_stress, &
          stack_fields, unstack_fields, iceTmask, iceUmask, &
          fld2, fld3, fld4, fld2c, dxhy, dyhx, cxp, cyp, cxm, cym
      use ice_flux, only: rdg_conv, strairxT, strairyT, &
          strairxU, strairyU, uocn, vocn, ss_tltx, ss_tlty, fmU, &
          strtltxU, strtltyU, strocnxU, strocnyU, strintxU, strintyU, taubxU, taubyU, &
          strax, stray, &
          TbU, hwater, &
          stressp_1, stressp_2, stressp_3, stressp_4, &
          stressm_1, stressm_2, stressm_3, stressm_4, &
          stress12_1, stress12_2, stress12_3, stress12_4
      use ice_grid, only: tmask, umask, dxT, dyT, dxU, dyU, &
          tarear, uarear, grid_average_X2Y, &
          grid_atm_dynu, grid_atm_dynv, grid_ocn_dynu, grid_ocn_dynv
      use ice_restoring, only: ice_restoring_halo, restore_ice
      use ice_state, only: aice, aiU, vice, vsno, uvel, vvel, divu, shear, vort, &
          aice_init, aice0, aicen, vicen, strength
      use ice_timers, only: timer_dynamics, timer_bound, &
          ice_timer_start, ice_timer_stop

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) :: &
         ksub       , & ! subcycle step
         iblk       , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, ij

      integer (kind=int_kind), dimension(max_blocks) :: &
         icellT     , & ! no. of cells where iceTmask = .true.
         icellU         ! no. of cells where iceUmask = .true.

      integer (kind=int_kind), dimension (nx_block*ny_block, max_blocks) :: &
         indxTi     , & ! compressed index in i-direction
         indxTj     , & ! compressed index in j-direction
         indxUi     , & ! compressed index in i-direction
         indxUj         ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         uocnU      , & ! i ocean current (m/s)
         vocnU      , & ! j ocean current (m/s)
         ss_tltxU   , & ! sea surface slope, x-direction (m/m)
         ss_tltyU   , & ! sea surface slope, y-direction (m/m)
         cdn_ocnU       ! ocn drag coefficient

      ! the i and j dimensions are flattened: element (i,j,iblk) of the block
      ! layout lives at ((j-1)*nx_block + i, iblk)
      real (kind=dbl_kind), dimension (nx_block*ny_block,max_blocks) :: &
         waterxU  , & ! for ocean stress calculation, x (m/s)
         wateryU  , & ! for ocean stress calculation, y (m/s)
         forcexU  , & ! work array: combined atm stress and ocn tilt, x
         forceyU      ! work array: combined atm stress and ocn tilt, y

      ! the i and j dimensions are flattened: element (i,j,iblk) of the block
      ! layout lives at ((j-1)*nx_block + i, iblk)
      real (kind=dbl_kind), dimension (nx_block*ny_block,max_blocks), target :: &
         umass      , & ! total mass of ice and snow (u grid)
         umassdti       ! mass of U-cell/dte (kg/m^2 s)

      ! the i and j dimensions are flattened: element (i,j,iblk) of the block
      ! layout lives at ((j-1)*nx_block + i, iblk)
      real (kind=dbl_kind), dimension (nx_block*ny_block,max_blocks), target :: &
         tmass          ! total mass of ice and snow (kg/m^2)

      ! compressed, indexed by block-local slot 1..cnt; the LAST element is
      ! the "no such neighbour" dead slot the stepu stencil falls back on
      real (kind=dbl_kind), dimension(dynMaxCnt+1,8):: &
         strtmp         ! stress combinations for momentum equation

      logical (kind=log_kind) :: &
         calc_strair

      integer (kind=int_kind), dimension (nx_block,ny_block,max_blocks) :: &
         halomask       ! ice mask for halo update

      type (ice_halo) :: &
         halo_info_mask !  ghost cell update info for masked halo

      type (block) :: &
         this_block     ! block information for current block

      character(len=*), parameter :: subname = '(eap)'

      call ice_timer_start(timer_dynamics) ! dynamics

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

       ! This call is needed only if dt changes during runtime.
!      call set_evp_parameters (dt)

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block) SCHEDULE(runtime)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            rdg_conv (i,j,iblk) = c0
            rdg_shear(i,j,iblk) = c0 ! always zero. Could be moved
            divu (i,j,iblk) = c0
            shear(i,j,iblk) = c0
            vort(i,j,iblk)  = c0
            e11(i,j,iblk) = c0
            e12(i,j,iblk) = c0
            e22(i,j,iblk) = c0
            s11(i,j,iblk) = c0
            s12(i,j,iblk) = c0
            s22(i,j,iblk) = c0
            yieldstress11(i,j,iblk) = c0
            yieldstress12(i,j,iblk) = c0
            yieldstress22(i,j,iblk) = c0
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

      !----------------------------------------------------------------
      ! Set wind stress to values supplied via NEMO or other forcing
      ! This wind stress is rotated on u grid and multiplied by aice
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
         ! Initialize structure tensor
         !-----------------------------------------------------------------

         do j = 1, ny_block
         do i = 1, nx_block
            if (.not.iceTmask((j-1)*nx_block+i,iblk)) then
               if (tmask((j-1)*nx_block+i,iblk)) then
                  ! structure tensor
                  a11_1(i,j,iblk) = p5
                  a11_2(i,j,iblk) = p5
                  a11_3(i,j,iblk) = p5
                  a11_4(i,j,iblk) = p5
                else
                  a11_1(i,j,iblk) = c0
                  a11_2(i,j,iblk) = c0
                  a11_3(i,j,iblk) = c0
                  a11_4(i,j,iblk) = c0
               endif
               a12_1(i,j,iblk) = c0
               a12_2(i,j,iblk) = c0
               a12_3(i,j,iblk) = c0
               a12_4(i,j,iblk) = c0
            endif                  ! iceTmask
         enddo                     ! i
         enddo                     ! j

         !-----------------------------------------------------------------
         ! ice strength
         ! New strength used in Ukita Moritz rheology
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
      enddo  ! iblk
      !$OMP END PARALLEL DO

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (strength,           halo_info, &
                           field_loc_center,   field_type_scalar)
      ! velocities may have changed in dyn_prep2
      call stack_fields(uvel, vvel, fld2)
      call ice_HaloUpdate (fld2,               halo_info, &
                           field_loc_NEcorner, field_type_vector)
      call unstack_fields(fld2, uvel, vvel)
      call ice_timer_stop(timer_bound)

      if (maskhalo_dyn) then
         call ice_timer_start(timer_bound)
         halomask = 0
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               if (iceUmask((j-1)*nx_block+i,iblk)) halomask(i,j,iblk) = 1
            enddo
            enddo
         enddo
         call ice_HaloUpdate (halomask,          halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_timer_stop(timer_bound)
         call ice_HaloMask(halo_info_mask, halo_info, halomask)
      endif

      !-----------------------------------------------------------------
      ! seabed stress factor TbU (TbU is part of Cb coefficient)
      !-----------------------------------------------------------------

      if (seabed_stress) then
         if ( seabed_stress_method == 'LKD' ) then
            !$OMP PARALLEL DO PRIVATE(iblk)  SCHEDULE(runtime)
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
               call seabed_stress_factor_prob (nx_block         , ny_block         , &
                                               icellT(iblk), indxTi(:,iblk), indxTj(:,iblk), &
                                               icellU(iblk), indxUi(:,iblk), indxUj(:,iblk), &
                                               aicen(:,:,:,iblk), vicen(:,:,:,iblk), &
                                               hwater (:,iblk), TbU    (:,iblk))
            enddo
            !$OMP END PARALLEL DO
         endif
      endif

      ! Pack the velocities once, subcycle compressed, unpack once.  uvel and
      ! vvel stay in the block layout for everything outside the dynamics.
      !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
      do iblk = 1, nblocks
         call pack_to_slots(uvel(:,:,iblk), dynMap(:,iblk), uvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(vvel(:,:,iblk), dynMap(:,iblk), vvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         ! constant through the subcycling, so packed once and never unpacked
         call pack_to_slots(aiU (:,iblk), dynMap(:,iblk), aiU_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(fmU (:,iblk), dynMap(:,iblk), fmU_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(TbU (:,iblk), dynMap(:,iblk), TbU_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(umassdti(:,iblk), dynMap(:,iblk), umassdti_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(waterxU(:,iblk), dynMap(:,iblk), waterxU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(wateryU(:,iblk), dynMap(:,iblk), wateryU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(forcexU(:,iblk), dynMap(:,iblk), forcexU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(forceyU(:,iblk), dynMap(:,iblk), forceyU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(uocnU(:,:,iblk), dynMap(:,iblk), uocnU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(vocnU(:,:,iblk), dynMap(:,iblk), vocnU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(Cdn_ocnU(:,:,iblk), dynMap(:,iblk), cdn_ocnU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))

         ! State the subcycle carries between stress_eap and stepa.  Packed
         ! here and unpacked after the loop; the block arrays are what restart
         ! and history read, so they stay authoritative outside the subcycle.
         call pack_to_slots(stressp_1(:,iblk), dynMap(:,iblk), stressp_1_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(stressp_2(:,iblk), dynMap(:,iblk), stressp_2_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(stressp_3(:,iblk), dynMap(:,iblk), stressp_3_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(stressp_4(:,iblk), dynMap(:,iblk), stressp_4_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(stressm_1(:,iblk), dynMap(:,iblk), stressm_1_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(stressm_2(:,iblk), dynMap(:,iblk), stressm_2_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(stressm_3(:,iblk), dynMap(:,iblk), stressm_3_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(stressm_4(:,iblk), dynMap(:,iblk), stressm_4_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(stress12_1(:,iblk), dynMap(:,iblk), stress12_1_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(stress12_2(:,iblk), dynMap(:,iblk), stress12_2_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(stress12_3(:,iblk), dynMap(:,iblk), stress12_3_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(stress12_4(:,iblk), dynMap(:,iblk), stress12_4_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(a11_1(:,:,iblk), dynMap(:,iblk), a11_1_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(a11_2(:,:,iblk), dynMap(:,iblk), a11_2_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(a11_3(:,:,iblk), dynMap(:,iblk), a11_3_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(a11_4(:,:,iblk), dynMap(:,iblk), a11_4_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(a12_1(:,:,iblk), dynMap(:,iblk), a12_1_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(a12_2(:,:,iblk), dynMap(:,iblk), a12_2_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(a12_3(:,:,iblk), dynMap(:,iblk), a12_3_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(a12_4(:,:,iblk), dynMap(:,iblk), a12_4_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
      enddo
      !$OMP END PARALLEL DO

      do ksub = 1,ndte        ! subcycling

         !-----------------------------------------------------------------
         ! stress tensor equation, total surface stress
         !-----------------------------------------------------------------

         !$OMP PARALLEL DO PRIVATE(iblk,strtmp) SCHEDULE(runtime)
         do iblk = 1, nblocks

!            call ice_timer_start(timer_tmp1,iblk)
            call stress_eap  (nx_block,             ny_block,             &
                              ksub,                 ndte,                 &
                              icellT        (iblk),                       &
                              indxTi      (:,iblk), indxTj      (:,iblk), &
                              arlx1i,               denom1,               &
                              uvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              vvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              dynNbr(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk),:), &
                              dxT       (:, iblk), dyT       (:, iblk), &
                              dxU       (:, iblk), dyU       (:, iblk), &
                              dxhy      (:,iblk), dyhx      (:,iblk), &
                              cxp       (:,iblk), cyp       (:,iblk), &
                              cxm       (:,iblk), cym       (:,iblk), &
                              tarear    (:,  iblk), strength  (:,:,iblk), &
                              a11_1_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              a11_2_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              a11_3_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              a11_4_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              a12_1_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              a12_2_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              a12_3_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              a12_4_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              stressp_1_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              stressp_2_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              stressp_3_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              stressp_4_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              stressm_1_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              stressm_2_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              stressm_3_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              stressm_4_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              stress12_1_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              stress12_2_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              stress12_3_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              stress12_4_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                              shear     (:,:,iblk), divu      (:,:,iblk), &
                              vort      (:,:,iblk),                       &
                              e11       (:,:,iblk), e12       (:,:,iblk), &
                              e22       (:,:,iblk),                       &
                              s11       (:,:,iblk), s12       (:,:,iblk), &
                              s22       (:,:,iblk),                       &
                              yieldstress11                   (:,:,iblk), &
                              yieldstress12                   (:,:,iblk), &
                              yieldstress22                   (:,:,iblk), &
!                              rdg_conv  (:,:,iblk), rdg_shear (:,:,iblk), &
                              rdg_conv  (:,:,iblk), &
                              dynMap    (:,iblk), strtmp (:,:))
!            call ice_timer_stop(timer_tmp1,iblk)

            !-----------------------------------------------------------------
            ! momentum equation
            !-----------------------------------------------------------------

!            call ice_timer_start(timer_tmp2,iblk)
            call stepu (nx_block,            ny_block,            &
                        dynCnt       (iblk), ksub,                &
                        cdn_ocnU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),  &
                        skipU(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)), &
                        dynIJ(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)), &
                        aiU_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), strtmp   (:,:),      &
                        uocnU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), vocnU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                        waterxU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), wateryU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                        forcexU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), forceyU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                        umassdti_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), fmU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                        uarear_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),    &
                        strintxU (:,:,iblk), strintyU (:,:,iblk), &
                        taubxU   (:,iblk)  , taubyU   (:,iblk)  , &
                        uvel_init(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),   vvel_init(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),   &
                        uvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), vvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                        dynNbr(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk),:), &
                        TbU_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
!            call ice_timer_stop(timer_tmp2,iblk)

            !-----------------------------------------------------------------
            ! evolution of structure tensor A
            !-----------------------------------------------------------------

!            call ice_timer_start(timer_tmp3,iblk)
            if (mod(ksub,10) == 1) then ! only called every 10th timestep
            call stepa (nx_block            , ny_block            , &
                        dtei                , icellT        (iblk), &
                        indxTi      (:,iblk), indxTj      (:,iblk), &
                        dynMap      (:,iblk),                       &
                        a11       (:,:,iblk), a12       (:,:,iblk), &
                        a11_1_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                        a11_2_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                        a11_3_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                        a11_4_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                        a12_1_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                        a12_2_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                        a12_3_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                        a12_4_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                        stressp_1_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                        stressp_2_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                        stressp_3_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                        stressp_4_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                        stressm_1_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                        stressm_2_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                        stressm_3_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                        stressm_4_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                        stress12_1_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                        stress12_2_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                        stress12_3_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                        stress12_4_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            endif
!            call ice_timer_stop(timer_tmp3,iblk)
         enddo
         !$OMP END PARALLEL DO

         ! exchanged compressed and bundled, so the pair costs one message
         ! rather than two -- the compressed counterpart of stack_fields
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

         ! ice_restoring_halo works on the block layout and returns at once
         ! unless restoring is on, so convert around it only in that case
         if (restore_ice) then
            !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
            do iblk = 1, nblocks
               call unpack_from_slots(uvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), uvel(:,:,iblk))
               call unpack_from_slots(vvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), vvel(:,:,iblk))
            enddo
            !$OMP END PARALLEL DO
            call ice_restoring_halo(setfld='velocity')
            !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
            do iblk = 1, nblocks
               call pack_to_slots(uvel(:,:,iblk), dynMap(:,iblk), uvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
               call pack_to_slots(vvel(:,:,iblk), dynMap(:,iblk), vvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            enddo
            !$OMP END PARALLEL DO
         endif

      enddo                     ! subcycling

      !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
      do iblk = 1, nblocks
         call unpack_from_slots(uvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), uvel(:,:,iblk))
         call unpack_from_slots(vvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), vvel(:,:,iblk))

         ! hand the subcycle state back to the block layout for restart,
         ! history and the diagnostics below
         call unpack_from_slots(stressp_1_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), stressp_1(:,iblk))
         call unpack_from_slots(stressp_2_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), stressp_2(:,iblk))
         call unpack_from_slots(stressp_3_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), stressp_3(:,iblk))
         call unpack_from_slots(stressp_4_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), stressp_4(:,iblk))
         call unpack_from_slots(stressm_1_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), stressm_1(:,iblk))
         call unpack_from_slots(stressm_2_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), stressm_2(:,iblk))
         call unpack_from_slots(stressm_3_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), stressm_3(:,iblk))
         call unpack_from_slots(stressm_4_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), stressm_4(:,iblk))
         call unpack_from_slots(stress12_1_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), stress12_1(:,iblk))
         call unpack_from_slots(stress12_2_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), stress12_2(:,iblk))
         call unpack_from_slots(stress12_3_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), stress12_3(:,iblk))
         call unpack_from_slots(stress12_4_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), stress12_4(:,iblk))
         call unpack_from_slots(a11_1_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), a11_1(:,:,iblk))
         call unpack_from_slots(a11_2_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), a11_2(:,:,iblk))
         call unpack_from_slots(a11_3_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), a11_3(:,:,iblk))
         call unpack_from_slots(a11_4_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), a11_4(:,:,iblk))
         call unpack_from_slots(a12_1_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), a12_1(:,:,iblk))
         call unpack_from_slots(a12_2_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), a12_2(:,:,iblk))
         call unpack_from_slots(a12_3_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), a12_3(:,:,iblk))
         call unpack_from_slots(a12_4_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), a12_4(:,:,iblk))
      enddo
      !$OMP END PARALLEL DO

      if (maskhalo_dyn) call ice_HaloDestroy(halo_info_mask)

      !-----------------------------------------------------------------
      ! ice-ocean stress
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
      do iblk = 1, nblocks

         call dyn_finish                               &
              (nx_block,           ny_block,           &
               icellU      (iblk), Cdn_ocnU(:,:,iblk), &
               indxUi    (:,iblk), indxUj    (:,iblk), &
               uvel    (:,:,iblk), vvel    (:,:,iblk), &
               uocnU   (:,:,iblk), vocnU   (:,:,iblk), &
               aiU     (:,iblk),   fmU     (:,iblk), &
               strocnxU(:,:,iblk), strocnyU(:,:,iblk))

      enddo
      !$OMP END PARALLEL DO

      call ice_timer_stop(timer_dynamics)    ! dynamics

      end subroutine eap

!=======================================================================
! Initialize parameters and variables needed for the eap dynamics
! (based on init_dyn)

      subroutine init_eap

      use ice_blocks, only: nx_block, ny_block
      use ice_domain, only: nblocks
      use ice_calendar, only: dt_dyn
      use ice_dyn_shared, only: init_dyn_shared, dynTot

      ! local variables

      integer (kind=int_kind) :: &
         i, j, &
         iblk          ! block index

      real (kind=dbl_kind), parameter :: &
         eps6 = 1.0e-6_dbl_kind

      integer (kind=int_kind) :: &
         ix, iy, iz, ia, ierr

      integer (kind=int_kind), parameter :: &
         nz = 100

      real (kind=dbl_kind) :: &
         ainit, xinit, yinit, zinit, &
         da, dx, dy, dz, &
         phi

      real (kind=dbl_kind) :: invstressconviso

      character(len=*), parameter :: subname = '(init_eap)'

      call icepack_query_parameters(puny_out=puny, &
         pi_out=pi, pi2_out=pi2, piq_out=piq, pih_out=pih)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      phi = pi/c12 ! diamond shaped floe smaller angle (default phi = 30 deg)

      call init_dyn_shared(dt_dyn)

      allocate( a11_1        (nx_block,ny_block,max_blocks), &
                a11_2        (nx_block,ny_block,max_blocks), &
                a11_3        (nx_block,ny_block,max_blocks), &
                a11_4        (nx_block,ny_block,max_blocks), &
                a12_1        (nx_block,ny_block,max_blocks), &
                a12_2        (nx_block,ny_block,max_blocks), &
                a12_3        (nx_block,ny_block,max_blocks), &
                a12_4        (nx_block,ny_block,max_blocks), &
                e11          (nx_block,ny_block,max_blocks), &
                e12          (nx_block,ny_block,max_blocks), &
                e22          (nx_block,ny_block,max_blocks), &
                yieldstress11(nx_block,ny_block,max_blocks), &
                yieldstress12(nx_block,ny_block,max_blocks), &
                yieldstress22(nx_block,ny_block,max_blocks), &
                s11          (nx_block,ny_block,max_blocks), &
                s12          (nx_block,ny_block,max_blocks), &
                s22          (nx_block,ny_block,max_blocks), &
                a11          (nx_block,ny_block,max_blocks), &
                a12          (nx_block,ny_block,max_blocks), &
                stat=ierr)
      if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory')

      ! init_dyn_shared above has built the compressed layout, so dynTot is set
      allocate( a11_1_slot(dynTot), a11_2_slot(dynTot),      &
                a11_3_slot(dynTot), a11_4_slot(dynTot),      &
                a12_1_slot(dynTot), a12_2_slot(dynTot),      &
                a12_3_slot(dynTot), a12_4_slot(dynTot),      &
                stressp_1_slot (dynTot), stressp_2_slot (dynTot), &
                stressp_3_slot (dynTot), stressp_4_slot (dynTot), &
                stressm_1_slot (dynTot), stressm_2_slot (dynTot), &
                stressm_3_slot (dynTot), stressm_4_slot (dynTot), &
                stress12_1_slot(dynTot), stress12_2_slot(dynTot), &
                stress12_3_slot(dynTot), stress12_4_slot(dynTot), &
                stat=ierr)
      if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory')
      a11_1_slot = c0 ; a11_2_slot = c0 ; a11_3_slot = c0 ; a11_4_slot = c0
      a12_1_slot = c0 ; a12_2_slot = c0 ; a12_3_slot = c0 ; a12_4_slot = c0
      stressp_1_slot  = c0 ; stressp_2_slot  = c0
      stressp_3_slot  = c0 ; stressp_4_slot  = c0
      stressm_1_slot  = c0 ; stressm_2_slot  = c0
      stressm_3_slot  = c0 ; stressm_4_slot  = c0
      stress12_1_slot = c0 ; stress12_2_slot = c0
      stress12_3_slot = c0 ; stress12_4_slot = c0


      !$OMP PARALLEL DO PRIVATE(iblk,i,j) SCHEDULE(runtime)
      do iblk = 1, nblocks
      do j = 1, ny_block
      do i = 1, nx_block
         e11          (i,j,iblk) = c0
         e12          (i,j,iblk) = c0
         e22          (i,j,iblk) = c0
         s11          (i,j,iblk) = c0
         s12          (i,j,iblk) = c0
         s22          (i,j,iblk) = c0
         yieldstress11(i,j,iblk) = c0
         yieldstress12(i,j,iblk) = c0
         yieldstress22(i,j,iblk) = c0
         a11_1        (i,j,iblk) = p5
         a11_2        (i,j,iblk) = p5
         a11_3        (i,j,iblk) = p5
         a11_4        (i,j,iblk) = p5
         a12_1        (i,j,iblk) = c0
         a12_2        (i,j,iblk) = c0
         a12_3        (i,j,iblk) = c0
         a12_4        (i,j,iblk) = c0
         rdg_shear    (i,j,iblk) = c0
      enddo                     ! i
      enddo                     ! j
      enddo                     ! iblk
      !$OMP END PARALLEL DO

      !-----------------------------------------------------------------
      ! create lookup table for eap dynamics (see Appendix A1)
      !-----------------------------------------------------------------

      da = p5/real(na_yield-1,kind=dbl_kind)
      ainit = p5 - da
      dx = pi/real(nx_yield-1,kind=dbl_kind)
      xinit = pi + piq - dx
      dz = pi/real(nz,kind=dbl_kind)
      zinit = -pih
      dy = pi/real(ny_yield-1,kind=dbl_kind)
      yinit = -dy
      invdx = c1/dx
      invdy = c1/dy
      invda = c1/da

      do ia=1,na_yield
      do ix=1,nx_yield
      do iy=1,ny_yield
         s11r(ix,iy,ia) = c0
         s12r(ix,iy,ia) = c0
         s22r(ix,iy,ia) = c0
         s11s(ix,iy,ia) = c0
         s12s(ix,iy,ia) = c0
         s22s(ix,iy,ia) = c0
         if (ia <= na_yield-1) then
             do iz=1,nz
                s11r(ix,iy,ia) = s11r(ix,iy,ia) + 1*w1(ainit+ia*da)* &
                   exp(-w2(ainit+ia*da)*(zinit+iz*dz)*(zinit+iz*dz))* &
                s11kr(xinit+ix*dx,yinit+iy*dy,zinit+iz*dz,phi)*dz/sin(c2*phi)
                s12r(ix,iy,ia) = s12r(ix,iy,ia) + 1*w1(ainit+ia*da)* &
                   exp(-w2(ainit+ia*da)*(zinit+iz*dz)*(zinit+iz*dz))* &
                s12kr(xinit+ix*dx,yinit+iy*dy,zinit+iz*dz,phi)*dz/sin(c2*phi)
                s22r(ix,iy,ia) = s22r(ix,iy,ia) + 1*w1(ainit+ia*da)* &
                   exp(-w2(ainit+ia*da)*(zinit+iz*dz)*(zinit+iz*dz))* &
                s22kr(xinit+ix*dx,yinit+iy*dy,zinit+iz*dz,phi)*dz/sin(c2*phi)
                s11s(ix,iy,ia) = s11s(ix,iy,ia) + 1*w1(ainit+ia*da)* &
                   exp(-w2(ainit+ia*da)*(zinit+iz*dz)*(zinit+iz*dz))* &
                s11ks(xinit+ix*dx,yinit+iy*dy,zinit+iz*dz,phi)*dz/sin(c2*phi)
                s12s(ix,iy,ia) = s12s(ix,iy,ia) + 1*w1(ainit+ia*da)* &
                   exp(-w2(ainit+ia*da)*(zinit+iz*dz)*(zinit+iz*dz))* &
                s12ks(xinit+ix*dx,yinit+iy*dy,zinit+iz*dz,phi)*dz/sin(c2*phi)
                s22s(ix,iy,ia) = s22s(ix,iy,ia) + 1*w1(ainit+ia*da)* &
                   exp(-w2(ainit+ia*da)*(zinit+iz*dz)*(zinit+iz*dz))* &
                s22ks(xinit+ix*dx,yinit+iy*dy,zinit+iz*dz,phi)*dz/sin(c2*phi)
             enddo
             if (abs(s11r(ix,iy,ia)) < eps6) s11r(ix,iy,ia) = c0
             if (abs(s12r(ix,iy,ia)) < eps6) s12r(ix,iy,ia) = c0
             if (abs(s22r(ix,iy,ia)) < eps6) s22r(ix,iy,ia) = c0
             if (abs(s11s(ix,iy,ia)) < eps6) s11s(ix,iy,ia) = c0
             if (abs(s12s(ix,iy,ia)) < eps6) s12s(ix,iy,ia) = c0
             if (abs(s22s(ix,iy,ia)) < eps6) s22s(ix,iy,ia) = c0
         else
             s11r(ix,iy,ia) = p5*s11kr(xinit+ix*dx,yinit+iy*dy,c0,phi)/sin(c2*phi)
             s12r(ix,iy,ia) = p5*s12kr(xinit+ix*dx,yinit+iy*dy,c0,phi)/sin(c2*phi)
             s22r(ix,iy,ia) = p5*s22kr(xinit+ix*dx,yinit+iy*dy,c0,phi)/sin(c2*phi)
             s11s(ix,iy,ia) = p5*s11ks(xinit+ix*dx,yinit+iy*dy,c0,phi)/sin(c2*phi)
             s12s(ix,iy,ia) = p5*s12ks(xinit+ix*dx,yinit+iy*dy,c0,phi)/sin(c2*phi)
             s22s(ix,iy,ia) = p5*s22ks(xinit+ix*dx,yinit+iy*dy,c0,phi)/sin(c2*phi)
             if (abs(s11r(ix,iy,ia)) < eps6) s11r(ix,iy,ia) = c0
             if (abs(s12r(ix,iy,ia)) < eps6) s12r(ix,iy,ia) = c0
             if (abs(s22r(ix,iy,ia)) < eps6) s22r(ix,iy,ia) = c0
             if (abs(s11s(ix,iy,ia)) < eps6) s11s(ix,iy,ia) = c0
             if (abs(s12s(ix,iy,ia)) < eps6) s12s(ix,iy,ia) = c0
             if (abs(s22s(ix,iy,ia)) < eps6) s22s(ix,iy,ia) = c0
         endif
      enddo
      enddo
      enddo

      ! Factor to maintain the same stress as in EVP (see Section 3)
      ! Can be set to 1 otherwise

      invstressconviso = c1/(c1+kfriction*kfriction)
      invsin = c1/sin(pi2/c12) * invstressconviso

      end subroutine init_eap

!=======================================================================
! Function : w1 (see Gaussian function psi in Tsamados et al 2013)

      FUNCTION w1(a)
      double precision, intent(in) :: a

      real (kind=dbl_kind) :: w1
      character(len=*), parameter :: subname = '(w1)'

      w1 = - 223.87569446_dbl_kind &
           + 2361.2198663_dbl_kind*a &
           - 10606.56079975_dbl_kind*a*a &
           + 26315.50025642_dbl_kind*a*a*a &
           - 38948.30444297_dbl_kind*a*a*a*a &
           + 34397.72407466_dbl_kind*a*a*a*a*a &
           - 16789.98003081_dbl_kind*a*a*a*a*a*a &
           + 3495.82839237_dbl_kind*a*a*a*a*a*a*a

      end FUNCTION w1

!=======================================================================
! Function : w2 (see Gaussian function psi in Tsamados et al 2013)

      FUNCTION w2(a)
      double precision, intent(in) :: a

      real (kind=dbl_kind) :: w2
      character(len=*), parameter :: subname = '(w2)'

      w2 = - 6670.68911883_dbl_kind &
           + 70222.33061536_dbl_kind*a &
           - 314871.71525448_dbl_kind*a*a &
           + 779570.02793492_dbl_kind*a*a*a &
           - 1151098.82436864_dbl_kind*a*a*a*a &
           + 1013896.59464498_dbl_kind*a*a*a*a*a &
           - 493379.44906738_dbl_kind*a*a*a*a*a*a &
           + 102356.551518_dbl_kind*a*a*a*a*a*a*a

      end FUNCTION w2

!=======================================================================
! Function : s11kr

      FUNCTION s11kr(x,y,z,phi)

      real (kind=dbl_kind), intent(in) :: &
        x,y,z,phi

      real (kind=dbl_kind) :: &
         s11kr, p

      real (kind=dbl_kind) :: &
         n1t2i11, n1t2i12, n1t2i21, n1t2i22, &
         n2t1i11, n2t1i12, n2t1i21, n2t1i22, &
!        t1t2i11, t1t2i12, t1t2i21, t1t2i22, &
!        t2t1i11, t2t1i12, t2t1i21, t2t1i22, &
         d11, d12, d22, &
         IIn1t2, IIn2t1, &
!        IIt1t2, &
         Hen1t2, Hen2t1

      character(len=*), parameter :: subname = '(s11kr)'

      p = phi

      n1t2i11 = cos(z+pih-p) * cos(z+p)
      n1t2i12 = cos(z+pih-p) * sin(z+p)
      n1t2i21 = sin(z+pih-p) * cos(z+p)
      n1t2i22 = sin(z+pih-p) * sin(z+p)
      n2t1i11 = cos(z-pih+p) * cos(z-p)
      n2t1i12 = cos(z-pih+p) * sin(z-p)
      n2t1i21 = sin(z-pih+p) * cos(z-p)
      n2t1i22 = sin(z-pih+p) * sin(z-p)
!     t1t2i11 = cos(z-p) * cos(z+p)
!     t1t2i12 = cos(z-p) * sin(z+p)
!     t1t2i21 = sin(z-p) * cos(z+p)
!     t1t2i22 = sin(z-p) * sin(z+p)
!     t2t1i11 = cos(z+p) * cos(z-p)
!     t2t1i12 = cos(z+p) * sin(z-p)
!     t2t1i21 = sin(z+p) * cos(z-p)
!     t2t1i22 = sin(z+p) * sin(z-p)
! In expression of tensor d, with this formulatin d(x)=-d(x+pi)
! Solution, when diagonalizing always check sgn(a11-a22) if > then keep x else x=x-pi/2
      d11 = cos(y)*cos(y)*(cos(x)+sin(x)*tan(y)*tan(y))
      d12 = cos(y)*cos(y)*tan(y)*(-cos(x)+sin(x))
      d22 = cos(y)*cos(y)*(sin(x)+cos(x)*tan(y)*tan(y))
      IIn1t2 = n1t2i11 * d11 + (n1t2i12 + n1t2i21) * d12 + n1t2i22 * d22
      IIn2t1 = n2t1i11 * d11 + (n2t1i12 + n2t1i21) * d12 + n2t1i22 * d22
!     IIt1t2 = t1t2i11 * d11 + (t1t2i12 + t1t2i21) * d12 + t1t2i22 * d22

      if (-IIn1t2>=puny) then
         Hen1t2 = c1
      else
         Hen1t2 = c0
      endif

      if (-IIn2t1>=puny) then
         Hen2t1 = c1
      else
         Hen2t1 = c0
      endif

      s11kr = (- Hen1t2 * n1t2i11 - Hen2t1 * n2t1i11)

      end FUNCTION s11kr

!=======================================================================
! Function : s12kr

      FUNCTION s12kr(x,y,z,phi)

      real (kind=dbl_kind), intent(in) :: &
         x,y,z,phi

      real (kind=dbl_kind) :: &
         s12kr, s12r0, s21r0, p

      real (kind=dbl_kind) :: &
         n1t2i11, n1t2i12, n1t2i21, n1t2i22, &
         n2t1i11, n2t1i12, n2t1i21, n2t1i22, &
!        t1t2i11, t1t2i12, t1t2i21, t1t2i22, &
!        t2t1i11, t2t1i12, t2t1i21, t2t1i22, &
         d11, d12, d22, &
         IIn1t2, IIn2t1, &
!        IIt1t2, &
         Hen1t2, Hen2t1

      character(len=*), parameter :: subname = '(s12kr)'

      p = phi

      n1t2i11 = cos(z+pih-p) * cos(z+p)
      n1t2i12 = cos(z+pih-p) * sin(z+p)
      n1t2i21 = sin(z+pih-p) * cos(z+p)
      n1t2i22 = sin(z+pih-p) * sin(z+p)
      n2t1i11 = cos(z-pih+p) * cos(z-p)
      n2t1i12 = cos(z-pih+p) * sin(z-p)
      n2t1i21 = sin(z-pih+p) * cos(z-p)
      n2t1i22 = sin(z-pih+p) * sin(z-p)
!     t1t2i11 = cos(z-p) * cos(z+p)
!     t1t2i12 = cos(z-p) * sin(z+p)
!     t1t2i21 = sin(z-p) * cos(z+p)
!     t1t2i22 = sin(z-p) * sin(z+p)
!     t2t1i11 = cos(z+p) * cos(z-p)
!     t2t1i12 = cos(z+p) * sin(z-p)
!     t2t1i21 = sin(z+p) * cos(z-p)
!     t2t1i22 = sin(z+p) * sin(z-p)
      d11 = cos(y)*cos(y)*(cos(x)+sin(x)*tan(y)*tan(y))
      d12 = cos(y)*cos(y)*tan(y)*(-cos(x)+sin(x))
      d22 = cos(y)*cos(y)*(sin(x)+cos(x)*tan(y)*tan(y))
      IIn1t2 = n1t2i11 * d11 + (n1t2i12 + n1t2i21) * d12 + n1t2i22 * d22
      IIn2t1 = n2t1i11 * d11 + (n2t1i12 + n2t1i21) * d12 + n2t1i22 * d22
!     IIt1t2 = t1t2i11 * d11 + (t1t2i12 + t1t2i21) * d12 + t1t2i22 * d22

      if (-IIn1t2>=puny) then
         Hen1t2 = c1
      else
         Hen1t2 = c0
      endif

      if (-IIn2t1>=puny) then
         Hen2t1 = c1
      else
         Hen2t1 = c0
      endif

      s12r0 = (- Hen1t2 * n1t2i12 - Hen2t1 * n2t1i12)
      s21r0 = (- Hen1t2 * n1t2i21 - Hen2t1 * n2t1i21)
      s12kr=p5*(s12r0+s21r0)

      end FUNCTION s12kr

!=======================================================================
! Function : s22r

      FUNCTION s22kr(x,y,z,phi)

      real (kind=dbl_kind), intent(in) :: &
         x,y,z,phi

      real (kind=dbl_kind) :: &
         s22kr, p

      real (kind=dbl_kind) :: &
         n1t2i11, n1t2i12, n1t2i21, n1t2i22, &
         n2t1i11, n2t1i12, n2t1i21, n2t1i22, &
!        t1t2i11, t1t2i12, t1t2i21, t1t2i22, &
!        t2t1i11, t2t1i12, t2t1i21, t2t1i22, &
         d11, d12, d22, &
         IIn1t2, IIn2t1, &
!        IIt1t2, &
         Hen1t2, Hen2t1

      character(len=*), parameter :: subname = '(s22kr)'

      p = phi

      n1t2i11 = cos(z+pih-p) * cos(z+p)
      n1t2i12 = cos(z+pih-p) * sin(z+p)
      n1t2i21 = sin(z+pih-p) * cos(z+p)
      n1t2i22 = sin(z+pih-p) * sin(z+p)
      n2t1i11 = cos(z-pih+p) * cos(z-p)
      n2t1i12 = cos(z-pih+p) * sin(z-p)
      n2t1i21 = sin(z-pih+p) * cos(z-p)
      n2t1i22 = sin(z-pih+p) * sin(z-p)
!     t1t2i11 = cos(z-p) * cos(z+p)
!     t1t2i12 = cos(z-p) * sin(z+p)
!     t1t2i21 = sin(z-p) * cos(z+p)
!     t1t2i22 = sin(z-p) * sin(z+p)
!     t2t1i11 = cos(z+p) * cos(z-p)
!     t2t1i12 = cos(z+p) * sin(z-p)
!     t2t1i21 = sin(z+p) * cos(z-p)
!     t2t1i22 = sin(z+p) * sin(z-p)
      d11 = cos(y)*cos(y)*(cos(x)+sin(x)*tan(y)*tan(y))
      d12 = cos(y)*cos(y)*tan(y)*(-cos(x)+sin(x))
      d22 = cos(y)*cos(y)*(sin(x)+cos(x)*tan(y)*tan(y))
      IIn1t2 = n1t2i11 * d11 + (n1t2i12 + n1t2i21) * d12 + n1t2i22 * d22
      IIn2t1 = n2t1i11 * d11 + (n2t1i12 + n2t1i21) * d12 + n2t1i22 * d22
!     IIt1t2 = t1t2i11 * d11 + (t1t2i12 + t1t2i21) * d12 + t1t2i22 * d22

      if (-IIn1t2>=puny) then
         Hen1t2 = c1
      else
         Hen1t2 = c0
      endif

      if (-IIn2t1>=puny) then
         Hen2t1 = c1
      else
         Hen2t1 = c0
      endif

      s22kr = (- Hen1t2 * n1t2i22 - Hen2t1 * n2t1i22)

      end FUNCTION s22kr

!=======================================================================
! Function : s11ks

      FUNCTION s11ks(x,y,z,phi)

      real (kind=dbl_kind), intent(in):: &
         x,y,z,phi

      real (kind=dbl_kind) :: &
         s11ks, p

      real (kind=dbl_kind) :: &
         n1t2i11, n1t2i12, n1t2i21, n1t2i22, &
         n2t1i11, n2t1i12, n2t1i21, n2t1i22, &
         t1t2i11, &
         t1t2i12, t1t2i21, t1t2i22, &
         t2t1i11, &
!        t2t1i12, t2t1i21, t2t1i22, &
         d11, d12, d22, &
         IIn1t2, IIn2t1, IIt1t2, &
         Hen1t2, Hen2t1

      character(len=*), parameter :: subname = '(s11ks)'

      p = phi

      n1t2i11 = cos(z+pih-p) * cos(z+p)
      n1t2i12 = cos(z+pih-p) * sin(z+p)
      n1t2i21 = sin(z+pih-p) * cos(z+p)
      n1t2i22 = sin(z+pih-p) * sin(z+p)
      n2t1i11 = cos(z-pih+p) * cos(z-p)
      n2t1i12 = cos(z-pih+p) * sin(z-p)
      n2t1i21 = sin(z-pih+p) * cos(z-p)
      n2t1i22 = sin(z-pih+p) * sin(z-p)
      t1t2i11 = cos(z-p) * cos(z+p)
      t1t2i12 = cos(z-p) * sin(z+p)
      t1t2i21 = sin(z-p) * cos(z+p)
      t1t2i22 = sin(z-p) * sin(z+p)
      t2t1i11 = cos(z+p) * cos(z-p)
!     t2t1i12 = cos(z+p) * sin(z-p)
!     t2t1i21 = sin(z+p) * cos(z-p)
!     t2t1i22 = sin(z+p) * sin(z-p)
      d11 = cos(y)*cos(y)*(cos(x)+sin(x)*tan(y)*tan(y))
      d12 = cos(y)*cos(y)*tan(y)*(-cos(x)+sin(x))
      d22 = cos(y)*cos(y)*(sin(x)+cos(x)*tan(y)*tan(y))
      IIn1t2 = n1t2i11 * d11 + (n1t2i12 + n1t2i21) * d12 + n1t2i22 * d22
      IIn2t1 = n2t1i11 * d11 + (n2t1i12 + n2t1i21) * d12 + n2t1i22 * d22
      IIt1t2 = t1t2i11 * d11 + (t1t2i12 + t1t2i21) * d12 + t1t2i22 * d22

      if (-IIn1t2>=puny) then
         Hen1t2 = c1
      else
         Hen1t2 = c0
      endif

      if (-IIn2t1>=puny) then
         Hen2t1 = c1
      else
         Hen2t1 = c0
      endif

      s11ks = sign(c1,IIt1t2+puny)*(Hen1t2 * t1t2i11 + Hen2t1 * t2t1i11)

      end FUNCTION s11ks

!=======================================================================
! Function : s12ks

      FUNCTION s12ks(x,y,z,phi)

      real (kind=dbl_kind), intent(in) :: &
         x,y,z,phi

      real (kind=dbl_kind) :: &
         s12ks,s12s0,s21s0,p

      real (kind=dbl_kind) :: &
         n1t2i11, n1t2i12, n1t2i21, n1t2i22, &
         n2t1i11, n2t1i12, n2t1i21, n2t1i22, &
         t1t2i11, t1t2i12, t1t2i21, t1t2i22, &
!        t2t1i11, t2t1i22, &
         t2t1i12, t2t1i21, &
         d11, d12, d22, &
         IIn1t2, IIn2t1, IIt1t2, &
         Hen1t2, Hen2t1

      character(len=*), parameter :: subname = '(s12ks)'

      p =phi

      n1t2i11 = cos(z+pih-p) * cos(z+p)
      n1t2i12 = cos(z+pih-p) * sin(z+p)
      n1t2i21 = sin(z+pih-p) * cos(z+p)
      n1t2i22 = sin(z+pih-p) * sin(z+p)
      n2t1i11 = cos(z-pih+p) * cos(z-p)
      n2t1i12 = cos(z-pih+p) * sin(z-p)
      n2t1i21 = sin(z-pih+p) * cos(z-p)
      n2t1i22 = sin(z-pih+p) * sin(z-p)
      t1t2i11 = cos(z-p) * cos(z+p)
      t1t2i12 = cos(z-p) * sin(z+p)
      t1t2i21 = sin(z-p) * cos(z+p)
      t1t2i22 = sin(z-p) * sin(z+p)
!     t2t1i11 = cos(z+p) * cos(z-p)
      t2t1i12 = cos(z+p) * sin(z-p)
      t2t1i21 = sin(z+p) * cos(z-p)
!     t2t1i22 = sin(z+p) * sin(z-p)
      d11 = cos(y)*cos(y)*(cos(x)+sin(x)*tan(y)*tan(y))
      d12 = cos(y)*cos(y)*tan(y)*(-cos(x)+sin(x))
      d22 = cos(y)*cos(y)*(sin(x)+cos(x)*tan(y)*tan(y))
      IIn1t2 = n1t2i11 * d11 + (n1t2i12 + n1t2i21) * d12 + n1t2i22 * d22
      IIn2t1 = n2t1i11 * d11 + (n2t1i12 + n2t1i21) * d12 + n2t1i22 * d22
      IIt1t2 = t1t2i11 * d11 + (t1t2i12 + t1t2i21) * d12 + t1t2i22 * d22

      if (-IIn1t2>=puny) then
         Hen1t2 = c1
      else
         Hen1t2 = c0
      endif

      if (-IIn2t1>=puny) then
         Hen2t1 = c1
      else
         Hen2t1 = c0
      endif

      s12s0 = sign(c1,IIt1t2+puny)*(Hen1t2 * t1t2i12 + Hen2t1 * t2t1i12)
      s21s0 = sign(c1,IIt1t2+puny)*(Hen1t2 * t1t2i21 + Hen2t1 * t2t1i21)
      s12ks=p5*(s12s0+s21s0)

      end FUNCTION s12ks

!=======================================================================
! Function : s22ks

      FUNCTION s22ks(x,y,z,phi)

      real (kind=dbl_kind), intent(in) :: &
         x,y,z,phi

      real (kind=dbl_kind) :: &
         s22ks,p

      real (kind=dbl_kind) :: &
         n1t2i11, n1t2i12, n1t2i21, n1t2i22, &
         n2t1i11, n2t1i12, n2t1i21, n2t1i22, &
         t1t2i11, t1t2i12, t1t2i21, t1t2i22, &
!        t2t1i11, t2t1i12, t2t1i21, &
         t2t1i22, &
         d11, d12, d22, &
         IIn1t2, IIn2t1, IIt1t2, &
         Hen1t2, Hen2t1

      character(len=*), parameter :: subname = '(s22ks)'

      p = phi

      n1t2i11 = cos(z+pih-p) * cos(z+p)
      n1t2i12 = cos(z+pih-p) * sin(z+p)
      n1t2i21 = sin(z+pih-p) * cos(z+p)
      n1t2i22 = sin(z+pih-p) * sin(z+p)
      n2t1i11 = cos(z-pih+p) * cos(z-p)
      n2t1i12 = cos(z-pih+p) * sin(z-p)
      n2t1i21 = sin(z-pih+p) * cos(z-p)
      n2t1i22 = sin(z-pih+p) * sin(z-p)
      t1t2i11 = cos(z-p) * cos(z+p)
      t1t2i12 = cos(z-p) * sin(z+p)
      t1t2i21 = sin(z-p) * cos(z+p)
      t1t2i22 = sin(z-p) * sin(z+p)
!     t2t1i11 = cos(z+p) * cos(z-p)
!     t2t1i12 = cos(z+p) * sin(z-p)
!     t2t1i21 = sin(z+p) * cos(z-p)
      t2t1i22 = sin(z+p) * sin(z-p)
      d11 = cos(y)*cos(y)*(cos(x)+sin(x)*tan(y)*tan(y))
      d12 = cos(y)*cos(y)*tan(y)*(-cos(x)+sin(x))
      d22 = cos(y)*cos(y)*(sin(x)+cos(x)*tan(y)*tan(y))
      IIn1t2 = n1t2i11 * d11 + (n1t2i12 + n1t2i21) * d12 + n1t2i22 * d22
      IIn2t1 = n2t1i11 * d11 + (n2t1i12 + n2t1i21) * d12 + n2t1i22 * d22
      IIt1t2 = t1t2i11 * d11 + (t1t2i12 + t1t2i21) * d12 + t1t2i22 * d22

      if (-IIn1t2>=puny) then
         Hen1t2 = c1
      else
         Hen1t2 = c0
      endif

      if (-IIn2t1>=puny) then
         Hen2t1 = c1
      else
         Hen2t1 = c0
      endif

      s22ks = sign(c1,IIt1t2+puny)*(Hen1t2 * t1t2i22 + Hen2t1 * t2t1i22)

      end FUNCTION s22ks

!=======================================================================
! Computes the rates of strain and internal stress components for
! each of the four corners on each T-grid cell.
! Computes stress terms for the momentum equation
! (based on subroutine stress)

      subroutine stress_eap  (nx_block,   ny_block,       &
                              ksub,       ndte,           &
                              icellT,                     &
                              indxTi,     indxTj,         &
                              arlx1i,     denom1,         &
                              uvel,       vvel,     nbr,  &
                              dxT,        dyT,            &
                              dxU,        dyU,            &
                              dxhy,       dyhx,           &
                              cxp,        cyp,            &
                              cxm,        cym,            &
                              tarear,     strength,       &
                              a11_1, a11_2, a11_3, a11_4, &
                              a12_1, a12_2, a12_3, a12_4, &
                              stressp_1,  stressp_2,      &
                              stressp_3,  stressp_4,      &
                              stressm_1,  stressm_2,      &
                              stressm_3,  stressm_4,      &
                              stress12_1, stress12_2,     &
                              stress12_3, stress12_4,     &
                              shear,      divu,           &
                              vort,                       &
                              e11,        e12,            &
                              e22,                        &
                              s11,        s12,            &
                              s22,                        &
                              yieldstress11,              &
                              yieldstress12,              &
                              yieldstress22,              &
!                             rdg_conv,   rdg_shear,      &
                              rdg_conv, &
                              Xmap,     strtmp)

      use ice_dyn_shared, only: nbr_W, nbr_S, nbr_SW

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ksub              , & ! subcycling step
         ndte              , & ! number of subcycles
         icellT                ! no. of cells where iceTmask = .true.

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxTi   , & ! compressed index in i-direction
         indxTj       ! compressed index in j-direction

      real (kind=dbl_kind), intent(in) :: &
         arlx1i   , & ! dte/2T (original) or 1/alpha1 (revised)
         denom1       ! constant for stress equation

      ! uvel and vvel are COMPRESSED: this block's slice, 1-based, so cells
      ! are 1..cnt and the last element is the block's dead slot.  They must
      ! be the compressed
      ! velocities, not the block-layout ones: stepu subcycles on the slots
      ! and the block arrays are stale until the subcycling ends.
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         uvel     , & ! x-component of velocity (m/s)
         vvel         ! y-component of velocity (m/s)

      ! nbr(is,d) is the block-local slot one step in direction d from slot is
      integer (kind=int_kind), dimension(:,:), intent(in) :: &
         nbr          ! neighbour slots

      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(in) :: &
         strength , & ! ice strength (N/m)
         dxT      , & ! width of T-cell through the middle (m)
         dyT      , & ! height of T-cell through the middle (m)
         dxU      , & ! width of U-cell through the middle (m)
         dyU      , & ! height of U-cell through the middle (m)
         dxhy     , & ! 0.5*(HTE - HTW)
         dyhx     , & ! 0.5*(HTN - HTS)
         cyp      , & ! 1.5*HTE - 0.5*HTW
         cxp      , & ! 1.5*HTN - 0.5*HTS
         cym      , & ! 0.5*HTE - 1.5*HTW
         cxm          ! 0.5*HTN - 1.5*HTS

      ! tarear is assumed shape, 1D over the block in Fortran storage order,
      ! i.e. element (i,j) is stored at (j-1)*nx_block+i
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         tarear       ! 1/tarea

      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(inout) :: &
         stressp_1, stressp_2, stressp_3, stressp_4, & ! sigma11+sigma22
         stressm_1, stressm_2, stressm_3, stressm_4, & ! sigma11-sigma22
         stress12_1,stress12_2,stress12_3,stress12_4   ! sigma12

      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(inout) :: &
         a11_1, a11_2, a11_3, a11_4, & ! structure tensor
         a12_1, a12_2, a12_3, a12_4              ! structure tensor

      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(inout) :: &
         shear    , & ! strain rate II component (1/s)
         divu     , & ! strain rate I component, velocity divergence (1/s)
         vort     , & ! vorticity (1/s)
         e11      , & ! components of strain rate tensor (1/s)
         e12      , & !
         e22      , & !
         s11      , & ! components of stress tensor (kg/s^2)
         s12      , & !
         s22      , & !
         yieldstress11, & ! components of yield stress tensor (kg/s^2)
         yieldstress12, &
         yieldstress22, &
         rdg_conv     ! convergence term for ridging (1/s)
!        rdg_shear    ! shear term for ridging (1/s)

      ! COMPRESSED, indexed by block-local slot
      real (kind=dbl_kind), dimension (:,:), intent(out), contiguous :: &
         strtmp       ! stress combinations

      ! Xmap holds each cell's block-local slot, see ice_dyn_shared
      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         Xmap

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij, ijb, &
         isW, isS, isSW, & ! W, S and SW neighbour slots
         is             ! block-local compressed slot of (i,j)

      real (kind=dbl_kind) :: &
         stressptmp_1, stressptmp_2, stressptmp_3, stressptmp_4, & ! sigma11+sigma22
         stressmtmp_1, stressmtmp_2, stressmtmp_3, stressmtmp_4, & ! sigma11-sigma22
         stress12tmp_1,stress12tmp_2,stress12tmp_3,stress12tmp_4   ! sigma12

      real (kind=dbl_kind) :: &
         divune, divunw, divuse, divusw            , & ! divergence
         tensionne, tensionnw, tensionse, tensionsw, & ! tension
         shearne, shearnw, shearse, shearsw        , & ! shearing
         u00, uW, uS, uSW                          , & ! velocity stencil,
         v00, vW, vS, vSW                          , & ! self and W, S, SW
         dvdxn, dvdxs, dudye, dudyw                , & ! for vorticity calc
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
         alpharne, alpharnw, alpharsw, alpharse,     &
         alphasne, alphasnw, alphassw, alphasse

      character(len=*), parameter :: subname = '(stress_eap)'

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      strtmp(:,:) = c0

      do ij = 1, icellT
         i = indxTi(ij)
         j = indxTj(ij)
         ijb = (j-1)*nx_block + i
         is  = Xmap(ijb)

         ! the velocity stencil, read from the compressed slots.  A missing
         ! neighbour is the block's dead slot, the last one, held at zero --
         ! which is what the block layout carried at land.
         ! each neighbour slot was fetched twice, once for uvel and once for
         ! vvel; read them once instead
         isW  = nbr(is,nbr_W )
         isS  = nbr(is,nbr_S )
         isSW = nbr(is,nbr_SW)

         u00 = uvel(is)
         uW  = uvel(isW)
         uS  = uvel(isS)
         uSW = uvel(isSW)
         v00 = vvel(is)
         vW  = vvel(isW)
         vS  = vvel(isS)
         vSW = vvel(isSW)

         !-----------------------------------------------------------------
         ! strain rates
         ! NOTE these are actually strain rates * area  (m^2/s)
         !-----------------------------------------------------------------

         ! divergence  =  e_11 + e_22
         divune    = cyp(ijb)*u00 - dyT(ijb)*uW &
                   + cxp(ijb)*v00 - dxT(ijb)*vS
         divunw    = cym(ijb)*uW + dyT(ijb)*u00 &
                   + cxp(ijb)*vW - dxT(ijb)*vSW
         divusw    = cym(ijb)*uSW + dyT(ijb)*uS &
                   + cxm(ijb)*vSW + dxT(ijb)*vW
         divuse    = cyp(ijb)*uS - dyT(ijb)*uSW &
                   + cxm(ijb)*vS + dxT(ijb)*v00

         ! tension strain rate  =  e_11 - e_22
         tensionne = -cym(ijb)*u00 - dyT(ijb)*uW &
                   +  cxm(ijb)*v00 + dxT(ijb)*vS
         tensionnw = -cyp(ijb)*uW + dyT(ijb)*u00 &
                   +  cxm(ijb)*vW + dxT(ijb)*vSW
         tensionsw = -cyp(ijb)*uSW + dyT(ijb)*uS &
                   +  cxp(ijb)*vSW - dxT(ijb)*vW
         tensionse = -cym(ijb)*uS - dyT(ijb)*uSW &
                   +  cxp(ijb)*vS - dxT(ijb)*v00

         ! shearing strain rate  =  2*e_12
         shearne = -cym(ijb)*v00 - dyT(ijb)*vW &
                 -  cxm(ijb)*u00 - dxT(ijb)*uS
         shearnw = -cyp(ijb)*vW + dyT(ijb)*v00 &
                 -  cxm(ijb)*uW - dxT(ijb)*uSW
         shearsw = -cyp(ijb)*vSW + dyT(ijb)*vS &
                 -  cxp(ijb)*uSW + dxT(ijb)*uW
         shearse = -cym(ijb)*vS - dyT(ijb)*vSW &
                 -  cxp(ijb)*uS + dxT(ijb)*u00

         !-----------------------------------------------------------------
         ! Stress updated depending on strain rate and structure tensor
         !-----------------------------------------------------------------

         ! ne
         call update_stress_rdg (ksub, ndte, divune, tensionne, &
                                 shearne, a11_1(is), a12_1(is), &
                                 stressptmp_1, stressmtmp_1, &
                                 stress12tmp_1, strength(ijb), &
                                 alpharne, alphasne)
         ! nw
         call update_stress_rdg (ksub, ndte, divunw, tensionnw, &
                                 shearnw, a11_2(is), a12_2(is), &
                                 stressptmp_2, stressmtmp_2, &
                                 stress12tmp_2, strength(ijb), &
                                 alpharnw, alphasnw)
         ! sw
         call update_stress_rdg (ksub, ndte, divusw, tensionsw, &
                                 shearsw, a11_3(is), a12_3(is), &
                                 stressptmp_3, stressmtmp_3, &
                                 stress12tmp_3, strength(ijb), &
                                 alpharsw, alphassw)
         ! se
         call update_stress_rdg (ksub, ndte, divuse, tensionse, &
                                 shearse, a11_4(is), a12_4(is), &
                                 stressptmp_4, stressmtmp_4, &
                                 stress12tmp_4, strength(ijb), &
                                 alpharse, alphasse)

         !-----------------------------------------------------------------
         ! on last subcycle, save quantities for mechanical redistribution
         !-----------------------------------------------------------------

         if (ksub == ndte) then

            ! diagnostic only
            ! shear = sqrt(tension**2 + shearing**2)
            shear(ijb) = p25*tarear(ijb)*sqrt( &
                 (tensionne + tensionnw + tensionse + tensionsw)**2 &
                +  (shearne +   shearnw +   shearse +   shearsw)**2)

            ! vorticity
            dvdxn = dyU(ijb)*v00 - dyU(ijb-1)*vW
            dvdxs = dyU(ijb-nx_block)*vS - dyU(ijb-nx_block-1)*vSW
            dudye = dxU(ijb)*u00 - dxU(ijb-nx_block)*uS
            dudyw = dxU(ijb-1)*uW - dxU(ijb-nx_block-1)*uSW
            vort(ijb) = p5*tarear(ijb)*(dvdxn + dvdxs - dudye - dudyw)

            divu(ijb) = p25*(divune + divunw + divuse + divusw) * tarear(ijb)
            rdg_conv(ijb)  = -min(p25*(alpharne + alpharnw &
                                     + alpharsw + alpharse),c0) * tarear(ijb)
            !rdg_shear=0 for computing closing_net in ridge_prep
            !rdg_shear(i,j) = p25*(alphasne + alphasnw &
            !                    + alphassw + alphasse) * tarear(ijb)
         endif

         e11(ijb) = p5*p25*(divune + divunw + divuse + divusw + &
                    tensionne + tensionnw + tensionse + tensionsw) * tarear(ijb)

         e12(ijb) = p5*p25*(shearne + shearnw + shearse + shearsw) * tarear(ijb)

         e22(ijb) = p5*p25*(divune + divunw + divuse + divusw - &
                    tensionne - tensionnw - tensionse - tensionsw) * tarear(ijb)

         !-----------------------------------------------------------------
         ! elastic relaxation, see Eq. A12-A14
         !-----------------------------------------------------------------

         stressp_1(is) = (stressp_1(is) + stressptmp_1*arlx1i) &
                          * denom1
         stressp_2(is) = (stressp_2(is) + stressptmp_2*arlx1i) &
                          * denom1
         stressp_3(is) = (stressp_3(is) + stressptmp_3*arlx1i) &
                          * denom1
         stressp_4(is) = (stressp_4(is) + stressptmp_4*arlx1i) &
                          * denom1

         stressm_1(is) = (stressm_1(is) + stressmtmp_1*arlx1i) &
                          * denom1
         stressm_2(is) = (stressm_2(is) + stressmtmp_2*arlx1i) &
                          * denom1
         stressm_3(is) = (stressm_3(is) + stressmtmp_3*arlx1i) &
                          * denom1
         stressm_4(is) = (stressm_4(is) + stressmtmp_4*arlx1i) &
                          * denom1

         stress12_1(is) = (stress12_1(is) + stress12tmp_1*arlx1i) &
                          * denom1
         stress12_2(is) = (stress12_2(is) + stress12tmp_2*arlx1i) &
                          * denom1
         stress12_3(is) = (stress12_3(is) + stress12tmp_3*arlx1i) &
                          * denom1
         stress12_4(is) = (stress12_4(is) + stress12tmp_4*arlx1i) &
                          * denom1

         s11(ijb) = p5 * p25 * (stressp_1(is) + stressp_2(is) &
                              + stressp_3(is) + stressp_4(is) &
                              + stressm_1(is) + stressm_2(is) &
                              + stressm_3(is) + stressm_4(is))
         s22(ijb) = p5 * p25 * (stressp_1(is) + stressp_2(is) &
                              + stressp_3(is) + stressp_4(is) &
                              - stressm_1(is) - stressm_2(is) &
                              - stressm_3(is) - stressm_4(is))
         s12(ijb) =      p25 * (stress12_1(is) + stress12_2(is) &
                              + stress12_3(is) + stress12_4(is))

         yieldstress11(ijb) = p5 * p25 * (stressptmp_1 + stressptmp_2 &
                                        + stressptmp_3 + stressptmp_4 &
                                        + stressmtmp_1 + stressmtmp_2 &
                                        + stressmtmp_3 + stressmtmp_4)
         yieldstress22(ijb) = p5 * p25 * (stressptmp_1 + stressptmp_2 &
                                        + stressptmp_3 + stressptmp_4 &
                                        - stressmtmp_1 - stressmtmp_2 &
                                        - stressmtmp_3 - stressmtmp_4)
         yieldstress12(ijb) =      p25 * (stress12tmp_1 + stress12tmp_2 &
                                        + stress12tmp_3 + stress12tmp_4)

         !-----------------------------------------------------------------
         ! Eliminate underflows.
         ! The following code is commented out because it is relatively
         ! expensive and most compilers include a flag that accomplishes
         ! the same thing more efficiently.  This code is cheaper than
         ! handling underflows if the compiler lacks a flag; uncomment
         ! it in that case.  The compiler flag is often described with the
         ! phrase "flush to zero".
         !-----------------------------------------------------------------

!         stressp_1(is) = sign(max(abs(stressp_1(is)),puny),stressp_1(is))
!         stressp_2(is) = sign(max(abs(stressp_2(is)),puny),stressp_2(is))
!         stressp_3(is) = sign(max(abs(stressp_3(is)),puny),stressp_3(is))
!         stressp_4(is) = sign(max(abs(stressp_4(is)),puny),stressp_4(is))

!         stressm_1(is) = sign(max(abs(stressm_1(is)),puny),stressm_1(is))
!         stressm_2(is) = sign(max(abs(stressm_2(is)),puny),stressm_2(is))
!         stressm_3(is) = sign(max(abs(stressm_3(is)),puny),stressm_3(is))
!         stressm_4(is) = sign(max(abs(stressm_4(is)),puny),stressm_4(is))

!         stress12_1(is) = sign(max(abs(stress12_1(is)),puny),stress12_1(is))
!         stress12_2(is) = sign(max(abs(stress12_2(is)),puny),stress12_2(is))
!         stress12_3(is) = sign(max(abs(stress12_3(is)),puny),stress12_3(is))
!         stress12_4(is) = sign(max(abs(stress12_4(is)),puny),stress12_4(is))

         !-----------------------------------------------------------------
         ! combinations of the stresses for the momentum equation ! kg/s^2
         !-----------------------------------------------------------------

         ssigpn  = stressp_1(is) + stressp_2(is)
         ssigps  = stressp_3(is) + stressp_4(is)
         ssigpe  = stressp_1(is) + stressp_4(is)
         ssigpw  = stressp_2(is) + stressp_3(is)
         ssigp1  =(stressp_1(is) + stressp_3(is))*p055
         ssigp2  =(stressp_2(is) + stressp_4(is))*p055

         ssigmn  = stressm_1(is) + stressm_2(is)
         ssigms  = stressm_3(is) + stressm_4(is)
         ssigme  = stressm_1(is) + stressm_4(is)
         ssigmw  = stressm_2(is) + stressm_3(is)
         ssigm1  =(stressm_1(is) + stressm_3(is))*p055
         ssigm2  =(stressm_2(is) + stressm_4(is))*p055

         ssig12n = stress12_1(is) + stress12_2(is)
         ssig12s = stress12_3(is) + stress12_4(is)
         ssig12e = stress12_1(is) + stress12_4(is)
         ssig12w = stress12_2(is) + stress12_3(is)
         ssig121 =(stress12_1(is) + stress12_3(is))*p111
         ssig122 =(stress12_2(is) + stress12_4(is))*p111

         csigpne = p111*stressp_1(is) + ssigp2 + p027*stressp_3(is)
         csigpnw = p111*stressp_2(is) + ssigp1 + p027*stressp_4(is)
         csigpsw = p111*stressp_3(is) + ssigp2 + p027*stressp_1(is)
         csigpse = p111*stressp_4(is) + ssigp1 + p027*stressp_2(is)

         csigmne = p111*stressm_1(is) + ssigm2 + p027*stressm_3(is)
         csigmnw = p111*stressm_2(is) + ssigm1 + p027*stressm_4(is)
         csigmsw = p111*stressm_3(is) + ssigm2 + p027*stressm_1(is)
         csigmse = p111*stressm_4(is) + ssigm1 + p027*stressm_2(is)

         csig12ne = p222*stress12_1(is) + ssig122 &
                  + p055*stress12_3(is)
         csig12nw = p222*stress12_2(is) + ssig121 &
                  + p055*stress12_4(is)
         csig12sw = p222*stress12_3(is) + ssig122 &
                  + p055*stress12_1(is)
         csig12se = p222*stress12_4(is) + ssig121 &
                  + p055*stress12_2(is)

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
         strtmp(is,1) = -strp_tmp - strm_tmp - str12ew &
              + dxhy(ijb)*(-csigpne + csigmne) + dyhx(ijb)*csig12ne

         ! northwest (i+1,j)
         strtmp(is,2) = strp_tmp + strm_tmp - str12we &
              + dxhy(ijb)*(-csigpnw + csigmnw) + dyhx(ijb)*csig12nw

         strp_tmp  = p25*dyT(ijb)*(p333*ssigps  + p166*ssigpn)
         strm_tmp  = p25*dyT(ijb)*(p333*ssigms  + p166*ssigmn)

         ! southeast (i,j+1)
         strtmp(is,3) = -strp_tmp - strm_tmp + str12ew &
              + dxhy(ijb)*(-csigpse + csigmse) + dyhx(ijb)*csig12se

         ! southwest (i+1,j+1)
         strtmp(is,4) = strp_tmp + strm_tmp + str12we &
              + dxhy(ijb)*(-csigpsw + csigmsw) + dyhx(ijb)*csig12sw

         !-----------------------------------------------------------------
         ! for dF/dy (v momentum)
         !-----------------------------------------------------------------

         strp_tmp  = p25*dxT(ijb)*(p333*ssigpe  + p166*ssigpw)
         strm_tmp  = p25*dxT(ijb)*(p333*ssigme  + p166*ssigmw)

         ! northeast (i,j)
         strtmp(is,5) = -strp_tmp + strm_tmp - str12ns &
              - dyhx(ijb)*(csigpne + csigmne) + dxhy(ijb)*csig12ne

         ! southeast (i,j+1)
         strtmp(is,6) = strp_tmp - strm_tmp - str12sn &
              - dyhx(ijb)*(csigpse + csigmse) + dxhy(ijb)*csig12se

         strp_tmp  = p25*dxT(ijb)*(p333*ssigpw  + p166*ssigpe)
         strm_tmp  = p25*dxT(ijb)*(p333*ssigmw  + p166*ssigme)

         ! northwest (i+1,j)
         strtmp(is,7) = -strp_tmp + strm_tmp + str12ns &
              - dyhx(ijb)*(csigpnw + csigmnw) + dxhy(ijb)*csig12nw

         ! southwest (i+1,j+1)
         strtmp(is,8) = strp_tmp - strm_tmp + str12sn &
              - dyhx(ijb)*(csigpsw + csigmsw) + dxhy(ijb)*csig12sw

      enddo                     ! ij

      end subroutine stress_eap

!=======================================================================
! Updates the stress depending on values of strain rate and structure
! tensor and for ksub=ndte it computes closing and sliding rate

      subroutine update_stress_rdg (ksub, ndte, divu, tension, &
                                   shear, a11, a12, &
                                   stressp,  stressm, &
                                   stress12, strength, &
                                   alphar, alphas)

      integer (kind=int_kind), intent(in) :: &
         ksub, &
         ndte

      real (kind=dbl_kind), intent(in) :: &
         a11,     a12,                   &
         divu,    tension, shear,        &
         strength

      real (kind=dbl_kind), intent(out) :: &
         stressp, stressm, stress12, &
         alphar, alphas

      ! local variables

      integer (kind=int_kind) :: &
         kx ,ky, ka

      real (kind=dbl_kind) :: &
         stemp11r, stemp12r, stemp22r,   &
         stemp11s, stemp12s, stemp22s, &
         a22, Qd11Qd11, Qd11Qd12, Qd12Qd12, &
         Q11Q11, Q11Q12, Q12Q12, &
         dtemp11, dtemp12, dtemp22, &
         rotstemp11r, rotstemp12r, rotstemp22r,   &
         rotstemp11s, rotstemp12s, rotstemp22s, &
         sig11, sig12, sig22, &
         sgprm11, sgprm12, sgprm22, &
         Angle_denom_gamma,  Angle_denom_alpha, &
         Tany_1, Tany_2, &
         x, y, dx, &
         dtemp1, dtemp2, atempprime, &
         kxw, kyw, kaw

      ! tcraig, temporary, should be moved to namelist
      ! turns on interpolation in stress_rdg
      logical(kind=log_kind), parameter :: &
         interpolate_stress_rdg = .false.

      character(len=*), parameter :: subname = '(update_stress_rdg)'

      ! compute eigenvalues, eigenvectors and angles for structure tensor, strain rates

      ! 1) structure tensor

      a22 = c1-a11

      ! gamma: angle between general coordinates and principal axis of A
      ! here Tan2gamma = 2 a12 / (a11 - a22)

      Q11Q11 = c1
      Q12Q12 = puny
      Q11Q12 = puny

      if ((ABS(a11 - a22) > puny).or.(ABS(a12) > puny)) then
         Angle_denom_gamma = sqrt( ( a11 - a22 )*( a11 - a22) + &
                                     c4*a12*a12 )

         Q11Q11 = p5 + ( a11 - a22 )*p5/Angle_denom_gamma   !Cos^2
         Q12Q12 = p5 - ( a11 - a22 )*p5/Angle_denom_gamma   !Sin^2
         Q11Q12 = a12/Angle_denom_gamma                             !CosSin
      endif

      ! rotation Q*atemp*Q^T
      atempprime = Q11Q11*a11 + c2*Q11Q12*a12 + Q12Q12*a22

      ! make first principal value the largest
      atempprime = max(atempprime, c1 - atempprime)

      ! 2) strain rate

      dtemp11 = p5*(divu + tension)
      dtemp12 = shear*p5
      dtemp22 = p5*(divu - tension)

      ! here Tan2alpha = 2 dtemp12 / (dtemp11 - dtemp22)

      Qd11Qd11 = c1
      Qd12Qd12 = puny
      Qd11Qd12 = puny

      if ((ABS( dtemp11 - dtemp22) > puny).or.(ABS(dtemp12) > puny)) then
         Angle_denom_alpha = sqrt( ( dtemp11 - dtemp22 )* &
                ( dtemp11 - dtemp22 ) + c4*dtemp12*dtemp12)

         Qd11Qd11 = p5 + ( dtemp11 - dtemp22 )*p5/Angle_denom_alpha   !Cos^2
         Qd12Qd12 = p5 - ( dtemp11 - dtemp22 )*p5/Angle_denom_alpha   !Sin^2
         Qd11Qd12 = dtemp12/Angle_denom_alpha                         !CosSin
      endif

      dtemp1 = Qd11Qd11*dtemp11 + c2*Qd11Qd12*dtemp12 + Qd12Qd12*dtemp22
      dtemp2 = Qd12Qd12*dtemp11 - c2*Qd11Qd12*dtemp12 + Qd11Qd11*dtemp22

      ! In cos and sin values
      x = c0

      if ((ABS(dtemp1) > puny).or.(ABS(dtemp2) > puny)) then
!         invleng = c1/sqrt(dtemp1*dtemp1 + dtemp2*dtemp2) ! not sure if this is neccessary
!         dtemp1 = dtemp1*invleng
!         dtemp2 = dtemp2*invleng
         if (dtemp1 == c0) then
            x = pih
         else
            x = atan2(dtemp2,dtemp1)
         endif
      endif

      !echmod to ensure the angle lies between pi/4 and 9 pi/4
      if (x < piq) x = x + pi2
      !echmod require 0 <= x < (nx_yield-1)*dx = 2 pi
!      x = mod(x+pi2, pi2)
      ! y: angle between major principal axis of strain rate and structure tensor
      ! y = gamma - alpha
      ! Expressesed componently with
      ! Tany = (Singamma*Cosgamma - Sinalpha*Cosgamma)/(Cos^2gamma - Sin^alpha)

      Tany_1 = Q11Q12 - Qd11Qd12
      Tany_2 = Q11Q11 - Qd12Qd12

      y = c0

      if ((ABS(Tany_1) > puny).or.(ABS(Tany_2) > puny)) then
!         invleng = c1/sqrt(Tany_1*Tany_1 + Tany_2*Tany_2) ! not sure if this is neccessary
!         Tany_1 = Tany_1*invleng
!         Tany_2 = Tany_2*invleng
         if (Tany_2 == c0) then
            y = pih
         else
            y = atan2(Tany_1,Tany_2)
         endif
      endif

      ! to make sure y is between 0 and pi

      if (y > pi) y = y - pi
      if (y < 0)  y = y + pi

      if (interpolate_stress_rdg) then

         ! Interpolated lookup

         ! if (x>=9*pi/4) x=9*pi/4-puny; end
         ! if (y>=pi/2)   y=pi/2-puny; end
         ! if (atempprime>=1.0), atempprime=1.0-puny; end

         ! % need 8 coords and 8 weights
         ! % range in kx

         kx  = int((x-piq-pi)*invdx) + 1
         kxw = c1 - ((x-piq-pi)*invdx - (kx-1))

         ky  = int(y*invdy) + 1
         kyw = c1 - (y*invdy - (ky-1))

         ka  = int((atempprime-p5)*invda) + 1
         kaw = c1 - ((atempprime-p5)*invda - (ka-1))

         ! % Determine sigma_r(A1,Zeta,y) and sigma_s (see Section A1)

         stemp11r =  kxw* kyw      * kaw      * s11r(kx  ,ky  ,ka  ) &
             + (c1-kxw) * kyw      * kaw      * s11r(kx+1,ky  ,ka  ) &
             + kxw      * (c1-kyw) * kaw      * s11r(kx  ,ky+1,ka  ) &
             + kxw      * kyw      * (c1-kaw) * s11r(kx  ,ky  ,ka+1) &
             + (c1-kxw) * (c1-kyw) * kaw      * s11r(kx+1,ky+1,ka  ) &
             + (c1-kxw) * kyw      * (c1-kaw) * s11r(kx+1,ky  ,ka+1) &
             + kxw      * (c1-kyw) * (c1-kaw) * s11r(kx  ,ky+1,ka+1) &
             + (c1-kxw) * (c1-kyw) * (c1-kaw) * s11r(kx+1,ky+1,ka+1)

         stemp12r =  kxw* kyw      * kaw      * s12r(kx  ,ky  ,ka  ) &
             + (c1-kxw) * kyw      * kaw      * s12r(kx+1,ky  ,ka  ) &
             + kxw      * (c1-kyw) * kaw      * s12r(kx  ,ky+1,ka  ) &
             + kxw      * kyw      * (c1-kaw) * s12r(kx  ,ky  ,ka+1) &
             + (c1-kxw) * (c1-kyw) * kaw      * s12r(kx+1,ky+1,ka  ) &
             + (c1-kxw) * kyw      * (c1-kaw) * s12r(kx+1,ky  ,ka+1) &
             + kxw      * (c1-kyw) * (c1-kaw) * s12r(kx  ,ky+1,ka+1) &
             + (c1-kxw) * (c1-kyw) * (c1-kaw) * s12r(kx+1,ky+1,ka+1)

         stemp22r = kxw * kyw      * kaw      * s22r(kx  ,ky  ,ka  ) &
             + (c1-kxw) * kyw      * kaw      * s22r(kx+1,ky  ,ka  ) &
             + kxw      * (c1-kyw) * kaw      * s22r(kx  ,ky+1,ka  ) &
             + kxw      * kyw      * (c1-kaw) * s22r(kx  ,ky  ,ka+1) &
             + (c1-kxw) * (c1-kyw) * kaw      * s22r(kx+1,ky+1,ka  ) &
             + (c1-kxw) * kyw      * (c1-kaw) * s22r(kx+1,ky  ,ka+1) &
             + kxw      * (c1-kyw) * (c1-kaw) * s22r(kx  ,ky+1,ka+1) &
             + (c1-kxw) * (c1-kyw) * (c1-kaw) * s22r(kx+1,ky+1,ka+1)

         stemp11s =  kxw* kyw      * kaw      * s11s(kx  ,ky  ,ka  ) &
             + (c1-kxw) * kyw      * kaw      * s11s(kx+1,ky  ,ka  ) &
             + kxw      * (c1-kyw) * kaw      * s11s(kx  ,ky+1,ka  ) &
             + kxw      * kyw      * (c1-kaw) * s11s(kx  ,ky  ,ka+1) &
             + (c1-kxw) * (c1-kyw) * kaw      * s11s(kx+1,ky+1,ka  ) &
             + (c1-kxw) * kyw      * (c1-kaw) * s11s(kx+1,ky  ,ka+1) &
             + kxw      * (c1-kyw) * (c1-kaw) * s11s(kx  ,ky+1,ka+1) &
             + (c1-kxw) * (c1-kyw) * (c1-kaw) * s11s(kx+1,ky+1,ka+1)

         stemp12s =  kxw* kyw      * kaw      * s12s(kx  ,ky  ,ka  ) &
             + (c1-kxw) * kyw      * kaw      * s12s(kx+1,ky  ,ka  ) &
             + kxw      * (c1-kyw) * kaw      * s12s(kx  ,ky+1,ka  ) &
             + kxw      * kyw      * (c1-kaw) * s12s(kx  ,ky  ,ka+1) &
             + (c1-kxw) * (c1-kyw) * kaw      * s12s(kx+1,ky+1,ka  ) &
             + (c1-kxw) * kyw      * (c1-kaw) * s12s(kx+1,ky  ,ka+1) &
             + kxw      * (c1-kyw) * (c1-kaw) * s12s(kx  ,ky+1,ka+1) &
             + (c1-kxw) * (c1-kyw) * (c1-kaw) * s12s(kx+1,ky+1,ka+1)

         stemp22s =  kxw* kyw      * kaw      * s22s(kx  ,ky  ,ka  ) &
             + (c1-kxw) * kyw      * kaw      * s22s(kx+1,ky  ,ka  ) &
             + kxw      * (c1-kyw) * kaw      * s22s(kx  ,ky+1,ka  ) &
             + kxw      * kyw      * (c1-kaw) * s22s(kx  ,ky  ,ka+1) &
             + (c1-kxw) * (c1-kyw) * kaw      * s22s(kx+1,ky+1,ka  ) &
             + (c1-kxw) * kyw      * (c1-kaw) * s22s(kx+1,ky  ,ka+1) &
             + kxw      * (c1-kyw) * (c1-kaw) * s22s(kx  ,ky+1,ka+1) &
             + (c1-kxw) * (c1-kyw) * (c1-kaw) * s22s(kx+1,ky+1,ka+1)

      else

         kx = int((x-piq-pi)*invdx) + 1
         ky = int(y*invdy) + 1
         ka = int((atempprime-p5)*invda) + 1

         ! Determine sigma_r(A1,Zeta,y) and sigma_s (see Section A1)

         stemp11r = s11r(kx,ky,ka)
         stemp12r = s12r(kx,ky,ka)
         stemp22r = s22r(kx,ky,ka)

         stemp11s = s11s(kx,ky,ka)
         stemp12s = s12s(kx,ky,ka)
         stemp22s = s22s(kx,ky,ka)

      endif

      ! Calculate mean ice stress over a collection of floes (Equation 3)

      stressp  = strength*(stemp11r + kfriction*stemp11s &
                         + stemp22r + kfriction*stemp22s) * invsin
      stress12 = strength*(stemp12r + kfriction*stemp12s) * invsin
      stressm  = strength*(stemp11r + kfriction*stemp11s &
                         - stemp22r - kfriction*stemp22s) * invsin

      ! Back - rotation of the stress from principal axes into general coordinates

      ! Update stress

      sig11 = p5*(stressp + stressm)
      sig12 = stress12
      sig22 = p5*(stressp - stressm)

      sgprm11 = Q11Q11*sig11 + Q12Q12*sig22 -        c2*Q11Q12 *sig12
      sgprm12 = Q11Q12*sig11 - Q11Q12*sig22 + (Q11Q11 - Q12Q12)*sig12
      sgprm22 = Q12Q12*sig11 + Q11Q11*sig22 +        c2*Q11Q12 *sig12

      stressp  = sgprm11 + sgprm22
      stress12 = sgprm12
      stressm  = sgprm11 - sgprm22

      ! Compute ridging and sliding functions in general coordinates (Equation 11)

      if (ksub == ndte) then
         rotstemp11r = Q11Q11*stemp11r - c2*Q11Q12* stemp12r &
                     + Q12Q12*stemp22r
         rotstemp12r = Q11Q11*stemp12r +    Q11Q12*(stemp11r-stemp22r) &
                     - Q12Q12*stemp12r
         rotstemp22r = Q12Q12*stemp11r + c2*Q11Q12* stemp12r &
                     + Q11Q11*stemp22r

         rotstemp11s = Q11Q11*stemp11s - c2*Q11Q12* stemp12s &
                     + Q12Q12*stemp22s
         rotstemp12s = Q11Q11*stemp12s +    Q11Q12*(stemp11s-stemp22s) &
                     - Q12Q12*stemp12s
         rotstemp22s = Q12Q12*stemp11s + c2*Q11Q12* stemp12s &
                     + Q11Q11*stemp22s

         alphar =  rotstemp11r*dtemp11 + c2*rotstemp12r*dtemp12 &
                 + rotstemp22r*dtemp22
         alphas =  rotstemp11s*dtemp11 + c2*rotstemp12s*dtemp12 &
                 + rotstemp22s*dtemp22
      endif

      end subroutine update_stress_rdg

!=======================================================================
! Solves evolution equation for structure tensor (A19, A20)

      subroutine stepa  (nx_block,   ny_block,       &
                         dtei,       icellT,         &
                         indxTi,     indxTj,         &
                         Xmap,                       &
                         a11, a12,                   &
                         a11_1, a11_2, a11_3, a11_4, &
                         a12_1, a12_2, a12_3, a12_4, &
                         stressp_1,  stressp_2,      &
                         stressp_3,  stressp_4,      &
                         stressm_1,  stressm_2,      &
                         stressm_3,  stressm_4,      &
                         stress12_1, stress12_2,     &
                         stress12_3, stress12_4)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellT                ! no. of cells where iceTmask = .true.

      real (kind=dbl_kind), intent(in) :: &
         dtei        ! 1/dte, where dte is subcycling timestep (1/s)

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxTi  , & ! compressed index in i-direction
         indxTj      ! compressed index in j-direction

      ! Xmap is 1D over the block in Fortran storage order: element (i,j) at
      ! (j-1)*nx_block+i, holding that cell's compressed slot
      integer (kind=int_kind), dimension (:), intent(in) :: &
         Xmap        ! compressed slot of each cell

      ! COMPRESSED: this block's slice, indexed with is.  These are the state
      ! stress_eap and stepa pass between them every subcycle, so both read and
      ! write the same slots and the block arrays are converted only either
      ! side of the loop, where restart and history want them.
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         ! ice stress tensor (kg/s^2) in each corner of T cell
          stressp_1,  stressp_2,  stressp_3,  stressp_4, & ! sigma11+sigma22
          stressm_1,  stressm_2,  stressm_3,  stressm_4, & ! sigma11-sigma22
         stress12_1, stress12_2, stress12_3, stress12_4    ! sigma12

      real (kind=dbl_kind), dimension (:), intent(inout), contiguous :: &
         ! structure tensor () in each corner of T cell, compressed
         a11_1, a11_2, a11_3, a11_4, & ! components of
         a12_1, a12_2, a12_3, a12_4    ! structure tensor ()

      ! a11 and a12 are the averages this routine produces.  stress_eap reads
      ! a11_1..4 and a12_1..4, never these, so they are outputs rather than
      ! subcycle state and stay in the block layout.
      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         a11, a12

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij, &
         ijb      , & ! 1D index of (i,j) in the block
         is           ! compressed slot of that cell

      real (kind=dbl_kind) :: &
         mresult11, mresult12, &
         dteikth, p5kth

      real (kind=dbl_kind), parameter :: &
         kth  = p2*p001

      character(len=*), parameter :: subname = '(stepa)'

      dteikth = c1 / (dtei + kth)
      p5kth = p5 * kth

      do ij = 1, icellT
         i = indxTi(ij)
         j = indxTj(ij)
         ijb = (j-1)*nx_block + i
         is  = Xmap(ijb)

         ! ne
         call calc_ffrac(stressp_1(is), stressm_1(is),  &
                         stress12_1(is),                 &
                         a11_1(is), a12_1(is),          &
                         mresult11, mresult12)

         a11_1(is) = (a11_1(is)*dtei + p5kth - mresult11) * dteikth ! implicit
         a12_1(is) = (a12_1(is)*dtei - mresult12) * dteikth ! implicit


         ! nw
         call calc_ffrac(stressp_2(is), stressm_2(is),  &
                         stress12_2(is),                 &
                         a11_2(is), a12_2(is),          &
                         mresult11, mresult12)

         a11_2(is) = (a11_2(is)*dtei + p5kth - mresult11) * dteikth ! implicit
         a12_2(is) = (a12_2(is)*dtei - mresult12) * dteikth ! implicit

         ! sw
         call calc_ffrac(stressp_3(is), stressm_3(is),  &
                         stress12_3(is),                 &
                         a11_3(is), a12_3(is),          &
                         mresult11, mresult12)

         a11_3(is) = (a11_3(is)*dtei + p5kth - mresult11) * dteikth ! implicit
         a12_3(is) = (a12_3(is)*dtei - mresult12) * dteikth ! implicit

         ! se
         call calc_ffrac(stressp_4(is), stressm_4(is),  &
                         stress12_4(is),                 &
                         a11_4(is), a12_4(is),          &
                         mresult11, mresult12)

         a11_4(is) = (a11_4(is)*dtei + p5kth - mresult11) * dteikth ! implicit
         a12_4(is) = (a12_4(is)*dtei - mresult12) * dteikth ! implicit

         ! average structure tensor

         a11(i,j) = p25*(a11_1(is) + a11_2(is) + a11_3(is) + a11_4(is))
         a12(i,j) = p25*(a12_1(is) + a12_2(is) + a12_3(is) + a12_4(is))

      enddo                     ! ij

      end subroutine stepa

!=======================================================================
! computes term in evolution equation for structure tensor which determines
! the ice floe re-orientation due to fracture
! Eq. 7: Ffrac = -kf(A-S) or = 0 depending on sigma_1 and sigma_2


      subroutine calc_ffrac (stressp, stressm, &
                             stress12,         &
                             a1x,  a2x,        &
                             mresult1, mresult2)

      real (kind=dbl_kind), intent(in) :: &
         stressp, stressm, stress12, a1x, a2x

      real (kind=dbl_kind), intent(out) :: &
         mresult1, mresult2

      ! local variables

      real (kind=dbl_kind) :: &
         sigma11, sigma12, sigma22, &
         gamma, sigma_1, sigma_2, &
         Q11, Q12, Q11Q11, Q11Q12, Q12Q12

      real (kind=dbl_kind), parameter :: &
         kfrac = p001, &
         threshold = c3*p1

      character(len=*), parameter :: subname = '(calc_ffrac)'

      sigma11 = p5*(stressp+stressm)
      sigma12 = stress12
      sigma22 = p5*(stressp-stressm)

!      if ((sigma11-sigma22) == c0) then sigma11-sigma22 == 0 => stressn ==0
      if (stressm == c0) then
         gamma = p5*(pih)
      else
         gamma = p5*atan2((c2*sigma12),(sigma11-sigma22))
      endif

      ! rotate tensor to get into sigma principal axis

      Q11 = cos(gamma)
      Q12 = sin(gamma)

      Q11Q11 = Q11*Q11
      Q11Q12 = Q11*Q12
      Q12Q12 = Q12*Q12

      sigma_1 = Q11Q11*sigma11 + c2*Q11Q12*sigma12  &
              + Q12Q12*sigma22 ! S(1,1)
      sigma_2 = Q12Q12*sigma11 - c2*Q11Q12*sigma12  &
              + Q11Q11*sigma22 ! S(2,2)

      ! Pure divergence
      if ((sigma_1 >= c0).and.(sigma_2 >= c0))  then
         mresult1 = c0
         mresult2 = c0

      ! Unconfined compression: cracking of blocks not along the axial splitting direction
      ! which leads to the loss of their shape, so we again model it through diffusion
      elseif ((sigma_1 >= c0).and.(sigma_2 < c0))  then
         mresult1 = kfrac * (a1x - Q12Q12)
         mresult2 = kfrac * (a2x + Q11Q12)

      ! Shear faulting
      elseif (sigma_2 == c0) then
         mresult1 = c0
         mresult2 = c0
      elseif ((sigma_1 <= c0).and.(sigma_1/sigma_2 <= threshold)) then
         mresult1 = kfrac * (a1x - Q12Q12)
         mresult2 = kfrac * (a2x + Q11Q12)

      ! Horizontal spalling
      else
         mresult1 = c0
         mresult2 = c0
      endif

      end subroutine calc_ffrac

!=======================================================================
!---! these subroutines write/read Fortran unformatted data files ..
!=======================================================================
! Dumps all values needed for a restart

      subroutine write_restart_eap ()

      use ice_restart, only: write_restart_field

      ! local variables

      logical (kind=log_kind) :: diag

      character(len=*), parameter :: subname = '(write_restart_eap)'

      diag = .true.

      !-----------------------------------------------------------------
      ! structure tensor
      !-----------------------------------------------------------------

      call write_restart_field(nu_dump_eap,0,a11_1,'ruf8','a11_1',1,diag)
      call write_restart_field(nu_dump_eap,0,a11_3,'ruf8','a11_3',1,diag)
      call write_restart_field(nu_dump_eap,0,a11_2,'ruf8','a11_2',1,diag)
      call write_restart_field(nu_dump_eap,0,a11_4,'ruf8','a11_4',1,diag)

      call write_restart_field(nu_dump_eap,0,a12_1,'ruf8','a12_1',1,diag)
      call write_restart_field(nu_dump_eap,0,a12_3,'ruf8','a12_3',1,diag)
      call write_restart_field(nu_dump_eap,0,a12_2,'ruf8','a12_2',1,diag)
      call write_restart_field(nu_dump_eap,0,a12_4,'ruf8','a12_4',1,diag)

      end subroutine write_restart_eap

!=======================================================================
! Reads all values needed for elastic anisotropic plastic dynamics restart

      subroutine read_restart_eap()

      use ice_blocks, only: nghost, ns_boundary_type
      use ice_boundary, only: ice_HaloUpdate_stress
      use ice_constants, only:  &
          field_loc_center, field_type_scalar
      use ice_domain, only: nblocks, halo_info
      use ice_restart, only: read_restart_field

      ! local variables

      integer (kind=int_kind) :: &
         i, j, iblk

      logical (kind=log_kind) :: &
         diag

      character(len=*), parameter :: subname = '(read_restart_eap)'

      diag = .true.

      !-----------------------------------------------------------------
      ! Structure tensor must be read and scattered in pairs in order
      ! to properly match corner values across a tripole grid cut.
      !-----------------------------------------------------------------
      if (my_task == master_task) write(nu_diag,*) &
         'structure tensor restart data'

      call read_restart_field(nu_restart_eap,0,a11_1,'ruf8', &
           'a11_1',1,diag,field_loc_center,field_type_scalar) ! a11_1
      call read_restart_field(nu_restart_eap,0,a11_3,'ruf8', &
           'a11_3',1,diag,field_loc_center,field_type_scalar) ! a11_3
      call read_restart_field(nu_restart_eap,0,a11_2,'ruf8', &
           'a11_2',1,diag,field_loc_center,field_type_scalar) ! a11_2
      call read_restart_field(nu_restart_eap,0,a11_4,'ruf8', &
           'a11_4',1,diag,field_loc_center,field_type_scalar) ! a11_4

      call read_restart_field(nu_restart_eap,0,a12_1,'ruf8', &
           'a12_1',1,diag,field_loc_center,field_type_scalar) ! a12_1
      call read_restart_field(nu_restart_eap,0,a12_3,'ruf8', &
           'a12_3',1,diag,field_loc_center,field_type_scalar) ! a12_3
      call read_restart_field(nu_restart_eap,0,a12_2,'ruf8', &
           'a12_2',1,diag,field_loc_center,field_type_scalar) ! a12_2
      call read_restart_field(nu_restart_eap,0,a12_4,'ruf8', &
           'a12_4',1,diag,field_loc_center,field_type_scalar) ! a12_4

      if (trim(ns_boundary_type) == 'tripole' .or. &
          trim(ns_boundary_type) == 'tripoleT') then

         call ice_HaloUpdate_stress(a11_1, a11_3,      halo_info, &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(a11_3, a11_1,      halo_info, &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(a11_2, a11_4,      halo_info, &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(a11_4, a11_2,      halo_info, &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(a12_1, a12_3,      halo_info, &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(a12_3, a12_1,      halo_info, &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(a12_2, a12_4,      halo_info, &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(a12_4, a12_2,      halo_info, &
                                    field_loc_center,  field_type_scalar)

      endif

      !-----------------------------------------------------------------
      ! Ensure unused values in west and south ghost cells are 0
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk,i,j) SCHEDULE(runtime)
      do iblk = 1, nblocks
         do j = 1, nghost
         do i = 1, nx_block
            a11_1 (i,j,iblk) = c0
            a11_2 (i,j,iblk) = c0
            a11_3 (i,j,iblk) = c0
            a11_4 (i,j,iblk) = c0
            a12_1 (i,j,iblk) = c0
            a12_2 (i,j,iblk) = c0
            a12_3 (i,j,iblk) = c0
            a12_4 (i,j,iblk) = c0
         enddo
         enddo
         do j = 1, ny_block
         do i = 1, nghost
            a11_1 (i,j,iblk) = c0
            a11_2 (i,j,iblk) = c0
            a11_3 (i,j,iblk) = c0
            a11_4 (i,j,iblk) = c0
            a12_1 (i,j,iblk) = c0
            a12_2 (i,j,iblk) = c0
            a12_3 (i,j,iblk) = c0
            a12_4 (i,j,iblk) = c0
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine read_restart_eap

!=======================================================================

      end module ice_dyn_eap

!=======================================================================
