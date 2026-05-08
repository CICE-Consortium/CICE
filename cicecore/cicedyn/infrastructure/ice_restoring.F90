!=======================================================================
!
! Reads and interpolates forcing data for atmosphere and ocean quantities.
!
! authors: Elizabeth C. Hunke, LANL

   module ice_restoring

#if (1 == 1)
!      use ice_arrays_column, only:ffracn,dhsn,oceanmixed_ice
      use ice_kinds_mod
      use ice_blocks, only: block, get_block, nx_block, ny_block, nghost, &
          nblocks_x, nblocks_y
      use ice_communicate, only: my_task, master_task
      use ice_constants, only: c0, c1, c2, p2, p5, c4
      use ice_domain_size, only: ncat, max_blocks, nilyr, nslyr
      use ice_domain, only: nblocks, blocks_ice, bdy_origin, &
          ew_boundary_type, ns_boundary_type, &
          max_set_boundary_flds, num_set_boundary_flds, set_boundary_flds
!! tcraig, iceUmask is a logical, can't restore it
!!      use ice_dyn_shared, only: iceUmask
      use ice_exit, only: abort_ice
      use ice_fileunits, only: nu_diag
!      use ice_flux, only: stressp_1, stressp_2, stressp_3, &
!          stressp_4, stressm_1, stressm_2, stressm_3, &
!          stressm_4, stress12_1, stress12_2, stress12_3, stress12_4, &
!          swvdr, swvdf, swidr, swidf, strocnyT_iavg, strocnxT_iavg, &
!          scale_factor, frz_onset, fsnow, frzmlt,  sst
      use ice_forcing, only: trestore, get_forcing_bry, &
          aicen_bry, vicen_bry, vsnon_bry, alvl_bry, qsno_bry, ffrac_bry, &
          Tsfc_bry, Tinz_bry, Sinz_bry, vlvl_bry, FY_bry, dhs_bry, &
          apnd_bry, hpnd_bry, ipnd_bry, iage_bry, uvel_bry, &
          vvel_bry, scale_factor_bry, swvdr_bry, swvdf_bry, swidr_bry, swidf_bry, &
          strocnxT_bry, strocnyT_bry, stressp_1_bry, stressp_2_bry, &
          stressp_3_bry, stressp_4_bry, stressm_1_bry, &
          stressm_2_bry, stressm_3_bry, stressm_4_bry, &
          stress12_1_bry, stress12_2_bry, stress12_3_bry, &
          stress12_4_bry, iceumask_bry, frz_onset_bry, fsnow_bry,  &
          frzmelt_bry, sst_bry
!      use ice_grid, only: tmask, hm
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
!      use icepack_intfc, only: icepack_init_trcr
      use icepack_intfc, only: icepack_query_parameters, &
          icepack_query_tracer_sizes, icepack_query_tracer_flags, &
          icepack_query_tracer_indices
      use ice_state, only: aicen, vicen, vsnon, trcrn, uvel, vvel
!          aice_init, aice0, aice, vice, vsno, trcr, &
!          trcr_depend, &
!          divu, shear, strength
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bound
#endif

      implicit none
      private
      public :: ice_HaloRestore_init, ice_HaloRestore, ice_HaloRestore_getbdy

      logical (kind=log_kind), public :: &
         restore_ice                 ! restore ice state if true

      character (len=7), parameter :: &
!         restore_ic = 'defined' ! restore to internally defined ice state
         restore_ic = 'initial' ! restore to initial ice state

      !-----------------------------------------------------------------
      ! state of the ice for each category
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable :: &
         aicen_rest , & ! concentration of ice
         vicen_rest , & ! volume per unit area of ice          (m)
         vsnon_rest ,&  ! volume per unit area of snow         (m)
         dhs_rest ,&
         ffrac_rest

      real (kind=dbl_kind), dimension (:,:,:), allocatable :: &
         uvel_rest, &
         vvel_rest, &
         scale_factor_rest, &
         swvdr_rest, &
         swvdf_rest, &
         swidr_rest, &
         swidf_rest, &
         strocnxT_rest, &
         strocnyT_rest, &
         stressp_1_rest, &
         stressp_2_rest, &
         stressp_3_rest, &
         stressp_4_rest, &
         stressm_1_rest, &
         stressm_2_rest, &
         stressm_3_rest, &
         stressm_4_rest, &
         stress12_1_rest, &
         stress12_2_rest, &
         stress12_3_rest, &
         stress12_4_rest, &
!         iceumask_rest, &
         frz_onset_rest, &
         fsnow_rest, &
         sst_rest, &
         frzmelt_rest

      real (kind=dbl_kind), dimension (:,:,:,:,:), allocatable :: &
         trcrn_rest     ! tracers

      integer(kind=int_kind), parameter :: &
         nfact=0        ! how far to restore into grid, 0=just halo

!=======================================================================

      contains

!=======================================================================

!  Allocates and initializes arrays needed for restoring the ice state
!  in cells surrounding the grid.


   subroutine ice_HaloRestore_init

#if (1 == 0)
!      use ice_flux, only: Tf, Tair, salinz, Tmltz
#endif

      integer (int_kind) :: &
         i,j,iblk,nt,n,k,    &! dummy loop indices
         ilo,ihi,jlo,jhi,    &! beginning and end of physical domain
         iglob(nx_block),    &! global indices
         jglob(ny_block),    &! global indices
         iblock, jblock,     &! block indices
         ntrcr                !

      type (block) :: &
         this_block  ! block info for current block

      character(len=*), parameter :: subname = '(ice_HaloRestore_init)'

      if (.not. restore_ice) return

      call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      allocate (aicen_rest(nx_block,ny_block,ncat,max_blocks), &
                vicen_rest(nx_block,ny_block,ncat,max_blocks), &
                vsnon_rest(nx_block,ny_block,ncat,max_blocks), &
                trcrn_rest(nx_block,ny_block,ntrcr,ncat,max_blocks))

      aicen_rest(:,:,:,:) = c0
      vicen_rest(:,:,:,:) = c0
      vsnon_rest(:,:,:,:) = c0
      trcrn_rest(:,:,:,:,:) = c0

      !-----------------------------------------------------------------------
      ! initialize
      ! halo cells have to be filled manually at this stage
      !-----------------------------------------------------------------------

#if (1 == 0)
      if (trim(restore_ic) == 'defined') then

         ! restore to defined ice state
         !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block, &
         !$OMP                     iglob,jglob,iblock,jblock)
         do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            iglob = this_block%i_glob
            jglob = this_block%j_glob
            iblock = this_block%iblock
            jblock = this_block%jblock

            call set_restore_var (nx_block,            ny_block,            &
                                  ilo, ihi,            jlo, jhi,            &
                                  iglob,               jglob,               &
                                  iblock,              jblock,              &
                                  Tair (:,:,    iblk), &
                                  Tf   (:,:,    iblk),                      &
                                  salinz(:,:,:, iblk), Tmltz(:,:,:,  iblk), &
                                  tmask(:,:,    iblk),                      &
                                  aicen_rest(:,:,  :,iblk), &
                                  trcrn_rest(:,:,:,:,iblk), ntrcr,         &
                                  vicen_rest(:,:,  :,iblk), &
                                  vsnon_rest(:,:,  :,iblk))
         enddo ! iblk
         !$OMP END PARALLEL DO

      else  ! restore_ic

         ! restore to initial ice state

         aicen_rest(:,:,:,:) = aicen(:,:,:,:)
         vicen_rest(:,:,:,:) = vicen(:,:,:,:)
         vsnon_rest(:,:,:,:) = vsnon(:,:,:,:)
         trcrn_rest(:,:,:,:,:) = trcrn(:,:,:,:,:)

      endif ! restore_ic

      !-----------------------------------------------------------------
      ! Impose land mask
      !-----------------------------------------------------------------

      do iblk = 1, nblocks
         do n = 1, ncat
            do j = 1, ny_block
            do i = 1, nx_block
               aicen_rest(i,j,n,iblk) = aicen_rest(i,j,n,iblk) * hm(i,j,iblk)
               vicen_rest(i,j,n,iblk) = vicen_rest(i,j,n,iblk) * hm(i,j,iblk)
               vsnon_rest(i,j,n,iblk) = vsnon_rest(i,j,n,iblk) * hm(i,j,iblk)
               do nt = 1, ntrcr
                  trcrn_rest(i,j,nt,n,iblk) = trcrn_rest(i,j,nt,n,iblk) &
                                                               * hm(i,j,iblk)
               enddo
            enddo
            enddo
         enddo
      enddo

      if (bdy_origin == 'restart_f') then
#endif
         allocate (uvel_rest        (nx_block,ny_block,max_blocks), &
                   vvel_rest        (nx_block,ny_block,max_blocks), &
                   scale_factor_rest(nx_block,ny_block,max_blocks), &
                   swvdr_rest       (nx_block,ny_block,max_blocks), &
                   swvdf_rest       (nx_block,ny_block,max_blocks), &
                   swidr_rest       (nx_block,ny_block,max_blocks), &
                   swidf_rest       (nx_block,ny_block,max_blocks), &
                   strocnxT_rest    (nx_block,ny_block,max_blocks), &
                   strocnyT_rest    (nx_block,ny_block,max_blocks), &
                   stressp_1_rest   (nx_block,ny_block,max_blocks), &
                   stressp_2_rest   (nx_block,ny_block,max_blocks), &
                   stressp_3_rest   (nx_block,ny_block,max_blocks), &
                   stressp_4_rest   (nx_block,ny_block,max_blocks), &
                   stressm_1_rest   (nx_block,ny_block,max_blocks), &
                   stressm_2_rest   (nx_block,ny_block,max_blocks), &
                   stressm_3_rest   (nx_block,ny_block,max_blocks), &
                   stressm_4_rest   (nx_block,ny_block,max_blocks), &
                   stress12_1_rest  (nx_block,ny_block,max_blocks), &
                   stress12_2_rest  (nx_block,ny_block,max_blocks), &
                   stress12_3_rest  (nx_block,ny_block,max_blocks), &
                   stress12_4_rest  (nx_block,ny_block,max_blocks), &
!                   iceumask_rest    (nx_block,ny_block,max_blocks), &
                   frz_onset_rest   (nx_block,ny_block,max_blocks), &
                   fsnow_rest       (nx_block,ny_block,max_blocks),  &
                   ffrac_rest       (nx_block,ny_block,ncat,max_blocks),  &
                   dhs_rest         (nx_block,ny_block,ncat,max_blocks),  &
                   sst_rest         (nx_block,ny_block,max_blocks),  &
                   frzmelt_rest     (nx_block,ny_block,max_blocks))

         uvel_rest(:,:,:) = c0
         vvel_rest(:,:,:) = c0
         scale_factor_rest(:,:,:) = c0
         swvdr_rest(:,:,:) = c0
         swvdf_rest(:,:,:) = c0
         swidr_rest(:,:,:) = c0
         swidf_rest(:,:,:) = c0
         strocnxT_rest(:,:,:) = c0
         strocnyT_rest(:,:,:) = c0
         stressp_1_rest(:,:,:) = c0
         stressp_2_rest(:,:,:) = c0
         stressp_3_rest(:,:,:) = c0
         stressp_4_rest(:,:,:) = c0
         stressm_1_rest(:,:,:) = c0
         stressm_2_rest(:,:,:) = c0
         stressm_3_rest(:,:,:) = c0
         stressm_4_rest(:,:,:) = c0
         stress12_1_rest(:,:,:) = c0
         stress12_2_rest(:,:,:) = c0
         stress12_3_rest(:,:,:) = c0
         stress12_4_rest(:,:,:) = c0
!         iceumask_rest(:,:,:) = c0
         frz_onset_rest(:,:,:) = c0
         fsnow_rest(:,:,:) = c0
         sst_rest(:,:,:) = c0
         frzmelt_rest(:,:,:) = c0
         ffrac_rest(:,:,:,:) = c0
         dhs_rest(:,:,:,:) = c0

!1 == 0      endif

   end subroutine ice_HaloRestore_init

!=======================================================================
!
   subroutine ice_HaloRestore_getbdy

      integer(kind=int_kind) :: &
         k            ! dummy arguments

      integer (kind=int_kind) :: ntrcr
      logical (kind=log_kind) :: tr_iage, tr_FY, tr_lvl, tr_pond, tr_aero, tr_fsd
      logical (kind=log_kind) :: tr_snow, tr_brine
      logical (kind=log_kind) :: tr_iso, tr_pond_lvl, tr_pond_topo, tr_pond_sealvl
      integer (kind=int_kind) :: nt_Tsfc, nt_sice, nt_qice, nt_qsno, nt_iage, nt_FY
      integer (kind=int_kind) :: nt_alvl, nt_vlvl, nt_apnd, nt_hpnd, nt_ipnd, nt_aero
      integer (kind=int_kind) :: nt_fsd, nt_isosno, nt_isoice, nt_fbri
      integer (kind=int_kind) :: nt_smice, nt_smliq, nt_rhos, nt_rsnw

      character(len=*), parameter :: subname = '(ice_HaloRestore_getbdy)'

      call icepack_query_tracer_flags(tr_iage_out=tr_iage, tr_FY_out=tr_FY, &
         tr_lvl_out=tr_lvl, tr_aero_out=tr_aero, tr_pond_out=tr_pond, &
         tr_pond_lvl_out=tr_pond_lvl, tr_pond_sealvl_out=tr_pond_sealvl, &
         tr_pond_topo_out=tr_pond_topo, tr_brine_out=tr_brine, tr_fsd_out=tr_fsd, &
         tr_snow_out=tr_snow, tr_iso_out=tr_iso)
      call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_sice_out=nt_sice, &
         nt_qice_out=nt_qice, nt_qsno_out=nt_qsno, nt_iage_out=nt_iage, nt_fy_out=nt_fy, &
         nt_alvl_out=nt_alvl, nt_vlvl_out=nt_vlvl, nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, &
         nt_ipnd_out=nt_ipnd, nt_fsd_out=nt_fsd, nt_aero_out=nt_aero, &
         nt_smice_out=nt_smice, nt_smliq_out=nt_smliq, nt_rhos_out=nt_rhos, nt_rsnw_out=nt_rsnw, &
         nt_isosno_out=nt_isosno,     nt_isoice_out=nt_isoice,       nt_fbri_out=nt_fbri)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__, line=__LINE__)

      call get_forcing_bry()

      trcrn_rest(:,:,nt_Tsfc,:,:) = Tsfc_bry(:,:,:,:)
      aicen_rest(:,:,:,:) = aicen_bry(:,:,:,:)
      vicen_rest(:,:,:,:) = vicen_bry(:,:,:,:)
      vsnon_rest(:,:,:,:) = vsnon_bry(:,:,:,:)

      do k = 1,nilyr
         trcrn_rest(:,:,nt_sice+k-1,:,:) = Sinz_bry(:,:,k,:,:)
         trcrn_rest(:,:,nt_qice+k-1,:,:) = Tinz_bry(:,:,k,:,:)
      enddo
      do k = 1,nslyr
         trcrn_rest(:,:,nt_qsno+k-1,:,:) = qsno_bry(:,:,k,:,:)
      enddo

      if (tr_pond_lvl) then ! bdy from file
         fsnow_rest(:,:,:) = fsnow_bry(:,:,:)
         trcrn_rest(:,:,nt_apnd,:,:) = apnd_bry(:,:,:,:)
         trcrn_rest(:,:,nt_ipnd,:,:) = ipnd_bry(:,:,:,:)
         trcrn_rest(:,:,nt_hpnd,:,:) = hpnd_bry(:,:,:,:)
!         dhs_rest(:,:,:,:) = dhs_bry(:,:,:,:)
!         ffrac_rest(:,:,:,:) = ffrac_bry(:,:,:,:)
      endif !tr_pond_lvl

!      if (oceanmixed_ice) then
!         frzmelt_rest(:,:,:) = frzmelt_bry(:,:,:)
!         sst_rest(:,:,:) = sst_bry(:,:,:)
!      endif !oceanmixed_ice
                   
      if (tr_FY) then
!         frz_onset_rest(:,:,:) = frz_onset_bry(:,:,:)
         trcrn_rest(:,:,nt_FY,:,:) = FY_bry(:,:,:,:)
      endif !tr_FY

      if (tr_lvl) then
         trcrn_rest(:,:,nt_alvl,:,:) = alvl_bry(:,:,:,:)
         trcrn_rest(:,:,nt_vlvl,:,:) = vlvl_bry(:,:,:,:)
      endif !tr_lvl

      if (tr_iage) then
         trcrn_rest(:,:,nt_iage,:,:) = iage_bry(:,:,:,:)
      endif !tr_iage

      uvel_rest(:,:,:) = uvel_bry(:,:,:)
      vvel_rest(:,:,:) = vvel_bry(:,:,:)
!      scale_factor_rest(:,:,:) = scale_factor_bry(:,:,:)
!      swvdr_rest(:,:,:) = swvdr_bry(:,:,:)
!      swvdf_rest(:,:,:) = swvdf_bry(:,:,:)
!      swidr_rest(:,:,:) = swidr_bry(:,:,:)
!      swidf_rest(:,:,:) = swidf_bry(:,:,:)
!      strocnxT_rest(:,:,:) = strocnxT_bry(:,:,:)
!      strocnyT_rest(:,:,:) = strocnyT_bry(:,:,:)
!      stressp_1_rest(:,:,:) = stressp_1_bry(:,:,:)
!      stressp_2_rest(:,:,:) = stressp_2_bry(:,:,:)
!      stressp_3_rest(:,:,:) = stressp_3_bry(:,:,:)
!      stressp_4_rest(:,:,:) = stressp_4_bry(:,:,:)
!      stressm_1_rest(:,:,:) = stressm_1_bry(:,:,:)
!      stressm_2_rest(:,:,:) = stressm_2_bry(:,:,:)
!      stressm_3_rest(:,:,:) = stressm_3_bry(:,:,:)
!      stressm_4_rest(:,:,:) = stressm_4_bry(:,:,:)
!      stress12_1_rest(:,:,:) = stress12_1_bry(:,:,:)
!      stress12_2_rest(:,:,:) = stress12_2_bry(:,:,:)
!      stress12_3_rest(:,:,:) = stress12_3_bry(:,:,:)
!      stress12_4_rest(:,:,:) = stress12_4_bry(:,:,:)
!      iceumask_rest(:,:,:) = iceumask_bry(:,:,:)

   end subroutine ice_HaloRestore_getbdy

!=======================================================================
#if (1 == 0)
! initialize restoring variables, based on set_state_var
! this routine assumes boundaries are not cyclic

   subroutine set_restore_var (nx_block, ny_block, &
                               ilo, ihi, jlo, jhi, &
                               iglob,    jglob,    &
                               iblock,   jblock,   &
                               Tair, &
                               Tf,                 &
                               salinz,   Tmltz,    &
                               tmask,    aicen,    &
                               trcrn,    ntrcr,    &
                               vicen,    vsnon)

! authors: E. C. Hunke, LANL

!      use ice_arrays_column, only: hin_max

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ilo, ihi          , & ! physical domain indices
         jlo, jhi          , & !
         iglob(nx_block)   , & ! global indices
         jglob(ny_block)   , & !
         iblock            , & ! block indices
         jblock            , & !
         ntrcr                 ! number of tracers in use

      real (kind=dbl_kind), dimension (:,:), intent(in) :: & ! (nx_block,ny_block)
         Tair    , & ! air temperature  (K)
         Tf          ! freezing temperature (C)

      real (kind=dbl_kind), dimension (:,:,:), intent(in) :: & ! (nx_block,ny_block,nilyr)
         salinz  , & ! initial salinity profile
         Tmltz       ! initial melting temperature profile

      logical (kind=log_kind), dimension (:,:), intent(in) :: & ! (nx_block,ny_block)
         tmask      ! true for ice/ocean cells

      real (kind=dbl_kind), dimension (:,:,:), intent(out) :: & ! (nx_block,ny_block,ncat)
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (:,:,:,:), intent(out) :: & ! (nx_block,ny_block,ntrcr,ncat)
         trcrn     ! ice tracers
                   ! 1: surface temperature of ice/snow (C)

      ! local variables

      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij          , & ! horizontal index, combines i and j loops
         ibc         , & ! ghost cell column or row
         npad        , & ! padding column/row counter
         k           , & ! ice layer index
         n           , & ! thickness category index
         it          , & ! tracer index
         nt_Tsfc     , & !
         nt_fbri     , & !
         nt_qice     , & !
         nt_sice     , & !
         nt_qsno     , & !
         icells          ! number of cells initialized with ice

      logical (kind=log_kind) :: &
         tr_brine

      integer (kind=int_kind), dimension(nx_block*ny_block) :: &
         indxi, indxj    ! compressed indices for cells with restoring

      real (kind=dbl_kind) :: &
         Tsfc, hbar, &
         hsno_init       ! initial snow thickness

      real (kind=dbl_kind), dimension(ncat) :: &
         ainit, hinit    ! initial area, thickness

      real (kind=dbl_kind), dimension(nilyr) :: &
         qin             ! ice enthalpy (J/m3)

      real (kind=dbl_kind), dimension(nslyr) :: &
         qsn             ! snow enthalpy (J/m3)

      character(len=*), parameter :: subname = '(set_restore_var)'

      call icepack_query_tracer_flags(tr_brine_out=tr_brine)
      call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_fbri_out=nt_fbri, &
           nt_qice_out=nt_qice, nt_sice_out=nt_sice, nt_qsno_out=nt_qsno)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      indxi(:) = 0
      indxj(:) = 0

      !-----------------------------------------------------------------
      ! Initialize restoring variables everywhere on grid
      !-----------------------------------------------------------------

      do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block
            aicen(i,j,n) = c0
            vicen(i,j,n) = c0
            vsnon(i,j,n) = c0
            if (tmask(i,j)) then
               trcrn(i,j,nt_Tsfc,n) = Tf(i,j)  ! surface temperature
            else
               trcrn(i,j,nt_Tsfc,n) = c0  ! on land gridcells
            endif
            if (ntrcr >= 2) then
               do it = 2, ntrcr
                  trcrn(i,j,it,n) = c0
               enddo
            endif
            if (tr_brine) trcrn(i,j,nt_fbri,n) = c1
         enddo
         enddo
      enddo

      !-----------------------------------------------------------------
      ! initial area and thickness in ice-occupied restoring cells
      !-----------------------------------------------------------------

      hbar = c2  ! initial ice thickness
      hsno_init = 0.20_dbl_kind ! initial snow thickness (m)
      do n = 1, ncat
         hinit(n) = c0
         ainit(n) = c0
         if (hbar > hin_max(n-1) .and. hbar < hin_max(n)) then
            hinit(n) = hbar
            ainit(n) = 0.95_dbl_kind ! initial ice concentration
         endif
      enddo

      !-----------------------------------------------------------------
      ! Define cells where ice is placed (or other values are used)
      ! Edges using initial values (zero, above) are commented out
      !-----------------------------------------------------------------

      icells = 0
      if (iblock == 1) then              ! west edge
            do j = 1, ny_block
            do i = 1, ilo
               if (tmask(i,j)) then
!               icells = icells + 1
!               indxi(icells) = i
!               indxj(icells) = j
               endif
            enddo
            enddo
      endif

      if (iblock == nblocks_x) then      ! east edge
            ! locate ghost cell column (avoid padding)
            ibc = nx_block
            do i = nx_block, 1, -1
               npad = 0
               if (iglob(i) == 0) then
                  do j = 1, ny_block
                     npad = npad + jglob(j)
                  enddo
               endif
               if (npad /= 0) ibc = ibc - 1
            enddo

            do j = 1, ny_block
            do i = ihi, ibc
               if (tmask(i,j)) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
               endif
            enddo
            enddo
      endif

      if (jblock == 1) then              ! south edge
            do j = 1, jlo
            do i = 1, nx_block
               if (tmask(i,j)) then
!               icells = icells + 1
!               indxi(icells) = i
!               indxj(icells) = j
               endif
            enddo
            enddo
      endif

      if (jblock == nblocks_y) then      ! north edge
            ! locate ghost cell row (avoid padding)
            ibc = ny_block
            do j = ny_block, 1, -1
               npad = 0
               if (jglob(j) == 0) then
                  do i = 1, nx_block
                     npad = npad + iglob(i)
                  enddo
               endif
               if (npad /= 0) ibc = ibc - 1
            enddo

            do j = jhi, ibc
            do i = 1, nx_block
               if (tmask(i,j)) then
!               icells = icells + 1
!               indxi(icells) = i
!               indxj(icells) = j
               endif
            enddo
            enddo
      endif

      !-----------------------------------------------------------------
      ! Set restoring variables
      !-----------------------------------------------------------------

         do n = 1, ncat

            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               ! ice volume, snow volume
               aicen(i,j,n) = ainit(n)
               vicen(i,j,n) = hinit(n) * ainit(n) ! m
               vsnon(i,j,n) = min(aicen(i,j,n)*hsno_init,p2*vicen(i,j,n))

               call icepack_init_trcr(Tair=Tair(i,j),    Tf=Tf(i,j),  &
                                      Sprofile=salinz(i,j,:),         &
                                      Tprofile=Tmltz(i,j,:),          &
                                      Tsfc=Tsfc,                      &
                                      qin=qin(:),        qsn=qsn(:))

               ! surface temperature
               trcrn(i,j,nt_Tsfc,n) = Tsfc ! deg C
               ! ice enthalpy, salinity
               do k = 1, nilyr
                  trcrn(i,j,nt_qice+k-1,n) = qin(k)
                  trcrn(i,j,nt_sice+k-1,n) = salinz(i,j,k)
               enddo
               ! snow enthalpy
               do k = 1, nslyr
                  trcrn(i,j,nt_qsno+k-1,n) = qsn(k)
               enddo               ! nslyr

            enddo               ! ij
         enddo                  ! ncat

         call icepack_warnings_flush(nu_diag)
         if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
            file=__FILE__, line=__LINE__)

   end subroutine set_restore_var
#endif
#if (1 == 0)
!=======================================================================
!  This subroutine is intended for restoring the ice state to desired
!  values in halo cells surrounding the grid.

   subroutine ice_Restore(setfld)

!      use ice_calendar, only: dt

      character (len=*), intent(in), optional :: &
         setfld                ! fields to restore

      ! local variables

      integer (int_kind) :: &
         iblk,i,j,n,nt,      & ! dummy loop indices
         ilo,ihi,jlo,jhi,    & ! beginning and end of physical domain
         ntrcr                 ! number of tracers in use

      type (block) :: &
         this_block  ! block info for current block

      real (dbl_kind) :: &
         secday,             & !
         crestore,           & ! restoring value restoring term, dt/trest
         clovalue,           & ! local value restoring term, 1 - dt/trest
         puny

      character (len=char_len_long) :: &
         lsetfld               ! local fields name

      logical (kind=log_kind) :: &
         l_stop                ! if true, abort model

      integer (kind=int_kind) :: &
         istop, jstop, k       ! indices of grid cell where model aborts

      character(len=*), parameter :: subname = '(ice_Restore)'

      if (.not. restore_ice) return

      if (present(setfld)) then
         lsetfld = trim(setfld)
      else
         lsetfld = 'none'
      endif

      if (lsetfld /= 'all' .and. lsetfld /= 'state' .and. &
          lsetfld /= 'velocity' .and. lsetfld /= 'none') then
         call abort_ice(error_message=subname//' ERROR: setfld option unknown = '//trim(lsetfld), &
            file=__FILE__, line=__LINE__)
      endif

      l_stop = .false.
      call ice_timer_start(timer_bound)

      call icepack_query_parameters(secday_out=secday)
      call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      !-----------------------------------------------------------------------
      !  Initialize
      !-----------------------------------------------------------------------

      ! for now, use same restoring constant as for SST
      if (trestore == c0) then
         crestore = c1
         clovalue = c0
      else
         crestore = max(abs(dt/(trestore*secday)),c1)
         clovalue = c1 - crestore
      endif

      !-----------------------------------------------------------------------
      !  Restore values in cells surrounding the grid
      !-----------------------------------------------------------------------

      if (bdy_origin == 'intern') then

         do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            if (this_block%iblock == 1) then          ! west edge
            if (trim(ew_boundary_type) /= 'cyclic') then
               do n = 1, ncat
               do j = 1, ny_block
               do i = 1, ilo
                  aicen (i,j,n,iblk) = aicen_rest(i,j,n,iblk)*crestore + aicen (i,j,n,iblk)*clovalue
                  vicen (i,j,n,iblk) = vicen_rest(i,j,n,iblk)*crestore + vicen (i,j,n,iblk)*clovalue
                  vsnon (i,j,n,iblk) = vsnon_rest(i,j,n,iblk)*crestore + vsnon (i,j,n,iblk)*clovalue
                  do nt = 1, ntrcr
                     trcrn(i,j,nt,n,iblk) = trcrn_rest(i,j,nt,n,iblk)*crestore + trcrn(i,j,nt,n,iblk)*clovalue
                  enddo
               enddo
               enddo
               enddo
            endif
            endif

            if (this_block%iblock == nblocks_x) then  ! east edge
            if (trim(ew_boundary_type) /= 'cyclic') then
               do n = 1, ncat
               do j = 1, ny_block
               do i = ihi, ihi+nghost
                  aicen (i,j,n,iblk) = aicen_rest(i,j,n,iblk)*crestore + aicen (i,j,n,iblk)*clovalue
                  vicen (i,j,n,iblk) = vicen_rest(i,j,n,iblk)*crestore + vicen (i,j,n,iblk)*clovalue
                  vsnon (i,j,n,iblk) = vsnon_rest(i,j,n,iblk)*crestore + vsnon (i,j,n,iblk)*clovalue
                  do nt = 1, ntrcr
                     trcrn(i,j,nt,n,iblk) = trcrn_rest(i,j,nt,n,iblk)*crestore + trcrn(i,j,nt,n,iblk)*clovalue
                  enddo
               enddo
               enddo
               enddo
            endif
            endif

            if (this_block%jblock == 1) then          ! south edge
            if (trim(ns_boundary_type) /= 'cyclic') then
               do n = 1, ncat
               do j = 1, jlo
               do i = 1, nx_block
                  aicen (i,j,n,iblk) = aicen_rest(i,j,n,iblk)*crestore + aicen (i,j,n,iblk)*clovalue
                  vicen (i,j,n,iblk) = vicen_rest(i,j,n,iblk)*crestore + vicen (i,j,n,iblk)*clovalue
                  vsnon (i,j,n,iblk) = vsnon_rest(i,j,n,iblk)*crestore + vsnon (i,j,n,iblk)*clovalue
                  do nt = 1, ntrcr
                     trcrn(i,j,nt,n,iblk) = trcrn_rest(i,j,nt,n,iblk)*crestore + trcrn(i,j,nt,n,iblk)*clovalue
                  enddo
               enddo
               enddo
               enddo
            endif
            endif

            if (this_block%jblock == nblocks_y) then  ! north edge
            if (trim(ns_boundary_type) /= 'cyclic' .and. &
                trim(ns_boundary_type) /= 'tripole' .and. &
                trim(ns_boundary_type) /= 'tripoleT') then

               do n = 1, ncat
               do j = jhi, jhi+nghost
               do i = 1, nx_block
                  aicen (i,j,n,iblk) = aicen_rest(i,j,n,iblk)*crestore + aicen (i,j,n,iblk)*clovalue
                  vicen (i,j,n,iblk) = vicen_rest(i,j,n,iblk)*crestore + vicen (i,j,n,iblk)*clovalue
                  vsnon (i,j,n,iblk) = vsnon_rest(i,j,n,iblk)*crestore + vsnon (i,j,n,iblk)*clovalue
                  do nt = 1, ntrcr
                     trcrn(i,j,nt,n,iblk) = trcrn_rest(i,j,nt,n,iblk)*crestore + trcrn(i,j,nt,n,iblk)*clovalue
                  enddo
               enddo
               enddo
               enddo
            endif
            endif
         enddo ! iblk

      elseif (bdy_origin == 'restart_f') then

         do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            if (trim(ew_boundary_type) /= 'cyclic') then
               ! West Edge
               if (this_block%iblock == 1) then
                  call restore_cells(iblk, 1, ilo-1+nfact, 1, ny_block, &
                       east=.false.,north=.false.,crestore=crestore,setfld=lsetfld)
               endif
               ! East Edge
               if (this_block%iblock == nblocks_x) then
                  call restore_cells(iblk, ihi+1-nfact, ihi+nghost, 1, ny_block, &
                       east=.true.,north=.false.,crestore=crestore,setfld=lsetfld)
               endif
            endif

            ! South Edge
            if (trim(ns_boundary_type) /= 'cyclic') then
               if (this_block%jblock == 1) then
                  call restore_cells(iblk, 1, nx_block, 1, jlo-1+nfact, &
                       east=.false.,north=.false.,crestore=crestore,setfld=lsetfld)
               endif
            endif

            ! North Edge
            if (trim(ns_boundary_type) /= 'cyclic' .and. &
                trim(ns_boundary_type) /= 'tripole' .and. &
                trim(ns_boundary_type) /= 'tripoleT') then
               if (this_block%jblock == nblocks_y) then
                  call restore_cells(iblk, 1, nx_block, jhi+1-nfact, jhi+nghost, &
                       east=.false.,north=.true.,crestore=crestore,setfld=lsetfld)
               endif
            endif

         enddo ! iblk

      else

         call abort_ice(error_message=subname//' ERROR: bdy_origin unknown = '//trim(bdy_origin), &
            file=__FILE__, line=__LINE__)

      endif

      call ice_timer_stop(timer_bound)

   end subroutine ice_Restore
#endif
!=======================================================================
!  This subroutine is intended for restoring the ice state to desired
!  values in halo cells surrounding the grid.

   subroutine ice_HaloRestore(setfld)

      character (len=*), intent(in), optional :: &
         setfld                ! field to restore

      ! local variables

      integer (int_kind) :: &
         iblk,i,j,n,k,nt,    & ! dummy loop indices
         ilo,ihi,jlo,jhi,    & ! beginning and end of physical domain
         ntrcr                 ! number of tracers in use

      type (block) :: &
         this_block  ! block info for current block

      real (dbl_kind) :: &
         crestore              ! restoring value restoring term, dt/trest

      character (len=char_len) :: &
         lsetfld, &            ! local field name
         psetfld               ! passed field name

      logical (log_kind), save :: &
         fldflag(max_set_boundary_flds)  ! matched 

      logical (log_kind), save :: &
         first_call=.true.               ! initialize fldflag

      character(len=*), parameter :: subname = '(ice_HaloRestore)'

      !-----------------------------------------------------------------------
      ! Return conditions to speed up model
      !-----------------------------------------------------------------------

      ! return if no fields set by user
      if (num_set_boundary_flds == 0) return

      ! ignore cyclic and tripole bcs
      if (trim(ew_boundary_type) == 'cyclic' .and. &
          (trim(ns_boundary_type) == 'cyclic' .or. &
           trim(ns_boundary_type) == 'tripole' .or. &
           trim(ns_boundary_type) == 'tripoleT')) then
         return
      endif

      ! ignore setfld = '' or 'none'
      if (present(setfld)) then
         if (setfld == '') return
         lsetfld = trim(setfld)
      else
         return
      endif

      !-----------------------------------------------------------------------
      ! Manage and document usage
      !-----------------------------------------------------------------------

      if (first_call) then
         first_call = .false.
         do n = 1,max_set_boundary_flds
            if (set_boundary_flds(n) == '' .or. set_boundary_flds(n) == 'none') then
               fldflag(n) = .true.   ! don't check if '' or 'none'
            else
               fldflag(n) = .false.
            endif
         enddo
      endif

      !-----------------------------------------------------------------------
      ! Set outer boundary halo values
      !-----------------------------------------------------------------------

      ! set "restoring" to external value
      crestore = c1

      do n = 1,num_set_boundary_flds
         ! Match halo restore call flds to those set by user
         psetfld = 'none'
         if (lsetfld == set_boundary_flds(n)) then
            psetfld = lsetfld
         elseif (lsetfld == 'state' .and. set_boundary_flds(n) == 'aicen') then
            psetfld = 'aicen'
         elseif (lsetfld == 'state' .and. set_boundary_flds(n) == 'vicen') then
            psetfld = 'vicen'
         elseif (lsetfld == 'state' .and. set_boundary_flds(n) == 'vsnon') then
            psetfld = 'vsnon'
         elseif (lsetfld == 'state' .and. set_boundary_flds(n) == 'trcrn') then
            psetfld = 'trcrn'
         endif

         if (psetfld /= 'none') then

            if (.not.fldflag(n)) then
               if (my_task == master_task) &
                              write(nu_diag,*) subname,' setting halo field '//trim(psetfld)
               fldflag(n) = .true.
               do k = 1,num_set_boundary_flds
                  if (.not.fldflag(k)) then
                     if (my_task == master_task) &
                              write(nu_diag,*) subname,' waiting to set '//trim(set_boundary_flds(k))
                  endif
               enddo
            endif

            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi

               if (trim(ew_boundary_type) /= 'cyclic') then
                  ! West Edge
                  if (this_block%iblock == 1) then
                     call restore_cells(iblk, 1, ilo-1+nfact, 1, ny_block, &
                          east=.false.,north=.false.,crestore=crestore,setfld=psetfld)
                  endif
                  ! East Edge
                  if (this_block%iblock == nblocks_x) then
                     call restore_cells(iblk, ihi+1-nfact, ihi+nghost, 1, ny_block, &
                          east=.true.,north=.false.,crestore=crestore,setfld=psetfld)
                  endif
               endif

               ! South Edge
               if (trim(ns_boundary_type) /= 'cyclic') then
                  if (this_block%jblock == 1) then
                     call restore_cells(iblk, 1, nx_block, 1, jlo-1+nfact, &
                          east=.false.,north=.false.,crestore=crestore,setfld=psetfld)
                  endif
               endif

               ! North Edge
               if (trim(ns_boundary_type) /= 'cyclic' .and. &
                   trim(ns_boundary_type) /= 'tripole' .and. &
                   trim(ns_boundary_type) /= 'tripoleT') then
                  if (this_block%jblock == nblocks_y) then
                     call restore_cells(iblk, 1, nx_block, jhi+1-nfact, jhi+nghost, &
                          east=.false.,north=.true.,crestore=crestore,setfld=psetfld)
                  endif
               endif

            enddo ! iblk

         endif
      enddo ! n

   end subroutine ice_HaloRestore

!=======================================================================

   subroutine restore_cells(iblk,i1,i2,j1,j2,east,north,crestore,setfld)

      integer(kind=int_kind), intent(in) :: &
         iblk,      & ! block id
         i1, i2,    & ! i start and end indices
         j1, j2       ! j start and end indices

      logical (kind=log_kind), intent(in) :: &
         east, north  ! are these east or north edges

      real (kind=dbl_kind), intent(in) :: &
         crestore     ! restoring weight in restoring

      character(len=*), intent(in) :: &
         setfld       ! fields to restore

      ! local variable

      integer(kind=int_kind) :: &
         i, j, n, nt, ntrcr, & ! local indices
         i1v, i2v,  & ! i start and end indices for velocity points
         j1v, j2v     ! j start and end indices for velocity points

      real(kind=dbl_kind) :: &
         clovalue     ! local weight in restoring

      character(len=*), parameter :: subname = '(ice_HaloRestore)'

      if (crestore == c0 .or. setfld == 'none') return

      clovalue = c1 - crestore

      call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      i1v = i1
      i2v = i2
      j1v = j1
      j2v = j2
      if (east ) i1v = i1 - 1
      if (north) j1v = j1 - 1

      if (setfld == 'state') then

         ! center gridcell
         do i = i1,i2
         do j = j1,j2
            do n = 1, ncat
               aicen (i,j,n,iblk) = aicen_rest(i,j,n,iblk)*crestore + aicen (i,j,n,iblk)*clovalue
               vicen (i,j,n,iblk) = vicen_rest(i,j,n,iblk)*crestore + vicen (i,j,n,iblk)*clovalue
               vsnon (i,j,n,iblk) = vsnon_rest(i,j,n,iblk)*crestore + vsnon (i,j,n,iblk)*clovalue
               do nt = 1, ntrcr
                  trcrn(i,j,nt,n,iblk) = trcrn_rest(i,j,nt,n,iblk)*crestore + trcrn(i,j,nt,n,iblk)*clovalue
               enddo
            enddo
         enddo
         enddo

      elseif (setfld == 'aicen') then

         ! center gridcell
         do i = i1,i2
         do j = j1,j2
            do n = 1, ncat
               aicen (i,j,n,iblk) = aicen_rest(i,j,n,iblk)*crestore + aicen (i,j,n,iblk)*clovalue
            enddo
         enddo
         enddo

      elseif (setfld == 'vicen') then

         ! center gridcell
         do i = i1,i2
         do j = j1,j2
            do n = 1, ncat
               vicen (i,j,n,iblk) = vicen_rest(i,j,n,iblk)*crestore + vicen (i,j,n,iblk)*clovalue
            enddo
         enddo
         enddo

      elseif (setfld == 'vsnon') then

         ! center gridcell
         do i = i1,i2
         do j = j1,j2
            do n = 1, ncat
               vsnon (i,j,n,iblk) = vsnon_rest(i,j,n,iblk)*crestore + vsnon (i,j,n,iblk)*clovalue
            enddo
         enddo
         enddo

      elseif (setfld == 'trcrn') then

         ! center gridcell
         do i = i1,i2
         do j = j1,j2
            do n = 1, ncat
               do nt = 1, ntrcr
                  trcrn(i,j,nt,n,iblk) = trcrn_rest(i,j,nt,n,iblk)*crestore + trcrn(i,j,nt,n,iblk)*clovalue
               enddo
            enddo
         enddo
         enddo

      elseif (setfld == 'velocity') then

         ! NE gridcell
         do i = i1v, i2v
         do j = j1v, j2v
            uvel(i,j,iblk) = uvel_rest(i,j,iblk)*crestore + uvel(i,j,iblk)*clovalue
            vvel(i,j,iblk) = vvel_rest(i,j,iblk)*crestore + vvel(i,j,iblk)*clovalue
         enddo
         enddo

      elseif (setfld == 'all') then

         ! NE gridcell
         do i = i1v, i2v
         do j = j1v, j2v
            uvel(i,j,iblk) = uvel_rest(i,j,iblk)*crestore + uvel(i,j,iblk)*clovalue
            vvel(i,j,iblk) = vvel_rest(i,j,iblk)*crestore + vvel(i,j,iblk)*clovalue
         enddo
         enddo

         ! center gridcell
         do i = i1,i2
         do j = j1,j2
            do n = 1, ncat
               aicen (i,j,n,iblk) = aicen_rest(i,j,n,iblk)*crestore + aicen (i,j,n,iblk)*clovalue
               vicen (i,j,n,iblk) = vicen_rest(i,j,n,iblk)*crestore + vicen (i,j,n,iblk)*clovalue
               vsnon (i,j,n,iblk) = vsnon_rest(i,j,n,iblk)*crestore + vsnon (i,j,n,iblk)*clovalue
!               dhsn  (i,j,n,iblk) = dhs_rest  (i,j,n,iblk)*crestore + dhsn  (i,j,n,iblk)*clovalue
!               ffracn(i,j,n,iblk) = ffrac_rest(i,j,n,iblk)*crestore + ffracn(i,j,n,iblk)*clovalue
               do nt = 1, ntrcr
                  trcrn(i,j,nt,n,iblk) = trcrn_rest(i,j,nt,n,iblk)*crestore + trcrn(i,j,nt,n,iblk)*clovalue
               enddo
            enddo
!            swvdr(i,j,iblk) = swvdr_rest(i,j,iblk)*crestore + swvdr(i,j,iblk)*clovalue
!            swvdf(i,j,iblk) = swvdf_rest(i,j,iblk)*crestore + swvdf(i,j,iblk)*clovalue
!            swidr(i,j,iblk) = swidr_rest(i,j,iblk)*crestore + swidr(i,j,iblk)*clovalue
!            swidf(i,j,iblk) = swidf_rest(i,j,iblk)*crestore + swidf(i,j,iblk)*clovalue
!            scale_factor(i,j,iblk) = scale_factor_rest(i,j,iblk)*crestore + scale_factor(i,j,iblk)*clovalue
!            strocnxT_iavg(i,j,iblk) = strocnxT_rest(i,j,iblk)*crestore + strocnxT_iavg(i,j,iblk)*clovalue
!            strocnyT_iavg(i,j,iblk) = strocnyT_rest(i,j,iblk)*crestore + strocnyT_iavg(i,j,iblk)*clovalue
!            stressp_1 (i,j,iblk) = stressp_1_rest (i,j,iblk)*crestore + stressp_1 (i,j,iblk)*clovalue
!            stressp_2 (i,j,iblk) = stressp_2_rest (i,j,iblk)*crestore + stressp_2 (i,j,iblk)*clovalue
!            stressp_3 (i,j,iblk) = stressp_3_rest (i,j,iblk)*crestore + stressp_3 (i,j,iblk)*clovalue
!            stressp_4 (i,j,iblk) = stressp_4_rest (i,j,iblk)*crestore + stressp_4 (i,j,iblk)*clovalue
!            stressm_1 (i,j,iblk) = stressm_1_rest (i,j,iblk)*crestore + stressm_1 (i,j,iblk)*clovalue
!            stressm_2 (i,j,iblk) = stressm_2_rest (i,j,iblk)*crestore + stressm_2 (i,j,iblk)*clovalue
!            stressm_3 (i,j,iblk) = stressm_3_rest (i,j,iblk)*crestore + stressm_3 (i,j,iblk)*clovalue
!            stressm_4 (i,j,iblk) = stressm_4_rest (i,j,iblk)*crestore + stressm_4 (i,j,iblk)*clovalue
!            stress12_1(i,j,iblk) = stress12_1_rest(i,j,iblk)*crestore + stress12_1(i,j,iblk)*clovalue
!            stress12_2(i,j,iblk) = stress12_2_rest(i,j,iblk)*crestore + stress12_2(i,j,iblk)*clovalue
!            stress12_3(i,j,iblk) = stress12_3_rest(i,j,iblk)*crestore + stress12_3(i,j,iblk)*clovalue
!            stress12_4(i,j,iblk) = stress12_4_rest(i,j,iblk)*crestore + stress12_4(i,j,iblk)*clovalue
!   !         iceUmask  (i,j,iblk) = iceumask_rest  (i,j,iblk)*crestore + iceUmask  (i,j,iblk)*clovalue
!            frz_onset (i,j,iblk) = frz_onset_rest (i,j,iblk)*crestore + frz_onset (i,j,iblk)*clovalue
!            fsnow     (i,j,iblk) = fsnow_rest     (i,j,iblk)*crestore + fsnow     (i,j,iblk)*clovalue
!            sst       (i,j,iblk) = sst_rest       (i,j,iblk)*crestore + sst       (i,j,iblk)*clovalue
!            frzmlt    (i,j,iblk) = frzmelt_rest   (i,j,iblk)*crestore + frzmlt    (i,j,iblk)*clovalue
         enddo
         enddo

      else

         call abort_ice(error_message=subname//' ERROR: setfld unknown = '//trim(setfld), &
            file=__FILE__, line=__LINE__)

      endif

   end subroutine restore_cells
!=======================================================================

   end module ice_restoring

!=======================================================================
