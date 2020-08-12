!=======================================================================
!
!  Contains CICE component driver routines common to all drivers.
!
!  authors Elizabeth C. Hunke, LANL
!          Philip W. Jones, LANL
!          William H. Lipscomb, LANL
!
! 2008 ECH: created module by moving subroutines from drivers/cice4/
! 2014 ECH: created column package

      module ice_step_mod

      use ice_kinds_mod
      use ice_constants, only: c0, c1, c1000, c4
      use ice_exit, only: abort_ice
      use ice_fileunits, only: nu_diag
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_prep_radiation
      use icepack_intfc, only: icepack_step_therm1
      use icepack_intfc, only: icepack_step_therm2
      use icepack_intfc, only: icepack_aggregate
      use icepack_intfc, only: icepack_step_ridge
      use icepack_intfc, only: icepack_step_wavefracture
      use icepack_intfc, only: icepack_step_radiation
      use icepack_intfc, only: icepack_ocn_mixed_layer, icepack_atm_boundary
      use icepack_intfc, only: icepack_biogeochemistry, icepack_load_ocean_bio_array
      use icepack_intfc, only: icepack_max_algae, icepack_max_nbtrcr, icepack_max_don
      use icepack_intfc, only: icepack_max_doc, icepack_max_dic, icepack_max_aero
      use icepack_intfc, only: icepack_max_fe, icepack_max_iso
      use icepack_intfc, only: icepack_query_parameters
      use icepack_intfc, only: icepack_query_tracer_flags, icepack_query_tracer_sizes
      use icepack_intfc, only: icepack_query_tracer_indices

      implicit none
      private

      public :: step_therm1, step_therm2, step_dyn_horiz, step_dyn_ridge, &
                prep_radiation, step_radiation, ocean_mixed_layer, &
                update_state, biogeochemistry, save_init, step_dyn_wave

!=======================================================================

      contains

!=======================================================================

      subroutine save_init
! saves initial values for aice, aicen, vicen, vsnon

      use ice_state, only: aice, aicen, aice_init, aicen_init, &
          vicen, vicen_init, vsnon, vsnon_init

      !-----------------------------------------------------------------
      ! Save the ice area passed to the coupler (so that history fields
      !  can be made consistent with coupler fields).
      ! Save the initial ice area and volume in each category.
      !-----------------------------------------------------------------

          aice_init = aice
         aicen_init = aicen
         vicen_init = vicen
         vsnon_init = vsnon

      end subroutine save_init

!=======================================================================
!
! Scales radiation fields computed on the previous time step.
!
! authors: Elizabeth Hunke, LANL

      subroutine prep_radiation (iblk)

      use ice_blocks, only: block, get_block
      use ice_domain, only: blocks_ice
      use ice_domain_size, only: ncat, nilyr, nslyr
      use ice_flux, only: scale_factor, swvdr, swvdf, swidr, swidf, &
          alvdr_ai, alvdf_ai, alidr_ai, alidf_ai, &
          alvdr_init, alvdf_init, alidr_init, alidf_init
      use ice_arrays_column, only: fswsfcn, fswintn, &
           fswthrun, fswthrun_vdr, fswthrun_vdf, fswthrun_idr, fswthrun_idf, &
           fswpenln, Sswabsn, Iswabsn
      use ice_state, only: aice, aicen
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_sw

      integer (kind=int_kind), intent(in) :: &
         iblk    ! block index

      ! local variables

      integer (kind=int_kind) :: &
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j               ! horizontal indices

      type (block) :: &
         this_block      ! block information for current block

      character(len=*), parameter :: subname = '(prep_radiation)'

      call ice_timer_start(timer_sw)      ! shortwave

      alvdr_init(:,:,iblk) = c0
      alvdf_init(:,:,iblk) = c0
      alidr_init(:,:,iblk) = c0
      alidf_init(:,:,iblk) = c0

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

      !-----------------------------------------------------------------
      ! Compute netsw scaling factor (new netsw / old netsw)
      !-----------------------------------------------------------------

         do j = jlo, jhi
         do i = ilo, ihi

            alvdr_init(i,j,iblk) = alvdr_ai(i,j,iblk)
            alvdf_init(i,j,iblk) = alvdf_ai(i,j,iblk)
            alidr_init(i,j,iblk) = alidr_ai(i,j,iblk)
            alidf_init(i,j,iblk) = alidf_ai(i,j,iblk)

            call icepack_prep_radiation (ncat=ncat, nilyr=nilyr, nslyr=nslyr,                 &
                        scale_factor=scale_factor(i,j,iblk),                                  &
                        aice     = aice    (i,j,    iblk), aicen    = aicen   (i,j,  :,iblk), &
                        swvdr    = swvdr   (i,j,    iblk), swvdf    = swvdf   (i,j,    iblk), &
                        swidr    = swidr   (i,j,    iblk), swidf    = swidf   (i,j,    iblk), &
                        alvdr_ai = alvdr_ai(i,j,    iblk), alvdf_ai = alvdf_ai(i,j,    iblk), &
                        alidr_ai = alidr_ai(i,j,    iblk), alidf_ai = alidf_ai(i,j,    iblk), &
                        fswsfcn  = fswsfcn (i,j,  :,iblk), fswintn  = fswintn (i,j,  :,iblk), &
                        fswthrun = fswthrun(i,j,  :,iblk), &
                        fswthrun_vdr = fswthrun_vdr(i,j,  :,iblk), &
                        fswthrun_vdf = fswthrun_vdf(i,j,  :,iblk), &
                        fswthrun_idr = fswthrun_idr(i,j,  :,iblk), &
                        fswthrun_idf = fswthrun_idf(i,j,  :,iblk), &
                        fswpenln = fswpenln(i,j,:,:,iblk), &
                        Sswabsn  = Sswabsn (i,j,:,:,iblk), Iswabsn  = Iswabsn (i,j,:,:,iblk))

         enddo               ! i
         enddo               ! j

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      call ice_timer_stop(timer_sw)     ! shortwave

      end subroutine prep_radiation

!=======================================================================
!
! Driver for updating ice and snow internal temperatures and
! computing thermodynamic growth rates and coupler fluxes.
!
! authors: William H. Lipscomb, LANL

      subroutine step_therm1 (dt, iblk)

      use ice_arrays_column, only: ffracn, dhsn, &
          Cdn_ocn, Cdn_ocn_skin, Cdn_ocn_floe, Cdn_ocn_keel, Cdn_atm_ratio, &
          Cdn_atm, Cdn_atm_skin, Cdn_atm_floe, Cdn_atm_rdg, Cdn_atm_pond, &
          hfreebd, hdraft, hridge, distrdg, hkeel, dkeel, lfloe, dfloe, &
          fswsfcn, fswintn, Sswabsn, Iswabsn, &
          fswthrun, fswthrun_vdr, fswthrun_vdf, fswthrun_idr, fswthrun_idf
      use ice_blocks, only: block, get_block, nx_block, ny_block
      use ice_calendar, only: yday
      use ice_domain, only: blocks_ice
      use ice_domain_size, only: ncat, nilyr, nslyr, n_iso, n_aero
      use ice_flux, only: frzmlt, sst, Tf, strocnxT, strocnyT, rside, fbot, Tbot, Tsnice, &
          meltsn, melttn, meltbn, congeln, snoicen, uatm, vatm, fside, &
          wind, rhoa, potT, Qa, zlvl, strax, stray, flatn, fsensn, fsurfn, fcondtopn, &
          flw, fsnow, fpond, sss, mlt_onset, frz_onset, fcondbotn, fcondbot, &
          frain, Tair, strairxT, strairyT, fsurf, fcondtop, fsens, &
          flat, fswabs, flwout, evap, evaps, evapi, Tref, Qref, Uref, fresh, fsalt, fhocn, &
          fswthru, fswthru_vdr, fswthru_vdf, fswthru_idr, fswthru_idf, &
          meltt, melts, meltb, congel, snoice, &
          flatn_f, fsensn_f, fsurfn_f, fcondtopn_f, &
          send_i2x_per_cat, fswthrun_ai
      use ice_flux_bgc, only: dsnown, faero_atm, faero_ocn, fiso_atm, fiso_ocn, &
          Qa_iso, Qref_iso, fiso_evap, HDO_ocn, H2_16O_ocn, H2_18O_ocn
      use ice_grid, only: lmask_n, lmask_s, tmask
      use ice_state, only: aice, aicen, aice_init, aicen_init, vicen_init, &
          vice, vicen, vsno, vsnon, trcrn, uvel, vvel, vsnon_init

#ifdef CESMCOUPLED
      use ice_prescribed_mod, only: prescribed_ice
#else
      logical (kind=log_kind) :: & 
         prescribed_ice ! if .true., use prescribed ice instead of computed
#endif
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      integer (kind=int_kind), intent(in) :: &
         iblk    ! block index

      ! local variables
#ifdef CICE_IN_NEMO
      real (kind=dbl_kind)    :: & 
         raice              ! temporary reverse ice concentration
#endif
      integer (kind=int_kind) :: &
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j           , & ! horizontal indices
         n              , & ! thickness category index
         k, kk              ! indices for aerosols

      integer (kind=int_kind) :: &
         ntrcr, nt_apnd, nt_hpnd, nt_ipnd, nt_alvl, nt_vlvl, nt_Tsfc, &
         nt_iage, nt_FY, nt_qice, nt_sice, nt_aero, nt_qsno, &
         nt_isosno, nt_isoice

      logical (kind=log_kind) :: &
         tr_iage, tr_FY, tr_iso, tr_aero, tr_pond, tr_pond_cesm, &
         tr_pond_lvl, tr_pond_topo, calc_Tsfc

      real (kind=dbl_kind) :: &
         puny

      real (kind=dbl_kind), dimension(n_aero,2,ncat) :: &
         aerosno,  aeroice    ! kg/m^2

      real (kind=dbl_kind), dimension(n_iso,ncat) :: &
         isosno,  isoice      ! kg/m^2

      type (block) :: &
         this_block      ! block information for current block

      character(len=*), parameter :: subname = '(step_therm1)'

      call icepack_query_parameters(puny_out=puny)
      call icepack_query_parameters(calc_Tsfc_out=calc_Tsfc)
      call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
      call icepack_query_tracer_flags( &
         tr_iage_out=tr_iage, tr_FY_out=tr_FY, tr_iso_out=tr_iso, &
         tr_aero_out=tr_aero, tr_pond_out=tr_pond, tr_pond_cesm_out=tr_pond_cesm, &
         tr_pond_lvl_out=tr_pond_lvl, tr_pond_topo_out=tr_pond_topo)
      call icepack_query_tracer_indices( &
         nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, nt_ipnd_out=nt_ipnd, &
         nt_alvl_out=nt_alvl, nt_vlvl_out=nt_vlvl, nt_Tsfc_out=nt_Tsfc, &
         nt_iage_out=nt_iage, nt_FY_out=nt_FY, &
         nt_qice_out=nt_qice, nt_sice_out=nt_sice, &
         nt_aero_out=nt_aero, nt_qsno_out=nt_qsno, &
         nt_isosno_out=nt_isosno, nt_isoice_out=nt_isoice)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

#ifndef CESMCOUPLED
      prescribed_ice = .false.
#endif

      isosno (:,:)   = c0
      isoice (:,:)   = c0
      aerosno(:,:,:) = c0
      aeroice(:,:,:) = c0

#ifdef CICE_IN_NEMO
      do j = 1, ny_block
      do i = 1, nx_block

      !---------------------------------------------------------------
      ! Scale frain and fsnow by ice concentration as these fields
      ! are supplied by NEMO multiplied by ice concentration
      !---------------------------------------------------------------

         if (aice_init(i,j,iblk) > puny) then
            raice           = c1 / aice_init(i,j,iblk)
            frain(i,j,iblk) = frain(i,j,iblk)*raice
            fsnow(i,j,iblk) = fsnow(i,j,iblk)*raice
         else
            frain(i,j,iblk) = c0
            fsnow(i,j,iblk) = c0
         endif

      enddo ! i
      enddo ! j
#endif

      this_block = get_block(blocks_ice(iblk),iblk)         
      ilo = this_block%ilo
      ihi = this_block%ihi
      jlo = this_block%jlo
      jhi = this_block%jhi
      
      do j = jlo, jhi
      do i = ilo, ihi

         if (tr_iso) then ! trcrn(nt_iso*) has units kg/m^3
            do n=1,ncat
               do k=1,n_iso
                  isosno(k,n) = trcrn(i,j,nt_isosno+k-1,n,iblk) * vsnon_init(i,j,n,iblk)
                  isoice(k,n) = trcrn(i,j,nt_isoice+k-1,n,iblk) * vicen_init(i,j,n,iblk)
               enddo
            enddo
         endif ! tr_iso

         if (tr_aero) then ! trcrn(nt_aero) has units kg/m^3
            do n=1,ncat
               do k=1,n_aero
                  aerosno (k,:,n) = &
                     trcrn(i,j,nt_aero+(k-1)*4  :nt_aero+(k-1)*4+1,n,iblk) &
                                  * vsnon_init(i,j,n,iblk)
                  aeroice (k,:,n) = &
                     trcrn(i,j,nt_aero+(k-1)*4+2:nt_aero+(k-1)*4+3,n,iblk) &
                                  * vicen_init(i,j,n,iblk)
               enddo
            enddo
         endif ! tr_aero

         if (tmask(i,j,iblk)) then

         call icepack_step_therm1(dt=dt, ncat=ncat,            &
                      nilyr=nilyr, nslyr=nslyr,                &
                      aicen_init   = aicen_init  (i,j,:,iblk), &
                      vicen_init   = vicen_init  (i,j,:,iblk), &
                      vsnon_init   = vsnon_init  (i,j,:,iblk), &
                      aice         = aice        (i,j,  iblk), &
                      aicen        = aicen       (i,j,:,iblk), &
                      vice         = vice        (i,j,  iblk), &
                      vicen        = vicen       (i,j,:,iblk), &
                      vsno         = vsno        (i,j,  iblk), &
                      vsnon        = vsnon       (i,j,:,iblk), &
                      uvel         = uvel        (i,j,  iblk), &
                      vvel         = vvel        (i,j,  iblk), &
                      Tsfc         = trcrn       (i,j,nt_Tsfc,:,iblk),                   &
                      zqsn         = trcrn       (i,j,nt_qsno:nt_qsno+nslyr-1,:,iblk),   & 
                      zqin         = trcrn       (i,j,nt_qice:nt_qice+nilyr-1,:,iblk),   & 
                      zSin         = trcrn       (i,j,nt_sice:nt_sice+nilyr-1,:,iblk),   & 
                      alvl         = trcrn       (i,j,nt_alvl,:,iblk),                   & 
                      vlvl         = trcrn       (i,j,nt_vlvl,:,iblk),                   & 
                      apnd         = trcrn       (i,j,nt_apnd,:,iblk),                   & 
                      hpnd         = trcrn       (i,j,nt_hpnd,:,iblk),                   & 
                      ipnd         = trcrn       (i,j,nt_ipnd,:,iblk),                   & 
                      iage         = trcrn       (i,j,nt_iage,:,iblk),                   &
                      FY           = trcrn       (i,j,nt_FY  ,:,iblk),                   & 
                      aerosno      = aerosno     (:,:,:),      &
                      aeroice      = aeroice     (:,:,:),      &
                      isosno       = isosno      (:,:),        &
                      isoice       = isoice      (:,:),        &
                      uatm         = uatm        (i,j,  iblk), &
                      vatm         = vatm        (i,j,  iblk), &
                      wind         = wind        (i,j,  iblk), &
                      zlvl         = zlvl        (i,j,  iblk), &
                      Qa           = Qa          (i,j,  iblk), &
                      Qa_iso       = Qa_iso      (i,j,:,iblk), &
                      rhoa         = rhoa        (i,j,  iblk), &
                      Tair         = Tair        (i,j,  iblk), &
                      Tref         = Tref        (i,j,  iblk), &
                      Qref         = Qref        (i,j,  iblk), &
                      Qref_iso     = Qref_iso    (i,j,:,iblk), &
                      Uref         = Uref        (i,j,  iblk), &
                      Cdn_atm_ratio= Cdn_atm_ratio(i,j, iblk), &
                      Cdn_ocn      = Cdn_ocn     (i,j,  iblk), &
                      Cdn_ocn_skin = Cdn_ocn_skin(i,j,  iblk), &
                      Cdn_ocn_floe = Cdn_ocn_floe(i,j,  iblk), &
                      Cdn_ocn_keel = Cdn_ocn_keel(i,j,  iblk), &
                      Cdn_atm      = Cdn_atm     (i,j,  iblk), &
                      Cdn_atm_skin = Cdn_atm_skin(i,j,  iblk), &
                      Cdn_atm_floe = Cdn_atm_floe(i,j,  iblk), &
                      Cdn_atm_pond = Cdn_atm_pond(i,j,  iblk), &
                      Cdn_atm_rdg  = Cdn_atm_rdg (i,j,  iblk), &
                      hfreebd      = hfreebd     (i,j,  iblk), &
                      hdraft       = hdraft      (i,j,  iblk), &
                      hridge       = hridge      (i,j,  iblk), &
                      distrdg      = distrdg     (i,j,  iblk), &
                      hkeel        = hkeel       (i,j,  iblk), &
                      dkeel        = dkeel       (i,j,  iblk), &
                      lfloe        = lfloe       (i,j,  iblk), &
                      dfloe        = dfloe       (i,j,  iblk), &
                      strax        = strax       (i,j,  iblk), &
                      stray        = stray       (i,j,  iblk), &
                      strairxT     = strairxT    (i,j,  iblk), &
                      strairyT     = strairyT    (i,j,  iblk), &
                      potT         = potT        (i,j,  iblk), &
                      sst          = sst         (i,j,  iblk), &
                      sss          = sss         (i,j,  iblk), &
                      Tf           = Tf          (i,j,  iblk), &
                      strocnxT     = strocnxT    (i,j,  iblk), &
                      strocnyT     = strocnyT    (i,j,  iblk), &
                      fbot         = fbot        (i,j,  iblk), &
                      Tbot         = Tbot        (i,j,  iblk), &
                      Tsnice       = Tsnice       (i,j, iblk), &
                      frzmlt       = frzmlt      (i,j,  iblk), &
                      rside        = rside       (i,j,  iblk), &
                      fside        = fside       (i,j,  iblk), &
                      fsnow        = fsnow       (i,j,  iblk), &
                      frain        = frain       (i,j,  iblk), &
                      fpond        = fpond       (i,j,  iblk), &
                      fsurf        = fsurf       (i,j,  iblk), &
                      fsurfn       = fsurfn      (i,j,:,iblk), &
                      fcondtop     = fcondtop    (i,j,  iblk), &
                      fcondtopn    = fcondtopn   (i,j,:,iblk), &
                      fcondbot     = fcondbot    (i,j,  iblk), &
                      fcondbotn    = fcondbotn   (i,j,:,iblk), &
                      fswsfcn      = fswsfcn     (i,j,:,iblk), &
                      fswintn      = fswintn     (i,j,:,iblk), &
                      fswthrun     = fswthrun    (i,j,:,iblk), &
                      fswthrun_vdr = fswthrun_vdr (i,j,:,iblk),&
                      fswthrun_vdf = fswthrun_vdf (i,j,:,iblk),&
                      fswthrun_idr = fswthrun_idr (i,j,:,iblk),&
                      fswthrun_idf = fswthrun_idf (i,j,:,iblk),&
                      fswabs       = fswabs      (i,j,  iblk), &
                      flwout       = flwout      (i,j,  iblk), &
                      Sswabsn      = Sswabsn     (i,j,:,:,iblk), &
                      Iswabsn      = Iswabsn     (i,j,:,:,iblk), &
                      flw          = flw         (i,j,  iblk), &
                      fsens        = fsens       (i,j,  iblk), &
                      fsensn       = fsensn      (i,j,:,iblk), &
                      flat         = flat        (i,j,  iblk), &
                      flatn        = flatn       (i,j,:,iblk), &
                      evap         = evap        (i,j,  iblk), &
                      evaps        = evaps       (i,j,  iblk), &
                      evapi        = evapi       (i,j,  iblk), &
                      fresh        = fresh       (i,j,  iblk), &
                      fsalt        = fsalt       (i,j,  iblk), &
                      fhocn        = fhocn       (i,j,  iblk), &
                      fswthru      = fswthru     (i,j,  iblk), &
                      fswthru_vdr  = fswthru_vdr  (i,j,  iblk),&
                      fswthru_vdf  = fswthru_vdf  (i,j,  iblk),&
                      fswthru_idr  = fswthru_idr  (i,j,  iblk),&
                      fswthru_idf  = fswthru_idf  (i,j,  iblk),&
                      flatn_f      = flatn_f     (i,j,:,iblk), &
                      fsensn_f     = fsensn_f    (i,j,:,iblk), &
                      fsurfn_f     = fsurfn_f    (i,j,:,iblk), &
                      fcondtopn_f  = fcondtopn_f (i,j,:,iblk), &
                      faero_atm    = faero_atm   (i,j,1:n_aero,iblk), &
                      faero_ocn    = faero_ocn   (i,j,1:n_aero,iblk), &
                      fiso_atm     = fiso_atm    (i,j,:,iblk), &
                      fiso_ocn     = fiso_ocn    (i,j,:,iblk), &
                      fiso_evap    = fiso_evap   (i,j,:,iblk), &
                      HDO_ocn      = HDO_ocn     (i,j,  iblk), &
                      H2_16O_ocn   = H2_16O_ocn  (i,j,  iblk), &
                      H2_18O_ocn   = H2_18O_ocn  (i,j,  iblk), &
                      dhsn         = dhsn        (i,j,:,iblk), &
                      ffracn       = ffracn      (i,j,:,iblk), &
                      meltt        = meltt       (i,j,  iblk), &
                      melttn       = melttn      (i,j,:,iblk), &
                      meltb        = meltb       (i,j,  iblk), &
                      meltbn       = meltbn      (i,j,:,iblk), &
                      melts        = melts       (i,j,  iblk), &
                      meltsn       = meltsn      (i,j,:,iblk), &
                      congel       = congel      (i,j,  iblk), &
                      congeln      = congeln     (i,j,:,iblk), &
                      snoice       = snoice      (i,j,  iblk), &
                      snoicen      = snoicen     (i,j,:,iblk), &
                      dsnown       = dsnown      (i,j,:,iblk), &
                      lmask_n      = lmask_n     (i,j,  iblk), &
                      lmask_s      = lmask_s     (i,j,  iblk), &
                      mlt_onset    = mlt_onset   (i,j,  iblk), &
                      frz_onset    = frz_onset   (i,j,  iblk), &
                      yday=yday, prescribed_ice=prescribed_ice)

      !-----------------------------------------------------------------
      ! handle per-category i2x fields, no merging
      !-----------------------------------------------------------------

         if (send_i2x_per_cat) then
            do n = 1, ncat
               ! TODO (mvertens, 2018-12-22): do we need to add the band separated quantities
               ! for MOM6 here also?

               fswthrun_ai(i,j,n,iblk) = fswthrun(i,j,n,iblk)*aicen_init(i,j,n,iblk)
            enddo                  ! ncat
         endif

         endif

         if (tr_iso) then
            do n = 1, ncat
               if (vicen(i,j,n,iblk) > puny) &
                  isoice(:,n) = isoice(:,n)/vicen(i,j,n,iblk)
               if (vsnon(i,j,n,iblk) > puny) &
                  isosno(:,n) = isosno(:,n)/vsnon(i,j,n,iblk)
               do k = 1, n_iso
                  trcrn(i,j,nt_isosno+k-1,n,iblk) = isosno(k,n)
                  trcrn(i,j,nt_isoice+k-1,n,iblk) = isoice(k,n)
               enddo
            enddo
         endif ! tr_iso

         if (tr_aero) then
            do n = 1, ncat
               if (vicen(i,j,n,iblk) > puny) &
                  aeroice(:,:,n) = aeroice(:,:,n)/vicen(i,j,n,iblk)
               if (vsnon(i,j,n,iblk) > puny) &
                  aerosno(:,:,n) = aerosno(:,:,n)/vsnon(i,j,n,iblk)
               do k = 1, n_aero
                  do kk = 1, 2
                     trcrn(i,j,nt_aero+(k-1)*4+kk-1,n,iblk)=aerosno(k,kk,n)
                     trcrn(i,j,nt_aero+(k-1)*4+kk+1,n,iblk)=aeroice(k,kk,n)
                  enddo
               enddo
            enddo
         endif ! tr_aero

      enddo ! i
      enddo ! j

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      end subroutine step_therm1

!=======================================================================
! Driver for thermodynamic changes not needed for coupling:
! transport in thickness space, lateral growth and melting.
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL

      subroutine step_therm2 (dt, iblk)

      use ice_arrays_column, only: hin_max, fzsal, ocean_bio, wave_sig_ht, &
          wave_spectrum, wavefreq, dwavefreq, &
          first_ice, bgrid, cgrid, igrid, floe_rad_c, floe_binwidth, &
          d_afsd_latg, d_afsd_newi, d_afsd_latm, d_afsd_weld
      use ice_blocks, only: block, get_block
      use ice_calendar, only: yday
      use ice_domain, only: blocks_ice
      use ice_domain_size, only: ncat, nilyr, nslyr, n_aero, nblyr, nfsd
      use ice_flux, only: fresh, frain, fpond, frzmlt, frazil, frz_onset, &
          update_ocn_f, fsalt, Tf, sss, salinz, fhocn, rside, fside, &
          meltl, frazil_diag
      use ice_flux_bgc, only: flux_bio, faero_ocn, &
          fiso_ocn, HDO_ocn, H2_16O_ocn, H2_18O_ocn
      use ice_grid, only: tmask
      use ice_state, only: aice, aicen, aice0, trcr_depend, &
          aicen_init, vicen_init, trcrn, vicen, vsnon, &
          trcr_base, n_trcr_strata, nt_strata

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      integer (kind=int_kind), intent(in) :: &
         iblk    ! block index

      ! local variables

      integer (kind=int_kind) :: &
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j               ! horizontal indices

      integer (kind=int_kind) :: &
         ntrcr, nbtrcr, nltrcr

      logical (kind=log_kind) :: &
         tr_fsd,          & ! floe size distribution tracers
         z_tracers

      type (block) :: &
         this_block         ! block information for current block

      character(len=*), parameter :: subname = '(step_therm2)'

      call icepack_query_parameters(z_tracers_out=z_tracers)
      call icepack_query_tracer_sizes(ntrcr_out=ntrcr, nbtrcr_out=nbtrcr)
      call icepack_query_tracer_flags(tr_fsd_out=tr_fsd)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      ! tcraig, nltrcr used to be the number of zbgc tracers, but it's used as a zbgc flag in icepack
      if (z_tracers) then
         nltrcr = 1
      else
         nltrcr = 0
      endif

      this_block = get_block(blocks_ice(iblk),iblk)         
      ilo = this_block%ilo
      ihi = this_block%ihi
      jlo = this_block%jlo
      jhi = this_block%jhi

      do j = jlo, jhi
      do i = ilo, ihi

         if (tmask(i,j,iblk)) then

         ! significant wave height for FSD
         if (tr_fsd) &
         wave_sig_ht(i,j,iblk) = c4*SQRT(SUM(wave_spectrum(i,j,:,iblk)*dwavefreq(:)))

         call icepack_step_therm2(dt=dt, ncat=ncat, &
                      nltrcr=nltrcr, nilyr=nilyr, nslyr=nslyr, nblyr=nblyr, &
                      hin_max    = hin_max   (:),          &   
                      aicen      = aicen     (i,j,:,iblk), &
                      vicen      = vicen     (i,j,:,iblk), &
                      vsnon      = vsnon     (i,j,:,iblk), &
                      aicen_init = aicen_init(i,j,:,iblk), &
                      vicen_init = vicen_init(i,j,:,iblk), &
                      trcrn      = trcrn     (i,j,:,:,iblk), &
                      aice0      = aice0     (i,j,  iblk), &
                      aice       = aice      (i,j,  iblk), &
                      trcr_depend= trcr_depend(:),         &
                      trcr_base  = trcr_base(:,:),         &
                      n_trcr_strata = n_trcr_strata(:),    &
                      nt_strata  = nt_strata(:,:),         &
                      Tf         = Tf        (i,j,  iblk), &
                      sss        = sss       (i,j,  iblk), &
                      salinz     = salinz    (i,j,:,iblk), &
                      rside      = rside     (i,j,  iblk), &
                      meltl      = meltl     (i,j,  iblk), &
                      fside      = fside     (i,j,  iblk), &
                      frzmlt     = frzmlt    (i,j,  iblk), &
                      frazil     = frazil    (i,j,  iblk), &
                      frain      = frain     (i,j,  iblk), &
                      fpond      = fpond     (i,j,  iblk), &
                      fresh      = fresh     (i,j,  iblk), &
                      fsalt      = fsalt     (i,j,  iblk), &
                      fhocn      = fhocn     (i,j,  iblk), &
                      update_ocn_f = update_ocn_f,         &
                      bgrid      = bgrid,                  &
                      cgrid      = cgrid,                  &
                      igrid      = igrid,                  &
                      faero_ocn  = faero_ocn (i,j,:,iblk), &
                      first_ice  = first_ice (i,j,:,iblk), &
                      fzsal      = fzsal     (i,j,  iblk), &
                      flux_bio   = flux_bio  (i,j,1:nbtrcr,iblk), &
                      ocean_bio  = ocean_bio (i,j,1:nbtrcr,iblk), &
                      frazil_diag= frazil_diag(i,j,iblk),  &
                      frz_onset  = frz_onset (i,j,  iblk), &
                      yday       = yday,                   &
                      fiso_ocn   = fiso_ocn  (i,j,:,iblk), &
                      HDO_ocn    = HDO_ocn   (i,j,  iblk), &
                      H2_16O_ocn = H2_16O_ocn(i,j,  iblk), &
                      H2_18O_ocn = H2_18O_ocn(i,j,  iblk), &
                      nfsd       = nfsd,                   &
                      wave_sig_ht= wave_sig_ht(i,j,iblk),  &
                      wave_spectrum = wave_spectrum(i,j,:,iblk),  &
                      wavefreq   = wavefreq(:),            &
                      dwavefreq  = dwavefreq(:),           &
                      d_afsd_latg= d_afsd_latg(i,j,:,iblk),&
                      d_afsd_newi= d_afsd_newi(i,j,:,iblk),&
                      d_afsd_latm= d_afsd_latm(i,j,:,iblk),&
                      d_afsd_weld= d_afsd_weld(i,j,:,iblk),&
                      floe_rad_c = floe_rad_c(:),          &
                      floe_binwidth = floe_binwidth(:))
         endif ! tmask

      enddo                     ! i
      enddo                     ! j

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      end subroutine step_therm2

!=======================================================================
!
! finalize thermo updates
!
! authors: Elizabeth Hunke, LANL

      subroutine update_state (dt, daidt, dvidt, dagedt, offset)

      use ice_blocks, only: nx_block, ny_block
      use ice_domain, only: nblocks
      use ice_domain_size, only: ncat
!     use ice_grid, only: tmask
      use ice_state, only: aicen, trcrn, vicen, vsnon, &
                           aice,  trcr,  vice,  vsno, aice0, trcr_depend, &
                           bound_state, trcr_base, nt_strata, n_trcr_strata
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bound

      real (kind=dbl_kind), intent(in) :: &
         dt    , & ! time step
         offset    ! d(age)/dt time offset = dt for thermo, 0 for dyn

      real (kind=dbl_kind), dimension(:,:,:), intent(inout) :: &
          daidt, & ! change in ice area per time step
          dvidt, & ! change in ice volume per time step
          dagedt   ! change in ice age per time step

      integer (kind=int_kind) :: & 
         iblk,  & ! block index 
         i,j,   & ! horizontal indices
         ntrcr, & !
         nt_iage  !

      logical (kind=log_kind) :: &
         tr_iage  !

      character(len=*), parameter :: subname='(update_state)'

      call icepack_query_tracer_flags(tr_iage_out=tr_iage)
      call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
      call icepack_query_tracer_indices(nt_iage_out=nt_iage)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      !-------------------------------------------------------------------
      ! Ghost cell updates for state variables.
      !-------------------------------------------------------------------

      call ice_timer_start(timer_bound)
      call bound_state (aicen,        &
                        vicen, vsnon, &
                        ntrcr, trcrn)
      call ice_timer_stop(timer_bound)

      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block

      !-----------------------------------------------------------------
      ! Aggregate the updated state variables (includes ghost cells). 
      !----------------------------------------------------------------- 
 
!        if (tmask(i,j,iblk)) &
            call icepack_aggregate(ncat  = ncat,                  &
                                   aicen = aicen(i,j,:,iblk),     &
                                   trcrn = trcrn(i,j,:,:,iblk),   &
                                   vicen = vicen(i,j,:,iblk),     &
                                   vsnon = vsnon(i,j,:,iblk),     &
                                   aice  = aice (i,j,  iblk),     &
                                   trcr  = trcr (i,j,:,iblk),     &
                                   vice  = vice (i,j,  iblk),     &
                                   vsno  = vsno (i,j,  iblk),     &
                                   aice0 = aice0(i,j,  iblk),     &
                                   ntrcr = ntrcr,                 &
                                   trcr_depend   = trcr_depend(:),   &
                                   trcr_base     = trcr_base(:,:),   &
                                   n_trcr_strata = n_trcr_strata(:), &
                                   nt_strata     = nt_strata(:,:))

      !-----------------------------------------------------------------
      ! Compute thermodynamic area and volume tendencies.
      !-----------------------------------------------------------------

         daidt(i,j,iblk) = (aice(i,j,iblk) - daidt(i,j,iblk)) / dt
         dvidt(i,j,iblk) = (vice(i,j,iblk) - dvidt(i,j,iblk)) / dt
         if (tr_iage) then
            if (offset > c0) then                 ! thermo
               if (trcr(i,j,nt_iage,iblk) > c0) &
               dagedt(i,j,iblk) = (trcr(i,j,nt_iage,iblk) &
                                - dagedt(i,j,iblk) - offset) / dt
            else                                  ! dynamics
               dagedt(i,j,iblk) = (trcr(i,j,nt_iage,iblk) &
                                - dagedt(i,j,iblk)) / dt
            endif
         endif

         enddo ! i
         enddo ! j
      enddo    ! iblk
      !$OMP END PARALLEL DO

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      end subroutine update_state

!=======================================================================
!
! Run one time step of wave-fracturing the floe size distribution
!
! authors: Lettie Roach, NIWA
!          Elizabeth C. Hunke, LANL

      subroutine step_dyn_wave (dt)

      use ice_arrays_column, only: wave_spectrum, wave_sig_ht, &
          d_afsd_wave, floe_rad_l, floe_rad_c, wavefreq, dwavefreq
      use ice_blocks, only: block, get_block
      use ice_domain, only: blocks_ice, nblocks
      use ice_domain_size, only: ncat, nfsd, nfreq
      use ice_state, only: trcrn, aicen, aice, vice
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_column, &
          timer_fsd

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      type (block) :: &
         this_block      ! block information for current block

      integer (kind=int_kind) :: &
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         iblk,            & ! block index
         i, j,            & ! horizontal indices
         ntrcr,           & !
         nbtrcr             !

      character (len=char_len) :: wave_spec_type

      character(len=*), parameter :: subname = '(step_dyn_wave)'

      call ice_timer_start(timer_column)
      call ice_timer_start(timer_fsd)

      call icepack_query_parameters(wave_spec_type_out=wave_spec_type)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks

         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            d_afsd_wave(i,j,:,iblk) = c0
            call icepack_step_wavefracture (wave_spec_type, &
                                            dt, ncat, nfsd, nfreq,         &
                                            aice           (i,j,    iblk), &
                                            vice           (i,j,    iblk), &
                                            aicen          (i,j,:,  iblk), &
                                            floe_rad_l(:), floe_rad_c(:),  &
                                            wave_spectrum  (i,j,:,  iblk), &
                                            wavefreq(:),   dwavefreq(:),   &
                                            trcrn          (i,j,:,:,iblk), &
                                            d_afsd_wave    (i,j,:,  iblk))
         end do ! i
         end do ! j
      end do    ! iblk
      !$OMP END PARALLEL DO

      call ice_timer_stop(timer_fsd)
      call ice_timer_stop(timer_column)

      end subroutine step_dyn_wave

!=======================================================================
!
! Run one time step of dynamics and horizontal transport.
! NOTE: The evp and transport modules include boundary updates, so
!       they cannot be done inside a single block loop.
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL

      subroutine step_dyn_horiz (dt)

      use ice_dyn_evp, only: evp
      use ice_dyn_eap, only: eap
      use ice_dyn_shared, only: kdyn, ktransport
      use ice_flux, only: init_history_dyn
      use ice_transport_driver, only: advection, transport_upwind, transport_remap

      real (kind=dbl_kind), intent(in) :: &
         dt      ! dynamics time step

      character(len=*), parameter :: subname = '(step_dyn_horiz)'

      call init_history_dyn     ! initialize dynamic history variables

      !-----------------------------------------------------------------
      ! Elastic-viscous-plastic ice dynamics
      !-----------------------------------------------------------------

      if (kdyn == 1) call evp (dt)
      if (kdyn == 2) call eap (dt)

      !-----------------------------------------------------------------
      ! Horizontal ice transport
      !-----------------------------------------------------------------

      if (ktransport > 0) then
      if (advection == 'upwind') then
         call transport_upwind (dt)    ! upwind
      else
         call transport_remap (dt)     ! incremental remapping
      endif
      endif

      end subroutine step_dyn_horiz

!=======================================================================
!
! Run one time step of ridging.
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL

      subroutine step_dyn_ridge (dt, ndtd, iblk)

      use ice_arrays_column, only: hin_max, fzsal, first_ice
      use ice_blocks, only: block, get_block
      use ice_domain, only: blocks_ice
      use ice_domain_size, only: ncat, nilyr, nslyr, n_aero, nblyr
      use ice_flux, only: &
          rdg_conv, rdg_shear, dardg1dt, dardg2dt, &
          dvirdgdt, opening, fpond, fresh, fhocn, &
          aparticn, krdgn, aredistn, vredistn, dardg1ndt, dardg2ndt, &
          dvirdgndt, araftn, vraftn, fsalt
      use ice_flux_bgc, only: flux_bio, faero_ocn, fiso_ocn
      use ice_grid, only: tmask
      use ice_state, only: trcrn, vsnon, aicen, vicen, &
          aice, aice0, trcr_depend, n_trcr_strata, &
          trcr_base, nt_strata
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_column, &
          timer_ridge

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      integer (kind=int_kind), intent(in) :: &
         ndtd, & ! number of dynamics subcycles
         iblk    ! block index 

      ! local variables

      type (block) :: &
         this_block      ! block information for current block

      integer (kind=int_kind) :: & 
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j,            & ! horizontal indices
         ntrcr,           & !
         nbtrcr             !

      character(len=*), parameter :: subname = '(step_dyn_ridge)'

      !-----------------------------------------------------------------
      ! Ridging
      !-----------------------------------------------------------------

      call ice_timer_start(timer_column)
      call ice_timer_start(timer_ridge)

      call icepack_query_tracer_sizes(ntrcr_out=ntrcr, nbtrcr_out=nbtrcr)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      this_block = get_block(blocks_ice(iblk), iblk)
      ilo = this_block%ilo
      ihi = this_block%ihi
      jlo = this_block%jlo
      jhi = this_block%jhi

      do j = jlo, jhi
      do i = ilo, ihi

!echmod: this changes the answers, continue using tmask for now
!      call aggregate_area (ncat, aicen(i,j,:,iblk), atmp, atmp0)
!      if (atmp > c0) then

         if (tmask(i,j,iblk)) then

            call icepack_step_ridge (dt=dt, ndtd=ndtd,                 &
                         nilyr=nilyr, nslyr=nslyr, nblyr=nblyr,        &
                         ncat=ncat, n_aero=n_aero, hin_max=hin_max(:), &
                         trcr_depend   = trcr_depend  (:),   &
                         trcr_base     = trcr_base    (:,:), &
                         n_trcr_strata = n_trcr_strata(:),   &
                         nt_strata     = nt_strata    (:,:), &
                         trcrn     = trcrn    (i,j,:,:,iblk), &
                         rdg_conv  = rdg_conv (i,j,  iblk), &
                         rdg_shear = rdg_shear(i,j,  iblk), &
                         aicen     = aicen    (i,j,:,iblk), &
                         vicen     = vicen    (i,j,:,iblk), &
                         vsnon     = vsnon    (i,j,:,iblk), &
                         aice0     = aice0    (i,j,  iblk), &
                         dardg1dt  = dardg1dt (i,j,  iblk), &
                         dardg2dt  = dardg2dt (i,j,  iblk), &
                         dvirdgdt  = dvirdgdt (i,j,  iblk), &
                         opening   = opening  (i,j,  iblk), &
                         fpond     = fpond    (i,j,  iblk), &
                         fresh     = fresh    (i,j,  iblk), &
                         fhocn     = fhocn    (i,j,  iblk), &
                         faero_ocn = faero_ocn(i,j,:,iblk), &
                         fiso_ocn  = fiso_ocn (i,j,:,iblk), &
                         aparticn  = aparticn (i,j,:,iblk), &
                         krdgn     = krdgn    (i,j,:,iblk), &
                         aredistn  = aredistn (i,j,:,iblk), &
                         vredistn  = vredistn (i,j,:,iblk), &
                         dardg1ndt = dardg1ndt(i,j,:,iblk), &
                         dardg2ndt = dardg2ndt(i,j,:,iblk), &
                         dvirdgndt = dvirdgndt(i,j,:,iblk), &
                         araftn    = araftn   (i,j,:,iblk), &
                         vraftn    = vraftn   (i,j,:,iblk), &
                         aice      = aice     (i,j,  iblk), &
                         fsalt     = fsalt    (i,j,  iblk), &
                         first_ice = first_ice(i,j,:,iblk), &
                         fzsal     = fzsal    (i,j,  iblk), &
                         flux_bio  = flux_bio (i,j,1:nbtrcr,iblk))

         endif ! tmask

      enddo ! i
      enddo ! j

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      call ice_timer_stop(timer_ridge)
      call ice_timer_stop(timer_column)

      end subroutine step_dyn_ridge

!=======================================================================
!
! Computes radiation fields
!
! authors: William H. Lipscomb, LANL
!          David Bailey, NCAR
!          Elizabeth C. Hunke, LANL

      subroutine step_radiation (dt, iblk)

      use ice_arrays_column, only: ffracn, dhsn, &
          fswsfcn, fswintn, fswpenln, Sswabsn, Iswabsn, &
          fswthrun, fswthrun_vdr, fswthrun_vdf, fswthrun_idr, fswthrun_idf, &
          albicen, albsnon, albpndn, &
          alvdrn, alidrn, alvdfn, alidfn, apeffn, trcrn_sw, snowfracn, &
          kaer_tab, waer_tab, gaer_tab, kaer_bc_tab, waer_bc_tab, &
          gaer_bc_tab, bcenh, swgrid, igrid
      use ice_blocks, only: block, get_block
      use ice_calendar, only: calendar_type, days_per_year, nextsw_cday, yday, sec
      use ice_domain, only: blocks_ice
      use ice_domain_size, only: ncat, n_aero, nilyr, nslyr, n_zaero, n_algae, nblyr
      use ice_flux, only: swvdr, swvdf, swidr, swidf, coszen, fsnow
      use ice_grid, only: TLAT, TLON, tmask
      use ice_state, only: aicen, vicen, vsnon, trcrn
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_sw
      use ice_communicate, only: my_task
      use ice_diagnostics, only: npnt, print_points, pmloc, piloc, pjloc

      real (kind=dbl_kind), intent(in) :: &
         dt                 ! time step

      integer (kind=int_kind), intent(in) :: &
         iblk               ! block index

      ! local variables

      integer (kind=int_kind) :: &
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, n,   k,    & ! horizontal indices
         ipoint             ! index for print diagnostic

      type (block) :: &
         this_block         ! block information for current block

      integer (kind=int_kind) :: &
         nt_Tsfc, nt_alvl, &
         nt_apnd, nt_hpnd, nt_ipnd, nt_aero, nlt_chl_sw, &
         ntrcr, nbtrcr, nbtrcr_sw, nt_fbri

      integer (kind=int_kind), dimension(icepack_max_algae) :: &
         nt_bgc_N

      integer (kind=int_kind), dimension(icepack_max_aero) :: &
         nlt_zaero_sw, nt_zaero

      logical (kind=log_kind) :: &
         tr_bgc_N, tr_zaero, tr_brine, dEdd_algae, modal_aero

      real (kind=dbl_kind), dimension(ncat) :: &
         fbri                 ! brine height to ice thickness

      real(kind= dbl_kind), dimension(:,:), allocatable :: &
         ztrcr_sw

      logical (kind=log_kind) :: &
         debug, &           ! flag for printing debugging information
         l_print_point      ! flag for printing debugging information

      character(len=*), parameter :: subname = '(step_radiation)'

      call ice_timer_start(timer_sw)      ! shortwave

      call icepack_query_tracer_sizes(ntrcr_out=ntrcr, &
         nbtrcr_out=nbtrcr, nbtrcr_sw_out=nbtrcr_sw)
      call icepack_query_tracer_flags( &
         tr_brine_out=tr_brine, tr_bgc_N_out=tr_bgc_N, tr_zaero_out=tr_zaero)
      call icepack_query_tracer_indices( &
         nt_Tsfc_out=nt_Tsfc, nt_alvl_out=nt_alvl, &
         nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, nt_ipnd_out=nt_ipnd, nt_aero_out=nt_aero, &
         nlt_chl_sw_out=nlt_chl_sw, nlt_zaero_sw_out=nlt_zaero_sw, &
         nt_fbri_out=nt_fbri, nt_zaero_out=nt_zaero, nt_bgc_N_out=nt_bgc_N)
      call icepack_query_parameters(dEdd_algae_out=dEdd_algae, modal_aero_out=modal_aero)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      allocate(ztrcr_sw(nbtrcr_sw,ncat))

      this_block = get_block(blocks_ice(iblk),iblk)         
      ilo = this_block%ilo
      ihi = this_block%ihi
      jlo = this_block%jlo
      jhi = this_block%jhi

      do j = jlo, jhi
      do i = ilo, ihi

         l_print_point = .false.
         debug = .false.
         if (debug .and. print_points) then
            do ipoint = 1, npnt
               if (my_task == pmloc(ipoint) .and. &
                    i == piloc(ipoint) .and. &
                    j == pjloc(ipoint)) &
                    l_print_point = .true.
                    write (nu_diag, *) 'my_task = ',my_task
            enddo ! ipoint
         endif
         fbri(:) = c0
         ztrcr_sw(:,:) = c0
         do n = 1, ncat
           if (tr_brine)  fbri(n) = trcrn(i,j,nt_fbri,n,iblk)
         enddo

         if (tmask(i,j,iblk)) then

            call icepack_step_radiation (dt=dt,   ncat=ncat,                  &
                         nblyr=nblyr, nilyr=nilyr, nslyr=nslyr,               &
                         dEdd_algae=dEdd_algae,                               &
                         swgrid=swgrid(:),        igrid=igrid(:),             &
                         fbri=fbri(:),                                        &
                         aicen=aicen(i,j,        :,iblk),                     &
                         vicen=vicen(i,j,        :,iblk),                     &
                         vsnon=vsnon(i,j,        :,iblk),                     &
                         Tsfcn=trcrn(i,j,nt_Tsfc,:,iblk),                     &
                         alvln=trcrn(i,j,nt_alvl,:,iblk),                     &
                         apndn=trcrn(i,j,nt_apnd,:,iblk),                     &
                         hpndn=trcrn(i,j,nt_hpnd,:,iblk),                     &
                         ipndn=trcrn(i,j,nt_ipnd,:,iblk),                     &
                         aeron=trcrn(i,j,nt_aero:nt_aero+4*n_aero-1,:,iblk),  &
                         bgcNn=trcrn(i,j,nt_bgc_N(1):nt_bgc_N(1)+n_algae*(nblyr+3)-1,:,iblk), &
                         zaeron=trcrn(i,j,nt_zaero(1):nt_zaero(1)+n_zaero*(nblyr+3)-1,:,iblk), &
                         trcrn_bgcsw=ztrcr_sw,                                &
                         TLAT=TLAT(i,j,iblk),     TLON=TLON(i,j,iblk),        &
                         calendar_type=calendar_type,                         &
                         days_per_year=days_per_year,                         &
                         nextsw_cday=nextsw_cday, yday=yday,                  &
                         sec=sec,                                             &
                         kaer_tab=kaer_tab, kaer_bc_tab=kaer_bc_tab(:,:),     &
                         waer_tab=waer_tab, waer_bc_tab=waer_bc_tab(:,:),     &
                         gaer_tab=gaer_tab, gaer_bc_tab=gaer_bc_tab(:,:),     &
                         bcenh=bcenh(:,:,:),                                  &
                         modal_aero=modal_aero,                               &
                         swvdr    =swvdr    (i,j    ,iblk), swvdf   =swvdf   (i,j    ,iblk), &
                         swidr    =swidr    (i,j    ,iblk), swidf   =swidf   (i,j    ,iblk), &
                         coszen   =coszen   (i,j    ,iblk), fsnow   =fsnow   (i,j    ,iblk), &
                         alvdrn   =alvdrn   (i,j,:  ,iblk), alvdfn  =alvdfn  (i,j,:  ,iblk), &
                         alidrn   =alidrn   (i,j,:  ,iblk), alidfn  =alidfn  (i,j,:  ,iblk), &
                         fswsfcn  =fswsfcn  (i,j,:  ,iblk), fswintn =fswintn (i,j,:  ,iblk), &
                         fswthrun =fswthrun (i,j,:  ,iblk), &
                         fswthrun_vdr =fswthrun_vdr (i,j,:  ,iblk), &
                         fswthrun_vdf =fswthrun_vdf (i,j,:  ,iblk), &
                         fswthrun_idr =fswthrun_idr (i,j,:  ,iblk), &
                         fswthrun_idf =fswthrun_idf (i,j,:  ,iblk), &
                         fswpenln=fswpenln(i,j,:,:,iblk), &
                         Sswabsn  =Sswabsn  (i,j,:,:,iblk), Iswabsn =Iswabsn (i,j,:,:,iblk), &
                         albicen  =albicen  (i,j,:  ,iblk), albsnon =albsnon (i,j,:  ,iblk), &
                         albpndn  =albpndn  (i,j,:  ,iblk), apeffn  =apeffn  (i,j,:  ,iblk), &
                         snowfracn=snowfracn(i,j,:  ,iblk),                                  &
                         dhsn     =dhsn     (i,j,:  ,iblk), ffracn  =ffracn(i,j,:,iblk),     &
                         l_print_point=l_print_point)

         endif
         
         if (dEdd_algae .and. (tr_zaero .or. tr_bgc_N)) then
           do n = 1, ncat
              do k = 1, nbtrcr_sw
                 trcrn_sw(i,j,k,n,iblk) = ztrcr_sw(k,n)
              enddo
           enddo
         endif

      enddo ! i
      enddo ! j

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      deallocate(ztrcr_sw)

      call ice_timer_stop(timer_sw)     ! shortwave

      end subroutine step_radiation

!=======================================================================
! Ocean mixed layer calculation (internal to sea ice model).
! Allows heat storage in ocean for uncoupled runs.
!
! authors:   John Weatherly, CRREL
!            C.M. Bitz, UW
!            Elizabeth C. Hunke, LANL
!            Bruce P. Briegleb, NCAR
!            William H. Lipscomb, LANL

      subroutine ocean_mixed_layer (dt, iblk)

      use ice_arrays_column, only: Cdn_atm, Cdn_atm_ratio
      use ice_blocks, only: nx_block, ny_block
      use ice_flux, only: sst, Tf, Qa, uatm, vatm, wind, potT, rhoa, zlvl, &
           frzmlt, fhocn, fswthru, flw, flwout_ocn, fsens_ocn, flat_ocn, evap_ocn, &
           alvdr_ocn, alidr_ocn, alvdf_ocn, alidf_ocn, swidf, swvdf, swidr, swvdr, &
           qdp, hmix, strairx_ocn, strairy_ocn, Tref_ocn, Qref_ocn
      use ice_grid, only: tmask
      use ice_state, only: aice

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      integer (kind=int_kind), intent(in) :: &
         iblk    ! block index

      ! local variables

      real (kind=dbl_kind) :: albocn

      real (kind=dbl_kind), parameter :: &
         frzmlt_max = c1000   ! max magnitude of frzmlt (W/m^2)

      integer (kind=int_kind) :: &
         i, j           , & ! horizontal indices
         ij                 ! combined ij index

      real (kind=dbl_kind), dimension(nx_block,ny_block) :: &
         delt  , & ! potential temperature difference   (K)
         delq  , & ! specific humidity difference   (kg/kg)
         shcoef, & ! transfer coefficient for sensible heat
         lhcoef    ! transfer coefficient for latent heat

      integer (kind=int_kind) :: &
         icells    ! number of ocean cells

      integer (kind=int_kind), dimension(nx_block*ny_block) :: &
         indxi, indxj    ! compressed indices for ocean cells

      character(len=*), parameter :: subname = '(ocn_mixed_layer)'

      !-----------------------------------------------------------------

         call icepack_query_parameters(albocn_out=albocn)
         call icepack_warnings_flush(nu_diag)
         if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
            file=__FILE__, line=__LINE__)

      !-----------------------------------------------------------------
      ! Identify ocean cells.
      ! Set fluxes to zero in land cells.
      !-----------------------------------------------------------------

         icells = 0
         indxi(:) = 0
         indxj(:) = 0

         do j = 1, ny_block
         do i = 1, nx_block

            if (tmask(i,j,iblk)) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            else
               sst       (i,j,iblk) = c0
               frzmlt    (i,j,iblk) = c0
               flwout_ocn(i,j,iblk) = c0
               fsens_ocn (i,j,iblk) = c0
               flat_ocn  (i,j,iblk) = c0
               evap_ocn  (i,j,iblk) = c0
            endif
         enddo                  ! i
         enddo                  ! j

      !-----------------------------------------------------------------
      ! Compute boundary layer quantities
      !-----------------------------------------------------------------

         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            call icepack_atm_boundary(sfctype = 'ocn',    &
                         Tsf     = sst        (i,j,iblk), &    
                         potT    = potT       (i,j,iblk), &
                         uatm    = uatm       (i,j,iblk), &   
                         vatm    = vatm       (i,j,iblk), &   
                         wind    = wind       (i,j,iblk), &   
                         zlvl    = zlvl       (i,j,iblk), &   
                         Qa      = Qa         (i,j,iblk), &     
                         rhoa    = rhoa       (i,j,iblk), &
                         strx    = strairx_ocn(i,j,iblk), & 
                         stry    = strairy_ocn(i,j,iblk), & 
                         Tref    = Tref_ocn   (i,j,iblk), & 
                         Qref    = Qref_ocn   (i,j,iblk), & 
                         delt    = delt       (i,j),      &    
                         delq    = delq       (i,j),      &
                         lhcoef  = lhcoef     (i,j),      &
                         shcoef  = shcoef     (i,j),      &
                         Cdn_atm = Cdn_atm    (i,j,iblk), & 
                         Cdn_atm_ratio_n = Cdn_atm_ratio(i,j,iblk))    
         enddo ! ij

         call icepack_warnings_flush(nu_diag)
         if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
            file=__FILE__, line=__LINE__)

      !-----------------------------------------------------------------
      ! Ocean albedo
      ! For now, assume albedo = albocn in each spectral band.
      !-----------------------------------------------------------------

         alvdr_ocn(:,:,iblk) = albocn
         alidr_ocn(:,:,iblk) = albocn
         alvdf_ocn(:,:,iblk) = albocn
         alidf_ocn(:,:,iblk) = albocn

      !-----------------------------------------------------------------
      ! Compute ocean fluxes and update SST
      !-----------------------------------------------------------------

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         call icepack_ocn_mixed_layer(alvdr_ocn=alvdr_ocn(i,j,iblk), swvdr =swvdr (i,j,iblk), &
                                      alidr_ocn=alidr_ocn(i,j,iblk), swidr =swidr (i,j,iblk), &
                                      alvdf_ocn=alvdf_ocn(i,j,iblk), swvdf =swvdf (i,j,iblk), &
                                      alidf_ocn=alidf_ocn(i,j,iblk), swidf =swidf (i,j,iblk), &
                                      sst      =sst      (i,j,iblk), flwout_ocn=flwout_ocn(i,j,iblk), &
                                      fsens_ocn=fsens_ocn(i,j,iblk), shcoef=shcoef(i,j), &
                                      flat_ocn =flat_ocn (i,j,iblk), lhcoef=lhcoef(i,j), &
                                      evap_ocn =evap_ocn (i,j,iblk), flw   =flw   (i,j,iblk), &
                                      delt     =delt     (i,j),      delq  =delq  (i,j), &
                                      aice     =aice     (i,j,iblk), fhocn =fhocn (i,j,iblk), &
                                      fswthru  =fswthru  (i,j,iblk), hmix  =hmix  (i,j,iblk), &
                                      Tf       =Tf       (i,j,iblk), qdp   =qdp   (i,j,iblk), &
                                      frzmlt   =frzmlt   (i,j,iblk), dt    =dt)
      enddo                    ! ij

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      end subroutine ocean_mixed_layer

!=======================================================================

      subroutine biogeochemistry (dt, iblk)

      use ice_arrays_column, only: upNO, upNH, iDi, iki, zfswin, &
                           zsal_tot, darcy_V, grow_net,  &
                           PP_net, hbri,dhbr_bot, dhbr_top, Zoo,&
                           fbio_snoice, fbio_atmice, ocean_bio,  &
                           first_ice, fswpenln, bphi, bTiz, ice_bio_net,  &
                           snow_bio_net, fswthrun, Rayleigh_criteria, &
                           ocean_bio_all, sice_rho, fzsal, fzsal_g, &
                           bgrid, igrid, icgrid, cgrid
      use ice_blocks, only: block, get_block
      use ice_domain, only: blocks_ice
      use ice_domain_size, only: nblyr, nilyr, nslyr, n_algae, n_zaero, ncat, &
                                 n_doc, n_dic,  n_don, n_fed, n_fep
      use ice_flux, only: meltbn, melttn, congeln, snoicen, &
                          sst, sss, fsnow, meltsn
      use ice_flux_bgc, only: hin_old, flux_bio, flux_bio_atm, faero_atm, & 
          nit, amm, sil, dmsp, dms, algalN, doc, don, dic, fed, fep, zaeros, hum
      use ice_state, only: aicen_init, vicen_init, aicen, vicen, vsnon, &
          trcrn, vsnon_init, aice0                    
      use ice_timers, only: timer_bgc, ice_timer_start, ice_timer_stop

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      integer (kind=int_kind), intent(in) :: &
         iblk    ! block index

      ! local variables

      integer (kind=int_kind) :: &
         i, j           , & ! horizontal indices
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         mm              ! tracer index

      type (block) :: &
         this_block      ! block information for current block

      integer (kind=int_kind) :: &
         nbtrcr, ntrcr

      integer (kind=int_kind), dimension(icepack_max_aero) :: &
         nlt_zaero

      integer (kind=int_kind), dimension(icepack_max_nbtrcr) :: &
         bio_index_o

      logical (kind=log_kind) :: &
         skl_bgc, tr_brine, tr_zaero

      character(len=*), parameter :: subname='(biogeochemistry)'

      call icepack_query_tracer_flags(tr_brine_out=tr_brine)
      call icepack_query_parameters(skl_bgc_out=skl_bgc)
      call icepack_query_tracer_sizes(ntrcr_out=ntrcr, nbtrcr_out=nbtrcr)
      call icepack_query_tracer_flags(tr_zaero_out=tr_zaero)
      call icepack_query_tracer_indices(nlt_zaero_out=nlt_zaero)
      call icepack_query_tracer_indices(bio_index_o_out=bio_index_o)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (tr_brine .or. skl_bgc) then

      call ice_timer_start(timer_bgc) ! biogeochemistry

      this_block = get_block(blocks_ice(iblk),iblk)         
      ilo = this_block%ilo
      ihi = this_block%ihi
      jlo = this_block%jlo
      jhi = this_block%jhi

      ! Define ocean concentrations for tracers used in simulation
      do j = jlo, jhi
      do i = ilo, ihi    

         call icepack_load_ocean_bio_array(max_nbtrcr = icepack_max_nbtrcr, &
                max_algae = icepack_max_algae, max_don = icepack_max_don, &
                max_doc   = icepack_max_doc,   max_dic = icepack_max_dic, &
                max_aero  = icepack_max_aero,  max_fe  = icepack_max_fe,  &
                nit = nit(i,j,  iblk), amm    = amm   (i,j,  iblk), &
                sil = sil(i,j,  iblk), dmsp   = dmsp  (i,j,  iblk), &
                dms = dms(i,j,  iblk), algalN = algalN(i,j,:,iblk), &
                doc = doc(i,j,:,iblk), don    = don   (i,j,:,iblk), &
                dic = dic(i,j,:,iblk), fed    = fed   (i,j,:,iblk), &
                fep = fep(i,j,:,iblk), zaeros = zaeros(i,j,:,iblk), &
                hum = hum(i,j,  iblk),                              &
                ocean_bio_all = ocean_bio_all(i,j,:,iblk))

         do mm = 1,nbtrcr
            ocean_bio(i,j,mm,iblk) = ocean_bio_all(i,j,bio_index_o(mm),iblk)  
         enddo  ! mm    
         if (tr_zaero) then
            do mm = 1, n_zaero  ! update aerosols
               flux_bio_atm(i,j,nlt_zaero(mm),iblk) = faero_atm(i,j,mm,iblk)
            enddo  ! mm
         endif

         call icepack_biogeochemistry(dt=dt, ntrcr=ntrcr, nbtrcr=nbtrcr,&
                              bgrid=bgrid, igrid=igrid, icgrid=icgrid, cgrid=cgrid,             &
                              nblyr=nblyr, nilyr=nilyr, nslyr=nslyr, n_algae=n_algae, n_zaero=n_zaero,   &
                              ncat=ncat, n_doc=n_doc, n_dic=n_dic, n_don=n_don, n_fed=n_fed, n_fep=n_fep, &
                              upNO         = upNO        (i,j,          iblk), &
                              upNH         = upNH        (i,j,          iblk), &
                              iDi          = iDi         (i,j,:,:,      iblk), &
                              iki          = iki         (i,j,:,:,      iblk), &
                              zfswin       = zfswin      (i,j,:,:,      iblk), &
                              zsal_tot     = zsal_tot    (i,j,          iblk), &
                              darcy_V      = darcy_V     (i,j,:,        iblk), &
                              grow_net     = grow_net    (i,j,          iblk), &
                              PP_net       = PP_net      (i,j,          iblk), &
                              hbri         = hbri        (i,j,          iblk), &
                              dhbr_bot     = dhbr_bot    (i,j,:,        iblk), &
                              dhbr_top     = dhbr_top    (i,j,:,        iblk), &
                              Zoo          = Zoo         (i,j,:,:,      iblk), &
                              fbio_snoice  = fbio_snoice (i,j,:,        iblk), &
                              fbio_atmice  = fbio_atmice (i,j,:,        iblk), &
                              ocean_bio    = ocean_bio   (i,j,1:nbtrcr, iblk), &
                              first_ice    = first_ice   (i,j,:,        iblk), &
                              fswpenln     = fswpenln    (i,j,:,:,      iblk), &
                              bphi         = bphi        (i,j,:,:,      iblk), &
                              bTiz         = bTiz        (i,j,:,:,      iblk), &
                              ice_bio_net  = ice_bio_net (i,j,1:nbtrcr, iblk), &
                              snow_bio_net = snow_bio_net(i,j,1:nbtrcr, iblk), &
                              fswthrun     = fswthrun    (i,j,:,        iblk), &
                              sice_rho     = sice_rho    (i,j,:,        iblk), &
                              fzsal        = fzsal       (i,j,          iblk), &   
                              fzsal_g      = fzsal_g     (i,j,          iblk), &
                              meltbn       = meltbn      (i,j,:,        iblk), &
                              melttn       = melttn      (i,j,:,        iblk), &
                              congeln      = congeln     (i,j,:,        iblk), &
                              snoicen      = snoicen     (i,j,:,        iblk), & 
                              sst          = sst         (i,j,          iblk), &    
                              sss          = sss         (i,j,          iblk), &
                              fsnow        = fsnow       (i,j,          iblk), &
                              meltsn       = meltsn      (i,j,:,        iblk), &
                              hin_old      = hin_old     (i,j,:,        iblk), &
                              flux_bio     = flux_bio    (i,j,1:nbtrcr, iblk), &
                              flux_bio_atm = flux_bio_atm(i,j,1:nbtrcr, iblk), &
                              aicen_init   = aicen_init  (i,j,:,        iblk), &
                              vicen_init   = vicen_init  (i,j,:,        iblk), &
                              aicen        = aicen       (i,j,:,        iblk), &
                              vicen        = vicen       (i,j,:,        iblk), &
                              vsnon        = vsnon       (i,j,:,        iblk), &
                              aice0        = aice0       (i,j,          iblk), &
                              trcrn        = trcrn       (i,j,:,:,      iblk), &
                              vsnon_init   = vsnon_init  (i,j,:,        iblk), &
                              Rayleigh_criteria = Rayleigh_criteria(i,j,iblk), &
                              skl_bgc      = skl_bgc)

      enddo               ! i
      enddo               ! j

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      call ice_timer_stop(timer_bgc) ! biogeochemistry

      endif  ! tr_brine .or. skl_bgc

      end subroutine biogeochemistry

!=======================================================================

      end module ice_step_mod

!=======================================================================
