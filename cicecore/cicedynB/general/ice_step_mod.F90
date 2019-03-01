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
      use ice_constants, only: c0, c1000
      use ice_exit, only: abort_ice
      use ice_fileunits, only: nu_diag
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_prep_radiation
      use icepack_intfc, only: icepack_step_therm1
      use icepack_intfc, only: icepack_step_therm2
      use icepack_intfc, only: icepack_aggregate
      use icepack_intfc, only: icepack_step_ridge
      use icepack_intfc, only: icepack_step_radiation
      use icepack_intfc, only: icepack_ocn_mixed_layer, icepack_atm_boundary
      use icepack_intfc, only: icepack_biogeochemistry, icepack_init_OceanConcArray
      use icepack_intfc, only: icepack_max_algae, icepack_max_nbtrcr, icepack_max_don
      use icepack_intfc, only: icepack_max_doc, icepack_max_dic, icepack_max_aero
      use icepack_intfc, only: icepack_max_fe
      use icepack_intfc, only: icepack_query_parameters
      use icepack_intfc, only: icepack_query_tracer_flags, icepack_query_tracer_numbers
      use icepack_intfc, only: icepack_query_tracer_indices

      implicit none
      private

      public :: step_therm1, step_therm2, step_dyn_horiz, step_dyn_ridge, &
                prep_radiation, step_radiation, ocean_mixed_layer, &
                update_state, biogeochemistry, save_init

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
      use ice_arrays_column, only: fswsfcn, fswintn, fswthrun, &
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

      alvdr_init(:,:,:) = c0
      alvdf_init(:,:,:) = c0
      alidr_init(:,:,:) = c0
      alidf_init(:,:,:) = c0

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

            call icepack_prep_radiation (ncat, nilyr, nslyr,             &
                        aice    (i,j,    iblk), aicen   (i,j,  :,iblk), &
                        swvdr   (i,j,    iblk), swvdf   (i,j,    iblk), &
                        swidr   (i,j,    iblk), swidf   (i,j,    iblk), &
                        alvdr_ai(i,j,    iblk), alvdf_ai(i,j,    iblk), &
                        alidr_ai(i,j,    iblk), alidf_ai(i,j,    iblk), &
                        scale_factor(i,j,iblk),                         &
                        fswsfcn (i,j,  :,iblk), fswintn (i,j,  :,iblk), &
                        fswthrun(i,j,  :,iblk), fswpenln(i,j,:,:,iblk), &
                        Sswabsn (i,j,:,:,iblk), Iswabsn (i,j,:,:,iblk))

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
          fswsfcn, fswintn, fswthrun, Sswabsn, Iswabsn
      use ice_blocks, only: block, get_block, nx_block, ny_block
      use ice_calendar, only: yday
      use ice_domain, only: blocks_ice
      use ice_domain_size, only: ncat, nilyr, nslyr, n_aero
      use ice_flux, only: frzmlt, sst, Tf, strocnxT, strocnyT, rside, fbot, Tbot, Tsnice, &
          meltsn, melttn, meltbn, congeln, snoicen, uatm, vatm, &
          wind, rhoa, potT, Qa, zlvl, strax, stray, flatn, fsensn, fsurfn, fcondtopn, &
          flw, fsnow, fpond, sss, mlt_onset, frz_onset, fcondbotn, fcondbot, &
          frain, Tair, strairxT, strairyT, fsurf, fcondtop, fsens, &
          flat, fswabs, flwout, evap, evaps, evapi, Tref, Qref, Uref, fresh, fsalt, fhocn, &
          fswthru, meltt, melts, meltb, congel, snoice, &
          flatn_f, fsensn_f, fsurfn_f, fcondtopn_f
      use ice_flux_bgc, only: dsnown, faero_atm, faero_ocn
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

      integer (kind=int_kind) :: &
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j           , & ! horizontal indices
         n              , & ! thickness category index
         k, kk              ! indices for aerosols

      integer (kind=int_kind) :: &
         ntrcr, nt_apnd, nt_hpnd, nt_ipnd, nt_alvl, nt_vlvl, nt_Tsfc, &
         nt_iage, nt_FY, nt_qice, nt_sice, nt_aero, nt_qsno

      logical (kind=log_kind) :: &
         tr_iage, tr_FY, tr_aero, tr_pond, tr_pond_cesm, &
         tr_pond_lvl, tr_pond_topo, calc_Tsfc

      real (kind=dbl_kind) :: &
         puny

      real (kind=dbl_kind), dimension(n_aero,2,ncat) :: &
         aerosno,  aeroice    ! kg/m^2

      type (block) :: &
         this_block      ! block information for current block

      character(len=*), parameter :: subname = '(step_therm1)'

      call icepack_query_parameters(puny_out=puny)
      call icepack_query_parameters(calc_Tsfc_out=calc_Tsfc)
      call icepack_query_tracer_numbers(ntrcr_out=ntrcr)
      call icepack_query_tracer_flags( &
         tr_iage_out=tr_iage, tr_FY_out=tr_FY, &
         tr_aero_out=tr_aero, tr_pond_out=tr_pond, tr_pond_cesm_out=tr_pond_cesm, &
         tr_pond_lvl_out=tr_pond_lvl, tr_pond_topo_out=tr_pond_topo)
      call icepack_query_tracer_indices( &
         nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, nt_ipnd_out=nt_ipnd, &
         nt_alvl_out=nt_alvl, nt_vlvl_out=nt_vlvl, nt_Tsfc_out=nt_Tsfc, &
         nt_iage_out=nt_iage, nt_FY_out=nt_FY, &
         nt_qice_out=nt_qice, nt_sice_out=nt_sice, &
         nt_aero_out=nt_aero, nt_qsno_out=nt_qsno)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

#ifndef CESMCOUPLED
      prescribed_ice = .false.
#endif

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

         if (tr_aero) then
         ! trcrn(nt_aero) has units kg/m^3
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

         if (tmask(i,j,iblk)) &
         call icepack_step_therm1(dt, ncat, nilyr, nslyr, n_aero,                &
                            aicen_init  (i,j,:,iblk),                           &
                            vicen_init  (i,j,:,iblk), vsnon_init  (i,j,:,iblk), &
                            aice        (i,j,  iblk), aicen       (i,j,:,iblk), &
                            vice        (i,j,  iblk), vicen       (i,j,:,iblk), &
                            vsno        (i,j,  iblk), vsnon       (i,j,:,iblk), &
                            uvel        (i,j,  iblk), vvel        (i,j,  iblk), &
                            trcrn       (i,j,nt_Tsfc,:,iblk),                   &
                            trcrn       (i,j,nt_qsno:nt_qsno+nslyr-1,:,iblk),   & 
                            trcrn       (i,j,nt_qice:nt_qice+nilyr-1,:,iblk),   & 
                            trcrn       (i,j,nt_sice:nt_sice+nilyr-1,:,iblk),   & 
                            trcrn       (i,j,nt_alvl,:,iblk),                   & 
                            trcrn       (i,j,nt_vlvl,:,iblk),                   & 
                            trcrn       (i,j,nt_apnd,:,iblk),                   & 
                            trcrn       (i,j,nt_hpnd,:,iblk),                   & 
                            trcrn       (i,j,nt_ipnd,:,iblk),                   & 
                            trcrn       (i,j,nt_iage,:,iblk),                   &
                            trcrn       (i,j,nt_FY  ,:,iblk),                   & 
                            aerosno     (:,:,:),      aeroice     (:,:,:),      &
                            uatm        (i,j,  iblk), vatm        (i,j,  iblk), &
                            wind        (i,j,  iblk), zlvl        (i,j,  iblk), &
                            Qa          (i,j,  iblk), rhoa        (i,j,  iblk), &
                            Tair        (i,j,  iblk), Tref        (i,j,  iblk), &
                            Qref        (i,j,  iblk), Uref        (i,j,  iblk), &
                            Cdn_atm_ratio(i,j, iblk),                           &
                            Cdn_ocn     (i,j,  iblk), Cdn_ocn_skin(i,j,  iblk), &
                            Cdn_ocn_floe(i,j,  iblk), Cdn_ocn_keel(i,j,  iblk), &
                            Cdn_atm     (i,j,  iblk), Cdn_atm_skin(i,j,  iblk), &
                            Cdn_atm_floe(i,j,  iblk), Cdn_atm_pond(i,j,  iblk), &
                            Cdn_atm_rdg (i,j,  iblk), hfreebd     (i,j,  iblk), &
                            hdraft      (i,j,  iblk), hridge      (i,j,  iblk), &
                            distrdg     (i,j,  iblk), hkeel       (i,j,  iblk), &
                            dkeel       (i,j,  iblk), lfloe       (i,j,  iblk), &
                            dfloe       (i,j,  iblk),                           &
                            strax       (i,j,  iblk), stray       (i,j,  iblk), &
                            strairxT    (i,j,  iblk), strairyT    (i,j,  iblk), &
                            potT        (i,j,  iblk), sst         (i,j,  iblk), &
                            sss         (i,j,  iblk), Tf          (i,j,  iblk), &
                            strocnxT    (i,j,  iblk), strocnyT    (i,j,  iblk), &
                            fbot        (i,j,  iblk),                           &
                            Tbot        (i,j,  iblk), Tsnice       (i,j, iblk),  &
                            frzmlt      (i,j,  iblk), rside       (i,j,  iblk), &
                            fsnow       (i,j,  iblk), frain       (i,j,  iblk), &
                            fpond       (i,j,  iblk),                           &
                            fsurf       (i,j,  iblk), fsurfn      (i,j,:,iblk), &
                            fcondtop    (i,j,  iblk), fcondtopn   (i,j,:,iblk), &
                            fcondbot    (i,j,  iblk), fcondbotn   (i,j,:,iblk), &
                            fswsfcn     (i,j,:,iblk), fswintn     (i,j,:,iblk), &
                            fswthrun    (i,j,:,iblk), fswabs      (i,j,  iblk), &
                            flwout      (i,j,  iblk),                           &
                            Sswabsn   (i,j,:,:,iblk), Iswabsn   (i,j,:,:,iblk), &
                            flw         (i,j,  iblk), &
                            fsens       (i,j,  iblk), fsensn      (i,j,:,iblk), &
                            flat        (i,j,  iblk), flatn       (i,j,:,iblk), &
                            evap        (i,j,  iblk),                           &
                            evaps       (i,j,  iblk), evapi       (i,j,  iblk), &
                            fresh       (i,j,  iblk), fsalt       (i,j,  iblk), &
                            fhocn       (i,j,  iblk), fswthru     (i,j,  iblk), &
                            flatn_f     (i,j,:,iblk), fsensn_f    (i,j,:,iblk), &
                            fsurfn_f    (i,j,:,iblk), fcondtopn_f (i,j,:,iblk), &
                            faero_atm   (i,j,1:n_aero,iblk),                    &
                            faero_ocn   (i,j,1:n_aero,iblk),                    &
                            dhsn        (i,j,:,iblk), ffracn      (i,j,:,iblk), &
                            meltt       (i,j,  iblk), melttn      (i,j,:,iblk), &
                            meltb       (i,j,  iblk), meltbn      (i,j,:,iblk), &
                            melts       (i,j,  iblk), meltsn      (i,j,:,iblk), &
                            congel      (i,j,  iblk), congeln     (i,j,:,iblk), &
                            snoice      (i,j,  iblk), snoicen     (i,j,:,iblk), &
                            dsnown      (i,j,:,iblk), &
                            lmask_n     (i,j,  iblk), lmask_s     (i,j,  iblk), &
                            mlt_onset   (i,j,  iblk), frz_onset   (i,j,  iblk), &
                            yday,                     prescribed_ice)

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

      use ice_arrays_column, only: hin_max, fzsal, ocean_bio, &
          first_ice, bgrid, cgrid, igrid
      use ice_blocks, only: block, get_block
      use ice_calendar, only: yday
      use ice_domain, only: blocks_ice
      use ice_domain_size, only: ncat, nilyr, nslyr, n_aero, nblyr
      use ice_flux, only: fresh, frain, fpond, frzmlt, frazil, frz_onset, &
          update_ocn_f, fsalt, Tf, sss, salinz, fhocn, rside, &
          meltl, frazil_diag
      use ice_flux_bgc, only: flux_bio, faero_ocn 
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
         z_tracers

      type (block) :: &
         this_block      ! block information for current block

      character(len=*), parameter :: subname = '(step_therm2)'

      call icepack_query_parameters(z_tracers_out=z_tracers)
      call icepack_query_tracer_numbers(ntrcr_out=ntrcr, nbtrcr_out=nbtrcr)
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

         call icepack_step_therm2(dt, ncat, n_aero, nltrcr,                &
                           nilyr,                  nslyr,                  &
                           hin_max   (:),          nblyr,                  &   
                           aicen     (i,j,:,iblk),                         &
                           vicen     (i,j,:,iblk), vsnon     (i,j,:,iblk), &
                           aicen_init(i,j,:,iblk), vicen_init(i,j,:,iblk), &
                           trcrn     (i,j,:,:,iblk),                       &
                           aice0     (i,j,  iblk), aice      (i,j,  iblk), &
                           trcr_depend(:),         trcr_base(:,:),         &
                           n_trcr_strata(:),       nt_strata(:,:),         &
                           Tf        (i,j,  iblk), sss       (i,j,  iblk), &
                           salinz    (i,j,:,iblk),                         &
                           rside     (i,j,  iblk), meltl     (i,j,  iblk), &
                           frzmlt    (i,j,  iblk), frazil    (i,j,  iblk), &
                           frain     (i,j,  iblk), fpond     (i,j,  iblk), &
                           fresh     (i,j,  iblk), fsalt     (i,j,  iblk), &
                           fhocn     (i,j,  iblk), update_ocn_f,           &
                           bgrid,                  cgrid,                  &
                           igrid,                  faero_ocn (i,j,:,iblk), &
                           first_ice (i,j,:,iblk), fzsal     (i,j,  iblk), &
                           flux_bio  (i,j,1:nbtrcr,iblk),                  &
                           ocean_bio (i,j,1:nbtrcr,iblk),                  &
                           frazil_diag(i,j, iblk),                         &
                           frz_onset (i,j,  iblk), yday)

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
      call icepack_query_tracer_numbers(ntrcr_out=ntrcr)
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
         call icepack_aggregate (ncat,            aicen(i,j,:,iblk),    &
                               trcrn(i,j,:,:,iblk),                     &
                               vicen(i,j,:,iblk), vsnon(i,j,  :,iblk),  &
                               aice (i,j,  iblk),                       &
                               trcr (i,j,:,  iblk),                     &
                               vice (i,j,  iblk), vsno (i,j,    iblk),  &
                               aice0(i,j,  iblk),                       &
                               ntrcr,                                   &
                               trcr_depend(:),                          &
                               trcr_base    (:,:),                      &
                               n_trcr_strata(:),                        &
                               nt_strata    (:,:))

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
      use ice_flux_bgc, only: flux_bio, faero_ocn
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

      call icepack_query_tracer_numbers(ntrcr_out=ntrcr, nbtrcr_out=nbtrcr)
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

            call icepack_step_ridge (dt,        ndtd,                  &
                         nilyr,                 nslyr,                 &
                         nblyr,                                        &
                         ncat,                  hin_max  (:),          &
                         rdg_conv (i,j,  iblk), rdg_shear(i,j,  iblk), &
                         aicen    (i,j,:,iblk),                        &
                         trcrn    (i,j,:,:,iblk),                      &
                         vicen    (i,j,:,iblk), vsnon    (i,j,:,iblk), &
                         aice0    (i,j,  iblk), trcr_depend(:),        &
                         trcr_base(:,:),        n_trcr_strata(:),      &
                         nt_strata(:,:),                               &
                         dardg1dt (i,j,  iblk), dardg2dt (i,j,  iblk), &
                         dvirdgdt (i,j,  iblk), opening  (i,j,  iblk), &
                         fpond    (i,j,  iblk),                        &
                         fresh    (i,j,  iblk), fhocn    (i,j,  iblk), &
                         n_aero,                                       &
                         faero_ocn(i,j,:,iblk),                        &
                         aparticn (i,j,:,iblk), krdgn    (i,j,:,iblk), &
                         aredistn (i,j,:,iblk), vredistn (i,j,:,iblk), &
                         dardg1ndt(i,j,:,iblk), dardg2ndt(i,j,:,iblk), &
                         dvirdgndt(i,j,:,iblk),                        &
                         araftn   (i,j,:,iblk), vraftn   (i,j,:,iblk), &
                         aice     (i,j,  iblk), fsalt    (i,j,  iblk), &
                         first_ice(i,j,:,iblk), fzsal    (i,j,  iblk), &
                         flux_bio (i,j,1:nbtrcr,iblk)                  )

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
          fswsfcn, fswintn, fswthrun, fswpenln, Sswabsn, Iswabsn, &
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

      integer (kind=int_kind), dimension(icepack_max_aero) :: &
         nlt_zaero_sw, nt_zaero

      logical (kind=log_kind) :: &
         tr_bgc_N, tr_zaero, tr_brine, dEdd_algae, modal_aero

      real (kind=dbl_kind), dimension(ncat) :: &
         fbri                 ! brine height to ice thickness

      real(kind= dbl_kind), dimension(:,:), allocatable :: &
         ztrcr    , &
         ztrcr_sw

      logical (kind=log_kind) :: &
         debug, &           ! flag for printing debugging information
         l_print_point      ! flag for printing debugging information

      character(len=*), parameter :: subname = '(step_radiation)'

      call ice_timer_start(timer_sw)      ! shortwave

      call icepack_query_tracer_numbers(ntrcr_out=ntrcr, &
         nbtrcr_out=nbtrcr, nbtrcr_sw_out=nbtrcr_sw)
      call icepack_query_tracer_flags( &
         tr_brine_out=tr_brine, tr_bgc_N_out=tr_bgc_N, tr_zaero_out=tr_zaero)
      call icepack_query_tracer_indices( &
         nt_Tsfc_out=nt_Tsfc, nt_alvl_out=nt_alvl, &
         nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, nt_ipnd_out=nt_ipnd, nt_aero_out=nt_aero, &
         nlt_chl_sw_out=nlt_chl_sw, nlt_zaero_sw_out=nlt_zaero_sw, &
         nt_fbri_out=nt_fbri, nt_zaero_out=nt_zaero)
      call icepack_query_parameters(dEdd_algae_out=dEdd_algae, modal_aero_out=modal_aero)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      allocate(ztrcr(ntrcr,ncat))
      allocate(ztrcr_sw(ntrcr,ncat))

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
           do k = 1, ntrcr
             ztrcr(k,n) = trcrn(i,j,k,n,iblk)
           enddo
           if (tr_brine)  fbri(n) = trcrn(i,j,nt_fbri,n,iblk)
         enddo

         if (tmask(i,j,iblk)) then

            call icepack_step_radiation (dt,   ncat,                      &
                          n_algae,   tr_zaero, nblyr,                     &
                          ntrcr,     nbtrcr_sw,                           &
                          nilyr,    nslyr,       n_aero,                  &
                          n_zaero,  dEdd_algae,  nlt_chl_sw,              &
                          nlt_zaero_sw(:),                                &
                          swgrid(:),           igrid(:),                  &
                          fbri(:),                                        &
                          aicen(i,j,:,iblk),     vicen(i,j,:,iblk),       &
                          vsnon(i,j,:,iblk),                              &
                          trcrn(i,j,nt_Tsfc,:,iblk),                      &
                          trcrn(i,j,nt_alvl,:,iblk),                      &
                          trcrn(i,j,nt_apnd,:,iblk),                      &
                          trcrn(i,j,nt_hpnd,:,iblk),                      &
                          trcrn(i,j,nt_ipnd,:,iblk),                      &
                          trcrn(i,j,nt_aero:nt_aero+4*n_aero-1,:,iblk),   &
                          ztrcr_sw,                                       &
                          ztrcr,                                          &
                          TLAT(i,j,iblk),        TLON(i,j,iblk),          &
                          calendar_type,         days_per_year,           &
                          nextsw_cday,           yday,                    &
                          sec,                                            &
                          kaer_tab, waer_tab,                             &
                          gaer_tab,                                       &
                          kaer_bc_tab(:,:),      waer_bc_tab(:,:),        &
                          gaer_bc_tab(:,:),      bcenh(:,:,:),            &
                          modal_aero,                                     &
                          swvdr(i,j,iblk),       swvdf(i,j,iblk),         &
                          swidr(i,j,iblk),       swidf(i,j,iblk),         &
                          coszen(i,j,iblk),      fsnow(i,j,iblk),         &
                          alvdrn(i,j,:,iblk),    alvdfn(i,j,:,iblk),      &
                          alidrn(i,j,:,iblk),    alidfn(i,j,:,iblk),      &
                          fswsfcn(i,j,:,iblk),   fswintn(i,j,:,iblk),     &
                          fswthrun(i,j,:,iblk),  fswpenln(i,j,:,:,iblk),  &
                          Sswabsn(i,j,:,:,iblk), Iswabsn(i,j,:,:,iblk),   &
                          albicen(i,j,:,iblk),   albsnon(i,j,:,iblk),     &
                          albpndn(i,j,:,iblk),   apeffn(i,j,:,iblk),      &
                          snowfracn(i,j,:,iblk),                          &
                          dhsn(i,j,:,iblk),      ffracn(i,j,:,iblk),      &
                          l_print_point)

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

      deallocate(ztrcr)
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
      use ice_blocks, only: block, get_block
      use ice_domain, only: blocks_ice
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
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
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

      type (block) :: &
         this_block         ! block information for current block

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

            call icepack_atm_boundary( 'ocn',               &
                                     sst        (i,j,iblk), &    
                                     potT       (i,j,iblk), &
                                     uatm       (i,j,iblk), &   
                                     vatm       (i,j,iblk), &   
                                     wind       (i,j,iblk), &   
                                     zlvl       (i,j,iblk), &   
                                     Qa         (i,j,iblk), &     
                                     rhoa       (i,j,iblk), &
                                     strairx_ocn(i,j,iblk), & 
                                     strairy_ocn(i,j,iblk), & 
                                     Tref_ocn   (i,j,iblk), & 
                                     Qref_ocn   (i,j,iblk), & 
                                     delt       (i,j),      &    
                                     delq       (i,j),      &
                                     lhcoef     (i,j),      &
                                     shcoef     (i,j),      &
                                     Cdn_atm    (i,j,iblk), & 
                                     Cdn_atm_ratio(i,j,iblk))    
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

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         call icepack_ocn_mixed_layer (alvdr_ocn(i,j,iblk), swvdr     (i,j,iblk), &
                                      alidr_ocn(i,j,iblk), swidr     (i,j,iblk), &
                                      alvdf_ocn(i,j,iblk), swvdf     (i,j,iblk), &
                                      alidf_ocn(i,j,iblk), swidf     (i,j,iblk), &
                                      sst      (i,j,iblk), flwout_ocn(i,j,iblk), &
                                      fsens_ocn(i,j,iblk), shcoef    (i,j),      &
                                      flat_ocn (i,j,iblk), lhcoef    (i,j),      &
                                      evap_ocn (i,j,iblk), flw       (i,j,iblk), &
                                      delt     (i,j),      delq      (i,j),      &
                                      aice     (i,j,iblk), fhocn     (i,j,iblk), &
                                      fswthru  (i,j,iblk), hmix      (i,j,iblk), &
                                      Tf       (i,j,iblk), qdp       (i,j,iblk), &
                                      frzmlt   (i,j,iblk), dt)
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
      call icepack_query_tracer_numbers(ntrcr_out=ntrcr, nbtrcr_out=nbtrcr)
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

         call icepack_init_OceanConcArray(icepack_max_nbtrcr, &
                icepack_max_algae, icepack_max_don,  icepack_max_doc,        &
                icepack_max_dic,   icepack_max_aero, icepack_max_fe,         &
                nit(i,j,  iblk), amm   (i,j,  iblk), &
                sil(i,j,  iblk), dmsp  (i,j,  iblk), &
                dms(i,j,  iblk), algalN(i,j,:,iblk), &
                doc(i,j,:,iblk), don   (i,j,:,iblk), &
                dic(i,j,:,iblk), fed   (i,j,:,iblk), &
                fep(i,j,:,iblk), zaeros(i,j,:,iblk), &
                ocean_bio_all(i,j,:,iblk), &
                hum(i,j,  iblk))

         do mm = 1,nbtrcr
            ocean_bio(i,j,mm,iblk) = ocean_bio_all(i,j,bio_index_o(mm),iblk)  
         enddo  ! mm    
         if (tr_zaero) then
            do mm = 1, n_zaero  ! update aerosols
               flux_bio_atm(i,j,nlt_zaero(mm),iblk) = faero_atm(i,j,mm,iblk)
            enddo  ! mm
         endif

         call icepack_biogeochemistry(dt, ntrcr, nbtrcr,&
                              upNO        (i,j,          iblk),        &
                              upNH        (i,j,          iblk),        &
                              iDi         (i,j,:,:,      iblk),        &
                              iki         (i,j,:,:,      iblk),        &
                              zfswin      (i,j,:,:,      iblk),        &
                              zsal_tot    (i,j,          iblk),        &
                              darcy_V     (i,j,:,        iblk),        &
                              grow_net    (i,j,          iblk),        &
                              PP_net      (i,j,          iblk),        &
                              hbri        (i,j,          iblk),        &
                              dhbr_bot    (i,j,:,        iblk),        &
                              dhbr_top    (i,j,:,        iblk),        &
                              Zoo         (i,j,:,:,      iblk),        &
                              fbio_snoice (i,j,:,        iblk),        &
                              fbio_atmice (i,j,:,        iblk),        &
                              ocean_bio   (i,j,1:nbtrcr, iblk),        &
                              first_ice   (i,j,:,        iblk),        &
                              fswpenln    (i,j,:,:,      iblk),        &
                              bphi        (i,j,:,:,      iblk),        &
                              bTiz        (i,j,:,:,      iblk),        &
                              ice_bio_net (i,j,1:nbtrcr, iblk),        &
                              snow_bio_net(i,j,1:nbtrcr, iblk),        &
                              fswthrun    (i,j,:,        iblk),        &
                              Rayleigh_criteria(i,j,     iblk),        &
                              sice_rho    (i,j,:,        iblk),        &
                              fzsal       (i,j,          iblk),        &   
                              fzsal_g     (i,j,          iblk),        &
                              bgrid, igrid, icgrid, cgrid,             &
                              nblyr, nilyr, nslyr, n_algae, n_zaero,   &
                              ncat, n_doc, n_dic, n_don, n_fed, n_fep, &
                              meltbn      (i,j,:,        iblk),        &
                              melttn      (i,j,:,        iblk),        &
                              congeln     (i,j,:,        iblk),        &
                              snoicen     (i,j,:,        iblk),        & 
                              sst         (i,j,          iblk),        &    
                              sss         (i,j,          iblk),        &
                              fsnow       (i,j,          iblk),        &
                              meltsn      (i,j,:,        iblk),        &
                              hin_old     (i,j,:,        iblk),        &
                              flux_bio    (i,j,1:nbtrcr, iblk),        &
                              flux_bio_atm(i,j,1:nbtrcr, iblk),        &
                              aicen_init  (i,j,:,        iblk),        &
                              vicen_init  (i,j,:,        iblk),        &
                              aicen       (i,j,:,        iblk),        &
                              vicen       (i,j,:,        iblk),        &
                              vsnon       (i,j,:,        iblk),        &
                              aice0       (i,j,          iblk),        &
                              trcrn       (i,j,:,:,iblk),              &
                              vsnon_init  (i,j,:,        iblk),        &
                              skl_bgc)

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
