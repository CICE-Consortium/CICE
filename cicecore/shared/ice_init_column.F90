!=========================================================================
!
! Initialization routines for the column package.
!
! author: Elizabeth C. Hunke, LANL
!
      module ice_init_column

      use ice_kinds_mod
      use ice_constants
      use ice_communicate, only: my_task, master_task, ice_barrier
      use ice_domain_size, only: ncat, max_blocks 
      use ice_domain_size, only: nblyr, nilyr, nslyr
      use ice_domain_size, only: n_aero, n_zaero, n_algae
      use ice_domain_size, only: n_doc, n_dic, n_don
      use ice_domain_size, only: n_fed, n_fep
      use ice_fileunits, only: nu_diag
      use ice_fileunits, only: nu_nml, nml_filename, get_fileunit, &
                               release_fileunit, flush_fileunit
      use ice_exit, only: abort_ice
      use icepack_intfc, only: icepack_max_don, icepack_max_doc, icepack_max_dic
      use icepack_intfc, only: icepack_max_algae, icepack_max_aero, icepack_max_fe
      use icepack_intfc, only: icepack_max_nbtrcr
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_init_tracer_sizes, icepack_init_tracer_flags
      use icepack_intfc, only: icepack_init_tracer_indices
      use icepack_intfc, only: icepack_init_parameters
      use icepack_intfc, only: icepack_query_tracer_sizes, icepack_query_tracer_flags
      use icepack_intfc, only: icepack_query_tracer_indices, icepack_query_tracer_sizes
      use icepack_intfc, only: icepack_query_parameters
      use icepack_intfc, only: icepack_write_tracer_sizes, icepack_write_tracer_flags
      use icepack_intfc, only: icepack_write_tracer_indices, icepack_write_tracer_sizes
      use icepack_intfc, only: icepack_init_fsd, icepack_cleanup_fsd
      use icepack_intfc, only: icepack_init_zbgc
      use icepack_intfc, only: icepack_init_thermo
      use icepack_intfc, only: icepack_step_radiation, icepack_init_orbit
      use icepack_intfc, only: icepack_init_bgc, icepack_init_zsalinity
      use icepack_intfc, only: icepack_init_ocean_bio, icepack_load_ocean_bio_array
      use icepack_intfc, only: icepack_init_hbrine

      implicit none

      private
      public :: init_thermo_vertical, init_shortwave, &
                init_age, init_FY, init_lvl, init_fsd, &
                init_meltponds_cesm, init_meltponds_lvl, init_meltponds_topo, &
                init_aerosol, init_bgc, init_hbrine, init_zbgc, input_zbgc, &
                count_tracers, init_isotope

      ! namelist parameters needed locally

      real (kind=dbl_kind) :: &
         tau_min            , tau_max            , &
         nitratetype        , ammoniumtype       , silicatetype,  &
         dmspptype          , dmspdtype          , humtype

      real (kind=dbl_kind), dimension(icepack_max_dic) :: &
         dictype

      real (kind=dbl_kind), dimension(icepack_max_algae) :: &
         algaltype   ! tau_min for both retention and release

      real (kind=dbl_kind), dimension(icepack_max_doc) :: &
         doctype

      real (kind=dbl_kind), dimension(icepack_max_don) :: &
         dontype

      real (kind=dbl_kind), dimension(icepack_max_fe) :: &
         fedtype

      real (kind=dbl_kind), dimension(icepack_max_fe) :: &
         feptype

      real (kind=dbl_kind), dimension(icepack_max_aero) :: &
         zaerotype

      real (kind=dbl_kind) :: &
          grid_o, l_sk, grid_o_t, initbio_frac, &
          frazil_scav, grid_oS, l_skS, &
          phi_snow, &
          ratio_Si2N_diatoms , ratio_Si2N_sp      , ratio_Si2N_phaeo   ,  &
          ratio_S2N_diatoms  , ratio_S2N_sp       , ratio_S2N_phaeo    ,  &
          ratio_Fe2C_diatoms , ratio_Fe2C_sp      , ratio_Fe2C_phaeo   ,  &
          ratio_Fe2N_diatoms , ratio_Fe2N_sp      , ratio_Fe2N_phaeo   ,  &
          ratio_Fe2DON       , ratio_Fe2DOC_s     , ratio_Fe2DOC_l     ,  &
          fr_resp            , &
          algal_vel          , R_dFe2dust         , dustFe_sol         ,  &
          chlabs_diatoms     , chlabs_sp          , chlabs_phaeo       ,  &
          alpha2max_low_diatoms,alpha2max_low_sp  , alpha2max_low_phaeo,  &
          beta2max_diatoms   , beta2max_sp        , beta2max_phaeo     ,  &
          mu_max_diatoms     , mu_max_sp          , mu_max_phaeo       ,  &
          grow_Tdep_diatoms  , grow_Tdep_sp       , grow_Tdep_phaeo    ,  &
          fr_graze_diatoms   , fr_graze_sp        , fr_graze_phaeo     ,  &
          mort_pre_diatoms   , mort_pre_sp        , mort_pre_phaeo     ,  &
          mort_Tdep_diatoms  , mort_Tdep_sp       , mort_Tdep_phaeo    ,  &
          k_exude_diatoms    , k_exude_sp         , k_exude_phaeo      ,  &
          K_Nit_diatoms      , K_Nit_sp           , K_Nit_phaeo        ,  &
          K_Am_diatoms       , K_Am_sp            , K_Am_phaeo         ,  &
          K_Sil_diatoms      , K_Sil_sp           , K_Sil_phaeo        ,  &
          K_Fe_diatoms       , K_Fe_sp            , K_Fe_phaeo         ,  &
          f_don_protein      , kn_bac_protein     , f_don_Am_protein   ,  &
          f_doc_s            , f_doc_l            , f_exude_s          ,  &
          f_exude_l          , k_bac_s            , k_bac_l            ,  &
          T_max              , fsal               , op_dep_min         ,  &
          fr_graze_s         , fr_graze_e         , fr_mort2min        ,  &
          fr_dFe             , k_nitrif           , t_iron_conv        ,  &
          max_loss           , max_dfe_doc1       , fr_resp_s          ,  &
          y_sk_DMS           , t_sk_conv          , t_sk_ox            ,  &
          algaltype_diatoms  , algaltype_sp       , algaltype_phaeo    ,  &
          doctype_s          , doctype_l          , dontype_protein    ,  &
          fedtype_1          , feptype_1          , zaerotype_bc1      ,  &
          zaerotype_bc2      , zaerotype_dust1    , zaerotype_dust2    ,  &
          zaerotype_dust3    , zaerotype_dust4    , ratio_C2N_diatoms  ,  &
          ratio_C2N_sp       , ratio_C2N_phaeo    , ratio_chl2N_diatoms,  &
          ratio_chl2N_sp     , ratio_chl2N_phaeo  , F_abs_chl_diatoms  ,  &
          F_abs_chl_sp       , F_abs_chl_phaeo    , ratio_C2N_proteins

!=======================================================================

      contains

!=======================================================================
!
! Initialize the vertical profile of ice salinity and melting temperature.
!
! authors: C. M. Bitz, UW
!          William H. Lipscomb, LANL

      subroutine init_thermo_vertical

      use ice_blocks, only: nx_block, ny_block
      use ice_flux, only: salinz, Tmltz

      integer (kind=int_kind) :: &
         i, j, iblk, &  ! horizontal indices
         k              ! ice layer index

      real (kind=dbl_kind), dimension(nilyr+1) :: &
         sprofile                         ! vertical salinity profile

      real (kind=dbl_kind) :: &
         depressT

      character(len=*), parameter :: subname='(init_thermo_vertical)'

      !-----------------------------------------------------------------
      ! initialize heat_capacity, l_brine, and salinity profile
      !-----------------------------------------------------------------

      call icepack_query_parameters(depressT_out=depressT)
      call icepack_init_thermo(nilyr=nilyr, sprofile=sprofile)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__, line=__LINE__)

      !-----------------------------------------------------------------
      ! Prescibe vertical profile of salinity and melting temperature.
      ! Note this profile is only used for BL99 thermodynamics.
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,k)
      do iblk = 1,max_blocks
      do j = 1, ny_block
      do i = 1, nx_block
         do k = 1, nilyr+1
            salinz(i,j,k,iblk) = sprofile(k)
            Tmltz (i,j,k,iblk) = -salinz(i,j,k,iblk)*depressT
         enddo ! k
      enddo    ! i
      enddo    ! j
      enddo    ! iblk
      !$OMP END PARALLEL DO

      end subroutine init_thermo_vertical

!=======================================================================
!
!  Initialize shortwave

      subroutine init_shortwave

      use ice_arrays_column, only: fswpenln, Iswabsn, Sswabsn, albicen, &
          albsnon, alvdrn, alidrn, alvdfn, alidfn, fswsfcn, &
          fswthrun, fswthrun_vdr, fswthrun_vdf, fswthrun_idr, fswthrun_idf, &
          fswintn, albpndn, apeffn, trcrn_sw, dhsn, ffracn, snowfracn, &
          kaer_tab, waer_tab, gaer_tab, kaer_bc_tab, waer_bc_tab, gaer_bc_tab, bcenh, &
          swgrid, igrid
      use ice_blocks, only: block, get_block, nx_block, ny_block
      use ice_calendar, only: dt, calendar_type, &
          days_per_year, nextsw_cday, yday, sec
      use ice_diagnostics, only: npnt, print_points, pmloc, piloc, pjloc
      use ice_domain, only: nblocks, blocks_ice
      use ice_flux, only: alvdf, alidf, alvdr, alidr, &
                          alvdr_ai, alidr_ai, alvdf_ai, alidf_ai, &
                          swvdr, swvdf, swidr, swidf, scale_factor, snowfrac, &
                          albice, albsno, albpnd, apeff_ai, coszen, fsnow
      use ice_grid, only: tlat, tlon, tmask
      use ice_restart_shared, only: restart, runtype
      use ice_state, only: aicen, vicen, vsnon, trcrn

      integer (kind=int_kind) :: &
         i, j , k    , & ! horizontal indices
         iblk        , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         n                  ! thickness category index

      real (kind=dbl_kind) :: &
         netsw           ! flag for shortwave radiation presence

      type (block) :: &
         this_block      ! block information for current block

      logical (kind=log_kind) :: &
         l_print_point, & ! flag to print designated grid point diagnostics
         debug,         & ! if true, print diagnostics
         dEdd_algae,    & ! from icepack
         modal_aero       ! from icepack

      character (char_len) :: shortwave

      integer (kind=int_kind) :: &
         ipoint

      real (kind=dbl_kind), dimension(ncat) :: &
         fbri                 ! brine height to ice thickness

      real(kind=dbl_kind), allocatable :: &
         ztrcr_sw(:,:)    !

      logical (kind=log_kind) :: tr_brine, tr_zaero, tr_bgc_n
      integer (kind=int_kind) :: nt_alvl, nt_apnd, nt_hpnd, nt_ipnd, nt_aero, &
         nt_fbri, nt_tsfc, ntrcr, nbtrcr, nbtrcr_sw
      integer (kind=int_kind), dimension(icepack_max_algae) :: &
         nt_bgc_N
      integer (kind=int_kind), dimension(icepack_max_aero) :: &
         nt_zaero
      real (kind=dbl_kind) :: puny

      character(len=*), parameter :: subname='(init_shortwave)'

      call icepack_query_parameters(puny_out=puny)
      call icepack_query_parameters(shortwave_out=shortwave)
      call icepack_query_parameters(dEdd_algae_out=dEdd_algae)
      call icepack_query_parameters(modal_aero_out=modal_aero)
      call icepack_query_tracer_sizes(ntrcr_out=ntrcr, nbtrcr_out=nbtrcr, nbtrcr_sw_out=nbtrcr_sw)
      call icepack_query_tracer_flags(tr_brine_out=tr_brine, tr_zaero_out=tr_zaero, &
         tr_bgc_n_out=tr_bgc_n)
      call icepack_query_tracer_indices(nt_alvl_out=nt_alvl, nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, &
         nt_ipnd_out=nt_ipnd, nt_aero_out=nt_aero, nt_fbri_out=nt_fbri, nt_tsfc_out=nt_tsfc, &
         nt_bgc_N_out=nt_bgc_N, nt_zaero_out=nt_zaero)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__,line= __LINE__)

      allocate(ztrcr_sw(nbtrcr_sw, ncat))

      do iblk=1,nblocks

         ! Initialize
         fswpenln(:,:,:,:,iblk) = c0
         Iswabsn(:,:,:,:,iblk) = c0
         Sswabsn(:,:,:,:,iblk) = c0

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi
         
         do j = 1, ny_block ! can be jlo, jhi
         do i = 1, nx_block ! can be ilo, ihi

            l_print_point = .false.
            debug = .false.
            if (debug .and. print_points) then
               do ipoint = 1, npnt
                  if (my_task == pmloc(ipoint) .and. &
                       i == piloc(ipoint) .and. &
                       j == pjloc(ipoint)) &
                       l_print_point = .true.
                       write (nu_diag, *) 'my_task = ',my_task
               enddo ! n
            endif

            alvdf(i,j,iblk) = c0
            alidf(i,j,iblk) = c0
            alvdr(i,j,iblk) = c0
            alidr(i,j,iblk) = c0
            alvdr_ai(i,j,iblk) = c0
            alidr_ai(i,j,iblk) = c0
            alvdf_ai(i,j,iblk) = c0
            alidf_ai(i,j,iblk) = c0
            albice(i,j,iblk) = c0
            albsno(i,j,iblk) = c0
            albpnd(i,j,iblk) = c0
            snowfrac(i,j,iblk) = c0
            apeff_ai(i,j,iblk) = c0

            do n = 1, ncat
               alvdrn(i,j,n,iblk) = c0
               alidrn(i,j,n,iblk) = c0
               alvdfn(i,j,n,iblk) = c0
               alidfn(i,j,n,iblk) = c0
               fswsfcn(i,j,n,iblk) = c0
               fswintn(i,j,n,iblk) = c0
               fswthrun(i,j,n,iblk) = c0
               fswthrun_vdr(i,j,n,iblk) = c0
               fswthrun_vdf(i,j,n,iblk) = c0
               fswthrun_idr(i,j,n,iblk) = c0
               fswthrun_idf(i,j,n,iblk) = c0
            enddo   ! ncat

         enddo
         enddo
         do j = jlo, jhi
         do i = ilo, ihi

            if (trim(shortwave) == 'dEdd') then ! delta Eddington

#ifndef CESMCOUPLED
               ! initialize orbital parameters
               ! These come from the driver in the coupled model.
               call icepack_init_orbit()
               call icepack_warnings_flush(nu_diag)
               if (icepack_warnings_aborted()) call abort_ice(subname//' init_orbit', &
                  file=__FILE__, line=__LINE__)
#endif
            endif

            fbri(:) = c0
            ztrcr_sw(:,:) = c0
            do n = 1, ncat
              if (tr_brine)  fbri(n) = trcrn(i,j,nt_fbri,n,iblk)
            enddo

            if (tmask(i,j,iblk)) then
               call icepack_step_radiation (dt=dt, ncat=ncat,                  &
                          nblyr=nblyr, nilyr=nilyr, nslyr=nslyr,               &
                          dEdd_algae=dEdd_algae,                               &
                          swgrid=swgrid(:), igrid=igrid(:),                    &
                          fbri=fbri(:),                                        &
                          aicen=aicen(i,j,:,iblk),                             &
                          vicen=vicen(i,j,:,iblk),                             &
                          vsnon=vsnon(i,j,:,iblk),                             &
                          Tsfcn=trcrn(i,j,nt_Tsfc,:,iblk),                     &
                          alvln=trcrn(i,j,nt_alvl,:,iblk),                     &
                          apndn=trcrn(i,j,nt_apnd,:,iblk),                     &
                          hpndn=trcrn(i,j,nt_hpnd,:,iblk),                     &
                          ipndn=trcrn(i,j,nt_ipnd,:,iblk),                     &
                          aeron=trcrn(i,j,nt_aero:nt_aero+4*n_aero-1,:,iblk),  &
                          bgcNn=trcrn(i,j,nt_bgc_N(1):nt_bgc_N(1)+n_algae*(nblyr+3)-1,:,iblk), &
                          zaeron=trcrn(i,j,nt_zaero(1):nt_zaero(1)+n_zaero*(nblyr+3)-1,:,iblk), &
                          trcrn_bgcsw=ztrcr_sw,                                &
                          TLAT=TLAT(i,j,iblk), TLON=TLON(i,j,iblk),            &
                          calendar_type=calendar_type,                         &
                          days_per_year=days_per_year,                         &
                          nextsw_cday=nextsw_cday, yday=yday,                  &
                          sec=sec,                                             &
                          kaer_tab=kaer_tab, kaer_bc_tab=kaer_bc_tab(:,:),     &
                          waer_tab=waer_tab, waer_bc_tab=waer_bc_tab(:,:),     &
                          gaer_tab=gaer_tab, gaer_bc_tab=gaer_bc_tab(:,:),     &
                          bcenh=bcenh(:,:,:),                                  &
                          modal_aero=modal_aero,                               &
                          swvdr=swvdr(i,j,iblk),         swvdf=swvdf(i,j,iblk),&
                          swidr=swidr(i,j,iblk),         swidf=swidf(i,j,iblk),&
                          coszen=coszen(i,j,iblk),       fsnow=fsnow(i,j,iblk),&
                          alvdrn=alvdrn(i,j,:,iblk),     alvdfn=alvdfn(i,j,:,iblk), &
                          alidrn=alidrn(i,j,:,iblk),     alidfn=alidfn(i,j,:,iblk), &
                          fswsfcn=fswsfcn(i,j,:,iblk),   fswintn=fswintn(i,j,:,iblk), &
                          fswthrun=fswthrun(i,j,:,iblk),                       &
                          fswthrun_vdr=fswthrun_vdr(i,j,:,iblk),               &
                          fswthrun_vdf=fswthrun_vdf(i,j,:,iblk),               &
                          fswthrun_idr=fswthrun_idr(i,j,:,iblk),               &
                          fswthrun_idf=fswthrun_idf(i,j,:,iblk),               &
                          fswpenln=fswpenln(i,j,:,:,iblk),                     &
                          Sswabsn=Sswabsn(i,j,:,:,iblk), Iswabsn=Iswabsn(i,j,:,:,iblk), &
                          albicen=albicen(i,j,:,iblk),   albsnon=albsnon(i,j,:,iblk), &
                          albpndn=albpndn(i,j,:,iblk),   apeffn=apeffn(i,j,:,iblk), &
                          snowfracn=snowfracn(i,j,:,iblk),                     &
                          dhsn=dhsn(i,j,:,iblk),         ffracn=ffracn(i,j,:,iblk), &
                          l_print_point=l_print_point,                         &
                          initonly = .true.)
            endif
         
      !-----------------------------------------------------------------
      ! Define aerosol tracer on shortwave grid
      !-----------------------------------------------------------------

            if (dEdd_algae .and. (tr_zaero .or. tr_bgc_N)) then
              do n = 1, ncat
                 do k = 1, nbtrcr_sw
                    trcrn_sw(i,j,k,n,iblk) = ztrcr_sw(k,n)
                 enddo
              enddo
            endif

         enddo ! i
         enddo ! j

      !-----------------------------------------------------------------
      ! Aggregate albedos 
      ! Match loop order in coupling_prep for same order of operations
      !-----------------------------------------------------------------

         do n = 1, ncat
         do j = jlo, jhi
         do i = ilo, ihi

               if (aicen(i,j,n,iblk) > puny) then
                  
                  alvdf(i,j,iblk) = alvdf(i,j,iblk) &
                       + alvdfn(i,j,n,iblk)*aicen(i,j,n,iblk)
                  alidf(i,j,iblk) = alidf(i,j,iblk) &
                       + alidfn(i,j,n,iblk)*aicen(i,j,n,iblk)
                  alvdr(i,j,iblk) = alvdr(i,j,iblk) &
                       + alvdrn(i,j,n,iblk)*aicen(i,j,n,iblk)
                  alidr(i,j,iblk) = alidr(i,j,iblk) &
                       + alidrn(i,j,n,iblk)*aicen(i,j,n,iblk)
                  
                  netsw = swvdr(i,j,iblk) + swidr(i,j,iblk) &
                        + swvdf(i,j,iblk) + swidf(i,j,iblk)
                  if (netsw > puny) then ! sun above horizon
                     albice(i,j,iblk) = albice(i,j,iblk) &
                          + albicen(i,j,n,iblk)*aicen(i,j,n,iblk)
                     albsno(i,j,iblk) = albsno(i,j,iblk) &
                          + albsnon(i,j,n,iblk)*aicen(i,j,n,iblk)
                     albpnd(i,j,iblk) = albpnd(i,j,iblk) &
                          + albpndn(i,j,n,iblk)*aicen(i,j,n,iblk)
                  endif
                  
                  apeff_ai(i,j,iblk) = apeff_ai(i,j,iblk) &
                       + apeffn(i,j,n,iblk)*aicen(i,j,n,iblk)
                  snowfrac(i,j,iblk) = snowfrac(i,j,iblk) &
                       + snowfracn(i,j,n,iblk)*aicen(i,j,n,iblk)
               
               endif ! aicen > puny

         enddo ! i
         enddo ! j
         enddo ! ncat

         do j = 1, ny_block
         do i = 1, nx_block

      !----------------------------------------------------------------
      ! Store grid box mean albedos and fluxes before scaling by aice
      !----------------------------------------------------------------

            alvdf_ai  (i,j,iblk) = alvdf  (i,j,iblk)
            alidf_ai  (i,j,iblk) = alidf  (i,j,iblk)
            alvdr_ai  (i,j,iblk) = alvdr  (i,j,iblk)
            alidr_ai  (i,j,iblk) = alidr  (i,j,iblk)

            ! for history averaging
!echmod?            cszn = c0
!echmod            if (coszen(i,j,iblk) > puny) cszn = c1
!echmod            do n = 1, nstreams
!echmod               albcnt(i,j,iblk,n) = albcnt(i,j,iblk,n) + cszn
!echmod            enddo

      !----------------------------------------------------------------
      ! Save net shortwave for scaling factor in scale_factor
      !----------------------------------------------------------------
            if (runtype == 'initial' .and. .not. restart) then
               scale_factor(i,j,iblk) = &
                      swvdr(i,j,iblk)*(c1 - alvdr_ai(i,j,iblk)) &
                    + swvdf(i,j,iblk)*(c1 - alvdf_ai(i,j,iblk)) &
                    + swidr(i,j,iblk)*(c1 - alidr_ai(i,j,iblk)) &
                    + swidf(i,j,iblk)*(c1 - alidf_ai(i,j,iblk))
            endif

         enddo ! i
         enddo ! j
      enddo ! iblk

      deallocate(ztrcr_sw)

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      end subroutine init_shortwave

!=======================================================================

!  Initialize ice age tracer (call prior to reading restart data)

      subroutine init_age(iage)

      real(kind=dbl_kind), dimension(:,:,:), intent(out) :: iage
      character(len=*),parameter :: subname='(init_age)'

      iage(:,:,:) = c0

      end subroutine init_age

!=======================================================================

!  Initialize ice FY tracer (call prior to reading restart data)

      subroutine init_FY(firstyear)

      real(kind=dbl_kind), dimension(:,:,:), intent(out) :: firstyear
      character(len=*),parameter :: subname='(init_FY)'

      firstyear(:,:,:) = c0

      end subroutine init_FY

!=======================================================================

!  Initialize ice lvl tracers (call prior to reading restart data)

      subroutine init_lvl(iblk, alvl, vlvl) 

      use ice_constants, only: c0, c1
      use ice_arrays_column, only: ffracn, dhsn

      integer (kind=int_kind), intent(in)  :: iblk

      real (kind=dbl_kind), dimension(:,:,:), intent(out) :: &
         alvl , & ! level ice area fraction
         vlvl     ! level ice volume
      character(len=*),parameter :: subname='(init_lvl)'

      alvl(:,:,:) = c1 ! level ice area fraction
      vlvl(:,:,:) = c1 ! level ice volume
      ffracn(:,:,:,iblk) = c0
      dhsn(:,:,:,iblk) = c0

      end subroutine init_lvl

!=======================================================================

!  Initialize melt ponds.

      subroutine init_meltponds_cesm(apnd, hpnd)

      real(kind=dbl_kind), dimension(:,:,:), intent(out) :: &
         apnd , & ! melt pond area fraction
         hpnd     ! melt pond depth
      character(len=*),parameter :: subname='(init_meltponds_cesm)'

      apnd(:,:,:) = c0
      hpnd(:,:,:) = c0

      end subroutine init_meltponds_cesm

!=======================================================================

!  Initialize melt ponds.

      subroutine init_meltponds_lvl(apnd, hpnd, ipnd, dhsn)

      real(kind=dbl_kind), dimension(:,:,:), intent(out) :: &
         apnd , & ! melt pond area fraction
         hpnd , & ! melt pond depth
         ipnd , & ! melt pond refrozen lid thickness
         dhsn     ! depth difference for snow on sea ice and pond ice
      character(len=*),parameter :: subname='(init_meltponds_lvl)'

      apnd(:,:,:) = c0
      hpnd(:,:,:) = c0
      ipnd(:,:,:) = c0
      dhsn(:,:,:) = c0

      end subroutine init_meltponds_lvl

!=======================================================================

!  Initialize melt ponds.

      subroutine init_meltponds_topo(apnd, hpnd, ipnd)

      real(kind=dbl_kind), dimension(:,:,:), intent(out) :: &
         apnd , & ! melt pond area fraction
         hpnd , & ! melt pond depth
         ipnd     ! melt pond refrozen lid thickness
      character(len=*),parameter :: subname='(init_meltponds_topo)'

      apnd(:,:,:) = c0
      hpnd(:,:,:) = c0
      ipnd(:,:,:) = c0
        
      end subroutine init_meltponds_topo

!=======================================================================

!  Initialize floe size distribution tracer (call prior to reading restart data)

      subroutine init_fsd(floesize)

      use ice_arrays_column, only: floe_rad_c, floe_binwidth, &
         wavefreq, dwavefreq, wave_sig_ht, wave_spectrum, &
         d_afsd_newi, d_afsd_latg, d_afsd_latm, d_afsd_wave, d_afsd_weld
      use ice_blocks, only: nx_block, ny_block
      use ice_domain_size, only: ncat, max_blocks, nfsd
      use ice_init, only: ice_ic
      use ice_state, only: aicen

      real(kind=dbl_kind), dimension(:,:,:,:,:), intent(out) :: &
         floesize            ! floe size distribution tracer

      ! local variables

      real (kind=dbl_kind), dimension(nfsd) :: &
         afsd                ! floe size distribution "profile"

      real (kind=dbl_kind), dimension(nfsd,ncat) :: &
         afsdn               ! floe size distribution "profile"

      real (kind=dbl_kind) :: puny

      integer (kind=int_kind) :: &
         i, j, iblk     , &  ! horizontal indices
         n, k                ! category index

      logical (kind=log_kind) :: tr_fsd

      character(len=*), parameter :: subname='(init_fsd)'

      call icepack_query_parameters(puny_out=puny)

      wavefreq       (:)       = c0
      dwavefreq      (:)       = c0
      wave_sig_ht    (:,:,:)   = c0
      wave_spectrum  (:,:,:,:) = c0
      d_afsd_newi    (:,:,:,:) = c0
      d_afsd_latg    (:,:,:,:) = c0
      d_afsd_latm    (:,:,:,:) = c0
      d_afsd_wave    (:,:,:,:) = c0
      d_afsd_weld    (:,:,:,:) = c0

      ! default: floes occupy the smallest size category in all thickness categories
      afsdn(:,:) = c0
      afsdn(1,:) = c1
      floesize(:,:,:,:,:) = c0
      floesize(:,:,1,:,:) = c1

      call icepack_query_tracer_flags(tr_fsd_out=tr_fsd)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__,line= __LINE__)

      if (tr_fsd) then

         ! initialize floe size distribution the same in every column and category
         call icepack_init_fsd(nfsd, ice_ic, &
            floe_rad_c,    &  ! fsd size bin centre in m (radius)
            floe_binwidth, &  ! fsd size bin width in m (radius)
            afsd)             ! floe size distribution

         do iblk = 1, max_blocks
            do j = 1, ny_block
            do i = 1, nx_block
               do n = 1, ncat
               do k = 1, nfsd
                  if (aicen(i,j,n,iblk) > puny) afsdn(k,n) = afsd(k)
               enddo    ! k
               enddo    ! n

               call icepack_cleanup_fsd (ncat, nfsd, afsdn) ! renormalize

               do n = 1, ncat
               do k = 1, nfsd
                  floesize(i,j,k,n,iblk) = afsdn(k,n)
               enddo    ! k
               enddo    ! n
            enddo       ! i
            enddo       ! j
         enddo          ! iblk

         call icepack_warnings_flush(nu_diag)
         if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
            file=__FILE__, line=__LINE__)

      endif ! tr_fsd

      end subroutine init_fsd

!=======================================================================

!  Initialize isotope tracers (call prior to reading restart data)

      subroutine init_isotope(isosno, isoice)

      real(kind=dbl_kind), dimension(:,:,:,:), intent(out) :: &
         isosno, isoice
      character(len=*),parameter :: subname='(init_isotope)'

      isosno(:,:,:,:) = c0
      isoice(:,:,:,:) = c0

      end subroutine init_isotope

!=======================================================================

!  Initialize ice aerosol tracer (call prior to reading restart data)

      subroutine init_aerosol(aero)

      real(kind=dbl_kind), dimension(:,:,:,:), intent(out) :: &
         aero ! aerosol tracers
      character(len=*),parameter :: subname='(init_aerosol)'

      aero(:,:,:,:) = c0

      end subroutine init_aerosol

!=======================================================================

!  Initialize vertical profile for biogeochemistry

      subroutine init_bgc() 

      use ice_arrays_column, only: zfswin, trcrn_sw, &
          ocean_bio_all, ice_bio_net, snow_bio_net, &
          cgrid, igrid, bphi, iDi, bTiz, iki, &
          Rayleigh_criteria, Rayleigh_real
      use ice_blocks, only: block, get_block
      use ice_domain, only: nblocks, blocks_ice
      use ice_flux, only: sss
      use ice_flux_bgc, only: nit, amm, sil, dmsp, dms, algalN, &
          doc, don, dic, fed, fep, zaeros, hum
      use ice_forcing_bgc, only: init_bgc_data, get_forcing_bgc
      use ice_restart_column, only: restart_zsal, &
          read_restart_bgc, restart_bgc
      use ice_state, only: trcrn

      ! local variables

      integer (kind=int_kind) :: &
         i, j, iblk       , & ! horizontal indices
         ilo,ihi,jlo,jhi  , & ! beginning and end of physical domain
         k                , & ! vertical index
         n                    ! category index

      integer (kind=int_kind) :: &
         max_nbtrcr, max_algae, max_don, max_doc, max_dic, max_aero, max_fe

      logical (kind=log_kind) :: &
         RayleighC    , &
         solve_zsal
        
      type (block) :: &
         this_block      ! block information for current block

      real(kind=dbl_kind), allocatable :: &
         trcrn_bgc(:,:)
      
      real(kind=dbl_kind), dimension(nilyr,ncat) :: &
         sicen    

      real(kind=dbl_kind) :: &
         RayleighR

      integer (kind=int_kind) :: &
         nbtrcr, ntrcr, ntrcr_o, &
         nt_sice, nt_bgc_S

      character(len=*), parameter :: subname='(init_bgc)'

      ! Initialize

      call icepack_query_parameters(solve_zsal_out=solve_zsal)
      call icepack_query_tracer_sizes(nbtrcr_out=nbtrcr, ntrcr_out=ntrcr, ntrcr_o_out=ntrcr_o)
      call icepack_query_tracer_indices(nt_sice_out=nt_sice, nt_bgc_S_out=nt_bgc_S)
      call icepack_query_tracer_sizes(max_nbtrcr_out=max_nbtrcr,         &
           max_algae_out=max_algae, max_don_out=max_don, max_doc_out=max_doc,   &
           max_dic_out=max_dic, max_aero_out=max_aero, max_fe_out=max_fe)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__,line= __LINE__)

      allocate(trcrn_bgc(ntrcr,ncat))

      bphi(:,:,:,:,:) = c0   ! initial porosity for no ice 
      iDi (:,:,:,:,:) = c0   ! interface diffusivity
      bTiz(:,:,:,:,:) = c0   ! initial bio grid ice temperature
      iki (:,:,:,:,:) = c0   ! permeability

      ocean_bio_all(:,:,:,:)   = c0
      ice_bio_net  (:,:,:,:)   = c0 ! integrated ice tracer conc (mmol/m^2 or mg/m^2) 
      snow_bio_net (:,:,:,:)   = c0 ! integrated snow tracer conc (mmol/m^2 or mg/m^2)
      zfswin       (:,:,:,:,:) = c0 ! shortwave flux on bio grid
      trcrn_sw     (:,:,:,:,:) = c0 ! tracers active in the shortwave calculation
      trcrn_bgc    (:,:)       = c0
      RayleighR                = c0
      RayleighC                = .false.

      !-----------------------------------------------------------------
      ! zsalinity initialization
      !-----------------------------------------------------------------
      
      if (solve_zsal) then  ! default values 

         !$OMP PARALLEL DO PRIVATE(iblk,i,j,k,n,ilo,ihi,jlo,jhi,this_block,trcrn_bgc)
         do iblk = 1, nblocks

            this_block = get_block(blocks_ice(iblk),iblk)         
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            do j = jlo, jhi
            do i = ilo, ihi  
               call icepack_init_zsalinity(nblyr=nblyr, ntrcr_o=ntrcr_o, &
                            Rayleigh_criteria = RayleighC, &
                            Rayleigh_real     = RayleighR, &
                            trcrn_bgc         = trcrn_bgc, &
                            nt_bgc_S          = nt_bgc_S,  &
                            ncat              = ncat,      &
                            sss               = sss(i,j,iblk))
               if (.not. restart_zsal) then
                  Rayleigh_real    (i,j,iblk) = RayleighR
                  Rayleigh_criteria(i,j,iblk) = RayleighC
                  do n = 1,ncat
                     do k  = 1, nblyr
                        trcrn    (i,j,nt_bgc_S+k-1,        n,iblk) = &
                        trcrn_bgc(    nt_bgc_S+k-1-ntrcr_o,n)
                     enddo
                  enddo
               endif
            enddo      ! i
            enddo      ! j  
         enddo         ! iblk
         !$OMP END PARALLEL DO
         call icepack_warnings_flush(nu_diag)
         if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
            file=__FILE__, line=__LINE__)
      endif ! solve_zsal

      if (.not. solve_zsal) restart_zsal = .false.

      !-----------------------------------------------------------------
      ! biogeochemistry initialization
      !-----------------------------------------------------------------

      if (.not. restart_bgc) then       
     
      !-----------------------------------------------------------------
      ! Initial Ocean Values if not coupled to the ocean bgc
      !-----------------------------------------------------------------
         !$OMP PARALLEL DO PRIVATE(iblk,i,j,n,ilo,ihi,jlo,jhi,this_block)
         do iblk = 1, nblocks

            this_block = get_block(blocks_ice(iblk),iblk)         
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            do j = jlo, jhi
            do i = ilo, ihi
               call icepack_init_ocean_bio ( &
                    amm=amm      (i,j,  iblk), dmsp=dmsp(i,j,  iblk), dms=dms(i,j,  iblk), &
                    algalN=algalN(i,j,:,iblk), doc=doc  (i,j,:,iblk), dic=dic(i,j,:,iblk), &
                    don=don      (i,j,:,iblk), fed=fed  (i,j,:,iblk), fep=fep(i,j,:,iblk), &
                    hum=hum      (i,j,  iblk), nit=nit  (i,j,  iblk), sil=sil(i,j,  iblk), &
                    zaeros=zaeros(i,j,:,iblk), &
                    max_dic = icepack_max_dic, max_don  = icepack_max_don, &
                    max_fe  = icepack_max_fe,  max_aero = icepack_max_aero)
            enddo  ! i
            enddo  ! j

         enddo     ! iblk
         !$OMP END PARALLEL DO

         call icepack_warnings_flush(nu_diag)
         if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
            file=__FILE__, line=__LINE__)

         call init_bgc_data(fed(:,:,1,:),fep(:,:,1,:)) ! input dFe from file
         call get_forcing_bgc                          ! defines nit and sil

      endif     ! .not. restart

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,n,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi  

            call icepack_load_ocean_bio_array(max_nbtrcr=icepack_max_nbtrcr,     &
                         max_algae=icepack_max_algae, max_don=icepack_max_don,   &
                         max_doc=icepack_max_doc,     max_fe=icepack_max_fe,     &
                         max_dic=icepack_max_dic,     max_aero=icepack_max_aero, &
                         nit =nit (i,j,  iblk), amm=amm(i,j,  iblk), sil   =sil   (i,j,  iblk), &
                         dmsp=dmsp(i,j,  iblk), dms=dms(i,j,  iblk), algalN=algalN(i,j,:,iblk), &
                         doc =doc (i,j,:,iblk), don=don(i,j,:,iblk), dic   =dic   (i,j,:,iblk), &  
                         fed =fed (i,j,:,iblk), fep=fep(i,j,:,iblk), zaeros=zaeros(i,j,:,iblk), &
                         hum=hum  (i,j,  iblk), ocean_bio_all=ocean_bio_all(i,j,:,iblk))

         enddo  ! i
         enddo  ! j

      enddo     ! iblk
      !$OMP END PARALLEL DO

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (.not. restart_bgc) then       
         !$OMP PARALLEL DO PRIVATE(iblk,i,j,k,n,ilo,ihi,jlo,jhi,this_block,sicen,trcrn_bgc)
         do iblk = 1, nblocks

            this_block = get_block(blocks_ice(iblk),iblk)         
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            do j = jlo, jhi
            do i = ilo, ihi  
                do n = 1, ncat
                do k = 1, nilyr
                   sicen(k,n) = trcrn(i,j,nt_sice+k-1,n,iblk)
                enddo
                do k = ntrcr_o+1, ntrcr
                   trcrn_bgc(k-ntrcr_o,n) = trcrn(i,j,k,n,iblk)
                enddo
                enddo
            call icepack_init_bgc(ncat=ncat, nblyr=nblyr, nilyr=nilyr, ntrcr_o=ntrcr_o,  &
                         cgrid=cgrid, igrid=igrid, ntrcr=ntrcr, nbtrcr=nbtrcr,           &
                         sicen=sicen(:,:), trcrn=trcrn_bgc(:,:), sss=sss(i,j, iblk), &
                         ocean_bio_all=ocean_bio_all(i,j,:,iblk))
            enddo  ! i
            enddo  ! j
         enddo     ! iblk
         !$OMP END PARALLEL DO

         call icepack_warnings_flush(nu_diag)
         if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
            file=__FILE__, line=__LINE__)

      endif ! .not. restart

      !-----------------------------------------------------------------
      ! read restart to complete BGC initialization
      !-----------------------------------------------------------------

      if (restart_zsal .or. restart_bgc) call read_restart_bgc  

      deallocate(trcrn_bgc)

      end subroutine init_bgc

!=======================================================================

!  Initialize brine height tracer

      subroutine init_hbrine()

      use ice_arrays_column, only: first_ice, bgrid, igrid, cgrid, &
          icgrid, swgrid
      use ice_state, only: trcrn

      real (kind=dbl_kind) :: phi_snow
      integer (kind=int_kind) :: nt_fbri
      logical (kind=log_kind) :: tr_brine
      character(len=*), parameter :: subname='(init_hbrine)'

      call icepack_query_parameters(phi_snow_out=phi_snow)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__,line= __LINE__)

      call icepack_init_hbrine(bgrid=bgrid, igrid=igrid, cgrid=cgrid, icgrid=icgrid, &
            swgrid=swgrid, nblyr=nblyr, nilyr=nilyr, phi_snow=phi_snow)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__, line=__LINE__)

      call icepack_init_parameters(phi_snow_in=phi_snow)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__,line= __LINE__)

      call icepack_query_tracer_flags(tr_brine_out=tr_brine)
      call icepack_query_tracer_indices(nt_fbri_out=nt_fbri)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__,line= __LINE__)

      first_ice(:,:,:,:) = .true.            
      if (tr_brine) trcrn(:,:,nt_fbri,:,:) = c1

      end subroutine init_hbrine

!=======================================================================

! Namelist variables, set to default values; may be altered at run time
! 
! author Elizabeth C. Hunke, LANL
!        Nicole Jeffery, LANL

      subroutine input_zbgc

      use ice_arrays_column, only: restore_bgc
      use ice_broadcast, only: broadcast_scalar
      use ice_restart_column, only: restart_bgc, restart_zsal, &
          restart_hbrine
      use ice_restart_shared, only: restart

      character (len=char_len) :: &
         shortwave        ! from icepack

      logical (kind=log_kind) :: &
         tr_brine, &
         tr_bgc_Nit,    tr_bgc_Am,    tr_bgc_Sil,   &
         tr_bgc_DMS,    tr_bgc_PON,   &
         tr_bgc_N,      tr_bgc_C,     tr_bgc_chl,   &
         tr_bgc_DON,    tr_bgc_Fe,    tr_zaero,     &
         tr_bgc_hum,    tr_aero
 
      integer (kind=int_kind) :: &
         ktherm

      logical (kind=log_kind) :: &
         solve_zsal, skl_bgc, z_tracers, scale_bgc, solve_zbgc, dEdd_algae, &
         modal_aero

      character (char_len) :: &
         bgc_flux_type

      integer (kind=int_kind) :: &
         nml_error, & ! namelist i/o error flag
         abort_flag

      character(len=*), parameter :: subname='(input_zbgc)'

      !-----------------------------------------------------------------
      ! namelist variables
      !-----------------------------------------------------------------

      namelist /zbgc_nml/  &
        tr_brine, restart_hbrine, tr_zaero, modal_aero, skl_bgc, &
        z_tracers, dEdd_algae, solve_zbgc, bgc_flux_type, &
        restore_bgc, restart_bgc, scale_bgc, solve_zsal, restart_zsal, &
        tr_bgc_Nit, tr_bgc_C, tr_bgc_chl, tr_bgc_Am, tr_bgc_Sil, &
        tr_bgc_DMS, tr_bgc_PON, tr_bgc_hum, tr_bgc_DON, tr_bgc_Fe, &
        grid_o, grid_o_t, l_sk, grid_oS, &   
        l_skS, phi_snow,  initbio_frac, frazil_scav, &
        ratio_Si2N_diatoms , ratio_Si2N_sp      , ratio_Si2N_phaeo   ,  &
        ratio_S2N_diatoms  , ratio_S2N_sp       , ratio_S2N_phaeo    ,  &
        ratio_Fe2C_diatoms , ratio_Fe2C_sp      , ratio_Fe2C_phaeo   ,  &
        ratio_Fe2N_diatoms , ratio_Fe2N_sp      , ratio_Fe2N_phaeo   ,  &
        ratio_Fe2DON       , ratio_Fe2DOC_s     , ratio_Fe2DOC_l     ,  &
        fr_resp            , tau_min            , tau_max            ,  &
        algal_vel          , R_dFe2dust         , dustFe_sol         ,  &
        chlabs_diatoms     , chlabs_sp          , chlabs_phaeo       ,  &
        alpha2max_low_diatoms,alpha2max_low_sp  , alpha2max_low_phaeo,  &
        beta2max_diatoms   , beta2max_sp        , beta2max_phaeo     ,  &
        mu_max_diatoms     , mu_max_sp          , mu_max_phaeo       ,  &
        grow_Tdep_diatoms  , grow_Tdep_sp       , grow_Tdep_phaeo    ,  &
        fr_graze_diatoms   , fr_graze_sp        , fr_graze_phaeo     ,  &
        mort_pre_diatoms   , mort_pre_sp        , mort_pre_phaeo     ,  &
        mort_Tdep_diatoms  , mort_Tdep_sp       , mort_Tdep_phaeo    ,  &
        k_exude_diatoms    , k_exude_sp         , k_exude_phaeo      ,  &
        K_Nit_diatoms      , K_Nit_sp           , K_Nit_phaeo        ,  &
        K_Am_diatoms       , K_Am_sp            , K_Am_phaeo         ,  &
        K_Sil_diatoms      , K_Sil_sp           , K_Sil_phaeo        ,  &
        K_Fe_diatoms       , K_Fe_sp            , K_Fe_phaeo         ,  &
        f_don_protein      , kn_bac_protein     , f_don_Am_protein   ,  &
        f_doc_s            , f_doc_l            , f_exude_s          ,  &
        f_exude_l          , k_bac_s            , k_bac_l            ,  &
        T_max              , fsal               , op_dep_min         ,  &
        fr_graze_s         , fr_graze_e         , fr_mort2min        ,  &
        fr_dFe             , k_nitrif           , t_iron_conv        ,  &
        max_loss           , max_dfe_doc1       , fr_resp_s          ,  &
        y_sk_DMS           , t_sk_conv          , t_sk_ox            ,  &
        algaltype_diatoms  , algaltype_sp       , algaltype_phaeo    ,  &
        nitratetype        , ammoniumtype       , silicatetype       ,  &
        dmspptype          , dmspdtype          , humtype            ,  &
        doctype_s          , doctype_l          , dontype_protein    ,  &
        fedtype_1          , feptype_1          , zaerotype_bc1      ,  &
        zaerotype_bc2      , zaerotype_dust1    , zaerotype_dust2    ,  &
        zaerotype_dust3    , zaerotype_dust4    , ratio_C2N_diatoms  ,  &
        ratio_C2N_sp       , ratio_C2N_phaeo    , ratio_chl2N_diatoms,  & 
        ratio_chl2N_sp     , ratio_chl2N_phaeo  , F_abs_chl_diatoms  ,  &
        F_abs_chl_sp       , F_abs_chl_phaeo    , ratio_C2N_proteins 

      !-----------------------------------------------------------------

      abort_flag = 0

      call icepack_query_tracer_flags(tr_aero_out=tr_aero)
      call icepack_query_parameters(ktherm_out=ktherm, shortwave_out=shortwave)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__, line=__LINE__)

      !-----------------------------------------------------------------
      ! default values
      !-----------------------------------------------------------------
      tr_brine        = .false.  ! brine height differs from ice height
      tr_zaero        = .false.  ! z aerosol tracers
      modal_aero      = .false.  ! use modal aerosol treatment of aerosols
      restore_bgc     = .false.  ! restore bgc if true
      solve_zsal      = .false.  ! update salinity tracer profile from solve_S_dt
      restart_bgc     = .false.  ! biogeochemistry restart
      restart_zsal    = .false.  ! salinity restart
      restart_hbrine  = .false.  ! hbrine restart
      scale_bgc       = .false.  ! initial bgc tracers proportional to S  
      skl_bgc         = .false.  ! solve skeletal biochemistry 
      z_tracers       = .false.  ! solve vertically resolved tracers
      dEdd_algae      = .false.  ! dynamic algae contributes to shortwave absorption
                                 ! in delta-Eddington calculation
      solve_zbgc      = .false.  ! turn on z layer biochemistry 
      tr_bgc_PON      = .false.  !---------------------------------------------   
      tr_bgc_Nit      = .false.  ! biogeochemistry (skl or zbgc)
      tr_bgc_C        = .false.  ! if skl_bgc = .true. then skl
      tr_bgc_chl      = .false.  ! if z_tracers = .true. then vertically resolved
      tr_bgc_Sil      = .false.  ! if z_tracers + solve_zbgc = .true. then
      tr_bgc_Am       = .false.  ! vertically resolved with reactions  
      tr_bgc_DMS      = .false.  !------------------------------------------------
      tr_bgc_DON      = .false.  ! 
      tr_bgc_hum      = .false.  !
      tr_bgc_Fe       = .false.  ! 
      tr_bgc_N        = .true.   !

      ! brine height parameter
      phi_snow        = p5       ! snow porosity

      ! skl biology parameters
      bgc_flux_type   = 'Jin2006'! type of ocean-ice poston velocity ('constant')

      ! z biology parameters  
      grid_o          = c5           ! for bottom flux        
      grid_o_t        = c5           ! for top flux        
      l_sk            = 7.0_dbl_kind ! characteristic diffusive scale (m)   
      initbio_frac    = c1           ! fraction of ocean trcr concentration in bio trcrs
      frazil_scav     = c1           ! increase in initial bio tracer from ocean scavenging 
      ratio_Si2N_diatoms = 1.8_dbl_kind    ! algal Si to N (mol/mol)                       
      ratio_Si2N_sp      = c0              ! diatoms, small plankton, phaeocystis
      ratio_Si2N_phaeo   = c0
      ratio_S2N_diatoms  = 0.03_dbl_kind   ! algal S  to N (mol/mol)
      ratio_S2N_sp       = 0.03_dbl_kind 
      ratio_S2N_phaeo    = 0.03_dbl_kind
      ratio_Fe2C_diatoms = 0.0033_dbl_kind ! algal Fe to C  (umol/mol)
      ratio_Fe2C_sp      = 0.0033_dbl_kind
      ratio_Fe2C_phaeo   = p1
      ratio_Fe2N_diatoms = 0.023_dbl_kind  ! algal Fe to N  (umol/mol)
      ratio_Fe2N_sp      = 0.023_dbl_kind
      ratio_Fe2N_phaeo   = 0.7_dbl_kind
      ratio_Fe2DON       = 0.023_dbl_kind  ! Fe to N of DON (nmol/umol)
      ratio_Fe2DOC_s     = p1              ! Fe to C of DOC (nmol/umol) saccharids
      ratio_Fe2DOC_l     = 0.033_dbl_kind  ! Fe to C of DOC (nmol/umol) lipids
      fr_resp            = 0.05_dbl_kind   ! frac of algal growth lost due to respiration      
      tau_min            = 5200.0_dbl_kind ! rapid mobile to stationary exchanges (s)
      tau_max            = 1.73e5_dbl_kind ! long time mobile to stationary exchanges (s)
      algal_vel          = 1.11e-8_dbl_kind! 0.5 cm/d(m/s) Lavoie 2005  1.5 cm/day
      R_dFe2dust         = 0.035_dbl_kind  !  g/g (3.5% content) Tagliabue 2009
      dustFe_sol         = 0.005_dbl_kind  ! solubility fraction
      chlabs_diatoms     = 0.03_dbl_kind   ! chl absorption (1/m/(mg/m^3))
      chlabs_sp          = 0.01_dbl_kind
      chlabs_phaeo       = 0.05_dbl_kind
      alpha2max_low_diatoms = 0.8_dbl_kind ! light limitation (1/(W/m^2))  
      alpha2max_low_sp      = 0.67_dbl_kind
      alpha2max_low_phaeo   = 0.67_dbl_kind
      beta2max_diatoms   = 0.018_dbl_kind  ! light inhibition (1/(W/m^2))  
      beta2max_sp        = 0.0025_dbl_kind
      beta2max_phaeo     = 0.01_dbl_kind
      mu_max_diatoms     = 1.2_dbl_kind    ! maximum growth rate (1/day) 
      mu_max_sp          = 0.851_dbl_kind
      mu_max_phaeo       = 0.851_dbl_kind
      grow_Tdep_diatoms  = 0.06_dbl_kind ! Temperature dependence of growth (1/C)
      grow_Tdep_sp       = 0.06_dbl_kind
      grow_Tdep_phaeo    = 0.06_dbl_kind
      fr_graze_diatoms   = 0.01_dbl_kind ! Fraction grazed
      fr_graze_sp        = p1
      fr_graze_phaeo     = p1
      mort_pre_diatoms   = 0.007_dbl_kind! Mortality (1/day)
      mort_pre_sp        = 0.007_dbl_kind
      mort_pre_phaeo     = 0.007_dbl_kind
      mort_Tdep_diatoms  = 0.03_dbl_kind ! T dependence of mortality (1/C)
      mort_Tdep_sp       = 0.03_dbl_kind
      mort_Tdep_phaeo    = 0.03_dbl_kind
      k_exude_diatoms    = c0            ! algal exudation (1/d)
      k_exude_sp         = c0
      k_exude_phaeo      = c0
      K_Nit_diatoms      = c1            ! nitrate half saturation (mmol/m^3)
      K_Nit_sp           = c1
      K_Nit_phaeo        = c1
      K_Am_diatoms       = 0.3_dbl_kind  ! ammonium half saturation (mmol/m^3)
      K_Am_sp            = 0.3_dbl_kind
      K_Am_phaeo         = 0.3_dbl_kind
      K_Sil_diatoms      = 4.0_dbl_kind  ! silicate half saturation (mmol/m^3)
      K_Sil_sp           = c0
      K_Sil_phaeo        = c0
      K_Fe_diatoms       = c1            ! iron half saturation (nM)
      K_Fe_sp            = 0.2_dbl_kind
      K_Fe_phaeo         = p1
      f_don_protein      = 0.6_dbl_kind  ! fraction of spilled grazing to proteins           
      kn_bac_protein     = 0.03_dbl_kind ! Bacterial degredation of DON (1/d)                
      f_don_Am_protein   = 0.25_dbl_kind ! fraction of remineralized DON to ammonium         
      f_doc_s            = 0.4_dbl_kind  ! fraction of mortality to DOC 
      f_doc_l            = 0.4_dbl_kind
      f_exude_s          = c1            ! fraction of exudation to DOC
      f_exude_l          = c1
      k_bac_s            = 0.03_dbl_kind ! Bacterial degredation of DOC (1/d)
      k_bac_l            = 0.03_dbl_kind
      T_max              = c0            ! maximum temperature (C)
      fsal               = c1            ! Salinity limitation (ppt)
      op_dep_min         = p1            ! Light attenuates for optical depths exceeding min
      fr_graze_s         = p5            ! fraction of grazing spilled or slopped
      fr_graze_e         = p5            ! fraction of assimilation excreted 
      fr_mort2min        = p5            ! fractionation of mortality to Am
      fr_dFe             = 0.3_dbl_kind  ! fraction of remineralized nitrogen
                                         ! (in units of algal iron)
      k_nitrif           = c0            ! nitrification rate (1/day)           
      t_iron_conv        = 3065.0_dbl_kind ! desorption loss pFe to dFe (day)
      max_loss           = 0.9_dbl_kind ! restrict uptake to % of remaining value 
      max_dfe_doc1       = 0.2_dbl_kind ! max ratio of dFe to saccharides in the ice 
                                         !(nM Fe/muM C)    
      fr_resp_s          = 0.75_dbl_kind ! DMSPd fraction of respiration loss as DMSPd
      y_sk_DMS           = p5            ! fraction conversion given high yield
      t_sk_conv          = 3.0_dbl_kind  ! Stefels conversion time (d)
      t_sk_ox            = 10.0_dbl_kind ! DMS oxidation time (d)
      algaltype_diatoms  = c0            ! ------------------
      algaltype_sp       = p5            !
      algaltype_phaeo    = p5            !
      nitratetype        = -c1           ! mobility type between
      ammoniumtype       = c1            ! stationary <-->  mobile
      silicatetype       = -c1           !
      dmspptype          = p5            !
      dmspdtype          = -c1           !
      humtype            = c1            !
      doctype_s          = p5            !
      doctype_l          = p5            !
      dontype_protein    = p5            !
      fedtype_1          = p5            !
      feptype_1          = p5            !
      zaerotype_bc1      = c1            !
      zaerotype_bc2      = c1            !
      zaerotype_dust1    = c1            !
      zaerotype_dust2    = c1            !
      zaerotype_dust3    = c1            !
      zaerotype_dust4    = c1            !--------------------
      ratio_C2N_diatoms  = 7.0_dbl_kind  ! algal C to N ratio (mol/mol)
      ratio_C2N_sp       = 7.0_dbl_kind
      ratio_C2N_phaeo    = 7.0_dbl_kind
      ratio_chl2N_diatoms= 2.1_dbl_kind  ! algal chlorophyll to N ratio (mg/mmol)
      ratio_chl2N_sp     = 1.1_dbl_kind
      ratio_chl2N_phaeo  = 0.84_dbl_kind
      F_abs_chl_diatoms  = 2.0_dbl_kind  ! scales absorbed radiation for dEdd
      F_abs_chl_sp       = 4.0_dbl_kind
      F_abs_chl_phaeo    = 5.0
      ratio_C2N_proteins = 7.0_dbl_kind  ! ratio of C to N in proteins (mol/mol)       

      ! z salinity  parameters
      grid_oS         = c5            ! for bottom flux         
      l_skS           = 7.0_dbl_kind  ! characteristic diffusive scale (m)  

      !-----------------------------------------------------------------
      ! read from input file
      !-----------------------------------------------------------------

      call get_fileunit(nu_nml)

      if (my_task == master_task) then
         open (nu_nml, file=trim(nml_filename), status='old',iostat=nml_error)
         if (nml_error /= 0) then
            nml_error = -1
         else
            nml_error =  1
         endif 

         print*,'Reading zbgc_nml'
         do while (nml_error > 0)
            read(nu_nml, nml=zbgc_nml,iostat=nml_error)
         end do
         if (nml_error == 0) close(nu_nml)
      endif
      call broadcast_scalar(nml_error, master_task)
      if (nml_error /= 0) then
         call abort_ice(subname//'ERROR: reading zbgc namelist')
      endif
      call release_fileunit(nu_nml)

      !-----------------------------------------------------------------
      ! broadcast
      !-----------------------------------------------------------------

      call broadcast_scalar(solve_zsal,         master_task)  
      call broadcast_scalar(restart_zsal,       master_task)  
      call broadcast_scalar(tr_brine,           master_task)
      call broadcast_scalar(restart_hbrine,     master_task) 

      call broadcast_scalar(phi_snow,           master_task)
      call broadcast_scalar(grid_oS,            master_task)
      call broadcast_scalar(l_skS,              master_task)

      call broadcast_scalar(solve_zbgc,         master_task)
      call broadcast_scalar(skl_bgc,            master_task)
      call broadcast_scalar(restart_bgc,        master_task)
      call broadcast_scalar(bgc_flux_type,      master_task)
      call broadcast_scalar(restore_bgc,        master_task)
      call broadcast_scalar(tr_bgc_N,           master_task)
      call broadcast_scalar(tr_bgc_C,           master_task)
      call broadcast_scalar(tr_bgc_chl,         master_task)
      call broadcast_scalar(tr_bgc_Nit,         master_task)
      call broadcast_scalar(tr_bgc_Am,          master_task)
      call broadcast_scalar(tr_bgc_Sil,         master_task)
      call broadcast_scalar(tr_bgc_hum,         master_task)
      call broadcast_scalar(tr_bgc_DMS,         master_task) 
      call broadcast_scalar(tr_bgc_PON,         master_task) 
      call broadcast_scalar(tr_bgc_DON,         master_task) 
      call broadcast_scalar(tr_bgc_Fe,          master_task) 

      call broadcast_scalar(z_tracers,          master_task)
      call broadcast_scalar(tr_zaero,           master_task)
      call broadcast_scalar(dEdd_algae,         master_task) 
      call broadcast_scalar(modal_aero,         master_task)
      call broadcast_scalar(grid_o,             master_task)
      call broadcast_scalar(grid_o_t,           master_task)
      call broadcast_scalar(l_sk,               master_task)
      call broadcast_scalar(scale_bgc,          master_task)
      call broadcast_scalar(initbio_frac,       master_task)
      call broadcast_scalar(frazil_scav,        master_task)
      call broadcast_scalar(ratio_Si2N_diatoms, master_task)
      call broadcast_scalar(ratio_Si2N_sp,      master_task)
      call broadcast_scalar(ratio_Si2N_phaeo,   master_task)
      call broadcast_scalar(ratio_S2N_diatoms,  master_task)
      call broadcast_scalar(ratio_S2N_sp,       master_task)
      call broadcast_scalar(ratio_S2N_phaeo,    master_task)
      call broadcast_scalar(ratio_Fe2C_diatoms, master_task)
      call broadcast_scalar(ratio_Fe2C_sp,      master_task)
      call broadcast_scalar(ratio_Fe2C_phaeo,   master_task)
      call broadcast_scalar(ratio_Fe2N_diatoms, master_task)
      call broadcast_scalar(ratio_Fe2N_sp,      master_task)
      call broadcast_scalar(ratio_Fe2N_phaeo,   master_task)
      call broadcast_scalar(ratio_Fe2DON   ,    master_task)
      call broadcast_scalar(ratio_Fe2DOC_s ,    master_task)
      call broadcast_scalar(ratio_Fe2DOC_l ,    master_task)
      call broadcast_scalar(fr_resp     ,       master_task)
      call broadcast_scalar(tau_min  ,          master_task)
      call broadcast_scalar(tau_max  ,          master_task)
      call broadcast_scalar(algal_vel  ,        master_task)
      call broadcast_scalar(R_dFe2dust ,        master_task)
      call broadcast_scalar(dustFe_sol      ,   master_task)
      call broadcast_scalar(chlabs_diatoms ,  master_task)
      call broadcast_scalar(chlabs_sp      ,  master_task)
      call broadcast_scalar(chlabs_phaeo     ,  master_task)
      call broadcast_scalar(alpha2max_low_diatoms ,  master_task) 
      call broadcast_scalar(alpha2max_low_sp      ,  master_task)
      call broadcast_scalar(alpha2max_low_phaeo   ,  master_task)
      call broadcast_scalar(beta2max_diatoms ,  master_task) 
      call broadcast_scalar(beta2max_sp      ,  master_task) 
      call broadcast_scalar(beta2max_phaeo   ,  master_task) 
      call broadcast_scalar(mu_max_diatoms   ,  master_task) 
      call broadcast_scalar(mu_max_sp        ,  master_task) 
      call broadcast_scalar(mu_max_phaeo     ,  master_task) 
      call broadcast_scalar(grow_Tdep_diatoms,  master_task) 
      call broadcast_scalar(grow_Tdep_sp     ,  master_task) 
      call broadcast_scalar(grow_Tdep_phaeo  ,  master_task) 
      call broadcast_scalar(fr_graze_diatoms ,  master_task) 
      call broadcast_scalar(fr_graze_sp      ,  master_task) 
      call broadcast_scalar(fr_graze_phaeo   ,  master_task)  
      call broadcast_scalar(mort_pre_diatoms ,  master_task) 
      call broadcast_scalar(mort_pre_sp      ,  master_task)  
      call broadcast_scalar(mort_pre_phaeo   ,  master_task)  
      call broadcast_scalar(mort_Tdep_diatoms,  master_task) 
      call broadcast_scalar(mort_Tdep_sp     ,  master_task)   
      call broadcast_scalar(mort_Tdep_phaeo  ,  master_task)   
      call broadcast_scalar(k_exude_diatoms  ,  master_task) 
      call broadcast_scalar(k_exude_sp       ,  master_task)   
      call broadcast_scalar(k_exude_phaeo    ,  master_task)   
      call broadcast_scalar(K_Nit_diatoms    ,  master_task) 
      call broadcast_scalar(K_Nit_sp         ,  master_task)
      call broadcast_scalar(K_Nit_phaeo      ,  master_task)
      call broadcast_scalar(K_Am_diatoms     ,  master_task)
      call broadcast_scalar(K_Am_sp          ,  master_task)
      call broadcast_scalar(K_Am_phaeo       ,  master_task)
      call broadcast_scalar(K_Sil_diatoms    ,  master_task)
      call broadcast_scalar(K_Sil_sp         ,  master_task)
      call broadcast_scalar(K_Sil_phaeo      ,  master_task)
      call broadcast_scalar(K_Fe_diatoms     ,  master_task)
      call broadcast_scalar(K_Fe_sp          ,  master_task)
      call broadcast_scalar(K_Fe_phaeo       ,  master_task)
      call broadcast_scalar(f_don_protein    ,  master_task)
      call broadcast_scalar(kn_bac_protein   ,  master_task)
      call broadcast_scalar(f_don_Am_protein ,  master_task)
      call broadcast_scalar(f_doc_s          ,  master_task)
      call broadcast_scalar(f_doc_l          ,  master_task)  
      call broadcast_scalar(f_exude_s        ,  master_task)
      call broadcast_scalar(f_exude_l        ,  master_task)
      call broadcast_scalar(k_bac_s          ,  master_task) 
      call broadcast_scalar(k_bac_l          ,  master_task)
      call broadcast_scalar(T_max            ,  master_task)
      call broadcast_scalar(fsal             ,  master_task)
      call broadcast_scalar(op_dep_min       ,  master_task)
      call broadcast_scalar(fr_graze_s       ,  master_task) 
      call broadcast_scalar(fr_graze_e       ,  master_task) 
      call broadcast_scalar(fr_mort2min      ,  master_task) 
      call broadcast_scalar(fr_dFe           ,  master_task)
      call broadcast_scalar(k_nitrif         ,  master_task)
      call broadcast_scalar(t_iron_conv      ,  master_task)
      call broadcast_scalar(max_loss         ,  master_task)
      call broadcast_scalar(max_dfe_doc1     ,  master_task)
      call broadcast_scalar(fr_resp_s        ,  master_task)
      call broadcast_scalar(y_sk_DMS         ,  master_task)
      call broadcast_scalar(t_sk_conv        ,  master_task) 
      call broadcast_scalar(t_sk_ox          ,  master_task)
      call broadcast_scalar(algaltype_diatoms,  master_task)
      call broadcast_scalar(algaltype_sp       ,  master_task) 
      call broadcast_scalar(algaltype_phaeo    ,  master_task) 
      call broadcast_scalar(nitratetype        ,  master_task)
      call broadcast_scalar(ammoniumtype       ,  master_task)
      call broadcast_scalar(silicatetype       ,  master_task)
      call broadcast_scalar(dmspptype          ,  master_task)  
      call broadcast_scalar(dmspdtype          ,  master_task) 
      call broadcast_scalar(humtype            ,  master_task)
      call broadcast_scalar(doctype_s          ,  master_task) 
      call broadcast_scalar(doctype_l          ,  master_task)
      call broadcast_scalar(dontype_protein    ,  master_task)
      call broadcast_scalar(fedtype_1          ,  master_task)
      call broadcast_scalar(feptype_1          ,  master_task)
      call broadcast_scalar(zaerotype_bc1      ,  master_task)
      call broadcast_scalar(zaerotype_bc2      ,  master_task)
      call broadcast_scalar(zaerotype_dust1    ,  master_task)
      call broadcast_scalar(zaerotype_dust2    ,  master_task)
      call broadcast_scalar(zaerotype_dust3    ,  master_task)
      call broadcast_scalar(zaerotype_dust4    ,  master_task)
      call broadcast_scalar(ratio_C2N_diatoms  ,  master_task)
      call broadcast_scalar(ratio_C2N_sp       ,  master_task)
      call broadcast_scalar(ratio_C2N_phaeo    ,  master_task)
      call broadcast_scalar(ratio_chl2N_diatoms,  master_task)
      call broadcast_scalar(ratio_chl2N_sp     ,  master_task)
      call broadcast_scalar(ratio_chl2N_phaeo  ,  master_task)
      call broadcast_scalar(F_abs_chl_diatoms  ,  master_task)
      call broadcast_scalar(F_abs_chl_sp       ,  master_task)
      call broadcast_scalar(F_abs_chl_phaeo    ,  master_task)
      call broadcast_scalar(ratio_C2N_proteins ,  master_task) 

      !-----------------------------------------------------------------
      ! zsalinity and brine
      !-----------------------------------------------------------------

      if (.not.restart) then
         if (my_task == master_task) &
            write(nu_diag,*) subname//' WARNING: restart = false, setting bgc restart flags to false'
         restart_bgc =  .false.
         restart_hbrine =  .false.
         restart_zsal =  .false.
      endif

      if (solve_zsal .and. nblyr < 1)  then
         if (my_task == master_task) then
            write(nu_diag,*) subname,' ERROR: solve_zsal=T but 0 zsalinity tracers'
         endif
         abort_flag = 101
      endif 

      if (solve_zsal .and. ((.not. tr_brine) .or. (ktherm /= 1))) then
         if (my_task == master_task) then
            write(nu_diag,*) subname,' ERROR: solve_zsal needs tr_brine=T and ktherm=1'
         endif
         abort_flag = 102
      endif

      if (tr_brine .and. nblyr < 1 ) then
         if (my_task == master_task) then
            write(nu_diag,*) subname,' ERROR: tr_brine=T but no biology layers compiled'
         endif
         abort_flag = 103
      endif 

      !-----------------------------------------------------------------
      ! biogeochemistry
      !-----------------------------------------------------------------

      if (.not. tr_brine) then
         if (solve_zbgc) then
            if (my_task == master_task) then
               write(nu_diag,*) subname,' ERROR: tr_brine = F and solve_zbgc = T'
            endif
            abort_flag = 104
         endif
         if (tr_zaero) then
            if (my_task == master_task) then
               write(nu_diag,*) subname,' ERROR: tr_brine = F and tr_zaero = T'
            endif
            abort_flag = 105
         endif
      endif

      if ((skl_bgc .AND. solve_zbgc) .or. (skl_bgc .AND. z_tracers)) then
         if (my_task == master_task) then
            write(nu_diag,*) subname,' ERROR: skl_bgc and solve_zbgc or z_tracers are both true'
         endif
         abort_flag = 106
      endif

      if (skl_bgc .AND. tr_zaero) then
         if (my_task == master_task) then
            write(nu_diag,*) subname,' ERROR: skl_bgc does not use vertical tracers'
         endif
         abort_flag = 107
      endif

      if (dEdd_algae .AND. trim(shortwave) /= 'dEdd') then 
         if (my_task == master_task) then
            write(nu_diag,*) subname,' ERROR: dEdd_algae = T but shortwave /= dEdd'
         endif
         abort_flag = 108
      endif

      if (dEdd_algae .AND. (.NOT. tr_bgc_N) .AND. (.NOT. tr_zaero)) then 
         if (my_task == master_task) then
            write(nu_diag,*) subname,' ERROR: need tr_bgc_N or tr_zaero for dEdd_algae'
         endif
         abort_flag = 109
      endif

      if (modal_aero .AND. (.NOT. tr_zaero) .AND. (.NOT. tr_aero)) then
         if (my_task == master_task) then
            write(nu_diag,*) subname,' ERROR: modal_aero T with tr_zaero and tr_aero'
         endif
         abort_flag = 110
      endif
         
      if (modal_aero .AND. trim(shortwave) /= 'dEdd') then 
         if (my_task == master_task) then
            write(nu_diag,*) subname,' ERROR: modal_aero = T but shortwave /= dEdd'
         endif
         abort_flag = 111
      endif
      if (n_algae > icepack_max_algae) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//'ERROR: number of algal types exceeds icepack_max_algae'
         endif
         abort_flag = 112
      endif
      if (n_doc > icepack_max_doc) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//'ERROR: number of doc types exceeds icepack_max_doc'
         endif
         abort_flag = 113
      endif
      if (n_dic > icepack_max_doc) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//'ERROR: number of dic types exceeds icepack_max_dic'
         endif
         abort_flag = 114
      endif
      if (n_don > icepack_max_don) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//'ERROR: number of don types exceeds icepack_max_don'
         endif
         abort_flag = 115
      endif
      if (n_fed  > icepack_max_fe ) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//'ERROR: number of dissolved fe types exceeds icepack_max_fe '
         endif
         abort_flag = 116
      endif
      if (n_fep  > icepack_max_fe ) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//'ERROR: number of particulate fe types exceeds icepack_max_fe '
         endif
         abort_flag = 117
      endif

      if (n_algae == 0 .and. skl_bgc) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//'ERROR: skl_bgc=T but 0 bgc or algal tracers compiled'
         endif
         abort_flag = 118
      endif

      if (n_algae == 0 .and. solve_zbgc) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//'ERROR: solve_zbgc=T but 0 zbgc or algal tracers compiled'
         endif
         abort_flag = 119
      endif

      if (solve_zbgc .and. .not. z_tracers) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//'ERROR: solve_zbgc=T but not z_tracers'
         endif
         abort_flag = 120
      endif

      if (skl_bgc .or. solve_zbgc) then
         if (.not. tr_bgc_N) then
            if (my_task == master_task) then
               write(nu_diag,*) subname//'ERROR: tr_bgc_N must be on for bgc'
            endif
            abort_flag = 121
         endif
         if (.not. tr_bgc_Nit) then
            if (my_task == master_task) then
               write(nu_diag,*) subname//'ERROR: tr_bgc_Nit must be on for bgc'
            endif
            abort_flag = 122
         endif
      else
         ! tcraig, allow bgc to be turned off in this case?
         tr_bgc_N         = .false.
         tr_bgc_C         = .false.
         tr_bgc_chl       = .false.
         tr_bgc_Nit       = .false.
         tr_bgc_Am        = .false.
         tr_bgc_Sil       = .false.
         tr_bgc_hum       = .false.
         tr_bgc_DMS       = .false.
         tr_bgc_PON       = .false.
         tr_bgc_DON       = .false.
         tr_bgc_Fe        = .false.
      endif

      !-----------------------------------------------------------------
      ! z layer aerosols
      !-----------------------------------------------------------------
      if (tr_zaero .and. .not. z_tracers) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//'ERROR: tr_zaero and not z_tracers'
         endif
         abort_flag = 123
      endif

      if (n_zaero > icepack_max_aero) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//'ERROR: number of z aerosols exceeds icepack_max_aero'
         endif
         abort_flag = 124
      endif

      !-----------------------------------------------------------------
      ! output
      !-----------------------------------------------------------------

      if (my_task == master_task) then
         write(nu_diag,1010) ' tr_brine                  = ', tr_brine
         if (tr_brine) then
         write(nu_diag,1010) ' restart_hbrine            = ', restart_hbrine
         write(nu_diag,1005) ' phi_snow                  = ', phi_snow
         endif
         write(nu_diag,1010) ' solve_zsal                = ', solve_zsal
         if (solve_zsal) then
         write(nu_diag,1010) ' restart_zsal              = ', restart_zsal
         write(nu_diag,1000) ' grid_oS                   = ', grid_oS
         write(nu_diag,1005) ' l_skS                     = ', l_skS
         endif

         write(nu_diag,1010) ' skl_bgc                   = ', skl_bgc
         write(nu_diag,1010) ' restart_bgc               = ', restart_bgc
         write(nu_diag,1010) ' tr_bgc_N                  = ', tr_bgc_N
         write(nu_diag,1010) ' tr_bgc_C                  = ', tr_bgc_C
         write(nu_diag,1010) ' tr_bgc_chl                = ', tr_bgc_chl
         write(nu_diag,1010) ' tr_bgc_Nit                = ', tr_bgc_Nit
         write(nu_diag,1010) ' tr_bgc_Am                 = ', tr_bgc_Am
         write(nu_diag,1010) ' tr_bgc_Sil                = ', tr_bgc_Sil
         write(nu_diag,1010) ' tr_bgc_hum                = ', tr_bgc_hum
         write(nu_diag,1010) ' tr_bgc_DMS                = ', tr_bgc_DMS
         write(nu_diag,1010) ' tr_bgc_PON                = ', tr_bgc_PON
         write(nu_diag,1010) ' tr_bgc_DON                = ', tr_bgc_DON
         write(nu_diag,1010) ' tr_bgc_Fe                 = ', tr_bgc_Fe
         write(nu_diag,1020) ' n_aero                    = ', n_aero
         write(nu_diag,1020) ' n_zaero                   = ', n_zaero
         write(nu_diag,1020) ' n_algae                   = ', n_algae
         write(nu_diag,1020) ' n_doc                     = ', n_doc
         write(nu_diag,1020) ' n_dic                     = ', n_dic
         write(nu_diag,1020) ' n_don                     = ', n_don
         write(nu_diag,1020) ' n_fed                     = ', n_fed
         write(nu_diag,1020) ' n_fep                     = ', n_fep

        if (skl_bgc) then

         write(nu_diag,1030) ' bgc_flux_type             = ', bgc_flux_type
         write(nu_diag,1010) ' restore_bgc               = ', restore_bgc

        elseif (z_tracers) then

         write(nu_diag,1010) ' dEdd_algae                = ', dEdd_algae
         write(nu_diag,1010) ' modal_aero                = ', modal_aero
         write(nu_diag,1010) ' scale_bgc                 = ', scale_bgc
         write(nu_diag,1010) ' solve_zbgc                = ', solve_zbgc
         write(nu_diag,1010) ' tr_zaero                  = ', tr_zaero
         write(nu_diag,1020) ' number of aerosols        = ', n_zaero
         ! bio parameters
         write(nu_diag,1000) ' grid_o                    = ', grid_o
         write(nu_diag,1000) ' grid_o_t                  = ', grid_o_t
         write(nu_diag,1005) ' l_sk                      = ', l_sk
         write(nu_diag,1000) ' initbio_frac              = ', initbio_frac
         write(nu_diag,1000) ' frazil_scav               = ', frazil_scav

        endif  ! skl_bgc or solve_bgc
      endif

      !-----------------------------------------------------------------
      ! abort if abort flag is set
      !-----------------------------------------------------------------

      if (abort_flag /= 0) then
        call flush_fileunit(nu_diag)
      endif
      call ice_barrier()
      if (abort_flag /= 0) then
         write(nu_diag,*) subname,' ERROR: abort_flag=',abort_flag
         call abort_ice (subname//' ABORTING on input ERRORS', &
            file=__FILE__, line=__LINE__)
      endif

      !-----------------------------------------------------------------
      ! set values in icepack
      !-----------------------------------------------------------------

      call icepack_init_parameters( &
          ktherm_in=ktherm, shortwave_in=shortwave, solve_zsal_in=solve_zsal, &
          skl_bgc_in=skl_bgc, z_tracers_in=z_tracers, scale_bgc_in=scale_bgc, &
          dEdd_algae_in=dEdd_algae, &
          solve_zbgc_in=solve_zbgc, &
          bgc_flux_type_in=bgc_flux_type, grid_o_in=grid_o, l_sk_in=l_sk, &
          initbio_frac_in=initbio_frac, &
          grid_oS_in=grid_oS, l_skS_in=l_skS, &
          phi_snow_in=phi_snow, frazil_scav_in = frazil_scav, &
          modal_aero_in=modal_aero)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__, line=__LINE__)

      call icepack_init_tracer_flags(tr_brine_in=tr_brine, &
         tr_bgc_Nit_in=tr_bgc_Nit, tr_bgc_Am_in =tr_bgc_Am,  tr_bgc_Sil_in=tr_bgc_Sil,   &
         tr_bgc_DMS_in=tr_bgc_DMS, tr_bgc_PON_in=tr_bgc_PON, &
         tr_bgc_N_in  =tr_bgc_N,   tr_bgc_C_in  =tr_bgc_C,   tr_bgc_chl_in=tr_bgc_chl,   &
         tr_bgc_DON_in=tr_bgc_DON, tr_bgc_Fe_in =tr_bgc_Fe,  tr_zaero_in  =tr_zaero,     &
         tr_bgc_hum_in=tr_bgc_hum)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__, line=__LINE__)

 1000    format (a30,2x,f9.2)  ! a30 to align formatted, unformatted statements
 1005    format (a30,2x,f9.6)  ! float
 1010    format (a30,2x,l6)    ! logical
 1020    format (a30,2x,i6)    ! integer
 1030    format (a30,   a8)    ! character

      end subroutine input_zbgc

!=======================================================================

! Count and index tracers
!
! author Elizabeth C. Hunke, LANL

      subroutine count_tracers

      use ice_domain_size, only: nilyr, nslyr, nblyr, nfsd, n_iso, &
          n_aero, n_zaero, n_algae, n_doc, n_dic, n_don, n_fed, n_fep

      ! local variables

      integer (kind=int_kind) :: &
         k, mm    , & ! loop index
         nk       , & ! layer index
         nk_bgc       ! layer index

      integer (kind=int_kind) :: ntrcr
      logical (kind=log_kind) :: tr_iage, tr_FY, tr_lvl, tr_pond, tr_aero, tr_fsd
      logical (kind=log_kind) :: tr_iso, tr_pond_cesm, tr_pond_lvl, tr_pond_topo
      integer (kind=int_kind) :: nt_Tsfc, nt_sice, nt_qice, nt_qsno, nt_iage, nt_FY
      integer (kind=int_kind) :: nt_alvl, nt_vlvl, nt_apnd, nt_hpnd, nt_ipnd, nt_aero
      integer (kind=int_kind) :: nt_fsd, nt_isosno, nt_isoice

      integer (kind=int_kind) :: &
         nbtrcr,        nbtrcr_sw,     &
         ntrcr_o,       nt_fbri,       &
         nt_bgc_Nit,    nt_bgc_Am,     nt_bgc_Sil,   &
         nt_bgc_DMS,    nt_bgc_PON,    nt_bgc_S,     &
         nt_bgc_DMSPp,  nt_bgc_DMSPd,  &
         nt_zbgc_frac,  nlt_chl_sw,    &
         nlt_bgc_Nit,   nlt_bgc_Am,    nlt_bgc_Sil, &
         nlt_bgc_DMS,   nlt_bgc_DMSPp, nlt_bgc_DMSPd, &
         nlt_bgc_PON,   nt_bgc_hum,    nlt_bgc_hum

      integer (kind=int_kind), dimension(icepack_max_aero) :: &
         nlt_zaero_sw       ! points to aerosol in trcrn_sw

      integer (kind=int_kind), dimension(icepack_max_algae) :: &
         nlt_bgc_N      , & ! algae
         nlt_bgc_chl

      integer (kind=int_kind), dimension(icepack_max_doc) :: &
         nlt_bgc_DOC        ! disolved organic carbon

      integer (kind=int_kind), dimension(icepack_max_don) :: &
         nlt_bgc_DON        !

      integer (kind=int_kind), dimension(icepack_max_dic) :: &
         nlt_bgc_DIC        ! disolved inorganic carbon

      integer (kind=int_kind), dimension(icepack_max_fe) :: &
         nlt_bgc_Fed    , & !
         nlt_bgc_Fep        !

      integer (kind=int_kind), dimension(icepack_max_aero) :: &
         nlt_zaero          ! non-reacting layer aerosols

      integer (kind=int_kind), dimension(icepack_max_algae) :: &
         nt_bgc_N , & ! diatoms, phaeocystis, pico/small
         nt_bgc_chl   ! diatoms, phaeocystis, pico/small

      integer (kind=int_kind), dimension(icepack_max_doc) :: &
         nt_bgc_DOC      !  dissolved organic carbon

      integer (kind=int_kind), dimension(icepack_max_don) :: &
         nt_bgc_DON         !  dissolved organic nitrogen

      integer (kind=int_kind), dimension(icepack_max_dic) :: &
         nt_bgc_DIC         !  dissolved inorganic carbon

      integer (kind=int_kind), dimension(icepack_max_fe) :: &
         nt_bgc_Fed,     & !  dissolved iron
         nt_bgc_Fep        !  particulate iron

      integer (kind=int_kind), dimension(icepack_max_aero) :: &
         nt_zaero       !  black carbon and other aerosols

      logical (kind=log_kind) :: &
          tr_brine, &
          tr_bgc_Nit,    tr_bgc_Am,    tr_bgc_Sil,   &
          tr_bgc_DMS,    tr_bgc_PON,   &
          tr_bgc_N,      tr_bgc_C,     tr_bgc_chl,   &
          tr_bgc_DON,    tr_bgc_Fe,    tr_zaero,     &
          tr_bgc_hum
 
      logical (kind=log_kind) :: &
          solve_zsal, skl_bgc, z_tracers

      character(len=*), parameter :: subname='(count_tracers)'

      !-----------------------------------------------------------------

      call icepack_query_parameters( &
          solve_zsal_out=solve_zsal, &
          skl_bgc_out=skl_bgc, z_tracers_out=z_tracers)

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__, line=__LINE__)

      call icepack_query_tracer_flags(tr_iage_out=tr_iage, tr_FY_out=tr_FY, &
         tr_lvl_out=tr_lvl, tr_aero_out=tr_aero, tr_pond_out=tr_pond, &
         tr_pond_cesm_out=tr_pond_cesm, tr_pond_lvl_out=tr_pond_lvl, &
         tr_pond_topo_out=tr_pond_topo, tr_brine_out=tr_brine, tr_fsd_out=tr_fsd, &
         tr_iso_out=tr_iso, &
         tr_bgc_Nit_out=tr_bgc_Nit, tr_bgc_Am_out =tr_bgc_Am,  tr_bgc_Sil_out=tr_bgc_Sil,   &
         tr_bgc_DMS_out=tr_bgc_DMS, tr_bgc_PON_out=tr_bgc_PON, &
         tr_bgc_N_out  =tr_bgc_N,   tr_bgc_C_out  =tr_bgc_C,   tr_bgc_chl_out=tr_bgc_chl,   &
         tr_bgc_DON_out=tr_bgc_DON, tr_bgc_Fe_out =tr_bgc_Fe,  tr_zaero_out  =tr_zaero,     &
         tr_bgc_hum_out=tr_bgc_hum)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__, line=__LINE__)

      ntrcr = 0

      ntrcr = ntrcr + 1             ! count tracers, starting with Tsfc = 1
      nt_Tsfc = ntrcr               ! index tracers, starting with Tsfc = 1

      nt_qice = ntrcr + 1
      ntrcr = ntrcr + nilyr ! qice in nilyr layers

      nt_qsno = ntrcr + 1
      ntrcr = ntrcr + nslyr ! qsno in nslyr layers

      nt_sice = ntrcr + 1
      ntrcr = ntrcr + nilyr ! sice in nilyr layers

      nt_iage = 0
      if (tr_iage) then
          ntrcr = ntrcr + 1
          nt_iage = ntrcr   ! chronological ice age
      endif

      nt_FY = 0
      if (tr_FY) then
          ntrcr = ntrcr + 1
          nt_FY = ntrcr     ! area of first year ice
      endif

      nt_alvl = 0
      nt_vlvl = 0
      if (tr_lvl) then
          ntrcr = ntrcr + 1
          nt_alvl = ntrcr
          ntrcr = ntrcr + 1
          nt_vlvl = ntrcr
      endif

      nt_apnd = 0
      nt_hpnd = 0
      nt_ipnd = 0
      if (tr_pond) then            ! all explicit melt pond schemes
          ntrcr = ntrcr + 1
          nt_apnd = ntrcr
          ntrcr = ntrcr + 1
          nt_hpnd = ntrcr
          if (tr_pond_lvl) then
              ntrcr = ntrcr + 1    ! refrozen pond ice lid thickness
              nt_ipnd = ntrcr      ! on level-ice ponds (if frzpnd='hlid')
          endif
          if (tr_pond_topo) then
              ntrcr = ntrcr + 1    ! 
              nt_ipnd = ntrcr      ! refrozen pond ice lid thickness
          endif
      endif

      nt_fsd = 0
      if (tr_fsd) then
          nt_fsd = ntrcr + 1       ! floe size distribution
          ntrcr = ntrcr + nfsd
      endif

      nt_isosno = 0
      nt_isoice = 0
      if (tr_iso) then
          nt_isosno = ntrcr + 1    ! isotopes in snow
          ntrcr = ntrcr + n_iso
          nt_isoice = ntrcr + 1    ! isotopes in ice
          ntrcr = ntrcr + n_iso
      endif

      nt_aero = 0
      if (tr_aero) then
          nt_aero = ntrcr + 1
          ntrcr = ntrcr + 4*n_aero ! 4 dEdd layers, n_aero species
      else
!tcx, modify code so we don't have to reset n_aero here
          n_aero = 0       !echmod - this is not getting set correctly (overwritten later?)
      endif
              
      !-----------------------------------------------------------------
      ! initialize zbgc tracer indices
      !-----------------------------------------------------------------

      nbtrcr = 0
      nbtrcr_sw = 0
      nt_zbgc_frac = 0

      ! vectors of size icepack_max_algae
      nlt_bgc_N(:) = 0
      nlt_bgc_chl(:) = 0
      nt_bgc_N(:) = 0
      nt_bgc_chl(:) = 0

      ! vectors of size icepack_max_dic
      nlt_bgc_DIC(:) = 0
      nt_bgc_DIC(:) = 0

      ! vectors of size icepack_max_doc
      nlt_bgc_DOC(:) = 0
      nt_bgc_DOC(:) = 0

      ! vectors of size icepack_max_don
      nlt_bgc_DON(:) = 0
      nt_bgc_DON(:) = 0

      ! vectors of size icepack_max_fe
      nlt_bgc_Fed(:) = 0
      nlt_bgc_Fep(:) = 0
      nt_bgc_Fed(:) = 0
      nt_bgc_Fep(:) = 0

      ! vectors of size icepack_max_aero
      nlt_zaero(:) = 0
      nlt_zaero_sw(:) = 0
      nt_zaero(:) = 0

      nlt_bgc_Nit    = 0
      nlt_bgc_Am     = 0
      nlt_bgc_Sil    = 0
      nlt_bgc_DMSPp  = 0
      nlt_bgc_DMSPd  = 0
      nlt_bgc_DMS    = 0
      nlt_bgc_PON    = 0
      nlt_bgc_hum    = 0
!      nlt_bgc_C      = 0
      nlt_chl_sw     = 0

      nt_bgc_Nit    = 0
      nt_bgc_Am     = 0
      nt_bgc_Sil    = 0
      nt_bgc_DMSPp  = 0
      nt_bgc_DMSPd  = 0
      nt_bgc_DMS    = 0
      nt_bgc_PON    = 0
      nt_bgc_hum    = 0
!      nt_bgc_C      = 0

      ntrcr_o = ntrcr
      nt_fbri = 0
      if (tr_brine) then
          nt_fbri = ntrcr + 1   ! ice volume fraction with salt
          ntrcr = ntrcr + 1
      endif

      nt_bgc_S = 0
      if (solve_zsal) then       ! .true. only if tr_brine = .true.
          nt_bgc_S = ntrcr + 1
          ntrcr = ntrcr + nblyr
      endif

      if (skl_bgc .or. z_tracers) then

         if (skl_bgc) then
            nk = 1
         elseif (z_tracers) then ! defined on nblyr+1 in ice
                                 ! and 2 snow layers (snow surface + interior)
            nk = nblyr + 1
         endif ! skl_bgc or z_tracers
         nk_bgc = nk                 ! number of bgc layers in ice
         if (nk > 1) nk_bgc = nk + 2 ! number of bgc layers in ice and snow

         !-----------------------------------------------------------------
         ! count tracers and assign tracer indices
         !-----------------------------------------------------------------

         if (tr_bgc_N) then
            do mm = 1, n_algae
               nt_bgc_N(mm) = ntrcr + 1
               do k = 1, nk_bgc
                  ntrcr = ntrcr + 1
               enddo
               nbtrcr = nbtrcr + 1
               nlt_bgc_N(mm) = nbtrcr
            enddo   ! mm
         endif ! tr_bgc_N

         if (tr_bgc_Nit) then
            nt_bgc_Nit = ntrcr + 1
            do k = 1, nk_bgc
               ntrcr = ntrcr + 1
            enddo
            nbtrcr = nbtrcr + 1
            nlt_bgc_Nit = nbtrcr
         endif ! tr_bgc_Nit

         if (tr_bgc_C) then
          !
          ! Algal C is not yet distinct from algal N
          ! * Reqires exudation and/or changing C:N ratios
          ! for implementation
          !
          !  do mm = 1,n_algae
          !     nt_bgc_C(mm) = ntrcr + 1
          !     do k = 1, nk_bgc
          !        ntrcr = ntrcr + 1
          !     enddo
          !     nbtrcr = nbtrcr + 1
          !     nlt_bgc_C(mm) = nbtrcr
          !  enddo   ! mm

            do mm = 1, n_doc
               nt_bgc_DOC(mm) = ntrcr + 1
               do k = 1, nk_bgc
                  ntrcr = ntrcr + 1
               enddo
               nbtrcr = nbtrcr + 1
               nlt_bgc_DOC(mm) = nbtrcr
            enddo   ! mm
            do mm = 1, n_dic
               nt_bgc_DIC(mm) = ntrcr + 1
               do k = 1, nk_bgc
                  ntrcr = ntrcr + 1
               enddo
               nbtrcr = nbtrcr + 1
               nlt_bgc_DIC(mm) = nbtrcr
            enddo   ! mm
         endif      ! tr_bgc_C

         if (tr_bgc_chl) then
            do mm = 1, n_algae
               nt_bgc_chl(mm) = ntrcr + 1
               do k = 1, nk_bgc
                  ntrcr = ntrcr + 1
               enddo
               nbtrcr = nbtrcr + 1
               nlt_bgc_chl(mm) = nbtrcr
            enddo   ! mm
         endif      ! tr_bgc_chl

         if (tr_bgc_Am) then
            nt_bgc_Am = ntrcr + 1
            do k = 1, nk_bgc
               ntrcr = ntrcr + 1
            enddo
            nbtrcr = nbtrcr + 1
            nlt_bgc_Am = nbtrcr
         endif
         if (tr_bgc_Sil) then
            nt_bgc_Sil = ntrcr + 1
            do k = 1, nk_bgc
               ntrcr = ntrcr + 1
            enddo
            nbtrcr = nbtrcr + 1
            nlt_bgc_Sil = nbtrcr
         endif

         if (tr_bgc_DMS) then   ! all together
            nt_bgc_DMSPp = ntrcr + 1
            do k = 1, nk_bgc
               ntrcr = ntrcr + 1
            enddo
            nbtrcr = nbtrcr + 1
            nlt_bgc_DMSPp = nbtrcr

            nt_bgc_DMSPd = ntrcr + 1
            do k = 1, nk_bgc
               ntrcr = ntrcr + 1
            enddo
            nbtrcr = nbtrcr + 1
            nlt_bgc_DMSPd = nbtrcr

            nt_bgc_DMS = ntrcr + 1
            do k = 1, nk_bgc
               ntrcr = ntrcr + 1
            enddo
            nbtrcr = nbtrcr + 1
            nlt_bgc_DMS = nbtrcr
         endif

         if (tr_bgc_PON) then
            nt_bgc_PON = ntrcr + 1
            do k = 1, nk_bgc
               ntrcr = ntrcr + 1
            enddo
            nbtrcr = nbtrcr + 1
            nlt_bgc_PON = nbtrcr
         endif

         if (tr_bgc_DON) then
            do mm = 1, n_don
               nt_bgc_DON(mm) = ntrcr + 1
               do k = 1, nk_bgc
                  ntrcr = ntrcr + 1
               enddo
               nbtrcr = nbtrcr + 1
               nlt_bgc_DON(mm) = nbtrcr
            enddo   ! mm
         endif      ! tr_bgc_DON

         if (tr_bgc_Fe) then
            do mm = 1, n_fed
               nt_bgc_Fed(mm) = ntrcr + 1
               do k = 1, nk_bgc
                  ntrcr = ntrcr + 1
               enddo
               nbtrcr = nbtrcr + 1
               nlt_bgc_Fed(mm) = nbtrcr
            enddo   ! mm
            do mm = 1, n_fep
               nt_bgc_Fep(mm) = ntrcr + 1
               do k = 1, nk_bgc
                  ntrcr = ntrcr + 1
               enddo
               nbtrcr = nbtrcr + 1
               nlt_bgc_Fep(mm) = nbtrcr
            enddo   ! mm
         endif      ! tr_bgc_Fe

         if (tr_bgc_hum) then
            nt_bgc_hum = ntrcr + 1
            do k = 1, nk_bgc
               ntrcr = ntrcr + 1
            enddo
            nbtrcr = nbtrcr + 1
            nlt_bgc_hum = nbtrcr
         endif

      endif ! skl_bgc .or. z_tracers

      if (z_tracers) then ! defined on nblyr+1 in ice
                          ! and 2 snow layers (snow surface + interior)
         ! z layer aerosols
         if (tr_zaero) then
            do mm = 1, n_zaero
               nt_zaero(mm) = ntrcr + 1
               do k = 1, nk_bgc
                  ntrcr = ntrcr + 1
               enddo
               nbtrcr = nbtrcr + 1
               nlt_zaero(mm) = nbtrcr
            enddo   ! mm
         endif      ! tr_zaero

         if (nbtrcr > 0) then
            nt_zbgc_frac = ntrcr + 1
            ntrcr = ntrcr + nbtrcr
         endif
      endif ! z_tracers

!tcx, +1 here is the unused tracer, want to get rid of it
      ntrcr = ntrcr + 1

!tcx, reset unusaed tracer index, eventually get rid of it.
      if (nt_iage  <= 0) nt_iage  = ntrcr
      if (nt_FY    <= 0) nt_FY    = ntrcr
      if (nt_alvl  <= 0) nt_alvl  = ntrcr
      if (nt_vlvl  <= 0) nt_vlvl  = ntrcr
      if (nt_apnd  <= 0) nt_apnd  = ntrcr
      if (nt_hpnd  <= 0) nt_hpnd  = ntrcr
      if (nt_ipnd  <= 0) nt_ipnd  = ntrcr
      if (nt_fsd   <= 0) nt_fsd   = ntrcr
      if (nt_isosno<= 0) nt_isosno= ntrcr
      if (nt_isoice<= 0) nt_isoice= ntrcr
      if (nt_aero  <= 0) nt_aero  = ntrcr
      if (nt_fbri  <= 0) nt_fbri  = ntrcr
      if (nt_bgc_S <= 0) nt_bgc_S = ntrcr

      if (my_task == master_task) then
         write(nu_diag,*) ' '
         write(nu_diag,1020) ' ntrcr                     = ', ntrcr
         write(nu_diag,1020) ' nbtrcr                    = ', nbtrcr
         write(nu_diag,1020) ' nbtrcr_sw                 = ', nbtrcr_sw
         write(nu_diag,*) ' '
         write(nu_diag,1020) ' nt_sice                   = ', nt_sice
         write(nu_diag,1020) ' nt_qice                   = ', nt_qice
         write(nu_diag,1020) ' nt_qsno                   = ', nt_qsno
         write(nu_diag,*)' '
 1020    format (a30,2x,i6)     ! integer
         call flush_fileunit(nu_diag)
      endif                     ! my_task = master_task
      call icepack_init_tracer_sizes(ntrcr_in=ntrcr, &
         ntrcr_o_in=ntrcr_o, nbtrcr_in=nbtrcr, nbtrcr_sw_in=nbtrcr_sw)
      call icepack_init_tracer_indices(nt_Tsfc_in=nt_Tsfc, nt_sice_in=nt_sice, &
         nt_qice_in=nt_qice, nt_qsno_in=nt_qsno, nt_iage_in=nt_iage, nt_fy_in=nt_fy, &
         nt_alvl_in=nt_alvl, nt_vlvl_in=nt_vlvl, nt_apnd_in=nt_apnd, nt_hpnd_in=nt_hpnd, &
         nt_ipnd_in=nt_ipnd, nt_fsd_in=nt_fsd, nt_aero_in=nt_aero, &
         nt_isosno_in=nt_isosno,     nt_isoice_in=nt_isoice,       nt_fbri_in=nt_fbri,      &
         nt_bgc_Nit_in=nt_bgc_Nit,   nt_bgc_Am_in=nt_bgc_Am,       nt_bgc_Sil_in=nt_bgc_Sil,   &
         nt_bgc_DMS_in=nt_bgc_DMS,   nt_bgc_PON_in=nt_bgc_PON,     nt_bgc_S_in=nt_bgc_S,     &
         nt_bgc_N_in=nt_bgc_N,       nt_bgc_chl_in=nt_bgc_chl,   &
         nt_bgc_DOC_in=nt_bgc_DOC,   nt_bgc_DON_in=nt_bgc_DON,     nt_bgc_DIC_in=nt_bgc_DIC,   &
         nt_zaero_in=nt_zaero,       nt_bgc_DMSPp_in=nt_bgc_DMSPp, nt_bgc_DMSPd_in=nt_bgc_DMSPd, &
         nt_bgc_Fed_in=nt_bgc_Fed,   nt_bgc_Fep_in=nt_bgc_Fep,     nt_zbgc_frac_in=nt_zbgc_frac, &
         nlt_zaero_sw_in=nlt_zaero_sw,  nlt_chl_sw_in=nlt_chl_sw,  nlt_bgc_Sil_in=nlt_bgc_Sil, &
         nlt_bgc_N_in=nlt_bgc_N,     nlt_bgc_Nit_in=nlt_bgc_Nit,   nlt_bgc_Am_in=nlt_bgc_Am, &
         nlt_bgc_DMS_in=nlt_bgc_DMS, nlt_bgc_DMSPp_in=nlt_bgc_DMSPp, nlt_bgc_DMSPd_in=nlt_bgc_DMSPd, &
         nlt_zaero_in=nlt_zaero,     nlt_bgc_chl_in=nlt_bgc_chl, &
         nlt_bgc_DIC_in=nlt_bgc_DIC, nlt_bgc_DOC_in=nlt_bgc_DOC,   nlt_bgc_PON_in=nlt_bgc_PON, &
         nlt_bgc_DON_in=nlt_bgc_DON, nlt_bgc_Fed_in=nlt_bgc_Fed,   nlt_bgc_Fep_in=nlt_bgc_Fep, &
         nt_bgc_hum_in=nt_bgc_hum,   nlt_bgc_hum_in=nlt_bgc_hum)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname//' Icepack Abort2', &
         file=__FILE__, line=__LINE__)

      if (my_task == master_task) then
         call icepack_write_tracer_flags(nu_diag)
         call icepack_write_tracer_sizes(nu_diag)
         call icepack_write_tracer_indices(nu_diag)
      endif
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname//' Icepack Abort3', &
         file=__FILE__, line=__LINE__)

      end subroutine count_tracers

!=======================================================================

! Initialize vertical biogeochemistry
!
! author Elizabeth C. Hunke, LANL
!        Nicole Jeffery, LANL

      subroutine init_zbgc

      use ice_state, only: trcr_base, trcr_depend, n_trcr_strata, &
          nt_strata
      use ice_arrays_column, only: R_C2N, R_chl2N, R_C2N_DON, R_Si2N

      integer (kind=int_kind) :: &
         nbtrcr,        nbtrcr_sw,     nt_fbri,       &
         nt_bgc_Nit,    nt_bgc_Am,     nt_bgc_Sil,    &
         nt_bgc_DMS,    nt_bgc_PON,    nt_bgc_S,      &
         nt_bgc_DMSPp,  nt_bgc_DMSPd,                 &
         nt_zbgc_frac,  nlt_chl_sw,                   &
         nlt_bgc_Nit,   nlt_bgc_Am,    nlt_bgc_Sil,   &
         nlt_bgc_DMS,   nlt_bgc_DMSPp, nlt_bgc_DMSPd, &
         nlt_bgc_PON,   nt_bgc_hum,    nlt_bgc_hum

      integer (kind=int_kind), dimension(icepack_max_aero) :: &
         nlt_zaero_sw       ! points to aerosol in trcrn_sw

      integer (kind=int_kind), dimension(icepack_max_algae) :: &
         nlt_bgc_N      , & ! algae
         nlt_bgc_chl

      integer (kind=int_kind), dimension(icepack_max_doc) :: &
         nlt_bgc_DOC        ! disolved organic carbon

      integer (kind=int_kind), dimension(icepack_max_don) :: &
         nlt_bgc_DON        !

      integer (kind=int_kind), dimension(icepack_max_dic) :: &
         nlt_bgc_DIC        ! disolved inorganic carbon

      integer (kind=int_kind), dimension(icepack_max_fe) :: &
         nlt_bgc_Fed    , & !
         nlt_bgc_Fep        !

      integer (kind=int_kind), dimension(icepack_max_aero) :: &
         nlt_zaero          ! non-reacting layer aerosols

      integer (kind=int_kind), dimension(icepack_max_algae) :: &
         nt_bgc_N , & ! diatoms, phaeocystis, pico/small
         nt_bgc_chl   ! diatoms, phaeocystis, pico/small

      integer (kind=int_kind), dimension(icepack_max_doc) :: &
         nt_bgc_DOC      !  dissolved organic carbon

      integer (kind=int_kind), dimension(icepack_max_don) :: &
         nt_bgc_DON         !  dissolved organic nitrogen

      integer (kind=int_kind), dimension(icepack_max_dic) :: &
         nt_bgc_DIC         !  dissolved inorganic carbon

      integer (kind=int_kind), dimension(icepack_max_fe) :: &
         nt_bgc_Fed,     & !  dissolved iron
         nt_bgc_Fep        !  particulate iron

      integer (kind=int_kind), dimension(icepack_max_aero) :: &
         nt_zaero       !  black carbon and other aerosols

      integer (kind=int_kind), dimension(icepack_max_nbtrcr) :: &
         bio_index_o         ! relates nlt_bgc_NO to ocean concentration index

      integer (kind=int_kind), dimension(icepack_max_nbtrcr) :: &
         bio_index           ! relates bio indices, ie.  nlt_bgc_N to nt_bgc_N

      logical (kind=log_kind) :: &
         tr_brine, &
         tr_bgc_Nit,    tr_bgc_Am,    tr_bgc_Sil,   &
         tr_bgc_DMS,    tr_bgc_PON,   &
         tr_bgc_N,      tr_bgc_C,     tr_bgc_chl,   &
         tr_bgc_DON,    tr_bgc_Fe,    tr_zaero,     &
         tr_bgc_hum

      real (kind=dbl_kind) :: &
         initbio_frac, &
         frazil_scav

      real (kind=dbl_kind), dimension(icepack_max_nbtrcr) :: &
         zbgc_frac_init,&! initializes mobile fraction
         bgc_tracer_type ! described tracer in mobile or stationary phases
                         ! < 0 is purely mobile (eg. nitrate)
                         ! > 0 has timescales for transitions between
                         ! phases based on whether the ice is melting or growing

     real (kind=dbl_kind), dimension(icepack_max_nbtrcr) :: &
         zbgc_init_frac, &   ! fraction of ocean tracer  concentration in new ice
         tau_ret,        &   ! retention timescale  (s), mobile to stationary phase
         tau_rel             ! release timescale    (s), stationary to mobile phase

      logical (kind=log_kind) :: &
         skl_bgc, z_tracers, dEdd_algae, solve_zsal

      real (kind=dbl_kind), dimension(icepack_max_algae) :: &
         F_abs_chl          ! to scale absorption in Dedd

       real (kind=dbl_kind),  dimension(icepack_max_algae) :: &
         R_S2N      , & ! algal S to N (mole/mole)
         ! Marchetti et al 2006, 3 umol Fe/mol C for iron limited Pseudo-nitzschia
         R_Fe2C     , & ! algal Fe to carbon (umol/mmol)
         R_Fe2N         ! algal Fe to N (umol/mmol)

      real (kind=dbl_kind), dimension(icepack_max_don) :: &
         R_Fe2DON       ! Fe to N of DON (nmol/umol)

      real (kind=dbl_kind), dimension(icepack_max_doc) :: &
         R_Fe2DOC       ! Fe to C of DOC (nmol/umol)

      real (kind=dbl_kind), dimension(icepack_max_algae) :: &
         chlabs           , & ! chla absorption 1/m/(mg/m^3)
         alpha2max_low    , & ! light limitation (1/(W/m^2))
         beta2max         , & ! light inhibition (1/(W/m^2))
         mu_max           , & ! maximum growth rate (1/d)
         grow_Tdep        , & ! T dependence of growth (1/C)
         fr_graze         , & ! fraction of algae grazed
         mort_pre         , & ! mortality (1/day)
         mort_Tdep        , & ! T dependence of mortality (1/C)
         k_exude          , & ! algal carbon  exudation rate (1/d)
         K_Nit            , & ! nitrate half saturation (mmol/m^3)
         K_Am             , & ! ammonium half saturation (mmol/m^3)
         K_Sil            , & ! silicon half saturation (mmol/m^3)
         K_Fe                 ! iron half saturation  or micromol/m^3

      real (kind=dbl_kind), dimension(icepack_max_DON) :: &
         f_don            , & ! fraction of spilled grazing to DON
         kn_bac           , & ! Bacterial degredation of DON (1/d)
         f_don_Am             ! fraction of remineralized DON to Am

      real (kind=dbl_kind), dimension(icepack_max_DOC) :: &
         f_doc            , & ! fraction of mort_N that goes to each doc pool
         f_exude          , & ! fraction of exuded carbon to each DOC pool
         k_bac                ! Bacterial degredation of DOC (1/d)

      integer (kind=int_kind) :: &
         k, mm    , & ! loop index
         nk       , & ! layer index
         ierr

      integer (kind=int_kind) :: &
        ntd      , & ! for tracer dependency calculation
        nt_depend

      character(len=*), parameter :: subname='(init_zbgc)'

      !------------------------------------------------------------
      !        Tracers have mobile and stationary phases.
      ! ice growth allows for retention, ice melt facilitates mobility
      ! bgc_tracer_type defines the exchange timescales between these phases
      ! -1 : entirely in the mobile phase, no exchange  (this is the default)
      !  0 : retention time scale is tau_min, release time scale is tau_max
      !  1 : retention time scale is tau_max, release time scale is tau_min
      ! 0.5: retention time scale is tau_min, release time scale is tau_min
      !  2 : retention time scale is tau_max, release time scale is tau_max
      ! tau_min and tau_max are defined in icepack_intfc.f90
      !------------------------------------------------------------

      !-----------------------------------------------------------------
      ! get values from icepack
      !-----------------------------------------------------------------

      call icepack_query_parameters( &
          solve_zsal_out=solve_zsal, &
          skl_bgc_out=skl_bgc, z_tracers_out=z_tracers, &
          dEdd_algae_out=dEdd_algae, &
          grid_o_out=grid_o, l_sk_out=l_sk, &
          initbio_frac_out=initbio_frac, &
          grid_oS_out=grid_oS, l_skS_out=l_skS, &
          phi_snow_out=phi_snow, frazil_scav_out = frazil_scav)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__, line=__LINE__)

      call icepack_query_tracer_sizes( &
          nbtrcr_out=nbtrcr, nbtrcr_sw_out=nbtrcr_sw)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__, line=__LINE__)

      call icepack_query_tracer_flags( &
          tr_brine_out =tr_brine, &
          tr_bgc_Nit_out=tr_bgc_Nit, tr_bgc_Am_out=tr_bgc_Am,  tr_bgc_Sil_out=tr_bgc_Sil,   &
          tr_bgc_DMS_out=tr_bgc_DMS, tr_bgc_PON_out=tr_bgc_PON, &
          tr_bgc_N_out =tr_bgc_N,   tr_bgc_C_out =tr_bgc_C,   tr_bgc_chl_out=tr_bgc_chl,   &
          tr_bgc_DON_out=tr_bgc_DON, tr_bgc_Fe_out=tr_bgc_Fe,  tr_zaero_out =tr_zaero,     &
          tr_bgc_hum_out=tr_bgc_hum)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__, line=__LINE__)

      call icepack_query_tracer_indices( &
          nt_fbri_out=nt_fbri,         &
          nt_bgc_Nit_out=nt_bgc_Nit,   nt_bgc_Am_out=nt_bgc_Am,       nt_bgc_Sil_out=nt_bgc_Sil,   &
          nt_bgc_DMS_out=nt_bgc_DMS,   nt_bgc_PON_out=nt_bgc_PON,     nt_bgc_S_out=nt_bgc_S,     &
          nt_bgc_N_out=nt_bgc_N,       nt_bgc_chl_out=nt_bgc_chl,   &
          nt_bgc_DOC_out=nt_bgc_DOC,   nt_bgc_DON_out=nt_bgc_DON,     nt_bgc_DIC_out=nt_bgc_DIC,   &
          nt_zaero_out=nt_zaero,       nt_bgc_DMSPp_out=nt_bgc_DMSPp, nt_bgc_DMSPd_out=nt_bgc_DMSPd, &
          nt_bgc_Fed_out=nt_bgc_Fed,   nt_bgc_Fep_out=nt_bgc_Fep,     nt_zbgc_frac_out=nt_zbgc_frac, &
          nlt_zaero_sw_out=nlt_zaero_sw,  nlt_chl_sw_out=nlt_chl_sw,  nlt_bgc_Sil_out=nlt_bgc_Sil, &
          nlt_bgc_N_out=nlt_bgc_N,     nlt_bgc_Nit_out=nlt_bgc_Nit,   nlt_bgc_Am_out=nlt_bgc_Am, &
          nlt_bgc_DMS_out=nlt_bgc_DMS, nlt_bgc_DMSPp_out=nlt_bgc_DMSPp, nlt_bgc_DMSPd_out=nlt_bgc_DMSPd, &
          nlt_zaero_out=nlt_zaero,     nlt_bgc_chl_out=nlt_bgc_chl, &
          nlt_bgc_DIC_out=nlt_bgc_DIC, nlt_bgc_DOC_out=nlt_bgc_DOC,   nlt_bgc_PON_out=nlt_bgc_PON, &
          nlt_bgc_DON_out=nlt_bgc_DON, nlt_bgc_Fed_out=nlt_bgc_Fed,   nlt_bgc_Fep_out=nlt_bgc_Fep, &
          nt_bgc_hum_out=nt_bgc_hum,   nlt_bgc_hum_out=nlt_bgc_hum)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__, line=__LINE__)

      !-----------------------------------------------------------------
      ! Define array parameters
      !-----------------------------------------------------------------

      allocate(          &
         R_C2N_DON(icepack_max_don), & ! carbon to nitrogen mole ratio of DON pool
         R_C2N(icepack_max_algae),   & ! algal C to N (mole/mole)
         R_chl2N(icepack_max_algae), & ! 3 algal chlorophyll to N (mg/mmol)
         R_Si2N(icepack_max_algae),  & ! silica to nitrogen mole ratio for algal groups
         stat=ierr)
      if (ierr/=0) call abort_ice(subname//' Out of Memory')

      R_Si2N(1) = ratio_Si2N_diatoms
      R_Si2N(2) = ratio_Si2N_sp
      R_Si2N(3) = ratio_Si2N_phaeo

      R_S2N(1) = ratio_S2N_diatoms
      R_S2N(2) = ratio_S2N_sp
      R_S2N(3) = ratio_S2N_phaeo

      R_Fe2C(1) = ratio_Fe2C_diatoms
      R_Fe2C(2) = ratio_Fe2C_sp
      R_Fe2C(3) = ratio_Fe2C_phaeo

      R_Fe2N(1) = ratio_Fe2N_diatoms
      R_Fe2N(2) = ratio_Fe2N_sp
      R_Fe2N(3) = ratio_Fe2N_phaeo

      R_C2N(1) = ratio_C2N_diatoms
      R_C2N(2) = ratio_C2N_sp
      R_C2N(3) = ratio_C2N_phaeo

      R_chl2N(1) = ratio_chl2N_diatoms
      R_chl2N(2) = ratio_chl2N_sp
      R_chl2N(3) = ratio_chl2N_phaeo

      F_abs_chl(1) = F_abs_chl_diatoms
      F_abs_chl(2) = F_abs_chl_sp
      F_abs_chl(3) = F_abs_chl_phaeo

      R_Fe2DON(1) = ratio_Fe2DON
      R_C2N_DON(1) = ratio_C2N_proteins

      R_Fe2DOC(1) = ratio_Fe2DOC_s
      R_Fe2DOC(2) = ratio_Fe2DOC_l
      R_Fe2DOC(3) = c0

      chlabs(1) = chlabs_diatoms
      chlabs(2) = chlabs_sp
      chlabs(3) = chlabs_phaeo

      alpha2max_low(1) = alpha2max_low_diatoms
      alpha2max_low(2) = alpha2max_low_sp
      alpha2max_low(3) = alpha2max_low_phaeo

      beta2max(1) = beta2max_diatoms
      beta2max(2) = beta2max_sp
      beta2max(3) = beta2max_phaeo

      mu_max(1) = mu_max_diatoms
      mu_max(2) = mu_max_sp
      mu_max(3) = mu_max_phaeo

      grow_Tdep(1) = grow_Tdep_diatoms
      grow_Tdep(2) = grow_Tdep_sp
      grow_Tdep(3) = grow_Tdep_phaeo

      fr_graze(1) = fr_graze_diatoms
      fr_graze(2) = fr_graze_sp
      fr_graze(3) = fr_graze_phaeo

      mort_pre(1) = mort_pre_diatoms
      mort_pre(2) = mort_pre_sp
      mort_pre(3) = mort_pre_phaeo

      mort_Tdep(1) = mort_Tdep_diatoms
      mort_Tdep(2) = mort_Tdep_sp
      mort_Tdep(3) = mort_Tdep_phaeo

      k_exude(1) = k_exude_diatoms
      k_exude(2) = k_exude_sp
      k_exude(3) = k_exude_phaeo

      K_Nit(1) = K_Nit_diatoms
      K_Nit(2) = K_Nit_sp
      K_Nit(3) = K_Nit_phaeo

      K_Am(1) = K_Am_diatoms
      K_Am(2) = K_Am_sp
      K_Am(3) = K_Am_phaeo

      K_Sil(1) = K_Sil_diatoms
      K_Sil(2) = K_Sil_sp
      K_Sil(3) = K_Sil_phaeo

      K_Fe(1) = K_Fe_diatoms
      K_Fe(2) = K_Fe_sp
      K_Fe(3) = K_Fe_phaeo

      f_don(1) = f_don_protein
      kn_bac(1) = kn_bac_protein
      f_don_Am(1) = f_don_Am_protein

      f_doc(1) = f_doc_s
      f_doc(2) = f_doc_l

      f_exude(1) = f_exude_s
      f_exude(2) = f_exude_l
      k_bac(1) = k_bac_s
      k_bac(2) = k_bac_l

      dictype(:) = -c1

      algaltype(1) = algaltype_diatoms
      algaltype(2) = algaltype_sp
      algaltype(3) = algaltype_phaeo

      doctype(1) = doctype_s
      doctype(2) = doctype_l

      dontype(1) = dontype_protein

      fedtype(1) = fedtype_1
      feptype(1) = feptype_1

      zaerotype(1) = zaerotype_bc1
      zaerotype(2) = zaerotype_bc2
      zaerotype(3) = zaerotype_dust1
      zaerotype(4) = zaerotype_dust2
      zaerotype(5) = zaerotype_dust3
      zaerotype(6) = zaerotype_dust4

      call icepack_init_zbgc ( &
         R_S2N_in=R_S2N, R_Fe2C_in=R_Fe2C, R_Fe2N_in=R_Fe2N, R_C2N_in=R_C2N, &
         R_chl2N_in=R_chl2N, F_abs_chl_in=F_abs_chl, R_Fe2DON_in=R_Fe2DON, R_Fe2DOC_in=R_Fe2DOC, &
         mort_Tdep_in=mort_Tdep, k_exude_in=k_exude, &
         K_Nit_in=K_Nit, K_Am_in=K_Am, K_sil_in=K_Sil, K_Fe_in=K_Fe, &
         f_don_in=f_don, kn_bac_in=kn_bac, f_don_Am_in=f_don_Am, f_exude_in=f_exude, k_bac_in=k_bac, &
         fr_resp_in=fr_resp, algal_vel_in=algal_vel, R_dFe2dust_in=R_dFe2dust, &
         dustFe_sol_in=dustFe_sol, T_max_in=T_max, fr_mort2min_in=fr_mort2min, fr_dFe_in=fr_dFe, &
         op_dep_min_in=op_dep_min, fr_graze_s_in=fr_graze_s, fr_graze_e_in=fr_graze_e, &
         k_nitrif_in=k_nitrif, t_iron_conv_in=t_iron_conv, max_loss_in=max_loss, max_dfe_doc1_in=max_dfe_doc1, &
         fr_resp_s_in=fr_resp_s, y_sk_DMS_in=y_sk_DMS, t_sk_conv_in=t_sk_conv, t_sk_ox_in=t_sk_ox, &
         mu_max_in=mu_max, R_Si2N_in=R_Si2N, R_C2N_DON_in=R_C2N_DON, chlabs_in=chlabs, &
         alpha2max_low_in=alpha2max_low, beta2max_in=beta2max, grow_Tdep_in=grow_Tdep, &
         fr_graze_in=fr_graze, mort_pre_in=mort_pre, f_doc_in=f_doc,fsal_in=fsal)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__, line=__LINE__)

      !-----------------------------------------------------------------
      ! assign tracer dependencies
      ! bgc_tracer_type: < 0  purely mobile , >= 0 stationary
      !------------------------------------------------------------------

      if (tr_brine) then
          trcr_depend(nt_fbri)   = 1   ! volume-weighted
          trcr_base  (nt_fbri,1) = c0  ! volume-weighted
          trcr_base  (nt_fbri,2) = c1  ! volume-weighted
          trcr_base  (nt_fbri,3) = c0  ! volume-weighted
          n_trcr_strata(nt_fbri) = 0
          nt_strata  (nt_fbri,1) = 0
          nt_strata  (nt_fbri,2) = 0
      endif

      ntd = 0                    ! if nt_fbri /= 0 then use fbri dependency
      if (nt_fbri == 0) ntd = -1 ! otherwise make tracers depend on ice volume

      if (solve_zsal) then       ! .true. only if tr_brine = .true.
          do k = 1,nblyr
             trcr_depend(nt_bgc_S + k - 1) = 2 + nt_fbri + ntd
             trcr_base  (nt_bgc_S,1) = c0  ! default: ice area
             trcr_base  (nt_bgc_S,2) = c1
             trcr_base  (nt_bgc_S,3) = c0
             n_trcr_strata(nt_bgc_S) = 1
             nt_strata(nt_bgc_S,1) = nt_fbri
             nt_strata(nt_bgc_S,2) = 0
          enddo
      endif

      bio_index(:)   = 0
      bio_index_o(:) = 0

      if (skl_bgc) then
         nk = 1
         nt_depend = 0
      elseif (z_tracers) then ! defined on nblyr+1 in ice
                              ! and 2 snow layers (snow surface + interior)
         nk = nblyr + 1
         nt_depend = 2 + nt_fbri + ntd
      endif ! skl_bgc or z_tracers

      if (skl_bgc .or. z_tracers) then

      if (tr_bgc_N) then
         do mm = 1, n_algae      
            call init_bgc_trcr(nk,              nt_fbri,       &
                               nt_bgc_N(mm),    nlt_bgc_N(mm), &
                               algaltype(mm),   nt_depend,     &
                               bgc_tracer_type, trcr_depend,   &
                               trcr_base,       n_trcr_strata, &
                               nt_strata,       bio_index)
            bio_index_o(nlt_bgc_N(mm)) = mm
         enddo   ! mm
      endif ! tr_bgc_N

      if (tr_bgc_Nit) then
            call init_bgc_trcr(nk,              nt_fbri,       &
                               nt_bgc_Nit,      nlt_bgc_Nit,   &
                               nitratetype,     nt_depend,     &
                               bgc_tracer_type, trcr_depend,   &
                               trcr_base,       n_trcr_strata, &
                               nt_strata,       bio_index)
            bio_index_o(nlt_bgc_Nit) = icepack_max_algae + 1
      endif ! tr_bgc_Nit
         
      if (tr_bgc_C) then
       !
       ! Algal C is not yet distinct from algal N
       ! * Reqires exudation and/or changing C:N ratios
       ! for implementation
       !
       !  do mm = 1,n_algae      
       !     call init_bgc_trcr(nk,              nt_fbri,       &
       !                        nt_bgc_C(mm),    nlt_bgc_C(mm), &
       !                        algaltype(mm),   nt_depend,     &
       !                        bgc_tracer_type, trcr_depend,   &
       !                        trcr_base,       n_trcr_strata, &
       !                        nt_strata,       bio_index)
       !     bio_index_o(nlt_bgc_C(mm)) = icepack_max_algae + 1 + mm
       !  enddo   ! mm

         do mm = 1, n_doc
            call init_bgc_trcr(nk,              nt_fbri,       &
                               nt_bgc_DOC(mm),  nlt_bgc_DOC(mm), &
                               doctype(mm),     nt_depend,     &
                               bgc_tracer_type, trcr_depend,   &
                               trcr_base,       n_trcr_strata, &
                               nt_strata,       bio_index)
            bio_index_o(nlt_bgc_DOC(mm)) = icepack_max_algae + 1 + mm
         enddo   ! mm
         do mm = 1, n_dic
            call init_bgc_trcr(nk,              nt_fbri,       &
                               nt_bgc_DIC(mm),  nlt_bgc_DIC(mm), &
                               dictype(mm),     nt_depend,     &
                               bgc_tracer_type, trcr_depend,   &
                               trcr_base,       n_trcr_strata, &
                               nt_strata,       bio_index)
            bio_index_o(nlt_bgc_DIC(mm)) = icepack_max_algae + icepack_max_doc + 1 + mm
         enddo   ! mm
      endif      ! tr_bgc_C

      if (tr_bgc_chl) then
         do mm = 1, n_algae
            call init_bgc_trcr(nk,              nt_fbri,       &
                               nt_bgc_chl(mm),  nlt_bgc_chl(mm), &
                               algaltype(mm),   nt_depend,     &
                               bgc_tracer_type, trcr_depend,   &
                               trcr_base,       n_trcr_strata, &
                               nt_strata,       bio_index)
            bio_index_o(nlt_bgc_chl(mm)) = icepack_max_algae + 1 + icepack_max_doc + icepack_max_dic + mm
         enddo   ! mm
      endif      ! tr_bgc_chl

      if (tr_bgc_Am) then
            call init_bgc_trcr(nk,              nt_fbri,       &
                               nt_bgc_Am,       nlt_bgc_Am,    &
                               ammoniumtype,    nt_depend,     &
                               bgc_tracer_type, trcr_depend,   &
                               trcr_base,       n_trcr_strata, &
                               nt_strata,       bio_index)
            bio_index_o(nlt_bgc_Am) = 2*icepack_max_algae + icepack_max_doc + icepack_max_dic + 2
      endif    
      if (tr_bgc_Sil) then
            call init_bgc_trcr(nk,              nt_fbri,       &
                               nt_bgc_Sil,      nlt_bgc_Sil,   &
                               silicatetype,    nt_depend,     &
                               bgc_tracer_type, trcr_depend,   &
                               trcr_base,       n_trcr_strata, &
                               nt_strata,       bio_index)
            bio_index_o(nlt_bgc_Sil) = 2*icepack_max_algae + icepack_max_doc + icepack_max_dic + 3
      endif    
      if (tr_bgc_DMS) then   ! all together
            call init_bgc_trcr(nk,              nt_fbri,       &
                               nt_bgc_DMSPp,    nlt_bgc_DMSPp, &
                               dmspptype,       nt_depend,     &
                               bgc_tracer_type, trcr_depend,   &
                               trcr_base,       n_trcr_strata, &
                               nt_strata,       bio_index)
            bio_index_o(nlt_bgc_DMSPp) = 2*icepack_max_algae + icepack_max_doc + icepack_max_dic + 4

            call init_bgc_trcr(nk,              nt_fbri,       &
                               nt_bgc_DMSPd,    nlt_bgc_DMSPd, &
                               dmspdtype,       nt_depend,     &
                               bgc_tracer_type, trcr_depend,   &
                               trcr_base,       n_trcr_strata, &
                               nt_strata,       bio_index)
            bio_index_o(nlt_bgc_DMSPd) = 2*icepack_max_algae + icepack_max_doc + icepack_max_dic + 5

            call init_bgc_trcr(nk,              nt_fbri,       &
                               nt_bgc_DMS,      nlt_bgc_DMS,   &
                               dmspdtype,       nt_depend,     &
                               bgc_tracer_type, trcr_depend,   &
                               trcr_base,       n_trcr_strata, &
                               nt_strata,       bio_index)
            bio_index_o(nlt_bgc_DMS) = 2*icepack_max_algae + icepack_max_doc + icepack_max_dic + 6
      endif    
      if (tr_bgc_PON) then
            call init_bgc_trcr(nk,              nt_fbri,       &
                               nt_bgc_PON,      nlt_bgc_PON, &
                               nitratetype,     nt_depend,     &
                               bgc_tracer_type, trcr_depend,   &
                               trcr_base,       n_trcr_strata, &
                               nt_strata,       bio_index)
            bio_index_o(nlt_bgc_PON) =  2*icepack_max_algae + icepack_max_doc + icepack_max_dic + 7
      endif
      if (tr_bgc_DON) then
         do mm = 1, n_don
            call init_bgc_trcr(nk,              nt_fbri,       &
                               nt_bgc_DON(mm),  nlt_bgc_DON(mm), &
                               dontype(mm),     nt_depend,     &
                               bgc_tracer_type, trcr_depend,   &
                               trcr_base,       n_trcr_strata, &
                               nt_strata,       bio_index)
            bio_index_o(nlt_bgc_DON(mm)) = 2*icepack_max_algae + icepack_max_doc + icepack_max_dic + 7 + mm
         enddo   ! mm
      endif      ! tr_bgc_DON
      if (tr_bgc_Fe) then
         do mm = 1, n_fed
            call init_bgc_trcr(nk,              nt_fbri,       &
                               nt_bgc_Fed(mm),  nlt_bgc_Fed(mm), &
                               fedtype(mm),     nt_depend,     &
                               bgc_tracer_type, trcr_depend,   &
                               trcr_base,       n_trcr_strata, &
                               nt_strata,       bio_index)
            bio_index_o(nlt_bgc_Fed(mm)) = 2*icepack_max_algae + icepack_max_doc + icepack_max_dic &
                                         + icepack_max_don + 7 + mm
         enddo   ! mm
         do mm = 1, n_fep
            call init_bgc_trcr(nk,              nt_fbri,       &
                               nt_bgc_Fep(mm),  nlt_bgc_Fep(mm), &
                               feptype(mm),     nt_depend,     &
                               bgc_tracer_type, trcr_depend,   &
                               trcr_base,       n_trcr_strata, &
                               nt_strata,       bio_index)
            bio_index_o(nlt_bgc_Fep(mm)) = 2*icepack_max_algae + icepack_max_doc + icepack_max_dic &
                                         + icepack_max_don + icepack_max_fe + 7 + mm
         enddo   ! mm
      endif      ! tr_bgc_Fe 
  
      if (tr_bgc_hum) then
            call init_bgc_trcr(nk,              nt_fbri,       &
                               nt_bgc_hum,      nlt_bgc_hum,   &
                               humtype,         nt_depend,     &
                               bgc_tracer_type, trcr_depend,   &
                               trcr_base,       n_trcr_strata, &
                               nt_strata,       bio_index)
            bio_index_o(nlt_bgc_hum) =   2*icepack_max_algae + icepack_max_doc + 8 + icepack_max_dic &
                                         + icepack_max_don + 2*icepack_max_fe + icepack_max_aero 
      endif
      endif  ! skl_bgc or z_tracers

      if (skl_bgc) then
         if (dEdd_algae) then
           nlt_chl_sw = 1
           nbtrcr_sw = nilyr+nslyr+2  ! only the bottom layer will be nonzero
         endif

      elseif (z_tracers) then ! defined on nblyr+1 in ice
                              ! and 2 snow layers (snow surface + interior)
         if (tr_bgc_N) then
            if (dEdd_algae) then
               nlt_chl_sw = 1
               nbtrcr_sw =  nilyr+nslyr+2
            endif
         endif ! tr_bgc_N
      endif ! skl_bgc or z_tracers

      if (z_tracers) then ! defined on nblyr+1 in ice
                          ! and 2 snow layers (snow surface + interior)

         nk = nblyr + 1
         nt_depend = 2 + nt_fbri + ntd 

         ! z layer aerosols
         if (tr_zaero) then
            do mm = 1, n_zaero
               if (dEdd_algae) then
                  nlt_zaero_sw(mm) = nbtrcr_sw + 1
                  nbtrcr_sw = nbtrcr_sw + nilyr + nslyr+2
               endif
               call init_bgc_trcr(nk,              nt_fbri,       &
                                  nt_zaero(mm),    nlt_zaero(mm), &
                                  zaerotype(mm),   nt_depend,     &
                                  bgc_tracer_type, trcr_depend,   &
                                  trcr_base,       n_trcr_strata, &
                                  nt_strata,       bio_index)
               bio_index_o(nlt_zaero(mm)) = 2*icepack_max_algae + icepack_max_doc + icepack_max_dic &
                                          + icepack_max_don + 2*icepack_max_fe + 7 + mm
            enddo   ! mm
         endif      ! tr_zaero

         if (nbtrcr > 0) then
            do k = 1,nbtrcr 
               zbgc_frac_init(k) = c1   
               trcr_depend(nt_zbgc_frac+k-1) =  2+nt_fbri 
               trcr_base(nt_zbgc_frac+ k - 1,1)  = c0
               trcr_base(nt_zbgc_frac+ k - 1,2)  = c1
               trcr_base(nt_zbgc_frac+ k - 1,3)  = c0
               n_trcr_strata(nt_zbgc_frac+ k - 1)= 1  
               nt_strata(nt_zbgc_frac+ k - 1,1)  = nt_fbri
               nt_strata(nt_zbgc_frac+ k - 1,2)  = 0 
               tau_ret(k) = c1
               tau_rel(k) = c1
               if (bgc_tracer_type(k) >=  c0 .and. bgc_tracer_type(k) < p5) then
                  tau_ret(k) = tau_min
                  tau_rel(k) = tau_max
                  zbgc_frac_init(k) = c1
               elseif (bgc_tracer_type(k) >= p5 .and. bgc_tracer_type(k) < c1) then
                  tau_ret(k) = tau_min
                  tau_rel(k) = tau_min
                  zbgc_frac_init(k) = c1
               elseif (bgc_tracer_type(k) >= c1 .and. bgc_tracer_type(k) < c2) then
                  tau_ret(k) = tau_max
                  tau_rel(k) = tau_min
                  zbgc_frac_init(k) = c1
               elseif (bgc_tracer_type(k) >= c2 ) then
                  tau_ret(k) = tau_max
                  tau_rel(k) = tau_max
                  zbgc_frac_init(k) = c1
               endif
            enddo
         endif

      endif ! z_tracers

      do k = 1, nbtrcr
         zbgc_init_frac(k) = frazil_scav
         if (bgc_tracer_type(k) < c0)  zbgc_init_frac(k) = initbio_frac
      enddo  

      !-----------------------------------------------------------------
      ! set values in icepack
      !-----------------------------------------------------------------

      call icepack_init_zbgc( &
         zbgc_init_frac_in=zbgc_init_frac, tau_ret_in=tau_ret, tau_rel_in=tau_rel, &
         zbgc_frac_init_in=zbgc_frac_init, bgc_tracer_type_in=bgc_tracer_type)
      call icepack_init_tracer_indices( &
         bio_index_o_in=bio_index_o, bio_index_in=bio_index)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      !-----------------------------------------------------------------
      ! final consistency checks
      !----------------------------------------------------------------- 
      if (nbtrcr > icepack_max_nbtrcr) then
         write (nu_diag,*) subname,' '
         write (nu_diag,*) subname,'nbtrcr > icepack_max_nbtrcr'
         write (nu_diag,*) subname,'nbtrcr, icepack_max_nbtrcr:',nbtrcr, icepack_max_nbtrcr
         call abort_ice (subname//'ERROR: nbtrcr > icepack_max_nbtrcr')
      endif
      if (.NOT. dEdd_algae) nbtrcr_sw = 1

      !-----------------------------------------------------------------
      ! spew
      !-----------------------------------------------------------------
      if (my_task == master_task) then
      if (skl_bgc) then

         write(nu_diag,1020) ' number of bio tracers     = ', nbtrcr
         write(nu_diag,1020) ' number of Isw tracers     = ', nbtrcr_sw
        
      elseif (z_tracers) then
         
         write(nu_diag,1020) ' number of ztracers        = ', nbtrcr
         write(nu_diag,1020) ' number of Isw tracers     = ', nbtrcr_sw
         write(nu_diag,1000) ' initbio_frac              = ', initbio_frac
         write(nu_diag,1000) ' frazil_scav               = ', frazil_scav  

      endif  ! skl_bgc or solve_bgc
      endif  ! master_task

 1000    format (a30,2x,f9.2)  ! a30 to align formatted, unformatted statements
 1020    format (a30,2x,i6)    ! integer

      end subroutine init_zbgc

!=======================================================================

      subroutine init_bgc_trcr(nk,              nt_fbri,       &
                               nt_bgc,          nlt_bgc,       &
                               bgctype,         nt_depend,     &
                               bgc_tracer_type, trcr_depend,   &
                               trcr_base,       n_trcr_strata, &
                               nt_strata,       bio_index)

      integer (kind=int_kind), intent(in) :: &
         nk           , & ! counter
         nt_depend    , & ! tracer dependency index
         nt_bgc       , & ! tracer index
         nlt_bgc      , & ! bio tracer index
         nt_fbri

      integer (kind=int_kind), dimension(:), intent(inout) :: &
         trcr_depend  , & ! tracer dependencies
         n_trcr_strata, & ! number of underlying tracer layers
         bio_index        !

      integer (kind=int_kind), dimension(:,:), intent(inout) :: &
         nt_strata        ! indices of underlying tracer layers

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         trcr_base        ! = 0 or 1 depending on tracer dependency
                          ! argument 2:  (1) aice, (2) vice, (3) vsno

      real (kind=dbl_kind), intent(in) :: &
         bgctype          ! bio tracer transport type (mobile vs stationary)

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         bgc_tracer_type  ! bio tracer transport type array

      ! local variables

      integer (kind=int_kind) :: &
         k         , & ! loop index
         n_strata  , & ! temporary values
         nt_strata1, & ! 
         nt_strata2

      real (kind=dbl_kind) :: &
         trcr_base1, & ! temporary values
         trcr_base2, &
         trcr_base3

      character(len=*), parameter :: subname='(init_bgc_trcr)'

      !--------

      bgc_tracer_type(nlt_bgc) = bgctype
         
      if (nk > 1) then ! include vertical bgc in snow
         do k = nk, nk+1
            trcr_depend  (nt_bgc + k  ) = 2 ! snow volume
            trcr_base    (nt_bgc + k,1) = c0
            trcr_base    (nt_bgc + k,2) = c0
            trcr_base    (nt_bgc + k,3) = c1
            n_trcr_strata(nt_bgc + k  ) = 0
            nt_strata    (nt_bgc + k,1) = 0
            nt_strata    (nt_bgc + k,2) = 0
         enddo

         trcr_base1 = c0      
         trcr_base2 = c1     
         trcr_base3 = c0
         n_strata = 1    
         nt_strata1 = nt_fbri
         nt_strata2 = 0
      else  ! nk = 1
         trcr_base1 = c1
         trcr_base2 = c0
         trcr_base3 = c0
         n_strata = 0
         nt_strata1 = 0
         nt_strata2 = 0
      endif ! nk

      do k = 1, nk     ! in ice
         trcr_depend  (nt_bgc + k - 1  ) = nt_depend
         trcr_base    (nt_bgc + k - 1,1) = trcr_base1
         trcr_base    (nt_bgc + k - 1,2) = trcr_base2
         trcr_base    (nt_bgc + k - 1,3) = trcr_base3
         n_trcr_strata(nt_bgc + k - 1  ) = n_strata
         nt_strata    (nt_bgc + k - 1,1) = nt_strata1
         nt_strata    (nt_bgc + k - 1,2) = nt_strata2
      enddo

      bio_index (nlt_bgc) = nt_bgc

      end subroutine init_bgc_trcr

!=======================================================================

      end module ice_init_column

!=======================================================================
