!=======================================================================

! parameter and variable initializations
!
! authors Elizabeth C. Hunke and William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! 2004 WHL: Block structure added
! 2006 ECH: Added namelist variables, warnings.
!           Replaced old default initial ice conditions with 3.14 version.
!           Converted to free source form (F90).

      module ice_init

      use ice_kinds_mod
      use ice_communicate, only: my_task, master_task, ice_barrier
      use ice_constants, only: c0, c1, c2, c3, c5, c10, c12, &
          p001, p01, p2, p3, p5, p75, p166, cm_to_m
      use ice_exit, only: abort_ice
      use ice_fileunits, only: nu_nml, nu_diag, nml_filename, diag_type, &
          ice_stdout, get_fileunit, release_fileunit, bfbflag, flush_fileunit, &
          ice_IOUnitsMinUnit, ice_IOUnitsMaxUnit
#ifdef CESMCOUPLED
      use ice_fileunits, only: inst_suffix, nu_diag_set
#endif
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_aggregate
      use icepack_intfc, only: icepack_init_trcr
      use icepack_intfc, only: icepack_init_parameters, icepack_write_parameters
      use icepack_intfc, only: icepack_init_tracer_flags
      use icepack_intfc, only: icepack_init_tracer_sizes
      use icepack_intfc, only: icepack_query_tracer_flags
      use icepack_intfc, only: icepack_query_tracer_sizes
      use icepack_intfc, only: icepack_query_tracer_indices
      use icepack_intfc, only: icepack_query_parameters

      implicit none
      private

      character(len=char_len_long), public :: &
         ice_ic      ! method of ice cover initialization
                     ! 'internal' => set from ice_data_ namelist
                     ! 'none'     => no ice
                     ! filename   => read file

      public :: input_data, init_state, set_state_var

!=======================================================================

      contains

!=======================================================================

! Namelist variables, set to default values; may be altered
! at run time
!
! author Elizabeth C. Hunke, LANL

      subroutine input_data

      use ice_broadcast, only: broadcast_scalar, broadcast_array
      use ice_diagnostics, only: &
          diag_file, print_global, print_points, latpnt, lonpnt, &
          debug_model, debug_model_step, debug_model_task, &
          debug_model_i, debug_model_j, debug_model_iblk
      use ice_domain, only: close_boundaries
      use ice_domain_size, only: &
          ncat, nilyr, nslyr, nblyr, nfsd, nfreq, &
          n_iso, n_aero, n_zaero, n_algae, &
          n_doc, n_dic, n_don, n_fed, n_fep, &
          max_nstrm
      use ice_calendar, only: &
          year_init, month_init, day_init, sec_init, &
          istep0, histfreq, histfreq_n, histfreq_base, &
          dumpfreq, dumpfreq_n, diagfreq, dumpfreq_base, &
          npt, dt, ndtd, days_per_year, use_leap_years, &
          write_ic, dump_last, npt_unit
      use ice_arrays_column, only: oceanmixed_ice
      use ice_restart_column, only: &
          restart_age, restart_FY, restart_lvl, &
          restart_pond_lvl, restart_pond_topo, restart_pond_sealvl, &
          restart_aero, restart_fsd, restart_iso, restart_snow
      use ice_restart_shared, only: &
          restart, restart_ext, restart_coszen, use_restart_time, &
          runtype, restart_file, restart_dir, runid, pointer_file, &
          restart_format, restart_rearranger, restart_iotasks, restart_root, &
          restart_stride, restart_deflate, restart_chunksize, restart_mod
      use ice_history_shared, only: &
          history_precision, hist_avg, history_format, history_file, incond_file, &
          history_dir, incond_dir, version_name, history_rearranger, &
          hist_suffix, history_iotasks, history_root, history_stride, &
          history_deflate, history_chunksize, hist_time_axis
      use ice_flux, only: update_ocn_f, cpl_frazil, l_mpond_fresh
      use ice_flux, only: default_season
      use ice_flux_bgc, only: cpl_bgc
      use ice_forcing, only: &
          ycycle,          fyear_init,    debug_forcing, &
          atm_data_type,   atm_data_dir,  precip_units, rotate_wind, &
          atm_data_format, ocn_data_format, atm_data_version, &
          bgc_data_type, &
          ocn_data_type, ocn_data_dir, wave_spec_file,  &
          oceanmixed_file, restore_ocn, trestore, &
          ice_data_type, ice_data_conc, ice_data_dist, &
          snw_filename, &
          snw_tau_fname, snw_kappa_fname, snw_drdt0_fname, &
          snw_rhos_fname, snw_Tgrd_fname, snw_T_fname
      use ice_arrays_column, only: bgc_data_dir, fe_data_type
      use ice_grid, only: &
          grid_file, gridcpl_file, kmt_file, &
          bathymetry_file, use_bathymetry, &
          bathymetry_format, kmt_type, &
          grid_type, grid_format, grid_outfile, &
          grid_ice, grid_ice_thrm, grid_ice_dynu, grid_ice_dynv, &
          grid_ocn, grid_ocn_thrm, grid_ocn_dynu, grid_ocn_dynv, &
          grid_atm, grid_atm_thrm, grid_atm_dynu, grid_atm_dynv, &
          dxrect, dyrect, dxscale, dyscale, scale_dxdy, &
          lonrefrect, latrefrect, save_ghte_ghtn
      use ice_dyn_shared, only: &
          ndte, kdyn, revised_evp, yield_curve, &
          evp_algorithm, visc_method,     &
          seabed_stress, seabed_stress_method,  &
          k1, k2, alphab, threshold_hw, Ktens,  &
          e_yieldcurve, e_plasticpot, coriolis, &
          ssh_stress, kridge, brlx, arlx,       &
          deltaminEVP, deltaminVP, capping,     &
          elasticDamp, dyn_area_min, dyn_mass_min
      use ice_dyn_vp, only: &
          maxits_nonlin, precond, dim_fgmres, dim_pgmres, maxits_fgmres, &
          maxits_pgmres, monitor_nonlin, monitor_fgmres, &
          monitor_pgmres, reltol_nonlin, reltol_fgmres, reltol_pgmres, &
          algo_nonlin, fpfunc_andacc, dim_andacc, reltol_andacc, &
          damping_andacc, start_andacc, use_mean_vrel, ortho_type
      use ice_transport_driver, only: advection, conserv_check
      use ice_restoring, only: restore_ice
      use ice_timers, only: timer_stats
      use ice_memusage, only: memory_stats
      use ice_fileunits, only: goto_nml

#ifdef CESMCOUPLED
      use shr_file_mod, only: shr_file_setIO
#endif

      ! local variables

      integer (kind=int_kind) :: &
         nml_error, & ! namelist i/o error flag
         n            ! loop index

#ifdef CESMCOUPLED
      logical :: exists
#endif

      real (kind=dbl_kind) :: ustar_min, albicev, albicei, albsnowv, albsnowi, &
        ahmax, R_ice, R_pnd, R_snw, dT_mlt, rsnw_mlt, emissivity, hi_min, &
        mu_rdg, hs0, dpscale, rfracmin, rfracmax, pndaspect, apnd_sl, hs1, hp1, &
        a_rapid_mode, Rac_rapid_mode, aspect_rapid_mode, dSdt_slow_mode, &
        phi_c_slow_mode, phi_i_mushy, kalg, atmiter_conv, Pstar, Cstar, &
        sw_frac, sw_dtemp, floediam, hfrazilmin, iceruf, iceruf_ocn, &
        rsnw_fall, rsnw_tmax, rhosnew, rhosmin, rhosmax, Tliquidus_max, &
        windmin, drhosdwind, snwlvlfac, tscale_pnd_drain

      integer (kind=int_kind) :: ktherm, kstrength, krdg_partic, krdg_redist, natmiter, &
        kitd, kcatbound, ktransport

      character (len=char_len) :: shortwave, albedo_type, conduct, fbot_xfer_type, &
        tfrz_option, saltflux_option, frzpnd, atmbndy, wave_spec_type, snwredist, snw_aging_table, &
        congel_freeze, capping_method, snw_ssp_table

      logical (kind=log_kind) :: calc_Tsfc, formdrag, highfreq, calc_strair, wave_spec, &
        sw_redist, calc_dragio, use_smliq_pnd, snwgrain, semi_implicit_Tsfc, vapor_flux_correction

      logical (kind=log_kind) :: tr_iage, tr_FY, tr_lvl, tr_pond
      logical (kind=log_kind) :: tr_iso, tr_aero, tr_fsd, tr_snow
      logical (kind=log_kind) :: tr_pond_lvl, tr_pond_topo, tr_pond_sealvl
      integer (kind=int_kind) :: numin, numax  ! unit number limits
      logical (kind=log_kind) :: lcdf64  ! deprecated, backwards compatibility
      logical (kind=log_kind) :: orca_halogrid !deprecated

      integer (kind=int_kind) :: rplvl, rptopo, rpsealvl
      real (kind=dbl_kind)    :: Cf, ksno, puny, ice_ref_salinity, Tocnfrz

      character (len=char_len) :: abort_list
      character (len=char_len)      :: nml_name ! namelist name
      character (len=char_len_long) :: tmpstr2

      character(len=*), parameter :: subname='(input_data)'

      !-----------------------------------------------------------------
      ! Namelist variables
      !-----------------------------------------------------------------

      namelist /setup_nml/ &
        days_per_year,  use_leap_years, istep0,          npt_unit,      &
        dt,             npt,            ndtd,            numin,         &
        runtype,        runid,          bfbflag,         numax,         &
        ice_ic,         restart,        restart_dir,     restart_file,  &
        restart_ext,    use_restart_time, restart_format, lcdf64,       &
        restart_root,   restart_stride, restart_iotasks, restart_rearranger, &
        restart_deflate, restart_chunksize, restart_mod,                &
        pointer_file,   dumpfreq,       dumpfreq_n,      dump_last,     &
        diagfreq,       diag_type,      diag_file,       history_format,&
        history_root,   history_stride, history_iotasks, history_rearranger, &
        hist_time_axis,                                                 &
        print_global,   print_points,   latpnt,          lonpnt,        &
        debug_forcing,  histfreq,       histfreq_n,      hist_avg,      &
        hist_suffix, history_deflate, history_chunksize,                &
        history_dir,    history_file,   history_precision, cpl_bgc,     &
        histfreq_base,  dumpfreq_base,  timer_stats,     memory_stats,  &
        conserv_check,  debug_model,    debug_model_step,               &
        debug_model_i,  debug_model_j,  debug_model_iblk, debug_model_task, &
        year_init,      month_init,     day_init,        sec_init,      &
        write_ic,       incond_dir,     incond_file,     version_name

      namelist /grid_nml/ &
        grid_format,    grid_type,       grid_file,     kmt_file,       &
        bathymetry_file, use_bathymetry, nfsd,          bathymetry_format, &
        ncat,           nilyr,           nslyr,         nblyr,          &
        kcatbound,      gridcpl_file,    dxrect,        dyrect,         &
        dxscale,        dyscale,         lonrefrect,    latrefrect,     &
        scale_dxdy,     grid_outfile,                                   &
        close_boundaries, orca_halogrid, grid_ice,      kmt_type,       &
        grid_atm,       grid_ocn

      namelist /tracer_nml/                                             &
        tr_iage, restart_age,                                           &
        tr_FY, restart_FY,                                              &
        tr_lvl, restart_lvl,                                            &
        tr_pond_lvl, restart_pond_lvl,                                  &
        tr_pond_sealvl, restart_pond_sealvl,                            &
        tr_pond_topo, restart_pond_topo,                                &
        tr_snow, restart_snow,                                          &
        tr_iso, restart_iso,                                            &
        tr_aero, restart_aero,                                          &
        tr_fsd, restart_fsd,                                            &
        n_iso, n_aero, n_zaero, n_algae,                                &
        n_doc, n_dic, n_don, n_fed, n_fep

      namelist /thermo_nml/ &
        kitd,           ktherm,          conduct,     ksno,             &
        a_rapid_mode,   Rac_rapid_mode,  aspect_rapid_mode,             &
        dSdt_slow_mode, phi_c_slow_mode, phi_i_mushy,                   &
        floediam,       hfrazilmin,      Tliquidus_max,   hi_min,       &
        tscale_pnd_drain


      namelist /dynamics_nml/ &
        kdyn,           ndte,           revised_evp,    yield_curve,    &
        evp_algorithm,  elasticDamp,                                    &
        brlx,           arlx,           ssh_stress,                     &
        advection,      coriolis,       kridge,         ktransport,     &
        kstrength,      krdg_partic,    krdg_redist,    mu_rdg,         &
        e_yieldcurve,   e_plasticpot,   visc_method,                    &
        maxits_nonlin,  precond,        dim_fgmres,                     &
        dim_pgmres,     maxits_fgmres,  maxits_pgmres,  monitor_nonlin, &
        monitor_fgmres, monitor_pgmres, reltol_nonlin,  reltol_fgmres,  &
        reltol_pgmres,  algo_nonlin,    dim_andacc,     reltol_andacc,  &
        damping_andacc, start_andacc,   fpfunc_andacc,  use_mean_vrel,  &
        ortho_type,     seabed_stress,  seabed_stress_method,           &
        k1, k2,         alphab,         threshold_hw,                   &
        deltaminEVP,    deltaminVP,     capping_method,                 &
        Cf,             Pstar,          Cstar,          Ktens,          &
        dyn_area_min,   dyn_mass_min

      namelist /shortwave_nml/ &
        shortwave,      albedo_type,     snw_ssp_table,                 &
        albicev,        albicei,         albsnowv,      albsnowi,       &
        ahmax,          R_ice,           R_pnd,         R_snw,          &
        sw_redist,      sw_frac,         sw_dtemp,                      &
        dT_mlt,         rsnw_mlt,        kalg

      namelist /ponds_nml/ &
        hs0,            dpscale,         frzpnd,                        &
        rfracmin,       rfracmax,        pndaspect,     hs1,            &
        hp1,            apnd_sl

      namelist /snow_nml/ &
        snwredist,      snwgrain,        rsnw_fall,     rsnw_tmax,      &
        rhosnew,        rhosmin,         rhosmax,       snwlvlfac,      &
        windmin,        drhosdwind,      use_smliq_pnd, snw_aging_table,&
        snw_filename,   snw_rhos_fname,  snw_Tgrd_fname,snw_T_fname,    &
        snw_tau_fname,  snw_kappa_fname, snw_drdt0_fname

      namelist /forcing_nml/ &
        formdrag,       atmbndy,         calc_strair,   calc_Tsfc,      &
        highfreq,       natmiter,        atmiter_conv,  calc_dragio,    &
        ustar_min,      emissivity,      iceruf,        iceruf_ocn,     &
        fbot_xfer_type, update_ocn_f,    l_mpond_fresh, tfrz_option,    &
        saltflux_option,ice_ref_salinity,cpl_frazil,    congel_freeze,  &
        oceanmixed_ice, restore_ice,     restore_ocn,   trestore,       &
        precip_units,   default_season,  wave_spec_type,nfreq,          &
        atm_data_type,  ocn_data_type,   bgc_data_type, fe_data_type,   &
        ice_data_type,  ice_data_conc,   ice_data_dist,                 &
        fyear_init,     ycycle,          wave_spec_file,restart_coszen, &
        atm_data_dir,   ocn_data_dir,    bgc_data_dir,                  &
        atm_data_format, ocn_data_format, rotate_wind,                  &
        oceanmixed_file, atm_data_version,semi_implicit_Tsfc,           &
        vapor_flux_correction

      !-----------------------------------------------------------------
      ! default values
      !-----------------------------------------------------------------

      abort_list = ""

      call icepack_query_parameters(puny_out=puny,Tocnfrz_out=Tocnfrz)
! nu_diag not yet defined
!      call icepack_warnings_flush(nu_diag)
!      if (icepack_warnings_aborted()) call abort_ice(error_message=subname//'Icepack Abort0', &
!         file=__FILE__, line=__LINE__)

      days_per_year = 365    ! number of days in a year
      use_leap_years= .false.! if true, use leap years (Feb 29)
      year_init = 0          ! initial year
      month_init = 1         ! initial month
      day_init = 1           ! initial day
      sec_init = 0           ! initial second
      istep0 = 0             ! no. of steps taken in previous integrations,
                             ! real (dumped) or imagined (to set calendar)
#ifndef CESMCOUPLED
      dt = 3600.0_dbl_kind   ! time step, s
#endif
      numin = 11             ! min allowed unit number
      numax = 99             ! max allowed unit number
      npt = 99999            ! total number of time steps (dt)
      npt_unit = '1'         ! units of npt 'y', 'm', 'd', 's', '1'
      diagfreq = 24          ! how often diag output is written
      debug_model  = .false. ! debug output
      debug_model_step = 0   ! debug model after this step number
      debug_model_i = -1     ! debug model local i index
      debug_model_j = -1     ! debug model local j index
      debug_model_iblk = -1  ! debug model local iblk number
      debug_model_task = -1  ! debug model local task number
      print_points = .false. ! if true, print point data
      print_global = .true.  ! if true, print global diagnostic data
      timer_stats = .false.  ! if true, print out detailed timer statistics
      memory_stats = .false. ! if true, print out memory information
      bfbflag = 'off'        ! off = optimized
      diag_type = 'stdout'
      diag_file = 'ice_diag.d'
      histfreq(:) = 'x'
      histfreq(1) = '1'      ! output frequency option for different streams
      histfreq(2) = 'h'      ! output frequency option for different streams
      histfreq(3) = 'd'      ! output frequency option for different streams
      histfreq(4) = 'm'      ! output frequency option for different streams
      histfreq(5) = 'y'      ! output frequency option for different streams
      histfreq_n(:) = 1      ! output frequency
      histfreq_base(:) = 'zero' ! output frequency reference date
      hist_avg(:) = .true.   ! if true, write time-averages (not snapshots)
      hist_suffix(:) = 'x'   ! appended to 'history_file' in filename when not 'x'
      history_format = 'cdf1'! history file format
      history_root = -99     ! history iotasks, root, stride sets pes for pio
      history_stride = -99   ! history iotasks, root, stride sets pes for pio
      history_iotasks = -99  ! history iotasks, root, stride sets pes for pio
      history_rearranger = 'default' ! history rearranger for pio
      hist_time_axis = 'end' ! History file time axis averaging interval position
      history_dir  = './'    ! write to executable dir for default
      history_file = 'iceh'  ! history file name prefix
      history_precision = 4  ! precision of history files
      history_deflate = 0    ! compression level for netcdf4
      history_chunksize(:) = 0 ! chunksize for netcdf4
      write_ic = .false.     ! write out initial condition
      cpl_bgc = .false.      ! couple bgc thru driver
      incond_dir = history_dir ! write to history dir for default
      incond_file = 'iceh_ic'! file prefix
      dumpfreq(:) = 'x'         ! restart frequency option
      dumpfreq_n(:) = 1         ! restart frequency
      dumpfreq_base(:) = 'init' ! restart frequency reference date
      dumpfreq(1) = 'y'         ! restart frequency option
      dumpfreq_n(1) = 1         ! restart frequency
      dump_last = .false.    ! write restart on last time step
      restart_dir  = './'    ! write to executable dir for default
      restart_file = 'iced'  ! restart file name prefix
      restart_ext  = .false. ! if true, read/write ghost cells
      restart_coszen  = .false.   ! if true, read/write coszen
      pointer_file = 'ice.restart_file'
      restart_format = 'cdf1'     ! restart file format
      restart_root = -99     ! restart iotasks, root, stride sets pes for pio
      restart_stride = -99   ! restart iotasks, root, stride sets pes for pio
      restart_iotasks = -99  ! restart iotasks, root, stride sets pes for pio
      restart_rearranger = 'default'  ! restart rearranger for pio
      restart_deflate = 0    ! compression level for netcdf4
      restart_chunksize(:) = 0    ! chunksize for netcdf4
      lcdf64       = .false.      ! 64 bit offset for netCDF
      ice_ic       = 'default'    ! latitude and sst-dependent
      grid_format  = 'bin'        ! grid format
         ! ('bin'=binary or 'pop_nc'=pop netcdf or 'mom_nc'=mom netcdf)
      grid_type    = 'rectangular'! define rectangular grid internally
      grid_file    = 'unknown_grid_file'
      grid_ice     = 'B'          ! underlying grid system
      grid_atm     = 'A'          ! underlying atm forcing/coupling grid
      grid_ocn     = 'A'          ! underlying atm forcing/coupling grid
      gridcpl_file = 'unknown_gridcpl_file'
      grid_outfile = .false.      ! write out one-time grid history file
      orca_halogrid = .false.     ! orca haloed grid - deprecated
      bathymetry_file   = 'unknown_bathymetry_file'
      bathymetry_format = 'default'
      use_bathymetry    = .false.
      kmt_type     = 'file'
      kmt_file     = 'unknown_kmt_file'
      version_name = 'unknown_version_name'
      ncat  = 0          ! number of ice thickness categories
      nfsd  = 1          ! number of floe size categories (1 = default)
      nilyr = 0          ! number of vertical ice layers
      nslyr = 0          ! number of vertical snow layers
      nblyr = 0          ! number of bio layers

      kitd = 1           ! type of itd conversions (0 = delta, 1 = linear)
      kcatbound = 1      ! category boundary formula (0 = old, 1 = new, etc)
      kdyn = 1           ! type of dynamics (-1, 0 = off, 1 = evp, 2 = eap, 3 = vp)
      ndtd = 1           ! dynamic time steps per thermodynamic time step
      ndte = 120         ! subcycles per dynamics timestep:  ndte=dt_dyn/dte
      evp_algorithm = 'standard_2d'  ! EVP kernel (standard_2d=standard cice evp; shared_mem_1d=1d shared memory and no mpi
      elasticDamp = 0.36_dbl_kind    ! coefficient for calculating the parameter E
      save_ghte_ghtn = .false.       ! if true, save global hte and htn (global ext.)
      brlx   = 300.0_dbl_kind ! revised_evp values. Otherwise overwritten in ice_dyn_shared
      arlx   = 300.0_dbl_kind ! revised_evp values. Otherwise overwritten in ice_dyn_shared
      revised_evp = .false.   ! if true, use revised procedure for evp dynamics
      yield_curve = 'ellipse' ! yield curve
      kstrength = 1           ! 1 = Rothrock 75 strength, 0 = Hibler 79
      Pstar = 2.75e4_dbl_kind ! constant in Hibler strength formula (kstrength = 0)
      Cstar = 20._dbl_kind    ! constant in Hibler strength formula (kstrength = 0)
      dyn_area_min = p001     ! minimum ice area concentration to activate dynamics
      dyn_mass_min = p01      ! minimum ice mass to activate dynamics (kg/m^2)
      krdg_partic = 1         ! 1 = new participation, 0 = Thorndike et al 75
      krdg_redist = 1         ! 1 = new redistribution, 0 = Hibler 80
      mu_rdg = 3              ! e-folding scale of ridged ice, krdg_partic=1 (m^0.5)
      Cf = 17.0_dbl_kind      ! ratio of ridging work to PE change in ridging
      ksno = 0.3_dbl_kind     ! snow thermal conductivity
      dxrect = 0.0_dbl_kind   ! user defined grid spacing in cm in x direction
      dyrect = 0.0_dbl_kind   ! user defined grid spacing in cm in y direction
      lonrefrect = -156.50_dbl_kind  ! lower left corner lon for rectgrid
      latrefrect =   71.35_dbl_kind  ! lower left corner lat for rectgrid
      scale_dxdy = .false.    ! apply dxscale, dyscale to rectgrid
      dxscale = 1.0_dbl_kind   ! user defined rectgrid x-grid scale factor (e.g., 1.02)
      dyscale = 1.0_dbl_kind   ! user defined rectgrid y-grid scale factor (e.g., 1.02)
      close_boundaries = .false.   ! true = set land on edges of grid
      seabed_stress= .false.  ! if true, seabed stress for landfast is on
      seabed_stress_method  = 'LKD'! LKD = Lemieux et al 2015, probabilistic = Dupont et al. 2022
      k1 = 7.5_dbl_kind       ! 1st free parameter for landfast parameterization
      k2 = 15.0_dbl_kind      ! 2nd free parameter (N/m^3) for landfast parametrization
      alphab = 20.0_dbl_kind  ! alphab=Cb factor in Lemieux et al 2015
      threshold_hw = 30.0_dbl_kind ! max water depth for grounding
      Ktens = 0.0_dbl_kind    ! T=Ktens*P (tensile strength: see Konig and Holland, 2010)
      e_yieldcurve = 2.0_dbl_kind  ! VP aspect ratio of elliptical yield curve
      e_plasticpot = 2.0_dbl_kind  ! VP aspect ratio of elliptical plastic potential
      visc_method = 'avg_zeta' ! calc viscosities at U point: avg_strength, avg_zeta
      deltaminEVP = 1e-11_dbl_kind ! minimum delta for viscosities (EVP, Hunke 2001)
      deltaminVP  = 2e-9_dbl_kind  ! minimum delta for viscosities (VP, Hibler 1979)
      capping_method  = 'max'  ! method for capping of viscosities (max=Hibler 1979,sum=Kreyscher2000)
      maxits_nonlin = 10       ! max nb of iteration for nonlinear solver
      precond = 'pgmres'       ! preconditioner for fgmres: 'ident' (identity), 'diag' (diagonal),
                               ! 'pgmres' (Jacobi-preconditioned GMRES)
      dim_fgmres = 50          ! size of fgmres Krylov subspace
      dim_pgmres = 5           ! size of pgmres Krylov subspace
      maxits_fgmres = 50       ! max nb of iteration for fgmres
      maxits_pgmres = 5        ! max nb of iteration for pgmres
      monitor_nonlin = .false. ! print nonlinear residual norm
      monitor_fgmres = .false. ! print fgmres residual norm
      monitor_pgmres = .false. ! print pgmres residual norm
      ortho_type = 'mgs'       ! orthogonalization procedure 'cgs' or 'mgs'
      reltol_nonlin = 1e-8_dbl_kind ! nonlinear stopping criterion: reltol_nonlin*res(k=0)
      reltol_fgmres = 1e-1_dbl_kind ! fgmres stopping criterion: reltol_fgmres*res(k)
      reltol_pgmres = 1e-6_dbl_kind ! pgmres stopping criterion: reltol_pgmres*res(k)
      algo_nonlin = 'picard'        ! nonlinear algorithm: 'picard' (Picard iteration), 'anderson' (Anderson acceleration)
      fpfunc_andacc = 1        ! fixed point function for Anderson acceleration:
                               ! 1: g(x) = FMGRES(A(x),b(x)), 2: g(x) = x - A(x)x + b(x)
      dim_andacc = 5           ! size of Anderson minimization matrix (number of saved previous residuals)
      reltol_andacc = 1e-6_dbl_kind  ! relative tolerance for Anderson acceleration
      damping_andacc = 0       ! damping factor for Anderson acceleration
      start_andacc = 0         ! acceleration delay factor (acceleration starts at this iteration)
      use_mean_vrel = .true.   ! use mean of previous 2 iterates to compute vrel
      advection  = 'remap'     ! incremental remapping transport scheme
      conserv_check = .false.  ! tracer conservation check
      shortwave = 'ccsm3'      ! 'ccsm3' or 'dEdd' (delta-Eddington)
      snw_ssp_table = 'test'   ! 'test' or 'snicar' dEdd_snicar_ad table data
      albedo_type = 'ccsm3'    ! 'ccsm3' or 'constant'
      ktherm = 1               ! -1 = OFF, 1 = BL99, 2 = mushy thermo
      conduct = 'bubbly'       ! 'MU71' or 'bubbly' (Pringle et al 2007)
      coriolis = 'latitude'    ! latitude dependent, or 'constant'
      ssh_stress = 'geostrophic'  ! 'geostrophic' or 'coupled'
      kridge   = 1             ! -1 = off, 1 = on
      ktransport = 1           ! -1 = off, 1 = on
      calc_Tsfc = .true.       ! calculate surface temperature
      semi_implicit_Tsfc = .false.  ! surface temperature coupling option based on d(hf)/dTs
      vapor_flux_correction = .false.  ! mass/enthalpy correction for evaporation/sublimation
      update_ocn_f = .false.   ! include fresh water and salt fluxes for frazil
      cpl_frazil = 'fresh_ice_correction' ! type of coupling for frazil ice
      ustar_min = 0.005        ! minimum friction velocity for ocean heat flux (m/s)
      hi_min = p01             ! minimum ice thickness allowed (m)
      iceruf = 0.0005_dbl_kind ! ice surface roughness at atmosphere interface (m)
      iceruf_ocn = 0.03_dbl_kind ! under-ice roughness (m)
      calc_dragio = .false.    ! compute dragio from iceruf_ocn and thickness of first ocean level
      emissivity = 0.985       ! emissivity of snow and ice
      l_mpond_fresh = .false.  ! logical switch for including meltpond freshwater
                               ! flux feedback to ocean model
      fbot_xfer_type = 'constant' ! transfer coefficient type for ocn heat flux
      R_ice     = 0.00_dbl_kind   ! tuning parameter for sea ice
      R_pnd     = 0.00_dbl_kind   ! tuning parameter for ponded sea ice
      R_snw     = 1.50_dbl_kind   ! tuning parameter for snow over sea ice
      dT_mlt    = 1.5_dbl_kind    ! change in temp to give non-melt to melt change
                                  ! in snow grain radius
      rsnw_mlt  = 1500._dbl_kind  ! maximum melting snow grain radius
      kalg      = 0.60_dbl_kind   ! algae absorption coefficient for 0.5 m thick layer
                                  ! 0.5 m path of 75 mg Chl a / m2
      apnd_sl   = 0.27_dbl_kind   ! equilibrium pond fraction in sealvl ponds
      hp1       = 0.01_dbl_kind   ! critical pond lid thickness for topo ponds
      hs0       = 0.03_dbl_kind   ! snow depth for transition to bare sea ice (m)
      hs1       = 0.03_dbl_kind   ! snow depth for transition to bare pond ice (m)
      dpscale   = c1              ! alter e-folding time scale for flushing
      frzpnd    = 'cesm'          ! melt pond refreezing parameterization
      rfracmin  = 0.15_dbl_kind   ! minimum retained fraction of meltwater
      rfracmax  = 0.85_dbl_kind   ! maximum retained fraction of meltwater
      pndaspect = 0.8_dbl_kind    ! ratio of pond depth to area fraction
      tscale_pnd_drain = c10      ! mushy macroscopic drainage timescale (days)
      snwredist = 'none'          ! type of snow redistribution
      snw_aging_table = 'test'    ! snow aging lookup table
      snw_filename    = 'unknown' ! snowtable filename
      snw_tau_fname   = 'unknown' ! snowtable file tau fieldname
      snw_kappa_fname = 'unknown' ! snowtable file kappa fieldname
      snw_drdt0_fname = 'unknown' ! snowtable file drdt0 fieldname
      snw_rhos_fname  = 'unknown' ! snowtable file rhos fieldname
      snw_Tgrd_fname  = 'unknown' ! snowtable file Tgrd fieldname
      snw_T_fname     = 'unknown' ! snowtable file T fieldname
      snwgrain  = .false.         ! snow metamorphosis
      use_smliq_pnd = .false.     ! use liquid in snow for ponds
      rsnw_fall =  100.0_dbl_kind ! radius of new snow (10^-6 m) ! advanced snow physics: 54.526 x 10^-6 m
      rsnw_tmax = 1500.0_dbl_kind ! maximum snow radius (10^-6 m)
      rhosnew   =  100.0_dbl_kind ! new snow density (kg/m^3)
      rhosmin   =  100.0_dbl_kind ! minimum snow density (kg/m^3)
      rhosmax   =  450.0_dbl_kind ! maximum snow density (kg/m^3)
      windmin   =   10.0_dbl_kind ! minimum wind speed to compact snow (m/s)
      drhosdwind=   27.3_dbl_kind ! wind compaction factor for snow (kg s/m^4)
      snwlvlfac =    0.3_dbl_kind ! fractional increase in snow depth for bulk redistribution
      albicev   = 0.78_dbl_kind   ! visible ice albedo for h > ahmax
      albicei   = 0.36_dbl_kind   ! near-ir ice albedo for h > ahmax
      albsnowv  = 0.98_dbl_kind   ! cold snow albedo, visible
      albsnowi  = 0.70_dbl_kind   ! cold snow albedo, near IR
      ahmax     = 0.3_dbl_kind    ! thickness above which ice albedo is constant (m)
      atmbndy   = 'similarity'    ! Atm boundary layer: 'similarity', 'constant' or 'mixed'
      default_season  = 'winter'  ! default forcing data, if data is not read in
      fyear_init = 1900           ! first year of forcing cycle
      ycycle = 1                  ! number of years in forcing cycle
      atm_data_format = 'bin'     ! file format ('bin'=binary or 'nc'=netcdf)
      atm_data_type   = 'default'
      atm_data_dir    = ' '
      atm_data_version = '_undef'  ! date atm_data_file was generated.
      rotate_wind     = .true.    ! rotate wind/stress composants to computational grid orientation
      calc_strair     = .true.    ! calculate wind stress
      formdrag        = .false.   ! calculate form drag
      highfreq        = .false.   ! calculate high frequency RASM coupling
      natmiter        = 5         ! number of iterations for atm boundary layer calcs
      atmiter_conv    = c0        ! ustar convergence criteria
      precip_units    = 'mks'     ! 'mm_per_month' or
                                  ! 'mm_per_sec' = 'mks' = kg/m^2 s
      congel_freeze   = 'two-step'! congelation freezing method
      tfrz_option     = 'mushy'   ! freezing temp formulation
      saltflux_option = 'constant'    ! saltflux calculation
      ice_ref_salinity = 4.0_dbl_kind ! Ice reference salinity for coupling
      oceanmixed_ice  = .false.   ! if true, use internal ocean mixed layer
      wave_spec_type  = 'none'    ! type of wave spectrum forcing
      nfreq           = 25        ! number of wave frequencies
      wave_spec_file  = ' '       ! wave forcing file name
      ocn_data_format = 'bin'     ! file format ('bin'=binary or 'nc'=netcdf)
      bgc_data_type   = 'default'
      fe_data_type    = 'default'
      ice_data_type   = 'default' ! used by some tests to initialize ice state (overall type and mask)
      ice_data_conc   = 'default' ! used by some tests to initialize ice state (concentration)
      ice_data_dist   = 'default' ! used by some tests to initialize ice state (distribution)
      bgc_data_dir    = 'unknown_bgc_data_dir'
      ocn_data_type   = 'default'
      ocn_data_dir    = 'unknown_ocn_data_dir'
      oceanmixed_file = 'unknown_oceanmixed_file' ! ocean forcing data
      restore_ocn     = .false.   ! restore sst if true
      trestore        = 90        ! restoring timescale, days (0 instantaneous)
      restore_ice     = .false.   ! restore ice state on grid edges if true
      restart_mod     = 'none'    ! restart modification option
      debug_forcing   = .false.   ! true writes diagnostics for input forcing

      latpnt(1) =  90._dbl_kind   ! latitude of diagnostic point 1 (deg)
      lonpnt(1) =   0._dbl_kind   ! longitude of point 1 (deg)
      latpnt(2) = -65._dbl_kind   ! latitude of diagnostic point 2 (deg)
      lonpnt(2) = -45._dbl_kind   ! longitude of point 2 (deg)

#ifndef CESMCOUPLED
      runid   = 'unknown'   ! run ID used in CESM and for machine 'bering'
      runtype = 'initial'   ! run type: 'initial', 'continue'
      restart = .false.     ! if true, read ice state from restart file
      use_restart_time = .false.   ! if true, use time info written in file
#endif

      ! extra tracers
      tr_iage      = .false. ! ice age
      restart_age  = .false. ! ice age restart
      tr_FY        = .false. ! ice age
      restart_FY   = .false. ! ice age restart
      tr_lvl       = .false. ! level ice
      restart_lvl  = .false. ! level ice restart
      tr_pond_lvl  = .false. ! level-ice melt ponds
      restart_pond_lvl  = .false. ! melt ponds restart
      tr_pond_sealvl  = .false. ! Sea level melt ponds
      restart_pond_sealvl  = .false. ! Sea level melt ponds restart
      tr_pond_topo = .false. ! explicit melt ponds (topographic)
      restart_pond_topo = .false. ! melt ponds restart
      tr_snow      = .false. ! advanced snow physics
      restart_snow = .false. ! advanced snow physics restart
      tr_iso       = .false. ! isotopes
      restart_iso  = .false. ! isotopes restart
      tr_aero      = .false. ! aerosols
      restart_aero = .false. ! aerosols restart
      tr_fsd       = .false. ! floe size distribution
      restart_fsd  = .false. ! floe size distribution restart

      n_iso = 0
      n_aero = 0
      n_zaero = 0
      n_algae = 0
      n_doc = 0
      n_dic = 0
      n_don = 0
      n_fed = 0
      n_fep = 0

      ! mushy layer gravity drainage physics
      a_rapid_mode      =  0.5e-3_dbl_kind ! channel radius for rapid drainage mode (m)
      Rac_rapid_mode    =    10.0_dbl_kind ! critical Rayleigh number
      aspect_rapid_mode =     1.0_dbl_kind ! aspect ratio (larger is wider)
      dSdt_slow_mode    = -1.5e-7_dbl_kind ! slow mode drainage strength (m s-1 K-1)
      phi_c_slow_mode   =    0.05_dbl_kind ! critical liquid fraction porosity cutoff
      phi_i_mushy       =    0.85_dbl_kind ! liquid fraction of congelation ice
      Tliquidus_max     =    0.00_dbl_kind ! maximum liquidus temperature of mush (C)

      floediam          =   300.0_dbl_kind ! min thickness of new frazil ice (m)
      hfrazilmin        =    0.05_dbl_kind ! effective floe diameter (m)

      ! shortwave redistribution in the thermodynamics
      sw_redist = .false.
      sw_frac   = 0.9_dbl_kind
      sw_dtemp  = 0.02_dbl_kind

      !-----------------------------------------------------------------
      ! read from input file
      !-----------------------------------------------------------------

#ifdef CESMCOUPLED
      nml_filename  = 'ice_in'//trim(inst_suffix)
#endif

      if (my_task == master_task) then

         ! open namelist file
         call get_fileunit(nu_nml)
         open (nu_nml, file=trim(nml_filename), status='old',iostat=nml_error)
         if (nml_error /= 0) then
            call abort_ice(subname//'ERROR: open file '// &
               trim(nml_filename), &
               file=__FILE__, line=__LINE__)
         endif

         ! read setup_nml
         nml_name = 'setup_nml'
         write(nu_diag,*) subname,' Reading ', trim(nml_name)
         ! goto namelist in file
         call goto_nml(nu_nml,trim(nml_name),nml_error)
         if (nml_error /= 0) then
            call abort_ice(subname//'ERROR: searching for '// trim(nml_name), &
               file=__FILE__, line=__LINE__)
         endif

         ! read namelist
         nml_error =  1
         do while (nml_error > 0)
            read(nu_nml, nml=setup_nml,iostat=nml_error)
            ! check if error
            if (nml_error /= 0) then
               ! backspace and re-read erroneous line
               backspace(nu_nml)
               read(nu_nml,fmt='(A)') tmpstr2
               call abort_ice(subname//'ERROR: '//trim(nml_name)//' reading '// &
                    trim(tmpstr2), file=__FILE__, line=__LINE__)
            endif
         end do

         ! read grid_nml
         nml_name = 'grid_nml'
         write(nu_diag,*) subname,' Reading ', trim(nml_name)
         ! goto namelist in file
         call goto_nml(nu_nml,trim(nml_name),nml_error)
         if (nml_error /= 0) then
            call abort_ice(subname//'ERROR: searching for '// trim(nml_name), &
               file=__FILE__, line=__LINE__)
         endif

         ! read namelist
         nml_error =  1
         do while (nml_error > 0)
            read(nu_nml, nml=grid_nml,iostat=nml_error)
            ! check if error
            if (nml_error /= 0) then
               ! backspace and re-read erroneous line
               backspace(nu_nml)
               read(nu_nml,fmt='(A)') tmpstr2
               call abort_ice(subname//'ERROR: ' //trim(nml_name)//' reading '// &
                    trim(tmpstr2), file=__FILE__, line=__LINE__)
            endif
         end do

         ! read tracer_nml
         nml_name = 'tracer_nml'
         write(nu_diag,*) subname,' Reading ', trim(nml_name)
         ! goto namelist in file
         call goto_nml(nu_nml,trim(nml_name),nml_error)
         if (nml_error /= 0) then
            call abort_ice(subname//'ERROR: searching for '// trim(nml_name), &
               file=__FILE__, line=__LINE__)
         endif

         ! read namelist
         nml_error =  1
         do while (nml_error > 0)
            read(nu_nml, nml=tracer_nml,iostat=nml_error)
            ! check if error
            if (nml_error /= 0) then
               ! backspace and re-read erroneous line
               backspace(nu_nml)
               read(nu_nml,fmt='(A)') tmpstr2
               call abort_ice(subname//'ERROR: ' //trim(nml_name)//' reading '// &
                    trim(tmpstr2), file=__FILE__, line=__LINE__)
            endif
         end do

         ! read thermo_nml
         nml_name = 'thermo_nml'
         write(nu_diag,*) subname,' Reading ', trim(nml_name)
         ! goto namelist in file
         call goto_nml(nu_nml,trim(nml_name),nml_error)
         if (nml_error /= 0) then
            call abort_ice(subname//'ERROR: searching for '// trim(nml_name), &
               file=__FILE__, line=__LINE__)
         endif

         ! read namelist
         nml_error =  1
         do while (nml_error > 0)
            read(nu_nml, nml=thermo_nml,iostat=nml_error)
            ! check if error
            if (nml_error /= 0) then
               ! backspace and re-read erroneous line
               backspace(nu_nml)
               read(nu_nml,fmt='(A)') tmpstr2
               call abort_ice(subname//'ERROR: '//trim(nml_name)//' reading '// &
                    trim(tmpstr2), file=__FILE__, line=__LINE__)
            endif
         end do

         ! read dynamics_nml
         nml_name = 'dynamics_nml'
         write(nu_diag,*) subname,' Reading ', trim(nml_name)

         ! goto namelist in file
         call goto_nml(nu_nml,trim(nml_name),nml_error)
         if (nml_error /= 0) then
            call abort_ice(subname//'ERROR: searching for '// trim(nml_name), &
               file=__FILE__, line=__LINE__)
         endif

         ! read namelist
         nml_error =  1
         do while (nml_error > 0)
            read(nu_nml, nml=dynamics_nml,iostat=nml_error)
            ! check if error
            if (nml_error /= 0) then
               ! backspace and re-read erroneous line
               backspace(nu_nml)
               read(nu_nml,fmt='(A)') tmpstr2
               call abort_ice(subname//'ERROR: '//trim(nml_name)//' reading '// &
                    trim(tmpstr2), file=__FILE__, line=__LINE__)
            endif
         end do

         ! read shortwave_nml
         nml_name = 'shortwave_nml'
         write(nu_diag,*) subname,' Reading ', trim(nml_name)

         ! goto namelist in file
         call goto_nml(nu_nml,trim(nml_name),nml_error)
         if (nml_error /= 0) then
            call abort_ice(subname//'ERROR: searching for '// trim(nml_name), &
               file=__FILE__, line=__LINE__)
         endif

         ! read namelist
         nml_error =  1
         do while (nml_error > 0)
            read(nu_nml, nml=shortwave_nml,iostat=nml_error)
            ! check if error
            if (nml_error /= 0) then
               ! backspace and re-read erroneous line
               backspace(nu_nml)
               read(nu_nml,fmt='(A)') tmpstr2
               call abort_ice(subname//'ERROR: '//trim(nml_name)//' reading '//&
                    trim(tmpstr2), file=__FILE__, line=__LINE__)
            endif
         end do

         ! read ponds_nml
         nml_name = 'ponds_nml'
         write(nu_diag,*) subname,' Reading ', trim(nml_name)

         ! goto namelist in file
         call goto_nml(nu_nml,trim(nml_name),nml_error)
         if (nml_error /= 0) then
            call abort_ice(subname//'ERROR: searching for '// trim(nml_name), &
               file=__FILE__, line=__LINE__)
         endif

         ! read namelist
         nml_error =  1
         do while (nml_error > 0)
            read(nu_nml, nml=ponds_nml,iostat=nml_error)
            ! check if error
            if (nml_error /= 0) then
               ! backspace and re-read erroneous line
               backspace(nu_nml)
               read(nu_nml,fmt='(A)') tmpstr2
               call abort_ice(subname//'ERROR: '//trim(nml_name)//' reading '// &
                    trim(tmpstr2), file=__FILE__, line=__LINE__)
            endif
         end do

         ! read snow_nml
         nml_name = 'snow_nml'
         write(nu_diag,*) subname,' Reading ', trim(nml_name)

         ! goto namelist in file
         call goto_nml(nu_nml,trim(nml_name),nml_error)
         if (nml_error /= 0) then
            call abort_ice(subname//'ERROR: searching for '// trim(nml_name), &
               file=__FILE__, line=__LINE__)
         endif

         ! read  namelist
         nml_error =  1
         do while (nml_error > 0)
            read(nu_nml, nml=snow_nml,iostat=nml_error)
            ! check if error
            if (nml_error /= 0) then
               ! backspace and re-read erroneous line
               backspace(nu_nml)
               read(nu_nml,fmt='(A)') tmpstr2
               call abort_ice(subname//'ERROR: '//trim(nml_name)//' reading '// &
                    trim(tmpstr2), file=__FILE__, line=__LINE__)
            endif
         end do

         ! read forcing_nml
         nml_name = 'forcing_nml'
         write(nu_diag,*) subname,' Reading ', trim(nml_name)

         ! goto namelist in file
         call goto_nml(nu_nml,trim(nml_name),nml_error)
         if (nml_error /= 0) then
            call abort_ice(subname//'ERROR: searching for '// trim(nml_name), &
               file=__FILE__, line=__LINE__)
         endif

         ! read namelist
         nml_error =  1
         do while (nml_error > 0)
            read(nu_nml, nml=forcing_nml,iostat=nml_error)
            ! check if error
            if (nml_error /= 0) then
               ! backspace and re-read erroneous line
               backspace(nu_nml)
               read(nu_nml,fmt='(A)') tmpstr2
               call abort_ice(subname//'ERROR: '// trim(nml_name)//' reading '// &
                    trim(tmpstr2), file=__FILE__, line=__LINE__)
            endif
         end do

         ! done reading namelist.
         close(nu_nml)
         call release_fileunit(nu_nml)
      endif

      !-----------------------------------------------------------------
      ! set up diagnostics output and resolve conflicts
      !-----------------------------------------------------------------

#ifdef CESMCOUPLED
      ! Note in CESMCOUPLED mode diag_file is not utilized and
      ! runid and runtype are obtained from the driver, not from the namelist

      if (my_task == master_task) then
         history_file  = trim(runid) // ".cice" // trim(inst_suffix) //".h"
         restart_file  = trim(runid) // ".cice" // trim(inst_suffix) //".r"
         incond_file   = trim(runid) // ".cice" // trim(inst_suffix) //".i"
         ! Note by tcraig - this if test is needed because the nuopc cap sets
         ! nu_diag before this routine is called.  This creates a conflict.
         ! In addition, in the nuopc cap, shr_file_setIO will fail if the
         ! needed namelist is missing (which it is in the CIME nuopc implementation)
         if (.not. nu_diag_set) then
            inquire(file='ice_modelio.nml'//trim(inst_suffix),exist=exists)
            if (exists) then
               call get_fileUnit(nu_diag)
               call shr_file_setIO('ice_modelio.nml'//trim(inst_suffix),nu_diag)
            end if
         endif
      else
         ! each task gets unique ice log filename when if test is true, for debugging
         if (1 == 0) then
            call get_fileUnit(nu_diag)
            write(tmpstr2,'(a,i4.4)') "ice.log.task_",my_task
            open(nu_diag,file=tmpstr2)
         endif
      end if
      if (trim(ice_ic) /= 'default' .and. &
          trim(ice_ic) /= 'none'    .and. &
          trim(ice_ic) /= 'internal') then
         restart = .true.
      end if
#else
      if (trim(diag_type) == 'file') call get_fileunit(nu_diag)
#endif

      !-----------------------------------------------------------------
      ! broadcast namelist settings
      !-----------------------------------------------------------------

      call broadcast_scalar(numin,                master_task)
      call broadcast_scalar(numax,                master_task)
      call broadcast_scalar(days_per_year,        master_task)
      call broadcast_scalar(use_leap_years,       master_task)
      call broadcast_scalar(year_init,            master_task)
      call broadcast_scalar(month_init,           master_task)
      call broadcast_scalar(day_init,             master_task)
      call broadcast_scalar(sec_init,             master_task)
      call broadcast_scalar(istep0,               master_task)
      call broadcast_scalar(dt,                   master_task)
      call broadcast_scalar(npt,                  master_task)
      call broadcast_scalar(npt_unit,             master_task)
      call broadcast_scalar(diagfreq,             master_task)
      call broadcast_scalar(debug_model,          master_task)
      call broadcast_scalar(debug_model_step,     master_task)
      call broadcast_scalar(debug_model_i,        master_task)
      call broadcast_scalar(debug_model_j,        master_task)
      call broadcast_scalar(debug_model_iblk,     master_task)
      call broadcast_scalar(debug_model_task,     master_task)
      call broadcast_scalar(print_points,         master_task)
      call broadcast_scalar(print_global,         master_task)
      call broadcast_scalar(timer_stats,          master_task)
      call broadcast_scalar(memory_stats,         master_task)
      call broadcast_scalar(bfbflag,              master_task)
      call broadcast_scalar(diag_type,            master_task)
      call broadcast_scalar(diag_file,            master_task)
      do n = 1, max_nstrm
         call broadcast_scalar(histfreq(n),       master_task)
         call broadcast_scalar(histfreq_base(n),  master_task)
         call broadcast_scalar(dumpfreq(n),       master_task)
         call broadcast_scalar(dumpfreq_base(n),  master_task)
         call broadcast_scalar(hist_suffix(n),    master_task)
      enddo
      call broadcast_array(hist_avg,              master_task)
      call broadcast_array(histfreq_n,            master_task)
      call broadcast_array(dumpfreq_n,            master_task)
      call broadcast_scalar(history_dir,          master_task)
      call broadcast_scalar(history_file,         master_task)
      call broadcast_scalar(history_precision,    master_task)
      call broadcast_scalar(history_format,       master_task)
      call broadcast_scalar(history_iotasks,      master_task)
      call broadcast_scalar(history_root,         master_task)
      call broadcast_scalar(history_stride,       master_task)
      call broadcast_scalar(history_rearranger,   master_task)
      call broadcast_scalar(hist_time_axis,       master_task)
      call broadcast_scalar(history_deflate,      master_task)
      call broadcast_array(history_chunksize,     master_task)
      call broadcast_scalar(write_ic,             master_task)
      call broadcast_scalar(cpl_bgc,              master_task)
      call broadcast_scalar(incond_dir,           master_task)
      call broadcast_scalar(incond_file,          master_task)
      call broadcast_scalar(version_name,         master_task)
      call broadcast_scalar(dump_last,            master_task)
      call broadcast_scalar(restart_file,         master_task)
      call broadcast_scalar(restart,              master_task)
      call broadcast_scalar(restart_dir,          master_task)
      call broadcast_scalar(restart_ext,          master_task)
      call broadcast_scalar(restart_coszen,       master_task)
      call broadcast_scalar(use_restart_time,     master_task)
      call broadcast_scalar(restart_format,       master_task)
      call broadcast_scalar(restart_iotasks,      master_task)
      call broadcast_scalar(restart_root,         master_task)
      call broadcast_scalar(restart_stride,       master_task)
      call broadcast_scalar(restart_rearranger,   master_task)
      call broadcast_scalar(restart_deflate,      master_task)
      call broadcast_array(restart_chunksize,     master_task)
      call broadcast_scalar(restart_mod,          master_task)
      call broadcast_scalar(lcdf64,               master_task)
      call broadcast_scalar(pointer_file,         master_task)
      call broadcast_scalar(ice_ic,               master_task)
      call broadcast_scalar(grid_format,          master_task)
      call broadcast_scalar(dxrect,               master_task)
      call broadcast_scalar(dyrect,               master_task)
      call broadcast_scalar(scale_dxdy,           master_task)
      call broadcast_scalar(dxscale,              master_task)
      call broadcast_scalar(dyscale,              master_task)
      call broadcast_scalar(lonrefrect,           master_task)
      call broadcast_scalar(latrefrect,           master_task)
      call broadcast_scalar(close_boundaries,     master_task)
      call broadcast_scalar(grid_type,            master_task)
      call broadcast_scalar(grid_ice,             master_task)
      call broadcast_scalar(grid_ocn,             master_task)
      call broadcast_scalar(grid_atm,             master_task)
      call broadcast_scalar(grid_file,            master_task)
      call broadcast_scalar(gridcpl_file,         master_task)
      call broadcast_scalar(grid_outfile,         master_task)
      call broadcast_scalar(orca_halogrid,        master_task)
      call broadcast_scalar(bathymetry_file,      master_task)
      call broadcast_scalar(bathymetry_format,    master_task)
      call broadcast_scalar(use_bathymetry,       master_task)
      call broadcast_scalar(kmt_type,             master_task)
      call broadcast_scalar(kmt_file,             master_task)
      call broadcast_scalar(kitd,                 master_task)
      call broadcast_scalar(kcatbound,            master_task)
      call broadcast_scalar(kdyn,                 master_task)
      call broadcast_scalar(ndtd,                 master_task)
      call broadcast_scalar(ndte,                 master_task)
      call broadcast_scalar(evp_algorithm,        master_task)
      call broadcast_scalar(elasticDamp,          master_task)
      call broadcast_scalar(brlx,                 master_task)
      call broadcast_scalar(arlx,                 master_task)
      call broadcast_scalar(revised_evp,          master_task)
      call broadcast_scalar(yield_curve,          master_task)
      call broadcast_scalar(kstrength,            master_task)
      call broadcast_scalar(Pstar,                master_task)
      call broadcast_scalar(Cstar,                master_task)
      call broadcast_scalar(dyn_area_min,         master_task)
      call broadcast_scalar(dyn_mass_min,         master_task)
      call broadcast_scalar(krdg_partic,          master_task)
      call broadcast_scalar(krdg_redist,          master_task)
      call broadcast_scalar(mu_rdg,               master_task)
      call broadcast_scalar(Cf,                   master_task)
      call broadcast_scalar(ksno,                 master_task)
      call broadcast_scalar(seabed_stress,        master_task)
      call broadcast_scalar(seabed_stress_method, master_task)
      call broadcast_scalar(k1,                   master_task)
      call broadcast_scalar(k2,                   master_task)
      call broadcast_scalar(alphab,               master_task)
      call broadcast_scalar(threshold_hw,         master_task)
      call broadcast_scalar(Ktens,                master_task)
      call broadcast_scalar(e_yieldcurve,         master_task)
      call broadcast_scalar(e_plasticpot,         master_task)
      call broadcast_scalar(visc_method,          master_task)
      call broadcast_scalar(deltaminEVP,          master_task)
      call broadcast_scalar(deltaminVP,           master_task)
      call broadcast_scalar(capping_method,       master_task)
      call broadcast_scalar(advection,            master_task)
      call broadcast_scalar(conserv_check,        master_task)
      call broadcast_scalar(shortwave,            master_task)
      call broadcast_scalar(snw_ssp_table,        master_task)
      call broadcast_scalar(albedo_type,          master_task)
      call broadcast_scalar(ktherm,               master_task)
      call broadcast_scalar(coriolis,             master_task)
      call broadcast_scalar(ssh_stress,           master_task)
      call broadcast_scalar(kridge,               master_task)
      call broadcast_scalar(ktransport,           master_task)
      call broadcast_scalar(maxits_nonlin,        master_task)
      call broadcast_scalar(precond,              master_task)
      call broadcast_scalar(dim_fgmres,           master_task)
      call broadcast_scalar(dim_pgmres,           master_task)
      call broadcast_scalar(maxits_fgmres,        master_task)
      call broadcast_scalar(maxits_pgmres,        master_task)
      call broadcast_scalar(monitor_nonlin,       master_task)
      call broadcast_scalar(monitor_fgmres,       master_task)
      call broadcast_scalar(monitor_pgmres,       master_task)
      call broadcast_scalar(ortho_type,           master_task)
      call broadcast_scalar(reltol_nonlin,        master_task)
      call broadcast_scalar(reltol_fgmres,        master_task)
      call broadcast_scalar(reltol_pgmres,        master_task)
      call broadcast_scalar(algo_nonlin,          master_task)
      call broadcast_scalar(fpfunc_andacc,        master_task)
      call broadcast_scalar(dim_andacc,           master_task)
      call broadcast_scalar(reltol_andacc,        master_task)
      call broadcast_scalar(damping_andacc,       master_task)
      call broadcast_scalar(start_andacc,         master_task)
      call broadcast_scalar(use_mean_vrel,        master_task)
      call broadcast_scalar(conduct,              master_task)
      call broadcast_scalar(R_ice,                master_task)
      call broadcast_scalar(R_pnd,                master_task)
      call broadcast_scalar(R_snw,                master_task)
      call broadcast_scalar(dT_mlt,               master_task)
      call broadcast_scalar(rsnw_mlt,             master_task)
      call broadcast_scalar(kalg,                 master_task)
      call broadcast_scalar(apnd_sl,              master_task)
      call broadcast_scalar(hp1,                  master_task)
      call broadcast_scalar(hs0,                  master_task)
      call broadcast_scalar(hs1,                  master_task)
      call broadcast_scalar(dpscale,              master_task)
      call broadcast_scalar(frzpnd,               master_task)
      call broadcast_scalar(rfracmin,             master_task)
      call broadcast_scalar(rfracmax,             master_task)
      call broadcast_scalar(pndaspect,            master_task)
      call broadcast_scalar(tscale_pnd_drain,     master_task)
      call broadcast_scalar(snwredist,            master_task)
      call broadcast_scalar(snw_aging_table,      master_task)
      call broadcast_scalar(snw_filename,         master_task)
      call broadcast_scalar(snw_tau_fname,        master_task)
      call broadcast_scalar(snw_kappa_fname,      master_task)
      call broadcast_scalar(snw_drdt0_fname,      master_task)
      call broadcast_scalar(snw_rhos_fname,       master_task)
      call broadcast_scalar(snw_Tgrd_fname,       master_task)
      call broadcast_scalar(snw_T_fname,          master_task)
      call broadcast_scalar(snwgrain,             master_task)
      call broadcast_scalar(use_smliq_pnd,        master_task)
      call broadcast_scalar(rsnw_fall,            master_task)
      call broadcast_scalar(rsnw_tmax,            master_task)
      call broadcast_scalar(rhosnew,              master_task)
      call broadcast_scalar(rhosmin,              master_task)
      call broadcast_scalar(rhosmax,              master_task)
      call broadcast_scalar(windmin,              master_task)
      call broadcast_scalar(drhosdwind,           master_task)
      call broadcast_scalar(snwlvlfac,            master_task)
      call broadcast_scalar(albicev,              master_task)
      call broadcast_scalar(albicei,              master_task)
      call broadcast_scalar(albsnowv,             master_task)
      call broadcast_scalar(albsnowi,             master_task)
      call broadcast_scalar(ahmax,                master_task)
      call broadcast_scalar(atmbndy,              master_task)
      call broadcast_scalar(default_season,       master_task)
      call broadcast_scalar(fyear_init,           master_task)
      call broadcast_scalar(ycycle,               master_task)
      call broadcast_scalar(atm_data_format,      master_task)
      call broadcast_scalar(atm_data_type,        master_task)
      call broadcast_scalar(atm_data_dir,         master_task)
      call broadcast_scalar(atm_data_version,     master_task)
      call broadcast_scalar(rotate_wind,          master_task)
      call broadcast_scalar(calc_strair,          master_task)
      call broadcast_scalar(calc_Tsfc,            master_task)
      call broadcast_scalar(semi_implicit_Tsfc,   master_task)
      call broadcast_scalar(vapor_flux_correction,master_task)
      call broadcast_scalar(formdrag,             master_task)
      call broadcast_scalar(highfreq,             master_task)
      call broadcast_scalar(natmiter,             master_task)
      call broadcast_scalar(atmiter_conv,         master_task)
      call broadcast_scalar(update_ocn_f,         master_task)
      call broadcast_scalar(cpl_frazil,           master_task)
      call broadcast_scalar(l_mpond_fresh,        master_task)
      call broadcast_scalar(ustar_min,            master_task)
      call broadcast_scalar(hi_min,               master_task)
      call broadcast_scalar(iceruf,               master_task)
      call broadcast_scalar(iceruf_ocn,           master_task)
      call broadcast_scalar(calc_dragio,          master_task)
      call broadcast_scalar(emissivity,           master_task)
      call broadcast_scalar(fbot_xfer_type,       master_task)
      call broadcast_scalar(precip_units,         master_task)
      call broadcast_scalar(oceanmixed_ice,       master_task)
      call broadcast_scalar(wave_spec_type,       master_task)
      call broadcast_scalar(wave_spec_file,       master_task)
      call broadcast_scalar(nfreq,                master_task)
      call broadcast_scalar(congel_freeze,        master_task)
      call broadcast_scalar(tfrz_option,          master_task)
      call broadcast_scalar(saltflux_option,      master_task)
      call broadcast_scalar(ice_ref_salinity,     master_task)
      call broadcast_scalar(ocn_data_format,      master_task)
      call broadcast_scalar(bgc_data_type,        master_task)
      call broadcast_scalar(fe_data_type,         master_task)
      call broadcast_scalar(ice_data_type,        master_task)
      call broadcast_scalar(ice_data_conc,        master_task)
      call broadcast_scalar(ice_data_dist,        master_task)
      call broadcast_scalar(bgc_data_dir,         master_task)
      call broadcast_scalar(ocn_data_type,        master_task)
      call broadcast_scalar(ocn_data_dir,         master_task)
      call broadcast_scalar(oceanmixed_file,      master_task)
      call broadcast_scalar(restore_ocn,          master_task)
      call broadcast_scalar(trestore,             master_task)
      call broadcast_scalar(restore_ice,          master_task)
      call broadcast_scalar(debug_forcing,        master_task)
      call broadcast_array (latpnt(1:2),          master_task)
      call broadcast_array (lonpnt(1:2),          master_task)
      call broadcast_scalar(runid,                master_task)
      call broadcast_scalar(runtype,              master_task)
      !call broadcast_scalar(nu_diag,              master_task)

      ! tracers
      call broadcast_scalar(tr_iage,              master_task)
      call broadcast_scalar(restart_age,          master_task)
      call broadcast_scalar(tr_FY,                master_task)
      call broadcast_scalar(restart_FY,           master_task)
      call broadcast_scalar(tr_lvl,               master_task)
      call broadcast_scalar(restart_lvl,          master_task)
      call broadcast_scalar(tr_pond_lvl,          master_task)
      call broadcast_scalar(restart_pond_lvl,     master_task)
      call broadcast_scalar(tr_pond_sealvl,       master_task)
      call broadcast_scalar(restart_pond_sealvl,  master_task)
      call broadcast_scalar(tr_pond_topo,         master_task)
      call broadcast_scalar(restart_pond_topo,    master_task)
      call broadcast_scalar(tr_snow,              master_task)
      call broadcast_scalar(restart_snow,         master_task)
      call broadcast_scalar(tr_iso,               master_task)
      call broadcast_scalar(restart_iso,          master_task)
      call broadcast_scalar(tr_aero,              master_task)
      call broadcast_scalar(restart_aero,         master_task)
      call broadcast_scalar(tr_fsd,               master_task)
      call broadcast_scalar(restart_fsd,          master_task)
      call broadcast_scalar(ncat,                 master_task)
      call broadcast_scalar(nfsd,                 master_task)
      call broadcast_scalar(nilyr,                master_task)
      call broadcast_scalar(nslyr,                master_task)
      call broadcast_scalar(nblyr,                master_task)
      call broadcast_scalar(n_iso,                master_task)
      call broadcast_scalar(n_aero,               master_task)
      call broadcast_scalar(n_zaero,              master_task)
      call broadcast_scalar(n_algae,              master_task)
      call broadcast_scalar(n_doc,                master_task)
      call broadcast_scalar(n_dic,                master_task)
      call broadcast_scalar(n_don,                master_task)
      call broadcast_scalar(n_fed,                master_task)
      call broadcast_scalar(n_fep,                master_task)
      call broadcast_scalar(a_rapid_mode,         master_task)
      call broadcast_scalar(floediam,             master_task)
      call broadcast_scalar(hfrazilmin,           master_task)
      call broadcast_scalar(Rac_rapid_mode,       master_task)
      call broadcast_scalar(aspect_rapid_mode,    master_task)
      call broadcast_scalar(dSdt_slow_mode,       master_task)
      call broadcast_scalar(phi_c_slow_mode,      master_task)
      call broadcast_scalar(phi_i_mushy,          master_task)
      call broadcast_scalar(Tliquidus_max,        master_task)
      call broadcast_scalar(sw_redist,            master_task)
      call broadcast_scalar(sw_frac,              master_task)
      call broadcast_scalar(sw_dtemp,             master_task)

      !-----------------------------------------------------------------
      ! update defaults
      !-----------------------------------------------------------------

      if (trim(ice_ic)        == 'default') ice_ic        = 'internal'
      if (trim(ice_data_conc) == 'default') ice_data_conc = 'parabolic'
      if (trim(ice_data_dist) == 'default') ice_data_dist = 'uniform'
      if (trim(ice_data_type) == 'default') ice_data_type = 'latsst'

      ! For backward compatibility
      if (grid_format ==  'nc') grid_format = 'pop_nc'

      !-----------------------------------------------------------------
      ! verify inputs
      !-----------------------------------------------------------------

      if (my_task == master_task) then
         if (trim(diag_type) == 'file') then
            write(ice_stdout,*) 'Diagnostic output will be in file ',diag_file
            open (nu_diag, file=diag_file, status='unknown')
         endif
         write(nu_diag,*) '--------------------------------'
         write(nu_diag,*) '   ',subname
         write(nu_diag,*) '  CICE model diagnostic output  '
         write(nu_diag,*) '--------------------------------'
         write(nu_diag,*) ' '
      endif

      if (trim(runtype) == 'continue') then
         if (my_task == master_task) then
            write(nu_diag,*) subname//'NOTE: runtype=continue, setting restart=.true.'
            if (.not. use_restart_time) &
               write(nu_diag,*) subname//'NOTE: runtype=continue, setting use_restart_time=.true.'
            write(nu_diag,*) ' '
         endif
         restart = .true.
         use_restart_time = .true.
      elseif (trim(runtype) == 'initial') then
         if (ice_ic == 'none' .or. ice_ic == 'internal') then
            if (my_task == master_task) then
               write(nu_diag,*) subname//'NOTE: ice_ic = none or internal, setting restart flags to .false.'
               if (.not. use_restart_time) &
                  write(nu_diag,*) subname//'NOTE: ice_ic = none or internal, setting use_restart_time=.false.'
               write(nu_diag,*) ' '
            endif
            use_restart_time = .false.
            restart = .false.
            restart_iso =  .false.
            restart_aero =  .false.
            restart_fsd =  .false.
            restart_age =  .false.
            restart_fy =  .false.
            restart_lvl =  .false.
            restart_pond_lvl =  .false.
            restart_pond_topo =  .false.
            restart_snow = .false.
! tcraig, OK to leave as true, needed for boxrestore case
!            restart_ext =  .false.
         else
            if (my_task == master_task) then
               write(nu_diag,*) subname//'NOTE: ice_ic /= none or internal, setting restart=.true.'
               write(nu_diag,*) ' '
            endif
            restart = .true.
         endif
      else
         if (my_task == master_task) then
            write(nu_diag,*) subname//' ERROR: runtype unknown = ',trim(runtype)
         endif
         abort_list = trim(abort_list)//":1"
      endif

      if (history_format /= 'cdf1'        .and. &
          history_format /= 'cdf2'        .and. &
          history_format /= 'cdf5'        .and. &
          history_format /= 'hdf5'        .and. &
          history_format /= 'pnetcdf1'    .and. &
          history_format /= 'pnetcdf2'    .and. &
          history_format /= 'pnetcdf5'    .and. &
          history_format /= 'pio_netcdf'  .and. &  ! backwards compatibility
          history_format /= 'pio_pnetcdf' .and. &  ! backwards compatibility
          history_format /= 'binary'      .and. &
          history_format /= 'default')     then    ! backwards compatibility
         if (my_task == master_task) then
            write(nu_diag,*) subname//' ERROR: history_format unknown = ',trim(history_format)
         endif
         abort_list = trim(abort_list)//":50"
      endif

      if (restart_format /= 'cdf1'        .and. &
          restart_format /= 'cdf2'        .and. &
          restart_format /= 'cdf5'        .and. &
          restart_format /= 'hdf5'        .and. &
          restart_format /= 'pnetcdf1'    .and. &
          restart_format /= 'pnetcdf2'    .and. &
          restart_format /= 'pnetcdf5'    .and. &
          restart_format /= 'pio_netcdf'  .and. &  ! backwards compatibility
          restart_format /= 'pio_pnetcdf' .and. &  ! backwards compatibility
          restart_format /= 'binary'      .and. &
          restart_format /= 'default')     then    ! backwards compatibility
         if (my_task == master_task) then
            write(nu_diag,*) subname//' ERROR: restart_format unknown = ',trim(restart_format)
         endif
         abort_list = trim(abort_list)//":51"
      endif

      ! backwards compatibility for history and restart formats, lcdf64

      if (history_format == 'pio_pnetcdf' .or. history_format == 'pio_netcdf') then
         if (my_task == master_task) then
            write(nu_diag,*) subname//' WARNING: history_format='//trim(history_format)// &
                             ' is deprecated, please update namelist settings'
         endif
      endif
      if (restart_format == 'pio_pnetcdf' .or. restart_format == 'pio_netcdf') then
         if (my_task == master_task) then
            write(nu_diag,*) subname//' WARNING: restart_format='//trim(restart_format)// &
                             ' is deprecated, please update namelist settings'
         endif
      endif

      if (lcdf64) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//' WARNING: lcdf64 is deprecated, please update namelist settings'
         endif

         if (history_format == 'default' .or. history_format == 'pio_netcdf') then
            history_format = 'cdf2'
         elseif (history_format == 'pio_pnetcdf') then
            history_format = 'pnetcdf2'
         else
            if (my_task == master_task) then
               write(nu_diag,*) subname//' ERROR: lcdf64 is T and history_format not supported for '//trim(history_format)
            endif
            abort_list = trim(abort_list)//":52"
         endif

         if (restart_format == 'default' .or. restart_format == 'pio_netcdf') then
            restart_format = 'cdf2'
         elseif (restart_format == 'pio_pnetcdf') then
            restart_format = 'pnetcdf2'
         else
            if (my_task == master_task) then
               write(nu_diag,*) subname//' ERROR: lcdf64 is T and restart_format not supported for '//trim(restart_format)
            endif
            abort_list = trim(abort_list)//":53"
         endif
      else
         if (history_format == 'default' .or. history_format == 'pio_netcdf') then
            history_format = 'cdf1'
         elseif (history_format == 'pio_pnetcdf') then
            history_format = 'pnetcdf1'
         endif

         if (restart_format == 'default' .or. restart_format == 'pio_netcdf') then
            restart_format = 'cdf1'
         elseif (restart_format == 'pio_pnetcdf') then
            restart_format = 'pnetcdf1'
         endif
      endif

      if (ktransport <= 0) then
         advection = 'none'
      endif

      if (ktransport > 0 .and. advection /= 'remap' .and. advection /= 'upwind') then
         if (my_task == master_task) write(nu_diag,*) subname//' ERROR: invalid advection=',trim(advection)
         abort_list = trim(abort_list)//":3"
      endif

      if (ncat == 1 .and. kitd == 1) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//' ERROR: kitd incompatability: ncat=1 and kitd=1'
            write(nu_diag,*) subname//' ERROR:   Remapping the ITD is not allowed for ncat=1.'
            write(nu_diag,*) subname//' ERROR:   Use kitd = 0 (delta function ITD) with kcatbound = 0'
            write(nu_diag,*) subname//' ERROR:   or for column configurations use kcatbound = -1'
         endif
         abort_list = trim(abort_list)//":4"
      endif

      if (ncat /= 1 .and. kcatbound == -1) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//' ERROR: ITD required for ncat > 1'
            write(nu_diag,*) subname//' ERROR:   ncat=',ncat,' kcatbound=',kcatbound
            write(nu_diag,*) subname//' ERROR:   Please review user guide'
         endif
         abort_list = trim(abort_list)//":5"
      endif

      if (kdyn == 1 .and. evp_algorithm == 'shared_mem_1d') then
         save_ghte_ghtn = .true.
      endif

      if (kdyn == 2 .and. revised_evp) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//' WARNING: revised_evp = T with EAP dynamics'
            write(nu_diag,*) subname//' WARNING:   revised_evp is ignored'
         endif
         revised_evp = .false.
      endif

      if (kdyn > 3) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//' WARNING: kdyn out of range'
         endif
         abort_list = trim(abort_list)//":33"
      endif

      if (seabed_stress) then
         if (seabed_stress_method /= 'LKD' .and. seabed_stress_method /= 'probabilistic') then
            if (my_task == master_task) then
               write(nu_diag,*) subname//' ERROR: invalid seabed stress method'
               write(nu_diag,*) subname//' ERROR: seabed_stress_method should be LKD or probabilistic'
            endif
            abort_list = trim(abort_list)//":48"
         endif
      endif

      if (close_boundaries) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//' ERROR: close_boundaries deprecated, '// &
              'use ew_boundary_type=closed and/or ns_boundary_type=closed'
         endif
         abort_list = trim(abort_list)//":49"
      endif

      if (grid_ice == 'CD') then
         if (my_task == master_task) then
            write(nu_diag,*) subname//' ERROR: grid_ice = CD not supported yet'
         endif
         abort_list = trim(abort_list)//":47"
      elseif (grid_ice == 'C_override_D') then
         if (my_task == master_task) then
            write(nu_diag,*) subname//' WARNING: using grid_ice = CD, not supported'
         endif
         grid_ice = 'CD'
      endif

      if (grid_ice == 'C' .or. grid_ice == 'CD') then
         if (kdyn > 1 .or. (kdyn == 1 .and. evp_algorithm /= 'standard_2d')) then
            if (my_task == master_task) then
              write(nu_diag,*) subname//' ERROR: grid_ice = C | CD only supported with kdyn=1 and evp_algorithm=standard_2d'
              write(nu_diag,*) subname//' ERROR: kdyn and/or evp_algorithm and grid_ice inconsistency'
            endif
            abort_list = trim(abort_list)//":46"
         endif
         if (visc_method /= 'avg_zeta' .and. visc_method /= 'avg_strength') then
            if (my_task == master_task) then
               write(nu_diag,*) subname//' ERROR: invalid method for viscosities'
               write(nu_diag,*) subname//' ERROR: visc_method should be avg_zeta or avg_strength'
            endif
            abort_list = trim(abort_list)//":44"
         endif
      endif

      capping = -9.99e30
      if (kdyn == 1 .or. kdyn == 3) then
         if (capping_method == 'max') then
            capping = c1
         elseif (capping_method == 'sum') then
            capping = c0
         else
            if (my_task == master_task) then
               write(nu_diag,*) subname//' ERROR: invalid method for capping viscosities'
               write(nu_diag,*) subname//' ERROR: capping_method should be equal to max or sum'
            endif
            abort_list = trim(abort_list)//":45"
         endif
      endif

      rplvl  = 0
      rptopo = 0
      rpsealvl = 0
      if (tr_pond_lvl ) rplvl  = 1
      if (tr_pond_sealvl ) rpsealvl  = 1
      if (tr_pond_topo) rptopo = 1

      tr_pond = .false. ! explicit melt ponds
      if (rplvl + rptopo + rpsealvl > 0) tr_pond = .true.

      if (rplvl + rptopo + rpsealvl > 1) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//' ERROR: Must use only one melt pond scheme'
         endif
         abort_list = trim(abort_list)//":6"
      endif

      if (tr_pond_lvl .and. .not. tr_lvl) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//' ERROR: tr_pond_lvl=T but tr_lvl=F'
         endif
         abort_list = trim(abort_list)//":30"
      endif

! tcraig - this was originally implemented by resetting hs0=0. EH says it might be OK
! to not reset it but extra calculations are done and it might not be bfb.  In our
! testing, we should explicitly set hs0 to 0. when setting tr_pond_lvl=T, and otherwise
! this will abort (safest option until additional testing is done)
      if (tr_pond_lvl .and. abs(hs0) > puny) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//' ERROR: tr_pond_lvl=T and hs0 /= 0'
         endif
         abort_list = trim(abort_list)//":7"
      endif

      if (semi_implicit_Tsfc .and. tr_pond_topo) then
         if (my_task == master_task) then
            write(nu_diag,*)'ERROR: semi_implicit_Tsfc and tr_pond_topo not supported together'
         endif
         abort_list = trim(abort_list)//":57"
      endif

      if (shortwave(1:4) /= 'dEdd' .and. tr_pond .and. calc_tsfc) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//' ERROR: tr_pond=T, calc_tsfc=T, invalid shortwave'
            write(nu_diag,*) subname//' ERROR:   Must use shortwave=dEdd or dEdd_snicar_ad'
         endif
         abort_list = trim(abort_list)//":8"
      endif

      if (snwredist(1:3) == 'ITD' .and. .not. tr_snow) then
         if (my_task == master_task) then
            write (nu_diag,*) 'ERROR: snwredist on but tr_snow=F'
            write (nu_diag,*) 'ERROR: Use tr_snow=T for snow redistribution'
         endif
         abort_list = trim(abort_list)//":37"
      endif
      if (snwredist(1:4) == 'bulk' .and. .not. tr_lvl) then
         if (my_task == master_task) then
            write (nu_diag,*) 'ERROR: snwredist=bulk but tr_lvl=F'
            write (nu_diag,*) 'ERROR: Use tr_lvl=T for snow redistribution'
         endif
         abort_list = trim(abort_list)//":38"
      endif
      if (snwredist(1:6) == 'ITDrdg' .and. .not. tr_lvl) then
         if (my_task == master_task) then
            write (nu_diag,*) 'ERROR: snwredist=ITDrdg but tr_lvl=F'
            write (nu_diag,*) 'ERROR: Use tr_lvl=T for snow redistribution'
         endif
         abort_list = trim(abort_list)//":39"
      endif
      if (use_smliq_pnd .and. .not. snwgrain) then
         if (my_task == master_task) then
            write (nu_diag,*) 'ERROR: use_smliq_pnd = T but'
            write (nu_diag,*) 'ERROR: snow metamorphosis not used'
            write (nu_diag,*) 'ERROR: Use snwgrain=T with smliq for ponds'
         endif
         abort_list = trim(abort_list)//":40"
      endif
      if (use_smliq_pnd .and. .not. tr_snow) then
         if (my_task == master_task) then
            write (nu_diag,*) 'ERROR: use_smliq_pnd = T but'
            write (nu_diag,*) 'ERROR: snow tracers are not active'
            write (nu_diag,*) 'ERROR: Use tr_snow=T with smliq for ponds'
         endif
         abort_list = trim(abort_list)//":41"
      endif
      if (snwgrain .and. .not. tr_snow) then
         if (my_task == master_task) then
            write (nu_diag,*) 'ERROR: snwgrain=T but tr_snow=F'
            write (nu_diag,*) 'ERROR: Use tr_snow=T for snow metamorphosis'
         endif
         abort_list = trim(abort_list)//":42"
      endif
      if (trim(snw_aging_table) /= 'test' .and. &
          trim(snw_aging_table) /= 'snicar' .and. &
          trim(snw_aging_table) /= 'file') then
         if (my_task == master_task) then
            write (nu_diag,*) 'ERROR: unknown snw_aging_table = '//trim(snw_aging_table)
         endif
         abort_list = trim(abort_list)//":43"
      endif

      if (tr_iso .and. n_iso==0) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//' ERROR: isotopes activated but'
            write(nu_diag,*) subname//' ERROR:   not allocated in tracer array.'
            write(nu_diag,*) subname//' ERROR:   if tr_iso, n_iso must be > 0.'
         endif
         abort_list = trim(abort_list)//":31"
      endif

      if (tr_aero .and. n_aero==0) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//' ERROR: aerosols activated but'
            write(nu_diag,*) subname//' ERROR:   not allocated in tracer array.'
            write(nu_diag,*) subname//' ERROR:   if tr_aero, n_aero must be > 0.'
         endif
         abort_list = trim(abort_list)//":9"
      endif

      if (ncat < 1) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//' ERROR: ncat < 1'
         endif
         abort_list = trim(abort_list)//":32"
      endif

      if (nilyr < 1) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//' ERROR: nilyr < 1'
         endif
         abort_list = trim(abort_list)//":2"
      endif

      if (nslyr < 1) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//' ERROR: nslyr < 1'
         endif
         abort_list = trim(abort_list)//":34"
      endif

      if (nblyr < 1) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//' ERROR: nblyr < 1'
            write(nu_diag,*) subname//' ERROR:   not allowed due to history implementation.'
         endif
         abort_list = trim(abort_list)//":35"
      endif

      if (nfsd < 1) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//' ERROR: nfsd < 1'
            write(nu_diag,*) subname//' ERROR:   not allowed due to history implementation.'
         endif
         abort_list = trim(abort_list)//":36"
      endif

      if (shortwave(1:4) /= 'dEdd' .and. tr_aero) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//' ERROR: tr_aero=T, invalid shortwave'
            write(nu_diag,*) subname//' ERROR:   Must use shortwave=dEdd or dEdd_snicar_ad'
         endif
         abort_list = trim(abort_list)//":10"
      endif

      if (shortwave(1:4) /= 'dEdd' .and. snwgrain) then
         if (my_task == master_task) then
            write (nu_diag,*) subname//' ERROR: snow grain radius is activated'
            write (nu_diag,*) subname//' ERROR:   Must use shortwave=dEdd or dEdd_snicar_ad'
         endif
         abort_list = trim(abort_list)//":17"
      endif

      if ((rfracmin < -puny .or. rfracmin > c1+puny) .or. &
          (rfracmax < -puny .or. rfracmax > c1+puny) .or. &
          (rfracmin > rfracmax)) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//' ERROR: rfracmin, rfracmax must be between 0 and 1'
            write(nu_diag,*) subname//' ERROR:   and rfracmax >= rfracmin'
         endif
         abort_list = trim(abort_list)//":11"
      endif
      rfracmin = min(max(rfracmin,c0),c1)
      rfracmax = min(max(rfracmax,c0),c1)

      if (trim(atm_data_type) == 'monthly' .and. calc_strair) then
         if (my_task == master_task) write(nu_diag,*) subname//' ERROR: atm_data_type=monthly and calc_strair=T'
         abort_list = trim(abort_list)//":12"
      endif

      if (ktherm == 2 .and. .not. calc_Tsfc) then
         if (my_task == master_task) write(nu_diag,*) subname//' ERROR: ktherm = 2 and calc_Tsfc=F'
         abort_list = trim(abort_list)//":13"
      endif

! ech: allow inconsistency for testing sensitivities.  It's not recommended for science runs
      if (ktherm == 1 .and. trim(tfrz_option) /= 'linear_salt') then
         if (my_task == master_task) then
            write(nu_diag,*) subname//' WARNING: ktherm = 1 and tfrz_option = ',trim(tfrz_option)
            write(nu_diag,*) subname//' WARNING:   For consistency, set tfrz_option = linear_salt'
         endif
      endif
      if (ktherm == 2 .and. trim(tfrz_option) /= 'mushy') then
         if (my_task == master_task) then
            write(nu_diag,*) subname//' WARNING: ktherm = 2 and tfrz_option = ',trim(tfrz_option)
            write(nu_diag,*) subname//' WARNING:   For consistency, set tfrz_option = mushy'
         endif
      endif
      if (ktherm == 1 .and. trim(saltflux_option) /= 'constant') then
         if (my_task == master_task) then
            write(nu_diag,*) subname//' WARNING: ktherm = 1 and saltflux_option = ',trim(saltflux_option)
            write(nu_diag,*) subname//' WARNING:   For consistency, set saltflux_option = constant'
         endif
      endif
      if (ktherm == 1 .and. .not.sw_redist) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//' WARNING: ktherm = 1 and sw_redist = ',sw_redist
            write(nu_diag,*) subname//' WARNING:   For consistency, set sw_redist = .true.'
         endif
      endif

      if (trim(atmbndy) == 'default') then
         if (my_task == master_task) then
            write(nu_diag,*) subname//' WARNING: atmbndy = default is deprecated'
            write(nu_diag,*) subname//' WARNING:   setting atmbndy = similarity'
         endif
         atmbndy = 'similarity'
      endif

      if (formdrag) then
         if (trim(atmbndy) == 'constant') then
            if (my_task == master_task) write(nu_diag,*) subname//' ERROR: formdrag=T and atmbndy=constant'
            abort_list = trim(abort_list)//":14"
         endif

         if (.not. calc_strair) then
            if (my_task == master_task) write(nu_diag,*) subname//' ERROR: formdrag=T and calc_strair=F'
            abort_list = trim(abort_list)//":15"
         endif

         if (.not. tr_pond) then
            if (my_task == master_task) write(nu_diag,*) subname//' ERROR: formdrag=T and tr_pond=F'
            abort_list = trim(abort_list)//":16"
         endif

         if (.not. tr_lvl) then
            if (my_task == master_task) write(nu_diag,*) subname//' ERROR: formdrag=T and tr_lvl=F'
            abort_list = trim(abort_list)//":18"
         endif
      endif

      if (trim(fbot_xfer_type) == 'Cdn_ocn' .and. .not. formdrag)  then
         if (my_task == master_task) write(nu_diag,*) subname//' ERROR: formdrag=F and fbot_xfer_type=Cdn_ocn'
         abort_list = trim(abort_list)//":19"
      endif

      if (history_precision .ne. 4 .and. history_precision .ne. 8) then
         write (nu_diag,*) subname//' ERROR: bad value for history_precision, allowed values: 4, 8'
         abort_list = trim(abort_list)//":22"
      endif

      do n = 1,max_nstrm
         if (histfreq_base(n) /= 'init' .and. histfreq_base(n) /= 'zero') then
            write (nu_diag,*) subname//' ERROR: bad value for histfreq_base, allowed values: init, zero: '//trim(histfreq_base(n))
            abort_list = trim(abort_list)//":24"
         endif

         if (dumpfreq_base(n) /= 'init' .and. dumpfreq_base(n) /= 'zero') then
            write (nu_diag,*) subname//' ERROR: bad value for dumpfreq_base, allowed values: init, zero: '//trim(dumpfreq_base(n))
            abort_list = trim(abort_list)//":25"
         endif

         if (.not.(scan(dumpfreq(n)(1:1),'ymdhx1YMDHX') == 1 .and. (dumpfreq(n)(2:2) == '1' .or. dumpfreq(n)(2:2) == ' '))) then
            if (my_task == master_task) then
               write(nu_diag,*) subname//' WARNING: unrecognized dumpfreq=', trim(dumpfreq(n))
               write(nu_diag,*) subname//' WARNING:   No restarts files will be written for this stream'
               write(nu_diag,*) subname//' WARNING:   Allowed values : y,m,d,h,x,1 followed by an optional 1'
            endif
            dumpfreq(n) = 'x'
         endif
      enddo

      if (trim(hist_time_axis) /= 'begin' .and. trim(hist_time_axis) /= 'middle' .and. trim(hist_time_axis) /= 'end') then
         write (nu_diag,*) subname//' ERROR: hist_time_axis value not valid = '//trim(hist_time_axis)
         abort_list = trim(abort_list)//":29"
      endif

#ifdef USE_PIO1
      if (history_deflate/=0 .or. restart_deflate/=0 .or. &
          history_chunksize(1)/=0 .or. history_chunksize(2)/=0 .or. &
          restart_chunksize(1)/=0 .or. restart_chunksize(2)/=0) then
         if (my_task == master_task) write (nu_diag,*) subname//' ERROR: _deflate and _chunksize not compatible with PIO1'
         abort_list = trim(abort_list)//":54"
      endif
#else
#ifndef CESMCOUPLED
      ! history_format not used by nuopc driver
      if (history_format/='hdf5' .and. history_deflate/=0) then
         if (my_task == master_task) then
            write (nu_diag,*) subname//' WARNING: history_deflate not compatible with '//history_format
            write (nu_diag,*) subname//' WARNING: netcdf compression only possible with history_type="hdf5" '
         endif
      endif

      if (history_format/='hdf5' .and. (history_chunksize(1)/=0 .or. history_chunksize(2)/=0)) then
         if (my_task == master_task) then
            write (nu_diag,*) subname//' WARNING: history_chunksize not compatible with '//history_format
            write (nu_diag,*) subname//' WARNING: netcdf chunking only possible with history_type="hdf5" '
         endif
      endif

      if (restart_format/='hdf5' .and. restart_deflate/=0) then
         if (my_task == master_task) then
            write (nu_diag,*) subname//' WARNING: restart_deflate not compatible with '//restart_format
            write (nu_diag,*) subname//' WARNING: netcdf compression only possible with restart_type="hdf5" '
         endif
      endif

      if (restart_format/='hdf5' .and. (restart_chunksize(1)/=0 .or. restart_chunksize(2)/=0)) then
         if (my_task == master_task) then
            write (nu_diag,*) subname//' WARNING: restart_chunksize not compatible with '//restart_format
            write (nu_diag,*) subname//' WARNING: netcdf chunking only possible with restart_type="hdf5" '
         endif
      endif
#endif

      if (history_deflate<0 .or. history_deflate>9) then
         if (my_task == master_task) write (nu_diag,*) subname//&
            ' ERROR: history_deflate value not valid. Allowed range: integers from 0 to 9 '
         abort_list = trim(abort_list)//":55"
      endif

      if (restart_deflate<0 .or. restart_deflate>9) then
         if (my_task == master_task) write (nu_diag,*) subname//&
            ' ERROR: restart_deflate value not valid. Allowed range: integers from 0 to 9 '
         abort_list = trim(abort_list)//":56"
      endif
#endif

      ! Implicit solver input validation
      if (kdyn == 3) then
         if (.not. (trim(algo_nonlin) == 'picard' .or. trim(algo_nonlin) == 'anderson')) then
            if (my_task == master_task) then
               write(nu_diag,*) subname//' ERROR: unknown algo_nonlin: '//algo_nonlin
               write(nu_diag,*) subname//' ERROR:   allowed values: ''picard'', ''anderson'''
            endif
            abort_list = trim(abort_list)//":60"
         endif

         if (trim(algo_nonlin) == 'picard') then
            ! Picard solver is implemented in the Anderson solver; reset number of saved residuals to zero
            dim_andacc = 0
         endif

         if (.not. (trim(precond) == 'ident' .or. trim(precond) == 'diag' .or. trim(precond) == 'pgmres')) then
            if (my_task == master_task) then
               write(nu_diag,*) subname//' ERROR: unknown precond: '//precond
               write(nu_diag,*) subname//' ERROR:   allowed values: ''ident'', ''diag'', ''pgmres'''
            endif
            abort_list = trim(abort_list)//":61"
         endif

         if (.not. (trim(ortho_type) == 'cgs' .or. trim(ortho_type) == 'mgs')) then
            if (my_task == master_task) then
               write(nu_diag,*) subname//' ERROR: unknown ortho_type: '//ortho_type
               write(nu_diag,*) subname//' ERROR:   allowed values: ''cgs'', ''mgs'''
            endif
            abort_list = trim(abort_list)//":62"
         endif
      endif

      if (orca_halogrid) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//' ERROR: orca_halogrid has been deprecated'
         endif
         abort_list = trim(abort_list)//":63"
      endif

      if (trim(grid_type) == 'cpom_grid') then
         if (my_task == master_task) then
            write(nu_diag,*) subname//" ERROR: grid_type = 'cpom_grid' has been deprecated"
         endif
         abort_list = trim(abort_list)//":64"
      endif

      ice_IOUnitsMinUnit = numin
      ice_IOUnitsMaxUnit = numax

      call icepack_init_parameters(Cf_in=Cf)
      call icepack_init_parameters(ksno_in=ksno)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname//'Icepack Abort1', &
         file=__FILE__, line=__LINE__)

      wave_spec = .false.
      if (tr_fsd .and. (trim(wave_spec_type) /= 'none')) wave_spec = .true.
      if (tr_fsd .and. (trim(wave_spec_type) == 'none')) then
            if (my_task == master_task) then
               write(nu_diag,*) subname//' WARNING: tr_fsd=T but wave_spec=F - not recommended'
            endif
      end if

      ! compute grid locations for thermo, u and v fields

      grid_ice_thrm = 'T'
      if (grid_ice == 'A') then
         grid_ice_dynu = 'T'
         grid_ice_dynv = 'T'
      elseif (grid_ice == 'B') then
         grid_ice_dynu = 'U'
         grid_ice_dynv = 'U'
      elseif (grid_ice == 'C') then
         grid_ice_dynu = 'E'
         grid_ice_dynv = 'N'
      elseif (grid_ice == 'CD') then
         grid_ice_dynu = 'NE'
         grid_ice_dynv = 'NE'
      else
         if (my_task == master_task) then
            write(nu_diag,*) subname//' ERROR: unknown grid_ice: '//trim(grid_ice)
         endif
         abort_list = trim(abort_list)//":64"
      endif

      grid_atm_thrm = 'T'
      if (grid_atm == 'A') then
         grid_atm_dynu = 'T'
         grid_atm_dynv = 'T'
      elseif (grid_atm == 'B') then
         grid_atm_dynu = 'U'
         grid_atm_dynv = 'U'
      elseif (grid_atm == 'C') then
         grid_atm_dynu = 'E'
         grid_atm_dynv = 'N'
      elseif (grid_atm == 'CD') then
         grid_atm_dynu = 'NE'
         grid_atm_dynv = 'NE'
      else
         if (my_task == master_task) then
            write(nu_diag,*) subname//' ERROR: unknown grid_atm: '//trim(grid_atm)
         endif
         abort_list = trim(abort_list)//":65"
      endif

      grid_ocn_thrm = 'T'
      if (grid_ocn == 'A') then
         grid_ocn_dynu = 'T'
         grid_ocn_dynv = 'T'
      elseif (grid_ocn == 'B') then
         grid_ocn_dynu = 'U'
         grid_ocn_dynv = 'U'
      elseif (grid_ocn == 'C') then
         grid_ocn_dynu = 'E'
         grid_ocn_dynv = 'N'
      elseif (grid_ocn == 'CD') then
         grid_ocn_dynu = 'NE'
         grid_ocn_dynv = 'NE'
      else
         if (my_task == master_task) then
            write(nu_diag,*) subname//' ERROR: unknown grid_ocn: '//trim(grid_ocn)
         endif
         abort_list = trim(abort_list)//":66"
      endif

      !-----------------------------------------------------------------
      ! spew
      !-----------------------------------------------------------------

      if (my_task == master_task) then

         write(nu_diag,*) ' Overview of model configuration with relevant parameters'
         write(nu_diag,*) '========================================================='
         write(nu_diag,*) 'For details, compare namelist output below with the'
         write(nu_diag,*) 'Case Settings section in the model documentation.'
         write(nu_diag,*) ' '
         write(nu_diag,*) ' Calendar'
         write(nu_diag,*) '--------------------------------'
         write(nu_diag,1020) ' days_per_year    = ',days_per_year,' : number of days in a model year'
         if (use_leap_years) then
            tmpstr2 = ' : leap days are included'
         else
            tmpstr2 = ' : leap days are not included'
         endif
         write(nu_diag,1010) ' use_leap_years   = ',use_leap_years,trim(tmpstr2)
         write(nu_diag,1002) ' dt               = ', dt, ' : model time step'

         write(nu_diag,*) ' '
         write(nu_diag,*) ' Grid, Discretization'
         write(nu_diag,*) '--------------------------------'
         write(nu_diag,1030) ' grid_format      = ',trim(grid_format)
         tmpstr2 = ' '
         if (trim(grid_type) == 'rectangular')    tmpstr2 = ' : internally defined, rectangular grid'
         if (trim(grid_type) == 'regional')       tmpstr2 = ' : grid file, regional grid'
         if (trim(grid_type) == 'displaced_pole') tmpstr2 = ' : grid file with rotated north pole'
         if (trim(grid_type) == 'tripole')        tmpstr2 = ' : grid file with northern hemisphere zipper'
         if (trim(grid_type) == 'latlon')         tmpstr2 = ' : cesm latlon domain file'
         write(nu_diag,1030) ' grid_type        = ',trim(grid_type),trim(tmpstr2)
         write(nu_diag,1030) ' grid_ice         = ',trim(grid_ice)
         write(nu_diag,1030) '   grid_ice_thrm  = ',trim(grid_ice_thrm)
         write(nu_diag,1030) '   grid_ice_dynu  = ',trim(grid_ice_dynu)
         write(nu_diag,1030) '   grid_ice_dynv  = ',trim(grid_ice_dynv)
         write(nu_diag,1030) ' grid_atm         = ',trim(grid_atm)
         write(nu_diag,1030) '   grid_atm_thrm  = ',trim(grid_atm_thrm)
         write(nu_diag,1030) '   grid_atm_dynu  = ',trim(grid_atm_dynu)
         write(nu_diag,1030) '   grid_atm_dynv  = ',trim(grid_atm_dynv)
         write(nu_diag,1030) ' grid_ocn         = ',trim(grid_ocn)
         write(nu_diag,1030) '   grid_ocn_thrm  = ',trim(grid_ocn_thrm)
         write(nu_diag,1030) '   grid_ocn_dynu  = ',trim(grid_ocn_dynu)
         write(nu_diag,1030) '   grid_ocn_dynv  = ',trim(grid_ocn_dynv)
         write(nu_diag,1030) ' kmt_type         = ',trim(kmt_type)
         if (trim(grid_type) == 'rectangular') then
            write(nu_diag,1004) 'lon/lat refrect   = ',lonrefrect,latrefrect
            write(nu_diag,1004) 'dx/dy rect (cm)   = ',dxrect,dyrect
            write(nu_diag,1010) 'scale_dxdy        = ',scale_dxdy
            write(nu_diag,1004) 'dx/dy scale       = ',dxscale,dyscale
         else
            if (use_bathymetry) then
               tmpstr2 = ' : bathymetric input data is used'
            else
               tmpstr2 = ' : bathymetric input data is not used'
            endif
            write(nu_diag,1010) ' use_bathymetry   = ', use_bathymetry,trim(tmpstr2)
            write(nu_diag,1030) ' bathymetry_format= ', trim(bathymetry_format)
         endif
         write(nu_diag,1020) ' nilyr            = ', nilyr, ' : number of ice layers (equal thickness)'
         write(nu_diag,1020) ' nslyr            = ', nslyr, ' : number of snow layers (equal thickness)'
         write(nu_diag,1020) ' nblyr            = ', nblyr, ' : number of bio layers (equal thickness)'
         if (shortwave(1:4) == 'dEdd') &
         write(nu_diag,*) 'dEdd interior and sfc scattering layers are used in both ice, snow (unequal)'
         write(nu_diag,1020) ' ncat             = ', ncat,  ' : number of ice categories'
         if (kcatbound == 0) then
            tmpstr2 = ' : original ITD category bounds'
         elseif (kcatbound == 1) then
            tmpstr2 = ' : round-number category bounds'
         elseif (kcatbound == 2) then
            tmpstr2 = ' : WMO standard ITD categories'
         elseif (kcatbound == -1) then
            tmpstr2 = ' : one thickness category'
         else
            tmpstr2 = ' : unknown value'
         endif
         write(nu_diag,1020) ' kcatbound        = ', kcatbound,trim(tmpstr2)
         if (kitd==0) then
            tmpstr2 = ' : delta function ITD approx'
         else
            tmpstr2 = ' : linear remapping ITD approx'
         endif
         write(nu_diag,1020) ' kitd             = ', kitd,trim(tmpstr2)

         if (tr_fsd) then
            tmpstr2 = ' : floe size distribution is enabled'
         else
            tmpstr2 = ' : floe size distribution is disabled'
         endif
         write(nu_diag,1010) ' tr_fsd           = ', tr_fsd,trim(tmpstr2)
         write(nu_diag,1020) ' nfsd             = ', nfsd, ' : number of floe size categories'

         write(nu_diag,*) ' '
         write(nu_diag,*) ' Horizontal Dynamics'
         write(nu_diag,*) '--------------------------------'
         if (kdyn == 1) then
            tmpstr2 = ' : elastic-viscous-plastic dynamics'
         elseif (kdyn == 2) then
            tmpstr2 = ' : elastic-anisotropic-plastic dynamics'
         elseif (kdyn == 3) then
            tmpstr2 = ' : viscous-plastic dynamics'
         elseif (kdyn < 1) then
            tmpstr2 = ' : dynamics disabled'
         else
            tmpstr2 = ' : unknown value'
         endif
         write(nu_diag,1020) ' kdyn             = ', kdyn,trim(tmpstr2)
         write(nu_diag,1003) ' dyn_area_min     = ', dyn_area_min,' : min ice area concentration to activate dynamics'
         write(nu_diag,1003) ' dyn_mass_min     = ', dyn_mass_min,' : min ice mass to activate dynamics (kg/m2)'
         if (kdyn >= 1) then
            if (kdyn == 1 .or. kdyn == 2) then
               if (revised_evp) then
                  tmpstr2 = ' : revised EVP formulation used'
                  write(nu_diag,1002) ' arlx             = ', arlx, ' : stress equation factor alpha'
                  write(nu_diag,1002) ' brlx             = ', brlx, ' : stress equation factor beta'
               else
                  tmpstr2 = ' : revised EVP formulation not used'
               endif
               write(nu_diag,1010) ' revised_evp      = ', revised_evp,trim(tmpstr2)

               if (evp_algorithm == 'standard_2d') then
                  tmpstr2 = ' : standard 2d EVP solver'
               elseif (evp_algorithm == 'shared_mem_1d') then
                  tmpstr2 = ' : vectorized 1d EVP solver'
               else
                  tmpstr2 = ' : unknown value'
               endif
               write(nu_diag,1031) ' evp_algorithm    = ', trim(evp_algorithm),trim(tmpstr2)
               write(nu_diag,1020) ' ndtd             = ', ndtd, ' : number of dynamics/advection/ridging/steps per thermo timestep'
               write(nu_diag,1020) ' ndte             = ', ndte, ' : number of EVP or EAP subcycles'
            endif

            if (kdyn == 1 .or. kdyn == 3) then
               write(nu_diag,1030) ' yield_curve      = ', trim(yield_curve), ' : yield curve'
               if (trim(yield_curve) == 'ellipse') &
                  write(nu_diag,1002) ' e_yieldcurve     = ', e_yieldcurve, ' : aspect ratio of yield curve'
                  write(nu_diag,1002) ' e_plasticpot     = ', e_plasticpot, ' : aspect ratio of plastic potential'
            endif

            if (kdyn == 1) then
               write(nu_diag,1003) ' deltamin     = ', deltaminEVP, ' : minimum delta for viscosities'
               write(nu_diag,1030) ' capping_meth = ', trim(capping_method), ' : capping method for viscosities'
            elseif (kdyn == 3) then
               write(nu_diag,1003) ' deltamin     = ', deltaminVP, ' : minimum delta for viscosities'
               write(nu_diag,1030) ' capping_meth = ', trim(capping_method), ' : capping method for viscosities'
            endif
            !write(nu_diag,1002) ' capping      = ', capping, ' : capping value for viscosities'

            write(nu_diag,1002) ' elasticDamp  = ', elasticDamp, ' : coefficient for calculating the parameter E'

            if (trim(coriolis) == 'latitude') then
               tmpstr2 = ' : latitude-dependent Coriolis parameter'
            elseif (trim(coriolis) == 'contant') then
               tmpstr2 = ' : = 1.46e-4/s'
            elseif (trim(coriolis) == 'zero') then
               tmpstr2 = ' : = 0.0'
            else
               tmpstr2 = ': unknown value'
            endif
            write(nu_diag,1030) ' coriolis         = ',trim(coriolis),trim(tmpstr2)

            if (trim(ssh_stress) == 'geostrophic') then
               tmpstr2 = ' : from ocean velocity'
            elseif (trim(ssh_stress) == 'coupled') then
               tmpstr2 = ' : from coupled sea surface height gradients'
            else
               tmpstr2 = ' : unknown value'
            endif
            write(nu_diag,1030) ' ssh_stress       = ',trim(ssh_stress),trim(tmpstr2)

            if (trim(advection) == 'remap') then
               tmpstr2 = ' : linear remapping advection'
            elseif (trim(advection) == 'upwind') then
               tmpstr2 = ' : donor cell (upwind) advection'
            elseif (trim(advection) == 'none') then
               tmpstr2 = ' : advection disabled by ktransport namelist'
            else
               tmpstr2 = ' : unknown value'
            endif
            write(nu_diag,1030) ' advection        = ', trim(advection),trim(tmpstr2)

            if (seabed_stress) then
               tmpstr2 = ' : use seabed stress parameterization for landfast ice'
            else
               tmpstr2 = ' : no seabed stress parameterization'
            endif
            write(nu_diag,1010) ' seabed_stress    = ', seabed_stress,trim(tmpstr2)
            if (seabed_stress) then
               write(nu_diag,1030) ' seabed method    = ',trim(seabed_stress_method)
               if (seabed_stress_method == 'LKD') then
                  write(nu_diag,1002) ' k1               = ', k1, ' : free parameter for landfast ice'
                  write(nu_diag,1002) ' k2               = ', k2, ' : free parameter for landfast ice'
                  write(nu_diag,1002) ' alphab           = ', alphab, ' : factor for landfast ice'
                  write(nu_diag,1002) ' threshold_hw     = ', threshold_hw, ' : max water depth for grounding ice'
               elseif (seabed_stress_method == 'probabilistic') then
                  write(nu_diag,1002) ' alphab           = ', alphab, ' : factor for landfast ice'
               endif
            endif
            if (grid_ice == 'C' .or. grid_ice == 'CD') then
               write(nu_diag,1030) ' visc_method= ', trim(visc_method),' : viscosities method (U point)'
            endif

            write(nu_diag,1002) ' Ktens            = ', Ktens, ' : tensile strength factor'

            if (kdyn == 3) then
               write(nu_diag,1020) ' maxits_nonlin    = ', maxits_nonlin,' : max nb of iteration for nonlinear solver'
               write(nu_diag,1030) ' precond          = ', trim(precond),' : preconditioner for FGMRES'
               write(nu_diag,1020) ' dim_fgmres       = ', dim_fgmres,' : size of FGMRES Krylov subspace'
               write(nu_diag,1020) ' dim_pgmres       = ', dim_pgmres,' : size of PGMRES Krylov subspace'
               write(nu_diag,1020) ' maxits_fgmres    = ', maxits_fgmres,' : max nb of iteration for FGMRES'
               write(nu_diag,1020) ' maxits_pgmres    = ', maxits_pgmres,' : max nb of iteration for PGMRES'
               write(nu_diag,1010) ' monitor_nonlin   = ', monitor_nonlin,' : print nonlinear residual norm'
               write(nu_diag,1010) ' monitor_fgmres   = ', monitor_fgmres,' : print FGMRES residual norm'
               write(nu_diag,1010) ' monitor_pgmres   = ', monitor_pgmres,' : print PGMRES residual norm'
               write(nu_diag,1030) ' ortho_type       = ', trim(ortho_type),' : type of orthogonalization for FGMRES'
               write(nu_diag,1009) ' reltol_nonlin    = ', reltol_nonlin,' : nonlinear stopping criterion'
               write(nu_diag,1009) ' reltol_fgmres    = ', reltol_fgmres,' : FGMRES stopping criterion'
               write(nu_diag,1009) ' reltol_pgmres    = ', reltol_pgmres,' : PGMRES stopping criterion'
               write(nu_diag,1030) ' algo_nonlin      = ', trim(algo_nonlin),' : nonlinear algorithm'
               write(nu_diag,1010) ' use_mean_vrel    = ', use_mean_vrel,' : use mean of previous 2 iterates to compute vrel'
               if (algo_nonlin == 'anderson') then
                  write(nu_diag,1020) ' fpfunc_andacc    = ', fpfunc_andacc,' : fixed point function for Anderson acceleration'
                  write(nu_diag,1020) ' dim_andacc       = ', dim_andacc,' : size of Anderson minimization matrix'
                  write(nu_diag,1009) ' reltol_andacc    = ', reltol_andacc,' : relative tolerance for Anderson acceleration'
                  write(nu_diag,1000) ' damping_andacc   = ', damping_andacc,' : damping factor for Anderson acceleration'
                  write(nu_diag,1020) ' start_andacc     = ', start_andacc,' : nonlinear iteration at which acceleration starts'
               endif
            endif

         endif ! kdyn enabled

         write(nu_diag,*) ' '
         write(nu_diag,*) ' Mechanical Deformation (Ridging) and Ice Strength'
         write(nu_diag,*) '--------------------------------------------------'
         if (kridge == 1) then
            tmpstr2 = ' : ridging enabled'
         else
            tmpstr2 = ' : ridging disabled'
         endif
         write(nu_diag,1010) ' tr_lvl           = ', tr_lvl,' : ridging related tracers'
         write(nu_diag,1020) ' kridge           = ', kridge,trim(tmpstr2)
         if (kridge == 1) then
            if (krdg_partic == 1) then
               tmpstr2 = ' : new participation function'
            else
               tmpstr2 = ' : old participation function'
            endif
            write(nu_diag,1020) ' krdg_partic      = ', krdg_partic,trim(tmpstr2)
            if (krdg_partic == 1) &
            write(nu_diag,1002) ' mu_rdg           = ', mu_rdg,' : e-folding scale of ridged ice'
            if (krdg_redist == 1) then
               tmpstr2 = ' : new redistribution function'
            else
               tmpstr2 = ' : old redistribution function'
            endif
            write(nu_diag,1020) ' krdg_redist      = ', krdg_redist,trim(tmpstr2)
         endif

         if (kstrength == 0) then
            tmpstr2 = ' : Hibler (1979)'
         elseif (kstrength == 1) then
            tmpstr2 = ' : Rothrock (1975)'
         else
            tmpstr2 = ' : unknown value'
         endif
         write(nu_diag,1020) ' kstrength        = ', kstrength,trim(tmpstr2)
         if (kstrength == 0) then
            write(nu_diag,1009) ' Pstar            = ', Pstar, ' : P* strength factor'
            write(nu_diag,1002) ' Cstar            = ', Cstar, ' : C* strength exponent factor'
         elseif (kstrength == 1) then
            write(nu_diag,1002) ' Cf               = ', Cf, ' : ratio of ridging work to PE change'
         endif

         write(nu_diag,*) ' '
         write(nu_diag,*) ' Thermodynamics'
         write(nu_diag,*) '--------------------------------'

         if (ktherm == 1) then
            tmpstr2 = ' : Bitz and Lipscomb 1999 thermo'
         elseif (ktherm == 2) then
            tmpstr2 = ' : mushy-layer thermo'
         elseif (ktherm < 0) then
            tmpstr2 = ' : Thermodynamics disabled'
         else
            tmpstr2 = ' : unknown value'
         endif
         write(nu_diag,1020) ' ktherm           = ', ktherm,trim(tmpstr2)
         if (ktherm >= 0) then
            write(nu_diag,1002) ' dt               = ', dt, ' : thermodynamic time step'
            write(nu_diag,1002) ' ksno             = ', ksno,' : snow thermal conductivity'
            if (ktherm == 1) &
            write(nu_diag,1030) ' conduct          = ', trim(conduct),' : ice thermal conductivity'
            if (ktherm == 2) then
               write(nu_diag,1002) ' a_rapid_mode     = ', a_rapid_mode,' : brine channel diameter'
               write(nu_diag,1002) ' Rac_rapid_mode   = ', Rac_rapid_mode,' : critical Rayleigh number'
               write(nu_diag,1002) ' aspect_rapid_mode= ', aspect_rapid_mode,' : brine convection aspect ratio'
               write(nu_diag,1009) ' dSdt_slow_mode   = ', dSdt_slow_mode,' : drainage strength parameter'
               write(nu_diag,1002) ' phi_c_slow_mode  = ', phi_c_slow_mode,' : critical liquid fraction'
               write(nu_diag,1002) ' phi_i_mushy      = ', phi_i_mushy,' : solid fraction at lower boundary'
               write(nu_diag,1002) ' Tliquidus_max    = ', Tliquidus_max,' : max mush liquidus temperature'
            endif
            write(nu_diag,1002) ' hfrazilmin       = ', hfrazilmin,' : minimum new frazil ice thickness'

            write(nu_diag,*) ' '
            write(nu_diag,*) ' Radiation'
            write(nu_diag,*) '--------------------------------'
            if (trim(shortwave) == 'dEdd') then
               tmpstr2 = ' : delta-Eddington multiple-scattering method'
            elseif (trim(shortwave) == 'dEdd_snicar_ad') then
               tmpstr2 = ' : delta-Eddington multiple-scattering method with SNICAR AD'
            elseif (trim(shortwave) == 'ccsm3') then
               tmpstr2 = ' : NCAR CCSM3 distribution method'
            else
               tmpstr2 = ' : unknown value'
            endif
            write(nu_diag,1030) ' shortwave        = ', trim(shortwave),trim(tmpstr2)
            if (shortwave(1:4) == 'dEdd') then
               write(nu_diag,1002) ' R_ice            = ', R_ice,' : tuning parameter for sea ice albedo'
               write(nu_diag,1002) ' R_pnd            = ', R_pnd,' : tuning parameter for ponded sea ice albedo'
               write(nu_diag,1002) ' R_snw            = ', R_snw,' : tuning parameter for snow broadband albedo'
               write(nu_diag,1002) ' dT_mlt           = ', dT_mlt,' : change in temperature per change in snow grain radius'
               write(nu_diag,1002) ' rsnw_mlt         = ', rsnw_mlt,' : maximum melting snow grain radius'
               write(nu_diag,1002) ' kalg             = ', kalg,' : absorption coefficient for algae'
             if (trim(shortwave) == 'dEdd_snicar_ad') then
               write(nu_diag,1030) ' snw_ssp_table    = ', trim(snw_ssp_table)
             endif
            else
               if (trim(albedo_type) == 'ccsm3') then
                  tmpstr2 = ' : NCAR CCSM3 albedos'
               elseif (trim(albedo_type) == 'constant') then
                  tmpstr2 = ' : four constant albedos'
               else
                  tmpstr2 = ' : unknown value'
                  abort_list = trim(abort_list)//":23"
               endif
               write(nu_diag,1030) ' albedo_type     = ', trim(albedo_type),trim(tmpstr2)
               if (trim(albedo_type) == 'ccsm3') then
                  write(nu_diag,1002) ' albicev          = ', albicev,' : visible  ice albedo for thicker ice'
                  write(nu_diag,1002) ' albicei          = ', albicei,' : near infrared ice albedo for thicker ice'
                  write(nu_diag,1002) ' albsnowv         = ', albsnowv,' : visible, cold snow albedo'
                  write(nu_diag,1002) ' albsnowi         = ', albsnowi,' : near infrared, cold snow albedo'
                  write(nu_diag,1002) ' ahmax            = ', ahmax,' : albedo is constant above this thickness'
                  write(nu_diag,1002) ' ahmax            = ', ahmax,' : albedo is constant above this thickness'
               endif
            endif
            write(nu_diag,1000) ' emissivity       = ', emissivity,' : emissivity of snow and ice'
            write(nu_diag,1010) ' sw_redist        = ', sw_redist,' : redistribute internal shortwave to surface'
            if (sw_redist) then
               write(nu_diag,1002) ' sw_frac          = ', sw_frac,' : fraction redistributed'
               write(nu_diag,1002) ' sw_dtemp         = ', sw_dtemp,' : temperature difference from freezing to redistribute'
            endif
         endif

         write(nu_diag,*) ' '
         write(nu_diag,*) ' Atmospheric Forcing / Coupling'
         write(nu_diag,*) '--------------------------------'
         write(nu_diag,1010) ' calc_Tsfc        = ', calc_Tsfc,' : calculate surface temperature as part of thermo'
         write(nu_diag,1010) ' calc_strair      = ', calc_strair,' : calculate wind stress and speed'
         write(nu_diag,1010) ' rotate_wind      = ', rotate_wind,' : rotate wind/stress to computational grid'
         write(nu_diag,1010) ' formdrag         = ', formdrag,' : use form drag parameterization'
         write(nu_diag,1000) ' iceruf           = ', iceruf, ' : ice surface roughness at atmosphere interface (m)'
         write(nu_diag,1010) ' semi_implicit_Tsfc    = ', semi_implicit_Tsfc,' : surface temperature coupling option based on d(hf)/dTs'
         write(nu_diag,1010) ' vapor_flux_correction = ', vapor_flux_correction,' : mass/enthalpy correction for evaporation/sublimation'
         if (trim(atmbndy) == 'constant') then
            tmpstr2 = ' : constant-based boundary layer'
         elseif (trim(atmbndy) == 'similarity' .or. &
                 trim(atmbndy) == 'mixed') then
            write(nu_diag,1010) ' highfreq         = ', highfreq,' : high-frequency atmospheric coupling'
            write(nu_diag,1020) ' natmiter         = ', natmiter,' : number of atmo boundary layer iterations'
            write(nu_diag,1002) ' atmiter_conv     = ', atmiter_conv,' : convergence criterion for ustar'
            if (trim(atmbndy) == 'similarity') then
               tmpstr2 = ' : stability-based boundary layer'
            else
               tmpstr2 = ' : stability-based boundary layer for wind stress, constant-based for sensible+latent heat fluxes'
            endif
         else
            tmpstr2 = ' : unknown value'
         endif
         write(nu_diag,1030) ' atmbndy          = ', trim(atmbndy),trim(tmpstr2)

         write(nu_diag,*) ' '
         write(nu_diag,*) ' Oceanic Forcing / Coupling'
         write(nu_diag,*) '--------------------------------'
         if (oceanmixed_ice) then
            tmpstr2 = ' : ocean mixed layer calculation (SST) enabled'
         else
            tmpstr2 = ' : ocean mixed layer calculation (SST) disabled'
         endif
         write(nu_diag,1010) ' oceanmixed_ice   = ', oceanmixed_ice,trim(tmpstr2)
         if (oceanmixed_ice) then
            write(nu_diag,*) '     WARNING: ocean mixed layer ON'
            write(nu_diag,*) '     WARNING: will impact ocean forcing interaction'
            write(nu_diag,*) '     WARNING: coupled forcing will be modified by mixed layer routine'
         endif
         write(nu_diag,1030) ' saltflux_option  = ', trim(saltflux_option)
         if (trim(saltflux_option) == 'constant') then
            write(nu_diag,1002) ' ice_ref_salinity = ',ice_ref_salinity
         endif
         if (trim(tfrz_option) == 'constant') then
            tmpstr2 = ' : constant ocean freezing temperature (Tocnfrz)'
         elseif (trim(tfrz_option) == 'minus1p8') then
            tmpstr2 = ' : constant ocean freezing temperature (-1.8C) (to be deprecated)'
         elseif (trim(tfrz_option) == 'linear_salt') then
            tmpstr2 = ' : linear function of salinity (use with ktherm=1)'
         elseif (trim(tfrz_option) == 'mushy') then
            tmpstr2 = ' : Assur (1958) as in mushy-layer thermo (ktherm=2)'
         else
            tmpstr2 = ' : unknown value'
         endif
         write(nu_diag,1030) ' tfrz_option      = ', trim(tfrz_option),trim(tmpstr2)
         if (trim(tfrz_option) == 'constant') then
            write(nu_diag,1002) ' Tocnfrz          = ', Tocnfrz
         endif
         write(nu_diag,1030) ' congel_freeze    = ', trim(congel_freeze)
         if (update_ocn_f) then
            tmpstr2 = ' : frazil water/salt fluxes included in ocean fluxes'
         else
            tmpstr2 = ' : frazil water/salt fluxes not included in ocean fluxes'
         endif
         write(nu_diag,1010) ' update_ocn_f     = ', update_ocn_f,trim(tmpstr2)
         write(nu_diag,1030) ' cpl_frazil       = ', trim(cpl_frazil)
         if (l_mpond_fresh .and. tr_pond_topo) then
            tmpstr2 = ' : retain (topo) pond water until ponds drain'
         else
            tmpstr2 = ' : pond water not retained on ice (virtual only)'
         endif
         write(nu_diag,1010) ' l_mpond_fresh    = ', l_mpond_fresh,trim(tmpstr2)
         if (trim(fbot_xfer_type) == 'constant') then
            tmpstr2 = ' : ocean heat transfer coefficient is constant'
         elseif (trim(fbot_xfer_type) == 'Cdn_ocn') then
            tmpstr2 = ' : variable ocean heat transfer coefficient'  ! only used with form_drag=T?
         else
            tmpstr2 = ' : unknown value'
         endif
         write(nu_diag,1030) ' fbot_xfer_type   = ', trim(fbot_xfer_type),trim(tmpstr2)
         write(nu_diag,1000) ' ustar_min        = ', ustar_min,' : minimum value of ocean friction velocity'
         write(nu_diag,1000) ' hi_min           = ', hi_min,' : minimum ice thickness allowed (m)'
         if (calc_dragio) then
            tmpstr2 = ' : dragio computed from iceruf_ocn'
         else
            tmpstr2 = ' : dragio hard-coded'
         endif
         write(nu_diag,1010) ' calc_dragio   = ', calc_dragio,trim(tmpstr2)
         if (calc_dragio) then
            write(nu_diag,1002) ' iceruf_ocn       = ', iceruf_ocn,' : under-ice roughness length'
         endif

         if (tr_fsd) then
            write(nu_diag,1002) ' floediam         = ', floediam, ' constant floe diameter'
            if (wave_spec) then
               tmpstr2 = ' : use wave spectrum for floe size distribution'
            else
               tmpstr2 = 'WARNING : floe size distribution does not use wave spectrum'
            endif
            write(nu_diag,1010) ' wave_spec          = ', wave_spec,trim(tmpstr2)
            if (wave_spec) then
               if (trim(wave_spec_type) == 'none') then
                  tmpstr2 = ' : no wave data provided, no wave-ice interactions'
               elseif (trim(wave_spec_type) == 'profile') then
                  tmpstr2 = ' : use fixed dummy wave spectrum for testing, sea surface height generated '// &
                            'using constant phase (1 iteration of wave fracture)'
               elseif (trim(wave_spec_type) == 'constant') then
                  tmpstr2 = ' : wave spectrum data file provided, sea surface height generated '// &
                            'using constant phase (1 iteration of wave fracture)'
               elseif (trim(wave_spec_type) == 'random') then
                  tmpstr2 = ' : wave spectrum data file provided, sea surface height generated using '// &
                            'random number (multiple iterations of wave fracture to convergence)'
               else
                  tmpstr2 = ' : unknown value'
               endif
               write(nu_diag,1030) ' wave_spec_type   = ', trim(wave_spec_type),trim(tmpstr2)
            endif
            write(nu_diag,1020) ' nfreq            = ', nfreq,' : number of wave spectral forcing frequencies'
         endif

         write(nu_diag,*) ' '
         write(nu_diag,*) ' Age related tracers'
         write(nu_diag,*) '--------------------------------'
         write(nu_diag,1010) ' tr_iage          = ', tr_iage,' : chronological ice age'
         write(nu_diag,1010) ' tr_FY            = ', tr_FY,' : first-year ice area'

         write(nu_diag,*) ' '
         write(nu_diag,*) ' Melt ponds'
         write(nu_diag,*) '--------------------------------'
         if (tr_pond_lvl) then
            write(nu_diag,1010) ' tr_pond_lvl      = ', tr_pond_lvl,' : level-ice pond formulation'
            write(nu_diag,1002) ' pndaspect        = ', pndaspect,' : ratio of pond depth to area fraction'
            write(nu_diag,1000) ' dpscale          = ', dpscale,' : time scale for flushing in permeable ice'
            if (trim(frzpnd) == 'hlid') then
               tmpstr2 = ' : Stefan refreezing with pond ice thickness'
            elseif (trim(frzpnd) == 'cesm') then
               tmpstr2 = ' : CESM refreezing empirical formula'
            else
               tmpstr2 = ' : unknown value'
            endif
            write(nu_diag,1030) ' frzpnd           = ', trim(frzpnd),trim(tmpstr2)
            write(nu_diag,1002) ' hs1              = ', hs1,' : snow depth of transition to pond ice'
         elseif (tr_pond_topo) then
            write(nu_diag,1010) ' tr_pond_topo     = ', tr_pond_topo,' : topo pond formulation'
            write(nu_diag,*) '     WARNING: dpnd history fields are turned off for topo ponds'
            write(nu_diag,1002) ' hp1              = ', hp1,' : critical ice lid thickness for topo ponds'
         elseif (tr_pond_sealvl) then
            write(nu_diag,1010) ' tr_pond_sealvl   = ', tr_pond_sealvl,' : sealvl pond formulation'
            write(nu_diag,1002) ' apnd_sl          = ', apnd_sl,' : equilibrium pond fraction'
         elseif (trim(shortwave) == 'ccsm3') then
            write(nu_diag,*) 'Pond effects on radiation are treated implicitly in the ccsm3 shortwave scheme'
         else
            write(nu_diag,*) 'Using default dEdd melt pond scheme for testing only'
         endif

         if (shortwave(1:4) == 'dEdd') then
            write(nu_diag,1002) ' hs0              = ', hs0,' : snow depth of transition to bare sea ice'
         endif

         write(nu_diag,1002) ' rfracmin         = ', rfracmin,' : minimum fraction of melt water added to ponds'
         write(nu_diag,1002) ' rfracmax         = ', rfracmax,' : maximum fraction of melt water added to ponds'

         write(nu_diag,*) ' '
         write(nu_diag,*) ' Snow redistribution/metamorphism tracers'
         write(nu_diag,*) '-----------------------------------------'
         if (tr_snow) then
            write(nu_diag,1010) ' tr_snow          = ', tr_snow, &
                                ' : advanced snow physics'
            if (snwredist(1:4) == 'none') then
               write(nu_diag,1030) ' snwredist        = ', trim(snwredist), &
                                   ' : Snow redistribution scheme turned off'
            else
               if (snwredist(1:4) == 'bulk') then
                  write(nu_diag,1030) ' snwredist        = ', trim(snwredist), &
                                      ' : Using bulk snow redistribution scheme'
               elseif (snwredist(1:6) == 'ITDrdg') then
                  write(nu_diag,1030) ' snwredist        = ', trim(snwredist), &
                                      ' : Using ridging based snow redistribution scheme'
                  write(nu_diag,1002) ' rhosnew          = ', rhosnew, &
                                      ' : new snow density (kg/m^3)'
                  write(nu_diag,1002) ' rhosmin          = ', rhosmin, &
                                      ' : minimum snow density (kg/m^3)'
                  write(nu_diag,1002) ' rhosmax          = ', rhosmax, &
                                      ' : maximum snow density (kg/m^3)'
                  write(nu_diag,1002) ' windmin          = ', windmin, &
                                      ' : minimum wind speed to compact snow (m/s)'
                  write(nu_diag,1002) ' drhosdwind       = ', drhosdwind, &
                                      ' : wind compaction factor (kg s/m^4)'
               endif
               write(nu_diag,1002) ' snwlvlfac        = ', snwlvlfac, &
                                   ' : fractional increase in snow depth for redistribution on ridges'
            endif
            if (.not. snwgrain) then
               write(nu_diag,1010) ' snwgrain         = ', snwgrain, &
                                   ' : Snow metamorphosis turned off'
            else
               write(nu_diag,1010) ' snwgrain         = ', snwgrain, &
                                   ' : Using snow metamorphosis scheme'
               write(nu_diag,1002) ' rsnw_tmax        = ', rsnw_tmax, &
                                   ' : maximum snow radius (10^-6 m)'
            endif
            write(nu_diag,1002) ' rsnw_fall        = ', rsnw_fall, &
                                ' : radius of new snow (10^-6 m)'
            if (snwgrain) then
               if (use_smliq_pnd) then
                  tmpstr2 = ' : Using liquid water in snow for melt ponds'
               else
                  tmpstr2 = ' : NOT using liquid water in snow for melt ponds'
               endif
                  write(nu_diag,1010) ' use_smliq_pnd    = ', use_smliq_pnd, trim(tmpstr2)
               if (snw_aging_table == 'test') then
                  tmpstr2 = ' : Using 5x5x1 test matrix of internallly defined snow aging parameters'
                  write(nu_diag,1030) ' snw_aging_table  = ', trim(snw_aging_table),trim(tmpstr2)
               elseif (snw_aging_table == 'snicar') then
                  tmpstr2 = ' : Reading 3D snow aging parameters from SNICAR file'
                  write(nu_diag,1030) ' snw_aging_table  = ', trim(snw_aging_table),trim(tmpstr2)
                  write(nu_diag,1031) ' snw_filename     = ',trim(snw_filename)
                  write(nu_diag,1031) ' snw_tau_fname    = ',trim(snw_tau_fname)
                  write(nu_diag,1031) ' snw_kappa_fname  = ',trim(snw_kappa_fname)
                  write(nu_diag,1031) ' snw_drdt0_fname  = ',trim(snw_drdt0_fname)
               elseif (snw_aging_table == 'file') then
                  tmpstr2 = ' : Reading 1D and 3D snow aging dimensions and parameters from external file'
                  write(nu_diag,1030) ' snw_aging_table  = ', trim(snw_aging_table),trim(tmpstr2)
                  write(nu_diag,1031) ' snw_filename     = ',trim(snw_filename)
                  write(nu_diag,1031) ' snw_rhos_fname   = ',trim(snw_rhos_fname)
                  write(nu_diag,1031) ' snw_Tgrd_fname   = ',trim(snw_Tgrd_fname)
                  write(nu_diag,1031) ' snw_T_fname      = ',trim(snw_T_fname)
                  write(nu_diag,1031) ' snw_tau_fname    = ',trim(snw_tau_fname)
                  write(nu_diag,1031) ' snw_kappa_fname  = ',trim(snw_kappa_fname)
                  write(nu_diag,1031) ' snw_drdt0_fname  = ',trim(snw_drdt0_fname)
               endif
            endif
         endif

         write(nu_diag,*) ' '
         write(nu_diag,*) ' Primary state variables, tracers'
         write(nu_diag,*) '   (excluding biogeochemistry)'
         write(nu_diag,*) '---------------------------------'
         write(nu_diag,*) 'Conserved properties (all tracers are conserved):'
         write(nu_diag,*) 'ice concentration, volume and enthalpy'
         write(nu_diag,*) 'snow volume and enthalpy'
         if (ktherm == 2)  write(nu_diag,1030) ' ice salinity'
         if (tr_fsd)       write(nu_diag,1010) ' tr_fsd           = ', tr_fsd,' : floe size distribution'
         if (tr_lvl)       write(nu_diag,1010) ' tr_lvl           = ', tr_lvl,' : ridging related tracers'
         if (tr_pond_lvl)  write(nu_diag,1010) ' tr_pond_lvl      = ', tr_pond_lvl,' : level-ice pond formulation'
         if (tr_pond_sealvl) write(nu_diag,1010) ' tr_pond_sealvl   = ', tr_pond_sealvl,' : sea level pond formulation'
         if (tr_pond_topo) write(nu_diag,1010) ' tr_pond_topo     = ', tr_pond_topo,' : topo pond formulation'
         if (tr_snow)      write(nu_diag,1010) ' tr_snow          = ', tr_snow,' : advanced snow physics'
         if (tr_iage)      write(nu_diag,1010) ' tr_iage          = ', tr_iage,' : chronological ice age'
         if (tr_FY)        write(nu_diag,1010) ' tr_FY            = ', tr_FY,' : first-year ice area'
         if (tr_iso)       write(nu_diag,1010) ' tr_iso           = ', tr_iso,' : diagnostic isotope tracers'
         if (tr_aero)      write(nu_diag,1010) ' tr_aero          = ', tr_aero,' : CESM aerosol tracers'
         write(nu_diag,*) 'Non-conserved properties:'
         write(nu_diag,*) 'ice surface temperature'
         write(nu_diag,*) 'ice velocity components and internal stress'

         write(nu_diag,*) ' '
         write(nu_diag,*) ' Other ice_in namelist parameters:'
         write(nu_diag,*) '===================================== '
         if (trim(runid) /= 'unknown') &
         write(nu_diag,1031) ' runid            = ', trim(runid)
         write(nu_diag,1031) ' version_name     = ', trim(version_name)
         write(nu_diag,1031) ' runtype          = ', trim(runtype)
         write(nu_diag,1021) ' year_init        = ', year_init
         write(nu_diag,1021) ' month_init       = ', month_init
         write(nu_diag,1021) ' day_init         = ', day_init
         write(nu_diag,1021) ' sec_init         = ', sec_init
         write(nu_diag,1021) ' istep0           = ', istep0
         write(nu_diag,1031) ' npt_unit         = ', trim(npt_unit)
         write(nu_diag,1021) ' npt              = ', npt
         write(nu_diag,1021) ' diagfreq         = ', diagfreq
         write(nu_diag,1011) ' print_global     = ', print_global
         write(nu_diag,1011) ' print_points     = ', print_points
         write(nu_diag,1011) ' debug_model      = ', debug_model
         write(nu_diag,1022) ' debug_model_step = ', debug_model_step
         write(nu_diag,1021) ' debug_model_i    = ', debug_model_i
         write(nu_diag,1021) ' debug_model_i    = ', debug_model_j
         write(nu_diag,1021) ' debug_model_iblk = ', debug_model_iblk
         write(nu_diag,1021) ' debug_model_task = ', debug_model_task
         write(nu_diag,1011) ' timer_stats      = ', timer_stats
         write(nu_diag,1011) ' memory_stats     = ', memory_stats
         write(nu_diag,1031) ' bfbflag          = ', trim(bfbflag)
         write(nu_diag,1021) ' numin            = ', numin
         write(nu_diag,1021) ' numax            = ', numax
         write(nu_diag,1011) ' grid_outfile     = ', grid_outfile
         write(nu_diag,1033) ' histfreq         = ', histfreq(1:max_nstrm-1)
         write(nu_diag,1023) ' histfreq_n       = ', histfreq_n(1:max_nstrm-1)
         write(nu_diag,1033) ' histfreq_base    = ', histfreq_base(1:max_nstrm-1)
         write(nu_diag,1013) ' hist_avg         = ', hist_avg(1:max_nstrm-1)
         write(nu_diag,1033) ' hist_suffix      = ', hist_suffix(1:max_nstrm-1)
         write(nu_diag,1031) ' history_dir      = ', trim(history_dir)
         write(nu_diag,1031) ' history_file     = ', trim(history_file)
         write(nu_diag,1021) ' history_precision= ', history_precision
         write(nu_diag,1031) ' history_format   = ', trim(history_format)
         write(nu_diag,1031) ' history_rearranger = ', trim(history_rearranger)
         write(nu_diag,1021) ' history_iotasks  = ', history_iotasks
         write(nu_diag,1021) ' history_root     = ', history_root
         write(nu_diag,1021) ' history_stride   = ', history_stride
         write(nu_diag,1031) ' hist_time_axis   = ', trim(hist_time_axis)
         write(nu_diag,1021) ' history_deflate  = ', history_deflate
         write(nu_diag,1023) ' history_chunksize= ', history_chunksize
         if (write_ic) then
            write(nu_diag,1039) ' Initial condition will be written in ', &
                               trim(incond_dir)
         endif
         write(nu_diag,1033) ' dumpfreq         = ', dumpfreq(1:max_nstrm-1)
         write(nu_diag,1023) ' dumpfreq_n       = ', dumpfreq_n(1:max_nstrm-1)
         write(nu_diag,1033) ' dumpfreq_base    = ', dumpfreq_base(1:max_nstrm-1)
         write(nu_diag,1011) ' dump_last        = ', dump_last
         write(nu_diag,1011) ' restart          = ', restart
         write(nu_diag,1031) ' restart_dir      = ', trim(restart_dir)
         write(nu_diag,1011) ' restart_ext      = ', restart_ext
         write(nu_diag,1031) ' restart_mod      = ', trim(restart_mod)
         write(nu_diag,1011) ' restart_coszen   = ', restart_coszen
         write(nu_diag,1031) ' restart_format   = ', trim(restart_format)
         write(nu_diag,1021) ' restart_deflate  = ', restart_deflate
         write(nu_diag,1023) ' restart_chunksize= ', restart_chunksize
!         write(nu_diag,1011) ' lcdf64           = ', lcdf64   ! deprecated
         write(nu_diag,1031) ' restart_rearranger = ', trim(restart_rearranger)
         write(nu_diag,1021) ' restart_iotasks  = ', restart_iotasks
         write(nu_diag,1021) ' restart_root     = ', restart_root
         write(nu_diag,1021) ' restart_stride   = ', restart_stride
         write(nu_diag,1031) ' restart_file     = ', trim(restart_file)
         write(nu_diag,1031) ' pointer_file     = ', trim(pointer_file)
         write(nu_diag,1011) ' use_restart_time = ', use_restart_time
         write(nu_diag,1031) ' ice_ic           = ', trim(ice_ic)
         if (trim(grid_type) /= 'rectangular' .or. &
             trim(grid_type) /= 'column') then
            write(nu_diag,1031) ' grid_file        = ', trim(grid_file)
            write(nu_diag,1031) ' gridcpl_file     = ', trim(gridcpl_file)
            write(nu_diag,1031) ' bathymetry_file  = ', trim(bathymetry_file)
            if (trim(kmt_type) == 'file') &
               write(nu_diag,1031) ' kmt_file         = ', trim(kmt_file)
         endif

         write(nu_diag,1011) ' conserv_check    = ', conserv_check

         write(nu_diag,1021) ' fyear_init       = ', fyear_init
         write(nu_diag,1021) ' ycycle           = ', ycycle
         write(nu_diag,1031) ' atm_data_type    = ', trim(atm_data_type)
         write(nu_diag,1031) ' atm_data_version = ', trim(atm_data_version)

         if (trim(atm_data_type) /= 'default') then
            write(nu_diag,1031) ' atm_data_dir     = ', trim(atm_data_dir)
            write(nu_diag,1031) ' precip_units     = ', trim(precip_units)
         elseif (trim(atm_data_type) == 'default') then
            write(nu_diag,1031) ' default_season   = ', trim(default_season)
         endif

         if (wave_spec) then
            write(nu_diag,1031) ' wave_spec_file   = ', trim(wave_spec_file)
         endif
         if (trim(bgc_data_type) == 'ncar' .or. &
             trim(ocn_data_type) == 'ncar') then
            write(nu_diag,1031) ' oceanmixed_file  = ', trim(oceanmixed_file)
         endif
         if (cpl_bgc) then
             write(nu_diag,*) 'BGC coupling is switched ON'
         else
             write(nu_diag,*) 'BGC coupling is switched OFF'
         endif
         write(nu_diag,1031) ' bgc_data_type    = ', trim(bgc_data_type)
         write(nu_diag,1031) ' fe_data_type     = ', trim(fe_data_type)
         write(nu_diag,1031) ' ice_data_type    = ', trim(ice_data_type)
         write(nu_diag,1031) ' ice_data_conc    = ', trim(ice_data_conc)
         write(nu_diag,1031) ' ice_data_dist    = ', trim(ice_data_dist)
         write(nu_diag,1031) ' bgc_data_dir     = ', trim(bgc_data_dir)
         write(nu_diag,1031) ' ocn_data_type    = ', trim(ocn_data_type)
         if (trim(bgc_data_type) /= 'default' .or. &
             trim(ocn_data_type) /= 'default') then
            write(nu_diag,1031) ' ocn_data_dir     = ', trim(ocn_data_dir)
            write(nu_diag,1011) ' restore_ocn      = ', restore_ocn
         endif
         write(nu_diag,1011) ' restore_ice      = ', restore_ice
         if (restore_ice .or. restore_ocn) &
         write(nu_diag,1021) ' trestore         = ', trestore

         write(nu_diag,*) ' '
         write(nu_diag,'(a31,2f8.2)') 'Diagnostic point 1: lat, lon =', &
                            latpnt(1), lonpnt(1)
         write(nu_diag,'(a31,2f8.2)') 'Diagnostic point 2: lat, lon =', &
                            latpnt(2), lonpnt(2)
         write(nu_diag,*) ' '

         ! tracer restarts
         write(nu_diag,1011) ' restart_age      = ', restart_age
         write(nu_diag,1011) ' restart_FY       = ', restart_FY
         write(nu_diag,1011) ' restart_lvl      = ', restart_lvl
         write(nu_diag,1011) ' restart_pond_lvl = ', restart_pond_lvl
         write(nu_diag,1011) ' restart_pond_sealvl = ', restart_pond_sealvl
         write(nu_diag,1011) ' restart_pond_topo= ', restart_pond_topo
         write(nu_diag,1011) ' restart_snow     = ', restart_snow
         write(nu_diag,1011) ' restart_iso      = ', restart_iso
         write(nu_diag,1011) ' restart_aero     = ', restart_aero
         write(nu_diag,1011) ' restart_fsd      = ', restart_fsd

         write(nu_diag,1021) ' n_iso            = ', n_iso
         write(nu_diag,1021) ' n_aero           = ', n_aero
         write(nu_diag,1021) ' n_zaero          = ', n_zaero
         write(nu_diag,1021) ' n_algae          = ', n_algae
         write(nu_diag,1021) ' n_doc            = ', n_doc
         write(nu_diag,1021) ' n_dic            = ', n_dic
         write(nu_diag,1021) ' n_don            = ', n_don
         write(nu_diag,1021) ' n_fed            = ', n_fed
         write(nu_diag,1021) ' n_fep            = ', n_fep
         write(nu_diag,*)    ' '

      endif                     ! my_task = master_task

      if (grid_format /=  'pop_nc'        .and. &
          grid_format /=  'mom_nc'        .and. &
          grid_format /=  'geosnc'        .and. &
          grid_format /=  'meshnc'        .and. &
          grid_format /=  'bin' ) then
         if (my_task == master_task) write(nu_diag,*) subname//' ERROR: unknown grid_format=',trim(grid_type)
         abort_list = trim(abort_list)//":67"
      endif

      if (grid_type  /=  'displaced_pole' .and. &
          grid_type  /=  'tripole'        .and. &
          grid_type  /=  'column'         .and. &
          grid_type  /=  'rectangular'    .and. &
          grid_type  /=  'regional'       .and. &
          grid_type  /=  'latlon') then
         if (my_task == master_task) write(nu_diag,*) subname//' ERROR: unknown grid_type=',trim(grid_type)
         abort_list = trim(abort_list)//":20"
      endif

      if (grid_ice /=  'B' .and. &
          grid_ice /=  'C' .and. &
          grid_ice /=  'CD' ) then
         if (my_task == master_task) write(nu_diag,*) subname//' ERROR: unknown grid_ice=',trim(grid_ice)
         abort_list = trim(abort_list)//":26"
      endif

      if (kmt_type  /=  'file'    .and. &
          kmt_type  /=  'channel' .and. &
          kmt_type  /=  'channel_oneeast' .and. &
          kmt_type  /=  'channel_onenorth' .and. &
          kmt_type  /=  'wall'    .and. &
          kmt_type  /=  'default' .and. &
          kmt_type  /=  'boxislands'.and. &
          kmt_type  /=  'none' ) then
         if (my_task == master_task) write(nu_diag,*) subname//' ERROR: unknown kmt_type=',trim(kmt_type)
         abort_list = trim(abort_list)//":27"
      endif

      if (grid_type  /=  'column'      .and. &
          grid_type  /=  'rectangular' .and. &
          kmt_type   /=  'file' .and. &
          kmt_type   /=  'none') then
         if (my_task == master_task) write(nu_diag,*) subname//' ERROR: need kmt file, kmt_type=',trim(kmt_type)
         abort_list = trim(abort_list)//":28"
      endif

      if (kdyn         == 1                .and. &
          evp_algorithm /= 'standard_2d'   .and. &
          evp_algorithm /= 'shared_mem_1d') then
         if (my_task == master_task) write(nu_diag,*) subname//' ERROR: unknown evp_algorithm=',trim(evp_algorithm)
             abort_list = trim(abort_list)//":21"
      endif

      if (abort_list /= "") then
         call flush_fileunit(nu_diag)
      endif
      call ice_barrier()
      if (abort_list /=  "") then
         write(nu_diag,*) subname,' ERROR: abort_list = ',trim(abort_list)
         call abort_ice (subname//' ABORTING on input ERRORS', &
            file=__FILE__, line=__LINE__)
      endif

      call icepack_init_parameters(ustar_min_in=ustar_min, albicev_in=albicev, albicei_in=albicei, &
         albsnowv_in=albsnowv, albsnowi_in=albsnowi, natmiter_in=natmiter, atmiter_conv_in=atmiter_conv, &
         emissivity_in=emissivity, snw_ssp_table_in=snw_ssp_table, hi_min_in=hi_min, &
         ahmax_in=ahmax, shortwave_in=shortwave, albedo_type_in=albedo_type, R_ice_in=R_ice, R_pnd_in=R_pnd, &
         R_snw_in=R_snw, dT_mlt_in=dT_mlt, rsnw_mlt_in=rsnw_mlt, &
         kstrength_in=kstrength, krdg_partic_in=krdg_partic, krdg_redist_in=krdg_redist, mu_rdg_in=mu_rdg, &
         atmbndy_in=atmbndy, calc_strair_in=calc_strair, formdrag_in=formdrag, highfreq_in=highfreq, &
         kitd_in=kitd, kcatbound_in=kcatbound, hs0_in=hs0, dpscale_in=dpscale, frzpnd_in=frzpnd, &
         rfracmin_in=rfracmin, rfracmax_in=rfracmax, pndaspect_in=pndaspect, hs1_in=hs1, hp1_in=hp1, &
         apnd_sl_in=apnd_sl, &
         ktherm_in=ktherm, calc_Tsfc_in=calc_Tsfc, conduct_in=conduct, semi_implicit_Tsfc_in=semi_implicit_Tsfc, &
         a_rapid_mode_in=a_rapid_mode, Rac_rapid_mode_in=Rac_rapid_mode, vapor_flux_correction_in=vapor_flux_correction, &
         floediam_in=floediam, hfrazilmin_in=hfrazilmin, Tliquidus_max_in=Tliquidus_max, &
         aspect_rapid_mode_in=aspect_rapid_mode, dSdt_slow_mode_in=dSdt_slow_mode, &
         phi_c_slow_mode_in=phi_c_slow_mode, phi_i_mushy_in=phi_i_mushy, conserv_check_in=conserv_check, &
         wave_spec_type_in = wave_spec_type, wave_spec_in=wave_spec, nfreq_in=nfreq, &
         update_ocn_f_in=update_ocn_f, cpl_frazil_in=cpl_frazil, congel_freeze_in=congel_freeze, &
         tfrz_option_in=tfrz_option, kalg_in=kalg, fbot_xfer_type_in=fbot_xfer_type, &
         saltflux_option_in=saltflux_option, ice_ref_salinity_in=ice_ref_salinity, &
         Pstar_in=Pstar, Cstar_in=Cstar, iceruf_in=iceruf, iceruf_ocn_in=iceruf_ocn, calc_dragio_in=calc_dragio, &
         windmin_in=windmin, drhosdwind_in=drhosdwind, &
         rsnw_fall_in=rsnw_fall, rsnw_tmax_in=rsnw_tmax, rhosnew_in=rhosnew, &
         snwlvlfac_in=snwlvlfac, rhosmin_in=rhosmin, rhosmax_in=rhosmax, &
         snwredist_in=snwredist, snwgrain_in=snwgrain, snw_aging_table_in=trim(snw_aging_table), &
         sw_redist_in=sw_redist, sw_frac_in=sw_frac, sw_dtemp_in=sw_dtemp, &
         tscale_pnd_drain_in=tscale_pnd_drain)
      call icepack_init_tracer_flags(tr_iage_in=tr_iage, tr_FY_in=tr_FY, &
         tr_lvl_in=tr_lvl, tr_iso_in=tr_iso, tr_aero_in=tr_aero, &
         tr_fsd_in=tr_fsd, tr_snow_in=tr_snow, tr_pond_in=tr_pond, &
         tr_pond_lvl_in=tr_pond_lvl, tr_pond_sealvl_in=tr_pond_sealvl, tr_pond_topo_in=tr_pond_topo)
      call icepack_init_tracer_sizes(ncat_in=ncat, nilyr_in=nilyr, nslyr_in=nslyr, nblyr_in=nblyr, &
         nfsd_in=nfsd, n_algae_in=n_algae, n_iso_in=n_iso, n_aero_in=n_aero, &
         n_DOC_in=n_DOC, n_DON_in=n_DON, &
         n_DIC_in=n_DIC, n_fed_in=n_fed, n_fep_in=n_fep, n_zaero_in=n_zaero)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

 1000    format (a20,1x,f13.6,1x,a) ! float
 1002    format (a20,5x,f9.2,1x,a)
 1003    format (a20,1x,g13.4,1x,a)
 1004    format (a20,1x,2g13.4,1x,a)
 1009    format (a20,1x,d13.6,1x,a)
 1010    format (a20,8x,l6,1x,a)  ! logical
 1011    format (a20,1x,l6)
 1013    format (a20,1x,6l3)
 1020    format (a20,8x,i6,1x,a)  ! integer
 1021    format (a20,1x,i6)
 1022    format (a20,1x,i12)
 1023    format (a20,1x,6i6)
 1030    format (a20,a14,1x,a)    ! character
 1031    format (a20,1x,a,a)
 1033    format (a20,1x,6a6)
 1039    format (a,1x,a,1x,a,1x,a)

      end subroutine input_data

!=======================================================================

! Initialize state for the itd model
!
! authors: C. M. Bitz, UW
!          William H. Lipscomb, LANL

      subroutine init_state

      use ice_blocks, only: block, get_block, nx_block, ny_block
      use ice_domain, only: nblocks, blocks_ice, halo_info
      use ice_domain_size, only: ncat, nilyr, nslyr, n_iso, n_aero, nfsd
      use ice_flux, only: sst, Tf, Tair, salinz, Tmltz
      use ice_grid, only: tmask, umask, ULON, TLAT, grid_ice, grid_average_X2Y
      use ice_boundary, only: ice_HaloUpdate
      use ice_constants, only: field_loc_Nface, field_loc_Eface, field_type_scalar
      use ice_state, only: trcr_depend, aicen, trcrn, vicen, vsnon, &
          aice0, aice, vice, vsno, trcr, aice_init, bound_state, &
          n_trcr_strata, nt_strata, trcr_base, uvel, vvel, &
          uvelN, vvelN, uvelE, vvelE

      integer (kind=int_kind) :: &
         ilo, ihi    , & ! physical domain indices
         jlo, jhi    , & ! physical domain indices
         iglob(nx_block), & ! global indices
         jglob(ny_block), & ! global indices
         i, j        , & ! horizontal indices
         k           , & ! vertical index
         it          , & ! tracer index
         iblk            ! block index


      integer (kind=int_kind) :: ntrcr
      logical (kind=log_kind) :: tr_iage, tr_FY, tr_lvl, tr_iso, tr_aero
      logical (kind=log_kind) :: tr_pond_lvl, tr_pond_topo, tr_pond_sealvl
      logical (kind=log_kind) :: tr_snow, tr_fsd
      integer (kind=int_kind) :: nt_Tsfc, nt_sice, nt_qice, nt_qsno, nt_iage, nt_FY
      integer (kind=int_kind) :: nt_alvl, nt_vlvl, nt_apnd, nt_hpnd, nt_ipnd
      integer (kind=int_kind) :: nt_smice, nt_smliq, nt_rhos, nt_rsnw
      integer (kind=int_kind) :: nt_isosno, nt_isoice, nt_aero, nt_fsd

      type (block) :: &
         this_block           ! block information for current block

      character(len=*), parameter :: subname='(init_state)'

      !-----------------------------------------------------------------

      call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
      call icepack_query_tracer_flags(tr_iage_out=tr_iage, tr_FY_out=tr_FY, &
        tr_lvl_out=tr_lvl, tr_iso_out=tr_iso, tr_aero_out=tr_aero, &
        tr_pond_lvl_out=tr_pond_lvl, tr_pond_topo_out=tr_pond_topo, &
        tr_pond_sealvl_out=tr_pond_sealvl, &
        tr_snow_out=tr_snow, tr_fsd_out=tr_fsd)
      call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_sice_out=nt_sice, &
        nt_qice_out=nt_qice, nt_qsno_out=nt_qsno, nt_iage_out=nt_iage, nt_fy_out=nt_fy, &
        nt_alvl_out=nt_alvl, nt_vlvl_out=nt_vlvl, nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, &
        nt_ipnd_out=nt_ipnd, nt_aero_out=nt_aero, nt_fsd_out=nt_fsd, &
        nt_smice_out=nt_smice, nt_smliq_out=nt_smliq, &
        nt_rhos_out=nt_rhos, nt_rsnw_out=nt_rsnw, &
        nt_isosno_out=nt_isosno, nt_isoice_out=nt_isoice)

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      !-----------------------------------------------------------------
      ! Check number of layers in ice and snow.
      !-----------------------------------------------------------------

      if (my_task == master_task) then

         if (nilyr < 1) then
            write(nu_diag,*) subname//' ERROR: Must have at least one ice layer'
            write(nu_diag,*) subname//' ERROR:   nilyr =', nilyr
            call abort_ice (error_message=subname//' Not enough ice layers', &
               file=__FILE__, line=__LINE__)
         endif

         if (nslyr < 1) then
            write(nu_diag,*) subname//' ERROR: Must have at least one snow layer'
            write(nu_diag,*) subname//' ERROR:   nslyr =', nslyr
            call abort_ice(error_message=subname//' Not enough snow layers', &
               file=__FILE__, line=__LINE__)
         endif

      endif      ! my_task

      !-----------------------------------------------------------------
      ! Set tracer types
      !-----------------------------------------------------------------

      trcr_depend(nt_Tsfc) = 0 ! ice/snow surface temperature
      do k = 1, nilyr
         trcr_depend(nt_sice + k - 1) = 1 ! volume-weighted ice salinity
         trcr_depend(nt_qice + k - 1) = 1 ! volume-weighted ice enthalpy
      enddo
      do k = 1, nslyr
         trcr_depend(nt_qsno + k - 1) = 2 ! volume-weighted snow enthalpy
      enddo
      if (tr_iage) trcr_depend(nt_iage)  = 1   ! volume-weighted ice age
      if (tr_FY)   trcr_depend(nt_FY)    = 0   ! area-weighted first-year ice area
      if (tr_lvl)  trcr_depend(nt_alvl)  = 0   ! level ice area
      if (tr_lvl)  trcr_depend(nt_vlvl)  = 1   ! level ice volume
      if (tr_pond_lvl) then
                   trcr_depend(nt_apnd)  = 2+nt_alvl   ! melt pond area
                   trcr_depend(nt_hpnd)  = 2+nt_apnd   ! melt pond depth
                   trcr_depend(nt_ipnd)  = 2+nt_apnd   ! refrozen pond lid
      endif
      if (tr_pond_topo .or. tr_pond_sealvl) then
                   trcr_depend(nt_apnd)  = 0           ! melt pond area
                   trcr_depend(nt_hpnd)  = 2+nt_apnd   ! melt pond depth
                   trcr_depend(nt_ipnd)  = 2+nt_apnd   ! refrozen pond lid
      endif
      if (tr_snow) then                                ! snow-volume-weighted snow tracers
         do k = 1, nslyr
            trcr_depend(nt_smice + k - 1) = 2          ! ice mass in snow
            trcr_depend(nt_smliq + k - 1) = 2          ! liquid mass in snow
            trcr_depend(nt_rhos  + k - 1) = 2          ! effective snow density
            trcr_depend(nt_rsnw  + k - 1) = 2          ! snow radius
         enddo
      endif
      if (tr_fsd) then
         do it = 1, nfsd
            trcr_depend(nt_fsd + it - 1) = 0    ! area-weighted floe size distribution
         enddo
      endif
      if (tr_iso) then  ! isotopes
         do it = 1, n_iso
            trcr_depend(nt_isosno+it-1) = 2     ! snow
            trcr_depend(nt_isoice+it-1) = 1     ! ice
         enddo
      endif
      if (tr_aero) then ! volume-weighted aerosols
         do it = 1, n_aero
            trcr_depend(nt_aero+(it-1)*4  ) = 2 ! snow
            trcr_depend(nt_aero+(it-1)*4+1) = 2 ! snow
            trcr_depend(nt_aero+(it-1)*4+2) = 1 ! ice
            trcr_depend(nt_aero+(it-1)*4+3) = 1 ! ice
         enddo
      endif

      trcr_base = c0

      do it = 1, ntrcr
         ! mask for base quantity on which tracers are carried
         if (trcr_depend(it) == 0) then      ! area
            trcr_base(it,1) = c1
         elseif (trcr_depend(it) == 1) then  ! ice volume
            trcr_base(it,2) = c1
         elseif (trcr_depend(it) == 2) then  ! snow volume
            trcr_base(it,3) = c1
         else
            trcr_base(it,1) = c1    ! default: ice area
            trcr_base(it,2) = c0
            trcr_base(it,3) = c0
         endif

         ! initialize number of underlying tracer layers
         n_trcr_strata(it) = 0
         ! default indices of underlying tracer layers
         nt_strata   (it,1) = 0
         nt_strata   (it,2) = 0
      enddo

      if (tr_pond_lvl) then
         n_trcr_strata(nt_apnd)   = 1       ! melt pond area
         nt_strata    (nt_apnd,1) = nt_alvl ! on level ice area
         n_trcr_strata(nt_hpnd)   = 2       ! melt pond depth
         nt_strata    (nt_hpnd,2) = nt_apnd ! on melt pond area
         nt_strata    (nt_hpnd,1) = nt_alvl ! on level ice area
         n_trcr_strata(nt_ipnd)   = 2       ! refrozen pond lid
         nt_strata    (nt_ipnd,2) = nt_apnd ! on melt pond area
         nt_strata    (nt_ipnd,1) = nt_alvl ! on level ice area
      endif
      if (tr_pond_topo .or. tr_pond_sealvl) then
         n_trcr_strata(nt_hpnd)   = 1       ! melt pond depth
         nt_strata    (nt_hpnd,1) = nt_apnd ! on melt pond area
         n_trcr_strata(nt_ipnd)   = 1       ! refrozen pond lid
         nt_strata    (nt_ipnd,1) = nt_apnd ! on melt pond area
      endif

      !-----------------------------------------------------------------
      ! Set state variables
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block, &
      !$OMP                     iglob,jglob)
      do iblk = 1, nblocks

         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi
         iglob = this_block%i_glob
         jglob = this_block%j_glob

         call set_state_var (nx_block,            ny_block,            &
                             ilo, ihi,            jlo, jhi,            &
                             iglob,               jglob,               &
                             ice_ic,              tmask(:,:,    iblk), &
                             umask(:,:,    iblk), &
                             ULON (:,:,    iblk), &
                             TLAT (:,:,    iblk), &
                             Tair (:,:,    iblk), sst  (:,:,    iblk), &
                             Tf   (:,:,    iblk),                      &
                             salinz(:,:,:, iblk), Tmltz(:,:,:,  iblk), &
                             aicen(:,:,  :,iblk), trcrn(:,:,:,:,iblk), &
                             vicen(:,:,  :,iblk), vsnon(:,:,  :,iblk), &
                             uvel (:,:,    iblk), vvel (:,:,    iblk))

      enddo                     ! iblk
      !$OMP END PARALLEL DO

      !-----------------------------------------------------------------
      ! ghost cell updates
      !-----------------------------------------------------------------

      call bound_state (aicen,        &
                        vicen, vsnon, &
                        ntrcr, trcrn)

      if (grid_ice == 'CD' .or. grid_ice == 'C') then

         call grid_average_X2Y('A',uvel,'U',uvelN,'N')
         call grid_average_X2Y('A',vvel,'U',vvelN,'N')
         call grid_average_X2Y('A',uvel,'U',uvelE,'E')
         call grid_average_X2Y('A',vvel,'U',vvelE,'E')

         ! Halo update on North, East faces
         call ice_HaloUpdate(uvelN, halo_info, &
                             field_loc_Nface, field_type_scalar, fillvalue=c0)
         call ice_HaloUpdate(vvelN, halo_info, &
                             field_loc_Nface, field_type_scalar, fillvalue=c0)

         call ice_HaloUpdate(uvelE, halo_info, &
                             field_loc_Eface, field_type_scalar, fillvalue=c0)
         call ice_HaloUpdate(vvelE, halo_info, &
                             field_loc_Eface, field_type_scalar, fillvalue=c0)

      endif

      !-----------------------------------------------------------------
      ! compute aggregate ice state and open water area
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk,it,i,j)
      do iblk = 1, nblocks

      do j = 1, ny_block
      do i = 1, nx_block
         aice(i,j,iblk) = c0
         vice(i,j,iblk) = c0
         vsno(i,j,iblk) = c0
         do it = 1, ntrcr
            trcr(i,j,it,iblk) = c0
         enddo

         if (tmask(i,j,iblk)) &
            call icepack_aggregate(aicen = aicen(i,j,:,iblk),     &
                                   trcrn = trcrn(i,j,:,:,iblk),   &
                                   vicen = vicen(i,j,:,iblk),     &
                                   vsnon = vsnon(i,j,:,iblk),     &
                                   aice  = aice (i,j,  iblk),     &
                                   trcr  = trcr (i,j,:,iblk),     &
                                   vice  = vice (i,j,  iblk),     &
                                   vsno  = vsno (i,j,  iblk),     &
                                   aice0 = aice0(i,j,  iblk),     &
                                   trcr_depend   = trcr_depend(:),   &
                                   trcr_base     = trcr_base(:,:),   &
                                   n_trcr_strata = n_trcr_strata(:), &
                                   nt_strata     = nt_strata(:,:),   &
                                   Tf            = Tf(i,j,iblk))

         aice_init(i,j,iblk) = aice(i,j,iblk)

      enddo
      enddo

      enddo                     ! iblk
      !$OMP END PARALLEL DO

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      end subroutine init_state

!=======================================================================

! Initialize state in each ice thickness category
!
! authors: C. M. Bitz
!          William H. Lipscomb, LANL

      subroutine set_state_var (nx_block, ny_block, &
                                ilo, ihi, jlo, jhi, &
                                iglob,    jglob,    &
                                ice_ic,   tmask,    &
                                umask, &
                                ULON,  &
                                TLAT,  &
                                Tair,     sst,  &
                                Tf,       &
                                salinz,   Tmltz, &
                                aicen,    trcrn, &
                                vicen,    vsnon, &
                                uvel,     vvel)


      use ice_arrays_column, only: hin_max
      use ice_domain_size, only: nilyr, nslyr, nx_global, ny_global, ncat
      use ice_grid, only: dxrect, dyrect
      use ice_forcing, only: ice_data_type, ice_data_conc, ice_data_dist

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ilo, ihi          , & ! physical domain indices
         jlo, jhi          , & !
         iglob(nx_block)   , & ! global indices
         jglob(ny_block)       !

      character(len=char_len_long), intent(in) :: &
         ice_ic      ! method of ice cover initialization

      logical (kind=log_kind), dimension (nx_block,ny_block), intent(in) :: &
         tmask  , & ! true for ice/ocean cells
         umask      ! for U points

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         ULON   , & ! longitude of velocity pts (radians)
         TLAT       ! latitude of temperature pts (radians)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         Tair    , & ! air temperature  (K)
         Tf      , & ! freezing temperature (C)
         sst         ! sea surface temperature (C)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr), intent(in) :: &
         salinz  , & ! initial salinity profile
         Tmltz       ! initial melting temperature profile

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), intent(out) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), intent(out), dimension (:,:,:,:) :: & ! (nx_block,ny_block,ntrcr,ncat)
         trcrn     ! ice tracers
                   ! 1: surface temperature of ice/snow (C)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out) :: &
         uvel    , & ! ice velocity B grid
         vvel        !

      ! local variables
      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij          , & ! horizontal index, combines i and j loops
         k           , & ! ice layer index
         n           , & ! thickness category index
         it          , & ! tracer index
         iedge       , & ! edge around big block
         jedge       , & ! edge around big block
         icells          ! number of cells initialized with ice

      logical (kind=log_kind) :: &
         in_slot, in_cyl ! boxslotcyl flags

      real (kind=dbl_kind) :: &  ! boxslotcyl parameters
         diam    , & ! cylinder diameter
         radius  , & ! cylinder radius
         center_x, & ! cylinder center
         center_y, &
         width   , & ! slot width
         length      ! slot height

      integer (kind=int_kind), dimension(nx_block*ny_block) :: &
         indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind) :: &
         Tsfc, asum, hbar, abar, puny, rhos, Lfresh, rad_to_deg, rsnw_fall, dist_ratio, Tffresh

      real (kind=dbl_kind), dimension(ncat) :: &
         ainit, hinit    ! initial area, thickness

      real (kind=dbl_kind), dimension(nilyr) :: &
         qin             ! ice enthalpy (J/m3)

      real (kind=dbl_kind), dimension(nslyr) :: &
         qsn             ! snow enthalpy (J/m3)

      real (kind=dbl_kind), parameter :: &
         hsno_init = 0.20_dbl_kind   , & ! initial snow thickness (m)
         edge_init_nh =  70._dbl_kind, & ! initial ice edge, N.Hem. (deg)
         edge_init_sh = -60._dbl_kind    ! initial ice edge, S.Hem. (deg)

      real (kind=dbl_kind) :: &  ! boxslotcyl
         pi             , & ! pi
         secday         , & ! seconds per day
         max_vel        , & ! max velocity
         domain_length  , & ! physical domain length
         period             ! rotational period

      logical (kind=log_kind) :: tr_brine, tr_lvl, tr_snow
      integer (kind=int_kind) :: ntrcr
      integer (kind=int_kind) :: nt_Tsfc, nt_qice, nt_qsno, nt_sice
      integer (kind=int_kind) :: nt_fbri, nt_alvl, nt_vlvl
      integer (kind=int_kind) :: nt_smice, nt_smliq, nt_rhos, nt_rsnw

      character(len=*), parameter :: subname='(set_state_var)'

      !-----------------------------------------------------------------

      call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
      call icepack_query_tracer_flags(tr_brine_out=tr_brine, tr_lvl_out=tr_lvl, &
        tr_snow_out=tr_snow)
      call icepack_query_tracer_indices( nt_Tsfc_out=nt_Tsfc, nt_qice_out=nt_qice, &
        nt_qsno_out=nt_qsno, nt_sice_out=nt_sice, &
        nt_fbri_out=nt_fbri, nt_alvl_out=nt_alvl, nt_vlvl_out=nt_vlvl, &
        nt_smice_out=nt_smice, nt_smliq_out=nt_smliq, &
        nt_rhos_out=nt_rhos, nt_rsnw_out=nt_rsnw)
      call icepack_query_parameters(rhos_out=rhos, Lfresh_out=Lfresh, puny_out=puny, &
        rad_to_deg_out=rad_to_deg, rsnw_fall_out=rsnw_fall, Tffresh_out=Tffresh)
      call icepack_query_parameters(secday_out=secday, pi_out=pi)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      indxi(:) = 0
      indxj(:) = 0

      ! Initialize state variables.
      ! If restarting, these values are overwritten.

      do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block
            aicen(i,j,n) = c0
            vicen(i,j,n) = c0
            vsnon(i,j,n) = c0
            if (tmask(i,j)) then
               trcrn(i,j,nt_Tsfc,n) = Tf(i,j)  ! surface temperature
            else
               trcrn(i,j,nt_Tsfc,n) = c0       ! at land grid cells (for clean history/restart files)
            endif
            if (ntrcr >= 2) then
               do it = 2, ntrcr
                  trcrn(i,j,it,n) = c0
               enddo
            endif
            if (tr_lvl)   trcrn(i,j,nt_alvl,n) = c1
            if (tr_lvl)   trcrn(i,j,nt_vlvl,n) = c1
            if (tr_brine) trcrn(i,j,nt_fbri,n) = c1
            do k = 1, nilyr
               trcrn(i,j,nt_sice+k-1,n) = salinz(i,j,k)
            enddo
            do k = 1, nslyr
               trcrn(i,j,nt_qsno+k-1,n) = -rhos * Lfresh
            enddo
            if (tr_snow) then
               do k = 1, nslyr
                  trcrn(i,j,nt_rsnw +k-1,n) = rsnw_fall
                  trcrn(i,j,nt_rhos +k-1,n) = rhos
                  trcrn(i,j,nt_smice+k-1,n) = rhos
                  trcrn(i,j,nt_smliq+k-1,n) = c0
               enddo               ! nslyr
            endif
         enddo
         enddo
      enddo

      if (trim(ice_ic) == 'internal') then

         !---------------------------------------------------------
         ! ice concentration/thickness
         !---------------------------------------------------------

         if (trim(ice_data_conc) == 'p5' .or. &
             trim(ice_data_conc) == 'p8' .or. &
             trim(ice_data_conc) == 'p9' .or. &
             trim(ice_data_conc) == 'c1' .or. &
             trim(ice_data_conc) == 'box2001') then

            if (trim(ice_data_conc) == 'p5') then
               hbar = c2  ! initial ice thickness
               abar = p5  ! initial ice concentration
            elseif (trim(ice_data_conc) == 'p8') then
               hbar = c1  ! initial ice thickness
               abar = 0.8_dbl_kind  ! initial ice concentration
            elseif (trim(ice_data_conc) == 'p9') then
               hbar = c1  ! initial ice thickness
               abar = 0.9_dbl_kind  ! initial ice concentration
            elseif (trim(ice_data_conc) == 'c1') then
               hbar = c1  ! initial ice thickness
               abar = c1  ! initial ice concentration
            elseif (trim(ice_data_conc) == 'box2001') then
               hbar = c2  ! initial ice thickness
               abar = p5  ! initial ice concentration
            endif

            do n = 1, ncat
               hinit(n) = c0
               ainit(n) = c0
               if (hbar > hin_max(n-1) .and. hbar <= hin_max(n)) then
                  hinit(n) = hbar
                  ainit(n) = abar
               endif
            enddo

         elseif (trim(ice_data_conc) == 'parabolic') then

            ! initial category areas in cells with ice
            hbar = c3  ! initial ice thickness with greatest area
                       ! Note: the resulting average ice thickness
                       ! tends to be less than hbar due to the
                       ! nonlinear distribution of ice thicknesses
            asum = c0
            do n = 1, ncat
               if (n < ncat) then
                  hinit(n) = p5*(hin_max(n-1) + hin_max(n)) ! m
               else                ! n=ncat
                  hinit(n) = (hin_max(n-1) + c1) ! m
               endif
               ! parabola, max at h=hbar, zero at h=0, 2*hbar
               ainit(n) = max(c0, (c2*hbar*hinit(n) - hinit(n)**2))
               asum = asum + ainit(n)
            enddo
            do n = 1, ncat
               ainit(n) = ainit(n) / (asum + puny/ncat) ! normalize
            enddo

         else

            call abort_ice(subname//'ERROR: ice_data_conc setting = '//trim(ice_data_conc), &
               file=__FILE__, line=__LINE__)

         endif ! ice_data_conc

         !---------------------------------------------------------
         ! location of ice
         !---------------------------------------------------------

         icells = 0

         if (trim(ice_data_type) == 'boxslotcyl') then

            ! Geometric configuration of the slotted cylinder
            diam     = p3 *dxrect*(nx_global-1)
            center_x = p5 *dxrect*(nx_global-1)
            center_y = p75*dyrect*(ny_global-1)
            radius   = p5*diam
            width    = p166*diam
            length   = c5*p166*diam

            icells = 0
            do j = jlo, jhi
            do i = ilo, ihi
               if (tmask(i,j)) then
                  ! check if grid point is inside slotted cylinder
                  in_slot = (dxrect*real(iglob(i)-1, kind=dbl_kind) >= center_x - width/c2) .and. &
                            (dxrect*real(iglob(i)-1, kind=dbl_kind) <= center_x + width/c2) .and. &
                            (dyrect*real(jglob(j)-1, kind=dbl_kind) >= center_y - radius) .and. &
                            (dyrect*real(jglob(j)-1, kind=dbl_kind) <= center_y + (length - radius))

                  in_cyl  = sqrt((dxrect*real(iglob(i)-1, kind=dbl_kind) - center_x)**c2 + &
                                 (dyrect*real(jglob(j)-1, kind=dbl_kind) - center_y)**c2) <= radius

                  if (in_cyl .and. .not. in_slot) then
                     icells = icells + 1
                     indxi(icells) = i
                     indxj(icells) = j
                  endif
               endif
            enddo
            enddo

         elseif (trim(ice_data_type) == 'uniform' .or. trim(ice_data_type) == 'box2001') then
            ! all cells not land mask are ice
            ! box2001 used to have a check for west of 50W, this was changed, so now box2001 is
            ! the same as uniform.  keep box2001 option for backwards compatibility.
            icells = 0
            do j = jlo, jhi
            do i = ilo, ihi
               if (tmask(i,j)) then
                  icells = icells + 1
                  indxi(icells) = i
                  indxj(icells) = j
               endif
            enddo
            enddo

         elseif (ice_data_type(1:7) == 'channel') then
            ! channel ice in center of domain in i direction
            icells = 0
            do j = jlo, jhi
            do i = ilo, ihi
               if (jglob(j) > ny_global/4 .and. jglob(j) < 3*nx_global/4) then
                  icells = icells + 1
                  indxi(icells) = i
                  indxj(icells) = j
               endif
            enddo
            enddo

         elseif (trim(ice_data_type) == 'block') then
            ! ice in 50% of domain, not at edges
            icells = 0
            iedge = int(real(nx_global,kind=dbl_kind) * 0.25) + 1
            jedge = int(real(ny_global,kind=dbl_kind) * 0.25) + 1
            do j = jlo, jhi
            do i = ilo, ihi
               if ((iglob(i) > iedge .and. iglob(i) < nx_global-iedge+1) .and. &
                   (jglob(j) > jedge .and. jglob(j) < ny_global-jedge+1)) then
                  icells = icells + 1
                  indxi(icells) = i
                  indxj(icells) = j
               endif
            enddo
            enddo

         elseif (trim(ice_data_type) == 'eastblock') then
            ! block on east half of domain in center of domain
            icells = 0
            do j = jlo, jhi
            do i = ilo, ihi
               if (jglob(j) > ny_global/4 .and. jglob(j) < 3*nx_global/4 .and. &
                   iglob(i) >= nx_global/2) then
                  icells = icells + 1
                  indxi(icells) = i
                  indxj(icells) = j
               endif
            enddo
            enddo

         elseif (trim(ice_data_type) == 'latsst') then

            !-----------------------------------------------------------------
            ! Place ice where ocean surface is cold.
            ! Note: If SST is not read from a file, then the ocean is assumed
            !       to be at its freezing point everywhere, and ice will
            !       extend to the prescribed edges.
            !-----------------------------------------------------------------

            icells = 0
            do j = jlo, jhi
            do i = ilo, ihi
               if (tmask(i,j)) then
                  ! place ice in high latitudes where ocean sfc is cold
#ifdef CESMCOUPLED
                  ! Option to use Tair instead.
                  if ( (Tair (i,j) <= Tffresh) .and. &
#else
                  if ( (sst (i,j) <= Tf(i,j)+p2) .and. &
#endif
                       (TLAT(i,j) < edge_init_sh/rad_to_deg .or. &
                        TLAT(i,j) > edge_init_nh/rad_to_deg) ) then
                     icells = icells + 1
                     indxi(icells) = i
                     indxj(icells) = j
                  endif            ! cold surface
               endif               ! tmask
            enddo                  ! i
            enddo                  ! j

         else

            call abort_ice(subname//'ERROR: ice_data_type setting = '//trim(ice_data_type), &
               file=__FILE__, line=__LINE__)

         endif                     ! ice_data_type

         !---------------------------------------------------------
         ! ice distribution
         !---------------------------------------------------------

         do n = 1, ncat

            ! ice volume, snow volume
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               aicen(i,j,n) = ainit(n)

               if (trim(ice_data_dist) == 'box2001') then
                  if (hinit(n) > c0) then
!                  ! varies linearly from 0 to 1 in x direction
                     aicen(i,j,n) = (real(iglob(i), kind=dbl_kind)-p5) &
                                  / (real(nx_global,kind=dbl_kind))
!                  ! constant slope from 0 to 0.5 in x direction
!                     aicen(i,j,n) = (real(iglob(i), kind=dbl_kind)-p5) &
!                                  / (real(nx_global,kind=dbl_kind)) * p5
!                  ! quadratic
!                     aicen(i,j,n) = max(c0,(real(iglob(i), kind=dbl_kind)-p5) &
!                                         / (real(nx_global,kind=dbl_kind)) &
!                                         * (real(jglob(j), kind=dbl_kind)-p5) &
!                                         / (real(ny_global,kind=dbl_kind)) * p5)
!                     aicen(i,j,n) = max(c0,(real(nx_global, kind=dbl_kind) &
!                                         -  real(iglob(i), kind=dbl_kind)-p5) &
!                                         / (real(nx_global,kind=dbl_kind)) &
!                                         * (real(ny_global, kind=dbl_kind) &
!                                         -  real(jglob(j), kind=dbl_kind)-p5) &
!                                         / (real(ny_global,kind=dbl_kind)) * p5)
                  endif

               elseif (trim(ice_data_dist) == 'gauss') then
                  if (hinit(n) > c0) then
                     dist_ratio = 8._dbl_kind * &
                                  sqrt((real(iglob(i),kind=dbl_kind)-real(nx_global+1,kind=dbl_kind)/c2)**2 + &
                                       (real(jglob(j),kind=dbl_kind)-real(ny_global+1,kind=dbl_kind)/c2)**2) / &
                                  sqrt((real(nx_global,kind=dbl_kind))**2 + &
                                       (real(ny_global,kind=dbl_kind))**2)
                     aicen(i,j,n) = ainit(n) * exp(-dist_ratio)
                  endif

               elseif (trim(ice_data_dist) == 'uniform') then

                  ! nothing extra to do

               else

                  call abort_ice(subname//'ERROR: ice_data_dist setting = '//trim(ice_data_dist), &
                  file=__FILE__, line=__LINE__)

               endif  ! ice_data_dist

               vicen(i,j,n) = hinit(n) * aicen(i,j,n) ! m
               vsnon(i,j,n) = min(aicen(i,j,n)*hsno_init,p2*vicen(i,j,n))

               call icepack_init_trcr(Tair  = Tair(i,j), Tf = Tf(i,j),  &
                                      Sprofile = salinz(i,j,:),         &
                                      Tprofile = Tmltz(i,j,:),          &
                                      Tsfc  = Tsfc,                     &
                                      qin   = qin(:),    qsn = qsn(:))

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
               ! brine fraction
               if (tr_brine) trcrn(i,j,nt_fbri,n) = c1

            enddo               ! ij
         enddo                  ! ncat

         !---------------------------------------------------------
         ! ice velocity
         ! these velocites are defined on B-grid
         !---------------------------------------------------------

         if (trim(ice_data_type) == 'boxslotcyl') then
            domain_length = dxrect*cm_to_m*nx_global
            period        = c12*secday               ! 12 days rotational period
            max_vel       = pi*domain_length/period

            do j = 1, ny_block
            do i = 1, nx_block

               if (umask(i,j)) then
                  uvel(i,j) =  c2*max_vel*(real(jglob(j), kind=dbl_kind) - p5) &
                            / real(ny_global - 1, kind=dbl_kind) - max_vel
                  vvel(i,j) = -c2*max_vel*(real(iglob(i), kind=dbl_kind) - p5) &
                            / real(nx_global - 1, kind=dbl_kind) + max_vel
               else
                  uvel(i,j) = c0
                  vvel(i,j) = c0
               endif
            enddo               ! j
            enddo               ! i
         else
            uvel = c0
            vvel = c0
         endif

      endif                     ! ice_ic

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      end subroutine set_state_var

!=======================================================================

      end module ice_init

!=======================================================================
