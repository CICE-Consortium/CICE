!  SVN:$Id: ice_init.F90 1228 2017-05-23 21:33:34Z tcraig $
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
      use ice_constants, only: c0, c1, c2, c3, p2, p5
      use ice_exit, only: abort_ice
      use ice_fileunits, only: nu_nml, nu_diag, nml_filename, diag_type, &
          ice_stdout, get_fileunit, release_fileunit, bfbflag, flush_fileunit
      use ice_fileunits, only: inst_suffix
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_aggregate
      use icepack_intfc, only: icepack_init_trcr
      use icepack_intfc, only: icepack_init_parameters
      use icepack_intfc, only: icepack_init_tracer_flags
      use icepack_intfc, only: icepack_init_tracer_numbers
      use icepack_intfc, only: icepack_init_tracer_indices
      use icepack_intfc, only: icepack_query_tracer_flags
      use icepack_intfc, only: icepack_query_tracer_numbers
      use icepack_intfc, only: icepack_query_tracer_indices
      use icepack_intfc, only: icepack_query_parameters

      implicit none

      character(len=char_len_long) :: &
         ice_ic      ! method of ice cover initialization
                     ! 'default'  => latitude and sst dependent
                     ! 'none'     => no ice
                     ! note:  restart = .true. overwrites

!=======================================================================

      contains

!=======================================================================

! Namelist variables, set to default values; may be altered
! at run time
!
! author Elizabeth C. Hunke, LANL

      subroutine input_data

      use ice_broadcast, only: broadcast_scalar, broadcast_array
      use ice_diagnostics, only: diag_file, print_global, print_points, latpnt, lonpnt
      use ice_domain_size, only: max_nstrm, nilyr, nslyr, max_ntrcr, ncat, n_aero
      use ice_calendar, only: year_init, istep0, histfreq, histfreq_n, &
                              dumpfreq, dumpfreq_n, diagfreq, nstreams, &
                              npt, dt, ndtd, days_per_year, use_leap_years, &
                              write_ic, dump_last
      use ice_arrays_column, only: oceanmixed_ice
      use ice_restart_column, only: restart_age, restart_FY, restart_lvl, &
          restart_pond_cesm, restart_pond_lvl, restart_pond_topo, restart_aero
      use ice_restart_shared, only: &
          restart, restart_ext, restart_dir, restart_file, pointer_file, &
          runid, runtype, use_restart_time, restart_format, lcdf64
      use ice_history_shared, only: hist_avg, history_dir, history_file, &
                             incond_dir, incond_file, version_name
      use ice_flux, only: update_ocn_f, l_mpond_fresh
      use ice_flux_bgc, only: cpl_bgc
      use ice_forcing, only: &
          ycycle,          fyear_init,    dbug, &
          atm_data_type,   atm_data_dir,  precip_units, &
          atm_data_format, ocn_data_format, &
          sss_data_type,   sst_data_type, ocn_data_dir, &
          oceanmixed_file, restore_sst,   trestore
      use ice_grid, only: grid_file, gridcpl_file, kmt_file, grid_type, grid_format
      use ice_dyn_shared, only: ndte, kdyn, revised_evp, yield_curve, &
                                basalstress, Ktens, e_ratio
      use ice_transport_driver, only: advection
      use ice_restoring, only: restore_ice
#ifdef CESMCOUPLED
      use shr_file_mod, only: shr_file_setIO
#endif

      ! local variables

      integer (kind=int_kind) :: &
        nml_error, & ! namelist i/o error flag
        n            ! loop index

      character (len=6) :: chartmp
      character (len=32) :: str

      logical :: exists

      real (kind=dbl_kind) :: ustar_min, albicev, albicei, albsnowv, albsnowi, &
        ahmax, R_ice, R_pnd, R_snw, dT_mlt, rsnw_mlt, emissivity, &
        mu_rdg, hs0, dpscale, rfracmin, rfracmax, pndaspect, hs1, hp1, &
        a_rapid_mode, Rac_rapid_mode, aspect_rapid_mode, dSdt_slow_mode, &
        phi_c_slow_mode, phi_i_mushy, kalg

      integer (kind=int_kind) :: ktherm, kstrength, krdg_partic, krdg_redist, natmiter, &
        kitd, kcatbound

      character (len=char_len) :: shortwave, albedo_type, conduct, fbot_xfer_type, &
        tfrz_option, frzpnd, atmbndy

      logical (kind=log_kind) :: calc_Tsfc, formdrag, highfreq, calc_strair

      integer (kind=int_kind) :: ntrcr
      logical (kind=log_kind) :: tr_iage, tr_FY, tr_lvl, tr_pond, tr_aero
      logical (kind=log_kind) :: tr_pond_cesm, tr_pond_lvl, tr_pond_topo
      integer (kind=int_kind) :: nt_Tsfc, nt_sice, nt_qice, nt_qsno, nt_iage, nt_FY
      integer (kind=int_kind) :: nt_alvl, nt_vlvl, nt_apnd, nt_hpnd, nt_ipnd, nt_aero

      real (kind=real_kind) :: rpcesm, rplvl, rptopo 
      real (kind=dbl_kind) :: Cf, puny
      integer :: abort_flag

      character(len=*), parameter :: subname='(input_data)'

      !-----------------------------------------------------------------
      ! Namelist variables.
      !-----------------------------------------------------------------

      namelist /setup_nml/ &
        days_per_year,  use_leap_years, year_init,       istep0,        &
        dt,             npt,            ndtd,                           &
        runtype,        runid,          bfbflag,                        &
        ice_ic,         restart,        restart_dir,     restart_file,  &
        restart_ext,    use_restart_time, restart_format, lcdf64,       &
        pointer_file,   dumpfreq,       dumpfreq_n,      dump_last,     &
        diagfreq,       diag_type,      diag_file,                      &
        print_global,   print_points,   latpnt,          lonpnt,        &
        dbug,           histfreq,       histfreq_n,      hist_avg,      &
        history_dir,    history_file,   cpl_bgc,                        &
        write_ic,       incond_dir,     incond_file,     version_name

      namelist /grid_nml/ &
        grid_format,    grid_type,       grid_file,     kmt_file,       &
        kcatbound,      gridcpl_file

      namelist /thermo_nml/ &
        kitd,           ktherm,          conduct,                       &
        a_rapid_mode,   Rac_rapid_mode,  aspect_rapid_mode,             &
        dSdt_slow_mode, phi_c_slow_mode, phi_i_mushy

      namelist /dynamics_nml/ &
        kdyn,           ndte,           revised_evp,    yield_curve,    &
        advection,                                                      &
        kstrength,      krdg_partic,    krdg_redist,    mu_rdg,         &
        e_ratio,        Ktens,          Cf,             basalstress

      namelist /shortwave_nml/ &
        shortwave,      albedo_type,                                    &
        albicev,        albicei,         albsnowv,      albsnowi,       &
        ahmax,          R_ice,           R_pnd,         R_snw,          &
        dT_mlt,         rsnw_mlt,        kalg

      namelist /ponds_nml/ &
        hs0,            dpscale,         frzpnd,                        &
        rfracmin,       rfracmax,        pndaspect,     hs1,            &
        hp1

      namelist /forcing_nml/ &
        atmbndy,        fyear_init,      ycycle,        atm_data_format,&
        atm_data_type,  atm_data_dir,    calc_strair,   calc_Tsfc,      &
        precip_units,   update_ocn_f,    l_mpond_fresh, ustar_min,      &
        fbot_xfer_type, emissivity,                                     &
        oceanmixed_ice, ocn_data_format, sss_data_type, sst_data_type,  &
        ocn_data_dir,   oceanmixed_file, restore_sst,   trestore,       &
        restore_ice,    formdrag,        highfreq,      natmiter,       &
        tfrz_option

      namelist /tracer_nml/   &
        tr_iage, restart_age, &
        tr_FY, restart_FY, &
        tr_lvl, restart_lvl, &
        tr_pond_cesm, restart_pond_cesm, &
        tr_pond_lvl, restart_pond_lvl, &
        tr_pond_topo, restart_pond_topo, &
        tr_aero, restart_aero

      !-----------------------------------------------------------------
      ! default values
      !-----------------------------------------------------------------

      abort_flag = 0

      call icepack_query_parameters(puny_out=puny)
! nu_diag not yet defined
!      call icepack_warnings_flush(nu_diag)
!      if (icepack_warnings_aborted()) call abort_ice(error_message=subname//'Icepack Abort0', &
!         file=__FILE__, line=__LINE__)

      days_per_year = 365    ! number of days in a year
      use_leap_years= .false.! if true, use leap years (Feb 29)
      year_init = 0          ! initial year
      istep0 = 0             ! no. of steps taken in previous integrations,
                             ! real (dumped) or imagined (to set calendar)
#ifndef CESMCOUPLED
      dt = 3600.0_dbl_kind   ! time step, s      
#endif
      npt = 99999            ! total number of time steps (dt) 
      diagfreq = 24          ! how often diag output is written
      print_points = .false. ! if true, print point data
      print_global = .true.  ! if true, print global diagnostic data
      bfbflag = .false.      ! if true, do bit-for-bit computations
      diag_type = 'stdout'
      diag_file = 'ice_diag.d'
      histfreq(1) = '1'      ! output frequency option for different streams
      histfreq(2) = 'h'      ! output frequency option for different streams
      histfreq(3) = 'd'      ! output frequency option for different streams
      histfreq(4) = 'm'      ! output frequency option for different streams
      histfreq(5) = 'y'      ! output frequency option for different streams
      histfreq_n(:) = 1      ! output frequency 
      hist_avg = .true.      ! if true, write time-averages (not snapshots)
      history_dir  = './'    ! write to executable dir for default
      history_file = 'iceh'  ! history file name prefix
      write_ic = .false.     ! write out initial condition
      cpl_bgc = .false.      ! history file name prefix
      incond_dir = history_dir ! write to history dir for default
      incond_file = 'iceh_ic'! file prefix
      dumpfreq='y'           ! restart frequency option
      dumpfreq_n = 1         ! restart frequency
      dump_last = .false.    ! write restart on last time step
      restart = .false.      ! if true, read restart files for initialization
      restart_dir  = './'     ! write to executable dir for default
      restart_file = 'iced'  ! restart file name prefix
      restart_ext  = .false. ! if true, read/write ghost cells
      use_restart_time = .true.     ! if true, use time info written in file
      pointer_file = 'ice.restart_file'
      restart_format = 'nc'  ! file format ('bin'=binary or 'nc'=netcdf or 'pio')
      lcdf64       = .false. ! 64 bit offset for netCDF
      ice_ic       = 'default'      ! latitude and sst-dependent
      grid_format  = 'bin'          ! file format ('bin'=binary or 'nc'=netcdf)
      grid_type    = 'rectangular'  ! define rectangular grid internally
      grid_file    = 'unknown_grid_file'
      gridcpl_file = 'unknown_gridcpl_file'
      kmt_file     = 'unknown_kmt_file'
      version_name = 'unknown_version_name'

      kitd = 1           ! type of itd conversions (0 = delta, 1 = linear)
      kcatbound = 1      ! category boundary formula (0 = old, 1 = new, etc)
      kdyn = 1           ! type of dynamics (1 = evp, 2 = eap)
      ndtd = 1           ! dynamic time steps per thermodynamic time step
      ndte = 120         ! subcycles per dynamics timestep:  ndte=dt_dyn/dte
      revised_evp = .false.  ! if true, use revised procedure for evp dynamics
      yield_curve = 'ellipse'
      kstrength = 1          ! 1 = Rothrock 75 strength, 0 = Hibler 79
      krdg_partic = 1        ! 1 = new participation, 0 = Thorndike et al 75
      krdg_redist = 1        ! 1 = new redistribution, 0 = Hibler 80
      mu_rdg = 3             ! e-folding scale of ridged ice, krdg_partic=1 (m^0.5)
      Cf = 17.0_dbl_kind     ! ratio of ridging work to PE change in ridging 
      basalstress= .false.   ! if true, basal stress for landfast is on
      Ktens = 0.0_dbl_kind   ! T=Ktens*P (tensile strength: see Konig and Holland, 2010)
      e_ratio = 2.0_dbl_kind ! EVP ellipse aspect ratio
      advection  = 'remap'   ! incremental remapping transport scheme
      shortwave = 'default'  ! 'default' or 'dEdd' (delta-Eddington)
      albedo_type = 'default'! or 'constant'
      ktherm = 1             ! 0 = 0-layer, 1 = BL99, 2 = mushy thermo
      conduct = 'bubbly'     ! 'MU71' or 'bubbly' (Pringle et al 2007)
      calc_Tsfc = .true.     ! calculate surface temperature
      update_ocn_f = .false. ! include fresh water and salt fluxes for frazil
      ustar_min = 0.005      ! minimum friction velocity for ocean heat flux (m/s)
      emissivity = 0.95      ! emissivity of snow and ice
      l_mpond_fresh = .false.     ! logical switch for including meltpond freshwater
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
      hp1       = 0.01_dbl_kind   ! critical pond lid thickness for topo ponds
      hs0       = 0.03_dbl_kind   ! snow depth for transition to bare sea ice (m)
      hs1       = 0.03_dbl_kind   ! snow depth for transition to bare pond ice (m)
      dpscale   = c1              ! alter e-folding time scale for flushing 
      frzpnd    = 'cesm'          ! melt pond refreezing parameterization
      rfracmin  = 0.15_dbl_kind   ! minimum retained fraction of meltwater
      rfracmax  = 0.85_dbl_kind   ! maximum retained fraction of meltwater
      pndaspect = 0.8_dbl_kind    ! ratio of pond depth to area fraction
      albicev   = 0.78_dbl_kind   ! visible ice albedo for h > ahmax
      albicei   = 0.36_dbl_kind   ! near-ir ice albedo for h > ahmax
      albsnowv  = 0.98_dbl_kind   ! cold snow albedo, visible
      albsnowi  = 0.70_dbl_kind   ! cold snow albedo, near IR
      ahmax     = 0.3_dbl_kind    ! thickness above which ice albedo is constant (m)
      atmbndy   = 'default'       ! or 'constant'

      fyear_init = 1900           ! first year of forcing cycle
      ycycle = 1                  ! number of years in forcing cycle
      atm_data_format = 'bin'     ! file format ('bin'=binary or 'nc'=netcdf)
      atm_data_type   = 'default'
      atm_data_dir    = ' '
      calc_strair     = .true.    ! calculate wind stress
      formdrag        = .false.   ! calculate form drag
      highfreq        = .false.   ! calculate high frequency RASM coupling
      natmiter        = 5         ! number of iterations for atm boundary layer calcs
      precip_units    = 'mks'     ! 'mm_per_month' or
                                  ! 'mm_per_sec' = 'mks' = kg/m^2 s
      tfrz_option     = 'mushy'   ! freezing temp formulation
      oceanmixed_ice  = .false.   ! if true, use internal ocean mixed layer
      ocn_data_format = 'bin'     ! file format ('bin'=binary or 'nc'=netcdf)
      sss_data_type   = 'default'
      sst_data_type   = 'default'
      ocn_data_dir    = ' '
      oceanmixed_file = 'unknown_oceanmixed_file' ! ocean forcing data
      restore_sst     = .false.   ! restore sst if true
      trestore        = 90        ! restoring timescale, days (0 instantaneous)
      restore_ice     = .false.   ! restore ice state on grid edges if true
      dbug      = .false.         ! true writes diagnostics for input forcing

      latpnt(1) =  90._dbl_kind   ! latitude of diagnostic point 1 (deg)
      lonpnt(1) =   0._dbl_kind   ! longitude of point 1 (deg)
      latpnt(2) = -65._dbl_kind   ! latitude of diagnostic point 2 (deg)
      lonpnt(2) = -45._dbl_kind   ! longitude of point 2 (deg)

#ifndef CESMCOUPLED
      runid   = 'unknown'   ! run ID used in CESM and for machine 'bering'
      runtype = 'initial'   ! run type: 'initial', 'continue'
#endif

      ! extra tracers
      tr_iage      = .false. ! ice age
      restart_age  = .false. ! ice age restart
      tr_FY        = .false. ! ice age
      restart_FY   = .false. ! ice age restart
      tr_lvl       = .false. ! level ice 
      restart_lvl  = .false. ! level ice restart
      tr_pond_cesm = .false. ! CESM melt ponds
      restart_pond_cesm = .false. ! melt ponds restart
      tr_pond_lvl  = .false. ! level-ice melt ponds
      restart_pond_lvl  = .false. ! melt ponds restart
      tr_pond_topo = .false. ! explicit melt ponds (topographic)
      restart_pond_topo = .false. ! melt ponds restart
      tr_aero      = .false. ! aerosols
      restart_aero = .false. ! aerosols restart

      ! mushy layer gravity drainage physics
      a_rapid_mode      =  0.5e-3_dbl_kind ! channel radius for rapid drainage mode (m)
      Rac_rapid_mode    =    10.0_dbl_kind ! critical Rayleigh number
      aspect_rapid_mode =     1.0_dbl_kind ! aspect ratio (larger is wider)
      dSdt_slow_mode    = -1.5e-7_dbl_kind ! slow mode drainage strength (m s-1 K-1)
      phi_c_slow_mode   =    0.05_dbl_kind ! critical liquid fraction porosity cutoff
      phi_i_mushy       =    0.85_dbl_kind ! liquid fraction of congelation ice

      !-----------------------------------------------------------------
      ! read from input file
      !-----------------------------------------------------------------

#ifdef CESMCOUPLED
      nml_filename  = 'ice_in'//trim(inst_suffix)
#endif

      call get_fileunit(nu_nml)

      if (my_task == master_task) then
         open (nu_nml, file=nml_filename, status='old',iostat=nml_error)
         if (nml_error /= 0) then
            nml_error = -1
         else
            nml_error =  1
         endif 

         do while (nml_error > 0)
            print*,'Reading setup_nml'
               read(nu_nml, nml=setup_nml,iostat=nml_error)
               if (nml_error /= 0) exit
            print*,'Reading grid_nml'
               read(nu_nml, nml=grid_nml,iostat=nml_error)
               if (nml_error /= 0) exit
            print*,'Reading tracer_nml'
               read(nu_nml, nml=tracer_nml,iostat=nml_error)
               if (nml_error /= 0) exit
            print*,'Reading thermo_nml'
               read(nu_nml, nml=thermo_nml,iostat=nml_error)
               if (nml_error /= 0) exit
            print*,'Reading dynamics_nml'
               read(nu_nml, nml=dynamics_nml,iostat=nml_error)
               if (nml_error /= 0) exit
            print*,'Reading shortwave_nml'
               read(nu_nml, nml=shortwave_nml,iostat=nml_error)
               if (nml_error /= 0) exit
            print*,'Reading ponds_nml'
               read(nu_nml, nml=ponds_nml,iostat=nml_error)
               if (nml_error /= 0) exit
            print*,'Reading forcing_nml'
               read(nu_nml, nml=forcing_nml,iostat=nml_error)
               if (nml_error /= 0) exit
         end do
         if (nml_error == 0) close(nu_nml)
      endif
      call broadcast_scalar(nml_error, master_task)
      if (nml_error /= 0) then
         call abort_ice(subname//'ERROR: reading namelist', &
            file=__FILE__, line=__LINE__)
      endif
      call release_fileunit(nu_nml)

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
         inquire(file='ice_modelio.nml'//trim(inst_suffix),exist=exists)
         if (exists) then
            call get_fileUnit(nu_diag)
            call shr_file_setIO('ice_modelio.nml'//trim(inst_suffix),nu_diag)
         end if
      else
         ! each task gets unique ice log filename when if test is true, for debugging
         if (1 == 0) then
            call get_fileUnit(nu_diag)
            write(str,'(a,i4.4)') "ice.log.task_",my_task
            open(nu_diag,file=str)
         endif
      end if
      if (trim(ice_ic) /= 'default' .and. trim(ice_ic) /= 'none') then
         restart = .true.
      end if
#else
      if (trim(diag_type) == 'file') call get_fileunit(nu_diag)
#endif

      !-----------------------------------------------------------------
      ! broadcast namelist settings
      !-----------------------------------------------------------------

      call broadcast_scalar(days_per_year,      master_task)
      call broadcast_scalar(use_leap_years,     master_task)
      call broadcast_scalar(year_init,          master_task)
      call broadcast_scalar(istep0,             master_task)
      call broadcast_scalar(dt,                 master_task)
      call broadcast_scalar(npt,                master_task)
      call broadcast_scalar(diagfreq,           master_task)
      call broadcast_scalar(print_points,       master_task)
      call broadcast_scalar(print_global,       master_task)
      call broadcast_scalar(bfbflag,            master_task)
      call broadcast_scalar(diag_type,          master_task)
      call broadcast_scalar(diag_file,          master_task)
      do n = 1, max_nstrm
         call broadcast_scalar(histfreq(n),     master_task)
      enddo  
      call broadcast_array(histfreq_n,          master_task)
      call broadcast_scalar(hist_avg,           master_task)
      call broadcast_scalar(history_dir,        master_task)
      call broadcast_scalar(history_file,       master_task)
      call broadcast_scalar(write_ic,           master_task)
      call broadcast_scalar(cpl_bgc,            master_task)
      call broadcast_scalar(incond_dir,         master_task)
      call broadcast_scalar(incond_file,        master_task)
      call broadcast_scalar(dumpfreq,           master_task)
      call broadcast_scalar(dumpfreq_n,         master_task)
      call broadcast_scalar(dump_last,          master_task)
      call broadcast_scalar(restart_file,       master_task)
      call broadcast_scalar(restart,            master_task)
      call broadcast_scalar(restart_dir,        master_task)
      call broadcast_scalar(restart_ext,        master_task)
      call broadcast_scalar(use_restart_time,   master_task)
      call broadcast_scalar(restart_format,     master_task)
      call broadcast_scalar(lcdf64,             master_task)
      call broadcast_scalar(pointer_file,       master_task)
      call broadcast_scalar(ice_ic,             master_task)
      call broadcast_scalar(grid_format,        master_task)
      call broadcast_scalar(grid_type,          master_task)
      call broadcast_scalar(grid_file,          master_task)
      call broadcast_scalar(gridcpl_file,       master_task)
      call broadcast_scalar(kmt_file,           master_task)
      call broadcast_scalar(kitd,               master_task)
      call broadcast_scalar(kcatbound,          master_task)
      call broadcast_scalar(kdyn,               master_task)
      call broadcast_scalar(ndtd,               master_task)
      call broadcast_scalar(ndte,               master_task)
      call broadcast_scalar(revised_evp,        master_task)
      call broadcast_scalar(yield_curve,        master_task)
      call broadcast_scalar(kstrength,          master_task)
      call broadcast_scalar(krdg_partic,        master_task)
      call broadcast_scalar(krdg_redist,        master_task)
      call broadcast_scalar(mu_rdg,             master_task)
      call broadcast_scalar(Cf,                 master_task)
      call broadcast_scalar(basalstress,        master_task)
      call broadcast_scalar(Ktens,              master_task)
      call broadcast_scalar(e_ratio,            master_task)
      call broadcast_scalar(advection,          master_task)
      call broadcast_scalar(shortwave,          master_task)
      call broadcast_scalar(albedo_type,        master_task)
      call broadcast_scalar(ktherm,             master_task)
      call broadcast_scalar(conduct,            master_task)
      call broadcast_scalar(R_ice,              master_task)
      call broadcast_scalar(R_pnd,              master_task)
      call broadcast_scalar(R_snw,              master_task)
      call broadcast_scalar(dT_mlt,             master_task)
      call broadcast_scalar(rsnw_mlt,           master_task)
      call broadcast_scalar(kalg,               master_task)
      call broadcast_scalar(hp1,                master_task)
      call broadcast_scalar(hs0,                master_task)
      call broadcast_scalar(hs1,                master_task)
      call broadcast_scalar(dpscale,            master_task)
      call broadcast_scalar(frzpnd,             master_task)
      call broadcast_scalar(rfracmin,           master_task)
      call broadcast_scalar(rfracmax,           master_task)
      call broadcast_scalar(pndaspect,          master_task)
      call broadcast_scalar(albicev,            master_task)
      call broadcast_scalar(albicei,            master_task)
      call broadcast_scalar(albsnowv,           master_task)
      call broadcast_scalar(albsnowi,           master_task)
      call broadcast_scalar(ahmax,              master_task)
      call broadcast_scalar(atmbndy,            master_task)
      call broadcast_scalar(fyear_init,         master_task)
      call broadcast_scalar(ycycle,             master_task)
      call broadcast_scalar(atm_data_format,    master_task)
      call broadcast_scalar(atm_data_type,      master_task)
      call broadcast_scalar(atm_data_dir,       master_task)
      call broadcast_scalar(calc_strair,        master_task)
      call broadcast_scalar(calc_Tsfc,          master_task)
      call broadcast_scalar(formdrag,           master_task)
      call broadcast_scalar(highfreq,           master_task)
      call broadcast_scalar(natmiter,           master_task)
      call broadcast_scalar(update_ocn_f,       master_task)
      call broadcast_scalar(l_mpond_fresh,      master_task)
      call broadcast_scalar(ustar_min,          master_task)
      call broadcast_scalar(emissivity,         master_task)
      call broadcast_scalar(fbot_xfer_type,     master_task)
      call broadcast_scalar(precip_units,       master_task)
      call broadcast_scalar(oceanmixed_ice,     master_task)
      call broadcast_scalar(tfrz_option,        master_task)
      call broadcast_scalar(ocn_data_format,    master_task)
      call broadcast_scalar(sss_data_type,      master_task)
      call broadcast_scalar(sst_data_type,      master_task)
      call broadcast_scalar(ocn_data_dir,       master_task)
      call broadcast_scalar(oceanmixed_file,    master_task)
      call broadcast_scalar(restore_sst,        master_task)
      call broadcast_scalar(trestore,           master_task)
      call broadcast_scalar(restore_ice,        master_task)
      call broadcast_scalar(dbug,               master_task)
      call broadcast_array (latpnt(1:2),        master_task)
      call broadcast_array (lonpnt(1:2),        master_task)
      call broadcast_scalar(runid,              master_task)
      call broadcast_scalar(runtype,            master_task)

      if (dbug) & ! else only master_task writes to file
      call broadcast_scalar(nu_diag,            master_task)

      ! tracers
      call broadcast_scalar(tr_iage,            master_task)
      call broadcast_scalar(restart_age,        master_task)
      call broadcast_scalar(tr_FY,              master_task)
      call broadcast_scalar(restart_FY,         master_task)
      call broadcast_scalar(tr_lvl,             master_task)
      call broadcast_scalar(restart_lvl,        master_task)
      call broadcast_scalar(tr_pond_cesm,       master_task)
      call broadcast_scalar(restart_pond_cesm,  master_task)
      call broadcast_scalar(tr_pond_lvl,        master_task)
      call broadcast_scalar(restart_pond_lvl,   master_task)
      call broadcast_scalar(tr_pond_topo,       master_task)
      call broadcast_scalar(restart_pond_topo,  master_task)
      call broadcast_scalar(tr_pond,            master_task)
      call broadcast_scalar(tr_aero,            master_task)
      call broadcast_scalar(restart_aero,       master_task)
      call broadcast_scalar(a_rapid_mode,       master_task)
      call broadcast_scalar(Rac_rapid_mode,     master_task)
      call broadcast_scalar(aspect_rapid_mode,  master_task)
      call broadcast_scalar(dSdt_slow_mode,     master_task)
      call broadcast_scalar(phi_c_slow_mode,    master_task)
      call broadcast_scalar(phi_i_mushy,        master_task)

#ifdef CESMCOUPLED
      pointer_file = trim(pointer_file) // trim(inst_suffix)
#endif

      !-----------------------------------------------------------------
      ! verify inputs
      !-----------------------------------------------------------------

      if (my_task == master_task) then
         if (trim(diag_type) == 'file') then
            write(ice_stdout,*) 'Diagnostic output will be in file ',diag_file
            open (nu_diag, file=diag_file, status='unknown')
         endif
         write(nu_diag,*) '--------------------------------'
         write(nu_diag,*) '  CICE model diagnostic output  '
         write(nu_diag,*) '--------------------------------'
         write(nu_diag,*) ' '
      endif

      if (trim(runtype) == 'continue' .and. .not.restart) then
         if (my_task == master_task) &
            write(nu_diag,*) 'WARNING: runtype=continue, setting restart=.true.'
         restart = .true.
      endif

      if (trim(runtype) /= 'continue' .and. restart .and. &
         (ice_ic == 'none' .or. ice_ic == 'default')) then
         if (my_task == master_task) &
            write(nu_diag,*) 'WARNING: runtype ne continue and ice_ic=none|default, setting restart=.false.'
         restart = .false.
      endif

      if (trim(runtype) /= 'continue' .and. (ice_ic == 'none' .or. ice_ic == 'default')) then
         if (my_task == master_task) &
            write(nu_diag,*) 'WARNING: ice_ic = none or default, setting restart flags to .false.'
         restart = .false.
         restart_aero =  .false. 
         restart_age =  .false. 
         restart_fy =  .false. 
         restart_lvl =  .false. 
         restart_pond_cesm =  .false. 
         restart_pond_lvl =  .false. 
         restart_pond_topo =  .false. 
! tcraig, probably needs to be uncommented when we can test bgc
!         restart_bgc =  .false. 
!         restart_hbrine =  .false. 
!         restart_zsal =  .false. 
! tcraig, OK to leave as true, needed for boxrestore case
!         restart_ext =  .false. 
      endif

      if (trim(runtype) == 'initial' .and. .not.(restart) .and. &
          ice_ic /= 'none' .and. ice_ic /= 'default') then
         if (my_task == master_task) then
            write(nu_diag,*) 'ERROR: runtype, restart, ice_ic are inconsistent:'
            write(nu_diag,*) 'ERROR:   runtype=',trim(runtype), 'restart=',restart, 'ice_ic=',trim(ice_ic)
            write(nu_diag,*) 'ERROR:   Please review user guide'
         endif
         abort_flag = 1
      endif

#ifndef ncdf
      if (grid_format /= 'bin' .or. atm_data_format /= 'bin' .or. ocn_data_format /= 'bin') then
         if (my_task == master_task) then
            write(nu_diag,*)'ERROR: ncdf CPP flag unset, data formats must be bin'
            write(nu_diag,*)'ERROR:   check grid_format, atm_data_format, ocn_data_format or set ncdf CPP'
         endif
         abort_flag = 2
      endif
#endif

      if (advection /= 'remap' .and. advection /= 'upwind' .and. advection /= 'none') then
         if (my_task == master_task) write(nu_diag,*)'ERROR: invalid advection=',trim(advection)
         abort_flag = 3
      endif

      if (ncat == 1 .and. kitd == 1) then
         if (my_task == master_task) then
            write(nu_diag,*) 'ERROR: kitd incompatability: ncat=1 and kitd=1'
            write(nu_diag,*) 'ERROR:   Remapping the ITD is not allowed for ncat=1.'
            write(nu_diag,*) 'ERROR:   Use kitd = 0 (delta function ITD) with kcatbound = 0'
            write(nu_diag,*) 'ERROR:   or for column configurations use kcatbound = -1'
         endif
         abort_flag = 4
      endif

      if (ncat /= 1 .and. kcatbound == -1) then
         if (my_task == master_task) then
            write(nu_diag,*) 'ERROR: ITD required for ncat > 1'
            write(nu_diag,*) 'ERROR:   ncat=',ncat,' kcatbound=',kcatbound
            write(nu_diag,*) 'ERROR:   Please review user guide'
         endif
         abort_flag = 5
      endif

      if (kdyn == 2 .and. revised_evp) then
         if (my_task == master_task) then
            write(nu_diag,*) 'WARNING: revised_evp = T with EAP dynamics'
            write(nu_diag,*) 'WARNING:   revised_evp is ignored'
         endif
         revised_evp = .false.
      endif

      rpcesm = c0
      rplvl  = c0
      rptopo = c0
      if (tr_pond_cesm) rpcesm = c1
      if (tr_pond_lvl ) rplvl  = c1
      if (tr_pond_topo) rptopo = c1

      tr_pond = .false. ! explicit melt ponds
      if (rpcesm + rplvl + rptopo > puny) tr_pond = .true.

      if (rpcesm + rplvl + rptopo > c1 + puny) then
         if (my_task == master_task) then
            write(nu_diag,*) 'ERROR: Must use only one melt pond scheme'
         endif
         abort_flag = 6
      endif

      if (tr_pond_lvl .and. .not. tr_lvl) then
         if (my_task == master_task) then
            write(nu_diag,*) 'WARNING: tr_pond_lvl=T but tr_lvl=F'
            write(nu_diag,*) 'WARNING: Setting tr_lvl=T'
         endif
         tr_lvl = .true.
      endif

! tcraig - this was originally implemented by resetting hs0=0. EH says it might be OK
! to not reset it but extra calculations are done and it might not be bfb.  In our
! testing, we should explicitly set hs0 to 0. when setting tr_pond_lvl=T, and otherwise
! this will abort (safest option until additional testing is done)
      if (tr_pond_lvl .and. abs(hs0) > puny) then
         if (my_task == master_task) then
            write(nu_diag,*) 'ERROR: tr_pond_lvl=T and hs0 /= 0'
         endif
         abort_flag = 7
      endif

      if (trim(shortwave) /= 'dEdd' .and. tr_pond .and. calc_tsfc) then
         if (my_task == master_task) then
            write(nu_diag,*) 'ERROR: tr_pond=T, calc_tsfc=T, invalid shortwave'
            write(nu_diag,*) 'ERROR:   Must use shortwave=dEdd'
         endif
         abort_flag = 8
      endif

      if (tr_aero .and. n_aero==0) then
         if (my_task == master_task) then
            write(nu_diag,*) 'ERROR: aerosols activated but'
            write(nu_diag,*) 'ERROR:   not allocated in tracer array.'
            write(nu_diag,*) 'ERROR:   Activate in compilation script.'
         endif
         abort_flag = 9
      endif

      if (trim(shortwave) /= 'dEdd' .and. tr_aero) then
         if (my_task == master_task) then
            write(nu_diag,*) 'ERROR: tr_aero=T, invalid shortwave'
            write(nu_diag,*) 'ERROR:   Must use shortwave=dEdd'
         endif
         abort_flag = 10
      endif

      if ((rfracmin < -puny .or. rfracmin > c1+puny) .or. &
          (rfracmax < -puny .or. rfracmax > c1+puny) .or. &
          (rfracmin > rfracmax)) then
         if (my_task == master_task) then
            write(nu_diag,*) 'ERROR: rfracmin, rfracmax must be between 0 and 1'
            write(nu_diag,*) 'ERROR:   and rfracmax >= rfracmin'
         endif
         abort_flag = 11
      endif
      rfracmin = min(max(rfracmin,c0),c1)
      rfracmax = min(max(rfracmax,c0),c1)

      if (trim(atm_data_type) == 'monthly' .and. calc_strair) then
         if (my_task == master_task) write(nu_diag,*)'ERROR: atm_data_type=monthly and calc_strair=T'
         abort_flag = 12
      endif

      if (ktherm == 2 .and. .not. calc_Tsfc) then
         if (my_task == master_task) write(nu_diag,*) 'ERROR: ktherm = 2 and calc_Tsfc=F'
         abort_flag = 13
      endif

! tcraig, is it really OK for users to run inconsistently?
      if (ktherm == 1 .and. trim(tfrz_option) /= 'linear_salt') then
         if (my_task == master_task) then
            write(nu_diag,*) 'WARNING: ktherm = 1 and tfrz_option = ',trim(tfrz_option)
            write(nu_diag,*) 'WARNING:   For consistency, set tfrz_option = linear_salt'
         endif
      endif
      if (ktherm == 2 .and. trim(tfrz_option) /= 'mushy') then
         if (my_task == master_task) then
            write(nu_diag,*) 'WARNING: ktherm = 2 and tfrz_option = ',trim(tfrz_option)
            write(nu_diag,*) 'WARNING:   For consistency, set tfrz_option = mushy'
         endif
      endif
!tcraig

      if (formdrag) then
         if (trim(atmbndy) == 'constant') then
            if (my_task == master_task) write(nu_diag,*) 'ERROR: formdrag=T and atmbndy=constant'
            abort_flag = 14
         endif

         if (.not. calc_strair) then
            if (my_task == master_task) write(nu_diag,*) 'ERROR: formdrag=T and calc_strair=F'
            abort_flag = 15
         endif

         if (tr_pond_cesm) then
            if (my_task == master_task) write(nu_diag,*)'ERROR: formdrag=T and frzpnd=cesm'
            abort_flag = 16
         endif

         if (.not. tr_lvl) then
            if (my_task == master_task) write(nu_diag,*) 'ERROR: formdrag=T and tr_lvl=F'
            abort_flag = 17
         endif
      endif

      if (trim(fbot_xfer_type) == 'Cdn_ocn' .and. .not. formdrag)  then
         if (my_task == master_task) write(nu_diag,*) 'ERROR: formdrag=F and fbot_xfer_type=Cdn_ocn'
         abort_flag = 18
      endif

      call icepack_init_parameters(Cf_in=Cf)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname//'Icepack Abort1', &
         file=__FILE__, line=__LINE__)

      !-----------------------------------------------------------------
      ! spew
      !-----------------------------------------------------------------

      if (my_task == master_task) then

         write(nu_diag,*) ' Document ice_in namelist parameters:'
         write(nu_diag,*) ' ==================================== '
         write(nu_diag,*) ' '
         if (trim(runid) /= 'unknown') &
         write(nu_diag,*)    ' runid                     = ', &
                               trim(runid)
         write(nu_diag,1030) ' runtype                   = ', &
                               trim(runtype)
         write(nu_diag,1020) ' days_per_year             = ', days_per_year
         write(nu_diag,1010) ' use_leap_years            = ', use_leap_years
         write(nu_diag,1020) ' year_init                 = ', year_init
         write(nu_diag,1020) ' istep0                    = ', istep0
         write(nu_diag,1000) ' dt                        = ', dt
         write(nu_diag,1020) ' npt                       = ', npt
         write(nu_diag,1020) ' diagfreq                  = ', diagfreq
         write(nu_diag,1010) ' print_global              = ', print_global
         write(nu_diag,1010) ' print_points              = ', print_points
         write(nu_diag,1010) ' bfbflag                   = ', bfbflag
         write(nu_diag,1050) ' histfreq                  = ', histfreq(:)
         write(nu_diag,1040) ' histfreq_n                = ', histfreq_n(:)
         write(nu_diag,1010) ' hist_avg                  = ', hist_avg
         if (.not. hist_avg) write(nu_diag,*) 'History data will be snapshots'
         write(nu_diag,*)    ' history_dir               = ', &
                               trim(history_dir)
         write(nu_diag,*)    ' history_file              = ', &
                               trim(history_file)
         if (write_ic) then
            write(nu_diag,*) 'Initial condition will be written in ', &
                               trim(incond_dir)
         endif
         write(nu_diag,1030) ' dumpfreq                  = ', &
                               trim(dumpfreq)
         write(nu_diag,1020) ' dumpfreq_n                = ', dumpfreq_n
         write(nu_diag,1010) ' dump_last                 = ', dump_last
         write(nu_diag,1010) ' restart                   = ', restart
         write(nu_diag,*)    ' restart_dir               = ', &
                               trim(restart_dir)
         write(nu_diag,*)    ' restart_ext               = ', restart_ext
         write(nu_diag,*)    ' restart_format            = ', &
                               trim(restart_format)
         write(nu_diag,*)    ' lcdf64                    = ', &
                               lcdf64
         write(nu_diag,*)    ' restart_file              = ', &
                               trim(restart_file)
         write(nu_diag,*)    ' pointer_file              = ', &
                               trim(pointer_file)
         write(nu_diag,*)    ' use_restart_time          = ', use_restart_time
         write(nu_diag,*)    ' ice_ic                    = ', &
                               trim(ice_ic)
         write(nu_diag,*)    ' grid_type                 = ', &
                               trim(grid_type)
         if (trim(grid_type) /= 'rectangular' .or. &
             trim(grid_type) /= 'column') then
            write(nu_diag,*) ' grid_file                 = ', &
                               trim(grid_file)
            write(nu_diag,*) ' gridcpl_file              = ', &
                               trim(gridcpl_file)
            write(nu_diag,*) ' kmt_file                  = ', &
                               trim(kmt_file)
         endif
         write(nu_diag,1020) ' kitd                      = ', kitd
         write(nu_diag,1020) ' kcatbound                 = ', &
                               kcatbound
         if (kdyn == 1) then
           write(nu_diag,1021) ' kdyn                      = ','evp ', kdyn
         elseif (kdyn == 2) then
           write(nu_diag,1021) ' kdyn                      = ','eap ', kdyn
         else
           write(nu_diag,1020) ' kdyn                      = ', kdyn
         endif
         write(nu_diag,1020) ' ndtd                      = ', ndtd
         write(nu_diag,1020) ' ndte                      = ', ndte
         write(nu_diag,1010) ' revised_evp               = ', &
                               revised_evp
         if (kdyn == 1) &
         write(nu_diag,*)    ' yield_curve               = ', &
                               trim(yield_curve)
         write(nu_diag,1020) ' kstrength                 = ', kstrength
         write(nu_diag,1020) ' krdg_partic               = ', &
                               krdg_partic
         write(nu_diag,1020) ' krdg_redist               = ', &
                               krdg_redist
         if (krdg_redist == 1) &
         write(nu_diag,1000) ' mu_rdg                    = ', mu_rdg
         if (kstrength == 1) &
         write(nu_diag,1000) ' Cf                        = ', Cf
         write(nu_diag,1010) ' basalstress               = ', basalstress
         write(nu_diag,1005) ' Ktens                     = ', Ktens
         write(nu_diag,1005) ' e_ratio                   = ', e_ratio    
         write(nu_diag,1030) ' advection                 = ', &
                               trim(advection)
         write(nu_diag,1030) ' shortwave                 = ', &
                               trim(shortwave)
         if (cpl_bgc) then
             write(nu_diag,1000) ' BGC coupling is switched ON'
         else
             write(nu_diag,1000) ' BGC coupling is switched OFF'
          endif

         if (trim(shortwave) == 'dEdd') then
         write(nu_diag,1000) ' R_ice                     = ', R_ice
         write(nu_diag,1000) ' R_pnd                     = ', R_pnd
         write(nu_diag,1000) ' R_snw                     = ', R_snw
         write(nu_diag,1000) ' dT_mlt                    = ', dT_mlt
         write(nu_diag,1000) ' rsnw_mlt                  = ', rsnw_mlt
         write(nu_diag,1000) ' kalg                      = ', kalg
         write(nu_diag,1000) ' hp1                       = ', hp1
         write(nu_diag,1000) ' hs0                       = ', hs0
         else
         write(nu_diag,1030) ' albedo_type               = ', &
                               trim(albedo_type)
         write(nu_diag,1000) ' albicev                   = ', albicev
         write(nu_diag,1000) ' albicei                   = ', albicei
         write(nu_diag,1000) ' albsnowv                  = ', albsnowv
         write(nu_diag,1000) ' albsnowi                  = ', albsnowi
         write(nu_diag,1000) ' ahmax                     = ', ahmax
         endif

         write(nu_diag,1000) ' rfracmin                  = ', rfracmin
         write(nu_diag,1000) ' rfracmax                  = ', rfracmax
         if (tr_pond_lvl) then
         write(nu_diag,1000) ' hs1                       = ', hs1
         write(nu_diag,1000) ' dpscale                   = ', dpscale
         write(nu_diag,1030) ' frzpnd                    = ', trim(frzpnd)
         endif
         if (tr_pond .and. .not. tr_pond_lvl) &
         write(nu_diag,1000) ' pndaspect                 = ', pndaspect

         write(nu_diag,1020) ' ktherm                    = ', ktherm
         if (ktherm == 1) &
         write(nu_diag,1030) ' conduct                   = ', conduct
         if (ktherm == 2) then
         write(nu_diag,1005) ' a_rapid_mode              = ', a_rapid_mode
         write(nu_diag,1005) ' Rac_rapid_mode            = ', Rac_rapid_mode
         write(nu_diag,1005) ' aspect_rapid_mode         = ', aspect_rapid_mode
         write(nu_diag,1005) ' dSdt_slow_mode            = ', dSdt_slow_mode
         write(nu_diag,1005) ' phi_c_slow_mode           = ', phi_c_slow_mode
         write(nu_diag,1005) ' phi_i_mushy               = ', phi_i_mushy
         endif

         write(nu_diag,1030) ' atmbndy                   = ', &
                               trim(atmbndy)
         write(nu_diag,1010) ' formdrag                  = ', formdrag
         write(nu_diag,1010) ' highfreq                  = ', highfreq
         write(nu_diag,1020) ' natmiter                  = ', natmiter
         write(nu_diag,1010) ' calc_strair               = ', calc_strair
         write(nu_diag,1010) ' calc_Tsfc                 = ', calc_Tsfc

         write(nu_diag,1020) ' fyear_init                = ', &
                               fyear_init
         write(nu_diag,1020) ' ycycle                    = ', ycycle
         write(nu_diag,*)    ' atm_data_type             = ', &
                               trim(atm_data_type)
         if (trim(atm_data_type) /= 'default') then
            write(nu_diag,*) ' atm_data_dir              = ', &
                               trim(atm_data_dir)
            write(nu_diag,*) ' precip_units              = ', &
                               trim(precip_units)
         endif 

         write(nu_diag,1010) ' update_ocn_f              = ', update_ocn_f
         write(nu_diag,1010) ' l_mpond_fresh             = ', l_mpond_fresh
         write(nu_diag,1005) ' ustar_min                 = ', ustar_min
         write(nu_diag,1005) ' emissivity                = ', emissivity
         write(nu_diag, *)   ' fbot_xfer_type            = ', &
                               trim(fbot_xfer_type)
         write(nu_diag,1010) ' oceanmixed_ice            = ', &
                               oceanmixed_ice
         write(nu_diag,*)    ' tfrz_option               = ', &
                               trim(tfrz_option)
         if (trim(sss_data_type) == 'ncar' .or. &
             trim(sst_data_type) == 'ncar') then
            write(nu_diag,*) ' oceanmixed_file           = ', &
                               trim(oceanmixed_file)
         endif
         write(nu_diag,*)    ' sss_data_type             = ', &
                               trim(sss_data_type)
         write(nu_diag,*)    ' sst_data_type             = ', &
                               trim(sst_data_type)
         if (trim(sss_data_type) /= 'default' .or. &
             trim(sst_data_type) /= 'default') then
            write(nu_diag,*) ' ocn_data_dir              = ', &
                               trim(ocn_data_dir)
            write(nu_diag,1010) ' restore_sst               = ', &
                               restore_sst
         endif 
         write(nu_diag,1010) ' restore_ice               = ', &
                               restore_ice
         if (restore_ice .or. restore_sst) &
         write(nu_diag,1020) ' trestore                  = ', trestore
 
#ifdef coupled
         if( oceanmixed_ice ) then
            write(nu_diag,*) 'WARNING ** WARNING ** WARNING ** WARNING '
            write(nu_diag,*) 'WARNING: coupled CPP and oceanmixed_ice namelist are BOTH ON'
            write(nu_diag,*) 'WARNING:   Ocean data received from coupler will'
            write(nu_diag,*) 'WARNING:   be altered by mixed layer routine!'
            write(nu_diag,*) 'WARNING ** WARNING ** WARNING ** WARNING '
            write(nu_diag,*) ' '
         endif
#endif

         write(nu_diag,*) ' '
         write(nu_diag,'(a30,2f8.2)') 'Diagnostic point 1: lat, lon =', &
                            latpnt(1), lonpnt(1)
         write(nu_diag,'(a30,2f8.2)') 'Diagnostic point 2: lat, lon =', &
                            latpnt(2), lonpnt(2)

         ! tracers
         write(nu_diag,1010) ' tr_iage                   = ', tr_iage
         write(nu_diag,1010) ' restart_age               = ', restart_age
         write(nu_diag,1010) ' tr_FY                     = ', tr_FY
         write(nu_diag,1010) ' restart_FY                = ', restart_FY
         write(nu_diag,1010) ' tr_lvl                    = ', tr_lvl
         write(nu_diag,1010) ' restart_lvl               = ', restart_lvl
         write(nu_diag,1010) ' tr_pond_cesm              = ', tr_pond_cesm
         write(nu_diag,1010) ' restart_pond_cesm         = ', restart_pond_cesm
         write(nu_diag,1010) ' tr_pond_lvl               = ', tr_pond_lvl
         write(nu_diag,1010) ' restart_pond_lvl          = ', restart_pond_lvl
         write(nu_diag,1010) ' tr_pond_topo              = ', tr_pond_topo
         write(nu_diag,1010) ' restart_pond_topo         = ', restart_pond_topo
         write(nu_diag,1010) ' tr_aero                   = ', tr_aero
         write(nu_diag,1010) ' restart_aero              = ', restart_aero

         nt_Tsfc = 1           ! index tracers, starting with Tsfc = 1
         ntrcr = 1             ! count tracers, starting with Tsfc = 1

         nt_qice = ntrcr + 1
         ntrcr = ntrcr + nilyr ! qice in nilyr layers

         nt_qsno = ntrcr + 1
         ntrcr = ntrcr + nslyr ! qsno in nslyr layers

         nt_sice = ntrcr + 1
         ntrcr = ntrcr + nilyr ! sice in nilyr layers

         nt_iage = max_ntrcr
         if (tr_iage) then
             ntrcr = ntrcr + 1
             nt_iage = ntrcr   ! chronological ice age
         endif

         nt_FY = max_ntrcr
         if (tr_FY) then
             ntrcr = ntrcr + 1
             nt_FY = ntrcr     ! area of first year ice
         endif

         nt_alvl = max_ntrcr
         nt_vlvl = max_ntrcr
         if (tr_lvl) then
             ntrcr = ntrcr + 1
             nt_alvl = ntrcr
             ntrcr = ntrcr + 1
             nt_vlvl = ntrcr
         endif

         nt_apnd = max_ntrcr
         nt_hpnd = max_ntrcr
         nt_ipnd = max_ntrcr
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

         ! tcraig, tcx, this is a BAD kludge, NTRAERO should be 0 if tr_aero is false
         nt_aero = max_ntrcr - 4*n_aero
         if (tr_aero) then
             nt_aero = ntrcr + 1
             ntrcr = ntrcr + 4*n_aero ! 4 dEdd layers, n_aero species
         endif
              
         if (ntrcr > max_ntrcr-1) then
            if (my_task == master_task) then
               write(nu_diag,*) 'ERROR: max_ntrcr-1 < number of namelist tracers'
               write(nu_diag,*) 'ERROR:   max_ntrcr-1 = ',max_ntrcr-1,' ntrcr = ',ntrcr
            endif
            abort_flag = 19
         endif                               

         write(nu_diag,*) ' '
         write(nu_diag,1020) 'ntrcr = ', ntrcr
         write(nu_diag,*) ' '
         write(nu_diag,1020)'nt_sice = ', nt_sice
         write(nu_diag,1020)'nt_qice = ', nt_qice
         write(nu_diag,1020)'nt_qsno = ', nt_qsno
         write(nu_diag,*)' '
         write(nu_diag,1020)'nilyr', nilyr
         write(nu_diag,1020)'nslyr', nslyr
         write(nu_diag,*)' '

 1000    format (a30,2x,f9.2)  ! a30 to align formatted, unformatted statements
 1005    format (a30,2x,f9.6)  ! float
 1010    format (a30,2x,l6)    ! logical
 1020    format (a30,2x,i6)    ! integer
 1021    format (a30,2x,a8,i6) ! char, int
 1030    format (a30,   a8)    ! character
 1040    format (a30,2x,6i6)   ! integer
 1050    format (a30,2x,6a6)   ! character

         write(nu_diag,*) ' '
         if (grid_type  /=  'displaced_pole' .and. &
             grid_type  /=  'tripole'        .and. &
             grid_type  /=  'column'         .and. &
             grid_type  /=  'rectangular'    .and. &
             grid_type  /=  'cpom_grid'      .and. &
             grid_type  /=  'regional'       .and. &
             grid_type  /=  'latlon' ) then 
            if (my_task == master_task) write(nu_diag,*)'ERROR: unknown grid_type=',trim(grid_type)
            abort_flag = 20
         endif

      endif                     ! my_task = master_task

      call broadcast_scalar(ntrcr,    master_task)
      call broadcast_scalar(nt_Tsfc,  master_task)
      call broadcast_scalar(nt_sice,  master_task)
      call broadcast_scalar(nt_qice,  master_task)
      call broadcast_scalar(nt_qsno,  master_task)
      call broadcast_scalar(nt_iage,  master_task)
      call broadcast_scalar(nt_FY,    master_task)
      call broadcast_scalar(nt_alvl,  master_task)
      call broadcast_scalar(nt_vlvl,  master_task)
      call broadcast_scalar(nt_apnd,  master_task)
      call broadcast_scalar(nt_hpnd,  master_task)
      call broadcast_scalar(nt_ipnd,  master_task)
      call broadcast_scalar(nt_aero,  master_task)

      if (formdrag) then
         if (nt_apnd==0) then
            if (my_task == master_task) write(nu_diag,*)'ERROR: formdrag=T, nt_apnd=',nt_apnd
            abort_flag = 21
         elseif (nt_hpnd==0) then
            if (my_task == master_task) write(nu_diag,*)'ERROR: formdrag=T, nt_hpnd=',nt_hpnd
            abort_flag = 22
         elseif (nt_ipnd==0) then
            if (my_task == master_task) write(nu_diag,*)'ERROR: formdrag=T, nt_ipnd=',nt_ipnd
            abort_flag = 23
         elseif (nt_alvl==0) then
            if (my_task == master_task) write(nu_diag,*)'ERROR: formdrag=T, nt_alvl=',nt_alvl
            abort_flag = 24
         elseif (nt_vlvl==0) then
            if (my_task == master_task) write(nu_diag,*)'ERROR: formdrag=T, nt_vlvl=',nt_vlvl
            abort_flag = 25
         endif
      endif

      call flush_fileunit(nu_diag)
      call icepack_init_parameters(ustar_min_in=ustar_min, albicev_in=albicev, albicei_in=albicei, &
         albsnowv_in=albsnowv, albsnowi_in=albsnowi, natmiter_in=natmiter, emissivity_in=emissivity, &
         ahmax_in=ahmax, shortwave_in=shortwave, albedo_type_in=albedo_type, R_ice_in=R_ice, R_pnd_in=R_pnd, &
         R_snw_in=R_snw, dT_mlt_in=dT_mlt, rsnw_mlt_in=rsnw_mlt, &
         kstrength_in=kstrength, krdg_partic_in=krdg_partic, krdg_redist_in=krdg_redist, mu_rdg_in=mu_rdg, &
         atmbndy_in=atmbndy, calc_strair_in=calc_strair, formdrag_in=formdrag, highfreq_in=highfreq, &
         kitd_in=kitd, kcatbound_in=kcatbound, hs0_in=hs0, dpscale_in=dpscale, frzpnd_in=frzpnd, &
         rfracmin_in=rfracmin, rfracmax_in=rfracmax, pndaspect_in=pndaspect, hs1_in=hs1, hp1_in=hp1, &
         ktherm_in=ktherm, calc_Tsfc_in=calc_Tsfc, conduct_in=conduct, &
         a_rapid_mode_in=a_rapid_mode, Rac_rapid_mode_in=Rac_rapid_mode, &
         aspect_rapid_mode_in=aspect_rapid_mode, dSdt_slow_mode_in=dSdt_slow_mode, &
         phi_c_slow_mode_in=phi_c_slow_mode, phi_i_mushy_in=phi_i_mushy, &
         tfrz_option_in=tfrz_option, kalg_in=kalg, fbot_xfer_type_in=fbot_xfer_type)
      call icepack_init_tracer_numbers(ntrcr_in=ntrcr)
      call icepack_init_tracer_flags(tr_iage_in=tr_iage, tr_FY_in=tr_FY, &
         tr_lvl_in=tr_lvl, tr_aero_in=tr_aero, tr_pond_in=tr_pond, &
         tr_pond_cesm_in=tr_pond_cesm, tr_pond_lvl_in=tr_pond_lvl, tr_pond_topo_in=tr_pond_topo)
      call icepack_init_tracer_indices(nt_Tsfc_in=nt_Tsfc, nt_sice_in=nt_sice, &
         nt_qice_in=nt_qice, nt_qsno_in=nt_qsno, nt_iage_in=nt_iage, nt_fy_in=nt_fy, &
         nt_alvl_in=nt_alvl, nt_vlvl_in=nt_vlvl, nt_apnd_in=nt_apnd, nt_hpnd_in=nt_hpnd, &
         nt_ipnd_in=nt_ipnd, nt_aero_in=nt_aero)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname//' Icepack Abort2', &
         file=__FILE__, line=__LINE__)

      call flush_fileunit(nu_diag)
      call ice_barrier()
      if (abort_flag /= 0) then
         write(nu_diag,*) subname,' ERROR: abort_flag=',abort_flag
         call abort_ice (subname//' ABORTING on input ERRORS', &
            file=__FILE__, line=__LINE__)
      endif

      end subroutine input_data

!=======================================================================

! Initialize state for the itd model
!
! authors: C. M. Bitz, UW
!          William H. Lipscomb, LANL

      subroutine init_state

      use ice_blocks, only: block, get_block, nx_block, ny_block
      use ice_domain, only: nblocks, blocks_ice
      use ice_domain_size, only: ncat, nilyr, nslyr, max_ntrcr, n_aero
      use ice_flux, only: sst, Tf, Tair, salinz, Tmltz
      use ice_grid, only: tmask, ULON, ULAT, TLON, TLAT
      use ice_state, only: trcr_depend, aicen, trcrn, vicen, vsnon, &
          aice0, aice, vice, vsno, trcr, aice_init, bound_state, &
          n_trcr_strata, nt_strata, trcr_base

      integer (kind=int_kind) :: &
         ilo, ihi    , & ! physical domain indices
         jlo, jhi    , & ! physical domain indices
         iglob(nx_block), & ! global indices
         jglob(ny_block), & ! global indices
         i, j        , & ! horizontal indices
         k           , & ! vertical index
         it          , & ! tracer index
         iblk            ! block index

      logical (kind=log_kind) :: &
         heat_capacity   ! from icepack

      integer (kind=int_kind) :: ntrcr
      logical (kind=log_kind) :: tr_iage, tr_FY, tr_lvl, tr_aero
      logical (kind=log_kind) :: tr_pond_cesm, tr_pond_lvl, tr_pond_topo
      integer (kind=int_kind) :: nt_Tsfc, nt_sice, nt_qice, nt_qsno, nt_iage, nt_fy
      integer (kind=int_kind) :: nt_alvl, nt_vlvl, nt_apnd, nt_hpnd, nt_ipnd, nt_aero

      type (block) :: &
         this_block           ! block information for current block

      character(len=*), parameter :: subname='(init_state)'

      call icepack_query_parameters(heat_capacity_out=heat_capacity)
      call icepack_query_tracer_numbers(ntrcr_out=ntrcr)
      call icepack_query_tracer_flags(tr_iage_out=tr_iage, tr_FY_out=tr_FY, &
        tr_lvl_out=tr_lvl, tr_aero_out=tr_aero, &
        tr_pond_cesm_out=tr_pond_cesm, tr_pond_lvl_out=tr_pond_lvl, tr_pond_topo_out=tr_pond_topo)
      call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_sice_out=nt_sice, &
        nt_qice_out=nt_qice, nt_qsno_out=nt_qsno, nt_iage_out=nt_iage, nt_fy_out=nt_fy, &
        nt_alvl_out=nt_alvl, nt_vlvl_out=nt_vlvl, nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, &
        nt_ipnd_out=nt_ipnd, nt_aero_out=nt_aero)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      !-----------------------------------------------------------------
      ! Check number of layers in ice and snow.
      !-----------------------------------------------------------------

      if (my_task == master_task) then
 
         if (nilyr < 1) then
            write(nu_diag,*) 'ERROR: Must have at least one ice layer'
            write(nu_diag,*) 'ERROR:   nilyr =', nilyr
            call abort_ice (error_message=subname//' Not enough ice layers', &
               file=__FILE__, line=__LINE__)
         endif

         if (nslyr < 1) then
            write(nu_diag,*) 'ERROR: Must have at least one snow layer'
            write(nu_diag,*) 'ERROR:   nslyr =', nslyr
            call abort_ice(error_message=subname//' Not enough snow layers', &
               file=__FILE__, line=__LINE__)
         endif

         if (.not.heat_capacity) then

            if (nilyr > 1) then
               write(nu_diag,*) 'ERROR: Must have nilyr = 1 if heat_capacity=F'
               write(nu_diag,*) 'ERROR:   nilyr =', nilyr
               call abort_ice(error_message=subname//' Too many ice layers', &
                  file=__FILE__, line=__LINE__)
            endif

            if (nslyr > 1) then
               write(nu_diag,*) 'ERROR: Must have nslyr = 1 if heat_capacity=F'
               write(nu_diag,*) 'ERROR:  nslyr =', nslyr
               call abort_ice(error_message=subname//' Too many snow layers', &
                  file=__FILE__, line=__LINE__)
            endif

         endif   ! heat_capacity = F

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
      if (tr_pond_cesm) then
                   trcr_depend(nt_apnd)  = 0           ! melt pond area
                   trcr_depend(nt_hpnd)  = 2+nt_apnd   ! melt pond depth
      endif
      if (tr_pond_lvl) then
                   trcr_depend(nt_apnd)  = 2+nt_alvl   ! melt pond area
                   trcr_depend(nt_hpnd)  = 2+nt_apnd   ! melt pond depth
                   trcr_depend(nt_ipnd)  = 2+nt_apnd   ! refrozen pond lid
      endif
      if (tr_pond_topo) then
                   trcr_depend(nt_apnd)  = 0           ! melt pond area
                   trcr_depend(nt_hpnd)  = 2+nt_apnd   ! melt pond depth
                   trcr_depend(nt_ipnd)  = 2+nt_apnd   ! refrozen pond lid
      endif
      if (tr_aero) then ! volume-weighted aerosols
         do it = 1, n_aero
            trcr_depend(nt_aero+(it-1)*4  ) = 2 ! snow
            trcr_depend(nt_aero+(it-1)*4+1) = 2 ! snow
            trcr_depend(nt_aero+(it-1)*4+2) = 1 ! ice
            trcr_depend(nt_aero+(it-1)*4+3) = 1 ! ice
         enddo
      endif

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

      if (tr_pond_cesm) then
         n_trcr_strata(nt_hpnd)   = 1       ! melt pond depth
         nt_strata    (nt_hpnd,1) = nt_apnd ! on melt pond area
      endif
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
      if (tr_pond_topo) then
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
                             ULON (:,:,    iblk), ULAT (:,:,    iblk), &
                             TLON (:,:,    iblk), TLAT (:,:,    iblk), &
                             Tair (:,:,    iblk), sst  (:,:,    iblk), &
                             Tf   (:,:,    iblk),                      &
                             salinz(:,:,:, iblk), Tmltz(:,:,:,  iblk), &
                             aicen(:,:,  :,iblk), trcrn(:,:,:,:,iblk), &
                             vicen(:,:,  :,iblk), vsnon(:,:,  :,iblk))

      enddo                     ! iblk
      !$OMP END PARALLEL DO

      !-----------------------------------------------------------------
      ! ghost cell updates
      !-----------------------------------------------------------------

      call bound_state (aicen,        &
                        vicen, vsnon, &
                        ntrcr, trcrn)

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
         do it = 1, max_ntrcr
            trcr(i,j,it,iblk) = c0
         enddo

         if (tmask(i,j,iblk)) &
         call icepack_aggregate (ncat,               &
                                aicen(i,j,:,iblk),  &
                                trcrn(i,j,1:ntrcr,:,iblk), &
                                vicen(i,j,:,iblk),   &
                                vsnon(i,j,:,iblk),   &
                                aice (i,j,  iblk),   &
                                trcr (i,j,1:ntrcr,iblk),   &
                                vice (i,j,  iblk),   &
                                vsno (i,j,  iblk),   &
                                aice0(i,j,  iblk),   &
                                ntrcr,               &
                                trcr_depend  (1:ntrcr),&
                                trcr_base    (1:ntrcr,:),&
                                n_trcr_strata(1:ntrcr),&
                                nt_strata    (1:ntrcr,:))

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
                                ULON,     ULAT, &
                                TLON,     TLAT, &
                                Tair,     sst,  &
                                Tf,       &
                                salinz,   Tmltz, &
                                aicen,    trcrn, &
                                vicen,    vsnon)

      use ice_arrays_column, only: hin_max
      use ice_domain_size, only: nilyr, nslyr, nx_global, ny_global, max_ntrcr, ncat
      use ice_grid, only: grid_type
      use ice_forcing, only: atm_data_type

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ilo, ihi          , & ! physical domain indices
         jlo, jhi          , & !
         iglob(nx_block)   , & ! global indices
         jglob(ny_block)       !

      character(len=char_len_long), intent(in) :: & 
         ice_ic      ! method of ice cover initialization

      logical (kind=log_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         tmask      ! true for ice/ocean cells

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         ULON   , & ! longitude of velocity pts (radians)
         ULAT   , & ! latitude of velocity pts (radians)
         TLON   , & ! longitude of temperature pts (radians)
         TLAT       ! latitude of temperature pts (radians)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         Tair    , & ! air temperature  (K)
         Tf      , & ! freezing temperature (C) 
         sst         ! sea surface temperature (C) 

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr), &
         intent(in) :: &
         salinz  , & ! initial salinity profile
         Tmltz       ! initial melting temperature profile

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), &
         intent(out) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr,ncat), &
         intent(out) :: &
         trcrn     ! ice tracers
                   ! 1: surface temperature of ice/snow (C)

      ! local variables

      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij          , & ! horizontal index, combines i and j loops
         k           , & ! ice layer index
         n           , & ! thickness category index
         it          , & ! tracer index
         icells          ! number of cells initialized with ice

      integer (kind=int_kind), dimension(nx_block*ny_block) :: &
         indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind) :: &
         Tsfc, sum, hbar, puny, rhos, Lfresh, rad_to_deg

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

      logical (kind=log_kind) :: tr_brine, tr_lvl
      integer (kind=int_kind) :: nt_Tsfc, nt_qice, nt_qsno, nt_sice
      integer (kind=int_kind) :: nt_fbri, nt_alvl, nt_vlvl

      character(len=*), parameter :: subname='(set_state_var)'

      call icepack_query_tracer_flags(tr_brine_out=tr_brine, tr_lvl_out=tr_lvl)
      call icepack_query_tracer_indices( nt_Tsfc_out=nt_Tsfc, nt_qice_out=nt_qice, &
        nt_qsno_out=nt_qsno, nt_sice_out=nt_sice, &
        nt_fbri_out=nt_fbri, nt_alvl_out=nt_alvl, nt_vlvl_out=nt_vlvl)
      call icepack_query_parameters(rhos_out=rhos, Lfresh_out=Lfresh, puny_out=puny, &
        rad_to_deg_out=rad_to_deg)
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
            trcrn(i,j,nt_Tsfc,n) = Tf(i,j)  ! surface temperature 
            if (max_ntrcr >= 2) then
               do it = 2, max_ntrcr
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
         enddo
         enddo
      enddo

      if (trim(ice_ic) == 'default') then

      !-----------------------------------------------------------------
      ! Place ice where ocean surface is cold.
      ! Note: If SST is not read from a file, then the ocean is assumed
      !       to be at its freezing point everywhere, and ice will
      !       extend to the prescribed edges.
      !-----------------------------------------------------------------

         if (trim(atm_data_type) == 'box') then

            hbar = c2  ! initial ice thickness
            do n = 1, ncat
               hinit(n) = c0
               ainit(n) = c0
               if (hbar > hin_max(n-1) .and. hbar < hin_max(n)) then
                  hinit(n) = hbar
                  ainit(n) = 0.50 !echmod symm
               endif
            enddo

         else

      ! initial category areas in cells with ice
         hbar = c3  ! initial ice thickness with greatest area
                    ! Note: the resulting average ice thickness 
                    ! tends to be less than hbar due to the
                    ! nonlinear distribution of ice thicknesses 
         sum = c0
         do n = 1, ncat
            if (n < ncat) then
               hinit(n) = p5*(hin_max(n-1) + hin_max(n)) ! m
            else                ! n=ncat
               hinit(n) = (hin_max(n-1) + c1) ! m
            endif
            ! parabola, max at h=hbar, zero at h=0, 2*hbar
            ainit(n) = max(c0, (c2*hbar*hinit(n) - hinit(n)**2))
            sum = sum + ainit(n)
         enddo
         do n = 1, ncat
            ainit(n) = ainit(n) / (sum + puny/ncat) ! normalize
         enddo

         endif ! atm_data_type

         if (trim(grid_type) == 'rectangular') then

         ! place ice on left side of domain
         icells = 0
         do j = jlo, jhi
         do i = ilo, ihi
            if (tmask(i,j)) then
               if (ULON(i,j) < -50./rad_to_deg) then
                  icells = icells + 1
                  indxi(icells) = i
                  indxj(icells) = j
               endif            ! ULON
            endif               ! tmask
         enddo                  ! i
         enddo                  ! j

         else

         ! place ice at high latitudes where ocean sfc is cold
         icells = 0
         do j = jlo, jhi
         do i = ilo, ihi
            if (tmask(i,j)) then
               ! place ice in high latitudes where ocean sfc is cold
               if ( (sst (i,j) <= Tf(i,j)+p2) .and. &
                    (TLAT(i,j) < edge_init_sh/rad_to_deg .or. &
                     TLAT(i,j) > edge_init_nh/rad_to_deg) ) then
                  icells = icells + 1
                  indxi(icells) = i
                  indxj(icells) = j
               endif            ! cold surface
            endif               ! tmask
         enddo                  ! i
         enddo                  ! j

         endif                  ! rectgrid

         do n = 1, ncat

            ! ice volume, snow volume
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               aicen(i,j,n) = ainit(n)

               if (trim(atm_data_type) == 'box') then
                  if (hinit(n) > c0) then
!                  ! constant slope from 0 to 1 in x direction
!                     aicen(i,j,n) = (real(iglob(i), kind=dbl_kind)-p5) &
!                                  / (real(nx_global,kind=dbl_kind))
!                  ! constant slope from 0 to 0.5 in x direction
!                     aicen(i,j,n) = (real(iglob(i), kind=dbl_kind)-p5) &
!                                  / (real(nx_global,kind=dbl_kind)) * p5
                  ! quadratic
!                     aicen(i,j,n) = max(c0,(real(iglob(i), kind=dbl_kind)-p5) &
!                                         / (real(nx_global,kind=dbl_kind)) &
!                                         * (real(jglob(j), kind=dbl_kind)-p5) &
!                                         / (real(ny_global,kind=dbl_kind)) * p5)
                     aicen(i,j,n) = max(c0,(real(nx_global, kind=dbl_kind) &
                                         -  real(iglob(i), kind=dbl_kind)-p5) &
                                         / (real(nx_global,kind=dbl_kind)) &
                                         * (real(ny_global, kind=dbl_kind) &
                                         -  real(jglob(j), kind=dbl_kind)-p5) &
                                         / (real(ny_global,kind=dbl_kind)) * p5)
                  endif
                  vicen(i,j,n) = hinit(n) * aicen(i,j,n) ! m
               else
                  vicen(i,j,n) = hinit(n) * ainit(n) ! m
               endif
               vsnon(i,j,n) = min(aicen(i,j,n)*hsno_init,p2*vicen(i,j,n))

               call icepack_init_trcr(Tair(i,j),     Tf(i,j),      &
                                     salinz(i,j,:), Tmltz(i,j,:), &
                                     Tsfc,                        &
                                     nilyr,         nslyr,        &
                                     qin(:),        qsn(:))

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
      endif                     ! ice_ic

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      end subroutine set_state_var

!=======================================================================

      end module ice_init

!=======================================================================
