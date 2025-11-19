module ice_comp_nuopc

  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for CICE
  !----------------------------------------------------------------------------

  use ESMF
  use NUOPC              , only : NUOPC_CompDerive, NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize
  use NUOPC              , only : NUOPC_CompFilterPhaseMap, NUOPC_IsUpdated, NUOPC_IsAtTime
  use NUOPC              , only : NUOPC_CompAttributeGet, NUOPC_Advertise
  use NUOPC              , only : NUOPC_SetAttribute, NUOPC_CompAttributeGet, NUOPC_CompAttributeSet
  use NUOPC_Model        , only : model_routine_SS           => SetServices
  use NUOPC_Model        , only : model_label_Advance        => label_Advance
  use NUOPC_Model        , only : model_label_DataInitialize => label_DataInitialize
  use NUOPC_Model        , only : model_label_SetRunClock    => label_SetRunClock
  use NUOPC_Model        , only : model_label_Finalize       => label_Finalize
  use NUOPC_Model        , only : NUOPC_ModelGet, SetVM
  use ice_constants      , only : ice_init_constants, c0
  use ice_shr_methods    , only : chkerr, state_setscalar, state_getscalar, state_diagnose, alarmInit
  use ice_shr_methods    , only : get_component_instance, state_flddebug

  use ice_import_export  , only : ice_import, ice_export, ice_advertise_fields, ice_realize_fields
  use ice_domain_size    , only : nx_global, ny_global
  use ice_grid           , only : grid_format, init_grid2
  use ice_communicate    , only : init_communicate, my_task, master_task, mpi_comm_ice
  use ice_calendar       , only : force_restart_now, write_ic
  use ice_calendar       , only : idate, idate0,  mday, mmonth, myear, year_init, month_init, day_init
  use ice_calendar       , only : msec, dt, calendar, calendar_type, nextsw_cday, istep, use_leap_years
  use ice_calendar       , only : ice_calendar_noleap, ice_calendar_proleptic_gregorian, ice_calendar_gregorian
  use ice_kinds_mod      , only : dbl_kind, int_kind, char_len, char_len_long
  use ice_fileunits      , only : nu_diag, nu_diag_set, inst_index, inst_name
  use ice_fileunits      , only : inst_suffix, release_all_fileunits, flush_fileunit
  use ice_restart_shared , only : runid, runtype, restart, use_restart_time, restart_dir, restart_file, &
                                  restart_format, restart_chunksize, pointer_date
  use ice_history        , only : accum_hist
  use ice_history_shared , only : history_format, history_chunksize
  use ice_exit           , only : abort_ice
  use icepack_intfc      , only : icepack_warnings_flush, icepack_warnings_aborted
  use icepack_intfc      , only : icepack_init_orbit, icepack_init_parameters, icepack_query_orbit
  use icepack_intfc      , only : icepack_query_tracer_flags, icepack_query_parameters
  use cice_wrapper_mod   , only : t_startf, t_stopf, t_barrierf
  use cice_wrapper_mod   , only : shr_log_getlogunit, shr_log_setlogunit
  use cice_wrapper_mod   , only : ufs_settimer, ufs_logtimer, ufs_file_setlogunit, wtime
#ifdef CESMCOUPLED
  use shr_const_mod
  use shr_orb_mod        , only : shr_orb_decl, shr_orb_params, SHR_ORB_UNDEF_REAL, SHR_ORB_UNDEF_INT
  use ice_scam           , only : scmlat, scmlon, scol_mask, scol_frac, scol_ni, scol_nj, scol_area
  use nuopc_shr_methods  , only : set_component_logging
#else
  use ice_shr_methods    , only : set_component_logging
#endif
  use ice_timers
  use CICE_InitMod       , only : cice_init1, cice_init2
  use CICE_RunMod        , only : cice_run
  use ice_mesh_mod       , only : ice_mesh_set_distgrid, ice_mesh_setmask_from_maskfile, ice_mesh_check
  use ice_mesh_mod       , only : ice_mesh_init_tlon_tlat_area_hm, ice_mesh_create_scolumn
  use ice_prescribed_mod , only : ice_prescribed_init
  use ice_scam           , only : scol_valid, single_column
#ifndef CESMCOUPLED
  use shr_is_restart_fh_mod, only : init_is_restart_fh, is_restart_fh, is_restart_fh_type
#endif

  implicit none
  private

  public  :: SetServices
  public  :: SetVM

  private :: InitializeP0
  private :: InitializeAdvertise
  private :: InitializeRealize
  private :: ModelAdvance
  private :: ModelSetRunClock
  private :: ModelFinalize
  private :: ice_orbital_init ! only valid for cesm

  character(len=char_len_long) :: flds_scalar_name = ''
  integer                      :: flds_scalar_num = 0
  integer                      :: flds_scalar_index_nx = 0
  integer                      :: flds_scalar_index_ny = 0
  integer                      :: flds_scalar_index_nextsw_cday = 0

  character(len=char_len_long) :: orb_mode        ! attribute - orbital mode
  integer                      :: orb_iyear       ! attribute - orbital year
  integer                      :: orb_iyear_align ! attribute - associated with model year
  real(dbl_kind)               :: orb_obliq       ! attribute - obliquity in degrees
  real(dbl_kind)               :: orb_mvelp       ! attribute - moving vernal equinox longitude
  real(dbl_kind)               :: orb_eccen       ! attribute and update-  orbital eccentricity

  character(len=*) , parameter :: orb_fixed_year       = 'fixed_year'
  character(len=*) , parameter :: orb_variable_year    = 'variable_year'
  character(len=*) , parameter :: orb_fixed_parameters = 'fixed_parameters'

  type(ESMF_Mesh)              :: ice_mesh

  integer                      :: nthrds       ! Number of threads to use in this component
  integer                      :: nu_timer = 6 ! Simple timer log, unused except by UFS
  integer                      :: dbug = 0
  logical                      :: profile_memory = .false.
  logical                      :: mastertask
  logical                      :: runtimelog = .false.
  logical                      :: restart_eor = .false. !End of run restart flag
#ifndef CESMCOUPLED
  type(is_restart_fh_type)     :: restartfh_info     ! For flexible restarts in UFS
#endif
  integer                      :: start_ymd          ! Start date (YYYYMMDD)
  integer                      :: start_tod          ! start time of day (s)
  integer                      :: curr_ymd           ! Current date (YYYYMMDD)
  integer                      :: curr_tod           ! Current time of day (s)
  integer                      :: stop_ymd           ! stop date (YYYYMMDD)
  integer                      :: stop_tod           ! stop time of day (sec)
  integer                      :: ref_ymd            ! Reference date (YYYYMMDD)
  integer                      :: ref_tod            ! reference time of day (s)
  integer     , parameter      :: debug_import = 0 ! internal debug level
  integer     , parameter      :: debug_export = 0 ! internal debug level
  character(*), parameter      :: modName =  "(ice_comp_nuopc)"
  character(*), parameter      :: u_FILE_u = &
       __FILE__

!=======================================================================
contains
!===============================================================================

  subroutine SetServices(gcomp, rc)

    ! Arguments
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Local variables
    character(len=*),parameter  :: subname=trim(modName)//':(SetServices) '
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! the NUOPC gcomp component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! switching to IPD versions
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         userRoutine=InitializeP0, phase=0, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p3"/), userRoutine=InitializeRealize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! attach specializing method(s)
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
         specRoutine=ModelAdvance, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MethodRemove(gcomp, label=model_label_SetRunClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_SetRunClock, &
         specRoutine=ModelSetRunClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
         specRoutine=ModelFinalize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine SetServices

  !===============================================================================

  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)

    ! Arguments
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc

    logical                      :: isPresent, isSet
    character(len=64)            :: value
    character(len=char_len_long) :: logmsg
    !--------------------------------

    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, &
         acceptStringList=(/"IPDv01p"/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    profile_memory = .false.
    call NUOPC_CompAttributeGet(gcomp, name="ProfileMemory", value=value, &
         isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) profile_memory=(trim(value)=="true")
    write(logmsg,*) profile_memory
    call ESMF_LogWrite('CICE_cap:ProfileMemory = '//trim(logmsg), ESMF_LOGMSG_INFO)

  end subroutine InitializeP0

  !===============================================================================

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)

    ! Arguments
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! Local variables
    character(len=char_len_long) :: cvalue
    character(len=char_len_long) :: ice_meshfile
    character(len=char_len_long) :: ice_maskfile
    character(len=char_len_long) :: errmsg
    logical                      :: isPresent, isSet
    real(dbl_kind)               :: eccen, obliqr, lambm0, mvelpp
    type(ESMF_DistGrid)          :: ice_distGrid
    real(kind=dbl_kind)          :: atmiter_conv
    real(kind=dbl_kind)          :: atmiter_conv_driver
    integer (kind=int_kind)      :: natmiter
    integer (kind=int_kind)      :: natmiter_driver
    integer                      :: localPet
    integer                      :: npes
    type(ESMF_VM)                :: vm
    integer                      :: lmpicom            ! local communicator
    type(ESMF_Time)              :: currTime           ! Current time
    type(ESMF_Time)              :: startTime          ! Start time
    type(ESMF_Time)              :: stopTime           ! Stop time
    type(ESMF_Time)              :: refTime            ! Ref time
    type(ESMF_TimeInterval)      :: timeStep           ! Model timestep
    type(ESMF_CalKind_Flag)      :: esmf_caltype       ! esmf calendar type
    integer                      :: yy,mm,dd           ! Temporaries for time query
    integer                      :: dtime              ! time step
    integer                      :: shrlogunit         ! original log unit
    character(len=char_len)      :: starttype          ! infodata start type
    integer                      :: lsize              ! local size of coupling array
    integer                      :: n,c,g,i,j,m        ! indices
    integer                      :: iblk, jblk         ! indices
    integer                      :: ig, jg             ! indices
    integer                      :: ilo, ihi, jlo, jhi ! beginning and end of physical domain

    character(len=char_len_long) :: diag_filename = 'unset'
    character(len=char_len_long) :: logmsg
    character(len=char_len_long) :: single_column_lnd_domainfile
    real(dbl_kind)               :: scol_lon
    real(dbl_kind)               :: scol_lat
    real(dbl_kind)               :: scol_spval
    character(len=char_len)      :: tfrz_option    ! tfrz_option from cice namelist
    character(len=char_len)      :: tfrz_option_driver    ! tfrz_option from cice namelist
    character(len=*), parameter :: subname=trim(modName)//':(InitializeAdvertise) '
    !--------------------------------

    call ufs_settimer(wtime)

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       flds_scalar_name = trim(cvalue)
       call ESMF_LogWrite(trim(subname)//' flds_scalar_name = '//trim(flds_scalar_name), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call abort_ice(subname//'Need to set attribute ScalarFieldName')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldCount", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue, *) flds_scalar_num
       write(logmsg,*) flds_scalar_num
       call ESMF_LogWrite(trim(subname)//' flds_scalar_num = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call abort_ice(subname//'Need to set attribute ScalarFieldCount')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNX", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_nx
       write(logmsg,*) flds_scalar_index_nx
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_nx = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call abort_ice(subname//'Need to set attribute ScalarFieldIdxGridNX')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNY", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_ny
       write(logmsg,*) flds_scalar_index_ny
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_ny = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call abort_ice(subname//'Need to set attribute ScalarFieldIdxGridNY')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxNextSwCday", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_nextsw_cday
       write(logmsg,*) flds_scalar_index_nextsw_cday
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_nextsw_cday = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    call NUOPC_CompAttributeGet(gcomp, name='dbug_flag', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
     read(cvalue,*) dbug
    end if
    write(logmsg,'(i6)') dbug
    call ESMF_LogWrite('CICE_cap: dbug = '//trim(logmsg), ESMF_LOGMSG_INFO)

    call NUOPC_CompAttributeGet(gcomp, name="RunTimeLog", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) runtimelog=(trim(cvalue)=="true")
    write(logmsg,*) runtimelog
    call ESMF_LogWrite('CICE_cap:RunTimeLog = '//trim(logmsg), ESMF_LOGMSG_INFO)

    call NUOPC_CompAttributeGet(gcomp, name="write_restart_at_endofrun", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       if (trim(cvalue) .eq. '.true.') restart_eor = .true.
    endif

#ifdef CESMCOUPLED
    pointer_date = .true.
#endif

    ! set CICE internal pointer_date variable based on nuopc settings
    ! this appends a datestamp to the "rpointer" file
    call NUOPC_CompAttributeGet(gcomp, name="restart_pointer_append_date", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) pointer_date = (trim(cvalue) .eq. ".true.")
    !----------------------------------------------------------------------------
    ! generate local mpi comm
    !----------------------------------------------------------------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, mpiCommunicator=lmpicom, localPet=localPet, PetCount=npes, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

#ifdef CESMCOUPLED
    call ESMF_VMGet(vm, pet=localPet, peCount=nthrds, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (nthrds==1) then
       call NUOPC_CompAttributeGet(gcomp, "nthreads", value=cvalue, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
       read(cvalue,*) nthrds
    endif
!$  call omp_set_num_threads(nthrds)
#endif

    !----------------------------------------------------------------------------
    ! Initialize cice communicators
    !----------------------------------------------------------------------------

    call init_communicate(lmpicom)     ! initial setup for message passing
    mastertask = .false.
    if (my_task == master_task) mastertask = .true.

    !----------------------------------------------------------------------------
    ! determine instance information
    !----------------------------------------------------------------------------

    call get_component_instance(gcomp, inst_suffix, inst_index, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

!   inst_name = "ICE"//trim(inst_suffix)
    inst_name = "ICE"

    !----------------------------------------------------------------------------
    ! start cice timers
    !----------------------------------------------------------------------------

    call t_startf ('cice_init_total')

    !----------------------------------------------------------------------------
    ! Initialize constants
    !----------------------------------------------------------------------------

#ifdef CESMCOUPLED
    call ice_init_constants(omega_in=SHR_CONST_OMEGA, radius_in=SHR_CONST_REARTH, &
       spval_dbl_in=SHR_CONST_SPVAL)

    ! TODO: get tfrz_option from driver

    call icepack_init_parameters( &
       secday_in           = SHR_CONST_CDAY,                  &
       rhoi_in             = SHR_CONST_RHOICE,                &
       rhow_in             = SHR_CONST_RHOSW,                 &
       cp_air_in           = SHR_CONST_CPDAIR,                &
       cp_ice_in           = SHR_CONST_CPICE,                 &
       cp_ocn_in           = SHR_CONST_CPSW,                  &
       gravit_in           = SHR_CONST_G,                     &
       rhofresh_in         = SHR_CONST_RHOFW,                 &
       zvir_in             = SHR_CONST_ZVIR,                  &
       vonkar_in           = SHR_CONST_KARMAN,                &
       cp_wv_in            = SHR_CONST_CPWV,                  &
       stefan_boltzmann_in = SHR_CONST_STEBOL,                &
       Tffresh_in          = SHR_CONST_TKFRZ,                 &
       Lsub_in             = SHR_CONST_LATSUB,                &
       Lvap_in             = SHR_CONST_LATVAP,                &
      !Lfresh_in           = SHR_CONST_LATICE,                & ! computed in init_parameters as Lsub-Lvap
       Timelt_in           = SHR_CONST_TKFRZ-SHR_CONST_TKFRZ, &
       Tsmelt_in           = SHR_CONST_TKFRZ-SHR_CONST_TKFRZ, &
       ice_ref_salinity_in = SHR_CONST_ICE_REF_SAL,           &
       depressT_in         = 0.054_dbl_kind,                  &
       Tocnfrz_in          = -34.0_dbl_kind*0.054_dbl_kind,   &
       pi_in               = SHR_CONST_PI,                    &
       snowpatch_in        = 0.005_dbl_kind)

    call icepack_warnings_flush(nu_diag)
    if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
        file=__FILE__, line=__LINE__)
#endif

    !----------------------------------------------------------------------------
    ! Determine attributes - also needed in realize phase to get grid information
    !----------------------------------------------------------------------------


    ! Get orbital values
    ! Note that these values are obtained in a call to init_orbit in ice_shortwave.F90
    ! if CESMCOUPLED is not defined

    call ice_orbital_init(gcomp, clock, nu_diag, mastertask, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine runtype and possibly nextsw_cday
    call NUOPC_CompAttributeGet(gcomp, name='start_type', value=cvalue, isPresent=isPresent, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       read(cvalue,*) starttype
       if (trim(starttype) == trim('startup')) then
          runtype = "initial"
       else if (trim(starttype) == trim('continue') ) then
          runtype = "continue"
          restart = .true.
          use_restart_time = .true.
       else if (trim(starttype) == trim('branch')) then
          runtype = "continue"
          restart = .true.
          use_restart_time = .true.
       else
          call abort_ice( subname//' ERROR: unknown starttype' )
       end if

       ! We assume here that on startup - nextsw_cday is just the current time
       ! TOOD (mvertens, 2019-03-21): need to get the perpetual run working
       if (trim(runtype) /= 'initial') then
          ! Set nextsw_cday to -1 (this will skip an orbital calculation on initialization
          nextsw_cday = -1.0_dbl_kind
       else
          call ESMF_ClockGet( clock, currTime=currTime, rc=rc )
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_TimeGet( currTime, dayOfYear_r8=nextsw_cday, rc=rc )
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    else
       runtype = 'initial' ! determined from the namelist in ice_init if CESMCOUPLED is not defined
    end if

    ! Determine runid
    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (isPresent .and. isSet) then
       read(cvalue,*) runid
    else
       ! read in from the namelist in ice_init.F90 if this is not an attribute passed from the driver
       runid = 'unknown'
    end if

    ! Get clock information before call to cice_init
    call ESMF_ClockGet( clock, &
         currTime=currTime, startTime=startTime, stopTime=stopTime, refTime=RefTime, &
         timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet( currTime, yy=yy, mm=mm, dd=dd, s=curr_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ice_cal_ymd2date(yy,mm,dd,curr_ymd)

    call ESMF_TimeGet( startTime, yy=yy, mm=mm, dd=dd, s=start_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ice_cal_ymd2date(yy,mm,dd,start_ymd)

    call ESMF_TimeGet( stopTime, yy=yy, mm=mm, dd=dd, s=stop_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ice_cal_ymd2date(yy,mm,dd,stop_ymd)

    call ESMF_TimeGet( refTime, yy=yy, mm=mm, dd=dd, s=ref_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ice_cal_ymd2date(yy,mm,dd,ref_ymd)

    call ESMF_TimeIntervalGet( timeStep, s=dtime, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    dt = real(dtime)

    call ESMF_TimeGet( currTime, calkindflag=esmf_caltype, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (esmf_caltype == ESMF_CALKIND_NOLEAP) then
       calendar_type = ice_calendar_noleap
    else if (esmf_caltype == ESMF_CALKIND_GREGORIAN) then
       calendar_type = ice_calendar_proleptic_gregorian
    else
       call abort_ice( subname//'ERROR:: bad calendar for ESMF' )
    end if

    !----------------------------------------------------------------------------
    ! Set cice logging
    !----------------------------------------------------------------------------
    ! Note - this must be done AFTER the communicators are set
    ! Note that sets the nu_diag module variable in ice_fileunits
    ! Set the nu_diag_set flag so it's not reset later

    call shr_log_setLogUnit (shrlogunit)
    call ufs_file_setLogUnit('./log.ice.timer',nu_timer,runtimelog)

    call NUOPC_CompAttributeGet(gcomp, name="diro", value=cvalue, &
         isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       diag_filename = trim(cvalue)
    end if
    call NUOPC_CompAttributeGet(gcomp, name="logfile", value=cvalue, &
         isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       diag_filename = trim(diag_filename) // '/' // trim(cvalue)
    end if

    if (trim(diag_filename) /= 'unset') then
       call set_component_logging(gcomp, mastertask, nu_diag, shrlogunit, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       nu_diag_set = .true.
    end if

    !----------------------------------------------------------------------------
    ! First cice initialization phase - before initializing grid info
    !----------------------------------------------------------------------------

#ifdef CESMCOUPLED
    ! Determine if single column

    call NUOPC_CompAttributeGet(gcomp, name='scol_lon', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scmlon
    call NUOPC_CompAttributeGet(gcomp, name='scol_lat', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scmlat
    call NUOPC_CompAttributeGet(gcomp, name='scol_spval', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scol_spval

    if (scmlon > scol_spval .and. scmlat > scol_spval) then
       call NUOPC_CompAttributeGet(gcomp, name='single_column_lnd_domainfile', &
            value=single_column_lnd_domainfile, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (trim(single_column_lnd_domainfile) /= 'UNSET') then
          single_column = .true.
       else
          call abort_ice('single_column_domainfile cannot be null for single column mode')
       end if
       call NUOPC_CompAttributeGet(gcomp, name='scol_ocnmask', value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) scol_mask
       call NUOPC_CompAttributeGet(gcomp, name='scol_ocnfrac', value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) scol_frac
       call NUOPC_CompAttributeGet(gcomp, name='scol_ni', value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) scol_ni
       call NUOPC_CompAttributeGet(gcomp, name='scol_nj', value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) scol_nj
       call NUOPC_CompAttributeGet(gcomp, name='scol_area', value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) scol_area

       call ice_mesh_create_scolumn(scmlon, scmlat, ice_mesh, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       scol_valid = (scol_mask == 1)
       if (.not. scol_valid) then
          ! Read the cice namelist as part of the call to cice_init1
          ! Note that if single_column is true and scol_valid is not - will never get here
          call t_startf ('cice_init1')
          call cice_init1
          call t_stopf ('cice_init1')
          ! Advertise fields
          call ice_advertise_fields(gcomp, importState, exportState, flds_scalar_name, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          call t_stopf ('cice_init_total')

          ! *******************
          ! *** RETURN HERE ***
          ! *******************
          RETURN
       end if
    end if
    ! Read the cice namelist as part of the call to cice_init1
    ! Note that if single_column is true and scol_valid is not - will never get here
    call t_startf ('cice_init1')
    call cice_init1
    call t_stopf ('cice_init1')

    !-----------------------------------------------------------------
    ! Advertise fields
    !-----------------------------------------------------------------
    call ice_advertise_fields(gcomp, importState, exportState, flds_scalar_name, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return


    ! Form of ocean freezing temperature
    ! 'minus1p8' = -1.8 C
    ! 'linear_salt' = -depressT * sss
    ! 'mushy' conforms with ktherm=2
    call NUOPC_CompAttributeGet(gcomp, name="tfreeze_option", value=tfrz_option_driver, &
         isPresent=isPresent, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (.not. isPresent) then
       tfrz_option_driver = 'linear_salt'
    end if
    call icepack_query_parameters( tfrz_option_out=tfrz_option)
    if (tfrz_option_driver  /= tfrz_option) then
       write(errmsg,'(a)') trim(subname)//'WARNING: tfrz_option from driver '//trim(tfrz_option_driver)//&
            ' is overwriting tfrz_option from cice namelist '//trim(tfrz_option)
       if (mastertask) write(nu_diag,*) trim(errmsg)
       call icepack_warnings_flush(nu_diag)
       call icepack_init_parameters(tfrz_option_in=tfrz_option_driver)
    endif

    ! Flux convergence tolerance - always use the driver attribute value
    call NUOPC_CompAttributeGet(gcomp, name="flux_convergence", value=cvalue, &
         isPresent=isPresent, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       read(cvalue,*) atmiter_conv_driver
       call icepack_query_parameters( atmiter_conv_out=atmiter_conv)
       if (atmiter_conv_driver /= atmiter_conv) then
          write(errmsg,'(a,d13.5,a,d13.5)') trim(subname)//'WARNING: atmiter_ from driver ',&
               atmiter_conv_driver,' is overwritting atmiter_conv from cice namelist ',atmiter_conv
          if(mastertask) write(nu_diag,*) trim(errmsg)
          call icepack_warnings_flush(nu_diag)
          call icepack_init_parameters(atmiter_conv_in=atmiter_conv_driver)
       end if
    end if

    ! Number of iterations for boundary layer calculations
    call NUOPC_CompAttributeGet(gcomp, name="flux_max_iteration", value=cvalue, isPresent=isPresent, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       read(cvalue,*) natmiter_driver
    else
       natmiter_driver = 5
    end if
    call icepack_query_parameters( natmiter_out=natmiter)
    if (natmiter_driver  /= natmiter) then
       write(errmsg,'(a,i8,a,i8)') trim(subname)//'error: natmiter_driver ',natmiter_driver, &
            ' must be the same as natmiter from cice namelist ',natmiter
       call abort_ice(trim(errmsg))
    endif

    ! Netcdf output created by PIO
    call NUOPC_CompAttributeGet(gcomp, name="pio_typename", value=cvalue, &
         isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
      if (trim(history_format)/='cdf1' .and. mastertask) then
         write(nu_diag,*) trim(subname)//history_format//'WARNING: history_format from cice_namelist ignored'
         write(nu_diag,*) trim(subname)//'WARNING: using '//trim(cvalue)//' from ICE_modelio'
      endif
      if (trim(restart_format)/='cdf1' .and. mastertask) then
         write(nu_diag,*) trim(subname)//restart_format//'WARNING: restart_format from cice_namelist ignored'
         write(nu_diag,*) trim(subname)//'WARNING: using '//trim(cvalue)//' from ICE_modelio'
      endif

      ! The only reason to set these is to detect in ice_history_write if the chunk/deflate settings are ok.
      select case (trim(cvalue))
      case ('netcdf4p')
         history_format='hdf5'
         restart_format='hdf5'
      case ('netcdf4c')
         if (mastertask) write(nu_diag,*) trim(subname)//'WARNING: pio_typename = netcdf4c is superseded, use netcdf4p'
         history_format='hdf5'
         restart_format='hdf5'
      case default !pio_typename=netcdf or pnetcdf
         ! do nothing
      end select
    else
      if(mastertask) write(nu_diag,*) trim(subname)//'WARNING: pio_typename from driver needs to be set for netcdf output to work'
    end if

#else

    ! Read the cice namelist as part of the call to cice_init1
    call t_startf ('cice_init1')
    call cice_init1
    call t_stopf ('cice_init1')

    !-----------------------------------------------------------------
    ! Advertise fields
    !-----------------------------------------------------------------
    call ice_advertise_fields(gcomp, importState, exportState, flds_scalar_name, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

#endif

    !----------------------------------------------------------------------------
    ! Initialize grid info
    !----------------------------------------------------------------------------

    if (single_column .and. scol_valid) then
       call ice_mesh_init_tlon_tlat_area_hm()
    else
       ! Determine mesh input file
       call NUOPC_CompAttributeGet(gcomp, name='mesh_ice', value=ice_meshfile, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Determine mask input file
       call NUOPC_CompAttributeGet(gcomp, name='mesh_mask', value=cvalue, isPresent=isPresent, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (isPresent .and. isSet) then
          ice_maskfile = trim(cvalue)
       else
          ice_maskfile = ice_meshfile
       end if
       if (my_task == master_task) then
          write(nu_diag,*)'mesh file for cice domain is ',trim(ice_meshfile)
          write(nu_diag,*)'mask file for cice domain is ',trim(ice_maskfile)
       end if

       ! Determine the model distgrid using the decomposition obtained in
       ! call to init_grid1 called from cice_init1
       call ice_mesh_set_distgrid(localpet, npes, ice_distgrid, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Read in the ice mesh on the cice distribution
       ice_mesh = ESMF_MeshCreate(filename=trim(ice_meshfile), fileformat=ESMF_FILEFORMAT_ESMFMESH, &
            elementDistGrid=ice_distgrid, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Initialize the cice mesh and the cice mask
       if (trim(grid_format) == 'meshnc') then
          ! In this case cap code determines the mask file
          call ice_mesh_setmask_from_maskfile(ice_maskfile, ice_mesh, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ice_mesh_init_tlon_tlat_area_hm()
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          ! In this case init_grid2 will initialize tlon, tlat, area and hm
          call init_grid2()
          call ice_mesh_check(gcomp,ice_mesh, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end if

    call t_stopf ('cice_init_total')
    if (mastertask) call ufs_logtimer(nu_timer,msec,'InitializeAdvertise time: ',runtimelog,wtime)
  end subroutine InitializeAdvertise

  !===============================================================================

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)

    ! Arguments
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! Local variables
    integer                                :: n
    integer                                :: fieldcount
    type(ESMF_Field)                       :: lfield
    character(len=char_len_long)           :: cvalue
    real(dbl_kind)                         :: scol_lon
    real(dbl_kind)                         :: scol_lat
    real(dbl_kind)                         :: scol_spval
    real(dbl_kind), pointer                :: fldptr1d(:)
    real(dbl_kind), pointer                :: fldptr2d(:,:)
    integer                                :: rank
    character(len=char_len)      :: tfrz_option    ! tfrz_option from cice namelist
    integer(int_kind)            :: ktherm

    character(len=char_len_long)           :: single_column_lnd_domainfile
    character(len=char_len_long) , pointer :: lfieldnamelist(:) => null()
    character(len=*), parameter            :: subname=trim(modName)//':(InitializeRealize) '
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    call ufs_settimer(wtime)
    !----------------------------------------------------------------------------
    ! Second cice initialization phase -after initializing grid info
    !----------------------------------------------------------------------------
    ! Note that cice_init2 also sets time manager info as well as mpi communicator info,
    ! including master_task and my_task
    ! Note that cice_init2 calls ice_init() which in turn calls icepack_init_parameters
    ! which sets the tfrz_option
    call t_startf ('cice_init2')
    call cice_init2()
    call t_stopf ('cice_init2')
    !---------------------------------------------------------------------------
    ! use EClock to reset calendar information
    !---------------------------------------------------------------------------

    ! - on initial run
    !   - iyear, month and mday obtained from sync clock
    !   - time determined from myear, month and mday
    !   - istep0 and istep1 are set to 0
    ! - on restart run
    !   - istep0, time and time_forc are read from restart file
    !   - istep1 is set to istep0
    !   - idate is determined from time via the call to calendar (see below)

    if (runtype == 'initial') then
       if (ref_ymd /= start_ymd .or. ref_tod /= start_tod) then
          if (my_task == master_task) then
             write(nu_diag,*) trim(subname),': ref_ymd ',ref_ymd, ' must equal start_ymd ',start_ymd
             write(nu_diag,*) trim(subname),': ref_tod',ref_tod, ' must equal start_tod ',start_tod
          end if
       end if

       if (my_task == master_task) then
          write(nu_diag,*) trim(subname),' idate from sync clock = ', start_ymd
          write(nu_diag,*) trim(subname),'   tod from sync clock = ', start_tod
          write(nu_diag,*) trim(subname),' resetting idate to match sync clock'
       end if
       idate = curr_ymd

       if (idate < 0) then
          if (my_task == master_task) then
             write(nu_diag,*) trim(subname),' ERROR curr_ymd,year_init =',curr_ymd,year_init
             write(nu_diag,*) trim(subname),' ERROR idate lt zero',idate
          end if
          call abort_ice(subname//' :: ERROR idate lt zero')
       endif
       myear = (idate/10000)                     ! integer year of basedate
       mmonth= (idate-myear*10000)/100           ! integer month of basedate
       mday  =  idate-myear*10000-mmonth*100     ! day of month of basedate
       msec  = start_tod                         ! start from basedate

       if (my_task == master_task) then
          write(nu_diag,*) trim(subname),' curr_ymd = ',curr_ymd
          write(nu_diag,*) trim(subname),' cice year_init = ',year_init
          write(nu_diag,*) trim(subname),' cice start date = ',idate
          write(nu_diag,*) trim(subname),' cice start ymds = ',myear,mmonth,mday,start_tod
          write(nu_diag,*) trim(subname),' cice calendar_type = ',trim(calendar_type)
       endif

    end if

    !  - start time from ESMF clock. Used to set history time units
    idate0    = start_ymd
    year_init = (idate0/10000)
    month_init= (idate0-year_init*10000)/100           ! integer month of basedate
    day_init  = idate0-year_init*10000-month_init*100

    !  - Set use_leap_years based on calendar (as some CICE calls use this instead of the calendar type)
    if (calendar_type == ice_calendar_proleptic_gregorian) then
      use_leap_years = .true.
    else
      use_leap_years = .false. ! no_leap calendars
    endif

    call calendar()     ! update calendar info

    !----------------------------------------------------------------------------
    ! reset shr logging to my log file
    !----------------------------------------------------------------------------

    call icepack_query_parameters(ktherm_out=ktherm)
    call icepack_query_parameters(tfrz_option_out=tfrz_option)
    call icepack_warnings_flush(nu_diag)
    if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
        file=__FILE__, line=__LINE__)

    ! Now write output to nu_diag - this must happen AFTER call to cice_init
    if (mastertask) then
       write(nu_diag,'(a,d21.14)') trim(subname)//' cice init nextsw_cday = ',nextsw_cday
       write(nu_diag,'(a)') trim(subname)//' tfrz_option = '//trim(tfrz_option)
       if (ktherm == 2 .and. trim(tfrz_option) /= 'mushy') then
          write(nu_diag,*) trim(subname),' Warning: Using ktherm = 2 and tfrz_option = ', trim(tfrz_option)
       endif
       write(nu_diag,'(a    )') trim(subname)//' inst_name   = '//trim(inst_name)
       write(nu_diag,'(a,i8 )') trim(subname)//' inst_index  = ',inst_index
       write(nu_diag,'(a    )') trim(subname)//' inst_suffix = ',trim(inst_suffix)
    endif


    if (write_ic) then
       call accum_hist(dt)  ! write initial conditions
    end if

    !-----------------------------------------------------------------
    ! Prescribed ice initialization
    !-----------------------------------------------------------------

    call ice_prescribed_init(gcomp, clock, ice_mesh, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

#ifdef CESMCOUPLED
    ! if single column is not valid - set all export state fields to zero and return
    if (single_column .and. .not. scol_valid) then
       write(nu_diag,'(a)')' (ice_comp_nuopc) single column mode point does not contain any ocn/ice '&
            //' - setting all export data to 0'
       call ice_realize_fields(gcomp, mesh=ice_mesh, &
            flds_scalar_name=flds_scalar_name, flds_scalar_num=flds_scalar_num, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_StateGet(exportState, itemCount=fieldCount, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       allocate(lfieldnamelist(fieldCount))
       call ESMF_StateGet(exportState, itemNameList=lfieldnamelist, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       do n = 1, fieldCount
          if (trim(lfieldnamelist(n)) /= flds_scalar_name) then
             call ESMF_StateGet(exportState, itemName=trim(lfieldnamelist(n)), field=lfield, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call ESMF_FieldGet(lfield, rank=rank, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             if (rank == 2) then
                call ESMF_FieldGet(lfield, farrayPtr=fldptr2d, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                fldptr2d(:,:) = 0._dbl_kind
             else
                call ESMF_FieldGet(lfield, farrayPtr=fldptr1d, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                fldptr1d(:) = 0._dbl_kind
             end if
          end if
       enddo
       deallocate(lfieldnamelist)
       call State_SetScalar(dble(0), flds_scalar_index_nx, exportState, &
            flds_scalar_name, flds_scalar_num, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call State_SetScalar(dble(0), flds_scalar_index_ny, exportState, &
            flds_scalar_name, flds_scalar_num, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! *******************
       ! *** RETURN HERE ***
       ! *******************
       RETURN
    else if(single_column) then
       write(nu_diag,'(a,3(f10.5,2x))')' (ice_comp_nuopc) single column mode lon/lat/frac is ',&
            scmlon,scmlat,scol_frac
    end if
#endif

    !-----------------------------------------------------------------
    ! Realize the actively coupled fields
    !-----------------------------------------------------------------

    call ice_realize_fields(gcomp, mesh=ice_mesh, &
         flds_scalar_name=flds_scalar_name, flds_scalar_num=flds_scalar_num, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !-----------------------------------------------------------------
    ! Create cice export state
    !-----------------------------------------------------------------

    call ice_export (exportstate, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call State_SetScalar(dble(nx_global), flds_scalar_index_nx, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_SetScalar(dble(ny_global), flds_scalar_index_ny, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! diagnostics
    !--------------------------------

    ! TODO (mvertens, 2018-12-21): fill in iceberg_prognostic as .false.
    if (debug_export > 0 .and. my_task==master_task) then
       call State_fldDebug(exportState, flds_scalar_name, 'cice_export:', &
            idate, msec, nu_diag, rc=rc)
    end if

    if (dbug > 0) then
       call state_diagnose(exportState,subname//':ES',rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

    call flush_fileunit(nu_diag)

    if (mastertask) call ufs_logtimer(nu_timer,msec,'InitializeRealize time: ',runtimelog,wtime)
  end subroutine InitializeRealize

  !===============================================================================

  subroutine ModelAdvance(gcomp, rc)

    !---------------------------------------------------------------------------
    ! Run CICE
    !---------------------------------------------------------------------------

    ! Arguments
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Local variables
    type(ESMF_Clock)           :: clock
    type(ESMF_Alarm)           :: alarm
    type(ESMF_Time)            :: startTime
    type(ESMF_Time)            :: currTime
    type(ESMF_Time)            :: nextTime
    type(ESMF_TimeInterval)    :: timeStep
    type(ESMF_State)           :: importState, exportState
    character(ESMF_MAXSTR)     :: cvalue
    real(dbl_kind)             :: eccen, obliqr, lambm0, mvelpp
    integer                    :: shrlogunit ! original log unit
    integer                    :: k,n        ! index
    logical                    :: stop_now   ! .true. ==> stop at the end of this run phase
    integer                    :: ymd        ! Current date (YYYYMMDD)
    integer                    :: tod        ! Current time of day (sec)
    integer                    :: curr_ymd   ! Current date (YYYYMMDD)
    integer                    :: curr_tod   ! Current time of day (s)
    integer                    :: yy,mm,dd   ! year, month, day, time of day
    integer                    :: ymd_sync   ! Sync date (YYYYMMDD)
    integer                    :: yr_sync    ! Sync current year
    integer                    :: mon_sync   ! Sync current month
    integer                    :: day_sync   ! Sync current day
    integer                    :: tod_sync   ! Sync current time of day (sec)
    character(char_len_long)   :: restart_date
    character(char_len_long)   :: restart_filename
    logical                    :: isPresent, isSet
#ifndef CESMCOUPLED
    logical                    :: write_restartfh
#endif
    character(len=*),parameter :: subname=trim(modName)//':(ModelAdvance) '
    character(char_len_long)   :: msgString
    !--------------------------------

    rc = ESMF_SUCCESS
    if (mastertask) call ufs_logtimer(nu_timer,msec,'ModelAdvance time since last step: ',runtimelog,wtime)
    call ufs_settimer(wtime)

    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    if (single_column .and. .not. scol_valid) then
       ! *******************
       ! *** RETURN HERE ***
       ! *******************
       RETURN
    end if

    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockPrint(clock, options="currTime", &
      preString="------>Advancing ICE from: ", unit=msgString, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(subname//trim(msgString), ESMF_LOGMSG_INFO)

    call ESMF_ClockGet(clock, startTime=startTime, currTime=currTime, &
      timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimePrint(currTime + timeStep, &
      preString="--------------------------------> to: ", unit=msgString, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)

    !--------------------------------
    ! Turn on timers
    !--------------------------------

    call ice_timer_start(timer_total) ! time entire run
    call t_barrierf('cice_run_total_BARRIER',mpi_comm_ice)
    call t_startf ('cice_run_total')

    !--------------------------------
    ! Reset shr logging to my log file
    !--------------------------------

    call shr_log_getLogUnit (shrlogunit)
    call shr_log_setLogUnit (nu_diag)

    !--------------------------------
    ! Query the Component for its clock, importState and exportState
    !--------------------------------

    call NUOPC_ModelGet(gcomp, modelClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Determine time of next atmospheric shortwave calculation
    !--------------------------------

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxNextSwCday", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       call State_GetScalar(importState, flds_scalar_index_nextsw_cday, nextsw_cday, &
            flds_scalar_name, flds_scalar_num, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call ESMF_ClockGetNextTime(clock, nextTime, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeGet(nextTime, dayOfYear_r8=nextsw_cday, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    if (my_task == master_task) then
       write(nu_diag,'(a,2x,i8,2x,d24.14)') trim(subname)//' cice istep, nextsw_cday = ',istep, nextsw_cday
    end if

    !--------------------------------
    ! Obtain orbital values
    !--------------------------------
    call ice_orbital_init(gcomp, clock, nu_diag, my_task==master_task, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! check that cice internal time is in sync with master clock before timestep update
    !--------------------------------

    ! cice clock
    tod = msec
    ymd = idate

    ! model clock
    call ESMF_ClockGet( clock, currTime=currTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet( currTime, yy=yr_sync, mm=mon_sync, dd=day_sync, s=tod_sync, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ice_cal_ymd2date(yr_sync, mon_sync, day_sync, ymd_sync)

    ! error check
    if ( (ymd /= ymd_sync) .or. (tod /= tod_sync) ) then
       if (my_task == master_task) then
          write(nu_diag,*)' cice ymd=',ymd     ,'  cice tod= ',tod
          write(nu_diag,*)' sync ymd=',ymd_sync,'  sync tod= ',tod_sync
       end if
       call ESMF_LogWrite(subname//" CICE clock not in sync with ESMF model clock",ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return
    end if

    !--------------------------------
    ! Determine if time to write restart
    !--------------------------------

    ! Note this logic triggers off of the component clock rather than the internal cice time
    ! The component clock does not get advanced until the end of the loop - not at the beginning

    call ESMF_ClockGetAlarm(clock, alarmname='alarm_restart', alarm=alarm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    force_restart_now = .false.

    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       force_restart_now = .true.
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_ClockGetNextTime(clock, nextTime=nextTime, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeGet(nexttime, yy=yy, mm=mm, dd=dd, s=tod, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       write(restart_date,"(i4.4,a,i2.2,a,i2.2,a,i5.5)") yy, '-', mm, '-',dd,'-',tod
       write(restart_filename,'(4a)') trim(restart_dir), trim(restart_file), '.', trim(restart_date)
    endif

    ! Handle end of run restart
    if (restart_eor) then
       call ESMF_ClockGetAlarm(clock, alarmname='alarm_stop', alarm=alarm, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          force_restart_now = .true.
          call ESMF_AlarmRingerOff( alarm, rc=rc )
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       endif
    endif

#ifndef CESMCOUPLED
    call is_restart_fh(clock, restartfh_info, write_restartfh)
    if (write_restartfh) force_restart_now = .true.
#endif

    !--------------------------------
    ! Unpack import state
    !--------------------------------

    call t_barrierf('cice_run_import_BARRIER',mpi_comm_ice)
    call t_startf ('cice_run_import')
    call ice_import(importState, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf ('cice_run_import')

    ! write Debug output
    if (debug_import  > 0 .and. my_task==master_task) then
       call State_fldDebug(importState, flds_scalar_name, 'cice_import:', &
            idate, msec, nu_diag, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    if (dbug > 0) then
       call state_diagnose(importState,subname//':IS',rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !--------------------------------
    ! Advance cice and timestep update
    !--------------------------------

    if(profile_memory) call ESMF_VMLogMemInfo("Entering CICE_Run : ")
    call CICE_Run()
    if(profile_memory) call ESMF_VMLogMemInfo("Leaving CICE_Run : ")

    !--------------------------------
    ! Create export state
    !--------------------------------

    call t_barrierf('cice_run_export_BARRIER',mpi_comm_ice)
    call t_startf ('cice_run_export')
    call ice_export(exportState, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf ('cice_run_export')

    ! write Debug output
    if (debug_export > 0 .and. my_task==master_task) then
       call State_fldDebug(exportState, flds_scalar_name, 'cice_export:', &
            idate, msec, nu_diag, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    if (dbug > 0) then
       call state_diagnose(exportState,subname//':ES',rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! reset shr logging to my original values
    call shr_log_setLogUnit (shrlogunit)

    !--------------------------------
    ! stop timers and print timer info
    !--------------------------------
    ! Need to have this logic here instead of in finalize phase
    ! since the finalize phase will still be called even in aqua-planet mode

    !--------------------------------
    ! Determine if time to stop
    !--------------------------------

    call ESMF_ClockGetAlarm(clock, alarmname='alarm_stop', alarm=alarm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       stop_now = .true.
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       stop_now = .false.
    endif

    call t_stopf ('cice_run_total')

    ! Need to stop this at the end of every run phase in a coupled run.
    call ice_timer_stop(timer_total)
    if (stop_now) then
       call ice_timer_print_all(stats=.true.) ! print timing information
       call release_all_fileunits
    endif

  105  format( A, 2i8, A, f10.2, A, f10.2, A)

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

    if (mastertask) call ufs_logtimer(nu_timer,msec,'ModelAdvance time: ',runtimelog,wtime)
    call ufs_settimer(wtime)

  end subroutine ModelAdvance

  !===============================================================================

  subroutine ModelSetRunClock(gcomp, rc)

    ! intput/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)         :: mclock, dclock
    type(ESMF_Time)          :: mcurrtime, dcurrtime
    type(ESMF_Time)          :: mstoptime
    type(ESMF_TimeInterval)  :: mtimestep, dtimestep
    character(len=256)       :: cvalue
    character(len=256)       :: restart_option ! Restart option units
    integer                  :: restart_n      ! Number until restart interval
    integer                  :: restart_ymd    ! Restart date (YYYYMMDD)
    type(ESMF_ALARM)         :: restart_alarm
    character(len=256)       :: stop_option    ! Stop option units
    integer                  :: stop_n         ! Number until stop interval
    integer                  :: stop_ymd       ! Stop date (YYYYMMDD)
    type(ESMF_ALARM)         :: stop_alarm
    character(len=128)       :: name
    integer                  :: alarmcount
#ifndef CESMCOUPLED
    integer                  :: dtime
#endif
    character(len=*),parameter :: subname=trim(modName)//':(ModelSetRunClock) '
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! query the Component for its clocks
    call NUOPC_ModelGet(gcomp, driverClock=dclock, modelClock=mclock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(dclock, currTime=dcurrtime, timeStep=dtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(mclock, currTime=mcurrtime, timeStep=mtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! force model clock currtime and timestep to match driver and set stoptime
    !--------------------------------

    mstoptime = mcurrtime + dtimestep
    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! set restart and stop alarms
    !--------------------------------

    call ESMF_ClockGetAlarmList(mclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmCount=alarmCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (alarmCount == 0) then

       call ESMF_GridCompGet(gcomp, name=name, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite(subname//'setting alarms for ' // trim(name), ESMF_LOGMSG_INFO)

       !----------------
       ! Stop alarm
       !----------------
       call NUOPC_CompAttributeGet(gcomp, name="stop_option", value=stop_option, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="stop_n", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) stop_n

       call NUOPC_CompAttributeGet(gcomp, name="stop_ymd", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) stop_ymd

       call alarmInit(mclock, stop_alarm, stop_option, &
            opt_n   = stop_n,           &
            opt_ymd = stop_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_stop', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(stop_alarm, clock=mclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !----------------
       ! Restart alarm
       !----------------
       call NUOPC_CompAttributeGet(gcomp, name="restart_option", value=restart_option, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="restart_n", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_n

       call NUOPC_CompAttributeGet(gcomp, name="restart_ymd", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_ymd

       call alarmInit(mclock, restart_alarm, restart_option, &
            opt_n   = restart_n,           &
            opt_ymd = restart_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_restart', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(restart_alarm, clock=mclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

#ifndef CESMCOUPLED
       call ESMF_TimeIntervalGet( dtimestep, s=dtime, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call init_is_restart_fh(mcurrTime, dtime, my_task == master_task, restartfh_info)
#endif
    end if

    !--------------------------------
    ! Advance model clock to trigger alarms then reset model clock back to currtime
    !--------------------------------

    call ESMF_ClockAdvance(mclock,rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelSetRunClock

  !===============================================================================

  subroutine ModelFinalize(gcomp, rc)

    !--------------------------------
    ! Finalize routine
    !--------------------------------

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(*), parameter :: F91 = "('(ice_comp_nuopc) ',73('-'))"
    character(len=*),parameter  :: subname=trim(modName)//':(ModelFinalize) '
    !--------------------------------

    rc = ESMF_SUCCESS
    call ufs_settimer(wtime)
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
    if (my_task == master_task) then
       write(nu_diag,F91)
       write(nu_diag,'(a)') 'CICE: end of main integration loop'
       write(nu_diag,F91)
    end if
    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

    if(mastertask) call ufs_logtimer(nu_timer,msec,'ModelFinalize time: ',runtimelog,wtime)

  end subroutine ModelFinalize

  !===============================================================================

  subroutine ice_orbital_init(gcomp, clock, logunit, mastertask, rc)

    !----------------------------------------------------------
    ! Initialize orbital related values for cesm coupled
    !----------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp)                 :: gcomp
    type(ESMF_Clock)    , intent(in)    :: clock
    integer             , intent(in)    :: logunit
    logical             , intent(in)    :: mastertask
    integer             , intent(out)   :: rc              ! output error

    ! local variables
    real(dbl_kind)               :: eccen, obliqr, lambm0, mvelpp
    character(len=char_len_long) :: msgstr   ! temporary
    character(len=char_len_long) :: cvalue   ! temporary
    type(ESMF_Time)              :: CurrTime ! current time
    integer                      :: year     ! model year at current time
    integer                      :: orb_year ! orbital year for current orbital computation
    integer, save                :: prev_orb_year=0 ! orbital year for previous orbital computation
    logical                      :: lprint
    logical, save                :: first_time = .true.
    character(len=*) , parameter :: subname = "(cice_orbital_init)"
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

#ifndef CESMCOUPLED
    return
#else
    if (first_time) then

       ! Determine orbital attributes from input
       call NUOPC_CompAttributeGet(gcomp, name='orb_mode', value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) orb_mode
       call NUOPC_CompAttributeGet(gcomp, name="orb_iyear", value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) orb_iyear
       call NUOPC_CompAttributeGet(gcomp, name="orb_iyear_align", value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) orb_iyear_align
       call NUOPC_CompAttributeGet(gcomp, name="orb_obliq", value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) orb_obliq
       call NUOPC_CompAttributeGet(gcomp, name="orb_eccen", value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) orb_eccen
       call NUOPC_CompAttributeGet(gcomp, name="orb_mvelp", value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) orb_mvelp

       ! Error checks
       if (trim(orb_mode) == trim(orb_fixed_year)) then
          orb_obliq = SHR_ORB_UNDEF_REAL
          orb_eccen = SHR_ORB_UNDEF_REAL
          orb_mvelp = SHR_ORB_UNDEF_REAL
          if (orb_iyear == SHR_ORB_UNDEF_INT) then
             if (mastertask) then
                write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(orb_mode)
                write(logunit,*) trim(subname),' ERROR: fixed_year settings = ',orb_iyear
                write (msgstr, *) ' ERROR: invalid settings for orb_mode '//trim(orb_mode)
             end if
             call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
             return  ! bail out
          endif
       elseif (trim(orb_mode) == trim(orb_variable_year)) then
          orb_obliq = SHR_ORB_UNDEF_REAL
          orb_eccen = SHR_ORB_UNDEF_REAL
          orb_mvelp = SHR_ORB_UNDEF_REAL
          if (orb_iyear == SHR_ORB_UNDEF_INT .or. orb_iyear_align == SHR_ORB_UNDEF_INT) then
             if (mastertask) then
                write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(orb_mode)
                write(logunit,*) trim(subname),' ERROR: variable_year settings = ',orb_iyear, orb_iyear_align
                write (msgstr, *) subname//' ERROR: invalid settings for orb_mode '//trim(orb_mode)
             end if
             call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
             return  ! bail out
          endif
       elseif (trim(orb_mode) == trim(orb_fixed_parameters)) then
          !-- force orb_iyear to undef to make sure shr_orb_params works properly
          orb_iyear = SHR_ORB_UNDEF_INT
          orb_iyear_align = SHR_ORB_UNDEF_INT
          if (orb_eccen == SHR_ORB_UNDEF_REAL .or. &
               orb_obliq == SHR_ORB_UNDEF_REAL .or. &
               orb_mvelp == SHR_ORB_UNDEF_REAL) then
             if (mastertask) then
                write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(orb_mode)
                write(logunit,*) trim(subname),' ERROR: orb_eccen = ',orb_eccen
                write(logunit,*) trim(subname),' ERROR: orb_obliq = ',orb_obliq
                write(logunit,*) trim(subname),' ERROR: orb_mvelp = ',orb_mvelp
                write (msgstr, *) subname//' ERROR: invalid settings for orb_mode '//trim(orb_mode)
             end if
             call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
             return  ! bail out
          endif
       else
          write (msgstr, *) subname//' ERROR: invalid orb_mode '//trim(orb_mode)
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
          rc = ESMF_FAILURE
          return  ! bail out
       endif
    end if
    lprint = .false.
    if (trim(orb_mode) == trim(orb_variable_year)) then
       call ESMF_ClockGet(clock, CurrTime=CurrTime, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeGet(CurrTime, yy=year, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       orb_year = orb_iyear + (year - orb_iyear_align)
    else
       orb_year = orb_iyear
    end if

    if (orb_year .ne. prev_orb_year) then
       lprint = mastertask
       ! this prevents the orbital print happening before the log file is opened.
       if (.not. first_time) prev_orb_year = orb_year
    endif
    eccen = orb_eccen

    call shr_orb_params(orb_year, eccen, orb_obliq, orb_mvelp, obliqr, lambm0, mvelpp, lprint)

    if ( eccen  == SHR_ORB_UNDEF_REAL .or. obliqr == SHR_ORB_UNDEF_REAL .or. &
         mvelpp == SHR_ORB_UNDEF_REAL .or. lambm0 == SHR_ORB_UNDEF_REAL) then
       write (msgstr, *) subname//' ERROR: orb params incorrect'
       call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
       return  ! bail out
    endif

    call icepack_init_orbit(eccen_in=eccen, mvelpp_in=mvelpp, lambm0_in=lambm0, obliqr_in=obliqr)
    call icepack_warnings_flush(nu_diag)
    if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
        file=__FILE__, line=__LINE__)

    first_time = .false.
#endif
  end subroutine ice_orbital_init

  !===============================================================================
  subroutine ice_cal_ymd2date(year, month, day, date)

    ! input/output parameters:
    integer,intent(in ) :: year,month,day  ! calendar year,month,day
    integer,intent(out) :: date            ! coded (yyyymmdd) calendar date

    !--- local ---
    character(*),parameter :: subName = "(ice_cal_ymd2date)"
    !-------------------------------------------------------------------------------
    ! NOTE:
    !   this calendar has a year zero (but no day or month zero)
    !-------------------------------------------------------------------------------

    date = abs(year)*10000 + month*100 + day  ! coded calendar date
    if (year < 0) date = -date

  end subroutine ice_cal_ymd2date

end module ice_comp_nuopc
