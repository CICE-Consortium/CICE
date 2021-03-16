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
  use ice_constants      , only : ice_init_constants
  use ice_shr_methods    , only : chkerr, state_setscalar, state_getscalar, state_diagnose, alarmInit
  use ice_shr_methods    , only : set_component_logging, get_component_instance
  use ice_shr_methods    , only : state_flddebug
  use ice_import_export  , only : ice_import, ice_export
  use ice_import_export  , only : ice_advertise_fields, ice_realize_fields
  use ice_domain_size    , only : nx_global, ny_global
  use ice_domain         , only : nblocks, blocks_ice, distrb_info
  use ice_blocks         , only : block, get_block, nx_block, ny_block, nblocks_x, nblocks_y
  use ice_blocks         , only : nblocks_tot, get_block_parameter
  use ice_distribution   , only : ice_distributiongetblockloc
  use ice_grid           , only : tlon, tlat, hm, tarea, ULON, ULAT
  use ice_communicate    , only : init_communicate, my_task, master_task, mpi_comm_ice
  use ice_calendar       , only : force_restart_now, write_ic
  use ice_calendar       , only : idate, mday, time, month, daycal, time2sec, year_init
  use ice_calendar       , only : sec, dt, calendar, calendar_type, nextsw_cday, istep
  use ice_kinds_mod      , only : dbl_kind, int_kind, char_len, char_len_long
  use ice_scam           , only : scmlat, scmlon, single_column
  use ice_fileunits      , only : nu_diag, nu_diag_set, inst_index, inst_name
  use ice_fileunits      , only : inst_suffix, release_all_fileunits, flush_fileunit
  use ice_restart_shared , only : runid, runtype, restart_dir, restart_file
  use ice_history        , only : accum_hist
  use CICE_InitMod       , only : cice_init
  use CICE_RunMod        , only : cice_run
  use ice_exit           , only : abort_ice
  use icepack_intfc      , only : icepack_warnings_flush, icepack_warnings_aborted
  use icepack_intfc      , only : icepack_init_orbit, icepack_init_parameters, icepack_query_orbit
  use icepack_intfc      , only : icepack_query_tracer_flags, icepack_query_parameters
  use cice_wrapper_mod   , only : t_startf, t_stopf, t_barrierf
  use cice_wrapper_mod   , only : shr_file_getlogunit, shr_file_setlogunit
#ifdef CESMCOUPLED
  use shr_const_mod
  use shr_orb_mod        , only : shr_orb_decl, shr_orb_params, SHR_ORB_UNDEF_REAL, SHR_ORB_UNDEF_INT
#endif
  use ice_timers
  use ice_prescribed_mod , only : ice_prescribed_init

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

  character(len=*),parameter   :: shr_cal_noleap    = 'NO_LEAP'
  character(len=*),parameter   :: shr_cal_gregorian = 'GREGORIAN'

  integer                      :: dbug = 0
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
    !--------------------------------

    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, &
         acceptStringList=(/"IPDv01p"/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

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
    character(len=char_len_long) :: logmsg
    logical                      :: isPresent, isSet
    character(len=*), parameter :: subname=trim(modName)//':(InitializeAdvertise) '
    !--------------------------------

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

    call ice_advertise_fields(gcomp, importState, exportState, flds_scalar_name, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

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
    real(dbl_kind)               :: eccen, obliqr, lambm0, mvelpp
    type(ESMF_DistGrid)          :: distGrid
    type(ESMF_Mesh)              :: Emesh, EmeshTemp
    integer                      :: spatialDim
    integer                      :: numOwnedElements
    real(dbl_kind), pointer      :: ownedElemCoords(:)
    real(dbl_kind), pointer      :: lat(:), latMesh(:)
    real(dbl_kind), pointer      :: lon(:), lonMesh(:)
    integer , allocatable        :: gindex_ice(:)
    integer , allocatable        :: gindex_elim(:)
    integer , allocatable        :: gindex(:)
    integer                      :: globalID
    character(ESMF_MAXSTR)       :: cvalue
    character(len=char_len)      :: tfrz_option
    character(ESMF_MAXSTR)       :: convCIM, purpComp
    type(ESMF_VM)                :: vm
    type(ESMF_Time)              :: currTime           ! Current time
    type(ESMF_Time)              :: startTime          ! Start time
    type(ESMF_Time)              :: stopTime           ! Stop time
    type(ESMF_Time)              :: refTime            ! Ref time
    type(ESMF_TimeInterval)      :: timeStep           ! Model timestep
    type(ESMF_Calendar)          :: esmf_calendar      ! esmf calendar
    type(ESMF_CalKind_Flag)      :: esmf_caltype       ! esmf calendar type
    integer                      :: start_ymd          ! Start date (YYYYMMDD)
    integer                      :: start_tod          ! start time of day (s)
    integer                      :: curr_ymd           ! Current date (YYYYMMDD)
    integer                      :: curr_tod           ! Current time of day (s)
    integer                      :: stop_ymd           ! stop date (YYYYMMDD)
    integer                      :: stop_tod           ! stop time of day (sec)
    integer                      :: ref_ymd            ! Reference date (YYYYMMDD)
    integer                      :: ref_tod            ! reference time of day (s)
    integer                      :: yy,mm,dd           ! Temporaries for time query
    integer                      :: iyear              ! yyyy
    integer                      :: dtime              ! time step
    integer                      :: lmpicom
    integer                      :: shrlogunit         ! original log unit
    character(len=char_len)      :: starttype          ! infodata start type
    integer                      :: lsize              ! local size of coupling array
    logical                      :: isPresent
    logical                      :: isSet
    integer                      :: localPet
    integer                      :: n,c,g,i,j,m        ! indices
    integer                      :: iblk, jblk         ! indices
    integer                      :: ig, jg             ! indices
    integer                      :: ilo, ihi, jlo, jhi ! beginning and end of physical domain
    type(block)                  :: this_block         ! block information for current block
    integer                      :: compid             ! component id
    character(len=char_len_long) :: tempc1,tempc2
    real(dbl_kind)               :: diff_lon
    integer                      :: npes
    integer                      :: num_elim_global
    integer                      :: num_elim_local
    integer                      :: num_elim
    integer                      :: num_ice
    integer                      :: num_elim_gcells    ! local number of eliminated gridcells
    integer                      :: num_elim_blocks    ! local number of eliminated blocks
    integer                      :: num_total_blocks
    integer                      :: my_elim_start, my_elim_end
    real(dbl_kind)               :: rad_to_deg
    integer(int_kind)            :: ktherm
    logical                      :: mastertask
    character(len=char_len_long) :: diag_filename = 'unset'
    character(len=*), parameter  :: F00   = "('(ice_comp_nuopc) ',2a,1x,d21.14)"
    character(len=*), parameter  :: subname=trim(modName)//':(InitializeRealize) '
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !----------------------------------------------------------------------------
    ! generate local mpi comm
    !----------------------------------------------------------------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, mpiCommunicator=lmpicom, localPet=localPet, PetCount=npes, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

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

    inst_name = "ICE"//trim(inst_suffix)

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
       snowpatch_in        = 0.005_dbl_kind,                  &
       dragio_in           = 0.00962_dbl_kind)

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
       else if (trim(starttype) == trim('branch')) then
          runtype = "continue"
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

    ! Determine if single column
    call NUOPC_CompAttributeGet(gcomp, name='single_column', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) single_column
       if (single_column) then
          call NUOPC_CompAttributeGet(gcomp, name='scmlon', value=cvalue, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          read(cvalue,*) scmlon
          call NUOPC_CompAttributeGet(gcomp, name='scmlat', value=cvalue, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          read(cvalue,*) scmlat
       end if
    else
       single_column = .false.
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
       calendar_type = shr_cal_noleap
    else if (esmf_caltype == ESMF_CALKIND_GREGORIAN) then
       calendar_type = shr_cal_gregorian
    else
       call abort_ice( subname//'ERROR:: bad calendar for ESMF' )
    end if

    !----------------------------------------------------------------------------
    ! Set cice logging
    !----------------------------------------------------------------------------
    ! Note - this must be done AFTER the communicators are set
    ! Note that sets the nu_diag module variable in ice_fileunits
    ! Set the nu_diag_set flag so it's not reset later

    call shr_file_setLogUnit (shrlogunit)

    call NUOPC_CompAttributeGet(gcomp, name="diro", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       diag_filename = trim(cvalue)
    end if
    call NUOPC_CompAttributeGet(gcomp, name="logfile", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
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
    ! Initialize cice
    !----------------------------------------------------------------------------

    ! Note that cice_init also sets time manager info as well as mpi communicator info,
    ! including master_task and my_task

    call t_startf ('cice_init')
    call cice_init
    call t_stopf ('cice_init')

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
       write(nu_diag,F00) trim(subname),' cice init nextsw_cday = ',nextsw_cday
       write(nu_diag,*) trim(subname),' tfrz_option = ',trim(tfrz_option)
       if (ktherm == 2 .and. trim(tfrz_option) /= 'mushy') then
          write(nu_diag,*) trim(subname),' Warning: Using ktherm = 2 and tfrz_option = ', trim(tfrz_option)
       endif
       write(nu_diag,*) trim(subname),' inst_name   = ',trim(inst_name)
       write(nu_diag,*) trim(subname),' inst_index  = ',inst_index
       write(nu_diag,*) trim(subname),' inst_suffix = ',trim(inst_suffix)
    endif

    !---------------------------------------------------------------------------
    ! use EClock to reset calendar information on initial start
    !---------------------------------------------------------------------------

    ! - on initial run
    !   - iyear, month and mday obtained from sync clock
    !   - time determined from iyear, month and mday
    !   - istep0 and istep1 are set to 0
    ! - on restart run
    !   - istep0, time and time_forc are read from restart file
    !   - istep1 is set to istep0
    !   - idate is determined from time via the call to calendar (see below)

    if (runtype == 'initial') then
       if (ref_ymd /= start_ymd .or. ref_tod /= start_tod) then
          if (my_task == master_task) then
             write(nu_diag,*) trim(subname),': ref_ymd ',ref_ymd, ' must equal start_ymd ',start_ymd
             write(nu_diag,*) trim(subname),': ref_ymd ',ref_tod, ' must equal start_ymd ',start_tod
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
       iyear = (idate/10000)                     ! integer year of basedate
       month = (idate-iyear*10000)/100           ! integer month of basedate
       mday  =  idate-iyear*10000-month*100      ! day of month of basedate

       if (my_task == master_task) then
          write(nu_diag,*) trim(subname),' curr_ymd = ',curr_ymd
          write(nu_diag,*) trim(subname),' cice year_init = ',year_init
          write(nu_diag,*) trim(subname),' cice start date = ',idate
          write(nu_diag,*) trim(subname),' cice start ymds = ',iyear,month,mday,start_tod
          write(nu_diag,*) trim(subname),' cice calendar_type = ',trim(calendar_type)
       endif

#ifdef CESMCOUPLED
       if (calendar_type == "GREGORIAN" .or. &
           calendar_type == "Gregorian" .or. &
           calendar_type == "gregorian") then
          call time2sec(iyear-(year_init-1),month,mday,time)
       else
          call time2sec(iyear-year_init,month,mday,time)
       endif
#endif
       time = time+start_tod
    end if

    call calendar(time)     ! update calendar info
    if (write_ic) then
       call accum_hist(dt)  ! write initial conditions
    end if

    !---------------------------------------------------------------------------
    ! Determine the global index space needed for the distgrid
    !---------------------------------------------------------------------------

    ! number the local grid to get allocation size for gindex_ice
    lsize = 0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       do j = jlo, jhi
          do i = ilo, ihi
             lsize = lsize + 1
          enddo
       enddo
    enddo

    ! set global index array
    allocate(gindex_ice(lsize))
    n = 0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       do j = jlo, jhi
          do i = ilo, ihi
             n = n+1
             ig = this_block%i_glob(i)
             jg = this_block%j_glob(j)
             gindex_ice(n) = (jg-1)*nx_global + ig
          enddo
       enddo
    enddo

    ! Determine total number of eliminated blocks globally
    globalID = 0
    num_elim_global = 0  ! number of eliminated blocks
    num_total_blocks = 0
    do jblk=1,nblocks_y
       do iblk=1,nblocks_x
          globalID = globalID + 1
          num_total_blocks = num_total_blocks + 1
          if (distrb_info%blockLocation(globalID) == 0) then
             num_elim_global = num_elim_global + 1
          end if
       end do
    end do

    if (num_elim_global > 0) then

       ! Distribute the eliminated blocks in a round robin fashion amoung processors
       num_elim_local = num_elim_global / npes
       my_elim_start = num_elim_local*localPet + min(localPet, mod(num_elim_global, npes)) + 1
       if (localPet < mod(num_elim_global, npes)) then
          num_elim_local = num_elim_local + 1
       end if
       my_elim_end = my_elim_start + num_elim_local - 1

       ! Determine the number of eliminated gridcells locally
       globalID = 0
       num_elim_blocks = 0  ! local number of eliminated blocks
       num_elim_gcells = 0
       do jblk=1,nblocks_y
          do iblk=1,nblocks_x
             globalID = globalID + 1
             if (distrb_info%blockLocation(globalID) == 0) then
                num_elim_blocks = num_elim_blocks + 1
                if (num_elim_blocks >= my_elim_start .and. num_elim_blocks <= my_elim_end) then
                   this_block = get_block(globalID, globalID)
                   num_elim_gcells = num_elim_gcells + &
                        (this_block%jhi-this_block%jlo+1) * (this_block%ihi-this_block%ilo+1)
                end if
             end if
          end do
       end do

       ! Determine the global index space of the eliminated gridcells
       allocate(gindex_elim(num_elim_gcells))
       globalID = 0
       num_elim_gcells = 0  ! local number of eliminated gridcells
       num_elim_blocks = 0  ! local number of eliminated blocks
       do jblk=1,nblocks_y
          do iblk=1,nblocks_x
             globalID = globalID + 1
             if (distrb_info%blockLocation(globalID) == 0) then
                this_block = get_block(globalID, globalID)
                num_elim_blocks = num_elim_blocks + 1
                if (num_elim_blocks >= my_elim_start .and. num_elim_blocks <= my_elim_end) then
                   do j=this_block%jlo,this_block%jhi
                      do i=this_block%ilo,this_block%ihi
                         num_elim_gcells = num_elim_gcells + 1
                         ig = this_block%i_glob(i)
                         jg = this_block%j_glob(j)
                         gindex_elim(num_elim_gcells) = (jg-1)*nx_global + ig
                      end do
                   end do
                end if
             end if
          end do
       end do

       ! create a global index that includes both active and eliminated gridcells
       num_ice  = size(gindex_ice)
       num_elim = size(gindex_elim)
       allocate(gindex(num_elim + num_ice))
       do n = 1,num_ice
          gindex(n) = gindex_ice(n)
       end do
       do n = num_ice+1,num_ice+num_elim
          gindex(n) = gindex_elim(n-num_ice)
       end do

       deallocate(gindex_elim)

    else

       ! No eliminated land blocks
       num_ice = size(gindex_ice)
       allocate(gindex(num_ice))
       do n = 1,num_ice
          gindex(n) = gindex_ice(n)
       end do

    end if

    !---------------------------------------------------------------------------
    ! Create distGrid from global index array
    !---------------------------------------------------------------------------

    DistGrid = ESMF_DistGridCreate(arbSeqIndexList=gindex, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------------------------------------------
    ! Create the CICE mesh
    !---------------------------------------------------------------------------

    ! read in the mesh
    call NUOPC_CompAttributeGet(gcomp, name='mesh_ice', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    EMeshTemp = ESMF_MeshCreate(filename=trim(cvalue), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (my_task == master_task) then
       write(nu_diag,*)'mesh file for cice domain is ',trim(cvalue)
    end if

    ! recreate the mesh using the above distGrid
    EMesh = ESMF_MeshCreate(EMeshTemp, elementDistgrid=Distgrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! obtain mesh lats and lons
    call ESMF_MeshGet(Emesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(ownedElemCoords(spatialDim*numOwnedElements))
    allocate(lonMesh(numOwnedElements), latMesh(numOwnedElements))
    call ESMF_MeshGet(Emesh, ownedElemCoords=ownedElemCoords)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    do n = 1,numOwnedElements
       lonMesh(n) = ownedElemCoords(2*n-1)
       latMesh(n) = ownedElemCoords(2*n)
    end do

    call icepack_query_parameters(rad_to_deg_out=rad_to_deg)
    call icepack_warnings_flush(nu_diag)
    if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
        file=__FILE__, line=__LINE__)

    ! obtain internally generated cice lats and lons for error checks
    allocate(lon(lsize))
    allocate(lat(lsize))
    n = 0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       do j = jlo, jhi
          do i = ilo, ihi
             n = n+1
             lon(n) = tlon(i,j,iblk)*rad_to_deg
             lat(n) = tlat(i,j,iblk)*rad_to_deg
          enddo
       enddo
    enddo

    ! error check differences between internally generated lons and those read in
    do n = 1,lsize
       diff_lon = abs(lonMesh(n) - lon(n))
       if ( (diff_lon > 1.e2  .and. abs(diff_lon - 360_dbl_kind) > 1.e-1) .or.&
            (diff_lon > 1.e-3 .and. diff_lon < 1._dbl_kind) ) then
          !write(6,100)n,lonMesh(n),lon(n), diff_lon
100       format('ERROR: CICE  n, lonmesh(n), lon(n), diff_lon = ',i6,2(f21.13,3x),d21.5)
          !call abort_ice()
       end if
       if (abs(latMesh(n) - lat(n)) > 1.e-1) then
          !write(6,101)n,latMesh(n),lat(n), abs(latMesh(n)-lat(n))
101       format('ERROR: CICE n, latmesh(n), lat(n), diff_lat = ',i6,2(f21.13,3x),d21.5)
          !call abort_ice()
       end if
    end do

    ! deallocate memory
    deallocate(ownedElemCoords)
    deallocate(lon, lonMesh)
    deallocate(lat, latMesh)

    !-----------------------------------------------------------------
    ! Realize the actively coupled fields
    !-----------------------------------------------------------------

    call ice_realize_fields(gcomp, mesh=Emesh, &
         flds_scalar_name=flds_scalar_name, flds_scalar_num=flds_scalar_num, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !-----------------------------------------------------------------
    ! Prescribed ice initialization - first get compid
    !-----------------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='MCTID', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) compid  ! convert from string to integer
    else
       compid = 0
    end if
    call ice_prescribed_init(lmpicom, compid, gindex_ice)

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

    ! TODO (mvertens, 2018-12-21): fill in iceberg_prognostic as .false.
    if (debug_export > 0 .and. my_task==master_task) then
       call State_fldDebug(exportState, flds_scalar_name, 'cice_export:', &
            idate, sec, nu_diag, rc=rc)
    end if

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (dbug > 0) then
       call state_diagnose(exportState,subname//':ES',rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

    call t_stopf ('cice_init_total')

    deallocate(gindex_ice)
    deallocate(gindex)

    call flush_fileunit(nu_diag)

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
    character(*)   , parameter :: F00   = "('(ice_comp_nuopc) ',2a,i8,d21.14)"
    character(len=*),parameter :: subname=trim(modName)//':(ModelAdvance) '
    character(char_len_long)   :: msgString
    !--------------------------------

    rc = ESMF_SUCCESS
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

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_setLogUnit (nu_diag)

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
       write(nu_diag,F00) trim(subname),' cice istep, nextsw_cday = ',istep, nextsw_cday
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
    tod = sec
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
    else
       force_restart_now = .false.
    endif

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
            idate, sec, nu_diag, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    if (dbug > 0) then
       call state_diagnose(importState,subname//':IS',rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !--------------------------------
    ! Advance cice and timestep update
    !--------------------------------

    call CICE_Run()

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
            idate, sec, nu_diag, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    if (dbug > 0) then
       call state_diagnose(exportState,subname//':ES',rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! reset shr logging to my original values
    call shr_file_setLogUnit (shrlogunit)

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
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(*), parameter :: F00   = "('(ice_comp_nuopc) ',8a)"
    character(*), parameter :: F91   = "('(ice_comp_nuopc) ',73('-'))"
    character(len=*),parameter  :: subname=trim(modName)//':(ModelFinalize) '
    !--------------------------------

    !--------------------------------
    ! Finalize routine
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    if (my_task == master_task) then
       write(nu_diag,F91)
       write(nu_diag,F00) 'CICE: end of main integration loop'
       write(nu_diag,F91)
    end if

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

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
    logical                      :: lprint
    logical                      :: first_time = .true.
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

    if (trim(orb_mode) == trim(orb_variable_year)) then
       call ESMF_ClockGet(clock, CurrTime=CurrTime, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeGet(CurrTime, yy=year, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       orb_year = orb_iyear + (year - orb_iyear_align)
       lprint = mastertask
    else
       orb_year = orb_iyear
       if (first_time) then
          lprint = mastertask
       else
          lprint = .false.
       end if
    end if

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

  !===============================================================================


end module ice_comp_nuopc
