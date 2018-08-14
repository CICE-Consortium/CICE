module ice_comp_esmf

#ifdef ESMF_INTERFACE

!---------------------------------------------------------------------------
!BOP
!
! !MODULE: ice_comp_esmf
!
! !DESCRIPTION:
! CICE interface routine for the ccsm cpl7 mct system
!
! !USES:

  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod,  only : shr_sys_abort, shr_sys_flush
! use shr_mem_mod,  only : shr_get_memusage, shr_init_memusage
  use shr_file_mod, only : shr_file_getlogunit, shr_file_getloglevel,  &
		           shr_file_setloglevel, shr_file_setlogunit
  use mct_mod
#ifdef USE_ESMF_LIB
  use esmf
#else
  use esmf, only: ESMF_clock
#endif

  use seq_flds_mod
  use seq_infodata_mod,only : seq_infodata_start_type_cont, &
		              seq_infodata_start_type_brnch, seq_infodata_start_type_start
  use seq_timemgr_mod, only : seq_timemgr_eclockgetdata, &
                              seq_timemgr_restartalarmison, &
		              seq_timemgr_eclockdateinsync, &
                              seq_timemgr_stopalarmison
  use seq_comm_mct,    only : seq_comm_suffix, seq_comm_inst, seq_comm_name
  use perf_mod,        only : t_startf, t_stopf, t_barrierf
  use esmfshr_util_mod, only : esmfshr_util_StateArrayDestroy
  use esmfshr_util_mod, only : esmfshr_util_ArrayGetIndex
  use esmf2mct_mod,     only : esmf2mct_copy, esmf2mct_init

  use ice_cpl_indices
  use ice_import_export
  use ice_state,       only : aice
  use ice_domain_size, only : nx_global, ny_global, block_size_x, block_size_y, max_blocks
  use ice_domain,      only : nblocks, blocks_ice, halo_info, distrb_info
  use ice_blocks,      only : block, get_block, nx_block, ny_block
  use ice_grid,        only : tlon, tlat, tarea, tmask, anglet, hm, &
 		              grid_type, t2ugrid_vector, gridcpl_file, ocn_gridcell_frac
  use ice_constants,   only : c0, c1, spval_dbl, rad_to_deg, radius, secday
  use ice_communicate, only : my_task, master_task, MPI_COMM_ICE
  use ice_calendar,    only : istep, istep1, force_restart_now, write_ic,&
                              idate, idate0, mday, time, month, daycal,  &
		              sec, dt, dt_dyn, calendar,                 &
                              calendar_type, nextsw_cday, days_per_year, &
                              nyr, new_year, time2sec, year_init
  use icepack_orbital,     only : eccen, obliqr, lambm0, mvelpp
  use ice_timers

  use ice_kinds_mod,   only : int_kind, dbl_kind, char_len_long, log_kind
  use ice_boundary,    only : ice_HaloUpdate 
  use ice_scam,        only : scmlat, scmlon, single_column
  use ice_fileunits,   only : nu_diag, inst_index, inst_name, inst_suffix, &
                              release_all_fileunits
  use ice_prescribed_mod
  use ice_step_mod
  use ice_global_reductions
  use ice_broadcast
  use CICE_RunMod
  use ice_exit, only: abort_ice
  use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted

! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  public :: ice_register_esmf
  public :: ice_init_esmf
  public :: ice_run_esmf
  public :: ice_final_esmf
  private                              ! By default make data private
!
! ! PUBLIC DATA:
!
! !REVISION HISTORY:
! Author: Jacob Sewall, Mariana Vertenstein
!
!EOP
! !PRIVATE MEMBER FUNCTIONS:
  private :: ice_distgrid_esmf
  private :: ice_domain_esmf

!
! !PRIVATE VARIABLES

  integer (kind=int_kind) :: ICEID       
  type(mct_gGrid) :: dom_i
  type(mct_gsMap) :: gsMap_i

  !--- for coupling on other grid from gridcpl_file ---
  !type(mct_gsMap) :: gsMap_iloc  ! local gsmaps
  !type(mct_gGrid) :: dom_iloc                 ! local domain
  !type(mct_aVect) :: x2i_iloc, i2x_iloc
  !type(mct_rearr) :: rearr_ice2iloc
  !type(mct_rearr) :: rearr_iloc2ice
  !integer         :: nxcpl, nycpl  ! size of coupling grid
  logical         :: other_cplgrid    ! using different coupling grid

!=======================================================================

contains

!=======================================================================
subroutine ice_register_esmf(comp, rc)
    type(ESMF_GridComp)  :: comp
    integer, intent(out) :: rc
    character(len=*), parameter :: subname = '(ice_register_esmf)'

    rc = ESMF_SUCCESS
    print *, "In ice register routine"
    ! Register the callback routines.

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
      ice_init_esmf, phase=1, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, &
      ice_run_esmf, phase=1, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, &
      ice_final_esmf, phase=1, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

end subroutine

!=======================================================================
!BOP
!
! !IROUTINE: ice_init_esmf
!
! !INTERFACE:
  subroutine ice_init_esmf(comp, import_state, export_state, EClock, rc)
!
! !DESCRIPTION:
! Initialize thermodynamic ice model and obtain relevant atmospheric model
! arrays back from driver 
!
! !USES:

    use CICE_InitMod
    use ice_restart_shared, only: runid, runtype, restart_dir, restart_format
    use ice_history,        only: accum_hist
    use ice_history_shared, only: history_dir, history_file
    use icepack_intfc, only: tr_aero, tr_zaero
!
! !ARGUMENTS:
    type(ESMF_GridComp)          :: comp
    type(ESMF_State)             :: import_state
    type(ESMF_State)             :: export_state
    type(ESMF_Clock)             :: EClock
    integer, intent(out)         :: rc
!
! !LOCAL VARIABLES:
!
    type(ESMF_ArraySpec)   :: arrayspec
    type(ESMF_DistGrid)    :: distgrid
    type(ESMF_Array)       :: i2x, x2i, dom
    type(ESMF_VM)          :: vm

    integer                               :: lsize,lsize_loc
    integer                               :: xoff,yoff
    integer                               :: nxg,nyg
    integer                               :: k, iblk
 
    character(len=256) :: drvarchdir         ! driver archive directory
    character(len=32)  :: starttype          ! infodata start type
    integer            :: start_ymd          ! Start date (YYYYMMDD)
    integer            :: start_tod          ! start time of day (s)
    integer            :: curr_ymd           ! Current date (YYYYMMDD)
    integer            :: curr_tod           ! Current time of day (s)
    integer            :: ref_ymd            ! Reference date (YYYYMMDD)
    integer            :: ref_tod            ! reference time of day (s)
    integer            :: iyear              ! yyyy
    integer            :: nyrp               ! yyyy
    integer            :: dtime              ! time step
    integer            :: shrlogunit,shrloglev ! old values
    integer            :: iam,ierr
    integer            :: lbnum
    integer            :: daycal(13)  !number of cumulative days per month
    integer            :: nleaps      ! number of leap days before current year
    integer            :: mpicom_loc, mpicom_vm, gsize
    integer            :: nfields
    logical (kind=log_kind) :: atm_aero
    real(r8), pointer      :: fptr(:,:)
    character(ESMF_MAXSTR) :: convCIM, purpComp
    real(r8) :: mrss, mrss0,msize,msize0
    character(len=*), parameter :: subname = '(ice_init_esmf)'
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!EOP
!-----------------------------------------------------------------------

    call t_barrierf('cice_init_total_BARRIER',MPI_COMM_ICE)
    call t_startf ('cice_init_total')

    rc = ESMF_SUCCESS

    !--------------------------------------------------------------------------
    ! Determine attribute vector indices
    !--------------------------------------------------------------------------

    call ice_cpl_indices_set()

    ! duplicate the mpi communicator from the current VM 
    call ESMF_VMGetCurrent(vm, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_VMGet(vm, mpiCommunicator=mpicom_vm, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call MPI_Comm_dup(mpicom_vm, mpicom_loc, rc)
    if(rc /= 0) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Initialize cice id
   
    call ESMF_AttributeGet(export_state, name="ID", value=ICEID, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Determine time of next atmospheric shortwave calculation
    call ESMF_AttributeGet(export_state, name="nextsw_cday", value=nextsw_cday, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Determine if aerosols are coming through the coupler
    call ESMF_AttributeGet(export_state, name="atm_aero", value=atm_aero, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="ID", value=ICEID, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Determine orbital parameters
    call ESMF_AttributeGet(export_state, name="orb_eccen", value=eccen, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="orb_obliqr", value=obliqr, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="orb_lambm0", value=lambm0, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="orb_mvelpp", value=mvelpp, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !   call shr_init_memusage()

    !---------------------------------------------------------------------------
    ! use infodata to determine type of run
    !---------------------------------------------------------------------------

    ! Preset single column values

    single_column = .false.
    scmlat = -999.
    scmlon = -999.

    call ESMF_AttributeGet(export_state, name="case_name", value=runid, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="single_column", value=single_column, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="scmlat", value=scmlat, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="scmlon", value=scmlon, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="start_type", value=starttype, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
       runtype = "initial"
    else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
       runtype = "continue"
    else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
       runtype = "continue"
    else
       write(nu_diag,*) trim(subname),' ERROR: unknown starttype'
       call shr_sys_abort()
    end if

    ! Set nextsw_cday to -1 for continue and branch runs.

    if (trim(runtype) /= 'initial') nextsw_cday = -1

    !=============================================================
    ! Set ice dtime to ice coupling frequency
    !=============================================================

    call seq_timemgr_EClockGetData(EClock, dtime=dtime, calendar=calendar_type)
    dt = real(dtime)

    !=============================================================
    ! Initialize cice because grid information is needed for
    ! creation of gsmap_i.  cice_init also sets time manager info
    !=============================================================

    inst_name   = seq_comm_name(ICEID)
    inst_index  = seq_comm_inst(ICEID)
    inst_suffix = seq_comm_suffix(ICEID)

    write(nu_diag,*) trim(subname),'inst_name   = ',trim(inst_name)
    write(nu_diag,*) trim(subname),'inst_index  = ',inst_index
    write(nu_diag,*) trim(subname),'inst_suffix = ',trim(inst_suffix)

    call t_startf ('cice_init')
    call cice_init( mpicom_loc )
    call t_stopf ('cice_init')

    !---------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !---------------------------------------------------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (nu_diag)

    !---------------------------------------------------------------------------
    ! use EClock to reset calendar information on initial start
    !---------------------------------------------------------------------------

    ! - on restart run 
    !   - istep0, time and time_forc are read from restart file
    !   - istep1 is set to istep0
    !   - idate is determined from time via the call to calendar (see below)
    ! - on initial run 
    !   - iyear, month and mday obtained from sync clock
    !   - time determined from iyear, month and mday
    !   - istep0 and istep1 are set to 0 

    call seq_timemgr_EClockGetData(EClock,               &
         start_ymd=start_ymd, start_tod=start_tod,       &
         curr_ymd=curr_ymd,   curr_tod=curr_tod,         &
         ref_ymd=ref_ymd,     ref_tod=ref_tod)

    if (runtype == 'initial') then
       if (ref_ymd /= start_ymd .or. ref_tod /= start_tod) then
          if (my_task == master_task) then
             write(nu_diag,*) trim(subname),' ref_ymd ',ref_ymd, &
                  ' must equal start_ymd ',start_ymd
             write(nu_diag,*) trim(subname),' ref_ymd ',ref_tod, &
                  ' must equal start_ymd ',start_tod
          end if
       end if

       if (my_task == master_task) then
          write(nu_diag,*) trim(subname)//' idate from sync clock = ', &
               start_ymd
          write(nu_diag,*) trim(subname)//' tod from sync clock = ', &
               start_tod
          write(nu_diag,*) trim(subname)//' resetting idate to match sync clock'
       end if

       idate0 = curr_ymd - (year_init*10000)
       idate = curr_ymd - (year_init*10000)

       if (idate < 0) then
          write(nu_diag,*) trim(subname),' ERROR curr_ymd,year_init =',curr_ymd,year_init
          write(nu_diag,*) trim(subname),' ERROR idate lt zero',idate
          call shr_sys_abort(subname//' :: ERROR idate lt zero')
       endif

       iyear = (idate/10000)                     ! integer year of basedate
       month = (idate-iyear*10000)/100           ! integer month of basedate
       mday  =  idate-iyear*10000-month*100      ! day of month of basedate

       if (my_task == master_task) then
          write(nu_diag,*) trim(subname),' curr_ymd = ',curr_ymd
          write(nu_diag,*) trim(subname),' cice year_init = ',year_init
          write(nu_diag,*) trim(subname),' cice start date = ',idate
          write(nu_diag,*) trim(subname),' cice start ymds = ',iyear,month,mday,start_tod
       endif

       call time2sec(iyear,month,mday,time)
       time = time+start_tod

       call shr_sys_flush(nu_diag)
    end if

    call calendar(time)     ! update calendar info
    if (write_ic) call accum_hist(dt) ! write initial conditions
 
    !---------------------------------------------------------------------------
    ! Initialize MCT attribute vectors and indices
    !---------------------------------------------------------------------------

    call t_startf ('cice_esmf_init')

    ! check gridcpl_file

    other_cplgrid = .true.
    if (trim(gridcpl_file) == 'unknown_gridcpl_file') then
       other_cplgrid = .false.
    else
       call shr_sys_abort(subname//' :: gridcpl_file not supported with esmf interface yet')
    endif

    !-----------------------------------------
    ! Initialize distgrid and gsmap_i
    ! (gsmap_i is needed for prescribed_ice)
    !-----------------------------------------

    distgrid = ice_distgrid_esmf(gsize)

    call ESMF_AttributeSet(export_state, name="gsize", value=gsize, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !-----------------------------------------
    !  Set arrayspec for dom, l2x and x2l
    !-----------------------------------------
    
    call ESMF_ArraySpecSet(arrayspec, rank=2, typekind=ESMF_TYPEKIND_R8, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !-----------------------------------------
    ! Create dom 
    !-----------------------------------------

    nfields = shr_string_listGetNum(trim(seq_flds_dom_fields))

    dom = ESMF_ArrayCreate(distgrid=distgrid, arrayspec=arrayspec, distgridToArrayMap=(/2/), &
         undistLBound=(/1/), undistUBound=(/nfields/), name="domain", rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(dom, name="mct_names", value=trim(seq_flds_dom_fields), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Set values of dom 
    call ice_domain_esmf(dom)

    !----------------------------------------- 
    !  Create i2x 
    !-----------------------------------------

    ! 1d undistributed index of fields, 2d is packed data

    nfields = shr_string_listGetNum(trim(seq_flds_i2x_fields))

    i2x = ESMF_ArrayCreate(distgrid=distgrid, arrayspec=arrayspec, distgridToArrayMap=(/2/), &
         undistLBound=(/1/), undistUBound=(/nfields/), name="d2x", rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(i2x, name="mct_names", value=trim(seq_flds_i2x_fields), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
 
    !----------------------------------------- 
    !  Create x2i 
    !-----------------------------------------

    nfields = shr_string_listGetNum(trim(seq_flds_x2i_fields))

    x2i = ESMF_ArrayCreate(distgrid=distgrid, arrayspec=arrayspec, distgridToArrayMap=(/2/), &
         undistLBound=(/1/), undistUBound=(/nfields/), name="x2d", rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(x2i, name="mct_names", value=trim(seq_flds_x2i_fields), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !----------------------------------------- 
    ! Add esmf arrays to import and export state 
    !-----------------------------------------
 
    call ESMF_StateAdd(export_state, (/dom/), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_StateAdd(export_state, (/i2x/), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
 
    call ESMF_StateAdd(import_state, (/x2i/), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !-----------------------------------------------------------------
    ! Second phase of prescribed ice initialization
    ! Need to create gsmap_i and dom_i (module variables)
    !-----------------------------------------------------------------

    call esmf2mct_init(distgrid, ICEID, gsmap_i, MPI_COMM_ICE, gsize=gsize, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call esmf2mct_init(dom, dom_i, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call esmf2mct_copy(dom, dom_i%data, rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !-----------------------------------------------------------------
    ! Prescribed ice initialization
    !-----------------------------------------------------------------

    if (other_cplgrid) then
       call shr_sys_abort(subname//' :: gridcpl_file not supported with esmf interface yet')
    else
       call ice_prescribed_init(ICEID, gsmap_i, dom_i)
    endif

    !---------------------------------------------------------------------------
    ! create ice export state
    !---------------------------------------------------------------------------

    call ESMF_ArrayGet(i2x, localDe=0, farrayPtr=fptr, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    if (other_cplgrid) then
       call shr_sys_abort(subname//' :: gridcpl_file not supported with esmf interface yet')
    else
       call ice_export (fptr)  !Send initial state to driver
    endif

    call ESMF_AttributeSet(export_state, name="ice_prognostic", value=.true., rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

! tcraig: iceberg_prognostic is false by default in cesm1.3
! not explicitly setting it here makes cice5 work in cesm1.1
!    call ESMF_AttributeSet(export_state, name="iceberg_prognostic", value=.false., rc=rc)
!    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(export_state, name="ice_nx", value=nx_global, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(export_state, name="ice_ny", value=ny_global, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

#ifdef USE_ESMF_METADATA
    convCIM  = "CIM"
    purpComp = "Model Component Simulation Description"

    call ESMF_AttributeAdd(comp,  &
                           convention=convCIM, purpose=purpComp, rc=rc)

    call ESMF_AttributeSet(comp, "ShortName", "CICE", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "LongName", &
                           "Community Ice CodE", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "Description", &
                  "CICE4 is the latest version of the Los Alamos Sea Ice " // &
                  "Model, sometimes referred to as the Community Ice " // &
                  "CodE.  It is the result of a community effort to  " // &
                  "develop a portable, efficient sea ice model that can " // &
                  "be run coupled in a global climate model or uncoupled " // &
                  "as a stand-alone ice model.  It has been released as " // &
                  "the sea ice component of the Community Earth System " // &
                  "Model (CESM), a fully-coupled global climate model " // &
                  "that provides simulations of the earths past, present " // &
                  "and future climate states.  CICE4 is supported on " // &
                  "high- and low-resolution Greenland Pole and tripole " // &
                  "grids, which are identical to those used by the " // &
                  "Parallel Ocean Program (POP) ocean model.  The high " // &
                  "resolution version is best suited for simulating " // &
                  "present-day and future climate scenarios while the low " // &
                  "resolution option is used for paleoclimate simulations " // &
                  "and debugging.", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ReleaseDate", "2010", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ModelType", "Sea Ice", &
                           convention=convCIM, purpose=purpComp, rc=rc)

!    call ESMF_AttributeSet(comp, "Name", "someone", &
!                           convention=convCIM, purpose=purpComp, rc=rc)
!    call ESMF_AttributeSet(comp, "EmailAddress", &
!                           "someone@someplace", &
!                           convention=convCIM, purpose=purpComp, rc=rc)
!    call ESMF_AttributeSet(comp, "ResponsiblePartyRole", "contact", &
!                           convention=convCIM, purpose=purpComp, rc=rc)

#endif

    call t_stopf ('cice_esmf_init')

    ! Error check
    if ((tr_aero .and. .not. atm_aero) .or. (tr_zaero .and. .not. atm_aero)) then
       write(nu_diag,*) 'ice_import ERROR: atm_aero must be set for tr_aero or tr_zaero' 
       call shr_sys_abort()
    end if

    !---------------------------------------------------------------------------
    ! Reset shr logging to original values
    !---------------------------------------------------------------------------

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

    call t_stopf ('cice_init_total')

    !   call ice_timer_stop(timer_total) ! time entire run
    !   call shr_get_memusage(msize,mrss)
    !   call shr_mpi_max(mrss, mrss0, MPI_COMM_ICE,trim(subname)//' mrss0')
    !   call shr_mpi_max(msize,msize0,MPI_COMM_ICE,trim(subname)//' msize0')
    !   if(my_task == 0) then
    !   write(shrlogunit,105) trim(subname)//' memory_write: model date = ',start_ymd,start_tod, &
    !           ' memory = ',msize0,' MB (highwater)    ',mrss0,' MB (usage)'
    !   endif
 
  105  format( A, 2i8, A, f10.2, A, f10.2, A)

  end subroutine ice_init_esmf

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ice_run_esmf
!
! !INTERFACE:
  subroutine ice_run_esmf(comp, import_state, export_state, EClock, rc)
!
! !DESCRIPTION:
! Run thermodynamic CICE
!
! !USES:
    use ice_history
    use ice_restart
    use ice_diagnostics
    use ice_restoring , only : restore_ice, ice_HaloRestore
    use icepack_shortwave , only : init_shortwave

! !ARGUMENTS:
    type(ESMF_GridComp)          :: comp
    type(ESMF_State)             :: import_state
    type(ESMF_State)             :: export_state
    type(ESMF_Clock)             :: EClock
    integer, intent(out)         :: rc

! !LOCAL VARIABLES:
    integer :: k             ! index
    logical :: stop_now      ! .true. ==> stop at the end of this run phase
    integer :: ymd           ! Current date (YYYYMMDD)
    integer :: tod           ! Current time of day (sec)
    integer :: yr_sync       ! Sync current year
    integer :: mon_sync      ! Sync current month
    integer :: day_sync      ! Sync current day
    integer :: tod_sync      ! Sync current time of day (sec)
    integer :: ymd_sync      ! Current year of sync clock
    integer :: curr_ymd           ! Current date (YYYYMMDD)
    integer :: curr_tod           ! Current time of day (s)
    integer :: shrlogunit,shrloglev ! old values
    integer :: lbnum
    integer :: n, nyrp
    type(ESMF_Array)             :: i2x, x2i
    real(R8), pointer            :: fptr(:,:)
    character(len=*), parameter :: subname = '(ice_run_esmf)'

    real(r8) :: mrss, mrss0,msize,msize0
    logical, save :: first_time = .true.

!
! !REVISION HISTORY:
! Author: Jacob Sewall, Mariana Vertenstein
!
!EOP
!---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call ice_timer_start(timer_total) ! time entire run
    call t_barrierf('cice_run_total_BARRIER',MPI_COMM_ICE)
    call t_startf ('cice_run_total')

    !---------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !---------------------------------------------------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (nu_diag)
   
    ! Determine time of next atmospheric shortwave calculation

    call ESMF_AttributeGet(export_state, name="nextsw_cday", value=nextsw_cday, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Determine orbital parameters
    call ESMF_AttributeGet(export_state, name="orb_eccen", value=eccen, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="orb_obliqr", value=obliqr, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="orb_lambm0", value=lambm0, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="orb_mvelpp", value=mvelpp, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Get clock information
    call seq_timemgr_EClockGetData(EClock,               &
         curr_ymd=curr_ymd, curr_tod=curr_tod)

    force_restart_now = seq_timemgr_RestartAlarmIsOn(EClock)

    if (calendar_type .eq. "GREGORIAN") then 	
       nyrp = nyr
       nyr = (curr_ymd/10000)+1           ! integer year of basedate
       if (nyr /= nyrp) then
          new_year = .true.
       else
          new_year = .false.
       end if
    end if

    !-------------------------------------------------------------------
    ! get import state
    !-------------------------------------------------------------------
    
    call t_barrierf('cice_run_import_BARRIER',MPI_COMM_ICE)
    call t_startf ('cice_run_import')

    call ESMF_StateGet(import_state, itemName="x2d", array=x2i, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_ArrayGet(x2i, localDe=0, farrayPtr=fptr, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ice_timer_start(timer_cplrecv)
    if (other_cplgrid) then
       call shr_sys_abort(subname//' :: gridcpl_file not supported with esmf interface yet')
    else
       call ice_import( fptr )
    endif
    call ice_timer_stop(timer_cplrecv)
    call t_stopf ('cice_run_import')
 
    !--------------------------------------------------------------------
    ! timestep update
    !--------------------------------------------------------------------

    call CICE_Run()

    !-----------------------------------------------------------------
    ! send export state to driver 
    !-----------------------------------------------------------------
    
    call t_barrierf('cice_run_export_BARRIER',MPI_COMM_ICE)
    call t_startf ('cice_run_export')
    call ice_timer_start(timer_cplsend)

    call ESMF_StateGet(export_state, itemName="d2x", array=i2x, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_ArrayGet(i2x, localDe=0, farrayPtr=fptr, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    if (other_cplgrid) then
       call shr_sys_abort(subname//' :: gridcpl_file not supported with esmf interface yet')
    else
       call ice_export (fptr)
    endif
    call ice_timer_stop(timer_cplsend)
    call t_stopf ('cice_run_export')
    
    !--------------------------------------------------------------------
    ! check that internal clock is in sync with master clock
    !--------------------------------------------------------------------

    tod = sec
    ymd = idate
    if (.not. seq_timemgr_EClockDateInSync( EClock, ymd, tod )) then
       call seq_timemgr_EClockGetData( EClock, curr_ymd=ymd_sync, &
          curr_tod=tod_sync )
       write(nu_diag,*)' cice ymd=',ymd     ,'  cice tod= ',tod
       write(nu_diag,*)' sync ymd=',ymd_sync,'  sync tod= ',tod_sync
       call shr_sys_abort( SubName// &
          ":: Internal sea-ice clock not in sync with Sync Clock")
    end if
   
    ! reset shr logging to my original values

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

    !-------------------------------------------------------------------
    ! stop timers and print timer info
    !-------------------------------------------------------------------
    ! Need to have this logic here instead of in ice_final_esmf since 
    ! the ice_final_esmf.F90 will still be called even in aqua-planet mode
    ! Could put this logic in the driver - but it seems easier here 

    ! Need to stop this at the end of every run phase in a coupled run.
    call ice_timer_stop(timer_total)        ! stop timing

    stop_now = seq_timemgr_StopAlarmIsOn( EClock )
    if (stop_now) then
       call ice_timer_print_all(stats=.true.) ! print timing information
       call release_all_fileunits
    end if
    
!   if(tod == 0) then
!      call shr_get_memusage(msize,mrss)
!      call shr_mpi_max(mrss, mrss0, MPI_COMM_ICE,trim(subname)//' mrss0')
!      call shr_mpi_max(msize,msize0,MPI_COMM_ICE,trim(subname)//' msize0')
!      if(my_task == 0 ) then
!          write(shrlogunit,105) trim(subname)//': memory_write: model date = ',ymd,tod, &
!               ' memory = ',msize0,' MB (highwater)    ',mrss0,' MB (usage)'
!      endif
!   endif
    call t_stopf ('cice_run_total')
 
  105  format( A, 2i8, A, f10.2, A, f10.2, A)

  end subroutine ice_run_esmf

!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ice_final_esmf
!
! !INTERFACE:
  subroutine ice_final_esmf(comp, import_state, export_state, EClock, rc)
!
! !DESCRIPTION:
! Finalize CICE
!
! !USES:
!
!------------------------------------------------------------------------------
!BOP
!
! !ARGUMENTS:
   type(ESMF_GridComp)          :: comp
   type(ESMF_State)             :: import_state
   type(ESMF_State)             :: export_state
   type(ESMF_Clock)             :: EClock
   integer, intent(out)         :: rc

   character(len=*), parameter :: subname = '(ice_final_esmf)'
!
! !REVISION HISTORY:
!
!EOP
!---------------------------------------------------------------------------

    ! Note that restart for final timestep was written in run phase.
    rc = ESMF_SUCCESS

    ! Destroy ESMF objects

    call esmfshr_util_StateArrayDestroy(export_state,"d2x",rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call esmfshr_util_StateArrayDestroy(export_state,"domain",rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call esmfshr_util_StateArrayDestroy(import_state,"x2d",rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

  end subroutine ice_final_esmf

!===============================================================================

  function ice_distgrid_esmf(gsize)

    implicit none
    !-------------------------------------------------------------------
    !
    ! Arguments
    !
    integer, intent(out)    :: gsize
    !
    ! Return
    type(esmf_distgrid)     :: ice_distgrid_esmf
    !
    ! Local variables
    !
    integer,allocatable :: gindex(:)
    integer     :: lat
    integer     :: lon
    integer     :: i, j, iblk, n, gi
    integer     :: lsize
    integer     :: rc
    integer     :: ilo, ihi, jlo, jhi ! beginning and end of physical domain
    type(block) :: this_block         ! block information for current block
    character(len=*), parameter :: subname = '(ice_distgrid_esmf)'
    !-------------------------------------------------------------------

    ! number the local grid

    n=0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)         
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       
       do j = jlo, jhi
          do i = ilo, ihi
             n = n+1
          enddo !i
       enddo    !j
    enddo        !iblk
    lsize = n

    ! not valid for padded decomps
    !    lsize = block_size_x*block_size_y*nblocks
    gsize = nx_global*ny_global

    allocate(gindex(lsize))
    n=0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)         
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       
       do j = jlo, jhi
          do i = ilo, ihi
             n = n+1
             lon = this_block%i_glob(i)
             lat = this_block%j_glob(j)
             gi = (lat-1)*nx_global + lon
             gindex(n) = gi
          enddo !i
       enddo    !j
    enddo        !iblk
   
    ice_distgrid_esmf = ESMF_DistGridCreate(arbSeqIndexList=gindex, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    deallocate(gindex)

  end function ice_DistGrid_esmf

!====================================================================================

  subroutine ice_domain_esmf( dom )

    !-------------------------------------------------------------------
    !
    ! Arguments
    !
    type(ESMF_Array), intent(inout) :: dom
    !
    ! Local Variables
    !
    integer           :: i, j, iblk, n               ! indices
    integer           :: ilo, ihi, jlo, jhi          ! beginning and end of physical domain
    integer           :: klon,klat,karea,kmask,kfrac ! domain fields
    type(block)       :: this_block                  ! block information for current block
    real(R8), pointer :: fptr (:,:)
    integer           :: rc
    character(len=*), parameter :: subname = '(ice_domain_esmf)'
    !-------------------------------------------------------------------

    call ESMF_ArrayGet(dom, localDe=0, farrayPtr=fptr, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Fill in correct values for domain components
    klon  = esmfshr_util_ArrayGetIndex(dom,'lon ',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    klat  = esmfshr_util_ArrayGetIndex(dom,'lat ',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    karea = esmfshr_util_ArrayGetIndex(dom,'area',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    kmask = esmfshr_util_ArrayGetIndex(dom,'mask',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    kfrac = esmfshr_util_ArrayGetIndex(dom,'frac',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    fptr(:,:) = -9999.0_R8
    n=0
    do iblk = 1, nblocks
       this_block = get_block(blocks_ice(iblk),iblk)         
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       
       do j = jlo, jhi
       do i = ilo, ihi
          n = n+1
          fptr(klon, n)  = TLON(i,j,iblk)*rad_to_deg 
          fptr(klat, n)  = TLAT(i,j,iblk)*rad_to_deg 
          fptr(karea, n) = tarea(i,j,iblk)/(radius*radius)
          fptr(kmask, n) = real(nint(hm(i,j,iblk)),kind=dbl_kind)
          if (trim(grid_type) == 'latlon') then
             fptr(kfrac, n) = ocn_gridcell_frac(i,j,iblk)
          else
             fptr(kfrac, n) = real(nint(hm(i,j,iblk)),kind=dbl_kind)
          end if
       enddo    !i
       enddo    !j
    enddo       !iblk

  end subroutine ice_domain_esmf

!=======================================================================
#endif

end module ice_comp_esmf

