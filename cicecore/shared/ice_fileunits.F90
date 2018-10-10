!  SVN:$Id: ice_fileunits.F90 1228 2017-05-23 21:33:34Z tcraig $
!=======================================================================
!
!  This module contains an I/O unit manager for tracking, assigning
!  and reserving I/O unit numbers.
!
!  There are three reserved I/O units set as parameters in this
!  module.  The default units for standard input (stdin), standard
!  output (stdout) and standard error (stderr).  These are currently
!  set as units 5,6,6, respectively as that is the most commonly
!  used among vendors. However, the user may change these if those
!  default units are conflicting with other models or if the
!  vendor is using different values.
!
!  The maximum number of I/O units per node is currently set by
!  the parameter ice\_IOMaxUnit.
!
! author: Elizabeth C. Hunke, LANL
! 2006: ECH converted to free source form (F90)
! 2007: ECH added dynamic file units, modified from POP_IOUnitsMod.F90

      module ice_fileunits

      use ice_kinds_mod
#ifdef CESMCOUPLED
      use shr_file_mod, only : shr_file_getunit, shr_file_freeunit
#endif

      implicit none
      private
      public :: init_fileunits, get_fileunit, flush_fileunit, &
                release_fileunit, release_all_fileunits

      character (len=char_len), public :: &
         diag_type               ! 'stdout' or 'file'

      logical (log_kind), public :: &
         bfbflag                 ! logical for bit-for-bit computations

      integer (kind=int_kind), public :: &
         nu_grid       , &  ! grid file
         nu_kmt        , &  ! land mask file
         nu_nml        , &  ! namelist input file
         nu_forcing    , &  ! forcing data file
         nu_dump       , &  ! dump file for restarting
         nu_restart    , &  ! restart input file
         nu_dump_age   , &  ! dump file for restarting ice age tracer
         nu_restart_age, &  ! restart input file for ice age tracer
         nu_dump_FY    , &  ! dump file for restarting first-year area tracer
         nu_restart_FY , &  ! restart input file for first-year area tracer
         nu_dump_lvl   , &  ! dump file for restarting level ice tracers
         nu_restart_lvl, &  ! restart input file for level ice tracers
         nu_dump_pond  , &  ! dump file for restarting melt pond tracer
         nu_restart_pond,&  ! restart input file for melt pond tracer
         nu_dump_aero  , &  ! dump file for restarting aerosol tracer
         nu_restart_aero,&  ! restart input file for aerosol tracer
         nu_dump_bgc   , &  ! dump file for restarting bgc
         nu_restart_bgc, &  ! restart input file for bgc
         nu_dump_hbrine, &  ! dump file for restarting hbrine
         nu_restart_hbrine, &  ! restart input file for hbrine
         nu_dump_eap   , &  ! dump file for restarting eap dynamics
         nu_restart_eap, &  ! restart input file for eap dynamics
         nu_rst_pointer, &  ! pointer to latest restart file
         nu_history    , &  ! binary history output file
         nu_hdr        , &  ! header file for binary history output
         nu_diag            ! diagnostics output file

      character (32), public :: &
         nml_filename = 'ice_in' ! namelist input file name

      integer (kind=int_kind), parameter, public :: &
         ice_stdin  =  5, & ! reserved unit for standard input
         ice_stdout =  6, & ! reserved unit for standard output
         ice_stderr =  6    ! reserved unit for standard error

      integer (kind=int_kind), public :: &
         ice_IOUnitsMinUnit = 11, & ! do not use unit numbers below 
         ice_IOUnitsMaxUnit = 99    ! or above, set by setup_nml

      logical (kind=log_kind), dimension(:), allocatable :: &
         ice_IOUnitsInUse   ! flag=.true. if unit currently open

      ! instance control
      integer (kind=int_kind), public :: inst_index
      character(len=16)      , public :: inst_name
      character(len=16)      , public :: inst_suffix

!=======================================================================

      contains

!=======================================================================

!  This routine grabs needed unit numbers. 
!  nu_diag is set to 6 (stdout) but may be reset later by the namelist. 
!  nu_nml is obtained separately.

      subroutine init_fileunits

         character(len=*),parameter :: subname='(init_fileunits)'

         nu_diag = ice_stdout  ! default

         allocate(ice_IOUnitsInUse(ice_IOUnitsMaxUnit))
         ice_IOUnitsInUse = .false.
         ice_IOUnitsInUse(ice_stdin)  = .true. ! reserve unit 5
         ice_IOUnitsInUse(ice_stdout) = .true. ! reserve unit 6
         ice_IOUnitsInUse(ice_stderr) = .true.

         call get_fileunit(nu_grid)
         call get_fileunit(nu_kmt)
         call get_fileunit(nu_forcing)
         call get_fileunit(nu_dump)
         call get_fileunit(nu_restart)
         call get_fileunit(nu_dump_age)
         call get_fileunit(nu_restart_age)
         call get_fileunit(nu_dump_FY)
         call get_fileunit(nu_restart_FY)
         call get_fileunit(nu_dump_lvl)
         call get_fileunit(nu_restart_lvl)
         call get_fileunit(nu_dump_pond)
         call get_fileunit(nu_restart_pond)
         call get_fileunit(nu_dump_aero)
         call get_fileunit(nu_restart_aero)
         call get_fileunit(nu_dump_bgc)
         call get_fileunit(nu_restart_bgc)
         call get_fileunit(nu_dump_hbrine)
         call get_fileunit(nu_restart_hbrine)
         call get_fileunit(nu_dump_eap)
         call get_fileunit(nu_restart_eap)
         call get_fileunit(nu_rst_pointer)
         call get_fileunit(nu_history)
         call get_fileunit(nu_hdr)

      end subroutine init_fileunits

!=======================================================================

!  This routine returns the next available I/O unit and marks it as
!  in use to prevent any later use.
!  Note that {\em all} processors must call this routine even if only
!  the master task is doing the I/O.  This is necessary insure that
!  the units remain synchronized for other parallel I/O functions.

      subroutine get_fileunit(iunit)

         integer (kind=int_kind), intent(out) :: &
            iunit                     ! next free I/O unit

         ! local variables

#ifndef CESMCOUPLED
         integer (kind=int_kind) :: n  ! dummy loop index
         logical (kind=log_kind) :: alreadyInUse
#endif

         character(len=*),parameter :: subname='(get_fileunit)'

#ifdef CESMCOUPLED
         iunit = shr_file_getUnit()
#else

         srch_units: do n=ice_IOUnitsMinUnit, ice_IOUnitsMaxUnit
            if (.not. ice_IOUnitsInUse(n)) then   ! I found one, I found one

               !*** make sure not in use by library or calling routines
               INQUIRE (unit=n,OPENED=alreadyInUse)

               if (.not. alreadyInUse) then
                  iunit = n        ! return the free unit number
                  ice_IOUnitsInUse(iunit) = .true.  ! mark iunit as being in use
                  exit srch_units
               else
                  !*** if inquire shows this unit in use, mark it as
                  !***    in use to prevent further queries
                  ice_IOUnitsInUse(n) = .true.
               endif
            endif
         end do srch_units

         if (iunit > ice_IOUnitsMaxUnit) stop 'ice_IOUnitsGet: No free units'

#endif

      end subroutine get_fileunit

!=======================================================================

!  This routine releases unit numbers at the end of a run. 

      subroutine release_all_fileunits

         character(len=*),parameter :: subname='(release_all_fileunits)'

         call release_fileunit(nu_grid)
         call release_fileunit(nu_kmt)
         call release_fileunit(nu_forcing)
         call release_fileunit(nu_dump)
         call release_fileunit(nu_restart)
         call release_fileunit(nu_dump_age)
         call release_fileunit(nu_restart_age)
         call release_fileunit(nu_dump_FY)
         call release_fileunit(nu_restart_FY)
         call release_fileunit(nu_dump_lvl)
         call release_fileunit(nu_restart_lvl)
         call release_fileunit(nu_dump_pond)
         call release_fileunit(nu_restart_pond)
         call release_fileunit(nu_dump_aero)
         call release_fileunit(nu_restart_aero)
         call release_fileunit(nu_dump_bgc)
         call release_fileunit(nu_restart_bgc)
         call release_fileunit(nu_dump_hbrine)
         call release_fileunit(nu_restart_hbrine)
         call release_fileunit(nu_dump_eap)
         call release_fileunit(nu_restart_eap)
         call release_fileunit(nu_rst_pointer)
         call release_fileunit(nu_history)
         call release_fileunit(nu_hdr)
         if (nu_diag /= ice_stdout) call release_fileunit(nu_diag)

      end subroutine release_all_fileunits

!=======================================================================

!  This routine releases an I/O unit (marks it as available).
!  Note that {\em all} processors must call this routine even if only
!  the master task is doing the I/O.  This is necessary insure that
!  the units remain synchronized for other parallel I/O functions.

      subroutine release_fileunit(iunit)

         integer (kind=int_kind), intent(in) :: &
            iunit                    ! I/O unit to be released

         character(len=*),parameter :: subname='(release_fileunit)'

#ifdef CESMCOUPLED
         call shr_file_freeUnit(iunit)
#else
!  check for proper unit number
         if (iunit < 1 .or. iunit > ice_IOUnitsMaxUnit) then
            stop 'release_fileunit: bad unit'
         endif

!  mark the unit as not in use
         ice_IOUnitsInUse(iunit) = .false.  !  that was easy...
#endif

      end subroutine release_fileunit

!=======================================================================


!  This routine enables a user to flush the output from an IO unit
!  (typically stdout) to force output when the system is buffering
!  such output.  Because this system function is system dependent,
!  we only support this wrapper and users are welcome to insert the
!  code relevant to their local machine.  In the case where the CESM
!  libraries are available, the shared routine for sys flush can be
!  used (and is provided here under a preprocessor option).

      subroutine flush_fileunit(iunit)

#ifdef CESMCOUPLED
         use shr_sys_mod, only : shr_sys_flush
#endif

         integer (kind=int_kind), intent(in) :: &
            iunit                    ! I/O unit to be flushed

         character(len=*),parameter :: subname='(flush_fileunit)'

!-----------------------------------------------------------------------
!
!  insert your system code here
!
!-----------------------------------------------------------------------

#ifdef CESMCOUPLED
         call shr_sys_flush(iunit)
#else
#ifndef NO_F2003
         flush(iunit)
#else
! Place holder for old call.
#endif
#endif

      end subroutine flush_fileunit

!=======================================================================

      end module ice_fileunits

!=======================================================================
