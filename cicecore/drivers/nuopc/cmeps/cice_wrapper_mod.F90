module cice_wrapper_mod

#ifdef CESMCOUPLED
  use perf_mod      , only : t_startf, t_stopf, t_barrierf
  use shr_file_mod  , only : shr_file_getlogunit, shr_file_setlogunit

  use ice_kinds_mod , only : dbl_kind, int_kind, char_len, char_len_long

  implicit none

  real(dbl_kind) :: wtime = 0.0
contains
  ! Define stub routines that do nothing - they are just here to avoid
  ! having cppdefs in the main program
  subroutine ufs_settimer(timevalue)
    real(dbl_kind),    intent(inout) :: timevalue
  end subroutine ufs_settimer
  subroutine ufs_logtimer(nunit,elapsedsecs,string,runtimelog,time0)
    integer,           intent(in)  :: nunit
    integer(int_kind), intent(in)  :: elapsedsecs
    character(len=*),  intent(in)  :: string
    logical,           intent(in)  :: runtimelog
    real(dbl_kind),    intent(in)  :: time0
  end subroutine ufs_logtimer
  subroutine ufs_file_setLogUnit(filename,nunit,runtimelog)
    character(len=*),  intent(in)  :: filename
    logical,           intent(in)  :: runtimelog
    integer,           intent(out) :: nunit
  end subroutine ufs_file_setLogUnit
  subroutine ufs_logfhour(msg,hour)
    character(len=*),  intent(in)  :: msg
    real(dbl_kind),    intent(in)  :: hour
  end subroutine ufs_logfhour
#else

  use ice_kinds_mod , only : dbl_kind, int_kind, char_len, char_len_long

  implicit none

  real(dbl_kind) :: wtime = 0.0
contains
  subroutine ufs_settimer(timevalue)
    real(dbl_kind),    intent(inout) :: timevalue
    real(dbl_kind)                   :: MPI_Wtime
    timevalue = MPI_Wtime()
  end subroutine ufs_settimer

  subroutine ufs_logtimer(nunit,elapsedsecs,string,runtimelog,time0)
    integer,           intent(in)    :: nunit
    integer(int_kind), intent(in)    :: elapsedsecs
    character(len=*),  intent(in)    :: string
    logical,           intent(in)    :: runtimelog
    real(dbl_kind),    intent(in)    :: time0
    real(dbl_kind)                   :: MPI_Wtime, timevalue
    if (.not. runtimelog) return
    if (time0 > 0.) then
       timevalue = MPI_Wtime()-time0
       write(nunit,*)elapsedsecs,' CICE '//trim(string),timevalue
    end if
  end subroutine ufs_logtimer

  subroutine ufs_file_setLogUnit(filename,nunit,runtimelog)
    character(len=*),  intent(in)    :: filename
    logical,           intent(in)    :: runtimelog
    integer,           intent(out)   :: nunit
    if (.not. runtimelog) return
    open (newunit=nunit, file=trim(filename))
  end subroutine ufs_file_setLogUnit

  subroutine ufs_logfhour(msg,hour)
    character(len=*),  intent(in)    :: msg
    real(dbl_kind),    intent(in)    :: hour
    character(len=char_len)          :: filename
    integer(int_kind)                :: nunit
    write(filename,'(a,i3.3)')'log.ice.f',int(hour)
    open(newunit=nunit,file=trim(filename))
    write(nunit,'(a)')'completed: cice'
    write(nunit,'(a,f10.3)')'forecast hour:',hour
    write(nunit,'(a)')'valid time: '//trim(msg)
    close(nunit)
  end subroutine ufs_logfhour

  ! Define stub routines that do nothing - they are just here to avoid
  ! having cppdefs in the main program
  subroutine shr_file_setLogUnit(nunit)
    integer, intent(in) :: nunit
  end subroutine shr_file_setLogUnit
  subroutine shr_file_getLogUnit(nunit)
    integer, intent(in) :: nunit
  end subroutine shr_file_getLogUnit
  subroutine t_startf(string)
    character(len=*) :: string
  end subroutine t_startf
  subroutine t_stopf(string)
    character(len=*) :: string
  end subroutine t_stopf
  subroutine t_barrierf(string, comm)
    character(len=*) :: string
    integer:: comm
  end subroutine t_barrierf
#endif

end module cice_wrapper_mod
