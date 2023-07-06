module cice_wrapper_mod

#ifdef CESMCOUPLED
  use perf_mod     , only : t_startf, t_stopf, t_barrierf
  use shr_file_mod , only : shr_file_getlogunit, shr_file_setlogunit

#else
  use ice_kinds_mod , only : dbl_kind, int_kind, char_len, char_len_long

  implicit none

  real(dbl_kind) :: timeiads = 0.0, timeirls = 0.0, timeadv = 0.0, timefs = 0.0
#endif
contains
  ! Define stub routines that do nothing - they are just here to avoid
  ! having cppdefs in the main program
#ifndef CESMCOUPLED
  subroutine ufs_settimer(timevalue)
    real(dbl_kind), intent(inout) :: timevalue
    real(dbl_kind)                :: MPI_Wtime
    timevalue = MPI_Wtime()
  end subroutine ufs_settimer
  subroutine ufs_logtimer(nunit,etime,string,time0)
    integer,          intent(in)  :: nunit
    integer(int_kind), intent(in) :: etime
    character(len=*), intent(in)  :: string
    real(dbl_kind),   intent(in)  :: time0
    real(dbl_kind)                :: MPI_Wtime, timevalue
    timevalue = MPI_Wtime()-time0
    write(nunit,*)etime,' CICE '//trim(string),timevalue
  end subroutine ufs_logtimer

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
#else
  subroutine ufs_settimer(timevalue)
    real(dbl_kind), intent(out) :: timevalue
  end subroutine ufs_settimer
  subroutine ufs_logtimer(nunit,etime,string,time0)
    character(len=*), intent(in)  :: string
    integer(int_kind), intent(in) :: etime
    real(dbl_kind),   intent(in)  :: time0
    integer,          intent(in)  :: nunit
  end subroutine ufs_logtimer
#endif

end module cice_wrapper_mod
