module ice_shr_methods

  use ESMF         , only : ESMF_LOGERR_PASSTHRU, ESMF_LogFoundError, ESMF_LOGMSG_ERROR, ESMF_MAXSTR
  use ESMF         , only : ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_FAILURE

  implicit none
  private

  public  :: chkerr


!===============================================================================
contains
!===============================================================================


  logical function chkerr(rc, line, file)

    integer, intent(in) :: rc
    integer, intent(in) :: line
    character(len=*), intent(in) :: file

    integer :: lrc
    character(len=*), parameter :: subname='(chkerr)'

    chkerr = .false.
    lrc = rc
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=line, file=file)) then
       chkerr = .true.
    endif
  end function chkerr

end module ice_shr_methods
