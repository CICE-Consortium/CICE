!=======================================================

! module to help reading namelists

module ice_namelist_mod

  use ice_kinds_mod

  public :: goto_nml_group
  
!=======================================================
contains
!=======================================================

  subroutine goto_nml_group(iunit, nml, status)
    
    ! Search to namelist group within ice_in file.
    ! for compilers that do not allow optional namelists
    integer(kind=int_kind), intent(in) :: &
         iunit ! namelist file unit
    
    character(len=*), intent(in) :: &
         nml ! namelist to search for
    
    integer(kind=int_kind), intent(out) :: &
         status ! status of subrouine
 
    ! local variables
    character(len=char_len) :: &
         nml_str, & ! string in file
         test_str   ! string to test
    
    integer(kind=int_kind) :: &
         i, n ! dummy integers
    
    
    ! rewind file
    rewind(iunit)
    
    ! define test string with ampersand
    test_str = '&' // trim(adjustl(nml))
    
    ! search for the record containing the namelist group we're looking for
    do
       read(iunit, '(a)', iostat=status) nml_str
       if (status /= 0) then
          exit ! e.g. end of file
       else
          if (index(adjustl(nml_str), test_str) == 1) then
             exit ! i.e. found record we're looking for
          end if
       end if
    end do
    
    backspace(iunit)
    
  end subroutine goto_nml_group
  
!=======================================================================
end module ice_namelist_mod
      
!======================================================================
