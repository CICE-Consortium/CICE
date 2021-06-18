!=======================================================================
!
! Exit the model.
!
! authors William H. Lipscomb (LANL)
!         Elizabeth C. Hunke (LANL)
! 2006 ECH: separated serial and mpi functionality

      module ice_exit

      use ice_kinds_mod
      use ice_fileunits, only: nu_diag, flush_fileunit
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
#ifdef CESMCOUPLED
      use shr_sys_mod
#endif

      implicit none
      public

!=======================================================================

      contains

!=======================================================================

      subroutine abort_ice(error_message,file,line,doabort)

!  This routine aborts the ice model and prints an error message.

      character (len=*), intent(in),optional :: error_message  ! error message
      character (len=*), intent(in),optional :: file           ! file
      integer (kind=int_kind), intent(in), optional :: line    ! line number
      logical (kind=log_kind), intent(in), optional :: doabort ! abort flag

      logical (kind=log_kind) :: ldoabort   ! local doabort
      character(len=*), parameter :: subname='(abort_ice)'

      ldoabort = .true.
      if (present(doabort)) ldoabort = doabort

#ifdef CESMCOUPLED
      call icepack_warnings_flush(nu_diag)
      write(nu_diag,*) ' '
      write(nu_diag,*) subname, 'ABORTED: '
      if (present(file))   write (nu_diag,*) subname,' called from ',trim(file)
      if (present(line))   write (nu_diag,*) subname,' line number ',line
      if (present(error_message)) write (nu_diag,*) subname,' error = ',trim(error_message)
      if (ldoabort) call shr_sys_abort(subname//trim(error_message))
#else
      call icepack_warnings_flush(nu_diag)
      write(nu_diag,*) ' '
      write(nu_diag,*) subname, 'ABORTED: '
      if (present(file))   write (nu_diag,*) subname,' called from ',trim(file)
      if (present(line))   write (nu_diag,*) subname,' line number ',line
      if (present(error_message)) write (nu_diag,*) subname,' error = ',trim(error_message)
      call flush_fileunit(nu_diag)
      if (ldoabort) stop
#endif

      end subroutine abort_ice

!=======================================================================

      subroutine end_run

      character(len=*), parameter :: subname = '(end_run)'

! Ends parallel run by calling MPI_FINALIZE.
! Does nothing in serial runs.

      end subroutine end_run

!=======================================================================

      end module ice_exit

!=======================================================================
