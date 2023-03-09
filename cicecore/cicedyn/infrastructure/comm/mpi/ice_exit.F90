!=======================================================================
!
! Exit the model. 
! authors William H. Lipscomb (LANL)
!         Elizabeth C. Hunke (LANL)
! 2006 ECH: separated serial and mpi functionality

      module ice_exit

      use ice_kinds_mod
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted

      implicit none
      public

!=======================================================================

      contains

!=======================================================================

      subroutine abort_ice(error_message, file, line, doabort)

!  This routine aborts the ice model and prints an error message.

#if (defined CESMCOUPLED)
      use ice_fileunits, only: nu_diag, flush_fileunit
      use shr_sys_mod
#else
      use ice_fileunits, only: nu_diag, ice_stderr, flush_fileunit
      use mpi   ! MPI Fortran module
#endif

      character (len=*), intent(in),optional :: error_message  ! error message
      character (len=*), intent(in),optional :: file           ! file
      integer (kind=int_kind), intent(in), optional :: line    ! line number
      logical (kind=log_kind), intent(in), optional :: doabort ! abort flag

      ! local variables

#ifndef CESMCOUPLED
      integer (int_kind) :: &
         ierr,       & ! MPI error flag
         error_code    ! return code
#endif
      logical (log_kind) :: ldoabort   ! local doabort flag
      character(len=*), parameter :: subname='(abort_ice)'

      ldoabort = .true.
      if (present(doabort)) ldoabort = doabort

#if (defined CESMCOUPLED)
      call flush_fileunit(nu_diag)
      call icepack_warnings_flush(nu_diag)
      write(nu_diag,*) ' '
      write(nu_diag,*) subname, 'ABORTED: '
      if (present(file))   write (nu_diag,*) subname,' called from ',trim(file)
      if (present(line))   write (nu_diag,*) subname,' line number ',line
      if (present(error_message)) write (nu_diag,*) subname,' error = ',trim(error_message)
      call flush_fileunit(nu_diag)
      if (ldoabort) call shr_sys_abort(subname//trim(error_message))
#else
      call flush_fileunit(nu_diag)
      call icepack_warnings_flush(nu_diag)
      write(ice_stderr,*) ' '
      write(ice_stderr,*) subname, 'ABORTED: '
      if (present(file))   write (ice_stderr,*) subname,' called from ',trim(file)
      if (present(line))   write (ice_stderr,*) subname,' line number ',line
      if (present(error_message)) write (ice_stderr,*) subname,' error = ',trim(error_message)
      call flush_fileunit(ice_stderr)
      error_code = 128
      if (ldoabort) then
         call MPI_ABORT(MPI_COMM_WORLD, error_code, ierr)
         stop
      endif
#endif

      end subroutine abort_ice

!=======================================================================

      subroutine end_run

! Ends run by calling MPI_FINALIZE.

      integer (int_kind) :: ierr ! MPI error flag
      character(len=*), parameter :: subname = '(end_run)'

      call MPI_FINALIZE(ierr)

      end subroutine end_run

!=======================================================================

      end module ice_exit

!=======================================================================
