
!=======================================================================
!
! Exit the model.
! authors William H. Lipscomb (LANL)
!         Elizabeth C. Hunke (LANL)
! 2006 ECH: separated serial and mpi functionality

      module ice_exit

      use ice_kinds_mod
      use ice_fileunits, only: nu_diag, ice_stderr, flush_fileunit
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
#if (defined CESMCOUPLED)
      use shr_sys_mod
#else
#ifndef SERIAL_REMOVE_MPI
      use mpi   ! MPI Fortran module
#endif
#endif

      implicit none
      public

!=======================================================================

      contains

!=======================================================================

      subroutine abort_ice(error_message, file, line, doabort)

!  This routine aborts the ice model and prints an error message.

      character (len=*), intent(in),optional :: error_message  ! error message
      character (len=*), intent(in),optional :: file           ! file
      integer (kind=int_kind), intent(in), optional :: line    ! line number
      logical (kind=log_kind), intent(in), optional :: doabort ! abort flag

      ! local variables

      integer (int_kind) :: &
         ierr,       & ! MPI error flag
         outunit,    & ! output unit
         error_code    ! return code
      logical (log_kind) :: ldoabort   ! local doabort flag
      character(len=*), parameter :: subname='(abort_ice)'

      ldoabort = .true.
      if (present(doabort)) ldoabort = doabort

#if (defined CESMCOUPLED)
      outunit = nu_diag
#else
      outunit = ice_stderr
#endif

      call flush_fileunit(nu_diag)
      call icepack_warnings_flush(nu_diag)
      write(outunit,*) ' '
      write(outunit,*) subname, 'ABORTED: '
      if (present(file))   write (outunit,*) subname,' called from ',trim(file)
      if (present(line))   write (outunit,*) subname,' line number ',line
      if (present(error_message)) write (outunit,*) subname,' error = ',trim(error_message)
      call flush_fileunit(outunit)

      if (ldoabort) then
#if (defined CESMCOUPLED)
         call shr_sys_abort(subname//trim(error_message))
#else
#ifndef SERIAL_REMOVE_MPI
         error_code = 128
         call MPI_ABORT(MPI_COMM_WORLD, error_code, ierr)
#endif
         stop
#endif
      endif

      end subroutine abort_ice

!=======================================================================

      subroutine end_run

! Ends run by calling MPI_FINALIZE
! Does nothing in serial runs

      integer (int_kind) :: ierr ! MPI error flag
      character(len=*), parameter :: subname = '(end_run)'

#ifndef SERIAL_REMOVE_MPI
      call MPI_FINALIZE(ierr)
#endif

      end subroutine end_run

!=======================================================================

      end module ice_exit

!=======================================================================
