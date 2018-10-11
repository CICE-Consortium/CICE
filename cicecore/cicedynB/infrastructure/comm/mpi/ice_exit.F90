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

      subroutine abort_ice(error_message, file, line)

!  This routine aborts the ice model and prints an error message.

#if (defined CESMCOUPLED)
      use ice_fileunits, only: nu_diag, flush_fileunit
      use shr_sys_mod
#else
      use ice_fileunits, only: nu_diag, ice_stderr, flush_fileunit
      include 'mpif.h'   ! MPI Fortran include file
#endif

      character (len=*), intent(in),optional :: error_message
      character (len=*), intent(in),optional :: file
      integer (kind=int_kind), intent(in), optional :: &
         line       ! line number

      ! local variables

#ifndef CESMCOUPLED
      integer (int_kind) :: ierr ! MPI error flag
#endif
      character(len=*), parameter :: subname='(abort_ice)'

#if (defined CESMCOUPLED)
      call flush_fileunit(nu_diag)
      call icepack_warnings_flush(nu_diag)
      write(nu_diag,*) ' '
      write(nu_diag,*) subname, 'ABORTED: '
      if (present(file))   write (nu_diag,*) subname,' called from ',trim(file)
      if (present(line))   write (nu_diag,*) subname,' line number ',line
      if (present(error_message)) write (nu_diag,*) subname,' error = ',trim(error_message)
      call flush_fileunit(nu_diag)
      call shr_sys_abort(subname//trim(error_message))
#else
      call flush_fileunit(nu_diag)
      call icepack_warnings_flush(nu_diag)
      write(ice_stderr,*) ' '
      write(ice_stderr,*) subname, 'ABORTED: '
      if (present(file))   write (ice_stderr,*) subname,' called from ',trim(file)
      if (present(line))   write (ice_stderr,*) subname,' line number ',line
      if (present(error_message)) write (ice_stderr,*) subname,' error = ',trim(error_message)
      call flush_fileunit(ice_stderr)
      call MPI_ABORT(MPI_COMM_WORLD, ierr)
      stop
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
