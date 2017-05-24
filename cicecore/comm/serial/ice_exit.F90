!  SVN:$Id: ice_exit.F90 700 2013-08-15 19:17:39Z eclare $
!=======================================================================
!
! Exit the model.
!
! authors William H. Lipscomb (LANL)
!         Elizabeth C. Hunke (LANL)
! 2006 ECH: separated serial and mpi functionality

      module ice_exit

      implicit none
      public

!=======================================================================

      contains

!=======================================================================

      subroutine abort_ice(error_message)

!  This routine aborts the ice model and prints an error message.

      use ice_fileunits, only: nu_diag, flush_fileunit
#ifdef CCSM
      use shr_sys_mod
#endif

      character (len=*), intent(in) :: error_message

#ifdef CCSM
      call shr_sys_abort(error_message)
#else
      write (nu_diag,*) error_message
      call flush_fileunit(nu_diag)
      stop
#endif

      end subroutine abort_ice

!=======================================================================

      subroutine end_run

! Ends parallel run by calling MPI_FINALIZE.
! Does nothing in serial runs.

      end subroutine end_run

!=======================================================================

      end module ice_exit

!=======================================================================
