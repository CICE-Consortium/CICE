!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module ice_communicate

!  This module contains the necessary routines and variables for
!  communicating between processors.  This instance of the module
!  is for serial execution so not much is done.
!
! author: Phil Jones, LANL
! Oct. 2004: Adapted from POP version by William H. Lipscomb, LANL

   use ice_kinds_mod
   use ice_exit, only: abort_ice
   use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted

   implicit none
   private

   public  :: init_communicate,          &
              get_num_procs,             &
              ice_barrier,               &
              create_communicator

   integer (int_kind), public :: &
      MPI_COMM_ICE,             &! MPI communicator for ice comms
      mpi_dbl,                  &! MPI type for dbl_kind
      my_task,                  &! MPI task number for this task
      master_task                ! task number of master task

   logical (log_kind), public :: &
      add_mpi_barriers      = .false. ! turn on mpi barriers for throttling

!***********************************************************************

 contains

!***********************************************************************

 subroutine init_communicate

!  This routine sets up MPI environment and defines ice communicator.

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   character(len=*), parameter :: subname = '(init_communicate)'

!-----------------------------------------------------------------------
!
!  initiate mpi environment and create communicator for internal
!  ice communications
!
!-----------------------------------------------------------------------

   my_task = 0
   master_task = 0

!-----------------------------------------------------------------------

 end subroutine init_communicate

!***********************************************************************

 function get_num_procs()

!  This function returns the number of processors assigned to
!  the ice model.

   integer (int_kind) :: get_num_procs

   character(len=*), parameter :: subname = '(get_num_procs)'

!-----------------------------------------------------------------------
!
!  serial execution, must be only 1
!
!-----------------------------------------------------------------------

   get_num_procs = 1

!-----------------------------------------------------------------------

 end function get_num_procs

!***********************************************************************

 subroutine ice_barrier()

!  This function is an MPI_BARRIER on the MPI side

   character(len=*), parameter :: subname = '(ice_barrier)'

!-----------------------------------------------------------------------
!
!  serial execution, no-op
!
!-----------------------------------------------------------------------

   ! do nothing

!-----------------------------------------------------------------------

 end subroutine ice_barrier

!***********************************************************************

 subroutine create_communicator(new_comm, num_procs)

!  This routine creates a separate communicator for a subset of
!  processors under default ice communicator.
!
!  this routine should be called from init_domain1 when the
!  domain configuration (e.g. nprocs_btrop) has been determined

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      num_procs         ! num of procs in new distribution

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) :: &
      new_comm          ! new communicator for this distribution

   new_comm = MPI_COMM_ICE

!-----------------------------------------------------------------------

 end subroutine create_communicator

!***********************************************************************

 end module ice_communicate

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
