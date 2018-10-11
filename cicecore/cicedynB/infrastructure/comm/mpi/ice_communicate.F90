!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module ice_communicate

!  This module contains the necessary routines and variables for
!  communicating between processors.
!
! author: Phil Jones, LANL
! Oct. 2004: Adapted from POP version by William H. Lipscomb, LANL

   use ice_kinds_mod
   use ice_exit, only: abort_ice
   use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted

#if defined key_oasis3 || key_oasis3mct
   use cpl_oasis3
#endif

#if defined key_oasis4
   use cpl_oasis4
#endif

#if defined key_iomput
   use lib_mpp, only:   mpi_comm_opa      ! MPP library
#endif

   implicit none
   private

   public  :: init_communicate,          &
              get_num_procs,             &
              ice_barrier,               &
              create_communicator

   integer (int_kind), public :: &
      MPI_COMM_ICE,             &! MPI communicator for ice comms
      mpiR16,                   &! MPI type for r16_kind
      mpiR8,                    &! MPI type for dbl_kind
      mpiR4,                    &! MPI type for real_kind
      my_task,                  &! MPI task number for this task
      master_task                ! task number of master task

   integer (int_kind), parameter, public :: &
      mpitagHalo            = 1,    &! MPI tags for various
      mpitag_gs             = 1000   ! communication patterns

!***********************************************************************

 contains

!***********************************************************************

 subroutine init_communicate(mpicom)

!  This routine sets up MPI environment and defines ice
!  communicator.

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   include 'mpif.h'   ! MPI Fortran include file

   integer (kind=int_kind), optional, intent(in) :: mpicom ! specified communicator

   integer (int_kind) :: ierr  ! MPI error flag
   logical            :: flag  ! MPI logical flag
   integer (int_kind) :: ice_comm

   character(len=*), parameter :: subname = '(init_communicate)'

!-----------------------------------------------------------------------
!
!  initiate mpi environment and create communicator for internal
!  ice communications
!
!-----------------------------------------------------------------------

   if (present(mpicom)) then
     ice_comm = mpicom
   else
#if (defined key_oasis3 || defined key_oasis3mct || defined key_oasis4)
     ice_comm = localComm       ! communicator from NEMO/OASISn 
#elif defined key_iomput
     ice_comm = mpi_comm_opa    ! communicator from NEMO/XIOS
#else
     ice_comm = MPI_COMM_WORLD  ! Global communicator 
#endif 
   endif

   call MPI_INITIALIZED(flag,ierr)
   if (.not.flag) call MPI_INIT(ierr)

   call MPI_BARRIER (ice_comm, ierr)
   call MPI_COMM_DUP(ice_comm, MPI_COMM_ICE, ierr)

   master_task = 0
   call MPI_COMM_RANK  (MPI_COMM_ICE, my_task, ierr)

   mpiR16 = MPI_REAL16
   mpiR8  = MPI_REAL8
   mpiR4  = MPI_REAL4

!-----------------------------------------------------------------------

 end subroutine init_communicate

!***********************************************************************

 function get_num_procs()

!  This function returns the number of processor assigned to
!  MPI_COMM_ICE

   integer (int_kind) :: get_num_procs

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: ierr
   character(len=*), parameter :: subname = '(get_num_procs)'

!-----------------------------------------------------------------------

   call MPI_COMM_SIZE(MPI_COMM_ICE, get_num_procs, ierr)

!-----------------------------------------------------------------------

 end function get_num_procs

!***********************************************************************

 subroutine ice_barrier()

!  This function calls an MPI_BARRIER

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: ierr
   character(len=*), parameter :: subname = '(ice_barrier)'

!-----------------------------------------------------------------------

   call MPI_BARRIER(MPI_COMM_ICE, ierr)

!-----------------------------------------------------------------------

 end subroutine ice_barrier

!***********************************************************************

 subroutine create_communicator(new_comm, num_procs)

!  This routine creates a separate communicator for a subset of
!  processors under default ice communicator.
!
!  this routine should be called from init_domain1 when the
!  domain configuration (e.g. nprocs_btrop) has been determined

   include 'mpif.h'

   integer (int_kind), intent(in) :: &
      num_procs         ! num of procs in new distribution

   integer (int_kind), intent(out) :: &
      new_comm          ! new communicator for this distribution

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     MPI_GROUP_ICE,         &! group of processors assigned to ice
     MPI_GROUP_NEW           ! group of processors assigned to new dist

   integer (int_kind) :: &
     ierr                    ! error flag for MPI comms

   integer (int_kind), dimension(3) :: &
     range                   ! range of tasks assigned to new dist
                             !  (assumed 0,num_procs-1)

   character(len=*), parameter :: subname = '(create_communicator)'

!-----------------------------------------------------------------------
!
!  determine group of processes assigned to distribution
!
!-----------------------------------------------------------------------

   call MPI_COMM_GROUP (MPI_COMM_ICE, MPI_GROUP_ICE, ierr)

   range(1) = 0
   range(2) = num_procs-1
   range(3) = 1

!-----------------------------------------------------------------------
!
!  create subroup and communicator for new distribution
!  note: MPI_COMM_CREATE must be called by all procs in MPI_COMM_ICE
!
!-----------------------------------------------------------------------

   call MPI_GROUP_RANGE_INCL(MPI_GROUP_ICE, 1, range, &
                             MPI_GROUP_NEW, ierr)

   call MPI_COMM_CREATE (MPI_COMM_ICE, MPI_GROUP_NEW,  &
                         new_comm, ierr)

!-----------------------------------------------------------------------

 end subroutine create_communicator

!***********************************************************************

 end module ice_communicate

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
