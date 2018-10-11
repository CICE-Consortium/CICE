!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module ice_broadcast

!  This module contains all the broadcast routines.  This
!  particular version contains MPI versions of these routines.
!
! author: Phil Jones, LANL
! Oct. 2004: Adapted from POP version by William H. Lipscomb, LANL

   use ice_kinds_mod
   use ice_communicate, only: mpiR8, mpir4, MPI_COMM_ICE
   use ice_exit, only: abort_ice
   use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted

   implicit none
   private

   public  :: broadcast_scalar,         &
              broadcast_array

!-----------------------------------------------------------------------
!
!  generic interfaces for module procedures
!
!-----------------------------------------------------------------------

   interface broadcast_scalar
     module procedure broadcast_scalar_dbl,  &
                      broadcast_scalar_real, &
                      broadcast_scalar_int,  &
                      broadcast_scalar_log,  &
                      broadcast_scalar_char
   end interface

   interface broadcast_array
     module procedure broadcast_array_dbl_1d,  &
                      broadcast_array_real_1d, &
                      broadcast_array_int_1d,  &
                      broadcast_array_log_1d,  &
                      broadcast_array_dbl_2d,  &
                      broadcast_array_real_2d, &
                      broadcast_array_int_2d,  &
                      broadcast_array_log_2d,  &
                      broadcast_array_dbl_3d,  &
                      broadcast_array_real_3d, &
                      broadcast_array_int_3d,  &
                      broadcast_array_log_3d
   end interface

!***********************************************************************

 contains

!***********************************************************************

 subroutine broadcast_scalar_dbl(scalar, root_pe)

!  Broadcasts a scalar dbl variable from one processor (root_pe)
!  to all other processors. This is a specific instance of the generic
!  broadcast\_scalar interface.

   include 'mpif.h'  ! MPI Fortran include file

   integer (int_kind), intent(in) :: &
      root_pe              ! processor number to broadcast from

   real (dbl_kind), intent(inout) :: &
      scalar               ! scalar to be broadcast

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: ierr  ! local MPI error flag
   character(len=*), parameter :: subname = '(broadcast_scalar_dbl)'

!-----------------------------------------------------------------------

   call MPI_BCAST(scalar, 1, mpiR8, root_pe, MPI_COMM_ICE, ierr)
   call MPI_BARRIER(MPI_COMM_ICE, ierr)

!-----------------------------------------------------------------------

end subroutine broadcast_scalar_dbl

!***********************************************************************

subroutine broadcast_scalar_real(scalar, root_pe)

!  Broadcasts a scalar real variable from one processor (root_pe)
!  to all other processors. This is a specific instance of the generic
!  broadcast\_scalar interface.

   include 'mpif.h'  ! MPI Fortran include file

   integer (int_kind), intent(in) :: &
      root_pe              ! processor number to broadcast from

   real (real_kind), intent(inout) :: &
      scalar               ! scalar to be broadcast

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: ierr  ! local MPI error flag
   character(len=*), parameter :: subname = '(broadcast_scalar_real)'

!-----------------------------------------------------------------------

   call MPI_BCAST(scalar, 1, mpiR4, root_pe, MPI_COMM_ICE, ierr)
   call MPI_BARRIER(MPI_COMM_ICE, ierr)

!-----------------------------------------------------------------------

 end subroutine broadcast_scalar_real

!***********************************************************************

subroutine broadcast_scalar_int(scalar, root_pe)

!  Broadcasts a scalar integer variable from one processor (root_pe)
!  to all other processors. This is a specific instance of the generic
!  broadcast\_scalar interface.

   include 'mpif.h'  ! MPI Fortran include file

   integer (int_kind), intent(in) :: &
      root_pe              ! processor number to broadcast from

   integer (int_kind), intent(inout) :: &
      scalar                ! scalar to be broadcast

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: ierr  ! local MPI error flag
   character(len=*), parameter :: subname = '(broadcast_scalar_int)'

!-----------------------------------------------------------------------

   call MPI_BCAST(scalar, 1, MPI_INTEGER, root_pe, MPI_COMM_ICE,ierr)
   call MPI_BARRIER(MPI_COMM_ICE, ierr)

!-----------------------------------------------------------------------

 end subroutine broadcast_scalar_int

!***********************************************************************

subroutine broadcast_scalar_log(scalar, root_pe)

!  Broadcasts a scalar logical variable from one processor (root_pe)
!  to all other processors. This is a specific instance of the generic
!  broadcast\_scalar interface.

   include 'mpif.h'  ! MPI Fortran include file

   integer (int_kind), intent(in) :: &
     root_pe              ! processor number to broadcast from

   logical (log_kind), intent(inout) :: &
     scalar               ! scalar to be broadcast

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     itmp,               &! local temporary
     ierr                 ! MPI error flag
   character(len=*), parameter :: subname = '(broadcast_scalar_log)'

!-----------------------------------------------------------------------

   if (scalar) then
     itmp = 1
   else
     itmp = 0
   endif

   call MPI_BCAST(itmp, 1, MPI_INTEGER, root_pe, MPI_COMM_ICE, ierr)
   call MPI_BARRIER(MPI_COMM_ICE, ierr)

   if (itmp == 1) then
     scalar = .true.
   else
     scalar = .false.
   endif

!-----------------------------------------------------------------------

 end subroutine broadcast_scalar_log

!***********************************************************************

subroutine broadcast_scalar_char(scalar, root_pe)

!  Broadcasts a scalar character variable from one processor (root_pe)
!  to all other processors. This is a specific instance of the generic
!  broadcast\_scalar interface.

   include 'mpif.h'  ! MPI Fortran include file

   integer (int_kind), intent(in) :: &
     root_pe              ! processor number to broadcast from

   character (*), intent(inout) :: &
     scalar               ! scalar to be broadcast

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     clength,            &! length of character
     ierr                 ! MPI error flag
   character(len=*), parameter :: subname = '(broadcast_scalar_char)'

!-----------------------------------------------------------------------

   clength = len(scalar)

   call MPI_BCAST(scalar, clength, MPI_CHARACTER, root_pe, MPI_COMM_ICE, ierr)
   call MPI_BARRIER(MPI_COMM_ICE, ierr)

!--------------------------------------------------------------------

 end subroutine broadcast_scalar_char

!***********************************************************************

subroutine broadcast_array_dbl_1d(array, root_pe)

!  Broadcasts a vector dbl variable from one processor (root_pe)
!  to all other processors. This is a specific instance of the generic
!  broadcast\_array interface.

   include 'mpif.h'  ! MPI Fortran include file

   integer (int_kind), intent(in) :: &
     root_pe           ! processor number to broadcast from

   real (dbl_kind), dimension(:), intent(inout) :: &
     array             ! array to be broadcast

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     nelements,       &! size of array
     ierr              ! local MPI error flag
   character(len=*), parameter :: subname = '(broadcast_array_dbl_1d)'

!-----------------------------------------------------------------------

   nelements = size(array)

   call MPI_BCAST(array, nelements, mpiR8, root_pe, MPI_COMM_ICE, ierr)
   call MPI_BARRIER(MPI_COMM_ICE, ierr)

!-----------------------------------------------------------------------

 end subroutine broadcast_array_dbl_1d

!***********************************************************************

subroutine broadcast_array_real_1d(array, root_pe)

!  Broadcasts a real vector from one processor (root_pe)
!  to all other processors. This is a specific instance of the generic
!  broadcast\_array interface.

   include 'mpif.h'  ! MPI Fortran include file

   integer (int_kind), intent(in) :: &
     root_pe              ! processor number to broadcast from

   real (real_kind), dimension(:), intent(inout) :: &
     array                ! array to be broadcast

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag
   character(len=*), parameter :: subname = '(broadcast_array_real_1d)'

!-----------------------------------------------------------------------

   nelements = size(array)

   call MPI_BCAST(array, nelements, mpiR4, root_pe, MPI_COMM_ICE, ierr)
   call MPI_BARRIER(MPI_COMM_ICE, ierr)

!-----------------------------------------------------------------------

 end subroutine broadcast_array_real_1d

!***********************************************************************

subroutine broadcast_array_int_1d(array, root_pe)

!  Broadcasts an integer vector from one processor (root_pe)
!  to all other processors. This is a specific instance of the generic
!  broadcast\_array interface.

   include 'mpif.h'  ! MPI Fortran include file

   integer (int_kind), intent(in) :: &
     root_pe              ! processor number to broadcast from

   integer (int_kind), dimension(:), intent(inout) :: &
       array              ! array to be broadcast

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag
   character(len=*), parameter :: subname = '(broadcast_array_int_1d)'

!-----------------------------------------------------------------------

   nelements = size(array)

   call MPI_BCAST(array, nelements, MPI_INTEGER, root_pe, MPI_COMM_ICE, ierr)
   call MPI_BARRIER(MPI_COMM_ICE, ierr)

!-----------------------------------------------------------------------

 end subroutine broadcast_array_int_1d

!***********************************************************************

subroutine broadcast_array_log_1d(array, root_pe)

!  Broadcasts a logical vector from one processor (root_pe)
!  to all other processors. This is a specific instance of the generic
!  broadcast\_array interface.

   include 'mpif.h'  ! MPI Fortran include file

   integer (int_kind), intent(in) :: &
     root_pe              ! processor number to broadcast from

   logical (log_kind), dimension(:), intent(inout) :: &
     array                ! array to be broadcast

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind), dimension(:), allocatable :: &
      array_int            ! temporary array for MPI bcast

   integer (int_kind) :: &
      nelements,          &! size of array to be broadcast
      ierr                 ! local MPI error flag

   character(len=*), parameter :: subname = '(broadcast_array_log_1d)'

!-----------------------------------------------------------------------

   nelements = size(array)
   allocate(array_int(nelements))

   where (array)
     array_int = 1
   elsewhere
     array_int = 0
   end where

   call MPI_BCAST(array_int, nelements, MPI_INTEGER, root_pe, &
                  MPI_COMM_ICE, ierr)
   call MPI_BARRIER(MPI_COMM_ICE, ierr)

   where (array_int == 1)
     array = .true.
   elsewhere
     array = .false.
   end where

   deallocate(array_int)

!-----------------------------------------------------------------------

 end subroutine broadcast_array_log_1d

!***********************************************************************

 subroutine broadcast_array_dbl_2d(array, root_pe)

!  Broadcasts a dbl 2d array from one processor (root_pe)
!  to all other processors. This is a specific instance of the generic
!  broadcast\_array interface.

   include 'mpif.h'  ! MPI Fortran include file

   integer (int_kind), intent(in) :: &
     root_pe           ! processor number to broadcast from

   real (dbl_kind), dimension(:,:), intent(inout) :: &
     array             ! array to be broadcast

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      nelements,         &! size of array
      ierr                ! local MPI error flag
   character(len=*), parameter :: subname = '(broadcast_array_dbl_2d)'

!-----------------------------------------------------------------------

   nelements = size(array)

   call MPI_BCAST(array, nelements, mpiR8, root_pe, MPI_COMM_ICE, ierr)
   call MPI_BARRIER(MPI_COMM_ICE, ierr)

!-----------------------------------------------------------------------

 end subroutine broadcast_array_dbl_2d

!***********************************************************************

 subroutine broadcast_array_real_2d(array, root_pe)

!  Broadcasts a real 2d array from one processor (root_pe)
!  to all other processors. This is a specific instance of the generic
!  broadcast\_array interface.

   include 'mpif.h'  ! MPI Fortran include file

   integer (int_kind), intent(in) :: &
     root_pe              ! processor number to broadcast from

   real (real_kind), dimension(:,:), intent(inout) :: &
     array                ! array to be broadcast

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag
   character(len=*), parameter :: subname = '(broadcast_array_real_2d)'

!-----------------------------------------------------------------------

   nelements = size(array)

   call MPI_BCAST(array, nelements, mpiR4, root_pe, MPI_COMM_ICE, ierr)
   call MPI_BARRIER(MPI_COMM_ICE, ierr)

!-----------------------------------------------------------------------

 end subroutine broadcast_array_real_2d

!***********************************************************************

 subroutine broadcast_array_int_2d(array, root_pe)

!  Broadcasts a 2d integer array from one processor (root_pe)
!  to all other processors. This is a specific instance of the generic
!  broadcast\_array interface.

   include 'mpif.h'  ! MPI Fortran include file

   integer (int_kind), intent(in) :: &
     root_pe              ! processor number to broadcast from

   integer (int_kind), dimension(:,:), intent(inout) :: &
       array              ! array to be broadcast

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag
   character(len=*), parameter :: subname = '(broadcast_array_int_2d)'

!-----------------------------------------------------------------------

   nelements = size(array)

   call MPI_BCAST(array, nelements, MPI_INTEGER, root_pe, MPI_COMM_ICE, ierr)
   call MPI_BARRIER(MPI_COMM_ICE, ierr)

!-----------------------------------------------------------------------

 end subroutine broadcast_array_int_2d

!***********************************************************************

 subroutine broadcast_array_log_2d(array, root_pe)

!  Broadcasts a logical 2d array from one processor (root_pe)
!  to all other processors. This is a specific instance of the generic
!  broadcast\_array interface.

   include 'mpif.h'  ! MPI Fortran include file

   integer (int_kind), intent(in) :: &
     root_pe              ! processor number to broadcast from

   logical (log_kind), dimension(:,:), intent(inout) :: &
     array                ! array to be broadcast

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind), dimension(:,:), allocatable :: &
     array_int            ! temporary array for MPI bcast

   integer (int_kind) :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag

   character(len=*), parameter :: subname = '(broadcast_array_log_2d)'

!-----------------------------------------------------------------------

   nelements = size(array)
   allocate(array_int(size(array,dim=1),size(array,dim=2)))

   where (array)
     array_int = 1
   elsewhere
     array_int = 0
   end where

   call MPI_BCAST(array_int, nelements, MPI_INTEGER, root_pe, &
                  MPI_COMM_ICE, ierr)
   call MPI_BARRIER(MPI_COMM_ICE, ierr)

   where (array_int == 1)
     array = .true.
   elsewhere
     array = .false.
   end where

   deallocate(array_int)

!-----------------------------------------------------------------------

 end subroutine broadcast_array_log_2d

!***********************************************************************

 subroutine broadcast_array_dbl_3d(array, root_pe)

!  Broadcasts a double 3d array from one processor (root_pe)
!  to all other processors. This is a specific instance of the generic
!  broadcast\_array interface.

   include 'mpif.h'  ! MPI Fortran include file

   integer (int_kind), intent(in) :: &
     root_pe           ! processor number to broadcast from

   real (dbl_kind), dimension(:,:,:), intent(inout) :: &
     array             ! array to be broadcast

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     nelements,       &! size of array
     ierr              ! local MPI error flag
   character(len=*), parameter :: subname = '(broadcast_array_dbl_3d)'

!-----------------------------------------------------------------------

   nelements = size(array)

   call MPI_BCAST(array, nelements, mpiR8, root_pe, MPI_COMM_ICE, ierr)
   call MPI_BARRIER(MPI_COMM_ICE, ierr)

!-----------------------------------------------------------------------

 end subroutine broadcast_array_dbl_3d

!***********************************************************************

 subroutine broadcast_array_real_3d(array, root_pe)

!  Broadcasts a real 3d array from one processor (root_pe)
!  to all other processors. This is a specific instance of the generic
!  broadcast\_array interface.

   include 'mpif.h'  ! MPI Fortran include file

   integer (int_kind), intent(in) :: &
     root_pe              ! processor number to broadcast from

   real (real_kind), dimension(:,:,:), intent(inout) :: &
     array                ! array to be broadcast

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag
   character(len=*), parameter :: subname = '(broadcast_array_real_3d)'

!-----------------------------------------------------------------------

   nelements = size(array)

   call MPI_BCAST(array, nelements, mpiR4, root_pe, MPI_COMM_ICE, ierr)
   call MPI_BARRIER(MPI_COMM_ICE, ierr)

!-----------------------------------------------------------------------

 end subroutine broadcast_array_real_3d

!***********************************************************************

 subroutine broadcast_array_int_3d(array, root_pe)

!  Broadcasts an integer 3d array from one processor (root_pe)
!  to all other processors. This is a specific instance of the generic
!  broadcast\_array interface.

   include 'mpif.h'  ! MPI Fortran include file

   integer (int_kind), intent(in) :: &
     root_pe              ! processor number to broadcast from

   integer (int_kind), dimension(:,:,:), intent(inout) :: &
       array              ! array to be broadcast

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag
   character(len=*), parameter :: subname = '(broadcast_array_int_3d)'

!-----------------------------------------------------------------------

   nelements = size(array)

   call MPI_BCAST(array, nelements, MPI_INTEGER, root_pe, MPI_COMM_ICE, ierr)
   call MPI_BARRIER(MPI_COMM_ICE, ierr)

!-----------------------------------------------------------------------

 end subroutine broadcast_array_int_3d

!***********************************************************************

 subroutine broadcast_array_log_3d(array, root_pe)

!  Broadcasts a logical 3d array from one processor (root_pe)
!  to all other processors. This is a specific instance of the generic
!  broadcast\_array interface.

   include 'mpif.h'  ! MPI Fortran include file

   integer (int_kind), intent(in) :: &
     root_pe              ! processor number to broadcast from

   logical (log_kind), dimension(:,:,:), intent(inout) :: &
     array                ! array to be broadcast

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind), dimension(:,:,:), allocatable :: &
     array_int            ! temporary array for MPI bcast

   integer (int_kind) :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag

   character(len=*), parameter :: subname = '(broadcast_array_log_3d)'

!-----------------------------------------------------------------------

   nelements = size(array)
   allocate(array_int(size(array,dim=1), &
                      size(array,dim=2), &
                      size(array,dim=3)))

   where (array)
     array_int = 1
   elsewhere
     array_int = 0
   end where

   call MPI_BCAST(array_int, nelements, MPI_INTEGER, root_pe, &
                  MPI_COMM_ICE, ierr)
   call MPI_BARRIER(MPI_COMM_ICE, ierr)

   where (array_int == 1)
     array = .true.
   elsewhere
     array = .false.
   end where

   deallocate(array_int)

!-----------------------------------------------------------------------

 end subroutine broadcast_array_log_3d

!***********************************************************************

 end module ice_broadcast

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
