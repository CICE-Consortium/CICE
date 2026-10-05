! Stub `mpi` module, for gfortran -fsyntax-only sweeps ONLY.
!
! There is no MPI on this machine, and the flattened comm/ directory uses
! `use mpi` unconditionally, so without this the sweep stops at
! ice_communicate and everything downstream goes unchecked -- which is most of
! the code this refactor touches.
!
! Only names the CICE source actually references are declared; the values are
! arbitrary.  Nothing here is ever linked or run.
module mpi

   implicit none
   public
   save

   integer, parameter :: MPI_STATUS_SIZE = 8
   integer, parameter :: MPI_SUCCESS     = 0
   integer, parameter :: MPI_COMM_WORLD  = 1
   integer, parameter :: MPI_CHARACTER   = 10
   integer, parameter :: MPI_INTEGER     = 11
   integer, parameter :: MPI_INTEGER8    = 12
   integer, parameter :: MPI_REAL        = 13
   integer, parameter :: MPI_REAL4       = 14
   integer, parameter :: MPI_REAL8       = 15
   integer, parameter :: MPI_REAL16      = 16
   integer, parameter :: MPI_COMPLEX16   = 17
   integer, parameter :: MPI_DOUBLE      = 18
   integer, parameter :: MPI_SUM         = 20
   integer, parameter :: MPI_MIN         = 21
   integer, parameter :: MPI_MAX         = 22

contains

   ! Only the ones used as functions need a body; everything else is called as
   ! an external subroutine and needs no interface.
   double precision function MPI_Wtime()
      MPI_Wtime = 0.0d0
   end function MPI_Wtime

end module mpi
