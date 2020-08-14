!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#define SERIAL_REMOVE_MPI

 module ice_global_reductions

!  This module contains all the routines for performing global
!  reductions like global sums, minvals, maxvals, etc.
!
! author: Phil Jones, LANL
! Oct. 2004: Adapted from POP version by William H. Lipscomb, LANL
! Feb. 2008: Updated from POP version by Elizabeth C. Hunke, LANL
! Aug. 2014: Added bit-for-bit reproducible options for global_sum_dbl
!            and global_sum_prod_dbl by T Craig NCAR
! Mar. 2019: Refactored bit-for-bit option, T Craig

#ifndef SERIAL_REMOVE_MPI
   use mpi   ! MPI Fortran module
#endif
   use ice_kinds_mod
   use ice_blocks, only: block, get_block, nx_block, ny_block
#ifdef SERIAL_REMOVE_MPI
   use ice_communicate, only: my_task, master_task
#else
   use ice_communicate, only: my_task, mpiR16, mpiR8, mpiR4, master_task
#endif
   use ice_constants, only: field_loc_Nface, field_loc_NEcorner, c0
   use ice_fileunits, only: bfbflag
   use ice_exit, only: abort_ice
   use ice_distribution, only: distrb, ice_distributionGet, &
       ice_distributionGetBlockID
   use ice_domain_size, only: nx_global, ny_global, max_blocks
   use ice_gather_scatter, only: gather_global
   use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
   use ice_reprosum, only: ice_reprosum_calc

   implicit none
   private

   public :: global_sum,      &
             global_allreduce_sum, &
             global_sum_prod, &
             global_maxval,   &
             global_minval

!-----------------------------------------------------------------------
!
!  generic interfaces for module procedures
!
!-----------------------------------------------------------------------

   interface global_sum
     module procedure global_sum_dbl,              &
                      global_sum_real,             &
                      global_sum_int,              &
                      global_sum_scalar_dbl,       &
                      global_sum_scalar_real,      &
                      global_sum_scalar_int
   end interface

   interface global_allreduce_sum
     module procedure global_allreduce_sum_vector_dbl!,   &
     ! module procedure global_allreduce_sum_vector_real, & ! not yet implemented
     ! module procedure global_allreduce_sum_vector_int     ! not yet implemented
   end interface

   interface global_sum_prod
     module procedure global_sum_prod_dbl,         &
                      global_sum_prod_real,        &
                      global_sum_prod_int
   end interface

   interface global_maxval
     module procedure global_maxval_dbl,           &
                      global_maxval_real,          &
                      global_maxval_int,           &
                      global_maxval_scalar_dbl,    &
                      global_maxval_scalar_real,   &
                      global_maxval_scalar_int
   end interface

   interface global_minval
     module procedure global_minval_dbl,           &
                      global_minval_real,          &
                      global_minval_int,           &
                      global_minval_scalar_dbl,    &
                      global_minval_scalar_real,   &
                      global_minval_scalar_int
   end interface

!***********************************************************************

 contains

!***********************************************************************

 function global_sum_dbl(array, dist, field_loc, mMask, lMask) &
          result(globalSum)

!  Computes the global sum of the physical domain of a 2-d array.
!
!  This is actually the specific interface for the generic global_sum
!  function corresponding to double precision arrays.  The generic
!  interface is identical but will handle real and integer 2-d slabs
!  and real, integer, and double precision scalars.

   real (dbl_kind), dimension(:,:,:), intent(in) :: &
      array                ! array to be summed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   integer (int_kind), intent(in) :: &
      field_loc            ! location of field on staggered grid

   real (dbl_kind), dimension(:,:,:), intent(in), optional :: &
      mMask                ! optional multiplicative mask

   logical (log_kind), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

   real (dbl_kind) :: &
      globalSum            ! resulting global sum

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,iblock,n, &! local counters
      ib,ie,jb,je,  &! beg,end of physical domain
      blockID,      &! block location
      numProcs,     &! number of processor participating
      numBlocks,    &! number of local blocks
      communicator, &! communicator for this distribution
      maxiglob       ! maximum non-redundant value of i_global

   logical (log_kind) :: &
      Nrow           ! this field is on a N row (a velocity row)

   real (dbl_kind), dimension(:,:), allocatable :: &
      work           ! temporary local array

   real (dbl_kind), dimension(:), allocatable :: &
      sums           ! array of sums

   type (block) :: &
      this_block     ! holds local block information

   character(len=*), parameter :: subname = '(global_sum_dbl)'

!-----------------------------------------------------------------------


   call ice_distributionGet(dist,          &
                            numLocalBlocks = numBlocks, &
                            nprocs = numProcs,       &
                            communicator = communicator)

   if (numBlocks == 0) then
      allocate(work(1,1))
   else
      allocate(work(nx_block*ny_block*numBlocks,1))
   endif
   allocate(sums(1))
   work = 0.0_dbl_kind
   sums = 0.0_dbl_kind

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

      !*** if this row along or beyond tripole boundary
      !*** must eliminate redundant points from global sum

      maxiglob = -1
      if (this_block%tripole) then
         Nrow=(field_loc == field_loc_Nface .or. &
               field_loc == field_loc_NEcorner)
         if (Nrow .and. this_block%tripoleTFlag) then
            maxiglob = 0 ! entire u-row on T-fold grid
         elseif (Nrow .or. this_block%tripoleTFlag) then
            maxiglob = nx_global/2 ! half T-row on T-fold and u-row on u-fold
         else
            maxiglob = -1 ! nothing to do for T-row on u-fold
         endif
      endif

      n = (iblock-1)*nx_block*ny_block

      do j=jb,je
      do i=ib,ie
         n = n + 1
         ! eliminate redundant points
         if (maxiglob > 0 .and. j == je .and. this_block%i_glob(i) > maxiglob) then
            work(n,1) = 0._dbl_kind
         else
            if (present(mMask)) then
               work(n,1) = array(i,j,iblock)*mMask(i,j,iblock)
            else if (present(lMask)) then
               if (lMask(i,j,iblock)) then
                  work(n,1) = array(i,j,iblock)
               endif
            else
               work(n,1) = array(i,j,iblock)
            endif
         endif
      end do
      end do
   end do

   call compute_sums_dbl(work,sums,communicator,numProcs)

   globalSum = sums(1)

   deallocate(work)
   deallocate(sums)

!-----------------------------------------------------------------------

 end function global_sum_dbl

!***********************************************************************

 function global_sum_real(array, dist, field_loc, mMask, lMask) &
          result(globalSum)

!  Computes the global sum of the physical domain of a 2-d array.
!
!  This is actually the specific interface for the generic global_sum
!  function corresponding to single precision arrays.  The generic
!  interface is identical but will handle double, real and integer 2-d slabs
!  and real, integer, and double precision scalars.

   real (real_kind), dimension(:,:,:), intent(in) :: &
      array                ! array to be summed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   integer (int_kind), intent(in) :: &
      field_loc            ! location of field on staggered grid

   real (real_kind), dimension(:,:,:), intent(in), optional :: &
      mMask                ! optional multiplicative mask

   logical (log_kind), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

   real (real_kind) :: &
      globalSum            ! resulting global sum

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,iblock,n, &! local counters
      ib,ie,jb,je,  &! beg,end of physical domain
      blockID,      &! block location
      numProcs,     &! number of processor participating
      numBlocks,    &! number of local blocks
      communicator, &! communicator for this distribution
      maxiglob       ! maximum non-redundant value of i_global

   logical (log_kind) :: &
      Nrow           ! this field is on a N row (a velocity row)

   real (dbl_kind), dimension(:,:), allocatable :: &
      work           ! temporary local array

   real (dbl_kind), dimension(:), allocatable :: &
      sums           ! array of sums

   type (block) :: &
      this_block     ! holds local block information

   character(len=*), parameter :: subname = '(global_sum_real)'

!-----------------------------------------------------------------------


   call ice_distributionGet(dist,          &
                            numLocalBlocks = numBlocks, &
                            nprocs = numProcs,       &
                            communicator = communicator)

   if (numBlocks == 0) then
      allocate(work(1,1))
   else
      allocate(work(nx_block*ny_block*numBlocks,1))
   endif
   allocate(sums(1))
   work = 0.0_dbl_kind
   sums = 0.0_dbl_kind

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

      !*** if this row along or beyond tripole boundary
      !*** must eliminate redundant points from global sum

      maxiglob = -1
      if (this_block%tripole) then
         Nrow=(field_loc == field_loc_Nface .or. &
               field_loc == field_loc_NEcorner)
         if (Nrow .and. this_block%tripoleTFlag) then
            maxiglob = 0 ! entire u-row on T-fold grid
         elseif (Nrow .or. this_block%tripoleTFlag) then
            maxiglob = nx_global/2 ! half T-row on T-fold and u-row on u-fold
         else
            maxiglob = -1 ! nothing to do for T-row on u-fold
         endif
      endif

      n = (iblock-1)*nx_block*ny_block

      do j=jb,je
      do i=ib,ie
         n = n + 1
         ! eliminate redundant points
         if (maxiglob > 0 .and. j == je .and. this_block%i_glob(i) > maxiglob) then
            work(n,1) = 0._dbl_kind
         else
            if (present(mMask)) then
               work(n,1) = real(array(i,j,iblock)*mMask(i,j,iblock),dbl_kind)
            else if (present(lMask)) then
               if (lMask(i,j,iblock)) then
                  work(n,1) = real(array(i,j,iblock),dbl_kind)
               endif
            else
               work(n,1) = real(array(i,j,iblock),dbl_kind)
            endif
         endif
      end do
      end do
   end do

   call compute_sums_dbl(work,sums,communicator,numProcs)

   globalSum = real(sums(1),real_kind)

   deallocate(work)
   deallocate(sums)

!-----------------------------------------------------------------------

 end function global_sum_real

!***********************************************************************

 function global_sum_int(array, dist, field_loc, mMask, lMask) &
          result(globalSum)

!  Computes the global sum of the physical domain of a 2-d array.
!
!  This is actually the specific interface for the generic global_sum
!  function corresponding to integer arrays.  The generic
!  interface is identical but will handle real and integer 2-d slabs
!  and real, integer, and double precision scalars.

   integer (int_kind), dimension(:,:,:), intent(in) :: &
      array                ! array to be summed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   integer (int_kind), intent(in) :: &
      field_loc            ! location of field on staggered grid

   integer (int_kind), dimension(:,:,:), intent(in), optional :: &
      mMask                ! optional multiplicative mask

   logical (log_kind), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

   integer (int_kind) :: &
      globalSum            ! resulting global sum

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      blockSum,     &! sum of local block domain
      localSum       ! sum of all local block domains

   integer (int_kind) :: &
      i,j,iblock,   &! local counters
      ib,ie,jb,je,  &! beg,end of physical domain
      ierr,         &! mpi error flag
      blockID,      &! block location
      numProcs,     &! number of processor participating
      numBlocks,    &! number of local blocks
      communicator, &! communicator for this distribution
      maxiglob       ! maximum non-redundant value of i_global

   logical (log_kind) :: &
      Nrow           ! this field is on a N row (a velocity row)

   type (block) :: &
      this_block     ! holds local block information

   character(len=*), parameter :: subname = '(global_sum_int)'

!-----------------------------------------------------------------------

   localSum  = 0_int_kind
   globalSum = 0_int_kind

   call ice_distributionGet(dist,          &
                            numLocalBlocks = numBlocks, &
                            nprocs = numProcs,       &
                            communicator = communicator)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

      !*** if this row along or beyond tripole boundary
      !*** must eliminate redundant points from global sum

      maxiglob = -1
      if (this_block%tripole) then
         Nrow=(field_loc == field_loc_Nface .or. &
               field_loc == field_loc_NEcorner)
         if (Nrow .and. this_block%tripoleTFlag) then
            maxiglob = 0 ! entire u-row on T-fold grid
         elseif (Nrow .or. this_block%tripoleTFlag) then
            maxiglob = nx_global/2 ! half T-row on T-fold and u-row on u-fold
         else
            maxiglob = -1 ! nothing to do for T-row on u-fold
         endif
      endif

      blockSum = 0_int_kind

      do j=jb,je
      do i=ib,ie
         ! eliminate redundant points
         if (maxiglob > 0 .and. j == je .and. this_block%i_glob(i) > maxiglob) then
!            blockSum = blockSum + 0_int_kind
         else
            if (present(mMask)) then
               blockSum = blockSum + array(i,j,iblock)*mMask(i,j,iblock)
            else if (present(lMask)) then
               if (lMask(i,j,iblock)) then
                  blockSum = blockSum + array(i,j,iblock)
              endif
            else
               blockSum = blockSum + array(i,j,iblock)
            endif
         endif
      end do
      end do

      !*** now add block sum to global sum

      localSum = localSum + blockSum

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local sum to global sum
!
!-----------------------------------------------------------------------

#ifdef SERIAL_REMOVE_MPI
   globalSum = localSum
#else
   if (my_task < numProcs) then
      call MPI_ALLREDUCE(localSum, globalSum, 1, &
                         MPI_INTEGER, MPI_SUM, communicator, ierr)
   endif
#endif

!-----------------------------------------------------------------------

 end function global_sum_int

!***********************************************************************

 function global_sum_scalar_dbl(scalar, dist) &
          result(globalSum)

!  Computes the global sum of a set of scalars distributed across
!  a parallel machine.
!
!  This is actually the specific interface for the generic global_sum
!  function corresponding to double precision scalars.  The generic
!  interface is identical but will handle real and integer 2-d slabs
!  and real, integer, and double precision scalars.

   real (dbl_kind), intent(in) :: &
      scalar               ! scalar to be summed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   real (dbl_kind) :: &
      globalSum            ! resulting global sum

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      ierr,         &! mpi error flag
      numProcs,     &! number of processor participating
      numBlocks,    &! number of local blocks
      communicator   ! communicator for this distribution

   real (dbl_kind), dimension(:,:), allocatable :: &
      work           ! temporary local array

   real (dbl_kind), dimension(:), allocatable :: &
      sums           ! array of sums

   character(len=*), parameter :: subname = '(global_sum_scalar_dbl)'

!-----------------------------------------------------------------------
!
!  get communicator for MPI calls
!
!-----------------------------------------------------------------------

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks, &
                            nprocs = numProcs,        &
                            communicator = communicator)


   allocate(work(1,1))
   allocate(sums(1))
   work(1,1) = scalar
   sums = 0.0_dbl_kind

   call compute_sums_dbl(work,sums,communicator,numProcs)

   globalSum = sums(1)

   deallocate(work)
   deallocate(sums)

!-----------------------------------------------------------------------

 end function global_sum_scalar_dbl

!***********************************************************************

 function global_sum_scalar_real(scalar, dist) &
          result(globalSum)

!  Computes the global sum of a set of scalars distributed across
!  a parallel machine.
!
!  This is actually the specific interface for the generic global_sum
!  function corresponding to real scalars.  The generic
!  interface is identical but will handle real and integer 2-d slabs
!  and real, integer, and double precision scalars.

   real (real_kind), intent(in) :: &
      scalar               ! scalar to be summed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   real (real_kind) :: &
      globalSum            ! resulting global sum

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      ierr,         &! mpi error flag
      numProcs,     &! number of processor participating
      numBlocks,    &! number of local blocks
      communicator   ! communicator for this distribution

   real (dbl_kind), dimension(:,:), allocatable :: &
      work           ! temporary local array

   real (dbl_kind), dimension(:), allocatable :: &
      sums           ! array of sums

   character(len=*), parameter :: subname = '(global_sum_scalar_real)'

!-----------------------------------------------------------------------
!
!  get communicator for MPI calls
!
!-----------------------------------------------------------------------

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks, &
                            nprocs = numProcs,        &
                            communicator = communicator)

   allocate(work(1,1))
   allocate(sums(1))
   work(1,1) = real(scalar,dbl_kind)
   sums = 0.0_dbl_kind

   call compute_sums_dbl(work,sums,communicator,numProcs)

   globalSum = real(sums(1),real_kind)

   deallocate(work)
   deallocate(sums)

!-----------------------------------------------------------------------

 end function global_sum_scalar_real

!***********************************************************************

 function global_sum_scalar_int(scalar, dist) &
          result(globalSum)

!  Computes the global sum of a set of scalars distributed across
!  a parallel machine.
!
!  This is actually the specific interface for the generic global_sum
!  function corresponding to integer scalars.  The generic
!  interface is identical but will handle real and integer 2-d slabs
!  and real, integer, and double precision scalars.

   integer (int_kind), intent(in) :: &
      scalar               ! scalar to be summed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   integer (int_kind) :: &
      globalSum            ! resulting global sum

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      ierr,         &! mpi error flag
      numProcs,     &! number of processor participating
      numBlocks,    &! number of local blocks
      communicator   ! communicator for this distribution

   character(len=*), parameter :: subname = '(global_sum_scalar_int)'

!-----------------------------------------------------------------------
!
!  get communicator for MPI calls
!
!-----------------------------------------------------------------------

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks, &
                            nprocs = numProcs,        &
                            communicator = communicator)

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local sum to global sum
!
!-----------------------------------------------------------------------

#ifdef SERIAL_REMOVE_MPI
   globalSum = scalar
#else
   if (my_task < numProcs) then
      call MPI_ALLREDUCE(scalar, globalSum, 1, &
                         MPI_INTEGER, MPI_SUM, communicator, ierr)
   endif
#endif

!-----------------------------------------------------------------------

 end function global_sum_scalar_int

!***********************************************************************

 function global_allreduce_sum_vector_dbl(vector, dist) &
          result(globalSums)

!  Computes the global sums of sets of scalars (elements of 'vector') 
!  distributed across a parallel machine.
!
!  This is actually the specific interface for the generic global_allreduce_sum
!  function corresponding to double precision vectors.  The generic
!  interface is identical but will handle real and integer vectors.

   real (dbl_kind), dimension(:), intent(in) :: &
      vector               ! vector whose components are to be summed

   type (distrb), intent(in) :: &
      dist                 ! block distribution

   real (dbl_kind), dimension(size(vector)) :: &
      globalSums           ! resulting array of global sums

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      ierr,         &! mpi error flag
      numProcs,     &! number of processor participating
      numBlocks,    &! number of local blocks
      communicator, &! communicator for this distribution
      numElem        ! number of elements in vector

   real (dbl_kind), dimension(:,:), allocatable :: &
      work           ! temporary local array

   character(len=*), parameter :: subname = '(global_allreduce_sum_vector_dbl)'

!-----------------------------------------------------------------------
!
!  get communicator for MPI calls
!
!-----------------------------------------------------------------------

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks, &
                            nprocs = numProcs,        &
                            communicator = communicator)

   numElem = size(vector)
   allocate(work(1,numElem))
   work(1,:) = vector
   globalSums = c0

   call compute_sums_dbl(work,globalSums,communicator,numProcs)

   deallocate(work)

!-----------------------------------------------------------------------

 end function global_allreduce_sum_vector_dbl

!***********************************************************************

 function global_sum_prod_dbl (array1, array2, dist, field_loc, &
                               mMask, lMask) &
          result(globalSum)

!  Computes the global sum of the physical domain of a product of
!  two 2-d arrays.
!
!  This is actually the specific interface for the generic 
!  global_sum_prod function corresponding to double precision arrays.
!  The generic interface is identical but will handle real and integer 
!  2-d slabs.

   real (dbl_kind), dimension(:,:,:), intent(in) :: &
      array1, array2       ! arrays whose product is to be summed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for arrays

   integer (int_kind), intent(in) :: &
      field_loc            ! location of field on staggered grid

   real (dbl_kind), dimension(:,:,:), intent(in), optional :: &
      mMask                ! optional multiplicative mask

   logical (log_kind), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

   real (dbl_kind) :: &
      globalSum            ! resulting global sum

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,iblock,n, &! local counters
      ib,ie,jb,je,  &! beg,end of physical domain
      blockID,      &! block location
      numProcs,     &! number of processor participating
      numBlocks,    &! number of local blocks
      communicator, &! communicator for this distribution
      maxiglob       ! maximum non-redundant value of i_global

   logical (log_kind) :: &
      Nrow           ! this field is on a N row (a velocity row)

   real (dbl_kind), dimension(:,:), allocatable :: &
      work           ! temporary local array

   real (dbl_kind), dimension(:), allocatable :: &
      sums           ! array of sums

   type (block) :: &
      this_block     ! holds local block information

   character(len=*), parameter :: subname = '(global_sum_prod_dbl)'

!-----------------------------------------------------------------------


   call ice_distributionGet(dist,          &
                            numLocalBlocks = numBlocks, &
                            nprocs = numProcs,       &
                            communicator = communicator)

   if (numBlocks == 0) then
      allocate(work(1,1))
   else
      allocate(work(nx_block*ny_block*numBlocks,1))
   endif
   allocate(sums(1))
   work = 0.0_dbl_kind
   sums = 0.0_dbl_kind

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

      !*** if this row along or beyond tripole boundary
      !*** must eliminate redundant points from global sum

      maxiglob = -1
      if (this_block%tripole) then
         Nrow=(field_loc == field_loc_Nface .or. &
               field_loc == field_loc_NEcorner)
         if (Nrow .and. this_block%tripoleTFlag) then
            maxiglob = 0 ! entire u-row on T-fold grid
         elseif (Nrow .or. this_block%tripoleTFlag) then
            maxiglob = nx_global/2 ! half T-row on T-fold and u-row on u-fold
         else
            maxiglob = -1 ! nothing to do for T-row on u-fold
         endif
      endif

      n = (iblock-1)*nx_block*ny_block

      do j=jb,je
      do i=ib,ie
         n = n + 1
         ! eliminate redundant points
         if (maxiglob > 0 .and. j == je .and. this_block%i_glob(i) > maxiglob) then
            work(n,1) = 0._dbl_kind
         else
            if (present(mMask)) then
               work(n,1) = array1(i,j,iblock)*array2(i,j,iblock)*mMask(i,j,iblock)
            else if (present(lMask)) then
               if (lMask(i,j,iblock)) then
                  work(n,1) = array1(i,j,iblock)*array2(i,j,iblock)
               endif
            else
               work(n,1) = array1(i,j,iblock)*array2(i,j,iblock)
            endif
         endif
      end do
      end do
   end do

   call compute_sums_dbl(work,sums,communicator,numProcs)

   globalSum = sums(1)

   deallocate(work)
   deallocate(sums)

!-----------------------------------------------------------------------

 end function global_sum_prod_dbl

!***********************************************************************

 function global_sum_prod_real (array1, array2, dist, field_loc, &
                                mMask, lMask) &
          result(globalSum)

!  Computes the global sum of the physical domain of a product of
!  two 2-d arrays.
!
!  This is actually the specific interface for the generic 
!  global_sum_prod function corresponding to single precision arrays.
!  The generic interface is identical but will handle real and integer 
!  2-d slabs.

   real (real_kind), dimension(:,:,:), intent(in) :: &
      array1, array2       ! arrays whose product is to be summed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for arrays

   integer (int_kind), intent(in) :: &
      field_loc            ! location of field on staggered grid

   real (real_kind), dimension(:,:,:), intent(in), optional :: &
      mMask                ! optional multiplicative mask

   logical (log_kind), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

   real (real_kind) :: &
      globalSum            ! resulting global sum

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,iblock,n, &! local counters
      ib,ie,jb,je,  &! beg,end of physical domain
      blockID,      &! block location
      numProcs,     &! number of processor participating
      numBlocks,    &! number of local blocks
      communicator, &! communicator for this distribution
      maxiglob       ! maximum non-redundant value of i_global

   logical (log_kind) :: &
      Nrow           ! this field is on a N row (a velocity row)

   real (dbl_kind), dimension(:,:), allocatable :: &
      work           ! temporary local array

   real (dbl_kind), dimension(:), allocatable :: &
      sums           ! array of sums

   type (block) :: &
      this_block     ! holds local block information

   character(len=*), parameter :: subname = '(global_sum_prod_dbl)'

!-----------------------------------------------------------------------


   call ice_distributionGet(dist,          &
                            numLocalBlocks = numBlocks, &
                            nprocs = numProcs,       &
                            communicator = communicator)

   if (numBlocks == 0) then
      allocate(work(1,1))
   else
      allocate(work(nx_block*ny_block*numBlocks,1))
   endif
   allocate(sums(1))
   work = 0.0_dbl_kind
   sums = 0.0_dbl_kind

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

      !*** if this row along or beyond tripole boundary
      !*** must eliminate redundant points from global sum

      maxiglob = -1
      if (this_block%tripole) then
         Nrow=(field_loc == field_loc_Nface .or. &
               field_loc == field_loc_NEcorner)
         if (Nrow .and. this_block%tripoleTFlag) then
            maxiglob = 0 ! entire u-row on T-fold grid
         elseif (Nrow .or. this_block%tripoleTFlag) then
            maxiglob = nx_global/2 ! half T-row on T-fold and u-row on u-fold
         else
            maxiglob = -1 ! nothing to do for T-row on u-fold
         endif
      endif

      n = (iblock-1)*nx_block*ny_block

      do j=jb,je
      do i=ib,ie
         n = n + 1
         ! eliminate redundant points
         if (maxiglob > 0 .and. j == je .and. this_block%i_glob(i) > maxiglob) then
            work(n,1) = 0._dbl_kind
         else
            if (present(mMask)) then
               work(n,1) = real(array1(i,j,iblock)*array2(i,j,iblock)*mMask(i,j,iblock),dbl_kind)
            else if (present(lMask)) then
               if (lMask(i,j,iblock)) then
                  work(n,1) = real(array1(i,j,iblock)*array2(i,j,iblock),dbl_kind)
               endif
            else
               work(n,1) = real(array1(i,j,iblock)*array2(i,j,iblock),dbl_kind)
            endif
         endif
      end do
      end do
   end do

   call compute_sums_dbl(work,sums,communicator,numProcs)

   globalSum = real(sums(1),real_kind)

   deallocate(work)
   deallocate(sums)

!-----------------------------------------------------------------------

 end function global_sum_prod_real

!***********************************************************************

 function global_sum_prod_int (array1, array2, dist, field_loc, &
                               mMask, lMask) &
          result(globalSum)

!  Computes the global sum of the physical domain of a product of
!  two 2-d arrays.
!
!  This is actually the specific interface for the generic 
!  global_sum_prod function corresponding to integer arrays.
!  The generic interface is identical but will handle real and integer 
!  2-d slabs.

   integer (int_kind), dimension(:,:,:), intent(in) :: &
      array1, array2       ! arrays whose product is to be summed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for arrays

   integer (int_kind), intent(in) :: &
      field_loc            ! location of field on staggered grid

   integer (int_kind), dimension(:,:,:), intent(in), optional :: &
      mMask                ! optional multiplicative mask

   logical (log_kind), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

   integer (int_kind) :: &
      globalSum            ! resulting global sum

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      blockSum,     &! sum of local block domain
      localSum       ! sum of all local block domains

   integer (int_kind) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      blockID,         &! block location
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator,    &! communicator for this distribution
      maxiglob          ! maximum non-redundant value of i_global

   logical (log_kind) :: &
      Nrow           ! this field is on a N row (a velocity row)

   type (block) :: &
      this_block          ! holds local block information

   character(len=*), parameter :: subname = '(global_sum_prod_int)'

!-----------------------------------------------------------------------

   localSum  = 0_int_kind
   globalSum = 0_int_kind

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks, &
                            nprocs = numProcs,        &
                            communicator = communicator)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

      !*** if this row along or beyond tripole boundary
      !*** must eliminate redundant points from global sum

      maxiglob = -1
      if (this_block%tripole) then
         Nrow=(field_loc == field_loc_Nface .or. &
               field_loc == field_loc_NEcorner)
         if (Nrow .and. this_block%tripoleTFlag) then
            maxiglob = 0 ! entire u-row on T-fold grid
         elseif (Nrow .or. this_block%tripoleTFlag) then
            maxiglob = nx_global/2 ! half T-row on T-fold and u-row on u-fold
         else
            maxiglob = -1 ! nothing to do for T-row on u-fold
         endif
      endif

      blockSum = 0_int_kind

      do j=jb,je
      do i=ib,ie
         ! eliminate redundant points
         if (maxiglob > 0 .and. j == je .and. this_block%i_glob(i) > maxiglob) then
!            blockSum = blockSum + 0_int_kind
         else
            if (present(mMask)) then
               blockSum = blockSum + array1(i,j,iblock)*array2(i,j,iblock)*mMask(i,j,iblock)
            else if (present(lMask)) then
               if (lMask(i,j,iblock)) then
                  blockSum = blockSum + array1(i,j,iblock)*array2(i,j,iblock)
               endif
            else
               blockSum = blockSum + array1(i,j,iblock)*array2(i,j,iblock)
            endif
         endif
      end do
      end do

      !*** now add block sum to global sum

      localSum = localSum + blockSum

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local sum to global sum
!
!-----------------------------------------------------------------------

#ifdef SERIAL_REMOVE_MPI
   globalSum = localSum
#else
   if (my_task < numProcs) then
      call MPI_ALLREDUCE(localSum, globalSum, 1, &
                         MPI_INTEGER, MPI_SUM, communicator, ierr)
   endif
#endif

!-----------------------------------------------------------------------

 end function global_sum_prod_int

!***********************************************************************

 function global_maxval_dbl (array, dist, lMask) &
          result(globalMaxval)

!  Computes the global maximum value of the physical domain of a 2-d field
!
!  This is actually the specific interface for the generic global_maxval
!  function corresponding to double precision arrays.  

   real (dbl_kind), dimension(:,:,:), intent(in) :: &
      array                ! array for which max value needed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   logical (log_kind), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

   real (dbl_kind) :: &
      globalMaxval         ! resulting maximum value of array

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (dbl_kind) ::    &
      blockMaxval,     &! sum of local block domain
      localMaxval       ! sum of all local block domains

   integer (int_kind) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator,    &! communicator for this distribution
      blockID           ! block location

   type (block) :: &
      this_block          ! holds local block information

   character(len=*), parameter :: subname = '(global_maxval_dbl)'

!-----------------------------------------------------------------------

   localMaxval  = -HUGE(0.0_dbl_kind)
   globalMaxval = -HUGE(0.0_dbl_kind)

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks, &
                            nprocs = numProcs,        &
                            communicator = communicator)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

      blockMaxval = -HUGE(0.0_dbl_kind)

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               blockMaxval = max(blockMaxval,array(i,j,iblock))
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            blockMaxval = max(blockMaxval,array(i,j,iblock))
         end do
         end do
      endif

      localMaxval = max(localMaxval,blockMaxval)

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local maxval to global maxval
!
!-----------------------------------------------------------------------

#ifdef SERIAL_REMOVE_MPI
   globalMaxval = localMaxval
#else
   if (my_task < numProcs) then
      call MPI_ALLREDUCE(localMaxval, globalMaxval, 1, &
                         mpiR8, MPI_MAX, communicator, ierr)
   endif
#endif

!-----------------------------------------------------------------------

 end function global_maxval_dbl

!***********************************************************************

 function global_maxval_real (array, dist, lMask) &
          result(globalMaxval)

!  Computes the global maximum value of the physical domain of a 2-d field
!
!  This is actually the specific interface for the generic global_maxval
!  function corresponding to single precision arrays.  

   real (real_kind), dimension(:,:,:), intent(in) :: &
      array                ! array for which max value needed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   logical (log_kind), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

   real (real_kind) :: &
      globalMaxval         ! resulting maximum value of array

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (real_kind) ::    &
      blockMaxval,     &! sum of local block domain
      localMaxval       ! sum of all local block domains

   integer (int_kind) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator,    &! communicator for this distribution
      blockID           ! block location

   type (block) :: &
      this_block          ! holds local block information

   character(len=*), parameter :: subname = '(global_maxval_real)'

!-----------------------------------------------------------------------

   localMaxval  = -HUGE(0.0_real_kind)
   globalMaxval = -HUGE(0.0_real_kind)

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks, &
                            nprocs = numProcs,        &
                            communicator = communicator)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

      blockMaxval = -HUGE(0.0_real_kind)

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               blockMaxval = max(blockMaxval,array(i,j,iblock))
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            blockMaxval = max(blockMaxval,array(i,j,iblock))
         end do
         end do
      endif

      localMaxval = max(localMaxval,blockMaxval)

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local maxval to global maxval
!
!-----------------------------------------------------------------------

#ifdef SERIAL_REMOVE_MPI
   globalMaxval = localMaxval
#else
   if (my_task < numProcs) then
      call MPI_ALLREDUCE(localMaxval, globalMaxval, 1, &
                         mpiR4, MPI_MAX, communicator, ierr)
   endif
#endif

!-----------------------------------------------------------------------

 end function global_maxval_real

!***********************************************************************

 function global_maxval_int (array, dist, lMask) &
          result(globalMaxval)

!  Computes the global maximum value of the physical domain of a 2-d field
!
!  This is actually the specific interface for the generic global_maxval
!  function corresponding to integer arrays.  

   integer (int_kind), dimension(:,:,:), intent(in) :: &
      array                ! array for which max value needed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   logical (log_kind), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

   integer (int_kind) :: &
      globalMaxval         ! resulting maximum value of array

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::    &
      blockMaxval,     &! sum of local block domain
      localMaxval       ! sum of all local block domains

   integer (int_kind) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator,    &! communicator for this distribution
      blockID           ! block location

   type (block) :: &
      this_block          ! holds local block information

   character(len=*), parameter :: subname = '(global_maxval_int)'

!-----------------------------------------------------------------------

   localMaxval  = -HUGE(0_int_kind)
   globalMaxval = -HUGE(0_int_kind)

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks, &
                            nprocs = numProcs,        &
                            communicator = communicator)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

      blockMaxval = -HUGE(0_int_kind)

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               blockMaxval = max(blockMaxval,array(i,j,iblock))
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            blockMaxval = max(blockMaxval,array(i,j,iblock))
         end do
         end do
      endif

      localMaxval = max(localMaxval,blockMaxval)

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local maxval to global maxval
!
!-----------------------------------------------------------------------

#ifdef SERIAL_REMOVE_MPI
   globalMaxval = localMaxval
#else
   if (my_task < numProcs) then
      call MPI_ALLREDUCE(localMaxval, globalMaxval, 1, &
                         MPI_INTEGER, MPI_MAX, communicator, ierr)
   endif
#endif

!-----------------------------------------------------------------------

 end function global_maxval_int

!***********************************************************************

 function global_maxval_scalar_dbl (scalar, dist) &
          result(globalMaxval)

!  Computes the global maximum value of a scalar value across
!  a distributed machine.
!
!  This is actually the specific interface for the generic global_maxval
!  function corresponding to double precision scalars.  

   real (dbl_kind), intent(in) :: &
      scalar               ! scalar for which max value needed

   type (distrb), intent(in) :: &
      dist                 ! block distribution

   real (dbl_kind) :: &
      globalMaxval         ! resulting maximum value

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      ierr,            &! mpi error flag
      numProcs,        &! number of processor participating
      communicator      ! communicator for this distribution

   character(len=*), parameter :: subname = '(global_maxval_scalar_dbl)'

!-----------------------------------------------------------------------

   call ice_distributionGet(dist, &
                            nprocs = numProcs,        &
                            communicator = communicator)

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local maxval to global maxval
!
!-----------------------------------------------------------------------

#ifdef SERIAL_REMOVE_MPI
   globalMaxval = scalar
#else
   if (my_task < numProcs) then
      call MPI_ALLREDUCE(scalar, globalMaxval, 1, &
                         mpiR8, MPI_MAX, communicator, ierr)
   endif
#endif

!-----------------------------------------------------------------------

 end function global_maxval_scalar_dbl

!***********************************************************************

 function global_maxval_scalar_real (scalar, dist) &
          result(globalMaxval)

!  Computes the global maximum value of a scalar value across
!  a distributed machine.
!
!  This is actually the specific interface for the generic global_maxval
!  function corresponding to single precision scalars.  

   real (real_kind), intent(in) :: &
      scalar               ! scalar for which max value needed

   type (distrb), intent(in) :: &
      dist                 ! block distribution

   real (real_kind) :: &
      globalMaxval         ! resulting maximum value

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      ierr,            &! mpi error flag
      numProcs,        &! number of processor participating
      communicator      ! communicator for this distribution

   character(len=*), parameter :: subname = '(global_maxval_scalar_real)'

!-----------------------------------------------------------------------

   call ice_distributionGet(dist, &
                            nprocs = numProcs,        &
                            communicator = communicator)

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local maxval to global maxval
!
!-----------------------------------------------------------------------

#ifdef SERIAL_REMOVE_MPI
   globalMaxval = scalar
#else
   if (my_task < numProcs) then
      call MPI_ALLREDUCE(scalar, globalMaxval, 1, &
                         mpiR4, MPI_MAX, communicator, ierr)
   endif
#endif

!-----------------------------------------------------------------------

 end function global_maxval_scalar_real

!***********************************************************************

 function global_maxval_scalar_int (scalar, dist) &
          result(globalMaxval)

!  Computes the global maximum value of a scalar value across
!  a distributed machine.
!
!  This is actually the specific interface for the generic global_maxval
!  function corresponding to single precision scalars.  

   integer (int_kind), intent(in) :: &
      scalar               ! scalar for which max value needed

   type (distrb), intent(in) :: &
      dist                 ! block distribution

   integer (int_kind) :: &
      globalMaxval         ! resulting maximum value

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      ierr,            &! mpi error flag
      numProcs,        &! number of processor participating
      communicator      ! communicator for this distribution

   character(len=*), parameter :: subname = '(global_maxval_scalar_int)'

!-----------------------------------------------------------------------

   call ice_distributionGet(dist, &
                            nprocs = numProcs,        &
                            communicator = communicator)

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local maxval to global maxval
!
!-----------------------------------------------------------------------

#ifdef SERIAL_REMOVE_MPI
   globalMaxval = scalar
#else
   if (my_task < numProcs) then
      call MPI_ALLREDUCE(scalar, globalMaxval, 1, &
                         MPI_INTEGER, MPI_MAX, communicator, ierr)
   endif
#endif

!-----------------------------------------------------------------------

 end function global_maxval_scalar_int

!***********************************************************************

 function global_minval_dbl (array, dist, lMask) &
          result(globalMinval)

!  Computes the global minimum value of the physical domain of a 2-d field
!
!  This is actually the specific interface for the generic global_minval
!  function corresponding to double precision arrays.  

   real (dbl_kind), dimension(:,:,:), intent(in) :: &
      array                ! array for which min value needed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   logical (log_kind), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

   real (dbl_kind) :: &
      globalMinval         ! resulting minimum value of array

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (dbl_kind) ::    &
      blockMinval,     &! sum of local block domain
      localMinval       ! sum of all local block domains

   integer (int_kind) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator,    &! communicator for this distribution
      blockID           ! block location

   type (block) :: &
      this_block          ! holds local block information

   character(len=*), parameter :: subname = '(global_minval_dbl)'

!-----------------------------------------------------------------------

   localMinval  = HUGE(0.0_dbl_kind)
   globalMinval = HUGE(0.0_dbl_kind)

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks, &
                            nprocs = numProcs,        &
                            communicator = communicator)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

      blockMinval = HUGE(0.0_dbl_kind)

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               blockMinval = min(blockMinval,array(i,j,iblock))
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            blockMinval = min(blockMinval,array(i,j,iblock))
         end do
         end do
      endif

      localMinval = min(localMinval,blockMinval)

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local minval to global minval
!
!-----------------------------------------------------------------------

#ifdef SERIAL_REMOVE_MPI
   globalMinval = localMinval
#else
   if (my_task < numProcs) then
      call MPI_ALLREDUCE(localMinval, globalMinval, 1, &
                         mpiR8, MPI_MIN, communicator, ierr)
   endif
#endif

!-----------------------------------------------------------------------

 end function global_minval_dbl

!***********************************************************************

 function global_minval_real (array, dist, lMask) &
          result(globalMinval)

!  Computes the global minimum value of the physical domain of a 2-d field
!
!  This is actually the specific interface for the generic global_minval
!  function corresponding to single precision arrays.  

   real (real_kind), dimension(:,:,:), intent(in) :: &
      array                ! array for which min value needed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   logical (log_kind), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

   real (real_kind) :: &
      globalMinval         ! resulting minimum value of array

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (real_kind) ::    &
      blockMinval,     &! sum of local block domain
      localMinval       ! sum of all local block domains

   integer (int_kind) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator,    &! communicator for this distribution
      blockID           ! block location

   type (block) :: &
      this_block          ! holds local block information

   character(len=*), parameter :: subname = '(global_minval_real)'

!-----------------------------------------------------------------------

   localMinval  = HUGE(0.0_real_kind)
   globalMinval = HUGE(0.0_real_kind)

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks, &
                            nprocs = numProcs,        &
                            communicator = communicator)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

      blockMinval = HUGE(0.0_real_kind)

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               blockMinval = min(blockMinval,array(i,j,iblock))
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            blockMinval = min(blockMinval,array(i,j,iblock))
         end do
         end do
      endif

      localMinval = min(localMinval,blockMinval)

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local minval to global minval
!
!-----------------------------------------------------------------------

#ifdef SERIAL_REMOVE_MPI
   globalMinval = localMinval
#else
   if (my_task < numProcs) then
      call MPI_ALLREDUCE(localMinval, globalMinval, 1, &
                         mpiR4, MPI_MIN, communicator, ierr)
   endif
#endif

!-----------------------------------------------------------------------

 end function global_minval_real

!***********************************************************************

 function global_minval_int (array, dist, lMask) &
          result(globalMinval)

!  Computes the global minimum value of the physical domain of a 2-d field
!
!  This is actually the specific interface for the generic global_minval
!  function corresponding to integer arrays.  

   integer (int_kind), dimension(:,:,:), intent(in) :: &
      array                ! array for which min value needed

   type (distrb), intent(in) :: &
      dist                 ! block distribution for array X

   logical (log_kind), dimension(:,:,:), intent(in), optional :: &
      lMask                ! optional logical mask

   integer (int_kind) :: &
      globalMinval         ! resulting minimum value of array

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::    &
      blockMinval,     &! sum of local block domain
      localMinval       ! sum of all local block domains

   integer (int_kind) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      ierr,            &! mpi error flag
      numBlocks,       &! number of local blocks
      numProcs,        &! number of processor participating
      communicator,    &! communicator for this distribution
      blockID           ! block location

   type (block) :: &
      this_block          ! holds local block information

   character(len=*), parameter :: subname = '(global_minval_int)'

!-----------------------------------------------------------------------

   localMinval  = HUGE(0_int_kind)
   globalMinval = HUGE(0_int_kind)

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks, &
                            nprocs = numProcs,        &
                            communicator = communicator)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

      blockMinval = HUGE(0_int_kind)

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               blockMinval = min(blockMinval,array(i,j,iblock))
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            blockMinval = min(blockMinval,array(i,j,iblock))
         end do
         end do
      endif

      localMinval = min(localMinval,blockMinval)

   end do

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local minval to global minval
!
!-----------------------------------------------------------------------

#ifdef SERIAL_REMOVE_MPI
   globalMinval = localMinval
#else
   if (my_task < numProcs) then
      call MPI_ALLREDUCE(localMinval, globalMinval, 1, &
                         MPI_INTEGER, MPI_MIN, communicator, ierr)
   endif
#endif

!-----------------------------------------------------------------------

 end function global_minval_int

!***********************************************************************

 function global_minval_scalar_dbl (scalar, dist) &
          result(globalMinval)

!  Computes the global minimum value of a scalar value across
!  a distributed machine.
!
!  This is actually the specific interface for the generic global_minval
!  function corresponding to double precision scalars.  

   real (dbl_kind), intent(in) :: &
      scalar               ! scalar for which min value needed

   type (distrb), intent(in) :: &
      dist                 ! block distribution

   real (dbl_kind) :: &
      globalMinval         ! resulting minimum value

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      ierr,            &! mpi error flag
      numProcs,        &! number of processor participating
      communicator      ! communicator for this distribution

   character(len=*), parameter :: subname = '(global_minval_scalar_dbl)'

!-----------------------------------------------------------------------

   call ice_distributionGet(dist, &
                            nprocs = numProcs,        &
                            communicator = communicator)

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local minval to global minval
!
!-----------------------------------------------------------------------

#ifdef SERIAL_REMOVE_MPI
   globalMinval = scalar
#else
   if (my_task < numProcs) then
      call MPI_ALLREDUCE(scalar, globalMinval, 1, &
                         mpiR8, MPI_MIN, communicator, ierr)
   endif
#endif

!-----------------------------------------------------------------------

 end function global_minval_scalar_dbl

!***********************************************************************

 function global_minval_scalar_real (scalar, dist) &
          result(globalMinval)

!  Computes the global minimum value of a scalar value across
!  a distributed machine.
!
!  This is actually the specific interface for the generic global_minval
!  function corresponding to single precision scalars.  

   real (real_kind), intent(in) :: &
      scalar               ! scalar for which min value needed

   type (distrb), intent(in) :: &
      dist                 ! block distribution

   real (real_kind) :: &
      globalMinval         ! resulting minimum value

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      ierr,            &! mpi error flag
      numProcs,        &! number of processor participating
      communicator      ! communicator for this distribution

   character(len=*), parameter :: subname = '(global_minval_scalar_real)'

!-----------------------------------------------------------------------

   call ice_distributionGet(dist, &
                            nprocs = numProcs,        &
                            communicator = communicator)

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local minval to global minval
!
!-----------------------------------------------------------------------

#ifdef SERIAL_REMOVE_MPI
   globalMinval = scalar
#else
   if (my_task < numProcs) then
      call MPI_ALLREDUCE(scalar, globalMinval, 1, &
                         mpiR4, MPI_MIN, communicator, ierr)
   endif
#endif

!-----------------------------------------------------------------------

 end function global_minval_scalar_real

!***********************************************************************

 function global_minval_scalar_int (scalar, dist) &
          result(globalMinval)

!  Computes the global minimum value of a scalar value across
!  a distributed machine.
!
!  This is actually the specific interface for the generic global_minval
!  function corresponding to single precision scalars.  

   integer (int_kind), intent(in) :: &
      scalar               ! scalar for which min value needed

   type (distrb), intent(in) :: &
      dist                 ! block distribution

   integer (int_kind) :: &
      globalMinval         ! resulting minimum value

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      ierr,            &! mpi error flag
      numProcs,        &! number of processor participating
      communicator      ! communicator for this distribution

   character(len=*), parameter :: subname = '(global_minval_scalar_int)'

!-----------------------------------------------------------------------

   call ice_distributionGet(dist, &
                            nprocs = numProcs,        &
                            communicator = communicator)

!-----------------------------------------------------------------------
!
!  now use MPI global reduction to reduce local minval to global minval
!
!-----------------------------------------------------------------------

#ifdef SERIAL_REMOVE_MPI
   globalMinval = scalar
#else
   if (my_task < numProcs) then
      call MPI_ALLREDUCE(scalar, globalMinval, 1, &
                         MPI_INTEGER, MPI_MIN, communicator, ierr)
   endif
#endif

!-----------------------------------------------------------------------

 end function global_minval_scalar_int

!***********************************************************************
!***********************************************************************

subroutine compute_sums_dbl(array2,sums8,mpicomm,numprocs)

! Computes the global sum of a 2-d array over fields
! with first dimension values and second dimension fields
!
! Several different options are supported.
! lsum4 = local sum with real*4 and scalar mpi allreduce, unlikely to be bfb
! lsum8 = local sum with real*8 and scalar mpi allreduce, unlikely to be bfb
! lsum16 = local sum with real*16 and scalar mpi allreduce, likely to be bfb
!    WARNING: this does not work in several compilers and mpi
!    implementations due to support for quad precision and consistency
!    between underlying datatype in fortran and c.  The source code
!    can be turned off with a cpp NO_R16.  Otherwise, it is recommended
!    that the results be validated on any platform where it might be used.
! reprosum = fixed point method based on ordered double integer sums.
!    that requires two scalar reductions per global sum.
!    This is extremely likely to be bfb.
!    (See Mirin and Worley, 2012, IJHPCA, 26, 1730, 
!    https://journals.sagepub.com/doi/10.1177/1094342011412630)
! ddpdd = parallel double-double algorithm using single scalar reduction.
!    This is very likely to be bfb.
!    (See He and Ding, 2001, Journal of Supercomputing, 18, 259,
!    https://link.springer.com/article/10.1023%2FA%3A1008153532043)

   real (dbl_kind), dimension(:,:), intent(in) :: &
      array2               ! array to be summed

   real (dbl_kind), dimension(:), intent(inout) :: &
      sums8                 ! resulting global sum

   integer(int_kind), intent(in) :: &
      mpicomm

   integer(int_kind), intent(in) :: &
      numprocs

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (real_kind), allocatable :: psums4(:)
   real (real_kind), allocatable :: sums4(:)
   real (dbl_kind) , allocatable :: psums8(:)
#ifndef NO_R16
   real (r16_kind) , allocatable :: psums16(:)
   real (r16_kind) , allocatable :: sums16(:)
#endif

   integer (int_kind) :: ns,nf,i,j, ierr

   character(len=*), parameter :: subname = '(compute_sums_dbl)'

!-----------------------------------------------------------------------

   sums8 = 0._dbl_kind
   ns = size(array2,dim=1)
   nf = size(array2,dim=2)

   if (bfbflag == 'off' .or. bfbflag == 'lsum8') then
      allocate(psums8(nf))
      psums8(:) = 0._dbl_kind

      do j = 1, nf
      do i = 1, ns
         psums8(j) = psums8(j) + array2(i,j)
      enddo
      enddo

#ifdef SERIAL_REMOVE_MPI
      sums8 = psums8
#else
      if (my_task < numProcs) then
         call MPI_ALLREDUCE(psums8, sums8, nf, mpiR8, MPI_SUM, mpicomm, ierr)
      endif
#endif

      deallocate(psums8)

#ifndef NO_R16
   elseif (bfbflag == 'lsum16') then
      allocate(psums16(nf))
      psums16(:) = 0._r16_kind
      allocate(sums16(nf))
      sums16(:) = 0._r16_kind

      do j = 1, nf
      do i = 1, ns
         psums16(j) = psums16(j) + real(array2(i,j),r16_kind)
      enddo
      enddo

#ifdef SERIAL_REMOVE_MPI
      sums16 = psums16
#else
      if (my_task < numProcs) then
         call MPI_ALLREDUCE(psums16, sums16, nf, mpiR16, MPI_SUM, mpicomm, ierr)
      endif
#endif
      sums8 = real(sums16,dbl_kind)

      deallocate(psums16,sums16)
#endif

   elseif (bfbflag == 'lsum4') then
      allocate(psums4(nf))
      psums4(:) = 0._real_kind
      allocate(sums4(nf))
      sums4(:) = 0._real_kind

      do j = 1, nf
      do i = 1, ns
         psums4(j) = psums4(j) + real(array2(i,j),real_kind)
      enddo
      enddo

#ifdef SERIAL_REMOVE_MPI
      sums4 = psums4
#else
      if (my_task < numProcs) then
         call MPI_ALLREDUCE(psums4, sums4, nf, mpiR4, MPI_SUM, mpicomm, ierr)
      endif
#endif
      sums8 = real(sums4,dbl_kind)

      deallocate(psums4,sums4)

   elseif (bfbflag == 'ddpdd') then
      if (my_task < numProcs) then
         call ice_reprosum_calc(array2,sums8,ns,ns,nf,ddpdd_sum=.true.,commid=mpicomm)
      endif

   elseif (bfbflag == 'reprosum') then
      if (my_task < numProcs) then
         call ice_reprosum_calc(array2,sums8,ns,ns,nf,ddpdd_sum=.false.,commid=mpicomm)
      endif

   else
      call abort_ice(subname//'ERROR: bfbflag unknown '//trim(bfbflag))
   endif

end subroutine compute_sums_dbl

!***********************************************************************

end module ice_global_reductions

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
