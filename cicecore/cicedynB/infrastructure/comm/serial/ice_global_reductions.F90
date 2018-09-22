!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module ice_global_reductions

!  This module contains all the routines for performing global
!  reductions like global sums, minvals, maxvals, etc.
!
! author: Phil Jones, LANL
! Oct. 2004: Adapted from POP version by William H. Lipscomb, LANL
! Feb. 2008: Updated from POP version by Elizabeth C. Hunke, LANL

   use ice_kinds_mod
   use ice_constants, only: field_loc_Nface, field_loc_NEcorner
   use ice_blocks, only: block, get_block
   use ice_distribution, only: distrb, ice_distributionGet, &
       ice_distributionGetBlockID
   use ice_domain_size, only: nx_global
   use ice_exit, only: abort_ice
   use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted

   implicit none
   private

   public :: global_sum,      &
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
      dist                 ! block distribution for array 

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

#ifdef REPRODUCIBLE
   real (r16_kind) :: &
      localSum,     &! sum of local block domain
      globalSumTmp   ! higher precision global sum
#else
   real (dbl_kind) :: &
      localSum       ! sum of local block domain
#endif

   integer (int_kind) :: &
      i,j,iblock,   &! local counters
      ib,ie,jb,je,  &! beg,end of physical domain
      blockID,      &! block location
      numBlocks,    &! number of local blocks
      maxiglob       ! maximum non-redundant value of i_global

   logical (log_kind) :: &
      Nrow           ! this field is on a N row (a velocity row)

   type (block) :: &
      this_block     ! holds local block information

   character(len=*), parameter :: subname = '(global_sum_dbl)'

!-----------------------------------------------------------------------

#ifdef REPRODUCIBLE
   globalSumTmp  = 0.0_r16_kind
#else
   globalSum = 0.0_dbl_kind
#endif

   call ice_distributionGet(dist,          &
                            numLocalBlocks = numBlocks)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

#ifdef REPRODUCIBLE
      localSum = 0.0_r16_kind
#else
      localSum = 0.0_dbl_kind
#endif

      if (present(mMask)) then
         do j=jb,je
         do i=ib,ie
            localSum = &
            localSum + array(i,j,iblock)*mMask(i,j,iblock)
         end do
         end do
      else if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               localSum = &
               localSum + array(i,j,iblock)
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            localSum = localSum + array(i,j,iblock)
         end do
         end do
      endif

      !*** if this row along or beyond tripole boundary
      !*** must eliminate redundant points from global sum

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

         if (maxiglob > 0) then

            j = je

            if (present(mMask)) then
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     localSum = &
                     localSum - array(i,j,iblock)*mMask(i,j,iblock)
                  endif
               end do
            else if (present(lMask)) then
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     if (lMask(i,j,iblock)) &
                     localSum = localSum - array(i,j,iblock)
                  endif
               end do
            else
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     localSum = localSum - array(i,j,iblock)
                  endif
               end do
            endif

         endif
      endif

      !*** now add block sum to global sum

#ifdef REPRODUCIBLE
      globalSumTmp = globalSumTmp + localSum
#else
      globalSum = globalSum + localSum
#endif

   end do

#ifdef REPRODUCIBLE
   globalSum = globalSumTmp
#endif

!-----------------------------------------------------------------------

 end function global_sum_dbl

!***********************************************************************

 function global_sum_real(array, dist, field_loc, mMask, lMask) &
          result(globalSum)

!  Computes the global sum of the physical domain of a 2-d array.
!
!  This is actually the specific interface for the generic global_sum
!  function corresponding to real arrays.  The generic
!  interface is identical but will handle real and integer 2-d slabs
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

#ifdef REPRODUCIBLE
   real (dbl_kind) :: &
      localSum,     &! sum of local block domain
      globalSumTmp   ! higher precision global sum
#else
   real (real_kind) :: &
      localSum       ! sum of local block domain
#endif

   integer (int_kind) :: &
      i,j,iblock,   &! local counters
      ib,ie,jb,je,  &! beg,end of physical domain
      blockID,      &! block location
      numBlocks,    &! number of local blocks
      maxiglob       ! maximum non-redundant value of i_global

   logical (log_kind) :: &
      Nrow           ! this field is on a N row (a velocity row)

   type (block) :: &
      this_block     ! holds local block information

   character(len=*), parameter :: subname = '(global_sum_real)'

!-----------------------------------------------------------------------

#ifdef REPRODUCIBLE
   globalSumTmp  = 0.0_dbl_kind
#else
   globalSum  = 0.0_real_kind
#endif

   call ice_distributionGet(dist,          &
                            numLocalBlocks = numBlocks)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

#ifdef REPRODUCIBLE
      localSum = 0.0_dbl_kind
#else
      localSum = 0.0_real_kind
#endif

      if (present(mMask)) then
         do j=jb,je
         do i=ib,ie
            localSum = &
            localSum + array(i,j,iblock)*mMask(i,j,iblock)
         end do
         end do
      else if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               localSum = &
               localSum + array(i,j,iblock)
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            localSum = localSum + array(i,j,iblock)
         end do
         end do
      endif

      !*** if this row along or beyond tripole boundary
      !*** must eliminate redundant points from global sum

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

         if (maxiglob > 0) then

            j = je

            if (present(mMask)) then
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     localSum = &
                     localSum - array(i,j,iblock)*mMask(i,j,iblock)
                  endif
               end do
            else if (present(lMask)) then
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     if (lMask(i,j,iblock)) &
                     localSum = localSum - array(i,j,iblock)
                  endif
               end do
            else
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     localSum = localSum - array(i,j,iblock)
                  endif
               end do
            endif

         endif
      endif

      !*** now add block sum to global sum

#ifdef REPRODUCIBLE
      globalSumTmp = globalSumTmp + localSum
#else
      globalSum = globalSum + localSum
#endif

   end do

#ifdef REPRODUCIBLE
   globalSum = globalSumTmp
#endif

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
      localSum       ! sum of all local block domains

   integer (int_kind) :: &
      i,j,iblock,   &! local counters
      ib,ie,jb,je,  &! beg,end of physical domain
      blockID,      &! block location
      numBlocks,    &! number of local blocks
      maxiglob       ! maximum non-redundant value of i_global

   logical (log_kind) :: &
      Nrow           ! this field is on a N row (a velocity row)

   type (block) :: &
      this_block     ! holds local block information

   character(len=*), parameter :: subname = '(global_sum_int)'

!-----------------------------------------------------------------------

   globalSum = 0_int_kind

   call ice_distributionGet(dist,          &
                            numLocalBlocks = numBlocks)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

      localSum = 0

      if (present(mMask)) then
         do j=jb,je
         do i=ib,ie
            localSum = &
            localSum + array(i,j,iblock)*mMask(i,j,iblock)
         end do
         end do
      else if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               localSum = &
               localSum + array(i,j,iblock)
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            localSum = localSum + array(i,j,iblock)
         end do
         end do
      endif

      !*** if this row along or beyond tripole boundary
      !*** must eliminate redundant points from global sum

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

         if (maxiglob > 0) then

            j = je

            if (present(mMask)) then
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     localSum = &
                     localSum - array(i,j,iblock)*mMask(i,j,iblock)
                  endif
               end do
            else if (present(lMask)) then
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     if (lMask(i,j,iblock)) &
                     localSum = localSum - array(i,j,iblock)
                  endif
               end do
            else
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     localSum = localSum - array(i,j,iblock)
                  endif
               end do
            endif

         endif
      endif

      !*** now add block sum to global sum

      globalSum = globalSum + localSum

   end do

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
      dist                 ! block distribution 

   real (dbl_kind) :: &
      globalSum            ! resulting global sum

   character(len=*), parameter :: subname = '(global_sum_scalar_dbl)'

!-----------------------------------------------------------------------
!
!  no operation needed for serial execution
!
!-----------------------------------------------------------------------

   globalSum = scalar

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
      dist                 ! block distribution 

   real (real_kind) :: &
      globalSum            ! resulting global sum

   character(len=*), parameter :: subname = '(global_sum_scalar_real)'

!-----------------------------------------------------------------------
!
!  no operation needed for serial execution
!
!-----------------------------------------------------------------------

   globalSum = scalar

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
      dist                 ! block distribution 

   integer (int_kind) :: &
      globalSum            ! resulting global sum

   character(len=*), parameter :: subname = '(global_sum_scalar_int)'

!-----------------------------------------------------------------------
!
!  no operation needed for serial execution
!
!-----------------------------------------------------------------------

   globalSum = scalar

!-----------------------------------------------------------------------

 end function global_sum_scalar_int

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

#ifdef REPRODUCIBLE
   real (r16_kind) :: &
      localSum,      &! sum of local block domain
      globalSumTmp    ! higher precision for reproducibility
#else
   real (dbl_kind) :: &
      localSum        ! sum of local block domain
#endif

   integer (int_kind) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      blockID,         &! block location
      numBlocks,       &! number of local blocks
      maxiglob          ! maximum non-redundant value of i_global

   logical (log_kind) :: &
      Nrow           ! this field is on a N row (a velocity row)

   type (block) :: &
      this_block          ! holds local block information

   character(len=*), parameter :: subname = '(global_sum_prod_dbl)'

!-----------------------------------------------------------------------

#ifdef REPRODUCIBLE
   globalSum  = 0.0_r16_kind
#else
   globalSum  = 0.0_dbl_kind
#endif

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

#ifdef REPRODUCIBLE
      localSum = 0.0_r16_kind
#else
      localSum = 0.0_dbl_kind
#endif

      if (present(mMask)) then
         do j=jb,je
         do i=ib,ie
            localSum = &
            localSum + array1(i,j,iblock)*array2(i,j,iblock)* &
                       mMask(i,j,iblock)
         end do
         end do
      else if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               localSum = &
               localSum + array1(i,j,iblock)*array2(i,j,iblock)
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            localSum = localSum + array1(i,j,iblock)*array2(i,j,iblock)
         end do
         end do
      endif

      !*** if this row along or beyond tripole boundary
      !*** must eliminate redundant points from global sum

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

         if (maxiglob > 0) then

            j = je

            if (present(mMask)) then
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     localSum = &
                     localSum - array1(i,j,iblock)*array2(i,j,iblock)* &
                                mMask(i,j,iblock)
                  endif
               end do
            else if (present(lMask)) then
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     if (lMask(i,j,iblock)) &
                        localSum = localSum - &
                                   array1(i,j,iblock)*array2(i,j,iblock)
                  endif
               end do
            else
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     localSum = localSum - &
                                array1(i,j,iblock)*array2(i,j,iblock)
                  endif
               end do
            endif

         endif
      endif

      !*** now add block sum to global sum

#ifdef REPRODUCIBLE
      globalSumTmp = globalSumTmp + localSum
#else
      globalSum = globalSum + localSum
#endif

   end do

#ifdef REPRODUCIBLE
   globalSum = globalSumTmp
#endif

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

#ifdef REPRODUCIBLE
   real (dbl_kind) :: &
      localSum,      &! sum of local block domain
      globalSumTmp    ! higher precision for reproducibility
#else
   real (real_kind) :: &
      localSum        ! sum of local block domain
#endif

   integer (int_kind) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      blockID,         &! block location
      numBlocks,       &! number of local blocks
      maxiglob          ! maximum non-redundant value of i_global

   logical (log_kind) :: &
      Nrow           ! this field is on a N row (a velocity row)

   type (block) :: &
      this_block          ! holds local block information

   character(len=*), parameter :: subname = '(global_sum_prod_real)'

!-----------------------------------------------------------------------

#ifdef REPRODUCIBLE
   globalSumTmp = 0.0_dbl_kind
#else
   globalSum = 0.0_real_kind
#endif

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

#ifdef REPRODUCIBLE
      localSum = 0.0_dbl_kind
#else
      localSum = 0.0_real_kind
#endif

      if (present(mMask)) then
         do j=jb,je
         do i=ib,ie
            localSum = &
            localSum + array1(i,j,iblock)*array2(i,j,iblock)* &
                       mMask(i,j,iblock)
         end do
         end do
      else if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               localSum = &
               localSum + array1(i,j,iblock)*array2(i,j,iblock)
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            localSum = localSum + array1(i,j,iblock)*array2(i,j,iblock)
         end do
         end do
      endif

      !*** if this row along or beyond tripole boundary
      !*** must eliminate redundant points from global sum

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

         if (maxiglob > 0) then

            j = je

            if (present(mMask)) then
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     localSum = &
                     localSum - array1(i,j,iblock)*array2(i,j,iblock)* &
                                mMask(i,j,iblock)
                  endif
               end do
            else if (present(lMask)) then
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     if (lMask(i,j,iblock)) &
                        localSum = localSum - &
                                   array1(i,j,iblock)*array2(i,j,iblock)
                  endif
               end do
            else
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     localSum = localSum - &
                                array1(i,j,iblock)*array2(i,j,iblock)
                  endif
               end do
            endif

         endif
      endif

      !*** now add block sum to global sum

#ifdef REPRODUCIBLE
      globalSumTmp = globalSumTmp + localSum
#else
      globalSum = globalSum + localSum
#endif

   end do

#ifdef REPRODUCIBLE
   globalSum = globalSumTmp
#endif

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
      localSum          ! sum of local block domain

   integer (int_kind) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      blockID,         &! block location
      numBlocks,       &! number of local blocks
      maxiglob          ! maximum non-redundant value of i_global

   logical (log_kind) :: &
      Nrow           ! this field is on a N row (a velocity row)

   type (block) :: &
      this_block          ! holds local block information

   character(len=*), parameter :: subname = '(global_sum_prod_int)'

!-----------------------------------------------------------------------

   globalSum = 0_int_kind

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

      localSum = 0

      if (present(mMask)) then
         do j=jb,je
         do i=ib,ie
            localSum = &
            localSum + array1(i,j,iblock)*array2(i,j,iblock)* &
                       mMask(i,j,iblock)
         end do
         end do
      else if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               localSum = &
               localSum + array1(i,j,iblock)*array2(i,j,iblock)
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            localSum = localSum + array1(i,j,iblock)*array2(i,j,iblock)
         end do
         end do
      endif

      !*** if this row along or beyond tripole boundary
      !*** must eliminate redundant points from global sum

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

         if (maxiglob > 0) then

            j = je

            if (present(mMask)) then
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     localSum = &
                     localSum - array1(i,j,iblock)*array2(i,j,iblock)* &
                                mMask(i,j,iblock)
                  endif
               end do
            else if (present(lMask)) then
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     if (lMask(i,j,iblock)) &
                        localSum = localSum - &
                                   array1(i,j,iblock)*array2(i,j,iblock)
                  endif
               end do
            else
               do i=ib,ie
                  if (this_block%i_glob(i) > maxiglob) then
                     localSum = localSum - &
                                array1(i,j,iblock)*array2(i,j,iblock)
                  endif
               end do
            endif

         endif
      endif

      !*** now add block sum to global sum

      globalSum = globalSum + localSum

   end do

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
      localMaxval       ! sum of local block domain

   integer (int_kind) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      numBlocks,       &! number of local blocks
      blockID           ! block location

   type (block) :: &
      this_block          ! holds local block information

   character(len=*), parameter :: subname = '(global_maxval_dbl)'

!-----------------------------------------------------------------------

   globalMaxval = -HUGE(0.0_dbl_kind)

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

      localMaxval = -HUGE(0.0_dbl_kind)

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               localMaxval = max(localMaxval,array(i,j,iblock))
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            localMaxval = max(localMaxval,array(i,j,iblock))
         end do
         end do
      endif

      globalMaxval = max(globalMaxval,localMaxval)

   end do

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
      localMaxval       ! sum of local block domain

   integer (int_kind) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      numBlocks,       &! number of local blocks
      blockID           ! block location

   type (block) :: &
      this_block          ! holds local block information

   character(len=*), parameter :: subname = '(global_maxval_real)'

!-----------------------------------------------------------------------

   globalMaxval = -HUGE(0.0_real_kind)

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

      localMaxval = -HUGE(0.0_real_kind)

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               localMaxval = max(localMaxval,array(i,j,iblock))
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            localMaxval = max(localMaxval,array(i,j,iblock))
         end do
         end do
      endif

      globalMaxval = max(globalMaxval,localMaxval)

   end do

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
      localMaxval       ! sum of local block domain

   integer (int_kind) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      numBlocks,       &! number of local blocks
      blockID           ! block location

   type (block) :: &
      this_block          ! holds local block information

   character(len=*), parameter :: subname = '(global_maxval_int)'

!-----------------------------------------------------------------------

   globalMaxval = -HUGE(0_int_kind)

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

      localMaxval = -HUGE(0_int_kind)

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               localMaxval = max(localMaxval,array(i,j,iblock))
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            localMaxval = max(localMaxval,array(i,j,iblock))
         end do
         end do
      endif

      globalMaxval = max(globalMaxval,localMaxval)

   end do

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

   character(len=*), parameter :: subname = '(global_maxval_scalar_dbl)'

!-----------------------------------------------------------------------
!
!  no operations required for serial execution
!
!-----------------------------------------------------------------------

   globalMaxval = scalar

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

   character(len=*), parameter :: subname = '(global_maxval_scalar_real)'

!-----------------------------------------------------------------------
!
!  no operations required for serial execution
!
!-----------------------------------------------------------------------

   globalMaxval = scalar

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

   character(len=*), parameter :: subname = '(global_maxval_scalar_int)'

!-----------------------------------------------------------------------
!
!  no operations required for serial execution
!
!-----------------------------------------------------------------------

   globalMaxval = scalar

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
      localMinval       ! sum of local block domain

   integer (int_kind) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      numBlocks,       &! number of local blocks
      blockID           ! block location

   type (block) :: &
      this_block          ! holds local block information

   character(len=*), parameter :: subname = '(global_minval_dbl)'

!-----------------------------------------------------------------------

   globalMinval = HUGE(0.0_dbl_kind)

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

      localMinval = HUGE(0.0_dbl_kind)

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               localMinval = min(localMinval,array(i,j,iblock))
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            localMinval = min(localMinval,array(i,j,iblock))
         end do
         end do
      endif

      globalMinval = min(globalMinval,localMinval)

   end do

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
      localMinval       ! sum of local block domain

   integer (int_kind) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      numBlocks,       &! number of local blocks
      blockID           ! block location

   type (block) :: &
      this_block          ! holds local block information

   character(len=*), parameter :: subname = '(global_minval_real)'

!-----------------------------------------------------------------------

   globalMinval = HUGE(0.0_real_kind)

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

      localMinval = HUGE(0.0_real_kind)

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               localMinval = min(localMinval,array(i,j,iblock))
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            localMinval = min(localMinval,array(i,j,iblock))
         end do
         end do
      endif

      globalMinval = min(globalMinval,localMinval)

   end do

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
      localMinval       ! sum of local block domain

   integer (int_kind) :: &
      i,j,iblock,      &! local counters
      ib,ie,jb,je,     &! beg,end of physical domain
      numBlocks,       &! number of local blocks
      blockID           ! block location

   type (block) :: &
      this_block          ! holds local block information

   character(len=*), parameter :: subname = '(global_minval_int)'

!-----------------------------------------------------------------------

   globalMinval = HUGE(0_int_kind)

   call ice_distributionGet(dist, &
                            numLocalBlocks = numBlocks)

   do iblock=1,numBlocks
      call ice_distributionGetBlockID(dist, iblock, blockID)

      this_block = get_block(blockID, blockID)

      ib = this_block%ilo
      ie = this_block%ihi
      jb = this_block%jlo
      je = this_block%jhi

      localMinval = HUGE(0_int_kind)

      if (present(lMask)) then
         do j=jb,je
         do i=ib,ie
            if (lMask(i,j,iblock)) then
               localMinval = min(localMinval,array(i,j,iblock))
            endif
         end do
         end do
      else
         do j=jb,je
         do i=ib,ie
            localMinval = min(localMinval,array(i,j,iblock))
         end do
         end do
      endif

      globalMinval = min(globalMinval,localMinval)

   end do

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

   character(len=*), parameter :: subname = '(global_minval_scalar_dbl)'

!-----------------------------------------------------------------------
!
!  no operations required for serial execution
!
!-----------------------------------------------------------------------

   globalMinval = scalar

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

   character(len=*), parameter :: subname = '(global_minval_scalar_real)'

!-----------------------------------------------------------------------
!
!  no operations required for serial execution
!
!-----------------------------------------------------------------------

   globalMinval = scalar

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

   character(len=*), parameter :: subname = '(global_minval_scalar_int)'

!-----------------------------------------------------------------------
!
!  no operations required for serial execution
!
!-----------------------------------------------------------------------

   globalMinval = scalar

!-----------------------------------------------------------------------

 end function global_minval_scalar_int

!***********************************************************************

 end module ice_global_reductions

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
