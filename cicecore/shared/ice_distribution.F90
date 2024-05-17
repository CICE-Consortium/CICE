!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module ice_distribution

!  This module provides data types and routines for distributing
!  blocks across processors.
!
! author: Phil Jones, LANL
! Oct. 2004: Adapted from POP by William H. Lipscomb, LANL
! Jan. 2008: Elizabeth Hunke updated to new POP infrastructure

   use ice_kinds_mod
   use ice_domain_size, only: max_blocks
   use ice_communicate, only: my_task, master_task, create_communicator
   use ice_blocks, only: nblocks_x, nblocks_y, nblocks_tot, debug_blocks
   use ice_exit, only: abort_ice
   use ice_fileunits, only: nu_diag
   use ice_memusage, only: ice_memusage_allocErr

   implicit none
   private
   save

   type, public :: distrb  ! distribution data type
      integer (int_kind) :: &
         nprocs            ,&! number of processors in this dist
         communicator      ,&! communicator to use in this dist
         numLocalBlocks      ! number of blocks distributed to this
                             !   local processor

      integer (int_kind), dimension(:), pointer :: &
         blockLocation     ,&! processor location for all blocks
         blockLocalID      ,&! local  block id for all blocks
         blockGlobalID       ! global block id for each local block

      integer (int_kind), dimension(:), pointer ::  blockCnt
   end type

   public :: create_distribution, &
             ice_distributionGet,         &
             ice_distributionGetBlockLoc, &
             ice_distributionGetBlockID, &
             create_local_block_ids, &
             proc_decomposition 

   character (char_len), public :: &
       processor_shape       ! 'square-pop' (approx) POP default config
                             ! 'square-ice' like square-pop but better for ice
                             ! 'slenderX1' (NPX x 1)
                             ! 'slenderX2' (NPX x 2)

!***********************************************************************

 contains

!***********************************************************************

 function create_distribution(dist_type, nprocs, work_per_block)

!  This routine determines the distribution of blocks across processors
!  by call the appropriate subroutine based on distribution type
!  requested.  Currently three distributions are supported:
!  2-d Cartesian distribution (cartesian), a load-balanced
!  distribution using a rake algorithm based on an input amount of work
!  per block, and a space-filling-curve algorithm.

   character (*), intent(in) :: &
      dist_type             ! method for distributing blocks
                            !  either cartesian or rake

   integer (int_kind), intent(in) :: &
      nprocs                ! number of processors in this distribution

   integer (int_kind), dimension(:), intent(in) :: &
      work_per_block        ! amount of work per block

   type (distrb) :: &
      create_distribution   ! resulting structure describing
                            !  distribution of blocks

   character(len=*),parameter :: subname='(create_distribution)'

!----------------------------------------------------------------------
!
!  select the appropriate distribution type
!
!----------------------------------------------------------------------

   select case (trim(dist_type))

   case('cartesian')

      create_distribution = create_distrb_cart(nprocs, work_per_block)

   case('rake')

      create_distribution = create_distrb_rake(nprocs, work_per_block)

   case('roundrobin')

      create_distribution = create_distrb_roundrobin(nprocs, work_per_block)

   case('spiralcenter')

      create_distribution = create_distrb_spiralcenter(nprocs, work_per_block)

   case('wghtfile')

      create_distribution = create_distrb_wghtfile(nprocs, work_per_block)

   case('sectrobin')

      create_distribution = create_distrb_sectrobin(nprocs, work_per_block)

   case('sectcart')

      create_distribution = create_distrb_sectcart(nprocs, work_per_block)

   case('spacecurve')

      create_distribution = create_distrb_spacecurve(nprocs, work_per_block)

   case default

      call abort_ice(subname//'ERROR: ice distribution: unknown distribution type', &
         file=__FILE__, line=__LINE__)

   end select

!-----------------------------------------------------------------------

 end function create_distribution

!***********************************************************************

 subroutine create_local_block_ids(block_ids, distribution)

!  This routine determines which blocks in an input distribution are
!  located on the local processor and creates an array of block ids
!  for all local blocks.

   type (distrb), intent(in) :: &
      distribution           ! input distribution for which local
                             !  blocks required

   integer (int_kind), dimension(:), pointer :: &
      block_ids              ! array of block ids for every block
                             ! that resides on the local processor
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      n, bcount,         &! dummy counters
      istat               ! status flag for deallocate

   character(len=*),parameter :: subname='(create_local_block_ids)'

!-----------------------------------------------------------------------
!
!  first determine number of local blocks to allocate array
!
!-----------------------------------------------------------------------

   bcount = 0
   do n=1,size(distribution%blockLocation)
      if (distribution%blockLocation(n) == my_task+1) bcount = bcount + 1
   end do

!-----------------------------------------------------------------------
!
!  now fill array with proper block ids
!
!-----------------------------------------------------------------------

   if (bcount > 0) then
      allocate(block_ids(bcount), stat=istat)
      if (ice_memusage_allocErr(istat,subname//'alloc block_ids')) return
      do n=1,size(distribution%blockLocation)
         if (distribution%blockLocation(n) == my_task+1) then
            block_ids(distribution%blockLocalID(n)) = n
         endif
      end do
   endif

 end subroutine create_local_block_ids

!***********************************************************************

 subroutine proc_decomposition(nprocs, nprocs_x, nprocs_y)

!  This subroutine attempts to find an optimal (nearly square)
!  2d processor decomposition for a given number of processors.

   integer (int_kind), intent(in) :: &
      nprocs                       ! total number or processors

   integer (int_kind), intent(out) :: &
      nprocs_x, nprocs_y           ! number of procs in each dimension

!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      iguess, jguess               ! guesses for nproc_x,y

   real (real_kind) :: &
      square                       ! square root of nprocs

   character(len=*),parameter :: subname='(proc_decomposition)'

!----------------------------------------------------------------------
!
!  start with an initial guess
!
!----------------------------------------------------------------------

   square = sqrt(real(nprocs,kind=real_kind))
   nprocs_x = 0
   nprocs_y = 0

   if (processor_shape == 'square-pop') then ! make as square as possible
      iguess = nint(square)
      jguess = nprocs/iguess
   elseif (processor_shape == 'square-ice') then ! better for bipolar ice
      jguess = nint(square)
      iguess = nprocs/jguess
   elseif (processor_shape == 'slenderX1') then ! 1 proc in y direction
      jguess = 1
      iguess = nprocs/jguess
   else                                  ! 2 processors in y direction
      jguess = min(2, nprocs)
      iguess = nprocs/jguess
   endif

!----------------------------------------------------------------------
!
!  try various decompositions to find the best
!
!----------------------------------------------------------------------

   proc_loop: do
   if (processor_shape == 'square-pop') then
      jguess = nprocs/iguess
   else
      iguess = nprocs/jguess
   endif

      if (iguess*jguess == nprocs) then ! valid decomp

         !*** if the blocks can be evenly distributed, it is a
         !*** good decomposition
         if (mod(nblocks_x,iguess) == 0 .and. &
             mod(nblocks_y,jguess) == 0) then
            nprocs_x = iguess
            nprocs_y = jguess
            exit proc_loop

         !*** if the blocks can be evenly distributed in a
         !*** transposed direction, it is a good decomposition
         else if (mod(nblocks_x,jguess) == 0 .and. &
                mod(nblocks_y,iguess) == 0) then
            nprocs_x = jguess
            nprocs_y = iguess
            exit proc_loop

         !*** A valid decomposition, but keep searching for
         !***  a better one
         else
            if (nprocs_x == 0) then
               nprocs_x = iguess
               nprocs_y = jguess
            endif
            if (processor_shape == 'square-pop') then
               iguess = iguess - 1
               if (iguess == 0) then
                  exit proc_loop
               else
                  cycle proc_loop
               endif
            else
               jguess = jguess - 1
               if (jguess == 0) then
                  exit proc_loop
               else
                  cycle proc_loop
               endif
            endif
         endif

      else ! invalid decomp - keep trying

         if (processor_shape == 'square-pop') then
            iguess = iguess - 1
            if (iguess == 0) then
               exit proc_loop
            else
               cycle proc_loop
            endif
         else
            jguess = jguess - 1
            if (jguess == 0) then
               exit proc_loop
            else
               cycle proc_loop
            endif
         endif
      endif

   end do proc_loop

   if (nprocs_x == 0) then
      call abort_ice(subname//'ERROR: Unable to find 2d processor config', &
         file=__FILE__, line=__LINE__)
   endif

   if (my_task == master_task) then
     write(nu_diag,'(a,a23,i4,a3,i4)') subname,'  Processors (X x Y) = ', &
                                        nprocs_x,' x ',nprocs_y
   endif

!----------------------------------------------------------------------

 end subroutine proc_decomposition

!**********************************************************************

 subroutine ice_distributionDestroy(distribution)

!  This routine destroys a defined distribution by deallocating
!  all memory associated with the distribution.

   type (distrb), intent(inout) :: &
      distribution          ! distribution to destroy

!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: istat  ! status flag for deallocate

   character(len=*),parameter :: subname='(ice_distributionDestroy)'

!----------------------------------------------------------------------
!
!  reset scalars
!
!----------------------------------------------------------------------

   distribution%nprocs         = 0
   distribution%communicator   = 0
   distribution%numLocalBlocks = 0

!----------------------------------------------------------------------
!
!  deallocate arrays
!
!----------------------------------------------------------------------

   deallocate(distribution%blockLocation, stat=istat)
   if (ice_memusage_allocErr(istat,subname//'dealloc blockLocation')) return

   deallocate(distribution%blockLocalID , stat=istat)
   if (ice_memusage_allocErr(istat,subname//'dealloc blockLocalID')) return

   deallocate(distribution%blockGlobalID, stat=istat)
   if (ice_memusage_allocErr(istat,subname//'dealloc blockGlobalID')) return

   deallocate(distribution%blockCnt     , stat=istat)
   if (ice_memusage_allocErr(istat,subname//'dealloc blockCnt')) return

!-----------------------------------------------------------------------

 end subroutine ice_distributionDestroy

!***********************************************************************

 subroutine ice_distributionGet(distribution,&
                            nprocs, communicator, numLocalBlocks, &
                            blockLocation, blockLocalID, blockGlobalID)

!  This routine extracts information from a distribution.

   type (distrb), intent(in) :: &
      distribution        ! input distribution for which information
                          !  is requested

   integer (int_kind), intent(out), optional ::   &
      nprocs            ,&! number of processors in this dist
      communicator      ,&! communicator to use in this dist
      numLocalBlocks      ! number of blocks distributed to this
                          !   local processor

   integer (int_kind), dimension(:), optional :: &
      blockLocation     ,&! processor location for all blocks
      blockLocalID      ,&! local  block id for all blocks
      blockGlobalID       ! global block id for each local block

   character(len=*),parameter :: subname='(ice_distributionGet)'

!-----------------------------------------------------------------------
!
!  depending on which optional arguments are present, extract the
!  requested info
!
!-----------------------------------------------------------------------

   if (present(nprocs))       nprocs       = distribution%nprocs
   if (present(communicator))   communicator   = distribution%communicator
   if (present(numLocalBlocks)) numLocalBlocks = distribution%numLocalBlocks

   if (present(blockLocation)) then
      if (associated(distribution%blockLocation)) then
         blockLocation = distribution%blockLocation
      else
         call abort_ice(subname//'ERROR: blockLocation not allocated', &
            file=__FILE__, line=__LINE__)
         return
      endif
   endif

   if (present(blockLocalID)) then
      if (associated(distribution%blockLocalID)) then
         blockLocalID = distribution%blockLocalID
      else
         call abort_ice(subname//'ERROR: blockLocalID not allocated', &
            file=__FILE__, line=__LINE__)
         return
      endif
   endif

   if (present(blockGlobalID)) then
      if (associated(distribution%blockGlobalID)) then
         blockGlobalID = distribution%blockGlobalID
      else
         call abort_ice(subname//'ERROR: blockGlobalID not allocated', &
            file=__FILE__, line=__LINE__)
         return
      endif
   endif

!-----------------------------------------------------------------------

 end subroutine ice_distributionGet

!***********************************************************************

 subroutine ice_distributionGetBlockLoc(distribution, blockID, &
                                        processor, localID)

!  Given a distribution of blocks and a global block ID, return
!  the processor and local index for the block.  A zero for both
!  is returned in the case that the block has been eliminated from
!  the distribution (i.e. has no active points).

   type (distrb), intent(in) :: &
      distribution           ! input distribution for which information
                             !  is requested

   integer (int_kind), intent(in) :: &
      blockID                ! global block id for which location requested

   integer (int_kind), intent(out) ::  &
      processor,            &! processor on which block resides
      localID                ! local index for this block on this proc

   character(len=*),parameter :: subname='(ice_distributionGetBlockLoc)'

!-----------------------------------------------------------------------
!
!  check for valid blockID
!
!-----------------------------------------------------------------------

   if (blockID < 0 .or. blockID > nblocks_tot) then
      call abort_ice(subname//'ERROR: invalid block id', &
         file=__FILE__, line=__LINE__)
      return
   endif

!-----------------------------------------------------------------------
!
!  extract the location from the distribution data structure
!
!-----------------------------------------------------------------------

   processor = distribution%blockLocation(blockID)
   localID   = distribution%blockLocalID (blockID)

!-----------------------------------------------------------------------

 end subroutine ice_distributionGetBlockLoc

!***********************************************************************

 subroutine ice_distributionGetBlockID(distribution, localID, &
                                       blockID)

!  Given a distribution of blocks and a local block index, return
!  the global block id for the block.

   type (distrb), intent(in) :: &
      distribution           ! input distribution for which information
                             !  is requested

   integer (int_kind), intent(in) ::  &
      localID                ! local index for this block on this proc

   integer (int_kind), intent(out) :: &
      blockID                ! global block id for this local block

   character(len=*),parameter :: subname='(ice_distributionGetBlockID)'

!-----------------------------------------------------------------------
!
!  check for valid localID
!
!-----------------------------------------------------------------------

   if (localID < 0 .or. localID > distribution%numLocalBlocks) then
      call abort_ice(subname//'ERROR: invalid local id', &
         file=__FILE__, line=__LINE__)
      return
   endif

!-----------------------------------------------------------------------
!
!  extract the global ID from the distribution data structure
!
!-----------------------------------------------------------------------

   blockID   = distribution%blockGlobalID (localID)

!-----------------------------------------------------------------------

 end subroutine ice_distributionGetBlockID

!***********************************************************************

 function create_distrb_cart(nprocs, workPerBlock, max_blocks_calc) result(newDistrb)

!  This function creates a distribution of blocks across processors
!  using a 2-d Cartesian distribution.

   integer (int_kind), intent(in) :: &
      nprocs            ! number of processors in this distribution

   integer (int_kind), dimension(:), intent(in) :: &
      workPerBlock      ! amount of work per block

   logical (log_kind), optional :: &
      max_blocks_calc   ! compute max_blocks (default true)

   type (distrb) :: &
      newDistrb         ! resulting structure describing Cartesian
                        !  distribution of blocks

!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      i, j, n,               &! dummy loop indices
      istat,                 &! status flag for allocation
      iblock, jblock,        &!
      is, ie, js, je,        &! start, end block indices for each proc
      processor,             &! processor position in cartesian decomp
      globalID,              &! global block ID
      localID,               &! block location on this processor
      nprocsX,               &! num of procs in x for global domain
      nprocsY,               &! num of procs in y for global domain
      numBlocksXPerProc,     &! num of blocks per processor in x
      numBlocksYPerProc,     &! num of blocks per processor in y
      numBlocksPerProc        ! required number of blocks per processor

   logical (log_kind) :: &
      lmax_blocks_calc        ! local max_blocks_calc setting

   character(len=*),parameter :: subname='(create_distrb_cart)'

!----------------------------------------------------------------------

   lmax_blocks_calc = .true.
   if (present(max_blocks_calc)) then
      lmax_blocks_calc = max_blocks_calc
   endif

!----------------------------------------------------------------------
!
!  create communicator for this distribution
!
!----------------------------------------------------------------------

   call create_communicator(newDistrb%communicator, nprocs)

!----------------------------------------------------------------------
!
!  try to find best processor arrangement
!
!----------------------------------------------------------------------

   newDistrb%nprocs = nprocs

   call proc_decomposition(nprocs, nprocsX, nprocsY)

!----------------------------------------------------------------------
!
!  allocate space for decomposition
!
!----------------------------------------------------------------------

   allocate(newDistrb%blockLocation(nblocks_tot), &
            newDistrb%blockLocalID (nblocks_tot), stat=istat)
   if (ice_memusage_allocErr(istat,subname//'alloc blockLocation or blockLocalID')) return

   allocate(newDistrb%blockCnt(nprocs), stat=istat)
   if (ice_memusage_allocErr(istat,subname//'alloc blockLocation or blockLocalID')) return

!----------------------------------------------------------------------
!
!  distribute blocks linearly across processors in each direction
!
!----------------------------------------------------------------------

   numBlocksXPerProc = (nblocks_x-1)/nprocsX + 1
   numBlocksYPerProc = (nblocks_y-1)/nprocsY + 1

   newDistrb%blockCnt(:) = 0
   do j=1,nprocsY
   do i=1,nprocsX
      processor = (j-1)*nprocsX + i      ! number the processors
                                         ! left to right, bot to top

      is = (i-1)*numBlocksXPerProc + 1   ! starting block in i
      ie =  i   *numBlocksXPerProc       ! ending   block in i
      if (ie > nblocks_x) ie = nblocks_x
      js = (j-1)*numBlocksYPerProc + 1   ! starting block in j
      je =  j   *numBlocksYPerProc       ! ending   block in j
      if (je > nblocks_y) je = nblocks_y

      do jblock = js,je
      do iblock = is,ie
         globalID = (jblock - 1)*nblocks_x + iblock
         if (workPerBlock(globalID) /= 0) then
            newDistrb%blockCnt(processor) = newDistrb%blockCnt(processor) + 1
            localID = newDistrb%blockCnt(processor)
            newDistrb%blockLocation(globalID) = processor
            newDistrb%blockLocalID (globalID) = localID
         else  ! no work - eliminate block from distribution
            newDistrb%blockLocation(globalID) = 0
            newDistrb%blockLocalID (globalID) = 0
         endif
      end do
      end do

   end do
   end do

   newDistrb%numLocalBlocks = newDistrb%blockCnt(my_task+1)

   ! set local blockGlobalID array
   allocate(newDistrb%blockGlobalID(newDistrb%numLocalBlocks), stat=istat)
   if (ice_memusage_allocErr(istat,subname//'alloc numLocalBlocks')) return
   do n = 1,nblocks_tot
      if (my_task+1 == newDistrb%blockLocation(n)) then
         localID = newDistrb%blockLocalID(n)
         newDistrb%blockGlobalID (localID) = n
      endif
   enddo

   ! set/check max_blocks
   if (lmax_blocks_calc) then
      if (max_blocks < 0) then
         max_blocks = newDistrb%numLocalBlocks
      endif
   endif

!----------------------------------------------------------------------

 end function create_distrb_cart

!**********************************************************************

 function create_distrb_rake(nprocs, workPerBlock) result(newDistrb)

!  This  function distributes blocks across processors in a
!  load-balanced manner based on the amount of work per block.
!  A rake algorithm is used in which the blocks are first distributed
!  in a Cartesian distribution and then a rake is applied in each
!  Cartesian direction.

   integer (int_kind), intent(in) :: &
      nprocs                ! number of processors in this distribution

   integer (int_kind), dimension(:), intent(in) :: &
      workPerBlock        ! amount of work per block

   type (distrb) :: &
      newDistrb           ! resulting structure describing
                          ! load-balanced distribution of blocks

!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) ::  &
      i, j, n,            &! dummy loop indices
      processor,          &! dummy for processor id
      istat,              &! status flag for allocates
      globalID,           &! global block ID
      localID,            &! block location on this processor
      numOcnBlocks,       &! number of ocean blocks
      maxWork,            &! max amount of work in any block
      nprocsX,            &! num of procs in x for global domain
      nprocsY              ! num of procs in y for global domain

   integer (int_kind), dimension(:), allocatable :: &
      priority,           &! priority for moving blocks
      workTmp,            &! work per row or column for rake algrthm
      procTmp              ! temp processor id for rake algrthm

   type (distrb) :: dist   ! temp hold distribution

   character(len=*),parameter :: subname='(create_distrb_rake)'

!----------------------------------------------------------------------
!
!  first set up as Cartesian distribution
!
!----------------------------------------------------------------------

   ! ignore max_block calc in create_distrb_cart and recompute below
   dist = create_distrb_cart(nprocs, workPerBlock, max_blocks_calc=.false.)

!----------------------------------------------------------------------
!
!  if the number of blocks is close to the number of processors,
!  only do a 1-d rake on the entire distribution
!
!----------------------------------------------------------------------

   numOcnBlocks = count(workPerBlock /= 0)
   maxWork = maxval(workPerBlock)

   if (numOcnBlocks <= 2*nprocs) then
      if (my_task == master_task) &
         write(nu_diag,*) subname,' 1d rake on entire distribution'

      allocate(priority(nblocks_tot), stat=istat)
      if (ice_memusage_allocErr(istat,subname//'alloc priority')) return

      !*** initialize priority array

      do j=1,nblocks_y
      do i=1,nblocks_x
         n=(j-1)*nblocks_x + i
         if (workPerBlock(n) > 0) then
            priority(n) = maxWork + n - workPerBlock(n)
         else
            priority(n) = 0
         endif
      end do
      end do

      allocate(workTmp(nprocs), procTmp(nprocs), stat=istat)
      if (ice_memusage_allocErr(istat,subname//'alloc procTmp')) return

      workTmp(:) = 0
      do i=1,nprocs
         procTmp(i) = i
         do n=1,nblocks_tot
            if (dist%blockLocation(n) == i) then
               workTmp(i) = workTmp(i) + workPerBlock(n)
            endif
         end do
      end do

      call ice_distributionRake (workTmp, procTmp, workPerBlock, &
                                 priority, dist)

      deallocate(workTmp, procTmp, stat=istat)
      if (ice_memusage_allocErr(istat,subname//'dealloc procTmp')) return

!----------------------------------------------------------------------
!
!  otherwise re-distribute blocks using a rake in each direction
!
!----------------------------------------------------------------------

   else
      if (my_task == master_task) &
         write(nu_diag,*) subname,' rake in each direction'

      call proc_decomposition(dist%nprocs, nprocsX, nprocsY)

!----------------------------------------------------------------------
!
!     load-balance using a rake algorithm in the x-direction first
!
!----------------------------------------------------------------------

      allocate(priority(nblocks_tot), stat=istat)
      if (ice_memusage_allocErr(istat,subname//'alloc priority')) return

      !*** set highest priority such that eastern-most blocks
      !*** and blocks with the least amount of work are
      !*** moved first

      do j=1,nblocks_y
      do i=1,nblocks_x
         n=(j-1)*nblocks_x + i
         if (workPerBlock(n) > 0) then
            priority(n) = (maxWork + 1)*(nblocks_x + i) - &
                          workPerBlock(n)
         else
            priority(n) = 0
         endif
      end do
      end do

      allocate(workTmp(nprocsX), procTmp(nprocsX), stat=istat)
      if (ice_memusage_allocErr(istat,subname//'alloc procTmp')) return

      do j=1,nprocsY

         workTmp(:) = 0
         do i=1,nprocsX
            processor = (j-1)*nprocsX + i
            procTmp(i) = processor
            do n=1,nblocks_tot
               if (dist%blockLocation(n) == processor) then
                  workTmp(i) = workTmp(i) + workPerBlock(n)
               endif
            end do
         end do

         call ice_distributionRake (workTmp, procTmp, workPerBlock, &
                                    priority, dist)
      end do

      deallocate(workTmp, procTmp, stat=istat)
      if (ice_memusage_allocErr(istat,subname//'dealloc procTmp')) return

!----------------------------------------------------------------------
!
!     use a rake algorithm in the y-direction now
!
!----------------------------------------------------------------------

      !*** set highest priority for northern-most blocks

      do j=1,nblocks_y
      do i=1,nblocks_x
         n=(j-1)*nblocks_x + i
         if (workPerBlock(n) > 0) then
            priority(n) = (maxWork + 1)*(nblocks_y + j) - &
                          workPerBlock(n)
         else
            priority(n) = 0
         endif
      end do
      end do

      allocate(workTmp(nprocsY), procTmp(nprocsY), stat=istat)
      if (ice_memusage_allocErr(istat,subname//'alloc procTmp')) return

      do i=1,nprocsX

         workTmp(:) = 0
         do j=1,nprocsY
            processor = (j-1)*nprocsX + i
            procTmp(j) = processor
            do n=1,nblocks_tot
               if (dist%blockLocation(n) == processor) then
                  workTmp(j) = workTmp(j) + workPerBlock(n)
               endif
            end do
         end do

         call ice_distributionRake (workTmp, procTmp, workPerBlock, &
                                    priority, dist)

      end do

      deallocate(workTmp, procTmp, priority, stat=istat)
      if (ice_memusage_allocErr(istat,subname//'dealloc procTmp')) return

   endif  ! 1d or 2d rake

!----------------------------------------------------------------------
!
!  create new distribution with info extracted from the temporary
!  distribution
!
!----------------------------------------------------------------------

   newDistrb%nprocs     = nprocs
   newDistrb%communicator = dist%communicator

   allocate(newDistrb%blockLocation(nblocks_tot), &
            newDistrb%blockLocalID(nblocks_tot), stat=istat)
   if (ice_memusage_allocErr(istat,subname//'alloc blockLocation or blockLocalID')) return

   allocate(newDistrb%blockCnt(nprocs), stat=istat)
   if (ice_memusage_allocErr(istat,subname//'alloc blockCnt')) return
   newDistrb%blockCnt(:) = 0

   do n=1,nblocks_tot
      globalID = n
      processor = dist%blockLocation(globalID)  ! processor id
      newDistrb%blockLocation(globalID) = processor

      if (processor > 0) then
         newDistrb%blockCnt(processor) = newDistrb%blockCnt(processor) + 1
         localID = newDistrb%blockCnt(processor)
         newDistrb%blockLocation(globalID) = processor
         newDistrb%blockLocalID (globalID) = localID
      else
         newDistrb%blockLocation(globalID) = 0
         newDistrb%blockLocalID (globalID) = 0
      endif
   end do

   newDistrb%numLocalBlocks = newDistrb%blockCnt(my_task+1)

   ! set local blockGlobalID array
   allocate(newDistrb%blockGlobalID(newDistrb%numLocalBlocks), stat=istat)
   if (ice_memusage_allocErr(istat,subname//'alloc numLocalBlocks')) return
   do n = 1,nblocks_tot
      if (my_task+1 == newDistrb%blockLocation(n)) then
         localID = newDistrb%blockLocalID(n)
         newDistrb%blockGlobalID (localID) = n
      endif
   enddo

   ! set/check max_blocks
   if (max_blocks < 0) then
      max_blocks = newDistrb%numLocalBlocks
   endif

   ! destroy cart distribution
   call ice_distributionDestroy(dist)

!----------------------------------------------------------------------

 end function create_distrb_rake

!***********************************************************************

 function create_distrb_roundrobin(nprocs, workPerBlock) result(newDistrb)

!  This function creates a distribution of blocks across processors
!  using a simple roundrobin algorithm. Mean for prescribed ice or
!  standalone CAM mode.

   integer (int_kind), intent(in) :: &
      nprocs              ! number of processors in this distribution

   integer (int_kind), dimension(:), intent(in) :: &
      workPerBlock        ! amount of work per block

   type (distrb) :: &
      newDistrb           ! resulting structure describing Cartesian
                          !  distribution of blocks

!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      i, j, n,               &! dummy loop indices
      istat,                 &! status flag for allocation
      processor,             &! processor position in cartesian decomp
      globalID,              &! global block ID
      localID                 ! block location on this processor

   character(len=*),parameter :: subname='(create_distrb_roundrobin)'

!----------------------------------------------------------------------
!
!  create communicator for this distribution
!
!----------------------------------------------------------------------

   call create_communicator(newDistrb%communicator, nprocs)

!----------------------------------------------------------------------
!
!  try to find best processor arrangement
!
!----------------------------------------------------------------------

   newDistrb%nprocs = nprocs

!----------------------------------------------------------------------
!
!  allocate space for decomposition
!
!----------------------------------------------------------------------

   allocate(newDistrb%blockLocation(nblocks_tot), &
            newDistrb%blockLocalID (nblocks_tot), stat=istat)
   if (ice_memusage_allocErr(istat,subname//'alloc blockLocation or blockLocalID')) return

   allocate(newDistrb%blockCnt(nprocs), stat=istat)
   if (ice_memusage_allocErr(istat,subname//'alloc blockCnt')) return

!----------------------------------------------------------------------
!
!  distribute blocks across processors, one block per proc until used
!
!----------------------------------------------------------------------

   processor = 0
   globalID = 0
   newDistrb%numLocalBlocks = 0
   newDistrb%blockCnt(:) = 0

   ! compute decomposition
   do j=1,nblocks_y
   do i=1,nblocks_x
      globalID = globalID + 1
      if (workPerBlock(globalID) /= 0) then
         processor = mod(processor,nprocs) + 1
         newDistrb%blockCnt(processor) = newDistrb%blockCnt(processor) + 1
         localID = newDistrb%blockCnt(processor)
         newDistrb%blockLocation(globalID) = processor
         newDistrb%blockLocalID (globalID) = localID
      else  ! no work - eliminate block from distribution
         newDistrb%blockLocation(globalID) = 0
         newDistrb%blockLocalID (globalID) = 0
      endif
   enddo
   enddo
   newDistrb%numLocalBlocks = newDistrb%blockCnt(my_task+1)

   ! set local blockGlobalID array
   allocate(newDistrb%blockGlobalID(newDistrb%numLocalBlocks), stat=istat)
   if (ice_memusage_allocErr(istat,subname//'alloc numLocalBlocks')) return
   do n = 1,nblocks_tot
      if (my_task+1 == newDistrb%blockLocation(n)) then
         localID = newDistrb%blockLocalID(n)
         newDistrb%blockGlobalID (localID) = n
      endif
   enddo

   ! set/check max_blocks
   if (max_blocks < 0) then
      max_blocks = newDistrb%numLocalBlocks
   endif

!----------------------------------------------------------------------

 end function create_distrb_roundrobin

!***********************************************************************

 function create_distrb_spiralcenter(nprocs, workPerBlock) result(newDistrb)

!  This function creates a distribution of blocks across processors
!  using a simple spiralcenter algorithm. Mean for prescribed ice or
!  standalone CAM mode.

   integer (int_kind), intent(in) :: &
      nprocs              ! number of processors in this distribution

   integer (int_kind), dimension(:), intent(in) :: &
      workPerBlock        ! amount of work per block

   type (distrb) :: &
      newDistrb           ! resulting structure describing Cartesian
                          !  distribution of blocks

!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      n, i, j, ic, jc, id, jd, cnt,  &! dummy loop indices
      istat,            &! status flag for allocation
      processor,        &! processor position in cartesian decomp
      nblocklist,       &! number of blocks in blocklist
      globalID,         &! global block ID
      localID            ! block location on this processor

   integer (int_kind), dimension(:), allocatable :: &
      blocklist          ! temp block ordered list
   integer (int_kind), dimension(:,:), allocatable :: &
      blockchk           ! temp block check array

   character(len=*),parameter :: subname='(create_distrb_spiralcenter)'

!----------------------------------------------------------------------
!
!  create communicator for this distribution
!
!----------------------------------------------------------------------

   call create_communicator(newDistrb%communicator, nprocs)

!----------------------------------------------------------------------
!
!  try to find best processor arrangement
!
!----------------------------------------------------------------------

   newDistrb%nprocs = nprocs

!----------------------------------------------------------------------
!
!  allocate space for decomposition
!
!----------------------------------------------------------------------

   allocate(newDistrb%blockLocation(nblocks_tot), &
            newDistrb%blockLocalID (nblocks_tot), stat=istat)
   if (ice_memusage_allocErr(istat,subname//'alloc blockLocation or blockLocalID')) return

   allocate(newDistrb%blockCnt(nprocs), stat=istat)
   if (ice_memusage_allocErr(istat,subname//'alloc blockCnt')) return

!----------------------------------------------------------------------
!
!  create list of blocks starting from center in spiral
!  pattern is start in center, right 1, up 1, left 2, down 2,
!  right 3, up 3, left 4, down 4, right 5, up 5, etc.
!  until all blocks have been accounted for just once.
!  cnt tracks the up, left, down, right counts and is the emergency
!  stop
!
!----------------------------------------------------------------------

   allocate(blocklist(nblocks_tot), stat=istat)
   if (ice_memusage_allocErr(istat,subname//'alloc blocklist')) return
   allocate(blockchk(nblocks_x,nblocks_y), stat=istat)
   if (ice_memusage_allocErr(istat,subname//'alloc blockchk')) return
   nblocklist = 0
   blocklist = 0
   blockchk = 0
   processor = 0
   globalID = 0

   jc = nblocks_y/2
   ic = nblocks_x/2

   ! center block
   cnt = 0
   j = jc
   i = ic
   globalID = (j-1)*nblocks_x + i
   nblocklist = nblocklist + 1
   blocklist(nblocklist) = globalID
   blockchk(i,j) = 1

   do while (minval(blocklist) < 1 .and. cnt < max(nblocks_x,nblocks_y) )

     cnt = cnt + 1

     ! right, j held constant
     ic = i
     do id = ic+1,ic+cnt
       i = max(min(id,nblocks_x),1)
       if (blockchk(i,j) == 0) then
         globalID = (j-1)*nblocks_x + i
         nblocklist = nblocklist + 1
         blocklist(nblocklist) = globalID
         blockchk(i,j) = 1
       endif
     enddo

     ! up, i held constant
     jc = j
     do jd = jc+1,jc+cnt
       j = max(min(jd,nblocks_y),1)
       if (blockchk(i,j) == 0) then
         globalID = (j-1)*nblocks_x + i
         nblocklist = nblocklist + 1
         blocklist(nblocklist) = globalID
         blockchk(i,j) = 1
       endif
     enddo

     cnt = cnt + 1

     ! left, j held constant
     ic = i
     do id = ic-1,ic-cnt,-1
       i = max(min(id,nblocks_x),1)
       if (blockchk(i,j) == 0) then
         globalID = (j-1)*nblocks_x + i
         nblocklist = nblocklist + 1
         blocklist(nblocklist) = globalID
         blockchk(i,j) = 1
       endif
     enddo

     ! down, i held constant
     jc = j
     do jd = jc-1,jc-cnt,-1
       j = max(min(jd,nblocks_y),1)
       if (blockchk(i,j) == 0) then
         globalID = (j-1)*nblocks_x + i
         nblocklist = nblocklist + 1
         blocklist(nblocklist) = globalID
         blockchk(i,j) = 1
       endif
     enddo

   enddo

   if (nblocklist /= nblocks_x*nblocks_y .or. &
       maxval(blockchk) /= 1 .or. minval(blockchk) /= 1) then
     call abort_ice(subname//'ERROR: blockchk invalid', &
         file=__FILE__, line=__LINE__)
     return
   endif
   deallocate(blockchk, stat=istat)
   if (ice_memusage_allocErr(istat,subname//'dealloc blockchk')) return

!----------------------------------------------------------------------
!
!  now distribute the blocks in the blocklist roundrobin
!
!----------------------------------------------------------------------

   newDistrb%numLocalBlocks = 0
   newDistrb%blockCnt(:) = 0

   do n = 1,nblocklist
      globalID = blocklist(n)
      if (workPerBlock(globalID) /= 0) then
         processor = mod(processor,nprocs) + 1
         newDistrb%blockCnt(processor) = newDistrb%blockCnt(processor) + 1
         localID = newDistrb%blockCnt(processor)
         newDistrb%blockLocation(globalID) = processor
         newDistrb%blockLocalID (globalID) = localID
      else  ! no work - eliminate block from distribution
         newDistrb%blockLocation(globalID) = 0
         newDistrb%blockLocalID (globalID) = 0
      endif
   end do
   newDistrb%numLocalBlocks = newDistrb%blockCnt(my_task+1)

   ! set local blockGlobalID array
   allocate(newDistrb%blockGlobalID(newDistrb%numLocalBlocks), stat=istat)
   if (ice_memusage_allocErr(istat,subname//'alloc numLocalBlocks')) return
   do n = 1,nblocks_tot
      if (my_task+1 == newDistrb%blockLocation(n)) then
         localID = newDistrb%blockLocalID(n)
         newDistrb%blockGlobalID (localID) = n
      endif
   enddo

   ! set/check max_blocks
   if (max_blocks < 0) then
      max_blocks = newDistrb%numLocalBlocks
   endif

   deallocate(blocklist, stat=istat)
   if (ice_memusage_allocErr(istat,subname//'dealloc blocklist')) return

!----------------------------------------------------------------------

 end function create_distrb_spiralcenter

!***********************************************************************

 function create_distrb_wghtfile(nprocs, workPerBlock) result(newDistrb)

!  This function creates a distribution of blocks across processors
!  using a simple wghtfile algorithm. Meant for prescribed ice or
!  standalone CAM mode.

   integer (int_kind), intent(in) :: &
      nprocs              ! number of processors in this distribution

   integer (int_kind), dimension(:), intent(in) :: &
      workPerBlock        ! amount of work per block

   type (distrb) :: &
      newDistrb           ! resulting structure describing Cartesian
                          !  distribution of blocks

!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      i, j, n,               &! dummy loop indices
      cnt,                   &! counter
      istat,                 &! status flag for allocation
      processor,             &! processor position in cartesian decomp
      globalID,              &! global block ID
      localID                 ! block location on this processor

   logical (log_kind) ::  up   ! direction of pe counting

   character(len=*),parameter :: subname='(create_distrb_wghtfile)'

!----------------------------------------------------------------------
!
!  create communicator for this distribution
!
!----------------------------------------------------------------------

   call create_communicator(newDistrb%communicator, nprocs)

!----------------------------------------------------------------------
!
!  try to find best processor arrangement
!
!----------------------------------------------------------------------

   newDistrb%nprocs = nprocs

!----------------------------------------------------------------------
!
!  allocate space for decomposition
!
!----------------------------------------------------------------------

   allocate(newDistrb%blockLocation(nblocks_tot), &
            newDistrb%blockLocalID (nblocks_tot), stat=istat)
   if (ice_memusage_allocErr(istat,subname//'alloc blockLocation or blockLocalID')) return

   allocate(newDistrb%blockCnt(nprocs), stat=istat)
   if (ice_memusage_allocErr(istat,subname//'alloc blockCnt')) return

!----------------------------------------------------------------------
!
!  distribute blocks across processors, one block per proc until used
!  work from most expensive workPerBlock to least and go up/down/up/down
!  in terms of the pe index to try to get better load balance.
!
!----------------------------------------------------------------------

   processor = 0
   newDistrb%numLocalBlocks = 0
   newDistrb%blockCnt(:) = 0
   up = .true.

   if (my_task == master_task) &
      write(nu_diag,*) subname,' workPerBlock = ',minval(workPerBlock),maxval(workPerBlock)
   if (minval(workPerBlock) < 0 .or. maxval(workPerBlock) > 12) then
      write(nu_diag,*) subname,' workPerBlock = ',minval(workPerBlock),maxval(workPerBlock)
      call abort_ice(subname//'ERROR: workPerBlock incorrect', &
         file=__FILE__, line=__LINE__)
      return
   endif

   ! do not distribution blocks with work=0
   do n = maxval(workPerBlock),1,-1
      cnt = 0
      do j=1,nblocks_y
      do i=1,nblocks_x
         if (mod(j,2) == 1) then
            globalID = (j-1)*nblocks_x + i
         else
            globalID = (j-1)*nblocks_x + nblocks_x - i + 1
         endif
         if (workPerBlock(globalID) == 0) then  ! no work - eliminate block from distribution
            newDistrb%blockLocation(globalID) = 0
            newDistrb%blockLocalID (globalID) = 0
         elseif (workPerBlock(globalID) == n) then
            cnt = cnt + 1
            if (up) then
               processor = processor + 1
            else
               processor = processor - 1
            endif
            if (processor > nprocs) then
               up = .false.
               processor = nprocs
            elseif (processor < 1) then
               up = .true.
               processor = 1
            endif
            newDistrb%blockCnt(processor) = newDistrb%blockCnt(processor) + 1
            localID = newDistrb%blockCnt(processor)
            newDistrb%blockLocation(globalID) = processor
            newDistrb%blockLocalID (globalID) = localID
         endif
      end do
      end do
!      write(nu_diag,*) subname,'n cnt = ',n,cnt
   end do
   newDistrb%numLocalBlocks = newDistrb%blockCnt(my_task+1)

   ! set local blockGlobalID array
   allocate(newDistrb%blockGlobalID(newDistrb%numLocalBlocks), stat=istat)
   if (ice_memusage_allocErr(istat,subname//'alloc numLocalBlocks')) return
   do n = 1,nblocks_tot
      if (my_task+1 == newDistrb%blockLocation(n)) then
         localID = newDistrb%blockLocalID(n)
         newDistrb%blockGlobalID (localID) = n
      endif
   enddo

   ! set/check max_blocks
   if (max_blocks < 0) then
      max_blocks = newDistrb%numLocalBlocks
   endif

!   write(nu_diag,*) subname,'my_task,newDistrb%numLocalBlocks',&
!      my_task,newDistrb%numLocalBlocks

!----------------------------------------------------------------------

 end function create_distrb_wghtfile

!***********************************************************************

 function create_distrb_sectrobin(nprocs, workPerBlock) result(newDistrb)

!  This function creates a distribution of blocks across processors
!  using a simple sectrobin algorithm. Mean for prescribed ice or
!  standalone CAM mode.

   integer (int_kind), intent(in) :: &
      nprocs              ! number of processors in this distribution

   integer (int_kind), dimension(:), intent(in) :: &
      workPerBlock        ! amount of work per block

   type (distrb) :: &
      newDistrb           ! resulting structure describing Cartesian
                          !  distribution of blocks

!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      i, j, n,               &! dummy loop indices
      istat,                 &! status flag for allocation
      mblocks,               &! estimate of max blocks per pe
      processor,             &! processor position in cartesian decomp
      globalID,              &! global block ID
      localID                 ! block location on this processor

   logical (log_kind), dimension(:), allocatable :: &
      bfree              ! map of assigned blocks, true = free

   integer (int_kind) :: cnt, blktogether, i2
   integer (int_kind) :: totblocks, nchunks
   logical (log_kind) :: keepgoing

   character(len=*),parameter :: subname='(create_distrb_sectrobin)'

!----------------------------------------------------------------------
!
!  create communicator for this distribution
!
!----------------------------------------------------------------------

   call create_communicator(newDistrb%communicator, nprocs)

!----------------------------------------------------------------------
!
!  try to find best processor arrangement
!
!----------------------------------------------------------------------

   newDistrb%nprocs = nprocs

!----------------------------------------------------------------------
!
!  allocate space for decomposition
!
!----------------------------------------------------------------------

   allocate(newDistrb%blockLocation(nblocks_tot), &
            newDistrb%blockLocalID (nblocks_tot), stat=istat)
   if (ice_memusage_allocErr(istat,subname//'alloc blockLocation or blockLocalID')) return

   allocate(newDistrb%blockCnt(nprocs), stat=istat)
   if (ice_memusage_allocErr(istat,subname//'alloc blockCnt')) return

!----------------------------------------------------------------------
!
!  distribute groups of blocks across processors, one per proc until used
!
!----------------------------------------------------------------------

   processor = 0
   globalID = 0
   newDistrb%numLocalBlocks = 0
   newDistrb%blockCnt(:) = 0
   allocate(bfree(nblocks_x*nblocks_y), stat=istat)
   if (ice_memusage_allocErr(istat,subname//'alloc bfree')) return
   bfree=.true.

   totblocks = 0
   do j=1,nblocks_y
   do i=1,nblocks_x
      globalID = (j-1)*nblocks_x + i
      if (workPerBlock(globalID) /= 0) then
         totblocks=totblocks+1
      else  ! no work - eliminate block from distribution
         bfree(globalID) = .false.
         newDistrb%blockLocation(globalID) = 0
         newDistrb%blockLocalID (globalID) = 0
      endif
   enddo
   enddo

   mblocks = totblocks/nprocs
   if (mod(totblocks,nprocs) > 0) mblocks=mblocks+1

   blktogether = max(1,nint(float(totblocks)/float(6*nprocs)))

!   write(nu_diag,*) subname,'totblocks = ',totblocks,nblocks_y*nblocks_x

   !------------------------------
   ! southern group of blocks
   !   weave back and forth in i vs j
   !   go south to north, low - high pes
   !   keepgoing to false to stop distribution
   !------------------------------

   processor=1
   cnt = 0
   keepgoing = .true.
   do j=1,nblocks_y
   do i=1,nblocks_x
      if (mod(j,2) == 0) then
         i2 = nblocks_x - i + 1
      else
         i2 = i
      endif
      globalID = (j-1)*nblocks_x + i2
      if (cnt >= blktogether) then
         processor = mod(processor,nprocs) + 1
         cnt = 0
         if (processor == 1) keepgoing = .false.
      endif
!      write(nu_diag,'(a,6i7,l2)') subname,i,j,globalID,cnt,blktogether,processor,keepgoing

      if (keepgoing) then
         if (bfree(globalID)) then
         if (workPerBlock(globalID) /= 0) then
            newDistrb%blockCnt(processor) = newDistrb%blockCnt(processor) + 1
            localID = newDistrb%blockCnt(processor)
            newDistrb%blockLocation(globalID) = processor
            newDistrb%blockLocalID (globalID) = localID
            cnt = cnt + 1
            totblocks = totblocks-1
            bfree(globalID) = .false.
         else  ! no work - eliminate block from distribution
            bfree(globalID) = .false.
            newDistrb%blockLocation(globalID) = 0
            newDistrb%blockLocalID (globalID) = 0
         endif
         endif  ! bfree
      endif
   end do
   end do

!   write(nu_diag,*) subname,'totblocks left after southern = ',totblocks

   !------------------------------
   ! northern group of blocks
   !   weave back and forth in i vs j
   !   go north to south, high - low pes
   !   keepgoing to false to stop distribution
   !------------------------------

   processor=nprocs
   cnt = 0
   keepgoing = .true.
   do j=nblocks_y,1,-1
   do i=1,nblocks_x
      if (mod(j,2) == 1) then
         i2 = nblocks_x - i + 1
      else
         i2 = i
      endif
      globalID = (j-1)*nblocks_x + i2
      if (cnt >= blktogether) then
         processor = mod(processor+nprocs-2,nprocs) + 1
         cnt = 0
         if (processor == nprocs) keepgoing = .false.
      endif

      if (keepgoing) then
         if (bfree(globalID)) then
         if (workPerBlock(globalID) /= 0) then
            newDistrb%blockCnt(processor) = newDistrb%blockCnt(processor) + 1
            localID = newDistrb%blockCnt(processor)
            newDistrb%blockLocation(globalID) = processor
            newDistrb%blockLocalID (globalID) = localID
            cnt = cnt + 1
            totblocks = totblocks-1
            bfree(globalID) = .false.
         else  ! no work - eliminate block from distribution
            bfree(globalID) = .false.
            newDistrb%blockLocation(globalID) = 0
            newDistrb%blockLocalID (globalID) = 0
         endif
         endif  ! bfree
      endif
   end do
   end do

!   write(nu_diag,*) subname,'totblocks left after northern = ',totblocks

   !------------------------------
   ! central group of blocks
   !   weave back and forth in i vs j
   !   go north to south, low - high / low - high pes
   !   distribute rest of blocks in 2 chunks per proc
   !------------------------------

   nchunks = 2*nprocs
   blktogether = max(1,nint(float(totblocks)/float(nchunks)))
   processor=1
   cnt = 0
   do j=nblocks_y,1,-1
   do i=1,nblocks_x
      if (mod(j,2) == 1) then
         i2 = nblocks_x - i + 1
      else
         i2 = i
      endif
      globalID = (j-1)*nblocks_x + i2
      if (totblocks > 0) then
         do while (newDistrb%blockCnt(processor) >= mblocks .or. cnt >= blktogether)
            nchunks = nchunks - 1
            if (nchunks == 0) then
               blktogether = 1
            else
               blktogether = max(1,nint(float(totblocks)/float(nchunks)))
            endif
            cnt = 0
            processor = mod(processor,nprocs) + 1
         enddo
      endif

!      write(nu_diag,*) subname,'central ',i,j,totblocks,cnt,nchunks,blktogether,processor

      if (bfree(globalID)) then
      if (workPerBlock(globalID) /= 0) then
         newDistrb%blockCnt(processor) = newDistrb%blockCnt(processor) + 1
         localID = newDistrb%blockCnt(processor)
         newDistrb%blockLocation(globalID) = processor
         newDistrb%blockLocalID (globalID) = localID
         cnt = cnt + 1
         totblocks = totblocks-1
         bfree(globalID) = .false.
      else  ! no work - eliminate block from distribution
         bfree(globalID) = .false.
         newDistrb%blockLocation(globalID) = 0
         newDistrb%blockLocalID (globalID) = 0
      endif
      endif  ! bfree
   end do
   end do

   newDistrb%numLocalBlocks = newDistrb%blockCnt(my_task+1)

   ! set local blockGlobalID array
   allocate(newDistrb%blockGlobalID(newDistrb%numLocalBlocks), stat=istat)
   if (ice_memusage_allocErr(istat,subname//'alloc numLocalBlocks')) return
   do n = 1,nblocks_tot
      if (my_task+1 == newDistrb%blockLocation(n)) then
         localID = newDistrb%blockLocalID(n)
         newDistrb%blockGlobalID (localID) = n
      endif
   enddo

   ! set/check max_blocks
   if (max_blocks < 0) then
      max_blocks = newDistrb%numLocalBlocks
   endif

   deallocate(bfree, stat=istat)
   if (ice_memusage_allocErr(istat,subname//'dealloc bfree')) return

!----------------------------------------------------------------------

 end function create_distrb_sectrobin

!***********************************************************************

 function create_distrb_sectcart(nprocs, workPerBlock) result(newDistrb)

!  This function creates a distribution of blocks across processors
!  using a simple sectcart algorithm. Mean for prescribed ice or
!  standalone CAM mode.

   integer (int_kind), intent(in) :: &
      nprocs              ! number of processors in this distribution

   integer (int_kind), dimension(:), intent(in) :: &
      workPerBlock        ! amount of work per block

   type (distrb) :: &
      newDistrb           ! resulting structure describing Cartesian
                          !  distribution of blocks

!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      i, j, i2, j2,          &! dummy loop indices
      istat,                 &! status flag for allocation
      processor,             &! processor position in cartesian decomp
      globalID,              &! global block ID
      localID,               &! block location on this processor
      blktogether,           &! number of blocks together
      cnt                     ! counter

   integer (int_kind) :: n

   character(len=*),parameter :: subname='(create_distrb_sectcart)'

!----------------------------------------------------------------------
!
!  create communicator for this distribution
!
!----------------------------------------------------------------------

   call create_communicator(newDistrb%communicator, nprocs)

!----------------------------------------------------------------------
!
!  try to find best processor arrangement
!
!----------------------------------------------------------------------

   newDistrb%nprocs = nprocs

!----------------------------------------------------------------------
!
!  allocate space for decomposition
!
!----------------------------------------------------------------------

   allocate(newDistrb%blockLocation(nblocks_tot), &
            newDistrb%blockLocalID (nblocks_tot), stat=istat)
   if (ice_memusage_allocErr(istat,subname//'alloc blockLocation or blockLocalID')) return

   allocate(newDistrb%blockCnt(nprocs), stat=istat)
   if (ice_memusage_allocErr(istat,subname//'alloc blockCnt')) return

!----------------------------------------------------------------------
!
!  distribute blocks linearly across processors in quadrants
!
!----------------------------------------------------------------------

   blktogether = max(1,nint(float(nblocks_x*nblocks_y)/float(4*nprocs)))

   ! --- two phases, reset processor and cnt for each phase
   ! --- phase 1 is south to north, east to west on the left half of the domain
   ! --- phase 2 is north to south, east to west on the right half of the domain

   if (mod(nblocks_x,2) /= 0) then
      call abort_ice(subname//'ERROR: nblocks_x not divisible by 2', &
         file=__FILE__, line=__LINE__)
      return
   endif

   newDistrb%numLocalBlocks = 0
   newDistrb%blockCnt(:) = 0

   do n=1,2
   processor = 1
   cnt = 0
   do j2=1,nblocks_y
   do i2=1,nblocks_x/2

      if (n == 1) then
         i = i2
         j = j2
      else
         i = nblocks_x/2 + i2
         j = nblocks_y - j2 + 1
      endif

      globalID = (j-1)*nblocks_x + i
      if (cnt >= blktogether) then
         processor = mod(processor,nprocs) + 1
         cnt = 0
      endif
      cnt = cnt + 1

      if (workPerBlock(globalID) /= 0) then
         newDistrb%blockCnt(processor) = newDistrb%blockCnt(processor) + 1
         localID = newDistrb%blockCnt(processor)
         newDistrb%blockLocation(globalID) = processor
         newDistrb%blockLocalID (globalID) = localID
      else  ! no work - eliminate block from distribution
         newDistrb%blockLocation(globalID) = 0
         newDistrb%blockLocalID (globalID) = 0
      endif

   end do
   end do
   end do
   newDistrb%numLocalBlocks = newDistrb%blockCnt(my_task+1)

   ! set local blockGlobalID array
   allocate(newDistrb%blockGlobalID(newDistrb%numLocalBlocks), stat=istat)
   if (ice_memusage_allocErr(istat,subname//'alloc numLocalBlocks')) return
   do n = 1,nblocks_tot
      if (my_task+1 == newDistrb%blockLocation(n)) then
         localID = newDistrb%blockLocalID(n)
         newDistrb%blockGlobalID (localID) = n
      endif
   enddo

   ! set/check max_blocks
   if (max_blocks < 0) then
      max_blocks = newDistrb%numLocalBlocks
   endif

!----------------------------------------------------------------------

 end function create_distrb_sectcart

!**********************************************************************

 function create_distrb_spacecurve(nprocs,work_per_block) result(newDistrb)

!  This function distributes blocks across processors in a
!  load-balanced manner using space-filling curves
!  added by J. Dennis 3/10/06

   use ice_spacecurve

   integer (int_kind), intent(in) :: &
      nprocs            ! number of processors in this distribution

   integer (int_kind), dimension(:), intent(in) :: &
      work_per_block    ! amount of work per block

   type (distrb) :: &
      newDistrb         ! resulting structure describing Cartesian
                        !  distribution of blocks

!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      i, j, n,            &! dummy loop indices
      istat,              &! status flag for allocation
      processor,          &! processor position in cartesian decomp
      globalID,           &! global block ID
      localID              ! local block position on processor

   integer (int_kind), dimension(:),allocatable :: &
        idxT_i,idxT_j      ! Temporary indices for SFC

   integer (int_kind), dimension(:,:),allocatable :: &
        Mesh,             &!   !arrays to hold Space-filling curve
        Mesh2,            &!
        Mesh3              !

   integer (int_kind) :: &
        nblocksL,nblocks, &! Number of blocks local and total
        ii,extra,tmp1,    &! loop tempories used for
        s1,ig              ! partitioning curve

   type (factor_t) :: xdim,ydim

   integer (int_kind) :: it,jj,i2,j2
   integer (int_kind) :: curveSize,sb_x,sb_y,itmp,numfac
   integer (int_kind) :: subNum, sfcNum
   logical            :: foundx

   character(len=*),parameter :: subname='(create_distrb_spacecurve)'

   !------------------------------------------------------
   ! Space filling curves only work if:
   !
   !    nblocks_x = nblocks_y
   !       nblocks_x = 2^m 3^n 5^p where m,n,p are integers
   !------------------------------------------------------

   if((.not. IsFactorable(nblocks_y)) .or. (.not. IsFactorable(nblocks_x))) then
     newDistrb = create_distrb_cart(nprocs, work_per_block)
     return
   endif

!----------------------------------------------------------------------
!
!  create communicator for this distribution
!
!----------------------------------------------------------------------

   call create_communicator(newDistrb%communicator, nprocs)

!----------------------------------------------------------------------
!
!  try to find best processor arrangement
!
!----------------------------------------------------------------------

   newDistrb%nprocs = nprocs

!----------------------------------------------------------------------
!
!  allocate space for decomposition
!
!----------------------------------------------------------------------

   allocate(newDistrb%blockLocation(nblocks_tot), &
            newDistrb%blockLocalID (nblocks_tot), stat=istat)
   if (ice_memusage_allocErr(istat,subname//'alloc blockLocation or blockLocalID')) return

   allocate(newDistrb%blockCnt(nprocs), stat=istat)
   if (ice_memusage_allocErr(istat,subname//'alloc blockCnt')) return

   !-----------------------------------------------
   ! Factor the numbers of blocks in each dimension
   !-----------------------------------------------

   xdim = Factor(nblocks_x)
   ydim = Factor(nblocks_y)
   numfac = xdim%numfact

   !---------------------------------------------
   ! Match the common factors to create SFC curve
   !---------------------------------------------

   curveSize=1
   do it=1,numfac
      call MatchFactor(xdim,ydim,itmp,foundX)
      curveSize = itmp*curveSize
   enddo

   !--------------------------------------
   ! determine the size of the sub-blocks
   ! within the space-filling curve
   !--------------------------------------

   sb_x = ProdFactor(xdim)
   sb_y = ProdFactor(ydim)

   !----------------------------------------------------------------------
   !  Create the array to hold the SFC and indices into it
   !----------------------------------------------------------------------

   allocate(Mesh(curveSize,curveSize),  &
            Mesh2(nblocks_x,nblocks_y), &
            Mesh3(nblocks_x,nblocks_y), &
            idxT_i(nblocks_tot),        &
            idxT_j(nblocks_tot), stat=istat)
   if (ice_memusage_allocErr(istat,subname//'alloc meshes')) return

   Mesh  = 0
   Mesh2 = 0
   Mesh3 = 0

   !----------------------------------------------------------------------
   !  Generate the space-filling curve
   !----------------------------------------------------------------------

   call GenSpaceCurve(Mesh)
   Mesh = Mesh + 1 ! make it 1-based indexing
!   if (debug_blocks) then
!     if (my_task == master_task) call PrintCurve(Mesh)
!   endif

   !-----------------------------------------------
   ! Reindex the SFC to address internal sub-blocks
   !-----------------------------------------------

   do j=1,curveSize
   do i=1,curveSize
      sfcNum = (Mesh(i,j) - 1)*(sb_x*sb_y) + 1
      do jj=1,sb_y
      do ii=1,sb_x
         subNum = (jj-1)*sb_x + (ii-1)
         i2 = (i-1)*sb_x + ii
         j2 = (j-1)*sb_y + jj
         Mesh2(i2,j2) = sfcNum + subNum
      enddo
      enddo
   enddo
   enddo

   !------------------------------------------------
   ! create a linear array of i,j coordinates of SFC
   !------------------------------------------------

   idxT_i=0;idxT_j=0
   do j=1,nblocks_y
     do i=1,nblocks_x
        n = (j-1)*nblocks_x + i
        ig = Mesh2(i,j)
        if(work_per_block(n) /= 0) then
            idxT_i(ig)=i;idxT_j(ig)=j
        endif
     enddo
   enddo

   !-----------------------------
   ! Compress out the land blocks
   !-----------------------------

   ii=0
   do i=1,nblocks_tot
      if(IdxT_i(i) .gt. 0) then
         ii=ii+1
         Mesh3(idxT_i(i),idxT_j(i)) = ii
      endif
   enddo
   nblocks=ii
   if (debug_blocks) then
     if (my_task == master_task) call PrintCurve(Mesh3)
   endif

   !----------------------------------------------------
   ! Compute the partitioning of the space-filling curve
   !----------------------------------------------------

   nblocksL = nblocks/nprocs
   ! every cpu gets nblocksL blocks, but the first 'extra' get nblocksL+1
   extra = mod(nblocks,nprocs)
   s1 = extra*(nblocksL+1)
   ! split curve into two curves:
   ! 1 ... s1  s2 ... nblocks
   !
   !  s1 = extra*(nblocksL+1)         (count be 0)
   !  s2 = s1+1
   !
   ! First region gets nblocksL+1 blocks per partition
   ! Second region gets nblocksL blocks per partition
!   if(debug_blocks) write(nu_diag,*) subname,'nprocs,extra,nblocks,nblocksL,s1: ', &
!                nprocs,extra,nblocks,nblocksL,s1

   !-----------------------------------------------------------
   ! Use the SFC to partition the blocks across processors
   !-----------------------------------------------------------

   do j=1,nblocks_y
   do i=1,nblocks_x
      n = (j-1)*nblocks_x + i
      ii = Mesh3(i,j)
      if(ii>0) then
        if(ii<=s1) then
           ! ------------------------------------
           ! If on the first region of curve
           ! all processes get nblocksL+1 blocks
           ! ------------------------------------
           ii=ii-1
           tmp1 = ii/(nblocksL+1)
           newDistrb%blockLocation(n) = tmp1+1
        else
           ! ------------------------------------
           ! If on the second region of curve
           ! all processes get nblocksL blocks
           ! ------------------------------------
           ii=ii-s1-1
           tmp1 = ii/nblocksL
           newDistrb%blockLocation(n) = extra + tmp1 + 1
        endif
      endif
   enddo
   enddo

   !----------------------------------------------------------------------
   !  Reset the dist data structure
   !----------------------------------------------------------------------

   globalID = 0
   newDistrb%numLocalBlocks = 0
   newDistrb%blockCnt(:) = 0

   do n=1,nblocks_tot
      globalID = n
      processor = newDistrb%blockLocation(globalID)
      if (processor > 0) then
         newDistrb%blockCnt(processor) = newDistrb%blockCnt(processor) + 1
         localID = newDistrb%blockCnt(processor)
         newDistrb%blockLocalID (globalID) = localID
      else  ! no work - eliminate block from distribution
         newDistrb%blockLocation(globalID) = 0
         newDistrb%blockLocalID (globalID) = 0
      endif
   enddo

   newDistrb%numLocalBlocks = newDistrb%blockCnt(my_task+1)

   ! set local blockGlobalID array
   allocate(newDistrb%blockGlobalID(newDistrb%numLocalBlocks), stat=istat)
   if (ice_memusage_allocErr(istat,subname//'alloc numLocalBlocks')) return
   do n = 1,nblocks_tot
      if (my_task+1 == newDistrb%blockLocation(n)) then
         localID = newDistrb%blockLocalID(n)
         newDistrb%blockGlobalID (localID) = n
      endif
   enddo

   ! set/check max_blocks
   if (max_blocks < 0) then
      max_blocks = newDistrb%numLocalBlocks
   endif

!   if (debug_blocks) then
!      if (my_task == master_task) write(nu_diag,*) subname,'dist%blockLocation:= ',dist%blockLocation
!      write(nu_diag,*) subname,'IAM: ',my_task,' SpaceCurve: Number of blocks {total,local} :=', &
!                nblocks_tot,nblocks,newDistrb%numLocalBlocks
!   endif

   !---------------------------------
   ! Deallocate temporary arrays
   !---------------------------------

   deallocate(Mesh,Mesh2,Mesh3,idxT_i,idxT_j, stat=istat)
   if (ice_memusage_allocErr(istat,subname//'dealloc meshes')) return

!----------------------------------------------------------------------

 end function create_distrb_spacecurve

!**********************************************************************

 subroutine ice_distributionRake (procWork, procID, blockWork, &
                                  priority, distribution)

!  This subroutine performs a rake algorithm to distribute the work
!  along a vector of processors.  In the rake algorithm, a work
!  threshold is first set.  Then, moving from left to right, work
!  above that threshold is raked to the next processor in line.
!  The process continues until the end of the vector is reached
!  and then the threshold is reduced by one for a second rake pass.
!  In this implementation, a priority for moving blocks is defined
!  such that the rake algorithm chooses the highest priority
!  block to be moved to the next processor.  This can be used
!  for example to always choose the eastern-most block or to
!  ensure a block does not stray too far from its neighbors.

   integer (int_kind), intent(in), dimension(:) :: &
      blockWork,          &! amount of work per block
      procID               ! global processor number

   integer (int_kind), intent(inout), dimension(:) :: &
      procWork,           &! amount of work per processor
      priority             ! priority for moving a given block

   type (distrb), intent(inout) :: &
      distribution         ! distribution to change

!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      i, n,                  &! dummy loop indices
      np1,                   &! n+1 corrected for cyclical wrap
      iproc, inext,          &! processor ids for current and next
      nprocs, numBlocks,     &! number of blocks, processors
      lastPriority,          &! priority for most recent block
      minPriority,           &! minimum priority
      lastLoc,               &! location for most recent block
      meanWork, maxWork,     &! mean,max work per processor
      diffWork,              &! work differences
      numTransfers            ! counter for number of block transfers

   character(len=*),parameter :: subname='(ice_distributionRake)'

!----------------------------------------------------------------------
!
!  initialization
!
!----------------------------------------------------------------------

   nprocs  = size(procWork)
   numBlocks = size(blockWork)

   !*** compute mean,max work per processor

   meanWork = sum(procWork)/nprocs + 1
   maxWork  = maxval(procWork)
!  residual = mod(meanWork,nprocs)

   minPriority = 1000000
   do n=1,nprocs
      iproc = procID(n)
      do i=1,numBlocks
         if (distribution%blockLocation(i) == iproc) then
            minPriority = min(minPriority,priority(i))
         endif
      end do
   end do

!----------------------------------------------------------------------
!
!  do two sets of transfers
!
!----------------------------------------------------------------------

   transferLoop: do

!----------------------------------------------------------------------
!
!     do rake across the processors
!
!----------------------------------------------------------------------

      numTransfers = 0
      do n=1,nprocs
         if (n < nprocs) then
            np1   = n+1
         else
            np1   = 1
         endif
         iproc = procID(n)
         inext = procID(np1)

         if (procWork(n) > meanWork) then !*** pass work to next

            diffWork = procWork(n) - meanWork

            rake1: do while (diffWork > 1)

               !*** attempt to find a block with the required
               !*** amount of work and with the highest priority
               !*** for moving (eg boundary blocks first)

               lastPriority = 0
               lastLoc = 0

               do i=1,numBlocks
                  if (distribution%blockLocation(i) == iproc) then
                     if (priority(i) > lastPriority ) then
                        lastPriority = priority(i)
                        lastLoc = i
                     endif
                  endif
               end do
               if (lastLoc == 0) exit rake1 ! could not shift work

               numTransfers = numTransfers + 1
               distribution%blockLocation(lastLoc) = inext
               if (np1 == 1) priority(lastLoc) = minPriority
               diffWork = diffWork - blockWork(lastLoc)

               procWork(n  ) = procWork(n  )-blockWork(lastLoc)
               procWork(np1) = procWork(np1)+blockWork(lastLoc)
            end do rake1
         endif

      end do

!----------------------------------------------------------------------
!
!     increment meanWork by one and repeat
!
!----------------------------------------------------------------------

      meanWork = meanWork + 1
      if (numTransfers == 0 .or. meanWork > maxWork) exit transferLoop

   end do transferLoop

!----------------------------------------------------------------------

end subroutine ice_distributionRake

!***********************************************************************

end module ice_distribution

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
