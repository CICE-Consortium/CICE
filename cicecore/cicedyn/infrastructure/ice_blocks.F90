!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module ice_blocks

!  This module contains data types and tools for decomposing a global
!  horizontal domain into a set of blocks.  It contains a data type
!  for describing each block and contains routines for creating and
!  querying the block decomposition for a global domain.
!
! author: Phil Jones, LANL
! Oct. 2004: Adapted from POP by William H. Lipscomb, LANL

   use ice_kinds_mod
   use ice_domain_size, only: block_size_x, block_size_y
   use ice_fileunits, only: nu_diag
   use ice_exit, only: abort_ice
   use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted

   implicit none
   private

   type, public :: block   ! block data type
      integer (int_kind) :: &
         block_id           ,&! global block number
         local_id           ,&! local address of block in current distrib
         ilo, ihi, jlo, jhi ,&! begin, end indices for physical domain
         iblock, jblock       ! cartesian i,j position for block

      logical (log_kind) :: &
         tripole,           & ! flag is true if block is at tripole bndy
         tripoleTFlag         ! tripole boundary is a T-fold

      integer (int_kind), dimension(:), pointer :: &
         i_glob, j_glob       ! global domain location for each point.
                              ! valid values between 1:nx_global, 1:ny_global.
                              ! outside that range may occur in the halo with
                              ! open or closed bcs or on the tripole.
                              ! by definition, tripole is only on the north
                              ! boundary and in that case, the j_glob values
                              ! will be valid j_glob values with minus sign.
   end type

   public :: create_blocks       ,&
             get_block           ,&
             get_block_parameter ,&
             ice_blocksGetNbrID

   integer (int_kind), parameter, public :: &
      nghost = 1       ! number of ghost cells around each block

   integer (int_kind), public :: &! size of block domain in
      nx_block, ny_block          !  x,y dir including ghost

   ! predefined directions for neighbor id routine
   ! Note: the directions that are commented out are implemented in
   !       POP but not in CICE.  If the tripole cut were in the south
   !       instead of the north, these would need to be used (and also
   !       implemented in ice_boundary.F90).
   integer (int_kind), parameter, public :: &
      ice_blocksNorth          =  1,      & ! (i  ,j+1)
      ice_blocksSouth          =  2,      & ! (i  ,j-1)
      ice_blocksEast           =  3,      & ! (i+1,j  )
      ice_blocksWest           =  4,      & ! (i-1,j  )
      ice_blocksNorthEast      =  5,      & ! (i+1,j+1)
      ice_blocksNorthWest      =  6,      & ! (i-1,j+1)
      ice_blocksSouthEast      =  7,      & ! (i+1,j-1)
      ice_blocksSouthWest      =  8         ! (i-1,j-1)
   integer (int_kind), parameter, public :: &
!      ice_blocksNorth2         =  9,      & ! (i  ,j+2)
!      ice_blocksSouth2         = 10,      & ! (i  ,j-2)
      ice_blocksEast2          = 11,      & ! (i+2,j  )
      ice_blocksWest2          = 12         ! (i-2,j  )
!      ice_blocksNorthEast2     = 13,      & ! (i+2,j+2)
!      ice_blocksNorthWest2     = 14,      & ! (i-2,j+2)
!      ice_blocksSouthEast2     = 15,      & ! (i+2,j-2)
!      ice_blocksSouthWest2     = 16         ! (i-2,j-2)
   integer (int_kind), parameter, public :: &
      ice_blocksEastNorthEast  = 17,      & ! (i+2,j+1)
!      ice_blocksEastSouthEast  = 18,      & ! (i+2,j-1)
      ice_blocksWestNorthWest  = 19         ! (i-2,j+1)
!      ice_blocksWestSouthWest  = 20,      & ! (i-2,j-1)
!      ice_blocksNorthNorthEast = 21,      & ! (i+1,j-2)
!      ice_blocksSouthSouthEast = 22,      & ! (i+1,j-2)
!      ice_blocksNorthNorthWest = 23,      & ! (i-1,j+2)
!      ice_blocksSouthSouthWest = 24         ! (i-1,j-2)

   integer (int_kind), public :: &
      nblocks_tot      ,&! total number of blocks in decomposition
      nblocks_x        ,&! tot num blocks in i direction
      nblocks_y          ! tot num blocks in j direction

   logical (kind=log_kind), public :: &
      debug_blocks       ! print verbose block information

!-----------------------------------------------------------------------
!
!  module private data
!
!-----------------------------------------------------------------------

   type (block), dimension(:), allocatable, public :: &
      all_blocks         ! block information for all blocks in domain

   integer (int_kind), dimension(:,:),allocatable, public :: &
      all_blocks_ij      ! block index stored in Cartesian order
                         !   useful for determining block index
                         !   of neighbor blocks

   integer (int_kind), dimension(:,:), allocatable, target, public :: &
      i_global,         &! global i index for each point in each block
      j_global           ! global j index for each point in each block

!***********************************************************************

contains

!***********************************************************************

 subroutine create_blocks(nx_global, ny_global, ew_boundary_type, &
                                                ns_boundary_type)

!  This subroutine decomposes the global domain into blocks and
!  fills the data structures with all the necessary block information.

   use ice_communicate, only: my_task, master_task

   integer (int_kind), intent(in) :: &
      nx_global, ny_global           ! global domain size in x,y

   character (*), intent(in) :: &
      ew_boundary_type,  &! type of boundary in logical east-west dir
      ns_boundary_type    ! type of boundary in logical north-south dir

!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      i, j, n              ,&! loop indices
      iblock, jblock       ,&! block loop indices
      is, ie, js, je         ! temp start, end indices

   character(len=*), parameter :: subname = '(create_blocks)'

!----------------------------------------------------------------------
!
!  Compute number of blocks and cartesian decomposition.
!  If the requested block size does not divide the global domain
!  size evenly, add additional block space to accomodate padding.
!
!  Compute the global indices for each block including on the halo.
!  The global indices go from 1:nx_global and 1:ny_global for
!  most of the domain including the halo that's in the internal part
!  of the domain.  On the outer boundaries, the global indices will
!  be wrapped around for the 'cyclic' option and will be given a
!  negative value on the north tripole.  Padded gridcells will be
!  given a global index of zero (0).  All other cases will extrapolate
!  the global index outside of 1:nx_global, 1:ny_global.  That means
!  the global index will go from -nghost+1:0 on the lower boundary
!  and n*_global+1:n*_global+nghost on the upper boundary and the
!  haloUpdate and scatter, for instance, will not fill those values
!  in those cases.  Other boundary condition methods will fill the 
!  outer halo values in cases where ice exists on those boundaries.
!
!----------------------------------------------------------------------

   nx_block    = block_size_x + 2*nghost ! size of block domain in x,y dir
   ny_block    = block_size_y + 2*nghost !  including ghost cells
   nblocks_x   = (nx_global-1)/block_size_x + 1
   nblocks_y   = (ny_global-1)/block_size_y + 1
   nblocks_tot = nblocks_x*nblocks_y

!----------------------------------------------------------------------
!
!  allocate block arrays
!
!----------------------------------------------------------------------

   if (.not.allocated(all_blocks)) allocate(all_blocks(nblocks_tot))
   if (.not.allocated(i_global)) allocate(i_global(nx_block,nblocks_tot))
   if (.not.allocated(j_global)) allocate(j_global(ny_block,nblocks_tot))
   if (.not.allocated(all_blocks_ij)) allocate(all_blocks_ij(nblocks_x,nblocks_y))

!----------------------------------------------------------------------
!
!  fill block data structures for all blocks in domain
!
!----------------------------------------------------------------------

   n = 0
   do jblock=1,nblocks_y
      js = (jblock-1)*block_size_y + 1
      if (js > ny_global) call abort_ice(subname// &
            ' ERROR: Bad block decomp: ny_block too large?')
      je = js + block_size_y - 1
      if (je > ny_global) je = ny_global ! pad array

      do iblock=1,nblocks_x
         n = n + 1  ! global block id

         is = (iblock-1)*block_size_x + 1
         if (is > nx_global) call abort_ice(subname// &
            ' ERROR: Bad block decomp: nx_block too large?')
         ie = is + block_size_x - 1
         if (ie > nx_global) ie = nx_global

         all_blocks(n)%block_id = n
         all_blocks(n)%iblock   = iblock
         all_blocks(n)%jblock   = jblock
         all_blocks(n)%ilo      = nghost + 1
         all_blocks(n)%jlo      = nghost + 1
         all_blocks(n)%ihi      = nx_block - nghost ! default value
         all_blocks(n)%jhi      = ny_block - nghost ! default value

         if (jblock == nblocks_y .and. &
             (ns_boundary_type == 'tripole' .or. &
             ns_boundary_type == 'tripoleT')) then
             all_blocks(n)%tripole = .true.
         else
             all_blocks(n)%tripole = .false.
         endif
         all_blocks(n)%tripoleTFlag = (ns_boundary_type == 'tripoleT')

         all_blocks_ij(iblock,jblock) = n

         do j=1,ny_block
            j_global(j,n) = js - nghost + j - 1  ! simple lower to upper counting

            !*** southern ghost cells

            if (j_global(j,n) < 1) then
               select case (ns_boundary_type)
               case ('cyclic')
                  j_global(j,n) = j_global(j,n) + ny_global
               case ('open')
                  ! lower to upper
               case ('closed')
                  ! lower to upper
               case ('tripole')
                  ! lower to upper
               case ('tripoleT')
                  ! lower to upper
               case default
                  call abort_ice(subname//' ERROR: unknown n-s bndy type')
               end select
            endif

            !*** padding required

            if (j_global(j,n) > ny_global + nghost) then
               j_global(j,n) = 0   ! padding

            !*** northern ghost cells

            else if (j_global(j,n) > ny_global) then
               select case (ns_boundary_type)
               case ('cyclic')
                  j_global(j,n) = j_global(j,n) - ny_global
               case ('open')
                  ! lower to upper
               case ('closed')
                  ! lower to upper
               case ('tripole')
                  j_global(j,n) = -j_global(j,n)  ! negative
               case ('tripoleT')
                  j_global(j,n) = -j_global(j,n)  ! negative
               case default
                  call abort_ice(subname//' ERROR: unknown n-s bndy type')
               end select

            !*** set last physical point if padded domain

            else if (j_global(j,n) == ny_global .and. &
                     j >= all_blocks(n)%jlo      .and. &
                     j <  all_blocks(n)%jhi) then
               all_blocks(n)%jhi = j   ! last physical point in padded domain
            endif
         end do

         all_blocks(n)%j_glob => j_global(:,n)

         do i=1,nx_block
            i_global(i,n) = is - nghost + i - 1  ! left to right counting

            !*** western ghost cells

            if (i_global(i,n) < 1) then
               select case (ew_boundary_type)
               case ('cyclic')
                  i_global(i,n) = i_global(i,n) + nx_global
               case ('open')
                  ! left to right
               case ('closed')
                  ! left to right
               case default
                  call abort_ice(subname//' ERROR: unknown e-w bndy type')
               end select
            endif

            !*** padded domain - fill padded region with zero

            if (i_global(i,n) > nx_global + nghost) then
               i_global(i,n) = 0

            !*** eastern ghost cells

            else if (i_global(i,n) > nx_global) then
               select case (ew_boundary_type)
               case ('cyclic')
                  i_global(i,n) = i_global(i,n) - nx_global
               case ('open')
                  ! left to right
               case ('closed')
                  ! left to right
               case default
                  call abort_ice(subname//' ERROR: unknown e-w bndy type')
               end select

            !*** last physical point in padded domain

            else if (i_global(i,n) == nx_global .and. &
                     i >= all_blocks(n)%ilo      .and. &
                     i <  all_blocks(n)%ihi) then
               all_blocks(n)%ihi = i
            endif
         end do

         all_blocks(n)%i_glob => i_global(:,n)

      end do
   end do

   if (debug_blocks) then
      if (my_task == master_task) then
      write(nu_diag,*) ' '
      write(nu_diag,'(2a)') subname,' block ID, iblock, jblock Locations:'
      do n = 1, nblocks_tot
         write(nu_diag,'(2a,3i8,l4)') subname,' global block ID, iblock, jblock, tripole:', &
         all_blocks(n)%block_id, &
         all_blocks(n)%iblock,   &
         all_blocks(n)%jblock,   &
         all_blocks(n)%tripole
      enddo
      endif
   endif

!----------------------------------------------------------------------

end subroutine create_blocks

!***********************************************************************

 function ice_blocksGetNbrID(blockID, direction, iBoundary, jBoundary) &
                             result (nbrID)

!  This function returns the block id of a neighboring block in a
!  requested direction.  Directions:
!      ice\_blocksNorth             (i  ,j+1)
!      ice\_blocksSouth             (i  ,j-1)
!      ice\_blocksEast              (i+1,j  )
!      ice\_blocksWest              (i-1,j  )
!      ice\_blocksNorthEast         (i+1,j+1)
!      ice\_blocksNorthWest         (i-1,j+1)
!      ice\_blocksSouthEast         (i  ,j-1)
!      ice\_blocksSouthWest         (i-1,j-1)
!      ice\_blocksNorth2            (i  ,j+2)
!      ice\_blocksSouth2            (i  ,j-2)
!      ice\_blocksEast2             (i+2,j  )
!      ice\_blocksWest2             (i-2,j  )
!      ice\_blocksNorthEast2        (i+2,j+2)
!      ice\_blocksNorthWest2        (i-2,j+2)
!      ice\_blocksSouthEast2        (i+2,j-2)
!      ice\_blocksSouthWest2        (i-2,j-2)
!      ice\_blocksEastNorthEast     (i+2,j+1)
!      ice\_blocksEastSouthEast     (i+2,j-1)
!      ice\_blocksWestNorthWest     (i-2,j+1)
!      ice\_blocksWestSouthWest     (i-2,j-1)
!      ice\_blocksNorthNorthEast    (i+1,j-2)
!      ice\_blocksSouthSouthEast    (i+1,j-2)
!      ice\_blocksNorthNorthWest    (i-1,j+2)
!      ice\_blocksSouthSouthWest    (i-1,j-2)
!

   integer (int_kind), intent(in)  :: &
      blockID,       &! id of block for which neighbor id requested
      direction       ! direction for which to look for neighbor -
                      !   must be one of the predefined module
                      !   variables for block direction

   character (*), intent(in) :: &
      iBoundary,     &! determines what to do at edges of domain
      jBoundary       !  options are - open, closed, cyclic, tripole, tripoleT

   integer (int_kind) :: &
      nbrID           ! block ID of neighbor in requested dir

!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      iBlock, jBlock,  &! i,j block location of current block
      inbr,   jnbr      ! i,j block location of neighboring block

   character(len=*), parameter :: subname = '(ice_blocksGetNbrID)'

!----------------------------------------------------------------------
!
!  retrieve info for current block
!
!----------------------------------------------------------------------

   call get_block_parameter(blockID, iblock=iBlock, jblock=jBlock)
   nbrID = 0   ! initial default

!----------------------------------------------------------------------
!
!  compute i,j block location of neighbor
!
!----------------------------------------------------------------------

   select case(direction)

   case (ice_blocksNorth)

      inbr = iBlock
      jnbr = jBlock + 1
      if (jnbr > nblocks_y) then
         select case(jBoundary)
         case ('open')
            jnbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = 1
         case ('tripole':'tripoleT')
            !*** return negative j value to flag tripole
            !*** i index of main northern neighbor across the
            !*** tripole cut - may also need i+1,i-1 to get
            !*** other points if there has been padding or
            !*** if the block size does not divide the domain
            !*** evenly
            inbr =  nblocks_x - iBlock + 1
            jnbr = -jBlock
         case default
            call abort_ice(subname//' ERROR: unknown north boundary')
         end select
      endif

   case (ice_blocksSouth)

      inbr = iBlock
      jnbr = jBlock - 1
      if (jnbr < 1) then
         select case(jBoundary)
         case ('open')
            jnbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = nblocks_y
         case ('tripole')
            jnbr = 0 ! do not write into the neighbor's ghost cells
         case ('tripoleT')
            jnbr = 0 ! do not write into the neighbor's ghost cells
         case default
            call abort_ice(subname//' ERROR: unknown south boundary')
         end select
      endif

   case (ice_blocksEast )

      inbr = iBlock + 1
      jnbr = jBlock
      if (inbr > nblocks_x) then
         select case(iBoundary)
         case ('open')
            inbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = 1
         case default
            call abort_ice(subname//' ERROR: unknown east boundary')
         end select
      endif

   case (ice_blocksWest )

      jnbr = jBlock
      inbr = iBlock - 1
      if (inbr < 1) then
         select case(iBoundary)
         case ('open')
            inbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = nblocks_x
         case default
            call abort_ice(subname//' ERROR: unknown west boundary')
         end select
      endif

   case (ice_blocksNorthEast)

      inbr = iBlock + 1
      jnbr = jBlock + 1
      if (inbr > nblocks_x) then
         select case(iBoundary)
         case ('open')
            inbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = 1
         case default
            call abort_ice(subname//' ERROR: unknown east boundary')
         end select
      endif
      if (jnbr > nblocks_y) then
         select case(jBoundary)
         case ('open')
            jnbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = 1
         case ('tripole':'tripoleT')
            !*** return negative j value to flag tripole
            !*** i index of main northern neighbor across the
            !*** tripole cut - may also need i+1,i-1 to get
            !*** other points if there has been padding or
            !*** if the block size does not divide the domain
            !*** evenly
            inbr =  nblocks_x - iBlock
            if (inbr == 0) inbr = nblocks_x
            jnbr = -jBlock
         case default
            call abort_ice(subname//' ERROR: unknown north boundary')
         end select
      endif

   case (ice_blocksNorthWest)

      inbr = iBlock - 1
      jnbr = jBlock + 1
      if (inbr < 1) then
         select case(iBoundary)
         case ('open')
            inbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = nblocks_x
         case default
            call abort_ice(subname//' ERROR: unknown west boundary')
         end select
      endif
      if (jnbr > nblocks_y) then
         select case(jBoundary)
         case ('open')
            jnbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = 1
         case ('tripole':'tripoleT')
            !*** return negative j value to flag tripole
            !*** i index of main northern neighbor across the
            !*** tripole cut - may also need i+1,i-1 to get
            !*** other points if there has been padding or
            !*** if the block size does not divide the domain
            !*** evenly
            inbr =  nblocks_x - iBlock + 2
            if (inbr > nblocks_x) inbr = 1
            jnbr = -jBlock
         case default
            call abort_ice(subname//' ERROR: unknown north boundary')
         end select
      endif

   case (ice_blocksSouthEast )

      inbr = iBlock + 1
      jnbr = jBlock - 1
      if (inbr > nblocks_x) then
         select case(iBoundary)
         case ('open')
            inbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = 1
         case default
            call abort_ice(subname//' ERROR: unknown east boundary')
         end select
      endif
      if (jnbr < 1) then
         select case(jBoundary)
         case ('open')
            inbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = nblocks_y
         case ('tripole')
            jnbr = 0 ! do not write into the neighbor's ghost cells
         case ('tripoleT')
            jnbr = 0 ! do not write into the neighbor's ghost cells
         case default
            call abort_ice(subname//' ERROR: unknown south boundary')
         end select
      endif

   case (ice_blocksSouthWest )
      inbr = iBlock - 1
      jnbr = jBlock - 1
      if (inbr < 1) then
         select case(iBoundary)
         case ('open')
            inbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = nblocks_x
         case default
            call abort_ice(subname//' ERROR: unknown west boundary')
         end select
      endif
      if (jnbr < 1) then
         select case(jBoundary)
         case ('open')
            jnbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = nblocks_y
         case ('tripole')
            jnbr = 0 ! do not write into the neighbor's ghost cells
         case ('tripoleT')
            jnbr = 0 ! do not write into the neighbor's ghost cells
         case default
            call abort_ice(subname//' ERROR: unknown south boundary')
         end select
      endif

   case (ice_blocksEast2)

      inbr = iBlock + 2
      jnbr = jBlock
      if (inbr > nblocks_x) then
         select case(iBoundary)
         case ('open')
            inbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = inbr - nblocks_x
         case default
            call abort_ice(subname//' ERROR: unknown east boundary')
         end select
      endif

   case (ice_blocksWest2)
      jnbr = jBlock
      inbr = iBlock - 2
      if (inbr < 1) then
         select case(iBoundary)
         case ('open')
            inbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = nblocks_x + inbr
         case default
            call abort_ice(subname//' ERROR: unknown west boundary')
         end select
      endif

   case (ice_blocksEastNorthEast)

      inbr = iBlock + 2
      jnbr = jBlock + 1
      if (inbr > nblocks_x) then
         select case(iBoundary)
         case ('open')
            inbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = inbr - nblocks_x
         case default
            call abort_ice(subname//' ERROR: unknown east boundary')
         end select
      endif
      if (jnbr > nblocks_y) then
         select case(jBoundary)
         case ('open')
            jnbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = jnbr - nblocks_y
         case ('tripole':'tripoleT')
            !*** return negative j value to flag tripole
            !*** i index of main northern neighbor across the
            !*** tripole cut - may also need i+1,i-1 to get
            !*** other points if there has been padding or
            !*** if the block size does not divide the domain
            !*** evenly
            inbr =  nblocks_x - iBlock - 1
            if (inbr <= 0) inbr = inbr + nblocks_x
            jnbr = -jBlock
         case default
            call abort_ice(subname//' ERROR: unknown north boundary')
         end select
      endif

   case (ice_blocksWestNorthWest)

      inbr = iBlock - 2
      jnbr = jBlock + 1
      if (inbr < 1) then
         select case(iBoundary)
         case ('open')
            inbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            inbr = 0
         case ('cyclic')
            inbr = nblocks_x + inbr
         case default
            call abort_ice(subname//' ERROR: unknown west boundary')
         end select
      endif
      if (jnbr > nblocks_y) then
         select case(jBoundary)
         case ('open')
            jnbr = 0 ! do not write into the neighbor's ghost cells
         case ('closed')
            jnbr = 0
         case ('cyclic')
            jnbr = jnbr + nblocks_y
         case ('tripole':'tripoleT')
            !*** return negative j value to flag tripole
            !*** i index of main northern neighbor across the
            !*** tripole cut - may also need i+1,i-1 to get
            !*** other points if there has been padding or
            !*** if the block size does not divide the domain
            !*** evenly
            inbr =  nblocks_x - iBlock + 3
            if (inbr > nblocks_x) inbr = inbr - nblocks_x
            jnbr = -jBlock
         case default
            call abort_ice(subname//' ERROR: unknown north boundary')
         end select
      endif

   case default

      call abort_ice(subname//' ERROR: unknown direction')
      return

   end select

!----------------------------------------------------------------------
!
!  now get block id for this neighbor block
!
!----------------------------------------------------------------------

   if (inbr > 0 .and. jnbr > 0) then
      nbrID = all_blocks_ij(inbr,jnbr)
   else if (inbr > 0 .and. jnbr < 0) then  ! tripole upper boundary
      !*** return negative value to flag tripole
      nbrID = -all_blocks_ij(inbr,abs(jnbr))
   else
      nbrID = 0   ! neighbor outside domain
   endif

!----------------------------------------------------------------------

 end function ice_blocksGetNbrID

!**********************************************************************

 function get_block(block_id,local_id)

!  This function returns the block data structure for the block
!  associated with the input block id.

   integer (int_kind), intent(in) :: &
      block_id,   &! global block id for requested block info
      local_id     ! local  block id to assign to this block

   type (block) :: &
      get_block    ! block information returned for requested block

   character(len=*), parameter :: subname = '(get_block)'

!----------------------------------------------------------------------
!
!  check for valid id.  if valid, return block info for requested block
!
!----------------------------------------------------------------------

   if (block_id < 1 .or. block_id > nblocks_tot) then
      call abort_ice(subname//' ERROR: invalid block_id')
   endif

   get_block = all_blocks(block_id)
   get_block%local_id = local_id

!----------------------------------------------------------------------

 end function get_block

!**********************************************************************

 subroutine get_block_parameter(block_id, local_id,           &
                                ilo, ihi, jlo, jhi,           &
                                iblock, jblock, tripole,      &
                                i_glob, j_glob)

!  This routine returns requested parts of the block data type
!  for the block associated with the input block id

   integer (int_kind), intent(in) :: &
      block_id   ! global block id for which parameters are requested

   !(optional) parts of block data type to extract if requested

   integer (int_kind), intent(out), optional :: &
      local_id           ,&! local id assigned to block in current distrb
      ilo, ihi, jlo, jhi ,&! begin,end indices for physical domain
      iblock, jblock       ! cartesian i,j position for bloc

   logical (log_kind), intent(out), optional :: &
      tripole              ! flag is true if block on tripole bndy

   integer (int_kind), dimension(:), pointer, optional :: &
      i_glob, j_glob     ! global domain location for each point

   character(len=*), parameter :: subname = '(get_block_parameter)'

!----------------------------------------------------------------------
!
!  extract each component of data type if requested
!
!----------------------------------------------------------------------

   if (block_id < 1 .or. block_id > nblocks_tot) then
      call abort_ice(subname//' ERROR: invalid block_id')
   endif

   if (present(local_id)) local_id = all_blocks(block_id)%local_id
   if (present(ilo     )) ilo      = all_blocks(block_id)%ilo
   if (present(ihi     )) ihi      = all_blocks(block_id)%ihi
   if (present(jlo     )) jlo      = all_blocks(block_id)%jlo
   if (present(jhi     )) jhi      = all_blocks(block_id)%jhi
   if (present(iblock  )) iblock   = all_blocks(block_id)%iblock
   if (present(jblock  )) jblock   = all_blocks(block_id)%jblock
   if (present(i_glob  )) i_glob   => all_blocks(block_id)%i_glob
   if (present(j_glob  )) j_glob   => all_blocks(block_id)%j_glob
   if (present(tripole )) tripole  = all_blocks(block_id)%tripole

!----------------------------------------------------------------------

 end subroutine get_block_parameter

!**********************************************************************

 end module ice_blocks

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
