!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module ice_gather_scatter

!  This module contains routines that mimic the behavior of the mpi
!  version in the case of a single processor:  gathering data to a single
!  processor from a distributed array, and scattering data from a
!  single processor to a distributed array.
!
! author: Phil Jones, LANL
! Oct. 2004: Adapted from POP version by William H. Lipscomb, LANL
! Jan. 2008: Elizabeth Hunke replaced old routines with new POP
!              infrastructure, added specialized routine scatter_global_stress

   use ice_kinds_mod
   use ice_constants
   use ice_blocks, only: block, get_block, nghost, nx_block, ny_block, &
       nblocks_x, nblocks_y, nblocks_tot
   use ice_distribution, only: distrb
   use ice_domain_size, only: nx_global, ny_global
   use ice_exit, only: abort_ice
   use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted

   implicit none
   private

   public :: gather_global,      &
             gather_global_ext,  &
             scatter_global,     &
             scatter_global_ext, &
             scatter_global_stress

!-----------------------------------------------------------------------
!
!  overload module functions
!
!-----------------------------------------------------------------------

   interface gather_global_ext
     module procedure gather_global_ext_dbl,  &
                      gather_global_ext_int,  &
                      gather_global_ext_log
   end interface

   interface gather_global
     module procedure gather_global_dbl,  &
                      gather_global_real, &
                      gather_global_int
   end interface

   interface scatter_global
     module procedure scatter_global_dbl,  &
                      scatter_global_real, &
                      scatter_global_int
   end interface

!-----------------------------------------------------------------------
!
!  module variables
!
!-----------------------------------------------------------------------

!***********************************************************************

 contains

!***********************************************************************

 subroutine gather_global_dbl(ARRAY_G, ARRAY, dst_task, src_dist, spc_val)

!  This subroutine gathers a distributed array to a global-sized
!  array on the processor dst_task.
!
!  This is the specific inteface for double precision arrays
!  corresponding to the generic interface gather_global.  It is shown
!  to provide information on the generic interface (the generic
!  interface is identical, but chooses a specific inteface based
!  on the data type of the input argument).

   integer (int_kind), intent(in) :: &
     dst_task   ! task to which array should be gathered

   type (distrb), intent(in) :: &
     src_dist   ! distribution of blocks in the source array

   real (dbl_kind), dimension(:,:,:), intent(in) :: &
     ARRAY      ! array containing horizontal slab of distributed field

   real (dbl_kind), intent(in), optional :: &
     spc_val

   real (dbl_kind), dimension(:,:), intent(inout) :: &
     ARRAY_G    ! array containing global horizontal field on dst_task

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n       ! dummy loop counters

   real (dbl_kind) :: &
     special_value

   type (block) :: &
     this_block  ! block info for current block

   character(len=*), parameter :: subname = '(gather_global_dbl)'

   if (present(spc_val)) then
      special_value = spc_val
   else
      special_value = spval_dbl
   endif

!-----------------------------------------------------------------------
!
!  copy local array into block decomposition
!
!-----------------------------------------------------------------------

   do n=1,nblocks_tot

      this_block = get_block(n,n)

      !*** copy local blocks

      if (src_dist%blockLocation(n) /= 0) then

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
            ARRAY_G(this_block%i_glob(i), &
                    this_block%j_glob(j)) = &
                  ARRAY(i,j,src_dist%blockLocalID(n))
         end do
         end do

      else !*** fill land blocks with special values

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
            ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = special_value
         end do
         end do

      endif

   end do

!-----------------------------------------------------------------------

 end subroutine gather_global_dbl

!***********************************************************************

 subroutine gather_global_real(ARRAY_G, ARRAY, dst_task, src_dist)

!-----------------------------------------------------------------------
!
!  This subroutine gathers a distributed array to a global-sized
!  array on the processor dst_task.
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!  input variables
!
!-----------------------------------------------------------------------

   integer (int_kind), intent(in) :: &
     dst_task       ! task to which array should be gathered

   type (distrb), intent(in) :: &
     src_dist       ! distribution of blocks in the source array

   real (real_kind), dimension(:,:,:), intent(in) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  output variables
!
!-----------------------------------------------------------------------

   real (real_kind), dimension(:,:), intent(inout) :: &
     ARRAY_G        ! array containing global field on dst_task

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,n        ! dummy loop counters

   type (block) :: &
      this_block   ! block info for current block

   character(len=*), parameter :: subname = '(gather_global_real)'

!-----------------------------------------------------------------------
!
!  copy local array into block decomposition
!
!-----------------------------------------------------------------------

   do n=1,nblocks_tot

      this_block = get_block(n,n)

      !*** copy local blocks

      if (src_dist%blockLocation(n) /= 0) then

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
            ARRAY_G(this_block%i_glob(i), &
                    this_block%j_glob(j)) = &
                  ARRAY(i,j,src_dist%blockLocalID(n))
         end do
         end do

      else !*** fill land blocks with zeroes

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
            ARRAY_G(this_block%i_glob(i), &
                    this_block%j_glob(j)) = 0._real_kind
         end do
         end do

      endif

   end do

!-----------------------------------------------------------------------

 end subroutine gather_global_real

!***********************************************************************

 subroutine gather_global_int(ARRAY_G, ARRAY, dst_task, src_dist)

!-----------------------------------------------------------------------
!
!  This subroutine gathers a distributed array to a global-sized
!  array on the processor dst_task.
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!  input variables
!
!-----------------------------------------------------------------------

   integer (int_kind), intent(in) :: &
     dst_task       ! task to which array should be gathered

   type (distrb), intent(in) :: &
     src_dist       ! distribution of blocks in the source array

   integer (int_kind), dimension(:,:,:), intent(in) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  output variables
!
!-----------------------------------------------------------------------

   integer (int_kind), dimension(:,:), intent(inout) :: &
     ARRAY_G        ! array containing global field on dst_task

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,n         ! dummy loop counters

   type (block) :: &
      this_block    ! block info for current block

   character(len=*), parameter :: subname = '(gather_global_int)'

!-----------------------------------------------------------------------
!
!  copy local array into block decomposition
!
!-----------------------------------------------------------------------

   do n=1,nblocks_tot

      this_block = get_block(n,n)

      !*** copy local blocks

      if (src_dist%blockLocation(n) /= 0) then

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
            ARRAY_G(this_block%i_glob(i), &
                    this_block%j_glob(j)) = &
                  ARRAY(i,j,src_dist%blockLocalID(n))
         end do
         end do

      else !*** fill land blocks with zeroes

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
            ARRAY_G(this_block%i_glob(i), &
                    this_block%j_glob(j)) = 0
         end do
         end do

      endif

   end do

!-----------------------------------------------------------------------

 end subroutine gather_global_int

!***********************************************************************

 subroutine gather_global_ext_dbl(ARRAY_G, ARRAY, dst_task, src_dist, spc_val)

!  This subroutine gathers a distributed array to a global-sized
!  array on the processor dst_task, including ghost cells.

   integer (int_kind), intent(in) :: &
     dst_task   ! task to which array should be gathered

   type (distrb), intent(in) :: &
     src_dist   ! distribution of blocks in the source array

   real (dbl_kind), dimension(:,:,:), intent(in) :: &
     ARRAY      ! array containing horizontal slab of distributed field

   real (dbl_kind), optional :: &
     spc_val

   real (dbl_kind), dimension(:,:), intent(inout) :: &
     ARRAY_G    ! array containing global horizontal field on dst_task

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n          ,&! dummy loop counters
     nx, ny           ! global dimensions

   real (dbl_kind) :: &
     special_value

   type (block) :: &
     this_block  ! block info for current block

   character(len=*), parameter :: subname = '(gather_global_ext_dbl)'

   if (present(spc_val)) then
      special_value = spc_val
   else
      special_value = spval_dbl
   endif
   ARRAY_G = special_value

   nx = nx_global + 2*nghost
   ny = ny_global + 2*nghost

!-----------------------------------------------------------------------
!
!  copy local array into block decomposition
!
!-----------------------------------------------------------------------

   do n=1,nblocks_tot

      this_block = get_block(n,n)

      !*** copy local blocks

      if (src_dist%blockLocation(n) /= 0) then

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
            ARRAY_G(this_block%i_glob(i)+nghost, &
                    this_block%j_glob(j)+nghost) = &
                  ARRAY(i,j,src_dist%blockLocalID(n))
         end do
         end do

         ! fill ghost cells
         if (this_block%jblock == 1) then
            ! south block
            do j=1, nghost
            do i=this_block%ilo,this_block%ihi
              ARRAY_G(this_block%i_glob(i)+nghost,j) = &
              ARRAY  (i,j,src_dist%blockLocalID(n))
            end do
            end do
            if (this_block%iblock == 1) then
               ! southwest corner
               do j=1, nghost
               do i=1, nghost
                 ARRAY_G(i,j) = &
                 ARRAY  (i,j,src_dist%blockLocalID(n))
               end do
               end do
            endif
         endif
         if (this_block%jblock == nblocks_y) then
            ! north block
            do j=1, nghost
            do i=this_block%ilo,this_block%ihi
              ARRAY_G(this_block%i_glob(i)+nghost, &
                      ny_global + nghost + j) = &
              ARRAY  (i,this_block%jhi+nghost-j+1,src_dist%blockLocalID(n))
            end do
            end do
            if (this_block%iblock == nblocks_x) then
               ! northeast corner
               do j=1, nghost
               do i=1, nghost
                 ARRAY_G(nx-i+1, ny-j+1) = &
                 ARRAY  (this_block%ihi+nghost-i+1, &
                         this_block%jhi+nghost-j+1, &
                         src_dist%blockLocalID(n))
               end do
               end do
            endif
         endif
         if (this_block%iblock == 1) then
            ! west block
            do j=this_block%jlo,this_block%jhi
            do i=1, nghost
              ARRAY_G(i,this_block%j_glob(j)+nghost) = &
              ARRAY  (i,j,src_dist%blockLocalID(n))
            end do
            end do
            if (this_block%jblock == nblocks_y) then
               ! northwest corner
               do j=1, nghost
               do i=1, nghost
                 ARRAY_G(i,                   ny-j+1) = &
                 ARRAY  (i,this_block%jhi+nghost-j+1,src_dist%blockLocalID(n))
               end do
               end do
            endif
         endif
         if (this_block%iblock == nblocks_x) then
            ! east block
            do j=this_block%jlo,this_block%jhi
            do i=1, nghost
              ARRAY_G(nx_global + nghost + i, &
                      this_block%j_glob(j)+nghost) = &
              ARRAY  (this_block%ihi+nghost-i+1,j,src_dist%blockLocalID(n))
            end do
            end do
            if (this_block%jblock == 1) then
               ! southeast corner
               do j=1, nghost
               do i=1, nghost
                 ARRAY_G(                   nx-i+1,j) = &
                 ARRAY  (this_block%ihi+nghost-i+1,j,src_dist%blockLocalID(n))
               end do
               end do
            endif
         endif

      endif  ! src_dist%blockLocation

   end do

!-----------------------------------------------------------------------

 end subroutine gather_global_ext_dbl

!***********************************************************************

 subroutine gather_global_ext_int(ARRAY_G, ARRAY, dst_task, src_dist, spc_val)

!  This subroutine gathers a distributed array to a global-sized
!  array on the processor dst_task, including ghost cells.

   integer (int_kind), intent(in) :: &
     dst_task   ! task to which array should be gathered

   type (distrb), intent(in) :: &
     src_dist   ! distribution of blocks in the source array

   integer (int_kind), dimension(:,:,:), intent(in) :: &
     ARRAY      ! array containing horizontal slab of distributed field

   integer (int_kind), optional :: &
     spc_val

   integer (int_kind), dimension(:,:), intent(inout) :: &
     ARRAY_G    ! array containing global horizontal field on dst_task

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n          ,&! dummy loop counters
     nx, ny           ! global dimensions

   integer (int_kind) :: &
     special_value

   type (block) :: &
     this_block  ! block info for current block

   character(len=*), parameter :: subname = '(gather_global_ext_int)'

   if (present(spc_val)) then
      special_value = spc_val
   else
      special_value = -9999
   endif
   ARRAY_G = special_value

   nx = nx_global + 2*nghost
   ny = ny_global + 2*nghost

!-----------------------------------------------------------------------
!
!  copy local array into block decomposition
!
!-----------------------------------------------------------------------

   do n=1,nblocks_tot

      this_block = get_block(n,n)

      !*** copy local blocks

      if (src_dist%blockLocation(n) /= 0) then

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
            ARRAY_G(this_block%i_glob(i)+nghost, &
                    this_block%j_glob(j)+nghost) = &
                  ARRAY(i,j,src_dist%blockLocalID(n))
         end do
         end do

         ! fill ghost cells
         if (this_block%jblock == 1) then
            ! south block
            do j=1, nghost
            do i=this_block%ilo,this_block%ihi
              ARRAY_G(this_block%i_glob(i)+nghost,j) = &
              ARRAY  (i,j,src_dist%blockLocalID(n))
            end do
            end do
            if (this_block%iblock == 1) then
               ! southwest corner
               do j=1, nghost
               do i=1, nghost
                 ARRAY_G(i,j) = &
                 ARRAY  (i,j,src_dist%blockLocalID(n))
               end do
               end do
            endif
         endif
         if (this_block%jblock == nblocks_y) then
            ! north block
            do j=1, nghost
            do i=this_block%ilo,this_block%ihi
              ARRAY_G(this_block%i_glob(i)+nghost, &
                      ny_global + nghost + j) = &
              ARRAY  (i,this_block%jhi+nghost-j+1,src_dist%blockLocalID(n))
            end do
            end do
            if (this_block%iblock == nblocks_x) then
               ! northeast corner
               do j=1, nghost
               do i=1, nghost
                 ARRAY_G(nx-i+1, ny-j+1) = &
                 ARRAY  (this_block%ihi+nghost-i+1, &
                         this_block%jhi+nghost-j+1, &
                         src_dist%blockLocalID(n))
               end do
               end do
            endif
         endif
         if (this_block%iblock == 1) then
            ! west block
            do j=this_block%jlo,this_block%jhi
            do i=1, nghost
              ARRAY_G(i,this_block%j_glob(j)+nghost) = &
              ARRAY  (i,j,src_dist%blockLocalID(n))
            end do
            end do
            if (this_block%jblock == nblocks_y) then
               ! northwest corner
               do j=1, nghost
               do i=1, nghost
                 ARRAY_G(i,                   ny-j+1) = &
                 ARRAY  (i,this_block%jhi+nghost-j+1,src_dist%blockLocalID(n))
               end do
               end do
            endif
         endif
         if (this_block%iblock == nblocks_x) then
            ! east block
            do j=this_block%jlo,this_block%jhi
            do i=1, nghost
              ARRAY_G(nx_global + nghost + i, &
                      this_block%j_glob(j)+nghost) = &
              ARRAY  (this_block%ihi+nghost-i+1,j,src_dist%blockLocalID(n))
            end do
            end do
            if (this_block%jblock == 1) then
               ! southeast corner
               do j=1, nghost
               do i=1, nghost
                 ARRAY_G(                   nx-i+1,j) = &
                 ARRAY  (this_block%ihi+nghost-i+1,j,src_dist%blockLocalID(n))
               end do
               end do
            endif
         endif

      endif  ! src_dist%blockLocation

   end do

!-----------------------------------------------------------------------

 end subroutine gather_global_ext_int

!***********************************************************************

 subroutine gather_global_ext_log(ARRAY_G, ARRAY, dst_task, src_dist, spc_val)

!  This subroutine gathers a distributed array to a global-sized
!  array on the processor dst_task, including ghost cells.

   integer (int_kind), intent(in) :: &
     dst_task   ! task to which array should be gathered

   type (distrb), intent(in) :: &
     src_dist   ! distribution of blocks in the source array

   logical (log_kind), dimension(:,:,:), intent(in) :: &
     ARRAY      ! array containing horizontal slab of distributed field

   logical (log_kind), optional :: &
     spc_val

   logical (log_kind), dimension(:,:), intent(inout) :: &
     ARRAY_G    ! array containing global horizontal field on dst_task

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n          ,&! dummy loop counters
     nx, ny           ! global dimensions

   logical (log_kind) :: &
     special_value

   type (block) :: &
     this_block  ! block info for current block

   character(len=*), parameter :: subname = '(gather_global_ext_log)'

   if (present(spc_val)) then
      special_value = spc_val
   else
      special_value = .false.
   endif
   ARRAY_G = special_value

   nx = nx_global + 2*nghost
   ny = ny_global + 2*nghost

!-----------------------------------------------------------------------
!
!  copy local array into block decomposition
!
!-----------------------------------------------------------------------

   do n=1,nblocks_tot

      this_block = get_block(n,n)

      !*** copy local blocks

      if (src_dist%blockLocation(n) /= 0) then

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
            ARRAY_G(this_block%i_glob(i)+nghost, &
                    this_block%j_glob(j)+nghost) = &
                  ARRAY(i,j,src_dist%blockLocalID(n))
         end do
         end do

         ! fill ghost cells
         if (this_block%jblock == 1) then
            ! south block
            do j=1, nghost
            do i=this_block%ilo,this_block%ihi
              ARRAY_G(this_block%i_glob(i)+nghost,j) = &
              ARRAY  (i,j,src_dist%blockLocalID(n))
            end do
            end do
            if (this_block%iblock == 1) then
               ! southwest corner
               do j=1, nghost
               do i=1, nghost
                 ARRAY_G(i,j) = &
                 ARRAY  (i,j,src_dist%blockLocalID(n))
               end do
               end do
            endif
         endif
         if (this_block%jblock == nblocks_y) then
            ! north block
            do j=1, nghost
            do i=this_block%ilo,this_block%ihi
              ARRAY_G(this_block%i_glob(i)+nghost, &
                      ny_global + nghost + j) = &
              ARRAY  (i,this_block%jhi+nghost-j+1,src_dist%blockLocalID(n))
            end do
            end do
            if (this_block%iblock == nblocks_x) then
               ! northeast corner
               do j=1, nghost
               do i=1, nghost
                 ARRAY_G(nx-i+1, ny-j+1) = &
                 ARRAY  (this_block%ihi+nghost-i+1, &
                         this_block%jhi+nghost-j+1, &
                         src_dist%blockLocalID(n))
               end do
               end do
            endif
         endif
         if (this_block%iblock == 1) then
            ! west block
            do j=this_block%jlo,this_block%jhi
            do i=1, nghost
              ARRAY_G(i,this_block%j_glob(j)+nghost) = &
              ARRAY  (i,j,src_dist%blockLocalID(n))
            end do
            end do
            if (this_block%jblock == nblocks_y) then
               ! northwest corner
               do j=1, nghost
               do i=1, nghost
                 ARRAY_G(i,                   ny-j+1) = &
                 ARRAY  (i,this_block%jhi+nghost-j+1,src_dist%blockLocalID(n))
               end do
               end do
            endif
         endif
         if (this_block%iblock == nblocks_x) then
            ! east block
            do j=this_block%jlo,this_block%jhi
            do i=1, nghost
              ARRAY_G(nx_global + nghost + i, &
                      this_block%j_glob(j)+nghost) = &
              ARRAY  (this_block%ihi+nghost-i+1,j,src_dist%blockLocalID(n))
            end do
            end do
            if (this_block%jblock == 1) then
               ! southeast corner
               do j=1, nghost
               do i=1, nghost
                 ARRAY_G(                   nx-i+1,j) = &
                 ARRAY  (this_block%ihi+nghost-i+1,j,src_dist%blockLocalID(n))
               end do
               end do
            endif
         endif

      endif  ! src_dist%blockLocation

   end do

!-----------------------------------------------------------------------

 end subroutine gather_global_ext_log

!***********************************************************************

 subroutine scatter_global_dbl(ARRAY, ARRAY_G, src_task, dst_dist, &
                               field_loc, field_type)

!  This subroutine scatters a global-sized array on the processor
!  src\_task to a distribution of blocks given by dst\_dist.
!
!  This is the specific interface for double precision arrays
!  corresponding to the generic interface scatter_global.  It is shown
!  to provide information on the generic interface (the generic
!  interface is identical, but chooses a specific interface based
!  on the data type of the input argument).

   integer (int_kind), intent(in) :: &
     src_task       ! task from which array should be scattered

   type (distrb), intent(in) :: &
     dst_dist       ! distribution of resulting blocks

   real (dbl_kind), dimension(:,:), intent(in) :: &
     ARRAY_G        ! array containing global field on src_task

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      field_loc                  ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

   real (dbl_kind), dimension(:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n,              &! dummy loop indices
     isrc, jsrc,         &! source addresses
     xoffset, yoffset,   &! offsets for tripole boundary conditions
     yoffset2,           &!
     isign,              &! sign factor for tripole boundary conditions
     dst_block            ! local block index in dest distribution

   type (block) :: &
     this_block  ! block info for current block

   character(len=*), parameter :: subname = '(scatter_global_dbl)'

!-----------------------------------------------------------------------
!
!  initialize return array to zero and set up tripole quantities
!
!-----------------------------------------------------------------------

   ARRAY = c0

   this_block = get_block(1,1) ! for the tripoleTflag - all blocks have it
   if (this_block%tripoleTFlag) then
     select case (field_loc)
     case (field_loc_center)   ! cell center location
        xoffset = 2
        yoffset = 0
     case (field_loc_NEcorner) ! cell corner (velocity) location
        xoffset = 1
        yoffset = -1
     case (field_loc_Eface)    ! cell face location
        xoffset = 1
        yoffset = 0
     case (field_loc_Nface)    ! cell face location
        xoffset = 2
        yoffset = -1
     case (field_loc_noupdate) ! ghost cells never used - use cell center
        xoffset = 1
        yoffset = 1
     end select
   else
     select case (field_loc)
     case (field_loc_center)   ! cell center location
        xoffset = 1
        yoffset = 1
     case (field_loc_NEcorner) ! cell corner (velocity) location
        xoffset = 0
        yoffset = 0
     case (field_loc_Eface)    ! cell face location
        xoffset = 0
        yoffset = 1
     case (field_loc_Nface)    ! cell face location
        xoffset = 1
        yoffset = 0
     case (field_loc_noupdate) ! ghost cells never used - use cell center
        xoffset = 1
        yoffset = 1
     end select
   endif

   select case (field_type)
   case (field_type_scalar)
      isign =  1
   case (field_type_vector)
      isign = -1
   case (field_type_angle)
      isign = -1
   case (field_type_noupdate) ! ghost cells not needed - use cell center
      isign =  1
   case default
      call abort_ice(subname//'ERROR: Unknown field type in scatter')
   end select

!-----------------------------------------------------------------------
!
!  copy blocks of global array into local block distribution
!
!-----------------------------------------------------------------------

   do n=1,nblocks_tot

      if (dst_dist%blockLocation(n) /= 0) then

         this_block = get_block(n,n)
         dst_block  = dst_dist%blockLocalID(n)

         do j=1,ny_block
            if (this_block%jblock == nblocks_y .and. this_block%j_glob(j) < 0) then
               ! tripole is top block with j_glob < 0
               ! for yoffset=0 or 1, yoffset2=0,0
               ! for yoffset=-1, yoffset2=0,1, for u-rows on T-fold grid
               do yoffset2=0,max(yoffset,0)-yoffset
                 jsrc = ny_global + yoffset + yoffset2 + &
                      (this_block%j_glob(j) + ny_global)
                 do i=1,nx_block
                   if (this_block%i_glob(i) /= 0) then
                      isrc = nx_global + xoffset - this_block%i_glob(i)
                      if (isrc < 1) isrc = isrc + nx_global
                      if (isrc > nx_global) isrc = isrc - nx_global
                      ARRAY(i,j-yoffset2,dst_block) &
                        = isign * ARRAY_G(isrc,jsrc)
                   endif
                 end do
               end do
            else
               ! normal block
               do i=1,nx_block
                  isrc = this_block%i_glob(i)
                  jsrc = this_block%j_glob(j)
                  if (isrc >=1 .and. isrc <= nx_global .and. &
                      jsrc >=1 .and. jsrc <= ny_global) &
                     ARRAY(i,j,dst_block) = ARRAY_G(isrc,jsrc)
               end do
            endif
         end do
      endif ! dst block not land
   end do  ! block loop

   !-----------------------------------------------------------------
   ! Set ghost cell values to 0 for noupdate
   !-----------------------------------------------------------------

   if (field_loc == field_loc_noupdate) then
     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) /= 0 .and. &
           dst_dist%blockLocalID(n) > 0) then

         dst_block = dst_dist%blockLocalID(n)
         this_block = get_block(n,n)
         ! north edge
         do j = this_block%jhi+1,ny_block
         do i = 1, nx_block
            ARRAY (i,j,dst_block) = c0
         enddo
         enddo
         ! east edge
         do j = 1, ny_block
         do i = this_block%ihi+1,nx_block
            ARRAY (i,j,dst_block) = c0
         enddo
         enddo
         ! south edge
         do j = 1, this_block%jlo-1
         do i = 1, nx_block
            ARRAY (i,j,dst_block) = c0
         enddo
         enddo
         ! west edge
         do j = 1, ny_block
         do i = 1, this_block%ilo-1
            ARRAY (i,j,dst_block) = c0
         enddo
         enddo

       endif
     enddo
   endif

!-----------------------------------------------------------------------

 end subroutine scatter_global_dbl

!***********************************************************************

 subroutine scatter_global_real(ARRAY, ARRAY_G, src_task, dst_dist, &
                                field_loc, field_type)

!-----------------------------------------------------------------------
!
!  This subroutine scatters a global-sized array on the processor
!  src\_task to a distribution of blocks given by dst\_dist.
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!  input variables
!
!-----------------------------------------------------------------------

   integer (int_kind), intent(in) :: &
     src_task       ! task from which array should be scattered

   type (distrb), intent(in) :: &
     dst_dist       ! distribution of resulting blocks

   real (real_kind), dimension(:,:), intent(in) :: &
     ARRAY_G        ! array containing global field on src_task

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      field_loc                  ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

!-----------------------------------------------------------------------
!
!  output variables
!
!-----------------------------------------------------------------------

   real (real_kind), dimension(:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n,              &! dummy loop indices
     isrc, jsrc,         &! source addresses
     xoffset, yoffset,   &! offsets for tripole boundary conditions
     yoffset2,           &!
     isign,              &! sign factor for tripole boundary conditions
     dst_block            ! local block index in dest distribution

   type (block) :: &
     this_block  ! block info for current block

   character(len=*), parameter :: subname = '(scatter_global_real)'

!-----------------------------------------------------------------------
!
!  initialize return array to zero and set up tripole quantities
!
!-----------------------------------------------------------------------

   ARRAY = 0._real_kind

   this_block = get_block(1,1) ! for the tripoleTflag - all blocks have it
   if (this_block%tripoleTFlag) then
     select case (field_loc)
     case (field_loc_center)   ! cell center location
        xoffset = 2
        yoffset = 0
     case (field_loc_NEcorner) ! cell corner (velocity) location
        xoffset = 1
        yoffset = 1
     case (field_loc_Eface)    ! cell face location
        xoffset = 1
        yoffset = 0
     case (field_loc_Nface)    ! cell face location
        xoffset = 2
        yoffset = 1
     case (field_loc_noupdate) ! ghost cells never used - use cell center
        xoffset = 1
        yoffset = 1
     end select
   else
     select case (field_loc)
     case (field_loc_center)   ! cell center location
        xoffset = 1
        yoffset = 1
     case (field_loc_NEcorner) ! cell corner (velocity) location
        xoffset = 0
        yoffset = 0
     case (field_loc_Eface)    ! cell face location
        xoffset = 0
        yoffset = 1
     case (field_loc_Nface)    ! cell face location
        xoffset = 1
        yoffset = 0
     case (field_loc_noupdate) ! ghost cells never used - use cell center
        xoffset = 1
        yoffset = 1
     end select
   endif

   select case (field_type)
   case (field_type_scalar)
      isign =  1
   case (field_type_vector)
      isign = -1
   case (field_type_angle)
      isign = -1
   case (field_type_noupdate) ! ghost cells not needed - use cell center
      isign =  1
   case default
      call abort_ice(subname//'ERROR: Unknown field type in scatter')
   end select

!-----------------------------------------------------------------------
!
!  copy blocks of global array into local block distribution
!
!-----------------------------------------------------------------------

   do n=1,nblocks_tot

      if (dst_dist%blockLocation(n) /= 0) then

         this_block = get_block(n,n)
         dst_block  = dst_dist%blockLocalID(n)

         do j=1,ny_block
            if (this_block%jblock == nblocks_y .and. this_block%j_glob(j) < 0) then
               ! tripole is top block with j_glob < 0
               ! for yoffset=0 or 1, yoffset2=0,0
               ! for yoffset=-1, yoffset2=0,1, for u-rows on T-fold grid
               do yoffset2=0,max(yoffset,0)-yoffset
                 jsrc = ny_global + yoffset + yoffset2 + &
                      (this_block%j_glob(j) + ny_global)
                 do i=1,nx_block
                   if (this_block%i_glob(i) /= 0) then
                      isrc = nx_global + xoffset - this_block%i_glob(i)
                      if (isrc < 1) isrc = isrc + nx_global
                      if (isrc > nx_global) isrc = isrc - nx_global
                      ARRAY(i,j-yoffset2,dst_block) &
                        = isign * ARRAY_G(isrc,jsrc)
                   endif
                 end do
               end do
            else
               ! normal block
               do i=1,nx_block
                  isrc = this_block%i_glob(i)
                  jsrc = this_block%j_glob(j)
                  if (isrc >=1 .and. isrc <= nx_global .and. &
                      jsrc >=1 .and. jsrc <= ny_global) &
                     ARRAY(i,j,dst_block) = ARRAY_G(isrc,jsrc)
               end do
            endif
         end do
      endif ! dst block not land
   end do  ! block loop

   !-----------------------------------------------------------------
   ! Set ghost cell values to 0 for noupdate
   !-----------------------------------------------------------------

   if (field_loc == field_loc_noupdate) then
     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) /= 0 .and. &
           dst_dist%blockLocalID(n) > 0) then

         dst_block = dst_dist%blockLocalID(n)
         this_block = get_block(n,n)
         ! north edge
         do j = this_block%jhi+1,ny_block
         do i = 1, nx_block
            ARRAY (i,j,dst_block) = 0._real_kind
         enddo
         enddo
         ! east edge
         do j = 1, ny_block
         do i = this_block%ihi+1,nx_block
            ARRAY (i,j,dst_block) = 0._real_kind
         enddo
         enddo
         ! south edge
         do j = 1, this_block%jlo-1
         do i = 1, nx_block
            ARRAY (i,j,dst_block) = 0._real_kind
         enddo
         enddo
         ! west edge
         do j = 1, ny_block
         do i = 1, this_block%ilo-1
            ARRAY (i,j,dst_block) = 0._real_kind
         enddo
         enddo

       endif
     enddo
   endif

!-----------------------------------------------------------------------

 end subroutine scatter_global_real

!***********************************************************************

 subroutine scatter_global_int(ARRAY, ARRAY_G, src_task, dst_dist, &
                               field_loc, field_type)

!-----------------------------------------------------------------------
!
!  This subroutine scatters a global-sized array on the processor
!  src\_task to a distribution of blocks given by dst\_dist.
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!  input variables
!
!-----------------------------------------------------------------------

   integer (int_kind), intent(in) :: &
     src_task       ! task from which array should be scattered

   type (distrb), intent(in) :: &
     dst_dist       ! distribution of resulting blocks

   integer (int_kind), dimension(:,:), intent(in) :: &
     ARRAY_G        ! array containing global field on src_task

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      field_loc                  ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

!-----------------------------------------------------------------------
!
!  output variables
!
!-----------------------------------------------------------------------

   integer (int_kind), dimension(:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,n,             &! dummy loop indices
      isrc, jsrc,        &! source addresses
      xoffset, yoffset,  &! offsets for tripole boundary conditions
      isign,             &! sign factor for tripole boundary conditions
      yoffset2,          &!
      dst_block           ! local block index in dest distribution

   type (block) :: &
      this_block  ! block info for current block

   character(len=*), parameter :: subname = '(scatter_global_int)'

!-----------------------------------------------------------------------
!
!  initialize return array to zero and set up tripole quantities
!
!-----------------------------------------------------------------------

   ARRAY = 0

   this_block = get_block(1,1) ! for the tripoleTflag - all blocks have it
   if (this_block%tripoleTFlag) then
     select case (field_loc)
     case (field_loc_center)   ! cell center location
        xoffset = 2
        yoffset = 0
     case (field_loc_NEcorner) ! cell corner (velocity) location
        xoffset = 1
        yoffset = 1
     case (field_loc_Eface)    ! cell face location
        xoffset = 1
        yoffset = 0
     case (field_loc_Nface)    ! cell face location
        xoffset = 2
        yoffset = 1
     case (field_loc_noupdate) ! ghost cells never used - use cell center
        xoffset = 1
        yoffset = 1
     end select
   else
     select case (field_loc)
     case (field_loc_center)   ! cell center location
        xoffset = 1
        yoffset = 1
     case (field_loc_NEcorner) ! cell corner (velocity) location
        xoffset = 0
        yoffset = 0
     case (field_loc_Eface)    ! cell face location
        xoffset = 0
        yoffset = 1
     case (field_loc_Nface)    ! cell face location
        xoffset = 1
        yoffset = 0
     case (field_loc_noupdate) ! ghost cells never used - use cell center
        xoffset = 1
        yoffset = 1
     end select
   endif

   select case (field_type)
   case (field_type_scalar)
      isign =  1
   case (field_type_vector)
      isign = -1
   case (field_type_angle)
      isign = -1
   case (field_type_noupdate) ! ghost cells not needed - use cell center
      isign =  1
   case default
      call abort_ice(subname//'ERROR: Unknown field type in scatter')
   end select

!-----------------------------------------------------------------------
!
!  copy blocks of global array into local block distribution
!
!-----------------------------------------------------------------------

   do n=1,nblocks_tot

      if (dst_dist%blockLocation(n) /= 0) then

         this_block = get_block(n,n)
         dst_block  = dst_dist%blockLocalID(n)

         do j=1,ny_block
            if (this_block%jblock == nblocks_y .and. this_block%j_glob(j) < 0) then
               ! tripole is top block with j_glob < 0
               ! for yoffset=0 or 1, yoffset2=0,0
               ! for yoffset=-1, yoffset2=0,1, for u-rows on T-fold grid
               do yoffset2=0,max(yoffset,0)-yoffset
                 jsrc = ny_global + yoffset + yoffset2 + &
                      (this_block%j_glob(j) + ny_global)
                 do i=1,nx_block
                   if (this_block%i_glob(i) /= 0) then
                      isrc = nx_global + xoffset - this_block%i_glob(i)
                      if (isrc < 1) isrc = isrc + nx_global
                      if (isrc > nx_global) isrc = isrc - nx_global
                      ARRAY(i,j-yoffset2,dst_block) &
                        = isign * ARRAY_G(isrc,jsrc)
                   endif
                 end do
               end do
            else
               ! normal block
               do i=1,nx_block
                  isrc = this_block%i_glob(i)
                  jsrc = this_block%j_glob(j)
                  if (isrc >=1 .and. isrc <= nx_global .and. &
                      jsrc >=1 .and. jsrc <= ny_global) &
                     ARRAY(i,j,dst_block) = ARRAY_G(isrc,jsrc)
               end do
            endif
         end do
      endif ! dst block not land
   end do  ! block loop

   !-----------------------------------------------------------------
   ! Set ghost cell values to 0 for noupdate
   !-----------------------------------------------------------------

   if (field_loc == field_loc_noupdate) then
     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) /= 0 .and. &
           dst_dist%blockLocalID(n) > 0) then

         dst_block = dst_dist%blockLocalID(n)
         this_block = get_block(n,n)
         ! north edge
         do j = this_block%jhi+1,ny_block
         do i = 1, nx_block
            ARRAY (i,j,dst_block) = 0
         enddo
         enddo
         ! east edge
         do j = 1, ny_block
         do i = this_block%ihi+1,nx_block
            ARRAY (i,j,dst_block) = 0
         enddo
         enddo
         ! south edge
         do j = 1, this_block%jlo-1
         do i = 1, nx_block
            ARRAY (i,j,dst_block) = 0
         enddo
         enddo
         ! west edge
         do j = 1, ny_block
         do i = 1, this_block%ilo-1
            ARRAY (i,j,dst_block) = 0
         enddo
         enddo

         endif
      enddo
   endif

!-----------------------------------------------------------------------

 end subroutine scatter_global_int

!***********************************************************************

 subroutine scatter_global_ext(ARRAY, ARRAY_G, src_task, dst_dist)

!  This subroutine scatters a global-sized array on the processor
!  src\_task to a distribution of blocks given by dst\_dist.
!
!  This is the specific interface for double precision arrays
!  corresponding to the generic interface scatter_global.  It is shown
!  to provide information on the generic interface (the generic
!  interface is identical, but chooses a specific interface based
!  on the data type of the input argument).

   integer (int_kind), intent(in) :: &
     src_task       ! task from which array should be scattered

   type (distrb), intent(in) :: &
     dst_dist       ! distribution of resulting blocks

   real (dbl_kind), dimension(:,:), intent(in) :: &
     ARRAY_G        ! array containing global field on src_task

   real (dbl_kind), dimension(:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n,              &! dummy loop indices
     isrc, jsrc,         &! source addresses
     iblk, jblk,         &! source addresses
     dst_block            ! local block index in dest distribution

   type (block) :: &
     this_block  ! block info for current block

   character(len=*), parameter :: subname = '(scatter_global_ext)'

!-----------------------------------------------------------------------
!
!  initialize return array to zero
!
!-----------------------------------------------------------------------

   ARRAY = c0

!-----------------------------------------------------------------------
!
!  copy blocks of global array into local block distribution
!
!-----------------------------------------------------------------------

   do n=1,nblocks_tot

      if (dst_dist%blockLocation(n) /= 0) then

         this_block = get_block(n,n)
         dst_block  = dst_dist%blockLocalID(n)

         ! interior
         do j = 1, ny_block
         do i = 1, nx_block
            ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i)+nghost,&
                                           this_block%j_glob(j)+nghost)
         end do
         end do

         if (this_block%jblock == 1) then
            ! south edge
            do j = 1, nghost
            do i = this_block%ilo,this_block%ihi
               ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i)+nghost,j)
            enddo
               do i = 1, nghost
               ! southwest corner
                  iblk = i
                  jblk = j
                  isrc = this_block%i_glob(this_block%ilo)+i-1
                  jsrc = j
                  ARRAY(iblk,jblk,dst_block) = ARRAY_G(isrc,jsrc)
               ! southeast corner
                  iblk = this_block%ihi+i
                  isrc = this_block%i_glob(this_block%ihi)+nghost+i
                  ARRAY(iblk,jblk,dst_block) = ARRAY_G(isrc,jsrc)
               enddo
            enddo
         endif
         if (this_block%jblock == nblocks_y) then
            ! north edge
            do j = 1, nghost
            do i = this_block%ilo,this_block%ihi
               ARRAY(i,this_block%jhi+j,dst_block) = ARRAY_G(this_block%i_glob(i)+nghost,&
                                                             ny_global+nghost+j)
            enddo
               do i = 1, nghost
               ! northwest corner
                  iblk = i
                  jblk = this_block%jhi+j
                  isrc = this_block%i_glob(this_block%ilo)+i-1
                  jsrc = ny_global+nghost+j
                  ARRAY(iblk,jblk,dst_block) = ARRAY_G(isrc,jsrc)
               ! northeast corner
                  iblk = this_block%ihi+i
                  isrc = this_block%i_glob(this_block%ihi)+nghost+i
                  ARRAY(iblk,jblk,dst_block) = ARRAY_G(isrc,jsrc)
               enddo
            enddo
         endif
         if (this_block%iblock == 1) then
            ! west edge
            do j = this_block%jlo,this_block%jhi
            do i = 1, nghost
               ARRAY(i,j,dst_block) = ARRAY_G(i,this_block%j_glob(j)+nghost)
            enddo
            enddo
               do j = 1, nghost
               do i = 1, nghost
               ! northwest corner
                  iblk = i
                  jblk = this_block%jhi+j
                  isrc = i
                  jsrc = this_block%j_glob(this_block%jhi)+nghost+j
                  ARRAY(iblk,jblk,dst_block) = ARRAY_G(isrc,jsrc)
               ! southwest corner
                  jblk = j
                  jsrc = this_block%j_glob(this_block%jlo)+j-1
                  ARRAY(iblk,jblk,dst_block) = ARRAY_G(isrc,jsrc)
               enddo
               enddo
         endif
         if (this_block%iblock == nblocks_x) then
            ! east edge
            do j = this_block%jlo,this_block%jhi
            do i = 1, nghost
               ARRAY(this_block%ihi+i,j,dst_block) = ARRAY_G(nx_global+nghost+i, &
                                                             this_block%j_glob(j)+nghost)
            enddo
            enddo
               do j = 1, nghost
               do i = 1, nghost
               ! northeast corner
                  iblk = this_block%ihi+i
                  jblk = this_block%jhi+j
                  isrc = nx_global+nghost+i
                  jsrc = this_block%j_glob(this_block%jhi)+nghost+j
                  ARRAY(iblk,jblk,dst_block) = ARRAY_G(isrc,jsrc)
               ! southeast corner
                  jblk = j
                  jsrc = this_block%j_glob(this_block%jlo)+j-1
                  ARRAY(iblk,jblk,dst_block) = ARRAY_G(isrc,jsrc)
               enddo
               enddo
         endif
      endif ! dst block not land
   end do  ! block loop

!-----------------------------------------------------------------------

 end subroutine scatter_global_ext

!***********************************************************************

 subroutine scatter_global_stress(ARRAY, ARRAY_G1, ARRAY_G2, &
                                  src_task, dst_dist)

!  This subroutine scatters global stresses to a distributed array.
!
!  Ghost cells in the stress tensor must be handled separately on tripole
!  grids, because matching the corner values requires 2 different arrays.

   integer (int_kind), intent(in) :: &
     src_task       ! task from which array should be scattered

   type (distrb), intent(in) :: &
     dst_dist       ! distribution of resulting blocks

   real (dbl_kind), dimension(:,:), intent(in) :: &
     ARRAY_G1,     &! array containing global field on src_task
     ARRAY_G2       ! array containing global field on src_task

   real (dbl_kind), dimension(:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n,              &! dummy loop indices
     isrc, jsrc,         &! source addresses
     xoffset, yoffset,   &! offsets for tripole boundary conditions
     yoffset2,           &!
     isign,              &! sign factor for tripole boundary conditions
     dst_block            ! local block index in dest distribution

   type (block) :: &
     this_block  ! block info for current block

   character(len=*), parameter :: subname = '(scatter_global_stress)'

!-----------------------------------------------------------------------
!
!  initialize return array to zero and set up tripole quantities
!
!-----------------------------------------------------------------------

   ARRAY = c0

   this_block = get_block(1,1) ! for the tripoleTflag - all blocks have it
   if (this_block%tripoleTFlag) then
     xoffset = 2  ! treat stresses as cell-centered scalars (they are not
     yoffset = 0  ! shared with neighboring grid cells)
   else
     xoffset = 1  ! treat stresses as cell-centered scalars (they are not
     yoffset = 1  ! shared with neighboring grid cells)
   endif
   isign   = 1

!-----------------------------------------------------------------------
!
!  copy blocks of global array into local block distribution
!
!-----------------------------------------------------------------------

   do n=1,nblocks_tot

      if (dst_dist%blockLocation(n) /= 0) then

         this_block = get_block(n,n)
         dst_block  = dst_dist%blockLocalID(n)

         do j=1,ny_block
            if (this_block%jblock == nblocks_y .and. this_block%j_glob(j) < 0) then
               ! tripole is top block with j_glob < 0
               ! for yoffset=0 or 1, yoffset2=0,0
               ! for yoffset=-1, yoffset2=0,1, for u-rows on T-fold grid
               do yoffset2=0,max(yoffset,0)-yoffset
                 jsrc = ny_global + yoffset + yoffset2 + &
                      (this_block%j_glob(j) + ny_global)
                 do i=1,nx_block
                   if (this_block%i_glob(i) /= 0) then
                      isrc = nx_global + xoffset - this_block%i_glob(i)
                      if (isrc < 1) isrc = isrc + nx_global
                      if (isrc > nx_global) isrc = isrc - nx_global
                      ARRAY(i,j-yoffset2,dst_block) &
                        = isign * ARRAY_G2(isrc,jsrc)
                   endif
                 end do
               end do
            else
               ! normal block
               do i=1,nx_block
                  isrc = this_block%i_glob(i)
                  jsrc = this_block%j_glob(j)
                  if (isrc >=1 .and. isrc <= nx_global .and. &
                      jsrc >=1 .and. jsrc <= ny_global) &
                     ARRAY(i,j,dst_block) = ARRAY_G1(isrc,jsrc)
               end do
            endif
         end do
      endif ! dst block not land
   end do  ! block loop

!-----------------------------------------------------------------------

 end subroutine scatter_global_stress

!***********************************************************************

 end module ice_gather_scatter

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
