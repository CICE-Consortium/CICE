!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module ice_gather_scatter

!  This module contains routines for gathering data to a single
!  processor from a distributed array, and scattering data from a
!  single processor to a distributed array.
!
!  NOTE: The arrays gathered and scattered are assumed to have
!        horizontal dimensions (nx_block, ny_block).
!
! author: Phil Jones, LANL
! Oct. 2004: Adapted from POP version by William H. Lipscomb, LANL
! Jan. 2008: Elizabeth Hunke replaced old routines with new POP
!              infrastructure, added specialized routine scatter_global_stress

   use mpi   ! MPI Fortran module
   use ice_kinds_mod
   use ice_communicate, only: my_task, mpiR8, mpiR4, mpitag_gs, MPI_COMM_ICE, &
       ice_barrier, add_mpi_barriers
   use ice_constants, only: spval_dbl, c0, &
       field_loc_center, field_loc_NEcorner, field_loc_Nface, field_loc_Eface, &
       field_loc_noupdate, &
       field_type_scalar, field_type_vector, field_type_angle, &
       field_type_noupdate
   use ice_blocks, only: block, nx_block, ny_block, nblocks_tot, get_block, &
       nblocks_x, nblocks_y, nghost
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
     i,j,n          ,&! dummy loop counters
     nsends         ,&! number of actual sends
     src_block      ,&! block locator for send
     ierr             ! MPI error flag

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status

   integer (int_kind), dimension(:), allocatable :: &
     snd_request

   integer (int_kind), dimension(:,:), allocatable :: &
     snd_status

   real (dbl_kind), dimension(:,:), allocatable :: &
     msg_buffer

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
!  if this task is the dst_task, copy local blocks into the global
!  array and post receives for non-local blocks.
!
!-----------------------------------------------------------------------

   if (my_task == dst_task) then

     do n=1,nblocks_tot

       !*** copy local blocks

       if (src_dist%blockLocation(n) == my_task+1) then

         this_block = get_block(n,n)

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = &
                  ARRAY(i,j,src_dist%blockLocalID(n))
         end do
         end do

       !*** fill land blocks with special values

       else if (src_dist%blockLocation(n) == 0) then

         this_block = get_block(n,n)

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = special_value
         end do
         end do
       endif

     end do

     !*** receive blocks to fill up the rest

     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (src_dist%blockLocation(n) > 0 .and. &
           src_dist%blockLocation(n) /= my_task+1) then

         this_block = get_block(n,n)

         call MPI_RECV(msg_buffer, size(msg_buffer), &
                       mpiR8, src_dist%blockLocation(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, status, ierr)

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = msg_buffer(i,j)
         end do
         end do
       endif
     end do

     deallocate(msg_buffer)

!-----------------------------------------------------------------------
!
!  otherwise send data to dst_task
!
!-----------------------------------------------------------------------

   else

     allocate(snd_request(nblocks_tot), &
              snd_status (MPI_STATUS_SIZE, nblocks_tot))

     nsends = 0
     do n=1,nblocks_tot
       if (src_dist%blockLocation(n) == my_task+1) then

         nsends = nsends + 1
         src_block = src_dist%blockLocalID(n)
         call MPI_ISEND(ARRAY(1,1,src_block), nx_block*ny_block, &
                     mpiR8, dst_task, 3*mpitag_gs+n, &
                     MPI_COMM_ICE, snd_request(nsends), ierr)
       endif
     end do

     if (nsends > 0) &
       call MPI_WAITALL(nsends, snd_request, snd_status, ierr)
     deallocate(snd_request, snd_status)

   endif

   if (add_mpi_barriers) then
     call ice_barrier()
   endif

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
     i,j,n          ,&! dummy loop counters
     nsends         ,&! number of actual sends
     src_block      ,&! block locator for send
     ierr             ! MPI error flag

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status

   integer (int_kind), dimension(:), allocatable :: &
     snd_request

   integer (int_kind), dimension(:,:), allocatable :: &
     snd_status

   real (real_kind), dimension(:,:), allocatable :: &
     msg_buffer

   type (block) :: &
     this_block  ! block info for current block

   character(len=*), parameter :: subname = '(gather_global_real)'

!-----------------------------------------------------------------------
!
!  if this task is the dst_task, copy local blocks into the global
!  array and post receives for non-local blocks.
!
!-----------------------------------------------------------------------

   if (my_task == dst_task) then

     do n=1,nblocks_tot

       !*** copy local blocks

       if (src_dist%blockLocation(n) == my_task+1) then

         this_block = get_block(n,n)

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = &
                  ARRAY(i,j,src_dist%blockLocalID(n))
         end do
         end do

       !*** fill land blocks with zeroes

       else if (src_dist%blockLocation(n) == 0) then

         this_block = get_block(n,n)

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = 0._real_kind
         end do
         end do
       endif

     end do

     !*** receive blocks to fill up the rest

     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (src_dist%blockLocation(n) > 0 .and. &
           src_dist%blockLocation(n) /= my_task+1) then

         this_block = get_block(n,n)

         call MPI_RECV(msg_buffer, size(msg_buffer), &
                       mpiR4, src_dist%blockLocation(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, status, ierr)

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = msg_buffer(i,j)
         end do
         end do
       endif
     end do

     deallocate(msg_buffer)

!-----------------------------------------------------------------------
!
!  otherwise send data to dst_task
!
!-----------------------------------------------------------------------

   else

     allocate(snd_request(nblocks_tot), &
              snd_status (MPI_STATUS_SIZE, nblocks_tot))

     nsends = 0
     do n=1,nblocks_tot
       if (src_dist%blockLocation(n) == my_task+1) then

         nsends = nsends + 1
         src_block = src_dist%blockLocalID(n)
         call MPI_ISEND(ARRAY(1,1,src_block), nx_block*ny_block, &
                     mpiR4, dst_task, 3*mpitag_gs+n, &
                     MPI_COMM_ICE, snd_request(nsends), ierr)
       endif
     end do

     if (nsends > 0) &
       call MPI_WAITALL(nsends, snd_request, snd_status, ierr)
     deallocate(snd_request, snd_status)

   endif

   if (add_mpi_barriers) then
     call ice_barrier()
   endif

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
     i,j,n          ,&! dummy loop counters
     nsends         ,&! number of actual sends
     src_block      ,&! block locator for send
     ierr             ! MPI error flag

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status

   integer (int_kind), dimension(:), allocatable :: &
     snd_request

   integer (int_kind), dimension(:,:), allocatable :: &
     snd_status

   integer (int_kind), dimension(:,:), allocatable :: &
     msg_buffer

   type (block) :: &
     this_block  ! block info for current block

   character(len=*), parameter :: subname = '(gather_global_int)'

!-----------------------------------------------------------------------
!
!  if this task is the dst_task, copy local blocks into the global
!  array and post receives for non-local blocks.
!
!-----------------------------------------------------------------------

   if (my_task == dst_task) then

     do n=1,nblocks_tot

       !*** copy local blocks

       if (src_dist%blockLocation(n) == my_task+1) then

         this_block = get_block(n,n)

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = &
                  ARRAY(i,j,src_dist%blockLocalID(n))
         end do
         end do

       !*** fill land blocks with zeroes

       else if (src_dist%blockLocation(n) == 0) then

         this_block = get_block(n,n)

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = 0
         end do
         end do
       endif

     end do

     !*** receive blocks to fill up the rest

     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (src_dist%blockLocation(n) > 0 .and. &
           src_dist%blockLocation(n) /= my_task+1) then

         this_block = get_block(n,n)

         call MPI_RECV(msg_buffer, size(msg_buffer), &
                       mpi_integer, src_dist%blockLocation(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, status, ierr)

         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = msg_buffer(i,j)
         end do
         end do
       endif
     end do

     deallocate(msg_buffer)

!-----------------------------------------------------------------------
!
!  otherwise send data to dst_task
!
!-----------------------------------------------------------------------

   else

     allocate(snd_request(nblocks_tot), &
              snd_status (MPI_STATUS_SIZE, nblocks_tot))

     nsends = 0
     do n=1,nblocks_tot
       if (src_dist%blockLocation(n) == my_task+1) then

         nsends = nsends + 1
         src_block = src_dist%blockLocalID(n)
         call MPI_ISEND(ARRAY(1,1,src_block), nx_block*ny_block, &
                     mpi_integer, dst_task, 3*mpitag_gs+n, &
                     MPI_COMM_ICE, snd_request(nsends), ierr)
       endif
     end do

     if (nsends > 0) &
       call MPI_WAITALL(nsends, snd_request, snd_status, ierr)
     deallocate(snd_request, snd_status)

   endif

   if (add_mpi_barriers) then
     call ice_barrier()
   endif

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

   real (dbl_kind), dimension(:,:), intent(inout) :: &
     ARRAY_G    ! array containing global horizontal field on dst_task

   real (dbl_kind), intent(in), optional :: &
     spc_val

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n          ,&! dummy loop counters
     nx, ny         ,&! global dimensions
     nsends         ,&! number of actual sends
     src_block      ,&! block locator for send
     ierr             ! MPI error flag

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status

   integer (int_kind), dimension(:), allocatable :: &
     snd_request

   integer (int_kind), dimension(:,:), allocatable :: &
     snd_status

   real (dbl_kind), dimension(:,:), allocatable :: &
     msg_buffer

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
!  if this task is the dst_task, copy local blocks into the global
!  array and post receives for non-local blocks.
!
!-----------------------------------------------------------------------

   if (my_task == dst_task) then

     do n=1,nblocks_tot

       !*** copy local blocks

       if (src_dist%blockLocation(n) == my_task+1) then

         this_block = get_block(n,n)

         ! interior
         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
           ARRAY_G(this_block%i_glob(i)+nghost, &
                   this_block%j_glob(j)+nghost) = &
           ARRAY  (i,j,src_dist%blockLocalID(n))
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

     !*** receive blocks to fill up the rest

     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (src_dist%blockLocation(n) > 0 .and. &
           src_dist%blockLocation(n) /= my_task+1) then

         this_block = get_block(n,n)

         call MPI_RECV(msg_buffer, size(msg_buffer), &
                       mpiR8, src_dist%blockLocation(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, status, ierr)

         ! block interior
         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
           ARRAY_G(this_block%i_glob(i)+nghost, &
                   this_block%j_glob(j)+nghost) = msg_buffer(i,j)
         end do
         end do
         if (this_block%jblock == 1) then
            ! south block
            do j=1, nghost
            do i=this_block%ilo,this_block%ihi
              ARRAY_G   (this_block%i_glob(i)+nghost,j) = &
              msg_buffer(i,j)
            end do
            end do
            if (this_block%iblock == 1) then
               ! southwest corner
               do j=1, nghost
               do i=1, nghost
                 ARRAY_G(i,j) = msg_buffer(i,j)
               end do
               end do
            endif
         endif
         if (this_block%jblock == nblocks_y) then
            ! north block
            do j=1, nghost
            do i=this_block%ilo,this_block%ihi
              ARRAY_G   (this_block%i_glob(i)+nghost, &
                         ny_global + nghost + j) = &
              msg_buffer(i, this_block%jhi+j)
            end do
            end do
            if (this_block%iblock == nblocks_x) then
               ! northeast corner
               do j=1, nghost
               do i=1, nghost
                 ARRAY_G   (nx-i+1, ny-j+1) = &
                 msg_buffer(this_block%ihi+nghost-i+1,&
                            this_block%jhi+nghost-j+1)
               end do
               end do
            endif
         endif
         if (this_block%iblock == 1) then
            ! west block
            do j=this_block%jlo,this_block%jhi
            do i=1, nghost
              ARRAY_G   (i, this_block%j_glob(j)+nghost) = &
              msg_buffer(i, j)
            end do
            end do
            if (this_block%jblock == nblocks_y) then
               ! northwest corner
               do j=1, nghost
               do i=1, nghost
                 ARRAY_G   (i, ny-j+1) = &
                 msg_buffer(i, this_block%jhi+nghost-j+1)
               end do
               end do
            endif
         endif
         if (this_block%iblock == nblocks_x) then
            ! east block
            do j=this_block%jlo,this_block%jhi
            do i=1, nghost
              ARRAY_G   (nx_global+nghost+i, &
                         this_block%j_glob(j)+nghost) = &
              msg_buffer(this_block%ihi+i, j)
            end do
            end do
            if (this_block%jblock == 1) then
               ! southeast corner
               do j=1, nghost
               do i=1, nghost
                 ARRAY_G   (nx-i+1, j) = &
                 msg_buffer(this_block%ihi+nghost-i+1, j)
               end do
               end do
            endif
         endif
       endif
     end do

     deallocate(msg_buffer)

!-----------------------------------------------------------------------
!
!  otherwise send data to dst_task
!
!-----------------------------------------------------------------------

   else  ! master task

     allocate(snd_request(nblocks_tot), &
              snd_status (MPI_STATUS_SIZE, nblocks_tot))

     nsends = 0
     do n=1,nblocks_tot
       if (src_dist%blockLocation(n) == my_task+1) then

         nsends = nsends + 1
         src_block = src_dist%blockLocalID(n)
         call MPI_ISEND(ARRAY(1,1,src_block), nx_block*ny_block, &
                     mpiR8, dst_task, 3*mpitag_gs+n, &
                     MPI_COMM_ICE, snd_request(nsends), ierr)
       endif
     end do

     if (nsends > 0) &
       call MPI_WAITALL(nsends, snd_request, snd_status, ierr)
     deallocate(snd_request, snd_status)

   endif  ! master task

   if (add_mpi_barriers) then
     call ice_barrier()
   endif

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

   integer (int_kind), dimension(:,:), intent(inout) :: &
     ARRAY_G    ! array containing global horizontal field on dst_task

   integer (int_kind), intent(in), optional :: &
     spc_val

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n          ,&! dummy loop counters
     nx, ny         ,&! global dimensions
     nsends         ,&! number of actual sends
     src_block      ,&! block locator for send
     ierr             ! MPI error flag

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status

   integer (int_kind), dimension(:), allocatable :: &
     snd_request

   integer (int_kind), dimension(:,:), allocatable :: &
     snd_status

   integer (int_kind), dimension(:,:), allocatable :: &
     msg_buffer

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
!  if this task is the dst_task, copy local blocks into the global
!  array and post receives for non-local blocks.
!
!-----------------------------------------------------------------------

   if (my_task == dst_task) then

     do n=1,nblocks_tot

       !*** copy local blocks

       if (src_dist%blockLocation(n) == my_task+1) then

         this_block = get_block(n,n)

         ! interior
         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
           ARRAY_G(this_block%i_glob(i)+nghost, &
                   this_block%j_glob(j)+nghost) = &
           ARRAY  (i,j,src_dist%blockLocalID(n))
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

     !*** receive blocks to fill up the rest

     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (src_dist%blockLocation(n) > 0 .and. &
           src_dist%blockLocation(n) /= my_task+1) then

         this_block = get_block(n,n)

         call MPI_RECV(msg_buffer, size(msg_buffer), &
                       MPI_INTEGER, src_dist%blockLocation(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, status, ierr)

         ! block interior
         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
           ARRAY_G(this_block%i_glob(i)+nghost, &
                   this_block%j_glob(j)+nghost) = msg_buffer(i,j)
         end do
         end do
         if (this_block%jblock == 1) then
            ! south block
            do j=1, nghost
            do i=this_block%ilo,this_block%ihi
              ARRAY_G   (this_block%i_glob(i)+nghost,j) = &
              msg_buffer(i,j)
            end do
            end do
            if (this_block%iblock == 1) then
               ! southwest corner
               do j=1, nghost
               do i=1, nghost
                 ARRAY_G(i,j) = msg_buffer(i,j)
               end do
               end do
            endif
         endif
         if (this_block%jblock == nblocks_y) then
            ! north block
            do j=1, nghost
            do i=this_block%ilo,this_block%ihi
              ARRAY_G   (this_block%i_glob(i)+nghost, &
                         ny_global + nghost + j) = &
              msg_buffer(i, this_block%jhi+j)
            end do
            end do
            if (this_block%iblock == nblocks_x) then
               ! northeast corner
               do j=1, nghost
               do i=1, nghost
                 ARRAY_G   (nx-i+1, ny-j+1) = &
                 msg_buffer(this_block%ihi+nghost-i+1,&
                            this_block%jhi+nghost-j+1)
               end do
               end do
            endif
         endif
         if (this_block%iblock == 1) then
            ! west block
            do j=this_block%jlo,this_block%jhi
            do i=1, nghost
              ARRAY_G   (i, this_block%j_glob(j)+nghost) = &
              msg_buffer(i, j)
            end do
            end do
            if (this_block%jblock == nblocks_y) then
               ! northwest corner
               do j=1, nghost
               do i=1, nghost
                 ARRAY_G   (i, ny-j+1) = &
                 msg_buffer(i, this_block%jhi+nghost-j+1)
               end do
               end do
            endif
         endif
         if (this_block%iblock == nblocks_x) then
            ! east block
            do j=this_block%jlo,this_block%jhi
            do i=1, nghost
              ARRAY_G   (nx_global+nghost+i, &
                         this_block%j_glob(j)+nghost) = &
              msg_buffer(this_block%ihi+i, j)
            end do
            end do
            if (this_block%jblock == 1) then
               ! southeast corner
               do j=1, nghost
               do i=1, nghost
                 ARRAY_G   (nx-i+1, j) = &
                 msg_buffer(this_block%ihi+nghost-i+1, j)
               end do
               end do
            endif
         endif
       endif
     end do

     deallocate(msg_buffer)

!-----------------------------------------------------------------------
!
!  otherwise send data to dst_task
!
!-----------------------------------------------------------------------

   else  ! master task

     allocate(snd_request(nblocks_tot), &
              snd_status (MPI_STATUS_SIZE, nblocks_tot))

     nsends = 0
     do n=1,nblocks_tot
       if (src_dist%blockLocation(n) == my_task+1) then

         nsends = nsends + 1
         src_block = src_dist%blockLocalID(n)
         call MPI_ISEND(ARRAY(1,1,src_block), nx_block*ny_block, &
                     MPI_INTEGER, dst_task, 3*mpitag_gs+n, &
                     MPI_COMM_ICE, snd_request(nsends), ierr)
       endif
     end do

     if (nsends > 0) &
       call MPI_WAITALL(nsends, snd_request, snd_status, ierr)
     deallocate(snd_request, snd_status)

   endif  ! master task

   if (add_mpi_barriers) then
     call ice_barrier()
   endif

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

   logical (log_kind), dimension(:,:), intent(inout) :: &
     ARRAY_G    ! array containing global horizontal field on dst_task

   logical (log_kind), intent(in), optional :: &
     spc_val

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n          ,&! dummy loop counters
     nx, ny         ,&! global dimensions
     nsends         ,&! number of actual sends
     src_block      ,&! block locator for send
     ierr             ! MPI error flag

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status

   integer (int_kind), dimension(:), allocatable :: &
     snd_request

   integer (int_kind), dimension(:,:), allocatable :: &
     snd_status

   logical (log_kind), dimension(:,:), allocatable :: &
     msg_buffer

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
!  if this task is the dst_task, copy local blocks into the global
!  array and post receives for non-local blocks.
!
!-----------------------------------------------------------------------

   if (my_task == dst_task) then

     do n=1,nblocks_tot

       !*** copy local blocks

       if (src_dist%blockLocation(n) == my_task+1) then

         this_block = get_block(n,n)

         ! interior
         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
           ARRAY_G(this_block%i_glob(i)+nghost, &
                   this_block%j_glob(j)+nghost) = &
           ARRAY  (i,j,src_dist%blockLocalID(n))
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

     !*** receive blocks to fill up the rest

     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (src_dist%blockLocation(n) > 0 .and. &
           src_dist%blockLocation(n) /= my_task+1) then

         this_block = get_block(n,n)

         call MPI_RECV(msg_buffer, size(msg_buffer), &
                       MPI_LOGICAL, src_dist%blockLocation(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, status, ierr)

         ! block interior
         do j=this_block%jlo,this_block%jhi
         do i=this_block%ilo,this_block%ihi
           ARRAY_G(this_block%i_glob(i)+nghost, &
                   this_block%j_glob(j)+nghost) = msg_buffer(i,j)
         end do
         end do
         if (this_block%jblock == 1) then
            ! south block
            do j=1, nghost
            do i=this_block%ilo,this_block%ihi
              ARRAY_G   (this_block%i_glob(i)+nghost,j) = &
              msg_buffer(i,j)
            end do
            end do
            if (this_block%iblock == 1) then
               ! southwest corner
               do j=1, nghost
               do i=1, nghost
                 ARRAY_G(i,j) = msg_buffer(i,j)
               end do
               end do
            endif
         endif
         if (this_block%jblock == nblocks_y) then
            ! north block
            do j=1, nghost
            do i=this_block%ilo,this_block%ihi
              ARRAY_G   (this_block%i_glob(i)+nghost, &
                         ny_global + nghost + j) = &
              msg_buffer(i, this_block%jhi+j)
            end do
            end do
            if (this_block%iblock == nblocks_x) then
               ! northeast corner
               do j=1, nghost
               do i=1, nghost
                 ARRAY_G   (nx-i+1, ny-j+1) = &
                 msg_buffer(this_block%ihi+nghost-i+1,&
                            this_block%jhi+nghost-j+1)
               end do
               end do
            endif
         endif
         if (this_block%iblock == 1) then
            ! west block
            do j=this_block%jlo,this_block%jhi
            do i=1, nghost
              ARRAY_G   (i, this_block%j_glob(j)+nghost) = &
              msg_buffer(i, j)
            end do
            end do
            if (this_block%jblock == nblocks_y) then
               ! northwest corner
               do j=1, nghost
               do i=1, nghost
                 ARRAY_G   (i, ny-j+1) = &
                 msg_buffer(i, this_block%jhi+nghost-j+1)
               end do
               end do
            endif
         endif
         if (this_block%iblock == nblocks_x) then
            ! east block
            do j=this_block%jlo,this_block%jhi
            do i=1, nghost
              ARRAY_G   (nx_global+nghost+i, &
                         this_block%j_glob(j)+nghost) = &
              msg_buffer(this_block%ihi+i, j)
            end do
            end do
            if (this_block%jblock == 1) then
               ! southeast corner
               do j=1, nghost
               do i=1, nghost
                 ARRAY_G   (nx-i+1, j) = &
                 msg_buffer(this_block%ihi+nghost-i+1, j)
               end do
               end do
            endif
         endif
       endif
     end do

     deallocate(msg_buffer)

!-----------------------------------------------------------------------
!
!  otherwise send data to dst_task
!
!-----------------------------------------------------------------------

   else  ! master task

     allocate(snd_request(nblocks_tot), &
              snd_status (MPI_STATUS_SIZE, nblocks_tot))

     nsends = 0
     do n=1,nblocks_tot
       if (src_dist%blockLocation(n) == my_task+1) then

         nsends = nsends + 1
         src_block = src_dist%blockLocalID(n)
         call MPI_ISEND(ARRAY(1,1,src_block), nx_block*ny_block, &
                     MPI_LOGICAL, dst_task, 3*mpitag_gs+n, &
                     MPI_COMM_ICE, snd_request(nsends), ierr)
       endif
     end do

     if (nsends > 0) &
       call MPI_WAITALL(nsends, snd_request, snd_status, ierr)
     deallocate(snd_request, snd_status)

   endif  ! master task

   if (add_mpi_barriers) then
     call ice_barrier()
   endif

!-----------------------------------------------------------------------

 end subroutine gather_global_ext_log

!***********************************************************************

 subroutine scatter_global_dbl(ARRAY, ARRAY_G, src_task, dst_dist, &
                               field_loc, field_type)

!  This subroutine scatters a global-sized array to a distributed array.
!
!  This is the specific interface for double precision arrays
!  corresponding to the generic interface scatter_global.

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
     nrecvs,             &! actual number of messages received
     isrc, jsrc,         &! source addresses
     dst_block,          &! location of block in dst array
     xoffset, yoffset,   &! offsets for tripole boundary conditions
     yoffset2,           &!
     isign,              &! sign factor for tripole boundary conditions
     ierr                 ! MPI error flag

   type (block) :: &
     this_block  ! block info for current block

   integer (int_kind), dimension(:), allocatable :: &
     rcv_request     ! request array for receives

   integer (int_kind), dimension(:,:), allocatable :: &
     rcv_status      ! status array for receives

   real (dbl_kind), dimension(:,:), allocatable :: &
     msg_buffer      ! buffer for sending blocks

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
   case (field_type_noupdate) ! ghost cells never used - use cell center
      isign =  1
   case default
      call abort_ice(subname//'ERROR: Unknown field type in scatter')
   end select

!-----------------------------------------------------------------------
!
!  if this task is the src_task, copy blocks of global array into
!  message buffer and send to other processors. also copy local blocks
!
!-----------------------------------------------------------------------

   if (my_task == src_task) then

     !*** send non-local blocks away

     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) > 0 .and. &
           dst_dist%blockLocation(n)-1 /= my_task) then

         msg_buffer = c0
         this_block = get_block(n,n)

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
                      msg_buffer(i,j-yoffset2) &
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
                     msg_buffer(i,j) = ARRAY_G(isrc,jsrc)
               end do
            endif
         end do

         call MPI_SEND(msg_buffer, nx_block*ny_block, &
                       mpiR8, dst_dist%blockLocation(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, ierr)

       endif
     end do

     deallocate(msg_buffer)

     !*** copy any local blocks

     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) == my_task+1) then
         dst_block = dst_dist%blockLocalID(n)
         this_block = get_block(n,n)

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
       endif
     end do

!-----------------------------------------------------------------------
!
!  otherwise receive data from src_task
!
!-----------------------------------------------------------------------

   else

     allocate (rcv_request(nblocks_tot), &
               rcv_status(MPI_STATUS_SIZE, nblocks_tot))

     rcv_request = 0
     rcv_status  = 0

     nrecvs = 0
     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) == my_task+1) then
         nrecvs = nrecvs + 1
         dst_block = dst_dist%blockLocalID(n)
         call MPI_IRECV(ARRAY(1,1,dst_block), nx_block*ny_block, &
                       mpiR8, src_task, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, rcv_request(nrecvs), ierr)
       endif
     end do

     if (nrecvs > 0) &
       call MPI_WAITALL(nrecvs, rcv_request, rcv_status, ierr)

     deallocate(rcv_request, rcv_status)
   endif

   !-----------------------------------------------------------------
   ! Set ghost cell values to 0 for noupdate
   !-----------------------------------------------------------------

   if (field_loc == field_loc_noupdate) then
     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) == my_task+1 .and. &
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

   if (add_mpi_barriers) then
     call ice_barrier()
   endif

!-----------------------------------------------------------------------

 end subroutine scatter_global_dbl

!***********************************************************************

 subroutine scatter_global_real(ARRAY, ARRAY_G, src_task, dst_dist, &
                               field_loc, field_type)

!-----------------------------------------------------------------------
!
!  This subroutine scatters a global-sized array to a distributed array.
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
     nrecvs,             &! actual number of messages received
     isrc, jsrc,         &! source addresses
     dst_block,          &! location of block in dst array
     xoffset, yoffset,   &! offsets for tripole boundary conditions
     yoffset2,           &!
     isign,              &! sign factor for tripole boundary conditions
     ierr                 ! MPI error flag

   type (block) :: &
     this_block  ! block info for current block

   integer (int_kind), dimension(:), allocatable :: &
     rcv_request     ! request array for receives

   integer (int_kind), dimension(:,:), allocatable :: &
     rcv_status      ! status array for receives

   real (real_kind), dimension(:,:), allocatable :: &
     msg_buffer      ! buffer for sending blocks

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
   case (field_type_noupdate) ! ghost cells never used - use cell center
      isign =  1
   case default
      call abort_ice(subname//'ERROR: Unknown field type in scatter')
   end select

!-----------------------------------------------------------------------
!
!  if this task is the src_task, copy blocks of global array into
!  message buffer and send to other processors. also copy local blocks
!
!-----------------------------------------------------------------------

   if (my_task == src_task) then

     !*** send non-local blocks away

     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) > 0 .and. &
           dst_dist%blockLocation(n)-1 /= my_task) then

         msg_buffer = 0._real_kind
         this_block = get_block(n,n)

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
                      msg_buffer(i,j-yoffset2) &
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
                     msg_buffer(i,j) = ARRAY_G(isrc,jsrc)
               end do
            endif
         end do

         call MPI_SEND(msg_buffer, nx_block*ny_block, &
                       mpiR4, dst_dist%blockLocation(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, ierr)

       endif
     end do

     deallocate(msg_buffer)

     !*** copy any local blocks

     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) == my_task+1) then
         dst_block = dst_dist%blockLocalID(n)
         this_block = get_block(n,n)

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
       endif
     end do

!-----------------------------------------------------------------------
!
!  otherwise receive data from src_task
!
!-----------------------------------------------------------------------

   else

     allocate (rcv_request(nblocks_tot), &
               rcv_status(MPI_STATUS_SIZE, nblocks_tot))

     rcv_request = 0
     rcv_status  = 0

     nrecvs = 0
     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) == my_task+1) then
         nrecvs = nrecvs + 1
         dst_block = dst_dist%blockLocalID(n)
         call MPI_IRECV(ARRAY(1,1,dst_block), nx_block*ny_block, &
                       mpiR4, src_task, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, rcv_request(nrecvs), ierr)
       endif
     end do

     if (nrecvs > 0) &
       call MPI_WAITALL(nrecvs, rcv_request, rcv_status, ierr)

     deallocate(rcv_request, rcv_status)
   endif

   !-----------------------------------------------------------------
   ! Set ghost cell values to 0 for noupdate
   !-----------------------------------------------------------------

   if (field_loc == field_loc_noupdate) then
     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) == my_task+1 .and. &
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

   if (add_mpi_barriers) then
     call ice_barrier()
   endif

!-----------------------------------------------------------------------

 end subroutine scatter_global_real

!***********************************************************************

 subroutine scatter_global_int(ARRAY, ARRAY_G, src_task, dst_dist, &
                               field_loc, field_type)

!-----------------------------------------------------------------------
!
!  This subroutine scatters a global-sized array to a distributed array.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  input variables
!
!-----------------------------------------------------------------------

   integer (int_kind), intent(in) :: &
     src_task       ! task from which array should be scattered

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      field_loc                  ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

   type (distrb), intent(in) :: &
     dst_dist       ! distribution of resulting blocks

   integer (int_kind), dimension(:,:), intent(in) :: &
     ARRAY_G        ! array containing global field on src_task

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
     i,j,n,              &! dummy loop indices
     nrecvs,             &! actual number of messages received
     isrc, jsrc,         &! source addresses
     dst_block,          &! location of block in dst array
     xoffset, yoffset,   &! offsets for tripole boundary conditions
     yoffset2,           &!
     isign,              &! sign factor for tripole boundary conditions
     ierr                 ! MPI error flag

   type (block) :: &
     this_block  ! block info for current block

   integer (int_kind), dimension(:), allocatable :: &
     rcv_request     ! request array for receives

   integer (int_kind), dimension(:,:), allocatable :: &
     rcv_status      ! status array for receives

   integer (int_kind), dimension(:,:), allocatable :: &
     msg_buffer      ! buffer for sending blocks

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
   case (field_type_noupdate) ! ghost cells never used - use cell center
      isign =  1
   case default
      call abort_ice(subname//'ERROR: Unknown field type in scatter')
   end select

!-----------------------------------------------------------------------
!
!  if this task is the src_task, copy blocks of global array into
!  message buffer and send to other processors. also copy local blocks
!
!-----------------------------------------------------------------------

   if (my_task == src_task) then

     !*** send non-local blocks away

     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) > 0 .and. &
           dst_dist%blockLocation(n)-1 /= my_task) then

         msg_buffer = 0
         this_block = get_block(n,n)

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
                      msg_buffer(i,j-yoffset2) &
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
                     msg_buffer(i,j) = ARRAY_G(isrc,jsrc)
               end do
            endif
         end do

         call MPI_SEND(msg_buffer, nx_block*ny_block, &
                       mpi_integer, dst_dist%blockLocation(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, ierr)

       endif
     end do

     deallocate(msg_buffer)

     !*** copy any local blocks

     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) == my_task+1) then
         dst_block = dst_dist%blockLocalID(n)
         this_block = get_block(n,n)

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
       endif
     end do

!-----------------------------------------------------------------------
!
!  otherwise receive data from src_task
!
!-----------------------------------------------------------------------

   else

     allocate (rcv_request(nblocks_tot), &
               rcv_status(MPI_STATUS_SIZE, nblocks_tot))

     rcv_request = 0
     rcv_status  = 0

     nrecvs = 0
     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) == my_task+1) then
         nrecvs = nrecvs + 1
         dst_block = dst_dist%blockLocalID(n)
         call MPI_IRECV(ARRAY(1,1,dst_block), nx_block*ny_block, &
                       mpi_integer, src_task, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, rcv_request(nrecvs), ierr)
       endif
     end do

     if (nrecvs > 0) &
       call MPI_WAITALL(nrecvs, rcv_request, rcv_status, ierr)

     deallocate(rcv_request, rcv_status)
   endif

   !-----------------------------------------------------------------
   ! Set ghost cell values to 0 for noupdate
   !-----------------------------------------------------------------

   if (field_loc == field_loc_noupdate) then
     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) == my_task+1 .and. &
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

   if (add_mpi_barriers) then
     call ice_barrier()
   endif

!-----------------------------------------------------------------------

 end subroutine scatter_global_int

!***********************************************************************

 subroutine scatter_global_ext(ARRAY, ARRAY_G, src_task, dst_dist)

!  This subroutine scatters a global-sized array to a distributed array.
!
!  This is the specific interface for double precision arrays
!  corresponding to the generic interface scatter_global.

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
     iblk, jblk,         &! block indices
     isrc, jsrc,         &! global indices
     nrecvs,             &! actual number of messages received
     dst_block,          &! location of block in dst array
     ierr                 ! MPI error flag

   type (block) :: &
     this_block  ! block info for current block

   integer (int_kind), dimension(:), allocatable :: &
     rcv_request     ! request array for receives

   integer (int_kind), dimension(:,:), allocatable :: &
     rcv_status      ! status array for receives

   real (dbl_kind), dimension(:,:), allocatable :: &
     msg_buffer      ! buffer for sending blocks

   character(len=*), parameter :: subname = '(scatter_global_ext)'

!-----------------------------------------------------------------------
!
!  initialize return array to zero
!
!-----------------------------------------------------------------------

   ARRAY = c0

!-----------------------------------------------------------------------
!
!  if this task is the src_task, copy blocks of global array into
!  message buffer and send to other processors. also copy local blocks
!
!-----------------------------------------------------------------------

   if (my_task == src_task) then

     !*** send non-local blocks away

     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) > 0 .and. &
           dst_dist%blockLocation(n)-1 /= my_task) then

         msg_buffer = c0
         this_block = get_block(n,n)

         ! interior
         do j = 1, ny_block
         do i = 1, nx_block
            if (this_block%j_glob(j) > 0) then
               msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i)+nghost,&
                                         this_block%j_glob(j)+nghost)
            endif
         end do
         end do

         if (this_block%jblock == 1) then
            ! south edge
            do j = 1, nghost
            do i = this_block%ilo,this_block%ihi
               msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i)+nghost,j)
            enddo
               do i = 1, nghost
               ! southwest corner
                  iblk = i
                  jblk = j
                  isrc = this_block%i_glob(this_block%ilo)+i-1
                  jsrc = j
                  msg_buffer(iblk,jblk) = ARRAY_G(isrc,jsrc)
               ! southeast corner
                  iblk = this_block%ihi+i
                  isrc = this_block%i_glob(this_block%ihi)+nghost+i
                  msg_buffer(iblk,jblk) = ARRAY_G(isrc,jsrc)
               enddo
            enddo
         endif
         if (this_block%jblock == nblocks_y) then
            ! north edge
            do j = 1, nghost
            do i = this_block%ilo,this_block%ihi
               msg_buffer(i,this_block%jhi+j) = ARRAY_G(this_block%i_glob(i)+nghost,&
                                                        ny_global+nghost+j)
            enddo
               do i = 1, nghost
               ! northwest corner
                  iblk = i
                  jblk = this_block%jhi+j
                  isrc = this_block%i_glob(this_block%ilo)+i-1
                  jsrc = ny_global+nghost+j
                  msg_buffer(iblk,jblk) = ARRAY_G(isrc,jsrc)
               ! northeast corner
                  iblk = this_block%ihi+i
                  isrc = this_block%i_glob(this_block%ihi)+nghost+i
                  msg_buffer(iblk,jblk) = ARRAY_G(isrc,jsrc)
               enddo
            enddo
         endif
         if (this_block%iblock == 1) then
            ! west edge
            do j = this_block%jlo,this_block%jhi
            do i = 1, nghost
               msg_buffer(i,j) = ARRAY_G(i,this_block%j_glob(j)+nghost)
            enddo
            enddo
               do j = 1, nghost
               do i = 1, nghost
               ! northwest corner
                  iblk = i
                  jblk = this_block%jhi+j
                  isrc = i
                  jsrc = this_block%j_glob(this_block%jhi)+nghost+j
                  msg_buffer(iblk,jblk) = ARRAY_G(isrc,jsrc)
               ! southwest corner
                  jblk = j
                  jsrc = this_block%j_glob(this_block%jlo)+j-1
                  msg_buffer(iblk,jblk) = ARRAY_G(isrc,jsrc)
               enddo
               enddo
         endif
         if (this_block%iblock == nblocks_x) then
            ! east edge
            do j = this_block%jlo,this_block%jhi
            do i = 1, nghost
               msg_buffer(this_block%ihi+i,j) = ARRAY_G(nx_global+nghost+i, &
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
                  msg_buffer(iblk,jblk) = ARRAY_G(isrc,jsrc)
               ! southeast corner
                  jblk = j
                  jsrc = this_block%j_glob(this_block%jlo)+j-1
                  msg_buffer(iblk,jblk) = ARRAY_G(isrc,jsrc)
               enddo
               enddo
         endif

         call MPI_SEND(msg_buffer, nx_block*ny_block, &
                       mpiR8, dst_dist%blockLocation(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, ierr)

       endif
     end do

     deallocate(msg_buffer)

     !*** copy any local blocks

     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) == my_task+1) then
         dst_block = dst_dist%blockLocalID(n)
         this_block = get_block(n,n)

         ! interior
         do j = 1, ny_block
         do i = 1, nx_block
            if (this_block%j_glob(j) > 0) then
               ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i)+nghost,&
                                              this_block%j_glob(j)+nghost)
            endif
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
       endif
     end do

!-----------------------------------------------------------------------
!
!  otherwise receive data from src_task
!
!-----------------------------------------------------------------------

   else

     allocate (rcv_request(nblocks_tot), &
               rcv_status(MPI_STATUS_SIZE, nblocks_tot))

     rcv_request = 0
     rcv_status  = 0

     nrecvs = 0
     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) == my_task+1) then
         nrecvs = nrecvs + 1
         dst_block = dst_dist%blockLocalID(n)
         call MPI_IRECV(ARRAY(1,1,dst_block), nx_block*ny_block, &
                       mpiR8, src_task, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, rcv_request(nrecvs), ierr)
       endif
     end do

     if (nrecvs > 0) &
       call MPI_WAITALL(nrecvs, rcv_request, rcv_status, ierr)

     deallocate(rcv_request, rcv_status)
   endif

   if (add_mpi_barriers) then
     call ice_barrier()
   endif

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
     nrecvs,             &! actual number of messages received
     isrc, jsrc,         &! source addresses
     dst_block,          &! location of block in dst array
     xoffset, yoffset,   &! offsets for tripole boundary conditions
     yoffset2,           &!
     isign,              &! sign factor for tripole boundary conditions
     ierr                 ! MPI error flag

   type (block) :: &
     this_block  ! block info for current block

   integer (int_kind), dimension(:), allocatable :: &
     rcv_request     ! request array for receives

   integer (int_kind), dimension(:,:), allocatable :: &
     rcv_status      ! status array for receives

   real (dbl_kind), dimension(:,:), allocatable :: &
     msg_buffer      ! buffer for sending blocks

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
!  if this task is the src_task, copy blocks of global array into
!  message buffer and send to other processors. also copy local blocks
!
!-----------------------------------------------------------------------

   if (my_task == src_task) then

     !*** send non-local blocks away

     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) > 0 .and. &
           dst_dist%blockLocation(n)-1 /= my_task) then

         msg_buffer = c0
         this_block = get_block(n,n)

         do j=1,ny_block
            if (this_block%jblock == nblocks_y .and. this_block%j_glob(j) < 0) then
               ! tripole is top block with j_glob < 0
               jsrc = ny_global + yoffset + &
                      (this_block%j_glob(j) + ny_global)
               do i=1,nx_block
                 if (this_block%i_glob(i) /= 0) then
                    isrc = nx_global + xoffset - this_block%i_glob(i)
                    if (isrc < 1) isrc = isrc + nx_global
                    if (isrc > nx_global) isrc = isrc - nx_global
                    msg_buffer(i,j) = isign * ARRAY_G2(isrc,jsrc)
                   endif
               end do
            else
               ! normal block
               do i=1,nx_block
                  isrc = this_block%i_glob(i)
                  jsrc = this_block%j_glob(j)
                  if (isrc >=1 .and. isrc <= nx_global .and. &
                      jsrc >=1 .and. jsrc <= ny_global) &
                     msg_buffer(i,j) = ARRAY_G1(isrc,jsrc)
               end do
            endif
         end do

         call MPI_SEND(msg_buffer, nx_block*ny_block, &
                       mpiR8, dst_dist%blockLocation(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, ierr)

       endif
     end do

     deallocate(msg_buffer)

     !*** copy any local blocks

     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) == my_task+1) then
         dst_block = dst_dist%blockLocalID(n)
         this_block = get_block(n,n)

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
       endif
     end do

!-----------------------------------------------------------------------
!
!  otherwise receive data from src_task
!
!-----------------------------------------------------------------------

   else

     allocate (rcv_request(nblocks_tot), &
               rcv_status(MPI_STATUS_SIZE, nblocks_tot))

     rcv_request = 0
     rcv_status  = 0

     nrecvs = 0
     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) == my_task+1) then
         nrecvs = nrecvs + 1
         dst_block = dst_dist%blockLocalID(n)
         call MPI_IRECV(ARRAY(1,1,dst_block), nx_block*ny_block, &
                       mpiR8, src_task, 3*mpitag_gs+n, &
                       MPI_COMM_ICE, rcv_request(nrecvs), ierr)
       endif
     end do

     if (nrecvs > 0) &
       call MPI_WAITALL(nrecvs, rcv_request, rcv_status, ierr)

     deallocate(rcv_request, rcv_status)
   endif

   if (add_mpi_barriers) then
     call ice_barrier()
   endif

!-----------------------------------------------------------------------

 end subroutine scatter_global_stress

!***********************************************************************

 end module ice_gather_scatter

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
