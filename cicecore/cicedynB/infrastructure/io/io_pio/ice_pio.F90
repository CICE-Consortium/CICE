!============================================================================
!  Writes netcdf files
!    Created by Mariana Vertenstein, June 2009

  module ice_pio

  use shr_kind_mod, only: r8 => shr_kind_r8, in=>shr_kind_in
  use shr_kind_mod, only: cl => shr_kind_cl
  use shr_sys_mod , only: shr_sys_flush
  use ice_kinds_mod
  use ice_blocks
  use ice_broadcast
  use ice_communicate
  use ice_domain, only : nblocks, blocks_ice
  use ice_domain_size
  use ice_fileunits  
  use ice_exit, only: abort_ice
  use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
  use pio

  implicit none
  private

  interface ice_pio_initdecomp
     module procedure ice_pio_initdecomp_2d
     module procedure ice_pio_initdecomp_3d
     module procedure ice_pio_initdecomp_4d
     module procedure ice_pio_initdecomp_3d_inner
  end interface

  public ice_pio_init
  public ice_pio_initdecomp

  type(iosystem_desc_t), pointer, public :: ice_pio_subsystem

!===============================================================================

  contains

!===============================================================================

!    Initialize the io subsystem
!    2009-Feb-17 - J. Edwards - initial version

   subroutine ice_pio_init(mode, filename, File, clobber, cdf64)

   use shr_pio_mod, only: shr_pio_getiosys, shr_pio_getiotype
     
   implicit none
   character(len=*)     , intent(in),    optional :: mode
   character(len=*)     , intent(in),    optional :: filename
   type(file_desc_t)    , intent(inout), optional :: File
   logical              , intent(in),    optional :: clobber
   logical              , intent(in),    optional :: cdf64

   ! local variables

   integer (int_kind) :: &
      nml_error          ! namelist read error flag

   integer :: pio_iotype
   logical :: exists
   logical :: lclobber
   logical :: lcdf64
   integer :: status
   integer :: nmode
   character(len=*), parameter :: subname = '(ice_pio_init)'
   logical, save :: first_call = .true.

   ice_pio_subsystem => shr_pio_getiosys(inst_name)
   pio_iotype =  shr_pio_getiotype(inst_name)

   if (present(mode) .and. present(filename) .and. present(File)) then
      
      if (trim(mode) == 'write') then
         lclobber = .false.
         if (present(clobber)) lclobber=clobber
         
         lcdf64 = .false.
         if (present(cdf64)) lcdf64=cdf64
         
         if (File%fh<0) then
            ! filename not open
            inquire(file=trim(filename),exist=exists)
            if (exists) then
               if (lclobber) then
                  nmode = pio_clobber
                  if (lcdf64) nmode = ior(nmode,PIO_64BIT_OFFSET)
                  status = pio_createfile(ice_pio_subsystem, File, pio_iotype, trim(filename), nmode)
                  if (my_task == master_task) then
                     write(nu_diag,*) subname,' create file ',trim(filename)
                  end if
               else
                  nmode = pio_write
                  status = pio_openfile(ice_pio_subsystem, File, pio_iotype, trim(filename), nmode)
                  if (my_task == master_task) then
                     write(nu_diag,*) subname,' open file ',trim(filename)
                  end if
               endif
            else
               nmode = pio_noclobber
               if (lcdf64) nmode = ior(nmode,PIO_64BIT_OFFSET)
               status = pio_createfile(ice_pio_subsystem, File, pio_iotype, trim(filename), nmode)
               if (my_task == master_task) then
                  write(nu_diag,*) subname,' create file ',trim(filename)
               end if
            endif
         else
            ! filename is already open, just return
         endif
      end if
      
      if (trim(mode) == 'read') then
         inquire(file=trim(filename),exist=exists)
         if (exists) then
            status = pio_openfile(ice_pio_subsystem, File, pio_iotype, trim(filename), pio_nowrite)
         else
            if(my_task==master_task) then
               write(nu_diag,*) 'ice_pio_ropen ERROR: file invalid ',trim(filename)
            end if
            call abort_ice(subname//'ERROR: aborting with invalid file')
         endif
      end if

   end if

   end subroutine ice_pio_init

!================================================================================

   subroutine ice_pio_initdecomp_2d(iodesc)

      type(io_desc_t), intent(out) :: iodesc

      integer (kind=int_kind) :: &
          iblk,ilo,ihi,jlo,jhi,lon,lat,i,j,n,k

      type(block) :: this_block 

      integer(kind=int_kind), pointer :: dof2d(:)
      character(len=*), parameter :: subname = '(ice_pio_initdecomp_2d)'

      allocate(dof2d(nx_block*ny_block*nblocks))

      n=0
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi
         
         do j=1,ny_block
         do i=1,nx_block
            n = n+1
            if (j < jlo .or. j>jhi) then
               dof2d(n) = 0
            else if (i < ilo .or. i > ihi) then
               dof2d(n) = 0
            else
               lon = this_block%i_glob(i)
               lat = this_block%j_glob(j)
               dof2d(n) = (lat-1)*nx_global + lon
            endif
         enddo !i
         enddo !j
      end do

      call pio_initdecomp(ice_pio_subsystem, pio_double, (/nx_global,ny_global/), &
           dof2d, iodesc)

      deallocate(dof2d)
 
   end subroutine ice_pio_initdecomp_2d

!================================================================================

   subroutine ice_pio_initdecomp_3d (ndim3, iodesc, remap)

      integer(kind=int_kind), intent(in) :: ndim3
      type(io_desc_t), intent(out) :: iodesc
      logical, optional :: remap
      integer (kind=int_kind) :: &
          iblk,ilo,ihi,jlo,jhi,lon,lat,i,j,n,k 

      type(block) :: this_block 
      logical :: lremap
      integer(kind=int_kind), pointer :: dof3d(:)
      character(len=*), parameter :: subname = '(ice_pio_initdecomp_2d)'

      allocate(dof3d(nx_block*ny_block*nblocks*ndim3))
      lremap=.false.
      if (present(remap)) lremap=remap
      if (lremap) then
         ! Reorder the ndim3 and nblocks loops to avoid a temporary array in restart read/write
         n=0
         do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)         
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            do k=1,ndim3         
               do j=1,ny_block
                  do i=1,nx_block
                     n = n+1
                     if (j < jlo .or. j>jhi) then
                        dof3d(n)=0
                     else if (i < ilo .or. i > ihi) then
                        dof3d(n) = 0
                     else
                        lon = this_block%i_glob(i)
                        lat = this_block%j_glob(j)
                        dof3d(n) = ((lat-1)*nx_global + lon) + (k-1)*nx_global*ny_global 
                     endif
                  enddo !i
               enddo !j
            enddo !ndim3
         enddo ! iblk
   else
         n=0
         do k=1,ndim3         
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)         
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j=1,ny_block
                  do i=1,nx_block
                     n = n+1
                     if (j < jlo .or. j>jhi) then
                        dof3d(n)=0
                     else if (i < ilo .or. i > ihi) then
                        dof3d(n) = 0
                     else
                        lon = this_block%i_glob(i)
                        lat = this_block%j_glob(j)
                        dof3d(n) = ((lat-1)*nx_global + lon) + (k-1)*nx_global*ny_global 
                     endif
                  enddo !i
               enddo !j
            enddo ! iblk
         enddo !ndim3
      endif

      call pio_initdecomp(ice_pio_subsystem, pio_double, (/nx_global,ny_global,ndim3/), &
           dof3d, iodesc)

      deallocate(dof3d)

   end subroutine ice_pio_initdecomp_3d

!================================================================================

   subroutine ice_pio_initdecomp_3d_inner(ndim3, inner_dim, iodesc)

      integer(kind=int_kind), intent(in) :: ndim3
      logical, intent(in) :: inner_dim
      type(io_desc_t), intent(out) :: iodesc

      integer (kind=int_kind) :: &
          iblk,ilo,ihi,jlo,jhi,lon,lat,i,j,n,k 

      type(block) :: this_block 

      integer(kind=int_kind), pointer :: dof3d(:)

      character(len=*), parameter :: subname = '(ice_pio_initdecomp_3d_inner)'

      allocate(dof3d(nx_block*ny_block*nblocks*ndim3))

      n=0
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi
         
         do j=1,ny_block
         do i=1,nx_block
         do k=1,ndim3
            n = n+1
            if (j < jlo .or. j>jhi) then
               dof3d(n) = 0
            else if (i < ilo .or. i > ihi) then
               dof3d(n) = 0
            else
               lon = this_block%i_glob(i)
               lat = this_block%j_glob(j)
               dof3d(n) = k + ((lon-1) + (lat-1)*nx_global)*ndim3
            endif
         end do !ndim3
         enddo  !i
         enddo  !j
      end do    !iblk

      call pio_initdecomp(ice_pio_subsystem, pio_double, (/ndim3,nx_global,ny_global/), &
           dof3d, iodesc)

      deallocate(dof3d)

   end subroutine ice_pio_initdecomp_3d_inner

!================================================================================

   subroutine ice_pio_initdecomp_4d (ndim3, ndim4, iodesc)

      integer(kind=int_kind), intent(in) :: ndim3, ndim4
      type(io_desc_t), intent(out) :: iodesc

      integer (kind=int_kind) :: &
          iblk,ilo,ihi,jlo,jhi,lon,lat,i,j,n,k,l 

      type(block) :: this_block 

      integer(kind=int_kind), pointer :: dof4d(:)

      character(len=*), parameter :: subname = '(ice_pio_initdecomp_4d)'

      allocate(dof4d(nx_block*ny_block*nblocks*ndim3*ndim4))

      n=0
      do l=1,ndim4
      do k=1,ndim3
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi
         
         do j=1,ny_block
         do i=1,nx_block
            n = n+1
            if (j < jlo .or. j>jhi) then
               dof4d(n)=0
            else if (i < ilo .or. i > ihi) then
               dof4d(n) = 0
            else
               lon = this_block%i_glob(i)
               lat = this_block%j_glob(j)
               dof4d(n) = ((lat-1)*nx_global + lon) &
                        + (k-1)*nx_global*ny_global & 
                        + (l-1)*nx_global*ny_global*ndim3 
            endif
         enddo !i
         enddo !j
      enddo ! iblk
      enddo !ndim3
      enddo !ndim4

      call pio_initdecomp(ice_pio_subsystem, pio_double, &
          (/nx_global,ny_global,ndim3,ndim4/), dof4d, iodesc)

      deallocate(dof4d)

   end subroutine ice_pio_initdecomp_4d
   
!================================================================================

  end module ice_pio

!================================================================================
