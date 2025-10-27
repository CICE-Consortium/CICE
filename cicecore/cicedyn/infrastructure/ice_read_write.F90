#ifdef ncdf
#define USE_NETCDF
#endif
!=======================================================================

! Routines for opening, reading and writing external files
!
! author: Tony Craig, NCAR
!
! 2004: Block structure added by William Lipscomb, LANL
! 2006: Converted to free source form (F90) by Elizabeth Hunke
! 2007: netcdf versions added by Alison McLaren & Ann Keen, Met Office

      module ice_read_write

      use,intrinsic :: ieee_arithmetic
      use ice_kinds_mod
      use ice_constants, only: c0, spval_dbl, &
          field_loc_noupdate, field_type_noupdate
      use ice_communicate, only: my_task, master_task
      use ice_broadcast, only: broadcast_scalar
      use ice_domain, only: distrb_info
      use ice_domain_size, only: max_blocks, nx_global, ny_global, ncat
      use ice_blocks, only: nx_block, ny_block, nghost
      use ice_exit, only: abort_ice
      use ice_fileunits, only: nu_diag

#ifdef USE_NETCDF
      use netcdf
#endif

      implicit none

      private

      integer (kind=int_kind), parameter, private :: &
         bits_per_byte = 8 ! number of bits per byte.
                           ! used to determine RecSize in ice_open

      public :: ice_open,           &
                ice_open_ext,       &
                ice_open_nc,        &
                ice_read,           &
                ice_read_ext,       &
                ice_read_nc,        &
                ice_read_global,    &
                ice_read_global_nc, &
                ice_read_nc_uv,     &
                ice_read_nc_xyf,    &
                ice_write,          &
                ice_write_nc,       &
                ice_write_ext,      &
                ice_read_vec_nc,    &
                ice_get_ncvarsize,  &
                ice_check_nc,       &
                ice_close_nc

      interface ice_write
         module procedure ice_write_xyt,  &
                          ice_write_xyzt
      end interface

      interface ice_read
         module procedure ice_read_xyt,  &
                          ice_read_xyzt
      end interface

      interface ice_read_nc
         module procedure ice_read_nc_xy,  &
                          ice_read_nc_xyz, &
                          !ice_read_nc_xyf, &
                          ice_read_nc_point, &
                          ice_read_nc_1D,  &
                          ice_read_nc_2D,  &
                          ice_read_nc_3D,  &
                          ice_read_nc_z
      end interface

      interface ice_write_nc
         module procedure ice_write_nc_xy,  &
                          ice_write_nc_xyz
      end interface

!=======================================================================

      contains

!=======================================================================

! Opens an unformatted file for reading.
! nbits indicates whether the file is sequential or direct access.
!
! author: Tony Craig, NCAR

      subroutine ice_open(nu, filename, nbits, algn)

      integer (kind=int_kind), intent(in) :: &
         nu        , & ! unit number
         nbits         ! no. of bits per variable (0 for sequential access)

      integer (kind=int_kind), intent(in), optional :: algn
      integer (kind=int_kind) :: RecSize, Remnant, nbytes

      character (*) :: filename

      character(len=*), parameter :: subname = '(ice_open)'

      if (my_task == master_task) then

         if (nbits == 0) then   ! sequential access

            open(nu,file=filename,form='unformatted')

         else                   ! direct access

            ! use nbytes to compute RecSize.
            ! this prevents integer overflow with large global grids using nbits
            ! nx*ny*nbits > 2^31 -1 (i.e., global grid 9000x7054x64)
            nbytes = nbits/bits_per_byte
            RecSize = nx_global*ny_global*nbytes

            if (present(algn)) then
              ! If data is keept in blocks using given sizes (=algn)
              !  Used in eg. HYCOM binary files, which are stored as "blocks" dividable by 16384 bit (=algn)
              if (algn /= 0) then
                Remnant = modulo(RecSize,algn)
                if (Remnant /= 0) then
                  RecSize = RecSize + (algn - Remnant)
                endif
              endif
            endif
            open(nu,file=filename,recl=RecSize, &
                  form='unformatted',access='direct')
         endif                   ! nbits = 0

      endif                      ! my_task = master_task

      end subroutine ice_open

!=======================================================================

! Opens an unformatted file for reading, incl ghost cells (direct access).
! nbits indicates whether the file is sequential or direct access.
!
! authors: Tony Craig, NCAR
!          David Hebert, NRLSSC

      subroutine ice_open_ext(nu, filename, nbits)

      integer (kind=int_kind), intent(in) :: &
         nu        , & ! unit number
         nbits         ! no. of bits per variable (0 for sequential access)

      integer (kind=int_kind) :: RecSize, nbytes

      character (*) :: filename

      integer (kind=int_kind) :: &
         nx, ny        ! grid dimensions including ghost cells

      character(len=*), parameter :: subname = '(ice_open_ext)'

      if (my_task == master_task) then

         if (nbits == 0) then   ! sequential access

            open(nu,file=filename,form='unformatted')

         else                   ! direct access

            nx = nx_global + 2*nghost
            ny = ny_global + 2*nghost

            ! use nbytes to compute RecSize.
            ! this prevents integer overflow with large global grids using nbits
            ! nx*ny*nbits > 2^31 -1 (i.e., global grid 9000x7054x64)
            nbytes = nbits/bits_per_byte
            RecSize = nx*ny*nbytes
            open(nu,file=filename,recl=RecSize, &
                  form='unformatted',access='direct')
         endif                   ! nbits = 0

      endif                      ! my_task = master_task

      end subroutine ice_open_ext

!=======================================================================

! Read an unformatted file and scatter to processors.
! work is a real array, atype indicates the format of the data.
! If the optional variables field_loc and field_type are present,
! the ghost cells are filled using values from the global array.
! This prevents them from being filled with zeroes in land cells
! (subroutine ice_HaloUpdate need not be called).
!
! author: Tony Craig, NCAR

      subroutine ice_read_xyt(nu, nrec, work, atype, diag, &
                          field_loc, field_type, &
                          ignore_eof, hit_eof)

      use ice_gather_scatter, only: scatter_global

      integer (kind=int_kind), intent(in) :: &
         nu            , & ! unit number
         nrec              ! record number (0 for sequential access)

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks), intent(out) :: &
         work              ! output array (real, 8-byte)

      character (len=4), intent(in) :: &
         atype             ! format for input array
                           ! (real/integer, 4-byte/8-byte)

      logical (kind=log_kind), intent(in) :: &
         diag              ! if true, write diagnostic output

      integer (kind=int_kind), optional, intent(in) :: &
         field_loc, &      ! location of field on staggered grid
         field_type        ! type of field (scalar, vector, angle)

      logical (kind=log_kind), optional, intent(in)  :: ignore_eof
      logical (kind=log_kind), optional, intent(out) :: hit_eof

      ! local variables

      integer (kind=int_kind) :: i, j, ios

      real (kind=dbl_kind) :: &
         amin, amax, asum  ! min, max values and sum of input array

      logical (kind=log_kind) :: ignore_eof_use

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g1

      real (kind=real_kind), dimension(:,:), allocatable :: &
         work_gr

      integer(kind=int_kind), dimension(:,:), allocatable :: &
         work_gi4

      integer(selected_int_kind(13)), dimension(:,:), allocatable :: &
         work_gi8

      character(len=*), parameter :: subname = '(ice_read_xyt)'

      if (my_task == master_task) then
         allocate(work_g1(nx_global,ny_global))
      else
         allocate(work_g1(1,1))   ! to save memory
      endif

      if (my_task == master_task) then

         !-------------------------------------------------------------------
         ! Read global array according to format atype
         !-------------------------------------------------------------------

         if (present(hit_eof)) hit_eof = .false.

         if (atype == 'ida4') then
            allocate(work_gi4(nx_global,ny_global))
            read(nu,rec=nrec) work_gi4
            work_g1 = real(work_gi4,kind=dbl_kind)
            deallocate(work_gi4)
         elseif (atype == 'ida8') then
            allocate(work_gi8(nx_global,ny_global))
            read(nu,rec=nrec) work_gi8
            work_g1 = real(work_gi8,kind=dbl_kind)
            deallocate(work_gi8)
         elseif (atype == 'rda4') then
            allocate(work_gr(nx_global,ny_global))
            read(nu,rec=nrec) work_gr
            work_g1 = work_gr
            deallocate(work_gr)
         elseif (atype == 'rda8') then
            read(nu,rec=nrec) work_g1
         elseif (atype == 'ruf8') then
            if (present(ignore_eof)) then
               ignore_eof_use = ignore_eof
            else
               ignore_eof_use = .false.
            endif
            if (ignore_eof_use) then
               ! Read line from file, checking for end-of-file
               read(nu, iostat=ios) ((work_g1(i,j),i=1,nx_global), &
                                                   j=1,ny_global)
               if (present(hit_eof)) hit_eof = ios < 0
            else
               read(nu) ((work_g1(i,j),i=1,nx_global),j=1,ny_global)
            endif
         else
            write(nu_diag,*) subname,' ERROR: reading unknown atype ',atype
         endif
      endif                     ! my_task = master_task

      if (present(hit_eof)) then
         call broadcast_scalar(hit_eof,master_task)
         if (hit_eof) then
            deallocate(work_g1)
            return
         endif
      endif

      !-------------------------------------------------------------------
      ! optional diagnostics
      !-------------------------------------------------------------------

      if (my_task==master_task .and. diag) then
         amin = minval(work_g1)
         amax = maxval(work_g1, mask = work_g1 /= spval_dbl)
         asum = sum(work_g1, mask = work_g1 /= spval_dbl)
         write(nu_diag,*) subname,' read_global ',nu, nrec, amin, amax, asum
      endif

      !-------------------------------------------------------------------
      ! Scatter data to individual processors.
      ! NOTE: Ghost cells are not updated unless field_loc is present.
      !-------------------------------------------------------------------

      if (present(field_loc)) then
         call scatter_global(work, work_g1, master_task, distrb_info, &
                             field_loc, field_type)
      else
         call scatter_global(work, work_g1, master_task, distrb_info, &
                             field_loc_noupdate, field_type_noupdate)
      endif

      deallocate(work_g1)

      end subroutine ice_read_xyt

!=======================================================================
! Read an unformatted file and scatter to processors.
! work is a real array, atype indicates the format of the data.
! If the optional variables field_loc and field_type are present,
! the ghost cells are filled using values from the global array.
! This prevents them from being filled with zeroes in land cells
! (subroutine ice_HaloUpdate need not be called).
!
! author: Tony Craig, NCAR

      subroutine ice_read_xyzt(nu,  nrec,  work, atype, diag, &
                          field_loc, field_type, &
                          ignore_eof, hit_eof)

      use ice_gather_scatter, only: scatter_global
      use ice_domain_size, only: nblyr

      integer (kind=int_kind), intent(in) :: &
         nu            , & ! unit number
         nrec              ! record number (0 for sequential access)

      real (kind=dbl_kind), dimension(nx_block,ny_block,nblyr+2,max_blocks), intent(out) :: &
         work              ! output array (real, 8-byte)

      character (len=4), intent(in) :: &
         atype             ! format for input array
                           ! (real/integer, 4-byte/8-byte)

      logical (kind=log_kind), intent(in) :: &
         diag              ! if true, write diagnostic output

      integer (kind=int_kind), optional, intent(in) :: &
         field_loc, &      ! location of field on staggered grid
         field_type        ! type of field (scalar, vector, angle)

      logical (kind=log_kind), optional, intent(in)  :: ignore_eof
      logical (kind=log_kind), optional, intent(out) :: hit_eof

      ! local variables

      integer (kind=int_kind) :: i, j, k, ios

      real (kind=dbl_kind) :: &
         amin, amax, asum  ! min, max values and sum of input array

      logical (kind=log_kind) :: ignore_eof_use


      real (kind=dbl_kind), dimension(:,:,:), allocatable :: &
         work_g4

      integer(kind=int_kind), dimension(:,:,:), allocatable :: &
         work_gi5

      integer(selected_int_kind(13)), dimension(:,:,:), allocatable :: &
         work_gi9

      real (kind=real_kind), dimension(:,:,:), allocatable :: &
         work_gr3

      character(len=*), parameter :: subname = '(ice_read_xyzt)'

      if (my_task == master_task) then
         allocate(work_g4(nx_global,ny_global,nblyr+2))
      else
         allocate(work_g4(1,1,nblyr+2))   ! to save memory
      endif

      if (my_task == master_task) then

         !-------------------------------------------------------------------
         ! Read global array according to format atype
         !-------------------------------------------------------------------

         if (present(hit_eof)) hit_eof = .false.

         if (atype == 'ida4') then
            allocate(work_gi5(nx_global,ny_global,nblyr+2))
            read(nu,rec=nrec) work_gi5
            work_g4 = real(work_gi5,kind=dbl_kind)
            deallocate(work_gi5)
         elseif (atype == 'ida8') then
            allocate(work_gi9(nx_global,ny_global,nblyr+2))
            read(nu,rec=nrec) work_gi9
            work_g4 = real(work_gi9,kind=dbl_kind)
            deallocate(work_gi9)
         elseif (atype == 'rda4') then
            allocate(work_gr3(nx_global,ny_global,nblyr+2))
            read(nu,rec=nrec) work_gr3
            work_g4 = work_gr3
            deallocate(work_gr3)
         elseif (atype == 'rda8') then
            read(nu,rec=nrec) work_g4
         elseif (atype == 'ruf8') then
            if (present(ignore_eof)) then
               ignore_eof_use = ignore_eof
            else
               ignore_eof_use = .false.
            endif
            if (ignore_eof_use) then
               ! Read line from file, checking for end-of-file
               read(nu, iostat=ios) (((work_g4(i,j,k),i=1,nx_global), &
                                                      j=1,ny_global), &
                                                      k=1,nblyr+2)
               if (present(hit_eof)) hit_eof = ios < 0
            else
               read(nu) (((work_g4(i,j,k),i=1,nx_global),j=1,ny_global),&
                                                         k=1,nblyr+2)
            endif
         else
            write(nu_diag,*) subname,' ERROR: reading unknown atype ',atype
         endif
      endif                     ! my_task = master_task

      if (present(hit_eof)) then
         call broadcast_scalar(hit_eof,master_task)
         if (hit_eof) then
            deallocate(work_g4)
            return
         endif
      endif

      !-------------------------------------------------------------------
      ! optional diagnostics
      !-------------------------------------------------------------------

      if (my_task==master_task .and. diag) then
         amin = minval(work_g4)
         amax = maxval(work_g4, mask = work_g4 /= spval_dbl)
         asum = sum   (work_g4, mask = work_g4 /= spval_dbl)
         write(nu_diag,*) subname,' read_global ',nu, nrec, amin, amax, asum
      endif

      !-------------------------------------------------------------------
      ! Scatter data to individual processors.
      ! NOTE: Ghost cells are not updated unless field_loc is present.
      !-------------------------------------------------------------------

      do k = 1, nblyr+2

         if (present(field_loc)) then
            call scatter_global(work(:,:,k,:), work_g4(:,:,k), master_task, distrb_info, &
                                field_loc, field_type)

         else

            call scatter_global(work(:,:,k,:), work_g4(:,:,k), master_task, distrb_info, &
                                field_loc_noupdate, field_type_noupdate)
         endif

      enddo   !k
      deallocate(work_g4)

      end subroutine ice_read_xyzt

!=======================================================================

! Read an unformatted file
! Just like ice_read except that it returns a global array.
! work_g is a real array, atype indicates the format of the data
!
! Adapted by William Lipscomb, LANL, from ice_read

      subroutine ice_read_global (nu,  nrec,  work_g, atype, diag, &
                                  ignore_eof, hit_eof)

      integer (kind=int_kind), intent(in) :: &
         nu            , & ! unit number
         nrec              ! record number (0 for sequential access)

      real (kind=dbl_kind), dimension(nx_global,ny_global), intent(out) :: &
         work_g            ! output array (real, 8-byte)

      character (len=4) :: &
         atype             ! format for input array
                           ! (real/integer, 4-byte/8-byte)

      logical (kind=log_kind) :: &
         diag              ! if true, write diagnostic output

      logical (kind=log_kind), optional, intent(in)  :: ignore_eof
      logical (kind=log_kind), optional, intent(out) :: hit_eof

      ! local variables

      integer (kind=int_kind) :: i, j, ios

      real (kind=dbl_kind) :: &
         amin, amax, asum  ! min, max values and sum of input array

      logical (kind=log_kind) :: ignore_eof_use

      real (kind=real_kind), dimension(:,:), allocatable :: &
         work_gr

      integer(kind=int_kind), dimension(:,:), allocatable :: &
         work_gi4

      integer(selected_int_kind(13)), dimension(:,:), allocatable :: &
         work_gi8

      character(len=*), parameter :: subname = '(ice_read_global)'

      work_g(:,:) = c0

      if (my_task == master_task) then

         !-------------------------------------------------------------------
         ! Read global array according to format atype
         !-------------------------------------------------------------------

         if (present(hit_eof)) hit_eof = .false.

         if (atype == 'ida4') then
            allocate(work_gi4(nx_global,ny_global))
            read(nu,rec=nrec) work_gi4
            work_g = real(work_gi4,kind=dbl_kind)
            deallocate(work_gi4)
         elseif (atype == 'ida8') then
            allocate(work_gi8(nx_global,ny_global))
            read(nu,rec=nrec) work_gi8
            work_g = real(work_gi8,kind=dbl_kind)
            deallocate(work_gi8)
         elseif (atype == 'rda4') then
            allocate(work_gr(nx_global,ny_global))
            read(nu,rec=nrec) work_gr
            work_g = work_gr
            deallocate(work_gr)
         elseif (atype == 'rda8') then
            read(nu,rec=nrec) work_g
         elseif (atype == 'ruf8') then
            if (present(ignore_eof)) then
               ignore_eof_use = ignore_eof
            else
               ignore_eof_use = .false.
            endif
            if (ignore_eof_use) then
               ! Read line from file, checking for end-of-file
               read(nu, iostat=ios) ((work_g(i,j),i=1,nx_global), &
                                                  j=1,ny_global)
               if (present(hit_eof)) hit_eof = ios < 0
            else
               read(nu) ((work_g(i,j),i=1,nx_global),j=1,ny_global)
            endif
         else
            write(nu_diag,*) subname,' ERROR: reading unknown atype ',atype
         endif
      endif                     ! my_task = master_task

      if (present(hit_eof)) then
         call broadcast_scalar(hit_eof,master_task)
         if (hit_eof) return
      endif

      !-------------------------------------------------------------------
      ! optional diagnostics
      !-------------------------------------------------------------------

      if (my_task == master_task .and. diag) then
         amin = minval(work_g)
         amax = maxval(work_g, mask = work_g /= spval_dbl)
         asum = sum   (work_g, mask = work_g /= spval_dbl)
         write(nu_diag,*) subname,' read_global ',nu, nrec, amin, amax,asum
      endif

      end subroutine ice_read_global

!=======================================================================

! Read an unformatted file and scatter to processors, incl ghost cells.
! work is a real array, atype indicates the format of the data.
! (subroutine ice_HaloUpdate need not be called).

      subroutine ice_read_ext(nu,  nrec,  work, atype, diag, &
                          ignore_eof, hit_eof)

      use ice_gather_scatter, only: scatter_global_ext

      integer (kind=int_kind), intent(in) :: &
         nu            , & ! unit number
         nrec              ! record number (0 for sequential access)

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks), intent(out) :: &
         work              ! output array (real, 8-byte)

      character (len=4), intent(in) :: &
         atype             ! format for input array
                           ! (real/integer, 4-byte/8-byte)

      logical (kind=log_kind), intent(in) :: &
         diag              ! if true, write diagnostic output

      logical (kind=log_kind), optional, intent(in)  :: ignore_eof
      logical (kind=log_kind), optional, intent(out) :: hit_eof

      ! local variables

      integer (kind=int_kind) :: i, j, ios, nx, ny

      real (kind=dbl_kind) :: &
         amin, amax, asum  ! min, max values and sum of input array

      logical (kind=log_kind) :: ignore_eof_use

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g1

      real (kind=real_kind), dimension(:,:), allocatable :: &
         work_gr

      integer(kind=int_kind), dimension(:,:), allocatable :: &
         work_gi4

      integer(selected_int_kind(13)), dimension(:,:), allocatable :: &
         work_gi8

      character(len=*), parameter :: subname = '(ice_read_ext)'

      nx = nx_global + 2*nghost
      ny = ny_global + 2*nghost

      if (my_task == master_task) then
         allocate(work_g1(nx,ny))
      else
         allocate(work_g1(1,1))   ! to save memory
      endif

      if (my_task == master_task) then

         !-------------------------------------------------------------------
         ! Read global array according to format atype
         !-------------------------------------------------------------------

         if (present(hit_eof)) hit_eof = .false.

         if (atype == 'ida4') then
            allocate(work_gi4(nx,ny))
            read(nu,rec=nrec) work_gi4
            work_g1 = real(work_gi4,kind=dbl_kind)
            deallocate(work_gi4)
         elseif (atype == 'ida8') then
            allocate(work_gi8(nx,ny))
            read(nu,rec=nrec) work_gi8
            work_g1 = real(work_gi8,kind=dbl_kind)
            deallocate(work_gi8)
         elseif (atype == 'rda4') then
            allocate(work_gr(nx,ny))
            read(nu,rec=nrec) work_gr
            work_g1 = work_gr
            deallocate(work_gr)
         elseif (atype == 'rda8') then
            read(nu,rec=nrec) work_g1
         elseif (atype == 'ruf8') then
            if (present(ignore_eof)) then
               ignore_eof_use = ignore_eof
            else
               ignore_eof_use = .false.
            endif
            if (ignore_eof_use) then
               ! Read line from file, checking for end-of-file
               read(nu, iostat=ios) ((work_g1(i,j),i=1,nx), &
                                                   j=1,ny)
               if (present(hit_eof)) hit_eof = ios < 0
            else
               read(nu) ((work_g1(i,j),i=1,nx),j=1,ny)
            endif
         else
            write(nu_diag,*) subname,' ERROR: reading unknown atype ',atype
         endif
      endif                     ! my_task = master_task

      if (present(hit_eof)) then
         call broadcast_scalar(hit_eof,master_task)
         if (hit_eof) then
            deallocate(work_g1)
            return
         endif
      endif

      !-------------------------------------------------------------------
      ! optional diagnostics
      !-------------------------------------------------------------------

      if (my_task==master_task .and. diag) then
         amin = minval(work_g1)
         amax = maxval(work_g1, mask = work_g1 /= spval_dbl)
         asum = sum   (work_g1, mask = work_g1 /= spval_dbl)
         write(nu_diag,*) subname,' read_global ',nu, nrec, amin, amax, asum
      endif

      !-------------------------------------------------------------------
      ! Scatter data to individual processors.
      ! NOTE: Ghost cells are always updated
      !-------------------------------------------------------------------

      call scatter_global_ext(work, work_g1, master_task, distrb_info)

      deallocate(work_g1)

      end subroutine ice_read_ext

!=======================================================================

! Writes an unformatted file
! work is a real array, atype indicates the format of the data

      subroutine ice_write_xyt(nu, nrec, work, atype, diag)

      use ice_gather_scatter, only: gather_global

      integer (kind=int_kind), intent(in) :: &
         nu            , & ! unit number
         nrec              ! record number (0 for sequential access)

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks), intent(in) :: &
         work              ! input array (real, 8-byte)

      character (len=4), intent(in) :: &
         atype             ! format for output array
                           ! (real/integer, 4-byte/8-byte)

      logical (kind=log_kind), intent(in) :: &
         diag              ! if true, write diagnostic output

      ! local variables

      integer (kind=int_kind) :: i, j

      real (kind=dbl_kind) :: &
         amin, amax, asum  ! min, max values and sum of input array

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g1

      real (kind=real_kind), dimension(:,:), allocatable :: &
         work_gr

      integer(kind=int_kind), dimension(:,:), allocatable :: &
         work_gi4

      integer(selected_int_kind(13)), dimension(:,:), allocatable :: &
         work_gi8

      character(len=*), parameter :: subname = '(ice_write_xyt)'

      !-------------------------------------------------------------------
      ! Gather data from individual processors
      !-------------------------------------------------------------------

      if (my_task == master_task) then
         allocate(work_g1(nx_global,ny_global))
      else
         allocate(work_g1(1,1)) ! to save memory
      endif

      call gather_global(work_g1, work, master_task, distrb_info, spc_val=c0)

      if (my_task == master_task) then

         !-------------------------------------------------------------------
         ! Write global array according to format atype
         !-------------------------------------------------------------------

         if (atype == 'ida4') then
            allocate(work_gi4(nx_global,ny_global))
            work_gi4 = nint(work_g1)
            write(nu,rec=nrec) work_gi4
            deallocate(work_gi4)
         elseif (atype == 'ida8') then
            allocate(work_gi8(nx_global,ny_global))
            work_gi8 = nint(work_g1)
            write(nu,rec=nrec) work_gi8
            deallocate(work_gi8)
         elseif (atype == 'rda4') then
            allocate(work_gr(nx_global,ny_global))
            work_gr = real(work_g1,real_kind)
            write(nu,rec=nrec) work_gr
            deallocate(work_gr)
         elseif (atype == 'rda8') then
            write(nu,rec=nrec) work_g1
         elseif (atype == 'ruf8') then
            write(nu) ((work_g1(i,j),i=1,nx_global),j=1,ny_global)
         else
            write(nu_diag,*) subname,' ERROR: writing unknown atype ',atype
         endif

         !-------------------------------------------------------------------
         ! diagnostics
         !-------------------------------------------------------------------

         if (diag) then
            amin = minval(work_g1)
            amax = maxval(work_g1, mask = work_g1 /= spval_dbl)
            asum = sum   (work_g1, mask = work_g1 /= spval_dbl)
            write(nu_diag,*) subname,' write_global ', nu, nrec, amin, amax, asum
         endif

      endif                     ! my_task = master_task

      deallocate(work_g1)

      end subroutine ice_write_xyt

!=======================================================================

! Writes an unformatted file
! work is a real array, atype indicates the format of the data

      subroutine ice_write_xyzt(nu, nrec, work, atype, diag)

      use ice_gather_scatter, only: gather_global
      use ice_domain_size, only: nblyr

      integer (kind=int_kind), intent(in) :: &
         nu            , & ! unit number
         nrec              ! record number (0 for sequential access)

      real (kind=dbl_kind), dimension(nx_block,ny_block,nblyr+2,max_blocks), intent(in) :: &
         work              ! input array (real, 8-byte)

      character (len=4), intent(in) :: &
         atype             ! format for output array
                           ! (real/integer, 4-byte/8-byte)

      logical (kind=log_kind), intent(in) :: &
         diag              ! if true, write diagnostic output

      ! local variables

      integer (kind=int_kind) :: i, j, k

      real (kind=dbl_kind) :: &
         amin, amax, asum  ! min, max values and sum of input array

      real (kind=dbl_kind), dimension(:,:,:), allocatable :: &
         work_g4

      real (kind=real_kind), dimension(:,:,:), allocatable :: &
         work_gr3

      integer(kind=int_kind), dimension(:,:,:), allocatable :: &
         work_gi5

      integer(selected_int_kind(13)), dimension(:,:,:), allocatable :: &
         work_gi9

      character(len=*), parameter :: subname = '(ice_write_xyzt)'

      !-------------------------------------------------------------------
      ! Gather data from individual processors
      !-------------------------------------------------------------------

      if (my_task == master_task) then
         allocate(work_g4(nx_global,ny_global,nblyr+2))
      else
         allocate(work_g4(1,1,nblyr+2)) ! to save memory
      endif
      do k = 1,nblyr+2
         call gather_global(work_g4(:,:,k), work(:,:,k,:), master_task, &
                            distrb_info, spc_val=c0)
      enddo   !k

      if (my_task == master_task) then

         !-------------------------------------------------------------------
         ! Write global array according to format atype
         !-------------------------------------------------------------------

         if (atype == 'ida4') then
            allocate(work_gi5(nx_global,ny_global,nblyr+2))
            work_gi5 = nint(work_g4)
            write(nu,rec=nrec) work_gi5
            deallocate(work_gi5)
         elseif (atype == 'ida8') then
            allocate(work_gi9(nx_global,ny_global,nblyr+2))
            work_gi9 = nint(work_g4)
            write(nu,rec=nrec) work_gi9
            deallocate(work_gi9)
         elseif (atype == 'rda4') then
            allocate(work_gr3(nx_global,ny_global,nblyr+2))
            work_gr3 = real(work_g4,real_kind)
            write(nu,rec=nrec) work_gr3
            deallocate(work_gr3)
         elseif (atype == 'rda8') then
            write(nu,rec=nrec) work_g4
         elseif (atype == 'ruf8') then
            write(nu)(((work_g4(i,j,k),i=1,nx_global),j=1,ny_global), &
                                       k=1,nblyr+2)
         else
            write(nu_diag,*) subname,' ERROR: writing unknown atype ',atype
         endif

         !-------------------------------------------------------------------
         ! diagnostics
         !-------------------------------------------------------------------

         if (diag) then
            amin = minval(work_g4)
            amax = maxval(work_g4, mask = work_g4 /= spval_dbl)
            asum = sum   (work_g4, mask = work_g4 /= spval_dbl)
            write(nu_diag,*) subname,' write_global ', nu, nrec, amin, amax, asum
         endif

      endif                     ! my_task = master_task

      deallocate(work_g4)

      end subroutine ice_write_xyzt

!=======================================================================
!
! Writes an unformatted file, including ghost cells
! work is a real array, atype indicates the format of the data
!
! author: Tony Craig, NCAR

      subroutine ice_write_ext(nu, nrec, work, atype, diag)

      use ice_gather_scatter, only: gather_global_ext

      integer (kind=int_kind), intent(in) :: &
         nu            , & ! unit number
         nrec              ! record number (0 for sequential access)

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks), intent(in) :: &
         work              ! input array (real, 8-byte)

      character (len=4), intent(in) :: &
         atype             ! format for output array
                           ! (real/integer, 4-byte/8-byte)

      logical (kind=log_kind), intent(in) :: &
         diag              ! if true, write diagnostic output

      ! local variables

      integer (kind=int_kind) :: i, j, nx, ny

      real (kind=dbl_kind) :: &
         amin, amax, asum  ! min, max values and sum of input array

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g1

      real (kind=real_kind), dimension(:,:), allocatable :: &
         work_gr

      integer(kind=int_kind), dimension(:,:), allocatable :: &
         work_gi4

      integer(selected_int_kind(13)), dimension(:,:), allocatable :: &
         work_gi8

      character(len=*), parameter :: subname = '(ice_write_ext)'

      !-------------------------------------------------------------------
      ! Gather data from individual processors
      !-------------------------------------------------------------------

      nx = nx_global + 2*nghost
      ny = ny_global + 2*nghost

      if (my_task == master_task) then
         allocate(work_g1(nx,ny))
      else
         allocate(work_g1(1,1)) ! to save memory
      endif

      call gather_global_ext(work_g1, work, master_task, distrb_info, spc_val=c0)

      if (my_task == master_task) then

         !-------------------------------------------------------------------
         ! Write global array according to format atype
         !-------------------------------------------------------------------

         if (atype == 'ida4') then
            allocate(work_gi4(nx,ny))
            work_gi4 = nint(work_g1)
            write(nu,rec=nrec) work_gi4
            deallocate(work_gi4)
         elseif (atype == 'ida8') then
            allocate(work_gi8(nx,ny))
            work_gi8 = nint(work_g1)
            write(nu,rec=nrec) work_gi8
            deallocate(work_gi8)
         elseif (atype == 'rda4') then
            allocate(work_gr(nx,ny))
            work_gr = real(work_g1,real_kind)
            write(nu,rec=nrec) work_gr
            deallocate(work_gr)
         elseif (atype == 'rda8') then
            write(nu,rec=nrec) work_g1
         elseif (atype == 'ruf8') then
            write(nu) ((work_g1(i,j),i=1,nx),j=1,ny)
         else
            write(nu_diag,*) subname,' ERROR: writing unknown atype ',atype
         endif

         !-------------------------------------------------------------------
         ! diagnostics
         !-------------------------------------------------------------------

         if (diag) then
            amin = minval(work_g1)
            amax = maxval(work_g1, mask = work_g1 /= spval_dbl)
            asum = sum   (work_g1, mask = work_g1 /= spval_dbl)
            write(nu_diag,*) subname,' write_global ', nu, nrec, amin, amax, asum
         endif

      endif                     ! my_task = master_task

      deallocate(work_g1)

      end subroutine ice_write_ext

!=======================================================================

! Opens a netCDF file for reading
! Adapted by Alison McLaren, Met Office from ice_open

      subroutine ice_open_nc(filename, fid)

      character (char_len_long), intent(in) :: &
         filename      ! netCDF filename

      integer (kind=int_kind), intent(out) :: &
         fid           ! unit number

      ! local variables

      character(len=*), parameter :: subname = '(ice_open_nc)'

#ifdef USE_NETCDF
      integer (kind=int_kind) :: &
         status        ! status variable from netCDF routine

      if (my_task == master_task) then

          status = nf90_open(filename, NF90_NOWRITE, fid)
          call ice_check_nc(status, subname//' ERROR: Cannot open '//trim(filename), &
                            file=__FILE__, line=__LINE__)

      endif                      ! my_task = master_task

#else
      call abort_ice(subname//' ERROR: USE_NETCDF cpp not defined for '//trim(filename), &
         file=__FILE__, line=__LINE__)
      fid = -999 ! to satisfy intent(out) attribute
#endif
      end subroutine ice_open_nc

!=======================================================================

! Read a netCDF file and scatter to processors.
! If the optional variables field_loc and field_type are present,
! the ghost cells are filled using values from the global array.
! This prevents them from being filled with zeroes in land cells
! (subroutine ice_HaloUpdate need not be called).
!
! Adapted by Alison McLaren, Met Office from ice_read

      subroutine ice_read_nc_xy(fid,  nrec,  varname, work,  diag, &
                             field_loc, field_type, restart_ext)

      use ice_gather_scatter, only: scatter_global, scatter_global_ext

      integer (kind=int_kind), intent(in) :: &
         fid           , & ! file id
         nrec              ! record number

      logical (kind=log_kind), intent(in) :: &
         diag              ! if true, write diagnostic output

      character (len=*), intent(in) :: &
         varname           ! field name in netcdf file

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks), intent(out) :: &
         work              ! output array (real, 8-byte)

      logical (kind=log_kind), optional, intent(in) :: &
         restart_ext       ! if true, read extended grid

      integer (kind=int_kind), optional, intent(in) :: &
         field_loc, &      ! location of field on staggered grid
         field_type        ! type of field (scalar, vector, angle)

      ! local variables

      character(len=*), parameter :: subname = '(ice_read_nc_xy)'

#ifdef USE_NETCDF
! netCDF file diagnostics:
      integer (kind=int_kind) :: &
         varid         , & ! variable id
         status        , & ! status output from netcdf routines
         ndims         , & ! number of dimensions
         dimlen            ! dimension size

      integer (kind=int_kind), dimension(10) :: &
         dimids            ! generic size dimids

      real (kind=dbl_kind) :: &
         missingvalue, &
         amin, amax, asum  ! min, max values and sum of input array

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g1

      logical, dimension(:,:), allocatable :: mask

      integer (kind=int_kind) :: nx, ny

      integer (kind=int_kind) :: lnrec       ! local value of nrec

      logical (kind=log_kind) :: lrestart_ext  ! local value of restart_ext

      lnrec = nrec

      work = c0 ! to satisfy intent(out) attribute

      lrestart_ext = .false.
      if (present(restart_ext)) then
         lrestart_ext = restart_ext
      endif

      if (lrestart_ext) then
         nx = nx_global + 2*nghost
         ny = ny_global + 2*nghost
      else
         nx = nx_global
         ny = ny_global
      endif

      if (my_task == master_task) then
         allocate(work_g1(nx,ny))
      else
         allocate(work_g1(1,1))   ! to save memory
      endif

      if (my_task == master_task) then

         !-------------------------------------------------------------
         ! Find out ID of required variable
         !-------------------------------------------------------------

         status = nf90_inq_varid(fid, trim(varname), varid)
         call ice_check_nc(status, subname//' ERROR: Cannot find variable '//trim(varname), &
                           file=__FILE__, line=__LINE__)

         !-------------------------------------------------------------
         ! Check nrec axis size
         !-------------------------------------------------------------

         status = nf90_inquire_variable(fid, varid, ndims=ndims, dimids=dimids)
         call ice_check_nc(status, subname//' ERROR: inquire variable dimids '//trim(varname), &
                           file=__FILE__, line=__LINE__)
         if (ndims > 2) then
            status = nf90_inquire_dimension(fid, dimids(3), len=dimlen)
            call ice_check_nc(status, subname//' ERROR: inquire dimension size 3 '//trim(varname), &
                              file=__FILE__, line=__LINE__)
            if (lnrec > dimlen) then
               write(nu_diag,*) subname,' ERROR: not enough records, ',trim(varname),lnrec,dimlen
               call abort_ice(subname//' ERROR: not enough records '//trim(varname), &
                              file=__FILE__, line=__LINE__)
            endif
         endif

         !--------------------------------------------------------------
         ! Read global array
         !--------------------------------------------------------------

         status = nf90_get_var( fid, varid, work_g1, &
                                 start=(/1,1,lnrec/), count=(/nx,ny,1/))
         call ice_check_nc(status, subname//' ERROR: Cannot get variable '//trim(varname), &
                           file=__FILE__, line=__LINE__)

      endif                     ! my_task = master_task

      !-------------------------------------------------------------------
      ! optional diagnostics
      !-------------------------------------------------------------------

      if (my_task==master_task .and. diag) then
           write(nu_diag,'(2a,i8,a,i8,2a)') &
             subname,' fid= ',fid, ', lnrec = ',lnrec, &
             ', varname = ',trim(varname)
!          status = nf90_inquire(fid, nDimensions=ndim, nVariables=nvar)
!          write(nu_diag,*) subname,' ndim= ',ndim,', nvar= ',nvar
!          do id=1,ndim
!            status = nf90_inquire_dimension(fid,id,name=dimname,len=dimlen)
!            write(nu_diag,*) subname,' Dim name = ',trim(dimname),', size = ',dimlen
!         enddo
         ! optional
         missingvalue = spval_dbl
         status = nf90_get_att(fid, varid, "_FillValue", missingvalue)
!          call ice_check_nc(status, subname//' ERROR: Missing _FillValue', &
!                            file=__FILE__, line=__LINE__)
!         write(nu_diag,*) subname,' missingvalue= ',missingvalue
         allocate(mask(nx,ny))
         if ( ieee_is_nan(missingvalue) ) then
            mask = ieee_is_nan(work_g1)
         else
            mask = work_g1 /= missingvalue
         endif
         amin = minval(work_g1, mask = mask )
         amax = maxval(work_g1, mask = mask )
         asum = sum   (work_g1, mask = mask )
         write(nu_diag,*) subname,' min, max, sum =', amin, amax, asum, trim(varname)
         deallocate(mask)
      endif

      !-------------------------------------------------------------------
      ! Scatter data to individual processors.
      ! NOTE: Ghost cells are not updated unless field_loc is present.
      !-------------------------------------------------------------------

      if (lrestart_ext) then
         call scatter_global_ext(work, work_g1, master_task, distrb_info)
      else
         if (present(field_loc)) then
            call scatter_global(work, work_g1, master_task, distrb_info, &
                                field_loc, field_type)
         else
            call scatter_global(work, work_g1, master_task, distrb_info, &
                                field_loc_noupdate, field_type_noupdate)
         endif
      endif

      deallocate(work_g1)

! echmod:  this should not be necessary if fill/missing are only on land
      where (work > 1.0e+30_dbl_kind) work = c0

#else
      work = c0 ! to satisfy intent(out) attribute
      call abort_ice(subname//' ERROR: USE_NETCDF cpp not defined', &
         file=__FILE__, line=__LINE__)
#endif
      end subroutine ice_read_nc_xy

!=======================================================================

! Read a netCDF file and scatter to processors.
! If the optional variables field_loc and field_type are present,
! the ghost cells are filled using values from the global array.
! This prevents them from being filled with zeroes in land cells
! (subroutine ice_HaloUpdate need not be called).
!
! Adapted by David Bailey, NCAR from ice_read_nc_xy

      subroutine ice_read_nc_xyz(fid,  nrec,  varname, work,  diag, &
                                 field_loc, field_type, restart_ext)

      use ice_gather_scatter, only: scatter_global, scatter_global_ext

      integer (kind=int_kind), intent(in) :: &
         fid           , & ! file id
         nrec              ! record number

      character (len=*), intent(in) :: &
         varname           ! field name in netcdf file

      logical (kind=log_kind), intent(in) :: &
         diag              ! if true, write diagnostic output

      real (kind=dbl_kind), dimension(nx_block,ny_block,ncat,max_blocks), intent(out) :: &
         work              ! output array (real, 8-byte)

      logical (kind=log_kind), optional, intent(in) :: &
         restart_ext       ! if true, read extended grid

      integer (kind=int_kind), optional, intent(in) :: &
         field_loc, &      ! location of field on staggered grid
         field_type        ! type of field (scalar, vector, angle)

      ! local variables

      character(len=*), parameter :: subname = '(ice_read_nc_xyz)'

#ifdef USE_NETCDF
! netCDF file diagnostics:
      integer (kind=int_kind) :: &
         n,              & ! ncat index
         varid         , & ! variable id
         status        , & ! status output from netcdf routines
         ndims         , & ! number of dimensions
         dimlen            ! dimension size

      integer (kind=int_kind), dimension(10) :: &
         dimids            ! generic size dimids

      real (kind=dbl_kind) :: &
         missingvalue,   & ! missing value
         amin, amax, asum  ! min, max values and sum of input array

!     character (char_len) :: &
!        dimname           ! dimension name

      real (kind=dbl_kind), dimension(:,:,:), allocatable :: &
         work_g1

      logical, dimension(:,:), allocatable :: mask

      integer (kind=int_kind) :: nx, ny

      integer (kind=int_kind) :: lnrec       ! local value of nrec

      logical (kind=log_kind) :: lrestart_ext  ! local value of restart_ext

      lnrec = nrec

      lrestart_ext = .false.
      if (present(restart_ext)) then
         lrestart_ext = restart_ext
      endif

      if (lrestart_ext) then
         nx = nx_global + 2*nghost
         ny = ny_global + 2*nghost
      else
         nx = nx_global
         ny = ny_global
      endif

      if (my_task == master_task) then
         allocate(work_g1(nx,ny,ncat))
      else
         allocate(work_g1(1,1,ncat))   ! to save memory
      endif

      if (my_task == master_task) then

         !-------------------------------------------------------------
         ! Find out ID of required variable
         !-------------------------------------------------------------

         status = nf90_inq_varid(fid, trim(varname), varid)
         call ice_check_nc(status, subname//' ERROR: Cannot find variable '//trim(varname), &
                           file=__FILE__, line=__LINE__)

         !-------------------------------------------------------------
         ! Check nrec axis size
         !-------------------------------------------------------------

         status = nf90_inquire_variable(fid, varid, ndims=ndims, dimids=dimids)
         call ice_check_nc(status, subname//' ERROR: inquire variable dimids '//trim(varname), &
                           file=__FILE__, line=__LINE__)
         if (ndims > 3) then
            status = nf90_inquire_dimension(fid, dimids(4), len=dimlen)
            call ice_check_nc(status, subname//' ERROR: inquire dimension size 4 '//trim(varname), &
                              file=__FILE__, line=__LINE__)
            if (lnrec > dimlen) then
               write(nu_diag,*) subname,' ERROR: not enough records, ',trim(varname),lnrec,dimlen
               call abort_ice(subname//' ERROR: not enough records '//trim(varname), &
                              file=__FILE__, line=__LINE__)
            endif
         endif

         !--------------------------------------------------------------
         ! Read global array
         !--------------------------------------------------------------

         status = nf90_get_var( fid, varid, work_g1, &
                                 start=(/1,1,1,lnrec/), count=(/nx,ny,ncat,1/))
         call ice_check_nc(status, subname//' ERROR: Cannot get variable '//trim(varname), &
                           file=__FILE__, line=__LINE__)

      endif                     ! my_task = master_task

      !-------------------------------------------------------------------
      ! optional diagnostics
      !-------------------------------------------------------------------

      if (my_task==master_task .and. diag) then
           write(nu_diag,'(2a,i8,a,i8,2a)') &
              subname,' fid= ',fid, ', lnrec = ',lnrec, &
             ', varname = ',trim(varname)
!          status = nf90_inquire(fid, nDimensions=ndim, nVariables=nvar)
!          write(nu_diag,*) subname,' ndim= ',ndim,', nvar= ',nvar
!          do id=1,ndim
!            status = nf90_inquire_dimension(fid,id,name=dimname,len=dimlen)
!            write(nu_diag,*) subname,' Dim name = ',trim(dimname),', size = ',dimlen
!         enddo
         ! optional
         missingvalue = spval_dbl
         status = nf90_get_att(fid, varid, "_FillValue", missingvalue)
!          call ice_check_nc(status, subname//' ERROR: Missing _FillValue', &
!                            file=__FILE__, line=__LINE__)
         allocate(mask(nx,ny))
         do n=1,ncat
            if ( ieee_is_nan(missingvalue) ) then
               mask = ieee_is_nan(work_g1(:,:,n))
            else
               mask = work_g1(:,:,n) /= missingvalue
            endif
            amin = minval(work_g1(:,:,n), mask = mask )
            amax = maxval(work_g1(:,:,n), mask = mask )
            asum = sum   (work_g1(:,:,n), mask = mask )
            write(nu_diag,*) subname,' min, max, sum =', amin, amax, asum, trim(varname)
         enddo
         deallocate(mask)
      endif

      !-------------------------------------------------------------------
      ! Scatter data to individual processors.
      ! NOTE: Ghost cells are not updated unless field_loc is present.
      !-------------------------------------------------------------------

      if (lrestart_ext) then
         do n=1,ncat
            call scatter_global_ext(work(:,:,n,:), work_g1(:,:,n), &
                                    master_task, distrb_info)
         enddo
      else
         if (present(field_loc)) then
            do n=1,ncat
               call scatter_global(work(:,:,n,:), work_g1(:,:,n), master_task, &
                    distrb_info, field_loc, field_type)
            enddo
         else
            do n=1,ncat
               call scatter_global(work(:,:,n,:), work_g1(:,:,n), master_task, &
                    distrb_info, field_loc_noupdate, field_type_noupdate)
            enddo
         endif
      endif

      deallocate(work_g1)

#else
      call abort_ice(subname//' ERROR: USE_NETCDF cpp not defined', &
         file=__FILE__, line=__LINE__)
      work = c0 ! to satisfy intent(out) attribute
#endif
      end subroutine ice_read_nc_xyz

!=======================================================================

! Read a netCDF file and scatter to processors.
! If the optional variables field_loc and field_type are present,
! the ghost cells are filled using values from the global array.
! This prevents them from being filled with zeroes in land cells
! (subroutine ice_HaloUpdate need not be called).
!
! Adapted by David Bailey, NCAR from ice_read_nc_xy
! Adapted by Lettie Roach, NIWA to read nfreq
! by changing all occurrences of ncat to nfreq

      subroutine ice_read_nc_xyf(fid,  nrec,  varname, work,  diag, &
                                 field_loc, field_type, restart_ext)

      use ice_fileunits, only: nu_diag
      use ice_domain_size, only: nfreq
      use ice_gather_scatter, only: scatter_global, scatter_global_ext

      integer (kind=int_kind), intent(in) :: &
         fid           , & ! file id
         nrec              ! record number

      character (len=*), intent(in) :: &
         varname           ! field name in netcdf file

      logical (kind=log_kind), intent(in) :: &
         diag              ! if true, write diagnostic output

      real (kind=dbl_kind), dimension(nx_block,ny_block,nfreq,1,max_blocks), intent(out) :: &
         work              ! output array (real, 8-byte)

      logical (kind=log_kind), optional, intent(in) :: &
         restart_ext       ! if true, read extended grid

      integer (kind=int_kind), optional, intent(in) :: &
         field_loc, &      ! location of field on staggered grid
         field_type        ! type of field (scalar, vector, angle)

      ! local variables

! netCDF file diagnostics:
      integer (kind=int_kind) :: &
         varid,          & ! variable id
         status,         & ! status output from netcdf routines
         ndim, nvar,     & ! sizes of netcdf file
         id,             & ! dimension index
         n,              & ! ncat index
         ndims,          & ! number of dimensions
         dimlen            ! dimension size

      integer (kind=int_kind), dimension(10) :: &
         dimids            ! generic size dimids

      real (kind=dbl_kind) :: &
         missingvalue,   & ! missing value
         amin, amax, asum  ! min, max values and sum of input array

      character (char_len) :: &
         dimname           ! dimension name

      real (kind=dbl_kind), dimension(:,:,:), allocatable :: &
         work_g1

      logical, dimension(:,:), allocatable :: mask

      integer (kind=int_kind) :: nx, ny

      integer (kind=int_kind) :: lnrec       ! local value of nrec

      logical (kind=log_kind) :: lrestart_ext  ! local value of restart_ext

      character(len=*), parameter :: subname = '(ice_read_nc_xyf)'

#ifdef USE_NETCDF

      lnrec = nrec

      lrestart_ext = .false.
      if (present(restart_ext)) then
         lrestart_ext = restart_ext
      endif

      if (lrestart_ext) then
         nx = nx_global + 2*nghost
         ny = ny_global + 2*nghost
      else
         nx = nx_global
         ny = ny_global
      endif

      if (my_task == master_task) then
         allocate(work_g1(nx,ny,nfreq))
      else
         allocate(work_g1(1,1,nfreq))   ! to save memory
      endif

      if (my_task == master_task) then

         !-------------------------------------------------------------
         ! Find out ID of required variable
         !-------------------------------------------------------------

         status = nf90_inq_varid(fid, trim(varname), varid)
         call ice_check_nc(status, subname//' ERROR: Cannot find variable '//trim(varname), &
                           file=__FILE__, line=__LINE__)

         !-------------------------------------------------------------
         ! Check nrec axis size
         !-------------------------------------------------------------

         status = nf90_inquire_variable(fid, varid, ndims=ndims, dimids=dimids)
         call ice_check_nc(status, subname//' ERROR: inquire variable dimids '//trim(varname), &
                           file=__FILE__, line=__LINE__)
         if (ndims > 3) then
            status = nf90_inquire_dimension(fid, dimids(4), len=dimlen)
            call ice_check_nc(status, subname//' ERROR: inquire dimension size 4 '//trim(varname), &
                              file=__FILE__, line=__LINE__)
            if (lnrec > dimlen) then
               write(nu_diag,*) subname,' ERROR: not enough records, ',trim(varname),lnrec,dimlen
               call abort_ice(subname//' ERROR: not enough records '//trim(varname), &
                              file=__FILE__, line=__LINE__)
            endif
         endif

         !--------------------------------------------------------------
         ! Read global array
         !--------------------------------------------------------------

         status = nf90_get_var( fid, varid, work_g1, &
                                 start=(/1,1,1,lnrec/), count=(/nx,ny,nfreq,1/))
         call ice_check_nc(status, subname//' ERROR: Cannot get variable '//trim(varname), &
                           file=__FILE__, line=__LINE__)

      endif                     ! my_task = master_task

      !-------------------------------------------------------------------
      ! optional diagnostics
      !-------------------------------------------------------------------

      if (my_task==master_task .and. diag) then
         write(nu_diag,'(2a,i8,a,i8,2a)') &
            subname,' fid= ',fid, ', lnrec = ',lnrec, &
           ', varname = ',trim(varname)
         status = nf90_inquire(fid, nDimensions=ndim, nVariables=nvar)
         write(nu_diag,*) subname,' ndim= ',ndim,', nvar= ',nvar
         do id=1,ndim
            status = nf90_inquire_dimension(fid,id,name=dimname,len=dimlen)
            write(nu_diag,*) subname,' Dim name = ',trim(dimname),', size = ',dimlen
         enddo
         ! optional
         missingvalue = spval_dbl
         status = nf90_get_att(fid, varid, "_FillValue", missingvalue)
!          call ice_check_nc(status, subname//' ERROR: Missing _FillValue', &
!                            file=__FILE__, line=__LINE__)
         allocate(mask(nx,ny))
         do n=1,ncat
            if ( ieee_is_nan(missingvalue) ) then
               mask = ieee_is_nan(work_g1(:,:,n))
            else
               mask = work_g1(:,:,n) /= missingvalue
            endif
            amin = minval(work_g1(:,:,n), mask = mask )
            amax = maxval(work_g1(:,:,n), mask = mask )
            asum = sum   (work_g1(:,:,n), mask = mask )
            write(nu_diag,*) subname,' min, max, sum =', amin, amax, asum, trim(varname)
         enddo
         deallocate(mask)
      endif

      !-------------------------------------------------------------------
      ! Scatter data to individual processors.
      ! NOTE: Ghost cells are not updated unless field_loc is present.
      !-------------------------------------------------------------------

      if (lrestart_ext) then
         do n = 1, nfreq
            call scatter_global_ext(work(:,:,n,1,:), work_g1(:,:,n), &
                                    master_task, distrb_info)
         enddo
      else
         if (present(field_loc)) then
            do n = 1, nfreq
               call scatter_global(work(:,:,n,1,:), work_g1(:,:,n), master_task, &
                    distrb_info, field_loc, field_type)
            enddo
         else
            do n = 1, nfreq
               call scatter_global(work(:,:,n,1,:), work_g1(:,:,n), master_task, &
                    distrb_info, field_loc_noupdate, field_type_noupdate)
            enddo
         endif
      endif

! echmod:  this should not be necessary if fill/missing are only on land
      where (work > 1.0e+30_dbl_kind) work = c0

      deallocate(work_g1)

#else
      call abort_ice(subname//' ERROR: USE_NETCDF cpp not defined', &
         file=__FILE__, line=__LINE__)
      work = c0 ! to satisfy intent(out) attribute
#endif

     end subroutine ice_read_nc_xyf

!=======================================================================

! Read a netCDF file
! Adapted by Alison McLaren, Met Office from ice_read

      subroutine ice_read_nc_point(fid,  nrec,  varname, work, diag, &
                                   field_loc, field_type)

      integer (kind=int_kind), intent(in) :: &
         fid           , & ! file id
         nrec              ! record number

      logical (kind=log_kind), intent(in) :: &
         diag              ! if true, write diagnostic output

      character (char_len), intent(in) :: &
         varname           ! field name in netcdf file

      integer (kind=int_kind), optional, intent(in) :: &
         field_loc, &      ! location of field on staggered grid
         field_type        ! type of field (scalar, vector, angle)

      real (kind=dbl_kind), intent(out) :: &
         work              ! output variable (real, 8-byte)

      ! local variables

      character(len=*), parameter :: subname = '(ice_read_nc_point)'

#ifdef USE_NETCDF
! netCDF file diagnostics:
      integer (kind=int_kind) :: &
         varid,          & ! netcdf id for field
         status,         & ! status output from netcdf routines
         ndim, nvar,     & ! sizes of netcdf file
         id,             & ! dimension index
         ndims,          & ! number of dimensions
         dimlen            ! dimension size

      integer (kind=int_kind), dimension(10) :: &
         dimids            ! generic size dimids

      real (kind=dbl_kind), dimension(1) :: &
         workg             ! temporary work variable

      integer (kind=int_kind) :: lnrec       ! local value of nrec

      character (char_len) :: &
         dimname           ! dimension name

      lnrec = nrec

      if (my_task == master_task) then

         !-------------------------------------------------------------
         ! Find out ID of required variable
         !-------------------------------------------------------------

         status = nf90_inq_varid(fid, trim(varname), varid)
         call ice_check_nc(status, subname//' ERROR: Cannot find variable '//trim(varname), &
                           file=__FILE__, line=__LINE__)

         !-------------------------------------------------------------
         ! Check nrec axis size
         !-------------------------------------------------------------

         status = nf90_inquire_variable(fid, varid, ndims=ndims, dimids=dimids)
         call ice_check_nc(status, subname//' ERROR: inquire variable dimids '//trim(varname), &
                           file=__FILE__, line=__LINE__)
         if (ndims > 0) then
            status = nf90_inquire_dimension(fid, dimids(1), len=dimlen)
            call ice_check_nc(status, subname//' ERROR: inquire dimension size 1 '//trim(varname), &
                              file=__FILE__, line=__LINE__)
            if (lnrec > dimlen) then
               write(nu_diag,*) subname,' ERROR: not enough records, ',trim(varname),lnrec,dimlen
               call abort_ice(subname//' ERROR: not enough records '//trim(varname), &
                              file=__FILE__, line=__LINE__)
            endif
         endif

         !--------------------------------------------------------------
         ! Read point variable
         !--------------------------------------------------------------

         status = nf90_get_var(fid, varid, workg, &
                               start= (/ lnrec /), count=(/ 1 /))
         call ice_check_nc(status, subname//' ERROR: Cannot get variable '//trim(varname), &
                           file=__FILE__, line=__LINE__)
      endif                     ! my_task = master_task

      !-------------------------------------------------------------------
      ! optional diagnostics
      !-------------------------------------------------------------------

      if (my_task==master_task .and. diag) then
          write(nu_diag,'(2a,i8,a,i8,2a)') &
             subname,' fid= ',fid, ', lnrec = ',lnrec, &
            ', varname = ',trim(varname)
          status = nf90_inquire(fid, nDimensions=ndim, nVariables=nvar)
          write(nu_diag,*) subname,' ndim= ',ndim,', nvar= ',nvar
          do id=1,ndim
            status = nf90_inquire_dimension(fid,id,name=dimname,len=dimlen)
            write(nu_diag,*) subname,' Dim name = ',trim(dimname),', size = ',dimlen
         enddo
      endif

      work = workg(1)

#else
      call abort_ice(subname//' ERROR: USE_NETCDF cpp not defined', &
         file=__FILE__, line=__LINE__)
      work = c0 ! to satisfy intent(out) attribute
#endif
      end subroutine ice_read_nc_point

!=======================================================================

! Written by T. Craig

      subroutine ice_read_nc_1D(fid,  varname, work, diag, &
                                xdim)

      use ice_fileunits, only: nu_diag

      integer (kind=int_kind), intent(in) :: &
         fid           , & ! file id
         xdim              ! field dimensions

      logical (kind=log_kind), intent(in) :: &
         diag              ! if true, write diagnostic output

      character (char_len), intent(in) :: &
         varname           ! field name in netcdf file

      real (kind=dbl_kind), dimension(:), intent(out) :: &
         work              ! output array

      ! local variables

      character(len=*), parameter :: subname = '(ice_read_nc_1D)'

#ifdef USE_NETCDF
! netCDF file diagnostics:
      integer (kind=int_kind) :: &
         varid,          & ! netcdf id for field
         status,         & ! status output from netcdf routines
         ndim, nvar        ! sizes of netcdf file

      real (kind=dbl_kind), dimension(xdim) :: &
         workg             ! output array (real, 8-byte)

      !--------------------------------------------------------------

      if (my_task == master_task) then

         if (size(work,dim=1) < xdim) then
            write(nu_diag,*) subname,' work, dim=1 ',size(work,dim=1),xdim
            call abort_ice (subname//' ERROR: work array wrong size '//trim(varname), &
                            file=__FILE__, line=__LINE__ )
         endif

         !-------------------------------------------------------------
         ! Find out ID of required variable
         !-------------------------------------------------------------

         status = nf90_inq_varid(fid, trim(varname), varid)
         call ice_check_nc(status, subname//' ERROR: Cannot find variable '//trim(varname), &
                           file=__FILE__, line=__LINE__ )

         !--------------------------------------------------------------
         ! Read  array
         !--------------------------------------------------------------

         status = nf90_get_var( fid, varid, workg, &
                                start=(/1/), count=(/xdim/) )
         call ice_check_nc(status, subname//' ERROR: Cannot get variable '//trim(varname), &
                           file=__FILE__, line=__LINE__ )
         work(1:xdim) = workg(1:xdim)

         !-------------------------------------------------------------------
         ! optional diagnostics
         !-------------------------------------------------------------------

         if (diag) then
            write(nu_diag,*) subname, &
              ' fid= ',fid, ', xdim = ',xdim, &
              ' varname = ',trim(varname)
            status = nf90_inquire(fid, nDimensions=ndim, nVariables=nvar)
            write(nu_diag,*) 'ndim= ',ndim,', nvar= ',nvar
         endif
      endif
#else
      call abort_ice(subname//' ERROR: USE_NETCDF cpp not defined', &
                     file=__FILE__, line=__LINE__)
      work = c0 ! to satisfy intent(out) attribute
#endif

      end subroutine ice_read_nc_1D

!=======================================================================

! Written by T. Craig

      subroutine ice_read_nc_2D(fid,  varname, work, diag, &
                                xdim, ydim)

      use ice_fileunits, only: nu_diag

      integer (kind=int_kind), intent(in) :: &
         fid           , & ! file id
         xdim, ydim        ! field dimensions

      logical (kind=log_kind), intent(in) :: &
         diag              ! if true, write diagnostic output

      character (char_len), intent(in) :: &
         varname           ! field name in netcdf file

      real (kind=dbl_kind), dimension(:,:), intent(out) :: &
         work              ! output array

      ! local variables

      character(len=*), parameter :: subname = '(ice_read_nc_2D)'

#ifdef USE_NETCDF
! netCDF file diagnostics:
      integer (kind=int_kind) :: &
         varid,          & ! netcdf id for field
         status,         & ! status output from netcdf routines
         ndim, nvar        ! sizes of netcdf file

      real (kind=dbl_kind), dimension(xdim,ydim) :: &
         workg             ! output array (real, 8-byte)

      !--------------------------------------------------------------

      if (my_task == master_task) then

         if (size(work,dim=1) < xdim .or. &
             size(work,dim=2) < ydim) then
            write(nu_diag,*) subname,' work, dim=1 ',size(work,dim=1),xdim
            write(nu_diag,*) subname,' work, dim=2 ',size(work,dim=2),ydim
            call abort_ice (subname//' ERROR: work array wrong size '//trim(varname), &
                            file=__FILE__, line=__LINE__ )
         endif

         !-------------------------------------------------------------
         ! Find out ID of required variable
         !-------------------------------------------------------------

         status = nf90_inq_varid(fid, trim(varname), varid)
         call ice_check_nc(status,subname//' ERROR: Cannot find variable '//trim(varname), &
                           file=__FILE__, line=__LINE__ )

         !--------------------------------------------------------------
         ! Read  array
         !--------------------------------------------------------------

         status = nf90_get_var( fid, varid, workg, &
                                start=(/1,1/), count=(/xdim,ydim/) )
         call ice_check_nc(status, subname//' ERROR: Cannot get variable '//trim(varname), &
                           file=__FILE__, line=__LINE__ )
         work(1:xdim,1:ydim) = workg(1:xdim, 1:ydim)

         !-------------------------------------------------------------------
         ! optional diagnostics
         !-------------------------------------------------------------------

         if (diag) then
            write(nu_diag,*) subname, &
              ' fid= ',fid, ', xdim = ',xdim, &
              ' ydim= ', ydim, ' varname = ',trim(varname)
            status = nf90_inquire(fid, nDimensions=ndim, nVariables=nvar)
            write(nu_diag,*) 'ndim= ',ndim,', nvar= ',nvar
         endif
      endif
#else
      call abort_ice(subname//' ERROR: USE_NETCDF cpp not defined', &
                     file=__FILE__, line=__LINE__)
      work = c0 ! to satisfy intent(out) attribute
#endif

      end subroutine ice_read_nc_2D

!=======================================================================

! Written by T. Craig

      subroutine ice_read_nc_3D(fid,  varname, work, diag, &
                                xdim, ydim, zdim)

      use ice_fileunits, only: nu_diag

      integer (kind=int_kind), intent(in) :: &
         fid           , & ! file id
         xdim, ydim,zdim   ! field dimensions

      logical (kind=log_kind), intent(in) :: &
         diag              ! if true, write diagnostic output

      character (char_len), intent(in) :: &
         varname           ! field name in netcdf file

      real (kind=dbl_kind), dimension(:,:,:), intent(out) :: &
         work              ! output array

      ! local variables

      character(len=*), parameter :: subname = '(ice_read_nc_3D)'

#ifdef USE_NETCDF
! netCDF file diagnostics:
      integer (kind=int_kind) :: &
         varid,          & ! netcdf id for field
         status,         & ! status output from netcdf routines
         ndim, nvar        ! sizes of netcdf file

      real (kind=dbl_kind), dimension(xdim,ydim,zdim) :: &
         workg             ! output array (real, 8-byte)

      !--------------------------------------------------------------

      if (my_task == master_task) then

         if (size(work,dim=1) < xdim .or. &
             size(work,dim=2) < ydim .or. &
             size(work,dim=3) < zdim ) then
            write(nu_diag,*) subname,' work, dim=1 ',size(work,dim=1),xdim
            write(nu_diag,*) subname,' work, dim=2 ',size(work,dim=2),ydim
            write(nu_diag,*) subname,' work, dim=3 ',size(work,dim=3),zdim
            call abort_ice (subname//' ERROR: work array wrong size '//trim(varname), &
                            file=__FILE__, line=__LINE__ )
         endif

         !-------------------------------------------------------------
         ! Find out ID of required variable
         !-------------------------------------------------------------

         status = nf90_inq_varid(fid, trim(varname), varid)
         call ice_check_nc(status, subname//' ERROR: Cannot find variable '//trim(varname), &
                           file=__FILE__, line=__LINE__ )

         !--------------------------------------------------------------
         ! Read  array
         !--------------------------------------------------------------

         status = nf90_get_var( fid, varid, workg, &
                                start=(/1,1,1/), count=(/xdim,ydim,zdim/) )
         call ice_check_nc(status, subname//' ERROR: Cannot get variable '//trim(varname), &
                           file=__FILE__, line=__LINE__ )
         work(1:xdim,1:ydim,1:zdim) = workg(1:xdim, 1:ydim, 1:zdim)

         !-------------------------------------------------------------------
         ! optional diagnostics
         !-------------------------------------------------------------------

         if (diag) then
            write(nu_diag,*) subname, &
              ' fid= ',fid, ', xdim = ',xdim, &
              ' ydim= ', ydim,' zdim = ',zdim, ' varname = ',trim(varname)
            status = nf90_inquire(fid, nDimensions=ndim, nVariables=nvar)
            write(nu_diag,*) 'ndim= ',ndim,', nvar= ',nvar
         endif
      endif
#else
      call abort_ice(subname//' ERROR: USE_NETCDF cpp not defined', &
                     file=__FILE__, line=__LINE__)
      work = c0 ! to satisfy intent(out) attribute
#endif

      end subroutine ice_read_nc_3D

!=======================================================================

! Adapted by Nicole Jeffery, LANL

      subroutine ice_read_nc_z(fid,  nrec,  varname, work,  diag, &
                               field_loc, field_type)

      use ice_domain_size, only: nilyr

      integer (kind=int_kind), intent(in) :: &
         fid           , & ! file id
         nrec              ! record number

      logical (kind=log_kind), intent(in) :: &
         diag              ! if true, write diagnostic output

      character (char_len), intent(in) :: &
         varname           ! field name in netcdf file

      integer (kind=int_kind), optional, intent(in) :: &
         field_loc, &      ! location of field on staggered grid
         field_type        ! type of field (scalar, vector, angle)

      real (kind=dbl_kind), dimension(nilyr), intent(out) :: &
         work              ! output array (real, 8-byte)

      ! local variables

#ifdef USE_NETCDF
      real (kind=dbl_kind), dimension(:), allocatable :: &
         work_z

! netCDF file diagnostics:
      integer (kind=int_kind) :: &
         varid,          & ! netcdf id for field
         status,         & ! status output from netcdf routines
         ndim, nvar,     & ! sizes of netcdf file
         id,             & ! dimension index
         ndims,          & ! number of dimensions
         dimlen            ! dimension size

      integer (kind=int_kind), dimension(10) :: &
         dimids            ! generic size dimids

      character (char_len) :: &
         dimname           ! dimension name

      integer (kind=int_kind) :: lnrec       ! local value of nrec

#endif

      character(len=*), parameter :: subname = '(ice_read_nc_z)'

#ifdef USE_NETCDF

      lnrec = nrec

      allocate(work_z(nilyr))

      if (my_task == master_task) then

         !-------------------------------------------------------------
         ! Find out ID of required variable
         !-------------------------------------------------------------

         status = nf90_inq_varid(fid, trim(varname), varid)
         call ice_check_nc(status, subname//' ERROR: Cannot find variable '//trim(varname), &
                           file=__FILE__, line=__LINE__)

         !-------------------------------------------------------------
         ! Check nrec axis size
         !-------------------------------------------------------------

         status = nf90_inquire_variable(fid, varid, ndims=ndims, dimids=dimids)
         call ice_check_nc(status, subname//' ERROR: inquire variable dimids '//trim(varname), &
                           file=__FILE__, line=__LINE__)
         if (ndims > 1) then
            status = nf90_inquire_dimension(fid, dimids(2), len=dimlen)
            call ice_check_nc(status, subname//' ERROR: inquire dimension size 2 '//trim(varname), &
                              file=__FILE__, line=__LINE__)
            if (lnrec > dimlen) then
               write(nu_diag,*) subname,' ERROR: not enough records, ',trim(varname),lnrec,dimlen
               call abort_ice(subname//' ERROR: not enough records '//trim(varname), &
                              file=__FILE__, line=__LINE__)
            endif
         endif

         !--------------------------------------------------------------
         ! Read global array
         !--------------------------------------------------------------

         status = nf90_get_var( fid, varid, work_z, &
                                start=(/1,lnrec/), count=(/nilyr,1/))
         call ice_check_nc(status, subname//' ERROR: Cannot get variable '//trim(varname), &
                           file=__FILE__, line=__LINE__)
      endif                     ! my_task = master_task

      !-------------------------------------------------------------------
      ! optional diagnostics
      !-------------------------------------------------------------------

      if (my_task==master_task .and. diag) then
          write(nu_diag,'(2a,i8,a,i8,2a)') &
             subname,' fid= ',fid, ', lnrec = ',lnrec, &
            ', varname = ',trim(varname)
          status = nf90_inquire(fid, nDimensions=ndim, nVariables=nvar)
          write(nu_diag,*) subname,' ndim= ',ndim,', nvar= ',nvar
          do id=1,ndim
            status = nf90_inquire_dimension(fid,id,name=dimname,len=dimlen)
            write(nu_diag,*) subname,' Dim name = ',trim(dimname),', size = ',dimlen
         enddo
      endif

      work(:) = work_z(:)
      deallocate(work_z)

#else
      call abort_ice(subname//' ERROR: USE_NETCDF cpp not defined', &
         file=__FILE__, line=__LINE__)
      work = c0 ! to satisfy intent(out) attribute
#endif
      end subroutine ice_read_nc_z

!=======================================================================

! Write a netCDF file.
!
! Adapted by David Bailey, NCAR

      subroutine ice_write_nc_xy(fid,  nrec,  varid, work,  diag, &
                                 restart_ext, varname)

      use ice_gather_scatter, only: gather_global, gather_global_ext

      integer (kind=int_kind), intent(in) :: &
         fid           , & ! file id
         varid         , & ! variable id
         nrec              ! record number

      logical (kind=log_kind), intent(in) :: &
         diag              ! if true, write diagnostic output

      logical (kind=log_kind), optional, intent(in) :: &
         restart_ext       ! if true, write extended grid

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks), intent(in) :: &
         work              ! output array (real, 8-byte)

      character (len=*), optional, intent(in) :: &
         varname           ! variable name

      ! local variables

      character(len=*), parameter :: subname = '(ice_write_nc_xy)'

#ifdef USE_NETCDF
! netCDF file diagnostics:
      integer (kind=int_kind) :: &
         status            ! status output from netcdf routines
!        ndim, nvar,     & ! sizes of netcdf file
!        id,             & ! dimension index
!        dimlen            ! size of dimension

      real (kind=dbl_kind) :: &
         amin, amax, asum  ! min, max values and sum of input array

      character (char_len) :: &
         lvarname          ! variable name
!        dimname           ! dimension name

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g1

      integer (kind=int_kind) :: nx, ny

      logical (kind=log_kind) :: lrestart_ext  ! local value of restart_ext

      lrestart_ext = .false.
      if (present(restart_ext)) then
         lrestart_ext = restart_ext
      endif

      if (lrestart_ext) then
         nx = nx_global + 2*nghost
         ny = ny_global + 2*nghost
      else
         nx = nx_global
         ny = ny_global
      endif

      if (present(varname)) then
         lvarname = trim(varname)
      else
         lvarname = ' '
      endif

      if (my_task == master_task) then
         allocate(work_g1(nx,ny))
      else
         allocate(work_g1(1,1))   ! to save memory
      endif

      if (lrestart_ext) then
         call gather_global_ext(work_g1, work, master_task, distrb_info, spc_val=c0)
      else
         call gather_global(work_g1, work, master_task, distrb_info, spc_val=c0)
      endif

      if (my_task == master_task) then

         !--------------------------------------------------------------
         ! Write global array
         !--------------------------------------------------------------

         status = nf90_put_var( fid, varid, work_g1, &
                                start=(/1,1,nrec/), count=(/nx,ny,1/))
         call ice_check_nc(status, subname//' ERROR: Cannot put variable ', &
                           file=__FILE__, line=__LINE__ )
      endif                     ! my_task = master_task

      !-------------------------------------------------------------------
      ! optional diagnostics
      !-------------------------------------------------------------------

      if (my_task==master_task .and. diag) then
!          write(nu_diag,*) &
!             subname,' fid= ',fid, ', nrec = ',nrec, &
!            ', varid = ',varid
!          status = nf90_inquire(fid, nDimensions=ndim, nVariables=nvar)
!          write(nu_diag,*) subname,' ndim= ',ndim,', nvar= ',nvar
!          do id=1,ndim
!            status = nf90_inquire_dimension(fid,id,name=dimname,len=dimlen)
!            write(nu_diag,*) subname,' Dim name = ',trim(dimname),', size = ',dimlen
!         enddo
         amin = minval(work_g1)
         amax = maxval(work_g1, mask = work_g1 /= spval_dbl)
         asum = sum   (work_g1, mask = work_g1 /= spval_dbl)
         write(nu_diag,*) subname,' min, max, sum =', amin, amax, asum, trim(lvarname)
      endif

      deallocate(work_g1)

#else
      call abort_ice(subname//' ERROR: USE_NETCDF cpp not defined', &
         file=__FILE__, line=__LINE__)
#endif

      end subroutine ice_write_nc_xy

!=======================================================================

! Write a netCDF file.
!
! Adapted by David Bailey, NCAR

      subroutine ice_write_nc_xyz(fid,  nrec,  varid, work,  diag, &
                                  restart_ext, varname)

      use ice_gather_scatter, only: gather_global, gather_global_ext

      integer (kind=int_kind), intent(in) :: &
         fid           , & ! file id
         varid         , & ! variable id
         nrec              ! record number

      logical (kind=log_kind), intent(in) :: &
         diag              ! if true, write diagnostic output

      logical (kind=log_kind), optional, intent(in) :: &
         restart_ext       ! if true, read extended grid

      real (kind=dbl_kind), dimension(nx_block,ny_block,ncat,max_blocks), intent(in) :: &
         work              ! output array (real, 8-byte)

      character (len=*), optional, intent(in) :: &
         varname           ! variable name

      ! local variables

      character(len=*), parameter :: subname = '(ice_write_nc_xyz)'

#ifdef USE_NETCDF
! netCDF file diagnostics:
      integer (kind=int_kind) :: &
         n,              & ! ncat index
         status            ! status output from netcdf routines
!        ndim, nvar,     & ! sizes of netcdf file
!        id,             & ! dimension index
!        dimlen            ! size of dimension

      real (kind=dbl_kind) :: &
         amin, amax, asum  ! min, max values and sum of input array

      character (char_len) :: &
         lvarname          ! variable name
!        dimname           ! dimension name

      real (kind=dbl_kind), dimension(:,:,:), allocatable :: &
         work_g1

      integer (kind=int_kind) :: nx, ny

      logical (kind=log_kind) :: lrestart_ext  ! local value of restart_ext

      lrestart_ext = .false.
      if (present(restart_ext)) then
         lrestart_ext = restart_ext
      endif

      if (lrestart_ext) then
         nx = nx_global + 2*nghost
         ny = ny_global + 2*nghost
      else
         nx = nx_global
         ny = ny_global
      endif

      if (my_task == master_task) then
         allocate(work_g1(nx,ny,ncat))
      else
         allocate(work_g1(1,1,ncat))   ! to save memory
      endif

      if (lrestart_ext) then
         do n=1,ncat
            call gather_global_ext(work_g1(:,:,n), work(:,:,n,:), &
                 master_task, distrb_info, spc_val=c0)
         enddo
      else
         do n=1,ncat
            call gather_global(work_g1(:,:,n), work(:,:,n,:), &
                    master_task, distrb_info, spc_val=c0)
         enddo
      endif

      if (present(varname)) then
         lvarname = trim(varname)
      else
         lvarname = ' '
      endif

      if (my_task == master_task) then

         !--------------------------------------------------------------
         ! Write global array
         !--------------------------------------------------------------

         status = nf90_put_var( fid, varid, work_g1, &
                                start=(/1,1,1,nrec/), count=(/nx,ny,ncat,1/))
         call ice_check_nc(status, subname//' ERROR: Cannot put variable ', &
                           file=__FILE__, line=__LINE__ )
      endif                     ! my_task = master_task

      !-------------------------------------------------------------------
      ! optional diagnostics
      !-------------------------------------------------------------------

      if (my_task==master_task .and. diag) then
!          write(nu_diag,*) &
!             subname,' fid= ',fid, ', nrec = ',nrec, &
!            ', varid = ',varid
!          status = nf90_inquire(fid, nDimensions=ndim, nVariables=nvar)
!          write(nu_diag,*) subname,' ndim= ',ndim,', nvar= ',nvar
!          do id=1,ndim
!            status = nf90_inquire_dimension(fid,id,name=dimname,len=dimlen)
!            write(nu_diag,*) subname,' Dim name = ',trim(dimname),', size = ',dimlen
!         enddo
         amin =  10000._dbl_kind
         amax = -10000._dbl_kind
         do n=1,ncat
            amin = minval(work_g1(:,:,n))
            amax = maxval(work_g1(:,:,n), mask = work_g1(:,:,n) /= spval_dbl)
            asum = sum   (work_g1(:,:,n), mask = work_g1(:,:,n) /= spval_dbl)
            write(nu_diag,*) subname,' min, max, sum =', amin, amax, asum, trim(lvarname)
         enddo
      endif

      deallocate(work_g1)

#else
      call abort_ice(subname//' ERROR: USE_NETCDF cpp not defined', &
         file=__FILE__, line=__LINE__)
#endif

      end subroutine ice_write_nc_xyz

!=======================================================================

! Read a netcdf file.
! Just like ice_read_nc except that it returns a global array.
! work_g is a real array
!
! Adapted by William Lipscomb, LANL, from ice_read
! Adapted by Ann Keen, Met Office, to read from a netcdf file

      subroutine ice_read_global_nc (fid,  nrec, varname, work_g, diag)

      integer (kind=int_kind), intent(in) :: &
         fid           , & ! file id
         nrec              ! record number

      character (char_len), intent(in) :: &
         varname           ! field name in netcdf file

      real (kind=dbl_kind), dimension(:,:), intent(out) :: &
         work_g            ! output array (real, 8-byte)

      logical (kind=log_kind) :: &
         diag              ! if true, write diagnostic output

      ! local variables

      character(len=*), parameter :: subname = '(ice_read_global_nc)'

#ifdef USE_NETCDF
      ! netCDF file diagnostics:
      integer (kind=int_kind) :: &
         varid,  &         ! netcdf id for field
         status, &         ! status output from netcdf routines
         ndim,  &         ! number of variable dimensions
         dimids(NF90_MAX_VAR_DIMS) , & !ids of dimensions
         dimlen            ! size of dimension
!        nvar,     & ! sizes of netcdf file
!        id,             & ! dimension index

      real (kind=dbl_kind) :: &
         amin, amax, asum  ! min, max values and sum of input array

!    character (char_len) :: &
!        dimname           ! dimension name
!

      work_g(:,:) = c0

      if (my_task == master_task) then

         !-------------------------------------------------------------
         ! Find out ID of required variable
         !-------------------------------------------------------------

         status = nf90_inq_varid(fid, trim(varname), varid)
         call ice_check_nc(status, subname//' ERROR: Cannot find variable '//trim(varname), &
                           file=__FILE__, line=__LINE__)

         !--------------------------------------------------------------
         ! Read global array
         !--------------------------------------------------------------

         ! Check var size : is var 2d ?
         status = nf90_inquire_variable(fid, varid, ndims=ndim, dimids=dimids)
         call ice_check_nc(status, subname//' ERROR: Cannot check variable '//trim(varname), &
                           file=__FILE__, line=__LINE__)
         if ( ndim > 2 ) then
            call abort_ice(subname//' ERROR: '//trim(varname)//' cannot have more than 2 dimensions', &
                           file=__FILE__, line=__LINE__)
         endif
         ! Is work_g the same size as the variable?
         status = nf90_inquire_dimension(fid, dimids(1), len=dimlen)
         call ice_check_nc(status, subname//' ERROR: Cannot check variable '//trim(varname), &
                           file=__FILE__, line=__LINE__)
         if ( dimlen /= size(work_g,1) ) then
            call abort_ice(subname//' ERROR: x dim of '//trim(varname)//' wrong size, check nx_global', &
                           file=__FILE__, line=__LINE__)
         endif
         status = nf90_inquire_dimension(fid, dimids(2), len=dimlen)
         call ice_check_nc(status, subname//' ERROR: Cannot check variable '//trim(varname), &
                           file=__FILE__, line=__LINE__)
         if ( dimlen /= size(work_g,2) ) then
            call abort_ice(subname//' ERROR: y dim of '//trim(varname)//' wrong size, check ny_global', &
                           file=__FILE__, line=__LINE__)
         endif

         ! Get the data
         status = nf90_get_var( fid, varid, work_g, start=(/1,1,nrec/))
         call ice_check_nc(status, subname//' ERROR: Cannot get variable '//trim(varname), &
                           file=__FILE__, line=__LINE__)

      endif                     ! my_task = master_task

      !-------------------------------------------------------------------
      ! optional diagnostics
      !-------------------------------------------------------------------

      if (my_task == master_task .and. diag) then
!          write(nu_diag,*) &
!             subname,' fid= ',fid, ', nrec = ',nrec, &
!            ', varname = ',trim(varname)
!          status = nf90_inquire(fid, nDimensions=ndim, nVariables=nvar)
!          write(nu_diag,*) subname,' ndim= ',ndim,', nvar= ',nvar
!          do id=1,ndim
!            status = nf90_inquire_dimension(fid,id,name=dimname,len=dimlen)
!            write(nu_diag,*) subname,' Dim name = ',trim(dimname),', size = ',dimlen
!         enddo
         amin = minval(work_g)
         amax = maxval(work_g, mask = work_g /= spval_dbl)
         asum = sum   (work_g, mask = work_g /= spval_dbl)
         write(nu_diag,*) subname,' min, max, sum = ', amin, amax, asum, trim(varname)
      endif

#else
      call abort_ice(subname//' ERROR: USE_NETCDF cpp not defined', &
         file=__FILE__, line=__LINE__)
      work_g = c0 ! to satisfy intent(out) attribute
#endif

      end subroutine ice_read_global_nc

!=======================================================================

! Report a netcdf error
! author: T. Craig

      subroutine ice_check_nc(status, abort_msg, file, line)
      integer(kind=int_kind), intent (in) :: status
      character (len=*)     , intent (in) :: abort_msg
      character (len=*)     , intent (in), optional :: file
      integer(kind=int_kind), intent (in), optional :: line

      ! local variables

      character(len=*), parameter :: subname = '(ice_check_nc)'

#ifdef USE_NETCDF
      if (status /= nf90_noerr) then
         if (present(file) .and. present(line)) then
            call abort_ice(subname//' '//trim(nf90_strerror(status))//', '//trim(abort_msg), &
                           file=file, line=line)
         elseif (present(file)) then
            call abort_ice(subname//' '//trim(nf90_strerror(status))//', '//trim(abort_msg), &
                           file=file)
         else
            call abort_ice(subname//' '//trim(nf90_strerror(status))//', '//trim(abort_msg))
         endif
      endif
#else
      call abort_ice(subname//' ERROR: USE_NETCDF cpp not defined', &
         file=__FILE__, line=__LINE__)
#endif

      end subroutine ice_check_nc

!=======================================================================

! Closes a netCDF file
! author: Alison McLaren, Met Office

      subroutine ice_close_nc(fid)

      integer (kind=int_kind), intent(in) :: &
         fid          ! unit number

      ! local variables

      character(len=*), parameter :: subname = '(ice_close_nc)'

#ifdef USE_NETCDF
      integer (kind=int_kind) :: &
        status        ! status variable from netCDF routine

      if (my_task == master_task) then
         status = nf90_close(fid)
         call ice_check_nc(status, subname//' ERROR: Cannot close file ', &
                           file=__FILE__, line=__LINE__ )
      endif
#else
      call abort_ice(subname//' ERROR: USE_NETCDF cpp not defined', &
         file=__FILE__, line=__LINE__)
#endif

      end subroutine ice_close_nc

!=======================================================================

! Read a netCDF file and scatter to processors.
! If the optional variables field_loc and field_type are present,
! the ghost cells are filled using values from the global array.
! This prevents them from being filled with zeroes in land cells
! (subroutine ice_HaloUpdate need not be called).
!
! Adapted by Elizabeth Hunke for reading 3D ocean currents

      subroutine ice_read_nc_uv(fid,  nrec, nzlev,  varname, work,  diag, &
                             field_loc, field_type, restart_ext)

      use ice_gather_scatter, only: scatter_global, scatter_global_ext

      integer (kind=int_kind), intent(in) :: &
         fid           , & ! file id
         nrec          , & ! record number
         nzlev             ! z level

      logical (kind=log_kind), intent(in) :: &
         diag              ! if true, write diagnostic output

      character (len=*), intent(in) :: &
         varname           ! field name in netcdf file

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks), intent(out) :: &
         work              ! output array (real, 8-byte)

      logical (kind=log_kind), optional, intent(in) :: &
         restart_ext       ! if true, read extended grid

      integer (kind=int_kind), optional, intent(in) :: &
         field_loc, &      ! location of field on staggered grid
         field_type        ! type of field (scalar, vector, angle)

      ! local variables

      character(len=*), parameter :: subname = '(ice_read_nc_uv)'

#ifdef USE_NETCDF
! netCDF file diagnostics:
      integer (kind=int_kind) :: &
         varid         , & ! variable id
         status            ! status output from netcdf routines
!        ndim, nvar    , & ! sizes of netcdf file
!        id,             & ! dimension index
!        dimlen            ! size of dimension

      real (kind=dbl_kind) :: &
         amin, amax, asum  ! min, max values and sum of input array

!     character (char_len) :: &
!        dimname           ! dimension name

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g1

      integer (kind=int_kind) :: nx, ny

      logical (kind=log_kind) :: lrestart_ext  ! local value of restart_ext

      lrestart_ext = .false.
      if (present(restart_ext)) then
         lrestart_ext = restart_ext
      endif

      if (lrestart_ext) then
         nx = nx_global + 2*nghost
         ny = ny_global + 2*nghost
      else
         nx = nx_global
         ny = ny_global
      endif

      if (my_task == master_task) then
         allocate(work_g1(nx,ny))
      else
         allocate(work_g1(1,1))   ! to save memory
      endif

      if (my_task == master_task) then

         !-------------------------------------------------------------
         ! Find out ID of required variable
         !-------------------------------------------------------------

         status = nf90_inq_varid(fid, trim(varname), varid)
         call ice_check_nc(status, subname//' ERROR: Cannot find variable '//trim(varname), &
                           file=__FILE__, line=__LINE__)

         !--------------------------------------------------------------
         ! Read global array
         !--------------------------------------------------------------

         status = nf90_get_var( fid, varid, work_g1, &
                                start=(/1,1,nzlev,nrec/), count=(/nx,ny,1,1/))
         call ice_check_nc(status, subname//' ERROR: Cannot get variable '//trim(varname), &
                           file=__FILE__, line=__LINE__)

      endif                     ! my_task = master_task

      !-------------------------------------------------------------------
      ! optional diagnostics
      !-------------------------------------------------------------------

      if (my_task==master_task .and. diag) then
         amin = minval(work_g1)
         amax = maxval(work_g1, mask = work_g1 /= spval_dbl)
         asum = sum   (work_g1, mask = work_g1 /= spval_dbl)
         write(nu_diag,*) subname,' min, max, sum =', amin, amax, asum, trim(varname)
      endif

      !-------------------------------------------------------------------
      ! Scatter data to individual processors.
      ! NOTE: Ghost cells are not updated unless field_loc is present.
      !-------------------------------------------------------------------

      if (lrestart_ext) then
         call scatter_global_ext(work, work_g1, master_task, distrb_info)
      else
         if (present(field_loc)) then
            call scatter_global(work, work_g1, master_task, distrb_info, &
                                field_loc, field_type)
         else
            call scatter_global(work, work_g1, master_task, distrb_info, &
                                field_loc_noupdate, field_type_noupdate)
         endif
      endif

      deallocate(work_g1)

#else
      call abort_ice(subname//' ERROR: USE_NETCDF cpp not defined', &
         file=__FILE__, line=__LINE__)
      work = c0 ! to satisfy intent(out) attribute
#endif

      end subroutine ice_read_nc_uv

!=======================================================================
! Read a vector in a netcdf file.
! Just like ice_read_global_nc except that it returns a vector.
! work_g is a real vector
!
! Adapted by William Lipscomb, LANL, from ice_read
! Adapted by Ann Keen, Met Office, to read from a netcdf file

      subroutine ice_read_vec_nc (fid,  nrec, varname, work_g, diag)

      integer (kind=int_kind), intent(in) :: &
         fid           , & ! file id
         nrec              ! record number

      character (char_len), intent(in) :: &
         varname           ! field name in netcdf file

      real (kind=dbl_kind), dimension(nrec), intent(out) :: &
         work_g            ! output array (real, 8-byte)

      logical (kind=log_kind) :: &
         diag              ! if true, write diagnostic output

      ! local variables

      character(len=*), parameter :: subname = '(ice_read_vec_nc)'

#ifdef USE_NETCDF
! netCDF file diagnostics:
      integer (kind=int_kind) :: &
         varid,          & ! netcdf id for field
         status            ! status output from netcdf routines

      real (kind=dbl_kind) :: &
         amin, amax        ! min, max values of input vector

      work_g(:) = c0

      if (my_task == master_task) then

         !-------------------------------------------------------------
         ! Find out ID of required variable
         !-------------------------------------------------------------

         status = nf90_inq_varid(fid, trim(varname), varid)
         call ice_check_nc(status, subname//' ERROR: Cannot find variable '//trim(varname), &
                           file=__FILE__, line=__LINE__)

         !--------------------------------------------------------------
         ! Read global array
         !--------------------------------------------------------------

         status = nf90_get_var( fid, varid, work_g, &
                                start=(/1/), count=(/nrec/))
         call ice_check_nc(status, subname//' ERROR: Cannot get variable '//trim(varname), &
                           file=__FILE__, line=__LINE__)

      endif                     ! my_task = master_task

      !-------------------------------------------------------------------
      ! optional diagnostics
      !-------------------------------------------------------------------

      if (my_task == master_task .and. diag) then
         amin = minval(work_g)
         amax = maxval(work_g)
         write(nu_diag,*) subname,' min, max, nrec = ', amin, amax, nrec
      endif

#else
      call abort_ice(subname//' ERROR: USE_NETCDF cpp not defined', &
         file=__FILE__, line=__LINE__)
      work_g = c0 ! to satisfy intent(out) attribute
#endif

      end subroutine ice_read_vec_nc

!=======================================================================
! Get number of variables of a given variable
      subroutine ice_get_ncvarsize(fid,varname,recsize)

      integer (kind=int_kind), intent(in) :: &
         fid                 ! file id
      character (char_len), intent(in) :: &
         varname             ! field name in netcdf file
      integer (kind=int_kind), intent(out) :: &
         recsize             ! Number of records in file

      ! local variables

#ifdef USE_NETCDF
      integer (kind=int_kind) :: &
         ndims, i, status
      character (char_len) :: &
         cvar
#endif
      character(len=*), parameter :: subname = '(ice_get_ncvarsize)'

#ifdef USE_NETCDF
      if (my_task ==  master_task) then
         status=nf90_inquire(fid, nDimensions = nDims)
         call ice_check_nc(status, subname//' ERROR: inquire nDimensions', &
                           file=__FILE__, line=__LINE__ )
         do i=1,nDims
            status = nf90_inquire_dimension(fid,i,name=cvar,len=recsize)
            call ice_check_nc(status, subname//' ERROR: inquire len for variable '//trim(cvar), &
                              file=__FILE__, line=__LINE__)
            if (trim(cvar) == trim(varname)) exit
         enddo
         if (trim(cvar) .ne. trim(varname)) then
            call abort_ice(subname//' ERROR: Did not find variable '//trim(varname), &
                           file=__FILE__, line=__LINE__)
         endif
      endif
#else
      call abort_ice(subname//' ERROR: USE_NETCDF cpp not defined', &
                     file=__FILE__, line=__LINE__)
      recsize = 0 ! to satisfy intent(out) attribute
#endif

      end subroutine ice_get_ncvarsize

!=======================================================================

      end module ice_read_write

!=======================================================================
