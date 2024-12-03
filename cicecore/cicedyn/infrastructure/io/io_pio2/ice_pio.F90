!============================================================================
!  Writes netcdf files
!    Created by Mariana Vertenstein, June 2009

  module ice_pio

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
  public ice_pio_finalize
  public ice_pio_check

#ifdef CESMCOUPLED
  type(iosystem_desc_t), pointer :: ice_pio_subsystem
#else
  type(iosystem_desc_t)          :: ice_pio_subsystem
#endif

!===============================================================================

  contains

!===============================================================================

!    Initialize the io subsystem
!    2009-Feb-17 - J. Edwards - initial version

   subroutine ice_pio_init(mode, filename, File, clobber, fformat, &
                           rearr, iotasks, root, stride, debug)

#ifdef CESMCOUPLED
   use shr_pio_mod, only: shr_pio_getiosys, shr_pio_getiotype, shr_pio_getioformat
#else
#ifdef GPTL
   use perf_mod, only : t_initf
#endif
#endif

   implicit none
   character(len=*)     , intent(in),    optional :: mode
   character(len=*)     , intent(in),    optional :: filename
   type(file_desc_t)    , intent(inout), optional :: File
   logical              , intent(in),    optional :: clobber
   character(len=*)     , intent(in),    optional :: fformat
   character(len=*)     , intent(in),    optional :: rearr
   integer              , intent(in),    optional :: iotasks
   integer              , intent(in),    optional :: root
   integer              , intent(in),    optional :: stride
   logical              , intent(in),    optional :: debug

   ! local variables

   integer (int_kind) :: &
      nml_error          ! namelist read error flag

   integer :: nprocs , lstride, lroot, liotasks, rearranger
   integer :: pio_iotype, status, nmode0, nmode
   logical :: lclobber, exists, ldebug
   character(len=*), parameter :: subname = '(ice_pio_init)'

#ifdef CESMCOUPLED
   ice_pio_subsystem => shr_pio_getiosys(inst_name)
   pio_iotype =  shr_pio_getiotype(inst_name)
   if ((pio_iotype==PIO_IOTYPE_NETCDF).or.(pio_iotype==PIO_IOTYPE_PNETCDF)) then
      nmode0 = shr_pio_getioformat(inst_name)
   else
      nmode0 = 0
   endif

   call pio_seterrorhandling(ice_pio_subsystem, PIO_RETURN_ERROR)
#else

#ifdef GPTL
   !--- initialize gptl
   call t_initf('undefined_NLFileName', LogPrint=.false., mpicom=MPI_COMM_ICE, &
         MasterTask=.true.)
#endif

   !--- initialize type of io
   ldebug = .false.
   if (present(debug)) then
      ldebug = debug
   endif

   if (present(fformat)) then
      if (fformat(1:3) == 'cdf') then
         pio_iotype = PIO_IOTYPE_NETCDF
      elseif (fformat(1:3) == 'hdf') then
         pio_iotype = PIO_IOTYPE_NETCDF4P
      elseif (fformat(1:7) == 'pnetcdf') then
         pio_iotype = PIO_IOTYPE_PNETCDF
      else
         call abort_ice(subname//' ERROR: format not allowed for '//trim(fformat), &
            file=__FILE__, line=__LINE__)
      endif

      if (fformat == 'cdf2' .or. fformat == 'pnetcdf2') then
         nmode0 = PIO_64BIT_OFFSET
      elseif (fformat == 'cdf5' .or. fformat == 'pnetcdf5') then
         nmode0 = PIO_64BIT_DATA
      else
         nmode0 = 0
      endif
   else
      pio_iotype = PIO_IOTYPE_NETCDF
      nmode0 = 0
   endif

   if (present(rearr)) then
      if (rearr == 'box' .or. rearr == 'default') then
         rearranger = PIO_REARR_BOX
      elseif (rearr == 'subset') then
         rearranger = PIO_REARR_SUBSET
      else
         call abort_ice(subname//' ERROR: rearr not allowed for '//trim(rearr), &
            file=__FILE__, line=__LINE__)
      endif
   else
      rearranger = PIO_REARR_BOX
   endif

   nprocs = get_num_procs()
   lstride = 4
   lroot = min(1,nprocs-1)
!  Adjustments for PIO2 iotask issue, https://github.com/NCAR/ParallelIO/issues/1986
!   liotasks = max(1,(nprocs-lroot)/lstride)  ! very conservative
   liotasks = max(1,nprocs/lstride - lroot/lstride)  ! less conservative (note integer math)
!   liotasks = 1 + (nprocs-lroot-1)/lstride   ! optimal

   if (present(iotasks)) then
      if (iotasks /= -99) liotasks=iotasks
   endif
   if (present(root)) then
      if (root /= -99) lroot=root
   endif
   if (present(stride)) then
      if (stride /= -99) lstride=stride
   endif

   if (liotasks < 1 .or. lroot < 0 .or. lstride < 1) then
      call abort_ice(subname//' ERROR: iotasks, root, stride incorrect ', &
         file=__FILE__, line=__LINE__)
   endif

   ! adjust to fit in nprocs, preserve root and stride as much as possible
   lroot = min(lroot,nprocs-1)   ! lroot <= nprocs-1
!  Adjustments for PIO2 iotask issue, https://github.com/NCAR/ParallelIO/issues/1986
!   liotasks = max(1,min(liotasks, (nprocs-lroot)/lstride))  ! very conservative
   liotasks = max(1,min(liotasks,nprocs/lstride - lroot/lstride))  ! less conservative (note integer math)
!   liotasks = max(1,min(liotasks, 1 + (nprocs-lroot-1)/lstride))  ! optimal

   !--- initialize ice_pio_subsystem

   if (ldebug .and. my_task == master_task) then
      write(nu_diag,*) subname,' nprocs     = ',nprocs
      write(nu_diag,*) subname,' pio_iotype = ',pio_iotype
      write(nu_diag,*) subname,' iotasks    = ',liotasks
      write(nu_diag,*) subname,' baseroot   = ',lroot
      write(nu_diag,*) subname,' stride     = ',lstride
      write(nu_diag,*) subname,' nmode      = ',nmode0
   end if

   call pio_init(my_task, MPI_COMM_ICE, liotasks, master_task, lstride, &
                 rearranger, ice_pio_subsystem, base=lroot)

   call pio_seterrorhandling(ice_pio_subsystem, PIO_RETURN_ERROR)

#endif

   if (present(mode) .and. present(filename) .and. present(File)) then

      if (trim(mode) == 'write') then

         lclobber = .false.
         if (present(clobber)) then
            lclobber=clobber
         endif

         if (File%fh<0) then
            ! filename not open
            inquire(file=trim(filename),exist=exists)
            if (exists) then
               if (lclobber) then
                  nmode = ior(PIO_CLOBBER,nmode0)
                  status = pio_createfile(ice_pio_subsystem, File, pio_iotype, trim(filename), nmode)
                  call ice_pio_check(status, subname//' ERROR: Failed to overwrite file '//trim(filename), &
                       file=__FILE__,line=__LINE__)
                  if (my_task == master_task) then
                     write(nu_diag,*) subname,' create file ',trim(filename)
                  end if
               else
                  nmode = pio_write
                  status = pio_openfile(ice_pio_subsystem, File, pio_iotype, trim(filename), nmode)
                  call ice_pio_check( status,  subname//' ERROR: Failed to open file '//trim(filename), &
                       file=__FILE__,line=__LINE__)
                  if (my_task == master_task) then
                     write(nu_diag,*) subname,' open file ',trim(filename)
                  end if
               endif
            else
               nmode = ior(PIO_NOCLOBBER,nmode0)
               status = pio_createfile(ice_pio_subsystem, File, pio_iotype, trim(filename), nmode)
               call ice_pio_check( status, subname//' ERROR: Failed to create file '//trim(filename), &
                    file=__FILE__,line=__LINE__)
               if (my_task == master_task) then
                  write(nu_diag,*) subname,' create file ',trim(filename)
               end if
            endif
         ! else: filename is already open, just return
         endif
      end if

      if (trim(mode) == 'read') then
         inquire(file=trim(filename),exist=exists)
         if (exists) then
            if (my_task == master_task) then
               write(nu_diag,*) subname//' opening file for reading '//trim(filename)
            endif
            status = pio_openfile(ice_pio_subsystem, File, pio_iotype, trim(filename), pio_nowrite)
            if (status /= PIO_NOERR) then
               if (my_task == master_task) then
                  write(nu_diag,*) subname//' opening '//trim(filename)//' as type '//trim(fformat)//' failed, retrying as type cdf1'
               endif
               status = pio_openfile(ice_pio_subsystem, File, PIO_IOTYPE_NETCDF, trim(filename), pio_nowrite)
               call ice_pio_check( status, subname//' ERROR: Failed to open file '//trim(filename), &
               file=__FILE__,line=__LINE__)
            endif
         else
            if(my_task==master_task) then
               write(nu_diag,*) subname//' ERROR: file not found '//trim(filename)
            end if
            call abort_ice(subname//' ERROR: aborting with invalid file '//trim(filename))
         endif
      end if

   end if

   call pio_seterrorhandling(ice_pio_subsystem, PIO_INTERNAL_ERROR)

   end subroutine ice_pio_init

!================================================================================

   subroutine ice_pio_initdecomp_2d(iodesc, precision)

      type(io_desc_t), intent(out) :: iodesc
      integer(kind=int_kind), optional, intent(in) :: precision

      integer (kind=int_kind) :: &
          iblk,ilo,ihi,jlo,jhi,lon,lat,i,j,n,k

      type(block) :: this_block

      integer(kind=int_kind), pointer :: dof2d(:)
      integer(kind=int_kind) :: lprecision
      character(len=*), parameter :: subname = '(ice_pio_initdecomp_2d)'

      lprecision = 8
      if (present(precision)) lprecision = precision

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

      if (lprecision == 8) then
         call pio_initdecomp(ice_pio_subsystem, pio_double, (/nx_global,ny_global/), &
              dof2d, iodesc)
      else
         call pio_initdecomp(ice_pio_subsystem, pio_real, (/nx_global,ny_global/), &
              dof2d, iodesc)
      endif

      deallocate(dof2d)

   end subroutine ice_pio_initdecomp_2d

!================================================================================

   subroutine ice_pio_initdecomp_3d (ndim3, iodesc, remap, precision)

      integer(kind=int_kind), intent(in) :: ndim3
      type(io_desc_t), intent(out) :: iodesc
      logical, optional :: remap
      integer(kind=int_kind), optional, intent(in) :: precision
      integer (kind=int_kind) :: &
          iblk,ilo,ihi,jlo,jhi,lon,lat,i,j,n,k

      type(block) :: this_block
      logical :: lremap
      integer(kind=int_kind), pointer :: dof3d(:)
      integer(kind=int_kind) :: lprecision
      character(len=*), parameter :: subname = '(ice_pio_initdecomp_2d)'

      lprecision = 8
      if (present(precision)) lprecision = precision

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

      if (lprecision == 8) then
         call pio_initdecomp(ice_pio_subsystem, pio_double, (/nx_global,ny_global,ndim3/), &
              dof3d, iodesc)
      else
         call pio_initdecomp(ice_pio_subsystem, pio_real, (/nx_global,ny_global,ndim3/), &
              dof3d, iodesc)
      endif

      deallocate(dof3d)

   end subroutine ice_pio_initdecomp_3d

!================================================================================

   subroutine ice_pio_initdecomp_3d_inner(ndim3, inner_dim, iodesc, precision)

      integer(kind=int_kind), intent(in) :: ndim3
      logical, intent(in) :: inner_dim
      type(io_desc_t), intent(out) :: iodesc
      integer(kind=int_kind), optional, intent(in) :: precision

      integer (kind=int_kind) :: &
          iblk,ilo,ihi,jlo,jhi,lon,lat,i,j,n,k

      type(block) :: this_block

      integer(kind=int_kind), pointer :: dof3d(:)
      integer(kind=int_kind) :: lprecision
      character(len=*), parameter :: subname = '(ice_pio_initdecomp_3d_inner)'

      lprecision = 8
      if (present(precision)) lprecision = precision

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

      if (lprecision == 8) then
         call pio_initdecomp(ice_pio_subsystem, pio_double, (/ndim3,nx_global,ny_global/), &
              dof3d, iodesc)
      else
         call pio_initdecomp(ice_pio_subsystem, pio_real, (/ndim3,nx_global,ny_global/), &
              dof3d, iodesc)
      endif

      deallocate(dof3d)

   end subroutine ice_pio_initdecomp_3d_inner

!================================================================================

   subroutine ice_pio_initdecomp_4d (ndim3, ndim4, iodesc, precision)

      integer(kind=int_kind), intent(in) :: ndim3, ndim4
      type(io_desc_t), intent(out) :: iodesc
      integer(kind=int_kind), optional, intent(in) :: precision

      integer (kind=int_kind) :: &
          iblk,ilo,ihi,jlo,jhi,lon,lat,i,j,n,k,l

      type(block) :: this_block

      integer(kind=int_kind), pointer :: dof4d(:)
      integer(kind=int_kind) :: lprecision
      character(len=*), parameter :: subname = '(ice_pio_initdecomp_4d)'

      lprecision = 8
      if (present(precision)) lprecision = precision

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

      if (lprecision == 8) then
         call pio_initdecomp(ice_pio_subsystem, pio_double, &
             (/nx_global,ny_global,ndim3,ndim4/), dof4d, iodesc)
      else
         call pio_initdecomp(ice_pio_subsystem, pio_real, &
             (/nx_global,ny_global,ndim3,ndim4/), dof4d, iodesc)
      endif

      deallocate(dof4d)

   end subroutine ice_pio_initdecomp_4d


!================================================================================

   ! PIO Finalize

   subroutine ice_pio_finalize()

      integer(kind=int_kind)      :: status
      character(len=*), parameter :: subname = '(ice_pio_finalize)'

      status = PIO_NOERR
#ifndef CESMCOUPLED
      call pio_seterrorhandling(ice_pio_subsystem, PIO_RETURN_ERROR)
      call pio_finalize(ice_pio_subsystem,status)
      call ice_pio_check( status, subname//' ERROR: Failed to finalize ', &
         file=__FILE__,line=__LINE__)
! do not call this, ice_pio_subsystem does not exist anymore
!      call pio_seterrorhandling(ice_pio_subsystem, PIO_INTERNAL_ERROR)
#endif

   end subroutine ice_pio_finalize

!================================================================================

   ! PIO Error handling
   ! Author: Anton Steketee, ACCESS-NRI

   subroutine ice_pio_check(status, abort_msg, file, line)
      integer(kind=int_kind), intent (in) :: status
      character (len=*)     , intent (in) :: abort_msg
      character (len=*)     , intent (in), optional :: file
      integer(kind=int_kind), intent (in), optional :: line

      ! local variables

      character(len=pio_max_name) :: err_msg
      integer(kind=int_kind)      :: strerror_status
      character(len=*), parameter :: subname = '(ice_pio_check)'

      if (status /= PIO_NOERR) then
#ifdef USE_PIO1
         err_msg = ''
#else
         strerror_status = pio_strerror(status, err_msg)
#endif
         if (present(file) .and. present(line)) then
            call abort_ice(subname//trim(err_msg)//', '//trim(abort_msg), file=file, line=line)
         elseif (present(file)) then
            call abort_ice(subname//trim(err_msg)//', '//trim(abort_msg), file=file)
         else
            call abort_ice(subname//trim(err_msg)//', '//trim(abort_msg))
         endif
      endif
   end subroutine ice_pio_check

!================================================================================

  end module ice_pio

!================================================================================
