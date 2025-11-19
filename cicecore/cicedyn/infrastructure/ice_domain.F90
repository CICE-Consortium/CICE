#ifdef ncdf
#define USE_NETCDF
#endif
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module ice_domain

!  This module contains the model domain and routines for initializing
!  the domain.  It also initializes the decompositions and
!  distributions across processors/threads by calling relevant
!  routines in the block, distribution modules.
!
! author: Phil Jones, LANL
! Oct. 2004: Adapted from POP by William H. Lipscomb, LANL
! Feb. 2007: E. Hunke removed NE and SW boundary options (they were buggy
!  and not used anyhow).

   use ice_kinds_mod
   use ice_constants, only: shlat, nhlat
   use ice_communicate, only: my_task, master_task, get_num_procs, &
       add_mpi_barriers, ice_barrier
   use ice_broadcast, only: broadcast_scalar, broadcast_array
   use ice_blocks, only: block, get_block, create_blocks, nghost, &
       nblocks_x, nblocks_y, nblocks_tot, nx_block, ny_block, debug_blocks
   use ice_distribution, only: distrb
   use ice_boundary, only: ice_halo
   use ice_exit, only: abort_ice
   use ice_fileunits, only: nu_nml, nml_filename, nu_diag, &
       get_fileunit, release_fileunit, flush_fileunit
   use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
   use icepack_intfc, only: icepack_query_parameters

#ifdef USE_NETCDF
   use netcdf
#endif

   implicit none
   private

   public  :: init_domain_blocks ,&
              init_domain_distribution

   integer (int_kind), public :: &
      nblocks         ! actual number of blocks on this processor

   logical (kind=log_kind), public :: &
      close_boundaries   ! deprecated Nov, 2025

   integer (int_kind), dimension(:), pointer, public :: &
      blocks_ice => null()        ! block ids for local blocks

   type (distrb), public :: &
      distrb_info        ! block distribution info

   type (ice_halo), public :: &
      halo_info          ! ghost cell update info

   character (char_len), public :: &
      ew_boundary_type, &! type of domain bndy in each logical
      ns_boundary_type   !    direction (ew is i, ns is j)

   logical (kind=log_kind), public :: &
      maskhalo_dyn   , & ! if true, use masked halo updates for dynamics
      maskhalo_remap , & ! if true, use masked halo updates for transport
      maskhalo_bound , & ! if true, use masked halo updates for bound_state
      halo_dynbundle , & ! if true, bundle halo update in dynamics
      landblockelim      ! if true, land block elimination is on

!-----------------------------------------------------------------------
!
!   module private variables - for the most part these appear as
!   module variables to facilitate sharing info between init_domain1
!   and init_domain2.
!
!-----------------------------------------------------------------------

    character (char_len) :: &
       distribution_type,   &! method to use for distributing blocks
                             ! 'cartesian', 'roundrobin', 'sectrobin', 'sectcart'
                             ! 'rake', 'spacecurve', etc
       distribution_wght     ! method for weighting work per block
                             ! 'block' = block weighted method with land block elimination
                             ! 'blockall' = block method with NO land block elimination and minimum weight given to land blocks
                             ! 'blockfull'= block method with NO land block elimination and full weight given to land blocks
                             ! 'latitude' = no. ocean points * |lat|
                             ! 'file' = read distribution_wgth_file
    character (char_len_long) :: &
       distribution_wght_file  ! file for distribution_wght=file

    integer (int_kind) :: &
       nprocs                ! num of processors

!***********************************************************************

 contains

!***********************************************************************

 subroutine init_domain_blocks

!  This routine reads in domain information and calls the routine
!  to set up the block decomposition.

   use ice_distribution, only: processor_shape, proc_decomposition
   use ice_domain_size, only: ncat, nilyr, nslyr, max_blocks, &
       nx_global, ny_global, block_size_x, block_size_y
   use ice_fileunits, only: goto_nml

!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind) :: &
      nml_error, &          ! namelist read error flag
      nprocs_x, nprocs_y    ! procs decomposed into blocks

   character(len=char_len)      :: nml_name ! text namelist name
   character(len=char_len_long) :: tmpstr2 ! for namelist check

   character(len=*), parameter :: subname = '(init_domain_blocks)'

!----------------------------------------------------------------------
!
!  input namelists
!
!----------------------------------------------------------------------

   namelist /domain_nml/ nprocs, &
                         max_blocks,   &
                         block_size_x, &
                         block_size_y, &
                         nx_global,    &
                         ny_global,    &
                         processor_shape,   &
                         distribution_type, &
                         distribution_wght, &
                         distribution_wght_file, &
                         ew_boundary_type,  &
                         ns_boundary_type,  &
                         maskhalo_dyn,      &
                         maskhalo_remap,    &
                         maskhalo_bound,    &
                         add_mpi_barriers,  &
                         debug_blocks

!----------------------------------------------------------------------
!
!  read domain information from namelist input
!
!----------------------------------------------------------------------

   nprocs = -1
   processor_shape   = 'slenderX2'
   distribution_type = 'cartesian'
   distribution_wght = 'latitude'
   distribution_wght_file = 'unknown'
   ew_boundary_type  = 'cyclic'
   ns_boundary_type  = 'open'
   maskhalo_dyn      = .false.     ! if true, use masked halos for dynamics
   maskhalo_remap    = .false.     ! if true, use masked halos for transport
   maskhalo_bound    = .false.     ! if true, use masked halos for bound_state
   halo_dynbundle    = .true.      ! if true, bundle halo updates in dynamics
   add_mpi_barriers  = .false.     ! if true, throttle communication
   debug_blocks      = .false.     ! if true, print verbose block information
   max_blocks        = -1          ! max number of blocks per processor
   block_size_x      = -1          ! size of block in first horiz dimension
   block_size_y      = -1          ! size of block in second horiz dimension
   nx_global         = -1          ! NXGLOB,  i-axis size
   ny_global         = -1          ! NYGLOB,  j-axis size
   landblockelim     = .true.      ! on by default

   if (my_task == master_task) then
      nml_name = 'domain_nml'
      write(nu_diag,*) subname,' Reading ', trim(nml_name)

      call get_fileunit(nu_nml)
      open (nu_nml, file=trim(nml_filename), status='old',iostat=nml_error)
      if (nml_error /= 0) then
         call abort_ice(subname//' ERROR: domain_nml open file '// &
              trim(nml_filename), file=__FILE__, line=__LINE__)
      endif

      call goto_nml(nu_nml,trim(nml_name),nml_error)
      if (nml_error /= 0) then
         call abort_ice(subname//' ERROR: searching for '// trim(nml_name), &
              file=__FILE__, line=__LINE__)
      endif

      nml_error =  1
      do while (nml_error > 0)
         read(nu_nml, nml=domain_nml,iostat=nml_error)
         ! check if error
         if (nml_error /= 0) then
            ! backspace and re-read erroneous line
            backspace(nu_nml)
            read(nu_nml,fmt='(A)') tmpstr2
            call abort_ice(subname//' ERROR: ' // trim(nml_name) // ' reading ' // &
                 trim(tmpstr2), file=__FILE__, line=__LINE__)
         endif
      end do

      close(nu_nml)
      call release_fileunit(nu_nml)

   endif

   call broadcast_scalar(nprocs,            master_task)
   call broadcast_scalar(processor_shape,   master_task)
   call broadcast_scalar(distribution_type, master_task)
   call broadcast_scalar(distribution_wght, master_task)
   call broadcast_scalar(distribution_wght_file, master_task)
   call broadcast_scalar(ew_boundary_type,  master_task)
   call broadcast_scalar(ns_boundary_type,  master_task)
   call broadcast_scalar(maskhalo_dyn,      master_task)
   call broadcast_scalar(maskhalo_remap,    master_task)
   call broadcast_scalar(maskhalo_bound,    master_task)
   call broadcast_scalar(add_mpi_barriers,  master_task)
   call broadcast_scalar(debug_blocks,      master_task)
   call broadcast_scalar(max_blocks,        master_task)
   call broadcast_scalar(block_size_x,      master_task)
   call broadcast_scalar(block_size_y,      master_task)
   call broadcast_scalar(nx_global,         master_task)
   call broadcast_scalar(ny_global,         master_task)

!----------------------------------------------------------------------
!
! Set nprocs if not explicitly set to valid value in namelist
!
!----------------------------------------------------------------------

#ifdef CESMCOUPLED
   nprocs = get_num_procs()
#else
   if (nprocs < 0) then
      nprocs = get_num_procs()
   else if (nprocs /= get_num_procs()) then
      write(nu_diag,*) subname,' ERROR: nprocs, get_num_procs = ',nprocs,get_num_procs()
      call abort_ice(subname//' ERROR: Input nprocs not same as system (e.g MPI) request', file=__FILE__, line=__LINE__)
   endif
#endif

!----------------------------------------------------------------------
!
!  perform some basic checks on domain
!
!----------------------------------------------------------------------

   if (nx_global < 1 .or. ny_global < 1 .or. ncat < 1) then
      !***
      !*** domain size zero or negative
      !***
      call abort_ice(subname//' ERROR: Invalid domain: size < 1', file=__FILE__, line=__LINE__) ! no domain
   else if (nghost < 1) then
      !***
      !*** must have at least 1 layer of ghost cells
      !***
      call abort_ice(subname//' ERROR: Not enough ghost cells allocated', file=__FILE__, line=__LINE__)
   endif

!----------------------------------------------------------------------
!
!  compute block decomposition and details
!
!----------------------------------------------------------------------

   call create_blocks(nx_global, ny_global, trim(ew_boundary_type), &
                                            trim(ns_boundary_type))

!----------------------------------------------------------------------
!
!  Now we need grid info before proceeding further
!  Print some domain information
!
!----------------------------------------------------------------------

   if (my_task == master_task) then
     write(nu_diag,'(/,a18,/)')'Domain Information'
     write(nu_diag,'(a,i6)')  '  Horizontal domain: nx = ', nx_global
     write(nu_diag,'(a,i6)')  '                     ny = ', ny_global
     write(nu_diag,'(a,i6)')  '  No. of categories: nc = ', ncat
     write(nu_diag,'(a,i6)')  '  No. of ice layers: ni = ', nilyr
     write(nu_diag,'(a,i6)')  '  No. of snow layers:ns = ', nslyr
     write(nu_diag,'(a,i6)')  '  Processors:  total    = ', nprocs
     write(nu_diag,'(a,a)')   '  Processor shape       = ', trim(processor_shape)
     write(nu_diag,'(a,a)')   '  Distribution type     = ', trim(distribution_type)
     write(nu_diag,'(a,a)')   '  Distribution weight   = ', trim(distribution_wght)
     write(nu_diag,'(a,a)')   '  Distribution wght file= ', trim(distribution_wght_file)
     write(nu_diag,'(a,a)')   '  ew_boundary_type      = ', trim(ew_boundary_type)
     write(nu_diag,'(a,a)')   '  ns_boundary_type      = ', trim(ns_boundary_type)
     write(nu_diag,'(a,l6)')  '  maskhalo_dyn          = ', maskhalo_dyn
     write(nu_diag,'(a,l6)')  '  maskhalo_remap        = ', maskhalo_remap
     write(nu_diag,'(a,l6)')  '  maskhalo_bound        = ', maskhalo_bound
     write(nu_diag,'(a,l6)')  '  add_mpi_barriers      = ', add_mpi_barriers
     write(nu_diag,'(a,l6)')  '  debug_blocks          = ', debug_blocks
     write(nu_diag,'(a,2i6)') '  block_size_x,_y       = ', block_size_x, block_size_y
     write(nu_diag,'(a,i6)')  '  max_blocks            = ', max_blocks
     write(nu_diag,'(a,i6,/)')'  Number of ghost cells = ', nghost
   endif

!----------------------------------------------------------------------

 end subroutine init_domain_blocks

!***********************************************************************

 subroutine init_domain_distribution(KMTG,ULATG,grid_ice)

!  This routine calls appropriate setup routines to distribute blocks
!  across processors and defines arrays with block ids for any local
!  blocks. Information about ghost cell update routines is also
!  initialized here through calls to the appropriate boundary routines.

   use ice_boundary, only: ice_HaloCreate
   use ice_distribution, only: create_distribution, create_local_block_ids, ice_distributionGet
   use ice_domain_size, only: max_blocks, nx_global, ny_global
   use ice_global_reductions, only: global_sum, global_maxval

   real (dbl_kind), dimension(nx_global,ny_global), intent(in) :: &
      KMTG           ,&! global topography
      ULATG            ! global latitude field (radians)

   character(len=*), intent(in) :: &
      grid_ice         ! grid_ice, B, C, CD, etc

!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (int_kind), dimension (nx_global, ny_global) :: &
      flat                 ! latitude-dependent scaling factor

   character (char_len) :: outstring

   integer (int_kind), parameter :: &
      max_work_unit=10    ! quantize the work into values from 1,max

   integer (int_kind) :: &
      i,j,n              ,&! dummy loop indices
      ig,jg              ,&! global indices
      igm1,igp1,jgm1,jgp1,&! global indices
      ninfo              ,&! ice_distributionGet check
      np, nlb, m         ,&! debug blocks temporaries
      work_unit          ,&! size of quantized work unit
#ifdef USE_NETCDF
      fid                ,&! file id
      varid              ,&! var id
      status             ,&! netcdf return code
#endif
      tblocks_tmp        ,&! total number of blocks
      nblocks_tmp        ,&! temporary value of nblocks
      nblocks_max          ! max blocks on proc

   real (dbl_kind) :: &
      puny, &              ! puny limit
      rad_to_deg           ! radians to degrees

   integer (int_kind), dimension(:), allocatable :: &
      blkinfo            ,&! ice_distributionGet check
      nocn               ,&! number of ocean points per block
      work_per_block       ! number of work units per block

   type (block) :: &
      this_block           ! block information for current block

   real (dbl_kind), dimension(:,:), allocatable :: &
      wght                 ! wghts from file

   character(len=*), parameter :: subname = '(init_domain_distribution)'

   call icepack_query_parameters(puny_out=puny, rad_to_deg_out=rad_to_deg)
   call icepack_warnings_flush(nu_diag)
   if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
      file=__FILE__, line=__LINE__)

!----------------------------------------------------------------------
!
!  estimate the amount of work per processor using the topography
!  and latitude
!
!----------------------------------------------------------------------

   if (distribution_wght == 'latitude') then
       flat = max(NINT(abs(ULATG*rad_to_deg), int_kind),1)  ! linear function
   else
       flat = 1
   endif

   if (distribution_wght == 'blockall' ) landblockelim = .false.
   if (distribution_wght == 'blockfull') landblockelim = .false.

   allocate(nocn(nblocks_tot))

   if (distribution_wght == 'file') then
      allocate(wght(nx_global,ny_global))
      if (my_task == master_task) then
         ! cannot use ice_read_write due to circular dependency
#ifdef USE_NETCDF
         status = nf90_open(distribution_wght_file, NF90_NOWRITE, fid)
         if (status /= nf90_noerr) then
            call abort_ice(subname//' ERROR: Cannot open '//trim(distribution_wght_file), &
                           file=__FILE__, line=__LINE__)
         endif
         status = nf90_inq_varid(fid, 'wght', varid)
         if (status /= nf90_noerr) then
            call abort_ice(subname//' ERROR: Cannot find wght '//trim(distribution_wght_file), &
                           file=__FILE__, line=__LINE__)
         endif
         status = nf90_get_var(fid, varid, wght)
         if (status /= nf90_noerr) then
            call abort_ice(subname//' ERROR: Cannot get wght '//trim(distribution_wght_file), &
                           file=__FILE__, line=__LINE__)
         endif
         status = nf90_close(fid)
         if (status /= nf90_noerr) then
            call abort_ice(subname//' ERROR: Cannot close '//trim(distribution_wght_file), &
                           file=__FILE__, line=__LINE__)
         endif
         write(nu_diag,*) 'read ',trim(distribution_wght_file),minval(wght),maxval(wght)
#else
         call abort_ice(subname//' ERROR: USE_NETCDF cpp not defined', &
             file=__FILE__, line=__LINE__)
#endif
      endif
      call broadcast_array(wght, master_task)
      nocn = 0
      do n=1,nblocks_tot
         this_block = get_block(n,n)
         do j=this_block%jlo,this_block%jhi
            jg = this_block%j_glob(j)
            if (jg > 0) then
               do i=this_block%ilo,this_block%ihi
                  ig = this_block%i_glob(i)
                  if (ig > 0) then
!                     if (KMTG(ig,jg) > puny) &
!                        nocn(n) = max(nocn(n),nint(wght(ig,jg)+1.0_dbl_kind))
                     if (KMTG(ig,jg) > puny) then
                        if (wght(ig,jg) > 0.00001_dbl_kind) then
                           nocn(n) = nocn(n)+nint(wght(ig,jg))
                        else
                           nocn(n) = max(nocn(n),1)
                        endif
                     endif
                  endif
               end do
            endif
         end do
      enddo
      deallocate(wght)
   else
      nocn = 0
      do n=1,nblocks_tot
         this_block = get_block(n,n)
         do j=this_block%jlo,this_block%jhi
            jg = this_block%j_glob(j)
            if (jg > 0) then
               do i=this_block%ilo,this_block%ihi
                  ig = this_block%i_glob(i)
                  if (ig > 0) then
                     if (grid_ice == 'C' .or. grid_ice == 'CD') then
                        ! Have to be careful about block elimination with C/CD
                        ! Use a bigger stencil
                        igm1 = mod(ig-2+nx_global,nx_global)+1
                        igp1 = mod(ig,nx_global)+1
                        jgm1 = max(jg-1,1)
                        jgp1 = min(jg+1,ny_global)
                        if ((KMTG(ig  ,jg  ) > puny .or.                             &
                             KMTG(igm1,jg  ) > puny .or. KMTG(igp1,jg  ) > puny .or. &
                             KMTG(ig  ,jgp1) > puny .or. KMTG(ig  ,jgm1) > puny) .and. &
                            (ULATG(ig,jg) < shlat/rad_to_deg .or.          &
                             ULATG(ig,jg) > nhlat/rad_to_deg) )            &
                             nocn(n) = nocn(n) + flat(ig,jg)
                     else
                        if (KMTG(ig,jg) > puny .and.                      &
                           (ULATG(ig,jg) < shlat/rad_to_deg .or.          &
                            ULATG(ig,jg) > nhlat/rad_to_deg) )            &
                             nocn(n) = nocn(n) + flat(ig,jg)
                     endif
                  endif
               end do
            endif
         end do

         !*** with array syntax, we actually do work on non-ocean
         !*** points, so where the block is not completely land,
         !*** reset nocn to be the full size of the block

         ! use processor_shape = 'square-pop' and distribution_wght = 'block'
         ! to make CICE and POP decompositions/distributions identical.

#ifdef CICE_IN_NEMO
         ! Keep all blocks even the ones only containing land points
         ! tcraig, use 'blockfull', get rid of the CPP, keep for backwards compatibility for now
         if (distribution_wght == 'block') nocn(n) = nx_block*ny_block
#else
         if (distribution_wght == 'blockfull') nocn(n) = nx_block*ny_block
         if (distribution_wght == 'block' .and. nocn(n) > 0) nocn(n) = nx_block*ny_block
         if (.not. landblockelim) nocn(n) = max(nocn(n),1)
#endif
      end do
   endif  ! distribution_wght = file

   work_unit = maxval(nocn)/max_work_unit + 1

   !*** find number of work units per block

   allocate(work_per_block(nblocks_tot))

   where (nocn > 1)
      work_per_block = nocn/work_unit + 2
   elsewhere (nocn == 1)
      work_per_block = nocn/work_unit + 1
   elsewhere
      work_per_block = 0
   end where
   if (my_task == master_task) then
      write(nu_diag,'(2a,4i9)') subname,' work_unit      = ',work_unit, max_work_unit
      write(nu_diag,'(2a,4i9)') subname,' nocn           = ',minval(nocn),maxval(nocn),sum(nocn)
      write(nu_diag,'(2a,4i9)') subname,' work_per_block = ',minval(work_per_block),maxval(work_per_block),sum(work_per_block)
   endif
   deallocate(nocn)


!----------------------------------------------------------------------
!
!  determine the distribution of blocks across processors
!
!----------------------------------------------------------------------

   distrb_info = create_distribution(distribution_type, &
                                     nprocs, work_per_block)

   deallocate(work_per_block)

!----------------------------------------------------------------------
!
!  allocate and determine block id for any local blocks
!
!----------------------------------------------------------------------

   call create_local_block_ids(blocks_ice, distrb_info)

!----------------------------------------------------------------------
!
!  check block sizes and max_blocks
!
!----------------------------------------------------------------------

   if (associated(blocks_ice)) then
      nblocks = size(blocks_ice)
   else
      nblocks = 0
   endif

   tblocks_tmp = global_sum(nblocks, distrb_info)
   nblocks_max = global_maxval(nblocks, distrb_info)

   if (my_task == master_task) then
      write(nu_diag,'(2a,i8)') subname,' total number of blocks is', tblocks_tmp
   endif

   if (nblocks > max_blocks) then
      write(nu_diag,'(2a,2i6)') subname,' ERROR: nblocks, max_blocks = ',nblocks,max_blocks
      write(nu_diag,'(2a,2i6)') subname,' ERROR: max_blocks too small: increase to', nblocks_max
      call abort_ice(subname//' ERROR max_blocks too small', file=__FILE__, line=__LINE__)
   else if (nblocks_max < max_blocks) then
      if (my_task == master_task) then
          write(nu_diag,'(2a,2i6)') subname,' NOTE: max_blocks too large: decrease to', nblocks_max
      endif
   endif

!----------------------------------------------------------------------
!
! write out block distribution
! internal check of icedistributionGet as part of verification process
!
!----------------------------------------------------------------------

   if (debug_blocks) then

      call flush_fileunit(nu_diag)
      call ice_barrier()
      if (my_task == master_task) then
         write(nu_diag,*) ' '
         write(nu_diag,'(2a)') subname, ' Blocks by Proc:'
      endif
      call ice_distributionGet(distrb_info, nprocs=np, numLocalBlocks=nlb)
      do m = 1, np
         if (m == my_task+1) then
            do n=1,nlb
               write(nu_diag,'(2a,3i8)') &
                     subname,' my_task, local block ID, global block ID: ', &
                     my_task, n, distrb_info%blockGlobalID(n)
            enddo
            call flush_fileunit(nu_diag)
         endif
         call ice_barrier()
      enddo

      if (my_task == master_task) then
         write(nu_diag,*) ' '
         write(nu_diag,'(2a)') subname, ' Blocks by Global Block ID:'
         do m = 1, nblocks_tot
            write(nu_diag,'(2a,3i8)') &
                  subname,' global block id, proc, local block ID: ', &
                  m, distrb_info%blockLocation(m), distrb_info%blockLocalID(m)
         enddo
         call flush_fileunit(nu_diag)
      endif
      call ice_barrier()

      call ice_distributionGet(distrb_info, nprocs=ninfo)
      if (ninfo /= distrb_info%nprocs) &
         call abort_ice(subname//' ice_distributionGet nprocs ERROR', file=__FILE__, line=__LINE__)

      call ice_distributionGet(distrb_info, communicator=ninfo)
      if (ninfo /= distrb_info%communicator) &
         call abort_ice(subname//' ice_distributionGet communicator ERROR', file=__FILE__, line=__LINE__)

      call ice_distributionGet(distrb_info, numLocalBlocks=ninfo)
      if (ninfo /= distrb_info%numLocalBlocks) &
         call abort_ice(subname//' ice_distributionGet numLocalBlocks ERROR', file=__FILE__, line=__LINE__)

      allocate(blkinfo(ninfo))

      call ice_distributionGet(distrb_info, blockGlobalID = blkinfo)
      do n = 1, ninfo
         if (blkinfo(n) /= distrb_info%blockGlobalID(n)) &
            call abort_ice(subname//' ice_distributionGet blockGlobalID ERROR', file=__FILE__, line=__LINE__)
      enddo

      deallocate(blkinfo)
      allocate(blkinfo(nblocks_tot))

      call ice_distributionGet(distrb_info, blockLocation = blkinfo)
      do n = 1, nblocks_tot
         if (blkinfo(n) /= distrb_info%blockLocation(n)) &
            call abort_ice(subname//' ice_distributionGet blockLocation ERROR', file=__FILE__, line=__LINE__)
      enddo

      call ice_distributionGet(distrb_info, blockLocalID = blkinfo)
      do n = 1, nblocks_tot
         if (blkinfo(n) /= distrb_info%blockLocalID(n)) &
            call abort_ice(subname//' ice_distributionGet blockLocalID ERROR', file=__FILE__, line=__LINE__)
      enddo

      deallocate(blkinfo)

      if (my_task == master_task) then
         write(nu_diag,*) ' '
         write(nu_diag,'(2a)') subname,' ice_distributionGet checks pass'
         write(nu_diag,*) ' '
      endif
   endif

!----------------------------------------------------------------------
!
!  Set up ghost cell updates for each distribution.
!  Boundary types are cyclic, closed, tripole or tripoleT.
!
!----------------------------------------------------------------------

   ! update ghost cells on all four boundaries
   halo_info = ice_HaloCreate(distrb_info,     &
                        trim(ns_boundary_type),     &
                        trim(ew_boundary_type),     &
                        nx_global)

!----------------------------------------------------------------------

 end subroutine init_domain_distribution

!***********************************************************************

 end module ice_domain

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
