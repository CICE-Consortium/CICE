#ifdef ncdf
#define USE_NETCDF
#endif
!=======================================================================

! Spatial grids, masks, and boundary conditions
!
! authors: Elizabeth C. Hunke and William H. Lipscomb, LANL
!          Tony Craig, NCAR
!
! 2004: Block structure added by William Lipscomb
!       init_grid split into two parts as in POP 2.0
!       Boundary update routines replaced by POP versions
! 2006: Converted to free source form (F90) by Elizabeth Hunke
! 2007: Option to read from netcdf files (A. Keen, Met Office)
!       Grid reading routines reworked by E. Hunke for boundary values
! 2021: Add N (center of north face) and E (center of east face) grids
!       to support C and CD solvers.  Defining T at center of cells, U at
!       NE corner, N at center of top face, E at center of right face.
!       All cells are quadrilaterals with NE, E, and N associated with
!       directions relative to logical grid.

      module ice_grid

      use ice_kinds_mod
      use ice_broadcast, only: broadcast_scalar, broadcast_array
      use ice_boundary, only: ice_HaloUpdate, ice_HaloExtrapolate
      use ice_constants, only: c0, c1, c1p5, c2, c4, c20, c180, c360, &
          p5, p25, radius, cm_to_m, &
          field_loc_center, field_loc_NEcorner, field_loc_Nface, field_loc_Eface, &
          field_type_scalar, field_type_vector, field_type_angle
      use ice_communicate, only: my_task, master_task
      use ice_blocks, only: block, get_block, nx_block, ny_block, nghost
      use ice_domain_size, only: nx_global, ny_global, max_blocks
      use ice_domain, only: blocks_ice, nblocks, halo_info, distrb_info, &
          ew_boundary_type, ns_boundary_type, init_domain_distribution, &
          close_boundaries
      use ice_fileunits, only: nu_diag, nu_grid, nu_kmt, &
          get_fileunit, release_fileunit, flush_fileunit
      use ice_gather_scatter, only: gather_global, scatter_global
      use ice_read_write, only: ice_read, ice_read_nc, ice_read_global, &
          ice_read_global_nc, ice_open, ice_open_nc, ice_close_nc, ice_check_nc
      use ice_timers, only: timer_bound, ice_timer_start, ice_timer_stop
      use ice_exit, only: abort_ice
      use ice_global_reductions, only: global_minval, global_maxval
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_parameters, icepack_init_parameters

      implicit none
      private
      public :: init_grid1, init_grid2, grid_average_X2Y, makemask, &
                alloc_grid, dealloc_grid, &
                grid_neighbor_min, grid_neighbor_max

      character (len=char_len_long), public :: &
         grid_format  , & ! file format ('bin'=binary or 'pop_nc'= pop netcdf or 'mom_nc'=mom (supergrid) netcdf)
         gridcpl_file , & !  input file for POP coupling grid info
         grid_file    , & !  input file for POP grid info
         kmt_file     , & !  input file for POP grid info
         kmt_type     , & !  options are file, default, boxislands
         bathymetry_file, & !  input bathymetry for seabed stress
         bathymetry_format, & ! bathymetry file format (default or pop)
         grid_spacing , & !  default of 30.e3m or set by user in namelist
         grid_ice  , & !  Underlying model grid structure (A, B, C, CD)
         grid_ice_thrm, & !  ocean forcing grid for thermo fields (T, U, N, E)
         grid_ice_dynu, & !  ocean forcing grid for dyn U fields  (T, U, N, E)
         grid_ice_dynv, & !  ocean forcing grid for dyn V fields  (T, U, N, E)
         grid_atm     , & !  atmos forcing grid structure (A, B, C, CD)
         grid_atm_thrm, & !  atmos forcing grid for thermo fields (T, U, N, E)
         grid_atm_dynu, & !  atmos forcing grid for dyn U fields  (T, U, N, E)
         grid_atm_dynv, & !  atmos forcing grid for dyn V fields  (T, U, N, E)
         grid_ocn     , & !  ocean forcing grid structure (A B, C, CD)
         grid_ocn_thrm, & !  ocean forcing grid for thermo fields (T, U, N, E)
         grid_ocn_dynu, & !  ocean forcing grid for dyn U fields  (T, U, N, E)
         grid_ocn_dynv, & !  ocean forcing grid for dyn V fields  (T, U, N, E)
         grid_type        !  current options are rectangular (default),
                          !  displaced_pole, tripole, regional

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         dxT    , & ! width of T-cell through the middle (m)
         dyT    , & ! height of T-cell through the middle (m)
         dxU    , & ! width of U-cell through the middle (m)
         dyU    , & ! height of U-cell through the middle (m)
         dxN    , & ! width of N-cell through the middle (m)
         dyN    , & ! height of N-cell through the middle (m)
         dxE    , & ! width of E-cell through the middle (m)
         dyE    , & ! height of E-cell through the middle (m)
         HTE    , & ! length of eastern edge of T-cell (m)
         HTN    , & ! length of northern edge of T-cell (m)
         tarea  , & ! area of T-cell (m^2), valid in halo
         uarea  , & ! area of U-cell (m^2), valid in halo
         narea  , & ! area of N-cell (m^2), valid in halo
         earea  , & ! area of E-cell (m^2), valid in halo
         tarear , & ! 1/tarea, valid in halo
         uarear , & ! 1/uarea, valid in halo
         narear , & ! 1/narea, valid in halo
         earear , & ! 1/earea, valid in halo
         tarean , & ! area of NH T-cells
         tareas , & ! area of SH T-cells
         ULON   , & ! longitude of velocity pts, NE corner of T pts (radians)
         ULAT   , & ! latitude of velocity pts, NE corner of T pts (radians)
         TLON   , & ! longitude of temp (T) pts (radians)
         TLAT   , & ! latitude of temp (T) pts (radians)
         NLON   , & ! longitude of center of north face of T pts (radians)
         NLAT   , & ! latitude of center of north face of T pts (radians)
         ELON   , & ! longitude of center of east face of T pts (radians)
         ELAT   , & ! latitude of center of east face of T pts (radians)
         ANGLE  , & ! for conversions between POP grid and lat/lon
         ANGLET , & ! ANGLE converted to T-cells, valid in halo
         bathymetry      , & ! ocean depth, for grounding keels and bergs (m)
         ocn_gridcell_frac   ! ocean gridcell fraction
                             ! gridcell value of [1 - (land fraction)] (T-cell)

      real (kind=dbl_kind), dimension (:,:), allocatable, public :: &
         G_HTE  , & ! length of eastern edge of T-cell (global ext.)
         G_HTN      ! length of northern edge of T-cell (global ext.)

      ! grid dimensions for rectangular grid
      real (kind=dbl_kind), public ::  &
         dxrect, & !  user_specified spacing (cm) in x-direction (uniform HTN)
         dyrect    !  user_specified spacing (cm) in y-direction (uniform HTE)

      ! growth factor for variable spaced grid
      real (kind=dbl_kind), public ::  &
         dxscale, & !  scale factor for grid spacing in x direction (e.g., 1.02)
         dyscale    !  scale factor for gird spacing in y direction (e.g., 1.02)

      real (kind=dbl_kind), public :: &
         lonrefrect, & ! lower left lon for rectgrid
         latrefrect    ! lower left lat for rectgrid

      ! Corners of grid boxes for history output
      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         lont_bounds, & ! longitude of gridbox corners for T point
         latt_bounds, & ! latitude of gridbox corners for T point
         lonu_bounds, & ! longitude of gridbox corners for U point
         latu_bounds, & ! latitude of gridbox corners for U point
         lonn_bounds, & ! longitude of gridbox corners for N point
         latn_bounds, & ! latitude of gridbox corners for N point
         lone_bounds, & ! longitude of gridbox corners for E point
         late_bounds    ! latitude of gridbox corners for E point

      ! masks
      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         hm     , & ! land/boundary mask, thickness (T-cell)
         bm     , & ! task/block id
         uvm    , & ! land/boundary mask (U-cell)
         npm    , & ! land/boundary mask (N-cell)
         epm    , & ! land/boundary mask (E-cell)
         kmt        ! ocean topography mask for bathymetry (T-cell)

      logical (kind=log_kind), public :: &
         grid_outfile,   & ! flag to write out one-time grid history file
         use_bathymetry, & ! flag for reading in bathymetry_file
         save_ghte_ghtn, & ! flag for saving global hte and htn during initialization
         scale_dxdy        ! flag to apply scale factor to vary dx/dy in rectgrid

      logical (kind=log_kind), dimension (:,:,:), allocatable, public :: &
         tmask  , & ! land/boundary mask, thickness (T-cell)
         umask  , & ! land/boundary mask  (U-cell) (1 if all surrounding T cells are ocean)
         umaskCD, & ! land/boundary mask  (U-cell) (1 if at least two surrounding T cells are ocean)
         nmask  , & ! land/boundary mask, (N-cell)
         emask  , & ! land/boundary mask, (E-cell)
         opmask , & ! land/boundary orphan mask, ocean cells in atmosphere but not ocean/ice
         lmask_n, & ! northern hemisphere mask
         lmask_s    ! southern hemisphere mask

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         rndex_global       ! global index for local subdomain (dbl)

      logical (kind=log_kind), private :: &
         l_readCenter ! If anglet exist in grid file read it otherwise calculate it

      character (len=char_len), private :: &
         mask_fieldname !field/var name for the mask variable (in nc files)

      interface grid_average_X2Y
         module procedure grid_average_X2Y_base , &
                          grid_average_X2Y_userwghts, &
                          grid_average_X2Y_NEversion
      end interface

!=======================================================================

      contains

!=======================================================================
!
! Allocate space for all variables
!
      subroutine alloc_grid

      integer (int_kind) :: ierr

      character(len=*), parameter :: subname = '(alloc_grid)'

      allocate( &
         dxT      (nx_block,ny_block,max_blocks), & ! width of T-cell through the middle (m)
         dyT      (nx_block,ny_block,max_blocks), & ! height of T-cell through the middle (m)
         dxU      (nx_block,ny_block,max_blocks), & ! width of U-cell through the middle (m)
         dyU      (nx_block,ny_block,max_blocks), & ! height of U-cell through the middle (m)
         dxN      (nx_block,ny_block,max_blocks), & ! width of N-cell through the middle (m)
         dyN      (nx_block,ny_block,max_blocks), & ! height of N-cell through the middle (m)
         dxE      (nx_block,ny_block,max_blocks), & ! width of E-cell through the middle (m)
         dyE      (nx_block,ny_block,max_blocks), & ! height of E-cell through the middle (m)
         HTE      (nx_block,ny_block,max_blocks), & ! length of eastern edge of T-cell (m)
         HTN      (nx_block,ny_block,max_blocks), & ! length of northern edge of T-cell (m)
         tarea    (nx_block,ny_block,max_blocks), & ! area of T-cell (m^2)
         uarea    (nx_block,ny_block,max_blocks), & ! area of U-cell (m^2)
         narea    (nx_block,ny_block,max_blocks), & ! area of N-cell (m^2)
         earea    (nx_block,ny_block,max_blocks), & ! area of E-cell (m^2)
         tarear   (nx_block,ny_block,max_blocks), & ! 1/tarea
         uarear   (nx_block,ny_block,max_blocks), & ! 1/uarea
         narear   (nx_block,ny_block,max_blocks), & ! 1/narea
         earear   (nx_block,ny_block,max_blocks), & ! 1/earea
         tarean   (nx_block,ny_block,max_blocks), & ! area of NH T-cells
         tareas   (nx_block,ny_block,max_blocks), & ! area of SH T-cells
         ULON     (nx_block,ny_block,max_blocks), & ! longitude of U pts, NE corner (radians)
         ULAT     (nx_block,ny_block,max_blocks), & ! latitude of U pts, NE corner (radians)
         TLON     (nx_block,ny_block,max_blocks), & ! longitude of T pts (radians)
         TLAT     (nx_block,ny_block,max_blocks), & ! latitude of T pts (radians)
         NLON     (nx_block,ny_block,max_blocks), & ! longitude of N pts, N face (radians)
         NLAT     (nx_block,ny_block,max_blocks), & ! latitude of N pts, N face (radians)
         ELON     (nx_block,ny_block,max_blocks), & ! longitude of E pts, E face (radians)
         ELAT     (nx_block,ny_block,max_blocks), & ! latitude of E pts, E face (radians)
         ANGLE    (nx_block,ny_block,max_blocks), & ! for conversions between POP grid and lat/lon
         ANGLET   (nx_block,ny_block,max_blocks), & ! ANGLE converted to T-cells
         bathymetry(nx_block,ny_block,max_blocks),& ! ocean depth, for grounding keels and bergs (m)
         ocn_gridcell_frac(nx_block,ny_block,max_blocks),& ! only relevant for lat-lon grids
         hm       (nx_block,ny_block,max_blocks), & ! land/boundary mask, thickness (T-cell)
         bm       (nx_block,ny_block,max_blocks), & ! task/block id
         uvm      (nx_block,ny_block,max_blocks), & ! land/boundary mask, velocity (U-cell)
         npm      (nx_block,ny_block,max_blocks), & ! land/boundary mask (N-cell)
         epm      (nx_block,ny_block,max_blocks), & ! land/boundary mask (E-cell)
         kmt      (nx_block,ny_block,max_blocks), & ! ocean topography mask for bathymetry (T-cell)
         tmask    (nx_block,ny_block,max_blocks), & ! land/boundary mask, thickness (T-cell)
         umask    (nx_block,ny_block,max_blocks), & ! land/boundary mask, velocity (U-cell)
         umaskCD  (nx_block,ny_block,max_blocks), & ! land/boundary mask, velocity (U-cell)
         nmask    (nx_block,ny_block,max_blocks), & ! land/boundary mask (N-cell)
         emask    (nx_block,ny_block,max_blocks), & ! land/boundary mask (E-cell)
         opmask   (nx_block,ny_block,max_blocks), & ! land/boundary orphan mask (atm ocean/ice cell)
         lmask_n  (nx_block,ny_block,max_blocks), & ! northern hemisphere mask
         lmask_s  (nx_block,ny_block,max_blocks), & ! southern hemisphere mask
         rndex_global(nx_block,ny_block,max_blocks), & ! global index for local subdomain (dbl)
         lont_bounds(4,nx_block,ny_block,max_blocks), & ! longitude of gridbox corners for T point
         latt_bounds(4,nx_block,ny_block,max_blocks), & ! latitude of gridbox corners for T point
         lonu_bounds(4,nx_block,ny_block,max_blocks), & ! longitude of gridbox corners for U point
         latu_bounds(4,nx_block,ny_block,max_blocks), & ! latitude of gridbox corners for U point
         lonn_bounds(4,nx_block,ny_block,max_blocks), & ! longitude of gridbox corners for N point
         latn_bounds(4,nx_block,ny_block,max_blocks), & ! latitude of gridbox corners for N point
         lone_bounds(4,nx_block,ny_block,max_blocks), & ! longitude of gridbox corners for E point
         late_bounds(4,nx_block,ny_block,max_blocks), & ! latitude of gridbox corners for E point
         stat=ierr)
      if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory1', file=__FILE__, line=__LINE__)

      ocn_gridcell_frac(:,:,:) = -c1   ! special value to start, will be ignored unless set elsewhere

      if (save_ghte_ghtn) then
         if (my_task == master_task) then
            allocate( &
               G_HTE(nx_global+2*nghost, ny_global+2*nghost), & ! length of eastern edge of T-cell (global ext.)
               G_HTN(nx_global+2*nghost, ny_global+2*nghost), & ! length of northern edge of T-cell (global ext.)
               stat=ierr)
         else
            allocate( &
               G_HTE(1,1), & ! needed for debug checks
               G_HTN(1,1), & ! never used in code
               stat=ierr)
         endif
         if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory3', file=__FILE__, line=__LINE__)
      endif

      end subroutine alloc_grid

!=======================================================================
!
! DeAllocate space for variables no longer needed after initialization
!
      subroutine dealloc_grid

      integer (int_kind) :: ierr

      character(len=*), parameter :: subname = '(dealloc_grid)'

      if (save_ghte_ghtn) then
         deallocate(G_HTE, G_HTN, stat=ierr)
         if (ierr/=0) call abort_ice(subname//' ERROR: Dealloc error1', file=__FILE__, line=__LINE__)
      endif

      end subroutine dealloc_grid

!=======================================================================
! Distribute blocks across processors.  The distribution is optimized
! based on latitude and topography, contained in the ULAT and KMT arrays.
!
! authors: William Lipscomb and Phil Jones, LANL

      subroutine init_grid1

#ifdef USE_NETCDF
      use netcdf, only: nf90_inq_varid , nf90_noerr
      integer (kind=int_kind) :: status, varid
#endif

      integer (kind=int_kind) :: &
         fid_grid, &    ! file id for netCDF grid file
         fid_kmt        ! file id for netCDF kmt file

      character (char_len) :: &
         fieldname       ! field name in netCDF file

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g1, work_g2, work_mom

      integer (kind=int_kind) :: &
         max_blocks_min, & ! min value of max_blocks across procs
         max_blocks_max, &    ! max value of max_blocks across procs
         i, j, im, jm, ierr

      real (kind=dbl_kind) :: &
         rad_to_deg

      character(len=*), parameter :: subname = '(init_grid1)'

      !-----------------------------------------------------------------
      ! Get global ULAT and KMT arrays used for block decomposition.
      !-----------------------------------------------------------------

      call icepack_query_parameters(rad_to_deg_out=rad_to_deg)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      allocate( &
         work_g1(nx_global,ny_global), &
         work_g2(nx_global,ny_global), &
         stat=ierr &
      )
      if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory', file=__FILE__, line=__LINE__)

      ! check tripole flags here
      ! can't check in init_data because ns_boundary_type is not yet read
      ! can't check in init_domain_blocks because grid_type is not accessible due to circular logic

      if (grid_type == 'tripole' .and. ns_boundary_type /= 'tripole' .and. &
          ns_boundary_type /= 'tripoleT') then
         call abort_ice(subname//' ERROR: grid_type tripole needs tripole ns_boundary_type', &
                        file=__FILE__, line=__LINE__)
      endif

      if (grid_type == 'tripole' .and. (mod(nx_global,2)/=0)) then
         call abort_ice(subname//' ERROR: grid_type tripole requires even nx_global number', &
                        file=__FILE__, line=__LINE__)
      endif

      if (grid_format == 'mom_nc' .and. ns_boundary_type == 'tripoleT') then
         call abort_ice(subname//" ERROR: ns_boundary_type='tripoleT' not implemented "// &
                        "for grid_format='mom_nc'. Use 'tripole' instead.", &
                        file=__FILE__, line=__LINE__)
      endif

      if (trim(grid_type) == 'displaced_pole' .or. &
          trim(grid_type) == 'tripole' .or. &
          trim(grid_type) == 'regional') then

         ! Fill ULAT
         select case(trim(grid_format))
            case ('mom_nc')

               if (my_task == master_task) then
                  allocate(work_mom(nx_global*2+1, ny_global*2+1), stat=ierr)
                  if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory', file=__FILE__, line=__LINE__)

                  fieldname='y'                ! use mom y field to fill cice ULAT
                  call ice_open_nc(grid_file,fid_grid)
                  call ice_read_global_nc(fid_grid,1,fieldname,work_mom,.true.)
                  call ice_close_nc(fid_grid)
                  im = 3
                  do i = 1, nx_global
                      jm = 3
                      do j = 1, ny_global
                         work_g1(i,j) = work_mom(im, jm)
                         jm = jm + 2
                      enddo
                      im = im + 2
                  enddo

                  deallocate(work_mom, stat=ierr)
                  if (ierr/=0) call abort_ice(subname//' ERROR: Dealloc error', file=__FILE__, line=__LINE__)

               endif

            case('pop_nc', 'geosnc')

               fieldname='ulat'
               call ice_open_nc(grid_file,fid_grid)
               call ice_read_global_nc(fid_grid,1,fieldname,work_g1,.true.)
               call ice_close_nc(fid_grid)

            case default

               call ice_open(nu_grid,grid_file,64)
               call ice_read_global(nu_grid,1,work_g1,'rda8',.true.)
               if (my_task == master_task) close (nu_grid)

         end select

      else   ! rectangular grid
         work_g1(:,:) = 75._dbl_kind/rad_to_deg  ! arbitrary polar latitude
      endif

      ! Fill kmt
      if (trim(kmt_type) =='file') then
         select case(trim(grid_format))
            case ('mom_nc', 'pop_nc', 'geosnc')

               ! mask variable name might be kmt or mask, check both
               call ice_open_nc(kmt_file,fid_kmt)
#ifdef USE_NETCDF
               if ( my_task==master_task ) then
                  status = nf90_inq_varid(fid_kmt, 'kmt', varid)
                  if (status == nf90_noerr) then
                     mask_fieldname = 'kmt'
                  else
                     status = nf90_inq_varid(fid_kmt, 'mask', varid)
                     call ice_check_nc(status, subname//' ERROR: does '//trim(kmt_file)//&
                                       ' contain "kmt" or "mask" variable?', file=__FILE__, line=__LINE__)
                     mask_fieldname = 'mask'
                  endif
               endif
#endif
               call broadcast_scalar(mask_fieldname, master_task)

               call ice_read_global_nc(fid_kmt,1,mask_fieldname,work_g2,.true.)
               call ice_close_nc(fid_kmt)

            case default

               call ice_open(nu_kmt, kmt_file, 32) ! KMT
               call ice_read_global(nu_kmt, 1,work_g2,'ida4',.true.)  ! KMT
               if (my_task == master_task) close (nu_kmt)

         end select

      else
         work_g2(:,:) = c1
      endif

      call broadcast_array(work_g1, master_task)   ! ULAT
      call broadcast_array(work_g2, master_task)   ! KMT

      !-----------------------------------------------------------------
      ! distribute blocks among processors
      !-----------------------------------------------------------------

      call init_domain_distribution(work_g2, work_g1, grid_ice)  ! KMT, ULAT

      deallocate(work_g1, work_g2, stat = ierr)
      if (ierr/=0) call abort_ice(subname//' ERROR: Dealloc error', file=__FILE__, line=__LINE__)

      !-----------------------------------------------------------------
      ! write additional domain information
      !-----------------------------------------------------------------

      max_blocks_min = global_minval(max_blocks, distrb_info)
      max_blocks_max = global_maxval(max_blocks, distrb_info)
      if (my_task == master_task) then
        write(nu_diag,*        ) ''
        write(nu_diag,'(2a)'   ) subname,' Block size:'
        write(nu_diag,'(2a,i8)') subname,'   nx_block        = ',nx_block
        write(nu_diag,'(2a,i8)') subname,'   ny_block        = ',ny_block
        write(nu_diag,'(2a,i8)') subname,'   min(max_blocks) = ',max_blocks_min
        write(nu_diag,'(2a,i8)') subname,'   max(max_blocks) = ',max_blocks_max
      endif

      end subroutine init_grid1

!=======================================================================
! Horizontal grid initialization:
!
!     U{LAT,LONG} = true {latitude,longitude} of U points
!     HT{N,E} = cell widths on {N,E} sides of T cell
!     ANGLE = angle between local x direction and true east
!     hm = land mask (c1 for ocean points, c0 for land points)
!     D{X,Y}{T,U} = {x,y} spacing centered at {T,U} points
!     T-grid and ghost cell values
!     Various grid quantities needed for dynamics and transport
!
! author: Elizabeth C. Hunke, LANL

      subroutine init_grid2

#if defined (_OPENMP)
      use OMP_LIB
#endif

      integer (kind=int_kind) :: &
         i, j, iblk,      &
         ilo,ihi,jlo,jhi, &      ! beginning and end of physical domain
         ierr

      real (kind=dbl_kind) :: &
         angle_0, angle_w, angle_s, angle_sw, pi

      logical (kind=log_kind), dimension(nx_block,ny_block,max_blocks):: &
         out_of_range

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g1

      type (block) :: &
         this_block           ! block information for current block

#if defined (_OPENMP)
      integer(kind=omp_sched_kind) :: ompsk  ! openmp schedule
      integer(kind=int_kind) :: ompcs        ! openmp schedule count
#endif

      character(len=*), parameter :: subname = '(init_grid2)'

      !-----------------------------------------------------------------
      ! lat, lon, cell widths, angle, land mask
      !-----------------------------------------------------------------

      l_readCenter = .false.
      call icepack_query_parameters(pi_out=pi)

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (trim(grid_type) == 'displaced_pole' .or. &
          trim(grid_type) == 'tripole' .or. &
          trim(grid_type) == 'regional'      ) then
         select case (trim(grid_format))
            case('mom_nc')
               call mom_grid        ! derive cice grid from MOM supergrid nc file
            case ('pop_nc')
               call popgrid_nc      ! read POP grid lengths from nc file
            case ('geosnc')
               call geosgrid_nc     ! read GEOS MOM grid used from nc file
            case default
               call popgrid         ! read POP grid lengths directly
         end select
#ifdef CESMCOUPLED
      elseif (trim(grid_type) == 'latlon') then
         call latlongrid        ! lat lon grid for sequential CESM (CAM mode)
         return
#endif
      else
         call rectgrid          ! regular rectangular grid
      endif

      if (trim(kmt_type) =='none') then
         kmt(:,:,:) = c1
         hm(:,:,:)  = c1
      else if (trim(kmt_type) =='file') then
         select case (trim(grid_format))
            case('mom_nc', 'pop_nc' ,'geosnc')
               call kmtmask('nc')
            case default
               call kmtmask('bin')
         end select
      endif ! the other types are handled by rectgrid

      !-----------------------------------------------------------------
      ! Diagnose OpenMP thread schedule, force order in output
      !-----------------------------------------------------------------

#if defined (_OPENMP)
       !$OMP PARALLEL DO ORDERED PRIVATE(iblk) SCHEDULE(runtime)
       do iblk = 1, nblocks
          if (my_task == master_task) then
             !$OMP ORDERED
             if (iblk == 1) then
                call omp_get_schedule(ompsk,ompcs)
!               write(nu_diag,*) ''
                write(nu_diag,*) subname,' OpenMP runtime thread schedule:'
                write(nu_diag,*) subname,'  omp schedule = ',ompsk,ompcs
             endif
             write(nu_diag,*) subname,' block, thread = ',iblk,OMP_GET_THREAD_NUM()
             !$OMP END ORDERED
          endif
       enddo
       !$OMP END PARALLEL DO
       call flush_fileunit(nu_diag)
#endif

      !-----------------------------------------------------------------
      ! T-grid cell and U-grid cell quantities
      ! Fill halo data locally where possible to avoid missing
      ! data associated with land block elimination
      ! Note: HTN, HTE, dx*, dy* are all defined from global arrays
      ! at halos.
      !-----------------------------------------------------------------

      if (trim(grid_format) /= 'mom_nc') then
         !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
         do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            do j = 1,ny_block
            do i = 1,nx_block
               tarea(i,j,iblk) = dxT(i,j,iblk)*dyT(i,j,iblk)
               uarea(i,j,iblk) = dxU(i,j,iblk)*dyU(i,j,iblk)
               narea(i,j,iblk) = dxN(i,j,iblk)*dyN(i,j,iblk)
               earea(i,j,iblk) = dxE(i,j,iblk)*dyE(i,j,iblk)
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO
      endif

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = 1,ny_block
         do i = 1,nx_block
            if (tarea(i,j,iblk) > c0) then
               tarear(i,j,iblk) = c1/tarea(i,j,iblk)
            else
               tarear(i,j,iblk) = c0 ! possible on boundaries
            endif
            if (uarea(i,j,iblk) > c0) then
               uarear(i,j,iblk) = c1/uarea(i,j,iblk)
            else
               uarear(i,j,iblk) = c0 ! possible on boundaries
            endif
            if (narea(i,j,iblk) > c0) then
               narear(i,j,iblk) = c1/narea(i,j,iblk)
            else
               narear(i,j,iblk) = c0 ! possible on boundaries
            endif
            if (earea(i,j,iblk) > c0) then
               earear(i,j,iblk) = c1/earea(i,j,iblk)
            else
               earear(i,j,iblk) = c0 ! possible on boundaries
            endif

         enddo
         enddo

      enddo                     ! iblk
      !$OMP END PARALLEL DO

      !-----------------------------------------------------------------
      ! Ghost cell updates
      ! On the tripole grid, one must be careful with updates of
      !  quantities that involve a difference of cell lengths.
      ! For example, dyhx and dxhy are cell-centered vector components.
      ! Also note that on the tripole grid, cxp and cxm would swap places,
      !  as would cyp and cym.  These quantities are computed only
      !  in north and east ghost cells (above), not south and west.
      !-----------------------------------------------------------------

      call ice_timer_start(timer_bound)

      ! Update just on the tripole seam to ensure bit-for-bit symmetry across seam
      call ice_HaloUpdate (tarea,              halo_info, &
                           field_loc_center,   field_type_scalar, &
                           fillValue=c1,       tripoleOnly=.true.)
      call ice_HaloUpdate (uarea,              halo_info, &
                           field_loc_NEcorner, field_type_scalar, &
                           fillValue=c1,       tripoleOnly=.true.)
      call ice_HaloUpdate (narea,              halo_info, &
                           field_loc_Nface,    field_type_scalar, &
                           fillValue=c1,       tripoleOnly=.true.)
      call ice_HaloUpdate (earea,              halo_info, &
                           field_loc_Eface,    field_type_scalar, &
                           fillValue=c1,       tripoleOnly=.true.)
      call ice_HaloUpdate (tarear,             halo_info, &
                           field_loc_center,   field_type_scalar, &
                           fillValue=c1,       tripoleOnly=.true.)
      call ice_HaloUpdate (uarear,             halo_info, &
                           field_loc_NEcorner, field_type_scalar, &
                           fillValue=c1,       tripoleOnly=.true.)
      call ice_HaloUpdate (narear,             halo_info, &
                           field_loc_Nface,    field_type_scalar, &
                           fillValue=c1,       tripoleOnly=.true.)
      call ice_HaloUpdate (earear,             halo_info, &
                           field_loc_Eface,    field_type_scalar, &
                           fillValue=c1,       tripoleOnly=.true.)

      call ice_timer_stop(timer_bound)

      !-----------------------------------------------------------------
      ! Ensure that -pi <= ANGLE <= pi
      !-----------------------------------------------------------------

      out_of_range = .false.
      where (ANGLE < -pi .or. ANGLE > pi) out_of_range = .true.
      if (count(out_of_range) > 0) then
         write(nu_diag,*) subname,' angle = ',minval(ANGLE),maxval(ANGLE),count(out_of_range)
         call abort_ice (subname//' ANGLE out of expected range', &
             file=__FILE__, line=__LINE__)
      endif

      if (l_readCenter) then
         out_of_range = .false.
         where (ANGLET < -pi .or. ANGLET > pi) out_of_range = .true.
         if (count(out_of_range) > 0) then
            write(nu_diag,*) subname,' angle = ',minval(ANGLET),maxval(ANGLET),count(out_of_range)
            call abort_ice (subname//' ANGLET out of expected range', &
               file=__FILE__, line=__LINE__)
         endif
      endif

      !-----------------------------------------------------------------
      ! Compute ANGLE on T-grid
      !-----------------------------------------------------------------
      if (.not. (l_readCenter)) then
         ANGLET = c0

         !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block, &
         !$OMP                     angle_0,angle_w,angle_s,angle_sw)
         do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            do j = jlo, jhi
            do i = ilo, ihi
               angle_0  = ANGLE(i  ,j  ,iblk) !   w----0
               angle_w  = ANGLE(i-1,j  ,iblk) !   |    |
               angle_s  = ANGLE(i,  j-1,iblk) !   |    |
               angle_sw = ANGLE(i-1,j-1,iblk) !   sw---s
               ANGLET(i,j,iblk) = atan2(p25*(sin(angle_0)+ &
                                             sin(angle_w)+ &
                                             sin(angle_s)+ &
                                             sin(angle_sw)),&
                                        p25*(cos(angle_0)+ &
                                             cos(angle_w)+ &
                                             cos(angle_s)+ &
                                             cos(angle_sw)))
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO
      endif

      if ((trim(grid_type) == 'regional' .or. &
           trim(grid_type) == 'rectangular') .and. &
          (.not. (l_readCenter))) then
         ! for W boundary extrapolate from interior
         !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
         do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            i = ilo
            if (this_block%i_glob(i) == 1) then
               do j = jlo, jhi
                  ANGLET(i,j,iblk) = c2*ANGLET(i+1,j,iblk)-ANGLET(i+2,j,iblk)
               enddo
            endif
         enddo
         !$OMP END PARALLEL DO
      endif  ! regional

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (ANGLET,           halo_info, &
                           field_loc_center, field_type_angle, &
                           fillValue=c1)
      call ice_timer_stop(timer_bound)

      call makemask          ! velocity mask, hemisphere masks

      !----------------------------------------------------------------
      ! Coordinates for all T/N/E cells
      !----------------------------------------------------------------

      if (trim(grid_format) /= 'mom_nc') then
         if (.not. (l_readCenter)) then
            call Tlatlon        ! get lat, lon on the T grid
         endif
         call NElatlon          ! get lat, lon on the N, E grid

         ! corners for CF-compliant output
         call gridbox_corners
         call gridbox_edges
      endif

      !-----------------------------------------------------------------
      ! bathymetry
      !-----------------------------------------------------------------

      if (trim(bathymetry_format) == 'default') then
         call get_bathymetry
      elseif (trim(bathymetry_format) == 'pop') then
         call get_bathymetry_popfile
      else
         call abort_ice(subname//' ERROR: bathymetry_format value must be default or pop', &
            file=__FILE__, line=__LINE__)
      endif

      !-----------------------------------------------------------------
      ! Compute global index (used for unpacking messages from coupler)
      !-----------------------------------------------------------------

      if (my_task==master_task) then
         allocate(work_g1(nx_global,ny_global), stat=ierr)
         if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory', file=__FILE__, line=__LINE__)
         do j=1,ny_global
         do i=1,nx_global
            work_g1(i,j) = real((j-1)*nx_global + i,kind=dbl_kind)
         enddo
         enddo
      else
         allocate(work_g1(1,1), stat=ierr) ! to save memory
         if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory', file=__FILE__, line=__LINE__)
      endif

      call scatter_global(rndex_global, work_g1,  &
                          master_task,  distrb_info, &
                          field_loc_center, field_type_scalar)

      deallocate(work_g1, stat=ierr)
      if (ierr/=0) call abort_ice(subname//' ERROR: Dealloc error', file=__FILE__, line=__LINE__)

      end subroutine init_grid2

!=======================================================================
! POP land mask
! Land mask record number and field is (1) KMT.

      subroutine kmtmask(filetype)

      character(len=*), intent(in) :: &
         filetype        ! 'nc' or 'bin'

      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi ! beginning and end of physical domain

      integer (kind=int_kind) :: &
         fid_kmt         ! file id for netCDF kmt file

      logical (kind=log_kind) :: diag

      real (kind=dbl_kind) :: &
         puny

      type (block) :: &
         this_block       ! block information for current block

      character(len=*), parameter :: subname = '(kmtmask)'

      call icepack_query_parameters(puny_out=puny)

      diag = .true.       ! write diagnostic info

      kmt(:,:,:) = c0
      hm (:,:,:) = c0

      if (filetype == 'bin') then
         call ice_open(nu_kmt,kmt_file,32)
         call ice_read(nu_kmt,1,kmt,'ida4',diag, &
                       field_loc=field_loc_center, &
                       field_type=field_type_scalar)
         if (my_task == master_task) then
            close (nu_kmt)
         endif
      elseif (filetype == 'nc') then
         call ice_open_nc(kmt_file,fid_kmt)
         call ice_read_nc(fid_kmt,1,mask_fieldname,kmt,diag, &
                           field_loc=field_loc_center, &
                           field_type=field_type_scalar)
         call ice_close_nc(fid_kmt)
      else
         call abort_ice(subname//' ERROR: invalid filetype='//trim(filetype), file=__FILE__, line=__LINE__)
      endif

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            ! force grid cells to land if ocn_gridcell_frac is defined
            if (ocn_gridcell_frac(i,j,iblk) >= c0 .and. &
                ocn_gridcell_frac(i,j,iblk) < puny) then
               kmt(i,j,iblk)  = c0
            endif
            if (kmt(i,j,iblk) >= p5) hm(i,j,iblk) = c1
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine kmtmask

!=======================================================================
! POP displaced pole grid (or tripole).
! Grid record number, field and units are: \\
! (1) ULAT  (radians)    \\
! (2) ULON  (radians)    \\
! (3) HTN   (cm)         \\
! (4) HTE   (cm)         \\
! (5) HUS   (cm)         \\
! (6) HUW   (cm)         \\
! (7) ANGLE (radians)
!
! author: Elizabeth C. Hunke, LANL

      subroutine popgrid

      logical (kind=log_kind) :: diag

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g1

      integer (int_kind) :: ierr

      character(len=*), parameter :: subname = '(popgrid)'

      call ice_open(nu_grid,grid_file,64)

      diag = .true.       ! write diagnostic info

      !-----------------------------------------------------------------
      ! lat, lon, angle
      !-----------------------------------------------------------------

      allocate(work_g1(nx_global,ny_global), stat=ierr)
      if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory', file=__FILE__, line=__LINE__)

      call ice_read_global(nu_grid,1,work_g1,'rda8',.true.)   ! ULAT
      call gridbox_verts(work_g1,latt_bounds)
      call scatter_global(ULAT, work_g1, master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      call ice_HaloExtrapolate(ULAT, distrb_info, &
                               ew_boundary_type, ns_boundary_type)

      call ice_read_global(nu_grid,2,work_g1,'rda8',.true.)   ! ULON
      call gridbox_verts(work_g1,lont_bounds)
      call scatter_global(ULON, work_g1, master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      call ice_HaloExtrapolate(ULON, distrb_info, &
                               ew_boundary_type, ns_boundary_type)

      call ice_read_global(nu_grid,7,work_g1,'rda8',.true.)   ! ANGLE
      call scatter_global(ANGLE, work_g1, master_task, distrb_info, &
                          field_loc_NEcorner, field_type_angle)

      !-----------------------------------------------------------------
      ! cell dimensions
      ! calculate derived quantities from global arrays to preserve
      ! information on boundaries
      !-----------------------------------------------------------------

      call ice_read_global(nu_grid,3,work_g1,'rda8',.true.)   ! HTN
      call primary_grid_lengths_HTN(work_g1)                  ! dxU, dxT, dxN, dxE

      call ice_read_global(nu_grid,4,work_g1,'rda8',.true.)   ! HTE
      call primary_grid_lengths_HTE(work_g1)                  ! dyU, dyT, dyN, dyE

      deallocate(work_g1, stat=ierr)
      if (ierr/=0) call abort_ice(subname//' ERROR: Dealloc error', file=__FILE__, line=__LINE__)

      if (my_task == master_task) then
         close (nu_grid)
      endif

      end subroutine popgrid

!=======================================================================
! POP displaced pole grid and land mask.
! Grid record number, field and units are: \\
! (1) ULAT  (radians)    \\
! (2) ULON  (radians)    \\
! (3) HTN   (cm)         \\
! (4) HTE   (cm)         \\
! (5) HUS   (cm)         \\
! (6) HUW   (cm)         \\
! (7) ANGLE (radians)
!
! author: Elizabeth C. Hunke, LANL
! Revised for netcdf input: Ann Keen, Met Office, May 2007

      subroutine popgrid_nc

#ifdef USE_NETCDF
      use netcdf, only : nf90_inq_varid , nf90_inq_dimid, &
                         nf90_inquire_dimension, nf90_get_var,  nf90_noerr
#endif

      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi, &     ! beginning and end of physical domain
         fid_grid , &           ! file id for netCDF grid file
         ierr

      logical (kind=log_kind) :: diag

      character (char_len) :: &
         fieldname              ! field name in netCDF file

      real (kind=dbl_kind) :: &
         pi

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g1

      integer(kind=int_kind) :: &
         varid, status

      character(len=*), parameter :: subname = '(popgrid_nc)'

#ifdef USE_NETCDF
      call icepack_query_parameters(pi_out=pi)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      call ice_open_nc(grid_file,fid_grid)

      diag = .true.       ! write diagnostic info

      !-----------------------------------------------------------------
      ! lat, lon, angle
      !-----------------------------------------------------------------

      allocate(work_g1(nx_global,ny_global), stat=ierr)
      if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory', file=__FILE__, line=__LINE__)

      fieldname='ulat'
      call ice_read_global_nc(fid_grid,1,fieldname,work_g1,diag) ! ULAT
      call gridbox_verts(work_g1,latt_bounds)
      call scatter_global(ULAT, work_g1, master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      call ice_HaloExtrapolate(ULAT, distrb_info, &
                               ew_boundary_type, ns_boundary_type)

      fieldname='ulon'
      call ice_read_global_nc(fid_grid,1,fieldname,work_g1,diag) ! ULON
      call gridbox_verts(work_g1,lont_bounds)
      call scatter_global(ULON, work_g1, master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      call ice_HaloExtrapolate(ULON, distrb_info, &
                               ew_boundary_type, ns_boundary_type)

      fieldname='angle'
      call ice_read_global_nc(fid_grid,1,fieldname,work_g1,diag) ! ANGLE
      call scatter_global(ANGLE, work_g1, master_task, distrb_info, &
                          field_loc_NEcorner, field_type_angle)
      ! fix ANGLE: roundoff error due to single precision
      where (ANGLE >  pi) ANGLE =  pi
      where (ANGLE < -pi) ANGLE = -pi

      ! if grid file includes anglet then read instead
      fieldname='anglet'
      if (my_task == master_task) then
         status = nf90_inq_varid(fid_grid, fieldname , varid)
         if (status /= nf90_noerr) then
            write(nu_diag,*) subname//' CICE will calculate angleT, TLON and TLAT'
         else
            write(nu_diag,*) subname//' angleT, TLON and TLAT is read from grid file'
            l_readCenter = .true.
         endif
      endif
      call broadcast_scalar(l_readCenter,master_task)
      if (l_readCenter) then
         call ice_read_global_nc(fid_grid,1,fieldname,work_g1,diag)
         call scatter_global(ANGLET, work_g1, master_task, distrb_info, &
                             field_loc_center, field_type_angle)
         where (ANGLET >  pi) ANGLET =  pi
         where (ANGLET < -pi) ANGLET = -pi
         fieldname="tlon"
         call ice_read_global_nc(fid_grid,1,fieldname,work_g1,diag)
         call scatter_global(TLON, work_g1, master_task, distrb_info, &
                             field_loc_center, field_type_scalar)
         fieldname="tlat"
         call ice_read_global_nc(fid_grid,1,fieldname,work_g1,diag)
         call scatter_global(TLAT, work_g1, master_task, distrb_info, &
                             field_loc_center, field_type_scalar)
      endif
      !-----------------------------------------------------------------
      ! cell dimensions
      ! calculate derived quantities from global arrays to preserve
      ! information on boundaries
      !-----------------------------------------------------------------

      fieldname='htn'
      call ice_read_global_nc(fid_grid,1,fieldname,work_g1,diag) ! HTN
      call primary_grid_lengths_HTN(work_g1)                  ! dxU, dxT, dxN, dxE
      fieldname='hte'
      call ice_read_global_nc(fid_grid,1,fieldname,work_g1,diag) ! HTE
      call primary_grid_lengths_HTE(work_g1)                  ! dyU, dyT, dyN, dyE

      deallocate(work_g1, stat=ierr)
      if (ierr/=0) call abort_ice(subname//' ERROR: Dealloc error', file=__FILE__, line=__LINE__)

      call ice_close_nc(fid_grid)

#else
      call abort_ice(subname//' ERROR: USE_NETCDF cpp not defined', &
          file=__FILE__, line=__LINE__)
#endif

      end subroutine popgrid_nc

#ifdef CESMCOUPLED
!=======================================================================
! Read in kmt file that matches CAM lat-lon grid and has single column
! functionality
! author: Mariana Vertenstein
! 2007: Elizabeth Hunke upgraded to netcdf90 and cice ncdf calls

      subroutine latlongrid

      use ice_scam, only : scmlat, scmlon, single_column
#ifdef USE_NETCDF
      use netcdf, only : nf90_inq_varid , nf90_inq_dimid, &
                         nf90_inquire_dimension, nf90_get_var
#endif

      integer (kind=int_kind) :: &
         i, j, iblk

      integer (kind=int_kind) :: &
         ni, nj, ncid, dimid, varid, ier

      type (block) :: &
         this_block           ! block information for current block

      integer (kind=int_kind) :: &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      real (kind=dbl_kind) :: &
           closelat, &        ! Single-column latitude value
           closelon, &        ! Single-column longitude value
           closelatidx, &     ! Single-column latitude index to retrieve
           closelonidx        ! Single-column longitude index to retrieve

      integer (kind=int_kind) :: &
           start(2), &        ! Start index to read in
           count(2)           ! Number of points to read in

      integer (kind=int_kind) :: &
           start3(3), &        ! Start index to read in
           count3(3)           ! Number of points to read in

      integer (kind=int_kind) :: &
        status                ! status flag

      real (kind=dbl_kind), allocatable :: &
           lats(:),lons(:),pos_lons(:), glob_grid(:,:)  ! temporaries

      real (kind=dbl_kind) :: &
         pos_scmlon,&         ! temporary
         pi, &
         scamdata             ! temporary

      character(len=*), parameter :: subname = '(lonlatgrid)'

#ifdef USE_NETCDF
      !-----------------------------------------------------------------
      ! - kmt file is actually clm fractional land file
      ! - Determine consistency of dimensions
      ! - Read in lon/lat centers in degrees from kmt file
      ! - Read in ocean from "kmt" file (1 for ocean, 0 for land)
      !-----------------------------------------------------------------

      call icepack_query_parameters(pi_out=pi)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      ! Determine dimension of domain file and check for consistency

      if (my_task == master_task) then
         call ice_open_nc(kmt_file, ncid)

         status = nf90_inq_dimid (ncid, 'ni', dimid)
         call ice_check_nc(status, subname//' ERROR: inq_dimid ni', file=__FILE__, line=__LINE__)
         status = nf90_inquire_dimension(ncid, dimid, len=ni)
         call ice_check_nc(status, subname//' ERROR: inq dim ni', file=__FILE__, line=__LINE__)
         status = nf90_inq_dimid (ncid, 'nj', dimid)
         call ice_check_nc(status, subname//' ERROR: inq_dimid nj', file=__FILE__, line=__LINE__)
         status = nf90_inquire_dimension(ncid, dimid, len=nj)
         call ice_check_nc(status, subname//' ERROR: inq dim nj', file=__FILE__, line=__LINE__)
      end if

      ! Determine start/count to read in for either single column or global lat-lon grid
      ! If single_column, then assume that only master_task is used since there is only one task

      if (single_column) then
         ! Check for consistency
         if (my_task == master_task) then
            if ((nx_global /= 1).or. (ny_global /= 1)) then
               write(nu_diag,*) 'Because you have selected the column model flag'
               write(nu_diag,*) 'Please set nx_global=ny_global=1 in file'
               write(nu_diag,*) 'ice_domain_size.F and recompile'
               call abort_ice (subname//' ERROR: check nx_global, ny_global', file=__FILE__, line=__LINE__)
            endif
         end if

         ! Read in domain file for single column
         allocate(lats(nj))
         allocate(lons(ni))
         allocate(pos_lons(ni))
         allocate(glob_grid(ni,nj))

         start3=(/1,1,1/)
         count3=(/ni,nj,1/)
         status = nf90_inq_varid(ncid, 'xc' , varid)
         call ice_check_nc(status, subname//' ERROR: inq_varid xc', file=__FILE__, line=__LINE__)
         status = nf90_get_var(ncid, varid, glob_grid, start3, count3)
         call ice_check_nc(status, subname//' ERROR: get_var xc', file=__FILE__, line=__LINE__)
         do i = 1,ni
            lons(i) = glob_grid(i,1)
         end do

         status = nf90_inq_varid(ncid, 'yc' , varid)
         call ice_check_nc(status, subname//' ERROR: inq_varid yc', file=__FILE__, line=__LINE__)
         status = nf90_get_var(ncid, varid, glob_grid, start3, count3)
         call ice_check_nc(status, subname//' ERROR: get_var yc', file=__FILE__, line=__LINE__)
         do j = 1,nj
            lats(j) = glob_grid(1,j)
         end do

         ! convert lons array and scmlon to 0,360 and find index of value closest to 0
         ! and obtain single-column longitude/latitude indices to retrieve

         pos_lons(:)= mod(lons(:) + 360._dbl_kind,360._dbl_kind)
         pos_scmlon = mod(scmlon  + 360._dbl_kind,360._dbl_kind)
         start(1) = (MINLOC(abs(pos_lons-pos_scmlon),dim=1))
         start(2) = (MINLOC(abs(lats    -scmlat    ),dim=1))

         deallocate(lats)
         deallocate(lons)
         deallocate(pos_lons)
         deallocate(glob_grid)

         status = nf90_inq_varid(ncid, 'xc' , varid)
         call ice_check_nc(status, subname//' ERROR: inq_varid xc', file=__FILE__, line=__LINE__)
         status = nf90_get_var(ncid, varid, scamdata, start)
         call ice_check_nc(status, subname//' ERROR: get_var xc', file=__FILE__, line=__LINE__)
         TLON = scamdata
         status = nf90_inq_varid(ncid, 'yc' , varid)
         call ice_check_nc(status, subname//' ERROR: inq_varid yc', file=__FILE__, line=__LINE__)
         status = nf90_get_var(ncid, varid, scamdata, start)
         call ice_check_nc(status, subname//' ERROR: get_var yc', file=__FILE__, line=__LINE__)
         TLAT = scamdata
         status = nf90_inq_varid(ncid, 'area' , varid)
         call ice_check_nc(status, subname//' ERROR: inq_varid area', file=__FILE__, line=__LINE__)
         status = nf90_get_var(ncid, varid, scamdata, start)
         call ice_check_nc(status, subname//' ERROR: get_var are', file=__FILE__, line=__LINE__)
         tarea = scamdata
         status = nf90_inq_varid(ncid, 'mask' , varid)
         call ice_check_nc(status, subname//' ERROR: inq_varid mask', file=__FILE__, line=__LINE__)
         status = nf90_get_var(ncid, varid, scamdata, start)
         call ice_check_nc(status, subname//' ERROR: get_var mask', file=__FILE__, line=__LINE__)
         hm = scamdata
         status = nf90_inq_varid(ncid, 'frac' , varid)
         call ice_check_nc(status, subname//' ERROR: inq_varid frac', file=__FILE__, line=__LINE__)
         status = nf90_get_var(ncid, varid, scamdata, start)
         call ice_check_nc(status, subname//' ERROR: get_var frac', file=__FILE__, line=__LINE__)
         ocn_gridcell_frac = scamdata
      else
         ! Check for consistency
         if (my_task == master_task) then
            if (nx_global /= ni .and. ny_global /= nj) then
              write(nu_diag,*) 'latlongrid: ni,nj = ',ni,nj
              write(nu_diag,*) 'latlongrid: nx_g,ny_g = ',nx_global, ny_global
              call abort_ice (subname//' ERROR: ni,nj not equal to nx_global,ny_global', &
                              file=__FILE__, line=__LINE__)
            end if
         end if

         ! Read in domain file for global lat-lon grid
         call ice_read_nc(ncid, 1, 'xc'  , TLON             , diag=.true.)
         call ice_read_nc(ncid, 1, 'yc'  , TLAT             , diag=.true.)
         call ice_read_nc(ncid, 1, 'area', tarea            , diag=.true., &
            field_loc=field_loc_center,field_type=field_type_scalar)
         call ice_read_nc(ncid, 1, 'mask', hm               , diag=.true.)
         call ice_read_nc(ncid, 1, 'frac', ocn_gridcell_frac, diag=.true.)
      end if

      if (my_task == master_task) then
         call ice_close_nc(ncid)
      end if

     !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,i,j)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            ! Convert from degrees to radians
            TLON(i,j,iblk) = pi*TLON(i,j,iblk)/180._dbl_kind

            ! Convert from degrees to radians
            TLAT(i,j,iblk) = pi*TLAT(i,j,iblk)/180._dbl_kind

            ! Convert from radians^2 to m^2
            ! (area in domain file is in radians^2 and tarea is in m^2)
            tarea(i,j,iblk) = tarea(i,j,iblk) * (radius*radius)
         end do
         end do
      end do
      !$OMP END PARALLEL DO

      !-----------------------------------------------------------------
      ! Calculate various geometric 2d arrays
      ! The U grid (velocity) is not used when run with sequential CAM
      ! because we only use thermodynamic sea ice.  However, ULAT is used
      ! in the default initialization of CICE so we calculate it here as
      ! a "dummy" so that CICE will initialize with ice.  If a no ice
      ! initialization is OK (or desired) this can be commented out and
      ! ULAT will remain 0 as specified above.  ULAT is located at the
      ! NE corner of the grid cell, TLAT at the center, so here ULAT is
      ! hacked by adding half the latitudinal spacing (in radians) to
      ! TLAT.
      !-----------------------------------------------------------------

     !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,i,j)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi

            if (ny_global == 1) then
               uarea(i,j,iblk)  = tarea(i,j,  iblk)
            else
               uarea(i,j,iblk)  = p25*  &
                                 (tarea(i,j,  iblk) + tarea(i+1,j,  iblk) &
                                + tarea(i,j+1,iblk) + tarea(i+1,j+1,iblk))
            endif
            tarear(i,j,iblk)   = c1/tarea(i,j,iblk)
            uarear(i,j,iblk)   = c1/uarea(i,j,iblk)

            if (single_column) then
               ULAT  (i,j,iblk) = TLAT(i,j,iblk)+(pi/nj)
            else
               if (ny_global == 1) then
                  ULAT  (i,j,iblk) = TLAT(i,j,iblk)
               else
                  ULAT  (i,j,iblk) = TLAT(i,j,iblk)+(pi/ny_global)
               endif
            endif
            ULON  (i,j,iblk) = c0
            NLON  (i,j,iblk) = c0
            NLAT  (i,j,iblk) = c0
            ELON  (i,j,iblk) = c0
            ELAT  (i,j,iblk) = c0
            ANGLE (i,j,iblk) = c0

            ANGLET(i,j,iblk) = c0
            HTN   (i,j,iblk) = 1.e36_dbl_kind
            HTE   (i,j,iblk) = 1.e36_dbl_kind
            dxT   (i,j,iblk) = 1.e36_dbl_kind
            dyT   (i,j,iblk) = 1.e36_dbl_kind
            dxU   (i,j,iblk) = 1.e36_dbl_kind
            dyU   (i,j,iblk) = 1.e36_dbl_kind
            dxN   (i,j,iblk) = 1.e36_dbl_kind
            dyN   (i,j,iblk) = 1.e36_dbl_kind
            dxE   (i,j,iblk) = 1.e36_dbl_kind
            dyE   (i,j,iblk) = 1.e36_dbl_kind
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      call makemask
#else
      call abort_ice(subname//' ERROR: USE_NETCDF cpp not defined', &
          file=__FILE__, line=__LINE__)
#endif

      end subroutine latlongrid
#endif

!=======================================================================
! Create the CICE grid from the MOM supergrid netcdf file.
! CICE fields and units are:
! ULAT, ULON, TLAT, TLON, ELAT, ELON, NLAT, NLON (radians)
! HTN, HTE   (m)
! dxT, dyT, dxU, dyU, dxN, dyN, dxE, dyE,   (m)
! ANGLE, ANGLET (radians)
! tarea, uarea, narea, earea (m^2)
! lont_bounds, latt_bounds, etc (degrees)

      subroutine mom_grid

      integer (kind=int_kind) :: &
         fid_grid, &            ! file id for netCDF grid file
         varid, &               ! netcdf varid
         ierr

      logical (kind=log_kind) :: diag

      character (char_len) :: &
         fieldname              ! field name in netCDF file

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         G_TLON, work_gE, G_ULON, work_gN, work_mom, G_ULAT, G_TLAT, work_g1

      character(len=*), parameter :: subname = '(mom_grid)'

      call ice_open_nc(grid_file,fid_grid)

      !-----------------------------------------------------------------
      ! lat, lon, angle
      !-----------------------------------------------------------------

      if (my_task == master_task) then
         allocate( &
            work_mom(nx_global*2+1, ny_global*2+1), &
            work_gE(nx_global+1,ny_global+1)      , &
            work_gN(nx_global+1,ny_global+1)      , &
            G_ULAT(nx_global+1,ny_global+1)       , & !include left and bottom
            G_TLAT(nx_global+1,ny_global+1)       , & !include top and right
            G_TLON(nx_global+1,ny_global+1)       , & !include left and bottom
            G_ULON(nx_global+1,ny_global+1)       , & !include top and right
            stat = ierr &
         )
      else
         allocate(work_mom(1,1), work_gE(1,1), work_gN(1,1), &
            G_ULAT(1,1), G_TLAT(1,1), G_TLON(1,1), G_ULON(1,1), &
            stat=ierr)
      endif
      if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory', file=__FILE__, line=__LINE__)

      ! populate all LAT fields
      fieldname='y'
      call ice_read_global_nc(fid_grid,1,fieldname,work_mom,diag)
      call mom_corners_global(work_mom, G_ULAT, G_TLAT, work_gE, work_gN)
      ! create bounds fields for cf-compliant output
      call mom_bounds(G_ULAT, latt_bounds) ! u points define corners for t-cells
      call mom_bounds(G_TLAT, latu_bounds)
      call mom_bounds(work_gN, late_bounds)
      call mom_bounds(work_gE, latn_bounds)
      !distribute global array to local
      call mom_corners_scatter(G_ULAT, G_TLAT, work_gE, work_gN, &
                          ULAT, TLAT, ELAT, NLAT)

      ! populate all LON fields
      fieldname='x'
      call ice_read_global_nc(fid_grid,1,fieldname,work_mom,diag)
      call mom_corners_global(work_mom, G_ULON, G_TLON, work_gE, work_gN)
      call mom_bounds(G_ULON, lont_bounds)
      call mom_bounds(G_TLON, lonu_bounds)
      call mom_bounds(work_gN, lone_bounds)
      call mom_bounds(work_gE, lonn_bounds)
      call mom_corners_scatter(G_ULON, G_TLON, work_gE, work_gN, &
                                  ULON, TLON, ELON, NLON)

      deallocate(work_gE, work_gN, stat=ierr)
      if (ierr/=0) call abort_ice(subname//' ERROR: Dealloc error', file=__FILE__, line=__LINE__)
      if (my_task == master_task) then
         allocate(work_g1(nx_global, ny_global), stat=ierr)       !array for angle field
      else
         allocate(work_g1(1, 1), stat=ierr)
      endif
      if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory', file=__FILE__, line=__LINE__)

      ! populate angle fields, angle is u-points, angleT is t-points
      ! even though mom supergrid files contain angle_dx, mom6 calculates internally
      call mom_grid_rotation_angle(G_ULON, G_ULAT, G_TLON(1:nx_global,1:ny_global), work_g1) ! anglet
      call scatter_global(ANGLET, work_g1, master_task, distrb_info, &
                           field_loc_center, field_type_angle)
      call mom_grid_rotation_angle(G_TLON, G_TLAT, G_ULON(2:nx_global+1,2:ny_global+1), work_g1) ! angle
      call scatter_global(ANGLE, work_g1, master_task, distrb_info, &
                           field_loc_NEcorner, field_type_angle)

      deallocate(work_g1, G_ULAT, G_TLAT, G_TLON, G_ULON, stat=ierr)
      if (ierr/=0) call abort_ice(subname//' ERROR: Dealloc error', file=__FILE__, line=__LINE__)

      !-----------------------------------------------------------------
      ! cell dimensions
      !-----------------------------------------------------------------
      fieldname='dx'
      ! dx uses the cells in x, edges in y, reallocate work_mom to this size
      deallocate(work_mom, stat=ierr)
      if (ierr/=0) call abort_ice(subname//' ERROR: Dealloc error', file=__FILE__, line=__LINE__)
      if (my_task == master_task) then
         allocate(work_mom(nx_global*2, ny_global*2+1), stat=ierr)
      else
         allocate(work_mom(1, 1), stat=ierr)
      endif
      if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory', file=__FILE__, line=__LINE__)

      call ice_read_global_nc(fid_grid,1,fieldname,work_mom,diag)
      call mom_dx(work_mom)
      deallocate(work_mom, stat=ierr)
      if (ierr/=0) call abort_ice(subname//' ERROR: Dealloc error', file=__FILE__, line=__LINE__)

      fieldname='dy'
      ! dy uses the edges in x, cells in y, reallocate work_mom to this size
      if (my_task == master_task) then
         allocate(work_mom(nx_global*2+1, ny_global*2), stat=ierr)
      else
         allocate(work_mom(1, 1), stat=ierr)
      endif
      if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory', file=__FILE__, line=__LINE__)

      call ice_read_global_nc(fid_grid,1,fieldname,work_mom,diag)
      call mom_dy(work_mom)
      deallocate(work_mom, stat=ierr)
      if (ierr/=0) call abort_ice(subname//' ERROR: Dealloc error', file=__FILE__, line=__LINE__)


      !-----------------------------------------------------------------
      ! cell areas
      !-----------------------------------------------------------------
      fieldname = 'area'
      if (my_task == master_task) then
         allocate(work_mom(nx_global*2, ny_global*2), stat=ierr)
      else
         allocate(work_mom(1, 1), stat=ierr)
      endif
      if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory', file=__FILE__, line=__LINE__)

      call ice_read_global_nc(fid_grid,1,fieldname,work_mom,diag)
      call mom_area(work_mom)
      deallocate(work_mom, stat=ierr)
      if (ierr/=0) call abort_ice(subname//' ERROR: Dealloc', file=__FILE__, line=__LINE__)

      !-----------------------------------------------------------------
      ! fin
      !-----------------------------------------------------------------
      call ice_close_nc(fid_grid)
      l_readCenter = .true.       ! we have read t quantities

      end subroutine mom_grid

!=======================================================================

      subroutine mom_corners_global(work_mom, G_U, G_T, G_E, G_N)

      ! mom supergrid has four cells for every model cell
      ! we need to select the correct edges to get lat & lon for a model cell
      ! we include left/bottom edges for U-points, and top/right edges for T-points
      ! and close per ew_boundary_type & ns_boundary_type

      real (kind=dbl_kind), dimension(:,:), intent(in) :: work_mom
         ! supergrid array of x or y

      real (kind=dbl_kind), dimension(:,:), intent(out) :: G_U, G_T, G_E, G_N
         ! global grids

      integer (kind=int_kind) :: &
         i, j, &
         im1, im2, jm1, jm2  ! i & j for mom supergrid

      character(len=*), parameter :: subname = '(mom_corners_global)'

      if (my_task == master_task) then

         im1 = 1 ; im2 = 2  ! lh , middle  hand edge of first col
         do i = 1, nx_global
            jm1 = 1; jm2 = 2 ! bottom, middle of first row
            do j = 1, ny_global
               G_U(i,j) = work_mom(im1, jm1)     ! ULAT/LON
               G_N(i,j) = work_mom(im2, jm1)     ! NLAT/LON
               G_E(i,j) = work_mom(im1, jm2)     ! ELAT/LON
               G_T(i,j) = work_mom(im2, jm2)     ! TLAT/LON
               jm1 = jm1 + 2 ; jm2 = jm2 + 2
            enddo
            im1 = im1 + 2 ; im2 = im2 + 2
         enddo

         ! fill last col
         jm1 = 1; jm2 = 2 ! bottom, middle of first row
         do j = 1, ny_global
            G_U(nx_global+1,j) = work_mom(2*nx_global+1, jm1)
            G_E(nx_global+1,j) = work_mom(2*nx_global+1, jm2)
            jm1 = jm1 + 2 ; jm2 = jm2 + 2
         enddo
         select case (trim(ew_boundary_type))
            case('cyclic')
               G_T(nx_global+1,:) = G_T(1,:)
               G_N(nx_global+1,:) = G_N(1,:)
            case('open')
               do j=1, ny_global+1
                  G_T(nx_global+1,j) = 2 * G_T(nx_global, j) - G_T(nx_global-1, j)
                  G_N(nx_global+1,j) = 2 * G_N(nx_global, j) - G_N(nx_global-1, j)
               enddo
         end select

         ! fill last row
         im1 = 1 ; im2 = 2
         do i = 1, nx_global+1
            G_U(i,ny_global + 1) = work_mom(im1, 2*ny_global+1)
            G_N(i,ny_global + 1) = work_mom(im2, 2*ny_global+1)
            im1 = im1 + 2
         enddo
         select case (trim(ns_boundary_type))
            case ('tripole')
               do i = 1, nx_global+1
                  G_T(i,ny_global+1) = G_T(nx_global+1-i, ny_global)
                  G_E(i,ny_global+1) = G_E(nx_global+1-i, ny_global)
               enddo
            case ('cyclic')
               G_T(:,ny_global+1) = G_T(:,1)
               G_E(:,ny_global+1) = G_E(:,1)
            case ('open')
               do i = 1, nx_global+1
                  G_T(i,ny_global+1) = 2 * G_T(i, ny_global) - G_T(i, ny_global-1)
                  G_E(i,ny_global+1) = 2 * G_E(i, ny_global) - G_E(i, ny_global-1)
               enddo
         end select

      endif

      end subroutine mom_corners_global

!=======================================================================

      subroutine mom_bounds(G_corners, bounds)

      ! with an global array of corner points, subset and distribute
      ! into a cice bounds variables
      ! e.g. The tracer coordinates are the corner of the u-cells,
      ! so use mom_bounds(G_TLON, lonu_bounds)

      real (kind=dbl_kind), dimension(:,:), intent(in) :: G_corners
      real (kind=dbl_kind), dimension(:,:,:,:), intent(out) :: bounds

      ! local vars

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work_bounds

      character(len=*), parameter :: subname = '(mom_bounds)'

      ! Get bounds of grid boxes for each block as follows:
      ! (1) SW corner, (2) SE corner, (3) NE corner, (4) NW corner
      call scatter_global(work_bounds, G_corners(1:nx_global, 1:ny_global), &
                         master_task, distrb_info, &
                         field_loc_NEcorner, field_type_scalar)
      bounds(1,:,:,:) = work_bounds(:,:,:)
      call scatter_global(work_bounds, G_corners(2:nx_global+1, 1:ny_global), &
                         master_task, distrb_info, &
                         field_loc_NEcorner, field_type_scalar)
      bounds(2,:,:,:) = work_bounds(:,:,:)
      call scatter_global(work_bounds, G_corners(2:nx_global+1, 2:ny_global+1), &
                         master_task, distrb_info, &
                         field_loc_NEcorner, field_type_scalar)
      bounds(3,:,:,:) = work_bounds(:,:,:)
      call scatter_global(work_bounds, G_corners(1:nx_global, 2:ny_global+1), &
                         master_task, distrb_info, &
                         field_loc_NEcorner, field_type_scalar)
      bounds(4,:,:,:) = work_bounds(:,:,:)

      end subroutine mom_bounds

!=======================================================================

      subroutine mom_corners_scatter(G_U, G_T, G_E, G_N, U, T, E, N )

      ! with a global array of corner points in degrees, convert to rad and scatter to workers

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: G_U, G_T, G_E, G_N
         ! global grids

      real (kind=dbl_kind), dimension(:,:,:), intent(out) :: U, T, E, N ! local grids

      real (kind=dbl_kind) :: deg_to_rad , pi

      character(len=*), parameter :: subname = '(mom_corners_scatter)'

      call icepack_query_parameters(pi_out=pi)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      deg_to_rad = pi/c180

      ! convert to rad
      G_T = G_T * deg_to_rad
      G_U = G_U * deg_to_rad
      G_N = G_N * deg_to_rad
      G_E = G_E * deg_to_rad

      ! distribute to processors
      ! subset G_T to active cells by dropping right/top halo
      call scatter_global(T, G_T(1:nx_global, 1:ny_global), &
                          master_task, distrb_info, &
                          field_loc_center, field_type_scalar)
      call ice_HaloExtrapolate(T, distrb_info, &
                               ew_boundary_type, ns_boundary_type)
      ! subset G_U/G_E/G_N to active cells by dropping left/bottom edge
      call scatter_global(U, G_U(2:nx_global+1, 2:ny_global+1), &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      call ice_HaloExtrapolate(U, distrb_info, &
                               ew_boundary_type, ns_boundary_type)
      call scatter_global(N, G_N(1:nx_global, 2:ny_global+1), master_task, distrb_info, &
                         field_loc_Nface, field_type_scalar)
      call ice_HaloExtrapolate(N, distrb_info, &
                               ew_boundary_type, ns_boundary_type)
      call scatter_global(E, G_E(2:nx_global+1, 1:ny_global), master_task, distrb_info, &
                          field_loc_Eface, field_type_scalar)
      call ice_HaloExtrapolate(E, distrb_info, &
                               ew_boundary_type, ns_boundary_type)

      end subroutine mom_corners_scatter

!=======================================================================

      subroutine mom_dx(work_mom)

      ! mom supergrid has four cells for every model cell, sum the sidelengths to get model dx

      real (kind=dbl_kind), dimension(:,:) :: work_mom

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         G_dxT, G_dxN, G_dxE, G_dxU

      integer (kind=int_kind) :: &
         i, j , &
         im1, im2, jm1, jm2, im3, jm3 , &  ! i & j for mom supergrid
         ierr

      character(len=*), parameter :: subname = '(mom_dx)'

      if (my_task == master_task) then
         allocate( &
            G_dxT(nx_global,ny_global), &
            G_dxN(nx_global,ny_global), &
            G_dxE(nx_global,ny_global), &
            G_dxU(nx_global,ny_global), &
            stat=ierr &
         )
      else
         allocate(G_dxT(1,1), G_dxE(1,1), G_dxU(1,1), G_dxN(1,1), stat=ierr)
      endif
      if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory', file=__FILE__, line=__LINE__)

      if (my_task == master_task) then
         im1 = 1 ; im2 = 2 ! left ; center column of first t-cell
         im3 = 3 ! left column of second T-cell, (ie right column of U-cell)
         do i = 1, nx_global - 1
            jm1 = 2 ; jm2 = 3 ! middle , top of first row
            do j = 1, ny_global
               G_dxT(i,j) = work_mom(im1, jm1) + work_mom(im2, jm1)     !dxT
               G_dxN(i,j) = work_mom(im1, jm2) + work_mom(im2, jm2)     !dxN
               G_dxE(i,j) = work_mom(im2, jm1) + work_mom(im3, jm1)     !dxE
               G_dxU(i,j) = work_mom(im2, jm2) + work_mom(im3, jm2)     !dxU
               jm1 = jm1 + 2 ; jm2 = jm2 + 2
            enddo
            im1 = im1 + 2 ; im2 = im2 + 2 ; im3 = im3 + 2
         enddo

         ! fill the last col
         jm1 = 2 ; jm2 = 3 ! middle , top of first row
         do j = 1, ny_global
            G_dxT(nx_global,j) = work_mom(2*nx_global - 1, jm1) + work_mom(2*nx_global, jm1)     !dxT
            G_dxN(nx_global,j) = work_mom(2*nx_global - 1, jm2) + work_mom(2*nx_global, jm2)     !dxN
            jm1 = jm1 + 2 ; jm2 = jm2 + 2
         enddo
         jm1 = 2 ; jm2 = 3 ! middle , top of first row
         if (trim(ew_boundary_type) == 'cyclic') then
            do j = 1, ny_global
               G_dxE(nx_global,j) = work_mom(2*nx_global, jm1) + work_mom(1, jm1)     !dxE
               G_dxU(nx_global,j) = work_mom(2*nx_global, jm2) + work_mom(1, jm2)     !dxU
               jm1 = jm1 + 2 ; jm2 = jm2 + 2
            enddo
         else if (trim(ew_boundary_type) == 'open') then
            do j = 1, ny_global
               G_dxE(nx_global,j) = 4*work_mom(2*nx_global, jm1) - 2*work_mom(2*nx_global-1, jm1)     !dxE
               G_dxU(nx_global,j) = 4*work_mom(2*nx_global, jm2) - 2*work_mom(2*nx_global-1, jm2)     !dxU
               jm1 = jm1 + 2 ; jm2 = jm2 + 2
            enddo
         endif

         if (save_ghte_ghtn) then
            do j = 1, ny_global
               do i = 1, nx_global
                  G_HTN(i+nghost,j+nghost) = G_dxN(i,j)
               enddo
            enddo
            call global_ext_halo(G_HTN)
         endif
      endif

      call scatter_global(dxT, G_dxT, master_task, distrb_info, &
                           field_loc_center, field_type_scalar)
      call scatter_global(HTN, G_dxN, master_task, distrb_info, &
                           field_loc_Nface, field_type_scalar)
      dxN(:,:,:) = HTN(:,:,:)
      call scatter_global(dxE, G_dxE, master_task, distrb_info, &
                           field_loc_center, field_type_scalar)
      call scatter_global(dxU, G_dxU, master_task, distrb_info, &
                           field_loc_NEcorner, field_type_scalar)

      deallocate(G_dxT, G_dxE, G_dxU, G_dxN, stat=ierr)
      if (ierr/=0) call abort_ice(subname//' ERROR: Dealloc error', file=__FILE__, line=__LINE__)

      end subroutine mom_dx

!=======================================================================

      subroutine mom_dy(work_mom)

      ! mom supergrid has four cells for every model cell, sum the sidelengths to get model dy

      real (kind=dbl_kind), dimension(:,:) :: work_mom

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         G_dyT, G_dyN, G_dyE, G_dyU

      integer (kind=int_kind) :: &
         i, j, &
         im1, im2, jm1, jm2, im3, jm3 , &  ! i & j for mom supergrid
         ierr

      character(len=*), parameter :: subname = '(mom_dy)'

      if (my_task == master_task) then
         allocate( &
            G_dyT(nx_global,ny_global), &
            G_dyN(nx_global,ny_global), &
            G_dyE(nx_global,ny_global), &
            G_dyU(nx_global,ny_global), &
            stat=ierr &
         )
      else
         allocate(G_dyT(1,1), G_dyE(1,1), G_dyU(1,1), G_dyN(1,1), stat=ierr)
      endif
      if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory', file=__FILE__, line=__LINE__)

      if (my_task == master_task) then
         im1 = 2 ; im2 = 3 ! middle , right edge of first T-cell
         do i = 1, nx_global
            jm1 = 1 ; jm2 = 2 ; jm3 = 3
            do j = 1, ny_global - 1
               G_dyT(i,j) = work_mom(im1, jm1) + work_mom(im1, jm2)     !dyT
               G_dyN(i,j) = work_mom(im1, jm2) + work_mom(im1, jm3)     !dyN
               G_dyE(i,j) = work_mom(im2, jm1) + work_mom(im2, jm2)     !dyE
               G_dyU(i,j) = work_mom(im2, jm2) + work_mom(im2, jm3)     !dyU
               jm1 = jm1 + 2 ; jm2 = jm2 + 2 ; jm3 = jm3 + 2
            enddo
            im1 = im1 + 2 ; im2 = im2 + 2
         enddo

         ! fill the top row
         im1 = 2 ; im2 = 3 ! middle , right edge of first column
         do i = 1, nx_global
            G_dyT(i,ny_global) = work_mom(im1, 2*ny_global - 1) + work_mom(im1, 2*ny_global)                   !dyT
            G_dyE(i,ny_global) = work_mom(im2, 2*ny_global - 1) + work_mom(im2, 2*ny_global)                   !dyE
            im1 = im1 + 2 ; im2 = im2 + 2
         enddo
         im1 = 2 ; im2 = 3
         if (trim(ns_boundary_type)  == 'tripole') then
            do i = 1, nx_global
               G_dyN(i,ny_global) = work_mom(im1, 2*ny_global) + work_mom(2*nx_global+2-im1, 2*ny_global)      !dyN
               G_dyU(i,ny_global) = work_mom(im2, 2*ny_global) + work_mom(2*nx_global+2-im2, 2*ny_global)      !dyU
               im1 = im1 + 2 ; im2 = im2 + 2
            enddo
         else if (trim(ns_boundary_type) == 'cyclic') then
            do i = 1, nx_global
               G_dyN(i,ny_global) = work_mom(im1, 2*ny_global) + work_mom(im1, 1)                              !dyN
               G_dyU(i,ny_global) = work_mom(im2, 2*ny_global) + work_mom(im2, 1)                              !dyU
               im1 = im1 + 2 ; im2 = im2 + 2
            enddo
         else if (trim(ns_boundary_type) == 'open') then
            do i = 1, nx_global
               G_dyN(i,ny_global) = 4*work_mom(im1, 2*ny_global) - 2*work_mom(im1, 2*ny_global-1)               !dyN
               G_dyU(i,ny_global) = 4*work_mom(im2, 2*ny_global) - 2*work_mom(im2, 2*ny_global-1)               !dyU
               im1 = im1 + 2 ; im2 = im2 + 2
            enddo
         endif

         if (save_ghte_ghtn) then
            do j = 1, ny_global
               do i = 1, nx_global
                  G_HTE(i+nghost,j+nghost) = G_dyE(i,j)
               enddo
            enddo
            call global_ext_halo(G_HTE)
         endif
      endif

      call scatter_global(dyT, G_dyT, master_task, distrb_info, &
      field_loc_center, field_type_scalar)
      call scatter_global(dyN, G_dyN, master_task, distrb_info, &
            field_loc_Nface, field_type_scalar)
      call scatter_global(HTE, G_dyE, master_task, distrb_info, &
            field_loc_center, field_type_scalar)
      dyE(:,:,:) = HTE(:,:,:)
      call scatter_global(dyU, G_dyU, master_task, distrb_info, &
            field_loc_NEcorner, field_type_scalar)

      deallocate(G_dyT, G_dyN, G_dyE, G_dyU)
      if (ierr/=0) call abort_ice(subname//' ERROR: Dealloc error', file=__FILE__, line=__LINE__)

      end subroutine mom_dy

!=======================================================================

      subroutine mom_area(work_mom)

      ! mom supergrid has four cells for every model cell, sum these
      ! to get uarea and tarea
      ! earea and narea are calculated from dx & dy - see https://github.com/NOAA-GFDL/MOM6/issues/740

      real (kind=dbl_kind), dimension(:,:), intent(in) :: work_mom

      integer (kind=int_kind) :: &
         i, j, iblk, &
         im1, im2, jm1, jm2, im3, jm3 , & ! i & j for mom supergrid
         ilo,ihi,jlo,jhi , &    ! beginning and end of physical domain
         ierr

      type (block) :: &
         this_block           ! block information for current block

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         G_tarea, G_uarea

      character(len=*), parameter :: subname = '(mom_area)'

      ! calculate narea and earea
      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = 1,ny_block
         do i = 1,nx_block
            narea(i,j,iblk) = dxN(i,j,iblk)*dyN(i,j,iblk)
            earea(i,j,iblk) = dxE(i,j,iblk)*dyE(i,j,iblk)
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      if (my_task == master_task) then
         allocate( &
            G_tarea(nx_global,ny_global), &
            G_uarea(nx_global,ny_global), &
            stat=ierr &
         )
      else
         allocate(G_tarea(1,1), G_uarea(1,1), stat=ierr )
      endif
      if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory', file=__FILE__, line=__LINE__)

      ! load tarea and uarea
      if (my_task == master_task) then
         im1 = 1 ; im2 = 2 ! left/right -half of first column
         im3 = 3 ! right of first U - cell
         do i = 1, nx_global - 1
            jm1 = 1 ; jm2 = 2 ! bottom/top -half of first row
            jm3 = 3 ! top of first U - cell
            do j = 1, ny_global - 1
               G_tarea(i,j) = work_mom(im1, jm1) + work_mom(im1, jm2) &
                              + work_mom(im2, jm1) + work_mom(im2, jm2)
               G_uarea(i,j) = work_mom(im2, jm2) + work_mom(im2, jm3) &
                              + work_mom(im3, jm2) + work_mom(im3, jm3)
               jm1 = jm1 + 2 ; jm2 = jm2 + 2 ; jm3 = jm3 + 2
            enddo
            im1 = im1 + 2 ; im2 = im2 + 2 ; im3 = im3 + 2
         enddo

         ! fill last column
         jm1 = 1 ; jm2 = 2 ; jm3 = 3
         im1 = 2*nx_global - 1 ; im2 = 2*nx_global ; im3 = 1
         do j = 1, ny_global - 1
            G_tarea(nx_global,j) = work_mom(im1, jm1) + work_mom(im1, jm2) &
                                 + work_mom(im2, jm1) + work_mom(im2, jm2)
            if (trim(ew_boundary_type) == 'cyclic') then
               G_uarea(nx_global,j) = work_mom(im2, jm2) + work_mom(im2, jm3) &
                                    + work_mom(im3, jm2) + work_mom(im3, jm3)
            else if (trim(ew_boundary_type) == 'open') then
               G_uarea(nx_global,j) = 4*work_mom(im2, jm2) + 4*work_mom(im2, jm3) &
                                    - 2*work_mom(im1, jm2) - 2*work_mom(im1, jm3)
            endif
            jm1 = jm1 + 2 ; jm2 = jm2 + 2 ; jm3 = jm3 + 2
         enddo

         ! fill last row
         jm1 = ny_global*2 - 1 ; jm2 = ny_global*2
         im1 = 1 ; im2 = 2 ; im3 = 3
         do i = 1, nx_global -1
            G_tarea(i,ny_global) = work_mom(im1, jm1) + work_mom(im1, jm2) &
                                 + work_mom(im2, jm1) + work_mom(im2, jm2)
            if (trim(ns_boundary_type) == 'tripole') then
               G_uarea(i,ny_global) = work_mom(im2, jm2) + work_mom(2*nx_global+1-im2, jm2) &
                                    + work_mom(im3, jm2) + work_mom(2*nx_global+1-im3, jm2)
            else if (trim(ns_boundary_type) == 'cyclic') then
               G_uarea(i,ny_global) = work_mom(im2, jm2) + work_mom(im2, jm3) &
                                    + work_mom(im3, jm2) + work_mom(im3, jm3)
            else if (trim(ns_boundary_type) == 'open') then
               G_uarea(i,ny_global) = 4*work_mom(im2, jm2) + 4*work_mom(im3, jm2) &
                                    - 2*work_mom(im2, jm1) - 2*work_mom(im3, jm1)
            endif
            im1 = im1 + 2 ; im2 = im2 + 2 ; im3 = im3 + 2
         enddo

         ! the top right corner
         im1 = nx_global*2-1 ; im2 = nx_global*2
         jm1 = ny_global*2-1 ; jm2 = ny_global*2
         G_tarea(nx_global,ny_global) = work_mom(im1, jm1) + work_mom(im1, jm2) &
                                       + work_mom(im2, jm1) + work_mom(im2, jm2)
         if (trim(ns_boundary_type) == 'tripole') then
            G_uarea(nx_global,ny_global) = 2*(work_mom(im2, jm2) + work_mom(1, jm2))
         else if (trim(ns_boundary_type) == 'cyclic' &
                  .and. trim(ew_boundary_type) == 'cyclic') then
            G_uarea(nx_global,ny_global) = work_mom(im2, jm2) + work_mom(1, jm2) &
                                          + work_mom(im2, 1) + work_mom(1, 1)
         else if (trim(ns_boundary_type) == 'cyclic' &
                  .and. trim(ew_boundary_type) == 'open') then
            G_uarea(nx_global,ny_global) = 4*work_mom(im2, jm2) + 4*work_mom(im2, 1) &
                                          - 2*work_mom(im1, jm2) - 2*work_mom(im1, 1)
         else if (trim(ns_boundary_type) == 'open' &
                  .and. trim(ew_boundary_type) == 'cyclic') then
            G_uarea(nx_global,ny_global) = 4*work_mom(im2, jm2) + 4*work_mom(1, jm2) &
                                          - 2*work_mom(im2, jm1) - 2*work_mom(1, jm1)
         else if (trim(ns_boundary_type) == 'open' &
                  .and. trim(ew_boundary_type) == 'open') then
            G_uarea(nx_global,ny_global) = 8*work_mom(im2, jm2) &
                                 - 2*work_mom(im2, jm1) - 2*work_mom(im1, jm2)
         endif

      endif

      call scatter_global(tarea, G_tarea, master_task, distrb_info, &
                         field_loc_center, field_type_scalar)
      call scatter_global(uarea, G_uarea, master_task, distrb_info, &
                         field_loc_NEcorner, field_type_scalar)
      deallocate(G_tarea, G_uarea, stat=ierr)
      if (ierr/=0) call abort_ice(subname//' ERROR: Dealloc error', file=__FILE__, line=__LINE__)

      end subroutine mom_area

!=======================================================================

      subroutine mom_grid_rotation_angle(lon_cnr, lat_cnr, lon_cen, angle)
      !  create angles in the same way mom6 creates the angle
      !  based on https://github.com/mom-ocean/MOM6/blob/129e1bda02d454fb280819d1d87ae16347fd044c/src/initialization/MOM_shared_initialization.F90#L535
      !  the angle is between logical north on the grid and true north.

      ! global lat/lons/angles
      real (kind=dbl_kind), dimension(:,:), intent(in) :: &
         lon_cnr, & ! array of lon corner points
         lat_cnr, & ! array of lat corner points
         lon_cen ! array of lon centre points (i.e. the location the angle is calculated for)
      real (kind=dbl_kind), dimension(:,:), intent(out) :: angle

      ! local vars
      real (kind=dbl_kind)   :: &
         lon_scale, &  ! The trigonometric scaling factor converting changes in longitude to equivalent distances in latitudes [nondim]
         len_lon, &
         lon_adj, &
         lonB(2,2)
      integer (kind=int_kind) :: i, j, m, n

      character(len=*), parameter :: subname = '(mom_grid_rotation_angle)'

      if (my_task == master_task) then
         len_lon = maxval(lon_cnr)-minval(lon_cnr)  ! The periodic range of longitudes, usually 2pi.

         do j=1,ny_global
            do i=1,nx_global
               lon_adj = lon_cen(i,j)-p5*len_lon
               do n=1,2 ; do m=1,2
                  ! shift 4 lon corner points to be similar range to centre point
                  ! e.g. upper limit of 0 might be shifted to 2*pi
                  lonB(m,n) = modulo(lon_cnr(i+m-1,j+n-1)-lon_adj, len_lon) &
                                    + lon_adj
               enddo ; enddo
               lon_scale = cos(p25*(lat_cnr(I,J) + lat_cnr(I+1,J+1) + lat_cnr(I+1,J) + lat_cnr(I,J+1)))
               angle(i,j) = atan2(lon_scale*((lonB(1,2) - lonB(2,1) + lonB(2,2) - lonB(1,1))), &
                           (lat_cnr(I,J+1) - lat_cnr(I+1,J) + lat_cnr(I+1,J+1) - lat_cnr(I,J)) )
            enddo
         enddo
      endif

      end subroutine mom_grid_rotation_angle

!=======================================================================
! GEOS MOM grid
! Grid record number, field and units are: \\
! (1) ULAT  (radians)    \\
! (2) ULON  (radians)    \\
! (3) ANGLE (radians)    \\
! (4) ANGLET (radians)   \\
! (5) HTN   (cm)         \\
! (6) HTE   (cm)         \\
!
! Land mask record number and field is (1) KMT.
!

      subroutine geosgrid_nc

      use ice_blocks, only: nx_block, ny_block
      use ice_constants, only: c0, c1, &
          field_loc_center, field_loc_NEcorner, &
          field_type_scalar, field_type_angle
      use ice_domain_size, only: max_blocks
#ifdef USE_NETCDF
      use netcdf
#endif

      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi, &    ! beginning and end of physical domain
         fid_grid              ! file id for netCDF grid file

      logical (kind=log_kind) :: diag

      character (char_len) :: &
         fieldname             ! field name in netCDF file

      real (kind=dbl_kind) :: &
         pi

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g1

      type (block) :: &
         this_block            ! block information for current block

      integer(kind=int_kind) :: &
         varid
      integer (kind=int_kind) :: &
         status                ! status flag

      character(len=*), parameter :: subname = '(geosgrid_nc)'

#ifdef USE_NETCDF
      call icepack_query_parameters(pi_out=pi)
      call icepack_warnings_flush(nu_diag)
      if   (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      call ice_open_nc(grid_file,fid_grid)

      diag = .true.       ! write diagnostic info
      l_readCenter = .false.

      !-----------------------------------------------------------------
      ! lat, lon, angle
      !-----------------------------------------------------------------

      allocate(work_g1(nx_global,ny_global))

      fieldname='ulat'
      call ice_read_global_nc(fid_grid,1,fieldname,work_g1,diag) ! ULAT
      call gridbox_verts(work_g1,latt_bounds)
      call scatter_global(ULAT, work_g1, master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      call ice_HaloExtrapolate(ULAT, distrb_info, &
                               ew_boundary_type, ns_boundary_type)

      fieldname='ulon'
      call ice_read_global_nc(fid_grid,1,fieldname,work_g1,diag) ! ULON
      call gridbox_verts(work_g1,lont_bounds)
      call scatter_global(ULON, work_g1, master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      call ice_HaloExtrapolate(ULON, distrb_info, &
                               ew_boundary_type, ns_boundary_type)

      fieldname='angle'
      call ice_read_global_nc(fid_grid,1,fieldname,work_g1,diag) ! ANGLE
      call scatter_global(ANGLE, work_g1, master_task, distrb_info, &
                          field_loc_NEcorner, field_type_angle)
      ! fix ANGLE: roundoff error due to single precision
      where (ANGLE >  pi) ANGLE =  pi
      where (ANGLE < -pi) ANGLE = -pi

      ! if grid file includes anglet then read instead
      fieldname='anglet'
      if (my_task == master_task) then
         status = nf90_inq_varid(fid_grid, trim(fieldname) , varid)
         if (status /= nf90_noerr) then
            write(nu_diag,*) subname//' CICE will calculate angleT, TLON and TLAT'
         else
            write(nu_diag,*) subname//' angleT, TLON and TLAT is read from grid file'
            l_readCenter = .true.
         endif
      endif
      call broadcast_scalar(l_readCenter,master_task)
      if (l_readCenter) then
         call ice_read_global_nc(fid_grid,1,fieldname,work_g1,diag)
         call scatter_global(ANGLET, work_g1, master_task, distrb_info, &
                             field_loc_center, field_type_angle)
         where (ANGLET >  pi) ANGLET =  pi
         where (ANGLET < -pi) ANGLET = -pi
         fieldname="tlon"
         call ice_read_global_nc(fid_grid,1,fieldname,work_g1,diag)
         call scatter_global(TLON, work_g1, master_task, distrb_info, &
                             field_loc_center, field_type_scalar)
         fieldname="tlat"
         call ice_read_global_nc(fid_grid,1,fieldname,work_g1,diag)
         call scatter_global(TLAT, work_g1, master_task, distrb_info, &
                             field_loc_center, field_type_scalar)
      endif
      !-----------------------------------------------------------------
      ! cell dimensions
      ! calculate derived quantities from global arrays to preserve
      ! information on boundaries
      !-----------------------------------------------------------------

      fieldname='htn'
      call ice_read_global_nc(fid_grid,1,fieldname,work_g1,diag) ! HTN
      call primary_grid_lengths_HTN(work_g1)                  ! dxu, dxt
      fieldname='hte'
      call ice_read_global_nc(fid_grid,1,fieldname,work_g1,diag) ! HTE
      call primary_grid_lengths_HTE(work_g1)                  ! dyu, dyt

      deallocate(work_g1)

      if (my_task == master_task) then
         call ice_close_nc(fid_grid)
      endif
#else
      call abort_ice(subname//'ERROR: USE_NETCDF cpp not defined', &
          file=__FILE__, line=__LINE__)
#endif

      end subroutine geosgrid_nc

!=======================================================================
! Regular rectangular grid and mask
!
! author: Elizabeth C. Hunke, LANL

      subroutine rectgrid

      integer (kind=int_kind) :: &
         i, j, &
         imid, jmid

      real (kind=dbl_kind) :: &
         length,  &
         rad_to_deg

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g1

      character(len=*), parameter :: subname = '(rectgrid)'

      !-----------------------------------------------------------------
      ! Calculate various geometric 2d arrays
      !-----------------------------------------------------------------

      call icepack_query_parameters(rad_to_deg_out=rad_to_deg)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      hm (:,:,:) = c0
      kmt(:,:,:) = c0
      angle(:,:,:) = c0   ! "square with the world"

      allocate(work_g1(nx_global,ny_global))

      if (scale_dxdy) then
         ! scale grid spacing from center outward.
         ! this different than original method in it
         ! needs to define grid spacing before lat/lon.
         ! original rectgrid defines latlon first
         call rectgrid_scale_dxdy
      else
         ! rectgrid no grid spacing.
         ! original method with addition to use namelist lat/lon reference

         if (my_task == master_task) then
            work_g1 = c0
            length = dxrect*cm_to_m/radius*rad_to_deg

            work_g1(1,:) = lonrefrect ! reference lon from namelist

            do j = 1, ny_global
            do i = 2, nx_global
               work_g1(i,j) = work_g1(i-1,j) + length   ! ULON
            enddo
            enddo
            work_g1(:,:) = work_g1(:,:) / rad_to_deg
         endif
         call scatter_global(ULON, work_g1, master_task, distrb_info, &
                             field_loc_NEcorner, field_type_scalar)
         call ice_HaloExtrapolate(ULON, distrb_info, &
                                  ew_boundary_type, ns_boundary_type)

         if (my_task == master_task) then
            work_g1 = c0
            length = dyrect*cm_to_m/radius*rad_to_deg

            work_g1(:,1) = latrefrect ! reference latitude from namelist

            do i = 1, nx_global
            do j = 2, ny_global
               work_g1(i,j) = work_g1(i,j-1) + length   ! ULAT
            enddo
            enddo
            work_g1(:,:) = work_g1(:,:) / rad_to_deg
         endif
         call scatter_global(ULAT, work_g1, master_task, distrb_info, &
                             field_loc_NEcorner, field_type_scalar)
         call ice_HaloExtrapolate(ULAT, distrb_info, &
                                  ew_boundary_type, ns_boundary_type)

         if (my_task == master_task) then
            do j = 1, ny_global
            do i = 1, nx_global
               work_g1(i,j) = dxrect             ! HTN
            enddo
            enddo
         endif
         call primary_grid_lengths_HTN(work_g1)  ! dxU, dxT, dxN, dxE

         if (my_task == master_task) then
            do j = 1, ny_global
            do i = 1, nx_global
               work_g1(i,j) = dyrect             ! HTE
            enddo
            enddo
         endif
         call primary_grid_lengths_HTE(work_g1)  ! dyU, dyT, dyN, dyE

      endif ! scale_dxdy

      !-----------------------------------------------------------------
      ! Construct T-cell land mask
      ! Keyed on ew_boundary_type; ns_boundary_type should be 'open'.
      !-----------------------------------------------------------------

      if (my_task == master_task) then
         work_g1(:,:) = c0      ! initialize hm as land

         if (trim(kmt_type) == 'boxislands') then

            call grid_boxislands_kmt(work_g1)

         elseif (trim(kmt_type) == 'channel') then

            do j = 3,ny_global-2     ! closed top and bottom
            do i = 1,nx_global       ! open sides
               work_g1(i,j) = c1     ! NOTE nx_global > 5
            enddo
            enddo

         elseif (trim(kmt_type) == 'channel_oneeast') then

            do j = ny_global/2,ny_global/2    ! one channel wide
            do i = 1,nx_global       ! open sides
               work_g1(i,j) = c1     ! NOTE nx_global > 5
            enddo
            enddo

         elseif (trim(kmt_type) == 'channel_onenorth') then

            do j = 1,ny_global       ! open sides
            do i = nx_global/2,nx_global/2    ! one channel wide
               work_g1(i,j) = c1     ! NOTE nx_global > 5
            enddo
            enddo

         elseif (trim(kmt_type) == 'wall') then

            do j = 1,ny_global       ! open except
            do i = 1,nx_global-2     ! closed east edge
               work_g1(i,j) = c1
            enddo
            enddo

         elseif (trim(kmt_type) == 'default') then

            ! land in the upper left and lower right corners,
            ! otherwise open boundaries
            imid = nint(aint(real(nx_global)/c2))
            jmid = nint(aint(real(ny_global)/c2))

            do j = 3,ny_global-2
            do i = 3,nx_global-2
               work_g1(i,j) = c1    ! open central domain
            enddo
            enddo

            if (nx_global > 5 .and. ny_global > 5) then

               do j = 1, jmid+2
               do i = 1, imid+2
                  work_g1(i,j) = c1    ! open lower left corner
               enddo
               enddo

               do j = max(jmid-2,1), ny_global
               do i = max(imid-2,1), nx_global
                  work_g1(i,j) = c1    ! open upper right corner
               enddo
               enddo

            endif ! > 5x5 grid

         else

            call abort_ice(subname//' ERROR: unknown kmt_type '//trim(kmt_type), &
                 file=__FILE__, line=__LINE__)

         endif ! kmt_type

         if (close_boundaries) then
            work_g1(:, 1:2) = c0
            work_g1(:, ny_global-1:ny_global) = c0
            work_g1(1:2, :) = c0
            work_g1(nx_global-1:nx_global, :) = c0
         endif

      endif

      call scatter_global(hm, work_g1, master_task, distrb_info, &
                          field_loc_center, field_type_scalar)

      deallocate(work_g1)

      end subroutine rectgrid

!=======================================================================

      subroutine rectgrid_scale_dxdy

        ! generate a variable spaced rectangluar grid.
        ! extend spacing from center of grid outward.

        integer (kind=int_kind) :: &
             i, j, iblk, &
             imid, jmid, &
             center1, center2 ! array centers for expanding dx, dy

        real (kind=dbl_kind) :: &
             length,  &
             rad_to_deg

        real (kind=dbl_kind), dimension(:,:), allocatable :: &
             work_g1

        character(len=*), parameter :: subname = '(rectgrid_scale_dxdy)'

        call icepack_query_parameters(rad_to_deg_out=rad_to_deg)

        allocate(work_g1(nx_global,ny_global))

        ! determine dx spacing
        ! strategy: initialize with dxrect.
        ! if want to scale the grid, work from center outwards,
        ! multplying neighbor cell by scale factor.
        ! this assumes dx varies in x direction only.
        ! (i.e, dx is the same across same y location)
        if (my_task == master_task) then

           ! initialize with initial dxrect
           work_g1(:,:) = dxrect

           ! check if nx is even or odd
           ! if even, middle 2 columns are center
           ! of odd,  middle 1 column is center
           if (mod(nx_global,2) == 0) then ! nx_global is even

              ! with even number of x locatons,
              ! the center two y columns are center
              center1 = nx_global/2  ! integer math
              center2 = center1 + 1  ! integer math

           else ! nx_global = odd
              ! only one center index. set center2=center1
              center1 = ceiling(real(nx_global/2),int_kind)
              center2 = center1
           endif

           ! note loop over only half the x grid points (center1)-1
           ! working from the center outward.
           do j = 1, ny_global
           do i = 1, center1-1
              ! work from center1 to left
              work_g1(center1-i,j) = dxscale*work_g1(center1-i+1,j)

              ! work from center2 to right
              work_g1(center2+i,j) = dxscale*work_g1(center2+i-1,j)
           enddo ! i
           enddo ! j

        endif       ! my_task == master_task


        ! note work_g1 is converted to meters in primary_grid_lengths_HTN
        call primary_grid_lengths_HTN(work_g1)  ! dxU, dxT, dxN, dxE

        ! make ULON array
        if (my_task == master_task) then

           ! make first column reference lon in radians.
           ! the remaining work_g1 is still dx in meters
           work_g1(1,:) = lonrefrect/rad_to_deg ! radians

           ! loop over remaining points and add spacing to successive
           ! x locations
           do j = 1, ny_global
           do i = 2, nx_global ! start from i=2. i=1 is lonrefrect
              length = work_g1(i,j)/radius             ! grid spacing in radians
              work_g1(i,j) = work_g1(i-1,j) + length   ! ULON
           enddo ! i
           enddo ! j
        endif    ! mytask == master_task
        call scatter_global(ULON, work_g1, master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
        call ice_HaloExtrapolate(ULON, distrb_info, &
                                 ew_boundary_type, ns_boundary_type)

        ! determine dy spacing
        ! strategy: initialize with dyrect.
        ! if want to scale the grid, work from center outwards,
        ! multplying neighbor cell by scale factor.
        ! this assumes dy varies in y direction only.
        ! (i.e, dy is the same across same x location)
        if (my_task == master_task) then

           ! initialize with initial dxrect
           work_g1(:,:) = dyrect

           ! check if ny is even or odd
           ! if even, middle 2 rows are center
           ! of odd,  middle 1 row is center
           if (mod(ny_global,2) == 0) then ! ny_global is even

              ! with even number of x locatons,
              ! the center two y columns are center
              center1 = ny_global/2  ! integer math
              center2 = center1 + 1  ! integer math

           else ! ny_global = odd
              ! only one center index. set center2=center1
              center1 = ceiling(real(ny_global/2),int_kind)
              center2 = center1
           endif

           ! note loop over only half the y grid points (center1)-1
           ! working from the center outward.
           do i = 1, nx_global
           do j = 1, center1-1
              ! work from center1 to bottom
              work_g1(i,center1-j) = dyscale*work_g1(i,center1-j+1)

              ! work from center2 to top
              work_g1(i,center2+j) = dyscale*work_g1(i,center2+j-1)
           enddo ! i
           enddo ! j
        endif    ! mytask == master_task
        ! note work_g1 is converted to meters primary_grid_lengths_HTE
        call primary_grid_lengths_HTE(work_g1)  ! dyU, dyT, dyN, dyE

        ! make ULAT array
        if (my_task == master_task) then

           ! make first row reference lat in radians.
           ! the remaining work_g1 is still dy in meters
           work_g1(:,1) = latrefrect/rad_to_deg ! radians


           ! loop over remaining points and add spacing to successive
           ! x locations
           do j = 2, ny_global ! start from j=2. j=1 is latrefrect
           do i = 1, nx_global
              length = work_g1(i,j)/radius             ! grid spacing in radians
              work_g1(i,j) = work_g1(i,j-1) + length   ! ULAT
           enddo ! i
           enddo ! j
        endif    ! mytask == master_task
        call scatter_global(ULAT, work_g1, master_task, distrb_info, &
                            field_loc_NEcorner, field_type_scalar)
        call ice_HaloExtrapolate(ULAT, distrb_info, &
                                 ew_boundary_type, ns_boundary_type)


        deallocate(work_g1)

      end subroutine rectgrid_scale_dxdy

!=======================================================================
      ! Complex land mask for testing box cases
      ! Requires nx_global, ny_global > 20
      ! Assumes work array has been initialized to 1 (ocean) and north and
      ! south land boundaries have been applied (ew_boundary_type='cyclic')

      subroutine grid_boxislands_kmt (work)

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: work

      integer (kind=int_kind) :: &
         i, j, k, & ! indices
         nxb, nyb   ! convenient cell-block sizes for building the mask

      character(len=*), parameter :: subname = '(grid_boxislands_kmt)'

      ! number of cells in 5% of global grid x and y lengths
      nxb = int(real(nx_global, dbl_kind) / c20, int_kind)
      nyb = int(real(ny_global, dbl_kind) / c20, int_kind)

      if (nxb < 1 .or. nyb < 1) &
         call abort_ice(subname//' ERROR: requires larger grid size', &
              file=__FILE__, line=__LINE__)

      ! initialize work area as all ocean (c1).
      work(:,:) = c1

      ! now add land points (c0)
      ! northeast triangle
      k = 0
      do j = ny_global, ny_global-3*nyb, -1
         k = k+1
         do i = nx_global-3*nxb+k, nx_global
            work(i,j) = c0
         enddo
      enddo

      ! northwest docks
      do j = ny_global-3*nyb, ny_global
         do i = 1, 1
            work(i,j) = c0
         enddo
      enddo
      do i = 1, 2*nxb
         do j = ny_global-3*nyb, ny_global-nyb-2
            work(i,j) = c0
         enddo
         do j = ny_global-nyb, ny_global-nyb+1
            work(i,j) = c0
         enddo
      enddo

      ! southwest docks
      do j = 2*nyb, 3*nyb
         do i = 1, 1
            work(i,j) = c0
         enddo
      enddo
      do j = 1, 2*nyb
         do i = 2, nxb
            work(i,j) = c0
         enddo
         do i = 2*nxb-1, 2*nxb
            work(i,j) = c0
         enddo
         do i = 2*nxb+2,4*nxb
            work(i,j) = c0
         enddo
      enddo

      ! tiny island
      do j = 14*nyb, 14*nyb+1
         do i = 14*nxb, 14*nxb+1
            work(i,j) = c0
         enddo
      enddo

      ! X islands
      ! left triangle
      k = 0
      do i = 2*nxb, 4*nxb
         k=k+1
         do j = 10*nyb+k, 14*nyb-k
            work(i,j) = c0
         enddo
      enddo
      ! upper triangle
      k = 0
      do j = 14*nyb, 12*nyb, -1
         k=k+1
         do i = 2*nxb+2+k, 6*nxb-2-k
            work(i,j) = c0
         enddo
      enddo
      ! diagonal
      k = 0
      do j = 10*nyb, 14*nyb
         k=k+1
         do i = 2*nxb+4+k, 2*nxb+6+k
            work(i,j) = c0
         enddo
      enddo
      ! lower right triangle
      k = 0
      do j = 12*nyb, 10*nyb, -1
         k=k+1
         do i = 5*nxb+k, 8*nxb
            work(i,j) = c0
         enddo
      enddo

      ! bar islands
      do i = 10*nxb, 16*nxb
         do j = 4*nyb, 5*nyb
            work(i,j) = c0
         enddo
         do j = 6*nyb+2, 8*nyb
            work(i,j) = c0
         enddo
         do j = 8*nyb+2, 8*nyb+3
            work(i,j) = c0
         enddo
      enddo

      end subroutine grid_boxislands_kmt


!=======================================================================
! Calculate dxU and dxT from HTN on the global grid, to preserve
! ghost cell and/or land values that might otherwise be lost. Scatter
! dxU, dxT and HTN to all processors.
!
! author: Elizabeth C. Hunke, LANL

      subroutine primary_grid_lengths_HTN(work_g)

      real (kind=dbl_kind), dimension(:,:) :: work_g ! global array holding HTN

      ! local variables

      integer (kind=int_kind) :: &
         i, j, &
         ip1 , &     ! i+1
         ierr

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g2

      character(len=*), parameter :: subname = '(primary_grid_lengths_HTN)'

      if (my_task == master_task) then
         allocate(work_g2(nx_global,ny_global), stat=ierr)
      else
         allocate(work_g2(1,1), stat=ierr)
      endif
      if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory', file=__FILE__, line=__LINE__)

      ! HTN, dxU = average of 2 neighbor HTNs in i

      if (my_task == master_task) then
         do j = 1, ny_global
         do i = 1, nx_global
            work_g(i,j) = work_g(i,j) * cm_to_m                ! HTN
         enddo
         enddo
         do j = 1, ny_global
         do i = 1, nx_global
            ! assume cyclic; noncyclic will be handled during scatter
            ip1 = i+1
            if (i == nx_global) ip1 = 1
            work_g2(i,j) = p5*(work_g(i,j) + work_g(ip1,j))    ! dxU
         enddo
         enddo
         if (save_ghte_ghtn) then
            do j = 1, ny_global
            do i = 1,nx_global
               G_HTN(i+nghost,j+nghost) = work_g(i,j)
            enddo
            enddo
            call global_ext_halo(G_HTN)
         endif
      endif
      call scatter_global(HTN, work_g, master_task, distrb_info, &
                          field_loc_Nface, field_type_scalar)
      call scatter_global(dxU, work_g2, master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)

      ! dxT = average of 2 neighbor HTNs in j

      if (my_task == master_task) then
         do j = 2, ny_global
         do i = 1, nx_global
            work_g2(i,j) = p5*(work_g(i,j) + work_g(i,j-1)) ! dxT
         enddo
         enddo
         ! extrapolate to obtain dxT along j=1
         do i = 1, nx_global
            work_g2(i,1) = c2*work_g(i,2) - work_g(i,3) ! dxT
         enddo
      endif
      call scatter_global(dxT, work_g2, master_task, distrb_info, &
                          field_loc_center, field_type_scalar)

      ! dxN = HTN

      dxN(:,:,:) = HTN(:,:,:)   ! dxN

      ! dxE = average of 4 surrounding HTNs

      if (my_task == master_task) then
         do j = 2, ny_global
         do i = 1, nx_global
            ! assume cyclic; noncyclic will be handled during scatter
            ip1 = i+1
            if (i == nx_global) ip1 = 1
            work_g2(i,j) = p25*(work_g(i,j)+work_g(ip1,j)+work_g(i,j-1)+work_g(ip1,j-1))   ! dxE
         enddo
         enddo
         ! extrapolate to obtain dxT along j=1
         do i = 1, nx_global
            ! assume cyclic; noncyclic will be handled during scatter
            ip1 = i+1
            if (i == nx_global) ip1 = 1
            work_g2(i,1) = p5*(c2*work_g(i  ,2) - work_g(i  ,3) + &
                               c2*work_g(ip1,2) - work_g(ip1,3))      ! dxE
         enddo
      endif
      call scatter_global(dxE, work_g2, master_task, distrb_info, &
                          field_loc_center, field_type_scalar)

      deallocate(work_g2, stat=ierr)
      if (ierr/=0) call abort_ice(subname//' ERROR: Dealloc error', file=__FILE__, line=__LINE__)

      end subroutine primary_grid_lengths_HTN

!=======================================================================
! Calculate dyU and dyT from HTE on the global grid, to preserve
! ghost cell and/or land values that might otherwise be lost. Scatter
! dyU, dyT and HTE to all processors.
!
! author: Elizabeth C. Hunke, LANL

      subroutine primary_grid_lengths_HTE(work_g)

      real (kind=dbl_kind), dimension(:,:) :: work_g ! global array holding HTE

      ! local variables

      integer (kind=int_kind) :: &
         i, j, &
         im1, &     ! i-1
         ierr

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g2

      character(len=*), parameter :: subname = '(primary_grid_lengths_HTE)'

      if (my_task == master_task) then
         allocate(work_g2(nx_global,ny_global), stat=ierr)
      else
         allocate(work_g2(1,1), stat=ierr)
      endif
      if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory', file=__FILE__, line=__LINE__)

      ! HTE, dyU = average of 2 neighbor HTE in j

      if (my_task == master_task) then
         do j = 1, ny_global
         do i = 1, nx_global
            work_g(i,j) = work_g(i,j) * cm_to_m                ! HTE
         enddo
         enddo
         do j = 1, ny_global-1
         do i = 1, nx_global
            work_g2(i,j) = p5*(work_g(i,j) + work_g(i,j+1)) ! dyU
         enddo
         enddo
         ! extrapolate to obtain dyU along j=ny_global
         if (ny_global > 1) then
            do i = 1, nx_global
               work_g2(i,ny_global) = c2*work_g(i,ny_global-1) - work_g(i,ny_global-2)  ! dyU
            enddo
         endif
         if (save_ghte_ghtn) then
            do j = 1, ny_global
            do i = 1, nx_global
               G_HTE(i+nghost,j+nghost) = work_g(i,j)
            enddo
            enddo
            call global_ext_halo(G_HTE)
         endif
      endif
      call scatter_global(HTE, work_g, master_task, distrb_info, &
                          field_loc_Eface, field_type_scalar)
      call scatter_global(dyU, work_g2, master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)

      ! dyT = average of 2 neighbor HTE in i

      if (my_task == master_task) then
         do j = 1, ny_global
         do i = 1, nx_global
            ! assume cyclic; noncyclic will be handled during scatter
            im1 = i-1
            if (i == 1) im1 = nx_global
            work_g2(i,j) = p5*(work_g(i,j) + work_g(im1,j))    ! dyT
         enddo
         enddo
      endif
      call scatter_global(dyT, work_g2, master_task, distrb_info, &
                          field_loc_center, field_type_scalar)

      ! dyN = average of 4 neighbor HTEs

      if (my_task == master_task) then
         do j = 1, ny_global-1
         do i = 1, nx_global
            ! assume cyclic; noncyclic will be handled during scatter
            im1 = i-1
            if (i == 1) im1 = nx_global
            work_g2(i,j) = p25*(work_g(i,j) + work_g(im1,j) + work_g(i,j+1) + work_g(im1,j+1))   ! dyN
         enddo
         enddo
         ! extrapolate to obtain dyN along j=ny_global
         if (ny_global > 1) then
            do i = 1, nx_global
               ! assume cyclic; noncyclic will be handled during scatter
               im1 = i-1
               if (i == 1) im1 = nx_global
               work_g2(i,ny_global) = p5*(c2*work_g(i  ,ny_global-1) - work_g(i  ,ny_global-2) + &
                                          c2*work_g(im1,ny_global-1) - work_g(im1,ny_global-2))     ! dyN
            enddo
         endif
      endif
      call scatter_global(dyN, work_g2, master_task, distrb_info, &
                          field_loc_center, field_type_scalar)

      ! dyE = HTE

      dyE(:,:,:) = HTE(:,:,:)

      deallocate(work_g2, stat=ierr)
      if (ierr/=0) call abort_ice(subname//' ERROR: Dealloc errro', file=__FILE__, line=__LINE__)

      end subroutine primary_grid_lengths_HTE

!=======================================================================
!  This subroutine fills ghost cells in global extended grid

      subroutine global_ext_halo(array)

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         array   ! extended global grid size nx+2*nghost, ny+2*nghost
                 ! nghost+1:nghost+nx_global and nghost+1:nghost+ny_global filled on entry

      integer (kind=int_kind) :: n

      character(len=*), parameter :: subname = '(global_ext_halo)'

      do n = 1,nghost
         if (ns_boundary_type =='cyclic') then
            array(:,n)                  = array(:,ny_global+n)
            array(:,ny_global+nghost+n) = array(:,nghost+n)
         elseif (ns_boundary_type == 'open') then
            array(:,n)                  = array(:,nghost+1)
            array(:,ny_global+nghost+n) = array(:,ny_global+nghost)
         else
            array(:,n)                  = c0
            array(:,ny_global+nghost+n) = c0
         endif
      enddo

      do n = 1,nghost
         if (ew_boundary_type =='cyclic') then
            array(n                 ,:) = array(nx_global+n,:)
            array(nx_global+nghost+n,:) = array(nghost+n   ,:)
         elseif (ew_boundary_type == 'open') then
            array(n                 ,:) = array(nghost+1        ,:)
            array(nx_global+nghost+n,:) = array(nx_global+nghost,:)
         else
            array(n                 ,:) = c0
            array(nx_global+nghost+n,:) = c0
         endif
      enddo

      end subroutine global_ext_halo

!=======================================================================
! Sets the boundary values for the T cell land mask (hm) and
! makes the logical land masks for T and U cells (tmask, umask)
! and N and E cells (nmask, emask).
! Also creates hemisphere masks (mask-n northern, mask-s southern)
!
! author: Elizabeth C. Hunke, LANL

      subroutine makemask

      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi, &      ! beginning and end of physical domain
         ierr

      real (kind=dbl_kind) :: &
         puny

      real (kind=dbl_kind), dimension(:,:,:), allocatable :: &
            uvmCD

      type (block) :: &
         this_block           ! block information for current block

      character(len=*), parameter :: subname = '(makemask)'

      call icepack_query_parameters(puny_out=puny)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (kmt,              halo_info, &
                           field_loc_center, field_type_scalar)
      call ice_HaloUpdate (hm,               halo_info, &
                           field_loc_center, field_type_scalar)
      call ice_timer_stop(timer_bound)

      !-----------------------------------------------------------------
      ! construct T-cell and U-cell masks
      !-----------------------------------------------------------------

      bm = c0
      allocate(uvmCD(nx_block,ny_block,max_blocks), stat=ierr)
      if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory', file=__FILE__, line=__LINE__)

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            uvm(i,j,iblk) = min (hm(i,j,  iblk), hm(i+1,j,  iblk), &
                                 hm(i,j+1,iblk), hm(i+1,j+1,iblk))
            npm(i,j,iblk) = min (hm(i,j,  iblk), hm(i,j+1,iblk))
            epm(i,j,iblk) = min (hm(i,j,  iblk), hm(i+1,j,iblk))
            bm(i,j,iblk) = my_task + iblk/100.0_dbl_kind
            uvmCD(i,j,iblk) = (hm(i,j,  iblk)+hm(i+1,j,  iblk) &
                            +  hm(i,j+1,iblk)+hm(i+1,j+1,iblk))
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (uvm,                halo_info, &
                           field_loc_NEcorner, field_type_scalar)
      call ice_HaloUpdate (uvmCD,              halo_info, &
                           field_loc_NEcorner, field_type_scalar)
      call ice_HaloUpdate (npm,                halo_info, &
                           field_loc_Nface,    field_type_scalar)
      call ice_HaloUpdate (epm,                halo_info, &
                           field_loc_Eface,    field_type_scalar)
      call ice_HaloUpdate (bm,                 halo_info, &
                           field_loc_center,   field_type_scalar)
      call ice_HaloUpdate (ocn_gridcell_frac,  halo_info, &
                           field_loc_center,   field_type_scalar)
      call ice_timer_stop(timer_bound)

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         ! needs to cover halo (no halo update for logicals)
         tmask(:,:,iblk)   = .false.
         umask(:,:,iblk)   = .false.
         umaskCD(:,:,iblk) = .false.
         nmask(:,:,iblk)   = .false.
         emask(:,:,iblk)   = .false.
         opmask(:,:,iblk)  = .false.
         do j = jlo-nghost, jhi+nghost
         do i = ilo-nghost, ihi+nghost
            if ( hm(i,j,iblk)   > p5  ) tmask  (i,j,iblk)   = .true.
            if (uvm(i,j,iblk)   > p5  ) umask  (i,j,iblk)   = .true.
            if (uvmCD(i,j,iblk) > c1p5) umaskCD(i,j,iblk)   = .true.
            if (npm(i,j,iblk)   > p5  ) nmask  (i,j,iblk)   = .true.
            if (epm(i,j,iblk)   > p5  ) emask  (i,j,iblk)   = .true.
            if (ocn_gridcell_frac(i,j,iblk) > puny .and. .not. tmask(i,j,iblk)) &
                opmask(i,j,iblk) = .true.
         enddo
         enddo

      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

      !-----------------------------------------------------------------
      ! create hemisphere masks
      !-----------------------------------------------------------------

         lmask_n(:,:,iblk) = .false.
         lmask_s(:,:,iblk) = .false.

         tarean(:,:,iblk) = c0
         tareas(:,:,iblk) = c0

         do j = jlo,jhi
         do i = ilo,ihi

            if (ULAT(i,j,iblk) >= -puny) then
               lmask_n(i,j,iblk) = .true. ! N. Hem.
            else
               lmask_s(i,j,iblk) = .true. ! S. Hem.
            endif

            ! N hemisphere area mask (m^2)
            if (lmask_n(i,j,iblk)) tarean(i,j,iblk) = tarea(i,j,iblk) &
                                                    * hm(i,j,iblk)

            ! S hemisphere area mask (m^2)
            if (lmask_s(i,j,iblk)) tareas(i,j,iblk) = tarea(i,j,iblk) &
                                                    * hm(i,j,iblk)

         enddo
         enddo

      enddo  ! iblk
      !$OMP END PARALLEL DO

      deallocate(uvmCD, stat=ierr)
      if (ierr/=0) call abort_ice(subname//' ERROR: Dealloc errro', file=__FILE__, line=__LINE__)

      end subroutine makemask

!=======================================================================
! Initializes latitude and longitude on T grid
!
! author: Elizabeth C. Hunke, LANL; code originally based on POP grid
! generation routine

      subroutine Tlatlon

      integer (kind=int_kind) :: &
           i, j, iblk       , & ! horizontal indices
           ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      real (kind=dbl_kind) :: &
           z1,x1,y1,z2,x2,y2,z3,x3,y3,z4,x4,y4,tx,ty,tz,da, &
           rad_to_deg

      type (block) :: &
           this_block           ! block information for current block

      character(len=*), parameter :: subname = '(Tlatlon)'

      if (my_task==master_task) then
         write(nu_diag,*) subname,' called'
      endif

      call icepack_query_parameters(rad_to_deg_out=rad_to_deg)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      TLAT(:,:,:) = c0
      TLON(:,:,:) = c0

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block, &
      !$OMP                     x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4, &
      !$OMP                     tx,ty,tz,da)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi

            z1 = cos(ULAT(i-1,j-1,iblk))
            x1 = cos(ULON(i-1,j-1,iblk))*z1
            y1 = sin(ULON(i-1,j-1,iblk))*z1
            z1 = sin(ULAT(i-1,j-1,iblk))

            z2 = cos(ULAT(i,j-1,iblk))
            x2 = cos(ULON(i,j-1,iblk))*z2
            y2 = sin(ULON(i,j-1,iblk))*z2
            z2 = sin(ULAT(i,j-1,iblk))

            z3 = cos(ULAT(i-1,j,iblk))
            x3 = cos(ULON(i-1,j,iblk))*z3
            y3 = sin(ULON(i-1,j,iblk))*z3
            z3 = sin(ULAT(i-1,j,iblk))

            z4 = cos(ULAT(i,j,iblk))
            x4 = cos(ULON(i,j,iblk))*z4
            y4 = sin(ULON(i,j,iblk))*z4
            z4 = sin(ULAT(i,j,iblk))

            ! ---------
            ! TLON/TLAT 4 pt computation (pts 1, 2, 3, 4)
            ! ---------

            tx = (x1+x2+x3+x4)/c4
            ty = (y1+y2+y3+y4)/c4
            tz = (z1+z2+z3+z4)/c4
            da = sqrt(tx**2+ty**2+tz**2)

            tz = tz/da

            ! TLON in radians East
            TLON(i,j,iblk) = c0
            if (tx /= c0 .or. ty /= c0) TLON(i,j,iblk) = atan2(ty,tx)

            ! TLAT in radians North
            TLAT(i,j,iblk) = asin(tz)

         enddo                  ! i
         enddo                  ! j
      enddo                     ! iblk
      !$OMP END PARALLEL DO

      if (trim(grid_type) == 'regional' .or. &
          trim(grid_type) == 'rectangular') then
         ! for W boundary extrapolate from interior
         !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
         do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            i = ilo
            if (this_block%i_glob(i) == 1) then
               do j = jlo, jhi
                  TLON(i,j,iblk) = c2*TLON(i+1,j,iblk) - &
                                      TLON(i+2,j,iblk)
                  TLAT(i,j,iblk) = c2*TLAT(i+1,j,iblk) - &
                                      TLAT(i+2,j,iblk)
               enddo
            endif
         enddo
         !$OMP END PARALLEL DO
      endif   ! regional

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (TLON,             halo_info, &
                           field_loc_center, field_type_scalar, &
                           fillValue=c1)
      call ice_HaloUpdate (TLAT,             halo_info, &
                           field_loc_center, field_type_scalar, &
                           fillValue=c1)
      call ice_HaloExtrapolate(TLON, distrb_info, &
                               ew_boundary_type, ns_boundary_type)
      call ice_HaloExtrapolate(TLAT, distrb_info, &
                               ew_boundary_type, ns_boundary_type)

      end subroutine Tlatlon

!=======================================================================
! Initializes latitude and longitude on N, E grid
!
! author: T. Craig from Tlatlon

      subroutine NElatlon

      use ice_constants, only: c0, c1, c1p5, c2, c4, p5, &
          field_loc_center, field_loc_Nface, field_loc_Eface, &
          field_type_scalar

      integer (kind=int_kind) :: &
           i, j, iblk       , & ! horizontal indices
           ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      real (kind=dbl_kind) :: &
           z1,x1,y1,z2,x2,y2,z3,x3,y3,z4,x4,y4,tx,ty,tz,da, &
           rad_to_deg

      type (block) :: &
           this_block           ! block information for current block

      character(len=*), parameter :: subname = '(NElatlon)'

      if (my_task==master_task) then
         write(nu_diag,*) subname,' called'
      endif

      call icepack_query_parameters(rad_to_deg_out=rad_to_deg)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      NLAT(:,:,:) = c0
      NLON(:,:,:) = c0
      ELAT(:,:,:) = c0
      ELON(:,:,:) = c0

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block, &
      !$OMP                     x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4, &
      !$OMP                     tx,ty,tz,da)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi

            z1 = cos(ULAT(i-1,j-1,iblk))
            x1 = cos(ULON(i-1,j-1,iblk))*z1
            y1 = sin(ULON(i-1,j-1,iblk))*z1
            z1 = sin(ULAT(i-1,j-1,iblk))

            z2 = cos(ULAT(i,j-1,iblk))
            x2 = cos(ULON(i,j-1,iblk))*z2
            y2 = sin(ULON(i,j-1,iblk))*z2
            z2 = sin(ULAT(i,j-1,iblk))

            z3 = cos(ULAT(i-1,j,iblk))
            x3 = cos(ULON(i-1,j,iblk))*z3
            y3 = sin(ULON(i-1,j,iblk))*z3
            z3 = sin(ULAT(i-1,j,iblk))

            z4 = cos(ULAT(i,j,iblk))
            x4 = cos(ULON(i,j,iblk))*z4
            y4 = sin(ULON(i,j,iblk))*z4
            z4 = sin(ULAT(i,j,iblk))

            ! ---------
            ! NLON/NLAT 2 pt computation (pts 3, 4)
            ! ---------

            tx = (x3+x4)/c2
            ty = (y3+y4)/c2
            tz = (z3+z4)/c2
            da = sqrt(tx**2+ty**2+tz**2)

            tz = tz/da

            ! NLON in radians East
            NLON(i,j,iblk) = c0
            if (tx /= c0 .or. ty /= c0) NLON(i,j,iblk) = atan2(ty,tx)

            ! NLAT in radians North
            NLAT(i,j,iblk) = asin(tz)

            ! ---------
            ! ELON/ELAT 2 pt computation (pts 2, 4)
            ! ---------

            tx = (x2+x4)/c2
            ty = (y2+y4)/c2
            tz = (z2+z4)/c2
            da = sqrt(tx**2+ty**2+tz**2)

            tz = tz/da

            ! ELON in radians East
            ELON(i,j,iblk) = c0
            if (tx /= c0 .or. ty /= c0) ELON(i,j,iblk) = atan2(ty,tx)

            ! ELAT in radians North
            ELAT(i,j,iblk) = asin(tz)

         enddo                  ! i
         enddo                  ! j
      enddo                     ! iblk
      !$OMP END PARALLEL DO

      if (trim(grid_type) == 'regional' .or. &
          trim(grid_type) == 'rectangular') then
         ! for W boundary extrapolate from interior
         !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
         do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            i = ilo
            if (this_block%i_glob(i) == 1) then
               do j = jlo, jhi
                  NLON(i,j,iblk) = c1p5*TLON(i+1,j,iblk) - &
                                     p5*TLON(i+2,j,iblk)
                  NLAT(i,j,iblk) = c1p5*TLAT(i+1,j,iblk) - &
                                     p5*TLAT(i+2,j,iblk)
               enddo
            endif
         enddo
         !$OMP END PARALLEL DO
      endif   ! regional

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (NLON,             halo_info, &
                           field_loc_Nface,  field_type_scalar, &
                           fillValue=c1)
      call ice_HaloUpdate (NLAT,             halo_info, &
                           field_loc_Nface,  field_type_scalar, &
                           fillValue=c1)
      call ice_HaloUpdate (ELON,             halo_info, &
                           field_loc_Eface,  field_type_scalar, &
                           fillValue=c1)
      call ice_HaloUpdate (ELAT,             halo_info, &
                           field_loc_Eface,  field_type_scalar, &
                           fillValue=c1)
      call ice_HaloExtrapolate(NLON, distrb_info, &
                               ew_boundary_type, ns_boundary_type)
      call ice_HaloExtrapolate(NLAT, distrb_info, &
                               ew_boundary_type, ns_boundary_type)
      call ice_HaloExtrapolate(ELON, distrb_info, &
                               ew_boundary_type, ns_boundary_type)
      call ice_HaloExtrapolate(ELAT, distrb_info, &
                               ew_boundary_type, ns_boundary_type)
      call ice_timer_stop(timer_bound)

      x1 = global_minval(TLON, distrb_info, tmask)
      x2 = global_maxval(TLON, distrb_info, tmask)
      x3 = global_minval(TLAT, distrb_info, tmask)
      x4 = global_maxval(TLAT, distrb_info, tmask)

      y1 = global_minval(ULON, distrb_info, umask)
      y2 = global_maxval(ULON, distrb_info, umask)
      y3 = global_minval(ULAT, distrb_info, umask)
      y4 = global_maxval(ULAT, distrb_info, umask)

      if (my_task==master_task) then
         write(nu_diag,*) ' '
         write(nu_diag,*) subname,' min/max ULON:', y1*rad_to_deg, y2*rad_to_deg
         write(nu_diag,*) subname,' min/max ULAT:', y3*rad_to_deg, y4*rad_to_deg
         write(nu_diag,*) subname,' min/max TLON:', x1*rad_to_deg, x2*rad_to_deg
         write(nu_diag,*) subname,' min/max TLAT:', x3*rad_to_deg, x4*rad_to_deg
      endif                     ! my_task

      x1 = global_minval(NLON, distrb_info, nmask)
      x2 = global_maxval(NLON, distrb_info, nmask)
      x3 = global_minval(NLAT, distrb_info, nmask)
      x4 = global_maxval(NLAT, distrb_info, nmask)

      y1 = global_minval(ELON, distrb_info, emask)
      y2 = global_maxval(ELON, distrb_info, emask)
      y3 = global_minval(ELAT, distrb_info, emask)
      y4 = global_maxval(ELAT, distrb_info, emask)

      if (my_task==master_task) then
         write(nu_diag,*) ' '
         write(nu_diag,*) subname,' min/max NLON:', x1*rad_to_deg, x2*rad_to_deg
         write(nu_diag,*) subname,' min/max NLAT:', x3*rad_to_deg, x4*rad_to_deg
         write(nu_diag,*) subname,' min/max ELON:', y1*rad_to_deg, y2*rad_to_deg
         write(nu_diag,*) subname,' min/max ELAT:', y3*rad_to_deg, y4*rad_to_deg
      endif                     ! my_task

      end subroutine NElatlon

!=======================================================================
! Shifts quantities from one grid to another
! Constructs the shift based on the grid
! NOTE: Input array includes ghost cells that must be updated before
!       calling this routine.
!
! author: T. Craig

      subroutine grid_average_X2Y_base(type,work1,grid1,work2,grid2)

      character(len=*) , intent(in) :: &
         type, grid1, grid2

      real (kind=dbl_kind), intent(in) :: &
         work1(:,:,:)

      real (kind=dbl_kind), intent(out)   :: &
         work2(:,:,:)

      ! local variables

      character(len=16) :: X2Y

      character(len=*), parameter :: subname = '(grid_average_X2Y_base)'

      if (trim(grid1) == trim(grid2)) then
         work2 = work1
      else
         X2Y = trim(grid1)//'2'//trim(grid2)//trim(type)
         call grid_average_X2Y_1(X2Y,work1,work2)
      endif

      end subroutine grid_average_X2Y_base

!=======================================================================
! Shifts quantities from one grid to another
! NOTE: Input array includes ghost cells that must be updated before
!       calling this routine.
!
! author: T. Craig

      subroutine grid_average_X2Y_userwghts(type,work1,grid1,wght1,mask1,work2,grid2)

      character(len=*) , intent(in) :: &
         type, grid1, grid2

      real (kind=dbl_kind), intent(in) :: &
         work1(:,:,:), &
         wght1(:,:,:), &
         mask1(:,:,:)

      real (kind=dbl_kind), intent(out)   :: &
         work2(:,:,:)

      ! local variables

      character(len=16) :: X2Y

      character(len=*), parameter :: subname = '(grid_average_X2Y_userwghts)'

      if (trim(grid1) == trim(grid2)) then
         work2 = work1
      else
         X2Y = trim(grid1)//'2'//trim(grid2)//trim(type)
         call grid_average_X2Y_1f(X2Y,work1,wght1,mask1,work2)
      endif

      end subroutine grid_average_X2Y_userwghts

!=======================================================================
! Shifts quantities from one grid to another
! NOTE: Input array includes ghost cells that must be updated before
!       calling this routine.
!
! author: T. Craig

      subroutine grid_average_X2Y_NEversion(type,work1a,grid1a,work1b,grid1b,work2,grid2)

      character(len=*) , intent(in) :: &
         type, grid1a, grid1b, grid2

      real (kind=dbl_kind), intent(in) :: &
         work1a(:,:,:), work1b(:,:,:)

      real (kind=dbl_kind), intent(out)   :: &
         work2(:,:,:)

      ! local variables

      character(len=16) :: X2Y

      character(len=*), parameter :: subname = '(grid_average_X2Y_NEversion)'

      X2Y = trim(grid1a)//trim(grid1b)//'2'//trim(grid2)//trim(type)

      select case (trim(X2Y))

         ! state masked
         case('NE2US')
            call grid_average_X2Y_2('NE2US',work1a,narea,npm,work1b,earea,epm,work2)
         case('EN2US')
            call grid_average_X2Y_2('NE2US',work1b,narea,npm,work1a,earea,epm,work2)
         case('NE2TS')
            call grid_average_X2Y_2('NE2TS',work1a,narea,npm,work1b,earea,epm,work2)
         case('EN2TS')
            call grid_average_X2Y_2('NE2TS',work1b,narea,npm,work1a,earea,epm,work2)

         ! state unmasked
         case('NE2UA')
            call grid_average_X2Y_2('NE2UA',work1a,narea,npm,work1b,earea,epm,work2)
         case('EN2UA')
            call grid_average_X2Y_2('NE2UA',work1b,narea,npm,work1a,earea,epm,work2)
         case('NE2TA')
            call grid_average_X2Y_2('NE2TA',work1a,narea,npm,work1b,earea,epm,work2)
         case('EN2TA')
            call grid_average_X2Y_2('NE2TA',work1b,narea,npm,work1a,earea,epm,work2)

         case default
            call abort_ice(subname//' ERROR: unknown X2Y '//trim(X2Y), file=__FILE__, line=__LINE__)
      end select

      end subroutine grid_average_X2Y_NEversion

!=======================================================================
! Shifts quantities from one grid to another
! NOTE: Input array includes ghost cells that must be updated before
!       calling this routine.
!
! author: T. Craig

      subroutine grid_average_X2Y_1(X2Y,work1,work2)

      character(len=*) , intent(in) :: &
         X2Y

      real (kind=dbl_kind), intent(in) :: &
         work1(:,:,:)

      real (kind=dbl_kind), intent(out)   :: &
         work2(:,:,:)

      ! local variables

      character(len=*), parameter :: subname = '(grid_average_X2Y_1)'

      select case (trim(X2Y))

         ! flux unmasked
         case('T2UF')
            call grid_average_X2YF('NE',work1,tarea,work2,uarea)
         case('T2EF')
            call grid_average_X2YF('E' ,work1,tarea,work2,earea)
         case('T2NF')
            call grid_average_X2YF('N' ,work1,tarea,work2,narea)
         case('U2TF')
            call grid_average_X2YF('SW',work1,uarea,work2,tarea)
         case('U2EF')
            call grid_average_X2YF('S' ,work1,uarea,work2,earea)
         case('U2NF')
            call grid_average_X2YF('W' ,work1,uarea,work2,narea)
         case('E2TF')
            call grid_average_X2YF('W' ,work1,earea,work2,tarea)
         case('E2UF')
            call grid_average_X2YF('N' ,work1,earea,work2,uarea)
         case('E2NF')
            call grid_average_X2YF('NW',work1,earea,work2,narea)
         case('N2TF')
            call grid_average_X2YF('S' ,work1,narea,work2,tarea)
         case('N2UF')
            call grid_average_X2YF('E' ,work1,narea,work2,uarea)
         case('N2EF')
            call grid_average_X2YF('SE',work1,narea,work2,earea)

         ! state masked
         case('T2US')
            call grid_average_X2YS('NE',work1,tarea,hm ,work2)
         case('T2ES')
            call grid_average_X2YS('E' ,work1,tarea,hm ,work2)
         case('T2NS')
            call grid_average_X2YS('N' ,work1,tarea,hm ,work2)
         case('U2TS')
            call grid_average_X2YS('SW',work1,uarea,uvm,work2)
         case('U2ES')
            call grid_average_X2YS('S' ,work1,uarea,uvm,work2)
         case('U2NS')
            call grid_average_X2YS('W' ,work1,uarea,uvm,work2)
         case('E2TS')
            call grid_average_X2YS('W' ,work1,earea,epm,work2)
         case('E2US')
            call grid_average_X2YS('N' ,work1,earea,epm,work2)
         case('E2NS')
            call grid_average_X2YS('NW',work1,earea,epm,work2)
         case('N2TS')
            call grid_average_X2YS('S' ,work1,narea,npm,work2)
         case('N2US')
            call grid_average_X2YS('E' ,work1,narea,npm,work2)
         case('N2ES')
            call grid_average_X2YS('SE',work1,narea,npm,work2)

         ! state unmasked
         case('T2UA')
            call grid_average_X2YA('NE',work1,tarea,work2)
         case('T2EA')
            call grid_average_X2YA('E' ,work1,tarea,work2)
         case('T2NA')
            call grid_average_X2YA('N' ,work1,tarea,work2)
         case('U2TA')
            call grid_average_X2YA('SW',work1,uarea,work2)
         case('U2EA')
            call grid_average_X2YA('S' ,work1,uarea,work2)
         case('U2NA')
            call grid_average_X2YA('W' ,work1,uarea,work2)
         case('E2TA')
            call grid_average_X2YA('W' ,work1,earea,work2)
         case('E2UA')
            call grid_average_X2YA('N' ,work1,earea,work2)
         case('E2NA')
            call grid_average_X2YA('NW',work1,earea,work2)
         case('N2TA')
            call grid_average_X2YA('S' ,work1,narea,work2)
         case('N2UA')
            call grid_average_X2YA('E' ,work1,narea,work2)
         case('N2EA')
            call grid_average_X2YA('SE',work1,narea,work2)

         case default
            call abort_ice(subname//' ERROR: unknown X2Y '//trim(X2Y), file=__FILE__, line=__LINE__)
      end select

      end subroutine grid_average_X2Y_1

!=======================================================================
! Shifts quantities from one grid to another
! NOTE: Input array includes ghost cells that must be updated before
!       calling this routine.
!
! author: T. Craig

      subroutine grid_average_X2Y_1f(X2Y,work1,wght1,mask1,work2)

      character(len=*) , intent(in) :: &
         X2Y

      real (kind=dbl_kind), intent(in) :: &
         work1(:,:,:), &
         wght1(:,:,:), &
         mask1(:,:,:)

      real (kind=dbl_kind), intent(out)   :: &
         work2(:,:,:)

      ! local variables

      character(len=*), parameter :: subname = '(grid_average_X2Y_1f)'

      select case (trim(X2Y))

! don't support these for now, requires extra destination wght
!         ! flux unmasked
!         case('T2UF')
!            call grid_average_X2YF('NE',work1,tarea,work2,uarea)
!         case('T2EF')
!            call grid_average_X2YF('E' ,work1,tarea,work2,earea)
!         case('T2NF')
!            call grid_average_X2YF('N' ,work1,tarea,work2,narea)
!         case('U2TF')
!            call grid_average_X2YF('SW',work1,uarea,work2,tarea)
!         case('U2EF')
!            call grid_average_X2YF('S' ,work1,uarea,work2,earea)
!         case('U2NF')
!            call grid_average_X2YF('W' ,work1,uarea,work2,narea)
!         case('E2TF')
!            call grid_average_X2YF('W' ,work1,earea,work2,tarea)
!         case('E2UF')
!            call grid_average_X2YF('N' ,work1,earea,work2,uarea)
!         case('E2NF')
!            call grid_average_X2YF('NW',work1,earea,work2,narea)
!         case('N2TF')
!            call grid_average_X2YF('S' ,work1,narea,work2,tarea)
!         case('N2UF')
!            call grid_average_X2YF('E' ,work1,narea,work2,uarea)
!         case('N2EF')
!            call grid_average_X2YF('SE',work1,narea,work2,earea)

         ! state masked
         case('T2US')
            call grid_average_X2YS('NE',work1,wght1,mask1,work2)
         case('T2ES')
            call grid_average_X2YS('E' ,work1,wght1,mask1,work2)
         case('T2NS')
            call grid_average_X2YS('N' ,work1,wght1,mask1,work2)
         case('U2TS')
            call grid_average_X2YS('SW',work1,wght1,mask1,work2)
         case('U2ES')
            call grid_average_X2YS('S' ,work1,wght1,mask1,work2)
         case('U2NS')
            call grid_average_X2YS('W' ,work1,wght1,mask1,work2)
         case('E2TS')
            call grid_average_X2YS('W' ,work1,wght1,mask1,work2)
         case('E2US')
            call grid_average_X2YS('N' ,work1,wght1,mask1,work2)
         case('E2NS')
            call grid_average_X2YS('NW',work1,wght1,mask1,work2)
         case('N2TS')
            call grid_average_X2YS('S' ,work1,wght1,mask1,work2)
         case('N2US')
            call grid_average_X2YS('E' ,work1,wght1,mask1,work2)
         case('N2ES')
            call grid_average_X2YS('SE',work1,wght1,mask1,work2)

         ! state unmasked
         case('T2UA')
            call grid_average_X2YA('NE',work1,wght1,work2)
         case('T2EA')
            call grid_average_X2YA('E' ,work1,wght1,work2)
         case('T2NA')
            call grid_average_X2YA('N' ,work1,wght1,work2)
         case('U2TA')
            call grid_average_X2YA('SW',work1,wght1,work2)
         case('U2EA')
            call grid_average_X2YA('S' ,work1,wght1,work2)
         case('U2NA')
            call grid_average_X2YA('W' ,work1,wght1,work2)
         case('E2TA')
            call grid_average_X2YA('W' ,work1,wght1,work2)
         case('E2UA')
            call grid_average_X2YA('N' ,work1,wght1,work2)
         case('E2NA')
            call grid_average_X2YA('NW',work1,wght1,work2)
         case('N2TA')
            call grid_average_X2YA('S' ,work1,wght1,work2)
         case('N2UA')
            call grid_average_X2YA('E' ,work1,wght1,work2)
         case('N2EA')
            call grid_average_X2YA('SE',work1,wght1,work2)

         case default
            call abort_ice(subname//' ERROR: unknown X2Y '//trim(X2Y), file=__FILE__, line=__LINE__)
      end select

      end subroutine grid_average_X2Y_1f

!=======================================================================
! Shifts quantities from one grid to another
! State masked version, simple area weighted averager
! NOTE: Input array includes ghost cells that must be updated before
!       calling this routine.
!
! author: T. Craig

      subroutine grid_average_X2YS(dir,work1,wght1,mask1,work2)

      character(len=*) , intent(in) :: &
         dir

      real (kind=dbl_kind), intent(in) :: &
         work1(:,:,:), &
         wght1(:,:,:), &
         mask1(:,:,:)

      real (kind=dbl_kind), intent(out) :: &
         work2(:,:,:)

      ! local variables

      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      real (kind=dbl_kind) :: &
         wtmp

      type (block) :: &
         this_block           ! block information for current block

      character(len=*), parameter :: subname = '(grid_average_X2YS)'

      work2(:,:,:) = c0

      select case (trim(dir))

         case('NE')
            !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo, jhi
               do i = ilo, ihi
                  wtmp = (mask1(i  ,j  ,iblk)*wght1(i  ,j  ,iblk)  &
                        + mask1(i+1,j  ,iblk)*wght1(i+1,j  ,iblk)  &
                        + mask1(i  ,j+1,iblk)*wght1(i  ,j+1,iblk)  &
                        + mask1(i+1,j+1,iblk)*wght1(i+1,j+1,iblk))
                  if (wtmp /= c0) &
                  work2(i,j,iblk) = (mask1(i  ,j  ,iblk)*work1(i  ,j  ,iblk)*wght1(i  ,j  ,iblk)  &
                                   + mask1(i+1,j  ,iblk)*work1(i+1,j  ,iblk)*wght1(i+1,j  ,iblk)  &
                                   + mask1(i  ,j+1,iblk)*work1(i  ,j+1,iblk)*wght1(i  ,j+1,iblk)  &
                                   + mask1(i+1,j+1,iblk)*work1(i+1,j+1,iblk)*wght1(i+1,j+1,iblk)) &
                                   / wtmp
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

         case('SW')
            !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo, jhi
               do i = ilo, ihi
                  wtmp = (mask1(i  ,j  ,iblk)*wght1(i  ,j  ,iblk)  &
                        + mask1(i-1,j  ,iblk)*wght1(i-1,j  ,iblk)  &
                        + mask1(i  ,j-1,iblk)*wght1(i  ,j-1,iblk)  &
                        + mask1(i-1,j-1,iblk)*wght1(i-1,j-1,iblk))
                  if (wtmp /= c0) &
                  work2(i,j,iblk) = (mask1(i  ,j  ,iblk)*work1(i  ,j  ,iblk)*wght1(i  ,j  ,iblk)  &
                                   + mask1(i-1,j  ,iblk)*work1(i-1,j  ,iblk)*wght1(i-1,j  ,iblk)  &
                                   + mask1(i  ,j-1,iblk)*work1(i  ,j-1,iblk)*wght1(i  ,j-1,iblk)  &
                                   + mask1(i-1,j-1,iblk)*work1(i-1,j-1,iblk)*wght1(i-1,j-1,iblk)) &
                                   / wtmp
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

         case('NW')
            !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo, jhi
               do i = ilo, ihi
                  wtmp = (mask1(i-1,j  ,iblk)*wght1(i-1,j  ,iblk)  &
                        + mask1(i  ,j  ,iblk)*wght1(i  ,j  ,iblk)  &
                        + mask1(i-1,j+1,iblk)*wght1(i-1,j+1,iblk)  &
                        + mask1(i  ,j+1,iblk)*wght1(i  ,j+1,iblk))
                  if (wtmp /= c0) &
                  work2(i,j,iblk) = (mask1(i-1,j  ,iblk)*work1(i-1,j  ,iblk)*wght1(i-1,j  ,iblk)  &
                                   + mask1(i  ,j  ,iblk)*work1(i  ,j  ,iblk)*wght1(i  ,j  ,iblk)  &
                                   + mask1(i-1,j+1,iblk)*work1(i-1,j+1,iblk)*wght1(i-1,j+1,iblk)  &
                                   + mask1(i  ,j+1,iblk)*work1(i  ,j+1,iblk)*wght1(i  ,j+1,iblk)) &
                                   / wtmp
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

         case('SE')
            !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo, jhi
               do i = ilo, ihi
                  wtmp = (mask1(i  ,j-1,iblk)*wght1(i  ,j-1,iblk)  &
                        + mask1(i+1,j-1,iblk)*wght1(i+1,j-1,iblk)  &
                        + mask1(i  ,j  ,iblk)*wght1(i  ,j  ,iblk)  &
                        + mask1(i+1,j  ,iblk)*wght1(i+1,j  ,iblk))
                  if (wtmp /= c0) &
                  work2(i,j,iblk) = (mask1(i  ,j-1,iblk)*work1(i  ,j-1,iblk)*wght1(i  ,j-1,iblk)  &
                                   + mask1(i+1,j-1,iblk)*work1(i+1,j-1,iblk)*wght1(i+1,j-1,iblk)  &
                                   + mask1(i  ,j  ,iblk)*work1(i  ,j  ,iblk)*wght1(i  ,j  ,iblk)  &
                                   + mask1(i+1,j  ,iblk)*work1(i+1,j  ,iblk)*wght1(i+1,j  ,iblk)) &
                                   / wtmp
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

         case('E')
            !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo, jhi
               do i = ilo, ihi
                  wtmp = (mask1(i  ,j,iblk)*wght1(i  ,j,iblk)  &
                        + mask1(i+1,j,iblk)*wght1(i+1,j,iblk))
                  if (wtmp /= c0) &
                  work2(i,j,iblk) = (mask1(i  ,j,iblk)*work1(i  ,j,iblk)*wght1(i  ,j,iblk)  &
                                   + mask1(i+1,j,iblk)*work1(i+1,j,iblk)*wght1(i+1,j,iblk)) &
                                   / wtmp
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

         case('W')
            !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo, jhi
               do i = ilo, ihi
                  wtmp = (mask1(i-1,j,iblk)*wght1(i-1,j,iblk)  &
                        + mask1(i  ,j,iblk)*wght1(i  ,j,iblk))
                  if (wtmp /= c0) &
                  work2(i,j,iblk) = (mask1(i-1,j,iblk)*work1(i-1,j,iblk)*wght1(i-1,j,iblk)  &
                                   + mask1(i  ,j,iblk)*work1(i  ,j,iblk)*wght1(i  ,j,iblk)) &
                                   / wtmp
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

         case('N')
            !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo, jhi
               do i = ilo, ihi
                  wtmp = (mask1(i,j  ,iblk)*wght1(i,j  ,iblk)  &
                        + mask1(i,j+1,iblk)*wght1(i,j+1,iblk))
                  if (wtmp /= c0) &
                  work2(i,j,iblk) = (mask1(i,j  ,iblk)*work1(i,j  ,iblk)*wght1(i,j  ,iblk)  &
                                   + mask1(i,j+1,iblk)*work1(i,j+1,iblk)*wght1(i,j+1,iblk)) &
                                   / wtmp
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

         case('S')
            !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo, jhi
               do i = ilo, ihi
                  wtmp = (mask1(i,j-1,iblk)*wght1(i,j-1,iblk)  &
                        + mask1(i,j  ,iblk)*wght1(i,j  ,iblk))
                  if (wtmp /= c0) &
                  work2(i,j,iblk) = (mask1(i,j-1,iblk)*work1(i,j-1,iblk)*wght1(i,j-1,iblk)  &
                                   + mask1(i,j  ,iblk)*work1(i,j  ,iblk)*wght1(i,j  ,iblk)) &
                                   / wtmp
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

         case default
            call abort_ice(subname//' ERROR: unknown option '//trim(dir), file=__FILE__, line=__LINE__)
      end select

      end subroutine grid_average_X2YS

!=======================================================================
! Shifts quantities from one grid to another
! State unmasked version, simple weighted averager
! NOTE: Input array includes ghost cells that must be updated before
!       calling this routine.
!
! author: T. Craig

      subroutine grid_average_X2YA(dir,work1,wght1,work2)

      character(len=*) , intent(in) :: &
         dir

      real (kind=dbl_kind), intent(in) :: &
         work1(:,:,:), &
         wght1(:,:,:)

      real (kind=dbl_kind), intent(out) :: &
         work2(:,:,:)

      ! local variables

      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      real (kind=dbl_kind) :: &
         wtmp

      type (block) :: &
         this_block           ! block information for current block

      character(len=*), parameter :: subname = '(grid_average_X2YA)'

      work2(:,:,:) = c0

      select case (trim(dir))

         case('NE')
            !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo, jhi
               do i = ilo, ihi
                  wtmp = (wght1(i  ,j  ,iblk)  &
                        + wght1(i+1,j  ,iblk)  &
                        + wght1(i  ,j+1,iblk)  &
                        + wght1(i+1,j+1,iblk))
                  if (wtmp /= c0) &
                  work2(i,j,iblk) = (work1(i  ,j  ,iblk)*wght1(i  ,j  ,iblk)  &
                                   + work1(i+1,j  ,iblk)*wght1(i+1,j  ,iblk)  &
                                   + work1(i  ,j+1,iblk)*wght1(i  ,j+1,iblk)  &
                                   + work1(i+1,j+1,iblk)*wght1(i+1,j+1,iblk)) &
                                   / wtmp
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

         case('SW')
            !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo, jhi
               do i = ilo, ihi
                  wtmp = (wght1(i  ,j  ,iblk)  &
                        + wght1(i-1,j  ,iblk)  &
                        + wght1(i  ,j-1,iblk)  &
                        + wght1(i-1,j-1,iblk))
                  if (wtmp /= c0) &
                  work2(i,j,iblk) = (work1(i  ,j  ,iblk)*wght1(i  ,j  ,iblk)  &
                                   + work1(i-1,j  ,iblk)*wght1(i-1,j  ,iblk)  &
                                   + work1(i  ,j-1,iblk)*wght1(i  ,j-1,iblk)  &
                                   + work1(i-1,j-1,iblk)*wght1(i-1,j-1,iblk)) &
                                   / wtmp
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

         case('NW')
            !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo, jhi
               do i = ilo, ihi
                  wtmp = (wght1(i-1,j  ,iblk)  &
                        + wght1(i  ,j  ,iblk)  &
                        + wght1(i-1,j+1,iblk)  &
                        + wght1(i  ,j+1,iblk))
                  if (wtmp /= c0) &
                  work2(i,j,iblk) = (work1(i-1,j  ,iblk)*wght1(i-1,j  ,iblk)  &
                                   + work1(i  ,j  ,iblk)*wght1(i  ,j  ,iblk)  &
                                   + work1(i-1,j+1,iblk)*wght1(i-1,j+1,iblk)  &
                                   + work1(i  ,j+1,iblk)*wght1(i  ,j+1,iblk)) &
                                   / wtmp
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

         case('SE')
            !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo, jhi
               do i = ilo, ihi
                  wtmp = (wght1(i  ,j-1,iblk)  &
                        + wght1(i+1,j-1,iblk)  &
                        + wght1(i  ,j  ,iblk)  &
                        + wght1(i+1,j  ,iblk))
                  if (wtmp /= c0) &
                  work2(i,j,iblk) = (work1(i  ,j-1,iblk)*wght1(i  ,j-1,iblk)  &
                                   + work1(i+1,j-1,iblk)*wght1(i+1,j-1,iblk)  &
                                   + work1(i  ,j  ,iblk)*wght1(i  ,j  ,iblk)  &
                                   + work1(i+1,j  ,iblk)*wght1(i+1,j  ,iblk)) &
                                   / wtmp
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

         case('E')
            !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo, jhi
               do i = ilo, ihi
                  wtmp = (wght1(i  ,j,iblk)  &
                        + wght1(i+1,j,iblk))
                  if (wtmp /= c0) &
                  work2(i,j,iblk) = (work1(i  ,j,iblk)*wght1(i  ,j,iblk)  &
                                   + work1(i+1,j,iblk)*wght1(i+1,j,iblk)) &
                                   / wtmp
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

         case('W')
            !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo, jhi
               do i = ilo, ihi
                  wtmp = (wght1(i-1,j,iblk)  &
                        + wght1(i  ,j,iblk))
                  if (wtmp /= c0) &
                  work2(i,j,iblk) = (work1(i-1,j,iblk)*wght1(i-1,j,iblk)  &
                                   + work1(i  ,j,iblk)*wght1(i  ,j,iblk)) &
                                   / wtmp
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

         case('N')
            !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo, jhi
               do i = ilo, ihi
                  wtmp = (wght1(i,j  ,iblk)  &
                        + wght1(i,j+1,iblk))
                  if (wtmp /= c0) &
                  work2(i,j,iblk) = (work1(i,j  ,iblk)*wght1(i,j  ,iblk)  &
                                   + work1(i,j+1,iblk)*wght1(i,j+1,iblk)) &
                                   / wtmp
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

         case('S')
            !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo, jhi
               do i = ilo, ihi
                  wtmp = (wght1(i,j-1,iblk)  &
                        + wght1(i,j  ,iblk))
                  if (wtmp /= c0) &
                  work2(i,j,iblk) = (work1(i,j-1,iblk)*wght1(i,j-1,iblk)  &
                                   + work1(i,j  ,iblk)*wght1(i,j  ,iblk)) &
                                   / wtmp
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

         case default
            call abort_ice(subname//' ERROR: unknown option '//trim(dir), file=__FILE__, line=__LINE__)
      end select

      end subroutine grid_average_X2YA

!=======================================================================
! Shifts quantities from one grid to another
! Flux masked, original implementation based on earlier t2u and u2t versions
! NOTE: Input array includes ghost cells that must be updated before
!       calling this routine.
!
! author: T. Craig

      subroutine grid_average_X2YF(dir,work1,wght1,work2,wght2)

      character(len=*) , intent(in) :: &
         dir

      real (kind=dbl_kind), intent(in) :: &
         work1(:,:,:), &
         wght1(:,:,:), &
         wght2(:,:,:)

      real (kind=dbl_kind), intent(out) :: &
         work2(:,:,:)

      ! local variables

      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      type (block) :: &
         this_block           ! block information for current block

      character(len=*), parameter :: subname = '(grid_average_X2YF)'

      work2(:,:,:) = c0

      select case (trim(dir))

         case('NE')
            !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo, jhi
               do i = ilo, ihi
                  work2(i,j,iblk) = p25 * &
                                    (work1(i  ,j  ,iblk)*wght1(i  ,j  ,iblk)  &
                                   + work1(i+1,j  ,iblk)*wght1(i+1,j  ,iblk)  &
                                   + work1(i  ,j+1,iblk)*wght1(i  ,j+1,iblk)  &
                                   + work1(i+1,j+1,iblk)*wght1(i+1,j+1,iblk)) &
                                   / wght2(i  ,j  ,iblk)
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

         case('SW')
            !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo, jhi
               do i = ilo, ihi
                  work2(i,j,iblk) = p25 *  &
                                   (work1(i  ,j  ,iblk)*wght1(i  ,j  ,iblk)  &
                                  + work1(i-1,j  ,iblk)*wght1(i-1,j  ,iblk)  &
                                  + work1(i  ,j-1,iblk)*wght1(i  ,j-1,iblk)  &
                                  + work1(i-1,j-1,iblk)*wght1(i-1,j-1,iblk)) &
                                  / wght2(i  ,j  ,iblk)
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

         case('NW')
            !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo, jhi
               do i = ilo, ihi
                  work2(i,j,iblk) = p25 * &
                                    (work1(i-1,j  ,iblk)*wght1(i-1,j  ,iblk)  &
                                   + work1(i  ,j  ,iblk)*wght1(i  ,j  ,iblk)  &
                                   + work1(i-1,j+1,iblk)*wght1(i-1,j+1,iblk)  &
                                   + work1(i  ,j+1,iblk)*wght1(i  ,j+1,iblk)) &
                                   / wght2(i  ,j  ,iblk)
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

         case('SE')
            !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo, jhi
               do i = ilo, ihi
                  work2(i,j,iblk) = p25 *  &
                                   (work1(i  ,j-1,iblk)*wght1(i  ,j-1,iblk)  &
                                  + work1(i+1,j-1,iblk)*wght1(i+1,j-1,iblk)  &
                                  + work1(i  ,j  ,iblk)*wght1(i  ,j  ,iblk)  &
                                  + work1(i+1,j  ,iblk)*wght1(i+1,j  ,iblk)) &
                                  / wght2(i  ,j  ,iblk)
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

         case('E')
            !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo, jhi
               do i = ilo, ihi
                  work2(i,j,iblk) = p5 * &
                                    (work1(i  ,j,iblk)*wght1(i  ,j,iblk)  &
                                   + work1(i+1,j,iblk)*wght1(i+1,j,iblk)) &
                                   / wght2(i  ,j,iblk)
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

         case('W')
            !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo, jhi
               do i = ilo, ihi
                  work2(i,j,iblk) = p5 * &
                                    (work1(i-1,j,iblk)*wght1(i-1,j,iblk)  &
                                   + work1(i  ,j,iblk)*wght1(i  ,j,iblk)) &
                                   / wght2(i  ,j,iblk)
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

         case('N')
            !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo, jhi
               do i = ilo, ihi
                  work2(i,j,iblk) = p5 * &
                                    (work1(i,j  ,iblk)*wght1(i,j  ,iblk)  &
                                   + work1(i,j+1,iblk)*wght1(i,j+1,iblk)) &
                                   / wght2(i  ,j,iblk)
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

         case('S')
            !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo, jhi
               do i = ilo, ihi
                  work2(i,j,iblk) = p5 * &
                                    (work1(i,j-1,iblk)*wght1(i,j-1,iblk)  &
                                   + work1(i,j  ,iblk)*wght1(i,j  ,iblk)) &
                                   / wght2(i  ,j,iblk)
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

         case default
            call abort_ice(subname//' ERROR: unknown option '//trim(dir), file=__FILE__, line=__LINE__)
      end select

      end subroutine grid_average_X2YF

!=======================================================================
! Shifts quantities from one grid to another
! State masked version, simple weighted averager
! NOTE: Input array includes ghost cells that must be updated before
!       calling this routine.
!
! author: T. Craig

      subroutine grid_average_X2Y_2(dir,work1a,wght1a,mask1a,work1b,wght1b,mask1b,work2)

      character(len=*) , intent(in) :: &
         dir

      real (kind=dbl_kind), intent(in) :: &
         work1a(:,:,:), work1b(:,:,:), &
         wght1a(:,:,:), wght1b(:,:,:), &
         mask1a(:,:,:), mask1b(:,:,:)

      real (kind=dbl_kind), intent(out) :: &
         work2(:,:,:)

      ! local variables

      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      real (kind=dbl_kind) :: &
         wtmp

      type (block) :: &
         this_block           ! block information for current block

      character(len=*), parameter :: subname = '(grid_average_X2Y_2)'

      work2(:,:,:) = c0

      select case (trim(dir))

         case('NE2US')
            !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo, jhi
               do i = ilo, ihi
                  wtmp = (mask1a(i  ,j  ,iblk)*wght1a(i  ,j  ,iblk)  &
                        + mask1a(i+1,j  ,iblk)*wght1a(i+1,j  ,iblk)  &
                        + mask1b(i  ,j  ,iblk)*wght1b(i  ,j  ,iblk)  &
                        + mask1b(i  ,j+1,iblk)*wght1b(i  ,j+1,iblk))
                  if (wtmp /= c0) &
                  work2(i,j,iblk) = (mask1a(i  ,j  ,iblk)*work1a(i  ,j  ,iblk)*wght1a(i  ,j  ,iblk)  &
                                   + mask1a(i+1,j  ,iblk)*work1a(i+1,j  ,iblk)*wght1a(i+1,j  ,iblk)  &
                                   + mask1b(i  ,j  ,iblk)*work1b(i  ,j  ,iblk)*wght1b(i  ,j  ,iblk)  &
                                   + mask1b(i  ,j+1,iblk)*work1b(i  ,j+1,iblk)*wght1b(i  ,j+1,iblk)) &
                                   / wtmp
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

         case('NE2TS')
            !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo, jhi
               do i = ilo, ihi
                  wtmp = (mask1a(i  ,j-1,iblk)*wght1a(i  ,j-1,iblk)  &
                        + mask1a(i  ,j  ,iblk)*wght1a(i  ,j  ,iblk)  &
                        + mask1b(i-1,j  ,iblk)*wght1b(i-1,j  ,iblk)  &
                        + mask1b(i  ,j  ,iblk)*wght1b(i  ,j  ,iblk))
                  if (wtmp /= c0) &
                  work2(i,j,iblk) = (mask1a(i  ,j-1,iblk)*work1a(i  ,j-1,iblk)*wght1a(i  ,j-1,iblk)  &
                                   + mask1a(i  ,j  ,iblk)*work1a(i  ,j  ,iblk)*wght1a(i  ,j  ,iblk)  &
                                   + mask1b(i-1,j  ,iblk)*work1b(i-1,j  ,iblk)*wght1b(i-1,j  ,iblk)  &
                                   + mask1b(i  ,j  ,iblk)*work1b(i  ,j  ,iblk)*wght1b(i  ,j  ,iblk)) &
                                   / wtmp
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

         case('NE2UA')
            !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo, jhi
               do i = ilo, ihi
                  wtmp = (wght1a(i  ,j  ,iblk)  &
                        + wght1a(i+1,j  ,iblk)  &
                        + wght1b(i  ,j  ,iblk)  &
                        + wght1b(i  ,j+1,iblk))
                  if (wtmp /= c0) &
                  work2(i,j,iblk) = (work1a(i  ,j  ,iblk)*wght1a(i  ,j  ,iblk)  &
                                   + work1a(i+1,j  ,iblk)*wght1a(i+1,j  ,iblk)  &
                                   + work1b(i  ,j  ,iblk)*wght1b(i  ,j  ,iblk)  &
                                   + work1b(i  ,j+1,iblk)*wght1b(i  ,j+1,iblk)) &
                                   / wtmp
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

         case('NE2TA')
            !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block,wtmp)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo, jhi
               do i = ilo, ihi
                  wtmp = (wght1a(i  ,j-1,iblk)  &
                        + wght1a(i  ,j  ,iblk)  &
                        + wght1b(i-1,j  ,iblk)  &
                        + wght1b(i  ,j  ,iblk))
                  if (wtmp /= c0) &
                  work2(i,j,iblk) = (work1a(i  ,j-1,iblk)*wght1a(i  ,j-1,iblk)  &
                                   + work1a(i  ,j  ,iblk)*wght1a(i  ,j  ,iblk)  &
                                   + work1b(i-1,j  ,iblk)*wght1b(i-1,j  ,iblk)  &
                                   + work1b(i  ,j  ,iblk)*wght1b(i  ,j  ,iblk)) &
                                   / wtmp
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

         case default
            call abort_ice(subname//' ERROR: unknown option '//trim(dir), file=__FILE__, line=__LINE__)
      end select

      end subroutine grid_average_X2Y_2

!=======================================================================
! Compute the minimum of adjacent values of a field at specific indices,
! depending on the grid location (U, E, N)
!
      real(kind=dbl_kind) function grid_neighbor_min(field, i, j, grid_location) result(mini)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         field    ! field defined at T point

      integer (kind=int_kind), intent(in) :: &
         i, j

      character(len=*), intent(in) :: &
         grid_location ! grid location at which to compute the minumum (U, E, N)

      character(len=*), parameter :: subname = '(grid_neighbor_min)'

      select case (trim(grid_location))
         case('U')
            mini = min(field(i,j), field(i+1,j), field(i,j+1), field(i+1,j+1))
         case('E')
            mini = min(field(i,j), field(i+1,j))
         case('N')
            mini = min(field(i,j), field(i,j+1))
         case default
            call abort_ice(subname // ' unknown grid_location: ' // grid_location, file=__FILE__, line=__LINE__)
      end select

      end function grid_neighbor_min

!=======================================================================
! Compute the maximum of adjacent values of a field at specific indices,
! depending on the grid location (U, E, N)
!
      real(kind=dbl_kind) function grid_neighbor_max(field, i, j, grid_location) result(maxi)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         field    ! field defined at T point

      integer (kind=int_kind), intent(in) :: &
         i, j

      character(len=*), intent(in) :: &
         grid_location ! grid location at which to compute the maximum (U, E, N)


      character(len=*), parameter :: subname = '(grid_neighbor_max)'

      select case (trim(grid_location))
         case('U')
            maxi = max(field(i,j), field(i+1,j), field(i,j+1), field(i+1,j+1))
         case('E')
            maxi = max(field(i,j), field(i+1,j))
         case('N')
            maxi = max(field(i,j), field(i,j+1))
         case default
            call abort_ice(subname // ' unknown grid_location: ' // grid_location, file=__FILE__, line=__LINE__)
      end select

      end function grid_neighbor_max

!=======================================================================
! The following code is used for obtaining the coordinates of the grid
! vertices for CF-compliant netCDF history output. Approximate!
!=======================================================================

! These fields are only used for netcdf history output, and the
! ghost cell values are not needed.
! NOTE:  Extrapolations were used: these fields are approximate!
!
! authors:   A. McLaren, Met Office
!            E. Hunke, LANL

      subroutine gridbox_corners

      integer (kind=int_kind) :: &
         i,j,iblk,icorner,& ! index counters
         ilo,ihi,jlo,jhi, &    ! beginning and end of physical domain
         ierr

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g2

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1

      real (kind=dbl_kind) :: &
         rad_to_deg

      type (block) :: &
         this_block           ! block information for current block

      character(len=*), parameter :: subname = '(gridbox_corners)'

      call icepack_query_parameters(rad_to_deg_out=rad_to_deg)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      !-------------------------------------------------------------
      ! Get coordinates of grid boxes for each block as follows:
      ! (1) SW corner, (2) SE corner, (3) NE corner, (4) NW corner
      !-------------------------------------------------------------

      latu_bounds(:,:,:,:) = c0
      lonu_bounds(:,:,:,:) = c0

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi

            latu_bounds(1,i,j,iblk)=TLAT(i  ,j  ,iblk)*rad_to_deg
            latu_bounds(2,i,j,iblk)=TLAT(i+1,j  ,iblk)*rad_to_deg
            latu_bounds(3,i,j,iblk)=TLAT(i+1,j+1,iblk)*rad_to_deg
            latu_bounds(4,i,j,iblk)=TLAT(i  ,j+1,iblk)*rad_to_deg

            lonu_bounds(1,i,j,iblk)=TLON(i  ,j  ,iblk)*rad_to_deg
            lonu_bounds(2,i,j,iblk)=TLON(i+1,j  ,iblk)*rad_to_deg
            lonu_bounds(3,i,j,iblk)=TLON(i+1,j+1,iblk)*rad_to_deg
            lonu_bounds(4,i,j,iblk)=TLON(i  ,j+1,iblk)*rad_to_deg

         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !----------------------------------------------------------------
      ! extrapolate on global grid to get edge values
      !----------------------------------------------------------------

      if (my_task == master_task) then
         allocate(work_g2(nx_global,ny_global), stat=ierr)
      else
         allocate(work_g2(1,1), stat=ierr)
      endif
      if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory', file=__FILE__, line=__LINE__)

      work1(:,:,:) = latu_bounds(2,:,:,:)
!     work_g2 = c0

      call gather_global(work_g2, work1, master_task, distrb_info)
      if (my_task == master_task) then
         do j = 1, ny_global
            work_g2(nx_global,j) = c2*work_g2(nx_global-1,j) &
                                    - work_g2(nx_global-2,j)
         enddo
      endif
      call scatter_global(work1, work_g2, &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      latu_bounds(2,:,:,:) = work1(:,:,:)

      work1(:,:,:) = latu_bounds(3,:,:,:)
      call gather_global(work_g2, work1, master_task, distrb_info)
      if (my_task == master_task) then
         do i = 1, nx_global
            work_g2(i,ny_global) = c2*work_g2(i,ny_global-1) &
                                    - work_g2(i,ny_global-2)
         enddo
         do j = 1, ny_global
            work_g2(nx_global,j) = c2*work_g2(nx_global-1,j) &
                                    - work_g2(nx_global-2,j)
         enddo
      endif
      call scatter_global(work1, work_g2, &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      latu_bounds(3,:,:,:) = work1(:,:,:)

      work1(:,:,:) = latu_bounds(4,:,:,:)
      call gather_global(work_g2, work1, master_task, distrb_info)
      if (my_task == master_task) then
         do i = 1, nx_global
            work_g2(i,ny_global) = c2*work_g2(i,ny_global-1) &
                                    - work_g2(i,ny_global-2)
         enddo
      endif
      call scatter_global(work1, work_g2, &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      latu_bounds(4,:,:,:) = work1(:,:,:)

      work1(:,:,:) = lonu_bounds(2,:,:,:)
      call gather_global(work_g2, work1, master_task, distrb_info)
      if (my_task == master_task) then
         do j = 1, ny_global
            work_g2(nx_global,j) = c2*work_g2(nx_global-1,j) &
                                    - work_g2(nx_global-2,j)
         enddo
      endif
      call scatter_global(work1, work_g2, &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      lonu_bounds(2,:,:,:) = work1(:,:,:)

      work1(:,:,:) = lonu_bounds(3,:,:,:)
      call gather_global(work_g2, work1, master_task, distrb_info)
      if (my_task == master_task) then
         do i = 1, nx_global
            work_g2(i,ny_global) = c2*work_g2(i,ny_global-1) &
                                    - work_g2(i,ny_global-2)
         enddo
         do j = 1, ny_global
            work_g2(nx_global,j) = c2*work_g2(nx_global-1,j) &
                                    - work_g2(nx_global-2,j)
         enddo
      endif
      call scatter_global(work1, work_g2, &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      lonu_bounds(3,:,:,:) = work1(:,:,:)

      work1(:,:,:) = lonu_bounds(4,:,:,:)
      call gather_global(work_g2, work1, master_task, distrb_info)
      if (my_task == master_task) then
         do i = 1, nx_global
            work_g2(i,ny_global) = c2*work_g2(i,ny_global-1) &
                                    - work_g2(i,ny_global-2)
         enddo
      endif
      call scatter_global(work1, work_g2, &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      lonu_bounds(4,:,:,:) = work1(:,:,:)

      deallocate(work_g2, stat=ierr)
      if (ierr/=0) call abort_ice(subname//' ERROR: Dealloc errro', file=__FILE__, line=__LINE__)

      !----------------------------------------------------------------
      ! Convert longitude to Degrees East >0 for history output
      !----------------------------------------------------------------

      allocate(work_g2(nx_block,ny_block), stat=ierr)  ! not used as global here
      if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory', file=__FILE__, line=__LINE__)
      !OMP fails in this loop
      do iblk = 1, nblocks
         do icorner = 1, 4
            work_g2(:,:) = lont_bounds(icorner,:,:,iblk) + c360
            where (work_g2 > c360) work_g2 = work_g2 - c360
            where (work_g2 < c0 )  work_g2 = work_g2 + c360
            lont_bounds(icorner,:,:,iblk) = work_g2(:,:)
            work_g2(:,:) = lonu_bounds(icorner,:,:,iblk) + c360
            where (work_g2 > c360) work_g2 = work_g2 - c360
            where (work_g2 < c0 )  work_g2 = work_g2 + c360
            lonu_bounds(icorner,:,:,iblk) = work_g2(:,:)
         enddo
      enddo
      deallocate(work_g2, stat=ierr)
      if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory', file=__FILE__, line=__LINE__)

      end subroutine gridbox_corners

!=======================================================================
! The following code is used for obtaining the coordinates of the grid
! vertices for CF-compliant netCDF history output. Approximate!
!=======================================================================

! These fields are only used for netcdf history output, and the
! ghost cell values are not needed.
! NOTE:  Extrapolations were used: these fields are approximate!
!

      subroutine gridbox_edges

      integer (kind=int_kind) :: &
         i,j,iblk,icorner,& ! index counters
         ilo,ihi,jlo,jhi , &   ! beginning and end of physical domain
         ierr

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g2

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1

      real (kind=dbl_kind) :: &
         rad_to_deg

      type (block) :: &
         this_block           ! block information for current block

      character(len=*), parameter :: subname = '(gridbox_edges)'

      call icepack_query_parameters(rad_to_deg_out=rad_to_deg)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      !-------------------------------------------------------------
      ! Get coordinates of grid boxes for each block as follows:
      ! for N pt: (1) W edge, (2) E edge, (3) E edge j+1, (4) W edge j+1
      ! for E pt: (1) S edge, (2) S edge i+1, (3) N edge, i+1 (4) N edge
      !-------------------------------------------------------------

      latn_bounds(:,:,:,:) = c0
      lonn_bounds(:,:,:,:) = c0
      late_bounds(:,:,:,:) = c0
      lone_bounds(:,:,:,:) = c0

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi

            latn_bounds(1,i,j,iblk)=ELAT(i-1,j  ,iblk)*rad_to_deg
            latn_bounds(2,i,j,iblk)=ELAT(i  ,j  ,iblk)*rad_to_deg
            latn_bounds(3,i,j,iblk)=ELAT(i  ,j+1,iblk)*rad_to_deg
            latn_bounds(4,i,j,iblk)=ELAT(i-1,j+1,iblk)*rad_to_deg

            lonn_bounds(1,i,j,iblk)=ELON(i-1,j  ,iblk)*rad_to_deg
            lonn_bounds(2,i,j,iblk)=ELON(i  ,j  ,iblk)*rad_to_deg
            lonn_bounds(3,i,j,iblk)=ELON(i  ,j+1,iblk)*rad_to_deg
            lonn_bounds(4,i,j,iblk)=ELON(i-1,j+1,iblk)*rad_to_deg

            late_bounds(1,i,j,iblk)=NLAT(i  ,j-1,iblk)*rad_to_deg
            late_bounds(2,i,j,iblk)=NLAT(i+1,j-1,iblk)*rad_to_deg
            late_bounds(3,i,j,iblk)=NLAT(i+1,j  ,iblk)*rad_to_deg
            late_bounds(4,i,j,iblk)=NLAT(i  ,j  ,iblk)*rad_to_deg

            lone_bounds(1,i,j,iblk)=NLON(i  ,j-1,iblk)*rad_to_deg
            lone_bounds(2,i,j,iblk)=NLON(i+1,j-1,iblk)*rad_to_deg
            lone_bounds(3,i,j,iblk)=NLON(i+1,j  ,iblk)*rad_to_deg
            lone_bounds(4,i,j,iblk)=NLON(i  ,j  ,iblk)*rad_to_deg

         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !----------------------------------------------------------------
      ! extrapolate on global grid to get edge values
      !----------------------------------------------------------------

      if (my_task == master_task) then
         allocate(work_g2(nx_global,ny_global), stat=ierr)
      else
         allocate(work_g2(1,1), stat=ierr)
      endif
      if (ierr/=0) call abort_ice(subname//' ERROR: Dealloc error', file=__FILE__, line=__LINE__)

      ! latn_bounds

      work1(:,:,:) = latn_bounds(1,:,:,:)
      call gather_global(work_g2, work1, master_task, distrb_info)
      if (my_task == master_task) then
         do j = 1, ny_global
            work_g2(1,j) = c2*work_g2(2,j) &
                            - work_g2(3,j)
         enddo
      endif
      call scatter_global(work1, work_g2, &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      latn_bounds(1,:,:,:) = work1(:,:,:)

      work1(:,:,:) = latn_bounds(3,:,:,:)
      call gather_global(work_g2, work1, master_task, distrb_info)
      if (my_task == master_task) then
         do i = 1, nx_global
            work_g2(i,ny_global) = c2*work_g2(i,ny_global-1) &
                                    - work_g2(i,ny_global-2)
         enddo
      endif
      call scatter_global(work1, work_g2, &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      latn_bounds(3,:,:,:) = work1(:,:,:)

      work1(:,:,:) = latn_bounds(4,:,:,:)
      call gather_global(work_g2, work1, master_task, distrb_info)
      if (my_task == master_task) then
         do i = 1, nx_global
            work_g2(i,ny_global) = c2*work_g2(i,ny_global-1) &
                                    - work_g2(i,ny_global-2)
         enddo
         do j = 1, ny_global
            work_g2(1,j) = c2*work_g2(2,j) &
                            - work_g2(3,j)
         enddo
      endif
      call scatter_global(work1, work_g2, &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      latn_bounds(4,:,:,:) = work1(:,:,:)

      ! lonn_bounds

      work1(:,:,:) = lonn_bounds(1,:,:,:)
      call gather_global(work_g2, work1, master_task, distrb_info)
      if (my_task == master_task) then
         do j = 1, ny_global
            work_g2(1,j) = c2*work_g2(2,j) &
                            - work_g2(3,j)
         enddo
      endif
      call scatter_global(work1, work_g2, &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      lonn_bounds(1,:,:,:) = work1(:,:,:)

      work1(:,:,:) = lonn_bounds(3,:,:,:)
      call gather_global(work_g2, work1, master_task, distrb_info)
      if (my_task == master_task) then
         do i = 1, nx_global
            work_g2(i,ny_global) = c2*work_g2(i,ny_global-1) &
                                    - work_g2(i,ny_global-2)
         enddo
      endif
      call scatter_global(work1, work_g2, &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      lonn_bounds(3,:,:,:) = work1(:,:,:)

      work1(:,:,:) = lonn_bounds(4,:,:,:)
      call gather_global(work_g2, work1, master_task, distrb_info)
      if (my_task == master_task) then
         do i = 1, nx_global
            work_g2(i,ny_global) = c2*work_g2(i,ny_global-1) &
                                    - work_g2(i,ny_global-2)
         enddo
         do j = 1, ny_global
            work_g2(1,j) = c2*work_g2(2,j) &
                            - work_g2(3,j)
         enddo
      endif
      call scatter_global(work1, work_g2, &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      lonn_bounds(4,:,:,:) = work1(:,:,:)

      ! late_bounds

      work1(:,:,:) = late_bounds(1,:,:,:)
      call gather_global(work_g2, work1, master_task, distrb_info)
      if (my_task == master_task) then
         do i = 1, nx_global
            work_g2(i,1) = c2*work_g2(i,2) &
                            - work_g2(i,3)
         enddo
      endif
      call scatter_global(work1, work_g2, &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      late_bounds(1,:,:,:) = work1(:,:,:)

      work1(:,:,:) = late_bounds(2,:,:,:)
      call gather_global(work_g2, work1, master_task, distrb_info)
      if (my_task == master_task) then
         do i = 1, nx_global
            work_g2(i,1) = c2*work_g2(i,2) &
                            - work_g2(i,3)
         enddo
         do j = 1, ny_global
            work_g2(nx_global,j) = c2*work_g2(nx_global-1,j) &
                                    - work_g2(nx_global-2,j)
         enddo
      endif
      call scatter_global(work1, work_g2, &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      late_bounds(2,:,:,:) = work1(:,:,:)

      work1(:,:,:) = late_bounds(3,:,:,:)
      call gather_global(work_g2, work1, master_task, distrb_info)
      if (my_task == master_task) then
         do j = 1, ny_global
            work_g2(nx_global,j) = c2*work_g2(nx_global-1,j) &
                                    - work_g2(nx_global-2,j)
         enddo
      endif
      call scatter_global(work1, work_g2, &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      late_bounds(3,:,:,:) = work1(:,:,:)

      ! lone_bounds

      work1(:,:,:) = lone_bounds(1,:,:,:)
      call gather_global(work_g2, work1, master_task, distrb_info)
      if (my_task == master_task) then
         do i = 1, nx_global
            work_g2(i,1) = c2*work_g2(i,2) &
                            - work_g2(i,3)
         enddo
      endif
      call scatter_global(work1, work_g2, &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      lone_bounds(1,:,:,:) = work1(:,:,:)

      work1(:,:,:) = lone_bounds(2,:,:,:)
      call gather_global(work_g2, work1, master_task, distrb_info)
      if (my_task == master_task) then
         do i = 1, nx_global
            work_g2(i,1) = c2*work_g2(i,2) &
                            - work_g2(i,3)
         enddo
         do j = 1, ny_global
            work_g2(nx_global,j) = c2*work_g2(nx_global-1,j) &
                                    - work_g2(nx_global-2,j)
         enddo
      endif
      call scatter_global(work1, work_g2, &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      lone_bounds(2,:,:,:) = work1(:,:,:)

      work1(:,:,:) = lone_bounds(3,:,:,:)
      call gather_global(work_g2, work1, master_task, distrb_info)
      if (my_task == master_task) then
         do j = 1, ny_global
            work_g2(nx_global,j) = c2*work_g2(nx_global-1,j) &
                                    - work_g2(nx_global-2,j)
         enddo
      endif
      call scatter_global(work1, work_g2, &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      lone_bounds(3,:,:,:) = work1(:,:,:)

      deallocate(work_g2, stat=ierr)
      if (ierr/=0) call abort_ice(subname//' ERROR: Dealloc error', file=__FILE__, line=__LINE__)

      !----------------------------------------------------------------
      ! Convert longitude to Degrees East >0 for history output
      !----------------------------------------------------------------

      allocate(work_g2(nx_block,ny_block), stat=ierr)  ! not used as global here
      if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory', file=__FILE__, line=__LINE__)
      !OMP fails in this loop
      do iblk = 1, nblocks
         do icorner = 1, 4
            work_g2(:,:) = lonn_bounds(icorner,:,:,iblk) + c360
            where (work_g2 > c360) work_g2 = work_g2 - c360
            where (work_g2 < c0 )  work_g2 = work_g2 + c360
            lonn_bounds(icorner,:,:,iblk) = work_g2(:,:)
            work_g2(:,:) = lone_bounds(icorner,:,:,iblk) + c360
            where (work_g2 > c360) work_g2 = work_g2 - c360
            where (work_g2 < c0 )  work_g2 = work_g2 + c360
            lone_bounds(icorner,:,:,iblk) = work_g2(:,:)
         enddo
      enddo
      deallocate(work_g2, stat=ierr)

      end subroutine gridbox_edges

!=======================================================================
! NOTE:  Boundary conditions for fields on NW, SW, SE corners
!        have not been implemented; using NE corner location for all.
!        Extrapolations are also used: these fields are approximate!
!
! authors:   A. McLaren, Met Office
!            E. Hunke, LANL

      subroutine gridbox_verts(work_g,vbounds)

      real (kind=dbl_kind), dimension(:,:), intent(in) :: &
          work_g

      real (kind=dbl_kind), dimension(4,nx_block,ny_block,max_blocks), intent(out) :: &
          vbounds

      integer (kind=int_kind) :: &
         i,j , &                ! index counters
         ierr

      real (kind=dbl_kind) :: &
          rad_to_deg

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g2

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1

      character(len=*), parameter :: subname = '(gridbox_verts)'

      call icepack_query_parameters(rad_to_deg_out=rad_to_deg)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (my_task == master_task) then
         allocate(work_g2(nx_global,ny_global), stat=ierr)
      else
         allocate(work_g2(1,1), stat=ierr)
      endif
      if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory', file=__FILE__, line=__LINE__)

      !-------------------------------------------------------------
      ! Get coordinates of grid boxes for each block as follows:
      ! (1) SW corner, (2) SE corner, (3) NE corner, (4) NW corner
      !-------------------------------------------------------------

      work_g2(:,:) = c0
      if (my_task == master_task) then
         do j = 2, ny_global
         do i = 2, nx_global
            work_g2(i,j) = work_g(i-1,j-1) * rad_to_deg
         enddo
         enddo
         ! extrapolate
         do j = 1, ny_global
            work_g2(1,j) = c2*work_g2(2,j) - work_g2(3,j)
         enddo
         do i = 1, nx_global
            work_g2(i,1) = c2*work_g2(i,2) - work_g2(i,3)
         enddo
      endif
      call scatter_global(work1, work_g2, &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      vbounds(1,:,:,:) = work1(:,:,:)

      work_g2(:,:) = c0
      if (my_task == master_task) then
         do j = 2, ny_global
         do i = 1, nx_global
            work_g2(i,j) = work_g(i,j-1) * rad_to_deg
         enddo
         enddo
         ! extrapolate
         do i = 1, nx_global
            work_g2(i,1) = (c2*work_g2(i,2) - work_g2(i,3))
         enddo
      endif
      call scatter_global(work1, work_g2, &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      vbounds(2,:,:,:) = work1(:,:,:)

      work_g2(:,:) = c0
      if (my_task == master_task) then
         do j = 1, ny_global
         do i = 1, nx_global
            work_g2(i,j) = work_g(i,j) * rad_to_deg
         enddo
         enddo
      endif
      call scatter_global(work1, work_g2, &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      vbounds(3,:,:,:) = work1(:,:,:)

      work_g2(:,:) = c0
      if (my_task == master_task) then
         do j = 1, ny_global
         do i = 2, nx_global
            work_g2(i,j) = work_g(i-1,j  ) * rad_to_deg
         enddo
         enddo
         ! extrapolate
         do j = 1, ny_global
            work_g2(1,j) = c2*work_g2(2,j) - work_g2(3,j)
         enddo
      endif
      call scatter_global(work1, work_g2, &
                          master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)
      vbounds(4,:,:,:) = work1(:,:,:)

      deallocate (work_g2, stat=ierr)
      if (ierr/=0) call abort_ice(subname//' ERROR: Dealloc error', file=__FILE__, line=__LINE__)

      end subroutine gridbox_verts

!=======================================================================
! ocean bathymetry for grounded sea ice (seabed stress) or icebergs
! currently hardwired for 40 levels (gx3, gx1 grids)
! should be read from a file instead (see subroutine read_seabedstress_bathy)

      subroutine get_bathymetry

      integer (kind=int_kind) :: &
         i, j, k, iblk      ! loop indices

      integer (kind=int_kind), parameter :: &
         nlevel = 40        ! number of layers (gx3 grid)

      real (kind=dbl_kind), dimension(nlevel) :: &
         depth              ! total depth, m

      logical (kind=log_kind) :: &
         calc_dragio

      real (kind=dbl_kind), dimension(nlevel), parameter :: &
         thick  = (/ &                        ! ocean layer thickness, m
            10.01244_dbl_kind,  10.11258_dbl_kind,  10.31682_dbl_kind, &
            10.63330_dbl_kind,  11.07512_dbl_kind,  11.66145_dbl_kind, &
            12.41928_dbl_kind,  13.38612_dbl_kind,  14.61401_dbl_kind, &
            16.17561_dbl_kind,  18.17368_dbl_kind,  20.75558_dbl_kind, &
            24.13680_dbl_kind,  28.63821_dbl_kind,  34.74644_dbl_kind, &
            43.20857_dbl_kind,  55.16812_dbl_kind,  72.30458_dbl_kind, &
            96.74901_dbl_kind,  130.0392_dbl_kind,  170.0489_dbl_kind, &
            207.9933_dbl_kind,  233.5694_dbl_kind,  245.2719_dbl_kind, &
            248.9804_dbl_kind,  249.8322_dbl_kind,  249.9787_dbl_kind, &
            249.9979_dbl_kind,  249.9998_dbl_kind,  250.0000_dbl_kind, &
            250.0000_dbl_kind,  250.0000_dbl_kind,  250.0000_dbl_kind, &
            250.0000_dbl_kind,  250.0000_dbl_kind,  250.0000_dbl_kind, &
            250.0000_dbl_kind,  250.0000_dbl_kind,  250.0000_dbl_kind, &
            250.0000_dbl_kind   /)

      character(len=*), parameter :: subname = '(get_bathymetry)'

      call icepack_query_parameters(calc_dragio_out=calc_dragio)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (use_bathymetry) then

         call read_seabedstress_bathy

      else

         ! convert to total depth
         depth(1) = thick(1)
         do k = 2, nlevel
            depth(k) = depth(k-1) + thick(k)
         enddo

         bathymetry = c0
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               k = min(nint(kmt(i,j,iblk)),nlevel)
               if (k > nlevel) call abort_ice(subname//' kmt gt nlevel error', &
                                    file=__FILE__, line=__LINE__)
               if (k > 0) bathymetry(i,j,iblk) = depth(k)
            enddo
            enddo
         enddo

         ! For consistency, set thickness_ocn_layer1 in Icepack if 'calc_dragio' is active
         if (calc_dragio) then
            call icepack_init_parameters(thickness_ocn_layer1_in=thick(1))
            call icepack_warnings_flush(nu_diag)
            if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
               file=__FILE__, line=__LINE__)
         endif

      endif ! bathymetry_file

      end subroutine get_bathymetry

!=======================================================================
! with use_bathymetry = false, vertical depth profile generated for max KMT
! with use_bathymetry = true, expects to read in pop vert_grid file

      subroutine get_bathymetry_popfile

      integer (kind=int_kind) :: &
         i, j, k, iblk      ! loop indices

      integer (kind=int_kind) :: &
         ntmp, nlevel   , & ! number of levels (max KMT)
         k1             , & ! levels
         ierr           , & ! error tag
         fid                ! fid unit number

      real (kind=dbl_kind), dimension(:),allocatable :: &
         depth          , & ! total depth, m
         thick              ! layer thickness, cm -> m

      logical (kind=log_kind) :: &
         calc_dragio

      character(len=*), parameter :: subname = '(get_bathymetry_popfile)'

      ntmp = maxval(nint(KMT))
      nlevel = global_maxval(ntmp,distrb_info)

      if (my_task==master_task) then
         write(nu_diag,*) subname,' KMT max = ',nlevel
      endif

      allocate(depth(nlevel),thick(nlevel), stat=ierr)
      if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory', file=__FILE__, line=__LINE__)
      thick = -999999.
      depth = -999999.

      if (use_bathymetry) then

         write (nu_diag,*) subname,' Bathymetry file = ', trim(bathymetry_file)
         if (my_task == master_task) then
            call get_fileunit(fid)
            open(fid,file=bathymetry_file,form='formatted',iostat=ierr)
            if (ierr/=0) call abort_ice(subname//' open error', file=__FILE__, line=__LINE__)
            do k = 1,nlevel
               read(fid,*,iostat=ierr) thick(k)
               if (ierr/=0) call abort_ice(subname//' read error', file=__FILE__, line=__LINE__)
            enddo
            call release_fileunit(fid)
         endif

         call broadcast_array(thick,master_task)

      else

         ! create thickness profile
         k1 = min(5,nlevel)
         do k = 1,k1
            thick(k) = max(10000._dbl_kind/float(nlevel),500._dbl_kind)
         enddo
         do k = k1+1,nlevel
            thick(k) = min(thick(k-1)*1.2_dbl_kind,20000._dbl_kind)
         enddo

      endif

      ! convert thick from cm to m
      thick = thick / 100._dbl_kind

      ! convert to total depth
      depth(1) = thick(1)
      do k = 2, nlevel
         depth(k) = depth(k-1) + thick(k)
         if (depth(k) < 0.) call abort_ice(subname//' negative depth error', file=__FILE__, line=__LINE__)
      enddo

      if (my_task==master_task) then
         do k = 1,nlevel
           write(nu_diag,'(2a,i6,2f13.7)') subname,'   k, thick(m), depth(m) = ',k,thick(k),depth(k)
         enddo
      endif

      bathymetry = 0._dbl_kind
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            k = nint(kmt(i,j,iblk))
            if (k > nlevel) call abort_ice(subname//' kmt gt nlevel error', file=__FILE__, line=__LINE__)
            if (k > 0) bathymetry(i,j,iblk) = depth(k)
         enddo
         enddo
      enddo

      ! For consistency, set thickness_ocn_layer1 in Icepack if 'calc_dragio' is active
      call icepack_query_parameters(calc_dragio_out=calc_dragio)
      if (calc_dragio) then
         call icepack_init_parameters(thickness_ocn_layer1_in=thick(1))
      endif
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      deallocate(depth,thick, stat=ierr)
      if (ierr/=0) call abort_ice(subname//' ERROR: Dealloc error', file=__FILE__, line=__LINE__)

      end subroutine get_bathymetry_popfile

!=======================================================================
! Read bathymetry data for seabed stress calculation (grounding scheme for
! landfast ice) in CICE stand-alone mode. When CICE is in coupled mode
! (e.g. CICE-NEMO), hwater should be uptated at each time level so that
! it varies with ocean dynamics.
!
! author: Fred Dupont, CMC

      subroutine read_seabedstress_bathy

      ! use module
      use ice_read_write

      ! local variables
      integer (kind=int_kind) :: &
         fid_init        ! file id for netCDF init file

      character (char_len_long) :: &        ! input data file names
         fieldname

      logical (kind=log_kind) :: diag=.true.

      character(len=*), parameter :: subname = '(read_seabedstress_bathy)'

      if (my_task == master_task) then
          write (nu_diag,*) ' '
          write (nu_diag,*) 'Bathymetry file: ', trim(bathymetry_file)
          call icepack_warnings_flush(nu_diag)
      endif

      call ice_open_nc(bathymetry_file,fid_init)

      fieldname='Bathymetry'

      if (my_task == master_task) then
         write(nu_diag,*) subname,' reading ',TRIM(fieldname)
         call icepack_warnings_flush(nu_diag)
      endif
      call ice_read_nc(fid_init,1,fieldname,bathymetry,diag, &
                    field_loc=field_loc_center, &
                    field_type=field_type_scalar)

      call ice_close_nc(fid_init)

      if (my_task == master_task) then
         write(nu_diag,*) subname,' closing file ',TRIM(bathymetry_file)
         call icepack_warnings_flush(nu_diag)
      endif

      end subroutine read_seabedstress_bathy

!=======================================================================

      end module ice_grid

!=======================================================================
