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

      module ice_grid

      use ice_kinds_mod
      use ice_broadcast, only: broadcast_scalar, broadcast_array
      use ice_boundary, only: ice_HaloUpdate, ice_HaloExtrapolate
      use ice_communicate, only: my_task, master_task
      use ice_blocks, only: block, get_block, nx_block, ny_block, nghost
      use ice_domain_size, only: nx_global, ny_global, max_blocks
      use ice_domain, only: blocks_ice, nblocks, halo_info, distrb_info, &
          ew_boundary_type, ns_boundary_type, init_domain_distribution
      use ice_fileunits, only: nu_diag, nu_grid, nu_kmt, &
          get_fileunit, release_fileunit, flush_fileunit
      use ice_gather_scatter, only: gather_global, scatter_global
      use ice_read_write, only: ice_read, ice_read_nc, ice_read_global, &
          ice_read_global_nc, ice_open, ice_open_nc, ice_close_nc
      use ice_timers, only: timer_bound, ice_timer_start, ice_timer_stop
      use ice_exit, only: abort_ice
      use ice_global_reductions, only: global_minval, global_maxval
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_parameters

      implicit none
      private
      public :: init_grid1, init_grid2, &
                t2ugrid_vector, u2tgrid_vector, &
                to_ugrid, to_tgrid, alloc_grid

      character (len=char_len_long), public :: &
         grid_format  , & ! file format ('bin'=binary or 'nc'=netcdf)
         gridcpl_file , & !  input file for POP coupling grid info
         grid_file    , & !  input file for POP grid info
         kmt_file     , & !  input file for POP grid info
         bathymetry_file, & !  input bathymetry for seabed stress
         bathymetry_format, & ! bathymetry file format (default or pop)
         grid_spacing , & !  default of 30.e3m or set by user in namelist 
         grid_type        !  current options are rectangular (default),
                          !  displaced_pole, tripole, regional

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         dxt    , & ! width of T-cell through the middle (m)
         dyt    , & ! height of T-cell through the middle (m)
         dxu    , & ! width of U-cell through the middle (m)
         dyu    , & ! height of U-cell through the middle (m)
         HTE    , & ! length of eastern edge of T-cell (m)
         HTN    , & ! length of northern edge of T-cell (m)
         tarea  , & ! area of T-cell (m^2)
         uarea  , & ! area of U-cell (m^2)
         tarear , & ! 1/tarea
         uarear , & ! 1/uarea
         tinyarea,& ! puny*tarea
         tarean , & ! area of NH T-cells
         tareas , & ! area of SH T-cells
         ULON   , & ! longitude of velocity pts (radians)
         ULAT   , & ! latitude of velocity pts (radians)
         TLON   , & ! longitude of temp pts (radians)
         TLAT   , & ! latitude of temp pts (radians)
         ANGLE  , & ! for conversions between POP grid and lat/lon
         ANGLET , & ! ANGLE converted to T-cells
         bathymetry      , & ! ocean depth, for grounding keels and bergs (m)
         ocn_gridcell_frac   ! only relevant for lat-lon grids
                             ! gridcell value of [1 - (land fraction)] (T-cell)

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         cyp    , & ! 1.5*HTE - 0.5*HTE
         cxp    , & ! 1.5*HTN - 0.5*HTN
         cym    , & ! 0.5*HTE - 1.5*HTE
         cxm    , & ! 0.5*HTN - 1.5*HTN
         dxhy   , & ! 0.5*(HTE - HTE)
         dyhx       ! 0.5*(HTN - HTN)

      ! grid dimensions for rectangular grid
      real (kind=dbl_kind), public ::  &
         dxrect, & !  user_specified spacing (cm) in x-direction (uniform HTN)
         dyrect    !  user_specified spacing (cm) in y-direction (uniform HTE)

      ! Corners of grid boxes for history output
      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         lont_bounds, & ! longitude of gridbox corners for T point
         latt_bounds, & ! latitude of gridbox corners for T point
         lonu_bounds, & ! longitude of gridbox corners for U point
         latu_bounds    ! latitude of gridbox corners for U point       

      ! geometric quantities used for remapping transport
      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         xav  , & ! mean T-cell value of x
         yav  , & ! mean T-cell value of y
         xxav , & ! mean T-cell value of xx
!         xyav , & ! mean T-cell value of xy
!         yyav , & ! mean T-cell value of yy
         yyav     ! mean T-cell value of yy
!         xxxav, & ! mean T-cell value of xxx
!         xxyav, & ! mean T-cell value of xxy
!         xyyav, & ! mean T-cell value of xyy
!         yyyav    ! mean T-cell value of yyy

      real (kind=dbl_kind), &
         dimension (:,:,:,:,:), allocatable, public :: &
         mne, & ! matrices used for coordinate transformations in remapping
         mnw, & ! ne = northeast corner, nw = northwest, etc.
         mse, & 
         msw

      ! masks
      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         hm     , & ! land/boundary mask, thickness (T-cell)
         bm     , & ! task/block id
         uvm    , & ! land/boundary mask, velocity (U-cell)
         kmt        ! ocean topography mask for bathymetry (T-cell)

      logical (kind=log_kind), public :: &
         use_bathymetry     ! flag for reading in bathymetry_file

      logical (kind=log_kind), &
         dimension (:,:,:), allocatable, public :: &
         tmask  , & ! land/boundary mask, thickness (T-cell)
         umask  , & ! land/boundary mask, velocity (U-cell)
         lmask_n, & ! northern hemisphere mask
         lmask_s    ! southern hemisphere mask

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         rndex_global       ! global index for local subdomain (dbl)

      logical (kind=log_kind), private :: &
         l_readCenter ! If anglet exist in grid file read it otherwise calculate it


!=======================================================================

      contains

!=======================================================================
!
! Allocate space for all variables 
!
      subroutine alloc_grid

      integer (int_kind) :: ierr

      allocate( &
         dxt      (nx_block,ny_block,max_blocks), & ! width of T-cell through the middle (m)
         dyt      (nx_block,ny_block,max_blocks), & ! height of T-cell through the middle (m)
         dxu      (nx_block,ny_block,max_blocks), & ! width of U-cell through the middle (m)
         dyu      (nx_block,ny_block,max_blocks), & ! height of U-cell through the middle (m)
         HTE      (nx_block,ny_block,max_blocks), & ! length of eastern edge of T-cell (m)
         HTN      (nx_block,ny_block,max_blocks), & ! length of northern edge of T-cell (m)
         tarea    (nx_block,ny_block,max_blocks), & ! area of T-cell (m^2)
         uarea    (nx_block,ny_block,max_blocks), & ! area of U-cell (m^2)
         tarear   (nx_block,ny_block,max_blocks), & ! 1/tarea
         uarear   (nx_block,ny_block,max_blocks), & ! 1/uarea
         tinyarea (nx_block,ny_block,max_blocks), & ! puny*tarea
         tarean   (nx_block,ny_block,max_blocks), & ! area of NH T-cells
         tareas   (nx_block,ny_block,max_blocks), & ! area of SH T-cells
         ULON     (nx_block,ny_block,max_blocks), & ! longitude of velocity pts (radians)
         ULAT     (nx_block,ny_block,max_blocks), & ! latitude of velocity pts (radians)
         TLON     (nx_block,ny_block,max_blocks), & ! longitude of temp pts (radians)
         TLAT     (nx_block,ny_block,max_blocks), & ! latitude of temp pts (radians)
         ANGLE    (nx_block,ny_block,max_blocks), & ! for conversions between POP grid and lat/lon
         ANGLET   (nx_block,ny_block,max_blocks), & ! ANGLE converted to T-cells
         bathymetry(nx_block,ny_block,max_blocks),& ! ocean depth, for grounding keels and bergs (m)
         ocn_gridcell_frac(nx_block,ny_block,max_blocks),& ! only relevant for lat-lon grids
         cyp      (nx_block,ny_block,max_blocks), & ! 1.5*HTE - 0.5*HTE
         cxp      (nx_block,ny_block,max_blocks), & ! 1.5*HTN - 0.5*HTN
         cym      (nx_block,ny_block,max_blocks), & ! 0.5*HTE - 1.5*HTE
         cxm      (nx_block,ny_block,max_blocks), & ! 0.5*HTN - 1.5*HTN
         dxhy     (nx_block,ny_block,max_blocks), & ! 0.5*(HTE - HTE)
         dyhx     (nx_block,ny_block,max_blocks), & ! 0.5*(HTN - HTN)
         xav      (nx_block,ny_block,max_blocks), & ! mean T-cell value of x
         yav      (nx_block,ny_block,max_blocks), & ! mean T-cell value of y
         xxav     (nx_block,ny_block,max_blocks), & ! mean T-cell value of xx
         yyav     (nx_block,ny_block,max_blocks), & ! mean T-cell value of yy
         hm       (nx_block,ny_block,max_blocks), & ! land/boundary mask, thickness (T-cell)
         bm       (nx_block,ny_block,max_blocks), & ! task/block id
         uvm      (nx_block,ny_block,max_blocks), & ! land/boundary mask, velocity (U-cell)
         kmt      (nx_block,ny_block,max_blocks), & ! ocean topography mask for bathymetry (T-cell)
         tmask    (nx_block,ny_block,max_blocks), & ! land/boundary mask, thickness (T-cell)
         umask    (nx_block,ny_block,max_blocks), & ! land/boundary mask, velocity (U-cell)
         lmask_n  (nx_block,ny_block,max_blocks), & ! northern hemisphere mask
         lmask_s  (nx_block,ny_block,max_blocks), & ! southern hemisphere mask
         rndex_global(nx_block,ny_block,max_blocks), & ! global index for local subdomain (dbl)
         lont_bounds(4,nx_block,ny_block,max_blocks), & ! longitude of gridbox corners for T point
         latt_bounds(4,nx_block,ny_block,max_blocks), & ! latitude of gridbox corners for T point
         lonu_bounds(4,nx_block,ny_block,max_blocks), & ! longitude of gridbox corners for U point
         latu_bounds(4,nx_block,ny_block,max_blocks), & ! latitude of gridbox corners for U point       
         mne  (2,2,nx_block,ny_block,max_blocks), & ! matrices used for coordinate transformations in remapping
         mnw  (2,2,nx_block,ny_block,max_blocks), & ! ne = northeast corner, nw = northwest, etc.
         mse  (2,2,nx_block,ny_block,max_blocks), &
         msw  (2,2,nx_block,ny_block,max_blocks), &
         stat=ierr)
      if (ierr/=0) call abort_ice('(alloc_grid): Out of memory')

      end subroutine alloc_grid

!=======================================================================

! Distribute blocks across processors.  The distribution is optimized
! based on latitude and topography, contained in the ULAT and KMT arrays. 
!
! authors: William Lipscomb and Phil Jones, LANL

      subroutine init_grid1 

      use ice_blocks, only: nx_block, ny_block
      use ice_broadcast, only: broadcast_array
      use ice_constants, only: c1

      integer (kind=int_kind) :: &
         fid_grid, &     ! file id for netCDF grid file
         fid_kmt         ! file id for netCDF kmt file

      character (char_len) :: &
         fieldname       ! field name in netCDF file

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g1, work_g2

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

      allocate(work_g1(nx_global,ny_global))
      allocate(work_g2(nx_global,ny_global))

      if (trim(grid_type) == 'displaced_pole' .or. &
          trim(grid_type) == 'tripole' .or. &
          trim(grid_type) == 'regional'     ) then

         if (trim(grid_format) == 'nc') then

            call ice_open_nc(grid_file,fid_grid)
            call ice_open_nc(kmt_file,fid_kmt)

            fieldname='ulat'
            call ice_read_global_nc(fid_grid,1,fieldname,work_g1,.true.)
            fieldname='kmt'
            call ice_read_global_nc(fid_kmt,1,fieldname,work_g2,.true.)

            if (my_task == master_task) then
               call ice_close_nc(fid_grid)
               call ice_close_nc(fid_kmt)
            endif

         else

            call ice_open(nu_grid,grid_file,64) ! ULAT
            call ice_open(nu_kmt, kmt_file, 32) ! KMT

            call ice_read_global(nu_grid,1,work_g1,'rda8',.true.)  ! ULAT
            call ice_read_global(nu_kmt, 1,work_g2,'ida4',.true.)  ! KMT

            if (my_task == master_task) then
               close (nu_grid)
               close (nu_kmt)
            endif

         endif

      else   ! rectangular grid

         work_g1(:,:) = 75._dbl_kind/rad_to_deg  ! arbitrary polar latitude
         work_g2(:,:) = c1

      endif

      call broadcast_array(work_g1, master_task)   ! ULAT
      call broadcast_array(work_g2, master_task)   ! KMT

      !-----------------------------------------------------------------
      ! distribute blocks among processors
      !-----------------------------------------------------------------

      call init_domain_distribution(work_g2, work_g1)  ! KMT, ULAT

      deallocate(work_g1)
      deallocate(work_g2)

      !-----------------------------------------------------------------
      ! write additional domain information
      !-----------------------------------------------------------------

      if (my_task == master_task) then
        write(nu_diag,'(a26,i6)') '  Block size:  nx_block = ',nx_block
        write(nu_diag,'(a26,i6)') '               ny_block = ',ny_block
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

      use ice_blocks, only: get_block, block, nx_block, ny_block
      use ice_constants, only: c0, c1, c2, p5, p25, c1p5, &
          field_loc_center, field_loc_NEcorner, &
          field_type_scalar, field_type_vector, field_type_angle
      use ice_domain_size, only: max_blocks

      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      real (kind=dbl_kind) :: &
         angle_0, angle_w, angle_s, angle_sw, &
         pi, pi2, puny

      logical (kind=log_kind), dimension(nx_block,ny_block,max_blocks):: &
         out_of_range

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g1

      type (block) :: &
         this_block           ! block information for current block

      character(len=*), parameter :: subname = '(init_grid2)'

      !-----------------------------------------------------------------
      ! lat, lon, cell widths, angle, land mask
      !-----------------------------------------------------------------

      call icepack_query_parameters(pi_out=pi, pi2_out=pi2, puny_out=puny)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (trim(grid_type) == 'displaced_pole' .or. &
          trim(grid_type) == 'tripole' .or. &
          trim(grid_type) == 'regional'      ) then
         if (trim(grid_format) == 'nc') then
            call popgrid_nc     ! read POP grid lengths from nc file
         else
            call popgrid        ! read POP grid lengths directly
         endif 
#ifdef CESMCOUPLED
      elseif (trim(grid_type) == 'latlon') then
         call latlongrid        ! lat lon grid for sequential CESM (CAM mode)
         return
#endif
      elseif (trim(grid_type) == 'cpom_grid') then
         call cpomgrid          ! cpom model orca1 type grid
      else
         call rectgrid          ! regular rectangular grid
      endif

      !-----------------------------------------------------------------
      ! T-grid cell and U-grid cell quantities
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            tarea(i,j,iblk) = dxt(i,j,iblk)*dyt(i,j,iblk)
            uarea(i,j,iblk) = dxu(i,j,iblk)*dyu(i,j,iblk)
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
            tinyarea(i,j,iblk) = puny*tarea(i,j,iblk)

            dxhy(i,j,iblk) = p5*(HTE(i,j,iblk) - HTE(i-1,j,iblk))
            dyhx(i,j,iblk) = p5*(HTN(i,j,iblk) - HTN(i,j-1,iblk))
         enddo
         enddo

         do j = jlo, jhi+1
         do i = ilo, ihi+1
            cyp(i,j,iblk) = (c1p5*HTE(i,j,iblk) - p5*HTE(i-1,j,iblk))
            cxp(i,j,iblk) = (c1p5*HTN(i,j,iblk) - p5*HTN(i,j-1,iblk))
            ! match order of operations in cyp, cxp for tripole grids
            cym(i,j,iblk) = -(c1p5*HTE(i-1,j,iblk) - p5*HTE(i,j,iblk)) 
            cxm(i,j,iblk) = -(c1p5*HTN(i,j-1,iblk) - p5*HTN(i,j,iblk)) 
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
      call ice_HaloUpdate (tarea,              halo_info, &
                           field_loc_center,   field_type_scalar, &
                           fillValue=c1)
      call ice_HaloUpdate (uarea,              halo_info, &
                           field_loc_NEcorner, field_type_scalar, &
                           fillValue=c1)
      call ice_HaloUpdate (tarear,             halo_info, &
                           field_loc_center,   field_type_scalar, &
                           fillValue=c1)
      call ice_HaloUpdate (uarear,             halo_info, &
                           field_loc_NEcorner, field_type_scalar, &
                           fillValue=c1)
      call ice_HaloUpdate (tinyarea,           halo_info, &
                           field_loc_center,   field_type_scalar, &
                           fillValue=c1)
      call ice_HaloUpdate (dxhy,               halo_info, &
                           field_loc_center,   field_type_vector, &
                           fillValue=c1)
      call ice_HaloUpdate (dyhx,               halo_info, &
                           field_loc_center,   field_type_vector, &
                           fillValue=c1)
      call ice_timer_stop(timer_bound)

      !-----------------------------------------------------------------
      ! Calculate ANGLET to be compatible with POP ocean model
      ! First, ensure that -pi <= ANGLE <= pi
      !-----------------------------------------------------------------

      out_of_range = .false.
      where (ANGLE < -pi .or. ANGLE > pi) out_of_range = .true.
      if (count(out_of_range) > 0) then
         write(nu_diag,*) subname,' angle = ',minval(ANGLE),maxval(ANGLE),count(out_of_range)
         call abort_ice (subname//' ANGLE out of expected range', &
             file=__FILE__, line=__LINE__)
      endif

      !-----------------------------------------------------------------
      ! Compute ANGLE on T-grid
      !-----------------------------------------------------------------
      if (trim(grid_type) == 'cpom_grid') then
         ANGLET(:,:,:) = ANGLE(:,:,:)
      else if (.not. (l_readCenter)) then
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
      endif ! cpom_grid
      if (trim(grid_type) == 'regional' .and. &
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
      if (.not. (l_readCenter)) then
         call Tlatlon           ! get lat, lon on the T grid
      endif
      !-----------------------------------------------------------------
      ! bathymetry
      !-----------------------------------------------------------------

      if (trim(bathymetry_format) == 'default') then
         call get_bathymetry
      elseif (trim(bathymetry_format) == 'pop') then
         call get_bathymetry_popfile
      else
         call abort_ice(subname//'ERROR: bathymetry_format value must be default or pop', &
            file=__FILE__, line=__LINE__)
      endif

      !----------------------------------------------------------------
      ! Corner coordinates for CF compliant history files
      !----------------------------------------------------------------

      call gridbox_corners

      !-----------------------------------------------------------------
      ! Compute global index (used for unpacking messages from coupler)
      !-----------------------------------------------------------------

      if (my_task==master_task) then
         allocate(work_g1(nx_global,ny_global))
         do j=1,ny_global
         do i=1,nx_global
            work_g1(i,j) = real((j-1)*nx_global + i,kind=dbl_kind)
         enddo
         enddo
      else
         allocate(work_g1(1,1)) ! to save memory
      endif

      call scatter_global(rndex_global, work_g1,  &
                          master_task,  distrb_info, &
                          field_loc_center, field_type_scalar)

      deallocate(work_g1)

      end subroutine init_grid2

!=======================================================================

! POP displaced pole grid and land mask (or tripole). 
! Grid record number, field and units are: \\
! (1) ULAT  (radians)    \\
! (2) ULON  (radians)    \\
! (3) HTN   (cm)         \\
! (4) HTE   (cm)         \\
! (5) HUS   (cm)         \\
! (6) HUW   (cm)         \\
! (7) ANGLE (radians)   
!
! Land mask record number and field is (1) KMT.
!
! author: Elizabeth C. Hunke, LANL

      subroutine popgrid

      use ice_blocks, only: nx_block, ny_block
      use ice_constants, only: c0, c1, &
          field_loc_center, field_loc_NEcorner, &
          field_type_scalar, field_type_angle
      use ice_domain_size, only: max_blocks

      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      logical (kind=log_kind) :: diag

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g1

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1

      type (block) :: &
         this_block           ! block information for current block

      character(len=*), parameter :: subname = '(popgrid)'

      call ice_open(nu_grid,grid_file,64)
      call ice_open(nu_kmt,kmt_file,32)

      diag = .true.       ! write diagnostic info

      !-----------------------------------------------------------------
      ! topography
      !-----------------------------------------------------------------

      call ice_read(nu_kmt,1,work1,'ida4',diag, &
                    field_loc=field_loc_center, & 
                    field_type=field_type_scalar)

      hm (:,:,:) = c0
      kmt(:,:,:) = c0
      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            kmt(i,j,iblk) = work1(i,j,iblk)
            if (kmt(i,j,iblk) >= c1) hm(i,j,iblk) = c1
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !-----------------------------------------------------------------
      ! lat, lon, angle
      !-----------------------------------------------------------------

      allocate(work_g1(nx_global,ny_global))

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
      call primary_grid_lengths_HTN(work_g1)                  ! dxu, dxt

      call ice_read_global(nu_grid,4,work_g1,'rda8',.true.)   ! HTE
      call primary_grid_lengths_HTE(work_g1)                  ! dyu, dyt

      deallocate(work_g1)

      if (my_task == master_task) then
         close (nu_grid)
         close (nu_kmt)
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
! Land mask record number and field is (1) KMT.
!
! author: Elizabeth C. Hunke, LANL
! Revised for netcdf input: Ann Keen, Met Office, May 2007

      subroutine popgrid_nc

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
         ilo,ihi,jlo,jhi, &     ! beginning and end of physical domain
         fid_grid, &            ! file id for netCDF grid file
         fid_kmt                ! file id for netCDF kmt file

      logical (kind=log_kind) :: diag

      character (char_len) :: &
         fieldname              ! field name in netCDF file

      real (kind=dbl_kind) :: &
         pi

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g1

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1

      type (block) :: &
         this_block           ! block information for current block
      
      integer(kind=int_kind) :: &
         varid
      integer (kind=int_kind) :: &
         status                ! status flag


      character(len=*), parameter :: subname = '(popgrid_nc)'

#ifdef USE_NETCDF
      call icepack_query_parameters(pi_out=pi)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      call ice_open_nc(grid_file,fid_grid)
      call ice_open_nc(kmt_file,fid_kmt)

      diag = .true.       ! write diagnostic info
      l_readCenter = .false.
      !-----------------------------------------------------------------
      ! topography
      !-----------------------------------------------------------------

      fieldname='kmt'
      call ice_read_nc(fid_kmt,1,fieldname,work1,diag, &
                       field_loc=field_loc_center, & 
                       field_type=field_type_scalar)

      hm (:,:,:) = c0
      kmt(:,:,:) = c0
      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            kmt(i,j,iblk) = work1(i,j,iblk)
            if (kmt(i,j,iblk) >= c1) hm(i,j,iblk) = c1
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

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
         call ice_close_nc(fid_kmt)
      endif
#else
      call abort_ice(subname//'ERROR: USE_NETCDF cpp not defined', &
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

!     use ice_boundary
      use ice_domain_size
      use ice_scam, only : scmlat, scmlon, single_column
      use ice_constants, only: c0, c1, p5, p25, &
          field_loc_center, field_type_scalar, radius
#ifdef USE_NETCDF
      use netcdf
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
         puny, &
         scamdata             ! temporary

      character(len=*), parameter :: subname = '(lonlatgrid)'

#ifdef USE_NETCDF
      !-----------------------------------------------------------------
      ! - kmt file is actually clm fractional land file
      ! - Determine consistency of dimensions
      ! - Read in lon/lat centers in degrees from kmt file
      ! - Read in ocean from "kmt" file (1 for ocean, 0 for land)
      !-----------------------------------------------------------------

      call icepack_query_parameters(pi_out=pi, puny_out=puny)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      ! Determine dimension of domain file and check for consistency

      if (my_task == master_task) then
         call ice_open_nc(kmt_file, ncid)

         status = nf90_inq_dimid (ncid, 'ni', dimid)
         status = nf90_inquire_dimension(ncid, dimid, len=ni)
         status = nf90_inq_dimid (ncid, 'nj', dimid)
         status = nf90_inquire_dimension(ncid, dimid, len=nj)
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
               call abort_ice (subname//'ERROR: check nx_global, ny_global')
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
         if (status /= nf90_noerr) call abort_ice (subname//' inq_varid xc')
         status = nf90_get_var(ncid, varid, glob_grid, start3, count3)
         if (status /= nf90_noerr) call abort_ice (subname//' get_var xc')
         do i = 1,ni
            lons(i) = glob_grid(i,1)
         end do

         status = nf90_inq_varid(ncid, 'yc' , varid)
         if (status /= nf90_noerr) call abort_ice (subname//' inq_varid yc')
         status = nf90_get_var(ncid, varid, glob_grid, start3, count3)
         if (status /= nf90_noerr) call abort_ice (subname//' get_var yc')
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
         if (status /= nf90_noerr) call abort_ice (subname//' inq_varid xc')
         status = nf90_get_var(ncid, varid, scamdata, start)
         if (status /= nf90_noerr) call abort_ice (subname//' get_var xc')
         TLON = scamdata
         status = nf90_inq_varid(ncid, 'yc' , varid)
         if (status /= nf90_noerr) call abort_ice (subname//' inq_varid yc')
         status = nf90_get_var(ncid, varid, scamdata, start)
         if (status /= nf90_noerr) call abort_ice (subname//' get_var yc')
         TLAT = scamdata
         status = nf90_inq_varid(ncid, 'area' , varid)
         if (status /= nf90_noerr) call abort_ice (subname//' inq_varid area')
         status = nf90_get_var(ncid, varid, scamdata, start)
         if (status /= nf90_noerr) call abort_ice (subname//' get_var are')
         tarea = scamdata
         status = nf90_inq_varid(ncid, 'mask' , varid)
         if (status /= nf90_noerr) call abort_ice (subname//' inq_varid mask')
         status = nf90_get_var(ncid, varid, scamdata, start)
         if (status /= nf90_noerr) call abort_ice (subname//' get_var mask')
         hm = scamdata
         status = nf90_inq_varid(ncid, 'frac' , varid)
         if (status /= nf90_noerr) call abort_ice (subname//' inq_varid frac')
         status = nf90_get_var(ncid, varid, scamdata, start)
         if (status /= nf90_noerr) call abort_ice (subname//' get_var frac')
         ocn_gridcell_frac = scamdata
      else
         ! Check for consistency
         if (my_task == master_task) then
            if (nx_global /= ni .and. ny_global /= nj) then
              write(nu_diag,*) 'latlongrid: ni,nj = ',ni,nj
              write(nu_diag,*) 'latlongrid: nx_g,ny_g = ',nx_global, ny_global
              call abort_ice (subname//'ERROR: ni,nj not equal to nx_global,ny_global')
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
            tinyarea(i,j,iblk) = puny*tarea(i,j,iblk)

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
            ANGLE (i,j,iblk) = c0                             

            ANGLET(i,j,iblk) = c0                             
            HTN   (i,j,iblk) = 1.e36_dbl_kind
            HTE   (i,j,iblk) = 1.e36_dbl_kind
            dxt   (i,j,iblk) = 1.e36_dbl_kind
            dyt   (i,j,iblk) = 1.e36_dbl_kind
            dxu   (i,j,iblk) = 1.e36_dbl_kind
            dyu   (i,j,iblk) = 1.e36_dbl_kind
            dxhy  (i,j,iblk) = 1.e36_dbl_kind
            dyhx  (i,j,iblk) = 1.e36_dbl_kind
            cyp   (i,j,iblk) = 1.e36_dbl_kind
            cxp   (i,j,iblk) = 1.e36_dbl_kind
            cym   (i,j,iblk) = 1.e36_dbl_kind
            cxm   (i,j,iblk) = 1.e36_dbl_kind
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      call makemask
#else
      call abort_ice(subname//'ERROR: USE_NETCDF cpp not defined', &
          file=__FILE__, line=__LINE__)
#endif

      end subroutine latlongrid
#endif

!=======================================================================

! Regular rectangular grid and mask
!
! author: Elizabeth C. Hunke, LANL

      subroutine rectgrid

      use ice_blocks, only: nx_block, ny_block
      use ice_constants, only: c0, c1, c2, radius, cm_to_m, &
          field_loc_center, field_loc_NEcorner, field_type_scalar
      use ice_domain, only: close_boundaries

      integer (kind=int_kind) :: &
         i, j, iblk, &
         imid, jmid

      real (kind=dbl_kind) :: &
         length, &
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

      ! Weddell Sea
      ! lower left corner of grid is 55W, 75S

      ! Barrow AK
      ! lower left corner of grid is 156.5W, 71.35N

      if (my_task == master_task) then
         work_g1 = c0
         length = dxrect*cm_to_m/radius*rad_to_deg

!         work_g1(1,:) = -55._dbl_kind   ! Weddell Sea
         work_g1(1,:) = -156.5_dbl_kind ! Barrow AK

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

!         work_g1(:,1) = -75._dbl_kind ! Weddell Sea
         work_g1(:,1) = 71.35_dbl_kind ! Barrow AK

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
      call primary_grid_lengths_HTN(work_g1)  ! dxu, dxt

      if (my_task == master_task) then
         do j = 1, ny_global
         do i = 1, nx_global
            work_g1(i,j) = dyrect             ! HTE
         enddo
         enddo
      endif
      call primary_grid_lengths_HTE(work_g1)  ! dyu, dyt

      !-----------------------------------------------------------------
      ! Construct T-cell land mask
      ! Keyed on ew_boundary_type; ns_boundary_type should be 'open'.
      !-----------------------------------------------------------------

      if (my_task == master_task) then
         work_g1(:,:) = c0      ! initialize hm as land

         if (trim(ew_boundary_type) == 'cyclic') then

            do j = 3,ny_global-2      ! closed top and bottom
            do i = 1,nx_global        ! open sides
               work_g1(i,j) = c1    ! NOTE nx_global > 5
            enddo
            enddo

         elseif (trim(ew_boundary_type) == 'open') then

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

            endif

            if (close_boundaries) then
              work_g1(:, 1:2) = c0
              work_g1(:, ny_global-1:ny_global) = c0
              work_g1(1:2, :) = c0
              work_g1(nx_global-1:nx_global, :) = c0
            endif

         elseif (trim(ew_boundary_type) == 'closed') then

            call abort_ice(subname//'ERROR: closed boundaries not available')

         endif
      endif

      call scatter_global(hm, work_g1, master_task, distrb_info, &
                          field_loc_center, field_type_scalar)

      deallocate(work_g1)

      end subroutine rectgrid

!=======================================================================

! CPOM displaced pole grid and land mask. \\
! Grid record number, field and units are: \\
! (1) ULAT  (degrees)    \\
! (2) ULON  (degrees)    \\
! (3) HTN   (m)          \\
! (4) HTE   (m)          \\
! (7) ANGLE (radians)    \\
!
! Land mask record number and field is (1) KMT.
!
! author: Adrian K. Turner, CPOM, UCL, 09/08/06

      subroutine cpomgrid

      use ice_blocks, only: nx_block, ny_block
      use ice_constants, only: c0, c1, m_to_cm, &
          field_loc_NEcorner, field_type_scalar
      use ice_domain_size, only: max_blocks

      integer (kind=int_kind) :: &
           i, j, iblk,           &
           ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      logical (kind=log_kind) :: diag

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g1

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1

      real (kind=dbl_kind) :: &
         rad_to_deg

      type (block) :: &
           this_block           ! block information for current block

      character(len=*), parameter :: subname = '(cpomgrid)'

      call icepack_query_parameters(rad_to_deg_out=rad_to_deg)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      call ice_open(nu_grid,grid_file,64)
      call ice_open(nu_kmt,kmt_file,32)

      diag = .true.       ! write diagnostic info

      ! topography
      call ice_read(nu_kmt,1,work1,'ida4',diag)

      hm (:,:,:) = c0
      kmt(:,:,:) = c0
      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi
            kmt(i,j,iblk) = work1(i,j,iblk)
            if (kmt(i,j,iblk) >= c1) hm(i,j,iblk) = c1
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      allocate(work_g1(nx_global,ny_global))

      ! lat, lon, cell dimensions, angles
      call ice_read_global(nu_grid,1,work_g1, 'rda8',diag)
      call scatter_global(ULAT, work_g1, master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)

      call ice_read_global(nu_grid,2,work_g1, 'rda8',diag)
      call scatter_global(ULON, work_g1, master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)

      call ice_read_global(nu_grid,3,work_g1,  'rda8',diag)
      work_g1 = work_g1 * m_to_cm
      call primary_grid_lengths_HTN(work_g1)  ! dxu, dxt

      call ice_read_global(nu_grid,4,work_g1,  'rda8',diag)
      work_g1 = work_g1 * m_to_cm
      call primary_grid_lengths_HTE(work_g1)  ! dyu, dyt

      call ice_read_global(nu_grid,7,work_g1,'rda8',diag)
      call scatter_global(ANGLE, work_g1, master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)

      ! fix units
      ULAT  = ULAT  / rad_to_deg
      ULON  = ULON  / rad_to_deg

      deallocate(work_g1)

      if (my_task == master_task) then
         close (nu_grid)
         close (nu_kmt)
      endif

      write(nu_diag,*) "min/max HTN: ", minval(HTN), maxval(HTN)
      write(nu_diag,*) "min/max HTE: ", minval(HTE), maxval(HTE)

      end subroutine cpomgrid

!=======================================================================

! Calculate dxu and dxt from HTN on the global grid, to preserve
! ghost cell and/or land values that might otherwise be lost. Scatter
! dxu, dxt and HTN to all processors.
!
! author: Elizabeth C. Hunke, LANL

      subroutine primary_grid_lengths_HTN(work_g)

      use ice_constants, only: p5, c2, cm_to_m, &
          field_loc_center, field_loc_NEcorner, &
          field_loc_Nface, field_type_scalar

      real (kind=dbl_kind), dimension(:,:) :: work_g ! global array holding HTN

      ! local variables

      integer (kind=int_kind) :: &
         i, j, &
         ip1     ! i+1

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g2

      character(len=*), parameter :: subname = '(primary_grid_lengths_HTN)'

      if (my_task == master_task) then
         allocate(work_g2(nx_global,ny_global))
      else
         allocate(work_g2(1,1))
      endif

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
         work_g2(i,j) = p5*(work_g(i,j) + work_g(ip1,j))    ! dxu
      enddo
      enddo
      endif
      call scatter_global(HTN, work_g, master_task, distrb_info, &
                          field_loc_Nface, field_type_scalar)
      call scatter_global(dxu, work_g2, master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)

      if (my_task == master_task) then
      do j = 2, ny_global
         do i = 1, nx_global
            work_g2(i,j) = p5*(work_g(i,j) + work_g(i,j-1)) ! dxt
         enddo
      enddo
      ! extrapolate to obtain dxt along j=1
      do i = 1, nx_global
         work_g2(i,1) = c2*work_g(i,2) - work_g(i,3) ! dxt
      enddo
      endif
      call scatter_global(dxt, work_g2, master_task, distrb_info, &
                          field_loc_center, field_type_scalar)

      deallocate(work_g2)

      end subroutine primary_grid_lengths_HTN

!=======================================================================
! Calculate dyu and dyt from HTE on the global grid, to preserve
! ghost cell and/or land values that might otherwise be lost. Scatter
! dyu, dyt and HTE to all processors.
!
! author: Elizabeth C. Hunke, LANL

      subroutine primary_grid_lengths_HTE(work_g)

      use ice_constants, only: p5, c2, cm_to_m, &
          field_loc_center, field_loc_NEcorner, &
          field_loc_Eface, field_type_scalar

      real (kind=dbl_kind), dimension(:,:) :: work_g ! global array holding HTE

      ! local variables

      integer (kind=int_kind) :: &
         i, j, &
         im1     ! i-1

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g2

      character(len=*), parameter :: subname = '(primary_grid_lengths_HTE)'

      if (my_task == master_task) then
         allocate(work_g2(nx_global,ny_global))
      else
         allocate(work_g2(1,1))
      endif

      if (my_task == master_task) then
         do j = 1, ny_global
         do i = 1, nx_global
            work_g(i,j) = work_g(i,j) * cm_to_m                ! HTE
         enddo
         enddo
         do j = 1, ny_global-1
         do i = 1, nx_global
            work_g2(i,j) = p5*(work_g(i,j) + work_g(i,j+1)) ! dyu
         enddo
         enddo
         ! extrapolate to obtain dyu along j=ny_global
         if (ny_global > 1) then
            do i = 1, nx_global
               work_g2(i,ny_global) = c2*work_g(i,ny_global-1) &
                                       - work_g(i,ny_global-2) ! dyu
            enddo
         endif
      endif
      call scatter_global(HTE, work_g, master_task, distrb_info, &
                          field_loc_Eface, field_type_scalar)
      call scatter_global(dyu, work_g2, master_task, distrb_info, &
                          field_loc_NEcorner, field_type_scalar)

      if (my_task == master_task) then
      do j = 1, ny_global
      do i = 1, nx_global
         ! assume cyclic; noncyclic will be handled during scatter
         im1 = i-1
         if (i == 1) im1 = nx_global 
         work_g2(i,j) = p5*(work_g(i,j) + work_g(im1,j))    ! dyt
      enddo
      enddo
      endif
      call scatter_global(dyt, work_g2, master_task, distrb_info, &
                          field_loc_center, field_type_scalar)

      deallocate(work_g2)

      end subroutine primary_grid_lengths_HTE

!=======================================================================

! Sets the boundary values for the T cell land mask (hm) and
! makes the logical land masks for T and U cells (tmask, umask).
! Also creates hemisphere masks (mask-n northern, mask-s southern)
!
! author: Elizabeth C. Hunke, LANL

      subroutine makemask

      use ice_constants, only: c0, p5, &
          field_loc_center, field_loc_NEcorner, field_type_scalar

      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      real (kind=dbl_kind) :: &
         puny

      type (block) :: &
         this_block           ! block information for current block

      character(len=*), parameter :: subname = '(makemask)'

      call icepack_query_parameters(puny_out=puny)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (kmt,               halo_info, &
                           field_loc_center, field_type_scalar)
      call ice_HaloUpdate (hm,               halo_info, &
                           field_loc_center, field_type_scalar)
      call ice_timer_stop(timer_bound)

      !-----------------------------------------------------------------
      ! construct T-cell and U-cell masks
      !-----------------------------------------------------------------

      bm = c0

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
            bm(i,j,iblk) = my_task + iblk/100.0_dbl_kind
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (uvm,                halo_info, &
                           field_loc_NEcorner, field_type_scalar)
      call ice_HaloUpdate (bm,               halo_info, &
                           field_loc_center, field_type_scalar)
      call ice_timer_stop(timer_bound)

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         ! needs to cover halo (no halo update for logicals)
         tmask(:,:,iblk) = .false.
         umask(:,:,iblk) = .false.
         do j = jlo-nghost, jhi+nghost
         do i = ilo-nghost, ihi+nghost
            if ( hm(i,j,iblk) > p5) tmask(i,j,iblk) = .true.
            if (uvm(i,j,iblk) > p5) umask(i,j,iblk) = .true.
         enddo
         enddo

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

      end subroutine makemask

!=======================================================================

! Initializes latitude and longitude on T grid
!
! author: Elizabeth C. Hunke, LANL; code originally based on POP grid
! generation routine

      subroutine Tlatlon

      use ice_constants, only: c0, c1, c2, c4, &
          field_loc_center, field_type_scalar

      integer (kind=int_kind) :: &
           i, j, iblk       , & ! horizontal indices
           ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      real (kind=dbl_kind) :: &
           z1,x1,y1,z2,x2,y2,z3,x3,y3,z4,x4,y4,tx,ty,tz,da, &
           rad_to_deg

      type (block) :: &
           this_block           ! block information for current block

      character(len=*), parameter :: subname = '(Tlatlon)'

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
      if (trim(grid_type) == 'regional') then
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
         if (nx_block > 5+2*nghost .and. ny_block > 5+2*nghost) then
         write(nu_diag,*) 'min/max ULON:', y1*rad_to_deg, y2*rad_to_deg
         write(nu_diag,*) 'min/max ULAT:', y3*rad_to_deg, y4*rad_to_deg
         endif
         write(nu_diag,*) 'min/max TLON:', x1*rad_to_deg, x2*rad_to_deg
         write(nu_diag,*) 'min/max TLAT:', x3*rad_to_deg, x4*rad_to_deg
      endif                     ! my_task

      end subroutine Tlatlon

!=======================================================================

! Transfer vector component from T-cell centers to U-cell centers.
!
! author: Elizabeth C. Hunke, LANL

      subroutine t2ugrid_vector (work)

      use ice_blocks, only: nx_block, ny_block
      use ice_constants, only: field_loc_center, field_type_vector
      use ice_domain_size, only: max_blocks

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks), intent(inout) :: & 
           work

      ! local variables

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1

      character(len=*), parameter :: subname = '(t2ugrid_vector)'

      work1(:,:,:) = work(:,:,:)

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (work1,            halo_info, &
                           field_loc_center, field_type_vector)
      call ice_timer_stop(timer_bound)

      call to_ugrid(work1,work)

      end subroutine t2ugrid_vector

!=======================================================================

! Shifts quantities from the T-cell midpoint (work1) to the U-cell
! midpoint (work2)
! NOTE: Input array includes ghost cells that must be updated before
!       calling this routine.
!
! author: Elizabeth C. Hunke, LANL

      subroutine to_ugrid(work1,work2)

      use ice_constants, only: c0, p25

      real (kind=dbl_kind), intent(in) :: &
         work1(nx_block,ny_block,max_blocks)

      real (kind=dbl_kind), intent(out) :: &
         work2(nx_block,ny_block,max_blocks)

      type (block) :: &
         this_block           ! block information for current block

      character(len=*), parameter :: subname = '(to_ugrid)'

      ! local variables

      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      work2(:,:,:) = c0

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
                              (work1(i,  j,  iblk)*tarea(i,  j,  iblk)  &
                             + work1(i+1,j,  iblk)*tarea(i+1,j,  iblk)  &
                             + work1(i,  j+1,iblk)*tarea(i,  j+1,iblk)  &
                             + work1(i+1,j+1,iblk)*tarea(i+1,j+1,iblk)) &
                             / uarea(i,  j,  iblk)
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine to_ugrid

!=======================================================================

! Transfer from U-cell centers to T-cell centers. Writes work into
! another array that has ghost cells
! NOTE: Input array is dimensioned only over physical cells.
!
! author: Elizabeth C. Hunke, LANL

      subroutine u2tgrid_vector (work)

      use ice_blocks, only: nx_block, ny_block
      use ice_constants, only: field_loc_NEcorner, field_type_vector
      use ice_domain_size, only: max_blocks

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), intent(inout) :: &
         work

      ! local variables

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1

      character(len=*), parameter :: subname = '(u2tgrid_vector)'

      work1(:,:,:) = work(:,:,:)

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (work1,              halo_info, &
                           field_loc_NEcorner, field_type_vector)
      call ice_timer_stop(timer_bound)

      call to_tgrid(work1,work)

      end subroutine u2tgrid_vector

!=======================================================================

! Shifts quantities from the U-cell midpoint (work1) to the T-cell
! midpoint (work2)
! NOTE: Input array includes ghost cells that must be updated before
!       calling this routine.
!
! author: Elizabeth C. Hunke, LANL

      subroutine to_tgrid(work1, work2)

      use ice_constants, only: p25

      real (kind=dbl_kind) :: work1(nx_block,ny_block,max_blocks), &
                              work2(nx_block,ny_block,max_blocks)

      ! local variables

      integer (kind=int_kind) :: &
         i, j, iblk, &
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      type (block) :: &
         this_block           ! block information for current block
      
      character(len=*), parameter :: subname = '(to_tgrid)'

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
                             (work1(i,  j  ,iblk) * uarea(i,  j,  iblk)  &
                            + work1(i-1,j  ,iblk) * uarea(i-1,j,  iblk)  &
                            + work1(i,  j-1,iblk) * uarea(i,  j-1,iblk)  & 
                            + work1(i-1,j-1,iblk) * uarea(i-1,j-1,iblk)) &
                            / tarea(i,  j,  iblk)
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine to_tgrid

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

      use ice_blocks, only: nx_block, ny_block
      use ice_constants, only: c0,  c2, c360, &
          field_loc_NEcorner, field_type_scalar
      use ice_domain_size, only: max_blocks

      integer (kind=int_kind) :: &
          i,j,iblk,icorner,& ! index counters
          ilo,ihi,jlo,jhi    ! beginning and end of physical domain

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
         allocate(work_g2(nx_global,ny_global))
      else
         allocate(work_g2(1,1))
      endif

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

      deallocate(work_g2)

      !----------------------------------------------------------------
      ! Convert longitude to Degrees East >0 for history output
      !----------------------------------------------------------------

      allocate(work_g2(nx_block,ny_block))  ! not used as global here
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
      deallocate(work_g2)

      end subroutine gridbox_corners

!=======================================================================

! NOTE:  Boundary conditions for fields on NW, SW, SE corners
!        have not been implemented; using NE corner location for all.
!        Extrapolations are also used: these fields are approximate!
!
! authors:   A. McLaren, Met Office
!            E. Hunke, LANL

      subroutine gridbox_verts(work_g,vbounds)

      use ice_blocks, only: nx_block, ny_block
      use ice_constants, only: c0, c2, &
          field_loc_NEcorner, field_type_scalar
      use ice_domain_size, only: max_blocks

      real (kind=dbl_kind), dimension(:,:), intent(in) :: &
          work_g

      real (kind=dbl_kind), dimension(4,nx_block,ny_block,max_blocks), intent(out) :: &
          vbounds

      integer (kind=int_kind) :: &
          i,j                 ! index counters

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
         allocate(work_g2(nx_global,ny_global))
      else
         allocate(work_g2(1,1))
      endif

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

      deallocate (work_g2)

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

      real (kind=dbl_kind) :: &
         puny

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

      call icepack_query_parameters(puny_out=puny)
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

         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               k = min(nint(kmt(i,j,iblk)),nlevel)
               if (k > nlevel) call abort_ice(subname//' kmt gt nlevel error')
               if (k > 0) bathymetry(i,j,iblk) = depth(k)
            enddo
            enddo
         enddo

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

      character(len=*), parameter :: subname = '(get_bathymetry_popfile)'

      ntmp = maxval(nint(KMT))
      nlevel = global_maxval(ntmp,distrb_info)

      if (my_task==master_task) then
         write(nu_diag,*) subname,' KMT max = ',nlevel
      endif

      allocate(depth(nlevel),thick(nlevel))
      thick = -999999.
      depth = -999999.

      if (use_bathymetry) then

         write (nu_diag,*) subname,' Bathymetry file = ', trim(bathymetry_file)
         if (my_task == master_task) then
            call get_fileunit(fid)
            open(fid,file=bathymetry_file,form='formatted',iostat=ierr)
            if (ierr/=0) call abort_ice(subname//' open error')
            do k = 1,nlevel
               read(fid,*,iostat=ierr) thick(k)
               if (ierr/=0) call abort_ice(subname//' read error')
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
         if (depth(k) < 0.) call abort_ice(subname//' negative depth error')
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
            if (k > nlevel) call abort_ice(subname//' kmt gt nlevel error')
            if (k > 0) bathymetry(i,j,iblk) = depth(k)
         enddo
         enddo
      enddo

      deallocate(depth,thick)

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
      use ice_constants, only: field_loc_center, field_type_scalar

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
         write(nu_diag,*) 'reading ',TRIM(fieldname)
         call icepack_warnings_flush(nu_diag)
      endif
      call ice_read_nc(fid_init,1,fieldname,bathymetry,diag, &
                    field_loc=field_loc_center, &
                    field_type=field_type_scalar)

      call ice_close_nc(fid_init)

      if (my_task == master_task) then
         write(nu_diag,*) 'closing file ',TRIM(bathymetry_file)
         call icepack_warnings_flush(nu_diag)
      endif

      end subroutine read_seabedstress_bathy
      
!=======================================================================

      end module ice_grid

!=======================================================================
