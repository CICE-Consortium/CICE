!=======================================================================

! Hemispheric diagnostics written to a netcdf file (diag_netcdf), as time
! averages over each diagfreq interval (diag_netcdf_avg) or instantaneous values.
!
! The diagnostics are global sums or maxima of per grid cell values, using the
! same global reductions (and bfbflag) as the diagnostic log.  For time averages
! the per grid cell values are accumulated every timestep, so the only global
! reductions are when a record is written.
!
! authors: Anton Steketee, ACCESS-NRI

      module ice_diagnostics_netcdf

      use ice_kinds_mod
      use ice_blocks, only: nx_block, ny_block
      use ice_communicate, only: my_task, master_task
      use ice_constants, only: c0, c1
      use ice_diagnostics, only: aice_extmin, maxval_spval, &
          diag_kinetic_energy, diag_ice_speed
      use ice_exit, only: abort_ice
      use ice_fileunits, only: nu_diag
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_parameters

      implicit none
      private
      public :: diag_netcdf_step

      logical (kind=log_kind), public :: &
         diag_netcdf      , & ! if true, write hemispheric diagnostics to a netcdf file
         diag_netcdf_avg      ! if true, time averages, otherwise instantaneous values

      !-----------------------------------------------------------------
      ! Hemispheric netcdf diagnostics (diag_netcdf)
      !
      ! Each diagnostic is a per grid cell value that is reduced to a value
      ! for each hemisphere, either an area weighted global sum or a global
      ! maximum.  The same per grid cell values are used for instantaneous
      ! output and, accumulated every timestep, for time averages.
      !
      ! To add a diagnostic:
      !   1. add an index ncd_* below and increase n_ncdfld
      !   2. set its names, units and reduction in init_ncdiag
      !   3. compute its per grid cell value in ncdiag_cell_values
      !-----------------------------------------------------------------

      integer (kind=int_kind), parameter :: &
         ncd_aice   = 1   , & ! sea-ice area
         ncd_extent = 2   , & ! sea-ice extent
         ncd_vice   = 3   , & ! sea-ice volume
         ncd_snwm   = 4   , & ! snow mass
         ncd_ke     = 5   , & ! kinetic energy
         ncd_albn   = 6   , & ! albedo numerator, for derived albedo
         ncd_albd   = 7   , & ! albedo denominator, for derived albedo
         ncd_vmax   = 8   , & ! max ice volume per unit area
         ncd_umax   = 9   , & ! max ice speed
         ncd_pmax   = 10  , & ! max ice strength
         n_ncdfld   = 10      ! number of per grid cell diagnostics

      type ncdiag_field
         character (len=char_len) :: &
            vname_n     , & ! netcdf variable name, northern hemisphere
            vname_s     , & ! netcdf variable name, southern hemisphere
            lname       , & ! long name
            units       , & ! units after scale is applied
            sname           ! standard_name, blank for none
         character (len=4) :: &
            reduce          ! 'sum' area weighted global sum or 'max' global maximum
         real (kind=dbl_kind) :: &
            scale           ! unit conversion applied to the hemispheric value
         logical (kind=log_kind) :: &
            lwrite          ! written to file, false if only used for derived diagnostics
      end type ncdiag_field

      type (ncdiag_field), dimension(n_ncdfld) :: &
         ncdfld               ! per grid cell diagnostics

      logical (kind=log_kind) :: &
         ncdfld_init = .false. ! ncdfld has been initialized

      logical (kind=log_kind), dimension(2) :: &
         ncdiag_hemi          ! hemisphere has grid cells, otherwise values are missing

      ! accumulation for time averages
      real (kind=dbl_kind), dimension(:,:,:,:), allocatable :: &
         acc_ncdiag           ! accumulated per grid cell diagnostics

      integer (kind=int_kind) :: &
         nacc_steps = 0       ! number of accumulated timesteps

      real (kind=dbl_kind) :: &
         acc_time_beg     , & ! start of current accumulation interval (s)
         ncdiag_time_beg      ! start of the interval being written (s)

      ! output list, the values written in each record
      integer (kind=int_kind), parameter :: &
         max_ncdiag = 64      ! max number of output variables per hemisphere

      integer (kind=int_kind) :: &
         n_ncdiag    = 0  , & ! number of output variables per hemisphere
         nrec_ncdiag = 0      ! number of records in the netcdf diagnostics file

      logical (kind=log_kind) :: &
         ncdiag_opened = .false. ! netcdf diagnostics file has been opened this run

      character (len=char_len), dimension(max_ncdiag,2) :: &
         ncdiag_name          ! netcdf variable name, northern and southern hemisphere

      character (len=char_len), dimension(max_ncdiag) :: &
         ncdiag_lname     , & ! netcdf long name
         ncdiag_units     , & ! netcdf units
         ncdiag_sname     , & ! netcdf standard_name, blank for none
         ncdiag_cmethods      ! netcdf cell_methods, blank for instantaneous

      real (kind=dbl_kind), dimension(max_ncdiag) :: &
         ncdiag_valn      , & ! northern hemisphere value
         ncdiag_vals          ! southern hemisphere value

      character (len=char_len_long) :: &
         ncdiag_file          ! netcdf diagnostics filename

!=======================================================================

      contains

!=======================================================================

! Hemispheric netcdf diagnostics for this timestep, called every timestep
! from the drivers.  Accumulates time averages and writes a record when
! the other diagnostics are written (write_diag).

      subroutine diag_netcdf_step

      use ice_calendar, only: write_diag

      character(len=*), parameter :: subname = '(diag_netcdf_step)'

      if (.not. diag_netcdf) return

      if (diag_netcdf_avg) call accum_ncdiag()

      if (write_diag) then
         call compute_ncdiag()     ! global reductions on all tasks
         call write_ncdiag()       ! master task
      endif

      end subroutine diag_netcdf_step

!=======================================================================

! Define the per grid cell hemispheric netcdf diagnostics, see the
! comments with ncdfld in the module header

      subroutine init_ncdiag

      use ice_constants, only: field_loc_center
      use ice_domain, only: distrb_info
      use ice_domain_size, only: max_blocks
      use ice_global_reductions, only: global_sum
      use ice_grid, only: tarean, tareas

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1             ! temporary

      character(len=*), parameter :: subname = '(init_ncdiag)'

      ! a hemisphere without grid cells (e.g. regional domains) is written as missing
      work1(:,:,:) = c1
      ncdiag_hemi(1) = global_sum(work1, distrb_info, field_loc_center, tarean) > c0
      ncdiag_hemi(2) = global_sum(work1, distrb_info, field_loc_center, tareas) > c0

      ! CMIP sea-ice diagnostics, in CMIP units
      call def_ncdfld(ncd_aice, 'siarean'  , 'siareas'  , 'sea-ice area'  , '1e6 km2', 'sum', &
                      1.0e-12_dbl_kind, sname='sea_ice_area')
      call def_ncdfld(ncd_extent, 'siextentn', 'siextents', 'sea-ice extent', '1e6 km2', 'sum', &
                      1.0e-12_dbl_kind, sname='sea_ice_extent')
      call def_ncdfld(ncd_vice, 'sivoln'   , 'sivols'   , 'sea-ice volume', '1e3 km3', 'sum', &
                      1.0e-12_dbl_kind, sname='sea_ice_volume')
      call def_ncdfld(ncd_snwm, 'sisnmassn', 'sisnmasss', 'snow mass on sea ice', 'kg', 'sum', c1)

      ! other diagnostics, the same quantities as in the log
      call def_ncdfld(ncd_ke  , 'kinetic_energy_nh'  , 'kinetic_energy_sh'  , &
                      'total kinetic energy', 'J', 'sum', c1)
      call def_ncdfld(ncd_albn, 'albedo_num_nh'      , 'albedo_num_sh'      , &
                      'albedo numerator', '1', 'sum', c1, lwrite=.false.)
      call def_ncdfld(ncd_albd, 'albedo_den_nh'      , 'albedo_den_sh'      , &
                      'albedo denominator', '1', 'sum', c1, lwrite=.false.)
      call def_ncdfld(ncd_vmax, 'ice_volume_max_nh'  , 'ice_volume_max_sh'  , &
                      'max ice volume per unit area', 'm', 'max', c1)
      call def_ncdfld(ncd_umax, 'ice_speed_max_nh'   , 'ice_speed_max_sh'   , &
                      'max ice speed', 'm/s', 'max', c1)
      call def_ncdfld(ncd_pmax, 'ice_strength_max_nh', 'ice_strength_max_sh', &
                      'max ice strength', 'kN/m', 'max', 1.0e-3_dbl_kind)

      ncdfld_init = .true.

      end subroutine init_ncdiag

!=======================================================================

! Set the definition of one per grid cell hemispheric netcdf diagnostic

      subroutine def_ncdfld(k, vname_n, vname_s, lname, units, reduce, scale, sname, lwrite)

      integer (kind=int_kind), intent(in) :: &
         k               ! index in ncdfld

      character (len=*), intent(in) :: &
         vname_n     , & ! netcdf variable name, northern hemisphere
         vname_s     , & ! netcdf variable name, southern hemisphere
         lname       , & ! long name
         units       , & ! units after scale is applied
         reduce          ! 'sum' area weighted global sum or 'max' global maximum

      real (kind=dbl_kind), intent(in) :: &
         scale           ! unit conversion applied to the hemispheric value

      character (len=*), intent(in), optional :: &
         sname           ! standard_name

      logical (kind=log_kind), intent(in), optional :: &
         lwrite          ! written to file, default true

      character(len=*), parameter :: subname = '(def_ncdfld)'

      if (reduce /= 'sum' .and. reduce /= 'max') then
         call abort_ice(subname//' ERROR: reduce must be sum or max for '//trim(vname_n), &
            file=__FILE__, line=__LINE__)
      endif

      ncdfld(k)%vname_n = vname_n
      ncdfld(k)%vname_s = vname_s
      ncdfld(k)%lname   = lname
      ncdfld(k)%units   = units
      ncdfld(k)%reduce  = reduce
      ncdfld(k)%scale   = scale
      ncdfld(k)%sname   = ''
      if (present(sname)) ncdfld(k)%sname = sname
      ncdfld(k)%lwrite  = .true.
      if (present(lwrite)) ncdfld(k)%lwrite = lwrite

      end subroutine def_ncdfld

!=======================================================================

! Compute the per grid cell values of the hemispheric netcdf diagnostics.
! Values to be summed are per unit area; they are multiplied by the cell
! area when reduced.

      subroutine ncdiag_cell_values(fld)

      use ice_domain, only: nblocks
      use ice_domain_size, only: max_blocks
      use ice_flux, only: alvdr, alidr, alvdf, alidf
      use ice_state, only: aice, vice, vsno, strength

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks,n_ncdfld), intent(out) :: &
         fld               ! per grid cell values

      ! local variables

      integer (kind=int_kind) :: &
         i, j, iblk

      real (kind=dbl_kind) :: &
         puny, rhos, awtvdr, awtidr, awtvdf, awtidf, alb

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1, work2      ! temporary

      character(len=*), parameter :: subname = '(ncdiag_cell_values)'

      call icepack_query_parameters(puny_out=puny, rhos_out=rhos, &
           awtvdr_out=awtvdr, awtidr_out=awtidr, awtvdf_out=awtvdf, awtidf_out=awtidf)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      call diag_kinetic_energy(work1)
      call diag_ice_speed(work2)

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,alb)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            fld(i,j,iblk,ncd_aice) = aice(i,j,iblk)
            fld(i,j,iblk,ncd_extent) = c0
            if (aice(i,j,iblk) >= aice_extmin) fld(i,j,iblk,ncd_extent) = c1
            fld(i,j,iblk,ncd_vice) = vice(i,j,iblk)
            fld(i,j,iblk,ncd_snwm) = rhos*vsno(i,j,iblk)
            fld(i,j,iblk,ncd_ke  ) = work1(i,j,iblk)
            ! ice area weighted albedo, masked where the sun is below the horizon
            alb = alvdr(i,j,iblk)*awtvdr + alidr(i,j,iblk)*awtidr &
                + alvdf(i,j,iblk)*awtvdf + alidf(i,j,iblk)*awtidf
            fld(i,j,iblk,ncd_albn) = c0
            fld(i,j,iblk,ncd_albd) = c0
            if (alb > puny) then
               fld(i,j,iblk,ncd_albn) = aice(i,j,iblk)*alb
               fld(i,j,iblk,ncd_albd) = aice(i,j,iblk)
            endif
            fld(i,j,iblk,ncd_vmax) = vice(i,j,iblk)
            fld(i,j,iblk,ncd_umax) = work2(i,j,iblk)
            fld(i,j,iblk,ncd_pmax) = strength(i,j,iblk)
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine ncdiag_cell_values

!=======================================================================

! Accumulate the per grid cell hemispheric netcdf diagnostics every
! timestep for time averages.  This is local work only, the global
! reductions are done once per write in compute_ncdiag.

      subroutine accum_ncdiag

      use ice_calendar, only: timesecs, dt
      use ice_domain_size, only: max_blocks

      ! local variables

      integer (kind=int_kind) :: &
         k

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable :: &
         fld               ! per grid cell values

      character(len=*), parameter :: subname = '(accum_ncdiag)'

      if (.not. ncdfld_init) call init_ncdiag()

      if (.not. allocated(acc_ncdiag)) then
         allocate(acc_ncdiag(nx_block,ny_block,max_blocks,n_ncdfld))
         call reset_ncdiag()
      endif

      ! start of the averaging interval is the start of the first step
      if (nacc_steps == 0) acc_time_beg = timesecs - dt
      nacc_steps = nacc_steps + 1

      allocate(fld(nx_block,ny_block,max_blocks,n_ncdfld))
      call ncdiag_cell_values(fld)

      do k = 1, n_ncdfld
         if (ncdfld(k)%reduce == 'max') then
            acc_ncdiag(:,:,:,k) = max(acc_ncdiag(:,:,:,k), fld(:,:,:,k))
         else
            acc_ncdiag(:,:,:,k) = acc_ncdiag(:,:,:,k) + fld(:,:,:,k)
         endif
      enddo

      deallocate(fld)

      end subroutine accum_ncdiag

!=======================================================================

! Reset the accumulation for time averages

      subroutine reset_ncdiag

      integer (kind=int_kind) :: &
         k

      do k = 1, n_ncdfld
         if (ncdfld(k)%reduce == 'max') then
            acc_ncdiag(:,:,:,k) = -huge(c1)
         else
            acc_ncdiag(:,:,:,k) = c0
         endif
      enddo
      nacc_steps = 0

      end subroutine reset_ncdiag

!=======================================================================

! Compute the hemispheric netcdf diagnostics for this record, time averages
! from the accumulation or instantaneous values, and fill the output list.
! Called on all tasks as it does global reductions.

      subroutine compute_ncdiag

      use ice_calendar, only: timesecs
      use ice_constants, only: c2, field_loc_center
      use ice_domain, only: distrb_info
      use ice_domain_size, only: max_blocks
      use ice_global_reductions, only: global_sum, global_maxval
      use ice_grid, only: lmask_n, lmask_s, tarean, tareas

      ! local variables

      integer (kind=int_kind) :: &
         k

      real (kind=dbl_kind) :: &
         rn, puny, rhoi, urmsn, urmss, albn, albs

      real (kind=dbl_kind), dimension(n_ncdfld) :: &
         valn, vals        ! hemispheric values, before scale is applied

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable :: &
         fld               ! per grid cell values

      character (len=char_len) :: &
         cm_sum, cm_max    ! cell_methods for sums and maxima

      character(len=*), parameter :: subname = '(compute_ncdiag)'

      call icepack_query_parameters(puny_out=puny, rhoi_out=rhoi)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (.not. ncdfld_init) call init_ncdiag()

      allocate(fld(nx_block,ny_block,max_blocks,n_ncdfld))
      if (diag_netcdf_avg) then
         if (nacc_steps == 0) then
            call abort_ice(subname//' ERROR: no accumulated diagnostics', &
               file=__FILE__, line=__LINE__)
         endif
         fld(:,:,:,:) = acc_ncdiag(:,:,:,:)
         rn = c1/real(nacc_steps,kind=dbl_kind)
         cm_sum = 'time: mean'
         cm_max = 'time: maximum'
      else
         call ncdiag_cell_values(fld)
         rn = c1
         cm_sum = ''
         cm_max = ''
      endif

      do k = 1, n_ncdfld
         if (ncdfld(k)%reduce == 'max') then
            valn(k) = global_maxval(fld(:,:,:,k), distrb_info, lmask_n)
            vals(k) = global_maxval(fld(:,:,:,k), distrb_info, lmask_s)
            if (valn(k) < maxval_spval) valn(k) = c0
            if (vals(k) < maxval_spval) vals(k) = c0
         else
            valn(k) = global_sum(fld(:,:,:,k), distrb_info, field_loc_center, tarean)*rn
            vals(k) = global_sum(fld(:,:,:,k), distrb_info, field_loc_center, tareas)*rn
         endif
      enddo
      deallocate(fld)

      n_ncdiag = 0
      do k = 1, n_ncdfld
         if (ncdfld(k)%lwrite) then
            if (ncdfld(k)%reduce == 'max') then
               call add_ncdiag(ncdfld(k)%vname_n, ncdfld(k)%vname_s, ncdfld(k)%lname, &
                               ncdfld(k)%units, valn(k)*ncdfld(k)%scale, vals(k)*ncdfld(k)%scale, &
                               ncdfld(k)%sname, cm_max)
            else
               call add_ncdiag(ncdfld(k)%vname_n, ncdfld(k)%vname_s, ncdfld(k)%lname, &
                               ncdfld(k)%units, valn(k)*ncdfld(k)%scale, vals(k)*ncdfld(k)%scale, &
                               ncdfld(k)%sname, cm_sum)
            endif
         endif
      enddo

      ! derived diagnostics

      ! rms ice speed from kinetic energy and ice and snow mass
      urmsn = c2*valn(ncd_ke)/(rhoi*valn(ncd_vice) + valn(ncd_snwm) + puny)
      urmss = c2*vals(ncd_ke)/(rhoi*vals(ncd_vice) + vals(ncd_snwm) + puny)
      if (urmsn > puny) then
         urmsn = sqrt(urmsn)
      else
         urmsn = c0
      endif
      if (urmss > puny) then
         urmss = sqrt(urmss)
      else
         urmss = c0
      endif
      call add_ncdiag('ice_speed_rms_nh', 'ice_speed_rms_sh', 'rms ice speed', 'm/s', &
                      urmsn, urmss, '', cm_sum)

      ! ice area weighted albedo over sunlit points
      albn = c0
      albs = c0
      if (valn(ncd_albd) > c0) albn = valn(ncd_albn)/valn(ncd_albd)
      if (vals(ncd_albd) > c0) albs = vals(ncd_albn)/vals(ncd_albd)
      call add_ncdiag('albedo_nh', 'albedo_sh', 'average albedo', '1', &
                      albn, albs, '', cm_sum)

      if (diag_netcdf_avg) then
         ncdiag_time_beg = acc_time_beg
         call reset_ncdiag()
      endif

      end subroutine compute_ncdiag

!=======================================================================

! Add a hemispheric diagnostic to the list written to the netcdf
! diagnostics file

      subroutine add_ncdiag(vname_n, vname_s, lname, units, valn, vals, sname, cmethods)

      character (len=*), intent(in) :: &
         vname_n     , & ! variable name, northern hemisphere
         vname_s     , & ! variable name, southern hemisphere
         lname       , & ! long name
         units       , & ! units
         sname       , & ! standard_name, blank for none
         cmethods        ! cell_methods, blank for instantaneous

      real (kind=dbl_kind), intent(in) :: &
         valn        , & ! northern hemisphere value
         vals            ! southern hemisphere value

      character(len=*), parameter :: subname = '(add_ncdiag)'

      n_ncdiag = n_ncdiag + 1
      if (n_ncdiag > max_ncdiag) then
         call abort_ice(subname//' ERROR: too many diagnostics, increase max_ncdiag', &
            file=__FILE__, line=__LINE__)
      endif

      ncdiag_name(n_ncdiag,1) = vname_n
      ncdiag_name(n_ncdiag,2) = vname_s
      ncdiag_lname   (n_ncdiag) = lname
      ncdiag_units   (n_ncdiag) = units
      ncdiag_sname   (n_ncdiag) = sname
      ncdiag_cmethods(n_ncdiag) = cmethods
      ncdiag_valn    (n_ncdiag) = valn
      ncdiag_vals    (n_ncdiag) = vals

      end subroutine add_ncdiag

!=======================================================================

! Write hemispheric diagnostics to the netcdf diagnostics file, master task only.
! There is one file for each of time averages and instantaneous values.  It is
! created on the first write, or appended to if it already exists, and reopened
! for each record so it is complete if the run stops early.

      subroutine write_ncdiag

#ifdef USE_NETCDF
      use netcdf
#endif
      use ice_calendar, only: timesecs, idate0, hh_init, mm_init, ss_init, &
          days_per_year, use_leap_years, seconds_per_day, dt
      use ice_history_shared, only: history_dir, history_file, history_format, &
          history_deflate, history_precision
      use ice_constants, only: spval, spval_dbl, p5
      use ice_read_write, only: ice_check_nc

      ! local variables

      integer (kind=int_kind) :: &
         n              , & ! diagnostic index
         hemi           , & ! hemisphere index, 1 north, 2 south
         status         , & ! netcdf status
         iflag          , & ! netcdf create mode
         lprecision     , & ! netcdf precision of diagnostic variables
         ncid           , & ! netcdf file id
         timid          , & ! netcdf time dimension id
         bndid          , & ! netcdf bounds dimension id
         varid              ! netcdf variable id

      character (len=8) :: &
         cdate          , & ! initial date string
         current_date       ! wall clock date

      character (len=10) :: &
         current_time       ! wall clock time

      character (len=char_len) :: &
         vname          , & ! netcdf variable name
         cal_units      , & ! time units
         cal_att        , & ! calendar attribute
         hist_att           ! history attribute

      real (kind=dbl_kind) :: &
         val            , & ! value to write
         tcur               ! half a timestep before the current time (days)

      real (kind=dbl_kind), dimension(:), allocatable :: &
         ftime              ! times already in the file (days)

      integer (kind=int_kind) :: &
         ntime              ! number of records already in the file

      logical (kind=log_kind) :: &
         lexist             ! file exists

      character (len=19), dimension(2), parameter :: &
         hemi_lname = (/'northern hemisphere','southern hemisphere'/)

      character(len=*), parameter :: subname = '(write_ncdiag)'

#ifdef USE_NETCDF
      if (my_task /= master_task) return

      ! filename, _mean for time averages or _inst for instantaneous values,
      ! as for history files
      if (diag_netcdf_avg) then
         ncdiag_file = trim(history_dir)//trim(history_file)//'_diag_mean.nc'
      else
         ncdiag_file = trim(history_dir)//trim(history_file)//'_diag_inst.nc'
      endif
      lexist = .true.
      if (.not. ncdiag_opened) inquire(file=trim(ncdiag_file), exist=lexist)

      ! create the file if it does not exist, otherwise open it to append a record
      if (.not. lexist) then

         write(nu_diag,*) subname,' creating ',trim(ncdiag_file)

         ! file format consistent with history files, written serially
         if (history_format == 'cdf2' .or. history_format == 'pnetcdf2') then
#ifdef NO_CDF2
            call abort_ice(subname//' ERROR: history_format cdf2 not available ', &
                           file=__FILE__, line=__LINE__)
#else
            iflag = ior(nf90_clobber,nf90_64bit_offset)
#endif
         elseif (history_format == 'cdf5' .or. history_format == 'pnetcdf5') then
#ifdef NO_CDF5
            call abort_ice(subname//' ERROR: history_format cdf5 not available ', &
                           file=__FILE__, line=__LINE__)
#else
            iflag = ior(nf90_clobber,nf90_64bit_data)
#endif
         elseif (history_format == 'hdf5') then
#ifdef NO_HDF5
            call abort_ice(subname//' ERROR: history_format hdf5 not available ', &
                           file=__FILE__, line=__LINE__)
#else
            iflag = ior(nf90_clobber,nf90_netcdf4)
#endif
         else
            iflag = nf90_clobber
         endif

         ! precision consistent with history files, time is always double
         lprecision = nf90_float
         if (history_precision == 8) lprecision = nf90_double

         ! create the file
         status = nf90_create(trim(ncdiag_file), iflag, ncid)
         call ice_check_nc(status, subname//' ERROR: creating '//trim(ncdiag_file), &
                           file=__FILE__, line=__LINE__)

         ! define the time coordinate
         status = nf90_def_dim(ncid, 'time', nf90_unlimited, timid)
         call ice_check_nc(status, subname//' ERROR: defining dim time', &
                           file=__FILE__, line=__LINE__)
         if (diag_netcdf_avg) then
            status = nf90_def_dim(ncid, 'nbnd', 2, bndid)
            call ice_check_nc(status, subname//' ERROR: defining dim nbnd', &
                              file=__FILE__, line=__LINE__)
         endif

         write(cdate,'(i8.8)') idate0
         write(cal_units,'(a,a4,a1,a2,a1,a2,a1,i2.2,a1,i2.2,a1,i2.2)') 'days since ', &
               cdate(1:4),'-',cdate(5:6),'-',cdate(7:8),' ', &
               hh_init,':',mm_init,':',ss_init

         if (days_per_year == 360) then
            cal_att = '360_day'
         elseif (days_per_year == 365 .and. .not.use_leap_years ) then
            cal_att = 'noleap'
         elseif (use_leap_years) then
            cal_att = 'proleptic_gregorian'
         else
            call abort_ice(subname//' ERROR: invalid calendar settings', &
                           file=__FILE__, line=__LINE__)
         endif

         status = nf90_def_var(ncid, 'time', nf90_double, (/timid/), varid)
         call ice_check_nc(status, subname//' ERROR: defining var time', &
                           file=__FILE__, line=__LINE__)
         status = nf90_put_att(ncid, varid, 'long_name', 'time')
         call ice_check_nc(status, subname//' ERROR: defining att long_name time', &
                           file=__FILE__, line=__LINE__)
         status = nf90_put_att(ncid, varid, 'standard_name', 'time')
         call ice_check_nc(status, subname//' ERROR: defining att standard_name time', &
                           file=__FILE__, line=__LINE__)
         status = nf90_put_att(ncid, varid, 'axis', 'T')
         call ice_check_nc(status, subname//' ERROR: defining att axis time', &
                           file=__FILE__, line=__LINE__)
         status = nf90_put_att(ncid, varid, 'units', trim(cal_units))
         call ice_check_nc(status, subname//' ERROR: defining att units time', &
                           file=__FILE__, line=__LINE__)
         status = nf90_put_att(ncid, varid, 'calendar', trim(cal_att))
         call ice_check_nc(status, subname//' ERROR: defining att calendar time', &
                           file=__FILE__, line=__LINE__)
         if (use_leap_years) then
            ! model time has no leap seconds (CF-1.11)
            status = nf90_put_att(ncid, varid, 'units_metadata', 'leap_seconds: none')
            call ice_check_nc(status, subname//' ERROR: defining att units_metadata time', &
                              file=__FILE__, line=__LINE__)
         endif

         if (diag_netcdf_avg) then
            status = nf90_put_att(ncid, varid, 'bounds', 'time_bounds')
            call ice_check_nc(status, subname//' ERROR: defining att bounds time', &
                              file=__FILE__, line=__LINE__)
            status = nf90_def_var(ncid, 'time_bounds', nf90_double, (/bndid,timid/), varid)
            call ice_check_nc(status, subname//' ERROR: defining var time_bounds', &
                              file=__FILE__, line=__LINE__)
            ! CF bounds variables inherit their attributes from time
         endif

         ! define each variable for the hemispheric diagnostics, with attributes
         do n = 1, n_ncdiag
            do hemi = 1, 2
               vname = ncdiag_name(n,hemi)
               status = nf90_def_var(ncid, trim(vname), lprecision, (/timid/), varid)
               call ice_check_nc(status, subname//' ERROR: defining var '//trim(vname), &
                                 file=__FILE__, line=__LINE__)
   #ifndef NO_HDF5
               if (history_format == 'hdf5' .and. history_deflate /= 0) then
                  status = nf90_def_var_deflate(ncid, varid, shuffle=0, deflate=1, &
                                                deflate_level=history_deflate)
                  call ice_check_nc(status, subname//' ERROR: deflating var '//trim(vname), &
                                    file=__FILE__, line=__LINE__)
               endif
   #endif
               status = nf90_put_att(ncid, varid, 'long_name', &
                                    trim(ncdiag_lname(n))//', '//hemi_lname(hemi))
               call ice_check_nc(status, subname//' ERROR: defining att long_name '//trim(vname), &
                                 file=__FILE__, line=__LINE__)
               status = nf90_put_att(ncid, varid, 'units', trim(ncdiag_units(n)))
               call ice_check_nc(status, subname//' ERROR: defining att units '//trim(vname), &
                                 file=__FILE__, line=__LINE__)
               if (lprecision == nf90_double) then
                  status = nf90_put_att(ncid, varid, '_FillValue', spval_dbl)
               else
                  status = nf90_put_att(ncid, varid, '_FillValue', spval)
               endif
               call ice_check_nc(status, subname//' ERROR: defining att _FillValue '//trim(vname), &
                                 file=__FILE__, line=__LINE__)
               if (len_trim(ncdiag_sname(n)) > 0) then
                  status = nf90_put_att(ncid, varid, 'standard_name', trim(ncdiag_sname(n)))
                  call ice_check_nc(status, subname//' ERROR: defining att standard_name '//trim(vname), &
                                    file=__FILE__, line=__LINE__)
               endif
               if (len_trim(ncdiag_cmethods(n)) > 0) then
                  status = nf90_put_att(ncid, varid, 'cell_methods', trim(ncdiag_cmethods(n)))
                  call ice_check_nc(status, subname//' ERROR: defining att cell_methods '//trim(vname), &
                                    file=__FILE__, line=__LINE__)
               endif
            enddo
         enddo

         ! global attributes, consistent with history files
         status = nf90_put_att(ncid, nf90_global, 'title', 'CICE hemispheric diagnostics')
         call ice_check_nc(status, subname//' ERROR: defining global att title', &
                           file=__FILE__, line=__LINE__)
         if (diag_netcdf_avg) then
            status = nf90_put_att(ncid, nf90_global, 'comment', &
                                  'time averages (or maxima) over each diagfreq interval')
         else
            status = nf90_put_att(ncid, nf90_global, 'comment', &
                                  'instantaneous values at diagfreq, as written to the diagnostic log')
         endif
         call ice_check_nc(status, subname//' ERROR: defining global att comment', &
                           file=__FILE__, line=__LINE__)
         status = nf90_put_att(ncid, nf90_global, 'Conventions', 'CF-1.8')
         call ice_check_nc(status, subname//' ERROR: defining global att Conventions', &
                           file=__FILE__, line=__LINE__)

         call date_and_time(date=current_date, time=current_time)
         write(hist_att,'(a,a4,a1,a2,a1,a2,a4,a2,a1,a2,a1,a2)') 'This dataset was created on ', &
               current_date(1:4),'-',current_date(5:6),'-',current_date(7:8),' at ', &
               current_time(1:2),':',current_time(3:4),':',current_time(5:6)
         status = nf90_put_att(ncid, nf90_global, 'history', trim(hist_att))
         call ice_check_nc(status, subname//' ERROR: defining global att history', &
                           file=__FILE__, line=__LINE__)

         status = nf90_enddef(ncid)
         call ice_check_nc(status, subname//' ERROR: in nf90_enddef', &
                           file=__FILE__, line=__LINE__)

      else ! file already exists, open it to append a record

         status = nf90_open(trim(ncdiag_file), nf90_write, ncid)
         call ice_check_nc(status, subname//' ERROR: opening '//trim(ncdiag_file), &
                           file=__FILE__, line=__LINE__)

         if (.not. ncdiag_opened) then
            ! first write of this run to an existing file, append after the last
            ! record before the current time, so records from a rerun are overwritten
            write(nu_diag,*) subname,' appending to ',trim(ncdiag_file)
            status = nf90_inq_dimid(ncid, 'time', timid)
            call ice_check_nc(status, subname//' ERROR: getting dimid for time', &
                              file=__FILE__, line=__LINE__)
            status = nf90_inquire_dimension(ncid, timid, len=ntime)
            call ice_check_nc(status, subname//' ERROR: getting length of time', &
                              file=__FILE__, line=__LINE__)
            nrec_ncdiag = ntime
            if (ntime > 0) then
               allocate(ftime(ntime))
               status = nf90_inq_varid(ncid, 'time', varid)
               call ice_check_nc(status, subname//' ERROR: getting varid for time', &
                                 file=__FILE__, line=__LINE__)
               status = nf90_get_var(ncid, varid, ftime)
               call ice_check_nc(status, subname//' ERROR: reading time', &
                                 file=__FILE__, line=__LINE__)
               ! records at or after the current time, within half a timestep
               tcur = (timesecs - p5*dt)/real(seconds_per_day,kind=dbl_kind)
               do n = 1, ntime
                  if (ftime(n) > tcur) then
                     nrec_ncdiag = n - 1
                     exit
                  endif
               enddo
               deallocate(ftime)
               if (nrec_ncdiag < ntime) then
                  write(nu_diag,*) subname,' WARNING: overwriting records from ', &
                                   nrec_ncdiag+1,' of ',ntime
               endif
            endif
         endif

      endif

      ncdiag_opened = .true.

      ! append the record, starting with the time
      nrec_ncdiag = nrec_ncdiag + 1

      status = nf90_inq_varid(ncid, 'time', varid)
      call ice_check_nc(status, subname//' ERROR: getting varid for time', &
                        file=__FILE__, line=__LINE__)
      status = nf90_put_var(ncid, varid, (/timesecs/real(seconds_per_day,kind=dbl_kind)/), &
                            start=(/nrec_ncdiag/), count=(/1/))
      call ice_check_nc(status, subname//' ERROR: writing time', &
                        file=__FILE__, line=__LINE__)

      if (diag_netcdf_avg) then
         status = nf90_inq_varid(ncid, 'time_bounds', varid)
         call ice_check_nc(status, subname//' ERROR: getting varid for time_bounds', &
                           file=__FILE__, line=__LINE__)
         status = nf90_put_var(ncid, varid, &
                               (/ncdiag_time_beg, timesecs/)/real(seconds_per_day,kind=dbl_kind), &
                               start=(/1,nrec_ncdiag/), count=(/2,1/))
         call ice_check_nc(status, subname//' ERROR: writing time_bounds', &
                           file=__FILE__, line=__LINE__)
      endif

      ! write the hemispheric diagnostics
      do n = 1, n_ncdiag
         do hemi = 1, 2
            vname = ncdiag_name(n,hemi)
            if (hemi == 1) then
               val = ncdiag_valn(n)
            else
               val = ncdiag_vals(n)
            endif
            if (.not. ncdiag_hemi(hemi)) val = spval_dbl
            status = nf90_inq_varid(ncid, trim(vname), varid)
            call ice_check_nc(status, subname//' ERROR: getting varid for '//trim(vname), &
                              file=__FILE__, line=__LINE__)
            status = nf90_put_var(ncid, varid, (/val/), start=(/nrec_ncdiag/), count=(/1/))
            call ice_check_nc(status, subname//' ERROR: writing '//trim(vname), &
                              file=__FILE__, line=__LINE__)
         enddo
      enddo

      status = nf90_close(ncid)
      call ice_check_nc(status, subname//' ERROR: closing '//trim(ncdiag_file), &
                        file=__FILE__, line=__LINE__)
#else
      call abort_ice(subname//' ERROR: USE_NETCDF cpp not defined', &
                     file=__FILE__, line=__LINE__)
#endif

      end subroutine write_ncdiag

!=======================================================================

      end module ice_diagnostics_netcdf

!=======================================================================
