!=======================================================================

! Hemispheric diagnostics written to a netcdf file at diag_netcdf_freq, as time
! averages over each interval (diag_netcdf_avg) or instantaneous values.
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
      use ice_diagnostics, only: hdiag, n_hdiag, hdiag_initialized, hdiag_hemi, &
          init_hdiag, hdiag_cell_values, hdiag_reduce, hdiag_derived, &
          hdiag_global, hdiag_global_values, print_global, &
          hd_ptot, hd_brine, hd_spond, hd_fisoa, hd_faeroa
      use ice_exit, only: abort_ice
      use ice_fileunits, only: nu_diag
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_parameters

      implicit none
      private
      public :: diag_netcdf_step

      logical (kind=log_kind), public :: &
         diag_netcdf_avg      ! if true, time averages, otherwise instantaneous values

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
         max_ncdiag = 256     ! max number of output variables per hemisphere

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
! from the drivers.  Accumulates time averages and writes a record at
! diag_netcdf_freq (write_diag_netcdf).

      subroutine diag_netcdf_step

      use ice_calendar, only: write_diag_netcdf, diag_netcdf_freq, write_histrest, write_restart

      character(len=*), parameter :: subname = '(diag_netcdf_step)'

      if (diag_netcdf_freq == 'x') return

      if (diag_netcdf_avg) call accum_ncdiag()

      if (write_diag_netcdf) then
         call compute_ncdiag()     ! global reductions on all tasks
         call write_ncdiag()       ! master task
      endif

      ! restart the averaging, with the history restarts
      if (diag_netcdf_avg .and. write_histrest .and. write_restart == 1) &
         call write_ncdiag_restart()

      end subroutine diag_netcdf_step

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

      if (.not. hdiag_initialized) call init_hdiag()

      if (.not. allocated(acc_ncdiag)) then
         allocate(acc_ncdiag(nx_block,ny_block,max_blocks,n_hdiag))
         call reset_ncdiag()
         ! continue the averaging from a restart, if there is one
         call read_ncdiag_restart()
      endif

      ! start of the averaging interval is the start of the first step
      if (nacc_steps == 0) acc_time_beg = timesecs - dt
      nacc_steps = nacc_steps + 1

      allocate(fld(nx_block,ny_block,max_blocks,n_hdiag))
      call hdiag_cell_values(fld)

      do k = 1, n_hdiag
         if (hdiag(k)%reduce == 'max') then
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

      do k = 1, n_hdiag
         if (hdiag(k)%reduce == 'max') then
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

      use ice_calendar, only: timesecs, dt
      use ice_domain_size, only: n_iso, n_aero
      use ice_domain_size, only: max_blocks

      ! local variables

      integer (kind=int_kind) :: &
         k

      real (kind=dbl_kind) :: &
         rn, urmsn, urmss, albn, albs, rdt

      type (hdiag_global_values) :: &
         gn, gs            ! print_global diagnostics, northern and southern hemisphere

      character (len=char_len) :: &
         cn                ! tracer index string

      integer (kind=int_kind) :: &
         n                 ! tracer index

      real (kind=dbl_kind), dimension(n_hdiag) :: &
         valn, vals        ! hemispheric values, before scale is applied

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable :: &
         fld               ! per grid cell values

      character (len=char_len) :: &
         cm_sum, cm_max    ! cell_methods for sums and maxima

      character(len=*), parameter :: subname = '(compute_ncdiag)'

      if (.not. hdiag_initialized) call init_hdiag()

      allocate(fld(nx_block,ny_block,max_blocks,n_hdiag))
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
         call hdiag_cell_values(fld)
         rn = c1
         cm_sum = ''
         cm_max = ''
      endif

      call hdiag_reduce(fld, rn, valn, vals)
      deallocate(fld)

      n_ncdiag = 0
      do k = 1, n_hdiag
         if (hdiag(k)%lwrite) then
            if (hdiag(k)%reduce == 'max') then
               call add_ncdiag(hdiag(k)%vname_n, hdiag(k)%vname_s, hdiag(k)%lname, &
                               hdiag(k)%units, valn(k)*hdiag(k)%scale, vals(k)*hdiag(k)%scale, &
                               hdiag(k)%sname, cm_max)
            else
               call add_ncdiag(hdiag(k)%vname_n, hdiag(k)%vname_s, hdiag(k)%lname, &
                               hdiag(k)%units, valn(k)*hdiag(k)%scale, vals(k)*hdiag(k)%scale, &
                               hdiag(k)%sname, cm_sum)
            endif
         endif
      enddo

      ! derived diagnostics

      ! rms ice speed, and ice area weighted albedo over sunlit points
      call hdiag_derived(valn, urmsn, albn)
      call hdiag_derived(vals, urmss, albs)
      call add_ncdiag('ice_speed_rms_nh', 'ice_speed_rms_sh', 'rms ice speed', 'm/s', &
                      urmsn, urmss, '', cm_sum)

      call add_ncdiag('albedo_nh', 'albedo_sh', 'average albedo', '1', &
                      albn, albs, '', cm_sum)

      ! print_global diagnostics, the same as the log.  Water, salt and tracer
      ! fluxes are rates.  The changes and errors are over the last timestep,
      ! so are only written as instantaneous values.
      if (print_global) then
         call hdiag_global(valn, 1, dt, gn)
         call hdiag_global(vals, 2, dt, gs)
         rdt = c1/dt

         ! mass, energy and salt
         call add_ncdiag('ice_mass_nh', 'ice_mass_sh', 'total ice mass', 'kg', &
                         gn%mice, gs%mice, '', cm_sum)
         if (hd_ptot > 0) &
         call add_ncdiag('pond_mass_nh', 'pond_mass_sh', 'total pond mass', 'kg', &
                         gn%mpnd, gs%mpnd, '', cm_sum)
         call add_ncdiag('total_mass_nh', 'total_mass_sh', 'total ice, snow and pond mass', 'kg', &
                         gn%mtot, gs%mtot, '', cm_sum)
         call add_ncdiag('total_energy_nh', 'total_energy_sh', 'total ice and snow energy', 'J', &
                         gn%etot, gs%etot, '', cm_sum)
         call add_ncdiag('salt_mass_nh', 'salt_mass_sh', 'total salt mass', 'kg', &
                         gn%mslt, gs%mslt, '', cm_sum)
         if (hd_brine > 0) &
         call add_ncdiag('brine_volume_nh', 'brine_volume_sh', 'total brine tracer volume', 'm^3', &
                         gn%shmaxt, gs%shmaxt, '', cm_sum)

         ! water and salt fluxes
         call add_ncdiag('rain_flux_nh', 'rain_flux_sh', 'rain on sea ice', 'kg s-1', &
                         gn%rn*rdt, gs%rn*rdt, '', cm_sum)
         call add_ncdiag('snow_flux_nh', 'snow_flux_sh', 'snow on sea ice', 'kg s-1', &
                         gn%sn*rdt, gs%sn*rdt, '', cm_sum)
         call add_ncdiag('evap_flux_nh', 'evap_flux_sh', 'evaporation and sublimation', 'kg s-1', &
                         gn%evp*rdt, gs%evp*rdt, '', cm_sum)
         call add_ncdiag('frazil_flux_nh', 'frazil_flux_sh', 'frazil ice formation', 'kg s-1', &
                         gn%frz*rdt, gs%frz*rdt, '', cm_sum)
         if (hd_spond > 0) &
         call add_ncdiag('pond_flux_nh', 'pond_flux_sh', 'pond water flux', 'kg s-1', &
                         gn%spond*rdt, gs%spond*rdt, '', cm_sum)
         call add_ncdiag('fresh_flux_nh', 'fresh_flux_sh', 'fresh water flux to ocean', 'kg s-1', &
                         gn%sfresh*rdt, gs%sfresh*rdt, '', cm_sum)
         call add_ncdiag('water_flux_nh', 'water_flux_sh', 'net water flux into sea ice', 'kg s-1', &
                         gn%flux*rdt, gs%flux*rdt, '', cm_sum)
         call add_ncdiag('salt_flux_nh', 'salt_flux_sh', 'salt flux to ocean', 'kg s-1', &
                         gn%sfsalt*rdt, gs%sfsalt*rdt, '', cm_sum)

         ! heat fluxes
         call add_ncdiag('atm_heat_flux_nh', 'atm_heat_flux_sh', 'atmosphere heat flux into sea ice', 'W', &
                         gn%fhatm, gs%fhatm, '', cm_sum)
         call add_ncdiag('ocn_heat_flux_nh', 'ocn_heat_flux_sh', 'ocean heat flux used by sea ice', 'W', &
                         gn%fhocn, gs%fhocn, '', cm_sum)
         call add_ncdiag('frazil_heat_flux_nh', 'frazil_heat_flux_sh', 'frazil heat flux', 'W', &
                         gn%fhfrz, gs%fhfrz, '', cm_sum)
         call add_ncdiag('net_heat_flux_nh', 'net_heat_flux_sh', 'net heat flux into sea ice', 'W', &
                         gn%hnet*rdt, gs%hnet*rdt, '', cm_sum)
         call add_ncdiag('sw_down_flux_nh', 'sw_down_flux_sh', 'incoming shortwave on sea ice', 'W', &
                         gn%fswdn, gs%fswdn, '', cm_sum)
         call add_ncdiag('sw_absorbed_flux_nh', 'sw_absorbed_flux_sh', 'absorbed shortwave', 'W', &
                         gn%fswnet, gs%fswnet, '', cm_sum)

         ! isotopes and aerosols
         do n = 1, n_iso
            if (hd_fisoa(n) == 0) cycle
            write(cn,'(i0)') n
            call add_ncdiag('isotope'//trim(cn)//'_atm_flux_nh', 'isotope'//trim(cn)//'_atm_flux_sh', &
                            'isotope flux from atmosphere', 'kg s-1', &
                            gn%fisoa(n)*rdt, gs%fisoa(n)*rdt, '', cm_sum)
            call add_ncdiag('isotope'//trim(cn)//'_ocn_flux_nh', 'isotope'//trim(cn)//'_ocn_flux_sh', &
                            'isotope flux to ocean', 'kg s-1', &
                            gn%fisoo(n)*rdt, gs%fisoo(n)*rdt, '', cm_sum)
            call add_ncdiag('isotope'//trim(cn)//'_total_nh', 'isotope'//trim(cn)//'_total_sh', &
                            'total isotope', 'kg', gn%isotot(n), gs%isotot(n), '', cm_sum)
            call add_ncdiag('isotope'//trim(cn)//'_max_nh', 'isotope'//trim(cn)//'_max_sh', &
                            'maximum isotope per unit area', 'kg m-2', gn%isomx(n), gs%isomx(n), '', cm_max)
         enddo
         do n = 1, n_aero
            if (hd_faeroa(n) == 0) cycle
            write(cn,'(i0)') n
            call add_ncdiag('aerosol'//trim(cn)//'_atm_flux_nh', 'aerosol'//trim(cn)//'_atm_flux_sh', &
                            'aerosol flux from atmosphere', 'kg s-1', &
                            gn%faeroa(n)*rdt, gs%faeroa(n)*rdt, '', cm_sum)
            call add_ncdiag('aerosol'//trim(cn)//'_ocn_flux_nh', 'aerosol'//trim(cn)//'_ocn_flux_sh', &
                            'aerosol flux to ocean', 'kg s-1', &
                            gn%faeroo(n)*rdt, gs%faeroo(n)*rdt, '', cm_sum)
            call add_ncdiag('aerosol'//trim(cn)//'_total_nh', 'aerosol'//trim(cn)//'_total_sh', &
                            'total aerosol', 'kg', gn%aerotot(n), gs%aerotot(n), '', cm_sum)
            call add_ncdiag('aerosol'//trim(cn)//'_max_nh', 'aerosol'//trim(cn)//'_max_sh', &
                            'maximum aerosol per unit area', 'kg m-2', gn%aeromx(n), gs%aeromx(n), '', cm_max)
         enddo

         ! changes and errors over the last timestep, instantaneous only
         if (.not. diag_netcdf_avg) then
            call add_ncdiag('total_mass_change_nh', 'total_mass_change_sh', &
                            'change in total mass over the timestep', 'kg', gn%delmi, gs%delmi, '', '')
            call add_ncdiag('water_flux_err_nh', 'water_flux_err_sh', &
                            'water conservation error over the timestep, fraction of total mass', '1', &
                            gn%werr, gs%werr, '', '')
            call add_ncdiag('total_energy_change_nh', 'total_energy_change_sh', &
                            'change in total energy over the timestep', 'J', gn%delei, gs%delei, '', '')
            call add_ncdiag('heat_err_nh', 'heat_err_sh', &
                            'heat conservation error over the timestep, fraction of total energy', '1', &
                            gn%herr, gs%herr, '', '')
            call add_ncdiag('sw_err_nh', 'sw_err_sh', &
                            'shortwave error, fraction of absorbed shortwave', '1', &
                            gn%swerr, gs%swerr, '', '')
            call add_ncdiag('salt_mass_change_nh', 'salt_mass_change_sh', &
                            'change in salt mass over the timestep', 'kg', gn%delmslt, gs%delmslt, '', '')
            call add_ncdiag('salt_flux_err_nh', 'salt_flux_err_sh', &
                            'salt conservation error over the timestep, fraction of salt mass', '1', &
                            gn%serr, gs%serr, '', '')
            do n = 1, n_iso
               if (hd_fisoa(n) == 0) cycle
               write(cn,'(i0)') n
               call add_ncdiag('isotope'//trim(cn)//'_err_nh', 'isotope'//trim(cn)//'_err_sh', &
                               'isotope conservation error over the timestep, fraction of total', '1', &
                               gn%isor(n), gs%isor(n), '', '')
            enddo
            do n = 1, n_aero
               if (hd_faeroa(n) == 0) cycle
               write(cn,'(i0)') n
               call add_ncdiag('aerosol'//trim(cn)//'_err_nh', 'aerosol'//trim(cn)//'_err_sh', &
                               'aerosol conservation error over the timestep, fraction of total', '1', &
                               gn%aerr(n), gs%aerr(n), '', '')
            enddo
         endif
      endif

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

      ! filename, _avg for time averages or _inst for instantaneous values,
      ! as for history files
      if (diag_netcdf_avg) then
         ncdiag_file = trim(history_dir)//trim(history_file)//'_diag_avg.nc'
      else
         ncdiag_file = trim(history_dir)//trim(history_file)//'_diag_inst.nc'
      endif
      lexist = .true.
      if (.not. ncdiag_opened) inquire(file=trim(ncdiag_file), exist=lexist)

      ! create the file if it does not exist, otherwise open it to append a record
      if (.not. lexist) then

         write(nu_diag,*) subname,' creating ',trim(ncdiag_file)

         ! file format consistent with history files, written serially
         iflag = ncdiag_create_mode()

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
            if (.not. hdiag_hemi(hemi)) val = spval_dbl
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

! netcdf create mode for the diagnostics files, consistent with history_format

      integer function ncdiag_create_mode()

#ifdef USE_NETCDF
      use netcdf
#endif
      use ice_history_shared, only: history_format

      character(len=*), parameter :: subname = '(ncdiag_create_mode)'

      ncdiag_create_mode = 0
#ifdef USE_NETCDF
      if (history_format == 'cdf2' .or. history_format == 'pnetcdf2') then
#ifdef NO_CDF2
         call abort_ice(subname//' ERROR: history_format cdf2 not available ', &
                        file=__FILE__, line=__LINE__)
#else
         ncdiag_create_mode = ior(nf90_clobber,nf90_64bit_offset)
#endif
      elseif (history_format == 'cdf5' .or. history_format == 'pnetcdf5') then
#ifdef NO_CDF5
         call abort_ice(subname//' ERROR: history_format cdf5 not available ', &
                        file=__FILE__, line=__LINE__)
#else
         ncdiag_create_mode = ior(nf90_clobber,nf90_64bit_data)
#endif
      elseif (history_format == 'hdf5') then
#ifdef NO_HDF5
         call abort_ice(subname//' ERROR: history_format hdf5 not available ', &
                        file=__FILE__, line=__LINE__)
#else
         ncdiag_create_mode = ior(nf90_clobber,nf90_netcdf4)
#endif
      else
         ncdiag_create_mode = nf90_clobber
      endif
#endif

      end function ncdiag_create_mode

!=======================================================================

! Restart file name for the time averaged netcdf diagnostics at a date,
! in restart_dir like the history restarts

      subroutine ncdiag_restart_name(fname, ayear, amon, aday, asec)

      use ice_history_shared, only: history_file
      use ice_restart_shared, only: restart_dir

      character (len=char_len_long), intent(out) :: &
         fname           ! restart file name

      integer (kind=int_kind), intent(in) :: &
         ayear, amon, aday, asec   ! date

      write(fname,'(a,i4.4,a,i2.2,a,i2.2,a,i5.5,a)') &
         trim(restart_dir)//trim(history_file)//'_diag_r.', &
         ayear,'-',amon,'-',aday,'-',asec,'.nc'

      end subroutine ncdiag_restart_name

!=======================================================================

! Write the accumulation for time averages to a restart file, so the
! averaging continues across restarts.  Written with the history restarts.

      subroutine write_ncdiag_restart

#ifdef USE_NETCDF
      use netcdf
#endif
      use ice_calendar, only: myear, mmonth, mday, msec
      use ice_domain, only: distrb_info
      use ice_domain_size, only: nx_global, ny_global
      use ice_gather_scatter, only: gather_global
      use ice_read_write, only: ice_check_nc

      ! local variables

      integer (kind=int_kind) :: &
         k, status, ncid, dimid(2), varid

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g          ! global field

      character (len=char_len_long) :: &
         fname           ! restart file name

      character (len=char_len) :: &
         vname           ! variable name

      character(len=*), parameter :: subname = '(write_ncdiag_restart)'

#ifdef USE_NETCDF
      if (.not. allocated(acc_ncdiag)) return

      if (my_task == master_task) then
         allocate(work_g(nx_global,ny_global))
         call ncdiag_restart_name(fname, myear, mmonth, mday, msec)
         write(nu_diag,*) subname,' writing ',trim(fname)
         status = nf90_create(trim(fname), ncdiag_create_mode(), ncid)
         call ice_check_nc(status, subname//' ERROR: creating '//trim(fname), &
                           file=__FILE__, line=__LINE__)
         status = nf90_def_dim(ncid, 'ni', nx_global, dimid(1))
         call ice_check_nc(status, subname//' ERROR: defining dim ni', &
                           file=__FILE__, line=__LINE__)
         status = nf90_def_dim(ncid, 'nj', ny_global, dimid(2))
         call ice_check_nc(status, subname//' ERROR: defining dim nj', &
                           file=__FILE__, line=__LINE__)
         status = nf90_def_var(ncid, 'nacc_steps', nf90_int, varid)
         call ice_check_nc(status, subname//' ERROR: defining var nacc_steps', &
                           file=__FILE__, line=__LINE__)
         status = nf90_def_var(ncid, 'acc_time_beg', nf90_double, varid)
         call ice_check_nc(status, subname//' ERROR: defining var acc_time_beg', &
                           file=__FILE__, line=__LINE__)
         do k = 1, n_hdiag
            write(vname,'(a,i3.3)') 'acc_ncdiag_',k
            status = nf90_def_var(ncid, trim(vname), nf90_double, dimid, varid)
            call ice_check_nc(status, subname//' ERROR: defining var '//trim(vname), &
                              file=__FILE__, line=__LINE__)
            status = nf90_put_att(ncid, varid, 'long_name', 'accumulated '//trim(hdiag(k)%lname))
            call ice_check_nc(status, subname//' ERROR: defining att long_name '//trim(vname), &
                              file=__FILE__, line=__LINE__)
         enddo
         status = nf90_enddef(ncid)
         call ice_check_nc(status, subname//' ERROR: in nf90_enddef', &
                           file=__FILE__, line=__LINE__)
         status = nf90_inq_varid(ncid, 'nacc_steps', varid)
         status = nf90_put_var(ncid, varid, nacc_steps)
         call ice_check_nc(status, subname//' ERROR: writing nacc_steps', &
                           file=__FILE__, line=__LINE__)
         status = nf90_inq_varid(ncid, 'acc_time_beg', varid)
         status = nf90_put_var(ncid, varid, acc_time_beg)
         call ice_check_nc(status, subname//' ERROR: writing acc_time_beg', &
                           file=__FILE__, line=__LINE__)
      else
         allocate(work_g(1,1))
      endif

      do k = 1, n_hdiag
         call gather_global(work_g, acc_ncdiag(:,:,:,k), master_task, distrb_info)
         if (my_task == master_task) then
            write(vname,'(a,i3.3)') 'acc_ncdiag_',k
            status = nf90_inq_varid(ncid, trim(vname), varid)
            call ice_check_nc(status, subname//' ERROR: getting varid for '//trim(vname), &
                              file=__FILE__, line=__LINE__)
            status = nf90_put_var(ncid, varid, work_g)
            call ice_check_nc(status, subname//' ERROR: writing '//trim(vname), &
                              file=__FILE__, line=__LINE__)
         endif
      enddo

      if (my_task == master_task) then
         status = nf90_close(ncid)
         call ice_check_nc(status, subname//' ERROR: closing '//trim(fname), &
                           file=__FILE__, line=__LINE__)
      endif
      deallocate(work_g)
#endif

      end subroutine write_ncdiag_restart

!=======================================================================

! Read the accumulation for time averages from a restart file at the start
! of the run, if it exists, like the history restarts

      subroutine read_ncdiag_restart

#ifdef USE_NETCDF
      use netcdf
#endif
      use ice_broadcast, only: broadcast_scalar
      use ice_calendar, only: timesecs, dt, calendar_time2date
      use ice_constants, only: field_loc_noupdate, field_type_noupdate
      use ice_domain, only: distrb_info
      use ice_domain_size, only: nx_global, ny_global
      use ice_gather_scatter, only: scatter_global
      use ice_read_write, only: ice_check_nc

      ! local variables

      integer (kind=int_kind) :: &
         k, status, ncid, varid, &
         syear, smon, sday, ssec  ! run start date

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g          ! global field

      character (len=char_len_long) :: &
         fname           ! restart file name

      character (len=char_len) :: &
         vname           ! variable name

      character(len=*), parameter :: subname = '(read_ncdiag_restart)'

#ifdef USE_NETCDF
      if (my_task == master_task) then
         ! the restart is from the start of the first step of this run
         call calendar_time2date(timesecs - dt, syear, smon, sday, ssec)
         call ncdiag_restart_name(fname, syear, smon, sday, ssec)
         status = nf90_open(trim(fname), nf90_nowrite, ncid)
         if (status /= nf90_noerr) then
            write(nu_diag,*) subname,' file not found ',trim(fname)
         else
            write(nu_diag,*) subname,' reading ',trim(fname)
         endif
      endif
      call broadcast_scalar(status, master_task)
      if (status /= nf90_noerr) return

      if (my_task == master_task) then
         allocate(work_g(nx_global,ny_global))
         status = nf90_inq_varid(ncid, 'nacc_steps', varid)
         if (status == nf90_noerr) status = nf90_get_var(ncid, varid, nacc_steps)
         call ice_check_nc(status, subname//' ERROR: reading nacc_steps', &
                           file=__FILE__, line=__LINE__)
         status = nf90_inq_varid(ncid, 'acc_time_beg', varid)
         if (status == nf90_noerr) status = nf90_get_var(ncid, varid, acc_time_beg)
         call ice_check_nc(status, subname//' ERROR: reading acc_time_beg', &
                           file=__FILE__, line=__LINE__)
      else
         allocate(work_g(1,1))
      endif
      call broadcast_scalar(nacc_steps, master_task)
      call broadcast_scalar(acc_time_beg, master_task)

      do k = 1, n_hdiag
         if (my_task == master_task) then
            write(vname,'(a,i3.3)') 'acc_ncdiag_',k
            status = nf90_inq_varid(ncid, trim(vname), varid)
            call ice_check_nc(status, subname//' ERROR: getting varid for '//trim(vname), &
                              file=__FILE__, line=__LINE__)
            status = nf90_get_var(ncid, varid, work_g)
            call ice_check_nc(status, subname//' ERROR: reading '//trim(vname), &
                              file=__FILE__, line=__LINE__)
         endif
         call scatter_global(acc_ncdiag(:,:,:,k), work_g, master_task, distrb_info, &
                             field_loc_noupdate, field_type_noupdate)
      enddo

      if (my_task == master_task) then
         status = nf90_close(ncid)
         call ice_check_nc(status, subname//' ERROR: closing '//trim(fname), &
                           file=__FILE__, line=__LINE__)
      endif
      deallocate(work_g)
#endif

      end subroutine read_ncdiag_restart

!=======================================================================

      end module ice_diagnostics_netcdf

!=======================================================================
