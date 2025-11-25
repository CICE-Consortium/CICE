#ifdef ncdf
#define USE_NETCDF
#endif
!=======================================================================
!
! Reads and interpolates forcing data for atmosphere and ocean quantities.
!
! authors: Elizabeth C. Hunke and William H. Lipscomb, LANL
!
! 2004 WHL: Block structure added
! 2005 WHL: ECMWF option added
! 2006 ECH: LY option added
! 2006 WHL: Module name changed from ice_flux_in
! 2006 ECH: Fixed bugs, rearranged routines, edited comments, etc.
!           Added NCAR ocean forcing file
!           Converted to free source form (F90)
! 2007: netcdf version of read_data added by Alison McLaren, Met Office
!
      module ice_forcing

      use ice_kinds_mod
      use ice_boundary, only: ice_HaloUpdate
      use ice_blocks, only: nx_block, ny_block
      use ice_domain, only: halo_info
      use ice_domain_size, only: ncat, max_blocks, nx_global, ny_global, nfreq
      use ice_communicate, only: my_task, master_task
      use ice_calendar, only: istep, istep1, &
                              msec, mday, mmonth, myear, yday, daycal, &
                              daymo, days_per_year, compute_days_between
      use ice_fileunits, only: nu_diag, nu_forcing
      use ice_exit, only: abort_ice
      use ice_read_write, only: ice_open, ice_read, ice_check_nc, &
                                ice_get_ncvarsize, ice_read_vec_nc, &
                                ice_open_nc, ice_read_nc, ice_close_nc
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_readwrite, &
                            timer_bound, timer_forcing
      use ice_arrays_column, only: oceanmixed_ice, restore_bgc
      use ice_constants, only: c0, c1, c2, c3, c4, c5, c8, c10, c12, c15, c20, &
                               c180, c360, c365, c1000, c3600
      use ice_constants, only: p001, p01, p1, p2, p25, p5, p6
      use ice_constants, only: cm_to_m
      use ice_constants, only: field_loc_center, field_type_scalar, &
                               field_type_vector, field_loc_NEcorner
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_sea_freezing_temperature
      use icepack_intfc, only: icepack_init_wave, icepack_init_parameters
      use icepack_intfc, only: icepack_query_tracer_indices, icepack_query_parameters

      implicit none
      private
      public :: init_forcing_atmo, init_forcing_ocn, alloc_forcing, &
                get_forcing_atmo, get_forcing_ocn, get_wave_spec, &
                read_clim_data, read_clim_data_nc, &
                interpolate_data, interp_coeff_monthly, &
                read_data_nc_point, interp_coeff, &
                init_snowtable

      integer (kind=int_kind), public :: &
         ycycle          , & ! number of years in forcing cycle, set by namelist
         fyear_init      , & ! first year of data in forcing cycle, set by namelist
         fyear           , & ! current year in forcing cycle, varying during the run
         fyear_final         ! last year in cycle, computed at init

      character (char_len_long) :: &        ! input data file names
          uwind_file, &  ! this is also used a generic file containing all fields for JRA55
          vwind_file, &
           wind_file, &
          strax_file, &
          stray_file, &
           tair_file, &
          humid_file, &
           rhoa_file, &
            fsw_file, &
            flw_file, &
           rain_file, &
            sst_file, &
            sss_file, &
         sublim_file, &
           snow_file

      character (char_len_long), dimension(:), allocatable, public :: &  ! input data file names
        topmelt_file, &
        botmelt_file

      real (kind=dbl_kind), public  :: &
           c1intp, c2intp     ! interpolation coefficients

      integer (kind=int_kind) :: &
           oldrecnum = 0  , & ! old record number (save between steps)
           oldrecnum4X = 0    !

      real (kind=dbl_kind), dimension(:,:,:), allocatable, public :: &
          cldf                ! cloud fraction

      real (kind=dbl_kind), dimension(:,:,:,:), allocatable, public :: &
            fsw_data, & ! field values at 2 temporal data points
           cldf_data, &
          fsnow_data, &
           Tair_data, &
           uatm_data, &
           vatm_data, &
           wind_data, &
          strax_data, &
          stray_data, &
             Qa_data, &
           rhoa_data, &
            flw_data, &
            sst_data, &
            sss_data, &
           uocn_data, &
           vocn_data, &
         sublim_data, &
          frain_data

      real (kind=dbl_kind), dimension(:,:,:,:,:), allocatable, public :: &
        topmelt_data, &
        botmelt_data

      real (kind=dbl_kind), dimension(:,:,:,:,:), allocatable :: &
          wave_spectrum_data ! field values at 2 temporal data points

      character(char_len), public :: &
         atm_data_format   , & ! 'bin'=binary or 'nc'=netcdf
         ocn_data_format   , & ! 'bin'=binary or 'nc'=netcdf
         atm_data_type     , & ! 'default', 'monthly', 'ncar', 'box2001'
                               ! 'hadgem', 'oned', 'calm', 'uniform'
                               ! 'JRA55' or 'JRA55do'
         atm_data_version  , & ! date of atm_forcing file creation
         bgc_data_type     , & ! 'default', 'clim'
         ocn_data_type     , & ! 'default', 'clim', 'ncar', 'oned', 'calm', 'box2001'
                               ! 'hadgem_sst' or 'hadgem_sst_uvocn', 'uniform'
         ice_data_type     , & ! 'latsst', 'box2001', 'boxslotcyl', etc
         ice_data_conc     , & ! 'p5','p8','p9','c1','parabolic', 'box2001', etc
         ice_data_dist     , & ! 'box2001','gauss', 'uniform', etc
         precip_units          ! 'mm_per_month', 'mm_per_sec', 'mks','m_per_sec'

      logical (kind=log_kind), public :: &
         rotate_wind      ! rotate wind/stress to computational grid from true north directed

      character(char_len_long), public :: &
         atm_data_dir , & ! top directory for atmospheric data
         ocn_data_dir , & ! top directory for ocean data
         wave_spec_dir, & ! dir name for wave spectrum
         wave_spec_file,& ! file name for wave spectrum
         oceanmixed_file  ! file name for ocean forcing data

      integer (kind=int_kind), parameter :: &
         nfld = 8   ! number of fields to search for in forcing file

      ! as in the dummy atm (latm)
      real (kind=dbl_kind), parameter, public :: &
         frcvdr = 0.28_dbl_kind, & ! frac of incoming sw in vis direct band
         frcvdf = 0.24_dbl_kind, & ! frac of incoming sw in vis diffuse band
         frcidr = 0.31_dbl_kind, & ! frac of incoming sw in near IR direct band
         frcidf = 0.17_dbl_kind    ! frac of incoming sw in near IR diffuse band

      real (kind=dbl_kind), dimension (:,:,:,:,:), allocatable, public :: &
         ocn_frc_m   ! ocn data for 12 months

      logical (kind=log_kind), public :: &
         restore_ocn                 ! restore sst if true

      integer (kind=int_kind), public :: &
         trestore                    ! restoring time scale (days)

      real (kind=dbl_kind), public :: &
         trest                       ! restoring time scale (sec)

      logical (kind=log_kind), public :: &
         debug_forcing               ! prints forcing debugging output if true

      real (dbl_kind), dimension(:), allocatable, public :: &
         jday_atm  ! jday time vector from atm forcing files

      integer (kind=int_kind), public :: &
         Njday_atm ! Number of atm forcing timesteps

      character (len=char_len_long), public :: &
         snw_filename        ! filename for snow lookup table

      character (char_len), public :: &
         snw_rhos_fname , &  ! snow table 1d rhos field name
         snw_Tgrd_fname , &  ! snow table 1d Tgrd field name
         snw_T_fname    , &  ! snow table 1d T field name
         snw_tau_fname  , &  ! snow table 3d tau field name
         snw_kappa_fname, &  ! snow table 3d kappa field name
         snw_drdt0_fname     ! snow table 3d drdt0 field name

      ! PRIVATE:

      real (dbl_kind), parameter :: &
         mixed_layer_depth_default = c20  ! default mixed layer depth in m

      logical (kind=log_kind), parameter :: &
         local_debug = .false.   ! local debug flag

!=======================================================================

      contains

!=======================================================================
!
! Allocate space for all variables
!
      subroutine alloc_forcing
      integer (int_kind) :: ierr
      character(len=*), parameter :: subname = '(alloc_forcing)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      allocate ( &
         cldf        (nx_block,ny_block,  max_blocks), & ! cloud fraction
         fsw_data    (nx_block,ny_block,2,max_blocks), & ! field values at 2 temporal data points
         cldf_data   (nx_block,ny_block,2,max_blocks), &
         fsnow_data  (nx_block,ny_block,2,max_blocks), &
         Tair_data   (nx_block,ny_block,2,max_blocks), &
         uatm_data   (nx_block,ny_block,2,max_blocks), &
         vatm_data   (nx_block,ny_block,2,max_blocks), &
         wind_data   (nx_block,ny_block,2,max_blocks), &
         strax_data  (nx_block,ny_block,2,max_blocks), &
         stray_data  (nx_block,ny_block,2,max_blocks), &
         Qa_data     (nx_block,ny_block,2,max_blocks), &
         rhoa_data   (nx_block,ny_block,2,max_blocks), &
         flw_data    (nx_block,ny_block,2,max_blocks), &
         sst_data    (nx_block,ny_block,2,max_blocks), &
         sss_data    (nx_block,ny_block,2,max_blocks), &
         uocn_data   (nx_block,ny_block,2,max_blocks), &
         vocn_data   (nx_block,ny_block,2,max_blocks), &
         sublim_data (nx_block,ny_block,2,max_blocks), &
         frain_data  (nx_block,ny_block,2,max_blocks), &
         topmelt_data(nx_block,ny_block,2,max_blocks,ncat), &
         botmelt_data(nx_block,ny_block,2,max_blocks,ncat), &
         ocn_frc_m   (nx_block,ny_block,  max_blocks,nfld,12), & ! ocn data for 12 months
         topmelt_file(ncat), &
         botmelt_file(ncat), &
         wave_spectrum_data(nx_block,ny_block,nfreq,2,max_blocks), &
         stat=ierr)
      if (ierr/=0) call abort_ice('(alloc_forcing): Out of Memory')

      cldf         = c0
      fsw_data     = c0
      cldf_data    = c0
      fsnow_data   = c0
      Tair_data    = c0
      uatm_data    = c0
      vatm_data    = c0
      wind_data    = c0
      strax_data   = c0
      stray_data   = c0
      Qa_data      = c0
      rhoa_data    = c0
      flw_data     = c0
      sst_data     = c0
      sss_data     = c0
      uocn_data    = c0
      vocn_data    = c0
      sublim_data  = c0
      frain_data   = c0
      topmelt_data = c0
      botmelt_data = c0
      ocn_frc_m    = c0
      topmelt_file = ''
      botmelt_file = ''
      wave_spectrum_data = c0

      end subroutine alloc_forcing

!=======================================================================

      subroutine init_forcing_atmo

! Determine the current and final year of the forcing cycle based on
! namelist input; initialize the atmospheric forcing data filenames.

      use ice_calendar, only: use_leap_years

      integer (kind=int_kind) :: modadj   ! adjustment for mod function
      character(len=*), parameter :: subname = '(init_forcing_atmo)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      modadj      = abs((min(0,myear-fyear_init)/ycycle+1)*ycycle)
      fyear       = fyear_init + mod(myear-fyear_init+modadj,ycycle)
      fyear_final = fyear_init + ycycle - 1 ! last year in forcing cycle

      if (local_debug .and. my_task == master_task) then
         write(nu_diag,*) subname,'fdbg fyear = ',fyear,fyear_init,fyear_final
         write(nu_diag,*) subname,'fdbg atm_data_type = ',trim(atm_data_type)
      endif

      if (trim(atm_data_type) /= 'default' .and. &
                          my_task == master_task) then
         write (nu_diag,*) ' Initial forcing data year = ',fyear_init
         write (nu_diag,*) ' Final   forcing data year = ',fyear_final
      endif

      if (trim(atm_data_type) == 'hadgem' .and. &
          trim(precip_units) /= 'mks') then
         if (my_task == master_task) then
            write (nu_diag,*) 'WARNING: HadGEM atmospheric data chosen with wrong precip_units'
            write (nu_diag,*) 'WARNING:   Changing precip_units to mks (i.e. kg/m2 s).'
         endif
         call abort_ice(error_message=subname//' HadGEM precip_units error', &
            file=__FILE__, line=__LINE__)
      endif

      if (use_leap_years .and. (index(trim(atm_data_type),'JRA55') == 0 .and. &
                                trim(atm_data_type) /= 'hycom'       .and. &
                                trim(atm_data_type) /= 'box2001'))   then
         write(nu_diag,*) 'use_leap_years option is currently only supported for'
         write(nu_diag,*) 'JRA55, JRA55do, default , and box2001 atmospheric data'
         call abort_ice(error_message=subname, file=__FILE__, line=__LINE__)
      endif

    !-------------------------------------------------------------------
    ! Get filenames for input forcing data
    !-------------------------------------------------------------------

      ! default forcing values from init_flux_atm
      if (trim(atm_data_type) == 'ncar') then
         call NCAR_files(fyear)
      elseif (index(trim(atm_data_type),'JRA55') > 0) then
         call JRA55_files(fyear)
      elseif (trim(atm_data_type) == 'hadgem') then
         call hadgem_files(fyear)
      elseif (trim(atm_data_type) == 'monthly') then
         call monthly_files(fyear)
      elseif (trim(atm_data_type) == 'oned') then
         call oned_files
      elseif (trim(atm_data_type) == 'ISPOL') then
         call ISPOL_files
      elseif (trim(atm_data_type) == 'box2001') then
         call box2001_data_atm
      elseif (trim(atm_data_type) == 'uniform_northeast') then
         call uniform_data_atm('NE')
      elseif (trim(atm_data_type) == 'uniform_north') then
         call uniform_data_atm('N')
      elseif (trim(atm_data_type) == 'uniform_east') then
         call uniform_data_atm('E')
      elseif (trim(atm_data_type) == 'uniform_south') then
         call uniform_data_atm('S')
      elseif (trim(atm_data_type) == 'uniform_west') then
         call uniform_data_atm('W')
      elseif (trim(atm_data_type) == 'calm') then
         call uniform_data_atm('N',c0) ! direction does not matter when c0
      elseif (trim(atm_data_type) == 'hycom') then
         call hycom_atm_files
      elseif (trim(atm_data_type) == 'default') then
         ! don't need to do anything more
      else
        call abort_ice (error_message=subname//' ERROR atm_data_type unknown = '// &
                        trim(atm_data_type), file=__FILE__, line=__LINE__)
      endif

      end subroutine init_forcing_atmo

!=======================================================================

      subroutine init_forcing_ocn(dt)

! Set sea surface salinity and freezing temperature to annual mean value
!  using a 12-month climatology.
! Read sst data for current month, and adjust sst based on freezing
! temperature.  No interpolation in time.

! Note: SST is subsequently prognosed if CICE is run
! with a mixed layer ocean (oceanmixed_ice = T), and can be
! restored to data (restore_ocn = T).

      use ice_blocks, only: nx_block, ny_block
      use ice_domain, only: nblocks
      use ice_domain_size, only: max_blocks
      use ice_flux, only: sss, sst, Tf

      real (kind=dbl_kind), intent(in) :: &
         dt                   ! time step

      ! local variables

      integer (kind=int_kind) :: &
         i, j, iblk       , & ! horizontal indices
         k                , & ! month index
         fid              , & ! file id for netCDF file
         nbits

      logical (kind=log_kind) :: diag

      real (kind=dbl_kind) :: secday

      character (char_len) :: &
         fieldname            ! field name in netcdf file

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1

      character(len=*), parameter :: subname = '(init_forcing_ocn)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      call icepack_query_parameters(secday_out=secday)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      call alloc_forcing()

      sst_data(:,:,:,:) = c0
      sss_data(:,:,:,:) = c0
      uocn_data(:,:,:,:) = c0
      vocn_data(:,:,:,:) = c0

      nbits = 64              ! double precision data

      if (restore_ocn .or. restore_bgc) then
         if (trestore == 0) then
            trest = dt        ! use data instantaneously
         else
            trest = real(trestore,kind=dbl_kind) * secday ! seconds
         endif
      endif

    !-------------------------------------------------------------------
    ! Sea surface salinity (SSS)
    ! initialize to annual climatology created from monthly data
    !-------------------------------------------------------------------

      if (trim(ocn_data_type) == 'clim') then

         sss_file = trim(ocn_data_dir)//'/sss.mm.100x116.da' ! gx3 only

         if (my_task == master_task) then
            write (nu_diag,*) ' '
            write (nu_diag,*) 'SSS climatology computed from:'
            write (nu_diag,*) trim(sss_file)
         endif

         if (my_task == master_task) &
              call ice_open (nu_forcing, sss_file, nbits)

         sss(:,:,:) = c0

         do k = 1,12            ! loop over 12 months
            call ice_read (nu_forcing, k, work1, 'rda8', debug_forcing, &
                           field_loc_center, field_type_scalar)
            !$OMP PARALLEL DO PRIVATE(iblk,i,j)
            do iblk = 1, nblocks
               do j = 1, ny_block
               do i = 1, nx_block
                  sss(i,j,iblk) = sss(i,j,iblk) + work1(i,j,iblk)
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO
         enddo                  ! k

         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               sss(i,j,iblk) = sss(i,j,iblk) / c12   ! annual average
               sss(i,j,iblk) = max(sss(i,j,iblk),c0)
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO

         call ocn_freezing_temperature

         if (my_task == master_task) close(nu_forcing)

    !-------------------------------------------------------------------
    ! Sea surface temperature (SST)
    ! initialize to data for current month
    !-------------------------------------------------------------------

         if (nx_global == 320) then ! gx1
            sst_file = trim(ocn_data_dir)//'/sst_clim_hurrell.dat'
         else                   ! gx3
            sst_file = trim(ocn_data_dir)//'/sst.mm.100x116.da'
         endif

         if (my_task == master_task) then
            write (nu_diag,*) ' '
            write (nu_diag,*) 'Initial SST file:', trim(sst_file)
         endif

         if (my_task == master_task) &
              call ice_open (nu_forcing, sst_file, nbits)

         call ice_read (nu_forcing, mmonth, sst, 'rda8', debug_forcing, &
                        field_loc_center, field_type_scalar)

         if (my_task == master_task) close(nu_forcing)

         ! Make sure sst is not less than freezing temperature Tf
         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               sst(i,j,iblk) = max(sst(i,j,iblk),Tf(i,j,iblk))
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO

      elseif (trim(ocn_data_type) == 'hadgem_sst' .or.  &
          trim(ocn_data_type) == 'hadgem_sst_uvocn') then

         diag = .true.   ! write diagnostic information

         sst_file = trim (ocn_data_dir)//'/MONTHLY/sst.1997.nc'

         if (my_task == master_task) then

             write (nu_diag,*) ' '
             write (nu_diag,*) 'Initial SST file:', trim(sst_file)

             call ice_open_nc(sst_file,fid)

         endif

         fieldname='sst'
         call ice_read_nc(fid,mmonth,fieldname,sst,diag)

         if (my_task == master_task) call ice_close_nc(fid)

         ! Make sure sst is not less than freezing temperature Tf
         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               sst(i,j,iblk) = max(sst(i,j,iblk),Tf(i,j,iblk))
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO

      elseif (trim(ocn_data_type) == 'ncar') then
         call ocn_data_ncar_init
!        call ocn_data_ncar_init_3D

      elseif (trim(ocn_data_type) == 'hycom') then
         call ocn_data_hycom_init

      elseif (trim(ocn_data_type) == 'box2001') then
         call box2001_data_ocn

      ! uniform forcing options
      elseif (trim(ocn_data_type) == 'uniform_northeast') then
         call uniform_data_ocn('NE',p1)
      elseif (trim(ocn_data_type) == 'uniform_east') then
         call uniform_data_ocn('E',p1)
      elseif (trim(ocn_data_type) == 'uniform_north') then
         call uniform_data_ocn('N',p1)
      elseif (trim(ocn_data_type) == 'calm') then
         call uniform_data_ocn('N',c0) ! directon does not matter for c0
      elseif (trim(ocn_data_type) == 'default') then
         ! don't need to do anything more
      else
         call abort_ice (error_message=subname//' ERROR ocn_data_type unknown = '// &
                         trim(ocn_data_type), file=__FILE__, line=__LINE__)
      endif

      end subroutine init_forcing_ocn

!=======================================================================

      subroutine ocn_freezing_temperature

 ! Compute ocean freezing temperature Tf based on tfrz_option
 ! 'minus1p8'         Tf = -1.8 C (default)
 ! 'linear_salt'      Tf = -depressT * sss
 ! 'mushy'            Tf conforms with mushy layer thermo (ktherm=2)

      use ice_blocks, only: nx_block, ny_block
      use ice_domain, only: nblocks
      use ice_flux, only: sss, Tf

      ! local variables

      integer (kind=int_kind) :: &
         i, j, iblk           ! horizontal indices

      character(len=*), parameter :: subname = '(ocn_freezing_temperature)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            Tf(i,j,iblk) = icepack_sea_freezing_temperature(sss(i,j,iblk))
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      end subroutine ocn_freezing_temperature

!=======================================================================

      subroutine get_forcing_atmo

! Get atmospheric forcing data and interpolate as necessary

      use ice_blocks, only: block, get_block
      use ice_domain, only: nblocks, blocks_ice
      use ice_flux, only: Tair, fsw, flw, frain, fsnow, Qa, rhoa, &
          uatm, vatm, strax, stray, zlvl, wind, swvdr, swvdf, swidr, swidf, &
          potT, sst
      use ice_state, only: aice, trcr
      use ice_grid, only: ANGLET, hm

      integer (kind=int_kind) :: &
         iblk, &              ! block index
         ilo,ihi,jlo,jhi, &   ! beginning and end of physical domain
         modadj, &            ! adjustment to make mod a postive number
         fyear_old, &         ! fyear setting on last timestep
         nt_Tsfc

      type (block) :: &
         this_block           ! block information for current block

      character(len=*), parameter :: subname = '(get_forcing_atmo)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      call ice_timer_start(timer_forcing)

      fyear_old = fyear
      modadj    = abs((min(0,myear-fyear_init)/ycycle+1)*ycycle)
      fyear     = fyear_init + mod(myear-fyear_init+modadj,ycycle)
      if (trim(atm_data_type) /= 'default' .and. &
          (istep <= 1 .or. fyear /= fyear_old)) then
         if (my_task == master_task) then
            write (nu_diag,*) ' Set current forcing data year = ',fyear
         endif
      endif

      call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      !-------------------------------------------------------------------
      ! Read and interpolate atmospheric data
      !-------------------------------------------------------------------

      if (local_debug .and. my_task == master_task) then
         write(nu_diag,*) subname,'fdbg fyear = ',fyear
         write(nu_diag,*) subname,'fdbg atm_data_type = ',trim(atm_data_type)
      endif

      if (trim(atm_data_type) == 'ncar') then
         call ncar_data
      elseif (index(trim(atm_data_type),'JRA55') > 0) then
         call JRA55_data
      elseif (trim(atm_data_type) == 'hadgem') then
         call hadgem_data
      elseif (trim(atm_data_type) == 'monthly') then
         call monthly_data
      elseif (trim(atm_data_type) == 'oned') then
         call oned_data
      elseif (trim(atm_data_type) == 'box2001') then
         call box2001_data_atm
      elseif (trim(atm_data_type) == 'uniform_northeast') then
         call uniform_data_atm('NE')
      elseif (trim(atm_data_type) == 'uniform_north') then
         call uniform_data_atm('N')
      elseif (trim(atm_data_type) == 'uniform_east') then
         call uniform_data_atm('E')
      elseif (trim(atm_data_type) == 'uniform_south') then
         call uniform_data_atm('S')
      elseif (trim(atm_data_type) == 'uniform_west') then
         call uniform_data_atm('W')
      elseif (trim(atm_data_type) == 'calm') then
         call uniform_data_atm('N',c0) ! direction does not matter when c0
      elseif (trim(atm_data_type) == 'hycom') then
         call hycom_atm_data
      !elseif (trim(atm_data_type) == 'uniform_northeast') then
      !elseif (trim(atm_data_type) == 'uniform_east') then
      !elseif (trim(atm_data_type) == 'uniform_north') then
      else    ! default values set in init_flux
         return
      endif

      !-------------------------------------------------------------------
      ! Convert forcing data to fields needed by ice model
      !-------------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks

         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         call prepare_forcing (nx_block, ny_block, &
                               ilo, ihi, jlo, jhi, &
                               hm    (:,:,iblk),   &
                               Tair  (:,:,iblk),   &
                               fsw   (:,:,iblk),   &
                               cldf  (:,:,iblk),   &
                               flw   (:,:,iblk),   &
                               frain (:,:,iblk),   &
                               fsnow (:,:,iblk),   &
                               Qa    (:,:,iblk),   &
                               rhoa  (:,:,iblk),   &
                               uatm  (:,:,iblk),   &
                               vatm  (:,:,iblk),   &
                               strax (:,:,iblk),   &
                               stray (:,:,iblk),   &
                               zlvl  (:,:,iblk),   &
                               wind  (:,:,iblk),   &
                               swvdr (:,:,iblk),   &
                               swvdf (:,:,iblk),   &
                               swidr (:,:,iblk),   &
                               swidf (:,:,iblk),   &
                               potT  (:,:,iblk),   &
                               ANGLET(:,:,iblk),   &
                               trcr  (:,:,nt_Tsfc,iblk), &
                               sst   (:,:,iblk),   &
                               aice  (:,:,iblk) )

      enddo                     ! iblk
      !$OMP END PARALLEL DO

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (swvdr,             halo_info, &
                           field_loc_center,  field_type_scalar, fillvalue=c0)
      call ice_HaloUpdate (swvdf,             halo_info, &
                           field_loc_center,  field_type_scalar, fillvalue=c0)
      call ice_HaloUpdate (swidr,             halo_info, &
                           field_loc_center,  field_type_scalar, fillvalue=c0)
      call ice_HaloUpdate (swidf,             halo_info, &
                           field_loc_center,  field_type_scalar, fillvalue=c0)
      call ice_timer_stop(timer_bound)

      call ice_timer_stop(timer_forcing)

      end subroutine get_forcing_atmo

!=======================================================================

      subroutine get_forcing_ocn (dt)

! Read and interpolate annual climatologies of SSS and SST.
! Restore model SST to data if desired.
! Interpolate ocean fields to U grid if necessary.

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      character(len=*), parameter :: subname = '(get_forcing_ocn)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      call ice_timer_start(timer_forcing)

      if (local_debug .and. my_task == master_task) then
         write(nu_diag,*) subname,'fdbg fyear = ',fyear
         write(nu_diag,*) subname,'fdbg ocn_data_type = ',trim(ocn_data_type)
      endif

      if (trim(ocn_data_type) == 'clim') then
         call ocn_data_clim(dt)
      elseif (trim(ocn_data_type) == 'ncar' .or.  &
              trim(ocn_data_type) == 'ISPOL') then
         call ocn_data_ncar(dt)
      elseif (trim(ocn_data_type) == 'hadgem_sst' .or.  &
              trim(ocn_data_type) == 'hadgem_sst_uvocn') then
         call ocn_data_hadgem(dt)
      elseif (trim(ocn_data_type) == 'oned') then
         call ocn_data_oned
      elseif (trim(ocn_data_type) == 'hycom') then
!         call ocn_data_hycom(dt)
!MHRI: NOT IMPLEMENTED YET
      elseif (trim(ocn_data_type) == 'box2001') then
         call box2001_data_ocn
      ! uniform forcing options
      elseif (trim(ocn_data_type) == 'uniform_northeast') then
! tcraig, not time varying
         call uniform_data_ocn('NE',p1)
      elseif (trim(ocn_data_type) == 'uniform_east') then
         call uniform_data_ocn('E',p1)
      elseif (trim(ocn_data_type) == 'uniform_north') then
         call uniform_data_ocn('N',p1)
      elseif (trim(ocn_data_type) == 'calm') then
         call uniform_data_ocn('N',c0) ! directon does not matter for c0
      endif

      call ice_timer_stop(timer_forcing)

      end subroutine get_forcing_ocn

!=======================================================================

      subroutine read_data (flag, recd, yr, ixm, ixx, ixp, &
                            maxrec, data_file, field_data, &
                            field_loc, field_type)

! If data is at the beginning of a one-year record, get data from
!  the previous year.
! If data is at the end of a one-year record, get data from the
!  following year.
! If no earlier data exists (beginning of fyear_init), then
!  (1) For monthly data, get data from the end of fyear_final.
!  (2) For more frequent data, let the ixm value equal the
!      first value of the year.
! If no later data exists (end of fyear_final), then
!  (1) For monthly data, get data from the beginning of fyear_init.
!  (2) For more frequent data, let the ixp value
!      equal the last value of the year.
! In other words, we assume persistence when daily or 6-hourly
!   data is missing, and we assume periodicity when monthly data
!   is missing.

      use ice_diagnostics, only: debug_model_step

      logical (kind=log_kind), intent(in) :: flag

      integer (kind=int_kind), intent(in) :: &
         recd                , & ! baseline record number
         yr                  , & ! year of forcing data
         ixm, ixx, ixp       , & ! record numbers of 3 data values
                                 ! relative to recd
         maxrec                  ! maximum record value

      real (kind=dbl_kind), dimension(nx_block,ny_block,2,max_blocks), intent(inout) :: &
         field_data              ! 2 values needed for interpolation

      integer (kind=int_kind), intent(in) :: &
           field_loc, &      ! location of field on staggered grid
           field_type        ! type of field (scalar, vector, angle)

      ! local variables

      character (char_len_long) :: &
         data_file               ! data file to be read

      integer (kind=int_kind) :: &
         nbits            , & ! = 32 for single precision, 64 for double
         nrec             , & ! record number to read
         n2, n4           , & ! like ixm and ixp, but
                              ! adjusted at beginning and end of data
         arg                  ! value of time argument in field_data

      character(len=*), parameter :: subname = '(read_data)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      call ice_timer_start(timer_readwrite)  ! reading/writing

      nbits = 64              ! double precision data

      if (istep1 > debug_model_step) debug_forcing = .true.  !! debugging

      if (my_task==master_task .and. (debug_forcing)) then
         write(nu_diag,*) '  ', trim(data_file)
      endif

      if (flag) then

      !-----------------------------------------------------------------
      ! Initialize record counters
      ! (n2, n4 will change only at the very beginning or end of
      !  a forcing cycle.)
      !-----------------------------------------------------------------
         n2 = ixm
         n4 = ixp
         arg = 0

      !-----------------------------------------------------------------
      ! read data
      !-----------------------------------------------------------------

         if (ixm /= -99) then
         ! currently in first half of data interval
            if (ixx <= 1) then
               if (yr > fyear_init) then ! get data from previous year
                  call file_year (data_file, yr-1)
               else             ! yr = fyear_init, no prior data exists
                  if (maxrec > 12) then ! extrapolate from first record
                     if (ixx == 1) n2 = ixx
                  else          ! go to end of fyear_final
                     call file_year (data_file, fyear_final)
                  endif
               endif            ! yr > fyear_init
            endif               ! ixx <= 1

            call ice_open (nu_forcing, data_file, nbits)

            arg = 1
            nrec = recd + n2
            call ice_read (nu_forcing, nrec, field_data(:,:,arg,:), &
                           'rda8', debug_forcing, field_loc, field_type)

            if (ixx==1 .and. my_task == master_task) close(nu_forcing)
         endif                  ! ixm ne -99

         ! always read ixx data from data file for current year
         call file_year (data_file, yr)
         call ice_open (nu_forcing, data_file, nbits)

         arg = arg + 1
         nrec = recd + ixx
         call ice_read (nu_forcing, nrec, field_data(:,:,arg,:), &
                        'rda8', debug_forcing, field_loc, field_type)

         if (ixp /= -99) then
         ! currently in latter half of data interval
            if (ixx==maxrec) then
               if (yr < fyear_final) then ! get data from following year
                  if (my_task == master_task) close(nu_forcing)
                  call file_year (data_file, yr+1)
                  call ice_open (nu_forcing, data_file, nbits)
               else             ! yr = fyear_final, no more data exists
                  if (maxrec > 12) then ! extrapolate from ixx
                     n4 = ixx
                  else          ! go to beginning of fyear_init
                     if (my_task == master_task) close(nu_forcing)
                     call file_year (data_file, fyear_init)

                     call ice_open (nu_forcing, data_file, nbits)

                  endif
               endif            ! yr < fyear_final
            endif               ! ixx = maxrec

            arg = arg + 1
            nrec = recd + n4
            call ice_read (nu_forcing, nrec, field_data(:,:,arg,:), &
                           'rda8', debug_forcing, field_loc, field_type)
         endif                  ! ixp /= -99

         if (my_task == master_task) close(nu_forcing)

      endif                     ! flag

      call ice_timer_stop(timer_readwrite)  ! reading/writing

      end subroutine read_data

!=======================================================================

      subroutine read_data_nc (flag, recd, yr, ixm, ixx, ixp, &
                            maxrec, data_file, fieldname, field_data, &
                            field_loc, field_type)

! If data is at the beginning of a one-year record, get data from
!  the previous year.
! If data is at the end of a one-year record, get data from the
!  following year.
! If no earlier data exists (beginning of fyear_init), then
!  (1) For monthly data, get data from the end of fyear_final.
!  (2) For more frequent data, let the ixm value equal the
!      first value of the year.
! If no later data exists (end of fyear_final), then
!  (1) For monthly data, get data from the beginning of fyear_init.
!  (2) For more frequent data, let the ixp value
!      equal the last value of the year.
! In other words, we assume persistence when daily or 6-hourly
!   data is missing, and we assume periodicity when monthly data
!   is missing.
!
! Adapted by Alison McLaren, Met Office from read_data

      use ice_diagnostics, only: debug_model_step

      logical (kind=log_kind), intent(in) :: flag

      integer (kind=int_kind), intent(in) :: &
         recd                , & ! baseline record number
         yr                  , & ! year of forcing data
         ixm, ixx, ixp       , & ! record numbers of 3 data values
                                 ! relative to recd
         maxrec                  ! maximum record value

      character (char_len_long) :: &
         data_file               ! data file to be read

      character (char_len), intent(in) :: &
         fieldname               ! field name in netCDF file

      integer (kind=int_kind), intent(in) :: &
         field_loc, &            ! location of field on staggered grid
         field_type              ! type of field (scalar, vector, angle)

      real (kind=dbl_kind), dimension(nx_block,ny_block,2,max_blocks), intent(out) :: &
         field_data              ! 2 values needed for interpolation

      ! local variables

      integer (kind=int_kind) :: &
         nrec             , & ! record number to read
         n2, n4           , & ! like ixm and ixp, but
                              ! adjusted at beginning and end of data
         arg              , & ! value of time argument in field_data
         fid                  ! file id for netCDF routines

      character(len=*), parameter :: subname = '(read_data_nc)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      call ice_timer_start(timer_readwrite)  ! reading/writing

      if (istep1 > debug_model_step) debug_forcing = .true.  !! debugging

      if (my_task==master_task .and. (debug_forcing)) then
         write(nu_diag,*) '  ', trim(data_file)
      endif

      if (flag) then

      !-----------------------------------------------------------------
      ! Initialize record counters
      ! (n2, n4 will change only at the very beginning or end of
      !  a forcing cycle.)
      !-----------------------------------------------------------------
         n2 = ixm
         n4 = ixp
         arg = 0

      !-----------------------------------------------------------------
      ! read data
      !-----------------------------------------------------------------

         if (ixm /= -99) then
         ! currently in first half of data interval
            if (ixx <= 1) then
               if (yr > fyear_init) then ! get data from previous year
                  call file_year (data_file, yr-1)
               else             ! yr = fyear_init, no prior data exists
                  if (maxrec > 12) then ! extrapolate from first record
                     if (ixx == 1) n2 = ixx
                  else          ! go to end of fyear_final
                     call file_year (data_file, fyear_final)
                  endif
               endif            ! yr > fyear_init
            endif               ! ixx <= 1

            call ice_open_nc (data_file, fid)

            arg = 1
            nrec = recd + n2

            call ice_read_nc &
                 (fid, nrec, fieldname, field_data(:,:,arg,:), debug_forcing, &
                  field_loc, field_type)

            if (ixx==1) call ice_close_nc(fid)
         endif                  ! ixm ne -99

         ! always read ixx data from data file for current year
         call file_year (data_file, yr)
         call ice_open_nc (data_file, fid)

         arg = arg + 1
         nrec = recd + ixx

         call ice_read_nc &
              (fid, nrec, fieldname, field_data(:,:,arg,:), debug_forcing, &
               field_loc, field_type)

         if (ixp /= -99) then
         ! currently in latter half of data interval
            if (ixx==maxrec) then
               if (yr < fyear_final) then ! get data from following year
                  call ice_close_nc(fid)
                  call file_year (data_file, yr+1)
                  call ice_open_nc (data_file, fid)
               else             ! yr = fyear_final, no more data exists
                  if (maxrec > 12) then ! extrapolate from ixx
                     n4 = ixx
                  else          ! go to beginning of fyear_init
                     call ice_close_nc(fid)
                     call file_year (data_file, fyear_init)
                     call ice_open_nc (data_file, fid)

                  endif
               endif            ! yr < fyear_final
            endif               ! ixx = maxrec

            arg = arg + 1
            nrec = recd + n4

            call ice_read_nc &
                 (fid, nrec, fieldname, field_data(:,:,arg,:), debug_forcing, &
                  field_loc, field_type)
         endif                  ! ixp /= -99

         call ice_close_nc(fid)

      endif                     ! flag

      call ice_timer_stop(timer_readwrite)  ! reading/writing

      end subroutine read_data_nc

!=======================================================================

      subroutine read_data_nc_hycom (flag, recd,  &
                            data_file, fieldname, field_data, &
                            field_loc, field_type)

!  Data is assumed to cover the entire time period of simulation.
!  It is not bounded by start of year nor end of year
!  Data must be accesible both before and after (or on) the point in time
!    Assume increasing timeaxis within the forcing files, but they do not
!      have to be equal spaced. Read time vector from "MT" in "init_hycom"
!
! Adapted by Mads Hvid Ribergaard, DMI from read_data_nc

      use ice_diagnostics, only: debug_model_step
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_readwrite

      logical (kind=log_kind), intent(in) :: flag

      integer (kind=int_kind), intent(in) :: &
         recd                    ! baseline record number

      character (char_len_long) :: &
         data_file               ! data file to be read
      character (char_len), intent(in) :: &
         fieldname               ! field name in netCDF file

      integer (kind=int_kind), intent(in) :: &
           field_loc, &      ! location of field on staggered grid
           field_type        ! type of field (scalar, vector, angle)

      real (kind=dbl_kind), dimension(nx_block,ny_block,2,max_blocks), &
         intent(out) :: &
         field_data              ! 2 values needed for interpolation

      ! local variables
      integer (kind=int_kind) :: &
         fid                  ! file id for netCDF routines

      character(len=*), parameter :: subname = '(read_data_nc_hycom)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      call ice_timer_start(timer_readwrite)  ! reading/writing

      if (istep1 > debug_model_step) debug_forcing = .true.  !! debugging

      if (my_task==master_task .and. (debug_forcing)) then
         write(nu_diag,*) '  ', trim(data_file)
      endif

      if (flag) then

         call ice_open_nc (data_file, fid)
      !-----------------------------------------------------------------
      ! read data
      !-----------------------------------------------------------------
         call ice_read_nc &
               (fid, recd , fieldname, field_data(:,:,1,:), debug_forcing, &
                field_loc, field_type)

         call ice_read_nc &
              (fid, recd+1, fieldname, field_data(:,:,2,:), debug_forcing, &
               field_loc, field_type)

         call ice_close_nc(fid)

      endif                     ! flag

      call ice_timer_stop(timer_readwrite)  ! reading/writing

      end subroutine read_data_nc_hycom

!=======================================================================

      subroutine read_clim_data (readflag, recd, ixm, ixx, ixp, &
                                 data_file, field_data, &
                                 field_loc, field_type)

! Read data needed for interpolation, as in read_data.
! Assume a one-year cycle of climatological data, so that there is
!  no need to get data from other years or to extrapolate data beyond
!  the forcing time period.

      use ice_diagnostics, only: debug_model_step

      logical (kind=log_kind),intent(in) :: readflag

      integer (kind=int_kind), intent(in) :: &
        recd            , & ! baseline record number
        ixm,ixx,ixp         ! record numbers of 3 data values
                            ! relative to recd

      character (char_len_long), intent(in) ::  data_file

      integer (kind=int_kind), intent(in) :: &
           field_loc, &      ! location of field on staggered grid
           field_type        ! type of field (scalar, vector, angle)

      real (kind=dbl_kind), dimension(nx_block,ny_block,2,max_blocks), intent(inout) :: &
        field_data         ! 2 values needed for interpolation

      ! local variables

      integer (kind=int_kind) :: &
        nbits          , & ! = 32 for single precision, 64 for double
        nrec           , & ! record number to read
        arg                ! value of time argument in field_data

      character(len=*), parameter :: subname = '(read_clim_data)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      call ice_timer_start(timer_readwrite)  ! reading/writing

      nbits = 64                ! double precision data

      if (istep1 > debug_model_step) debug_forcing = .true.  !! debugging

      if (my_task==master_task .and. (debug_forcing)) &
        write(nu_diag,*) '  ', trim(data_file)

      if (readflag) then

      !-----------------------------------------------------------------
      ! read data
      !-----------------------------------------------------------------

         call ice_open (nu_forcing, data_file, nbits)

         arg = 0
         if (ixm /= -99) then
            arg = 1
            nrec = recd + ixm
            call ice_read (nu_forcing, nrec, field_data(:,:,arg,:), &
                           'rda8', debug_forcing, field_loc, field_type)
         endif

         arg = arg + 1
         nrec = recd + ixx
         call ice_read (nu_forcing, nrec, field_data(:,:,arg,:), &
                        'rda8', debug_forcing, field_loc, field_type)

         if (ixp /= -99) then
            arg = arg + 1
            nrec = recd + ixp
            call ice_read (nu_forcing, nrec, field_data(:,:,arg,:), &
                           'rda8', debug_forcing, field_loc, field_type)
         endif

         if (my_task == master_task) close (nu_forcing)
      endif                     ! readflag

      call ice_timer_stop(timer_readwrite)  ! reading/writing

      end subroutine read_clim_data

!=======================================================================

      subroutine read_clim_data_nc (readflag, recd, ixm, ixx, ixp, &
                                 data_file, fieldname, field_data, &
                                 field_loc, field_type)

! Read data needed for interpolation, as in read_data.
! Assume a one-year cycle of climatological data, so that there is
!  no need to get data from other years or to extrapolate data beyond
!  the forcing time period.

      use ice_diagnostics, only: debug_model_step

      logical (kind=log_kind),intent(in) :: readflag

      integer (kind=int_kind), intent(in) :: &
        recd            , & ! baseline record number
        ixm,ixx,ixp         ! record numbers of 3 data values
                            ! relative to recd

      character (char_len_long), intent(in) ::  data_file

      character (char_len), intent(in) :: &
         fieldname               ! field name in netCDF file

      integer (kind=int_kind), intent(in) :: &
           field_loc, &      ! location of field on staggered grid
           field_type        ! type of field (scalar, vector, angle)

      real (kind=dbl_kind), dimension(nx_block,ny_block,2,max_blocks), intent(out) :: &
        field_data         ! 2 values needed for interpolation

      ! local variables

      integer (kind=int_kind) :: &
        nrec           , & ! record number to read
        arg            , & ! value of time argument in field_data
        fid                ! file id for netCDF routines

      character(len=*), parameter :: subname = '(read_clim_data_nc)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      call ice_timer_start(timer_readwrite)  ! reading/writing

      if (istep1 > debug_model_step) debug_forcing = .true.  !! debugging

      if (my_task==master_task .and. (debug_forcing)) &
        write(nu_diag,*) '  ', trim(data_file)

      if (readflag) then

      !-----------------------------------------------------------------
      ! read data
      !-----------------------------------------------------------------

         call ice_open_nc (data_file, fid)

         arg = 0
         if (ixm /= -99) then
            arg = 1
            nrec = recd + ixm
            call ice_read_nc &
                 (fid, nrec, fieldname, field_data(:,:,arg,:), &
                  debug_forcing, field_loc, field_type)
         endif

         arg = arg + 1
         nrec = recd + ixx
         call ice_read_nc &
                 (fid, nrec, fieldname, field_data(:,:,arg,:), &
                  debug_forcing, field_loc, field_type)

         if (ixp /= -99) then
            arg = arg + 1
            nrec = recd + ixp
            call ice_read_nc &
                 (fid, nrec, fieldname, field_data(:,:,arg,:), &
                  debug_forcing, field_loc, field_type)
         endif

         if (my_task == master_task) call ice_close_nc (fid)
      endif                     ! readflag

      call ice_timer_stop(timer_readwrite)  ! reading/writing

      end subroutine read_clim_data_nc

!=======================================================================

      subroutine interp_coeff_monthly (recslot)

! Compute coefficients for interpolating monthly data to current time step.

      integer (kind=int_kind), intent(in) :: &
          recslot         ! slot (1 or 2) for current record

      ! local variables

      real (kind=dbl_kind) :: &
          secday       , & ! seconds in day
          tt           , & ! days elapsed in current year
          t1, t2           ! days elapsed at month midpoint

      real (kind=dbl_kind) :: &
          daymid(0:13)     ! month mid-points

      character(len=*), parameter :: subname = '(interp_coeff_monthly)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      call icepack_query_parameters(secday_out=secday)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      daymid(1:13) = 14._dbl_kind   ! time frame ends 0 sec into day 15
      daymid(0)    = 14._dbl_kind - daymo(12)  ! Dec 15, 0 sec

      ! compute days since Jan 1, 00h, yday is the day counter for the year
      tt = real(yday-1,kind=dbl_kind) + real(msec,kind=dbl_kind)/secday

      ! Find neighboring times

      if (recslot==2) then      ! first half of month
        t2 = daycal(mmonth) + daymid(mmonth)   ! midpoint, current month
        if (mmonth == 1) then
          t1 = daymid(0)                 ! Dec 15 (0 sec)
        else
          t1 = daycal(mmonth-1) + daymid(mmonth-1) ! midpoint, previous month
        endif
      else                      ! second half of month
        t1 = daycal(mmonth) + daymid(mmonth)    ! midpoint, current month
        t2 = daycal(mmonth+1) + daymid(mmonth+1)! day 15 of next month (0 sec)
      endif

      if (tt < t1 .or. tt > t2) then
        write(nu_diag,*) subname,' ERROR in tt',tt,t1,t2
        call abort_ice (error_message=subname//' ERROR in tt', &
           file=__FILE__, line=__LINE__)
      endif

      ! Compute coefficients
      c1intp = (t2 - tt) / (t2 - t1)
      c2intp =  c1 - c1intp

      end subroutine interp_coeff_monthly

!=======================================================================

      subroutine interp_coeff (recnum, recslot, secint, dataloc)

! Compute coefficients for interpolating data to current time step.
! Works for any data interval that divides evenly into a
!  year (daily, 6-hourly, etc.)
! Use interp_coef_monthly for monthly data.

      integer (kind=int_kind), intent(in) :: &
          recnum      , & ! record number for current data value
          recslot     , & ! spline slot for current record
          dataloc         ! = 1 for data located in middle of time interval
                          ! = 2 for date located at end of time interval

      real (kind=dbl_kind), intent(in) :: &
          secint                    ! seconds in data interval

      ! local variables

      real (kind=dbl_kind) :: &
          secday           ! seconds in a day

      real (kind=dbl_kind) :: &
          tt           , & ! seconds elapsed in current year
          t1, t2       , & ! seconds elapsed at data points
          rcnum            ! recnum => dbl_kind

      character(len=*), parameter :: subname = '(interp_coeff)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      call icepack_query_parameters(secday_out=secday)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      ! compute seconds since Jan 1, 00h, yday is the day counter for the year
      tt = real(yday-1,kind=dbl_kind)*secday + real(msec,kind=dbl_kind)

      ! Find neighboring times
      rcnum = real(recnum,kind=dbl_kind)
      if (recslot==2) then           ! current record goes in slot 2
         if (dataloc==1) then        ! data located at middle of interval
            t2 = (rcnum-p5)*secint
         else                        !  data located at end of interval
            t2 = rcnum*secint
         endif
         t1 = t2 - secint            !  - 1 interval
      else                           ! recslot = 1
         if (dataloc==1) then        ! data located at middle of interval
            t1 = (rcnum-p5)*secint
         else
            t1 = rcnum*secint        ! data located at end of interval
         endif
         t2 = t1 + secint            !  + 1 interval
      endif

      ! Compute coefficients
      c1intp =  abs((t2 - tt) / (t2 - t1))
      c2intp =  c1 - c1intp

      if (local_debug .and. my_task == master_task) then
         write(nu_diag,*) subname,'fdbg yday,sec = ',yday,msec
         write(nu_diag,*) subname,'fdbg tt = ',tt
         write(nu_diag,*) subname,'fdbg c12intp = ',c1intp,c2intp
      endif

      end subroutine interp_coeff

!=======================================================================

      subroutine interp_coeff2 (tt, t1, t2)

! Compute coefficients for interpolating data to current time step.
! Works for any data interval using decimal daynumbers

      ! local variables
      real (kind=dbl_kind), intent(in) :: &
          tt      , &  ! current decimal daynumber
          t1, t2       ! first+last decimal daynumber
      character(len=*), parameter :: subname = '(interp_coeff2)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      ! Compute coefficients
      c1intp =  abs((t2 - tt) / (t2 - t1))
      c2intp =  c1 - c1intp

      end subroutine interp_coeff2

!=======================================================================

      subroutine interpolate_data (field_data, field)

! Linear interpolation

! author: Elizabeth C. Hunke, LANL

      use ice_domain, only: nblocks

      real (kind=dbl_kind), dimension(nx_block,ny_block,2,max_blocks), intent(in) :: &
        field_data    ! 2 values used for interpolation

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks), intent(out) :: &
        field         ! interpolated field

      ! local variables

      integer (kind=int_kind) :: i,j, iblk

      character(len=*), parameter :: subname = '(interpolate data)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            field(i,j,iblk) = c1intp * field_data(i,j,1,iblk) &
                            + c2intp * field_data(i,j,2,iblk)
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine interpolate_data

!=======================================================================

      subroutine interpolate_wavespec_data (field_data, field)

! Linear interpolation

! author: Elizabeth C. Hunke, LANL

      use ice_domain, only: nblocks

      real (kind=dbl_kind), dimension(nx_block,ny_block,nfreq,2,max_blocks), intent(in) :: &
        field_data    ! 2 values used for interpolation

      real (kind=dbl_kind), dimension(nx_block,ny_block,nfreq,max_blocks), intent(out) :: &
        field         ! interpolated field

      ! local variables

      integer (kind=int_kind) :: i,j, iblk, freq

      character(len=*), parameter :: subname = '(interpolate data)'

      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
         do freq = 1, nfreq
            field(i,j,freq,iblk) = c1intp * field_data(i,j,freq,1,iblk) &
                            + c2intp * field_data(i,j,freq,2,iblk)
         enddo
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine interpolate_wavespec_data


!=======================================================================

      subroutine file_year (data_file, yr)

! Construct the correct name of the atmospheric data file
! to be read, given the year and assuming the naming convention
! that filenames end with 'yyyy.dat' or 'yyyy.r' or 'yyyy.nc'.

      character (char_len_long), intent(inout) ::  data_file

      integer (kind=int_kind), intent(in) :: yr

      character (char_len_long) :: tmpname

      integer (kind=int_kind) :: i

      character(len=*), parameter :: subname = '(file_year)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      if (trim(atm_data_type) == 'hadgem') then ! netcdf
         i = index(data_file,'.nc') - 5
         tmpname = data_file
         write(data_file,'(a,i4.4,a)') tmpname(1:i), yr, '.nc'
      elseif (index(trim(atm_data_type),'JRA55') > 0) then ! netcdf
         i = index(data_file,'.nc') - 5
         tmpname = data_file
         write(data_file,'(a,i4.4,a)') tmpname(1:i), yr, '.nc'
      else                                     ! LANL/NCAR naming convention
         i = index(data_file,'.dat') - 5
         tmpname = data_file
         write(data_file,'(a,i4.4,a)') tmpname(1:i), yr, '.dat'
      endif

      end subroutine file_year

!=======================================================================

      subroutine prepare_forcing (nx_block, ny_block, &
                                  ilo, ihi, jlo, jhi, &
                                  hm,                 &
                                  Tair,     fsw,      &
                                  cldf,     flw,      &
                                  frain,    fsnow,    &
                                  Qa,       rhoa,     &
                                  uatm,     vatm,     &
                                  strax,    stray,    &
                                  zlvl,     wind,     &
                                  swvdr,    swvdf,    &
                                  swidr,    swidf,    &
                                  potT,     ANGLET,   &
                                  Tsfc,     sst,      &
                                  aice)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ilo,ihi,jlo,jhi       ! beginning and end of physical domain

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) :: &
         ANGLET  , & ! ANGLE converted to T-cells
         Tsfc    , & ! ice skin temperature
         sst     , & ! sea surface temperature
         aice    , & ! ice area fraction
         hm          ! land mask

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(inout) :: &
         fsw     , & ! incoming shortwave radiation (W/m^2)
         cldf    , & ! cloud fraction
         frain   , & ! rainfall rate (kg/m^2 s)
         fsnow   , & ! snowfall rate (kg/m^2 s)
         Tair    , & ! air temperature  (K)
         Qa      , & ! specific humidity (kg/kg)
         rhoa    , & ! air density (kg/m^3)
         uatm    , & ! wind velocity components (m/s)
         vatm    , &
         strax   , & ! wind stress components (N/m^2)
         stray   , &
         zlvl    , & ! atm level height (m)
         wind    , & ! wind speed (m/s)
         flw     , & ! incoming longwave radiation (W/m^2)
         swvdr   , & ! sw down, visible, direct  (W/m^2)
         swvdf   , & ! sw down, visible, diffuse (W/m^2)
         swidr   , & ! sw down, near IR, direct  (W/m^2)
         swidf   , & ! sw down, near IR, diffuse (W/m^2)
         potT        ! air potential temperature  (K)

      ! local variables

      integer (kind=int_kind) :: &
         i, j

      real (kind=dbl_kind) :: workx, worky, &
         precip_factor, zlvl0, secday, Tffresh, puny

      logical (kind=log_kind) :: calc_strair

      character(len=*), parameter :: subname = '(prepare_forcing)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      call icepack_query_parameters(Tffresh_out=Tffresh, puny_out=puny)
      call icepack_query_parameters(secday_out=secday)
      call icepack_query_parameters(calc_strair_out=calc_strair)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      do j = jlo, jhi
      do i = ilo, ihi

         zlvl0 = c10 ! default

      !-----------------------------------------------------------------
      ! make sure interpolated values are physically realistic
      !-----------------------------------------------------------------
         cldf (i,j) = max(min(cldf(i,j),c1),c0)
         fsw  (i,j) = max(fsw(i,j),c0)
         fsnow(i,j) = max(fsnow(i,j),c0)
         rhoa (i,j) = max(rhoa(i,j),c0)
         Qa   (i,j) = max(Qa(i,j),c0)

!        if (rhoa(i,j) .lt. puny) rhoa(i,j) = 1.3_dbl_kind
!        if (Tair(i,j) .lt. puny) Tair(i,j) = Tffresh
!        if (Qa(i,j) .lt. puny) Qa(i,j) = 0.0035_dbl_kind
      enddo                     ! i
      enddo                     ! j

      !-----------------------------------------------------------------
      ! calculations specific to datasets
      !-----------------------------------------------------------------

      if (trim(atm_data_type) == 'ncar') then

         ! precip is in mm/month

         zlvl0 = c10

         do j = jlo, jhi
         do i = ilo, ihi
            ! correct known biases in NCAR data (as in CESM latm)
            Qa (i,j) = Qa (i,j) * 0.94_dbl_kind
            fsw(i,j) = fsw(i,j) * 0.92_dbl_kind

            ! downward longwave as in Parkinson and Washington (1979)
            call longwave_parkinson_washington(Tair(i,j), cldf(i,j), &
                                              flw(i,j))
         enddo
         enddo

      elseif (trim(atm_data_type) == 'oned') then  ! rectangular grid

         ! precip is in kg/m^2/s

         zlvl0 = c10

         do j = jlo, jhi
         do i = ilo, ihi

      !-----------------------------------------------------------------
      ! compute downward longwave as in Parkinson and Washington (1979)
      !-----------------------------------------------------------------

            ! downward longwave as in Parkinson and Washington (1979)
            call longwave_parkinson_washington(Tair(i,j), cldf(i,j), &
                                               flw(i,j))

            ! longwave based on Rosati and Miyakoda, JPO 18, p. 1607 (1988)
!            call longwave_rosati_miyakoda(cldf(i,j), Tsfc(i,j), &
!                                          aice(i,j), sst(i,j),  &
!                                          Qa(i,j),   Tair(i,j), &
!                                          hm(i,j),   flw(i,j))
         enddo
         enddo

      endif                     ! atm_data_type

      !-----------------------------------------------------------------
      ! Compute other fields needed by model
      !-----------------------------------------------------------------

      ! convert precipitation units to kg/m^2 s
      if (trim(precip_units) == 'mm_per_month') then
         precip_factor = c12/(secday*real(days_per_year,kind=dbl_kind))
      elseif (trim(precip_units) == 'mm_per_day') then
         precip_factor = c1/secday
      elseif (trim(precip_units) == 'mm_per_sec' .or. &
              trim(precip_units) == 'mks') then
         precip_factor = c1    ! mm/sec = kg/m^2 s
      elseif (trim(precip_units) == 'm_per_sec') then
         precip_factor = c1000
      endif

      do j = jlo, jhi
      do i = ilo, ihi

         zlvl(i,j) = zlvl0
         potT(i,j) = Tair(i,j)

        ! divide shortwave into spectral bands
         swvdr(i,j) = fsw(i,j)*frcvdr        ! visible direct
         swvdf(i,j) = fsw(i,j)*frcvdf        ! visible diffuse
         swidr(i,j) = fsw(i,j)*frcidr        ! near IR direct
         swidf(i,j) = fsw(i,j)*frcidf        ! near IR diffuse

        ! convert precipitation units to kg/m^2 s
         fsnow(i,j) = fsnow(i,j) * precip_factor
      enddo                     ! i
      enddo                     ! j

      ! determine whether precip is rain or snow
      ! HadGEM forcing provides separate snowfall and rainfall rather
      ! than total precipitation
      if (trim(atm_data_type) /= 'hadgem') then

        do j = jlo, jhi
        do i = ilo, ihi
           frain(i,j) = c0
           if (Tair(i,j) >= Tffresh) then
               frain(i,j) = fsnow(i,j)
               fsnow(i,j) = c0
           endif
        enddo                     ! i
        enddo                     ! j

      endif

      if (calc_strair) then

        if (rotate_wind) then
          do j = jlo, jhi
          do i = ilo, ihi
             wind(i,j) = sqrt(uatm(i,j)**2 + vatm(i,j)**2)
      !-----------------------------------------------------------------
      ! Rotate zonal/meridional vectors to local coordinates.
      ! Velocity comes in on T grid, but is oriented geographically ---
      ! need to rotate to pop-grid FIRST using ANGLET
      ! then interpolate to the U-cell centers  (otherwise we
      ! interpolate across the pole).
      ! Use ANGLET which is on the T grid !
      ! Atmo variables are needed in T cell centers in subroutine
      ! atmo_boundary_layer, and are interpolated to the U grid later as
      ! necessary.
      !-----------------------------------------------------------------
             workx      = uatm(i,j) ! wind velocity, m/s
             worky      = vatm(i,j)
             uatm (i,j) = workx*cos(ANGLET(i,j)) & ! convert to POP grid
                        + worky*sin(ANGLET(i,j))   ! note uatm, vatm, wind
             vatm (i,j) = worky*cos(ANGLET(i,j)) & !  are on the T-grid here
                        - workx*sin(ANGLET(i,j))
          enddo                     ! i
          enddo                     ! j
        else ! not rotated
          do j = jlo, jhi
          do i = ilo, ihi
             wind(i,j) = sqrt(uatm(i,j)**2 + vatm(i,j)**2)
          enddo                     ! i
          enddo                     ! j
        endif ! rotated

      else ! strax, stray, wind are read from files

        if (rotate_wind) then
          do j = jlo, jhi
          do i = ilo, ihi
             workx      = strax(i,j) ! wind stress
             worky      = stray(i,j)
             strax(i,j) = workx*cos(ANGLET(i,j)) & ! convert to POP grid
                        + worky*sin(ANGLET(i,j))   ! note strax, stray, wind
             stray(i,j) = worky*cos(ANGLET(i,j)) & !  are on the T-grid here
                        - workx*sin(ANGLET(i,j))
          enddo                     ! i
          enddo                     ! j
        else ! not rotated
          ! wind (speed) is already read from file, so all is in place
        endif ! rotated

      endif                   ! calc_strair

      end subroutine prepare_forcing

!=======================================================================

      subroutine longwave_parkinson_washington(Tair, cldf, flw)

      ! compute downward longwave as in Parkinson and Washington (1979)
      ! (for now)
      ! Parkinson, C. L. and W. M. Washington (1979),
      ! Large-scale numerical-model of sea ice,
      ! JGR, 84, 311-337, doi:10.1029/JC084iC01p00311

      real(kind=dbl_kind), intent(in) :: &
           Tair , & ! air temperature  (K)
           cldf     ! cloud fraction

      real(kind=dbl_kind), intent(out) :: &
           flw      ! incoming longwave radiation (W/m^2)

      real(kind=dbl_kind) :: &
           Tffresh, stefan_boltzmann

      character(len=*), parameter :: subname = '(longwave_parkinson_washington)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      call icepack_query_parameters(Tffresh_out=Tffresh, &
           stefan_boltzmann_out=stefan_boltzmann)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      flw = stefan_boltzmann*Tair**4 &
             * (c1 - 0.261_dbl_kind &
             * exp(-7.77e-4_dbl_kind*(Tffresh - Tair)**2)) &
             * (c1 + 0.275_dbl_kind*cldf)

      end subroutine longwave_parkinson_washington

!=======================================================================

      subroutine longwave_rosati_miyakoda(cldf, Tsfc, &
                                          aice, sst,  &
                                          Qa,   Tair, &
                                          hm,   flw)

      ! based on
      ! Rosati, A. and K. Miyakoda (1988),
      ! A general-circulation model for upper ocean simulation,
      ! J. Physical Oceanography, 18, 1601-1626,
      ! doi:10.1175/1520-0485(1988)018<1601:AGCMFU>2.0.CO;2

      real(kind=dbl_kind), intent(in) :: &
           cldf , & ! cloud fraction
           Tsfc , & ! ice skin temperature
           aice , & ! ice area fraction
           sst  , & ! sea surface temperature
           Qa   , & ! specific humidity (kg/kg)
           Tair , & ! air temperature  (K)
           hm       ! land mask

      real(kind=dbl_kind), intent(out) :: &
           flw      ! incoming longwave radiation (W/m^2)

      real(kind=dbl_kind) :: &
           fcc  , & ! cloudiness modification
           sstk , & ! ice/ocean surface temperature (K)
           rtea , & ! square root of the vapour pressure
           ptem , & ! potential air temperature (K)
           qlwm

      real(kind=dbl_kind) :: &
           Tffresh, stefan_boltzmann, emissivity

      character(len=*), parameter :: subname = '(longwave_rosati_miyakoda)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      call icepack_query_parameters(Tffresh_out=Tffresh, &
           stefan_boltzmann_out=stefan_boltzmann, &
           emissivity_out=emissivity)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      fcc = c1 - 0.8_dbl_kind * cldf
      sstk = (Tsfc * aice &
           + sst * (c1 - aice)) + Tffresh
      rtea = sqrt(c1000*Qa /  &
           (0.622_dbl_kind+0.378_dbl_kind*Qa))
      ptem = Tair    ! get this from stability?
      qlwm = ptem * ptem * ptem  &
                 * ( ptem*(0.39_dbl_kind-0.05_dbl_kind*rtea)*fcc  &
                 + c4*(sstk-ptem) )
      flw = emissivity*stefan_boltzmann * ( sstk**4 - qlwm )
      flw = flw * hm ! land mask

      end subroutine longwave_rosati_miyakoda

!=======================================================================
! NCAR atmospheric forcing
!=======================================================================

      subroutine ncar_files (yr)

      ! Construct filenames based on the LANL naming conventions for NCAR data.
      ! Edit for other directory structures or filenames.
      ! Note: The year number in these filenames does not matter, because
      !       subroutine file_year will insert the correct year.
      ! Note: atm_data_dir may have NCAR_bulk or not
      !
      ! atm_data_type should be 'ncar'
      ! atm_dat_dir should be ${CICE_DATA_root}/forcing/$grid/[NCAR_bulk,'']
      !    atm_data_dir should be set to ${CICE_DATA_root}/forcing/$grid/[JRA55,JRA55do,'']
      !       NCAR_bulk at the end of the atm_data_dir is optional to provide backwards
      !       compatibility and if not included, will be appended automaticaly.
      !       The grid is typically gx1, gx3, tx1, or similar.

      integer (kind=int_kind), intent(in) :: &
           yr                   ! current forcing year

      character (char_len_long) :: &
           atm_data_dir_extra   ! atm_dat_dir extra if needed

      integer (kind=int_kind) :: &
           strind   ! string index

      character(len=*), parameter :: subname = '(ncar_files)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      ! decide whether NCAR_bulk is part of atm_data_dir and set atm_data_dir_extra
      atm_data_dir_extra = '/NCAR_bulk'
      strind = index(trim(atm_data_dir),'NCAR_bulk')
      if (strind > 0) then
         atm_data_dir_extra = ''
      endif

      fsw_file   = trim(atm_data_dir)//trim(atm_data_dir_extra)//'/MONTHLY/swdn.1996.dat'
      call file_year(fsw_file,yr)

      flw_file   = trim(atm_data_dir)//trim(atm_data_dir_extra)//'/MONTHLY/cldf.1996.dat'
      call file_year(flw_file,yr)

      rain_file  = trim(atm_data_dir)//trim(atm_data_dir_extra)//'/MONTHLY/prec.1996.dat'
      call file_year(rain_file,yr)

      uwind_file = trim(atm_data_dir)//trim(atm_data_dir_extra)//'/4XDAILY/u_10.1996.dat'
      call file_year(uwind_file,yr)

      vwind_file = trim(atm_data_dir)//trim(atm_data_dir_extra)//'/4XDAILY/v_10.1996.dat'
      call file_year(vwind_file,yr)

      tair_file  = trim(atm_data_dir)//trim(atm_data_dir_extra)//'/4XDAILY/t_10.1996.dat'
      call file_year(tair_file,yr)

      humid_file = trim(atm_data_dir)//trim(atm_data_dir_extra)//'/4XDAILY/q_10.1996.dat'
      call file_year(humid_file,yr)

      rhoa_file  = trim(atm_data_dir)//trim(atm_data_dir_extra)//'/4XDAILY/dn10.1996.dat'
      call file_year(rhoa_file,yr)

      if (my_task == master_task) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'Forcing data year =', fyear
         write (nu_diag,*) 'Atmospheric data files:'
         write (nu_diag,'(3a)') trim(fsw_file)
         write (nu_diag,'(3a)') trim(flw_file)
         write (nu_diag,'(3a)') trim(rain_file)
         write (nu_diag,'(3a)') trim(uwind_file)
         write (nu_diag,'(3a)') trim(vwind_file)
         write (nu_diag,'(3a)') trim(tair_file)
         write (nu_diag,'(3a)') trim(humid_file)
         write (nu_diag,'(3a)') trim(rhoa_file)
      endif                     ! master_task

      end subroutine ncar_files

!=======================================================================

      subroutine ncar_data

      use ice_flux, only: fsw, fsnow, Tair, uatm, vatm, rhoa, Qa

      integer (kind=int_kind) :: &
          ixm,ixx,ixp , & ! record numbers for neighboring months
          recnum      , & ! record number
          maxrec      , & ! maximum record number
          recslot     , & ! spline slot for current record
          dataloc     , & ! = 1 for data located in middle of time interval
                          ! = 2 for date located at end of time interval
          midmonth        ! middle day of month

      real (kind=dbl_kind) :: &
          secday, &           ! number of seconds in day
          sec6hr              ! number of seconds in 6 hours

      logical (kind=log_kind) :: readm, read6

      character(len=*), parameter :: subname = '(ncar_data)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      call icepack_query_parameters(secday_out=secday)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

    !-------------------------------------------------------------------
    ! monthly data
    !
    ! Assume that monthly data values are located in the middle of the
    ! month.
    !-------------------------------------------------------------------

      midmonth = 15  ! data is given on 15th of every month
!      midmonth = fix(p5 * real(daymo(mmonth)))  ! exact middle

      ! Compute record numbers for surrounding months
      maxrec = 12
      ixm  = mod(mmonth+maxrec-2,maxrec) + 1
      ixp  = mod(mmonth,         maxrec) + 1
      if (mday >= midmonth) ixm = -99  ! other two points will be used
      if (mday <  midmonth) ixp = -99

      ! Determine whether interpolation will use values 1:2 or 2:3
      ! recslot = 2 means we use values 1:2, with the current value (2)
      !  in the second slot
      ! recslot = 1 means we use values 2:3, with the current value (2)
      !  in the first slot
      recslot = 1                             ! latter half of month
      if (mday < midmonth) recslot = 2        ! first half of month

      ! Find interpolation coefficients
      call interp_coeff_monthly (recslot)

      ! Read 2 monthly values
      readm = .false.
      if (istep==1 .or. (mday==midmonth .and. msec==0)) readm = .true.

      if (trim(atm_data_format) == 'bin') then
         call read_data (readm, 0, fyear, ixm, mmonth, ixp, &
                         maxrec, fsw_file, fsw_data, &
                         field_loc_center, field_type_scalar)
         call read_data (readm, 0, fyear, ixm, mmonth, ixp, &
                         maxrec, flw_file, cldf_data, &
                         field_loc_center, field_type_scalar)
         call read_data (readm, 0, fyear, ixm, mmonth, ixp, &
                         maxrec, rain_file, fsnow_data, &
                         field_loc_center, field_type_scalar)
      else
         call abort_ice (error_message=subname//'nonbinary atm_data_format unavailable', &
            file=__FILE__, line=__LINE__)
!        The routine exists, for example:
!         call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, &
!                            maxrec, fsw_file, 'fsw', fsw_data, &
!                            field_loc_center, field_type_scalar)
!         call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, &
!                            maxrec, flw_file, 'cldf',cldf_data, &
!                            field_loc_center, field_type_scalar)
!         call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, &
!                            maxrec, rain_file,'prec',fsnow_data, &
!                            field_loc_center, field_type_scalar)
      endif

      ! Interpolate to current time step
      call interpolate_data (fsw_data,   fsw)
      call interpolate_data (cldf_data,  cldf)
      call interpolate_data (fsnow_data, fsnow)

    !-------------------------------------------------------------------
    ! 6-hourly data
    !
    ! Assume that the 6-hourly value is located at the end of the
    !  6-hour period.  This is the convention for NCEP reanalysis data.
    !  E.g. record 1 gives conditions at 6 am GMT on 1 January.
    !-------------------------------------------------------------------

      dataloc = 2               ! data located at end of interval
      sec6hr = secday/c4        ! seconds in 6 hours
      maxrec = 1460             ! 365*4

      ! current record number
      recnum = 4*int(yday) - 3 + int(real(msec,kind=dbl_kind)/sec6hr)

      ! Compute record numbers for surrounding data

      ixm = mod(recnum+maxrec-2,maxrec) + 1
      ixx = mod(recnum-1,       maxrec) + 1
!      ixp = mod(recnum,         maxrec) + 1

      ! Compute interpolation coefficients
      ! If data is located at the end of the time interval, then the
      !  data value for the current record always goes in slot 2.

      recslot = 2
      ixp = -99
      call interp_coeff (recnum, recslot, sec6hr, dataloc)

      ! Read
      read6 = .false.
      if (istep==1 .or. oldrecnum /= recnum) read6 = .true.

      if (trim(atm_data_format) == 'bin') then
         call read_data (read6, 0, fyear, ixm, ixx, ixp, &
                         maxrec, tair_file, Tair_data, &
                         field_loc_center, field_type_scalar)
         call read_data (read6, 0, fyear, ixm, ixx, ixp, &
                         maxrec, uwind_file, uatm_data, &
                         field_loc_center, field_type_vector)
         call read_data (read6, 0, fyear, ixm, ixx, ixp, &
                         maxrec, vwind_file, vatm_data, &
                         field_loc_center, field_type_vector)
         call read_data (read6, 0, fyear, ixm, ixx, ixp, &
                         maxrec, rhoa_file, rhoa_data, &
                         field_loc_center, field_type_scalar)
         call read_data (read6, 0, fyear, ixm, ixx, ixp, &
                         maxrec, humid_file, Qa_data, &
                         field_loc_center, field_type_scalar)
      else
         call abort_ice (error_message=subname//'nonbinary atm_data_format unavailable', &
            file=__FILE__, line=__LINE__)
      endif

      ! Interpolate
      call interpolate_data (Tair_data, Tair)
      call interpolate_data (uatm_data, uatm)
      call interpolate_data (vatm_data, vatm)
      call interpolate_data (rhoa_data, rhoa)
      call interpolate_data (Qa_data,   Qa)

      ! Save record number for next time step
      oldrecnum = recnum

      end subroutine ncar_data

!=======================================================================

      subroutine JRA55_files(yr)

      ! find the JRA55 files:
      ! This subroutine finds the JRA55 atm forcing files based on settings
      ! in atm_data_type and atm_data_dir.  Because the filenames are not
      ! entirely consistent, we need a flexible method.
      !
      ! atm_data_type could be JRA55 or JRA55do with/without _grid appended
      ! atm_data_dir could contain JRA55 or JRA55do or not
      ! actual files could have grid in name in two location or not at all
      !
      ! The files will generally be of the format
      !    $atm_data_type/[JRA55,JRA55do,'']/8XDAILY/[JRA55,JRA55do][_$grid,'']_03hr_forcing[_$grid,'']_$year.nc
      ! The options defined by cnt try several versions of paths/filenames
      ! As a user,
      !    atm_data_type should be set to JRA55, JRA55do, JRA55_xxx, or JRA55do_xxx
      !       where xxx can be any set of characters.  The _xxx if included will be ignored.
      !       Historically, these were set to JRA55_gx1 and so forth but the _gx1 is no longer needed
      !       but this is still allowed for backwards compatibility.  atm_data_type_prefix
      !       is atm_data_type with _ and everything after _ removed.
      !    atm_data_dir should be set to ${CICE_DATA_root}/forcing/$grid/[JRA55,JRA55do,'']
      !       The [JRA55,JRA55do] at the end of the atm_data_dir is optional to provide backwards
      !       compatibility and if not included, will be appended automaticaly using
      !       the atm_data_type_prefix value.  The grid is typically gx1, gx3, tx1, or similar.
      ! In general, we recommend using the following format
      !    atm_data_type = [JRA55,JRA55do]
      !    atm_data_dir = ${CICE_DATA_root}/forcing/$grid

      integer (kind=int_kind), intent(in) :: &
           yr         ! current forcing year

      ! local variables
      character(len=16) :: &
           grd        ! gx3, gx1, tx1

      character(len=64) :: &
           atm_data_type_prefix  ! atm_data_type prefix

      integer (kind=int_kind) :: &
           cnt    , & ! search for files
           strind     ! string index

      logical :: &
           exists     ! file existance

      character(len=*), parameter :: subname = '(JRA55_files)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      ! this could be JRA55[do] or JRA55[do]_grid, drop the _grid if set
      atm_data_type_prefix = trim(atm_data_type)
      strind = index(trim(atm_data_type),'_')
      if (strind > 0) then
         atm_data_type_prefix = atm_data_type(1:strind-1)
      endif

      ! check for grid version using fortran INDEX intrinsic
      if (index(trim(atm_data_dir),'gx1') > 0) then
         grd = 'gx1'
      else if (index(trim(atm_data_dir),'gx3') > 0) then
         grd = 'gx3'
      else if (index(trim(atm_data_dir),'tx1') > 0) then
         grd = 'tx1'
      else
         call abort_ice(error_message=subname//' unknown grid type')
      endif

      ! cnt represents the possible file format options and steps thru them until one is found
      exists = .false.
      cnt = 1
      do while (.not.exists .and. cnt <= 6)

         if (cnt == 1) uwind_file = trim(atm_data_dir)//'/'//trim(atm_data_type_prefix)//     &
                                    '/8XDAILY/'//trim(atm_data_type_prefix)//'_'//trim(grd)// &
                                    '_03hr_forcing'//trim(atm_data_version)//'_2005.nc'

         if (cnt == 2) uwind_file = trim(atm_data_dir)//'/'//trim(atm_data_type_prefix)//                  &
                                    '/8XDAILY/'//trim(atm_data_type_prefix)//'_03hr_forcing_'//trim(grd)// &
                                    trim(atm_data_version)//'_2005.nc'

         if (cnt == 3) uwind_file = trim(atm_data_dir)//'/'//trim(atm_data_type_prefix)// &
                                    '/8XDAILY/'//trim(atm_data_type_prefix)//             &
                                    '_03hr_forcing'//trim(atm_data_version)//'_2005.nc'

         if (cnt == 4) uwind_file = trim(atm_data_dir)//                                      &
                                    '/8XDAILY/'//trim(atm_data_type_prefix)//'_'//trim(grd)// &
                                    '_03hr_forcing'//trim(atm_data_version)//'_2005.nc'

         if (cnt == 5) uwind_file = trim(atm_data_dir)//                                                   &
                                    '/8XDAILY/'//trim(atm_data_type_prefix)//'_03hr_forcing_'//trim(grd)// &
                                    trim(atm_data_version)//'_2005.nc'

         if (cnt == 6) uwind_file = trim(atm_data_dir)//                                  &
                                    '/8XDAILY/'//trim(atm_data_type_prefix)//             &
                                    '_03hr_forcing'//trim(atm_data_version)//'_2005.nc'


         call file_year(uwind_file,yr)
         INQUIRE(FILE=uwind_file,EXIST=exists)

         if (debug_forcing .and. (my_task == master_task)) then
            write(nu_diag,*) subname,cnt,exists,trim(uwind_file)
         endif

         cnt = cnt + 1
      enddo

      if (.not.exists) then
         write(nu_diag,*) subname,' atm_data_dir = ',trim(atm_data_dir)
         write(nu_diag,*) subname,' atm_data_type_prefix = ',trim(atm_data_type_prefix)
         write(nu_diag,*) subname,' atm_data_version = ',trim(atm_data_version)
         call abort_ice(error_message=subname//' could not find forcing file')
      endif

      if (my_task == master_task) then
         write (nu_diag,'(2a)') ' '
         write (nu_diag,'(2a)') subname,'Atmospheric data files:'
         write (nu_diag,'(2a)') subname,trim(uwind_file)
      endif

    end subroutine JRA55_files

!=======================================================================

      subroutine JRA55_data

      use ice_blocks, only: block, get_block
      use ice_global_reductions, only: global_minval, global_maxval
      use ice_domain, only: nblocks, distrb_info
      use ice_flux, only: fsnow, Tair, uatm, vatm, Qa, fsw, flw
      use ice_grid, only: hm, tmask, umask
      use ice_state, only: aice
      use ice_calendar, only: days_per_year

      integer (kind=int_kind) :: &
          ncid        , & ! netcdf file id
          i, j, n1    , &
          lfyear      , & ! local year value
          recnum      , & ! record number
          maxrec      , & ! maximum record number
          iblk            ! block index

      integer (kind=int_kind), save :: &
          frec_info(2,2) = -99    ! remember prior values to reduce reading
                                  ! first dim is yr, recnum
                                  ! second dim is data1 data2

      real (kind=dbl_kind) :: &
          sec3hr          , & ! number of seconds in 3 hours
          secday          , & ! number of seconds in day
          eps, tt         , & ! for interpolation coefficients
          Tffresh         , &
          vmin, vmax

      character(len=64) :: fieldname !netcdf field name
      character (char_len_long) :: uwind_file_old
      character(len=*), parameter :: subname = '(JRA55_data)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      call icepack_query_parameters(Tffresh_out=Tffresh)
      call icepack_query_parameters(secday_out=secday)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      sec3hr = secday/c8        ! seconds in 3 hours
      maxrec = days_per_year * 8

      if (local_debug .and. my_task == master_task) then
         write(nu_diag,*) subname,'fdbg dpy, maxrec = ',days_per_year,maxrec
      endif

      !-------------------------------------------------------------------
      ! 3-hourly data
      ! states are instantaneous, 1st record is 00z Jan 1
      ! fluxes are 3 hour averages, 1st record is 00z-03z Jan 1
      ! interpolate states, do not interpolate fluxes
      !-------------------------------------------------------------------
      ! File is NETCDF with winds in NORTH and EAST direction
      ! file variable names are:
      ! glbrad   (shortwave W/m^2), 3 hr average
      ! dlwsfc   (longwave W/m^2), 3 hr average
      ! wndewd   (eastward wind m/s), instantaneous
      ! wndnwd   (northward wind m/s), instantaneous
      ! airtmp   (air temperature K), instantaneous
      ! spchmd   (specific humidity kg/kg), instantaneous
      ! ttlpcp   (precipitation kg/m s-1), 3 hr average
      !-------------------------------------------------------------------

      uwind_file_old = uwind_file
      if (uwind_file /= uwind_file_old .and. my_task == master_task) then
         write(nu_diag,'(2a)') subname,' reading forcing file = ',trim(uwind_file)
      endif

      call ice_open_nc(uwind_file,ncid)

      do n1 = 1, 2

         lfyear = fyear
         call file_year(uwind_file,lfyear)
         if (n1 == 1) then
            recnum = 8*int(yday) - 7 + int(real(msec,kind=dbl_kind)/sec3hr)
            if (my_task == master_task .and. (recnum <= 2 .or. recnum >= maxrec-1)) then
               write(nu_diag,'(3a)') subname,' reading forcing file 1st ts = ',trim(uwind_file)
            endif
         elseif (n1 == 2) then
            recnum = 8*int(yday) - 7 + int(real(msec,kind=dbl_kind)/sec3hr) + 1
            if (recnum > maxrec) then
               lfyear = fyear + 1  ! next year
               if (lfyear > fyear_final) lfyear = fyear_init
               recnum = 1
               call file_year(uwind_file,lfyear)
               if (my_task == master_task) then
                  write(nu_diag,'(3a)') subname,' reading forcing file 2nd ts = ',trim(uwind_file)
               endif
               call ice_close_nc(ncid)
               call ice_open_nc(uwind_file,ncid)
            endif
         endif

         if (local_debug .and. my_task == master_task) then
            write(nu_diag,*) subname,'fdbg read recnum = ',recnum,n1
         endif

         ! to reduce reading, check whether it's the same data as last read

         if (lfyear /= frec_info(1,n1) .or. recnum /= frec_info(2,n1)) then

            ! check whether we can copy values from 2 to 1, should be faster than reading
            ! can only do this from 2 to 1 or 1 to 2 without setting up a temporary
            ! it's more likely that the values from data2 when time advances are needed in data1
            ! compare n1=1 year/record with data from last timestep at n1=2

            if (n1 == 1 .and. lfyear == frec_info(1,2) .and. recnum == frec_info(2,2)) then
                Tair_data(:,:,1,:) =  Tair_data(:,:,2,:)
                uatm_data(:,:,1,:) =  uatm_data(:,:,2,:)
                vatm_data(:,:,1,:) =  vatm_data(:,:,2,:)
                  Qa_data(:,:,1,:) =    Qa_data(:,:,2,:)
                 fsw_data(:,:,1,:) =   fsw_data(:,:,2,:)
                 flw_data(:,:,1,:) =   flw_data(:,:,2,:)
               fsnow_data(:,:,1,:) = fsnow_data(:,:,2,:)
            else

               fieldname = 'airtmp'
               call ice_read_nc(ncid,recnum,fieldname,Tair_data(:,:,n1,:),local_debug, &
                    field_loc=field_loc_center, &
                    field_type=field_type_scalar)

               fieldname = 'wndewd'
               call ice_read_nc(ncid,recnum,fieldname,uatm_data(:,:,n1,:),local_debug, &
                    field_loc=field_loc_center, &
                    field_type=field_type_scalar)

               fieldname = 'wndnwd'
               call ice_read_nc(ncid,recnum,fieldname,vatm_data(:,:,n1,:),local_debug, &
                    field_loc=field_loc_center, &
                    field_type=field_type_scalar)

               fieldname = 'spchmd'
               call ice_read_nc(ncid,recnum,fieldname,Qa_data(:,:,n1,:),local_debug, &
                    field_loc=field_loc_center, &
                    field_type=field_type_scalar)

               fieldname = 'glbrad'
               call ice_read_nc(ncid,recnum,fieldname,fsw_data(:,:,n1,:),local_debug, &
                    field_loc=field_loc_center, &
                    field_type=field_type_scalar)

               fieldname = 'dlwsfc'
               call ice_read_nc(ncid,recnum,fieldname,flw_data(:,:,n1,:),local_debug, &
                    field_loc=field_loc_center, &
                    field_type=field_type_scalar)

               fieldname = 'ttlpcp'
               call ice_read_nc(ncid,recnum,fieldname,fsnow_data(:,:,n1,:),local_debug, &
                    field_loc=field_loc_center, &
                    field_type=field_type_scalar)
            endif  ! copy data from n1=2 from last timestep to n1=1
         endif  ! input data is same as last timestep

         frec_info(1,n1) = lfyear
         frec_info(2,n1) = recnum

      enddo  ! n1

      call ice_close_nc(ncid)

      ! reset uwind_file to original year
      call file_year(uwind_file,fyear)

      ! Compute interpolation coefficients
      eps = 1.0e-6
      tt = real(mod(msec,nint(sec3hr)),kind=dbl_kind)
      c2intp = tt / sec3hr
      if (c2intp < c0 .and. c2intp > c0-eps) c2intp = c0
      if (c2intp > c1 .and. c2intp < c1+eps) c2intp = c1
      c1intp = 1.0_dbl_kind - c2intp
      if (c2intp < c0 .or. c2intp > c1) then
         write(nu_diag,*) subname,' ERROR: c2intp = ',c2intp
         call abort_ice (error_message=subname//' ERROR: c2intp out of range', &
            file=__FILE__, line=__LINE__)
      endif
      if (local_debug .and. my_task == master_task) then
         write(nu_diag,*) subname,'fdbg c12intp = ',c1intp,c2intp
      endif

      ! Interpolate
      call interpolate_data (Tair_data, Tair)
      call interpolate_data (uatm_data, uatm)
      call interpolate_data (vatm_data, vatm)
      call interpolate_data (Qa_data, Qa)
      ! use 3 hr average for heat flux and precip fields, no interpolation
!      call interpolate_data (fsw_data, fsw)
!      call interpolate_data (flw_data, flw)
!      call interpolate_data (fsnow_data, fsnow)
      fsw(:,:,:) = fsw_data(:,:,1,:)
      flw(:,:,:) = flw_data(:,:,1,:)
      fsnow(:,:,:) = fsnow_data(:,:,1,:)

      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
        ! limit summer Tair values where ice is present
        do j = 1, ny_block
          do i = 1, nx_block
            if (aice(i,j,iblk) > p1) Tair(i,j,iblk) = min(Tair(i,j,iblk), Tffresh+p1)
          enddo
        enddo

        do j = 1, ny_block
          do i = 1, nx_block
            Qa  (i,j,iblk) = Qa  (i,j,iblk) * hm(i,j,iblk)
            Tair(i,j,iblk) = Tair(i,j,iblk) * hm(i,j,iblk)
            uatm(i,j,iblk) = uatm(i,j,iblk) * hm(i,j,iblk)
            vatm(i,j,iblk) = vatm(i,j,iblk) * hm(i,j,iblk)
            fsw (i,j,iblk) = fsw (i,j,iblk) * hm(i,j,iblk)
            flw (i,j,iblk) = flw (i,j,iblk) * hm(i,j,iblk)
            fsnow(i,j,iblk) = fsnow (i,j,iblk) * hm(i,j,iblk)
          enddo
        enddo

      enddo  ! iblk
      !$OMP END PARALLEL DO

      if (debug_forcing .or. local_debug) then
         if (my_task.eq.master_task) write (nu_diag,*) subname,'fdbg JRA55_bulk_data'
         vmin = global_minval(fsw,distrb_info,tmask)
         vmax = global_maxval(fsw,distrb_info,tmask)
         if (my_task.eq.master_task) write (nu_diag,*) subname,'fdbg fsw',vmin,vmax
         vmin = global_minval(flw,distrb_info,tmask)
         vmax = global_maxval(flw,distrb_info,tmask)
         if (my_task.eq.master_task) write (nu_diag,*) subname,'fdbg flw',vmin,vmax
         vmin =global_minval(fsnow,distrb_info,tmask)
         vmax =global_maxval(fsnow,distrb_info,tmask)
         if (my_task.eq.master_task) write (nu_diag,*) subname,'fdbg fsnow',vmin,vmax
         vmin = global_minval(Tair,distrb_info,tmask)
         vmax = global_maxval(Tair,distrb_info,tmask)
         if (my_task.eq.master_task) write (nu_diag,*) subname,'fdbg Tair',vmin,vmax
         vmin = global_minval(uatm,distrb_info,umask)
         vmax = global_maxval(uatm,distrb_info,umask)
         if (my_task.eq.master_task) write (nu_diag,*) subname,'fdbg uatm',vmin,vmax
         vmin = global_minval(vatm,distrb_info,umask)
         vmax = global_maxval(vatm,distrb_info,umask)
         if (my_task.eq.master_task) write (nu_diag,*) subname,'fdbg vatm',vmin,vmax
         vmin = global_minval(Qa,distrb_info,tmask)
         vmax = global_maxval(Qa,distrb_info,tmask)
         if (my_task.eq.master_task) write (nu_diag,*) subname,'fdbg Qa',vmin,vmax
      endif                   ! debug_forcing

      end subroutine JRA55_data

!=======================================================================
!
! AOMIP shortwave forcing
! standard calculation using solar declination angle
! then shortwave is reduced using a function of cloud fraction

      subroutine compute_shortwave(nx_block,  ny_block, &
                                   ilo, ihi, jlo, jhi, &
                                   TLON, TLAT, hm, Qa, cldf, fsw)

!---!-------------------------------------------------------------------
!---!-------------------------------------------------------------------

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ilo,ihi,jlo,jhi       ! beginning and end of physical domain

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) :: &
         TLON, TLAT     , & ! longitude, latitude
         Qa             , & ! specific humidity
         cldf           , & ! cloud fraction
         hm                 ! land mask

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(inout) :: &
         fsw                ! shortwave

      real (kind=dbl_kind) :: &
         hour_angle, &
         solar_time, &
         declin    , &
         cosZ      , &
         e, d      , &
         sw0       , &
         secday    , &
         pi        , &
         lontmp    , &
         deg2rad

      integer (kind=int_kind) :: &
         i, j

      character(len=*), parameter :: subname = '(compute_shortwave)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      call icepack_query_parameters(secday_out=secday, pi_out=pi)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      do j=jlo,jhi
       do i=ilo,ihi
        deg2rad = pi/c180
!       solar_time = mod(real(msec,kind=dbl_kind),secday)/c3600 &
!                  + c12*sin(p5*TLON(i,j))

!       Convert longitude to range of -180 to 180 for LST calculation

        lontmp = mod(TLON(i,j)/deg2rad,c360)
        if (lontmp .gt. c180) lontmp = lontmp - c360
        if (lontmp .lt. -c180) lontmp = lontmp + c360

        solar_time = mod(real(msec,kind=dbl_kind),secday)/c3600 &
                   + lontmp/c15
        if (solar_time .ge. 24._dbl_kind) solar_time = solar_time - 24._dbl_kind
        hour_angle = (c12 - solar_time)*pi/c12
        declin = 23.44_dbl_kind*cos((172._dbl_kind-yday) &
                 * c2*pi/c365)*deg2rad     ! use dayyr instead of c365???
        cosZ = sin(TLAT(i,j))*sin(declin) &
             + cos(TLAT(i,j))*cos(declin)*cos(hour_angle)
        cosZ = max(cosZ,c0)
        e = 1.e5*Qa(i,j)/(0.622_dbl_kind + 0.378_dbl_kind*Qa(i,j))
        d = (cosZ+2.7_dbl_kind)*e*1.e-5_dbl_kind+1.085_dbl_kind*cosZ+p1
        sw0 = 1353._dbl_kind*cosZ**2/d
        sw0 = max(sw0,c0)

        ! total downward shortwave for cice
        Fsw(i,j) = sw0*(c1-p6*cldf(i,j)**3)
        Fsw(i,j) = Fsw(i,j)*hm(i,j)
       enddo
      enddo

      end subroutine compute_shortwave

!=======================================================================
!
! prevents humidity from being super-saturated

      subroutine Qa_fixLY(nx_block, ny_block, Tair, Qa)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block ! block dimensions

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) :: &
         Tair               ! air temperature

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(inout) :: &
         Qa                 ! specific humidity

      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         worka

      real (kind=dbl_kind) :: &
         Tffresh, puny

      character(len=*), parameter :: subname = '(Qa_fixLY)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      call icepack_query_parameters(Tffresh_out=Tffresh, puny_out=puny)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      worka = Tair - Tffresh
      worka = c2 + (0.7859_dbl_kind + 0.03477_dbl_kind*worka) &
                     /(c1 + 0.00412_dbl_kind*worka) & ! 2+ converts ea mb -> Pa
                + 0.00422_dbl_kind*worka              ! for ice
      ! vapor pressure
      worka = (c10**worka)      ! saturated
      worka = max(worka,puny)   ! puny over land to prevent division by zero
      ! specific humidity
      worka = 0.622_dbl_kind*worka/(1.e5_dbl_kind-0.378_dbl_kind*worka)

      Qa = min(Qa, worka)

      end subroutine Qa_fixLY

!=======================================================================
! HadGEM or HadGAM atmospheric forcing
!=======================================================================

      subroutine hadgem_files (yr)

! Construct filenames based on selected model options
!
! Note: The year number in these filenames does not matter, because
!       subroutine file_year will insert the correct year.
!
! author: Alison McLaren, Met Office

      integer (kind=int_kind), intent(in) :: &
           yr                   ! current forcing year

      integer (kind=int_kind) :: &
           n           ! thickness category index

      logical (kind=log_kind) :: calc_strair, calc_Tsfc

      character(len=*), parameter :: subname = '(hadgem_files)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      call icepack_query_parameters(calc_strair_out=calc_strair, &
           calc_Tsfc_out=calc_Tsfc)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      ! -----------------------------------------------------------
      ! Rainfall and snowfall
      ! -----------------------------------------------------------

      snow_file = &
           trim(atm_data_dir)//'/MONTHLY/snowfall.1996.nc'
           call file_year(snow_file,yr)

      rain_file = &
           trim(atm_data_dir)//'/MONTHLY/rainfall.1996.nc'
           call file_year(rain_file,yr)

      if (my_task == master_task) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'Atmospheric data files:'
         write (nu_diag,*) trim(rain_file)
         write (nu_diag,*) trim(snow_file)
      endif

      if (calc_strair) then

         ! --------------------------------------------------------
         ! Wind velocity
         ! --------------------------------------------------------

         uwind_file = &
           trim(atm_data_dir)//'/MONTHLY/u_10.1996.nc'
           call file_year(uwind_file,yr)

         vwind_file = &
           trim(atm_data_dir)//'/MONTHLY/v_10.1996.nc'
           call file_year(vwind_file,yr)

         if (my_task == master_task) then
            write (nu_diag,*) trim(uwind_file)
            write (nu_diag,*) trim(vwind_file)
         endif

      else

         ! --------------------------------------------------------
         ! Wind stress
         ! --------------------------------------------------------

         strax_file = &
              trim(atm_data_dir)//'/MONTHLY/taux.1996.nc'
         call file_year(strax_file,yr)

         stray_file = &
              trim(atm_data_dir)//'/MONTHLY/tauy.1996.nc'
         call file_year(stray_file,yr)

         if (my_task == master_task) then
            write (nu_diag,*) trim(strax_file)
            write (nu_diag,*) trim(stray_file)
         endif

         if (calc_Tsfc .or. oceanmixed_ice) then

            ! --------------------------------------------------
            ! Wind speed
            ! --------------------------------------------------

            wind_file = &
               trim(atm_data_dir)//'/MONTHLY/wind_10.1996.nc'
            call file_year(wind_file,yr)

            if (my_task == master_task) then
               write (nu_diag,*) trim(wind_file)
            endif

         endif   ! calc_Tsfc or oceanmixed_ice

      endif  ! calc_strair

      ! --------------------------------------------------------------
      ! Atmosphere properties.  Even if these fields are not
      ! being used to force the ice (i.e. calc_Tsfc=.false.), they
      ! are still needed to generate forcing for mixed layer model or
      ! to calculate wind stress
      ! --------------------------------------------------------------

       if (calc_Tsfc .or. oceanmixed_ice .or. calc_strair) then

         fsw_file = &
           trim(atm_data_dir)//'/MONTHLY/SW_incoming.1996.nc'
           call file_year(fsw_file,yr)

         flw_file = &
           trim(atm_data_dir)//'/MONTHLY/LW_incoming.1996.nc'
           call file_year(flw_file,yr)

         tair_file = &
           trim(atm_data_dir)//'/MONTHLY/t_10.1996.nc'
           call file_year(tair_file,yr)

         humid_file = &
           trim(atm_data_dir)//'/MONTHLY/q_10.1996.nc'
           call file_year(humid_file,yr)

         rhoa_file = &
           trim(atm_data_dir)//'/MONTHLY/rho_10.1996.nc'
           call file_year(rhoa_file,yr)

         if (my_task == master_task) then
            write (nu_diag,*) trim(fsw_file)
            write (nu_diag,*) trim(flw_file)
            write (nu_diag,*) trim(tair_file)
            write (nu_diag,*) trim(humid_file)
            write (nu_diag,*) trim(rhoa_file)
         endif                     ! master_task

      endif ! calc_Tsfc or oceanmixed_ice  or calc_strair

      if (.not. calc_Tsfc) then

         ! ------------------------------------------------------
         ! Sublimation, topmelt and botmelt
         ! ------------------------------------------------------

         do n = 1, ncat

            ! 'topmelt' = fsurf - fcondtop.
            write(topmelt_file(n), '(a,i1,a)')  &
              trim(atm_data_dir)//'/MONTHLY/topmeltn',n,'.1996.nc'
              call file_year(topmelt_file(n),yr)

            ! 'botmelt' = fcondtop.
            write(botmelt_file(n), '(a,i1,a)')  &
              trim(atm_data_dir)//'/MONTHLY/botmeltn',n,'.1996.nc'
              call file_year(botmelt_file(n),yr)

         enddo

         ! 'sublim' = - flat / Lsub.
         sublim_file = &
           trim(atm_data_dir)//'/MONTHLY/sublim.1996.nc'
           call file_year(sublim_file,yr)

         if (my_task == master_task) then
            do n = 1, ncat
               write (nu_diag,*) trim(topmelt_file(n))
               write (nu_diag,*) trim(botmelt_file(n))
            enddo
            write (nu_diag,*) trim(sublim_file)

         endif

      endif  ! .not. calc_Tsfc

      end subroutine hadgem_files

!=======================================================================

! read HadGEM or HadGAM atmospheric data

      subroutine hadgem_data

! authors: Alison McLaren, Met Office

      use ice_domain, only: nblocks
      use ice_flux, only: fsnow, frain, uatm, vatm, strax, stray, wind, &
          fsw, flw, Tair, rhoa, Qa, fcondtopn_f, fsurfn_f, flatn_f

      integer (kind=int_kind) :: &
          i, j        , & ! horizontal indices
          n           , & ! thickness category index
          iblk        , & ! block index
          ixm,ixp     , & ! record numbers for neighboring months
          maxrec      , & ! maximum record number
          recslot     , & ! spline slot for current record
          midmonth        ! middle day of month

      logical (kind=log_kind) :: readm

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
            topmelt, & ! temporary fields
            botmelt, &
            sublim

      character (char_len) :: &
            fieldname    ! field name in netcdf file

      real (kind=dbl_kind) :: &
            Lsub

      logical (kind=log_kind) :: &
            calc_strair, &
            calc_Tsfc

      character(len=*), parameter :: subname = '(hadgem_data)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      call icepack_query_parameters(Lsub_out=Lsub)
      call icepack_query_parameters(calc_strair_out=calc_strair, &
           calc_Tsfc_out=calc_Tsfc)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

    !-------------------------------------------------------------------
    ! monthly data
    !
    ! Assume that monthly data values are located in the middle of the
    ! month.
    !-------------------------------------------------------------------

      midmonth = 15  ! data is given on 15th of every month
!      midmonth = fix(p5 * real(daymo(mmonth)))  ! exact middle

      ! Compute record numbers for surrounding months
      maxrec = 12
      ixm  = mod(mmonth+maxrec-2,maxrec) + 1
      ixp  = mod(mmonth,         maxrec) + 1
      if (mday >= midmonth) ixm = -99  ! other two points will be used
      if (mday <  midmonth) ixp = -99

      ! Determine whether interpolation will use values 1:2 or 2:3
      ! recslot = 2 means we use values 1:2, with the current value (2)
      !  in the second slot
      ! recslot = 1 means we use values 2:3, with the current value (2)
      !  in the first slot
      recslot = 1                             ! latter half of month
      if (mday < midmonth) recslot = 2        ! first half of month

      ! Find interpolation coefficients
      call interp_coeff_monthly (recslot)

      ! Read 2 monthly values
      readm = .false.
      if (istep==1 .or. (mday==midmonth .and. msec==0)) readm = .true.

      ! -----------------------------------------------------------
      ! Rainfall and snowfall
      ! -----------------------------------------------------------

      fieldname='rainfall'
      call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, &
                      maxrec, rain_file, fieldname, frain_data, &
                      field_loc_center, field_type_scalar)
      fieldname='snowfall'
      call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, &
                      maxrec, snow_file, fieldname, fsnow_data, &
                      field_loc_center, field_type_scalar)

      ! Interpolate to current time step
      call interpolate_data (fsnow_data, fsnow)
      call interpolate_data (frain_data, frain)

      if (calc_strair) then

         ! --------------------------------------------------------
         ! Wind velocity
         ! --------------------------------------------------------

         fieldname='u_10'
         call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, &
                      maxrec, uwind_file, fieldname, uatm_data, &
                      field_loc_center, field_type_vector)
         fieldname='v_10'
         call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, &
                      maxrec, vwind_file, fieldname, vatm_data, &
                      field_loc_center, field_type_vector)

         ! Interpolate to current time step
         call interpolate_data (uatm_data, uatm)
         call interpolate_data (vatm_data, vatm)

      else

         ! --------------------------------------------------------
         ! Wind stress
         ! --------------------------------------------------------

         fieldname='taux'
         call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, &
                      maxrec, strax_file, fieldname, strax_data, &
                      field_loc_center, field_type_vector)
         fieldname='tauy'
         call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, &
                      maxrec, stray_file, fieldname, stray_data, &
                      field_loc_center, field_type_vector)

         ! Interpolate to current time step
         call interpolate_data (strax_data, strax)
         call interpolate_data (stray_data, stray)

         if (calc_Tsfc .or. oceanmixed_ice) then

            ! --------------------------------------------------
            ! Wind speed
            ! --------------------------------------------------

            fieldname='wind_10'
            call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, &
                      maxrec, wind_file, fieldname, wind_data, &
                      field_loc_center, field_type_scalar)

            ! Interpolate to current time step
            call interpolate_data (wind_data, wind)

         endif   ! calc_Tsfc or oceanmixed_ice

      endif      ! calc_strair

      ! -----------------------------------------------------------
      ! SW incoming, LW incoming, air temperature, density and
      ! humidity at 10m.
      !
      ! Even if these fields are not being used to force the ice
      ! (i.e. calc_Tsfc=.false.), they are still needed to generate
      ! forcing for mixed layer model or to calculate wind stress
      ! -----------------------------------------------------------

      if (calc_Tsfc .or. oceanmixed_ice .or. calc_strair) then

         fieldname='SW_incoming'
         call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, &
                      maxrec, fsw_file, fieldname, fsw_data, &
                      field_loc_center, field_type_scalar)
         fieldname='LW_incoming'
         call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, &
                      maxrec, flw_file, fieldname, flw_data, &
                      field_loc_center, field_type_scalar)
         fieldname='t_10'
         call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, &
                      maxrec, tair_file, fieldname, Tair_data, &
                      field_loc_center, field_type_scalar)
         fieldname='rho_10'
         call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, &
                      maxrec, rhoa_file, fieldname, rhoa_data, &
                      field_loc_center, field_type_scalar)
         fieldname='q_10'
         call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, &
                      maxrec, humid_file, fieldname, Qa_data, &
                      field_loc_center, field_type_scalar)

         ! Interpolate onto current timestep

         call interpolate_data (fsw_data,   fsw)
         call interpolate_data (flw_data,  flw)
         call interpolate_data (Tair_data, Tair)
         call interpolate_data (rhoa_data, rhoa)
         call interpolate_data (Qa_data,   Qa)

      endif       ! calc_Tsfc or oceanmixed_ice or calc_strair

      if (.not. calc_Tsfc) then

         ! ------------------------------------------------------
         ! Sublimation, topmelt and botmelt
         ! ------------------------------------------------------

         fieldname='sublim'
         call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, &
                      maxrec, sublim_file, fieldname, sublim_data, &
                      field_loc_center, field_type_scalar)

         ! Interpolate to current time step
         call interpolate_data (sublim_data, sublim)

         do n = 1, ncat
            write(fieldname, '(a,i1)') 'topmeltn',n
            call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, &
              maxrec, topmelt_file(n), fieldname, topmelt_data(:,:,:,:,n), &
              field_loc_center, field_type_scalar)

            write(fieldname, '(a,i1)') 'botmeltn',n
            call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, &
              maxrec, botmelt_file(n), fieldname, botmelt_data(:,:,:,:,n), &
              field_loc_center, field_type_scalar)

            call interpolate_data (topmelt_data(:,:,:,:,n), topmelt)
            call interpolate_data (botmelt_data(:,:,:,:,n), botmelt)

            !--------------------------------------------------------
            ! Convert from UM variables to CICE variables
            !  topmelt = fsurf - fcondtop
            !  botmelt = fcondtop  (as zero layer)
            !
            ! Convert UM sublimation data into CICE LH flux
            ! (sublim = - flatn / Lsub) and have same value for all
            ! categories
            !--------------------------------------------------------

            !$OMP PARALLEL DO PRIVATE(iblk,i,j)
            do iblk = 1, nblocks
               do j = 1, ny_block
               do i = 1, nx_block
                  fcondtopn_f(i,j,n,iblk) = botmelt(i,j,iblk)
                  fsurfn_f(i,j,n,iblk)    = topmelt(i,j,iblk) &
                                            + botmelt(i,j,iblk)
                  flatn_f(i,j,n,iblk)    = - sublim(i,j,iblk)*Lsub
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO

         enddo  ! ncat

      endif   ! .not. calc_Tsfc

      end subroutine hadgem_data

!=======================================================================
! monthly forcing
!=======================================================================

      subroutine monthly_files (yr)

! Construct filenames based on the LANL naming conventions for NCAR data.
! Edit for other directory structures or filenames.
! Note: The year number in these filenames does not matter, because
!       subroutine file_year will insert the correct year.

! author: Elizabeth C. Hunke, LANL

      integer (kind=int_kind), intent(in) :: &
           yr                   ! current forcing year

      character(len=*), parameter :: subname = '(monthly_files)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      flw_file = &
           trim(atm_data_dir)//'/MONTHLY/cldf.omip.dat'

      rain_file = &
           trim(atm_data_dir)//'/MONTHLY/prec.nmyr.dat'

      tair_file = &
           trim(atm_data_dir)//'/MONTHLY/t_10.1996.dat'
      call file_year(tair_file,yr)

      humid_file = &
           trim(atm_data_dir)//'/MONTHLY/q_10.1996.dat'
      call file_year(humid_file,yr)

      ! stress/speed is used instead of wind components
      strax_file = &
           trim(atm_data_dir)//'/MONTHLY/strx.1996.dat'
      call file_year(strax_file,yr)

      stray_file = &
           trim(atm_data_dir)//'/MONTHLY/stry.1996.dat'
      call file_year(stray_file,yr)

      wind_file = &
           trim(atm_data_dir)//'/MONTHLY/wind.1996.dat'
      call file_year(wind_file,yr)

      if (my_task == master_task) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'Forcing data year = ', fyear
         write (nu_diag,*) 'Atmospheric data files:'
         write (nu_diag,*) trim(flw_file)
         write (nu_diag,*) trim(rain_file)
         write (nu_diag,*) trim(tair_file)
         write (nu_diag,*) trim(humid_file)
         write (nu_diag,*) trim(uwind_file)
         write (nu_diag,*) trim(vwind_file)
      endif                     ! master_task

      end subroutine monthly_files

!=======================================================================
! read monthly atmospheric data

      subroutine monthly_data

      use ice_blocks, only: block, get_block
      use ice_global_reductions, only: global_minval, global_maxval
      use ice_domain, only: nblocks, distrb_info, blocks_ice
      use ice_flux, only: fsnow, Tair, Qa, wind, strax, stray, fsw
      use ice_grid, only: hm, tlon, tlat, tmask, umask

      integer (kind=int_kind) :: &
          i, j        , &
          ixm,ixp     , & ! record numbers for neighboring months
          maxrec      , & ! maximum record number
          recslot     , & ! spline slot for current record
          midmonth    , & ! middle day of month
          iblk        , & ! block index
          ilo,ihi,jlo,jhi ! beginning and end of physical domain

      real (kind=dbl_kind) :: &
          vmin, vmax

      logical (kind=log_kind) :: readm

      type (block) :: &
         this_block           ! block information for current block

      character(len=*), parameter :: subname = '(monthly_data)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

    !-------------------------------------------------------------------
    ! monthly data
    !
    ! Assume that monthly data values are located in the middle of the
    ! month.
    !-------------------------------------------------------------------

      midmonth = 15  ! data is given on 15th of every month
!      midmonth = fix(p5 * real(daymo(mmonth)))  ! exact middle

      ! Compute record numbers for surrounding months
      maxrec = 12
      ixm  = mod(mmonth+maxrec-2,maxrec) + 1
      ixp  = mod(mmonth,         maxrec) + 1
      if (mday >= midmonth) ixm = -99  ! other two points will be used
      if (mday <  midmonth) ixp = -99

      ! Determine whether interpolation will use values 1:2 or 2:3
      ! recslot = 2 means we use values 1:2, with the current value (2)
      !  in the second slot
      ! recslot = 1 means we use values 2:3, with the current value (2)
      !  in the first slot
      recslot = 1                             ! latter half of month
      if (mday < midmonth) recslot = 2        ! first half of month

      ! Find interpolation coefficients
      call interp_coeff_monthly (recslot)

      ! Read 2 monthly values
      readm = .false.
      if (istep==1 .or. (mday==midmonth .and. msec==0)) readm = .true.

      call read_clim_data (readm, 0, ixm, mmonth, ixp,  &
             flw_file, cldf_data, &
             field_loc_center, field_type_scalar)
      call read_clim_data (readm, 0, ixm, mmonth, ixp,  &
             rain_file, fsnow_data, &
             field_loc_center, field_type_scalar)
      call read_clim_data (readm, 0, ixm, mmonth, ixp,  &
             tair_file, Tair_data, &
             field_loc_center, field_type_scalar)
      call read_clim_data (readm, 0, ixm, mmonth, ixp,  &
             humid_file, Qa_data, &
             field_loc_center, field_type_scalar)
      call read_clim_data (readm, 0, ixm, mmonth, ixp,  &
             wind_file, wind_data, &
             field_loc_center, field_type_scalar)
      call read_clim_data (readm, 0, ixm, mmonth, ixp,  &
             strax_file, strax_data, &
             field_loc_center, field_type_vector)
      call read_clim_data (readm, 0, ixm, mmonth, ixp,  &
             stray_file, stray_data, &
             field_loc_center, field_type_vector)

      call interpolate_data (cldf_data, cldf)
      call interpolate_data (fsnow_data, fsnow)  ! units mm/s = kg/m^2/s
      call interpolate_data (Tair_data, Tair)
      call interpolate_data (Qa_data, Qa)
      call interpolate_data (wind_data, wind)
      call interpolate_data (strax_data, strax)
      call interpolate_data (stray_data, stray)

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks
        call Qa_fixLY(nx_block,  ny_block, &
                                 Tair (:,:,iblk), &
                                 Qa   (:,:,iblk))

        do j = 1, ny_block
          do i = 1, nx_block
            Qa   (i,j,iblk) = Qa   (i,j,iblk) * hm(i,j,iblk)
            Tair (i,j,iblk) = Tair (i,j,iblk) * hm(i,j,iblk)
            wind (i,j,iblk) = wind (i,j,iblk) * hm(i,j,iblk)
            strax(i,j,iblk) = strax(i,j,iblk) * hm(i,j,iblk)
            stray(i,j,iblk) = stray(i,j,iblk) * hm(i,j,iblk)
          enddo
        enddo

      ! AOMIP
      this_block = get_block(blocks_ice(iblk),iblk)
      ilo = this_block%ilo
      ihi = this_block%ihi
      jlo = this_block%jlo
      jhi = this_block%jhi

      call compute_shortwave(nx_block, ny_block, &
                             ilo, ihi, jlo, jhi, &
                             TLON (:,:,iblk), &
                             TLAT (:,:,iblk), &
                             hm   (:,:,iblk), &
                             Qa   (:,:,iblk), &
                             cldf (:,:,iblk), &
                             fsw  (:,:,iblk))

      enddo  ! iblk
      !$OMP END PARALLEL DO

         if (debug_forcing) then
           if (my_task == master_task) write (nu_diag,*) 'LY_bulk_data'
           vmin = global_minval(fsw,distrb_info,tmask)
           vmax = global_maxval(fsw,distrb_info,tmask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'fsw',vmin,vmax
           vmin = global_minval(cldf,distrb_info,tmask)
           vmax = global_maxval(cldf,distrb_info,tmask)
           if (my_task.eq.master_task) &
               write (nu_diag,*) 'cldf',vmin,vmax
           vmin =global_minval(fsnow,distrb_info,tmask)
           vmax =global_maxval(fsnow,distrb_info,tmask)
           if (my_task.eq.master_task) &
               write (nu_diag,*) 'fsnow',vmin,vmax
           vmin = global_minval(Tair,distrb_info,tmask)
           vmax = global_maxval(Tair,distrb_info,tmask)
           if (my_task.eq.master_task) &
               write (nu_diag,*) 'Tair',vmin,vmax
           vmin = global_minval(wind,distrb_info,umask)
           vmax = global_maxval(wind,distrb_info,umask)
           if (my_task.eq.master_task) &
               write (nu_diag,*) 'wind',vmin,vmax
           vmin = global_minval(strax,distrb_info,umask)
           vmax = global_maxval(strax,distrb_info,umask)
           if (my_task.eq.master_task) &
               write (nu_diag,*) 'strax',vmin,vmax
           vmin = global_minval(stray,distrb_info,umask)
           vmax = global_maxval(stray,distrb_info,umask)
           if (my_task.eq.master_task) &
               write (nu_diag,*) 'stray',vmin,vmax
           vmin = global_minval(Qa,distrb_info,tmask)
           vmax = global_maxval(Qa,distrb_info,tmask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'Qa',vmin,vmax

        endif                   ! debug_forcing

      end subroutine monthly_data

!=======================================================================
! Oned atmospheric data
!=======================================================================

      subroutine oned_data

      use ice_flux, only: uatm, vatm, Tair, fsw, fsnow, Qa, rhoa, frain

      ! local parameters

      character (char_len_long) :: &
         met_file,   &    ! netcdf filename
         fieldname        ! field name in netcdf file

      integer (kind=int_kind) :: &
         fid              ! file id for netCDF file

      real (kind=dbl_kind):: &
         work             ! temporary variable

      logical (kind=log_kind) :: diag

      integer (kind=int_kind) :: &
         status           ! status flag

      real (kind=dbl_kind) :: & ! used to determine specific humidity
         Temp               , & ! air temperature (K)
         rh                 , & ! relative humidity (%)
         Psat               , & ! saturation vapour pressure (hPa)
         ws                     ! saturation mixing ratio

      real (kind=dbl_kind), parameter :: & ! coefficients for Hyland-Wexler Qa
         ps1 = 0.58002206e4_dbl_kind,    & ! (K)
         ps2 = 1.3914993_dbl_kind,       & !
         ps3 = 0.48640239e-1_dbl_kind,   & ! (K^-1)
         ps4 = 0.41764768e-4_dbl_kind,   & ! (K^-2)
         ps5 = 0.14452093e-7_dbl_kind,   & ! (K^-3)
         ps6 = 6.5459673_dbl_kind,       & !
         ws1 = 621.97_dbl_kind,          & ! for saturation mixing ratio
         Pair = 1020._dbl_kind             ! Sea level pressure (hPa)

      character(len=*), parameter :: subname = '(oned_data)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      diag = .false.   ! write diagnostic information

      if (trim(atm_data_format) == 'nc') then     ! read nc file

        ! hourly data beginning Jan 1, 1989, 01:00
        ! HARDWIRED for dt = 1 hour!
        met_file = uwind_file
        call ice_open_nc(met_file,fid)

        fieldname='Uatm'
        call ice_read_nc(fid,istep1,fieldname,work,diag)
        uatm(:,:,:) = work

        fieldname='Vatm'
        call ice_read_nc(fid,istep1,fieldname,work,diag)
        vatm(:,:,:) = work

        fieldname='Tair'
        call ice_read_nc(fid,istep1,fieldname,work,diag)
        Temp = work
        Tair(:,:,:) = Temp

        call ice_close_nc(fid)

        ! hourly solar data beginning Jan 1, 1989, 01:00
        met_file = fsw_file
        call ice_open_nc(met_file,fid)

        fieldname='fsw'
        call ice_read_nc(fid,istep1,fieldname,work,diag)
        fsw(:,:,:) = work

        call ice_close_nc(fid)

        ! hourly interpolated monthly  data beginning Jan 1, 1989, 01:00
        met_file = humid_file
        call ice_open_nc(met_file,fid)

        fieldname='rh'
        call ice_read_nc(fid,istep1,fieldname,work,diag)
        rh = work

        fieldname='fsnow'
        call ice_read_nc(fid,istep1,fieldname,work,diag)
        fsnow(:,:,:) = work

        call ice_close_nc(fid)

      !-------------------------------------------------------------------
      ! Find specific humidity using Hyland-Wexler formulation
      ! Hyland, R.W. and A. Wexler, Formulations for the Thermodynamic
      ! Properties of the saturated phases of H20 from 173.15K to 473.15K,
      ! ASHRAE Trans, 89(2A), 500-519, 1983
      !-------------------------------------------------------------------

        Psat = exp(-ps1/Temp + ps2 - ps3*Temp + ps4*Temp**2 - ps5 * Temp**3  &
              + ps6 * log(Temp))*p01          ! saturation vapour pressure
        ws = ws1 * Psat/(Pair - Psat)         ! saturation mixing ratio
        Qa(:,:,:) = rh * ws * p01/(c1 + rh * ws * p01) * p001
                                              ! specific humidity (kg/kg)
      endif ! atm_data_format

      ! flw calculated in prepare_forcing
        rhoa (:,:,:) = 1.3_dbl_kind ! air density (kg/m^3)
        cldf (:,:,:) = p25          ! cloud fraction
        frain(:,:,:) = c0           ! this is available in hourlymet_rh file

      end subroutine oned_data

!=======================================================================

      subroutine oned_files

      character(len=*), parameter :: subname = '(oned_files)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      fsw_file = &
           trim(atm_data_dir)//'/hourlysolar_brw1989_5yr.nc'

      rain_file = &
           trim(atm_data_dir)//'/hourlymet_rh_5yr.nc'

      uwind_file = &
           trim(atm_data_dir)//'/hourlymet_brw1989_5yr.nc'

      vwind_file = &
           trim(atm_data_dir)//'/hourlymet_brw1989_5yr.nc'

      tair_file = &
           trim(atm_data_dir)//'/hourlymet_brw1989_5yr.nc'

      humid_file = &
           trim(atm_data_dir)//'/hourlymet_rh_5yr.nc'

      if (my_task == master_task) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'Atmospheric data files:'
         write (nu_diag,*) trim(fsw_file)
         write (nu_diag,*) trim(rain_file)
         write (nu_diag,*) trim(uwind_file)
         write (nu_diag,*) trim(vwind_file)
         write (nu_diag,*) trim(tair_file)
         write (nu_diag,*) trim(humid_file)
      endif                     ! master_task

      end subroutine oned_files

!=======================================================================
! Climatological ocean forcing
!=======================================================================

      subroutine ocn_data_clim (dt)

! Interpolate monthly sss, sst data to timestep.
! Restore prognostic sst to data.
! Interpolate fields from U grid to T grid if necessary.

! author: Elizabeth C. Hunke and William H. Lipscomb, LANL

      use ice_domain, only: nblocks
      use ice_flux, only: sss, sst

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) :: &
          i, j, iblk  , & ! horizontal indices
          ixm,ixp     , & ! record numbers for neighboring months
          maxrec      , & ! maximum record number
          recslot     , & ! spline slot for current record
          midmonth        ! middle day of month

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
          sstdat              ! data value toward which SST is restored

      logical (kind=log_kind) :: readm

      character(len=*), parameter :: subname = '(ocn_data_clim)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      if (my_task == master_task .and. istep == 1) then
         if (trim(ocn_data_type)=='clim') then
            write (nu_diag,*) ' '
            write (nu_diag,*) 'SSS data interpolated to timestep:'
            write (nu_diag,*) trim(sss_file)
            write (nu_diag,*) ' '
            write (nu_diag,*) 'SST data interpolated to timestep:'
            write (nu_diag,*) trim(sst_file)
            if (restore_ocn) write (nu_diag,*) &
              'SST restoring timescale (days) =', trestore
         endif
      endif                     ! my_task, istep

    !-------------------------------------------------------------------
    ! monthly data
    !
    ! Assume that monthly data values are located in the middle of the
    ! month.
    !-------------------------------------------------------------------

      if (trim(ocn_data_type)=='clim') then

         midmonth = 15          ! data is given on 15th of every month
!!!      midmonth = fix(p5 * real(daymo(mmonth)))  ! exact middle

         ! Compute record numbers for surrounding months
         maxrec = 12
         ixm  = mod(mmonth+maxrec-2,maxrec) + 1
         ixp  = mod(mmonth,         maxrec) + 1
         if (mday >= midmonth) ixm = -99 ! other two points will be used
         if (mday <  midmonth) ixp = -99

         ! Determine whether interpolation will use values 1:2 or 2:3
         ! recslot = 2 means we use values 1:2, with the current value (2)
         !  in the second slot
         ! recslot = 1 means we use values 2:3, with the current value (2)
         !  in the first slot
         recslot = 1            ! latter half of month
         if (mday < midmonth) recslot = 2 ! first half of month

         ! Find interpolation coefficients
         call interp_coeff_monthly (recslot)

         readm = .false.
         if (istep==1 .or. (mday==midmonth .and. msec==0)) readm = .true.

    !-------------------------------------------------------------------
    ! Read two monthly SSS values and interpolate.
    ! Note: SSS is restored instantaneously to data.
    !-------------------------------------------------------------------

         call read_clim_data (readm, 0, ixm, mmonth, ixp, &
                              sss_file, sss_data, &
                              field_loc_center, field_type_scalar)
         call interpolate_data (sss_data, sss)

         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               sss(i,j,iblk) = max(sss(i,j,iblk), c0)
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO

         call ocn_freezing_temperature
      endif

    !-------------------------------------------------------------------
    ! Read two monthly SST values and interpolate.
    ! Restore toward interpolated value.
    !-------------------------------------------------------------------

      if (trim(ocn_data_type)=='clim') then
         call read_clim_data (readm, 0, ixm, mmonth, ixp, &
                              sst_file, sst_data, &
                              field_loc_center, field_type_scalar)
         call interpolate_data (sst_data, sstdat)

         if (restore_ocn) then
         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               sst(i,j,iblk) = sst(i,j,iblk)  &
                         + (sstdat(i,j,iblk)-sst(i,j,iblk))*dt/trest
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO
         endif
      endif

      end subroutine ocn_data_clim

!=======================================================================
! NCAR CESM M-configuration (AIO) ocean forcing
!=======================================================================

      subroutine ocn_data_ncar_init

! Reads NCAR pop ocean forcing data set 'pop_frc_gx1v3_010815.nc'
!
! List of ocean forcing fields: Note that order is important!
! (order is determined by field list in vname).
!
! For ocean mixed layer-----------------------------units
!
! 1  sst------temperature---------------------------(C)
! 2  sss------salinity------------------------------(ppt)
! 3  hbl------depth---------------------------------(m)
! 4  u--------surface u current---------------------(m/s)
! 5  v--------surface v current---------------------(m/s)
! 6  dhdx-----surface tilt x direction--------------(m/m)
! 7  dhdy-----surface tilt y direction--------------(m/m)
! 8  qdp------ocean sub-mixed layer heat flux-------(W/m2)
!
! Fields 4, 5, 6, 7 are on the U-grid; 1, 2, 3, and 8 are
! on the T-grid.

! authors: Bruce Briegleb, NCAR
!          Elizabeth Hunke, LANL

      use ice_blocks, only: nx_block, ny_block
      use ice_domain_size, only: max_blocks
#ifdef USE_NETCDF
      use netcdf
#endif

      integer (kind=int_kind) :: &
        n   , & ! field index
        m   , & ! month index
        nrec, & ! record number for direct access
        nbits

      character(char_len) :: &
        vname(nfld) ! variable names to search for in file
      data vname /  &
           'T',      'S',      'hblt',  'U',     'V', &
           'dhdx',   'dhdy',   'qdp' /

      integer (kind=int_kind) :: &
        status  , & ! status flag
        fid     , & ! file id
        dimid   , & ! dimension id
        nlat    , & ! number of longitudes of data
        nlon        ! number of latitudes  of data

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1

      character(len=*), parameter :: subname = '(ocn_data_ncar_init)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      if (my_task == master_task) then

         write (nu_diag,*) 'WARNING: evp_prep calculates surface tilt'
         write (nu_diag,*) 'WARNING: stress from geostrophic currents,'
         write (nu_diag,*) 'WARNING: not data from ocean forcing file.'
         write (nu_diag,*) 'WARNING: Alter ice_dyn_evp.F90 if desired.'

         if (restore_ocn) write (nu_diag,*)  &
             'SST restoring timescale = ',trestore,' days'

         sst_file = trim(ocn_data_dir)//'/'//trim(oceanmixed_file) ! not just sst

        !---------------------------------------------------------------
        ! Read in ocean forcing data from an existing file
        !---------------------------------------------------------------
        write (nu_diag,*) 'ocean mixed layer forcing data file = ', &
                           trim(sst_file)

      endif ! master_task

      if (trim(ocn_data_format) == 'nc') then
#ifdef USE_NETCDF
        if (my_task == master_task) then
          call ice_open_nc(sst_file, fid)

!          status = nf90_inq_dimid(fid,'nlon',dimid)
          status = nf90_inq_dimid(fid,'ni',dimid)
          call ice_check_nc(status, subname//' ERROR: inq dimid ni', file=__FILE__, line=__LINE__)
          status = nf90_inquire_dimension(fid,dimid,len=nlon)
          call ice_check_nc(status, subname//' ERROR: inq dim ni', file=__FILE__, line=__LINE__)

!          status = nf90_inq_dimid(fid,'nlat',dimid)
          status = nf90_inq_dimid(fid,'nj',dimid)
          call ice_check_nc(status, subname//' ERROR: inq dimid nj', file=__FILE__, line=__LINE__)
          status = nf90_inquire_dimension(fid,dimid,len=nlat)
          call ice_check_nc(status, subname//' ERROR: inq dim nj', file=__FILE__, line=__LINE__)

          if( nlon .ne. nx_global ) then
            call abort_ice (error_message=subname//'ice: ocn frc file nlon ne nx_global', &
               file=__FILE__, line=__LINE__)
          endif
          if( nlat .ne. ny_global ) then
            call abort_ice (error_message=subname//'ice: ocn frc file nlat ne ny_global', &
               file=__FILE__, line=__LINE__)
          endif

        endif ! master_task

        ! Read in ocean forcing data for all 12 months
        do n=1,nfld
          do m=1,12

            ! Note: netCDF does single to double conversion if necessary
!           if (n >= 4 .and. n <= 7) then
!              call ice_read_nc(fid, m, vname(n), work1, debug_forcing, &
!                               field_loc_NEcorner, field_type_vector)
!           else
               call ice_read_nc(fid, m, vname(n), work1, debug_forcing, &
                                field_loc_center, field_type_scalar)
!           endif

            ocn_frc_m(:,:,:,n,m) = work1(:,:,:)

          enddo               ! month loop
        enddo               ! field loop

        if (my_task == master_task) call ice_close_nc(fid)
#else
      call abort_ice(subname//'ERROR: USE_NETCDF cpp not defined for '//trim(sst_file), &
          file=__FILE__, line=__LINE__)
#endif

      else  ! binary format

        nbits = 64
        call ice_open (nu_forcing, sst_file, nbits)

        nrec = 0
        do n=1,nfld
           do m=1,12
              nrec = nrec + 1
              if (n >= 4 .and. n <= 7) then
                call ice_read (nu_forcing, nrec, work1, 'rda8', debug_forcing, &
                               field_loc_NEcorner, field_type_vector)
              else
                call ice_read (nu_forcing, nrec, work1, 'rda8', debug_forcing, &
                               field_loc_center, field_type_scalar)
              endif
              ocn_frc_m(:,:,:,n,m) = work1(:,:,:)
           enddo               ! month loop
        enddo               ! field loop
        close (nu_forcing)

      endif

!echmod - currents cause Fram outflow to be too large
!             ocn_frc_m(:,:,:,4,:) = c0
!             ocn_frc_m(:,:,:,5,:) = c0
!echmod

      end subroutine ocn_data_ncar_init

!=======================================================================

      subroutine ocn_data_ncar_init_3D

! Reads NCAR pop ocean forcing data set 'oceanmixed_ice_depth.nc'
!
! List of ocean forcing fields: Note that order is important!
! (order is determined by field list in vname).
!
! For ocean mixed layer-----------------------------units
!
! 1  sst------temperature---------------------------(C)
! 2  sss------salinity------------------------------(ppt)
! 3  hbl------depth---------------------------------(m)
! 4  u--------surface u current---------------------(m/s)
! 5  v--------surface v current---------------------(m/s)
! 6  dhdx-----surface tilt x direction--------------(m/m)
! 7  dhdy-----surface tilt y direction--------------(m/m)
! 8  qdp------ocean sub-mixed layer heat flux-------(W/m2)
!
! All fields are on the T-grid.
!
! authors: Bruce Briegleb, NCAR
!          Elizabeth Hunke, LANL

      use ice_blocks, only: nx_block, ny_block
      use ice_domain_size, only: max_blocks
      use ice_grid, only: grid_average_X2Y, ANGLET
      use ice_read_write, only: ice_read_nc_uv
#ifdef USE_NETCDF
      use netcdf
#endif

#ifdef USE_NETCDF
      integer (kind=int_kind) :: &
        n   , & ! field index
        m   , & ! month index
        nzlev   ! z level of currents

      character(char_len) :: &
        vname(nfld) ! variable names to search for in file
      data vname /  &
           'T',      'S',      'hblt',  'U',     'V', &
           'dhdx',   'dhdy',   'qdp' /

      integer (kind=int_kind) :: &
        fid        , & ! file id
        dimid          ! dimension id

      integer (kind=int_kind) :: &
        status  , & ! status flag
        nlat    , & ! number of longitudes of data
        nlon        ! number of latitudes  of data

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1, work2
#endif

      character(len=*), parameter :: subname = '(ocn_data_ncar_init_3D)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      if (my_task == master_task) then

         write (nu_diag,*) 'WARNING: evp_prep calculates surface tilt'
         write (nu_diag,*) 'WARNING: stress from geostrophic currents,'
         write (nu_diag,*) 'WARNING: not data from ocean forcing file.'
         write (nu_diag,*) 'WARNING: Alter ice_dyn_evp.F if desired.'

         if (restore_ocn) write (nu_diag,*)  &
             'SST restoring timescale = ',trestore,' days'

         sst_file = trim(ocn_data_dir)//'/'//trim(oceanmixed_file) ! not just sst

        !---------------------------------------------------------------
        ! Read in ocean forcing data from an existing file
        !---------------------------------------------------------------
        write (nu_diag,*) 'ocean mixed layer forcing data file = ', &
                           trim(sst_file)
        write (nu_diag,*)

      endif ! master_task

      if (trim(ocn_data_format) == 'nc') then
#ifdef USE_NETCDF
        if (my_task == master_task) then
          call ice_open_nc(sst_file, fid)

!          status = nf90_inq_dimid(fid,'nlon',dimid)
          status = nf90_inq_dimid(fid,'ni',dimid)
          call ice_check_nc(status, subname//' ERROR: inq dimid ni', file=__FILE__, line=__LINE__)
          status = nf90_inquire_dimension(fid,dimid,len=nlon)
          call ice_check_nc(status, subname//' ERROR: inq dim ni', file=__FILE__, line=__LINE__)

!          status = nf90_inq_dimid(fid,'nlat',dimid)
          status = nf90_inq_dimid(fid,'nj',dimid)
          call ice_check_nc(status, subname//' ERROR: inq dimid nj', file=__FILE__, line=__LINE__)
          status = nf90_inquire_dimension(fid,dimid,len=nlat)
          call ice_check_nc(status, subname//' ERROR: inq dim nj', file=__FILE__, line=__LINE__)

          if( nlon .ne. nx_global ) then
            call abort_ice (error_message=subname//'ice: ocn frc file nlon ne nx_global', &
               file=__FILE__, line=__LINE__)
          endif
          if( nlat .ne. ny_global ) then
            call abort_ice (error_message=subname//'ice: ocn frc file nlat ne ny_global', &
               file=__FILE__, line=__LINE__)
          endif

        endif ! master_task

        ! Read in ocean forcing data for all 12 months
        do n=1,nfld
          do m=1,12

            ! Note: netCDF does single to double conversion if necessary
            if (n == 4 .or. n == 5) then ! 3D currents
               nzlev = 1                 ! surface currents
               call ice_read_nc_uv(fid, m, nzlev, vname(n), work1, debug_forcing, &
                                field_loc_center, field_type_scalar)
            else
               call ice_read_nc(fid, m, vname(n), work1, debug_forcing, &
                                field_loc_center, field_type_scalar)
            endif

            ! the land mask used in ocean_mixed_depth.nc does not
            ! match our gx1v3 mask (hm)
            where (work1(:,:,:) < -900.) work1(:,:,:) = c0

            ocn_frc_m(:,:,:,n,m) = work1(:,:,:)

          enddo               ! month loop
        enddo               ! field loop

        if (my_task == master_task) call ice_close_nc(fid)

        ! Rotate vector quantities and shift to U-grid
        do n=4,6,2
          do m=1,12

             work1(:,:,:) = ocn_frc_m(:,:,:,n  ,m)
             work2(:,:,:) = ocn_frc_m(:,:,:,n+1,m)
             ocn_frc_m(:,:,:,n  ,m) = work1(:,:,:)*cos(ANGLET(:,:,:)) &
                                    + work2(:,:,:)*sin(ANGLET(:,:,:))
             ocn_frc_m(:,:,:,n+1,m) = work2(:,:,:)*cos(ANGLET(:,:,:)) &
                                    - work1(:,:,:)*sin(ANGLET(:,:,:))

             work1(:,:,:) = ocn_frc_m(:,:,:,n  ,m)
             work2(:,:,:) = ocn_frc_m(:,:,:,n+1,m)
             call grid_average_X2Y('A',work1,'T',ocn_frc_m(:,:,:,n  ,m),'U')
             call grid_average_X2Y('A',work2,'T',ocn_frc_m(:,:,:,n+1,m),'U')

          enddo               ! month loop
        enddo               ! field loop

#else
        call abort_ice(subname//'ERROR: USE_NETCDF cpp not defined', &
            file=__FILE__, line=__LINE__)
#endif

      else  ! binary format

        call abort_ice (error_message=subname//'new ocean forcing is netcdf only', &
           file=__FILE__, line=__LINE__)

      endif

      end subroutine ocn_data_ncar_init_3D

!=======================================================================

      subroutine ocn_data_ncar(dt)

! Interpolate monthly ocean data to timestep.
! Restore sst if desired. sst is updated with surface fluxes in ice_ocean.F.

      use ice_blocks, only: nx_block, ny_block
      use ice_global_reductions, only: global_minval, global_maxval
      use ice_domain, only: nblocks, distrb_info
      use ice_domain_size, only: max_blocks
      use ice_flux, only: sss, sst, Tf, uocn, vocn, ss_tltx, ss_tlty, &
            qdp, hmix
      use ice_restart_shared, only: restart
      use ice_grid, only: hm, tmask, umask

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      integer (kind=int_kind) :: &
          i, j, n, iblk   , &
          ixm,ixp         , & ! record numbers for neighboring months
          maxrec          , & ! maximum record number
          recslot         , & ! spline slot for current record
          midmonth            ! middle day of month

      real (kind=dbl_kind) :: &
          vmin, vmax

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1

      character(len=*), parameter :: subname = '(ocn_data_ncar)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

    !-------------------------------------------------------------------
    ! monthly data
    !
    ! Assume that monthly data values are located in the middle of the
    ! month.
    !-------------------------------------------------------------------

      midmonth = 15  ! data is given on 15th of every month
!      midmonth = fix(p5 * real(daymo(mmonth),kind=dbl_kind))  ! exact middle

      ! Compute record numbers for surrounding months
      maxrec = 12
      ixm  = mod(mmonth+maxrec-2,maxrec) + 1
      ixp  = mod(mmonth,         maxrec) + 1
      if (mday >= midmonth) ixm = -99  ! other two points will be used
      if (mday <  midmonth) ixp = -99

      ! Determine whether interpolation will use values 1:2 or 2:3
      ! recslot = 2 means we use values 1:2, with the current value (2)
      !  in the second slot
      ! recslot = 1 means we use values 2:3, with the current value (2)
      !  in the first slot
      recslot = 1                             ! latter half of month
      if (mday < midmonth) recslot = 2        ! first half of month

      ! Find interpolation coefficients
      call interp_coeff_monthly (recslot)

      sst_data(:,:,:,:) = c0
      do n = nfld, 1, -1
        do iblk = 1, nblocks
        ! use sst_data arrays as temporary work space until n=1
        if (ixm /= -99) then  ! first half of month
          sst_data(:,:,1,iblk) = ocn_frc_m(:,:,iblk,n,ixm)
          sst_data(:,:,2,iblk) = ocn_frc_m(:,:,iblk,n,mmonth)
        else                 ! second half of month
          sst_data(:,:,1,iblk) = ocn_frc_m(:,:,iblk,n,mmonth)
          sst_data(:,:,2,iblk) = ocn_frc_m(:,:,iblk,n,ixp)
        endif
        enddo

        call interpolate_data (sst_data,work1)
        ! masking by hm is necessary due to NaNs in the data file
        do j = 1, ny_block
          do i = 1, nx_block
            if (n == 2) sss    (i,j,:) = c0
            if (n == 3) hmix   (i,j,:) = c0
            if (n == 4) uocn   (i,j,:) = c0
            if (n == 5) vocn   (i,j,:) = c0
            if (n == 6) ss_tltx(i,j,:) = c0
            if (n == 7) ss_tlty(i,j,:) = c0
            if (n == 8) qdp    (i,j,:) = c0
            do iblk = 1, nblocks
              if (hm(i,j,iblk) == c1) then
                if (n == 2) sss    (i,j,iblk) = work1(i,j,iblk)
                if (n == 3) hmix   (i,j,iblk) = max(mixed_layer_depth_default,work1(i,j,iblk))
                if (n == 4) uocn   (i,j,iblk) = work1(i,j,iblk)
                if (n == 5) vocn   (i,j,iblk) = work1(i,j,iblk)
                if (n == 6) ss_tltx(i,j,iblk) = work1(i,j,iblk)
                if (n == 7) ss_tlty(i,j,iblk) = work1(i,j,iblk)
                if (n == 8) qdp    (i,j,iblk) = work1(i,j,iblk)
              endif
            enddo
          enddo
        enddo
      enddo

      do j = 1, ny_block
         do i = 1, nx_block
            sss (i,j,:) = max (sss(i,j,:), c0)
            hmix(i,j,:) = max(hmix(i,j,:), c0)
         enddo
      enddo

      call ocn_freezing_temperature

      if (restore_ocn) then
        do j = 1, ny_block
         do i = 1, nx_block
           sst(i,j,:) = sst(i,j,:) + (work1(i,j,:)-sst(i,j,:))*dt/trest
         enddo
        enddo
!     else sst is only updated in ice_ocean.F
      endif

      ! initialize sst properly on first step
      if (istep1 <= 1 .and. .not. (restart)) then
        call interpolate_data (sst_data,sst)
        !$OMP PARALLEL DO PRIVATE(iblk,i,j)
        do iblk = 1, nblocks
         do j = 1, ny_block
          do i = 1, nx_block
            if (hm(i,j,iblk) == c1) then
              sst(i,j,iblk) =  max (sst(i,j,iblk), Tf(i,j,iblk))
            else
              sst(i,j,iblk) = c0
            endif
          enddo
         enddo
        enddo
        !$OMP END PARALLEL DO
      endif

      if (debug_forcing) then
         if (my_task == master_task)  &
               write (nu_diag,*) 'ocn_data_ncar'
           vmin = global_minval(Tf,distrb_info,tmask)
           vmax = global_maxval(Tf,distrb_info,tmask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'Tf',vmin,vmax
           vmin = global_minval(sst,distrb_info,tmask)
           vmax = global_maxval(sst,distrb_info,tmask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'sst',vmin,vmax
           vmin = global_minval(sss,distrb_info,tmask)
           vmax = global_maxval(sss,distrb_info,tmask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'sss',vmin,vmax
           vmin = global_minval(hmix,distrb_info,tmask)
           vmax = global_maxval(hmix,distrb_info,tmask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'hmix',vmin,vmax
           vmin = global_minval(uocn,distrb_info,umask)
           vmax = global_maxval(uocn,distrb_info,umask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'uocn',vmin,vmax
           vmin = global_minval(vocn,distrb_info,umask)
           vmax = global_maxval(vocn,distrb_info,umask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'vocn',vmin,vmax
           vmin = global_minval(ss_tltx,distrb_info,umask)
           vmax = global_maxval(ss_tltx,distrb_info,umask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'ss_tltx',vmin,vmax
           vmin = global_minval(ss_tlty,distrb_info,umask)
           vmax = global_maxval(ss_tlty,distrb_info,umask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'ss_tlty',vmin,vmax
           vmin = global_minval(qdp,distrb_info,tmask)
           vmax = global_maxval(qdp,distrb_info,tmask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'qdp',vmin,vmax
      endif

      end subroutine ocn_data_ncar

!=======================================================================
! ocean data for oned configuration
! Current (released) values are the same as the defaults (ice_flux.F90)

      subroutine ocn_data_oned

      use ice_flux, only: sss, sst, Tf, uocn, vocn, ss_tltx, ss_tlty, &
            qdp, hmix, frzmlt

      character(len=*), parameter :: subname = '(ocn_data_oned)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      sss    (:,:,:) = 34.0_dbl_kind   ! sea surface salinity (ppt)

      call ocn_freezing_temperature

      sst    (:,:,:) = Tf(:,:,:)       ! sea surface temp (C)
      uocn   (:,:,:) = c0              ! surface ocean currents (m/s)
      vocn   (:,:,:) = c0
      ss_tltx(:,:,:) = c0              ! sea surface tilt (m/m)
      ss_tlty(:,:,:) = c0
      frzmlt (:,:,:) = c0              ! freezing/melting potential (W/m^2)
      qdp    (:,:,:) = c0              ! deep ocean heat flux (W/m^2)
      hmix   (:,:,:) = mixed_layer_depth_default   ! ocean mixed layer depth

      end subroutine ocn_data_oned

!=======================================================================

      subroutine ocn_data_hadgem(dt)

!  Reads in HadGEM ocean forcing data as required from netCDF files
!  Current options (selected by ocn_data_type)
!  hadgem_sst: read in sst only
!  hadgem_sst_uvocn: read in sst plus uocn and vocn

! authors: Ann Keen, Met Office

      use ice_domain, only: nblocks
      use ice_domain_size, only: max_blocks
      use ice_flux, only: sst, uocn, vocn
      use ice_grid, only: ANGLET

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      integer (kind=int_kind) :: &
          i, j, iblk  , &
          ixm,ixp     , & ! record numbers for neighboring months
          maxrec      , & ! maximum record number
          recslot     , & ! spline slot for current record
          midmonth        ! middle day of month

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
          sstdat              ! data value toward which SST is restored

      real (kind=dbl_kind) :: workx, worky

      logical (kind=log_kind) :: readm

      character (char_len) :: &
            fieldname     ! field name in netcdf file

      character (char_len_long) :: &
            filename      ! name of netCDF file

      character(len=*), parameter :: subname = '(ocn_data_hadgem)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

    !-------------------------------------------------------------------
    ! monthly data
    !
    ! Assume that monthly data values are located in the middle of the
    ! month.
    !-------------------------------------------------------------------

      midmonth = 15  ! data is given on 15th of every month
!      midmonth = fix(p5 * real(daymo(mmonth)))  ! exact middle

      ! Compute record numbers for surrounding months
      maxrec = 12
      ixm  = mod(mmonth+maxrec-2,maxrec) + 1
      ixp  = mod(mmonth,         maxrec) + 1
      if (mday >= midmonth) ixm = -99  ! other two points will be used
      if (mday <  midmonth) ixp = -99

      ! Determine whether interpolation will use values 1:2 or 2:3
      ! recslot = 2 means we use values 1:2, with the current value (2)
      !  in the second slot
      ! recslot = 1 means we use values 2:3, with the current value (2)
      !  in the first slot
      recslot = 1                             ! latter half of month
      if (mday < midmonth) recslot = 2        ! first half of month

      ! Find interpolation coefficients
      call interp_coeff_monthly (recslot)

      ! Read 2 monthly values
      readm = .false.
      if (istep==1 .or. (mday==midmonth .and. msec==0)) readm = .true.

      if (my_task == master_task .and. istep == 1) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'SST data interpolated to timestep:'
         write (nu_diag,*) trim(ocn_data_dir)//'/MONTHLY/sst.1997.nc'
         if (restore_ocn) write (nu_diag,*) &
              'SST restoring timescale (days) =', trestore
         if (trim(ocn_data_type)=='hadgem_sst_uvocn') then
            write (nu_diag,*) ' '
            write (nu_diag,*) 'uocn and vocn interpolated to timestep:'
            write (nu_diag,*) trim(ocn_data_dir)//'/MONTHLY/uocn.1997.nc'
            write (nu_diag,*) trim(ocn_data_dir)//'/MONTHLY/vocn.1997.nc'
         endif
      endif                     ! my_task, istep

      ! -----------------------------------------------------------
      ! SST
      ! -----------------------------------------------------------
      sst_file = trim(ocn_data_dir)//'/MONTHLY/sst.1997.nc'
      fieldname='sst'
      call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, &
                      maxrec, sst_file, fieldname, sst_data, &
                      field_loc_center, field_type_scalar)

      ! Interpolate to current time step
      call interpolate_data (sst_data, sstdat)

      ! Restore SSTs if required
        if (restore_ocn) then
         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               sst(i,j,iblk) = sst(i,j,iblk)  &
                         + (sstdat(i,j,iblk)-sst(i,j,iblk))*dt/trest
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO
         endif

      ! -----------------------------------------------------------
      ! Ocean currents
      ! --------------
      ! Values read in are on T grid and oriented geographically, hence
      ! vectors need to be rotated to model grid and then interpolated
      ! to U grid.
      ! Also need to be converted from cm s-1 (UM) to m s-1 (CICE)
      ! -----------------------------------------------------------

      if (trim(ocn_data_type)=='hadgem_sst_uvocn') then

        filename = trim(ocn_data_dir)//'/MONTHLY/uocn.1997.nc'
        fieldname='uocn'
        call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, &
                      maxrec, filename, fieldname, uocn_data, &
                      field_loc_center, field_type_vector)

        ! Interpolate to current time step
        call interpolate_data (uocn_data, uocn)

        filename = trim(ocn_data_dir)//'/MONTHLY/vocn.1997.nc'
        fieldname='vocn'
        call read_data_nc (readm, 0, fyear, ixm, mmonth, ixp, &
                      maxrec, filename, fieldname, vocn_data, &
                      field_loc_center, field_type_vector)

        ! Interpolate to current time step
        call interpolate_data (vocn_data, vocn)

     !-----------------------------------------------------------------
     ! Rotate zonal/meridional vectors to local coordinates,
     ! and change  units
     !-----------------------------------------------------------------

         !$OMP PARALLEL DO PRIVATE(iblk,i,j,workx,worky)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block

               workx      = uocn(i,j,iblk)
               worky      = vocn(i,j,iblk)
               uocn(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) &
                                  + worky*sin(ANGLET(i,j,iblk))
               vocn(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) &
                                  - workx*sin(ANGLET(i,j,iblk))

               uocn(i,j,iblk) = uocn(i,j,iblk) * cm_to_m
               vocn(i,j,iblk) = vocn(i,j,iblk) * cm_to_m

            enddo   ! i
            enddo   ! j
         enddo      ! nblocks
         !$OMP END PARALLEL DO

     !-----------------------------------------------------------------
     ! Interpolate to U grid
     !-----------------------------------------------------------------

         ! tcraig, this is now computed in dynamics for consistency

     endif    !   ocn_data_type = hadgem_sst_uvocn

     end subroutine ocn_data_hadgem

!=======================================================================

      subroutine ocn_data_hycom_init
        ! Read SSS+SST from a HYCOM file converted to NetCDF format.
        ! HYCOM binary2NetCDF: hcdata2ncdf2d (or hcdata2ncdf3z)
        !   + rename/link file
        use ice_blocks, only: nx_block, ny_block
        use ice_domain, only: nblocks
        use ice_flux, only: sss, sst, Tf

        integer (kind=int_kind) :: &
           i, j, iblk       , & ! horizontal indices
           fid                  ! file id for netCDF file

        character (char_len) :: &
           fieldname            ! field name in netcdf file

        character(len=*), parameter :: subname = '(ocn_data_hycom_init)'

        if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

        if (trim(ocn_data_type) == 'hycom') then
           sss_file = trim(ocn_data_dir)//'ice.restart.surf.nc'

           if (my_task == master_task) then
             write (nu_diag,*)' '
             write (nu_diag,*)'Initial SSS file: ',trim(sss_file)
           endif

           fieldname = 'sss'
           call ice_open_nc (sss_file, fid)
           call ice_read_nc (fid, 1 , fieldname, sss, debug_forcing, &
                             field_loc_center, field_type_scalar)
           call ice_close_nc(fid)

           call ocn_freezing_temperature

           sst_file = trim(ocn_data_dir)//'ice.restart.surf.nc'

           if (my_task == master_task) then
             write (nu_diag,*)' '
             write (nu_diag,*)'Initial SST file: ',trim(sst_file)
           endif

           fieldname = 'sst'
           call ice_open_nc (sst_file, fid)
           call ice_read_nc (fid, 1 , fieldname, sst, debug_forcing, &
                                field_loc_center, field_type_scalar)
           call ice_close_nc(fid)

           ! Make sure sst is not less than freezing temperature Tf
           !$OMP PARALLEL DO PRIVATE(iblk,i,j)
           do iblk = 1, nblocks
              do j = 1, ny_block
              do i = 1, nx_block
                 sst(i,j,iblk) = max(sst(i,j,iblk),Tf(i,j,iblk))
              enddo
              enddo
           enddo
           !$OMP END PARALLEL DO
        endif

      end subroutine ocn_data_hycom_init

!=======================================================================

      subroutine hycom_atm_files

      use ice_broadcast, only: broadcast_array, broadcast_scalar

      integer (kind = int_kind) :: &
            fid          ! File id
      character (char_len) :: &
            varname      ! variable name in netcdf file
      character(len=*), parameter :: subname = '(hycom_atm_files)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      fsw_file   = trim(atm_data_dir)//'/forcing.shwflx.nc'
      flw_file   = trim(atm_data_dir)//'/forcing.radflx.nc'
      rain_file  = trim(atm_data_dir)//'/forcing.precip.nc'
      uwind_file = trim(atm_data_dir)//'/forcing.wndewd.nc'
      vwind_file = trim(atm_data_dir)//'/forcing.wndnwd.nc'
      tair_file  = trim(atm_data_dir)//'/forcing.airtmp.nc'
      humid_file = trim(atm_data_dir)//'/forcing.vapmix.nc'

      ! Read time vector from "tair_file"
      call ice_open_nc(tair_file, fid)
      varname='MT'
      call ice_get_ncvarsize(fid,varname,Njday_atm)
      call broadcast_scalar(Njday_atm,master_task)
      allocate(jday_atm(Njday_atm))
      call ice_read_vec_nc(fid,Njday_atm, varname,jday_atm, .true.)
      call ice_close_nc(fid)
      call broadcast_array(jday_atm ,master_task)

      ! Write diag info
      if (my_task == master_task) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'CICE: Atm. (hycomdate) Start = ',jday_atm(1)
         write (nu_diag,*) 'CICE: Atm. (hycomdate) End   = ',jday_atm(Njday_atm)
         write (nu_diag,*) 'CICE: Total Atm timesteps    = ',Njday_atm
         write (nu_diag,*) 'CICE: Atmospheric forcing files:'
         write (nu_diag,*) trim(fsw_file)
         write (nu_diag,*) trim(flw_file)
         write (nu_diag,*) trim(rain_file)
         write (nu_diag,*) trim(uwind_file)
         write (nu_diag,*) trim(vwind_file)
         write (nu_diag,*) trim(tair_file)
         write (nu_diag,*) trim(humid_file)
      endif                     ! master_task

      end subroutine hycom_atm_files

!=======================================================================

      subroutine hycom_atm_data

      use ice_flux, only: fsw, fsnow, Tair, uatm, vatm, Qa, flw
      use ice_domain, only: nblocks

      integer (kind=int_kind) :: &
          recnum       ! record number

      real (kind=dbl_kind) :: &
          hcdate         ! current time in HYCOM jday units

      logical (kind=log_kind) :: read6

      character (char_len) :: &
            fieldname    ! field name in netcdf file

      integer (kind=int_kind) :: &
         i, j, iblk      ! horizontal indices

      real (kind=dbl_kind) :: Tffresh, secday

      character(len=*), parameter :: subname = '(hycom_atm_data)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      call icepack_query_parameters(Tffresh_out=Tffresh)
      call icepack_query_parameters(secday_out=secday)

      ! current time in HYCOM jday units (HYCOM ref year: 1900,12,31,000000)
      hcdate = real(compute_days_between(1900,12,31,myear,mmonth,mday)) + msec/secday

      ! Init recnum try
      recnum=min(max(oldrecnum,1),Njday_atm-1)

      ! Find correct time in ATM data ... assume cont. incr. time-axis
      do while ( recnum<Njday_atm )
        if ( hcdate>=jday_atm(recnum) .and. &
             hcdate<=jday_atm(recnum+1) ) exit
        if ( abs(hcdate-jday_atm(recnum))<p001 ) exit  ! Accept within tolerance = 0.001 days
        if ( abs(hcdate-jday_atm(recnum+1))<p001 ) exit  ! Accept within tolerance = 0.001 days
        recnum=recnum+1
      enddo

      ! Last Atm date might be the same as last CICE date.
      recnum=min(recnum,Njday_atm-1)

      ! Check if current time do not exceed last forcing time
      ! + check forcing is available before (or at) current forcing time
      if ( hcdate>jday_atm(recnum+1)+p001 .or. hcdate<jday_atm(recnum)-p001) then
         write (nu_diag,*) &
         'ERROR: CICE: Atm forcing not available at hcdate =',hcdate
         write (nu_diag,*) &
         'ERROR: CICE: myear, yday ,msec = ',myear, yday, msec
         call abort_ice ('ERROR: CICE stopped')
      endif

      ! Compute interpolation coefficients
      call interp_coeff2 (hcdate, jday_atm(recnum), jday_atm(recnum+1) )


      ! Read
      read6 = .false.
      if (istep==1 .or. oldrecnum /= recnum) read6 = .true.

      if (trim(atm_data_format) == 'nc') then

         if (read6 .and. my_task == master_task) write(nu_diag,*) &
           'CICE: Atm. read: = ',jday_atm(recnum), jday_atm(recnum+1)

         fieldname = 'airtmp'
         call read_data_nc_hycom (read6, recnum, &
                          tair_file, fieldname, Tair_data, &
                          field_loc_center, field_type_scalar)
         fieldname = 'wndewd'
         call read_data_nc_hycom (read6, recnum, &
                          uwind_file, fieldname, uatm_data, &
                          field_loc_center, field_type_vector)
         fieldname = 'wndnwd'
         call read_data_nc_hycom (read6, recnum, &
                          vwind_file, fieldname, vatm_data, &
                          field_loc_center, field_type_vector)
         fieldname = 'vapmix'
         call read_data_nc_hycom (read6, recnum, &
                          humid_file, fieldname, Qa_data,  &
                          field_loc_center, field_type_scalar)
         fieldname = 'shwflx'
         call read_data_nc_hycom (read6, recnum, &
                          fsw_file, fieldname, fsw_data,   &
                          field_loc_center, field_type_scalar)
         fieldname = 'radflx'
         call read_data_nc_hycom (read6, recnum, &
                          flw_file, fieldname, flw_data,   &
                          field_loc_center, field_type_scalar)
         fieldname = 'precip'
         call read_data_nc_hycom (read6, recnum, &
                          rain_file, fieldname, fsnow_data,&
                          field_loc_center, field_type_scalar)

      else
         call abort_ice(subname//'ERROR: atm_data_format unavailable for hycom')
      endif

      ! Interpolate
      if (debug_forcing) then
        if (my_task == master_task) then
           write(nu_diag,*)'CICE: Atm. interpolate: = ',&
                                  hcdate,c1intp,c2intp
        endif
      endif
      call interpolate_data (Tair_data, Tair)
      call interpolate_data (uatm_data, uatm)
      call interpolate_data (vatm_data, vatm)
      call interpolate_data ( fsw_data, fsw)
      call interpolate_data ( flw_data, flw)
      call interpolate_data (  Qa_data, Qa)
      call interpolate_data (fsnow_data, fsnow)

      ! Adjust data forcing to CICE units
      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            ! Air temperature: Degrees --> Kelvin
            Tair(i,j,iblk) = Tair(i,j,iblk) + Tffresh
         enddo  ! i
         enddo  ! j
      enddo     ! nblocks
      !$OMP END PARALLEL DO

      ! Save record number for next time step
      oldrecnum = recnum

      end subroutine hycom_atm_data

!=======================================================================
!
      subroutine read_data_nc_point (flag, recd, yr, ixm, ixx, ixp, &
                            maxrec, data_file, fieldname, field_data, &
                            field_loc, field_type)
!
! If data is at the beginning of a one-year record, get data from
!  the previous year.
! If data is at the end of a one-year record, get data from the
!  following year.
! If no earlier data exists (beginning of fyear_init), then
!  (1) For monthly data, get data from the end of fyear_final.
!  (2) For more frequent data, let the ixm value equal the
!      first value of the year.
! If no later data exists (end of fyear_final), then
!  (1) For monthly data, get data from the beginning of fyear_init.
!  (2) For more frequent data, let the ixp value
!      equal the last value of the year.
! In other words, we assume persistence when daily or 6-hourly
!   data is missing, and we assume periodicity when monthly data
!   is missing.
!
      use ice_diagnostics, only: debug_model_step

      logical (kind=log_kind), intent(in) :: flag

      integer (kind=int_kind), intent(in) :: &
         recd                , & ! baseline record number
         yr                  , & ! year of forcing data
         ixm, ixx, ixp       , & ! record numbers of 3 data values
                                 ! relative to recd
         maxrec                  ! maximum record value

      character (char_len_long), intent(in) :: &
         data_file               ! data file to be read

      character (char_len), intent(in) :: &
         fieldname               ! field name in netCDF file

      integer (kind=int_kind), intent(in) :: &
           field_loc, &      ! location of field on staggered grid
           field_type        ! type of field (scalar, vector, angle)

      real (kind=dbl_kind), dimension(2), intent(inout) :: &
         field_data              ! 2 values needed for interpolation

      integer (kind=int_kind) :: &
         nrec             , & ! record number to read
         n2, n4           , & ! like ixm and ixp, but
                              ! adjusted at beginning and end of data
         arg              , & ! value of time argument in field_data
         fid                  ! file id for netCDF routines

      character(len=*), parameter :: subname = '(read_data_nc_point)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      call ice_timer_start(timer_readwrite)  ! reading/writing

      field_data = c0 ! to satisfy intent(out) attribute

      if (istep1 > debug_model_step) debug_forcing = .true.  !! debugging

      if (my_task==master_task .and. (debug_forcing)) then
         write(nu_diag,*) '  ', trim(data_file)
      endif

      if (flag) then

      !-----------------------------------------------------------------
      ! Initialize record counters
      ! (n2, n4 will change only at the very beginning or end of
      !  a forcing cycle.)
      !-----------------------------------------------------------------
         n2 = ixm
         n4 = ixp
         arg = 0

      !-----------------------------------------------------------------
      ! read data
      !-----------------------------------------------------------------

         if (ixm /= -99) then
         ! currently in first half of data interval
            if (ixx <= 1) then
               if (yr > fyear_init) then ! get data from previous year
                  !call file_year (data_file, yr-1)
               else             ! yr = fyear_init, no prior data exists
                  if (maxrec > 12) then ! extrapolate from first record
                     if (ixx == 1) n2 = ixx
                  else          ! go to end of fyear_final
                    ! call file_year (data_file, fyear_final)
                  endif
               endif            ! yr > fyear_init
            endif               ! ixx <= 1

      ! write(nu_diag,*) '!! read_data_nc !!!', trim(data_file)
      ! write(nu_diag,*) 'istep  ', istep
      ! write(nu_diag,*) 'fyear_final  ', fyear_final
      ! write(nu_diag,*) 'fyear_init  ', fyear_init
      ! write(nu_diag,*) 'ixm, ixx, ixp  ', ixm, ixx, ixp
      ! write(nu_diag,*) 'maxrec ', maxrec
      ! write(nu_diag,*) 'fieldname  ', fieldname

            call ice_open_nc (data_file, fid)

            arg = 1
            nrec = recd + n2

            call ice_read_nc &
                 (fid, nrec, fieldname, field_data(arg), debug_forcing, &
                  field_loc, field_type)

            !if (ixx==1) call ice_close_nc(fid)
            call ice_close_nc(fid)
         endif                  ! ixm ne -99

        ! always read ixx data from data file for current year
        ! call file_year (data_file, yr)
         call ice_open_nc (data_file, fid)

         arg = arg + 1
         nrec = recd + ixx

         call ice_read_nc &
              (fid, nrec, fieldname, field_data(arg), debug_forcing, &
               field_loc, field_type)

         if (ixp /= -99) then
         ! currently in latter half of data interval
            if (ixx==maxrec) then
               if (yr < fyear_final) then ! get data from following year
                  call ice_close_nc(fid)
                  !call file_year (data_file, yr+1)
                  call ice_open_nc (data_file, fid)
               else             ! yr = fyear_final, no more data exists
                  if (maxrec > 12) then ! extrapolate from ixx
                     n4 = ixx
                  else          ! go to beginning of fyear_init
                     call ice_close_nc(fid)
                    ! call file_year (data_file, fyear_init)
                     call ice_open_nc (data_file, fid)

                  endif
               endif            ! yr < fyear_final
            endif               ! ixx = maxrec

            arg = arg + 1
            nrec = recd + n4

            call ice_read_nc &
                 (fid, nrec, fieldname, field_data(arg), debug_forcing, &
                  field_loc, field_type)
         endif                  ! ixp /= -99

         call ice_close_nc(fid)

      endif                     ! flag

      call ice_timer_stop(timer_readwrite)  ! reading/writing

      end subroutine read_data_nc_point

!=======================================================================

      subroutine ISPOL_files

      character(len=*), parameter :: subname = '(ISPOL_files)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      fsw_file = &
           trim(atm_data_dir)//'/fsw_sfc_4Xdaily.nc'

      flw_file = &
           trim(atm_data_dir)//'/flw_sfc_4Xdaily.nc'

      rain_file = &
           trim(atm_data_dir)//'/fsnow_sfc_daily_mod3.nc'

      uwind_file = &
           trim(atm_data_dir)//'/uatm_10m_daily.nc'

      vwind_file = &
           trim(atm_data_dir)//'/vatm_10m_daily.nc'

      tair_file = &
           trim(atm_data_dir)//'/Tair_2m_daily.nc'

      humid_file = &
           trim(atm_data_dir)//'/Qa_2m_daily.nc'

      if (my_task == master_task) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'Atmospheric data files:'
         write (nu_diag,*) trim(fsw_file)
         write (nu_diag,*) trim(flw_file)
         write (nu_diag,*) trim(rain_file)
         write (nu_diag,*) trim(uwind_file)
         write (nu_diag,*) trim(vwind_file)
         write (nu_diag,*) trim(tair_file)
         write (nu_diag,*) trim(humid_file)
      endif                     ! master_task

      end subroutine ISPOL_files

!=======================================================================

      subroutine ISPOL_data

! Defines atmospheric data fields for Antarctic Weddell sea location

! authors: Nicole Jeffery, LANL
!
      use ice_flux, only: uatm, vatm, Tair, fsw,  Qa, rhoa, &
          frain, fsnow, flw

!local parameters

      character (char_len_long) :: &
         met_file,   &    ! netcdf filename
         fieldname        ! field name in netcdf file

      real (kind=dbl_kind), dimension(2), save :: &
         Tair_data_p      , &      ! air temperature (K) for interpolation
         Qa_data_p,  fsnow_data_p, &
         fsw_data_p, flw_data_p, &
         uatm_data_p, vatm_data_p

      real (kind=dbl_kind), parameter :: & ! coefficients for Hyland-Wexler Qa
         ps1 = 0.58002206e4_dbl_kind,    & ! (K)
         ps2 = 1.3914993_dbl_kind,       & !
         ps3 = 0.48640239e-1_dbl_kind,   & ! (K^-1)
         ps4 = 0.41764768e-4_dbl_kind,   & ! (K^-2)
         ps5 = 0.14452093e-7_dbl_kind,   & ! (K^-3)
         ps6 = 6.5459673_dbl_kind,       & !
         ws1 = 621.97_dbl_kind,          & ! for saturation mixing ratio
         Pair = 1020._dbl_kind,          & ! Sea level pressure (hPa)
         lapse_rate = 0.0065_dbl_kind      ! (K/m) lapse rate over sea level

      ! for interpolation of hourly data
      integer (kind=int_kind) :: &
         ixm,ixx,ixp , &  ! record numbers for neighboring months
         maxrec      , &  ! maximum record number
         recslot     , &  ! spline slot for current record
         dataloc          ! = 1 for data located in middle of time interval
                          ! = 2 for date located at end of time interval
      real (kind=dbl_kind) :: &
         secday    , &
         Qa_pnt

      real (kind=dbl_kind) :: &
         sec1hr           ! number of seconds in 1 hour

      logical (kind=log_kind) :: read1

      integer (kind=int_kind) :: &
          recnum      , & ! record number
          recnum4X        ! record number

      character(len=*), parameter :: subname = '(ISPOL_data)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      call icepack_query_parameters(secday_out=secday)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (trim(atm_data_format) == 'nc') then     ! read nc file

     !-------------------------------------------------------------------
     ! data from NCEP_DOE Reanalysis 2 and Bareiss et al 2008
     ! daily data located at the end of the 24-hour period.
     !-------------------------------------------------------------------

      dataloc = 2                          ! data located at end of interval
      sec1hr = secday                      ! seconds in day
      maxrec = 366                         !

      ! current record number
      recnum = int(yday)

      ! Compute record numbers for surrounding data (2 on each side)
      ixm = mod(recnum+maxrec-2,maxrec) + 1
      ixx = mod(recnum-1,       maxrec) + 1
!     ixp = mod(recnum,         maxrec) + 1

      ! Compute interpolation coefficients
      ! If data is located at the end of the time interval, then the
      !  data value for the current record goes in slot 2

      recslot = 2
      ixp = -99
      call interp_coeff (recnum, recslot, sec1hr, dataloc)

      read1 = .false.
      if (istep==1 .or. oldrecnum .ne. recnum) read1 = .true.

      ! Daily 2m Air temperature 1991

        met_file = tair_file
        fieldname='Tair'

        call read_data_nc_point(read1, 0, fyear, ixm, ixx, ixp, &
                    maxrec, met_file, fieldname, Tair_data_p, &
                    field_loc_center, field_type_scalar)

        Tair(:,:,:) =  c1intp * Tair_data_p(1) &
                       + c2intp * Tair_data_p(2) &
                     - lapse_rate*8.0_dbl_kind

        met_file = humid_file
        fieldname='Qa'

        call read_data_nc_point(read1, 0, fyear, ixm, ixx, ixp, &
                    maxrec, met_file, fieldname, Qa_data_p, &
                    field_loc_center, field_type_scalar)

        Qa_pnt= c1intp * Qa_data_p(1) &
                          + c2intp * Qa_data_p(2)
        Qa(:,:,:) = Qa_pnt

        met_file = uwind_file
        fieldname='uatm'

        call read_data_nc_point(read1, 0, fyear, ixm, ixx, ixp, &
                    maxrec, met_file, fieldname, uatm_data_p, &
                    field_loc_center, field_type_scalar)

        uatm(:,:,:) =  c1intp * uatm_data_p(1) &
                          + c2intp * uatm_data_p(2)

        met_file = vwind_file
        fieldname='vatm'

        call read_data_nc_point(read1, 0, fyear, ixm, ixx, ixp, &
                    maxrec, met_file, fieldname, vatm_data_p, &
                    field_loc_center, field_type_scalar)

        vatm(:,:,:) =  c1intp * vatm_data_p(1) &
                          + c2intp * vatm_data_p(2)

        met_file = rain_file
        fieldname='fsnow'

        call read_data_nc_point(read1, 0, fyear, ixm, ixx, ixp, &
                    maxrec, met_file, fieldname, fsnow_data_p, &
                    field_loc_center, field_type_scalar)

        fsnow(:,:,:) =  (c1intp * fsnow_data_p(1) + &
                         c2intp * fsnow_data_p(2))

        !-----------------------------
        !fsw and flw are every 6 hours
        !------------------------------
        dataloc = 2                          ! data located at end of interval
        sec1hr = secday/c4                   ! seconds in 6 hours
        maxrec = 1460                        ! 366*4

      ! current record number
        recnum4X = 4*int(yday) - 3 + int(real(msec,kind=dbl_kind)/sec1hr)

      ! Compute record numbers for surrounding data (2 on each side)
      ixm = mod(recnum4X+maxrec-2,maxrec) + 1
      ixx = mod(recnum4X-1,       maxrec) + 1

      ! Compute interpolation coefficients
      ! If data is located at the end of the time interval, then the
      !  data value for the current record goes in slot 2

      recslot = 2
      ixp = -99
      call interp_coeff (recnum4X, recslot, sec1hr, dataloc)

      read1 = .false.
      if (istep==1 .or. oldrecnum4X .ne. recnum4X) read1 = .true.

        met_file = fsw_file
        fieldname='fsw'

        call read_data_nc_point(read1, 0, fyear, ixm, ixx, ixp, &
                    maxrec, met_file, fieldname, fsw_data_p, &
                    field_loc_center, field_type_scalar)

        fsw(:,:,:) =  c1intp * fsw_data_p(1) &
                          + c2intp * fsw_data_p(2)

        met_file = flw_file
        fieldname='flw'

        call read_data_nc_point(read1, 0, fyear, ixm, ixx, ixp, &
                    maxrec, met_file, fieldname, flw_data_p, &
                    field_loc_center, field_type_scalar)

        flw(:,:,:) =  c1intp * flw_data_p(1) &
                          + c2intp * flw_data_p(2)
     endif  !nc

      !flw   given cldf and Tair  calculated in prepare_forcing

      !-----------------------------
      ! fixed data
      ! May not be needed
      !-----------------------------
      rhoa (:,:,:) = 1.3_dbl_kind ! air density (kg/m^3)
      cldf(:,:,:) =  c1  !0.25_dbl_kind ! cloud fraction
      frain(:,:,:) = c0            ! this is available in hourlymet_rh file

      ! Save record number for next time step
      oldrecnum = recnum
      oldrecnum4X = recnum4X

      end subroutine ISPOL_data

!=======================================================================

      subroutine ocn_data_ispol_init

! Reads NCAR pop ocean forcing data set 'pop_frc_gx1v3_010815.nc'
! at the ISPOL location -67.4677N, 310.4375E
!
! For ocean mixed layer-----------------------------units
!
! 1  sst------temperature---------------------------(C)
! 2  sss------salinity------------------------------(ppt)
! 3  hbl------depth---------------------------------(m)
! 4  u--------surface u current---------------------(m/s)
! 5  v--------surface v current---------------------(m/s)
! 6  dhdx-----surface tilt x direction--------------(m/m)
! 7  dhdy-----surface tilt y direction--------------(m/m)
! 8  qdp------ocean sub-mixed layer heat flux-------(W/m2)
!
! Fields 4, 5, 6, 7 are on the U-grid; 1, 2, 3, and 8 are
! on the T-grid.
!
! authors: Nicole Jeffery, LANL
!
      use ice_gather_scatter
      use ice_read_write

      integer (kind=int_kind) :: &
        n   , & ! field index
        m       ! month index

      character(char_len) :: &
        vname(nfld) ! variable names to search for in file
      data vname /  &
           'T',      'S',      'hblt',  'U',     'V', &
           'dhdx',   'dhdy',   'qdp' /

      real (kind=dbl_kind) :: &
        work

      integer (kind=int_kind) :: &
        fid         ! file id

      character(len=*), parameter :: subname = '(ocn_data_ispol_init)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      if (my_task == master_task) then

         if (restore_ocn) write (nu_diag,*)  &
             'SST restoring timescale = ',trestore,' days'

         sst_file = trim(ocn_data_dir)//'/'//trim(oceanmixed_file) ! not just sst

        !---------------------------------------------------------------
        ! Read in ocean forcing data from an existing file
        !---------------------------------------------------------------
        write (nu_diag,*) 'ocean mixed layer forcing data file = ', &
                           sst_file

      endif ! master_task

      if (trim(ocn_data_format) == 'nc') then
        if (my_task == master_task) then
          call ice_open_nc(sst_file, fid)
        endif ! master_task

        ! Read in ocean forcing data for all 12 months
        do n=1,nfld
          do m=1,12
            ! Note: netCDF does single to double conversion if necessary
            if (n >= 4 .and. n <= 7) then
               call ice_read_nc(fid, m, vname(n), work, debug_forcing, &
                                field_loc_NEcorner, field_type_vector)
            else
               call ice_read_nc(fid, m, vname(n), work, debug_forcing, &
                                field_loc_center, field_type_scalar)
            endif
            ocn_frc_m(:,:,:,n,m) = work
          enddo               ! month loop
        enddo               ! field loop

        if (my_task == master_task) call ice_close_nc(fid)

      else  ! binary format
         call abort_ice (error_message=subname//'new ocean forcing is netcdf only', &
            file=__FILE__, line=__LINE__)
      endif

!echmod - currents cause Fram outflow to be too large
              ocn_frc_m(:,:,:,4,:) = c0
              ocn_frc_m(:,:,:,5,:) = c0
!echmod

      end subroutine ocn_data_ispol_init

!=======================================================================
!
      subroutine box2001_data_atm

! wind fields as in Hunke, JCP 2001
! these are defined at the u point
! authors: Elizabeth Hunke, LANL

      use ice_domain, only: nblocks, blocks_ice
      use ice_calendar, only: timesecs
      use ice_blocks, only: block, get_block, nx_block, ny_block, nghost
      use ice_flux, only: uatm, vatm, wind, rhoa, strax, stray
      use ice_state, only: aice

      ! local parameters

      integer (kind=int_kind) :: &
         iblk, i,j           ! loop indices

      integer (kind=int_kind) :: &
         iglob(nx_block), & ! global indices
         jglob(ny_block)    ! global indices

      type (block) :: &
         this_block           ! block information for current block

      real (kind=dbl_kind) :: &
         secday, pi , puny, period, pi2, tau

      character(len=*), parameter :: subname = '(box2001_data_atm)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      call icepack_query_parameters(pi_out=pi, pi2_out=pi2, puny_out=puny)
      call icepack_query_parameters(secday_out=secday)

      period = c4*secday

      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block

         this_block = get_block(blocks_ice(iblk),iblk)
         iglob = this_block%i_glob
         jglob = this_block%j_glob

         ! wind components
         uatm(i,j,iblk) = c5 + (sin(pi2*timesecs/period)-c3) &
                              * sin(pi2*real(iglob(i), kind=dbl_kind)  &
                                       /real(nx_global,kind=dbl_kind)) &
                              * sin(pi *real(jglob(j), kind=dbl_kind)  &
                                       /real(ny_global,kind=dbl_kind))
         vatm(i,j,iblk) = c5 + (sin(pi2*timesecs/period)-c3) &
                              * sin(pi *real(iglob(i), kind=dbl_kind)  &
                                       /real(nx_global,kind=dbl_kind)) &
                              * sin(pi2*real(jglob(j), kind=dbl_kind)  &
                                       /real(ny_global,kind=dbl_kind))
         ! wind stress
         wind(i,j,iblk) = sqrt(uatm(i,j,iblk)**2 + vatm(i,j,iblk)**2)
         tau = rhoa(i,j,iblk) * 0.0012_dbl_kind * wind(i,j,iblk)

         strax(i,j,iblk) = aice(i,j,iblk) * tau * uatm(i,j,iblk)
         stray(i,j,iblk) = aice(i,j,iblk) * tau * vatm(i,j,iblk)

! initialization test
       ! Diagonal wind vectors 1
         !uatm(i,j,iblk) = c1 *real(j-nghost, kind=dbl_kind) &
         !                   / real(ny_global,kind=dbl_kind)
         !vatm(i,j,iblk) = c1 *real(j-nghost, kind=dbl_kind) &
         !                   / real(ny_global,kind=dbl_kind)

       ! Diagonal wind vectors 2
         !uatm(i,j,iblk) = c1 *real(i-nghost, kind=dbl_kind) &
         !                   / real(nx_global,kind=dbl_kind)
         !vatm(i,j,iblk) = -c1 *real(i-nghost, kind=dbl_kind) &
         !                   / real(nx_global,kind=dbl_kind)

       ! Wind in x direction
        ! uatm(i,j,iblk) = c1 *real(i-nghost, kind=dbl_kind) &
        !                    / real(nx_global,kind=dbl_kind)
        ! vatm(i,j,iblk) = c0

       ! Wind in y direction
        ! uatm(i,j,iblk) = c0
        ! vatm(i,j,iblk) = c1 *real(j-nghost, kind=dbl_kind) &
        !                    / real(ny_global,kind=dbl_kind)
! initialization test

         enddo
         enddo
      enddo ! nblocks

      end subroutine box2001_data_atm

!=======================================================================
!
      subroutine box2001_data_ocn

! current fields as in Hunke, JCP 2001
! these are defined at the u point
! authors: Elizabeth Hunke, LANL

      use ice_domain, only: nblocks, blocks_ice
      use ice_blocks, only: block, get_block, nx_block, ny_block, nghost
      use ice_flux, only: uocn, vocn

      ! local parameters

      integer (kind=int_kind) :: &
         iblk, i,j           ! loop indices

      integer (kind=int_kind) :: &
         iglob(nx_block), & ! global indices
         jglob(ny_block)    ! global indices

      type (block) :: &
         this_block           ! block information for current block

      character(len=*), parameter :: subname = '(box2001_data_ocn)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block

         this_block = get_block(blocks_ice(iblk),iblk)
         iglob = this_block%i_glob
         jglob = this_block%j_glob

         ! ocean current
         ! constant in time, could be initialized in ice_flux.F90
         uocn(i,j,iblk) =  p2*real(jglob(j), kind=dbl_kind) &
                            / real(ny_global,kind=dbl_kind) - p1
         vocn(i,j,iblk) = -p2*real(iglob(i), kind=dbl_kind) &
                            / real(nx_global,kind=dbl_kind) + p1

         enddo
         enddo
      enddo ! nblocks

      end subroutine box2001_data_ocn

!=======================================================================
!
      subroutine uniform_data_atm(dir,spd)
!     uniform wind fields in some direction

      use ice_domain, only: nblocks
      use ice_blocks, only: nx_block, ny_block, nghost
      use ice_flux, only: uatm, vatm, wind, rhoa, strax, stray
      use ice_state, only: aice

      character(len=*), intent(in) :: dir
      real(kind=dbl_kind), intent(in), optional :: spd ! velocity

      ! local parameters

      integer (kind=int_kind) :: &
         iblk, i,j           ! loop indices

      real (kind=dbl_kind) :: &
         tau, &
         atm_val ! value to use for atm speed

      character(len=*), parameter :: subname = '(uniform_data_atm)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      ! check for optional spd
      if (present(spd)) then
         atm_val = spd
      else
         atm_val = c5 ! default
      endif

      ! wind components
      if (dir == 'NE') then
         uatm = atm_val
         vatm = atm_val
      elseif (dir == 'N') then
         uatm = c0
         vatm = atm_val
      elseif (dir == 'E') then
         uatm = atm_val
         vatm = c0
      elseif (dir == 'S') then
         uatm = c0
         vatm = -atm_val
      elseif (dir == 'W') then
         uatm = -atm_val
         vatm = c0
      else
         call abort_ice (subname//'ERROR: dir unknown, dir = '//trim(dir), &
              file=__FILE__, line=__LINE__)
      endif

      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block

            ! wind stress
            wind(i,j,iblk) = sqrt(uatm(i,j,iblk)**2 + vatm(i,j,iblk)**2)
            tau = rhoa(i,j,iblk) * 0.0012_dbl_kind * wind(i,j,iblk)
            strax(i,j,iblk) = aice(i,j,iblk) * tau * uatm(i,j,iblk)
            stray(i,j,iblk) = aice(i,j,iblk) * tau * vatm(i,j,iblk)

         enddo
         enddo
      enddo ! nblocks

      end subroutine uniform_data_atm
!=======================================================================

!
      subroutine uniform_data_ocn(dir,spd)

!     uniform current fields in some direction

      use ice_flux, only: uocn, vocn

      character(len=*), intent(in) :: dir

      real(kind=dbl_kind), intent(in), optional :: spd ! velocity

      ! local parameters

      real(kind=dbl_kind) :: &
           ocn_val ! value to use for ocean currents

      character(len=*), parameter :: subname = '(uniform_data_ocn)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      if (present(spd)) then
         ocn_val = spd
      else
         ocn_val = p1 ! default
      endif

      ! ocn components
      if (dir == 'NE') then
         uocn = ocn_val
         vocn = ocn_val
      elseif (dir == 'N') then
         uocn = c0
         vocn = ocn_val
      elseif (dir == 'E') then
         uocn = ocn_val
         vocn = c0
      else
         call abort_ice (subname//'ERROR: dir unknown, dir = '//trim(dir), &
              file=__FILE__, line=__LINE__)
      endif

      end subroutine uniform_data_ocn
!=======================================================================

      subroutine get_wave_spec

      use ice_read_write, only: ice_read_nc_xyf
      use ice_arrays_column, only: wave_spectrum, &
                                   dwavefreq, wavefreq
      use ice_constants, only: c0
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_fsd

      ! local variables
      integer (kind=int_kind) :: &
         fid                    ! file id for netCDF routines

      real(kind=dbl_kind), dimension(nfreq) :: &
         wave_spectrum_profile  ! wave spectrum

      character(char_len) :: wave_spec_type
      logical (kind=log_kind) :: wave_spec
      character(len=*), parameter :: subname = '(get_wave_spec)'

      if (local_debug .and. my_task == master_task) write(nu_diag,*) subname,'fdbg start'

      call ice_timer_start(timer_fsd)

      call icepack_query_parameters(wave_spec_out=wave_spec, &
                                    wave_spec_type_out=wave_spec_type)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      ! if no wave data is provided, wave_spectrum is zero everywhere
      wave_spectrum(:,:,:,:) = c0
      debug_forcing = .false.

      ! wave spectrum and frequencies
      if (wave_spec) then
      ! get hardwired frequency bin info and a dummy wave spectrum profile
      ! the latter is used if wave_spec_type == profile
         call icepack_init_wave(nfreq     = nfreq,    &
                                wave_spectrum_profile = wave_spectrum_profile, &
                                wavefreq  = wavefreq, &
                                dwavefreq = dwavefreq)

         ! read more realistic data from a file
         if ((trim(wave_spec_type) == 'constant').OR.(trim(wave_spec_type) == 'random')) then
            if (trim(wave_spec_file(1:4)) == 'unkn') then
               call abort_ice (subname//'ERROR: wave_spec_file '//trim(wave_spec_file), &
                  file=__FILE__, line=__LINE__)
            else
#ifdef USE_NETCDF
               call wave_spec_data
#else
               write (nu_diag,*) "wave spectrum file not available, requires cpp USE_NETCDF"
               write (nu_diag,*) "wave spectrum file not available, using default profile"
               call abort_ice (subname//'ERROR: wave_spec_file '//trim(wave_spec_file), &
                  file=__FILE__, line=__LINE__)
#endif
            endif
         endif
      endif

      call ice_timer_stop(timer_fsd)

      end subroutine get_wave_spec

!=======================================================================
!
!   Read in wave spectrum forcing as a function of time. 6 hourly
!   LR started working from JRA55_data routine
!   Changed fields, and changed 3 hourly to 6 hourly
!
      subroutine wave_spec_data

      use ice_blocks, only: block, get_block
      use ice_global_reductions, only: global_minval, global_maxval
      use ice_domain, only: nblocks, distrb_info, blocks_ice
      use ice_arrays_column, only: wave_spectrum, &
                                   dwavefreq, wavefreq
      use ice_read_write, only: ice_read_nc_xyf
      use ice_grid, only: hm, tlon, tlat, tmask, umask
      use ice_calendar, only: days_per_year, use_leap_years

      integer (kind=int_kind) :: &
          ncid        , & ! netcdf file id
          i, j, freq  , &
          ixm,ixx,ixp , & ! record numbers for neighboring months
          recnum      , & ! record number
          maxrec      , & ! maximum record number
          recslot     , & ! spline slot for current record
          midmonth    , & ! middle day of month
          dataloc     , & ! = 1 for data located in middle of time interval
                          ! = 2 for date located at end of time interval
          iblk        , & ! block index
          ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
          yr              ! current forcing year

      real (kind=dbl_kind) :: &
          sec6hr          , & ! number of seconds in 3 hours
          secday          , & ! number of seconds in day
          vmin, vmax

      logical (kind=log_kind) :: readm, read6,debug_n_d

      type (block) :: &
         this_block           ! block information for current block

      real(kind=dbl_kind), dimension(nfreq) :: &
         wave_spectrum_profile  ! wave spectrum

      character(len=64) :: fieldname !netcdf field name
      character(char_len_long) :: spec_file
      character(char_len) :: wave_spec_type
      logical (kind=log_kind) :: wave_spec
      character(len=*), parameter :: subname = '(wave_spec_data)'

      debug_n_d = .false.  !usually false

      call icepack_query_parameters(secday_out=secday)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

         call icepack_init_wave(nfreq     = nfreq,    &
                                wave_spectrum_profile = wave_spectrum_profile, &
                                wavefreq  = wavefreq, &
                                dwavefreq = dwavefreq)

      !spec_file = trim(ocn_data_dir)//'/'//trim(wave_spec_file)
      spec_file = trim(wave_spec_file)
      wave_spectrum_data = c0
      wave_spectrum = c0
      yr = fyear  ! current year
    !-------------------------------------------------------------------
    ! 6-hourly data
    !
    ! Assume that the 6-hourly value is located at the end of the
    !  6-hour period.  This is the convention for NCEP reanalysis data.
    !  E.g. record 1 gives conditions at 6 am GMT on 1 January.
    !-------------------------------------------------------------------

      dataloc = 2               ! data located at end of interval
      sec6hr = secday/c4        ! seconds in 6 hours
      !maxrec = 2920            ! 365*8; for leap years = 366*8

      if (use_leap_years) days_per_year = 366 !overrides setting of 365 in ice_calendar
      maxrec = days_per_year*4

      if(days_per_year == 365 .and. (mod(yr,  4) == 0)) then
      call abort_ice('days_per_year should be set to 366 for leap years')
      end if

      ! current record number
      recnum = 4*int(yday) - 3 + int(real(msec,kind=dbl_kind)/sec6hr)

      ! Compute record numbers for surrounding data (2 on each side)

      ixm = mod(recnum+maxrec-2,maxrec) + 1
      ixx = mod(recnum-1,       maxrec) + 1

      ! Compute interpolation coefficients
      ! If data is located at the end of the time interval, then the
      !  data value for the current record goes in slot 2

      recslot = 2
      ixp = -99
      call interp_coeff (recnum, recslot, sec6hr, dataloc)

      ! Read
      read6 = .false.
      if (istep==1 .or. oldrecnum .ne. recnum) read6 = .true.
         !-------------------------------------------------------------------
         ! File is NETCDF
         ! file variable names are:
         ! efreq   (wave spectrum, energy as a function of wave frequency UNITS)
         !-------------------------------------------------------------------
         call ice_open_nc(spec_file,ncid)

         call ice_read_nc_xyf(ncid,recnum,'efreq',wave_spectrum_data(:,:,:,1,:),debug_n_d, &
              field_loc=field_loc_center, &
              field_type=field_type_scalar)
         call ice_read_nc_xyf(ncid,recnum,'efreq',wave_spectrum_data(:,:,:,2,:),debug_n_d, &
              field_loc=field_loc_center, &
              field_type=field_type_scalar)
         call ice_close_nc(ncid)


      ! Interpolate
      call interpolate_wavespec_data (wave_spectrum_data, wave_spectrum)

      ! Save record number
      oldrecnum = recnum

         if (local_debug) then
           if (my_task == master_task) write (nu_diag,*) &
              'wave_spec_data ',spec_file
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'maxrec',maxrec
               write (nu_diag,*) 'days_per_year', days_per_year

        endif                   ! local debug

      end subroutine wave_spec_data

!=======================================================================

! initial snow aging lookup table
!
! Dry snow metamorphism table
! snicar_drdt_bst_fit_60_c070416.nc
! Flanner (file metadata units mislabelled)
! drdsdt0 (10^-6 m/hr) tau (10^-6 m)
!
     subroutine init_snowtable

      use ice_broadcast, only: broadcast_array, broadcast_scalar
      integer (kind=int_kind) :: &
          idx_T_max   , &  ! Table dimensions
          idx_rhos_max, &
          idx_Tgrd_max
      real (kind=dbl_kind), allocatable :: &
          snowage_rhos (:), &
          snowage_Tgrd (:), &
          snowage_T    (:), &
          snowage_tau  (:,:,:), &
          snowage_kappa(:,:,:), &
          snowage_drdt0(:,:,:)

      ! local variables

      logical (kind=log_kind) :: diag = .false.

      integer (kind=int_kind) :: &
         fid                  ! file id for netCDF file

      character (char_len) :: &
         snw_aging_table, &   ! aging table setting
         fieldname            ! field name in netcdf file

      character(len=*), parameter :: subname = '(init_snowtable)'

      !-----------------------------------------------------------------
      ! read table of snow aging parameters
      !-----------------------------------------------------------------

      call icepack_query_parameters(snw_aging_table_out=snw_aging_table, &
         isnw_rhos_out=idx_rhos_max, isnw_Tgrd_out=idx_Tgrd_max, isnw_T_out=idx_T_max)

      if (my_task == master_task) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'Snow aging file:', trim(snw_filename)
      endif

      if (snw_aging_table == 'snicar') then
         ! just read the 3d data and pass it in

         call ice_open_nc(snw_filename,fid)

         allocate(snowage_tau  (idx_rhos_max, idx_Tgrd_max, idx_T_max))
         allocate(snowage_kappa(idx_rhos_max, idx_Tgrd_max, idx_T_max))
         allocate(snowage_drdt0(idx_rhos_max, idx_Tgrd_max, idx_T_max))

         fieldname = trim(snw_tau_fname)
         call ice_read_nc(fid,fieldname,snowage_tau,  diag, &
                          idx_rhos_max,idx_Tgrd_max,idx_T_max)
         fieldname = trim(snw_kappa_fname)
         call ice_read_nc(fid,fieldname,snowage_kappa,diag, &
                          idx_rhos_max,idx_Tgrd_max,idx_T_max)
         fieldname = trim(snw_drdt0_fname)
         call ice_read_nc(fid,fieldname,snowage_drdt0,diag, &
                          idx_rhos_max,idx_Tgrd_max,idx_T_max)

         call ice_close_nc(fid)

         call broadcast_array(snowage_tau  , master_task)
         call broadcast_array(snowage_kappa, master_task)
         call broadcast_array(snowage_drdt0, master_task)

         if (my_task == master_task) then
            write(nu_diag,*) subname,'  '
            write(nu_diag,*) subname,' Successfully read snow aging properties:'
            write(nu_diag,*) subname,' snw_aging_table = ',trim(snw_aging_table)
            write(nu_diag,*) subname,' idx_rhos_max = ',idx_rhos_max
            write(nu_diag,*) subname,' idx_Tgrd_max = ',idx_Tgrd_max
            write(nu_diag,*) subname,' idx_T_max    = ',idx_T_max
            write(nu_diag,*) subname,' Data at rhos, Tgrd, T at first index '
            write(nu_diag,*) subname,' snoage_tau (1,1,1)         = ',snowage_tau  (1,1,1)
            write(nu_diag,*) subname,' snoage_kappa (1,1,1)       = ',snowage_kappa(1,1,1)
            write(nu_diag,*) subname,' snoage_drdt0 (1,1,1)       = ',snowage_drdt0(1,1,1)
            write(nu_diag,*) subname,' Data at rhos, Tgrd, T at max index'
            write(nu_diag,*) subname,' snoage_tau (max,max,max)   = ',snowage_tau  (idx_rhos_max, idx_Tgrd_max, idx_T_max)
            write(nu_diag,*) subname,' snoage_kappa (max,max,max) = ',snowage_kappa(idx_rhos_max, idx_Tgrd_max, idx_T_max)
            write(nu_diag,*) subname,' snoage_drdt0 (max,max,max) = ',snowage_drdt0(idx_rhos_max, idx_Tgrd_max, idx_T_max)
         endif

         call icepack_init_parameters(        &
            snowage_tau_in   = snowage_tau,   &
            snowage_kappa_in = snowage_kappa, &
            snowage_drdt0_in = snowage_drdt0 )

         deallocate(snowage_tau)
         deallocate(snowage_kappa)
         deallocate(snowage_drdt0)

      else
         ! read everything and pass it in

         call ice_open_nc(snw_filename,fid)

         fieldname = trim(snw_rhos_fname)
         call ice_get_ncvarsize(fid,fieldname,idx_rhos_max)
         fieldname = trim(snw_Tgrd_fname)
         call ice_get_ncvarsize(fid,fieldname,idx_Tgrd_max)
         fieldname = trim(snw_T_fname)
         call ice_get_ncvarsize(fid,fieldname,idx_T_max)

         call broadcast_scalar(idx_rhos_max, master_task)
         call broadcast_scalar(idx_Tgrd_max, master_task)
         call broadcast_scalar(idx_T_max   , master_task)

         allocate(snowage_rhos (idx_rhos_max))
         allocate(snowage_Tgrd (idx_Tgrd_max))
         allocate(snowage_T    (idx_T_max))
         allocate(snowage_tau  (idx_rhos_max, idx_Tgrd_max, idx_T_max))
         allocate(snowage_kappa(idx_rhos_max, idx_Tgrd_max, idx_T_max))
         allocate(snowage_drdt0(idx_rhos_max, idx_Tgrd_max, idx_T_max))

         fieldname = trim(snw_rhos_fname)
         call ice_read_nc(fid,fieldname,snowage_rhos,  diag, &
                          idx_rhos_max)
         fieldname = trim(snw_Tgrd_fname)
         call ice_read_nc(fid,fieldname,snowage_Tgrd,  diag, &
                          idx_Tgrd_max)
         fieldname = trim(snw_T_fname)
         call ice_read_nc(fid,fieldname,snowage_T,  diag, &
                          idx_T_max)

         fieldname = trim(snw_tau_fname)
         call ice_read_nc(fid,fieldname,snowage_tau,  diag, &
                          idx_rhos_max,idx_Tgrd_max,idx_T_max)
         fieldname = trim(snw_kappa_fname)
         call ice_read_nc(fid,fieldname,snowage_kappa,diag, &
                          idx_rhos_max,idx_Tgrd_max,idx_T_max)
         fieldname = trim(snw_drdt0_fname)
         call ice_read_nc(fid,fieldname,snowage_drdt0,diag, &
                          idx_rhos_max,idx_Tgrd_max,idx_T_max)

         call ice_close_nc(fid)

         call broadcast_array(snowage_rhos , master_task)
         call broadcast_array(snowage_Tgrd , master_task)
         call broadcast_array(snowage_T    , master_task)
         call broadcast_array(snowage_tau  , master_task)
         call broadcast_array(snowage_kappa, master_task)
         call broadcast_array(snowage_drdt0, master_task)

         if (my_task == master_task) then
            write(nu_diag,*) subname,'  '
            write(nu_diag,*) subname,' Successfully read snow aging properties:'
            write(nu_diag,*) subname,' idx_rhos_max = ',idx_rhos_max
            write(nu_diag,*) subname,' idx_Tgrd_max = ',idx_Tgrd_max
            write(nu_diag,*) subname,' idx_T_max    = ',idx_T_max
            write(nu_diag,*) subname,' Data at rhos, Tgrd, T = ',snowage_rhos(1),snowage_Tgrd(1),snowage_T(1)
            write(nu_diag,*) subname,' snoage_tau (1,1,1)         = ',snowage_tau  (1,1,1)
            write(nu_diag,*) subname,' snoage_kappa (1,1,1)       = ',snowage_kappa(1,1,1)
            write(nu_diag,*) subname,' snoage_drdt0 (1,1,1)       = ',snowage_drdt0(1,1,1)
            write(nu_diag,*) subname,' Data at rhos, Tgrd, T = ', &
                                       snowage_rhos(idx_rhos_max),snowage_Tgrd(idx_Tgrd_max),snowage_T(idx_T_max)
            write(nu_diag,*) subname,' snoage_tau (max,max,max)   = ',snowage_tau  (idx_rhos_max, idx_Tgrd_max, idx_T_max)
            write(nu_diag,*) subname,' snoage_kappa (max,max,max) = ',snowage_kappa(idx_rhos_max, idx_Tgrd_max, idx_T_max)
            write(nu_diag,*) subname,' snoage_drdt0 (max,max,max) = ',snowage_drdt0(idx_rhos_max, idx_Tgrd_max, idx_T_max)
         endif

         call icepack_init_parameters(        &
            isnw_t_in        = idx_T_max,     &
            isnw_Tgrd_in     = idx_Tgrd_max,  &
            isnw_rhos_in     = idx_rhos_max,  &
            snowage_rhos_in  = snowage_rhos,  &
            snowage_Tgrd_in  = snowage_Tgrd,  &
            snowage_T_in     = snowage_T,     &
            snowage_tau_in   = snowage_tau,   &
            snowage_kappa_in = snowage_kappa, &
            snowage_drdt0_in = snowage_drdt0 )

         deallocate(snowage_rhos)
         deallocate(snowage_Tgrd)
         deallocate(snowage_T)
         deallocate(snowage_tau)
         deallocate(snowage_kappa)
         deallocate(snowage_drdt0)

      endif

      end subroutine init_snowtable

!=======================================================================

      end module ice_forcing

!=======================================================================
