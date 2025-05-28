#ifdef ncdf
#define USE_NETCDF
#endif
!=======================================================================
!
! Writes history in netCDF format
!
! authors Tony Craig and Bruce Briegleb, NCAR
!         Elizabeth C. Hunke and William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! 2004 WHL: Block structure added
! 2006 ECH: Accepted some CESM code into mainstream CICE
!           Added ice_present, aicen, vicen; removed aice1...10, vice1...1.
!           Added histfreq_n and histfreq='h' options, removed histfreq='w'
!           Converted to free source form (F90)
!           Added option for binary output instead of netCDF
! 2009 D Bailey and ECH: Generalized for multiple frequency output
! 2010 Alison McLaren and ECH: Added 3D capability
! 2013 ECH split from ice_history.F90

      module ice_history_write

      use ice_constants, only: c0, c360, p5, spval, spval_dbl
      use ice_fileunits, only: nu_diag
      use ice_exit, only: abort_ice
      use ice_read_write, only: ice_check_nc
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_parameters
      use ice_kinds_mod, only: int_kind
#ifdef USE_NETCDF
      use netcdf
#endif

      implicit none
      private

      TYPE coord_attributes         ! netcdf coordinate attributes
        character (len=11)   :: short_name
        character (len=45)   :: long_name
        character (len=30)   :: units
        character (len=8)    :: axis
      END TYPE coord_attributes

      TYPE req_attributes         ! req'd netcdf attributes
        type (coord_attributes) :: req
        character (len=20)   :: coordinates
      END TYPE req_attributes

      public :: ice_write_hist

      integer (kind=int_kind) :: imtid,jmtid

!=======================================================================

      contains

!=======================================================================
!
! write average ice quantities or snapshots
!
! author:   Elizabeth C. Hunke, LANL

      subroutine ice_write_hist (ns)

      use ice_kinds_mod
      use ice_arrays_column, only: hin_max, floe_rad_c
      use ice_blocks, only: nx_block, ny_block
      use ice_broadcast, only: broadcast_scalar
      use ice_calendar, only: msec, timesecs, idate, idate0, write_ic, &
          histfreq, histfreq_n, days_per_year, use_leap_years, dayyr, &
          hh_init, mm_init, ss_init
      use ice_communicate, only: my_task, master_task
      use ice_domain, only: distrb_info
      use ice_domain_size, only: nx_global, ny_global, max_nstrm, max_blocks
      use ice_gather_scatter, only: gather_global
      use ice_grid, only: TLON, TLAT, ULON, ULAT, NLON, NLAT, ELON, ELAT, &
          hm, uvm, npm, epm, bm, tarea, uarea, narea, earea, &
          dxU, dxT, dyU, dyT, dxN, dyN, dxE, dyE, HTN, HTE, ANGLE, ANGLET, &
          lont_bounds, latt_bounds, lonu_bounds, latu_bounds, &
          lonn_bounds, latn_bounds, lone_bounds, late_bounds
      use ice_history_shared
#ifdef CESMCOUPLED
      use ice_restart_shared, only: runid
#endif

      integer (kind=int_kind), intent(in) :: ns

      ! local variables

      real (kind=dbl_kind), dimension(:,:),   allocatable :: work_g1
      real (kind=dbl_kind), dimension(:,:,:), allocatable :: work1_3
      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: work1

      integer (kind=int_kind) :: i,k,ic,n,nn, &
         ncid,status,kmtidi,kmtids,kmtidb, cmtid,timid,varid, &
         nvertexid,ivertex,kmtida,iflag, fmtid
      integer (kind=int_kind), dimension(3) :: dimid
      integer (kind=int_kind), dimension(4) :: dimidz
      integer (kind=int_kind), dimension(5) :: dimidcz
      integer (kind=int_kind), dimension(3) :: dimid_nverts
      integer (kind=int_kind), dimension(6) :: dimidex
      real (kind=dbl_kind)  :: ltime2
      character (char_len) :: title, cal_units, cal_att
      character (char_len) :: time_period_freq = 'none'
      character (char_len_long) :: ncfile
      character (len=512) :: extvars
      real (kind=dbl_kind)  :: secday, rad_to_deg

      integer (kind=int_kind) :: ind,boundid, lprecision

      character (char_len) :: start_time,current_date,current_time
      character (len=8) :: cdate

      ! time coord
      TYPE(coord_attributes) :: time_coord

      ! 4 vertices in each grid cell
      INTEGER (kind=int_kind), PARAMETER :: nverts = 4

      ! 8 variables describe T, U grid boundaries:
      ! lont_bounds, latt_bounds, lonu_bounds, latu_bounds
      INTEGER (kind=int_kind), PARAMETER :: nvar_verts = 8

      TYPE(req_attributes), dimension(nvar_grd) :: var_grd
      TYPE(coord_attributes), dimension(ncoord) :: var_coord
      TYPE(coord_attributes), dimension(nvar_verts) :: var_nverts
      TYPE(coord_attributes), dimension(nvar_grdz) :: var_grdz
      CHARACTER (char_len), dimension(ncoord) :: coord_bounds

      character(len=*), parameter :: subname = '(ice_write_hist)'

#ifdef USE_NETCDF
      call icepack_query_parameters(secday_out=secday, rad_to_deg_out=rad_to_deg)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__, line=__LINE__)

      extvars = ''
      lprecision = nf90_float
      if (history_precision == 8) lprecision = nf90_double

      if (my_task == master_task) then

         call construct_filename(ncfile,'nc',ns)

         ! add local directory path name to ncfile
         if (write_ic) then
            ncfile = trim(incond_dir)//ncfile
         else
            ncfile = trim(history_dir)//ncfile
         endif

         ! create file
         if (history_format == 'cdf1') then
           iflag = nf90_clobber
         elseif (history_format == 'cdf2') then
#ifdef NO_CDF2
           call abort_ice(subname//' ERROR: history_format cdf2 not available ', &
              file=__FILE__, line=__LINE__)
#else
           iflag = ior(nf90_clobber,nf90_64bit_offset)
#endif
         elseif (history_format == 'cdf5') then
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
           call abort_ice(subname//' ERROR: history_format not allowed for '//trim(history_format), &
              file=__FILE__, line=__LINE__)
         endif
         status = nf90_create(ncfile, iflag, ncid)
         call ice_check_nc(status, subname// ' ERROR: creating history ncfile '//ncfile, &
                           file=__FILE__, line=__LINE__)

         !-----------------------------------------------------------------
         ! define dimensions
         !-----------------------------------------------------------------

         if (hist_avg(ns) .and. .not. write_ic) then
            status = nf90_def_dim(ncid,'nbnd',2,boundid)
            call ice_check_nc(status, subname// ' ERROR: defining dim nbnd', &
                              file=__FILE__, line=__LINE__)
         endif

         status = nf90_def_dim(ncid,'ni',nx_global,imtid)
         call ice_check_nc(status, subname// ' ERROR: defining dim ni', &
                           file=__FILE__, line=__LINE__)

         status = nf90_def_dim(ncid,'nj',ny_global,jmtid)
         call ice_check_nc(status, subname// ' ERROR: defining dim nj', &
                           file=__FILE__, line=__LINE__)

         status = nf90_def_dim(ncid,'nc',ncat_hist,cmtid)
         call ice_check_nc(status, subname// ' ERROR: defining dim nc', &
                           file=__FILE__, line=__LINE__)

         status = nf90_def_dim(ncid,'nkice',nzilyr,kmtidi)
         call ice_check_nc(status, subname// ' ERROR: defining dim nkice', &
                           file=__FILE__, line=__LINE__)

         status = nf90_def_dim(ncid,'nksnow',nzslyr,kmtids)
         call ice_check_nc(status, subname// ' ERROR: defining dim nksnow', &
                           file=__FILE__, line=__LINE__)

         status = nf90_def_dim(ncid,'nkbio',nzblyr,kmtidb)
         call ice_check_nc(status, subname// ' ERROR: defining dim nkbio', &
                           file=__FILE__, line=__LINE__)

         status = nf90_def_dim(ncid,'nkaer',nzalyr,kmtida)
         call ice_check_nc(status, subname// ' ERROR: defining dim nkaer', &
                           file=__FILE__, line=__LINE__)

         ! do not write time axis on grid output file
         timid = -99
         if (histfreq(ns)/='g') then
            status = nf90_def_dim(ncid,'time',NF90_UNLIMITED,timid)
            call ice_check_nc(status, subname// ' ERROR: defining dim time', &
                              file=__FILE__, line=__LINE__)
         endif

         status = nf90_def_dim(ncid,'nvertices',nverts,nvertexid)
         call ice_check_nc(status, subname// ' ERROR: defining dim nvertices', &
                           file=__FILE__, line=__LINE__)

         status = nf90_def_dim(ncid,'nf',nfsd_hist,fmtid)
         call ice_check_nc(status, subname// ' ERROR: defining dim nf', &
                           file=__FILE__, line=__LINE__)

         !-----------------------------------------------------------------
         ! define coordinate variables: time, time_bounds
         !-----------------------------------------------------------------

         ! do not write time axis on grid output file
         if (histfreq(ns)/='g') then

            write(cdate,'(i8.8)') idate0
            write(cal_units,'(a,a4,a1,a2,a1,a2,a1,i2.2,a1,i2.2,a1,i2.2)') 'days since ', &
                  cdate(1:4),'-',cdate(5:6),'-',cdate(7:8),' ', &
                  hh_init,':',mm_init,':',ss_init

            if (days_per_year == 360) then
               cal_att='360_day'
            elseif (days_per_year == 365 .and. .not.use_leap_years ) then
               cal_att='noleap'
            elseif (use_leap_years) then
               cal_att='proleptic_gregorian'
            else
               call abort_ice(subname//' ERROR: invalid calendar settings', file=__FILE__, line=__LINE__)
            endif

            time_coord = coord_attributes('time', 'time', trim(cal_units),'T')
            call ice_hist_coord_def(ncid, time_coord, nf90_double, (/timid/), varid)

            status = nf90_put_att(ncid,varid,'calendar',cal_att) !extra attribute
            call ice_check_nc(status,  subname//' ERROR: defining att calendar: '//cal_att,file=__FILE__,line=__LINE__)
            if (hist_avg(ns) .and. .not. write_ic) then
               status = nf90_put_att(ncid,varid,'bounds','time_bounds')
               call ice_check_nc(status, subname//' ERROR: defining att bounds time_bounds',file=__FILE__,line=__LINE__)
            endif

            ! Define coord time_bounds if hist_avg is true
            ! bounds inherit attributes
            if (hist_avg(ns) .and. .not. write_ic) then
               time_coord = coord_attributes('time_bounds', 'undefined', 'undefined', 'undefined')

               dimid(1) = boundid
               dimid(2) = timid

               call ice_hist_coord_def(ncid, time_coord, nf90_double, dimid(1:2), varid)
            endif

         endif  ! histfreq(ns)/='g'

         !-----------------------------------------------------------------
         ! define information for required time-invariant variables
         !-----------------------------------------------------------------

         do ind = 1, ncoord
            select case (ind)
               case(n_tlon)
                  var_coord(ind) = coord_attributes('TLON', &
                                   'T grid center longitude', 'degrees_east', 'X')
                  coord_bounds(ind) = 'lont_bounds'
               case(n_tlat)
                  var_coord(ind) = coord_attributes('TLAT', &
                                   'T grid center latitude',  'degrees_north', 'Y')
                  coord_bounds(ind) = 'latt_bounds'
               case(n_ulon)
                  var_coord(ind) = coord_attributes('ULON', &
                                   'U grid center longitude', 'degrees_east', 'X')
                  coord_bounds(ind) = 'lonu_bounds'
               case(n_ulat)
                  var_coord(ind) = coord_attributes('ULAT', &
                                   'U grid center latitude',  'degrees_north', 'Y')
                  coord_bounds(ind) = 'latu_bounds'
               case(n_nlon)
                  var_coord(ind) = coord_attributes('NLON', &
                                   'N grid center longitude', 'degrees_east', 'X')
                  coord_bounds(ind) = 'lonn_bounds'
               case(n_nlat)
                  var_coord(ind) = coord_attributes('NLAT', &
                                   'N grid center latitude',  'degrees_north', 'Y')
                  coord_bounds(ind) = 'latn_bounds'
               case(n_elon)
                  var_coord(ind) = coord_attributes('ELON', &
                                   'E grid center longitude', 'degrees_east', 'X')
                  coord_bounds(ind) = 'lone_bounds'
               case(n_elat)
                  var_coord(ind) = coord_attributes('ELAT', &
                                   'E grid center latitude',  'degrees_north', 'Y')
                  coord_bounds(ind) = 'late_bounds'
            end select
         end do

         var_grdz(1) = coord_attributes('NCAT', 'category maximum thickness', 'm', 'undefined')
         var_grdz(2) = coord_attributes('VGRDi', 'vertical ice levels', '1', 'undefined')
         var_grdz(3) = coord_attributes('VGRDs', 'vertical snow levels', '1', 'undefined')
         var_grdz(4) = coord_attributes('VGRDb', 'vertical ice-bio levels', '1', 'undefined')
         var_grdz(5) = coord_attributes('VGRDa', 'vertical snow-ice-bio levels', '1', 'undefined')
         var_grdz(6) = coord_attributes('NFSD', 'category floe size (center)', 'm', 'undefined')

         !-----------------------------------------------------------------
         ! define information for optional time-invariant variables
         !-----------------------------------------------------------------

         var_grd(n_tmask)%req = coord_attributes('tmask', &
                     'mask of T grid cells, 0 = land, 1 = ocean', '1', 'undefined')
         var_grd(n_tmask)%coordinates = 'TLON TLAT'
         var_grd(n_umask)%req = coord_attributes('umask', &
                     'mask of U grid cells, 0 = land, 1 = ocean', '1', 'undefined')
         var_grd(n_umask)%coordinates = 'ULON ULAT'
         var_grd(n_nmask)%req = coord_attributes('nmask', &
                     'mask of N grid cells, 0 = land, 1 = ocean', '1', 'undefined')
         var_grd(n_nmask)%coordinates = 'NLON NLAT'
         var_grd(n_emask)%req = coord_attributes('emask', &
                     'mask of E grid cells, 0 = land, 1 = ocean', '1', 'undefined')
         var_grd(n_emask)%coordinates = 'ELON ELAT'

         var_grd(n_tarea)%req = coord_attributes('tarea', &
                     'area of T grid cells', 'm^2', 'undefined')
         var_grd(n_tarea)%coordinates = 'TLON TLAT'
         var_grd(n_uarea)%req = coord_attributes('uarea', &
                     'area of U grid cells', 'm^2', 'undefined')
         var_grd(n_uarea)%coordinates = 'ULON ULAT'
         var_grd(n_narea)%req = coord_attributes('narea', &
                     'area of N grid cells', 'm^2', 'undefined')
         var_grd(n_narea)%coordinates = 'NLON NLAT'
         var_grd(n_earea)%req = coord_attributes('earea', &
                     'area of E grid cells', 'm^2', 'undefined')
         var_grd(n_earea)%coordinates = 'ELON ELAT'

         var_grd(n_blkmask)%req = coord_attributes('blkmask', &
                     'block id of T grid cells, mytask + iblk/100', '1', 'undefined')
         var_grd(n_blkmask)%coordinates = 'TLON TLAT'

         var_grd(n_dxt)%req = coord_attributes('dxt', &
                     'T cell width through middle', 'm', 'undefined')
         var_grd(n_dxt)%coordinates = 'TLON TLAT'
         var_grd(n_dyt)%req = coord_attributes('dyt', &
                     'T cell height through middle', 'm', 'undefined')
         var_grd(n_dyt)%coordinates = 'TLON TLAT'
         var_grd(n_dxu)%req = coord_attributes('dxu', &
                     'U cell width through middle', 'm', 'undefined')
         var_grd(n_dxu)%coordinates = 'ULON ULAT'
         var_grd(n_dyu)%req = coord_attributes('dyu', &
                     'U cell height through middle', 'm', 'undefined')
         var_grd(n_dyu)%coordinates = 'ULON ULAT'
         var_grd(n_dxn)%req = coord_attributes('dxn', &
                     'N cell width through middle', 'm', 'undefined')
         var_grd(n_dxn)%coordinates = 'NLON NLAT'
         var_grd(n_dyn)%req = coord_attributes('dyn', &
                     'N cell height through middle', 'm', 'undefined')
         var_grd(n_dyn)%coordinates = 'NLON NLAT'
         var_grd(n_dxe)%req = coord_attributes('dxe', &
                     'E cell width through middle', 'm', 'undefined')
         var_grd(n_dxe)%coordinates = 'ELON ELAT'
         var_grd(n_dye)%req = coord_attributes('dye', &
                     'E cell height through middle', 'm', 'undefined')
         var_grd(n_dye)%coordinates = 'ELON ELAT'

         var_grd(n_HTN)%req = coord_attributes('HTN', &
                     'T cell width on North side','m', 'undefined')
         var_grd(n_HTN)%coordinates = 'TLON TLAT'
         var_grd(n_HTE)%req = coord_attributes('HTE', &
                     'T cell width on East side', 'm', 'undefined')
         var_grd(n_HTE)%coordinates = 'TLON TLAT'
         var_grd(n_ANGLE)%req = coord_attributes('ANGLE', &
                     'angle grid makes with latitude line on U grid', &
                     'radians', 'undefined')
         var_grd(n_ANGLE)%coordinates = 'ULON ULAT'
         var_grd(n_ANGLET)%req = coord_attributes('ANGLET', &
                     'angle grid makes with latitude line on T grid', &
                     'radians', 'undefined')
         var_grd(n_ANGLET)%coordinates = 'TLON TLAT'

         ! bounds fields are required for CF compliance
         ! dimensions (nx,ny,nverts)
         ! bounds inherit attributes
         var_nverts(n_lont_bnds) = coord_attributes('lont_bounds','und','und','und')
         var_nverts(n_latt_bnds) = coord_attributes('latt_bounds','und','und','und')
         var_nverts(n_lonu_bnds) = coord_attributes('lonu_bounds','und','und','und')
         var_nverts(n_latu_bnds) = coord_attributes('latu_bounds','und','und','und')
         var_nverts(n_lonn_bnds) = coord_attributes('lonn_bounds','und','und','und')
         var_nverts(n_latn_bnds) = coord_attributes('latn_bounds','und','und','und')
         var_nverts(n_lone_bnds) = coord_attributes('lone_bounds','und','und','und')
         var_nverts(n_late_bnds) = coord_attributes('late_bounds','und','und','und')

         !-----------------------------------------------------------------
         ! define attributes for time-invariant variables
         !-----------------------------------------------------------------

         dimid(1) = imtid
         dimid(2) = jmtid
         dimid(3) = timid

         do i = 1, ncoord
            if(icoord(i) .or. histfreq(ns)=='g') then
               call ice_hist_coord_def(ncid, var_coord(i), lprecision, dimid(1:2), varid)
               call ice_write_hist_fill(ncid,varid,var_coord(i)%short_name,history_precision)
               if (var_coord(i)%short_name == 'ULAT') then
                  status = nf90_put_att(ncid,varid,'comment', &
                       'Latitude of NE corner of T grid cell')
                  call ice_check_nc(status, subname// ' ERROR: defining comment for '//var_coord(i)%short_name, &
                                    file=__FILE__, line=__LINE__)
               endif
               if (f_bounds .or. histfreq(ns)=='g') then
                  status = nf90_put_att(ncid, varid, 'bounds', coord_bounds(i))
                  call ice_check_nc(status, subname// ' ERROR: defining bounds for '//var_coord(i)%short_name, &
                                    file=__FILE__, line=__LINE__)
               endif
            else
               extvars = trim(extvars)//' '//trim(var_coord(i)%short_name)
            endif
         enddo

         ! Extra dimensions (NCAT, NZILYR, NZSLYR, NZBLYR, NZALYR, NFSD)
         dimidex(1)=cmtid
         dimidex(2)=kmtidi
         dimidex(3)=kmtids
         dimidex(4)=kmtidb
         dimidex(5)=kmtida
         dimidex(6)=fmtid

         do i = 1, nvar_grdz
            if (igrdz(i) .or. histfreq(ns)=='g') then
               call ice_hist_coord_def(ncid, var_grdz(i), lprecision, dimidex(i:i), varid)
            else
               extvars = trim(extvars)//' '//trim(var_grdz(i)%short_name)
            endif
         enddo

         do i = 1, nvar_grd
            if (igrd(i) .or. histfreq(ns)=='g') then
               call ice_hist_coord_def(ncid, var_grd(i)%req, lprecision, dimid(1:2), varid)
               status = nf90_put_att(ncid, varid, 'coordinates', var_grd(i)%coordinates)
               call ice_check_nc(status, subname// ' ERROR: defining coordinates for '//var_grd(i)%req%short_name, &
                                 file=__FILE__, line=__LINE__)
               call ice_write_hist_fill(ncid,varid,var_grd(i)%req%short_name,history_precision)
            else
               extvars = trim(extvars)//' '//trim(var_grd(i)%req%short_name)
            endif
         enddo

         ! bounds fields with dimensions (nverts,nx,ny)
         ! bounds inherits attributes
         dimid_nverts(1) = nvertexid
         dimid_nverts(2) = imtid
         dimid_nverts(3) = jmtid
         do i = 1, nvar_verts
            if (f_bounds .or. histfreq(ns)=='g') then
               call ice_hist_coord_def(ncid, var_nverts(i), lprecision, dimid_nverts, varid)
            endif
         enddo

         !-----------------------------------------------------------------
         ! define attributes for time-variant variables
         !-----------------------------------------------------------------

         do n=1,num_avail_hist_fields_2D
            if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
               call ice_hist_field_def(ncid, avail_hist_fields(n),lprecision, dimid,ns)
            endif
         enddo  ! num_avail_hist_fields_2D

         dimidz(1) = imtid
         dimidz(2) = jmtid
         dimidz(3) = cmtid
         dimidz(4) = timid

         do n = n2D + 1, n3Dccum
            if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
               call ice_hist_field_def(ncid, avail_hist_fields(n),lprecision, dimidz,ns)
            endif
         enddo  ! num_avail_hist_fields_3Dc

         dimidz(1) = imtid
         dimidz(2) = jmtid
         dimidz(3) = kmtidi
         dimidz(4) = timid

         do n = n3Dccum + 1, n3Dzcum
            if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
               call ice_hist_field_def(ncid, avail_hist_fields(n),lprecision, dimidz,ns)
            endif
         enddo  ! num_avail_hist_fields_3Dz

         dimidz(1) = imtid
         dimidz(2) = jmtid
         dimidz(3) = kmtidb
         dimidz(4) = timid

         do n = n3Dzcum + 1, n3Dbcum
            if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
               call ice_hist_field_def(ncid, avail_hist_fields(n),lprecision, dimidz,ns)
            endif
         enddo  ! num_avail_hist_fields_3Db

         dimidz(1) = imtid
         dimidz(2) = jmtid
         dimidz(3) = kmtida
         dimidz(4) = timid

         do n = n3Dbcum + 1, n3Dacum
            if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
               call ice_hist_field_def(ncid, avail_hist_fields(n),lprecision, dimidz,ns)
            endif
         enddo  ! num_avail_hist_fields_3Da

         dimidz(1) = imtid
         dimidz(2) = jmtid
         dimidz(3) = fmtid
         dimidz(4) = timid

         do n = n3Dacum + 1, n3Dfcum
            if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
               call ice_hist_field_def(ncid, avail_hist_fields(n),lprecision, dimidz,ns)
            endif
         enddo  ! num_avail_hist_fields_3Df

         dimidcz(1) = imtid
         dimidcz(2) = jmtid
         dimidcz(3) = kmtidi
         dimidcz(4) = cmtid
         dimidcz(5) = timid

         do n = n3Dfcum + 1, n4Dicum
            if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
               call ice_hist_field_def(ncid, avail_hist_fields(n),lprecision, dimidcz,ns)
            endif
         enddo  ! num_avail_hist_fields_4Di

         dimidcz(1) = imtid
         dimidcz(2) = jmtid
         dimidcz(3) = kmtids
         dimidcz(4) = cmtid
         dimidcz(5) = timid

         do n = n4Dicum + 1, n4Dscum
            if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
               call ice_hist_field_def(ncid, avail_hist_fields(n),lprecision, dimidcz,ns)
            endif
         enddo  ! num_avail_hist_fields_4Ds

         dimidcz(1) = imtid
         dimidcz(2) = jmtid
         dimidcz(3) = fmtid
         dimidcz(4) = cmtid
         dimidcz(5) = timid

         do n = n4Dscum + 1, n4Dfcum
            if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
               call ice_hist_field_def(ncid, avail_hist_fields(n),lprecision, &
                  ! dimidcz, ns)
                  dimidcz(1:4),ns) ! ferret
             endif
         enddo  ! num_avail_hist_fields_4Df

         !-----------------------------------------------------------------
         ! global attributes
         !-----------------------------------------------------------------
         ! ... the user should change these to something useful ...
         !-----------------------------------------------------------------
#ifdef CESMCOUPLED
         status = nf90_put_att(ncid,nf90_global,'title',runid)
         call ice_check_nc(status, subname// ' ERROR: in global attribute title', &
                           file=__FILE__, line=__LINE__)
#else
         title  = 'sea ice model output for CICE'
         status = nf90_put_att(ncid,nf90_global,'title',title)
         call ice_check_nc(status, subname// ' ERROR: in global attribute title', &
                           file=__FILE__, line=__LINE__)
#endif
         title = 'Diagnostic and Prognostic Variables'
         status = nf90_put_att(ncid,nf90_global,'contents',title)
         call ice_check_nc(status, subname// ' ERROR: global attribute contents', &
                           file=__FILE__, line=__LINE__)

         write(title,'(2a)') 'CICE Sea Ice Model, ', trim(version_name)
         status = nf90_put_att(ncid,nf90_global,'source',title)
         call ice_check_nc(status, subname// ' ERROR: global attribute source', &
                           file=__FILE__, line=__LINE__)

         if (use_leap_years) then
           write(title,'(a,i3,a)') 'This year has ',dayyr,' days'
         else
           write(title,'(a,i3,a)') 'All years have exactly ',dayyr,' days'
         endif
         status = nf90_put_att(ncid,nf90_global,'comment',title)
         call ice_check_nc(status, subname// ' ERROR: global attribute comment', &
                           file=__FILE__, line=__LINE__)

         write(title,'(a,i8.8)') 'File written on model date ',idate
         status = nf90_put_att(ncid,nf90_global,'comment2',title)
         call ice_check_nc(status, subname// ' ERROR: global attribute date1', &
                           file=__FILE__, line=__LINE__)

         write(title,'(a,i6)') 'seconds elapsed into model date: ',msec
         status = nf90_put_att(ncid,nf90_global,'comment3',title)
         call ice_check_nc(status, subname// ' ERROR: global attribute date2', &
                           file=__FILE__, line=__LINE__)

         select case (histfreq(ns))
            case ("y", "Y")
               write(time_period_freq,'(a,i0)') 'year_',histfreq_n(ns)
            case ("m", "M")
               write(time_period_freq,'(a,i0)') 'month_',histfreq_n(ns)
            case ("d", "D")
               write(time_period_freq,'(a,i0)') 'day_',histfreq_n(ns)
            case ("h", "H")
               write(time_period_freq,'(a,i0)') 'hour_',histfreq_n(ns)
            case ("1")
               write(time_period_freq,'(a,i0)') 'step_',histfreq_n(ns)
         end select

         if (.not.write_ic .and. trim(time_period_freq) /= 'none') then
            status = nf90_put_att(ncid,nf90_global,'time_period_freq',trim(time_period_freq))
            call ice_check_nc(status, subname// ' ERROR: global attribute time_period_freq', &
                              file=__FILE__, line=__LINE__)
         endif

         if (hist_avg(ns)) then
            status = nf90_put_att(ncid,nf90_global,'time_axis_position',trim(hist_time_axis))
            call ice_check_nc(status, subname// ' ERROR: global attribute time axis position', &
                              file=__FILE__, line=__LINE__)
         endif

         title = 'CF-1.8'
         status =  nf90_put_att(ncid,nf90_global,'Conventions',title)
         call ice_check_nc(status, subname// ' ERROR: in global attribute conventions', &
                           file=__FILE__, line=__LINE__)

         status =  nf90_put_att(ncid,nf90_global,'external_variables',trim(extvars))
         call ice_check_nc(status, subname// ' ERROR: in global attribute external_variables', &
                           file=__FILE__, line=__LINE__)

         call date_and_time(date=current_date, time=current_time)
         write(start_time,1000) current_date(1:4), current_date(5:6), &
                                current_date(7:8), current_time(1:2), &
                                current_time(3:4), current_time(5:8)
1000     format('This dataset was created on ', &
                 a,'-',a,'-',a,' at ',a,':',a,':',a)

         status = nf90_put_att(ncid,nf90_global,'history',start_time)
         call ice_check_nc(status, subname// ' ERROR: global attribute history', &
                           file=__FILE__, line=__LINE__)

         write(start_time,1001) current_date(1:4), current_date(5:6), &
                                current_date(7:8), current_time(1:2), &
                                current_time(3:4), current_time(5:8)
1001     format(a,'-',a,'-',a,' ',a,':',a,':',a)

         status = nf90_put_att(ncid,nf90_global,'date_created',start_time)
         call ice_check_nc(status, subname// ' ERROR: global attribute date_created', &
                           file=__FILE__, line=__LINE__)

         status = nf90_put_att(ncid,nf90_global,'io_flavor','io_netcdf '//trim(history_format))
         call ice_check_nc(status, subname// ' ERROR: global attribute io_flavor', &
                           file=__FILE__, line=__LINE__)

         !-----------------------------------------------------------------
         ! end define mode
         !-----------------------------------------------------------------

         status = nf90_enddef(ncid)
         call ice_check_nc(status, subname// ' ERROR: in nf90_enddef', &
                           file=__FILE__, line=__LINE__)

         !-----------------------------------------------------------------
         ! write time variable and time bounds info
         !-----------------------------------------------------------------

         ! do not write time axis on grid output file
         if (histfreq(ns)/='g') then

            ltime2 = timesecs/secday ! hist_time_axis = 'end' (default)

            ! Some coupled models require the time axis "stamp" to be in the middle
            ! or even beginning of averaging interval.
            if (hist_avg(ns)) then
               if (trim(hist_time_axis) == "begin" ) ltime2 = time_beg(ns)
               if (trim(hist_time_axis) == "middle") ltime2 = p5*(time_beg(ns)+time_end(ns))
            endif

            status = nf90_inq_varid(ncid,'time',varid)
            call ice_check_nc(status, subname// ' ERROR: getting time varid', &
                              file=__FILE__, line=__LINE__)
            status = nf90_put_var(ncid,varid,ltime2)
            call ice_check_nc(status, subname// ' ERROR: writing time variable', &
                              file=__FILE__, line=__LINE__)

            if (hist_avg(ns) .and. .not. write_ic) then
               status = nf90_inq_varid(ncid,'time_bounds',varid)
               call ice_check_nc(status, subname// ' ERROR: getting time_bounds id', &
                                 file=__FILE__, line=__LINE__)
               status = nf90_put_var(ncid,varid,time_beg(ns),start=(/1/))
               call ice_check_nc(status, subname// ' ERROR: writing time_beg', &
                                 file=__FILE__, line=__LINE__)
               status = nf90_put_var(ncid,varid,time_end(ns),start=(/2/))
               call ice_check_nc(status, subname// ' ERROR: writing time_end', &
                                 file=__FILE__, line=__LINE__)
            endif

         endif  ! histfreq(ns)/='g'

      endif                     ! master_task

      if (my_task==master_task) then
         allocate(work_g1(nx_global,ny_global))
      else
         allocate(work_g1(1,1))
      endif

      work_g1(:,:) = c0

      !-----------------------------------------------------------------
      ! write coordinate variables
      !-----------------------------------------------------------------

      do i = 1,ncoord
         if(icoord(i) .or. histfreq(ns)=='g') then
            call broadcast_scalar(var_coord(i)%short_name,master_task)
            SELECT CASE (var_coord(i)%short_name)
               CASE ('TLON')
                  ! Convert T grid longitude from -180 -> 180 to 0 to 360
                  work1 = TLON*rad_to_deg + c360
                  where (work1 > c360) work1 = work1 - c360
                  where (work1 < c0 )  work1 = work1 + c360
                  call gather_global(work_g1,work1,master_task,distrb_info)
               CASE ('TLAT')
                  work1 = TLAT*rad_to_deg
                  call gather_global(work_g1,work1,master_task,distrb_info)
               CASE ('ULON')
                  work1 = ULON*rad_to_deg
                  call gather_global(work_g1,work1,master_task,distrb_info)
               CASE ('ULAT')
                  work1 = ULAT*rad_to_deg
                  call gather_global(work_g1,work1,master_task,distrb_info)
               CASE ('NLON')
                  work1 = NLON*rad_to_deg
                  call gather_global(work_g1,work1,master_task,distrb_info)
               CASE ('NLAT')
                  work1 = NLAT*rad_to_deg
                  call gather_global(work_g1,work1,master_task,distrb_info)
               CASE ('ELON')
                  work1 = ELON*rad_to_deg
                  call gather_global(work_g1,work1,master_task,distrb_info)
               CASE ('ELAT')
                  work1 = ELAT*rad_to_deg
                  call gather_global(work_g1,work1,master_task,distrb_info)
            END SELECT

            if (my_task == master_task) then
               status = nf90_inq_varid(ncid, var_coord(i)%short_name, varid)
               call ice_check_nc(status, subname// ' ERROR: getting varid for '//var_coord(i)%short_name, &
                                 file=__FILE__, line=__LINE__)
               status = nf90_put_var(ncid,varid,work_g1)
               call ice_check_nc(status, subname// ' ERROR: writing'//var_coord(i)%short_name, &
                                 file=__FILE__, line=__LINE__)
            endif
         endif
      enddo

      ! Extra dimensions (NCAT, NFSD, VGRD*)

      do i = 1, nvar_grdz
         if (igrdz(i) .or. histfreq(ns)=='g') then
            call broadcast_scalar(var_grdz(i)%short_name,master_task)
            if (my_task == master_task) then
               status = nf90_inq_varid(ncid, var_grdz(i)%short_name, varid)
               call ice_check_nc(status, subname// ' ERROR: getting varid for '//var_grdz(i)%short_name, &
                                 file=__FILE__, line=__LINE__)
               SELECT CASE (var_grdz(i)%short_name)
                  CASE ('NCAT')
                     status = nf90_put_var(ncid,varid,hin_max(1:ncat_hist))
                  CASE ('NFSD')
                     status = nf90_put_var(ncid,varid,floe_rad_c(1:nfsd_hist))
                  CASE ('VGRDi') ! index - needed for Met Office analysis code
                     status = nf90_put_var(ncid,varid,(/(k, k=1,nzilyr)/))
                  CASE ('VGRDs') ! index - needed for Met Office analysis code
                     status = nf90_put_var(ncid,varid,(/(k, k=1,nzslyr)/))
                  CASE ('VGRDb')
                     status = nf90_put_var(ncid,varid,(/(k, k=1,nzblyr)/))
                  CASE ('VGRDa')
                     status = nf90_put_var(ncid,varid,(/(k, k=1,nzalyr)/))
               END SELECT
               call ice_check_nc(status, subname// ' ERROR: put var '//var_grdz(i)%short_name, &
                                 file=__FILE__, line=__LINE__)
            endif
         endif
      enddo

      !-----------------------------------------------------------------
      ! write grid masks, area and rotation angle
      !-----------------------------------------------------------------

      do i = 1, nvar_grd
         if (igrd(i) .or. histfreq(ns)=='g') then
            call broadcast_scalar(var_grd(i)%req%short_name,master_task)
            SELECT CASE (var_grd(i)%req%short_name)
               CASE ('tmask')
                  call gather_global(work_g1,    hm, master_task, distrb_info)
               CASE ('umask')
                  call gather_global(work_g1,   uvm, master_task, distrb_info)
               CASE ('nmask')
                  call gather_global(work_g1,   npm, master_task, distrb_info)
               CASE ('emask')
                  call gather_global(work_g1,   epm, master_task, distrb_info)
               CASE ('tarea')
                  call gather_global(work_g1, tarea, master_task, distrb_info)
               CASE ('uarea')
                  call gather_global(work_g1, uarea, master_task, distrb_info)
               CASE ('narea')
                  call gather_global(work_g1, narea, master_task, distrb_info)
               CASE ('earea')
                  call gather_global(work_g1, earea, master_task, distrb_info)
               CASE ('blkmask')
                  call gather_global(work_g1,    bm, master_task, distrb_info)
               CASE ('dxu')
                  call gather_global(work_g1,   dxU, master_task, distrb_info)
               CASE ('dyu')
                  call gather_global(work_g1,   dyU, master_task, distrb_info)
               CASE ('dxt')
                  call gather_global(work_g1,   dxT, master_task, distrb_info)
               CASE ('dyt')
                  call gather_global(work_g1,   dyT, master_task, distrb_info)
               CASE ('dxn')
                  call gather_global(work_g1,   dxN, master_task, distrb_info)
               CASE ('dyn')
                  call gather_global(work_g1,   dyN, master_task, distrb_info)
               CASE ('dxe')
                  call gather_global(work_g1,   dxE, master_task, distrb_info)
               CASE ('dye')
                  call gather_global(work_g1,   dyE, master_task, distrb_info)
               CASE ('HTN')
                  call gather_global(work_g1,   HTN, master_task, distrb_info)
               CASE ('HTE')
                  call gather_global(work_g1,   HTE, master_task, distrb_info)
               CASE ('ANGLE')
                  call gather_global(work_g1, ANGLE, master_task, distrb_info)
               CASE ('ANGLET')
                  call gather_global(work_g1, ANGLET,master_task, distrb_info)
            END SELECT

            if (my_task == master_task) then
               status = nf90_inq_varid(ncid, var_grd(i)%req%short_name, varid)
               call ice_check_nc(status, subname// ' ERROR: getting varid for '//var_grd(i)%req%short_name, &
                                 file=__FILE__, line=__LINE__)
               status = nf90_put_var(ncid,varid,work_g1)
               call ice_check_nc(status, subname// ' ERROR: writing variable '//var_grd(i)%req%short_name, &
                                 file=__FILE__, line=__LINE__)
            endif
         endif
      enddo

      !----------------------------------------------------------------
      ! Write coordinates of grid box vertices
      !----------------------------------------------------------------

      if (f_bounds .or. histfreq(ns)=='g') then
         if (my_task==master_task) then
            allocate(work1_3(nverts,nx_global,ny_global))
         else
            allocate(work1_3(1,1,1))   ! to save memory
         endif

         work1_3(:,:,:) = c0
         work1  (:,:,:) = c0

         do i = 1, nvar_verts
            call broadcast_scalar(var_nverts(i)%short_name,master_task)
            SELECT CASE (var_nverts(i)%short_name)
               CASE ('lont_bounds')
                  do ivertex = 1, nverts
                     work1(:,:,:) = lont_bounds(ivertex,:,:,:)
                     call gather_global(work_g1, work1, master_task, distrb_info)
                     if (my_task == master_task) work1_3(ivertex,:,:) = work_g1(:,:)
                  enddo
               CASE ('latt_bounds')
                  do ivertex = 1, nverts
                     work1(:,:,:) = latt_bounds(ivertex,:,:,:)
                     call gather_global(work_g1, work1, master_task, distrb_info)
                     if (my_task == master_task) work1_3(ivertex,:,:) = work_g1(:,:)
                  enddo
               CASE ('lonu_bounds')
                  do ivertex = 1, nverts
                     work1(:,:,:) = lonu_bounds(ivertex,:,:,:)
                     call gather_global(work_g1, work1, master_task, distrb_info)
                     if (my_task == master_task) work1_3(ivertex,:,:) = work_g1(:,:)
                  enddo
               CASE ('latu_bounds')
                  do ivertex = 1, nverts
                     work1(:,:,:) = latu_bounds(ivertex,:,:,:)
                     call gather_global(work_g1, work1, master_task, distrb_info)
                     if (my_task == master_task) work1_3(ivertex,:,:) = work_g1(:,:)
                  enddo
               CASE ('lonn_bounds')
                  do ivertex = 1, nverts
                     work1(:,:,:) = lonn_bounds(ivertex,:,:,:)
                     call gather_global(work_g1, work1, master_task, distrb_info)
                     if (my_task == master_task) work1_3(ivertex,:,:) = work_g1(:,:)
                  enddo
               CASE ('latn_bounds')
                  do ivertex = 1, nverts
                     work1(:,:,:) = latn_bounds(ivertex,:,:,:)
                     call gather_global(work_g1, work1, master_task, distrb_info)
                     if (my_task == master_task) work1_3(ivertex,:,:) = work_g1(:,:)
                  enddo
               CASE ('lone_bounds')
                  do ivertex = 1, nverts
                     work1(:,:,:) = lone_bounds(ivertex,:,:,:)
                     call gather_global(work_g1, work1, master_task, distrb_info)
                     if (my_task == master_task) work1_3(ivertex,:,:) = work_g1(:,:)
                  enddo
               CASE ('late_bounds')
                  do ivertex = 1, nverts
                     work1(:,:,:) = late_bounds(ivertex,:,:,:)
                     call gather_global(work_g1, work1, master_task, distrb_info)
                     if (my_task == master_task) work1_3(ivertex,:,:) = work_g1(:,:)
                  enddo
            END SELECT

            if (my_task == master_task) then
               status = nf90_inq_varid(ncid, var_nverts(i)%short_name, varid)
               call ice_check_nc(status, subname// ' ERROR: getting varid for '//var_nverts(i)%short_name, &
                                 file=__FILE__, line=__LINE__)
               status = nf90_put_var(ncid,varid,work1_3)
               call ice_check_nc(status, subname// ' ERROR: writing variable '//var_nverts(i)%short_name, &
                                 file=__FILE__, line=__LINE__)
            endif
         enddo
         deallocate(work1_3)
      endif

      !-----------------------------------------------------------------
      ! write variable data
      !-----------------------------------------------------------------

      work_g1(:,:) = c0

      do n=1,num_avail_hist_fields_2D
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            call gather_global(work_g1, a2D(:,:,n,:), &
                               master_task, distrb_info)
            if (my_task == master_task) then
               status  = nf90_inq_varid(ncid,avail_hist_fields(n)%vname,varid)
               call ice_check_nc(status, subname// ' ERROR: getting varid for '//avail_hist_fields(n)%vname, &
                                 file=__FILE__, line=__LINE__)
               status  = nf90_put_var(ncid,varid,work_g1, &
                                      count=(/nx_global,ny_global/))
               call ice_check_nc(status, subname// ' ERROR: writing variable '//avail_hist_fields(n)%vname, &
                                 file=__FILE__, line=__LINE__)
            endif

         endif
      enddo ! num_avail_hist_fields_2D

      work_g1(:,:) = c0

      do n = n2D + 1, n3Dccum
         nn = n - n2D
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            if (my_task == master_task) then
               status  = nf90_inq_varid(ncid,avail_hist_fields(n)%vname,varid)
               call ice_check_nc(status, subname// ' ERROR: getting varid for '//avail_hist_fields(n)%vname, &
                                 file=__FILE__, line=__LINE__)
            endif
            do k = 1, ncat_hist
               call gather_global(work_g1, a3Dc(:,:,k,nn,:), &
                                  master_task, distrb_info)

               if (my_task == master_task) then
                  status  = nf90_inq_varid(ncid,avail_hist_fields(n)%vname,varid)
                  call ice_check_nc(status, subname// ' ERROR: getting varid for '//avail_hist_fields(n)%vname, &
                                    file=__FILE__, line=__LINE__)
                  status  = nf90_put_var(ncid,varid,work_g1, &
                                         start=(/        1,        1,k/), &
                                         count=(/nx_global,ny_global,1/))
                  call ice_check_nc(status, subname// ' ERROR: writing variable '//avail_hist_fields(n)%vname, &
                                    file=__FILE__, line=__LINE__)
               endif
            enddo ! k
         endif
      enddo ! num_avail_hist_fields_3Dc

      work_g1(:,:) = c0

      do n = n3Dccum+1, n3Dzcum
         nn = n - n3Dccum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            if (my_task == master_task) then
               status  = nf90_inq_varid(ncid,avail_hist_fields(n)%vname,varid)
               call ice_check_nc(status, subname// ' ERROR: getting varid for '//avail_hist_fields(n)%vname, &
                                 file=__FILE__, line=__LINE__)
            endif
            do k = 1, nzilyr
               call gather_global(work_g1, a3Dz(:,:,k,nn,:), &
                                  master_task, distrb_info)

               if (my_task == master_task) then
                  status  = nf90_put_var(ncid,varid,work_g1, &
                                         start=(/        1,        1,k/), &
                                         count=(/nx_global,ny_global,1/))
                  call ice_check_nc(status, subname// ' ERROR: writing variable '//avail_hist_fields(n)%vname, &
                                    file=__FILE__, line=__LINE__)
               endif
            enddo ! k
         endif
      enddo ! num_avail_hist_fields_3Dz

      work_g1(:,:) = c0

      do n = n3Dzcum+1, n3Dbcum
         nn = n - n3Dzcum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            if (my_task == master_task) then
               status  = nf90_inq_varid(ncid,avail_hist_fields(n)%vname,varid)
               call ice_check_nc(status, subname// ' ERROR: getting varid for '//avail_hist_fields(n)%vname, &
                                 file=__FILE__, line=__LINE__)
            endif
            do k = 1, nzblyr
               call gather_global(work_g1, a3Db(:,:,k,nn,:), &
                                  master_task, distrb_info)

               if (my_task == master_task) then
                  status  = nf90_put_var(ncid,varid,work_g1, &
                                         start=(/        1,        1,k/), &
                                         count=(/nx_global,ny_global,1/))
                  call ice_check_nc(status, subname// ' ERROR: writing variable '//avail_hist_fields(n)%vname, &
                                    file=__FILE__, line=__LINE__)
               endif
            enddo ! k
         endif
      enddo ! num_avail_hist_fields_3Db

      work_g1(:,:) = c0

      do n = n3Dbcum+1, n3Dacum
         nn = n - n3Dbcum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            if (my_task == master_task) then
               status  = nf90_inq_varid(ncid,avail_hist_fields(n)%vname,varid)
               call ice_check_nc(status, subname// ' ERROR: getting varid for '//avail_hist_fields(n)%vname, &
                                 file=__FILE__, line=__LINE__)
            endif
            do k = 1, nzalyr
               call gather_global(work_g1, a3Da(:,:,k,nn,:), &
                                  master_task, distrb_info)

               if (my_task == master_task) then
                  status  = nf90_put_var(ncid,varid,work_g1, &
                                         start=(/        1,        1,k/), &
                                         count=(/nx_global,ny_global,1/))
                  call ice_check_nc(status, subname// ' ERROR: writing variable '//avail_hist_fields(n)%vname, &
                                    file=__FILE__, line=__LINE__)
               endif
            enddo ! k
         endif
      enddo ! num_avail_hist_fields_3Da

      work_g1(:,:) = c0

      do n = n3Dacum+1, n3Dfcum
         nn = n - n3Dacum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            if (my_task == master_task) then
               status  = nf90_inq_varid(ncid,avail_hist_fields(n)%vname,varid)
               call ice_check_nc(status, subname// ' ERROR: getting varid for '//avail_hist_fields(n)%vname, &
                                 file=__FILE__, line=__LINE__)
            endif
            do k = 1, nfsd_hist
               call gather_global(work_g1, a3Df(:,:,k,nn,:), &
                                  master_task, distrb_info)
               if (my_task == master_task) then
                  status  = nf90_put_var(ncid,varid,work_g1, &
                                         start=(/        1,        1,k/), &
                                         count=(/nx_global,ny_global,1/))
                  call ice_check_nc(status, subname// ' ERROR: writing variable '//avail_hist_fields(n)%vname, &
                                    file=__FILE__, line=__LINE__)
               endif
            enddo ! k
         endif
      enddo ! num_avail_hist_fields_3Df

      work_g1(:,:) = c0

      do n = n3Dfcum+1, n4Dicum
         nn = n - n3Dfcum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            if (my_task == master_task) then
               status  = nf90_inq_varid(ncid,avail_hist_fields(n)%vname,varid)
               call ice_check_nc(status, subname// ' ERROR: getting varid for '//avail_hist_fields(n)%vname, &
                                 file=__FILE__, line=__LINE__)
            endif
            do ic = 1, ncat_hist
            do k = 1, nzilyr
               call gather_global(work_g1, a4Di(:,:,k,ic,nn,:), &
                                  master_task, distrb_info)
               if (my_task == master_task) then
                  status  = nf90_put_var(ncid,varid,work_g1, &
                                         start=(/        1,        1,k,ic/), &
                                         count=(/nx_global,ny_global,1, 1/))
                  call ice_check_nc(status, subname// ' ERROR: writing variable '//avail_hist_fields(n)%vname, &
                                    file=__FILE__, line=__LINE__)
               endif
            enddo ! k
            enddo ! ic
         endif
      enddo ! num_avail_hist_fields_4Di

      work_g1(:,:) = c0

      do n = n4Dicum+1, n4Dscum
         nn = n - n4Dicum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            if (my_task == master_task) then
               status  = nf90_inq_varid(ncid,avail_hist_fields(n)%vname,varid)
               call ice_check_nc(status, subname// ' ERROR: getting varid for '//avail_hist_fields(n)%vname, &
                                 file=__FILE__, line=__LINE__)
            endif
            do ic = 1, ncat_hist
            do k = 1, nzslyr
               call gather_global(work_g1, a4Ds(:,:,k,ic,nn,:), &
                                  master_task, distrb_info)
               if (my_task == master_task) then
                  status  = nf90_put_var(ncid,varid,work_g1, &
                                         start=(/        1,        1,k,ic/), &
                                         count=(/nx_global,ny_global,1, 1/))
                  call ice_check_nc(status, subname// ' ERROR: writing variable '//avail_hist_fields(n)%vname, &
                                    file=__FILE__, line=__LINE__)
               endif
            enddo ! k
            enddo ! ic
         endif
      enddo ! num_avail_hist_fields_4Ds

      do n = n4Dscum+1, n4Dfcum
         nn = n - n4Dscum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            if (my_task == master_task) then
               status  = nf90_inq_varid(ncid,avail_hist_fields(n)%vname,varid)
               call ice_check_nc(status, subname// ' ERROR: getting varid for '//avail_hist_fields(n)%vname, &
                                 file=__FILE__, line=__LINE__)
            endif
            do ic = 1, ncat_hist
            do k = 1, nfsd_hist
               call gather_global(work_g1, a4Df(:,:,k,ic,nn,:), &
                                  master_task, distrb_info)
               if (my_task == master_task) then
                  status  = nf90_put_var(ncid,varid,work_g1, &
                                         start=(/        1,        1,k,ic/), &
                                         count=(/nx_global,ny_global,1, 1/))
                  call ice_check_nc(status, subname// ' ERROR: writing variable '//avail_hist_fields(n)%vname, &
                                    file=__FILE__, line=__LINE__)
               endif
            enddo ! k
            enddo ! ic
         endif
      enddo ! num_avail_hist_fields_4Df

      deallocate(work_g1)

      !-----------------------------------------------------------------
      ! close output dataset
      !-----------------------------------------------------------------

      if (my_task == master_task) then
         status = nf90_close(ncid)
         call ice_check_nc(status, subname// ' ERROR: closing netCDF history file', &
                           file=__FILE__, line=__LINE__)
         write(nu_diag,*) ' '
         write(nu_diag,*) 'Finished writing ',trim(ncfile)
      endif

#else
      call abort_ice(subname//' ERROR: USE_NETCDF cpp not defined', file=__FILE__, line=__LINE__)
#endif

      end subroutine ice_write_hist

!=======================================================================
! Defines a (time-dependent) history var in the history file
! variables have short_name, long_name and units, coordiantes and cell_measures attributes,
!  and are compressed and chunked for 'hdf5'

      subroutine ice_hist_field_def(ncid, hfield, lprecision, dimids, ns)

      use ice_history_shared, only: history_deflate, history_chunksize, history_format, ice_hist_field, &
         history_precision, hist_avg
      use ice_calendar, only: histfreq, histfreq_n, write_ic

      integer(kind=int_kind), intent(in) :: ncid, dimids(:), lprecision, ns
      type(ice_hist_field), intent(in) :: hfield

      !local vars
      integer(kind=int_kind) :: chunks(size(dimids)), i, status, varid

      character(len=*), parameter :: subname = '(ice_hist_field_def)'

#ifdef USE_NETCDF
      status = nf90_def_var(ncid, hfield%vname, lprecision, dimids, varid)
      call ice_check_nc(status, subname//' ERROR: defining var '//trim(hfield%vname),file=__FILE__,line=__LINE__)

#ifdef NO_HDF5
      if (history_format=='hdf5') then
          call abort_ice(subname//' ERROR: history_format hdf5 not available ', &
              file=__FILE__, line=__LINE__)
      endif
#else
      if (history_format=='hdf5' .and. size(dimids)>1) then
         if (dimids(1)==imtid .and. dimids(2)==jmtid) then
            chunks(1)=history_chunksize(1)
            chunks(2)=history_chunksize(2)
            do i = 3, size(dimids)
               chunks(i) = 0
            enddo
            status = nf90_def_var_chunking(ncid, varid, NF90_CHUNKED, chunksizes=chunks)
            call ice_check_nc(status, subname//' ERROR chunking var '//trim(hfield%vname), file=__FILE__, line=__LINE__)
         endif
      endif

      if (history_format=='hdf5' .and. history_deflate/=0) then
         status = nf90_def_var_deflate(ncid, varid, shuffle=0, deflate=1, deflate_level=history_deflate)
         call ice_check_nc(status, subname//' ERROR deflating var '//trim(hfield%vname), file=__FILE__, line=__LINE__)
      endif
#endif

      ! add attributes
      status = nf90_put_att(ncid,varid,'units', hfield%vunit)
      call ice_check_nc(status, subname// ' ERROR: defining units for '//hfield%vname, &
                        file=__FILE__, line=__LINE__)

      status = nf90_put_att(ncid,varid, 'long_name', hfield%vdesc)
      call ice_check_nc(status, subname// ' ERROR: defining long_name for '//hfield%vname, &
                        file=__FILE__, line=__LINE__)

      status = nf90_put_att(ncid,varid,'coordinates', hfield%vcoord)
      call ice_check_nc(status, subname// ' ERROR: defining coordinates for '//hfield%vname, &
                        file=__FILE__, line=__LINE__)

      status = nf90_put_att(ncid,varid,'cell_measures', hfield%vcellmeas)
      call ice_check_nc(status, subname// ' ERROR: defining cell measures for '//hfield%vname, &
                        file=__FILE__, line=__LINE__)

      if (hfield%vcomment /= "none") then
         status = nf90_put_att(ncid,varid,'comment', hfield%vcomment)
         call ice_check_nc(status, subname// ' ERROR: defining comment for '//hfield%vname, &
                           file=__FILE__, line=__LINE__)
      endif

      call ice_write_hist_fill(ncid,varid,hfield%vname,history_precision)

      ! Add cell_methods attribute to variables if averaged
      if (hist_avg(ns) .and. .not. write_ic) then
         if    (TRIM(hfield%vname(1:4))/='sig1' &
           .and.TRIM(hfield%vname(1:4))/='sig2' &
           .and.TRIM(hfield%vname(1:9))/='sistreave' &
           .and.TRIM(hfield%vname(1:9))/='sistremax' &
           .and.TRIM(hfield%vname(1:4))/='sigP') then
             status = nf90_put_att(ncid,varid,'cell_methods','time: mean')
             call ice_check_nc(status, subname// ' ERROR: defining cell methods for '//hfield%vname, &
                               file=__FILE__, line=__LINE__)
         endif
      endif

      if ((histfreq(ns) == '1' .and. histfreq_n(ns) == 1) &
          .or..not. hist_avg(ns)                          &
          .or. write_ic                                   &
          .or.TRIM(hfield%vname(1:4))=='divu' &
          .or.TRIM(hfield%vname(1:5))=='shear' &
          .or.TRIM(hfield%vname(1:4))=='vort' &
          .or.TRIM(hfield%vname(1:4))=='sig1' &
          .or.TRIM(hfield%vname(1:4))=='sig2' &
          .or.TRIM(hfield%vname(1:4))=='sigP' &
          .or.TRIM(hfield%vname(1:5))=='trsig' &
          .or.TRIM(hfield%vname(1:9))=='sistreave' &
          .or.TRIM(hfield%vname(1:9))=='sistremax' &
          .or.TRIM(hfield%vname(1:9))=='mlt_onset' &
          .or.TRIM(hfield%vname(1:9))=='frz_onset' &
          .or.TRIM(hfield%vname(1:6))=='hisnap' &
          .or.TRIM(hfield%vname(1:6))=='aisnap') then
         status = nf90_put_att(ncid,varid,'time_rep','instantaneous')
      else
         status = nf90_put_att(ncid,varid,'time_rep','averaged')
      endif
      call ice_check_nc(status, subname// ' ERROR: defining time rep for '//hfield%vname, &
                        file=__FILE__, line=__LINE__)

#else
      call abort_ice(subname//' ERROR: USE_NETCDF cpp not defined', file=__FILE__, line=__LINE__)
#endif

      end subroutine ice_hist_field_def

!=======================================================================
! Defines missing_value and _FillValue attributes

      subroutine ice_write_hist_fill(ncid,varid,vname,precision)

      integer (kind=int_kind), intent(in) :: ncid   ! netcdf file id
      integer (kind=int_kind), intent(in) :: varid  ! netcdf var id
      character(len=*),        intent(in) :: vname  ! var name
      integer (kind=int_kind), intent(in) :: precision   ! precision

      ! local variables

      integer (kind=int_kind) :: status
      character(len=*), parameter :: subname = '(ice_write_hist_fill)'

#ifdef USE_NETCDF
      if (precision == 8) then
         status = nf90_put_att(ncid,varid,'missing_value',spval_dbl)
      else
         status = nf90_put_att(ncid,varid,'missing_value',spval)
      endif
      call ice_check_nc(status, subname// ' ERROR: defining missing_value for '//trim(vname), &
                        file=__FILE__, line=__LINE__)

      if (precision == 8) then
         status = nf90_put_att(ncid,varid,'_FillValue',spval_dbl)
      else
         status = nf90_put_att(ncid,varid,'_FillValue',spval)
      endif
      call ice_check_nc(status, subname// ' ERROR: defining _FillValue for '//trim(vname), &
                        file=__FILE__, line=__LINE__)
#else
      call abort_ice(subname//' ERROR : USE_NETCDF cpp not defined', file=__FILE__, line=__LINE__)
#endif

      end subroutine ice_write_hist_fill

!=======================================================================
! Defines a coordinate var in the history file
! coordinates have short_name, long_name and units attributes,
!  and are compressed for 'hdf5' when more than one dimensional

      subroutine ice_hist_coord_def(ncid, coord, lprecision, dimids, varid)

      use ice_history_shared, only: history_deflate, history_format, history_chunksize

      integer(kind=int_kind), intent(in) :: ncid, dimids(:), lprecision
      type(coord_attributes), intent(in) :: coord
      integer(kind=int_kind), intent(inout) :: varid

      !local vars
      integer(kind=int_kind) ::chunks(size(dimids)), i, status

      character(len=*), parameter :: subname = '(ice_hist_coord_def)'

#ifdef USE_NETCDF
      status = nf90_def_var(ncid, coord%short_name, lprecision, dimids, varid)
      call ice_check_nc(status, subname//' ERROR: defining coord '//coord%short_name,file=__FILE__,line=__LINE__)

#ifdef NO_HDF5
      if (history_format=='hdf5') then
           call abort_ice(subname//' ERROR: history_format hdf5 not available ', &
              file=__FILE__, line=__LINE__)
      endif
#else
      if (history_format=='hdf5' .and. size(dimids)>1) then
         if (dimids(1)==imtid .and. dimids(2)==jmtid) then
            chunks(1)=history_chunksize(1)
            chunks(2)=history_chunksize(2)
            do i = 3, size(dimids)
               chunks(i) = 0
            enddo
            status = nf90_def_var_chunking(ncid, varid, NF90_CHUNKED, chunksizes=chunks)
            call ice_check_nc(status, subname//' ERROR chunking var '//trim(coord%short_name), file=__FILE__, line=__LINE__)
         endif
      endif

      if (history_format=='hdf5' .and. history_deflate/=0) then
         status=nf90_def_var_deflate(ncid, varid, shuffle=0, deflate=1, deflate_level=history_deflate)
         call ice_check_nc(status, subname//' ERROR deflating var '//trim(coord%short_name), file=__FILE__, line=__LINE__)
      endif
#endif

      if (coord%long_name(1:3) /= 'und') then
         status = nf90_put_att(ncid,varid,'long_name',trim(coord%long_name))
         call ice_check_nc(status, subname// ' ERROR: defining long_name for '//coord%short_name, &
                           file=__FILE__, line=__LINE__)
      endif
      if (coord%units(1:3) /= 'und') then
         status = nf90_put_att(ncid, varid, 'units', trim(coord%units))
         call ice_check_nc(status, subname// ' ERROR: defining units for '//coord%short_name, &
                           file=__FILE__, line=__LINE__)
      endif
      if (coord%axis(1:3) /= 'und') then
         status = nf90_put_att(ncid, varid, 'axis', trim(coord%axis))
         call ice_check_nc(status, subname// ' ERROR: defining axis for '//coord%short_name, &
                           file=__FILE__, line=__LINE__)
      endif

#else
      call abort_ice(subname//' ERROR: USE_NETCDF cpp not defined', &
            file=__FILE__, line=__LINE__)
#endif

      end subroutine ice_hist_coord_def

!=======================================================================

      end module ice_history_write

!=======================================================================
