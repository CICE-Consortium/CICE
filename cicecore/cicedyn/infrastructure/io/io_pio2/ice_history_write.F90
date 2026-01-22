!=======================================================================
!
! Writes history in netCDF format using NCAR ParallelIO library
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
! 2025 T Craig: Add history restart capability
!
      module ice_history_write

      use ice_kinds_mod
      use ice_constants, only: c0, c360, p5, spval, spval_dbl
      use ice_fileunits, only: nu_diag
      use ice_exit, only: abort_ice
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_parameters
      use ice_calendar, only: write_ic, histfreq
      use ice_pio

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

      public :: ice_write_hist, ice_read_hist

      integer (kind=int_kind) :: imtid,jmtid

!=======================================================================

      contains

!=======================================================================
!
! write average ice quantities or snapshots
! supports history output, write_ic, and history restarts
!
! author:   Elizabeth C. Hunke, LANL

      subroutine ice_write_hist (ns)

      use ice_blocks, only: nx_block, ny_block
      use ice_broadcast, only: broadcast_scalar
      use ice_calendar, only: msec, timesecs, idate, idate0, &
          histfreq_n, days_per_year, use_leap_years, dayyr, &
          hh_init, mm_init, ss_init
      use ice_communicate, only: my_task, master_task
      use ice_domain, only: distrb_info, nblocks
      use ice_domain_size, only: nx_global, ny_global, max_blocks
      use ice_flux, only: albcnt, snwcnt
      use ice_gather_scatter, only: gather_global
      use ice_grid, only: TLON, TLAT, ULON, ULAT, NLON, NLAT, ELON, ELAT, &
          hm, bm, uvm, npm, epm, &
          dxU, dxT, dyU, dyT, dxN, dyN, dxE, dyE, HTN, HTE, ANGLE, ANGLET, &
          tarea, uarea, narea, earea, tmask, umask, nmask, emask, &
          lont_bounds, latt_bounds, lonu_bounds, latu_bounds, &
          lonn_bounds, latn_bounds, lone_bounds, late_bounds
      use ice_history_shared
      use ice_arrays_column, only: hin_max, floe_rad_c
      use ice_restart_shared, only: runid, restart_dir
      use pio

      integer (kind=int_kind), intent(in) :: ns

      ! local variables

      integer (kind=int_kind) :: i,j,k,ic,n,nn, &
         ncid,status,kmtidi,kmtids,kmtidb,cmtid,timid, &
         length,nvertexid,ivertex,kmtida,fmtid,lhistprec
      integer (kind=int_kind), dimension(2) :: dimid2
      integer (kind=int_kind), dimension(3) :: dimid3
      integer (kind=int_kind), dimension(4) :: dimidz
      integer (kind=int_kind), dimension(5) :: dimidcz
      integer (kind=int_kind), dimension(3) :: dimid_nverts
      integer (kind=int_kind), dimension(6) :: dimidex
      real (kind= dbl_kind) :: ltime2
      character (len=8) :: cdate
      character (len=1) :: cns
      character (len=char_len_long) :: title, cal_units, cal_att
      character (len=char_len) :: time_period_freq = 'none'
      character (len=char_len_long) :: ncfile
      character (len=512) :: extvars

      integer (kind=int_kind) :: icategory,ind,i_aice,boundid, lprecision

      character (len=char_len) :: start_time,current_date,current_time
      character (len=16) :: c_aice

      type(file_desc_t)     :: File
      type(io_desc_t)       :: iodesc2d, &
                               iodesc3dc, iodesc3dv, iodesc3di, iodesc3db, iodesc3da, &
                               iodesc3df, &
                               iodesc4di, iodesc4ds, iodesc4df
      type(var_desc_t)      :: varid

      ! time coord
      TYPE(coord_attributes) :: time_coord

      ! 4 vertices in each grid cell
      INTEGER (kind=int_kind), PARAMETER :: nverts = 4

      ! 8 variables describe T, U, N, E grid boundaries:
      ! lont_bounds, latt_bounds, lonu_bounds, latu_bounds
      ! lonn_bounds, latn_bounds, lone_bounds, late_bounds
      INTEGER (kind=int_kind), PARAMETER :: nvar_verts = 8

      TYPE(req_attributes), dimension(nvar_grd) :: var_grd
      TYPE(coord_attributes), dimension(ncoord) :: var_coord
      TYPE(coord_attributes), dimension(nvar_verts) :: var_nverts
      TYPE(coord_attributes), dimension(nvar_grdz) :: var_grdz
      CHARACTER (char_len), dimension(ncoord) :: coord_bounds

      real (kind=dbl_kind) , allocatable :: workd2(:,:,:)
      real (kind=dbl_kind) , allocatable :: workd3(:,:,:,:)
      real (kind=dbl_kind) , allocatable :: workd4(:,:,:,:,:)
      real (kind=dbl_kind) , allocatable :: workd3v(:,:,:,:)

      real (kind=real_kind), allocatable :: workr2(:,:,:)
      real (kind=real_kind), allocatable :: workr3(:,:,:,:)
      real (kind=real_kind), allocatable :: workr4(:,:,:,:,:)
      real (kind=real_kind), allocatable :: workr3v(:,:,:,:)

      character(len=char_len_long) ::  filename

      integer (kind=int_kind), dimension(1) ::  &
         tim_start,tim_length          ! dimension quantities for netCDF

      integer (kind=int_kind), dimension(2) ::  &
         bnd_start,bnd_length          ! dimension quantities for netCDF

      real (kind=dbl_kind) :: secday
      real (kind=dbl_kind) :: rad_to_deg

      logical (kind=log_kind), save :: first_call = .true.

      character(len=*), parameter :: subname = '(ice_write_hist)'

      call icepack_query_parameters(secday_out=secday)
      call icepack_query_parameters(rad_to_deg_out=rad_to_deg)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__, line=__LINE__)

      extvars = ''

      write(cns,'(i1.1)') ns

      ! modify history restart output
      lhistprec = history_precision
      if (write_histrest_now) then
         history_precision = 8
      endif

      if (my_task == master_task) then
        if (write_histrest_now) then
           call construct_filename(ncfile,'nc',ns,option='histrest')
        else
           call construct_filename(ncfile,'nc',ns)
        endif

        ! add local directory path name to ncfile
        if (write_ic) then
          ncfile = trim(incond_dir)//ncfile
        elseif (write_histrest_now) then
          ncfile = trim(restart_dir)//ncfile
        else
          ncfile = trim(history_dir)//ncfile
        endif
        filename = ncfile
      end if
      call broadcast_scalar(filename, master_task)

      ! create file
      File%fh=-1
      call ice_pio_init(mode='write', filename=trim(filename), File=File, &
           clobber=.true., fformat=trim(history_format), rearr=trim(history_rearranger), &
           iotasks=history_iotasks, root=history_root, stride=history_stride, debug=first_call)

      call ice_pio_initdecomp(iodesc=iodesc2d, precision=history_precision)
      call ice_pio_initdecomp(ndim3=ncat_hist, iodesc=iodesc3dc, precision=history_precision)
      call ice_pio_initdecomp(ndim3=nzilyr,    iodesc=iodesc3di, precision=history_precision)
      call ice_pio_initdecomp(ndim3=nzblyr,    iodesc=iodesc3db, precision=history_precision)
      call ice_pio_initdecomp(ndim3=nzalyr,    iodesc=iodesc3da, precision=history_precision)
      call ice_pio_initdecomp(ndim3=nfsd_hist, iodesc=iodesc3df, precision=history_precision)
      call ice_pio_initdecomp(ndim3=nverts,    iodesc=iodesc3dv, inner_dim=.true., precision=history_precision)
      call ice_pio_initdecomp(ndim3=nzilyr,    ndim4=ncat_hist, iodesc=iodesc4di, precision=history_precision)
      call ice_pio_initdecomp(ndim3=nzslyr,    ndim4=ncat_hist, iodesc=iodesc4ds, precision=history_precision)
      call ice_pio_initdecomp(ndim3=nfsd_hist, ndim4=ncat_hist, iodesc=iodesc4df, precision=history_precision)

      ! option of turning on double precision history files
      lprecision = pio_real
      if (history_precision == 8) lprecision = pio_double

      !-----------------------------------------------------------------
      ! define dimensions
      !-----------------------------------------------------------------
      call pio_seterrorhandling(File, PIO_RETURN_ERROR)

      if (hist_avg(ns) .and. .not. write_ic) then
         call ice_pio_check(pio_def_dim(File,'nbnd',2,boundid), &
              subname//' ERROR: defining dim nbnd with len 2',file=__FILE__,line=__LINE__)
      endif

      call ice_pio_check(pio_def_dim(File,'ni',nx_global,imtid), &
           subname//' ERROR: defining dim ni',file=__FILE__,line=__LINE__)

      call ice_pio_check(pio_def_dim(File,'nj',ny_global,jmtid), &
           subname//' ERROR: defining dim nj',file=__FILE__,line=__LINE__)

      call ice_pio_check(pio_def_dim(File,'nc',ncat_hist,cmtid), &
           subname//' ERROR: defining dim nc',file=__FILE__,line=__LINE__)

      call ice_pio_check(pio_def_dim(File,'nkice',nzilyr,kmtidi), &
           subname//' ERROR: defining dim nkice',file=__FILE__,line=__LINE__)

      call ice_pio_check(pio_def_dim(File,'nksnow',nzslyr,kmtids), &
           subname//' ERROR: defining dim nksnow',file=__FILE__,line=__LINE__)

      call ice_pio_check(pio_def_dim(File,'nkbio',nzblyr,kmtidb), &
           subname//' ERROR: defining dim nkbio',file=__FILE__,line=__LINE__)

      call ice_pio_check(pio_def_dim(File,'nkaer',nzalyr,kmtida), &
           subname//' ERROR: defining dim nkaer',file=__FILE__,line=__LINE__)

      ! do not write time axis on grid output file
      timid = -99
      if (histfreq(ns)/='g') then
         call ice_pio_check(pio_def_dim(File,'time',PIO_UNLIMITED,timid), &
              subname//' ERROR: defining dim time',file=__FILE__,line=__LINE__)
      endif

      call ice_pio_check(pio_def_dim(File,'nvertices',nverts,nvertexid), &
           subname//' ERROR: defining dim nvertices',file=__FILE__,line=__LINE__)

      call ice_pio_check(pio_def_dim(File,'nf',nfsd_hist,fmtid), &
           subname//' ERROR: defining dim nf',file=__FILE__,line=__LINE__)

      !-----------------------------------------------------------------
      ! define coordinate variables:  time, time_bounds
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
            call abort_ice(subname//' ERROR: invalid calendar settings')
         endif

         time_coord = coord_attributes('time', 'time', trim(cal_units), 'T')
         call ice_hist_coord_def(File, time_coord, pio_double, (/timid/), varid)
         call ice_pio_check(pio_put_att(File,varid,'calendar',cal_att), &
                 subname//' ERROR: defining att calendar: '//cal_att,file=__FILE__,line=__LINE__)
         if (hist_avg(ns) .and. .not. write_ic) then
            call ice_pio_check(pio_put_att(File,varid,'bounds','time_bounds'), &
                 subname//' ERROR: defining att bounds time_bounds',file=__FILE__,line=__LINE__)
         endif

         ! Define coord time_bounds if hist_avg is true
         if (hist_avg(ns) .and. .not. write_ic) then
            time_coord = coord_attributes('time_bounds', 'time interval bounds', trim(cal_units), 'undefined')

            dimid2(1) = boundid
            dimid2(2) = timid

            call ice_hist_coord_def(File, time_coord, pio_double, dimid2, varid)
            call ice_pio_check(pio_put_att(File,varid,'calendar',cal_att), &
                 subname//' ERROR: defining att calendar: '//cal_att,file=__FILE__,line=__LINE__)
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

      var_grd(n_blkmask)%req = coord_attributes('blkmask', &
                  'block id of T grid cells, mytask + iblk/100', '1', 'undefined')
      var_grd(n_blkmask)%coordinates = 'TLON TLAT'

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
      var_nverts(n_lont_bnds) = coord_attributes('lont_bounds','longitude bounds (T-cell)','degrees_east','und')
      var_nverts(n_latt_bnds) = coord_attributes('latt_bounds','latitude bounds (T-cell)','degrees_north','und')
      var_nverts(n_lonu_bnds) = coord_attributes('lonu_bounds','longitude bounds (U-cell)','degrees_east','und')
      var_nverts(n_latu_bnds) = coord_attributes('latu_bounds','latitude bounds (U-cell)','degrees_north','und')
      var_nverts(n_lonn_bnds) = coord_attributes('lonn_bounds','longitude bounds (N-cell)','degrees_east','und')
      var_nverts(n_latn_bnds) = coord_attributes('latn_bounds','latitude bounds (N-cell)','degrees_north','und')
      var_nverts(n_lone_bnds) = coord_attributes('lone_bounds','longitude bounds (E-cell)','degrees_east','und')
      var_nverts(n_late_bnds) = coord_attributes('late_bounds','latitude bounds (E-cell)','degrees_north','und')

      !-----------------------------------------------------------------
      ! define attributes for time-invariant variables
      !-----------------------------------------------------------------

      dimid2(1) = imtid
      dimid2(2) = jmtid

      do i = 1, ncoord
         if (icoord(i) .or. histfreq(ns)=='g') then
            call ice_hist_coord_def(File, var_coord(i), lprecision, dimid2, varid)
            call ice_write_hist_fill(File,varid,var_coord(i)%short_name,history_precision)
            if (var_coord(i)%short_name == 'ULAT') then
               call ice_pio_check(pio_put_att(File,varid,'comment', &
                    trim('Latitude of NE corner of T grid cell')), &
                    subname//' ERROR: defining att comment',file=__FILE__,line=__LINE__)
            endif
            if (f_bounds .or. histfreq(ns)=='g') then
               call ice_pio_check(pio_put_att(File, varid, 'bounds', trim(coord_bounds(i))), &
                    subname//' ERROR: defining att bounds '//trim(coord_bounds(i)),file=__FILE__,line=__LINE__)
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
            call ice_hist_coord_def(File, var_grdz(i), lprecision, dimidex(i:i), varid)
         else
            extvars = trim(extvars)//' '//trim(var_grdz(i)%short_name)
         endif
      enddo

      do i = 1, nvar_grd
         if (igrd(i) .or. histfreq(ns)=='g') then
            call ice_hist_coord_def(File, var_grd(i)%req, lprecision, dimid2, varid)
            call ice_pio_check(pio_put_att(File, varid, 'coordinates', trim(var_grd(i)%coordinates)), &
                 subname//' ERROR: defining att coordinates '//trim(var_grd(i)%coordinates),file=__FILE__,line=__LINE__)
            call ice_write_hist_fill(File,varid,var_grd(i)%req%short_name,history_precision)
         else
            extvars = trim(extvars)//' '//trim(var_grd(i)%req%short_name)
         endif
      enddo

      ! bounds fields with dimensions (nverts,nx,ny)
      dimid_nverts(1) = nvertexid
      dimid_nverts(2) = imtid
      dimid_nverts(3) = jmtid
      do i = 1, nvar_verts
         if (f_bounds .or. histfreq(ns)=='g') then
            call ice_hist_coord_def(File, var_nverts(i), lprecision, dimid_nverts, varid)
         endif
      enddo

      !-----------------------------------------------------------------
      ! define attributes for time-variant variables
      !-----------------------------------------------------------------

      ! 2D
      dimid3(1) = imtid
      dimid3(2) = jmtid
      dimid3(3) = timid

      if (write_histrest_now) then
         status = pio_def_var(File, 'time_beg', lprecision, varid)
         status = pio_def_var(File, 'avgct', lprecision, varid)
         status = pio_def_var(File, 'albcnt'//cns, lprecision, dimid3, varid)
         status = pio_def_var(File, 'snwcnt'//cns, lprecision, dimid3, varid)
      endif

      do n=1,num_avail_hist_fields_2D
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            call ice_hist_field_def(File, avail_hist_fields(n),lprecision, dimid3, ns)
         endif
      enddo

      ! 3D (category)
      dimidz(1) = imtid
      dimidz(2) = jmtid
      dimidz(3) = cmtid
      dimidz(4) = timid

      do n = n2D + 1, n3Dccum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            call ice_hist_field_def(File, avail_hist_fields(n),lprecision, dimidz,ns)
         endif
      enddo  ! num_avail_hist_fields_3Dc

      ! 3D (ice layers)
      dimidz(1) = imtid
      dimidz(2) = jmtid
      dimidz(3) = kmtidi
      dimidz(4) = timid

      do n = n3Dccum + 1, n3Dzcum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            call ice_hist_field_def(File, avail_hist_fields(n),lprecision, dimidz,ns)
         endif
      enddo  ! num_avail_hist_fields_3Dz

      ! 3D (biology ice layers)
      dimidz(1) = imtid
      dimidz(2) = jmtid
      dimidz(3) = kmtidb
      dimidz(4) = timid

      do n = n3Dzcum + 1, n3Dbcum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            call ice_hist_field_def(File, avail_hist_fields(n),lprecision, dimidz,ns)
         endif
      enddo  ! num_avail_hist_fields_3Db

      ! 3D (biology snow layers)
      dimidz(1) = imtid
      dimidz(2) = jmtid
      dimidz(3) = kmtida
      dimidz(4) = timid

      do n = n3Dbcum + 1, n3Dacum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            call ice_hist_field_def(File, avail_hist_fields(n),lprecision, dimidz,ns)
         endif
      enddo  ! num_avail_hist_fields_3Da

      ! 3D (fsd)
      dimidz(1) = imtid
      dimidz(2) = jmtid
      dimidz(3) = fmtid
      dimidz(4) = timid

      do n = n3Dacum + 1, n3Dfcum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            call ice_hist_field_def(File, avail_hist_fields(n),lprecision, dimidz,ns)
         endif
      enddo  ! num_avail_hist_fields_3Df

      ! 4D (ice categories)
      dimidcz(1) = imtid
      dimidcz(2) = jmtid
      dimidcz(3) = kmtidi
      dimidcz(4) = cmtid
      dimidcz(5) = timid

      do n = n3Dfcum + 1, n4Dicum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            call ice_hist_field_def(File, avail_hist_fields(n),lprecision, dimidcz,ns)
         endif
      enddo  ! num_avail_hist_fields_4Di

      ! 4D (snow layers)
      dimidcz(1) = imtid
      dimidcz(2) = jmtid
      dimidcz(3) = kmtids
      dimidcz(4) = cmtid
      dimidcz(5) = timid

      do n = n4Dicum + 1, n4Dscum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            call ice_hist_field_def(File, avail_hist_fields(n),lprecision, dimidcz,ns)
         endif
      enddo  ! num_avail_hist_fields_4Ds

      ! 4D (fsd layers)
      dimidcz(1) = imtid
      dimidcz(2) = jmtid
      dimidcz(3) = fmtid
      dimidcz(4) = cmtid
      dimidcz(5) = timid

      do n = n4Dscum + 1, n4Dfcum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            call ice_hist_field_def(File, avail_hist_fields(n),lprecision, dimidcz,ns)
         endif
      enddo  ! num_avail_hist_fields_4Df

      !-----------------------------------------------------------------
      ! global attributes
      !-----------------------------------------------------------------
      ! ... the user should change these to something useful ...
      !-----------------------------------------------------------------
#ifdef CESMCOUPLED
      call ice_pio_check(pio_put_att(File,pio_global,'title',runid), &
           subname//' ERROR: defining att title '//runid,file=__FILE__,line=__LINE__)
#else
      title  = 'sea ice model output for CICE'
      call ice_pio_check(pio_put_att(File,pio_global,'title',trim(title)), &
           subname//' ERROR: defining att title '//trim(title),file=__FILE__,line=__LINE__)
#endif
      title = 'Diagnostic and Prognostic Variables'
      call ice_pio_check(pio_put_att(File,pio_global,'contents',trim(title)), &
           subname//' ERROR: defining att contents '//trim(title),file=__FILE__,line=__LINE__)

      write(title,'(2a)') 'CICE Sea Ice Model, ', trim(version_name)
      call ice_pio_check(pio_put_att(File,pio_global,'source',trim(title)), &
           subname//' ERROR: defining att source '//trim(title),file=__FILE__,line=__LINE__)

      if (use_leap_years) then
         write(title,'(a,i3,a)') 'This year has ',dayyr,' days'
      else
         write(title,'(a,i3,a)') 'All years have exactly ',dayyr,' days'
      endif
      call ice_pio_check(pio_put_att(File,pio_global,'comment',trim(title)), &
           subname//' ERROR: defining att comment '//trim(title),file=__FILE__,line=__LINE__)

      write(title,'(a,i8.8)') 'File written on model date ',idate
      call ice_pio_check(pio_put_att(File,pio_global,'comment2',trim(title)), &
           subname//' ERROR: defining att comment2 '//trim(title),file=__FILE__,line=__LINE__)

      write(title,'(a,i6)') 'seconds elapsed into model date: ',msec
      call ice_pio_check(pio_put_att(File,pio_global,'comment3',trim(title)), &
           subname//' ERROR: defining att comment3 '//trim(title),file=__FILE__,line=__LINE__)

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
         call ice_pio_check(pio_put_att(File,pio_global,'time_period_freq',trim(time_period_freq)), &
              subname//' ERROR: defining att time_period_freq '//trim(time_period_freq),file=__FILE__,line=__LINE__)
      endif

      if (hist_avg(ns)) &
         call ice_pio_check(pio_put_att(File,pio_global,'time_axis_position',trim(hist_time_axis)), &
              subname//' ERROR: defining att time_axis_position '//trim(hist_time_axis),file=__FILE__,line=__LINE__)

      title = 'CF-1.8'
      call ice_pio_check(pio_put_att(File,pio_global,'Conventions',trim(title)), &
           subname//' ERROR: defining att conventions '//trim(title),file=__FILE__,line=__LINE__)

      call ice_pio_check(pio_put_att(File,pio_global,'external_variables',trim(extvars)), &
           subname//' ERROR: defining att external_variables '//trim(extvars),file=__FILE__,line=__LINE__)

      call date_and_time(date=current_date, time=current_time)
      write(start_time,1000) current_date(1:4), current_date(5:6), &
                             current_date(7:8), current_time(1:2), &
                             current_time(3:4)
1000  format('This dataset was created on ', &
              a,'-',a,'-',a,' at ',a,':',a)
      call ice_pio_check(pio_put_att(File,pio_global,'history',trim(start_time)), &
           subname//' ERROR: defining att history '//trim(start_time),file=__FILE__,line=__LINE__)

      write(start_time,1001) current_date(1:4), current_date(5:6), &
                             current_date(7:8), current_time(1:2), &
                             current_time(3:4)
1001  format(a,'-',a,'-',a,' ',a,':',a)
      call ice_pio_check(pio_put_att(File,pio_global,'date_created',trim(start_time)), &
           subname//' ERROR: defining att date_created '//trim(start_time),file=__FILE__,line=__LINE__)

#ifdef USE_PIO1
      call ice_pio_check(pio_put_att(File,pio_global,'io_flavor','io_pio1 '//trim(history_format)), &
           subname//' ERROR: defining att io_flavor',file=__FILE__,line=__LINE__)
#else
      call ice_pio_check(pio_put_att(File,pio_global,'io_flavor','io_pio2 '//trim(history_format)), &
           subname//' ERROR: defining att io_flavor',file=__FILE__,line=__LINE__)
#endif

      !-----------------------------------------------------------------
      ! end define mode
      !-----------------------------------------------------------------

      call ice_pio_check(pio_enddef(File), &
           subname//' ERROR: ending pio definitions',file=__FILE__,line=__LINE__)

      !-----------------------------------------------------------------
      ! write time and time bounds info
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

         call ice_pio_check(pio_inq_varid(File,'time',varid), &
              subname//' ERROR: getting var time',file=__FILE__,line=__LINE__)
         call ice_pio_check(pio_put_var(File,varid,(/1/),ltime2), &
              subname//' ERROR: setting var time',file=__FILE__,line=__LINE__)

         if (hist_avg(ns) .and. .not. write_ic) then
            call ice_pio_check(pio_inq_varid(File,'time_bounds',varid), &
                 subname//' ERROR: getting time_bounds' ,file=__FILE__,line=__LINE__)
            time_bounds=(/time_beg(ns),time_end(ns)/)
            bnd_start  = (/1,1/)
            bnd_length = (/2,1/)
            call ice_pio_check(pio_put_var(File,varid,ival=time_bounds,start=bnd_start(:),count=bnd_length(:)), &
                 subname//' ERROR: setting time_bounds' ,file=__FILE__,line=__LINE__)
         endif
      endif  ! histfreq(ns)/='g'

      !-----------------------------------------------------------------
      ! write coordinate variables
      !-----------------------------------------------------------------

      allocate(workd2(nx_block,ny_block,nblocks))
      allocate(workr2(nx_block,ny_block,nblocks))

      do i = 1,ncoord
         if(icoord(i) .or. histfreq(ns)=='g') then
            call ice_pio_check(pio_inq_varid(File, var_coord(i)%short_name, varid), &
                 subname//' ERROR: getting '//var_coord(i)%short_name ,file=__FILE__,line=__LINE__)
            SELECT CASE (var_coord(i)%short_name)
               CASE ('TLON')
                 ! Convert T grid longitude from -180 -> 180 to 0 to 360
                    workd2(:,:,:) = mod(tlon(:,:,1:nblocks)*rad_to_deg + c360, c360)
               CASE ('TLAT')
                 workd2(:,:,:) = tlat(:,:,1:nblocks)*rad_to_deg
               CASE ('ULON')
                 workd2(:,:,:) = ulon(:,:,1:nblocks)*rad_to_deg
               CASE ('ULAT')
                 workd2(:,:,:) = ulat(:,:,1:nblocks)*rad_to_deg
               CASE ('NLON')
                 workd2(:,:,:) = nlon(:,:,1:nblocks)*rad_to_deg
               CASE ('NLAT')
                 workd2(:,:,:) = nlat(:,:,1:nblocks)*rad_to_deg
               CASE ('ELON')
                 workd2(:,:,:) = elon(:,:,1:nblocks)*rad_to_deg
               CASE ('ELAT')
                 workd2(:,:,:) = elat(:,:,1:nblocks)*rad_to_deg
            END SELECT
            if (history_precision == 8) then
               call pio_write_darray(File, varid, iodesc2d, &
                    workd2, status, fillval=spval_dbl)
            else
               workr2 = workd2
               call pio_write_darray(File, varid, iodesc2d, &
                    workr2, status, fillval=spval)
            endif

            call ice_pio_check(status,subname//' ERROR: writing '//avail_hist_fields(n)%vname, &
                               file=__FILE__,line=__LINE__)
         endif
      enddo

      ! Extra dimensions (NCAT, NFSD, VGRD*)

      do i = 1, nvar_grdz
         if (igrdz(i) .or. histfreq(ns)=='g') then
            call ice_pio_check(pio_inq_varid(File, var_grdz(i)%short_name, varid), &
                 subname//' ERROR: getting '//var_grdz(i)%short_name,file=__FILE__,line=__LINE__)
            SELECT CASE (var_grdz(i)%short_name)
               CASE ('NCAT')
                  call ice_pio_check(pio_put_var(File, varid, hin_max(1:ncat_hist)), &
                       subname//' ERROR: setting '//var_grdz(i)%short_name,file=__FILE__,line=__LINE__)
               CASE ('NFSD')
                  call ice_pio_check(pio_put_var(File, varid, floe_rad_c(1:nfsd_hist)), &
                       subname//' ERROR: setting '//var_grdz(i)%short_name,file=__FILE__,line=__LINE__)
               CASE ('VGRDi')
                  call ice_pio_check(pio_put_var(File, varid, (/(k, k=1,nzilyr)/)), &
                       subname//' ERROR: setting '//var_grdz(i)%short_name,file=__FILE__,line=__LINE__)
               CASE ('VGRDs')
                  call ice_pio_check(pio_put_var(File, varid, (/(k, k=1,nzslyr)/)), &
                       subname//' ERROR: setting '//var_grdz(i)%short_name,file=__FILE__,line=__LINE__)
               CASE ('VGRDb')
                  call ice_pio_check(pio_put_var(File, varid, (/(k, k=1,nzblyr)/)), &
                       subname//' ERROR: setting '//var_grdz(i)%short_name,file=__FILE__,line=__LINE__)
               CASE ('VGRDa')
                  call ice_pio_check(pio_put_var(File, varid, (/(k, k=1,nzalyr)/)), &
                       subname//' ERROR: setting '//var_grdz(i)%short_name,file=__FILE__,line=__LINE__)
            END SELECT
         endif
      enddo

      !-----------------------------------------------------------------
      ! write grid masks, area and rotation angle
      !-----------------------------------------------------------------

      do i = 1, nvar_grd
         if (igrd(i) .or. histfreq(ns)=='g') then
            SELECT CASE (var_grd(i)%req%short_name)
               CASE ('tmask')
                  workd2 = hm(:,:,1:nblocks)
               CASE ('umask')
                  workd2 = uvm(:,:,1:nblocks)
               CASE ('nmask')
                  workd2 = npm(:,:,1:nblocks)
               CASE ('emask')
                  workd2 = epm(:,:,1:nblocks)
               CASE ('blkmask')
                  workd2 = bm(:,:,1:nblocks)
               CASE ('tarea')
                  workd2 = tarea(:,:,1:nblocks)
               CASE ('uarea')
                  workd2 = uarea(:,:,1:nblocks)
               CASE ('narea')
                  workd2 = narea(:,:,1:nblocks)
               CASE ('earea')
                  workd2 = earea(:,:,1:nblocks)
               CASE ('dxt')
                  workd2 = dxT(:,:,1:nblocks)
               CASE ('dyt')
                  workd2 = dyT(:,:,1:nblocks)
               CASE ('dxu')
                  workd2 = dxU(:,:,1:nblocks)
               CASE ('dyu')
                  workd2 = dyU(:,:,1:nblocks)
               CASE ('dxn')
                  workd2 = dxN(:,:,1:nblocks)
               CASE ('dyn')
                  workd2 = dyN(:,:,1:nblocks)
               CASE ('dxe')
                  workd2 = dxE(:,:,1:nblocks)
               CASE ('dye')
                  workd2 = dyE(:,:,1:nblocks)
               CASE ('HTN')
                  workd2 = HTN(:,:,1:nblocks)
               CASE ('HTE')
                  workd2 = HTE(:,:,1:nblocks)
               CASE ('ANGLE')
                  workd2 = ANGLE(:,:,1:nblocks)
               CASE ('ANGLET')
                  workd2 = ANGLET(:,:,1:nblocks)
            END SELECT
            call ice_pio_check(pio_inq_varid(File, var_grd(i)%req%short_name, varid), &
                 subname//' ERROR: getting '//var_grd(i)%req%short_name,file=__FILE__,line=__LINE__)
            if (history_precision == 8) then
               call pio_write_darray(File, varid, iodesc2d, &
                    workd2, status, fillval=spval_dbl)
            else
               workr2 = workd2
               call pio_write_darray(File, varid, iodesc2d, &
                    workr2, status, fillval=spval)
            endif

            call ice_pio_check(status,subname//' ERROR: writing '//avail_hist_fields(n)%vname, &
                               file=__FILE__,line=__LINE__)
         endif
      enddo

      !----------------------------------------------------------------
      ! Write coordinates of grid box vertices
      !----------------------------------------------------------------

      if (f_bounds .or. histfreq(ns)=='g') then
         allocate(workd3v(nverts,nx_block,ny_block,nblocks))
         allocate(workr3v(nverts,nx_block,ny_block,nblocks))
         workd3v (:,:,:,:) = c0
         do i = 1, nvar_verts
            SELECT CASE (var_nverts(i)%short_name)
               CASE ('lont_bounds')
                  do ivertex = 1, nverts
                     workd3v(ivertex,:,:,:) = lont_bounds(ivertex,:,:,1:nblocks)
                  enddo
               CASE ('latt_bounds')
                  do ivertex = 1, nverts
                     workd3v(ivertex,:,:,:) = latt_bounds(ivertex,:,:,1:nblocks)
                  enddo
               CASE ('lonu_bounds')
                  do ivertex = 1, nverts
                     workd3v(ivertex,:,:,:) = lonu_bounds(ivertex,:,:,1:nblocks)
                  enddo
               CASE ('latu_bounds')
                  do ivertex = 1, nverts
                     workd3v(ivertex,:,:,:) = latu_bounds(ivertex,:,:,1:nblocks)
                  enddo
               CASE ('lonn_bounds')
                  do ivertex = 1, nverts
                     workd3v(ivertex,:,:,:) = lonn_bounds(ivertex,:,:,1:nblocks)
                  enddo
               CASE ('latn_bounds')
                  do ivertex = 1, nverts
                     workd3v(ivertex,:,:,:) = latn_bounds(ivertex,:,:,1:nblocks)
                  enddo
               CASE ('lone_bounds')
                  do ivertex = 1, nverts
                     workd3v(ivertex,:,:,:) = lone_bounds(ivertex,:,:,1:nblocks)
                  enddo
               CASE ('late_bounds')
                  do ivertex = 1, nverts
                     workd3v(ivertex,:,:,:) = late_bounds(ivertex,:,:,1:nblocks)
                  enddo
            END SELECT

            call ice_pio_check(pio_inq_varid(File, var_nverts(i)%short_name, varid), &
                 subname//' ERROR: getting '//var_nverts(i)%short_name,file=__FILE__,line=__LINE__)
            if (history_precision == 8) then
               call pio_write_darray(File, varid, iodesc3dv, &
                                   workd3v, status, fillval=spval_dbl)
            else
               workr3v = workd3v
               call pio_write_darray(File, varid, iodesc3dv, &
                                     workr3v, status, fillval=spval)
            endif

            call ice_pio_check(status,subname//' ERROR: writing '//avail_hist_fields(n)%vname, &
                               file=__FILE__,line=__LINE__)
         enddo
         deallocate(workd3v)
         deallocate(workr3v)
      endif  ! f_bounds

      !-----------------------------------------------------------------
      ! write variable data
      !-----------------------------------------------------------------

      if (write_histrest_now) then
         call ice_pio_check(pio_inq_varid(File,'time_beg',varid), &
                            subname// ' ERROR: getting varid for '//'time_beg', &
                            file=__FILE__, line=__LINE__)
         call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)
#ifdef CESM1_PIO
         call pio_setframe(varid, int(1,kind=PIO_OFFSET))
#else
         call pio_setframe(File, varid, int(1,kind=PIO_OFFSET_KIND))
#endif
         call pio_seterrorhandling(File, PIO_RETURN_ERROR)
         call ice_pio_check(pio_put_var(File,varid,(/1/),time_beg(ns)), &
                            subname// ' ERROR: writing variable '//'time_beg', &
                            file=__FILE__, line=__LINE__)

         call ice_pio_check(pio_inq_varid(File,'avgct',varid), &
                            subname// ' ERROR: getting varid for '//'avgct', &
                            file=__FILE__, line=__LINE__)
         call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)
#ifdef CESM1_PIO
         call pio_setframe(varid, int(1,kind=PIO_OFFSET))
#else
         call pio_setframe(File, varid, int(1,kind=PIO_OFFSET_KIND))
#endif
         call pio_seterrorhandling(File, PIO_RETURN_ERROR)
         call ice_pio_check(pio_put_var(File,varid,(/1/),avgct(ns)), &
                            subname// ' ERROR: writing variable '//'avgct', &
                            file=__FILE__, line=__LINE__)

         call ice_pio_check(pio_inq_varid(File,'albcnt'//cns,varid), &
              subname//' ERROR: getting varid for '//'albcnt'//cns,file=__FILE__,line=__LINE__)
         workd2(:,:,:) = albcnt(:,:,1:nblocks,ns)
         call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)
#ifdef CESM1_PIO
         call pio_setframe(varid, int(1,kind=PIO_OFFSET))
#else
         call pio_setframe(File, varid, int(1,kind=PIO_OFFSET_KIND))
#endif
         call pio_seterrorhandling(File, PIO_RETURN_ERROR)
         if (history_precision == 8) then
            call pio_write_darray(File, varid, iodesc2d,&
                                  workd2, status, fillval=spval_dbl)
         else
            workr2 = workd2
            call pio_write_darray(File, varid, iodesc2d,&
                                  workr2, status, fillval=spval)
         endif
         call ice_pio_check(status,subname//' ERROR: writing '//'albcnt'//cns, &
                            file=__FILE__,line=__LINE__)

         call ice_pio_check(pio_inq_varid(File,'snwcnt'//cns,varid), &
              subname//' ERROR: getting varid for '//'snwcnt'//cns,file=__FILE__,line=__LINE__)
         workd2(:,:,:) = snwcnt(:,:,1:nblocks,ns)
         call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)
#ifdef CESM1_PIO
         call pio_setframe(varid, int(1,kind=PIO_OFFSET))
#else
         call pio_setframe(File, varid, int(1,kind=PIO_OFFSET_KIND))
#endif
         call pio_seterrorhandling(File, PIO_RETURN_ERROR)
         if (history_precision == 8) then
            call pio_write_darray(File, varid, iodesc2d,&
                                  workd2, status, fillval=spval_dbl)
         else
            workr2 = workd2
            call pio_write_darray(File, varid, iodesc2d,&
                                  workr2, status, fillval=spval)
         endif
         call ice_pio_check(status,subname//' ERROR: writing '//'snwcnt'//cns, &
                            file=__FILE__,line=__LINE__)

      endif

      ! 2D
      do n=1,num_avail_hist_fields_2D
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            call ice_pio_check(pio_inq_varid(File,avail_hist_fields(n)%vname,varid), &
                 subname//' ERROR: getting varid for '//avail_hist_fields(n)%vname,file=__FILE__,line=__LINE__)
            workd2(:,:,:) = a2D(:,:,n,1:nblocks)
            call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)
#ifdef CESM1_PIO
            call pio_setframe(varid, int(1,kind=PIO_OFFSET))
#else
            call pio_setframe(File, varid, int(1,kind=PIO_OFFSET_KIND))
#endif
            call pio_seterrorhandling(File, PIO_RETURN_ERROR)
            if (history_precision == 8) then
               call pio_write_darray(File, varid, iodesc2d,&
                                     workd2, status, fillval=spval_dbl)
            else
               workr2 = workd2
               call pio_write_darray(File, varid, iodesc2d,&
                                     workr2, status, fillval=spval)
            endif

            call ice_pio_check(status,subname//' ERROR: writing '//avail_hist_fields(n)%vname, &
                               file=__FILE__,line=__LINE__)
         endif
      enddo ! num_avail_hist_fields_2D

      deallocate(workd2)
      deallocate(workr2)

      ! 3D (category)
      allocate(workd3(nx_block,ny_block,nblocks,ncat_hist))
      allocate(workr3(nx_block,ny_block,nblocks,ncat_hist))
      do n = n2D + 1, n3Dccum
         nn = n - n2D
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            call ice_pio_check(pio_inq_varid(File,avail_hist_fields(n)%vname,varid), &
                 subname//' ERROR: getting varid for '//avail_hist_fields(n)%vname,file=__FILE__,line=__LINE__)
            do j = 1, nblocks
            do i = 1, ncat_hist
               workd3(:,:,j,i) = a3Dc(:,:,i,nn,j)
            enddo
            enddo
            call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)
#ifdef CESM1_PIO
            call pio_setframe(varid, int(1,kind=PIO_OFFSET))
#else
            call pio_setframe(File, varid, int(1,kind=PIO_OFFSET_KIND))
#endif
            call pio_seterrorhandling(File, PIO_RETURN_ERROR)
            if (history_precision == 8) then
               call pio_write_darray(File, varid, iodesc3dc,&
                                     workd3, status, fillval=spval_dbl)
            else
               workr3 = workd3
               call pio_write_darray(File, varid, iodesc3dc,&
                                     workr3, status, fillval=spval)
            endif

            call ice_pio_check(status,subname//' ERROR: writing '//avail_hist_fields(n)%vname, &
                               file=__FILE__,line=__LINE__)
         endif
      enddo ! num_avail_hist_fields_3Dc
      deallocate(workd3)
      deallocate(workr3)

      ! 3D (vertical ice)
      allocate(workd3(nx_block,ny_block,nblocks,nzilyr))
      allocate(workr3(nx_block,ny_block,nblocks,nzilyr))
      do n = n3Dccum+1, n3Dzcum
         nn = n - n3Dccum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            call ice_pio_check(pio_inq_varid(File,avail_hist_fields(n)%vname,varid), &
                 subname//' ERROR: getting varid for '//avail_hist_fields(n)%vname,file=__FILE__,line=__LINE__)
            do j = 1, nblocks
            do i = 1, nzilyr
               workd3(:,:,j,i) = a3Dz(:,:,i,nn,j)
            enddo
            enddo
            call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)
#ifdef CESM1_PIO
            call pio_setframe(varid, int(1,kind=PIO_OFFSET))
#else
            call pio_setframe(File, varid, int(1,kind=PIO_OFFSET_KIND))
#endif
            call pio_seterrorhandling(File, PIO_RETURN_ERROR)
            if (history_precision == 8) then
               call pio_write_darray(File, varid, iodesc3di,&
                                     workd3, status, fillval=spval_dbl)
            else
               workr3 = workd3
               call pio_write_darray(File, varid, iodesc3di,&
                                     workr3, status, fillval=spval)
            endif

            call ice_pio_check(status,subname//' ERROR: writing '//avail_hist_fields(n)%vname, &
                               file=__FILE__,line=__LINE__)
         endif
      enddo ! num_avail_hist_fields_3Dz
      deallocate(workd3)
      deallocate(workr3)

      ! 3D (vertical ice biology)
      allocate(workd3(nx_block,ny_block,nblocks,nzblyr))
      allocate(workr3(nx_block,ny_block,nblocks,nzblyr))
      do n = n3Dzcum+1, n3Dbcum
         nn = n - n3Dzcum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            call ice_pio_check(pio_inq_varid(File,avail_hist_fields(n)%vname,varid), &
                 subname//' ERROR: getting varid for '//avail_hist_fields(n)%vname,file=__FILE__,line=__LINE__)
            do j = 1, nblocks
            do i = 1, nzblyr
               workd3(:,:,j,i) = a3Db(:,:,i,nn,j)
            enddo
            enddo
            call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)
#ifdef CESM1_PIO
            call pio_setframe(varid, int(1,kind=PIO_OFFSET))
#else
            call pio_setframe(File, varid, int(1,kind=PIO_OFFSET_KIND))
#endif
            call pio_seterrorhandling(File, PIO_RETURN_ERROR)
            if (history_precision == 8) then
               call pio_write_darray(File, varid, iodesc3db,&
                                     workd3, status, fillval=spval_dbl)
            else
               workr3 = workd3
               call pio_write_darray(File, varid, iodesc3db,&
                                     workr3, status, fillval=spval)
            endif

            call ice_pio_check(status,subname//' ERROR: writing '//avail_hist_fields(n)%vname, &
                               file=__FILE__,line=__LINE__)
         endif
      enddo ! num_avail_hist_fields_3Db
      deallocate(workd3)
      deallocate(workr3)

      ! 3D (vertical snow biology)
      allocate(workd3(nx_block,ny_block,nblocks,nzalyr))
      allocate(workr3(nx_block,ny_block,nblocks,nzalyr))
      do n = n3Dbcum+1, n3Dacum
         nn = n - n3Dbcum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            call ice_pio_check(pio_inq_varid(File,avail_hist_fields(n)%vname,varid), &
                 subname//' ERROR: getting varid for '//avail_hist_fields(n)%vname,file=__FILE__,line=__LINE__)
            do j = 1, nblocks
            do i = 1, nzalyr
               workd3(:,:,j,i) = a3Da(:,:,i,nn,j)
            enddo
            enddo
            call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)
#ifdef CESM1_PIO
            call pio_setframe(varid, int(1,kind=PIO_OFFSET))
#else
            call pio_setframe(File, varid, int(1,kind=PIO_OFFSET_KIND))
#endif
            call pio_seterrorhandling(File, PIO_RETURN_ERROR)
            if (history_precision == 8) then
               call pio_write_darray(File, varid, iodesc3da,&
                                     workd3, status, fillval=spval_dbl)
            else
               workr3 = workd3
               call pio_write_darray(File, varid, iodesc3da,&
                                     workr3, status, fillval=spval)
            endif

            call ice_pio_check(status,subname//' ERROR: writing '//avail_hist_fields(n)%vname, &
                               file=__FILE__,line=__LINE__)
         endif
      enddo ! num_avail_hist_fields_3Db
      deallocate(workd3)
      deallocate(workr3)

      ! 3D (fsd)
      allocate(workd3(nx_block,ny_block,nblocks,nfsd_hist))
      allocate(workr3(nx_block,ny_block,nblocks,nfsd_hist))
      do n = n3Dacum+1, n3Dfcum
         nn = n - n3Dacum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            call ice_pio_check(pio_inq_varid(File,avail_hist_fields(n)%vname,varid), &
                 subname//' ERROR: getting varid for '//avail_hist_fields(n)%vname,file=__FILE__,line=__LINE__)
            do j = 1, nblocks
            do i = 1, nfsd_hist
               workd3(:,:,j,i) = a3Df(:,:,i,nn,j)
            enddo
            enddo
            call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)
#ifdef CESM1_PIO
            call pio_setframe(varid, int(1,kind=PIO_OFFSET))
#else
            call pio_setframe(File, varid, int(1,kind=PIO_OFFSET_KIND))
#endif
            call pio_seterrorhandling(File, PIO_RETURN_ERROR)
            if (history_precision == 8) then
               call pio_write_darray(File, varid, iodesc3df,&
                                     workd3, status, fillval=spval_dbl)
            else
               workr3 = workd3
               call pio_write_darray(File, varid, iodesc3df,&
                                     workr3, status, fillval=spval)
            endif

            call ice_pio_check(status,subname//' ERROR: writing '//avail_hist_fields(n)%vname, &
                               file=__FILE__,line=__LINE__)
         endif
      enddo ! num_avail_hist_fields_3Df
      deallocate(workd3)
      deallocate(workr3)

      allocate(workd4(nx_block,ny_block,nblocks,nzilyr,ncat_hist))
      allocate(workr4(nx_block,ny_block,nblocks,nzilyr,ncat_hist))
      ! 4D (categories, vertical ice)
      do n = n3Dfcum+1, n4Dicum
         nn = n - n3Dfcum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            call ice_pio_check(pio_inq_varid(File,avail_hist_fields(n)%vname,varid), &
                 subname//' ERROR: getting varid for '//avail_hist_fields(n)%vname,file=__FILE__,line=__LINE__)
            do j = 1, nblocks
            do i = 1, ncat_hist
            do k = 1, nzilyr
               workd4(:,:,j,k,i) = a4Di(:,:,k,i,nn,j)
            enddo ! k
            enddo ! i
            enddo ! j
            call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)
#ifdef CESM1_PIO
            call pio_setframe(varid, int(1,kind=PIO_OFFSET))
#else
            call pio_setframe(File, varid, int(1,kind=PIO_OFFSET_KIND))
#endif
            call pio_seterrorhandling(File, PIO_RETURN_ERROR)
            if (history_precision == 8) then
               call pio_write_darray(File, varid, iodesc4di,&
                                     workd4, status, fillval=spval_dbl)
            else
               workr4 = workd4
               call pio_write_darray(File, varid, iodesc4di,&
                                     workr4, status, fillval=spval)
            endif
            call ice_pio_check(status,subname//' ERROR: writing '//avail_hist_fields(n)%vname, &
                               file=__FILE__,line=__LINE__)
         endif
      enddo ! num_avail_hist_fields_4Di
      deallocate(workd4)
      deallocate(workr4)

      allocate(workd4(nx_block,ny_block,nblocks,nzslyr,ncat_hist))
      allocate(workr4(nx_block,ny_block,nblocks,nzslyr,ncat_hist))
      ! 4D (categories, vertical snow)
      do n = n4Dicum+1, n4Dscum
         nn = n - n4Dicum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            call ice_pio_check(pio_inq_varid(File,avail_hist_fields(n)%vname,varid), &
                 subname//' ERROR: getting varid for '//avail_hist_fields(n)%vname,file=__FILE__,line=__LINE__)
            do j = 1, nblocks
            do i = 1, ncat_hist
            do k = 1, nzslyr
               workd4(:,:,j,k,i) = a4Ds(:,:,k,i,nn,j)
            enddo ! k
            enddo ! i
            enddo ! j
            call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)
#ifdef CESM1_PIO
            call pio_setframe(varid, int(1,kind=PIO_OFFSET))
#else
            call pio_setframe(File, varid, int(1,kind=PIO_OFFSET_KIND))
#endif
            call pio_seterrorhandling(File, PIO_RETURN_ERROR)
            if (history_precision == 8) then
               call pio_write_darray(File, varid, iodesc4ds,&
                                     workd4, status, fillval=spval_dbl)
            else
               workr4 = workd4
               call pio_write_darray(File, varid, iodesc4ds,&
                                     workr4, status, fillval=spval)
            endif

            call ice_pio_check(status,subname//' ERROR: writing '//avail_hist_fields(n)%vname, &
                               file=__FILE__,line=__LINE__)
         endif
      enddo ! num_avail_hist_fields_4Ds
      deallocate(workd4)
      deallocate(workr4)

      allocate(workd4(nx_block,ny_block,nblocks,nfsd_hist,ncat_hist))
      allocate(workr4(nx_block,ny_block,nblocks,nfsd_hist,ncat_hist))
      ! 4D (categories, fsd)
      do n = n4Dscum+1, n4Dfcum
         nn = n - n4Dscum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            call ice_pio_check(pio_inq_varid(File,avail_hist_fields(n)%vname,varid), &
                 subname//' ERROR: getting varid for '//avail_hist_fields(n)%vname,file=__FILE__,line=__LINE__)
            do j = 1, nblocks
            do i = 1, ncat_hist
            do k = 1, nfsd_hist
               workd4(:,:,j,k,i) = a4Df(:,:,k,i,nn,j)
            enddo ! k
            enddo ! i
            enddo ! j
            call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)
#ifdef CESM1_PIO
            call pio_setframe(varid, int(1,kind=PIO_OFFSET))
#else
            call pio_setframe(File, varid, int(1,kind=PIO_OFFSET_KIND))
#endif
            call pio_seterrorhandling(File, PIO_RETURN_ERROR)
            if (history_precision == 8) then
               call pio_write_darray(File, varid, iodesc4df,&
                                     workd4, status, fillval=spval_dbl)
            else
               workr4 = workd4
               call pio_write_darray(File, varid, iodesc4df,&
                                     workr4, status, fillval=spval)
            endif
            call ice_pio_check(status,subname//' ERROR: writing '//avail_hist_fields(n)%vname, &
                               file=__FILE__,line=__LINE__)
         endif
      enddo ! num_avail_hist_fields_4Df
      deallocate(workd4)
      deallocate(workr4)

!     similarly for num_avail_hist_fields_4Db (define workd4b, iodesc4db)

      !-----------------------------------------------------------------
      ! clean-up PIO descriptors
      !-----------------------------------------------------------------
      call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)

      call pio_freedecomp(File,iodesc2d)
      call pio_freedecomp(File,iodesc3dv)
      call pio_freedecomp(File,iodesc3dc)
      call pio_freedecomp(File,iodesc3di)
      call pio_freedecomp(File,iodesc3db)
      call pio_freedecomp(File,iodesc3da)
      call pio_freedecomp(File,iodesc3df)
      call pio_freedecomp(File,iodesc4di)
      call pio_freedecomp(File,iodesc4ds)
      call pio_freedecomp(File,iodesc4df)

      !-----------------------------------------------------------------
      ! close output dataset
      !-----------------------------------------------------------------

      call pio_closefile(File)
      if (my_task == master_task) then
         write(nu_diag,*) subname,' Finished writing ',trim(ncfile)
      endif

      !-----------------------------------------------------------------
      ! clean up PIO
      !-----------------------------------------------------------------

      call ice_pio_finalize()

      ! reset history parameters
      if (write_histrest_now) then
         history_precision = lhistprec
      endif

      first_call = .false.

      end subroutine ice_write_hist


!=======================================================================
!
! read history restarts, only called for history restarts
!
! author:   T. Craig Nov 2025

      subroutine ice_read_hist

      use ice_kinds_mod
      use ice_blocks, only: nx_block, ny_block
      use ice_broadcast, only: broadcast_scalar
      use ice_calendar, only: nstreams, histfreq
      use ice_communicate, only: my_task, master_task
      use ice_domain, only: nblocks
      use ice_domain_size, only: max_blocks
      use ice_flux, only: albcnt, snwcnt
      use ice_history_shared
      use ice_restart_shared, only: restart_dir
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_readwrite
      use pio

      ! local variables

      real (kind=dbl_kind), dimension(:,:,:),     allocatable :: work2
      real (kind=dbl_kind), dimension(:,:,:,:),   allocatable :: work3
      real (kind=dbl_kind), dimension(:,:,:,:,:), allocatable :: work4

      integer (kind=int_kind) :: i,j,k,n,nn,ns,ncid,status
      logical (kind=log_kind) :: exists
      character (char_len_long) :: ncfile
      character (len=1) :: cns
      character (len=32) :: readstr
      character (len=*), parameter :: readstrT = ' read ok:'
      character (len=*), parameter :: readstrF = ' DID NOT READ:'

      integer (kind=int_kind), parameter :: histprec=8  ! hardwired to double

      type(file_desc_t)     :: File
      type(io_desc_t)       :: iodesc2d, &
                               iodesc3dc, iodesc3di, iodesc3db, iodesc3da, iodesc3df, &
                               iodesc4di, iodesc4ds, iodesc4df
      type(var_desc_t)      :: varid

      logical (kind=log_kind), save :: first_call = .true.

      character(len=*), parameter :: subname = '(ice_read_hist)'

      call ice_timer_start(timer_readwrite)  ! reading/writing
      do ns = 1,nstreams
      if (hist_avg(ns)) then

         write(cns,'(i1.1)') ns

         if (my_task == master_task) then
            call construct_filename(ncfile,'nc',ns, option='histrest')
            ncfile = trim(restart_dir)//ncfile
            write(nu_diag,*) subname,' reading file ',trim(ncfile)
         endif
         call broadcast_scalar(ncfile, master_task)

         ! open file
         inquire(file=trim(ncfile),exist=exists)
         if (exists) then
            File%fh=-1
            call ice_pio_init(mode='read', filename=trim(ncfile), File=File, &
              fformat=trim(history_format), rearr=trim(history_rearranger), &
              iotasks=history_iotasks, root=history_root, stride=history_stride, debug=first_call)

            call pio_seterrorhandling(File, PIO_RETURN_ERROR)

            call ice_pio_initdecomp(iodesc=iodesc2d, precision=histprec)
            call ice_pio_initdecomp(ndim3=ncat_hist, iodesc=iodesc3dc, precision=histprec)
            call ice_pio_initdecomp(ndim3=nzilyr,    iodesc=iodesc3di, precision=histprec)
            call ice_pio_initdecomp(ndim3=nzblyr,    iodesc=iodesc3db, precision=histprec)
            call ice_pio_initdecomp(ndim3=nzalyr,    iodesc=iodesc3da, precision=histprec)
            call ice_pio_initdecomp(ndim3=nfsd_hist, iodesc=iodesc3df, precision=histprec)
            call ice_pio_initdecomp(ndim3=nzilyr,    ndim4=ncat_hist, iodesc=iodesc4di, precision=histprec)
            call ice_pio_initdecomp(ndim3=nzslyr,    ndim4=ncat_hist, iodesc=iodesc4ds, precision=histprec)
            call ice_pio_initdecomp(ndim3=nfsd_hist, ndim4=ncat_hist, iodesc=iodesc4df, precision=histprec)

            !-----------------------------------------------------------------
            ! read variable data
            !-----------------------------------------------------------------

            readstr = readstrF
            status = pio_inq_varid(File, 'time_beg', varid)
            if (status == PIO_NOERR) status  = pio_get_var(File,varid,(/1/),time_beg(ns))
            if (status == PIO_NOERR) readstr = readstrT
            if (my_task == master_task) then
               write(nu_diag,*) subname,trim(readstr),' time_beg'
            endif

            readstr = readstrF
            status = pio_inq_varid(File, 'avgct', varid)
            if (status == PIO_NOERR) status  = pio_get_var(File,varid,(/1/),avgct(ns))
            if (status == PIO_NOERR) readstr = readstrT
            if (my_task == master_task) then
               write(nu_diag,*) subname,trim(readstr),' avgct'
            endif

            allocate(work2(nx_block,ny_block,max_blocks))

            work2(:,:,:) = c0
            readstr = readstrF
            status = pio_inq_varid(File, 'albcnt'//cns, varid)
            if (status == PIO_NOERR) call pio_read_darray(File, varid, iodesc2d, work2, status)
            if (status == PIO_NOERR) then
               readstr = readstrT
               albcnt(:,:,:,ns) = work2(:,:,:)
            endif
            if (my_task == master_task) then
               write(nu_diag,*) subname,trim(readstr),' albcnt'//cns
            endif

            work2(:,:,:) = c0
            readstr = readstrF
            status = pio_inq_varid(File, 'snwcnt'//cns, varid)
            if (status == PIO_NOERR) call pio_read_darray(File, varid, iodesc2d, work2, status)
            if (status == PIO_NOERR) then
               readstr = readstrT
               snwcnt(:,:,:,ns) = work2(:,:,:)
            endif
            if (my_task == master_task) then
               write(nu_diag,*) subname,trim(readstr),' snwcnt'//cns
            endif

            do n=1,num_avail_hist_fields_2D
               if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) then
                  readstr = readstrF
                  status = pio_inq_varid(File,avail_hist_fields(n)%vname,varid)
                  work2(:,:,:) = c0
                  if (status == PIO_NOERR) call pio_read_darray(File, varid, iodesc2d, work2, status)
                  if (status == PIO_NOERR) then
                     readstr = readstrT
                     a2D(:,:,n,:) = work2(:,:,:)
                  endif
                  if (my_task == master_task) then
                     write(nu_diag,*) subname,trim(readstr),trim(avail_hist_fields(n)%vname)
                  endif
               endif
            enddo

            deallocate(work2)
            allocate(work3(nx_block,ny_block,max_blocks,ncat_hist))

            ! 2D
            do n = n2D + 1, n3Dccum
               nn = n - n2D
               if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) then
                  readstr = readstrF
                  status = pio_inq_varid(File,avail_hist_fields(n)%vname,varid)
                  work3(:,:,:,:) = c0
                  if (status == PIO_NOERR) call pio_read_darray(File, varid, iodesc3dc, work3, status)
                  if (status == PIO_NOERR) then
                     readstr = readstrT
                     do j = 1, nblocks
                     do i = 1, ncat_hist
                        a3Dc(:,:,i,nn,j) = work3(:,:,j,i)
                     enddo
                     enddo
                  endif
                  if (my_task == master_task) then
                     write(nu_diag,*) subname,trim(readstr),trim(avail_hist_fields(n)%vname)
                  endif
               endif
            enddo

            deallocate(work3)
            allocate(work3(nx_block,ny_block,max_blocks,ncat_hist))

            ! 3D (category)
            do n = n2D + 1, n3Dccum
               nn = n - n2D
               if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) then
                  readstr = readstrF
                  status = pio_inq_varid(File,avail_hist_fields(n)%vname,varid)
                  work3(:,:,:,:) = c0
                  if (status == PIO_NOERR) call pio_read_darray(File, varid, iodesc3dc, work3, status)
                  if (status == PIO_NOERR) then
                     readstr = readstrT
                     do j = 1, nblocks
                     do i = 1, ncat_hist
                        a3Dc(:,:,i,nn,j) = work3(:,:,j,i)
                     enddo
                     enddo
                  endif
                  if (my_task == master_task) then
                     write(nu_diag,*) subname,trim(readstr),trim(avail_hist_fields(n)%vname)
                  endif
               endif
            enddo

            deallocate(work3)
            allocate(work3(nx_block,ny_block,max_blocks,nzilyr))

            ! 3D (vertical ice)
            do n = n3Dccum+1, n3Dzcum
               nn = n - n3Dccum
               if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) then
                  readstr = readstrF
                  status = pio_inq_varid(File,avail_hist_fields(n)%vname,varid)
                  work3(:,:,:,:) = c0
                  if (status == PIO_NOERR) call pio_read_darray(File, varid, iodesc3di, work3, status)
                  if (status == PIO_NOERR) then
                     readstr = readstrT
                     do j = 1, nblocks
                     do i = 1, nzilyr
                        a3Dz(:,:,i,nn,j) = work3(:,:,j,i)
                     enddo
                     enddo
                  endif
                  if (my_task == master_task) then
                     write(nu_diag,*) subname,trim(readstr),trim(avail_hist_fields(n)%vname)
                  endif
               endif
            enddo

            deallocate(work3)
            allocate(work3(nx_block,ny_block,max_blocks,nzblyr))

            ! 3D (vertical ice biology)
            do n = n3Dzcum+1, n3Dbcum
               nn = n - n3Dzcum
               if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) then
                  readstr = readstrF
                  status = pio_inq_varid(File,avail_hist_fields(n)%vname,varid)
                  work3(:,:,:,:) = c0
                  if (status == PIO_NOERR) call pio_read_darray(File, varid, iodesc3db, work3, status)
                  if (status == PIO_NOERR) then
                     readstr = readstrT
                     do j = 1, nblocks
                     do i = 1, nzblyr
                        a3Db(:,:,i,nn,j) = work3(:,:,j,i)
                     enddo
                     enddo
                  endif
                  if (my_task == master_task) then
                     write(nu_diag,*) subname,trim(readstr),trim(avail_hist_fields(n)%vname)
                  endif
               endif
            enddo

            deallocate(work3)
            allocate(work3(nx_block,ny_block,max_blocks,nzalyr))

            ! 3D (vertical snow biology)
            do n = n3Dbcum+1, n3Dacum
               nn = n - n3Dbcum
               if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) then
                  readstr = readstrF
                  status = pio_inq_varid(File,avail_hist_fields(n)%vname,varid)
                  work3(:,:,:,:) = c0
                  if (status == PIO_NOERR) call pio_read_darray(File, varid, iodesc3da, work3, status)
                  if (status == PIO_NOERR) then
                     readstr = readstrT
                     do j = 1, nblocks
                     do i = 1, nzalyr
                        a3Da(:,:,i,nn,j) = work3(:,:,j,i)
                     enddo
                     enddo
                  endif
                  if (my_task == master_task) then
                     write(nu_diag,*) subname,trim(readstr),trim(avail_hist_fields(n)%vname)
                  endif
               endif
            enddo

            deallocate(work3)
            allocate(work3(nx_block,ny_block,max_blocks,nfsd_hist))

            ! 3D (fsd)
            do n = n3Dacum+1, n3Dfcum
               nn = n - n3Dacum
               if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) then
                  readstr = readstrF
                  status = pio_inq_varid(File,avail_hist_fields(n)%vname,varid)
                  work3(:,:,:,:) = c0
                  if (status == PIO_NOERR) call pio_read_darray(File, varid, iodesc3df, work3, status)
                  if (status == PIO_NOERR) then
                     readstr = readstrT
                     do j = 1, nblocks
                     do i = 1, nfsd_hist
                        a3Df(:,:,i,nn,j) = work3(:,:,j,i)
                     enddo
                     enddo
                  endif
                  if (my_task == master_task) then
                     write(nu_diag,*) subname,trim(readstr),trim(avail_hist_fields(n)%vname)
                  endif
               endif
            enddo

            deallocate(work3)
            allocate(work4(nx_block,ny_block,max_blocks,nzilyr,ncat_hist))

            ! 4D (categories, vertical ice)
            do n = n3Dfcum+1, n4Dicum
               nn = n - n3Dfcum
               if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) then
                  readstr = readstrF
                  status = pio_inq_varid(File,avail_hist_fields(n)%vname,varid)
                  work4(:,:,:,:,:) = c0
                  if (status == PIO_NOERR) call pio_read_darray(File, varid, iodesc4di, work4, status)
                  if (status == PIO_NOERR) then
                     readstr = readstrT
                     do j = 1, nblocks
                     do i = 1, ncat_hist
                     do k = 1, nzilyr
                        a4Di(:,:,k,i,nn,j) = work4(:,:,j,k,i)
                     enddo
                     enddo
                     enddo
                  endif
                  if (my_task == master_task) then
                     write(nu_diag,*) subname,trim(readstr),trim(avail_hist_fields(n)%vname)
                  endif
               endif
            enddo

            deallocate(work4)
            allocate(work4(nx_block,ny_block,max_blocks,nzslyr,ncat_hist))

            ! 4D (categories, vertical snow)
            do n = n4Dicum+1, n4Dscum
               nn = n - n4Dicum
               if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) then
                  readstr = readstrF
                  status = pio_inq_varid(File,avail_hist_fields(n)%vname,varid)
                  work4(:,:,:,:,:) = c0
                  if (status == PIO_NOERR) call pio_read_darray(File, varid, iodesc4ds, work4, status)
                  if (status == PIO_NOERR) then
                     readstr = readstrT
                     do j = 1, nblocks
                     do i = 1, ncat_hist
                     do k = 1, nzslyr
                        a4Ds(:,:,k,i,nn,j) = work4(:,:,j,k,i)
                     enddo
                     enddo
                     enddo
                  endif
                  if (my_task == master_task) then
                     write(nu_diag,*) subname,trim(readstr),trim(avail_hist_fields(n)%vname)
                  endif
               endif
            enddo

            deallocate(work4)
            allocate(work4(nx_block,ny_block,max_blocks,nfsd_hist,ncat_hist))

            ! 4D (categories, fsd)
            do n = n4Dscum+1, n4Dfcum
               nn = n - n4Dscum
               if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) then
                  readstr = readstrF
                  status = pio_inq_varid(File,avail_hist_fields(n)%vname,varid)
                  work4(:,:,:,:,:) = c0
                  if (status == PIO_NOERR) call pio_read_darray(File, varid, iodesc4df, work4, status)
                  if (status == PIO_NOERR) then
                     readstr = readstrT
                     do j = 1, nblocks
                     do i = 1, ncat_hist
                     do k = 1, nfsd_hist
                        a4Df(:,:,k,i,nn,j) = work4(:,:,j,k,i)
                     enddo
                     enddo
                     enddo
                  endif
                  if (my_task == master_task) then
                     write(nu_diag,*) subname,trim(readstr),trim(avail_hist_fields(n)%vname)
                  endif
               endif
            enddo

            deallocate(work4)

            !-----------------------------------------------------------------
            ! clean-up PIO descriptors
            !-----------------------------------------------------------------
            call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)

            call pio_freedecomp(File,iodesc2d)
            call pio_freedecomp(File,iodesc3dc)
            call pio_freedecomp(File,iodesc3di)
            call pio_freedecomp(File,iodesc3db)
            call pio_freedecomp(File,iodesc3da)
            call pio_freedecomp(File,iodesc3df)
            call pio_freedecomp(File,iodesc4di)
            call pio_freedecomp(File,iodesc4ds)
            call pio_freedecomp(File,iodesc4df)

            !-----------------------------------------------------------------
            ! close output dataset
            !-----------------------------------------------------------------

            call pio_closefile(File)
            if (my_task == master_task) then
               write(nu_diag,*) subname,' Finished reading ',trim(ncfile)
            endif

            !-----------------------------------------------------------------
            ! clean up PIO
            !-----------------------------------------------------------------

            call ice_pio_finalize()

            first_call = .false.

         endif ! open file success
      endif ! hist_avg
      enddo ! nstreams
      call ice_timer_stop(timer_readwrite)  ! reading/writing

      end subroutine ice_read_hist

!=======================================================================
! Defines a coordinate var in the history file
! coordinates have short_name, long_name and units attributes,
!  and are compressed for 'hdf5' when more than one dimensional

      subroutine ice_hist_coord_def(File, coord,lprecision, dimids,varid)

      use pio, only: file_desc_t, var_desc_t, pio_def_var,  pio_put_att
#ifndef USE_PIO1
      use pio, only: pio_def_var_deflate
      use pio_nf, only: pio_def_var_chunking !This is missing from pio module <2.6.0
      use netcdf, only: NF90_CHUNKED
      use ice_history_shared, only: history_deflate, history_chunksize, history_format
#endif

      type(file_desc_t),   intent(inout) :: File
      type(coord_attributes), intent(in) :: coord
      integer(kind=int_kind), intent(in) :: dimids(:), lprecision
      type(var_desc_t),    intent(inout) :: varid

      ! local vars
      integer(kind=int_kind) :: chunks(size(dimids)), i, status

      character(len=*), parameter :: subname = '(ice_hist_coord_def)'

      !define var, set deflate, long_name and units
      status = pio_def_var(File, coord%short_name, lprecision, dimids, varid)
      call ice_pio_check(status, &
         subname//' ERROR: defining coord '//coord%short_name,file=__FILE__,line=__LINE__)
#ifndef USE_PIO1
      if (history_deflate/=0 .and. history_format=='hdf5') then
         status = pio_def_var_deflate(File, varid, shuffle=0, deflate=1, deflate_level=history_deflate)
         call ice_pio_check(status, &
            subname//' ERROR: deflating coord '//coord%short_name,file=__FILE__,line=__LINE__)
      endif

      if (history_format=='hdf5' .and. size(dimids)>1) then
         if (dimids(1)==imtid .and. dimids(2)==jmtid) then
            chunks(1)=history_chunksize(1)
            chunks(2)=history_chunksize(2)
            do i = 3, size(dimids)
               chunks(i) = 0
            enddo
            status = pio_def_var_chunking(File, varid, NF90_CHUNKED, chunks)
            call ice_pio_check(status, &
               subname//' ERROR: chunking coord '//coord%short_name,file=__FILE__,line=__LINE__)
         endif
      endif
#endif
      if (coord%long_name(1:3) /= 'und') then
         call ice_pio_check(pio_put_att(File,varid,'long_name',trim(coord%long_name)), &
               subname//' ERROR: defining att long_name '//coord%long_name,file=__FILE__,line=__LINE__)
      endif
      if (coord%units(1:3) /= 'und') then
         call ice_pio_check(pio_put_att(File, varid, 'units', trim(coord%units)), &
               subname//' ERROR: defining att units '//coord%units,file=__FILE__,line=__LINE__)
      endif
      if (coord%axis(1:3) /= 'und') then
         call ice_pio_check(pio_put_att(File, varid, 'axis', trim(coord%axis)), &
               subname//' ERROR: defining att axis '//coord%units,file=__FILE__,line=__LINE__)
      endif

      end subroutine ice_hist_coord_def

!=======================================================================
! Defines a (time-dependent) history var in the history file
! variables have short_name, long_name and units, coordinates and cell_measures attributes,
!  and are compressed and chunked for 'hdf5'

      subroutine ice_hist_field_def(File, hfield,lprecision, dimids, ns)

      use pio, only: file_desc_t , var_desc_t, pio_def_var, pio_put_att
#ifndef USE_PIO1
      use pio, only: pio_def_var_deflate
      use pio_nf, only: pio_def_var_chunking !This is missing from pio module <2.6.0
      use netcdf, only: NF90_CHUNKED
      use ice_history_shared, only: history_deflate, history_chunksize, history_format
#endif
      use ice_history_shared, only: ice_hist_field, history_precision, hist_avg
      use ice_calendar, only: histfreq, histfreq_n, write_ic

      type(file_desc_t),   intent(inout) :: File
      type(ice_hist_field)  , intent(in) :: hfield
      integer(kind=int_kind), intent(in) :: dimids(:), lprecision, ns

      ! local vars
      type(var_desc_t) :: varid
      integer(kind=int_kind) :: chunks(size(dimids)), i, status

      character(len=*), parameter :: subname = '(ice_hist_field_def)'

      status = pio_def_var(File, hfield%vname, lprecision, dimids, varid)
      call ice_pio_check(status, &
         subname//' ERROR: defining var '//hfield%vname,file=__FILE__,line=__LINE__)

#ifndef USE_PIO1
      if (history_deflate/=0 .and. history_format=='hdf5') then
         status = pio_def_var_deflate(File, varid, shuffle=0, deflate=1, deflate_level=history_deflate)
         call ice_pio_check(status, &
            subname//' ERROR: deflating var '//hfield%vname,file=__FILE__,line=__LINE__)
      endif

      if (history_format=='hdf5' .and. size(dimids)>1) then
         if (dimids(1)==imtid .and. dimids(2)==jmtid) then
            chunks(1)=history_chunksize(1)
            chunks(2)=history_chunksize(2)
            do i = 3, size(dimids)
               chunks(i) = 0
            enddo
            status = pio_def_var_chunking(File, varid, NF90_CHUNKED, chunks)
            call ice_pio_check(status, subname//' ERROR: chunking var '//hfield%vname,file=__FILE__,line=__LINE__)
         endif
      endif
#endif

      !var attributes

      call ice_pio_check(pio_put_att(File,varid,'units', trim(hfield%vunit)), &
           subname//' ERROR: defining att units '//trim(hfield%vunit),file=__FILE__,line=__LINE__)

      call ice_pio_check(pio_put_att(File,varid, 'long_name', trim(hfield%vdesc)), &
           subname//' ERROR: defining att long_name '//trim(hfield%vdesc),file=__FILE__,line=__LINE__)

      call ice_pio_check(pio_put_att(File,varid,'coordinates', trim(hfield%vcoord)), &
           subname//' ERROR: defining att coordinates '//trim(hfield%vdesc),file=__FILE__,line=__LINE__)

      call ice_pio_check(pio_put_att(File,varid,'cell_measures',trim(hfield%vcellmeas)), &
           subname//' ERROR: defining att cell_measures '//trim(hfield%vcoord),file=__FILE__,line=__LINE__)

      if (hfield%vcomment /= "none") then
         call ice_pio_check(pio_put_att(File,varid,'comment', trim(hfield%vcomment)), &
              subname//' ERROR: defining att comment '//trim(hfield%vcomment),file=__FILE__,line=__LINE__)
      endif

      call ice_write_hist_fill(File,varid,hfield%vname,history_precision)

      ! Add cell_methods attribute to variables if averaged
      if (hist_avg(ns) .and. .not. write_ic) then
         if    (TRIM(hfield%vname(1:4))/='sig1' &
           .and.TRIM(hfield%vname(1:4))/='sig2' &
           .and.TRIM(hfield%vname(1:5))/='trsig' &
           .and.TRIM(hfield%vname(1:4))/='divu' &
           .and.TRIM(hfield%vname(1:5))/='shear' &
           .and.TRIM(hfield%vname(1:4))/='vort' &
           .and.TRIM(hfield%vname(1:9))/='frz_onset' &
           .and.TRIM(hfield%vname(1:9))/='mlt_onset' &
           .and.TRIM(hfield%vname(1:6))/='aisnap' &
           .and.TRIM(hfield%vname(1:6))/='hisnap' &
           .and.TRIM(hfield%vname(1:8))/='sidivvel' &
           .and.TRIM(hfield%vname(1:10))/='sishearvel' &
           .and.TRIM(hfield%vname(1:11))/='sistressave' &
           .and.TRIM(hfield%vname(1:11))/='sistressmax' &
           .and.TRIM(hfield%vname(1:4))/='sigP') then
            if (trim(hfield%avg_ice_present) /= 'none') then
               call ice_pio_check(pio_put_att(File,varid,'cell_methods', &
                    'area: time: mean where sea ice (mask=siconc)'), &
                    subname//' ERROR: defining att cell_methods',file=__FILE__,line=__LINE__)
            else
               call ice_pio_check(pio_put_att(File,varid,'cell_methods','time: mean'), &
                    subname//' ERROR: defining att cell_methods',file=__FILE__,line=__LINE__)
            endif
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
          .or.TRIM(hfield%vname(1:8))=='sidivvel' &
          .or.TRIM(hfield%vname(1:10))=='sishearvel' &
          .or.TRIM(hfield%vname(1:11))=='sistressave' &
          .or.TRIM(hfield%vname(1:11))=='sistressmax' &
          .or.TRIM(hfield%vname(1:9))=='mlt_onset' &
          .or.TRIM(hfield%vname(1:9))=='frz_onset' &
          .or.TRIM(hfield%vname(1:6))=='hisnap' &
          .or.TRIM(hfield%vname(1:6))=='aisnap') then
         call ice_pio_check(pio_put_att(File,varid,'time_rep','instantaneous'), &
              subname//' ERROR: defining att time_rep i',file=__FILE__,line=__LINE__)
      else
         call ice_pio_check(pio_put_att(File,varid,'time_rep','averaged'), &
              subname//' ERROR: defining att time_rep a',file=__FILE__,line=__LINE__)
      endif

      end subroutine ice_hist_field_def

!=======================================================================
! Defines missing_value and _FillValue attributes

      subroutine ice_write_hist_fill(File,varid,vname,precision)

      use pio, only: pio_put_att, file_desc_t, var_desc_t

      type(file_desc_t),    intent(inout) :: File
      type(var_desc_t),        intent(in) :: varid
      character(len=*),        intent(in) :: vname
      integer (kind=int_kind), intent(in) :: precision

      ! local variables

      integer (kind=int_kind) :: status
      character(len=*), parameter :: subname = '(ice_write_hist_fill)'

      if (precision == 8) then
         call ice_pio_check(pio_put_att(File, varid, 'missing_value', spval_dbl), &
              subname//' ERROR: defining att missing_value',file=__FILE__,line=__LINE__)
         call ice_pio_check(pio_put_att(File, varid,'_FillValue',spval_dbl), &
              subname//' ERROR: defining att _FillValue',file=__FILE__,line=__LINE__)
      else
         call ice_pio_check(pio_put_att(File, varid, 'missing_value', spval), &
              subname//' ERROR: defining att missing_value',file=__FILE__,line=__LINE__)
         call ice_pio_check(pio_put_att(File, varid,'_FillValue',spval), &
              subname//' ERROR: defining att _FillValue',file=__FILE__,line=__LINE__)
      endif

      end subroutine ice_write_hist_fill

!=======================================================================

      end module ice_history_write

!=======================================================================
