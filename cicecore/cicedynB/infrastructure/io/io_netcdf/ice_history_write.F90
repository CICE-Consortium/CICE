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

      use ice_constants, only: c0, c360, spval
      use ice_fileunits, only: nu_diag
      use ice_exit, only: abort_ice
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_parameters

      implicit none
      private
      public :: ice_write_hist
      
!=======================================================================

      contains

!=======================================================================
!
! write average ice quantities or snapshots
!
! author:   Elizabeth C. Hunke, LANL

      subroutine ice_write_hist (ns)

      use ice_kinds_mod
#ifdef ncdf
      use ice_arrays_column, only: hin_max
      use ice_blocks, only: nx_block, ny_block
      use ice_broadcast, only: broadcast_scalar
      use ice_calendar, only: time, sec, idate, idate0, write_ic, &
          histfreq, dayyr, days_per_year, use_leap_years
      use ice_communicate, only: my_task, master_task
      use ice_domain, only: distrb_info
      use ice_domain_size, only: nx_global, ny_global, max_nstrm, max_blocks
      use ice_gather_scatter, only: gather_global
      use ice_grid, only: TLON, TLAT, ULON, ULAT, hm, bm, tarea, uarea, &
          dxu, dxt, dyu, dyt, HTN, HTE, ANGLE, ANGLET, &
          lont_bounds, latt_bounds, lonu_bounds, latu_bounds
      use ice_history_shared
      use ice_restart_shared, only: runid, lcdf64
      use netcdf
#endif

      integer (kind=int_kind), intent(in) :: ns

      ! local variables

#ifdef ncdf
      real (kind=dbl_kind),  dimension(:,:),   allocatable :: work_g1
      real (kind=real_kind), dimension(:,:),   allocatable :: work_gr
      real (kind=real_kind), dimension(:,:,:), allocatable :: work_gr3
      real (kind=dbl_kind),  dimension(nx_block,ny_block,max_blocks) :: &
         work1

      integer (kind=int_kind) :: i,k,ic,n,nn, &
         ncid,status,imtid,jmtid,kmtidi,kmtids,kmtidb, cmtid,timid,varid, &
         nvertexid,ivertex,kmtida,iflag
      integer (kind=int_kind), dimension(3) :: dimid
      integer (kind=int_kind), dimension(4) :: dimidz
      integer (kind=int_kind), dimension(5) :: dimidcz
      integer (kind=int_kind), dimension(3) :: dimid_nverts
      integer (kind=int_kind), dimension(5) :: dimidex
!     real (kind=real_kind) :: ltime
      real (kind=dbl_kind)  :: ltime2
      character (char_len) :: title
      character (char_len_long) :: ncfile(max_nstrm)
      real (kind=dbl_kind)  :: secday, rad_to_deg

      integer (kind=int_kind) :: ind,boundid

      character (char_len) :: start_time,current_date,current_time
      character (len=8) :: cdate

      ! 4 coordinate variables: TLON, TLAT, ULON, ULAT
      INTEGER (kind=int_kind), PARAMETER :: ncoord = 4

      ! 4 vertices in each grid cell
      INTEGER (kind=int_kind), PARAMETER :: nverts = 4

      ! 4 variables describe T, U grid boundaries:
      ! lont_bounds, latt_bounds, lonu_bounds, latu_bounds
      INTEGER (kind=int_kind), PARAMETER :: nvar_verts = 4

      TYPE coord_attributes         ! netcdf coordinate attributes
        character (len=11)   :: short_name
        character (len=45)   :: long_name
        character (len=20)   :: units
      END TYPE coord_attributes

      TYPE req_attributes         ! req'd netcdf attributes
        type (coord_attributes) :: req
        character (len=20)   :: coordinates
      END TYPE req_attributes

      TYPE(req_attributes), dimension(nvar) :: var
      TYPE(coord_attributes), dimension(ncoord) :: coord_var
      TYPE(coord_attributes), dimension(nvar_verts) :: var_nverts
      TYPE(coord_attributes), dimension(nvarz) :: var_nz
      CHARACTER (char_len), dimension(ncoord) :: coord_bounds

      character(len=*), parameter :: subname = '(ice_write_hist)'

      call icepack_query_parameters(secday_out=secday, rad_to_deg_out=rad_to_deg)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__, line=__LINE__)

      if (my_task == master_task) then

!       ltime=time/int(secday)
        ltime2=time/int(secday)

        call construct_filename(ncfile(ns),'nc',ns)

        ! add local directory path name to ncfile
        if (write_ic) then
          ncfile(ns) = trim(incond_dir)//ncfile(ns)
        else
          ncfile(ns) = trim(history_dir)//ncfile(ns)
        endif

        ! create file
        iflag = nf90_clobber
        if (lcdf64) iflag = ior(iflag,nf90_64bit_offset)
        status = nf90_create(ncfile(ns), iflag, ncid)
        if (status /= nf90_noerr) call abort_ice(subname// &
           'ERROR: creating history ncfile '//ncfile(ns))

      !-----------------------------------------------------------------
      ! define dimensions
      !-----------------------------------------------------------------

        if (hist_avg) then
          status = nf90_def_dim(ncid,'d2',2,boundid)
          if (status /= nf90_noerr) call abort_ice(subname// &
             'ERROR: defining dim d2')
        endif

        status = nf90_def_dim(ncid,'ni',nx_global,imtid)
        if (status /= nf90_noerr) call abort_ice(subname// &
                      'ERROR: defining dim ni')

        status = nf90_def_dim(ncid,'nj',ny_global,jmtid)
        if (status /= nf90_noerr) call abort_ice(subname// &
                      'ERROR: defining dim nj')

        status = nf90_def_dim(ncid,'nc',ncat_hist,cmtid)
        if (status /= nf90_noerr) call abort_ice(subname// &
                      'ERROR: defining dim nc')

        status = nf90_def_dim(ncid,'nkice',nzilyr,kmtidi)
        if (status /= nf90_noerr) call abort_ice(subname// &
                      'ERROR: defining dim nki')

        status = nf90_def_dim(ncid,'nksnow',nzslyr,kmtids)
        if (status /= nf90_noerr) call abort_ice(subname// &
                      'ERROR: defining dim nks')

        status = nf90_def_dim(ncid,'nkbio',nzblyr,kmtidb)
        if (status /= nf90_noerr) call abort_ice(subname// &
                      'ERROR: defining dim nkb')

        status = nf90_def_dim(ncid,'nkaer',nzalyr,kmtida)
        if (status /= nf90_noerr) call abort_ice(subname// &
                      'ERROR: defining dim nka')

        status = nf90_def_dim(ncid,'time',NF90_UNLIMITED,timid)
        if (status /= nf90_noerr) call abort_ice(subname// &
                      'ERROR: defining dim time')

        status = nf90_def_dim(ncid,'nvertices',nverts,nvertexid)
        if (status /= nf90_noerr) call abort_ice(subname// &
                      'ERROR: defining dim nverts')

      !-----------------------------------------------------------------
      ! define coordinate variables
      !-----------------------------------------------------------------

!sgl        status = nf90_def_var(ncid,'time',nf90_float,timid,varid)
        status = nf90_def_var(ncid,'time',nf90_double,timid,varid)
        if (status /= nf90_noerr) call abort_ice(subname// &
                      'ERROR: defining var time')

        status = nf90_put_att(ncid,varid,'long_name','model time')
        if (status /= nf90_noerr) call abort_ice(subname// &
                      'ice Error: time long_name')

        write(cdate,'(i8.8)') idate0
        write(title,'(a,a,a,a,a,a,a,a)') 'days since ', &
              cdate(1:4),'-',cdate(5:6),'-',cdate(7:8),' 00:00:00'
        status = nf90_put_att(ncid,varid,'units',title)
        if (status /= nf90_noerr) call abort_ice(subname// &
                      'ERROR: time units')

        if (days_per_year == 360) then
           status = nf90_put_att(ncid,varid,'calendar','360_day')
           if (status /= nf90_noerr) call abort_ice(subname// &
                         'ERROR: time calendar')
        elseif (days_per_year == 365 .and. .not.use_leap_years ) then
           status = nf90_put_att(ncid,varid,'calendar','NoLeap')
           if (status /= nf90_noerr) call abort_ice(subname// &
                         'ERROR: time calendar')
        elseif (use_leap_years) then
           status = nf90_put_att(ncid,varid,'calendar','Gregorian')
           if (status /= nf90_noerr) call abort_ice(subname// &
                         'ERROR: time calendar')
        else
           call abort_ice(subname//'ERROR: invalid calendar settings')
        endif

        if (hist_avg) then
          status = nf90_put_att(ncid,varid,'bounds','time_bounds')
          if (status /= nf90_noerr) call abort_ice(subname// &
                      'ERROR: time bounds')
        endif

      !-----------------------------------------------------------------
      ! Define attributes for time bounds if hist_avg is true
      !-----------------------------------------------------------------

        if (hist_avg) then
          dimid(1) = boundid
          dimid(2) = timid
          status = nf90_def_var(ncid,'time_bounds',nf90_float,dimid(1:2),varid)
          if (status /= nf90_noerr) call abort_ice(subname// &
                        'ERROR: defining var time_bounds')
          status = nf90_put_att(ncid,varid,'long_name', &
                                'boundaries for time-averaging interval')
          if (status /= nf90_noerr) call abort_ice(subname// &
                        'ERROR: time_bounds long_name')
          write(cdate,'(i8.8)') idate0
          write(title,'(a,a,a,a,a,a,a,a)') 'days since ', &
                cdate(1:4),'-',cdate(5:6),'-',cdate(7:8),' 00:00:00'
          status = nf90_put_att(ncid,varid,'units',title)
          if (status /= nf90_noerr) call abort_ice(subname// &
                        'ERROR: time_bounds units')
        endif

      !-----------------------------------------------------------------
      ! define information for required time-invariant variables
      !-----------------------------------------------------------------

      ind = 0
      ind = ind + 1
      coord_var(ind) = coord_attributes('TLON', &
                       'T grid center longitude', 'degrees_east')
      coord_bounds(ind) = 'lont_bounds'
      ind = ind + 1
      coord_var(ind) = coord_attributes('TLAT', &
                       'T grid center latitude',  'degrees_north')
      coord_bounds(ind) = 'latt_bounds'
      ind = ind + 1
      coord_var(ind) = coord_attributes('ULON', &
                       'U grid center longitude', 'degrees_east')
      coord_bounds(ind) = 'lonu_bounds'
      ind = ind + 1
      coord_var(ind) = coord_attributes('ULAT', &
                       'U grid center latitude',  'degrees_north')
      coord_bounds(ind) = 'latu_bounds'

      var_nz(1) = coord_attributes('NCAT', 'category maximum thickness', 'm')
      var_nz(2) = coord_attributes('VGRDi', 'vertical ice levels', '1')
      var_nz(3) = coord_attributes('VGRDs', 'vertical snow levels', '1')
      var_nz(4) = coord_attributes('VGRDb', 'vertical ice-bio levels', '1')
      var_nz(5) = coord_attributes('VGRDa', 'vertical snow-ice-bio levels', '1')

      !-----------------------------------------------------------------
      ! define information for optional time-invariant variables
      !-----------------------------------------------------------------

      var(n_tarea)%req = coord_attributes('tarea', &
                  'area of T grid cells', 'm^2')
      var(n_tarea)%coordinates = 'TLON TLAT'
      var(n_uarea)%req = coord_attributes('uarea', &
                  'area of U grid cells', 'm^2')
      var(n_uarea)%coordinates = 'ULON ULAT'
      var(n_dxt)%req = coord_attributes('dxt', &
                  'T cell width through middle', 'm')
      var(n_dxt)%coordinates = 'TLON TLAT'
      var(n_dyt)%req = coord_attributes('dyt', &
                  'T cell height through middle', 'm')
      var(n_dyt)%coordinates = 'TLON TLAT'
      var(n_dxu)%req = coord_attributes('dxu', &
                  'U cell width through middle', 'm')
      var(n_dxu)%coordinates = 'ULON ULAT'
      var(n_dyu)%req = coord_attributes('dyu', &
                  'U cell height through middle', 'm')
      var(n_dyu)%coordinates = 'ULON ULAT'
      var(n_HTN)%req = coord_attributes('HTN', &
                  'T cell width on North side','m')
      var(n_HTN)%coordinates = 'TLON TLAT'
      var(n_HTE)%req = coord_attributes('HTE', &
                  'T cell width on East side', 'm')
      var(n_HTE)%coordinates = 'TLON TLAT'
      var(n_ANGLE)%req = coord_attributes('ANGLE', &
                  'angle grid makes with latitude line on U grid', &
                  'radians')
      var(n_ANGLE)%coordinates = 'ULON ULAT'
      var(n_ANGLET)%req = coord_attributes('ANGLET', &
                  'angle grid makes with latitude line on T grid', &
                  'radians')
      var(n_ANGLET)%coordinates = 'TLON TLAT'

      ! These fields are required for CF compliance
      ! dimensions (nx,ny,nverts)
      var_nverts(n_lont_bnds) = coord_attributes('lont_bounds', &
                  'longitude boundaries of T cells', 'degrees_east')
      var_nverts(n_latt_bnds) = coord_attributes('latt_bounds', &
                  'latitude boundaries of T cells', 'degrees_north')
      var_nverts(n_lonu_bnds) = coord_attributes('lonu_bounds', &
                  'longitude boundaries of U cells', 'degrees_east')
      var_nverts(n_latu_bnds) = coord_attributes('latu_bounds', &
                  'latitude boundaries of U cells', 'degrees_north')

      !-----------------------------------------------------------------
      ! define attributes for time-invariant variables
      !-----------------------------------------------------------------

        dimid(1) = imtid
        dimid(2) = jmtid
        dimid(3) = timid

        do i = 1, ncoord
          status = nf90_def_var(ncid, coord_var(i)%short_name, nf90_float, &
                                dimid(1:2), varid)
          if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining short_name for '//coord_var(i)%short_name)
          status = nf90_put_att(ncid,varid,'long_name',coord_var(i)%long_name)
          if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining long_name for '//coord_var(i)%short_name)
          status = nf90_put_att(ncid, varid, 'units', coord_var(i)%units)
          if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining units for '//coord_var(i)%short_name)
          status = nf90_put_att(ncid,varid,'missing_value',spval)
          if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining missing_value for '//coord_var(i)%short_name)
          status = nf90_put_att(ncid,varid,'_FillValue',spval)
          if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining _FillValue for '//coord_var(i)%short_name)
          if (coord_var(i)%short_name == 'ULAT') then
             status = nf90_put_att(ncid,varid,'comment', &
                  'Latitude of NE corner of T grid cell')
             if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining comment for '//coord_var(i)%short_name)
          endif
          if (f_bounds) then
             status = nf90_put_att(ncid, varid, 'bounds', coord_bounds(i))
             if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining bounds for '//coord_var(i)%short_name)
          endif          
        enddo

        ! Extra dimensions (NCAT, NZILYR, NZSLYR, NZBLYR, NZALYR)       
          dimidex(1)=cmtid
          dimidex(2)=kmtidi
          dimidex(3)=kmtids
          dimidex(4)=kmtidb
          dimidex(5)=kmtida
        
        do i = 1, nvarz
           if (igrdz(i)) then
             status = nf90_def_var(ncid, var_nz(i)%short_name, &
                                   nf90_float, dimidex(i), varid)
             if (status /= nf90_noerr) call abort_ice(subname// &
                'ERROR: defining short_name for '//var_nz(i)%short_name)
             status = nf90_put_att(ncid,varid,'long_name',var_nz(i)%long_name)
             if (status /= nf90_noerr) call abort_ice(subname// &
                'ERROR: defining long_name for '//var_nz(i)%short_name)
             status = nf90_put_att(ncid, varid, 'units', var_nz(i)%units)
             if (Status /= nf90_noerr) call abort_ice(subname// &
                'ERROR: defining units for '//var_nz(i)%short_name)
           endif
        enddo

        ! Attributes for tmask, blkmask defined separately, since they have no units
        if (igrd(n_tmask)) then
           status = nf90_def_var(ncid, 'tmask', nf90_float, dimid(1:2), varid)
           if (status /= nf90_noerr) call abort_ice(subname//'ERROR: defining var tmask')
           status = nf90_put_att(ncid,varid, 'long_name', 'ocean grid mask') 
           if (status /= nf90_noerr) call abort_ice(subname//'ERROR: tmask long_name') 
           status = nf90_put_att(ncid, varid, 'coordinates', 'TLON TLAT')
           if (status /= nf90_noerr) call abort_ice(subname//'ERROR: tmask units') 
           status = nf90_put_att(ncid,varid,'comment', '0 = land, 1 = ocean')
           if (status /= nf90_noerr) call abort_ice(subname//'ERROR: tmask comment') 
           status = nf90_put_att(ncid,varid,'missing_value',spval)
           if (status /= nf90_noerr) call abort_ice(subname//'ERROR: defining missing_value for tmask')
           status = nf90_put_att(ncid,varid,'_FillValue',spval)
           if (status /= nf90_noerr) call abort_ice(subname//'ERROR: defining _FillValue for tmask')
        endif

        if (igrd(n_blkmask)) then
           status = nf90_def_var(ncid, 'blkmask', nf90_float, dimid(1:2), varid)
           if (status /= nf90_noerr) call abort_ice(subname//'ERROR: defining var blkmask')
           status = nf90_put_att(ncid,varid, 'long_name', 'ice grid block mask') 
           if (status /= nf90_noerr) call abort_ice(subname//'ERROR: blkmask long_name') 
           status = nf90_put_att(ncid, varid, 'coordinates', 'TLON TLAT')
           if (status /= nf90_noerr) call abort_ice(subname//'ERROR: blkmask units') 
           status = nf90_put_att(ncid,varid,'comment', 'mytask + iblk/100')
           if (status /= nf90_noerr) call abort_ice(subname//'ERROR: blkmask comment') 
           status = nf90_put_att(ncid,varid,'missing_value',spval)
           if (status /= nf90_noerr) call abort_ice(subname//'ERROR: defining missing_value for blkmask')
           status = nf90_put_att(ncid,varid,'_FillValue',spval)
           if (status /= nf90_noerr) call abort_ice(subname//'ERROR: defining _FillValue for blkmask')
        endif

        do i = 3, nvar      ! note n_tmask=1, n_blkmask=2
          if (igrd(i)) then
             status = nf90_def_var(ncid, var(i)%req%short_name, &
                                   nf90_float, dimid(1:2), varid)
             if (status /= nf90_noerr) call abort_ice(subname// &
                  'ERROR: defining variable '//var(i)%req%short_name)
             status = nf90_put_att(ncid,varid, 'long_name', var(i)%req%long_name)
             if (status /= nf90_noerr) call abort_ice(subname// &
                  'ERROR: defining long_name for '//var(i)%req%short_name)
             status = nf90_put_att(ncid, varid, 'units', var(i)%req%units)
             if (status /= nf90_noerr) call abort_ice(subname// &
                  'ERROR: defining units for '//var(i)%req%short_name)
             status = nf90_put_att(ncid, varid, 'coordinates', var(i)%coordinates)
             if (status /= nf90_noerr) call abort_ice(subname// &
                  'ERROR: defining coordinates for '//var(i)%req%short_name)
             status = nf90_put_att(ncid,varid,'missing_value',spval)
             if (status /= nf90_noerr) call abort_ice(subname// &
                  'ERROR: defining missing_value for '//var(i)%req%short_name)
             status = nf90_put_att(ncid,varid,'_FillValue',spval)
             if (status /= nf90_noerr) call abort_ice(subname// &
                  'ERROR: defining _FillValue for '//var(i)%req%short_name)
          endif
        enddo

        ! Fields with dimensions (nverts,nx,ny)
        dimid_nverts(1) = nvertexid
        dimid_nverts(2) = imtid
        dimid_nverts(3) = jmtid
        do i = 1, nvar_verts
          if (f_bounds) then
             status = nf90_def_var(ncid, var_nverts(i)%short_name, &
                                   nf90_float,dimid_nverts, varid)
             if (status /= nf90_noerr) call abort_ice(subname// &
                  'ERROR: defining variable '//var_nverts(i)%short_name)
             status = nf90_put_att(ncid,varid, 'long_name', var_nverts(i)%long_name)
             if (status /= nf90_noerr) call abort_ice(subname// &
                  'ERROR: defining long_name for '//var_nverts(i)%short_name)
             status = nf90_put_att(ncid, varid, 'units', var_nverts(i)%units)
             if (status /= nf90_noerr) call abort_ice(subname// &
                  'ERROR: defining units for '//var_nverts(i)%short_name)
             status = nf90_put_att(ncid,varid,'missing_value',spval)
             if (status /= nf90_noerr) call abort_ice(subname// &
                  'ERROR: defining missing_value for '//var_nverts(i)%short_name)
             status = nf90_put_att(ncid,varid,'_FillValue',spval)
             if (status /= nf90_noerr) call abort_ice(subname// &
                  'ERROR: defining _FillValue for '//var_nverts(i)%short_name)
          endif
        enddo

        do n=1,num_avail_hist_fields_2D
          if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            status  = nf90_def_var(ncid, avail_hist_fields(n)%vname, &
                         nf90_float, dimid, varid)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining variable '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'units', &
                        avail_hist_fields(n)%vunit)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining units for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid, 'long_name', &
                        avail_hist_fields(n)%vdesc)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining long_name for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'coordinates', &
                        avail_hist_fields(n)%vcoord)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining coordinates for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'cell_measures', &
                        avail_hist_fields(n)%vcellmeas)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining cell measures for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'missing_value',spval)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining missing_value for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'_FillValue',spval)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining _FillValue for '//avail_hist_fields(n)%vname)

      !-----------------------------------------------------------------
      ! Add cell_methods attribute to variables if averaged
      !-----------------------------------------------------------------
            if (hist_avg) then
              if (TRIM(avail_hist_fields(n)%vname)/='sig1' &
              .or.TRIM(avail_hist_fields(n)%vname)/='sig2' & 
              .or.TRIM(avail_hist_fields(n)%vname)/='sistreave' & 
              .or.TRIM(avail_hist_fields(n)%vname)/='sistremax' & 
              .or.TRIM(avail_hist_fields(n)%vname)/='sigP') then
                status = nf90_put_att(ncid,varid,'cell_methods','time: mean')
                if (status /= nf90_noerr) call abort_ice(subname// &
                 'ERROR: defining cell methods for '//avail_hist_fields(n)%vname)
              endif
            endif

            if (histfreq(ns) == '1' .or. .not. hist_avg         &
                .or. n==n_divu(ns)      .or. n==n_shear(ns)     &  ! snapshots
                .or. n==n_sig1(ns)      .or. n==n_sig2(ns)      &
                .or. n==n_sigP(ns)      .or. n==n_trsig(ns)     &
                .or. n==n_sistreave(ns) .or. n==n_sistremax(ns) &
                .or. n==n_mlt_onset(ns) .or. n==n_frz_onset(ns) &
                .or. n==n_hisnap(ns)    .or. n==n_aisnap(ns)) then
               status = nf90_put_att(ncid,varid,'time_rep','instantaneous')
            else
               status = nf90_put_att(ncid,varid,'time_rep','averaged')
            endif
          endif
        enddo  ! num_avail_hist_fields_2D

        dimidz(1) = imtid
        dimidz(2) = jmtid
        dimidz(3) = cmtid
        dimidz(4) = timid

        do n = n2D + 1, n3Dccum
          if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            status  = nf90_def_var(ncid, avail_hist_fields(n)%vname, &
                         nf90_float, dimidz, varid)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining variable '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'units', &
                        avail_hist_fields(n)%vunit)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining units for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid, 'long_name', &
                        avail_hist_fields(n)%vdesc)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining long_name for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'coordinates', &
                        avail_hist_fields(n)%vcoord)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining coordinates for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'cell_measures', &
                        avail_hist_fields(n)%vcellmeas)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining cell measures for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'missing_value',spval)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining missing_value for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'_FillValue',spval)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining _FillValue for '//avail_hist_fields(n)%vname)

      !-----------------------------------------------------------------
      ! Add cell_methods attribute to variables if averaged
      !-----------------------------------------------------------------
            if (hist_avg) then
                status = nf90_put_att(ncid,varid,'cell_methods','time: mean')
                if (status /= nf90_noerr) call abort_ice(subname// &
                 'ERROR: defining cell methods for '//avail_hist_fields(n)%vname)
            endif

            if (histfreq(ns) == '1' .or. .not. hist_avg) then
               status = nf90_put_att(ncid,varid,'time_rep','instantaneous')
            else
               status = nf90_put_att(ncid,varid,'time_rep','averaged')
            endif
          endif
        enddo  ! num_avail_hist_fields_3Dc

        dimidz(1) = imtid
        dimidz(2) = jmtid
        dimidz(3) = kmtidi
        dimidz(4) = timid

        do n = n3Dccum + 1, n3Dzcum
          if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            status  = nf90_def_var(ncid, avail_hist_fields(n)%vname, &
                         nf90_float, dimidz, varid)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining variable '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'units', &
                        avail_hist_fields(n)%vunit)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining units for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid, 'long_name', &
                        avail_hist_fields(n)%vdesc)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining long_name for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'coordinates', &
                        avail_hist_fields(n)%vcoord)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining coordinates for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'cell_measures', &
                        avail_hist_fields(n)%vcellmeas)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining cell measures for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'missing_value',spval)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining missing_value for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'_FillValue',spval)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining _FillValue for '//avail_hist_fields(n)%vname)

          endif
        enddo  ! num_avail_hist_fields_3Dz
        
        dimidz(1) = imtid
        dimidz(2) = jmtid
        dimidz(3) = kmtidb
        dimidz(4) = timid

        do n = n3Dzcum + 1, n3Dbcum
          if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            status  = nf90_def_var(ncid, avail_hist_fields(n)%vname, &
                         nf90_float, dimidz, varid)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining variable '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'units', &
                        avail_hist_fields(n)%vunit)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining units for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid, 'long_name', &
                        avail_hist_fields(n)%vdesc)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining long_name for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'coordinates', &
                        avail_hist_fields(n)%vcoord)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining coordinates for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'cell_measures', &
                        avail_hist_fields(n)%vcellmeas)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining cell measures for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'missing_value',spval)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining missing_value for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'_FillValue',spval)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining _FillValue for '//avail_hist_fields(n)%vname)

          endif
        enddo  ! num_avail_hist_fields_3Db

        dimidz(1) = imtid
        dimidz(2) = jmtid
        dimidz(3) = kmtida
        dimidz(4) = timid

        do n = n3Dbcum + 1, n3Dacum
          if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            status  = nf90_def_var(ncid, avail_hist_fields(n)%vname, &
                         nf90_float, dimidz, varid)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining variable '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'units', &
                        avail_hist_fields(n)%vunit)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining units for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid, 'long_name', &
                        avail_hist_fields(n)%vdesc)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining long_name for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'coordinates', &
                        avail_hist_fields(n)%vcoord)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining coordinates for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'cell_measures', &
                        avail_hist_fields(n)%vcellmeas)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining cell measures for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'missing_value',spval)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining missing_value for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'_FillValue',spval)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining _FillValue for '//avail_hist_fields(n)%vname)

          endif
        enddo  ! num_avail_hist_fields_3Da
      
        dimidcz(1) = imtid
        dimidcz(2) = jmtid
        dimidcz(3) = kmtidi
        dimidcz(4) = cmtid
        dimidcz(5) = timid

        do n = n3Dacum + 1, n4Dicum
          if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            status  = nf90_def_var(ncid, avail_hist_fields(n)%vname, &
!                             nf90_float, dimidcz, varid)
                             nf90_float, dimidcz(1:4), varid) ! ferret    
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining variable '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'units', &
                        avail_hist_fields(n)%vunit)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining units for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid, 'long_name', &
                        avail_hist_fields(n)%vdesc)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining long_name for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'coordinates', &
                        avail_hist_fields(n)%vcoord)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining coordinates for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'cell_measures', &
                        avail_hist_fields(n)%vcellmeas)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining cell measures for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'missing_value',spval)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining missing_value for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'_FillValue',spval)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining _FillValue for '//avail_hist_fields(n)%vname)

      !-----------------------------------------------------------------
      ! Add cell_methods attribute to variables if averaged
      !-----------------------------------------------------------------
            if (hist_avg) then
                status = nf90_put_att(ncid,varid,'cell_methods','time: mean')
                if (status /= nf90_noerr) call abort_ice(subname// &
                 'ERROR: defining cell methods for '//avail_hist_fields(n)%vname)
            endif

            if (histfreq(ns) == '1' .or. .not. hist_avg) then
               status = nf90_put_att(ncid,varid,'time_rep','instantaneous')
            else
               status = nf90_put_att(ncid,varid,'time_rep','averaged')
            endif
          endif
        enddo  ! num_avail_hist_fields_4Di

        dimidcz(1) = imtid
        dimidcz(2) = jmtid
        dimidcz(3) = kmtids
        dimidcz(4) = cmtid
        dimidcz(5) = timid

        do n = n4Dicum + 1, n4Dscum
          if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            status  = nf90_def_var(ncid, avail_hist_fields(n)%vname, &
!                             nf90_float, dimidcz, varid)
                             nf90_float, dimidcz(1:4), varid) ! ferret    
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining variable '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'units', &
                        avail_hist_fields(n)%vunit)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining units for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid, 'long_name', &
                        avail_hist_fields(n)%vdesc)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining long_name for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'coordinates', &
                        avail_hist_fields(n)%vcoord)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining coordinates for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'cell_measures', &
                        avail_hist_fields(n)%vcellmeas)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining cell measures for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'missing_value',spval)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining missing_value for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'_FillValue',spval)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: defining _FillValue for '//avail_hist_fields(n)%vname)

      !-----------------------------------------------------------------
      ! Add cell_methods attribute to variables if averaged
      !-----------------------------------------------------------------
            if (hist_avg) then
                status = nf90_put_att(ncid,varid,'cell_methods','time: mean')
                if (status /= nf90_noerr) call abort_ice(subname// &
                 'ERROR: defining cell methods for '//avail_hist_fields(n)%vname)
            endif

            if (histfreq(ns) == '1' .or. .not. hist_avg) then
               status = nf90_put_att(ncid,varid,'time_rep','instantaneous')
            else
               status = nf90_put_att(ncid,varid,'time_rep','averaged')
            endif
          endif
        enddo  ! num_avail_hist_fields_4Ds

      !-----------------------------------------------------------------
      ! global attributes
      !-----------------------------------------------------------------
      ! ... the user should change these to something useful ...
      !-----------------------------------------------------------------
#ifdef CESMCOUPLED
        status = nf90_put_att(ncid,nf90_global,'title',runid)
        if (status /= nf90_noerr) call abort_ice(subname// &
                      'ERROR: in global attribute title')
#else
        title  = 'sea ice model output for CICE'
        status = nf90_put_att(ncid,nf90_global,'title',title)
        if (status /= nf90_noerr) call abort_ice(subname// &
                      'ERROR: in global attribute title')
#endif
        title = 'Diagnostic and Prognostic Variables'
        status = nf90_put_att(ncid,nf90_global,'contents',title)
        if (status /= nf90_noerr) call abort_ice(subname// &
                      'ERROR: global attribute contents')

        write(title,'(2a)') 'Los Alamos Sea Ice Model, ', trim(version_name)
        status = nf90_put_att(ncid,nf90_global,'source',title)
        if (status /= nf90_noerr) call abort_ice(subname// &
                      'ERROR: global attribute source')

        if (use_leap_years) then
          write(title,'(a,i3,a)') 'This year has ',int(dayyr),' days'
        else
          write(title,'(a,i3,a)') 'All years have exactly ',int(dayyr),' days'
        endif
        status = nf90_put_att(ncid,nf90_global,'comment',title)
        if (status /= nf90_noerr) call abort_ice(subname// &
                      'ERROR: global attribute comment')

        write(title,'(a,i8.8)') 'File written on model date ',idate
        status = nf90_put_att(ncid,nf90_global,'comment2',title)
        if (status /= nf90_noerr) call abort_ice(subname// &
                      'ERROR: global attribute date1')

        write(title,'(a,i6)') 'seconds elapsed into model date: ',sec
        status = nf90_put_att(ncid,nf90_global,'comment3',title)
        if (status /= nf90_noerr) call abort_ice(subname// &
                      'ERROR: global attribute date2')

        title = 'CF-1.0'
        status =  &
             nf90_put_att(ncid,nf90_global,'conventions',title)
        if (status /= nf90_noerr) call abort_ice(subname// &
             'ERROR: in global attribute conventions')

        call date_and_time(date=current_date, time=current_time)
        write(start_time,1000) current_date(1:4), current_date(5:6), &
                               current_date(7:8), current_time(1:2), &
                               current_time(3:4), current_time(5:8)
1000    format('This dataset was created on ', &
                a,'-',a,'-',a,' at ',a,':',a,':',a)

        status = nf90_put_att(ncid,nf90_global,'history',start_time)
        if (status /= nf90_noerr) call abort_ice(subname// &
                      'ERROR: global attribute history')

        status = nf90_put_att(ncid,nf90_global,'io_flavor','io_netcdf')
        if (status /= nf90_noerr) call abort_ice(subname// &
                      'ERROR: global attribute io_flavor')

      !-----------------------------------------------------------------
      ! end define mode
      !-----------------------------------------------------------------

        status = nf90_enddef(ncid)
        if (status /= nf90_noerr) call abort_ice(subname//'ERROR in nf90_enddef')

      !-----------------------------------------------------------------
      ! write time variable
      !-----------------------------------------------------------------

        status = nf90_inq_varid(ncid,'time',varid)
        if (status /= nf90_noerr) call abort_ice(subname// &
                      'ERROR: getting time varid')
!sgl        status = nf90_put_var(ncid,varid,ltime)
        status = nf90_put_var(ncid,varid,ltime2)
        if (status /= nf90_noerr) call abort_ice(subname// &
                      'ERROR: writing time variable')

      !-----------------------------------------------------------------
      ! write time_bounds info
      !-----------------------------------------------------------------

        if (hist_avg) then
          status = nf90_inq_varid(ncid,'time_bounds',varid)
          if (status /= nf90_noerr) call abort_ice(subname// &
                        'ERROR: getting time_bounds id')
          status = nf90_put_var(ncid,varid,time_beg(ns),start=(/1/))
          if (status /= nf90_noerr) call abort_ice(subname// &
                        'ERROR: writing time_beg')
          status = nf90_put_var(ncid,varid,time_end(ns),start=(/2/))
          if (status /= nf90_noerr) call abort_ice(subname// &
                        'ERROR: writing time_end')
        endif

      endif                     ! master_task

      if (my_task==master_task) then
         allocate(work_g1(nx_global,ny_global))
         allocate(work_gr(nx_global,ny_global))
      else
         allocate(work_gr(1,1))   ! to save memory
         allocate(work_g1(1,1))
      endif

      work_g1(:,:) = c0

      !-----------------------------------------------------------------
      ! write coordinate variables
      !-----------------------------------------------------------------

        do i = 1,ncoord
          call broadcast_scalar(coord_var(i)%short_name,master_task)
          SELECT CASE (coord_var(i)%short_name)
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
          END SELECT
          
          if (my_task == master_task) then
             work_gr = work_g1
             status = nf90_inq_varid(ncid, coord_var(i)%short_name, varid)
             if (status /= nf90_noerr) call abort_ice(subname// &
                  'ERROR: getting varid for '//coord_var(i)%short_name)
             status = nf90_put_var(ncid,varid,work_gr)
             if (status /= nf90_noerr) call abort_ice(subname// &
                  'ERROR: writing'//coord_var(i)%short_name)
          endif
        enddo

        ! Extra dimensions (NCAT, VGRD*)

        do i = 1, nvarz
          if (igrdz(i)) then
          call broadcast_scalar(var_nz(i)%short_name,master_task)
          if (my_task == master_task) then
             status = nf90_inq_varid(ncid, var_nz(i)%short_name, varid)
             if (status /= nf90_noerr) call abort_ice(subname// &
                  'ERROR: getting varid for '//var_nz(i)%short_name)
             SELECT CASE (var_nz(i)%short_name)
               CASE ('NCAT') 
                 status = nf90_put_var(ncid,varid,hin_max(1:ncat_hist))
               CASE ('VGRDi') ! index - needed for Met Office analysis code
                 status = nf90_put_var(ncid,varid,(/(k, k=1,nzilyr)/))
               CASE ('VGRDs') ! index - needed for Met Office analysis code
                 status = nf90_put_var(ncid,varid,(/(k, k=1,nzslyr)/))
               CASE ('VGRDb')
                 status = nf90_put_var(ncid,varid,(/(k, k=1,nzblyr)/))
               CASE ('VGRDa') 
                 status = nf90_put_var(ncid,varid,(/(k, k=1,nzalyr)/))
             END SELECT
             if (status /= nf90_noerr) call abort_ice(subname// &
                           'ERROR: writing'//var_nz(i)%short_name)
          endif
          endif
        enddo

      !-----------------------------------------------------------------
      ! write grid masks, area and rotation angle
      !-----------------------------------------------------------------

      if (igrd(n_tmask)) then
      call gather_global(work_g1, hm, master_task, distrb_info)
      if (my_task == master_task) then
        work_gr=work_g1
        status = nf90_inq_varid(ncid, 'tmask', varid)
        if (status /= nf90_noerr) call abort_ice(subname// &
                      'ERROR: getting varid for tmask')
        status = nf90_put_var(ncid,varid,work_gr)
        if (status /= nf90_noerr) call abort_ice(subname// &
                      'ERROR: writing variable tmask')
      endif
      endif

      if (igrd(n_blkmask)) then
      call gather_global(work_g1, bm, master_task, distrb_info)
      if (my_task == master_task) then
        work_gr=work_g1
        status = nf90_inq_varid(ncid, 'blkmask', varid)
        if (status /= nf90_noerr) call abort_ice(subname// &
                      'ERROR: getting varid for blkmask')
        status = nf90_put_var(ncid,varid,work_gr)
        if (status /= nf90_noerr) call abort_ice(subname// &
                      'ERROR: writing variable blkmask')
      endif
      endif

      do i = 3, nvar      ! note n_tmask=1, n_blkmask=2
        if (igrd(i)) then
        call broadcast_scalar(var(i)%req%short_name,master_task)
        SELECT CASE (var(i)%req%short_name)
          CASE ('tarea')
            call gather_global(work_g1, tarea, master_task, distrb_info)
          CASE ('uarea')
            call gather_global(work_g1, uarea, master_task, distrb_info)
          CASE ('dxu')
            call gather_global(work_g1,   dxu, master_task, distrb_info)
          CASE ('dyu')
            call gather_global(work_g1,   dyu, master_task, distrb_info)
          CASE ('dxt')
            call gather_global(work_g1,   dxt, master_task, distrb_info)
          CASE ('dyt')
            call gather_global(work_g1,   dyt, master_task, distrb_info)
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
          work_gr=work_g1
          status = nf90_inq_varid(ncid, var(i)%req%short_name, varid)
          if (status /= nf90_noerr) call abort_ice(subname// &
                        'ERROR: getting varid for '//var(i)%req%short_name)
          status = nf90_put_var(ncid,varid,work_gr)
          if (status /= nf90_noerr) call abort_ice(subname// &
                        'ERROR: writing variable '//var(i)%req%short_name)
        endif
        endif
      enddo

      deallocate(work_gr)

      !----------------------------------------------------------------
      ! Write coordinates of grid box vertices
      !----------------------------------------------------------------

      if (f_bounds) then
      if (my_task==master_task) then
         allocate(work_gr3(nverts,nx_global,ny_global))
      else
         allocate(work_gr3(1,1,1))   ! to save memory
      endif

      work_gr3(:,:,:) = c0
      work1   (:,:,:) = c0

      do i = 1, nvar_verts
        call broadcast_scalar(var_nverts(i)%short_name,master_task)
        SELECT CASE (var_nverts(i)%short_name)
        CASE ('lont_bounds')
        do ivertex = 1, nverts 
           work1(:,:,:) = lont_bounds(ivertex,:,:,:)
           call gather_global(work_g1, work1, master_task, distrb_info)
           if (my_task == master_task) work_gr3(ivertex,:,:) = work_g1(:,:)
        enddo
        CASE ('latt_bounds')
        do ivertex = 1, nverts
           work1(:,:,:) = latt_bounds(ivertex,:,:,:)
           call gather_global(work_g1, work1, master_task, distrb_info)
           if (my_task == master_task) work_gr3(ivertex,:,:) = work_g1(:,:)
        enddo
        CASE ('lonu_bounds')
        do ivertex = 1, nverts
           work1(:,:,:) = lonu_bounds(ivertex,:,:,:)
           call gather_global(work_g1, work1, master_task, distrb_info)
           if (my_task == master_task) work_gr3(ivertex,:,:) = work_g1(:,:)
        enddo
        CASE ('latu_bounds')
        do ivertex = 1, nverts
           work1(:,:,:) = latu_bounds(ivertex,:,:,:)
           call gather_global(work_g1, work1, master_task, distrb_info)
           if (my_task == master_task) work_gr3(ivertex,:,:) = work_g1(:,:)
        enddo
        END SELECT

        if (my_task == master_task) then
          status = nf90_inq_varid(ncid, var_nverts(i)%short_name, varid)
          if (status /= nf90_noerr) call abort_ice(subname// &
             'ERROR: getting varid for '//var_nverts(i)%short_name)
          status = nf90_put_var(ncid,varid,work_gr3)
          if (status /= nf90_noerr) call abort_ice(subname// &
             'ERROR: writing variable '//var_nverts(i)%short_name)
        endif
      enddo
      deallocate(work_gr3)
      endif

      !-----------------------------------------------------------------
      ! write variable data
      !-----------------------------------------------------------------

      if (my_task==master_task) then
         allocate(work_gr(nx_global,ny_global))
      else
         allocate(work_gr(1,1))     ! to save memory
      endif
      work_gr(:,:) = c0
      work_g1(:,:) = c0

      do n=1,num_avail_hist_fields_2D
        if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
          call gather_global(work_g1, a2D(:,:,n,:), &
                             master_task, distrb_info)
          if (my_task == master_task) then
            work_gr(:,:) = work_g1(:,:)
            status  = nf90_inq_varid(ncid,avail_hist_fields(n)%vname,varid)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: getting varid for '//avail_hist_fields(n)%vname)
            status  = nf90_put_var(ncid,varid,work_gr(:,:), &
                                   count=(/nx_global,ny_global/))
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: writing variable '//avail_hist_fields(n)%vname)
          endif
        endif
      enddo ! num_avail_hist_fields_2D

      work_gr(:,:) = c0
      work_g1(:,:) = c0

      do n = n2D + 1, n3Dccum
        nn = n - n2D
        if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
          if (my_task == master_task) then
            status  = nf90_inq_varid(ncid,avail_hist_fields(n)%vname,varid)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: getting varid for '//avail_hist_fields(n)%vname)
          endif
          do k = 1, ncat_hist
             call gather_global(work_g1, a3Dc(:,:,k,nn,:), &
                                master_task, distrb_info)
             work_gr(:,:) = work_g1(:,:)

             if (my_task == master_task) then
             status  = nf90_inq_varid(ncid,avail_hist_fields(n)%vname,varid)
             if (status /= nf90_noerr) call abort_ice(subname// &
                'ERROR: getting varid for '//avail_hist_fields(n)%vname)
             status  = nf90_put_var(ncid,varid,work_gr(:,:), &
                                    start=(/        1,        1,k/), &
                                    count=(/nx_global,ny_global,1/))
             if (status /= nf90_noerr) call abort_ice(subname// &
                'ERROR: writing variable '//avail_hist_fields(n)%vname)
             endif
          enddo ! k
        endif
      enddo ! num_avail_hist_fields_3Dc

      work_gr(:,:) = c0
      work_g1(:,:) = c0

      do n = n3Dccum+1, n3Dzcum
        nn = n - n3Dccum
        if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
          if (my_task == master_task) then
            status  = nf90_inq_varid(ncid,avail_hist_fields(n)%vname,varid)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: getting varid for '//avail_hist_fields(n)%vname)
          endif
          do k = 1, nzilyr
             call gather_global(work_g1, a3Dz(:,:,k,nn,:), &
                                master_task, distrb_info)
             work_gr(:,:) = work_g1(:,:)

             if (my_task == master_task) then
             status  = nf90_put_var(ncid,varid,work_gr(:,:), &
                                    start=(/        1,        1,k/), &
                                    count=(/nx_global,ny_global,1/))
             if (status /= nf90_noerr) call abort_ice(subname// &
                'ERROR: writing variable '//avail_hist_fields(n)%vname)
           endif
           enddo ! k
        endif
      enddo ! num_avail_hist_fields_3Dz

      work_gr(:,:) = c0
      work_g1(:,:) = c0

     do n = n3Dzcum+1, n3Dbcum
        nn = n - n3Dzcum
        if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
          if (my_task == master_task) then
            status  = nf90_inq_varid(ncid,avail_hist_fields(n)%vname,varid)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: getting varid for '//avail_hist_fields(n)%vname)
          endif
          do k = 1, nzblyr
             call gather_global(work_g1, a3Db(:,:,k,nn,:), &
                                master_task, distrb_info)
             work_gr(:,:) = work_g1(:,:)

             if (my_task == master_task) then
             status  = nf90_put_var(ncid,varid,work_gr(:,:), &
                                    start=(/        1,        1,k/), &
                                    count=(/nx_global,ny_global,1/))
             if (status /= nf90_noerr) call abort_ice(subname// &
                'ERROR: writing variable '//avail_hist_fields(n)%vname)
           endif
           enddo ! k
        endif
      enddo ! num_avail_hist_fields_3Db

      work_gr(:,:) = c0
      work_g1(:,:) = c0

      do n = n3Dbcum+1, n3Dacum
        nn = n - n3Dbcum
        if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
          if (my_task == master_task) then
            status  = nf90_inq_varid(ncid,avail_hist_fields(n)%vname,varid)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: getting varid for '//avail_hist_fields(n)%vname)
          endif
          do k = 1, nzalyr
             call gather_global(work_g1, a3Da(:,:,k,nn,:), &
                                master_task, distrb_info)
             work_gr(:,:) = work_g1(:,:)

             if (my_task == master_task) then
             status  = nf90_put_var(ncid,varid,work_gr(:,:), &
                                    start=(/        1,        1,k/), &
                                    count=(/nx_global,ny_global,1/))
             if (status /= nf90_noerr) call abort_ice(subname// &
                'ERROR: writing variable '//avail_hist_fields(n)%vname)
           endif
           enddo ! k
        endif
      enddo ! num_avail_hist_fields_3Da

      work_gr(:,:) = c0
      work_g1(:,:) = c0

      do n = n3Dacum+1, n4Dicum
        nn = n - n3Dacum
        if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
          if (my_task == master_task) then
            status  = nf90_inq_varid(ncid,avail_hist_fields(n)%vname,varid)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: getting varid for '//avail_hist_fields(n)%vname)
          endif
          do ic = 1, ncat_hist
             do k = 1, nzilyr
                call gather_global(work_g1, a4Di(:,:,k,ic,nn,:), &
                                master_task, distrb_info)
                work_gr(:,:) = work_g1(:,:)
                if (my_task == master_task) then
                  status  = nf90_put_var(ncid,varid,work_gr(:,:), &
                                         start=(/        1,        1,k,ic/), &
                                         count=(/nx_global,ny_global,1, 1/))
                  if (status /= nf90_noerr) call abort_ice(subname// &
                     'ERROR: writing variable '//avail_hist_fields(n)%vname)
                endif
             enddo ! k
          enddo ! ic
        endif
      enddo ! num_avail_hist_fields_4Di

      work_gr(:,:) = c0
      work_g1(:,:) = c0

      do n = n4Dicum+1, n4Dscum
        nn = n - n4Dicum
        if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
          if (my_task == master_task) then
            status  = nf90_inq_varid(ncid,avail_hist_fields(n)%vname,varid)
            if (status /= nf90_noerr) call abort_ice(subname// &
               'ERROR: getting varid for '//avail_hist_fields(n)%vname)
          endif
          do ic = 1, ncat_hist
             do k = 1, nzslyr
                call gather_global(work_g1, a4Ds(:,:,k,ic,nn,:), &
                                master_task, distrb_info)
                work_gr(:,:) = work_g1(:,:)
                if (my_task == master_task) then
                  status  = nf90_put_var(ncid,varid,work_gr(:,:), &
                                         start=(/        1,        1,k,ic/), &
                                         count=(/nx_global,ny_global,1, 1/))
                  if (status /= nf90_noerr) call abort_ice(subname// &
                     'ERROR: writing variable '//avail_hist_fields(n)%vname)
                endif
             enddo ! k
          enddo ! ic
        endif
      enddo ! num_avail_hist_fields_4Ds

      deallocate(work_gr)
      deallocate(work_g1)

      !-----------------------------------------------------------------
      ! close output dataset
      !-----------------------------------------------------------------

      if (my_task == master_task) then
         status = nf90_close(ncid)
         if (status /= nf90_noerr) call abort_ice(subname// &
                       'ERROR: closing netCDF history file')
         write(nu_diag,*) ' '
         write(nu_diag,*) 'Finished writing ',trim(ncfile(ns))
      endif
#endif

      end subroutine ice_write_hist

!=======================================================================

      end module ice_history_write

!=======================================================================
