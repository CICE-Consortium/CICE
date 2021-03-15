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
!
      module ice_history_write

      use ice_kinds_mod
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

      use ice_blocks, only: nx_block, ny_block
      use ice_broadcast, only: broadcast_scalar
      use ice_calendar, only: msec, timesecs, idate, idate0, write_ic, &
          histfreq, days_per_year, use_leap_years, dayyr
      use ice_communicate, only: my_task, master_task
      use ice_constants, only: c0, c360, spval, spval_dbl
      use ice_domain, only: distrb_info, nblocks
      use ice_domain_size, only: nx_global, ny_global, max_blocks, max_nstrm
      use ice_gather_scatter, only: gather_global
      use ice_grid, only: TLON, TLAT, ULON, ULAT, hm, bm, tarea, uarea, &
          dxu, dxt, dyu, dyt, HTN, HTE, ANGLE, ANGLET, tmask, &
          lont_bounds, latt_bounds, lonu_bounds, latu_bounds
      use ice_history_shared
      use ice_arrays_column, only: hin_max, floe_rad_c
      use ice_restart_shared, only: runid, lcdf64
      use ice_pio
      use pio

      integer (kind=int_kind), intent(in) :: ns

      ! local variables

      integer (kind=int_kind) :: i,j,k,ic,n,nn, &
         ncid,status,imtid,jmtid,kmtidi,kmtids,kmtidb, cmtid,timid, &
         length,nvertexid,ivertex,kmtida,fmtid
      integer (kind=int_kind), dimension(2) :: dimid2
      integer (kind=int_kind), dimension(3) :: dimid3
      integer (kind=int_kind), dimension(4) :: dimidz
      integer (kind=int_kind), dimension(5) :: dimidcz
      integer (kind=int_kind), dimension(3) :: dimid_nverts
      integer (kind=int_kind), dimension(6) :: dimidex
      real (kind=real_kind) :: ltime
      real (kind= dbl_kind) :: ltime2
      character (char_len) :: title
      character (char_len_long) :: ncfile(max_nstrm)
      integer (kind=int_kind) :: iotype

      integer (kind=int_kind) :: icategory,ind,i_aice,boundid

      character (char_len) :: start_time,current_date,current_time
      character (len=16) :: c_aice
      character (len=8) :: cdate

      type(file_desc_t)     :: File
      type(io_desc_t)       :: iodesc2d, &
                               iodesc3dc, iodesc3dv, iodesc3di, iodesc3db, iodesc3da, &
                               iodesc3df, &
                               iodesc4di, iodesc4ds, iodesc4df
      type(var_desc_t)      :: varid

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

      real (kind=dbl_kind), allocatable :: workr2(:,:,:)
      real (kind=dbl_kind), allocatable :: workr3(:,:,:,:)
      real (kind=dbl_kind), allocatable :: workr4(:,:,:,:,:)
      real (kind=dbl_kind), allocatable :: workr3v(:,:,:,:)

      character(len=char_len_long) :: &
           filename

      integer (kind=int_kind), dimension(1) ::  &
         tim_start,tim_length          ! dimension quantities for netCDF

      integer (kind=int_kind), dimension(2) ::  &
         bnd_start,bnd_length          ! dimension quantities for netCDF

      real (kind=dbl_kind) :: secday
      real (kind=dbl_kind) :: rad_to_deg

      integer (kind=int_kind) :: lprecision

      character(len=*), parameter :: subname = '(ice_write_hist)'

      call icepack_query_parameters(secday_out=secday)
      call icepack_query_parameters(rad_to_deg_out=rad_to_deg)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__, line=__LINE__)

      if (my_task == master_task) then
        call construct_filename(ncfile(ns),'nc',ns)

        ! add local directory path name to ncfile
        if (write_ic) then
          ncfile(ns) = trim(incond_dir)//ncfile(ns)
        else
          ncfile(ns) = trim(history_dir)//ncfile(ns)
        endif
        filename = ncfile(ns)
      end if
      call broadcast_scalar(filename, master_task)

      ! create file

      iotype = PIO_IOTYPE_NETCDF
      if (history_format == 'pio_pnetcdf') iotype = PIO_IOTYPE_PNETCDF
      File%fh=-1
      call ice_pio_init(mode='write', filename=trim(filename), File=File, &
        clobber=.true., cdf64=lcdf64, iotype=iotype)

      call ice_pio_initdecomp(iodesc=iodesc2d)
      call ice_pio_initdecomp(ndim3=ncat_hist, iodesc=iodesc3dc)
      call ice_pio_initdecomp(ndim3=nzilyr,    iodesc=iodesc3di)
      call ice_pio_initdecomp(ndim3=nzblyr,    iodesc=iodesc3db)
      call ice_pio_initdecomp(ndim3=nzalyr,    iodesc=iodesc3da)
      call ice_pio_initdecomp(ndim3=nfsd_hist, iodesc=iodesc3df)
      call ice_pio_initdecomp(ndim3=nverts,    iodesc=iodesc3dv, inner_dim=.true.)
      call ice_pio_initdecomp(ndim3=nzilyr,    ndim4=ncat_hist, iodesc=iodesc4di)
      call ice_pio_initdecomp(ndim3=nzslyr,    ndim4=ncat_hist, iodesc=iodesc4ds)
      call ice_pio_initdecomp(ndim3=nfsd_hist, ndim4=ncat_hist, iodesc=iodesc4df)

      ltime2 = timesecs/secday
      ltime  = real(timesecs/secday,kind=real_kind)

      ! option of turning on double precision history files
      lprecision = pio_real
      if (history_precision == 8) lprecision = pio_double

      !-----------------------------------------------------------------
      ! define dimensions
      !-----------------------------------------------------------------

       if (hist_avg .and. histfreq(ns) /= '1') then
          status = pio_def_dim(File,'d2',2,boundid)
        endif

        status = pio_def_dim(File,'ni',nx_global,imtid)
        status = pio_def_dim(File,'nj',ny_global,jmtid)
        status = pio_def_dim(File,'nc',ncat_hist,cmtid)
        status = pio_def_dim(File,'nkice',nzilyr,kmtidi)
        status = pio_def_dim(File,'nksnow',nzslyr,kmtids)
        status = pio_def_dim(File,'nkbio',nzblyr,kmtidb)
        status = pio_def_dim(File,'nkaer',nzalyr,kmtida)
        status = pio_def_dim(File,'time',PIO_UNLIMITED,timid)
        status = pio_def_dim(File,'nvertices',nverts,nvertexid)
        status = pio_def_dim(File,'nf',nfsd_hist,fmtid)

      !-----------------------------------------------------------------
      ! define coordinate variables:  time, time_bounds
      !-----------------------------------------------------------------

!sgl        status = pio_def_var(File,'time',pio_real,(/timid/),varid)
        status = pio_def_var(File,'time',pio_double,(/timid/),varid)
        status = pio_put_att(File,varid,'long_name','model time')

        write(cdate,'(i8.8)') idate0
        write(title,'(a,a,a,a,a,a,a)') 'days since ', &
              cdate(1:4),'-',cdate(5:6),'-',cdate(7:8),' 00:00:00'
        status = pio_put_att(File,varid,'units',trim(title))

        if (days_per_year == 360) then
           status = pio_put_att(File,varid,'calendar','360_day')
        elseif (days_per_year == 365 .and. .not.use_leap_years ) then
           status = pio_put_att(File,varid,'calendar','NoLeap')
        elseif (use_leap_years) then
           status = pio_put_att(File,varid,'calendar','Gregorian')
        else
           call abort_ice(subname//'ERROR: invalid calendar settings')
        endif

        if (hist_avg .and. histfreq(ns) /= '1') then
          status = pio_put_att(File,varid,'bounds','time_bounds')
        endif

        ! Define attributes for time_bounds if hist_avg is true
        if (hist_avg .and. histfreq(ns) /= '1') then
          dimid2(1) = boundid
          dimid2(2) = timid
!sgl          status = pio_def_var(File,'time_bounds',pio_real,dimid2,varid)
          status = pio_def_var(File,'time_bounds',pio_double,dimid2,varid)
          status = pio_put_att(File,varid,'long_name', &
                                'boundaries for time-averaging interval')
          write(cdate,'(i8.8)') idate0
          write(title,'(a,a,a,a,a,a,a,a)') 'days since ', &
                cdate(1:4),'-',cdate(5:6),'-',cdate(7:8),' 00:00:00'
          status = pio_put_att(File,varid,'units',trim(title))
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
      var_nz(6) = coord_attributes('NFSD', 'category floe size (center)', 'm')

      !-----------------------------------------------------------------
      ! define information for optional time-invariant variables
      !-----------------------------------------------------------------

      var(n_tmask)%req = coord_attributes('tmask', &
                  'ocean grid mask', ' ')
      var(n_tmask)%coordinates = 'TLON TLAT'

      var(n_blkmask)%req = coord_attributes('blkmask', &
                  'ice grid block mask', ' ')
      var(n_blkmask)%coordinates = 'TLON TLAT'

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

        dimid2(1) = imtid
        dimid2(2) = jmtid

        do i = 1, ncoord
          status = pio_def_var(File, trim(coord_var(i)%short_name), lprecision, &
                                dimid2, varid)
          status = pio_put_att(File,varid,'long_name',trim(coord_var(i)%long_name))
          status = pio_put_att(File, varid, 'units', trim(coord_var(i)%units))
          if (lprecision == pio_real) then
             status = pio_put_att(File, varid, 'missing_value', spval)
             status = pio_put_att(File, varid,'_FillValue',spval)
          else
             status = pio_put_att(File, varid, 'missing_value', spval_dbl)
             status = pio_put_att(File, varid,'_FillValue',spval_dbl)
          endif
          if (coord_var(i)%short_name == 'ULAT') then
             status = pio_put_att(File,varid,'comment', &
                  trim('Latitude of NE corner of T grid cell'))
          endif
          if (f_bounds) then
              status = pio_put_att(File, varid, 'bounds', trim(coord_bounds(i)))
          endif          
        enddo

        ! Extra dimensions (NCAT, NZILYR, NZSLYR, NZBLYR, NZALYR, NFSD)
          dimidex(1)=cmtid
          dimidex(2)=kmtidi
          dimidex(3)=kmtids
          dimidex(4)=kmtidb
          dimidex(5)=kmtida
          dimidex(6)=fmtid

        do i = 1, nvarz
           if (igrdz(i)) then
              status = pio_def_var(File, trim(var_nz(i)%short_name), lprecision, &
                                   (/dimidex(i)/), varid)
              status = pio_put_att(File, varid, 'long_name', var_nz(i)%long_name)
              status = pio_put_att(File, varid, 'units'    , var_nz(i)%units)
           endif
        enddo

        ! Attributes for tmask defined separately, since it has no units
        if (igrd(n_tmask)) then
           status = pio_def_var(File, 'tmask', lprecision, dimid2, varid)
           status = pio_put_att(File,varid, 'long_name', 'ocean grid mask') 
           status = pio_put_att(File, varid, 'coordinates', 'TLON TLAT')
           if (lprecision == pio_real) then
              status = pio_put_att(File, varid, 'missing_value', spval)
              status = pio_put_att(File, varid,'_FillValue',spval)
           else
              status = pio_put_att(File, varid, 'missing_value', spval_dbl)
              status = pio_put_att(File, varid,'_FillValue',spval_dbl)
           endif
           status = pio_put_att(File,varid,'comment', '0 = land, 1 = ocean')
        endif
        if (igrd(n_blkmask)) then
           status = pio_def_var(File, 'blkmask', lprecision, dimid2, varid)
           status = pio_put_att(File,varid, 'long_name', 'ice grid block mask') 
           status = pio_put_att(File, varid, 'coordinates', 'TLON TLAT')
           status = pio_put_att(File,varid,'comment', 'mytask + iblk/100')
           if (lprecision == pio_real) then
              status = pio_put_att(File, varid, 'missing_value', spval)
              status = pio_put_att(File, varid,'_FillValue',spval)
           else
              status = pio_put_att(File, varid, 'missing_value', spval_dbl)
              status = pio_put_att(File, varid,'_FillValue',spval_dbl)
           endif
        endif

        do i = 3, nvar       ! note: n_tmask=1, n_blkmask=2
          if (igrd(i)) then
             status = pio_def_var(File, trim(var(i)%req%short_name), &
                                   lprecision, dimid2, varid)
             status = pio_put_att(File,varid, 'long_name', trim(var(i)%req%long_name))
             status = pio_put_att(File, varid, 'units', trim(var(i)%req%units))
             status = pio_put_att(File, varid, 'coordinates', trim(var(i)%coordinates))
             if (lprecision == pio_real) then
                status = pio_put_att(File, varid, 'missing_value', spval)
                status = pio_put_att(File, varid,'_FillValue',spval)
             else
                status = pio_put_att(File, varid, 'missing_value', spval_dbl)
                status = pio_put_att(File, varid,'_FillValue',spval_dbl)
             endif
          endif
        enddo

        ! Fields with dimensions (nverts,nx,ny)
        dimid_nverts(1) = nvertexid
        dimid_nverts(2) = imtid
        dimid_nverts(3) = jmtid
        do i = 1, nvar_verts
          if (f_bounds) then
             status = pio_def_var(File, trim(var_nverts(i)%short_name), &
                                   lprecision,dimid_nverts, varid)
             status = & 
             pio_put_att(File,varid, 'long_name', trim(var_nverts(i)%long_name))
             status = &
             pio_put_att(File, varid, 'units', trim(var_nverts(i)%units))
             if (lprecision == pio_real) then
                status = pio_put_att(File, varid, 'missing_value', spval)
                status = pio_put_att(File, varid,'_FillValue',spval)
             else
                status = pio_put_att(File, varid, 'missing_value', spval_dbl)
                status = pio_put_att(File, varid,'_FillValue',spval_dbl)
             endif
          endif
        enddo
 
      !-----------------------------------------------------------------
      ! define attributes for time-variant variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! 2D
      !-----------------------------------------------------------------

        dimid3(1) = imtid
        dimid3(2) = jmtid
        dimid3(3) = timid

        do n=1,num_avail_hist_fields_2D
          if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            status  = pio_def_var(File, trim(avail_hist_fields(n)%vname), &
                         lprecision, dimid3, varid)
            status = pio_put_att(File,varid,'units', &
                        trim(avail_hist_fields(n)%vunit))
            status = pio_put_att(File,varid, 'long_name', &
                        trim(avail_hist_fields(n)%vdesc))
            status = pio_put_att(File,varid,'coordinates', &
                        trim(avail_hist_fields(n)%vcoord))
            status = pio_put_att(File,varid,'cell_measures', &
                        trim(avail_hist_fields(n)%vcellmeas))
            if (lprecision == pio_real) then
               status = pio_put_att(File, varid, 'missing_value', spval)
               status = pio_put_att(File, varid,'_FillValue',spval)
            else
               status = pio_put_att(File, varid, 'missing_value', spval_dbl)
               status = pio_put_att(File, varid,'_FillValue',spval_dbl)
            endif

            ! Add cell_methods attribute to variables if averaged
            if (hist_avg .and. histfreq(ns) /= '1') then
              if (TRIM(avail_hist_fields(n)%vname)/='sig1' &
              .or.TRIM(avail_hist_fields(n)%vname)/='sig2' &
              .or.TRIM(avail_hist_fields(n)%vname)/='sistreave' &
              .or.TRIM(avail_hist_fields(n)%vname)/='sistremax' &
              .or.TRIM(avail_hist_fields(n)%vname)/='sigP') then
                status = pio_put_att(File,varid,'cell_methods','time: mean')
              endif
            endif

            if (histfreq(ns) == '1' .or. .not. hist_avg         &
                .or. n==n_divu(ns)      .or. n==n_shear(ns)     &  ! snapshots
                .or. n==n_sig1(ns)      .or. n==n_sig2(ns)      &
                .or. n==n_sigP(ns)      .or. n==n_trsig(ns)     &
                .or. n==n_sistreave(ns) .or. n==n_sistremax(ns) &
                .or. n==n_mlt_onset(ns) .or. n==n_frz_onset(ns) &
                .or. n==n_hisnap(ns)    .or. n==n_aisnap(ns)) then
               status = pio_put_att(File,varid,'time_rep','instantaneous')
            else
               status = pio_put_att(File,varid,'time_rep','averaged')
            endif
          endif
        enddo  ! num_avail_hist_fields_2D

      !-----------------------------------------------------------------
      ! 3D (category)
      !-----------------------------------------------------------------

        dimidz(1) = imtid
        dimidz(2) = jmtid
        dimidz(3) = cmtid
        dimidz(4) = timid

        do n = n2D + 1, n3Dccum
          if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            status  = pio_def_var(File, trim(avail_hist_fields(n)%vname), &
                         lprecision, dimidz, varid)
            status = pio_put_att(File,varid,'units', &
                        trim(avail_hist_fields(n)%vunit))
            status = pio_put_att(File,varid, 'long_name', &
                        trim(avail_hist_fields(n)%vdesc))
            status = pio_put_att(File,varid,'coordinates', &
                        trim(avail_hist_fields(n)%vcoord))
            status = pio_put_att(File,varid,'cell_measures', &
                        trim(avail_hist_fields(n)%vcellmeas))
            if (lprecision == pio_real) then
               status = pio_put_att(File, varid, 'missing_value', spval)
               status = pio_put_att(File, varid,'_FillValue',spval)
            else
               status = pio_put_att(File, varid, 'missing_value', spval_dbl)
               status = pio_put_att(File, varid,'_FillValue',spval_dbl)
            endif

            ! Add cell_methods attribute to variables if averaged
            if (hist_avg .and. histfreq(ns) /= '1') then
                status = pio_put_att(File,varid,'cell_methods','time: mean')
            endif

            if (histfreq(ns) == '1' .or. .not. hist_avg) then
               status = pio_put_att(File,varid,'time_rep','instantaneous')
            else
               status = pio_put_att(File,varid,'time_rep','averaged')
            endif
          endif
        enddo  ! num_avail_hist_fields_3Dc

      !-----------------------------------------------------------------
      ! 3D (ice layers)
      !-----------------------------------------------------------------

        dimidz(1) = imtid
        dimidz(2) = jmtid
        dimidz(3) = kmtidi
        dimidz(4) = timid

        do n = n3Dccum + 1, n3Dzcum
          if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            status  = pio_def_var(File, trim(avail_hist_fields(n)%vname), &
                         lprecision, dimidz, varid)
            status = pio_put_att(File,varid,'units', &
                        trim(avail_hist_fields(n)%vunit))
            status = pio_put_att(File,varid, 'long_name', &
                        trim(avail_hist_fields(n)%vdesc))
            status = pio_put_att(File,varid,'coordinates', &
                        trim(avail_hist_fields(n)%vcoord))
            status = pio_put_att(File,varid,'cell_measures', &
                        trim(avail_hist_fields(n)%vcellmeas))
            if (lprecision == pio_real) then
               status = pio_put_att(File, varid, 'missing_value', spval)
               status = pio_put_att(File, varid,'_FillValue',spval)
            else
               status = pio_put_att(File, varid, 'missing_value', spval_dbl)
               status = pio_put_att(File, varid,'_FillValue',spval_dbl)
            endif

            ! Add cell_methods attribute to variables if averaged
            if (hist_avg .and. histfreq(ns) /= '1') then
                status = pio_put_att(File,varid,'cell_methods','time: mean')
            endif

            if (histfreq(ns) == '1' .or. .not. hist_avg) then
               status = pio_put_att(File,varid,'time_rep','instantaneous')
            else
               status = pio_put_att(File,varid,'time_rep','averaged')
            endif
          endif
        enddo  ! num_avail_hist_fields_3Dz
        
      !-----------------------------------------------------------------
      ! 3D (biology ice layers)
      !-----------------------------------------------------------------

        dimidz(1) = imtid
        dimidz(2) = jmtid
        dimidz(3) = kmtidb
        dimidz(4) = timid

        do n = n3Dzcum + 1, n3Dbcum
          if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            status  = pio_def_var(File, trim(avail_hist_fields(n)%vname), &
                         lprecision, dimidz, varid)
            status = pio_put_att(File,varid,'units', &
                        trim(avail_hist_fields(n)%vunit))
            status = pio_put_att(File,varid, 'long_name', &
                        trim(avail_hist_fields(n)%vdesc))
            status = pio_put_att(File,varid,'coordinates', &
                        trim(avail_hist_fields(n)%vcoord))
            status = pio_put_att(File,varid,'cell_measures', &
                        trim(avail_hist_fields(n)%vcellmeas))
            if (lprecision == pio_real) then
               status = pio_put_att(File, varid, 'missing_value', spval)
               status = pio_put_att(File, varid,'_FillValue',spval)
            else
               status = pio_put_att(File, varid, 'missing_value', spval_dbl)
               status = pio_put_att(File, varid,'_FillValue',spval_dbl)
            endif

            ! Add cell_methods attribute to variables if averaged
            if (hist_avg .and. histfreq(ns) /= '1') then
                status = pio_put_att(File,varid,'cell_methods','time: mean')
            endif

            if (histfreq(ns) == '1' .or. .not. hist_avg) then
               status = pio_put_att(File,varid,'time_rep','instantaneous')
            else
               status = pio_put_att(File,varid,'time_rep','averaged')
            endif
          endif
        enddo  ! num_avail_hist_fields_3Db

      !-----------------------------------------------------------------
      ! 3D (biology snow layers)
      !-----------------------------------------------------------------

        dimidz(1) = imtid
        dimidz(2) = jmtid
        dimidz(3) = kmtida
        dimidz(4) = timid

        do n = n3Dbcum + 1, n3Dacum
          if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            status  = pio_def_var(File, trim(avail_hist_fields(n)%vname), &
                         lprecision, dimidz, varid)
            status = pio_put_att(File,varid,'units', &
                        trim(avail_hist_fields(n)%vunit))
            status = pio_put_att(File,varid, 'long_name', &
                        trim(avail_hist_fields(n)%vdesc))
            status = pio_put_att(File,varid,'coordinates', &
                        trim(avail_hist_fields(n)%vcoord))
            status = pio_put_att(File,varid,'cell_measures', &
                        trim(avail_hist_fields(n)%vcellmeas))
            if (lprecision == pio_real) then
               status = pio_put_att(File, varid, 'missing_value', spval)
               status = pio_put_att(File, varid,'_FillValue',spval)
            else
               status = pio_put_att(File, varid, 'missing_value', spval_dbl)
               status = pio_put_att(File, varid,'_FillValue',spval_dbl)
            endif

            ! Add cell_methods attribute to variables if averaged
            if (hist_avg .and. histfreq(ns) /= '1') then
                status = pio_put_att(File,varid,'cell_methods','time: mean')
            endif

            if (histfreq(ns) == '1' .or. .not. hist_avg) then
               status = pio_put_att(File,varid,'time_rep','instantaneous')
            else
               status = pio_put_att(File,varid,'time_rep','averaged')
            endif
          endif
        enddo  ! num_avail_hist_fields_3Da

      !-----------------------------------------------------------------
      ! 3D (fsd)
      !-----------------------------------------------------------------

        dimidz(1) = imtid
        dimidz(2) = jmtid
        dimidz(3) = fmtid
        dimidz(4) = timid

        do n = n3Dacum + 1, n3Dfcum
          if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            status  = pio_def_var(File, trim(avail_hist_fields(n)%vname), &
                         lprecision, dimidz, varid)
            status = pio_put_att(File,varid,'units', &
                        trim(avail_hist_fields(n)%vunit))
            status = pio_put_att(File,varid, 'long_name', &
                        trim(avail_hist_fields(n)%vdesc))
            status = pio_put_att(File,varid,'coordinates', &
                        trim(avail_hist_fields(n)%vcoord))
            status = pio_put_att(File,varid,'cell_measures', &
                        trim(avail_hist_fields(n)%vcellmeas))
            if (lprecision == pio_real) then
               status = pio_put_att(File, varid, 'missing_value', spval)
               status = pio_put_att(File, varid,'_FillValue',spval)
            else
               status = pio_put_att(File, varid, 'missing_value', spval_dbl)
               status = pio_put_att(File, varid,'_FillValue',spval_dbl)
            endif

            ! Add cell_methods attribute to variables if averaged
            if (hist_avg .and. histfreq(ns) /= '1') then
                status = pio_put_att(File,varid,'cell_methods','time: mean')
            endif

            if (histfreq(ns) == '1' .or. .not. hist_avg) then
               status = pio_put_att(File,varid,'time_rep','instantaneous')
            else
               status = pio_put_att(File,varid,'time_rep','averaged')
            endif
          endif
        enddo  ! num_avail_hist_fields_3Df

      !-----------------------------------------------------------------
      ! define attributes for 4D variables
      ! time coordinate is dropped
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! 4D (ice categories)
      !-----------------------------------------------------------------

        dimidcz(1) = imtid
        dimidcz(2) = jmtid
        dimidcz(3) = kmtidi
        dimidcz(4) = cmtid
        dimidcz(5) = timid

        do n = n3Dfcum + 1, n4Dicum
          if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            status  = pio_def_var(File, trim(avail_hist_fields(n)%vname), &
                             lprecision, dimidcz, varid)
            status = pio_put_att(File,varid,'units', &
                        trim(avail_hist_fields(n)%vunit))
            status = pio_put_att(File,varid, 'long_name', &
                        trim(avail_hist_fields(n)%vdesc))
            status = pio_put_att(File,varid,'coordinates', &
                        trim(avail_hist_fields(n)%vcoord))
            status = pio_put_att(File,varid,'cell_measures', &
                        trim(avail_hist_fields(n)%vcellmeas))
            if (lprecision == pio_real) then
               status = pio_put_att(File, varid, 'missing_value', spval)
               status = pio_put_att(File, varid,'_FillValue',spval)
            else
               status = pio_put_att(File, varid, 'missing_value', spval_dbl)
               status = pio_put_att(File, varid,'_FillValue',spval_dbl)
            endif

            ! Add cell_methods attribute to variables if averaged
            if (hist_avg .and. histfreq(ns) /= '1') then
                status = pio_put_att(File,varid,'cell_methods','time: mean')
            endif

            if (histfreq(ns) == '1' .or. .not. hist_avg) then
               status = pio_put_att(File,varid,'time_rep','instantaneous')
            else
               status = pio_put_att(File,varid,'time_rep','averaged')
            endif
          endif
        enddo  ! num_avail_hist_fields_4Di

      !-----------------------------------------------------------------
      ! 4D (snow layers)
      !-----------------------------------------------------------------

        dimidcz(1) = imtid
        dimidcz(2) = jmtid
        dimidcz(3) = kmtids
        dimidcz(4) = cmtid
        dimidcz(5) = timid

        do n = n4Dicum + 1, n4Dscum
          if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            status  = pio_def_var(File, trim(avail_hist_fields(n)%vname), &
                             lprecision, dimidcz, varid)
            status = pio_put_att(File,varid,'units', &
                        trim(avail_hist_fields(n)%vunit))
            status = pio_put_att(File,varid, 'long_name', &
                        trim(avail_hist_fields(n)%vdesc))
            status = pio_put_att(File,varid,'coordinates', &
                        trim(avail_hist_fields(n)%vcoord))
            status = pio_put_att(File,varid,'cell_measures', &
                        trim(avail_hist_fields(n)%vcellmeas))
            if (lprecision == pio_real) then
               status = pio_put_att(File, varid, 'missing_value', spval)
               status = pio_put_att(File, varid,'_FillValue',spval)
            else
               status = pio_put_att(File, varid, 'missing_value', spval_dbl)
               status = pio_put_att(File, varid,'_FillValue',spval_dbl)
            endif

            ! Add cell_methods attribute to variables if averaged
            if (hist_avg .and. histfreq(ns) /= '1') then
                status = pio_put_att(File,varid,'cell_methods','time: mean')
            endif

            if (histfreq(ns) == '1' .or. .not. hist_avg) then
               status = pio_put_att(File,varid,'time_rep','instantaneous')
            else
               status = pio_put_att(File,varid,'time_rep','averaged')
            endif
          endif
        enddo  ! num_avail_hist_fields_4Ds


      !-----------------------------------------------------------------
      ! 4D (fsd layers)
      !-----------------------------------------------------------------

        dimidcz(1) = imtid
        dimidcz(2) = jmtid
        dimidcz(3) = fmtid
        dimidcz(4) = cmtid
        dimidcz(5) = timid

        do n = n4Dscum + 1, n4Dfcum
          if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            status  = pio_def_var(File, trim(avail_hist_fields(n)%vname), &
                             lprecision, dimidcz, varid)
            status = pio_put_att(File,varid,'units', &
                        trim(avail_hist_fields(n)%vunit))
            status = pio_put_att(File,varid, 'long_name', &
                        trim(avail_hist_fields(n)%vdesc))
            status = pio_put_att(File,varid,'coordinates', &
                        trim(avail_hist_fields(n)%vcoord))
            status = pio_put_att(File,varid,'cell_measures', &
                        trim(avail_hist_fields(n)%vcellmeas))
            if (lprecision == pio_real) then
               status = pio_put_att(File, varid, 'missing_value', spval)
               status = pio_put_att(File, varid,'_FillValue',spval)
            else
               status = pio_put_att(File, varid, 'missing_value', spval_dbl)
               status = pio_put_att(File, varid,'_FillValue',spval_dbl)
            endif

            ! Add cell_methods attribute to variables if averaged
            if (hist_avg .and. histfreq(ns) /= '1') then
                status = pio_put_att(File,varid,'cell_methods','time: mean')
            endif

            if (histfreq(ns) == '1' .or. .not. hist_avg) then
               status = pio_put_att(File,varid,'time_rep','instantaneous')
            else
               status = pio_put_att(File,varid,'time_rep','averaged')
            endif
          endif
        enddo  ! num_avail_hist_fields_4Df

      !-----------------------------------------------------------------
      ! global attributes
      !-----------------------------------------------------------------
      ! ... the user should change these to something useful ...
      !-----------------------------------------------------------------
#ifdef CESMCOUPLED
        status = pio_put_att(File,pio_global,'title',runid)
#else
        title  = 'sea ice model output for CICE'
        status = pio_put_att(File,pio_global,'title',trim(title))
#endif
        title = 'Diagnostic and Prognostic Variables'
        status = pio_put_att(File,pio_global,'contents',trim(title))

        write(title,'(2a)') 'Los Alamos Sea Ice Model, ', trim(version_name)
        status = pio_put_att(File,pio_global,'source',trim(title))

        if (use_leap_years) then
          write(title,'(a,i3,a)') 'This year has ',dayyr,' days'
        else
          write(title,'(a,i3,a)') 'All years have exactly ',dayyr,' days'
        endif
        status = pio_put_att(File,pio_global,'comment',trim(title))

        write(title,'(a,i8.8)') 'File written on model date ',idate
        status = pio_put_att(File,pio_global,'comment2',trim(title))

        write(title,'(a,i6)') 'seconds elapsed into model date: ',msec
        status = pio_put_att(File,pio_global,'comment3',trim(title))

        title = 'CF-1.0'
        status =  &
             pio_put_att(File,pio_global,'conventions',trim(title))

        call date_and_time(date=current_date, time=current_time)
        write(start_time,1000) current_date(1:4), current_date(5:6), &
                               current_date(7:8), current_time(1:2), &
                               current_time(3:4)
1000    format('This dataset was created on ', &
                a,'-',a,'-',a,' at ',a,':',a)
        status = pio_put_att(File,pio_global,'history',trim(start_time))

        if (history_format == 'pio_pnetcdf') then
           status = pio_put_att(File,pio_global,'io_flavor','io_pio pnetcdf')
        else
           status = pio_put_att(File,pio_global,'io_flavor','io_pio netcdf')
        endif

      !-----------------------------------------------------------------
      ! end define mode
      !-----------------------------------------------------------------

        status = pio_enddef(File)

      !-----------------------------------------------------------------
      ! write time variable
      !-----------------------------------------------------------------

        status = pio_inq_varid(File,'time',varid)
!sgl        status = pio_put_var(File,varid,(/1/),ltime)
        status = pio_put_var(File,varid,(/1/),ltime2)

      !-----------------------------------------------------------------
      ! write time_bounds info
      !-----------------------------------------------------------------

        if (hist_avg .and. histfreq(ns) /= '1') then
          status = pio_inq_varid(File,'time_bounds',varid)
          time_bounds=(/time_beg(ns),time_end(ns)/)
          bnd_start  = (/1,1/)
          bnd_length = (/2,1/)
          status = pio_put_var(File,varid,ival=time_bounds, &
                   start=bnd_start(:),count=bnd_length(:)) 
        endif

      !-----------------------------------------------------------------
      ! write coordinate variables
      !-----------------------------------------------------------------

        allocate(workr2(nx_block,ny_block,nblocks))

        do i = 1,ncoord
          status = pio_inq_varid(File, coord_var(i)%short_name, varid)
          SELECT CASE (coord_var(i)%short_name)
            CASE ('TLON')
              ! Convert T grid longitude from -180 -> 180 to 0 to 360
                 workr2(:,:,:) = mod(tlon(:,:,1:nblocks)*rad_to_deg + c360, c360)
            CASE ('TLAT')
              workr2(:,:,:) = tlat(:,:,1:nblocks)*rad_to_deg
            CASE ('ULON')
              workr2(:,:,:) = ulon(:,:,1:nblocks)*rad_to_deg
            CASE ('ULAT')
              workr2(:,:,:) = ulat(:,:,1:nblocks)*rad_to_deg
          END SELECT
          call pio_write_darray(File, varid, iodesc2d, &
               workr2, status, fillval=spval_dbl)
        enddo

        ! Extra dimensions (NCAT, NFSD, VGRD*)

        do i = 1, nvarz
          if (igrdz(i)) then
            status = pio_inq_varid(File, var_nz(i)%short_name, varid)
            SELECT CASE (var_nz(i)%short_name)
              CASE ('NCAT')
                status = pio_put_var(File, varid, hin_max(1:ncat_hist)) 
              CASE ('NFSD')
                status = pio_put_var(File, varid, floe_rad_c(1:nfsd_hist))
              CASE ('VGRDi')
                status = pio_put_var(File, varid, (/(k, k=1,nzilyr)/))
              CASE ('VGRDs')
                status = pio_put_var(File, varid, (/(k, k=1,nzslyr)/))
              CASE ('VGRDb')
                status = pio_put_var(File, varid, (/(k, k=1,nzblyr)/))
              CASE ('VGRDa')
                status = pio_put_var(File, varid, (/(k, k=1,nzalyr)/))
             END SELECT
           endif
        enddo

      !-----------------------------------------------------------------
      ! write grid masks, area and rotation angle
      !-----------------------------------------------------------------

!      if (igrd(n_tmask)) then
!        status = pio_inq_varid(File, 'tmask', varid)
!        call pio_write_darray(File, varid, iodesc2d, &
!                              hm(:,:,1:nblocks), status, fillval=spval_dbl)
!      endif
!      if (igrd(n_blkmask)) then
!        status = pio_inq_varid(File, 'blkmask', varid)
!        call pio_write_darray(File, varid, iodesc2d, &
!                              bm(:,:,1:nblocks), status, fillval=spval_dbl)
!      endif
        
      do i = 1, nvar       ! note: n_tmask=1, n_blkmask=2
         if (igrd(i)) then
            SELECT CASE (var(i)%req%short_name)
            CASE ('tmask')
               workr2 = hm(:,:,1:nblocks)
            CASE ('blkmask')
               workr2 = bm(:,:,1:nblocks)
            CASE ('tarea')
               workr2 = tarea(:,:,1:nblocks)
            CASE ('uarea')
               workr2 = uarea(:,:,1:nblocks)
            CASE ('dxu')
               workr2 = dxu(:,:,1:nblocks)
            CASE ('dyu')
               workr2 = dyu(:,:,1:nblocks)
            CASE ('dxt')
               workr2 = dxt(:,:,1:nblocks)
            CASE ('dyt')
               workr2 = dyt(:,:,1:nblocks)
            CASE ('HTN')
               workr2 = HTN(:,:,1:nblocks)
            CASE ('HTE')
               workr2 = HTE(:,:,1:nblocks)
            CASE ('ANGLE')
               workr2 = ANGLE(:,:,1:nblocks)
            CASE ('ANGLET')
               workr2 = ANGLET(:,:,1:nblocks)
            END SELECT
            status = pio_inq_varid(File, var(i)%req%short_name, varid)
            call pio_write_darray(File, varid, iodesc2d, &
                 workr2, status, fillval=spval_dbl)
         endif
      enddo

      !----------------------------------------------------------------
      ! Write coordinates of grid box vertices
      !----------------------------------------------------------------

      if (f_bounds) then
      allocate(workr3v(nverts,nx_block,ny_block,nblocks))
      workr3v (:,:,:,:) = c0
      do i = 1, nvar_verts
        SELECT CASE (var_nverts(i)%short_name)
        CASE ('lont_bounds')
           do ivertex = 1, nverts 
              workr3v(ivertex,:,:,:) = lont_bounds(ivertex,:,:,1:nblocks)
           enddo
        CASE ('latt_bounds')
           do ivertex = 1, nverts 
              workr3v(ivertex,:,:,:) = latt_bounds(ivertex,:,:,1:nblocks)
           enddo
        CASE ('lonu_bounds')
           do ivertex = 1, nverts 
              workr3v(ivertex,:,:,:) = lonu_bounds(ivertex,:,:,1:nblocks)
           enddo
        CASE ('latu_bounds')
           do ivertex = 1, nverts 
              workr3v(ivertex,:,:,:) = latu_bounds(ivertex,:,:,1:nblocks)
           enddo
        END SELECT

          status = pio_inq_varid(File, var_nverts(i)%short_name, varid)
          call pio_write_darray(File, varid, iodesc3dv, &
                                workr3v, status, fillval=spval_dbl)
      enddo
      deallocate(workr3v)
      endif  ! f_bounds


      !-----------------------------------------------------------------
      ! write variable data
      !-----------------------------------------------------------------

      ! 2D
      do n=1,num_avail_hist_fields_2D
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            status  = pio_inq_varid(File,avail_hist_fields(n)%vname,varid)
            if (status /= pio_noerr) call abort_ice(subname// &
               'ERROR getting varid for '//avail_hist_fields(n)%vname)
            workr2(:,:,:) = a2D(:,:,n,1:nblocks)
#ifdef CESM1_PIO
            call pio_setframe(varid, int(1,kind=PIO_OFFSET))
#else
            call pio_setframe(File, varid, int(1,kind=PIO_OFFSET_KIND))
#endif
            call pio_write_darray(File, varid, iodesc2d,&
                                  workr2, status, fillval=spval_dbl)
         endif
      enddo ! num_avail_hist_fields_2D

      deallocate(workr2)

      ! 3D (category)
      allocate(workr3(nx_block,ny_block,nblocks,ncat_hist))
      do n = n2D + 1, n3Dccum
         nn = n - n2D
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            status  = pio_inq_varid(File,avail_hist_fields(n)%vname,varid)
            if (status /= pio_noerr) call abort_ice(subname// &
               'ERROR: getting varid for '//avail_hist_fields(n)%vname)
            do j = 1, nblocks
            do i = 1, ncat_hist
               workr3(:,:,j,i) = a3Dc(:,:,i,nn,j)
            enddo
            enddo
#ifdef CESM1_PIO
            call pio_setframe(varid, int(1,kind=PIO_OFFSET))
#else
            call pio_setframe(File, varid, int(1,kind=PIO_OFFSET_KIND))
#endif
            call pio_write_darray(File, varid, iodesc3dc,&
                                  workr3, status, fillval=spval_dbl)
         endif
      enddo ! num_avail_hist_fields_3Dc
      deallocate(workr3)

      ! 3D (vertical ice)
      allocate(workr3(nx_block,ny_block,nblocks,nzilyr))
      do n = n3Dccum+1, n3Dzcum
         nn = n - n3Dccum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            status  = pio_inq_varid(File,avail_hist_fields(n)%vname,varid)
            if (status /= pio_noerr) call abort_ice(subname// &
               'ERROR: getting varid for '//avail_hist_fields(n)%vname)
            do j = 1, nblocks
            do i = 1, nzilyr
               workr3(:,:,j,i) = a3Dz(:,:,i,nn,j)
            enddo
            enddo
#ifdef CESM1_PIO
            call pio_setframe(varid, int(1,kind=PIO_OFFSET))
#else
            call pio_setframe(File, varid, int(1,kind=PIO_OFFSET_KIND))
#endif
            call pio_write_darray(File, varid, iodesc3di,&
                                  workr3, status, fillval=spval_dbl)
         endif
      enddo ! num_avail_hist_fields_3Dz
      deallocate(workr3)

      ! 3D (vertical ice biology)
      allocate(workr3(nx_block,ny_block,nblocks,nzblyr))
      do n = n3Dzcum+1, n3Dbcum
         nn = n - n3Dzcum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            status  = pio_inq_varid(File,avail_hist_fields(n)%vname,varid)
            if (status /= pio_noerr) call abort_ice(subname// &
               'ERROR: getting varid for '//avail_hist_fields(n)%vname)
            do j = 1, nblocks
            do i = 1, nzblyr
               workr3(:,:,j,i) = a3Db(:,:,i,nn,j)
            enddo
            enddo
#ifdef CESM1_PIO
            call pio_setframe(varid, int(1,kind=PIO_OFFSET))
#else
            call pio_setframe(File, varid, int(1,kind=PIO_OFFSET_KIND))
#endif
            call pio_write_darray(File, varid, iodesc3db,&
                                  workr3, status, fillval=spval_dbl)
         endif
      enddo ! num_avail_hist_fields_3Db
      deallocate(workr3)

      ! 3D (vertical snow biology)
      allocate(workr3(nx_block,ny_block,nblocks,nzalyr))
      do n = n3Dbcum+1, n3Dacum
         nn = n - n3Dbcum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            status  = pio_inq_varid(File,avail_hist_fields(n)%vname,varid)
            if (status /= pio_noerr) call abort_ice(subname// &
               'ERROR: getting varid for '//avail_hist_fields(n)%vname)
            do j = 1, nblocks
            do i = 1, nzalyr
               workr3(:,:,j,i) = a3Da(:,:,i,nn,j)
            enddo
            enddo
#ifdef CESM1_PIO
            call pio_setframe(varid, int(1,kind=PIO_OFFSET))
#else
            call pio_setframe(File, varid, int(1,kind=PIO_OFFSET_KIND))
#endif
            call pio_write_darray(File, varid, iodesc3da,&
                                  workr3, status, fillval=spval_dbl)
         endif
      enddo ! num_avail_hist_fields_3Db
      deallocate(workr3)

      ! 3D (fsd)
      allocate(workr3(nx_block,ny_block,nblocks,nfsd_hist))
      do n = n3Dacum+1, n3Dfcum
         nn = n - n3Dacum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            status  = pio_inq_varid(File,avail_hist_fields(n)%vname,varid)
            if (status /= pio_noerr) call abort_ice(subname// &
               'ERROR: getting varid for '//avail_hist_fields(n)%vname)
            do j = 1, nblocks
            do i = 1, nfsd_hist
               workr3(:,:,j,i) = a3Df(:,:,i,nn,j)
            enddo
            enddo
#ifdef CESM1_PIO
            call pio_setframe(varid, int(1,kind=PIO_OFFSET))
#else
            call pio_setframe(File, varid, int(1,kind=PIO_OFFSET_KIND))
#endif
            call pio_write_darray(File, varid, iodesc3df,&
                                  workr3, status, fillval=spval_dbl)
         endif
      enddo ! num_avail_hist_fields_3Df
      deallocate(workr3)

      allocate(workr4(nx_block,ny_block,nblocks,ncat_hist,nzilyr))
      ! 4D (categories, fsd)
      do n = n3Dfcum+1, n4Dicum
         nn = n - n3Dfcum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            status  = pio_inq_varid(File,avail_hist_fields(n)%vname,varid)
            if (status /= pio_noerr) call abort_ice(subname// &
               'ERROR: getting varid for '//avail_hist_fields(n)%vname)
            do j = 1, nblocks
            do i = 1, ncat_hist
            do k = 1, nzilyr
               workr4(:,:,j,i,k) = a4Di(:,:,k,i,nn,j)
            enddo ! k
            enddo ! i
            enddo ! j
#ifdef CESM1_PIO
            call pio_setframe(varid, int(1,kind=PIO_OFFSET))
#else
            call pio_setframe(File, varid, int(1,kind=PIO_OFFSET_KIND))
#endif
            call pio_write_darray(File, varid, iodesc4di,&
                                  workr4, status, fillval=spval_dbl)
         endif
      enddo ! num_avail_hist_fields_4Di
      deallocate(workr4)

      allocate(workr4(nx_block,ny_block,nblocks,ncat_hist,nzslyr))
      ! 4D (categories, vertical ice)
      do n = n4Dicum+1, n4Dscum
         nn = n - n4Dicum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            status  = pio_inq_varid(File,avail_hist_fields(n)%vname,varid)
            if (status /= pio_noerr) call abort_ice(subname// &
               'ERROR: getting varid for '//avail_hist_fields(n)%vname)
            do j = 1, nblocks
            do i = 1, ncat_hist
            do k = 1, nzslyr
               workr4(:,:,j,i,k) = a4Ds(:,:,k,i,nn,j)
            enddo ! k
            enddo ! i
            enddo ! j
#ifdef CESM1_PIO
            call pio_setframe(varid, int(1,kind=PIO_OFFSET))
#else
            call pio_setframe(File, varid, int(1,kind=PIO_OFFSET_KIND))
#endif
            call pio_write_darray(File, varid, iodesc4ds,&
                                  workr4, status, fillval=spval_dbl)
         endif
      enddo ! num_avail_hist_fields_4Ds
      deallocate(workr4)

      allocate(workr4(nx_block,ny_block,nblocks,ncat_hist,nfsd_hist))
      ! 4D (categories, vertical ice)
      do n = n4Dscum+1, n4Dfcum
         nn = n - n4Dscum
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            status  = pio_inq_varid(File,avail_hist_fields(n)%vname,varid)
            if (status /= pio_noerr) call abort_ice(subname// &
               'ERROR: getting varid for '//avail_hist_fields(n)%vname)
            do j = 1, nblocks
            do i = 1, ncat_hist
            do k = 1, nfsd_hist
               workr4(:,:,j,i,k) = a4Df(:,:,k,i,nn,j)
            enddo ! k
            enddo ! i
            enddo ! j
#ifdef CESM1_PIO
            call pio_setframe(varid, int(1,kind=PIO_OFFSET))
#else
            call pio_setframe(File, varid, int(1,kind=PIO_OFFSET_KIND))
#endif
            call pio_write_darray(File, varid, iodesc4df,&
                                  workr4, status, fillval=spval_dbl)
         endif
      enddo ! num_avail_hist_fields_4Df
      deallocate(workr4)

!     similarly for num_avail_hist_fields_4Db (define workr4b, iodesc4db)


      !-----------------------------------------------------------------
      ! clean-up PIO descriptors
      !-----------------------------------------------------------------

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
         write(nu_diag,*) ' '
         write(nu_diag,*) 'Finished writing ',trim(ncfile(ns))
      endif

      end subroutine ice_write_hist

!=======================================================================

      end module ice_history_write

!=======================================================================
