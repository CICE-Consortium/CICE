#ifdef ncdf
#define USE_NETCDF
#endif
!=======================================================================

! Read and write ice model restart files using netCDF or binary
! interfaces.
! authors David A Bailey, NCAR

      module ice_restart

      use ice_broadcast
      use ice_constants, only: c0
      use ice_communicate, only: my_task, master_task
      use ice_kinds_mod
#ifdef USE_NETCDF
      use netcdf
#endif
      use ice_read_write, only: ice_check_nc
      use ice_restart_shared
      use ice_fileunits, only: nu_diag, nu_rst_pointer
      use ice_exit, only: abort_ice
      use icepack_intfc, only: icepack_query_parameters
      use icepack_intfc, only: icepack_query_tracer_sizes
      use icepack_intfc, only: icepack_query_tracer_flags
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted

      implicit none
      private
      public :: init_restart_write, init_restart_read, &
         read_restart_field, write_restart_field, final_restart, &
         query_field

      integer (kind=int_kind) :: ncid , &
         dimid_ni,   & ! netCDF identifiers
         dimid_nj

!=======================================================================

      contains

!=======================================================================

! Sets up restart file for reading.
! author David A Bailey, NCAR

      subroutine init_restart_read(ice_ic)

      use ice_calendar, only: msec, mmonth, mday, myear, &
                             istep0, istep1, npt

      character(len=char_len_long), intent(in), optional :: ice_ic

      ! local variables

      character(len=char_len_long) :: &
         filename, filename0

      integer (kind=int_kind) :: status

      character(len=*), parameter :: subname = '(init_restart_read)'

#ifdef USE_NETCDF
      if (present(ice_ic)) then
         filename = trim(ice_ic)
      else
         if (my_task == master_task) then
            open(nu_rst_pointer,file=pointer_file, status='old')
            read(nu_rst_pointer,'(a)') filename0
            filename = trim(filename0)
            close(nu_rst_pointer)
            write(nu_diag,*) 'Read ',pointer_file(1:lenstr(pointer_file))
         endif
         call broadcast_scalar(filename, master_task)
      endif

      if (my_task == master_task) then
         write(nu_diag,*) 'Using restart dump=', trim(filename)

         status = nf90_open(trim(filename), nf90_nowrite, ncid)
         call ice_check_nc(status, subname//' ERROR: open '//trim(filename), file=__FILE__, line=__LINE__)

         if (use_restart_time) then
            ! for backwards compatibility, check nyr, month, and sec as well
            status = nf90_get_att(ncid, nf90_global, 'istep1', istep0)
            call ice_check_nc(status, subname//" ERROR: reading restart step ",file=__FILE__,line=__LINE__)

            status = nf90_get_att(ncid, nf90_global, 'myear', myear)
            if (status /= nf90_noerr) then
               status = nf90_get_att(ncid, nf90_global, 'nyr', myear)
               call ice_check_nc(status, subname//" ERROR: reading restart year ",file=__FILE__,line=__LINE__)
            endif

            status = nf90_get_att(ncid, nf90_global, 'mmonth', mmonth)
            if (status /= nf90_noerr) then
               status = nf90_get_att(ncid, nf90_global, 'month', mmonth)
               call ice_check_nc(status, subname//" ERROR: reading restart month ",file=__FILE__,line=__LINE__)
            endif

            status = nf90_get_att(ncid, nf90_global, 'mday', mday)
            call ice_check_nc(status, subname//" ERROR: reading restart day ",file=__FILE__,line=__LINE__)

            status = nf90_get_att(ncid, nf90_global, 'msec', msec)
            if (status /= nf90_noerr) then
               status = nf90_get_att(ncid, nf90_global, 'sec', msec)
               call ice_check_nc(status, subname//" ERROR: reading restart sec ",file=__FILE__,line=__LINE__)
            endif

         endif ! use namelist values if use_restart_time = F

      endif

      call broadcast_scalar(istep0,master_task)
      call broadcast_scalar(myear,master_task)
      call broadcast_scalar(mmonth,master_task)
      call broadcast_scalar(mday,master_task)
      call broadcast_scalar(msec,master_task)

      istep1 = istep0

      ! if runid is bering then need to correct npt for istep0
      if (trim(runid) == 'bering') then
         npt = npt - istep0
      endif
#else
      call abort_ice(subname//' ERROR: USE_NETCDF cpp not defined for '//trim(ice_ic), &
          file=__FILE__, line=__LINE__)
#endif

      end subroutine init_restart_read

!=======================================================================

! Sets up restart file for writing.
! author David A Bailey, NCAR

      subroutine init_restart_write(filename_spec)

      use ice_blocks, only: nghost
      use ice_calendar, only: msec, mmonth, mday, myear, istep1
      use ice_domain_size, only: nx_global, ny_global, ncat, nilyr, nslyr, &
                                 n_iso, n_aero, nblyr, n_zaero, n_algae, n_doc,   &
                                 n_dic, n_don, n_fed, n_fep, nfsd
      use ice_arrays_column, only: oceanmixed_ice
      use ice_dyn_shared, only: kdyn
      use ice_grid, only: grid_ice

      character(len=*), intent(in), optional :: filename_spec

      ! local variables

      logical (kind=log_kind) :: &
         skl_bgc, z_tracers, tr_fsd, &
         tr_iage, tr_FY, tr_lvl, tr_iso, tr_aero, &
         tr_pond_topo, tr_pond_lvl, tr_pond_sealvl, tr_brine, tr_snow, &
         tr_bgc_N, tr_bgc_C, tr_bgc_Nit, &
         tr_bgc_Sil, tr_bgc_DMS, &
         tr_bgc_chl, tr_bgc_Am,  &
         tr_bgc_PON, tr_bgc_DON, &
         tr_zaero,   tr_bgc_Fe,  &
         tr_bgc_hum

      integer (kind=int_kind) :: &
         k, n,                 & ! index
         nx, ny,               & ! global array size
         nbtrcr                  ! number of bgc tracers

      character(len=char_len_long) :: filename
      character(len=char_len_long) :: lpointer_file

      integer (kind=int_kind), allocatable :: dims(:)

      integer (kind=int_kind) :: &

        dimid_ncat, & !
        iflag,      & ! netCDF creation flag
        status        ! status variable from netCDF routine

      character (len=3) :: nchar, ncharb

      character(len=*), parameter :: subname = '(init_restart_write)'

#ifdef USE_NETCDF
      call icepack_query_parameters( &
         skl_bgc_out=skl_bgc, z_tracers_out=z_tracers)
      call icepack_query_tracer_sizes( &
         nbtrcr_out=nbtrcr)
      call icepack_query_tracer_flags( &
         tr_iage_out=tr_iage, tr_FY_out=tr_FY, tr_lvl_out=tr_lvl, tr_fsd_out=tr_fsd, &
         tr_iso_out=tr_iso, tr_aero_out=tr_aero, &
         tr_pond_topo_out=tr_pond_topo, tr_pond_lvl_out=tr_pond_lvl, &
         tr_pond_sealvl_out=tr_pond_sealvl, &
         tr_snow_out=tr_snow, tr_brine_out=tr_brine, &
         tr_bgc_N_out=tr_bgc_N, tr_bgc_C_out=tr_bgc_C, tr_bgc_Nit_out=tr_bgc_Nit, &
         tr_bgc_Sil_out=tr_bgc_Sil, tr_bgc_DMS_out=tr_bgc_DMS, &
         tr_bgc_chl_out=tr_bgc_chl, tr_bgc_Am_out=tr_bgc_Am, &
         tr_bgc_PON_out=tr_bgc_PON, tr_bgc_DON_out=tr_bgc_DON, &
         tr_zaero_out=tr_zaero,   tr_bgc_Fe_out=tr_bgc_Fe, &
         tr_bgc_hum_out=tr_bgc_hum)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      ! construct path/file
      if (present(filename_spec)) then
         filename = trim(filename_spec)
      else
         write(filename,'(a,a,a,i4.4,a,i2.2,a,i2.2,a,i5.5)') &
              restart_dir(1:lenstr(restart_dir)), &
              restart_file(1:lenstr(restart_file)),'.', &
              myear,'-',mmonth,'-',mday,'-',msec
      end if

      ! write pointer (path/file)
      if (my_task == master_task) then
         filename = trim(filename) // '.nc'
         lpointer_file = pointer_file
         if (pointer_date) then
            ! append date to pointer filename
            write(lpointer_file,'(a,i4.4,a,i2.2,a,i2.2,a,i5.5)') &
               trim(lpointer_file)//'.',myear,'-',mmonth,'-',mday,'-',msec
         end if
         open(nu_rst_pointer,file=lpointer_file)
         write(nu_rst_pointer,'(a)') filename
         close(nu_rst_pointer)

         if (restart_format == 'cdf1') then
           iflag = nf90_clobber
         elseif (restart_format == 'cdf2') then
#ifdef NO_CDF2
           call abort_ice(subname//' ERROR: restart_format cdf2 not available ', &
              file=__FILE__, line=__LINE__)
#else
           iflag = ior(nf90_clobber,nf90_64bit_offset)
#endif
         elseif (restart_format == 'cdf5') then
#ifdef NO_CDF5
           call abort_ice(subname//' ERROR: restart_format cdf5 not available ', &
              file=__FILE__, line=__LINE__)
#else
           iflag = ior(nf90_clobber,nf90_64bit_data)
#endif
         elseif (restart_format == 'hdf5') then
#ifdef NO_HDF5
           call abort_ice(subname//' ERROR: restart_format hdf5 not available ', &
              file=__FILE__, line=__LINE__)
#else
           iflag = ior(nf90_clobber,nf90_netcdf4)
#endif
         else
           call abort_ice(subname//' ERROR: restart_format not allowed for '//trim(restart_format), &
              file=__FILE__, line=__LINE__)
         endif
         status = nf90_create(trim(filename), iflag, ncid)
         call ice_check_nc(status, subname//' ERROR: creating '//trim(filename), file=__FILE__, line=__LINE__)

         status = nf90_put_att(ncid,nf90_global,'istep1',istep1)
         call ice_check_nc(status, subname//' ERROR: writing att istep', file=__FILE__, line=__LINE__)
         status = nf90_put_att(ncid,nf90_global,'myear',myear)
         call ice_check_nc(status, subname//' ERROR: writing att year', file=__FILE__, line=__LINE__)
         status = nf90_put_att(ncid,nf90_global,'mmonth',mmonth)
         call ice_check_nc(status, subname//' ERROR: writing att month', file=__FILE__, line=__LINE__)
         status = nf90_put_att(ncid,nf90_global,'mday',mday)
         call ice_check_nc(status, subname//' ERROR: writing att day', file=__FILE__, line=__LINE__)
         status = nf90_put_att(ncid,nf90_global,'msec',msec)
         call ice_check_nc(status, subname//' ERROR: writing att sec', file=__FILE__, line=__LINE__)

         nx = nx_global
         ny = ny_global
         if (restart_ext) then
            nx = nx_global + 2*nghost
            ny = ny_global + 2*nghost
         endif
         status = nf90_def_dim(ncid,'ni',nx,dimid_ni)
         call ice_check_nc(status, subname//' ERROR: writing dim ni', file=__FILE__, line=__LINE__)
         status = nf90_def_dim(ncid,'nj',ny,dimid_nj)
         call ice_check_nc(status, subname//' ERROR: writing dim nj', file=__FILE__, line=__LINE__)

         status = nf90_def_dim(ncid,'ncat',ncat,dimid_ncat)
         call ice_check_nc(status, subname//' ERROR: writing dim ncat', file=__FILE__, line=__LINE__)

         !-----------------------------------------------------------------
         ! 2D restart fields
         !-----------------------------------------------------------------

         allocate(dims(2))

         dims(1) = dimid_ni
         dims(2) = dimid_nj

         call define_rest_field(ncid,'uvel',dims)
         call define_rest_field(ncid,'vvel',dims)

         if (grid_ice == 'CD') then
            call define_rest_field(ncid,'uvelE',dims)
            call define_rest_field(ncid,'vvelE',dims)
            call define_rest_field(ncid,'uvelN',dims)
            call define_rest_field(ncid,'vvelN',dims)
         endif

         if (grid_ice == 'C') then
            call define_rest_field(ncid,'uvelE',dims)
            call define_rest_field(ncid,'vvelN',dims)
         endif

         if (restart_coszen) call define_rest_field(ncid,'coszen',dims)

         call define_rest_field(ncid,'scale_factor',dims)
         call define_rest_field(ncid,'swvdr',dims)
         call define_rest_field(ncid,'swvdf',dims)
         call define_rest_field(ncid,'swidr',dims)
         call define_rest_field(ncid,'swidf',dims)

         call define_rest_field(ncid,'strocnxT',dims)
         call define_rest_field(ncid,'strocnyT',dims)

         call define_rest_field(ncid,'stressp_1',dims)
         call define_rest_field(ncid,'stressp_2',dims)
         call define_rest_field(ncid,'stressp_3',dims)
         call define_rest_field(ncid,'stressp_4',dims)

         call define_rest_field(ncid,'stressm_1',dims)
         call define_rest_field(ncid,'stressm_2',dims)
         call define_rest_field(ncid,'stressm_3',dims)
         call define_rest_field(ncid,'stressm_4',dims)

         call define_rest_field(ncid,'stress12_1',dims)
         call define_rest_field(ncid,'stress12_2',dims)
         call define_rest_field(ncid,'stress12_3',dims)
         call define_rest_field(ncid,'stress12_4',dims)

         call define_rest_field(ncid,'iceumask',dims)

         if (grid_ice == 'CD' .or. grid_ice == 'C') then
            call define_rest_field(ncid,'stresspT' ,dims)
            call define_rest_field(ncid,'stressmT' ,dims)
            call define_rest_field(ncid,'stress12T',dims)
            call define_rest_field(ncid,'stresspU' ,dims)
            call define_rest_field(ncid,'stressmU' ,dims)
            call define_rest_field(ncid,'stress12U',dims)
            call define_rest_field(ncid,'icenmask',dims)
            call define_rest_field(ncid,'iceemask',dims)
         endif


         if (oceanmixed_ice) then
            call define_rest_field(ncid,'sst',dims)
            call define_rest_field(ncid,'frzmlt',dims)
         endif

         if (tr_FY) then
            call define_rest_field(ncid,'frz_onset',dims)
         endif

         if (kdyn == 2) then
            call define_rest_field(ncid,'a11_1',dims)
            call define_rest_field(ncid,'a11_2',dims)
            call define_rest_field(ncid,'a11_3',dims)
            call define_rest_field(ncid,'a11_4',dims)
            call define_rest_field(ncid,'a12_1',dims)
            call define_rest_field(ncid,'a12_2',dims)
            call define_rest_field(ncid,'a12_3',dims)
            call define_rest_field(ncid,'a12_4',dims)
         endif

         if (tr_pond_lvl .or. tr_pond_sealvl) then
            call define_rest_field(ncid,'fsnow',dims)
         endif

         if (nbtrcr > 0) then
            if (tr_bgc_N) then
            do k=1,n_algae
               write(nchar,'(i3.3)') k
               call define_rest_field(ncid,'algalN'//trim(nchar),dims)
            enddo
            endif
            if (tr_bgc_C) then
            do k=1,n_doc
               write(nchar,'(i3.3)') k
               call define_rest_field(ncid,'doc'//trim(nchar),dims)
            enddo
            do k=1,n_dic
               write(nchar,'(i3.3)') k
               call define_rest_field(ncid,'dic'//trim(nchar),dims)
            enddo
            endif
            call define_rest_field(ncid,'nit'   ,dims)
            if (tr_bgc_Am) &
            call define_rest_field(ncid,'amm'   ,dims)
            if (tr_bgc_Sil) &
            call define_rest_field(ncid,'sil'   ,dims)
            if (tr_bgc_hum) &
            call define_rest_field(ncid,'hum'   ,dims)
            if (tr_bgc_DMS) then
              call define_rest_field(ncid,'dmsp'  ,dims)
              call define_rest_field(ncid,'dms'   ,dims)
            endif
            if (tr_bgc_DON) then
            do k=1,n_don
               write(nchar,'(i3.3)') k
               call define_rest_field(ncid,'don'//trim(nchar),dims)
            enddo
            endif
            if (tr_bgc_Fe ) then
            do k=1,n_fed
               write(nchar,'(i3.3)') k
               call define_rest_field(ncid,'fed'//trim(nchar),dims)
            enddo
            do k=1,n_fep
               write(nchar,'(i3.3)') k
               call define_rest_field(ncid,'fep'//trim(nchar),dims)
            enddo
            endif
            if (tr_zaero) then
            do k=1,n_zaero
               write(nchar,'(i3.3)') k
               call define_rest_field(ncid,'zaeros'//trim(nchar),dims)
            enddo
            endif
         endif  !nbtrcr

         deallocate(dims)

         !-----------------------------------------------------------------
         ! 3D restart fields (ncat)
         !-----------------------------------------------------------------

         allocate(dims(3))

         dims(1) = dimid_ni
         dims(2) = dimid_nj
         dims(3) = dimid_ncat

         call define_rest_field(ncid,'aicen',dims)
         call define_rest_field(ncid,'vicen',dims)
         call define_rest_field(ncid,'vsnon',dims)
         call define_rest_field(ncid,'Tsfcn',dims)

         if (tr_iage) then
            call define_rest_field(ncid,'iage',dims)
         end if

         if (tr_FY) then
            call define_rest_field(ncid,'FY',dims)
         end if

         if (tr_lvl) then
            call define_rest_field(ncid,'alvl',dims)
            call define_rest_field(ncid,'vlvl',dims)
         end if

         if (tr_pond_topo) then
            call define_rest_field(ncid,'apnd',dims)
            call define_rest_field(ncid,'hpnd',dims)
            call define_rest_field(ncid,'ipnd',dims)
         end if

         if (tr_pond_lvl .or. tr_pond_sealvl) then
            call define_rest_field(ncid,'apnd',dims)
            call define_rest_field(ncid,'hpnd',dims)
            call define_rest_field(ncid,'ipnd',dims)
            call define_rest_field(ncid,'dhs',dims)
            call define_rest_field(ncid,'ffrac',dims)
         end if

         if (tr_brine) then
            call define_rest_field(ncid,'fbrn',dims)
            call define_rest_field(ncid,'first_ice',dims)
         endif

         if (skl_bgc) then
            do k = 1, n_algae
               write(nchar,'(i3.3)') k
               call define_rest_field(ncid,'bgc_N'//trim(nchar)    ,dims)
            enddo
            if (tr_bgc_C) then
            !  do k = 1, n_algae
            !     write(nchar,'(i3.3)') k
            !     call define_rest_field(ncid,'bgc_C'//trim(nchar)    ,dims)
            !  enddo
              do k = 1, n_doc
                 write(nchar,'(i3.3)') k
                 call define_rest_field(ncid,'bgc_DOC'//trim(nchar)    ,dims)
              enddo
              do k = 1, n_dic
                 write(nchar,'(i3.3)') k
                 call define_rest_field(ncid,'bgc_DIC'//trim(nchar)    ,dims)
              enddo
            endif
            if (tr_bgc_chl) then
              do k = 1, n_algae
                 write(nchar,'(i3.3)') k
                 call define_rest_field(ncid,'bgc_chl'//trim(nchar)    ,dims)
              enddo
            endif
            call define_rest_field(ncid,'bgc_Nit'  ,dims)
            if (tr_bgc_Am) &
            call define_rest_field(ncid,'bgc_Am'   ,dims)
            if (tr_bgc_Sil) &
            call define_rest_field(ncid,'bgc_Sil'  ,dims)
            if (tr_bgc_hum) &
            call define_rest_field(ncid,'bgc_hum'  ,dims)
            if (tr_bgc_DMS) then
              call define_rest_field(ncid,'bgc_DMSPp',dims)
              call define_rest_field(ncid,'bgc_DMSPd',dims)
              call define_rest_field(ncid,'bgc_DMS'  ,dims)
            endif
            if (tr_bgc_PON) &
            call define_rest_field(ncid,'bgc_PON'  ,dims)
            if (tr_bgc_DON) then
              do k = 1, n_don
                 write(nchar,'(i3.3)') k
                 call define_rest_field(ncid,'bgc_DON'//trim(nchar)    ,dims)
              enddo
            endif
            if (tr_bgc_Fe ) then
              do k = 1, n_fed
                 write(nchar,'(i3.3)') k
                 call define_rest_field(ncid,'bgc_Fed'//trim(nchar)    ,dims)
              enddo
              do k = 1, n_fep
                 write(nchar,'(i3.3)') k
                 call define_rest_field(ncid,'bgc_Fep'//trim(nchar)    ,dims)
              enddo
            endif
         endif   !skl_bgc

         !-----------------------------------------------------------------
         ! 4D restart fields, written as layers of 3D
         !-----------------------------------------------------------------

         do k=1,nilyr
            write(nchar,'(i3.3)') k
            call define_rest_field(ncid,'sice'//trim(nchar),dims)
            call define_rest_field(ncid,'qice'//trim(nchar),dims)
         enddo

         do k=1,nslyr
            write(nchar,'(i3.3)') k
            call define_rest_field(ncid,'qsno'//trim(nchar),dims)
         enddo

         if (tr_snow) then
            do k=1,nslyr
               write(nchar,'(i3.3)') k
               call define_rest_field(ncid,'smice'//trim(nchar),dims)
               call define_rest_field(ncid,'smliq'//trim(nchar),dims)
               call define_rest_field(ncid, 'rhos'//trim(nchar),dims)
               call define_rest_field(ncid, 'rsnw'//trim(nchar),dims)
            enddo
         endif

         if (tr_fsd) then
            do k=1,nfsd
               write(nchar,'(i3.3)') k
               call define_rest_field(ncid,'fsd'//trim(nchar),dims)
            enddo
         endif

         if (tr_iso) then
            do k=1,n_iso
               write(nchar,'(i3.3)') k
               call define_rest_field(ncid,'isosno'//trim(nchar),dims)
               call define_rest_field(ncid,'isoice'//trim(nchar),dims)
            enddo
         endif

         if (tr_aero) then
            do k=1,n_aero
               write(nchar,'(i3.3)') k
               call define_rest_field(ncid,'aerosnossl'//trim(nchar),dims)
               call define_rest_field(ncid,'aerosnoint'//trim(nchar),dims)
               call define_rest_field(ncid,'aeroicessl'//trim(nchar),dims)
               call define_rest_field(ncid,'aeroiceint'//trim(nchar),dims)
            enddo
         endif

         if (z_tracers) then
            if (tr_zaero) then
               do n = 1, n_zaero
                  write(ncharb,'(i3.3)') n
                  do k = 1, nblyr+3
                     write(nchar,'(i3.3)') k
                     call define_rest_field(ncid,'zaero'//trim(ncharb)//trim(nchar),dims)
                 enddo !k
              enddo  !n
            endif   !tr_zaero
            if (tr_bgc_Nit) then
               do k = 1, nblyr+3
                  write(nchar,'(i3.3)') k
                  call define_rest_field(ncid,'bgc_Nit'//trim(nchar),dims)
               enddo
            endif
            if (tr_bgc_N) then
               do n = 1, n_algae
                  write(ncharb,'(i3.3)') n
                  do k = 1, nblyr+3
                     write(nchar,'(i3.3)') k
                     call define_rest_field(ncid,'bgc_N'//trim(ncharb)//trim(nchar),dims)
                  enddo
               enddo
            endif
            if (tr_bgc_C) then
            !   do n = 1, n_algae
            !      write(ncharb,'(i3.3)') n
            !      do k = 1, nblyr+3
            !         write(nchar,'(i3.3)') k
            !         call define_rest_field(ncid,'bgc_C'//trim(ncharb)//trim(nchar),dims)
            !      enddo
            !   enddo
               do n = 1, n_doc
                  write(ncharb,'(i3.3)') n
                  do k = 1, nblyr+3
                     write(nchar,'(i3.3)') k
                     call define_rest_field(ncid,'bgc_DOC'//trim(ncharb)//trim(nchar),dims)
                  enddo
               enddo
               do n = 1, n_dic
                  write(ncharb,'(i3.3)') n
                  do k = 1, nblyr+3
                     write(nchar,'(i3.3)') k
                     call define_rest_field(ncid,'bgc_DIC'//trim(ncharb)//trim(nchar),dims)
                  enddo
               enddo
            endif
            if (tr_bgc_chl) then
               do n = 1, n_algae
                  write(ncharb,'(i3.3)') n
                  do k = 1, nblyr+3
                     write(nchar,'(i3.3)') k
                     call define_rest_field(ncid,'bgc_chl'//trim(ncharb)//trim(nchar),dims)
                  enddo
               enddo
            endif
            if (tr_bgc_Am) then
               do k = 1, nblyr+3
                  write(nchar,'(i3.3)') k
                  call define_rest_field(ncid,'bgc_Am'//trim(nchar),dims)
               enddo
            endif
            if (tr_bgc_Sil) then
               do k = 1, nblyr+3
                  write(nchar,'(i3.3)') k
                  call define_rest_field(ncid,'bgc_Sil'//trim(nchar),dims)
               enddo
            endif
            if (tr_bgc_hum) then
               do k = 1, nblyr+3
                  write(nchar,'(i3.3)') k
                  call define_rest_field(ncid,'bgc_hum'//trim(nchar),dims)
               enddo
            endif
            if (tr_bgc_DMS) then
               do k = 1, nblyr+3
                  write(nchar,'(i3.3)') k
                  call define_rest_field(ncid,'bgc_DMSPp'//trim(nchar),dims)
                  call define_rest_field(ncid,'bgc_DMSPd'//trim(nchar),dims)
                  call define_rest_field(ncid,'bgc_DMS'//trim(nchar),dims)
               enddo
            endif
            if (tr_bgc_PON) then
               do k = 1, nblyr+3
                  write(nchar,'(i3.3)') k
                  call define_rest_field(ncid,'bgc_PON'//trim(nchar),dims)
               enddo
            endif
            if (tr_bgc_DON) then
               do n = 1, n_don
                  write(ncharb,'(i3.3)') n
                  do k = 1, nblyr+3
                     write(nchar,'(i3.3)') k
                     call define_rest_field(ncid,'bgc_DON'//trim(ncharb)//trim(nchar),dims)
                  enddo
               enddo
            endif
            if (tr_bgc_Fe ) then
               do n = 1, n_fed
                  write(ncharb,'(i3.3)') n
                  do k = 1, nblyr+3
                     write(nchar,'(i3.3)') k
                     call define_rest_field(ncid,'bgc_Fed'//trim(ncharb)//trim(nchar),dims)
                  enddo
               enddo
               do n = 1, n_fep
                  write(ncharb,'(i3.3)') n
                  do k = 1, nblyr+3
                     write(nchar,'(i3.3)') k
                     call define_rest_field(ncid,'bgc_Fep'//trim(ncharb)//trim(nchar),dims)
                  enddo
               enddo
            endif
            do k = 1, nbtrcr
               write(nchar,'(i3.3)') k
               call define_rest_field(ncid,'zbgc_frac'//trim(nchar),dims)
            enddo
         endif   !z_tracers

         deallocate(dims)
         status = nf90_enddef(ncid)
         call ice_check_nc(status, subname//' ERROR: enddef', file=__FILE__, line=__LINE__)

         write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
      endif ! master_task

#else
      call abort_ice(subname//' ERROR: USE_NETCDF cpp not defined for '//trim(filename_spec), &
          file=__FILE__, line=__LINE__)
#endif

      end subroutine init_restart_write

!=======================================================================

! Reads a single restart field
! author David A Bailey, NCAR

      subroutine read_restart_field(nu,nrec,work,atype,vname,ndim3, &
                                    diag, field_loc, field_type)

      use ice_blocks, only: nx_block, ny_block
      use ice_domain_size, only: max_blocks, ncat
      use ice_read_write, only: ice_read_nc

      integer (kind=int_kind), intent(in) :: &
         nu            , & ! unit number (not used for netcdf)
         ndim3         , & ! third dimension
         nrec              ! record number (0 for sequential access)

      real (kind=dbl_kind), dimension(nx_block,ny_block,ndim3,max_blocks), intent(inout) :: &
         work              ! input array (real, 8-byte)

      character (len=4), intent(in) :: &
         atype             ! format for output array
                           ! (real/integer, 4-byte/8-byte)

      logical (kind=log_kind), intent(in) :: &
         diag              ! if true, write diagnostic output

      character (len=*), intent(in)  :: vname

      integer (kind=int_kind), optional, intent(in) :: &
         field_loc, &      ! location of field on staggered grid
         field_type        ! type of field (scalar, vector, angle)

      ! local variables

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
         work2              ! input array (real, 8-byte)

      character(len=*), parameter :: subname = '(read_restart_field)'

      work (:,:,:,:) = c0
      work2(:,:,:)   = c0
#ifdef USE_NETCDF
      if (present(field_loc)) then
         if (ndim3 == ncat) then
            call ice_read_nc(ncid,1,vname,work,diag, &
               field_loc=field_loc,field_type=field_type,restart_ext=restart_ext)
         elseif (ndim3 == 1) then
            call ice_read_nc(ncid,1,vname,work2,diag, &
               field_loc=field_loc,field_type=field_type,restart_ext=restart_ext)
            work(:,:,1,:) = work2(:,:,:)
         else
            write(nu_diag,*) 'ndim3 not supported ',ndim3
         endif
      else
         if (ndim3 == ncat) then
            call ice_read_nc(ncid, 1, vname, work, diag, restart_ext=restart_ext)
         elseif (ndim3 == 1) then
            call ice_read_nc(ncid, 1, vname, work2, diag, restart_ext=restart_ext)
            work(:,:,1,:) = work2(:,:,:)
         else
            write(nu_diag,*) 'ndim3 not supported ',ndim3
         endif
      endif

#else
      call abort_ice(subname//' ERROR: USE_NETCDF cpp not defined', &
          file=__FILE__, line=__LINE__)
#endif

      end subroutine read_restart_field

!=======================================================================

! Writes a single restart field.
! author David A Bailey, NCAR

      subroutine write_restart_field(nu,nrec,work,atype,vname,ndim3,diag)

      use ice_blocks, only: nx_block, ny_block
      use ice_domain_size, only: max_blocks, ncat
      use ice_read_write, only: ice_write_nc

      integer (kind=int_kind), intent(in) :: &
         nu            , & ! unit number
         ndim3         , & ! third dimension
         nrec              ! record number (0 for sequential access)

      real (kind=dbl_kind), dimension(nx_block,ny_block,ndim3,max_blocks), intent(in) :: &
         work              ! input array (real, 8-byte)

      character (len=4), intent(in) :: &
         atype             ! format for output array
                           ! (real/integer, 4-byte/8-byte)

      logical (kind=log_kind), intent(in) :: &
         diag              ! if true, write diagnostic output

      character (len=*), intent(in)  :: vname

      ! local variables

      integer (kind=int_kind) :: &
         varid         , & ! variable id
         status            ! status variable from netCDF routine

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
         work2             ! input array (real, 8-byte)

      character(len=*), parameter :: subname = '(write_restart_field)'

#ifdef USE_NETCDF
      varid = -99
      if (my_task == master_task) then
         ! ncid is only valid on master
         status = nf90_inq_varid(ncid,trim(vname),varid)
         call ice_check_nc(status, subname//' ERROR: inq varid '//trim(vname), file=__FILE__, line=__LINE__)
      endif
      if (ndim3 == ncat) then
         call ice_write_nc(ncid, 1, varid, work, diag, restart_ext, varname=trim(vname))
      elseif (ndim3 == 1) then
         work2(:,:,:) = work(:,:,1,:)
         call ice_write_nc(ncid, 1, varid, work2, diag, restart_ext, varname=trim(vname))
      else
         write(nu_diag,*) 'ndim3 not supported',ndim3
      endif

#else
      call abort_ice(subname//' ERROR: USE_NETCDF cpp not defined', &
          file=__FILE__, line=__LINE__)
#endif

      end subroutine write_restart_field

!=======================================================================

! Finalize the restart file.
! author David A Bailey, NCAR

      subroutine final_restart()

      use ice_calendar, only: istep1, myear, mmonth, mday, msec

      integer (kind=int_kind) :: status

      character(len=*), parameter :: subname = '(final_restart)'

#ifdef USE_NETCDF
      if (my_task == master_task) then
         ! ncid is only valid on master
         status = nf90_close(ncid)
         call ice_check_nc(status, subname//' ERROR: closing', file=__FILE__, line=__LINE__)
         write(nu_diag,'(a,i8,4x,i4.4,a,i2.2,a,i2.2,a,i5.5)') &
            'Restart read/written ',istep1,myear,'-',mmonth,'-',mday,'-',msec
      endif
#else
      call abort_ice(subname//' ERROR: USE_NETCDF cpp not defined', &
          file=__FILE__, line=__LINE__)
#endif

      end subroutine final_restart

!=======================================================================

! Defines a restart field
! author David A Bailey, NCAR

      subroutine define_rest_field(ncid, vname, dims)

      character (len=*)      , intent(in)  :: vname
      integer (kind=int_kind), intent(in)  :: dims(:)
      integer (kind=int_kind), intent(in)  :: ncid

      integer (kind=int_kind) :: varid

      integer (kind=int_kind) :: chunks(size(dims)), status, i

      character(len=*), parameter :: subname = '(define_rest_field)'

#ifdef USE_NETCDF

      status = nf90_def_var(ncid,trim(vname),nf90_double,dims,varid)
      call ice_check_nc(status, subname//' ERROR: def var '//trim(vname), file=__FILE__, line=__LINE__)

#ifdef NO_HDF5
      if (restart_format=='hdf5') then
         call abort_ice(subname//' ERROR: restart_format hdf5 not available ', &
            file=__FILE__, line=__LINE__)
      endif
#else
      if (restart_format=='hdf5' .and. size(dims)>1) then
         if (dims(1)==dimid_ni .and. dims(2)==dimid_nj) then
            chunks(1)=restart_chunksize(1)
            chunks(2)=restart_chunksize(2)
            do i = 3, size(dims)
               chunks(i) = 0
            enddo
            status = nf90_def_var_chunking(ncid, varid, NF90_CHUNKED, chunksizes=chunks)
            call ice_check_nc(status, subname//' ERROR: chunking var '//trim(vname), file=__FILE__, line=__LINE__)
         endif
      endif

      if (restart_format=='hdf5' .and. restart_deflate/=0) then
         status=nf90_def_var_deflate(ncid, varid, shuffle=0, deflate=1, deflate_level=restart_deflate)
         call ice_check_nc(status, subname//' ERROR deflating var '//trim(vname), file=__FILE__, line=__LINE__)
      endif
#endif

#else
      call abort_ice(subname//' ERROR: USE_NETCDF cpp not defined', &
          file=__FILE__, line=__LINE__)
#endif

      end subroutine define_rest_field

!=======================================================================

! Inquire field existance
! author T. Craig

      logical function query_field(nu,vname)

      integer (kind=int_kind), intent(in) :: nu     ! unit number
      character (len=*)      , intent(in) :: vname  ! variable name

      ! local variables

      integer (kind=int_kind) :: status, varid
      character(len=*), parameter :: subname = '(query_field)'

      query_field = .false.
#ifdef USE_NETCDF
      if (my_task == master_task) then
         status = nf90_inq_varid(ncid,trim(vname),varid)
         if (status == nf90_noerr) query_field = .true.
      endif
      call broadcast_scalar(query_field,master_task)
#else
      call abort_ice(subname//' ERROR: USE_NETCDF cpp not defined', &
          file=__FILE__, line=__LINE__)
#endif

      end function query_field

!=======================================================================

      end module ice_restart

!=======================================================================
