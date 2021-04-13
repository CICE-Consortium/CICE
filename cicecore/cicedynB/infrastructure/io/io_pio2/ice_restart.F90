!=======================================================================
!
! Read and write ice model restart files using pio interfaces.
! authors David A Bailey, NCAR

      module ice_restart

      use ice_broadcast
      use ice_exit, only: abort_ice
      use ice_fileunits, only: nu_diag, nu_restart, nu_rst_pointer
      use ice_kinds_mod
      use ice_restart_shared, only: &
          restart, restart_ext, restart_dir, restart_file, pointer_file, &
          runid, runtype, use_restart_time, restart_format, lcdf64, lenstr, &
          restart_coszen
      use ice_pio
      use pio
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_tracer_flags, icepack_query_tracer_indices
      use icepack_intfc, only: icepack_query_parameters, icepack_query_tracer_sizes

      implicit none
      private
      public :: init_restart_write, init_restart_read, &
                read_restart_field, write_restart_field, final_restart

      type(file_desc_t)     :: File
      type(var_desc_t)      :: vardesc

      type(io_desc_t)       :: iodesc2d
      type(io_desc_t)       :: iodesc3d_ncat

!=======================================================================

      contains

!=======================================================================

! Sets up restart file for reading.
! author David A Bailey, NCAR

      subroutine init_restart_read(ice_ic)

      use ice_calendar, only: istep0, istep1, myear, mmonth, &
                              mday, msec, npt
      use ice_communicate, only: my_task, master_task
      use ice_domain_size, only: ncat
      use ice_read_write, only: ice_open

      character(len=char_len_long), intent(in), optional :: ice_ic

      ! local variables

      character(len=char_len_long) :: &
         filename, filename0

      integer (kind=int_kind) :: status, status1

      integer (kind=int_kind) :: iotype

      character(len=*), parameter :: subname = '(init_restart_read)'

      if (present(ice_ic)) then 
         filename = trim(ice_ic)
      else
         if (my_task == master_task) then
            open(nu_rst_pointer,file=pointer_file)
            read(nu_rst_pointer,'(a)') filename0
            filename = trim(filename0)
            close(nu_rst_pointer)
            write(nu_diag,*) 'Read ',pointer_file(1:lenstr(pointer_file))
         endif
         call broadcast_scalar(filename, master_task)
      endif

      if (my_task == master_task) then
         write(nu_diag,*) 'Using restart dump=', trim(filename)
      end if

!     if (restart_format(1:3) == 'pio') then
         iotype = PIO_IOTYPE_NETCDF
         if (restart_format == 'pio_pnetcdf') iotype = PIO_IOTYPE_PNETCDF
         File%fh=-1
         call ice_pio_init(mode='read', filename=trim(filename), File=File, iotype=iotype)
      
         call ice_pio_initdecomp(iodesc=iodesc2d)
         call ice_pio_initdecomp(ndim3=ncat  , iodesc=iodesc3d_ncat,remap=.true.)

         if (use_restart_time) then
            status1 = PIO_noerr
            status = pio_get_att(File, pio_global, 'istep1', istep0)
!            status = pio_get_att(File, pio_global, 'time', time)
!            status = pio_get_att(File, pio_global, 'time_forc', time_forc)
            call pio_seterrorhandling(File, PIO_BCAST_ERROR)
            status = pio_get_att(File, pio_global, 'myear', myear)
            if (status /= PIO_noerr) status = pio_get_att(File, pio_global, 'nyr', myear)
            if (status /= PIO_noerr) status1 = status
            status = pio_get_att(File, pio_global, 'mmonth', mmonth)
            if (status /= PIO_noerr) status = pio_get_att(File, pio_global, 'month', mmonth)
            if (status /= PIO_noerr) status1 = status
            status = pio_get_att(File, pio_global, 'mday', mday)
            if (status /= PIO_noerr) status1 = status
            status = pio_get_att(File, pio_global, 'msec', msec)
            if (status /= PIO_noerr) status = pio_get_att(File, pio_global, 'sec', msec)
            if (status /= PIO_noerr) status1 = status
            if (status1 /= PIO_noerr) &
               call abort_ice(subname//"ERROR: reading restart time ")
            call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)
         endif ! use namelist values if use_restart_time = F
!     endif

      if (my_task == master_task) then
         write(nu_diag,*) 'Restart read at istep=',istep0,myear,mmonth,mday,msec
      endif

      call broadcast_scalar(istep0,master_task)
      call broadcast_scalar(myear,master_task)
      call broadcast_scalar(mmonth,master_task)
      call broadcast_scalar(mday,master_task)
      call broadcast_scalar(msec,master_task)
!      call broadcast_scalar(time,master_task)
!      call broadcast_scalar(time_forc,master_task)
      call broadcast_scalar(myear,master_task)
      
      istep1 = istep0

      ! if runid is bering then need to correct npt for istep0
      if (trim(runid) == 'bering') then
         npt = npt - istep0
      endif

      end subroutine init_restart_read

!=======================================================================

! Sets up restart file for writing.
! author David A Bailey, NCAR

      subroutine init_restart_write(filename_spec)

      use ice_calendar, only: msec, mmonth, mday, myear, istep1
      use ice_communicate, only: my_task, master_task
      use ice_domain_size, only: nx_global, ny_global, ncat, nilyr, nslyr, &
                                 n_iso, n_aero, nblyr, n_zaero, n_algae, n_doc,   &
                                 n_dic, n_don, n_fed, n_fep, nfsd
      use ice_dyn_shared, only: kdyn
      use ice_arrays_column, only: oceanmixed_ice

      logical (kind=log_kind) :: &
          solve_zsal, skl_bgc, z_tracers

      logical (kind=log_kind) :: &
          tr_iage, tr_FY, tr_lvl, tr_iso, tr_aero, tr_pond_cesm, &
          tr_pond_topo, tr_pond_lvl, tr_brine, &
          tr_bgc_N, tr_bgc_C, tr_bgc_Nit, &
          tr_bgc_Sil, tr_bgc_DMS, &
          tr_bgc_chl, tr_bgc_Am,  &
          tr_bgc_PON, tr_bgc_DON, &
          tr_zaero,   tr_bgc_Fe,  &
          tr_bgc_hum, tr_fsd

      integer (kind=int_kind) :: &
          nbtrcr

      character(len=char_len_long), intent(in), optional :: filename_spec

      ! local variables

      character(len=char_len_long) :: filename

      integer (kind=int_kind) :: dimid_ni, dimid_nj, dimid_ncat, &
                                 dimid_nilyr, dimid_nslyr, dimid_naero

      integer (kind=int_kind), allocatable :: dims(:)

      integer (kind=int_kind) :: iotype

      integer (kind=int_kind) :: &
        k,    n,    & ! loop index
        status        ! status variable from netCDF routine

      character (len=3) :: nchar, ncharb

      character(len=*), parameter :: subname = '(init_restart_write)'

      call icepack_query_tracer_sizes(nbtrcr_out=nbtrcr)
      call icepack_query_tracer_flags( &
          tr_iage_out=tr_iage, tr_FY_out=tr_FY, tr_lvl_out=tr_lvl, &
          tr_iso_out=tr_iso, tr_aero_out=tr_aero, tr_pond_cesm_out=tr_pond_cesm, &
          tr_pond_topo_out=tr_pond_topo, tr_pond_lvl_out=tr_pond_lvl, tr_brine_out=tr_brine, &
          tr_bgc_N_out=tr_bgc_N, tr_bgc_C_out=tr_bgc_C, tr_bgc_Nit_out=tr_bgc_Nit, &
          tr_bgc_Sil_out=tr_bgc_Sil, tr_bgc_DMS_out=tr_bgc_DMS, &
          tr_bgc_chl_out=tr_bgc_chl,  tr_bgc_Am_out=tr_bgc_Am, &
          tr_bgc_PON_out=tr_bgc_PON, tr_bgc_DON_out=tr_bgc_DON, &
          tr_zaero_out=tr_zaero,    tr_bgc_Fe_out=tr_bgc_Fe, &
          tr_bgc_hum_out=tr_bgc_hum, tr_fsd_out=tr_fsd)
      call icepack_query_parameters(solve_zsal_out=solve_zsal, skl_bgc_out=skl_bgc, &
          z_tracers_out=z_tracers)
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
        
      if (restart_format(1:3) /= 'bin') filename = trim(filename) // '.nc'

      ! write pointer (path/file)
      if (my_task == master_task) then
         open(nu_rst_pointer,file=pointer_file)
         write(nu_rst_pointer,'(a)') filename
         close(nu_rst_pointer)
      endif

!     if (restart_format(1:3) == 'pio') then
      
         iotype = PIO_IOTYPE_NETCDF
         if (restart_format == 'pio_pnetcdf') iotype = PIO_IOTYPE_PNETCDF
         File%fh=-1
         call ice_pio_init(mode='write',filename=trim(filename), File=File, &
              clobber=.true., cdf64=lcdf64, iotype=iotype)

         status = pio_put_att(File,pio_global,'istep1',istep1)
!         status = pio_put_att(File,pio_global,'time',time)
!         status = pio_put_att(File,pio_global,'time_forc',time_forc)
         status = pio_put_att(File,pio_global,'myear',myear)
         status = pio_put_att(File,pio_global,'mmonth',mmonth)
         status = pio_put_att(File,pio_global,'mday',mday)
         status = pio_put_att(File,pio_global,'msec',msec)

         status = pio_def_dim(File,'ni',nx_global,dimid_ni)
         status = pio_def_dim(File,'nj',ny_global,dimid_nj)
         status = pio_def_dim(File,'ncat',ncat,dimid_ncat)

      !-----------------------------------------------------------------
      ! 2D restart fields
      !-----------------------------------------------------------------

         allocate(dims(2))

         dims(1) = dimid_ni
         dims(2) = dimid_nj

         call define_rest_field(File,'uvel',dims)
         call define_rest_field(File,'vvel',dims)
         if (restart_coszen) call define_rest_field(File,'coszen',dims)
         call define_rest_field(File,'scale_factor',dims)
         call define_rest_field(File,'swvdr',dims)
         call define_rest_field(File,'swvdf',dims)
         call define_rest_field(File,'swidr',dims)
         call define_rest_field(File,'swidf',dims)

         call define_rest_field(File,'strocnxT',dims)
         call define_rest_field(File,'strocnyT',dims)

         call define_rest_field(File,'stressp_1',dims)
         call define_rest_field(File,'stressp_2',dims)
         call define_rest_field(File,'stressp_3',dims)
         call define_rest_field(File,'stressp_4',dims)

         call define_rest_field(File,'stressm_1',dims)
         call define_rest_field(File,'stressm_2',dims)
         call define_rest_field(File,'stressm_3',dims)
         call define_rest_field(File,'stressm_4',dims)

         call define_rest_field(File,'stress12_1',dims)
         call define_rest_field(File,'stress12_2',dims)
         call define_rest_field(File,'stress12_3',dims)
         call define_rest_field(File,'stress12_4',dims)

         call define_rest_field(File,'iceumask',dims)

         if (oceanmixed_ice) then
            call define_rest_field(File,'sst',dims)
            call define_rest_field(File,'frzmlt',dims)
         endif

         if (tr_FY) then
            call define_rest_field(File,'frz_onset',dims)
         end if

         if (kdyn == 2) then
            call define_rest_field(File,'a11_1',dims)
            call define_rest_field(File,'a11_2',dims)
            call define_rest_field(File,'a11_3',dims)
            call define_rest_field(File,'a11_4',dims)
            call define_rest_field(File,'a12_1',dims)
            call define_rest_field(File,'a12_2',dims)
            call define_rest_field(File,'a12_3',dims)
            call define_rest_field(File,'a12_4',dims)
         endif

         if (tr_pond_lvl) then
            call define_rest_field(File,'fsnow',dims)
         endif

         if (nbtrcr > 0) then
            if (tr_bgc_N) then
            do k=1,n_algae
               write(nchar,'(i3.3)') k
               call define_rest_field(File,'algalN'//trim(nchar),dims)
            enddo
            endif
            if (tr_bgc_C) then
            do k=1,n_doc
               write(nchar,'(i3.3)') k
               call define_rest_field(File,'doc'//trim(nchar),dims)
            enddo
            do k=1,n_dic
               write(nchar,'(i3.3)') k
               call define_rest_field(File,'dic'//trim(nchar),dims)
            enddo
            endif
            call define_rest_field(File,'nit'   ,dims)
            if (tr_bgc_Am) &
            call define_rest_field(File,'amm'   ,dims)
            if (tr_bgc_Sil) &
            call define_rest_field(File,'sil'   ,dims)
            if (tr_bgc_hum) &
            call define_rest_field(File,'hum'   ,dims)
            if (tr_bgc_DMS) then
              call define_rest_field(File,'dmsp'  ,dims)
              call define_rest_field(File,'dms'   ,dims)
            endif
            if (tr_bgc_DON) then
            do k=1,n_don
               write(nchar,'(i3.3)') k
               call define_rest_field(File,'don'//trim(nchar),dims)
            enddo
            endif
            if (tr_bgc_Fe ) then
            do k=1,n_fed
               write(nchar,'(i3.3)') k
               call define_rest_field(File,'fed'//trim(nchar),dims)
            enddo
            do k=1,n_fep
               write(nchar,'(i3.3)') k
               call define_rest_field(File,'fep'//trim(nchar),dims)
            enddo
            endif
            if (tr_zaero) then
            do k=1,n_zaero
               write(nchar,'(i3.3)') k
               call define_rest_field(File,'zaeros'//trim(nchar),dims)
            enddo
            endif
         endif  !nbtrcr

         if (solve_zsal) call define_rest_field(File,'sss',dims)

         deallocate(dims)

      !-----------------------------------------------------------------
      ! 3D restart fields (ncat)
      !-----------------------------------------------------------------

         allocate(dims(3))

         dims(1) = dimid_ni
         dims(2) = dimid_nj
         dims(3) = dimid_ncat

         call define_rest_field(File,'aicen',dims)
         call define_rest_field(File,'vicen',dims)
         call define_rest_field(File,'vsnon',dims)
         call define_rest_field(File,'Tsfcn',dims)

         if (tr_iage) then
            call define_rest_field(File,'iage',dims)
         end if

         if (tr_FY) then
            call define_rest_field(File,'FY',dims)
         end if

         if (tr_lvl) then
            call define_rest_field(File,'alvl',dims)
            call define_rest_field(File,'vlvl',dims)
         end if

         if (tr_pond_cesm) then
            call define_rest_field(File,'apnd',dims)
            call define_rest_field(File,'hpnd',dims)
         end if

         if (tr_pond_topo) then
            call define_rest_field(File,'apnd',dims)
            call define_rest_field(File,'hpnd',dims)
            call define_rest_field(File,'ipnd',dims)
         end if

         if (tr_pond_lvl) then
            call define_rest_field(File,'apnd',dims)
            call define_rest_field(File,'hpnd',dims)
            call define_rest_field(File,'ipnd',dims)
            call define_rest_field(File,'dhs',dims)
            call define_rest_field(File,'ffrac',dims)
         end if

         if (tr_brine) then
            call define_rest_field(File,'fbrn',dims)
            call define_rest_field(File,'first_ice',dims)
         endif

         if (skl_bgc) then
            do k = 1, n_algae
               write(nchar,'(i3.3)') k
               call define_rest_field(File,'bgc_N'//trim(nchar)    ,dims)
            enddo
            if (tr_bgc_C) then
            !  do k = 1, n_algae
            !     write(nchar,'(i3.3)') k
            !     call define_rest_field(File,'bgc_C'//trim(nchar)    ,dims)
            !  enddo
              do k = 1, n_doc
                 write(nchar,'(i3.3)') k
                 call define_rest_field(File,'bgc_DOC'//trim(nchar)    ,dims)
              enddo
              do k = 1, n_dic
                 write(nchar,'(i3.3)') k
                 call define_rest_field(File,'bgc_DIC'//trim(nchar)    ,dims)
              enddo
            endif
            if (tr_bgc_chl) then
              do k = 1, n_algae
                 write(nchar,'(i3.3)') k
                 call define_rest_field(File,'bgc_chl'//trim(nchar)    ,dims)
              enddo
            endif
            call define_rest_field(File,'bgc_Nit'  ,dims)
            if (tr_bgc_Am) &
            call define_rest_field(File,'bgc_Am'   ,dims)
            if (tr_bgc_Sil) &
            call define_rest_field(File,'bgc_Sil'  ,dims)
            if (tr_bgc_hum) &
            call define_rest_field(File,'bgc_hum'  ,dims)
            if (tr_bgc_DMS) then
              call define_rest_field(File,'bgc_DMSPp',dims)
              call define_rest_field(File,'bgc_DMSPd',dims)
              call define_rest_field(File,'bgc_DMS'  ,dims)
            endif
            if (tr_bgc_PON) &
            call define_rest_field(File,'bgc_PON'  ,dims)
            if (tr_bgc_DON) then
              do k = 1, n_don
                 write(nchar,'(i3.3)') k
                 call define_rest_field(File,'bgc_DON'//trim(nchar)    ,dims)
              enddo
            endif
            if (tr_bgc_Fe ) then
              do k = 1, n_fed
                 write(nchar,'(i3.3)') k
                 call define_rest_field(File,'bgc_Fed'//trim(nchar)    ,dims)
              enddo
              do k = 1, n_fep
                 write(nchar,'(i3.3)') k
                 call define_rest_field(File,'bgc_Fep'//trim(nchar)    ,dims)
              enddo
            endif
         endif   !skl_bgc
         if (solve_zsal) &
            call define_rest_field(File,'Rayleigh',dims)

      !-----------------------------------------------------------------
      ! 4D restart fields, written as layers of 3D
      !-----------------------------------------------------------------

         do k=1,nilyr
            write(nchar,'(i3.3)') k
            call define_rest_field(File,'sice'//trim(nchar),dims)
            call define_rest_field(File,'qice'//trim(nchar),dims)
         enddo

         do k=1,nslyr
            write(nchar,'(i3.3)') k
            call define_rest_field(File,'qsno'//trim(nchar),dims)
         enddo

         if (tr_fsd) then
            do k=1,nfsd
               write(nchar,'(i3.3)') k
               call define_rest_field(File,'fsd'//trim(nchar),dims)
            enddo
         endif

         if (tr_iso) then
            do k=1,n_iso
               write(nchar,'(i3.3)') k
               call define_rest_field(File,'isosno'//nchar, dims)
               call define_rest_field(File,'isoice'//nchar, dims)
            enddo
         endif

         if (tr_aero) then
            do k=1,n_aero
               write(nchar,'(i3.3)') k
               call define_rest_field(File,'aerosnossl'//nchar, dims)
               call define_rest_field(File,'aerosnoint'//nchar, dims)
               call define_rest_field(File,'aeroicessl'//nchar, dims)
               call define_rest_field(File,'aeroiceint'//nchar, dims)
            enddo
         endif

         if (solve_zsal) then
         do k = 1, nblyr
            write(nchar,'(i3.3)') k
            call define_rest_field(File,'zSalinity'//trim(nchar),dims)
         enddo
         endif
         if (z_tracers) then
            if (tr_zaero) then
             do n = 1, n_zaero
              write(ncharb,'(i3.3)') n
              do k = 1, nblyr+3
               write(nchar,'(i3.3)') k
               call define_rest_field(File,'zaero'//trim(ncharb)//trim(nchar),dims)
              enddo !k
             enddo  !n
            endif   !tr_zaero
            if (tr_bgc_Nit) then
              do k = 1, nblyr+3
               write(nchar,'(i3.3)') k
               call define_rest_field(File,'bgc_Nit'//trim(nchar),dims)
              enddo
            endif
            if (tr_bgc_N) then
             do n = 1, n_algae
              write(ncharb,'(i3.3)') n
              do k = 1, nblyr+3
               write(nchar,'(i3.3)') k
               call define_rest_field(File,'bgc_N'//trim(ncharb)//trim(nchar),dims)
              enddo
             enddo
            endif
            if (tr_bgc_C) then
            ! do n = 1, n_algae
            !  write(ncharb,'(i3.3)') n
            !  do k = 1, nblyr+3
            !     write(nchar,'(i3.3)') k
            !     call
            !     define_rest_field(File,'bgc_C'//trim(ncharb)//trim(nchar),dims)
            !  enddo
            ! enddo
             do n = 1, n_doc
              write(ncharb,'(i3.3)') n
              do k = 1, nblyr+3
                 write(nchar,'(i3.3)') k
                 call define_rest_field(File,'bgc_DOC'//trim(ncharb)//trim(nchar),dims)
              enddo
             enddo
             do n = 1, n_dic
              write(ncharb,'(i3.3)') n
              do k = 1, nblyr+3
                 write(nchar,'(i3.3)') k
                 call define_rest_field(File,'bgc_DIC'//trim(ncharb)//trim(nchar),dims)
              enddo
             enddo
            endif
            if (tr_bgc_chl) then
             do n = 1, n_algae
              write(ncharb,'(i3.3)') n
              do k = 1, nblyr+3
                 write(nchar,'(i3.3)') k
                 call define_rest_field(File,'bgc_chl'//trim(ncharb)//trim(nchar),dims)
              enddo
             enddo
            endif
            if (tr_bgc_Am) then
              do k = 1, nblyr+3
                 write(nchar,'(i3.3)') k
                 call define_rest_field(File,'bgc_Am'//trim(nchar),dims)
              enddo
            endif
            if (tr_bgc_Sil) then
              do k = 1, nblyr+3
                 write(nchar,'(i3.3)') k
                 call define_rest_field(File,'bgc_Sil'//trim(nchar),dims)
              enddo
            endif
            if (tr_bgc_hum) then
              do k = 1, nblyr+3
                 write(nchar,'(i3.3)') k
                 call define_rest_field(File,'bgc_hum'//trim(nchar),dims)
              enddo
            endif
            if (tr_bgc_DMS) then
              do k = 1, nblyr+3
                 write(nchar,'(i3.3)') k
                 call define_rest_field(File,'bgc_DMSPp'//trim(nchar),dims)
                 call define_rest_field(File,'bgc_DMSPd'//trim(nchar),dims)
                 call define_rest_field(File,'bgc_DMS'//trim(nchar),dims)
              enddo
            endif
            if (tr_bgc_PON) then
              do k = 1, nblyr+3
                 write(nchar,'(i3.3)') k
                 call define_rest_field(File,'bgc_PON'//trim(nchar),dims)
              enddo
            endif
            if (tr_bgc_DON) then
             do n = 1, n_don
              write(ncharb,'(i3.3)') n
              do k = 1, nblyr+3
                 write(nchar,'(i3.3)') k
                 call define_rest_field(File,'bgc_DON'//trim(ncharb)//trim(nchar),dims)
              enddo
             enddo
            endif
            if (tr_bgc_Fe ) then
             do n = 1, n_fed
              write(ncharb,'(i3.3)') n
              do k = 1, nblyr+3
                 write(nchar,'(i3.3)') k
                 call define_rest_field(File,'bgc_Fed'//trim(ncharb)//trim(nchar),dims)
              enddo
             enddo
             do n = 1, n_fep
              write(ncharb,'(i3.3)') n
              do k = 1, nblyr+3
                 write(nchar,'(i3.3)') k
                 call define_rest_field(File,'bgc_Fep'//trim(ncharb)//trim(nchar),dims)
              enddo
             enddo
            endif
            do k = 1, nbtrcr
               write(nchar,'(i3.3)') k
               call define_rest_field(File,'zbgc_frac'//trim(nchar),dims)
            enddo
         endif   !z_tracers

         deallocate(dims)
         status = pio_enddef(File)

         call ice_pio_initdecomp(iodesc=iodesc2d)
         call ice_pio_initdecomp(ndim3=ncat  , iodesc=iodesc3d_ncat, remap=.true.)

!     endif  ! restart_format

      if (my_task == master_task) then
         write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
      endif

      end subroutine init_restart_write

!=======================================================================

! Reads a single restart field
! author David A Bailey, NCAR

      subroutine read_restart_field(nu,nrec,work,atype,vname,ndim3,diag, &
                                    field_loc, field_type)

      use ice_blocks, only: nx_block, ny_block
      use ice_communicate, only: my_task, master_task
      use ice_constants, only: c0, field_loc_center
      use ice_boundary, only: ice_HaloUpdate
      use ice_domain, only: halo_info, distrb_info, nblocks
      use ice_domain_size, only: max_blocks, ncat
      use ice_global_reductions, only: global_minval, global_maxval, global_sum

      integer (kind=int_kind), intent(in) :: &
           nu            , & ! unit number
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

      integer (kind=int_kind) :: &
        j,     &      ! dimension counter
        n,     &      ! number of dimensions for variable
        ndims, &      ! number of variable dimensions
        status        ! status variable from netCDF routine

      real (kind=dbl_kind) :: amin,amax,asum

      character(len=*), parameter :: subname = '(read_restart_field)'

!     if (restart_format(1:3) == "pio") then
         if (my_task == master_task) &
            write(nu_diag,*)'Parallel restart file read: ',vname

         call pio_seterrorhandling(File, PIO_BCAST_ERROR)

         status = pio_inq_varid(File,trim(vname),vardesc)

         if (status /= PIO_noerr) then
            call abort_ice(subname//"ERROR: CICE restart? Missing variable: "//trim(vname))
         endif

         status = pio_inq_varndims(File, vardesc, ndims)

         call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)

!         if (ndim3 == ncat .and. ncat>1) then
         if (ndim3 == ncat .and. ndims == 3) then
            call pio_read_darray(File, vardesc, iodesc3d_ncat, work, status)
            if (present(field_loc)) then
               do n=1,ndim3
                  call ice_HaloUpdate (work(:,:,n,:), halo_info, &
                                       field_loc, field_type)
               enddo
            endif
!         elseif (ndim3 == 1) then
         elseif (ndim3 == 1 .and. ndims == 2) then
            call pio_read_darray(File, vardesc, iodesc2d, work, status)
            if (present(field_loc)) then
               call ice_HaloUpdate (work(:,:,1,:), halo_info, &
                                    field_loc, field_type)
            endif
         else
            write(nu_diag,*) "ndim3 not supported ",ndim3
         endif

         if (diag) then
            if (ndim3 > 1) then
               do n=1,ndim3
                  amin = global_minval(work(:,:,n,:),distrb_info)
                  amax = global_maxval(work(:,:,n,:),distrb_info)
                  asum = global_sum(work(:,:,n,:), distrb_info, field_loc_center)
                  if (my_task == master_task) then
                     write(nu_diag,*) ' min and max =', amin, amax
                     write(nu_diag,*) ' sum =',asum
                  endif
               enddo
            else
               amin = global_minval(work(:,:,1,:),distrb_info)
               amax = global_maxval(work(:,:,1,:),distrb_info)
               asum = global_sum(work(:,:,1,:), distrb_info, field_loc_center)
               if (my_task == master_task) then
                  write(nu_diag,*) ' min and max =', amin, amax
                  write(nu_diag,*) ' sum =',asum
                  write(nu_diag,*) ''
               endif
            endif
         
         endif
!     else
!        call abort_ice(subname//"ERROR: Invalid restart_format: "//trim(restart_format))
!     endif  ! restart_format

      end subroutine read_restart_field
      
!=======================================================================

! Writes a single restart field.
! author David A Bailey, NCAR

      subroutine write_restart_field(nu,nrec,work,atype,vname,ndim3,diag)

      use ice_blocks, only: nx_block, ny_block
      use ice_communicate, only: my_task, master_task
      use ice_constants, only: c0, field_loc_center
      use ice_domain, only: distrb_info, nblocks
      use ice_domain_size, only: max_blocks, ncat
      use ice_global_reductions, only: global_minval, global_maxval, global_sum

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
        j,     &      ! dimension counter
        n,     &      ! dimension counter
        ndims, &  ! number of variable dimensions
        status        ! status variable from netCDF routine

      real (kind=dbl_kind) :: amin,amax,asum

      character(len=*), parameter :: subname = '(write_restart_field)'

!      if (restart_format(1:3) == "pio") then
         if (my_task == master_task) &
            write(nu_diag,*)'Parallel restart file write: ',vname

         status = pio_inq_varid(File,trim(vname),vardesc)
         
         status = pio_inq_varndims(File, vardesc, ndims)

         if (ndims==3) then 
            call pio_write_darray(File, vardesc, iodesc3d_ncat,work(:,:,:,1:nblocks), &
                 status, fillval=c0)
         elseif (ndims == 2) then
            call pio_write_darray(File, vardesc, iodesc2d, work(:,:,1,1:nblocks), &
                 status, fillval=c0)
         else
            write(nu_diag,*) "ndims not supported",ndims,ndim3
         endif

         if (diag) then
            if (ndim3 > 1) then
               do n=1,ndim3
                  amin = global_minval(work(:,:,n,:),distrb_info)
                  amax = global_maxval(work(:,:,n,:),distrb_info)
                  asum = global_sum(work(:,:,n,:), distrb_info, field_loc_center)
                  if (my_task == master_task) then
                     write(nu_diag,*) ' min and max =', amin, amax
                     write(nu_diag,*) ' sum =',asum
                  endif
               enddo
            else
               amin = global_minval(work(:,:,1,:),distrb_info)
               amax = global_maxval(work(:,:,1,:),distrb_info)
               asum = global_sum(work(:,:,1,:), distrb_info, field_loc_center)
               if (my_task == master_task) then
                  write(nu_diag,*) ' min and max =', amin, amax
                  write(nu_diag,*) ' sum =',asum
               endif
            endif
         endif
!     else
!        call abort_ice(subname//"ERROR: Invalid restart_format: "//trim(restart_format))
!     endif

      end subroutine write_restart_field

!=======================================================================

! Finalize the restart file.
! author David A Bailey, NCAR

      subroutine final_restart()

      use ice_calendar, only: istep1, idate, msec
      use ice_communicate, only: my_task, master_task

      character(len=*), parameter :: subname = '(final_restart)'

      call PIO_freeDecomp(File,iodesc2d)
      call PIO_freeDecomp(File,iodesc3d_ncat)
      call pio_closefile(File)

      if (my_task == master_task) &
         write(nu_diag,*) 'Restart read/written ',istep1,idate,msec

      end subroutine final_restart

!=======================================================================

! Defines a restart field
! author David A Bailey, NCAR

      subroutine define_rest_field(File, vname, dims)

      type(file_desc_t)      , intent(in)  :: File
      character (len=*)      , intent(in)  :: vname
      integer (kind=int_kind), intent(in)  :: dims(:)

      integer (kind=int_kind) :: &
        status        ! status variable from netCDF routine

      character(len=*), parameter :: subname = '(define_rest_field)'

      status = pio_def_var(File,trim(vname),pio_double,dims,vardesc)
        
      end subroutine define_rest_field

!=======================================================================

      end module ice_restart

!=======================================================================
