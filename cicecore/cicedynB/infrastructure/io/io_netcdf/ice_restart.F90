!=======================================================================

! Read and write ice model restart files using netCDF or binary
! interfaces.
! authors David A Bailey, NCAR

      module ice_restart

      use ice_broadcast
      use ice_kinds_mod
      use netcdf
      use ice_restart_shared, only: &
          restart_ext, restart_dir, restart_file, pointer_file, &
          runid, use_restart_time, lcdf64, lenstr
      use ice_fileunits, only: nu_diag, nu_rst_pointer
      use ice_exit, only: abort_ice
      use icepack_intfc, only: icepack_query_parameters
      use icepack_intfc, only: icepack_query_tracer_numbers
      use icepack_intfc, only: icepack_query_tracer_flags
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted

      implicit none
      private
      public :: init_restart_write, init_restart_read, &
                read_restart_field, write_restart_field, final_restart

      integer (kind=int_kind) :: ncid

!=======================================================================

      contains

!=======================================================================

! Sets up restart file for reading.
! author David A Bailey, NCAR

      subroutine init_restart_read(ice_ic)

      use ice_calendar, only: sec, month, mday, nyr, istep0, istep1, &
                              time, time_forc, npt
      use ice_communicate, only: my_task, master_task

      character(len=char_len_long), intent(in), optional :: ice_ic

      ! local variables

      character(len=char_len_long) :: &
         filename, filename0

      integer (kind=int_kind) :: status

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

         status = nf90_open(trim(filename), nf90_nowrite, ncid)
         if (status /= nf90_noerr) call abort_ice(subname// &
            'ERROR: reading restart ncfile '//trim(filename))
      
         if (use_restart_time) then
         status = nf90_get_att(ncid, nf90_global, 'istep1', istep0)
         status = nf90_get_att(ncid, nf90_global, 'time', time)
         status = nf90_get_att(ncid, nf90_global, 'time_forc', time_forc)
         status = nf90_get_att(ncid, nf90_global, 'nyr', nyr)
         if (status == nf90_noerr) then
            status = nf90_get_att(ncid, nf90_global, 'month', month)
            status = nf90_get_att(ncid, nf90_global, 'mday', mday)
            status = nf90_get_att(ncid, nf90_global, 'sec', sec)
         endif
         endif ! use namelist values if use_restart_time = F

         write(nu_diag,*) 'Restart read at istep=',istep0,time,time_forc
      endif

      call broadcast_scalar(istep0,master_task)
      call broadcast_scalar(time,master_task)
      call broadcast_scalar(time_forc,master_task)
      call broadcast_scalar(nyr,master_task)
      
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

      use ice_blocks, only: nghost
      use ice_calendar, only: sec, month, mday, nyr, istep1, &
                              time, time_forc, year_init
      use ice_communicate, only: my_task, master_task
      use ice_domain_size, only: nx_global, ny_global, ncat, nilyr, nslyr, &
                                 n_aero, nblyr, n_zaero, n_algae, n_doc,   &
                                 n_dic, n_don, n_fed, n_fep
      use ice_arrays_column, only: oceanmixed_ice
      use ice_dyn_shared, only: kdyn

      character(len=char_len_long), intent(in), optional :: filename_spec

      ! local variables

      logical (kind=log_kind) :: &
         solve_zsal, skl_bgc, z_tracers, &
         tr_iage, tr_FY, tr_lvl, tr_aero, tr_pond_cesm, &
         tr_pond_topo, tr_pond_lvl, tr_brine, &
         tr_bgc_N, tr_bgc_C, tr_bgc_Nit, &
         tr_bgc_Sil, tr_bgc_DMS, &
         tr_bgc_chl,  tr_bgc_Am, &
         tr_bgc_PON, tr_bgc_DON, &
         tr_zaero,    tr_bgc_Fe, &
         tr_bgc_hum

      integer (kind=int_kind) :: &
         k,  n,                & ! index
         nx, ny,               & ! global array size
         iyear, imonth, iday,  & ! year, month, day
         nbtrcr                  ! number of bgc tracers

      character(len=char_len_long) :: filename

      integer (kind=int_kind), allocatable :: dims(:)

      integer (kind=int_kind) :: &
        dimid_ni,   & ! netCDF identifiers
        dimid_nj,   & !
        dimid_ncat, & !
        iflag,      & ! netCDF creation flag
        status        ! status variable from netCDF routine

      character (len=3) :: nchar, ncharb

      character(len=*), parameter :: subname = '(init_restart_write)'

      call icepack_query_parameters( &
         solve_zsal_out=solve_zsal, skl_bgc_out=skl_bgc, z_tracers_out=z_tracers)
      call icepack_query_tracer_numbers( &
         nbtrcr_out=nbtrcr)
      call icepack_query_tracer_flags( &
         tr_iage_out=tr_iage, tr_FY_out=tr_FY, tr_lvl_out=tr_lvl, &
         tr_aero_out=tr_aero, tr_pond_cesm_out=tr_pond_cesm, &
         tr_pond_topo_out=tr_pond_topo, tr_pond_lvl_out=tr_pond_lvl, tr_brine_out=tr_brine, &
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
         iyear = nyr + year_init - 1
      
         write(filename,'(a,a,a,i4.4,a,i2.2,a,i2.2,a,i5.5)') &
              restart_dir(1:lenstr(restart_dir)), &
              restart_file(1:lenstr(restart_file)),'.', &
              iyear,'-',month,'-',mday,'-',sec
      end if

      ! write pointer (path/file)
      if (my_task == master_task) then
         filename = trim(filename) // '.nc'
         open(nu_rst_pointer,file=pointer_file)
         write(nu_rst_pointer,'(a)') filename
         close(nu_rst_pointer)

         iflag = 0
         if (lcdf64) iflag = nf90_64bit_offset
         status = nf90_create(trim(filename), iflag, ncid)
         if (status /= nf90_noerr) call abort_ice(subname// &
            'ERROR: creating restart ncfile '//trim(filename))

         status = nf90_put_att(ncid,nf90_global,'istep1',istep1)
         status = nf90_put_att(ncid,nf90_global,'time',time)
         status = nf90_put_att(ncid,nf90_global,'time_forc',time_forc)
         status = nf90_put_att(ncid,nf90_global,'nyr',nyr)
         status = nf90_put_att(ncid,nf90_global,'month',month)
         status = nf90_put_att(ncid,nf90_global,'mday',mday)
         status = nf90_put_att(ncid,nf90_global,'sec',sec)

         nx = nx_global
         ny = ny_global
         if (restart_ext) then
            nx = nx_global + 2*nghost
            ny = ny_global + 2*nghost
         endif
         status = nf90_def_dim(ncid,'ni',nx,dimid_ni)
         status = nf90_def_dim(ncid,'nj',ny,dimid_nj)

         status = nf90_def_dim(ncid,'ncat',ncat,dimid_ncat)

      !-----------------------------------------------------------------
      ! 2D restart fields
      !-----------------------------------------------------------------

         allocate(dims(2))

         dims(1) = dimid_ni
         dims(2) = dimid_nj

         call define_rest_field(ncid,'uvel',dims)
         call define_rest_field(ncid,'vvel',dims)

#ifdef CESMCOUPLED
         call define_rest_field(ncid,'coszen',dims)
#endif
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

         if (tr_pond_lvl) then
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

         if (solve_zsal) call define_rest_field(ncid,'sss',dims)

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

         if (tr_pond_cesm) then
            call define_rest_field(ncid,'apnd',dims)
            call define_rest_field(ncid,'hpnd',dims)
         end if

         if (tr_pond_topo) then
            call define_rest_field(ncid,'apnd',dims)
            call define_rest_field(ncid,'hpnd',dims)
            call define_rest_field(ncid,'ipnd',dims)
         end if

         if (tr_pond_lvl) then
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
         if (solve_zsal) &
            call define_rest_field(ncid,'Rayleigh',dims)

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

         if (tr_aero) then
            do k=1,n_aero
               write(nchar,'(i3.3)') k
               call define_rest_field(ncid,'aerosnossl'//trim(nchar),dims)
               call define_rest_field(ncid,'aerosnoint'//trim(nchar),dims)
               call define_rest_field(ncid,'aeroicessl'//trim(nchar),dims)
               call define_rest_field(ncid,'aeroiceint'//trim(nchar),dims)
            enddo
         endif

         if (solve_zsal) then
         do k = 1, nblyr
            write(nchar,'(i3.3)') k
            call define_rest_field(ncid,'zSalinity'//trim(nchar),dims)
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
            ! do n = 1, n_algae
            !  write(ncharb,'(i3.3)') n
            !  do k = 1, nblyr+3
            !     write(nchar,'(i3.3)') k
            !     call define_rest_field(ncid,'bgc_C'//trim(ncharb)//trim(nchar),dims)
            !  enddo
            ! enddo
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

         write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
      endif ! master_task

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

      integer (kind=int_kind) :: &
        n,     &      ! number of dimensions for variable
        varid, &      ! variable id
        status        ! status variable from netCDF routine

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
           work2              ! input array (real, 8-byte)

      character(len=*), parameter :: subname = '(read_restart_field)'

         if (present(field_loc)) then
            if (ndim3 == ncat) then
               if (restart_ext) then
                  call ice_read_nc(ncid,1,vname,work,diag, &
                     field_loc=field_loc,field_type=field_type,restart_ext=restart_ext)
               else
                  call ice_read_nc(ncid,1,vname,work,diag,field_loc,field_type)
               endif
            elseif (ndim3 == 1) then
               if (restart_ext) then
                  call ice_read_nc(ncid,1,vname,work2,diag, &
                     field_loc=field_loc,field_type=field_type,restart_ext=restart_ext)
               else
                  call ice_read_nc(ncid,1,vname,work2,diag,field_loc,field_type)
               endif
               work(:,:,1,:) = work2(:,:,:)
            else
               write(nu_diag,*) 'ndim3 not supported ',ndim3
            endif
         else
            if (ndim3 == ncat) then
               if (restart_ext) then
                  call ice_read_nc(ncid, 1, vname, work, diag, restart_ext=restart_ext)
               else
                  call ice_read_nc(ncid, 1, vname, work, diag)
               endif
            elseif (ndim3 == 1) then
               if (restart_ext) then
                  call ice_read_nc(ncid, 1, vname, work2, diag, restart_ext=restart_ext)
               else
                  call ice_read_nc(ncid, 1, vname, work2, diag)
               endif
               work(:,:,1,:) = work2(:,:,:)
            else
               write(nu_diag,*) 'ndim3 not supported ',ndim3
            endif
         endif

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
        varid, &      ! variable id
        status        ! status variable from netCDF routine

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
           work2              ! input array (real, 8-byte)

      character(len=*), parameter :: subname = '(write_restart_field)'

         status = nf90_inq_varid(ncid,trim(vname),varid)
         if (ndim3 == ncat) then 
            if (restart_ext) then
               call ice_write_nc(ncid, 1, varid, work, diag, restart_ext, varname=trim(vname))
            else
               call ice_write_nc(ncid, 1, varid, work, diag, varname=trim(vname))
            endif
         elseif (ndim3 == 1) then
            work2(:,:,:) = work(:,:,1,:)
            if (restart_ext) then
               call ice_write_nc(ncid, 1, varid, work2, diag, restart_ext, varname=trim(vname))
            else
               call ice_write_nc(ncid, 1, varid, work2, diag, varname=trim(vname))
            endif
         else
            write(nu_diag,*) 'ndim3 not supported',ndim3
         endif

      end subroutine write_restart_field

!=======================================================================

! Finalize the restart file.
! author David A Bailey, NCAR

      subroutine final_restart()

      use ice_calendar, only: istep1, time, time_forc
      use ice_communicate, only: my_task, master_task

      integer (kind=int_kind) :: status

      character(len=*), parameter :: subname = '(final_restart)'

      status = nf90_close(ncid)

      if (my_task == master_task) &
         write(nu_diag,*) 'Restart read/written ',istep1,time,time_forc

      end subroutine final_restart

!=======================================================================

! Defines a restart field
! author David A Bailey, NCAR

      subroutine define_rest_field(ncid, vname, dims)

      character (len=*)      , intent(in)  :: vname
      integer (kind=int_kind), intent(in)  :: dims(:)
      integer (kind=int_kind), intent(in)  :: ncid

      integer (kind=int_kind) :: varid

      integer (kind=int_kind) :: &
        status        ! status variable from netCDF routine

      character(len=*), parameter :: subname = '(define_rest_field)'

      status = nf90_def_var(ncid,trim(vname),nf90_double,dims,varid)
        
      end subroutine define_rest_field

!=======================================================================

      end module ice_restart

!=======================================================================
