!=======================================================================

! Read and write ice model restart files using netCDF or binary
! interfaces.

! authors David A Bailey, NCAR

      module ice_restart

      use ice_broadcast
      use ice_kinds_mod
      use ice_restart_shared, only: &
          restart, restart_ext, restart_dir, restart_file, pointer_file, &
          runid, runtype, use_restart_time, restart_format, lenstr
      use ice_fileunits, only: nu_diag, nu_rst_pointer
      use ice_fileunits, only: nu_dump, nu_dump_eap, nu_dump_FY, nu_dump_age
      use ice_fileunits, only: nu_dump_lvl, nu_dump_pond, nu_dump_hbrine
      use ice_fileunits, only: nu_dump_bgc, nu_dump_aero, nu_dump_age
      use ice_fileunits, only: nu_restart, nu_restart_eap, nu_restart_FY, nu_restart_age 
      use ice_fileunits, only: nu_restart_lvl, nu_restart_pond, nu_restart_hbrine
      use ice_fileunits, only: nu_restart_bgc, nu_restart_aero, nu_restart_age
      use ice_exit, only: abort_ice
      use icepack_intfc, only: icepack_query_parameters
      use icepack_intfc, only: icepack_query_tracer_numbers
      use icepack_intfc, only: icepack_query_tracer_flags
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted

      implicit none
      private
      public :: init_restart_write, init_restart_read, &
                read_restart_field, write_restart_field, final_restart

!=======================================================================

      contains

!=======================================================================

! Sets up restart file for reading.
! author David A Bailey, NCAR

      subroutine init_restart_read(ice_ic)

      use ice_calendar, only: istep0, istep1, time, time_forc, npt, nyr
      use ice_communicate, only: my_task, master_task
      use ice_dyn_shared, only: kdyn
      use ice_read_write, only: ice_open, ice_open_ext

      character(len=char_len_long), intent(in), optional :: ice_ic

      ! local variables

      logical (kind=log_kind) :: &
         solve_zsal, &
         tr_iage, tr_FY, tr_lvl, tr_aero, tr_pond_cesm, &
         tr_pond_topo, tr_pond_lvl, tr_brine

      character(len=char_len_long) :: &
         filename, filename0

      integer (kind=int_kind) :: &
         n, &                    ! loop indices
         nbtrcr, &               ! number of bgc tracers
         iignore                 ! dummy variable

      real (kind=real_kind) :: &
         rignore                 ! dummy variable

      character(len=char_len_long) :: &
         string1, string2

      character(len=*), parameter :: subname = '(init_restart_read)'

      call icepack_query_parameters( &
         solve_zsal_out=solve_zsal)
      call icepack_query_tracer_numbers( &
         nbtrcr_out=nbtrcr)
      call icepack_query_tracer_flags( &
         tr_iage_out=tr_iage, tr_FY_out=tr_FY, tr_lvl_out=tr_lvl, &
         tr_aero_out=tr_aero, tr_pond_cesm_out=tr_pond_cesm, &
         tr_pond_topo_out=tr_pond_topo, tr_pond_lvl_out=tr_pond_lvl, tr_brine_out=tr_brine)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

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
         if (restart_ext) then
            call ice_open_ext(nu_restart,trim(filename),0)
         else
            call ice_open(nu_restart,trim(filename),0)
         endif
         if (use_restart_time) then
            read (nu_restart) istep0,time,time_forc,nyr
         else
            read (nu_restart) iignore,rignore,rignore ! use namelist values
         endif
         write(nu_diag,*) 'Restart read at istep=',istep0,time,time_forc
      endif

      call broadcast_scalar(istep0,master_task)
      call broadcast_scalar(time,master_task)
      call broadcast_scalar(time_forc,master_task)
      call broadcast_scalar(nyr,master_task)
      
      istep1 = istep0

      ! Supplemental restart files

      if (kdyn == 2) then
         if (my_task == master_task) then
            n = index(filename0,trim(restart_file))
            if (n == 0) call abort_ice(subname//'ERROR: eap restart: filename discrepancy')
            string1 = trim(filename0(1:n-1))
            string2 = trim(filename0(n+lenstr(restart_file):lenstr(filename0)))
            write(filename,'(a,a,a,a)') &
               string1(1:lenstr(string1)), &
               restart_file(1:lenstr(restart_file)),'.eap', &
               string2(1:lenstr(string2))
            if (restart_ext) then
               call ice_open_ext(nu_restart_eap,filename,0)
            else
               call ice_open(nu_restart_eap,filename,0)
            endif
            read (nu_restart_eap) iignore,rignore,rignore
            write(nu_diag,*) 'Reading ',filename(1:lenstr(filename))
         endif
      endif

      if (tr_iage) then
         if (my_task == master_task) then
            n = index(filename0,trim(restart_file))
            if (n == 0) call abort_ice(subname//'ERROR: iage restart: filename discrepancy')
            string1 = trim(filename0(1:n-1))
            string2 = trim(filename0(n+lenstr(restart_file):lenstr(filename0)))
            write(filename,'(a,a,a,a)') &
               string1(1:lenstr(string1)), &
               restart_file(1:lenstr(restart_file)),'.iage', &
               string2(1:lenstr(string2))
            if (restart_ext) then
               call ice_open_ext(nu_restart_age,filename,0)
            else
               call ice_open(nu_restart_age,filename,0)
            endif
            read (nu_restart_age) iignore,rignore,rignore
            write(nu_diag,*) 'Reading ',filename(1:lenstr(filename))
         endif
      endif

      if (tr_FY) then
         if (my_task == master_task) then
            n = index(filename0,trim(restart_file))
            if (n == 0) call abort_ice(subname//'ERROR: FY restart: filename discrepancy')
            string1 = trim(filename0(1:n-1))
            string2 = trim(filename0(n+lenstr(restart_file):lenstr(filename0)))
            write(filename,'(a,a,a,a)') &
               string1(1:lenstr(string1)), &
               restart_file(1:lenstr(restart_file)),'.FY', &
               string2(1:lenstr(string2))
            if (restart_ext) then
               call ice_open_ext(nu_restart_FY,filename,0)
            else
               call ice_open(nu_restart_FY,filename,0)
            endif
            read (nu_restart_FY) iignore,rignore,rignore
            write(nu_diag,*) 'Reading ',filename(1:lenstr(filename))
         endif
      endif

      if (tr_lvl) then
         if (my_task == master_task) then
            n = index(filename0,trim(restart_file))
            if (n == 0) call abort_ice(subname//'ERROR: lvl restart: filename discrepancy')
            string1 = trim(filename0(1:n-1))
            string2 = trim(filename0(n+lenstr(restart_file):lenstr(filename0)))
            write(filename,'(a,a,a,a)') &
               string1(1:lenstr(string1)), &
               restart_file(1:lenstr(restart_file)),'.lvl', &
               string2(1:lenstr(string2))
            if (restart_ext) then
               call ice_open_ext(nu_restart_lvl,filename,0)
            else
               call ice_open(nu_restart_lvl,filename,0)
            endif
            read (nu_restart_lvl) iignore,rignore,rignore
            write(nu_diag,*) 'Reading ',filename(1:lenstr(filename))
         endif
      endif

      if (tr_pond_cesm) then
         if (my_task == master_task) then
            n = index(filename0,trim(restart_file))
            if (n == 0) call abort_ice(subname//'ERROR: pond_cesm restart: filename discrepancy')
            string1 = trim(filename0(1:n-1))
            string2 = trim(filename0(n+lenstr(restart_file):lenstr(filename0)))
            write(filename,'(a,a,a,a)') &
               string1(1:lenstr(string1)), &
               restart_file(1:lenstr(restart_file)),'.pond_cesm', &
               string2(1:lenstr(string2))
            if (restart_ext) then
               call ice_open_ext(nu_restart_pond,filename,0)
            else
               call ice_open(nu_restart_pond,filename,0)
            endif
            read (nu_restart_pond) iignore,rignore,rignore
            write(nu_diag,*) 'Reading ',filename(1:lenstr(filename))
         endif
      endif

      if (tr_pond_lvl) then
         if (my_task == master_task) then
            n = index(filename0,trim(restart_file))
            if (n == 0) call abort_ice(subname//'ERROR:pond_lvl restart: filename discrepancy')
            string1 = trim(filename0(1:n-1))
            string2 = trim(filename0(n+lenstr(restart_file):lenstr(filename0)))
            write(filename,'(a,a,a,a)') &
               string1(1:lenstr(string1)), &
               restart_file(1:lenstr(restart_file)),'.pond_lvl', &
               string2(1:lenstr(string2))
            if (restart_ext) then
               call ice_open_ext(nu_restart_pond,filename,0)
            else
               call ice_open(nu_restart_pond,filename,0)
            endif
            read (nu_restart_pond) iignore,rignore,rignore
            write(nu_diag,*) 'Reading ',filename(1:lenstr(filename))
         endif
      endif

      if (tr_pond_topo) then
         if (my_task == master_task) then
            n = index(filename0,trim(restart_file))
            if (n == 0) call abort_ice(subname//'ERROR: pond_topo restart: filename discrepancy')
            string1 = trim(filename0(1:n-1))
            string2 = trim(filename0(n+lenstr(restart_file):lenstr(filename0)))
            write(filename,'(a,a,a,a)') &
               string1(1:lenstr(string1)), &
               restart_file(1:lenstr(restart_file)),'.pond_topo', &
               string2(1:lenstr(string2))
            if (restart_ext) then
               call ice_open_ext(nu_restart_pond,filename,0)
            else
               call ice_open(nu_restart_pond,filename,0)
            endif
            read (nu_restart_pond) iignore,rignore,rignore
            write(nu_diag,*) 'Reading ',filename(1:lenstr(filename))
         endif
      endif

      if (tr_brine) then
         if (my_task == master_task) then
            n = index(filename0,trim(restart_file))
            if (n == 0) call abort_ice(subname//'ERROR: brine restart: filename discrepancy')
            string1 = trim(filename0(1:n-1))
            string2 = trim(filename0(n+lenstr(restart_file):lenstr(filename0)))
            write(filename,'(a,a,a,a)') &
               string1(1:lenstr(string1)), &
               restart_file(1:lenstr(restart_file)),'.brine', &
               string2(1:lenstr(string2))
            if (restart_ext) then
               call ice_open_ext(nu_restart_hbrine,filename,0)
            else
               call ice_open(nu_restart_hbrine,filename,0)
            endif
            read (nu_restart_hbrine) iignore,rignore,rignore
            write(nu_diag,*) 'Reading ',filename(1:lenstr(filename))
         endif
      endif

      if (solve_zsal .or. nbtrcr > 0) then
         if (my_task == master_task) then
            n = index(filename0,trim(restart_file))
            if (n == 0) call abort_ice(subname//'ERROR: bgc restart: filename discrepancy')
            string1 = trim(filename0(1:n-1))
            string2 = trim(filename0(n+lenstr(restart_file):lenstr(filename0)))
            write(filename,'(a,a,a,a)') &
               string1(1:lenstr(string1)), &
               restart_file(1:lenstr(restart_file)),'.bgc', &
               string2(1:lenstr(string2))
            if (restart_ext) then
               call ice_open_ext(nu_restart_bgc,filename,0)
            else
               call ice_open(nu_restart_bgc,filename,0)
            endif
            read (nu_restart_bgc) iignore,rignore,rignore
            write(nu_diag,*) 'Reading ',filename(1:lenstr(filename))
         endif
      endif

      if (tr_aero) then
         if (my_task == master_task) then
            n = index(filename0,trim(restart_file))
            if (n == 0) call abort_ice(subname//'ERROR: aero restart: filename discrepancy')
            string1 = trim(filename0(1:n-1))
            string2 = trim(filename0(n+lenstr(restart_file):lenstr(filename0)))
            write(filename,'(a,a,a,a)') &
               string1(1:lenstr(string1)), &
               restart_file(1:lenstr(restart_file)),'.aero', &
               string2(1:lenstr(string2))
            if (restart_ext) then
               call ice_open_ext(nu_restart_aero,filename,0)
            else
               call ice_open(nu_restart_aero,filename,0)
            endif
            read (nu_restart_aero) iignore,rignore,rignore
            write(nu_diag,*) 'Reading ',filename(1:lenstr(filename))
         endif
      endif

      ! if runid is bering then need to correct npt for istep0
      if (trim(runid) == 'bering') then
         npt = npt - istep0
      endif

      end subroutine init_restart_read

!=======================================================================

! Sets up restart file for writing.
! author David A Bailey, NCAR

      subroutine init_restart_write(filename_spec)

      use ice_calendar, only: sec, month, mday, nyr, istep1, &
                              time, time_forc, year_init
      use ice_communicate, only: my_task, master_task
      use ice_dyn_shared, only: kdyn
      use ice_read_write, only: ice_open, ice_open_ext

      character(len=char_len_long), intent(in), optional :: filename_spec

      ! local variables

      logical (kind=log_kind) :: &
         solve_zsal, &
         tr_iage, tr_FY, tr_lvl, tr_aero, tr_pond_cesm, &
         tr_pond_topo, tr_pond_lvl, tr_brine

      integer (kind=int_kind) :: &
         nbtrcr, &               ! number of bgc tracers
         iyear, imonth, iday     ! year, month, day

      character(len=char_len_long) :: filename

      character(len=*), parameter :: subname = '(init_restart_write)'

      call icepack_query_parameters( &
         solve_zsal_out=solve_zsal)
      call icepack_query_tracer_numbers( &
         nbtrcr_out=nbtrcr)
      call icepack_query_tracer_flags( &
         tr_iage_out=tr_iage, tr_FY_out=tr_FY, tr_lvl_out=tr_lvl, &
         tr_aero_out=tr_aero, tr_pond_cesm_out=tr_pond_cesm, &
         tr_pond_topo_out=tr_pond_topo, tr_pond_lvl_out=tr_pond_lvl, tr_brine_out=tr_brine)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      ! construct path/file
      if (present(filename_spec)) then
         filename = trim(filename_spec)
      else
         iyear = nyr + year_init - 1
         imonth = month
         iday = mday
      
         write(filename,'(a,a,a,i4.4,a,i2.2,a,i2.2,a,i5.5)') &
              restart_dir(1:lenstr(restart_dir)), &
              restart_file(1:lenstr(restart_file)),'.', &
              iyear,'-',month,'-',mday,'-',sec
      end if
        
      ! write pointer (path/file)
      if (my_task == master_task) then
         open(nu_rst_pointer,file=pointer_file)
         write(nu_rst_pointer,'(a)') filename
         close(nu_rst_pointer)
         if (restart_ext) then
            call ice_open_ext(nu_dump,filename,0)
         else
            call ice_open(nu_dump,filename,0)
         endif
         write(nu_dump) istep1,time,time_forc,nyr
         write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
      endif

      ! begin writing restart data

      if (kdyn == 2) then

         write(filename,'(a,a,a,i4.4,a,i2.2,a,i2.2,a,i5.5)') &
              restart_dir(1:lenstr(restart_dir)), &
              restart_file(1:lenstr(restart_file)),'.eap.', &
              iyear,'-',month,'-',mday,'-',sec

         if (restart_ext) then
            call ice_open_ext(nu_dump_eap,filename,0)
         else
            call ice_open(nu_dump_eap,filename,0)
         endif

         if (my_task == master_task) then
           write(nu_dump_eap) istep1,time,time_forc
           write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
         endif

      endif

      if (tr_FY) then

         write(filename,'(a,a,a,i4.4,a,i2.2,a,i2.2,a,i5.5)') &
              restart_dir(1:lenstr(restart_dir)), &
              restart_file(1:lenstr(restart_file)),'.FY.', &
              iyear,'-',month,'-',mday,'-',sec

         if (restart_ext) then
            call ice_open_ext(nu_dump_FY,filename,0)
         else
            call ice_open(nu_dump_FY,filename,0)
         endif

         if (my_task == master_task) then
           write(nu_dump_FY) istep1,time,time_forc
           write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
         endif

      endif

      if (tr_iage) then

         write(filename,'(a,a,a,i4.4,a,i2.2,a,i2.2,a,i5.5)') &
              restart_dir(1:lenstr(restart_dir)), &
              restart_file(1:lenstr(restart_file)),'.iage.', &
              iyear,'-',month,'-',mday,'-',sec

         if (restart_ext) then
            call ice_open_ext(nu_dump_age,filename,0)
         else
            call ice_open(nu_dump_age,filename,0)
         endif

         if (my_task == master_task) then
           write(nu_dump_age) istep1,time,time_forc
           write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
         endif

      endif

      if (tr_lvl) then

         write(filename,'(a,a,a,i4.4,a,i2.2,a,i2.2,a,i5.5)') &
              restart_dir(1:lenstr(restart_dir)), &
              restart_file(1:lenstr(restart_file)),'.lvl.', &
              iyear,'-',month,'-',mday,'-',sec

         if (restart_ext) then
            call ice_open_ext(nu_dump_lvl,filename,0)
         else
            call ice_open(nu_dump_lvl,filename,0)
         endif

         if (my_task == master_task) then
           write(nu_dump_lvl) istep1,time,time_forc
           write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
         endif

      endif

      if (tr_pond_cesm) then

         write(filename,'(a,a,a,i4.4,a,i2.2,a,i2.2,a,i5.5)') &
              restart_dir(1:lenstr(restart_dir)), &
              restart_file(1:lenstr(restart_file)),'.pond_cesm.', &
              iyear,'-',month,'-',mday,'-',sec

         if (restart_ext) then
            call ice_open_ext(nu_dump_pond,filename,0)
         else
            call ice_open(nu_dump_pond,filename,0)
         endif

         if (my_task == master_task) then
           write(nu_dump_pond) istep1,time,time_forc
           write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
         endif

      endif

      if (tr_pond_lvl) then

         write(filename,'(a,a,a,i4.4,a,i2.2,a,i2.2,a,i5.5)') &
              restart_dir(1:lenstr(restart_dir)), &
              restart_file(1:lenstr(restart_file)),'.pond_lvl.', &
              iyear,'-',month,'-',mday,'-',sec

         if (restart_ext) then
            call ice_open_ext(nu_dump_pond,filename,0)
         else
            call ice_open(nu_dump_pond,filename,0)
         endif

         if (my_task == master_task) then
           write(nu_dump_pond) istep1,time,time_forc
           write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
         endif

      endif

      if (tr_pond_topo) then

         write(filename,'(a,a,a,i4.4,a,i2.2,a,i2.2,a,i5.5)') &
              restart_dir(1:lenstr(restart_dir)), &
              restart_file(1:lenstr(restart_file)),'.pond_topo.', &
              iyear,'-',month,'-',mday,'-',sec

         if (restart_ext) then
            call ice_open_ext(nu_dump_pond,filename,0)
         else
            call ice_open(nu_dump_pond,filename,0)
         endif

         if (my_task == master_task) then
           write(nu_dump_pond) istep1,time,time_forc
           write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
         endif

      endif

      if (tr_brine) then

         write(filename,'(a,a,a,i4.4,a,i2.2,a,i2.2,a,i5.5)') &
              restart_dir(1:lenstr(restart_dir)), &
              restart_file(1:lenstr(restart_file)),'.brine.', &
              iyear,'-',month,'-',mday,'-',sec

         if (restart_ext) then
            call ice_open_ext(nu_dump_hbrine,filename,0)
         else
            call ice_open(nu_dump_hbrine,filename,0)
         endif

         if (my_task == master_task) then
           write(nu_dump_hbrine) istep1,time,time_forc
           write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
         endif

      endif

      if (solve_zsal .or. nbtrcr > 0) then

         write(filename,'(a,a,a,i4.4,a,i2.2,a,i2.2,a,i5.5)') &
              restart_dir(1:lenstr(restart_dir)), &
              restart_file(1:lenstr(restart_file)),'.bgc.', &
              iyear,'-',month,'-',mday,'-',sec

         if (restart_ext) then
            call ice_open_ext(nu_dump_bgc,filename,0)
         else
            call ice_open(nu_dump_bgc,filename,0)
         endif

         if (my_task == master_task) then
           write(nu_dump_bgc) istep1,time,time_forc
           write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
         endif
      endif

      if (tr_aero) then

         write(filename,'(a,a,a,i4.4,a,i2.2,a,i2.2,a,i5.5)') &
              restart_dir(1:lenstr(restart_dir)), &
              restart_file(1:lenstr(restart_file)),'.aero.', &
              iyear,'-',month,'-',mday,'-',sec

         if (restart_ext) then
            call ice_open_ext(nu_dump_aero,filename,0)
         else
            call ice_open(nu_dump_aero,filename,0)
         endif

         if (my_task == master_task) then
           write(nu_dump_aero) istep1,time,time_forc
           write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
         endif

      endif

      end subroutine init_restart_write

!=======================================================================

! Reads a single restart field
! author David A Bailey, NCAR

      subroutine read_restart_field(nu,nrec,work,atype,vname,ndim3, &
                                    diag, field_loc, field_type)

      use ice_blocks, only: nx_block, ny_block
      use ice_domain_size, only: max_blocks
      use ice_read_write, only: ice_read, ice_read_ext

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

      character (len=*), intent(in) :: vname

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

         write(nu_diag,*) 'vname ',trim(vname)
         if (present(field_loc)) then
            do n=1,ndim3
               if (restart_ext) then
                  call ice_read_ext(nu,nrec,work2,atype,diag)
               else
                  call ice_read(nu,nrec,work2,atype,diag,field_loc,field_type)
               endif
               work(:,:,n,:) = work2(:,:,:)
            enddo
         else
            do n=1,ndim3
               if (restart_ext) then
                  call ice_read_ext(nu,nrec,work2,atype,diag)
               else
                  call ice_read(nu,nrec,work2,atype,diag)
               endif
               work(:,:,n,:) = work2(:,:,:)
            enddo
         endif

      end subroutine read_restart_field
      
!=======================================================================

! Writes a single restart field.
! author David A Bailey, NCAR

      subroutine write_restart_field(nu,nrec,work,atype,vname,ndim3,diag)

      use ice_blocks, only: nx_block, ny_block
      use ice_domain_size, only: max_blocks
      use ice_read_write, only: ice_write, ice_write_ext

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
        n,     &      ! dimension counter
        varid, &      ! variable id
        status        ! status variable from netCDF routine

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
           work2              ! input array (real, 8-byte)

      character(len=*), parameter :: subname = '(write_restart_field)'

         do n=1,ndim3
            work2(:,:,:) = work(:,:,n,:)
            if (restart_ext) then
               call ice_write_ext(nu,nrec,work2,atype,diag)
            else
               call ice_write(nu,nrec,work2,atype,diag)
            endif
         enddo

      end subroutine write_restart_field

!=======================================================================

! Finalize the restart file.
! author David A Bailey, NCAR

      subroutine final_restart()

      use ice_calendar, only: istep1, time, time_forc
      use ice_communicate, only: my_task, master_task

      logical (kind=log_kind) :: &
         solve_zsal, &
         tr_iage, tr_FY, tr_lvl, tr_aero, tr_pond_cesm, &
         tr_pond_topo, tr_pond_lvl, tr_brine

      integer (kind=int_kind) :: &
         nbtrcr, &               ! number of bgc tracers
         status

      character(len=*), parameter :: subname = '(final_restart)'

      call icepack_query_parameters( &
         solve_zsal_out=solve_zsal)
      call icepack_query_tracer_numbers( &
         nbtrcr_out=nbtrcr)
      call icepack_query_tracer_flags( &
         tr_iage_out=tr_iage, tr_FY_out=tr_FY, tr_lvl_out=tr_lvl, &
         tr_aero_out=tr_aero, tr_pond_cesm_out=tr_pond_cesm, &
         tr_pond_topo_out=tr_pond_topo, tr_pond_lvl_out=tr_pond_lvl, tr_brine_out=tr_brine)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (my_task == master_task) then
         close(nu_dump)

         if (tr_aero)      close(nu_dump_aero)
         if (tr_iage)      close(nu_dump_age)
         if (tr_FY)        close(nu_dump_FY)
         if (tr_lvl)       close(nu_dump_lvl)
         if (tr_pond_cesm) close(nu_dump_pond)
         if (tr_pond_lvl)  close(nu_dump_pond)
         if (tr_pond_topo) close(nu_dump_pond)
         if (tr_brine)     close(nu_dump_hbrine)
         if (solve_zsal .or. nbtrcr > 0) &
                           close(nu_dump_bgc)

         write(nu_diag,*) 'Restart read/written ',istep1,time,time_forc
      endif

      end subroutine final_restart

!=======================================================================

      end module ice_restart

!=======================================================================
