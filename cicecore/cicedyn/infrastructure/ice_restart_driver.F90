!=======================================================================

! Read and write ice model restart files
!
! authors Elizabeth C. Hunke, LANL
!         William H. Lipscomb LANL
!         David Bailey, NCAR
!
! 2004-05: Block structure added by William Lipscomb
!          Restart module separated from history module
! 2006 ECH: Accepted some CESM code into mainstream CICE
!           Converted to free source form (F90)
! 2008 ECH: Rearranged order in which internal stresses are written and read
! 2010 ECH: Changed eice, esno to qice, qsno
! 2012 ECH: Added routines for reading/writing extended grid
! 2013 DAB: Added generic interfaces for writing restart fields.

      module ice_restart_driver

      use ice_kinds_mod
      use ice_arrays_column, only: oceanmixed_ice
      use ice_communicate, only: my_task, master_task
      use ice_constants, only: c0, c1, p5, &
          field_loc_center, field_loc_NEcorner, &
          field_loc_Eface, field_loc_Nface, &
          field_type_scalar, field_type_vector
      use ice_restart_shared, only: restart_dir, pointer_file, &
          runid, use_restart_time, lenstr, restart_coszen, restart_mod
      use ice_restart
      use ice_exit, only: abort_ice
      use ice_fileunits, only: nu_diag, nu_rst_pointer, nu_restart, nu_dump
      use ice_grid, only: tmask, opmask, grid_ice, grid_average_X2Y
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_aggregate
      use icepack_intfc, only: icepack_query_tracer_indices, icepack_query_tracer_sizes
      use icepack_intfc, only: icepack_query_parameters
      use icepack_intfc, only: icepack_query_tracer_flags

      implicit none
      private
      public :: dumpfile, restartfile, restartfile_v4

!=======================================================================

      contains

!=======================================================================

!=======================================================================
!---! these subroutines write/read Fortran unformatted data files ..
!=======================================================================

! Dumps all values needed for a restart
! author Elizabeth C. Hunke, LANL

      subroutine dumpfile(filename_spec)

      use ice_blocks, only: nx_block, ny_block
      use ice_domain, only: nblocks
      use ice_domain_size, only: nilyr, nslyr, ncat, max_blocks
      use ice_dyn_shared, only: iceUmask, iceEmask, iceNmask, kdyn
      use ice_flux, only: scale_factor, swvdr, swvdf, swidr, swidf, &
          strocnxT_iavg, strocnyT_iavg, sst, frzmlt, &
          stressp_1, stressp_2, stressp_3, stressp_4, &
          stressm_1, stressm_2, stressm_3, stressm_4, &
          stress12_1, stress12_2, stress12_3, stress12_4, &
          stresspT, stressmT, stress12T, &
          stresspU, stressmU, stress12U
      use ice_flux, only: coszen
      use ice_state, only: aicen, vicen, vsnon, trcrn, uvel, vvel, &
                           uvelE, vvelE, uvelN, vvelN

      character(len=*), intent(in), optional :: filename_spec

      ! local variables

      integer (kind=int_kind) :: &
          i, j, k, iblk, &     ! counting indices
          nt_Tsfc, nt_sice, nt_qice, nt_qsno

      logical (kind=log_kind) :: diag

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1

      character (len=3) :: nchar

      character(len=*), parameter :: subname = '(dumpfile)'

      call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_sice_out=nt_sice, &
           nt_qice_out=nt_qice, nt_qsno_out=nt_qsno)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (present(filename_spec)) then
         call init_restart_write(filename_spec=filename_spec)
      else
         call init_restart_write
      endif

      diag = .true.

      !-----------------------------------------------------------------
      ! Zero out tracers over land
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            if (.not.tmask(i,j,iblk) .and. .not.opmask(i,j,iblk)) trcrn(i,j,:,:,iblk) = c0
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !-----------------------------------------------------------------
      ! state variables
      ! Tsfc is the only tracer written to binary files.  All other
      ! tracers are written to their own dump/restart binary files.
      !-----------------------------------------------------------------

      call write_restart_field(nu_dump,0,aicen(:,:,:,:),'ruf8','aicen',ncat,diag)
      call write_restart_field(nu_dump,0,vicen(:,:,:,:),'ruf8','vicen',ncat,diag)
      call write_restart_field(nu_dump,0,vsnon(:,:,:,:),'ruf8','vsnon',ncat,diag)
      call write_restart_field(nu_dump,0,trcrn(:,:,nt_Tsfc,:,:),'ruf8','Tsfcn',ncat,diag)

      do k=1,nilyr
         write(nchar,'(i3.3)') k
         call write_restart_field(nu_dump,0,trcrn(:,:,nt_sice+k-1,:,:),'ruf8', &
                                 'sice'//trim(nchar),ncat,diag)
      enddo

      do k=1,nilyr
         write(nchar,'(i3.3)') k
         call write_restart_field(nu_dump,0,trcrn(:,:,nt_qice+k-1,:,:),'ruf8', &
                                 'qice'//trim(nchar),ncat,diag)
      enddo

      do k=1,nslyr
         write(nchar,'(i3.3)') k
         call write_restart_field(nu_dump,0,trcrn(:,:,nt_qsno+k-1,:,:),'ruf8', &
                                 'qsno'//trim(nchar),ncat,diag)
      enddo

      !-----------------------------------------------------------------
      ! velocity
      !-----------------------------------------------------------------
      call write_restart_field(nu_dump,0,uvel,'ruf8','uvel',1,diag)
      call write_restart_field(nu_dump,0,vvel,'ruf8','vvel',1,diag)

      if (grid_ice == 'CD') then
         call write_restart_field(nu_dump,0,uvelE,'ruf8','uvelE',1,diag)
         call write_restart_field(nu_dump,0,vvelE,'ruf8','vvelE',1,diag)
         call write_restart_field(nu_dump,0,uvelN,'ruf8','uvelN',1,diag)
         call write_restart_field(nu_dump,0,vvelN,'ruf8','vvelN',1,diag)
      endif

      if (grid_ice == 'C') then
         call write_restart_field(nu_dump,0,uvelE,'ruf8','uvelE',1,diag)
         call write_restart_field(nu_dump,0,vvelN,'ruf8','vvelN',1,diag)
      endif

      !-----------------------------------------------------------------
      ! radiation fields
      !-----------------------------------------------------------------

      if (restart_coszen) call write_restart_field(nu_dump,0,coszen,'ruf8','coszen',1,diag)

      call write_restart_field(nu_dump,0,scale_factor,'ruf8','scale_factor',1,diag)

      call write_restart_field(nu_dump,0,swvdr,'ruf8','swvdr',1,diag)
      call write_restart_field(nu_dump,0,swvdf,'ruf8','swvdf',1,diag)
      call write_restart_field(nu_dump,0,swidr,'ruf8','swidr',1,diag)
      call write_restart_field(nu_dump,0,swidf,'ruf8','swidf',1,diag)

      !-----------------------------------------------------------------
      ! ocean stress (for bottom heat flux in thermo)
      !-----------------------------------------------------------------
      call write_restart_field(nu_dump,0,strocnxT_iavg,'ruf8','strocnxT',1,diag)
      call write_restart_field(nu_dump,0,strocnyT_iavg,'ruf8','strocnyT',1,diag)

      !-----------------------------------------------------------------
      ! internal stress
      !-----------------------------------------------------------------
      call write_restart_field(nu_dump,0,stressp_1,'ruf8','stressp_1',1,diag)
      call write_restart_field(nu_dump,0,stressp_3,'ruf8','stressp_3',1,diag)
      call write_restart_field(nu_dump,0,stressp_2,'ruf8','stressp_2',1,diag)
      call write_restart_field(nu_dump,0,stressp_4,'ruf8','stressp_4',1,diag)

      call write_restart_field(nu_dump,0,stressm_1,'ruf8','stressm_1',1,diag)
      call write_restart_field(nu_dump,0,stressm_3,'ruf8','stressm_3',1,diag)
      call write_restart_field(nu_dump,0,stressm_2,'ruf8','stressm_2',1,diag)
      call write_restart_field(nu_dump,0,stressm_4,'ruf8','stressm_4',1,diag)

      call write_restart_field(nu_dump,0,stress12_1,'ruf8','stress12_1',1,diag)
      call write_restart_field(nu_dump,0,stress12_3,'ruf8','stress12_3',1,diag)
      call write_restart_field(nu_dump,0,stress12_2,'ruf8','stress12_2',1,diag)
      call write_restart_field(nu_dump,0,stress12_4,'ruf8','stress12_4',1,diag)

      if (grid_ice == 'CD') then
         call write_restart_field(nu_dump,0,stresspT ,'ruf8','stresspT' ,1,diag)
         call write_restart_field(nu_dump,0,stressmT ,'ruf8','stressmT' ,1,diag)
         call write_restart_field(nu_dump,0,stress12T,'ruf8','stress12T',1,diag)
         call write_restart_field(nu_dump,0,stresspU ,'ruf8','stresspU' ,1,diag)
         call write_restart_field(nu_dump,0,stressmU ,'ruf8','stressmU' ,1,diag)
         call write_restart_field(nu_dump,0,stress12U,'ruf8','stress12U',1,diag)
      endif

      if (grid_ice == 'C') then
         call write_restart_field(nu_dump,0,stresspT ,'ruf8','stresspT' ,1,diag)
         call write_restart_field(nu_dump,0,stressmT ,'ruf8','stressmT' ,1,diag)
         call write_restart_field(nu_dump,0,stress12U,'ruf8','stress12U',1,diag)
      endif

      !-----------------------------------------------------------------
      ! ice mask for dynamics
      !-----------------------------------------------------------------
      if (kdyn > 0) then
         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               work1(i,j,iblk) = c0
               if (iceUmask(i,j,iblk)) work1(i,j,iblk) = c1
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO
         call write_restart_field(nu_dump,0,work1,'ruf8','iceumask',1,diag)

         if (grid_ice == 'CD' .or. grid_ice == 'C') then

            !$OMP PARALLEL DO PRIVATE(iblk,i,j)
            do iblk = 1, nblocks
               do j = 1, ny_block
               do i = 1, nx_block
                  work1(i,j,iblk) = c0
                  if (iceNmask(i,j,iblk)) work1(i,j,iblk) = c1
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO
            call write_restart_field(nu_dump,0,work1,'ruf8','icenmask',1,diag)

            !$OMP PARALLEL DO PRIVATE(iblk,i,j)
            do iblk = 1, nblocks
               do j = 1, ny_block
               do i = 1, nx_block
                  work1(i,j,iblk) = c0
                  if (iceEmask(i,j,iblk)) work1(i,j,iblk) = c1
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO
            call write_restart_field(nu_dump,0,work1,'ruf8','iceemask',1,diag)
         endif
      else
         work1(:,:,:) = c0
         call write_restart_field(nu_dump,0,work1,'ruf8','iceumask',1,diag)
         if (grid_ice == 'CD' .or. grid_ice == 'C') then
            call write_restart_field(nu_dump,0,work1,'ruf8','icenmask',1,diag)
            call write_restart_field(nu_dump,0,work1,'ruf8','iceemask',1,diag)
         endif
      endif

      ! for mixed layer model
      if (oceanmixed_ice) then
         call write_restart_field(nu_dump,0,sst,'ruf8','sst',1,diag)
         call write_restart_field(nu_dump,0,frzmlt,'ruf8','frzmlt',1,diag)
      endif

      end subroutine dumpfile

!=======================================================================

! Restarts from a dump
! author Elizabeth C. Hunke, LANL

      subroutine restartfile (ice_ic)

      use ice_boundary, only: ice_HaloUpdate_stress
      use ice_blocks, only: nghost, nx_block, ny_block
      use ice_calendar, only: istep0, npt, calendar
      use ice_domain, only: nblocks, halo_info, ns_boundary_type
      use ice_domain_size, only: nilyr, nslyr, ncat, &
          max_blocks
      use ice_dyn_shared, only: iceUmask, iceEmask, iceNmask,kdyn
      use ice_flux, only: scale_factor, swvdr, swvdf, swidr, swidf, &
          strocnxT_iavg, strocnyT_iavg, sst, frzmlt, &
          stressp_1, stressp_2, stressp_3, stressp_4, &
          stressm_1, stressm_2, stressm_3, stressm_4, &
          stress12_1, stress12_2, stress12_3, stress12_4, &
          stresspT, stressmT, stress12T, &
          stresspU, stressmU, stress12U
      use ice_flux, only: coszen, Tf
      use ice_state, only: trcr_depend, aice, vice, vsno, trcr, &
          aice0, aicen, vicen, vsnon, trcrn, aice_init, uvel, vvel, &
          uvelE, vvelE, uvelN, vvelN, &
          trcr_base, nt_strata, n_trcr_strata
      use icepack_itd, only: cleanup_itd  !for restart_mod
      use ice_arrays_column, only: first_ice, hin_max
      use ice_flux, only: fpond, fresh, fsalt, fhocn
      use ice_flux_bgc, only: faero_ocn, fiso_ocn, flux_bio
      use ice_calendar, only: dt

      character (*), optional :: ice_ic

      ! local variables

      integer (kind=int_kind) :: &
         i, j, k, iblk, &     ! counting indices
         ntrcr, &             ! number of tracers
         nt_Tsfc, nt_sice, nt_qice, nt_qsno

      logical (kind=log_kind) :: &
         diag, &
         tr_aero, tr_pond_topo

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1

      character (len=3) :: nchar

      character(len=*), parameter :: subname = '(restartfile)'

      call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__, line=__LINE__)

      call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_sice_out=nt_sice, &
           nt_qice_out=nt_qice, nt_qsno_out=nt_qsno)
      call icepack_query_tracer_flags(tr_aero_out=tr_aero, tr_pond_topo_out=tr_pond_topo)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      call init_restart_read(ice_ic)
      call calendar()

      diag = .true.

      !-----------------------------------------------------------------
      ! state variables
      ! Tsfc is the only tracer read in this file.  All other
      ! tracers are in their own dump/restart files.
      !-----------------------------------------------------------------
      if (my_task == master_task) &
           write(nu_diag,*) ' min/max area, vol ice, vol snow, Tsfc'

      call read_restart_field(nu_restart,0,aicen,'ruf8', &
              'aicen',ncat,diag,field_loc_center, field_type_scalar)
      call read_restart_field(nu_restart,0,vicen,'ruf8', &
              'vicen',ncat,diag,field_loc_center, field_type_scalar)
      call read_restart_field(nu_restart,0,vsnon,'ruf8', &
              'vsnon',ncat,diag,field_loc_center, field_type_scalar)
      call read_restart_field(nu_restart,0,trcrn(:,:,nt_Tsfc,:,:),'ruf8', &
              'Tsfcn',ncat,diag,field_loc_center, field_type_scalar)

      if (my_task == master_task) &
         write(nu_diag,*) 'min/max sice for each layer'
      do k=1,nilyr
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart,0,trcrn(:,:,nt_sice+k-1,:,:),'ruf8', &
              'sice'//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo

      if (my_task == master_task) &
         write(nu_diag,*) 'min/max qice for each layer'
      do k=1,nilyr
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart,0,trcrn(:,:,nt_qice+k-1,:,:),'ruf8', &
              'qice'//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo

      if (my_task == master_task) &
         write(nu_diag,*) 'min/max qsno for each layer'
      do k=1,nslyr
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart,0,trcrn(:,:,nt_qsno+k-1,:,:),'ruf8', &
              'qsno'//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo

      !-----------------------------------------------------------------
      ! velocity
      !-----------------------------------------------------------------
      if (my_task == master_task) &
           write(nu_diag,*) 'min/max velocity components'

      call read_restart_field(nu_restart,0,uvel,'ruf8', &
           'uvel',1,diag,field_loc_NEcorner, field_type_vector)
      call read_restart_field(nu_restart,0,vvel,'ruf8', &
           'vvel',1,diag,field_loc_NEcorner, field_type_vector)

      if (grid_ice == 'CD') then
         if (query_field(nu_restart,'uvelE')) &
            call read_restart_field(nu_restart,0,uvelE,'ruf8', &
                 'uvelE',1,diag,field_loc_Eface, field_type_vector)
         if (query_field(nu_restart,'vvelE')) &
            call read_restart_field(nu_restart,0,vvelE,'ruf8', &
                 'vvelE',1,diag,field_loc_Eface, field_type_vector)
         if (query_field(nu_restart,'uvelN')) &
            call read_restart_field(nu_restart,0,uvelN,'ruf8', &
                 'uvelN',1,diag,field_loc_Nface, field_type_vector)
         if (query_field(nu_restart,'vvelN')) &
            call read_restart_field(nu_restart,0,vvelN,'ruf8', &
                 'vvelN',1,diag,field_loc_Nface, field_type_vector)
      endif

      if (grid_ice == 'C') then
         if (query_field(nu_restart,'uvelE')) &
            call read_restart_field(nu_restart,0,uvelE,'ruf8', &
                 'uvelE',1,diag,field_loc_Eface, field_type_vector)
         if (query_field(nu_restart,'vvelN')) &
            call read_restart_field(nu_restart,0,vvelN,'ruf8', &
                 'vvelN',1,diag,field_loc_Nface, field_type_vector)
      endif

      !-----------------------------------------------------------------
      ! radiation fields
      !-----------------------------------------------------------------

      if (my_task == master_task) &
         write(nu_diag,*) 'radiation fields'

      if (restart_coszen) call read_restart_field(nu_restart,0,coszen,'ruf8', &
           'coszen',1,diag)
      call read_restart_field(nu_restart,0,scale_factor,'ruf8', &
           'scale_factor',1,diag, field_loc_center, field_type_scalar)
      call read_restart_field(nu_restart,0,swvdr,'ruf8', &
           'swvdr',1,diag,field_loc_center, field_type_scalar)
      call read_restart_field(nu_restart,0,swvdf,'ruf8', &
           'swvdf',1,diag,field_loc_center, field_type_scalar)
      call read_restart_field(nu_restart,0,swidr,'ruf8', &
           'swidr',1,diag,field_loc_center, field_type_scalar)
      call read_restart_field(nu_restart,0,swidf,'ruf8', &
           'swidf',1,diag,field_loc_center, field_type_scalar)

      !-----------------------------------------------------------------
      ! ocean stress
      !-----------------------------------------------------------------
      if (my_task == master_task) &
           write(nu_diag,*) 'min/max ocean stress components'

      call read_restart_field(nu_restart,0,strocnxT_iavg,'ruf8', &
           'strocnxT',1,diag,field_loc_center, field_type_vector)
      call read_restart_field(nu_restart,0,strocnyT_iavg,'ruf8', &
           'strocnyT',1,diag,field_loc_center, field_type_vector)

      !-----------------------------------------------------------------
      ! internal stress
      ! The stress tensor must be read and scattered in pairs in order
      ! to properly match corner values across a tripole grid cut.
      !-----------------------------------------------------------------
      if (my_task == master_task) write(nu_diag,*) &
           'internal stress components'

      call read_restart_field(nu_restart,0,stressp_1,'ruf8', &
           'stressp_1',1,diag,field_loc_center,field_type_scalar) ! stressp_1
      call read_restart_field(nu_restart,0,stressp_3,'ruf8', &
           'stressp_3',1,diag,field_loc_center,field_type_scalar) ! stressp_3
      call read_restart_field(nu_restart,0,stressp_2,'ruf8', &
           'stressp_2',1,diag,field_loc_center,field_type_scalar) ! stressp_2
      call read_restart_field(nu_restart,0,stressp_4,'ruf8', &
           'stressp_4',1,diag,field_loc_center,field_type_scalar) ! stressp_4

      call read_restart_field(nu_restart,0,stressm_1,'ruf8', &
           'stressm_1',1,diag,field_loc_center,field_type_scalar) ! stressm_1
      call read_restart_field(nu_restart,0,stressm_3,'ruf8', &
           'stressm_3',1,diag,field_loc_center,field_type_scalar) ! stressm_3
      call read_restart_field(nu_restart,0,stressm_2,'ruf8', &
           'stressm_2',1,diag,field_loc_center,field_type_scalar) ! stressm_2
      call read_restart_field(nu_restart,0,stressm_4,'ruf8', &
           'stressm_4',1,diag,field_loc_center,field_type_scalar) ! stressm_4

      call read_restart_field(nu_restart,0,stress12_1,'ruf8', &
           'stress12_1',1,diag,field_loc_center,field_type_scalar) ! stress12_1
      call read_restart_field(nu_restart,0,stress12_3,'ruf8', &
           'stress12_3',1,diag,field_loc_center,field_type_scalar) ! stress12_1

      call read_restart_field(nu_restart,0,stress12_2,'ruf8', &
           'stress12_2',1,diag,field_loc_center,field_type_scalar) ! stress12_2
      call read_restart_field(nu_restart,0,stress12_4,'ruf8', &
           'stress12_4',1,diag,field_loc_center,field_type_scalar) ! stress12_4

      if (grid_ice == 'CD' .or. grid_ice == 'C') then
         if (query_field(nu_restart,'stresspT')) &
            call read_restart_field(nu_restart,0,stresspT,'ruf8', &
               'stresspT' ,1,diag,field_loc_center,field_type_scalar) ! stresspT
         if (query_field(nu_restart,'stressmT')) &
            call read_restart_field(nu_restart,0,stressmT,'ruf8', &
               'stressmT' ,1,diag,field_loc_center,field_type_scalar) ! stressmT
         if (query_field(nu_restart,'stress12T')) &
            call read_restart_field(nu_restart,0,stress12T,'ruf8', &
               'stress12T',1,diag,field_loc_center,field_type_scalar) ! stress12T
         if (query_field(nu_restart,'stresspU')) &
            call read_restart_field(nu_restart,0,stresspU,'ruf8', &
               'stresspU' ,1,diag,field_loc_NEcorner,field_type_scalar) ! stresspU
         if (query_field(nu_restart,'stressmU')) &
            call read_restart_field(nu_restart,0,stressmU,'ruf8', &
               'stressmU' ,1,diag,field_loc_NEcorner,field_type_scalar) ! stressmU
         if (query_field(nu_restart,'stress12U')) &
            call read_restart_field(nu_restart,0,stress12U,'ruf8', &
               'stress12U',1,diag,field_loc_NEcorner,field_type_scalar) ! stress12U
      endif

      if (trim(ns_boundary_type) == 'tripole' .or. &
          trim(ns_boundary_type) == 'tripoleT') then
         call ice_HaloUpdate_stress(stressp_1, stressp_3, halo_info, &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_3, stressp_1, halo_info, &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_2, stressp_4, halo_info, &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_4, stressp_2, halo_info, &
                                    field_loc_center,  field_type_scalar)

         call ice_HaloUpdate_stress(stressm_1, stressm_3, halo_info, &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_3, stressm_1, halo_info, &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_2, stressm_4, halo_info, &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_4, stressm_2, halo_info, &
                                    field_loc_center,  field_type_scalar)

         call ice_HaloUpdate_stress(stress12_1, stress12_3, halo_info, &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_3, stress12_1, halo_info, &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_2, stress12_4, halo_info, &
                                    field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_4, stress12_2, halo_info, &
                                    field_loc_center,  field_type_scalar)
         ! TODO: CD-grid
      endif

      !-----------------------------------------------------------------
      ! ice mask for dynamics
      !-----------------------------------------------------------------
      if (kdyn > 0) then

         if (my_task == master_task) &
             write(nu_diag,*) 'ice mask for dynamics'
         if (query_field(nu_restart,'iceumask')) then
            call read_restart_field(nu_restart,0,work1,'ruf8', &
                'iceumask',1,diag,field_loc_center, field_type_scalar)

            iceUmask(:,:,:) = .false.
            !$OMP PARALLEL DO PRIVATE(iblk,i,j)
            do iblk = 1, nblocks
               do j = 1, ny_block
               do i = 1, nx_block
                  if (work1(i,j,iblk) > p5) iceUmask(i,j,iblk) = .true.
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO
         endif
         if (grid_ice == 'CD' .or. grid_ice == 'C') then

            if (query_field(nu_restart,'icenmask')) then
               call read_restart_field(nu_restart,0,work1,'ruf8', &
                    'icenmask',1,diag,field_loc_center, field_type_scalar)

               iceNmask(:,:,:) = .false.
               !$OMP PARALLEL DO PRIVATE(iblk,i,j)
               do iblk = 1, nblocks
                  do j = 1, ny_block
                  do i = 1, nx_block
                     if (work1(i,j,iblk) > p5) iceNmask(i,j,iblk) = .true.
                  enddo
                  enddo
               enddo
               !$OMP END PARALLEL DO
            endif

            if (query_field(nu_restart,'iceemask')) then
               call read_restart_field(nu_restart,0,work1,'ruf8', &
                    'iceemask',1,diag,field_loc_center, field_type_scalar)

               iceEmask(:,:,:) = .false.
               !$OMP PARALLEL DO PRIVATE(iblk,i,j)
               do iblk = 1, nblocks
                  do j = 1, ny_block
                  do i = 1, nx_block
                     if (work1(i,j,iblk) > p5) iceEmask(i,j,iblk) = .true.
                  enddo
                  enddo
               enddo
               !$OMP END PARALLEL DO
            endif
         endif
      else
         if (my_task == master_task) &
            write(nu_diag,*) 'ice mask for dynamics - not used, however mandatory to read in binary files'
         if (query_field(nu_restart,'iceumask')) then
            call read_restart_field(nu_restart,0,work1,'ruf8', &
               'iceumask',1,diag,field_loc_center, field_type_scalar)
         endif
         if (grid_ice == 'CD' .or. grid_ice == 'C') then
            if (query_field(nu_restart,'icenmask')) then
               call read_restart_field(nu_restart,0,work1,'ruf8', &
                    'icenmask',1,diag,field_loc_center, field_type_scalar)
            endif
            if (query_field(nu_restart,'iceemask')) then
               call read_restart_field(nu_restart,0,work1,'ruf8', &
                    'iceemask',1,diag,field_loc_center, field_type_scalar)
            endif
         endif
      endif

      ! set Tsfcn to c0 on land
      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            if (.not.tmask(i,j,iblk) .and. .not.opmask(i,j,iblk)) trcrn(i,j,nt_Tsfc,:,iblk) = c0
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      ! for mixed layer model
      if (oceanmixed_ice) then

         if (my_task == master_task) &
              write(nu_diag,*) 'min/max sst, frzmlt'

         call read_restart_field(nu_restart,0,sst,'ruf8', &
              'sst',1,diag,field_loc_center, field_type_scalar)
         call read_restart_field(nu_restart,0,frzmlt,'ruf8', &
              'frzmlt',1,diag,field_loc_center, field_type_scalar)
      endif

      !-----------------------------------------------------------------
      ! Ensure unused stress values in west and south ghost cells are 0
      !-----------------------------------------------------------------
      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
         do j = 1, nghost
         do i = 1, nx_block
            stressp_1 (i,j,iblk) = c0
            stressp_2 (i,j,iblk) = c0
            stressp_3 (i,j,iblk) = c0
            stressp_4 (i,j,iblk) = c0
            stressm_1 (i,j,iblk) = c0
            stressm_2 (i,j,iblk) = c0
            stressm_3 (i,j,iblk) = c0
            stressm_4 (i,j,iblk) = c0
            stress12_1(i,j,iblk) = c0
            stress12_2(i,j,iblk) = c0
            stress12_3(i,j,iblk) = c0
            stress12_4(i,j,iblk) = c0
         enddo
         enddo
         do j = 1, ny_block
         do i = 1, nghost
            stressp_1 (i,j,iblk) = c0
            stressp_2 (i,j,iblk) = c0
            stressp_3 (i,j,iblk) = c0
            stressp_4 (i,j,iblk) = c0
            stressm_1 (i,j,iblk) = c0
            stressm_2 (i,j,iblk) = c0
            stressm_3 (i,j,iblk) = c0
            stressm_4 (i,j,iblk) = c0
            stress12_1(i,j,iblk) = c0
            stress12_2(i,j,iblk) = c0
            stress12_3(i,j,iblk) = c0
            stress12_4(i,j,iblk) = c0
         enddo
         enddo
         ! TODO: CD-grid ?
      enddo
      !$OMP END PARALLEL DO

      !-----------------------------------------------------------------
      ! Ensure ice is binned in correct categories
      ! (should not be necessary unless restarting from a run with
      !  different category boundaries).
      !
      ! If called, this subroutine does not give exact restart.
      !-----------------------------------------------------------------
!!!      call cleanup_itd

      !-----------------------------------------------------------------
      ! compute aggregate ice state and open water area
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks

      do j = 1, ny_block
      do i = 1, nx_block
         if (tmask(i,j,iblk) .or. opmask(i,j,iblk)) &
            call icepack_aggregate(aicen = aicen(i,j,:,iblk),     &
                                   trcrn = trcrn(i,j,:,:,iblk),   &
                                   vicen = vicen(i,j,:,iblk),     &
                                   vsnon = vsnon(i,j,:,iblk),     &
                                   aice  = aice (i,j,  iblk),     &
                                   trcr  = trcr (i,j,:,iblk),     &
                                   vice  = vice (i,j,  iblk),     &
                                   vsno  = vsno (i,j,  iblk),     &
                                   aice0 = aice0(i,j,  iblk),     &
                                   trcr_depend   = trcr_depend,   &
                                   trcr_base     = trcr_base,     &
                                   n_trcr_strata = n_trcr_strata, &
                                   nt_strata     = nt_strata,     &
                                   Tf            = Tf(i,j,iblk))

         aice_init(i,j,iblk) = aice(i,j,iblk)
      enddo
      enddo

      enddo
      !$OMP END PARALLEL DO

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      ! if runid is bering then need to correct npt for istep0
      if (trim(runid) == 'bering') then
         npt = npt - istep0
      endif

      !-----------------------------------------------------------------
      ! update concentration from a file
      !-----------------------------------------------------------------
      if (restart_mod /= "none") then

         select case (trim(restart_mod))

         case('adjust_aice')
            call direct_adjust_aice

         case('adjust_aice_test')
            call direct_adjust_aice(test=.true.)

         case default
            call abort_ice(subname//'ERROR: unsupported restart_mod='//trim(restart_mod), &
               file=__FILE__, line=__LINE__)

         end select

         !-----------------------------------------------------------------
         ! Ensure ice is binned in correct categories
         !-----------------------------------------------------------------
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               if (tmask(i,j,iblk)) then

                  call cleanup_itd(dt,    hin_max,             &
                       aicen(i,j,:,iblk), trcrn(i,j,:,:,iblk), &
                       vicen(i,j,:,iblk), vsnon(i,j,:,  iblk), &
                       aice0(i,j,  iblk),  aice(i,j,    iblk), &
                       tr_aero,           tr_pond_topo,        &
                       first_ice(i,j,:,iblk),                  &
                       trcr_depend,       trcr_base,           &
                       n_trcr_strata,     nt_strata,           &
                       fpond     =     fpond(i,j,  iblk),      &
                       fresh     =     fresh(i,j,  iblk),      &
                       fsalt     =     fsalt(i,j,  iblk),      &
                       fhocn     =     fhocn(i,j,  iblk),      &
                       faero_ocn = faero_ocn(i,j,:,iblk),      &
                       fiso_ocn  =  fiso_ocn(i,j,:,iblk),      &
                       flux_bio  =  flux_bio(i,j,:,iblk),      &
                       Tf        =        Tf(i,j,  iblk),      &
                       limit_aice = .true. )

                  call icepack_aggregate( &
                                   aicen = aicen(i,j,:,iblk),     &
                                   trcrn = trcrn(i,j,:,:,iblk),   &
                                   vicen = vicen(i,j,:,iblk),     &
                                   vsnon = vsnon(i,j,:,iblk),     &
                                   aice  = aice (i,j,  iblk),     &
                                   trcr  = trcr (i,j,:,iblk),     &
                                   vice  = vice (i,j,  iblk),     &
                                   vsno  = vsno (i,j,  iblk),     &
                                   aice0 = aice0(i,j,  iblk),     &
                                   trcr_depend   = trcr_depend,   &
                                   trcr_base     = trcr_base,     &
                                   n_trcr_strata = n_trcr_strata, &
                                   nt_strata     = nt_strata,     &
                                   Tf            = Tf(i,j,iblk))

                  aice_init(i,j,iblk) = aice(i,j,iblk)

               endif  ! tmask
            enddo     ! i
            enddo     ! j
         enddo        ! iblk

      endif !restart_mod

      end subroutine restartfile

!=======================================================================

! Restarts from a CICE v4.1 (binary) dump
! author Elizabeth C. Hunke, LANL

      subroutine restartfile_v4 (ice_ic)

      use ice_broadcast, only: broadcast_scalar
      use ice_blocks, only: nghost, nx_block, ny_block
      use ice_calendar, only: istep0, istep1, timesecs, calendar, npt, &
          set_date_from_timesecs
      use ice_domain, only: nblocks, distrb_info
      use ice_domain_size, only: nilyr, nslyr, ncat, nx_global, ny_global, &
          max_blocks
      use ice_dyn_shared, only: iceUmask
      use ice_flux, only: scale_factor, swvdr, swvdf, swidr, swidf, &
          strocnxT_iavg, strocnyT_iavg, sst, frzmlt, Tf, &
          stressp_1, stressp_2, stressp_3, stressp_4, &
          stressm_1, stressm_2, stressm_3, stressm_4, &
          stress12_1, stress12_2, stress12_3, stress12_4
      use ice_gather_scatter, only: scatter_global_stress
      use ice_read_write, only: ice_open, ice_read, ice_read_global
      use ice_state, only: trcr_depend, aice, vice, vsno, trcr, &
          aice0, aicen, vicen, vsnon, trcrn, aice_init, uvel, vvel, &
          trcr_base, nt_strata, n_trcr_strata

      character (*), optional :: ice_ic

      ! local variables

      integer (kind=int_kind) :: &
         i, j, k, n, iblk, &     ! counting indices
         ntrcr, &                ! number of tracers
         nt_Tsfc, nt_sice, nt_qice, nt_qsno, &
         iignore                 ! dummy variable

      real (kind=real_kind) :: &
         rignore                 ! dummy variable

      character(len=char_len_long) :: &
         filename, filename0

      logical (kind=log_kind) :: &
         diag

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1

      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         work_g1, work_g2

      real (kind=dbl_kind) :: &
         time_forc      ! historic, now local

      character(len=*), parameter :: subname = '(restartfile_v4)'

      call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__, line=__LINE__)

      call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_sice_out=nt_sice, &
           nt_qice_out=nt_qice, nt_qsno_out=nt_qsno)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (present(ice_ic)) then
         filename = ice_ic
      elseif (my_task == master_task) then
         open(nu_rst_pointer,file=pointer_file)
         read(nu_rst_pointer,'(a)') filename0
         filename = trim(filename0)
         close(nu_rst_pointer)
         write(nu_diag,*) 'Read ',pointer_file(1:lenstr(pointer_file))
      endif

      call ice_open(nu_restart,filename,0)

      if (my_task == master_task) &
         write(nu_diag,*) 'Using restart dump=', trim(filename)

      if (use_restart_time) then

         if (my_task == master_task) then
            read (nu_restart) istep0,timesecs,time_forc
            write(nu_diag,*) 'Restart read at istep=',istep0,timesecs
         endif
         call broadcast_scalar(istep0,master_task)
         istep1 = istep0
         call broadcast_scalar(timesecs,master_task)
!         call broadcast_scalar(time_forc,master_task)
         call set_date_from_timesecs(timesecs)
         call calendar()

      else

         if (my_task == master_task) &
            read (nu_restart) iignore,rignore,rignore

      endif

      diag = .true.     ! write min/max diagnostics for field

      !-----------------------------------------------------------------
      ! state variables
      ! Tsfc is the only tracer read in this file.  All other
      ! tracers are in their own dump/restart files.
      !-----------------------------------------------------------------
      do n=1,ncat
         if (my_task == master_task) &
              write(nu_diag,*) 'cat ',n, &
                               ' min/max area, vol ice, vol snow, Tsfc'

         call ice_read(nu_restart,0,aicen(:,:,n,:),'ruf8',diag, &
                       field_loc_center, field_type_scalar)
         call ice_read(nu_restart,0,vicen(:,:,n,:),'ruf8',diag, &
                       field_loc_center, field_type_scalar)
         call ice_read(nu_restart,0,vsnon(:,:,n,:),'ruf8',diag, &
                       field_loc_center, field_type_scalar)
         call ice_read(nu_restart,0,trcrn(:,:,nt_Tsfc,n,:),'ruf8',diag, &
                       field_loc_center, field_type_scalar)

         if (my_task == master_task) &
              write(nu_diag,*) 'cat ',n, 'min/max sice for each layer'
         do k=1,nilyr
            call ice_read(nu_restart,0,trcrn(:,:,nt_sice+k-1,n,:),'ruf8',diag, &
                          field_loc_center, field_type_scalar)
         enddo

         if (my_task == master_task) &
              write(nu_diag,*) 'cat ',n, 'min/max qice for each layer'
         do k=1,nilyr
            call ice_read(nu_restart,0,trcrn(:,:,nt_qice+k-1,n,:),'ruf8',diag, &
                          field_loc_center, field_type_scalar)
         enddo

         if (my_task == master_task) &
              write(nu_diag,*) 'cat ',n, 'min/max qsno for each layer'
         do k=1,nslyr
            call ice_read(nu_restart,0,trcrn(:,:,nt_qsno+k-1,n,:),'ruf8',diag, &
                          field_loc_center, field_type_scalar)
         enddo
      enddo ! ncat

      !-----------------------------------------------------------------
      ! velocity
      !-----------------------------------------------------------------
      if (my_task == master_task) &
           write(nu_diag,*) 'min/max velocity components'

      call ice_read(nu_restart,0,uvel,'ruf8',diag, &
                       field_loc_NEcorner, field_type_vector)
      call ice_read(nu_restart,0,vvel,'ruf8',diag, &
                       field_loc_NEcorner, field_type_vector)

      !-----------------------------------------------------------------
      ! radiation fields
      !-----------------------------------------------------------------

      if (my_task == master_task) &
         write(nu_diag,*) 'radiation fields'

      call ice_read(nu_restart,0,scale_factor,'ruf8',diag, &
                    field_loc_center, field_type_scalar)
      call ice_read(nu_restart,0,swvdr,'ruf8',diag, &
                    field_loc_center, field_type_scalar)
      call ice_read(nu_restart,0,swvdf,'ruf8',diag, &
                    field_loc_center, field_type_scalar)
      call ice_read(nu_restart,0,swidr,'ruf8',diag, &
                    field_loc_center, field_type_scalar)
      call ice_read(nu_restart,0,swidf,'ruf8',diag, &
                    field_loc_center, field_type_scalar)

      !-----------------------------------------------------------------
      ! ocean stress
      !-----------------------------------------------------------------
      if (my_task == master_task) &
           write(nu_diag,*) 'min/max ocean stress components'

      call ice_read(nu_restart,0,strocnxT_iavg,'ruf8',diag, &
                       field_loc_center, field_type_vector)
      call ice_read(nu_restart,0,strocnyT_iavg,'ruf8',diag, &
                       field_loc_center, field_type_vector)

      !-----------------------------------------------------------------
      ! internal stress
      ! The stress tensor must be read and scattered in pairs in order
      ! to properly match corner values across a tripole grid cut.
      !-----------------------------------------------------------------
      if (my_task == master_task) write(nu_diag,*) &
           'internal stress components'

      allocate (work_g1(nx_global,ny_global), &
                work_g2(nx_global,ny_global))

      call ice_read_global(nu_restart,0,work_g1,'ruf8',diag) ! stressp_1
      call ice_read_global(nu_restart,0,work_g2,'ruf8',diag) ! stressp_3
      call scatter_global_stress(stressp_1, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stressp_3, work_g2, work_g1, &
                                 master_task, distrb_info)

      call ice_read_global(nu_restart,0,work_g1,'ruf8',diag) ! stressp_2
      call ice_read_global(nu_restart,0,work_g2,'ruf8',diag) ! stressp_4
      call scatter_global_stress(stressp_2, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stressp_4, work_g2, work_g1, &
                                 master_task, distrb_info)

      call ice_read_global(nu_restart,0,work_g1,'ruf8',diag) ! stressm_1
      call ice_read_global(nu_restart,0,work_g2,'ruf8',diag) ! stressm_3
      call scatter_global_stress(stressm_1, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stressm_3, work_g2, work_g1, &
                                 master_task, distrb_info)

      call ice_read_global(nu_restart,0,work_g1,'ruf8',diag) ! stressm_2
      call ice_read_global(nu_restart,0,work_g2,'ruf8',diag) ! stressm_4
      call scatter_global_stress(stressm_2, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stressm_4, work_g2, work_g1, &
                                 master_task, distrb_info)

      call ice_read_global(nu_restart,0,work_g1,'ruf8',diag) ! stress12_1
      call ice_read_global(nu_restart,0,work_g2,'ruf8',diag) ! stress12_3
      call scatter_global_stress(stress12_1, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stress12_3, work_g2, work_g1, &
                                 master_task, distrb_info)

      call ice_read_global(nu_restart,0,work_g1,'ruf8',diag) ! stress12_2
      call ice_read_global(nu_restart,0,work_g2,'ruf8',diag) ! stress12_4
      call scatter_global_stress(stress12_2, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stress12_4, work_g2, work_g1, &
                                 master_task, distrb_info)

      deallocate (work_g1, work_g2)

      !-----------------------------------------------------------------
      ! ice mask for dynamics
      !-----------------------------------------------------------------
      if (my_task == master_task) &
           write(nu_diag,*) 'ice mask for dynamics'

      call ice_read(nu_restart,0,work1,'ruf8',diag, &
                    field_loc_center, field_type_scalar)

      iceUmask(:,:,:) = .false.
      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            if (work1(i,j,iblk) > p5) iceUmask(i,j,iblk) = .true.
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      ! for mixed layer model
      if (oceanmixed_ice) then

         if (my_task == master_task) &
              write(nu_diag,*) 'min/max sst, frzmlt'

         call ice_read(nu_restart,0,sst,'ruf8',diag, &
                       field_loc_center, field_type_scalar)
         call ice_read(nu_restart,0,frzmlt,'ruf8',diag, &
                       field_loc_center, field_type_scalar)
      endif

      if (my_task == master_task) close(nu_restart)

      !-----------------------------------------------------------------
      ! Ensure unused stress values in west and south ghost cells are 0
      !-----------------------------------------------------------------
      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
         do j = 1, nghost
         do i = 1, nx_block
            stressp_1 (i,j,iblk) = c0
            stressp_2 (i,j,iblk) = c0
            stressp_3 (i,j,iblk) = c0
            stressp_4 (i,j,iblk) = c0
            stressm_1 (i,j,iblk) = c0
            stressm_2 (i,j,iblk) = c0
            stressm_3 (i,j,iblk) = c0
            stressm_4 (i,j,iblk) = c0
            stress12_1(i,j,iblk) = c0
            stress12_2(i,j,iblk) = c0
            stress12_3(i,j,iblk) = c0
            stress12_4(i,j,iblk) = c0
         enddo
         enddo
         do j = 1, ny_block
         do i = 1, nghost
            stressp_1 (i,j,iblk) = c0
            stressp_2 (i,j,iblk) = c0
            stressp_3 (i,j,iblk) = c0
            stressp_4 (i,j,iblk) = c0
            stressm_1 (i,j,iblk) = c0
            stressm_2 (i,j,iblk) = c0
            stressm_3 (i,j,iblk) = c0
            stressm_4 (i,j,iblk) = c0
            stress12_1(i,j,iblk) = c0
            stress12_2(i,j,iblk) = c0
            stress12_3(i,j,iblk) = c0
            stress12_4(i,j,iblk) = c0
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !-----------------------------------------------------------------
      ! Ensure ice is binned in correct categories
      ! (should not be necessary unless restarting from a run with
      !  different category boundaries).
      !
      ! If called, this subroutine does not give exact restart.
      !-----------------------------------------------------------------
!!!      call cleanup_itd

      !-----------------------------------------------------------------
      ! compute aggregate ice state and open water area
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks

      do j = 1, ny_block
      do i = 1, nx_block
         if (tmask(i,j,iblk) .or. opmask(i,j,iblk)) &
            call icepack_aggregate(aicen = aicen(i,j,:,iblk),     &
                                   trcrn = trcrn(i,j,:,:,iblk),   &
                                   vicen = vicen(i,j,:,iblk),     &
                                   vsnon = vsnon(i,j,:,iblk),     &
                                   aice  = aice (i,j,  iblk),     &
                                   trcr  = trcr (i,j,:,iblk),     &
                                   vice  = vice (i,j,  iblk),     &
                                   vsno  = vsno (i,j,  iblk),     &
                                   aice0 = aice0(i,j,  iblk),     &
                                   trcr_depend   = trcr_depend,   &
                                   trcr_base     = trcr_base,     &
                                   n_trcr_strata = n_trcr_strata, &
                                   nt_strata     = nt_strata,     &
                                   Tf            = Tf(i,j,iblk))

         aice_init(i,j,iblk) = aice(i,j,iblk)
      enddo
      enddo

      enddo
      !$OMP END PARALLEL DO

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      ! creates new file
      filename = trim(restart_dir) // '/iced.converted'
      call dumpfile(filename)
      call final_restart
      ! stop

      ! if runid is bering then need to correct npt for istep0
      if (trim(runid) == 'bering') then
         npt = npt - istep0
      endif

      end subroutine restartfile_v4

!=======================================================================

!=======================================================================
!  Direct insertion of ice concentration read from file.
!
!  Posey, et. al. 2015: Improving Arctic sea ice edge forecasts by
!  assimilating high horizontal resolution sea ice concentration
!  data into the US Navy's ice forecast systems.
!  The Cryosphere.  doi:10.5194/tc-9-1735-2015
!
!  Alan J. Wallcraft, COAPS/FSU, Nov 2024

      subroutine direct_adjust_aice(test)

        use ice_blocks, only: nghost, nx_block, ny_block
        use ice_domain, only: nblocks
        use ice_domain_size, only: nilyr, nslyr, ncat, max_blocks
        use ice_grid, only: tmask
        use ice_communicate, only: my_task, master_task
        use ice_constants, only: c0, c1, c4, c20, c100, &
             p5, p2, p1, p01, p001, &
             field_loc_center, field_loc_NEcorner, &
             field_type_scalar, field_type_vector
        use ice_fileunits, only: nu_diag
        use ice_flux, only:  &
             Tair, Tf, salinz, Tmltz, sst,                  &
             stressp_1, stressp_2, stressp_3, stressp_4,    &
             stressm_1, stressm_2, stressm_3, stressm_4,    &
             stress12_1, stress12_2, stress12_3, stress12_4
        use ice_state, only:  &
             aice, aicen, vicen, vsnon, trcrn
        use ice_read_write, only: ice_check_nc, ice_read_nc, &
            ice_open_nc, ice_close_nc
        use ice_arrays_column, only: hin_max
        ! use icepack_mushy_physics, only: enthalpy_mush
        use icepack_intfc, only: icepack_init_trcr

        logical(kind=log_kind), optional, intent(in) :: &
             test         ! use internally generated aice

        ! --- local variables
        real(kind=dbl_kind) :: &
             q     ,    & ! scale factor
             aice_m,    & ! model aice
             aice_o,    & ! observation aice
             aice_t,    & ! target aice
             aice_i,    & ! insert ice
             slope,     & ! used to compute surf Temp
             Ti,        & ! target surface temperature
             edge_om,   & ! nominal ice edge zone
             diff_om,   & ! allowed model vs obs difference
             hin_om,    & ! new ice thickness
             aicen_old, & ! old value of aice to check when adding ice
             vsnon_old, & ! old value of snow volume to check when adding ice
             Tsfc         ! surface temp.
        integer (kind=int_kind) :: &
             fid          ! file id for netCDF file
        integer (kind=int_kind) :: &
             i, j, k, n, iblk ! counting  indices
        logical (kind=log_kind) :: &
             diag,      & ! diagnostic output
             ltest        ! local value of test argument
        real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
             work1
        real (kind=dbl_kind), dimension(nilyr) :: &
             qin          ! ice enthalpy (J/m3)
        real (kind=dbl_kind), dimension(nslyr) :: &
             qsn          ! snow enthalpy (J/m3)
        character(len=char_len_long) :: &
             aice_filename, &! filename to read in
             aice_fldname    ! fieldname to read in

        ! parameters from icepack
        real (kind=dbl_kind) :: &
             puny, Tffresh, Tsmelt, Lfresh, cp_ice, cp_ocn, &
             rhos, rhoi
        integer (kind=int_kind) :: &
             nt_Tsfc, nt_sice, nt_qice, nt_qsno, &
             ktherm
        character(len=*), parameter :: subname = '(direct_adjust_aice)'

        diag = .true.
        ltest = .false.
        if (present(test)) then
           ltest = test
        endif
        aice_filename = trim(restart_dir)//'/sic.nc'
        aice_fldname = 'sic'

        ! get parameters from icepack
        call icepack_query_parameters( &
             puny_out=puny,            &
             Tffresh_out=Tffresh,      &
             Tsmelt_out=Tsmelt,        &
             Lfresh_out=Lfresh,        &
             cp_ice_out=cp_ice,        &
             cp_ocn_out=cp_ocn,        &
             rhos_out=rhos,            &
             rhoi_out=rhoi,            &
             ktherm_out=ktherm         )

        call icepack_query_tracer_indices( &
             nt_Tsfc_out=nt_Tsfc,          &
             nt_sice_out=nt_sice,          &
             nt_qice_out=nt_qice,          &
             nt_qsno_out=nt_qsno           )

        call icepack_warnings_flush(nu_diag)
        if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
             file=__FILE__, line=__LINE__)

        if (ltest) then
           if (my_task == master_task) then
              write(nu_diag,*) subname//" direct_adjust_aice rounding to nearest 1/20th"
           endif
           work1 = nint(aice*c20)/c20   ! round to nearest 5/100th
        else
           if (my_task == master_task) then
              write(nu_diag,*) subname//" direct_adjust_aice from "//trim(aice_filename)
           endif

           call ice_open_nc(trim(aice_filename), fid)
           call ice_read_nc(fid,1,trim(aice_fldname),work1,diag, &
                            field_loc=field_loc_center, &
                            field_type=field_type_scalar)
           call ice_close_nc(fid)
        endif

        edge_om = p2  ! nominal ice edge zone
        diff_om = p1  ! allowed model vs obs difference
        hin_om  = hin_max(1)*0.9_dbl_kind  !new ice thickness

        do iblk = 1, nblocks
           do j = 1, ny_block
           do i = 1, nx_block

              aice_o = work1(i,j,iblk) ! obs.  ice concentration
              aice_m = aice(i,j,iblk)  ! model ice concentration

              if     (.not.tmask(i,j,iblk)) then
                 ! land - do nothing
              elseif (aice_o.gt.p01 .and. &
                   abs(aice_o-aice_m).le.p01) then
                 ! model and obs are very close - do nothing
              elseif (min(aice_o,aice_m).ge.edge_om .and. &
                   abs(aice_o-aice_m).le.diff_om) then
                 ! model and obs are close enough - do nothing
              elseif (aice_o.eq.aice_m) then
              elseif (aice_o.lt.aice_m) then
                 if (aice_o.lt.p01)then
                    ! --- remove all ice ---
                    ! warm sst so the ice won't grow immediately
                    sst(i,j,iblk) = sst(i,j,iblk) + p2
                    do n=1,ncat
                       aicen(i,j,n,iblk) = c0
                       vicen(i,j,n,iblk) = c0
                       vsnon(i,j,n,iblk) = c0
                       call icepack_init_trcr( &
                            Tair     =   Tair(i,j,  iblk), &
                            Tf       =     Tf(i,j,  iblk), &
                            Sprofile = salinz(i,j,:,iblk), &
                            Tprofile =  Tmltz(i,j,:,iblk), &
                            Tsfc     = Tsfc,               &
                            qin      = qin(:),             &
                            qsn      = qsn(:)              )
                       ! surface temperature
                       trcrn(i,j,nt_Tsfc,n,iblk) = Tsfc ! deg C
                       ! ice enthalpy, salinity
                       do k = 1, nilyr
                          trcrn(i,j,nt_qice+k-1,n,iblk) = qin(k)
                          trcrn(i,j,nt_sice+k-1,n,iblk) = salinz(i,j,k,iblk)
                       enddo  ! nilyr
                       ! snow enthalpy
                       do k = 1, nslyr
                          trcrn(i,j,nt_qsno+k-1,n,iblk) = qsn(k)
                       enddo  ! nslyr
                    enddo !n
                 else !aice_o.ge.p01
                    if     (aice_o.lt.edge_om) then
                       ! --- target ice conc. is obs.
                       aice_t = aice_o
                    else !aice_m-aice_o.gt.diff_om
                       ! --- target ice conc. is obs.+diff_om
                       aice_t = aice_o + diff_om
                    endif
                    ! --- reduce ice to the target concentration,
                    !     completely exhasting ice categories in order ---
                    aice_i = aice_m - aice_t   !>=0.0
                    do n=1,ncat
                       if     (aice_i.le.p001) then
                          exit
                       elseif (aice_i.ge.aicen(i,j,n,iblk)) then
                          ! --- remove all of this category
                          aice_i = aice_i - aicen(i,j,n,iblk)
                          aicen(i,j,n,iblk) = c0
                          vicen(i,j,n,iblk) = c0
                          vsnon(i,j,n,iblk) = c0
                          call icepack_init_trcr( &
                               Tair     =   Tair(i,j,  iblk), &
                               Tf       =     Tf(i,j,  iblk), &
                               Sprofile = salinz(i,j,:,iblk), &
                               Tprofile =  Tmltz(i,j,:,iblk), &
                               Tsfc     = Tsfc,               &
                               qin      = qin(:),             &
                               qsn      = qsn(:)              )
                          ! surface temperature
                          trcrn(i,j,nt_Tsfc,n,iblk) = Tsfc ! deg C
                          ! ice enthalpy, salinity
                          do k = 1, nilyr
                             trcrn(i,j,nt_qice+k-1,n,iblk) = qin(k)
                             trcrn(i,j,nt_sice+k-1,n,iblk) = salinz(i,j,k,iblk)
                          enddo  ! nilyr
                          ! snow enthalpy
                          do k = 1, nslyr
                             trcrn(i,j,nt_qsno+k-1,n,iblk) = qsn(k)
                          enddo  ! nslyr
                       else  !aice_i.lt.aicen(i,j,n,iblk)
                          ! --- remove part of this category
                          q = (aicen(i,j,n,iblk) - aice_i) &
                               /aicen(i,j,n,iblk)              !<1
                          aice_i = c0

                          ! reduce aicen, vicen, vsnon by q
                          ! do not alter Tsfc since there is already
                          ! ice here.
                          aicen(i,j,n,iblk) = q*aicen(i,j,n,iblk)
                          vicen(i,j,n,iblk) = q*vicen(i,j,n,iblk)
                          vsnon(i,j,n,iblk) = q*vsnon(i,j,n,iblk)
                       endif    ! aice_i.gt.p001 and aice_i.lt.aicen
                    enddo       ! n
                 endif          ! aice_o.lt.p01
              elseif (aice_o.gt.p01) then  ! .and. aice_o.gt.aicen
                 if     (aice_m.lt.edge_om) then
                    ! --- target ice conc. is obs.
                    aice_t = aice_o
                 else !aice_o-aice_m.gt.diff_om
                    ! --- target ice conc. is obs.-diff_om
                    aice_t = aice_o - diff_om
                 endif
                 q = (aice_t-aice_m)
                 ! --- add ice to the target concentration,
                 ! --- with all new ice in category 1
                 ! --- cool sst so the ice won't melt immediately
                 sst(  i,j,  iblk) = sst(  i,j,  iblk) - q  ! 0 <= q <= 1
                 aicen_old         = aicen(i,j,1,iblk)  ! store to check for zero ice later
                 vsnon_old         = vsnon(i,j,1,iblk)  ! store to check for zero snow later
                 aicen(i,j,1,iblk) = aicen(i,j,1,iblk) + q
                 vicen(i,j,1,iblk) = vicen(i,j,1,iblk) + q*hin_om
                 vsnon(i,j,1,iblk) = vsnon(i,j,1,iblk) + q*hin_om*p2

                 ! ------------------------------------------------------
                 ! check for zero snow in 1st category.
                 ! It is possible that there was ice
                 ! but no snow. This would skip the loop below and an
                 ! error in snow thermo would occur. If snow was zero
                 ! specify enthalpy here
                 ! ------------------------------------------------------
                 if (vsnon_old < puny) then
                    do n=1,1           ! only do 1st category
                       ! --- snow layers
                       trcrn(i,j,nt_Tsfc,n,iblk) =   &      ! Tsfc
                            min(Tsmelt,Tair(i,j,iblk) - Tffresh)
                       Ti = min(c0,trcrn(i,j,nt_Tsfc,n,iblk))
                       do k=1,nslyr
                          trcrn(i,j,nt_qsno+k-1,n,iblk) = -rhos*(Lfresh - cp_ice*Ti)
                       enddo ! k
                    enddo ! n = 1,1
                 endif

                 ! ------------------------------------------------------
                 ! check for zero aice in 1st category.
                 ! if adding to an initially zero ice, we must define
                 ! qice, qsno, sice so thermo does not blow up.
                 ! ------------------------------------------------------
                 if (aicen_old < puny) then
                    do n =1,1  ! only do 1st category
                       call icepack_init_trcr( &
                            Tair     =   Tair(i,j,  iblk), &
                            Tf       =     Tf(i,j,  iblk), &
                            Sprofile = salinz(i,j,:,iblk), &
                            Tprofile =  Tmltz(i,j,:,iblk), &
                            Tsfc     = Tsfc,               &
                            qin      = qin(:),             &
                            qsn      = qsn(:)              )
                       ! surface temperature
                       trcrn(i,j,nt_Tsfc,n,iblk) = Tsfc ! deg C
                       ! ice enthalpy, salinity
                       do k = 1, nilyr
                          trcrn(i,j,nt_qice+k-1,n,iblk) = qin(k)
                          trcrn(i,j,nt_sice+k-1,n,iblk) = salinz(i,j,k,iblk)
                       enddo
                       ! snow enthalpy
                       do k = 1, nslyr
                          trcrn(i,j,nt_qsno+k-1,n,iblk) = qsn(k)
                       enddo               ! nslyr
                    enddo       ! n
                 endif          ! qice == c0
              endif             ! aice_o vs aice_m or tmask
           enddo                ! j
           enddo                ! i
        enddo                   ! iblk

      end subroutine direct_adjust_aice

!=======================================================================

      end module ice_restart_driver

!=======================================================================
