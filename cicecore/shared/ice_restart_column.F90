!=========================================================================
!
! Restart routines for the column package.
!
! author: Elizabeth C. Hunke, LANL
!
! 2014: Moved subroutines from column package modules

      module ice_restart_column

      use ice_kinds_mod
      use ice_communicate, only: my_task, master_task
      use ice_constants, only: c0, c1, p5
      use ice_constants, only: field_loc_center, field_type_scalar
      use ice_domain_size, only: ncat, nblyr
      use ice_restart,only: read_restart_field, write_restart_field
      use ice_exit, only: abort_ice
      use ice_fileunits, only: nu_diag
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_max_algae, icepack_max_doc, &
          icepack_max_don, icepack_max_dic, icepack_max_fe, icepack_max_aero
      use icepack_intfc, only: icepack_query_parameters, &
          icepack_query_tracer_numbers, icepack_query_tracer_flags, &
          icepack_query_tracer_indices

      implicit none

      private
      public ::  write_restart_age,       read_restart_age, &
                 write_restart_FY,        read_restart_FY, &
                 write_restart_lvl,       read_restart_lvl, &
                 write_restart_pond_cesm, read_restart_pond_cesm, &
                 write_restart_pond_lvl,  read_restart_pond_lvl, &
                 write_restart_pond_topo, read_restart_pond_topo, &
                 write_restart_aero,      read_restart_aero, &
                 write_restart_bgc,       read_restart_bgc,  &
                 write_restart_hbrine,    read_restart_hbrine

      logical (kind=log_kind), public :: & 
         restart_age      , & ! if .true., read age tracer restart file
         restart_FY       , & ! if .true., read FY tracer restart file
         restart_lvl      , & ! if .true., read lvl tracer restart file
         restart_pond_cesm, & ! if .true., read meltponds restart file
         restart_pond_lvl , & ! if .true., read meltponds restart file
         restart_pond_topo, & ! if .true., read meltponds restart file
         restart_aero     , & ! if .true., read aerosol tracer restart file
         restart_zsal     , & ! if .true., read Salinity from restart file 
         restart_hbrine   , & ! if .true., read hbrine from restart file
         restart_bgc          ! if .true., read bgc restart file

!=======================================================================

      contains

!=======================================================================

! Dumps all values needed for restarting
! author Elizabeth C. Hunke, LANL

      subroutine write_restart_age()

      use ice_fileunits, only: nu_dump_age
      use ice_state, only: trcrn

      ! local variables

      logical (kind=log_kind) :: diag
      integer (kind=int_kind) :: nt_iage
      character(len=*),parameter :: subname='(write_restart_age)'

      call icepack_query_tracer_indices(nt_iage_out=nt_iage)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      diag = .true.

      !-----------------------------------------------------------------

      call write_restart_field(nu_dump_age,0,trcrn(:,:,nt_iage,:,:),'ruf8', &
                               'iage',ncat,diag)

      end subroutine write_restart_age

!=======================================================================

! Reads all values needed for an ice age restart
! author Elizabeth C. Hunke, LANL

      subroutine read_restart_age()

      use ice_fileunits, only: nu_restart_age
      use ice_state, only: trcrn

      ! local variables

      logical (kind=log_kind) :: &
         diag
      integer (kind=int_kind) :: nt_iage
      character(len=*),parameter :: subname='(read_restart_age)'

      call icepack_query_tracer_indices(nt_iage_out=nt_iage)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      diag = .true.

      if (my_task == master_task) write(nu_diag,*) subname,'min/max age (s)'

      call read_restart_field(nu_restart_age,0,trcrn(:,:,nt_iage,:,:),'ruf8', &
                       'iage',ncat,diag,field_loc_center,field_type_scalar)

      end subroutine read_restart_age

!=======================================================================

! Dumps all values needed for restarting
! author Elizabeth C. Hunke, LANL

      subroutine write_restart_FY()

      use ice_fileunits, only: nu_dump_FY
      use ice_flux, only: frz_onset
      use ice_state, only: trcrn

      ! local variables

      logical (kind=log_kind) :: diag
      integer (kind=int_kind) :: nt_FY
      character(len=*),parameter :: subname='(write_restart_FY)'

      call icepack_query_tracer_indices(nt_FY_out=nt_FY)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      diag = .true.

      !-----------------------------------------------------------------

      call write_restart_field(nu_dump_FY,0,trcrn(:,:,nt_FY,:,:),'ruf8', &
                               'FY',ncat,diag)
      call write_restart_field(nu_dump_FY,0,frz_onset,'ruf8', &
                               'frz_onset',1,diag)

      end subroutine write_restart_FY

!=======================================================================

! Reads all values needed for an ice FY restart
! author Elizabeth C. Hunke, LANL

      subroutine read_restart_FY()

      use ice_fileunits, only: nu_restart_FY
      use ice_flux, only: frz_onset
      use ice_state, only: trcrn

      ! local variables

      logical (kind=log_kind) :: &
         diag
      integer (kind=int_kind) :: nt_FY
      character(len=*),parameter :: subname='(read_restart_FY)'

      call icepack_query_tracer_indices(nt_FY_out=nt_FY)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      diag = .true.

      if (my_task == master_task) write(nu_diag,*) subname,'min/max first-year ice area'

      call read_restart_field(nu_restart_FY,0,trcrn(:,:,nt_FY,:,:),'ruf8', &
                      'FY',ncat,diag,field_loc_center,field_type_scalar)

      if (my_task == master_task) write(nu_diag,*) subname,'min/max frz_onset'

      call read_restart_field(nu_restart_FY,0,frz_onset,'ruf8', &
                  'frz_onset',1,diag,field_loc_center,field_type_scalar)

      end subroutine read_restart_FY

!=======================================================================

! Dumps all values needed for restarting
!
! author Elizabeth C. Hunke, LANL

      subroutine write_restart_lvl()

      use ice_fileunits, only: nu_dump_lvl
      use ice_state, only: trcrn

      ! local variables

      logical (kind=log_kind) :: diag
      integer (kind=int_kind) :: nt_alvl, nt_vlvl
      character(len=*),parameter :: subname='(write_restart_lvl)'

      call icepack_query_tracer_indices(nt_alvl_out=nt_alvl, nt_vlvl_out=nt_vlvl)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      diag = .true.

      !-----------------------------------------------------------------

      call write_restart_field(nu_dump_lvl,0,trcrn(:,:,nt_alvl,:,:),'ruf8', &
                               'alvl',ncat,diag)
      call write_restart_field(nu_dump_lvl,0,trcrn(:,:,nt_vlvl,:,:),'ruf8', &
                               'vlvl',ncat,diag)

      end subroutine write_restart_lvl

!=======================================================================

! Reads all values needed for an ice lvl restart
!
! author Elizabeth C. Hunke, LANL

      subroutine read_restart_lvl()

      use ice_fileunits, only: nu_restart_lvl
      use ice_state, only: trcrn

      ! local variables

      logical (kind=log_kind) :: &
         diag
      integer (kind=int_kind) :: nt_alvl, nt_vlvl
      character(len=*),parameter :: subname='(read_restart_lvl)'

      call icepack_query_tracer_indices(nt_alvl_out=nt_alvl, nt_vlvl_out=nt_vlvl)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      diag = .true.

      if (my_task == master_task) write(nu_diag,*) subname,'min/max level ice area, volume'

      call read_restart_field(nu_restart_lvl,0,trcrn(:,:,nt_alvl,:,:),'ruf8', &
                       'alvl',ncat,diag,field_loc_center,field_type_scalar)
      call read_restart_field(nu_restart_lvl,0,trcrn(:,:,nt_vlvl,:,:),'ruf8', &
                       'vlvl',ncat,diag,field_loc_center,field_type_scalar)

      end subroutine read_restart_lvl

!=======================================================================
!
! Dumps all values needed for restarting
!
! authors Elizabeth C. Hunke, LANL
!         David A. Bailey, NCAR

      subroutine write_restart_pond_cesm()

      use ice_fileunits, only: nu_dump_pond
      use ice_state, only: trcrn

      ! local variables

      logical (kind=log_kind) :: diag
      integer (kind=int_kind) :: nt_apnd, nt_hpnd
      character(len=*),parameter :: subname='(write_restart_pond_cesm)'

      call icepack_query_tracer_indices(nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      diag = .true.

      call write_restart_field(nu_dump_pond,0,trcrn(:,:,nt_apnd,:,:),'ruf8', &
                               'apnd',ncat,diag)
      call write_restart_field(nu_dump_pond,0,trcrn(:,:,nt_hpnd,:,:),'ruf8', &
                               'hpnd',ncat,diag)

      end subroutine write_restart_pond_cesm

!=======================================================================

! Reads all values needed for a meltpond volume restart
!
! authors Elizabeth C. Hunke, LANL
!         David A. Bailey, NCAR

      subroutine read_restart_pond_cesm()

      use ice_fileunits, only: nu_restart_pond 
      use ice_state, only: trcrn

      ! local variables

      logical (kind=log_kind) :: &
         diag
      integer (kind=int_kind) :: nt_apnd, nt_hpnd
      character(len=*),parameter :: subname='(read_restart_pond_cesm)'

      call icepack_query_tracer_indices(nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      diag = .true.

      if (my_task == master_task) write(nu_diag,*) subname,'min/max cesm ponds'

      call read_restart_field(nu_restart_pond,0,trcrn(:,:,nt_apnd,:,:),'ruf8', &
                              'apnd',ncat,diag,field_loc_center,field_type_scalar)
      call read_restart_field(nu_restart_pond,0,trcrn(:,:,nt_hpnd,:,:),'ruf8', &
                              'hpnd',ncat,diag,field_loc_center,field_type_scalar)

      end subroutine read_restart_pond_cesm

!=======================================================================
!
! Dumps all values needed for restarting
!
! authors Elizabeth C. Hunke, LANL

      subroutine write_restart_pond_lvl()

      use ice_arrays_column, only: dhsn, ffracn
      use ice_fileunits, only: nu_dump_pond
      use ice_flux, only: fsnow
      use ice_state, only: trcrn

      ! local variables

      logical (kind=log_kind) :: diag
      integer (kind=int_kind) :: nt_apnd, nt_hpnd, nt_ipnd
      character(len=*),parameter :: subname='(write_restart_pond_lvl)'

      call icepack_query_tracer_indices(nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, &
           nt_ipnd_out=nt_ipnd)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      diag = .true.

      call write_restart_field(nu_dump_pond,0, trcrn(:,:,nt_apnd,:,:),'ruf8', &
                               'apnd',ncat,diag)
      call write_restart_field(nu_dump_pond,0, trcrn(:,:,nt_hpnd,:,:),'ruf8', &
                               'hpnd',ncat,diag)
      call write_restart_field(nu_dump_pond,0, trcrn(:,:,nt_ipnd,:,:),'ruf8', &
                               'ipnd',ncat,diag)
      call write_restart_field(nu_dump_pond,0, fsnow(:,:,          :),'ruf8', &
                               'fsnow',1,diag)
      call write_restart_field(nu_dump_pond,0,  dhsn(:,:,        :,:),'ruf8', &
                               'dhs',ncat,diag)
      call write_restart_field(nu_dump_pond,0,ffracn(:,:,        :,:),'ruf8', &
                               'ffrac',ncat,diag)

      end subroutine write_restart_pond_lvl

!=======================================================================

! Reads all values needed for a meltpond volume restart
!
! authors Elizabeth C. Hunke, LANL

      subroutine read_restart_pond_lvl()

      use ice_arrays_column, only: dhsn, ffracn
      use ice_fileunits, only: nu_restart_pond 
      use ice_flux, only: fsnow
      use ice_state, only: trcrn

      ! local variables

      logical (kind=log_kind) :: &
         diag
      integer (kind=int_kind) :: nt_apnd, nt_hpnd, nt_ipnd
      character(len=*),parameter :: subname='(read_restart_pond_lvl)'

      call icepack_query_tracer_indices(nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, &
           nt_ipnd_out=nt_ipnd)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      diag = .true.

      if (my_task == master_task) write(nu_diag,*) subname,'min/max level-ice ponds'

      call read_restart_field(nu_restart_pond,0, trcrn(:,:,nt_apnd,:,:),'ruf8', &
                              'apnd',ncat,diag,field_loc_center,field_type_scalar)
      call read_restart_field(nu_restart_pond,0, trcrn(:,:,nt_hpnd,:,:),'ruf8', &
                              'hpnd',ncat,diag,field_loc_center,field_type_scalar)
      call read_restart_field(nu_restart_pond,0, trcrn(:,:,nt_ipnd,:,:),'ruf8', &
                              'ipnd',ncat,diag,field_loc_center,field_type_scalar)
      call read_restart_field(nu_restart_pond,0, fsnow(:,:,          :),'ruf8', &
                              'fsnow',1,diag,field_loc_center,field_type_scalar)
      call read_restart_field(nu_restart_pond,0,  dhsn(:,:,        :,:),'ruf8', &
                              'dhs',ncat,diag,field_loc_center,field_type_scalar)
      call read_restart_field(nu_restart_pond,0,ffracn(:,:,        :,:),'ruf8', &
                              'ffrac',ncat,diag,field_loc_center,field_type_scalar)

      end subroutine read_restart_pond_lvl

!=======================================================================

! Dumps all values needed for restarting
!
! authors Elizabeth C. Hunke, LANL
!         David A. Bailey, NCAR

      subroutine write_restart_pond_topo()

      use ice_fileunits, only: nu_dump_pond
      use ice_state, only: trcrn

      ! local variables

      logical (kind=log_kind) :: diag
      integer (kind=int_kind) :: nt_apnd, nt_hpnd, nt_ipnd
      character(len=*),parameter :: subname='(write_restart_pond_topo)'

      call icepack_query_tracer_indices(nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, &
           nt_ipnd_out=nt_ipnd)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      diag = .true.

      call write_restart_field(nu_dump_pond,0,trcrn(:,:,nt_apnd,:,:),'ruf8', &
                               'apnd',ncat,diag)
      call write_restart_field(nu_dump_pond,0,trcrn(:,:,nt_hpnd,:,:),'ruf8', &
                               'hpnd',ncat,diag)
      call write_restart_field(nu_dump_pond,0,trcrn(:,:,nt_ipnd,:,:),'ruf8', &
                               'ipnd',ncat,diag)

      end subroutine write_restart_pond_topo

!=======================================================================

! Reads all values needed for a meltpond volume restart
!
! authors Elizabeth C. Hunke, LANL
!         David A. Bailey, NCAR

      subroutine read_restart_pond_topo()

      use ice_fileunits, only: nu_restart_pond 
      use ice_state, only: trcrn

      ! local variables

      logical (kind=log_kind) :: &
         diag
      integer (kind=int_kind) :: nt_apnd, nt_hpnd, nt_ipnd
      character(len=*),parameter :: subname='(read_restart_pond_topo)'

      call icepack_query_tracer_indices(nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, &
           nt_ipnd_out=nt_ipnd)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      diag = .true.

      if (my_task == master_task) write(nu_diag,*) subname,'min/max topo ponds'

      call read_restart_field(nu_restart_pond,0,trcrn(:,:,nt_apnd,:,:),'ruf8', &
                              'apnd',ncat,diag,field_loc_center,field_type_scalar)
      call read_restart_field(nu_restart_pond,0,trcrn(:,:,nt_hpnd,:,:),'ruf8', &
                              'hpnd',ncat,diag,field_loc_center,field_type_scalar)
      call read_restart_field(nu_restart_pond,0,trcrn(:,:,nt_ipnd,:,:),'ruf8', &
                              'ipnd',ncat,diag,field_loc_center,field_type_scalar)

      end subroutine read_restart_pond_topo

!=======================================================================

! Dumps all values needed for restarting
!
! authors Elizabeth Hunke, LANL (original version)
!         David Bailey, NCAR
!         Marika Holland, NCAR

      subroutine write_restart_aero()

      use ice_domain_size, only: n_aero
      use ice_state, only: trcrn
      use ice_fileunits, only: nu_dump_aero

      ! local variables

      integer (kind=int_kind) :: &
         k                    ! loop indices

      logical (kind=log_kind) :: diag

      character (len=3)       :: nchar
      integer (kind=int_kind) :: nt_aero
      character(len=*),parameter :: subname='(write_restart_aero)'

      !-----------------------------------------------------------------

      if (my_task == master_task) write(nu_diag,*) subname,'aerosols'

      call icepack_query_tracer_indices(nt_aero_out=nt_aero)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      diag = .true.

      do k = 1, n_aero
       write(nchar,'(i3.3)') k
       call write_restart_field(nu_dump_aero,0, &
            trcrn(:,:,nt_aero  +(k-1)*4,:,:),'ruf8','aerosnossl'//nchar, &
            ncat,diag)
       call write_restart_field(nu_dump_aero,0, &
            trcrn(:,:,nt_aero+1+(k-1)*4,:,:),'ruf8','aerosnoint'//nchar, &
            ncat,diag)
       call write_restart_field(nu_dump_aero,0, &
            trcrn(:,:,nt_aero+2+(k-1)*4,:,:),'ruf8','aeroicessl'//nchar, &
            ncat,diag)
       call write_restart_field(nu_dump_aero,0, &
            trcrn(:,:,nt_aero+3+(k-1)*4,:,:),'ruf8','aeroiceint'//nchar, &
            ncat,diag)
      enddo

      end subroutine write_restart_aero

!=======================================================================

! Reads all values needed for an ice aerosol restart
!
! authors Elizabeth Hunke, LANL (original version)
!         David Bailey, NCAR
!         Marika Holland, NCAR

      subroutine read_restart_aero()

      use ice_domain_size, only: n_aero
      use ice_state, only: trcrn
      use ice_fileunits, only: nu_restart_aero

      ! local variables

      integer (kind=int_kind) :: &
         k                    ! loop indices

      logical (kind=log_kind) :: &
         diag
      integer (kind=int_kind) :: nt_aero

      character (len=3)       :: nchar
      character(len=*),parameter :: subname='(read_restart_aero)'

      !-----------------------------------------------------------------

      if (my_task == master_task) write(nu_diag,*) subname,'aerosols'

      call icepack_query_tracer_indices(nt_aero_out=nt_aero)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      diag = .true.

      do k = 1, n_aero
       write(nchar,'(i3.3)') k
       call read_restart_field(nu_restart_aero,0, &
            trcrn(:,:,nt_aero  +(k-1)*4,:,:),'ruf8','aerosnossl'//trim(nchar), &
            ncat,diag,field_type=field_type_scalar,field_loc=field_loc_center)
       call read_restart_field(nu_restart_aero,0, &
            trcrn(:,:,nt_aero+1+(k-1)*4,:,:),'ruf8','aerosnoint'//trim(nchar), &
            ncat,diag,field_type=field_type_scalar,field_loc=field_loc_center)
       call read_restart_field(nu_restart_aero,0, &
            trcrn(:,:,nt_aero+2+(k-1)*4,:,:),'ruf8','aeroicessl'//trim(nchar), &
            ncat,diag,field_type=field_type_scalar,field_loc=field_loc_center)
       call read_restart_field(nu_restart_aero,0, &
            trcrn(:,:,nt_aero+3+(k-1)*4,:,:),'ruf8','aeroiceint'//trim(nchar), &
            ncat,diag,field_type=field_type_scalar,field_loc=field_loc_center)
      enddo

      end subroutine read_restart_aero

!=======================================================================

      subroutine read_restart_hbrine()

! Reads all values needed for hbrine
! author Elizabeth C. Hunke, LANL

      use ice_arrays_column, only: first_ice_real, first_ice
      use ice_blocks, only: block, get_block
      use ice_communicate, only: my_task, master_task
      use ice_domain, only: nblocks, blocks_ice
      use ice_fileunits, only: nu_restart_hbrine
      use ice_state, only: trcrn
      use ice_restart,only: read_restart_field

      ! local variables

      integer (kind=int_kind) :: &
         i, j, n, iblk    , & ! horizontal indices
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      type (block) :: &
         this_block      ! block information for current block

      logical (kind=log_kind) :: &
         diag

      integer (kind=int_kind) :: nt_fbri

      character(len=*),parameter :: subname='(read_restart_hbrine)'

      call icepack_query_tracer_indices(nt_fbri_out=nt_fbri)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      diag = .true.

      if (my_task == master_task) write(nu_diag,*) subname,'brine restart'

      call read_restart_field(nu_restart_hbrine,0,trcrn(:,:,nt_fbri,:,:),'ruf8', &
                              'fbrn',ncat,diag,field_loc_center,field_type_scalar)
      call read_restart_field(nu_restart_hbrine,0,first_ice_real(:,:,:,:),'ruf8', &
                              'first_ice',ncat,diag,field_loc_center,field_type_scalar)

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,n,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi  
            do n = 1, ncat
               if (first_ice_real(i,j,n,iblk) >= p5) then
                   first_ice     (i,j,n,iblk) = .true.
               else
                   first_ice     (i,j,n,iblk) = .false.
               endif
            enddo ! ncat
         enddo    ! i 
         enddo    ! j
      enddo       ! iblk

      end subroutine read_restart_hbrine

!=======================================================================

      subroutine write_restart_hbrine()

! Dumps all values needed for a hbrine restart
! author Elizabeth C. Hunke, LANL

      use ice_arrays_column, only: first_ice, first_ice_real
      use ice_blocks, only: block, get_block
      use ice_domain, only: nblocks, blocks_ice
      use ice_fileunits, only: nu_dump_hbrine
      use ice_state, only: trcrn
      use ice_restart,only: write_restart_field

      ! local variables

      integer (kind=int_kind) :: &
         i, j, n, iblk    , & ! horizontal indices
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      logical (kind=log_kind) :: diag

      integer (kind=int_kind) :: nt_fbri

      type (block) :: &
         this_block      ! block information for current block
      character(len=*),parameter :: subname='(write_restart_hbrine)'

      call icepack_query_tracer_indices(nt_fbri_out=nt_fbri)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      diag = .true.

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,n,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks

        this_block = get_block(blocks_ice(iblk),iblk)         
        ilo = this_block%ilo
        ihi = this_block%ihi
        jlo = this_block%jlo
        jhi = this_block%jhi

        do j = jlo, jhi
        do i = ilo, ihi  
           do n = 1, ncat
              if (first_ice     (i,j,n,iblk)) then
                  first_ice_real(i,j,n,iblk) = c1
              else
                  first_ice_real(i,j,n,iblk) = c0
              endif
           enddo ! n
        enddo    ! i
        enddo    ! j
      enddo      ! iblk

      call write_restart_field(nu_dump_hbrine,0,trcrn(:,:,nt_fbri,:,:),'ruf8', &
                               'fbrn',ncat,diag)
      call write_restart_field(nu_dump_hbrine,0,first_ice_real(:,:,:,:),'ruf8', &
                               'first_ice',ncat,diag)

      end subroutine write_restart_hbrine

!=======================================================================
!
! Dumps all values needed for a bgc restart
!
! author Elizabeth C. Hunke, LANL

      subroutine write_restart_bgc()

      use ice_arrays_column, only: Rayleigh_criteria, Rayleigh_real
      use ice_blocks, only: block, get_block
      use ice_domain, only: nblocks, blocks_ice
      use ice_domain_size, only: ncat, n_algae, n_doc, n_dic, &
          n_don, n_zaero, n_fed, n_fep
      use ice_fileunits, only: nu_dump_bgc
      use ice_flux_bgc, only: nit, amm, sil, dmsp, dms, algalN, &
          doc, don, dic, fed, fep, zaeros, hum
      use ice_state, only: trcrn
      use ice_flux, only: sss  
      use ice_restart, only:  write_restart_field

      ! local variables

      integer (kind=int_kind) :: &
       i, j, k, iblk    , & ! horizontal, vertical and block indices
       mm               , & ! n_algae
       ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      logical (kind=log_kind) :: diag

      character (len=3) :: nchar, ncharb

      integer (kind=int_kind) :: nt_bgc_S, nt_bgc_Am, &
         nt_bgc_DMS, nt_bgc_DMSPd, &
         nt_bgc_DMSPp, nt_bgc_Nit, nt_bgc_Sil, &
         nt_bgc_PON, nt_zbgc_frac, nt_bgc_hum, nbtrcr

      integer (kind=int_kind), dimension(icepack_max_algae) :: &  
         nt_bgc_N , & ! diatoms, phaeocystis, pico/small   
         nt_bgc_C , & ! diatoms, phaeocystis, pico/small   
         nt_bgc_chl   ! diatoms, phaeocystis, pico/small 

      integer (kind=int_kind), dimension(icepack_max_doc) :: &  
         nt_bgc_DOC      !  dissolved organic carbon

      integer (kind=int_kind), dimension(icepack_max_don) :: & 
         nt_bgc_DON         !  dissolved organic nitrogen

      integer (kind=int_kind), dimension(icepack_max_dic) :: &  
         nt_bgc_DIC         !  dissolved inorganic carbon

      integer (kind=int_kind), dimension(icepack_max_fe) :: & 
         nt_bgc_Fed,     & !  dissolved iron
         nt_bgc_Fep        !  particulate iron

      integer (kind=int_kind), dimension(icepack_max_aero) :: &  
         nt_zaero       !  black carbon and other aerosols
      
      logical (kind=log_kind) :: tr_bgc_Nit, tr_bgc_Am, tr_bgc_Sil,&
         tr_bgc_DMS, tr_bgc_PON, tr_bgc_N, tr_bgc_C, &
         tr_bgc_DON, tr_bgc_Fe,  tr_zaero , tr_bgc_chl, &
         tr_bgc_hum

      logical (kind=log_kind) :: skl_bgc, solve_zsal

      type (block) :: &
         this_block      ! block information for current block

      integer (kind=int_kind) :: &
         ipoint

      character(len=*),parameter :: subname='(write_restart_bgc)'

      call icepack_query_parameters(skl_bgc_out=skl_bgc, solve_zsal_out=solve_zsal)
      call icepack_query_tracer_numbers(nbtrcr_out=nbtrcr)
      call icepack_query_tracer_flags(tr_bgc_Nit_out=tr_bgc_Nit, &
          tr_bgc_Am_out=tr_bgc_Am, tr_bgc_Sil_out=tr_bgc_Sil, &
          tr_bgc_DMS_out=tr_bgc_DMS, tr_bgc_PON_out=tr_bgc_PON, &
          tr_bgc_N_out=tr_bgc_N, tr_bgc_C_out=tr_bgc_C, &
          tr_bgc_DON_out=tr_bgc_DON, tr_bgc_Fe_out=tr_bgc_Fe,  tr_zaero_out=tr_zaero, &
          tr_bgc_chl_out=tr_bgc_chl, tr_bgc_hum_out=tr_bgc_hum)
      call icepack_query_tracer_indices(nt_bgc_S_out=nt_bgc_S, nt_bgc_Am_out=nt_bgc_Am, &
          nt_bgc_DMS_out=nt_bgc_DMS, nt_bgc_DMSPd_out=nt_bgc_DMSPd, &
          nt_bgc_C_out=nt_bgc_C, nt_bgc_chl_out=nt_bgc_chl, &
          nt_bgc_DMSPp_out=nt_bgc_DMSPp, nt_bgc_Nit_out=nt_bgc_Nit, &
          nt_bgc_Sil_out=nt_bgc_Sil, nt_bgc_PON_out=nt_bgc_PON, &
          nt_bgc_DON_out=nt_bgc_DON, nt_bgc_DOC_out=nt_bgc_DOC, nt_bgc_DIC_out=nt_bgc_DIC, &
          nt_bgc_N_out=nt_bgc_N, nt_zaero_out=nt_zaero, nt_bgc_Fed_out=nt_bgc_Fed, &
          nt_bgc_hum_out=nt_bgc_hum, nt_bgc_Fep_out=nt_bgc_Fep, nt_zbgc_frac_out=nt_zbgc_frac)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      diag = .true.

      !-----------------------------------------------------------------
      ! Salinity and extras
      !-----------------------------------------------------------------
      if (solve_zsal) then 

      do k = 1,nblyr
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_S+k-1,:,:),'ruf8', &
                   'zSalinity'//trim(nchar),ncat,diag)
      enddo
    
      call write_restart_field(nu_dump_bgc,0,sss,'ruf8','sss',1,diag)

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi
         do j = jlo, jhi
         do i = ilo, ihi
            if (Rayleigh_criteria(i,j,iblk)) then
                Rayleigh_real    (i,j,iblk) = c1
            elseif (.NOT. Rayleigh_criteria(i,j,iblk)) then
                Rayleigh_real    (i,j,iblk) = c0
            endif
         enddo
         enddo
      enddo

      call write_restart_field(nu_dump_bgc,0,Rayleigh_real,'ruf8','Rayleigh',1,diag)

      endif ! solve_zsal

      !-----------------------------------------------------------------
      ! Skeletal layer BGC
      !-----------------------------------------------------------------

      if (skl_bgc) then
         do k = 1, n_algae
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_N(k),:,:), &
                                  'ruf8','bgc_N'//trim(nchar),ncat,diag)
            if (tr_bgc_chl) &
            call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_chl(k),:,:), &
                                  'ruf8','bgc_chl'//trim(nchar),ncat,diag)
         enddo
        if (tr_bgc_C)  then
          ! do k = 1, n_algae
          !  write(nchar,'(i3.3)') k
          !  call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_C(k),:,:), &
          !                        'ruf8','bgc_C'//trim(nchar),ncat,diag)
          ! enddo
           do k = 1, n_doc
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_DOC(k),:,:), &
                                  'ruf8','bgc_DOC'//trim(nchar),ncat,diag)
           enddo
           do k = 1, n_dic
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_DIC(k),:,:), &
                                  'ruf8','bgc_DIC'//trim(nchar),ncat,diag)
           enddo
         endif
         if (tr_bgc_Nit) &
         call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_Nit,:,:), &
                                  'ruf8','bgc_Nit',ncat,diag)
         if (tr_bgc_Am) &
         call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_Am,:,:), &
                                  'ruf8','bgc_Am',ncat,diag)
         if (tr_bgc_Sil) &
         call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_Sil,:,:), &
                                  'ruf8','bgc_Sil',ncat,diag)
         if (tr_bgc_hum) &
         call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_hum,:,:), &
                                  'ruf8','bgc_hum',ncat,diag)
         if (tr_bgc_DMS) then
           call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_DMSPp,:,:), &
                                  'ruf8','bgc_DMSPp',ncat,diag)
           call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_DMSPd,:,:), &
                                  'ruf8','bgc_DMSPd',ncat,diag)
           call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_DMS,:,:), &
                                  'ruf8','bgc_DMS',ncat,diag)
         endif
         if (tr_bgc_PON) &
         call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_PON,:,:), &
                                  'ruf8','bgc_PON',ncat,diag)
      
        if (tr_bgc_DON)  then
           do k = 1, n_don
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_DON(k),:,:), &
                                  'ruf8','bgc_DON'//trim(nchar),ncat,diag)
           enddo
         endif
        if (tr_bgc_Fe )  then
           do k = 1, n_fed 
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_Fed (k),:,:), &
                                  'ruf8','bgc_Fed'//trim(nchar),ncat,diag)
           enddo
           do k = 1, n_fep 
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_Fep (k),:,:), &
                                  'ruf8','bgc_Fep'//trim(nchar),ncat,diag)
           enddo
         endif

      else 

      !-----------------------------------------------------------------
      ! Z layer BGC
      !-----------------------------------------------------------------

         if (tr_bgc_Nit) then
         do k = 1,nblyr+3
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_Nit+k-1,:,:),'ruf8', &
                                 'bgc_Nit'//trim(nchar),ncat,diag)
         enddo
         endif
         if (tr_bgc_N) then
         do mm = 1,n_algae
         write(ncharb, '(i3.3)') mm
         do k = 1,nblyr+3
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0, &
                                  trcrn(:,:,nt_bgc_N(mm)+k-1,:,:),'ruf8', &
                                 'bgc_N'//trim(ncharb)//trim(nchar),ncat,diag)
         enddo
         if (tr_bgc_chl) then
         do k = 1,nblyr+3
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,  &
                                  trcrn(:,:,nt_bgc_chl(mm)+k-1,:,:),'ruf8', &
                                 'bgc_chl'//trim(ncharb)//trim(nchar),ncat,diag)
         enddo
         endif
         enddo   !n_algae
         endif   ! tr_bgc_N
         if (tr_bgc_C) then
        ! do mm = 1,n_algae
        ! write(ncharb, '(i3.3)') mm
        ! do k = 1,nblyr+3
        !    write(nchar,'(i3.3)') k
        !    call write_restart_field(nu_dump_bgc,0,  &
        !                          trcrn(:,:,nt_bgc_C(mm)+k-1,:,:),'ruf8', &
        !                         'bgc_C'//trim(ncharb)//trim(nchar),ncat,diag)
        ! enddo
        ! enddo
         do mm = 1,n_doc
         write(ncharb, '(i3.3)') mm
         do k = 1,nblyr+3
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,  &
                                  trcrn(:,:,nt_bgc_DOC(mm)+k-1,:,:),'ruf8', &
                                 'bgc_DOC'//trim(ncharb)//trim(nchar),ncat,diag)
         enddo
         enddo
         do mm = 1,n_dic
         write(ncharb, '(i3.3)') mm
         do k = 1,nblyr+3
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,  &
                                  trcrn(:,:,nt_bgc_DIC(mm)+k-1,:,:),'ruf8', &
                                 'bgc_DIC'//trim(ncharb)//trim(nchar),ncat,diag)
         enddo
         enddo
         endif  !tr_bgc_C
         if (tr_bgc_Am) then
         do k = 1,nblyr+3
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_Am+k-1,:,:),'ruf8', &
                                 'bgc_Am'//trim(nchar),ncat,diag)
         enddo
         endif
         if (tr_bgc_Sil) then
         do k = 1,nblyr+3
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_Sil+k-1,:,:),'ruf8', &
                                 'bgc_Sil'//trim(nchar),ncat,diag)
         enddo
         endif
         if (tr_bgc_hum) then
         do k = 1,nblyr+3
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_hum+k-1,:,:),'ruf8', &
                                 'bgc_hum'//trim(nchar),ncat,diag)
         enddo
         endif
         if (tr_bgc_DMS) then
         do k = 1,nblyr+3
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_DMSPp+k-1,:,:),'ruf8', &
                                 'bgc_DMSPp'//trim(nchar),ncat,diag)
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_DMSPd+k-1,:,:),'ruf8', &
                                 'bgc_DMSPd'//trim(nchar),ncat,diag)
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_DMS+k-1,:,:),'ruf8', &
                                 'bgc_DMS'//trim(nchar),ncat,diag)
         enddo
         endif
         if (tr_bgc_PON) then
         do k = 1,nblyr+3
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_PON+k-1,:,:),'ruf8', &
                                 'bgc_PON'//trim(nchar),ncat,diag)
         enddo
         endif
         if (tr_bgc_DON) then
         do mm = 1,n_don
         write(ncharb, '(i3.3)') mm
         do k = 1,nblyr+3
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,  &
                                  trcrn(:,:,nt_bgc_DON(mm)+k-1,:,:),'ruf8', &
                                 'bgc_DON'//trim(ncharb)//trim(nchar),ncat,diag)
         enddo
         enddo
         endif
         if (tr_bgc_Fe ) then
         do mm = 1,n_fed
         write(ncharb, '(i3.3)') mm
         do k = 1,nblyr+3
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,  &
                                  trcrn(:,:,nt_bgc_Fed(mm)+k-1,:,:),'ruf8', &
                                 'bgc_Fed'//trim(ncharb)//trim(nchar),ncat,diag)
         enddo
         enddo
         do mm = 1,n_fep
         write(ncharb, '(i3.3)') mm
         do k = 1,nblyr+3
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,  &
                                  trcrn(:,:,nt_bgc_Fep(mm)+k-1,:,:),'ruf8', &
                                 'bgc_Fep'//trim(ncharb)//trim(nchar),ncat,diag)
         enddo
         enddo
         endif
         if (tr_zaero) then
         do mm = 1,n_zaero
         write(ncharb, '(i3.3)') mm
         do k = 1,nblyr+3
            write(nchar,'(i3.3)') k
            call write_restart_field(nu_dump_bgc,0,  &
                                  trcrn(:,:,nt_zaero(mm)+k-1,:,:),'ruf8', &
                                 'zaero'//trim(ncharb)//trim(nchar),ncat,diag)
         enddo
         enddo
         endif
         do mm = 1,nbtrcr
          write(nchar,'(i3.3)') mm
          call write_restart_field(nu_dump_bgc,0,   &
                                trcrn(:,:,nt_zbgc_frac+mm-1,:,:),'ruf8', &
                                'zbgc_frac'//trim(nchar),ncat,diag)
         enddo
      endif

      !-----------------------------------------------------------------
      ! Ocean BGC
      !-----------------------------------------------------------------

      if (tr_bgc_N) then
      do k = 1,n_algae
          write(nchar,'(i3.3)') k
          call write_restart_field(nu_dump_bgc,0,algalN(:,:,k,:),'ruf8','algalN'//trim(nchar),1,diag)
      enddo  !k
      endif
      if (tr_bgc_C) then
      do k = 1,n_doc
          write(nchar,'(i3.3)') k
          call write_restart_field(nu_dump_bgc,0,doc(:,:,k,:),'ruf8','doc'//trim(nchar),1,diag)
      enddo  !k
      do k = 1,n_dic
          write(nchar,'(i3.3)') k
          call write_restart_field(nu_dump_bgc,0,dic(:,:,k,:),'ruf8','dic'//trim(nchar),1,diag)
      enddo  !k
      endif      
      if (tr_bgc_Nit) &
      call write_restart_field(nu_dump_bgc,0,nit,   'ruf8','nit',   1,diag)
      if (tr_bgc_Am) &
      call write_restart_field(nu_dump_bgc,0,amm,   'ruf8','amm',   1,diag)
      if (tr_bgc_Sil) &
      call write_restart_field(nu_dump_bgc,0,sil,   'ruf8','sil',   1,diag)
      if (tr_bgc_hum) &
      call write_restart_field(nu_dump_bgc,0,hum,   'ruf8','hum',   1,diag)
      if (tr_bgc_DMS) then
        call write_restart_field(nu_dump_bgc,0,dmsp,  'ruf8','dmsp',  1,diag)
        call write_restart_field(nu_dump_bgc,0,dms,   'ruf8','dms',   1,diag)
      endif
      if (tr_bgc_DON) then
      do k = 1,n_don
          write(nchar,'(i3.3)') k
          call write_restart_field(nu_dump_bgc,0,don(:,:,k,:),'ruf8','don'//trim(nchar),1,diag)
      enddo  !k
      endif
      if (tr_bgc_Fe ) then
      do k = 1,n_fed
          write(nchar,'(i3.3)') k
          call write_restart_field(nu_dump_bgc,0,fed(:,:,k,:),'ruf8','fed'//trim(nchar),1,diag)
      enddo  !k
      do k = 1,n_fep
          write(nchar,'(i3.3)') k
          call write_restart_field(nu_dump_bgc,0,fep(:,:,k,:),'ruf8','fep'//trim(nchar),1,diag)
      enddo  !k
      endif
      if (tr_zaero) then
      do k = 1,n_zaero
          write(nchar,'(i3.3)') k
          call write_restart_field(nu_dump_bgc,0,zaeros(:,:,k,:),'ruf8','zaeros'//trim(nchar),1,diag)
      enddo  !k
      endif

      end subroutine write_restart_bgc

!=======================================================================
!
! Reads all values needed for a bgc restart
!
! author Elizabeth C. Hunke, LANL

      subroutine read_restart_bgc()

      use ice_arrays_column, only: Rayleigh_real, Rayleigh_criteria
      use ice_blocks, only: block, get_block
      use ice_communicate, only: my_task, master_task
      use ice_domain, only: nblocks, blocks_ice
      use ice_domain_size, only: ncat, n_algae, n_doc, n_dic,&
          n_don, n_zaero, n_fed, n_fep
      use ice_fileunits, only: nu_restart_bgc
      use ice_flux, only: sss  
      use ice_flux_bgc, only: nit, amm, sil, dmsp, dms, algalN, &
          doc, don, dic, fed, fep, zaeros, hum
      use ice_state, only: trcrn
      use ice_restart, only: read_restart_field

      ! local variables

      type (block) :: &
         this_block      ! block information for current block

      integer (kind=int_kind) :: &
         i, j, k, iblk, & ! indices
         mm           , & ! n_algae
         ilo,ihi,jlo,jhi  ! beginning and end of physical domain

      logical (kind=log_kind) :: diag

      integer (kind=int_kind) :: nt_bgc_S, nt_bgc_Am, &
         nt_bgc_DMS, nt_bgc_DMSPd, &
         nt_bgc_DMSPp, nt_bgc_Nit, nt_bgc_Sil, &
         nt_bgc_PON, nt_zbgc_frac, nt_bgc_hum, nbtrcr

      integer (kind=int_kind), dimension(icepack_max_algae) :: &  
         nt_bgc_N , & ! diatoms, phaeocystis, pico/small   
         nt_bgc_C , & ! diatoms, phaeocystis, pico/small   
         nt_bgc_chl   ! diatoms, phaeocystis, pico/small 

      integer (kind=int_kind), dimension(icepack_max_doc) :: &  
         nt_bgc_DOC      !  dissolved organic carbon

      integer (kind=int_kind), dimension(icepack_max_don) :: & 
         nt_bgc_DON         !  dissolved organic nitrogen

      integer (kind=int_kind), dimension(icepack_max_dic) :: &  
         nt_bgc_DIC         !  dissolved inorganic carbon

      integer (kind=int_kind), dimension(icepack_max_fe) :: & 
         nt_bgc_Fed,     & !  dissolved iron
         nt_bgc_Fep        !  particulate iron

      integer (kind=int_kind), dimension(icepack_max_aero) :: &  
         nt_zaero       !  black carbon and other aerosols
      
      logical (kind=log_kind) :: tr_bgc_Nit, tr_bgc_Am, tr_bgc_Sil,&
         tr_bgc_DMS, tr_bgc_PON, tr_bgc_N, tr_bgc_C, &
         tr_bgc_DON, tr_bgc_Fe,  tr_zaero , tr_bgc_chl, &
         tr_bgc_hum

      logical (kind=log_kind) :: skl_bgc, solve_zsal

      character (len=3) :: nchar, ncharb

      character(len=*),parameter :: subname='(read_restart_bgc)'

      call icepack_query_parameters(skl_bgc_out=skl_bgc, solve_zsal_out=solve_zsal)
      call icepack_query_tracer_numbers(nbtrcr_out=nbtrcr)
      call icepack_query_tracer_flags(tr_bgc_Nit_out=tr_bgc_Nit, &
          tr_bgc_Am_out=tr_bgc_Am, tr_bgc_Sil_out=tr_bgc_Sil, &
          tr_bgc_DMS_out=tr_bgc_DMS, tr_bgc_PON_out=tr_bgc_PON, &
          tr_bgc_N_out=tr_bgc_N, tr_bgc_C_out=tr_bgc_C, &
          tr_bgc_DON_out=tr_bgc_DON, tr_bgc_Fe_out=tr_bgc_Fe,  tr_zaero_out=tr_zaero, &
          tr_bgc_chl_out=tr_bgc_chl, tr_bgc_hum_out=tr_bgc_hum)
      call icepack_query_tracer_indices(nt_bgc_S_out=nt_bgc_S, nt_bgc_Am_out=nt_bgc_Am, &
          nt_bgc_DMS_out=nt_bgc_DMS, nt_bgc_DMSPd_out=nt_bgc_DMSPd, &
          nt_bgc_C_out=nt_bgc_C, nt_bgc_chl_out=nt_bgc_chl, &
          nt_bgc_DMSPp_out=nt_bgc_DMSPp, nt_bgc_Nit_out=nt_bgc_Nit, &
          nt_bgc_Sil_out=nt_bgc_Sil, nt_bgc_PON_out=nt_bgc_PON, &
          nt_bgc_DON_out=nt_bgc_DON, nt_bgc_DOC_out=nt_bgc_DOC, nt_bgc_DIC_out=nt_bgc_DIC, &
          nt_bgc_N_out=nt_bgc_N, nt_zaero_out=nt_zaero, nt_bgc_Fed_out=nt_bgc_Fed, &
          nt_bgc_Fep_out=nt_bgc_Fep, nt_zbgc_frac_out=nt_zbgc_frac, nt_bgc_hum_out=nt_bgc_hum)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      diag = .true.

      !-----------------------------------------------------------------
      ! Salinity and extras
      !-----------------------------------------------------------------

      if (restart_zsal) then 

      if (my_task == master_task) write(nu_diag,*) subname,'zSalinity restart'
      do k = 1,nblyr
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_S+k-1,:,:),'ruf8', &
              'zSalinity'//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
 
      if (my_task == master_task) write(nu_diag,*) subname,'sea surface salinity'
      call read_restart_field(nu_restart_bgc,0,sss,'ruf8','sss',1,diag)
      call read_restart_field(nu_restart_bgc,0,Rayleigh_real,'ruf8','Rayleigh',1,diag)

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j = jlo, jhi
         do i = ilo, ihi  
            if (Rayleigh_real     (i,j,iblk) .GE. c1) then
                Rayleigh_criteria (i,j,iblk) = .true.
            elseif (Rayleigh_real (i,j,iblk) < c1) then
                Rayleigh_criteria (i,j,iblk) = .false.
            endif
         enddo
         enddo
      enddo
      endif ! restart_zsal

      !-----------------------------------------------------------------
      ! Skeletal Layer BGC
      !-----------------------------------------------------------------
      if (restart_bgc) then

      if (skl_bgc) then
       if (my_task == master_task) write(nu_diag,*) subname,'skl bgc restart'

       do k = 1, n_algae
          write(nchar,'(i3.3)') k
          call read_restart_field(nu_restart_bgc,0, &
                 trcrn(:,:,nt_bgc_N(k),:,:), &
                 'ruf8','bgc_N'//trim(nchar),ncat,diag)
          if (tr_bgc_chl) &
          call read_restart_field(nu_restart_bgc,0,  &
                 trcrn(:,:,nt_bgc_chl(k),:,:), &
                 'ruf8','bgc_chl'//trim(nchar),ncat,diag)
       enddo   !k
       if (tr_bgc_C) then
         ! do k = 1, n_algae
         !     write(nchar,'(i3.3)') k
         !     call read_restart_field(nu_restart_bgc,0,  &
         !        trcrn(:,:,nt_bgc_C(k),:,:), &
         !        'ruf8','bgc_C'//trim(nchar),ncat,diag)
         ! enddo
          do k = 1, n_doc
              write(nchar,'(i3.3)') k
              call read_restart_field(nu_restart_bgc,0,  &
                 trcrn(:,:,nt_bgc_DOC(k),:,:), &
                 'ruf8','bgc_DOC'//trim(nchar),ncat,diag)
          enddo
          do k = 1, n_dic
              write(nchar,'(i3.3)') k
              call read_restart_field(nu_restart_bgc,0,  &
                 trcrn(:,:,nt_bgc_DIC(k),:,:), &
                 'ruf8','bgc_DIC'//trim(nchar),ncat,diag)
          enddo
       endif
       if (tr_bgc_Nit) &
       call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_Nit,:,:), &
           'ruf8','bgc_Nit',ncat,diag)
       if (tr_bgc_Am) &
       call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_Am,:,:), &
           'ruf8','bgc_Am',ncat,diag)
       if (tr_bgc_Sil) &
       call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_Sil,:,:), &
           'ruf8','bgc_Sil',ncat,diag)
       if (tr_bgc_hum) &
       call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_hum,:,:), &
           'ruf8','bgc_hum',ncat,diag)
       if(tr_bgc_DMS) then
         call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_DMSPp,:,:), &
           'ruf8','bgc_DMSPp',ncat,diag)
         call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_DMSPd,:,:), &
           'ruf8','bgc_DMSPd',ncat,diag)
         call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_DMS,:,:), &
           'ruf8','bgc_DMS',ncat,diag)
       endif
       if (tr_bgc_PON) &
       call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_PON,:,:), &
           'ruf8','bgc_PON',ncat,diag)
       if (tr_bgc_DON) then
          do k = 1, n_don
              write(nchar,'(i3.3)') k
              call read_restart_field(nu_restart_bgc,0,  &
                 trcrn(:,:,nt_bgc_DON(k),:,:), &
                 'ruf8','bgc_DON'//trim(nchar),ncat,diag)
          enddo
       endif
       if (tr_bgc_Fe) then
          do k = 1, n_fed 
              write(nchar,'(i3.3)') k
              call read_restart_field(nu_restart_bgc,0,  &
                 trcrn(:,:,nt_bgc_Fed (k),:,:), &
                 'ruf8','bgc_Fed'//trim(nchar),ncat,diag)
          enddo
          do k = 1, n_fep 
              write(nchar,'(i3.3)') k
              call read_restart_field(nu_restart_bgc,0,  &
                 trcrn(:,:,nt_bgc_Fep (k),:,:), &
                 'ruf8','bgc_Fep'//trim(nchar),ncat,diag)
          enddo
       endif

      else

      !-----------------------------------------------------------------
      ! Z Layer BGC
      !-----------------------------------------------------------------

      if (tr_bgc_Nit) then
      if (my_task == master_task) write(nu_diag,*) subname,'z bgc restart: min/max Nitrate'
      do k=1,nblyr+3
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_Nit+k-1,:,:),'ruf8', &
              'bgc_Nit'//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      endif   !Nit
      if (tr_bgc_N) then
      do mm = 1,n_algae
      write(ncharb,'(i3.3)') mm
      if (my_task == master_task) write(nu_diag,*) subname,' min/max Algal N'
      do k=1,nblyr+3
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0, &
               trcrn(:,:,nt_bgc_N(mm)+k-1,:,:),'ruf8', &
              'bgc_N'//trim(ncharb)//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      if (tr_bgc_chl) then
      if (my_task == master_task) write(nu_diag,*) subname,' min/max Algal chla'
      do k=1,nblyr+3
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,  &
               trcrn(:,:,nt_bgc_chl(mm)+k-1,:,:),'ruf8', &
              'bgc_chl'//trim(ncharb)//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      endif   ! tr_bgc_chl
      enddo   !n_algae
      endif   ! tr_bgc_N
      if (tr_bgc_C) then
     ! do mm = 1,n_algae
     ! write(ncharb,'(i3.3)') mm
     ! if (my_task == master_task) write(nu_diag,*) subname,' min/max Algal C'
     ! do k=1,nblyr+3
     !    write(nchar,'(i3.3)') k
     !    call read_restart_field(nu_restart_bgc,0,  &
     !          trcrn(:,:,nt_bgc_C(mm)+k-1,:,:),'ruf8', &
     !         'bgc_C'//trim(ncharb)//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
     ! enddo
     ! enddo  !mm
      do mm = 1,n_doc
      write(ncharb,'(i3.3)') mm
      if (my_task == master_task) write(nu_diag,*) subname,' min/max DOC'
      do k=1,nblyr+3
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,  &
               trcrn(:,:,nt_bgc_DOC(mm)+k-1,:,:),'ruf8', &
              'bgc_DOC'//trim(ncharb)//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      enddo  !mm
      do mm = 1,n_dic
      write(ncharb,'(i3.3)') mm
      if (my_task == master_task) write(nu_diag,*) subname,' min/max DIC'
      do k=1,nblyr+3
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,  &
               trcrn(:,:,nt_bgc_DIC(mm)+k-1,:,:),'ruf8', &
              'bgc_DIC'//trim(ncharb)//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      enddo  !mm
      endif  ! tr_bgc_C
      if (tr_bgc_Am) then
      if (my_task == master_task) write(nu_diag,*) subname,' min/max ammonium'
      do k=1,nblyr+3
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_Am+k-1,:,:),'ruf8', &
              'bgc_Am'//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      endif
      if (tr_bgc_Sil) then
      if (my_task == master_task) write(nu_diag,*) subname,' min/max silicate'
      do k=1,nblyr+3
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_Sil+k-1,:,:),'ruf8', &
              'bgc_Sil'//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      endif
      if (tr_bgc_hum) then
      if (my_task == master_task) write(nu_diag,*) subname,' min/max humic material'
      do k=1,nblyr+3
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_hum+k-1,:,:),'ruf8', &
              'bgc_hum'//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      endif
      if (tr_bgc_DMS) then
      do k=1,nblyr+3
      if (my_task == master_task) write(nu_diag,*) subname,' min/max DMSPp'
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_DMSPp+k-1,:,:),'ruf8', &
              'bgc_DMSPp'//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      if (my_task == master_task) write(nu_diag,*) subname,' min/max DMSPd'
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_DMSPd+k-1,:,:),'ruf8', &
              'bgc_DMSPd'//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      if (my_task == master_task) write(nu_diag,*) subname,' min/max DMS'
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_DMS+k-1,:,:),'ruf8', &
              'bgc_DMS'//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      endif
      if (tr_bgc_PON) then
      if (my_task == master_task) write(nu_diag,*) subname,' min/max PON'
      do k=1,nblyr+3
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_PON+k-1,:,:),'ruf8', &
              'bgc_PON'//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      endif
      if (tr_bgc_DON) then
      do mm = 1,n_don
      write(ncharb,'(i3.3)') mm
      if (my_task == master_task) write(nu_diag,*) subname,' min/max DON'
      do k=1,nblyr+3
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,  &
               trcrn(:,:,nt_bgc_DON(mm)+k-1,:,:),'ruf8', &
              'bgc_DON'//trim(ncharb)//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      enddo  !mm
      endif
      if (tr_bgc_Fe) then
      do mm = 1,n_fed
      write(ncharb,'(i3.3)') mm
      if (my_task == master_task) write(nu_diag,*) subname,' min/max dFe '
      do k=1,nblyr+3
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,  &
               trcrn(:,:,nt_bgc_Fed (mm)+k-1,:,:),'ruf8', &
              'bgc_Fed'//trim(ncharb)//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      enddo  !mm
      do mm = 1,n_fep
      write(ncharb,'(i3.3)') mm
      if (my_task == master_task) write(nu_diag,*) subname,' min/max pFe '
      do k=1,nblyr+3
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,  &
               trcrn(:,:,nt_bgc_Fep (mm)+k-1,:,:),'ruf8', &
              'bgc_Fep'//trim(ncharb)//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      enddo  !mm
      endif
      if (tr_zaero) then
      do mm = 1,n_zaero
      write(ncharb,'(i3.3)') mm
      if (my_task == master_task) write(nu_diag,*) subname,' min/max z aerosols'
      do k=1,nblyr+3
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,  &
               trcrn(:,:,nt_zaero(mm)+k-1,:,:),'ruf8', &
              'zaero'//trim(ncharb)//trim(nchar),ncat,diag,field_loc_center,field_type_scalar)
      enddo
      enddo  !mm
      endif
      do mm = 1,nbtrcr
          write(nchar,'(i3.3)') mm
          call read_restart_field(nu_restart_bgc,0, &
               trcrn(:,:,nt_zbgc_frac+mm-1,:,:),'ruf8', &
               'zbgc_frac'//trim(nchar),ncat,diag)
      enddo
      endif

      !-----------------------------------------------------------------
      ! Ocean BGC
      !-----------------------------------------------------------------

      if (my_task == master_task) write(nu_diag,*) subname,'mixed layer ocean bgc restart'
      if (tr_bgc_N) then
      do k = 1,n_algae
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,algalN(:,:,k,:),'ruf8','algalN'//trim(nchar),1,diag)
      enddo  !k
      endif
      if (tr_bgc_C) then
      do k = 1,n_doc
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,doc(:,:,k,:),'ruf8','doc'//trim(nchar),1,diag)
      enddo  !k
      do k = 1,n_dic
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,dic(:,:,k,:),'ruf8','dic'//trim(nchar),1,diag)
      enddo  !k
      endif  !tr_bgc_C

      if (tr_bgc_Nit) &
      call read_restart_field(nu_restart_bgc,0,nit   ,'ruf8','nit'   ,&
                              1,diag,field_loc_center,field_type_scalar)
      if (tr_bgc_Am) &
      call read_restart_field(nu_restart_bgc,0,amm   ,'ruf8','amm'   ,&
                              1,diag,field_loc_center,field_type_scalar)
      if (tr_bgc_Sil) &
      call read_restart_field(nu_restart_bgc,0,sil   ,'ruf8','sil'   ,&
                              1,diag,field_loc_center,field_type_scalar)
      if (tr_bgc_hum) &
      call read_restart_field(nu_restart_bgc,0,hum   ,'ruf8','hum'   ,&
                              1,diag,field_loc_center,field_type_scalar)
      if (tr_bgc_DMS) then
        call read_restart_field(nu_restart_bgc,0,dmsp  ,'ruf8','dmsp'  ,1,diag)
        call read_restart_field(nu_restart_bgc,0,dms   ,'ruf8','dms'   ,1,diag)
      endif
      if (tr_bgc_DON) then
      do k = 1,n_don
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,don(:,:,k,:),'ruf8','don'//trim(nchar),1,diag)
      enddo  !k
      endif
      if (tr_bgc_Fe ) then
      do k = 1,n_fed
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,fed(:,:,k,:),'ruf8','fed'//trim(nchar),1,diag)
      enddo  !k
      do k = 1,n_fep
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,fep(:,:,k,:),'ruf8','fep'//trim(nchar),1,diag)
      enddo  !k
      endif
      if (tr_zaero) then
      do k = 1,n_zaero
         write(nchar,'(i3.3)') k
         call read_restart_field(nu_restart_bgc,0,zaeros(:,:,k,:),'ruf8','zaeros'//trim(nchar),1,diag)
      enddo  !k
      endif
      endif  ! restart_bgc
   
      end subroutine read_restart_bgc

!=======================================================================

      end module ice_restart_column

!=======================================================================
