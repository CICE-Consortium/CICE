!=======================================================================
!
! Reads and applies restoring data
! This handles both the halo restoring and the interior restoring
!   as two separate features.
! The halo restoring is on if set_boundary_flds is not set to 'none' or ''.
! The interior restoring is on if restore_ice is true and restore_flds is
!   not set to 'none' or ''.
! The restoring data for both halos and interior is one set of arrays,
! and the point by point restoring on the halo and interior is controlled
! by the mask_restoring and fval_restoring arrays.  The halo restoring
! is implemented so it's always set to the restoring data exactly.  The
! interior restoring is specified by the restore_mask, restore_width, and 
! restore_timescale namelist inputs.
!
! A few basic options are supported for reading restoring data through the
! restore_data namelist inputa and the subroutine ice_restoring_getdata.
! Users are encouraged to implement custom data reading to handle various
! file formats, fields, and possible spatial or temporal interpolation.
!
! authors: T. Craig

   module ice_restoring

      use ice_kinds_mod
      use ice_blocks, only: block, get_block, nx_block, ny_block, nghost, &
          nblocks_x, nblocks_y
      use ice_communicate, only: my_task, master_task
      use ice_constants, only: c0, c1, c2, p2, p5, c4
      use ice_domain_size, only: nx_global, ny_global, ncat, max_blocks, &
          nilyr, nslyr, nfsd, n_iso, n_aero
      use ice_domain, only: nblocks, blocks_ice, &
          ew_boundary_type, ns_boundary_type, &
          max_set_boundary_flds, num_set_boundary_flds, set_boundary_flds
      use ice_exit, only: abort_ice
      use ice_fileunits, only: nu_diag
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_parameters, &
          icepack_query_tracer_sizes, icepack_query_tracer_flags, &
          icepack_query_tracer_indices
      use ice_state, only: aicen, vicen, vsnon, trcrn, uvel, vvel
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_restore

      implicit none
      private
      public :: ice_restoring_init, ice_restoring_getdata, &
                ice_restoring_halo, ice_restoring_interior

      integer (kind=int_kind), parameter, public :: &
         max_restore_flds = 10

      integer (kind=int_kind), public :: &
         num_restore_flds

      !-----------------------------------------------------------------
      ! namelist inputs
      !-----------------------------------------------------------------

      logical (kind=log_kind), public :: &
         restore_ice                 ! restore ice state if true

      character (char_len_long), public :: &
         restore_data, &             ! restore data option
         restore_mask                ! restoring mask option

      integer (kind=int_kind), public :: &
         restore_width               ! restoring boundary depth (gridcells)

      real (kind=dbl_kind), public :: &
         restore_timescale           ! restoring timescale (days)

      character (char_len), dimension(max_restore_flds), public :: &
         restore_flds   ! set interior restoring for these fields

      !-----------------------------------------------------------------
      ! potential restoring fields
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable :: &
         aicen_restoring , & ! concentration of ice
         vicen_restoring , & ! volume per unit area of ice          (m)
         vsnon_restoring     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (:,:,:), allocatable :: &
         uvel_restoring  , & ! u velocity
         vvel_restoring  , & ! v velocity
         fval_restoring      ! restore forcing weight based on restoring time (0.-1.)

      real (kind=dbl_kind), dimension (:,:,:,:,:), allocatable :: &
         trcrn_restoring     ! tracers

      logical (kind=dbl_kind), dimension (:,:,:), allocatable :: &
         mask_restoring      ! mask where restoring is applied

!=======================================================================

      contains

!=======================================================================

!  Allocates and initializes arrays needed for restoring the ice state
!  in cells surrounding the grid.


   subroutine ice_restoring_init

      use ice_calendar, only: dt

      integer (kind=int_kind) :: &
         i,j,iblk          , & ! indices
         ilo, ihi, jlo, jhi, & ! physical domain indices
         igmin, igmax,       & ! restoring indices
         jgmin, jgmax,       & ! restoring indices
         ierr,               & ! error status
         minn,               & ! restoring index
         ntrcr                 ! number of tracers in use

      type (block) :: &
         this_block           ! block info for current block

      real (dbl_kind) :: &
         fval, &              ! fval computed from restore_timescale
         secday               ! seconds per day

      character(len=*), parameter :: subname = '(ice_restoring_init)'

      call icepack_query_parameters(secday_out=secday)
      call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      !-----------------------------------------------------------------------
      ! allocate and initialize
      !-----------------------------------------------------------------------

      allocate (aicen_restoring(nx_block,ny_block,ncat,max_blocks), &
                vicen_restoring(nx_block,ny_block,ncat,max_blocks), &
                vsnon_restoring(nx_block,ny_block,ncat,max_blocks), &
                trcrn_restoring(nx_block,ny_block,ntrcr,ncat,max_blocks), &
                stat=ierr)
      if (ierr/=0) call abort_ice(error_message=trim(subname)//' ERROR: Out of memory1', &
         file=__FILE__, line=__LINE__)

      allocate (uvel_restoring        (nx_block,ny_block,max_blocks), &
                vvel_restoring        (nx_block,ny_block,max_blocks), &
                fval_restoring        (nx_block,ny_block,max_blocks), &
                mask_restoring        (nx_block,ny_block,max_blocks), &
                stat=ierr)
      if (ierr/=0) call abort_ice(error_message=trim(subname)//' ERROR: Out of memory2', &
         file=__FILE__, line=__LINE__)

      aicen_restoring(:,:,:,:) = c0
      vicen_restoring(:,:,:,:) = c0
      vsnon_restoring(:,:,:,:) = c0
      trcrn_restoring(:,:,:,:,:) = c0

      uvel_restoring(:,:,:) = c0
      vvel_restoring(:,:,:) = c0

      ! default is restore everywhere instantaneously
      ! valid for halo

      mask_restoring(:,:,:) = .false.
      fval_restoring(:,:,:) = c0

      if (.not. restore_ice) return

      ! compute interior fval from restore_timescale

      if (restore_timescale == c0) then
         fval = c1
      else
         fval = min(abs(dt/(restore_timescale*secday)),c1)
      endif

      ! compute interior mask_restoring and fval_restoring from
      ! restore_mask and restore_width inputs

      igmin = restore_width
      igmax = nx_global-restore_width+1
      jgmin = restore_width
      jgmax = ny_global-restore_width+1
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi
         do j = jlo, jhi
         do i = ilo, ihi

            if (restore_mask == 'all') then
               ! set restoring everywhere
               mask_restoring(i,j,iblk) = .true.
               fval_restoring(i,j,iblk) = fval

            elseif (restore_mask == 'none') then
               ! turn off restoring everywhere
               mask_restoring(i,j,iblk) = .false.
               fval_restoring(i,j,iblk) = c0

            elseif (restore_mask == 'constant' .or. &
                    restore_mask == 'linear') then
               ! set restoring at interior gridcells associated with restore_width
               if (this_block%i_glob(i) > igmin .and. &
                   this_block%i_glob(i) < igmax .and. &
                   this_block%j_glob(j) > jgmin .and. &
                   this_block%j_glob(j) < jgmax) then
                  mask_restoring(i,j,iblk) = .false.
                  fval_restoring(i,j,iblk) = c0
               else
                  mask_restoring(i,j,iblk) = .true.
                  if (restore_mask == 'linear') then
                     minn = min(this_block%i_glob(i)-1, &
                                nx_global-this_block%i_glob(i), &
                                this_block%j_glob(j)-1, &
                                ny_global-this_block%j_glob(j))
                     fval_restoring(i,j,iblk) = (c1 - real(minn,dbl_kind) / real(restore_width,kind=dbl_kind)) * fval
                  else
                     fval_restoring(i,j,iblk) = fval
                  endif
!         write(100+my_task,'(2a,3i4,l4,e16.7)') subname,' fval check ',i,j,iblk,mask_restoring(i,j,iblk),fval_restoring(i,j,iblk)
               endif
!         write(100+my_task,'(2a,3i4,l4,e16.7)') subname,' fval check ',i,j,iblk,mask_restoring(i,j,iblk),fval_restoring(i,j,iblk)

            else
               call abort_ice(error_message=subname//' ERROR: restore_mask not supported '//trim(restore_mask), &
                  file=__FILE__, line=__LINE__)
            endif
         enddo
         enddo
      enddo

   end subroutine ice_restoring_init

!=======================================================================
! Sets restore arrays for outside bcs and internal restoring

   subroutine ice_restoring_getdata

      logical (log_kind), save :: &
         first_call = .true.                ! first call flag

      character(len=*), parameter :: subname = '(ice_restoring_getdata)'

      call ice_timer_start(timer_restore)

      if (restore_data == 'initial') then
         if (first_call) then
            ! set restore fields to initial ice state
            aicen_restoring(:,:,:,:) = aicen(:,:,:,:)
            vicen_restoring(:,:,:,:) = vicen(:,:,:,:)
            vsnon_restoring(:,:,:,:) = vsnon(:,:,:,:)
            trcrn_restoring(:,:,:,:,:) = trcrn(:,:,:,:,:)
            uvel_restoring(:,:,:) = uvel(:,:,:)
            vvel_restoring(:,:,:) = vvel(:,:,:)
         endif

      elseif (restore_data == 'defined') then
         if (first_call) then
            ! set restore fields to defined ice state
            call ice_restoring_data_define()
         endif

      elseif (restore_data == 'restartfiles') then
         call ice_restoring_data_restartfiles()

      else
         call abort_ice(error_message=subname//' ERROR: restore_data not supported '//trim(restore_data), &
            file=__FILE__, line=__LINE__)
      endif

      first_call = .false.

      call ice_timer_stop(timer_restore)

   end subroutine ice_restoring_getdata

!=======================================================================
! Sets restore arrays from individual restart files
! No time interpolation, one file per timestep

   subroutine ice_restoring_data_restartfiles

      use ice_calendar, only: msec, mday, mmonth, myear
      use ice_constants, only: field_loc_center, field_loc_necorner, &
          field_type_scalar, field_type_vector
      use ice_read_write, only: ice_read_nc, ice_open_nc, ice_close_nc

      integer(kind=int_kind) :: &
         k,           & ! dummy arguments
         nu_restoring     ! restore file id

      logical (kind=log_kind) :: skl_bgc, z_tracers, &
         tr_iage, tr_FY, tr_lvl, tr_pond, tr_aero, tr_fsd, &
         tr_snow, tr_brine, tr_iso

      integer (kind=int_kind) :: nt_Tsfc, nt_sice, nt_qice, nt_qsno, &
         nt_iage, nt_FY, nt_alvl, nt_vlvl, nt_apnd, nt_hpnd, nt_ipnd, nt_aero, &
         nt_fsd, nt_isosno, nt_isoice, nt_fbri, nt_smice, nt_smliq, nt_rhos, nt_rsnw

      character (char_len_long) :: &
         filebase, &         ! base filename
         data_file           ! data file to be read
      character (len=8) :: &
         nchar               ! temporary string
      logical (kind=log_kind) :: &
         diag                ! diagnose read call

      character(len=*), parameter :: subname = '(ice_restoring_data_restartfiles)'

      call icepack_query_parameters( &
         skl_bgc_out=skl_bgc, z_tracers_out=z_tracers)
      call icepack_query_tracer_flags(tr_iage_out=tr_iage, tr_FY_out=tr_FY, &
         tr_lvl_out=tr_lvl, tr_aero_out=tr_aero, tr_pond_out=tr_pond, &
         tr_brine_out=tr_brine, tr_fsd_out=tr_fsd, &
         tr_snow_out=tr_snow, tr_iso_out=tr_iso)
      call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_sice_out=nt_sice, &
         nt_qice_out=nt_qice, nt_qsno_out=nt_qsno, nt_iage_out=nt_iage, nt_fy_out=nt_fy, &
         nt_alvl_out=nt_alvl, nt_vlvl_out=nt_vlvl, nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, &
         nt_ipnd_out=nt_ipnd, nt_fsd_out=nt_fsd, nt_aero_out=nt_aero, &
         nt_smice_out=nt_smice, nt_smliq_out=nt_smliq, nt_rhos_out=nt_rhos, nt_rsnw_out=nt_rsnw, &
         nt_isosno_out=nt_isosno,     nt_isoice_out=nt_isoice,       nt_fbri_out=nt_fbri)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__, line=__LINE__)

      diag = .false.

      filebase = 'cice_bdy_restart_'
      write(data_file,'(a,i4.4,i2.2,i2.2,i5.5,a)') trim(filebase),myear,mmonth,mday,msec,'.nc'

      if (my_task == master_task) write (nu_diag,*) subname,' read ',trim(data_file)

      aicen_restoring(:,:,:,:) = c0
      vicen_restoring(:,:,:,:) = c0
      vsnon_restoring(:,:,:,:) = c0
      trcrn_restoring(:,:,:,:,:) = c0
      uvel_restoring(:,:,:) = c0
      vvel_restoring(:,:,:) = c0

      call ice_open_nc(data_file, nu_restoring)

      ! aicen, vicen, vsnon

      call ice_read_nc(nu_restoring, 1, 'aicen', aicen_restoring, diag, restart_ext=.true.)
      call ice_read_nc(nu_restoring, 1, 'vicen', vicen_restoring, diag, restart_ext=.true.)
      call ice_read_nc(nu_restoring, 1, 'vsnon', vsnon_restoring, diag, restart_ext=.true.)

      ! velocity

      call ice_read_nc(nu_restoring, 1, 'uvel', uvel_restoring, diag, restart_ext=.true.)
      call ice_read_nc(nu_restoring, 1, 'vvel', vvel_restoring, diag, restart_ext=.true.)

      ! Tsfcn

      call ice_read_nc(nu_restoring, 1, 'Tsfcn', trcrn_restoring(:,:,nt_Tsfc,:,:), diag, restart_ext=.true.)

      ! sice

      do k=1,nilyr
         write(nchar,'(i3.3)') k
         call ice_read_nc(nu_restoring, 1, 'sice'//trim(nchar), trcrn_restoring(:,:,nt_sice+k-1,:,:), diag, restart_ext=.true.)
      enddo

      ! qice

      do k=1,nilyr
         write(nchar,'(i3.3)') k
         call ice_read_nc(nu_restoring, 1, 'qice'//trim(nchar), trcrn_restoring(:,:,nt_qice+k-1,:,:), diag, restart_ext=.true.)
      enddo

      ! qsno

      do k=1,nslyr
         write(nchar,'(i3.3)') k
         call ice_read_nc(nu_restoring, 1, 'qsno'//trim(nchar), trcrn_restoring(:,:,nt_qsno+k-1,:,:), diag, restart_ext=.true.)
      enddo

      ! iage

      if (tr_iage) then
         call ice_read_nc(nu_restoring, 1, 'iage', trcrn_restoring(:,:,nt_iage,:,:), diag, restart_ext=.true.)
      endif

      ! FY

      if (tr_FY) then
         call ice_read_nc(nu_restoring, 1, 'FY', trcrn_restoring(:,:,nt_FY,:,:), diag, restart_ext=.true.)
      endif

      ! lvl

      if (tr_lvl) then
         call ice_read_nc(nu_restoring, 1, 'alvl', trcrn_restoring(:,:,nt_alvl,:,:), diag, restart_ext=.true.)
         call ice_read_nc(nu_restoring, 1, 'vlvl', trcrn_restoring(:,:,nt_vlvl,:,:), diag, restart_ext=.true.)
      endif

      ! pond lvl, pond topo, pond sealvl

      if (tr_pond) then
         call ice_read_nc(nu_restoring, 1, 'apnd',  trcrn_restoring(:,:,nt_apnd,:,:), diag, restart_ext=.true.)
         call ice_read_nc(nu_restoring, 1, 'hpnd',  trcrn_restoring(:,:,nt_hpnd,:,:), diag, restart_ext=.true.)
         call ice_read_nc(nu_restoring, 1, 'ipnd',  trcrn_restoring(:,:,nt_ipnd,:,:), diag, restart_ext=.true.)
      endif

      ! snow

      if (tr_snow) then
         do k=1,nslyr
            write(nchar,'(i3.3)') k
            call ice_read_nc(nu_restoring, 1, 'smice'//trim(nchar), trcrn_restoring(:,:,nt_smice+k-1,:,:), diag, restart_ext=.true.)
            call ice_read_nc(nu_restoring, 1, 'smliq'//trim(nchar), trcrn_restoring(:,:,nt_smliq+k-1,:,:), diag, restart_ext=.true.)
            call ice_read_nc(nu_restoring, 1, 'rhos'//trim(nchar), trcrn_restoring(:,:,nt_rhos+k-1,:,:), diag, restart_ext=.true.)
            call ice_read_nc(nu_restoring, 1, 'rsnw'//trim(nchar), trcrn_restoring(:,:,nt_rsnw+k-1,:,:), diag, restart_ext=.true.)
         enddo
      endif

      ! fsd

      if (tr_fsd) then
         do k=1,nfsd
            write(nchar,'(i3.3)') k
            call ice_read_nc(nu_restoring, 1, 'fsd'//trim(nchar), trcrn_restoring(:,:,nt_fsd+k-1,:,:), diag, restart_ext=.true.)
         enddo
      endif

      ! iso

      if (tr_iso) then
          do k = 1, n_iso
             write(nchar,'(i3.3)') k
             call ice_read_nc(nu_restoring, 1, 'isosno'//trim(nchar), trcrn_restoring(:,:,nt_isosno+k-1,:,:), diag, restart_ext=.true.)
             call ice_read_nc(nu_restoring, 1, 'isoice'//trim(nchar), trcrn_restoring(:,:,nt_isoice+k-1,:,:), diag, restart_ext=.true.)
          enddo
      endif

      ! aero

      if (tr_aero) then
         do k = 1, n_aero
            write(nchar,'(i3.3)') k
            call ice_read_nc(nu_restoring, 1, 'aerosnossl'//trim(nchar), trcrn_restoring(:,:,nt_aero  +(k-1)*4,:,:), diag, restart_ext=.true.)
            call ice_read_nc(nu_restoring, 1, 'aerosnoint'//trim(nchar), trcrn_restoring(:,:,nt_aero+1+(k-1)*4,:,:), diag, restart_ext=.true.)
            call ice_read_nc(nu_restoring, 1, 'aeroicessl'//trim(nchar), trcrn_restoring(:,:,nt_aero+2+(k-1)*4,:,:), diag, restart_ext=.true.)
            call ice_read_nc(nu_restoring, 1, 'aeroiceint'//trim(nchar), trcrn_restoring(:,:,nt_aero+3+(k-1)*4,:,:), diag, restart_ext=.true.)
         enddo
      endif

      ! brine

      if (tr_brine) then
         call ice_read_nc(nu_restoring, 1, 'fbrn', trcrn_restoring(:,:,nt_fbri,:,:), diag, restart_ext=.true.)
      endif

      ! bgc - To Be Done

      if (skl_bgc .or. z_tracers) then
         call abort_ice(error_message=subname//' ERROR: bgc not supported yet ', &
            file=__FILE__, line=__LINE__)
      endif

      call ice_close_nc(nu_restoring)

   end subroutine ice_restoring_data_restartfiles

!=======================================================================
! define restoring variables for halo and interior

   subroutine ice_restoring_data_define

      use ice_arrays_column, only: hin_max
      use ice_flux, only: Tf, Tair, salinz, Tmltz
      use ice_grid, only: tmask,hm
      use icepack_intfc, only: icepack_init_trcr

      ! local variables

      integer (kind=int_kind) :: &
         ilo, ihi          , & ! physical domain indices
         jlo, jhi          , & !
         iglob(nx_block)   , & ! global indices
         jglob(ny_block)   , & !
         iblock            , & ! block indices
         jblock            , & !
         ntrcr                 ! number of tracers in use

      integer (kind=int_kind) :: &
         i,j,k,n,nt,iblk   , & ! indices
         nt_Tsfc           , & !
         nt_fbri           , & !
         nt_qice           , & !
         nt_sice           , & !
         nt_qsno           , & !
         icells                ! number of cells initialized with ice

      type (block) :: &
         this_block           ! block info for current block

      logical (kind=log_kind) :: &
         tr_brine

      real (kind=dbl_kind) :: &
         Tsfc, hbar, &
         hsno_init       ! initial snow thickness

      real (kind=dbl_kind), dimension(ncat) :: &
         ainit, hinit    ! initial area, thickness

      real (kind=dbl_kind), dimension(nilyr) :: &
         qin             ! ice enthalpy (J/m3)

      real (kind=dbl_kind), dimension(nslyr) :: &
         qsn             ! snow enthalpy (J/m3)

      character(len=*), parameter :: subname = '(ice_restoring_data_define)'

      call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
      call icepack_query_tracer_flags(tr_brine_out=tr_brine)
      call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_fbri_out=nt_fbri, &
           nt_qice_out=nt_qice, nt_sice_out=nt_sice, nt_qsno_out=nt_qsno)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      !-----------------------------------------------------------------
      ! initial area and thickness in ice-occupied restoring cells
      !-----------------------------------------------------------------

      hbar = c2  ! initial ice thickness
      hsno_init = 0.20_dbl_kind ! initial snow thickness (m)
      do n = 1, ncat
         hinit(n) = c0
         ainit(n) = c0
         if (hbar > hin_max(n-1) .and. hbar < hin_max(n)) then
            hinit(n) = hbar
            ainit(n) = 0.95_dbl_kind ! initial ice concentration
         endif
      enddo

      !-----------------------------------------------------------------
      ! Initialize restoring variables everywhere on grid
      !-----------------------------------------------------------------

      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block
            aicen_restoring(i,j,n,iblk) = c0
            vicen_restoring(i,j,n,iblk) = c0
            vsnon_restoring(i,j,n,iblk) = c0
            if (tmask(i,j,iblk)) then
               trcrn_restoring(i,j,nt_Tsfc,n,iblk) = Tf(i,j,iblk)  ! surface temperature
            else
               trcrn_restoring(i,j,nt_Tsfc,n,iblk) = c0  ! on land gridcells
            endif
            if (ntrcr >= 2) then
               do nt = 2, ntrcr
                  trcrn_restoring(i,j,nt,n,iblk) = c0
               enddo
            endif
            if (tr_brine) trcrn_restoring(i,j,nt_fbri,n,iblk) = c1
         enddo
         enddo
         enddo
      enddo

      !-----------------------------------------------------------------
      ! Define cells where ice is placed (or other values are used)
      ! Edges using initial values (zero, above) are commented out
      !-----------------------------------------------------------------

      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi
         iglob = this_block%i_glob
         jglob = this_block%j_glob
         iblock = this_block%iblock
         jblock = this_block%jblock

         ! mask with iglob, jglob or tlon, tlat as needed

         do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block

               ! ice volume, snow volume
               aicen_restoring(i,j,n,iblk) = ainit(n)
               vicen_restoring(i,j,n,iblk) = hinit(n) * ainit(n) ! m
               vsnon_restoring(i,j,n,iblk) = min(aicen_restoring(i,j,n,iblk)*hsno_init, &
                                                 p2*vicen_restoring(i,j,n,iblk))

               call icepack_init_trcr(Tair=Tair(i,j,iblk), Tf=Tf(i,j,iblk), &
                                      Sprofile=salinz(i,j,:,iblk),          &
                                      Tprofile=Tmltz(i,j,:,iblk),           &
                                      Tsfc=Tsfc,                            &
                                      qin=qin(:),        qsn=qsn(:))

               ! surface temperature
               trcrn_restoring(i,j,nt_Tsfc,n,iblk) = Tsfc ! deg C
               ! ice enthalpy, salinity
               do k = 1, nilyr
                  trcrn_restoring(i,j,nt_qice+k-1,n,iblk) = qin(k)
                  trcrn_restoring(i,j,nt_sice+k-1,n,iblk) = salinz(i,j,k,iblk)
               enddo
               ! snow enthalpy
               do k = 1, nslyr
                  trcrn_restoring(i,j,nt_qsno+k-1,n,iblk) = qsn(k)
               enddo

         enddo  ! i
         enddo  ! j
         enddo  ! ncat
      enddo ! blocks

      !-----------------------------------------------------------------
      ! Impose land mask
      !-----------------------------------------------------------------

      do iblk = 1, nblocks
         do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block
            aicen_restoring(i,j,n,iblk) = aicen_restoring(i,j,n,iblk) * hm(i,j,iblk)
            vicen_restoring(i,j,n,iblk) = vicen_restoring(i,j,n,iblk) * hm(i,j,iblk)
            vsnon_restoring(i,j,n,iblk) = vsnon_restoring(i,j,n,iblk) * hm(i,j,iblk)
            do nt = 1, ntrcr
               trcrn_restoring(i,j,nt,n,iblk) = trcrn_restoring(i,j,nt,n,iblk) &
                                                            * hm(i,j,iblk)
            enddo
         enddo
         enddo
         enddo
      enddo

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

   end subroutine ice_restoring_data_define

!=======================================================================
!  This subroutine is intended for restoring the ice state to desired
!  values in the interior of the grid based on fval_restoring weights

   subroutine ice_restoring_interior(setfld)

      character (len=*), intent(in) :: &
         setfld                ! fields to restore

      ! local variables

      integer (int_kind) :: &
         iblk,n,k,           & ! dummy loop indices
         ilo,ihi,jlo,jhi,    & ! beginning and end of physical domain
         nvi, nvj              ! adjust restoring index for velocity (corner) field

      type (block) :: &
         this_block  ! block info for current block

      character (len=char_len) :: &
         psetfld               ! passed field name

      logical (log_kind), save :: &
         fldflag(max_restore_flds)  ! matched 

      logical (log_kind), save :: &
         first_call=.true.               ! initialize fldflag

      character(len=*), parameter :: subname = '(ice_restoring_interior)'

      if (.not. restore_ice) return

      ! return if no fields set by user
      if (num_restore_flds == 0) return

      ! ignore setfld = '' or 'none'
      if (setfld == '' .or. setfld == 'none') return

      call ice_timer_start(timer_restore)

      !-----------------------------------------------------------------------
      ! Manage and document usage
      !-----------------------------------------------------------------------

      if (first_call) then
         first_call = .false.
         do n = 1,max_restore_flds
            if (restore_flds(n) == '' .or. restore_flds(n) == 'none') then
               fldflag(n) = .true.   ! don't check if '' or 'none'
            else
               fldflag(n) = .false.
            endif
         enddo
      endif

      !-----------------------------------------------------------------------
      !  Restore interior
      !-----------------------------------------------------------------------

      do n = 1,num_restore_flds

         ! Match restore call flds to those set by user
         psetfld = 'none'
         if (setfld == restore_flds(n)) then
            psetfld = setfld
         elseif (setfld == 'state' .and. restore_flds(n) == 'aicen') then
            psetfld = 'aicen'
         elseif (setfld == 'state' .and. restore_flds(n) == 'vicen') then
            psetfld = 'vicen'
         elseif (setfld == 'state' .and. restore_flds(n) == 'vsnon') then
            psetfld = 'vsnon'
         elseif (setfld == 'state' .and. restore_flds(n) == 'trcrn') then
            psetfld = 'trcrn'
         endif

         if (psetfld /= 'none') then

            if (.not.fldflag(n)) then
               if (my_task == master_task) &
                              write(nu_diag,*) subname,' setting halo field '//trim(psetfld)
               fldflag(n) = .true.
               do k = 1,num_restore_flds
                  if (.not.fldflag(k)) then
                     if (my_task == master_task) &
                              write(nu_diag,*) subname,' waiting to set '//trim(restore_flds(k))
                  endif
               enddo
            endif

            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi

               nvi = 0
               nvj = 0
               ! reduce interior restoring by 1 gridcell on east and north edge
               if (psetfld == 'velocity') then
                  if (this_block%iblock == nblocks_x) nvi = -1
                  if (this_block%jblock == nblocks_y) nvj = -1
               endif

               call restore_cells(iblk, ilo, ihi+nvi, jlo, jhi+nvj, &
                    halo=.false., setfld=psetfld)

            enddo

         endif
      enddo

      call ice_timer_stop(timer_restore)

   end subroutine ice_restoring_interior

!=======================================================================
!  This subroutine is intended for restoring the ice state to desired
!  values in halo cells surrounding the grid.

   subroutine ice_restoring_halo(setfld)

      character (len=*), intent(in) :: &
         setfld                ! field to restore

      ! local variables

      integer (int_kind) :: &
         iblk,n,k,           & ! dummy loop indices
         ilo,ihi,jlo,jhi,    & ! beginning and end of physical domain
         nv                    ! adjust indexing for velocity (corner) field

      type (block) :: &
         this_block  ! block info for current block

      character (len=char_len) :: &
         psetfld               ! passed field name

      logical (log_kind), save :: &
         fldflag(max_set_boundary_flds)  ! matched 

      logical (log_kind), save :: &
         first_call=.true.               ! initialize fldflag

      character(len=*), parameter :: subname = '(ice_restoring_halo)'

      !-----------------------------------------------------------------------
      ! Return conditions to speed up model
      !-----------------------------------------------------------------------

      ! return if no fields set by user
      if (num_set_boundary_flds == 0) return

      ! ignore cyclic and tripole bcs
      if (trim(ew_boundary_type) == 'cyclic' .and. &
          (trim(ns_boundary_type) == 'cyclic' .or. &
           trim(ns_boundary_type) == 'tripole' .or. &
           trim(ns_boundary_type) == 'tripoleT')) then
         return
      endif

      ! ignore setfld = '' or 'none'
      if (setfld == '' .or. setfld == 'none') return

      call ice_timer_start(timer_restore)

      !-----------------------------------------------------------------------
      ! Manage and document usage
      !-----------------------------------------------------------------------

      if (first_call) then
         first_call = .false.
         do n = 1,max_set_boundary_flds
            if (set_boundary_flds(n) == '' .or. set_boundary_flds(n) == 'none') then
               fldflag(n) = .true.   ! don't check if '' or 'none'
            else
               fldflag(n) = .false.
            endif
         enddo
      endif

      !-----------------------------------------------------------------------
      ! Set outer boundary halo values
      !-----------------------------------------------------------------------

      do n = 1,num_set_boundary_flds
         ! Match halo restore call flds to those set by user
         psetfld = 'none'
         if (setfld == set_boundary_flds(n)) then
            psetfld = setfld
         elseif (setfld == 'state' .and. set_boundary_flds(n) == 'aicen') then
            psetfld = 'aicen'
         elseif (setfld == 'state' .and. set_boundary_flds(n) == 'vicen') then
            psetfld = 'vicen'
         elseif (setfld == 'state' .and. set_boundary_flds(n) == 'vsnon') then
            psetfld = 'vsnon'
         elseif (setfld == 'state' .and. set_boundary_flds(n) == 'trcrn') then
            psetfld = 'trcrn'
         endif

         nv = 0
         if (psetfld == 'velocity') nv = -1

         if (psetfld /= 'none') then

            if (.not.fldflag(n)) then
               if (my_task == master_task) &
                              write(nu_diag,*) subname,' setting halo field '//trim(psetfld)
               fldflag(n) = .true.
               do k = 1,num_set_boundary_flds
                  if (.not.fldflag(k)) then
                     if (my_task == master_task) &
                              write(nu_diag,*) subname,' waiting to set '//trim(set_boundary_flds(k))
                  endif
               enddo
            endif

            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi

               if (trim(ew_boundary_type) /= 'cyclic') then
                  ! West Edge
                  if (this_block%iblock == 1) then
                     call restore_cells(iblk, 1, ilo-1, 1, ny_block, &
                          halo=.true., setfld=psetfld)
                  endif
                  ! East Edge
                  if (this_block%iblock == nblocks_x) then
                     call restore_cells(iblk, ihi+1+nv, ihi+nghost, 1, ny_block, &
                          halo=.true., setfld=psetfld)
                  endif
               endif

               ! South Edge
               if (trim(ns_boundary_type) /= 'cyclic') then
                  if (this_block%jblock == 1) then
                     call restore_cells(iblk, 1, nx_block, 1, jlo-1, &
                          halo=.true., setfld=psetfld)
                  endif
               endif

               ! North Edge
               if (trim(ns_boundary_type) /= 'cyclic' .and. &
                   trim(ns_boundary_type) /= 'tripole' .and. &
                   trim(ns_boundary_type) /= 'tripoleT') then
                  if (this_block%jblock == nblocks_y) then
                     call restore_cells(iblk, 1, nx_block, jhi+1+nv, jhi+nghost, &
                          halo=.true., setfld=psetfld)
                  endif
               endif

            enddo ! iblk

         endif
      enddo ! n

      call ice_timer_stop(timer_restore)

   end subroutine ice_restoring_halo

!=======================================================================

   subroutine restore_cells(iblk,i1,i2,j1,j2,halo,setfld)

      integer(kind=int_kind), intent(in) :: &
         iblk,      & ! block id
         i1, i2,    & ! i start and end indices
         j1, j2       ! j start and end indices

      logical (kind=log_kind), intent(in) :: &
         halo         ! is this for the halo

      character(len=*), intent(in) :: &
         setfld       ! fields to restore

      ! local variable

      integer(kind=int_kind) :: &
         i, j, n, nt, ntrcr ! local indices

      real(kind=dbl_kind) :: &
         vr1, vv2     ! restoring weights

      character(len=*), parameter :: subname = '(restore_cells)'

      if (setfld == 'none') return

      call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (setfld == 'state') then

         if (halo) then   ! set values explicitly
            do i = i1,i2
            do j = j1,j2
               do n = 1, ncat
                  aicen (i,j,n,iblk) = aicen_restoring(i,j,n,iblk)
                  vicen (i,j,n,iblk) = vicen_restoring(i,j,n,iblk)
                  vsnon (i,j,n,iblk) = vsnon_restoring(i,j,n,iblk)
                  do nt = 1, ntrcr
                     trcrn(i,j,nt,n,iblk) = trcrn_restoring(i,j,nt,n,iblk)
                  enddo
               enddo
            enddo
            enddo
         else              ! restoring
            do i = i1,i2
            do j = j1,j2
               if (mask_restoring(i,j,iblk)) then
                  vr1 = fval_restoring(i,j,iblk)
                  vv2 = c1 - vr1
                  do n = 1, ncat
                     aicen (i,j,n,iblk) = aicen_restoring(i,j,n,iblk)*vr1 + aicen (i,j,n,iblk)*vv2
                     vicen (i,j,n,iblk) = vicen_restoring(i,j,n,iblk)*vr1 + vicen (i,j,n,iblk)*vv2
                     vsnon (i,j,n,iblk) = vsnon_restoring(i,j,n,iblk)*vr1 + vsnon (i,j,n,iblk)*vv2
                     do nt = 1, ntrcr
                        trcrn(i,j,nt,n,iblk) = trcrn_restoring(i,j,nt,n,iblk)*vr1 + trcrn(i,j,nt,n,iblk)*vv2
                     enddo
                  enddo
               endif
            enddo
            enddo
         endif

      elseif (setfld == 'aicen') then

         if (halo) then   ! set values explicitly
            do i = i1,i2
            do j = j1,j2
               do n = 1, ncat
                  aicen (i,j,n,iblk) = aicen_restoring(i,j,n,iblk)
               enddo
            enddo
            enddo
         else              ! restoring
            do i = i1,i2
            do j = j1,j2
               if (mask_restoring(i,j,iblk)) then
                  vr1 = fval_restoring(i,j,iblk)
                  vv2 = c1 - vr1
                  do n = 1, ncat
                     aicen (i,j,n,iblk) = aicen_restoring(i,j,n,iblk)*vr1 + aicen (i,j,n,iblk)*vv2
                  enddo
               endif
            enddo
            enddo
         endif

      elseif (setfld == 'vicen') then

         if (halo) then   ! set values explicitly
            do i = i1,i2
            do j = j1,j2
               do n = 1, ncat
                  vicen (i,j,n,iblk) = vicen_restoring(i,j,n,iblk)
               enddo
            enddo
            enddo
         else              ! restoring
            do i = i1,i2
            do j = j1,j2
               if (mask_restoring(i,j,iblk)) then
                  vr1 = fval_restoring(i,j,iblk)
                  vv2 = c1 - vr1
                  do n = 1, ncat
                     vicen (i,j,n,iblk) = vicen_restoring(i,j,n,iblk)*vr1 + vicen (i,j,n,iblk)*vv2
                  enddo
               endif
            enddo
            enddo
         endif

      elseif (setfld == 'vsnon') then

         ! center gridcell
         if (halo) then   ! set values explicitly
            do i = i1,i2
            do j = j1,j2
               do n = 1, ncat
                  vsnon (i,j,n,iblk) = vsnon_restoring(i,j,n,iblk)
               enddo
            enddo
            enddo
         else              ! restoring
            do i = i1,i2
            do j = j1,j2
               if (mask_restoring(i,j,iblk)) then
                  vr1 = fval_restoring(i,j,iblk)
                  vv2 = c1 - vr1
                  do n = 1, ncat
                     vsnon (i,j,n,iblk) = vsnon_restoring(i,j,n,iblk)*vr1 + vsnon (i,j,n,iblk)*vv2
                  enddo
               endif
            enddo
            enddo
         endif

      elseif (setfld == 'trcrn') then

         if (halo) then   ! set values explicitly
            do i = i1,i2
            do j = j1,j2
               do n = 1, ncat
                  do nt = 1, ntrcr
                     trcrn(i,j,nt,n,iblk) = trcrn_restoring(i,j,nt,n,iblk)
                  enddo
               enddo
            enddo
            enddo
         else              ! restoring
            do i = i1,i2
            do j = j1,j2
               if (mask_restoring(i,j,iblk)) then
                  vr1 = fval_restoring(i,j,iblk)
                  vv2 = c1 - vr1
                  do n = 1, ncat
                     do nt = 1, ntrcr
                        trcrn(i,j,nt,n,iblk) = trcrn_restoring(i,j,nt,n,iblk)*vr1 + trcrn(i,j,nt,n,iblk)*vv2
                     enddo
                  enddo
               endif
            enddo
            enddo
         endif

      elseif (setfld == 'velocity') then

         if (halo) then   ! set values explicitly
            do i = i1,i2
            do j = j1,j2
               uvel(i,j,iblk) = uvel_restoring(i,j,iblk)
               vvel(i,j,iblk) = vvel_restoring(i,j,iblk)
            enddo
            enddo
         else              ! restoring
            do i = i1, i2
            do j = j1, j2
               if (mask_restoring(i,j,iblk)) then
                  vr1 = fval_restoring(i,j,iblk)
                  vv2 = c1 - vr1
                  uvel(i,j,iblk) = uvel_restoring(i,j,iblk)*vr1 + uvel(i,j,iblk)*vv2
                  vvel(i,j,iblk) = vvel_restoring(i,j,iblk)*vr1 + vvel(i,j,iblk)*vv2
               endif
            enddo
            enddo
         endif

      else

         call abort_ice(error_message=subname//' ERROR: setfld unknown = '//trim(setfld), &
            file=__FILE__, line=__LINE__)

      endif

   end subroutine restore_cells
!=======================================================================

   end module ice_restoring

!=======================================================================
