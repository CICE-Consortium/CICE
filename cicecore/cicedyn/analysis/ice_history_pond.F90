!=======================================================================

! Melt pond history output
!
! 2012 Elizabeth Hunke split code from ice_history.F90

      module ice_history_pond

      use ice_kinds_mod
      use ice_domain_size, only: max_nstrm
      use ice_constants, only: c0, c1
      use ice_fileunits, only: nu_nml, nml_filename, &
          get_fileunit, release_fileunit
      use ice_fileunits, only: nu_diag
      use ice_exit, only: abort_ice
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_parameters, &
          icepack_query_tracer_flags, icepack_query_tracer_indices

      implicit none
      private
      public :: accum_hist_pond, init_hist_pond_2D, init_hist_pond_3Dc

      !---------------------------------------------------------------
      ! flags: write to output file if true or histfreq value
      !---------------------------------------------------------------

      character (len=max_nstrm), public :: &
           f_apondn      = 'm', f_apeffn       = 'm', &
           f_hpondn      = 'm',                     &
           f_apond       = 'x', f_apond_ai     = 'x', &
           f_hpond       = 'x', f_hpond_ai     = 'x', &
           f_ipond       = 'x', f_ipond_ai     = 'x', &
           f_apeff       = 'x', f_apeff_ai     = 'x', &
           f_dpnd_flush  = 'x', f_dpnd_expon   = 'x', &
           f_dpnd_freebd = 'x', f_dpnd_initial = 'x', &
           f_dpnd_dlid   = 'x', f_dpnd_melt    = 'x', &
           f_dpnd_ridge  = 'x',                     &
           f_dpnd_flushn = 'x', f_dpnd_exponn  = 'x', &
           f_dpnd_freebdn= 'x', f_dpnd_initialn= 'x', &
           f_dpnd_dlidn  = 'x'

      !---------------------------------------------------------------
      ! namelist variables
      !---------------------------------------------------------------

      namelist / icefields_pond_nml /     &
           f_apondn,    f_apeffn   , &
           f_hpondn,                 &
           f_apond,     f_apond_ai , &
           f_hpond,     f_hpond_ai , &
           f_ipond,     f_ipond_ai , &
           f_apeff,     f_apeff_ai , &
           f_dpnd_flush    , &
           f_dpnd_expon    , &
           f_dpnd_freebd   , &
           f_dpnd_initial  , &
           f_dpnd_dlid     , &
           f_dpnd_melt     , &
           f_dpnd_ridge    , &
           f_dpnd_flushn   , &
           f_dpnd_exponn   , &
           f_dpnd_freebdn  , &
           f_dpnd_initialn , &
           f_dpnd_dlidn

      !---------------------------------------------------------------
      ! field indices
      !---------------------------------------------------------------

      integer (kind=int_kind), dimension(max_nstrm) :: &
           n_apondn      , n_apeffn    , &
           n_hpondn      , &
           n_apond       , n_apond_ai, &
           n_hpond       , n_hpond_ai, &
           n_ipond       , n_ipond_ai, &
           n_apeff       , n_apeff_ai, &
           n_dpnd_flush  , n_dpnd_expon, &
           n_dpnd_freebd , n_dpnd_initial, &
           n_dpnd_dlid   , n_dpnd_melt, &
           n_dpnd_ridge  , &
           n_dpnd_flushn , n_dpnd_exponn, &
           n_dpnd_freebdn, n_dpnd_initialn, &
           n_dpnd_dlidn

!=======================================================================

      contains

!=======================================================================

      subroutine init_hist_pond_2D

      use ice_broadcast, only: broadcast_scalar
      use ice_calendar, only: nstreams, histfreq
      use ice_communicate, only: my_task, master_task
      use ice_history_shared, only: tstr2D, tcstr, define_hist_field
      use ice_fileunits, only: goto_nml

      integer (kind=int_kind) :: ns
      integer (kind=int_kind) :: nml_error ! namelist i/o error flag
      logical (kind=log_kind) :: tr_pond
      character(len=char_len_long) :: tmpstr2 ! for namelist check
      character(len=char_len)      :: nml_name ! text namelist name

      character(len=*), parameter :: subname = '(init_hist_pond_2D)'

      call icepack_query_tracer_flags(tr_pond_out=tr_pond)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      !-----------------------------------------------------------------
      ! read namelist
      !-----------------------------------------------------------------

      if (my_task == master_task) then
         nml_name = 'icefields_pond_nml'
         write(nu_diag,*) subname,' Reading ', trim(nml_name)

         ! open namelist file
         call get_fileunit(nu_nml)
         open (nu_nml, file=trim(nml_filename), status='old',iostat=nml_error)
         if (nml_error /= 0) then
            call abort_ice(subname//'ERROR: '//trim(nml_name)//' open file '// &
                 trim(nml_filename), &
                 file=__FILE__, line=__LINE__)
         endif

         ! goto this namelist in file
         call goto_nml(nu_nml,trim(nml_name),nml_error)
         if (nml_error /= 0) then
            call abort_ice(subname//'ERROR: searching for '// trim(nml_name), &
                 file=__FILE__, line=__LINE__)
         endif

         ! read namelist
         nml_error =  1
         do while (nml_error > 0)
            read(nu_nml, nml=icefields_pond_nml,iostat=nml_error)
            ! check if error
            if (nml_error /= 0) then
               ! backspace and re-read erroneous line
               backspace(nu_nml)
               read(nu_nml,fmt='(A)') tmpstr2
               call abort_ice(subname//'ERROR: ' // trim(nml_name) // ' reading ' // &
                    trim(tmpstr2), file=__FILE__, line=__LINE__)
            endif
         end do

         close(nu_nml)
         call release_fileunit(nu_nml)
      endif

      if (.not. tr_pond) then
          f_apondn    = 'x'
          f_hpondn    = 'x'
          f_apeffn    = 'x'
          f_apond     = 'x'
          f_hpond     = 'x'
          f_ipond     = 'x'
          f_apeff     = 'x'
          f_apond_ai  = 'x'
          f_hpond_ai  = 'x'
          f_ipond_ai  = 'x'
          f_apeff_ai  = 'x'
          f_dpnd_flush   = 'x'
          f_dpnd_expon   = 'x'
          f_dpnd_freebd  = 'x'
          f_dpnd_initial = 'x'
          f_dpnd_dlid    = 'x'
          f_dpnd_melt    = 'x'
          f_dpnd_ridge   = 'x'
          f_dpnd_flushn  = 'x'
          f_dpnd_exponn  = 'x'
          f_dpnd_freebdn = 'x'
          f_dpnd_initialn= 'x'
          f_dpnd_dlidn   = 'x'
      endif

      call broadcast_scalar (f_apondn, master_task)
      call broadcast_scalar (f_hpondn, master_task)
      call broadcast_scalar (f_apeffn, master_task)
      call broadcast_scalar (f_apond, master_task)
      call broadcast_scalar (f_hpond, master_task)
      call broadcast_scalar (f_ipond, master_task)
      call broadcast_scalar (f_apeff, master_task)
      call broadcast_scalar (f_apond_ai, master_task)
      call broadcast_scalar (f_hpond_ai, master_task)
      call broadcast_scalar (f_ipond_ai, master_task)
      call broadcast_scalar (f_apeff_ai, master_task)
      call broadcast_scalar (f_dpnd_flush   , master_task)
      call broadcast_scalar (f_dpnd_expon   , master_task)
      call broadcast_scalar (f_dpnd_freebd  , master_task)
      call broadcast_scalar (f_dpnd_initial , master_task)
      call broadcast_scalar (f_dpnd_dlid    , master_task)
      call broadcast_scalar (f_dpnd_melt    , master_task)
      call broadcast_scalar (f_dpnd_ridge   , master_task)
      call broadcast_scalar (f_dpnd_flushn  , master_task)
      call broadcast_scalar (f_dpnd_exponn  , master_task)
      call broadcast_scalar (f_dpnd_freebdn , master_task)
      call broadcast_scalar (f_dpnd_initialn, master_task)
      call broadcast_scalar (f_dpnd_dlidn   , master_task)

      if (tr_pond) then

      ! 2D variables
      do ns = 1, nstreams
      if (histfreq(ns) /= 'x') then

      if (f_apond(1:1) /= 'x') &
         call define_hist_field(n_apond,"apond","1",tstr2D, tcstr, &
             "melt pond fraction of sea ice",                      &
             "none", c1, c0,                                       &
             ns, f_apond)

      if (f_apond_ai(1:1) /= 'x') &
         call define_hist_field(n_apond_ai,"apond_ai","1",tstr2D, tcstr, &
             "melt pond fraction of grid cell",                    &
             "weighted by ice area", c1, c0,                       &
             ns, f_apond_ai)

      if (f_hpond(1:1) /= 'x') &
         call define_hist_field(n_hpond,"hpond","m",tstr2D, tcstr, &
             "mean melt pond depth over sea ice",                  &
             "none", c1, c0,                                       &
             ns, f_hpond)

      if (f_hpond_ai(1:1) /= 'x') &
         call define_hist_field(n_hpond_ai,"hpond_ai","m",tstr2D, tcstr, &
             "mean melt pond depth over grid cell",                &
             "weighted by ice area", c1, c0,                       &
             ns, f_hpond)

      if (f_ipond(1:1) /= 'x') &
         call define_hist_field(n_ipond,"ipond","m",tstr2D, tcstr, &
             "mean pond ice thickness over sea ice",               &
             "none", c1, c0,                                       &
             ns, f_ipond)

      if (f_ipond_ai(1:1) /= 'x') &
         call define_hist_field(n_ipond_ai,"ipond_ai","m",tstr2D, tcstr, &
             "mean pond ice thickness over grid cell",             &
             "weighted by ice area", c1, c0,                       &
             ns, f_ipond_ai)

      if (f_apeff(1:1) /= 'x') &
         call define_hist_field(n_apeff,"apeff","1",tstr2D, tcstr, &
             "radiation-effective pond area fraction of sea ice",  &
             "none", c1, c0,  &
             ns, f_apeff)

      if (f_apeff_ai(1:1) /= 'x') &
         call define_hist_field(n_apeff_ai,"apeff_ai","1",tstr2D, tcstr, &
             "radiation-effective pond area fraction over grid cell",  &
             "weighted by ice area", c1, c0,                       &
             ns, f_apeff_ai)

      if (f_dpnd_flush(1:1) /= 'x') &
         call define_hist_field(n_dpnd_flush,"dpnd_flush","m/s",tstr2D, tcstr, &
             "pond flushing rate due to ice permeability",           &
             "none", c1, c0,                                       &
             ns, f_dpnd_flush)

      if (f_dpnd_expon(1:1) /= 'x') &
         call define_hist_field(n_dpnd_expon,"dpnd_expon","m/s",tstr2D, tcstr, &
             "exponential pond drainage rate",                       &
             "none", c1, c0,                                       &
             ns, f_dpnd_expon)

      if (f_dpnd_freebd(1:1) /= 'x') &
         call define_hist_field(n_dpnd_freebd,"dpnd_freebd","m/s",tstr2D, tcstr, &
             "pond drainage rate due to freeboard constraint",       &
             "none", c1, c0,                                       &
             ns, f_dpnd_freebd)

      if (f_dpnd_initial(1:1) /= 'x') &
         call define_hist_field(n_dpnd_initial,"dpnd_initial","m/s",tstr2D, tcstr, &
             "runoff rate due to rfrac",                             &
             "none", c1, c0,                                       &
             ns, f_dpnd_initial)

      if (f_dpnd_dlid(1:1) /= 'x') &
         call define_hist_field(n_dpnd_dlid,"dpnd_dlid","m/s",tstr2D, tcstr, &
             "pond loss / gain to ice lid freezing / melting",         &
             "none", c1, c0,                                       &
             ns, f_dpnd_dlid)

      if (f_dpnd_melt(1:1) /= 'x') &
         call define_hist_field(n_dpnd_melt,"dpnd_melt","m/s",tstr2D, tcstr, &
             "pond drainage due to ice melting",                     &
             "none", c1, c0,                                       &
             ns, f_dpnd_melt)

      if (f_dpnd_ridge(1:1) /= 'x') &
         call define_hist_field(n_dpnd_ridge,"dpnd_ridge","m",tstr2D, tcstr, &
             "pond drainage due to ridging",                       &
             "none", c1, c0,                                       &
             ns, f_dpnd_ridge)

      endif ! histfreq(ns) /= 'x'
      enddo ! nstreams

      endif ! tr_pond

      end subroutine init_hist_pond_2D

!=======================================================================

      subroutine init_hist_pond_3Dc

      use ice_calendar, only: nstreams, histfreq
      use ice_history_shared, only: tstr3Dc, tcstr, define_hist_field

      integer (kind=int_kind) :: ns
      logical (kind=log_kind) :: tr_pond
      character(len=*), parameter :: subname = '(init_hist_pond_3Dc)'

      call icepack_query_tracer_flags(tr_pond_out=tr_pond)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (tr_pond) then

      ! 3D (category) variables must be looped separately
      do ns = 1, nstreams
        if (histfreq(ns) /= 'x') then

        if (f_apondn(1:1) /= 'x') &
           call define_hist_field(n_apondn,"apondn","1",tstr3Dc, tcstr, &
              "melt pond fraction, category","none", c1, c0,      &
              ns, f_apondn)

        if (f_hpondn(1:1) /= 'x') &
           call define_hist_field(n_hpondn,"hpondn","m",tstr3Dc, tcstr, &
              "melt pond depth, category","none", c1, c0,       &
              ns, f_hpondn)

        if (f_apeffn(1:1) /= 'x') &
           call define_hist_field(n_apeffn,"apeffn","1",tstr3Dc, tcstr, &
             "effective melt pond fraction, category",   &
             "none", c1, c0,                                  &
             ns, f_apeffn)

        if (f_dpnd_flushn(1:1) /= 'x') &
           call define_hist_field(n_dpnd_flushn,"dpnd_flushn","m/s",tstr3Dc, tcstr, &
               "category pond flushing rate due to ice permeability",     &
               "none", c1, c0,                                       &
               ns, f_dpnd_flushn)

        if (f_dpnd_exponn(1:1) /= 'x') &
           call define_hist_field(n_dpnd_exponn,"dpnd_exponn","m/s",tstr3Dc, tcstr, &
               "category exponential pond drainage rate",                 &
               "none", c1, c0,                                       &
               ns, f_dpnd_exponn)

        if (f_dpnd_freebdn(1:1) /= 'x') &
           call define_hist_field(n_dpnd_freebdn,"dpnd_freebdn","m/s",tstr3Dc, tcstr, &
               "category pond drainage rate due to freeboard constraint", &
               "none", c1, c0,                                       &
               ns, f_dpnd_freebdn)

        if (f_dpnd_initialn(1:1) /= 'x') &
           call define_hist_field(n_dpnd_initialn,"dpnd_initialn","m/s",tstr3Dc, tcstr, &
               "category runoff rate due to rfrac",                       &
               "none", c1, c0,                                       &
               ns, f_dpnd_initialn)

        if (f_dpnd_dlidn(1:1) /= 'x') &
           call define_hist_field(n_dpnd_dlidn,"dpnd_dlidn","m/s",tstr3Dc, tcstr, &
               "category pond loss / gain to ice lid freezing / melting",   &
               "none", c1, c0,                                       &
               ns, f_dpnd_dlidn)

        endif ! histfreq(ns) /= 'x'
      enddo ! ns

      endif ! tr_pond

      end subroutine init_hist_pond_3Dc

!=======================================================================

! accumulate average ice quantities or snapshots

      subroutine accum_hist_pond (iblk)

      use ice_arrays_column, only: apeffn
      use ice_blocks, only: block, get_block, nx_block, ny_block
      use ice_domain, only: blocks_ice
      use ice_flux, only: apeff_ai
      use ice_flux, only: dpnd_flush, dpnd_expon, dpnd_freebd, dpnd_initial
      use ice_flux, only: dpnd_dlid, dpnd_melt,dpnd_ridge, dpnd_dlidn
      use ice_flux, only: dpnd_flushn, dpnd_exponn, dpnd_freebdn, dpnd_initialn
      use ice_history_shared, only: n2D, a2D, a3Dc, ncat_hist, &
          accum_hist_field
      use ice_state, only: aice, trcr, trcrn

      integer (kind=int_kind), intent(in) :: &
           iblk                 ! block index

      ! local variables

      integer (kind=int_kind) :: &
           i,j, &
           ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         worka

      integer (kind=int_kind) :: &
         nt_apnd, nt_hpnd, nt_alvl, nt_ipnd
      logical (kind=log_kind) :: &
         tr_pond_lvl, tr_pond_sealvl, tr_pond_topo

      real (kind=dbl_kind) :: &
         puny

      type (block) :: &
         this_block           ! block information for current block

      character(len=*), parameter :: subname = '(accum_hist_pond)'

      !---------------------------------------------------------------
      ! increment field
      !---------------------------------------------------------------

         call icepack_query_parameters(puny_out=puny)
         call icepack_query_tracer_flags(tr_pond_lvl_out=tr_pond_lvl, &
              tr_pond_topo_out=tr_pond_topo, tr_pond_sealvl_out=tr_pond_sealvl)
         call icepack_query_tracer_indices(nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, &
              nt_alvl_out=nt_alvl, nt_ipnd_out=nt_ipnd)
         call icepack_warnings_flush(nu_diag)
         if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
            file=__FILE__, line=__LINE__)

         if (allocated(a2D)) then

         if (tr_pond_lvl) then

         if (f_apond(1:1)/= 'x') &
             call accum_hist_field(n_apond, iblk, &
                            trcr(:,:,nt_alvl,iblk) * trcr(:,:,nt_apnd,iblk), a2D)
         if (f_apond_ai(1:1)/= 'x') &
             call accum_hist_field(n_apond_ai, iblk, &
                            aice(:,:,iblk) &
                          * trcr(:,:,nt_alvl,iblk) * trcr(:,:,nt_apnd,iblk), a2D)
         if (f_hpond(1:1)/= 'x') &
             call accum_hist_field(n_hpond, iblk, &
                            trcr(:,:,nt_alvl,iblk) * trcr(:,:,nt_apnd,iblk) &
                                                   * trcr(:,:,nt_hpnd,iblk), a2D)
         if (f_hpond_ai(1:1)/= 'x') &
             call accum_hist_field(n_hpond_ai, iblk, &
                            aice(:,:,iblk) &
                          * trcr(:,:,nt_alvl,iblk) * trcr(:,:,nt_apnd,iblk) &
                                                   * trcr(:,:,nt_hpnd,iblk), a2D)
         if (f_ipond(1:1)/= 'x') &
             call accum_hist_field(n_ipond, iblk, &
                            trcr(:,:,nt_alvl,iblk) * trcr(:,:,nt_apnd,iblk) &
                                                   * trcr(:,:,nt_ipnd,iblk), a2D)
         if (f_ipond_ai(1:1)/= 'x') &
             call accum_hist_field(n_ipond_ai, iblk, &
                            aice(:,:,iblk) &
                          * trcr(:,:,nt_alvl,iblk) * trcr(:,:,nt_apnd,iblk) &
                                                   * trcr(:,:,nt_ipnd,iblk), a2D)

         elseif (tr_pond_topo .or. tr_pond_sealvl) then

         if (f_apond(1:1)/= 'x') &
             call accum_hist_field(n_apond, iblk, &
                                   trcr(:,:,nt_apnd,iblk), a2D)
         if (f_apond_ai(1:1)/= 'x') &
             call accum_hist_field(n_apond_ai, iblk, &
                                   aice(:,:,iblk) * trcr(:,:,nt_apnd,iblk), a2D)
         if (f_hpond(1:1)/= 'x') &
             call accum_hist_field(n_hpond, iblk, &
                                   trcr(:,:,nt_apnd,iblk) &
                                 * trcr(:,:,nt_hpnd,iblk), a2D)
         if (f_hpond_ai(1:1)/= 'x') &
             call accum_hist_field(n_hpond_ai, iblk, &
                                   aice(:,:,iblk) * trcr(:,:,nt_apnd,iblk) &
                                                  * trcr(:,:,nt_hpnd,iblk), a2D)
         if (f_ipond(1:1)/= 'x') &
             call accum_hist_field(n_ipond, iblk, &
                                   trcr(:,:,nt_apnd,iblk) &
                                 * trcr(:,:,nt_ipnd,iblk), a2D)
         if (f_ipond_ai(1:1)/= 'x') &
             call accum_hist_field(n_ipond_ai, iblk, &
                                   aice(:,:,iblk) * trcr(:,:,nt_apnd,iblk) &
                                                  * trcr(:,:,nt_ipnd,iblk), a2D)
         endif ! ponds

         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         if (f_apeff (1:1) /= 'x') then
             worka(:,:) = c0
             do j = jlo, jhi
             do i = ilo, ihi
                if (aice(i,j,iblk) > puny) worka(i,j) = apeff_ai(i,j,iblk) &
                                                      / aice(i,j,iblk)
             enddo
             enddo
             call accum_hist_field(n_apeff, iblk, worka(:,:), a2D)
         endif
         if (f_apeff_ai(1:1) /= 'x') &
             call accum_hist_field(n_apeff_ai, iblk, apeff_ai(:,:,iblk), a2D)

         if (f_dpnd_flush (1:1) /= 'x') &
             call accum_hist_field(n_dpnd_flush, iblk, dpnd_flush(:,:,iblk), a2D)
         if (f_dpnd_expon (1:1) /= 'x') &
             call accum_hist_field(n_dpnd_expon, iblk, dpnd_expon(:,:,iblk), a2D)
         if (f_dpnd_freebd (1:1) /= 'x') &
             call accum_hist_field(n_dpnd_freebd, iblk, dpnd_freebd(:,:,iblk), a2D)
         if (f_dpnd_initial (1:1) /= 'x') &
             call accum_hist_field(n_dpnd_initial, iblk, dpnd_initial(:,:,iblk), a2D)
         if (f_dpnd_dlid (1:1) /= 'x') &
             call accum_hist_field(n_dpnd_dlid, iblk, dpnd_dlid(:,:,iblk), a2D)
         if (f_dpnd_melt (1:1) /= 'x') &
             call accum_hist_field(n_dpnd_melt, iblk, dpnd_melt(:,:,iblk), a2D)
         if (f_dpnd_ridge (1:1) /= 'x') &
             call accum_hist_field(n_dpnd_ridge, iblk, dpnd_ridge(:,:,iblk), a2D)

         endif ! allocated(a2D)

         ! 3D category fields
         if (allocated(a3Dc)) then
         if (f_apondn   (1:1) /= 'x') &
             call accum_hist_field(n_apondn-n2D, iblk, ncat_hist, &
                  trcrn(:,:,nt_apnd,1:ncat_hist,iblk), a3Dc)
         if (f_apeffn (1:1) /= 'x') &
             call accum_hist_field(n_apeffn-n2D,  iblk, ncat_hist, &
                  apeffn(:,:,1:ncat_hist,iblk), a3Dc)
         if (f_hpondn   (1:1) /= 'x') &
             call accum_hist_field(n_hpondn-n2D, iblk, ncat_hist, &
                    trcrn(:,:,nt_apnd,1:ncat_hist,iblk) &
                  * trcrn(:,:,nt_hpnd,1:ncat_hist,iblk), a3Dc)

         if (f_dpnd_flushn (1:1) /= 'x') &
             call accum_hist_field(n_dpnd_flushn-n2D, iblk, ncat_hist, dpnd_flushn(:,:,:,iblk), a3Dc)
         if (f_dpnd_exponn (1:1) /= 'x') &
             call accum_hist_field(n_dpnd_exponn-n2D, iblk, ncat_hist, dpnd_exponn(:,:,:,iblk), a3Dc)
         if (f_dpnd_freebdn (1:1) /= 'x') &
             call accum_hist_field(n_dpnd_freebdn-n2D, iblk, ncat_hist, dpnd_freebdn(:,:,:,iblk), a3Dc)
         if (f_dpnd_initialn (1:1) /= 'x') &
             call accum_hist_field(n_dpnd_initialn-n2D, iblk, ncat_hist, dpnd_initialn(:,:,:,iblk), a3Dc)
         if (f_dpnd_dlidn (1:1) /= 'x') &
             call accum_hist_field(n_dpnd_dlidn-n2D, iblk, ncat_hist, dpnd_dlidn(:,:,:,iblk), a3Dc)

         endif ! allocated(a3Dc)

      end subroutine accum_hist_pond

!=======================================================================

      end module ice_history_pond

!=======================================================================
