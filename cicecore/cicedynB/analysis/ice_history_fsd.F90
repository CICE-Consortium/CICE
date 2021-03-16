!=======================================================================

! Floe size distribution history output
!
! 2019 Elizabeth Hunke

      module ice_history_fsd

      use ice_kinds_mod
      use ice_domain_size, only: max_nstrm, nfsd
      use ice_constants, only: c0, c1, c8, c100, mps_to_cmpdy
      use ice_fileunits, only: nu_nml, nml_filename, &
          get_fileunit, release_fileunit
      use ice_fileunits, only: nu_diag
      use ice_exit, only: abort_ice
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_parameters, &
          icepack_query_tracer_flags, icepack_query_tracer_indices

      implicit none
      private
      public :: accum_hist_fsd, init_hist_fsd_2D, init_hist_fsd_3Df, &
                init_hist_fsd_4Df
      
      !---------------------------------------------------------------
      ! flags: write to output file if true or histfreq value
      !---------------------------------------------------------------

      character (len=max_nstrm), public :: &
           f_afsd       = 'm', f_afsdn       = 'm', &
           f_dafsd_newi = 'm', f_dafsd_latg  = 'm', &
           f_dafsd_latm = 'm', f_dafsd_wave  = 'm', &
           f_dafsd_weld = 'm', f_wave_sig_ht = 'm', &
           f_aice_ww    = 'x', f_diam_ww     = 'x', &
           f_hice_ww    = 'x', f_fsdrad      = 'x', &
           f_fsdperim   = 'x'

      !---------------------------------------------------------------
      ! namelist variables
      !---------------------------------------------------------------

      namelist / icefields_fsd_nml /    &
           f_afsd      , f_afsdn      , &
           f_dafsd_newi, f_dafsd_latg , &
           f_dafsd_latm, f_dafsd_wave , &
           f_dafsd_weld, f_wave_sig_ht, &
           f_aice_ww   , f_diam_ww    , &
           f_hice_ww   , f_fsdrad     , &
           f_fsdperim

      !---------------------------------------------------------------
      ! field indices
      !---------------------------------------------------------------

      integer (kind=int_kind), dimension(max_nstrm) :: &
           n_afsd      , n_afsdn      , &
           n_dafsd_newi, n_dafsd_latg , &
           n_dafsd_latm, n_dafsd_wave , &
           n_dafsd_weld, n_wave_sig_ht, &
           n_aice_ww   , n_diam_ww    , &
           n_hice_ww   , n_fsdrad     , &
           n_fsdperim

!=======================================================================

      contains

!=======================================================================

! Initialize history files
! authors Elizabeth C. Hunke, LANL

      subroutine init_hist_fsd_2D

      use ice_broadcast, only: broadcast_scalar
      use ice_calendar, only: nstreams
      use ice_communicate, only: my_task, master_task
      use ice_history_shared, only: tstr2D, tcstr, define_hist_field

      integer (kind=int_kind) :: ns
      integer (kind=int_kind) :: nml_error ! namelist i/o error flag
      real    (kind=dbl_kind) :: secday
      logical (kind=log_kind) :: tr_fsd, wave_spec
      character(len=*), parameter :: subname = '(init_hist_fsd_2D)'

      call icepack_query_tracer_flags(tr_fsd_out=tr_fsd)
      call icepack_query_parameters(wave_spec_out=wave_spec)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (tr_fsd) then

      !-----------------------------------------------------------------
      ! read namelist
      !-----------------------------------------------------------------

      call get_fileunit(nu_nml)
      if (my_task == master_task) then
         open (nu_nml, file=nml_filename, status='old',iostat=nml_error)
         if (nml_error /= 0) then
            nml_error = -1
         else
            nml_error =  1
         endif
         do while (nml_error > 0)
            read(nu_nml, nml=icefields_fsd_nml,iostat=nml_error)
         end do
         if (nml_error == 0) close(nu_nml)
      endif
      call release_fileunit(nu_nml)

      call broadcast_scalar(nml_error, master_task)
      if (nml_error /= 0) then
         close (nu_nml)
         call abort_ice(subname//'ERROR: reading icefields_fsd_nml')
      endif

      call broadcast_scalar (f_afsd, master_task)
      call broadcast_scalar (f_afsdn, master_task)
      call broadcast_scalar (f_dafsd_newi, master_task)
      call broadcast_scalar (f_dafsd_latg, master_task)
      call broadcast_scalar (f_dafsd_latm, master_task)
      call broadcast_scalar (f_dafsd_wave, master_task)
      call broadcast_scalar (f_dafsd_weld, master_task)
      call broadcast_scalar (f_wave_sig_ht, master_task)
      call broadcast_scalar (f_aice_ww, master_task)
      call broadcast_scalar (f_diam_ww, master_task)
      call broadcast_scalar (f_hice_ww, master_task)
      call broadcast_scalar (f_fsdrad, master_task)
      call broadcast_scalar (f_fsdperim, master_task)

      ! 2D variables

      do ns = 1, nstreams

      if (f_wave_sig_ht(1:1) /= 'x') &
         call define_hist_field(n_wave_sig_ht,"wave_sig_ht","1",tstr2D, tcstr, &
             "significant height of wind and swell waves",  &
             "from attenuated spectrum in ice", c1, c0,     &
             ns, f_wave_sig_ht)
      if (f_aice_ww(1:1) /= 'x') &
         call define_hist_field(n_aice_ww,"aice_ww","1",tstr2D, tcstr, &
             "Concentration of floes > Dmin",               &
             "for waves", c1, c0,                           &
             ns, f_aice_ww)
      if (f_diam_ww(1:1) /= 'x') &
         call define_hist_field(n_diam_ww,"diam_ww","1",tstr2D, tcstr, &
             "Average (number) diameter of floes > Dmin",   &
             "for waves", c1, c0,                           &
             ns, f_diam_ww)
      if (f_hice_ww(1:1) /= 'x') &
         call define_hist_field(n_hice_ww,"hice_ww","m",tstr2D, tcstr, &
             "Thickness of floes > Dmin",                   &
             "for waves", c1, c0,                           &
             ns, f_hice_ww)
      if (f_fsdrad(1:1) /= 'x') &
         call define_hist_field(n_fsdrad,"fsdrad","m",tstr2D, tcstr, &
            "floe size distribution, representative radius",                  &
            " ", c1, c0, ns, f_fsdrad)
      if (f_fsdperim(1:1) /= 'x') &
         call define_hist_field(n_fsdperim,"fsdperim","1/m",tstr2D, tcstr, &
            "floe size distribution, perimeter",                  &
            "per unit ice area", c1, c0, ns, f_fsdperim)


      enddo ! nstreams

      else  ! tr_fsd

         f_afsd        = 'x'
         f_afsdn       = 'x'
         f_dafsd_newi  = 'x'
         f_dafsd_latg  = 'x'
         f_dafsd_latm  = 'x'
         f_dafsd_wave  = 'x'
         f_dafsd_weld  = 'x'
         f_wave_sig_ht = 'x'
         f_fsdrad      = 'x'
         f_fsdperim    = 'x'
         if (.not. wave_spec) then
            f_aice_ww  = 'x'
            f_diam_ww  = 'x'
            f_hice_ww  = 'x'
         endif

      endif ! tr_fsd

      end subroutine init_hist_fsd_2D

!=======================================================================

      subroutine init_hist_fsd_3Df

      use ice_calendar, only: nstreams, histfreq
      use ice_history_shared, only: tstr3Df, tcstr, define_hist_field

      logical (kind=log_kind) :: tr_fsd
      integer (kind=int_kind) :: ns
      character(len=*), parameter :: subname = '(init_hist_fsd_3Df)'

      call icepack_query_tracer_flags(tr_fsd_out=tr_fsd)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      !-----------------------------------------------------------------
      ! 3D (category) variables must be looped separately
      !-----------------------------------------------------------------

      if (tr_fsd) then

      do ns = 1, nstreams
         if (histfreq(ns) /= 'x') then

         if (f_afsd(1:1) /= 'x') &
            call define_hist_field(n_afsd,"afsd", "1", tstr3Df, tcstr, &
               "areal floe size distribution",                 &
               "per unit bin width ", c1, c0, ns, f_afsd)
         if (f_dafsd_newi(1:1) /= 'x') &
            call define_hist_field(n_dafsd_newi,"dafsd_newi","1",tstr3Df, tcstr, &
               "Change in fsd: new ice",                       &
               "Avg over freq period", c1, c0, ns, f_dafsd_newi)
         if (f_dafsd_latg(1:1) /= 'x') &
            call define_hist_field(n_dafsd_latg,"dafsd_latg","1",tstr3Df, tcstr, &
               "Change in fsd: lateral growth",                &
               "Avg over freq period", c1, c0, ns, f_dafsd_latg)
         if (f_dafsd_latm(1:1) /= 'x') &
            call define_hist_field(n_dafsd_latm,"dafsd_latm","1",tstr3Df, tcstr, &
               "Change in fsd: lateral melt",                  &
               "Avg over freq period", c1, c0, ns, f_dafsd_latm)
         if (f_dafsd_wave(1:1) /= 'x') &
            call define_hist_field(n_dafsd_wave,"dafsd_wave","1",tstr3Df, tcstr, &
               "Change in fsd: waves",                         &
               "Avg over freq period", c1, c0, ns, f_dafsd_wave)
         if (f_dafsd_weld(1:1) /= 'x') &
            call define_hist_field(n_dafsd_weld,"dafsd_weld","1",tstr3Df, tcstr, &
               "Change in fsd: welding",                       &
               "Avg over freq period", c1, c0, ns, f_dafsd_weld)
         endif ! if (histfreq(ns) /= 'x')
      enddo ! ns

      endif ! tr_fsd

      end subroutine init_hist_fsd_3Df

!=======================================================================

      subroutine init_hist_fsd_4Df

      use ice_calendar, only: nstreams, histfreq
      use ice_history_shared, only: tstr4Df, tcstr, define_hist_field

      logical (kind=log_kind) :: tr_fsd
      integer (kind=int_kind) :: ns
      character(len=*), parameter :: subname = '(init_hist_fsd_4Df)'

      call icepack_query_tracer_flags(tr_fsd_out=tr_fsd)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      !-----------------------------------------------------------------
      ! 4D (floe size, thickness category) variables
      !-----------------------------------------------------------------

      if (tr_fsd) then

      do ns = 1, nstreams
         if (histfreq(ns) /= 'x') then

         if (f_afsdn(1:1) /= 'x') &
            call define_hist_field(n_afsdn,"afsdn","1",tstr4Df, tcstr, & 
               "areal floe size and thickness distribution",    &
               "per unit bin width", c1, c0, ns, f_afsdn)

         endif ! if (histfreq(ns) /= 'x') then
      enddo ! ns 

      endif ! tr_fsd

      end subroutine init_hist_fsd_4Df

!=======================================================================

! accumulate average ice quantities or snapshots
! author:   Elizabeth C. Hunke, LANL

      subroutine accum_hist_fsd (iblk)

      use ice_blocks, only: nx_block, ny_block
      use ice_constants, only: c0, c1, c2, c4
      use ice_history_shared, only: a2D, a3Df, a4Df, nfsd_hist, &
         ncat_hist, accum_hist_field, n3Dacum, n4Dscum
      use ice_state, only: trcrn, aicen_init, vicen, aice_init
      use ice_arrays_column, only: wave_sig_ht, floe_rad_c, floe_binwidth, &
         d_afsd_newi, d_afsd_latg, d_afsd_latm, d_afsd_wave, d_afsd_weld

      integer (kind=int_kind), intent(in) :: &
           iblk                 ! block index

      ! local variables

      integer (kind=int_kind) :: &
         i, j, n, k, & ! loop indices
         nt_fsd        ! ! fsd tracer index
      logical (kind=log_kind) :: tr_fsd
      real (kind=dbl_kind) :: floeshape, puny

      real (kind=dbl_kind) :: workb, workc
      real (kind=dbl_kind), dimension(nx_block,ny_block) :: worka
      real (kind=dbl_kind), dimension(nx_block,ny_block,nfsd_hist) :: worke
      real (kind=dbl_kind), dimension(nx_block,ny_block,nfsd_hist,ncat_hist) :: workd

      character(len=*), parameter :: subname = '(accum_hist_fsd)'

      call icepack_query_parameters(floeshape_out=floeshape, puny_out=puny)
      call icepack_query_tracer_flags(tr_fsd_out=tr_fsd)
      call icepack_query_tracer_indices(nt_fsd_out=nt_fsd)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (tr_fsd) then

      !---------------------------------------------------------------
      ! increment field
      !---------------------------------------------------------------

      ! 2D fields
      if (allocated(a2D)) then

      if (f_wave_sig_ht(1:1)/= 'x') &
         call accum_hist_field(n_wave_sig_ht,   iblk, &
                               wave_sig_ht(:,:,iblk), a2D)

      if (f_aice_ww(1:1)/= 'x') then
         do j = 1, ny_block
         do i = 1, nx_block
            worka(i,j) = c0
            do n = 1, ncat_hist
            do k = 1, nfsd_hist
               worka(i,j) = worka(i,j) + aicen_init(i,j,n,iblk)*trcrn(i,j,nt_fsd+k-1,n,iblk)
            end do
            end do
         end do
         end do
         call accum_hist_field(n_aice_ww,   iblk, worka(:,:), a2D)
      endif

      if (f_diam_ww (1:1) /= 'x') then
         do j = 1, ny_block
         do i = 1, nx_block
            worka(i,j) = c0
            workb      = c0
            do n = 1, ncat_hist
            do k = 1, nfsd_hist
               workc = aicen_init(i,j,n,iblk)*trcrn(i,j,nt_fsd+k-1,n,iblk) &
                     / (c4*floeshape*floe_rad_c(k)**2)
               ! number-mean radius
               worka(i,j) = worka(i,j) + workc * floe_rad_c(k)
               ! normalization factor
               workb      = workb      + workc
            end do
            end do
            ! number-mean diameter, following Eqn. 5 in HT2017
            if (workb > c0) worka(i,j) = c2*worka(i,j) / workb
            worka(i,j) = MAX(c2*floe_rad_c(1),worka(i,j)) ! >= smallest resolved diameter
         end do
         end do
         call accum_hist_field(n_diam_ww, iblk, worka(:,:), a2D)
      endif

      if (f_hice_ww (1:1) /= 'x') then
         do j = 1, ny_block
         do i = 1, nx_block
            worka(i,j) = c0
            workb      = c0
            do n = 1, ncat_hist
            do k = 1, nfsd_hist
               workb = workb + aicen_init(i,j,n,iblk)*trcrn(i,j,nt_fsd+k-1,n,iblk)
               worka(i,j) = worka(i,j) + vicen(i,j,n,iblk)*trcrn(i,j,nt_fsd+k-1,n,iblk)
            end do
            end do
            if (workb > puny) then
               worka(i,j) = worka(i,j) / workb
            else
               worka(i,j) = c1
            endif
         end do
         end do
         call accum_hist_field(n_hice_ww, iblk, worka(:,:), a2D)
      endif

      if (f_fsdrad(1:1) /= 'x') then
         do j = 1, ny_block
         do i = 1, nx_block
            worka(i,j) = c0            
            if (aice_init(i,j,iblk) > puny) then
             do k = 1, nfsd_hist
                do n = 1, ncat_hist
                  worka(i,j) = worka(i,j) &
                               + (trcrn(i,j,nt_fsd+k-1,n,iblk) * floe_rad_c(k) &
                               * aicen_init(i,j,n,iblk)/aice_init(i,j,iblk))
                 end do
              end do
            endif
         end do
         end do
         call accum_hist_field(n_fsdrad, iblk, worka(:,:), a2D)
      endif

      if (f_fsdperim(1:1) /= 'x') then
         do j = 1, ny_block
         do i = 1, nx_block
            worka(i,j) = c0
            if (aice_init(i,j,iblk) > puny) then
             do k = 1, nfsd_hist
               do n = 1, ncat_hist
                  worka(i,j) = worka(i,j) &
                               + (c8*floeshape*trcrn(i,j,nt_fsd+k-1,n,iblk)*floe_rad_c(k) &
                                    *aicen_init(i,j,n,iblk)/(c4*floeshape*floe_rad_c(k)**2 *aice_init(i,j,iblk)))
               end do
              end do
            endif
         end do
         end do
         call accum_hist_field(n_fsdperim, iblk, worka, a2D)
      endif

      endif ! a2D allocated

      ! 3D category fields
      if (allocated(a3Df)) then

      if (f_afsd(1:1) /= 'x') then
         do j = 1, ny_block
         do i = 1, nx_block
            do k = 1, nfsd_hist
               worke(i,j,k)=c0
               do n = 1, ncat_hist
                  worke(i,j,k) = worke(i,j,k) + (trcrn(i,j,nt_fsd+k-1,n,iblk) &
                               * aicen_init(i,j,n,iblk)/floe_binwidth(k))
               end do
            end do
         end do
         end do
         call accum_hist_field(n_afsd-n3Dacum, iblk, nfsd_hist, worke, a3Df)
      endif
 
      if (f_dafsd_newi(1:1)/= 'x') &
             call accum_hist_field(n_dafsd_newi-n3Dacum, iblk, nfsd_hist, &
                                    d_afsd_newi(:,:,1:nfsd_hist,iblk), a3Df)
      if (f_dafsd_latg(1:1)/= 'x') &
             call accum_hist_field(n_dafsd_latg-n3Dacum, iblk, nfsd_hist, &
                                    d_afsd_latg(:,:,1:nfsd_hist,iblk), a3Df)
      if (f_dafsd_latm(1:1)/= 'x') &
             call accum_hist_field(n_dafsd_latm-n3Dacum, iblk, nfsd_hist, &
                                    d_afsd_latm(:,:,1:nfsd_hist,iblk), a3Df)
      if (f_dafsd_wave(1:1)/= 'x') &
             call accum_hist_field(n_dafsd_wave-n3Dacum, iblk, nfsd_hist, &
                                    d_afsd_wave(:,:,1:nfsd_hist,iblk), a3Df)
      if (f_dafsd_weld(1:1)/= 'x') &
             call accum_hist_field(n_dafsd_weld-n3Dacum, iblk, nfsd_hist, &
                                    d_afsd_weld(:,:,1:nfsd_hist,iblk), a3Df)
      endif ! a3Df allocated

      ! 4D floe size, thickness category fields
      if (allocated(a4Df)) then

      if (f_afsdn(1:1) /= 'x') then
         do n = 1, ncat_hist
         do k = 1, nfsd_hist 
         do j = 1, ny_block
         do i = 1, nx_block
            workd(i,j,k,n) = trcrn(i,j,nt_fsd+k-1,n,iblk) &
                           * aicen_init(i,j,n,iblk)/floe_binwidth(k)
         end do
         end do
         end do
         end do
         call accum_hist_field(n_afsdn-n4Dscum, iblk, nfsd_hist, ncat_hist, &
                               workd(:,:,1:nfsd_hist,1:ncat_hist), a4Df)
      endif

      endif ! a4Df allocated

      endif ! tr_fsd

      end subroutine accum_hist_fsd

!=======================================================================

      end module ice_history_fsd

!=======================================================================
