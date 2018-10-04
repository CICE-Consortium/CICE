!=======================================================================

! Mechanical redistribution history output
!
! 2012 Elizabeth Hunke split code from ice_history.F90

      module ice_history_mechred

      use ice_kinds_mod
      use ice_domain_size, only: max_nstrm
      use ice_constants, only: c0, c1, c100, mps_to_cmpdy
      use ice_fileunits, only: nu_nml, nml_filename, &
          get_fileunit, release_fileunit
      use ice_fileunits, only: nu_diag
      use ice_exit, only: abort_ice
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_parameters, &
          icepack_query_tracer_flags, icepack_query_tracer_indices

      implicit none
      private
      public :: accum_hist_mechred, init_hist_mechred_2D, init_hist_mechred_3Dc
      
      !---------------------------------------------------------------
      ! flags: write to output file if true or histfreq value
      !---------------------------------------------------------------

      character (len=max_nstrm), public :: &
           f_ardg      = 'm', f_vrdg       = 'm', &
           f_alvl      = 'm', f_vlvl       = 'm', &
           f_dardg1dt  = 'm', f_dardg2dt   = 'm', &
           f_dvirdgdt  = 'm', f_opening    = 'm', &
           f_ardgn     = 'x', f_vrdgn      = 'x', &
           f_dardg1ndt = 'x', f_dardg2ndt  = 'x', &
           f_dvirdgndt = 'x', &
           f_aparticn  = 'x', f_krdgn      = 'x', &
           f_aredistn  = 'x', f_vredistn   = 'x', &
           f_araftn    = 'x', f_vraftn     = 'x'

      !---------------------------------------------------------------
      ! namelist variables
      !---------------------------------------------------------------

      namelist / icefields_mechred_nml /     &
           f_ardg,      f_vrdg     , &
           f_alvl,      f_vlvl     , &
           f_dardg1dt,  f_dardg2dt , &
           f_dvirdgdt,  f_opening  , &
           f_ardgn,     f_vrdgn    , &
           f_dardg1ndt, f_dardg2ndt, &
           f_dvirdgndt, &
           f_aparticn,  f_krdgn    , &
           f_aredistn,  f_vredistn , &
           f_araftn,    f_vraftn

      !---------------------------------------------------------------
      ! field indices
      !---------------------------------------------------------------

      integer (kind=int_kind), dimension(max_nstrm) :: &
           n_ardg       , n_vrdg       , &
           n_alvl       , n_vlvl       , &
           n_dardg1dt   , n_dardg2dt   , &
           n_dvirdgdt   , n_opening    , &
           n_ardgn      , n_vrdgn      , &
           n_dardg1ndt  , n_dardg2ndt  , &
           n_dvirdgndt  , &
           n_aparticn   , n_krdgn      , &
           n_aredistn   , n_vredistn   , &
           n_araftn     , n_vraftn

!=======================================================================

      contains

!=======================================================================

! Initialize history files
! authors Elizabeth C. Hunke, LANL

      subroutine init_hist_mechred_2D

      use ice_broadcast, only: broadcast_scalar
      use ice_calendar, only: nstreams
      use ice_communicate, only: my_task, master_task
      use ice_history_shared, only: tstr2D, tcstr, define_hist_field

      integer (kind=int_kind) :: ns
      integer (kind=int_kind) :: nml_error ! namelist i/o error flag
      real    (kind=dbl_kind) :: secday
      logical (kind=log_kind) :: tr_lvl
      character(len=*), parameter :: subname = '(init_hist_mechred_2D)'

      call icepack_query_parameters(secday_out=secday)
      call icepack_query_tracer_flags(tr_lvl_out=tr_lvl)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

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
            read(nu_nml, nml=icefields_mechred_nml,iostat=nml_error)
         end do
         if (nml_error == 0) close(nu_nml)
      endif
      call release_fileunit(nu_nml)

      call broadcast_scalar(nml_error, master_task)
      if (nml_error /= 0) then
         close (nu_nml)
         call abort_ice(subname//'ERROR: reading icefields_mechred_nml')
      endif

      if (.not. tr_lvl) then
         f_ardg = 'x'
         f_vrdg = 'x'
         f_alvl = 'x'
         f_vlvl = 'x'
         f_ardgn = 'x'
         f_vrdgn = 'x'
         f_araftn = 'x'
         f_vraftn = 'x'
      endif
      if (f_araftn /= 'x' .or. f_vraftn /= 'x') f_ardgn = f_araftn

      call broadcast_scalar (f_ardg, master_task)
      call broadcast_scalar (f_vrdg, master_task)
      call broadcast_scalar (f_alvl, master_task)
      call broadcast_scalar (f_vlvl, master_task)
      call broadcast_scalar (f_dardg1dt, master_task)
      call broadcast_scalar (f_dardg2dt, master_task)
      call broadcast_scalar (f_dvirdgdt, master_task)
      call broadcast_scalar (f_opening, master_task)
      call broadcast_scalar (f_ardgn, master_task)
      call broadcast_scalar (f_vrdgn, master_task)
      call broadcast_scalar (f_dardg1ndt, master_task)
      call broadcast_scalar (f_dardg2ndt, master_task)
      call broadcast_scalar (f_dvirdgndt, master_task)
      call broadcast_scalar (f_krdgn, master_task)
      call broadcast_scalar (f_aparticn, master_task)
      call broadcast_scalar (f_aredistn, master_task)
      call broadcast_scalar (f_vredistn, master_task)
      call broadcast_scalar (f_araftn, master_task)
      call broadcast_scalar (f_vraftn, master_task)

      ! 2D variables

      do ns = 1, nstreams

      if (f_alvl(1:1) /= 'x') &
         call define_hist_field(n_alvl,"alvl","1",tstr2D, tcstr, &
             "level ice area fraction",                            &
             "none", c1, c0,                                       &
             ns, f_alvl)
      if (f_vlvl(1:1) /= 'x') &
         call define_hist_field(n_vlvl,"vlvl","m",tstr2D, tcstr, &
             "level ice volume",                           &
             "grid cell mean level ice thickness", c1, c0, &
             ns, f_vlvl)
      if (f_ardg(1:1) /= 'x') &
         call define_hist_field(n_ardg,"ardg","1",tstr2D, tcstr, &
             "ridged ice area fraction",                           &
             "none", c1, c0,                                       &
             ns, f_ardg)
      if (f_vrdg(1:1) /= 'x') &
         call define_hist_field(n_vrdg,"vrdg","m",tstr2D, tcstr, &
             "ridged ice volume",                          &
             "grid cell mean level ridged thickness", c1, c0, &
             ns, f_vrdg)

      if (f_dardg1dt(1:1) /= 'x') &
         call define_hist_field(n_dardg1dt,"dardg1dt","%/day",tstr2D, tcstr, &
             "ice area ridging rate",                                      &
             "none", secday*c100, c0,                                      &
             ns, f_dardg1dt)
      
      if (f_dardg2dt(1:1) /= 'x') &
         call define_hist_field(n_dardg2dt,"dardg2dt","%/day",tstr2D, tcstr, &
             "ridge area formation rate",                                  &
             "none", secday*c100, c0,                                      &
             ns, f_dardg2dt)
      
      if (f_dvirdgdt(1:1) /= 'x') &
         call define_hist_field(n_dvirdgdt,"dvirdgdt","cm/day",tstr2D, tcstr, &
             "ice volume ridging rate",                                     &
             "none", mps_to_cmpdy, c0,                                      &
             ns, f_dvirdgdt)

      if (f_opening(1:1) /= 'x') &
         call define_hist_field(n_opening,"opening","%/day",tstr2D, tcstr, &
             "lead area opening rate",                                   &
             "none", secday*c100, c0,                                    &
             ns, f_opening)

      enddo ! nstreams

      end subroutine init_hist_mechred_2D

!=======================================================================

      subroutine init_hist_mechred_3Dc

      use ice_calendar, only: nstreams
      use ice_history_shared, only: tstr3Dc, tcstr, define_hist_field

      integer (kind=int_kind) :: ns
      real (kind=dbl_kind) :: secday
      character(len=*), parameter :: subname = '(init_hist_mechred_3Dc)'

      !-----------------------------------------------------------------
      ! 3D (category) variables must be looped separately
      !-----------------------------------------------------------------

      call icepack_query_parameters(secday_out=secday)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      do ns = 1, nstreams

       if (f_ardgn(1:1) /= 'x') &
           call define_hist_field(n_ardgn,"ardgn","1",tstr3Dc, tcstr, &
             "ridged ice area fraction, category",                 &
             "none", c1, c0,                                       &
             ns, f_ardgn)

       if (f_vrdgn(1:1) /= 'x') &
           call define_hist_field(n_vrdgn,"vrdgn","m",tstr3Dc, tcstr, &
             "ridged ice volume, category",                &
             "grid cell mean ridged ice thickness", c1, c0, &
             ns, f_vrdgn)

       if (f_dardg1ndt(1:1) /= 'x') &
           call define_hist_field(n_dardg1ndt,"dardg1ndt","%/day",tstr3Dc, tcstr, &
             "ice area ridging rate, category",                            &
             "none", secday*c100, c0,                                      &
             ns, f_dardg1ndt)

       if (f_dardg2ndt(1:1) /= 'x') &
           call define_hist_field(n_dardg2ndt,"dardg2ndt","%/day",tstr3Dc, tcstr, &
             "ridge area formation rate, category",                        &
             "none", secday*c100, c0,                                      &
             ns, f_dardg2ndt)

       if (f_dvirdgndt(1:1) /= 'x') &
          call define_hist_field(n_dvirdgndt,"dvirdgndt","cm/day",tstr3Dc, tcstr, &
             "ice volume ridging rate, category",                          &
             "none", mps_to_cmpdy, c0,                                     &
             ns, f_dvirdgndt)

       if (f_krdgn(1:1) /= 'x') &
           call define_hist_field(n_krdgn,"krdgn","1",tstr3Dc, tcstr, &
             "ridging thickness factor, category",                    &
             "mean ridge thickness/thickness of ridging ice", c1, c0, &
             ns, f_krdgn)

       if (f_aparticn(1:1) /= 'x') &
           call define_hist_field(n_aparticn,"aparticn","1",tstr3Dc, tcstr, &
             "ridging ice participation function, category",       &
             "fraction of new ridge area added to cat", c1, c0,    &
             ns, f_aparticn)

       if (f_aredistn(1:1) /= 'x') &
           call define_hist_field(n_aredistn,"aredistn","1",tstr3Dc, tcstr, &
             "ridging ice area redistribution function, category",   &
             "fraction of new ridge volume added to cat", c1, c0,    &
             ns, f_aredistn)

       if (f_vredistn(1:1) /= 'x') &
           call define_hist_field(n_vredistn,"vredistn","1",tstr3Dc, tcstr, &
             "ridging ice volume redistribution function, category",       &
             "none", c1, c0,                                       &
             ns, f_vredistn)

       if (f_araftn(1:1) /= 'x') &
           call define_hist_field(n_araftn,"araftn","1",tstr3Dc, tcstr, &
             "rafted ice area fraction, category",                 &
             "none", c1, c0,                                       &
             ns, f_araftn)

       if (f_vraftn(1:1) /= 'x') &
           call define_hist_field(n_vraftn,"vraftn","1",tstr3Dc, tcstr, &
             "rafted ice volume, category",                 &
             "none", c1, c0,                                       &
             ns, f_vraftn)

      enddo ! ns

      end subroutine init_hist_mechred_3Dc

!=======================================================================

! accumulate average ice quantities or snapshots
! author:   Elizabeth C. Hunke, LANL

      subroutine accum_hist_mechred (iblk)

      use ice_history_shared, only: n2D, a2D, a3Dc, ncat_hist, &
          accum_hist_field
      use ice_state, only: aice, vice, trcr, aicen, vicen, trcrn
      use ice_flux, only: dardg1dt, dardg2dt, dvirdgdt, dardg1ndt,&
          dardg2ndt, dvirdgndt, krdgn, aparticn, aredistn, vredistn, &
          araftn, vraftn, opening

      integer (kind=int_kind), intent(in) :: &
           iblk                 ! block index

      ! local variables

      integer (kind=int_kind) :: &
           nt_alvl, nt_vlvl
      character(len=*), parameter :: subname = '(accum_hist_mechred)'

      call icepack_query_tracer_indices(nt_alvl_out=nt_alvl, nt_vlvl_out=nt_vlvl)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      !---------------------------------------------------------------
      ! increment field
      !---------------------------------------------------------------

         ! 2D fields

         if (f_alvl(1:1)/= 'x') &
             call accum_hist_field(n_alvl,   iblk, &
                                   aice(:,:,iblk) * trcr(:,:,nt_alvl,iblk), a2D)
         if (f_vlvl(1:1)/= 'x') &
             call accum_hist_field(n_vlvl,   iblk, &
                                   vice(:,:,iblk) * trcr(:,:,nt_vlvl,iblk), a2D)
         if (f_ardg(1:1)/= 'x') &
             call accum_hist_field(n_ardg,   iblk, &
                             aice(:,:,iblk) * (c1 - trcr(:,:,nt_alvl,iblk)), a2D)
         if (f_vrdg(1:1)/= 'x') &
             call accum_hist_field(n_vrdg,   iblk, &
                             vice(:,:,iblk) * (c1 - trcr(:,:,nt_vlvl,iblk)), a2D)
         if (f_dardg1dt(1:1)/= 'x') &
             call accum_hist_field(n_dardg1dt,iblk, dardg1dt(:,:,iblk), a2D)
         if (f_dardg2dt(1:1)/= 'x') &
             call accum_hist_field(n_dardg2dt,iblk, dardg2dt(:,:,iblk), a2D)
         if (f_dvirdgdt(1:1)/= 'x') &
             call accum_hist_field(n_dvirdgdt,iblk, dvirdgdt(:,:,iblk), a2D)
         if (f_opening(1:1) /= 'x') &
             call accum_hist_field(n_opening, iblk, opening(:,:,iblk), a2D)

         ! 3D category fields

         if (f_ardgn(1:1)/= 'x') &
             call accum_hist_field(n_ardgn-n2D, iblk, ncat_hist, &
                                   aicen(:,:,1:ncat_hist,iblk) &
                                 * (c1 - trcrn(:,:,nt_alvl,1:ncat_hist,iblk)), a3Dc)
         if (f_vrdgn(1:1)/= 'x') &
             call accum_hist_field(n_vrdgn-n2D, iblk, ncat_hist, &
                                   vicen(:,:,1:ncat_hist,iblk) &
                                 * (c1 - trcrn(:,:,nt_vlvl,1:ncat_hist,iblk)), a3Dc)
         if (f_dardg1ndt(1:1)/= 'x') &
             call accum_hist_field(n_dardg1ndt-n2D, iblk, ncat_hist, &
                                   dardg1ndt(:,:,1:ncat_hist,iblk), a3Dc)
         if (f_dardg2ndt(1:1)/= 'x') &
             call accum_hist_field(n_dardg2ndt-n2D, iblk, ncat_hist, &
                                   dardg2ndt(:,:,1:ncat_hist,iblk), a3Dc)
         if (f_dvirdgndt(1:1)/= 'x') &
             call accum_hist_field(n_dvirdgndt-n2D, iblk, ncat_hist, &
                                   dvirdgndt(:,:,1:ncat_hist,iblk), a3Dc)
         if (f_krdgn(1:1)/= 'x') &
             call accum_hist_field(n_krdgn-n2D, iblk, ncat_hist, &
                                   krdgn(:,:,1:ncat_hist,iblk), a3Dc)
         if (f_aparticn(1:1)/= 'x') &
             call accum_hist_field(n_aparticn-n2D, iblk, ncat_hist, &
                                   aparticn(:,:,1:ncat_hist,iblk), a3Dc)
         if (f_aredistn(1:1)/= 'x') &
             call accum_hist_field(n_aredistn-n2D, iblk, ncat_hist, &
                                   aredistn(:,:,1:ncat_hist,iblk), a3Dc)
         if (f_vredistn(1:1)/= 'x') &
             call accum_hist_field(n_vredistn-n2D, iblk, ncat_hist, &
                                   vredistn(:,:,1:ncat_hist,iblk), a3Dc)
         if (f_araftn(1:1)/= 'x') &
             call accum_hist_field(n_araftn-n2D, iblk, ncat_hist, &
                                   araftn(:,:,1:ncat_hist,iblk), a3Dc)
         if (f_vraftn(1:1)/= 'x') &
             call accum_hist_field(n_vraftn-n2D, iblk, ncat_hist, &
                                   vraftn(:,:,1:ncat_hist,iblk), a3Dc)

      end subroutine accum_hist_mechred

!=======================================================================

      end module ice_history_mechred

!=======================================================================
