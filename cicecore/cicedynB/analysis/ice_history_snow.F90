!=======================================================================

! Snow tracer history output

      module ice_history_snow

      use ice_kinds_mod
      use ice_constants, only: c0, c1, mps_to_cmpdy
      use ice_domain_size, only: max_nstrm, nslyr
      use ice_fileunits, only: nu_nml, nml_filename, &
          get_fileunit, release_fileunit
      use ice_fileunits, only: nu_diag
      use ice_exit, only: abort_ice
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_parameters, &
          icepack_query_tracer_flags, icepack_query_tracer_indices

      implicit none
      private
      public :: accum_hist_snow, init_hist_snow_2D, init_hist_snow_3Dc
      
      !---------------------------------------------------------------
      ! flags: write to output file if true or histfreq value
      !---------------------------------------------------------------

      character (len=max_nstrm), public :: &
           f_smice     = 'm', f_smicen      = 'x', &
           f_smliq     = 'm', f_smliqn      = 'x', &
           f_rhos      = 'm', f_rhosn       = 'x', &
           f_rsnw      = 'm', f_rsnwn       = 'x', &
           f_meltsliq  = 'm', f_fsloss      = 'x'

      !---------------------------------------------------------------
      ! namelist variables
      !---------------------------------------------------------------

      namelist / icefields_snow_nml / &
           f_smice,     f_smicen, &
           f_smliq,     f_smliqn, &
           f_rhos,      f_rhosn,  &
           f_rsnw,      f_rsnwn,  &
           f_meltsliq,  f_fsloss

      !---------------------------------------------------------------
      ! field indices
      !---------------------------------------------------------------

      integer (kind=int_kind), dimension(max_nstrm), public :: &
           n_smice,     n_smicen, &
           n_smliq,     n_smliqn, &
           n_rhos,      n_rhosn,  &
           n_rsnw,      n_rsnwn,  &
           n_meltsliq,  n_fsloss

!=======================================================================

      contains

!=======================================================================

      subroutine init_hist_snow_2D (dt)

      use ice_broadcast, only: broadcast_scalar
      use ice_calendar, only: nstreams, histfreq
      use ice_communicate, only: my_task, master_task
      use ice_history_shared, only: tstr2D, tcstr, define_hist_field
      use ice_fileunits, only: nu_nml, nml_filename, &
          get_fileunit, release_fileunit

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      integer (kind=int_kind) :: ns
      integer (kind=int_kind) :: nml_error ! namelist i/o error flag
      real    (kind=dbl_kind) :: rhofresh, secday
      logical (kind=log_kind) :: tr_snow
      character(len=*), parameter :: subname = '(init_hist_snow_2D)'

      call icepack_query_tracer_flags(tr_snow_out=tr_snow)
      call icepack_query_parameters(rhofresh_out=rhofresh,secday_out=secday)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (tr_snow) then

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
            read(nu_nml, nml=icefields_snow_nml,iostat=nml_error)
         end do
         if (nml_error == 0) close(nu_nml)
      endif
      call release_fileunit(nu_nml)

      call broadcast_scalar(nml_error, master_task)
      if (nml_error /= 0) then
         close (nu_nml)
         call abort_ice('ice: error reading icefields_snow_nml')
      endif

      else ! .not. tr_snow
          f_smice    = 'x'
          f_smliq    = 'x'
          f_rhos     = 'x'
          f_rsnw     = 'x'
          f_smicen   = 'x'
          f_smliqn   = 'x'
          f_rhosn    = 'x'
          f_rsnwn    = 'x'
          f_meltsliq = 'x'
          f_fsloss   = 'x'
      endif

      call broadcast_scalar (f_smice,    master_task)
      call broadcast_scalar (f_smliq,    master_task)
      call broadcast_scalar (f_rhos,     master_task)
      call broadcast_scalar (f_rsnw,     master_task)
      call broadcast_scalar (f_smicen,   master_task)
      call broadcast_scalar (f_smliqn,   master_task)
      call broadcast_scalar (f_rhosn,    master_task)
      call broadcast_scalar (f_rsnwn,    master_task)
      call broadcast_scalar (f_meltsliq, master_task)
      call broadcast_scalar (f_fsloss,   master_task)

      if (tr_snow) then

      ! 2D variables
      do ns = 1, nstreams
      if (histfreq(ns) /= 'x') then

      if (f_smice(1:1) /= 'x') &
         call define_hist_field(n_smice,"smice","kg/m^2",tstr2D, tcstr, &
             "ice mass per unit area in snow",                          &
             "none", c1, c0,                                            &
             ns, f_smice)

      if (f_smliq(1:1) /= 'x') &
         call define_hist_field(n_smliq,"smliq","kg/m^2",tstr2D, tcstr, &
             "liquid water mass per unit area in snow",                 &
             "none", c1, c0,                                            &
             ns, f_smliq)

      if (f_rhos(1:1) /= 'x') &
         call define_hist_field(n_rhos,"rhos","kg/m^3",tstr2D, tcstr, &
             "snow density: compaction",                              &
             "none", c1, c0,                                          &
             ns, f_rhos)

      if (f_rsnw(1:1) /= 'x') &
         call define_hist_field(n_rsnw,"rsnw","um",tstr2D, tcstr, &
             "average snow grain radius",                         &
             "none", c1, c0,                                      &
             ns, f_rsnw)

      if (f_meltsliq(1:1) /= 'x') &
         call define_hist_field(n_meltsliq,"meltsliq","kg/m^2/day",tstr2D, tcstr, &
             "snow liquid contribution to meltponds",                             &
             "none",secday/dt, c0,                                                &
             ns, f_meltsliq)

      if (f_fsloss(1:1) /= 'x') &
         call define_hist_field(n_fsloss,"fsloss","cm/day",tstr2D, tcstr, &
             "snow blown into leads (liquid)",                            &
             "none", mps_to_cmpdy/rhofresh, c0,                           &
             ns, f_fsloss)

      endif ! histfreq(ns) /= 'x'
      enddo ! nstreams
      endif ! tr_snow
      
      end subroutine init_hist_snow_2D

!=======================================================================

      subroutine init_hist_snow_3Dc

      use ice_calendar, only: nstreams, histfreq
      use ice_history_shared, only: tstr3Dc, tcstr, define_hist_field

      integer (kind=int_kind) :: ns
      logical (kind=log_kind) :: tr_snow
      character(len=*), parameter :: subname = '(init_hist_pond_3Dc)'
      
      call icepack_query_tracer_flags(tr_snow_out=tr_snow)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (tr_snow) then

      ! 3D (category) variables must be looped separately
      do ns = 1, nstreams
      if (histfreq(ns) /= 'x') then

      if (f_smicen(1:1) /= 'x') &
         call define_hist_field(n_smicen,"smicen","kg/m^2",tstr3Dc, tcstr, &
             "ice mass per unit area in snow, category",                   &
             "none", c1, c0,                                               &
             ns, f_smicen)

      if (f_smliqn(1:1) /= 'x') &
         call define_hist_field(n_smliqn,"smliqn","kg/m^2",tstr3Dc, tcstr, &
             "liquid water mass per unit area in snow, category",          &
             "none", c1, c0,                                               &
             ns, f_smliqn)

      if (f_rhosn(1:1) /= 'x') &
         call define_hist_field(n_rhosn,"rhosn","kg/m^3",tstr3Dc, tcstr, &
             "snow density: compaction, category",                       &
             "none", c1, c0,                                             &
             ns, f_rhosn)

      if (f_rsnwn(1:1) /= 'x') &
         call define_hist_field(n_rsnwn,"rsnwn","um",tstr3Dc, tcstr, &
             "average snow grain radius, category",                  &
             "none", c1, c0,                                         &
             ns, f_rsnwn)

      endif ! histfreq(ns) /= 'x'
      enddo ! ns

      endif ! tr_snow

      end subroutine init_hist_snow_3Dc

!=======================================================================

! accumulate average ice quantities or snapshots

      subroutine accum_hist_snow (iblk)

      use ice_arrays_column, only: meltsliq
      use ice_blocks, only: block, nx_block, ny_block
      use ice_domain, only: blocks_ice
      use ice_flux, only: fsloss
      use ice_history_shared, only: n2D, a2D, a3Dc, ncat_hist, &
          accum_hist_field, nzslyr
      use ice_state, only: vsno, vsnon, trcr, trcrn

      integer (kind=int_kind), intent(in) :: &
           iblk                 ! block index

      ! local variables

      integer (kind=int_kind) :: &
         i, j, k, n

      integer (kind=int_kind) :: &
         nt_smice, nt_smliq, nt_rhos, nt_rsnw

      logical (kind=log_kind) :: tr_snow

      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         worka

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat_hist) :: &
         workb

      character(len=*), parameter :: subname = '(accum_hist_snow)'

      !---------------------------------------------------------------
      ! increment field
      !---------------------------------------------------------------

      call icepack_query_tracer_flags(tr_snow_out=tr_snow)
      call icepack_query_tracer_indices(nt_smice_out=nt_smice, &
           nt_smliq_out=nt_smliq, nt_rhos_out=nt_rhos, nt_rsnw_out=nt_rsnw)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (allocated(a2D)) then
      if (tr_snow) then

         if (f_smice(1:1)/= 'x') then
            worka(:,:) = c0
            do k = 1, nzslyr
               worka(:,:) = worka(:,:) &
                          + trcr(:,:,nt_smice+k-1,iblk) &
                          * vsno(:,:,iblk) / real(nslyr,kind=dbl_kind)
            enddo
            call accum_hist_field(n_smice, iblk, worka, a2D) 
         endif
         if (f_smliq(1:1)/= 'x') then
            worka(:,:) = c0
            do k = 1, nzslyr
               worka(:,:) = worka(:,:) &
                          + trcr(:,:,nt_smliq+k-1,iblk) &
                          * vsno(:,:,iblk) / real(nslyr,kind=dbl_kind)
            enddo
            call accum_hist_field(n_smliq, iblk, worka, a2D) 
         endif
         if (f_rhos(1:1)/= 'x') then
            worka(:,:) = c0
            do k = 1, nzslyr
               worka(:,:) = worka(:,:) &
                          + trcr(:,:,nt_rhos+k-1,iblk) &
                          / real(nslyr,kind=dbl_kind)
            enddo
            call accum_hist_field(n_rhos, iblk, worka, a2D) 
         endif
         if (f_rsnw(1:1)/= 'x') then
            worka(:,:) = c0
            do k = 1, nzslyr
               worka(:,:) = worka(:,:) &
                          + trcr(:,:,nt_rsnw+k-1,iblk) &
                          / real(nslyr,kind=dbl_kind)
            enddo
            call accum_hist_field(n_rsnw, iblk, worka, a2D) 
         endif
         if (f_meltsliq(1:1)/= 'x') &
            call accum_hist_field(n_meltsliq, iblk, &
                 meltsliq(:,:,iblk), a2D) 
         if (f_fsloss(1:1)/= 'x') &
            call accum_hist_field(n_fsloss, iblk, &
                 fsloss(:,:,iblk), a2D) 
 
         endif ! allocated(a2D)

         ! 3D category fields
         if (allocated(a3Dc)) then
         if (f_smicen(1:1)/= 'x') then
            workb(:,:,:) = c0
            do k = 1, nzslyr
            do n = 1, ncat_hist
               workb(:,:,n) = workb(:,:,n) &
                            + trcrn(:,:,nt_smice+k-1,n,iblk) &
                            * vsnon(:,:,n,iblk) / real(nslyr,kind=dbl_kind)
            enddo
            enddo
            call accum_hist_field(n_smicen-n2D, iblk, ncat_hist, workb, a3Dc) 
         endif
         if (f_smliqn(1:1)/= 'x') then
            workb(:,:,:) = c0
            do k = 1, nzslyr
            do n = 1, ncat_hist
               workb(:,:,n) = workb(:,:,n) &
                            + trcrn(:,:,nt_smliq+k-1,n,iblk) &
                            * vsnon(:,:,n,iblk) / real(nslyr,kind=dbl_kind)
            enddo
            enddo
            call accum_hist_field(n_smliqn-n2D, iblk, ncat_hist, workb, a3Dc) 
         endif
         if (f_rhosn(1:1)/= 'x') then
            workb(:,:,:) = c0
            do k = 1, nzslyr
            do n = 1, ncat_hist
               workb(:,:,n) = workb(:,:,n) &
                            + trcrn(:,:,nt_rhos+k-1,n,iblk) &
                            / real(nslyr,kind=dbl_kind)
            enddo
            enddo
            call accum_hist_field(n_rhosn-n2D, iblk, ncat_hist, workb, a3Dc) 
         endif
         if (f_rsnwn(1:1)/= 'x') then
            workb(:,:,:) = c0
            do k = 1, nzslyr
            do n = 1, ncat_hist
               workb(:,:,n) = workb(:,:,n) &
                            + trcrn(:,:,nt_rsnw+k-1,n,iblk) &
                            / real(nslyr,kind=dbl_kind)
            enddo
            enddo
            call accum_hist_field(n_rsnwn-n2D, iblk, ncat_hist, workb, a3Dc) 
         endif
         endif ! allocated(a3Dc)

      endif ! tr_snow
      
      end subroutine accum_hist_snow

!=======================================================================

      end module ice_history_snow

!=======================================================================
