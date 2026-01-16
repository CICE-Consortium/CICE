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
           f_smassice  = 'm', f_smassicen   = 'x', &
           f_smassliq  = 'm', f_smassliqn   = 'x', &
           f_rhos_cmp  = 'm', f_rhos_cmpn   = 'x', &
           f_rhos_cnt  = 'm', f_rhos_cntn   = 'x', &
           f_rsnw      = 'm', f_rsnwn       = 'x', &
           f_meltsliq  = 'm', f_fsloss      = 'x'

      character (len=max_nstrm), public :: &
           f_sisndmasswind = 'm'

      !---------------------------------------------------------------
      ! namelist variables
      !---------------------------------------------------------------

      namelist / icefields_snow_nml / &
           f_smassice,  f_smassicen,  &
           f_smassliq,  f_smassliqn,  &
           f_rhos_cmp,  f_rhos_cmpn,  &
           f_rhos_cnt,  f_rhos_cntn,  &
           f_rsnw,      f_rsnwn,      &
           f_meltsliq,  f_fsloss,     &
           f_sisndmasswind

      !---------------------------------------------------------------
      ! field indices
      !---------------------------------------------------------------

      integer (kind=int_kind), dimension(max_nstrm), public :: &
           n_smassice,  n_smassicen, &
           n_smassliq,  n_smassliqn, &
           n_rhos_cmp,  n_rhos_cmpn, &
           n_rhos_cnt,  n_rhos_cntn, &
           n_rsnw,      n_rsnwn,     &
           n_meltsliq,  n_fsloss

      integer (kind=int_kind), dimension(max_nstrm), public :: &
           n_sisndmasswind

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
      use ice_fileunits, only: goto_nml

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      integer (kind=int_kind) :: ns
      integer (kind=int_kind) :: nml_error ! namelist i/o error flag
      real    (kind=dbl_kind) :: rhofresh, secday
      logical (kind=log_kind) :: tr_snow
      character(len=char_len_long) :: tmpstr2 ! for namelist check
      character(len=char_len)      :: nml_name ! for namelist check

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

         if (my_task == master_task) then
            nml_name = 'icefields_snow_nml'
            write(nu_diag,*) subname,' Reading ', trim(nml_name)

            ! open namelist file
            call get_fileunit(nu_nml)
            open (nu_nml, file=trim(nml_filename), status='old',iostat=nml_error)
            if (nml_error /= 0) then
               call abort_ice(subname//'ERROR: '//trim(nml_name)//' open file '// &
                  trim(nml_filename), &
                  file=__FILE__, line=__LINE__)
            endif

            ! goto namelist in file
            call goto_nml(nu_nml,trim(nml_name),nml_error)
            if (nml_error /= 0) then
               call abort_ice(subname//'ERROR: searching for '// trim(nml_name), &
                    file=__FILE__, line=__LINE__)
            endif

            ! read namelist
            nml_error =  1
            do while (nml_error > 0)
               read(nu_nml, nml=icefields_snow_nml,iostat=nml_error)
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

      else ! .not. tr_snow
          f_smassice = 'x'
          f_smassliq = 'x'
          f_rhos_cmp = 'x'
          f_rhos_cnt = 'x'
          f_rsnw     = 'x'
          f_smassicen= 'x'
          f_smassliqn= 'x'
          f_rhos_cmpn= 'x'
          f_rhos_cntn= 'x'
          f_rsnwn    = 'x'
          f_meltsliq = 'x'
          f_fsloss   = 'x'
          f_sisndmasswind = 'x'
      endif

      call broadcast_scalar (f_smassice, master_task)
      call broadcast_scalar (f_smassliq, master_task)
      call broadcast_scalar (f_rhos_cmp, master_task)
      call broadcast_scalar (f_rhos_cnt, master_task)
      call broadcast_scalar (f_rsnw,     master_task)
      call broadcast_scalar (f_smassicen,master_task)
      call broadcast_scalar (f_smassliqn,master_task)
      call broadcast_scalar (f_rhos_cmpn,master_task)
      call broadcast_scalar (f_rhos_cntn,master_task)
      call broadcast_scalar (f_rsnwn,    master_task)
      call broadcast_scalar (f_meltsliq, master_task)
      call broadcast_scalar (f_fsloss,   master_task)
      call broadcast_scalar (f_sisndmasswind,   master_task)

      if (tr_snow) then

      ! 2D variables
      do ns = 1, nstreams
      if (histfreq(ns) /= 'x') then

      if (f_smassice(1:1) /= 'x') &
         call define_hist_field(n_smassice,"smassice","kg/m^2",tstr2D, tcstr, &
             "ice mass per unit area in snow",                                &
             "none", c1, c0,                                                  &
             ns, f_smassice)

      if (f_smassliq(1:1) /= 'x') &
         call define_hist_field(n_smassliq,"smassliq","kg/m^2",tstr2D, tcstr, &
             "liquid mass per unit area in snow",                             &
             "none", c1, c0,                                                  &
             ns, f_smassliq)

      if (f_rhos_cmp(1:1) /= 'x') &
         call define_hist_field(n_rhos_cmp,"rhos_cmp","kg/m^3",tstr2D, tcstr, &
             "snow density: compaction",                                      &
             "none", c1, c0,                                                  &
             ns, f_rhos_cmp)

      if (f_rhos_cnt(1:1) /= 'x') &
         call define_hist_field(n_rhos_cnt,"rhos_cnt","kg/m^3",tstr2D, tcstr, &
             "snow density: content",                                         &
             "none", c1, c0,                                                  &
             ns, f_rhos_cnt)

      if (f_rsnw(1:1) /= 'x') &
         call define_hist_field(n_rsnw,"rsnw","10^-6 m",tstr2D, tcstr, &
             "average snow grain radius",                              &
             "none", c1, c0,                                           &
             ns, f_rsnw)

      if (f_meltsliq(1:1) /= 'x') &
         call define_hist_field(n_meltsliq,"meltsliq","kg/m^2/s",tstr2D, tcstr, &
             "snow liquid contribution to meltponds",                           &
             "none", c1/dt, c0,                                                 &
             ns, f_meltsliq)

      if (f_fsloss(1:1) /= 'x') &
         call define_hist_field(n_fsloss,"fsloss","kg/m^2/s",tstr2D, tcstr, &
             "rate of snow loss to leads (liquid)",                         &
             "none", c1, c0,                                                &
             ns, f_fsloss)

      if (f_sisndmasswind(1:1) /= 'x') &
         call define_hist_field(n_sisndmasswind,"sisndmasswind","kg/m^2/s",tstr2D, tcstr, &
             "snow mass rate of change through wind drift of snow",                       &
             "rate of change of snow mass due to wind-driven transport into the ocean",   &
             c1, c0,                                                                      &
             ns, f_sisndmasswind, avg_ice_present='none', mask_ice_free_points=.false.)

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

      if (f_smassicen(1:1) /= 'x') &
         call define_hist_field(n_smassicen,"smassicen","kg/m^2",tstr3Dc, tcstr, &
             "ice mass per unit area in snow, category",                         &
             "none", c1, c0,                                                     &
             ns, f_smassicen)

      if (f_smassliqn(1:1) /= 'x') &
         call define_hist_field(n_smassliqn,"smassliqn","kg/m^2",tstr3Dc, tcstr, &
             "liquid mass per unit area in snow, category",                      &
             "none", c1, c0,                                                     &
             ns, f_smassliqn)

      if (f_rhos_cmpn(1:1) /= 'x') &
         call define_hist_field(n_rhos_cmpn,"rhos_cmpn","kg/m^3",tstr3Dc, tcstr, &
             "snow density: compaction, category",                               &
             "none", c1, c0,                                                     &
             ns, f_rhos_cmpn)

      if (f_rhos_cntn(1:1) /= 'x') &
         call define_hist_field(n_rhos_cntn,"rhos_cntn","kg/m^3",tstr3Dc, tcstr, &
             "snow density: content, category",                                  &
             "none", c1, c0,                                                     &
             ns, f_rhos_cntn)

      if (f_rsnwn(1:1) /= 'x') &
         call define_hist_field(n_rsnwn,"rsnwn","10^-6 m",tstr3Dc, tcstr, &
             "average snow grain radius, category",                       &
             "none", c1, c0,                                              &
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
      use ice_flux, only: fsloss
      use ice_history_shared, only: n2D, a2D, a3Dc, ncat_hist, &
          accum_hist_field, nzslyr
      use ice_state, only: vsno, vsnon, trcr, trcrn

      integer (kind=int_kind), intent(in) :: &
           iblk                 ! block index

      ! local variables

      integer (kind=int_kind) :: &
         k, n

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

         if (f_smassice(1:1)/= 'x') then
            worka(:,:) = c0
            do k = 1, nzslyr
               worka(:,:) = worka(:,:) &
                          + trcr(:,:,nt_smice+k-1,iblk)
            enddo
            worka(:,:) = worka(:,:) * vsno(:,:,iblk) / real(nslyr,kind=dbl_kind)
            call accum_hist_field(n_smassice, iblk, worka, a2D)
         endif
         if (f_smassliq(1:1)/= 'x') then
            worka(:,:) = c0
            do k = 1, nzslyr
               worka(:,:) = worka(:,:) &
                          + trcr(:,:,nt_smliq+k-1,iblk)
            enddo
            worka(:,:) = worka(:,:) * vsno(:,:,iblk) / real(nslyr,kind=dbl_kind)
            call accum_hist_field(n_smassliq, iblk, worka, a2D)
         endif
         if (f_rhos_cmp(1:1)/= 'x') then
            worka(:,:) = c0
            do k = 1, nzslyr
               worka(:,:) = worka(:,:) &
                          + trcr(:,:,nt_rhos+k-1,iblk)
            enddo
            worka(:,:) = worka(:,:) / real(nslyr,kind=dbl_kind)
            call accum_hist_field(n_rhos_cmp, iblk, worka, a2D)
         endif
         if (f_rhos_cnt(1:1)/= 'x') then
            worka(:,:) = c0
            do k = 1, nzslyr
               worka(:,:) = worka(:,:) &
                          + trcr(:,:,nt_smice+k-1,iblk) &
                          + trcr(:,:,nt_smliq+k-1,iblk)
            enddo
            worka(:,:) = worka(:,:) / real(nslyr,kind=dbl_kind)
            call accum_hist_field(n_rhos_cnt, iblk, worka, a2D)
         endif
         if (f_rsnw(1:1)/= 'x') then
            worka(:,:) = c0
            do k = 1, nzslyr
               worka(:,:) = worka(:,:) &
                          + trcr(:,:,nt_rsnw+k-1,iblk)
            enddo
            worka(:,:) = worka(:,:) / real(nslyr,kind=dbl_kind)
            call accum_hist_field(n_rsnw, iblk, worka, a2D)
         endif
         if (f_meltsliq(1:1)/= 'x') &
            call accum_hist_field(n_meltsliq, iblk, &
                 meltsliq(:,:,iblk), a2D)
         if (f_fsloss(1:1)/= 'x') &
            call accum_hist_field(n_fsloss, iblk, &
                 fsloss(:,:,iblk), a2D)

         if (f_sisndmasswind(1:1)/= 'x') &
            call accum_hist_field(n_sisndmasswind, iblk, &
                 fsloss(:,:,iblk), a2D)

         endif ! allocated(a2D)

         ! 3D category fields
         if (allocated(a3Dc)) then
         if (f_smassicen(1:1)/= 'x') then
            workb(:,:,:) = c0
            do n = 1, ncat_hist
               do k = 1, nzslyr
                  workb(:,:,n) = workb(:,:,n) &
                               + trcrn(:,:,nt_smice+k-1,n,iblk)
               enddo
               workb(:,:,n) = workb(:,:,n) &
                            * vsnon(:,:,n,iblk) / real(nslyr,kind=dbl_kind)
            enddo
            call accum_hist_field(n_smassicen-n2D, iblk, ncat_hist, workb, a3Dc)
         endif
         if (f_smassliqn(1:1)/= 'x') then
            workb(:,:,:) = c0
            do n = 1, ncat_hist
               do k = 1, nzslyr
                  workb(:,:,n) = workb(:,:,n) &
                               + trcrn(:,:,nt_smliq+k-1,n,iblk)
               enddo
               workb(:,:,n) = workb(:,:,n) &
                            * vsnon(:,:,n,iblk) / real(nslyr,kind=dbl_kind)
            enddo
            call accum_hist_field(n_smassliqn-n2D, iblk, ncat_hist, workb, a3Dc)
         endif
         if (f_rhos_cmpn(1:1)/= 'x') then
            workb(:,:,:) = c0
            do n = 1, ncat_hist
               do k = 1, nzslyr
                  workb(:,:,n) = workb(:,:,n) &
                               + trcrn(:,:,nt_rhos+k-1,n,iblk)
               enddo
               workb(:,:,n) = workb(:,:,n) / real(nslyr,kind=dbl_kind)
            enddo
            call accum_hist_field(n_rhos_cmpn-n2D, iblk, ncat_hist, workb, a3Dc)
         endif
         if (f_rhos_cntn(1:1)/= 'x') then
            workb(:,:,:) = c0
            do n = 1, ncat_hist
               do k = 1, nzslyr
                  workb(:,:,n) = workb(:,:,n) &
                               + trcrn(:,:,nt_smice+k-1,n,iblk) &
                               + trcrn(:,:,nt_smliq+k-1,n,iblk)
               enddo
               workb(:,:,n) = workb(:,:,n) / real(nslyr,kind=dbl_kind)
            enddo
            call accum_hist_field(n_rhos_cntn-n2D, iblk, ncat_hist, workb, a3Dc)
         endif
         if (f_rsnwn(1:1)/= 'x') then
            workb(:,:,:) = c0
            do n = 1, ncat_hist
               do k = 1, nzslyr
                  workb(:,:,n) = workb(:,:,n) &
                               + trcrn(:,:,nt_rsnw+k-1,n,iblk)
               enddo
               workb(:,:,n) = workb(:,:,n) / real(nslyr,kind=dbl_kind)
            enddo
            call accum_hist_field(n_rsnwn-n2D, iblk, ncat_hist, workb, a3Dc)
         endif
         endif ! allocated(a3Dc)

      endif ! tr_snow

      end subroutine accum_hist_snow

!=======================================================================

      end module ice_history_snow

!=======================================================================
