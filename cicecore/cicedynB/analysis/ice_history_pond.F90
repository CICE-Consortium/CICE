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
           f_apondn    = 'm', f_apeffn     = 'm', &
           f_hpondn    = 'm',                     &
           f_apond     = 'x', f_apond_ai   = 'x', &
           f_hpond     = 'x', f_hpond_ai   = 'x', &
           f_ipond     = 'x', f_ipond_ai   = 'x', &
           f_apeff     = 'x', f_apeff_ai   = 'x'

      !---------------------------------------------------------------
      ! namelist variables
      !---------------------------------------------------------------

      namelist / icefields_pond_nml /     &
           f_apondn,    f_apeffn   , &
           f_hpondn,                 &
           f_apond,     f_apond_ai , &  
           f_hpond,     f_hpond_ai , &  
           f_ipond,     f_ipond_ai , &  
           f_apeff,     f_apeff_ai

      !---------------------------------------------------------------
      ! field indices
      !---------------------------------------------------------------

      integer (kind=int_kind), dimension(max_nstrm) :: &
           n_apondn      , n_apeffn    , & 
           n_hpondn      , &
           n_apond       , n_apond_ai, &
           n_hpond       , n_hpond_ai, &
           n_ipond       , n_ipond_ai, &
           n_apeff       , n_apeff_ai

!=======================================================================

      contains

!=======================================================================

      subroutine init_hist_pond_2D

      use ice_broadcast, only: broadcast_scalar
      use ice_calendar, only: nstreams
      use ice_communicate, only: my_task, master_task
      use ice_history_shared, only: tstr2D, tcstr, define_hist_field

      integer (kind=int_kind) :: ns
      integer (kind=int_kind) :: nml_error ! namelist i/o error flag
      logical (kind=log_kind) :: tr_pond
      character(len=*), parameter :: subname = '(init_hist_pond_2D)'

      !-----------------------------------------------------------------
      ! read namelist
      !-----------------------------------------------------------------

      call icepack_query_tracer_flags(tr_pond_out=tr_pond)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      call get_fileunit(nu_nml)
      if (my_task == master_task) then
         open (nu_nml, file=nml_filename, status='old',iostat=nml_error)
         if (nml_error /= 0) then
            nml_error = -1
         else
            nml_error =  1
         endif
         do while (nml_error > 0)
            read(nu_nml, nml=icefields_pond_nml,iostat=nml_error)
         end do
         if (nml_error == 0) close(nu_nml)
      endif
      call release_fileunit(nu_nml)

      call broadcast_scalar(nml_error, master_task)
      if (nml_error /= 0) then
         close (nu_nml)
         call abort_ice(subname//'ERROR: reading icefields_pond_nml')
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

      if (tr_pond) then

      ! 2D variables
      do ns = 1, nstreams

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

      enddo ! nstreams

      endif ! tr_pond
      
      end subroutine init_hist_pond_2D

!=======================================================================

      subroutine init_hist_pond_3Dc

      use ice_calendar, only: nstreams
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
         tr_pond_cesm, tr_pond_lvl, tr_pond_topo

      real (kind=dbl_kind) :: &
         puny

      type (block) :: &
         this_block           ! block information for current block

      character(len=*), parameter :: subname = '(accum_hist_pond)'

      !---------------------------------------------------------------
      ! increment field
      !---------------------------------------------------------------

         call icepack_query_parameters(puny_out=puny)
         call icepack_query_tracer_flags(tr_pond_cesm_out=tr_pond_cesm, &
              tr_pond_lvl_out=tr_pond_lvl, tr_pond_topo_out=tr_pond_topo)
         call icepack_query_tracer_indices(nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, &
              nt_alvl_out=nt_alvl, nt_ipnd_out=nt_ipnd)
         call icepack_warnings_flush(nu_diag)
         if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
            file=__FILE__, line=__LINE__)

         if (tr_pond_cesm) then
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

         elseif (tr_pond_lvl) then
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

         elseif (tr_pond_topo) then

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

         ! 3D category fields
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

      end subroutine accum_hist_pond

!=======================================================================

      end module ice_history_pond

!=======================================================================
