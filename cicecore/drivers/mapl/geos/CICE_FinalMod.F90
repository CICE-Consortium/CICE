!=======================================================================
!
!  This module contains routines for the final exit of the CICE model,
!  including final output and clean exit from any message passing
!  environments and frameworks.
!
!  authors: Philip W. Jones, LANL
!  2006: Converted to free source form (F90) by Elizabeth Hunke
!  2008: E. Hunke moved ESMF code to its own driver

      module CICE_FinalMod

      use ice_kinds_mod
      use ice_communicate, only: my_task, master_task
      use ice_exit, only: end_run, abort_ice
      use ice_fileunits, only: nu_diag, release_all_fileunits
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_parameters
      use icepack_intfc, only: icepack_query_tracer_flags, icepack_query_tracer_sizes

      implicit none
      private
      public :: CICE_Finalize, ice_checkpoint

!=======================================================================

      contains

!=======================================================================
!
!  This routine shuts down CICE by exiting all relevent environments.

      subroutine CICE_Finalize

      use ice_restart_shared, only: runid
      use ice_timers, only: ice_timer_stop, ice_timer_print_all, timer_total

      character(len=*), parameter :: subname = '(CICE_Finalize)'


      call ice_checkpoint

   !-------------------------------------------------------------------
   ! stop timers and print timer info
   !-------------------------------------------------------------------

      call ice_timer_stop(timer_total)        ! stop timing entire run
      call ice_timer_print_all(stats=.false.) ! print timing information

!echmod      if (nu_diag /= 6) close (nu_diag) ! diagnostic output
      call release_all_fileunits

   !-------------------------------------------------------------------
   ! quit MPI
   !-------------------------------------------------------------------

! standalone
!      call end_run       ! quit MPI

      end subroutine CICE_Finalize

!=======================================================================

!=======================================================================
      subroutine ice_checkpoint(time_stamp)

      use ice_boundary, only: ice_HaloUpdate
      use ice_calendar, only: dt, dt_dyn, ndtd, diagfreq, write_restart, istep
      use ice_calendar, only: idate, msec
      use ice_domain, only: halo_info, nblocks
      use ice_dyn_shared, only: kdyn, kridge
      use ice_dyn_eap, only: write_restart_eap
      use ice_restart, only: final_restart
      use ice_restart_shared, only: &
          restart_ext, restart_dir, restart_file, pointer_file, &
          runid, use_restart_time, lenstr, restart_coszen
      use ice_restart_column, only: write_restart_age, write_restart_FY, &
          write_restart_lvl, write_restart_pond_lvl, &
          write_restart_pond_topo, write_restart_aero, write_restart_fsd, &
          write_restart_iso, write_restart_bgc, write_restart_hbrine, &
          write_restart_snow
      use ice_restart_driver, only: dumpfile
      use ice_timers, only: ice_timer_start, ice_timer_stop, &
          timer_readwrite
      use ice_communicate, only: MPI_COMM_ICE


      character(len=*), intent(in), optional :: &
         time_stamp

      integer (kind=int_kind) :: &
         iblk        , & ! block index
         k           , & ! dynamics supercycling index
         ktherm          ! thermodynamics is off when ktherm = -1

      real (kind=dbl_kind) :: &
         offset          ! d(age)/dt time offset

      character(len=char_len_long) :: filename

      logical (kind=log_kind) :: &
          tr_iage, tr_FY, tr_lvl, tr_fsd, tr_snow, &
          tr_pond_lvl, tr_pond_topo, tr_brine, tr_iso, tr_aero, &
          calc_Tsfc, skl_bgc, z_tracers, wave_spec

      character(len=*), parameter :: subname = '(ice_checkpoint)'

      call icepack_query_parameters(calc_Tsfc_out=calc_Tsfc, skl_bgc_out=skl_bgc, &
           z_tracers_out=z_tracers, ktherm_out=ktherm, &
           wave_spec_out=wave_spec)
      call icepack_query_tracer_flags(tr_iage_out=tr_iage, tr_FY_out=tr_FY, &
           tr_lvl_out=tr_lvl, tr_pond_lvl_out=tr_pond_lvl, &
           tr_pond_topo_out=tr_pond_topo, tr_brine_out=tr_brine, tr_aero_out=tr_aero, &
           tr_iso_out=tr_iso, tr_fsd_out=tr_fsd, tr_snow_out=tr_snow)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      call ice_timer_start(timer_readwrite)  ! reading/writing

      if(present(time_stamp)) then
         filename = trim(restart_dir) // trim(restart_file) // '.' // trim(time_stamp)
      else
         filename = trim(restart_dir) // trim(restart_file)
      endif

      call dumpfile(filename_spec=trim(filename))   ! core variables for restarting
      if (tr_iage)      call write_restart_age
      if (tr_FY)        call write_restart_FY
      if (tr_lvl)       call write_restart_lvl
      !if (tr_pond_cesm) call write_restart_pond_cesm
      if (tr_pond_lvl)  call write_restart_pond_lvl
      if (tr_pond_topo) call write_restart_pond_topo
      if (tr_snow)      call write_restart_snow
      if (tr_fsd)       call write_restart_fsd
      if (tr_iso)       call write_restart_iso
      if (tr_aero)      call write_restart_aero
      if (skl_bgc .or. z_tracers) &
                        call write_restart_bgc
      if (tr_brine)     call write_restart_hbrine
      if (kdyn == 2)    call write_restart_eap
      call final_restart

      call ice_timer_stop(timer_readwrite)  ! reading/writing

      end subroutine ice_checkpoint

      end module CICE_FinalMod

!=======================================================================
