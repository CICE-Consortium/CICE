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

      implicit none
      private
      public :: CICE_Finalize

!=======================================================================

      contains

!=======================================================================
!
!  This routine shuts down CICE by exiting all relevent environments.

      subroutine CICE_Finalize

      use ice_restart_shared, only: runid
      use ice_timers, only: ice_timer_stop, ice_timer_print_all, timer_total

      character(len=*), parameter :: subname = '(CICE_Finalize)'

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

      end module CICE_FinalMod

!=======================================================================
