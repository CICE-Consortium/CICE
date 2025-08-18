!=======================================================================
! Copyright 1998-2025, Triad National Security, LLC
! All rights reserved.
!
! This program was produced under U.S. Government contract 89233218CNA000001
! for Los Alamos National Laboratory (LANL), which is operated by Triad
! National Security, LLC for the U.S. Department of Energy/National Nuclear
! Security Administration. All rights in the program are reserved by Triad
! National Security, LLC, and the U.S. Department of Energy/National Nuclear
! Security Administration. The Government is granted for itself and others
! acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
! license in this material to reproduce, prepare. derivative works,
! distribute copies to the public, perform publicly and display publicly,
! and to permit others to do so.
!
! The full license and distribution policy are available from
! https://github.com/CICE-Consortium
!
!=======================================================================
!
! Main driver routine for CICE.  Initializes and steps through the model.
! This program should be compiled if CICE is run as a separate executable,
!  but not if CICE subroutines are called from another program (e.g., CAM).
!
! authors Elizabeth C. Hunke and William H. Lipscomb, LANL
!
! 2006: Converted to free form source (F90) by Elizabeth Hunke
! 2008: E. Hunke moved ESMF code to its own driver
!
      program icemodel

      use CICE_InitMod
      use CICE_RunMod
      use CICE_FinalMod

      implicit none

      !-----------------------------------------------------------------
      ! Initialize CICE
      !-----------------------------------------------------------------

      call CICE_Initialize

      !-----------------------------------------------------------------
      ! Run CICE
      !-----------------------------------------------------------------

      call CICE_Run

      !-----------------------------------------------------------------
      ! Finalize CICE
      !-----------------------------------------------------------------

      call CICE_Finalize

      end program icemodel

!=======================================================================
