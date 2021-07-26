!=======================================================================
! Copyright (c) 2021, Triad National Security, LLC
! All rights reserved.
!                
! Copyright 2021. Triad National Security, LLC. This software was
! produced under U.S. Government contract DE-AC52-06NA25396 for Los 
! Alamos National Laboratory (LANL), which is operated by Triad
! National Security, LLC for the U.S. Department of Energy. The U.S.  
! Government has rights to use, reproduce, and distribute this software.  
! NEITHER THE GOVERNMENT NOR TRIAD NATIONAL SECURITY, LLC MAKES ANY  
! WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF
! THIS SOFTWARE. If software is modified to produce derivative works, 
! such modified software should be clearly marked, so as not to confuse 
! it with the version available from LANL.
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
      character(len=*), parameter :: subname='(icemodel)'

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
