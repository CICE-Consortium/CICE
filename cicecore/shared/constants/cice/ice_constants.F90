!  SVN:$Id: ice_constants.F90 1228 2017-05-23 21:33:34Z tcraig $
!=======================================================================
!
! This module defines a variety of physical and numerical constants
! used throughout the ice model 
!
! author Elizabeth C. Hunke, LANL

      module ice_constants

      use ice_kinds_mod

      implicit none

      !-----------------------------------------------------------------
      ! physical constants
      !-----------------------------------------------------------------

      real (kind=dbl_kind), parameter, public :: &
         omega     = 7.292e-5_dbl_kind   ,&! angular velocity of earth (rad/sec)
         radius    = 6.37e6_dbl_kind       ! earth radius (m)

      real (kind=dbl_kind), parameter, public :: &
         spval_dbl = 1.0e30_dbl_kind    ! special value (double precision)

      real (kind=real_kind), parameter, public :: &
         spval     = 1.0e30_real_kind   ! special value for netCDF output

      ! these are currently set so as to have no effect on the decomposition
      real (kind=dbl_kind), parameter, public :: &
         shlat  =  30.0_dbl_kind   ,&! artificial masking edge (deg)
         nhlat  = -30.0_dbl_kind     ! artificial masking edge (deg)
   
      !-----------------------------------------------------------------
      ! numbers used outside the column package
      !-----------------------------------------------------------------

      real (kind=dbl_kind), parameter, public :: &
        c0   = 0.0_dbl_kind, &
        c1   = 1.0_dbl_kind, &
        c1p5 = 1.5_dbl_kind, &
        c2   = 2.0_dbl_kind, &
        c3   = 3.0_dbl_kind, &
        c4   = 4.0_dbl_kind, &
        c5   = 5.0_dbl_kind, &
        c6   = 6.0_dbl_kind, &
        c8   = 8.0_dbl_kind, &
        c9   = 9.0_dbl_kind, &
        c10  = 10.0_dbl_kind, &
        c12  = 12.0_dbl_kind, &
        c15  = 15.0_dbl_kind, &
        c16  = 16.0_dbl_kind, &
        c20  = 20.0_dbl_kind, &
        c25  = 25.0_dbl_kind, &
        c30  = 30.0_dbl_kind, &
        c100 = 100.0_dbl_kind, &
        c180 = 180.0_dbl_kind, &
        c360 = 360.0_dbl_kind, &
        c365 = 365.0_dbl_kind, &
	c400 = 400.0_dbl_kind, &
        c1000= 1000.0_dbl_kind, &
        c3600= 3600.0_dbl_kind, &
        p001 = 0.001_dbl_kind, &
        p01  = 0.01_dbl_kind, &
        p025 = 0.025_dbl_kind, &
        p05  = 0.05_dbl_kind, &
        p1   = 0.1_dbl_kind, &
        p15  = 0.15_dbl_kind, &
        p2   = 0.2_dbl_kind, &
        p25  = 0.25_dbl_kind, &
        p4   = 0.4_dbl_kind, &
        p5   = 0.5_dbl_kind, &
        p6   = 0.6_dbl_kind, &
        p75  = 0.75_dbl_kind, &
        p111 = c1/c9, &
        p166 = c1/c6, &
        p222 = c2/c9, &
        p333 = c1/c3, &
        p666 = c2/c3, &
        p055 = p111*p5, &
        p027 = p055*p5, &
        eps13  = 1.0e-13_dbl_kind, &
        eps16  = 1.0e-16_dbl_kind

      !-----------------------------------------------------------------
      ! location of fields for staggered grids
      !-----------------------------------------------------------------

      integer (int_kind), parameter, public :: &   
        field_loc_unknown  =  0, & 
        field_loc_noupdate = -1, & 
        field_loc_center   =  1, & 
        field_loc_NEcorner =  2, & 
        field_loc_Nface    =  3, & 
        field_loc_Eface    =  4, &
        field_loc_Wface    =  5

      !-----------------------------------------------------------------
      ! field type attribute - necessary for handling
      ! changes of direction across tripole boundary
      !-----------------------------------------------------------------

      integer (int_kind), parameter, public :: &   
        field_type_unknown  =  0, & 
        field_type_noupdate = -1, & 
        field_type_scalar   =  1, & 
        field_type_vector   =  2, & 
        field_type_angle    =  3

      !-----------------------------------------------------------------
      ! conversion factors
      !-----------------------------------------------------------------

      real (kind=dbl_kind), parameter, public :: &
        cm_to_m       = 0.01_dbl_kind   ,&! cm to meters
        m_to_cm       = 100._dbl_kind   ,&! meters to cm
        m2_to_km2     = 1.e-6_dbl_kind  ,&! m^2 to km^2
        kg_to_g       = 1000._dbl_kind  ,&! kilograms to grams
        mps_to_cmpdy  = 8.64e6_dbl_kind   ! m per s to cm per day

!=======================================================================

      end module ice_constants

!=======================================================================
