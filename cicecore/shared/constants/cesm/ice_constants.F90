!  SVN:$Id: ice_constants.F90 890 2014-11-20 23:35:37Z eclare $
!=======================================================================
!
! This module defines a variety of physical and numerical constants
! used throughout the ice model 
!
! author Elizabeth C. Hunke, LANL

      module ice_constants

      use ice_kinds_mod
      use icepack_intfc ! all constants needed for column package

      implicit none
      save

      !-----------------------------------------------------------------
      ! physical constants
      !-----------------------------------------------------------------

      real (kind=dbl_kind), parameter, public :: &
         omega     = SHR_CONST_OMEGA ,&! angular velocity of earth (rad/sec)
         radius    = SHR_CONST_REARTH  ! earth radius (m)

      real (kind=dbl_kind), parameter, public :: &
         spval_dbl = SHR_CONST_SPVAL    ! special value

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
        c9   = 9.0_dbl_kind, &
        c12  = 12.0_dbl_kind, &
        c30  = 30.0_dbl_kind, &
        c180 = 180.0_dbl_kind, &
        c360 = 360.0_dbl_kind, &
        c365 = 365.0_dbl_kind, &
	c400 = 400.0_dbl_kind, &
        c3600= 3600.0_dbl_kind, &
        p025 = 0.025_dbl_kind, &
        p166 = c1/c6, &
        p111 = c1/c9, &
        p055 = p111*p5, &
        p027 = p055*p5, &
        p222 = c2/c9, &
        eps04  = 1.0e-4_dbl_kind, &
        eps11  = 1.0e-11_dbl_kind, &
        eps12  = 1.0e-12_dbl_kind, &
        eps13  = 1.0e-13_dbl_kind, &
        eps15  = 1.0e-15_dbl_kind, &
        eps16  = 1.0e-16_dbl_kind, &
        piq    = p5*pih, &
        pi2    = c2*pi

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
        mps_to_cmpdy  = 8.64e6_dbl_kind ,&! m per s to cm per day
        rad_to_deg    = 180._dbl_kind/pi  ! degree-radian conversion

#ifndef USE_ESMF
      integer (kind=int_kind), parameter :: &
         ESMF_SUCCESS = 0   ! otherwise ESMF defines this parameter
#endif

!=======================================================================

      end module ice_constants

!=======================================================================
