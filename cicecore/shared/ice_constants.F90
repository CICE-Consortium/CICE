!=======================================================================
!
! This module defines a variety of physical and numerical constants
! used throughout the ice model 
!
! author Elizabeth C. Hunke, LANL

      module ice_constants

      use ice_kinds_mod

      implicit none
      private

      public :: ice_init_constants
      public :: ice_query_constants

      !-----------------------------------------------------------------
      ! physical constants
      !-----------------------------------------------------------------

      real (kind=dbl_kind), public :: &
         omega     = 7.292e-5_dbl_kind   ,&! angular velocity of earth (rad/sec)
         radius    = 6.37e6_dbl_kind       ! earth radius (m)

      real (kind=dbl_kind), public :: &
         spval_dbl = 1.0e30_dbl_kind    ! special value (double precision)

      real (kind=real_kind), public :: &
         spval     = 1.0e30_real_kind   ! special value for netCDF output

      ! these are currently set so as to have no effect on the decomposition
      real (kind=dbl_kind), public :: &
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
        p3   = 0.3_dbl_kind, &
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
        eps04  = 1.0e-4_dbl_kind, &
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

      contains

!=======================================================================

! subroutine to set the cice constants

      subroutine ice_init_constants(   &
         omega_in, radius_in, spval_dbl_in, spval_in, shlat_in, nhlat_in)

      real (kind=dbl_kind), intent(in), optional :: &
         omega_in     , &   ! angular velocity of earth (rad/sec)
         radius_in    , &   ! earth radius (m)
         spval_dbl_in , &   ! special value (double precision)
         spval_in     , &   ! special value for netCDF output
         shlat_in     , &   ! artificial masking edge (deg)
         nhlat_in           ! artificial masking edge (deg)

      character(len=*),parameter :: subname='(ice_init_constants)'

      if (present(omega_in)) omega = omega_in
      if (present(radius_in)) radius = radius_in
      if (present(spval_dbl_in)) spval_dbl = spval_dbl_in
      if (present(spval_in)) spval = spval_in
      if (present(shlat_in)) shlat = shlat_in
      if (present(nhlat_in)) nhlat = nhlat_in

      end subroutine ice_init_constants

!=======================================================================

! subroutine to set the cice constants

      subroutine ice_query_constants(   &
         omega_out, radius_out, spval_dbl_out, spval_out, shlat_out, nhlat_out)

      real (kind=dbl_kind), intent(out), optional :: &
         omega_out     , &   ! angular velocity of earth (rad/sec)
         radius_out    , &   ! earth radius (m)
         spval_dbl_out , &   ! special value (double precision)
         spval_out     , &   ! special value for netCDF output
         shlat_out     , &   ! artificial masking edge (deg)
         nhlat_out           ! artificial masking edge (deg)

      character(len=*),parameter :: subname='(ice_query_constants)'

      if (present(omega_out)) omega_out = omega
      if (present(radius_out)) radius_out = radius
      if (present(spval_dbl_out)) spval_dbl_out = spval_dbl
      if (present(spval_out)) spval_out = spval
      if (present(shlat_out)) shlat_out = shlat
      if (present(nhlat_out)) nhlat_out = nhlat

      end subroutine ice_query_constants

!=======================================================================

      end module ice_constants

!=======================================================================
