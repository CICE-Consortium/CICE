!=======================================================================

! Flux variable declarations; these include fields sent from the coupler
! ("in"), sent to the coupler ("out"), written to diagnostic history files
! ("diagnostic"), and used internally ("internal").
!
! author Elizabeth C. Hunke, LANL
!
! 2004: Block structure added by William Lipscomb
!       Swappped, revised, and added some subroutines
! 2006: Converted to free source form (F90) by Elizabeth Hunke

      module ice_flux

      use ice_kinds_mod
      use ice_fileunits, only: nu_diag
      use ice_blocks, only: nx_block, ny_block
      use ice_domain_size, only: max_blocks, ncat, max_nstrm, nilyr
      use ice_constants, only: c0, c1, c5, c10, c20, c180
      use ice_exit, only: abort_ice
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_tracer_flags, icepack_query_tracer_indices
      use icepack_intfc, only: icepack_query_parameters
      use icepack_intfc, only: icepack_liquidus_temperature

      implicit none
      private
      public :: init_coupler_flux, init_history_therm, init_history_dyn, &
                init_flux_ocn, init_flux_atm, scale_fluxes, alloc_flux

      character (char_len), public :: &
         default_season ! seasonal default values for forcing

      !-----------------------------------------------------------------
      ! Dynamics component
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &

       ! in from atmos (if .not.calc_strair)  
         strax   , & ! wind stress components (N/m^2)
         stray   , & ! 

       ! in from ocean
         uocn    , & ! ocean current, x-direction (m/s)
         vocn    , & ! ocean current, y-direction (m/s)
         ss_tltx , & ! sea surface slope, x-direction (m/m)
         ss_tlty , & ! sea surface slope, y-direction
         hwater  , & ! water depth for seabed stress calc (landfast ice) 

       ! out to atmosphere
         strairxT, & ! stress on ice by air, x-direction
         strairyT, & ! stress on ice by air, y-direction

       ! out to ocean          T-cell (kg/m s^2)
       ! Note, CICE_IN_NEMO uses strocnx and strocny for coupling
         strocnxT, & ! ice-ocean stress, x-direction
         strocnyT    ! ice-ocean stress, y-direction

       ! diagnostic

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         sig1    , & ! normalized principal stress component
         sig2    , & ! normalized principal stress component
         sigP    , & ! internal ice pressure (N/m)
         taubx   , & ! seabed stress (x) (N/m^2)
         tauby   , & ! seabed stress (y) (N/m^2)
         strairx , & ! stress on ice by air, x-direction
         strairy , & ! stress on ice by air, y-direction
         strocnx , & ! ice-ocean stress, x-direction
         strocny , & ! ice-ocean stress, y-direction
         strtltx , & ! stress due to sea surface slope, x-direction
         strtlty , & ! stress due to sea surface slope, y-direction
         strintx , & ! divergence of internal ice stress, x (N/m^2)
         strinty , & ! divergence of internal ice stress, y (N/m^2)
         daidtd  , & ! ice area tendency due to transport   (1/s)
         dvidtd  , & ! ice volume tendency due to transport (m/s)
         dagedtd , & ! ice age tendency due to transport (s/s)
         dardg1dt, & ! rate of area loss by ridging ice (1/s)
         dardg2dt, & ! rate of area gain by new ridges (1/s)
         dvirdgdt, & ! rate of ice volume ridged (m/s)
         opening     ! rate of opening due to divergence/shear (1/s)

      real (kind=dbl_kind), & 
         dimension (:,:,:,:), allocatable, public :: &
       ! ridging diagnostics in categories
         dardg1ndt, & ! rate of area loss by ridging ice (1/s)
         dardg2ndt, & ! rate of area gain by new ridges (1/s)
         dvirdgndt, & ! rate of ice volume ridged (m/s)
         aparticn,  & ! participation function
         krdgn,     & ! mean ridge thickness/thickness of ridging ice
         ardgn,     & ! fractional area of ridged ice
         vrdgn,     & ! volume of ridged ice
         araftn,    & ! rafting ice area
         vraftn,    & ! rafting ice volume 
         aredistn,  & ! redistribution function: fraction of new ridge area
         vredistn     ! redistribution function: fraction of new ridge volume

       ! restart

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
       ! ice stress tensor in each corner of T cell (kg/s^2)
         stressp_1, stressp_2, stressp_3, stressp_4 , & ! sigma11+sigma22
         stressm_1, stressm_2, stressm_3, stressm_4 , & ! sigma11-sigma22
         stress12_1,stress12_2,stress12_3,stress12_4    ! sigma12

      logical (kind=log_kind), &
         dimension (:,:,:), allocatable, public :: &
         iceumask   ! ice extent mask (U-cell)

       ! internal

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         fm       , & ! Coriolis param. * mass in U-cell (kg/s)
         Tbu          ! factor for seabed stress (N/m^2)

      !-----------------------------------------------------------------
      ! Thermodynamic component
      !-----------------------------------------------------------------

       ! in from atmosphere (if calc_Tsfc)

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         zlvl    , & ! atm level height (m)
         uatm    , & ! wind velocity components (m/s)
         vatm    , &
         wind    , & ! wind speed (m/s)
         potT    , & ! air potential temperature  (K)
         Tair    , & ! air temperature  (K)
         Qa      , & ! specific humidity (kg/kg)
         rhoa    , & ! air density (kg/m^3)
         swvdr   , & ! sw down, visible, direct  (W/m^2)
         swvdf   , & ! sw down, visible, diffuse (W/m^2)
         swidr   , & ! sw down, near IR, direct  (W/m^2)
         swidf   , & ! sw down, near IR, diffuse (W/m^2)
         flw         ! incoming longwave radiation (W/m^2)

       ! in from atmosphere (if .not. Tsfc_calc)
       ! required for coupling to HadGEM3
       ! NOTE: when in CICE_IN_NEMO mode, these are gridbox mean fields,
       ! not per ice area. When in standalone mode, these are per ice area.

      real (kind=dbl_kind), & 
         dimension (:,:,:,:), allocatable, public :: &
         fsurfn_f   , & ! net flux to top surface, excluding fcondtop
         fcondtopn_f, & ! downward cond flux at top surface (W m-2)
         fsensn_f   , & ! sensible heat flux (W m-2)
         flatn_f        ! latent heat flux (W m-2)

       ! in from atmosphere

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         frain   , & ! rainfall rate (kg/m^2 s)
         fsnow       ! snowfall rate (kg/m^2 s)

       ! in from ocean

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         sss     , & ! sea surface salinity (ppt)
         sst     , & ! sea surface temperature (C)
         frzmlt  , & ! freezing/melting potential (W/m^2)
         frzmlt_init, & ! frzmlt used in current time step (W/m^2)
         Tf      , & ! freezing temperature (C)
         qdp     , & ! deep ocean heat flux (W/m^2), negative upward
         hmix    , & ! mixed layer depth (m)
         daice_da    ! data assimilation concentration increment rate 
                     ! (concentration s-1)(only used in hadgem drivers)

       ! out to atmosphere (if calc_Tsfc)
       ! note Tsfc is in ice_state.F

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         fsens   , & ! sensible heat flux (W/m^2)
         flat    , & ! latent heat flux   (W/m^2)
         fswabs  , & ! shortwave flux absorbed in ice and ocean (W/m^2)
         fswint_ai, & ! SW absorbed in ice interior below surface (W/m^2)
         flwout  , & ! outgoing longwave radiation (W/m^2)
         Tref    , & ! 2m atm reference temperature (K)
         Qref    , & ! 2m atm reference spec humidity (kg/kg)
         Uref    , & ! 10m atm reference wind speed (m/s)
         evap    , & ! evaporative water flux (kg/m^2/s)
         evaps   , & ! evaporative water flux over snow (kg/m^2/s)
         evapi       ! evaporative water flux over ice (kg/m^2/s)

       ! albedos aggregated over categories (if calc_Tsfc)
      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         alvdr   , & ! visible, direct   (fraction)
         alidr   , & ! near-ir, direct   (fraction)
         alvdf   , & ! visible, diffuse  (fraction)
         alidf   , & ! near-ir, diffuse  (fraction)
         ! grid-box-mean versions
         alvdr_ai, & ! visible, direct   (fraction)
         alidr_ai, & ! near-ir, direct   (fraction)
         alvdf_ai, & ! visible, diffuse  (fraction)
         alidf_ai, & ! near-ir, diffuse  (fraction)
         ! components for history
         albice    , & ! bare ice albedo
         albsno    , & ! snow albedo
         albpnd    , & ! melt pond albedo
         apeff_ai  , & ! effective pond area used for radiation calculation
         snowfrac  , & ! snow fraction used in radiation
         ! components for diagnostic
         alvdr_init, & ! visible, direct   (fraction)
         alidr_init, & ! near-ir, direct   (fraction)
         alvdf_init, & ! visible, diffuse  (fraction)
         alidf_init    ! near-ir, diffuse  (fraction)

      real (kind=dbl_kind), &
         dimension(:,:,:,:), allocatable, public :: &
         albcnt       ! counter for zenith angle

       ! out to ocean 
       ! (Note CICE_IN_NEMO does not use these for coupling.  
       !  It uses fresh_ai,fsalt_ai,fhocn_ai and fswthru_ai)
      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         fpond   , & ! fresh water flux to ponds (kg/m^2/s)
         fresh   , & ! fresh water flux to ocean (kg/m^2/s)
         fsalt   , & ! salt flux to ocean (kg/m^2/s)
         fhocn   , & ! net heat flux to ocean (W/m^2)
         fswthru , & ! shortwave penetrating to ocean (W/m^2)
         fswthru_vdr , & ! vis dir shortwave penetrating to ocean (W/m^2)
         fswthru_vdf , & ! vis dif shortwave penetrating to ocean (W/m^2)
         fswthru_idr , & ! nir dir shortwave penetrating to ocean (W/m^2)
         fswthru_idf     ! nir dif shortwave penetrating to ocean (W/m^2)

       ! internal

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         scale_factor! scaling factor for shortwave components

      logical (kind=log_kind), public :: &
         update_ocn_f, & ! if true, update fresh water and salt fluxes
         l_mpond_fresh   ! if true, include freshwater feedback from meltponds
                         ! when running in ice-ocean or coupled configuration

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         meltsn      , & ! snow melt in category n (m)
         melttn      , & ! top melt in category n (m)
         meltbn      , & ! bottom melt in category n (m)
         congeln     , & ! congelation ice formation in category n (m)
         snoicen         ! snow-ice formation in category n (m)

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         keffn_top       ! effective thermal conductivity of the top ice layer 
                         ! on categories (W/m^2/K)

      ! quantities passed from ocean mixed layer to atmosphere
      ! (for running with CAM)

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         strairx_ocn , & ! stress on ocean by air, x-direction
         strairy_ocn , & ! stress on ocean by air, y-direction
         fsens_ocn   , & ! sensible heat flux (W/m^2)
         flat_ocn    , & ! latent heat flux   (W/m^2)
         flwout_ocn  , & ! outgoing longwave radiation (W/m^2)
         evap_ocn    , & ! evaporative water flux (kg/m^2/s)
         alvdr_ocn   , & ! visible, direct   (fraction)
         alidr_ocn   , & ! near-ir, direct   (fraction)
         alvdf_ocn   , & ! visible, diffuse  (fraction)
         alidf_ocn   , & ! near-ir, diffuse  (fraction)
         Tref_ocn    , & ! 2m atm reference temperature (K)
         Qref_ocn        ! 2m atm reference spec humidity (kg/kg)

      ! diagnostic

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         fsurf , & ! net surface heat flux (excluding fcondtop)(W/m^2)
         fcondtop,&! top surface conductive flux        (W/m^2)
         fcondbot,&! bottom surface conductive flux     (W/m^2)
         fbot,   & ! heat flux at bottom surface of ice (excluding excess) (W/m^2)
         Tbot,   & ! temperature at bottom surface of ice (deg C)
         Tsnice,  & ! temperature at snow ice interface (deg C)
         congel, & ! basal ice growth         (m/step-->cm/day)
         frazil, & ! frazil ice growth        (m/step-->cm/day)
         snoice, & ! snow-ice formation       (m/step-->cm/day)
         meltt , & ! top ice melt             (m/step-->cm/day)
         melts , & ! snow melt                (m/step-->cm/day)
         meltb , & ! basal ice melt           (m/step-->cm/day)
         meltl , & ! lateral ice melt         (m/step-->cm/day)
         dsnow,  & ! change in snow thickness (m/step-->cm/day)
         daidtt, & ! ice area tendency thermo.   (s^-1)
         dvidtt, & ! ice volume tendency thermo. (m/s)
         dagedtt,& ! ice age tendency thermo.    (s/s)
         mlt_onset, &! day of year that sfc melting begins
         frz_onset, &! day of year that freezing begins (congel or frazil)
         frazil_diag ! frazil ice growth diagnostic (m/step-->cm/day)
         
      real (kind=dbl_kind), & 
         dimension (:,:,:,:), allocatable, public :: &
         fsurfn,   & ! category fsurf
         fcondtopn,& ! category fcondtop
         fcondbotn,& ! category fcondbot
         fsensn,   & ! category sensible heat flux
         flatn       ! category latent heat flux

      ! As above but these remain grid box mean values i.e. they are not
      ! divided by aice at end of ice_dynamics.  These are used in
      ! CICE_IN_NEMO for coupling and also for generating
      ! ice diagnostics and history files as these are more accurate. 
      ! (The others suffer from problem of incorrect values at grid boxes
      !  that change from an ice free state to an icy state.)
    
      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         fresh_ai, & ! fresh water flux to ocean (kg/m^2/s)
         fsalt_ai, & ! salt flux to ocean (kg/m^2/s)
         fhocn_ai, & ! net heat flux to ocean (W/m^2)
         fswthru_ai  ! shortwave penetrating to ocean (W/m^2)

      ! Used with data assimilation in hadgem drivers
      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         fresh_da, & ! fresh water flux to ocean due to data assim (kg/m^2/s)
         fsalt_da    ! salt flux to ocean due to data assimilation(kg/m^2/s)

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         fswthrun_ai  ! per-category fswthru * ai (W/m^2)
 
      logical (kind=log_kind), public :: send_i2x_per_cat = .false.

      !-----------------------------------------------------------------
      ! internal
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         rside   , & ! fraction of ice that melts laterally
         fside   , & ! lateral heat flux (W/m^2)
         fsw     , & ! incoming shortwave radiation (W/m^2)
         coszen  , & ! cosine solar zenith angle, < 0 for sun below horizon 
         rdg_conv, & ! convergence term for ridging (1/s)
         rdg_shear   ! shear term for ridging (1/s)
 
      real (kind=dbl_kind), dimension(:,:,:,:), allocatable, public :: &
         salinz    ,&   ! initial salinity  profile (ppt)   
         Tmltz          ! initial melting temperature (^oC)

!=======================================================================

      contains

!=======================================================================
!
! Allocate space for all variables 
!
      subroutine alloc_flux

      integer (int_kind) :: ierr

      allocate( &
         strax      (nx_block,ny_block,max_blocks), & ! wind stress components (N/m^2)
         stray      (nx_block,ny_block,max_blocks), & ! 
         uocn       (nx_block,ny_block,max_blocks), & ! ocean current, x-direction (m/s)
         vocn       (nx_block,ny_block,max_blocks), & ! ocean current, y-direction (m/s)
         ss_tltx    (nx_block,ny_block,max_blocks), & ! sea surface slope, x-direction (m/m)
         ss_tlty    (nx_block,ny_block,max_blocks), & ! sea surface slope, y-direction
         hwater     (nx_block,ny_block,max_blocks), & ! water depth for seabed stress calc (landfast ice) 
         strairxT   (nx_block,ny_block,max_blocks), & ! stress on ice by air, x-direction
         strairyT   (nx_block,ny_block,max_blocks), & ! stress on ice by air, y-direction
         strocnxT   (nx_block,ny_block,max_blocks), & ! ice-ocean stress, x-direction
         strocnyT   (nx_block,ny_block,max_blocks), & ! ice-ocean stress, y-direction
         sig1       (nx_block,ny_block,max_blocks), & ! normalized principal stress component
         sig2       (nx_block,ny_block,max_blocks), & ! normalized principal stress component
         sigP       (nx_block,ny_block,max_blocks), & ! internal ice pressure (N/m)
         taubx      (nx_block,ny_block,max_blocks), & ! seabed stress (x) (N/m^2)
         tauby      (nx_block,ny_block,max_blocks), & ! seabed stress (y) (N/m^2)
         strairx    (nx_block,ny_block,max_blocks), & ! stress on ice by air, x-direction
         strairy    (nx_block,ny_block,max_blocks), & ! stress on ice by air, y-direction
         strocnx    (nx_block,ny_block,max_blocks), & ! ice-ocean stress, x-direction
         strocny    (nx_block,ny_block,max_blocks), & ! ice-ocean stress, y-direction
         strtltx    (nx_block,ny_block,max_blocks), & ! stress due to sea surface slope, x-direction
         strtlty    (nx_block,ny_block,max_blocks), & ! stress due to sea surface slope, y-direction
         strintx    (nx_block,ny_block,max_blocks), & ! divergence of internal ice stress, x (N/m^2)
         strinty    (nx_block,ny_block,max_blocks), & ! divergence of internal ice stress, y (N/m^2)
         daidtd     (nx_block,ny_block,max_blocks), & ! ice area tendency due to transport   (1/s)
         dvidtd     (nx_block,ny_block,max_blocks), & ! ice volume tendency due to transport (m/s)
         dagedtd    (nx_block,ny_block,max_blocks), & ! ice age tendency due to transport (s/s)
         dardg1dt   (nx_block,ny_block,max_blocks), & ! rate of area loss by ridging ice (1/s)
         dardg2dt   (nx_block,ny_block,max_blocks), & ! rate of area gain by new ridges (1/s)
         dvirdgdt   (nx_block,ny_block,max_blocks), & ! rate of ice volume ridged (m/s)
         opening    (nx_block,ny_block,max_blocks), & ! rate of opening due to divergence/shear (1/s)
         stressp_1  (nx_block,ny_block,max_blocks), & ! sigma11+sigma22
         stressp_2  (nx_block,ny_block,max_blocks), & ! sigma11+sigma22
         stressp_3  (nx_block,ny_block,max_blocks), & ! sigma11+sigma22
         stressp_4  (nx_block,ny_block,max_blocks), & ! sigma11+sigma22
         stressm_1  (nx_block,ny_block,max_blocks), & ! sigma11-sigma22
         stressm_2  (nx_block,ny_block,max_blocks), & ! sigma11-sigma22
         stressm_3  (nx_block,ny_block,max_blocks), & ! sigma11-sigma22
         stressm_4  (nx_block,ny_block,max_blocks), & ! sigma11-sigma22
         stress12_1 (nx_block,ny_block,max_blocks), & ! sigma12
         stress12_2 (nx_block,ny_block,max_blocks), & ! sigma12
         stress12_3 (nx_block,ny_block,max_blocks), & ! sigma12
         stress12_4 (nx_block,ny_block,max_blocks), & ! sigma12
         iceumask   (nx_block,ny_block,max_blocks), & ! ice extent mask (U-cell)
         fm         (nx_block,ny_block,max_blocks), & ! Coriolis param. * mass in U-cell (kg/s)
         Tbu        (nx_block,ny_block,max_blocks), & ! factor for seabed stress (landfast ice)
         zlvl       (nx_block,ny_block,max_blocks), & ! atm level height (m)
         uatm       (nx_block,ny_block,max_blocks), & ! wind velocity components (m/s)
         vatm       (nx_block,ny_block,max_blocks), &
         wind       (nx_block,ny_block,max_blocks), & ! wind speed (m/s)
         potT       (nx_block,ny_block,max_blocks), & ! air potential temperature  (K)
         Tair       (nx_block,ny_block,max_blocks), & ! air temperature  (K)
         Qa         (nx_block,ny_block,max_blocks), & ! specific humidity (kg/kg)
         rhoa       (nx_block,ny_block,max_blocks), & ! air density (kg/m^3)
         swvdr      (nx_block,ny_block,max_blocks), & ! sw down, visible, direct  (W/m^2)
         swvdf      (nx_block,ny_block,max_blocks), & ! sw down, visible, diffuse (W/m^2)
         swidr      (nx_block,ny_block,max_blocks), & ! sw down, near IR, direct  (W/m^2)
         swidf      (nx_block,ny_block,max_blocks), & ! sw down, near IR, diffuse (W/m^2)
         flw        (nx_block,ny_block,max_blocks), & ! incoming longwave radiation (W/m^2)
         frain      (nx_block,ny_block,max_blocks), & ! rainfall rate (kg/m^2 s)
         fsnow      (nx_block,ny_block,max_blocks), & ! snowfall rate (kg/m^2 s)
         sss        (nx_block,ny_block,max_blocks), & ! sea surface salinity (ppt)
         sst        (nx_block,ny_block,max_blocks), & ! sea surface temperature (C)
         frzmlt     (nx_block,ny_block,max_blocks), & ! freezing/melting potential (W/m^2)
         frzmlt_init(nx_block,ny_block,max_blocks), & ! frzmlt used in current time step (W/m^2)
         Tf         (nx_block,ny_block,max_blocks), & ! freezing temperature (C)
         qdp        (nx_block,ny_block,max_blocks), & ! deep ocean heat flux (W/m^2), negative upward
         hmix       (nx_block,ny_block,max_blocks), & ! mixed layer depth (m)
         daice_da   (nx_block,ny_block,max_blocks), & ! data assimilation concentration increment rate (concentration s-1)(only used in hadgem drivers)
         fsens      (nx_block,ny_block,max_blocks), & ! sensible heat flux (W/m^2)
         flat       (nx_block,ny_block,max_blocks), & ! latent heat flux   (W/m^2)
         fswabs     (nx_block,ny_block,max_blocks), & ! shortwave flux absorbed in ice and ocean (W/m^2)
         fswint_ai  (nx_block,ny_block,max_blocks), & ! SW absorbed in ice interior below surface (W/m^2)
         flwout     (nx_block,ny_block,max_blocks), & ! outgoing longwave radiation (W/m^2)
         Tref       (nx_block,ny_block,max_blocks), & ! 2m atm reference temperature (K)
         Qref       (nx_block,ny_block,max_blocks), & ! 2m atm reference spec humidity (kg/kg)
         Uref       (nx_block,ny_block,max_blocks), & ! 10m atm reference wind speed (m/s)
         evap       (nx_block,ny_block,max_blocks), & ! evaporative water flux (kg/m^2/s)
         evaps      (nx_block,ny_block,max_blocks), & ! evaporative water flux over snow (kg/m^2/s)
         evapi      (nx_block,ny_block,max_blocks), & ! evaporative water flux over ice (kg/m^2/s)
         alvdr      (nx_block,ny_block,max_blocks), & ! visible, direct   (fraction)
         alidr      (nx_block,ny_block,max_blocks), & ! near-ir, direct   (fraction)
         alvdf      (nx_block,ny_block,max_blocks), & ! visible, diffuse  (fraction)
         alidf      (nx_block,ny_block,max_blocks), & ! near-ir, diffuse  (fraction)
         alvdr_ai   (nx_block,ny_block,max_blocks), & ! visible, direct   (fraction)
         alidr_ai   (nx_block,ny_block,max_blocks), & ! near-ir, direct   (fraction)
         alvdf_ai   (nx_block,ny_block,max_blocks), & ! visible, diffuse  (fraction)
         alidf_ai   (nx_block,ny_block,max_blocks), & ! near-ir, diffuse  (fraction)
         albice     (nx_block,ny_block,max_blocks), & ! bare ice albedo
         albsno     (nx_block,ny_block,max_blocks), & ! snow albedo
         albpnd     (nx_block,ny_block,max_blocks), & ! melt pond albedo
         apeff_ai   (nx_block,ny_block,max_blocks), & ! effective pond area used for radiation calculation
         snowfrac   (nx_block,ny_block,max_blocks), & ! snow fraction used in radiation
         alvdr_init (nx_block,ny_block,max_blocks), & ! visible, direct   (fraction)
         alidr_init (nx_block,ny_block,max_blocks), & ! near-ir, direct   (fraction)
         alvdf_init (nx_block,ny_block,max_blocks), & ! visible, diffuse  (fraction)
         alidf_init (nx_block,ny_block,max_blocks), & ! near-ir, diffuse  (fraction)
         fpond      (nx_block,ny_block,max_blocks), & ! fresh water flux to ponds (kg/m^2/s)
         fresh      (nx_block,ny_block,max_blocks), & ! fresh water flux to ocean (kg/m^2/s)
         fsalt      (nx_block,ny_block,max_blocks), & ! salt flux to ocean (kg/m^2/s)
         fhocn      (nx_block,ny_block,max_blocks), & ! net heat flux to ocean (W/m^2)
         fswthru    (nx_block,ny_block,max_blocks), & ! shortwave penetrating to ocean (W/m^2)
         fswthru_vdr (nx_block,ny_block,max_blocks), & ! vis dir shortwave penetrating to ocean (W/m^2)
         fswthru_vdf (nx_block,ny_block,max_blocks), & ! vis dif shortwave penetrating to ocean (W/m^2)
         fswthru_idr (nx_block,ny_block,max_blocks), & ! nir dir shortwave penetrating to ocean (W/m^2)
         fswthru_idf (nx_block,ny_block,max_blocks), & ! nir dif shortwave penetrating to ocean (W/m^2)
         scale_factor (nx_block,ny_block,max_blocks), & ! scaling factor for shortwave components
         strairx_ocn(nx_block,ny_block,max_blocks), & ! stress on ocean by air, x-direction
         strairy_ocn(nx_block,ny_block,max_blocks), & ! stress on ocean by air, y-direction
         fsens_ocn  (nx_block,ny_block,max_blocks), & ! sensible heat flux (W/m^2)
         flat_ocn   (nx_block,ny_block,max_blocks), & ! latent heat flux   (W/m^2)
         flwout_ocn (nx_block,ny_block,max_blocks), & ! outgoing longwave radiation (W/m^2)
         evap_ocn   (nx_block,ny_block,max_blocks), & ! evaporative water flux (kg/m^2/s)
         alvdr_ocn  (nx_block,ny_block,max_blocks), & ! visible, direct   (fraction)
         alidr_ocn  (nx_block,ny_block,max_blocks), & ! near-ir, direct   (fraction)
         alvdf_ocn  (nx_block,ny_block,max_blocks), & ! visible, diffuse  (fraction)
         alidf_ocn  (nx_block,ny_block,max_blocks), & ! near-ir, diffuse  (fraction)
         Tref_ocn   (nx_block,ny_block,max_blocks), & ! 2m atm reference temperature (K)
         Qref_ocn   (nx_block,ny_block,max_blocks), & ! 2m atm reference spec humidity (kg/kg)
         fsurf      (nx_block,ny_block,max_blocks), & ! net surface heat flux (excluding fcondtop)(W/m^2)
         fcondtop   (nx_block,ny_block,max_blocks), & ! top surface conductive flux (W/m^2)
         fcondbot   (nx_block,ny_block,max_blocks), & ! bottom surface conductive flux (W/m^2)
         fbot       (nx_block,ny_block,max_blocks), & ! heat flux at bottom surface of ice (excluding excess) (W/m^2)
         Tbot       (nx_block,ny_block,max_blocks), & ! temperature at bottom surface of ice (deg C)
         Tsnice     (nx_block,ny_block,max_blocks), & ! temperature at snow ice interface (deg C)
         congel     (nx_block,ny_block,max_blocks), & ! basal ice growth         (m/step-->cm/day)
         frazil     (nx_block,ny_block,max_blocks), & ! frazil ice growth        (m/step-->cm/day)
         snoice     (nx_block,ny_block,max_blocks), & ! snow-ice formation       (m/step-->cm/day)
         meltt      (nx_block,ny_block,max_blocks), & ! top ice melt             (m/step-->cm/day)
         melts      (nx_block,ny_block,max_blocks), & ! snow melt                (m/step-->cm/day)
         meltb      (nx_block,ny_block,max_blocks), & ! basal ice melt           (m/step-->cm/day)
         meltl      (nx_block,ny_block,max_blocks), & ! lateral ice melt         (m/step-->cm/day)
         dsnow      (nx_block,ny_block,max_blocks), & ! change in snow thickness (m/step-->cm/day)
         daidtt     (nx_block,ny_block,max_blocks), & ! ice area tendency thermo.   (s^-1)
         dvidtt     (nx_block,ny_block,max_blocks), & ! ice volume tendency thermo. (m/s)
         dagedtt    (nx_block,ny_block,max_blocks), & ! ice age tendency thermo.    (s/s)
         mlt_onset  (nx_block,ny_block,max_blocks), & ! day of year that sfc melting begins
         frz_onset  (nx_block,ny_block,max_blocks), & ! day of year that freezing begins (congel or frazil)
         frazil_diag(nx_block,ny_block,max_blocks), & ! frazil ice growth diagnostic (m/step-->cm/day)
         fresh_ai   (nx_block,ny_block,max_blocks), & ! fresh water flux to ocean (kg/m^2/s)
         fsalt_ai   (nx_block,ny_block,max_blocks), & ! salt flux to ocean (kg/m^2/s)
         fhocn_ai   (nx_block,ny_block,max_blocks), & ! net heat flux to ocean (W/m^2)
         fswthru_ai (nx_block,ny_block,max_blocks), &  ! shortwave penetrating to ocean (W/m^2)
         fresh_da   (nx_block,ny_block,max_blocks), & ! fresh water flux to ocean due to data assim (kg/m^2/s)
         fsalt_da   (nx_block,ny_block,max_blocks), & ! salt flux to ocean due to data assimilation(kg/m^2/s)
         rside      (nx_block,ny_block,max_blocks), & ! fraction of ice that melts laterally
         fside      (nx_block,ny_block,max_blocks), & ! lateral melt rate (W/m^2)
         fsw        (nx_block,ny_block,max_blocks), & ! incoming shortwave radiation (W/m^2)
         coszen     (nx_block,ny_block,max_blocks), & ! cosine solar zenith angle, < 0 for sun below horizon 
         rdg_conv   (nx_block,ny_block,max_blocks), & ! convergence term for ridging (1/s)
         rdg_shear  (nx_block,ny_block,max_blocks), & ! shear term for ridging (1/s)
         dardg1ndt  (nx_block,ny_block,ncat,max_blocks), & ! rate of area loss by ridging ice (1/s)
         dardg2ndt  (nx_block,ny_block,ncat,max_blocks), & ! rate of area gain by new ridges (1/s)
         dvirdgndt  (nx_block,ny_block,ncat,max_blocks), & ! rate of ice volume ridged (m/s)
         aparticn   (nx_block,ny_block,ncat,max_blocks), & ! participation function
         krdgn      (nx_block,ny_block,ncat,max_blocks), & ! mean ridge thickness/thickness of ridging ice
         ardgn      (nx_block,ny_block,ncat,max_blocks), & ! fractional area of ridged ice
         vrdgn      (nx_block,ny_block,ncat,max_blocks), & ! volume of ridged ice
         araftn     (nx_block,ny_block,ncat,max_blocks), & ! rafting ice area
         vraftn     (nx_block,ny_block,ncat,max_blocks), & ! rafting ice volume 
         aredistn   (nx_block,ny_block,ncat,max_blocks), & ! redistribution function: fraction of new ridge area
         vredistn   (nx_block,ny_block,ncat,max_blocks), & ! redistribution function: fraction of new ridge volume
         fsurfn_f   (nx_block,ny_block,ncat,max_blocks), & ! net flux to top surface, excluding fcondtop
         fcondtopn_f(nx_block,ny_block,ncat,max_blocks), & ! downward cond flux at top surface (W m-2)
         fsensn_f   (nx_block,ny_block,ncat,max_blocks), & ! sensible heat flux (W m-2)
         flatn_f    (nx_block,ny_block,ncat,max_blocks), & ! latent heat flux (W m-2)
         meltsn     (nx_block,ny_block,ncat,max_blocks), & ! snow melt in category n (m)
         melttn     (nx_block,ny_block,ncat,max_blocks), & ! top melt in category n (m)
         meltbn     (nx_block,ny_block,ncat,max_blocks), & ! bottom melt in category n (m)
         congeln    (nx_block,ny_block,ncat,max_blocks), & ! congelation ice formation in category n (m)
         snoicen    (nx_block,ny_block,ncat,max_blocks), & ! snow-ice formation in category n (m)
         keffn_top  (nx_block,ny_block,ncat,max_blocks), & ! effective thermal conductivity of the top ice layer
         fsurfn     (nx_block,ny_block,ncat,max_blocks), & ! category fsurf
         fcondtopn  (nx_block,ny_block,ncat,max_blocks), & ! category fcondtop
         fcondbotn  (nx_block,ny_block,ncat,max_blocks), & ! category fcondbot
         fsensn     (nx_block,ny_block,ncat,max_blocks), & ! category sensible heat flux
         flatn      (nx_block,ny_block,ncat,max_blocks), & ! category latent heat flux
         albcnt     (nx_block,ny_block,max_blocks,max_nstrm), & ! counter for zenith angle
         salinz     (nx_block,ny_block,nilyr+1,max_blocks), & ! initial salinity  profile (ppt)   
         Tmltz      (nx_block,ny_block,nilyr+1,max_blocks), & ! initial melting temperature (^oC)
         stat=ierr)
      if (ierr/=0) call abort_ice('(alloc_flux): Out of memory')

      end subroutine alloc_flux

!=======================================================================

! Initialize all fluxes exchanged with flux coupler
! and some data-derived fields
!
! author Elizabeth C. Hunke, LANL

      subroutine init_coupler_flux

      use ice_arrays_column, only: Cdn_atm
      use ice_flux_bgc, only: flux_bio_atm, flux_bio, faero_atm, fiso_atm, &
           fnit, famm, fsil, fdmsp, fdms, fhum, fdust, falgalN, &
           fdoc, fdon, fdic, ffed, ffep
      use ice_grid, only: bathymetry

      integer (kind=int_kind) :: n

      integer (kind=int_kind), parameter :: max_d = 6
      real (kind=dbl_kind) :: fcondtopn_d(max_d), fsurfn_d(max_d)
      real (kind=dbl_kind) :: stefan_boltzmann, Tffresh
      real (kind=dbl_kind) :: vonkar, zref, iceruf

      integer :: i, j, iblk

      character(len=*), parameter :: subname = '(init_coupler_flux)'

      data fcondtopn_d / -50.0_dbl_kind,-17.0_dbl_kind,-12.0_dbl_kind, &
                          -9.0_dbl_kind, -7.0_dbl_kind, -3.0_dbl_kind /
      data fsurfn_d    /  0.20_dbl_kind, 0.15_dbl_kind, 0.10_dbl_kind, &
                          0.05_dbl_kind, 0.01_dbl_kind, 0.01_dbl_kind /

      call icepack_query_parameters(stefan_boltzmann_out=stefan_boltzmann, &
         Tffresh_out=Tffresh, vonkar_out=vonkar, zref_out=zref, iceruf_out=iceruf)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      !-----------------------------------------------------------------
      ! fluxes received from atmosphere
      !-----------------------------------------------------------------
      zlvl  (:,:,:) = c10             ! atm level height (m)
      rhoa  (:,:,:) = 1.3_dbl_kind    ! air density (kg/m^3)
      uatm  (:,:,:) = c5              ! wind velocity    (m/s)
      vatm  (:,:,:) = c5
      strax (:,:,:) = 0.05_dbl_kind
      stray (:,:,:) = 0.05_dbl_kind
      fsnow (:,:,:) = c0              ! snowfall rate (kg/m2/s)
                                      ! fsnow must be 0 for exact restarts
      if (trim(default_season) == 'winter') then
         ! typical winter values
         potT  (:,:,:) = 253.0_dbl_kind  ! air potential temp (K)
         Tair  (:,:,:) = 253.0_dbl_kind  ! air temperature  (K)
         Qa    (:,:,:) = 0.0006_dbl_kind ! specific humidity (kg/kg)
         swvdr (:,:,:) = c0              ! shortwave radiation (W/m^2)
         swvdf (:,:,:) = c0              ! shortwave radiation (W/m^2)
         swidr (:,:,:) = c0              ! shortwave radiation (W/m^2)
         swidf (:,:,:) = c0              ! shortwave radiation (W/m^2)
         flw   (:,:,:) = c180            ! incoming longwave rad (W/m^2)
         frain (:,:,:) = c0              ! rainfall rate (kg/m2/s)
         do n = 1, ncat              ! conductive heat flux (W/m^2)
            fcondtopn_f(:,:,n,:) = fcondtopn_d(min(n,max_d))
         enddo
         fsurfn_f = fcondtopn_f      ! surface heat flux (W/m^2)
         flatn_f (:,:,:,:) = c0          ! latent heat flux (kg/m2/s)
         fsensn_f(:,:,:,:) = c0          ! sensible heat flux (W/m^2)
      elseif (trim(default_season) == 'summer') then
         ! typical summer values
         potT  (:,:,:) = 273.0_dbl_kind  ! air potential temp (K)
         Tair  (:,:,:) = 273.0_dbl_kind  ! air temperature  (K)
         Qa    (:,:,:) = 0.0035_dbl_kind ! specific humidity (kg/kg)
         swvdr (:,:,:) = 50._dbl_kind    ! shortwave radiation (W/m^2)
         swvdf (:,:,:) = 50._dbl_kind    ! shortwave radiation (W/m^2)
         swidr (:,:,:) = 50._dbl_kind    ! shortwave radiation (W/m^2)
         swidf (:,:,:) = 50._dbl_kind    ! shortwave radiation (W/m^2)
         flw   (:,:,:) = 280.0_dbl_kind  ! incoming longwave rad (W/m^2)
         frain (:,:,:) = c0              ! rainfall rate (kg/m2/s)
         do n = 1, ncat                   ! surface heat flux (W/m^2)
            fsurfn_f(:,:,n,:) = fsurfn_d(min(n,max_d))
         enddo
         fcondtopn_f(:,:,:,:) =  0.0_dbl_kind ! conductive heat flux (W/m^2)
         flatn_f    (:,:,:,:) = -2.0_dbl_kind ! latent heat flux (W/m^2)
         fsensn_f   (:,:,:,:) =  c0           ! sensible heat flux (W/m^2)
      else
         ! typical spring values
         potT  (:,:,:) = 263.15_dbl_kind ! air potential temp (K)
         Tair  (:,:,:) = 263.15_dbl_kind ! air temperature  (K)
         Qa    (:,:,:) = 0.001_dbl_kind  ! specific humidity (kg/kg)
         swvdr (:,:,:) = 25._dbl_kind    ! shortwave radiation (W/m^2)
         swvdf (:,:,:) = 25._dbl_kind    ! shortwave radiation (W/m^2)
         swidr (:,:,:) = 25._dbl_kind    ! shortwave radiation (W/m^2)
         swidf (:,:,:) = 25._dbl_kind    ! shortwave radiation (W/m^2)
         flw   (:,:,:) = 230.0_dbl_kind  ! incoming longwave rad (W/m^2)
         frain (:,:,:) = c0              ! rainfall rate (kg/m2/s)
         do n = 1, ncat                   ! surface heat flux (W/m^2)
            fsurfn_f(:,:,n,:) = fsurfn_d(min(n,max_d))
         enddo
         fcondtopn_f(:,:,:,:) =  c0           ! conductive heat flux (W/m^2)
         flatn_f    (:,:,:,:) = -1.0_dbl_kind ! latent heat flux (W/m^2)
         fsensn_f   (:,:,:,:) =  c0           ! sensible heat flux (W/m^2)
      endif !   

      fiso_atm  (:,:,:,:) = c0           ! isotope deposition rate (kg/m2/s)
      faero_atm (:,:,:,:) = c0           ! aerosol deposition rate (kg/m2/s)
      flux_bio_atm (:,:,:,:) = c0        ! zaero and bio deposition rate (kg/m2/s)

      !-----------------------------------------------------------------
      ! fluxes received from ocean
      !-----------------------------------------------------------------

      ss_tltx(:,:,:)= c0              ! sea surface tilt (m/m)
      ss_tlty(:,:,:)= c0
      uocn  (:,:,:) = c0              ! surface ocean currents (m/s)
      vocn  (:,:,:) = c0
      frzmlt(:,:,:) = c0              ! freezing/melting potential (W/m^2)
      frzmlt_init(:,:,:) = c0         ! freezing/melting potential (W/m^2)
      sss   (:,:,:) = 34.0_dbl_kind   ! sea surface salinity (ppt)

      do iblk = 1, size(Tf,3)
      do j = 1, size(Tf,2)
      do i = 1, size(Tf,1)
         Tf (i,j,iblk) = icepack_liquidus_temperature(sss(i,j,iblk)) ! freezing temp (C)
      enddo
      enddo
      enddo

      sst   (:,:,:) = Tf(:,:,:)       ! sea surface temp (C)
      qdp   (:,:,:) = c0              ! deep ocean heat flux (W/m^2)
      hmix  (:,:,:) = c20             ! ocean mixed layer depth (m)
      hwater(:,:,:) = bathymetry(:,:,:) ! ocean water depth (m)
      daice_da(:,:,:) = c0            ! data assimilation increment rate

      !-----------------------------------------------------------------
      ! fluxes sent to atmosphere
      !-----------------------------------------------------------------

      strairxT(:,:,:) = c0            ! wind stress, T grid
      strairyT(:,:,:) = c0

      fsens   (:,:,:) = c0
      flat    (:,:,:) = c0
      fswabs  (:,:,:) = c0
      fswint_ai(:,:,:) = c0
      flwout  (:,:,:) = -stefan_boltzmann*Tffresh**4   
                        ! in case atm model diagnoses Tsfc from flwout
      evap    (:,:,:) = c0
      evaps   (:,:,:) = c0
      evapi   (:,:,:) = c0
      Tref    (:,:,:) = c0
      Qref    (:,:,:) = c0
      Uref    (:,:,:) = c0
      alvdr   (:,:,:) = c0
      alidr   (:,:,:) = c0
      alvdf   (:,:,:) = c0
      alidf   (:,:,:) = c0

      !-----------------------------------------------------------------
      ! fluxes sent to ocean
      !-----------------------------------------------------------------

      strocnxT(:,:,:) = c0    ! ice-ocean stress, x-direction (T-cell)
      strocnyT(:,:,:) = c0    ! ice-ocean stress, y-direction (T-cell)
      fresh   (:,:,:) = c0
      fsalt   (:,:,:) = c0
      fpond   (:,:,:) = c0
      fhocn   (:,:,:) = c0
      fswthru (:,:,:) = c0
      fswthru_vdr (:,:,:) = c0
      fswthru_vdf (:,:,:) = c0
      fswthru_idr (:,:,:) = c0
      fswthru_idf (:,:,:) = c0
      fresh_da(:,:,:) = c0    ! data assimilation
      fsalt_da(:,:,:) = c0
      flux_bio (:,:,:,:) = c0 ! bgc
      fnit    (:,:,:) = c0
      fsil    (:,:,:) = c0
      famm    (:,:,:) = c0
      fdmsp   (:,:,:) = c0
      fdms    (:,:,:) = c0
      fhum    (:,:,:) = c0
      fdust   (:,:,:) = c0
      falgalN(:,:,:,:)= c0
      fdoc   (:,:,:,:)= c0
      fdic   (:,:,:,:)= c0
      fdon   (:,:,:,:)= c0
      ffep   (:,:,:,:)= c0
      ffed   (:,:,:,:)= c0
      
      if (send_i2x_per_cat) then
         allocate(fswthrun_ai(nx_block,ny_block,ncat,max_blocks))
         fswthrun_ai(:,:,:,:) = c0
      endif

      !-----------------------------------------------------------------
      ! derived or computed fields
      !-----------------------------------------------------------------

      coszen  (:,:,:) = c0            ! Cosine of the zenith angle
      fsw     (:,:,:) = c0            ! shortwave radiation (W/m^2)
      scale_factor(:,:,:) = c1        ! shortwave scaling factor 
      wind    (:,:,:) = sqrt(uatm(:,:,:)**2 &
                           + vatm(:,:,:)**2)  ! wind speed, (m/s)
      Cdn_atm(:,:,:) = (vonkar/log(zref/iceruf)) &
                     * (vonkar/log(zref/iceruf)) ! atmo drag for RASM
      alvdr_init(:,:,:) = c0
      alidr_init(:,:,:) = c0
      alvdf_init(:,:,:) = c0
      alidf_init(:,:,:) = c0


      end subroutine init_coupler_flux

!=======================================================================

! Initialize some fluxes sent to coupler for use by the atm model
!
! author: Elizabeth C. Hunke, LANL

      subroutine init_flux_atm

      use ice_flux_bgc, only: fiso_evap, Qref_iso, Qa_iso

      character(len=*), parameter :: subname = '(init_flux_atm)'

      !-----------------------------------------------------------------
      ! initialize albedo and fluxes
      !-----------------------------------------------------------------

      strairxT(:,:,:) = c0      ! wind stress, T grid
      strairyT(:,:,:) = c0
      ! for rectangular grid tests without thermo
      ! strairxT(:,:,:) = 0.15_dbl_kind
      ! strairyT(:,:,:) = 0.15_dbl_kind

      fsens   (:,:,:) = c0
      flat    (:,:,:) = c0
      fswabs  (:,:,:) = c0
      flwout  (:,:,:) = c0
      evap    (:,:,:) = c0
      Tref    (:,:,:) = c0
      Qref    (:,:,:) = c0
      Uref    (:,:,:) = c0

      fiso_evap(:,:,:,:) = c0
      Qref_iso (:,:,:,:) = c0
      Qa_iso   (:,:,:,:) = c0

      end subroutine init_flux_atm

!=======================================================================
! Initialize some fluxes sent to coupler for use by the ocean model
!
! NOTE: These fluxes should be initialized immediately after the
!       call to the coupler.  The atmospheric fluxes can be initialized
!       at the beginning of the following time step because they are
!       not modified by any subroutines between the call to_coupler
!       and the end of the time step.
!
! author: Elizabeth C. Hunke, LANL

      subroutine init_flux_ocn

      use ice_flux_bgc, only: faero_ocn, fiso_ocn, HDO_ocn, H2_16O_ocn, H2_18O_ocn

      character(len=*), parameter :: subname = '(init_flux_ocn)'

      !-----------------------------------------------------------------
      ! fluxes sent
      !-----------------------------------------------------------------

      fresh    (:,:,:)   = c0
      fsalt    (:,:,:)   = c0
      fpond    (:,:,:)   = c0
      fhocn    (:,:,:)   = c0
      fswthru  (:,:,:)   = c0
      fswthru_vdr  (:,:,:)   = c0
      fswthru_vdf  (:,:,:)   = c0
      fswthru_idr  (:,:,:)   = c0
      fswthru_idf  (:,:,:)   = c0

      faero_ocn (:,:,:,:) = c0
      fiso_ocn  (:,:,:,:) = c0
      HDO_ocn     (:,:,:) = c0
      H2_16O_ocn  (:,:,:) = c0
      H2_18O_ocn  (:,:,:) = c0

      if (send_i2x_per_cat) then
         fswthrun_ai(:,:,:,:) = c0
      endif

      end subroutine init_flux_ocn

!=======================================================================

! Initialize thermodynamic fields written to history files.
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL

      subroutine init_history_therm

      use ice_state, only: aice, vice, trcr
      use ice_arrays_column, only: &
          hfreebd, hdraft, hridge, distrdg, hkeel, dkeel, lfloe, dfloe, &
          Cdn_atm_skin, Cdn_atm_floe, Cdn_atm_pond, Cdn_atm_rdg, &
          Cdn_ocn_skin, Cdn_ocn_floe, Cdn_ocn_keel, Cdn_atm_ratio, &
          Cdn_atm, Cdn_ocn

      logical (kind=log_kind) :: &
          formdrag, &
          tr_iage

      integer (kind=int_kind) :: &
          nt_iage

      real (kind=dbl_kind) :: &
          dragio, &
          vonkar, &
          zref, &
          iceruf

      character(len=*), parameter :: subname = '(init_history_therm)'

      call icepack_query_parameters(formdrag_out=formdrag)
      call icepack_query_tracer_flags(tr_iage_out=tr_iage)
      call icepack_query_tracer_indices(nt_iage_out=nt_iage)
      call icepack_query_parameters( dragio_out=dragio, &
         vonkar_out=vonkar, zref_out=zref, iceruf_out=iceruf)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      fsurf  (:,:,:) = c0
      fcondtop(:,:,:)= c0
      fcondbot(:,:,:)= c0
      congel (:,:,:) = c0
      fbot   (:,:,:) = c0
      Tbot   (:,:,:) = c0
      Tsnice  (:,:,:) = c0
      frazil (:,:,:) = c0
      snoice (:,:,:) = c0
      dsnow  (:,:,:) = c0
      meltt  (:,:,:) = c0
      melts  (:,:,:) = c0
      meltb  (:,:,:) = c0
      meltl  (:,:,:) = c0
      daidtt (:,:,:) = aice(:,:,:) ! temporary initial area
      dvidtt (:,:,:) = vice(:,:,:) ! temporary initial volume
      if (tr_iage) then
         dagedtt(:,:,:) = trcr(:,:,nt_iage,:) ! temporary initial age
      else
         dagedtt(:,:,:) = c0
      endif
      fsurfn    (:,:,:,:) = c0
      fcondtopn (:,:,:,:) = c0
      fcondbotn (:,:,:,:) = c0
      flatn     (:,:,:,:) = c0
      fsensn    (:,:,:,:) = c0
      fresh_ai  (:,:,:) = c0
      fsalt_ai  (:,:,:) = c0
      fhocn_ai  (:,:,:) = c0
      fswthru_ai(:,:,:) = c0
      albice (:,:,:) = c0
      albsno (:,:,:) = c0
      albpnd (:,:,:) = c0
      apeff_ai (:,:,:) = c0
      snowfrac (:,:,:) = c0
      frazil_diag (:,:,:) = c0

      ! drag coefficients are computed prior to the atmo_boundary call, 
      ! during the thermodynamics section 
      Cdn_ocn(:,:,:) = dragio
      Cdn_atm(:,:,:) = (vonkar/log(zref/iceruf)) &
                     * (vonkar/log(zref/iceruf)) ! atmo drag for RASM
      Cdn_atm_ratio(:,:,:)= c0

      if (formdrag) then
        Cdn_atm_rdg (:,:,:) = c0
        Cdn_atm_floe(:,:,:) = c0
        Cdn_atm_pond(:,:,:) = c0
        Cdn_atm_skin(:,:,:) = c0
        Cdn_ocn_skin(:,:,:) = c0
        Cdn_ocn_keel(:,:,:) = c0
        Cdn_ocn_floe(:,:,:) = c0
        hfreebd     (:,:,:) = c0
        hdraft      (:,:,:) = c0
        hridge      (:,:,:) = c0
        distrdg     (:,:,:) = c0
        hkeel       (:,:,:) = c0
        dkeel       (:,:,:) = c0
        lfloe       (:,:,:) = c0
        dfloe       (:,:,:) = c0
      endif

      end subroutine init_history_therm

!=======================================================================

! Initialize dynamic fields written to history files.
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL

      subroutine init_history_dyn

      use ice_state, only: aice, vice, trcr, strength

      logical (kind=log_kind) :: &
          tr_iage

      integer (kind=int_kind) :: &
          nt_iage

      character(len=*), parameter :: subname = '(init_history_dyn)'

      call icepack_query_tracer_flags(tr_iage_out=tr_iage)
      call icepack_query_tracer_indices(nt_iage_out=nt_iage)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      sig1    (:,:,:) = c0
      sig2    (:,:,:) = c0
      taubx   (:,:,:) = c0
      tauby   (:,:,:) = c0
      strength (:,:,:) = c0
      strocnx (:,:,:) = c0
      strocny (:,:,:) = c0
      strairx (:,:,:) = c0
      strairy (:,:,:) = c0
      strtltx (:,:,:) = c0
      strtlty (:,:,:) = c0
      strintx (:,:,:) = c0
      strinty (:,:,:) = c0
      dardg1dt(:,:,:) = c0
      dardg2dt(:,:,:) = c0
      dvirdgdt(:,:,:) = c0
      opening (:,:,:) = c0
      daidtd  (:,:,:) = aice(:,:,:) ! temporary initial area
      dvidtd  (:,:,:) = vice(:,:,:) ! temporary initial volume
      if (tr_iage) &
         dagedtd (:,:,:) = trcr(:,:,nt_iage,:) ! temporary initial age
      fm      (:,:,:) = c0
      ardgn   (:,:,:,:) = c0
      vrdgn   (:,:,:,:) = c0
      krdgn   (:,:,:,:) = c1
      aparticn(:,:,:,:) = c0
      aredistn(:,:,:,:) = c0
      vredistn(:,:,:,:) = c0
      dardg1ndt(:,:,:,:) = c0
      dardg2ndt(:,:,:,:) = c0
      dvirdgndt(:,:,:,:) = c0
      araftn   (:,:,:,:) = c0
      vraftn   (:,:,:,:) = c0
      aredistn (:,:,:,:) = c0
      vredistn (:,:,:,:) = c0

      end subroutine init_history_dyn

!=======================================================================

!  Divide ice fluxes by ice area before sending them to the
!  coupler, since the coupler multiplies by ice area.
!
! authors: C.M.Bitz, William H. Lipscomb

      subroutine scale_fluxes (nx_block, ny_block, &
                               tmask,              &
                               nbtrcr,   max_aero, &
                               aice,     Tf,       &
                               Tair,     Qa,       &
                               strairxT, strairyT, &
                               fsens,    flat,     &
                               fswabs,   flwout,   &
                               evap,               &
                               Tref,     Qref,     &
                               fresh,    fsalt,    &
                               fhocn,    fswthru,  &
                               fswthru_vdr, fswthru_vdf, &
                               fswthru_idr, fswthru_idf, &
                               faero_ocn,          &
                               alvdr,    alidr,    &
                               alvdf,    alidf,    &
                               fzsal,    fzsal_g,  &
                               flux_bio,           &
                               fsurf,    fcondtop, &
                               Uref,     wind,     &
                               Qref_iso,           &
                               fiso_evap,fiso_ocn)

      use icepack_intfc, only: icepack_max_iso

      integer (kind=int_kind), intent(in) :: &
          nx_block, ny_block, &    ! block dimensions
          nbtrcr            , &    ! number of biology tracers
          max_aero                 ! maximum number of aerosols

      logical (kind=log_kind), dimension (nx_block,ny_block), intent(in) :: &
          tmask     ! land/boundary mask, thickness (T-cell)


      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) :: &
          aice    , & ! fractional ice area
          Tf      , & ! freezing temperature            (C)
          Tair    , & ! surface air temperature         (K)
          Qa          ! sfc air specific humidity       (kg/kg)

      real (kind=dbl_kind), dimension(nx_block,ny_block), optional, intent(in) :: &
          wind        ! wind speed                      (m/s)

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(inout) :: &
          strairxT, & ! air/ice zonal  stress           (N/m**2)
          strairyT, & ! air/ice merdnl stress           (N/m**2)
          fsens   , & ! sensible heat flx               (W/m**2)
          flat    , & ! latent   heat flx               (W/m**2)
          fswabs  , & ! shortwave absorbed heat flx     (W/m**2)
          flwout  , & ! upwd lw emitted heat flx        (W/m**2)
          evap    , & ! evaporation                     (kg/m2/s)
          Tref    , & ! air tmp reference level         (K)
          Qref    , & ! air sp hum reference level      (kg/kg)
          fresh   , & ! fresh water flux to ocean       (kg/m2/s)
          fsalt   , & ! salt flux to ocean              (kg/m2/s)
          fhocn   , & ! actual ocn/ice heat flx         (W/m**2)
          fswthru , & ! sw radiation through ice bot    (W/m**2)
          fswthru_vdr , & ! vis dir sw radiation through ice bot    (W/m**2)
          fswthru_vdf , & ! vis dif sw radiation through ice bot    (W/m**2)
          fswthru_idr , & ! nir dir sw radiation through ice bot    (W/m**2)
          fswthru_idf , & ! nir dif sw radiation through ice bot    (W/m**2)
          alvdr   , & ! visible, direct   (fraction)
          alidr   , & ! near-ir, direct   (fraction)
          alvdf   , & ! visible, diffuse  (fraction)
          alidf       ! near-ir, diffuse  (fraction)

      real (kind=dbl_kind), dimension(nx_block,ny_block), optional, intent(inout) :: &
          Uref        ! air speed reference level       (m/s)

      real (kind=dbl_kind), dimension(nx_block,ny_block,nbtrcr), intent(inout) :: &
          flux_bio    ! tracer flux to ocean from biology (mmol/m2/s)

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_aero), intent(inout) :: &
          faero_ocn   ! aerosol flux to ocean            (kg/m2/s)

      ! For hadgem drivers. Assumes either both fields are passed or neither
      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(inout), optional :: &
          fsurf   , & ! surface heat flux               (W/m**2)
          fcondtop    ! top surface conductive flux     (W/m**2)

      ! zsalinity fluxes
      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(inout) :: &
          fzsal   , & ! salt flux to ocean with prognositic salinity (kg/m2/s)  
          fzsal_g     ! Gravity drainage salt flux to ocean (kg/m2/s) 

      ! isotopes
      real (kind=dbl_kind), dimension(nx_block,ny_block,icepack_max_iso), &
          optional, intent(inout) :: &
          Qref_iso , & ! isotope air sp hum reference level      (kg/kg)
          fiso_evap, & ! isotope evaporation (kg/m2/s)
          fiso_ocn     ! isotope flux to ocean (kg/m2/s)

      ! local variables

      real (kind=dbl_kind) :: &
          ar, &   ! 1/aice
          stefan_boltzmann, &
          Tffresh

      integer (kind=int_kind) :: &
          i, j    ! horizontal indices

      character(len=*), parameter :: subname = '(scale_fluxes)'

      call icepack_query_parameters(stefan_boltzmann_out=stefan_boltzmann, &
         Tffresh_out=Tffresh)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      do j = 1, ny_block
      do i = 1, nx_block
         if (tmask(i,j) .and. aice(i,j) > c0) then
            ar = c1 / aice(i,j)
            strairxT(i,j) = strairxT(i,j) * ar
            strairyT(i,j) = strairyT(i,j) * ar
            fsens   (i,j) = fsens   (i,j) * ar
            flat    (i,j) = flat    (i,j) * ar
            fswabs  (i,j) = fswabs  (i,j) * ar
            flwout  (i,j) = flwout  (i,j) * ar
            evap    (i,j) = evap    (i,j) * ar
            Tref    (i,j) = Tref    (i,j) * ar
            Qref    (i,j) = Qref    (i,j) * ar
            if (present(Uref)) &
               Uref (i,j) = Uref    (i,j) * ar
            fresh   (i,j) = fresh   (i,j) * ar
            fsalt   (i,j) = fsalt   (i,j) * ar
            fhocn   (i,j) = fhocn   (i,j) * ar
            fswthru (i,j) = fswthru (i,j) * ar
            fswthru_vdr (i,j) = fswthru_vdr (i,j) * ar
            fswthru_vdf (i,j) = fswthru_vdf (i,j) * ar
            fswthru_idr (i,j) = fswthru_idr (i,j) * ar
            fswthru_idf (i,j) = fswthru_idf (i,j) * ar
            alvdr   (i,j) = alvdr   (i,j) * ar
            alidr   (i,j) = alidr   (i,j) * ar
            alvdf   (i,j) = alvdf   (i,j) * ar
            alidf   (i,j) = alidf   (i,j) * ar
            fzsal   (i,j) = fzsal   (i,j) * ar  
            fzsal_g (i,j) = fzsal_g (i,j) * ar  
            flux_bio (i,j,:) = flux_bio (i,j,:) * ar
            faero_ocn(i,j,:) = faero_ocn(i,j,:) * ar
            if (present(Qref_iso )) Qref_iso (i,j,:) = Qref_iso (i,j,:) * ar
            if (present(fiso_evap)) fiso_evap(i,j,:) = fiso_evap(i,j,:) * ar
            if (present(fiso_ocn )) fiso_ocn (i,j,:) = fiso_ocn (i,j,:) * ar
         else                   ! zero out fluxes
            strairxT(i,j) = c0
            strairyT(i,j) = c0
            fsens   (i,j) = c0
            flat    (i,j) = c0
            fswabs  (i,j) = c0
            flwout  (i,j) = -stefan_boltzmann *(Tf(i,j) + Tffresh)**4
               ! to make upward longwave over ocean reasonable for history file
            evap    (i,j) = c0
            Tref    (i,j) = Tair(i,j)
            Qref    (i,j) = Qa  (i,j)
            if (present(Uref) .and. present(wind)) &
               Uref (i,j) = wind(i,j)
            fresh   (i,j) = c0
            fsalt   (i,j) = c0
            fhocn   (i,j) = c0
            fswthru (i,j) = c0
            fswthru_vdr (i,j) = c0
            fswthru_vdf (i,j) = c0
            fswthru_idr (i,j) = c0
            fswthru_idf (i,j) = c0
            alvdr   (i,j) = c0  ! zero out albedo where ice is absent
            alidr   (i,j) = c0
            alvdf   (i,j) = c0 
            alidf   (i,j) = c0
            fzsal   (i,j) = c0  
            fzsal_g (i,j) = c0 
            flux_bio (i,j,:) = c0
            faero_ocn(i,j,:) = c0
            if (present(Qref_iso )) Qref_iso (i,j,:) = c0
            if (present(fiso_evap)) fiso_evap(i,j,:) = c0
            if (present(fiso_ocn )) fiso_ocn (i,j,:) = c0
         endif                  ! tmask and aice > 0
      enddo                     ! i
      enddo                     ! j

      ! Scale fluxes for history output
      if (present(fsurf) .and. present(fcondtop) ) then 
     
        do j = 1, ny_block
        do i = 1, nx_block
           if (tmask(i,j) .and. aice(i,j) > c0) then
              ar = c1 / aice(i,j)
              fsurf   (i,j) = fsurf   (i,j) * ar
              fcondtop(i,j) = fcondtop(i,j) * ar
           else                   ! zero out fluxes
              fsurf   (i,j) = c0
              fcondtop(i,j) = c0
           endif                  ! tmask and aice > 0
        enddo                     ! i
        enddo                     ! j
      
      endif                       ! present(fsurf & fcondtop)
      
      end subroutine scale_fluxes

!=======================================================================

      end module ice_flux

!=======================================================================
