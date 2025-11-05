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
      ! All variables are assumed to be on the atm or ocn thermodynamic
      ! grid except as noted
      !
      ! scale_fluxes divides several of these by aice "in place", so
      ! the state of some of these variables is not well defined.  In the
      ! future, we need to refactor and add "_iavg" versions of the
      ! fields to clearly differentiate fields that have been divided
      ! by aice and others that are not.  The challenge is that we need
      ! to go thru each field carefully to see which version is used.
      ! For instance, in diagnostics, there are places where these
      ! fields are multiplied by aice to compute things properly.
      ! strocn[x,y]T_iavg is the first field defined using _iavg.
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &

       ! in from atmos (if .not.calc_strair)
         strax   , & ! wind stress components (N/m^2), on grid_atm_dynu
         stray   , & !                                 on grid_atm_dynv

       ! in from ocean
         uocn    , & ! ocean current, x-direction (m/s),     on grid_ocn_dynu
         vocn    , & ! ocean current, y-direction (m/s),     on grid_ocn_dynv
         ss_tltx , & ! sea surface slope, x-direction (m/m), on grid_ocn_dynu
         ss_tlty , & ! sea surface slope, y-direction,       on grid_ocn_dynv
         hwater  , & ! water depth for seabed stress calc (landfast ice)

       ! out to atmosphere
         strairxT, & ! stress on ice by air, x-direction at T points, computed in icepack
         strairyT, & ! stress on ice by air, y-direction at T points, computed in icepack

       ! out to ocean          T-cell (kg/m s^2)
       ! Note, CICE_IN_NEMO uses strocnx and strocny for coupling
         strocnxT_iavg, & ! ice-ocean stress, x-direction at T points, per ice fraction (scaled flux)
         strocnyT_iavg    ! ice-ocean stress, y-direction at T points, per ice fraction (scaled flux)

       ! diagnostic

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         sig1    , & ! normalized principal stress component
         sig2    , & ! normalized principal stress component
         sigP    , & ! internal ice pressure (N/m)
         taubxU  , & ! seabed stress (x) (N/m^2)
         taubyU  , & ! seabed stress (y) (N/m^2)
         strairxU, & ! stress on ice by air, x-direction at U points
         strairyU, & ! stress on ice by air, y-direction at U points
         strocnxU, & ! ice-ocean stress, x-direction at U points, computed in dyn_finish
         strocnyU, & ! ice-ocean stress, y-direction at U points, computed in dyn_finish
         strtltxU, & ! stress due to sea surface slope, x-direction
         strtltyU, & ! stress due to sea surface slope, y-direction
         strintxU, & ! divergence of internal ice stress, x (N/m^2)
         strintyU, & ! divergence of internal ice stress, y (N/m^2)
         taubxN  , & ! seabed stress (x) at N points (N/m^2)
         taubyN  , & ! seabed stress (y) at N points (N/m^2)
         strairxN, & ! stress on ice by air, x-direction at N points
         strairyN, & ! stress on ice by air, y-direction at N points
         strocnxN, & ! ice-ocean stress, x-direction at N points, computed in dyn_finish
         strocnyN, & ! ice-ocean stress, y-direction at N points, computed in dyn_finish
         strtltxN, & ! stress due to sea surface slope, x-direction at N points
         strtltyN, & ! stress due to sea surface slope, y-direction at N points
         strintxN, & ! divergence of internal ice stress, x at N points (N/m^2)
         strintyN, & ! divergence of internal ice stress, y at N points (N/m^2)
         taubxE  , & ! seabed stress (x) at E points (N/m^2)
         taubyE  , & ! seabed stress (y) at E points (N/m^2)
         strairxE, & ! stress on ice by air, x-direction at E points
         strairyE, & ! stress on ice by air, y-direction at E points
         strocnxE, & ! ice-ocean stress, x-direction at E points, computed in dyn_finish
         strocnyE, & ! ice-ocean stress, y-direction at E points, computed in dyn_finish
         strtltxE, & ! stress due to sea surface slope, x-direction at E points
         strtltyE, & ! stress due to sea surface slope, y-direction at E points
         strintxE, & ! divergence of internal ice stress, x at E points (N/m^2)
         strintyE, & ! divergence of internal ice stress, y at E points (N/m^2)
         daidtd  , & ! ice area tendency due to transport   (1/s)
         dvidtd  , & ! ice volume tendency due to transport (m/s)
         dvsdtd  , & ! snow volume tendency due to transport (m/s)
         dagedtd , & ! ice age tendency due to transport (s/s)
         dardg1dt, & ! rate of area loss by ridging ice (1/s)
         dardg2dt, & ! rate of area gain by new ridges (1/s)
         dvirdgdt, & ! rate of ice volume ridged (m/s)
         opening     ! rate of opening due to divergence/shear (1/s)

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
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
         stress12_1,stress12_2,stress12_3,stress12_4, & ! sigma12
       ! ice stress tensor at U and T locations (grid_ice = 'C|CD') (kg/s^2)
         stresspT, stressmT, stress12T, & ! sigma11+sigma22, sigma11-sigma22, sigma12
         stresspU, stressmU, stress12U    ! "

       ! internal

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         fmU      , & ! Coriolis param. * mass in U-cell (kg/s)
         TbU      , & ! factor for seabed stress (N/m^2)
         fmE      , & ! Coriolis param. * mass in E-cell (kg/s)
         TbE      , & ! factor for seabed stress (N/m^2)
         fmN      , & ! Coriolis param. * mass in N-cell (kg/s)
         TbN          ! factor for seabed stress (N/m^2)

      !-----------------------------------------------------------------
      ! Thermodynamic component
      !-----------------------------------------------------------------

       ! in from atmosphere (if calc_Tsfc)

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         zlvl    , & ! atm level height (momentum) (m)
         zlvs    , & ! atm level height (scalar quantities) (m)
         uatm    , & ! wind velocity components (m/s), on grid_atm_dynu
         vatm    , & !                                 on grid_atm_dynv
         wind    , & ! wind speed (m/s)              , on grid_atm_dynu
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

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         fsurfn_f   , & ! net flux to top surface, excluding fcondtop
         fcondtopn_f, & ! downward cond flux at top surface (W m-2)
         fsensn_f   , & ! sensible heat flux (W m-2)
         flatn_f        ! latent heat flux (W m-2)

      ! in from atmosphere
      ! required for coupling in GEOS

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         evapn_f,      & ! evaporation/sublimation (kg m-2 s-1)
         dflatndTsfc_f,  & ! derivative of latent flux w.r.t. Tsfc
         dfsurfndTsfc_f    ! derivative of surface flux w.r.t. Tsfc

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         swuvrdr     , & !  vis uvr flux, direct (W m-2)
         swuvrdf     , & !  vis uvr flux, diffuse (W m-2)
         swpardr     , & !  vis par flux, direct (W m-2)
         swpardf         !  vis par flux, diffuse (W m-2)

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

      real (kind=dbl_kind), dimension(:,:,:,:), allocatable, public :: &
         albcnt       ! counter for zenith angle

       ! out to ocean
       ! (Note CICE_IN_NEMO does not use these for coupling.
       !  It uses fresh_ai,fsalt_ai,fhocn_ai and fswthru_ai)
      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         fpond   , & ! fresh water flux to ponds (kg/m^2/s)
         fresh   , & ! fresh water flux to ocean (kg/m^2/s)
         fsalt   , & ! salt flux to ocean (kg/m^2/s)
         fhocn   , & ! net heat flux to ocean (W/m^2)
         fsloss  , & ! rate of snow loss to leads (kg/m^2/s)
         fswthru , & ! shortwave penetrating to ocean (W/m^2)
         fswthru_vdr , & ! vis dir shortwave penetrating to ocean (W/m^2)
         fswthru_vdf , & ! vis dif shortwave penetrating to ocean (W/m^2)
         fswthru_idr , & ! nir dir shortwave penetrating to ocean (W/m^2)
         fswthru_idf , & ! nir dif shortwave penetrating to ocean (W/m^2)
         fswthru_uvrdr,& ! vis dir uvr SW penetrating to ocean (W/m^2)
         fswthru_uvrdf,& ! vis dif uvr SW penetrating to ocean (W/m^2)
         fswthru_pardr,& ! nir dir par SW penetrating to ocean (W/m^2)
         fswthru_pardf   ! nir dif par SW penetrating to ocean (W/m^2)

       ! internal

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         scale_factor! scaling factor for shortwave components

      logical (kind=log_kind), public :: &
         update_ocn_f, & ! if true, update fresh water and salt fluxes
         l_mpond_fresh   ! if true, include freshwater feedback from meltponds
                         ! when running in ice-ocean or coupled configuration

      character (char_len), public :: &
         cpl_frazil      ! type of coupling for frazil ice, 'fresh_ice_correction','internal','external'

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
         dvsdtt, & ! snow volume tendency thermo. (m/s)
         dagedtt,& ! ice age tendency thermo.    (s/s)
         mlt_onset, &! day of year that sfc melting begins
         frz_onset, &! day of year that freezing begins (congel or frazil)
         frazil_diag ! frazil ice growth diagnostic (m/step-->cm/day)

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         dpnd_flush,  & ! pond flushing rate due to ice permeability (m/step)
         dpnd_expon,  & ! exponential pond drainage rate (m/step)
         dpnd_freebd, & ! pond drainage rate due to freeboard constraint (m/step)
         dpnd_initial,& ! runoff rate due to rfrac (m/step)
         dpnd_dlid,   & ! pond loss/gain (+/-) to ice lid freezing/melting (m/step)
         dpnd_melt,   & ! pond 'drainage' due to ice melting (m / step)
         dpnd_ridge     ! pond 'drainage' due to ridging (m)

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         dpnd_flushn, & ! category pond flushing rate due to ice permeability (m/step)
         dpnd_exponn, & ! category exponential pond drainage rate (m/step)
         dpnd_freebdn,& ! category pond drainage rate due to freeboard constraint (m/step)
         dpnd_initialn,&! category runoff rate due to rfrac (m/step)
         dpnd_dlidn     ! category pond loss/gain (+/-) to ice lid freezing/melting (m/step)

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         fsurfn,   & ! category fsurf
         fcondtopn,& ! category fcondtop
         fcondbotn,& ! category fcondbot
         fsensn,   & ! category sensible heat flux
         flatn       ! category latent heat flux

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         snwcnt       ! counter for presence of snow

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
         uatmT   , & ! uatm on T grid (m/s)
         vatmT   , & ! vatm on T grid (m/s)
         wlat    , & ! lateral heat rate (m/s)
         fsw     , & ! incoming shortwave radiation (W/m^2)
         coszen  , & ! cosine solar zenith angle, < 0 for sun below horizon
         rdg_conv, & ! convergence term for ridging (1/s)
         rdg_shear   ! shear term for ridging (1/s)

      real (kind=dbl_kind), dimension(:,:,:,:), allocatable, public :: &
         rsiden    ,&   ! fraction of ice that melts laterally
         salinz    ,&   ! initial salinity  profile (ppt)
         Tmltz          ! initial melting temperature (^oC)

!=======================================================================

      contains

!=======================================================================
!
! Allocate space for all variables
!
      subroutine alloc_flux

      use ice_grid, only : grid_ice

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
         strocnxT_iavg(nx_block,ny_block,max_blocks), & ! ice-ocean stress, x-direction, per ice area
         strocnyT_iavg(nx_block,ny_block,max_blocks), & ! ice-ocean stress, y-direction, per ice area
         sig1       (nx_block,ny_block,max_blocks), & ! normalized principal stress component
         sig2       (nx_block,ny_block,max_blocks), & ! normalized principal stress component
         sigP       (nx_block,ny_block,max_blocks), & ! internal ice pressure (N/m)
         taubxU     (nx_block,ny_block,max_blocks), & ! seabed stress (x) (N/m^2)
         taubyU     (nx_block,ny_block,max_blocks), & ! seabed stress (y) (N/m^2)
         strairxU   (nx_block,ny_block,max_blocks), & ! stress on ice by air, x-direction
         strairyU   (nx_block,ny_block,max_blocks), & ! stress on ice by air, y-direction
         strocnxU   (nx_block,ny_block,max_blocks), & ! ice-ocean stress, x-direction
         strocnyU   (nx_block,ny_block,max_blocks), & ! ice-ocean stress, y-direction
         strtltxU   (nx_block,ny_block,max_blocks), & ! stress due to sea surface slope, x-direction
         strtltyU   (nx_block,ny_block,max_blocks), & ! stress due to sea surface slope, y-direction
         strintxU   (nx_block,ny_block,max_blocks), & ! divergence of internal ice stress, x (N/m^2)
         strintyU   (nx_block,ny_block,max_blocks), & ! divergence of internal ice stress, y (N/m^2)
         daidtd     (nx_block,ny_block,max_blocks), & ! ice area tendency due to transport   (1/s)
         dvidtd     (nx_block,ny_block,max_blocks), & ! ice volume tendency due to transport (m/s)
         dvsdtd     (nx_block,ny_block,max_blocks), & ! snow volume tendency due to transport (m/s)
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
         fmU        (nx_block,ny_block,max_blocks), & ! Coriolis param. * mass in U-cell (kg/s)
         TbU        (nx_block,ny_block,max_blocks), & ! factor for seabed stress (landfast ice)
         zlvl       (nx_block,ny_block,max_blocks), & ! atm level height (momentum) (m)
         zlvs       (nx_block,ny_block,max_blocks), & ! atm level height (scalar quantities) (m)
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
         swuvrdr    (nx_block,ny_block,max_blocks), & ! vis uvr flux, direct (W m-2)
         swuvrdf    (nx_block,ny_block,max_blocks), & ! vis uvr flux, diffuse (W m-2)
         swpardr    (nx_block,ny_block,max_blocks), & ! vis par flux, direct (W m-2)
         swpardf    (nx_block,ny_block,max_blocks), & ! vis par flux, diffuse (W m-2)
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
         daice_da   (nx_block,ny_block,max_blocks), & ! data assimilation concentration increment rate (concentration s-1)
                                                      ! (only used in hadgem drivers)
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
         fsloss     (nx_block,ny_block,max_blocks), & ! rate of snow loss to leads (kg/m^2/s)
         fswthru    (nx_block,ny_block,max_blocks), & ! shortwave penetrating to ocean (W/m^2)
         fswthru_vdr(nx_block,ny_block,max_blocks), & ! vis dir shortwave penetrating to ocean (W/m^2)
         fswthru_vdf(nx_block,ny_block,max_blocks), & ! vis dif shortwave penetrating to ocean (W/m^2)
         fswthru_idr(nx_block,ny_block,max_blocks), & ! nir dir shortwave penetrating to ocean (W/m^2)
         fswthru_idf(nx_block,ny_block,max_blocks), & ! nir dif shortwave penetrating to ocean (W/m^2)
         fswthru_uvrdr (nx_block,ny_block,max_blocks), & ! vis dir uvr SW penetrating to ocean (W/m^2)
         fswthru_uvrdf (nx_block,ny_block,max_blocks), & ! vis dir uvr SW penetrating to ocean (W/m^2)
         fswthru_pardr (nx_block,ny_block,max_blocks), & ! vis dir par SW penetrating to ocean (W/m^2)
         fswthru_pardf (nx_block,ny_block,max_blocks), & ! vis dir par SW penetrating to ocean (W/m^2)
         scale_factor  (nx_block,ny_block,max_blocks), & ! scaling factor for shortwave components
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
         dvsdtt     (nx_block,ny_block,max_blocks), & ! snow volume tendency thermo. (m/s)
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
         uatmT      (nx_block,ny_block,max_blocks), & ! uatm on T grid
         vatmT      (nx_block,ny_block,max_blocks), & ! vatm on T grid
         wlat       (nx_block,ny_block,max_blocks), & ! lateral melt rate (m/s)
         fsw        (nx_block,ny_block,max_blocks), & ! incoming shortwave radiation (W/m^2)
         coszen     (nx_block,ny_block,max_blocks), & ! cosine solar zenith angle, < 0 for sun below horizon
         rdg_conv   (nx_block,ny_block,max_blocks), & ! convergence term for ridging (1/s)
         rdg_shear  (nx_block,ny_block,max_blocks), & ! shear term for ridging (1/s)
         rsiden     (nx_block,ny_block,ncat,max_blocks), & ! fraction of ice that melts laterally
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
         evapn_f    (nx_block,ny_block,ncat,max_blocks), & ! evaporative water flux (kg/m^2/s) by atmosphere model
         dflatndTsfc_f (nx_block,ny_block,ncat,max_blocks), & ! derivative of flatn with respect to Tsfc
         dfsurfndTsfc_f(nx_block,ny_block,ncat,max_blocks), & ! derivative of fsurfn with respect to Tsfc
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
         snwcnt     (nx_block,ny_block,max_blocks,max_nstrm), & ! counter for snow
         salinz     (nx_block,ny_block,nilyr+1,max_blocks), & ! initial salinity  profile (ppt)
         Tmltz      (nx_block,ny_block,nilyr+1,max_blocks), & ! initial melting temperature (^oC)
         stat=ierr)
      if (ierr/=0) call abort_ice('(alloc_flux): Out of memory')

      strax      = c0
      stray      = c0
      uocn       = c0
      vocn       = c0
      ss_tltx    = c0
      ss_tlty    = c0
      hwater     = c0
      strairxT   = c0
      strairyT   = c0
      strocnxT_iavg= c0
      strocnyT_iavg= c0
      sig1       = c0
      sig2       = c0
      sigP       = c0
      taubxU     = c0
      taubyU     = c0
      strairxU   = c0
      strairyU   = c0
      strocnxU   = c0
      strocnyU   = c0
      strtltxU   = c0
      strtltyU   = c0
      strintxU   = c0
      strintyU   = c0
      daidtd     = c0
      dvidtd     = c0
      dvsdtd     = c0
      dagedtd    = c0
      dardg1dt   = c0
      dardg2dt   = c0
      dvirdgdt   = c0
      opening    = c0
      stressp_1  = c0
      stressp_2  = c0
      stressp_3  = c0
      stressp_4  = c0
      stressm_1  = c0
      stressm_2  = c0
      stressm_3  = c0
      stressm_4  = c0
      stress12_1 = c0
      stress12_2 = c0
      stress12_3 = c0
      stress12_4 = c0
      fmU        = c0
      TbU        = c0
      zlvl       = c0
      zlvs       = c0
      uatm       = c0
      vatm       = c0
      wind       = c0
      potT       = c0
      Tair       = c0
      Qa         = c0
      rhoa       = c0
      swvdr      = c0
      swvdf      = c0
      swidr      = c0
      swidf      = c0
      swuvrdr    = c0
      swuvrdf    = c0
      swpardr    = c0
      swpardf    = c0
      flw        = c0
      frain      = c0
      fsnow      = c0
      sss        = c0
      sst        = c0
      frzmlt     = c0
      frzmlt_init= c0
      Tf         = c0
      qdp        = c0
      hmix       = c0
      daice_da   = c0
      fsens      = c0
      flat       = c0
      fswabs     = c0
      fswint_ai  = c0
      flwout     = c0
      Tref       = c0
      Qref       = c0
      Uref       = c0
      evap       = c0
      evaps      = c0
      evapi      = c0
      alvdr      = c0
      alidr      = c0
      alvdf      = c0
      alidf      = c0
      alvdr_ai   = c0
      alidr_ai   = c0
      alvdf_ai   = c0
      alidf_ai   = c0
      albice     = c0
      albsno     = c0
      albpnd     = c0
      apeff_ai   = c0
      snowfrac   = c0
      alvdr_init = c0
      alidr_init = c0
      alvdf_init = c0
      alidf_init = c0
      fpond      = c0
      fresh      = c0
      fsalt      = c0
      fhocn      = c0
      fsloss     = c0
      fswthru    = c0
      fswthru_vdr= c0
      fswthru_vdf= c0
      fswthru_idr= c0
      fswthru_idf= c0
      fswthru_uvrdr = c0
      fswthru_uvrdf = c0
      fswthru_pardr = c0
      fswthru_pardf = c0
      scale_factor  = c0
      strairx_ocn= c0
      strairy_ocn= c0
      fsens_ocn  = c0
      flat_ocn   = c0
      flwout_ocn = c0
      evap_ocn   = c0
      alvdr_ocn  = c0
      alidr_ocn  = c0
      alvdf_ocn  = c0
      alidf_ocn  = c0
      Tref_ocn   = c0
      Qref_ocn   = c0
      fsurf      = c0
      fcondtop   = c0
      fcondbot   = c0
      fbot       = c0
      Tbot       = c0
      Tsnice     = c0
      congel     = c0
      frazil     = c0
      snoice     = c0
      meltt      = c0
      melts      = c0
      meltb      = c0
      meltl      = c0
      dsnow      = c0
      daidtt     = c0
      dvidtt     = c0
      dvsdtt     = c0
      dagedtt    = c0
      mlt_onset  = c0
      frz_onset  = c0
      frazil_diag= c0
      fresh_ai   = c0
      fsalt_ai   = c0
      fhocn_ai   = c0
      fswthru_ai = c0
      fresh_da   = c0
      fsalt_da   = c0
      uatmT      = c0
      vatmT      = c0
      wlat       = c0
      fsw        = c0
      coszen     = c0
      rdg_conv   = c0
      rdg_shear  = c0
      rsiden     = c0
      dardg1ndt  = c0
      dardg2ndt  = c0
      dvirdgndt  = c0
      aparticn   = c0
      krdgn      = c0
      ardgn      = c0
      vrdgn      = c0
      araftn     = c0
      vraftn     = c0
      aredistn   = c0
      vredistn   = c0
      fsurfn_f   = c0
      fcondtopn_f= c0
      fsensn_f   = c0
      flatn_f    = c0
      evapn_f    = c0
      dflatndTsfc_f = c0
      dfsurfndTsfc_f= c0
      meltsn     = c0
      melttn     = c0
      meltbn     = c0
      congeln    = c0
      snoicen    = c0
      keffn_top  = c0
      fsurfn     = c0
      fcondtopn  = c0
      fcondbotn  = c0
      fsensn     = c0
      flatn      = c0
      albcnt     = c0
      snwcnt     = c0
      salinz     = c0
      Tmltz      = c0

      if (grid_ice == "CD" .or. grid_ice == "C") then
         allocate( &
            taubxN     (nx_block,ny_block,max_blocks), & ! seabed stress (x) at N points (N/m^2)
            taubyN     (nx_block,ny_block,max_blocks), & ! seabed stress (y) at N points (N/m^2)
            strairxN   (nx_block,ny_block,max_blocks), & ! stress on ice by air, x-direction at N points
            strairyN   (nx_block,ny_block,max_blocks), & ! stress on ice by air, y-direction at N points
            strocnxN   (nx_block,ny_block,max_blocks), & ! ice-ocean stress, x-direction at N points
            strocnyN   (nx_block,ny_block,max_blocks), & ! ice-ocean stress, y-direction at N points
            strtltxN   (nx_block,ny_block,max_blocks), & ! stress due to sea surface slope, x-direction at N points
            strtltyN   (nx_block,ny_block,max_blocks), & ! stress due to sea surface slope, y-direction at N points
            strintxN   (nx_block,ny_block,max_blocks), & ! divergence of internal ice stress, x at N points (N/m^2)
            strintyN   (nx_block,ny_block,max_blocks), & ! divergence of internal ice stress, y at N points (N/m^2)
            fmN        (nx_block,ny_block,max_blocks), & ! Coriolis param. * mass in N-cell (kg/s)
            TbN        (nx_block,ny_block,max_blocks), & ! factor for seabed stress (landfast ice)
            taubxE     (nx_block,ny_block,max_blocks), & ! seabed stress (x) at E points (N/m^2)
            taubyE     (nx_block,ny_block,max_blocks), & ! seabed stress (y) at E points (N/m^2)
            strairxE   (nx_block,ny_block,max_blocks), & ! stress on ice by air, x-direction at E points
            strairyE   (nx_block,ny_block,max_blocks), & ! stress on ice by air, y-direction at E points
            strocnxE   (nx_block,ny_block,max_blocks), & ! ice-ocean stress, x-direction at E points
            strocnyE   (nx_block,ny_block,max_blocks), & ! ice-ocean stress, y-direction at E points
            strtltxE   (nx_block,ny_block,max_blocks), & ! stress due to sea surface slope, x-direction at E points
            strtltyE   (nx_block,ny_block,max_blocks), & ! stress due to sea surface slope, y-direction at E points
            strintxE   (nx_block,ny_block,max_blocks), & ! divergence of internal ice stress, x at E points (N/m^2)
            strintyE   (nx_block,ny_block,max_blocks), & ! divergence of internal ice stress, y at E points (N/m^2)
            fmE        (nx_block,ny_block,max_blocks), & ! Coriolis param. * mass in E-cell (kg/s)
            TbE        (nx_block,ny_block,max_blocks), & ! factor for seabed stress (landfast ice)
            stresspT   (nx_block,ny_block,max_blocks), & ! sigma11+sigma22
            stressmT   (nx_block,ny_block,max_blocks), & ! sigma11-sigma22
            stress12T  (nx_block,ny_block,max_blocks), & ! sigma12
            stresspU   (nx_block,ny_block,max_blocks), & ! sigma11+sigma22
            stressmU   (nx_block,ny_block,max_blocks), & ! sigma11-sigma22
            stress12U  (nx_block,ny_block,max_blocks), & ! sigma12
            stat=ierr)
         if (ierr/=0) call abort_ice('(alloc_flux): Out of memory (C or CD grid)')

         taubxN     = c0
         taubyN     = c0
         strairxN   = c0
         strairyN   = c0
         strocnxN   = c0
         strocnyN   = c0
         strtltxN   = c0
         strtltyN   = c0
         strintxN   = c0
         strintyN   = c0
         fmN        = c0
         TbN        = c0
         taubxE     = c0
         taubyE     = c0
         strairxE   = c0
         strairyE   = c0
         strocnxE   = c0
         strocnyE   = c0
         strtltxE   = c0
         strtltyE   = c0
         strintxE   = c0
         strintyE   = c0
         fmE        = c0
         TbE        = c0
         stresspT   = c0
         stressmT   = c0
         stress12T  = c0
         stresspU   = c0
         stressmU   = c0
         stress12U  = c0
      endif

      ! Pond diagnostics
      allocate( &
         dpnd_flush   (nx_block,ny_block,max_blocks), & ! pond flushing rate due to ice permeability (m/step)
         dpnd_expon   (nx_block,ny_block,max_blocks), & ! exponential pond drainage rate (m/step)
         dpnd_freebd  (nx_block,ny_block,max_blocks), & ! pond drainage rate due to freeboard constraint (m/step)
         dpnd_initial (nx_block,ny_block,max_blocks), & ! runoff rate due to rfrac (m/step)
         dpnd_dlid    (nx_block,ny_block,max_blocks), & ! pond loss/gain (+/-) to ice lid freezing/melting (m/step)
         dpnd_melt    (nx_block,ny_block,max_blocks), & ! pond 'drainage' due to ice melting (m / step)
         dpnd_ridge   (nx_block,ny_block,max_blocks), & ! pond 'drainage' due to ridging (m)
         dpnd_flushn  (nx_block,ny_block,ncat,max_blocks), & ! category pond flushing rate due to ice permeability (m/step)
         dpnd_exponn  (nx_block,ny_block,ncat,max_blocks), & ! category exponential pond drainage rate (m/step)
         dpnd_freebdn (nx_block,ny_block,ncat,max_blocks), & ! category pond drainage rate due to freeboard constraint (m/step)
         dpnd_initialn(nx_block,ny_block,ncat,max_blocks), & ! category runoff rate due to rfrac (m/step)
         dpnd_dlidn   (nx_block,ny_block,ncat,max_blocks), & ! category pond loss/gain (+/-) to ice lid freezing/melting (m/step)
         stat=ierr)
      if (ierr/=0) call abort_ice('(alloc_flux): Out of memory (ponds)')

      dpnd_flush   = c0
      dpnd_expon   = c0
      dpnd_freebd  = c0
      dpnd_initial = c0
      dpnd_dlid    = c0
      dpnd_melt    = c0
      dpnd_ridge   = c0
      dpnd_flushn  = c0
      dpnd_exponn  = c0
      dpnd_freebdn = c0
      dpnd_initialn= c0
      dpnd_dlidn   = c0

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
      zlvl  (:,:,:) = c10             ! atm level height (momentum) (m)
      zlvs  (:,:,:) = c10             ! atm level height (scalar quantities) (m)
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

      swuvrdr(:,:,:) = c0                ! visible uvr flux, direct (W/m^2)
      swuvrdf(:,:,:) = c0                ! visible uvr flux, diffuse (W/m^2)
      swpardr(:,:,:) = c0                ! visible par flux, direct (W/m^2)
      swpardf(:,:,:) = c0                ! visible par flux, diffuse (W/m^2)

      fiso_atm  (:,:,:,:) = c0           ! isotope deposition rate (kg/m2/s)
      faero_atm (:,:,:,:) = c0           ! aerosol deposition rate (kg/m2/s)
      flux_bio_atm (:,:,:,:) = c0        ! zaero and bio deposition rate (kg/m2/s)

      !-----------------------------------------------------------------
      ! fluxes received from ocean
      !-----------------------------------------------------------------

      ss_tltx (:,:,:) = c0              ! sea surface tilt (m/m)
      ss_tlty (:,:,:) = c0
      uocn    (:,:,:) = c0              ! surface ocean currents (m/s)
      vocn    (:,:,:) = c0
      frzmlt  (:,:,:) = c0              ! freezing/melting potential (W/m^2)
      frzmlt_init(:,:,:) = c0           ! freezing/melting potential (W/m^2)
      sss     (:,:,:) = 34.0_dbl_kind   ! sea surface salinity (ppt)

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

      strocnxT_iavg (:,:,:) = c0 ! ice-ocean stress, x-direction (T-cell)
      strocnyT_iavg (:,:,:) = c0 ! ice-ocean stress, y-direction (T-cell)
      fresh   (:,:,:) = c0
      fsalt   (:,:,:) = c0
      fpond   (:,:,:) = c0
      fhocn   (:,:,:) = c0
      fswthru (:,:,:) = c0
      fswthru_vdr (:,:,:) = c0
      fswthru_vdf (:,:,:) = c0
      fswthru_idr (:,:,:) = c0
      fswthru_idf (:,:,:) = c0
      fswthru_uvrdr (:,:,:) = c0
      fswthru_uvrdf (:,:,:) = c0
      fswthru_pardr (:,:,:) = c0
      fswthru_pardf (:,:,:) = c0
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

      allocate(fswthrun_ai(nx_block,ny_block,ncat,max_blocks))
      fswthrun_ai(:,:,:,:) = c0

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

      fsurf   (:,:,:) = c0
      fcondtop(:,:,:) = c0
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
      fswthru_uvrdr(:,:,:)   = c0
      fswthru_uvrdf(:,:,:)   = c0
      fswthru_pardr(:,:,:)   = c0
      fswthru_pardf(:,:,:)   = c0

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

      use ice_state, only: aice, vice, vsno, trcr
      use ice_arrays_column, only: &
          hfreebd, hdraft, hridge, distrdg, hkeel, dkeel, lfloe, dfloe, &
          Cdn_atm_skin, Cdn_atm_floe, Cdn_atm_pond, Cdn_atm_rdg, &
          Cdn_ocn_skin, Cdn_ocn_floe, Cdn_ocn_keel, Cdn_atm_ratio, &
          Cdn_atm, Cdn_ocn

      logical (kind=log_kind) :: &
          formdrag, &
          tr_pond,  &
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
      call icepack_query_tracer_flags(tr_pond_out=tr_pond)
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
      dvsdtt (:,:,:) = vsno(:,:,:) ! temporary initial volume
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

      ! Extra pond diagnostics
      dpnd_flush(:,:,:)   = c0
      dpnd_expon(:,:,:)   = c0
      dpnd_freebd(:,:,:)  = c0
      dpnd_initial(:,:,:) = c0
      dpnd_dlid(:,:,:)    = c0
      dpnd_melt(:,:,:)    = c0
      dpnd_ridge(:,:,:)   = c0
      dpnd_flushn(:,:,:,:)   = c0
      dpnd_exponn(:,:,:,:)   = c0
      dpnd_freebdn(:,:,:,:)  = c0
      dpnd_initialn(:,:,:,:) = c0
      dpnd_dlidn(:,:,:,:)    = c0

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

      use ice_state, only: aice, vice, vsno, trcr, strength, divu, shear, vort
      use ice_grid,  only: grid_ice

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
      divu    (:,:,:) = c0
      shear   (:,:,:) = c0
      vort    (:,:,:) = c0
      taubxU  (:,:,:) = c0
      taubyU  (:,:,:) = c0
      strength (:,:,:) = c0
      strocnxU(:,:,:) = c0
      strocnyU(:,:,:) = c0
      strairxU(:,:,:) = c0
      strairyU(:,:,:) = c0
      strtltxU(:,:,:) = c0
      strtltyU(:,:,:) = c0
      strintxU(:,:,:) = c0
      strintyU(:,:,:) = c0
      dardg1dt(:,:,:) = c0
      dardg2dt(:,:,:) = c0
      dvirdgdt(:,:,:) = c0
      opening (:,:,:) = c0
      daidtd  (:,:,:) = aice(:,:,:) ! temporary initial area
      dvidtd  (:,:,:) = vice(:,:,:) ! temporary initial volume
      dvsdtd  (:,:,:) = vsno(:,:,:) ! temporary initial volume
      if (tr_iage) &
         dagedtd (:,:,:) = trcr(:,:,nt_iage,:) ! temporary initial age
      fmU     (:,:,:) = c0
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

      if (grid_ice == "CD" .or. grid_ice == "C") then
         taubxE     (:,:,:) = c0
         taubyE     (:,:,:) = c0
         strocnxE   (:,:,:) = c0
         strocnyE   (:,:,:) = c0
         strairxE   (:,:,:) = c0
         strairyE   (:,:,:) = c0
         strtltxE   (:,:,:) = c0
         strtltyE   (:,:,:) = c0
         strintxE   (:,:,:) = c0
         strintyE   (:,:,:) = c0
         fmE        (:,:,:) = c0
         TbE        (:,:,:) = c0
         taubxN     (:,:,:) = c0
         taubyN     (:,:,:) = c0
         strocnxN   (:,:,:) = c0
         strocnyN   (:,:,:) = c0
         strairxN   (:,:,:) = c0
         strairyN   (:,:,:) = c0
         strtltxN   (:,:,:) = c0
         strtltyN   (:,:,:) = c0
         strintxN   (:,:,:) = c0
         strintyN   (:,:,:) = c0
         fmN        (:,:,:) = c0
         TbN        (:,:,:) = c0
      end if
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
          Tffresh, puny

      integer (kind=int_kind) :: &
          i, j    ! horizontal indices

      character(len=*), parameter :: subname = '(scale_fluxes)'

      call icepack_query_parameters(stefan_boltzmann_out=stefan_boltzmann, &
         Tffresh_out=Tffresh, puny_out=puny)
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
            ! Special case where aice_init was zero and aice > 0.
            if (flwout(i,j) > -puny) &
               flwout  (i,j) = -stefan_boltzmann *(Tf(i,j) + Tffresh)**4
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
