!=======================================================================

! Grid-dependent arrays needed for column package
! These were originally module variables in modules that became part of
! the column package

! author: Elizabeth C. Hunke, LANL

      module ice_arrays_column

      use ice_kinds_mod
      use ice_constants, only : c0
      use ice_fileunits, only: nu_diag
      use ice_blocks, only: nx_block, ny_block
      use ice_domain_size, only: max_blocks, ncat, nilyr, nslyr, &
          nblyr, nfsd, nfreq
      use icepack_intfc, only: icepack_nspint_3bd
      use icepack_intfc, only: icepack_query_tracer_sizes, icepack_query_parameters, &
          icepack_query_tracer_flags, &
          icepack_warnings_flush, icepack_warnings_aborted, icepack_query_tracer_sizes

      implicit none
      private

      public :: alloc_arrays_column

      ! icepack_atmo.F90
      ! Cdn variables on the T-grid
      real (kind=dbl_kind), public, dimension (:,:,:), allocatable :: &
         Cdn_atm     , & ! atm drag coefficient
         Cdn_ocn     , & ! ocn drag coefficient
                         ! form drag
         hfreebd,      & ! freeboard (m)
         hdraft,       & ! draft of ice + snow column (Stoessel1993)
         hridge,       & ! ridge height
         distrdg,      & ! distance between ridges
         hkeel,        & ! keel depth
         dkeel,        & ! distance between keels
         lfloe,        & ! floe length
         dfloe,        & ! distance between floes
         Cdn_atm_skin, & ! neutral skin drag coefficient
         Cdn_atm_floe, & ! neutral floe edge drag coefficient
         Cdn_atm_pond, & ! neutral pond edge drag coefficient
         Cdn_atm_rdg,  & ! neutral ridge drag coefficient
         Cdn_ocn_skin, & ! skin drag coefficient
         Cdn_ocn_floe, & ! floe edge drag coefficient
         Cdn_ocn_keel, & ! keel drag coefficient
         Cdn_atm_ratio   ! ratio drag atm / neutral drag atm

!-------------------------------------------------------------------
! a note regarding hi_min and hin_max(0):
! both represent a minimum ice thickness.  hin_max(0) is
! intended to be used for particular numerical implementations
! of category conversions in the ice thickness distribution.
! hi_min is a more general purpose parameter, but is specifically
! for maintaining stability in the thermodynamics.
! hin_max(0) = 0.1 m for the delta function itd
! hin_max(0) = 0.0 m for linear remapping
!
! Also note that the upper limit on the thickest category
! is only used for the linear remapping scheme
! and it is not a true upper limit on the thickness
!-------------------------------------------------------------------

      ! icepack_itd.F90
      real (kind=dbl_kind), public, allocatable :: &
         hin_max(:)   ! category limits (m)

      character (len=35), public, allocatable :: &
         c_hi_range(:)! string for history output

      ! icepack_snow.F90
      real (kind=dbl_kind), public, dimension (:,:,:), allocatable :: &
         meltsliq     ! snow melt mass (kg/m^2/step-->kg/m^2/day)

      real (kind=dbl_kind), public, dimension (:,:,:,:), allocatable :: &
         meltsliqn    ! snow melt mass in category n (kg/m^2)

      ! icepack_meltpond_lvl.F90
      real (kind=dbl_kind), public, dimension (:,:,:,:), allocatable :: &
         dhsn, &      ! depth difference for snow on sea ice and pond ice
         ffracn       ! fraction of fsurfn used to melt ipond

      ! icepack_shortwave.F90
      ! category albedos
      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         alvdrn, &    ! visible direct albedo           (fraction)
         alidrn, &    ! near-ir direct albedo           (fraction)
         alvdfn, &    ! visible diffuse albedo          (fraction)
         alidfn       ! near-ir diffuse albedo          (fraction)

      ! albedo components for history
      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         albicen, &   ! bare ice
         albsnon, &   ! snow
         albpndn, &   ! pond
         apeffn       ! effective pond area used for radiation calculation

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         snowfracn    ! Category snow fraction used in radiation

      ! shortwave components
      real (kind=dbl_kind), dimension (:,:,:,:,:), allocatable, public :: &
         Iswabsn      ! SW radiation absorbed in ice layers (W m-2)

      real (kind=dbl_kind), dimension (:,:,:,:,:), allocatable, public :: &
         Sswabsn      ! SW radiation absorbed in snow layers (W m-2)

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         fswsfcn      , & ! SW absorbed at ice/snow surface (W m-2)
         fswthrun     , & ! SW through ice to ocean            (W/m^2)
         fswthrun_vdr , & ! vis dir SW through ice to ocean            (W/m^2)
         fswthrun_vdf , & ! vis dif SW through ice to ocean            (W/m^2)
         fswthrun_idr , & ! nir dir SW through ice to ocean            (W/m^2)
         fswthrun_idf , & ! nir dif SW through ice to ocean            (W/m^2)
         fswthrun_uvrdr,& ! vis uvr dir SW through ice to ocean        (W/m^2)
         fswthrun_uvrdf,& ! vis uvr dif SW through ice to ocean        (W/m^2)
         fswthrun_pardr,& ! vis par dir SW through ice to ocean        (W/m^2)
         fswthrun_pardf,& ! vis par dif SW through ice to ocean        (W/m^2)
         fswintn          ! SW absorbed in ice interior, below surface (W m-2)

      real (kind=dbl_kind), dimension (:,:,:,:,:), allocatable, public :: &
         fswpenln     ! visible SW entering ice layers (W m-2)

      ! biogeochemistry components

      real (kind=dbl_kind), dimension (:), allocatable, public :: &
         bgrid            , &  ! biology nondimensional vertical grid points
         igrid            , &  ! biology vertical interface points
         cgrid            , &  ! CICE vertical coordinate
         icgrid           , &  ! interface grid for CICE (shortwave variable)
         swgrid                ! grid for ice tracers used in dEdd scheme

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         first_ice_real     ! .true. = c1, .false. = c0

      logical (kind=log_kind), dimension (:,:,:,:), allocatable, public :: &
         first_ice      ! distinguishes ice that disappears (e.g. melts)
                        ! and reappears (e.g. transport) in a grid cell
                        ! during a single time step from ice that was
                        ! there the entire time step (true until ice forms)

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         ocean_bio      ! contains all the ocean bgc tracer concentrations

      ! diagnostic fluxes
      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         fbio_snoice, & ! fluxes from snow to ice
         fbio_atmice    ! fluxes from atm to ice

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         ocean_bio_all      ! fixed order, all values even for tracers false
                            ! N(1:icepack_max_algae) = 1:icepack_max_algae
                            ! Nit = icepack_max_algae + 1
                            ! DOC(1:icepack_max_doc) = icepack_max_algae + 2: icepack_max_algae + icepack_max_doc + 1
                            ! DIC(1:icepack_max_dic) = icepack_max_algae + icepack_max_doc + 2: icepack_max_algae + icepack_max_doc + 1 + icepack_max_dic
                            ! chl(1:icepack_max_algae) =  icepack_max_algae + icepack_max_doc + 2 + icepack_max_dic: &
                            !                     2*icepack_max_algae + icepack_max_doc + 1 + icepack_max_dic
                            ! Am =  2*icepack_max_algae + icepack_max_doc + 2 + icepack_max_dic
                            ! Sil=  2*icepack_max_algae + icepack_max_doc + 3 + icepack_max_dic
                            ! DMSPp=  2*icepack_max_algae + icepack_max_doc + 4 + icepack_max_dic
                            ! DMSPd=  2*icepack_max_algae + icepack_max_doc + 5 + icepack_max_dic
                            ! DMS  =  2*icepack_max_algae + icepack_max_doc + 6 + icepack_max_dic
                            ! PON  =  2*icepack_max_algae + icepack_max_doc + 7 + icepack_max_dic
                            ! DON(1:icepack_max_don)  =  2*icepack_max_algae + icepack_max_doc + 8 + icepack_max_dic:
                            !                    2*icepack_max_algae + icepack_max_doc + 7 + icepack_max_dic + icepack_max_don
                            ! Fed(1:icepack_max_fe) = 2*icepack_max_algae + icepack_max_doc + 8 + icepack_max_dic + icepack_max_don:
                            !                2*icepack_max_algae + icepack_max_doc + 7 + icepack_max_dic + icepack_max_don + icepack_max_fe
                            ! Fep(1:icepack_max_fe) = 2*icepack_max_algae + icepack_max_doc + 8 + icepack_max_dic + icepack_max_don + icepack_max_fe:
                            !                2*icepack_max_algae + icepack_max_doc + 7 + icepack_max_dic + icepack_max_don + 2*icepack_max_fe
                            ! zaero(1:icepack_max_aero) = 2*icepack_max_algae + icepack_max_doc + 8 + icepack_max_dic + icepack_max_don + 2*icepack_max_fe:
                            !                     2*icepack_max_algae + icepack_max_doc + 7 + icepack_max_dic + icepack_max_don + 2*icepack_max_fe
                            !                     + icepack_max_aero
                            ! humic ==  2*icepack_max_algae + icepack_max_doc + 8 + icepack_max_dic + icepack_max_don + 2*icepack_max_fe
                            !                     + icepack_max_aero

      integer (kind=int_kind), dimension(:,:,:,:), allocatable, public :: &
        algal_peak          ! vertical location of algal maximum, 0 if no maximum

      real (kind=dbl_kind), dimension (:,:,:,:,:), allocatable, public :: &
         Zoo        ! N losses accumulated in timestep (ie. zooplankton/bacteria)
                    ! mmol/m^3

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         dhbr_top     , & ! brine top change
         dhbr_bot         ! brine bottom change

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         grow_net       , & ! Specific growth rate (/s) per grid cell
         PP_net         , & ! Total production (mg C/m^2/s) per grid cell
         hbri               ! brine height, area-averaged for comparison with hi (m)

      real (kind=dbl_kind), dimension (:,:,:,:,:), allocatable, public :: &
         bphi           , & ! porosity of layers
         bTiz               ! layer temperatures interpolated on bio grid (C)

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         darcy_V            ! darcy velocity positive up (m/s)

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         chl_net     , & ! Total chla (mg chla/m^2) per grid cell
         NO_net          ! Total nitrate per grid cell

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         sice_rho     ! avg sea ice density  (kg/m^3)  ! ech: diagnostic only?

      real (kind=dbl_kind), dimension (:,:,:,:,:), allocatable, public :: &
         zfswin       ! Shortwave flux into layers interpolated on bio grid  (W/m^2)

      real (kind=dbl_kind), dimension (:,:,:,:,:), allocatable, public :: &
         iDi      , & ! igrid Diffusivity (m^2/s)
         iki          ! Ice permeability (m^2)

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         upNO     , & ! nitrate uptake rate (mmol/m^2/d) times aice
         upNH         ! ammonium uptake rate (mmol/m^2/d) times aice

      real (kind=dbl_kind), dimension(:,:,:,:,:), allocatable, public :: &
         trcrn_sw        ! bgc tracers active in the delta-Eddington shortwave
                         ! calculation on the shortwave grid (swgrid)

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         ice_bio_net  , &   ! depth integrated tracer (mmol/m^2)
         snow_bio_net       ! depth integrated snow tracer (mmol/m^2)

      logical (kind=log_kind), public :: &
         oceanmixed_ice, &  ! if true, use internal ocean mixed layer
         restore_bgc        !

      character(char_len), public :: &
         fe_data_type   ! 'default', 'clim'

      character(char_len_long), public :: &
         bgc_data_dir   ! directory for biogeochemistry data

      real (kind=dbl_kind), dimension(:), allocatable, public :: &
         R_chl2N,       &  ! 3 algal chlorophyll to N (mg/mmol)
         R_C2N             ! algal C to N (mole/mole)

      ! floe size distribution
      real(kind=dbl_kind), dimension(:), allocatable, public ::  &
         floe_rad_l,    &  ! fsd size lower bound in m (radius)
         floe_rad_c,    &  ! fsd size bin centre in m (radius)
         floe_binwidth     ! fsd size bin width in m (radius)

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
        wave_sig_ht        ! significant height of waves (m)

      real (kind=dbl_kind), dimension (:), allocatable, public :: &
         wavefreq,      &  ! wave frequencies
         dwavefreq         ! wave frequency bin widths

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         wave_spectrum, &  ! wave spectrum
         ! change in floe size distribution due to processes
         d_afsd_newi, d_afsd_latg, d_afsd_latm, d_afsd_wave, d_afsd_weld

      character (len=35), public, allocatable :: c_fsd_range(:)

!=======================================================================

      contains

!=======================================================================

      subroutine alloc_arrays_column
        ! Allocate column arrays
        use ice_exit, only: abort_ice
        integer (int_kind) :: max_nbtrcr, max_algae, max_aero, &
           nmodal1, nmodal2, max_don
        integer (int_kind) :: ierr, ntrcr

      character(len=*),parameter :: subname='(alloc_arrays_column)'

      call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
      call icepack_query_tracer_sizes( max_nbtrcr_out=max_nbtrcr, &
         max_algae_out=max_algae, max_aero_out=max_aero, &
         nmodal1_out=nmodal1, nmodal2_out=nmodal2, max_don_out=max_don)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__,line= __LINE__)

      allocate(                                       &
         Cdn_atm      (nx_block,ny_block,max_blocks), & ! atm drag coefficient
         Cdn_ocn      (nx_block,ny_block,max_blocks), & ! ocn drag coefficient
         hfreebd      (nx_block,ny_block,max_blocks), & ! freeboard (m)
         hdraft       (nx_block,ny_block,max_blocks), & ! draft of ice + snow column (Stoessel1993)
         hridge       (nx_block,ny_block,max_blocks), & ! ridge height
         distrdg      (nx_block,ny_block,max_blocks), & ! distance between ridges
         hkeel        (nx_block,ny_block,max_blocks), & ! keel depth
         dkeel        (nx_block,ny_block,max_blocks), & ! distance between keels
         lfloe        (nx_block,ny_block,max_blocks), & ! floe length
         dfloe        (nx_block,ny_block,max_blocks), & ! distance between floes
         Cdn_atm_skin (nx_block,ny_block,max_blocks), & ! neutral skin drag coefficient
         Cdn_atm_floe (nx_block,ny_block,max_blocks), & ! neutral floe edge drag coefficient
         Cdn_atm_pond (nx_block,ny_block,max_blocks), & ! neutral pond edge drag coefficient
         Cdn_atm_rdg  (nx_block,ny_block,max_blocks), & ! neutral ridge drag coefficient
         Cdn_ocn_skin (nx_block,ny_block,max_blocks), & ! skin drag coefficient
         Cdn_ocn_floe (nx_block,ny_block,max_blocks), & ! floe edge drag coefficient
         Cdn_ocn_keel (nx_block,ny_block,max_blocks), & ! keel drag coefficient
         Cdn_atm_ratio(nx_block,ny_block,max_blocks), & ! ratio drag atm / neutral drag atm
         grow_net     (nx_block,ny_block,max_blocks), & ! Specific growth rate (/s) per grid cell
         PP_net       (nx_block,ny_block,max_blocks), & ! Total production (mg C/m^2/s) per grid cell
         hbri         (nx_block,ny_block,max_blocks), & ! brine height, area-averaged for comparison with hi (m)
         chl_net      (nx_block,ny_block,max_blocks), & ! Total chla (mg chla/m^2) per grid cell
         NO_net       (nx_block,ny_block,max_blocks), & ! Total nitrate per grid cell
         upNO         (nx_block,ny_block,max_blocks), & ! nitrate uptake rate (mmol/m^2/d) times aice
         upNH         (nx_block,ny_block,max_blocks), & ! ammonium uptake rate (mmol/m^2/d) times aice
         meltsliq     (nx_block,ny_block,max_blocks), & ! snow melt mass (kg/m^2)
         meltsliqn    (nx_block,ny_block,ncat,max_blocks), & ! snow melt mass in category n (kg/m^2)
         dhsn         (nx_block,ny_block,ncat,max_blocks), & ! depth difference for snow on sea ice and pond ice
         ffracn       (nx_block,ny_block,ncat,max_blocks), & ! fraction of fsurfn used to melt ipond
         alvdrn       (nx_block,ny_block,ncat,max_blocks), & ! visible direct albedo           (fraction)
         alidrn       (nx_block,ny_block,ncat,max_blocks), & ! near-ir direct albedo           (fraction)
         alvdfn       (nx_block,ny_block,ncat,max_blocks), & ! visible diffuse albedo          (fraction)
         alidfn       (nx_block,ny_block,ncat,max_blocks), & ! near-ir diffuse albedo          (fraction)
         albicen      (nx_block,ny_block,ncat,max_blocks), & ! bare ice
         albsnon      (nx_block,ny_block,ncat,max_blocks), & ! snow
         albpndn      (nx_block,ny_block,ncat,max_blocks), & ! pond
         apeffn       (nx_block,ny_block,ncat,max_blocks), & ! effective pond area used for radiation calculation
         snowfracn    (nx_block,ny_block,ncat,max_blocks), & ! Category snow fraction used in radiation
         fswsfcn      (nx_block,ny_block,ncat,max_blocks), & ! SW absorbed at ice/snow surface (W m-2)
         fswthrun     (nx_block,ny_block,ncat,max_blocks), & ! SW through ice to ocean            (W/m^2)
         fswthrun_vdr  (nx_block,ny_block,ncat,max_blocks), & ! vis dir SW through ice to ocean            (W/m^2)
         fswthrun_vdf  (nx_block,ny_block,ncat,max_blocks), & ! vis dif SW through ice to ocean            (W/m^2)
         fswthrun_idr  (nx_block,ny_block,ncat,max_blocks), & ! nir dir SW through ice to ocean            (W/m^2)
         fswthrun_idf  (nx_block,ny_block,ncat,max_blocks), & ! nir dif SW through ice to ocean            (W/m^2)
         fswthrun_uvrdr(nx_block,ny_block,ncat,max_blocks), & ! vis uvr dir SW uhrough ice to ocean        (W/m^2)
         fswthrun_uvrdf(nx_block,ny_block,ncat,max_blocks), & ! vis uvr dif SW through ice to ocean        (W/m^2)
         fswthrun_pardr(nx_block,ny_block,ncat,max_blocks), & ! vis par dir SW through ice to ocean        (W/m^2)
         fswthrun_pardf(nx_block,ny_block,ncat,max_blocks), & ! vis par dif SW through ice to ocean        (W/m^2)
         fswintn      (nx_block,ny_block,ncat,max_blocks), & ! SW absorbed in ice interior, below surface (W m-2)
         first_ice_real                                    &
                      (nx_block,ny_block,ncat,max_blocks), & ! .true. = c1, .false. = c0
         first_ice    (nx_block,ny_block,ncat,max_blocks), & ! distinguishes ice that disappears (melts) and reappears (transport)
         dhbr_top     (nx_block,ny_block,ncat,max_blocks), & ! brine top change
         dhbr_bot     (nx_block,ny_block,ncat,max_blocks), & ! brine bottom change
         darcy_V      (nx_block,ny_block,ncat,max_blocks), & ! darcy velocity positive up (m/s)
         sice_rho     (nx_block,ny_block,ncat,max_blocks), & ! avg sea ice density  (kg/m^3)  ! ech: diagnostic only?
         Iswabsn      (nx_block,ny_block,nilyr,ncat,max_blocks), & ! SW radiation absorbed in ice layers (W m-2)
         Sswabsn      (nx_block,ny_block,nslyr,ncat,max_blocks), & ! SW radiation absorbed in snow layers (W m-2)
         fswpenln     (nx_block,ny_block,nilyr+1,ncat,max_blocks), & ! visible SW entering ice layers (W m-2)
         Zoo          (nx_block,ny_block,nblyr+1,ncat,max_blocks), & ! N losses accumulated in timestep (ie. zooplankton/bacteria)
         zfswin       (nx_block,ny_block,nblyr+1,ncat,max_blocks), & ! Shortwave flux into layers interpolated on bio grid  (W/m^2)
         iDi          (nx_block,ny_block,nblyr+1,ncat,max_blocks), & ! igrid Diffusivity (m^2/s)
         iki          (nx_block,ny_block,nblyr+1,ncat,max_blocks), & ! Ice permeability (m^2)
         bphi         (nx_block,ny_block,nblyr+2,ncat,max_blocks), & ! porosity of layers
         bTiz         (nx_block,ny_block,nblyr+2,ncat,max_blocks), &    ! layer temperatures interpolated on bio grid (C)
         stat=ierr)
      if (ierr/=0) call abort_ice(subname//': Out of Memory1')

      Cdn_atm       = c0
      Cdn_ocn       = c0
      hfreebd       = c0
      hdraft        = c0
      hridge        = c0
      distrdg       = c0
      hkeel         = c0
      dkeel         = c0
      lfloe         = c0
      dfloe         = c0
      Cdn_atm_skin  = c0
      Cdn_atm_floe  = c0
      Cdn_atm_pond  = c0
      Cdn_atm_rdg   = c0
      Cdn_ocn_skin  = c0
      Cdn_ocn_floe  = c0
      Cdn_ocn_keel  = c0
      Cdn_atm_ratio = c0
      grow_net      = c0
      PP_net        = c0
      hbri          = c0
      chl_net       = c0
      NO_net        = c0
      upNO          = c0
      upNH          = c0
      meltsliq      = c0
      meltsliqn     = c0
      dhsn          = c0
      ffracn        = c0
      alvdrn        = c0
      alidrn        = c0
      alvdfn        = c0
      alidfn        = c0
      albicen       = c0
      albsnon       = c0
      albpndn       = c0
      apeffn        = c0
      snowfracn     = c0
      fswsfcn       = c0
      fswthrun      = c0
      fswthrun_vdr  = c0
      fswthrun_vdf  = c0
      fswthrun_idr  = c0
      fswthrun_idf  = c0
      fswthrun_uvrdr= c0
      fswthrun_uvrdf= c0
      fswthrun_pardr= c0
      fswthrun_pardf= c0
      fswintn       = c0
      first_ice_real= c0
      first_ice     = .false.
      dhbr_top      = c0
      dhbr_bot      = c0
      darcy_V       = c0
      sice_rho      = c0
      Iswabsn       = c0
      Sswabsn       = c0
      fswpenln      = c0
      Zoo           = c0
      zfswin        = c0
      iDi           = c0
      iki           = c0
      bphi          = c0
      bTiz          = c0

      allocate(                                       &
         ocean_bio    (nx_block,ny_block,max_nbtrcr,max_blocks), & ! contains all the ocean bgc tracer concentrations
         fbio_snoice  (nx_block,ny_block,max_nbtrcr,max_blocks), & ! fluxes from snow to ice
         fbio_atmice  (nx_block,ny_block,max_nbtrcr,max_blocks), & ! fluxes from atm to ice
         ocean_bio_all(nx_block,ny_block,max_nbtrcr,max_blocks), & ! fixed order, all values even for tracers false
         ice_bio_net  (nx_block,ny_block,max_nbtrcr,max_blocks), & ! depth integrated tracer (mmol/m^2)
         snow_bio_net (nx_block,ny_block,max_nbtrcr,max_blocks), & ! depth integrated snow tracer (mmol/m^2)
         algal_peak   (nx_block,ny_block,max_algae ,max_blocks), & ! vertical location of algal maximum, 0 if no maximum
         stat=ierr)
      if (ierr/=0) call abort_ice(subname//': Out of Memory2')

      ocean_bio    = c0
      fbio_snoice  = c0
      fbio_atmice  = c0
      ocean_bio_all= c0
      ice_bio_net  = c0
      snow_bio_net = c0
      algal_peak   = 0

      allocate(                                       &
         hin_max(0:ncat)            , & ! category limits (m)
         c_hi_range(ncat)           , & !
         bgrid(nblyr+2)             , & ! biology nondimensional vertical grid points
         igrid(nblyr+1)             , &  ! biology vertical interface points
         cgrid(nilyr+1)             , &  ! CICE vertical coordinate
         icgrid(nilyr+1)            , &  ! interface grid for CICE (shortwave variable)
         swgrid(nilyr+1)            , &  ! grid for ice tracers used in dEdd scheme
         stat=ierr)
      if (ierr/=0) call abort_ice(subname//' Out of Memory3')

      hin_max = c0
      c_hi_range = ''
      bgrid   = c0
      igrid   = c0
      cgrid   = c0
      icgrid  = c0
      swgrid  = c0

      ! floe size distribution
      allocate(                                                   &
         floe_rad_l     (nfsd)      , & ! fsd size lower bound in m (radius)
         floe_rad_c     (nfsd)      , & ! fsd size bin centre in m (radius)
         floe_binwidth  (nfsd)      , & ! fsd size bin width in m (radius)
         c_fsd_range    (nfsd)      , & ! fsd floe_rad bounds (m)
         wavefreq       (nfreq)     , & ! wave frequency
         dwavefreq      (nfreq)     , & ! wave frequency bin widths
         wave_sig_ht    (nx_block,ny_block,          max_blocks), & !
         wave_spectrum  (nx_block,ny_block,nfreq,    max_blocks), & !
         d_afsd_newi    (nx_block,ny_block,nfsd,     max_blocks), & !
         d_afsd_latg    (nx_block,ny_block,nfsd,     max_blocks), & !
         d_afsd_latm    (nx_block,ny_block,nfsd,     max_blocks), & !
         d_afsd_wave    (nx_block,ny_block,nfsd,     max_blocks), & !
         d_afsd_weld    (nx_block,ny_block,nfsd,     max_blocks), & !
         stat=ierr)
      if (ierr/=0) call abort_ice(subname//' Out of Memory5')

      floe_rad_l     = c0
      floe_rad_c     = c0
      floe_binwidth  = c0
      c_fsd_range    = ''
      wavefreq       = c0
      dwavefreq      = c0
      wave_sig_ht    = c0
      wave_spectrum  = c0
      d_afsd_newi    = c0
      d_afsd_latg    = c0
      d_afsd_latm    = c0
      d_afsd_wave    = c0
      d_afsd_weld    = c0

      end subroutine alloc_arrays_column

!=======================================================================

      end module ice_arrays_column

!=======================================================================
