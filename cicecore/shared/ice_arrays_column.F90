!=======================================================================

! Grid-dependent arrays needed for column package
! These were originally module variables in modules that became part of 
! the column package

! author: Elizabeth C. Hunke, LANL

      module ice_arrays_column

      use ice_kinds_mod
      use ice_blocks, only: nx_block, ny_block
      use icepack_intfc, only: icepack_nspint
      use ice_domain_size, only: max_blocks, ncat, nilyr, nslyr, &
           nblyr, max_nsw , max_ntrcr
      use icepack_intfc, only: icepack_max_nbtrcr, icepack_max_algae, icepack_max_aero, &
           icepack_nmodal1, icepack_nmodal2

      implicit none

      ! icepack_atmo.F90
      real (kind=dbl_kind), public, &
         dimension (nx_block,ny_block,max_blocks) :: &
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
      real (kind=dbl_kind), public :: &
         hin_max(0:ncat) ! category limits (m)

      character (len=35), public :: c_hi_range(ncat)

      ! icepack_meltpond_lvl.F90
      real (kind=dbl_kind), public, &
         dimension (nx_block,ny_block,ncat,max_blocks) :: &
         dhsn, &      ! depth difference for snow on sea ice and pond ice
         ffracn       ! fraction of fsurfn used to melt ipond

      ! icepack_shortwave.F90
      ! category albedos
      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,ncat,max_blocks), public :: &
         alvdrn      , & ! visible direct albedo           (fraction)
         alidrn      , & ! near-ir direct albedo           (fraction)
         alvdfn      , & ! visible diffuse albedo          (fraction)
         alidfn          ! near-ir diffuse albedo          (fraction)

      ! albedo components for history
      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,ncat,max_blocks), public :: &
         albicen, &   ! bare ice 
         albsnon, &   ! snow 
         albpndn, &   ! pond 
         apeffn       ! effective pond area used for radiation calculation

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat,max_blocks), &
         public :: &
         snowfracn    ! Category snow fraction used in radiation

      ! shortwave components
      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,nilyr,ncat,max_blocks), public :: &
         Iswabsn         ! SW radiation absorbed in ice layers (W m-2)

      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,nslyr,ncat,max_blocks), public :: &
         Sswabsn         ! SW radiation absorbed in snow layers (W m-2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat,max_blocks), &
         public :: &
         fswsfcn     , & ! SW absorbed at ice/snow surface (W m-2)
         fswthrun    , & ! SW through ice to ocean            (W/m^2)
         fswintn         ! SW absorbed in ice interior, below surface (W m-2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr+1,ncat,max_blocks), &
         public :: &
         fswpenln        ! visible SW entering ice layers (W m-2)

      ! aerosol optical properties   -> band  |
      !                                       v aerosol
      ! for combined dust category, use category 4 properties
      real (kind=dbl_kind), dimension(icepack_nspint,icepack_max_aero), public :: & 
         kaer_tab, & ! aerosol mass extinction cross section (m2/kg)
         waer_tab, & ! aerosol single scatter albedo (fraction)
         gaer_tab    ! aerosol asymmetry parameter (cos(theta))

      real (kind=dbl_kind), dimension(icepack_nspint,icepack_nmodal1), public :: & 
         kaer_bc_tab, & ! BC mass extinction cross section (m2/kg)
         waer_bc_tab, & ! BC single scatter albedo (fraction)
         gaer_bc_tab    ! BC aerosol asymmetry parameter (cos(theta))

      real (kind=dbl_kind), dimension (icepack_nspint,icepack_nmodal1,icepack_nmodal2), public :: &
          bcenh           ! BC absorption enhancement factor

      ! biogeochemistry components

      real (kind=dbl_kind), dimension (nblyr+2), public :: &
         bgrid              ! biology nondimensional vertical grid points

      real (kind=dbl_kind), dimension (nblyr+1), public :: &
         igrid              ! biology vertical interface points
 
      real (kind=dbl_kind), dimension (nilyr+1), public :: &
         cgrid            , &  ! CICE vertical coordinate   
         icgrid           , &  ! interface grid for CICE (shortwave variable)
         swgrid                ! grid for ice tracers used in dEdd scheme

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat,max_blocks), public :: &
         first_ice_real     ! .true. = c1, .false. = c0

      logical (kind=log_kind), &
         dimension (nx_block,ny_block,ncat,max_blocks), public :: &
         first_ice      ! distinguishes ice that disappears (e.g. melts)
                        ! and reappears (e.g. transport) in a grid cell
                        ! during a single time step from ice that was
                        ! there the entire time step (true until ice forms)

      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,icepack_max_nbtrcr,max_blocks), public :: &
         ocean_bio      ! contains all the ocean bgc tracer concentrations

      ! diagnostic fluxes
      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,icepack_max_nbtrcr,max_blocks), public :: &
         fbio_snoice, & ! fluxes from snow to ice
         fbio_atmice    ! fluxes from atm to ice

      real (kind=dbl_kind), dimension (nx_block,ny_block,icepack_max_nbtrcr, max_blocks), public :: &
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

      integer (kind=int_kind), dimension(nx_block, ny_block,icepack_max_algae, max_blocks), public :: &
        algal_peak          ! vertical location of algal maximum, 0 if no maximum 

      real (kind=dbl_kind), & 
         dimension (nx_block,ny_block,nblyr+1,ncat,max_blocks), public :: &
         Zoo        ! N losses accumulated in timestep (ie. zooplankton/bacteria)
                    ! mmol/m^3

      real (kind=dbl_kind), &  
         dimension (nx_block,ny_block,ncat,max_blocks), public :: &
         dhbr_top     , & ! brine top change
         dhbr_bot         ! brine bottom change

      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,max_blocks), public :: &
         grow_net       , & ! Specific growth rate (/s) per grid cell
         PP_net         , & ! Total production (mg C/m^2/s) per grid cell
         hbri               ! brine height, area-averaged for comparison with hi (m)

      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,nblyr+2,ncat,max_blocks), public :: &
         bphi           , & ! porosity of layers    
         bTiz               ! layer temperatures interpolated on bio grid (C)

      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,ncat,max_blocks), public :: &
         darcy_V            ! darcy velocity positive up (m/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), public :: &
         zsal_tot    , & ! Total ice salinity in per grid cell (g/m^2) 
         chl_net     , & ! Total chla (mg chla/m^2) per grid cell      
         NO_net          ! Total nitrate per grid cell  

      logical (kind=log_kind), dimension (nx_block,ny_block,max_blocks), public :: &
         Rayleigh_criteria    ! .true. means Ra_c was reached   

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), public :: &
         Rayleigh_real        ! .true. = c1, .false. = c0

      real (kind=dbl_kind), & 
         dimension (nx_block,ny_block,ncat,max_blocks), public :: &
         sice_rho     ! avg sea ice density  (kg/m^3)  ! ech: diagnostic only?

      real (kind=dbl_kind), & 
         dimension (nx_block,ny_block,ncat,max_blocks), public :: &
         fzsaln, &    ! category fzsal(kg/m^2/s) 
         fzsaln_g     ! salt flux from gravity drainage only

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), public :: &
         fzsal    , & ! Total flux  of salt to ocean at time step for conservation
         fzsal_g      ! Total gravity drainage flux

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+1,ncat,max_blocks), public :: &
         zfswin       ! Shortwave flux into layers interpolated on bio grid  (W/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nblyr+1,ncat,max_blocks), public :: &
         iDi      , & ! igrid Diffusivity (m^2/s)    
         iki          ! Ice permeability (m^2)     

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), public :: &
         upNO     , & ! nitrate uptake rate (mmol/m^2/d) times aice
         upNH         ! ammonium uptake rate (mmol/m^2/d) times aice
        
      real (kind=dbl_kind), &
         dimension(nx_block,ny_block,max_ntrcr, ncat, max_blocks), public :: &
         trcrn_sw        ! bgc tracers active in the delta-Eddington shortwave 
                         ! calculation on the shortwave grid (swgrid)

      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,icepack_max_nbtrcr, max_blocks), public :: &
         ice_bio_net  , &   ! depth integrated tracer (mmol/m^2) 
         snow_bio_net       ! depth integrated snow tracer (mmol/m^2)

      logical (kind=log_kind), public :: &
         oceanmixed_ice, &  ! if true, use internal ocean mixed layer
         restore_bgc        ! 

      character(char_len), public :: &
         sil_data_type  , & ! 'default', 'clim'
         nit_data_type  , & ! 'default', 'clim'
         fe_data_type   , & ! 'default', 'clim'
         bgc_data_dir   ! directory for biogeochemistry data

      real (kind=dbl_kind), dimension(icepack_max_algae) :: &
         R_C2N     ,      & ! algal C to N (mole/mole)
         R_chl2N            ! 3 algal chlorophyll to N (mg/mmol)

!=======================================================================

      end module ice_arrays_column

!=======================================================================
