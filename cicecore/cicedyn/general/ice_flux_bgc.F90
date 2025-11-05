!=======================================================================

! Flux variable declarations for biogeochemistry
!
! author Elizabeth C. Hunke, LANL
!
      module ice_flux_bgc

      use ice_kinds_mod
      use ice_constants, only: c0
      use ice_blocks, only: nx_block, ny_block
      use ice_domain_size, only: max_blocks, ncat
      use ice_fileunits, only: nu_diag
      use ice_exit, only: abort_ice
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_max_iso, icepack_max_aero, icepack_max_nbtrcr, &
          icepack_max_algae, icepack_max_doc, icepack_max_don, icepack_max_dic, icepack_max_fe, &
          icepack_query_tracer_indices, icepack_query_tracer_flags, icepack_query_parameters

      implicit none
      private

      public :: bgcflux_ice_to_ocn, alloc_flux_bgc

      ! in from atmosphere

      real (kind=dbl_kind), &   ! coupling variable for both tr_aero and tr_zaero
         dimension (:,:,:,:), allocatable, public :: &
         fiso_atm, & ! isotope deposition rate (kg/m^2 s)
         faero_atm   ! aerosol deposition rate (kg/m^2 s)

      real (kind=dbl_kind), &
         dimension (:,:,:,:), allocatable, public :: &
         flux_bio_atm  ! all bio fluxes to ice from atmosphere

      ! out to ocean

      real (kind=dbl_kind), &
         dimension (:,:,:,:), allocatable, public :: &
         fiso_ocn, & ! isotope flux to ocean  (kg/m^2/s)
         faero_ocn   ! aerosol flux to ocean  (kg/m^2/s)

      real (kind=dbl_kind), &
         dimension (:,:,:,:), allocatable, public :: &
         flux_bio   , & ! all bio fluxes to ocean
         flux_bio_ai    ! all bio fluxes to ocean, averaged over grid cell

      ! internal

      logical (kind=log_kind), public :: &
         cpl_bgc         ! switch to couple BGC via drivers

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         hin_old     , & ! old ice thickness
         dsnown          ! change in snow thickness in category n (m)

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         nit        , & ! ocean nitrate (mmol/m^3)
         amm        , & ! ammonia/um (mmol/m^3)
         sil        , & ! silicate (mmol/m^3)
         dmsp       , & ! dmsp (mmol/m^3)
         dms        , & ! dms (mmol/m^3)
         hum        , & ! humic material carbon (mmol/m^3)
         fnit       , & ! ice-ocean nitrate flux (mmol/m^2/s), positive to ocean
         famm       , & ! ice-ocean ammonia/um flux (mmol/m^2/s), positive to ocean
         fsil       , & ! ice-ocean silicate flux (mmol/m^2/s), positive to ocean
         fdmsp      , & ! ice-ocean dmsp (mmol/m^2/s), positive to ocean
         fdms       , & ! ice-ocean dms (mmol/m^2/s), positive to ocean
         fhum       , & ! ice-ocean humic material carbon (mmol/m^2/s), positive to ocean
         fdust          ! ice-ocean dust flux (kg/m^2/s), positive to ocean

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         algalN     , & ! ocean algal nitrogen (mmol/m^3) (diatoms, pico, phaeo)
         falgalN        ! ice-ocean algal nitrogen flux (mmol/m^2/s) (diatoms, pico, phaeo)

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         doc         , & ! ocean doc (mmol/m^3)  (saccharids, lipids, tbd )
         fdoc            ! ice-ocean doc flux (mmol/m^2/s)  (saccharids, lipids, tbd)

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         don         , & ! ocean don (mmol/m^3) (proteins and amino acids)
         fdon            ! ice-ocean don flux (mmol/m^2/s) (proteins and amino acids)

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         dic         , & ! ocean dic (mmol/m^3)
         fdic            ! ice-ocean dic flux (mmol/m^2/s)

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         fed, fep    , & ! ocean dissolved and particulate fe (nM)
         ffed, ffep      ! ice-ocean dissolved and particulate fe flux (umol/m^2/s)

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         zaeros          ! ocean aerosols (mmol/m^3)

      ! isotopes
      real (kind=dbl_kind), &   ! coupling variable for tr_iso
         dimension (:,:,:,:), allocatable, public :: &
         fiso_evap , & ! isotope evaporation rate (kg/m^2 s)
         Qa_iso    , & ! isotope specific humidity (kg/kg)
         Qref_iso      ! 2m atm reference isotope spec humidity (kg/kg)

      real (kind=dbl_kind), &   ! coupling variable for tr_iso
         dimension (:,:,:), allocatable, public :: &
         HDO_ocn   , & ! seawater concentration of HDO (kg/kg)
         H2_16O_ocn, & ! seawater concentration of H2_16O (kg/kg)
         H2_18O_ocn    ! seawater concentration of H2_18O (kg/kg)

!=======================================================================

      contains

!=======================================================================
!
! Allocate space for all variables
!
      subroutine alloc_flux_bgc

      integer (int_kind) :: ierr

      allocate( &
         nit         (nx_block,ny_block,max_blocks), & ! ocean nitrate (mmol/m^3)
         amm         (nx_block,ny_block,max_blocks), & ! ammonia/um (mmol/m^3)
         sil         (nx_block,ny_block,max_blocks), & ! silicate (mmol/m^3)
         dmsp        (nx_block,ny_block,max_blocks), & ! dmsp (mmol/m^3)
         dms         (nx_block,ny_block,max_blocks), & ! dms (mmol/m^3)
         hum         (nx_block,ny_block,max_blocks), & ! humic material carbon (mmol/m^3)
         fnit        (nx_block,ny_block,max_blocks), & ! ice-ocean nitrate flux (mmol/m^2/s), positive to ocean
         famm        (nx_block,ny_block,max_blocks), & ! ice-ocean ammonia/um flux (mmol/m^2/s), positive to ocean
         fsil        (nx_block,ny_block,max_blocks), & ! ice-ocean silicate flux (mmol/m^2/s), positive to ocean
         fdmsp       (nx_block,ny_block,max_blocks), & ! ice-ocean dmsp (mmol/m^2/s), positive to ocean
         fdms        (nx_block,ny_block,max_blocks), & ! ice-ocean dms (mmol/m^2/s), positive to ocean
         fhum        (nx_block,ny_block,max_blocks), & ! ice-ocean humic material carbon (mmol/m^2/s), positive to ocean
         fdust       (nx_block,ny_block,max_blocks), & ! ice-ocean dust flux (kg/m^2/s), positive to ocean
         hin_old     (nx_block,ny_block,ncat,max_blocks), & ! old ice thickness
         dsnown      (nx_block,ny_block,ncat,max_blocks), & ! change in snow thickness in category n (m)
         HDO_ocn     (nx_block,ny_block,max_blocks), & ! seawater concentration of HDO (kg/kg)
         H2_16O_ocn  (nx_block,ny_block,max_blocks), & ! seawater concentration of H2_16O (kg/kg)
         H2_18O_ocn  (nx_block,ny_block,max_blocks), & ! seawater concentration of H2_18O (kg/kg)
         Qa_iso      (nx_block,ny_block,icepack_max_iso,max_blocks), & ! isotope specific humidity (kg/kg)
         Qref_iso    (nx_block,ny_block,icepack_max_iso,max_blocks), & ! 2m atm reference isotope spec humidity (kg/kg)
         fiso_atm    (nx_block,ny_block,icepack_max_iso,max_blocks), & ! isotope deposition rate (kg/m^2 s)
         fiso_evap   (nx_block,ny_block,icepack_max_iso,max_blocks), & ! isotope evaporation rate (kg/m^2 s)
         fiso_ocn    (nx_block,ny_block,icepack_max_iso,max_blocks), & ! isotope flux to ocean  (kg/m^2/s)
         faero_atm   (nx_block,ny_block,icepack_max_aero,max_blocks), & ! aerosol deposition rate (kg/m^2 s)
         faero_ocn   (nx_block,ny_block,icepack_max_aero,max_blocks), & ! aerosol flux to ocean  (kg/m^2/s)
         zaeros      (nx_block,ny_block,icepack_max_aero,max_blocks), & ! ocean aerosols (mmol/m^3)
         flux_bio_atm(nx_block,ny_block,icepack_max_nbtrcr,max_blocks), & ! all bio fluxes to ice from atmosphere
         flux_bio    (nx_block,ny_block,icepack_max_nbtrcr,max_blocks), & ! all bio fluxes to ocean
         flux_bio_ai (nx_block,ny_block,icepack_max_nbtrcr,max_blocks), & ! all bio fluxes to ocean, averaged over grid cell
         algalN      (nx_block,ny_block,icepack_max_algae,max_blocks), & ! ocean algal nitrogen (mmol/m^3) (diatoms, pico, phaeo)
         falgalN     (nx_block,ny_block,icepack_max_algae,max_blocks), & ! ice-ocn algalN flux (mmol/m^2/s) (diatoms, pico, phaeo)
         doc         (nx_block,ny_block,icepack_max_doc,max_blocks), & ! ocean doc (mmol/m^3)  (saccharids, lipids, tbd )
         fdoc        (nx_block,ny_block,icepack_max_doc,max_blocks), & ! ice-ocean doc flux (mmol/m^2/s)  (saccharids, lipids, tbd)
         don         (nx_block,ny_block,icepack_max_don,max_blocks), & ! ocean don (mmol/m^3) (proteins and amino acids)
         fdon        (nx_block,ny_block,icepack_max_don,max_blocks), & ! ice-ocean don flux (mmol/m^2/s) (proteins and amino acids)
         dic         (nx_block,ny_block,icepack_max_dic,max_blocks), & ! ocean dic (mmol/m^3)
         fdic        (nx_block,ny_block,icepack_max_dic,max_blocks), & ! ice-ocean dic flux (mmol/m^2/s)
         fed         (nx_block,ny_block,icepack_max_fe, max_blocks), & ! ocean dissolved fe (nM)
         fep         (nx_block,ny_block,icepack_max_fe, max_blocks), & ! ocean particulate fe (nM)
         ffed        (nx_block,ny_block,icepack_max_fe, max_blocks), & ! ice-ocean dissolved fe flux (umol/m^2/s)
         ffep        (nx_block,ny_block,icepack_max_fe, max_blocks), & ! ice-ocean particulate fe flux (umol/m^2/s)
         stat=ierr)
      if (ierr/=0) call abort_ice('(alloc_flux_bgc): Out of memory')

      nit         = c0
      amm         = c0
      sil         = c0
      dmsp        = c0
      dms         = c0
      hum         = c0
      fnit        = c0
      famm        = c0
      fsil        = c0
      fdmsp       = c0
      fdms        = c0
      fhum        = c0
      fdust       = c0
      hin_old     = c0
      dsnown      = c0
      HDO_ocn     = c0
      H2_16O_ocn  = c0
      H2_18O_ocn  = c0
      Qa_iso      = c0
      Qref_iso    = c0
      fiso_atm    = c0
      fiso_evap   = c0
      fiso_ocn    = c0
      faero_atm   = c0
      faero_ocn   = c0
      zaeros      = c0
      flux_bio_atm= c0
      flux_bio    = c0
      flux_bio_ai = c0
      algalN      = c0
      falgalN     = c0
      doc         = c0
      fdoc        = c0
      don         = c0
      fdon        = c0
      dic         = c0
      fdic        = c0
      fed         = c0
      fep         = c0
      ffed        = c0
      ffep        = c0

      end subroutine alloc_flux_bgc

!=======================================================================

! Initialize some fluxes sent to coupler for use by the atm model
!
! author: Nicole Jeffery, LANL

      subroutine bgcflux_ice_to_ocn(nx_block,       &
                                  ny_block,         &
                                  flux_bio,         &
                                  f_nit,    f_sil,    &
                                  f_amm,    f_dmsp,   &
                                  f_dms,    f_hum,    &
                                  f_dust,   f_algalN, &
                                  f_doc,    f_dic,    &
                                  f_don,    f_fep,    &
                                  f_fed)

      use ice_constants, only: c0
      use ice_domain_size, only: n_zaero, n_algae, n_doc, n_dic, n_don, n_fed, n_fep

      real(kind=dbl_kind), dimension(:,:,:), intent(in) :: &
          flux_bio
      real(kind=dbl_kind), dimension(:,:), intent(out):: &
          f_nit,  &  ! nitrate flux mmol/m^2/s  positive to ocean
          f_sil,  &  ! silicate flux mmol/m^2/s
          f_amm,  &  ! ammonium flux mmol/m^2/s
          f_dmsp, &  ! DMSPd flux mmol/m^2/s
          f_dms,  &  ! DMS flux mmol/m^2/s
          f_hum,  &  ! humic flux mmol/m^2/s
          f_dust     ! dust flux kg/m^2/s

      real(kind=dbl_kind), dimension(:,:,:), intent(out):: &
          f_algalN, & ! algal nitrogen flux mmol/m^2/s
          f_doc,    & ! DOC flux mmol/m^2/s
          f_dic,    & ! DIC flux mmol/m^2/s
          f_don,    & ! DON flux mmol/m^2/s
          f_fep,    & ! particulate iron flux umol/m^2/s
          f_fed       ! dissolved iron flux umol/m^2/s

      integer (kind=int_kind), intent(in) :: &
          nx_block, &
          ny_block

      ! local variables

      integer (kind=int_kind) :: &
         i,j         , & ! horizontal indices
         k               ! tracer index

      logical (kind=log_kind) :: &
          skl_bgc, solve_zbgc, &
          tr_bgc_Nit, tr_bgc_N, &
          tr_bgc_DON, tr_bgc_C, tr_bgc_Am, tr_bgc_Sil, tr_bgc_DMS, tr_bgc_Fe, &
          tr_bgc_hum, tr_zaero

      integer (kind=int_kind) :: &
          nlt_bgc_Nit, nlt_bgc_Am, &
          nlt_bgc_Sil, nlt_bgc_DMSPd, nlt_bgc_DMS, nlt_bgc_hum

      integer (kind=int_kind), dimension(icepack_max_algae) :: &
          nlt_bgc_N, nlt_bgc_C   ! algae
      integer (kind=int_kind), dimension(icepack_max_doc) :: &
          nlt_bgc_DOC            ! disolved organic carbon
      integer (kind=int_kind), dimension(icepack_max_don) :: &
          nlt_bgc_DON            !
      integer (kind=int_kind), dimension(icepack_max_dic) :: &
          nlt_bgc_DIC            ! disolved inorganic carbon
      integer (kind=int_kind), dimension(icepack_max_fe) :: &
          nlt_bgc_Fed, nlt_bgc_Fep  !
      integer (kind=int_kind), dimension(icepack_max_aero) :: &
          nlt_zaero              ! non-reacting layer aerosols

      character(len=*), parameter :: subname = '(bgcflux_ice_to_ocn)'

      call icepack_query_parameters(skl_bgc_out=skl_bgc, solve_zbgc_out=solve_zbgc)
      call icepack_query_tracer_flags( &
          tr_bgc_Nit_out=tr_bgc_Nit, tr_bgc_N_out=tr_bgc_N, &
          tr_bgc_DON_out=tr_bgc_DON, tr_bgc_C_out=tr_bgc_C, tr_bgc_Am_out=tr_bgc_Am, &
          tr_bgc_Sil_out=tr_bgc_Sil, tr_bgc_DMS_out=tr_bgc_DMS, tr_bgc_Fe_out=tr_bgc_Fe, &
          tr_bgc_hum_out=tr_bgc_hum, tr_zaero_out=tr_zaero)
      call icepack_query_tracer_indices( &
          nlt_bgc_N_out=nlt_bgc_N, nlt_bgc_C_out=nlt_bgc_C, nlt_bgc_DOC_out=nlt_bgc_DOC, &
          nlt_bgc_DON_out=nlt_bgc_DON, nlt_bgc_DIC_out=nlt_bgc_DIC, &
          nlt_bgc_Fed_out=nlt_bgc_Fed, nlt_bgc_Fep_out=nlt_bgc_Fep, &
          nlt_zaero_out=nlt_zaero, nlt_bgc_Nit_out=nlt_bgc_Nit, nlt_bgc_Am_out=nlt_bgc_Am, &
          nlt_bgc_Sil_out=nlt_bgc_Sil, nlt_bgc_DMSPd_out=nlt_bgc_DMSPd, &
          nlt_bgc_DMS_out=nlt_bgc_DMS, nlt_bgc_hum_out=nlt_bgc_hum)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      f_nit    (:,:) = c0
      f_sil    (:,:) = c0
      f_amm    (:,:) = c0
      f_dmsp   (:,:) = c0
      f_dms    (:,:) = c0
      f_hum    (:,:) = c0
      f_dust   (:,:) = c0
      f_algalN(:,:,:)= c0
      f_doc   (:,:,:)= c0
      f_dic   (:,:,:)= c0
      f_don   (:,:,:)= c0
      f_fep   (:,:,:)= c0
      f_fed   (:,:,:)= c0

      do j = 1, ny_block
      do i = 1, nx_block
         if (skl_bgc .or. solve_zbgc) then
            do k = 1, n_algae
               f_algalN(i,j,k) = flux_bio(i,j,nlt_bgc_N(k))
            enddo
         endif
         if (tr_bgc_C) then
            do k = 1, n_doc
               f_doc(i,j,k) = flux_bio(i,j,nlt_bgc_DOC(k))
            enddo
            do k = 1, n_dic
               f_dic(i,j,k) = flux_bio(i,j,nlt_bgc_DIC(k))
            enddo
         endif
         if (tr_bgc_DON) then
            do k = 1, n_don
               f_don(i,j,k) = flux_bio(i,j,nlt_bgc_DON(k))
            enddo
         endif
         if (tr_bgc_Fe) then
            do k = 1, n_fep
               f_fep(i,j,k) = flux_bio(i,j,nlt_bgc_Fep(k))
            enddo
            do k = 1, n_fed
               f_fed(i,j,k) = flux_bio(i,j,nlt_bgc_Fed(k))
            enddo
         endif
         if (tr_bgc_Nit) f_nit(i,j)  = flux_bio(i,j,nlt_bgc_Nit)
         if (tr_bgc_Sil) f_sil(i,j)  = flux_bio(i,j,nlt_bgc_Sil)
         if (tr_bgc_Am)  f_amm(i,j)  = flux_bio(i,j,nlt_bgc_Am)
         if (tr_bgc_hum) f_hum(i,j)  = flux_bio(i,j,nlt_bgc_hum)
         if (tr_bgc_DMS) then
            f_dms(i,j) = flux_bio(i,j,nlt_bgc_DMS)
            f_dmsp(i,j) = flux_bio(i,j,nlt_bgc_DMSPd)
         endif
         if (tr_zaero) then
            do k = 3, n_zaero
                f_dust(i,j) = f_dust(i,j) + flux_bio(i,j,nlt_zaero(k))
            enddo
         endif
      enddo
      enddo

      end subroutine bgcflux_ice_to_ocn

!=======================================================================

      end module ice_flux_bgc

!=======================================================================
