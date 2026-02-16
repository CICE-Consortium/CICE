!=======================================================================
!
! Primary state variables in various configurations
! Note: other state variables are at the end of this...
! The primary state variable names are:
!-------------------------------------------------------------------
! for each category   aggregated over     units
!                       categories
!-------------------------------------------------------------------
! aicen(i,j,n)         aice(i,j)           ---
! vicen(i,j,n)         vice(i,j)           m
! vsnon(i,j,n)         vsno(i,j)           m
! trcrn(i,j,it,n)      trcr(i,j,it)
!
! Area is dimensionless because aice is the fractional area
! (normalized so that the sum over all categories, including open
! water, is 1.0).  That is why vice/vsno have units of m instead of m^3.
!
! Variable names follow these rules:
!
! (1) For 3D variables (indices i,j,n), write 'ice' or 'sno' or
!     'sfc' and put an 'n' at the end.
! (2) For 2D variables (indices i,j) aggregated over all categories,
!     write 'ice' or 'sno' or 'sfc' without the 'n'.
! (3) For 2D variables (indices i,j) associated with an individual
!     category, write 'i' or 's' instead of 'ice' or 'sno' and put an 'n'
!     at the end: e.g. hin, hsn.  These are not declared here
!     but in individual modules (e.g., icepack_therm_vertical).
!
! authors C. M. Bitz, UW
!         Elizabeth C. Hunke and William H. Lipscomb, LANL
!
! 2004: Block structure added by William Lipscomb
! 2006: Converted to free form source (F90) by Elizabeth Hunke

      module ice_state

      use ice_kinds_mod
      use ice_constants, only: field_loc_center, field_type_scalar, c0
      use ice_domain_size, only: max_blocks, ncat
      use ice_blocks, only: nx_block, ny_block
      use ice_exit, only: abort_ice
      use ice_fileunits, only: nu_diag
      use icepack_intfc, only: icepack_query_tracer_sizes
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted

      implicit none
      private
      public :: bound_state, alloc_state

      !-----------------------------------------------------------------
      ! state of the ice aggregated over all categories
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension(:,:,:), allocatable, public :: &
         aice  , & ! concentration of ice on T grid
         aiU   , & ! concentration of ice on U grid
         vice  , & ! volume per unit area of ice          (m)
         vsno      ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension(:,:,:,:), allocatable, public :: &
         trcr      ! ice tracers
                   ! 1: surface temperature of ice/snow (C)

      !-----------------------------------------------------------------
      ! state of the ice for each category
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public:: &
         aice0     ! concentration of open water

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (:,:,:,:,:), allocatable, public :: &
         trcrn     ! tracers
                   ! 1: surface temperature of ice/snow (C)

      !-----------------------------------------------------------------
      ! tracers infrastructure arrays
      !-----------------------------------------------------------------

      integer (kind=int_kind), dimension (:), allocatable, public :: &
         trcr_depend   ! = 0 for ice area tracers
                       ! = 1 for ice volume tracers
                       ! = 2 for snow volume tracers

      integer (kind=int_kind), dimension (:), allocatable, public :: &
         n_trcr_strata ! number of underlying tracer layers

      integer (kind=int_kind), dimension (:,:), allocatable, public :: &
         nt_strata     ! indices of underlying tracer layers

      real (kind=dbl_kind), dimension (:,:), allocatable, public :: &
         trcr_base     ! = 0 or 1 depending on tracer dependency
                       ! argument 2:  (1) aice, (2) vice, (3) vsno

      !-----------------------------------------------------------------
      ! dynamic variables closely related to the state of the ice
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension(:,:,:), allocatable, public :: &
         uvel     , & ! x-component of velocity on U grid (m/s)
         vvel     , & ! y-component of velocity on U grid (m/s)
         uvelE    , & ! x-component of velocity on E grid (m/s)
         vvelE    , & ! y-component of velocity on E grid (m/s)
         uvelN    , & ! x-component of velocity on N grid (m/s)
         vvelN    , & ! y-component of velocity on N grid (m/s)
         divu     , & ! strain rate I component, velocity divergence (1/s)
         shear    , & ! strain rate II component (1/s)
         vort     , & ! vorticity (1/s)
         strength     ! ice strength (N/m)

      !-----------------------------------------------------------------
      ! ice state at start of time step, saved for later in the step
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension(:,:,:), allocatable, public :: &
         aice_init       ! initial concentration of ice, for diagnostics

      real (kind=dbl_kind), dimension(:,:,:,:), allocatable, public :: &
         aicen_init  , & ! initial ice concentration, for linear ITD
         vicen_init  , & ! initial ice volume (m), for linear ITD
         vsnon_init  , & ! initial snow volume (m), for aerosol
         Tsfcn_init      ! initial ice surface temperature (degC)

!=======================================================================

      contains

!=======================================================================
!
! Allocate space for all state variables
!
      subroutine alloc_state
      integer (int_kind) :: ntrcr, ierr
      character(len=*),parameter :: subname='(alloc_state)'

      call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__, line=__LINE__)

      allocate ( &
         aice      (nx_block,ny_block,max_blocks) , & ! concentration of ice T grid
         aiU       (nx_block,ny_block,max_blocks) , & ! concentration of ice U grid
         vice      (nx_block,ny_block,max_blocks) , & ! volume per unit area of ice (m)
         vsno      (nx_block,ny_block,max_blocks) , & ! volume per unit area of snow (m)
         aice0     (nx_block,ny_block,max_blocks) , & ! concentration of open water
         uvel      (nx_block,ny_block,max_blocks) , & ! x-component of velocity on U grid (m/s)
         vvel      (nx_block,ny_block,max_blocks) , & ! y-component of velocity on U grid (m/s)
         uvelE     (nx_block,ny_block,max_blocks) , & ! x-component of velocity on E grid (m/s)
         vvelE     (nx_block,ny_block,max_blocks) , & ! y-component of velocity on E grid (m/s)
         uvelN     (nx_block,ny_block,max_blocks) , & ! x-component of velocity on N grid (m/s)
         vvelN     (nx_block,ny_block,max_blocks) , & ! y-component of velocity on N grid (m/s)
         divu      (nx_block,ny_block,max_blocks) , & ! strain rate I component, velocity divergence (1/s)
         shear     (nx_block,ny_block,max_blocks) , & ! strain rate II component (1/s)
         vort      (nx_block,ny_block,max_blocks) , & ! vorticity (1/s)
         strength  (nx_block,ny_block,max_blocks) , & ! ice strength (N/m)
         aice_init (nx_block,ny_block,max_blocks) , & ! initial concentration of ice, for diagnostics
         aicen     (nx_block,ny_block,ncat,max_blocks) , & ! concentration of ice
         vicen     (nx_block,ny_block,ncat,max_blocks) , & ! volume per unit area of ice (m)
         vsnon     (nx_block,ny_block,ncat,max_blocks) , & ! volume per unit area of snow (m)
         aicen_init(nx_block,ny_block,ncat,max_blocks) , & ! initial ice concentration, for linear ITD
         vicen_init(nx_block,ny_block,ncat,max_blocks) , & ! initial ice volume (m), for linear ITD
         vsnon_init(nx_block,ny_block,ncat,max_blocks) , & ! initial snow volume (m), for aerosol
         Tsfcn_init(nx_block,ny_block,ncat,max_blocks) , & ! initial snow/ice surface temperature(degC)
         trcr      (nx_block,ny_block,ntrcr,max_blocks) , & ! ice tracers: 1: surface temperature of ice/snow (C)
         trcrn     (nx_block,ny_block,ntrcr,ncat,max_blocks) , & ! tracers: 1: surface temperature of ice/snow (C)
         stat=ierr)
      if (ierr/=0) call abort_ice('(alloc_state): Out of memory1')

      aice  = c0
      aiU   = c0
      vice  = c0
      vsno  = c0
      aice0 = c0
      uvel  = c0
      vvel  = c0
      uvelE = c0
      vvelE = c0
      uvelN = c0
      vvelN = c0
      divu  = c0
      shear = c0
      vort  = c0
      strength   = c0
      aice_init  = c0
      aicen = c0
      vicen = c0
      vsnon = c0
      aicen_init = c0
      vicen_init = c0
      vsnon_init = c0
      Tsfcn_init = c0
      trcr  = c0
      trcrn = c0

      allocate ( &
         trcr_depend(ntrcr)   , & !
         n_trcr_strata(ntrcr) , & ! number of underlying tracer layers
         nt_strata(ntrcr,2)   , & ! indices of underlying tracer layers
         trcr_base(ntrcr,3)   , & ! = 0 or 1 depending on tracer dependency, (1) aice, (2) vice, (3) vsno
         stat=ierr)
      if (ierr/=0) call abort_ice('(alloc_state): Out of memory2')

      trcr_depend = 0
      n_trcr_strata = 0
      nt_strata = 0
      trcr_base = c0

      end subroutine alloc_state

!=======================================================================
!
! Get ghost cell values for ice state variables in each thickness category.
! NOTE: This subroutine cannot be called from inside a block loop!
!
! author: William H. Lipscomb, LANL

      subroutine bound_state (aicen,        &
                              vicen, vsnon, &
                              ntrcr, trcrn)

      use ice_boundary, only: ice_halo, ice_HaloMask, ice_HaloUpdate, &
          ice_HaloDestroy
      use ice_domain, only: halo_info, maskhalo_bound, nblocks

      integer (kind=int_kind), intent(in) :: &
         ntrcr     ! number of tracers in use

      real (kind=dbl_kind), dimension(nx_block,ny_block,ncat,max_blocks), intent(inout) :: &
         aicen , & ! fractional ice area
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), intent(inout), dimension(:,:,:,:,:) :: &  ! (nx_block,ny_block,ntrcr,ncat,max_blocks)
         trcrn     ! ice tracers

      ! local variables

      integer (kind=int_kind) :: i, j, n, iblk

      integer (kind=int_kind), &
         dimension(nx_block,ny_block,max_blocks) :: halomask

      type (ice_halo) :: halo_info_aicemask

      character(len=*), parameter :: subname = '(bound_state)'

      call ice_HaloUpdate (aicen,            halo_info, &
                           field_loc_center, field_type_scalar)

      if (maskhalo_bound) then
         halomask(:,:,:) = 0

         !$OMP PARALLEL DO PRIVATE(iblk,n,i,j)
         do iblk = 1, nblocks
         do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block
            if (aicen(i,j,n,iblk) > c0) halomask(i,j,iblk) = 1
         enddo
         enddo
         enddo
         enddo
         !$OMP END PARALLEL DO

         call ice_HaloMask(halo_info_aicemask, halo_info, halomask)

         call ice_HaloUpdate (trcrn(:,:,:,:,:), halo_info_aicemask, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (vicen,            halo_info_aicemask, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (vsnon,            halo_info_aicemask, &
                              field_loc_center, field_type_scalar)
         call ice_HaloDestroy(halo_info_aicemask)

      else
         call ice_HaloUpdate (trcrn(:,:,:,:,:), halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (vicen,            halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (vsnon,            halo_info, &
                              field_loc_center, field_type_scalar)
      endif

      end subroutine bound_state

!=======================================================================

      end module ice_state

!=======================================================================
