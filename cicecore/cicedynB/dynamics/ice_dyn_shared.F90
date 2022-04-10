!=======================================================================

! Elastic-viscous-plastic sea ice dynamics model code shared with other
! approaches
!
! author: Elizabeth C. Hunke, LANL
!
! 2013: Split from ice_dyn_evp.F90 by Elizabeth Hunke

      module ice_dyn_shared

      use ice_kinds_mod
      use ice_communicate, only: my_task, master_task, get_num_procs
      use ice_constants, only: c0, c1, c2, c3, c4, c6
      use ice_constants, only: omega, spval_dbl, p01, p001, p5
      use ice_blocks, only: nx_block, ny_block
      use ice_domain_size, only: max_blocks
      use ice_fileunits, only: nu_diag
      use ice_exit, only: abort_ice
      use ice_grid, only: grid_ice
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_parameters

      implicit none
      private
      public :: set_evp_parameters, stepu, stepuv_CD, stepu_C, stepv_C, &
                principal_stress, init_dyn, dyn_prep1, dyn_prep2, dyn_finish, &
                seabed_stress_factor_LKD, seabed_stress_factor_prob, &
                alloc_dyn_shared, &
                deformations, deformationsC_T, deformationsCD_T, &
                strain_rates, strain_rates_T, strain_rates_U, &
                visc_replpress, &
                dyn_haloUpdate, &
                stack_fields, unstack_fields

      ! namelist parameters

      integer (kind=int_kind), public :: &
         kdyn       , & ! type of dynamics ( -1, 0 = off, 1 = evp, 2 = eap )
         kridge     , & ! set to "-1" to turn off ridging
         ndte           ! number of subcycles:  ndte=dt/dte

      character (len=char_len), public :: &
         coriolis   , & ! 'constant', 'zero', or 'latitude'
         ssh_stress     ! 'geostrophic' or 'coupled'

      logical (kind=log_kind), public :: &
         revised_evp    ! if true, use revised evp procedure

      character (len=char_len), public :: &
         evp_algorithm  ! standard_2d = 2D org version (standard)
                        ! shared_mem_1d = 1d without mpi call and refactorization to 1d

      real (kind=dbl_kind), public :: &
         elasticDamp    ! coefficient for calculating the parameter E, elastic damping parameter

      ! other EVP parameters

      character (len=char_len), public :: &
         yield_curve      , & ! 'ellipse' ('teardrop' needs further testing)
         visc_method      , & ! method for viscosity calc at U points (C, CD grids)
         seabed_stress_method ! method for seabed stress calculation
                              ! LKD: Lemieux et al. 2015, probabilistic: Dupont et al. in prep.

      real (kind=dbl_kind), parameter, public :: &
         u0    = 5e-5_dbl_kind, & ! residual velocity for seabed stress (m/s)
         cosw  = c1           , & ! cos(ocean turning angle)  ! turning angle = 0
         sinw  = c0           , & ! sin(ocean turning angle)  ! turning angle = 0
         a_min = p001         , & ! minimum ice area
         m_min = p01              ! minimum ice mass (kg/m^2)

      real (kind=dbl_kind), public :: &
         revp        , & ! 0 for classic EVP, 1 for revised EVP
         e_yieldcurve, & ! VP aspect ratio of elliptical yield curve
         e_plasticpot, & ! VP aspect ratio of elliptical plastic potential
         epp2i       , & ! 1/(e_plasticpot)^2
         e_factor    , & ! (e_yieldcurve)^2/(e_plasticpot)^4
         ecci        , & ! temporary for 1d evp
         deltaminEVP , & ! minimum delta for viscosities (EVP)
         deltaminVP  , & ! minimum delta for viscosities (VP)
         capping     , & ! capping of viscosities (1=Hibler79, 0=Kreyscher2000)
         dtei        , & ! 1/dte, where dte is subcycling timestep (1/s)
!         dte2T       , & ! dte/2T
         denom1          ! constants for stress equation

      real (kind=dbl_kind), public :: & ! Bouillon et al relaxation constants
         arlx        , & ! alpha for stressp
         arlx1i      , & ! (inverse of alpha) for stressp
         brlx            ! beta   for momentum

      real (kind=dbl_kind), allocatable, public :: &
         fcor_blk(:,:,:)     ! Coriolis parameter (1/s)

      real (kind=dbl_kind), allocatable, public :: &
         fcorE_blk(:,:,:), & ! Coriolis parameter at E points (1/s)
         fcorN_blk(:,:,:)    ! Coriolis parameter at N points  (1/s)

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         uvel_init       , & ! x-component of velocity (m/s), beginning of timestep
         vvel_init           ! y-component of velocity (m/s), beginning of timestep

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         uvelN_init      , & ! x-component of velocity (m/s), beginning of timestep
         vvelN_init          ! y-component of velocity (m/s), beginning of timestep

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         uvelE_init      , & ! x-component of velocity (m/s), beginning of timestep
         vvelE_init          ! y-component of velocity (m/s), beginning of timestep

      real (kind=dbl_kind), allocatable, public :: &
         DminTarea(:,:,:)    ! deltamin * tarea (m^2/s)

      ! ice isotropic tensile strength parameter
      real (kind=dbl_kind), public :: &
         Ktens               ! T=Ktens*P (tensile strength: see Konig and Holland, 2010)

      ! seabed (basal) stress parameters and settings
      logical (kind=log_kind), public :: &
         seabed_stress       ! if true, seabed stress for landfast on

      real (kind=dbl_kind), public :: &
         k1              , & ! 1st free parameter for seabed1 grounding parameterization
         k2              , & ! second free parameter (N/m^3) for seabed1 grounding parametrization
         alphab          , & ! alphab=Cb factor in Lemieux et al 2015
         threshold_hw        ! max water depth for grounding
                             ! see keel data from Amundrud et al. 2004 (JGR)

      interface strain_rates_T
         module procedure strain_rates_Tdt
         module procedure strain_rates_Tdtsd
      end interface

      interface dyn_haloUpdate
         module procedure dyn_haloUpdate1
         module procedure dyn_haloUpdate2
         module procedure dyn_haloUpdate3
         module procedure dyn_haloUpdate4
         module procedure dyn_haloUpdate5
      end interface

      interface stack_fields
         module procedure stack_fields2
         module procedure stack_fields3
         module procedure stack_fields4
         module procedure stack_fields5
      end interface

      interface unstack_fields
         module procedure unstack_fields2
         module procedure unstack_fields3
         module procedure unstack_fields4
         module procedure unstack_fields5
      end interface

!=======================================================================

      contains

!=======================================================================
!
! Allocate space for all variables
!
      subroutine alloc_dyn_shared

      integer (int_kind) :: ierr

      character(len=*), parameter :: subname = '(alloc_dyn_shared)'

      allocate( &
         uvel_init (nx_block,ny_block,max_blocks), & ! x-component of velocity (m/s), beginning of timestep
         vvel_init (nx_block,ny_block,max_blocks), & ! y-component of velocity (m/s), beginning of timestep
         stat=ierr)
      if (ierr/=0) call abort_ice(subname//': Out of memory')

      if (grid_ice == 'CD' .or. grid_ice == 'C') then
         allocate( &
            uvelE_init (nx_block,ny_block,max_blocks), & ! x-component of velocity (m/s), beginning of timestep
            vvelE_init (nx_block,ny_block,max_blocks), & ! y-component of velocity (m/s), beginning of timestep
            uvelN_init (nx_block,ny_block,max_blocks), & ! x-component of velocity (m/s), beginning of timestep
            vvelN_init (nx_block,ny_block,max_blocks), & ! y-component of velocity (m/s), beginning of timestep
            stat=ierr)
         if (ierr/=0) call abort_ice(subname//': Out of memory')
      endif

      end subroutine alloc_dyn_shared

!=======================================================================
! Initialize parameters and variables needed for the dynamics
! author: Elizabeth C. Hunke, LANL

      subroutine init_dyn (dt)

      use ice_blocks, only: nx_block, ny_block
      use ice_domain, only: nblocks, halo_dynbundle
      use ice_domain_size, only: max_blocks
      use ice_flux, only: rdg_conv, rdg_shear, iceumask, &
          stressp_1, stressp_2, stressp_3, stressp_4, &
          stressm_1, stressm_2, stressm_3, stressm_4, &
          stress12_1, stress12_2, stress12_3, stress12_4, &
          stresspT, stressmT, stress12T, &
          stresspU, stressmU, stress12U
      use ice_state, only: uvel, vvel, uvelE, vvelE, uvelN, vvelN, divu, shear
      use ice_grid, only: ULAT, NLAT, ELAT, tarea

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) :: &
         i, j  , &  ! indices
         nprocs, &  ! number of processors
         iblk       ! block index

      character(len=*), parameter :: subname = '(init_dyn)'

      call set_evp_parameters (dt)

      ! Set halo_dynbundle, this is empirical at this point, could become namelist
      halo_dynbundle = .true.
      nprocs = get_num_procs()
      if (nx_block*ny_block/nprocs > 100) halo_dynbundle = .false.

      if (my_task == master_task) then
         write(nu_diag,*) 'dt  = ',dt
         write(nu_diag,*) 'dte = ',dt/real(ndte,kind=dbl_kind)
         write(nu_diag,*) 'tdamp =', elasticDamp * dt
         write(nu_diag,*) 'halo_dynbundle =', halo_dynbundle
      endif

      allocate(fcor_blk(nx_block,ny_block,max_blocks))
      allocate(DminTarea(nx_block,ny_block,max_blocks))

      if (grid_ice == 'CD' .or. grid_ice == 'C') then
         allocate(fcorE_blk(nx_block,ny_block,max_blocks))
         allocate(fcorN_blk(nx_block,ny_block,max_blocks))
      endif

      !$OMP PARALLEL DO PRIVATE(iblk,i,j) SCHEDULE(runtime)
      do iblk = 1, nblocks
      do j = 1, ny_block
      do i = 1, nx_block

         ! velocity
         uvel(i,j,iblk) = c0    ! m/s
         vvel(i,j,iblk) = c0    ! m/s
         if (grid_ice == 'CD' .or. grid_ice == 'C') then ! extra velocity variables
            uvelE(i,j,iblk) = c0
            vvelE(i,j,iblk) = c0
            uvelN(i,j,iblk) = c0
            vvelN(i,j,iblk) = c0
         endif

         ! strain rates
         divu (i,j,iblk) = c0
         shear(i,j,iblk) = c0
         rdg_conv (i,j,iblk) = c0
         rdg_shear(i,j,iblk) = c0

         ! Coriolis parameter
         if (trim(coriolis) == 'constant') then
            fcor_blk(i,j,iblk) = 1.46e-4_dbl_kind ! Hibler 1979, N. Hem; 1/s
         else if (trim(coriolis) == 'zero') then
            fcor_blk(i,j,iblk) = c0
         else
            fcor_blk(i,j,iblk) = c2*omega*sin(ULAT(i,j,iblk)) ! 1/s
         endif

         if (grid_ice == 'CD' .or. grid_ice == 'C') then

            if (trim(coriolis) == 'constant') then
               fcorE_blk(i,j,iblk) = 1.46e-4_dbl_kind ! Hibler 1979, N. Hem; 1/s
               fcorN_blk(i,j,iblk) = 1.46e-4_dbl_kind ! Hibler 1979, N. Hem; 1/s
            else if (trim(coriolis) == 'zero') then
               fcorE_blk(i,j,iblk) = c0
               fcorN_blk(i,j,iblk) = c0
            else
               fcorE_blk(i,j,iblk) = c2*omega*sin(ELAT(i,j,iblk)) ! 1/s
               fcorN_blk(i,j,iblk) = c2*omega*sin(NLAT(i,j,iblk)) ! 1/s
            endif

         endif

         ! stress tensor,  kg/s^2
         stressp_1 (i,j,iblk) = c0
         stressp_2 (i,j,iblk) = c0
         stressp_3 (i,j,iblk) = c0
         stressp_4 (i,j,iblk) = c0
         stressm_1 (i,j,iblk) = c0
         stressm_2 (i,j,iblk) = c0
         stressm_3 (i,j,iblk) = c0
         stressm_4 (i,j,iblk) = c0
         stress12_1(i,j,iblk) = c0
         stress12_2(i,j,iblk) = c0
         stress12_3(i,j,iblk) = c0
         stress12_4(i,j,iblk) = c0

         if (grid_ice == 'CD' .or. grid_ice == 'C') then
            stresspT  (i,j,iblk) = c0
            stressmT  (i,j,iblk) = c0
            stress12T (i,j,iblk) = c0
            stresspU  (i,j,iblk) = c0
            stressmU  (i,j,iblk) = c0
            stress12U (i,j,iblk) = c0
         endif

         if (kdyn == 1) then
            DminTarea(i,j,iblk) = deltaminEVP*tarea(i,j,iblk)
         elseif (kdyn == 3) then
            DminTarea(i,j,iblk) = deltaminVP*tarea(i,j,iblk)
         endif

         ! ice extent mask on velocity points
         iceumask(i,j,iblk) = .false.

      enddo                     ! i
      enddo                     ! j
      enddo                     ! iblk
      !$OMP END PARALLEL DO

      end subroutine init_dyn

!=======================================================================
! Set parameters needed for the evp dynamics.
! Note: This subroutine is currently called only during initialization.
!       If the dynamics time step can vary during runtime, it should
!        be called whenever the time step changes.
!
! author: Elizabeth C. Hunke, LANL

      subroutine set_evp_parameters (dt)

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      !real (kind=dbl_kind) :: &
         !dte         , & ! subcycling timestep for EVP dynamics, s
         !tdamp2          ! 2*(wave damping time scale T)

      character(len=*), parameter :: subname = '(set_evp_parameters)'

      ! elastic time step
      !dte = dt/real(ndte,kind=dbl_kind)        ! s
      !dtei = c1/dte              ! 1/s
      dtei = real(ndte,kind=dbl_kind)/dt

      ! variables for elliptical yield curve and plastic potential
      epp2i = c1/e_plasticpot**2
      e_factor = e_yieldcurve**2 / e_plasticpot**4
      ecci = c1/e_yieldcurve**2 ! temporary for 1d evp

      ! constants for stress equation
      !tdamp2 = c2 * elasticDamp * dt                 ! s
      !dte2T = dte/tdamp2    or c1/(c2*elasticDamp*real(ndte,kind=dbl_kind))               ! ellipse (unitless)

      if (revised_evp) then       ! Bouillon et al, Ocean Mod 2013
         revp   = c1
         denom1 = c1
         arlx1i = c1/arlx
      else                        ! Hunke, JCP 2013 with modified stress eq
         revp   = c0
         !arlx1i = dte2T
         !arlx   = c1/arlx1i
         !brlx   = dt*dtei
         arlx   = c2 * elasticDamp * real(ndte,kind=dbl_kind)
         arlx1i   = c1/arlx
         brlx   = real(ndte,kind=dbl_kind)
         denom1 = c1/(c1+arlx1i)
      endif
      if (my_task == master_task) then
         write (nu_diag,*) 'arlx, arlxi, brlx, denom1', &
                  arlx, arlx1i, brlx, denom1
      endif

      end subroutine set_evp_parameters

!=======================================================================
! Computes quantities needed in the stress tensor (sigma)
! and momentum (u) equations, but which do not change during
! the thermodynamics/transport time step:
! ice mass and ice extent masks
!
! author: Elizabeth C. Hunke, LANL

      subroutine dyn_prep1 (nx_block,  ny_block, &
                            ilo, ihi,  jlo, jhi, &
                            aice,      vice,     &
                            vsno,      tmask,    &
                            tmass,     icetmask)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ilo,ihi,jlo,jhi       ! beginning and end of physical domain

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         aice    , & ! concentration of ice
         vice    , & ! volume per unit area of ice          (m)
         vsno        ! volume per unit area of snow         (m)

      logical (kind=log_kind), dimension (nx_block,ny_block), intent(in) :: &
         tmask       ! land/boundary mask, thickness (T-cell)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out) :: &
         tmass       ! total mass of ice and snow (kg/m^2)

      integer (kind=int_kind), dimension (nx_block,ny_block), intent(out) :: &
         icetmask    ! ice extent mask (T-cell)

      ! local variables

      integer (kind=int_kind) :: &
         i, j

      real (kind=dbl_kind) :: &
         rhoi, rhos

      logical (kind=log_kind), dimension(nx_block,ny_block) :: &
         tmphm               ! temporary mask

      character(len=*), parameter :: subname = '(dyn_prep1)'

      call icepack_query_parameters(rhos_out=rhos, rhoi_out=rhoi)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      do j = 1, ny_block
      do i = 1, nx_block

         !-----------------------------------------------------------------
         ! total mass of ice and snow, centered in T-cell
         ! NOTE: vice and vsno must be up to date in all grid cells,
         !       including ghost cells
         !-----------------------------------------------------------------
         if (tmask(i,j)) then
            tmass(i,j) = (rhoi*vice(i,j) + rhos*vsno(i,j)) ! kg/m^2
         else
            tmass(i,j) = c0
         endif

         !-----------------------------------------------------------------
         ! ice extent mask (T-cells)
         !-----------------------------------------------------------------
         tmphm(i,j) = tmask(i,j) .and. (aice (i,j) > a_min) &
                                 .and. (tmass(i,j) > m_min)

         !-----------------------------------------------------------------
         ! augmented mask (land + open ocean)
         !-----------------------------------------------------------------
         icetmask (i,j) = 0

      enddo
      enddo

      do j = jlo, jhi
      do i = ilo, ihi

         ! extend ice extent mask (T-cells) to points around pack
         if (tmphm(i-1,j+1) .or. tmphm(i,j+1) .or. tmphm(i+1,j+1) .or. &
             tmphm(i-1,j)   .or. tmphm(i,j)   .or. tmphm(i+1,j)   .or. &
             tmphm(i-1,j-1) .or. tmphm(i,j-1) .or. tmphm(i+1,j-1) ) then
            icetmask(i,j) = 1
         endif

         if (.not.tmask(i,j)) icetmask(i,j) = 0

      enddo
      enddo

      end subroutine dyn_prep1

!=======================================================================
! Computes quantities needed in the stress tensor (sigma)
! and momentum (u) equations, but which do not change during
! the thermodynamics/transport time step:
! --wind stress shift to U grid,
! --ice mass and ice extent masks,
! initializes ice velocity for new points to ocean sfc current
!
! author: Elizabeth C. Hunke, LANL

      subroutine dyn_prep2 (nx_block,   ny_block,   &
                            ilo, ihi,   jlo, jhi,   &
                            icellt,     icellu,     &
                            indxti,     indxtj,     &
                            indxui,     indxuj,     &
                            aiu,        umass,      &
                            umassdti,   fcor,       &
                            umask,                  &
                            uocn,       vocn,       &
                            strairx,    strairy,    &
                            ss_tltx,    ss_tlty,    &
                            icetmask,   iceumask,   &
                            fm,         dt,         &
                            strtltx,    strtlty,    &
                            strocnx,    strocny,    &
                            strintx,    strinty,    &
                            taubx,      tauby,      &
                            waterx,     watery,     &
                            forcex,     forcey,     &
                            stressp_1,  stressp_2,  &
                            stressp_3,  stressp_4,  &
                            stressm_1,  stressm_2,  &
                            stressm_3,  stressm_4,  &
                            stress12_1, stress12_2, &
                            stress12_3, stress12_4, &
                            uvel_init,  vvel_init,  &
                            uvel,       vvel,       &
                            Tbu)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ilo,ihi,jlo,jhi       ! beginning and end of physical domain

      integer (kind=int_kind), intent(out) :: &
         icellt  , & ! no. of cells where icetmask = 1
         icellu      ! no. of cells where iceumask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(out) :: &
         indxti  , & ! compressed index in i-direction
         indxtj  , & ! compressed index in j-direction
         indxui  , & ! compressed index in i-direction
         indxuj      ! compressed index in j-direction

      logical (kind=log_kind), dimension (nx_block,ny_block), intent(in) :: &
         umask       ! land/boundary mask, thickness (U-cell)

      integer (kind=int_kind), dimension (nx_block,ny_block), intent(in) :: &
         icetmask    ! ice extent mask (T-cell)

      logical (kind=log_kind), dimension (nx_block,ny_block), intent(inout) :: &
         iceumask    ! ice extent mask (U-cell)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         aiu     , & ! ice fraction on u-grid
         umass   , & ! total mass of ice and snow (u grid)
         fcor    , & ! Coriolis parameter (1/s)
         strairx , & ! stress on ice by air, x-direction
         strairy , & ! stress on ice by air, y-direction
         uocn    , & ! ocean current, x-direction (m/s)
         vocn    , & ! ocean current, y-direction (m/s)
         ss_tltx , & ! sea surface slope, x-direction (m/m)
         ss_tlty     ! sea surface slope, y-direction

      real (kind=dbl_kind), intent(in) :: &
         dt          ! time step

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out) :: &
         Tbu,      & ! seabed stress factor (N/m^2)
         uvel_init,& ! x-component of velocity (m/s), beginning of time step
         vvel_init,& ! y-component of velocity (m/s), beginning of time step
         umassdti, & ! mass of U-cell/dt (kg/m^2 s)
         waterx  , & ! for ocean stress calculation, x (m/s)
         watery  , & ! for ocean stress calculation, y (m/s)
         forcex  , & ! work array: combined atm stress and ocn tilt, x
         forcey      ! work array: combined atm stress and ocn tilt, y

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         fm      , & ! Coriolis param. * mass in U-cell (kg/s)
         stressp_1, stressp_2, stressp_3, stressp_4 , & ! sigma11+sigma22
         stressm_1, stressm_2, stressm_3, stressm_4 , & ! sigma11-sigma22
         stress12_1,stress12_2,stress12_3,stress12_4, & ! sigma12
         uvel    , & ! x-component of velocity (m/s)
         vvel    , & ! y-component of velocity (m/s)
         strtltx , & ! stress due to sea surface slope, x-direction
         strtlty , & ! stress due to sea surface slope, y-direction
         strocnx , & ! ice-ocean stress, x-direction
         strocny , & ! ice-ocean stress, y-direction
         strintx , & ! divergence of internal ice stress, x (N/m^2)
         strinty , & ! divergence of internal ice stress, y (N/m^2)
         taubx   , & ! seabed stress, x-direction (N/m^2)
         tauby       ! seabed stress, y-direction (N/m^2)

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &
         gravit

      logical (kind=log_kind), dimension(nx_block,ny_block) :: &
         iceumask_old      ! old-time iceumask

      character(len=*), parameter :: subname = '(dyn_prep2)'

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      do j = 1, ny_block
      do i = 1, nx_block
         waterx   (i,j) = c0
         watery   (i,j) = c0
         forcex   (i,j) = c0
         forcey   (i,j) = c0
         umassdti (i,j) = c0
         Tbu      (i,j) = c0
         taubx    (i,j) = c0
         tauby    (i,j) = c0

         if (icetmask(i,j)==0) then
            stressp_1 (i,j) = c0
            stressp_2 (i,j) = c0
            stressp_3 (i,j) = c0
            stressp_4 (i,j) = c0
            stressm_1 (i,j) = c0
            stressm_2 (i,j) = c0
            stressm_3 (i,j) = c0
            stressm_4 (i,j) = c0
            stress12_1(i,j) = c0
            stress12_2(i,j) = c0
            stress12_3(i,j) = c0
            stress12_4(i,j) = c0
         endif
      enddo                     ! i
      enddo                     ! j

      !-----------------------------------------------------------------
      ! Identify cells where icetmask = 1
      ! Note: The icellt mask includes north and east ghost cells
      !       where stresses are needed.
      !-----------------------------------------------------------------

      icellt = 0
      do j = jlo, jhi+1
      do i = ilo, ihi+1
         if (icetmask(i,j) == 1) then
            icellt = icellt + 1
            indxti(icellt) = i
            indxtj(icellt) = j
         endif
      enddo
      enddo

      !-----------------------------------------------------------------
      ! Define iceumask
      ! Identify cells where iceumask is true
      ! Initialize velocity where needed
      !-----------------------------------------------------------------

      icellu = 0

      do j = jlo, jhi
      do i = ilo, ihi
         iceumask_old(i,j) = iceumask(i,j) ! save
!         if (grid_ice == 'B') then ! include ice mask.
         ! ice extent mask (U-cells)
         iceumask(i,j) = (umask(i,j)) .and. (aiu  (i,j) > a_min) &
                                         .and. (umass(i,j) > m_min)
!         else  ! ice mask shpuld be applied to cd grid. For now it is not implemented.
!            iceumask(i,j) = umask(i,j)
!         endif

         if (iceumask(i,j)) then
            icellu = icellu + 1
            indxui(icellu) = i
            indxuj(icellu) = j

            ! initialize velocity for new ice points to ocean sfc current
            if (.not. iceumask_old(i,j)) then
               uvel(i,j) = uocn(i,j)
               vvel(i,j) = vocn(i,j)
            endif
         else
            ! set velocity and stresses to zero for masked-out points
            uvel(i,j)    = c0
            vvel(i,j)    = c0
            strintx(i,j) = c0
            strinty(i,j) = c0
            strocnx(i,j) = c0
            strocny(i,j) = c0
         endif

         uvel_init(i,j) = uvel(i,j)
         vvel_init(i,j) = vvel(i,j)
      enddo
      enddo

      !-----------------------------------------------------------------
      ! Define variables for momentum equation
      !-----------------------------------------------------------------

      if (trim(ssh_stress) == 'coupled') then
         call icepack_query_parameters(gravit_out=gravit)
         call icepack_warnings_flush(nu_diag)
         if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
            file=__FILE__, line=__LINE__)
      endif

      do ij = 1, icellu
         i = indxui(ij)
         j = indxuj(ij)

         umassdti(i,j) = umass(i,j)/dt ! kg/m^2 s

         fm(i,j) = fcor(i,j)*umass(i,j)   ! Coriolis * mass

         ! for ocean stress
         waterx(i,j) = uocn(i,j)*cosw - vocn(i,j)*sinw*sign(c1,fm(i,j))
         watery(i,j) = vocn(i,j)*cosw + uocn(i,j)*sinw*sign(c1,fm(i,j))

         ! combine tilt with wind stress
         if (trim(ssh_stress) == 'geostrophic') then
            ! calculate tilt from geostrophic currents if needed
            strtltx(i,j) = -fm(i,j)*vocn(i,j)
            strtlty(i,j) =  fm(i,j)*uocn(i,j)
         elseif (trim(ssh_stress) == 'coupled') then
            strtltx(i,j) = -gravit*umass(i,j)*ss_tltx(i,j)
            strtlty(i,j) = -gravit*umass(i,j)*ss_tlty(i,j)
         else
            call abort_ice(subname//' ERROR: unknown ssh_stress='//trim(ssh_stress), &
               file=__FILE__, line=__LINE__)
         endif

         forcex(i,j) = strairx(i,j) + strtltx(i,j)
         forcey(i,j) = strairy(i,j) + strtlty(i,j)
      enddo

      end subroutine dyn_prep2

!=======================================================================
! Calculation of the surface stresses
! Integration of the momentum equation to find velocity (u,v)
!
! author: Elizabeth C. Hunke, LANL

      subroutine stepu (nx_block,   ny_block, &
                        icellu,     Cw,       &
                        indxui,     indxuj,   &
                        aiu,        str,      &
                        uocn,       vocn,     &
                        waterx,     watery,   &
                        forcex,     forcey,   &
                        umassdti,   fm,       &
                        uarear,               &
                        strintx,    strinty,  &
                        taubx,      tauby,    &
                        uvel_init,  vvel_init,&
                        uvel,       vvel,     &
                        Tbu)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellu                ! total count when iceumask is true

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxui  , & ! compressed index in i-direction
         indxuj      ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         Tbu,      & ! seabed stress factor (N/m^2)
         uvel_init,& ! x-component of velocity (m/s), beginning of timestep
         vvel_init,& ! y-component of velocity (m/s), beginning of timestep
         aiu     , & ! ice fraction on u-grid
         waterx  , & ! for ocean stress calculation, x (m/s)
         watery  , & ! for ocean stress calculation, y (m/s)
         forcex  , & ! work array: combined atm stress and ocn tilt, x
         forcey  , & ! work array: combined atm stress and ocn tilt, y
         umassdti, & ! mass of U-cell/dt (kg/m^2 s)
         uocn    , & ! ocean current, x-direction (m/s)
         vocn    , & ! ocean current, y-direction (m/s)
         fm      , & ! Coriolis param. * mass in U-cell (kg/s)
         uarear      ! 1/uarea

      real (kind=dbl_kind), dimension(nx_block,ny_block,8), intent(in) :: &
         str         ! temporary

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         uvel    , & ! x-component of velocity (m/s)
         vvel        ! y-component of velocity (m/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         strintx , & ! divergence of internal ice stress, x (N/m^2)
         strinty , & ! divergence of internal ice stress, y (N/m^2)
         taubx   , & ! seabed stress, x-direction (N/m^2)
         tauby       ! seabed stress, y-direction (N/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         Cw          ! ocean-ice neutral drag coefficient

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &
         uold, vold        , & ! old-time uvel, vvel
         vrel              , & ! relative ice-ocean velocity
         cca,ccb,ab2,cc1,cc2,& ! intermediate variables
         taux, tauy        , & ! part of ocean stress term
         Cb                , & ! complete seabed (basal) stress coeff
         rhow                  !

      character(len=*), parameter :: subname = '(stepu)'

      !-----------------------------------------------------------------
      ! integrate the momentum equation
      !-----------------------------------------------------------------

      call icepack_query_parameters(rhow_out=rhow)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      do ij =1, icellu
         i = indxui(ij)
         j = indxuj(ij)

         uold = uvel(i,j)
         vold = vvel(i,j)

         ! (magnitude of relative ocean current)*rhow*drag*aice
         vrel = aiu(i,j)*rhow*Cw(i,j)*sqrt((uocn(i,j) - uold)**2 + &
                                           (vocn(i,j) - vold)**2)  ! m/s
         ! ice/ocean stress
         taux = vrel*waterx(i,j) ! NOTE this is not the entire
         tauy = vrel*watery(i,j) ! ocn stress term

         Cb  = Tbu(i,j) / (sqrt(uold**2 + vold**2) + u0) ! for seabed stress
         ! revp = 0 for classic evp, 1 for revised evp
         cca = (brlx + revp)*umassdti(i,j) + vrel * cosw + Cb ! kg/m^2 s

         ccb = fm(i,j) + sign(c1,fm(i,j)) * vrel * sinw ! kg/m^2 s

         ab2 = cca**2 + ccb**2

         ! divergence of the internal stress tensor
         strintx(i,j) = uarear(i,j)* &
             (str(i,j,1) + str(i+1,j,2) + str(i,j+1,3) + str(i+1,j+1,4))
         strinty(i,j) = uarear(i,j)* &
             (str(i,j,5) + str(i,j+1,6) + str(i+1,j,7) + str(i+1,j+1,8))

         ! finally, the velocity components
         cc1 = strintx(i,j) + forcex(i,j) + taux &
             + umassdti(i,j)*(brlx*uold + revp*uvel_init(i,j))
         cc2 = strinty(i,j) + forcey(i,j) + tauy &
             + umassdti(i,j)*(brlx*vold + revp*vvel_init(i,j))

         uvel(i,j) = (cca*cc1 + ccb*cc2) / ab2 ! m/s
         vvel(i,j) = (cca*cc2 - ccb*cc1) / ab2

         ! calculate seabed stress component for outputs
         ! only needed on last iteration.
         taubx(i,j) = -uvel(i,j)*Cb
         tauby(i,j) = -vvel(i,j)*Cb
      enddo                     ! ij

      end subroutine stepu

!=======================================================================
! Integration of the momentum equation to find velocity (u,v) at E and N locations

      subroutine stepuv_CD (nx_block,   ny_block, &
                          icell,      Cw,       &
                          indxi,      indxj,    &
                                      aiu,      &
                          uocn,       vocn,     &
                          waterx,     watery,   &
                          forcex,     forcey,   &
                          massdti,    fm,       &
                          strintx,    strinty,  &
                          taubx,      tauby,    &
                          uvel_init,  vvel_init,&
                          uvel,       vvel,     &
                          Tb)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icell                 ! total count when ice[en]mask is true

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxi   , & ! compressed index in i-direction
         indxj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         Tb,       & ! seabed stress factor (N/m^2)
         uvel_init,& ! x-component of velocity (m/s), beginning of timestep
         vvel_init,& ! y-component of velocity (m/s), beginning of timestep
         aiu     , & ! ice fraction on [en]-grid
         waterx  , & ! for ocean stress calculation, x (m/s)
         watery  , & ! for ocean stress calculation, y (m/s)
         forcex  , & ! work array: combined atm stress and ocn tilt, x
         forcey  , & ! work array: combined atm stress and ocn tilt, y
         massdti , & ! mass of [EN]-cell/dt (kg/m^2 s)
         uocn    , & ! ocean current, x-direction (m/s)
         vocn    , & ! ocean current, y-direction (m/s)
         fm      , & ! Coriolis param. * mass in [EN]-cell (kg/s)
         strintx , & ! divergence of internal ice stress, x (N/m^2)
         strinty     ! divergence of internal ice stress, y (N/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         uvel    , & ! x-component of velocity (m/s)
         vvel        ! y-component of velocity (m/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         taubx   , & ! seabed stress, x-direction (N/m^2)
         tauby       ! seabed stress, y-direction (N/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         Cw          ! ocean-ice neutral drag coefficient

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &
         uold, vold         , & ! old-time uvel, vvel
         vrel               , & ! relative ice-ocean velocity
         cca,ccb,ccc,ab2    , & ! intermediate variables
         cc1,cc2            , & ! "
         taux, tauy         , & ! part of ocean stress term
         Cb                 , & ! complete seabed (basal) stress coeff
         rhow                   !

      character(len=*), parameter :: subname = '(stepuv_CD)'

      !-----------------------------------------------------------------
      ! integrate the momentum equation
      !-----------------------------------------------------------------

      call icepack_query_parameters(rhow_out=rhow)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      do ij =1, icell
         i = indxi(ij)
         j = indxj(ij)

         uold = uvel(i,j)
         vold = vvel(i,j)

         ! (magnitude of relative ocean current)*rhow*drag*aice
         vrel = aiu(i,j)*rhow*Cw(i,j)*sqrt((uocn(i,j) - uold)**2 + &
                                           (vocn(i,j) - vold)**2)  ! m/s
         ! ice/ocean stress
         taux = vrel*waterx(i,j) ! NOTE this is not the entire
         tauy = vrel*watery(i,j) ! ocn stress term

         ccc = sqrt(uold**2 + vold**2) + u0
         Cb  = Tb(i,j) / ccc ! for seabed stress
         ! revp = 0 for classic evp, 1 for revised evp
         cca = (brlx + revp)*massdti(i,j) + vrel * cosw + Cb ! kg/m^2 s

         ccb = fm(i,j) + sign(c1,fm(i,j)) * vrel * sinw ! kg/m^2 s

         ab2 = cca**2 + ccb**2

         ! compute the velocity components
         cc1 = strintx(i,j) + forcex(i,j) + taux &
             + massdti(i,j)*(brlx*uold + revp*uvel_init(i,j))
         cc2 = strinty(i,j) + forcey(i,j) + tauy &
             + massdti(i,j)*(brlx*vold + revp*vvel_init(i,j))
         uvel(i,j) = (cca*cc1 + ccb*cc2) / ab2 ! m/s
         vvel(i,j) = (cca*cc2 - ccb*cc1) / ab2

         ! calculate seabed stress component for outputs
         ! only needed on last iteration.
         taubx(i,j) = -uvel(i,j)*Cb
         tauby(i,j) = -vvel(i,j)*Cb

      enddo                     ! ij

      end subroutine stepuv_CD

!=======================================================================
! Integration of the momentum equation to find velocity u at E location on C grid

      subroutine stepu_C (nx_block,   ny_block, &
                           icell,      Cw,       &
                           indxi,      indxj,    &
                                       aiu,      &
                           uocn,       vocn,     &
                           waterx,     forcex,   &
                           massdti,    fm,       &
                           strintx,    taubx,    &
                           uvel_init,            &
                           uvel,       vvel,     &
                           Tb)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icell                 ! total count when ice[en]mask is true

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxi   , & ! compressed index in i-direction
         indxj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         Tb,       & ! seabed stress factor (N/m^2)
         uvel_init,& ! x-component of velocity (m/s), beginning of timestep
         aiu     , & ! ice fraction on [en]-grid
         waterx  , & ! for ocean stress calculation, x (m/s)
         forcex  , & ! work array: combined atm stress and ocn tilt, x
         massdti , & ! mass of e-cell/dt (kg/m^2 s)
         uocn    , & ! ocean current, x-direction (m/s)
         vocn    , & ! ocean current, y-direction (m/s)
         fm      , & ! Coriolis param. * mass in e-cell (kg/s)
         strintx , & ! divergence of internal ice stress, x (N/m^2)
         Cw      , & ! ocean-ice neutral drag coefficient
         vvel        ! y-component of velocity (m/s) interpolated to E location

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         uvel    , & ! x-component of velocity (m/s)
         taubx       ! seabed stress, x-direction (N/m^2)

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &
         uold, vold         , & ! old-time uvel, vvel
         vrel               , & ! relative ice-ocean velocity
         cca,ccb,ccc,cc1    , & ! intermediate variables
         taux               , & ! part of ocean stress term
         Cb                 , & ! complete seabed (basal) stress coeff
         rhow                   !

      character(len=*), parameter :: subname = '(stepu_C)'

      !-----------------------------------------------------------------
      ! integrate the momentum equation
      !-----------------------------------------------------------------

      call icepack_query_parameters(rhow_out=rhow)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      do ij =1, icell
         i = indxi(ij)
         j = indxj(ij)

         uold = uvel(i,j)
         vold = vvel(i,j)

         ! (magnitude of relative ocean current)*rhow*drag*aice
         vrel = aiu(i,j)*rhow*Cw(i,j)*sqrt((uocn(i,j) - uold)**2 + &
                                           (vocn(i,j) - vold)**2)  ! m/s
         ! ice/ocean stress
         taux = vrel*waterx(i,j) ! NOTE this is not the entire

         ccc = sqrt(uold**2 + vold**2) + u0
         Cb  = Tb(i,j) / ccc ! for seabed stress
         ! revp = 0 for classic evp, 1 for revised evp
         cca = (brlx + revp)*massdti(i,j) + vrel * cosw + Cb ! kg/m^2 s

         ccb = fm(i,j) + sign(c1,fm(i,j)) * vrel * sinw ! kg/m^2 s

         ! compute the velocity components
         cc1 = strintx(i,j) + forcex(i,j) + taux &
             + massdti(i,j)*(brlx*uold + revp*uvel_init(i,j))

         uvel(i,j) = (ccb*vold + cc1) / cca ! m/s

         ! calculate seabed stress component for outputs
         ! only needed on last iteration.
         taubx(i,j) = -uvel(i,j)*Cb

      enddo                     ! ij

      end subroutine stepu_C

!=======================================================================
! Integration of the momentum equation to find velocity v at N location on C grid

      subroutine stepv_C (nx_block,   ny_block, &
                           icell,      Cw,       &
                           indxi,      indxj,    &
                                       aiu,      &
                           uocn,       vocn,     &
                           watery,     forcey,   &
                           massdti,    fm,       &
                           strinty,    tauby,    &
                           vvel_init,            &
                           uvel,       vvel,     &
                           Tb)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icell                 ! total count when ice[en]mask is true

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxi   , & ! compressed index in i-direction
         indxj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         Tb,       & ! seabed stress factor (N/m^2)
         vvel_init,& ! y-component of velocity (m/s), beginning of timestep
         aiu     , & ! ice fraction on [en]-grid
         watery  , & ! for ocean stress calculation, y (m/s)
         forcey  , & ! work array: combined atm stress and ocn tilt, y
         massdti , & ! mass of n-cell/dt (kg/m^2 s)
         uocn    , & ! ocean current, x-direction (m/s)
         vocn    , & ! ocean current, y-direction (m/s)
         fm      , & ! Coriolis param. * mass in n-cell (kg/s)
         strinty , & ! divergence of internal ice stress, y (N/m^2)
         Cw      , & ! ocean-ice neutral drag coefficient
         uvel        ! x-component of velocity (m/s) interpolated to N location

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         vvel    , & ! y-component of velocity (m/s)
         tauby       ! seabed stress, y-direction (N/m^2)

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &
         uold, vold         , & ! old-time uvel, vvel
         vrel               , & ! relative ice-ocean velocity
         cca,ccb,ccc,cc2    , & ! intermediate variables
         tauy               , & ! part of ocean stress term
         Cb                 , & ! complete seabed (basal) stress coeff
         rhow                   !

      character(len=*), parameter :: subname = '(stepv_C)'

      !-----------------------------------------------------------------
      ! integrate the momentum equation
      !-----------------------------------------------------------------

      call icepack_query_parameters(rhow_out=rhow)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      do ij =1, icell
         i = indxi(ij)
         j = indxj(ij)

         uold = uvel(i,j)
         vold = vvel(i,j)

         ! (magnitude of relative ocean current)*rhow*drag*aice
         vrel = aiu(i,j)*rhow*Cw(i,j)*sqrt((uocn(i,j) - uold)**2 + &
                                           (vocn(i,j) - vold)**2)  ! m/s
         ! ice/ocean stress
         tauy = vrel*watery(i,j) ! NOTE this is not the entire ocn stress

         ccc = sqrt(uold**2 + vold**2) + u0
         Cb  = Tb(i,j) / ccc ! for seabed stress
         ! revp = 0 for classic evp, 1 for revised evp
         cca = (brlx + revp)*massdti(i,j) + vrel * cosw + Cb ! kg/m^2 s

         ccb = fm(i,j) + sign(c1,fm(i,j)) * vrel * sinw ! kg/m^2 s

         ! compute the velocity components
         cc2 = strinty(i,j) + forcey(i,j) + tauy &
             + massdti(i,j)*(brlx*vold + revp*vvel_init(i,j))

         vvel(i,j) = (-ccb*uold + cc2) / cca

         ! calculate seabed stress component for outputs
         ! only needed on last iteration.
         tauby(i,j) = -vvel(i,j)*Cb

      enddo                     ! ij

      end subroutine stepv_C

!=======================================================================
! Calculation of the ice-ocean stress.
! ...the sign will be reversed later...
!
! author: Elizabeth C. Hunke, LANL

      subroutine dyn_finish (nx_block, ny_block, &
                             icellu,   Cw,       &
                             indxui,   indxuj,   &
                             uvel,     vvel,     &
                             uocn,     vocn,     &
                             aiu,      fm,       &
                             strocnx,  strocny)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellu                ! total count when iceumask is true

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxui  , & ! compressed index in i-direction
         indxuj      ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         uvel    , & ! x-component of velocity (m/s)
         vvel    , & ! y-component of velocity (m/s)
         uocn    , & ! ocean current, x-direction (m/s)
         vocn    , & ! ocean current, y-direction (m/s)
         aiu     , & ! ice fraction on u-grid
         fm          ! Coriolis param. * mass in U-cell (kg/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         strocnx , & ! ice-ocean stress, x-direction
         strocny     ! ice-ocean stress, y-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         Cw          ! ocean-ice neutral drag coefficient

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &
         vrel    , & !
         rhow        !

      character(len=*), parameter :: subname = '(dyn_finish)'

      call icepack_query_parameters(rhow_out=rhow)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      ! ocean-ice stress for coupling
      do ij =1, icellu
         i = indxui(ij)
         j = indxuj(ij)

          vrel = rhow*Cw(i,j)*sqrt((uocn(i,j) - uvel(i,j))**2 + &
                 (vocn(i,j) - vvel(i,j))**2)  ! m/s

!        strocnx(i,j) = strocnx(i,j) &
!                     - vrel*(uvel(i,j)*cosw - vvel(i,j)*sinw) * aiu(i,j)
!        strocny(i,j) = strocny(i,j) &
!                     - vrel*(vvel(i,j)*cosw + uvel(i,j)*sinw) * aiu(i,j)

         ! update strocnx to most recent iterate and complete the term
         vrel = vrel * aiu(i,j)
         strocnx(i,j) = vrel*((uocn(i,j) - uvel(i,j))*cosw &
                            - (vocn(i,j) - vvel(i,j))*sinw*sign(c1,fm(i,j)))
         strocny(i,j) = vrel*((vocn(i,j) - vvel(i,j))*cosw &
                            + (uocn(i,j) - uvel(i,j))*sinw*sign(c1,fm(i,j)))

         ! Hibler/Bryan stress
         ! the sign is reversed later, therefore negative here
!         strocnx(i,j) = -(strairx(i,j) + strintx(i,j))
!         strocny(i,j) = -(strairy(i,j) + strinty(i,j))

      enddo

      end subroutine dyn_finish

!=======================================================================
! Computes seabed (basal) stress factor Tbu (landfast ice) based on mean
! thickness and bathymetry data. LKD refers to linear keel draft. This
! parameterization assumes that the largest keel draft varies linearly
! with the mean thickness.
!
! Lemieux, J. F., B. Tremblay, F. Dupont, M. Plante, G.C. Smith, D. Dumont (2015).
! A basal stress parameterization form modeling landfast ice, J. Geophys. Res.
! Oceans, 120, 3157-3173.
!
! Lemieux, J. F., F. Dupont, P. Blain, F. Roy, G.C. Smith, G.M. Flato (2016).
! Improving the simulation of landfast ice by combining tensile strength and a
! parameterization for grounded ridges, J. Geophys. Res. Oceans, 121, 7354-7368.
!
! author: JF Lemieux, Philippe Blain (ECCC)
!
! note1: Tbu is a part of the Cb as defined in Lemieux et al. 2015 and 2016.
! note2: Seabed stress (better name) was called basal stress in Lemieux et al. 2015

      subroutine seabed_stress_factor_LKD (nx_block, ny_block,         &
                                           icellu,                     &
                                           indxui,   indxuj,           &
                                           vice,     aice,             &
                                           hwater,   Tbu,              &
                                           grid_location)

      use ice_grid, only: grid_neighbor_min, grid_neighbor_max

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, &  ! block dimensions
         icellu                 ! no. of cells where ice[uen]mask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxui    , & ! compressed index in i-direction
         indxuj        ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         aice      , & ! concentration of ice at tracer location
         vice      , & ! volume per unit area of ice at tracer location (m)
         hwater        ! water depth at tracer location (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         Tbu           ! seabed stress factor at 'grid_location' (N/m^2)

      character(len=*), optional, intent(inout) :: &
         grid_location ! grid location (U, E, N), U assumed if not present

      real (kind=dbl_kind) :: &
         au        , & ! concentration of ice at u location
         hu        , & ! volume per unit area of ice at u location (mean thickness, m)
         hwu       , & ! water depth at u location (m)
         docalc_tbu, & ! logical as real (C0,C1) decides whether c0 is 0 or
         hcu           ! critical thickness at u location (m)

      integer (kind=int_kind) :: &
         i, j, ij

      character(len=char_len) :: &
         l_grid_location ! local version of 'grid_location'

      character(len=*), parameter :: subname = '(seabed_stress_factor_LKD)'

      ! Assume U location (NE corner) if grid_location not present
      if (.not. (present(grid_location))) then
         l_grid_location = 'U'
      else
         l_grid_location = grid_location
      endif

      do ij = 1, icellu
         i = indxui(ij)
         j = indxuj(ij)

         ! convert quantities to grid_location

         hwu = grid_neighbor_min(hwater, i, j, l_grid_location)

         docalc_tbu = merge(c1,c0,hwu < threshold_hw)


         au  = grid_neighbor_max(aice, i, j, l_grid_location)
         hu  = grid_neighbor_max(vice, i, j, l_grid_location)

         ! 1- calculate critical thickness
         hcu = au * hwu / k1

         ! 2- calculate seabed stress factor
         Tbu(i,j) = docalc_tbu*k2 * max(c0,(hu - hcu)) * exp(-alphab * (c1 - au))

      enddo                     ! ij

      end subroutine seabed_stress_factor_LKD

!=======================================================================
! Computes seabed (basal) stress factor Tbu (landfast ice) based on
! probability of contact between the ITD and the seabed. The water depth
! could take into account variations of the SSH. In the simplest
! formulation, hwater is simply the value of the bathymetry. To calculate
! the probability of contact, it is assumed that the bathymetry follows
! a normal distribution with sigma_b = 2.5d0. An improvement would
! be to provide the distribution based on high resolution data.
!
! Dupont, F. Dumont, D., Lemieux, J.F., Dumas-Lefebvre, E., Caya, A.
! in prep.
!
! authors: D. Dumont, J.F. Lemieux, E. Dumas-Lefebvre, F. Dupont
!
      subroutine seabed_stress_factor_prob (nx_block, ny_block,          &
                                            icellt, indxti,   indxtj,    &
                                            icellu, indxui,   indxuj,    &
                                            aicen,  vicen,               &
                                            hwater, Tbu,                 &
                                            TbE,    TbN,                 &
                                            icelle, indxei,   indxej,    &
                                            icelln, indxni,   indxnj)
! use modules

      use ice_arrays_column, only: hin_max
      use ice_domain_size, only: ncat
      use ice_grid, only: grid_neighbor_min, grid_neighbor_max

      integer (kind=int_kind), intent(in) :: &
           nx_block, ny_block, &  ! block dimensions
           icellt, icellu         ! no. of cells where ice[tu]mask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
           indxti  , & ! compressed index in i-direction
           indxtj  , & ! compressed index in j-direction
           indxui  , & ! compressed index in i-direction
           indxuj      ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
           hwater      ! water depth at tracer location (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), intent(in) :: &
           aicen,    & ! partial concentration for last thickness category in ITD
           vicen       ! partial volume for last thickness category in ITD (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
           Tbu         ! seabed stress factor at U location (N/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout), optional :: &
           TbE,      & ! seabed stress factor at E location (N/m^2)
           TbN         ! seabed stress factor at N location (N/m^2)

      integer (kind=int_kind), intent(in), optional :: &
           icelle, icelln ! no. of cells where ice[en]mask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in), optional :: &
           indxei  , & ! compressed index in i-direction
           indxej  , & ! compressed index in j-direction
           indxni  , & ! compressed index in i-direction
           indxnj      ! compressed index in j-direction

! local variables

      integer (kind=int_kind) :: &
           i, j, ij, ii, n

      integer, parameter :: &
           ncat_b = 100, &  ! number of bathymetry categories
           ncat_i = 100     ! number of ice thickness categories (log-normal)

      real (kind=dbl_kind), parameter :: &
           max_depth = 50.0_dbl_kind, & ! initial range of log-normal distribution
           mu_s = 0.1_dbl_kind, &       ! friction coefficient
           sigma_b = 2.5d0              ! Standard deviation of bathymetry

      real (kind=dbl_kind), dimension(ncat_i) :: & ! log-normal for ice thickness
           x_k, & ! center of thickness categories (m)
           g_k, & ! probability density function (thickness, 1/m)
           P_x    ! probability for each thickness category

      real (kind=dbl_kind), dimension(ncat_b) :: & ! normal dist for bathymetry
           y_n, & ! center of bathymetry categories (m)
           b_n, & ! probability density function (bathymetry, 1/m)
           P_y    ! probability for each bathymetry category

      real (kind=dbl_kind), dimension(ncat) :: &
           vcat, acat  ! vice, aice temporary arrays

      integer, dimension(ncat_b) :: &
           tmp    ! Temporary vector tmp = merge(1,0,gt)

      logical, dimension (ncat_b) :: &
           gt     !

      real (kind=dbl_kind) :: &
           wid_i, wid_b  , & ! parameters for PDFs
           mu_i, sigma_i , & !
           mu_b, m_i, v_i, & !
           atot, x_kmax  , & !
           cut           , & !
           rhoi, rhow    , & !
           gravit        , & !
           pi, puny          !

      real (kind=dbl_kind), dimension(ncat_i) :: &
           tb_tmp

      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
           Tbt    ! seabed stress factor at t point (N/m^2)

      character(len=*), parameter :: subname = '(seabed_stress_factor_prob)'

      call icepack_query_parameters(rhow_out=rhow, rhoi_out=rhoi)
      call icepack_query_parameters(gravit_out=gravit)
      call icepack_query_parameters(pi_out=pi)
      call icepack_query_parameters(puny_out=puny)

      Tbt=c0

      do ij = 1, icellt
         i = indxti(ij)
         j = indxtj(ij)

         atot = sum(aicen(i,j,1:ncat))

         if (atot > 0.05_dbl_kind .and. hwater(i,j) < max_depth) then

            mu_b = hwater(i,j)           ! mean of PDF (normal dist) bathymetry
            wid_i = max_depth/ncat_i     ! width of ice categories
            wid_b = c6*sigma_b/ncat_b    ! width of bathymetry categories (6 sigma_b = 2x3 sigma_b)

            x_k = (/( wid_i*( real(i,kind=dbl_kind) - p5 ), i=1, ncat_i )/)
            y_n = (/( ( mu_b-c3*sigma_b )+( real(i,kind=dbl_kind) - p5 )*( c6*sigma_b/ncat_b ), i=1, ncat_b )/)

            vcat(1:ncat) = vicen(i,j,1:ncat)
            acat(1:ncat) = aicen(i,j,1:ncat)

            m_i = sum(vcat)

            v_i=c0
            do n =1, ncat
               v_i = v_i + vcat(n)**2 / (max(acat(n), puny))
            enddo
            v_i = v_i - m_i**2

            mu_i    = log(m_i/sqrt(c1 + v_i/m_i**2)) ! parameters for the log-normal
            sigma_i = sqrt(log(c1 + v_i/m_i**2))

            ! max thickness associated with percentile of log-normal PDF
            ! x_kmax=x997 was obtained from an optimization procedure (Dupont et al.)

            x_kmax = exp(mu_i + sqrt(c2*sigma_i)*1.9430d0)

            ! Set x_kmax to hlev of the last category where there is ice
            ! when there is no ice in the last category
            cut = x_k(ncat_i)
            do n = ncat,-1,1
               if (acat(n) < puny) then
                  cut = hin_max(n-1)
               else
                  exit
               endif
            enddo
            x_kmax = min(cut, x_kmax)

            g_k = exp(-(log(x_k) - mu_i) ** 2 / (c2 * sigma_i ** 2)) / (x_k * sigma_i * sqrt(c2 * pi))

            b_n  = exp(-(y_n - mu_b) ** 2 / (c2 * sigma_b ** 2)) / (sigma_b * sqrt(c2 * pi))

            P_x = g_k*wid_i
            P_y = b_n*wid_b

            do n =1, ncat_i
               if (x_k(n) > x_kmax) P_x(n)=c0
            enddo

            ! calculate Tb factor at t-location
            do n=1, ncat_i
               gt = (y_n <= rhoi*x_k(n)/rhow)
               tmp = merge(1,0,gt)
               ii = sum(tmp)
               if (ii == 0) then
                  tb_tmp(n) = c0
               else
                  tb_tmp(n) = max(mu_s*gravit*P_x(n)*sum(P_y(1:ii)*(rhoi*x_k(n) - rhow*y_n(1:ii))),c0)
               endif
            enddo
            Tbt(i,j) = sum(tb_tmp)*exp(-alphab * (c1 - atot))
         endif
      enddo

      if (grid_ice == "B") then
         do ij = 1, icellu
            i = indxui(ij)
            j = indxuj(ij)
            ! convert quantities to U-location
            Tbu(i,j)  = grid_neighbor_max(Tbt, i, j, 'U')
         enddo                     ! ij
      elseif (grid_ice == "C" .or. grid_ice == "CD") then
         if (present(Tbe)    .and. present(TbN)    .and. &
             present(icelle) .and. present(icelln) .and. &
             present(indxei) .and. present(indxej) .and. &
             present(indxni) .and. present(indxnj)) then

            do ij = 1, icelle
               i = indxei(ij)
               j = indxej(ij)
               ! convert quantities to E-location
                  TbE(i,j)  = grid_neighbor_max(Tbt, i, j, 'E')
            enddo
            do ij = 1, icelln
               i = indxni(ij)
               j = indxnj(ij)
               ! convert quantities to N-location
               TbN(i,j)  = grid_neighbor_max(Tbt, i, j, 'N')
            enddo

         else
            call abort_ice(subname // ' insufficient number of arguments for grid_ice:' // grid_ice)
         endif
      endif

      end subroutine seabed_stress_factor_prob

!=======================================================================
! Computes principal stresses for comparison with the theoretical
! yield curve
!
! author: Elizabeth C. Hunke, LANL

      subroutine principal_stress(nx_block,   ny_block,  &
                                  stressp,    stressm,   &
                                  stress12,   strength,  &
                                  sig1,       sig2,      &
                                  sigP)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block  ! block dimensions

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         stressp   , & ! sigma11 + sigma22
         stressm   , & ! sigma11 - sigma22
         stress12  , & ! sigma12
         strength      ! for normalization of sig1 and sig2

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out):: &
         sig1    , & ! normalized principal stress component
         sig2    , & ! normalized principal stress component
         sigP        ! internal ice pressure (N/m)

      ! local variables

      integer (kind=int_kind) :: &
         i, j

      real (kind=dbl_kind) :: &
         puny

      character(len=*), parameter :: subname = '(principal_stress)'

      call icepack_query_parameters(puny_out=puny)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      do j = 1, ny_block
      do i = 1, nx_block
         if (strength(i,j) > puny) then
            ! ice internal pressure
            sigP(i,j) = -p5*stressp(i,j)

            ! normalized principal stresses
            sig1(i,j) = (p5*(stressp(i,j) &
                      + sqrt(stressm(i,j)**2+c4*stress12(i,j)**2))) &
                      / strength(i,j)
            sig2(i,j) = (p5*(stressp(i,j) &
                      - sqrt(stressm(i,j)**2+c4*stress12(i,j)**2))) &
                      / strength(i,j)
         else
            sig1(i,j) = spval_dbl
            sig2(i,j) = spval_dbl
            sigP(i,j) = spval_dbl
         endif
      enddo
      enddo

      end subroutine principal_stress

!=======================================================================
! Compute deformations for mechanical redistribution
!
! author: Elizabeth C. Hunke, LANL
!
! 2019: subroutine created by Philippe Blain, ECCC

      subroutine deformations (nx_block,   ny_block,   &
                               icellt,                 &
                               indxti,     indxtj,     &
                               uvel,       vvel,       &
                               dxt,        dyt,        &
                               cxp,        cyp,        &
                               cxm,        cym,        &
                               tarear,                 &
                               shear,      divu,       &
                               rdg_conv,   rdg_shear )

      use ice_constants, only: p25, p5

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellt                ! no. of cells where icetmask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxti   , & ! compressed index in i-direction
         indxtj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         uvel     , & ! x-component of velocity (m/s)
         vvel     , & ! y-component of velocity (m/s)
         dxt      , & ! width of T-cell through the middle (m)
         dyt      , & ! height of T-cell through the middle (m)
         cyp      , & ! 1.5*HTE - 0.5*HTW
         cxp      , & ! 1.5*HTN - 0.5*HTS
         cym      , & ! 0.5*HTE - 1.5*HTW
         cxm      , & ! 0.5*HTN - 1.5*HTS
         tarear       ! 1/tarea

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         shear    , & ! strain rate II component (1/s)
         divu     , & ! strain rate I component, velocity divergence (1/s)
         rdg_conv , & ! convergence term for ridging (1/s)
         rdg_shear    ! shear term for ridging (1/s)

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &                       ! at each corner :
        divune, divunw, divuse, divusw            , & ! divergence
        tensionne, tensionnw, tensionse, tensionsw, & ! tension
        shearne, shearnw, shearse, shearsw        , & ! shearing
        Deltane, Deltanw, Deltase, Deltasw        , & ! Delta
        tmp                                           ! useful combination

      character(len=*), parameter :: subname = '(deformations)'

      do ij = 1, icellt
         i = indxti(ij)
         j = indxtj(ij)

         !-----------------------------------------------------------------
         ! strain rates
         ! NOTE these are actually strain rates * area  (m^2/s)
         !-----------------------------------------------------------------
         call strain_rates (nx_block,   ny_block,   &
                            i,          j,          &
                            uvel,       vvel,       &
                            dxt,        dyt,        &
                            cxp,        cyp,        &
                            cxm,        cym,        &
                            divune,     divunw,     &
                            divuse,     divusw,     &
                            tensionne,  tensionnw,  &
                            tensionse,  tensionsw,  &
                            shearne,    shearnw,    &
                            shearse,    shearsw,    &
                            Deltane,    Deltanw,    &
                            Deltase,    Deltasw     )
         !-----------------------------------------------------------------
         ! deformations for mechanical redistribution
         !-----------------------------------------------------------------
         divu(i,j) = p25*(divune + divunw + divuse + divusw) * tarear(i,j)
         tmp = p25*(Deltane + Deltanw + Deltase + Deltasw)   * tarear(i,j)
         rdg_conv(i,j)  = -min(divu(i,j),c0)
         rdg_shear(i,j) = p5*(tmp-abs(divu(i,j)))

         ! diagnostic only
         ! shear = sqrt(tension**2 + shearing**2)
         shear(i,j) = p25*tarear(i,j)*sqrt( &
                      (tensionne + tensionnw + tensionse + tensionsw)**2 + &
                      (shearne   + shearnw   + shearse   + shearsw  )**2)

      enddo                     ! ij

      end subroutine deformations

!=======================================================================
! Compute deformations for mechanical redistribution at T point
!
! author: JF Lemieux, ECCC
! Nov 2021

      subroutine deformationsCD_T (nx_block,   ny_block,   &
                                   icellt,                 &
                                   indxti,     indxtj,     &
                                   uvelE,      vvelE,      &
                                   uvelN,      vvelN,      &
                                   dxN,        dyE,        &
                                   dxT,        dyT,        &
                                   tarear,                 &
                                   shear,      divu,       &
                                   rdg_conv,   rdg_shear )

      use ice_constants, only: p5

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellt                ! no. of cells where icetmask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxti   , & ! compressed index in i-direction
         indxtj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         uvelE    , & ! x-component of velocity (m/s) at the E point
         vvelE    , & ! y-component of velocity (m/s) at the N point
         uvelN    , & ! x-component of velocity (m/s) at the E point
         vvelN    , & ! y-component of velocity (m/s) at the N point
         dxN      , & ! width of N-cell through the middle (m)
         dyE      , & ! height of E-cell through the middle (m)
         dxT      , & ! width of T-cell through the middle (m)
         dyT      , & ! height of T-cell through the middle (m)
         tarear       ! 1/tarea

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         shear    , & ! strain rate II component (1/s)
         divu     , & ! strain rate I component, velocity divergence (1/s)
         rdg_conv , & ! convergence term for ridging (1/s)
         rdg_shear    ! shear term for ridging (1/s)

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
        divT      , & ! divergence at T point
        tensionT  , & ! tension at T point
        shearT    , & ! shear at T point
        DeltaT        ! delt at T point

      real (kind=dbl_kind) :: &
        tmp           ! useful combination

      character(len=*), parameter :: subname = '(deformations_T)'

      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------

      call strain_rates_T (nx_block   ,   ny_block   , &
                           icellt     ,                &
                           indxti(:)  , indxtj  (:)  , &
                           uvelE (:,:), vvelE   (:,:), &
                           uvelN (:,:), vvelN   (:,:), &
                           dxN   (:,:), dyE     (:,:), &
                           dxT   (:,:), dyT     (:,:), &
                           divT  (:,:), tensionT(:,:), &
                           shearT(:,:), DeltaT  (:,:)  )

      do ij = 1, icellt
         i = indxti(ij)
         j = indxtj(ij)

         !-----------------------------------------------------------------
         ! deformations for mechanical redistribution
         !-----------------------------------------------------------------
         divu(i,j) = divT(i,j) * tarear(i,j)
         tmp = Deltat(i,j) * tarear(i,j)
         rdg_conv(i,j)  = -min(divu(i,j),c0)
         rdg_shear(i,j) = p5*(tmp-abs(divu(i,j)))

         ! diagnostic only
         ! shear = sqrt(tension**2 + shearing**2)
         shear(i,j) = tarear(i,j)*sqrt( tensionT(i,j)**2 + shearT(i,j)**2 )

      enddo                     ! ij

    end subroutine deformationsCD_T


!=======================================================================
! Compute deformations for mechanical redistribution at T point
!
! author: JF Lemieux, ECCC
! Nov 2021

    subroutine deformationsC_T (nx_block,   ny_block,   &
                                icellt,                 &
                                indxti,     indxtj,     &
                                uvelE,      vvelE,      &
                                uvelN,      vvelN,      &
                                dxN,        dyE,        &
                                dxT,        dyT,        &
                                tarear,     uarea,      &
                                shearU,                 &
                                shear,      divu,       &
                                rdg_conv,   rdg_shear )

      use ice_constants, only: p5

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellt                ! no. of cells where icetmask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxti   , & ! compressed index in i-direction
         indxtj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         uvelE    , & ! x-component of velocity (m/s) at the E point
         vvelE    , & ! y-component of velocity (m/s) at the N point
         uvelN    , & ! x-component of velocity (m/s) at the E point
         vvelN    , & ! y-component of velocity (m/s) at the N point
         dxN      , & ! width of N-cell through the middle (m)
         dyE      , & ! height of E-cell through the middle (m)
         dxT      , & ! width of T-cell through the middle (m)
         dyT      , & ! height of T-cell through the middle (m)
         tarear   , & ! 1/tarea
         uarea    , & ! area of u cell
         shearU       ! shearU

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         shear    , & ! strain rate II component (1/s)
         divu     , & ! strain rate I component, velocity divergence (1/s)
         rdg_conv , & ! convergence term for ridging (1/s)
         rdg_shear    ! shear term for ridging (1/s)

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
        divT      , & ! divergence at T point
        tensionT  , & ! tension at T point
        shearT    , & ! shear at T point
        DeltaT        ! delt at T point

      real (kind=dbl_kind) :: &
        tmp       , & ! useful combination
        shearTsqr     ! strain rates squared at T point

      character(len=*), parameter :: subname = '(deformations_T2)'

      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------

      call strain_rates_T (nx_block   ,   ny_block   , &
                           icellt     ,                &
                           indxti(:)  , indxtj  (:)  , &
                           uvelE (:,:), vvelE   (:,:), &
                           uvelN (:,:), vvelN   (:,:), &
                           dxN   (:,:), dyE     (:,:), &
                           dxT   (:,:), dyT     (:,:), &
                           divT  (:,:), tensionT(:,:), &
                           shearT(:,:), DeltaT  (:,:)  )

      ! DeltaT is calc by strain_rates_T but replaced by calculation below.

      do ij = 1, icellt
         i = indxti(ij)
         j = indxtj(ij)

         !-----------------------------------------------------------------
         ! deformations for mechanical redistribution
         !-----------------------------------------------------------------
         
         shearTsqr = (shearU(i  ,j  )**2 * uarea(i  ,j  )  &
                    + shearU(i  ,j-1)**2 * uarea(i  ,j-1)  &
                    + shearU(i-1,j-1)**2 * uarea(i-1,j-1)  &
                    + shearU(i-1,j  )**2 * uarea(i-1,j  )) &
                    / (uarea(i,j)+uarea(i,j-1)+uarea(i-1,j-1)+uarea(i-1,j))

         DeltaT(i,j) = sqrt(divT(i,j)**2 + e_factor*(tensionT(i,j)**2 + shearTsqr))

         divu(i,j) = divT(i,j) * tarear(i,j)
         tmp = DeltaT(i,j) * tarear(i,j)
         rdg_conv(i,j)  = -min(divu(i,j),c0)
         rdg_shear(i,j) = p5*(tmp-abs(divu(i,j)))

         ! diagnostic only...maybe we dont want to use shearTsqr here????
         ! shear = sqrt(tension**2 + shearing**2)
         shear(i,j) = tarear(i,j)*sqrt( tensionT(i,j)**2 + shearT(i,j)**2 )

      enddo                     ! ij

    end subroutine deformationsC_T

!=======================================================================
! Compute strain rates
!
! author: Elizabeth C. Hunke, LANL
!
! 2019: subroutine created by Philippe Blain, ECCC

      subroutine strain_rates (nx_block,   ny_block,   &
                               i,          j,          &
                               uvel,       vvel,       &
                               dxt,        dyt,        &
                               cxp,        cyp,        &
                               cxm,        cym,        &
                               divune,     divunw,     &
                               divuse,     divusw,     &
                               tensionne,  tensionnw,  &
                               tensionse,  tensionsw,  &
                               shearne,    shearnw,    &
                               shearse,    shearsw,    &
                               Deltane,    Deltanw,    &
                               Deltase,    Deltasw     )

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block    ! block dimensions

      integer (kind=int_kind), intent(in) :: &
         i, j                  ! indices

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         uvel     , & ! x-component of velocity (m/s)
         vvel     , & ! y-component of velocity (m/s)
         dxt      , & ! width of T-cell through the middle (m)
         dyt      , & ! height of T-cell through the middle (m)
         cyp      , & ! 1.5*HTE - 0.5*HTW
         cxp      , & ! 1.5*HTN - 0.5*HTS
         cym      , & ! 0.5*HTE - 1.5*HTW
         cxm          ! 0.5*HTN - 1.5*HTS

      real (kind=dbl_kind), intent(out):: &           ! at each corner :
        divune, divunw, divuse, divusw            , & ! divergence
        tensionne, tensionnw, tensionse, tensionsw, & ! tension
        shearne, shearnw, shearse, shearsw        , & ! shearing
        Deltane, Deltanw, Deltase, Deltasw            ! Delta

      character(len=*), parameter :: subname = '(strain_rates)'

      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------

      ! divergence  =  e_11 + e_22
      divune    = cyp(i,j)*uvel(i  ,j  ) - dyt(i,j)*uvel(i-1,j  ) &
                + cxp(i,j)*vvel(i  ,j  ) - dxt(i,j)*vvel(i  ,j-1)
      divunw    = cym(i,j)*uvel(i-1,j  ) + dyt(i,j)*uvel(i  ,j  ) &
                + cxp(i,j)*vvel(i-1,j  ) - dxt(i,j)*vvel(i-1,j-1)
      divusw    = cym(i,j)*uvel(i-1,j-1) + dyt(i,j)*uvel(i  ,j-1) &
                + cxm(i,j)*vvel(i-1,j-1) + dxt(i,j)*vvel(i-1,j  )
      divuse    = cyp(i,j)*uvel(i  ,j-1) - dyt(i,j)*uvel(i-1,j-1) &
                + cxm(i,j)*vvel(i  ,j-1) + dxt(i,j)*vvel(i  ,j  )

      ! tension strain rate  =  e_11 - e_22
      tensionne = -cym(i,j)*uvel(i  ,j  ) - dyt(i,j)*uvel(i-1,j  ) &
                +  cxm(i,j)*vvel(i  ,j  ) + dxt(i,j)*vvel(i  ,j-1)
      tensionnw = -cyp(i,j)*uvel(i-1,j  ) + dyt(i,j)*uvel(i  ,j  ) &
                +  cxm(i,j)*vvel(i-1,j  ) + dxt(i,j)*vvel(i-1,j-1)
      tensionsw = -cyp(i,j)*uvel(i-1,j-1) + dyt(i,j)*uvel(i  ,j-1) &
                +  cxp(i,j)*vvel(i-1,j-1) - dxt(i,j)*vvel(i-1,j  )
      tensionse = -cym(i,j)*uvel(i  ,j-1) - dyt(i,j)*uvel(i-1,j-1) &
                +  cxp(i,j)*vvel(i  ,j-1) - dxt(i,j)*vvel(i  ,j  )

      ! shearing strain rate  =  2*e_12
      shearne = -cym(i,j)*vvel(i  ,j  ) - dyt(i,j)*vvel(i-1,j  ) &
              -  cxm(i,j)*uvel(i  ,j  ) - dxt(i,j)*uvel(i  ,j-1)
      shearnw = -cyp(i,j)*vvel(i-1,j  ) + dyt(i,j)*vvel(i  ,j  ) &
              -  cxm(i,j)*uvel(i-1,j  ) - dxt(i,j)*uvel(i-1,j-1)
      shearsw = -cyp(i,j)*vvel(i-1,j-1) + dyt(i,j)*vvel(i  ,j-1) &
              -  cxp(i,j)*uvel(i-1,j-1) + dxt(i,j)*uvel(i-1,j  )
      shearse = -cym(i,j)*vvel(i  ,j-1) - dyt(i,j)*vvel(i-1,j-1) &
              -  cxp(i,j)*uvel(i  ,j-1) + dxt(i,j)*uvel(i  ,j  )

      ! Delta (in the denominator of zeta, eta)
      Deltane = sqrt(divune**2 + e_factor*(tensionne**2 + shearne**2))
      Deltanw = sqrt(divunw**2 + e_factor*(tensionnw**2 + shearnw**2))
      Deltasw = sqrt(divusw**2 + e_factor*(tensionsw**2 + shearsw**2))
      Deltase = sqrt(divuse**2 + e_factor*(tensionse**2 + shearse**2))

      end subroutine strain_rates

!=======================================================================
! Compute dtsd (div, tension, shear, delta) strain rates at the T point
!
! author: JF Lemieux, ECCC
! Nov 2021

      subroutine strain_rates_Tdtsd (nx_block,   ny_block, &
                                 icellt,               &
                                 indxti,     indxtj,   &
                                 uvelE,      vvelE,    &
                                 uvelN,      vvelN,    &
                                 dxN,        dyE,      &
                                 dxT,        dyT,      &
                                 divT,       tensionT, &
                                 shearT,     DeltaT    )

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, &  ! block dimensions
         icellt

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxti   , & ! compressed index in i-direction
         indxtj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         uvelE    , & ! x-component of velocity (m/s) at the E point
         vvelE    , & ! y-component of velocity (m/s) at the N point
         uvelN    , & ! x-component of velocity (m/s) at the E point
         vvelN    , & ! y-component of velocity (m/s) at the N point
         dxN      , & ! width of N-cell through the middle (m)
         dyE      , & ! height of E-cell through the middle (m)
         dxT      , & ! width of T-cell through the middle (m)
         dyT          ! height of T-cell through the middle (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out):: &
         divT     , & ! divergence at T point
         tensionT , & ! tension at T point
         shearT   , & ! shear at T point
         DeltaT       ! strain rates at the T point

      ! local variables

      integer (kind=int_kind) :: &
         ij, i, j                  ! indices

      character(len=*), parameter :: subname = '(strain_rates_Tdtsd)'

      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------

      ! compute divT, tensionT
      call strain_rates_Tdt (nx_block,   ny_block, &
                             icellt,               &
                             indxti,     indxtj,   &
                             uvelE,      vvelE,    &
                             uvelN,      vvelN,    &
                             dxN,        dyE,      &
                             dxT,        dyT,      &
                             divT,       tensionT  )

      shearT  (:,:) = c0
      deltaT  (:,:) = c0

      do ij = 1, icellt
         i = indxti(ij)
         j = indxtj(ij)

         ! shearing strain rate  =  2*e_12
         shearT(i,j) = (dxT(i,j)**2)*(uvelN(i,j)/dxN(i,j) - uvelN(i,j-1)/dxN(i,j-1)) &
                     + (dyT(i,j)**2)*(vvelE(i,j)/dyE(i,j) - vvelE(i-1,j)/dyE(i-1,j))

         ! Delta (in the denominator of zeta, eta)
         DeltaT(i,j) = sqrt(divT(i,j)**2 + e_factor*(tensionT(i,j)**2 + shearT(i,j)**2))

      enddo

      end subroutine strain_rates_Tdtsd

!=======================================================================
! Compute the dt (div, tension) strain rates at the T point
!
! author: JF Lemieux, ECCC
! Nov 2021

      subroutine strain_rates_Tdt (nx_block,   ny_block, &
                                 icellt,               &
                                 indxti,     indxtj,   &
                                 uvelE,      vvelE,    &
                                 uvelN,      vvelN,    &
                                 dxN,        dyE,      &
                                 dxT,        dyT,      &
                                 divT,       tensionT  )

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, &  ! block dimensions
         icellt

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxti   , & ! compressed index in i-direction
         indxtj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         uvelE    , & ! x-component of velocity (m/s) at the E point
         vvelE    , & ! y-component of velocity (m/s) at the N point
         uvelN    , & ! x-component of velocity (m/s) at the E point
         vvelN    , & ! y-component of velocity (m/s) at the N point
         dxN      , & ! width of N-cell through the middle (m)
         dyE      , & ! height of E-cell through the middle (m)
         dxT      , & ! width of T-cell through the middle (m)
         dyT          ! height of T-cell through the middle (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out):: &
         divT     , & ! divergence at T point
         tensionT     ! tension at T point

      ! local variables

      integer (kind=int_kind) :: &
         ij, i, j                  ! indices

      character(len=*), parameter :: subname = '(strain_rates_Tdt)'

      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------

      divT    (:,:) = c0
      tensionT(:,:) = c0

      do ij = 1, icellt
         i = indxti(ij)
         j = indxtj(ij)

         ! divergence  =  e_11 + e_22
         divT    (i,j)= dyE(i,j)*uvelE(i  ,j  ) - dyE(i-1,j)*uvelE(i-1,j  ) &
                      + dxN(i,j)*vvelN(i  ,j  ) - dxN(i,j-1)*vvelN(i  ,j-1)

         ! tension strain rate  =  e_11 - e_22
         tensionT(i,j) = (dyT(i,j)**2)*(uvelE(i,j)/dyE(i,j) - uvelE(i-1,j)/dyE(i-1,j)) &
                       - (dxT(i,j)**2)*(vvelN(i,j)/dxN(i,j) - vvelN(i,j-1)/dxN(i,j-1))

      enddo

      end subroutine strain_rates_Tdt

!=======================================================================
! Compute strain rates at the U point including boundary conditions
!
! author: JF Lemieux, ECCC
! Nov 2021

      subroutine strain_rates_U (nx_block,   ny_block,  &
                                 icellu,                &
                                 indxui,     indxuj,    &
                                 uvelE,      vvelE,     &
                                 uvelN,      vvelN,     &
                                 uvelU,      vvelU,     &
                                 dxE,        dyN,       &
                                 dxU,        dyU,       &
                                 ratiodxN,   ratiodxNr, &
                                 ratiodyE,   ratiodyEr, &
                                 epm,        npm,       &
                                 divergU,    tensionU,  &
                                 shearU,     DeltaU     )

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellu

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxui   , & ! compressed index in i-direction
         indxuj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         uvelE    , & ! x-component of velocity (m/s) at the E point
         vvelE    , & ! y-component of velocity (m/s) at the N point
         uvelN    , & ! x-component of velocity (m/s) at the E point
         vvelN    , & ! y-component of velocity (m/s) at the N point
         uvelU    , & ! x-component of velocity (m/s) interp. at U point
         vvelU    , & ! y-component of velocity (m/s) interp. at U point
         dxE      , & ! width of E-cell through the middle (m)
         dyN      , & ! height of N-cell through the middle (m)
         dxU      , & ! width of U-cell through the middle (m)
         dyU      , & ! height of U-cell through the middle (m)
         ratiodxN , & ! -dxN(i+1,j)/dxN(i,j) for BCs
         ratiodxNr, & ! -dxN(i,j)/dxN(i+1,j) for BCs
         ratiodyE , & ! -dyE(i,j+1)/dyE(i,j) for BCs
         ratiodyEr, & ! -dyE(i,j)/dyE(i,j+1) for BCs
         epm      , & ! E-cell mask
         npm          ! N-cell mask

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out):: &
         divergU  , & ! divergence at U point
         tensionU , & ! tension at U point
         shearU   , & ! shear at U point
         DeltaU       ! delt at the U point

      ! local variables

      integer (kind=int_kind) :: &
         ij, i, j                  ! indices

      real (kind=dbl_kind) :: &
        uNip1j, uNij, vEijp1, vEij, uEijp1, uEij, vNip1j, vNij

      character(len=*), parameter :: subname = '(strain_rates_U)'

      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------

      divergU (:,:) = c0
      tensionU(:,:) = c0
      shearU  (:,:) = c0
      deltaU  (:,:) = c0

      do ij = 1, icellu
         i = indxui(ij)
         j = indxuj(ij)

         uNip1j = uvelN(i+1,j) * npm(i+1,j) &
                +(npm(i,j)-npm(i+1,j)) * npm(i,j)   * ratiodxN(i,j)  * uvelN(i,j)
         uNij   = uvelN(i,j) * npm(i,j) &
                +(npm(i+1,j)-npm(i,j)) * npm(i+1,j) * ratiodxNr(i,j) * uvelN(i+1,j)
         vEijp1 = vvelE(i,j+1) * epm(i,j+1) &
                +(epm(i,j)-epm(i,j+1)) * epm(i,j)   * ratiodyE(i,j)  * vvelE(i,j)
         vEij   = vvelE(i,j) * epm(i,j) &
                +(epm(i,j+1)-epm(i,j)) * epm(i,j+1) * ratiodyEr(i,j) * vvelE(i,j+1)

         ! divergence  =  e_11 + e_22
         divergU (i,j) = dyU(i,j) * ( uNip1j - uNij ) &
                       + uvelU(i,j) * ( dyN(i+1,j) - dyN(i,j) ) &
                       + dxU(i,j) * ( vEijp1 - vEij ) &
                       + vvelU(i,j) * ( dxE(i,j+1) - dxE(i,j) )

         ! tension strain rate  =  e_11 - e_22
         tensionU(i,j) = dyU(i,j) * ( uNip1j - uNij ) &
                       - uvelU(i,j) * ( dyN(i+1,j) - dyN(i,j) ) &
                       - dxU(i,j) * ( vEijp1 - vEij ) &
                       + vvelU(i,j) * ( dxE(i,j+1) - dxE(i,j) )

         uEijp1 = uvelE(i,j+1) * epm(i,j+1) &
                +(epm(i,j)-epm(i,j+1)) * epm(i,j)   * ratiodyE(i,j)  * uvelE(i,j)
         uEij   = uvelE(i,j) * epm(i,j) &
                +(epm(i,j+1)-epm(i,j)) * epm(i,j+1) * ratiodyEr(i,j) * uvelE(i,j+1)
         vNip1j = vvelN(i+1,j) * npm(i+1,j) &
                +(npm(i,j)-npm(i+1,j)) * npm(i,j)   * ratiodxN(i,j)  * vvelN(i,j)
         vNij   = vvelN(i,j) * npm(i,j) &
                +(npm(i+1,j)-npm(i,j)) * npm(i+1,j) * ratiodxNr(i,j) * vvelN(i+1,j)

         ! shearing strain rate  =  2*e_12
         shearU(i,j)   = dxU(i,j) * ( uEijp1 - uEij ) &
                       - uvelU(i,j) * ( dxE(i,j+1) - dxE(i,j) ) &
                       + dyU(i,j) * ( vNip1j - vNij ) &
                       - vvelU(i,j) * ( dyN(i+1,j) - dyN(i,j) )

         ! Delta (in the denominator of zeta, eta)
         DeltaU(i,j)   = sqrt(divergU(i,j)**2 + e_factor*(tensionU(i,j)**2 + shearU(i,j)**2))

      enddo

      end subroutine strain_rates_U

!=======================================================================
! Computes viscosities and replacement pressure for stress
! calculations. Note that tensile strength is included here.
!
! Hibler, W. D. (1979). A dynamic thermodynamic sea ice model. J. Phys.
! Oceanogr., 9, 817-846.
!
! Konig Beatty, C. and Holland, D. M.  (2010). Modeling landfast ice by
! adding tensile strength. J. Phys. Oceanogr. 40, 185-198.
!
! Lemieux, J. F. et al. (2016). Improving the simulation of landfast ice
! by combining tensile strength and a parameterization for grounded ridges.
! J. Geophys. Res. Oceans, 121, 7354-7368.

      subroutine visc_replpress(strength, DminArea, Delta, &
                                zetax2, etax2, rep_prs, capping)

      real (kind=dbl_kind), intent(in)::  &
         strength, & !
         DminArea    !

      real (kind=dbl_kind), intent(in)::  &
         Delta   , & !
         capping     !

      real (kind=dbl_kind), intent(out):: &
         zetax2  , & ! bulk viscosity
         etax2   , & ! shear viscosity
         rep_prs     ! replacement pressure

      ! local variables
      real (kind=dbl_kind) :: &
         tmpcalc     ! temporary

      character(len=*), parameter :: subname = '(visc_replpress)'

      ! NOTE: for comp. efficiency 2 x zeta and 2 x eta are used in the code

      tmpcalc =     capping *(strength/max(Delta,DminArea))+ &
                (c1-capping)*(strength/(Delta + DminArea))
      zetax2  = (c1+Ktens)*tmpcalc
      rep_prs = (c1-Ktens)*tmpcalc*Delta
      etax2   = epp2i*zetax2

      end subroutine visc_replpress

!=======================================================================
! Do a halo update on 1 field

      subroutine dyn_haloUpdate1(halo_info, halo_info_mask, field_loc, field_type, fld1)

      use ice_boundary, only: ice_halo, ice_HaloUpdate
      use ice_domain, only: nblocks, maskhalo_dyn, halo_dynbundle
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bound

      type (ice_halo), intent(in) :: &
         halo_info  , &  ! standard unmasked halo
         halo_info_mask  ! masked halo

      integer (kind=int_kind), intent(in) :: &
         field_loc  ,  & ! field loc
         field_type      ! field_type

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), intent(inout) :: &
         fld1            ! fields to halo

      ! local variables

      integer (kind=int_kind) :: &
         iblk            ! iblock

      real (kind=dbl_kind), dimension (nx_block,ny_block,1,max_blocks) :: &
         fldbundle       ! work array for boundary updates

      character(len=*), parameter :: subname = '(dyn_haloUpdate1)'

      call ice_timer_start(timer_bound)

         if (maskhalo_dyn) then
            call ice_HaloUpdate (fld1     , halo_info_mask, &
                                 field_loc, field_type)
         else
            call ice_HaloUpdate (fld1     , halo_info     , &
                                 field_loc, field_type)
         endif

      call ice_timer_stop(timer_bound)

      end subroutine dyn_haloUpdate1

!=======================================================================
! Do a halo update on 2 fields

      subroutine dyn_haloUpdate2(halo_info, halo_info_mask, field_loc, field_type, fld1, fld2)

      use ice_boundary, only: ice_halo, ice_HaloUpdate
      use ice_domain, only: nblocks, maskhalo_dyn, halo_dynbundle
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bound

      type (ice_halo), intent(in) :: &
         halo_info  , &  ! standard unmasked halo
         halo_info_mask  ! masked halo

      integer (kind=int_kind), intent(in) :: &
         field_loc  ,  & ! field loc
         field_type      ! field_type

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), intent(inout) :: &
         fld1        , & ! fields to halo
         fld2            !

      ! local variables

      integer (kind=int_kind) :: &
         iblk            ! iblock

      real (kind=dbl_kind), dimension (nx_block,ny_block,2,max_blocks) :: &
         fldbundle       ! work array for boundary updates

      character(len=*), parameter :: subname = '(dyn_haloUpdate2)'

      call ice_timer_start(timer_bound)
      ! single process performs better without bundling fields
      if (halo_dynbundle) then

         call stack_fields(fld1, fld2, fldbundle)
         if (maskhalo_dyn) then
            call ice_HaloUpdate (fldbundle, halo_info_mask, &
                                 field_loc, field_type)
         else
            call ice_HaloUpdate (fldbundle, halo_info     , &
                                 field_loc, field_type)
         endif
         call unstack_fields(fldbundle, fld1, fld2)

      else

         if (maskhalo_dyn) then
            call ice_HaloUpdate (fld1     , halo_info_mask, &
                                 field_loc, field_type)
            call ice_HaloUpdate (fld2     , halo_info_mask, &
                                 field_loc, field_type)
         else
            call ice_HaloUpdate (fld1     , halo_info     , &
                                 field_loc, field_type)
            call ice_HaloUpdate (fld2     , halo_info     , &
                                 field_loc, field_type)
         endif

      endif
      call ice_timer_stop(timer_bound)

      end subroutine dyn_haloUpdate2

!=======================================================================
! Do a halo update on 3 fields

      subroutine dyn_haloUpdate3(halo_info, halo_info_mask, field_loc, field_type, fld1, fld2, fld3)

      use ice_boundary, only: ice_halo, ice_HaloUpdate
      use ice_domain, only: nblocks, maskhalo_dyn, halo_dynbundle
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bound

      type (ice_halo), intent(in) :: &
         halo_info  , &  ! standard unmasked halo
         halo_info_mask  ! masked halo

      integer (kind=int_kind), intent(in) :: &
         field_loc  ,  & ! field loc
         field_type      ! field_type

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), intent(inout) :: &
         fld1        , & ! fields to halo
         fld2        , & !
         fld3            !

      ! local variables

      integer (kind=int_kind) :: &
         iblk            ! iblock

      real (kind=dbl_kind), dimension (nx_block,ny_block,3,max_blocks) :: &
         fldbundle       ! work array for boundary updates

      character(len=*), parameter :: subname = '(dyn_haloUpdate3)'

      call ice_timer_start(timer_bound)
      ! single process performs better without bundling fields
      if (halo_dynbundle) then

         call stack_fields(fld1, fld2, fld3, fldbundle)
         if (maskhalo_dyn) then
            call ice_HaloUpdate (fldbundle, halo_info_mask, &
                                 field_loc, field_type)
         else
            call ice_HaloUpdate (fldbundle, halo_info     , &
                                 field_loc, field_type)
         endif
         call unstack_fields(fldbundle, fld1, fld2, fld3)

      else

         if (maskhalo_dyn) then
            call ice_HaloUpdate (fld1     , halo_info_mask, &
                                 field_loc, field_type)
            call ice_HaloUpdate (fld2     , halo_info_mask, &
                                 field_loc, field_type)
            call ice_HaloUpdate (fld3     , halo_info_mask, &
                                 field_loc, field_type)
         else
            call ice_HaloUpdate (fld1     , halo_info     , &
                                 field_loc, field_type)
            call ice_HaloUpdate (fld2     , halo_info     , &
                                 field_loc, field_type)
            call ice_HaloUpdate (fld3     , halo_info     , &
                                 field_loc, field_type)
         endif

      endif
      call ice_timer_stop(timer_bound)

      end subroutine dyn_haloUpdate3

!=======================================================================
! Do a halo update on 4 fields

      subroutine dyn_haloUpdate4(halo_info, halo_info_mask, field_loc, field_type, fld1, fld2, fld3, fld4)

      use ice_boundary, only: ice_halo, ice_HaloUpdate
      use ice_domain, only: nblocks, maskhalo_dyn, halo_dynbundle
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bound

      type (ice_halo), intent(in) :: &
         halo_info  , &  ! standard unmasked halo
         halo_info_mask  ! masked halo

      integer (kind=int_kind), intent(in) :: &
         field_loc,    & ! field loc
         field_type      ! field_type

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), intent(inout) :: &
         fld1        , & ! fields to halo
         fld2        , & !
         fld3        , & !
         fld4            !

      ! local variables

      integer (kind=int_kind) :: &
         iblk            ! iblock

      real (kind=dbl_kind), dimension (nx_block,ny_block,4,max_blocks) :: &
         fldbundle       ! work array for boundary updates

      character(len=*), parameter :: subname = '(dyn_haloUpdate4)'

      call ice_timer_start(timer_bound)
      ! single process performs better without bundling fields
      if (halo_dynbundle) then

         call stack_fields(fld1, fld2, fld3, fld4, fldbundle)
         if (maskhalo_dyn) then
            call ice_HaloUpdate (fldbundle, halo_info_mask, &
                                 field_loc, field_type)
         else
            call ice_HaloUpdate (fldbundle, halo_info     , &
                                 field_loc, field_type)
         endif
         call unstack_fields(fldbundle, fld1, fld2, fld3, fld4)

      else

         if (maskhalo_dyn) then
            call ice_HaloUpdate (fld1     , halo_info_mask, &
                                 field_loc, field_type)
            call ice_HaloUpdate (fld2     , halo_info_mask, &
                                 field_loc, field_type)
            call ice_HaloUpdate (fld3     , halo_info_mask, &
                                 field_loc, field_type)
            call ice_HaloUpdate (fld4     , halo_info_mask, &
                                 field_loc, field_type)
         else
            call ice_HaloUpdate (fld1     , halo_info     , &
                                 field_loc, field_type)
            call ice_HaloUpdate (fld2     , halo_info     , &
                                 field_loc, field_type)
            call ice_HaloUpdate (fld3     , halo_info     , &
                                 field_loc, field_type)
            call ice_HaloUpdate (fld4     , halo_info     , &
                                 field_loc, field_type)
         endif

      endif
      call ice_timer_stop(timer_bound)

      end subroutine dyn_haloUpdate4

!=======================================================================
! Do a halo update on 5 fields

      subroutine dyn_haloUpdate5(halo_info, halo_info_mask, field_loc, field_type, fld1, fld2, fld3, fld4, fld5)

      use ice_boundary, only: ice_halo, ice_HaloUpdate
      use ice_domain, only: nblocks, maskhalo_dyn, halo_dynbundle
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bound

      type (ice_halo), intent(in) :: &
         halo_info  , &  ! standard unmasked halo
         halo_info_mask  ! masked halo

      integer (kind=int_kind), intent(in) :: &
         field_loc  ,  & ! field loc
         field_type      ! field_type

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks), intent(inout) :: &
         fld1        , & ! fields to halo
         fld2        , & !
         fld3        , & !
         fld4        , & !
         fld5            !

      ! local variables

      integer (kind=int_kind) :: &
         iblk            ! iblock

      real (kind=dbl_kind), dimension (nx_block,ny_block,5,max_blocks) :: &
         fldbundle       ! work array for boundary updates

      character(len=*), parameter :: subname = '(dyn_haloUpdate5)'

      call ice_timer_start(timer_bound)
      ! single process performs better without bundling fields
      if (halo_dynbundle) then

         call stack_fields(fld1, fld2, fld3, fld4, fld5, fldbundle)
         if (maskhalo_dyn) then
            call ice_HaloUpdate (fldbundle, halo_info_mask, &
                                 field_loc, field_type)
         else
            call ice_HaloUpdate (fldbundle, halo_info     , &
                                 field_loc, field_type)
         endif
         call unstack_fields(fldbundle, fld1, fld2, fld3, fld4, fld5)

      else

         if (maskhalo_dyn) then
            call ice_HaloUpdate (fld1     , halo_info_mask, &
                                 field_loc, field_type)
            call ice_HaloUpdate (fld2     , halo_info_mask, &
                                 field_loc, field_type)
            call ice_HaloUpdate (fld3     , halo_info_mask, &
                                 field_loc, field_type)
            call ice_HaloUpdate (fld4     , halo_info_mask, &
                                 field_loc, field_type)
            call ice_HaloUpdate (fld5     , halo_info_mask, &
                                 field_loc, field_type)
         else
            call ice_HaloUpdate (fld1     , halo_info     , &
                                 field_loc, field_type)
            call ice_HaloUpdate (fld2     , halo_info     , &
                                 field_loc, field_type)
            call ice_HaloUpdate (fld3     , halo_info     , &
                                 field_loc, field_type)
            call ice_HaloUpdate (fld4     , halo_info     , &
                                 field_loc, field_type)
            call ice_HaloUpdate (fld5     , halo_info     , &
                                 field_loc, field_type)
         endif

      endif
      call ice_timer_stop(timer_bound)

      end subroutine dyn_haloUpdate5

!=======================================================================
! Load fields into array for boundary updates

      subroutine stack_fields2(fld1, fld2, fldbundle)

      use ice_domain, only: nblocks
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bundbound

      real (kind=dbl_kind), dimension (:,:,:), intent(in) :: &
         fld1    , & ! fields to stack
         fld2        !

      real (kind=dbl_kind), dimension (:,:,:,:), intent(out) :: &
         fldbundle   ! work array for boundary updates (i,j,n,iblk)

      ! local variables

      integer (kind=int_kind) :: &
         iblk        ! block index

      character(len=*), parameter :: subname = '(stack_fields2)'

      call ice_timer_start(timer_bundbound)
      !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
      do iblk = 1, nblocks
         fldbundle(:,:,1,iblk) = fld1(:,:,iblk)
         fldbundle(:,:,2,iblk) = fld2(:,:,iblk)
      enddo
      !$OMP END PARALLEL DO
      call ice_timer_stop(timer_bundbound)

      end subroutine stack_fields2

!=======================================================================
! Load fields into array for boundary updates

      subroutine stack_fields3(fld1, fld2, fld3, fldbundle)

      use ice_domain, only: nblocks
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bundbound

      real (kind=dbl_kind), dimension (:,:,:), intent(in) :: &
         fld1    , & ! fields to stack
         fld2    , & !
         fld3        !

      real (kind=dbl_kind), dimension (:,:,:,:), intent(out) :: &
         fldbundle   ! work array for boundary updates (i,j,n,iblk)

      ! local variables

      integer (kind=int_kind) :: &
         iblk        ! block index

      character(len=*), parameter :: subname = '(stack_fields3)'

      call ice_timer_start(timer_bundbound)
      !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
      do iblk = 1, nblocks
         fldbundle(:,:,1,iblk) = fld1(:,:,iblk)
         fldbundle(:,:,2,iblk) = fld2(:,:,iblk)
         fldbundle(:,:,3,iblk) = fld3(:,:,iblk)
      enddo
      !$OMP END PARALLEL DO
      call ice_timer_stop(timer_bundbound)

      end subroutine stack_fields3

!=======================================================================
! Load fields into array for boundary updates

      subroutine stack_fields4(fld1, fld2, fld3, fld4, fldbundle)

      use ice_domain, only: nblocks
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bundbound

      real (kind=dbl_kind), dimension (:,:,:), intent(in) :: &
         fld1    , & ! fields to stack
         fld2    , & !
         fld3    , & !
         fld4        !

      real (kind=dbl_kind), dimension (:,:,:,:), intent(out) :: &
         fldbundle   ! work array for boundary updates (i,j,n,iblk)

      ! local variables

      integer (kind=int_kind) :: &
         iblk        ! block index

      character(len=*), parameter :: subname = '(stack_fields4)'

      call ice_timer_start(timer_bundbound)
      !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
      do iblk = 1, nblocks
         fldbundle(:,:,1,iblk) = fld1(:,:,iblk)
         fldbundle(:,:,2,iblk) = fld2(:,:,iblk)
         fldbundle(:,:,3,iblk) = fld3(:,:,iblk)
         fldbundle(:,:,4,iblk) = fld4(:,:,iblk)
      enddo
      !$OMP END PARALLEL DO
      call ice_timer_stop(timer_bundbound)

      end subroutine stack_fields4

!=======================================================================
! Load fields into array for boundary updates

      subroutine stack_fields5(fld1, fld2, fld3, fld4, fld5, fldbundle)

      use ice_domain, only: nblocks
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bundbound

      real (kind=dbl_kind), dimension (:,:,:), intent(in) :: &
         fld1    , & ! fields to stack
         fld2    , & !
         fld3    , & !
         fld4    , & !
         fld5        !

      real (kind=dbl_kind), dimension (:,:,:,:), intent(out) :: &
         fldbundle   ! work array for boundary updates (i,j,n,iblk)

      ! local variables

      integer (kind=int_kind) :: &
         iblk        ! block index

      character(len=*), parameter :: subname = '(stack_fields5)'

      call ice_timer_start(timer_bundbound)
      !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
      do iblk = 1, nblocks
         fldbundle(:,:,1,iblk) = fld1(:,:,iblk)
         fldbundle(:,:,2,iblk) = fld2(:,:,iblk)
         fldbundle(:,:,3,iblk) = fld3(:,:,iblk)
         fldbundle(:,:,4,iblk) = fld4(:,:,iblk)
         fldbundle(:,:,5,iblk) = fld5(:,:,iblk)
      enddo
      !$OMP END PARALLEL DO
      call ice_timer_stop(timer_bundbound)

      end subroutine stack_fields5

!=======================================================================
! Unload fields from array after boundary updates

      subroutine unstack_fields2(fldbundle, fld1, fld2)

      use ice_domain, only: nblocks
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bundbound

      real (kind=dbl_kind), dimension (:,:,:,:), intent(in) :: &
         fldbundle   ! work array for boundary updates (i,j,n,iblk)

      real (kind=dbl_kind), dimension (:,:,:), intent(out) :: &
         fld1    , & ! fields to unstack
         fld2        !

      ! local variables

      integer (kind=int_kind) :: &
         iblk        ! block index

      character(len=*), parameter :: subname = '(unstack_fields2)'

      call ice_timer_start(timer_bundbound)
      !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
      do iblk = 1, nblocks
         fld1(:,:,iblk) = fldbundle(:,:,1,iblk)
         fld2(:,:,iblk) = fldbundle(:,:,2,iblk)
      enddo
      !$OMP END PARALLEL DO
      call ice_timer_stop(timer_bundbound)

      end subroutine unstack_fields2

!=======================================================================
! Unload fields from array after boundary updates

      subroutine unstack_fields3(fldbundle, fld1, fld2, fld3)

      use ice_domain, only: nblocks
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bundbound

      real (kind=dbl_kind), dimension (:,:,:,:), intent(in) :: &
         fldbundle   ! work array for boundary updates (i,j,n,iblk)

      real (kind=dbl_kind), dimension (:,:,:), intent(out) :: &
         fld1    , & ! fields to unstack
         fld2    , & !
         fld3        !

      ! local variables

      integer (kind=int_kind) :: &
         iblk        ! block index

      character(len=*), parameter :: subname = '(unstack_fields3)'

      call ice_timer_start(timer_bundbound)
      !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
      do iblk = 1, nblocks
         fld1(:,:,iblk) = fldbundle(:,:,1,iblk)
         fld2(:,:,iblk) = fldbundle(:,:,2,iblk)
         fld3(:,:,iblk) = fldbundle(:,:,3,iblk)
      enddo
      !$OMP END PARALLEL DO
      call ice_timer_stop(timer_bundbound)

      end subroutine unstack_fields3

!=======================================================================
! Unload fields from array after boundary updates

      subroutine unstack_fields4(fldbundle, fld1, fld2, fld3, fld4)

      use ice_domain, only: nblocks
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bundbound

      real (kind=dbl_kind), dimension (:,:,:,:), intent(in) :: &
         fldbundle   ! work array for boundary updates (i,j,n,iblk)

      real (kind=dbl_kind), dimension (:,:,:), intent(out) :: &
         fld1    , & ! fields to unstack
         fld2    , & !
         fld3    , & !
         fld4        !

      ! local variables

      integer (kind=int_kind) :: &
         iblk        ! block index

      character(len=*), parameter :: subname = '(unstack_fields4)'

      call ice_timer_start(timer_bundbound)
      !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
      do iblk = 1, nblocks
         fld1(:,:,iblk) = fldbundle(:,:,1,iblk)
         fld2(:,:,iblk) = fldbundle(:,:,2,iblk)
         fld3(:,:,iblk) = fldbundle(:,:,3,iblk)
         fld4(:,:,iblk) = fldbundle(:,:,4,iblk)
      enddo
      !$OMP END PARALLEL DO
      call ice_timer_stop(timer_bundbound)

      end subroutine unstack_fields4

!=======================================================================
! Unload fields from array after boundary updates

      subroutine unstack_fields5(fldbundle, fld1, fld2, fld3, fld4, fld5)

      use ice_domain, only: nblocks
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bundbound

      real (kind=dbl_kind), dimension (:,:,:,:), intent(in) :: &
         fldbundle   ! work array for boundary updates (i,j,n,iblk)

      real (kind=dbl_kind), dimension (:,:,:), intent(out) :: &
         fld1    , & ! fields to unstack
         fld2    , & !
         fld3    , & !
         fld4    , & !
         fld5        !

      ! local variables

      integer (kind=int_kind) :: &
         iblk        ! block index

      character(len=*), parameter :: subname = '(unstack_fields5)'

      call ice_timer_start(timer_bundbound)
      !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
      do iblk = 1, nblocks
         fld1(:,:,iblk) = fldbundle(:,:,1,iblk)
         fld2(:,:,iblk) = fldbundle(:,:,2,iblk)
         fld3(:,:,iblk) = fldbundle(:,:,3,iblk)
         fld4(:,:,iblk) = fldbundle(:,:,4,iblk)
         fld5(:,:,iblk) = fldbundle(:,:,5,iblk)
      enddo
      !$OMP END PARALLEL DO
      call ice_timer_stop(timer_bundbound)

      end subroutine unstack_fields5

!=======================================================================

      end module ice_dyn_shared

!=======================================================================
