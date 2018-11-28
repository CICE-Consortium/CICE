!=======================================================================

! Elastic-viscous-plastic sea ice dynamics model code shared with other
! approaches
!
! author: Elizabeth C. Hunke, LANL
!
! 2013: Split from ice_dyn_evp.F90 by Elizabeth Hunke

      module ice_dyn_shared

      use ice_kinds_mod
      use ice_communicate, only: my_task, master_task
      use ice_constants, only: c0, c1, p01, p001
      use ice_blocks, only: nx_block, ny_block
      use ice_domain_size, only: max_blocks
      use ice_fileunits, only: nu_diag
      use ice_exit, only: abort_ice
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_parameters

      implicit none
      private
      public :: init_evp, set_evp_parameters, stepu, principal_stress, &
                dyn_prep1, dyn_prep2, dyn_finish, basal_stress_coeff,  &
                alloc_dyn_shared

      ! namelist parameters

      integer (kind=int_kind), public :: &
         kdyn       , & ! type of dynamics ( -1, 0 = off, 1 = evp, 2 = eap )
         kridge     , & ! set to "-1" to turn off ridging
         ktransport , & ! set to "-1" to turn off transport
         ndte         ! number of subcycles:  ndte=dt/dte

      character (len=char_len), public :: &
         coriolis     ! 'constant', 'zero', or 'latitude'

      logical (kind=log_kind), public :: &
         revised_evp  ! if true, use revised evp procedure

      ! other EVP parameters

      character (len=char_len), public :: & 
         yield_curve  ! 'ellipse' ('teardrop' needs further testing)
                                                                      ! 
      real (kind=dbl_kind), parameter, public :: &
         eyc = 0.36_dbl_kind, &
                         ! coefficient for calculating the parameter E
         cosw = c1   , & ! cos(ocean turning angle)  ! turning angle = 0
         sinw = c0   , & ! sin(ocean turning angle)  ! turning angle = 0
         a_min = p001, & ! minimum ice area
         m_min = p01     ! minimum ice mass (kg/m^2)

      real (kind=dbl_kind), public :: &
         revp     , & ! 0 for classic EVP, 1 for revised EVP
         e_ratio  , & ! e = EVP ellipse aspect ratio 
         ecci     , & ! 1/e^2
         dtei     , & ! 1/dte, where dte is subcycling timestep (1/s)
!         dte2T    , & ! dte/2T
         denom1       ! constants for stress equation

      real (kind=dbl_kind), public :: & ! Bouillon et al relaxation constants
         arlx     , & ! alpha for stressp
         arlx1i   , & ! (inverse of alpha) for stressp
         brlx         ! beta   for momentum

      real (kind=dbl_kind), allocatable, public :: & 
         fcor_blk(:,:,:)   ! Coriolis parameter (1/s)

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
         uvel_init, & ! x-component of velocity (m/s), beginning of timestep
         vvel_init    ! y-component of velocity (m/s), beginning of timestep
         
      ! ice isotropic tensile strength parameter
      real (kind=dbl_kind), public :: &
         Ktens         ! T=Ktens*P (tensile strength: see Konig and Holland, 2010)   

      logical (kind=log_kind), public :: &
         basalstress   ! if true, basal stress for landfast on

      real (kind=dbl_kind), public :: &
         k1            ! 1st free parameter for landfast parameterization

!=======================================================================

      contains

!=======================================================================
!
! Allocate space for all variables 
!
      subroutine alloc_dyn_shared

      integer (int_kind) :: ierr

      allocate( &
         uvel_init (nx_block,ny_block,max_blocks), & ! x-component of velocity (m/s), beginning of timestep
         vvel_init (nx_block,ny_block,max_blocks), & ! y-component of velocity (m/s), beginning of timestep
         stat=ierr)
      if (ierr/=0) call abort_ice('(alloc_dyn_shared): Out of memory')

      end subroutine alloc_dyn_shared

!=======================================================================

! Initialize parameters and variables needed for the evp dynamics
! author: Elizabeth C. Hunke, LANL

      subroutine init_evp (dt)

      use ice_blocks, only: nx_block, ny_block
      use ice_constants, only: c0, c2, omega
      use ice_domain, only: nblocks
      use ice_domain_size, only: max_blocks
      use ice_flux, only: rdg_conv, rdg_shear, iceumask, &
          stressp_1, stressp_2, stressp_3, stressp_4, &
          stressm_1, stressm_2, stressm_3, stressm_4, &
          stress12_1, stress12_2, stress12_3, stress12_4
      use ice_state, only: uvel, vvel, divu, shear
      use ice_grid, only: ULAT

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) :: &
         i, j, &
         iblk            ! block index

      character(len=*), parameter :: subname = '(init_evp)'

      call set_evp_parameters (dt)

      if (my_task == master_task) then
         write(nu_diag,*) 'dt  = ',dt
         write(nu_diag,*) 'dte = ',dt/real(ndte,kind=dbl_kind)
         write(nu_diag,*) 'tdamp =', eyc*dt
      endif

      allocate(fcor_blk(nx_block,ny_block,max_blocks))

      !$OMP PARALLEL DO PRIVATE(iblk,i,j)
      do iblk = 1, nblocks
      do j = 1, ny_block
      do i = 1, nx_block

         ! velocity
         uvel(i,j,iblk) = c0    ! m/s
         vvel(i,j,iblk) = c0    ! m/s

         ! strain rates
         divu (i,j,iblk) = c0
         shear(i,j,iblk) = c0
         rdg_conv (i,j,iblk) = c0
         rdg_shear(i,j,iblk) = c0

         ! Coriolis parameter
         if (trim(coriolis) == 'constant') then
            fcor_blk(i,j,iblk) = 1.46e-4_dbl_kind ! Hibler 1979, N. Hem; 1/s
         else if (trim(coriolis) == 'zero') then
            fcor_blk(i,j,iblk) = 0.0
         else
            fcor_blk(i,j,iblk) = c2*omega*sin(ULAT(i,j,iblk)) ! 1/s
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

         ! ice extent mask on velocity points
         iceumask(i,j,iblk) = .false.

      enddo                     ! i
      enddo                     ! j
      enddo                     ! iblk
      !$OMP END PARALLEL DO

      end subroutine init_evp

!=======================================================================

! Set parameters needed for the evp dynamics.
! Note: This subroutine is currently called only during initialization.
!       If the dynamics time step can vary during runtime, it should
!        be called whenever the time step changes.
!
! author: Elizabeth C. Hunke, LANL

      subroutine set_evp_parameters (dt)

      use ice_communicate, only: my_task, master_task
      use ice_constants, only: p25, c1, c2, p5
      use ice_domain, only: distrb_info
      use ice_global_reductions, only: global_minval
      use ice_grid, only: dxt, dyt, tmask

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      !real (kind=dbl_kind) :: &
         !dte         , & ! subcycling timestep for EVP dynamics, s
         !ecc         , & ! (ratio of major to minor ellipse axes)^2
         !tdamp2          ! 2*(wave damping time scale T)

      character(len=*), parameter :: subname = '(set_evp_parameters)'

      ! elastic time step
      !dte = dt/real(ndte,kind=dbl_kind)        ! s
      !dtei = c1/dte              ! 1/s
      dtei = real(ndte,kind=dbl_kind)/dt 

      ! major/minor axis length ratio, squared
      !ecc  = e_ratio**2
      !ecci = c1/ecc               ! 1/ecc
      ecci = c1/e_ratio**2               ! 1/ecc

      ! constants for stress equation
      !tdamp2 = c2*eyc*dt                    ! s
      !dte2T = dte/tdamp2    or c1/(c2*eyc*real(ndte,kind=dbl_kind))               ! ellipse (unitless)

      if (revised_evp) then       ! Bouillon et al, Ocean Mod 2013
         revp   = c1
         denom1 = c1
         arlx1i = c1/arlx
      else                        ! Hunke, JCP 2013 with modified stress eq
         revp   = c0
         !arlx1i = dte2T
         !arlx   = c1/arlx1i
         !brlx   = dt*dtei
         arlx   = c2*eyc*real(ndte,kind=dbl_kind)
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
                            strairxT,  strairyT, & 
                            strairx,   strairy,  & 
                            tmass,     icetmask)

      use ice_constants, only: c0

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ilo,ihi,jlo,jhi       ! beginning and end of physical domain

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         aice    , & ! concentration of ice
         vice    , & ! volume per unit area of ice          (m)
         vsno    , & ! volume per unit area of snow         (m)
         strairxT, & ! stress on ice by air, x-direction
         strairyT    ! stress on ice by air, y-direction

      logical (kind=log_kind), dimension (nx_block,ny_block), intent(in) :: &
         tmask       ! land/boundary mask, thickness (T-cell)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out) :: &
         strairx , & ! stress on ice by air, x-direction
         strairy , & ! stress on ice by air, y-direction
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
      ! prep to convert to U grid
      !-----------------------------------------------------------------
         ! these quantities include the factor of aice needed for 
         ! correct treatment of free drift
         strairx(i,j) = strairxT(i,j)
         strairy(i,j) = strairyT(i,j)

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

      use ice_constants, only: c0, c1

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ilo,ihi,jlo,jhi       ! beginning and end of physical domain

      integer (kind=int_kind), intent(out) :: &
         icellt   , & ! no. of cells where icetmask = 1
         icellu       ! no. of cells where iceumask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(out) :: &
         indxti   , & ! compressed index in i-direction
         indxtj   , & ! compressed index in j-direction
         indxui   , & ! compressed index in i-direction
         indxuj       ! compressed index in j-direction

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
         Tbu,      & ! coefficient for basal stress (N/m^2)
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
         taubx   , & ! basal stress, x-direction (N/m^2)
         tauby       ! basal stress, y-direction (N/m^2)

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

#ifdef coupled
      real (kind=dbl_kind) :: gravit
#endif

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

         if (revp==1) then               ! revised evp
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
         else if (icetmask(i,j)==0) then ! classic evp
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
         endif                  ! revp
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

         ! ice extent mask (U-cells)
         iceumask_old(i,j) = iceumask(i,j) ! save
         iceumask(i,j) = (umask(i,j)) .and. (aiu  (i,j) > a_min) & 
                                      .and. (umass(i,j) > m_min)

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

#ifdef coupled
      call icepack_query_parameters(gravit_out=gravit)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)
#endif

      do ij = 1, icellu
         i = indxui(ij)
         j = indxuj(ij)

         umassdti(i,j) = umass(i,j)/dt ! kg/m^2 s

         fm(i,j) = fcor(i,j)*umass(i,j)   ! Coriolis * mass

         ! for ocean stress
         waterx(i,j) = uocn(i,j)*cosw - vocn(i,j)*sinw*sign(c1,fm(i,j))
         watery(i,j) = vocn(i,j)*cosw + uocn(i,j)*sinw*sign(c1,fm(i,j))

         ! combine tilt with wind stress
#ifndef coupled
         ! calculate tilt from geostrophic currents if needed
         strtltx(i,j) = -fm(i,j)*vocn(i,j)
         strtlty(i,j) =  fm(i,j)*uocn(i,j)
#else
         strtltx(i,j) = -gravit*umass(i,j)*ss_tltx(i,j)
         strtlty(i,j) = -gravit*umass(i,j)*ss_tlty(i,j)
#endif
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
                        ksub,                 &
                        aiu,        str,      &
                        uocn,       vocn,     &
                        waterx,     watery,   &
                        forcex,     forcey,   &
                        umassdti,   fm,       &
                        uarear,               &
                        strocnx,    strocny,  &
                        strintx,    strinty,  &
                        taubx,      tauby,    &
                        uvel_init,  vvel_init,&
                        uvel,       vvel,     &
                        Tbu)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellu,             & ! total count when iceumask is true
         ksub                  ! subcycling iteration

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxui  , & ! compressed index in i-direction
         indxuj      ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         Tbu,      & ! coefficient for basal stress (N/m^2)
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
         str

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         uvel    , & ! x-component of velocity (m/s)
         vvel        ! y-component of velocity (m/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         strocnx , & ! ice-ocean stress, x-direction
         strocny , & ! ice-ocean stress, y-direction
         strintx , & ! divergence of internal ice stress, x (N/m^2)
         strinty , & ! divergence of internal ice stress, y (N/m^2)
         taubx   , & ! basal stress, x-direction (N/m^2)
         tauby       ! basal stress, y-direction (N/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         Cw                   ! ocean-ice neutral drag coefficient

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: &
         uold, vold        , & ! old-time uvel, vvel
         vrel              , & ! relative ice-ocean velocity
         cca,ccb,ab2,cc1,cc2,& ! intermediate variables
         taux, tauy        , & ! part of ocean stress term
         Cb                , & ! complete basal stress coeff
         rhow                  !

      real (kind=dbl_kind) :: &
         u0 = 5e-5_dbl_kind    ! residual velocity for basal stress (m/s)
         
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
      
         Cb  = Tbu(i,j) / (sqrt(uold**2 + vold**2) + u0) ! for basal stress
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

      !-----------------------------------------------------------------
      ! ocean-ice stress for coupling
      ! here, strocn includes the factor of aice
      !-----------------------------------------------------------------
         strocnx(i,j) = taux
         strocny(i,j) = tauy
         
      ! calculate basal stress component for outputs
         if (ksub == ndte) then ! on last subcycling iteration
          if ( basalstress ) then
           taubx(i,j) = -uvel(i,j)*Tbu(i,j) / (sqrt(uold**2 + vold**2) + u0)
           tauby(i,j) = -vvel(i,j)*Tbu(i,j) / (sqrt(uold**2 + vold**2) + u0)
          endif
         endif

      enddo                     ! ij

      end subroutine stepu

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
                             strintx,  strinty,  &
                             strairx,  strairy,  &
                             strocnx,  strocny,  &
                             strocnxT, strocnyT) 

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
         fm      , & ! Coriolis param. * mass in U-cell (kg/s)
         strintx , & ! divergence of internal ice stress, x (N/m^2)
         strinty , & ! divergence of internal ice stress, y (N/m^2)
         strairx , & ! stress on ice by air, x-direction
         strairy     ! stress on ice by air, y-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         strocnx , & ! ice-ocean stress, x-direction
         strocny , & ! ice-ocean stress, y-direction
         strocnxT, & ! ice-ocean stress, x-direction
         strocnyT    ! ice-ocean stress, y-direction

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij

      real (kind=dbl_kind) :: vrel, rhow
      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         Cw                   ! ocean-ice neutral drag coefficient 

      character(len=*), parameter :: subname = '(dyn_finish)'

      call icepack_query_parameters(rhow_out=rhow)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      do j = 1, ny_block
      do i = 1, nx_block
         strocnxT(i,j) = c0
         strocnyT(i,j) = c0
      enddo
      enddo

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

         ! Prepare to convert to T grid
         ! divide by aice for coupling
         strocnxT(i,j) = strocnx(i,j) / aiu(i,j)
         strocnyT(i,j) = strocny(i,j) / aiu(i,j)
      enddo

      end subroutine dyn_finish

!=======================================================================
! Computes basal stress Tbu coefficients (landfast ice)
!
! Lemieux, J. F., B. Tremblay, F. Dupont, M. Plante, G.C. Smith, D. Dumont (2015). 
! A basal stress parameterization form modeling landfast ice, J. Geophys. Res. 
! Oceans, 120, 3157-3173.
!
! Lemieux, J. F., F. Dupont, P. Blain, F. Roy, G.C. Smith, G.M. Flato (2016). 
! Improving the simulation of landfast ice by combining tensile strength and a
! parameterization for grounded ridges, J. Geophys. Res. Oceans, 121.
!
! author: JF Lemieux, Philippe Blain (ECCC)
!
! note: Tbu is a part of the Cb as defined in Lemieux et al. 2015 and 2016.
!
      subroutine basal_stress_coeff (nx_block, ny_block,         &
                                     icellu,                     &
                                     indxui,   indxuj,           &
                                     vice,     aice,             &
                                     hwater,   Tbu)

      use ice_constants, only: c0, c1                                     
                                     
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, &  ! block dimensions
         icellu                 ! no. of cells where icetmask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxui   , & ! compressed index in i-direction
         indxuj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         aice    , & ! concentration of ice at tracer location
         vice    , & ! volume per unit area of ice at tracer location
         hwater      ! water depth at tracer location

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         Tbu         ! coefficient for basal stress (N/m^2)

      real (kind=dbl_kind) :: &
         au,  & ! concentration of ice at u location
         hu,  & ! volume per unit area of ice at u location (mean thickness)
         hwu, & ! water depth at u location
         hcu, & ! critical thickness at u location
         k2 = 15.0_dbl_kind , &  ! second free parameter (N/m^3) for landfast parametrization 
         alphab = 20.0_dbl_kind  ! alphab=Cb factor in Lemieux et al 2015

      integer (kind=int_kind) :: &
         i, j, ij

      character(len=*), parameter :: subname = '(basal_stress_coeff)'

      do ij = 1, icellu
         i = indxui(ij)
         j = indxuj(ij)

         ! convert quantities to u-location
         au  = max(aice(i,j),aice(i+1,j),aice(i,j+1),aice(i+1,j+1))
         hwu = min(hwater(i,j),hwater(i+1,j),hwater(i,j+1),hwater(i+1,j+1))
         hu  = max(vice(i,j),vice(i+1,j),vice(i,j+1),vice(i+1,j+1))

         ! 1- calculate critical thickness
         hcu = au * hwu / k1

         ! 2- calculate basal stress factor                    
         Tbu(i,j) = k2 * max(c0,(hu - hcu)) * exp(-alphab * (c1 - au))

      enddo                     ! ij

      end subroutine basal_stress_coeff

!=======================================================================

! Computes principal stresses for comparison with the theoretical
! yield curve; northeast values
!
! author: Elizabeth C. Hunke, LANL

      subroutine principal_stress(nx_block,   ny_block,  &
                                  stressp_1,  stressm_1, &
                                  stress12_1, strength,  &
                                  sig1,       sig2,      &
                                  sigP)

      use ice_constants, only: spval_dbl, p5, c4

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block  ! block dimensions

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         stressp_1 , & ! sigma11 + sigma22
         stressm_1 , & ! sigma11 - sigma22
         stress12_1, & ! sigma12
         strength      ! for normalization of sig1 and sig2

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out):: &
         sig1    , & ! normalized principal stress component
         sig2    , & ! normalized principal stress component
         sigP        ! internal ice pressure (N/m)

      ! local variables

      integer (kind=int_kind) :: i, j

      real (kind=dbl_kind) :: puny

      character(len=*), parameter :: subname = '(principal_stress)'

      call icepack_query_parameters(puny_out=puny)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      do j = 1, ny_block
      do i = 1, nx_block
         if (strength(i,j) > puny) then
            ! ice internal pressure          
            sigP(i,j) = -p5*stressp_1(i,j) 
            
            ! normalized principal stresses
            sig1(i,j) = (p5*(stressp_1(i,j) &
                      + sqrt(stressm_1(i,j)**2+c4*stress12_1(i,j)**2))) &
                      / strength(i,j)
            sig2(i,j) = (p5*(stressp_1(i,j) &
                      - sqrt(stressm_1(i,j)**2+c4*stress12_1(i,j)**2))) &
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

      end module ice_dyn_shared

!=======================================================================
