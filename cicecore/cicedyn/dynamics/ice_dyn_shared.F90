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
      use ice_constants, only: c0, c1, c2, c3, c4, c6, c1p5
      use ice_constants, only: omega, spval_dbl, p01, p001, p5
      use ice_blocks, only: nx_block, ny_block
      use ice_domain_size, only: max_blocks
      use ice_fileunits, only: nu_diag
      use ice_exit, only: abort_ice
      ! The ragged layout lives in ice_grid now -- it is derived from the grid
      ! masks, and keeping it there lets the grid metrics that are compressed
      ! once per run be packed where they are created.  Re-exported here so the
      ! solvers that import it from ice_dyn_shared are unchanged.
      use ice_grid, only: grid_ice,                                   &
                          dynOff, dynCnt, dynTot, dynMaxCnt,          &
                          dynMap, dynNbr, dynIJ,                      &
                          set_ragged_layout, set_neighbours,          &
                          set_slot_cells, pack_to_slots,              &
                          unpack_from_slots,                          &
                          nbr_W, nbr_SW, nbr_S, nbr_E, nbr_NE, nbr_N, &
                          nbr_NW, nbr_SE, nbr_num,                     &
                          uarear_slot, dxN_slot, dyE_slot, uarea_slot, &
                          tarea_slot, earea_slot, narea_slot, hm_slot, &
                          dxU_slot, dyT_slot, dxT_slot, dyU_slot, &
                          dxE_slot, dyN_slot, epm_slot, npm_slot, uvm_slot, &
                          uareaavg_slot, wsumN_slot, wsumE_slot, &
                          wsumNW_slot, wsumSE_slot, wsumNE_slot
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_parameters

      implicit none
      private

      ! re-exported from ice_grid, which owns the ragged layout
      public :: dynOff, dynCnt, dynTot, dynMaxCnt, dynMap, dynNbr, dynIJ
      public :: nbr_W, nbr_SW, nbr_S, nbr_E, nbr_NE, nbr_N, nbr_NW, nbr_SE
      public :: nbr_num
      public :: uarear_slot, dxN_slot, dyE_slot, uarea_slot
      public :: tarea_slot, earea_slot, narea_slot, hm_slot
      public :: dxU_slot, dyT_slot, dxT_slot, dyU_slot
      public :: dxE_slot, dyN_slot, epm_slot, npm_slot, uvm_slot
      public :: uareaavg_slot, wsumN_slot, wsumE_slot
      public :: wsumNW_slot, wsumSE_slot, wsumNE_slot
      ! NB the *_slot arrays already carry the public attribute on their own
      ! declarations; repeating them here is an ACCESS conflict (and the ifort
      ! #6251 class that bit f1aab62).  Only the procedures need listing.
      public :: grid_average_X2YA_compressed, grid_average_X2YS_compressed
      public :: set_evp_parameters, stepu, stepuv_CD, stepu_C, stepv_C, &
                principal_stress, init_dyn_shared, dyn_prep1, dyn_prep2, dyn_finish, &
                seabed_stress_factor_LKD, seabed_stress_factor_prob, &
                deformations, deformationsC_T, deformationsCD_T, &
                strain_rates, strain_rates_T, strain_rates_U, &
                strain_rates_Tdt_compressed, strain_rates_Tdtsd_compressed, &
                visc_replpress, &
                dyn_haloUpdate, &
                pack_to_slots, unpack_from_slots, &
                stack_fields, unstack_fields

      ! namelist parameters

      integer (kind=int_kind), public :: &
         kdyn       , & ! type of dynamics ( -1, 0 = off, 1 = evp, 2 = eap , 3 = vp)
         kridge     , & ! set to "-1" to turn off ridging
         ndte           ! number of subcycles

      character (len=char_len), public :: &
         coriolis   , & ! 'constant', 'zero', or 'latitude'
         ssh_stress     ! 'geostrophic' or 'coupled'

      logical (kind=log_kind), public :: &
         revised_evp    ! if true, use revised evp procedure

      character (len=char_len), public :: &
         evp_algorithm  ! standard_2d = 2D org version (standard)
                        ! shared_mem_1d = 1d without mpi call and refactorization to 1d

      real (kind=dbl_kind), public :: &
         dyn_area_min,& ! minimum ice area concentration to activate dynamics
         dyn_mass_min,& ! minimum ice mass to activate dynamics (kg/m^2)
         elasticDamp    ! coefficient for calculating the parameter E, elastic damping parameter

      ! other dynamics parameters

      character (len=char_len), public :: &
         yield_curve      , & ! 'ellipse' ('teardrop' needs further testing)
         visc_method      , & ! method for viscosity calc at U points (C, CD grids)
         seabed_stress_method ! method for seabed stress calculation
                              ! LKD: Lemieux et al. 2015, probabilistic: Dupont et al. 2022

      real (kind=dbl_kind), parameter, public :: &
         rheo_area_min = 1e-3_dbl_kind, & ! minimum ice area concentration to activate rheology
         u0    = 5e-5_dbl_kind        , & ! residual velocity for seabed stress (m/s)
         cosw  = c1                   , & ! cos(ocean turning angle)  ! turning angle = 0
         sinw  = c0                       ! sin(ocean turning angle)  ! turning angle = 0

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
         dtei        , & ! ndte/dt, where dt/ndte is subcycling timestep (1/s)
         denom1          ! constants for stress equation

      real (kind=dbl_kind), public :: & ! Bouillon et al relaxation constants
         arlx        , & ! alpha for stressp
         arlx1i      , & ! (inverse of alpha) for stressp
         brlx            ! beta   for momentum

      real (kind=dbl_kind), allocatable, public :: &
         fcor_blk(:,:)       ! Coriolis parameter (1/s)

      real (kind=dbl_kind), allocatable, public :: &
         fcorE_blk(:,:)  , & ! Coriolis parameter at E points (1/s)
         fcorN_blk(:,:)      ! Coriolis parameter at N points  (1/s)

      real (kind=dbl_kind), allocatable, public :: &
         fld2(:,:,:,:), & ! 2 bundled fields
         fld3(:,:,:,:), & ! 3 bundled fields
         fld4(:,:,:,:)    ! 4 bundled fields

      ! bundle for COMPRESSED fields: (slot, field), slot 0 being each block's
      ! dead slot.  The compressed counterpart of fld2, so that the velocities
      ! cross the halo in one message instead of one each.
      real (kind=dbl_kind), allocatable, public :: &
         fld2c(:,:)   , & ! 2 bundled compressed fields
         fld3c(:,:)   , & ! 3 bundled compressed fields
         fld4c(:,:)       ! 4 bundled compressed fields

      ! ---------------------------------------------------------------------
      ! Compressed (ragged) storage on the static mask.
      !
      ! ONE list serves every grid location.  The static masks are set once in
      ! makemask and never change, so a cell's slot is fixed for the whole run:
      !
      !   dynCnt(iblk)    slots held by block iblk
      !   dynOff(iblk)    first slot of block iblk, minus one
      !   dynMap(ij,iblk) block-local slot of cell ij, zero if not stored
      !   dynTot          total slots over all local blocks
      !
      ! so cell (i,j,iblk) lives at dynOff(iblk) + dynMap((j-1)*nx_block+i,iblk),
      ! whether it is being used as a T, U, N or E quantity.
      !
      ! The list is built from the union of the masks the run needs, the way
      ! ice_dyn_evp1d builds its own from Tmask (or Tmask .or. Umask): a single
      ! list is what lets a U-point kernel read a T-point array, which is
      ! exactly what stepu does when it sums str over the four surrounding T
      ! cells.  Separate per-location lists cannot express that without
      ! cross-grid index tables.
      !
      ! Activity is NOT carried by the ordering but by the skip flags below,
      ! one per grid location over the same slots -- again as in evp1d, whose
      ! skipTcell/skipUcell both run over its single 1..nActive list.
      !
      ! The GHOST RING IS INCLUDED: stencil kernels and halo exchanges need the
      ! ghost cells to be addressable, and makemask sets the masks over the
      ! same range.
      ! ---------------------------------------------------------------------

      ! flattened: element (i,j,iblk) lives at ((j-1)*nx_block + i, iblk)

                            ! per-block scratch indexed by block-local slot

      ! set each timestep from the dynamic ice masks: .true. means the slot
      ! carries no ice of that kind this step, so the kernels leave it alone
      logical (kind=log_kind), dimension (:), allocatable, public :: &
         skipT, skipU, skipN, skipE

      ! ---------------------------------------------------------------------
      ! Neighbour slots.
      !
      ! Once compressed a neighbour is no longer at ijb+-1 or ijb+-nx_block,
      ! so stencil kernels need its slot looked up.  Xnbr(is,d) holds the slot
      ! of the cell one step in direction d from slot is, or 0 where that
      ! neighbour falls outside the block or is not in the static mask.
      ! Because the layout is static these are built once, at allocation.
      !
      ! These are the six directions the 1d kernels use.  ice_dyn_evp1d holds
      ! the same set under historical names whose spelling does not match the
      ! offset; the correspondence is
      !
      !     here     offset     ice_dyn_evp1d
      !     nbr_W    (-1, 0)    ee
      !     nbr_SW   (-1,-1)    ne
      !     nbr_S    ( 0,-1)    se
      !     nbr_E    (+1, 0)    nw
      !     nbr_NE   (+1,+1)    sw
      !     nbr_N    ( 0,+1)    sse
      !
      ! A missing neighbour is the block's DEAD SLOT, the last element of the
      ! block's span (cnt+1), held at zero -- which is what the block layout
      ! carried at land points.  So a kernel may read it without a branch.
      ! The slices are 1-BASED: ice_dyn_core1d's kernels take assumed-shape
      ! dummies, which rebase to 1 whatever the actual's bounds, so a 0-based
      ! slice would silently shift every neighbour index by one.
      ! ---------------------------------------------------------------------
      ! These are COMPASS directions.  Do not confuse nbr_NW below with the
      ! core1d dummy named nw, which is this table's nbr_E -- core1d names its
      ! neighbour arguments after the stress corner they feed, not the compass.
      !
      ! NW and SE exist for grid_average_X2Y on the C and CD grids: the E2NA
      ! shift averages over (0,0), N, W, NW and N2EA over (0,0), S, E, SE.
      ! Nothing on the B grid reads them, so they cost one unused column pair
      ! in dynNbr there -- two integers per slot, a few MB on a large grid,
      ! against the complexity of making nbr_num depend on grid_ice.

      ! shape (dynTot, nbr_num), so dynNbr(:,d) is contiguous and one
      ! direction can be passed as a contiguous rank-1 argument

      ! Reverse of dynMap: dynIJ(slot) is the cell (j-1)*nx_block+i that owns
      ! that slot.  Kernels loop over slots and need the block index only for
      ! the few block-layout outputs they still write.

      ! Compressed working copies of the velocities; each block's span runs
      ! 1..cnt with its dead slot last, so it is addressable without a branch.
      ! uvel/vvel themselves stay in
      ! the block layout: they are state, read by restart, history,
      ! diagnostics and the couplers.  These are packed once per timestep,
      ! subcycled on, and unpacked again.
      real (kind=dbl_kind), dimension (:), allocatable, public :: &
         uvel_slot       , & ! x-component of velocity (m/s), compressed
         vvel_slot           ! y-component of velocity (m/s), compressed

      ! Compressed working copies of the U-point coefficients stepu reads.
      ! Unlike the velocities these are CONSTANT through the subcycling --
      ! aiU comes from grid_average_X2Y, fmU from dyn_prep2 and TbU from the
      ! seabed stress, all before the loop -- so they are packed once and
      ! never unpacked, and need no halo.  They carry a dead slot like the
      ! velocities purely so that pack_to_slots and the per-block slicing are
      ! the same everywhere; no stencil reads them, so the dead
      ! slot is never actually addressed.  The block-layout
      ! originals stay untouched for the consumers that blocked compressing
      ! them outright: grid_average_X2Y and ice_step_mod read aiU over the
      ! whole block, and fmU feeds the strcorx/strcory history fields.
      real (kind=dbl_kind), dimension (:), allocatable, public :: &
         aiU_slot        , & ! ice fraction on U grid, compressed
         fmU_slot        , & ! Coriolis parameter * mass in U cell, compressed
         TbU_slot            ! seabed stress factor (N/m^2), compressed

      ! The dyn_prep2 work arrays, likewise compressed for stepu and likewise
      ! constant through the subcycling.  The block-layout originals must
      ! stay: they are handed WHOLE to dyn_evp1d_run, whose gather_dyn calls
      ! gather_global and needs the block layout.  These are additional
      ! copies, not a relayout of those arrays -- see the note on their
      ! dyn_prep2 declarations.
      ! uarear is a grid metric and never changes at all; it is repacked with
      ! the rest each timestep only to keep one code path, at the cost of one
      ! extra sweep per block.  uocnU/vocnU/cdn_ocnU are the first arrays
      ! compressed here that were still rank-3 rather than flat.
      real (kind=dbl_kind), dimension (:), allocatable, public :: &
         uocnU_slot      , & ! ocean current, x, compressed
         vocnU_slot      , & ! ocean current, y, compressed
         cdn_ocnU_slot       ! ocean-ice neutral drag coefficient, compressed

      ! C and CD grids.  Allocated only for those grids, beside skipN/skipE.
      !
      ! Which of these are compressed follows the same rule as the B grid:
      ! compress what the subcycling changes AND a stencil reads at a
      ! neighbour; everything else keeps block indexing through dynIJ(is).
      ! The one addition is a grid metric read at a NEIGHBOUR -- dxN, dyE and
      ! uarea -- because reaching those block-wise would need
      ! dynIJ(nbr(is,d)), a second indirection per access.  Those are constant,
      ! so they are packed once and never unpacked, like uarear_slot.
      !
      ! Stencils, all inside the six directions dynNbr already holds:
      !   strain_rates_T  S on uvelN,dxN   W on vvelE,dyE
      !   stressC_T       W, S, SW on shearU, uarea
      !   div_stress_E*   S, E             div_stress_N*  W, N
      !
      ! Copy, never relayout: the block-layout originals stay written in full.
      real (kind=dbl_kind), dimension (:), allocatable, public :: &
         uvelE_slot      , & ! x-component of E-point velocity, compressed
         vvelE_slot      , & ! y-component of E-point velocity, compressed
         uvelN_slot      , & ! x-component of N-point velocity, compressed
         vvelN_slot      , & ! y-component of N-point velocity, compressed
         shearU_slot     , & ! shear at the U point, compressed
         ! The T-point stresses stressC_T writes and div_stress_* reads one
         ! cell away (S/E for the E kernels, W/N for the N ones), plus the
         ! U-point shear stress stressC_U writes and the same kernels read.
         ! zetax2T, etax2T, divergU, tensionU, deltaU, strengthU, etax2U,
         ! strintxE and strintyN are all read POINTWISE, so they stay block
         ! indexed via dynIJ(is) and are deliberately absent here.
         stresspT_slot   , & ! sigma11+sigma22 at T, compressed
         stressmT_slot   , & ! sigma11-sigma22 at T, compressed
         stress12T_slot  , & ! sigma12 at T, compressed
         stress12U_slot  , & ! sigma12 at U, compressed
         ! CD only: div_stress_Ey and div_stress_Nx read these at a neighbour
         stresspU_slot   , & ! sigma11+sigma22 at U, compressed
         stressmU_slot   , & ! sigma11-sigma22 at U, compressed
         ! The viscosities and strength.  NOTHING outside the dynamics reads
         ! zetax2T, etax2T, zetax2U, etax2U or strengthU -- checked across the
         ! tree -- so these need no block form at all; they are compressed
         ! because grid_average_X2YS_compressed shifts T to U on them.
         zetax2T_slot    , & ! 2*zeta at T, compressed
         etax2T_slot     , & ! 2*eta  at T, compressed
         zetax2U_slot    , & ! 2*zeta at U, compressed
         etax2U_slot     , & ! 2*eta  at U, compressed
         strength_slot   , & ! ice strength at T, compressed
         strengthU_slot  , & ! ice strength at U, compressed
         ! div_stress_* writes these and the step kernels read them in the
         ! same pass, so they are pure subcycle state.  The block arrays of
         ! the same name are diagnostics -- history reads them through
         ! f_strintx/f_strinty -- and are written ONCE after the subcycle,
         ! the same treatment 136e861 gave the B grid's strintxU/strintyU.
         strintxE_slot   , & ! div of internal stress, x at E, compressed
         strintyE_slot   , & ! div of internal stress, y at E, compressed
         strintxN_slot   , & ! div of internal stress, x at N, compressed
         strintyN_slot   , & ! div of internal stress, y at N, compressed
         ! Seabed stress, written by the step kernels purely for output.
         ! Same treatment as strint*: slot while subcycling, unpacked once.
         taubxE_slot     , & ! seabed stress, x at E, compressed
         taubyE_slot     , & ! seabed stress, y at E, compressed
         taubxN_slot     , & ! seabed stress, x at N, compressed
         taubyN_slot     , & ! seabed stress, y at N, compressed
         ! The per-cell coefficients the step kernels read.  All are constant
         ! through the subcycling -- dyn_prep1/dyn_prep2 set them before the
         ! loop -- so the driver packs them once.  With these compressed there
         ! is no rank-2 read left in stepu_C/stepv_C/stepuv_CD, which is what
         ! lets their cell loop run on slots like the B grid's stepu.
         aiE_slot           , & ! ice fraction
         cdn_ocnE_slot      , & ! ocean-ice neutral drag coefficient
         uocnE_slot         , & ! ocean current, x (m/s)
         vocnE_slot         , & ! ocean current, y (m/s)
         waterxE_slot       , & ! for ocean stress, x (m/s)
         wateryE_slot       , & ! for ocean stress, y (m/s)
         forcexE_slot       , & ! atm stress + ocn tilt, x
         forceyE_slot       , & ! atm stress + ocn tilt, y
         emassdti_slot      , & ! mass of cell/dt (kg/m^2 s)
         fmE_slot           , & ! Coriolis param. * mass (kg/s)
         TbE_slot           , & ! seabed stress factor (N/m^2)
         aiN_slot           , & ! ice fraction
         cdn_ocnN_slot      , & ! ocean-ice neutral drag coefficient
         uocnN_slot         , & ! ocean current, x (m/s)
         vocnN_slot         , & ! ocean current, y (m/s)
         waterxN_slot       , & ! for ocean stress, x (m/s)
         wateryN_slot       , & ! for ocean stress, y (m/s)
         forcexN_slot       , & ! atm stress + ocn tilt, x
         forceyN_slot       , & ! atm stress + ocn tilt, y
         nmassdti_slot      , & ! mass of cell/dt (kg/m^2 s)
         fmN_slot           , & ! Coriolis param. * mass (kg/s)
         TbN_slot           , & ! seabed stress factor (N/m^2)
         strintxU_slot   , & ! div of internal stress, x at U, compressed
         strintyU_slot   , & ! div of internal stress, y at U, compressed
         taubxU_slot     , & ! seabed stress, x at U, compressed
         taubyU_slot         ! seabed stress, y at U, compressed

      real (kind=dbl_kind), dimension (:), allocatable, public :: &
         umassdti_slot   , & ! mass of U-cell/dt (kg/m^2 s), compressed
         waterxU_slot    , & ! ocean stress coefficient, x, compressed
         wateryU_slot    , & ! ocean stress coefficient, y, compressed
         forcexU_slot    , & ! atm stress + ocn tilt, x, compressed
         forceyU_slot        ! atm stress + ocn tilt, y, compressed

      real (kind=dbl_kind), dimension (:), allocatable, public :: &
         uvel_init       , & ! x-component of velocity (m/s), beginning of timestep
         vvel_init           ! y-component of velocity (m/s), beginning of timestep

      real (kind=dbl_kind), dimension (:), allocatable, public :: &
         uvelN_init      , & ! x-component of velocity (m/s), beginning of timestep
         vvelN_init          ! y-component of velocity (m/s), beginning of timestep

      real (kind=dbl_kind), dimension (:), allocatable, public :: &
         uvelE_init      , & ! x-component of velocity (m/s), beginning of timestep
         vvelE_init          ! y-component of velocity (m/s), beginning of timestep

      ! the i and j dimensions are flattened: element (i,j,iblk) of the block
      ! layout lives at ((j-1)*nx_block + i, iblk)
      logical (kind=log_kind), dimension (:,:), allocatable, public :: &
         iceUmask, &   ! ice extent mask (U-cell)
         iceNmask, &   ! ice extent mask (N-cell)
         iceEmask      ! ice extent mask (E-cell)

      ! the i and j dimensions are flattened: element (i,j,iblk) of the block
      ! layout lives at ((j-1)*nx_block + i, iblk)
      logical (kind=log_kind), dimension (:,:), allocatable, target, public :: &
         iceTmask      ! ice extent mask (T-cell)

      real (kind=dbl_kind), allocatable, public :: &
         DminTarea(:,:)      ! deltamin * tarea (m^2/s)

      ! the i and j dimensions are flattened: element (i,j,iblk) of the block
      ! layout lives at ((j-1)*nx_block + i, iblk)
      real (kind=dbl_kind), dimension (:,:), allocatable, public :: &
         cyp    , & ! 1.5*HTE(i,j)-0.5*HTW(i,j) = 1.5*HTE(i,j)-0.5*HTE(i-1,j)
         cxp    , & ! 1.5*HTN(i,j)-0.5*HTS(i,j) = 1.5*HTN(i,j)-0.5*HTN(i,j-1)
         cym    , & ! 0.5*HTE(i,j)-1.5*HTW(i,j) = 0.5*HTE(i,j)-1.5*HTE(i-1,j)
         cxm    , & ! 0.5*HTN(i,j)-1.5*HTS(i,j) = 0.5*HTN(i,j)-1.5*HTN(i,j-1)
         dxhy   , & ! 0.5*(HTE(i,j) - HTW(i,j)) = 0.5*(HTE(i,j) - HTE(i-1,j))
         dyhx       ! 0.5*(HTN(i,j) - HTS(i,j)) = 0.5*(HTN(i,j) - HTN(i,j-1))

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
         module procedure dyn_haloUpdate1_flat
         module procedure dyn_haloUpdate2_flat
         module procedure dyn_haloUpdate3_flat
         module procedure dyn_haloUpdate4_flat
      end interface

      interface stack_fields
         module procedure stack_fields2_compressed
         module procedure stack_fields3_compressed
         module procedure stack_fields4_compressed
         module procedure stack_fields2
         module procedure stack_fields3
         module procedure stack_fields3_flat
         module procedure stack_fields2_flat
         module procedure stack_fields4_flat
         module procedure stack_fields4
         module procedure stack_fields5
      end interface

      interface unstack_fields
         module procedure unstack_fields2_compressed
         module procedure unstack_fields3_compressed
         module procedure unstack_fields4_compressed
         module procedure unstack_fields2
         module procedure unstack_fields3
         module procedure unstack_fields3_flat
         module procedure unstack_fields2_flat
         module procedure unstack_fields4_flat
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

      use ice_grid, only: tmask, umask, umaskCD, nmask, emask

      integer (int_kind) :: ierr

      ! union of every static mask this run needs, so that one list can serve
      ! all grid locations

      character(len=*), parameter :: subname = '(alloc_dyn_shared)'

      ! The ragged layout itself is built by init_ragged_layout in ice_grid,
      ! at the end of init_grid2, because it is derived from the grid masks and
      ! because the grid metrics compressed once per run are packed there too.
      ! By the time any solver calls this, dynTot and dynMap are set.
      allocate(skipT(dynTot), skipU(dynTot), stat=ierr)
      if (ierr/=0) call abort_ice(subname//': Out of memory')
      skipT = .true.
      skipU = .true.

      allocate(uvel_slot(dynTot), vvel_slot(dynTot), stat=ierr)
      if (ierr/=0) call abort_ice(subname//': Out of memory')
      uvel_slot = c0
      vvel_slot = c0

      allocate(fld2c(dynTot,2), fld3c(dynTot,3), fld4c(dynTot,4), stat=ierr)
      if (ierr/=0) call abort_ice(subname//': Out of memory')
      fld2c = c0
      fld3c = c0
      fld4c = c0

      allocate(strintxU_slot(dynTot), strintyU_slot(dynTot), &
               taubxU_slot(dynTot), taubyU_slot(dynTot), stat=ierr)
      if (ierr/=0) call abort_ice(subname//': Out of memory')
      strintxU_slot = c0 ; strintyU_slot = c0
      taubxU_slot   = c0 ; taubyU_slot   = c0

      allocate(aiU_slot(dynTot), fmU_slot(dynTot), TbU_slot(dynTot), stat=ierr)
      if (ierr/=0) call abort_ice(subname//': Out of memory')
      aiU_slot = c0
      fmU_slot = c0
      TbU_slot = c0

      allocate(umassdti_slot(dynTot), waterxU_slot(dynTot), &
               wateryU_slot(dynTot), forcexU_slot(dynTot), &
               forceyU_slot(dynTot), stat=ierr)
      if (ierr/=0) call abort_ice(subname//': Out of memory')
      umassdti_slot = c0
      waterxU_slot  = c0
      wateryU_slot  = c0
      forcexU_slot  = c0
      forceyU_slot  = c0

      allocate(uocnU_slot(dynTot), vocnU_slot(dynTot), &
               cdn_ocnU_slot(dynTot), stat=ierr)
      if (ierr/=0) call abort_ice(subname//': Out of memory')
      uocnU_slot    = c0
      vocnU_slot    = c0
      cdn_ocnU_slot = c0

      allocate( &
         uvel_init (dynTot), & ! x-component of velocity (m/s), beginning of timestep
         vvel_init (dynTot), & ! y-component of velocity (m/s), beginning of timestep
         iceTmask  (nx_block*ny_block,max_blocks), & ! T mask for dynamics
         iceUmask  (nx_block*ny_block,max_blocks), & ! U mask for dynamics
         fcor_blk  (nx_block*ny_block,max_blocks), & ! Coriolis
         DminTarea (nx_block*ny_block,max_blocks), & !
         stat=ierr)
      if (ierr/=0) call abort_ice(subname//': Out of memory')

      uvel_init  = c0
      vvel_init  = c0
      iceTmask   = .false.
      iceUmask   = .false.
      fcor_blk   = c0
      DminTarea  = c0

      allocate( &
         fld2(nx_block,ny_block,2,max_blocks), &
         fld3(nx_block,ny_block,3,max_blocks), &
         fld4(nx_block,ny_block,4,max_blocks), &
         stat=ierr)
      if (ierr/=0) call abort_ice(subname//': Out of memory')

      fld2 = c0
      fld3 = c0
      fld4 = c0

      allocate( &
         cyp(nx_block*ny_block,max_blocks), & ! 1.5*HTE - 0.5*HTW
         cxp(nx_block*ny_block,max_blocks), & ! 1.5*HTN - 0.5*HTS
         cym(nx_block*ny_block,max_blocks), & ! 0.5*HTE - 1.5*HTW
         cxm(nx_block*ny_block,max_blocks), & ! 0.5*HTN - 1.5*HTS
         stat=ierr)
      if (ierr/=0) call abort_ice(subname//': Out of memory')

      cyp = c0
      cxp = c0
      cym = c0
      cxm = c0

      if (grid_ice == 'B' .and. evp_algorithm == "standard_2d") then
         allocate( &
            dxhy(nx_block*ny_block,max_blocks), & ! 0.5*(HTE - HTW)
            dyhx(nx_block*ny_block,max_blocks), & ! 0.5*(HTN - HTS)
            stat=ierr)
         if (ierr/=0) call abort_ice(subname//': Out of memory')
         dxhy = c0
         dyhx = c0
      endif

      if (grid_ice == 'CD' .or. grid_ice == 'C') then
         ! N and E share the one list; only their activity flags differ
         allocate(skipN(dynTot), skipE(dynTot), stat=ierr)
         if (ierr/=0) call abort_ice(subname//': Out of memory')
         skipN = .true.
         skipE = .true.

         allocate(uvelE_slot(dynTot), vvelE_slot(dynTot), &
                  uvelN_slot(dynTot), vvelN_slot(dynTot), &
                  shearU_slot(dynTot),                    &
                  stresspT_slot(dynTot), stressmT_slot(dynTot),  &
                  stress12T_slot(dynTot), stress12U_slot(dynTot), &
                  stresspU_slot(dynTot), stressmU_slot(dynTot), &
                  zetax2T_slot(dynTot), etax2T_slot(dynTot), &
                  zetax2U_slot(dynTot), etax2U_slot(dynTot), &
                  strength_slot(dynTot), strengthU_slot(dynTot), &
                  strintxE_slot(dynTot), strintyE_slot(dynTot), &
                  strintxN_slot(dynTot), strintyN_slot(dynTot), &
                  taubxE_slot(dynTot), taubyE_slot(dynTot), &
                  taubxN_slot(dynTot), taubyN_slot(dynTot), &
                  aiE_slot(dynTot), &
                  cdn_ocnE_slot(dynTot), &
                  uocnE_slot(dynTot), &
                  vocnE_slot(dynTot), &
                  waterxE_slot(dynTot), &
                  wateryE_slot(dynTot), &
                  forcexE_slot(dynTot), &
                  forceyE_slot(dynTot), &
                  emassdti_slot(dynTot), &
                  fmE_slot(dynTot), &
                  TbE_slot(dynTot), &
                  aiN_slot(dynTot), &
                  cdn_ocnN_slot(dynTot), &
                  uocnN_slot(dynTot), &
                  vocnN_slot(dynTot), &
                  waterxN_slot(dynTot), &
                  wateryN_slot(dynTot), &
                  forcexN_slot(dynTot), &
                  forceyN_slot(dynTot), &
                  nmassdti_slot(dynTot), &
                  fmN_slot(dynTot), &
                  TbN_slot(dynTot), &
                  stat=ierr)
         if (ierr/=0) call abort_ice(subname//': Out of memory')
         uvelE_slot  = c0
         vvelE_slot  = c0
         uvelN_slot  = c0
         vvelN_slot  = c0
         shearU_slot = c0
         stresspT_slot  = c0
         stressmT_slot  = c0
         stress12T_slot = c0
         stress12U_slot = c0
         stresspU_slot  = c0
         stressmU_slot  = c0
         zetax2T_slot   = c0 ; etax2T_slot    = c0
         zetax2U_slot   = c0 ; etax2U_slot    = c0
         strength_slot  = c0 ; strengthU_slot = c0
         strintxE_slot  = c0 ; strintyE_slot  = c0
         strintxN_slot  = c0 ; strintyN_slot  = c0
         taubxE_slot    = c0 ; taubyE_slot    = c0
         taubxN_slot    = c0 ; taubyN_slot    = c0
         aiE_slot = c0
         cdn_ocnE_slot = c0
         uocnE_slot = c0
         vocnE_slot = c0
         waterxE_slot = c0
         wateryE_slot = c0
         forcexE_slot = c0
         forceyE_slot = c0
         emassdti_slot = c0
         fmE_slot = c0
         TbE_slot = c0
         aiN_slot = c0
         cdn_ocnN_slot = c0
         uocnN_slot = c0
         vocnN_slot = c0
         waterxN_slot = c0
         wateryN_slot = c0
         forcexN_slot = c0
         forceyN_slot = c0
         nmassdti_slot = c0
         fmN_slot = c0
         TbN_slot = c0

         allocate( &
            uvelE_init (dynTot), & ! x-component of velocity (m/s), beginning of timestep
            vvelE_init (dynTot), & ! y-component of velocity (m/s), beginning of timestep
            uvelN_init (dynTot), & ! x-component of velocity (m/s), beginning of timestep
            vvelN_init (dynTot), & ! y-component of velocity (m/s), beginning of timestep
            iceEmask   (nx_block*ny_block,max_blocks), & ! T mask for dynamics
            iceNmask   (nx_block*ny_block,max_blocks), & ! U mask for dynamics
            fcorE_blk  (nx_block*ny_block,max_blocks), & ! Coriolis
            fcorN_blk  (nx_block*ny_block,max_blocks), &   ! Coriolis
            stat=ierr)
         if (ierr/=0) call abort_ice(subname//': Out of memory')
         uvelE_init  = c0
         vvelE_init  = c0
         uvelN_init  = c0
         vvelN_init  = c0
         iceEmask    = .false.
         iceNmask    = .false.
         fcorE_blk   = c0
         fcorN_blk   = c0
      endif

      end subroutine alloc_dyn_shared






!=======================================================================
! Compressed counterparts of grid_average_X2Y, for the C and CD subcycle.
!
! ADDITIVE, not a conversion.  grid_average_X2Y has some 200 call sites
! across the model -- ice_step_mod, ice_diagnostics, ice_forcing, vp, and the
! gridavgchk unit test -- so the block-layout family stays exactly as it is.
! These cover only the five shifts the C/CD subcycle performs:
!
!   T2US  NE   strength, etax2T -> U     masked   (X2YS)
!   E2UA  N    uvelE -> uvel             area     (X2YA)
!   N2UA  E    vvelN -> vvel             area     (X2YA)
!   E2NA  NW   uvelE -> uvelN            area     (X2YA)
!   N2EA  SE   vvelN -> vvelE            area     (X2YA)
!
! They live here rather than in ice_grid because they need dynMap, dynOff and
! dynNbr, and ice_dyn_shared already uses ice_grid -- the other way round
! would be circular.
!
! The block loop is kept rather than iterating slots: the originals run over
! the PHYSICAL domain only, ilo:ihi by jlo:jhi, and the compressed layout also
! holds ghost cells, so a plain slot sweep would cover cells the original never
! touches.  This is the same shape stress uses: keep the loop, derive the slot.
!
! SUMMATION ORDER IS PRESERVED EXACTLY.  It differs per direction -- NW sums
! W, self, NW, N; SE sums S, SE, self, E -- and reordering would cost
! bit-for-bit agreement with the block-layout routine.
!
! dynNbr holds BLOCK-LOCAL slots, so a global index is dynOff(iblk) + nbr.

      subroutine grid_average_X2YA_compressed(dir, work1, wght1, wsum, work2)

      use ice_blocks, only: block, get_block
      use ice_domain, only: nblocks, blocks_ice

      character(len=*), intent(in) :: &
         dir            ! shift direction: N, E, NW or SE

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         work1      , & ! field to shift, compressed
         wght1      , & ! weights, compressed
         wsum           ! precomputed weight sum -- see ice_grid

      real (kind=dbl_kind), dimension(:), intent(out) :: &
         work2          ! shifted field, compressed

      integer (kind=int_kind) :: &
         i, j, iblk, ijb, m, is, isN, isE, isW, isS, isNW, isSE, off
      integer (kind=int_kind) :: ilo,ihi,jlo,jhi
      type (block) :: this_block
      real (kind=dbl_kind) :: wtmp

      character(len=*), parameter :: subname = '(grid_average_X2YA_compressed)'

      work2(:) = c0

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ijb,m,is,isN,isE,isW,isS,isNW,isSE, &
      !$OMP                     off,ilo,ihi,jlo,jhi,this_block,wtmp) SCHEDULE(runtime)
      do iblk = 1, nblocks
         if (dynCnt(iblk) == 0) cycle
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi
         off = dynOff(iblk)

         do j = jlo, jhi
         do i = ilo, ihi
            ijb = (j-1)*nx_block + i
            m   = dynMap(ijb,iblk)
            if (m == 0) cycle
            is  = off + m

            select case (trim(dir))

            case('N')
               isN = off + dynNbr(is,nbr_N)
               wtmp = wsum(is)
               if (wtmp /= c0) &
               work2(is) = (work1(is )*wght1(is )  &
                          + work1(isN)*wght1(isN)) &
                          / wtmp

            case('E')
               isE = off + dynNbr(is,nbr_E)
               wtmp = wsum(is)
               if (wtmp /= c0) &
               work2(is) = (work1(is )*wght1(is )  &
                          + work1(isE)*wght1(isE)) &
                          / wtmp

            case('NW')
               isW  = off + dynNbr(is,nbr_W )
               isNW = off + dynNbr(is,nbr_NW)
               isN  = off + dynNbr(is,nbr_N )
               wtmp = wsum(is)
               if (wtmp /= c0) &
               work2(is) = (work1(isW )*wght1(isW )  &
                          + work1(is  )*wght1(is  )  &
                          + work1(isNW)*wght1(isNW)  &
                          + work1(isN )*wght1(isN )) &
                          / wtmp

            case('SE')
               isS  = off + dynNbr(is,nbr_S )
               isSE = off + dynNbr(is,nbr_SE)
               isE  = off + dynNbr(is,nbr_E )
               wtmp = wsum(is)
               if (wtmp /= c0) &
               work2(is) = (work1(isS )*wght1(isS )  &
                          + work1(isSE)*wght1(isSE)  &
                          + work1(is  )*wght1(is  )  &
                          + work1(isE )*wght1(isE )) &
                          / wtmp

            case default
               call abort_ice(subname//' ERROR: unsupported dir='//trim(dir), &
                  file=__FILE__, line=__LINE__)
            end select

         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine grid_average_X2YA_compressed

!=======================================================================
! Masked counterpart, for the T2US shift.  Only NE is needed.

      subroutine grid_average_X2YS_compressed(dir, work1, wght1, mask1, wsum, work2)

      use ice_blocks, only: block, get_block
      use ice_domain, only: nblocks, blocks_ice

      character(len=*), intent(in) :: &
         dir            ! shift direction: NE

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         work1      , & ! field to shift, compressed
         wght1      , & ! weights, compressed
         mask1      , & ! mask, compressed
         wsum           ! precomputed weight sum -- see ice_grid

      real (kind=dbl_kind), dimension(:), intent(out) :: &
         work2          ! shifted field, compressed

      integer (kind=int_kind) :: &
         i, j, iblk, ijb, m, is, isE, isN, isNE, off
      integer (kind=int_kind) :: ilo,ihi,jlo,jhi
      type (block) :: this_block
      real (kind=dbl_kind) :: wtmp

      character(len=*), parameter :: subname = '(grid_average_X2YS_compressed)'

      work2(:) = c0

      if (trim(dir) /= 'NE') then
         call abort_ice(subname//' ERROR: unsupported dir='//trim(dir), &
            file=__FILE__, line=__LINE__)
      endif

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ijb,m,is,isE,isN,isNE,off, &
      !$OMP                     ilo,ihi,jlo,jhi,this_block,wtmp) SCHEDULE(runtime)
      do iblk = 1, nblocks
         if (dynCnt(iblk) == 0) cycle
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi
         off = dynOff(iblk)

         do j = jlo, jhi
         do i = ilo, ihi
            ijb = (j-1)*nx_block + i
            m   = dynMap(ijb,iblk)
            if (m == 0) cycle
            is   = off + m
            isE  = off + dynNbr(is,nbr_E )
            isN  = off + dynNbr(is,nbr_N )
            isNE = off + dynNbr(is,nbr_NE)

            wtmp = wsum(is)
            if (wtmp /= c0) &
            work2(is) = (mask1(is  )*work1(is  )*wght1(is  )  &
                       + mask1(isE )*work1(isE )*wght1(isE )  &
                       + mask1(isN )*work1(isN )*wght1(isN )  &
                       + mask1(isNE)*work1(isNE)*wght1(isNE)) &
                       / wtmp
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine grid_average_X2YS_compressed

!=======================================================================
! Initialize parameters and variables needed for the dynamics
! author: Elizabeth C. Hunke, LANL

      subroutine init_dyn_shared (dt)

      use ice_blocks, only: block, get_block, ns_boundary_type, nghost
      use ice_boundary, only: ice_halo, ice_haloUpdate
      use ice_domain, only: nblocks, halo_dynbundle, blocks_ice, halo_info
      use ice_flux, only: &
          stressp_1, stressp_2, stressp_3, stressp_4, &
          stressm_1, stressm_2, stressm_3, stressm_4, &
          stress12_1, stress12_2, stress12_3, stress12_4, &
          stresspT, stressmT, stress12T, &
          stresspU, stressmU, stress12U
      use ice_state, only: uvel, vvel, uvelE, vvelE, uvelN, vvelN
      use ice_grid, only: ULAT, NLAT, ELAT, tarea, HTE, HTN
      use ice_constants, only: field_loc_center, field_type_vector

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) :: &
         i, j  ,             & ! indices
         ij    ,             & ! 1D index of (i,j) in the block
         ijb   ,             & ! 1D index of (i,j) in the flattened block layout
         ilo, ihi, jlo, jhi, & !min and max index for interior of blocks
         nprocs,             & ! number of processors
         iblk                  ! block index

      type (block) :: &
         this_block   ! block information for current block

      character(len=*), parameter :: subname = '(init_dyn_shared)'

      ! checks
      if (kdyn == 1 .and. evp_algorithm == 'shared_mem_1d' .and. &
          (ns_boundary_type == 'tripole' .or. ns_boundary_type == 'tripoleT')) then
         call abort_ice(subname//' ERROR: evp_alg shared mem 1d not supported with tripole', &
            file=__FILE__, line=__LINE__)
      endif

      ! calc_halo_parent builds the 1d halo assuming a single ghost layer: it
      ! pairs i==nx-nghost with ifind=1 and i==1+nghost with ifind=nx, which
      ! only spans one layer.  With nghost > 1 the outer layers would simply
      ! never be filled -- wrong answers, not a crash -- so refuse rather than
      ! run.  nghost is a parameter today, so this is a guard against someone
      ! raising it, matching the TODO in calc_halo_parent.
      if (kdyn == 1 .and. evp_algorithm == 'shared_mem_1d' .and. nghost > 1) then
         call abort_ice(subname//' ERROR: evp_alg shared mem 1d not supported with nghost > 1', &
            file=__FILE__, line=__LINE__)
      endif

      call set_evp_parameters (dt)
      ! allocate dyn shared (init_uvel,init_vvel)
      call alloc_dyn_shared
      ! Set halo_dynbundle, this is empirical at this point, could become namelist
      halo_dynbundle = .true.
      nprocs = get_num_procs()
      if (nx_block*ny_block/nprocs > 100) halo_dynbundle = .false.

      if (my_task == master_task) then
         write(nu_diag,*) 'dt  = ',dt
         write(nu_diag,*) 'dt_subcyle = ',dt/real(ndte,kind=dbl_kind)
         write(nu_diag,*) 'tdamp =', elasticDamp * dt
         write(nu_diag,*) 'halo_dynbundle =', halo_dynbundle
      endif

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ij) SCHEDULE(runtime)
      do iblk = 1, nblocks
      do j = 1, ny_block
      do i = 1, nx_block

         ij = (j-1)*nx_block + i

         ! velocity
         uvel(i,j,iblk) = c0    ! m/s
         vvel(i,j,iblk) = c0    ! m/s
         if (grid_ice == 'CD' .or. grid_ice == 'C') then ! extra velocity variables
            uvelE(i,j,iblk) = c0
            vvelE(i,j,iblk) = c0
            uvelN(i,j,iblk) = c0
            vvelN(i,j,iblk) = c0
         endif


         ! Coriolis parameter
         if (trim(coriolis) == 'constant') then
            fcor_blk(ij,iblk) = 1.46e-4_dbl_kind ! Hibler 1979, N. Hem; 1/s
         else if (trim(coriolis) == 'zero') then
            fcor_blk(ij,iblk) = c0
         else
            fcor_blk(ij,iblk) = c2*omega*sin(ULAT(ij,iblk)) ! 1/s
         endif

         if (grid_ice == 'CD' .or. grid_ice == 'C') then

            if (trim(coriolis) == 'constant') then
               fcorE_blk(ij,iblk) = 1.46e-4_dbl_kind ! Hibler 1979, N. Hem; 1/s
               fcorN_blk(ij,iblk) = 1.46e-4_dbl_kind ! Hibler 1979, N. Hem; 1/s
            else if (trim(coriolis) == 'zero') then
               fcorE_blk(ij,iblk) = c0
               fcorN_blk(ij,iblk) = c0
            else
               fcorE_blk(ij,iblk) = c2*omega*sin(ELAT(ij,iblk)) ! 1/s
               fcorN_blk(ij,iblk) = c2*omega*sin(NLAT(ij,iblk)) ! 1/s
            endif

         endif

         ! stress tensor,  kg/s^2
         stressp_1 (ij,iblk) = c0
         stressp_2 (ij,iblk) = c0
         stressp_3 (ij,iblk) = c0
         stressp_4 (ij,iblk) = c0
         stressm_1 (ij,iblk) = c0
         stressm_2 (ij,iblk) = c0
         stressm_3 (ij,iblk) = c0
         stressm_4 (ij,iblk) = c0
         stress12_1(ij,iblk) = c0
         stress12_2(ij,iblk) = c0
         stress12_3(ij,iblk) = c0
         stress12_4(ij,iblk) = c0

         if (grid_ice == 'CD' .or. grid_ice == 'C') then
            stresspT  (ij,iblk) = c0
            stressmT  (ij,iblk) = c0
            stress12T (ij,iblk) = c0
            stresspU  (ij,iblk) = c0
            stressmU  (ij,iblk) = c0
            stress12U (ij,iblk) = c0
         endif

         if (kdyn == 1) then
            DminTarea(ij,iblk) = deltaminEVP*tarea((j-1)*nx_block+i,iblk)
         elseif (kdyn == 3) then
            DminTarea(ij,iblk) = deltaminVP*tarea((j-1)*nx_block+i,iblk)
         endif

         ! ice extent mask on velocity points
         iceUmask(ij,iblk) = .false.
         if (grid_ice == 'CD' .or. grid_ice == 'C') then
            iceEmask(ij,iblk) = .false.
            iceNmask(ij,iblk) = .false.
         end if
      enddo                     ! i
      enddo                     ! j
      enddo                     ! iblk
      !$OMP END PARALLEL DO

      if (grid_ice == 'B' .and. evp_algorithm == "standard_2d") then
         do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            do j = jlo, jhi
            do i = ilo, ihi
               dxhy((j-1)*nx_block+i,iblk) = p5*(HTE((j-1)*nx_block+i,iblk) - HTE((j-1)*nx_block+i-1,iblk))
               dyhx((j-1)*nx_block+i,iblk) = p5*(HTN((j-1)*nx_block+i,iblk) - HTN((j-2)*nx_block+i,iblk))
            enddo
            enddo
         enddo

         call ice_HaloUpdate (dxhy,               halo_info, &
                              field_loc_center,   field_type_vector, &
                              fillValue=c1)
         call ice_HaloUpdate (dyhx,               halo_info, &
                              field_loc_center,   field_type_vector, &
                              fillValue=c1)

     endif

     do iblk = 1, nblocks
        this_block = get_block(blocks_ice(iblk),iblk)
        ilo = this_block%ilo
        ihi = this_block%ihi
        jlo = this_block%jlo
        jhi = this_block%jhi
        do j = jlo, jhi+1
        do i = ilo, ihi+1
           ijb = (j-1)*nx_block + i
           cyp(ijb,iblk) = (c1p5*HTE(ijb,iblk) - p5*HTE(ijb-1,iblk))
           cxp(ijb,iblk) = (c1p5*HTN(ijb,iblk) - p5*HTN(ijb-nx_block,iblk))
            ! match order of operations in cyp, cxp for tripole grids
           cym(ijb,iblk) = -(c1p5*HTE(ijb-1,iblk) - p5*HTE(ijb,iblk))
           cxm(ijb,iblk) = -(c1p5*HTN(ijb-nx_block,iblk) - p5*HTN(ijb,iblk))
        enddo
        enddo
     enddo

  end subroutine init_dyn_shared

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

      character(len=*), parameter :: subname = '(set_evp_parameters)'

      ! elastic time step
      dtei = real(ndte,kind=dbl_kind)/dt

      ! variables for elliptical yield curve and plastic potential
      epp2i = c1/e_plasticpot**2
      e_factor = e_yieldcurve**2 / e_plasticpot**4
      ecci = c1/e_yieldcurve**2 ! temporary for 1d evp

      if (revised_evp) then       ! Bouillon et al, Ocean Mod 2013
         revp   = c1
         denom1 = c1
         arlx1i = c1/arlx
      else                        ! Hunke, JCP 2013 with modified stress eq
         revp   = c0
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
                            vsno,      Tmask,    &
                            Tmass,     iceTmask)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ilo,ihi,jlo,jhi       ! beginning and end of physical domain

      ! all array arguments are assumed shape, so the actual arguments must be
      ! array sections of the right rank: (nx_block,ny_block) for the rank-2
      ! inputs and (nx_block*ny_block) for the flattened arguments, the latter
      ! in Fortran storage order so that element (i,j) is stored at
      ! (j-1)*nx_block+i

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         aice    , & ! concentration of ice
         vice    , & ! volume per unit area of ice          (m)
         vsno        ! volume per unit area of snow         (m)

      logical (kind=log_kind), dimension (:), intent(in) :: &
         Tmask       ! land/boundary mask, thickness (T-cell)

      real (kind=dbl_kind), dimension (:), intent(out) :: &
         Tmass       ! total mass of ice and snow (kg/m^2)

      logical (kind=log_kind), dimension (:), intent(out) :: &
         iceTmask    ! ice extent mask (T-cell)

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij       , & ! i, j and the matching 1D index into the block
         ijsub          , & ! 1D index over the physical domain only
         nisub, njsub       ! i and j extent of the physical domain

      real (kind=dbl_kind) :: &
         rhoi, rhos

      logical (kind=log_kind), dimension(nx_block*ny_block) :: &
         tmphm               ! temporary mask

      character(len=*), parameter :: subname = '(dyn_prep1)'

      call icepack_query_parameters(rhos_out=rhos, rhoi_out=rhoi)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      ! whole block, ghost cells included
      do ij = 1, nx_block*ny_block

         j = (ij-1)/nx_block + 1
         i = ij - (j-1)*nx_block

         !-----------------------------------------------------------------
         ! total mass of ice and snow, centered in T-cell
         ! NOTE: vice and vsno must be up to date in all grid cells,
         !       including ghost cells
         !-----------------------------------------------------------------
         if (Tmask(ij)) then
            Tmass(ij) = (rhoi*vice(i,j) + rhos*vsno(i,j)) ! kg/m^2
         else
            Tmass(ij) = c0
         endif

         !-----------------------------------------------------------------
         ! ice extent mask (T-cells)
         !-----------------------------------------------------------------
         tmphm(ij) = Tmask(ij) .and. (aice (i,j) > dyn_area_min) &
                                .and. (Tmass(ij)  > dyn_mass_min)

         !-----------------------------------------------------------------
         ! augmented mask (land + open ocean)
         !-----------------------------------------------------------------
         iceTmask (ij) = .false.

      enddo

      ! physical domain only: ilo:ihi within 1:nx_block, jlo:jhi within 1:ny_block,
      ! so ijsub runs over the subdomain and ij is its index in the full block
      nisub = ihi - ilo + 1
      njsub = jhi - jlo + 1

      do ijsub = 1, nisub*njsub

         j  = (ijsub-1)/nisub + jlo
         i  = ijsub - (j-jlo)*nisub + ilo - 1
         ij = (j-1)*nx_block + i

         ! extend ice extent mask (T-cells) to points around pack
         if (tmphm(ij+nx_block-1) .or. tmphm(ij+nx_block) .or. tmphm(ij+nx_block+1) .or. &
             tmphm(ij-1)          .or. tmphm(ij)          .or. tmphm(ij+1)          .or. &
             tmphm(ij-nx_block-1) .or. tmphm(ij-nx_block) .or. tmphm(ij-nx_block+1) ) then
            iceTmask(ij) = .true.
         endif

         if (.not.Tmask(ij)) iceTmask(ij) = .false.

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
                            icellT,     icellX,     &
                            indxTi,     indxTj,     &
                            indxXi,     indxXj,     &
                            aiX,        Xmass,      &
                            Xmassdti,   fcor,       &
                            Xmask,                  &
                            uocn,       vocn,       &
                            strairx,    strairy,    &
                            ss_tltx,    ss_tlty,    &
                            iceTmask,   iceXmask,   &
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
                            Xmap,       skipX,      &
                            TbU,        rheofactX,  &
                            Xmassdti_slot, waterx_slot, &
                            watery_slot,   forcex_slot, &
                            forcey_slot)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ilo,ihi,jlo,jhi       ! beginning and end of physical domain

      integer (kind=int_kind), intent(out) :: &
         icellT  , & ! no. of cells where iceTmask = .true.
         icellX      ! no. of cells where iceXmask = .true.

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(out) :: &
         indxTi  , & ! compressed index in i-direction on T grid
         indxTj  , & ! compressed index in j-direction
         indxXi  , & ! compressed index in i-direction on X grid, grid depends on call
         indxXj      ! compressed index in j-direction

      ! Xmask is assumed shape, 1D over the block in Fortran storage order,
      ! i.e. element (i,j) is stored at (j-1)*nx_block+i
      logical (kind=log_kind), dimension (:), intent(in) :: &
         Xmask       ! land/boundary mask, thickness (X-grid-cell)

      ! iceTmask is assumed shape, 1D over the block in Fortran storage order,
      ! i.e. element (i,j) is stored at (j-1)*nx_block+i
      logical (kind=log_kind), dimension (:), intent(in) :: &
         iceTmask    ! ice extent mask (T-cell)

      ! iceXmask is assumed shape, 1D over the block in Fortran storage order,
      ! i.e. element (i,j) is stored at (j-1)*nx_block+i
      logical (kind=log_kind), dimension (:), intent(inout) :: &
         iceXmask    ! ice extent mask (X-grid-cell)

      ! Xmass and Xmassdti are assumed shape, 1D over the block in Fortran
      ! storage order, i.e. element (i,j) is stored at (j-1)*nx_block+i
      real (kind=dbl_kind), dimension (:), intent(in) :: &
         aiX     , & ! ice fraction on u-grid (X grid)
         Xmass   , & ! total mass of ice and snow (X grid)
         fcor        ! Coriolis parameter (1/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         strairx , & ! stress on ice by air, x-direction
         strairy , & ! stress on ice by air, y-direction
         uocn    , & ! ocean current, x-direction (m/s)
         vocn    , & ! ocean current, y-direction (m/s)
         ss_tltx , & ! sea surface slope, x-direction (m/m)
         ss_tlty     ! sea surface slope, y-direction

      real (kind=dbl_kind), intent(in) :: &
         dt          ! time step

      ! uvel_init and vvel_init are COMPRESSED: cell ij of this block lives at
      ! slot Xmap(ij), see the layout description in this module's header
      real (kind=dbl_kind), dimension (:), intent(out) :: &
         uvel_init,& ! x-component of velocity (m/s), beginning of time step
         vvel_init   ! y-component of velocity (m/s), beginning of time step

      ! Xmap is assumed shape, 1D over the block in Fortran storage order,
      ! i.e. element (i,j) is at (j-1)*nx_block+i, holding that cell's static
      ! compressed slot, or 0 if the cell is not in the static mask
      integer (kind=int_kind), dimension (:), intent(in) :: &
         Xmap        ! static compressed slot of each cell

      ! skipX runs over this block's compressed slots: .true. where the slot
      ! carries no ice this step, so the kernels leave it alone
      logical (kind=log_kind), dimension (:), intent(out) :: &
         skipX       ! inactive-slot flag

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out) :: &
         TbU         ! seabed stress factor (N/m^2)

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
         strinty     ! divergence of internal ice stress, y (N/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout), &
           optional :: &
         rheofactX   ! mult. factor = 1, set to 0 if aiU <= rheo_area_min

      ! These five are in the BLOCK layout -- 1D over the block in Fortran
      ! storage order, element (i,j) at (j-1)*nx_block+i -- and must stay that
      ! way.  They are handed WHOLE to dyn_evp1d_run, whose gather_dyn calls
      ! gather_global and so requires the block layout.  Index THESE with ijb,
      ! never with a slot: they are allocated nx_block*ny_block in both evp and
      ! vp, so writing them at a compressed slot instead is in bounds and the
      ! compiler cannot see it -- it just puts every value at the wrong cell.
      ! That was the fff0dcd bug.
      !
      ! stepu does read them by slot, but through the separate umassdti_slot /
      ! waterxU_slot / wateryU_slot / forcexU_slot / forceyU_slot copies packed
      ! after this routine returns.  Those are additional arrays; they do not
      ! relayout these, and gather_dyn keeps seeing the block layout.
      real (kind=dbl_kind), dimension (:), intent(out) :: &
         Xmassdti, & ! mass of X-grid-cell/dt (kg/m^2 s)
         waterx  , & ! for ocean stress calculation, x (m/s)
         watery  , & ! for ocean stress calculation, y (m/s)
         forcex  , & ! work array: combined atm stress and ocn tilt, x
         forcey      ! work array: combined atm stress and ocn tilt, y

      ! Compressed copies of the five above, for the momentum kernel.  These
      ! are ADDITIONAL arrays: the block-layout originals are still written in
      ! full and still go whole to dyn_evp1d_run, so the rule stated above is
      ! not weakened.  Filling them here rather than packing afterwards saves a
      ! second traversal of the block -- the values are already to hand.
      !
      ! Each is this block's slice arr(off(iblk)+1:off(iblk)+cnt(iblk)+1), so
      ! element m is block-local slot m and the last element is the dead slot,
      ! exactly what pack_to_slots writes.  Optional because only the B-grid
      ! momentum kernel consumes them; the C/CD, eap and vp callers omit them.
      real (kind=dbl_kind), dimension (:), intent(out), optional :: &
         Xmassdti_slot, & ! Xmassdti at the cell's compressed slot
         waterx_slot  , & ! waterx   at the cell's compressed slot
         watery_slot  , & ! watery   at the cell's compressed slot
         forcex_slot  , & ! forcex   at the cell's compressed slot
         forcey_slot      ! forcey   at the cell's compressed slot

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij       , & ! ij is the compressed cell counter
         ijb                ! 1D index of (i,j) in the block

      real (kind=dbl_kind) :: &
         gravit

      logical (kind=log_kind), dimension(nx_block,ny_block) :: &
         iceXmask_old      ! old-time iceXmask

      ! the five compressed copies are supplied together or not at all, so one
      ! test taken once keeps the momentum loop free of present() checks
      logical (kind=log_kind) :: &
         do_slots          ! fill the compressed copies as well

      ! taubx and tauby are assumed shape, 1D over the block in Fortran
      ! storage order, i.e. element (i,j) is stored at (j-1)*nx_block+i
      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         taubx   , & ! seabed stress, x-direction (N/m^2)
         tauby       ! seabed stress, y-direction (N/m^2)

      character(len=*), parameter :: subname = '(dyn_prep2)'

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      ! compressed arrays are cleared over their own extent; every slot starts
      ! inactive and is switched on below only where there is ice
      waterx   (:) = c0
      watery   (:) = c0
      forcex   (:) = c0
      forcey   (:) = c0
      Xmassdti (:) = c0
      skipX    (:) = .true.

      ! The compressed copies are cleared over their own extent, dead slot
      ! included.  That reproduces pack_to_slots exactly: it sets the dead slot
      ! to zero and then copies every cell that owns a slot, and the block
      ! arrays it copies from have just been zeroed here, so the slots of cells
      ! without ice this step end up at zero either way.
      do_slots = present(Xmassdti_slot)
      if (do_slots) then
         Xmassdti_slot(:) = c0
         waterx_slot  (:) = c0
         watery_slot  (:) = c0
         forcex_slot  (:) = c0
         forcey_slot  (:) = c0
      endif

      do j = 1, ny_block
      do i = 1, nx_block
         ijb = (j-1)*nx_block + i
         TbU      (i,j) = c0
         taubx    (ijb) = c0
         tauby    (ijb) = c0

         if (.not.iceTmask(ijb)) then
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
      ! Identify cells where iceTmask = .true.
      ! Note: The icellT mask includes north and east ghost cells
      !       where stresses are needed.
      !-----------------------------------------------------------------

      icellT = 0
      do j = jlo, jhi+1
      do i = ilo, ihi+1
         ijb = (j-1)*nx_block + i
         if (iceTmask(ijb)) then
            icellT = icellT + 1
            indxTi(icellT) = i
            indxTj(icellT) = j
         endif
      enddo
      enddo

      !-----------------------------------------------------------------
      ! Define iceXmask
      ! Identify cells where iceXmask is true
      ! Initialize velocity where needed
      !-----------------------------------------------------------------

      icellX = 0

      do j = jlo, jhi
      do i = ilo, ihi
         ijb = (j-1)*nx_block + i
         iceXmask_old(i,j) = iceXmask(ijb) ! save
         ! ice extent mask (U-cells)
         iceXmask(ijb) = (Xmask(ijb)) .and. (aiX  (ijb) > dyn_area_min) &
                                      .and. (Xmass(ijb) > dyn_mass_min)

         if (iceXmask(ijb)) then
            icellX = icellX + 1
            indxXi(icellX) = i
            indxXj(icellX) = j

            ! initialize velocity for new ice points to ocean sfc current
            if (.not. iceXmask_old(i,j)) then
               uvel(i,j) = uocn(i,j)
               vvel(i,j) = vocn(i,j)
            endif

            ! compressed: this cell owns slot Xmap(ijb) for the whole run
            uvel_init(Xmap(ijb)) = uvel(i,j)
            vvel_init(Xmap(ijb)) = vvel(i,j)
            skipX    (Xmap(ijb)) = .false.
         else
            ! set velocity and stresses to zero for masked-out points
            uvel(i,j)    = c0
            vvel(i,j)    = c0
            strintx(i,j) = c0
            strinty(i,j) = c0
            strocnx(i,j) = c0
            strocny(i,j) = c0
         endif

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

      do ij = 1, icellX
         i = indxXi(ij)
         j = indxXj(ij)
         ijb = (j-1)*nx_block + i

         if (present(rheofactX)) then
            if ( aiX (ijb) > rheo_area_min ) then
               rheofactX(i,j) = c1
            else
               rheofactX(i,j) = c0
            endif
         endif

         Xmassdti(ijb) = Xmass(ijb)/dt ! kg/m^2 s

         fm(i,j) = fcor(ijb)*Xmass(ijb)   ! Coriolis * mass

         ! for ocean stress
         waterx(ijb) = uocn(i,j)*cosw - vocn(i,j)*sinw*sign(c1,fm(i,j))
         watery(ijb) = vocn(i,j)*cosw + uocn(i,j)*sinw*sign(c1,fm(i,j))

         ! combine tilt with wind stress
         if (trim(ssh_stress) == 'geostrophic') then
            ! calculate tilt from geostrophic currents if needed
            strtltx(i,j) = -fm(i,j)*vocn(i,j)
            strtlty(i,j) =  fm(i,j)*uocn(i,j)
         elseif (trim(ssh_stress) == 'coupled') then
            strtltx(i,j) = -gravit*Xmass(ijb)*ss_tltx(i,j)
            strtlty(i,j) = -gravit*Xmass(ijb)*ss_tlty(i,j)
         else
            call abort_ice(subname//' ERROR: unknown ssh_stress='//trim(ssh_stress), &
               file=__FILE__, line=__LINE__)
         endif

         forcex(ijb) = strairx(i,j) + strtltx(i,j)
         forcey(ijb) = strairy(i,j) + strtlty(i,j)

         ! copy, never relayout: the block arrays above keep every value, and
         ! these are read back from them so the two cannot drift apart
         if (do_slots) then
            Xmassdti_slot(Xmap(ijb)) = Xmassdti(ijb)
            waterx_slot  (Xmap(ijb)) = waterx  (ijb)
            watery_slot  (Xmap(ijb)) = watery  (ijb)
            forcex_slot  (Xmap(ijb)) = forcex  (ijb)
            forcey_slot  (Xmap(ijb)) = forcey  (ijb)
         endif
      enddo

      end subroutine dyn_prep2

!=======================================================================
! Calculation of the surface stresses
! Integration of the momentum equation to find velocity (u,v)
!
! author: Elizabeth C. Hunke, LANL

      subroutine stepu (&
                        ncells,               &
                        Cw,                   &
                        skipme,               &
                        aiX,        str,      &
                        uocn,       vocn,     &
                        waterx,     watery,   &
                        forcex,     forcey,   &
                        Umassdti,   fm,       &
                        uarear,               &
                        strintx_s,  strinty_s,&
                        taubx_s,    tauby_s,  &
                        uvel_init,  vvel_init,&
                        uvel,       vvel,     &
                        nbr,                  &
                        TbU)

      integer (kind=int_kind), intent(in) :: &
         ncells                ! slots held by this block

      ! skipme(is) is .true. where the slot carries no ice this step.
      ! Together they replace the icellU / indxUi / indxUj / Xmap chain: this
      ! kernel iterates over slots, the way ice_dyn_core1d's does.
      logical (kind=log_kind), dimension (:), intent(in), contiguous :: &
         skipme      ! inactive-slot flag

      ! Every coefficient this kernel reads is COMPRESSED: this block's slice,
      ! 1-based, cells 1..cnt with the dead slot last, matching the velocities.  All of them are constant
      ! through the subcycling -- the grid metric uarear is constant outright
      ! -- so the caller packs them once before the loop.  With these there is
      ! no rank-2 read left in the loop body, which is what allows the cell
      ! loop itself to run on slots.
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         aiX     , & ! ice fraction on X-grid
         fm      , & ! Coriolis param. * mass in U-cell (kg/s)
         TbU     , & ! seabed stress factor (N/m^2)
         Umassdti, & ! mass of U-cell/dt (kg/m^2 s)
         waterx  , & ! for ocean stress calculation, x (m/s)
         watery  , & ! for ocean stress calculation, y (m/s)
         forcex  , & ! work array: combined atm stress and ocn tilt, x
         forcey  , & ! work array: combined atm stress and ocn tilt, y
         uocn    , & ! ocean current, x-direction (m/s)
         vocn    , & ! ocean current, y-direction (m/s)
         Cw      , & ! ocean-ice neutral drag coefficient
         uarear      ! 1/uarea

      ! uvel_init/vvel_init are compressed but indexed from 1
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         uvel_init,& ! x-component of velocity (m/s), beginning of timestep
         vvel_init   ! y-component of velocity (m/s), beginning of timestep

      ! str is COMPRESSED and indexed by BLOCK-LOCAL slot.  Element 0 is the
      ! "no such neighbour" slot and is held at zero, which is what the block
      ! layout carried at land points, so the stencil below needs no branch.
      real (kind=dbl_kind), dimension(:,:), intent(in), contiguous :: &
         str         ! temporary

      ! nbr(is,d) is the block-local slot one step in direction d from slot is
      integer (kind=int_kind), dimension(:,:), intent(in) :: &
         nbr         ! neighbour slots

      ! uvel and vvel are COMPRESSED: this block's slice, 1-based, so cells
      ! are 1..cnt and the last element is the block's dead slot
      real (kind=dbl_kind), dimension (:), intent(inout), contiguous :: &
         uvel    , & ! x-component of velocity (m/s)
         vvel        ! y-component of velocity (m/s)

      ! compressed, output only -- see the store at the end of the loop
      real (kind=dbl_kind), dimension (:), intent(inout), contiguous :: &
         strintx_s , & ! divergence of internal ice stress, x (N/m^2)
         strinty_s     ! divergence of internal ice stress, y (N/m^2)

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij          , & ! ij is the compressed cell counter
         is                  , & ! compressed slot
         isE, isN, isNE          ! its E, N and NE neighbour slots

      real (kind=dbl_kind) :: &
         uold, vold        , & ! old-time uvel, vvel
         strintx_tmp       , & ! divergence of internal stress, x, this cell
         strinty_tmp       , & ! divergence of internal stress, y, this cell
         tmp_str2_E        , & ! the six neighbour stresses, loaded up front
         tmp_str3_N        , & !
         tmp_str4_NE       , & !
         tmp_str6_N        , & !
         tmp_str7_E        , & !
         tmp_str8_NE       , & !
         vrel              , & ! relative ice-ocean velocity
         cca,ccb,ab2,cc1,cc2,& ! intermediate variables
         taux, tauy        , & ! part of ocean stress term
         Cb                , & ! complete seabed (basal) stress coeff
         rhow                  !

      ! compressed, output only -- see the store at the end of the loop
      real (kind=dbl_kind), dimension (:), intent(inout), contiguous :: &
         taubx_s , & ! seabed stress, x-direction (N/m^2)
         tauby_s     ! seabed stress, y-direction (N/m^2)

      character(len=*), parameter :: subname = '(stepu)'

      !-----------------------------------------------------------------
      ! integrate the momentum equation
      !-----------------------------------------------------------------

      call icepack_query_parameters(rhow_out=rhow)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      do is = 1, ncells
         if (skipme(is)) cycle

         uold = uvel(is)
         vold = vvel(is)

         ! (magnitude of relative ocean current)*rhow*drag*aice
         vrel = aiX(is)*rhow*Cw(is)*sqrt((uocn(is) - uold)**2 + &
                                           (vocn(is) - vold)**2)  ! m/s
         ! ice/ocean stress
         taux = vrel*waterx(is) ! NOTE this is not the entire
         tauy = vrel*watery(is) ! ocn stress term

         Cb  = TbU(is) / (sqrt(uold**2 + vold**2) + u0) ! for seabed stress
         ! revp = 0 for classic evp, 1 for revised evp
         cca = (brlx + revp)*Umassdti(is) + vrel * cosw + Cb ! kg/m^2 s

         ccb = fm(is) + sign(c1,fm(is)) * vrel * sinw ! kg/m^2 s

         ab2 = cca**2 + ccb**2

         ! divergence of the internal stress tensor
         ! sum over the four T cells around this U point: (0,0), (+1,0),
         ! (0,+1), (+1,+1).  A missing neighbour resolves to slot 0, which
         ! holds zero, matching what the block layout held at land.
         ! The three neighbour slots are read once each and the six stresses
         ! loaded up front, as ice_dyn_core1d's stepu_1d does.  Each slot was
         ! fetched twice before -- str(nbr(is,nbr_E),2) and
         ! str(nbr(is,nbr_E),7) -- so this halves the nbr loads, and it lets
         ! the six dependent str loads issue before any is needed.
         isE  = nbr(is,nbr_E )
         isN  = nbr(is,nbr_N )
         isNE = nbr(is,nbr_NE)

         tmp_str2_E  = str(isE ,2)
         tmp_str3_N  = str(isN ,3)
         tmp_str4_NE = str(isNE,4)
         tmp_str6_N  = str(isN ,6)
         tmp_str7_E  = str(isE ,7)
         tmp_str8_NE = str(isNE,8)

         strintx_tmp = uarear(is)* &
             (str(is,1) + tmp_str2_E + tmp_str3_N + tmp_str4_NE)
         strinty_tmp = uarear(is)* &
             (str(is,5) + tmp_str6_N + tmp_str7_E + tmp_str8_NE)

         ! finally, the velocity components
         cc1 = strintx_tmp + forcex(is) + taux &
             + Umassdti(is)*(brlx*uold + revp*uvel_init(is))
         cc2 = strinty_tmp + forcey(is) + tauy &
             + Umassdti(is)*(brlx*vold + revp*vvel_init(is))

         uvel(is) = (cca*cc1 + ccb*cc2) / ab2 ! m/s
         vvel(is) = (cca*cc2 - ccb*cc1) / ab2

         ! These four are OUTPUT ONLY -- history reads them, nothing in the
         ! subcycle does.  They used to be written to the block arrays under
         ! "if (ksub == ndte)", which needed ksub and ndte passed in and cost
         ! a divide and a mod per cell to rebuild (i,j) from the slot.  Writing
         ! slots instead makes the guard pointless: the last pass's value is
         ! simply what is left there, and the driver unpacks once after the
         ! loop.  No branch, no index arithmetic.
         strintx_s(is) = strintx_tmp
         strinty_s(is) = strinty_tmp
         taubx_s  (is) = -uvel(is)*Cb
         tauby_s  (is) = -vvel(is)*Cb
      enddo                     ! ij

      end subroutine stepu

!=======================================================================
! Integration of the momentum equation to find velocity (u,v) at E and N locations

      subroutine stepuv_CD (ncells,     skipme,   &
                            Cw,                   &
                                        aiX,      &
                            uocn,       vocn,     &
                            waterx,     watery,   &
                            forcex,     forcey,   &
                            massdti,    fm,       &
                            strintx_s,  strinty_s,&
                            taubx_s,    tauby_s,  &
                            uvel_init,  vvel_init,&
                            uvel_s,     vvel_s,   &
                            Tb)

      integer (kind=int_kind), intent(in) :: &
         ncells                ! number of slots in this block

      ! skipme(is) is .true. where the slot carries no ice this step
      logical (kind=log_kind), dimension (:), intent(in), contiguous :: &
         skipme      ! inactive-slot flag

      ! Every coefficient below is COMPRESSED, this block's slice.  All are
      ! constant through the subcycling, so the caller packs them once; with
      ! them there is no rank-2 read left in the loop, which is what lets the
      ! cell loop run on slots exactly as the B grid's stepu does.
      ! aiX and massdti are compressed like the rest
      real (kind=dbl_kind), dimension (:), intent(in) :: &
         aiX     , & ! ice fraction on X-grid
         massdti , & ! mass of [EN]-cell/dt (kg/m^2 s)
         uvel_init,& ! x-component of velocity (m/s), beginning of timestep
         vvel_init,& ! y-component of velocity (m/s), beginning of timestep
         waterx  , & ! for ocean stress calculation, x (m/s)
         watery  , & ! for ocean stress calculation, y (m/s)
         forcex  , & ! work array: combined atm stress and ocn tilt, x
         forcey      ! work array: combined atm stress and ocn tilt, y

      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         Tb,       & ! seabed stress factor (N/m^2)
         uocn    , & ! ocean current, x-direction (m/s)
         vocn    , & ! ocean current, y-direction (m/s)
         fm          ! Coriolis param. * mass in [EN]-cell (kg/s)

      ! compressed, written by div_stress_* in the same pass
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         strintx_s , & ! divergence of internal ice stress, x (N/m^2)
         strinty_s     ! divergence of internal ice stress, y (N/m^2)

      ! Compressed: the velocities are what the subcycle carries between the
      ! step kernels, strain_rates_U and the stress kernels, so they stay in
      ! slots for the whole loop rather than round-tripping through blocks.
      real (kind=dbl_kind), dimension (:), intent(inout), contiguous :: &
         uvel_s  , & ! x-component of velocity (m/s)
         vvel_s      ! y-component of velocity (m/s)

      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         Cw          ! ocean-ice neutral drag coefficient

      ! local variables

      integer (kind=int_kind) :: &
         is                   ! compressed slot

      real (kind=dbl_kind) :: &
         uold, vold         , & ! old-time uvel, vvel
         vrel               , & ! relative ice-ocean velocity
         cca,ccb,ccc,ab2    , & ! intermediate variables
         cc1,cc2            , & ! "
         taux, tauy         , & ! part of ocean stress term
         Cb                 , & ! complete seabed (basal) stress coeff
         rhow                   !

      ! compressed, output only: written every pass, unpacked once after
      real (kind=dbl_kind), dimension (:), intent(inout), contiguous :: &
         taubx_s , & ! seabed stress, x-direction (N/m^2)
         tauby_s     ! seabed stress, y-direction (N/m^2)

      character(len=*), parameter :: subname = '(stepuv_CD)'

      !-----------------------------------------------------------------
      ! integrate the momentum equation
      !-----------------------------------------------------------------

      call icepack_query_parameters(rhow_out=rhow)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      do is = 1, ncells
         if (skipme(is)) cycle

         uold = uvel_s(is)
         vold = vvel_s(is)

         ! (magnitude of relative ocean current)*rhow*drag*aice
         vrel = aiX(is)*rhow*Cw(is)*sqrt((uocn(is) - uold)**2 + &
                                           (vocn(is) - vold)**2)  ! m/s
         ! ice/ocean stress
         taux = vrel*waterx(is) ! NOTE this is not the entire
         tauy = vrel*watery(is) ! ocn stress term

         ccc = sqrt(uold**2 + vold**2) + u0
         Cb  = Tb(is) / ccc ! for seabed stress
         ! revp = 0 for classic evp, 1 for revised evp
         cca = (brlx + revp)*massdti(is) + vrel * cosw + Cb ! kg/m^2 s

         ccb = fm(is) + sign(c1,fm(is)) * vrel * sinw ! kg/m^2 s

         ab2 = cca**2 + ccb**2

         ! compute the velocity components
         cc1 = strintx_s(is) + forcex(is) + taux &
             + massdti(is)*(brlx*uold + revp*uvel_init(is))
         cc2 = strinty_s(is) + forcey(is) + tauy &
             + massdti(is)*(brlx*vold + revp*vvel_init(is))
         uvel_s(is) = (cca*cc1 + ccb*cc2) / ab2 ! m/s
         vvel_s(is) = (cca*cc2 - ccb*cc1) / ab2

         ! calculate seabed stress component for outputs
         ! only needed on last iteration.
         taubx_s(is) = -uvel_s(is)*Cb
         tauby_s(is) = -vvel_s(is)*Cb

      enddo                     ! ij

      end subroutine stepuv_CD

!=======================================================================
! Integration of the momentum equation to find velocity u at E location on C grid

      subroutine stepu_C (ncells,     skipme,   &
                          Cw,                   &
                                      aiE,      &
                          uocn,       vocn,     &
                          waterx,     forcex,   &
                          massdti,    fm,       &
                          strintx_s,  taubx_s,  &
                          uvel_init,            &
                          uvel_s,     vvel_s,   &
                          Tb)

      integer (kind=int_kind), intent(in) :: &
         ncells                ! number of slots in this block

      ! skipme(is) is .true. where the slot carries no ice this step
      logical (kind=log_kind), dimension (:), intent(in), contiguous :: &
         skipme      ! inactive-slot flag

      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         aiE     , & ! ice fraction on E-grid
         massdti , & ! mass of e-cell/dt (kg/m^2 s)
         uvel_init,& ! x-component of velocity (m/s), beginning of timestep
         waterx  , & ! for ocean stress calculation, x (m/s)
         forcex      ! work array: combined atm stress and ocn tilt, x

      ! Every coefficient below is COMPRESSED, this block's slice.  All are
      ! constant through the subcycling, so the caller packs them once; with
      ! them there is no rank-2 read left in the loop, which is what lets the
      ! cell loop run on slots exactly as the B grid's stepu does.
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         Tb,       & ! seabed stress factor (N/m^2)
         uocn    , & ! ocean current, x-direction (m/s)
         vocn    , & ! ocean current, y-direction (m/s)
         fm      , & ! Coriolis param. * mass in e-cell (kg/s)
         Cw          ! ocean-ice neutral drag coefficient

      ! compressed, written by div_stress_Ex in the same pass
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         strintx_s   ! divergence of internal ice stress, x (N/m^2)

      ! Compressed: the velocities are what the subcycle carries between
      ! stepu_C, strain_rates_U and the stress kernels, so they stay in slots
      ! for the whole loop rather than round-tripping through the block form.
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         vvel_s      ! y-component of velocity (m/s) interpolated to E location

      ! compressed, output only: written every pass, unpacked once after the
      ! subcycle.  Nothing reads it in between, here or anywhere else.
      real (kind=dbl_kind), dimension (:), intent(inout), contiguous :: &
         taubx_s     ! seabed stress, x-direction (N/m^2)

      real (kind=dbl_kind), dimension (:), intent(inout), contiguous :: &
         uvel_s      ! x-component of velocity (m/s)

      ! local variables

      integer (kind=int_kind) :: &
         is                   ! compressed slot

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

      do is = 1, ncells
         if (skipme(is)) cycle

         uold = uvel_s(is)
         vold = vvel_s(is)

         ! (magnitude of relative ocean current)*rhow*drag*aice
         vrel = aiE(is)*rhow*Cw(is)*sqrt((uocn(is) - uold)**2 + &
                                           (vocn(is) - vold)**2)  ! m/s
         ! ice/ocean stress
         taux = vrel*waterx(is) ! NOTE this is not the entire

         ccc = sqrt(uold**2 + vold**2) + u0
         Cb  = Tb(is) / ccc ! for seabed stress
         ! revp = 0 for classic evp, 1 for revised evp
         cca = (brlx + revp)*massdti(is) + vrel * cosw + Cb ! kg/m^2 s

         ccb = fm(is) + sign(c1,fm(is)) * vrel * sinw ! kg/m^2 s

         ! compute the velocity components
         cc1 = strintx_s(is) + forcex(is) + taux &
             + massdti(is)*(brlx*uold + revp*uvel_init(is))

         uvel_s(is) = (ccb*vold + cc1) / cca ! m/s

         ! calculate seabed stress component for outputs
         ! only needed on last iteration.
         taubx_s(is) = -uvel_s(is)*Cb

      enddo                     ! ij

      end subroutine stepu_C

!=======================================================================
! Integration of the momentum equation to find velocity v at N location on C grid

      subroutine stepv_C (ncells,     skipme,   &
                          Cw,                   &
                                      aiN,      &
                          uocn,       vocn,     &
                          watery,     forcey,   &
                          massdti,    fm,       &
                          strinty_s,  tauby_s,  &
                          vvel_init,            &
                          uvel_s,     vvel_s,   &
                          Tb)

      integer (kind=int_kind), intent(in) :: &
         ncells                ! number of slots in this block

      ! skipme(is) is .true. where the slot carries no ice this step
      logical (kind=log_kind), dimension (:), intent(in), contiguous :: &
         skipme      ! inactive-slot flag

      ! Every coefficient below is COMPRESSED, this block's slice.  All are
      ! constant through the subcycling, so the caller packs them once; with
      ! them there is no rank-2 read left in the loop, which is what lets the
      ! cell loop run on slots exactly as the B grid's stepu does.
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         aiN     , & ! ice fraction on N-grid
         massdti , & ! mass of n-cell/dt (kg/m^2 s)
         vvel_init,& ! y-component of velocity (m/s), beginning of timestep
         watery  , & ! for ocean stress calculation, y (m/s)
         forcey      ! work array: combined atm stress and ocn tilt, y

      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         Tb,       & ! seabed stress factor (N/m^2)
         uocn    , & ! ocean current, x-direction (m/s)
         vocn    , & ! ocean current, y-direction (m/s)
         fm      , & ! Coriolis param. * mass in n-cell (kg/s)
         Cw          ! ocean-ice neutral drag coefficient

      ! compressed, written by div_stress_Ny in the same pass
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         strinty_s   ! divergence of internal ice stress, y (N/m^2)

      ! Compressed: the velocities are what the subcycle carries between the
      ! step kernels, strain_rates_U and the stress kernels, so they stay in
      ! slots for the whole loop rather than round-tripping through blocks.
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         uvel_s      ! x-component of velocity (m/s) interpolated to N location

      ! compressed, output only: written every pass, unpacked once after the
      ! subcycle.  Nothing reads it in between, here or anywhere else.
      real (kind=dbl_kind), dimension (:), intent(inout), contiguous :: &
         tauby_s     ! seabed stress, y-direction (N/m^2)

      real (kind=dbl_kind), dimension (:), intent(inout), contiguous :: &
         vvel_s      ! y-component of velocity (m/s)

      ! local variables

      integer (kind=int_kind) :: &
         is                   ! compressed slot

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

      do is = 1, ncells
         if (skipme(is)) cycle

         uold = uvel_s(is)
         vold = vvel_s(is)

         ! (magnitude of relative ocean current)*rhow*drag*aice
         vrel = aiN(is)*rhow*Cw(is)*sqrt((uocn(is) - uold)**2 + &
                                           (vocn(is) - vold)**2)  ! m/s
         ! ice/ocean stress
         tauy = vrel*watery(is) ! NOTE this is not the entire ocn stress

         ccc = sqrt(uold**2 + vold**2) + u0
         Cb  = Tb(is) / ccc ! for seabed stress
         ! revp = 0 for classic evp, 1 for revised evp
         cca = (brlx + revp)*massdti(is) + vrel * cosw + Cb ! kg/m^2 s

         ccb = fm(is) + sign(c1,fm(is)) * vrel * sinw ! kg/m^2 s

         ! compute the velocity components
         cc2 = strinty_s(is) + forcey(is) + tauy &
             + massdti(is)*(brlx*vold + revp*vvel_init(is))

         vvel_s(is) = (-ccb*uold + cc2) / cca

         ! calculate seabed stress component for outputs
         ! only needed on last iteration.
         tauby_s(is) = -vvel_s(is)*Cb

      enddo                     ! ij

      end subroutine stepv_C

!=======================================================================
! Calculation of the ice-ocean stress.
! ...the sign will be reversed later...
!
! author: Elizabeth C. Hunke, LANL

      subroutine dyn_finish (nx_block, ny_block, &
                             icellU,   Cw,       &
                             indxUi,   indxUj,   &
                             uvel,     vvel,     &
                             uocn,     vocn,     &
                             aiX,      fm,       &
                             strocnx,  strocny)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellU                ! total count when iceUmask is true

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxUi  , & ! compressed index in i-direction
         indxUj      ! compressed index in j-direction

      ! aiX is assumed shape, 1D over the block in Fortran storage order,
      ! i.e. element (i,j) is stored at (j-1)*nx_block+i
      real (kind=dbl_kind), dimension (:), intent(in) :: &
         aiX         ! ice fraction on X-grid

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         uvel    , & ! x-component of velocity (m/s)
         vvel    , & ! y-component of velocity (m/s)
         uocn    , & ! ocean current, x-direction (m/s)
         vocn    , & ! ocean current, y-direction (m/s)
         fm          ! Coriolis param. * mass in U-cell (kg/s)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         strocnx , & ! ice-ocean stress, x-direction
         strocny     ! ice-ocean stress, y-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         Cw          ! ocean-ice neutral drag coefficient

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij       , & ! ij is the compressed cell counter
         ijb                ! 1D index of (i,j) in the block

      real (kind=dbl_kind) :: &
         vrel    , & !
         rhow        !

      character(len=*), parameter :: subname = '(dyn_finish)'

      call icepack_query_parameters(rhow_out=rhow)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      ! ocean-ice stress for coupling
      do ij =1, icellU
         i = indxUi(ij)
         j = indxUj(ij)
         ijb = (j-1)*nx_block + i

          vrel = rhow*Cw(i,j)*sqrt((uocn(i,j) - uvel(i,j))**2 + &
                 (vocn(i,j) - vvel(i,j))**2)  ! m/s

!        strocnx(i,j) = strocnx(i,j) &
!                     - vrel*(uvel(i,j)*cosw - vvel(i,j)*sinw) * aiX(ijb)
!        strocny(i,j) = strocny(i,j) &
!                     - vrel*(vvel(i,j)*cosw + uvel(i,j)*sinw) * aiX(ijb)

         ! update strocnx to most recent iterate and complete the term
         vrel = vrel * aiX(ijb)
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
! Computes seabed (basal) stress factor TbU (landfast ice) based on mean
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
! note1: TbU is a part of the Cb as defined in Lemieux et al. 2015 and 2016.
! note2: Seabed stress (better name) was called basal stress in Lemieux et al. 2015

      subroutine seabed_stress_factor_LKD (nx_block, ny_block,         &
                                           icellU,                     &
                                           indxUi,   indxUj,           &
                                           vice,     aice,             &
                                           hwater,   TbU,              &
                                           grid_location)

      use ice_grid, only: grid_neighbor_min, grid_neighbor_max

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, &  ! block dimensions
         icellU                 ! no. of cells where ice[uen]mask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxUi    , & ! compressed index in i-direction
         indxUj        ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         aice      , & ! concentration of ice at tracer location
         vice      , & ! volume per unit area of ice at tracer location (m)
         hwater        ! water depth at tracer location (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         TbU           ! seabed stress factor at 'grid_location' (N/m^2)

      character(len=*), optional, intent(in) :: &
         grid_location ! grid location (U, E, N), U assumed if not present

      real (kind=dbl_kind) :: &
         au        , & ! concentration of ice at u location
         hu        , & ! volume per unit area of ice at u location (mean thickness, m)
         hwu       , & ! water depth at u location (m)
         docalc_tbu, & ! logical as real (C0,C1) decides whether c0 is 0 or
         hcu           ! critical thickness at u location (m)

      integer (kind=int_kind) :: &
         i, j, ij       , & ! ij is the compressed cell counter
         ijb                ! 1D index of (i,j) in the block

      character(len=char_len) :: &
         l_grid_location ! local version of 'grid_location'

      character(len=*), parameter :: subname = '(seabed_stress_factor_LKD)'

      ! Assume U location (NE corner) if grid_location not present
      if (.not. (present(grid_location))) then
         l_grid_location = 'U'
      else
         l_grid_location = grid_location
      endif

      do ij = 1, icellU
         i = indxUi(ij)
         j = indxUj(ij)

         ! convert quantities to grid_location

         hwu = grid_neighbor_min(hwater, i, j, l_grid_location)

         docalc_tbu = merge(c1,c0,hwu < threshold_hw)


         au  = grid_neighbor_max(aice, i, j, l_grid_location)
         hu  = grid_neighbor_max(vice, i, j, l_grid_location)

         ! 1- calculate critical thickness
         hcu = au * hwu / k1

         ! 2- calculate seabed stress factor
         TbU(i,j) = docalc_tbu*k2 * max(c0,(hu - hcu)) * exp(-alphab * (c1 - au))

      enddo                     ! ij

      end subroutine seabed_stress_factor_LKD

!=======================================================================
! Computes seabed (basal) stress factor TbU (landfast ice) based on
! probability of contact between the ITD and the seabed. The water depth
! could take into account variations of the SSH. In the simplest
! formulation, hwater is simply the value of the bathymetry. To calculate
! the probability of contact, it is assumed that the bathymetry follows
! a normal distribution with sigma_b = 2.5d0. An improvement would
! be to provide the distribution based on high resolution data.
!
! Dupont, F., D. Dumont, J.F. Lemieux, E. Dumas-Lefebvre, A. Caya (2022).
! A probabilistic seabed-ice keel interaction model, The Cryosphere, 16,
! 1963-1977.
!
! authors: D. Dumont, J.F. Lemieux, E. Dumas-Lefebvre, F. Dupont
!
      subroutine seabed_stress_factor_prob (nx_block, ny_block,          &
                                            icellT, indxTi,   indxTj,    &
                                            icellU, indxUi,   indxUj,    &
                                            aicen,  vicen,               &
                                            hwater, TbU,                 &
                                            TbE,    TbN,                 &
                                            icellE, indxEi,   indxEj,    &
                                            icellN, indxNi,   indxNj)
! use modules

      use ice_arrays_column, only: hin_max
      use ice_domain_size, only: ncat
      use ice_grid, only: grid_neighbor_min, grid_neighbor_max

      integer (kind=int_kind), intent(in) :: &
           nx_block, ny_block, &  ! block dimensions
           icellT, icellU         ! no. of cells where ice[tu]mask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
           indxTi  , & ! compressed index in i-direction
           indxTj  , & ! compressed index in j-direction
           indxUi  , & ! compressed index in i-direction
           indxUj      ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
           hwater      ! water depth at tracer location (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), intent(in) :: &
           aicen,    & ! partial concentration for last thickness category in ITD
           vicen       ! partial volume for last thickness category in ITD (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
           TbU         ! seabed stress factor at U location (N/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout), optional :: &
           TbE,      & ! seabed stress factor at E location (N/m^2)
           TbN         ! seabed stress factor at N location (N/m^2)

      integer (kind=int_kind), intent(in), optional :: &
           icellE, icellN ! no. of cells where ice[en]mask = 1

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in), optional :: &
           indxEi  , & ! compressed index in i-direction
           indxEj  , & ! compressed index in j-direction
           indxNi  , & ! compressed index in i-direction
           indxNj      ! compressed index in j-direction

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

      call icepack_query_parameters(rhow_out=rhow, rhoi_out=rhoi,gravit_out=gravit,pi_out=pi,puny_out=puny)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      Tbt=c0

      do ij = 1, icellT
         i = indxTi(ij)
         j = indxTj(ij)

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
            v_i = max((v_i - m_i**2), puny)

            mu_i    = log(m_i/sqrt(c1 + v_i/m_i**2)) ! parameters for the log-normal
            sigma_i = sqrt(log(c1 + v_i/m_i**2))

            ! max thickness associated with percentile of log-normal PDF
            ! x_kmax=x997 was obtained from an optimization procedure (Dupont et al. 2022)

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
         do ij = 1, icellU
            i = indxUi(ij)
            j = indxUj(ij)
            ! convert quantities to U-location
            TbU(i,j)  = grid_neighbor_max(Tbt, i, j, 'U')
         enddo                     ! ij
      elseif (grid_ice == "C" .or. grid_ice == "CD") then
         if (present(Tbe)    .and. present(TbN)    .and. &
             present(icellE) .and. present(icellN) .and. &
             present(indxEi) .and. present(indxEj) .and. &
             present(indxNi) .and. present(indxNj)) then

            do ij = 1, icellE
               i = indxEi(ij)
               j = indxEj(ij)
               ! convert quantities to E-location
                  TbE(i,j)  = grid_neighbor_max(Tbt, i, j, 'E')
            enddo
            do ij = 1, icellN
               i = indxNi(ij)
               j = indxNj(ij)
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
                               icellT,                 &
                               indxTi,     indxTj,     &
                               uvel,       vvel,       &
                               dxT,        dyT,        &
                               dxU,        dyU,        &
                               cxp,        cyp,        &
                               cxm,        cym,        &
                               tarear,     vort,       &
                               shear,      divu,       &
                               rdg_conv,   rdg_shear )

      use ice_constants, only: p25, p5

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellT                ! no. of cells where iceTmask = .true.

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxTi   , & ! compressed index in i-direction
         indxTj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(in) :: &
         uvel     , & ! x-component of velocity (m/s)
         vvel     , & ! y-component of velocity (m/s)
         dxT      , & ! width of T-cell through the middle (m)
         dyT      , & ! height of T-cell through the middle (m)
         dxU      , & ! width of U-cell through the middle (m)
         dyU      , & ! height of U-cell through the middle (m)
         cyp      , & ! 1.5*HTE - 0.5*HTW
         cxp      , & ! 1.5*HTN - 0.5*HTS
         cym      , & ! 0.5*HTE - 1.5*HTW
         cxm          ! 0.5*HTN - 1.5*HTS

      ! tarear is assumed shape, 1D over the block in Fortran storage order,
      ! i.e. element (i,j) is stored at (j-1)*nx_block+i
      real (kind=dbl_kind), dimension (:), intent(in) :: &
         tarear       ! 1/tarea

      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(inout) :: &
         vort     , & ! vorticity (1/s)
         shear    , & ! strain rate II component (1/s)
         divu     , & ! strain rate I component, velocity divergence (1/s)
         rdg_conv , & ! convergence term for ridging (1/s)
         rdg_shear    ! shear term for ridging (1/s)

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij, ijb

      real (kind=dbl_kind) :: &                       ! at each corner :
        divune, divunw, divuse, divusw            , & ! divergence
        tensionne, tensionnw, tensionse, tensionsw, & ! tension
        shearne, shearnw, shearse, shearsw        , & ! shearing
        Deltane, Deltanw, Deltase, Deltasw        , & ! Delta
        tmp                                           ! useful combination

      real (kind=dbl_kind) :: &                       ! at edges for vorticity calc :
         dvdxn, dvdxs, dudye, dudyw                   ! dvdx and dudy terms on edges

      character(len=*), parameter :: subname = '(deformations)'

      do ij = 1, icellT
         i = indxTi(ij)
         j = indxTj(ij)
         ijb = (j-1)*nx_block + i

         !-----------------------------------------------------------------
         ! strain rates
         ! NOTE these are actually strain rates * area  (m^2/s)
         !-----------------------------------------------------------------
         call strain_rates (nx_block,   ny_block,   &
                            ijb,                    &
                               uvel(ijb),  uvel(ijb-1),          &
                               uvel(ijb-nx_block), uvel(ijb-nx_block-1), &
                               vvel(ijb),  vvel(ijb-1),          &
                               vvel(ijb-nx_block), vvel(ijb-nx_block-1), &
                            dxT,        dyT,        &
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
         divu(ijb) = p25*(divune + divunw + divuse + divusw) * tarear(ijb)
         tmp = p25*(Deltane + Deltanw + Deltase + Deltasw)   * tarear(ijb)
         rdg_conv(ijb)  = -min(divu(ijb),c0)
         rdg_shear(ijb) = p5*(tmp-abs(divu(ijb)))

         ! diagnostic only
         ! shear = sqrt(tension**2 + shearing**2)
         shear(ijb) = p25*tarear(ijb)*sqrt( &
                      (tensionne + tensionnw + tensionse + tensionsw)**2 + &
                      (shearne   + shearnw   + shearse   + shearsw  )**2)

         ! vorticity
         dvdxn = dyU(ijb)*vvel(ijb) - dyU(ijb-1)*vvel(ijb-1)
         dvdxs = dyU(ijb-nx_block)*vvel(ijb-nx_block) - dyU(ijb-nx_block-1)*vvel(ijb-nx_block-1)
         dudye = dxU(ijb)*uvel(ijb) - dxU(ijb-nx_block)*uvel(ijb-nx_block)
         dudyw = dxU(ijb-1)*uvel(ijb-1) - dxU(ijb-nx_block-1)*uvel(ijb-nx_block-1)
         vort(ijb) = p5*tarear(ijb)*(dvdxn + dvdxs - dudye - dudyw)

      enddo                     ! ij

      end subroutine deformations

!=======================================================================
! Compute deformations for mechanical redistribution at T point
!
! author: JF Lemieux, ECCC
! Nov 2021

      subroutine deformationsCD_T (nx_block,   ny_block,   &
                                   icellT,                 &
                                   indxTi,     indxTj,     &
                                   uvelE,      vvelE,      &
                                   uvelN,      vvelN,      &
                                   dxN,        dyE,        &
                                   dxT,        dyT,        &
                                   tarear,     vort,       &
                                   shear,      divu,       &
                                   rdg_conv,   rdg_shear )

      use ice_constants, only: p5

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellT                ! no. of cells where iceTmask = .true.

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxTi   , & ! compressed index in i-direction
         indxTj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(in) :: &
         uvelE    , & ! x-component of velocity (m/s) at the E point
         vvelE    , & ! y-component of velocity (m/s) at the N point
         uvelN    , & ! x-component of velocity (m/s) at the E point
         vvelN    , & ! y-component of velocity (m/s) at the N point
         dxN      , & ! width of N-cell through the middle (m)
         dyE      , & ! height of E-cell through the middle (m)
         dxT      , & ! width of T-cell through the middle (m)
         dyT      , & ! height of T-cell through the middle (m)
         tarear       ! 1/tarea

      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(inout) :: &
         vort     , & ! vorticity (1/s)
         shear    , & ! strain rate II component (1/s)
         divu     , & ! strain rate I component, velocity divergence (1/s)
         rdg_conv , & ! convergence term for ridging (1/s)
         rdg_shear    ! shear term for ridging (1/s)

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij       , & ! ij is the compressed cell counter
         ijb                ! 1D index of (i,j) in the block

      real (kind=dbl_kind), dimension (nx_block*ny_block) :: &
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
                           icellT     ,                &
                           indxTi(:)  , indxTj  (:)  , &
                           uvelE (:)  , vvelE   (:)  , &
                           uvelN (:)  , vvelN   (:)  , &
                           dxN   (:)  , dyE     (:)  , &
                           dxT   (:)  , dyT     (:)  , &
                           divT  (:)  , tensionT(:)  , &
                           shearT(:)  , DeltaT  (:)    )

      do ij = 1, icellT
         i = indxTi(ij)
         j = indxTj(ij)
         ijb = (j-1)*nx_block + i

         !-----------------------------------------------------------------
         ! deformations for mechanical redistribution
         !-----------------------------------------------------------------
         divu(ijb) = divT(ijb) * tarear(ijb)
         tmp = Deltat(ijb) * tarear(ijb)
         rdg_conv(ijb)  = -min(divu(ijb),c0)
         rdg_shear(ijb) = p5*(tmp-abs(divu(ijb)))

         ! diagnostic only
         ! shear = sqrt(tension**2 + shearing**2)
         shear(ijb) = tarear(ijb)*sqrt( tensionT(ijb)**2 + shearT(ijb)**2 )
         ! vorticity
         vort (ijb) = tarear(ijb)*( ( dyE(ijb)*vvelE(ijb) - dyE(ijb-1)*vvelE(ijb-1) ) &
                                  - ( dxN(ijb)*uvelN(ijb) - dxN(ijb-nx_block)*uvelN(ijb-nx_block)) )

      enddo                     ! ij

    end subroutine deformationsCD_T


!=======================================================================
! Compute deformations for mechanical redistribution at T point
!
! author: JF Lemieux, ECCC
! Nov 2021

    subroutine deformationsC_T (nx_block,   ny_block,   &
                                icellT,                 &
                                indxTi,     indxTj,     &
                                uvelE,      vvelE,      &
                                uvelN,      vvelN,      &
                                dxN,        dyE,        &
                                dxT,        dyT,        &
                                tarear,     uarea,      &
                                shearU,     vort,       &
                                shear,      divu,       &
                                rdg_conv,   rdg_shear )

      use ice_constants, only: p5

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellT                ! no. of cells where iceTmask = .true.

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxTi   , & ! compressed index in i-direction
         indxTj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(in) :: &
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

      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(inout) :: &
         vort     , & ! vorticity (1/s)
         shear    , & ! strain rate II component (1/s)
         divu     , & ! strain rate I component, velocity divergence (1/s)
         rdg_conv , & ! convergence term for ridging (1/s)
         rdg_shear    ! shear term for ridging (1/s)

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij       , & ! ij is the compressed cell counter
         ijb                ! 1D index of (i,j) in the block

      real (kind=dbl_kind), dimension (nx_block*ny_block) :: &
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
                           icellT     ,                &
                           indxTi(:)  , indxTj  (:)  , &
                           uvelE (:)  , vvelE   (:)  , &
                           uvelN (:)  , vvelN   (:)  , &
                           dxN   (:)  , dyE     (:)  , &
                           dxT   (:)  , dyT     (:)  , &
                           divT  (:)  , tensionT(:)  , &
                           shearT(:)  , DeltaT  (:)    )

      ! DeltaT is calc by strain_rates_T but replaced by calculation below.

      do ij = 1, icellT
         i = indxTi(ij)
         j = indxTj(ij)
         ijb = (j-1)*nx_block + i

         !-----------------------------------------------------------------
         ! deformations for mechanical redistribution
         !-----------------------------------------------------------------

         shearTsqr = (shearU(ijb)**2 * uarea(ijb)  &
                    + shearU(ijb-nx_block)**2 * uarea(ijb-nx_block)  &
                    + shearU(ijb-nx_block-1)**2 * uarea(ijb-nx_block-1)  &
                    + shearU(ijb-1)**2 * uarea(ijb-1)) &
                    / (uarea(ijb)+uarea(ijb-nx_block)+uarea(ijb-nx_block-1)+uarea(ijb-1))

         DeltaT(ijb) = sqrt(divT(ijb)**2 + e_factor*(tensionT(ijb)**2 + shearTsqr))

         divu(ijb) = divT(ijb) * tarear(ijb)
         tmp = DeltaT(ijb) * tarear(ijb)
         rdg_conv(ijb)  = -min(divu(ijb),c0)
         rdg_shear(ijb) = p5*(tmp-abs(divu(ijb)))

         ! diagnostic only...maybe we dont want to use shearTsqr here????
         ! shear = sqrt(tension**2 + shearing**2)
         shear(ijb) = tarear(ijb)*sqrt( tensionT(ijb)**2 + shearT(ijb)**2 )
         ! vorticity
         vort (ijb) = tarear(ijb)*( ( dyE(ijb)*vvelE(ijb) - dyE(ijb-1)*vvelE(ijb-1) ) &
                                  - ( dxN(ijb)*uvelN(ijb) - dxN(ijb-nx_block)*uvelN(ijb-nx_block)) )

      enddo                     ! ij

    end subroutine deformationsC_T

!=======================================================================
! Compute strain rates
!
! author: Elizabeth C. Hunke, LANL
!
! 2019: subroutine created by Philippe Blain, ECCC

      subroutine strain_rates (nx_block,   ny_block,   &
                               ijb,                    &
                               uvel00,     uvelW,      &
                               uvelS,      uvelSW,     &
                               vvel00,     vvelW,      &
                               vvelS,      vvelSW,     &
                               dxT,        dyT,        &
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
         ijb                   ! 1D index of (i,j) in the block

      ! the i and j dimensions are flattened: element (i,j) is stored at
      ! (j-1)*nx_block + i, so neighbours are ijb+-1 (i+-1), ijb+-nx_block (j+-1)
      ! The velocity stencil arrives as values, not as arrays plus an index:
      ! this routine is called both with block-layout and with compressed
      ! velocities, and only the caller knows how to reach a neighbour.
      ! The four are the cell itself and its W (-1,0), S (0,-1) and
      ! SW (-1,-1) neighbours.
      real (kind=dbl_kind), intent(in) :: &
         uvel00   , & ! x-velocity at (i  ,j  ) (m/s)
         uvelW    , & ! x-velocity at (i-1,j  )
         uvelS    , & ! x-velocity at (i  ,j-1)
         uvelSW   , & ! x-velocity at (i-1,j-1)
         vvel00   , & ! y-velocity at (i  ,j  )
         vvelW    , & ! y-velocity at (i-1,j  )
         vvelS    , & ! y-velocity at (i  ,j-1)
         vvelSW       ! y-velocity at (i-1,j-1)

      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(in) :: &
         dxT      , & ! width of T-cell through the middle (m)
         dyT      , & ! height of T-cell through the middle (m)
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
      divune    = cyp(ijb)*uvel00 - dyT(ijb)*uvelW &
                + cxp(ijb)*vvel00 - dxT(ijb)*vvelS
      divunw    = cym(ijb)*uvelW + dyT(ijb)*uvel00 &
                + cxp(ijb)*vvelW - dxT(ijb)*vvelSW
      divusw    = cym(ijb)*uvelSW + dyT(ijb)*uvelS &
                + cxm(ijb)*vvelSW + dxT(ijb)*vvelW
      divuse    = cyp(ijb)*uvelS - dyT(ijb)*uvelSW &
                + cxm(ijb)*vvelS + dxT(ijb)*vvel00

      ! tension strain rate  =  e_11 - e_22
      tensionne = -cym(ijb)*uvel00 - dyT(ijb)*uvelW &
                +  cxm(ijb)*vvel00 + dxT(ijb)*vvelS
      tensionnw = -cyp(ijb)*uvelW + dyT(ijb)*uvel00 &
                +  cxm(ijb)*vvelW + dxT(ijb)*vvelSW
      tensionsw = -cyp(ijb)*uvelSW + dyT(ijb)*uvelS &
                +  cxp(ijb)*vvelSW - dxT(ijb)*vvelW
      tensionse = -cym(ijb)*uvelS - dyT(ijb)*uvelSW &
                +  cxp(ijb)*vvelS - dxT(ijb)*vvel00

      ! shearing strain rate  =  2*e_12
      shearne = -cym(ijb)*vvel00 - dyT(ijb)*vvelW &
              -  cxm(ijb)*uvel00 - dxT(ijb)*uvelS
      shearnw = -cyp(ijb)*vvelW + dyT(ijb)*vvel00 &
              -  cxm(ijb)*uvelW - dxT(ijb)*uvelSW
      shearsw = -cyp(ijb)*vvelSW + dyT(ijb)*vvelS &
              -  cxp(ijb)*uvelSW + dxT(ijb)*uvelW
      shearse = -cym(ijb)*vvelS - dyT(ijb)*vvelSW &
              -  cxp(ijb)*uvelS + dxT(ijb)*uvel00

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
                                     icellT,               &
                                     indxTi,     indxTj,   &
                                     uvelE,      vvelE,    &
                                     uvelN,      vvelN,    &
                                     dxN,        dyE,      &
                                     dxT,        dyT,      &
                                     divT,       tensionT, &
                                     shearT,     DeltaT    )

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, &  ! block dimensions
         icellT

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxTi   , & ! compressed index in i-direction
         indxTj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(in) :: &
         uvelE    , & ! x-component of velocity (m/s) at the E point
         vvelE    , & ! y-component of velocity (m/s) at the N point
         uvelN    , & ! x-component of velocity (m/s) at the E point
         vvelN    , & ! y-component of velocity (m/s) at the N point
         dxN      , & ! width of N-cell through the middle (m)
         dyE      , & ! height of E-cell through the middle (m)
         dxT      , & ! width of T-cell through the middle (m)
         dyT          ! height of T-cell through the middle (m)

      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(out):: &
         divT     , & ! divergence at T point
         tensionT , & ! tension at T point
         shearT   , & ! shear at T point
         DeltaT       ! strain rates at the T point

      ! local variables

      integer (kind=int_kind) :: &
         ij, i, j               , & ! ij is the compressed cell counter
         ijb                        ! 1D index of (i,j) in the block

      character(len=*), parameter :: subname = '(strain_rates_Tdtsd)'

      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------

      ! compute divT, tensionT
      call strain_rates_Tdt (nx_block,   ny_block, &
                             icellT,               &
                             indxTi,     indxTj,   &
                             uvelE,      vvelE,    &
                             uvelN,      vvelN,    &
                             dxN,        dyE,      &
                             dxT,        dyT,      &
                             divT,       tensionT  )

      shearT  (:) = c0
      deltaT  (:) = c0

      do ij = 1, icellT
         i = indxTi(ij)
         j = indxTj(ij)
         ijb = (j-1)*nx_block + i

         ! shearing strain rate  =  2*e_12
         shearT(ijb) = (dxT(ijb)**2)*(uvelN(ijb)/dxN(ijb) - uvelN(ijb-nx_block)/dxN(ijb-nx_block)) &
                     + (dyT(ijb)**2)*(vvelE(ijb)/dyE(ijb) - vvelE(ijb-1)/dyE(ijb-1))

         ! Delta (in the denominator of zeta, eta)
         DeltaT(ijb) = sqrt(divT(ijb)**2 + e_factor*(tensionT(ijb)**2 + shearT(ijb)**2))

      enddo

      end subroutine strain_rates_Tdtsd

!=======================================================================
! Compute the dt (div, tension) strain rates at the T point
!
! author: JF Lemieux, ECCC
! Nov 2021

      subroutine strain_rates_Tdt (nx_block,   ny_block, &
                                   icellT,               &
                                   indxTi,     indxTj,   &
                                   uvelE,      vvelE,    &
                                   uvelN,      vvelN,    &
                                   dxN,        dyE,      &
                                   dxT,        dyT,      &
                                   divT,       tensionT  )

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, &  ! block dimensions
         icellT

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxTi   , & ! compressed index in i-direction
         indxTj       ! compressed index in j-direction

      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(in) :: &
         uvelE    , & ! x-component of velocity (m/s) at the E point
         vvelE    , & ! y-component of velocity (m/s) at the N point
         uvelN    , & ! x-component of velocity (m/s) at the E point
         vvelN    , & ! y-component of velocity (m/s) at the N point
         dxN      , & ! width of N-cell through the middle (m)
         dyE      , & ! height of E-cell through the middle (m)
         dxT      , & ! width of T-cell through the middle (m)
         dyT          ! height of T-cell through the middle (m)

      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(out):: &
         divT     , & ! divergence at T point
         tensionT     ! tension at T point

      ! local variables

      integer (kind=int_kind) :: &
         ij, i, j               , & ! ij is the compressed cell counter
         ijb                        ! 1D index of (i,j) in the block

      character(len=*), parameter :: subname = '(strain_rates_Tdt)'

      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------

      divT    (:) = c0
      tensionT(:) = c0

      do ij = 1, icellT
         i = indxTi(ij)
         j = indxTj(ij)
         ijb = (j-1)*nx_block + i

         ! divergence  =  e_11 + e_22
         divT    (ijb)= dyE(ijb)*uvelE(ijb) - dyE(ijb-1)*uvelE(ijb-1) &
                      + dxN(ijb)*vvelN(ijb) - dxN(ijb-nx_block)*vvelN(ijb-nx_block)

         ! tension strain rate  =  e_11 - e_22
         tensionT(ijb) = (dyT(ijb)**2)*(uvelE(ijb)/dyE(ijb) - uvelE(ijb-1)/dyE(ijb-1)) &
                       - (dxT(ijb)**2)*(vvelN(ijb)/dxN(ijb) - vvelN(ijb-nx_block)/dxN(ijb-nx_block))

      enddo

      end subroutine strain_rates_Tdt

!=======================================================================
! Compressed form of strain_rates_Tdt, for the C/CD subcycle.
!
! ADDITIVE, like grid_average_X2Y*_compressed: the block version above is
! still what deformationsC_T and deformationsCD_T call, and those run after
! the subcycle on block arrays, so converting them is not needed to put the
! subcycle on slots.  ⚠ The two share their arithmetic -- change one and the
! other must change with it, or C/CD stops being bit-for-bit with itself.
!
! The stencil is W on uvelE/dyE and S on vvelN/dxN, both already in dynNbr.
! dxT and dyT are read POINTWISE and so stay block-indexed; divT and
! tensionT are consumed pointwise by the caller and stay block-indexed too.

      subroutine strain_rates_Tdt_compressed (nx_block,   ny_block, &
                                   icellT,               &
                                   indxTi,     indxTj,   &
                                   nbr,        Xmap,     &
                                   uvelE_s,    vvelN_s,  &
                                   dxN_s,      dyE_s,    &
                                   dxT,        dyT,      &
                                   divT,       tensionT  )

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, &  ! block dimensions
         icellT

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxTi   , & ! compressed index in i-direction
         indxTj       ! compressed index in j-direction

      ! nbr(is,d) is the block-local slot one step in direction d from slot is
      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nbr          ! neighbour slots

      ! Xmap is assumed shape, 1D over the block in Fortran storage order,
      ! i.e. element (i,j) is stored at (j-1)*nx_block+i
      integer (kind=int_kind), dimension (:), intent(in) :: &
         Xmap         ! static compressed slot of each cell

      ! Compressed, block-local slots 1..cnt with the dead slot last.
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         uvelE_s  , & ! x-component of velocity at the E point, compressed
         vvelN_s  , & ! y-component of velocity at the N point, compressed
         dxN_s    , & ! width of N-cell through the middle (m), compressed
         dyE_s        ! height of E-cell through the middle (m), compressed

      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(in) :: &
         dxT      , & ! width of T-cell through the middle (m)
         dyT          ! height of T-cell through the middle (m)

      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(out):: &
         divT     , & ! divergence at T point
         tensionT     ! tension at T point

      ! local variables

      integer (kind=int_kind) :: &
         ij, i, j               , & ! ij is the compressed cell counter
         ijb                    , & ! 1D index of (i,j) in the block
         is, isW, isS               ! own, W and S neighbour slots

      character(len=*), parameter :: subname = '(strain_rates_Tdt_compressed)'

      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------

      divT    (:) = c0
      tensionT(:) = c0

      do ij = 1, icellT
         i = indxTi(ij)
         j = indxTj(ij)
         ijb = (j-1)*nx_block + i

         is  = Xmap(ijb)
         isW = nbr(is,nbr_W)
         isS = nbr(is,nbr_S)

         ! divergence  =  e_11 + e_22
         divT    (ijb)= dyE_s(is)*uvelE_s(is) - dyE_s(isW)*uvelE_s(isW) &
                      + dxN_s(is)*vvelN_s(is) - dxN_s(isS)*vvelN_s(isS)

         ! tension strain rate  =  e_11 - e_22
         tensionT(ijb) = (dyT(ijb)**2)*(uvelE_s(is)/dyE_s(is) - uvelE_s(isW)/dyE_s(isW)) &
                       - (dxT(ijb)**2)*(vvelN_s(is)/dxN_s(is) - vvelN_s(isS)/dxN_s(isS))

      enddo

      end subroutine strain_rates_Tdt_compressed

!=======================================================================
! Compressed form of strain_rates_Tdtsd, for the CD subcycle.  Same additive
! reasoning and the same ⚠ as strain_rates_Tdt_compressed above.
!
! Note the stencil is the MIRROR of Tdt's: Tdt reads uvelE to the W and vvelN
! to the S, this reads uvelN to the S and vvelE to the W.  Getting those the
! wrong way round is type-correct and silent, which is what kernel_equiv.py
! exists to catch.

      subroutine strain_rates_Tdtsd_compressed (nx_block, ny_block, &
                                   icellT,               &
                                   indxTi,     indxTj,   &
                                   nbr,        Xmap,     &
                                   uvelE_s,    vvelE_s,  &
                                   uvelN_s,    vvelN_s,  &
                                   dxN_s,      dyE_s,    &
                                   dxT,        dyT,      &
                                   divT,       tensionT, &
                                   shearT,     deltaT    )

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, &  ! block dimensions
         icellT

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxTi   , & ! compressed index in i-direction
         indxTj       ! compressed index in j-direction

      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nbr          ! neighbour slots
      integer (kind=int_kind), dimension (:), intent(in) :: &
         Xmap         ! static compressed slot of each cell

      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         uvelE_s  , & ! x-component of velocity (m/s) at the E point
         vvelE_s  , & ! y-component of velocity (m/s) at the E point
         uvelN_s  , & ! x-component of velocity (m/s) at the N point
         vvelN_s  , & ! y-component of velocity (m/s) at the N point
         dxN_s    , & ! width of N-cell through the middle (m)
         dyE_s        ! height of E-cell through the middle (m)

      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(in) :: &
         dxT      , & ! width of T-cell through the middle (m)
         dyT          ! height of T-cell through the middle (m)

      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(out):: &
         divT     , & ! divergence at T point
         tensionT , & ! tension at T point
         shearT   , & ! shear at T point
         deltaT       ! delt at T point

      ! local variables

      integer (kind=int_kind) :: &
         ij, i, j               , & ! ij is the compressed cell counter
         ijb                    , & ! 1D index of (i,j) in the block
         is, isW, isS               ! own, W and S neighbour slots

      character(len=*), parameter :: subname = '(strain_rates_Tdtsd_compressed)'

      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------

      ! compute divT, tensionT
      call strain_rates_Tdt_compressed (nx_block, ny_block, &
                             icellT,               &
                             indxTi,     indxTj,   &
                             nbr,        Xmap,     &
                             uvelE_s,    vvelN_s,  &
                             dxN_s,      dyE_s,    &
                             dxT,        dyT,      &
                             divT,       tensionT  )

      shearT  (:) = c0
      deltaT  (:) = c0

      do ij = 1, icellT
         i = indxTi(ij)
         j = indxTj(ij)
         ijb = (j-1)*nx_block + i

         is  = Xmap(ijb)
         isW = nbr(is,nbr_W)
         isS = nbr(is,nbr_S)

         ! shearing strain rate  =  2*e_12
         shearT(ijb) = (dxT(ijb)**2)*(uvelN_s(is)/dxN_s(is) - uvelN_s(isS)/dxN_s(isS)) &
                     + (dyT(ijb)**2)*(vvelE_s(is)/dyE_s(is) - vvelE_s(isW)/dyE_s(isW))

         ! Delta (in the denominator of zeta, eta)
         DeltaT(ijb) = sqrt(divT(ijb)**2 + e_factor*(tensionT(ijb)**2 + shearT(ijb)**2))

      enddo

      end subroutine strain_rates_Tdtsd_compressed

!=======================================================================
! Compute strain rates at the U point including boundary conditions
!
! author: JF Lemieux, ECCC
! Nov 2021

!=======================================================================
! Compressed form of strain_rates_U, for the C-grid subcycle.
!
! ADDITIVE for now: the CD subcycle still calls the block version above, and
! converting CD is its own tranche.  ⚠ The two share their arithmetic -- when
! CD is converted this one becomes the only version and the block one goes.
      subroutine strain_rates_U (nx_block, ny_block,  &
                                 icellU,                &
                                 indxUi,     indxUj,    &
                                 nbr,        Xmap,      &
                                 uvelE_s,    vvelE_s,   &
                                 uvelN_s,    vvelN_s,   &
                                 uvelU,      vvelU,     &
                                 dxE_s,      dyN_s,     &
                                 dxU,        dyU,       &
                                 ratiodxN,   ratiodxNr, &
                                 ratiodyE,   ratiodyEr, &
                                 epm_s,      npm_s,     &
                                 divergU,    tensionU,  &
                                 shearU_s,   DeltaU     )

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellU

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxUi   , & ! compressed index in i-direction
         indxUj       ! compressed index in j-direction

      ! nbr(is,d) is the block-local slot one step in direction d from slot is
      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nbr          ! neighbour slots

      ! Xmap is assumed shape, 1D over the block in Fortran storage order,
      ! i.e. element (i,j) is stored at (j-1)*nx_block+i
      integer (kind=int_kind), dimension (:), intent(in) :: &
         Xmap         ! static compressed slot of each cell

      ! Compressed: this kernel reaches E for the N-point quantities and N for
      ! the E-point ones.  epm and npm are real 0/1 masks used arithmetically.
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         uvelE_s  , & ! x-component of velocity (m/s) at the E point
         vvelE_s  , & ! y-component of velocity (m/s) at the N point
         uvelN_s  , & ! x-component of velocity (m/s) at the E point
         vvelN_s  , & ! y-component of velocity (m/s) at the N point
         dxE_s    , & ! width of E-cell through the middle (m)
         dyN_s    , & ! height of N-cell through the middle (m)
         epm_s    , & ! E-cell mask
         npm_s        ! N-cell mask

      ! uvelU and vvelU are COMPRESSED -- the driver has passed uvel_slot
      ! and vvel_slot since the velocities stopped leaving slot form inside
      ! the subcycle.  They were left declared block-indexed and read as
      ! uvelU(is); an explicit-shape dummy is not bounds-checked against its
      ! actual, so that read ran off the end of the slice in silence and
      ! DeltaU came out NaN.
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         uvelU    , & ! x-component of velocity (m/s) interp. at U point
         vvelU        ! y-component of velocity (m/s) interp. at U point

      ! Pointwise, so block-indexed.
      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(in) :: &
         dxU      , & ! width of U-cell through the middle (m)
         dyU      , & ! height of U-cell through the middle (m)
         ratiodxN , & ! -dxN(i+1,j)/dxN(i,j) for BCs
         ratiodxNr, & ! -dxN(i,j)/dxN(i+1,j) for BCs
         ratiodyE , & ! -dyE(i,j+1)/dyE(i,j) for BCs
         ratiodyEr    ! -dyE(i,j)/dyE(i,j+1) for BCs

      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(out):: &
         divergU  , & ! divergence at U point
         tensionU , & ! tension at U point
         DeltaU       ! delt at the U point

      ! shearU is read at W, S and SW by stressC_T, so it is compressed.
      real (kind=dbl_kind), dimension (:), intent(out), contiguous :: &
         shearU_s     ! shear at U point

      ! local variables

      integer (kind=int_kind) :: &
         ij, i, j               , & ! ij is the compressed cell counter
         ijb                    , & ! 1D index of (i,j) in the block
         is, isE, isN               ! own, E and N neighbour slots

      real (kind=dbl_kind) :: &
        uNip1j, uNij, vEijp1, vEij, uEijp1, uEij, vNip1j, vNij

      character(len=*), parameter :: subname = '(strain_rates_U)'

      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------

      divergU (:) = c0
      tensionU(:) = c0
      shearU_s(:) = c0
      deltaU  (:) = c0

      do ij = 1, icellU
         i = indxUi(ij)
         j = indxUj(ij)
         ijb = (j-1)*nx_block + i

         is  = Xmap(ijb)
         isE = nbr(is,nbr_E)   ! (i+1,j  )
         isN = nbr(is,nbr_N)   ! (i  ,j+1)

         uNip1j = uvelN_s(isE) * npm_s(isE) &
                +(npm_s(is)-npm_s(isE)) * npm_s(is)   * ratiodxN(ijb)  * uvelN_s(is)
         uNij   = uvelN_s(is) * npm_s(is) &
                +(npm_s(isE)-npm_s(is)) * npm_s(isE) * ratiodxNr(ijb) * uvelN_s(isE)
         vEijp1 = vvelE_s(isN) * epm_s(isN) &
                +(epm_s(is)-epm_s(isN)) * epm_s(is)   * ratiodyE(ijb)  * vvelE_s(is)
         vEij   = vvelE_s(is) * epm_s(is) &
                +(epm_s(isN)-epm_s(is)) * epm_s(isN) * ratiodyEr(ijb) * vvelE_s(isN)

         ! divergence  =  e_11 + e_22
         divergU (ijb) = dyU(ijb) * ( uNip1j - uNij ) &
                       + uvelU(is) * ( dyN_s(isE) - dyN_s(is) ) &
                       + dxU(ijb) * ( vEijp1 - vEij ) &
                       + vvelU(is) * ( dxE_s(isN) - dxE_s(is) )

         ! tension strain rate  =  e_11 - e_22
         tensionU(ijb) = dyU(ijb) * ( uNip1j - uNij ) &
                       - uvelU(is) * ( dyN_s(isE) - dyN_s(is) ) &
                       - dxU(ijb) * ( vEijp1 - vEij ) &
                       + vvelU(is) * ( dxE_s(isN) - dxE_s(is) )

         uEijp1 = uvelE_s(isN) * epm_s(isN) &
                +(epm_s(is)-epm_s(isN)) * epm_s(is)   * ratiodyE(ijb)  * uvelE_s(is)
         uEij   = uvelE_s(is) * epm_s(is) &
                +(epm_s(isN)-epm_s(is)) * epm_s(isN) * ratiodyEr(ijb) * uvelE_s(isN)
         vNip1j = vvelN_s(isE) * npm_s(isE) &
                +(npm_s(is)-npm_s(isE)) * npm_s(is)   * ratiodxN(ijb)  * vvelN_s(is)
         vNij   = vvelN_s(is) * npm_s(is) &
                +(npm_s(isE)-npm_s(is)) * npm_s(isE) * ratiodxNr(ijb) * vvelN_s(isE)

         ! shearing strain rate  =  2*e_12
         shearU_s(is)   = dxU(ijb) * ( uEijp1 - uEij ) &
                       - uvelU(is) * ( dxE_s(isN) - dxE_s(is) ) &
                       + dyU(ijb) * ( vNip1j - vNij ) &
                       - vvelU(is) * ( dyN_s(isE) - dyN_s(is) )

         ! Delta (in the denominator of zeta, eta)
         DeltaU(ijb)   = sqrt(divergU(ijb)**2 + e_factor*(tensionU(ijb)**2 + shearU_s(is)**2))

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
                                zetax2, etax2, rep_prs)

      real (kind=dbl_kind), intent(in)::  &
         strength, & !
         DminArea    !

      real (kind=dbl_kind), intent(in)::  &
         Delta

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
      use ice_domain, only: maskhalo_dyn
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
      use ice_domain, only: maskhalo_dyn, halo_dynbundle
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
      use ice_domain, only: maskhalo_dyn, halo_dynbundle
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
      use ice_domain, only: maskhalo_dyn, halo_dynbundle
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
      use ice_domain, only: maskhalo_dyn, halo_dynbundle
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
! Load fields into array for boundary updates, where the fields have their
! i and j dimensions flattened, i.e. dimension(nx_block*ny_block,max_blocks),
! with element (i,j,iblk) at ((j-1)*nx_block + i, iblk)

      subroutine stack_fields3_flat(fld1, fld2, fld3, fldbundle)

      use ice_domain, only: nblocks
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bundbound

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         fld1    , & ! flattened fields to stack
         fld2    , & !
         fld3        !

      real (kind=dbl_kind), dimension (:,:,:,:), intent(out) :: &
         fldbundle   ! work array for boundary updates (i,j,n,iblk)

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij, & ! horizontal indices, ij is the 1D index in the block
         iblk        ! block index

      character(len=*), parameter :: subname = '(stack_fields3_flat)'

      call ice_timer_start(timer_bundbound)
      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ij) SCHEDULE(runtime)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            ij = (j-1)*nx_block + i
            fldbundle(i,j,1,iblk) = fld1(ij,iblk)
            fldbundle(i,j,2,iblk) = fld2(ij,iblk)
            fldbundle(i,j,3,iblk) = fld3(ij,iblk)
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      call ice_timer_stop(timer_bundbound)

      end subroutine stack_fields3_flat

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
! Unload fields from array after boundary updates, where the fields have
! their i and j dimensions flattened, i.e.
! dimension(nx_block*ny_block,max_blocks), with element (i,j,iblk) at
! ((j-1)*nx_block + i, iblk)

      subroutine unstack_fields3_flat(fldbundle, fld1, fld2, fld3)

      use ice_domain, only: nblocks
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bundbound

      real (kind=dbl_kind), dimension (:,:,:,:), intent(in) :: &
         fldbundle   ! work array for boundary updates (i,j,n,iblk)

      real (kind=dbl_kind), dimension (:,:), intent(out) :: &
         fld1    , & ! flattened fields to unstack
         fld2    , & !
         fld3        !

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij, & ! horizontal indices, ij is the 1D index in the block
         iblk        ! block index

      character(len=*), parameter :: subname = '(unstack_fields3_flat)'

      call ice_timer_start(timer_bundbound)
      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ij) SCHEDULE(runtime)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            ij = (j-1)*nx_block + i
            fld1(ij,iblk) = fldbundle(i,j,1,iblk)
            fld2(ij,iblk) = fldbundle(i,j,2,iblk)
            fld3(ij,iblk) = fldbundle(i,j,3,iblk)
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      call ice_timer_stop(timer_bundbound)

      end subroutine unstack_fields3_flat

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
! Do a halo update on 1 flattened field

      subroutine dyn_haloUpdate1_flat(halo_info, halo_info_mask, field_loc, field_type, fld1)

      use ice_boundary, only: ice_halo, ice_HaloUpdate
      use ice_domain, only: maskhalo_dyn
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bound

      type (ice_halo), intent(in) :: &
         halo_info  , &  ! standard unmasked halo
         halo_info_mask  ! masked halo

      integer (kind=int_kind), intent(in) :: &
         field_loc  ,  & ! field loc
         field_type      ! field_type

      ! flattened: element (i,j,iblk) lives at ((j-1)*nx_block + i, iblk)
      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         fld1            !

      ! local variables

      character(len=*), parameter :: subname = '(dyn_haloUpdate1_flat)'

      call ice_timer_start(timer_bound)
         if (maskhalo_dyn) then
            call ice_HaloUpdate (fld1     , halo_info_mask, &
                                 field_loc, field_type)
         else
            call ice_HaloUpdate (fld1     , halo_info     , &
                                 field_loc, field_type)
         endif
      call ice_timer_stop(timer_bound)

      end subroutine dyn_haloUpdate1_flat

!=======================================================================
! Do a halo update on 2 flattened fields

      subroutine dyn_haloUpdate2_flat(halo_info, halo_info_mask, field_loc, field_type, fld1, fld2)

      use ice_boundary, only: ice_halo, ice_HaloUpdate
      use ice_domain, only: maskhalo_dyn, halo_dynbundle
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bound

      type (ice_halo), intent(in) :: &
         halo_info  , &  ! standard unmasked halo
         halo_info_mask  ! masked halo

      integer (kind=int_kind), intent(in) :: &
         field_loc  ,  & ! field loc
         field_type      ! field_type

      ! flattened: element (i,j,iblk) lives at ((j-1)*nx_block + i, iblk)
      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         fld1        , & ! fields to halo
         fld2            !

      ! local variables

      real (kind=dbl_kind), dimension (nx_block,ny_block,2,max_blocks) :: &
         fldbundle       ! work array for boundary updates

      character(len=*), parameter :: subname = '(dyn_haloUpdate2_flat)'

      call ice_timer_start(timer_bound)
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

      end subroutine dyn_haloUpdate2_flat

!=======================================================================
! Do a halo update on 3 flattened fields

      subroutine dyn_haloUpdate3_flat(halo_info, halo_info_mask, field_loc, field_type, fld1, fld2, fld3)

      use ice_boundary, only: ice_halo, ice_HaloUpdate
      use ice_domain, only: maskhalo_dyn, halo_dynbundle
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bound

      type (ice_halo), intent(in) :: &
         halo_info  , &  ! standard unmasked halo
         halo_info_mask  ! masked halo

      integer (kind=int_kind), intent(in) :: &
         field_loc  ,  & ! field loc
         field_type      ! field_type

      ! flattened: element (i,j,iblk) lives at ((j-1)*nx_block + i, iblk)
      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         fld1        , & ! fields to halo
         fld2        , & ! fields to halo
         fld3            !

      ! local variables

      real (kind=dbl_kind), dimension (nx_block,ny_block,3,max_blocks) :: &
         fldbundle       ! work array for boundary updates

      character(len=*), parameter :: subname = '(dyn_haloUpdate3_flat)'

      call ice_timer_start(timer_bound)
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

      end subroutine dyn_haloUpdate3_flat

!=======================================================================
! Do a halo update on 4 flattened fields

      subroutine dyn_haloUpdate4_flat(halo_info, halo_info_mask, field_loc, field_type, fld1, fld2, fld3, fld4)

      use ice_boundary, only: ice_halo, ice_HaloUpdate
      use ice_domain, only: maskhalo_dyn, halo_dynbundle
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bound

      type (ice_halo), intent(in) :: &
         halo_info  , &  ! standard unmasked halo
         halo_info_mask  ! masked halo

      integer (kind=int_kind), intent(in) :: &
         field_loc  ,  & ! field loc
         field_type      ! field_type

      ! flattened: element (i,j,iblk) lives at ((j-1)*nx_block + i, iblk)
      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         fld1        , & ! fields to halo
         fld2        , & ! fields to halo
         fld3        , & ! fields to halo
         fld4            !

      ! local variables

      real (kind=dbl_kind), dimension (nx_block,ny_block,4,max_blocks) :: &
         fldbundle       ! work array for boundary updates

      character(len=*), parameter :: subname = '(dyn_haloUpdate4_flat)'

      call ice_timer_start(timer_bound)
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

      end subroutine dyn_haloUpdate4_flat

!=======================================================================
! Load 2 flattened fields into a bundle

!=======================================================================
! Bundle two COMPRESSED fields so they cross the halo in one message.
!
! The compressed counterpart of stack_fields2.  Both operands and the bundle
! are indexed by compressed slot, and the whole extent is copied including
! slot 0 -- each block's dead slot -- so the dead slots stay zero in the
! bundle and the stencils keep reading zero for a missing neighbour.
!
! No OpenMP here: this is a straight whole-array copy over the compressed
! extent, which has no block dimension to thread over.

      subroutine stack_fields3_compressed(fld1, fld2, fld3, fldbundle)

      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bundbound

      ! compressed: indexed by slot; each block's dead slot is its last
      real (kind=dbl_kind), dimension (:), intent(in) :: &
         fld1, fld2, fld3

      ! (slot, field)
      real (kind=dbl_kind), dimension (:,:), intent(out) :: &
         fldbundle   ! work array for boundary updates

      character(len=*), parameter :: subname = '(stack_fields3_compressed)'

      call ice_timer_start(timer_bundbound)
      fldbundle(:,1) = fld1(:)
      fldbundle(:,2) = fld2(:)
      fldbundle(:,3) = fld3(:)
      call ice_timer_stop(timer_bundbound)

      end subroutine stack_fields3_compressed

!=======================================================================
! Unbundle three COMPRESSED fields.  See stack_fields3_compressed.

      subroutine unstack_fields3_compressed(fldbundle, fld1, fld2, fld3)

      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bundbound

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         fldbundle   ! work array for boundary updates

      real (kind=dbl_kind), dimension (:), intent(out) :: &
         fld1, fld2, fld3

      character(len=*), parameter :: subname = '(unstack_fields3_compressed)'

      call ice_timer_start(timer_bundbound)
      fld1(:) = fldbundle(:,1)
      fld2(:) = fldbundle(:,2)
      fld3(:) = fldbundle(:,3)
      call ice_timer_stop(timer_bundbound)

      end subroutine unstack_fields3_compressed

!=======================================================================

      subroutine stack_fields4_compressed(fld1, fld2, fld3, fld4, fldbundle)

      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bundbound

      ! compressed: indexed by slot; each block's dead slot is its last
      real (kind=dbl_kind), dimension (:), intent(in) :: &
         fld1, fld2, fld3, fld4

      ! (slot, field)
      real (kind=dbl_kind), dimension (:,:), intent(out) :: &
         fldbundle   ! work array for boundary updates

      character(len=*), parameter :: subname = '(stack_fields4_compressed)'

      call ice_timer_start(timer_bundbound)
      fldbundle(:,1) = fld1(:)
      fldbundle(:,2) = fld2(:)
      fldbundle(:,3) = fld3(:)
      fldbundle(:,4) = fld4(:)
      call ice_timer_stop(timer_bundbound)

      end subroutine stack_fields4_compressed

!=======================================================================
! Unbundle four COMPRESSED fields.  See stack_fields4_compressed.

      subroutine unstack_fields4_compressed(fldbundle, fld1, fld2, fld3, fld4)

      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bundbound

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         fldbundle   ! work array for boundary updates

      real (kind=dbl_kind), dimension (:), intent(out) :: &
         fld1, fld2, fld3, fld4

      character(len=*), parameter :: subname = '(unstack_fields4_compressed)'

      call ice_timer_start(timer_bundbound)
      fld1(:) = fldbundle(:,1)
      fld2(:) = fldbundle(:,2)
      fld3(:) = fldbundle(:,3)
      fld4(:) = fldbundle(:,4)
      call ice_timer_stop(timer_bundbound)

      end subroutine unstack_fields4_compressed

!=======================================================================

      subroutine stack_fields2_compressed(fld1, fld2, fldbundle)

      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bundbound

      ! compressed: indexed by slot; each block's dead slot is its last
      real (kind=dbl_kind), dimension (:), intent(in) :: &
         fld1    , & !
         fld2        !

      ! (slot, field)
      real (kind=dbl_kind), dimension (:,:), intent(out) :: &
         fldbundle   ! work array for boundary updates

      character(len=*), parameter :: subname = '(stack_fields2_compressed)'

      call ice_timer_start(timer_bundbound)
      fldbundle(:,1) = fld1(:)
      fldbundle(:,2) = fld2(:)
      call ice_timer_stop(timer_bundbound)

      end subroutine stack_fields2_compressed

!=======================================================================
! Unbundle two COMPRESSED fields.  See stack_fields2_compressed.

      subroutine unstack_fields2_compressed(fldbundle, fld1, fld2)

      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bundbound

      ! (slot, field)
      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         fldbundle   ! work array for boundary updates

      ! compressed: indexed by slot; each block's dead slot is its last
      real (kind=dbl_kind), dimension (:), intent(out) :: &
         fld1    , & !
         fld2        !

      character(len=*), parameter :: subname = '(unstack_fields2_compressed)'

      call ice_timer_start(timer_bundbound)
      fld1(:) = fldbundle(:,1)
      fld2(:) = fldbundle(:,2)
      call ice_timer_stop(timer_bundbound)

      end subroutine unstack_fields2_compressed

!=======================================================================

      subroutine stack_fields2_flat(fld1, fld2, fldbundle)

      use ice_domain, only: nblocks
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bundbound

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         fld1    , & !
         fld2        !

      real (kind=dbl_kind), dimension (:,:,:,:), intent(out) :: &
         fldbundle   ! work array for boundary updates (i,j,n,iblk)

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij, & ! horizontal indices, ij is the 1D index in the block
         iblk        ! block index

      character(len=*), parameter :: subname = '(stack_fields2_flat)'

      call ice_timer_start(timer_bundbound)
      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ij) SCHEDULE(runtime)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            ij = (j-1)*nx_block + i
            fldbundle(i,j,1,iblk) = fld1(ij,iblk)
            fldbundle(i,j,2,iblk) = fld2(ij,iblk)
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      call ice_timer_stop(timer_bundbound)

      end subroutine stack_fields2_flat

!=======================================================================
! Load 4 flattened fields into a bundle

      subroutine stack_fields4_flat(fld1, fld2, fld3, fld4, fldbundle)

      use ice_domain, only: nblocks
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bundbound

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         fld1    , & !
         fld2    , & !
         fld3    , & !
         fld4        !

      real (kind=dbl_kind), dimension (:,:,:,:), intent(out) :: &
         fldbundle   ! work array for boundary updates (i,j,n,iblk)

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij, & ! horizontal indices, ij is the 1D index in the block
         iblk        ! block index

      character(len=*), parameter :: subname = '(stack_fields4_flat)'

      call ice_timer_start(timer_bundbound)
      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ij) SCHEDULE(runtime)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            ij = (j-1)*nx_block + i
            fldbundle(i,j,1,iblk) = fld1(ij,iblk)
            fldbundle(i,j,2,iblk) = fld2(ij,iblk)
            fldbundle(i,j,3,iblk) = fld3(ij,iblk)
            fldbundle(i,j,4,iblk) = fld4(ij,iblk)
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      call ice_timer_stop(timer_bundbound)

      end subroutine stack_fields4_flat

!=======================================================================
! Unload 2 flattened fields from a bundle

      subroutine unstack_fields2_flat(fldbundle, fld1, fld2)

      use ice_domain, only: nblocks
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bundbound

      real (kind=dbl_kind), dimension (:,:,:,:), intent(in) :: &
         fldbundle   ! work array for boundary updates (i,j,n,iblk)

      real (kind=dbl_kind), dimension (:,:), intent(out) :: &
         fld1    , & !
         fld2        !

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij, & ! horizontal indices, ij is the 1D index in the block
         iblk        ! block index

      character(len=*), parameter :: subname = '(unstack_fields2_flat)'

      call ice_timer_start(timer_bundbound)
      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ij) SCHEDULE(runtime)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            ij = (j-1)*nx_block + i
            fld1(ij,iblk) = fldbundle(i,j,1,iblk)
            fld2(ij,iblk) = fldbundle(i,j,2,iblk)
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      call ice_timer_stop(timer_bundbound)

      end subroutine unstack_fields2_flat

!=======================================================================
! Unload 4 flattened fields from a bundle

      subroutine unstack_fields4_flat(fldbundle, fld1, fld2, fld3, fld4)

      use ice_domain, only: nblocks
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bundbound

      real (kind=dbl_kind), dimension (:,:,:,:), intent(in) :: &
         fldbundle   ! work array for boundary updates (i,j,n,iblk)

      real (kind=dbl_kind), dimension (:,:), intent(out) :: &
         fld1    , & !
         fld2    , & !
         fld3    , & !
         fld4        !

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij, & ! horizontal indices, ij is the 1D index in the block
         iblk        ! block index

      character(len=*), parameter :: subname = '(unstack_fields4_flat)'

      call ice_timer_start(timer_bundbound)
      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ij) SCHEDULE(runtime)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            ij = (j-1)*nx_block + i
            fld1(ij,iblk) = fldbundle(i,j,1,iblk)
            fld2(ij,iblk) = fldbundle(i,j,2,iblk)
            fld3(ij,iblk) = fldbundle(i,j,3,iblk)
            fld4(ij,iblk) = fldbundle(i,j,4,iblk)
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      call ice_timer_stop(timer_bundbound)

      end subroutine unstack_fields4_flat

!=======================================================================

      end module ice_dyn_shared

!=======================================================================
