!=======================================================================
!
! Viscous-plastic sea ice dynamics model
! Computes ice velocity and deformation
!
! See:
!
! Lemieux, J.‐F., Tremblay, B., Thomas, S., Sedláček, J., and Mysak, L. A. (2008),
! Using the preconditioned Generalized Minimum RESidual (GMRES) method to solve
! the sea‐ice momentum equation, J. Geophys. Res., 113, C10004, doi:10.1029/2007JC004680.
!
! Hibler, W. D., and Ackley, S. F. (1983), Numerical simulation of the Weddell Sea pack ice,
! J. Geophys. Res., 88( C5), 2873– 2887, doi:10.1029/JC088iC05p02873.
!
! Y. Saad. A Flexible Inner-Outer Preconditioned GMRES Algorithm. SIAM J. Sci. Comput.,
! 14(2):461–469, 1993. URL: https://doi.org/10.1137/0914028, doi:10.1137/0914028.
!
! C. T. Kelley, Iterative Methods for Linear and Nonlinear Equations, SIAM, 1995.
! (https://www.siam.org/books/textbooks/fr16_book.pdf)
!
! Y. Saad, Iterative Methods for Sparse Linear Systems. SIAM, 2003.
! (http://www-users.cs.umn.edu/~saad/IterMethBook_2ndEd.pdf)
!
! Walker, H. F., & Ni, P. (2011). Anderson Acceleration for Fixed-Point Iterations.
! SIAM Journal on Numerical Analysis, 49(4), 1715–1735. https://doi.org/10.1137/10078356X
!
! Fang, H., & Saad, Y. (2009). Two classes of multisecant methods for nonlinear acceleration.
! Numerical Linear Algebra with Applications, 16(3), 197–221. https://doi.org/10.1002/nla.617
!
! Birken, P. (2015) Termination criteria for inexact fixed‐point schemes.
! Numer. Linear Algebra Appl., 22: 702– 716. doi: 10.1002/nla.1982.
!
! authors: JF Lemieux, ECCC, Philppe Blain, ECCC
!

      module ice_dyn_vp

      use ice_kinds_mod
      use ice_blocks, only: nx_block, ny_block
      use ice_boundary, only: ice_halo
      use ice_communicate, only: my_task, master_task, get_num_procs
      use ice_constants, only: field_loc_center, field_loc_NEcorner, &
          field_type_scalar, field_type_vector
      use ice_constants, only: c0, p027, p055, p111, p166, &
          p222, p25, p333, p5, c1
      use ice_domain, only: nblocks, distrb_info
      use ice_domain_size, only: max_blocks
      use ice_dyn_shared, only: dyn_prep1, dyn_prep2, dyn_finish, &
          cosw, sinw, fcor_blk, uvel_init, vvel_init, dynOff, dynCnt, dynMap, dynIJ, dynMaxCnt, skipT, skipU, &
          dynTot, pack_to_slots, unpack_from_slots, calc_seabed_stress_slot, &
          dynNbr, nbr_W, nbr_S, nbr_SW, nbr_E, nbr_N, nbr_NE, uarear_slot, &
          dxT_slot, dyT_slot, dxhy_slot, dyhx_slot, cxp_slot, cyp_slot, &
          uvel_slot, vvel_slot, tarear_slot, dxU_slot, dyU_slot, &
          cxm_slot, cym_slot, DminTarea_slot, strength_slot, &
          seabed_stress_factor_LKD, seabed_stress_factor_prob, seabed_stress_method, &
          seabed_stress, Ktens, stack_fields,  unstack_fields, fld2, fld3, fld4, fld2c
      use ice_fileunits, only: nu_diag
      use ice_flux, only: fmU
      use ice_global_reductions, only: global_sum
      use ice_grid, only: dxT, dyT, uarear
      use ice_exit, only: abort_ice
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_ice_strength, icepack_query_parameters

      implicit none
      private
      public :: implicit_solver, init_vp

      ! namelist parameters

      integer (kind=int_kind), public :: &
         maxits_nonlin  , & ! max nb of iteration for nonlinear solver
         dim_fgmres     , & ! size of fgmres Krylov subspace
         dim_pgmres     , & ! size of pgmres Krylov subspace
         maxits_fgmres  , & ! max nb of iteration for fgmres
         maxits_pgmres  , & ! max nb of iteration for pgmres
         fpfunc_andacc  , & ! fixed point function for Anderson acceleration:
                            ! 1: g(x) = FMGRES(A(x),b(x)), 2: g(x) = x - A(x)x + b(x)
         dim_andacc     , & ! size of Anderson minimization matrix (number of saved previous residuals)
         start_andacc       ! acceleration delay factor (acceleration starts at this iteration)

      ! VP's compressed working copies, indexed by slot over the whole run
      ! (dynTot).  Stage 2 converts the solver kernels to iterate slots the way
      ! EVP's do; until every consumer of a field is converted, the block-layout
      ! copy stays and the pack/unpack pair at the call site keeps the two in
      ! step.  Those pairs are transitional, not the destination.
      real (kind=dbl_kind), dimension(:), allocatable :: &
         Cw_s     , & ! ocean-ice neutral drag coefficient, compressed
         aiU_s    , & ! ice fraction on u-grid, compressed
         TbU_s    , & ! seabed stress factor (N/m^2), compressed
         uocn_s   , & ! ocean current, x-direction (m/s), compressed
         vocn_s   , & ! ocean current, y-direction (m/s), compressed
         ulin_s   , & ! uvel used to linearize vrel, compressed
         vlin_s   , & ! vvel used to linearize vrel, compressed
         vrel_s   , & ! coeff for tauw, compressed
         Cb_s     , & ! seabed stress coefficient, compressed
         taubx_s  , & ! seabed stress, x-direction (N/m^2), compressed
         tauby_s      ! seabed stress, y-direction (N/m^2), compressed

      ! The four-corner and eight-component T-cell fields.  Kept as (slot,n)
      ! rather than split into separate vectors the way EVP's str1..str8 are:
      ! those are read through the neighbour table with a different component
      ! at each neighbour, so splitting them separates the streams.  These are
      ! read with all components at one slot, which is the same number of
      ! streams either way, and (slot,n) keeps the slot index stride-1.
      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         zetax2_s , & ! zetax2 = 2*zeta (bulk viscosity), compressed
         etax2_s  , & ! etax2  = 2*eta  (shear viscosity), compressed
         repprs_s , & ! replacement pressure, compressed
         stPr_s       ! stress combinations from replacement pressure, compressed

      real (kind=dbl_kind), dimension(:), allocatable :: &
         waterx_s , & ! for ocean stress calculation, x (m/s), compressed
         watery_s , & ! for ocean stress calculation, y (m/s), compressed
         bxfix_s  , & ! part of bx constant during Picard, compressed
         byfix_s  , & ! part of by constant during Picard, compressed
         bx_s     , & ! b vector, x, compressed
         by_s     , & ! b vector, y, compressed
         mvu_s    , & ! matvec input velocity, x, compressed
         mvv_s    , & ! matvec input velocity, y, compressed
         mvAu_s   , & ! matvec output, x, compressed
         mvAv_s   , & ! matvec output, y, compressed
         umassdti_s,& ! mass of U-cell/dt (kg/m^2 s), compressed
         fmU_s    , & ! Coriolis param. * mass in U-cell (kg/s), compressed
         ! Krylov-side vectors.  Fx/Fy never leave the slot world; solx/soly,
         ! diagx/diagy and the fixed-point residual are packed per nonlinear
         ! iteration because arrays_to_vec still reads them block-shaped.
         Fx_s     , & ! x residual vector, compressed
         Fy_s     , & ! y residual vector, compressed
         fpx_s    , & ! fixed point residual, x, compressed
         fpy_s    , & ! fixed point residual, y, compressed
         diagx_s  , & ! diagonal of the system matrix, x, compressed
         diagy_s  , & ! diagonal of the system matrix, y, compressed
         solx_s   , & ! FGMRES solution, x, compressed
         soly_s       ! FGMRES solution, y, compressed

      logical (kind=log_kind), public :: &
         monitor_nonlin , & ! print nonlinear residual norm
         monitor_fgmres , & ! print fgmres residual norm
         monitor_pgmres , & ! print pgmres residual norm
         use_mean_vrel      ! use mean of previous 2 iterates to compute vrel (see Hibler and Ackley 1983)

      real (kind=dbl_kind), public :: &
         reltol_nonlin  , & ! nonlinear stopping criterion: reltol_nonlin*res(k=0)
         reltol_fgmres  , & ! fgmres stopping criterion: reltol_fgmres*res(k)
         reltol_pgmres  , & ! pgmres stopping criterion: reltol_pgmres*res(k)
         damping_andacc , & ! damping factor for Anderson acceleration
         reltol_andacc      ! relative tolerance for Anderson acceleration

      character (len=char_len), public :: &
         precond        , & ! preconditioner for fgmres: 'ident' (identity), 'diag' (diagonal),
                            ! 'pgmres' (Jacobi-preconditioned GMRES)
         algo_nonlin    , & ! nonlinear algorithm: 'picard' (Picard iteration), 'anderson' (Anderson acceleration)
         ortho_type         ! type of orthogonalization for FGMRES ('cgs' or 'mgs')

      ! module variables

      integer (kind=int_kind), allocatable :: &
         icellT(:)    , & ! no. of cells where iceTmask = .true.
         icellU(:)        ! no. of cells where iceUmask = .true.

      integer (kind=int_kind), allocatable :: &
         indxTi(:,:)  , & ! compressed index in i-direction
         indxTj(:,:)  , & ! compressed index in j-direction
         indxUi(:,:)  , & ! compressed index in i-direction
         indxUj(:,:)      ! compressed index in j-direction

!=======================================================================

      contains

!=======================================================================

! Initialize parameters and variables needed for the vp dynamics
! author: Philippe Blain, ECCC

      subroutine init_vp

      use ice_blocks, only: get_block, block
      use ice_boundary, only: ice_HaloUpdate
      use ice_calendar, only: dt_dyn
      use ice_dyn_shared, only: init_dyn_shared

      call init_dyn_shared(dt_dyn)

      ! Initialize module variables
      allocate(icellT(max_blocks), icellU(max_blocks))
      allocate(indxTi(nx_block*ny_block, max_blocks), &
               indxTj(nx_block*ny_block, max_blocks), &
               indxUi(nx_block*ny_block, max_blocks), &
               indxUj(nx_block*ny_block, max_blocks))

      ! init_dyn_shared above has built the compressed layout, so dynTot is set
      allocate(Cw_s  (dynTot), aiU_s (dynTot), TbU_s (dynTot), &
               uocn_s(dynTot), vocn_s(dynTot), ulin_s(dynTot), &
               vlin_s(dynTot), vrel_s(dynTot), Cb_s  (dynTot), &
               taubx_s(dynTot), tauby_s(dynTot))
      allocate(zetax2_s(dynTot,4), etax2_s(dynTot,4), &
               repprs_s(dynTot,4), stPr_s (dynTot,8))
      zetax2_s = c0 ; etax2_s = c0 ; repprs_s = c0 ; stPr_s = c0
      allocate(waterx_s(dynTot), watery_s(dynTot), bxfix_s(dynTot), &
               byfix_s (dynTot), bx_s    (dynTot), by_s   (dynTot))
      waterx_s = c0 ; watery_s = c0 ; bxfix_s = c0
      byfix_s  = c0 ; bx_s     = c0 ; by_s    = c0
      allocate(mvu_s(dynTot), mvv_s(dynTot), mvAu_s(dynTot), mvAv_s(dynTot), &
               umassdti_s(dynTot), fmU_s(dynTot))
      mvu_s = c0 ; mvv_s = c0 ; mvAu_s = c0 ; mvAv_s = c0
      umassdti_s = c0 ; fmU_s = c0
      ! Krylov-side vectors.  solx/soly and diagx/diagy are packed per
      ! nonlinear iteration because arrays_to_vec and the history still read
      ! them in the block layout; Fx/Fy never leave the slot world.
      allocate(Fx_s(dynTot), Fy_s(dynTot), fpx_s(dynTot), fpy_s(dynTot), &
               diagx_s(dynTot), diagy_s(dynTot), solx_s(dynTot), soly_s(dynTot))
      Fx_s = c0 ; Fy_s = c0 ; fpx_s = c0 ; fpy_s = c0
      diagx_s = c0 ; diagy_s = c0 ; solx_s = c0 ; soly_s = c0
      Cw_s   = c0 ; aiU_s  = c0 ; TbU_s = c0
      uocn_s = c0 ; vocn_s = c0 ; ulin_s = c0
      vlin_s = c0 ; vrel_s = c0 ; Cb_s  = c0
      taubx_s = c0 ; tauby_s = c0

      end subroutine init_vp

!=======================================================================

! Viscous-plastic dynamics driver
!
#ifdef CICE_IN_NEMO
! Wind stress is set during this routine from the values supplied
! via NEMO (unless calc_strair is true).  These values are supplied
! rotated on u grid and multiplied by aice.  strairxT = 0 in this
! case so operations in dyn_prep1 are pointless but carried out to
! minimise code changes.
#endif
!
! author: JF Lemieux, A. Qaddouri and F. Dupont ECCC

      subroutine implicit_solver (dt)

      use ice_arrays_column, only: Cdn_ocn
      use ice_boundary, only: ice_HaloMask, ice_HaloUpdate, &
          ice_HaloDestroy, ice_HaloUpdate_stress
      use ice_blocks, only: block, get_block, nx_block, ny_block, &
          ns_boundary_type
      use ice_domain, only: blocks_ice, halo_info, maskhalo_dyn
      use ice_domain_size, only: max_blocks
      use ice_dyn_shared, only: deformations, iceTmask, iceUmask, &
          cxp, cyp, cxm, cym
      use ice_flux, only: rdg_conv, rdg_shear, strairxT, strairyT,                    &
          strairxU, strairyU, uocn, vocn, ss_tltx, ss_tlty, fmU,                      &
          strtltxU, strtltyU, strocnxU, strocnyU, strintxU, strintyU, taubxU, taubyU, &
          strax, stray,                                                               &
          hwater,                                                                &
          stressp_1, stressp_2, stressp_3, stressp_4,                                 &
          stressm_1, stressm_2, stressm_3, stressm_4,                                 &
          stress12_1, stress12_2, stress12_3, stress12_4
      use ice_grid, only: tmask, umask, dxT, dyT, dxU, dyU, &
          tarear, grid_average_X2Y,                         &
          grid_atm_dynu, grid_atm_dynv, grid_ocn_dynu, grid_ocn_dynv
      use ice_state, only: aice, aiU, vice, vsno, uvel, vvel, divu, shear, vort, &
          aice_init, aice0, aicen, vicen, strength
      use ice_timers, only: timer_dynamics, timer_bound, &
          ice_timer_start, ice_timer_stop

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) :: &
         ntot           , & ! size of problem for Anderson
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, ij

      ! Flattened, as the group below already is: element (i,j,iblk) lives at
      ! ((j-1)*nx_block + i, iblk).  grid_average_X2Y resolves to its 3d-in,
      ! flat-out specific, and dyn_prep2 takes the slice directly.
      real (kind=dbl_kind), dimension (nx_block*ny_block,max_blocks) :: &
         ss_tltxU , & ! sea surface slope, x-direction (m/m)
         ss_tltyU     ! sea surface slope, y-direction (m/m)

      ! flattened like the group below: element (i,j,iblk) at
      ! ((j-1)*nx_block + i, iblk)
      real (kind=dbl_kind), dimension (nx_block*ny_block,max_blocks) :: &
         uocnU    , & ! i ocean current (m/s)
         vocnU    , & ! j ocean current (m/s)
         cdn_ocnU , & ! ocn drag coefficient
         bxfix    , & ! part of bx that is constant during Picard
         byfix    , & ! part of by that is constant during Picard
         Cb           ! seabed stress coefficient

      ! the i and j dimensions are flattened: element (i,j,iblk) of the block
      ! layout lives at ((j-1)*nx_block + i, iblk)
      real (kind=dbl_kind), dimension (nx_block*ny_block,max_blocks) :: &
         waterxU  , & ! for ocean stress calculation, x (m/s)
         wateryU  , & ! for ocean stress calculation, y (m/s)
         forcexU  , & ! work array: combined atm stress and ocn tilt, x
         forceyU      ! work array: combined atm stress and ocn tilt, y

      ! the i and j dimensions are flattened: element (i,j,iblk) of the block
      ! layout lives at ((j-1)*nx_block + i, iblk)
      real (kind=dbl_kind), dimension (nx_block*ny_block,max_blocks), target :: &
         umass    , & ! total mass of ice and snow (u grid)
         umassdti     ! mass of U-cell/dte (kg/m^2 s)

      ! the i and j dimensions are flattened: element (i,j,iblk) of the block
      ! layout lives at ((j-1)*nx_block + i, iblk)
      real (kind=dbl_kind), dimension (nx_block*ny_block,max_blocks), target :: &
         tmass        ! total mass of ice and snow (kg/m^2)

      logical (kind=log_kind) :: calc_strair

      integer (kind=int_kind), dimension (nx_block,ny_block,max_blocks) :: &
         halomask     ! generic halo mask

      type (ice_halo) :: &
         halo_info_mask !  ghost cell update info for masked halo

      type (block) :: &
         this_block           ! block information for current block

      real (kind=dbl_kind), allocatable :: &
         sol(:)          ! solution vector

      character(len=*), parameter :: subname = '(implicit_solver)'

      call ice_timer_start(timer_dynamics) ! dynamics

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

       ! This call is needed only if dt changes during runtime.
!      call set_evp_parameters (dt)

      !-----------------------------------------------------------------
      ! boundary updates
      ! commented out because the ghost cells are freshly
      ! updated after cleanup_itd
      !-----------------------------------------------------------------

!      call ice_timer_start(timer_bound)
!      call ice_HaloUpdate (aice,              halo_info, &
!                           field_loc_center,  field_type_scalar)
!      call ice_HaloUpdate (vice,              halo_info, &
!                           field_loc_center,  field_type_scalar)
!      call ice_HaloUpdate (vsno,              halo_info, &
!                           field_loc_center,  field_type_scalar)
!      call ice_timer_stop(timer_bound)

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block)
      do iblk = 1, nblocks

         do j = 1, ny_block
         do i = 1, nx_block
            rdg_conv ((j-1)*nx_block+i,iblk) = c0
            rdg_shear((j-1)*nx_block+i,iblk) = c0
            divu ((j-1)*nx_block+i,iblk) = c0
            shear((j-1)*nx_block+i,iblk) = c0
         enddo
         enddo

      !-----------------------------------------------------------------
      ! preparation for dynamics
      !-----------------------------------------------------------------

         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         call dyn_prep1 (nx_block,           ny_block,           &
                         ilo, ihi,           jlo, jhi,           &
                         aice    (:,:,iblk), vice    (:,:,iblk), &
                         vsno    (:,:,iblk), tmask   (:,  iblk), &
                         tmass   (:,iblk),   iceTmask(:,iblk))

      enddo                     ! iblk
      !$OMP END PARALLEL DO

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (iceTmask,          halo_info, &
                           field_loc_center,  field_type_scalar)
      call ice_timer_stop(timer_bound)

      !-----------------------------------------------------------------
      ! convert fields from T to U grid
      !-----------------------------------------------------------------

      call stack_fields(tmass, aice_init, cdn_ocn, fld3)
      call ice_HaloUpdate (fld3,             halo_info, &
                           field_loc_center, field_type_scalar)
      call stack_fields(uocn, vocn, ss_tltx, ss_tlty, fld4)
      call ice_HaloUpdate (fld4,             halo_info, &
                           field_loc_center, field_type_vector)
      call unstack_fields(fld3, tmass, aice_init, cdn_ocn)
      call unstack_fields(fld4, uocn, vocn, ss_tltx, ss_tlty)

      call grid_average_X2Y('S',tmass    , 'T'          , umass   , 'U')
      call grid_average_X2Y('S',aice_init, 'T'          , aiU     , 'U')
      call grid_average_X2Y('S',cdn_ocn  , 'T'          , cdn_ocnU, 'U')
      call grid_average_X2Y('S',uocn     , grid_ocn_dynu, uocnU   , 'U')
      call grid_average_X2Y('S',vocn     , grid_ocn_dynv, vocnU   , 'U')
      call grid_average_X2Y('S',ss_tltx  , grid_ocn_dynu, ss_tltxU, 'U')
      call grid_average_X2Y('S',ss_tlty  , grid_ocn_dynv, ss_tltyU, 'U')

      !----------------------------------------------------------------
      ! Set wind stress to values supplied via NEMO or other forcing
      ! This wind stress is rotated on u grid and multiplied by aice
      !----------------------------------------------------------------
      call icepack_query_parameters(calc_strair_out=calc_strair)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (.not. calc_strair) then
         call grid_average_X2Y('F', strax, grid_atm_dynu, strairxU, 'U')
         call grid_average_X2Y('F', stray, grid_atm_dynv, strairyU, 'U')
      else
         call ice_HaloUpdate (strairxT,         halo_info, &
                              field_loc_center, field_type_vector)
         call ice_HaloUpdate (strairyT,         halo_info, &
                              field_loc_center, field_type_vector)
         call grid_average_X2Y('F',strairxT,'T',strairxU,'U')
         call grid_average_X2Y('F',strairyT,'T',strairyU,'U')
      endif

      !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block,ij,i,j)
      do iblk = 1, nblocks

      !-----------------------------------------------------------------
      ! more preparation for dynamics
      !-----------------------------------------------------------------

         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         call dyn_prep2 (nx_block,             ny_block,                                                                               &
                         ilo, ihi,             jlo, jhi,                                                                               &
                         icellT(iblk),         icellU(iblk),                                                                           &
                         indxTi      (:,iblk), indxTj      (:,iblk),                                                                   &
                         indxUi      (:,iblk), indxUj      (:,iblk),                                                                   &
                         aiU       (:,iblk),   umass     (:,iblk),                                                                     &
                         umassdti  (:,iblk),   fcor_blk  (:,iblk)  ,                                                                   &
                         umask     (:,iblk)  ,                                                                                         &
                         uocnU     (:,iblk), vocnU     (:,iblk),                                                                       &
                         strairxU  (:,:,iblk), strairyU  (:,:,iblk),                                                                   &
                         ss_tltxU  (:,  iblk), ss_tltyU  (:,  iblk),                                                                   &
                         iceTmask  (:,iblk),   iceUmask  (:,iblk),                                                                     &
                         fmU       (:,iblk), dt,                                                                                       &
                         strtltxU  (:,:,iblk), strtltyU  (:,:,iblk),                                                                   &
                         strocnxU  (:,:,iblk), strocnyU  (:,:,iblk),                                                                   &
                         strintxU  (:,iblk), strintyU  (:,iblk),                                                                   &
                         taubxU    (:,iblk)  , taubyU    (:,iblk)  ,                                                                   &
                         waterxU   (:,iblk)  , wateryU   (:,iblk)  ,                                                                   &
                         forcexU   (:,iblk)  , forceyU   (:,iblk)  ,                                                                   &
                         stressp_1 (:,iblk), stressp_2 (:,iblk),                                                                       &
                         stressp_3 (:,iblk), stressp_4 (:,iblk),                                                                       &
                         stressm_1 (:,iblk), stressm_2 (:,iblk),                                                                       &
                         stressm_3 (:,iblk), stressm_4 (:,iblk),                                                                       &
                         stress12_1(:,iblk), stress12_2(:,iblk),                                                                       &
                         stress12_3(:,iblk), stress12_4(:,iblk),                                                                       &
                         uvel_init (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),   vvel_init (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)), &
                         uvel(:,iblk), vvel(:,iblk),                                                                                   &
                        dynMap(:,iblk), skipU(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),                                               &
                        skipT(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),                                                               &
                         TbU_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))

         call calc_bfix (nx_block            , ny_block            ,           &
                         icellU(iblk)        ,                                 &
                         indxUi      (:,iblk), indxUj      (:,iblk),           &
                         umassdti  (:,iblk),                                   &
                         forcexU   (:,iblk)  , forceyU   (:,iblk)  ,           &
                         uvel_init (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)), &
                         vvel_init (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)), &
                         dynMap      (:,iblk),                                 &
                         bxfix     (:,iblk), byfix     (:,iblk))

      !-----------------------------------------------------------------
      ! ice strength
      !-----------------------------------------------------------------

         strength(:,:,iblk) = c0  ! initialize
         do ij = 1, icellT(iblk)
            i = indxTi(ij, iblk)
            j = indxTj(ij, iblk)
            call icepack_ice_strength (aice  = aice       (i,j,  iblk), &
                                       vice  = vice       (i,j,  iblk), &
                                       aice0 = aice0      (i,j,  iblk), &
                                       aicen = aicen      (i,j,:,iblk), &
                                       vicen = vicen      (i,j,:,iblk), &
                                       strength = strength(i,j,  iblk))
         enddo  ! ij

      enddo  ! iblk
      !$OMP END PARALLEL DO

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (strength,           halo_info, &
                           field_loc_center,   field_type_scalar)
      ! velocities may have changed in dyn_prep2
      call stack_fields(uvel, vvel, fld2)
      call ice_HaloUpdate (fld2,               halo_info, &
                           field_loc_NEcorner, field_type_vector)
      call unstack_fields(fld2, uvel, vvel)
      call ice_timer_stop(timer_bound)

      if (maskhalo_dyn) then
         call ice_timer_start(timer_bound)
         halomask = 0
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               if (iceUmask((j-1)*nx_block+i,iblk)) halomask(i,j,iblk) = 1
            enddo
            enddo
         enddo
         call ice_HaloUpdate (halomask,          halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_timer_stop(timer_bound)
         call ice_HaloMask(halo_info_mask, halo_info, halomask)
      endif

      !-----------------------------------------------------------------
      ! seabed stress factor TbU (TbU is part of Cb coefficient)
      !-----------------------------------------------------------------
      if (seabed_stress) then
         if ( seabed_stress_method == 'LKD' ) then
            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks
               call seabed_stress_factor_LKD (nx_block,         ny_block,       &
                                              icellU  (iblk),                   &
                                              indxUi(:,iblk),   indxUj(:,iblk), &
                                              vice(:,:,iblk),   aice(:,:,iblk), &
                                              hwater(:,iblk), dynMap(:,iblk), TbU_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            enddo
            !$OMP END PARALLEL DO

         elseif ( seabed_stress_method == 'probabilistic' ) then
            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks

               call seabed_stress_factor_prob (nx_block,         ny_block,                   &
                                               icellT(iblk), indxTi(:,iblk), indxTj(:,iblk), &
                                               icellU(iblk), indxUi(:,iblk), indxUj(:,iblk), &
                                               aicen(:,:,:,iblk), vicen(:,:,:,iblk),         &
                                               hwater(:,iblk), dynMap(:,iblk), TbU_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            enddo
            !$OMP END PARALLEL DO

         endif
      endif


      !-----------------------------------------------------------------
      ! calc size of problem (ntot) and allocate solution vector
      !-----------------------------------------------------------------

      ntot = 0
      do iblk = 1, nblocks
         ntot = ntot + icellU(iblk)
      enddo
      ntot = 2 * ntot ! times 2 because of u and v

      allocate(sol(ntot))

      !-----------------------------------------------------------------
      ! Start of nonlinear iteration
      !-----------------------------------------------------------------
      call anderson_solver (icellU , &
                            indxUi  , indxUj , &
                            aiU     , ntot   , &
                            uocnU   , vocnU  , &
                            waterxU , wateryU, &
                            bxfix   , byfix  , &
                            umassdti, sol    , &
                            cdn_ocnU,          &
                            Cb, halo_info_mask)
      !-----------------------------------------------------------------
      ! End of nonlinear iteration
      !-----------------------------------------------------------------

      deallocate(sol)

      if (maskhalo_dyn) call ice_HaloDestroy(halo_info_mask)

      !-----------------------------------------------------------------
      ! Compute stresses
      !-----------------------------------------------------------------
      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks
         ! The nonlinear iteration has finished, so ulin_s and vlin_s are free
         ! to carry the converged velocity into the slot kernel.
         call pack_to_slots(uvel(:,iblk), dynMap(:,iblk), ulin_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(vvel(:,iblk), dynMap(:,iblk), vlin_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))

         call stress_vp (dynCnt(iblk), skipT(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)), &
                         ulin_s  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),          &
                         vlin_s  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),          &
                         dxT_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),          &
                         dyT_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),          &
                         cxp_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),          &
                         cyp_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),          &
                         cxm_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),          &
                         cym_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),          &
                         dynNbr  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk),:),          &
                         dynIJ   (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),            &
                         zetax2_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:),        &
                         etax2_s (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:),        &
                         repprs_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:),        &
                         stressp_1 (:,iblk), stressp_2 (:,iblk),                        &
                         stressp_3 (:,iblk), stressp_4 (:,iblk),                        &
                         stressm_1 (:,iblk), stressm_2 (:,iblk),                        &
                         stressm_3 (:,iblk), stressm_4 (:,iblk),                        &
                         stress12_1(:,iblk), stress12_2(:,iblk),                        &
                         stress12_3(:,iblk), stress12_4(:,iblk))
      enddo ! iblk
      !$OMP END PARALLEL DO

      !-----------------------------------------------------------------
      ! Compute deformations
      !-----------------------------------------------------------------
      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks
         ! VP's velocity is block-resident at this point (vec_to_arrays wrote
         ! it), so it is packed for the kernel.  EVP keeps uvel_slot current
         ! across its subcycles and needs no such pack.
         call pack_to_slots(uvel(:,iblk), dynMap(:,iblk), uvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(vvel(:,iblk), dynMap(:,iblk), vvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call deformations (                                                                                                              &
                            dynCnt(iblk), skipT(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),                                                &
                            dynNbr(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk),:), dynIJ(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),          &
                            uvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), vvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                            dxT_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dyT_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),   &
                            dxU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dyU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),   &
                            cxp_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), cyp_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),   &
                            cxm_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), cym_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),   &
                            tarear_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), vort(:,iblk),                                      &
                            shear     (:,iblk), divu      (:,iblk),                                                                   &
                            rdg_conv  (:,iblk), rdg_shear (:,iblk))
      enddo
      !$OMP END PARALLEL DO

      !-----------------------------------------------------------------
      ! Compute seabed stress (diagnostic)
      !-----------------------------------------------------------------
      if (seabed_stress) then
         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks
            ! The shared kernel EVP already uses.  It carries Cb in the
            ! taub array -- reads it, overwrites it with the stress -- so
            ! taubx_s is seeded with Cb rather than passed alongside it.
            call pack_to_slots(uvel(:,iblk), dynMap(:,iblk), ulin_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(vvel(:,iblk), dynMap(:,iblk), vlin_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            taubx_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1) = &
               Cb_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1)

            call calc_seabed_stress_slot(dynCnt(iblk), skipU(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)), &
                 ulin_s (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),                                   &
                 taubx_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),                                   &
                 vlin_s (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),                                   &
                 tauby_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))

            call unpack_from_slots(taubx_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), taubxU(:,iblk))
            call unpack_from_slots(tauby_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), taubyU(:,iblk))
         enddo
         !$OMP END PARALLEL DO
      endif

      ! Force symmetry across the tripole seam
      if (trim(ns_boundary_type) == 'tripole' .or. &
          trim(ns_boundary_type) == 'tripoleT') then
      if (maskhalo_dyn) then
         !-------------------------------------------------------
         ! set halomask to zero because ice_HaloMask always keeps
         ! local copies AND tripole zipper communication
         !-------------------------------------------------------
         halomask = 0
         call ice_HaloMask(halo_info_mask, halo_info, halomask)

         call ice_HaloUpdate_stress(stressp_1, stressp_3, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_3, stressp_1, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_2, stressp_4, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_4, stressp_2, halo_info_mask, &
                              field_loc_center,  field_type_scalar)

         call ice_HaloUpdate_stress(stressm_1, stressm_3, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_3, stressm_1, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_2, stressm_4, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_4, stressm_2, halo_info_mask, &
                              field_loc_center,  field_type_scalar)

         call ice_HaloUpdate_stress(stress12_1, stress12_3, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_3, stress12_1, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_2, stress12_4, halo_info_mask, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_4, stress12_2, halo_info_mask, &
                              field_loc_center,  field_type_scalar)

         call ice_HaloDestroy(halo_info_mask)
      else
         call ice_HaloUpdate_stress(stressp_1, stressp_3, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_3, stressp_1, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_2, stressp_4, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressp_4, stressp_2, halo_info, &
                              field_loc_center,  field_type_scalar)

         call ice_HaloUpdate_stress(stressm_1, stressm_3, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_3, stressm_1, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_2, stressm_4, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stressm_4, stressm_2, halo_info, &
                              field_loc_center,  field_type_scalar)

         call ice_HaloUpdate_stress(stress12_1, stress12_3, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_3, stress12_1, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_2, stress12_4, halo_info, &
                              field_loc_center,  field_type_scalar)
         call ice_HaloUpdate_stress(stress12_4, stress12_2, halo_info, &
                              field_loc_center,  field_type_scalar)
      endif   ! maskhalo
      endif   ! tripole

      !-----------------------------------------------------------------
      ! ice-ocean stress
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks

         call dyn_finish                                &
              (nx_block,           ny_block,            &
               icellU      (iblk), cdn_ocnU(:,iblk),    &
               indxUi    (:,iblk), indxUj    (:,iblk),  &
               uvel(:,iblk), vvel(:,iblk),              &
               uocnU   (:,iblk), vocnU   (:,iblk),      &
               aiU     (:,iblk),   fmU     (:,iblk),    &
!               strintxU(:,iblk), strintyU(:,iblk), &
!               strairxU(:,:,iblk), strairyU(:,:,iblk), &
               strocnxU(:,:,iblk), strocnyU(:,:,iblk))

      enddo
      !$OMP END PARALLEL DO

! shift velocity components from CD grid locations (N, E) to B grid location (U) for transport
! commented out in order to focus on EVP for now within the cdgrid
! should be used when routine is ready
!      if (grid_ice == 'CD' .or. grid_ice == 'C') then
!          call grid_average_X2Y('E2US',uvelE,uvel)
!          call grid_average_X2Y('N2US',vvelN,vvel)
!      endif
!end comment out
      call ice_timer_stop(timer_dynamics)    ! dynamics

      end subroutine implicit_solver

!=======================================================================

! Solve the nonlinear equation F(u,v) = 0, where
! F(u,v) := A(u,v) * (u,v) - b(u,v)
! using Anderson acceleration (accelerated fixed point (Picard) iteration)
!
! author: JF Lemieux, A. Qaddouri, F. Dupont and P. Blain ECCC
!
! Anderson algorithm adadpted from:
! H. F. Walker, “Anderson Acceleration: Algorithms and Implementations”
!   [Online]. Available: https://users.wpi.edu/~walker/Papers/anderson_accn_algs_imps.pdf

      subroutine anderson_solver (icellU , &
                                  indxUi  , indxUj , &
                                  aiU     , ntot   , &
                                  uocn    , vocn   , &
                                  waterxU , wateryU, &
                                  bxfix   , byfix  , &
                                  umassdti, sol    , &
                                  cdn_ocn ,          &
                                  Cb, halo_info_mask)

      use ice_blocks, only: nx_block, ny_block
      use ice_boundary, only: ice_HaloUpdate
      use ice_constants, only: c1
      use ice_domain, only: maskhalo_dyn, halo_info
      use ice_domain_size, only: max_blocks
      use ice_flux, only:   fmU
      use ice_grid, only: dxT, dyT, uarear
      use ice_dyn_shared, only: DminTarea, dxhy, dyhx, cxp, cyp, cxm, cym
      use ice_state, only: uvel, vvel, strength
      use ice_restoring, only: ice_restoring_halo
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bound

      integer (kind=int_kind), intent(in) :: &
         ntot         ! size of problem for Anderson

      integer (kind=int_kind), dimension(max_blocks), intent(in) :: &
         icellU       ! no. of cells where iceUmask = .true.

      integer (kind=int_kind), dimension (nx_block*ny_block, max_blocks), intent(in) :: &
         indxUi   , & ! compressed index in i-direction
         indxUj       ! compressed index in j-direction

      ! aiU and umassdti are assumed shape, with the i and j dimensions
      ! flattened, i.e. element (i,j,iblk) lives at ((j-1)*nx_block + i, iblk)
      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         aiU      , & ! ice fraction on u-grid
         umassdti     ! mass of U-cell/dte (kg/m^2 s)

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         waterxU  , & ! for ocean stress calculation, x (m/s)
         wateryU      ! for ocean stress calculation, y (m/s)

      real (kind=dbl_kind), dimension (nx_block*ny_block,max_blocks), intent(in) :: &
         uocn     , & ! i ocean current (m/s)
         vocn     , & ! j ocean current (m/s)
         cdn_ocn  , & ! ocn drag coefficient
         bxfix    , & ! part of bx that is constant during Picard
         byfix        ! part of by that is constant during Picard

      type (ice_halo), intent(in) :: &
         halo_info_mask !  ghost cell update info for masked halo

      real (kind=dbl_kind), dimension (nx_block*ny_block,max_blocks), intent(inout) :: &
         Cb           ! seabed stress coefficient

      real (kind=dbl_kind), dimension (ntot), intent(inout) :: &
         sol          ! current approximate solution

      ! local variables

      integer (kind=int_kind) :: &
         it_nl      , & ! nonlinear loop iteration index
         res_num    , & ! current number of stored residuals
         j          , & ! iteration index for QR update
         iblk       , & ! block index
         nbiter         ! number of FGMRES iterations performed

      integer (kind=int_kind), parameter :: &
         inc = 1        ! increment value for BLAS calls

      ! these four hold copies of uvel/vvel, so they carry its flattened
      ! shape: element (i,j,iblk) at ((j-1)*nx_block + i, iblk)
      real (kind=dbl_kind), dimension (nx_block*ny_block,max_blocks) :: &
         uprev_k  , & ! uvel at previous Picard iteration
         vprev_k  , & ! vvel at previous Picard iteration
         ulin     , & ! uvel to linearize vrel
         vlin         ! vvel to linearize vrel

      real (kind=dbl_kind), dimension (nx_block*ny_block,max_blocks) :: &
         vrel     , & ! coeff for tauw
         bx       , & ! b vector
         by       , & ! b vector
         Au       , & ! matvec, Fx = bx - Au
         Av       , & ! matvec, Fy = by - Av
         Fx       , & ! x residual vector, Fx = bx - Au
         Fy           ! y residual vector, Fy = by - Av

      ! Scratch for the progress-residual norm only: uvel - uprev_k is formed
      ! here and packed for global_norm.  Local, because nothing outside this
      ! routine ever reads it.
      real (kind=dbl_kind), dimension (nx_block*ny_block,max_blocks) :: &
         fpresx   , & ! fixed point residual, x: fx = uvel - uprev_k
         fpresy       ! fixed point residual, y: fy = vvel - vprev_k

      ! OMP-private scratch that never leaves the formDiag pair, so it is
      ! sized for the largest block rather than dynTot: step1 fills it one
      ! block at a time.  The +1 is the block's dead slot.
      real (kind=dbl_kind), dimension (dynMaxCnt+1,8):: &
         diag_rheo      ! contributions of the rhelogy term to the diagonal

      real (kind=dbl_kind), dimension (ntot) :: &
         res        , & ! current residual
         res_old    , & ! previous residual
         res_diff   , & ! difference between current and previous residuals
         fpfunc     , & ! current value of fixed point function
         fpfunc_old , & ! previous value of fixed point function
         tmp            ! temporary vector for BLAS calls

      real (kind=dbl_kind), dimension(ntot,dim_andacc) :: &
         Q        , & ! Q factor for QR factorization of F (residuals) matrix
         G_diff       ! Matrix containing the differences of g(x) (fixed point function) evaluations

      real (kind=dbl_kind), dimension(dim_andacc,dim_andacc) :: &
         R            ! R factor for QR factorization of F (residuals) matrix

      real (kind=dbl_kind), dimension(dim_andacc) :: &
         rhs_tri  , & ! right hand side vector for matrix-vector product
         coeffs       ! coeffs used to combine previous solutions

      real (kind=dbl_kind) :: &
         ! tol         , & ! tolerance for fixed point convergence: reltol_andacc * (initial fixed point residual norm)  [unused for now]
         tol_nl      , & ! tolerance for nonlinear convergence: reltol_nonlin * (initial nonlinear residual norm)
         fpres_norm  , & ! norm of current fixed point residual : f(x) = g(x) - x
         prog_norm   , & ! norm of difference between current and previous solution
         nlres_norm      ! norm of current nonlinear residual : F(x) = A(x)x -b(x)

#ifdef USE_LAPACK
      real (kind=dbl_kind) :: &
         ddot, dnrm2     ! external BLAS functions
#endif

      character(len=*), parameter :: subname = '(anderson_solver)'

      ! Initialization
      res_num = 0

      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks
         uprev_k(:,iblk) = uvel(:,iblk)
         vprev_k(:,iblk) = vvel(:,iblk)
      enddo
      !$OMP END PARALLEL DO

      ! Start iterations
      ! umassdti and fmU do not change during the solve, so they are compressed
      ! once here rather than at each of the five matvec call sites -- matvec
      ! is called from fgmres and pgmres as well as from the loop below.
      !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
      do iblk = 1, nblocks
         if (dynCnt(iblk) == 0) cycle
         call pack_to_slots(umassdti(:,iblk), dynMap(:,iblk), umassdti_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call pack_to_slots(fmU     (:,iblk), dynMap(:,iblk), fmU_s     (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
      enddo
      !$OMP END PARALLEL DO

      do it_nl = 0, maxits_nonlin        ! nonlinear iteration loop
         ! Compute quantities needed for computing PDE residual and g(x) (fixed point map)
         !-----------------------------------------------------------------
         ! Calc zetax2, etax2, dPr/dx, dPr/dy, Cb and vrel = f(uprev_k, vprev_k)
         !-----------------------------------------------------------------
         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks

            if (use_mean_vrel) then
               ulin(:,iblk) = p5*uprev_k(:,iblk) + p5*uvel(:,iblk)
               vlin(:,iblk) = p5*vprev_k(:,iblk) + p5*vvel(:,iblk)
            else
               ulin(:,iblk) = uvel(:,iblk)
               vlin(:,iblk) = vvel(:,iblk)
            endif
            uprev_k(:,iblk) = uvel(:,iblk)
            vprev_k(:,iblk) = vvel(:,iblk)

            ! The T-cell metrics are packed once per run in ice_dyn_shared
            ! and ice_grid, so only what changes per iteration is packed here.
            call pack_to_slots(uprev_k (:,iblk), dynMap(:,iblk), ulin_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(vprev_k (:,iblk), dynMap(:,iblk), vlin_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(strength(:,:,iblk), dynMap(:,iblk), strength_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))

            call calc_zeta_dPr (dynCnt(iblk), skipT(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)), &
                                ulin_s        (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),    &
                                vlin_s        (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),    &
                                dxT_slot      (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),    &
                                dyT_slot      (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),    &
                                dxhy_slot     (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),    &
                                dyhx_slot     (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),    &
                                cxp_slot      (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),    &
                                cyp_slot      (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),    &
                                cxm_slot      (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),    &
                                cym_slot      (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),    &
                                DminTarea_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),    &
                                strength_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),    &
                                dynNbr        (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk),:),    &
                                zetax2_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:),        &
                                etax2_s (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:),        &
                                repprs_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:),        &
                                stPr_s  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:))

            ! Transitional pack/unpack: every other consumer of vrel and Cb
            ! still reads the block layout, so the compressed copies are built
            ! here and written back.  The pair disappears as those consumers
            ! convert.
            call pack_to_slots(cdn_ocn(:,iblk), dynMap(:,iblk), Cw_s  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(aiU    (:,iblk), dynMap(:,iblk), aiU_s (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(uocn   (:,iblk), dynMap(:,iblk), uocn_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(vocn   (:,iblk), dynMap(:,iblk), vocn_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(ulin   (:,iblk), dynMap(:,iblk), ulin_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(vlin   (:,iblk), dynMap(:,iblk), vlin_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))

            call calc_vrel_Cb (dynCnt(iblk), skipU (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)), &
                               Cw_s  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),             &
                               aiU_s (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),             &
                               TbU_s (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),             &
                               uocn_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),             &
                               vocn_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),             &
                               ulin_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),             &
                               vlin_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),             &
                               vrel_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),             &
                               Cb_s  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))

            call unpack_from_slots(vrel_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), vrel(:,iblk))
            call unpack_from_slots(Cb_s  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), Cb  (:,iblk))

            ! prepare b vector (RHS)
            ! stPr and vrel are already compressed; uarear is compressed once
            ! per run.  Only the four still-block-layout inputs are packed.
            call pack_to_slots(waterxU(:,iblk), dynMap(:,iblk), waterx_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(wateryU(:,iblk), dynMap(:,iblk), watery_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(bxfix  (:,iblk), dynMap(:,iblk), bxfix_s (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(byfix  (:,iblk), dynMap(:,iblk), byfix_s (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))

            call calc_bvec (dynCnt(iblk), skipU(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)), &
                            stPr_s   (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:),       &
                            uarear_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),       &
                            waterx_s (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),         &
                            watery_s (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),         &
                            bxfix_s  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),         &
                            byfix_s  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),         &
                            bx_s     (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),         &
                            by_s     (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),         &
                            vrel_s   (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),         &
                            dynNbr   (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk),:))

            call unpack_from_slots(bx_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), bx(:,iblk))
            call unpack_from_slots(by_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), by(:,iblk))

            ! Compute nonlinear residual norm (PDE residual)
            call pack_to_slots(uprev_k(:,iblk), dynMap(:,iblk), mvu_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(vprev_k(:,iblk), dynMap(:,iblk), mvv_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call matvec (dynCnt(iblk), skipT(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),                                                &
                         skipU(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)), dynNbr(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk),:),          &
                         dxT_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dyT_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),   &
                         dxhy_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dyhx_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                         cxp_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), cyp_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),   &
                         cxm_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), cym_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),   &
                         mvu_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), mvv_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),         &
                         vrel_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), Cb_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),         &
                         zetax2_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:),                                                       &
                         etax2_s (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:),                                                       &
                         umassdti_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), fmU_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),    &
                         uarear_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),                                                      &
                         mvAu_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), mvAv_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call residual_vec (dynCnt(iblk), skipU(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),                                          &
                               bx_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), by_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),     &
                               mvAu_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), mvAv_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               Fx_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), Fy_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         enddo
         !$OMP END PARALLEL DO
         nlres_norm = global_norm(Fx_s, Fy_s)
         if (my_task == master_task .and. monitor_nonlin) then
            write(nu_diag, '(a,i4,a,d26.16)') "monitor_nonlin: iter_nonlin= ", it_nl, &
                                              " nonlin_res_L2norm= ", nlres_norm
         endif
         ! Compute relative tolerance at first iteration
         if (it_nl == 0) then
            tol_nl = reltol_nonlin*nlres_norm
         endif

         ! Check for nonlinear convergence
         if (nlres_norm < tol_nl) then
            exit
         endif

         ! Put initial guess for FGMRES in solx_s,soly_s and sol (needed for
         ! anderson).  The guess is uprev_k, so it is packed once here and
         ! then stays compressed through fgmres and arrays_to_vec alike.
         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks
            call pack_to_slots(uprev_k(:,iblk), dynMap(:,iblk), solx_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(vprev_k(:,iblk), dynMap(:,iblk), soly_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         enddo
         !$OMP END PARALLEL DO
         call arrays_to_vec (ntot, solx_s, soly_s, sol(:))

         ! Compute fixed point map g(x)
         if (fpfunc_andacc == 1) then
            ! g_1(x) = FGMRES(A(x), b(x))

            ! Prepare diagonal for preconditioner
            if (precond == 'diag' .or. precond == 'pgmres') then
               !$OMP PARALLEL DO PRIVATE(iblk,diag_rheo)
               do iblk = 1, nblocks
                  ! first compute diagonal contributions due to rheology term
                  call formDiag_step1 (dynCnt(iblk), skipU(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)), &
                                       dxT_slot   (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),       &
                                       dyT_slot   (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),       &
                                       dxhy_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),       &
                                       dyhx_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),       &
                                       cxp_slot   (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),       &
                                       cyp_slot   (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),       &
                                       cxm_slot   (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),       &
                                       cym_slot   (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),       &
                                       dynNbr     (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk),:),       &
                                       zetax2_s   (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:),     &
                                       etax2_s    (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:),     &
                                       diag_rheo(:,:))
                  ! second compute the full diagonal
                  call formDiag_step2 (dynCnt(iblk), skipU(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)), &
                                       diag_rheo(:,:),                                                &
                                       vrel_s     (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),       &
                                       umassdti_s (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),       &
                                       uarear_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),       &
                                       Cb_s       (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),       &
                                       diagx_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), diagy_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
               enddo
               !$OMP END PARALLEL DO
            endif

            ! FGMRES linear solver.  Everything it touches is compressed, so
            ! the solver never leaves the slot layout at all.
            call fgmres (Cb_s          , vrel_s    , &
                         umassdti_s    ,             &
                         halo_info_mask,             &
                         bx_s          , by_s      , &
                         diagx_s       , diagy_s   , &
                         reltol_fgmres , dim_fgmres, &
                         maxits_fgmres ,             &
                         solx_s        , soly_s    , &
                         nbiter)

            ! Put FGMRES solution in fpfunc vector (needed for Anderson).
            ! It is already compressed, so there is nothing to unpack.
            call arrays_to_vec (ntot, solx_s, soly_s, fpfunc(:))
         elseif (fpfunc_andacc == 2) then
            ! g_2(x) = x - A(x)x + b(x) = x - F(x)
            call abort_ice(error_message=subname // " Fixed point function g_2(x) not yet implemented (fpfunc_andacc = 2)" , &
               file=__FILE__, line=__LINE__)
         endif

         ! Compute fixed point residual f(x) = g(x) - x
         res = fpfunc - sol
#ifdef USE_LAPACK
         fpres_norm = global_sum(dnrm2(size(res), res, inc)**2, distrb_info)
#else
         ! The fixed-point residual is only wanted for its norm, so it goes
         ! straight into slots -- no block round trip at all.
         call vec_to_slots (ntot, res(:), fpx_s, fpy_s)
         fpres_norm = global_norm(fpx_s, fpy_s)
#endif
         if (my_task == master_task .and. monitor_nonlin) then
            write(nu_diag, '(a,i4,a,d26.16)') "monitor_nonlin: iter_nonlin= ", it_nl, &
                                              " fixed_point_res_L2norm= ", fpres_norm
         endif

         ! Not used for now (only nonlinear residual is checked)
         ! ! Store initial residual norm
         ! if (it_nl == 0) then
         !    tol = reltol_andacc*fpres_norm
         ! endif
         !
         ! ! Check residual
         ! if (fpres_norm < tol) then
         !    exit
         ! endif

         if (dim_andacc == 0 .or. it_nl < start_andacc) then
            ! Simple fixed point (Picard) iteration in this case
            sol = fpfunc
         else
#ifdef USE_LAPACK
            ! Begin Anderson acceleration
            if (get_num_procs() > 1) then
               ! Anderson solver is not yet parallelized; abort
               if (my_task == master_task) then
                  call abort_ice(error_message=subname // " Anderson solver (algo_nonlin = 'anderson') is not yet parallelized, and nprocs > 1 " , &
                     file=__FILE__, line=__LINE__)
               endif
            endif
            if (it_nl > start_andacc) then
               ! Update residual difference vector
               res_diff = res - res_old
               ! Update fixed point function difference matrix
               if (res_num < dim_andacc) then
                  ! Add column
                  G_diff(:,res_num+1) = fpfunc - fpfunc_old
               else
                  ! Delete first column and add column
                  G_diff(:,1:res_num-1) = G_diff(:,2:res_num)
                  G_diff(:,res_num) = fpfunc - fpfunc_old
               endif
               res_num = res_num + 1
            endif
            res_old = res
            fpfunc_old = fpfunc
            if (res_num == 0) then
               sol = fpfunc
            else
               if (res_num == 1) then
                  ! Initialize QR factorization
                  R(1,1) = dnrm2(size(res_diff), res_diff, inc)
                  Q(:,1) = res_diff/R(1,1)
               else
                  if (res_num > dim_andacc) then
                     ! Update factorization since 1st column was deleted
                     call qr_delete(Q,R)
                     res_num = res_num - 1
                  endif
                  ! Update QR factorization for new column
                  do j = 1, res_num - 1
                     R(j,res_num) = ddot(ntot, Q(:,j), inc, res_diff, inc)
                     res_diff = res_diff - R(j,res_num) * Q(:,j)
                  enddo
                  R(res_num, res_num) = dnrm2(size(res_diff) ,res_diff, inc)
                  Q(:,res_num) = res_diff / R(res_num, res_num)
               endif
               ! TODO: here, drop more columns to improve conditioning
               ! if (droptol) then

               ! endif
               ! Solve least square problem for coefficients
               ! 1. Compute rhs_tri = Q^T * res
               call dgemv ('t', size(Q,1), res_num, c1, Q(:,1:res_num), size(Q,1), res, inc, c0, rhs_tri, inc)
               ! 2. Solve R*coeffs = rhs_tri, put result in rhs_tri
               call dtrsv ('u', 'n', 'n', res_num, R(1:res_num,1:res_num), res_num, rhs_tri, inc)
               coeffs = rhs_tri
               ! Update approximate solution: x = fpfunc - G_diff*coeffs, put result in fpfunc
               call dgemv ('n', size(G_diff,1), res_num, -c1, G_diff(:,1:res_num), size(G_diff,1), coeffs, inc, c1, fpfunc, inc)
               sol = fpfunc
               ! Apply damping
               if (damping_andacc > 0 .and. damping_andacc /= 1) then
                  ! x = x - (1-beta) (res - Q*R*coeffs)

                  ! tmp = R*coeffs
                  call dgemv ('n', res_num, res_num, c1, R(1:res_num,1:res_num), res_num, coeffs, inc, c0, tmp, inc)
                  ! res = res - Q*tmp
                  call dgemv ('n', size(Q,1), res_num, -c1, Q(:,1:res_num), size(Q,1), tmp, inc, c1, res, inc)
                  ! x = x - (1-beta)*res
                  sol = sol - (1-damping_andacc)*res
               endif
            endif
#else
            ! Anderson solver is not usable without LAPACK; abort
            call abort_ice(error_message=subname // " CICE was not compiled with LAPACK, "// &
                           "and Anderson solver was chosen (algo_nonlin = 'anderson')" ,     &
               file=__FILE__, line=__LINE__)
#endif
         endif

         !-----------------------------------------------------------------------
         !     Put vector sol in uvel and vvel arrays
         !-----------------------------------------------------------------------
         call vec_to_arrays (nx_block    , ny_block    , &
                             nblocks     , max_blocks  , &
                             icellU   (:), ntot        , &
                             indxUi (:,:), indxUj (:,:), &
                             sol      (:),               &
                             uvel(:,:), vvel(:,:))

         ! Do halo update so that halo cells contain up to date info for advection
         call stack_fields(uvel, vvel, fld2)
         call ice_timer_start(timer_bound)
         if (maskhalo_dyn) then
            call ice_HaloUpdate (fld2,               halo_info_mask, &
                                 field_loc_NEcorner, field_type_vector)
         else
            call ice_HaloUpdate (fld2,               halo_info, &
                                 field_loc_NEcorner, field_type_vector)
         endif
         call ice_timer_stop(timer_bound)
         call unstack_fields(fld2, uvel, vvel)

         call ice_restoring_halo(setfld='velocity')

         ! Compute "progress" residual norm
         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks
            fpresx(:,iblk) = uvel(:,iblk) - uprev_k(:,iblk)
            fpresy(:,iblk) = vvel(:,iblk) - vprev_k(:,iblk)
            call pack_to_slots(fpresx(:,iblk), dynMap(:,iblk), fpx_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(fpresy(:,iblk), dynMap(:,iblk), fpy_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         enddo
         !$OMP END PARALLEL DO
         prog_norm = global_norm(fpx_s, fpy_s)
         if (my_task == master_task .and. monitor_nonlin) then
            write(nu_diag, '(a,i4,a,d26.16)') "monitor_nonlin: iter_nonlin= ", it_nl, &
                                              " progress_res_L2norm= ", prog_norm
         endif

      enddo ! nonlinear iteration loop

      end subroutine anderson_solver

!=======================================================================

! Computes the viscosities and dPr/dx, dPr/dy

      subroutine calc_zeta_dPr (ncells  , skipme  , &
                                uvel    , vvel    , &
                                dxT     , dyT     , &
                                dxhy    , dyhx    , &
                                cxp     , cyp     , &
                                cxm     , cym     , &
                                DminTarea,strength, &
                                nbr     ,           &
                                zetax2  , etax2   , &
                                rep_prs , stPr)

      use ice_dyn_shared, only: strain_rates, visc_replpress

      integer (kind=int_kind), intent(in) :: &
         ncells                ! slots held by this block

      logical (kind=log_kind), dimension (:), intent(in), contiguous :: &
         skipme       ! inactive-slot flag, T cells

      ! compressed: this block's slice, 1-based, cells 1..ncells, dead slot last
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         strength , & ! ice strength (N/m)
         uvel     , & ! x-component of velocity (m/s)
         vvel     , & ! y-component of velocity (m/s)
         dxT      , & ! width of T-cell through the middle (m)
         dyT      , & ! height of T-cell through the middle (m)
         dxhy     , & ! 0.5*(HTE - HTW)
         dyhx     , & ! 0.5*(HTN - HTS)
         cyp      , & ! 1.5*HTE - 0.5*HTW
         cxp      , & ! 1.5*HTN - 0.5*HTS
         cym      , & ! 0.5*HTE - 1.5*HTW
         cxm      , & ! 0.5*HTN - 1.5*HTS
         DminTarea    ! deltaminVP*tarea

      ! (slot, direction) -> neighbour slot.  A missing neighbour points at the
      ! block's dead slot, which holds zero, so the four corners below need no
      ! branch -- the same arrangement EVP's stress uses.
      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nbr

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         zetax2   , & ! zetax2 = 2*zeta (bulk viscosity)
         etax2    , & ! etax2  = 2*eta  (shear viscosity)
         rep_prs      ! replacement pressure

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         stPr          ! stress combinations from replacement pressure

      ! local variables

      integer (kind=int_kind) :: &
         is             , & ! slot counter
         isW, isS, isSW     ! the three neighbour slots this cell reads

      real (kind=dbl_kind) :: &
        divune, divunw, divuse, divusw                , & ! divergence
        tensionne, tensionnw, tensionse, tensionsw    , & ! tension
        shearne, shearnw, shearse, shearsw            , & ! shearing
        Deltane, Deltanw, Deltase, Deltasw            , & ! Delt
        ssigpn, ssigps, ssigpe, ssigpw, ssigp1, ssigp2, &
        csigpne, csigpnw, csigpsw, csigpse            , &
        stressp_1, stressp_2, stressp_3, stressp_4    , &
        strp_tmp

      character(len=*), parameter :: subname = '(calc_zeta_dPr)'

      ! Initialize stPr, zetax2 and etax2 to zero
      ! (for cells where iceTmask is false)
      stPr   = c0
      zetax2 = c0
      etax2  = c0

      do is = 1, ncells
         if (skipme(is)) cycle

         isW  = nbr(is,nbr_W )
         isS  = nbr(is,nbr_S )
         isSW = nbr(is,nbr_SW)

      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------
         call strain_rates (                            uvel(is), uvel(isW), &
                            uvel(isS), uvel(isSW), &
                            vvel(is), vvel(isW), &
                            vvel(isS), vvel(isSW), &
                            dxT(is) , dyT(is) , &
                            cxp(is) , cyp(is) , &
                            cxm(is) , cym(is) , &
                            divune   , divunw   , &
                            divuse   , divusw   , &
                            tensionne, tensionnw, &
                            tensionse, tensionsw, &
                            shearne  , shearnw  , &
                            shearse  , shearsw  , &
                            Deltane  , Deltanw  , &
                            Deltase  , Deltasw)

      !-----------------------------------------------------------------
      ! viscosities and replacement pressure
      !-----------------------------------------------------------------

         call visc_replpress (strength(is)  , DminTarea(is)  , &
                              Deltane        , zetax2   (is,1), &
                              etax2   (is,1), rep_prs  (is,1))
         call visc_replpress (strength(is)  , DminTarea(is)  , &
                              Deltanw        , zetax2   (is,2), &
                              etax2   (is,2), rep_prs  (is,2))
         call visc_replpress (strength(is)  , DminTarea(is)  , &
                              Deltasw        , zetax2   (is,3), &
                              etax2   (is,3), rep_prs  (is,3))
         call visc_replpress (strength(is)  , DminTarea(is)  , &
                              Deltase        , zetax2   (is,4), &
                              etax2   (is,4), rep_prs  (is,4))

      !-----------------------------------------------------------------
      ! the stresses                            ! kg/s^2
      ! (1) northeast, (2) northwest, (3) southwest, (4) southeast
      !-----------------------------------------------------------------

         stressp_1 = -rep_prs(is,1)
         stressp_2 = -rep_prs(is,2)
         stressp_3 = -rep_prs(is,3)
         stressp_4 = -rep_prs(is,4)

      !-----------------------------------------------------------------
      ! combinations of the Pr related stresses for the momentum equation ! kg/s^2
      !-----------------------------------------------------------------

         ssigpn  = stressp_1 + stressp_2
         ssigps  = stressp_3 + stressp_4
         ssigpe  = stressp_1 + stressp_4
         ssigpw  = stressp_2 + stressp_3
         ssigp1  =(stressp_1 + stressp_3)*p055
         ssigp2  =(stressp_2 + stressp_4)*p055

         csigpne = p111*stressp_1 + ssigp2 + p027*stressp_3
         csigpnw = p111*stressp_2 + ssigp1 + p027*stressp_4
         csigpsw = p111*stressp_3 + ssigp2 + p027*stressp_1
         csigpse = p111*stressp_4 + ssigp1 + p027*stressp_2

      !-----------------------------------------------------------------
      ! for dF/dx (u momentum)
      !-----------------------------------------------------------------
         strp_tmp  = p25*dyT(is)*(p333*ssigpn  + p166*ssigps)

         ! northeast(is)
         stPr(is,1) = -strp_tmp &
              + dxhy(is)*(-csigpne)

         ! northwest((j-1)*nx_block+i+1)
         stPr(is,2) = strp_tmp  &
              + dxhy(is)*(-csigpnw)

         strp_tmp  = p25*dyT(is)*(p333*ssigps  + p166*ssigpn)

         ! southeast((j+1-1)*nx_block+i)
         stPr(is,3) = -strp_tmp &
              + dxhy(is)*(-csigpse)

         ! southwest((j+1-1)*nx_block+i+1)
         stPr(is,4) = strp_tmp  &
              + dxhy(is)*(-csigpsw)

      !-----------------------------------------------------------------
      ! for dF/dy (v momentum)
      !-----------------------------------------------------------------
         strp_tmp  = p25*dxT(is)*(p333*ssigpe  + p166*ssigpw)

         ! northeast(is)
         stPr(is,5) = -strp_tmp &
              - dyhx(is)*(csigpne)

         ! southeast((j+1-1)*nx_block+i)
         stPr(is,6) = strp_tmp  &
              - dyhx(is)*(csigpse)

         strp_tmp  = p25*dxT(is)*(p333*ssigpw  + p166*ssigpe)

         ! northwest((j-1)*nx_block+i+1)
         stPr(is,7) = -strp_tmp &
              - dyhx(is)*(csigpnw)

         ! southwest((j+1-1)*nx_block+i+1)
         stPr(is,8) = strp_tmp  &
              - dyhx(is)*(csigpsw)

      enddo                     ! is

      end subroutine calc_zeta_dPr

!=======================================================================

! Computes the VP stresses (as diagnostic)

! Lemieux, J.-F., and Dupont, F. (2020), On the calculation of normalized
! viscous-plastic sea ice stresses, Geosci. Model Dev., 13, 1763–1769,

      subroutine stress_vp (ncells    , skipme    , &
                            uvel      , vvel      , &
                            dxT       , dyT       , &
                            cxp       , cyp       , &
                            cxm       , cym       , &
                            nbr       , ijslot    , &
                            zetax2    , etax2     , &
                            rep_prs   ,             &
                            stressp_1 , stressp_2 , &
                            stressp_3 , stressp_4 , &
                            stressm_1 , stressm_2 , &
                            stressm_3 , stressm_4 , &
                            stress12_1, stress12_2, &
                            stress12_3, stress12_4)

      use ice_dyn_shared, only: strain_rates

      integer (kind=int_kind), intent(in) :: &
         ncells                ! slots held by this block

      logical (kind=log_kind), dimension (:), intent(in), contiguous :: &
         skipme       ! inactive-slot flag, T cells

      ! compressed: this block's slice, 1-based, cells 1..ncells, dead slot last
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         uvel     , & ! x-component of velocity (m/s)
         vvel     , & ! y-component of velocity (m/s)
         dxT      , & ! width of T-cell through the middle (m)
         dyT      , & ! height of T-cell through the middle (m)
         cyp      , & ! 1.5*HTE - 0.5*HTW
         cxp      , & ! 1.5*HTN - 0.5*HTS
         cym      , & ! 0.5*HTE - 1.5*HTW
         cxm          ! 0.5*HTN - 1.5*HTS

      ! (slot, direction) -> neighbour slot.  A missing neighbour points at the
      ! block's dead slot, which holds zero, so the corners below need no
      ! branch -- the same arrangement calc_zeta_dPr and EVP's stress use.
      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nbr

      ! (j-1)*nx_block+i of the cell owning each slot.  The twelve stresses are
      ! restart and history fields, so they stay in the block layout; writing
      ! them through ijslot touches exactly the cells the icellT loop touched
      ! and leaves every other cell at its previous value, as before.  That is
      ! why there is no unpack here: an unpack would also write the cells the
      ! block-layout loop deliberately skipped.
      integer (kind=int_kind), dimension (:), intent(in) :: &
         ijslot

      real (kind=dbl_kind), dimension(:,:), intent(in) :: &
         zetax2   , & ! zetax2 = 2*zeta (bulk viscosity)
         etax2    , & ! etax2  = 2*eta  (shear viscosity)
         rep_prs

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         stressp_1, stressp_2, stressp_3, stressp_4 , & ! sigma11+sigma22
         stressm_1, stressm_2, stressm_3, stressm_4 , & ! sigma11-sigma22
         stress12_1,stress12_2,stress12_3,stress12_4    ! sigma12

      ! local variables

      integer (kind=int_kind) :: &
         is             , & ! slot counter
         isW, isS, isSW , & ! the three neighbour slots this cell reads
         ijb                ! 1D index of the owning cell in the block

      real (kind=dbl_kind) :: &
        divune, divunw, divuse, divusw            , & ! divergence
        tensionne, tensionnw, tensionse, tensionsw, & ! tension
        shearne, shearnw, shearse, shearsw        , & ! shearing
        Deltane, Deltanw, Deltase, Deltasw            ! Delt

      character(len=*), parameter :: subname = '(stress_vp)'

      do is = 1, ncells
         if (skipme(is)) cycle

         isW  = nbr(is,nbr_W )
         isS  = nbr(is,nbr_S )
         isSW = nbr(is,nbr_SW)
         ijb  = ijslot(is)

      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------
         call strain_rates (                            uvel(is), uvel(isW), &
                            uvel(isS), uvel(isSW),                           &
                            vvel(is), vvel(isW),                             &
                            vvel(isS), vvel(isSW),                           &
                            dxT(is) , dyT(is) ,                              &
                            cxp(is) , cyp(is) ,                              &
                            cxm(is) , cym(is) ,                              &
                            divune   , divunw   ,                            &
                            divuse   , divusw   ,                            &
                            tensionne, tensionnw,                            &
                            tensionse, tensionsw,                            &
                            shearne  , shearnw  ,                            &
                            shearse  , shearsw  ,                            &
                            Deltane  , Deltanw  ,                            &
                            Deltase  , Deltasw)

      !-----------------------------------------------------------------
      ! the stresses                            ! kg/s^2
      ! (1) northeast, (2) northwest, (3) southwest, (4) southeast
      !-----------------------------------------------------------------

         stressp_1(ijb) = zetax2(is,1)*divune - rep_prs(is,1)
         stressp_2(ijb) = zetax2(is,2)*divunw - rep_prs(is,2)
         stressp_3(ijb) = zetax2(is,3)*divusw - rep_prs(is,3)
         stressp_4(ijb) = zetax2(is,4)*divuse - rep_prs(is,4)

         stressm_1(ijb) = etax2(is,1)*tensionne
         stressm_2(ijb) = etax2(is,2)*tensionnw
         stressm_3(ijb) = etax2(is,3)*tensionsw
         stressm_4(ijb) = etax2(is,4)*tensionse

         stress12_1(ijb) = etax2(is,1)*shearne*p5
         stress12_2(ijb) = etax2(is,2)*shearnw*p5
         stress12_3(ijb) = etax2(is,3)*shearsw*p5
         stress12_4(ijb) = etax2(is,4)*shearse*p5

      enddo                     ! is

      end subroutine stress_vp

!=======================================================================

! Compute vrel and seabed stress coefficients

      subroutine calc_vrel_Cb (ncells  , skipme  , &
                               Cw      ,           &
                               aiU     , TbU     , &
                               uocn    , vocn    , &
                               uvel    , vvel    , &
                               vrel    , Cb)

      use ice_dyn_shared, only: u0 ! residual velocity for seabed stress (m/s)

      integer (kind=int_kind), intent(in) :: &
         ncells                ! slots held by this block

      ! skipme(is) is .true. where the slot carries no ice this step.  Together
      ! with ncells this replaces the icellU / indxUi / indxUj chain: the loop
      ! runs over slots, and the non-skipped ones are exactly the icellU cells
      ! in the same order, because dynMap numbers slots row-major over the
      ! block and dyn_prep2 scans the physical domain the same way.
      logical (kind=log_kind), dimension (:), intent(in), contiguous :: &
         skipme      ! inactive-slot flag

      ! every field below is COMPRESSED: this block's slice, 1-based, cells
      ! 1..ncells with the dead slot last
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         aiU     , & ! ice fraction on u-grid
         TbU     , & ! seabed stress factor (N/m^2)
         uocn    , & ! ocean current, x-direction (m/s)
         vocn    , & ! ocean current, y-direction (m/s)
         Cw      , & ! ocean-ice neutral drag coefficient
         uvel    , & ! x-component of velocity (m/s)
         vvel        ! y-component of velocity (m/s)

      real (kind=dbl_kind), dimension (:), intent(inout), contiguous :: &
         vrel      , & ! coeff for tauw
         Cb            ! seabed stress coeff

      ! local variables

      integer (kind=int_kind) :: &
         is                 ! slot counter

      real (kind=dbl_kind) :: &
         rhow                  !

      character(len=*), parameter :: subname = '(calc_vrel_Cb)'

      call icepack_query_parameters(rhow_out=rhow)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      do is = 1, ncells
         if (skipme(is)) cycle

         ! (magnitude of relative ocean current)*rhow*drag*aice
         vrel(is) = aiU(is)*rhow*Cw(is)*sqrt((uocn(is) - uvel(is))**2 + &
                                             (vocn(is) - vvel(is))**2)  ! m/s

         Cb(is)  = TbU(is) / (sqrt(uvel(is)**2 + vvel(is)**2) + u0) ! for seabed stress
      enddo                     ! is

      end subroutine calc_vrel_Cb

! Computes the matrix vector product A(u,v) * (u,v)
! Au = A(u,v)_[x] * uvel (x components of  A(u,v) * (u,v))
! Av = A(u,v)_[y] * vvel (y components of  A(u,v) * (u,v))

      subroutine matvec (ncells  , skipT   , &
                         skipU   , nbr     , &
                         dxT     , dyT     , &
                         dxhy    , dyhx    , &
                         cxp     , cyp     , &
                         cxm     , cym     , &
                         uvel    , vvel    , &
                         vrel    , Cb      , &
                         zetax2  , etax2   , &
                         umassdti, fmU     , &
                         uarear  ,           &
                         Au      , Av)

      use ice_dyn_shared, only: strain_rates

      integer (kind=int_kind), intent(in) :: &
         ncells                ! slots held by this block

      ! two flags, because this kernel has a T-cell pass and a U-cell pass
      logical (kind=log_kind), dimension (:), intent(in), contiguous :: &
         skipT   , & ! inactive-slot flag, T cells
         skipU       ! inactive-slot flag, U cells

      ! (slot, direction) -> neighbour slot; a missing one is the dead slot
      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nbr

      ! compressed: this block's slice, 1-based, cells 1..ncells, dead slot last
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         dxT      , & ! width of T-cell through the middle (m)
         dyT      , & ! height of T-cell through the middle (m)
         dxhy     , & ! 0.5*(HTE - HTW)
         dyhx     , & ! 0.5*(HTN - HTS)
         cyp      , & ! 1.5*HTE - 0.5*HTW
         cxp      , & ! 1.5*HTN - 0.5*HTS
         cym      , & ! 0.5*HTE - 1.5*HTW
         cxm      , & ! 0.5*HTN - 1.5*HTS
         umassdti , & ! mass of U-cell/dt (kg/m^2 s)
         uarear   , & ! 1/uarea
         uvel     , & ! x-component of velocity (m/s)
         vvel     , & ! y-component of velocity (m/s)
         vrel     , & ! coefficient for tauw
         Cb       , & ! coefficient for seabed stress
         fmU          ! Coriolis param. * mass in U-cell (kg/s)

      real (kind=dbl_kind), dimension(:,:), intent(in) :: &
         zetax2   , & ! zetax2 = 2*zeta (bulk viscosity)
         etax2        ! etax2  = 2*eta  (shear viscosity)

      real (kind=dbl_kind), dimension (:), intent(inout), contiguous :: &
         Au      , & ! matvec, Fx = bx - Au (N/m^2)
         Av          ! matvec, Fy = by - Av (N/m^2)

      ! local variables

      integer (kind=int_kind) :: &
         is             , & ! slot counter
         isW, isS, isSW , & ! neighbour slots read by the T pass
         isE, isN, isNE     ! neighbour slots read by the U pass

      ! Eight vectors rather than one (slot,8) array, matching EVP's
      ! str1_slot..str8_slot.  The U pass below reads a DIFFERENT component at
      ! each neighbour -- str1 at this slot, str2 at E, str3 at N, str4 at NE
      ! -- so one array per component gives each its own stream, where a
      ! (slot,8) array would stride across both dimensions at once.  ncells+1:
      ! the last entry is the block's dead slot, which every missing neighbour
      ! points at and which stays zero.
      real (kind=dbl_kind), dimension(ncells+1):: &
         str1, str2, str3, str4, str5, str6, str7, str8

      real (kind=dbl_kind) :: &
         ccaimp,ccb        , & ! intermediate variables
         tmp_str2_E, tmp_str3_N, tmp_str4_NE, & ! str read once per neighbour
         tmp_str6_N, tmp_str7_E, tmp_str8_NE, &
         strintx, strinty      ! divergence of the internal stress tensor

      real (kind=dbl_kind) :: &
        divune, divunw, divuse, divusw            , & ! divergence
        tensionne, tensionnw, tensionse, tensionsw, & ! tension
        shearne, shearnw, shearse, shearsw        , & ! shearing
        Deltane, Deltanw, Deltase, Deltasw        , & ! Delt
        ssigpn, ssigps, ssigpe, ssigpw            , &
        ssigmn, ssigms, ssigme, ssigmw            , &
        ssig12n, ssig12s, ssig12e, ssig12w        , &
        ssigp1, ssigp2, ssigm1, ssigm2, ssig121, ssig122, &
        csigpne, csigpnw, csigpse, csigpsw        , &
        csigmne, csigmnw, csigmse, csigmsw        , &
        csig12ne, csig12nw, csig12se, csig12sw    , &
        str12ew, str12we, str12ns, str12sn        , &
        strp_tmp, strm_tmp

      real (kind=dbl_kind) :: &
         stressp_1, stressp_2, stressp_3, stressp_4 , & ! sigma11+sigma22 (without Pr)
         stressm_1, stressm_2, stressm_3, stressm_4 , & ! sigma11-sigma22
         stress12_1,stress12_2,stress12_3,stress12_4    ! sigma12

      character(len=*), parameter :: subname = '(matvec)'

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      str1(:) = c0 ; str2(:) = c0 ; str3(:) = c0 ; str4(:) = c0
      str5(:) = c0 ; str6(:) = c0 ; str7(:) = c0 ; str8(:) = c0

      do is = 1, ncells
         if (skipT(is)) cycle

         isW  = nbr(is,nbr_W )
         isS  = nbr(is,nbr_S )
         isSW = nbr(is,nbr_SW)

      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------
         call strain_rates (                            uvel(is), uvel(isW), &
                            uvel(isS), uvel(isSW), &
                            vvel(is), vvel(isW), &
                            vvel(isS), vvel(isSW), &
                            dxT(is) , dyT(is) , &
                            cxp(is) , cyp(is) , &
                            cxm(is) , cym(is) , &
                            divune   , divunw   , &
                            divuse   , divusw   , &
                            tensionne, tensionnw, &
                            tensionse, tensionsw, &
                            shearne  , shearnw  , &
                            shearse  , shearsw  , &
                            Deltane  , Deltanw  , &
                            Deltase  , Deltasw)

      !-----------------------------------------------------------------
      ! the stresses                            ! kg/s^2
      ! (1) northeast, (2) northwest, (3) southwest, (4) southeast
      ! NOTE: commented part of stressp is from the replacement pressure Pr
      !-----------------------------------------------------------------

         stressp_1 = zetax2(is,1)*divune! - Deltane*(c1-Ktens))
         stressp_2 = zetax2(is,2)*divunw! - Deltanw*(c1-Ktens))
         stressp_3 = zetax2(is,3)*divusw! - Deltasw*(c1-Ktens))
         stressp_4 = zetax2(is,4)*divuse! - Deltase*(c1-Ktens))

         stressm_1 = etax2(is,1)*tensionne
         stressm_2 = etax2(is,2)*tensionnw
         stressm_3 = etax2(is,3)*tensionsw
         stressm_4 = etax2(is,4)*tensionse

         stress12_1 = etax2(is,1)*shearne*p5
         stress12_2 = etax2(is,2)*shearnw*p5
         stress12_3 = etax2(is,3)*shearsw*p5
         stress12_4 = etax2(is,4)*shearse*p5

      !-----------------------------------------------------------------
      ! combinations of the stresses for the momentum equation ! kg/s^2
      !-----------------------------------------------------------------

         ssigpn  = stressp_1 + stressp_2
         ssigps  = stressp_3 + stressp_4
         ssigpe  = stressp_1 + stressp_4
         ssigpw  = stressp_2 + stressp_3
         ssigp1  =(stressp_1 + stressp_3)*p055
         ssigp2  =(stressp_2 + stressp_4)*p055

         ssigmn  = stressm_1 + stressm_2
         ssigms  = stressm_3 + stressm_4
         ssigme  = stressm_1 + stressm_4
         ssigmw  = stressm_2 + stressm_3
         ssigm1  =(stressm_1 + stressm_3)*p055
         ssigm2  =(stressm_2 + stressm_4)*p055

         ssig12n = stress12_1 + stress12_2
         ssig12s = stress12_3 + stress12_4
         ssig12e = stress12_1 + stress12_4
         ssig12w = stress12_2 + stress12_3
         ssig121 =(stress12_1 + stress12_3)*p111
         ssig122 =(stress12_2 + stress12_4)*p111

         csigpne = p111*stressp_1 + ssigp2 + p027*stressp_3
         csigpnw = p111*stressp_2 + ssigp1 + p027*stressp_4
         csigpsw = p111*stressp_3 + ssigp2 + p027*stressp_1
         csigpse = p111*stressp_4 + ssigp1 + p027*stressp_2

         csigmne = p111*stressm_1 + ssigm2 + p027*stressm_3
         csigmnw = p111*stressm_2 + ssigm1 + p027*stressm_4
         csigmsw = p111*stressm_3 + ssigm2 + p027*stressm_1
         csigmse = p111*stressm_4 + ssigm1 + p027*stressm_2

         csig12ne = p222*stress12_1 + ssig122 &
                  + p055*stress12_3
         csig12nw = p222*stress12_2 + ssig121 &
                  + p055*stress12_4
         csig12sw = p222*stress12_3 + ssig122 &
                  + p055*stress12_1
         csig12se = p222*stress12_4 + ssig121 &
                  + p055*stress12_2

         str12ew = p5*dxT(is)*(p333*ssig12e + p166*ssig12w)
         str12we = p5*dxT(is)*(p333*ssig12w + p166*ssig12e)
         str12ns = p5*dyT(is)*(p333*ssig12n + p166*ssig12s)
         str12sn = p5*dyT(is)*(p333*ssig12s + p166*ssig12n)

      !-----------------------------------------------------------------
      ! for dF/dx (u momentum)
      !-----------------------------------------------------------------
         strp_tmp  = p25*dyT(is)*(p333*ssigpn  + p166*ssigps)
         strm_tmp  = p25*dyT(is)*(p333*ssigmn  + p166*ssigms)

         ! northeast(is)
         str1(is) = -strp_tmp - strm_tmp - str12ew &
              + dxhy(is)*(-csigpne + csigmne) + dyhx(is)*csig12ne

         ! northwest, read at this cell's E neighbour
         str2(is) = strp_tmp + strm_tmp - str12we &
              + dxhy(is)*(-csigpnw + csigmnw) + dyhx(is)*csig12nw

         strp_tmp  = p25*dyT(is)*(p333*ssigps  + p166*ssigpn)
         strm_tmp  = p25*dyT(is)*(p333*ssigms  + p166*ssigmn)

         ! southeast, read at this cell's N neighbour
         str3(is) = -strp_tmp - strm_tmp + str12ew &
              + dxhy(is)*(-csigpse + csigmse) + dyhx(is)*csig12se

         ! southwest, read at this cell's NE neighbour
         str4(is) = strp_tmp + strm_tmp + str12we &
              + dxhy(is)*(-csigpsw + csigmsw) + dyhx(is)*csig12sw

      !-----------------------------------------------------------------
      ! for dF/dy (v momentum)
      !-----------------------------------------------------------------
         strp_tmp  = p25*dxT(is)*(p333*ssigpe  + p166*ssigpw)
         strm_tmp  = p25*dxT(is)*(p333*ssigme  + p166*ssigmw)

         ! northeast(is)
         str5(is) = -strp_tmp + strm_tmp - str12ns &
              - dyhx(is)*(csigpne + csigmne) + dxhy(is)*csig12ne

         ! southeast, read at this cell's N neighbour
         str6(is) = strp_tmp - strm_tmp - str12sn &
              - dyhx(is)*(csigpse + csigmse) + dxhy(is)*csig12se

         strp_tmp  = p25*dxT(is)*(p333*ssigpw  + p166*ssigpe)
         strm_tmp  = p25*dxT(is)*(p333*ssigmw  + p166*ssigme)

         ! northwest, read at this cell's E neighbour
         str7(is) = -strp_tmp + strm_tmp + str12ns &
              - dyhx(is)*(csigpnw + csigmnw) + dxhy(is)*csig12nw

         ! southwest, read at this cell's NE neighbour
         str8(is) = strp_tmp - strm_tmp + str12sn &
              - dyhx(is)*(csigpsw + csigmsw) + dxhy(is)*csig12sw

      enddo ! ij - icellT

      !-----------------------------------------------------------------
      ! Form Au and Av
      !-----------------------------------------------------------------

      do is = 1, ncells
         if (skipU(is)) cycle

         isE  = nbr(is,nbr_E )
         isN  = nbr(is,nbr_N )
         isNE = nbr(is,nbr_NE)

         ! each neighbour slot read once, the six dependent loads issued
         ! together, as ice_dyn_core1d's stepu_1d does
         tmp_str2_E  = str2(isE )
         tmp_str3_N  = str3(isN )
         tmp_str4_NE = str4(isNE)
         tmp_str6_N  = str6(isN )
         tmp_str7_E  = str7(isE )
         tmp_str8_NE = str8(isNE)

         ccaimp = umassdti(is) + vrel(is) * cosw + Cb(is) ! kg/m^2 s

         ccb = fmU(is) + sign(c1,fmU(is)) * vrel(is) * sinw ! kg/m^2 s

         ! divergence of the internal stress tensor
         strintx = uarear(is)* &
             (str1(is) + tmp_str2_E + tmp_str3_N + tmp_str4_NE)
         strinty = uarear(is)* &
             (str5(is) + tmp_str6_N + tmp_str7_E + tmp_str8_NE)

         Au(is) = ccaimp*uvel(is) - ccb*vvel(is) - strintx
         Av(is) = ccaimp*vvel(is) + ccb*uvel(is) - strinty
      enddo ! ij - icellU

      end subroutine matvec

!=======================================================================

! Compute the constant component of b(u,v) i.e. the part of b(u,v) that
! does not depend on (u,v) and thus do not change during the nonlinear iteration

     subroutine calc_bfix  (nx_block , ny_block , &
                            icellU   ,            &
                            indxUi   , indxUj   , &
                            umassdti ,            &
                            forcexU  , forceyU  , &
                            uvel_init, vvel_init, &
                            Xmap,                 &
                            bxfix    , byfix)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icellU                ! no. of cells where iceUmask = .true.

      integer (kind=int_kind), dimension (nx_block*ny_block), intent(in) :: &
         indxUi   , & ! compressed index in i-direction
         indxUj       ! compressed index in j-direction

      ! umassdti is assumed shape, 1D over the block in Fortran storage order,
      ! i.e. element((j-1)*nx_block+i) is stored at (j-1)*nx_block+i
      real (kind=dbl_kind), dimension (:), intent(in) :: &
         umassdti, & ! mass of U-cell/dt (kg/m^2 s)
         uvel_init,& ! x-component of velocity (m/s), start of step, COMPRESSED
         vvel_init,& ! y-component of velocity (m/s), start of step, COMPRESSED
         forcexU , & ! work array: combined atm stress and ocn tilt, x
         forceyU     ! work array: combined atm stress and ocn tilt, y

      ! Xmap holds each cell's static compressed slot, see ice_dyn_shared
      integer (kind=int_kind), dimension (:), intent(in) :: &
         Xmap

      real (kind=dbl_kind), dimension (nx_block*ny_block), intent(out) :: &
         bxfix   , & ! bx = taux + bxfix
         byfix       ! by = tauy + byfix

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij       , & ! ij is the compressed cell counter
         ijb              , & ! 1D index of((j-1)*nx_block+i) in the block
         is                   ! compressed slot of((j-1)*nx_block+i)

      character(len=*), parameter :: subname = '(calc_bfix)'

      do ij = 1, icellU
         i = indxUi(ij)
         j = indxUj(ij)
         ijb = (j-1)*nx_block + i
         is  = Xmap(ijb)

         bxfix((j-1)*nx_block+i) = umassdti(ijb)*uvel_init(is) + forcexU(ijb)
         byfix((j-1)*nx_block+i) = umassdti(ijb)*vvel_init(is) + forceyU(ijb)
      enddo

      end subroutine calc_bfix

!=======================================================================

! Compute the vector b(u,v), i.e. the part of the nonlinear function F(u,v)
! that cannot be written as A(u,v)*(u,v), where A(u,v) is a matrix with entries
! depending on (u,v)

      subroutine calc_bvec (ncells  , skipme  , &
                            stPr    , uarear  , &
                            waterxU , wateryU , &
                            bxfix   , byfix   , &
                            bx      , by      , &
                            vrel    , nbr)

      integer (kind=int_kind), intent(in) :: &
         ncells                ! slots held by this block

      logical (kind=log_kind), dimension (:), intent(in), contiguous :: &
         skipme      ! inactive-slot flag

      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         waterxU , & ! for ocean stress calculation, x (m/s)
         wateryU , & ! for ocean stress calculation, y (m/s)
         uarear  , & ! 1/uarea
         bxfix   , & ! bx = taux + bxfix
         byfix   , & ! by = tauy + byfix
         vrel        ! relative ice-ocean velocity

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         stPr

      ! (slot, direction) -> neighbour slot; a missing one is the dead slot
      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nbr

      real (kind=dbl_kind), dimension (:), intent(inout), contiguous :: &
         bx      , & ! b vector, bx = taux + bxfix (N/m^2)
         by          ! b vector, by = tauy + byfix (N/m^2)

      ! local variables

      integer (kind=int_kind) :: &
         is             , & ! slot counter
         isE, isN, isNE     ! the three neighbour slots this cell reads

      real (kind=dbl_kind) :: &
         taux, tauy        , & ! part of ocean stress term
         strintx, strinty  , & ! divergence of the internal stress tensor (only Pr contributions)
         rhow                  !

      character(len=*), parameter :: subname = '(calc_bvec)'

      !-----------------------------------------------------------------
      ! calc b vector
      !-----------------------------------------------------------------

      call icepack_query_parameters(rhow_out=rhow)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      do is = 1, ncells
         if (skipme(is)) cycle

         ! the three T cells whose Pr stresses reach this U point, as EVP's
         ! stepu reads str: (+1,0), (0,+1) and (+1,+1)
         isE  = nbr(is,nbr_E )
         isN  = nbr(is,nbr_N )
         isNE = nbr(is,nbr_NE)

         ! ice/ocean stress
         taux = vrel(is)*waterxU(is) ! NOTE this is not the entire
         tauy = vrel(is)*wateryU(is) ! ocn stress term

         ! divergence of the internal stress tensor (only Pr part, i.e. dPr/dx, dPr/dy)
         strintx = uarear(is)* &
             (stPr(is,1) + stPr(isE,2) + stPr(isN,3) + stPr(isNE,4))
         strinty = uarear(is)* &
             (stPr(is,5) + stPr(isN,6) + stPr(isE,7) + stPr(isNE,8))

         bx(is) = bxfix(is) + taux + strintx
         by(is) = byfix(is) + tauy + strinty
      enddo                     ! is

      end subroutine calc_bvec

!=======================================================================

! Compute the non linear residual F(u,v) = b(u,v) - A(u,v) * (u,v),
! with Au, Av precomputed as
! Au = A(u,v)_[x] * uvel (x components of  A(u,v) * (u,v))
! Av = A(u,v)_[y] * vvel (y components of  A(u,v) * (u,v))

      subroutine residual_vec (ncells     , skipme  , &
                               bx         , by      , &
                               Au         , Av      , &
                               Fx         , Fy        )

      integer (kind=int_kind), intent(in) :: &
         ncells                ! slots held by this block

      logical (kind=log_kind), dimension (:), intent(in), contiguous :: &
         skipme       ! inactive-slot flag, U cells

      ! compressed: this block's slice, 1-based, cells 1..ncells, dead slot last
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         bx       , & ! b vector, bx = taux + bxfix (N/m^2)
         by       , & ! b vector, by = tauy + byfix (N/m^2)
         Au       , & ! matvec, Fx = bx - Au (N/m^2)
         Av           ! matvec, Fy = by - Av (N/m^2)

      real (kind=dbl_kind), dimension (:), intent(inout), contiguous :: &
         Fx      , & ! x residual vector, Fx = bx - Au (N/m^2)
         Fy          ! y residual vector, Fy = by - Av (N/m^2)

      ! local variables

      integer (kind=int_kind) :: &
         is                 ! slot counter

      character(len=*), parameter :: subname = '(residual_vec)'

      !-----------------------------------------------------------------
      ! compute residual
      !-----------------------------------------------------------------

      do is = 1, ncells
         if (skipme(is)) cycle

         Fx(is) = bx(is) - Au(is)
         Fy(is) = by(is) - Av(is)
      enddo                     ! is

      end subroutine residual_vec

!=======================================================================

! Form the diagonal of the matrix A(u,v) (first part of the computation)
! Part 1: compute the contributions to the diagonal from the rheology term

      subroutine formDiag_step1  (ncells  , skipme  , &
                                  dxT     , dyT     , &
                                  dxhy    , dyhx    , &
                                  cxp     , cyp     , &
                                  cxm     , cym     , &
                                  nbr     ,           &
                                  zetax2  , etax2   , &
                                  Drheo)

      integer (kind=int_kind), intent(in) :: &
         ncells                ! slots held by this block

      logical (kind=log_kind), dimension (:), intent(in), contiguous :: &
         skipme       ! inactive-slot flag, U cells

      ! compressed: this block's slice, 1-based, cells 1..ncells, dead slot last
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         dxT      , & ! width of T-cell through the middle (m)
         dyT      , & ! height of T-cell through the middle (m)
         dxhy     , & ! 0.5*(HTE - HTW)
         dyhx     , & ! 0.5*(HTN - HTS)
         cyp      , & ! 1.5*HTE - 0.5*HTW
         cxp      , & ! 1.5*HTN - 0.5*HTS
         cym      , & ! 0.5*HTE - 1.5*HTW
         cxm          ! 0.5*HTN - 1.5*HTS

      ! (slot, direction) -> neighbour slot.  Each of the eight terms below
      ! reads one of the four T cells around this U cell: its own, or one step
      ! E, N or NE.  A missing neighbour points at the block's dead slot,
      ! which holds zero, so no edge test is needed.
      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nbr

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         zetax2   , & ! zetax2 = 2*zeta (bulk viscosity)
         etax2        ! etax2  = 2*eta  (shear viscosity)

      real (kind=dbl_kind), dimension (:,:), intent(out) :: &
         Drheo          ! intermediate value for diagonal components of matrix A associated
                       ! with rheology term

      ! local variables

      integer (kind=int_kind) :: &
         is             , & ! U slot counter
         ist            , & ! the T slot this term reads (self, E, N or NE)
         tdir           , & ! direction from the U slot to that T slot
         cc                 ! which of the eight diagonal terms

      real (kind=dbl_kind) :: &
        divune, divunw, divuse, divusw            , & ! divergence
        tensionne, tensionnw, tensionse, tensionsw, & ! tension
        shearne, shearnw, shearse, shearsw        , & ! shearing
        uij,ui1j,uij1,ui1j1,vij,vi1j,vij1,vi1j1   , & ! == c0 or c1
        stressp_1, stressp_2, stressp_3, stressp_4, &
        stressm_1, stressm_2, stressm_3, stressm_4, &
        stress12_1,stress12_2,stress12_3,stress12_4,&
        ssigpn, ssigps, ssigpe, ssigpw            , &
        ssigmn, ssigms, ssigme, ssigmw            , &
        ssig12n, ssig12s, ssig12e, ssig12w        , &
        ssigp1, ssigp2, ssigm1, ssigm2, ssig121, ssig122, &
        csigpne, csigpnw, csigpse, csigpsw        , &
        csigmne, csigmnw, csigmse, csigmsw        , &
        csig12ne, csig12nw, csig12se, csig12sw    , &
        str12ew, str12we, str12ns, str12sn        , &
        strp_tmp, strm_tmp

      character(len=*), parameter :: subname = '(formDiag_step1)'

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      Drheo(:,:) = c0

      ! Be careful: Drheo contains 4 terms for u and 4 terms for v.
      ! These 8 terms come from the surrounding T cells but are all
      ! refrerenced to the i,j (u point) :

      ! Drheo(ist,1) corresponds to str(ist,1)
      ! Drheo(ist,2) corresponds to str((j-1)*nx_block+i+1,2)
      ! Drheo(ist,3) corresponds to str((j+1-1)*nx_block+i,3)
      ! Drheo(ist,4) corresponds to str((j+1-1)*nx_block+i+1,4))
      ! Drheo(ist,5) corresponds to str(ist,5)
      ! Drheo(ist,6) corresponds to str((j+1-1)*nx_block+i,6)
      ! Drheo(ist,7) corresponds to str((j-1)*nx_block+i+1,7)
      ! Drheo(ist,8) corresponds to str((j+1-1)*nx_block+i+1,8))

      do cc = 1, 8 ! 4 for u and 4 for v

         if (cc == 1) then     ! u comp, T cell i,j
            uij   = c1
            ui1j  = c0
            uij1  = c0
            ui1j1 = c0
            vij   = c0
            vi1j  = c0
            vij1  = c0
            vi1j1 = c0
            tdir  = 0       ! this U cell's own T slot
         elseif (cc == 2) then ! u comp, T cell i+1,j
            uij   = c0
            ui1j  = c1
            uij1  = c0
            ui1j1 = c0
            vij   = c0
            vi1j  = c0
            vij1  = c0
            vi1j1 = c0
            tdir  = nbr_E   ! T cell one step east
         elseif (cc == 3) then ! u comp, T cell i,j+1
            uij   = c0
            ui1j  = c0
            uij1  = c1
            ui1j1 = c0
            vij   = c0
            vi1j  = c0
            vij1  = c0
            vi1j1 = c0
            tdir  = nbr_N   ! T cell one step north
         elseif (cc == 4) then ! u comp, T cell i+1,j+1
            uij   = c0
            ui1j  = c0
            uij1  = c0
            ui1j1 = c1
            vij   = c0
            vi1j  = c0
            vij1  = c0
            vi1j1 = c0
            tdir  = nbr_NE  ! T cell one step north-east
         elseif (cc == 5) then ! v comp, T cell i,j
            uij   = c0
            ui1j  = c0
            uij1  = c0
            ui1j1 = c0
            vij   = c1
            vi1j  = c0
            vij1  = c0
            vi1j1 = c0
            tdir  = 0       ! this U cell's own T slot
         elseif (cc == 6) then ! v comp, T cell i,j+1
            uij   = c0
            ui1j  = c0
            uij1  = c0
            ui1j1 = c0
            vij   = c0
            vi1j  = c0
            vij1  = c1
            vi1j1 = c0
            tdir  = nbr_N   ! T cell one step north
         elseif (cc == 7) then ! v comp, T cell i+1,j
            uij   = c0
            ui1j  = c0
            uij1  = c0
            ui1j1 = c0
            vij   = c0
            vi1j  = c1
            vij1  = c0
            vi1j1 = c0
            tdir  = nbr_E   ! T cell one step east
         elseif (cc == 8) then ! v comp, T cell i+1,j+1
            uij   = c0
            ui1j  = c0
            uij1  = c0
            ui1j1 = c0
            vij   = c0
            vi1j  = c0
            vij1  = c0
            vi1j1 = c1
            tdir  = nbr_NE  ! T cell one step north-east
         endif

         do is = 1, ncells
            if (skipme(is)) cycle

            ! ist is the SHIFTED T slot this cc term reads, not the U slot;
            ! every metric and viscosity below is read there.  A missing
            ! neighbour points at the block's dead slot, which holds zero, so
            ! this needs no edge test -- the same arrangement matvec uses.
            if (tdir == 0) then
               ist = is
            else
               ist = nbr(is,tdir)
            endif

         !-----------------------------------------------------------------
         ! strain rates
         ! NOTE these are actually strain rates * area  (m^2/s)
         !-----------------------------------------------------------------
            ! divergence  =  e_11 + e_22
            divune    = cyp(ist)*uij   - dyT(ist)*ui1j  &
                      + cxp(ist)*vij   - dxT(ist)*vij1
            divunw    = cym(ist)*ui1j  + dyT(ist)*uij   &
                      + cxp(ist)*vi1j  - dxT(ist)*vi1j1
            divusw    = cym(ist)*ui1j1 + dyT(ist)*uij1  &
                      + cxm(ist)*vi1j1 + dxT(ist)*vi1j
            divuse    = cyp(ist)*uij1  - dyT(ist)*ui1j1 &
                      + cxm(ist)*vij1  + dxT(ist)*vij

            ! tension strain rate  =  e_11 - e_22
            tensionne = -cym(ist)*uij   - dyT(ist)*ui1j  &
                      +  cxm(ist)*vij   + dxT(ist)*vij1
            tensionnw = -cyp(ist)*ui1j  + dyT(ist)*uij   &
                      +  cxm(ist)*vi1j  + dxT(ist)*vi1j1
            tensionsw = -cyp(ist)*ui1j1 + dyT(ist)*uij1  &
                      +  cxp(ist)*vi1j1 - dxT(ist)*vi1j
            tensionse = -cym(ist)*uij1  - dyT(ist)*ui1j1 &
                      +  cxp(ist)*vij1  - dxT(ist)*vij

            ! shearing strain rate  =  2*e_12
            shearne = -cym(ist)*vij   - dyT(ist)*vi1j  &
                    -  cxm(ist)*uij   - dxT(ist)*uij1
            shearnw = -cyp(ist)*vi1j  + dyT(ist)*vij   &
                    -  cxm(ist)*ui1j  - dxT(ist)*ui1j1
            shearsw = -cyp(ist)*vi1j1 + dyT(ist)*vij1  &
                    -  cxp(ist)*ui1j1 + dxT(ist)*ui1j
            shearse = -cym(ist)*vij1  - dyT(ist)*vi1j1 &
                    -  cxp(ist)*uij1  + dxT(ist)*uij

         !-----------------------------------------------------------------
         ! the stresses                            ! kg/s^2
         ! (1) northeast, (2) northwest, (3) southwest, (4) southeast
         !-----------------------------------------------------------------

            stressp_1 = zetax2(ist,1)*divune
            stressp_2 = zetax2(ist,2)*divunw
            stressp_3 = zetax2(ist,3)*divusw
            stressp_4 = zetax2(ist,4)*divuse

            stressm_1 = etax2(ist,1)*tensionne
            stressm_2 = etax2(ist,2)*tensionnw
            stressm_3 = etax2(ist,3)*tensionsw
            stressm_4 = etax2(ist,4)*tensionse

            stress12_1 = etax2(ist,1)*shearne*p5
            stress12_2 = etax2(ist,2)*shearnw*p5
            stress12_3 = etax2(ist,3)*shearsw*p5
            stress12_4 = etax2(ist,4)*shearse*p5

         !-----------------------------------------------------------------
         ! combinations of the stresses for the momentum equation ! kg/s^2
         !-----------------------------------------------------------------

            ssigpn  = stressp_1 + stressp_2
            ssigps  = stressp_3 + stressp_4
            ssigpe  = stressp_1 + stressp_4
            ssigpw  = stressp_2 + stressp_3
            ssigp1  =(stressp_1 + stressp_3)*p055
            ssigp2  =(stressp_2 + stressp_4)*p055

            ssigmn  = stressm_1 + stressm_2
            ssigms  = stressm_3 + stressm_4
            ssigme  = stressm_1 + stressm_4
            ssigmw  = stressm_2 + stressm_3
            ssigm1  =(stressm_1 + stressm_3)*p055
            ssigm2  =(stressm_2 + stressm_4)*p055

            ssig12n = stress12_1 + stress12_2
            ssig12s = stress12_3 + stress12_4
            ssig12e = stress12_1 + stress12_4
            ssig12w = stress12_2 + stress12_3
            ssig121 =(stress12_1 + stress12_3)*p111
            ssig122 =(stress12_2 + stress12_4)*p111

            csigpne = p111*stressp_1 + ssigp2 + p027*stressp_3
            csigpnw = p111*stressp_2 + ssigp1 + p027*stressp_4
            csigpsw = p111*stressp_3 + ssigp2 + p027*stressp_1
            csigpse = p111*stressp_4 + ssigp1 + p027*stressp_2

            csigmne = p111*stressm_1 + ssigm2 + p027*stressm_3
            csigmnw = p111*stressm_2 + ssigm1 + p027*stressm_4
            csigmsw = p111*stressm_3 + ssigm2 + p027*stressm_1
            csigmse = p111*stressm_4 + ssigm1 + p027*stressm_2

            csig12ne = p222*stress12_1 + ssig122 &
                     + p055*stress12_3
            csig12nw = p222*stress12_2 + ssig121 &
                     + p055*stress12_4
            csig12sw = p222*stress12_3 + ssig122 &
                     + p055*stress12_1
            csig12se = p222*stress12_4 + ssig121 &
                     + p055*stress12_2

            str12ew = p5*dxT(ist)*(p333*ssig12e + p166*ssig12w)
            str12we = p5*dxT(ist)*(p333*ssig12w + p166*ssig12e)
            str12ns = p5*dyT(ist)*(p333*ssig12n + p166*ssig12s)
            str12sn = p5*dyT(ist)*(p333*ssig12s + p166*ssig12n)

         !-----------------------------------------------------------------
         ! for dF/dx (u momentum)
         !-----------------------------------------------------------------

            if (cc == 1) then ! T cell i,j

               strp_tmp  = p25*dyT(ist)*(p333*ssigpn  + p166*ssigps)
               strm_tmp  = p25*dyT(ist)*(p333*ssigmn  + p166*ssigms)

               ! northeast(ist)
               Drheo(is,1) = -strp_tmp - strm_tmp - str12ew &
                  + dxhy(ist)*(-csigpne + csigmne) + dyhx(ist)*csig12ne

            elseif (cc == 2) then ! T cell i+1,j

               strp_tmp  = p25*dyT(ist)*(p333*ssigpn  + p166*ssigps)
               strm_tmp  = p25*dyT(ist)*(p333*ssigmn  + p166*ssigms)

               ! northwest((j-1)*nx_block+i+1)
               Drheo(is,2) = strp_tmp + strm_tmp - str12we &
                  + dxhy(ist)*(-csigpnw + csigmnw) + dyhx(ist)*csig12nw

            elseif (cc == 3) then ! T cell i,j+1

               strp_tmp  = p25*dyT(ist)*(p333*ssigps  + p166*ssigpn)
               strm_tmp  = p25*dyT(ist)*(p333*ssigms  + p166*ssigmn)

               ! southeast((j+1-1)*nx_block+i)
               Drheo(is,3) = -strp_tmp - strm_tmp + str12ew &
                  + dxhy(ist)*(-csigpse + csigmse) + dyhx(ist)*csig12se

            elseif (cc == 4) then ! T cell i+1,j+1

               strp_tmp  = p25*dyT(ist)*(p333*ssigps  + p166*ssigpn)
               strm_tmp  = p25*dyT(ist)*(p333*ssigms  + p166*ssigmn)

               ! southwest((j+1-1)*nx_block+i+1)
               Drheo(is,4) = strp_tmp + strm_tmp + str12we &
                  + dxhy(ist)*(-csigpsw + csigmsw) + dyhx(ist)*csig12sw

         !-----------------------------------------------------------------
         ! for dF/dy (v momentum)
         !-----------------------------------------------------------------

            elseif (cc == 5) then ! T cell i,j

               strp_tmp  = p25*dxT(ist)*(p333*ssigpe  + p166*ssigpw)
               strm_tmp  = p25*dxT(ist)*(p333*ssigme  + p166*ssigmw)

               ! northeast(ist)
               Drheo(is,5) = -strp_tmp + strm_tmp - str12ns &
                  - dyhx(ist)*(csigpne + csigmne) + dxhy(ist)*csig12ne

            elseif (cc == 6) then ! T cell i,j+1

               strp_tmp  = p25*dxT(ist)*(p333*ssigpe  + p166*ssigpw)
               strm_tmp  = p25*dxT(ist)*(p333*ssigme  + p166*ssigmw)

               ! southeast((j+1-1)*nx_block+i)
               Drheo(is,6) = strp_tmp - strm_tmp - str12sn &
                  - dyhx(ist)*(csigpse + csigmse) + dxhy(ist)*csig12se

            elseif (cc == 7) then ! T cell i,j+1

               strp_tmp  = p25*dxT(ist)*(p333*ssigpw  + p166*ssigpe)
               strm_tmp  = p25*dxT(ist)*(p333*ssigmw  + p166*ssigme)

               ! northwest((j-1)*nx_block+i+1)
               Drheo(is,7) = -strp_tmp + strm_tmp + str12ns &
                  - dyhx(ist)*(csigpnw + csigmnw) + dxhy(ist)*csig12nw

            elseif (cc == 8) then ! T cell i+1,j+1

               strp_tmp  = p25*dxT(ist)*(p333*ssigpw  + p166*ssigpe)
               strm_tmp  = p25*dxT(ist)*(p333*ssigmw  + p166*ssigme)

               ! southwest((j+1-1)*nx_block+i+1)
               Drheo(is,8) = strp_tmp - strm_tmp + str12sn &
                  - dyhx(ist)*(csigpsw + csigmsw) + dxhy(ist)*csig12sw

            endif

         enddo                     ! ij

      enddo                     ! cc

      end subroutine formDiag_step1

!=======================================================================

! Form the diagonal of the matrix A(u,v) (second part of the computation)
! Part 2: compute diagonal

      subroutine formDiag_step2 (ncells  , skipme  , &
                                 Drheo   , vrel    , &
                                 umassdti,           &
                                 uarear  , Cb      , &
                                 diagx   , diagy)

      integer (kind=int_kind), intent(in) :: &
         ncells                ! slots held by this block

      logical (kind=log_kind), dimension (:), intent(in), contiguous :: &
         skipme       ! inactive-slot flag, U cells

      ! compressed: this block's slice, 1-based, cells 1..ncells, dead slot last
      real (kind=dbl_kind), dimension (:), intent(in) :: &
         umassdti, & ! mass of U-cell/dt (kg/m^2 s)
         uarear  , & ! 1/uarea
         vrel    , & ! coefficient for tauw
         Cb          ! coefficient for seabed stress

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         Drheo

      ! compressed as well: fgmres reads the diagonal in the slot layout, so
      ! there is no reason to build it block-shaped and pack it afterwards.
      real (kind=dbl_kind), dimension (:), intent(out) :: &
         diagx   , & ! Diagonal (x component) of the matrix A
         diagy       ! Diagonal (y component) of the matrix A

      ! local variables

      integer (kind=int_kind) :: &
         is                 ! slot counter

      real (kind=dbl_kind) :: &
         ccaimp             , & ! intermediate variables
         strintx, strinty       ! diagonal contributions to the divergence

      character(len=*), parameter :: subname = '(formDiag_step2)'

      !-----------------------------------------------------------------
      ! integrate the momentum equation
      !-----------------------------------------------------------------

      strintx = c0
      strinty = c0

      ! Be careful: Drheo contains 4 terms for u and 4 terms for v.
      ! These 8 terms come from the surrounding T cells but are all
      ! refrerenced to the i,j (u point) :

      ! Drheo(ijb,1) corresponds to str(ijb,1)
      ! Drheo(ijb,2) corresponds to str((j-1)*nx_block+i+1,2)
      ! Drheo(ijb,3) corresponds to str((j+1-1)*nx_block+i,3)
      ! Drheo(ijb,4) corresponds to str((j+1-1)*nx_block+i+1,4))
      ! Drheo(ijb,5) corresponds to str(ijb,5)
      ! Drheo(ijb,6) corresponds to str((j+1-1)*nx_block+i,6)
      ! Drheo(ijb,7) corresponds to str((j-1)*nx_block+i+1,7)
      ! Drheo(ijb,8) corresponds to str((j+1-1)*nx_block+i+1,8))

      do is = 1, ncells
         if (skipme(is)) cycle

         ccaimp = umassdti(is) + vrel(is) * cosw + Cb(is) ! kg/m^2 s

         strintx = uarear(is)* &
             (Drheo(is,1) + Drheo(is,2) + Drheo(is,3) + Drheo(is,4))
         strinty = uarear(is)* &
             (Drheo(is,5) + Drheo(is,6) + Drheo(is,7) + Drheo(is,8))

         diagx(is) = ccaimp - strintx
         diagy(is) = ccaimp - strinty
      enddo                     ! is

      end subroutine formDiag_step2

!=======================================================================

! Compute global l^2 norm of a vector field (field_x, field_y)

      function global_norm (field_x , field_y ) &
               result(norm)

      ! compressed: the whole local slot vector, every block
      real (kind=dbl_kind), dimension (dynTot), intent(in) :: &
         field_x     , & ! x-component of vector field
         field_y         ! y-component of vector field

      real (kind=dbl_kind) :: &
         norm      ! l^2 norm of vector field

      ! local variables

      character(len=*), parameter :: subname = '(global_norm)'

      norm = sqrt(global_dot_product (field_x   , field_y  , &
                                      field_x   , field_y  ))

      end function global_norm

!=======================================================================

! Compute global dot product of two grid vectors, each split into X and Y components

      function global_dot_product (vector1_x , vector1_y, &
                                   vector2_x , vector2_y) &
               result(dot_product)

      use ice_blocks, only: ns_boundary_type, nx_block, ny_block
      use ice_domain, only: distrb_info
      use ice_domain_size, only: max_blocks
      use ice_fileunits, only: bfbflag

      ! compressed: the whole local slot vector, every block
      real (kind=dbl_kind), dimension (dynTot), intent(in) :: &
         vector1_x     , & ! x-component of first vector
         vector1_y     , & ! y-component of first vector
         vector2_x     , & ! x-component of second vector
         vector2_y         ! y-component of second vector

      real (kind=dbl_kind) :: &
         dot_product      ! l^2 norm of vector field

      ! local variables

      integer (kind=int_kind) :: &
         is, iblk

      real (kind=dbl_kind) :: &
         p         ! this cell's contribution

      ! prod stays in the BLOCK layout deliberately.  The reprosum branch
      ! below hands the whole 2-D field to global_sum, which masks it by field
      ! location, so it cannot take a compressed vector.  Only the cells the
      ! old icellU loop wrote are set and the rest stay zero, so that branch
      ! receives a byte-identical array.
      real (kind=dbl_kind), dimension (nx_block*ny_block,max_blocks) :: &
         prod      ! temporary array

      real (kind=dbl_kind), dimension(max_blocks) :: &
         dot       ! temporary scalar for accumulating the result

      character(len=*), parameter :: subname = '(global_dot_product)'

      prod = c0
      dot  = c0

      ! Slot order is icellU order, so dot(iblk) accumulates in exactly the
      ! sequence the block-layout loop used.  The summation is unchanged, not
      ! merely equivalent -- which is what keeps both branches bit-for-bit.
      !$OMP PARALLEL DO PRIVATE(is, iblk, p)
      do iblk = 1, nblocks
         do is = dynOff(iblk)+1, dynOff(iblk)+dynCnt(iblk)
            if (skipU(is)) cycle
            p = vector1_x(is)*vector2_x(is) + vector1_y(is)*vector2_y(is)
            prod(dynIJ(is),iblk) = p
            dot(iblk) = dot(iblk) + p
         enddo ! is
      enddo
      !$OMP END PARALLEL DO

      ! Use faster local summation result for several bfbflag settings.
      ! The local implementation sums over each block, sums over local
      ! blocks, and calls global_sum on a scalar and should be just as accurate as
      ! bfbflag = 'off', 'lsum8', and 'lsum4' without the extra copies and overhead
      ! in the more general array global_sum.  But use the array global_sum
      ! if bfbflag is more strict or for tripole grids (requires special masking)
      if (ns_boundary_type /= 'tripole' .and. ns_boundary_type /= 'tripoleT' .and. &
         (bfbflag == 'off' .or. bfbflag == 'lsum8' .or. bfbflag == 'lsum4')) then
         dot_product = global_sum(sum(dot), distrb_info)
      else
         ! prod is flattened; global_sum has no flat specific the way
         ! global_maxval and global_minval do, so reshape it -- exact, since
         ! the storage order is the same.  A global_sum_dbl_flat would be the
         ! better fix if this branch ever matters for speed.
         dot_product = global_sum(reshape(prod,(/nx_block,ny_block,max_blocks/)), &
                                  distrb_info, field_loc_NEcorner)
      endif

      end function global_dot_product

!=======================================================================

! Convert a grid function (tpu,tpv) to a one dimensional vector

      subroutine arrays_to_vec (ntot                , &
                                tpu     , tpv       , &
                                outvec)

      integer (kind=int_kind), intent(in) :: &
         ntot                  ! size of problem for Anderson

      ! compressed: the whole local slot vector, every block
      real (kind=dbl_kind), dimension (dynTot), intent(in) :: &
         tpu     , & ! x-component of vector
         tpv         ! y-component of vector

      real (kind=dbl_kind), dimension (ntot), intent(out) :: &
         outvec      ! output 1D vector

      ! local variables

      integer (kind=int_kind) :: &
         iblk, tot, is

      character(len=*), parameter :: subname = '(arrays_to_vec)'

      !-----------------------------------------------------------------
      ! form vector (converts from max_blocks arrays to single vector)
      !-----------------------------------------------------------------

      outvec(:) = c0
      tot = 0

      ! Slot order is icellU order, so tot advances in exactly the sequence
      ! the block-layout loop used and outvec comes out byte-identical.
      do iblk = 1, nblocks
         do is = dynOff(iblk)+1, dynOff(iblk)+dynCnt(iblk)
            if (skipU(is)) cycle
            tot = tot + 1
            outvec(tot) = tpu(is)
            tot = tot + 1
            outvec(tot) = tpv(is)
         enddo
      enddo ! is

      end subroutine arrays_to_vec

!=======================================================================

! Convert one dimensional vector to a grid function (tpu,tpv)

      subroutine vec_to_arrays (nx_block, ny_block  , &
                                nblocks , max_blocks, &
                                icellU  , ntot      , &
                                indxUi  , indxUj    , &
                                invec   ,             &
                                tpu     , tpv)

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         nblocks,            & ! nb of blocks
         max_blocks,         & ! max nb of blocks
         ntot                  ! size of problem for Anderson

      integer (kind=int_kind), dimension (max_blocks), intent(in) :: &
         icellU

      integer (kind=int_kind), dimension (nx_block*ny_block, max_blocks), intent(in) :: &
         indxUi  , & ! compressed index in i-direction
         indxUj      ! compressed index in j-direction

      real (kind=dbl_kind), dimension (ntot), intent(in) :: &
         invec       ! input 1D vector

      real (kind=dbl_kind), dimension (nx_block*ny_block,max_blocks), intent(out) :: &
         tpu     , & ! x-component of vector
         tpv         ! y-component of vector

      ! local variables

      integer (kind=int_kind) :: &
         i, j, iblk, tot, ij

      character(len=*), parameter :: subname = '(vec_to_arrays)'

      !-----------------------------------------------------------------
      ! form arrays (converts from vector to the max_blocks arrays)
      !-----------------------------------------------------------------

      tpu(:,:) = c0
      tpv(:,:) = c0
      tot = 0

      do iblk = 1, nblocks
         do ij = 1, icellU(iblk)
            i = indxUi(ij, iblk)
            j = indxUj(ij, iblk)
            tot = tot + 1
            tpu((j-1)*nx_block+i,iblk) = invec(tot)
            tot = tot + 1
            tpv((j-1)*nx_block+i,iblk) = invec(tot)
         enddo
      enddo! ij

      end subroutine vec_to_arrays

!=======================================================================

! Compressed counterpart of vec_to_arrays: scatters the Anderson vector
! straight into slot vectors.  The traversal is identical, so invec is
! consumed in the same order and only the destination layout differs.
! vec_to_arrays itself stays, because uvel/vvel are read by the halo,
! advection, history and restart in the block layout.

      subroutine vec_to_slots (ntot, invec, tpu, tpv)

      integer (kind=int_kind), intent(in) :: &
         ntot        ! size of problem for Anderson

      real (kind=dbl_kind), dimension (ntot), intent(in) :: &
         invec       ! input 1D vector

      ! compressed: the whole local slot vector, every block
      real (kind=dbl_kind), dimension (dynTot), intent(out) :: &
         tpu     , & ! x-component of vector
         tpv         ! y-component of vector

      ! local variables

      integer (kind=int_kind) :: &
         iblk, tot, is

      character(len=*), parameter :: subname = '(vec_to_slots)'

      tpu(:) = c0
      tpv(:) = c0

      tot = 0
      do iblk = 1, nblocks
         do is = dynOff(iblk)+1, dynOff(iblk)+dynCnt(iblk)
            if (skipU(is)) cycle
            tot = tot + 1
            tpu(is) = invec(tot)
            tot = tot + 1
            tpv(is) = invec(tot)
         enddo
      enddo ! is

      end subroutine vec_to_slots

!=======================================================================

! Update Q and R factors after deletion of the 1st column of G_diff
!
! author: P. Blain ECCC
!
! adapted from :
! H. F. Walker, “Anderson Acceleration: Algorithms and Implementations”
!   [Online]. Available: https://users.wpi.edu/~walker/Papers/anderson_accn_algs_imps.pdf

      subroutine qr_delete(Q, R)

      real (kind=dbl_kind), intent(inout) :: &
         Q(:,:),  & ! Q factor
         R(:,:)     ! R factor

      ! local variables

      integer (kind=int_kind) :: &
         i, j, k, & ! loop indices
         m, n       ! size of Q matrix

      real (kind=dbl_kind) :: &
         temp, c, s

      character(len=*), parameter :: subname = '(qr_delete)'

      n = size(Q, 1)
      m = size(Q, 2)
      do i = 1, m-1
         temp = sqrt(R(i, i+1)**2 + R(i+1, i+1)**2)
         c = R(i  , i+1) / temp
         s = R(i+1, i+1) / temp
         R(i  , i+1) = temp
         R(i+1, i+1) = 0
         if (i < m-1) then
            do j = i+2, m
               temp      =  c*R(i, j) + s*R(i+1, j)
               R(i+1, j) = -s*R(i, j) + c*R(i+1, j)
               R(i  , j) = temp
            enddo
         endif
         do k = 1, n
            temp      =  c*Q(k, i) + s*Q(k, i+1);
            Q(k, i+1) = -s*Q(k, i) + c*Q(k, i+1);
            Q(k, i)   = temp
         enddo
      enddo
      R(:, 1:m-1) = R(:, 2:m)

      end subroutine qr_delete

!=======================================================================

! FGMRES: Flexible generalized minimum residual method (with restarts).
! Solves the linear system A x = b using GMRES with a varying (right) preconditioner
!
! authors: Stéphane Gaudreault, Abdessamad Qaddouri, Philippe Blain, ECCC

      subroutine fgmres (Cb       , vrel    , &
                         umassdti ,           &
                         halo_info_mask     , &
                         bx       , by      , &
                         diagx    , diagy   , &
                         tolerance, maxinner, &
                         maxouter ,           &
                         solx     , soly    , &
                         nbiter)

      use ice_boundary, only: ice_HaloUpdate3DR8_compressed
      use ice_domain, only: maskhalo_dyn, halo_info
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bound
      use ice_dyn_shared, only: dxhy, dyhx, cxp, cyp, cxm, cym

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         umassdti  ! mass of U-cell/dte (kg/m^2 s)

      real (kind=dbl_kind), dimension(dynTot), intent(in) :: &
         vrel  , & ! coefficient for tauw
         Cb        ! seabed stress coefficient

      type (ice_halo), intent(in) :: &
         halo_info_mask !  ghost cell update info for masked halo

      real (kind=dbl_kind), dimension(dynTot), intent(in) :: &
         bx       , & ! Right hand side of the linear system (x components)
         by       , & ! Right hand side of the linear system (y components)
         diagx    , & ! Diagonal of the system matrix (x components)
         diagy        ! Diagonal of the system matrix (y components)

      real (kind=dbl_kind), intent(in) :: &
         tolerance   ! Tolerance to achieve. The algorithm terminates when the relative
                     ! residual is below tolerance

      integer (kind=int_kind), intent(in) :: &
         maxinner, & ! Restart the method every maxinner inner (Arnoldi) iterations
         maxouter    ! Maximum number of outer (restarts) iterations
                     ! Iteration will stop after maxinner*maxouter Arnoldi steps
                     ! even if the specified tolerance has not been achieved

      real (kind=dbl_kind), dimension(dynTot), intent(inout) :: &
         solx     , & ! Initial guess on input, approximate solution on output (x components)
         soly         ! Initial guess on input, approximate solution on output (y components)

      integer (kind=int_kind), intent(out) :: &
         nbiter      ! Total number of Arnoldi iterations performed

      ! local variables

      integer (kind=int_kind) :: &
         iblk    , & ! block index
         is      , & ! slot index
         j           ! Givens-rotation / Hessenberg loop index

      real (kind=dbl_kind), dimension(dynTot) :: &
         workspace_x , & ! work vector (x components)
         workspace_y     ! work vector (y components)

      real (kind=dbl_kind), dimension (max_blocks) :: &
         norm_squared   ! array to accumulate squared norm of grid function over blocks

      ! Allocatable rather than automatic.  At dim_fgmres = 50 these four come
      ! to about 12 MB, which overflows a default 8 MB stack on entry to the
      ! routine.  The fault is then reported against the first executable
      ! statement, which reads like a code bug rather than a stack limit.
      ! Moving them to the heap changes no value and no operation order, so
      ! results are unaffected.  Local allocatables are deallocated
      ! automatically on return, so the early returns below need no cleanup.
      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         arnoldi_basis_x , & ! Arnoldi basis (x components)
         arnoldi_basis_y , & ! Arnoldi basis (y components)
         orig_basis_x    , & ! original basis (x components)
         orig_basis_y        ! original basis (y components)

      integer (kind=int_kind) :: &
         ierr                ! allocate status

      real (kind=dbl_kind) :: &
         norm_residual   , & ! current L^2 norm of residual vector
         inverse_norm    , & ! inverse of the norm of a vector
         norm_rhs        , & ! L^2 norm of right-hand-side vector
         nu, t               ! local temporary values

      integer (kind=int_kind) :: &
         initer          , & ! inner (Arnoldi) loop counter
         outiter         , & ! outer (restarts) loop counter
         nextit          , & ! nextit == initer+1
         it, k, ii, jj       ! reusable loop counters

      real (kind=dbl_kind), dimension(maxinner+1) :: &
         rot_cos         , & ! cosine elements of Givens rotations
         rot_sin         , & ! sine elements of Givens rotations
         rhs_hess            ! right hand side vector of the Hessenberg (least squares) system

      real (kind=dbl_kind), dimension(maxinner+1, maxinner) :: &
         hessenberg        ! system matrix of the Hessenberg (least squares) system

      character (len=char_len) :: &
         precond_type ! type of preconditioner

      real (kind=dbl_kind) :: &
         relative_tolerance  ! relative_tolerance, i.e. tolerance*norm(initial residual)

      character(len=*), parameter :: subname = '(fgmres)'

      ! Here we go !

      ! Initialize
      outiter = 0
      nbiter = 0

      norm_squared = c0
      precond_type = precond

      allocate(arnoldi_basis_x(dynTot, maxinner+1), &
               arnoldi_basis_y(dynTot, maxinner+1), &
               orig_basis_x   (dynTot, maxinner  ), &
               orig_basis_y   (dynTot, maxinner  ), &
               stat=ierr)
      if (ierr/=0) call abort_ice(subname//': Out of memory')

      ! Cells with no ice should be zero-initialized
      workspace_x = c0
      workspace_y = c0
      arnoldi_basis_x = c0
      arnoldi_basis_y = c0

      ! solution is zero if RHS is zero
      norm_rhs = global_norm(                             bx      , by        )
      if (norm_rhs == c0) then
         solx = bx
         soly = by
         return
      endif

      ! Residual of the initial iterate

      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks
         call matvec (dynCnt(iblk), skipT(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),                                                &
                      skipU(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)), dynNbr(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk),:),          &
                      dxT_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dyT_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),   &
                      dxhy_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dyhx_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                      cxp_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), cyp_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),   &
                      cxm_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), cym_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),   &
                      solx(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), soly(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),           &
                      vrel_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), Cb_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),         &
                      zetax2_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:),                                                       &
                      etax2_s (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:),                                                       &
                      umassdti_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), fmU_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),    &
                      uarear_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),                                                      &
                      workspace_x(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), workspace_y(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call residual_vec (dynCnt(iblk), skipU(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),                                                    &
                            bx(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), by(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),                   &
                            workspace_x(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), workspace_y(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                            arnoldi_basis_x(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,1),                                                    &
                            arnoldi_basis_y(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,1))
      enddo
      !$OMP END PARALLEL DO

      ! Start outer (restarts) loop
      do
         ! Compute norm of initial residual
         norm_residual = global_norm(                                     arnoldi_basis_x(:,1), &
                                     arnoldi_basis_y(:,1))

         if (my_task == master_task .and. monitor_fgmres) then
            write(nu_diag, '(a,i4,a,d26.16)') "monitor_fgmres: iter_fgmres= ", nbiter, &
                                              " fgmres_L2norm= ", norm_residual
         endif

         ! Current guess is a good enough solution TODO: reactivate and test this
         ! if (norm_residual < tolerance) then
         !    return
         ! end if

         ! Normalize the first Arnoldi vector
         inverse_norm = c1 / norm_residual
         !$OMP PARALLEL DO PRIVATE(iblk,is)
         do iblk = 1, nblocks
            do is = dynOff(iblk)+1, dynOff(iblk)+dynCnt(iblk)
                  if (skipU(is)) cycle

               arnoldi_basis_x(is,1) = arnoldi_basis_x(is,1) * inverse_norm
               arnoldi_basis_y(is,1) = arnoldi_basis_y(is,1) * inverse_norm
            enddo ! ij
         enddo
         !$OMP END PARALLEL DO

         if (outiter == 0) then
            relative_tolerance = tolerance * norm_residual
         end if

         ! Initialize 1-st term of RHS of Hessenberg system
         rhs_hess(1)  = norm_residual
         rhs_hess(2:) = c0

         initer = 0

         ! Start of inner (Arnoldi) loop
         do

            nbiter = nbiter + 1
            initer = initer + 1
            nextit = initer + 1
            ! precondition the current Arnoldi vector
            call precondition(Cb          , vrel           , &
                              umassdti    ,                  &
                              halo_info_mask,                &
                              arnoldi_basis_x(:,initer),     &
                              arnoldi_basis_y(:,initer),     &
                              diagx       , diagy          , &
                              precond_type,                  &
                              workspace_x , workspace_y)
            orig_basis_x(:,initer) = workspace_x
            orig_basis_y(:,initer) = workspace_y

            ! Update workspace with boundary values
            call stack_fields(workspace_x, workspace_y, fld2c)
            call ice_timer_start(timer_bound)
            if (maskhalo_dyn) then
               call ice_HaloUpdate3DR8_compressed(fld2c, halo_info_mask, &
                                    field_loc_NEcorner, field_type_vector, dynMap, dynOff)
            else
               call ice_HaloUpdate3DR8_compressed(fld2c, halo_info, &
                                    field_loc_NEcorner, field_type_vector, dynMap, dynOff)
            endif
            call ice_timer_stop(timer_bound)
            call unstack_fields(fld2c, workspace_x, workspace_y)

            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks
               call matvec (dynCnt(iblk), skipT(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),                                                    &
                            skipU(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)), dynNbr(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk),:),              &
                            dxT_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dyT_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),       &
                            dxhy_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dyhx_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),     &
                            cxp_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), cyp_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),       &
                            cxm_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), cym_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),       &
                            workspace_x(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), workspace_y(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                            vrel_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), Cb_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),             &
                            zetax2_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:),                                                           &
                            etax2_s (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:),                                                           &
                            umassdti_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), fmU_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),        &
                            uarear_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),                                                          &
                            arnoldi_basis_x(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,nextit),                                               &
                            arnoldi_basis_y(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,nextit))
            enddo
            !$OMP END PARALLEL DO

            ! Orthogonalize the new vector
            call orthogonalize(ortho_type     , initer         , &
                               nextit         , maxinner       , &
                               arnoldi_basis_x, arnoldi_basis_y, &
                               hessenberg)

            ! Compute norm of new Arnoldi vector and update Hessenberg matrix
            hessenberg(nextit,initer) = global_norm(                                                    arnoldi_basis_x(:,nextit), &
                                                    arnoldi_basis_y(:,nextit))

            ! Watch out for happy breakdown
            if (.not. almost_zero( hessenberg(nextit,initer) ) ) then
               ! Normalize next Arnoldi vector
               inverse_norm = c1 / hessenberg(nextit,initer)
               !$OMP PARALLEL DO PRIVATE(iblk,is)
               do iblk = 1, nblocks
                  do is = dynOff(iblk)+1, dynOff(iblk)+dynCnt(iblk)
                  if (skipU(is)) cycle

                     arnoldi_basis_x(is,nextit) = arnoldi_basis_x(is,nextit)*inverse_norm
                     arnoldi_basis_y(is,nextit) = arnoldi_basis_y(is,nextit)*inverse_norm
                  enddo ! ij
               enddo
               !$OMP END PARALLEL DO
            end if

            ! Apply previous Givens rotation to the last column of the Hessenberg matrix
            if (initer > 1) then
               do k = 2, initer
                  t = hessenberg(k-1, initer)
                  hessenberg(k-1, initer) =  rot_cos(k-1)*t + rot_sin(k-1)*hessenberg(k, initer)
                  hessenberg(k,   initer) = -rot_sin(k-1)*t + rot_cos(k-1)*hessenberg(k, initer)
               end do
            end if

            ! Compute and apply new Givens rotation
            nu = sqrt(hessenberg(initer,initer)**2 + hessenberg(nextit,initer)**2)
            if (.not. almost_zero(nu)) then
               rot_cos(initer) = hessenberg(initer,initer) / nu
               rot_sin(initer) = hessenberg(nextit,initer) / nu

               rhs_hess(nextit) = -rot_sin(initer) * rhs_hess(initer)
               rhs_hess(initer) =  rot_cos(initer) * rhs_hess(initer)

               hessenberg(initer,initer) = rot_cos(initer) * hessenberg(initer,initer) + rot_sin(initer) * hessenberg(nextit,initer)
            end if

            ! Check for convergence
            norm_residual = abs(rhs_hess(nextit))

            if (my_task == master_task .and. monitor_fgmres) then
               write(nu_diag, '(a,i4,a,d26.16)') "monitor_fgmres: iter_fgmres= ", nbiter, &
                                                 " fgmres_L2norm= ", norm_residual
            endif

             if ((initer >= maxinner) .or. (norm_residual <= relative_tolerance)) then
               exit
            endif

         end do ! end of inner (Arnoldi) loop

         ! At this point either the maximum number of inner iterations
         ! was reached or the absolute residual is below the scaled tolerance.

         ! Solve the (now upper triangular) system "hessenberg * sol_hess = rhs_hess"
         ! (sol_hess is stored in rhs_hess)
         rhs_hess(initer) = rhs_hess(initer) / hessenberg(initer,initer)
         do ii = 2, initer
            k  = initer - ii + 1
            t  = rhs_hess(k)
            do j = k + 1, initer
               t = t - hessenberg(k,j) * rhs_hess(j)
            end do
            rhs_hess(k) = t / hessenberg(k,k)
         end do

         ! Form linear combination to get new solution iterate
         do it = 1, initer
            t = rhs_hess(it)
            !$OMP PARALLEL DO PRIVATE(iblk,is)
            do iblk = 1, nblocks
               do is = dynOff(iblk)+1, dynOff(iblk)+dynCnt(iblk)
                  if (skipU(is)) cycle

                  solx(is) = solx(is) + t * orig_basis_x(is,it)
                  soly(is) = soly(is) + t * orig_basis_y(is,it)
               enddo ! ij
            enddo
            !$OMP END PARALLEL DO
         end do

         ! Increment outer loop counter and check for convergence
         outiter = outiter + 1
         if (norm_residual <= relative_tolerance .or. outiter >= maxouter) then
            return
         end if

         ! Solution is not convergent : compute residual vector and continue.

         ! The residual vector is computed here using (see Saad p. 177) :
         ! \begin{equation}
         !    r =  V_{m+1} * Q_m^T * (\gamma_{m+1} * e_{m+1})
         ! \end{equation}
         ! where :
         ! $r$ is the residual
         ! $V_{m+1}$ is a matrix whose columns are the Arnoldi vectors from 1 to nextit (m+1)
         ! $Q_m$ is the product of the Givens rotation : Q_m = G_m G_{m-1} ... G_1
         ! $gamma_{m+1}$ is the last element of rhs_hess
         ! $e_{m+1})$ is the unit vector (0, 0, ..., 1)^T \in \reals^{m+1}

         ! Apply the Givens rotation in reverse order to g := \gamma_{m+1} * e_{m+1},
         ! store the result in rhs_hess
         do it = 1, initer
            jj = nextit - it + 1
            rhs_hess(jj-1) = -rot_sin(jj-1) * rhs_hess(jj) ! + rot_cos(jj-1) * g(jj-1) (== 0)
            rhs_hess(jj)   =  rot_cos(jj-1) * rhs_hess(jj) ! + rot_sin(jj-1) * g(jj-1) (== 0)
         end do

         ! Compute the residual by multiplying V_{m+1} and rhs_hess
         workspace_x = c0
         workspace_y = c0
         do it = 1, nextit
            !$OMP PARALLEL DO PRIVATE(iblk,is)
            do iblk = 1, nblocks
               do is = dynOff(iblk)+1, dynOff(iblk)+dynCnt(iblk)
                  if (skipU(is)) cycle

                  workspace_x(is) = workspace_x(is) + rhs_hess(it) * arnoldi_basis_x(is,it)
                  workspace_y(is) = workspace_y(is) + rhs_hess(it) * arnoldi_basis_y(is,it)
               enddo ! ij
            enddo
            !$OMP END PARALLEL DO
            arnoldi_basis_x(:,1) = workspace_x
            arnoldi_basis_y(:,1) = workspace_y
         end do
      end do ! end of outer (restarts) loop

      end subroutine fgmres

!=======================================================================

! PGMRES: Right-preconditioned generalized minimum residual method (with restarts).
! Solves the linear A x = b using GMRES with a right preconditioner
!
! authors: Stéphane Gaudreault, Abdessamad Qaddouri, Philippe Blain, ECCC

      subroutine pgmres (Cb       , vrel    , &
                         umassdti ,           &
                         halo_info_mask,      &
                         bx       , by      , &
                         diagx    , diagy   , &
                         tolerance, maxinner, &
                         maxouter ,           &
                         solx     , soly    , &
                         nbiter)

      use ice_boundary, only: ice_HaloUpdate3DR8_compressed
      use ice_domain, only: maskhalo_dyn, halo_info
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bound
      use ice_dyn_shared, only: dyhx, dxhy, cxp, cyp, cxm, cym

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         umassdti  ! mass of U-cell/dte (kg/m^2 s)

      real (kind=dbl_kind), dimension(dynTot), intent(in) :: &
         vrel  , & ! coefficient for tauw
         Cb        ! seabed stress coefficient

      type (ice_halo), intent(in) :: &
         halo_info_mask !  ghost cell update info for masked halo

      real (kind=dbl_kind), dimension(dynTot), intent(in) :: &
         bx       , & ! Right hand side of the linear system (x components)
         by           ! Right hand side of the linear system (y components)

      real (kind=dbl_kind), dimension(dynTot), intent(in) :: &
         diagx    , & ! Diagonal of the system matrix (x components)
         diagy        ! Diagonal of the system matrix (y components)

      real (kind=dbl_kind), intent(in) :: &
         tolerance   ! Tolerance to achieve. The algorithm terminates when the relative
                     ! residual is below tolerance

      integer (kind=int_kind), intent(in) :: &
         maxinner, & ! Restart the method every maxinner inner (Arnoldi) iterations
         maxouter    ! Maximum number of outer (restarts) iterations
                     ! Iteration will stop after maxinner*maxouter Arnoldi steps
                     ! even if the specified tolerance has not been achieved

      real (kind=dbl_kind), dimension(dynTot), intent(inout) :: &
         solx     , & ! Initial guess on input, approximate solution on output (x components)
         soly         ! Initial guess on input, approximate solution on output (y components)

      integer (kind=int_kind), intent(out) :: &
         nbiter      ! Total number of Arnoldi iterations performed

      ! local variables

      integer (kind=int_kind) :: &
         iblk    , & ! block index
         is      , & ! slot index
         j           ! Givens-rotation / Hessenberg loop index

      real (kind=dbl_kind), dimension(dynTot) :: &
         workspace_x , & ! work vector (x components)
         workspace_y     ! work vector (y components)

      real (kind=dbl_kind), dimension (max_blocks) :: &
         norm_squared   ! array to accumulate squared norm of grid function over blocks

      ! Allocatable for the same reason as in fgmres.  dim_pgmres is small (5)
      ! so this pair is well under a megabyte today, but it is the same array
      ! with the same failure mode if that namelist value is ever raised.
      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         arnoldi_basis_x , & ! Arnoldi basis (x components)
         arnoldi_basis_y     ! Arnoldi basis (y components)

      integer (kind=int_kind) :: &
         ierr                ! allocate status

      real (kind=dbl_kind) :: &
         norm_residual   , & ! current L^2 norm of residual vector
         inverse_norm    , & ! inverse of the norm of a vector
         nu, t               ! local temporary values

      integer (kind=int_kind) :: &
         initer          , & ! inner (Arnoldi) loop counter
         outiter         , & ! outer (restarts) loop counter
         nextit          , & ! nextit == initer+1
         it, k, ii, jj       ! reusable loop counters

      real (kind=dbl_kind), dimension(maxinner+1) :: &
         rot_cos         , & ! cosine elements of Givens rotations
         rot_sin         , & ! sine elements of Givens rotations
         rhs_hess            ! right hand side vector of the Hessenberg (least squares) system

      real (kind=dbl_kind), dimension(maxinner+1, maxinner) :: &
         hessenberg          ! system matrix of the Hessenberg (least squares) system

      character(len=char_len) :: &
         precond_type    , & ! type of preconditioner
         ortho_type          ! type of orthogonalization

      real (kind=dbl_kind) :: &
         relative_tolerance  ! relative_tolerance, i.e. tolerance*norm(initial residual)

      character(len=*), parameter :: subname = '(pgmres)'

      ! Here we go !

      ! Initialize
      outiter = 0
      nbiter = 0

      norm_squared = c0
      precond_type = 'diag' ! Jacobi preconditioner
      ortho_type = 'cgs' ! classical gram-schmidt TODO: try with MGS

      allocate(arnoldi_basis_x(dynTot, maxinner+1), &
               arnoldi_basis_y(dynTot, maxinner+1), &
               stat=ierr)
      if (ierr/=0) call abort_ice(subname//': Out of memory')

      ! Cells with no ice should be zero-initialized
      workspace_x = c0
      workspace_y = c0
      arnoldi_basis_x = c0
      arnoldi_basis_y = c0

      ! Residual of the initial iterate

      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks
         call matvec (dynCnt(iblk), skipT(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),                                                &
                      skipU(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)), dynNbr(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk),:),          &
                      dxT_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dyT_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),   &
                      dxhy_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dyhx_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                      cxp_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), cyp_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),   &
                      cxm_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), cym_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),   &
                      solx(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), soly(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),           &
                      vrel_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), Cb_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),         &
                      zetax2_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:),                                                       &
                      etax2_s (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:),                                                       &
                      umassdti_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), fmU_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),    &
                      uarear_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),                                                      &
                      workspace_x(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), workspace_y(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         call residual_vec (dynCnt(iblk), skipU(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),                                                    &
                            bx(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), by(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),                   &
                            workspace_x(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), workspace_y(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                            arnoldi_basis_x(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,1),                                                    &
                            arnoldi_basis_y(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,1))
      enddo
      !$OMP END PARALLEL DO

      ! Start outer (restarts) loop
      do
         ! Compute norm of initial residual
         norm_residual = global_norm(                                     arnoldi_basis_x(:,1), &
                                     arnoldi_basis_y(:,1))

         if (my_task == master_task .and. monitor_pgmres) then
            write(nu_diag, '(a,i4,a,d26.16)') "monitor_pgmres: iter_pgmres= ", nbiter, &
                                              " pgmres_L2norm= ", norm_residual
         endif

         ! Current guess is a good enough solution
         ! if (norm_residual < tolerance) then
         !    return
         ! end if

         ! Normalize the first Arnoldi vector
         inverse_norm = c1 / norm_residual
         !$OMP PARALLEL DO PRIVATE(iblk,is)
         do iblk = 1, nblocks
            do is = dynOff(iblk)+1, dynOff(iblk)+dynCnt(iblk)
                  if (skipU(is)) cycle

               arnoldi_basis_x(is,1) = arnoldi_basis_x(is,1) * inverse_norm
               arnoldi_basis_y(is,1) = arnoldi_basis_y(is,1) * inverse_norm
            enddo ! ij
         enddo
         !$OMP END PARALLEL DO

         if (outiter == 0) then
            relative_tolerance = tolerance * norm_residual
         end if

         ! Initialize 1-st term of RHS of Hessenberg system
         rhs_hess(1)  = norm_residual
         rhs_hess(2:) = c0

         initer = 0

         ! Start of inner (Arnoldi) loop
         do

            nbiter = nbiter + 1
            initer = initer + 1
            nextit = initer + 1

            ! precondition the current Arnoldi vector
            call precondition(Cb          , vrel           , &
                              umassdti    ,                  &
                              halo_info_mask,                &
                              arnoldi_basis_x(:,initer),     &
                              arnoldi_basis_y(:,initer),     &
                              diagx       , diagy          , &
                              precond_type,                  &
                              workspace_x , workspace_y)

            ! Update workspace with boundary values
            call stack_fields(workspace_x, workspace_y, fld2c)
            call ice_timer_start(timer_bound)
            if (maskhalo_dyn) then
               call ice_HaloUpdate3DR8_compressed(fld2c, halo_info_mask, &
                                    field_loc_NEcorner, field_type_vector, dynMap, dynOff)
            else
               call ice_HaloUpdate3DR8_compressed(fld2c, halo_info, &
                                    field_loc_NEcorner, field_type_vector, dynMap, dynOff)
            endif
            call ice_timer_stop(timer_bound)
            call unstack_fields(fld2c, workspace_x, workspace_y)

            !$OMP PARALLEL DO PRIVATE(iblk)
            do iblk = 1, nblocks
               call matvec (dynCnt(iblk), skipT(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)),                                                    &
                            skipU(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)), dynNbr(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk),:),              &
                            dxT_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dyT_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),       &
                            dxhy_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dyhx_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),     &
                            cxp_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), cyp_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),       &
                            cxm_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), cym_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),       &
                            workspace_x(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), workspace_y(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                            vrel_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), Cb_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),             &
                            zetax2_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:),                                                           &
                            etax2_s (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:),                                                           &
                            umassdti_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), fmU_s(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),        &
                            uarear_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),                                                          &
                            arnoldi_basis_x(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,nextit),                                               &
                            arnoldi_basis_y(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,nextit))
            enddo
            !$OMP END PARALLEL DO

            ! Orthogonalize the new vector
            call orthogonalize(ortho_type     , initer         , &
                               nextit         , maxinner       , &
                               arnoldi_basis_x, arnoldi_basis_y, &
                               hessenberg)

            ! Compute norm of new Arnoldi vector and update Hessenberg matrix
            hessenberg(nextit,initer) = global_norm(                                                    arnoldi_basis_x(:,nextit), &
                                                    arnoldi_basis_y(:,nextit))

            ! Watch out for happy breakdown
            if (.not. almost_zero( hessenberg(nextit,initer) ) ) then
               ! Normalize next Arnoldi vector
               inverse_norm = c1 / hessenberg(nextit,initer)
               !$OMP PARALLEL DO PRIVATE(iblk,is)
               do iblk = 1, nblocks
                  do is = dynOff(iblk)+1, dynOff(iblk)+dynCnt(iblk)
                  if (skipU(is)) cycle

                     arnoldi_basis_x(is,nextit) = arnoldi_basis_x(is,nextit)*inverse_norm
                     arnoldi_basis_y(is,nextit) = arnoldi_basis_y(is,nextit)*inverse_norm
                  enddo ! ij
               enddo
               !$OMP END PARALLEL DO
            end if

            ! Apply previous Givens rotation to the last column of the Hessenberg matrix
            if (initer > 1) then
               do k = 2, initer
                  t = hessenberg(k-1, initer)
                  hessenberg(k-1, initer) =  rot_cos(k-1)*t + rot_sin(k-1)*hessenberg(k, initer)
                  hessenberg(k,   initer) = -rot_sin(k-1)*t + rot_cos(k-1)*hessenberg(k, initer)
               end do
            end if

            ! Compute and apply new Givens rotation
            nu = sqrt(hessenberg(initer,initer)**2 + hessenberg(nextit,initer)**2)
            if (.not. almost_zero(nu)) then
               rot_cos(initer) = hessenberg(initer,initer) / nu
               rot_sin(initer) = hessenberg(nextit,initer) / nu

               rhs_hess(nextit) = -rot_sin(initer) * rhs_hess(initer)
               rhs_hess(initer) =  rot_cos(initer) * rhs_hess(initer)

               hessenberg(initer,initer) = rot_cos(initer) * hessenberg(initer,initer) + rot_sin(initer) * hessenberg(nextit,initer)
            end if

            ! Check for convergence
            norm_residual = abs(rhs_hess(nextit))

            if (my_task == master_task .and. monitor_pgmres) then
               write(nu_diag, '(a,i4,a,d26.16)') "monitor_pgmres: iter_pgmres= ", nbiter, &
                                                 " pgmres_L2norm= ", norm_residual
            endif

             if ((initer >= maxinner) .or. (norm_residual <= relative_tolerance)) then
               exit
            endif

         end do ! end of inner (Arnoldi) loop

         ! At this point either the maximum number of inner iterations
         ! was reached or the absolute residual is below the scaled tolerance.

         ! Solve the (now upper triangular) system "hessenberg * sol_hess = rhs_hess"
         ! (sol_hess is stored in rhs_hess)
         rhs_hess(initer) = rhs_hess(initer) / hessenberg(initer,initer)
         do ii = 2, initer
            k  = initer - ii + 1
            t  = rhs_hess(k)
            do j = k + 1, initer
               t = t - hessenberg(k,j) * rhs_hess(j)
            end do
            rhs_hess(k) = t / hessenberg(k,k)
         end do

         ! Form linear combination to get new solution iterate
         workspace_x = c0
         workspace_y = c0
         do it = 1, initer
            t = rhs_hess(it)
            !$OMP PARALLEL DO PRIVATE(iblk,is)
            do iblk = 1, nblocks
               do is = dynOff(iblk)+1, dynOff(iblk)+dynCnt(iblk)
                  if (skipU(is)) cycle

                  workspace_x(is) = workspace_x(is) + t * arnoldi_basis_x(is,it)
                  workspace_y(is) = workspace_y(is) + t * arnoldi_basis_y(is,it)
               enddo ! ij
            enddo
            !$OMP END PARALLEL DO
         end do

         ! Call preconditioner
         call precondition(Cb          , vrel       , &
                           umassdti    ,              &
                           halo_info_mask,            &
                           workspace_x , workspace_y, &
                           diagx       , diagy      , &
                           precond_type,              &
                           workspace_x , workspace_y)

         solx = solx + workspace_x
         soly = soly + workspace_y

         ! Increment outer loop counter and check for convergence
         outiter = outiter + 1
         if (norm_residual <= relative_tolerance .or. outiter >= maxouter) then
            return
         end if

         ! Solution is not convergent : compute residual vector and continue.

         ! The residual vector is computed here using (see Saad p. 177) :
         ! \begin{equation}
         !    r =  V_{m+1} * Q_m^T * (\gamma_{m+1} * e_{m+1})
         ! \end{equation}
         ! where :
         ! $r$ is the residual
         ! $V_{m+1}$ is a matrix whose columns are the Arnoldi vectors from 1 to nextit (m+1)
         ! $Q_m$ is the product of the Givens rotation : Q_m = G_m G_{m-1} ... G_1
         ! $gamma_{m+1}$ is the last element of rhs_hess
         ! $e_{m+1})$ is the unit vector (0, 0, ..., 1)^T \in \reals^{m+1}

         ! Apply the Givens rotation in reverse order to g := \gamma_{m+1} * e_{m+1},
         ! store the result in rhs_hess
         do it = 1, initer
            jj = nextit - it + 1
            rhs_hess(jj-1) = -rot_sin(jj-1) * rhs_hess(jj) ! + rot_cos(jj-1) * g(jj-1) (== 0)
            rhs_hess(jj)   =  rot_cos(jj-1) * rhs_hess(jj) ! + rot_sin(jj-1) * g(jj-1) (== 0)
         end do

         ! Compute the residual by multiplying V_{m+1} and rhs_hess
         workspace_x = c0
         workspace_y = c0
         do it = 1, nextit
            !$OMP PARALLEL DO PRIVATE(iblk,is)
            do iblk = 1, nblocks
               do is = dynOff(iblk)+1, dynOff(iblk)+dynCnt(iblk)
                  if (skipU(is)) cycle

                  workspace_x(is) = workspace_x(is) + rhs_hess(it) * arnoldi_basis_x(is,it)
                  workspace_y(is) = workspace_y(is) + rhs_hess(it) * arnoldi_basis_y(is,it)
               enddo ! ij
            enddo
            !$OMP END PARALLEL DO
            arnoldi_basis_x(:,1) = workspace_x
            arnoldi_basis_y(:,1) = workspace_y
         end do
      end do ! end of outer (restarts) loop

      end subroutine pgmres

!=======================================================================

! Generic routine to precondition a vector
!
! authors: Philippe Blain, ECCC

      subroutine precondition(Cb          , vrel , &
                              umassdti    ,        &
                              halo_info_mask,      &
                              vx          , vy   , &
                              diagx       , diagy, &
                              precond_type,        &
                              wx          , wy)

      ! compressed: the whole local slot vector, every block
      real (kind=dbl_kind), dimension(dynTot), intent(in) :: &
         umassdti  ! mass of U-cell/dte (kg/m^2 s)

      real (kind=dbl_kind), dimension(dynTot), intent(in) :: &
         vrel  , & ! coefficient for tauw
         Cb        ! seabed stress coefficient

      type (ice_halo), intent(in) :: &
         halo_info_mask !  ghost cell update info for masked halo

      real (kind=dbl_kind), dimension (dynTot), intent(in) :: &
         vx       , & ! input vector (x components)
         vy           ! input vector (y components)

      real (kind=dbl_kind), dimension (dynTot), intent(in) :: &
         diagx    , & ! diagonal of the system matrix (x components)
         diagy        ! diagonal of the system matrix (y components)

      character (len=char_len), intent(in) :: &
         precond_type ! type of preconditioner

      real (kind=dbl_kind), dimension (dynTot), intent(inout) :: &
         wx       , & ! preconditionned vector (x components)
         wy           ! preconditionned vector (y components)

      ! local variables

      integer (kind=int_kind) :: &
         iblk    , & ! block index
         is      , & ! slot index
         j           ! Givens-rotation / Hessenberg loop index

      real (kind=dbl_kind) :: &
         tolerance   ! Tolerance for PGMRES

      integer (kind=int_kind) :: &
         maxinner    ! Restart parameter for PGMRES

      integer (kind=int_kind) :: &
         maxouter    ! Maximum number of outer iterations for PGMRES

      integer (kind=int_kind) :: &
         nbiter      ! Total number of iteration PGMRES performed

      character(len=*), parameter :: subname = '(precondition)'

      if     (precond_type == 'ident') then ! identity (no preconditioner)
         wx = vx
         wy = vy
      elseif (precond_type == 'diag') then ! Jacobi preconditioner (diagonal)
         !$OMP PARALLEL DO PRIVATE(iblk,is)
         do iblk = 1, nblocks
            do is = dynOff(iblk)+1, dynOff(iblk)+dynCnt(iblk)
               if (skipU(is)) cycle

               wx(is) = vx(is) / diagx(is)
               wy(is) = vy(is) / diagy(is)
            enddo ! is
         enddo
         !$OMP END PARALLEL DO
      elseif (precond_type == 'pgmres') then ! PGMRES (Jacobi-preconditioned GMRES)
         ! Initialize preconditioned vector to 0 ! TODO: try with wx = vx or vx/diagx
         wx = c0
         wy = c0
         tolerance = reltol_pgmres
         maxinner = dim_pgmres
         maxouter = maxits_pgmres
         call pgmres (Cb       , vrel    , &
                      umassdti ,           &
                      halo_info_mask     , &
                      vx       , vy      , &
                      diagx    , diagy   , &
                      tolerance, maxinner, &
                      maxouter ,           &
                      wx       , wy      , &
                      nbiter)
      else
         call abort_ice(error_message='wrong preconditioner in ' // subname, &
            file=__FILE__, line=__LINE__)
      endif
      end subroutine precondition

!=======================================================================

! Generic routine to orthogonalize a vector (arnoldi_basis_[xy](:, :, :, nextit))
! against a set of vectors (arnoldi_basis_[xy](:, :, :, 1:initer))
!
! authors: Philippe Blain, ECCC

      subroutine orthogonalize(ortho_type     , initer         , &
                               nextit         , maxinner       , &
                               arnoldi_basis_x, arnoldi_basis_y, &
                               hessenberg)

      character(len=*), intent(in) :: &
         ortho_type ! type of orthogonalization

      integer (kind=int_kind), intent(in) :: &
         initer  , & ! inner (Arnoldi) loop counter
         nextit  , & ! nextit == initer+1
         maxinner    ! Restart the method every maxinner inner iterations

      ! compressed: the whole local slot vector, every block
      real (kind=dbl_kind), dimension(dynTot, maxinner+1), intent(inout) :: &
         arnoldi_basis_x , & ! arnoldi basis (x components)
         arnoldi_basis_y     ! arnoldi basis (y components)

      real (kind=dbl_kind), dimension(maxinner+1, maxinner), intent(inout) :: &
         hessenberg        ! system matrix of the Hessenberg (least squares) system

      ! local variables

      integer (kind=int_kind) :: &
         it      , & ! reusable loop counter
         iblk    , & ! block index
         is          ! slot index

      character(len=*), parameter :: subname = '(orthogonalize)'

      if (trim(ortho_type) == 'cgs') then ! Classical Gram-Schmidt
         ! Classical Gram-Schmidt orthogonalisation process
         ! First loop of Gram-Schmidt (compute coefficients)
         do it = 1, initer
            hessenberg(it, initer) = global_dot_product(arnoldi_basis_x(:,it)    , &
                                                        arnoldi_basis_y(:,it)    , &
                                                        arnoldi_basis_x(:,nextit), &
                                                        arnoldi_basis_y(:,nextit))
         end do
         ! Second loop of Gram-Schmidt (orthonormalize)
         do it = 1, initer
            !$OMP PARALLEL DO PRIVATE(iblk,is)
            do iblk = 1, nblocks
               do is = dynOff(iblk)+1, dynOff(iblk)+dynCnt(iblk)
                  if (skipU(is)) cycle

                  arnoldi_basis_x(is,nextit) = arnoldi_basis_x(is,nextit) &
                                                        - hessenberg(it,initer) * arnoldi_basis_x(is,it)
                  arnoldi_basis_y(is,nextit) = arnoldi_basis_y(is,nextit) &
                                                        - hessenberg(it,initer) * arnoldi_basis_y(is,it)
               enddo ! is
            enddo
            !$OMP END PARALLEL DO
         end do
      elseif (trim(ortho_type) == 'mgs') then ! Modified Gram-Schmidt
         ! Modified Gram-Schmidt orthogonalisation process
         do it = 1, initer
            hessenberg(it, initer) = global_dot_product(arnoldi_basis_x(:,it)    , &
                                                        arnoldi_basis_y(:,it)    , &
                                                        arnoldi_basis_x(:,nextit), &
                                                        arnoldi_basis_y(:,nextit))

            !$OMP PARALLEL DO PRIVATE(iblk,is)
            do iblk = 1, nblocks
               do is = dynOff(iblk)+1, dynOff(iblk)+dynCnt(iblk)
                  if (skipU(is)) cycle

                  arnoldi_basis_x(is,nextit) = arnoldi_basis_x(is,nextit) &
                                                        - hessenberg(it,initer) * arnoldi_basis_x(is,it)
                  arnoldi_basis_y(is,nextit) = arnoldi_basis_y(is,nextit) &
                                                        - hessenberg(it,initer) * arnoldi_basis_y(is,it)
               enddo ! is
            enddo
            !$OMP END PARALLEL DO
         end do
      else
         call abort_ice(error_message='wrong orthonalization in ' // subname, &
         file=__FILE__, line=__LINE__)
      endif

   end subroutine orthogonalize

!=======================================================================

! Check if value A is close to zero, up to machine precision
!
!author
!     Stéphane Gaudreault, ECCC -- June 2014
!
!revision
!     v4-80 - Gaudreault S.         - gfortran compatibility
!     2019  - Philippe Blain, ECCC  - converted to CICE standards

      logical function almost_zero(A) result(retval)

      real (kind=dbl_kind), intent(in) :: A

      ! local variables

      character(len=*), parameter :: subname = '(almost_zero)'

      integer (kind=int8_kind) :: aBit
      integer (kind=int8_kind), parameter :: two_complement = int(Z'80000000', kind=int8_kind)
      aBit = 0
      aBit = transfer(A, aBit)
      if (aBit < 0) then
         aBit = two_complement - aBit
      end if
      ! lexicographic order test with a tolerance of 1 adjacent float
      retval = (abs(aBit) <= 1)

      end function almost_zero

!=======================================================================

      end module ice_dyn_vp

!=======================================================================
