! Module for 1d evp dynamics
! Mimics the 2d B grid solver
! functions in this module includes conversion from 1d to 2d and vice versa.
! cpp flag _OPENMP_TARGET is for gpu. Otherwize optimized for cpu
! FIXME: For now it allocates all water point, which in most cases could be avoided.
!===============================================================================
! Created by Till Rasmussen (DMI), Mads Hvid Ribergaard (DMI), and Jacob W. Poulsen, Intel

module ice_dyn_evp1d

  !- modules -------------------------------------------------------------------
  use ice_kinds_mod
  use ice_blocks, only: nx_block, ny_block, nghost
  use ice_constants
  use ice_communicate, only: my_task, master_task
  use ice_domain_size, only: max_blocks, nx_global, ny_global
  use ice_fileunits, only: nu_diag
  use ice_exit, only: abort_ice

  !- directives ----------------------------------------------------------------
  implicit none
  private

  !- public routines -----------------------------------------------------------
  public :: dyn_evp1d_init, dyn_evp1d_run, dyn_evp1d_finalize

  !- private routines ----------------------------------------------------------

  !- private vars --------------------------------------------------------------
  ! nx and ny are module variables for arrays after gather (G_*) Dimension according to CICE is
  ! nx_global+2*nghost, ny_global+2*nghost
  ! nactive are number of active points (both t and u). navel is number of active
  integer(kind=int_kind), save :: nx, ny, nActive, navel, nallocated

  ! indexes
  integer(kind=int_kind), allocatable, dimension(:,:) :: iwidx
  logical(kind=log_kind), allocatable, dimension(:)   :: skipTcell,skipUcell
  integer(kind=int_kind), allocatable, dimension(:)   :: ee,ne,se,nw,sw,sse ! arrays for neighbour points
  integer(kind=int_kind), allocatable, dimension(:)   :: indxti, indxtj, indxTij
  ! Scratch map over the gathered grid, linear index i+(j-1)*nx -> slot.
  ! -1 marks a cell the solver needs before slots are handed out, 0 a cell it
  ! does not. Built by calc_navel, consumed and released by convert_2d_1d_init,
  ! so it exists only during init and only on the master task.
  integer(kind=int_kind), allocatable, dimension(:)   :: ijslot

  ! 1D arrays to allocate

  ! Grid
  real   (kind=dbl_kind), allocatable, dimension(:)   :: &
     HTE_1d,HTN_1d, HTEm1_1d,HTNm1_1d, dxT_1d, dyT_1d, uarear_1d

  ! time varying
  real(kind=dbl_kind)   , allocatable, dimension(:)   ::                  &
    cdn_ocn,aiu,uocn,vocn,waterxU,wateryU,forcexU,forceyU,umassdti,fmU,   &
    strintxU,strintyU,uvel_init,vvel_init, strength, uvel, vvel,          &
    stressp_1, stressp_2, stressp_3, stressp_4, stressm_1, stressm_2,     &
    stressm_3, stressm_4, stress12_1, stress12_2, stress12_3, stress12_4, &
    str1, str2, str3, str4, str5, str6, str7, str8, Tbu, Cb

  ! Every boundary condition reduces to dst = s1 + w*(s1 - s2):
  !   cyclic         s1 is the matching cell on the opposite edge, w = 0
  !   zero_gradient  s1 is the cell just inside the edge,          w = 0
  !   linear_extrap  s1, s2 the two cells just inside,   w = the ghost depth
  ! One ordered list for all of them, built once in halo_sweep.
  integer(kind=int_kind), allocatable, dimension(:)   ::                       &
    halo_bc_dst, halo_bc_s1, halo_bc_s2
  real   (kind=dbl_kind), allocatable, dimension(:)   :: halo_bc_w
  integer(kind=int_kind)                              :: n_halo_bc

!=============================================================================
  contains
!=============================================================================
! module public subroutines
! In addition all water points are assumed to be active and allocated thereafter.
!=============================================================================

  subroutine dyn_evp1d_init

    use ice_grid, only: G_HTE, G_HTN

    implicit none

    ! local variables

    real(kind=dbl_kind)   , allocatable, dimension(:,:) :: G_dyT, G_dxT, G_uarear
    logical(kind=log_kind), allocatable, dimension(:,:) :: G_tmask

    integer(kind=int_kind) :: ierr

    character(len=*), parameter :: subname = '(dyn_evp1d_init)'

    nx=nx_global+2*nghost
    ny=ny_global+2*nghost

    allocate(G_dyT(nx,ny),G_dxT(nx,ny),G_uarear(nx,ny),G_tmask(nx,ny),stat=ierr)
    if (ierr/=0) then
       call abort_ice(subname//' ERROR: allocating', file=__FILE__, line=__LINE__)
    endif

    ! gather from blks to global
    call gather_static(G_uarear, G_dxT, G_dyT, G_tmask)

    ! calculate number of water points (T and U). Only needed for the static version
    ! tmask in ocean/ice
    if (my_task == master_task) then
      call calc_nActiveTU(G_tmask,nActive)
      call evp1d_alloc_static_na(nActive)
      call calc_2d_indices_init(nActive, G_tmask)
      call calc_navel(nActive, navel)
      call evp1d_alloc_static_navel(navel)
      call numainit(1,nActive,navel+1)
      call convert_2d_1d_init(nActive,G_HTE, G_HTN, G_uarear, G_dxT, G_dyT)
      call evp1d_alloc_static_halo()
      ! The halo geometry -- indxTij, na0, navel and the boundary types -- does
      ! not change, so the list is built once here.
      call halo_sweep(.false.)
      ! ijslot has served its purpose: slots are assigned and both halo lists
      ! are built.
      deallocate(ijslot, stat=ierr)
      if (ierr/=0) then
         call abort_ice(subname//' ERROR: deallocating', file=__FILE__, line=__LINE__)
      endif
    endif

    deallocate(G_dyT,G_dxT,G_uarear,G_tmask,stat=ierr)
    if (ierr/=0) then
       call abort_ice(subname//' ERROR: deallocating', file=__FILE__, line=__LINE__)
    endif

  end subroutine dyn_evp1d_init

!=============================================================================

  subroutine dyn_evp1d_run(L_stressp_1 , L_stressp_2 , L_stressp_3 , L_stressp_4 , &
                           L_stressm_1 , L_stressm_2 , L_stressm_3 , L_stressm_4 , &
                           L_stress12_1, L_stress12_2, L_stress12_3, L_stress12_4, &
                           L_strength,                                             &
                           L_cdn_ocn   , L_aiu       , L_uocn      , L_vocn      , &
                           L_waterxU   , L_wateryU   , L_forcexU   , L_forceyU   , &
                           L_umassdti  , L_fmU       , L_strintxU  , L_strintyU  , &
                           L_Tbu       , L_taubxU    , L_taubyU    , L_uvel      , &
                           L_vvel      , L_icetmask  , L_iceUmask)

    use ice_dyn_shared, only : ndte
    use ice_dyn_core1d, only : stress_1d, stepu_1d, calc_diag_1d
    use ice_timers    , only : ice_timer_start, ice_timer_stop, timer_evp1dcore

    use icepack_intfc , only : icepack_query_parameters, icepack_warnings_flush, &
      icepack_warnings_aborted

    implicit none

    ! nx_block, ny_block, max_blocks
    real(kind=dbl_kind)   , dimension(:,:,:), intent(inout) :: &
      L_stressp_1 , L_stressp_2 , L_stressp_3 , L_stressp_4 ,  &
      L_stressm_1 , L_stressm_2 , L_stressm_3 , L_stressm_4 ,  &
      L_stress12_1, L_stress12_2, L_stress12_3, L_stress12_4,  &
      L_strintxU  , L_strintyU  , L_uvel      , L_vvel      ,  &
      L_taubxU    , L_taubyU
    real(kind=dbl_kind)   , dimension(:,:,:), intent(in) ::    &
      L_strength  ,                                            &
      L_cdn_ocn   , L_aiu       , L_uocn     , L_vocn   ,      &
      L_waterxU   , L_wateryU   , L_forcexU  , L_forceyU,      &
      L_umassdti  , L_fmU       , L_Tbu
    logical(kind=log_kind), dimension(:,:,:), intent(in) ::    &
      L_iceUmask  , L_iceTmask

    ! local variables

    ! nx, ny
    real(kind=dbl_kind),    dimension(nx,ny) ::                &
      G_stressp_1 , G_stressp_2 , G_stressp_3 , G_stressp_4 ,  &
      G_stressm_1 , G_stressm_2 , G_stressm_3 , G_stressm_4 ,  &
      G_stress12_1, G_stress12_2, G_stress12_3, G_stress12_4,  &
      G_strength,                                              &
      G_cdn_ocn   , G_aiu       , G_uocn      , G_vocn      ,  &
      G_waterxU   , G_wateryU   , G_forcexU   , G_forceyU   ,  &
      G_umassdti  , G_fmU       , G_strintxU  , G_strintyU  ,  &
      G_Tbu       , G_uvel     , G_vvel       , G_taubxU    ,  &
      G_taubyU                  ! G_taubxU and G_taubyU are post processed from Cb
    logical(kind=log_kind), dimension (nx,ny) ::                 &
      G_iceUmask  , G_iceTmask

    character(len=*), parameter :: subname = '(dyn_evp1d_run)'

    integer(kind=int_kind) :: ksub

    real   (kind=dbl_kind) :: rhow

    ! From 3d to 2d on master task
    call gather_dyn(L_stressp_1 , L_stressp_2 , L_stressp_3 , L_stressp_4 , &
                    L_stressm_1 , L_stressm_2 , L_stressm_3 , L_stressm_4 , &
                    L_stress12_1, L_stress12_2, L_stress12_3, L_stress12_4, &
                    L_strength,                                             &
                    L_cdn_ocn   , L_aiu       , L_uocn      , L_vocn      , &
                    L_waterxU   , L_wateryU   , L_forcexU   , L_forceyU   , &
                    L_umassdti  , L_fmU       ,                             &
                    L_Tbu       , L_uvel      , L_vvel      ,               &
                    L_icetmask  , L_iceUmask  ,                             &
                    L_strintxU  , L_strintyU  , L_taubxU    , L_taubyU    , &
                    G_strintxU  , G_strintyU  , G_taubxU    , G_taubyU    , &
                    G_stressp_1 , G_stressp_2 , G_stressp_3 , G_stressp_4 , &
                    G_stressm_1 , G_stressm_2 , G_stressm_3 , G_stressm_4 , &
                    G_stress12_1, G_stress12_2, G_stress12_3, G_stress12_4, &
                    G_strength  ,                                           &
                    G_cdn_ocn   , G_aiu       , G_uocn      , G_vocn      , &
                    G_waterxU   , G_wateryU   , G_forcexU   , G_forceyU   , &
                    G_umassdti  , G_fmU       ,                             &
                    G_Tbu       , G_uvel      , G_vvel      ,               &
                    G_iceTmask,  G_iceUmask)

    if (my_task == master_task) then
       call set_skipMe(G_iceTmask, G_iceUmask,nActive)
       ! Map from 2d to 1d
       call convert_2d_1d_dyn(nActive,                                                    &
                              G_stressp_1 , G_stressp_2 , G_stressp_3 ,  G_stressp_4,     &
                              G_stressm_1 , G_stressm_2 , G_stressm_3 ,  G_stressm_4,     &
                              G_stress12_1, G_stress12_2, G_stress12_3, G_stress12_4,     &
                              G_strength,                                                 &
                              G_cdn_ocn   , G_aiu       , G_uocn     ,  G_vocn     ,      &
                              G_waterxU   , G_wateryU   , G_forcexU  , G_forceyU   ,      &
                              G_umassdti  , G_fmU       ,                                 &
                              G_Tbu       , G_uvel     , G_vvel)


       ! map from cpu to gpu (to) and back.
       ! This could be optimized considering which variables change from time step to time step
       ! and which are constant.
       ! in addition initialization of Cb and str1, str2, str3, str4, str5, str6, str7, str8
       call icepack_query_parameters(rhow_out=rhow)
       call icepack_warnings_flush(nu_diag)
       if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__, line=__LINE__)

       call ice_timer_start(timer_evp1dcore)
#ifdef _OPENMP_TARGET
       !$omp target data map(to: ee, ne, se, nw, sw, sse, skipUcell, skipTcell,&
       !$omp                 strength, dxT_1d, dyT_1d, HTE_1d,HTN_1d,HTEm1_1d, &
       !$omp                 HTNm1_1d,forcexU, forceyU, umassdti, fmU,         &
       !$omp                 uarear_1d,uvel_init, vvel_init, Tbu, Cb,          &
       !$omp                 str1, str2, str3, str4, str5, str6, str7, str8,   &
       !$omp                 cdn_ocn, aiu, uocn, vocn, waterxU, wateryU, rhow  &
       !$omp             map(tofrom: uvel,vvel,                                &
       !$omp                 stressp_1, stressp_2, stressp_3, stressp_4,       &
       !$omp                 stressm_1, stressm_2, stressm_3, stressm_4,       &
       !$omp                 stress12_1,stress12_2,stress12_3,stress12_4)
       !$omp target update to(arlx1i,denom1,capping,deltaminEVP,e_factor,epp2i,brlx)
#endif
       ! initialization of str? in order to avoid influence from old time steps
       str1(1:navel+1)=c0
       str2(1:navel+1)=c0
       str3(1:navel+1)=c0
       str4(1:navel+1)=c0
       str5(1:navel+1)=c0
       str6(1:navel+1)=c0
       str7(1:navel+1)=c0
       str8(1:navel+1)=c0

       do ksub = 1,ndte        ! subcycling
          call stress_1d (ee, ne, se, 1, nActive,                                    &
                          uvel, vvel, dxT_1d, dyT_1d, skipTcell, strength,           &
                          HTE_1d, HTN_1d, HTEm1_1d, HTNm1_1d,                        &
                          stressp_1,  stressp_2,  stressp_3,  stressp_4,             &
                          stressm_1,  stressm_2,  stressm_3,  stressm_4,             &
                          stress12_1, stress12_2, stress12_3, stress12_4,            &
                          str1, str2, str3, str4, str5, str6, str7, str8)

          call stepu_1d  (1, nActive, cdn_ocn, aiu, uocn, vocn,                         &
                          waterxU, wateryU, forcexU, forceyU, umassdti, fmU, uarear_1d, &
                          uvel_init, vvel_init, uvel, vvel,                             &
                          str1, str2, str3, str4, str5, str6, str7, str8,               &
                          nw, sw, sse, skipUcell, Tbu, Cb, rhow)
          call evp1d_halo_update()
       enddo
       ! This can be skipped if diagnostics of strintx and strinty is not needed
       ! They will either both be calculated or not.
       call calc_diag_1d(1        , nActive  , &
                         uarear_1d, skipUcell, &
                         str1     , str2     , &
                         str3     , str4     , &
                         str5     , str6     , &
                         str7     , str8     , &
                         nw       , sw       , &
                         sse      ,            &
                         strintxU, strintyU)

       call ice_timer_stop(timer_evp1dcore)

#ifdef _OPENMP_TARGET
       !$omp end target data
#endif
       ! Map results back to 2d
       call convert_1d_2d_dyn(nActive, navel,                                         &
                              G_stressp_1 , G_stressp_2 , G_stressp_3 , G_stressp_4 , &
                              G_stressm_1 , G_stressm_2 , G_stressm_3 , G_stressm_4 , &
                              G_stress12_1, G_stress12_2, G_stress12_3, G_stress12_4, &
                              G_strength  ,                                           &
                              G_cdn_ocn   , G_aiu       , G_uocn      , G_vocn      , &
                              G_waterxU   , G_wateryU   , G_forcexU   , G_forceyU   , &
                              G_umassdti  , G_fmU       , G_strintxU  , G_strintyU  , &
                              G_Tbu       , G_uvel      , G_vvel      , G_taubxU    , &
                              G_taubyU)

    endif ! master_task

    call scatter_dyn(L_stressp_1 , L_stressp_2 , L_stressp_3 , L_stressp_4 , &
                     L_stressm_1 , L_stressm_2 , L_stressm_3 , L_stressm_4 , &
                     L_stress12_1, L_stress12_2, L_stress12_3, L_stress12_4, &
                     L_strintxU  , L_strintyU  ,  L_uvel     , L_vvel      , &
                     L_taubxU    , L_taubyU    ,                             &
                     G_stressp_1 , G_stressp_2 , G_stressp_3 , G_stressp_4 , &
                     G_stressm_1 , G_stressm_2 , G_stressm_3 , G_stressm_4 , &
                     G_stress12_1, G_stress12_2, G_stress12_3, G_stress12_4, &
                     G_strintxU  , G_strintyU  , G_uvel      , G_vvel      , &
                     G_taubxU    , G_taubyU)
    ! calculate number of active points. allocate if initial or if array size should increase
    ! call calc_nActiveTU(iceTmask_log,nActive, iceUmask)
    ! if (nActiveold ==0) then ! first
    !     call evp_1d_alloc(nActive, nActive,nx,ny)
    !     nactiveold=nActive+buf1d ! allocate
    !     call init_unionTU(nx, ny, iceTmask_log,iceUmask)
    ! else if (nactiveold < nActive) then
    !     write(nu_diag,*) 'Warning nActive is bigger than old allocation. Need to re allocate'
    !     call evp_1d_dealloc() ! only deallocate if not first time step
    !     call evp_1d_alloc(nActive, nActive,nx,ny)
    !     nactiveold=nActive+buf1d ! allocate
    !     call init_unionTU(nx, ny, iceTmask_log,iceUmask)
    ! endif
    ! call cp_2dto1d(nActive)
    ! FIXME THIS IS THE LOGIC FOR RE ALLOCATION IF NEEDED
    ! call add_1d(nx, ny, natmp, iceTmask_log, iceUmask, ts)

  end subroutine dyn_evp1d_run

!=============================================================================

  subroutine dyn_evp1d_finalize()
    implicit none

    character(len=*), parameter :: subname = '(dyn_evp1d_finalize)'

    if (my_task == master_task) then
       write(nu_diag,*) 'Close evp 1d log'
    endif

  end subroutine dyn_evp1d_finalize

!=============================================================================

  subroutine evp1d_alloc_static_na(na0)
    implicit none

    integer(kind=int_kind), intent(in) :: na0
    integer(kind=int_kind) :: ierr
    character(len=*), parameter :: subname = '(evp1d_alloc_static_na)'

    allocate(skipTcell(1:na0),  &
             skipUcell(1:na0),  &
             iwidx(1:nx,1:ny),  &
             stat=ierr)

    if (ierr/=0) then
       call abort_ice(subname//' ERROR: allocating', file=__FILE__, line=__LINE__)
    endif

    allocate(indxTi(1:na0), &
             indxTj(1:na0), &
             stat=ierr)

    if (ierr/=0) then
       call abort_ice(subname//' ERROR: allocating', file=__FILE__, line=__LINE__)
    endif

    allocate(ee(1:na0) , &
             ne(1:na0) , &
             se(1:na0) , &
             nw(1:na0) , &
             sw(1:na0) , &
             sse(1:na0), &
             stat=ierr)

    if (ierr/=0) then
       call abort_ice(subname//' ERROR: allocating', file=__FILE__, line=__LINE__)
    endif

    allocate( HTE_1d    (1:na0), &
              HTN_1d    (1:na0), &
              HTEm1_1d  (1:na0), &
              HTNm1_1d  (1:na0), &
              dxT_1d    (1:na0), &
              dyT_1d    (1:na0), &
              strength  (1:na0), &
              stressp_1 (1:na0), &
              stressp_2 (1:na0), &
              stressp_3 (1:na0), &
              stressp_4 (1:na0), &
              stressm_1 (1:na0), &
              stressm_2 (1:na0), &
              stressm_3 (1:na0), &
              stressm_4 (1:na0), &
              stress12_1(1:na0), &
              stress12_2(1:na0), &
              stress12_3(1:na0), &
              stress12_4(1:na0), &
              stat=ierr)

    if (ierr/=0) then
       call abort_ice(subname//' ERROR: allocating', file=__FILE__, line=__LINE__)
    endif

    allocate(cdn_ocn  (1:na0), aiu      (1:na0), &
             uocn     (1:na0), vocn     (1:na0), &
             waterxU  (1:na0), wateryU  (1:na0), &
             forcexU  (1:na0), forceyU  (1:na0), &
             umassdti (1:na0), fmU      (1:na0), &
             uarear_1d(1:na0),                   &
             strintxU (1:na0), strintyU (1:na0), &
             Tbu      (1:na0), Cb       (1:na0), &
             uvel_init(1:na0), vvel_init(1:na0), &
             stat=ierr)

    if (ierr/=0) then
       call abort_ice(subname//' ERROR: allocating', file=__FILE__, line=__LINE__)
    endif

  end subroutine evp1d_alloc_static_na

!=============================================================================

  subroutine evp1d_alloc_static_navel(navel0)
    implicit none

    integer(kind=int_kind), intent(in) :: navel0
    integer(kind=int_kind) :: ierr
    character(len=*), parameter :: subname = '(evp1d_alloc_static_na)'

    ! navel0+1: the last slot is the dead slot.  Directions that leave the
    ! grid point at it, so the kernels can read a neighbour unconditionally
    ! and get zero -- the same convention the compressed 2-D layout uses.
    allocate(str1(1:navel0+1)   , str2(1:navel0+1), str3(1:navel0+1), &
             str4(1:navel0+1)   , str5(1:navel0+1), str6(1:navel0+1), &
             str7(1:navel0+1)   , str8(1:navel0+1),                   &
             indxTij(1:navel0+1), uvel(1:navel0+1), vvel(1:navel0+1), &
             stat=ierr)

    if (ierr/=0) then
       call abort_ice(subname//' ERROR: allocating', file=__FILE__, line=__LINE__)
    endif

  end subroutine evp1d_alloc_static_navel

!=============================================================================

  subroutine evp1d_alloc_static_halo()

    implicit none
    integer(kind=int_kind) :: ierr
    character(len=*), parameter :: subname = '(evp1d_alloc_static_halo)'

    ! One entry per ghost cell per edge, for the cyclic pass and again for
    ! the zero_gradient/linear_extrap pass.
    allocate(halo_bc_dst(4*nghost*(nx+ny)), halo_bc_s1(4*nghost*(nx+ny)), &
             halo_bc_s2(4*nghost*(nx+ny)), halo_bc_w (4*nghost*(nx+ny)), &
             stat=ierr)

    if (ierr/=0) then
       call abort_ice(subname//' ERROR: allocating', file=__FILE__, line=__LINE__)
    endif

  end subroutine evp1d_alloc_static_halo

!=============================================================================

  subroutine calc_nActiveTU(Tmask,na0, Umask)

    ! Calculate number of active points with a given mask.

    implicit none
    logical(kind=log_kind), intent(in) :: Tmask(:,:)
    logical(kind=log_kind), optional, intent(in) :: Umask(:,:)
    integer(kind=int_kind), intent(out)  :: na0
    integer(kind=int_kind)              :: i,j
    character(len=*), parameter :: subname = '(calc_nActivceTU)'

    na0=0
    if (present(Umask)) then
       do i=1+nghost,nx
       do j=1+nghost,ny
          if ((Tmask(i,j)) .or. (Umask(i,j))) then
             na0=na0+1
          endif
       enddo
       enddo
    else
       do i=1+nghost,nx
       do j=1+nghost,ny
          if (Tmask(i,j)) then
             na0=na0+1
          endif
       enddo
       enddo
    endif

  end subroutine calc_nActiveTU

!=============================================================================

  subroutine set_skipMe(iceTmask, iceUmask,na0)

    implicit none

    logical(kind=log_kind), intent(in) :: iceTmask(:,:), iceUmask(:,:)
    integer(kind=int_kind), intent(in) :: na0
    integer(kind=int_kind)              :: iw, i, j, niw
    character(len=*), parameter :: subname = '(set_skipMe)'

    skipUcell=.false.
    skipTcell=.false.
    niw=0
    ! first count
    do iw=1, na0
      i = indxti(iw)
      j = indxtj(iw)
      if ( iceTmask(i,j) .or. iceUmask(i,j)) then
         niw=niw+1
      endif
      if (.not. (iceTmask(i,j))) skipTcell(iw)=.true.
      if (.not. (iceUmask(i,j))) skipUcell(iw)=.true.
      if (i == nx)  skipUcell(iw)=.true.
      if (j == ny)  skipUcell(iw)=.true.
    enddo
    !    write(nu_diag,*) 'number of points and Active points', na0, niw

  end subroutine set_skipMe

!=============================================================================

  subroutine calc_2d_indices_init(na0, Tmask)
    ! All points are active. Need to find neighbors.
    ! This should include de selection of u points.

    implicit none

    integer(kind=int_kind), intent(in) :: na0
    ! nx, ny
    logical(kind=log_kind), dimension(:,:), intent(in) :: Tmask

    ! local variables

    integer(kind=int_kind) :: i, j, Nmaskt
    character(len=*), parameter :: subname = '(calc_2d_indices_init)'

    indxti(:) = 0
    indxtj(:) = 0
    Nmaskt = 0
    ! NOTE: T mask includes northern and eastern ghost cells
    do j = 1 + nghost, ny
    do i = 1 + nghost, nx
       if (Tmask(i,j)) then
          Nmaskt = Nmaskt + 1
          indxti(Nmaskt) = i
          indxtj(Nmaskt) = j
       end if
    end do
    end do

  end subroutine calc_2d_indices_init

  subroutine gather_static(G_uarear, G_dxT, G_dyT, G_Tmask)

     ! In standalone  distrb_info is an integer. Not needed anyway
     use ice_communicate, only : master_task
     use ice_gather_scatter, only : gather_global
     use ice_domain, only : distrb_info
     use ice_grid, only: dyT, dxT, uarear, tmask
     implicit none

     ! nx, ny
     real(kind=dbl_kind)   , dimension(:,:), intent(out) :: G_uarear, G_dxT, G_dyT
     logical(kind=log_kind), dimension(:,:), intent(out) :: G_Tmask

     character(len=*), parameter :: subname = '(gather_static)'

     G_uarear = c0
     G_dyT = c0
     G_dxT = c0
     G_tmask = .false.

     ! copy from distributed I_* to G_*
     call gather_global(G_uarear, uarear, master_task, distrb_info, grid_ext=.true.)
     call gather_global(G_dxT   , dxT   , master_task, distrb_info, grid_ext=.true.)
     call gather_global(G_dyT   , dyT   , master_task, distrb_info, grid_ext=.true.)
     call gather_global(G_Tmask , Tmask , master_task, distrb_info, grid_ext=.true.)

  end subroutine gather_static

!=============================================================================

  subroutine gather_dyn(L_stressp_1 , L_stressp_2 , L_stressp_3 , L_stressp_4 , &
                        L_stressm_1 , L_stressm_2 , L_stressm_3 , L_stressm_4 , &
                        L_stress12_1, L_stress12_2, L_stress12_3,L_stress12_4 , &
                        L_strength  ,                                           &
                        L_cdn_ocn   , L_aiu       , L_uocn      , L_vocn      , &
                        L_waterxU   , L_wateryU   , L_forcexU   , L_forceyU   , &
                        L_umassdti  , L_fmU       ,                             &
                        L_Tbu       , L_uvel      , L_vvel      ,               &
                        L_icetmask  , L_iceUmask  ,                             &
                        L_strintxU  , L_strintyU  , L_taubxU    , L_taubyU    , &
                        G_strintxU  , G_strintyU  , G_taubxU    , G_taubyU    , &
                        G_stressp_1 , G_stressp_2 , G_stressp_3 , G_stressp_4 , &
                        G_stressm_1 , G_stressm_2 , G_stressm_3 , G_stressm_4 , &
                        G_stress12_1, G_stress12_2, G_stress12_3, G_stress12_4, &
                        G_strength,                                             &
                        G_cdn_ocn   , G_aiu       , G_uocn      , G_vocn      , &
                        G_waterxU   , G_wateryU   , G_forcexU   , G_forceyU   , &
                        G_umassdti  , G_fmU       ,                             &
                        G_Tbu       , G_uvel      , G_vvel      ,               &
                        G_iceTmask,  G_iceUmask)

     use ice_communicate, only : master_task
     use ice_gather_scatter, only : gather_global
     use ice_domain, only : distrb_info
     implicit none

     ! nx_block, ny_block, max_blocks
     real(kind=dbl_kind)   , dimension(:,:,:), intent(in)  ::   &
        L_stressp_1 , L_stressp_2 , L_stressp_3 , L_stressp_4 , &
        L_stressm_1 , L_stressm_2 , L_stressm_3 , L_stressm_4 , &
        L_stress12_1, L_stress12_2, L_stress12_3, L_stress12_4, &
        L_strength  ,                                           &
        L_cdn_ocn   , L_aiu       , L_uocn      , L_vocn      , &
        L_waterxU   , L_wateryU   , L_forcexU   , L_forceyU   , &
        L_umassdti  , L_fmU       ,                             &
        L_Tbu       , L_uvel      , L_vvel
     logical(kind=log_kind), dimension(:,:,:), intent(in)  ::   &
        L_iceUmask  , L_iceTmask
     ! the four diagnostics, gathered so the round trip leaves alone any cell
     ! this solver does not own -- see the note in convert_1d_2d_dyn
     real(kind=dbl_kind)   , dimension(:,:,:), intent(in)  ::   &
        L_strintxU  , L_strintyU  , L_taubxU    , L_taubyU
     real(kind=dbl_kind)   , dimension(:,:), intent(inout) ::   &
        G_strintxU  , G_strintyU  , G_taubxU    , G_taubyU

     ! nx, ny
     real(kind=dbl_kind)   , dimension(:,:), intent(out) ::     &
        G_stressp_1 , G_stressp_2 , G_stressp_3 , G_stressp_4 , &
        G_stressm_1 , G_stressm_2 , G_stressm_3 , G_stressm_4 , &
        G_stress12_1, G_stress12_2, G_stress12_3, G_stress12_4, &
        G_strength,                                             &
        G_cdn_ocn   , G_aiu       , G_uocn      , G_vocn      , &
        G_waterxU   , G_wateryU   , G_forcexU   , G_forceyU   , &
        G_umassdti  , G_fmU       ,                             &
        G_Tbu       , G_uvel      , G_vvel
     logical(kind=log_kind), dimension(:,:), intent(out) ::     &
        G_iceUmask  , G_iceTmask

     character(len=*), parameter :: subname = '(gather_dyn)'

     ! copy from distributed I_* to G_*
     call gather_global(G_stressp_1 ,     L_stressp_1,     master_task, distrb_info,c0, grid_ext=.true.)
     call gather_global(G_stressp_2 ,     L_stressp_2,     master_task, distrb_info,c0, grid_ext=.true.)
     call gather_global(G_stressp_3 ,     L_stressp_3,     master_task, distrb_info,c0, grid_ext=.true.)
     call gather_global(G_stressp_4 ,     L_stressp_4,     master_task, distrb_info,c0, grid_ext=.true.)

     call gather_global(G_stressm_1 ,     L_stressm_1,     master_task, distrb_info,c0, grid_ext=.true.)
     call gather_global(G_stressm_2 ,     L_stressm_2,     master_task, distrb_info,c0, grid_ext=.true.)
     call gather_global(G_stressm_3 ,     L_stressm_3,     master_task, distrb_info,c0, grid_ext=.true.)
     call gather_global(G_stressm_4 ,     L_stressm_4,     master_task, distrb_info,c0, grid_ext=.true.)

     call gather_global(G_stress12_1,     L_stress12_1,    master_task, distrb_info,c0, grid_ext=.true.)
     call gather_global(G_stress12_2,     L_stress12_2,    master_task, distrb_info,c0, grid_ext=.true.)
     call gather_global(G_stress12_3,     L_stress12_3,    master_task, distrb_info,c0, grid_ext=.true.)
     call gather_global(G_stress12_4,     L_stress12_4,    master_task, distrb_info,c0, grid_ext=.true.)
     call gather_global(G_strength  ,     L_strength  ,    master_task, distrb_info,c0, grid_ext=.true.)

     call gather_global(G_cdn_ocn   ,     L_cdn_ocn   ,     master_task, distrb_info, grid_ext=.true.)
     call gather_global(G_aiu       ,     L_aiu       ,     master_task, distrb_info, grid_ext=.true.)
     call gather_global(G_uocn      ,     L_uocn      ,     master_task, distrb_info, grid_ext=.true.)
     call gather_global(G_vocn      ,     L_vocn      ,     master_task, distrb_info, grid_ext=.true.)

     call gather_global(G_waterxU   ,     L_waterxU   ,     master_task, distrb_info, grid_ext=.true.)
     call gather_global(G_wateryU   ,     L_wateryU   ,     master_task, distrb_info, grid_ext=.true.)
     call gather_global(G_forcexU   ,     L_forcexU   ,     master_task, distrb_info, grid_ext=.true.)
     call gather_global(G_forceyU   ,     L_forceyU   ,     master_task, distrb_info, grid_ext=.true.)

     call gather_global(G_umassdti  ,     L_umassdti  ,     master_task, distrb_info, grid_ext=.true.)
     call gather_global(G_fmU       ,     L_fmU       ,     master_task, distrb_info, grid_ext=.true.)

     call gather_global(G_Tbu       ,     L_Tbu       ,     master_task, distrb_info, grid_ext=.true.)
     call gather_global(G_uvel      ,     L_uvel      ,     master_task, distrb_info,c0, grid_ext=.true.)
     call gather_global(G_vvel      ,     L_vvel      ,     master_task, distrb_info,c0, grid_ext=.true.)
     ! Symmetric with the scatter in scatter_dyn, so the round trip is an
     ! identity at every cell this solver does not write.
     call gather_global(G_strintxU  ,     L_strintxU  ,     master_task, distrb_info,c0, grid_ext=.true.)
     call gather_global(G_strintyU  ,     L_strintyU  ,     master_task, distrb_info,c0, grid_ext=.true.)
     call gather_global(G_taubxU    ,     L_taubxU    ,     master_task, distrb_info,c0, grid_ext=.true.)
     call gather_global(G_taubyU    ,     L_taubyU    ,     master_task, distrb_info,c0, grid_ext=.true.)
     call gather_global(G_iceTmask  ,     L_iceTmask  ,     master_task, distrb_info, grid_ext=.true.)
     call gather_global(G_iceUmask  ,     L_iceUmask  ,     master_task, distrb_info, grid_ext=.true.)

  end subroutine gather_dyn

!=============================================================================

  subroutine scatter_dyn(L_stressp_1 , L_stressp_2 , L_stressp_3 , L_stressp_4 , &
                         L_stressm_1 , L_stressm_2 , L_stressm_3 , L_stressm_4 , &
                         L_stress12_1, L_stress12_2, L_stress12_3, L_stress12_4, &
                         L_strintxU  , L_strintyU  , L_uvel      , L_vvel      , &
                         L_taubxU    , L_taubyU    ,                             &
                         G_stressp_1 , G_stressp_2 , G_stressp_3 ,  G_stressp_4, &
                         G_stressm_1 , G_stressm_2 , G_stressm_3 ,  G_stressm_4, &
                         G_stress12_1, G_stress12_2, G_stress12_3, G_stress12_4, &
                         G_strintxU  , G_strintyU  , G_uvel      , G_vvel      , &
                         G_taubxU    , G_taubyU )

     use ice_communicate, only : master_task
     use ice_gather_scatter, only : scatter_global
     use ice_domain, only : distrb_info
     implicit none

     ! nx_block, ny_block, max_blocks
     real(kind=dbl_kind), dimension(:,:,:), intent(out) :: &
        L_stressp_1 , L_stressp_2 , L_stressp_3 , L_stressp_4 , &
        L_stressm_1 , L_stressm_2 , L_stressm_3 , L_stressm_4 , &
        L_stress12_1, L_stress12_2, L_stress12_3, L_stress12_4, &
        L_strintxU  , L_strintyU  , L_uvel      , L_vvel      , &
        L_taubxU    , L_taubyU

     ! nx, ny
     real(kind=dbl_kind), dimension(:,:), intent(in) ::         &
        G_stressp_1 , G_stressp_2 , G_stressp_3 , G_stressp_4 , &
        G_stressm_1 , G_stressm_2 , G_stressm_3 , G_stressm_4 , &
        G_stress12_1, G_stress12_2, G_stress12_3, G_stress12_4, &
        G_strintxU  , G_strintyU  , G_uvel      , G_vvel      , &
        G_taubxU    , G_taubyU

     character(len=*), parameter :: subname = '(scatter_dyn)'

     call scatter_global(L_stressp_1,  G_stressp_1,  master_task, distrb_info, grid_ext=.true.)
     call scatter_global(L_stressp_2,  G_stressp_2,  master_task, distrb_info, grid_ext=.true.)
     call scatter_global(L_stressp_3,  G_stressp_3,  master_task, distrb_info, grid_ext=.true.)
     call scatter_global(L_stressp_4,  G_stressp_4,  master_task, distrb_info, grid_ext=.true.)

     call scatter_global(L_stressm_1,  G_stressm_1,  master_task, distrb_info, grid_ext=.true.)
     call scatter_global(L_stressm_2,  G_stressm_2,  master_task, distrb_info, grid_ext=.true.)
     call scatter_global(L_stressm_3,  G_stressm_3,  master_task, distrb_info, grid_ext=.true.)
     call scatter_global(L_stressm_4,  G_stressm_4,  master_task, distrb_info, grid_ext=.true.)

     call scatter_global(L_stress12_1, G_stress12_1, master_task, distrb_info, grid_ext=.true.)
     call scatter_global(L_stress12_2, G_stress12_2, master_task, distrb_info, grid_ext=.true.)
     call scatter_global(L_stress12_3, G_stress12_3, master_task, distrb_info, grid_ext=.true.)
     call scatter_global(L_stress12_4, G_stress12_4, master_task, distrb_info, grid_ext=.true.)

     call scatter_global(L_strintxU  , G_strintxU  , master_task, distrb_info, grid_ext=.true.)
     call scatter_global(L_strintyU  , G_strintyU  , master_task, distrb_info, grid_ext=.true.)
     call scatter_global(L_uvel      , G_uvel      , master_task, distrb_info, grid_ext=.true.)
     call scatter_global(L_vvel      , G_vvel      , master_task, distrb_info, grid_ext=.true.)
     call scatter_global(L_taubxU    , G_taubxU    , master_task, distrb_info, grid_ext=.true.)
     call scatter_global(L_taubyU    , G_taubyU    , master_task, distrb_info, grid_ext=.true.)

  end subroutine scatter_dyn

!=============================================================================

  subroutine convert_2d_1d_init(na0, G_HTE, G_HTN, G_uarear,  G_dxT, G_dyT)

     implicit none

     integer(kind=int_kind), intent(in) ::  na0
     real (kind=dbl_kind), dimension(:, :), intent(in) :: G_HTE, G_HTN, G_uarear, G_dxT, G_dyT

     ! local variables

     integer(kind=int_kind) :: iw, lo, up, j, i, ij, islot, ierr

     character(len=*), parameter :: subname = '(convert_2d_1d_init)'

     ! Hand out the slots.  The layout is the one the kernels rely on: the
     ! active T and U points first, 1..na0, in the order calc_2d_indices_init
     ! found them, then the cells that are only ever read.  stress_1d and
     ! stepu_1d both run over 1..na0, so nothing downstream has to know where
     ! the halo tail begins beyond na0 itself.
     do iw = 1, na0
        ij = indxti(iw) + (indxtj(iw) - 1) * nx
        ijslot(ij)  = iw
        indxTij(iw) = ij
     end do

     ! The read-only cells take the slots after them, in ascending linear
     ! index -- the same order the sorted union used to produce.
     islot = na0
     do ij = 1, nx * ny
        if (ijslot(ij) == -1) then
           islot = islot + 1
           ijslot(ij)     = islot
           indxTij(islot) = ij
        endif
     end do
     if (islot /= navel) then
        call abort_ice(subname//' ERROR: slot count does not match navel', &
                       file=__FILE__, line=__LINE__)
     endif

     ! Neighbour slots.  A direction that leaves the grid takes the dead slot
     ! at navel+1, which holds zero and is never written.  That only happens
     ! for nw/sw/sse at i == nx or j == ny, where the cell has no U point to
     ! step -- skipUcell is already .true. there, matching the 2-D solver,
     ! which computes stress over jlo..jhi+1 but steps U only over jlo..jhi.
     ! ee/ne/se step towards i-1 and j-1 and active T cells start at
     ! 1+nghost, so for nghost >= 1 they always land inside the grid.
     do iw = 1, na0
        i = indxti(iw)
        j = indxtj(iw)
        ee (iw) = slot_of(i - 1, j    , navel + 1)
        ne (iw) = slot_of(i - 1, j - 1, navel + 1)
        se (iw) = slot_of(i    , j - 1, navel + 1)
        nw (iw) = slot_of(i + 1, j    , navel + 1)
        sw (iw) = slot_of(i + 1, j + 1, navel + 1)
        sse(iw) = slot_of(i    , j + 1, navel + 1)
     end do

     !tar      i$OMP PARALLEL PRIVATE(iw, lo, up, j, i)
     ! write 1D arrays from 2D arrays (target points)
     !tar      call domp_get_domain(1, na0, lo, up)
     lo=1
     up=na0
     do iw = 1, na0
        ! get 2D indices
        i = indxti(iw)
        j = indxtj(iw)
        ! map
        uarear_1d(iw)  = G_uarear(i, j)
        dxT_1d(iw)     = G_dxT(i, j)
        dyT_1d(iw)     = G_dyT(i, j)
        HTE_1d(iw)     = G_HTE(i, j)
        HTN_1d(iw)     = G_HTN(i, j)
        HTEm1_1d(iw)   = G_HTE(i - 1, j)
        HTNm1_1d(iw)   = G_HTN(i, j - 1)
     end do

  end subroutine convert_2d_1d_init

!=============================================================================

  subroutine convert_2d_1d_dyn(na0         ,                                           &
                               G_stressp_1 , G_stressp_2 , G_stressp_3 , G_stressp_4 , &
                               G_stressm_1 , G_stressm_2 , G_stressm_3 , G_stressm_4 , &
                               G_stress12_1, G_stress12_2, G_stress12_3, G_stress12_4, &
                               G_strength  , G_cdn_ocn   , G_aiu       , G_uocn      , &
                               G_vocn      , G_waterxU   , G_wateryU   , G_forcexU   , &
                               G_forceyU   , G_umassdti  , G_fmU       ,  G_Tbu      , &
                               G_uvel      , G_vvel      )

     implicit none

     integer(kind=int_kind), intent(in) ::  na0

     ! nx, ny
     real(kind=dbl_kind), dimension(:, :), intent(in) :: &
        G_stressp_1 , G_stressp_2 , G_stressp_3 , G_stressp_4, &
        G_stressm_1 , G_stressm_2 , G_stressm_3 , G_stressm_4, &
        G_stress12_1, G_stress12_2, G_stress12_3,G_stress12_4, &
        G_strength  , G_cdn_ocn   , G_aiu       , G_uocn     , &
        G_vocn      , G_waterxU   , G_wateryU   , G_forcexU  , &
        G_forceyU   , G_umassdti  , G_fmU       , G_Tbu      , &
        G_uvel      , G_vvel

     integer(kind=int_kind) ::  lo, up, iw, i, j
     character(len=*), parameter :: subname = '(convert_2d_1d_dyn)'

     lo=1
     up=na0
     do iw = 1, na0
        ! get 2D indices
        i = indxti(iw)
        j = indxtj(iw)
        ! map
        stressp_1(iw)  = G_stressp_1(i, j)
        stressp_2(iw)  = G_stressp_2(i, j)
        stressp_3(iw)  = G_stressp_3(i, j)
        stressp_4(iw)  = G_stressp_4(i, j)
        stressm_1(iw)  = G_stressm_1(i, j)
        stressm_2(iw)  = G_stressm_2(i, j)
        stressm_3(iw)  = G_stressm_3(i, j)
        stressm_4(iw)  = G_stressm_4(i, j)
        stress12_1(iw) = G_stress12_1(i, j)
        stress12_2(iw) = G_stress12_2(i, j)
        stress12_3(iw) = G_stress12_3(i, j)
        stress12_4(iw) = G_stress12_4(i, j)
        strength(iw)   = G_strength(i,j)
        cdn_ocn(iw)    = G_cdn_ocn(i, j)
        aiu(iw)        = G_aiu(i, j)
        uocn(iw)       = G_uocn(i, j)
        vocn(iw)       = G_vocn(i, j)
        waterxU(iw)    = G_waterxU(i, j)
        wateryU(iw)    = G_wateryU(i, j)
        forcexU(iw)    = G_forcexU(i, j)
        forceyU(iw)    = G_forceyU(i, j)
        umassdti(iw)   = G_umassdti(i, j)
        fmU(iw)        = G_fmU(i, j)
        strintxU(iw)   = C0
        strintyU(iw)   = C0
        Tbu(iw)        = G_Tbu(i, j)
        Cb(iw)         = c0
        uvel(iw)       = G_uvel(i,j)
        vvel(iw)       = G_vvel(i,j)
        uvel_init(iw)  = G_uvel(i,j)
        vvel_init(iw)  = G_vvel(i,j)
     end do

     ! Halos can potentially have values of u and v
     do iw=na0+1,navel
        j = int((indxTij(iw) - 1) / (nx)) + 1
        i =      indxTij(iw) - (j - 1) * nx
        uvel(iw)=G_uvel(i,j)
        vvel(iw)=G_vvel(i,j)
     end do

  end subroutine convert_2d_1d_dyn

!=============================================================================

  subroutine convert_1d_2d_dyn(na0         , navel0      ,                             &
                               G_stressp_1 , G_stressp_2 , G_stressp_3 , G_stressp_4 , &
                               G_stressm_1 , G_stressm_2 , G_stressm_3 , G_stressm_4 , &
                               G_stress12_1, G_stress12_2, G_stress12_3, G_stress12_4, &
                               G_strength,                                             &
                               G_cdn_ocn   , G_aiu       , G_uocn      , G_vocn      , &
                               G_waterxU   , G_wateryU   , G_forcexU   , G_forceyU   , &
                               G_umassdti  , G_fmU       , G_strintxU  , G_strintyU  , &
                               G_Tbu       , G_uvel      , G_vvel      , G_taubxU    , &
                               G_taubyU)

     implicit none

     integer(kind=int_kind), intent(in) ::  na0, navel0
     ! nx, ny
     real(kind=dbl_kind), dimension(:, :), intent(inout) ::     &
        G_stressp_1 , G_stressp_2 , G_stressp_3 , G_stressp_4 , &
        G_stressm_1 , G_stressm_2 , G_stressm_3 , G_stressm_4 , &
        G_stress12_1, G_stress12_2, G_stress12_3, G_stress12_4, &
        G_strength,                                             &
        G_cdn_ocn   , G_aiu       , G_uocn      , G_vocn      , &
        G_waterxU   , G_wateryU   , G_forcexU   , G_forceyU   , &
        G_umassdti  , G_fmU       , G_strintxU  , G_strintyU  , &
        G_Tbu       , G_uvel      , G_vvel      , G_taubxU    , &
        G_taubyU

     integer(kind=int_kind) ::  lo, up, iw, i, j
     character(len=*), parameter :: subname = '(convert_1d_2d_dyn)'

     G_stressp_1  = c0
     G_stressp_2  = c0
     G_stressp_3  = c0
     G_stressp_4  = c0
     G_stressm_1  = c0
     G_stressm_2  = c0
     G_stressm_3  = c0
     G_stressm_4  = c0
     G_stress12_1 = c0
     G_stress12_2 = c0
     G_stress12_3 = c0
     G_stress12_4 = c0
     G_strength   = c0
     G_cdn_ocn    = c0
     G_aiu        = c0
     G_uocn       = c0
     G_vocn       = c0
     G_waterxU    = c0
     G_wateryU    = c0
     G_forcexU    = c0
     G_forceyU    = c0
     G_umassdti   = c0
     G_fmU        = c0
     ! G_strintxU, G_strintyU, G_taubxU and G_taubyU are deliberately NOT
     ! zeroed, for the same reason as G_uvel and G_vvel below: gather_dyn has
     ! filled them with the block field, and scatter_dyn writes them straight
     ! back.  The 2-D solver's unpack_from_slots writes only the cells it has
     ! a slot for and leaves the rest, so zeroing here made evp1d clear
     ! diagnostics at cells it does not own where the 2-D solver does not.
     G_Tbu        = c0
     ! G_uvel and G_vvel are deliberately NOT zeroed.  Unlike the others they
     ! arrive holding the gathered block field -- gather_dyn fills them with
     ! grid_ext, ghosts included -- and they are scattered straight back.  A
     ! cell with no slot is one this solver does not own, typically a ghost
     ! whose interior neighbour is land; zeroing it here made the scatter
     ! overwrite the block value with zero, so evp1d silently cleared
     ! velocities that the 2-D solver leaves alone.  Leaving them means the
     ! scatter writes back what it read at every cell the solver does not
     ! touch.  Both loops below write every slot, so nothing stale survives
     ! where the solver does own the cell.

     lo=1
     up=na0
     do iw = lo, up
        ! get 2D indices
        i = indxti(iw)
        j = indxtj(iw)
        ! map to 2d
        G_stressp_1 (i,j) = stressp_1(iw)
        G_stressp_2 (i,j) = stressp_2(iw)
        G_stressp_3 (i,j) = stressp_3(iw)
        G_stressp_4 (i,j) = stressp_4(iw)
        G_stressm_1 (i,j) = stressm_1(iw)
        G_stressm_2 (i,j) = stressm_2(iw)
        G_stressm_3 (i,j) = stressm_3(iw)
        G_stressm_4 (i,j) = stressm_4(iw)
        G_stress12_1(i,j) = stress12_1(iw)
        G_stress12_2(i,j) = stress12_2(iw)
        G_stress12_3(i,j) = stress12_3(iw)
        G_stress12_4(i,j) = stress12_4(iw)
        G_strintxU(i,j)   = strintxU(iw)
        G_strintyU(i,j)   = strintyU (iw)
        G_taubxU(i,j)     = -uvel(iw)*Cb(iw)
        G_taubyU(i,j)     = -vvel(iw)*Cb(iw)
        G_uvel(i,j)       = uvel(iw)
        G_vvel(i,j)       = vvel(iw)
     end do

     do iw=na0+1,navel0
        j = int((indxTij(iw) - 1) / (nx)) + 1
        i = indxTij(iw) - (j - 1) * nx
        G_uvel(i,j)       = uvel(iw)
        G_vvel(i,j)       = vvel(iw)
     end do

  end subroutine convert_1d_2d_dyn

  subroutine calc_navel(na0, navel0)
     ! Calculate number of active points, including halo points.
     !
     ! Marks every cell the 1-D solver touches -- each active T cell, plus the
     ! six neighbours its kernels read -- and counts them.  The marks stay in
     ! ijslot for convert_2d_1d_init, which turns them into slot numbers.
     !
     ! This replaces a sorted-merge union of seven index vectors built from
     ! i + (j-1)*nx arithmetic.  That arithmetic has no way to say "there is no
     ! such neighbour": at i == nx the (+1,.) directions wrapped round into the
     ! next row, and at j == ny they ran off the end of the grid entirely, so
     ! the union carried indices naming cells that do not exist.  A scan with
     ! explicit bounds guards cannot generate them.  It is how set_neighbours
     ! builds dynNbr for the compressed 2-D layout, and the point of doing it
     ! the same way here is that the two solvers then agree by construction
     ! rather than by two separate pieces of index arithmetic agreeing.

     implicit none

     integer(kind=int_kind), intent(in) :: na0
     integer(kind=int_kind), intent(out) :: navel0

     ! local variables

     integer(kind=int_kind) :: iw, i, j, ij, ierr

     character(len=*), parameter :: subname = '(calc_navel)'

     allocate(ijslot(1:nx*ny), stat=ierr)
     if (ierr/=0) then
        call abort_ice(subname//' ERROR: allocating', file=__FILE__, line=__LINE__)
     endif

     ijslot(:) = 0

     do iw = 1, na0
        i = indxti(iw)
        j = indxtj(iw)
        ! No trailing backslashes in these comments: .F90 is preprocessed and
        ! cpp splices a line ending in one into the comment above it, which
        ! silently deletes the call.
        call mark_cell(i    , j    )   ! ( 0,  0) this cell
        ! the U points stress_1d reads
        call mark_cell(i - 1, j    )   ! (-1,  0)
        call mark_cell(i - 1, j - 1)   ! (-1, -1)
        call mark_cell(i    , j - 1)   ! ( 0, -1)
        ! the T cells stepu_1d reads
        call mark_cell(i + 1, j    )   ! (+1,  0)
        call mark_cell(i + 1, j + 1)   ! (+1, +1)
        call mark_cell(i    , j + 1)   ! ( 0, +1)
     end do

     ! the halo update reads and writes cells of its own, beyond anything the
     ! kernels reach -- see halo_sweep
     call halo_sweep(.true.)

     navel0 = 0
     do ij = 1, nx * ny
        if (ijslot(ij) /= 0) navel0 = navel0 + 1
     end do

  end subroutine calc_navel

!=======================================================================

  subroutine mark_cell(i, j)
     ! Mark (i,j) as a cell the solver needs, if it is inside the gathered
     ! grid.  Outside means the neighbour does not exist, which is precisely
     ! what the index arithmetic this replaces could not express.

     implicit none

     integer(kind=int_kind), intent(in) :: i, j

     character(len=*), parameter :: subname = '(mark_cell)'

     if (i < 1 .or. i > nx) return
     if (j < 1 .or. j > ny) return
     ijslot(i + (j - 1) * nx) = -1

  end subroutine mark_cell

!=======================================================================

  integer(kind=int_kind) function slot_of(i, j, dead)
     ! Slot holding cell (i,j), or the dead slot if it is off the grid.

     implicit none

     integer(kind=int_kind), intent(in) :: i, j, dead

     character(len=*), parameter :: subname = '(slot_of)'

     if (i < 1 .or. i > nx .or. j < 1 .or. j > ny) then
        slot_of = dead
     else
        slot_of = ijslot(i + (j - 1) * nx)
     endif

  end function slot_of

!=======================================================================

  subroutine numainit(lo,up,uu)

     implicit none
     integer(kind=int_kind),intent(in) :: lo,up,uu
     integer(kind=int_kind) :: iw
     character(len=*), parameter :: subname = '(numainit)'

     !$omp parallel do schedule(runtime) private(iw)
     do iw = lo,up
        skipTcell(iw)=.false.
        skipUcell(iw)=.false.
        ee(iw)=0
        ne(iw)=0
        se(iw)=0
        nw(iw)=0
        sw(iw)=0
        sse(iw)=0
        aiu(iw)=c0
        Cb(iw)=c0
        cdn_ocn(iw)=c0
        dxT_1d(iw)=c0
        dyT_1d(iw)=c0
        fmU(iw)=c0
        forcexU(iw)=c0
        forceyU(iw)=c0
        HTE_1d(iw)=c0
        HTEm1_1d(iw)=c0
        HTN_1d(iw)=c0
        HTNm1_1d(iw)=c0
        strength(iw)= c0
        stress12_1(iw)=c0
        stress12_2(iw)=c0
        stress12_3(iw)=c0
        stress12_4(iw)=c0
        stressm_1(iw)=c0
        stressm_2(iw)=c0
        stressm_3(iw)=c0
        stressm_4(iw)=c0
        stressp_1(iw)=c0
        stressp_2(iw)=c0
        stressp_3(iw)=c0
        stressp_4(iw)=c0
        strintxU(iw)= c0
        strintyU(iw)= c0
        Tbu(iw)=c0
        uarear_1d(iw)=c0
        umassdti(iw)=c0
        uocn(iw)=c0
        uvel_init(iw)=c0
        uvel(iw)=c0
        vocn(iw)=c0
        vvel_init(iw)=c0
        vvel(iw)=c0
        waterxU(iw)=c0
        wateryU(iw)=c0
     enddo
     !$omp end parallel do
     !$omp parallel do schedule(runtime) private(iw)
     do iw = lo,uu
        uvel(iw)=c0
        vvel(iw)=c0
        str1(iw)=c0
        str2(iw)=c0
        str3(iw)=c0
        str4(iw)=c0
        str5(iw)=c0
        str6(iw)=c0
        str7(iw)=c0
        str8(iw)=c0
     enddo
     !$omp end parallel do

  end subroutine numainit

!=======================================================================

  subroutine halo_sweep(domark)
     ! Walk every ghost point the boundary conditions fill, in the order the
     ! 2-D halo fills them.
     !
     ! One sweep, two uses.  With domark it marks the cells the halo update
     ! reads and writes, so calc_navel gives them slots; without it, it builds
     ! the list.  Sharing the loop is the point: a cell the build pass needs
     ! but the mark pass missed would drop out of the list with no sign of it.
     !
     ! Order follows ice_boundary.F90.  The cyclic exchange runs first in both
     ! directions, then zero_gradient/linear_extrap east/west over the full
     ! column and north/south over the full row.  Running each pass the full
     ! length of the edge is what fills the corners: the north/south pass
     ! reads a ghost the east/west pass has just written, so a bi-cyclic
     ! corner ends up with the diagonal cell, which is what the 2-D exchange
     ! gets from its diagonal neighbour.  The list is therefore ordered, and
     ! the loop applying it must stay sequential.

     use ice_blocks, only: ew_boundary_type, ns_boundary_type

     implicit none

     logical(kind=log_kind), intent(in) :: domark

     ! local variables

     integer(kind=int_kind) :: i, j, g, ilo, ihi, jlo, jhi
     real   (kind=dbl_kind) :: w

     character(len=*), parameter :: subname = '(halo_sweep)'

     if (.not. domark) n_halo_bc = 0

     ! the gathered field is one block, so its physical domain is the whole
     ! grid less the ghost rim
     ilo = 1  + nghost
     ihi = nx - nghost
     jlo = 1  + nghost
     jhi = ny - nghost

     ! the cyclic exchange, both directions
     if (trim(ew_boundary_type) == 'cyclic') then
        do j = 1, ny
        do g = 1, nghost
           call halo_point(g      , j, ihi-nghost+g, j, ihi-nghost+g, j, c0, domark)
           call halo_point(ihi + g, j, ilo+g-1     , j, ilo+g-1     , j, c0, domark)
        end do
        end do
     endif

     if (trim(ns_boundary_type) == 'cyclic') then
        do i = 1, nx
        do g = 1, nghost
           call halo_point(i, g      , i, jhi-nghost+g, i, jhi-nghost+g, c0, domark)
           call halo_point(i, jhi + g, i, jlo+g-1     , i, jlo+g-1     , c0, domark)
        end do
        end do
     endif

     ! then the extrapolating types, east/west before north/south
     if (trim(ew_boundary_type) == 'zero_gradient' .or. &
         trim(ew_boundary_type) == 'linear_extrap') then
        do j = 1, ny
        do g = 1, nghost
           w = c0
           if (trim(ew_boundary_type) == 'linear_extrap') w = real(nghost-g+1, dbl_kind)
           call halo_point(g      , j, ilo, j, ilo+1, j, w, domark)   ! west
           w = c0
           if (trim(ew_boundary_type) == 'linear_extrap') w = real(g, dbl_kind)
           call halo_point(ihi + g, j, ihi, j, ihi-1, j, w, domark)   ! east
        end do
        end do
     endif

     if (trim(ns_boundary_type) == 'zero_gradient' .or. &
         trim(ns_boundary_type) == 'linear_extrap') then
        do i = 1, nx
        do g = 1, nghost
           w = c0
           if (trim(ns_boundary_type) == 'linear_extrap') w = real(nghost-g+1, dbl_kind)
           call halo_point(i, g      , i, jlo, i, jlo+1, w, domark)   ! south
           w = c0
           if (trim(ns_boundary_type) == 'linear_extrap') w = real(g, dbl_kind)
           call halo_point(i, jhi + g, i, jhi, i, jhi-1, w, domark)   ! north
        end do
        end do
     endif

  end subroutine halo_sweep

!=======================================================================

  subroutine halo_point(id, jd, i1, j1, i2, j2, w, domark)
     ! Mark, or append, one halo point.
     !
     ! Marking is what lets the append succeed: these cells are the halo
     ! update's own inputs and outputs, and a ghost whose interior neighbour
     ! is land is reached by no kernel, so without marking it has no slot and
     ! the point would be dropped.  The 2-D halo fills it regardless, which is
     ! where evp1d used to differ.

     implicit none

     integer(kind=int_kind), intent(in) :: id, jd, i1, j1, i2, j2
     real   (kind=dbl_kind), intent(in) :: w
     logical(kind=log_kind), intent(in) :: domark

     ! local variables

     integer(kind=int_kind) :: sd, s1, s2

     character(len=*), parameter :: subname = '(halo_point)'

     if (domark) then
        call mark_cell(id, jd)
        call mark_cell(i1, j1)
        if (w /= c0) call mark_cell(i2, j2)
        return
     endif

     sd = ijslot(id + (jd - 1) * nx)
     s1 = ijslot(i1 + (j1 - 1) * nx)
     if (sd == 0 .or. s1 == 0) return

     s2 = ijslot(i2 + (j2 - 1) * nx)
     if (s2 == 0) then
        ! w == 0 does not read s2, so a missing one is no obstacle
        if (w /= c0) return
        s2 = s1
     endif

     n_halo_bc = n_halo_bc + 1
     halo_bc_dst(n_halo_bc) = sd
     halo_bc_s1 (n_halo_bc) = s1
     halo_bc_s2 (n_halo_bc) = s2
     halo_bc_w  (n_halo_bc) = w

  end subroutine halo_point

!=======================================================================

  subroutine evp1d_halo_update()

     implicit none
     integer(kind=int_kind) :: iw

     character(len=*), parameter :: subname = '(evp1d_halo_update)'

! One sequential pass over the list built in halo_sweep.  Not parallelised
! and not to be: north/south destinations are east/west sources at the
! corners, so the entries are order-dependent.  See halo_sweep.
     do iw = 1, n_halo_bc
        uvel(halo_bc_dst(iw)) = uvel(halo_bc_s1(iw)) + halo_bc_w(iw) *      &
                               (uvel(halo_bc_s1(iw)) - uvel(halo_bc_s2(iw)))
        vvel(halo_bc_dst(iw)) = vvel(halo_bc_s1(iw)) + halo_bc_w(iw) *      &
                               (vvel(halo_bc_s1(iw)) - vvel(halo_bc_s2(iw)))
     end do

  end subroutine evp1d_halo_update

end module ice_dyn_evp1d

