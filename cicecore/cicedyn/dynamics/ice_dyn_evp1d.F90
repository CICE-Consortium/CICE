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

  ! halo updates for circular domains
  integer(kind=int_kind), allocatable, dimension(:)   ::                       &
    halo_parent_outer_east , halo_parent_outer_west ,                          &
    halo_parent_outer_north, halo_parent_outer_south,                          &
    halo_inner_east        , halo_inner_west        ,                          &
    halo_inner_north       , halo_inner_south

  ! number of halo points (same for inner and outer)
  integer(kind=int_kind)                              ::                       &
    n_inner_east, n_inner_west, n_inner_north, n_inner_south

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

    integer(kind=int_kind) :: ios, ierr

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
      call numainit(1,nActive,navel)
      call convert_2d_1d_init(nActive,G_HTE, G_HTN, G_uarear, G_dxT, G_dyT)
      call evp1d_alloc_static_halo()
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

       call calc_halo_parent(Nactive,navel)

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
       str1(1:navel)=c0
       str2(1:navel)=c0
       str3(1:navel)=c0
       str4(1:navel)=c0
       str5(1:navel)=c0
       str6(1:navel)=c0
       str7(1:navel)=c0
       str8(1:navel)=c0

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

    allocate(str1(1:navel0)   , str2(1:navel0), str3(1:navel0), &
             str4(1:navel0)   , str5(1:navel0), str6(1:navel0), &
             str7(1:navel0)   , str8(1:navel0),                 &
             indxTij(1:navel0), uvel(1:navel0), vvel(1:navel0), &
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

    ! allocation of arrays to use for halo
    ! These are the size of one of the dimensions of the global grid but they could be
    ! reduced in size as only the number of active U points are used.
    ! Points to send data from are in the "inner" vectors. Data in outer points are named "outer"

    allocate(halo_inner_east (ny), halo_inner_west (ny), &
             halo_inner_north(nx), halo_inner_south(nx), &
             stat=ierr)

    if (ierr/=0) then
       call abort_ice(subname//' ERROR: allocating', file=__FILE__, line=__LINE__)
    endif

    allocate(halo_parent_outer_east (ny), halo_parent_outer_west (ny), &
             halo_parent_outer_north(nx), halo_parent_outer_south(nx), &
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

!=============================================================================

  subroutine union(x, y, xdim, ydim, xy, nxy)

    ! Find union (xy) of two sorted integer vectors (x and y), i.e.
    ! combined values of the two vectors with no repetitions
    implicit none
    integer(kind=int_kind), intent(in)  :: xdim, ydim
    integer(kind=int_kind), intent(in)  :: x(1:xdim), y(1:ydim)
    integer(kind=int_kind), intent(out) :: xy(1:xdim + ydim)
    integer(kind=int_kind), intent(out) :: nxy

    ! local variables

    integer(kind=int_kind) :: i, j, k
    character(len=*), parameter :: subname = '(union)'

    i = 1
    j = 1
    k = 1
    do while (i <= xdim .and. j <= ydim)
       if (x(i) < y(j)) then
          xy(k) = x(i)
          i = i + 1
       else if (x(i) > y(j)) then
          xy(k) = y(j)
          j = j + 1
       else
          xy(k) = x(i)
          i = i + 1
          j = j + 1
       endif
       k = k + 1
    enddo

    ! the rest
    do while (i <= xdim)
       xy(k) = x(i)
       i = i + 1
       k = k + 1
    enddo
    do while (j <= ydim)
      xy(k) = y(j)
      j = j + 1
      k = k + 1
    enddo
    nxy = k - 1

  end subroutine union

!=============================================================================

  subroutine gather_static(G_uarear, G_dxT, G_dyT, G_Tmask)

     ! In standalone  distrb_info is an integer. Not needed anyway
     use ice_communicate, only : master_task
     use ice_gather_scatter, only : gather_global_ext
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
     call gather_global_ext(G_uarear, uarear, master_task, distrb_info)
     call gather_global_ext(G_dxT   , dxT   , master_task, distrb_info)
     call gather_global_ext(G_dyT   , dyT   , master_task, distrb_info)
     call gather_global_ext(G_Tmask , Tmask , master_task, distrb_info)

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
     use ice_gather_scatter, only : gather_global_ext
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
     call gather_global_ext(G_stressp_1 ,     L_stressp_1,     master_task, distrb_info,c0)
     call gather_global_ext(G_stressp_2 ,     L_stressp_2,     master_task, distrb_info,c0)
     call gather_global_ext(G_stressp_3 ,     L_stressp_3,     master_task, distrb_info,c0)
     call gather_global_ext(G_stressp_4 ,     L_stressp_4,     master_task, distrb_info,c0)

     call gather_global_ext(G_stressm_1 ,     L_stressm_1,     master_task, distrb_info,c0)
     call gather_global_ext(G_stressm_2 ,     L_stressm_2,     master_task, distrb_info,c0)
     call gather_global_ext(G_stressm_3 ,     L_stressm_3,     master_task, distrb_info,c0)
     call gather_global_ext(G_stressm_4 ,     L_stressm_4,     master_task, distrb_info,c0)

     call gather_global_ext(G_stress12_1,     L_stress12_1,    master_task, distrb_info,c0)
     call gather_global_ext(G_stress12_2,     L_stress12_2,    master_task, distrb_info,c0)
     call gather_global_ext(G_stress12_3,     L_stress12_3,    master_task, distrb_info,c0)
     call gather_global_ext(G_stress12_4,     L_stress12_4,    master_task, distrb_info,c0)
     call gather_global_ext(G_strength  ,     L_strength  ,    master_task, distrb_info,c0)

     call gather_global_ext(G_cdn_ocn   ,     L_cdn_ocn   ,     master_task, distrb_info)
     call gather_global_ext(G_aiu       ,     L_aiu       ,     master_task, distrb_info)
     call gather_global_ext(G_uocn      ,     L_uocn      ,     master_task, distrb_info)
     call gather_global_ext(G_vocn      ,     L_vocn      ,     master_task, distrb_info)

     call gather_global_ext(G_waterxU   ,     L_waterxU   ,     master_task, distrb_info)
     call gather_global_ext(G_wateryU   ,     L_wateryU   ,     master_task, distrb_info)
     call gather_global_ext(G_forcexU   ,     L_forcexU   ,     master_task, distrb_info)
     call gather_global_ext(G_forceyU   ,     L_forceyU   ,     master_task, distrb_info)

     call gather_global_ext(G_umassdti  ,     L_umassdti  ,     master_task, distrb_info)
     call gather_global_ext(G_fmU       ,     L_fmU       ,     master_task, distrb_info)

     call gather_global_ext(G_Tbu       ,     L_Tbu       ,     master_task, distrb_info)
     call gather_global_ext(G_uvel      ,     L_uvel      ,     master_task, distrb_info,c0)
     call gather_global_ext(G_vvel      ,     L_vvel      ,     master_task, distrb_info,c0)
     call gather_global_ext(G_iceTmask  ,     L_iceTmask  ,     master_task, distrb_info)
     call gather_global_ext(G_iceUmask  ,     L_iceUmask  ,     master_task, distrb_info)

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
     use ice_gather_scatter, only : scatter_global_ext
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

     call scatter_global_ext(L_stressp_1,  G_stressp_1,  master_task, distrb_info)
     call scatter_global_ext(L_stressp_2,  G_stressp_2,  master_task, distrb_info)
     call scatter_global_ext(L_stressp_3,  G_stressp_3,  master_task, distrb_info)
     call scatter_global_ext(L_stressp_4,  G_stressp_4,  master_task, distrb_info)

     call scatter_global_ext(L_stressm_1,  G_stressm_1,  master_task, distrb_info)
     call scatter_global_ext(L_stressm_2,  G_stressm_2,  master_task, distrb_info)
     call scatter_global_ext(L_stressm_3,  G_stressm_3,  master_task, distrb_info)
     call scatter_global_ext(L_stressm_4,  G_stressm_4,  master_task, distrb_info)

     call scatter_global_ext(L_stress12_1, G_stress12_1, master_task, distrb_info)
     call scatter_global_ext(L_stress12_2, G_stress12_2, master_task, distrb_info)
     call scatter_global_ext(L_stress12_3, G_stress12_3, master_task, distrb_info)
     call scatter_global_ext(L_stress12_4, G_stress12_4, master_task, distrb_info)

     call scatter_global_ext(L_strintxU  , G_strintxU  , master_task, distrb_info)
     call scatter_global_ext(L_strintyU  , G_strintyU  , master_task, distrb_info)
     call scatter_global_ext(L_uvel      , G_uvel      , master_task, distrb_info)
     call scatter_global_ext(L_vvel      , G_vvel      , master_task, distrb_info)
     call scatter_global_ext(L_taubxU    , G_taubxU    , master_task, distrb_info)
     call scatter_global_ext(L_taubyU    , G_taubyU    , master_task, distrb_info)

  end subroutine scatter_dyn

!=============================================================================

  subroutine convert_2d_1d_init(na0, G_HTE, G_HTN, G_uarear,  G_dxT, G_dyT)

     implicit none

     integer(kind=int_kind), intent(in) ::  na0
     real (kind=dbl_kind), dimension(:, :), intent(in) :: G_HTE, G_HTN, G_uarear, G_dxT, G_dyT

     ! local variables

     integer(kind=int_kind) :: iw, lo, up, j, i
     integer(kind=int_kind), dimension(1:na0) :: &
       Iin, Iee, Ine, Ise, Inw, Isw, Isse

     integer(kind=int_kind), dimension(1:7 * na0) :: util1, util2

     character(len=*), parameter :: subname = '(convert_2d_1d_init)'

     ! calculate additional 1D indices used for finite differences
     do iw = 1, na0
     ! get 2D indices
        i = indxti(iw)
        j = indxtj(iw)
     ! calculate 1D indices
        Iin(iw)  = i     + (j - 1) * nx  ! ( 0, 0) target point
        Iee(iw)  = i - 1 + (j - 1) * nx  ! (-1, 0)
        Ine(iw)  = i - 1 + (j - 2) * nx  ! (-1,-1)
        Ise(iw)  = i     + (j - 2) * nx  ! ( 0,-1)
        Inw(iw)  = i + 1 + (j - 1) * nx  ! (+1, 0)
        Isw(iw)  = i + 1 + (j - 0) * nx  ! (+1,+1)
        Isse(iw) = i     + (j - 0) * nx  ! ( 0,+1)
     end do

     ! find number of points needed for finite difference calculations
     call union(Iin,   Iee,  na0, na0, util1,i    )
     call union(util1, Ine,  i,  na0, util2, j    )
     call union(util2, Ise,  j,  na0, util1, i    )
     call union(util1, Inw,  i,  na0, util2, j    )
     call union(util2, Isw,  j,  na0, util1, i    )
     call union(util1, Isse, i,  na0, util2, navel)

     ! index vector with sorted target points
     do iw = 1, na0
        indxTij(iw) = Iin(iw)
     end do
     ! sorted additional points
     call setdiff(util2, Iin, navel, na0, util1, j)
     do iw = na0 + 1, navel
         indxTij(iw) = util1(iw - na0)
     end do

     ! indices for additional points needed for uvel and vvel
     call findXinY(Iee,  indxTij, na0, navel, ee)
     call findXinY(Ine,  indxTij, na0, navel, ne)
     call findXinY(Ise,  indxTij, na0, navel, se)
     call findXinY(Inw,  indxTij, na0, navel, nw)
     call findXinY(Isw,  indxTij, na0, navel, sw)
     call findXinY(Isse, indxTij, na0, navel, sse)
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
     G_strintxU   = c0
     G_strintyU   = c0
     G_Tbu        = c0
     G_uvel       = c0
     G_vvel       = c0
     G_taubxU     = c0
     G_taubyU     = c0

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

!=======================================================================

  subroutine setdiff(x, y,  lvecx, lvecy,xy, nxy)
     ! Find element (xy) of two sorted integer vectors (x and y) that
     ! are in x, but not in y, or in y, but not in x

     implicit none

     integer(kind=int_kind), intent(in) :: lvecx,lvecy
     integer(kind=int_kind), intent(in) :: x(1:lvecx), y(1:lvecy)
     integer(kind=int_kind), intent(out) :: xy(1:lvecx + lvecy)
     integer(kind=int_kind), intent(out) :: nxy

     ! local variables

     integer(kind=int_kind) :: i, j, k

     character(len=*), parameter :: subname = '(setdiff)'

     i = 1
     j = 1
     k = 1
     do while (i <= lvecx .and. j <= lvecy)
        if (x(i) < y(j)) then
           xy(k) = x(i)
           i = i + 1
           k = k + 1
        else if (x(i) > y(j)) then
           xy(k) = y(j)
           j = j + 1
           k = k + 1
        else
           i = i + 1
           j = j + 1
        end if
     end do

     ! the rest
     do while (i <= lvecx)
        xy(k) = x(i)
        i = i + 1
        k = k + 1
     end do
     do while (j <= lvecy)
        xy(k) = y(j)
        j = j + 1
        k = k + 1
     end do
     nxy = k - 1

  end subroutine setdiff

!=======================================================================

  subroutine findXinY(x, y, lvecx, lvecy, indx)
     ! Find indx vector so that x(1:na) = y(indx(1:na))
     !
     !  Conditions:
     !   * EVERY item in x is found in y
     !   * x(1:lvecx) is a sorted integer vector
     !   * y(1:lvecy) consists of two sorted integer vectors:
     !        [y(1:lvecx); y(lvecy + 1:lvecx)]
     !   * lvecy >= lvecx

     implicit none

     integer (kind=int_kind), intent(in) :: lvecx, lvecy
     integer (kind=int_kind), intent(in) :: x(1:lvecx), y(1:lvecy)
     integer (kind=int_kind), intent(out) :: indx(1:lvecx)

      ! local variables

     integer (kind=int_kind) :: i, j1, j2

     character(len=*), parameter :: subname = '(findXinY)'

     i = 1
     j1 = 1
     j2 = lvecx + 1
     do while (i <= lvecx)
        if (x(i) == y(j1)) then
           indx(i) = j1
           i = i + 1
           j1 = j1 + 1
        else if (x(i) == y(j2)) then
           indx(i) = j2
           i = i + 1
           j2 = j2 + 1
        else if (x(i) > y(j1)) then
           j1 = j1 + 1
        else if (x(i) > y(j2)) then
           j2 = j2 + 1
        else
           stop
        end if
     end do

  end subroutine findXinY

!=======================================================================

  subroutine calc_navel(na0, navel0)
     ! Calculate number of active points, including halo points

     implicit none

     integer(kind=int_kind), intent(in) :: na0
     integer(kind=int_kind), intent(out) :: navel0

     ! local variables

     integer(kind=int_kind) :: iw, i, j
     integer(kind=int_kind), dimension(1:na0) :: &
        Iin, Iee, Ine, Ise, Inw, Isw, Isse, indi, indj

     integer(kind=int_kind), dimension(1:7 * na0) :: util1, util2

     character(len=*), parameter :: subname = '(calc_navel)'

     ! calculate additional 1D indices used for finite differences
     do iw = 1, na0
        ! get 2D indices
        i = indxti(iw)
        j = indxtj(iw)

        ! calculate 1D indices
        Iin(iw)  = i     + (j - 1) * nx  ! ( 0,  0) target point
        Iee(iw)  = i - 1 + (j - 1) * nx  ! (-1,  0)
        Ine(iw)  = i - 1 + (j - 2) * nx  ! (-1, -1)
        Ise(iw)  = i     + (j - 2) * nx  ! ( 0, -1)
        Inw(iw)  = i + 1 + (j - 1) * nx  ! (+1,  0)
        Isw(iw)  = i + 1 + (j - 0) * nx  ! (+1, +1)
        Isse(iw) = i +     (j - 0) * nx  ! ( 0, +1)
     end do

     ! find number of points needed for finite difference calculations
     call union(Iin  , Iee , na0, na0, util1, i     )
     call union(util1, Ine , i  , na0, util2, j     )
     call union(util2, Ise , j  , na0, util1, i     )
     call union(util1, Inw , i  , na0, util2, j     )
     call union(util2, Isw , j  , na0, util1, i     )
     call union(util1, Isse, i  , na0, util2, navel0)

  end subroutine calc_navel

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

  subroutine evp1d_halo_update()

     implicit none
     integer(kind=int_kind) :: iw

     character(len=*), parameter :: subname = '(evp1d_halo_update)'

!TILL    !$omp parallel do schedule(runtime) private(iw)
     do iw = 1, n_inner_east
        uvel(halo_parent_outer_east(iw)) = uvel(halo_inner_east(iw))
        vvel(halo_parent_outer_east(iw)) = vvel(halo_inner_east(iw))
     end do
! western halo
     do iw = 1, n_inner_west
        uvel(halo_parent_outer_west(iw)) = uvel(halo_inner_west(iw))
        vvel(halo_parent_outer_west(iw)) = vvel(halo_inner_west(iw))
     end do
     do iw = 1, n_inner_south
        uvel(halo_parent_outer_south(iw)) = uvel(halo_inner_south(iw))
        vvel(halo_parent_outer_south(iw)) = vvel(halo_inner_south(iw))
     end do
! western halo
     do iw = 1, n_inner_north
        uvel(halo_parent_outer_north(iw)) = uvel(halo_inner_north(iw))
        vvel(halo_parent_outer_north(iw)) = vvel(halo_inner_north(iw))
     end do

  end subroutine evp1d_halo_update

!=======================================================================

  subroutine calc_halo_parent(na0,navel0)
     ! splits the global domain in east and west boundary and find the inner (within) the domain and the outer (outside the domain)
     ! Implementation for circular boundaries. This means that mathes between the opposite directions must be found
     ! E.g. inner_west and outer_east
     ! Till Rasmussen, DMI 2023

     use ice_domain, only: ew_boundary_type, ns_boundary_type
     implicit none

     integer(kind=int_kind), intent(in) :: na0, navel0

     ! local variables

     ! Indexes, Directions are east, weast, north and south
     ! This is done to reduce the search windows.
     ! Iw runs from 1 to navel and the one to keep in the end
     ! Iw_inner_{direction} contains the indexes for

     integer(kind=int_kind) :: &
        iw, n_outer_east, n_outer_west, n_outer_south, n_outer_north

     integer(kind=int_kind) :: i, j, ifind, jfind !  2d index. ifind and jfind are points on the boundary

     integer(kind=int_kind), dimension(ny) :: &
        halo_outer_east, halo_outer_west,      &
        ind_inner_west , ind_inner_east

     integer(kind=int_kind), dimension(nx) :: &
        halo_outer_south, halo_outer_north,    &
        ind_inner_south , ind_inner_north

     character(len=*), parameter :: subname = '(calc_halo_parent)'

     !-----------------------------------------------------------------
     ! Indices for halo update:
     !     0: no halo point
     !    >0: index for halo point parent, related to indij vector
     !
     ! TODO: Implement for nghost > 1
     ! TODO: Implement for tripole grids
     !-----------------------------------------------------------------
     halo_inner_west(:) = 0
     halo_inner_east(:) = 0
     halo_inner_south(:) = 0
     halo_inner_north(:) = 0

     halo_outer_west(:) = 0
     halo_outer_east(:) = 0
     halo_outer_south(:) = 0
     halo_outer_north(:) = 0

     ind_inner_west(:)  = 0
     ind_inner_east(:)  = 0
     ind_inner_south(:)  = 0
     ind_inner_north(:)  = 0

     halo_parent_outer_east(:)=0
     halo_parent_outer_west(:)=0
     halo_parent_outer_north(:)=0
     halo_parent_outer_south(:)=0
     ! Index inner boundary
     n_inner_north=0
     n_inner_south=0
     n_inner_east=0
     n_inner_west=0
     ! Index outer boundary
     n_outer_east=0
     n_outer_west=0
     n_outer_north=0
     n_outer_south=0
     !TILL SHOULD CHANGE TO 1D
     do iw = 1, na0
        j = int((indxTij(iw) - 1) / (nx)) + 1
        i = indxTij(iw) - (j - 1) * nx
        ! All four boundaries find points internally that are within the domain and next to the boundary
        ! This can in principle be moved to previos loops that connects i and j to 1d index.
        ! ifind is i value on the halo to find.
        ! Some parts assume nghost = 1
        ! INNER EAST
        if (trim(ew_boundary_type) == 'cyclic') then
           if ((.not. skipUcell(iw)) .and. (i==nx-nghost)) then
              n_inner_east=n_inner_east+1
              ifind = 1
              ind_inner_east(n_inner_east)  = ifind     + (j - 1) * nx
              halo_inner_east(n_inner_east) = iw
           else if ((.not. skipUcell(iw)) .and. (i==1+nghost)) then
              n_inner_west=n_inner_west+1
              ifind = nx
              ind_inner_west(n_inner_west)  = ifind     + (j - 1) * nx
              halo_inner_west(n_inner_west) = iw
           endif
        endif
        if (trim(ns_boundary_type) == 'cyclic') then
           if ((.not. skipUcell(iw)) .and. (j==1+nghost)) then
              n_inner_south=n_inner_south+1
              jfind = ny
              ind_inner_south(n_inner_south)  = i     + (jfind - 1) * nx
              halo_inner_south(n_inner_south) = iw
           else if ((.not. skipUcell(iw)) .and. (j==ny-nghost)) then
              n_inner_north=n_inner_north+1
              jfind = 1
              ind_inner_north(n_inner_north)  = i     + (jfind - 1) * nx
              halo_inner_north(n_inner_north) = iw
           endif
         endif
         ! Finds all halos points on western halo WEST
         if (i == 1) then
            n_outer_west=n_outer_west+1
            halo_outer_west(n_outer_west)= iw
         endif
         ! Simiilar on East
         if (i == nx ) then
            n_outer_east=n_outer_east+1
            halo_outer_east(n_outer_east)=iw
         endif
         ! Finds all halos points on western halo WEST
         if (j == 1) then
            n_outer_south=n_outer_south+1
            halo_outer_south(n_outer_south)= iw
         endif
         ! Simiilar on East
         if (j == ny ) then
            n_outer_north=n_outer_north+1
            halo_outer_north(n_outer_north)=iw
         endif
     end do

     ! outer halo also needs points that are not active
     do iw = na0+1, navel0
        j = int((indxTij(iw) - 1) / (nx)) + 1
        i = indxTij(iw) - (j - 1) * nx
        ! outer halo west
         if (i == 1) then
            n_outer_west=n_outer_west+1
            halo_outer_west(n_outer_west)= iw
         endif
        ! outer halo east
         if (i == nx ) then
            n_outer_east=n_outer_east+1
            halo_outer_east(n_outer_east)=iw
         endif
        ! outer halo south
         if (j == 1) then
            n_outer_south=n_outer_south+1
            halo_outer_south(n_outer_south)= iw
         endif
        ! outer halo north
         if (j == ny ) then
            n_outer_north=n_outer_north+1
            halo_outer_north(n_outer_north)=iw
         endif
     end do
     ! Search is now reduced to a search between two reduced vectors for each boundary
     ! This runs through each boundary and matches
     ! number of active points for halo east and west (count of active u cells within the domain.
     ! reduce outer array to only match inner arrays
     ! East West
     if (trim(ew_boundary_type) == 'cyclic') then
        do i=1,n_inner_west
           do j=1,n_outer_east
              if (ind_inner_west(i) == indxTij(halo_outer_east(j))) then
                 halo_parent_outer_west(i)=halo_outer_east(j)
              endif
           end do
        end do

        do i=1,n_inner_east
           do j=1,n_outer_west
              if (ind_inner_east(i) == indxTij(halo_outer_west(j))) then
                 halo_parent_outer_east(i)=halo_outer_west(j)
              endif
           end do
        end do
     endif
     if (trim(ns_boundary_type) == 'cyclic') then
        do i=1,n_inner_south
           do j=1,n_outer_north
              if (ind_inner_south(i) == indxTij(halo_outer_north(j))) then
                 halo_parent_outer_south(i)=halo_outer_north(j)
              endif
           end do
        end do

        do i=1,n_inner_north
           do j=1,n_outer_south
              if (ind_inner_north(i) == indxTij(halo_outer_south(j))) then
                 halo_parent_outer_north(i)=halo_outer_south(j)
              endif
           end do
        end do
     endif

  end subroutine calc_halo_parent

!=======================================================================

end module ice_dyn_evp1d

