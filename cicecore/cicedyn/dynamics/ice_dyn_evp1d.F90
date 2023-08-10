!=============================================================================== 
module ice_dyn_evp1d
  !- modules -------------------------------------------------------------------
  use ice_kinds_mod
  use ice_constants
  !  none so far ...

  !- directives ----------------------------------------------------------------
  implicit none
  private

  !- interfaces ----------------------------------------------------------------
  !  none so far ...

  !- public routines -----------------------------------------------------------
  public :: dyn_evp1d_init, dyn_evp1d_run, dyn_evp1d_finalize, dyn_evp2d_dump

  !- private routines ----------------------------------------------------------
  private :: io_new_unit, evp_1d_nml_reader

  !- public vars ---------------------------------------------------------------
  !  none so far ...

  !- private vars --------------------------------------------------------------
  logical(kind=log_kind), save :: ldebug, linfo   ! output
  logical(kind=log_kind), save :: firsttime = .true.  ! is this called for the first time? Should be called from 
  integer(kind=int_kind), save :: ts  ! time step for debugging
  integer(kind=int_kind), save :: iu06 ! log file
  ! nx and ny are module variables for arrays after gather (G_*) Dimension according to CICE is nx_global+2*nghost,
  ! ny_global+2*nghost
  ! nactive are number of active points (both t and u). navel is number of active
  integer(kind=int_kind), save :: nx, ny, nActive, navel, nallocated, buf1d
  ! cice variables imported
  integer(kind=int_kind), save :: nghost, nx_block, ny_block, max_block
  integer(kind=int_kind), allocatable, dimension(:,:) :: iwidx
  logical(kind=log_kind), allocatable, dimension(:) :: unionTU
  logical(kind=log_kind), allocatable, dimension(:)   :: skipTcell,skipUcell
  real(kind=dbl_kind),    allocatable, dimension(:)   ::                       &
    HTE,HTN, HTEm1,HTNm1
  integer(kind=int_kind), allocatable, dimension(:) :: ee,ne,se,nw,sw,sse ! arrays for neighbour points
  integer(kind=int_kind), allocatable, dimension(:) ::                         &
    indxti, indxtj, indxtij, indxui, indxuj
! icepack - ocean density

    real(kind=dbl_kind), parameter :: rhow      = 1026.0_dbl_kind ! density of seawater (kg/m^3)
! 1D arrays to allocate
  real (kind=dbl_kind),   allocatable, dimension(:) ::                         &
    cdn_ocn,aiu,uocn,vocn,waterxU,wateryU,forcexU,forceyU,umassdti,fmU,uarear, &
    strintxU,strintyU,uvel_init,vvel_init,                                     &
!    cdn_ocn,aiu,uocn,vocn,forcexU,forceyU,umassdti,fmU,uarear,                 &
!    uvel_init,vvel_init,                                                       &
    strength, uvel, vvel, dxt, dyt,                                            &
    stressp_1, stressp_2, stressp_3, stressp_4, stressm_1, stressm_2,          &
    stressm_3, stressm_4, stress12_1, stress12_2, stress12_3, stress12_4,      &
    str1, str2, str3, str4, str5, str6, str7, str8, Tbu, Cb

  contains

  !=============================================================================
  ! module public subroutines 

  !=============================================================================
  ! This follows CICE naming
  ! In addition all water points are assumed to be active and allocated thereafter.
  subroutine dyn_evp1d_init(cnx, cny, cnx_block, cny_block, cmax_block, cnghost, &
                            L_dyT, L_dxT, L_uarear, L_tmask,                     &
                            G_HTE, G_HTN)
! Note that TMask is ocean/land
  use debug_evp1d, only : dump_init
    implicit none
! name within 
    integer(kind=int_kind), intent(in) :: cnx, cny, cnx_block, cny_block, cmax_block, cnghost
    real(kind=dbl_kind), dimension(:,:,:), intent(in) :: L_dyT, L_dxT, L_uarear
    logical(kind=log_kind), dimension(:,:,:), intent(in) :: L_tmask
    real(kind=dbl_kind), dimension(:,:)  , intent(in) :: G_HTE,G_HTN
    real(kind=dbl_kind), dimension(:,:), allocatable :: G_dyT, G_dxT, G_uarear
    logical(kind=log_kind), dimension(:,:), allocatable :: G_tmask
    integer(kind=int_kind) :: ios, ierr
    character(9), parameter :: logfilename='evp1d.log'
    iu06 = io_new_unit()
    open (iu06, file=logfilename,status='replace', iostat=ios)
    if (ios /= 0) then
      write  (*,*) 'Unable to open '//trim(logfilename)
      stop
    endif
    write(iu06,'(a19)') 'Initializing evp1d'
    call evp_1d_nml_reader()
! save as module vars
    ts=0
    nx_block=cnx_block
    ny_block=cny_block
    nghost=cnghost
    nx=cnx+2*nghost
    ny=cny+2*nghost
    max_block=cmax_block
    allocate(G_dyT(nx,ny),G_dxT(nx,ny),G_uarear(nx,ny),G_tmask(nx,ny),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
    write(iu06,*) nx, ny, nActive, nghost
    call gather_static(L_uarear,  L_dxT, L_dyT, L_Tmask, &
                       G_uarear,  G_dxT, G_dyT, G_Tmask)
    ! calculate number of water points (T and U). Only needed for the static version
    ! Tmask in ocean/ice
    call calc_nActiveTU(G_Tmask,nActive)
    call evp1d_alloc_static_na(nActive)
    call calc_2d_indices_init(nActive, G_Tmask)
    call calc_navel(nActive, navel)
    call evp1d_alloc_static_navel(navel)
    call convert_2d_1d_init(nActive,G_HTE, G_HTN, G_uarear, G_dxT, G_dyT)
    write(iu06,*) nx, ny, nActive
    ! activate dump
    call dump_init(iu06)
  end subroutine dyn_evp1d_init

  !=============================================================================
  subroutine dyn_evp1d_run(L_stressp_1 , L_stressp_2 , L_stressp_3, L_stressp_4,     &
                           L_stressm_1 , L_stressm_2 , L_stressm_3, L_stressm_4,     &
                           L_stress12_1, L_stress12_2, L_stress12_3,L_stress12_4,    &
                           L_cdn_ocn   , L_aiu       , L_uocn     , L_vocn     ,     &
                           L_waterxU   , L_wateryU   , L_forcexU  , L_forceyU  ,     &
                           L_umassdti  , L_fmU       , L_strintxU , L_strintyU ,     &
                           L_Tbu       , L_Cb        , L_uvel      , L_vvel,         &
                           L_icetmask , L_iceUmask                                   )
    use debug_evp1d, only : dumpall
    use ice_dyn_shared, only : ndte
! bench should be renamed to something better e.g evp1d_func
  use bench, only : stress, stepu
    implicit none
    real(kind=dbl_kind), dimension(:,:,:), intent(in) :: &
                           L_stressp_1 , L_stressp_2 , L_stressp_3, L_stressp_4,     &
                           L_stressm_1 , L_stressm_2 , L_stressm_3, L_stressm_4,     &
                           L_stress12_1, L_stress12_2, L_stress12_3,L_stress12_4,    &
                           L_cdn_ocn   , L_aiu       , L_uocn     , L_vocn     ,     &
                           L_waterxU   , L_wateryU   , L_forcexU  , L_forceyU  ,     &
                           L_umassdti  , L_fmU       , L_strintxU , L_strintyU ,     &
                           L_Tbu       , L_Cb        , L_uvel      , L_vvel
    logical, dimension (:,:,:), intent(in) :: L_iceUmask
    integer, dimension (:,:,:), intent(in) :: L_iceTmask
    real(kind=dbl_kind), dimension(nx,ny) :: &
                           G_stressp_1 , G_stressp_2 , G_stressp_3, G_stressp_4,     &
                           G_stressm_1 , G_stressm_2 , G_stressm_3, G_stressm_4,     &
                           G_stress12_1, G_stress12_2, G_stress12_3,G_stress12_4,    &
                           G_cdn_ocn   , G_aiu       , G_uocn     , G_vocn     ,     &
                           G_waterxU   , G_wateryU   , G_forcexU  , G_forceyU  ,     &
                           G_umassdti  , G_fmU       , G_strintxU , G_strintyU ,     &
                           G_Tbu       , G_Cb        , G_uvel     , G_vvel

    logical(kind=log_kind), dimension (nx,ny)  :: G_iceUmask, G_iceTmask_log !tar , iceTmask
    integer(kind=int_kind), dimension (nx,ny)  :: G_iceTmask
    integer(kind=int_kind) :: ksub
    character(10),parameter :: mydebugfile1='before1d'
    character(10),parameter :: mydebugfile2='after1d'
! FIXME NEEDED IN OLD CICE version that unittest are build on
! NOT NEEDED IN CONSTANT VERSION
    if (ldebug) then
      write(iu06,*) 'Handling time-step, number of Allocated points ', ts, nActive
    endif


   call gather_dyn(L_stressp_1 , L_stressp_2 , L_stressp_3, L_stressp_4,     &
                    L_stressm_1 , L_stressm_2 , L_stressm_3, L_stressm_4,     &
                    L_stress12_1, L_stress12_2, L_stress12_3,L_stress12_4,    &
                    L_cdn_ocn   , L_aiu       , L_uocn     , L_vocn     ,     &
                    L_waterxU   , L_wateryU   , L_forcexU  , L_forceyU  ,     &
                    L_umassdti  , L_fmU       , L_strintxU , L_strintyU ,     &
                    L_Tbu       , L_Cb        , L_uvel     , L_vvel     ,     &
                    L_icetmask  , L_iceUmask  ,                               &
                    G_stressp_1 , G_stressp_2 , G_stressp_3, G_stressp_4,     &
                    G_stressm_1 , G_stressm_2 , G_stressm_3, G_stressm_4,     &
                    G_stress12_1, G_stress12_2, G_stress12_3,G_stress12_4,    &
                    G_cdn_ocn   , G_aiu       , G_uocn     , G_vocn     ,     &
                    G_waterxU   , G_wateryU   , G_forcexU  , G_forceyU  ,     &
                    G_umassdti  , G_fmU       , G_strintxU , G_strintyU ,     &
                    G_Tbu       , G_Cb        , G_uvel     , G_vvel     ,     &
                    G_iceTmask,  G_iceUmask)
    call convert_Tint2log(G_iceTmask,G_iceTmask_log)
    call set_skipMe(G_iceTmask_log, G_iceUmask,nActive)
    call convert_2d_1d_dyn(nActive,                                                   &
                           G_stressp_1 , G_stressp_2 , G_stressp_3 ,  G_stressp_4,      &
                           G_stressm_1 , G_stressm_2 , G_stressm_3 ,  G_stressm_4,      &
                           G_stress12_1, G_stress12_2, G_stress12_3, G_stress12_4,      &
                           G_cdn_ocn   , G_aiu       , G_uocn     ,  G_vocn     ,      &
                           G_waterxU   , G_wateryU   , G_forcexU  , G_forceyU   ,      &
                           G_umassdti  , G_fmU       , G_strintxU , G_strintyU  ,      &
                           G_Tbu       , G_Cb        , G_uvel     , G_vvel)
    do ksub = 1,ndte        ! subcycling
       call stress (ee, ne, se, 1, nActive,                                      &
                   uvel, vvel, dxT, dyT, skipTcell, strength,                    &
                   HTE, HTN, HTEm1, HTNm1,                                       &
                   stressp_1,  stressp_2,  stressp_3,  stressp_4,                &
                   stressm_1,  stressm_2,  stressm_3,  stressm_4,                &
                   stress12_1, stress12_2, stress12_3, stress12_4,               &
                   str1, str2, str3, str4, str5, str6, str7, str8)

    !$OMP BARRIER

       call stepu (1, nActive, cdn_ocn, aiu, uocn, vocn,                          &
                  waterxU, wateryU, forcexU, forceyU, umassdti, fmU, uarear,      &
                  strintxU, strintyU, uvel_init, vvel_init, uvel, vvel,        &
                  str1, str2, str3, str4, str5, str6, str7, str8,            &
                   nw, sw, sse, skipUcell, Tbu, Cb)
    !$OMP BARRIER

    enddo
    call dumpall(mydebugfile1, ts, nx, ny, iu06,                             &
                       G_stressp_1 , G_stressp_2 , G_stressp_3, G_stressp_4, &
                       G_stressm_1 , G_stressm_2 , G_stressm_3, G_stressm_4, &
                       G_stress12_1, G_stress12_2, G_stress12_3,G_stress12_4,&
                       G_cdn_ocn   , G_aiu       , G_uocn     , G_vocn     , &
                       G_waterxU   , G_wateryU   , G_forcexU  , G_forceyU  , &
                       G_umassdti  , G_fmU       , G_strintxU , G_strintyU , &
                       G_Tbu       , G_Cb        , G_uvel     , G_vvel)

    call convert_1d_2d_dyn(nActive, &
                           G_stressp_1 , G_stressp_2 , G_stressp_3, G_stressp_4,     &
                           G_stressm_1 , G_stressm_2 , G_stressm_3, G_stressm_4,     &
                           G_stress12_1, G_stress12_2, G_stress12_3,G_stress12_4,    &
                           G_cdn_ocn   , G_aiu       , G_uocn     , G_vocn     ,     &
                           G_waterxU   , G_wateryU   , G_forcexU  , G_forceyU  ,     &
                           G_umassdti  , G_fmU       , G_strintxU , G_strintyU ,     &
                           G_Tbu       , G_Cb        , G_uvel     , G_vvel)

    call dumpall(mydebugfile2, ts, nx, ny, iu06,                              &
                       G_stressp_1 , G_stressp_2 , G_stressp_3, G_stressp_4, &
                       G_stressm_1 , G_stressm_2 , G_stressm_3, G_stressm_4, &
                       G_stress12_1, G_stress12_2, G_stress12_3,G_stress12_4,&
                       G_cdn_ocn   , G_aiu       , G_uocn     , G_vocn     , &
                       G_waterxU   , G_wateryU   , G_forcexU  , G_forceyU  , &
                       G_umassdti  , G_fmU       , G_strintxU , G_strintyU , &
                       G_Tbu       , G_Cb        , G_uvel     , G_vvel)

!call convert_Tint2log(iceTmask,iceTmask_log)
!call setSkipme(nActive, iceTmask_log, iceUmask)
! END FIXME
! calculate number of active points. allocate if initial or if array size should increase
!    call calc_nActiveTU(iceTmask_log,nActive, iceUmask)
!    if (nActiveold ==0) then ! first
!         call evp_1d_alloc(nActive, nActive,nx,ny)
!         nactiveold=nActive+buf1d ! allocate
!         call init_unionTU(nx, ny, iceTmask_log,iceUmask)
!    else if (nactiveold < nActive) then
!         write(iu06,*) 'Warning nActive is bigger than old allocation. Need to re allocate'
!         call evp_1d_dealloc() ! only deallocate if not first time step 
!         call evp_1d_alloc(nActive, nActive,nx,ny)
!         nactiveold=nActive+buf1d ! allocate
!         call init_unionTU(nx, ny, iceTmask_log,iceUmask)
!    endif
!#endif
!    call cp_2dto1d(nActive)   
! FIXME THIS IS THE LOGIC FOR RE ALLOCATION IF NEEDED
    ts=ts+1
!tar    call add_1d(nx, ny, natmp, iceTmask_log, iceUmask, ts)
  end subroutine dyn_evp1d_run

  !============================================================================
 subroutine dyn_evp2d_dump(L_stressp_1 , L_stressp_2 , L_stressp_3, L_stressp_4,     &
                            L_stressm_1 , L_stressm_2 , L_stressm_3, L_stressm_4,     &
                            L_stress12_1, L_stress12_2, L_stress12_3,L_stress12_4,    &
                            L_cdn_ocn   , L_aiu       , L_uocn     , L_vocn     ,     &
                            L_waterxU   , L_wateryU   , L_forcexU  , L_forceyU  ,     &
                            L_umassdti  , L_fmU       , L_strintxU , L_strintyU ,     &
                            L_Tbu       , L_Cb        , L_uvel     , L_vvel     ,     &
                            L_iceTmask  , L_iceUmask)
    use debug_evp1d, only : dumpall
    implicit none
    real(kind=dbl_kind), dimension(:,:,:), intent(in) :: &
                           L_stressp_1 , L_stressp_2 , L_stressp_3, L_stressp_4,     &
                           L_stressm_1 , L_stressm_2 , L_stressm_3, L_stressm_4,     &
                           L_stress12_1, L_stress12_2, L_stress12_3,L_stress12_4,    &
                           L_cdn_ocn   , L_aiu       , L_uocn     , L_vocn     ,     &
                           L_waterxU   , L_wateryU   , L_forcexU  , L_forceyU  ,     &
                           L_umassdti  , L_fmU       , L_strintxU , L_strintyU ,     &
                           L_Tbu       , L_Cb        , L_uvel     , L_vvel
    logical, dimension (:,:,:), intent(in) :: L_iceUmask
    integer, dimension (:,:,:), intent(in) :: L_iceTmask
    real(kind=dbl_kind), dimension(nx,ny) :: &
                           G_stressp_1 , G_stressp_2 , G_stressp_3, G_stressp_4,     &
                           G_stressm_1 , G_stressm_2 , G_stressm_3, G_stressm_4,     &
                           G_stress12_1, G_stress12_2, G_stress12_3,G_stress12_4,    &
                           G_cdn_ocn   , G_aiu       , G_uocn     , G_vocn     ,     &
                           G_waterxU   , G_wateryU   , G_forcexU  , G_forceyU  ,     &
                           G_umassdti  , G_fmU       , G_strintxU , G_strintyU ,     &
                           G_Tbu       , G_Cb        , G_uvel     , G_vvel
    character(10),parameter :: mydebugfile1='after2d'
    ! These are dummy here
    logical(kind=log_kind), dimension (nx,ny)  :: G_iceUmask
    integer(kind=int_kind), dimension (nx,ny)  :: G_iceTmask

    call gather_dyn(L_stressp_1 , L_stressp_2 , L_stressp_3, L_stressp_4,     &
                    L_stressm_1 , L_stressm_2 , L_stressm_3, L_stressm_4,     &
                    L_stress12_1, L_stress12_2, L_stress12_3,L_stress12_4,    &
                    L_cdn_ocn   , L_aiu       , L_uocn     , L_vocn     ,     &
                    L_waterxU   , L_wateryU   , L_forcexU  , L_forceyU  ,     &
                    L_umassdti  , L_fmU       , L_strintxU , L_strintyU ,     &
                    L_Tbu       , L_Cb        , L_uvel     , L_vvel     ,     &
                    L_iceTmask  , L_iceUmask  ,                               &
                    G_stressp_1 , G_stressp_2 , G_stressp_3, G_stressp_4,     &
                    G_stressm_1 , G_stressm_2 , G_stressm_3, G_stressm_4,     &
                    G_stress12_1, G_stress12_2, G_stress12_3,G_stress12_4,    &
                    G_cdn_ocn   , G_aiu       , G_uocn     , G_vocn     ,     &
                    G_waterxU   , G_wateryU   , G_forcexU  , G_forceyU  ,     &
                    G_umassdti  , G_fmU       , G_strintxU , G_strintyU ,     &
                    G_Tbu       , G_Cb        , G_uvel     , G_vvel     ,     &
                    G_iceTmask,  G_iceUmask)

! this is called after ts is updated
! this is called after ts is updated
    call dumpall(mydebugfile1, ts-1, nx, ny, iu06,                             &
                    G_stressp_1 , G_stressp_2 , G_stressp_3, G_stressp_4, &
                    G_stressm_1 , G_stressm_2 , G_stressm_3, G_stressm_4, &
                    G_stress12_1, G_stress12_2, G_stress12_3,G_stress12_4,&
                    G_cdn_ocn   , G_aiu       , G_uocn     , G_vocn     , &
                    G_waterxU   , G_wateryU   , G_forcexU  , G_forceyU  , &
                    G_umassdti  , G_fmU       , G_strintxU , G_strintyU , &
                    G_Tbu       , G_Cb        , G_uvel     , G_vvel)

    end subroutine dyn_evp2d_dump


  !=============================================================================
  subroutine dyn_evp1d_finalize()
    implicit none
    write(*,'(a21)') 'Time to die for evp1d'
    write(iu06,'(a21)') 'Time to die for evp1d'
    close (iu06)
  end subroutine dyn_evp1d_finalize

  !=============================================================================
!=============================================================================
  ! module private subroutines

  !=============================================================================
  subroutine evp_1d_nml_reader
    implicit none
    integer(kind=int_kind) :: nml_err, nmllun
    character(9), parameter :: logfilename='evp1d.log'
    namelist /evp1d_nml/ ldebug, linfo, buf1d
    ! default values
    ldebug=.true.       ! write all debug information
    linfo=.true.        ! write only info information
    buf1d=0             ! allow extra buffer for the number of 1d points that
                        ! As long as arrays are static buf1d should be 0.
                        ! may be added during the simulation
    nmllun = io_new_unit()
    open (nmllun, file='evp_1d.nml', status='old',iostat=nml_err)
    if (nml_err .ne. 0) then
      write(iu06,*) 'Read namelist evp1d_nml from file evp_1d.nml'
    endif
    read(nmllun,nml=evp1d_nml,iostat=nml_err)
    write (iu06,*) 'ldebug   = ',ldebug
    write (iu06,*) 'linfo    = ',linfo
    write (iu06,*) 'buf1d    = ',buf1d
  end subroutine evp_1d_nml_reader

  !=============================================================================
  integer(kind=int_kind) function io_new_unit ()
    implicit none
    integer(kind=int_kind) :: i
    logical    :: file_is_open
    integer(4), parameter :: io_min=7, io_max=130
    file_is_open = .true.
    i = io_min-1
    do while (file_is_open)
      i = i+1
      if (i > io_max) then
        stop 'ERROR: Could not find new I/O unit number'
      endif
      inquire(unit=i,opened=file_is_open)
    enddo
    io_new_unit = i
  end function io_new_unit

  subroutine evp1d_alloc_static_na(na0)
    implicit none
    integer(kind=int_kind), intent(in) :: na0
    integer(kind=int_kind) :: ierr

    if (linfo) then
      write(iu06,*) 'Handling 1D evp allocations on Active points'
      write(iu06,*) na0,buf1d,nx,ny
    endif

    allocate(skipTcell(1:na0+buf1d),                                           &
             skipUcell(1:na0+buf1d),                                           &
             iwidx(1:nx,1:ny),                                                 &
             unionTU(1:na0+buf1d),                                             &
             stat=ierr)
             ! FIXME add all 1D arrays
    if (ierr/=0) stop 'Error allocating'

             unionTU(:) = .false.

    allocate(indxTi(1:na0+buf1d),                                               &
             indxTj(1:na0+buf1d),                                               &
             stat=ierr)
    if (ierr/=0) stop 'Error allocating i and j indexes'

    allocate(ee(1:na0+buf1d),                                                  &
             ne(1:na0+buf1d),                                                  &
             se(1:na0+buf1d),                                                  &
             nw(1:na0+buf1d),                                                  &
             sw(1:na0+buf1d),                                                  &
             sse(1:na0+buf1d),                                                 &
             stat=ierr)
    if (ierr/=0) stop 'Error allocating indexes'

    allocate( HTE       (1:na0+buf1d   ),                                      &
              HTN       (1:na0+buf1d   ),                                      &
              HTEm1     (1:na0+buf1d   ),                                      &
              HTNm1     (1:na0+buf1d   ),                                      &
              dxt       (1:na0+buf1d   ),                                      &
              dyt       (1:na0+buf1d   ),                                      &
              strength  (1:na0+buf1d   ),                                      &
              stressp_1 (1:na0+buf1d   ),                                      &
              stressp_2 (1:na0+buf1d   ),                                      &
              stressp_3 (1:na0+buf1d   ),                                      &
              stressp_4 (1:na0+buf1d   ),                                      &
              stressm_1 (1:na0+buf1d   ),                                      &
              stressm_2 (1:na0+buf1d   ),                                      &
              stressm_3 (1:na0+buf1d   ),                                      &
              stressm_4 (1:na0+buf1d   ),                                      &
              stress12_1(1:na0+buf1d   ),                                      &
              stress12_2(1:na0+buf1d   ),                                      &
              stress12_3(1:na0+buf1d   ),                                      &
              stress12_4(1:na0+buf1d   ),stat=ierr)
    if (ierr/=0) stop 'Error allocating stress'

    allocate(cdn_ocn    (1:na0+buf1d), aiu        (1:na0+buf1d),               &
             uocn       (1:na0+buf1d), vocn       (1:na0+buf1d),               &
             waterxU    (1:na0+buf1d), wateryU    (1:na0+buf1d),               &
             forcexU    (1:na0+buf1d), forceyU    (1:na0+buf1d),               &
             umassdti   (1:na0+buf1d), fmU        (1:na0+buf1d),               &
             uarear     (1:na0+buf1d),                                         &
             strintxU   (1:na0+buf1d), strintyU   (1:na0+buf1d),               &
             Tbu        (1:na0+buf1d), Cb         (1:na0+buf1d),               &
             uvel_init  (1:na0+buf1d), vvel_init  (1:na0+buf1d),               &
             stat=ierr)
    if (ierr/=0) stop 'Error allocating in evp1d_alloc'
  end subroutine evp1d_alloc_static_na

  !=============================================================================
  subroutine evp1d_alloc_static_navel(navel0) !(na0, navel0, extrabuffer,...)
    implicit none
    integer(kind=int_kind), intent(in) :: navel0
    integer(kind=int_kind) :: ierr
    if (linfo) then
      write(iu06,*) 'Handling 1D evp allocations'
    endif

    allocate(str1(1:navel0+buf1d), str2(1:navel0+buf1d), str3(1:navel0+buf1d), &
             str4(1:navel0+buf1d), str5(1:navel0+buf1d), str6(1:navel0+buf1d), &
             str7(1:navel0+buf1d), str8(1:navel0+buf1d),                       &
             indxtij(1:navel0+buf1d),uvel(1:navel0+buf1d), vvel(1:navel0+buf1d),                       &
             stat=ierr)
    if (ierr/=0) stop 'Error allocating navel'
  end subroutine evp1d_alloc_static_navel

  subroutine evp1d_dealloc()
          ! FIXME
    implicit none
    integer(kind=int_kind) :: ierr
    if (ldebug) then
      write(iu06,*) 'Handling 1D evp de-allocations'
    endif
    ! FIXME
!    if (ierr/=0) stop 'Error allocating in evp_1d_alloc'
  end subroutine evp1d_dealloc

  subroutine update_1d(na_new, i, j)
    implicit none
    integer(kind=int_kind), intent(in)  :: na_new, i, j
    ! ensure that ALL 1D arrays get update values in index na_new, e.g. uvel(na_new) = uvel2d(i,j)
    ! FIXME 
  end subroutine update_1d

  subroutine convert_Tint2log(Tmaskint,Tmasklog)
    implicit none
    logical(kind=log_kind), intent(out) :: Tmasklog(nx,ny)
    integer(kind=int_kind), intent(in)  :: Tmaskint(nx,ny)
    integer(kind=int_kind)              :: i,j
    Tmasklog (:,:) = .false.
    do i=1,nx
      do j=1,ny
        if (Tmaskint(i,j)>0) Tmasklog(i,j)=.true.
      enddo
    enddo
  end subroutine

  subroutine calc_nActiveTU(Tmask,na0, Umask)
! Calculate number of active points with a given mask.
    implicit none
    logical(kind=log_kind), intent(in) :: Tmask(:,:)
    logical(kind=log_kind), optional, intent(in) :: Umask(:,:)
    integer(kind=int_kind), intent(out)  :: na0
    integer(kind=int_kind)              :: i,j
    na0=0
    if (present(Umask)) then
       do i=1,nx
         do j=1,ny
           if ((Tmask(i,j)) .or. (Umask(i,j))) then
             na0=na0+1
           endif
         enddo
       enddo
    else
       do i=1,nx
         do j=1,ny
           if (Tmask(i,j)) then
             na0=na0+1
           endif
         enddo
       enddo
    endif
  end subroutine calc_nActiveTU

  subroutine set_skipMe(iceTmask, iceUmask,na0)
    implicit none
    logical(kind=log_kind), intent(in) :: iceTmask(:,:), iceUmask(:,:)
    integer(kind=int_kind), intent(in) :: na0
    integer(kind=int_kind)              :: iw, i, j, niw
    skipUcell=.false.
    skipTcell=.false.
   niw=0
! first count
   do iw=1, na0
      i = indxti(iw)
      j = indxtj(iw)
      if ( iceTmask(i,j) .or. iceUmask(i,j)) niw=niw+1
      if (.not. (iceTmask(i,j))) skipTcell(iw)=.true.
      if (.not. (iceUmask(i,j))) skipUcell(iw)=.true.
   enddo
    write(iu06,*) 'number of Active points', niw

#ifdef dumjegskalnokslettes
    do i=1+nghost,nx
      do j=1+nghost,ny
        if ( iceTmask(i,j) .or. iceUmask(i,j)) niw=niw+1

      enddo
    enddo
    write(iu06,*) 'number of Active points', niw

    do i=1+nghost,nx
      do j=1+nghost,ny
        if (.not. (iceTmask(i,j))) skipTcell(iw)=.true.
      enddo
    enddo
    do i=1+nghost,nx-nghost
      do j=1+nghost,ny-nghost
        if (.not. iceUmask(i,j)) skipUcell(iw)=.true.
      enddo
    enddo
#endif
  end subroutine set_skipMe
    subroutine calc_2d_indices_init(na0, Tmask)
! All points are active. Need to find neighbors.
! This should include de selection of u points.

      implicit none

      integer(kind=int_kind), intent(in) :: na0
      logical(kind=log_kind), dimension(nx, ny), intent(in) :: &
         Tmask

      ! local variables

      integer(kind=int_kind) :: i, j, Nmaskt

      character(len=*), parameter :: subname = '(calc_2d_indices)'

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

  subroutine union(x, y, xdim, ydim, xy, nxy)

    ! Find union (xy) of two sorted integer vectors (x and y), i.e.
    ! combined values of the two vectors with no repetitions
    implicit none
    integer(4), intent(in)  :: xdim, ydim
    integer(4), intent(in)  :: x(1:xdim), y(1:ydim)
    integer(4), intent(out) :: xy(1:xdim + ydim)
    integer(4), intent(out) :: nxy
    ! local variables
    integer(kind=int_kind) :: i, j, k

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

  subroutine gather_static(L_uarear, L_dxT, L_dyT,L_Tmask, G_uarear,G_dxT,G_dyT, G_Tmask)
! In standalone  distrb_info is an integer. Not needed anyway
     use ice_communicate, only : master_task
     use ice_gather_scatter, only : gather_global_ext
     use ice_domain, only : distrb_info
     implicit none
     real(kind=dbl_kind), dimension(nx_block, ny_block, max_block), intent(in) :: &
         L_uarear,  L_dxT, L_dyT
     logical(kind=log_kind), dimension(nx_block, ny_block, max_block), intent(in) :: L_Tmask
     real(kind=dbl_kind), dimension(:, :), intent(inout) :: &
         G_uarear,G_dxT,G_dyT
     logical(kind=int_kind), dimension(nx, ny), intent(out) :: G_Tmask

! copy from distributed I_* to G_*
     call gather_global_ext(G_uarear,     L_uarear,     master_task, distrb_info    )
     call gather_global_ext(G_dxT,        L_dxT,        master_task, distrb_info    )
     call gather_global_ext(G_dyT,        L_dyT,        master_task, distrb_info    )
     call gather_global_ext(G_Tmask,      L_Tmask,        master_task, distrb_info    )
  end subroutine gather_static

  subroutine gather_dyn(L_stressp_1 , L_stressp_2 , L_stressp_3, L_stressp_4,     &
                        L_stressm_1 , L_stressm_2 , L_stressm_3, L_stressm_4,     &
                        L_stress12_1, L_stress12_2, L_stress12_3,L_stress12_4,    &
                        L_cdn_ocn   , L_aiu       , L_uocn     , L_vocn     ,     &
                        L_waterxU   , L_wateryU   , L_forcexU  , L_forceyU  ,     &
                        L_umassdti  , L_fmU       , L_strintxU , L_strintyU ,     &
                        L_Tbu       , L_Cb        , L_uvel     , L_vvel     ,     &
                        L_icetmask  , L_iceUmask  ,                               &
                        G_stressp_1 , G_stressp_2 , G_stressp_3, G_stressp_4,     &
                        G_stressm_1 , G_stressm_2 , G_stressm_3, G_stressm_4,     &
                        G_stress12_1, G_stress12_2, G_stress12_3,G_stress12_4,    &
                        G_cdn_ocn   , G_aiu       , G_uocn     , G_vocn     ,     &
                        G_waterxU   , G_wateryU   , G_forcexU  , G_forceyU  ,     &
                        G_umassdti  , G_fmU       , G_strintxU , G_strintyU ,     &
                        G_Tbu       , G_Cb        , G_uvel     , G_vvel     ,     &
                        G_iceTmask,  G_iceUmask)

     use ice_communicate, only : master_task
     use ice_gather_scatter, only : gather_global_ext
     use ice_domain, only : distrb_info
     implicit none
     real(kind=dbl_kind), dimension(nx_block, ny_block, max_block), intent(in) :: &
                            L_stressp_1 , L_stressp_2 , L_stressp_3, L_stressp_4,     &
                            L_stressm_1 , L_stressm_2 , L_stressm_3, L_stressm_4,     &
                            L_stress12_1, L_stress12_2, L_stress12_3,L_stress12_4,    &
                            L_cdn_ocn   , L_aiu       , L_uocn     , L_vocn     ,     &
                            L_waterxU   , L_wateryU   , L_forcexU  , L_forceyU  ,     &
                            L_umassdti  , L_fmU       , L_strintxU , L_strintyU ,     &
                            L_Tbu       , L_Cb        , L_uvel     , L_vvel

     logical(kind=log_kind), dimension(nx_block, ny_block, max_block), intent(in) ::  &
                            L_iceUmask

     integer(kind=int_kind), dimension(nx_block, ny_block, max_block), intent(in) ::  &
                           L_iceTmask

     real(kind=dbl_kind), dimension(nx, ny), intent(out) :: &
                            G_stressp_1 , G_stressp_2 , G_stressp_3, G_stressp_4,     &
                            G_stressm_1 , G_stressm_2 , G_stressm_3, G_stressm_4,     &
                            G_stress12_1, G_stress12_2, G_stress12_3,G_stress12_4,    &
                            G_cdn_ocn   , G_aiu       , G_uocn     , G_vocn     ,     &
                            G_waterxU   , G_wateryU   , G_forcexU  , G_forceyU  ,     &
                            G_umassdti  , G_fmU       , G_strintxU , G_strintyU ,     &
                            G_Tbu       , G_Cb        , G_uvel     , G_vvel
     integer(kind=int_kind), dimension(nx, ny), intent(out) :: G_iceTmask
     logical(kind=log_kind), dimension(nx, ny), intent(out) :: G_iceUmask

! copy from distributed I_* to G_*
    call gather_global_ext(G_stressp_1,     L_stressp_1,     master_task, distrb_info)
    call gather_global_ext(G_stressp_2,     L_stressp_2,     master_task, distrb_info)
    call gather_global_ext(G_stressp_3,     L_stressp_3,     master_task, distrb_info)
    call gather_global_ext(G_stressp_4,     L_stressp_4,     master_task, distrb_info)

    call gather_global_ext(G_stressm_1,     L_stressm_1,     master_task, distrb_info)
    call gather_global_ext(G_stressm_2,     L_stressm_2,     master_task, distrb_info)
    call gather_global_ext(G_stressm_3,     L_stressm_3,     master_task, distrb_info)
    call gather_global_ext(G_stressm_4,     L_stressm_4,     master_task, distrb_info)

    call gather_global_ext(G_stress12_1,     L_stress12_1,     master_task, distrb_info)
    call gather_global_ext(G_stress12_2,     L_stress12_2,     master_task, distrb_info)
    call gather_global_ext(G_stress12_3,     L_stress12_3,     master_task, distrb_info)
    call gather_global_ext(G_stress12_4,     L_stress12_4,     master_task, distrb_info)

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
    call gather_global_ext(G_strintxU  ,     L_strintxU  ,     master_task, distrb_info)
    call gather_global_ext(G_strintyU  ,     L_strintyU  ,     master_task, distrb_info)


    call gather_global_ext(G_Tbu       ,     L_Tbu       ,     master_task, distrb_info)
    call gather_global_ext(G_Cb        ,     L_Cb        ,     master_task, distrb_info)
    call gather_global_ext(G_uvel      ,     L_uvel      ,     master_task, distrb_info)
    call gather_global_ext(G_vvel      ,     L_vvel      ,     master_task, distrb_info)
    call gather_global_ext(G_iceTmask  ,     L_iceTmask  ,     master_task, distrb_info)
    call gather_global_ext(G_iceUmask  ,     L_iceUmask  ,     master_task, distrb_info)
 end subroutine gather_dyn

  subroutine convert_2d_1d_init(na0, G_HTE, G_HTN, &
      G_uarear,  G_dxT, G_dyT)

      implicit none

      integer(kind=int_kind), intent(in) ::  na0
      real (kind=dbl_kind), dimension(:, :), intent(in) :: &
        G_HTE, G_HTN, G_uarear, G_dxT, G_dyT

      ! local variables

      integer(kind=int_kind) :: iw, lo, up, j, i
      integer(kind=int_kind), dimension(1:na0) :: Iin, Iee, Ine, Ise, &
         Inw, Isw, Isse
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
         indxtij(iw) = Iin(iw)
      end do
      ! sorted additional points
      call setdiff(util2, Iin, navel, na0, util1, j)
      do iw = na0 + 1, navel
      ! FIXME = util2? or maybe just na0+1
         indxtij(iw) = util1(iw - na0)
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
      do iw = lo, up
         ! get 2D indices
         i = indxti(iw)
         j = indxtj(iw)
         ! map
         uarear(iw)     = G_uarear(i, j)
         dxT(iw)        = G_dxT(i, j)
         dyT(iw)        = G_dyT(i, j)
         HTE(iw)        = G_HTE(i, j)
         HTN(iw)        = G_HTN(i, j)
         HTEm1(iw)      = G_HTE(i - 1, j)
         HTNm1(iw)      = G_HTN(i, j - 1)
      end do
   end subroutine convert_2d_1d_init

   subroutine convert_2d_1d_dyn(na0, &
              G_stressp_1 , G_stressp_2 , G_stressp_3, G_stressp_4 ,  &
              G_stressm_1 , G_stressm_2 , G_stressm_3, G_stressm_4 ,  &
              G_stress12_1, G_stress12_2, G_stress12_3,G_stress12_4,  &
              G_cdn_ocn   , G_aiu       , G_uocn     , G_vocn      ,  &
              G_waterxU   , G_wateryU   , G_forcexU  , G_forceyU   ,  &
              G_umassdti  , G_fmU       , G_strintxU , G_strintyU  ,  &
              G_Tbu       , G_Cb        , G_uvel     , G_vvel      )


      implicit none

      integer(kind=int_kind), intent(in) ::  na0
      real(kind=dbl_kind), dimension(:, :), intent(in) :: &
                            G_stressp_1 , G_stressp_2 , G_stressp_3, G_stressp_4,     &
                            G_stressm_1 , G_stressm_2 , G_stressm_3, G_stressm_4,     &
                            G_stress12_1, G_stress12_2, G_stress12_3,G_stress12_4,    &
                            G_cdn_ocn   , G_aiu       , G_uocn     , G_vocn     ,     &
                            G_waterxU   , G_wateryU   , G_forcexU  , G_forceyU  ,     &
                            G_umassdti  , G_fmU       , G_strintxU , G_strintyU ,     &
                            G_Tbu       , G_Cb        , G_uvel     , G_vvel
      integer(kind=int_kind) ::  lo, up, iw, i, j

     lo=1
      up=na0
      do iw = lo, up
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
         strintxU(iw)   = G_strintxU(i, j)
         strintyU(iw)   = G_strintyU (i, j)
         Tbu(iw)        = G_Tbu(i, j)
         Cb(iw)         = G_Cb(i, j)
         uvel(iw)       = G_uvel(i,j)
         vvel(iw)       = G_vvel(i,j)
     end do
         uvel(na0+1:navel)=c0
         vvel(na0+1:navel)=c0
   end subroutine convert_2d_1d_dyn

subroutine convert_1d_2d_dyn(na0, &
              G_stressp_1 , G_stressp_2 , G_stressp_3, G_stressp_4,     &
              G_stressm_1 , G_stressm_2 , G_stressm_3, G_stressm_4,     &
              G_stress12_1, G_stress12_2, G_stress12_3,G_stress12_4,    &
              G_cdn_ocn   , G_aiu       , G_uocn     , G_vocn     ,     &
              G_waterxU   , G_wateryU   , G_forcexU  , G_forceyU  ,     &
              G_umassdti  , G_fmU       , G_strintxU , G_strintyU ,     &
              G_Tbu       , G_Cb        , G_uvel     , G_vvel     )

      implicit none

      integer(kind=int_kind), intent(in) ::  na0
      real(kind=dbl_kind), dimension(:, :), intent(inout) :: &
                            G_stressp_1 , G_stressp_2 , G_stressp_3, G_stressp_4,     &
                            G_stressm_1 , G_stressm_2 , G_stressm_3, G_stressm_4,     &
                            G_stress12_1, G_stress12_2, G_stress12_3,G_stress12_4,    &
                            G_cdn_ocn   , G_aiu       , G_uocn     , G_vocn     ,     &
                            G_waterxU   , G_wateryU   , G_forcexU  , G_forceyU  ,     &
                            G_umassdti  , G_fmU       , G_strintxU , G_strintyU ,     &
                            G_Tbu       , G_Cb        , G_uvel     , G_vvel
      integer(kind=int_kind) ::  lo, up, iw, i, j


      lo=1
      up=na0
      do iw = lo, up
         ! get 2D indices
         i = indxti(iw)
         j = indxtj(iw)
         ! map
         G_stressp_1 (i,j)  = stressp_1(iw)
         G_stressp_2 (i,j)  = stressp_2(iw)
         G_stressp_3 (i,j)  = stressp_3(iw)
         G_stressp_4 (i,j)  = stressp_4(iw)
         G_stressm_1 (i,j)  = stressm_1(iw)
         G_stressm_2 (i,j)  = stressm_2(iw)
         G_stressm_3 (i,j)  = stressm_3(iw)
         G_stressm_4 (i,j)  = stressm_4(iw)
         G_stress12_1(i,j) = stress12_1(iw)
         G_stress12_2(i,j) = stress12_2(iw)
         G_stress12_3(i,j) = stress12_3(iw)
         G_stress12_4(i,j) = stress12_4(iw)

         G_cdn_ocn(i,j)    = cdn_ocn(iw)
         G_aiu(i,j)        = aiu(iw)
         G_uocn(i,j)       = uocn(iw)
         G_vocn(i,j)       = vocn(iw)
         G_waterxU(i,j)    = waterxU(iw)
         G_wateryU(i,j)    = wateryU(iw)
         G_forcexU(i,j)    = forcexU(iw)
         G_forceyU(i,j)    = forceyU(iw)
         G_umassdti(i,j)   = umassdti(iw)
         G_fmU(i,j)        = fmU(iw)
         G_strintxU(i,j)   = strintxU(iw)
         G_strintyU(i,j)   = strintyU (iw)
         G_Tbu(i,j)        = Tbu(iw)
         G_Cb(i,j)         = Cb(iw)
         G_uvel(i,j)       = uvel(iw)
         G_vvel(i,j)       = vvel(iw)

     end do
   end subroutine convert_1d_2d_dyn

   !=======================================================================

   subroutine setdiff(x, y,  lvecx, lvecy,xy, nxy)
      ! Find element (xy) of two sorted integer vectors (x and y) that
      ! are in x, but not in y, or in y, but not in x

      implicit none

      integer(kind=int_kind), intent(in) :: x(1:lvecx), y(1:lvecy), lvecx,lvecy
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
!            call abort_ice(subname &
!               // ': ERROR: conditions not met')
         end if
      end do

   end subroutine findXinY

   subroutine calc_navel(na0, navel0)
      ! Calculate number of active points, including halo points

      implicit none

      integer(kind=int_kind), intent(in) :: na0
      integer(kind=int_kind), intent(out) :: navel0

      ! local variables

      integer(kind=int_kind) :: iw, i, j
      integer(kind=int_kind), dimension(1:na0) :: Iin, Iee, Ine, Ise, &
         Inw, Isw, Isse, indi, indj
      integer(kind=int_kind), dimension(1:7 * na0) :: util1, util2

      character(len=*), parameter :: subname = '(calc_navel)'

      ! calculate additional 1D indices used for finite differences
      do iw = 1, na0
         ! get 2D indices
         i = indi(iw)
         j = indj(iw)
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
      call union(Iin,   Iee,  na0, na0, util1, i    )
      call union(util1, Ine,  i,  na0, util2, j    )
      call union(util2, Ise,  j,  na0, util1, i    )
      call union(util1, Inw,  i,  na0, util2, j    )
      call union(util2, Isw,  j,  na0, util1, i    )
      call union(util1, Isse, i,  na0, util2, navel0)

   end subroutine calc_navel

end module ice_dyn_evp1d

