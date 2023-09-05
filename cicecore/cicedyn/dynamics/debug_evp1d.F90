!===============================================================================
module debug_evp1d
  !- modules -------------------------------------------------------------------
  use ice_kinds_mod
  !- directives ----------------------------------------------------------------
  implicit none
  integer(kind=int_kind) :: nts, nfield
  integer(kind=int_kind), allocatable :: lts(:)
  character(len=10), allocatable :: lfields(:)
  interface dumpfield
    module procedure dumpfield_double
    module procedure dumpfield_integer
    module procedure dumpfield_logical
  end interface
  public :: dumpfield, dumpall, dump_init, dumpall3d
  contains

  subroutine dump_init(iu06)
  implicit none
  integer(kind=int_kind) iu06
  integer(kind=int_kind) nmllun, nml_err
  namelist /dumpspec/ nts, nfield
  namelist /dumpts/ lts, lfields
  nmllun = 10
  open (nmllun, file='evp1d_debug.nml', status='old',iostat=nml_err)
  write (iu06,*) 'nts = ',nts
  if (nml_err .ne. 0) then
    write(iu06,*) 'Error open file evp1d_debug.nml'
  endif
  read(nmllun,nml=dumpspec,iostat=nml_err)
  if (nml_err .ne. 0) then
    write(iu06,*) 'Read namelist dumpspec from file evp1d_debug.nml'
  endif
  write (iu06,*) 'nts   = ',nts
  write (iu06,*) 'nfield   = ',nfield
  allocate(lts(nts))
  lts=1
  write (iu06,*) 'list   = ',lts
  read(nmllun,nml=dumpts,iostat=nml_err)
  if (nml_err .ne. 0) then
    write(iu06,*) 'Read namelist dumpts from file evp1d_debug.nml'
  endif
  write (iu06,*) 'list of ts = ',lts
  end subroutine dump_init

 subroutine dumpall(myname, ts, nx, ny, iu06,                               &
                       G_stressp_1 , G_stressp_2 , G_stressp_3, G_stressp_4, &
                       G_stressm_1 , G_stressm_2 , G_stressm_3, G_stressm_4, &
                       G_stress12_1, G_stress12_2, G_stress12_3,G_stress12_4,&
                       G_cdn_ocn   , G_aiu       , G_uocn     , G_vocn     , &
                       G_waterxU   , G_wateryU   , G_forcexU  , G_forceyU  , &
                       G_umassdti  , G_fmU       , G_strintxU , G_strintyU , &
                       G_Tbu       , G_uvel     , G_vvel     )
    implicit none
    integer(kind=int_kind), intent(in) :: ts, nx, ny, iu06
    real(kind=dbl_kind), dimension(:,:),intent(in) :: &
                       G_stressp_1 , G_stressp_2 , G_stressp_3, G_stressp_4, &
                       G_stressm_1 , G_stressm_2 , G_stressm_3, G_stressm_4, &
                       G_stress12_1, G_stress12_2, G_stress12_3,G_stress12_4,&
                       G_cdn_ocn   , G_aiu       , G_uocn     , G_vocn     , &
                       G_waterxU   , G_wateryU   , G_forcexU  , G_forceyU  , &
                       G_umassdti  , G_fmU       , G_strintxU , G_strintyU , &
                       G_Tbu       , G_uvel     , G_vvel

    character(len=10), intent(in) :: myname
    character(len=200) :: binfile
    character(4) :: ctmp
    integer(kind=int_kind) :: j
    integer(kind=int_kind) :: lun, ios
    logical(kind=log_kind) :: dumpthistime
    dumpthistime=.false.
    do j=1,nts
       if (ts==lts(j)) then
          write (iu06,*) 'Time to dump ',ts, j
          dumpthistime=.true.
          write (iu06,*) 'Dumping all field from 2D and 1D world into two 2D files'
          write(ctmp,'(i4.4)') ts
          binfile = trim(trim(myname)//'.double.ts-'//ctmp//'.bin')
          write (iu06,*) 'Opening file ', trim(binfile)
          lun=709
          open(lun,file=trim(binfile), form='unformatted', access='stream', action='write', status='replace', iostat=ios)
          if (ios .ne. 0) then
             stop ('Failed open file')
          endif
          write(lun,iostat=ios) &
                G_stressp_1 , G_stressp_2 , G_stressp_3, G_stressp_4,     &
                G_stressm_1 , G_stressm_2 , G_stressm_3, G_stressm_4,     &
                G_stress12_1, G_stress12_2, G_stress12_3,G_stress12_4,    &
                G_cdn_ocn   , G_aiu       , G_uocn     , G_vocn     ,     &
                G_waterxU   , G_wateryU   , G_forcexU  , G_forceyU  ,     &
                G_umassdti  , G_fmU       , G_strintxU , G_strintyU ,     &
                G_Tbu       , G_uvel     , G_vvel
          if (ios .ne. 0) then
             stop ('Failed write file')
          endif
          close(lun)
       endif
    enddo
    if (.not. dumpthistime) then
           write(iu06,*) 'No file dump at ts  ', ts
    endif
  end subroutine dumpall

 subroutine dumpfield_double(myname, ts, f)
    implicit none
    real(kind=dbl_kind), intent(in) :: f(:)
    integer(kind=int_kind), intent(in) :: ts
    character(len=10), intent(in) :: myname
    integer(kind=int_kind) :: lun, ios
    character(len=17) :: binfile
    character(4) :: ctmp
    write(ctmp,'(i4.4)') ts
    binfile = trim(trim(myname)//'.double.ts-'//ctmp//'.bin')
    write (*,*) 'Opening file ', trim(binfile)
    lun = 710
    open(lun,file=binfile, form='unformatted', access='stream', action='write', status='replace', iostat=ios)
    if (ios .ne. 0) then
      stop ('Failed open file')
    endif
    write(lun,iostat=ios) f
    if (ios .ne. 0) then
      stop ('Failed write file')
    endif
    close(lun)
  end subroutine dumpfield_double
  subroutine dumpfield_integer(myname, ts, f)
    implicit none
    integer(kind=int_kind), intent(in) :: f(:)
    integer(kind=int_kind), intent(in) :: ts
    character(len=10), intent(in) :: myname
    character(len=17) :: binfile
    integer(kind=int_kind) :: lun, ios
    character(4) :: ctmp
    write(ctmp,'(i4.4)') ts
    binfile = trim(trim(myname)//'.integer.ts-'//ctmp//'.bin')
    write (*,*) 'Opening file ', trim(binfile)
    ! FIXME
  end subroutine dumpfield_integer
  subroutine dumpfield_logical(myname, ts, f)
    implicit none
    logical(kind=log_kind), intent(in) :: f(:)
    integer(kind=int_kind), intent(in) :: ts
    character(len=10), intent(in) :: myname
    character(len=17) :: binfile
    character(4) :: ctmp
    write(ctmp,'(i4.4)') ts
    binfile = trim(trim(myname)//'.integer.ts-'//ctmp//'.bin')
    write (*,*) 'Opening file ', trim(binfile)
  end subroutine dumpfield_logical

 subroutine dumpall3d(myname, ts,                               &
                       G_stressp_1 , G_stressp_2 , G_stressp_3, G_stressp_4, &
                       G_stressm_1 , G_stressm_2 , G_stressm_3, G_stressm_4, &
                       G_stress12_1, G_stress12_2, G_stress12_3,G_stress12_4,&
                       G_cdn_ocn   , G_aiu       , G_uocn     , G_vocn     , &
                       G_waterxU   , G_wateryU   , G_forcexU  , G_forceyU  , &
                       G_umassdti  , G_fmU       , G_strintxU , G_strintyU , &
                       G_Tbu       , G_uvel     , G_vvel     )
    implicit none
    integer(kind=int_kind), intent(in) :: ts
    real(kind=dbl_kind), dimension(:,:,:),intent(in) :: &
                       G_stressp_1 , G_stressp_2 , G_stressp_3, G_stressp_4, &
                       G_stressm_1 , G_stressm_2 , G_stressm_3, G_stressm_4, &
                       G_stress12_1, G_stress12_2, G_stress12_3,G_stress12_4,&
                       G_cdn_ocn   , G_aiu       , G_uocn     , G_vocn     , &
                       G_waterxU   , G_wateryU   , G_forcexU  , G_forceyU  , &
                       G_umassdti  , G_fmU       , G_strintxU , G_strintyU , &
                       G_Tbu       , G_uvel     , G_vvel

    character(len=10), intent(in) :: myname
    character(len=200) :: binfile
    character(4) :: ctmp
    integer(kind=int_kind) :: lun, ios
    write (999,*) 'Time to dump ',ts
    write (999,*) 'Dumping all field from 2D and 1D world into two 2D files'
    write(ctmp,'(i4.4)') ts
    binfile = trim(trim(myname)//'.double.ts-'//ctmp//'.bin')
    write (999,*) 'Opening file ', trim(binfile)
    lun=709
    open(lun,file=trim(binfile), form='unformatted', access='stream', action='write', status='replace', iostat=ios)
    if (ios .ne. 0) then
        stop ('Failed open file')
    endif
    write(lun,iostat=ios) &
          G_stressp_1 , G_stressp_2 , G_stressp_3, G_stressp_4,     &
          G_stressm_1 , G_stressm_2 , G_stressm_3, G_stressm_4,     &
          G_stress12_1, G_stress12_2, G_stress12_3,G_stress12_4,    &
          G_cdn_ocn   , G_aiu       , G_uocn     , G_vocn     ,     &
          G_waterxU   , G_wateryU   , G_forcexU  , G_forceyU  ,     &
          G_umassdti  , G_fmU       , G_strintxU , G_strintyU ,     &
          G_Tbu       , G_uvel     , G_vvel
    if (ios .ne. 0) then
        stop ('Failed write file')
    endif
    close(lun)
  end subroutine dumpall3d

end module debug_evp1d

