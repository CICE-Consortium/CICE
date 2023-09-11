!===============================================================================
! Copyright (C) 2023, Intel Corporation
! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
! 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
! 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
! 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! Modified September 2023, Till Rasmussen DMI
! This version is only valid for v2. 
! It assumes that all variables are passed as function arguments
!===============================================================================
module numa
  use ice_kinds_mod
  use myomp, only    : domp_get_domain
  !- directives ----------------------------------------------------------------
  implicit none
  private
  !- interfaces ----------------------------------------------------------------
  interface numainit
    module procedure numainit_real
    module procedure numainit_int
    module procedure numainit_logical
  end interface 
  public  ::  numainit, numainit_all, numareinit
  private ::  numainit_real, numainit_int, numainit_logical
  contains
  subroutine numainit_real(l,u,a,tmpVa)
    real   (kind=dbl_kind), intent(in),    dimension(:), contiguous :: tmpVa
    real   (kind=dbl_kind), intent(inout), dimension(:), contiguous :: a
    integer(kind=int_kind), intent(in)                              :: l,u
    integer(kind=int_kind)                                          :: lo,up
    call domp_get_domain(l,u,lo,up)
    a(l:u)=tmpVa(l:u)
  end subroutine numainit_real
  subroutine numainit_int(l,u,i,tmpVi)
    integer(kind=int_kind), intent(in),    dimension(:), contiguous :: tmpVi
    integer(kind=int_kind), intent(inout), dimension(:), contiguous :: i
    integer(kind=int_kind), intent(in)                              :: l,u
    integer(kind=int_kind)                                          :: lo,up
    call domp_get_domain(l,u,lo,up)
    i(l:u)=tmpVi(l:u)
  end subroutine numainit_int
  subroutine numainit_logical(l,u,ll,tmpVll)
    logical(kind=log_kind), intent(in),    dimension(:), contiguous :: tmpVll
    logical(kind=log_kind), intent(inout), dimension(:), contiguous :: ll
    integer(kind=int_kind), intent(in)                              :: l,u
    integer(kind=int_kind)                                          :: lo,up
    call domp_get_domain(l,u,lo,up)
    ll(l:u)=tmpVll(l:u)
  end subroutine numainit_logical
  subroutine numainit_all(l,u,uu)
    !- modules -----------------------------------------------------------------
    use myomp,         only : domp_get_domain
    use ice_constants, only : c0, c1
    use vars,          only : uvel, vvel, dxT, dyT, strength,                  &
                              stressp_1, stressp_2, stressp_3, stressp_4,      &
                              stressm_1, stressm_2, stressm_3, stressm_4,      &
                              stress12_1, stress12_2, stress12_3, stress12_4,  &
                              cdn_ocn, aiu, uocn, vocn, waterxU, wateryU,      &
                              forcexU, forceyU, umassdti, fmU, uarear, Tbu,    &
                              strintxU, strintyU, uvel_init, vvel_init, Cb,    &
                              ee,ne,se,nw,sw,sse, skipTcell1d,  skipUcell1d,   &
                              str1, str2, str3, str4, str5, str6, str7, str8,  &
                              HTE1d,HTN1d, HTE1dm1,HTN1dm1
    implicit none
    integer(kind=int_kind),intent(in) :: l,u,uu
    integer(kind=int_kind) :: lo,up
    call domp_get_domain(l,u,lo,up)
    if ((lo > 0) .and. (up >= lo)) then
       skipTcell1d(lo:up)=.false.
       skipUcell1d(lo:up)=.false.
       ee(lo:up)=0
       ne(lo:up)=0
       se(lo:up)=0
       nw(lo:up)=0
       sw(lo:up)=0
       sse(lo:up)=0
       aiu(lo:up)=c0
       Cb(lo:up)=c0
       cdn_ocn(lo:up)=c0
       dxt(lo:up)=c0
       dyt(lo:up)=c0
       fmU(lo:up)=c0
       forcexU(lo:up)=c0
       forceyU(lo:up)=c0
       HTE1d(lo:up)=c0
       HTE1dm1(lo:up)=c0
       HTN1d(lo:up)=c0
       HTN1dm1(lo:up)=c0
       str1(lo:up)=c0
       str2(lo:up)=c0
       str3(lo:up)=c0
       str4(lo:up)=c0
       str5(lo:up)=c0
       str6(lo:up)=c0
       str7(lo:up)=c0
       str8(lo:up)=c0
       strength(lo:up)= c0
       stress12_1(lo:up)=c0
       stress12_2(lo:up)=c0
       stress12_3(lo:up)=c0
       stress12_4(lo:up)=c0
       stressm_1(lo:up)=c0
       stressm_2(lo:up)=c0
       stressm_3(lo:up)=c0
       stressm_4(lo:up)=c0
       stressp_1(lo:up)=c0
       stressp_2(lo:up)=c0
       stressp_3(lo:up)=c0
       stressp_4(lo:up)=c0
       strintxU(lo:up)= c0
       strintyU(lo:up)= c0
       Tbu(lo:up)=c0
       uarear(lo:up)=c0
       umassdti(lo:up)=c0
       uocn(lo:up)=c0
       uvel_init(lo:up)=c0
       uvel(lo:up)=c0
       vocn(lo:up)=c0
       vvel_init(lo:up)=c0
       vvel(lo:up)=c0
       waterxU(lo:up)=c0
       wateryU(lo:up)=c0
    endif
    call domp_get_domain(u+1,uu,lo,up)
    if ((lo > 0) .and. (up >= lo)) then
      uvel(lo:up)=c0
      vvel(lo:up)=c0
    endif
  end subroutine numainit_all
  subroutine numareinit()
    use vars,          only : uvel, vvel,                                      &
                              str1, str2, str3, str4, str5, str6, str7, str8,  &
                              dxT, dyT, strength,                              &
                              stressp_1, stressp_2, stressp_3, stressp_4,      &
                              stressm_1, stressm_2, stressm_3, stressm_4,      &
                              stress12_1, stress12_2, stress12_3, stress12_4,  &
                              cdn_ocn, aiu, uocn, vocn, waterxU, wateryU,      &
                              forcexU, forceyU, umassdti, fmU, uarear, Tbu,    &
                              strintxU, strintyU, uvel_init, vvel_init, Cb,    &
                              ee,ne,se,nw,sw,sse, skipTcell1d,  skipUcell1d,   &
                              HTE1d,HTN1d, HTE1dm1,HTN1dm1
    implicit none
    real   (kind=dbl_kind), dimension(:), allocatable :: tmpVa
    integer(kind=int_kind), dimension(:), allocatable :: tmpVi
    logical(kind=log_kind), dimension(:), allocatable :: tmpVll
    integer(kind=int_kind) :: nanew, navelnew
    integer(kind=int_kind) :: ierr
    navelnew = size(str1)
    ! navel: 
    ! real: uvel, vvel, str1...str8
    allocate(tmpVa(navelnew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'

    tmpVa=str1
    deallocate(str1)
    allocate(str1(navelnew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,navelnew,str1,tmpVa)
!$OMP END PARALLEL

    tmpVa=str2
    deallocate(str2)
    allocate(str2(navelnew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,navelnew,str2,tmpVa)
!$OMP END PARALLEL

    tmpVa=str3
    deallocate(str3)
    allocate(str3(navelnew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,navelnew,str3,tmpVa)
!$OMP END PARALLEL

    tmpVa=str4
    deallocate(str4)
    allocate(str4(navelnew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,navelnew,str4,tmpVa)
!$OMP END PARALLEL

    tmpVa=str5
    deallocate(str5)
    allocate(str5(navelnew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,navelnew,str5,tmpVa)
!$OMP END PARALLEL

    tmpVa=str6
    deallocate(str6)
    allocate(str6(navelnew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,navelnew,str6,tmpVa)
!$OMP END PARALLEL

    tmpVa=str7
    deallocate(str7)
    allocate(str7(navelnew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,navelnew,str7,tmpVa)
!$OMP END PARALLEL

    tmpVa=str8
    deallocate(str8)
    allocate(str8(navelnew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,navelnew,str8,tmpVa)
!$OMP END PARALLEL

    tmpVa=uvel
    deallocate(uvel)
    allocate(uvel(navelnew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,navelnew,uvel,tmpVa)
!$OMP END PARALLEL

    tmpVa=vvel
    deallocate(vvel)
    allocate(vvel(navelnew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,navelnew,vvel,tmpVa)
!$OMP END PARALLEL

    ! na:
    ! real: dxT, dyT, strength,
    !       stressp_1, stressp_2, stressp_3, stressp_4,      &
    !       stressm_1, stressm_2, stressm_3, stressm_4,      &
    !       stress12_1, stress12_2, stress12_3, stress12_4,  &
    !       cdn_ocn, aiu, uocn, vocn, waterxU, wateryU,      &
    !       forcexU, forceyU, umassdti, fmU, uarear, Tbu,    &
    !       strintxU, strintyU, uvel_init, vvel_init, Cb,    &
    !       HTE1d,HTN1d, HTE1dm1,HTN1dm1
    ! int:  ee,ne,se,nw,sw,sse
    ! logical: skipTcell1d,  skipUcell1d,   &
    deallocate(tmpVa)
    nanew = size(HTE1d)
    allocate(tmpVa(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'

    tmpVa=dxT
    deallocate(dxT)
    allocate(dxT(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,dxT,tmpVa)
!$OMP END PARALLEL

    tmpVa=dyT
    deallocate(dyT)
    allocate(dyT(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,dyT,tmpVa)
!$OMP END PARALLEL

    tmpVa=strength
    deallocate(strength)
    allocate(strength(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,strength,tmpVa)
!$OMP END PARALLEL

    tmpVa=stressp_1
    deallocate(stressp_1)
    allocate(stressp_1(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,stressp_1,tmpVa)
!$OMP END PARALLEL

    tmpVa=stressp_2
    deallocate(stressp_2)
    allocate(stressp_2(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,stressp_2,tmpVa)
!$OMP END PARALLEL

    tmpVa=stressp_3
    deallocate(stressp_3)
    allocate(stressp_3(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,stressp_3,tmpVa)
!$OMP END PARALLEL

    tmpVa=stressp_4
    deallocate(stressp_4)
    allocate(stressp_4(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,stressp_4,tmpVa)
!$OMP END PARALLEL

    tmpVa=stressm_1
    deallocate(stressm_1)
    allocate(stressm_1(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,stressm_1,tmpVa)
!$OMP END PARALLEL

    tmpVa=stressm_2
    deallocate(stressm_2)
    allocate(stressm_2(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,stressm_2,tmpVa)
!$OMP END PARALLEL

    tmpVa=stressm_3
    deallocate(stressm_3)
    allocate(stressm_3(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,stressm_3,tmpVa)
!$OMP END PARALLEL

    tmpVa=stressm_4
    deallocate(stressm_4)
    allocate(stressm_4(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,stressm_4,tmpVa)
!$OMP END PARALLEL

    tmpVa=stress12_1
    deallocate(stress12_1)
    allocate(stress12_1(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,stress12_1,tmpVa)
!$OMP END PARALLEL

    tmpVa=stress12_2
    deallocate(stress12_2)
    allocate(stress12_2(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,stress12_2,tmpVa)
!$OMP END PARALLEL

    tmpVa=stress12_3
    deallocate(stress12_3)
    allocate(stress12_3(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,stress12_3,tmpVa)
!$OMP END PARALLEL

    tmpVa=stress12_4
    deallocate(stress12_4)
    allocate(stress12_4(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,stress12_4,tmpVa)
!$OMP END PARALLEL

    tmpVa=cdn_ocn
    deallocate(cdn_ocn)
    allocate(cdn_ocn(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,cdn_ocn,tmpVa)
!$OMP END PARALLEL

    tmpVa=aiu
    deallocate(aiu)
    allocate(aiu(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,aiu,tmpVa)
!$OMP END PARALLEL

    tmpVa=uocn
    deallocate(uocn)
    allocate(uocn(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,uocn,tmpVa)
!$OMP END PARALLEL

    tmpVa=vocn
    deallocate(vocn)
    allocate(vocn(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,vocn,tmpVa)
!$OMP END PARALLEL

    tmpVa=waterxU
    deallocate(waterxU)
    allocate(waterxU(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,waterxU,tmpVa)
!$OMP END PARALLEL

    tmpVa=wateryU
    deallocate(wateryU)
    allocate(wateryU(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,wateryU,tmpVa)
!$OMP END PARALLEL

    tmpVa=forcexU
    deallocate(forcexU)
    allocate(forcexU(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,forcexU,tmpVa)
!$OMP END PARALLEL

    tmpVa=forceyU
    deallocate(forceyU)
    allocate(forceyU(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,forceyU,tmpVa)
!$OMP END PARALLEL

    tmpVa=umassdti
    deallocate(umassdti)
    allocate(umassdti(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,umassdti,tmpVa)
!$OMP END PARALLEL

    tmpVa=fmU
    deallocate(fmU)
    allocate(fmU(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,fmU,tmpVa)
!$OMP END PARALLEL

    tmpVa=uarear
    deallocate(uarear)
    allocate(uarear(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,uarear,tmpVa)
!$OMP END PARALLEL

    tmpVa=Tbu
    deallocate(Tbu)
    allocate(Tbu(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,Tbu,tmpVa)
!$OMP END PARALLEL

    tmpVa=strintxU
    deallocate(strintxU)
    allocate(strintxU(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,strintxU,tmpVa)
!$OMP END PARALLEL

    tmpVa=strintyU
    deallocate(strintyU)
    allocate(strintyU(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,strintyU,tmpVa)
!$OMP END PARALLEL

    tmpVa=uvel_init
    deallocate(uvel_init)
    allocate(uvel_init(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,uvel_init,tmpVa)
!$OMP END PARALLEL

    tmpVa=vvel_init
    deallocate(vvel_init)
    allocate(vvel_init(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,vvel_init,tmpVa)
!$OMP END PARALLEL

    tmpVa=Cb
    deallocate(Cb)
    allocate(Cb(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,Cb,tmpVa)
!$OMP END PARALLEL

    tmpVa=HTE1d
    deallocate(HTE1d)
    allocate(HTE1d(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,HTE1d,tmpVa)
!$OMP END PARALLEL

    tmpVa=HTN1d
    deallocate(HTN1d)
    allocate(HTN1d(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,HTN1d,tmpVa)
!$OMP END PARALLEL

    tmpVa=HTE1dm1
    deallocate(HTE1dm1)
    allocate(HTE1dm1(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,HTE1dm1,tmpVa)
!$OMP END PARALLEL

    tmpVa=HTN1dm1
    deallocate(HTN1dm1)
    allocate(HTN1dm1(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,HTN1dm1,tmpVa)
!$OMP END PARALLEL

    deallocate(tmpVa)

    ! int:  ee,ne,se,nw,sw,sse
    allocate(tmpVi(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'

    tmpVi=ee
    deallocate(ee)
    allocate(ee(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,ee,tmpVi)
!$OMP END PARALLEL

    tmpVi=ne
    deallocate(ne)
    allocate(ne(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,ne,tmpVi)
!$OMP END PARALLEL

    tmpVi=se
    deallocate(se)
    allocate(se(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,se,tmpVi)
!$OMP END PARALLEL

    tmpVi=nw
    deallocate(nw)
    allocate(nw(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,nw,tmpVi)
!$OMP END PARALLEL

    tmpVi=sw
    deallocate(sw)
    allocate(sw(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,sw,tmpVi)
!$OMP END PARALLEL

    tmpVi=sse
    deallocate(sse)
    allocate(sse(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,sse,tmpVi)
!$OMP END PARALLEL

    deallocate(tmpVi)

    allocate(tmpVll(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
    ! logical: skipTcell1d,  skipUcell1d

    tmpVll=skipTcell1d
    deallocate(skipTcell1d)
    allocate(skipTcell1d(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,skipTcell1d,tmpVll)
!$OMP END PARALLEL

    tmpVll=skipUcell1d
    deallocate(skipUcell1d)
    allocate(skipUcell1d(nanew),stat=ierr)
    if (ierr/=0) stop 'Error allocating'
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,nanew,skipUcell1d,tmpVll)
!$OMP END PARALLEL

    deallocate(tmpVll)
  end subroutine numareinit
end module numa
