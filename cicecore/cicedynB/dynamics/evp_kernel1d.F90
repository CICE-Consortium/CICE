! Modules used for:
!   * convert 2D arrays into 1D vectors
!   * Do stress/stepu interations
!   * convert 1D vectors into 2D matrices
!
!   Call from ice_dyn_evp.F90:
!     call evp_copyin(..)
!     call evp_kernel()
!     call evp_copyout(...)
!
! Jacob Weismann Poulsen (JWP), Mads Hvid Ribergaard (MHRI), DMI
!===============================================================================

!===============================================================================
module dmi_omp
  !- directives ----------------------------------------------------------------
  use ice_kinds_mod
  implicit none
  private
  INTEGER, PARAMETER :: JPIM = SELECTED_INT_KIND(9)

  real(kind=dbl_kind), parameter  :: zero = 0.0_dbl_kind,  &
     one = 1.0_dbl_kind, half = 0.5_dbl_kind

  !- interfaces ----------------------------------------------------------------
  interface domp_get_domain
    module procedure domp_get_domain_rlu
  end interface 

  integer(int_kind), private :: domp_iam, domp_nt

  private :: domp_get_domain_rlu

  !- public vars & methods -----------------------------------------------------
  public  :: domp_init, domp_get_domain, domp_get_thread_no

#if defined (_OPENMP)
  ! Please note, this constant will create a compiler info for a constant
  ! expression in IF statements:
  real(kind=dbl_kind), private :: rdomp_iam, rdomp_nt
  !$OMP THREADPRIVATE(domp_iam,domp_nt,rdomp_iam,rdomp_nt) 
#endif

contains
    
  ! ----------------------------------------------------------------------------

  subroutine domp_init(nt_out)

#if defined (_OPENMP)
    use omp_lib, only : omp_get_thread_num, omp_get_num_threads
#endif

    !- argument(s) -------------------------------------------------------------
    integer(int_kind), intent(out) :: nt_out

    !$OMP PARALLEL DEFAULT(none)
#if defined (_OPENMP)
    domp_iam  = omp_get_thread_num()
    rdomp_iam = real(domp_iam,8)
    domp_nt   = omp_get_num_threads()
    rdomp_nt  = real(domp_nt,8)
#else
    domp_iam  = 0
    domp_nt   = 1
#endif
    !$OMP END PARALLEL

#if defined (_OPENACC)
    write(*,'(a27)') 'Build with openACC support'
!#elif defined (_OPENMP)
!    write(*,'(a26)') 'Build with openMP support'
!#else
!    write(*,'(a41)') 'Build without openMP and openACC support'
#endif

    !- echo #threads:
    if (domp_nt > 1) then
      write(*,'(a20,i5,a8)') 'Running openMP with ', domp_nt, ' threads'
    else
#if defined (_OPENMP)
      write(*,'(a35)') 'Running openMP with a single thread'
#else
      write(*,'(a22)') 'Running without openMP'
#endif
    endif

    !- return value of #threads:
    nt_out = domp_nt

  end subroutine domp_init
 
  ! ----------------------------------------------------------------------------

  subroutine domp_get_domain_rlu(lower,upper,d_lower,d_upper)
#if defined (_OPENMP)
    use omp_lib,   only : omp_in_parallel
#endif

    !- arguments ---------------------------------------------------------------
    integer(KIND=JPIM), intent(in)  :: lower,upper
    integer(KIND=JPIM), intent(out) :: d_lower,d_upper

#if defined (_OPENMP)
    !  local variables ---------------------------------------------------------
    real(kind=dbl_kind)    :: dlen
    integer(int_kind) :: lr, ur
#endif

    ! proper action in "null" cases:
    if (upper <= 0 .or. upper < lower) then
      d_lower = 0
      d_upper = -1
      return
    endif

    ! proper action in serial sections
    d_lower = lower
    d_upper = upper

#if defined (_OPENMP)
    if (omp_in_parallel()) then
      dlen    = real(upper-lower+1, 8)
      d_lower = lower    + floor((rdomp_iam*dlen+half)/rdomp_nt, 4)
      d_upper = lower -1 + floor((rdomp_iam*dlen+dlen+half)/rdomp_nt, 4)
    endif
#endif

  if (.false.) then
    write(*,'(a14,i3,a24,i10,i10)') 'openMP thread ', domp_iam,               &
         ' handles range: ', d_lower, d_upper
  endif

  end subroutine domp_get_domain_rlu

  subroutine domp_get_thread_no (tnum)
    implicit none

    integer(int_kind), intent(out) :: tnum

    tnum = domp_iam + 1

  end subroutine domp_get_thread_no
  ! ----------------------------------------------------------------------------
end module dmi_omp
!===============================================================================
module bench_v1
  !- interfaces ----------------------------------------------------------------
  interface stress
    module procedure stress_i
    module procedure stress_l
  end interface 
  interface stepu
    module procedure stepu_iter
    module procedure stepu_last
  end interface 
  contains
    subroutine stress_i(NA_len, &
                     ee,ne,se,lb,ub,uvel,vvel,dxt,dyt,     & 
                     tinyarea,dxhy,dyhx,cxp,cyp,cxm,cym,   &
                     strength,stressp_1,stressp_2,stressp_3,stressp_4,  & 
                     stressm_1,stressm_2,stressm_3,stressm_4,stress12_1,         &
                     stress12_2,stress12_3,stress12_4,str1,str2,str3,str4,str5,  &
                     str6,str7,str8)
    !- modules -------------------------------------------------------------------
    use ice_kinds_mod
    use dmi_omp, only : domp_get_domain
    use ice_constants, only: p027, p055, p111, p166, p222, p25, p333, p5, c1p5
    use icepack_parameters, only: puny
    use ice_dyn_shared, only: ecci, denom1, arlx1i
    !- directives ----------------------------------------------------------------
    implicit none
    ! arguments ------------------------------------------------------------------
    integer (kind=int_kind), intent(in) :: NA_len
    integer (kind=int_kind), intent(in) :: lb,ub
    integer (kind=int_kind), dimension(:), intent(in),    contiguous :: ee,ne,se
    real    (kind=dbl_kind), dimension(:), intent(in),    contiguous ::          &
       strength, uvel, vvel, dxt, dyt,                                           &
       tinyarea,dxhy,dyhx,cxp,cyp,cxm,cym
    real    (kind=dbl_kind), dimension(:), intent(inout), contiguous ::          &
       stressp_1,stressp_2,  stressp_3, stressp_4, stressm_1,  stressm_2,        &
       stressm_3,stressm_4, stress12_1,stress12_2,stress12_3, stress12_4
    real    (kind=dbl_kind), dimension(:), intent(out),   contiguous ::          &
       str1,str2,str3,str4,str5,str6,str7,str8
    ! local variables ------------------------------------------------------------
    integer (kind=int_kind) :: iw,il,iu
    real    (kind=dbl_kind) ::                                                   &
      divune, divunw, divuse, divusw,tensionne, tensionnw, tensionse, tensionsw, & 
      shearne, shearnw, shearse, shearsw, Deltane, Deltanw, Deltase, Deltasw   , &
      c0ne, c0nw, c0se, c0sw, c1ne, c1nw, c1se, c1sw                           , &
      ssigpn, ssigps, ssigpe, ssigpw, ssigmn, ssigms, ssigme, ssigmw           , &
      ssig12n, ssig12s, ssig12e, ssig12w, ssigp1, ssigp2,ssigm1, ssigm2,ssig121, &
      ssig122, csigpne, csigpnw, csigpse, csigpsw,csigmne, csigmnw, csigmse    , &
      csigmsw, csig12ne, csig12nw, csig12se, csig12sw, str12ew, str12we,str12ns, &
      str12sn, strp_tmp, strm_tmp, tmp_uvel_ee, tmp_vvel_se, tmp_vvel_ee,        &
      tmp_vvel_ne, tmp_uvel_ne, tmp_uvel_se
  
!    real    (kind=dbl_kind) :: dxhy,dyhx,cxp,cyp,cxm,cym,tinyarea
   
#ifdef _OPENACC
    !$acc parallel                                                                 &
    !$acc present(ee,ne,se,strength,uvel,vvel,dxt,dyt,                             &
    !$acc         tinyarea,dxhy,dyhx,cxp,cyp,cxm,cym,                              &
    !$acc         str1,str2,str3,str4,str5,str6,str7,str8,                         &
    !$acc         stressp_1,stressp_2,stressp_3,stressp_4,                         &
    !$acc         stressm_1,stressm_2,stressm_3,stressm_4,                         &
    !$acc         stress12_1,stress12_2,stress12_3,stress12_4)
    !$acc loop 
    do iw = 1,NA_len
#else
    call domp_get_domain(lb,ub,il,iu)
    do iw = il, iu
#endif
!      tinyarea = puny*dxt(iw)*dyt(iw)
!      dxhy =     p5*(hte(iw)  - htem1(iw))
!      dyhx =     p5*(htn(iw)  - htnm1(iw))
!      cxp  =   c1p5*htn(iw)   - p5*htnm1(iw)
!      cyp  =   c1p5*hte(iw)   - p5*htem1(iw)
!      cxm  = -(c1p5*htnm1(iw) - p5*htn(iw))
!      cym  = -(c1p5*htem1(iw) - p5*hte(iw))
  
      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------
      ! divergence  =  e_11 + e_22
      tmp_uvel_ee = uvel(ee(iw))
      tmp_vvel_ee = vvel(ee(iw))
  
      tmp_vvel_se = vvel(se(iw))
      tmp_uvel_se = uvel(se(iw))
  
       ! ne
      divune    = cyp(iw)*uvel(iw) - dyt(iw)*tmp_uvel_ee                       &
                + cxp(iw)*vvel(iw) - dxt(iw)*tmp_vvel_se
      ! tension strain rate  =  e_11 - e_22
      tensionne = -cym(iw)*uvel(iw) - dyt(iw)*tmp_uvel_ee                      &
                +  cxm(iw)*vvel(iw) + dxt(iw)*tmp_vvel_se
      ! shearing strain rate  =  e_12
      shearne = -cym(iw)*vvel(iw) - dyt(iw)*tmp_vvel_ee                        &
              -  cxm(iw)*uvel(iw) - dxt(iw)*tmp_uvel_se
      ! Delta (in the denominator of zeta, eta)
      Deltane = sqrt(divune**2 + ecci*(tensionne**2 + shearne**2))
  
      ! These two can move after ne block
      !
      tmp_uvel_ne = uvel(ne(iw))
      tmp_vvel_ne = vvel(ne(iw))
  
      ! nw
      divunw    = cym(iw)*tmp_uvel_ee + dyt(iw)*uvel(iw)                       &
                + cxp(iw)*tmp_vvel_ee - dxt(iw)*tmp_vvel_ne
      tensionnw = -cyp(iw)*tmp_uvel_ee + dyt(iw)*uvel(iw)                      &
                +  cxm(iw)*tmp_vvel_ee + dxt(iw)*tmp_vvel_ne
      shearnw = -cyp(iw)*tmp_vvel_ee + dyt(iw)*vvel(iw)                        &
              -  cxm(iw)*tmp_uvel_ee - dxt(iw)*tmp_uvel_ne
      Deltanw = sqrt(divunw**2 + ecci*(tensionnw**2 + shearnw**2))
  
      ! sw
      divusw    = cym(iw)*tmp_uvel_ne + dyt(iw)*tmp_uvel_se                   &
                + cxm(iw)*tmp_vvel_ne + dxt(iw)*tmp_vvel_ee
      tensionsw = -cyp(iw)*tmp_uvel_ne + dyt(iw)*tmp_uvel_se                  &
                +  cxp(iw)*tmp_vvel_ne - dxt(iw)*tmp_vvel_ee
      shearsw = -cyp(iw)*tmp_vvel_ne + dyt(iw)*tmp_vvel_se                    &
              -  cxp(iw)*tmp_uvel_ne + dxt(iw)*tmp_uvel_ee
      Deltasw = sqrt(divusw**2 + ecci*(tensionsw**2 + shearsw**2))
  
      ! se
      divuse    = cyp(iw)*tmp_uvel_se - dyt(iw)*tmp_uvel_ne                   &
                + cxm(iw)*tmp_vvel_se + dxt(iw)*vvel(iw)
      tensionse = -cym(iw)*tmp_uvel_se - dyt(iw)*tmp_uvel_ne                  &
                +  cxp(iw)*tmp_vvel_se - dxt(iw)*vvel(iw)
      shearse = -cym(iw)*tmp_vvel_se - dyt(iw)*tmp_vvel_ne                    &
              -  cxp(iw)*tmp_uvel_se + dxt(iw)*uvel(iw)
      Deltase = sqrt(divuse**2 + ecci*(tensionse**2 + shearse**2))
  
    !-----------------------------------------------------------------
    ! replacement pressure/Delta                   ! kg/s
    ! save replacement pressure for principal stress calculation
    !-----------------------------------------------------------------
      c0ne = strength(iw)/max(Deltane,tinyarea(iw))
      c0nw = strength(iw)/max(Deltanw,tinyarea(iw))
      c0sw = strength(iw)/max(Deltasw,tinyarea(iw))
      c0se = strength(iw)/max(Deltase,tinyarea(iw))
  
      c1ne = c0ne*arlx1i
      c1nw = c0nw*arlx1i
      c1sw = c0sw*arlx1i
      c1se = c0se*arlx1i
  
      c0ne = c1ne*ecci
      c0nw = c1nw*ecci
      c0sw = c1sw*ecci
      c0se = c1se*ecci
  
    !-----------------------------------------------------------------
    ! the stresses                            ! kg/s^2
    ! (1) northeast, (2) northwest, (3) southwest, (4) southeast
    !-----------------------------------------------------------------
  
      stressp_1(iw) = (stressp_1(iw) + c1ne*(divune - Deltane)) * denom1
      stressp_2(iw) = (stressp_2(iw) + c1nw*(divunw - Deltanw)) * denom1
      stressp_3(iw) = (stressp_3(iw) + c1sw*(divusw - Deltasw)) * denom1
      stressp_4(iw) = (stressp_4(iw) + c1se*(divuse - Deltase)) * denom1
  
      stressm_1(iw) = (stressm_1(iw) + c0ne*tensionne) * denom1
      stressm_2(iw) = (stressm_2(iw) + c0nw*tensionnw) * denom1
      stressm_3(iw) = (stressm_3(iw) + c0sw*tensionsw) * denom1
      stressm_4(iw) = (stressm_4(iw) + c0se*tensionse) * denom1
      
      stress12_1(iw) = (stress12_1(iw) + c0ne*shearne*p5) * denom1
      stress12_2(iw) = (stress12_2(iw) + c0nw*shearnw*p5) * denom1
      stress12_3(iw) = (stress12_3(iw) + c0sw*shearsw*p5) * denom1
      stress12_4(iw) = (stress12_4(iw) + c0se*shearse*p5) * denom1
  
    !-----------------------------------------------------------------
    ! combinations of the stresses for the momentum equation ! kg/s^2
    !-----------------------------------------------------------------
  
      ssigpn  = stressp_1(iw) + stressp_2(iw)
      ssigps  = stressp_3(iw) + stressp_4(iw)
      ssigpe  = stressp_1(iw) + stressp_4(iw)
      ssigpw  = stressp_2(iw) + stressp_3(iw)
      ssigp1  =(stressp_1(iw) + stressp_3(iw))*p055
      ssigp2  =(stressp_2(iw) + stressp_4(iw))*p055
  
      ssigmn  = stressm_1(iw) + stressm_2(iw)
      ssigms  = stressm_3(iw) + stressm_4(iw)
      ssigme  = stressm_1(iw) + stressm_4(iw)
      ssigmw  = stressm_2(iw) + stressm_3(iw)
      ssigm1  =(stressm_1(iw) + stressm_3(iw))*p055
      ssigm2  =(stressm_2(iw) + stressm_4(iw))*p055
  
      ssig12n = stress12_1(iw) + stress12_2(iw)
      ssig12s = stress12_3(iw) + stress12_4(iw)
      ssig12e = stress12_1(iw) + stress12_4(iw)
      ssig12w = stress12_2(iw) + stress12_3(iw)
      ssig121 =(stress12_1(iw) + stress12_3(iw))*p111
      ssig122 =(stress12_2(iw) + stress12_4(iw))*p111
  
      csigpne = p111*stressp_1(iw) + ssigp2 + p027*stressp_3(iw)
      csigpnw = p111*stressp_2(iw) + ssigp1 + p027*stressp_4(iw)
      csigpsw = p111*stressp_3(iw) + ssigp2 + p027*stressp_1(iw)
      csigpse = p111*stressp_4(iw) + ssigp1 + p027*stressp_2(iw)
      
      csigmne = p111*stressm_1(iw) + ssigm2 + p027*stressm_3(iw)
      csigmnw = p111*stressm_2(iw) + ssigm1 + p027*stressm_4(iw)
      csigmsw = p111*stressm_3(iw) + ssigm2 + p027*stressm_1(iw)
      csigmse = p111*stressm_4(iw) + ssigm1 + p027*stressm_2(iw)
      
      csig12ne = p222*stress12_1(iw) + ssig122 + p055*stress12_3(iw)
      csig12nw = p222*stress12_2(iw) + ssig121 + p055*stress12_4(iw)
      csig12sw = p222*stress12_3(iw) + ssig122 + p055*stress12_1(iw)
      csig12se = p222*stress12_4(iw) + ssig121 + p055*stress12_2(iw)
  
      str12ew = p5*dxt(iw)*(p333*ssig12e + p166*ssig12w)
      str12we = p5*dxt(iw)*(p333*ssig12w + p166*ssig12e)
      str12ns = p5*dyt(iw)*(p333*ssig12n + p166*ssig12s)
      str12sn = p5*dyt(iw)*(p333*ssig12s + p166*ssig12n)
  
    !-----------------------------------------------------------------
    ! for dF/dx (u momentum)
    !-----------------------------------------------------------------
      strp_tmp  = p25*dyt(iw)*(p333*ssigpn  + p166*ssigps)
      strm_tmp  = p25*dyt(iw)*(p333*ssigmn  + p166*ssigms)
  
      ! northeast (iw)
      str1(iw) = -strp_tmp - strm_tmp - str12ew &
                 + dxhy(iw)*(-csigpne + csigmne) + dyhx(iw)*csig12ne
  
      ! northwest (i+1,j)
      str2(iw) = strp_tmp + strm_tmp - str12we &
                 + dxhy(iw)*(-csigpnw + csigmnw) + dyhx(iw)*csig12nw
  
      strp_tmp  = p25*dyt(iw)*(p333*ssigps  + p166*ssigpn)
      strm_tmp  = p25*dyt(iw)*(p333*ssigms  + p166*ssigmn)
  
      ! southeast (i,j+1)
      str3(iw) = -strp_tmp - strm_tmp + str12ew &
                 + dxhy(iw)*(-csigpse + csigmse) + dyhx(iw)*csig12se
  
      ! southwest (i+1,j+1)
      str4(iw) = strp_tmp + strm_tmp + str12we &
                 + dxhy(iw)*(-csigpsw + csigmsw) + dyhx(iw)*csig12sw
  
    !-----------------------------------------------------------------
    ! for dF/dy (v momentum)
    !-----------------------------------------------------------------
      strp_tmp  = p25*dxt(iw)*(p333*ssigpe  + p166*ssigpw)
      strm_tmp  = p25*dxt(iw)*(p333*ssigme  + p166*ssigmw)
  
      ! northeast (i,j)
      str5(iw) = -strp_tmp + strm_tmp - str12ns &
                 - dyhx(iw)*(csigpne + csigmne) + dxhy(iw)*csig12ne
  
      ! southeast (i,j+1)
      str6(iw) = strp_tmp - strm_tmp - str12sn &
                 - dyhx(iw)*(csigpse + csigmse) + dxhy(iw)*csig12se
  
      strp_tmp  = p25*dxt(iw)*(p333*ssigpw  + p166*ssigpe)
      strm_tmp  = p25*dxt(iw)*(p333*ssigmw  + p166*ssigme)
  
      ! northwest (i+1,j)
      str7(iw) = -strp_tmp + strm_tmp + str12ns &
                 - dyhx(iw)*(csigpnw + csigmnw) + dxhy(iw)*csig12nw
  
      ! southwest (i+1,j+1)
      str8(iw) = strp_tmp - strm_tmp + str12sn &
                 - dyhx(iw)*(csigpsw + csigmsw) + dxhy(iw)*csig12sw
    enddo   
    !$acc end parallel
  end subroutine stress_i
  
  subroutine stress_l(NA_len, tarear, &
                     ee,ne,se,lb,ub,uvel,vvel,dxt,dyt,                           & 
                     tinyarea,dxhy,dyhx,cxp,cyp,cxm,cym,                         &
                     strength,stressp_1,stressp_2,stressp_3,stressp_4,           & 
                     stressm_1,stressm_2,stressm_3,stressm_4,stress12_1,         &
                     stress12_2,stress12_3,stress12_4,                           &
                     divu,rdg_conv,rdg_shear,shear,                              &
                     str1,str2,str3,str4,str5,str6,str7,str8                     )
    !- modules -------------------------------------------------------------------
    use ice_kinds_mod
    use dmi_omp, only : domp_get_domain
    use ice_constants, only: p027, p055, p111, p166, p222, p25, p333, p5, c1p5, c0
    use icepack_parameters, only: puny
    use ice_dyn_shared, only: ecci, denom1, arlx1i
    !- directives ----------------------------------------------------------------
    implicit none
    ! arguments ------------------------------------------------------------------
    integer (kind=int_kind), intent(in) :: NA_len
    integer (kind=int_kind), intent(in) :: lb,ub
    integer (kind=int_kind), dimension(:), intent(in),    contiguous :: ee,ne,se
    real    (kind=dbl_kind), dimension(:), intent(in),    contiguous ::          &
       strength, uvel, vvel, dxt, dyt, tarear, tinyarea,dxhy,dyhx,cxp,cyp,cxm,cym
    real    (kind=dbl_kind), dimension(:), intent(inout), contiguous ::          &
       stressp_1,stressp_2,  stressp_3, stressp_4, stressm_1,  stressm_2,        &
       stressm_3,stressm_4, stress12_1,stress12_2,stress12_3, stress12_4
    real    (kind=dbl_kind), dimension(:), intent(out),   contiguous ::          &
            str1,str2,str3,str4,str5,str6,str7,str8,                             &
            divu,rdg_conv,rdg_shear,shear
    ! local variables ------------------------------------------------------------
    integer (kind=int_kind) :: iw,il,iu
    real    (kind=dbl_kind) ::                                                   &
      divune, divunw, divuse, divusw,tensionne, tensionnw, tensionse, tensionsw, & 
      shearne, shearnw, shearse, shearsw, Deltane, Deltanw, Deltase, Deltasw   , &
      c0ne, c0nw, c0se, c0sw, c1ne, c1nw, c1se, c1sw                           , &
      ssigpn, ssigps, ssigpe, ssigpw, ssigmn, ssigms, ssigme, ssigmw           , &
      ssig12n, ssig12s, ssig12e, ssig12w, ssigp1, ssigp2,ssigm1, ssigm2,ssig121, &
      ssig122, csigpne, csigpnw, csigpse, csigpsw,csigmne, csigmnw, csigmse    , &
      csigmsw, csig12ne, csig12nw, csig12se, csig12sw, str12ew, str12we,str12ns, &
      str12sn, strp_tmp, strm_tmp, tmp_uvel_ee, tmp_vvel_se, tmp_vvel_ee,        &
      tmp_vvel_ne, tmp_uvel_ne, tmp_uvel_se
  
!    real    (kind=dbl_kind) :: dxhy,dyhx,cxp,cyp,cxm,cym,tinyarea
  
#ifdef _OPENACC
    !$acc parallel                                                                 &
    !$acc present(ee,ne,se,strength,uvel,vvel,dxt,dyt,tarear,                      &
    !$acc         tinyarea,dxhy,dyhx,cxp,cyp,cxm,cym,                              &
    !$acc         str1,str2,str3,str4,str5,str6,str7,str8,                         &
    !$acc         stressp_1,stressp_2,stressp_3,stressp_4,                         &
    !$acc         stressm_1,stressm_2,stressm_3,stressm_4,                         &
    !$acc         stress12_1,stress12_2,stress12_3,stress12_4,                     &
    !$acc         divu,rdg_conv,rdg_shear,shear)
    !$acc loop 
    do iw = 1,NA_len
#else
    call domp_get_domain(lb,ub,il,iu)
    do iw = il, iu
#endif
!      tinyarea = puny*dxt(iw)*dyt(iw)
!      dxhy =     p5*(hte(iw)  - htem1(iw))
!      dyhx =     p5*(htn(iw)  - htnm1(iw))
!      cxp  =   c1p5*htn(iw)   - p5*htnm1(iw)
!      cyp  =   c1p5*hte(iw)   - p5*htem1(iw)
!      cxm  = -(c1p5*htnm1(iw) - p5*htn(iw))
!      cym  = -(c1p5*htem1(iw) - p5*hte(iw))
  
    !-----------------------------------------------------------------
    ! strain rates
    ! NOTE these are actually strain rates * area  (m^2/s)
    !-----------------------------------------------------------------
       ! divergence  =  e_11 + e_22
       tmp_uvel_ee = uvel(ee(iw))
       tmp_vvel_se = vvel(se(iw))
       tmp_vvel_ee = vvel(ee(iw))
       tmp_vvel_ne = vvel(ne(iw))
       tmp_uvel_ne = uvel(ne(iw))
       tmp_uvel_se = uvel(se(iw))
  
       divune    = cyp(iw)*uvel(iw) - dyt(iw)*tmp_uvel_ee                       &
                 + cxp(iw)*vvel(iw) - dxt(iw)*tmp_vvel_se
       divunw    = cym(iw)*tmp_uvel_ee + dyt(iw)*uvel(iw)                       &
                 + cxp(iw)*tmp_vvel_ee - dxt(iw)*tmp_vvel_ne
       divusw    = cym(iw)*tmp_uvel_ne + dyt(iw)*tmp_uvel_se                   &
                 + cxm(iw)*tmp_vvel_ne + dxt(iw)*tmp_vvel_ee
       divuse    = cyp(iw)*tmp_uvel_se - dyt(iw)*tmp_uvel_ne                   &
                 + cxm(iw)*tmp_vvel_se + dxt(iw)*vvel(iw)
  
       ! tension strain rate  =  e_11 - e_22
       tensionne = -cym(iw)*uvel(iw) - dyt(iw)*tmp_uvel_ee                      &
                 +  cxm(iw)*vvel(iw) + dxt(iw)*tmp_vvel_se
       tensionnw = -cyp(iw)*tmp_uvel_ee + dyt(iw)*uvel(iw)                      &
                 +  cxm(iw)*tmp_vvel_ee + dxt(iw)*tmp_vvel_ne
       tensionsw = -cyp(iw)*tmp_uvel_ne + dyt(iw)*tmp_uvel_se                  &
                 +  cxp(iw)*tmp_vvel_ne - dxt(iw)*tmp_vvel_ee
       tensionse = -cym(iw)*tmp_uvel_se - dyt(iw)*tmp_uvel_ne                  &
                 +  cxp(iw)*tmp_vvel_se - dxt(iw)*vvel(iw)
  
       ! shearing strain rate  =  e_12
       shearne = -cym(iw)*vvel(iw) - dyt(iw)*tmp_vvel_ee                        &
               -  cxm(iw)*uvel(iw) - dxt(iw)*tmp_uvel_se
       shearnw = -cyp(iw)*tmp_vvel_ee + dyt(iw)*vvel(iw)                        &
               -  cxm(iw)*tmp_uvel_ee - dxt(iw)*tmp_uvel_ne
       shearsw = -cyp(iw)*tmp_vvel_ne + dyt(iw)*tmp_vvel_se                    &
               -  cxp(iw)*tmp_uvel_ne + dxt(iw)*tmp_uvel_ee
       shearse = -cym(iw)*tmp_vvel_se - dyt(iw)*tmp_vvel_ne                    &
               -  cxp(iw)*tmp_uvel_se + dxt(iw)*uvel(iw)
       
       ! Delta (in the denominator of zeta, eta)
       Deltane = sqrt(divune**2 + ecci*(tensionne**2 + shearne**2))
       Deltanw = sqrt(divunw**2 + ecci*(tensionnw**2 + shearnw**2))
       Deltase = sqrt(divuse**2 + ecci*(tensionse**2 + shearse**2))
       Deltasw = sqrt(divusw**2 + ecci*(tensionsw**2 + shearsw**2))
  
!       !-----------------------------------------------------------------
!       ! on last subcycle, save quantities for mechanical redistribution
!       !-----------------------------------------------------------------
       divu(iw) = p25*(divune + divunw + divuse + divusw) * tarear(iw)
       rdg_conv(iw)  = -min(divu(iw),c0)                          ! Move outside the entire "kernel"
       rdg_shear(iw) = p5*( p25*(Deltane + Deltanw + Deltase + Deltasw) * tarear(iw) -abs(divu(iw)) )

       ! diagnostic only
       ! shear = sqrt(tension**2 + shearing**2)
       shear(iw) = p25*tarear(iw)*sqrt( &
                 (tensionne + tensionnw + tensionse + tensionsw)**2 &
                +  (shearne +   shearnw +   shearse +   shearsw)**2)

    !-----------------------------------------------------------------
    ! replacement pressure/Delta                   ! kg/s
    ! save replacement pressure for principal stress calculation
    !-----------------------------------------------------------------
       c0ne = strength(iw)/max(Deltane,tinyarea(iw))
       c0nw = strength(iw)/max(Deltanw,tinyarea(iw))
       c0sw = strength(iw)/max(Deltasw,tinyarea(iw))
       c0se = strength(iw)/max(Deltase,tinyarea(iw))
  
       c1ne = c0ne*arlx1i
       c1nw = c0nw*arlx1i
       c1sw = c0sw*arlx1i
       c1se = c0se*arlx1i
  
       c0ne = c1ne*ecci
       c0nw = c1nw*ecci
       c0sw = c1sw*ecci
       c0se = c1se*ecci
  
    !-----------------------------------------------------------------
    ! the stresses                            ! kg/s^2
    ! (1) northeast, (2) northwest, (3) southwest, (4) southeast
    !-----------------------------------------------------------------
  
       stressp_1(iw) = (stressp_1(iw) + c1ne*(divune - Deltane)) * denom1
       stressp_2(iw) = (stressp_2(iw) + c1nw*(divunw - Deltanw)) * denom1
       stressp_3(iw) = (stressp_3(iw) + c1sw*(divusw - Deltasw)) * denom1
       stressp_4(iw) = (stressp_4(iw) + c1se*(divuse - Deltase)) * denom1
  
       stressm_1(iw) = (stressm_1(iw) + c0ne*tensionne) * denom1
       stressm_2(iw) = (stressm_2(iw) + c0nw*tensionnw) * denom1
       stressm_3(iw) = (stressm_3(iw) + c0sw*tensionsw) * denom1
       stressm_4(iw) = (stressm_4(iw) + c0se*tensionse) * denom1
      
       stress12_1(iw) = (stress12_1(iw) + c0ne*shearne*p5) * denom1
       stress12_2(iw) = (stress12_2(iw) + c0nw*shearnw*p5) * denom1
       stress12_3(iw) = (stress12_3(iw) + c0sw*shearsw*p5) * denom1
       stress12_4(iw) = (stress12_4(iw) + c0se*shearse*p5) * denom1
  
    !-----------------------------------------------------------------
    ! combinations of the stresses for the momentum equation ! kg/s^2
    !-----------------------------------------------------------------
  
       ssigpn  = stressp_1(iw) + stressp_2(iw)
       ssigps  = stressp_3(iw) + stressp_4(iw)
       ssigpe  = stressp_1(iw) + stressp_4(iw)
       ssigpw  = stressp_2(iw) + stressp_3(iw)
       ssigp1  =(stressp_1(iw) + stressp_3(iw))*p055
       ssigp2  =(stressp_2(iw) + stressp_4(iw))*p055
  
       ssigmn  = stressm_1(iw) + stressm_2(iw)
       ssigms  = stressm_3(iw) + stressm_4(iw)
       ssigme  = stressm_1(iw) + stressm_4(iw)
       ssigmw  = stressm_2(iw) + stressm_3(iw)
       ssigm1  =(stressm_1(iw) + stressm_3(iw))*p055
       ssigm2  =(stressm_2(iw) + stressm_4(iw))*p055
  
       ssig12n = stress12_1(iw) + stress12_2(iw)
       ssig12s = stress12_3(iw) + stress12_4(iw)
       ssig12e = stress12_1(iw) + stress12_4(iw)
       ssig12w = stress12_2(iw) + stress12_3(iw)
       ssig121 =(stress12_1(iw) + stress12_3(iw))*p111
       ssig122 =(stress12_2(iw) + stress12_4(iw))*p111
  
       csigpne = p111*stressp_1(iw) + ssigp2 + p027*stressp_3(iw)
       csigpnw = p111*stressp_2(iw) + ssigp1 + p027*stressp_4(iw)
       csigpsw = p111*stressp_3(iw) + ssigp2 + p027*stressp_1(iw)
       csigpse = p111*stressp_4(iw) + ssigp1 + p027*stressp_2(iw)
       
       csigmne = p111*stressm_1(iw) + ssigm2 + p027*stressm_3(iw)
       csigmnw = p111*stressm_2(iw) + ssigm1 + p027*stressm_4(iw)
       csigmsw = p111*stressm_3(iw) + ssigm2 + p027*stressm_1(iw)
       csigmse = p111*stressm_4(iw) + ssigm1 + p027*stressm_2(iw)
       
       csig12ne = p222*stress12_1(iw) + ssig122 + p055*stress12_3(iw)
       csig12nw = p222*stress12_2(iw) + ssig121 + p055*stress12_4(iw)
       csig12sw = p222*stress12_3(iw) + ssig122 + p055*stress12_1(iw)
       csig12se = p222*stress12_4(iw) + ssig121 + p055*stress12_2(iw)
  
       str12ew = p5*dxt(iw)*(p333*ssig12e + p166*ssig12w)
       str12we = p5*dxt(iw)*(p333*ssig12w + p166*ssig12e)
       str12ns = p5*dyt(iw)*(p333*ssig12n + p166*ssig12s)
       str12sn = p5*dyt(iw)*(p333*ssig12s + p166*ssig12n)
  
    !-----------------------------------------------------------------
    ! for dF/dx (u momentum)
    !-----------------------------------------------------------------
       strp_tmp  = p25*dyt(iw)*(p333*ssigpn  + p166*ssigps)
       strm_tmp  = p25*dyt(iw)*(p333*ssigmn  + p166*ssigms)
  
       ! northeast (iw)
       str1(iw) = -strp_tmp - strm_tmp - str12ew &
                  + dxhy(iw)*(-csigpne + csigmne) + dyhx(iw)*csig12ne
  
       ! northwest (i+1,j)
       str2(iw) = strp_tmp + strm_tmp - str12we &
                  + dxhy(iw)*(-csigpnw + csigmnw) + dyhx(iw)*csig12nw
  
       strp_tmp  = p25*dyt(iw)*(p333*ssigps  + p166*ssigpn)
       strm_tmp  = p25*dyt(iw)*(p333*ssigms  + p166*ssigmn)
  
       ! southeast (i,j+1)
       str3(iw) = -strp_tmp - strm_tmp + str12ew &
                  + dxhy(iw)*(-csigpse + csigmse) + dyhx(iw)*csig12se
  
       ! southwest (i+1,j+1)
       str4(iw) = strp_tmp + strm_tmp + str12we &
                  + dxhy(iw)*(-csigpsw + csigmsw) + dyhx(iw)*csig12sw
  
    !-----------------------------------------------------------------
    ! for dF/dy (v momentum)
    !-----------------------------------------------------------------
       strp_tmp  = p25*dxt(iw)*(p333*ssigpe  + p166*ssigpw)
       strm_tmp  = p25*dxt(iw)*(p333*ssigme  + p166*ssigmw)
  
       ! northeast (i,j)
       str5(iw) = -strp_tmp + strm_tmp - str12ns &
                  - dyhx(iw)*(csigpne + csigmne) + dxhy(iw)*csig12ne
  
       ! southeast (i,j+1)
       str6(iw) = strp_tmp - strm_tmp - str12sn &
                  - dyhx(iw)*(csigpse + csigmse) + dxhy(iw)*csig12se
  
       strp_tmp  = p25*dxt(iw)*(p333*ssigpw  + p166*ssigpe)
       strm_tmp  = p25*dxt(iw)*(p333*ssigmw  + p166*ssigme)
  
       ! northwest (i+1,j)
       str7(iw) = -strp_tmp + strm_tmp + str12ns &
                  - dyhx(iw)*(csigpnw + csigmnw) + dxhy(iw)*csig12nw
  
       ! southwest (i+1,j+1)
       str8(iw) = strp_tmp - strm_tmp + str12sn &
                  - dyhx(iw)*(csigpsw + csigmsw) + dxhy(iw)*csig12sw
    enddo   
    !$acc end parallel
  end subroutine stress_l
  
  subroutine stepu_iter(NA_len,rhow, &
                    lb,ub,Cw,aiu,uocn,vocn,forcex,forcey, umassdti,fm,uarear,    &
                    uvel_init,vvel_init,uvel,vvel,                               &
                    str1,str2,str3,str4,str5,str6,str7,str8, nw,sw,se,skipme)
    !- modules -------------------------------------------------------------------
    use ice_kinds_mod
    use dmi_omp, only : domp_get_domain
    use ice_dyn_shared, only: ecci, denom1, arlx1i, brlx, revp
    use ice_constants, only: c0, c1
    !- directives ----------------------------------------------------------------
    implicit none
    ! arguments ------------------------------------------------------------------
    integer (kind=int_kind), intent(in) :: NA_len
    real (kind=dbl_kind), intent(in) :: rhow
    integer(kind=int_kind),intent(in)   :: lb,ub
    logical(kind=log_kind),intent(in), dimension(:)   :: skipme
    integer(kind=int_kind),dimension(:), intent(in), contiguous :: nw,sw,se
    real(kind=dbl_kind),dimension(:), intent(in), contiguous    ::               &
       uvel_init, vvel_init, aiu, forcex, forcey, umassdti,                      &
       uocn, vocn, fm, uarear,Cw,str1,str2,str3,str4,str5,str6,str7,str8
    real(kind=dbl_kind),dimension(:), intent(inout), contiguous ::               &
       uvel,vvel
    real (kind=dbl_kind), parameter :: &
           cosw = c1   , & ! cos(ocean turning angle)  ! turning angle = 0
           sinw = c0   
    ! local variables
    integer (kind=int_kind) :: iw,il,iu
    real    (kind=dbl_kind) :: uold, vold, vrel,cca,ccb,ab2,cc1,cc2,taux,tauy
    real    (kind=dbl_kind) :: tmp_str2_nw,tmp_str3_se,tmp_str4_sw, tmp_strintx
    real    (kind=dbl_kind) :: tmp_str6_se,tmp_str7_nw,tmp_str8_sw, tmp_strinty
    real    (kind=dbl_kind) :: waterx,watery
  
#ifdef _OPENACC
    !$acc parallel                                                                 &
    !$acc present(Cw,aiu,uocn,vocn,forcex,forcey,umassdti,                         &
    !$acc         fm,uarear,uvel_init,vvel_init,nw,sw,se,skipme,                   &
    !$acc         str1,str2,str3,str4,str5,str6,str7,str8,uvel,vvel)
    !$acc loop 
    do iw = 1,NA_len
#else
    call domp_get_domain(lb,ub,il,iu)
    do iw = il, iu
#endif
       if (skipme(iw)) cycle
       uold = uvel(iw)
       vold = vvel(iw)
       vrel = aiu(iw)*rhow*Cw(iw)*sqrt((uocn(iw)-uold)**2+(vocn(iw)-vold)**2)
       waterx = uocn(iw)*cosw - vocn(iw)*sinw*sign(c1,fm(iw))
       watery = vocn(iw)*cosw + uocn(iw)*sinw*sign(c1,fm(iw))
       taux = vrel*waterx 
       tauy = vrel*watery
       cca = (brlx + revp)*umassdti(iw) + vrel * cosw 
       ccb = fm(iw) + sign(c1,fm(iw)) * vrel * sinw 
       ab2 = cca**2 + ccb**2
       ! southeast(i,j+1)   = se
       ! northwest(i+1,j)   = nw
       ! southwest(i+1,j+1) = sw
       tmp_str2_nw = str2(nw(iw))
       tmp_str3_se = str3(se(iw))
       tmp_str4_sw = str4(sw(iw))
       tmp_str6_se = str6(se(iw))
       tmp_str7_nw = str7(nw(iw))
       tmp_str8_sw = str8(sw(iw))
  
       tmp_strintx = uarear(iw)*(str1(iw)+tmp_str2_nw+tmp_str3_se+tmp_str4_sw)
       tmp_strinty = uarear(iw)*(str5(iw)+tmp_str6_se+tmp_str7_nw+tmp_str8_sw)
       cc1 = tmp_strintx + forcex(iw) + taux &
           + umassdti(iw)*(brlx*uold + revp*uvel_init(iw))
       cc2 = tmp_strinty + forcey(iw) + tauy &
           + umassdti(iw)*(brlx*vold + revp*vvel_init(iw))
       uvel(iw) = (cca*cc1 + ccb*cc2) / ab2 
       vvel(iw) = (cca*cc2 - ccb*cc1) / ab2
    enddo
    !$acc end parallel
  end subroutine stepu_iter
  
  subroutine stepu_last(NA_len, rhow, &
                    lb,ub,Cw,aiu,uocn,vocn,forcex,forcey, umassdti,fm,uarear,    &
                    strocnx,strocny,strintx,strinty,                             &
                    uvel_init,vvel_init,uvel,vvel,                               &
                    str1,str2,str3,str4,str5,str6,str7,str8, nw,sw,se,skipme)
    !- modules -------------------------------------------------------------------
    use ice_kinds_mod
    use dmi_omp, only : domp_get_domain
    use ice_constants, only: c0, c1
    use icepack_parameters, only: puny
    use ice_dyn_shared, only: ecci, denom1, arlx1i, brlx, revp
    !- directives ----------------------------------------------------------------
    implicit none
    ! arguments ------------------------------------------------------------------
    integer (kind=int_kind), intent(in) :: NA_len
    real (kind=dbl_kind), intent(in) :: rhow
    logical(kind=log_kind),intent(in), dimension(:)   :: skipme
    integer(kind=int_kind),intent(in)   :: lb,ub
    integer(kind=int_kind),dimension(:), intent(in), contiguous :: nw,sw,se
    real(kind=dbl_kind),dimension(:), intent(in), contiguous    ::               &
       uvel_init, vvel_init, aiu, forcex, forcey, umassdti,      &
       uocn, vocn, fm, uarear,Cw,str1,str2,str3,str4,str5,str6,str7,str8
    real(kind=dbl_kind),dimension(:), intent(inout), contiguous ::               &
       uvel,vvel,strintx, strinty
    real(kind=dbl_kind),dimension(:), intent(out),   contiguous ::               &
       strocnx,strocny
    real (kind=dbl_kind), parameter :: &
           cosw = c1   , & ! cos(ocean turning angle)  ! turning angle = 0
           sinw = c0   
    ! local variables
    integer (kind=int_kind) :: iw,il,iu
    real    (kind=dbl_kind) :: uold, vold, vrel,cca,ccb,ab2,cc1,cc2,taux,tauy
    real    (kind=dbl_kind) :: tmp_str2_nw,tmp_str3_se,tmp_str4_sw
    real    (kind=dbl_kind) :: tmp_str6_se,tmp_str7_nw,tmp_str8_sw
    real    (kind=dbl_kind) :: waterx,watery
  
#ifdef _OPENACC
    !$acc parallel                                                                 &
    !$acc present(Cw,aiu,uocn,vocn,forcex,forcey,umassdti,                         &
    !$acc         fm,uarear,strintx,strinty,uvel_init,vvel_init,nw,sw,se,skipme,   &
    !$acc         str1,str2,str3,str4,str5,str6,str7,str8,strocnx,strocny,uvel,vvel)
    !$acc loop 
    do iw = 1,NA_len
#else
    call domp_get_domain(lb,ub,il,iu)
    do iw = il, iu
#endif
       if (skipme(iw)) cycle
       uold = uvel(iw)
       vold = vvel(iw)
       vrel = aiu(iw)*rhow*Cw(iw)*sqrt((uocn(iw)-uold)**2+(vocn(iw)-vold)**2)
       waterx = uocn(iw)*cosw - vocn(iw)*sinw*sign(c1,fm(iw))
       watery = vocn(iw)*cosw + uocn(iw)*sinw*sign(c1,fm(iw))
       taux = vrel*waterx 
       tauy = vrel*watery
       cca = (brlx + revp)*umassdti(iw) + vrel * cosw 
       ccb = fm(iw) + sign(c1,fm(iw)) * vrel * sinw 
       ab2 = cca**2 + ccb**2
       ! southeast(i,j+1)   = se
       ! northwest(i+1,j)   = nw
       ! southwest(i+1,j+1) = sw
       tmp_str2_nw = str2(nw(iw))
       tmp_str3_se = str3(se(iw))
       tmp_str4_sw = str4(sw(iw))
       tmp_str6_se = str6(se(iw))
       tmp_str7_nw = str7(nw(iw))
       tmp_str8_sw = str8(sw(iw))
  
       strintx(iw) = uarear(iw)*(str1(iw)+tmp_str2_nw+tmp_str3_se+tmp_str4_sw)
       strinty(iw) = uarear(iw)*(str5(iw)+tmp_str6_se+tmp_str7_nw+tmp_str8_sw)
       cc1 = strintx(iw) + forcex(iw) + taux &
           + umassdti(iw)*(brlx*uold + revp*uvel_init(iw))
       cc2 = strinty(iw) + forcey(iw) + tauy &
           + umassdti(iw)*(brlx*vold + revp*vvel_init(iw))
       uvel(iw) = (cca*cc1 + ccb*cc2) / ab2 
       vvel(iw) = (cca*cc2 - ccb*cc1) / ab2
       strocnx(iw) = taux
       strocny(iw) = tauy
    enddo
    !$acc end parallel
  end subroutine stepu_last

  subroutine halo_update(NAVEL_len,lb,ub,uvel,vvel, halo_parent)
    !- modules -------------------------------------------------------------------
    use ice_kinds_mod
    use dmi_omp, only : domp_get_domain
    !- directives ----------------------------------------------------------------
    implicit none
    ! arguments ------------------------------------------------------------------
    integer (kind=int_kind), intent(in) :: NAVEL_len
    integer(kind=int_kind),intent(in)   :: lb,ub
    integer(kind=int_kind),dimension(:), intent(in), contiguous :: halo_parent
    real(kind=dbl_kind),dimension(:), intent(inout), contiguous :: uvel,vvel
    ! local variables
    integer (kind=int_kind) :: iw,il,iu

#ifdef _OPENACC
    !$acc parallel                                   &
    !$acc present(uvel,vvel)                         &
    !$acc loop
    do iw = 1,NAVEL_len
#else
    call domp_get_domain(lb,ub,il,iu)
    do iw = il, iu
#endif
       if (halo_parent(iw)==0) cycle
       uvel(iw) = uvel(halo_parent(iw))
       vvel(iw) = vvel(halo_parent(iw))
    enddo
    !$acc end parallel
  end subroutine halo_update

end module bench_v1
!===============================================================================
module bench_v2
  !- interfaces ----------------------------------------------------------------
  interface stress
    module procedure stress_i
    module procedure stress_l
  end interface 
  interface stepu
    module procedure stepu_iter
    module procedure stepu_last
  end interface 
  contains
    subroutine stress_i(NA_len, &
                     ee,ne,se,lb,ub,uvel,vvel,dxt,dyt,                           & 
                     hte,htn,htem1,htnm1,                                        &
                     strength,stressp_1,stressp_2,stressp_3,stressp_4,           &
                     stressm_1,stressm_2,stressm_3,stressm_4,stress12_1,         &
                     stress12_2,stress12_3,stress12_4,str1,str2,str3,str4,str5,  &
                     str6,str7,str8)
    !- modules -------------------------------------------------------------------
    use ice_kinds_mod
    use dmi_omp, only : domp_get_domain
    use ice_constants, only: p027, p055, p111, p166, p222, p25, p333, p5, c1p5
    use icepack_parameters, only: puny
    use ice_dyn_shared, only: ecci, denom1, arlx1i
    !- directives ----------------------------------------------------------------
    implicit none
    ! arguments ------------------------------------------------------------------
    integer (kind=int_kind), intent(in) :: NA_len
    integer (kind=int_kind), intent(in) :: lb,ub
    integer (kind=int_kind), dimension(:), intent(in),    contiguous :: ee,ne,se
    real    (kind=dbl_kind), dimension(:), intent(in),    contiguous ::          &
       strength, uvel, vvel, dxt, dyt,                                           &
       hte,htn,htem1,htnm1
    real    (kind=dbl_kind), dimension(:), intent(inout), contiguous ::          &
       stressp_1,stressp_2,  stressp_3, stressp_4, stressm_1,  stressm_2,        &
       stressm_3,stressm_4, stress12_1,stress12_2,stress12_3, stress12_4
    real    (kind=dbl_kind), dimension(:), intent(out),   contiguous ::          &
       str1,str2,str3,str4,str5,str6,str7,str8
    ! local variables ------------------------------------------------------------
    integer (kind=int_kind) :: iw,il,iu
    real    (kind=dbl_kind) ::                                                   &
      divune, divunw, divuse, divusw,tensionne, tensionnw, tensionse, tensionsw, & 
      shearne, shearnw, shearse, shearsw, Deltane, Deltanw, Deltase, Deltasw   , &
      c0ne, c0nw, c0se, c0sw, c1ne, c1nw, c1se, c1sw                           , &
      ssigpn, ssigps, ssigpe, ssigpw, ssigmn, ssigms, ssigme, ssigmw           , &
      ssig12n, ssig12s, ssig12e, ssig12w, ssigp1, ssigp2,ssigm1, ssigm2,ssig121, &
      ssig122, csigpne, csigpnw, csigpse, csigpsw,csigmne, csigmnw, csigmse    , &
      csigmsw, csig12ne, csig12nw, csig12se, csig12sw, str12ew, str12we,str12ns, &
      str12sn, strp_tmp, strm_tmp, tmp_uvel_ee, tmp_vvel_se, tmp_vvel_ee,        &
      tmp_vvel_ne, tmp_uvel_ne, tmp_uvel_se
    real    (kind=dbl_kind) :: dxhy,dyhx,cxp,cyp,cxm,cym,tinyarea
   
#ifdef _OPENACC
    !$acc parallel                                                                 &
    !$acc present(ee,ne,se,strength,uvel,vvel,dxt,dyt,                             &
    !$acc         hte, htn, htem1, htnm1,                                          &
    !$acc         str1,str2,str3,str4,str5,str6,str7,str8,                         &
    !$acc         stressp_1,stressp_2,stressp_3,stressp_4,                         &
    !$acc         stressm_1,stressm_2,stressm_3,stressm_4,                         &
    !$acc         stress12_1,stress12_2,stress12_3,stress12_4)
    !$acc loop 
    do iw = 1,NA_len
#else
    call domp_get_domain(lb,ub,il,iu)
    do iw = il, iu
#endif
      tinyarea = puny*dxt(iw)*dyt(iw)
      dxhy =     p5*(hte(iw)  - htem1(iw))
      dyhx =     p5*(htn(iw)  - htnm1(iw))
      cxp  =   c1p5*htn(iw)   - p5*htnm1(iw)
      cyp  =   c1p5*hte(iw)   - p5*htem1(iw)
      cxm  = -(c1p5*htnm1(iw) - p5*htn(iw))
      cym  = -(c1p5*htem1(iw) - p5*hte(iw))
  
      !-----------------------------------------------------------------
      ! strain rates
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------
      ! divergence  =  e_11 + e_22
      tmp_uvel_ee = uvel(ee(iw))
      tmp_vvel_ee = vvel(ee(iw))
  
      tmp_vvel_se = vvel(se(iw))
      tmp_uvel_se = uvel(se(iw))
  
       ! ne
      divune    = cyp*uvel(iw) - dyt(iw)*tmp_uvel_ee                       &
                + cxp*vvel(iw) - dxt(iw)*tmp_vvel_se
      ! tension strain rate  =  e_11 - e_22
      tensionne = -cym*uvel(iw) - dyt(iw)*tmp_uvel_ee                      &
                +  cxm*vvel(iw) + dxt(iw)*tmp_vvel_se
      ! shearing strain rate  =  e_12
      shearne = -cym*vvel(iw) - dyt(iw)*tmp_vvel_ee                        &
              -  cxm*uvel(iw) - dxt(iw)*tmp_uvel_se
      ! Delta (in the denominator of zeta, eta)
      Deltane = sqrt(divune**2 + ecci*(tensionne**2 + shearne**2))
  
      ! These two can move after ne block
      !
      tmp_uvel_ne = uvel(ne(iw))
      tmp_vvel_ne = vvel(ne(iw))
  
      ! nw
      divunw    = cym*tmp_uvel_ee + dyt(iw)*uvel(iw)                       &
                + cxp*tmp_vvel_ee - dxt(iw)*tmp_vvel_ne
      tensionnw = -cyp*tmp_uvel_ee + dyt(iw)*uvel(iw)                      &
                +  cxm*tmp_vvel_ee + dxt(iw)*tmp_vvel_ne
      shearnw = -cyp*tmp_vvel_ee + dyt(iw)*vvel(iw)                        &
              -  cxm*tmp_uvel_ee - dxt(iw)*tmp_uvel_ne
      Deltanw = sqrt(divunw**2 + ecci*(tensionnw**2 + shearnw**2))
  
      ! sw
      divusw    = cym*tmp_uvel_ne + dyt(iw)*tmp_uvel_se                   &
                + cxm*tmp_vvel_ne + dxt(iw)*tmp_vvel_ee
      tensionsw = -cyp*tmp_uvel_ne + dyt(iw)*tmp_uvel_se                  &
                +  cxp*tmp_vvel_ne - dxt(iw)*tmp_vvel_ee
      shearsw = -cyp*tmp_vvel_ne + dyt(iw)*tmp_vvel_se                    &
              -  cxp*tmp_uvel_ne + dxt(iw)*tmp_uvel_ee
      Deltasw = sqrt(divusw**2 + ecci*(tensionsw**2 + shearsw**2))
  
      ! se
      divuse    = cyp*tmp_uvel_se - dyt(iw)*tmp_uvel_ne                   &
                + cxm*tmp_vvel_se + dxt(iw)*vvel(iw)
      tensionse = -cym*tmp_uvel_se - dyt(iw)*tmp_uvel_ne                  &
                +  cxp*tmp_vvel_se - dxt(iw)*vvel(iw)
      shearse = -cym*tmp_vvel_se - dyt(iw)*tmp_vvel_ne                    &
              -  cxp*tmp_uvel_se + dxt(iw)*uvel(iw)
      Deltase = sqrt(divuse**2 + ecci*(tensionse**2 + shearse**2))
  
    !-----------------------------------------------------------------
    ! replacement pressure/Delta                   ! kg/s
    ! save replacement pressure for principal stress calculation
    !-----------------------------------------------------------------
      c0ne = strength(iw)/max(Deltane,tinyarea)
      c0nw = strength(iw)/max(Deltanw,tinyarea)
      c0sw = strength(iw)/max(Deltasw,tinyarea)
      c0se = strength(iw)/max(Deltase,tinyarea)
  
      c1ne = c0ne*arlx1i
      c1nw = c0nw*arlx1i
      c1sw = c0sw*arlx1i
      c1se = c0se*arlx1i
  
      c0ne = c1ne*ecci
      c0nw = c1nw*ecci
      c0sw = c1sw*ecci
      c0se = c1se*ecci
  
    !-----------------------------------------------------------------
    ! the stresses                            ! kg/s^2
    ! (1) northeast, (2) northwest, (3) southwest, (4) southeast
    !-----------------------------------------------------------------
  
      stressp_1(iw) = (stressp_1(iw) + c1ne*(divune - Deltane)) * denom1
      stressp_2(iw) = (stressp_2(iw) + c1nw*(divunw - Deltanw)) * denom1
      stressp_3(iw) = (stressp_3(iw) + c1sw*(divusw - Deltasw)) * denom1
      stressp_4(iw) = (stressp_4(iw) + c1se*(divuse - Deltase)) * denom1
  
      stressm_1(iw) = (stressm_1(iw) + c0ne*tensionne) * denom1
      stressm_2(iw) = (stressm_2(iw) + c0nw*tensionnw) * denom1
      stressm_3(iw) = (stressm_3(iw) + c0sw*tensionsw) * denom1
      stressm_4(iw) = (stressm_4(iw) + c0se*tensionse) * denom1
      
      stress12_1(iw) = (stress12_1(iw) + c0ne*shearne*p5) * denom1
      stress12_2(iw) = (stress12_2(iw) + c0nw*shearnw*p5) * denom1
      stress12_3(iw) = (stress12_3(iw) + c0sw*shearsw*p5) * denom1
      stress12_4(iw) = (stress12_4(iw) + c0se*shearse*p5) * denom1
  
    !-----------------------------------------------------------------
    ! combinations of the stresses for the momentum equation ! kg/s^2
    !-----------------------------------------------------------------
  
      ssigpn  = stressp_1(iw) + stressp_2(iw)
      ssigps  = stressp_3(iw) + stressp_4(iw)
      ssigpe  = stressp_1(iw) + stressp_4(iw)
      ssigpw  = stressp_2(iw) + stressp_3(iw)
      ssigp1  =(stressp_1(iw) + stressp_3(iw))*p055
      ssigp2  =(stressp_2(iw) + stressp_4(iw))*p055
  
      ssigmn  = stressm_1(iw) + stressm_2(iw)
      ssigms  = stressm_3(iw) + stressm_4(iw)
      ssigme  = stressm_1(iw) + stressm_4(iw)
      ssigmw  = stressm_2(iw) + stressm_3(iw)
      ssigm1  =(stressm_1(iw) + stressm_3(iw))*p055
      ssigm2  =(stressm_2(iw) + stressm_4(iw))*p055
  
      ssig12n = stress12_1(iw) + stress12_2(iw)
      ssig12s = stress12_3(iw) + stress12_4(iw)
      ssig12e = stress12_1(iw) + stress12_4(iw)
      ssig12w = stress12_2(iw) + stress12_3(iw)
      ssig121 =(stress12_1(iw) + stress12_3(iw))*p111
      ssig122 =(stress12_2(iw) + stress12_4(iw))*p111
  
      csigpne = p111*stressp_1(iw) + ssigp2 + p027*stressp_3(iw)
      csigpnw = p111*stressp_2(iw) + ssigp1 + p027*stressp_4(iw)
      csigpsw = p111*stressp_3(iw) + ssigp2 + p027*stressp_1(iw)
      csigpse = p111*stressp_4(iw) + ssigp1 + p027*stressp_2(iw)
      
      csigmne = p111*stressm_1(iw) + ssigm2 + p027*stressm_3(iw)
      csigmnw = p111*stressm_2(iw) + ssigm1 + p027*stressm_4(iw)
      csigmsw = p111*stressm_3(iw) + ssigm2 + p027*stressm_1(iw)
      csigmse = p111*stressm_4(iw) + ssigm1 + p027*stressm_2(iw)
      
      csig12ne = p222*stress12_1(iw) + ssig122 + p055*stress12_3(iw)
      csig12nw = p222*stress12_2(iw) + ssig121 + p055*stress12_4(iw)
      csig12sw = p222*stress12_3(iw) + ssig122 + p055*stress12_1(iw)
      csig12se = p222*stress12_4(iw) + ssig121 + p055*stress12_2(iw)
  
      str12ew = p5*dxt(iw)*(p333*ssig12e + p166*ssig12w)
      str12we = p5*dxt(iw)*(p333*ssig12w + p166*ssig12e)
      str12ns = p5*dyt(iw)*(p333*ssig12n + p166*ssig12s)
      str12sn = p5*dyt(iw)*(p333*ssig12s + p166*ssig12n)
  
    !-----------------------------------------------------------------
    ! for dF/dx (u momentum)
    !-----------------------------------------------------------------
      strp_tmp  = p25*dyt(iw)*(p333*ssigpn  + p166*ssigps)
      strm_tmp  = p25*dyt(iw)*(p333*ssigmn  + p166*ssigms)
  
      ! northeast (iw)
      str1(iw) = -strp_tmp - strm_tmp - str12ew &
                 + dxhy*(-csigpne + csigmne) + dyhx*csig12ne
  
      ! northwest (i+1,j)
      str2(iw) = strp_tmp + strm_tmp - str12we &
                 + dxhy*(-csigpnw + csigmnw) + dyhx*csig12nw
  
      strp_tmp  = p25*dyt(iw)*(p333*ssigps  + p166*ssigpn)
      strm_tmp  = p25*dyt(iw)*(p333*ssigms  + p166*ssigmn)
  
      ! southeast (i,j+1)
      str3(iw) = -strp_tmp - strm_tmp + str12ew &
                 + dxhy*(-csigpse + csigmse) + dyhx*csig12se
  
      ! southwest (i+1,j+1)
      str4(iw) = strp_tmp + strm_tmp + str12we &
                 + dxhy*(-csigpsw + csigmsw) + dyhx*csig12sw
  
    !-----------------------------------------------------------------
    ! for dF/dy (v momentum)
    !-----------------------------------------------------------------
      strp_tmp  = p25*dxt(iw)*(p333*ssigpe  + p166*ssigpw)
      strm_tmp  = p25*dxt(iw)*(p333*ssigme  + p166*ssigmw)
  
      ! northeast (i,j)
      str5(iw) = -strp_tmp + strm_tmp - str12ns &
                 - dyhx*(csigpne + csigmne) + dxhy*csig12ne
  
      ! southeast (i,j+1)
      str6(iw) = strp_tmp - strm_tmp - str12sn &
                 - dyhx*(csigpse + csigmse) + dxhy*csig12se
  
      strp_tmp  = p25*dxt(iw)*(p333*ssigpw  + p166*ssigpe)
      strm_tmp  = p25*dxt(iw)*(p333*ssigmw  + p166*ssigme)
  
      ! northwest (i+1,j)
      str7(iw) = -strp_tmp + strm_tmp + str12ns &
                 - dyhx*(csigpnw + csigmnw) + dxhy*csig12nw
  
      ! southwest (i+1,j+1)
      str8(iw) = strp_tmp - strm_tmp + str12sn &
                 - dyhx*(csigpsw + csigmsw) + dxhy*csig12sw
    enddo   
    !$acc end parallel
  end subroutine stress_i
  
  subroutine stress_l(NA_len, tarear, &
                     ee,ne,se,lb,ub,uvel,vvel,dxt,dyt,                           & 
                     hte,htn,htem1,htnm1,                                        &
                     strength,stressp_1,stressp_2,stressp_3,stressp_4,           & 
                     stressm_1,stressm_2,stressm_3,stressm_4,stress12_1,         &
                     stress12_2,stress12_3,stress12_4,                           &
                     divu,rdg_conv,rdg_shear,shear,                              &
                     str1,str2,str3,str4,str5,str6,str7,str8                     )
    !- modules -------------------------------------------------------------------
    use ice_kinds_mod
    use dmi_omp, only : domp_get_domain
    use ice_constants, only: p027, p055, p111, p166, p222, p25, p333, p5, c1p5, c0
    use icepack_parameters, only: puny
    use ice_dyn_shared, only: ecci, denom1, arlx1i
    !- directives ----------------------------------------------------------------
    implicit none
    ! arguments ------------------------------------------------------------------
    integer (kind=int_kind), intent(in) :: NA_len
    integer (kind=int_kind), intent(in) :: lb,ub
    integer (kind=int_kind), dimension(:), intent(in),    contiguous :: ee,ne,se
    real    (kind=dbl_kind), dimension(:), intent(in),    contiguous ::          &
       strength, uvel, vvel, dxt, dyt, tarear,                                   &
       hte,htn,htem1,htnm1
    real    (kind=dbl_kind), dimension(:), intent(inout), contiguous ::          &
       stressp_1,stressp_2,  stressp_3, stressp_4, stressm_1,  stressm_2,        &
       stressm_3,stressm_4, stress12_1,stress12_2,stress12_3, stress12_4
    real    (kind=dbl_kind), dimension(:), intent(out),   contiguous ::          &
            str1,str2,str3,str4,str5,str6,str7,str8,                             &
            divu,rdg_conv,rdg_shear,shear
    ! local variables ------------------------------------------------------------
    integer (kind=int_kind) :: iw,il,iu
    real    (kind=dbl_kind) ::                                                   &
      divune, divunw, divuse, divusw,tensionne, tensionnw, tensionse, tensionsw, & 
      shearne, shearnw, shearse, shearsw, Deltane, Deltanw, Deltase, Deltasw   , &
      c0ne, c0nw, c0se, c0sw, c1ne, c1nw, c1se, c1sw                           , &
      ssigpn, ssigps, ssigpe, ssigpw, ssigmn, ssigms, ssigme, ssigmw           , &
      ssig12n, ssig12s, ssig12e, ssig12w, ssigp1, ssigp2,ssigm1, ssigm2,ssig121, &
      ssig122, csigpne, csigpnw, csigpse, csigpsw,csigmne, csigmnw, csigmse    , &
      csigmsw, csig12ne, csig12nw, csig12se, csig12sw, str12ew, str12we,str12ns, &
      str12sn, strp_tmp, strm_tmp, tmp_uvel_ee, tmp_vvel_se, tmp_vvel_ee,        &
      tmp_vvel_ne, tmp_uvel_ne, tmp_uvel_se
    real    (kind=dbl_kind) :: dxhy,dyhx,cxp,cyp,cxm,cym,tinyarea
  
#ifdef _OPENACC
    !$acc parallel                                                                 &
    !$acc present(ee,ne,se,strength,uvel,vvel,dxt,dyt,tarear,                      &
    !$acc         hte,htn,htem1,htnm1,                                             &
    !$acc         str1,str2,str3,str4,str5,str6,str7,str8,                         &
    !$acc         stressp_1,stressp_2,stressp_3,stressp_4,                         &
    !$acc         stressm_1,stressm_2,stressm_3,stressm_4,                         &
    !$acc         stress12_1,stress12_2,stress12_3,stress12_4,                     &
    !$acc         divu,rdg_conv,rdg_shear,shear)
    !$acc loop 
    do iw = 1,NA_len
#else
    call domp_get_domain(lb,ub,il,iu)
    do iw = il, iu
#endif
      tinyarea = puny*dxt(iw)*dyt(iw)
      dxhy =     p5*(hte(iw)  - htem1(iw))
      dyhx =     p5*(htn(iw)  - htnm1(iw))
      cxp  =   c1p5*htn(iw)   - p5*htnm1(iw)
      cyp  =   c1p5*hte(iw)   - p5*htem1(iw)
      cxm  = -(c1p5*htnm1(iw) - p5*htn(iw))
      cym  = -(c1p5*htem1(iw) - p5*hte(iw))
  
    !-----------------------------------------------------------------
    ! strain rates
    ! NOTE these are actually strain rates * area  (m^2/s)
    !-----------------------------------------------------------------
       ! divergence  =  e_11 + e_22
       tmp_uvel_ee = uvel(ee(iw))
       tmp_vvel_se = vvel(se(iw))
       tmp_vvel_ee = vvel(ee(iw))
       tmp_vvel_ne = vvel(ne(iw))
       tmp_uvel_ne = uvel(ne(iw))
       tmp_uvel_se = uvel(se(iw))
  
       divune    = cyp*uvel(iw) - dyt(iw)*tmp_uvel_ee                       &
                 + cxp*vvel(iw) - dxt(iw)*tmp_vvel_se
       divunw    = cym*tmp_uvel_ee + dyt(iw)*uvel(iw)                       &
                 + cxp*tmp_vvel_ee - dxt(iw)*tmp_vvel_ne
       divusw    = cym*tmp_uvel_ne + dyt(iw)*tmp_uvel_se                   &
                 + cxm*tmp_vvel_ne + dxt(iw)*tmp_vvel_ee
       divuse    = cyp*tmp_uvel_se - dyt(iw)*tmp_uvel_ne                   &
                 + cxm*tmp_vvel_se + dxt(iw)*vvel(iw)
  
       ! tension strain rate  =  e_11 - e_22
       tensionne = -cym*uvel(iw) - dyt(iw)*tmp_uvel_ee                      &
                 +  cxm*vvel(iw) + dxt(iw)*tmp_vvel_se
       tensionnw = -cyp*tmp_uvel_ee + dyt(iw)*uvel(iw)                      &
                 +  cxm*tmp_vvel_ee + dxt(iw)*tmp_vvel_ne
       tensionsw = -cyp*tmp_uvel_ne + dyt(iw)*tmp_uvel_se                  &
                 +  cxp*tmp_vvel_ne - dxt(iw)*tmp_vvel_ee
       tensionse = -cym*tmp_uvel_se - dyt(iw)*tmp_uvel_ne                  &
                 +  cxp*tmp_vvel_se - dxt(iw)*vvel(iw)
  
       ! shearing strain rate  =  e_12
       shearne = -cym*vvel(iw) - dyt(iw)*tmp_vvel_ee                        &
               -  cxm*uvel(iw) - dxt(iw)*tmp_uvel_se
       shearnw = -cyp*tmp_vvel_ee + dyt(iw)*vvel(iw)                        &
               -  cxm*tmp_uvel_ee - dxt(iw)*tmp_uvel_ne
       shearsw = -cyp*tmp_vvel_ne + dyt(iw)*tmp_vvel_se                    &
               -  cxp*tmp_uvel_ne + dxt(iw)*tmp_uvel_ee
       shearse = -cym*tmp_vvel_se - dyt(iw)*tmp_vvel_ne                    &
               -  cxp*tmp_uvel_se + dxt(iw)*uvel(iw)
       
       ! Delta (in the denominator of zeta, eta)
       Deltane = sqrt(divune**2 + ecci*(tensionne**2 + shearne**2))
       Deltanw = sqrt(divunw**2 + ecci*(tensionnw**2 + shearnw**2))
       Deltase = sqrt(divuse**2 + ecci*(tensionse**2 + shearse**2))
       Deltasw = sqrt(divusw**2 + ecci*(tensionsw**2 + shearsw**2))
  
!       !-----------------------------------------------------------------
!       ! on last subcycle, save quantities for mechanical redistribution
!       !-----------------------------------------------------------------
       divu(iw) = p25*(divune + divunw + divuse + divusw) * tarear(iw)
       rdg_conv(iw)  = -min(divu(iw),c0)                          ! Move outside the entire "kernel"
       rdg_shear(iw) = p5*( p25*(Deltane + Deltanw + Deltase + Deltasw) * tarear(iw) -abs(divu(iw)) )

       ! diagnostic only
       ! shear = sqrt(tension**2 + shearing**2)
       shear(iw) = p25*tarear(iw)*sqrt( &
                 (tensionne + tensionnw + tensionse + tensionsw)**2 &
                +  (shearne +   shearnw +   shearse +   shearsw)**2)

    !-----------------------------------------------------------------
    ! replacement pressure/Delta                   ! kg/s
    ! save replacement pressure for principal stress calculation
    !-----------------------------------------------------------------
       c0ne = strength(iw)/max(Deltane,tinyarea)
       c0nw = strength(iw)/max(Deltanw,tinyarea)
       c0sw = strength(iw)/max(Deltasw,tinyarea)
       c0se = strength(iw)/max(Deltase,tinyarea)
  
       c1ne = c0ne*arlx1i
       c1nw = c0nw*arlx1i
       c1sw = c0sw*arlx1i
       c1se = c0se*arlx1i
  
       c0ne = c1ne*ecci
       c0nw = c1nw*ecci
       c0sw = c1sw*ecci
       c0se = c1se*ecci
  
    !-----------------------------------------------------------------
    ! the stresses                            ! kg/s^2
    ! (1) northeast, (2) northwest, (3) southwest, (4) southeast
    !-----------------------------------------------------------------
  
       stressp_1(iw) = (stressp_1(iw) + c1ne*(divune - Deltane)) * denom1
       stressp_2(iw) = (stressp_2(iw) + c1nw*(divunw - Deltanw)) * denom1
       stressp_3(iw) = (stressp_3(iw) + c1sw*(divusw - Deltasw)) * denom1
       stressp_4(iw) = (stressp_4(iw) + c1se*(divuse - Deltase)) * denom1
  
       stressm_1(iw) = (stressm_1(iw) + c0ne*tensionne) * denom1
       stressm_2(iw) = (stressm_2(iw) + c0nw*tensionnw) * denom1
       stressm_3(iw) = (stressm_3(iw) + c0sw*tensionsw) * denom1
       stressm_4(iw) = (stressm_4(iw) + c0se*tensionse) * denom1
      
       stress12_1(iw) = (stress12_1(iw) + c0ne*shearne*p5) * denom1
       stress12_2(iw) = (stress12_2(iw) + c0nw*shearnw*p5) * denom1
       stress12_3(iw) = (stress12_3(iw) + c0sw*shearsw*p5) * denom1
       stress12_4(iw) = (stress12_4(iw) + c0se*shearse*p5) * denom1
  
    !-----------------------------------------------------------------
    ! combinations of the stresses for the momentum equation ! kg/s^2
    !-----------------------------------------------------------------
  
       ssigpn  = stressp_1(iw) + stressp_2(iw)
       ssigps  = stressp_3(iw) + stressp_4(iw)
       ssigpe  = stressp_1(iw) + stressp_4(iw)
       ssigpw  = stressp_2(iw) + stressp_3(iw)
       ssigp1  =(stressp_1(iw) + stressp_3(iw))*p055
       ssigp2  =(stressp_2(iw) + stressp_4(iw))*p055
  
       ssigmn  = stressm_1(iw) + stressm_2(iw)
       ssigms  = stressm_3(iw) + stressm_4(iw)
       ssigme  = stressm_1(iw) + stressm_4(iw)
       ssigmw  = stressm_2(iw) + stressm_3(iw)
       ssigm1  =(stressm_1(iw) + stressm_3(iw))*p055
       ssigm2  =(stressm_2(iw) + stressm_4(iw))*p055
  
       ssig12n = stress12_1(iw) + stress12_2(iw)
       ssig12s = stress12_3(iw) + stress12_4(iw)
       ssig12e = stress12_1(iw) + stress12_4(iw)
       ssig12w = stress12_2(iw) + stress12_3(iw)
       ssig121 =(stress12_1(iw) + stress12_3(iw))*p111
       ssig122 =(stress12_2(iw) + stress12_4(iw))*p111
  
       csigpne = p111*stressp_1(iw) + ssigp2 + p027*stressp_3(iw)
       csigpnw = p111*stressp_2(iw) + ssigp1 + p027*stressp_4(iw)
       csigpsw = p111*stressp_3(iw) + ssigp2 + p027*stressp_1(iw)
       csigpse = p111*stressp_4(iw) + ssigp1 + p027*stressp_2(iw)
       
       csigmne = p111*stressm_1(iw) + ssigm2 + p027*stressm_3(iw)
       csigmnw = p111*stressm_2(iw) + ssigm1 + p027*stressm_4(iw)
       csigmsw = p111*stressm_3(iw) + ssigm2 + p027*stressm_1(iw)
       csigmse = p111*stressm_4(iw) + ssigm1 + p027*stressm_2(iw)
       
       csig12ne = p222*stress12_1(iw) + ssig122 + p055*stress12_3(iw)
       csig12nw = p222*stress12_2(iw) + ssig121 + p055*stress12_4(iw)
       csig12sw = p222*stress12_3(iw) + ssig122 + p055*stress12_1(iw)
       csig12se = p222*stress12_4(iw) + ssig121 + p055*stress12_2(iw)
  
       str12ew = p5*dxt(iw)*(p333*ssig12e + p166*ssig12w)
       str12we = p5*dxt(iw)*(p333*ssig12w + p166*ssig12e)
       str12ns = p5*dyt(iw)*(p333*ssig12n + p166*ssig12s)
       str12sn = p5*dyt(iw)*(p333*ssig12s + p166*ssig12n)
  
    !-----------------------------------------------------------------
    ! for dF/dx (u momentum)
    !-----------------------------------------------------------------
       strp_tmp  = p25*dyt(iw)*(p333*ssigpn  + p166*ssigps)
       strm_tmp  = p25*dyt(iw)*(p333*ssigmn  + p166*ssigms)
  
       ! northeast (iw)
       str1(iw) = -strp_tmp - strm_tmp - str12ew &
                  + dxhy*(-csigpne + csigmne) + dyhx*csig12ne
  
       ! northwest (i+1,j)
       str2(iw) = strp_tmp + strm_tmp - str12we &
                  + dxhy*(-csigpnw + csigmnw) + dyhx*csig12nw
  
       strp_tmp  = p25*dyt(iw)*(p333*ssigps  + p166*ssigpn)
       strm_tmp  = p25*dyt(iw)*(p333*ssigms  + p166*ssigmn)
  
       ! southeast (i,j+1)
       str3(iw) = -strp_tmp - strm_tmp + str12ew &
                  + dxhy*(-csigpse + csigmse) + dyhx*csig12se
  
       ! southwest (i+1,j+1)
       str4(iw) = strp_tmp + strm_tmp + str12we &
                  + dxhy*(-csigpsw + csigmsw) + dyhx*csig12sw
  
    !-----------------------------------------------------------------
    ! for dF/dy (v momentum)
    !-----------------------------------------------------------------
       strp_tmp  = p25*dxt(iw)*(p333*ssigpe  + p166*ssigpw)
       strm_tmp  = p25*dxt(iw)*(p333*ssigme  + p166*ssigmw)
  
       ! northeast (i,j)
       str5(iw) = -strp_tmp + strm_tmp - str12ns &
                  - dyhx*(csigpne + csigmne) + dxhy*csig12ne
  
       ! southeast (i,j+1)
       str6(iw) = strp_tmp - strm_tmp - str12sn &
                  - dyhx*(csigpse + csigmse) + dxhy*csig12se
  
       strp_tmp  = p25*dxt(iw)*(p333*ssigpw  + p166*ssigpe)
       strm_tmp  = p25*dxt(iw)*(p333*ssigmw  + p166*ssigme)
  
       ! northwest (i+1,j)
       str7(iw) = -strp_tmp + strm_tmp + str12ns &
                  - dyhx*(csigpnw + csigmnw) + dxhy*csig12nw
  
       ! southwest (i+1,j+1)
       str8(iw) = strp_tmp - strm_tmp + str12sn &
                  - dyhx*(csigpsw + csigmsw) + dxhy*csig12sw
    enddo   
    !$acc end parallel
  end subroutine stress_l
  
  subroutine stepu_iter(NA_len,rhow, &
                    lb,ub,Cw,aiu,uocn,vocn,forcex,forcey,                        &
                    umassdti, fm,uarear,uvel_init,vvel_init,uvel,vvel,           &
                    str1,str2,str3,str4,str5,str6,str7,str8, nw,sw,se,skipme)
    !- modules -------------------------------------------------------------------
    use ice_kinds_mod
    use dmi_omp, only : domp_get_domain
    use ice_dyn_shared, only: ecci, denom1, arlx1i, brlx, revp
    use ice_constants, only: c0, c1
    !- directives ----------------------------------------------------------------
    implicit none
    ! arguments ------------------------------------------------------------------
    integer (kind=int_kind), intent(in) :: NA_len
    real (kind=dbl_kind), intent(in) :: rhow
    integer(kind=int_kind),intent(in)   :: lb,ub
    logical(kind=log_kind),intent(in), dimension(:)   :: skipme
    integer(kind=int_kind),dimension(:), intent(in), contiguous :: nw,sw,se
    real(kind=dbl_kind),dimension(:), intent(in), contiguous    ::               &
       uvel_init, vvel_init, aiu, forcex, forcey, umassdti,                      &
       uocn, vocn, fm, uarear,Cw,str1,str2,str3,str4,str5,str6,str7,str8
    real(kind=dbl_kind),dimension(:), intent(inout), contiguous ::               &
       uvel,vvel
    real (kind=dbl_kind), parameter :: &
           cosw = c1   , & ! cos(ocean turning angle)  ! turning angle = 0
           sinw = c0   
    ! local variables
    integer (kind=int_kind) :: iw,il,iu
    real    (kind=dbl_kind) :: uold, vold, vrel,cca,ccb,ab2,cc1,cc2,taux,tauy
    real    (kind=dbl_kind) :: tmp_str2_nw,tmp_str3_se,tmp_str4_sw, tmp_strintx
    real    (kind=dbl_kind) :: tmp_str6_se,tmp_str7_nw,tmp_str8_sw, tmp_strinty
    real    (kind=dbl_kind) :: waterx,watery
  
#ifdef _OPENACC
    !$acc parallel                                                                 &
    !$acc present(Cw,aiu,uocn,vocn,forcex,forcey,umassdti,                         &
    !$acc         fm,uarear,uvel_init,vvel_init,nw,sw,se,skipme,                   &
    !$acc         str1,str2,str3,str4,str5,str6,str7,str8,uvel,vvel)
    !$acc loop 
    do iw = 1,NA_len
#else
    call domp_get_domain(lb,ub,il,iu)
    do iw = il, iu
#endif
       if (skipme(iw)) cycle
       uold = uvel(iw)
       vold = vvel(iw)
       vrel = aiu(iw)*rhow*Cw(iw)*sqrt((uocn(iw)-uold)**2+(vocn(iw)-vold)**2)
       waterx = uocn(iw)*cosw - vocn(iw)*sinw*sign(c1,fm(iw))
       watery = vocn(iw)*cosw + uocn(iw)*sinw*sign(c1,fm(iw))
       taux = vrel*waterx 
       tauy = vrel*watery
       cca = (brlx + revp)*umassdti(iw) + vrel * cosw 
       ccb = fm(iw) + sign(c1,fm(iw)) * vrel * sinw 
       ab2 = cca**2 + ccb**2
       ! southeast(i,j+1)   = se
       ! northwest(i+1,j)   = nw
       ! southwest(i+1,j+1) = sw
       tmp_str2_nw = str2(nw(iw))
       tmp_str3_se = str3(se(iw))
       tmp_str4_sw = str4(sw(iw))
       tmp_str6_se = str6(se(iw))
       tmp_str7_nw = str7(nw(iw))
       tmp_str8_sw = str8(sw(iw))
  
       tmp_strintx = uarear(iw)*(str1(iw)+tmp_str2_nw+tmp_str3_se+tmp_str4_sw)
       tmp_strinty = uarear(iw)*(str5(iw)+tmp_str6_se+tmp_str7_nw+tmp_str8_sw)
       cc1 = tmp_strintx + forcex(iw) + taux &
           + umassdti(iw)*(brlx*uold + revp*uvel_init(iw))
       cc2 = tmp_strinty + forcey(iw) + tauy &
           + umassdti(iw)*(brlx*vold + revp*vvel_init(iw))
       uvel(iw) = (cca*cc1 + ccb*cc2) / ab2 
       vvel(iw) = (cca*cc2 - ccb*cc1) / ab2
    enddo
    !$acc end parallel

  end subroutine stepu_iter
  
  subroutine stepu_last(NA_len, rhow, &
                    lb,ub,Cw,aiu,uocn,vocn,forcex,forcey, umassdti,fm,uarear,    &
                    strocnx,strocny,strintx,strinty,                             &
                    uvel_init,vvel_init,uvel,vvel,                               &
                    str1,str2,str3,str4,str5,str6,str7,str8, nw,sw,se,skipme)
    !- modules -------------------------------------------------------------------
    use ice_kinds_mod
    use dmi_omp, only : domp_get_domain
    use ice_constants, only: c0, c1
    use icepack_parameters, only: puny
    use ice_dyn_shared, only: ecci, denom1, arlx1i, brlx, revp
    !- directives ----------------------------------------------------------------
    implicit none
    ! arguments ------------------------------------------------------------------
    integer (kind=int_kind), intent(in) :: NA_len
    real (kind=dbl_kind), intent(in) :: rhow
    logical(kind=log_kind),intent(in), dimension(:)   :: skipme
    integer(kind=int_kind),intent(in)   :: lb,ub
    integer(kind=int_kind),dimension(:), intent(in), contiguous :: nw,sw,se
    real(kind=dbl_kind),dimension(:), intent(in), contiguous    ::               &
       uvel_init, vvel_init, aiu, forcex, forcey, umassdti,      &
       uocn, vocn, fm, uarear,Cw,str1,str2,str3,str4,str5,str6,str7,str8
    real(kind=dbl_kind),dimension(:), intent(inout), contiguous ::               &
       uvel,vvel,strintx, strinty
    real(kind=dbl_kind),dimension(:), intent(out),   contiguous ::               &
       strocnx,strocny
    real (kind=dbl_kind), parameter :: &
           cosw = c1   , & ! cos(ocean turning angle)  ! turning angle = 0
           sinw = c0   
    ! local variables
    integer (kind=int_kind) :: iw,il,iu
    real    (kind=dbl_kind) :: uold, vold, vrel,cca,ccb,ab2,cc1,cc2,taux,tauy
    real    (kind=dbl_kind) :: tmp_str2_nw,tmp_str3_se,tmp_str4_sw
    real    (kind=dbl_kind) :: tmp_str6_se,tmp_str7_nw,tmp_str8_sw
    real    (kind=dbl_kind) :: waterx,watery
   
#ifdef _OPENACC
    !$acc parallel                                                                 &
    !$acc present(Cw,aiu,uocn,vocn,forcex,forcey,umassdti,                         &
    !$acc         fm,uarear,strintx,strinty,uvel_init,vvel_init,nw,sw,se,skipme,   &
    !$acc         str1,str2,str3,str4,str5,str6,str7,str8,strocnx,strocny,uvel,vvel)
    !$acc loop 
    do iw = 1,NA_len
#else
    call domp_get_domain(lb,ub,il,iu)
    do iw = il, iu
#endif
       if (skipme(iw)) cycle
       uold = uvel(iw)
       vold = vvel(iw)
       vrel = aiu(iw)*rhow*Cw(iw)*sqrt((uocn(iw)-uold)**2+(vocn(iw)-vold)**2)
       waterx = uocn(iw)*cosw - vocn(iw)*sinw*sign(c1,fm(iw))
       watery = vocn(iw)*cosw + uocn(iw)*sinw*sign(c1,fm(iw))
       taux = vrel*waterx 
       tauy = vrel*watery
       cca = (brlx + revp)*umassdti(iw) + vrel * cosw 
       ccb = fm(iw) + sign(c1,fm(iw)) * vrel * sinw 
       ab2 = cca**2 + ccb**2
       ! southeast(i,j+1)   = se
       ! northwest(i+1,j)   = nw
       ! southwest(i+1,j+1) = sw
       tmp_str2_nw = str2(nw(iw))
       tmp_str3_se = str3(se(iw))
       tmp_str4_sw = str4(sw(iw))
       tmp_str6_se = str6(se(iw))
       tmp_str7_nw = str7(nw(iw))
       tmp_str8_sw = str8(sw(iw))
  
       strintx(iw) = uarear(iw)*(str1(iw)+tmp_str2_nw+tmp_str3_se+tmp_str4_sw)
       strinty(iw) = uarear(iw)*(str5(iw)+tmp_str6_se+tmp_str7_nw+tmp_str8_sw)
       cc1 = strintx(iw) + forcex(iw) + taux &
           + umassdti(iw)*(brlx*uold + revp*uvel_init(iw))
       cc2 = strinty(iw) + forcey(iw) + tauy &
           + umassdti(iw)*(brlx*vold + revp*vvel_init(iw))
       uvel(iw) = (cca*cc1 + ccb*cc2) / ab2 
       vvel(iw) = (cca*cc2 - ccb*cc1) / ab2
       strocnx(iw) = taux
       strocny(iw) = tauy
    enddo
    !$acc end parallel
  end subroutine stepu_last

  subroutine halo_update(NAVEL_len,lb,ub,uvel,vvel, halo_parent)
    !- modules -------------------------------------------------------------------
    use ice_kinds_mod
    use dmi_omp, only : domp_get_domain
    !- directives ----------------------------------------------------------------
    implicit none
    ! arguments ------------------------------------------------------------------
    integer (kind=int_kind), intent(in) :: NAVEL_len
    integer(kind=int_kind),intent(in)   :: lb,ub
    integer(kind=int_kind),dimension(:), intent(in), contiguous :: halo_parent
    real(kind=dbl_kind),dimension(:), intent(inout), contiguous :: uvel,vvel
    ! local variables
    integer (kind=int_kind) :: iw,il,iu

#ifdef _OPENACC
    !$acc parallel                                   &
    !$acc present(uvel,vvel)                         &
    !$acc loop
    do iw = 1,NAVEL_len
#else
    call domp_get_domain(lb,ub,il,iu)
    do iw = il, iu
#endif
       if (halo_parent(iw)==0) cycle
       uvel(iw) = uvel(halo_parent(iw))
       vvel(iw) = vvel(halo_parent(iw))
    enddo
    !$acc end parallel
  end subroutine halo_update

end module bench_v2
  
!===============================================================================
  
!-- One dimension representation of EVP 2D arrays used for EVP kernels
module evp_kernel1d
  use ice_kinds_mod
  !-- BEGIN: specific for the KERNEL
  use ice_dyn_shared, only: revp, ecci, denom1, arlx1i, brlx
  !-- END: specific for the KERNEL
  implicit none
  !- interfaces
  interface evp_copyin
    module procedure evp_copyin_v0
!    module procedure evp_copyin_v1
!    module procedure evp_copyin_v2
  end interface 
  public :: evp_copyin, evp_copyout, evp_kernel_v1, evp_kernel_v2
  private
  save
  integer(kind=int_kind) ::                                        &
    NA_len, NAVEL_len
  logical(kind=log_kind), dimension(:), allocatable ::             &
    skipucell
  integer(kind=int_kind), dimension(:), allocatable ::             &
    ee,ne,se,nw,sw,sse,indi,indj,indij , halo_parent                  
  real (kind=dbl_kind), dimension(:), allocatable ::               &
    cdn_ocn,aiu,uocn,vocn,waterx,watery,forcex,forcey,tarear,      &
    umassdti,fm,uarear,strintx,strinty,uvel_init,vvel_init
  real (kind=dbl_kind), dimension(:), allocatable ::               &
    strength,uvel,vvel,dxt,dyt,dxhy,dyhx,cyp,cxp,cym,cxm,tinyarea, &
    stressp_1, stressp_2, stressp_3, stressp_4,                    &
    stressm_1, stressm_2, stressm_3, stressm_4,                    &
    stress12_1,stress12_2,stress12_3,stress12_4,                   &
    str1, str2, str3, str4, str5, str6, str7, str8,                &
    divu,rdg_conv,rdg_shear,shear
  real (kind=dbl_kind), dimension(:), allocatable ::               &
    HTE,HTN,                                                       &
    HTEm1,HTNm1
  real (kind=dbl_kind), dimension(:), allocatable ::               &
    strocnx,strocny
  contains
  subroutine alloc1d(na)
    implicit none
    integer(kind=int_kind),intent(in) :: na
    integer(kind=int_kind)            :: ierr,nb
    nb=na
    allocate(                                                      &
      ! U+T cells
      ! Helper index for neighbours
      indj(1:na),indi(1:na),                                       &
      ee(1:na),ne(1:na),se(1:na),                                  &
      nw(1:nb),sw(1:nb),sse(1:nb),                                 &
      skipucell(1:na),                             &
      ! Grid distances: HTE,HTN + "-1 neighbours" 
      HTE(1:na),HTN(1:na),                                         &
      HTEm1(1:na),HTNm1(1:na),                                     &
      ! T cells
      strength(1:na),dxt(1:na),dyt(1:na),dxhy(1:na),dyhx(1:na),           &
      cyp(1:na),cxp(1:na),cym(1:na),cxm(1:na),tinyarea(1:na),tarear(1:na),&
      stressp_1(1:na), stressp_2(1:na), stressp_3(1:na), stressp_4(1:na), &
      stressm_1(1:na), stressm_2(1:na), stressm_3(1:na), stressm_4(1:na), &
      stress12_1(1:na),stress12_2(1:na),stress12_3(1:na),stress12_4(1:na),&
      divu(1:na),rdg_conv(1:na),rdg_shear(1:na),shear(1:na),              &
      ! U cells
      cdn_ocn(1:nb),aiu(1:nb),uocn(1:nb),vocn(1:nb),                      &
      waterx(1:nb),watery(1:nb),forcex(1:nb),forcey(1:nb),                &
      umassdti(1:nb),fm(1:nb),uarear(1:nb),                               &
      strocnx(1:nb),strocny(1:nb),strintx(1:nb),strinty(1:nb),            &
      uvel_init(1:nb),vvel_init(1:nb),                                    &
      stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D'
  end subroutine alloc1d
  subroutine alloc1d_navel(navel)
    implicit none
    integer(kind=int_kind),intent(in) :: navel
    integer(kind=int_kind)            :: ierr
    allocate(                                                            &
      uvel(1:navel),vvel(1:navel), indij(1:navel), halo_parent(1:navel), &
      str1(1:navel),str2(1:navel),str3(1:navel),str4(1:navel),           &
      str5(1:navel),str6(1:navel),str7(1:navel),str8(1:navel),           &
      stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D navel'
  end subroutine alloc1d_navel
  subroutine dealloc1d
    implicit none
    integer(kind=int_kind) :: ierr
    deallocate(                                   &
      ! U+T cells
      ! Helper index for neighbours
      indj,indi,                                  &
      ee,ne,se,                                   &
      nw,sw,sse,                                  &
      skipucell,                                  &
      ! T cells
      strength,dxt,dyt,tarear,                    &
      stressp_1, stressp_2, stressp_3, stressp_4, &
      stressm_1, stressm_2, stressm_3, stressm_4, &
      stress12_1,stress12_2,stress12_3,stress12_4,&
      str1, str2,str3,str4,                       &
      str5, str6,str7,str8,                       &
      divu,rdg_conv,rdg_shear,shear,              &
      ! U cells
      cdn_ocn,aiu,uocn,vocn,                      &
      waterx,watery,forcex,forcey,                &
      umassdti,fm,uarear,                         &
      strocnx,strocny,strintx,strinty,            &
      uvel_init,vvel_init,                        &
      ! NAVEL 
      uvel,vvel, indij, halo_parent,              &
      stat=ierr)
    if (ierr/=0) stop 'Error de-allocating 1D'
    if (allocated(tinyarea)) then
      deallocate(                                 &
        dxhy,dyhx,cyp,cxp,cym,cxm,tinyarea,       &
        stat=ierr)
      if (ierr/=0) stop 'Error de-allocating 1D, v1'
    endif
    if (allocated(HTE)) then
      deallocate(                                 &
        ! Grid distances: HTE,HTN + "-1 neighbours" 
        HTE,HTN, HTEm1,HTNm1,                     &
        stat=ierr)
      if (ierr/=0) stop 'Error de-allocating 1D, v2'
    endif
  end subroutine dealloc1d
!===============================================================================
!===============================================================================
  subroutine evp_copyin_v0(nx,ny,nblk,nx_glob,ny_glob,                                        &
       I_HTE,I_HTN,I_dxhy,I_dyhx,I_cyp,I_cxp,I_cym,I_cxm,I_tinyarea,                          &
       I_icetmask,I_iceumask,                                                                 &
       I_cdn_ocn,I_aiu,I_uocn,I_vocn,I_waterx,I_watery,I_forcex,I_forcey,                     &
       I_umassdti,I_fm,I_uarear,I_tarear,I_strintx,I_strinty,I_uvel_init,I_vvel_init,         &
       I_strength,I_uvel,I_vvel,I_dxt,I_dyt,                                                  &
        I_stressp_1 ,I_stressp_2, I_stressp_3, I_stressp_4,                                   &
        I_stressm_1 ,I_stressm_2, I_stressm_3, I_stressm_4,                                   &
       I_stress12_1,I_stress12_2,I_stress12_3,I_stress12_4                                    )
    use ice_gather_scatter, only: gather_global_ext
    use ice_domain, only: distrb_info
    use ice_communicate, only: my_task, master_task
    use ice_constants, only: c0,c1,p5
    implicit none
    integer(int_kind), intent(in) :: nx, ny, nblk, nx_glob, ny_glob
    integer (kind=int_kind),dimension (nx,ny,nblk), intent(in) :: I_icetmask
    logical (kind=log_kind),dimension (nx,ny,nblk), intent(in) :: I_iceumask
    real (kind=dbl_kind), dimension(nx,ny,nblk), intent(in)  ::                               &
       I_HTE,I_HTN,I_dxhy,I_dyhx,I_cyp,I_cxp,I_cym,I_cxm,I_tinyarea,                          &
       I_cdn_ocn,I_aiu,I_uocn,I_vocn,I_waterx,I_watery,I_forcex,I_forcey,                     &
       I_umassdti,I_fm,I_uarear,I_tarear,I_strintx,I_strinty,I_uvel_init,I_vvel_init,         &
       I_strength,I_uvel,I_vvel,I_dxt,I_dyt,                                                  &
       I_stressp_1, I_stressp_2, I_stressp_3, I_stressp_4,                                    &
       I_stressm_1, I_stressm_2, I_stressm_3, I_stressm_4,                                    &
       I_stress12_1,I_stress12_2,I_stress12_3,I_stress12_4
    ! local variables
    integer (kind=int_kind),dimension (nx_glob,ny_glob) :: G_icetmask
    logical (kind=log_kind),dimension (nx_glob,ny_glob) :: G_iceumask
    real (kind=dbl_kind), dimension(nx_glob,ny_glob) ::                                       &
       G_HTE,G_HTN,G_dxhy,G_dyhx,G_cyp,G_cxp,G_cym,G_cxm,G_tinyarea,                          &
       G_cdn_ocn,G_aiu,G_uocn,G_vocn,G_waterx,G_watery,G_forcex,G_forcey,                     &
       G_umassdti,G_fm,G_uarear,G_tarear,G_strintx,G_strinty,G_uvel_init,G_vvel_init,         &
       G_strength,G_uvel,G_vvel,G_dxt,G_dyt,                                                  &
       G_stressp_1, G_stressp_2, G_stressp_3, G_stressp_4,                                    &
       G_stressm_1, G_stressm_2, G_stressm_3, G_stressm_4,                                    &
       G_stress12_1,G_stress12_2,G_stress12_3,G_stress12_4
    ! Temporary variables
    real (kind=dbl_kind), dimension(nx,ny,nblk) :: I_work
    real (kind=dbl_kind), dimension(nx_glob,ny_glob) :: G_work
    integer(int_kind) :: na, navel, i,j,n
    !-- Gather data into one single block --
    ! BEGIN: Gather data
    !-- icetmask: gather_global_ext only works for reals
    !$OMP PARALLEL DO PRIVATE(n,i,j)
    do n=1,nblk
      do j=1,ny
        do i=1,nx
          if (I_icetmask(i,j,n) == 1) then
            I_work(i,j,n)=c1
          else
            I_work(i,j,n)=c0
          endif
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    call gather_global_ext(G_work, I_work, master_task, distrb_info)
    !$OMP PARALLEL DO PRIVATE(i,j)
    do j=1,ny_glob
      do i=1,nx_glob
        if (G_work(i,j)>p5)  then
          G_icetmask(i,j)=1
        else
          G_icetmask(i,j)=0
        endif
      enddo
    enddo
    !$OMP END PARALLEL DO
    !-- iceumask: gather_global_ext only works for reals
    !$OMP PARALLEL DO PRIVATE(n,i,j)
    do n=1,nblk
      do j=1,ny
        do i=1,nx
          if (I_iceumask(i,j,n) .eqv. .true.) then
            I_work(i,j,n)=c1
          else
            I_work(i,j,n)=c0
          endif
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    call gather_global_ext(G_work, I_work, master_task, distrb_info)
    !$OMP PARALLEL DO PRIVATE(i,j)
    do j=1,ny_glob
      do i=1,nx_glob
        if (G_work(i,j)>p5)  then
          G_iceumask(i,j)=.true.
        else
          G_iceumask(i,j)=.false.
        endif
      enddo
    enddo
    !$OMP END PARALLEL DO
    !-- the rest
    call gather_global_ext(G_HTE, I_HTE, master_task, distrb_info)
    call gather_global_ext(G_HTN, I_HTN, master_task, distrb_info)
    call gather_global_ext(G_dxhy, I_dxhy, master_task, distrb_info)
    call gather_global_ext(G_dyhx, I_dyhx, master_task, distrb_info)
    call gather_global_ext(G_cyp, I_cyp, master_task, distrb_info)
    call gather_global_ext(G_cxp, I_cxp, master_task, distrb_info)
    call gather_global_ext(G_cym, I_cym, master_task, distrb_info)
    call gather_global_ext(G_cxm, I_cxm, master_task, distrb_info)
    call gather_global_ext(G_tinyarea, I_tinyarea, master_task, distrb_info)
    call gather_global_ext(G_cdn_ocn, I_cdn_ocn, master_task, distrb_info)
    call gather_global_ext(G_aiu, I_aiu, master_task, distrb_info)
    call gather_global_ext(G_uocn, I_uocn, master_task, distrb_info)
    call gather_global_ext(G_vocn, I_vocn, master_task, distrb_info)
    call gather_global_ext(G_waterx, I_waterx, master_task, distrb_info)
    call gather_global_ext(G_watery, I_watery, master_task, distrb_info)
    call gather_global_ext(G_forcex, I_forcex, master_task, distrb_info)
    call gather_global_ext(G_forcey, I_forcey, master_task, distrb_info)
    call gather_global_ext(G_umassdti, I_umassdti, master_task, distrb_info)
    call gather_global_ext(G_fm, I_fm, master_task, distrb_info)
    call gather_global_ext(G_uarear, I_uarear, master_task, distrb_info)
    call gather_global_ext(G_tarear, I_tarear, master_task, distrb_info)
    call gather_global_ext(G_strintx, I_strintx, master_task, distrb_info)
    call gather_global_ext(G_strinty, I_strinty, master_task, distrb_info)
    call gather_global_ext(G_uvel_init, I_uvel_init, master_task, distrb_info)
    call gather_global_ext(G_vvel_init, I_vvel_init, master_task, distrb_info)
    call gather_global_ext(G_strength, I_strength, master_task, distrb_info)
    call gather_global_ext(G_uvel, I_uvel, master_task, distrb_info)
    call gather_global_ext(G_vvel, I_vvel, master_task, distrb_info)
    call gather_global_ext(G_dxt, I_dxt, master_task, distrb_info)
    call gather_global_ext(G_dyt, I_dyt, master_task, distrb_info)
    call gather_global_ext(G_stressp_1, I_stressp_1, master_task, distrb_info)
    call gather_global_ext(G_stressp_2, I_stressp_2, master_task, distrb_info)
    call gather_global_ext(G_stressp_3, I_stressp_3, master_task, distrb_info)
    call gather_global_ext(G_stressp_4, I_stressp_4, master_task, distrb_info)
    call gather_global_ext(G_stressm_1, I_stressm_1, master_task, distrb_info)
    call gather_global_ext(G_stressm_2, I_stressm_2, master_task, distrb_info)
    call gather_global_ext(G_stressm_3, I_stressm_3, master_task, distrb_info)
    call gather_global_ext(G_stressm_4, I_stressm_4, master_task, distrb_info)
    call gather_global_ext(G_stress12_1, I_stress12_1, master_task, distrb_info)
    call gather_global_ext(G_stress12_2, I_stress12_2, master_task, distrb_info)
    call gather_global_ext(G_stress12_3, I_stress12_3, master_task, distrb_info)
    call gather_global_ext(G_stress12_4, I_stress12_4, master_task, distrb_info)
    ! END: Gather data

    !-- Find number of active points and allocate vectors --
    call calc_na(nx_glob,ny_glob,na,G_icetmask)
    call alloc1d(na)
    call calc_2d_indices(nx_glob,ny_glob,na, G_icetmask, G_iceumask)
    call calc_navel(nx_glob,ny_glob,na,navel)
    call alloc1d_navel(navel)
!MHRI    !$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,na,navel)
!MHRI    !$OMP END PARALLEL
    ! Remap 2d to 1d and fill in
    call convert_2d_1d(nx_glob,ny_glob,na,navel,                                        &
       G_HTE,G_HTN,G_dxhy,G_dyhx,G_cyp,G_cxp,G_cym,G_cxm,G_tinyarea,                    &
       G_cdn_ocn,G_aiu,G_uocn,G_vocn,G_waterx,G_watery,G_forcex,G_forcey,               &
       G_umassdti,G_fm,G_uarear,G_tarear,G_strintx,G_strinty,G_uvel_init,G_vvel_init,   &
       G_strength,G_uvel,G_vvel,G_dxt,G_dyt,                                            &
        G_stressp_1, G_stressp_2, G_stressp_3, G_stressp_4,                             &
        G_stressm_1, G_stressm_2, G_stressm_3, G_stressm_4,                             &
       G_stress12_1,G_stress12_2,G_stress12_3,G_stress12_4                              )
    call calc_halo_parent(nx_glob,ny_glob,na,navel, G_icetmask)
    NA_len=na
    NAVEL_len=navel
    !-- write check
!if (1 == 1) then
!  write(*,*)'MHRI: INDICES start: evp-copyin'
!  write(*,*) 'na,navel  ', na,navel
!  write(*,*) 'Min/max ee', minval(ee(1:na)), maxval(ee(1:na))
!  write(*,*) 'Min/max ne', minval(ne(1:na)), maxval(ne(1:na))
!  write(*,*) 'Min/max se', minval(se(1:na)), maxval(se(1:na))
!  write(*,*) 'Min/max nw', minval(nw(1:na)), maxval(nw(1:na))
!  write(*,*) 'Min/max sw', minval(sw(1:na)), maxval(sw(1:na))
!  write(*,*) 'Min/max sse', minval(sse(1:na)), maxval(sse(1:na))
!  write(*,*)'MHRI: INDICES end: evp-copyin'
!endif
  end subroutine evp_copyin_v0
  !===============================================================================
  subroutine evp_copyout(nx,ny,nblk,nx_glob,ny_glob,                    &
               I_uvel,I_vvel, I_strintx,I_strinty, I_strocnx,I_strocny, &
               I_stressp_1, I_stressp_2, I_stressp_3, I_stressp_4,      &
               I_stressm_1, I_stressm_2, I_stressm_3, I_stressm_4,      &
               I_stress12_1,I_stress12_2,I_stress12_3,I_stress12_4,     &
               I_divu,I_rdg_conv,I_rdg_shear,I_shear                    )
    use ice_constants, only : c0, field_loc_center, field_loc_NEcorner, &
                                  field_type_scalar, field_type_vector
    use ice_gather_scatter, only: scatter_global_ext, scatter_global
    use ice_domain, only: distrb_info
    use ice_communicate, only: my_task, master_task
    implicit none
    integer(int_kind), intent(in) :: nx,ny,nblk, nx_glob,ny_glob
    real(dbl_kind), dimension(nx,ny,nblk), intent(out) :: &
       I_uvel,I_vvel, I_strintx,I_strinty, I_strocnx,I_strocny, &
       I_stressp_1, I_stressp_2, I_stressp_3, I_stressp_4,      &
       I_stressm_1, I_stressm_2, I_stressm_3, I_stressm_4,      &
       I_stress12_1,I_stress12_2,I_stress12_3,I_stress12_4,     &
       I_divu,I_rdg_conv, I_rdg_shear,I_shear
    ! local variables
    real(dbl_kind), dimension(nx_glob,ny_glob) :: &
       G_uvel,G_vvel, G_strintx,G_strinty, G_strocnx,G_strocny, &
       G_stressp_1, G_stressp_2, G_stressp_3, G_stressp_4,      &
       G_stressm_1, G_stressm_2, G_stressm_3, G_stressm_4,      &
       G_stress12_1,G_stress12_2,G_stress12_3,G_stress12_4,     &
       G_divu,G_rdg_conv, G_rdg_shear,G_shear
    integer(int_kind) :: i,j,iw
    ! Remap 1d to 2d and fill in
    G_uvel       = c0
    G_vvel       = c0
    G_strintx    = c0
    G_strinty    = c0
    G_strocnx    = c0
    G_strocny    = c0
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
    G_divu       = c0
    G_rdg_conv   = c0
    G_rdg_shear  = c0
    G_shear      = c0
    !$OMP PARALLEL PRIVATE(iw,i,j)
    do iw=1,NA_len
      i=indi(iw)
      j=indj(iw)
      G_uvel(i,j)       = uvel(iw)
      G_vvel(i,j)       = vvel(iw)
      G_strintx(i,j)    = strintx(iw)
      G_strinty(i,j)    = strinty(iw)
      G_strocnx(i,j)    = strocnx(iw)
      G_strocny(i,j)    = strocny(iw)
      G_stressp_1(i,j)  = stressp_1(iw)
      G_stressp_2(i,j)  = stressp_2(iw)
      G_stressp_3(i,j)  = stressp_3(iw)
      G_stressp_4(i,j)  = stressp_4(iw)
      G_stressm_1(i,j)  = stressm_1(iw)
      G_stressm_2(i,j)  = stressm_2(iw)
      G_stressm_3(i,j)  = stressm_3(iw)
      G_stressm_4(i,j)  = stressm_4(iw)
      G_stress12_1(i,j) = stress12_1(iw)
      G_stress12_2(i,j) = stress12_2(iw)
      G_stress12_3(i,j) = stress12_3(iw)
      G_stress12_4(i,j) = stress12_4(iw)
      G_divu(i,j)       = divu(iw)
      G_rdg_conv(i,j)   = rdg_conv(iw)
      G_rdg_shear(i,j)  = rdg_shear(iw)
      G_shear(i,j)      = shear(iw)
    enddo
    !$OMP END PARALLEL
    call dealloc1d()
    !-- Scatter data into blocks --
    ! BEGIN: Scatter data
    call scatter_global_ext(I_uvel, G_uvel, master_task, distrb_info)
    call scatter_global_ext(I_vvel, G_vvel, master_task, distrb_info)
    call scatter_global_ext(I_strintx, G_strintx, master_task, distrb_info)
    call scatter_global_ext(I_strinty, G_strinty, master_task, distrb_info)
    call scatter_global_ext(I_strocnx, G_strocnx, master_task, distrb_info)
    call scatter_global_ext(I_strocny, G_strocny, master_task, distrb_info)
    call scatter_global_ext(I_stressp_1, G_stressp_1, master_task, distrb_info)
    call scatter_global_ext(I_stressp_2, G_stressp_2, master_task, distrb_info)
    call scatter_global_ext(I_stressp_3, G_stressp_3, master_task, distrb_info)
    call scatter_global_ext(I_stressp_4, G_stressp_4, master_task, distrb_info)
    call scatter_global_ext(I_stressm_1, G_stressm_1, master_task, distrb_info)
    call scatter_global_ext(I_stressm_2, G_stressm_2, master_task, distrb_info)
    call scatter_global_ext(I_stressm_3, G_stressm_3, master_task, distrb_info)
    call scatter_global_ext(I_stressm_4, G_stressm_4, master_task, distrb_info)
    call scatter_global_ext(I_stress12_1, G_stress12_1, master_task, distrb_info)
    call scatter_global_ext(I_stress12_2, G_stress12_2, master_task, distrb_info)
    call scatter_global_ext(I_stress12_3, G_stress12_3, master_task, distrb_info)
    call scatter_global_ext(I_stress12_4, G_stress12_4, master_task, distrb_info)
    call scatter_global_ext(I_divu, G_divu, master_task, distrb_info)
    call scatter_global_ext(I_rdg_conv, G_rdg_conv, master_task, distrb_info)
    call scatter_global_ext(I_rdg_shear, G_rdg_shear, master_task, distrb_info)
    call scatter_global_ext(I_shear, G_shear, master_task, distrb_info)

! -- Do not work ---
!    call scatter_global(I_uvel, G_uvel, master_task, distrb_info, field_loc_NEcorner, field_type_vector)
!    call scatter_global(I_vvel, G_vvel, master_task, distrb_info, field_loc_NEcorner, field_type_vector)
!    call scatter_global(I_strintx, G_strintx, master_task, distrb_info, field_loc_center, field_type_scalar)
!    call scatter_global(I_strinty, G_strinty, master_task, distrb_info, field_loc_center, field_type_scalar)
!    call scatter_global(I_strocnx, G_strocnx, master_task, distrb_info, field_loc_center, field_type_scalar)
!    call scatter_global(I_strocny, G_strocny, master_task, distrb_info, field_loc_center, field_type_scalar)
!    call scatter_global(I_stressp_1, G_stressp_1, master_task, distrb_info, field_loc_center, field_type_scalar)
!    call scatter_global(I_stressp_2, G_stressp_2, master_task, distrb_info, field_loc_center, field_type_scalar)
!    call scatter_global(I_stressp_3, G_stressp_3, master_task, distrb_info, field_loc_center, field_type_scalar)
!    call scatter_global(I_stressp_4, G_stressp_4, master_task, distrb_info, field_loc_center, field_type_scalar)
!    call scatter_global(I_stressm_1, G_stressm_1, master_task, distrb_info, field_loc_center, field_type_scalar)
!    call scatter_global(I_stressm_2, G_stressm_2, master_task, distrb_info, field_loc_center, field_type_scalar)
!    call scatter_global(I_stressm_3, G_stressm_3, master_task, distrb_info, field_loc_center, field_type_scalar)
!    call scatter_global(I_stressm_4, G_stressm_4, master_task, distrb_info, field_loc_center, field_type_scalar)
!    call scatter_global(I_stress12_1, G_stress12_1, master_task, distrb_info, field_loc_center, field_type_scalar)
!    call scatter_global(I_stress12_2, G_stress12_2, master_task, distrb_info, field_loc_center, field_type_scalar)
!    call scatter_global(I_stress12_3, G_stress12_3, master_task, distrb_info, field_loc_center, field_type_scalar)
!    call scatter_global(I_stress12_4, G_stress12_4, master_task, distrb_info, field_loc_center, field_type_scalar)
!    call scatter_global(I_divu, G_divu, master_task, distrb_info, field_loc_center, field_type_scalar)
!    call scatter_global(I_rdg_conv, G_rdg_conv, master_task, distrb_info, field_loc_center, field_type_scalar)
!    call scatter_global(I_rdg_shear, G_rdg_shear, master_task, distrb_info, field_loc_center, field_type_scalar)
!    call scatter_global(I_shear, G_shear, master_task, distrb_info, field_loc_center, field_type_scalar)
! -- Do not work ---

    ! END: Scatter data
  end subroutine evp_copyout
  !===============================================================================
  subroutine evp_kernel_v1
    use ice_constants, only : c0
    use ice_dyn_shared, only: ndte
    use bench_v1, only : stress, stepu, halo_update
    use dmi_omp, only : domp_init
    use icepack_intfc, only: icepack_query_parameters
    implicit none
    real(kind=dbl_kind) :: rhow
    integer (kind=int_kind) :: ierr, lun, i, nthreads
    integer (kind=int_kind) :: na,nb,navel
    !- Read constants...
    call icepack_query_parameters(rhow_out=rhow)
    na=NA_len
    nb=NA_len
    navel=NAVEL_len

    !- Initialize openmp ---------------------------------------------------------
    call domp_init(nthreads) ! ought to be called from main
  
    !- Initialize timers ---------------------------------------------------------
    str1=c0
    str2=c0
    str3=c0
    str4=c0
    str5=c0
    str6=c0
    str7=c0
    str8=c0
   
    if (ndte<2) STOP 'ndte must be 2 or higher for this kernel'
    !$OMP PARALLEL PRIVATE(i)
    do i = 1, ndte-1
      call stress (NA_len, &
                   ee,ne,se,1,na,uvel,vvel,dxt,dyt,                              & 
                   tinyarea,dxhy,dyhx,cxp,cyp,cxm,cym,                           &
                   strength,stressp_1,stressp_2,stressp_3,stressp_4,             & 
                   stressm_1,stressm_2,stressm_3,stressm_4,stress12_1,           &
                   stress12_2,stress12_3,stress12_4,str1,str2,str3,              &
                   str4,str5,str6,str7,str8)
      !$OMP BARRIER
      call stepu(NA_len, rhow, &
                   1,nb,cdn_ocn,aiu,uocn,vocn,forcex,forcey, umassdti,fm,uarear, &
                   uvel_init,vvel_init,uvel,vvel,                                &
                   str1,str2,str3,str4,str5,str6,str7,str8, nw,sw,sse,skipucell)
      !$OMP BARRIER
      call halo_update(NA_len,1,navel,uvel,vvel, halo_parent)
      !$OMP BARRIER
    enddo
    call stress   (NA_len, tarear,                                               &
                   ee,ne,se,1,na,uvel,vvel,dxt,dyt,                              & 
                   tinyarea,dxhy,dyhx,cxp,cyp,cxm,cym,                           &
                   strength,stressp_1,stressp_2,stressp_3,stressp_4,             & 
                   stressm_1,stressm_2,stressm_3,stressm_4,stress12_1,           &
                   stress12_2,stress12_3,stress12_4,                             &
                   divu,rdg_conv,rdg_shear,shear,                                &
                   str1,str2,str3,str4,str5,str6,str7,str8)
    !$OMP BARRIER
    call stepu     (NA_len, rhow, &
                   1,nb,cdn_ocn,aiu,uocn,vocn,forcex,forcey, umassdti,fm,uarear, &
                   strocnx,strocny,strintx,strinty,                              &
                   uvel_init,vvel_init,uvel,vvel,                                &
                   str1,str2,str3,str4,str5,str6,str7,str8, nw,sw,sse,skipucell)
    !$OMP BARRIER
    call halo_update(NA_len,1,navel,uvel,vvel, halo_parent)
    !$OMP END PARALLEL
  end subroutine evp_kernel_v1
  !===============================================================================
  subroutine evp_kernel_v2
    use ice_constants, only : c0
    use ice_dyn_shared, only: ndte
    use bench_v2, only : stress, stepu, halo_update
    use dmi_omp, only : domp_init
    use icepack_intfc, only: icepack_query_parameters
    implicit none
    real(kind=dbl_kind) :: rhow
    integer (kind=int_kind) :: ierr, lun, i, nthreads
    integer (kind=int_kind) :: na,nb,navel
    !- Read constants...
    call icepack_query_parameters(rhow_out=rhow)
    na=NA_len
    nb=NA_len
    navel=NAVEL_len

    !- Initialize openmp ---------------------------------------------------------
    call domp_init(nthreads) ! ought to be called from main
  
    !- Initialize timers ---------------------------------------------------------
    str1=c0
    str2=c0
    str3=c0
    str4=c0
    str5=c0
    str6=c0
    str7=c0
    str8=c0
   
    if (ndte<2) STOP 'ndte must be 2 or higher for this kernel'
    !$OMP PARALLEL PRIVATE(i)
    do i = 1, ndte-1
      call stress (NA_len, &
                   ee,ne,se,1,na,uvel,vvel,dxt,dyt,                              & 
                   hte,htn,htem1,htnm1,                                          &
                   strength,stressp_1,stressp_2,stressp_3,stressp_4,             & 
                   stressm_1,stressm_2,stressm_3,stressm_4,stress12_1,           &
                   stress12_2,stress12_3,stress12_4,str1,str2,str3,              &
                   str4,str5,str6,str7,str8)
      !$OMP BARRIER
      call stepu(NA_len, rhow, &
                   1,nb,cdn_ocn,aiu,uocn,vocn,forcex,forcey, umassdti,fm,uarear, &
                   uvel_init,vvel_init,uvel,vvel,                                &
                   str1,str2,str3,str4,str5,str6,str7,str8, nw,sw,sse,skipucell)
      !$OMP BARRIER
      call halo_update(NA_len,1,navel,uvel,vvel, halo_parent)
      !$OMP BARRIER
    enddo
    call stress   (NA_len, tarear,                                               &
                   ee,ne,se,1,na,uvel,vvel,dxt,dyt,                              & 
                   hte,htn,htem1,htnm1,                                          &
                   strength,stressp_1,stressp_2,stressp_3,stressp_4,             & 
                   stressm_1,stressm_2,stressm_3,stressm_4,stress12_1,           &
                   stress12_2,stress12_3,stress12_4,                             &
                   divu,rdg_conv,rdg_shear,shear,                                &
                   str1,str2,str3,str4,str5,str6,str7,str8)
    !$OMP BARRIER
    call stepu     (NA_len, rhow, &
                   1,nb,cdn_ocn,aiu,uocn,vocn,forcex,forcey, umassdti,fm,uarear, &
                   strocnx,strocny,strintx,strinty,                              &
                   uvel_init,vvel_init,uvel,vvel,                                &
                   str1,str2,str3,str4,str5,str6,str7,str8, nw,sw,sse,skipucell)
    !$OMP BARRIER
    call halo_update(NA_len,1,navel,uvel,vvel, halo_parent)
    !$OMP END PARALLEL
  end subroutine evp_kernel_v2
  !===============================================================================
  subroutine calc_na(nx,ny,na,icetmask)
    ! Calculate number of active points (na)
    use ice_blocks, only: nghost
    implicit none
    integer(int_kind),intent(in) :: nx,ny
    integer(int_kind),intent(out) :: na
    integer (kind=int_kind),dimension (nx,ny), intent(in) :: icetmask
    integer(int_kind) :: i,j
    na = 0
! Note: The icellt mask includes north and east ghost cells. (ice_dyn_shared.F90)
    do j = 1+nghost, ny ! -nghost
    do i = 1+nghost, nx ! -nghost
        if (icetmask(i,j)==1) then 
          na=na+1
        endif
    enddo
    enddo
  end subroutine calc_na
  subroutine calc_2d_indices(nx,ny,na,icetmask,iceumask)
    use ice_blocks, only: nghost
    implicit none
    integer(int_kind),intent(in) :: nx,ny,na
    integer (kind=int_kind),dimension (nx,ny), intent(in) :: icetmask
    logical (kind=log_kind),dimension (nx,ny), intent(in) :: iceumask
    integer(int_kind) :: i,j,Nmaskt
    skipucell(:)=.false.
    indi=0
    indj=0
    Nmaskt=0
! Note: The icellt mask includes north and east ghost cells. (ice_dyn_shared.F90)
    do j = 1+nghost, ny ! -nghost
    do i = 1+nghost, nx ! -nghost
        if (icetmask(i,j)==1) then
          Nmaskt=Nmaskt+1
          indi(Nmaskt) = i
          indj(Nmaskt) = j
          ! Umask do NOT include north/east ghost cells ... skip these as well
          if (iceumask(i,j) .eqv. .false. ) skipucell(Nmaskt) = .true.
          if (i==nx) skipucell(Nmaskt) = .true.
          if (j==ny) skipucell(Nmaskt) = .true.
        endif
    enddo
    enddo
    if (Nmaskt.ne.na) then
      write(*,*)'Nmaskt,na: ',Nmaskt,na
      stop 'Problem Nmaskt != na'
    endif
    if (Nmaskt==0) then
      write(*,*)'WARNING: NO ICE'
    endif
  end subroutine calc_2d_indices
  subroutine calc_navel(nx_block,ny_block,na,navel)
    ! Calculate number of active points including needed halo points (navel)
    implicit none
    integer(int_kind),intent(in)  :: nx_block,ny_block,na
    integer(int_kind),intent(out) :: navel
    integer(int_kind) :: iw,i,j
    integer(int_kind),dimension(1:na) :: Iin,Iee,Ine,Ise,Inw,Isw,Isse
    integer(int_kind),dimension(1:7*na) :: util1,util2
    ! Additional indices used for finite differences (FD)
    do iw=1,na
      i=indi(iw)
      j=indj(iw)
      Iin(iw) = i   + (j-1)*nx_block  ! ( 0, 0) Target point
      Iee(iw) = i-1 + (j-1)*nx_block  ! (-1, 0)
      Ine(iw) = i-1 + (j-2)*nx_block  ! (-1,-1)
      Ise(iw) = i   + (j-2)*nx_block  ! ( 0,-1)
      Inw(iw) = i+1 + (j-1)*nx_block  ! (+1, 0)
      Isw(iw) = i+1 + (j-0)*nx_block  ! (+1,+1)
      Isse(iw)= i   + (j-0)*nx_block  ! ( 0,+1)
    enddo
    !-- Find number of points needed for finite difference calculations
    call union(Iin,  Iee,na,na,util1,i)
    call union(util1,Ine, i,na,util2,j)
    call union(util2,Ise, j,na,util1,i)
    call union(util1,Inw, i,na,util2,j)
    call union(util2,Isw, j,na,util1,i)
    call union(util1,Isse,i,na,util2,navel)
    !-- Check bounds
    do iw=1,navel
      if (util2(iw)>nx_block*ny_block .or. util2(iw)<1) then
        write(*,*)'nx_block,ny_block,nx_block*ny_block: ',nx_block,ny_block,nx_block*ny_block
        write(*,*)'na,navel,iw,util2(iw): ',na,navel,iw,util2(iw)
        stop 'Problem with boundary. Check halo zone values'
      endif
    enddo
  end subroutine calc_navel
  subroutine convert_2d_1d(nx,ny, na,navel,                               &
       I_HTE,I_HTN,I_dxhy,I_dyhx,I_cyp,I_cxp,I_cym,I_cxm,I_tinyarea,      &
       I_cdn_ocn,I_aiu,I_uocn,I_vocn,I_waterx,I_watery,I_forcex,I_forcey, &
       I_umassdti,I_fm,I_uarear,I_tarear,I_strintx,I_strinty,             &
       I_uvel_init,I_vvel_init,                                           &
       I_strength,I_uvel,I_vvel,I_dxt,I_dyt,                              &
        I_stressp_1 ,I_stressp_2, I_stressp_3, I_stressp_4,               &
        I_stressm_1 ,I_stressm_2, I_stressm_3, I_stressm_4,               &
       I_stress12_1,I_stress12_2,I_stress12_3,I_stress12_4                )
    implicit none
    integer(int_kind),intent(in) :: nx,ny,na,navel
    real (kind=dbl_kind), dimension(nx,ny), intent(in)    ::              &
       I_HTE,I_HTN,I_dxhy,I_dyhx,I_cyp,I_cxp,I_cym,I_cxm,I_tinyarea,      &
       I_cdn_ocn,I_aiu,I_uocn,I_vocn,I_waterx,I_watery,I_forcex,I_forcey, &
       I_umassdti,I_fm,I_uarear,I_tarear,I_strintx,I_strinty,             &
       I_uvel_init,I_vvel_init,                                           &
       I_strength,I_uvel,I_vvel,I_dxt,I_dyt,                              &
        I_stressp_1 ,I_stressp_2, I_stressp_3, I_stressp_4,               &
        I_stressm_1 ,I_stressm_2, I_stressm_3, I_stressm_4,               &
       I_stress12_1,I_stress12_2,I_stress12_3,I_stress12_4
    integer(int_kind) :: iw,i,j, nx_block
    integer(int_kind),dimension(1:na) :: Iin,Iee,Ine,Ise,Inw,Isw,Isse
    integer(int_kind),dimension(1:7*na) :: util1,util2
    integer(int_kind) :: nachk
    ! Additional indices used for finite differences (FD)
    nx_block=nx  ! Total block size in x-dir 
    do iw=1,na
      i=indi(iw)
      j=indj(iw)
      Iin(iw) = i   + (j-1)*nx_block  ! ( 0, 0) Target point
      Iee(iw) = i-1 + (j-1)*nx_block  ! (-1, 0)
      Ine(iw) = i-1 + (j-2)*nx_block  ! (-1,-1)
      Ise(iw) = i   + (j-2)*nx_block  ! ( 0,-1)
      Inw(iw) = i+1 + (j-1)*nx_block  ! (+1, 0)
      Isw(iw) = i+1 + (j-0)*nx_block  ! (+1,+1)
      Isse(iw)= i   + (j-0)*nx_block  ! ( 0,+1)
    enddo
    !-- Find number of points needed for finite difference calculations
    call union(Iin,  Iee,na,na,util1,i)
    call union(util1,Ine, i,na,util2,j)
    call union(util2,Ise, j,na,util1,i)
    call union(util1,Inw, i,na,util2,j)
    call union(util2,Isw, j,na,util1,i)
    call union(util1,Isse,i,na,util2,nachk)
    if (nachk .ne. navel) then
      write(*,*)'ERROR: navel badly chosen: na,navel,nachk = ',na,navel,nachk
      stop
    endif

    ! indij: vector with target points (sorted) ...
    do iw=1,na
      indij(iw)=Iin(iw)
    enddo
    ! indij: ... followed by extra points (sorted)
    call setdiff(util2,Iin,navel,na,util1,j)
    do iw=na+1,navel
      indij(iw)=util1(iw-na)
    enddo

    !-- Create indices for additional points needed for uvel,vvel:
    call findXinY(Iee ,indij,na,navel, ee)
    call findXinY(Ine ,indij,na,navel, ne)
    call findXinY(Ise ,indij,na,navel, se)
    call findXinY(Inw ,indij,na,navel, nw)
    call findXinY(Isw ,indij,na,navel, sw)
    call findXinY(Isse,indij,na,navel,sse)

    !-- write check
!if (1 == 2) then
!    write(*,*)'MHRI: INDICES start: convert_2d_1d'
!    write(*,*) 'Min/max ee', minval(ee), maxval(ee)
!    write(*,*) 'Min/max ne', minval(ne), maxval(ne)
!    write(*,*) 'Min/max se', minval(se), maxval(se)
!    write(*,*) 'Min/max nw', minval(nw), maxval(nw)
!    write(*,*) 'Min/max sw', minval(sw), maxval(sw)
!    write(*,*) 'Min/max sse',minval(sse),maxval(sse)
!    write(*,*)'MHRI: INDICES end: convert_2d_1d'
!endif

    ! Write 1D data from 2D: Here only extra FD part, the rest follows...
    !$OMP PARALLEL DO PRIVATE(iw,i,j)
    do iw=na+1,navel
      j=int((indij(iw)-1)/(nx_block))+1
      i=indij(iw)-(j-1)*nx_block
      uvel(iw)=      I_uvel(i,j)
      vvel(iw)=      I_vvel(i,j)
    enddo
    !$OMP END PARALLEL DO
    
    ! Write 1D data from 2D
    !$OMP PARALLEL DO PRIVATE(iw,i,j)
    do iw=1,na
      i=indi(iw)
      j=indj(iw)
      uvel(iw)=      I_uvel(i,j)
      vvel(iw)=      I_vvel(i,j)
      cdn_ocn(iw)=   I_cdn_ocn(i,j)  
      aiu(iw)=       I_aiu(i,j)   
      uocn(iw)=      I_uocn(i,j)   
      vocn(iw)=      I_vocn(i,j)   
      waterx(iw)=    I_waterx(i,j)   
      watery(iw)=    I_watery(i,j)   
      forcex(iw)=    I_forcex(i,j)   
      forcey(iw)=    I_forcey(i,j)   
      umassdti(iw)=  I_umassdti(i,j)   
      fm(iw)=        I_fm(i,j)   
      tarear(iw)=    I_tarear(i,j)   
      uarear(iw)=    I_uarear(i,j)   
      strintx(iw)=   I_strintx(i,j)   
      strinty(iw)=   I_strinty(i,j)   
      uvel_init(iw)= I_uvel_init(i,j)  
      vvel_init(iw)= I_vvel_init(i,j)  
      strength(iw)=  I_strength(i,j)   
      dxt(iw)=       I_dxt(i,j)   
      dyt(iw)=       I_dyt(i,j)   
      stressp_1(iw)= I_stressp_1(i,j)  
      stressp_2(iw)= I_stressp_2(i,j)  
      stressp_3(iw)= I_stressp_3(i,j)  
      stressp_4(iw)= I_stressp_4(i,j)  
      stressm_1(iw)= I_stressm_1(i,j)  
      stressm_2(iw)= I_stressm_2(i,j)  
      stressm_3(iw)= I_stressm_3(i,j)  
      stressm_4(iw)= I_stressm_4(i,j)  
      stress12_1(iw)=I_stress12_1(i,j) 
      stress12_2(iw)=I_stress12_2(i,j) 
      stress12_3(iw)=I_stress12_3(i,j) 
      stress12_4(iw)=I_stress12_4(i,j)
      !
      dxhy(iw)=      I_dxhy(i,j)   
      dyhx(iw)=      I_dyhx(i,j)   
      cyp(iw)=       I_cyp(i,j)   
      cxp(iw)=       I_cxp(i,j)   
      cym(iw)=       I_cym(i,j)   
      cxm(iw)=       I_cxm(i,j)   
      tinyarea(iw)=  I_tinyarea(i,j)   
      ! Grid space
      HTE(iw)   = I_HTE(i,j)
      HTN(iw)   = I_HTN(i,j)
      HTEm1(iw) = I_HTE(i-1,j)
      HTNm1(iw) = I_HTN(i,j-1)
    enddo
    !$OMP END PARALLEL DO
  end subroutine convert_2d_1d
  subroutine calc_halo_parent(nx,ny,na,navel, I_icetmask)
    implicit none
    integer(int_kind),intent(in) :: nx,ny,na,navel
    integer(kind=int_kind), dimension(nx,ny), intent(in) :: I_icetmask
    integer(int_kind) :: iw,i,j !,masku,maskt
    integer(int_kind),dimension(1:navel-na) :: Ihalo
    ! Indices for halo update
    ! TODO: ONLY for nghost==1
    ! TODO: ONLY for circular grids - NOT tripole grids
    Ihalo(:)=0
    !$OMP PARALLEL DO PRIVATE(iw,i,j)
    do iw=na+1,navel
      j=int((indij(iw)-1)/(nx))+1
      i=indij(iw)-(j-1)*nx
      if (i==nx .and. I_icetmask(   2,j)==1) Ihalo(iw-na)=     2+ (j-1)*nx
      if (i==1  .and. I_icetmask(nx-1,j)==1) Ihalo(iw-na)=(nx-1)+ (j-1)*nx
      if (j==ny .and. I_icetmask(i,   2)==1) Ihalo(iw-na)=     i+       nx
      if (j==1  .and. I_icetmask(i,ny-1)==1) Ihalo(iw-na)=     i+(ny-2)*nx
    enddo
    !$OMP END PARALLEL DO
    call findXinY_halo(Ihalo,indij(1:na),navel-na,na,halo_parent(na+1:navel))
!    halo_parent(1:na)=0 ! Done in numainit
  end subroutine calc_halo_parent
  !=======================================================================
  subroutine union(x,y,nx,ny,xy,nxy)
    ! Find union (xy) of two sorted integer vectors (x and y)
    !   ie. Combined values of the two vectors with no repetitions.
    !use ice_kinds_mod
    implicit none
    integer (int_kind) :: i,j,k
    integer (int_kind),intent(in)  :: nx,ny
    integer (int_kind),intent(in)  :: x(1:nx),y(1:ny)
    integer (int_kind),intent(out) :: xy(1:nx+ny)
    integer (int_kind),intent(out) :: nxy
    i=1
    j=1
    k=1
    do while (i<=nx .and. j<=ny)
      if (x(i)<y(j)) then
        xy(k)=x(i)
        i=i+1
      else if (x(i)>y(j)) then
        xy(k)=y(j)
        j=j+1
      else !if (x(i)==y(j)) then
        xy(k)=x(i)
        i=i+1
        j=j+1
      endif
      k=k+1
    enddo
    ! The rest
    do while (i<=nx)
      xy(k)=x(i)
      i=i+1
      k=k+1
    enddo
    do while (j<=ny)
      xy(k)=y(j)
      j=j+1
      k=k+1
    enddo
    nxy=k-1
  end subroutine union
  !=======================================================================
  subroutine setdiff(x,y,nx,ny,xy,nxy)
    ! Find element (xy) of two sorted integer vectors (x and y)
    !   that are in x, but not in y ... or in y, but not in x
    !use ice_kinds_mod
    implicit none
    integer (int_kind) :: i,j,k
    integer (int_kind),intent(in)  :: nx,ny
    integer (int_kind),intent(in)  :: x(1:nx),y(1:ny)
    integer (int_kind),intent(out) :: xy(1:nx+ny)
    integer (int_kind),intent(out) :: nxy
    i=1
    j=1
    k=1
    do while (i<=nx .and. j<=ny)
      if (x(i)<y(j)) then
        xy(k)=x(i)
        i=i+1
        k=k+1
      else if (x(i)>y(j)) then
        xy(k)=y(j)
        j=j+1
        k=k+1
      else !if (x(i)==y(j)) then
        i=i+1
        j=j+1
      endif
    enddo
    ! The rest
    do while (i<=nx)
      xy(k)=x(i)
      i=i+1
      k=k+1
    enddo
    do while (j<=ny)
      xy(k)=y(j)
      j=j+1
      k=k+1
    enddo
    nxy=k-1
  end subroutine setdiff
  !=======================================================================
  subroutine findXinY(x,y,nx,ny,indx)
    ! Find indx vector so that x(1:na)=y(indx(1:na))
    !
    !  Conditions:
    !   * EVERY item in x is found in y.
    !   * x(1:nx) is a sorted integer vector.
    !   * y(1:ny) consists of two sorted integer vectors:
    !       [y(1:nx) ; y(nx+1:ny)]
    !   * ny>=nx
    !  Return: indx(1:na)
    !
    !use ice_kinds_mod
    implicit none
    integer (int_kind),intent(in)  :: nx,ny
    integer (int_kind),intent(in)  :: x(1:nx),y(1:ny)
    integer (int_kind),intent(out) :: indx(1:nx)
    integer (int_kind) :: i,j1,j2
    i=1
    j1=1
    j2=nx+1
    do while (i<=nx)
      if (x(i)==y(j1)) then
        indx(i)=j1
        i=i+1
        j1=j1+1
      else if (x(i)==y(j2)) then
        indx(i)=j2
        i=i+1
        j2=j2+1
      else if (x(i)>y(j1) ) then !.and. j1<nx) then
        j1=j1+1
      else if (x(i)>y(j2) ) then !.and. j2<ny) then
        j2=j2+1
      else
        write(*,*)'nx,ny: ',nx,ny
        write(*,*)'i,j1,j2: ',i,j1,j2
        write(*,*)'x(i),y(j1),y(j2): ',x(i),y(j1),y(j2)
        stop 'ERROR in findXinY'
      endif
    end do
  end subroutine findXinY
  subroutine findXinY_halo(x,y,nx,ny,indx)
    ! Find indx vector so that x(1:na)=y(indx(1:na))
    !
    !  Conditions:
    !   * EVERY item in x is found in y.
    !       Except for x==0, where indx=0 is returned
    !   * x(1:nx) is an integer vector. Not sorted.
    !   * y(1:ny) is a sorted integer vector
    !   * ny>=nx
    !  Return: indx(1:na)
    !
    !use ice_kinds_mod
    implicit none
    integer (int_kind),intent(in)  :: nx,ny
    integer (int_kind),intent(in)  :: x(1:nx),y(1:ny)
    integer (int_kind),intent(out) :: indx(1:nx)
    integer (int_kind) :: i,j1,nloop
    nloop=1
    i=1
    j1=int((ny+1)/2) ! initial guess in the middle
    do while (i<=nx)
      if (x(i)==0) then
        indx(i)=0
        i=i+1
        nloop=1
      else if (x(i)==y(j1)) then
        indx(i)=j1
        i=i+1
        j1=j1+1
        if (j1>ny) j1=int((ny+1)/2) ! initial guess in the middle
        nloop=1
      else if (x(i)<y(j1) ) then
        j1=1
      else if (x(i)>y(j1) ) then
        j1=j1+1
        if (j1>ny) then
          j1=1
          nloop=nloop+1
          if (nloop>2) then
            ! Stop for inf. loop. This check should not be necessary for halo
            write(*,*)'nx,ny: ',nx,ny
            write(*,*)'i,j1: ',i,j1
            write(*,*)'x(i),y(j1): ',x(i),y(j1)
            stop 'ERROR in findXinY_halo: too many loops'
          endif
        endif
      endif
    end do
  end subroutine findXinY_halo

  !=======================================================================
  subroutine numainit(l,u,uu)
    !- modules -----------------------------------------------------------------
    use dmi_omp, only  : domp_get_domain
    use ice_constants, only: c0
!    use ice_kinds_mod
!    use vars
    implicit none
    integer(kind=int_kind),intent(in) :: l,u,uu
    integer(kind=int_kind) :: lo,up
    call domp_get_domain(l,u,lo,up)
    ee(lo:up)=0
    ne(lo:up)=0
    se(lo:up)=0
    sse(lo:up)=0
    nw(lo:up)=0
    sw(lo:up)=0
    halo_parent(lo:up)=0
    strength(lo:up)=c0
    uvel(lo:up)=c0
    vvel(lo:up)=c0
    dxt(lo:up)=c0
    dyt(lo:up)=c0
    HTE(lo:up)=c0
    HTN(lo:up)=c0
    HTEm1(lo:up)=c0
    HTNm1(lo:up)=c0
    dxhy(lo:up)=c0
    dyhx(lo:up)=c0
    cyp(lo:up)=c0
    cxp(lo:up)=c0
    cym(lo:up)=c0
    cxm(lo:up)=c0
    tinyarea(lo:up)=c0
    stressp_1(lo:up)=c0
    stressp_2(lo:up)=c0
    stressp_3(lo:up)=c0
    stressp_4(lo:up)=c0
    stressm_1(lo:up)=c0
    stressm_2(lo:up)=c0
    stressm_3(lo:up)=c0
    stressm_4(lo:up)=c0
    stress12_1(lo:up)=c0
    stress12_2(lo:up)=c0
    stress12_3(lo:up)=c0
    stress12_4(lo:up)=c0
    tarear(lo:up)=c0
    divu(lo:up)=c0
    rdg_conv(lo:up)=c0
    rdg_shear(lo:up)=c0
    shear(lo:up)=c0
    str1(lo:up)=c0
    str2(lo:up)=c0
    str3(lo:up)=c0
    str4(lo:up)=c0
    str5(lo:up)=c0
    str6(lo:up)=c0
    str7(lo:up)=c0
    str8(lo:up)=c0
    call domp_get_domain(u+1,uu,lo,up)
    halo_parent(lo:up)=0
    uvel(lo:up)=c0
    vvel(lo:up)=c0
    str1(lo:up)=c0
    str2(lo:up)=c0
    str3(lo:up)=c0
    str4(lo:up)=c0
    str5(lo:up)=c0
    str6(lo:up)=c0
    str7(lo:up)=c0
    str8(lo:up)=c0
  end subroutine numainit
  
  !=======================================================================
end module evp_kernel1d

