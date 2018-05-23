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
module bench
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
       strength, uvel, vvel, dxt, dyt, tinyarea,dxhy,dyhx,cxp,cyp,cxm,cym
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
    do iw = 1,NA_len  ! FIXME hardcoding
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
!MHRI       str1,str2,str3,str4,str5,str6,str7,str8, prs_sig  
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
!MHRI       prs_sig(iw) = c0ne*Deltane ! northeast
  
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
                    lb,ub,Cw,aiu,uocn,vocn,forcex,forcey,                        &
                    umassdti, fm,uarear,uvel_init,                               &
                    vvel_init,uvel,vvel,str1,str2,str3,str4,str5,str6,str7,str8, &
                    nw,sw,se,skipme)
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
!    do iw = 1,144481  ! FIXME hardcoding
    do iw = 1,NA_len  ! FIXME hardcoding
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
                    lb,ub,Cw,aiu,uocn,vocn,forcex,forcey,umassdti, &
                    fm,uarear,strocnx,strocny,strintx,strinty,uvel_init,         &
                    vvel_init,uvel,vvel,str1,str2,str3,str4,str5,str6,str7,str8, &
                    nw,sw,se,skipme)
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
!    do iw = 1,144481  ! FIXME hardcoding
    do iw = 1,NA_len  ! FIXME hardcoding
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
end module bench
  
!===============================================================================
  
!-- One dimension representation of EVP 2D arrays used for EVP kernels
module evp_kernel1d
  use ice_kinds_mod
  !-- BEGIN: specific for the KERNEL
  use ice_dyn_shared, only: revp, ecci, denom1, arlx1i, brlx
  !-- END: specific for the KERNEL
  implicit none
  public :: evp_copyin, evp_copyout, evp_kernel
  private
  save
  integer(kind=int_kind) ::                                        &
    NA_len, NAVEL_len
  logical(kind=log_kind), dimension(:), allocatable ::             &
    skiptcell,skipucell
  integer(kind=int_kind), dimension(:), allocatable ::             &
    ee,ne,se,nw,sw,sse,indi,indj,indij                  
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
!MHRI    str1, str2, str3, str4, str5, str6, str7, str8, prs_sig
  real (kind=dbl_kind), dimension(:), allocatable ::               &
    HTE,HTN,                                                       &
!    HTE_m10,HTN_m01
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
      skiptcell(1:na),skipucell(1:na),                             &
      ! Grid distances: HTE,HTN + "-1 neighbours" 
      HTE(1:na),HTN(1:na),                                         &
      HTEm1(1:na),HTNm1(1:na),                                     &
      ! T cells
      strength(1:na),dxt(1:na),dyt(1:na),dxhy(1:na),dyhx(1:na),           &
      cyp(1:na),cxp(1:na),cym(1:na),cxm(1:na),tinyarea(1:na),tarear(1:na),&
      stressp_1(1:na), stressp_2(1:na), stressp_3(1:na), stressp_4(1:na), &
      stressm_1(1:na), stressm_2(1:na), stressm_3(1:na), stressm_4(1:na), &
      stress12_1(1:na),stress12_2(1:na),stress12_3(1:na),stress12_4(1:na),&
      str1(1:na), str2(1:na),str3(1:na),str4(1:na),                       &
      str5(1:na), str6(1:na),str7(1:na),str8(1:na),                       &
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
      uvel(1:navel),vvel(1:navel), indij(1:navel),                       &
!MHRI      prs_sig(1:navel),                                                  &
      stat=ierr)
    if (ierr/=0) stop 'Error allocating 1D navel'
  end subroutine alloc1d_navel
  subroutine dealloc1d
    implicit none
    integer(kind=int_kind) :: ierr
    deallocate(                                              &
      ! U+T cells
      ! Helper index for neighbours
      indj,indi,                                       &
      ee,ne,se,                                  &
      nw,sw,sse,                                 &
      skiptcell,skipucell,                             &
      ! Grid distances: HTE,HTN + "-1 neighbours" 
      HTE,HTN,                                         &
      HTEm1,HTNm1,                                     &
      ! T cells
      strength,dxt,dyt,dxhy,dyhx,           &
      cyp,cxp,cym,cxm,tinyarea,tarear,&
      stressp_1, stressp_2, stressp_3, stressp_4, &
      stressm_1, stressm_2, stressm_3, stressm_4, &
      stress12_1,stress12_2,stress12_3,stress12_4,&
      str1, str2,str3,str4,                       &
      str5, str6,str7,str8,                       &
      divu,rdg_conv,rdg_shear,shear,              &
      ! U cells
      cdn_ocn,aiu,uocn,vocn,                      &
      waterx,watery,forcex,forcey,                &
      umassdti,fm,uarear,                               &
      strocnx,strocny,strintx,strinty,            &
      uvel_init,vvel_init,                                    &
      ! NAVEL 
      uvel,vvel, indij,           &
!MHRI      prs_sig,                                      &
      stat=ierr)
    if (ierr/=0) stop 'Error de-allocating 1D'
  end subroutine dealloc1d
  subroutine write1d(na,navel)
    integer(kind=int_kind),intent(in) :: na,navel
    integer(kind=int_kind) :: lun, ios
    integer(kind=int_kind) :: nb
    nb=na
    lun=77
    open(lun,file='EVP_conv1D.bin', form='unformatted', access='stream',   &
             action='write', iostat=ios)
    if (ios/=0) stop 'Problem opening file EVP_conv1D.bin'
    write(lun,iostat=ios)                                                  &
      ! T cells
      strength(1:na),dxt(1:na),dyt(1:na),dxhy(1:na),dyhx(1:na),            &
      cyp(1:na),cym(1:na),cxm(1:na),tinyarea(1:na),stressp_1(1:na),        &
      stressp_3(1:na),stressp_4(1:na),stressm_1(1:na),stressm_2(1:na),     &
      stressm_3(1:na),stressm_4(1:na),stress12_1(1:na),                    &
      cxp(1:na),stressp_2(1:na),                                           &
      stress12_2(1:na),stress12_3(1:na),stress12_4(1:na),tarear(1:na),     &
      ! U cells
      cdn_ocn(1:nb),aiu(1:nb),uocn(1:nb),vocn(1:nb),                       &
      waterx(1:nb),watery(1:nb),forcex(1:nb),forcey(1:nb),                 &
      umassdti(1:nb),fm(1:nb),                                             &
      uarear(1:nb),strintx(1:nb),strinty(1:nb),                            &
      uvel_init(1:nb),vvel_init(1:nb),                                     &
      ! common
      uvel(1:navel),vvel(1:navel),                                         &
      ! Helper index for neighbours
!      indi(1:navel),indj(1:navel),                   &
      ee(1:na),ne(1:na),se(1:na),                                          &
      nw(1:nb),sw(1:nb),sse(1:nb), skiptcell(1:na),skipucell(1:na)
    if (ios/=0) stop 'Problem writing file EVP_conv1D.bin'
    close(lun)
  end subroutine write1d
!===============================================================================
!===============================================================================
  subroutine evp_copyin(nx,ny,na,                                                             &
       I_icetmask,I_iceumask, I_HTE,I_HTN,                                                    &
       I_cdn_ocn,I_aiu,I_uocn,I_vocn,I_waterx,I_watery,I_forcex,I_forcey,                     &
       I_umassdti,I_fm,I_uarear,I_tarear,I_strintx,I_strinty,I_uvel_init,I_vvel_init,         &
       I_strength,I_uvel,I_vvel,I_dxt,I_dyt,I_dxhy,I_dyhx,I_cyp,I_cxp,I_cym,I_cxm,I_tinyarea, &
        I_stressp_1 ,I_stressp_2, I_stressp_3, I_stressp_4,                                   &
        I_stressm_1 ,I_stressm_2, I_stressm_3, I_stressm_4,                                   &
       I_stress12_1,I_stress12_2,I_stress12_3,I_stress12_4                                    )
    implicit none
    integer(int_kind), intent(in) :: nx, ny, na
    integer (kind=int_kind),dimension (nx,ny), intent(in) :: I_icetmask
    logical (kind=log_kind),dimension (nx,ny), intent(in) :: I_iceumask
    real (kind=dbl_kind), dimension(nx,ny,1), intent(in)  ::                                  &
       I_HTE,I_HTN,                                                                           &
       I_cdn_ocn,I_aiu,I_uocn,I_vocn,I_waterx,I_watery,I_forcex,I_forcey,                     &
       I_umassdti,I_fm,I_uarear,I_tarear,I_strintx,I_strinty,I_uvel_init,I_vvel_init,         &
       I_strength,I_uvel,I_vvel,I_dxt,I_dyt,I_dxhy,I_dyhx,I_cyp,I_cxp,I_cym,I_cxm,I_tinyarea, &
       I_stressp_1, I_stressp_2, I_stressp_3, I_stressp_4,                                    &
       I_stressm_1, I_stressm_2, I_stressm_3, I_stressm_4,                                    &
       I_stress12_1,I_stress12_2,I_stress12_3,I_stress12_4
    integer(int_kind) :: navel
!    na=icellt(1)    ! MHRI/DMI: assume iblk==1
!    nx=nx_block    ! MHRI/DMI: assume iblk==1
!    ny=ny_block    ! MHRI/DMI: assume iblk==1
    call alloc1d(na)
    call calc_2d_indices(nx,ny,na, I_icetmask,I_iceumask)
    call calc_navel(nx,ny,na,navel)
    call alloc1d_navel(navel)
!$OMP PARALLEL DEFAULT(shared) 
    call numainit(1,na,navel)
!$OMP END PARALLEL
    ! Remap 2d to 1d and fill in
    call convert_2d_1d(nx,ny,na,navel, I_HTE,I_HTN,                                           &
       I_cdn_ocn,I_aiu,I_uocn,I_vocn,I_waterx,I_watery,I_forcex,I_forcey,                     &
       I_umassdti,I_fm,I_uarear,I_tarear,I_strintx,I_strinty,I_uvel_init,I_vvel_init,         &
       I_strength,I_uvel,I_vvel,I_dxt,I_dyt,I_dxhy,I_dyhx,I_cyp,I_cxp,I_cym,I_cxm,I_tinyarea, &
        I_stressp_1, I_stressp_2, I_stressp_3, I_stressp_4,                                   &
        I_stressm_1, I_stressm_2, I_stressm_3, I_stressm_4,                                   &
       I_stress12_1,I_stress12_2,I_stress12_3,I_stress12_4                                    )
    NA_len=na
    NAVEL_len=navel
    !-- write check
!if (1 == 2) then
!  write(*,*)'MHRI: INDICES start: evp-copyin'
!  write(*,*) 'Min/max ee', minval(ee(1:na)), maxval(ee(1:na))
!  write(*,*) 'Min/max ne', minval(ne(1:na)), maxval(ne(1:na))
!  write(*,*) 'Min/max se', minval(se(1:na)), maxval(se(1:na))
!  write(*,*) 'Min/max nw', minval(nw(1:na)), maxval(nw(1:na))
!  write(*,*) 'Min/max sw', minval(sw(1:na)), maxval(sw(1:na))
!  write(*,*) 'Min/max sse', minval(sse(1:na)), maxval(sse(1:na))
!  write(*,*)'MHRI: INDICES end: evp-copyin'
!endif
!MHRI chk
!    call write1d(na,navel)
!write(*,*)'HTE(1:10)',HTE(1:10)
!MHRI chk
  end subroutine evp_copyin
  subroutine evp_copyout(nx,ny,                                         &
               I_uvel,I_vvel, I_strintx,I_strinty, I_strocnx,I_strocny, &
               I_stressp_1, I_stressp_2, I_stressp_3, I_stressp_4,      &
               I_stressm_1, I_stressm_2, I_stressm_3, I_stressm_4,      &
               I_stress12_1,I_stress12_2,I_stress12_3,I_stress12_4,     &
               I_divu,I_rdg_conv,I_rdg_shear,I_shear                    )
    use ice_constants, only : c0
    implicit none
    integer(int_kind), intent(in) :: nx,ny
    real(dbl_kind), dimension(nx,ny,1), intent(out) :: &
       I_uvel,I_vvel, I_strintx,I_strinty, I_strocnx,I_strocny, &
       I_stressp_1, I_stressp_2, I_stressp_3, I_stressp_4,      &
       I_stressm_1, I_stressm_2, I_stressm_3, I_stressm_4,      &
       I_stress12_1,I_stress12_2,I_stress12_3,I_stress12_4,     &
       I_divu,I_rdg_conv, I_rdg_shear,I_shear
    integer(int_kind) :: i,j,iw
    ! Remap 1d to 2d and fill in
!    do iw=1,NAVEL_len
!      j=int((indij(iw)-1)/(nx))+1
!      i=indij(iw)-(j-1)*nx
    do iw=1,NA_len
      i=indi(iw)
      j=indj(iw)
      I_uvel(i,j,1)       = uvel(iw)
      I_vvel(i,j,1)       = vvel(iw)
      I_strintx(i,j,1)    = strintx(iw)
      I_strinty(i,j,1)    = strinty(iw)
      I_strocnx(i,j,1)    = strocnx(iw)
      I_strocny(i,j,1)    = strocny(iw)
      I_stressp_1(i,j,1)  = stressp_1(iw)
      I_stressp_2(i,j,1)  = stressp_2(iw)
      I_stressp_3(i,j,1)  = stressp_3(iw)
      I_stressp_4(i,j,1)  = stressp_4(iw)
      I_stressm_1(i,j,1)  = stressm_1(iw)
      I_stressm_2(i,j,1)  = stressm_2(iw)
      I_stressm_3(i,j,1)  = stressm_3(iw)
      I_stressm_4(i,j,1)  = stressm_4(iw)
      I_stress12_1(i,j,1) = stress12_1(iw)
      I_stress12_2(i,j,1) = stress12_2(iw)
      I_stress12_3(i,j,1) = stress12_3(iw)
      I_stress12_4(i,j,1) = stress12_4(iw)
      I_divu(i,j,1)       = divu(iw)
      I_rdg_conv(i,j,1)   = rdg_conv(iw)
      I_rdg_shear(i,j,1)  = rdg_shear(iw)
      I_shear(i,j,1)      = shear(iw)
    enddo
    call dealloc1d()
!   do j=1,ny
!     do i=1,nx
!       rdg_conv(i,j)  = -min(divu(i,j),c0)  ! From stress routine
!     enddo
!   enddo
  end subroutine evp_copyout
  !===============================================================================
  subroutine evp_kernel
    use ice_constants, only : c0
    use ice_dyn_shared, only: ndte
    use bench, only : stress, stepu
    use dmi_omp, only : domp_init
    use icepack_intfc, only: icepack_query_parameters
    implicit none
    real(kind=dbl_kind) :: rhow
    integer (kind=int_kind) :: ierr, lun, i, nthreads
    integer (kind=int_kind) :: na,nb,navel
    !- Read constants...
    call icepack_query_parameters(rhow_out=rhow)
    na=NA_len  ! TODO: full scale
    nb=NA_len  ! TODO: full scale
    navel=NAVEL_len  ! TODO: full scale
!write(*,*)'ndte,rhow',ndte,rhow
!write(*,*)'na,nb,navel',na,nb,navel

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
   
    !-- write check
!if (1 == 1) then
!  write(*,*)'MHRI: INDICES start: evp-kernel'
!  write(*,*) 'Min/max ee', minval(ee(1:na)), maxval(ee(1:na))
!  write(*,*) 'Min/max ne', minval(ne(1:na)), maxval(ne(1:na))
!  write(*,*) 'Min/max se', minval(se(1:na)), maxval(se(1:na))
!  write(*,*) 'Min/max nw', minval(nw(1:nb)), maxval(nw(1:nb))
!  write(*,*) 'Min/max sw', minval(sw(1:nb)), maxval(sw(1:nb))
!  write(*,*) 'Min/max sse', minval(sse(1:nb)), maxval(sse(1:nb))
!  write(*,*)'MHRI: INDICES end: evp-kernel'
!endif

!write(*,*)'HTE(1:10)',HTE(1:10)
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
      call stepu  (NA_len,rhow, &
                   1,nb,cdn_ocn,aiu,uocn,vocn,forcex,forcey,                     &
                   umassdti,fm,uarear,uvel_init,vvel_init,uvel,vvel,             &
                   str1,str2,str3,str4,str5,str6,str7,str8,nw,sw,sse,skipucell)
  !$OMP BARRIER
    enddo
    call stress   (NA_len, tarear,                                               &
                   ee,ne,se,1,na,uvel,vvel,dxt,dyt,                              & 
                   tinyarea,dxhy,dyhx,cxp,cyp,cxm,cym,                           &
                   strength,stressp_1,stressp_2,stressp_3,stressp_4,             & 
                   stressm_1,stressm_2,stressm_3,stressm_4,stress12_1,           &
                   stress12_2,stress12_3,stress12_4,                             &
!MHRI              prs_sig,                                                      &
                   divu,rdg_conv,rdg_shear,shear,                                &
                   str1,str2,str3,str4,str5,str6,str7,str8)
  !$OMP BARRIER
    call stepu    (NA_len,rhow, &
                   1,nb,cdn_ocn,aiu,uocn,vocn,forcex,                            &
                   forcey,umassdti,fm,uarear,strocnx,strocny,strintx,strinty,    & 
                   uvel_init,vvel_init,uvel,vvel,str1,str2,str3,                 &
                   str4,str5,str6,str7,str8,nw,sw,sse,skipucell)
  !$OMP END PARALLEL
  end subroutine evp_kernel
  !===============================================================================
  subroutine calc_2d_indices(nx,ny,na,icetmask,iceumask)
    implicit none
    integer(int_kind),intent(in) :: nx,ny,na
    integer (kind=int_kind),dimension (nx,ny), intent(in) :: icetmask
    logical (kind=log_kind),dimension (nx,ny), intent(in) :: iceumask
    integer(int_kind) :: i1d,i,j,masku,maskt
    skiptcell(:)=.false.
    skipucell(:)=.false.
    masku=0
    maskt=0
    indi=0
    indj=0
    i1d = 0
    do j = 1, ny ! FIXME check bounds
      do i = 1, nx ! FIXME check bounds
         if ((icetmask(i,j)==1)) then 
           if ((iceumask(i,j)) .eqv..true.) then
             i1d=i1d+1
             indi(i1d) = i
             indj(i1d) = j
           else
             i1d=i1d+1
             indi(i1d) = i
             indj(i1d) = j
             skipucell(i1d) = .true.
           endif
         endif
         if ((iceumask(i,j).eqv..true.)) then
           masku=masku+1
         endif
         if ((icetmask(i,j)==1)) then
           maskt=maskt+1
         endif
      enddo
    enddo
    if (i1d.ne.na) then
      write(*,*)'i1d,na: ',i1d,na
      stop 'Problem i1d != na'
    endif
    if (masku.ge.maskt) then
      if (maskt==0) then
        write(*,*)'WARNING: NO ICE: masku, maskt: ',masku,maskt
      else
        write(*,*)'masku, maskt: ',masku,maskt
        stop 'Problem masku>= maskt'   ! FIXME: JWP/MHRI: Are you sure this is not possible??
      endif
    endif
  end subroutine calc_2d_indices
  subroutine calc_navel(nx_block,ny_block,na,navel)
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
       I_HTE,I_HTN,                                                       &
       I_cdn_ocn,I_aiu,I_uocn,I_vocn,I_waterx,I_watery,I_forcex,I_forcey, &
       I_umassdti,I_fm,I_uarear,I_tarear,I_strintx,I_strinty,             &
       I_uvel_init,I_vvel_init,                                           &
       I_strength,I_uvel,I_vvel,I_dxt,I_dyt,I_dxhy,I_dyhx,                &
       I_cyp,I_cxp,I_cym,I_cxm,I_tinyarea,                                &
        I_stressp_1 ,I_stressp_2, I_stressp_3, I_stressp_4,               &
        I_stressm_1 ,I_stressm_2, I_stressm_3, I_stressm_4,               &
       I_stress12_1,I_stress12_2,I_stress12_3,I_stress12_4                )
    implicit none
    integer(int_kind),intent(in) :: nx,ny,na,navel
    real (kind=dbl_kind), dimension(nx,ny,1), intent(in)    ::            &
       I_HTE,I_HTN,                                                       &
       I_cdn_ocn,I_aiu,I_uocn,I_vocn,I_waterx,I_watery,I_forcex,I_forcey, &
       I_umassdti,I_fm,I_uarear,I_tarear,I_strintx,I_strinty,             &
       I_uvel_init,I_vvel_init,                                           &
       I_strength,I_uvel,I_vvel,I_dxt,I_dyt,I_dxhy,I_dyhx,                &
       I_cyp,I_cxp,I_cym,I_cxm,I_tinyarea,                                &
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
!if (1 == 1) then
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
    do iw=na+1,navel
      j=int((indij(iw)-1)/(nx_block))+1
      i=indij(iw)-(j-1)*nx_block
      uvel(iw)=      I_uvel(i,j,1)
      vvel(iw)=      I_vvel(i,j,1)
    enddo
    
    ! Write 1D data from 2D
    do iw=1,na
      i=indi(iw)
      j=indj(iw)
      uvel(iw)=      I_uvel(i,j,1)
      vvel(iw)=      I_vvel(i,j,1)
      cdn_ocn(iw)=   I_cdn_ocn(i,j,1)  
      aiu(iw)=       I_aiu(i,j,1)   
      uocn(iw)=      I_uocn(i,j,1)   
      vocn(iw)=      I_vocn(i,j,1)   
      waterx(iw)=    I_waterx(i,j,1)   
      watery(iw)=    I_watery(i,j,1)   
      forcex(iw)=    I_forcex(i,j,1)   
      forcey(iw)=    I_forcey(i,j,1)   
      umassdti(iw)=  I_umassdti(i,j,1)   
      fm(iw)=        I_fm(i,j,1)   
      tarear(iw)=    I_tarear(i,j,1)   
      uarear(iw)=    I_uarear(i,j,1)   
      strintx(iw)=   I_strintx(i,j,1)   
      strinty(iw)=   I_strinty(i,j,1)   
      uvel_init(iw)= I_uvel_init(i,j,1)  
      vvel_init(iw)= I_vvel_init(i,j,1)  
      strength(iw)=  I_strength(i,j,1)   
      dxt(iw)=       I_dxt(i,j,1)   
      dyt(iw)=       I_dyt(i,j,1)   
      dxhy(iw)=      I_dxhy(i,j,1)   
      dyhx(iw)=      I_dyhx(i,j,1)   
      cyp(iw)=       I_cyp(i,j,1)   
      cxp(iw)=       I_cxp(i,j,1)   
      cym(iw)=       I_cym(i,j,1)   
      cxm(iw)=       I_cxm(i,j,1)   
      tinyarea(iw)=  I_tinyarea(i,j,1)   
      stressp_1(iw)= I_stressp_1(i,j,1)  
      stressp_2(iw)= I_stressp_2(i,j,1)  
      stressp_3(iw)= I_stressp_3(i,j,1)  
      stressp_4(iw)= I_stressp_4(i,j,1)  
      stressm_1(iw)= I_stressm_1(i,j,1)  
      stressm_2(iw)= I_stressm_2(i,j,1)  
      stressm_3(iw)= I_stressm_3(i,j,1)  
      stressm_4(iw)= I_stressm_4(i,j,1)  
      stress12_1(iw)=I_stress12_1(i,j,1) 
      stress12_2(iw)=I_stress12_2(i,j,1) 
      stress12_3(iw)=I_stress12_3(i,j,1) 
      stress12_4(iw)=I_stress12_4(i,j,1)
      ! Grid space
      HTE(iw)     = I_HTE(i  ,j  ,1)
      HTN(iw)     = I_HTN(i  ,j  ,1)
      HTEm1(iw) = I_HTE(i-1,j  ,1)
      HTNm1(iw) = I_HTN(i  ,j-1,1)
    enddo
  end subroutine convert_2d_1d
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
!write(*,*)'i,j,k,x(i),y(j):',i,j,k,x(i),y(j)
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
!write(*,*)'I,k,x(i),y(j):',i,j-1,k,x(i),y(j-1)
      xy(k)=x(i)
      i=i+1
      k=k+1
    enddo
    do while (j<=ny)
!write(*,*)'i,J,k,x(i),y(j):',i-1,j,k,x(i-1),y(j)
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

  !=======================================================================
  subroutine numainit(l,u,uu)
    !- modules -----------------------------------------------------------------
    use dmi_omp, only  : domp_get_domain
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
    strength(lo:up)= 0.0
    uvel(lo:up)=0.0
    vvel(lo:up)=0.0
    dxt(lo:up)=0.0
    dyt(lo:up)=0.0
    dxhy(lo:up)=0.0
    dyhx(lo:up)=0.0
    cyp(lo:up)=0.0
    cxp(lo:up)=0.0
    cym(lo:up)=0.0
    cxm(lo:up)=0.0
    tinyarea(lo:up)=0.0
    stressp_1(lo:up)= 0.0
    stressp_2(lo:up)=0.0
    stressp_3(lo:up)=0.0
    stressp_4(lo:up)=0.0
    stressm_1(lo:up)=0.0
    stressm_2(lo:up)=0.0
    stressm_3(lo:up)=0.0
    stressm_4(lo:up)=0.0
    stress12_1(lo:up)=0.0
    stress12_2(lo:up)=0.0
    stress12_3(lo:up)=0.0
    stress12_4(lo:up)=0.0
!MHRI    prs_sig(lo:up)=0.0
    tarear(lo:up)=0.0
    divu(lo:up)=0.0
    rdg_conv(lo:up)=0.0
    rdg_shear(lo:up)=0.0
    shear(lo:up)=0.0
    str1(lo:up)=0.0
    str2(lo:up)=0.0
    str3(lo:up)=0.0
    str4(lo:up)=0.0
    str5(lo:up)=0.0
    str6(lo:up)=0.0
    str7(lo:up)=0.0
    str8(lo:up)=0.0
    call domp_get_domain(u+1,uu,lo,up)
    uvel(lo:up)=0.0
    vvel(lo:up)=0.0
  end subroutine numainit
  
  !=======================================================================
end module evp_kernel1d

