! ice_dyn_evp_1d
!
! Contained 3 Fortran modules,
!   * dmi_omp
!   * bench_v2
!   * ice_dyn_evp_1d
! These were merged into one module, ice_dyn_evp_1d to support some
!   coupled build systems.
!
! Modules used for:
!   * convert 2D arrays into 1D vectors
!   * Do stress/stepu/halo_update interations
!   * convert 1D vectors into 2D matrices
!
!   Call from ice_dyn_evp.F90:
!     call ice_dyn_evp_1d_copyin(...)
!     call ice_dyn_evp_1d_kernel()
!     call ice_dyn_evp_1d_copyout(...)
!
! * REAL4 internal version:
!     mv evp_kernel1d.F90 evp_kernel1d_r8.F90
!     cat evp_kernel1d_r8.F90 | sed s/DBL_KIND/REAL_KIND/g > evp_kernel1d.F90
!
! * !v1 : a) "dxhy,dyhx,cyp,cxp,cym,cxm,tinyarea" input variables is replaced by
!            "HTE,HTN"->"HTE,HTN,HTEm1,HTNm1" and variables are calculated in-line
!         b) "waterx,watery" is calculated using existing input "uocn,vocn"
!
! Jacob Weismann Poulsen (JWP), Mads Hvid Ribergaard (MHRI), DMI
!===============================================================================

!===============================================================================
  
!-- One dimension representation of EVP 2D arrays used for EVP kernels
module ice_dyn_evp_1d

  use ice_kinds_mod
  use ice_fileunits, only: nu_diag
  use ice_exit, only: abort_ice
  use icepack_intfc, only: icepack_query_parameters
  use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted

  implicit none
  private
  public :: ice_dyn_evp_1d_copyin, ice_dyn_evp_1d_copyout, ice_dyn_evp_1d_kernel

  interface ice_dyn_evp_1d_copyin
!    module procedure evp_copyin_v1
    module procedure evp_copyin_v2
  end interface 

  interface ice_dyn_evp_1d_kernel
!    module procedure evp_kernel_v1
    module procedure evp_kernel_v2
  end interface 

  interface ice_dyn_evp_1d_copyout
    module procedure evp_copyout
  end interface 

  interface convert_2d_1d
!    module procedure convert_2d_1d_v1
    module procedure convert_2d_1d_v2
  end interface 

  integer(kind=int_kind) ::                                        &
    NA_len, NAVEL_len
  logical(kind=log_kind), dimension(:), allocatable ::             &
    skipucell
  integer(kind=int_kind), dimension(:), allocatable ::             &
    ee,ne,se,nw,sw,sse,indi,indj,indij , halo_parent
  real (kind=dbl_kind), dimension(:), allocatable ::               &
    cdn_ocn,aiu,uocn,vocn,forcex,forcey,Tbu,tarear,                &
    umassdti,fm,uarear,strintx,strinty,uvel_init,vvel_init
  real (kind=dbl_kind), dimension(:), allocatable ::               &
    strength,uvel,vvel,dxt,dyt,                                    &
!v1    dxhy,dyhx,cyp,cxp,cym,cxm,tinyarea,                            &
!v1    waterx,watery,                                                 &
    stressp_1, stressp_2, stressp_3, stressp_4,                    &
    stressm_1, stressm_2, stressm_3, stressm_4,                    &
    stress12_1,stress12_2,stress12_3,stress12_4,                   &
    divu,rdg_conv,rdg_shear,shear,taubx,tauby
  real (kind=DBL_KIND), dimension(:), allocatable ::               &
    str1, str2, str3, str4, str5, str6, str7, str8
  real (kind=dbl_kind), dimension(:), allocatable ::               &
    HTE,HTN,                                                       &
    HTEm1,HTNm1
  logical(kind=log_kind),parameter :: dbug = .false.


!--- dmi_omp ---------------------------
  interface domp_get_domain
    module procedure domp_get_domain_rlu
  end interface 

  INTEGER, PARAMETER :: JPIM = SELECTED_INT_KIND(9)
  integer(int_kind) :: domp_iam, domp_nt

#if defined (_OPENMP)
  ! Please note, this constant will create a compiler info for a constant
  ! expression in IF statements:
  real(kind=dbl_kind) :: rdomp_iam, rdomp_nt
  !$OMP THREADPRIVATE(domp_iam,domp_nt,rdomp_iam,rdomp_nt) 
#endif
!--- dmi_omp ---------------------------

!--- bench_v2 --------------------------
  interface evp1d_stress
    module procedure stress_i
    module procedure stress_l
  end interface 
  interface evp1d_stepu
    module procedure stepu_iter
    module procedure stepu_last
  end interface 
!--- bench_v2 --------------------------

  contains

!===============================================================================
!former module dmi_omp

  subroutine domp_init(nt_out)

#if defined (_OPENMP)
    use omp_lib, only : omp_get_thread_num, omp_get_num_threads
#endif

    integer(int_kind), intent(out) :: nt_out

    character(len=*), parameter :: subname = '(domp_init)'
    !---------------------------------------

    !$OMP PARALLEL DEFAULT(none)
#if defined (_OPENMP)
    domp_iam  = omp_get_thread_num()
    rdomp_iam = real(domp_iam,dbl_kind)
    domp_nt   = omp_get_num_threads()
    rdomp_nt  = real(domp_nt,dbl_kind)
#else
    domp_iam  = 0
    domp_nt   = 1
#endif
    !$OMP END PARALLEL

    if (dbug) then
#if defined (_OPENACC)
      write(nu_diag,'(2a)') subname,' Build with openACC support'
!#elif defined (_OPENMP)
!      write(nu_diag,'(2a)') subname,' Build with openMP support'
!#else
!      write(nu_diag,'(2a)') subname,' Build without openMP and openACC support'
#endif

      !- echo #threads:
      if (domp_nt > 1) then
        write(nu_diag,'(2a,i5,a)') subname,' Running openMP with ', domp_nt, ' threads'
      else
#if defined (_OPENMP)
        write(nu_diag,'(2a)') subname,' Running openMP with a single thread'
#else
        write(nu_diag,'(2a)') subname,' Running without openMP'
#endif
      endif
    endif

    !- return value of #threads:
    nt_out = domp_nt

  end subroutine domp_init
 
!----------------------------------------------------------------------------

  subroutine domp_get_domain_rlu(lower,upper,d_lower,d_upper)

#if defined (_OPENMP)
    use omp_lib,   only : omp_in_parallel
    use ice_constants, only: p5
#endif

    integer(KIND=JPIM), intent(in)  :: lower,upper
    integer(KIND=JPIM), intent(out) :: d_lower,d_upper

#if defined (_OPENMP)
    !-- local variables
    real(kind=dbl_kind)    :: dlen
#endif

    character(len=*), parameter :: subname = '(domp_get_domain_rlu)'
    !---------------------------------------

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
      dlen    = real(upper-lower+1, dbl_kind)
      d_lower = lower    + floor((rdomp_iam*dlen+p5)/rdomp_nt, JPIM)
      d_upper = lower -1 + floor((rdomp_iam*dlen+dlen+p5)/rdomp_nt, JPIM)
    endif
#endif

  if (.false.) then
    write(nu_diag,'(2a,i3,a,2i10)') subname,' openMP thread ', domp_iam, &
         ' handles range: ', d_lower, d_upper
  endif

  end subroutine domp_get_domain_rlu

!----------------------------------------------------------------------------

  subroutine domp_get_thread_no (tnum)

    implicit none
    integer(int_kind), intent(out) :: tnum
    character(len=*), parameter :: subname = '(domp_get_thread_no)'

    tnum = domp_iam + 1

  end subroutine domp_get_thread_no

!----------------------------------------------------------------------------

!former end module dmi_omp

!===============================================================================

!former module bench_v2

!----------------------------------------------------------------------------

    subroutine stress_i(NA_len, &
                     ee,ne,se,lb,ub,uvel,vvel,dxt,dyt,                           & 
                     hte,htn,htem1,htnm1,                                        &
                     strength,stressp_1,stressp_2,stressp_3,stressp_4,           &
                     stressm_1,stressm_2,stressm_3,stressm_4,stress12_1,         &
                     stress12_2,stress12_3,stress12_4,str1,str2,str3,str4,str5,  &
                     str6,str7,str8)

    use ice_kinds_mod
    use ice_constants, only: p027, p055, p111, p166, p222, p25, p333, p5, c1p5, c1
    use ice_dyn_shared, only: ecci, denom1, arlx1i, Ktens, revp

    implicit none

    integer (kind=int_kind), intent(in) :: NA_len
    integer (kind=int_kind), intent(in) :: lb,ub
    integer (kind=int_kind), dimension(:), intent(in),    contiguous :: ee,ne,se
    real    (kind=dbl_kind), dimension(:), intent(in),    contiguous ::          &
       strength, uvel, vvel, dxt, dyt,                                           &
       hte,htn,htem1,htnm1
    real    (kind=dbl_kind), dimension(:), intent(inout), contiguous ::          &
       stressp_1,stressp_2,  stressp_3, stressp_4, stressm_1,  stressm_2,        &
       stressm_3,stressm_4, stress12_1,stress12_2,stress12_3, stress12_4
    real    (kind=DBL_KIND), dimension(:), intent(out),   contiguous ::          &
       str1,str2,str3,str4,str5,str6,str7,str8

    !-- local variables

    integer (kind=int_kind) :: iw,il,iu
    real    (kind=dbl_kind) ::                                                   &
      puny
    real    (kind=DBL_KIND) ::                                                   &
      divune, divunw, divuse, divusw,tensionne, tensionnw, tensionse, tensionsw, & 
      shearne, shearnw, shearse, shearsw, Deltane, Deltanw, Deltase, Deltasw   , &
      c0ne, c0nw, c0se, c0sw, c1ne, c1nw, c1se, c1sw                           , &
      ssigpn, ssigps, ssigpe, ssigpw, ssigmn, ssigms, ssigme, ssigmw           , &
      ssig12n, ssig12s, ssig12e, ssig12w, ssigp1, ssigp2,ssigm1, ssigm2,ssig121, &
      ssig122, csigpne, csigpnw, csigpse, csigpsw,csigmne, csigmnw, csigmse    , &
      csigmsw, csig12ne, csig12nw, csig12se, csig12sw, str12ew, str12we,str12ns, &
      str12sn, strp_tmp, strm_tmp, tmp_uvel_ee, tmp_vvel_se, tmp_vvel_ee,        &
      tmp_vvel_ne, tmp_uvel_ne, tmp_uvel_se
    real    (kind=DBL_KIND) :: dxhy,dyhx,cxp,cyp,cxm,cym,tinyarea
   
    character(len=*), parameter :: subname = '(stress_i)'
    !---------------------------------------

    call icepack_query_parameters(puny_out=puny)
    call icepack_warnings_flush(nu_diag)
    if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
        file=__FILE__, line=__LINE__)

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
      ! shearing strain rate  =  2*e_12
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
  
      stressp_1(iw) = (stressp_1(iw)*(c1-arlx1i*revp) + c1ne*(divune*(c1+Ktens) - Deltane*(c1-Ktens))) * denom1
      stressp_2(iw) = (stressp_2(iw)*(c1-arlx1i*revp) + c1nw*(divunw*(c1+Ktens) - Deltanw*(c1-Ktens))) * denom1
      stressp_3(iw) = (stressp_3(iw)*(c1-arlx1i*revp) + c1sw*(divusw*(c1+Ktens) - Deltasw*(c1-Ktens))) * denom1
      stressp_4(iw) = (stressp_4(iw)*(c1-arlx1i*revp) + c1se*(divuse*(c1+Ktens) - Deltase*(c1-Ktens))) * denom1
  
      stressm_1(iw) = (stressm_1(iw)*(c1-arlx1i*revp) + c0ne*tensionne*(c1+Ktens)) * denom1
      stressm_2(iw) = (stressm_2(iw)*(c1-arlx1i*revp) + c0nw*tensionnw*(c1+Ktens)) * denom1
      stressm_3(iw) = (stressm_3(iw)*(c1-arlx1i*revp) + c0sw*tensionsw*(c1+Ktens)) * denom1
      stressm_4(iw) = (stressm_4(iw)*(c1-arlx1i*revp) + c0se*tensionse*(c1+Ktens)) * denom1
      
      stress12_1(iw) = (stress12_1(iw)*(c1-arlx1i*revp) + c0ne*shearne*p5*(c1+Ktens)) * denom1
      stress12_2(iw) = (stress12_2(iw)*(c1-arlx1i*revp) + c0nw*shearnw*p5*(c1+Ktens)) * denom1
      stress12_3(iw) = (stress12_3(iw)*(c1-arlx1i*revp) + c0sw*shearsw*p5*(c1+Ktens)) * denom1
      stress12_4(iw) = (stress12_4(iw)*(c1-arlx1i*revp) + c0se*shearse*p5*(c1+Ktens)) * denom1

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
#ifdef _OPENACC
    !$acc end parallel
#endif

  end subroutine stress_i
  
!----------------------------------------------------------------------------

  subroutine stress_l(NA_len, tarear, &
                     ee,ne,se,lb,ub,uvel,vvel,dxt,dyt,                           & 
                     hte,htn,htem1,htnm1,                                        &
                     strength,stressp_1,stressp_2,stressp_3,stressp_4,           & 
                     stressm_1,stressm_2,stressm_3,stressm_4,stress12_1,         &
                     stress12_2,stress12_3,stress12_4,                           &
                     divu,rdg_conv,rdg_shear,shear,                              &
                     str1,str2,str3,str4,str5,str6,str7,str8                     )

    use ice_kinds_mod
    use ice_constants, only: p027, p055, p111, p166, p222, p25, p333, p5, c1p5, c0, c1
    use ice_dyn_shared, only: ecci, denom1, arlx1i, Ktens, revp

    implicit none

    integer (kind=int_kind), intent(in) :: NA_len
    integer (kind=int_kind), intent(in) :: lb,ub
    integer (kind=int_kind), dimension(:), intent(in),    contiguous :: ee,ne,se
    real    (kind=dbl_kind), dimension(:), intent(in),    contiguous ::          &
       strength, uvel, vvel, dxt, dyt, tarear,                                   &
       hte,htn,htem1,htnm1
    real    (kind=dbl_kind), dimension(:), intent(inout), contiguous ::          &
       stressp_1,stressp_2,  stressp_3, stressp_4, stressm_1,  stressm_2,        &
       stressm_3,stressm_4, stress12_1,stress12_2,stress12_3, stress12_4
    real    (kind=DBL_KIND), dimension(:), intent(out),   contiguous ::          &
       str1,str2,str3,str4,str5,str6,str7,str8
    real    (kind=dbl_kind), dimension(:), intent(out),   contiguous ::          &
       divu,rdg_conv,rdg_shear,shear

    !-- local variables

    integer (kind=int_kind) :: iw,il,iu
    real    (kind=dbl_kind) ::                                                   &
      puny
    real    (kind=DBL_KIND) ::                                                   &
      divune, divunw, divuse, divusw,tensionne, tensionnw, tensionse, tensionsw, & 
      shearne, shearnw, shearse, shearsw, Deltane, Deltanw, Deltase, Deltasw   , &
      c0ne, c0nw, c0se, c0sw, c1ne, c1nw, c1se, c1sw                           , &
      ssigpn, ssigps, ssigpe, ssigpw, ssigmn, ssigms, ssigme, ssigmw           , &
      ssig12n, ssig12s, ssig12e, ssig12w, ssigp1, ssigp2,ssigm1, ssigm2,ssig121, &
      ssig122, csigpne, csigpnw, csigpse, csigpsw,csigmne, csigmnw, csigmse    , &
      csigmsw, csig12ne, csig12nw, csig12se, csig12sw, str12ew, str12we,str12ns, &
      str12sn, strp_tmp, strm_tmp, tmp_uvel_ee, tmp_vvel_se, tmp_vvel_ee,        &
      tmp_vvel_ne, tmp_uvel_ne, tmp_uvel_se
    real    (kind=DBL_KIND) :: dxhy,dyhx,cxp,cyp,cxm,cym,tinyarea
  
    character(len=*), parameter :: subname = '(stress_l)'
    !---------------------------------------

    call icepack_query_parameters(puny_out=puny)
    call icepack_warnings_flush(nu_diag)
    if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
        file=__FILE__, line=__LINE__)

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
  
      ! shearing strain rate  =  2*e_12
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
  
      !-----------------------------------------------------------------
      ! on last subcycle, save quantities for mechanical redistribution
      !-----------------------------------------------------------------
      divu(iw) = p25*(divune + divunw + divuse + divusw) * tarear(iw)
      rdg_conv(iw)  = -min(divu(iw),c0)    ! Could move outside the entire "kernel"
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

      stressp_1(iw) = (stressp_1(iw)*(c1-arlx1i*revp) + c1ne*(divune*(c1+Ktens) - Deltane*(c1-Ktens))) * denom1
      stressp_2(iw) = (stressp_2(iw)*(c1-arlx1i*revp) + c1nw*(divunw*(c1+Ktens) - Deltanw*(c1-Ktens))) * denom1
      stressp_3(iw) = (stressp_3(iw)*(c1-arlx1i*revp) + c1sw*(divusw*(c1+Ktens) - Deltasw*(c1-Ktens))) * denom1
      stressp_4(iw) = (stressp_4(iw)*(c1-arlx1i*revp) + c1se*(divuse*(c1+Ktens) - Deltase*(c1-Ktens))) * denom1
 
      stressm_1(iw) = (stressm_1(iw)*(c1-arlx1i*revp) + c0ne*tensionne*(c1+Ktens)) * denom1
      stressm_2(iw) = (stressm_2(iw)*(c1-arlx1i*revp) + c0nw*tensionnw*(c1+Ktens)) * denom1
      stressm_3(iw) = (stressm_3(iw)*(c1-arlx1i*revp) + c0sw*tensionsw*(c1+Ktens)) * denom1
      stressm_4(iw) = (stressm_4(iw)*(c1-arlx1i*revp) + c0se*tensionse*(c1+Ktens)) * denom1

      stress12_1(iw) = (stress12_1(iw)*(c1-arlx1i*revp) + c0ne*shearne*p5*(c1+Ktens)) * denom1
      stress12_2(iw) = (stress12_2(iw)*(c1-arlx1i*revp) + c0nw*shearnw*p5*(c1+Ktens)) * denom1
      stress12_3(iw) = (stress12_3(iw)*(c1-arlx1i*revp) + c0sw*shearsw*p5*(c1+Ktens)) * denom1
      stress12_4(iw) = (stress12_4(iw)*(c1-arlx1i*revp) + c0se*shearse*p5*(c1+Ktens)) * denom1
 
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
#ifdef _OPENACC
    !$acc end parallel
#endif
  end subroutine stress_l
  
!----------------------------------------------------------------------------

  subroutine stepu_iter(NA_len,rhow, &
                    lb,ub,Cw,aiu,uocn,vocn,forcex,forcey,umassdti,fm,uarear,Tbu, &
                    uvel_init,vvel_init,uvel,vvel,                               &
                    str1,str2,str3,str4,str5,str6,str7,str8, nw,sw,se,skipme)

    use ice_kinds_mod
    use ice_dyn_shared, only: brlx, revp
    use ice_constants, only: c0, c1

    implicit none

    integer (kind=int_kind), intent(in) :: NA_len
    real (kind=dbl_kind), intent(in) :: rhow
    integer(kind=int_kind),intent(in)   :: lb,ub
    logical(kind=log_kind),intent(in), dimension(:)   :: skipme
    integer(kind=int_kind),dimension(:), intent(in), contiguous :: nw,sw,se
    real(kind=dbl_kind),dimension(:), intent(in), contiguous    ::               &
       uvel_init, vvel_init, aiu, forcex, forcey, umassdti, Tbu,                 &
       uocn, vocn, fm, uarear,Cw
    real(kind=DBL_KIND),dimension(:), intent(in), contiguous ::                  &
       str1,str2,str3,str4,str5,str6,str7,str8
    real(kind=dbl_kind),dimension(:), intent(inout), contiguous ::               &
       uvel,vvel
    real (kind=dbl_kind), parameter :: &
           cosw = c1   , & ! cos(ocean turning angle)  ! turning angle = 0
           sinw = c0   

    !-- local variables

    integer (kind=int_kind) :: iw,il,iu
    real    (kind=dbl_kind) :: uold, vold, vrel,cca,ccb,ab2,cc1,cc2,taux,tauy,Cb
    real    (kind=dbl_kind) :: tmp_str2_nw,tmp_str3_se,tmp_str4_sw, tmp_strintx
    real    (kind=dbl_kind) :: tmp_str6_se,tmp_str7_nw,tmp_str8_sw, tmp_strinty
    real    (kind=dbl_kind) :: waterx,watery
    real    (kind=dbl_kind) :: u0 = 5.e-5_dbl_kind    ! residual velocity for seabed stress (m/s)
  
    character(len=*), parameter :: subname = '(stepu_iter)'
    !---------------------------------------

#ifdef _OPENACC
    !$acc parallel                                                       &
    !$acc present(Cw,aiu,uocn,vocn,forcex,forcey,umassdti,fm,uarear,Tbu, &
    !$acc         uvel_init,vvel_init,nw,sw,se,skipme,                   &
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
       Cb  = Tbu(iw) / (sqrt(uold**2 + vold**2) + u0)
       cca = (brlx + revp)*umassdti(iw) + vrel * cosw + Cb 
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
#ifdef _OPENACC
    !$acc end parallel
#endif

  end subroutine stepu_iter
  
!----------------------------------------------------------------------------

  subroutine stepu_last(NA_len, rhow, &
                    lb,ub,Cw,aiu,uocn,vocn,forcex,forcey,umassdti,fm,uarear,Tbu, &
                    strintx,strinty,taubx,tauby,                                 &
                    uvel_init,vvel_init,uvel,vvel,                               &
                    str1,str2,str3,str4,str5,str6,str7,str8, nw,sw,se,skipme)

    use ice_kinds_mod
    use ice_constants, only: c0, c1
    use ice_dyn_shared, only: brlx, revp, seabed_stress

    implicit none

    integer (kind=int_kind), intent(in) :: NA_len
    real (kind=dbl_kind), intent(in) :: rhow
    logical(kind=log_kind),intent(in), dimension(:)   :: skipme
    integer(kind=int_kind),intent(in)   :: lb,ub
    integer(kind=int_kind),dimension(:), intent(in), contiguous :: nw,sw,se
    real(kind=dbl_kind),dimension(:), intent(in), contiguous    ::               &
       uvel_init, vvel_init, aiu, forcex, forcey, umassdti, Tbu,                 &
       uocn, vocn, fm, uarear,Cw
    real(kind=DBL_KIND),dimension(:), intent(in), contiguous ::                  &
       str1,str2,str3,str4,str5,str6,str7,str8
    real(kind=dbl_kind),dimension(:), intent(inout), contiguous ::               &
       uvel,vvel, strintx,strinty, taubx,tauby
    real (kind=dbl_kind), parameter :: &
           cosw = c1   , & ! cos(ocean turning angle)  ! turning angle = 0
           sinw = c0   

    !-- local variables

    integer (kind=int_kind) :: iw,il,iu
    real    (kind=dbl_kind) :: uold, vold, vrel,cca,ccb,ab2,cc1,cc2,taux,tauy,Cb
    real    (kind=dbl_kind) :: tmp_str2_nw,tmp_str3_se,tmp_str4_sw
    real    (kind=dbl_kind) :: tmp_str6_se,tmp_str7_nw,tmp_str8_sw
    real    (kind=dbl_kind) :: waterx,watery
    real    (kind=dbl_kind) :: u0 = 5.e-5_dbl_kind    ! residual velocity for seabed stress (m/s)
   
    character(len=*), parameter :: subname = '(stepu_last)'
    !---------------------------------------

#ifdef _OPENACC
    !$acc parallel                                                                 &
    !$acc present(Cw,aiu,uocn,vocn,forcex,forcey,umassdti,fm,uarear,Tbu,           &
    !$acc         strintx,strinty,taubx,tauby,uvel_init,vvel_init,nw,sw,se,skipme, &
    !$acc         str1,str2,str3,str4,str5,str6,str7,str8,uvel,vvel                )
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
       Cb  = Tbu(iw) / (sqrt(uold**2 + vold**2) + u0)
       cca = (brlx + revp)*umassdti(iw) + vrel * cosw + Cb
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
       ! calculate seabed stress component for outputs
       if ( seabed_stress ) then
         taubx(iw) = -uvel(iw)*Tbu(iw) / (sqrt(uold**2 + vold**2) + u0)
         tauby(iw) = -vvel(iw)*Tbu(iw) / (sqrt(uold**2 + vold**2) + u0)
       endif
    enddo
#ifdef _OPENACC
    !$acc end parallel
#endif

  end subroutine stepu_last

!----------------------------------------------------------------------------

  subroutine evp1d_halo_update(NAVEL_len,lb,ub,uvel,vvel, halo_parent)

    use ice_kinds_mod

    implicit none

    integer (kind=int_kind), intent(in) :: NAVEL_len
    integer(kind=int_kind),intent(in)   :: lb,ub
    integer(kind=int_kind),dimension(:), intent(in), contiguous :: halo_parent
    real(kind=dbl_kind),dimension(:), intent(inout), contiguous :: uvel,vvel

    !-- local variables

    integer (kind=int_kind) :: iw,il,iu

    character(len=*), parameter :: subname = '(evp1d_halo_update)'
    !---------------------------------------

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
#ifdef _OPENACC
    !$acc end parallel
#endif

  end subroutine evp1d_halo_update

!----------------------------------------------------------------------------

!former end module bench_v2
  
!===============================================================================
!----------------------------------------------------------------------------

  subroutine alloc1d(na)

    implicit none

    integer(kind=int_kind),intent(in) :: na
    integer(kind=int_kind)            :: ierr,nb

    character(len=*), parameter :: subname = '(alloc1d)'
    !---------------------------------------

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
!v1      dxhy(1:na),dyhx(1:na),cyp(1:na),cxp(1:na),cym(1:na),cxm(1:na),tinyarea(1:na),&
      strength(1:na),dxt(1:na),dyt(1:na), tarear(1:na),                   &
      stressp_1(1:na), stressp_2(1:na), stressp_3(1:na), stressp_4(1:na), &
      stressm_1(1:na), stressm_2(1:na), stressm_3(1:na), stressm_4(1:na), &
      stress12_1(1:na),stress12_2(1:na),stress12_3(1:na),stress12_4(1:na),&
      divu(1:na),rdg_conv(1:na),rdg_shear(1:na),shear(1:na),              &
      ! U cells
!v1      waterx(1:nb),watery(1:nb),                                          &
      cdn_ocn(1:nb),aiu(1:nb),uocn(1:nb),vocn(1:nb),                      &
      forcex(1:nb),forcey(1:nb),Tbu(1:nb),                                &
      umassdti(1:nb),fm(1:nb),uarear(1:nb),                               &
      strintx(1:nb),strinty(1:nb),                                        &
      uvel_init(1:nb),vvel_init(1:nb),                                    &
      taubx(1:nb),tauby(1:nb),                                            &
      stat=ierr)

    if (ierr/=0) call abort_ice(subname//': ERROR allocating 1D')

  end subroutine alloc1d

!----------------------------------------------------------------------------

  subroutine alloc1d_navel(navel)

    implicit none

    integer(kind=int_kind),intent(in) :: navel
    integer(kind=int_kind)            :: ierr

    character(len=*), parameter :: subname = '(alloc1d_navel)'
    !---------------------------------------

    allocate(                                                            &
      uvel(1:navel),vvel(1:navel), indij(1:navel), halo_parent(1:navel), &
      str1(1:navel),str2(1:navel),str3(1:navel),str4(1:navel),           &
      str5(1:navel),str6(1:navel),str7(1:navel),str8(1:navel),           &
      stat=ierr)
    if (ierr/=0) call abort_ice(subname// ': Error allocating 1D navel')

  end subroutine alloc1d_navel

!----------------------------------------------------------------------------

  subroutine dealloc1d

    implicit none

    integer(kind=int_kind) :: ierr

    character(len=*), parameter :: subname = '(dealloc1d)'
    !---------------------------------------

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
      forcex,forcey,Tbu,                          &
      umassdti,fm,uarear,                         &
      strintx,strinty,                            &
      uvel_init,vvel_init,                        &
      taubx,tauby,                                &
      ! NAVEL 
      uvel,vvel, indij, halo_parent,              &
      stat=ierr)

    if (ierr/=0) call abort_ice(subname//': Error de-allocating 1D')

!v1    if (allocated(tinyarea)) then
!v1      deallocate(                                 &
!v1        dxhy,dyhx,cyp,cxp,cym,cxm,tinyarea,       &
!v1        waterx,watery,                            &
!v1        stat=ierr)
!v1      if (ierr/=0) call abort_ice(subname//': Error de-allocating 1D, v1')
!v1    endif

    if (allocated(HTE)) then
      deallocate(                                 &
        ! Grid distances: HTE,HTN + "-1 neighbours" 
        HTE,HTN, HTEm1,HTNm1,                     &
        stat=ierr)
      if (ierr/=0) call abort_ice(subname//': Error de-allocating 1D, v2')
    endif

  end subroutine dealloc1d

!----------------------------------------------------------------------------

  subroutine evp_copyin_v2(nx,ny,nblk,nx_glob,ny_glob,                                &
       I_HTE,I_HTN,                                                                   &
!v1       I_dxhy,I_dyhx,I_cyp,I_cxp,I_cym,I_cxm,I_tinyarea,                              &
!v1       I_waterx,I_watery,                                                             &
       I_icetmask,I_iceumask,                                                         &
       I_cdn_ocn,I_aiu,I_uocn,I_vocn,I_forcex,I_forcey,I_Tbu,                         &
       I_umassdti,I_fm,I_uarear,I_tarear,I_strintx,I_strinty,I_uvel_init,I_vvel_init, &
       I_strength,I_uvel,I_vvel,I_dxt,I_dyt,                                          &
        I_stressp_1 ,I_stressp_2, I_stressp_3, I_stressp_4,                           &
        I_stressm_1 ,I_stressm_2, I_stressm_3, I_stressm_4,                           &
       I_stress12_1,I_stress12_2,I_stress12_3,I_stress12_4                            )

    use ice_gather_scatter, only: gather_global_ext
    use ice_domain, only: distrb_info
    use ice_communicate, only: my_task, master_task

    implicit none

    integer(int_kind), intent(in) :: nx, ny, nblk, nx_glob, ny_glob
    integer (kind=int_kind),dimension (nx,ny,nblk), intent(in) :: I_icetmask
    logical (kind=log_kind),dimension (nx,ny,nblk), intent(in) :: I_iceumask
    real (kind=dbl_kind), dimension(nx,ny,nblk), intent(in)  ::                       &
       I_HTE,I_HTN,                                                                   &
!v1       I_dxhy,I_dyhx,I_cyp,I_cxp,I_cym,I_cxm,I_tinyarea,                              &
!v1       I_waterx,I_watery,                                                             &
       I_cdn_ocn,I_aiu,I_uocn,I_vocn,I_forcex,I_forcey,I_Tbu,                         &
       I_umassdti,I_fm,I_uarear,I_tarear,I_strintx,I_strinty,I_uvel_init,I_vvel_init, &
       I_strength,I_uvel,I_vvel,I_dxt,I_dyt,                                          &
       I_stressp_1, I_stressp_2, I_stressp_3, I_stressp_4,                            &
       I_stressm_1, I_stressm_2, I_stressm_3, I_stressm_4,                            &
       I_stress12_1,I_stress12_2,I_stress12_3,I_stress12_4

    !-- local variables

    integer (kind=int_kind),dimension (nx_glob,ny_glob) :: G_icetmask
    logical (kind=log_kind),dimension (nx_glob,ny_glob) :: G_iceumask
    real (kind=dbl_kind), dimension(nx_glob,ny_glob) ::                               &
       G_HTE,G_HTN,                                                                   &
!v1       G_dxhy,G_dyhx,G_cyp,G_cxp,G_cym,G_cxm,G_tinyarea,                              &
!v1       G_waterx,G_watery,                                                             &
       G_cdn_ocn,G_aiu,G_uocn,G_vocn,G_forcex,G_forcey,G_Tbu,                         &
       G_umassdti,G_fm,G_uarear,G_tarear,G_strintx,G_strinty,G_uvel_init,G_vvel_init, &
       G_strength,G_uvel,G_vvel,G_dxt,G_dyt,                                          &
       G_stressp_1, G_stressp_2, G_stressp_3, G_stressp_4,                            &
       G_stressm_1, G_stressm_2, G_stressm_3, G_stressm_4,                            &
       G_stress12_1,G_stress12_2,G_stress12_3,G_stress12_4
    integer(kind=int_kind) :: na, navel

    character(len=*), parameter :: subname = '(evp_copyin_v2)'
    !---------------------------------------
    !-- Gather data into one single block --

    call gather_global_ext(G_icetmask, I_icetmask, master_task, distrb_info)
    call gather_global_ext(G_iceumask, I_iceumask, master_task, distrb_info)
    call gather_global_ext(G_HTE, I_HTE, master_task, distrb_info)
    call gather_global_ext(G_HTN, I_HTN, master_task, distrb_info)
!v1    call gather_global_ext(G_dxhy, I_dxhy, master_task, distrb_info)
!v1    call gather_global_ext(G_dyhx, I_dyhx, master_task, distrb_info)
!v1    call gather_global_ext(G_cyp, I_cyp, master_task, distrb_info)
!v1    call gather_global_ext(G_cxp, I_cxp, master_task, distrb_info)
!v1    call gather_global_ext(G_cym, I_cym, master_task, distrb_info)
!v1    call gather_global_ext(G_cxm, I_cxm, master_task, distrb_info)
!v1    call gather_global_ext(G_tinyarea, I_tinyarea, master_task, distrb_info)
!v1    call gather_global_ext(G_waterx, I_waterx, master_task, distrb_info)
!v1    call gather_global_ext(G_watery, I_watery, master_task, distrb_info)
    call gather_global_ext(G_cdn_ocn, I_cdn_ocn, master_task, distrb_info)
    call gather_global_ext(G_aiu, I_aiu, master_task, distrb_info)
    call gather_global_ext(G_uocn, I_uocn, master_task, distrb_info)
    call gather_global_ext(G_vocn, I_vocn, master_task, distrb_info)
    call gather_global_ext(G_forcex, I_forcex, master_task, distrb_info)
    call gather_global_ext(G_forcey, I_forcey, master_task, distrb_info)
    call gather_global_ext(G_Tbu, I_Tbu, master_task, distrb_info)
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

    !-- All calculations has to be done on the master-task --

    if (my_task == master_task) then
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
         G_HTE,G_HTN,                                                                     &
!v1         G_dxhy,G_dyhx,G_cyp,G_cxp,G_cym,G_cxm,G_tinyarea,                                &
!v1         G_waterx,G_watery,                                                               &
         G_cdn_ocn,G_aiu,G_uocn,G_vocn,G_forcex,G_forcey,G_Tbu,                           &
         G_umassdti,G_fm,G_uarear,G_tarear,G_strintx,G_strinty,G_uvel_init,G_vvel_init,   &
         G_strength,G_uvel,G_vvel,G_dxt,G_dyt,                                            &
          G_stressp_1, G_stressp_2, G_stressp_3, G_stressp_4,                             &
          G_stressm_1, G_stressm_2, G_stressm_3, G_stressm_4,                             &
         G_stress12_1,G_stress12_2,G_stress12_3,G_stress12_4                              )
      call calc_halo_parent(nx_glob,ny_glob,na,navel, G_icetmask)
      NA_len=na
      NAVEL_len=navel
    endif

    !-- write check
!if (1 == 1) then
!  write(nu_diag,*) subname,' MHRI: INDICES start:'
!  write(nu_diag,*) 'na,navel  ', na,navel
!  write(nu_diag,*) 'Min/max ee', minval(ee(1:na)), maxval(ee(1:na))
!  write(nu_diag,*) 'Min/max ne', minval(ne(1:na)), maxval(ne(1:na))
!  write(nu_diag,*) 'Min/max se', minval(se(1:na)), maxval(se(1:na))
!  write(nu_diag,*) 'Min/max nw', minval(nw(1:na)), maxval(nw(1:na))
!  write(nu_diag,*) 'Min/max sw', minval(sw(1:na)), maxval(sw(1:na))
!  write(nu_diag,*) 'Min/max sse', minval(sse(1:na)), maxval(sse(1:na))
!  write(nu_diag,*) subname,' MHRI: INDICES end:'
!endif

  end subroutine evp_copyin_v2

!----------------------------------------------------------------------------

  subroutine evp_copyout(nx,ny,nblk,nx_glob,ny_glob,                    &
               I_uvel,I_vvel, I_strintx,I_strinty,                      &
               I_stressp_1, I_stressp_2, I_stressp_3, I_stressp_4,      &
               I_stressm_1, I_stressm_2, I_stressm_3, I_stressm_4,      &
               I_stress12_1,I_stress12_2,I_stress12_3,I_stress12_4,     &
               I_divu,I_rdg_conv,I_rdg_shear,I_shear,I_taubx,I_tauby    )

    use ice_constants, only : c0
    use ice_gather_scatter, only: scatter_global_ext
    use ice_domain, only: distrb_info
    use ice_communicate, only: my_task, master_task

    implicit none

    integer(int_kind), intent(in) :: nx,ny,nblk, nx_glob,ny_glob
    real(dbl_kind), dimension(nx,ny,nblk), intent(out) ::       &
       I_uvel,I_vvel, I_strintx,I_strinty,                      &
       I_stressp_1, I_stressp_2, I_stressp_3, I_stressp_4,      &
       I_stressm_1, I_stressm_2, I_stressm_3, I_stressm_4,      &
       I_stress12_1,I_stress12_2,I_stress12_3,I_stress12_4,     &
       I_divu,I_rdg_conv, I_rdg_shear,I_shear, I_taubx,I_tauby

    !-- local variables

    real(dbl_kind), dimension(nx_glob,ny_glob) :: &
       G_uvel,G_vvel, G_strintx,G_strinty,                      &
       G_stressp_1, G_stressp_2, G_stressp_3, G_stressp_4,      &
       G_stressm_1, G_stressm_2, G_stressm_3, G_stressm_4,      &
       G_stress12_1,G_stress12_2,G_stress12_3,G_stress12_4,     &
       G_divu,G_rdg_conv, G_rdg_shear,G_shear, G_taubx,G_tauby
    integer(int_kind) :: i,j,iw, nx_block

    character(len=*), parameter :: subname = '(evp_copyout)'
    !---------------------------------------
    ! Remap 1d to 2d and fill in
    nx_block=nx_glob  ! Total block size in x-dir 

    if (my_task == master_task) then
      G_uvel       = c0
      G_vvel       = c0
      G_strintx    = c0
      G_strinty    = c0
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
      G_taubx      = c0
      G_tauby      = c0
      !$OMP PARALLEL PRIVATE(iw,i,j)
      do iw=1,NAVEL_len
        j=int((indij(iw)-1)/(nx_block))+1
        i=indij(iw)-(j-1)*nx_block
        G_uvel(i,j)       = uvel(iw)
        G_vvel(i,j)       = vvel(iw)
      enddo
      !$OMP END PARALLEL
      !$OMP PARALLEL PRIVATE(iw,i,j)
      do iw=1,NA_len
        i=indi(iw)
        j=indj(iw)
!        G_uvel(i,j)       = uvel(iw)  ! done above
!        G_vvel(i,j)       = vvel(iw)  ! done above
        G_strintx(i,j)    = strintx(iw)
        G_strinty(i,j)    = strinty(iw)
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
        G_taubx(i,j)      = taubx(iw)
        G_tauby(i,j)      = tauby(iw)
      enddo
      !$OMP END PARALLEL
      call dealloc1d()
    endif

    !-- Scatter data into blocks --
    !--   has to be done on all tasks --

    call scatter_global_ext(I_uvel, G_uvel, master_task, distrb_info)
    call scatter_global_ext(I_vvel, G_vvel, master_task, distrb_info)
    call scatter_global_ext(I_strintx, G_strintx, master_task, distrb_info)
    call scatter_global_ext(I_strinty, G_strinty, master_task, distrb_info)
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
    call scatter_global_ext(I_taubx, G_taubx, master_task, distrb_info)
    call scatter_global_ext(I_tauby, G_tauby, master_task, distrb_info)

  end subroutine evp_copyout

!----------------------------------------------------------------------------

  subroutine evp_kernel_v2

    use ice_constants, only : c0
    use ice_dyn_shared, only: ndte
    use ice_communicate, only: my_task, master_task
    implicit none

    real(kind=dbl_kind) :: rhow
    integer (kind=int_kind) :: i, nthreads
    integer (kind=int_kind) :: na,nb,navel

    character(len=*), parameter :: subname = '(evp_kernel_v2)'
    !---------------------------------------
    !-- All calculations has to be done on one single node (choose master-task) --

    if (my_task == master_task) then

      !- Read constants...
      call icepack_query_parameters(rhow_out=rhow)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__, line=__LINE__)
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
   
      if (ndte<2) call abort_ice(subname//' ERROR: ndte must be 2 or higher for this kernel')

      !$OMP PARALLEL PRIVATE(i)
      do i = 1, ndte-1
        call evp1d_stress(NA_len,                                                    &
                     ee,ne,se,1,na,uvel,vvel,dxt,dyt,                                &
                     hte,htn,htem1,htnm1,                                            &
                     strength,stressp_1,stressp_2,stressp_3,stressp_4,               &
                     stressm_1,stressm_2,stressm_3,stressm_4,stress12_1,             &
                     stress12_2,stress12_3,stress12_4,str1,str2,str3,                &
                     str4,str5,str6,str7,str8)
        !$OMP BARRIER
        call evp1d_stepu(NA_len, rhow,                                               &
                     1,nb,cdn_ocn,aiu,uocn,vocn,forcex,forcey,umassdti,fm,uarear,Tbu,&
                     uvel_init,vvel_init,uvel,vvel,                                  &
                     str1,str2,str3,str4,str5,str6,str7,str8, nw,sw,sse,skipucell)
        !$OMP BARRIER
        call evp1d_halo_update(NA_len,1,navel,uvel,vvel, halo_parent)
        !$OMP BARRIER
      enddo

      call evp1d_stress(NA_len, tarear,                                              &
                     ee,ne,se,1,na,uvel,vvel,dxt,dyt,                                &
                     hte,htn,htem1,htnm1,                                            &
                     strength,stressp_1,stressp_2,stressp_3,stressp_4,               &
                     stressm_1,stressm_2,stressm_3,stressm_4,stress12_1,             &
                     stress12_2,stress12_3,stress12_4,                               &
                     divu,rdg_conv,rdg_shear,shear,                                  &
                     str1,str2,str3,str4,str5,str6,str7,str8)
      !$OMP BARRIER
      call evp1d_stepu(NA_len, rhow,                                                 &
                     1,nb,cdn_ocn,aiu,uocn,vocn,forcex,forcey,umassdti,fm,uarear,Tbu,&
                     strintx,strinty,taubx,tauby,                                    &
                     uvel_init,vvel_init,uvel,vvel,                                  &
                     str1,str2,str3,str4,str5,str6,str7,str8, nw,sw,sse,skipucell)
      !$OMP BARRIER
      call evp1d_halo_update(NA_len,1,navel,uvel,vvel, halo_parent)
      !$OMP END PARALLEL

    endif

  end subroutine evp_kernel_v2

!----------------------------------------------------------------------------

  subroutine calc_na(nx,ny,na,icetmask)
    ! Calculate number of active points (na)
    use ice_blocks, only: nghost

    implicit none

    integer(int_kind),intent(in) :: nx,ny
    integer(int_kind),intent(out) :: na
    integer (kind=int_kind),dimension (nx,ny), intent(in) :: icetmask
    integer(int_kind) :: i,j

    character(len=*), parameter :: subname = '(calc_na)'
    !---------------------------------------

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

!----------------------------------------------------------------------------

  subroutine calc_2d_indices(nx,ny,na,icetmask,iceumask)

    use ice_blocks, only: nghost

    implicit none

    integer(int_kind),intent(in) :: nx,ny,na
    integer (kind=int_kind),dimension (nx,ny), intent(in) :: icetmask
    logical (kind=log_kind),dimension (nx,ny), intent(in) :: iceumask
    integer(int_kind) :: i,j,Nmaskt

    character(len=*), parameter :: subname = '(calc_2d_indices)'
    !---------------------------------------

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
      write(nu_diag,*) subname,' Nmaskt,na: ',Nmaskt,na
      call abort_ice(subname//': ERROR Problem Nmaskt != na')
    endif
    if (Nmaskt==0) then
      write(nu_diag,*) subname,' WARNING: NO ICE'
    endif

  end subroutine calc_2d_indices

!----------------------------------------------------------------------------

  subroutine calc_navel(nx_block,ny_block,na,navel)
    ! Calculate number of active points including needed halo points (navel)

    implicit none

    integer(int_kind),intent(in)  :: nx_block,ny_block,na
    integer(int_kind),intent(out) :: navel

    integer(int_kind) :: iw,i,j
    integer(int_kind),dimension(1:na) :: Iin,Iee,Ine,Ise,Inw,Isw,Isse
    integer(int_kind),dimension(1:7*na) :: util1,util2

    character(len=*), parameter :: subname = '(calc_navel)'

    !---------------------------------------
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
        write(nu_diag,*) subname,' nx_block,ny_block,nx_block*ny_block: ',nx_block,ny_block,nx_block*ny_block
        write(nu_diag,*) subname,' na,navel,iw,util2(iw): ',na,navel,iw,util2(iw)
        call abort_ice(subname//': Problem with boundary. Check halo zone values')
      endif
    enddo

  end subroutine calc_navel

!----------------------------------------------------------------------------

  subroutine convert_2d_1d_v2(nx,ny, na,navel,                &
       I_HTE,I_HTN,                                           &
!v1       I_dxhy,I_dyhx,I_cyp,I_cxp,I_cym,I_cxm,I_tinyarea,      &
!v1       I_waterx,I_watery,                                     &
       I_cdn_ocn,I_aiu,I_uocn,I_vocn,I_forcex,I_forcey,I_Tbu, &
       I_umassdti,I_fm,I_uarear,I_tarear,I_strintx,I_strinty, &
       I_uvel_init,I_vvel_init,                               &
       I_strength,I_uvel,I_vvel,I_dxt,I_dyt,                  &
        I_stressp_1 ,I_stressp_2, I_stressp_3, I_stressp_4,   &
        I_stressm_1 ,I_stressm_2, I_stressm_3, I_stressm_4,   &
       I_stress12_1,I_stress12_2,I_stress12_3,I_stress12_4    )

    implicit none

    integer(int_kind),intent(in) :: nx,ny,na,navel
    real (kind=dbl_kind), dimension(nx,ny), intent(in)    ::  &
       I_HTE,I_HTN,                                           &
!v1       I_dxhy,I_dyhx,I_cyp,I_cxp,I_cym,I_cxm,I_tinyarea,      &
!v1       I_waterx,I_watery,                                     &
       I_cdn_ocn,I_aiu,I_uocn,I_vocn,I_forcex,I_forcey,I_Tbu, &
       I_umassdti,I_fm,I_uarear,I_tarear,I_strintx,I_strinty, &
       I_uvel_init,I_vvel_init,                               &
       I_strength,I_uvel,I_vvel,I_dxt,I_dyt,                  &
        I_stressp_1 ,I_stressp_2, I_stressp_3, I_stressp_4,   &
        I_stressm_1 ,I_stressm_2, I_stressm_3, I_stressm_4,   &
       I_stress12_1,I_stress12_2,I_stress12_3,I_stress12_4

    integer(int_kind) :: iw,i,j, nx_block
    integer(int_kind),dimension(1:na) :: Iin,Iee,Ine,Ise,Inw,Isw,Isse
    integer(int_kind),dimension(1:7*na) :: util1,util2
    integer(int_kind) :: nachk

    character(len=*), parameter :: subname = '(convert_2d_1d_v2)'

    !---------------------------------------
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
      write(nu_diag,*) subname,' ERROR: navel badly chosen: na,navel,nachk = ',na,navel,nachk
      call abort_ice(subname//': ERROR: navel badly chosen')
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
!    write(nu_diag,*) subname,' MHRI: INDICES start:'
!    write(nu_diag,*) 'Min/max ee', minval(ee), maxval(ee)
!    write(nu_diag,*) 'Min/max ne', minval(ne), maxval(ne)
!    write(nu_diag,*) 'Min/max se', minval(se), maxval(se)
!    write(nu_diag,*) 'Min/max nw', minval(nw), maxval(nw)
!    write(nu_diag,*) 'Min/max sw', minval(sw), maxval(sw)
!    write(nu_diag,*) 'Min/max sse',minval(sse),maxval(sse)
!    write(nu_diag,*) subname,' MHRI: INDICES end:'
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
      forcex(iw)=    I_forcex(i,j)   
      forcey(iw)=    I_forcey(i,j)
      Tbu(iw)=       I_Tbu(i,j)
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
!v1      dxhy(iw)=      I_dxhy(i,j)   
!v1      dyhx(iw)=      I_dyhx(i,j)   
!v1      cyp(iw)=       I_cyp(i,j)   
!v1      cxp(iw)=       I_cxp(i,j)   
!v1      cym(iw)=       I_cym(i,j)   
!v1      cxm(iw)=       I_cxm(i,j)   
!v1      tinyarea(iw)=  I_tinyarea(i,j)   
!v1      waterx(iw)=    I_waterx(i,j)   
!v1      watery(iw)=    I_watery(i,j)   
      HTE(iw)   = I_HTE(i,j)
      HTN(iw)   = I_HTN(i,j)
      HTEm1(iw) = I_HTE(i-1,j)
      HTNm1(iw) = I_HTN(i,j-1)
    enddo
    !$OMP END PARALLEL DO

  end subroutine convert_2d_1d_v2

!----------------------------------------------------------------------------

  subroutine calc_halo_parent(nx,ny,na,navel, I_icetmask)

    implicit none

    integer(kind=int_kind),intent(in) :: nx,ny,na,navel
    integer(kind=int_kind), dimension(nx,ny), intent(in) :: I_icetmask

    integer(kind=int_kind) :: iw,i,j !,masku,maskt
    integer(kind=int_kind),dimension(1:navel) :: Ihalo

    character(len=*), parameter :: subname = '(calc_halo_parent)'

    !---------------------------------------
    ! Indices for halo update:
    !      0: no halo point
    !     >0: index for halo point parent. Finally related to indij vector
    ! TODO: ONLY for nghost==1
    ! TODO: ONLY for circular grids - NOT tripole grids

    Ihalo(:)=0
    halo_parent(:)=0

    !$OMP PARALLEL DO PRIVATE(iw,i,j)
    do iw=1,navel
      j=int((indij(iw)-1)/(nx))+1
      i=indij(iw)-(j-1)*nx
      ! If within ghost-zone:
      if (i==nx .and. I_icetmask(   2,j)==1) Ihalo(iw)=     2+ (j-1)*nx
      if (i==1  .and. I_icetmask(nx-1,j)==1) Ihalo(iw)=(nx-1)+ (j-1)*nx
      if (j==ny .and. I_icetmask(i,   2)==1) Ihalo(iw)=     i+       nx
      if (j==1  .and. I_icetmask(i,ny-1)==1) Ihalo(iw)=     i+(ny-2)*nx
    enddo
    !$OMP END PARALLEL DO

    ! Relate halo indices to indij vector
    call findXinY_halo(Ihalo,indij,navel,navel,halo_parent)

    !-- write check
!if (1 == 1) then
!    integer(kind=int_kind) :: iiw,ii,jj !,masku,maskt  MHRI
!    write(nu_diag,*) subname,' MHRI: halo boundary start:'
!    do iw=1,navel
!      if (halo_parent(iw)>0) then
!      iiw=halo_parent(iw)
!      j=int((indij(iiw)-1)/(nx))+1
!      i=indij(iiw)-(j-1)*nx
!      ii=i
!      jj=j
!      j=int((indij(iw)-1)/(nx))+1
!      i=indij(iw)-(j-1)*nx
!      write(nu_diag,*)iw,i,j,iiw,ii,jj
!      endif
!    enddo
!    write(nu_diag,*) subname,' MHRI: halo boundary end:'
!endif

  end subroutine calc_halo_parent

!----------------------------------------------------------------------------

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

    character(len=*), parameter :: subname = '(union)'

    !---------------------------------------

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

!----------------------------------------------------------------------------

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

    character(len=*), parameter :: subname = '(setdiff)'
    !---------------------------------------

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

!----------------------------------------------------------------------------

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

    character(len=*), parameter :: subname = '(findXinY)'
    !---------------------------------------

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
        write(nu_diag,*) subname,' nx,ny: ',nx,ny
        write(nu_diag,*) subname,' i,j1,j2: ',i,j1,j2
        write(nu_diag,*) subname,' x(i),y(j1),y(j2): ',x(i),y(j1),y(j2)
        call abort_ice(subname//': ERROR')
      endif
    end do

  end subroutine findXinY

!----------------------------------------------------------------------------

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

    character(len=*), parameter :: subname = '(findXinY_halo)'
    !---------------------------------------

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
            write(nu_diag,*) subname,' nx,ny: ',nx,ny
            write(nu_diag,*) subname,' i,j1: ',i,j1
            write(nu_diag,*) subname,' x(i),y(j1): ',x(i),y(j1)
            call abort_ice(subname//': ERROR too many loops')
          endif
        endif
      endif
    end do

  end subroutine findXinY_halo

!----------------------------------------------------------------------------

  subroutine numainit(l,u,uu)

    use ice_constants, only: c0

    implicit none

    integer(kind=int_kind),intent(in) :: l,u,uu

    integer(kind=int_kind) :: lo,up

    character(len=*), parameter :: subname = '(numainit)'
    !---------------------------------------

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
    uvel_init(lo:up)=c0
    vvel_init(lo:up)=c0
    uocn(lo:up)=c0
    vocn(lo:up)=c0
    dxt(lo:up)=c0
    dyt(lo:up)=c0
    HTE(lo:up)=c0
    HTN(lo:up)=c0
    HTEm1(lo:up)=c0
    HTNm1(lo:up)=c0
!v1    dxhy(lo:up)=c0
!v1    dyhx(lo:up)=c0
!v1    cyp(lo:up)=c0
!v1    cxp(lo:up)=c0
!v1    cym(lo:up)=c0
!v1    cxm(lo:up)=c0
!v1    tinyarea(lo:up)=c0
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
    Tbu(lo:up)=c0
    taubx(lo:up)=c0
    tauby(lo:up)=c0
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
  
!----------------------------------------------------------------------------
!===============================================================================

end module ice_dyn_evp_1d

