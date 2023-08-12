!===============================================================================
! Copyright (C) 2023, Intel Corporation
! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
! 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
! 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
! 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!===============================================================================

!===============================================================================
!
! Elastic-viscous-plastic sea ice dynamics model
! Computes ice velocity and deformation
!
! See:
!
! Hunke, E. C., and J. K. Dukowicz (1997). An elastic-viscous-plastic model
! for sea ice dynamics. J. Phys. Oceanogr., 27, 1849-1867.
!
! Hunke, E. C. (2001).  Viscous-Plastic Sea Ice Dynamics with the EVP Model:
! Linearization Issues. J. Comput. Phys., 170, 18-38.
!
! Hunke, E. C., and J. K. Dukowicz (2002).  The Elastic-Viscous-Plastic
! Sea Ice Dynamics Model in General Orthogonal Curvilinear Coordinates
! on a Sphere - Incorporation of Metric Terms. Mon. Weather Rev.,
! 130, 1848-1865.
!
! Hunke, E. C., and J. K. Dukowicz (2003).  The sea ice momentum
! equation in the free drift regime.  Los Alamos Tech. Rep. LA-UR-03-2219.
!
! Hibler, W. D. (1979). A dynamic thermodynamic sea ice model. J. Phys.
! Oceanogr., 9, 817-846.
!
! Bouillon, S., T. Fichefet, V. Legat and G. Madec (2013).  The
! elastic-viscous-plastic method revisited.  Ocean Model., 71, 2-12.
!
! author: Elizabeth C. Hunke, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb (LANL)
! 2004: Block structure added by William Lipscomb
! 2005: Removed boundary calls for stress arrays (WHL)
! 2006: Streamlined for efficiency by Elizabeth Hunke
!       Converted to free source form (F90)
!===============================================================================
! 2023: Intel
!       Refactored to support explicit inlining
!       Refactored to allow private subroutines in stress to become pure
!===============================================================================
!===============================================================================
! 2023: DMI
!       Updated to match requirements from CICE
!===============================================================================

module bench_v2b
contains
!DIR$ ATTRIBUTES FORCEINLINE :: stress
subroutine stress (ee, ne, se, lb, ub,                                         &
                   uvel, vvel, dxT, dyT, skipme, strength,                     &
                   hte, htn, htem1, htnm1,                                     &
                   stressp_1,  stressp_2,  stressp_3,  stressp_4,              &
                   stressm_1,  stressm_2,  stressm_3,  stressm_4,              &
                   stress12_1, stress12_2, stress12_3, stress12_4,             &
                   str1, str2, str3, str4, str5, str6, str7, str8)
  use ice_kinds_mod
  use myomp, only : domp_get_domain
  use ice_constants, only: c1, p027, p055, p111, p166, c1p5,        &
                           p2, p222, p25, p333, p5 
  use ice_dyn_shared, only: arlx1i, denom1, Ktens, revp, capping,    &
                           deltaminEVP, epp2i, e_factor
! 
  implicit none
  ! arguments ------------------------------------------------------------------
  integer (kind=int_kind), intent(in)                           :: lb,ub
  integer (kind=int_kind), dimension(:), intent(in), contiguous :: ee,ne,se
  logical (kind=log_kind), dimension(:), intent(in), contiguous :: skipme
  real    (kind=dbl_kind), dimension(:), intent(in), contiguous ::             &
           strength , & ! ice strength (N/m)
           uvel     , & ! x-component of velocity (m/s)
           vvel     , & ! y-component of velocity (m/s)
           dxT      , & ! width of T-cell through the middle (m)
           dyT      , & ! height of T-cell through the middle (m)
           hte      , &
           htn      , &
           htem1    , &
           htnm1
  real    (kind=dbl_kind), dimension(:), intent(inout), contiguous ::          &
           stressp_1, stressp_2, stressp_3, stressp_4 , & ! sigma11+sigma22
           stressm_1, stressm_2, stressm_3, stressm_4 , & ! sigma11-sigma22
           stress12_1,stress12_2,stress12_3,stress12_4    ! sigma12
  real    (kind=dbl_kind), dimension(:), intent(out), contiguous ::            &
           str1,str2,str3,str4,str5,str6,str7,str8
  ! local variables
  integer (kind=int_kind) :: iw,il,iu
  real    (kind=dbl_kind) ::                          &
    divune, divunw, divuse, divusw            ,       & ! divergence
    tensionne, tensionnw, tensionse, tensionsw,       & ! tension
    shearne, shearnw, shearse, shearsw        ,       & ! shearing
    Deltane, Deltanw, Deltase, Deltasw        ,       & ! Delt
    zetax2ne, zetax2nw, zetax2se, zetax2sw    ,       & ! 2 x zeta (bulk visc)
    etax2ne, etax2nw, etax2se, etax2sw        ,       & ! 2 x eta (shear visc)
    rep_prsne, rep_prsnw, rep_prsse, rep_prssw,       & ! replacement pressure
    ssigpn, ssigps, ssigpe, ssigpw            ,       &
    ssigmn, ssigms, ssigme, ssigmw            ,       &
    ssig12n, ssig12s, ssig12e, ssig12w        ,       &
    ssigp1, ssigp2, ssigm1, ssigm2, ssig121, ssig122, &
    csigpne, csigpnw, csigpse, csigpsw        ,       &
    csigmne, csigmnw, csigmse, csigmsw        ,       &
    csig12ne, csig12nw, csig12se, csig12sw    ,       &
    str12ew, str12we, str12ns, str12sn        ,       &
    strp_tmp, strm_tmp

  real    (kind=dbl_kind) ::                                                   &
    tmp_uvel_ee, tmp_vvel_se, tmp_vvel_ee,                                     &
    tmp_vvel_ne, tmp_uvel_ne, tmp_uvel_se,                                     &
    tmp_uvel_cc, tmp_vvel_cc, tmp_dxT, tmp_dyT,                                &
    tmp_cxp, tmp_cyp, tmp_cxm, tmp_cym,                                        &
    tmp_strength, tmp_DminTarea, tmparea,                                      &
    tmp_dxhy, tmp_dyhx

  call domp_get_domain(lb,ub,il,iu)
  do iw = il, iu
    if (skipme(iw)) cycle
    ! divergence  =  e_11 + e_22
    tmp_uvel_cc = uvel(iw)
    tmp_vvel_cc = vvel(iw)
    tmp_uvel_ee = uvel(ee(iw))
    tmp_vvel_se = vvel(se(iw))
    tmp_vvel_ee = vvel(ee(iw))
    tmp_vvel_ne = vvel(ne(iw))
    tmp_uvel_ne = uvel(ne(iw))
    tmp_uvel_se = uvel(se(iw))
    tmp_dxT=dxT(iw)
    tmp_dyT=dyT(iw)
    tmp_cxp       = c1p5 * htn(iw) - p5 * htnm1(iw)
    tmp_cyp       = c1p5 * hte(iw) - p5 * htem1(iw)
    tmp_cxm       = -(c1p5 * htnm1(iw) - p5 * htn(iw))
    tmp_cym       = -(c1p5 * htem1(iw) - p5 * hte(iw))
    tmp_strength=strength(iw)
    tmparea = dxT(iw) * dyT(iw) ! necessary to split calc of DminTarea. Otherwize not binary identical
    tmp_DminTarea = deltaminEVP * tmparea
    tmp_dxhy      = p5 * (hte(iw) - htem1(iw))
    tmp_dyhx      = p5 * (htn(iw) - htnm1(iw))

    !--------------------------------------------------------------------------
    ! strain rates - NOTE these are actually strain rates * area  (m^2/s)
    !--------------------------------------------------------------------------
    call strain_rates (tmp_uvel_cc, tmp_vvel_cc, &
                            tmp_uvel_ee, tmp_vvel_ee, &
                            tmp_uvel_se, tmp_vvel_se, &
                            tmp_uvel_ne, tmp_vvel_ne, &
                            tmp_dxT,     tmp_dyT,     &
                            tmp_cxp,     tmp_cyp,     &
                            tmp_cxm,     tmp_cym,     &
                            divune,     divunw,       &
                            divuse,     divusw,       &
                            tensionne,  tensionnw,    &
                            tensionse,  tensionsw,    &
                            shearne,    shearnw,      &
                            shearse,    shearsw,      &
                            Deltane,    Deltanw,      &
                            Deltase,    Deltasw, e_factor )

    !--------------------------------------------------------------------------
    ! viscosities and replacement pressure
    !--------------------------------------------------------------------------
    call visc_replpress (tmp_strength, tmp_DminTarea, Deltane, &
                              zetax2ne, etax2ne, rep_prsne, capping, epp2i)

    call visc_replpress (tmp_strength, tmp_DminTarea, Deltanw, &
                              zetax2nw, etax2nw, rep_prsnw, capping, epp2i)

    call visc_replpress (tmp_strength, tmp_DminTarea, Deltasw, &
                              zetax2sw, etax2sw, rep_prssw, capping, epp2i)

    call visc_replpress (tmp_strength, tmp_DminTarea, Deltase, &
                              zetax2se, etax2se, rep_prsse, capping, epp2i)

    !--------------------------------------------------------------------------
    ! the stresses                            ! kg/s^2
    ! (1) northeast, (2) northwest, (3) southwest, (4) southeast
    !--------------------------------------------------------------------------
    ! NOTE: for comp. efficiency 2 x zeta and 2 x eta are used in the code
    stressp_1 (iw) = (stressp_1 (iw)*(c1-arlx1i*revp) &
                     + arlx1i*(zetax2ne*divune - rep_prsne)) * denom1
    stressp_2 (iw) = (stressp_2 (iw)*(c1-arlx1i*revp) &
                    + arlx1i*(zetax2nw*divunw - rep_prsnw)) * denom1
    stressp_3 (iw) = (stressp_3 (iw)*(c1-arlx1i*revp) &
                          + arlx1i*(zetax2sw*divusw - rep_prssw)) * denom1
    stressp_4 (iw) = (stressp_4 (iw)*(c1-arlx1i*revp) &
                          + arlx1i*(zetax2se*divuse - rep_prsse)) * denom1

    stressm_1 (iw) = (stressm_1 (iw)*(c1-arlx1i*revp) &
                          + arlx1i*etax2ne*tensionne) * denom1
    stressm_2 (iw) = (stressm_2 (iw)*(c1-arlx1i*revp) &
                          + arlx1i*etax2nw*tensionnw) * denom1
    stressm_3 (iw) = (stressm_3 (iw)*(c1-arlx1i*revp) &
                          + arlx1i*etax2sw*tensionsw) * denom1
    stressm_4 (iw) = (stressm_4 (iw)*(c1-arlx1i*revp) &
                           + arlx1i*etax2se*tensionse) * denom1

    stress12_1(iw) = (stress12_1(iw)*(c1-arlx1i*revp) &
                          + arlx1i*p5*etax2ne*shearne) * denom1
    stress12_2(iw) = (stress12_2(iw)*(c1-arlx1i*revp) &
                          + arlx1i*p5*etax2nw*shearnw) * denom1
    stress12_3(iw) = (stress12_3(iw)*(c1-arlx1i*revp) &
                          + arlx1i*p5*etax2sw*shearsw) * denom1
    stress12_4(iw) = (stress12_4(iw)*(c1-arlx1i*revp) &
                          + arlx1i*p5*etax2se*shearse) * denom1

     !--------------------------------------------------------------------------
     ! combinations of the stresses for the momentum equation ! kg/s^2
     !--------------------------------------------------------------------------
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

    csig12ne = p222*stress12_1(iw) + ssig122 &
             + p055*stress12_3(iw)
    csig12nw = p222*stress12_2(iw) + ssig121 &
             + p055*stress12_4(iw)
    csig12sw = p222*stress12_3(iw) + ssig122 &
             + p055*stress12_1(iw)
    csig12se = p222*stress12_4(iw) + ssig121 &
             + p055*stress12_2(iw)

    str12ew = p5*tmp_dxt*(p333*ssig12e + p166*ssig12w)
    str12we = p5*tmp_dxt*(p333*ssig12w + p166*ssig12e)
    str12ns = p5*tmp_dyt*(p333*ssig12n + p166*ssig12s)
    str12sn = p5*tmp_dyt*(p333*ssig12s + p166*ssig12n)

    !--------------------------------------------------------------------------
    ! for dF/dx (u momentum)
    !--------------------------------------------------------------------------
    strp_tmp  = p25*tmp_dyT*(p333*ssigpn  + p166*ssigps)
    strm_tmp  = p25*tmp_dyT*(p333*ssigmn  + p166*ssigms)

    ! northeast (i,j)
    str1(iw) = -strp_tmp - strm_tmp - str12ew &
               +tmp_dxhy*(-csigpne + csigmne) + tmp_dyhx*csig12ne

    ! northwest (i+1,j)
    str2(iw) = strp_tmp + strm_tmp - str12we &
               +tmp_dxhy*(-csigpnw + csigmnw) + tmp_dyhx*csig12nw

    strp_tmp  = p25*tmp_dyT*(p333*ssigps  + p166*ssigpn)
    strm_tmp  = p25*tmp_dyT*(p333*ssigms  + p166*ssigmn)

    ! southeast (i,j+1)
    str3(iw) = -strp_tmp - strm_tmp + str12ew &
               +tmp_dxhy*(-csigpse + csigmse) + tmp_dyhx*csig12se

    ! southwest (i+1,j+1)
    str4(iw) = strp_tmp + strm_tmp + str12we &
               +tmp_dxhy*(-csigpsw + csigmsw) + tmp_dyhx*csig12sw

    !--------------------------------------------------------------------------
    ! for dF/dy (v momentum)
    !--------------------------------------------------------------------------
    strp_tmp  = p25*tmp_dxT*(p333*ssigpe  + p166*ssigpw)
    strm_tmp  = p25*tmp_dxT*(p333*ssigme  + p166*ssigmw)

    ! northeast (i,j)
    str5(iw) = -strp_tmp + strm_tmp - str12ns &
               -tmp_dyhx*(csigpne + csigmne) + tmp_dxhy*csig12ne

    ! southeast (i,j+1)
    str6(iw) =  strp_tmp - strm_tmp - str12sn &
               -tmp_dyhx*(csigpse + csigmse) + tmp_dxhy*csig12se

    strp_tmp  = p25*tmp_dxT*(p333*ssigpw  + p166*ssigpe)
    strm_tmp  = p25*tmp_dxT*(p333*ssigmw  + p166*ssigme)

     ! northwest (i+1,j)
    str7(iw) = -strp_tmp + strm_tmp + str12ns &
               -tmp_dyhx*(csigpnw + csigmnw) + tmp_dxhy*csig12nw

    ! southwest (i+1,j+1)
    str8(iw) =  strp_tmp - strm_tmp + str12sn &
               -tmp_dyhx*(csigpsw + csigmsw) + tmp_dxhy*csig12sw
  enddo
end subroutine stress

!===============================================================================
! Compute strain rates
!
! author: Elizabeth C. Hunke, LANL
!
! 2019: subroutine created by Philippe Blain, ECCC
!DIR$ ATTRIBUTES FORCEINLINE :: strain_rates
pure subroutine strain_rates (tmp_uvel_cc, tmp_vvel_cc,                        &
                         tmp_uvel_ee, tmp_vvel_ee,                             &
                         tmp_uvel_se, tmp_vvel_se,                             &
                         tmp_uvel_ne, tmp_vvel_ne,                             &
                         dxT,        dyT,                                      &
                         cxp,        cyp,                                      &
                         cxm,        cym,                                      &
                         divune,     divunw,                                   &
                         divuse,     divusw,                                   &
                         tensionne,  tensionnw,                                &
                         tensionse,  tensionsw,                                &
                         shearne,    shearnw,                                  &
                         shearse,    shearsw,                                  &
                         Deltane,    Deltanw,                                  &
                         Deltase,    Deltasw, e_factor )

  use ice_kinds_mod

  real (kind=dbl_kind), intent(in) ::                                          &
         tmp_uvel_ee, tmp_vvel_ee, tmp_uvel_se, tmp_vvel_se,                   &
         tmp_uvel_cc, tmp_vvel_cc, tmp_uvel_ne, tmp_vvel_ne, e_factor

  real (kind=dbl_kind), intent(in) ::                                          &
         dxT      , & ! width of T-cell through the middle (m)
         dyT      , & ! height of T-cell through the middle (m)
         cyp      , & ! 1.5*HTE - 0.5*HTW
         cxp      , & ! 1.5*HTN - 0.5*HTS
         cym      , & ! 0.5*HTE - 1.5*HTW
         cxm          ! 0.5*HTN - 1.5*HTS

  real (kind=dbl_kind), intent(out)::               & ! at each corner :
        divune, divunw, divuse, divusw            , & ! divergence
        tensionne, tensionnw, tensionse, tensionsw, & ! tension
        shearne, shearnw, shearse, shearsw        , & ! shearing
        Deltane, Deltanw, Deltase, Deltasw            ! Delta

  !-----------------------------------------------------------------------------
  ! strain rates
  ! NOTE these are actually strain rates * area  (m^2/s)
  !-----------------------------------------------------------------------------

  ! divergence  =  e_11 + e_22
  divune    = cyp*tmp_uvel_cc - dyT*tmp_uvel_ee &
                + cxp*tmp_vvel_cc - dxT*tmp_vvel_se
  divunw    = cym*tmp_uvel_ee + dyT*tmp_uvel_cc &
                + cxp*tmp_vvel_ee - dxT*tmp_vvel_ne
  divusw    = cym*tmp_uvel_ne + dyT*tmp_uvel_se &
                + cxm*tmp_vvel_ne + dxT*tmp_vvel_ee
  divuse    = cyp*tmp_uvel_se - dyT*tmp_uvel_ne &
                + cxm*tmp_vvel_se + dxT*tmp_vvel_cc

  ! tension strain rate  =  e_11 - e_22
  tensionne = -cym*tmp_uvel_cc - dyT*tmp_uvel_ee &
                +  cxm*tmp_vvel_cc + dxT*tmp_vvel_se
  tensionnw = -cyp*tmp_uvel_ee + dyT*tmp_uvel_cc&
                +  cxm*tmp_vvel_ee + dxT*tmp_vvel_ne
  tensionsw = -cyp*tmp_uvel_ne + dyT*tmp_uvel_se &
                +  cxp*tmp_vvel_ne - dxT*tmp_vvel_ee
  tensionse = -cym*tmp_uvel_se - dyT*tmp_uvel_ne &
                +  cxp*tmp_vvel_se - dxT*tmp_vvel_cc

  ! shearing strain rate  =  2*e_12
  shearne = -cym*tmp_vvel_cc - dyT*tmp_vvel_ee &
              -  cxm*tmp_uvel_cc - dxT*tmp_uvel_se
  shearnw = -cyp*tmp_vvel_ee + dyT*tmp_vvel_cc &
              -  cxm*tmp_uvel_ee - dxT*tmp_uvel_ne
  shearsw = -cyp*tmp_vvel_ne + dyT*tmp_vvel_se &
              -  cxp*tmp_uvel_ne + dxT*tmp_uvel_ee
  shearse = -cym*tmp_vvel_se - dyT*tmp_vvel_ne &
              -  cxp*tmp_uvel_se + dxT*tmp_uvel_cc
  ! Delta (in the denominator of zeta, eta)
  Deltane = sqrt(divune**2 + e_factor*(tensionne**2 + shearne**2))
  Deltanw = sqrt(divunw**2 + e_factor*(tensionnw**2 + shearnw**2))
  Deltasw = sqrt(divusw**2 + e_factor*(tensionsw**2 + shearsw**2))
  Deltase = sqrt(divuse**2 + e_factor*(tensionse**2 + shearse**2))

end subroutine strain_rates

!===============================================================================
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

!DIR$ ATTRIBUTES FORCEINLINE :: visc_replpress
pure subroutine visc_replpress(strength, DminArea, Delta,                      &
                          zetax2, etax2, rep_prs, capping, epp2i)
  use ice_kinds_mod
  use ice_constants, only: c1
  use ice_dyn_shared, only: Ktens

  real (kind=dbl_kind), intent(in)::  strength, DminArea
  real (kind=dbl_kind), intent(in)::  Delta, capping, epp2i
  real (kind=dbl_kind), intent(out)::                                          &
         zetax2  , & ! bulk viscosity
         etax2   , & ! shear viscosity
         rep_prs     ! replacement pressure

  ! local variables
  real (kind=dbl_kind) :: tmpcalc     ! temporary

  ! NOTE: for comp. efficiency 2 x zeta and 2 x eta are used in the code

  tmpcalc = capping *(strength/max(Delta,DminArea))+ &
            (c1-capping)*(strength/(Delta + DminArea))
  zetax2  = (c1+Ktens)*tmpcalc
  rep_prs = (c1-Ktens)*tmpcalc*Delta
  etax2   = epp2i*zetax2

end subroutine visc_replpress

!===============================================================================
! Calculation of the surface stresses
! Integration of the momentum equation to find velocity (u,v)
! author: Elizabeth C. Hunke, LANL
!DIR$ ATTRIBUTES FORCEINLINE :: stepu
subroutine stepu (lb, ub,                                                      &
                        Cw,                                                    &
                        aiX,                                                   &
                        uocn,       vocn,                                      &
                        waterx,     watery,                                    &
                        forcex,     forcey,                                    &
                        umassdti,   fm,                                        &
                        uarear,                                                &
                        strintx,    strinty,                                   &
                        uvel_init,  vvel_init,                                 &
                        uvel,       vvel,                                      &
                        str1,str2,str3,str4,                                   &
                        str5,str6,str7,str8,                                   &
                        nw,sw,sse,skipme,                                      &
                        Tbu, Cb)
  use ice_kinds_mod
  use myomp, only : domp_get_domain
  use ice_constants, only: c0, c1
  use ice_dyn_shared, only: brlx, revp, u0, cosw, sinw
  implicit none
  ! arguments ------------------------------------------------------------------
  integer(kind=int_kind), intent(in)                           :: lb,ub
  integer(kind=int_kind), intent(in), dimension(:), contiguous :: nw,sw,sse
  logical(kind=log_kind), intent(in), dimension(:), contiguous :: skipme
  real   (kind=dbl_kind), intent(in), dimension(:), contiguous ::              &
     Tbu,      & ! coefficient for basal stress (N/m^2)
     uvel_init,& ! x-component of velocity (m/s), beginning of timestep
     vvel_init,& ! y-component of velocity (m/s), beginning of timestep
     aiX,      & ! ice fraction on u-grid
     waterx,   & ! for ocean stress calculation, x (m/s)
     watery,   & ! for ocean stress calculation, y (m/s)
     forcex,   & ! work array: combined atm stress and ocn tilt, x
     forcey,   & ! work array: combined atm stress and ocn tilt, y
     Umassdti, & ! mass of U-cell/dt (kg/m^2 s)
     uocn,     & ! ocean current, x-direction (m/s)
     vocn,     & ! ocean current, y-direction (m/s)
     fm,       & ! Coriolis param. * mass in U-cell (kg/s)
     uarear,   & ! 1/uarea
     Cw

  real (kind=dbl_kind),dimension(:), intent(in),    contiguous ::              &
     str1,str2,str3,str4,str5,str6,str7,str8
  real (kind=dbl_kind),dimension(:), intent(inout), contiguous ::              &
     uvel, vvel, strintx, strinty
  ! basal stress coefficient
  real (kind=dbl_kind),dimension(:), intent(out), contiguous :: Cb
  real (kind=dbl_kind), parameter ::                                           &
         rhow =  1026._dbl_kind               
  ! local variables
  integer (kind=int_kind) :: iw,il,iu
  real (kind=dbl_kind) ::                                                      &
         uold, vold        , & ! old-time uvel, vvel
         vrel              , & ! relative ice-ocean velocity
         cca,ccb,ab2,cc1,cc2,& ! intermediate variables
         taux, tauy          ! part of ocean stress term

  real (kind=dbl_kind) :: tmp_str2_nw,tmp_str3_sse,tmp_str4_sw,                &
                          tmp_str6_sse,tmp_str7_nw,tmp_str8_sw

  !-----------------------------------------------------------------------------
  ! integrate the momentum equation
  !-----------------------------------------------------------------------------
  call domp_get_domain(lb,ub,il,iu)
  do iw = il, iu
    if (skipme(iw)) cycle

    uold = uvel(iw)
    vold = vvel(iw)

    ! (magnitude of relative ocean current)*rhow*drag*aice
    vrel = aiX(iw)*rhow*Cw(iw)*sqrt((uocn(iw)-uold)**2+(vocn(iw)-vold)**2)
    ! ice/ocean stress
    taux = vrel*waterx(iw) ! NOTE this is not the entire
    tauy = vrel*watery(iw) ! ocn stress term

    Cb(iw)  = Tbu(iw) / (sqrt(uold**2 + vold**2) + u0) ! for basal stress

    ! revp = 0 for classic evp, 1 for revised evp
    cca = (brlx + revp)*umassdti(iw) + vrel * cosw + Cb(iw) ! kg/m^2 s
    ccb = fm(iw) + sign(c1,fm(iw)) * vrel * sinw ! kg/m^2 s

    ab2 = cca**2 + ccb**2

    tmp_str2_nw = str2(nw(iw))
    tmp_str3_sse = str3(sse(iw))
    tmp_str4_sw = str4(sw(iw))
    tmp_str6_sse = str6(sse(iw))
    tmp_str7_nw = str7(nw(iw))
    tmp_str8_sw = str8(sw(iw))

    ! divergence of the internal stress tensor
    strintx(iw) = uarear(iw)*(str1(iw)+tmp_str2_nw+tmp_str3_sse+tmp_str4_sw)

    strinty(iw) = uarear(iw)*(str5(iw)+tmp_str6_sse+tmp_str7_nw+tmp_str8_sw)

    ! finally, the velocity components
    cc1 = strintx(iw) + forcex(iw) + taux &
          + umassdti(iw)*(brlx*uold + revp*uvel_init(iw))
    cc2 = strinty(iw) + forcey(iw) + tauy &
          + umassdti(iw)*(brlx*vold + revp*vvel_init(iw))
    uvel(iw) = (cca*cc1 + ccb*cc2) / ab2 ! m/s
    vvel(iw) = (cca*cc2 - ccb*cc1) / ab2
    ! calculate seabed stress component for outputs
    ! only needed on last iteration.
  enddo

end subroutine stepu

end module bench_v2b
