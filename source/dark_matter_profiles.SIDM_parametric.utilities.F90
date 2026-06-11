!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
!!    Andrew Benson <abenson@carnegiescience.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    Galacticus is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

!+    Contributions to this file made by: Niusha Ahvazi

!!{
Contains a module that implements the fitting functions of the parametric self-interacting dark matter model of
\cite{yang_parametric_2024}: the gravothermal evolution timescale and the maps between cold dark matter (NFW) and
self-interacting dark matter halo properties.
!!}

module SIDM_Parametric_Model
  !!{
  Implements the fitting functions of the parametric self-interacting dark matter model of \cite{yang_parametric_2024}.
  !!}
  implicit none
  private
  public :: get_tc   , dvmaxt , drmaxt   , Rmax_NFW, Vmax_NFW              , &
       &    r_s0     , rho_s0 , get_rho_s, get_r_s , get_r_c               , &
       &    radiusMaximumToScaleNFW      , densityScaleFactorNFW

  ! Constants related to NFW profiles.
  !! The ratio of the radius of the peak of the rotation curve to the scale radius.
  double precision, parameter :: radiusMaximumToScaleNFW=2.1625815870646120d0
  !! The factor, f, appearing in the expression:
  !    ρₛ = (Vₘₐₓ/f rₛ)^2 / G
  !  and which is equal to:
  !    f = √{4π(log(2)-½)} Vₛ/Vₘₐₓ
  double precision, parameter :: densityScaleFactorNFW  =1.6483500453640068d0

contains

  double precision function get_tc(darkMatterParticle_,C,Vmax,Rvmax,VmaxSIDM)
    !!{
    Evaluate the gravothermal evolution timescale $t_\mathrm{c}$ following eqn.~(2.2) of \cite{yang_parametric_2024}, given the
    self-interacting dark matter particle \mono{darkMatterParticle\_}, the dimensionless model parameter \mono{C}, the maximum
    circular velocity \mono{Vmax} and its radius \mono{Rvmax}, and the SIDM maximum circular velocity \mono{VmaxSIDM}.
    !!}
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal             , megaParsec, massSolar, MpcPerKmPerSToGyr
    use :: Numerical_Constants_Prefixes    , only : centi                                      , milli
    use :: Error                           , only : Error_Report
    use :: Dark_Matter_Particles           , only : darkMatterParticleClass                    , darkMatterParticleSelfInteractingDarkMatter
    implicit none
    class           (darkMatterParticleClass), intent(inout) :: darkMatterParticle_
    double precision                         , intent(in   ) :: C                             , Vmax    , &
         &                                                      Rvmax                         , VmaxSIDM
    ! Numerical coefficient in the gravothermal timescale (Yang et al. 2024; JCAP; 2; 32; eqn. 2.2).
    double precision                         , parameter     :: timescaleNormalization=150.0d0
    double precision                                         :: sigmaeff                      , reff    , &
         &                                                      rhoeff

    select type (darkMatterParticle_)
    class is (darkMatterParticleSelfInteractingDarkMatter)
       ! Get the effective cross section and convert to units of Mpc²/M☉.
       sigmaeff=+darkMatterParticle_%crossSectionEffective(VmaxSIDM) &
            &   *(centi     **2/milli    )                           &
            &   /(megaParsec**2/massSolar)
    class default
       sigmaeff=+0.0d0
       call Error_Report('unexpected class'//{introspection:location})
    end select
    reff   =+Rvmax                   &
         &  /radiusMaximumToScaleNFW
    rhoeff = rho_s0(reff,Vmax)
    get_tc =+timescaleNormalization                &
         &  /C                                     &
         &  /sigmaeff                              &
         &  /rhoeff                                &
         &  /reff                                  &
         &  /sqrt(                                 &
         &         +4.0d0                          &
         &         *Pi                             &
         &         *gravitationalConstant_internal &
         &         *rhoeff                         &
         &        )                                &
         &  *MpcPerKmPerSToGyr
    return
  end function get_tc

  double precision function dvmaxt(tau, Vmaxt)
    !!{
    Return the derivative $\mathrm{d}V_\mathrm{max}/\mathrm{d}\tau$ of the SIDM maximum circular
    velocity with respect to the dimensionless gravothermal time $\tau$, scaled by the CDM
    maximum circular velocity \mono{Vmaxt}, using the polynomial fit of
    \cite[][eqn.~2.4]{yang_parametric_2024}. The result is zero for $\tau>1$ (the core-collapsed
    regime).
    !!}
    implicit none
    double precision, intent(in   ) :: tau, Vmaxt

    if (tau > 1.0d0) then
       dvmaxt=+ 0.000000000000000d0
    else if (tau <= 1.0d0) then
       ! This is the derivative with respect to τ of equation 2.4 of Yang et al. 2024; JCAP; 2; 32. The exact values (more precise
       ! than those in the published paper) were taken from Daneng Yang's code.
       dvmaxt=+Vmaxt                         &
            & *(+ 0.177749020000000d0        &
            &   -13.195824689999998d0*tau**2 &
            &   +66.620926760000000d0*tau**3 &
            &   -94.337060499999990d0*tau**4 &
            &   +63.537661110000000d0*tau**6 &
            &   -21.925108889999997d0*tau**8 &
            & )
    end if
    return
  end function dvmaxt

  double precision function drmaxt(tau, Rmaxt)
    !!{
    Return the derivative $\mathrm{d}R_\mathrm{max}/\mathrm{d}\tau$ of the SIDM maximum-circular-velocity radius with respect to
    the dimensionless gravothermal time $\tau$, scaled by the CDM radius \mono{Rmaxt}, using the polynomial fit of
    \cite{yang_parametric_2024}. The result is zero for $\tau>1$ (the core-collapsed regime).
    !!}
    implicit none
    double precision, intent(in   ) :: tau, Rmaxt

    if (tau > 1.0d0) then
      drmaxt=+0.0d0
    else
       ! This is the derivative with respect to τ of equation 2.4 of Yang et al. 2024; JCAP; 2; 32. The exact values (more precise
       ! than those in the published paper) were taken from Daneng Yang's code.
       drmaxt=+Rmaxt                 &
            & *(                     &
            &   +0.00762288d0        &
            &   -1.43996392d0*tau    &
            &   +1.01282643d0*tau**2 &
            &   -0.55015288d0*tau**3 &
            &  )
    end if
    return
  end function drmaxt

  double precision function Rmax_NFW(RmaxSIDM, tau)
    !!{
    Map an SIDM maximum-circular-velocity radius \mono{RmaxSIDM} back to the equivalent CDM (NFW) value by
    inverting the $R_\mathrm{max}(\tau)$ evolution fit of \cite{yang_parametric_2024}. \mono{tau} is clamped to
    the range $[0,1]$.
    !!}
    implicit none
    double precision, intent(in) :: RmaxSIDM, tau
    double precision             :: tau_

    ! This is equation 2.4 of Yang et al. 2024; JCAP; 2; 32.
    tau_    =min(max(tau,0.0d0),1.0d0)
    Rmax_NFW=+RmaxSIDM             &
         &   /(                    &
         &     +1.000000d0         &
         &     +0.007623d0*tau_    &
         &     -0.720000d0*tau_**2 &
         &     +0.337600d0*tau_**3 &
         &     -0.137500d0*tau_**4 &
         &    )
    return
  end function Rmax_NFW

  double precision function Vmax_NFW(VmaxSIDM, tau)
    !!{
    Map an SIDM maximum circular velocity \mono{VmaxSIDM} back to the equivalent CDM (NFW) value by inverting the
    $V_\mathrm{max}(\tau)$ evolution fit of \cite{yang_parametric_2024}. \mono{tau} is clamped to the range
    $[0,1]$.
    !!}
    implicit none
    double precision, intent(in) :: VmaxSIDM, tau
    double precision             :: tau_

    ! This is equation 2.4 of Yang et al. 2024; JCAP; 2; 32.
    tau_    =min(max(tau,0.0d0),1.0d0)
    Vmax_NFW=+VmaxSIDM           &
         &   /(                  &
         &     +1.0000d0         &
         &     +0.1777d0*tau_    &
         &     -4.3990d0*tau_**3 &
         &     +16.660d0*tau_**4 &
         &     -18.870d0*tau_**5 &
         &     +9.0770d0*tau_**7 &
         &     -2.4360d0*tau_**9 &
         &    )
    return
  end function Vmax_NFW

  double precision function r_s0(Rmax)
    !!{
    Return the NFW scale radius corresponding to a maximum-circular-velocity radius \mono{Rmax}, using the
    standard NFW relation $r_\mathrm{s} \approx R_\mathrm{max}/2.163$.
    !!}
    implicit none
    double precision, intent(in) :: Rmax

    r_s0   =+Rmax                    &
         &  /radiusMaximumToScaleNFW
    return
  end function r_s0

  double precision function rho_s0(Rs, Vmax)
    !!{
    Return the NFW characteristic density corresponding to a given scale radius \mono{Rs} and maximum circular velocity
    \mono{Vmax}.
    !!}
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_astronomical, only : gravitationalConstant_internal
    implicit none
    double precision, intent(in) :: Rs, Vmax

    rho_s0=+(                              &
         &   +Vmax                         &
         &   /densityScaleFactorNFW        &
         &   /Rs                           &
         &  )**2                           &
         & /gravitationalConstant_internal
    return
  end function rho_s0

  double precision function get_rho_s(rho_s0, tau)
    !!{
    Return the SIDM parametric-profile characteristic density as a function of the dimensionless gravothermal time $\tau$,
    normalized to the initial NFW characteristic density \mono{rho\_s0}, using the fit of
    \cite{yang_parametric_2024}. \mono{tau} is clamped to the range $[0,1]$.
    !!}
    implicit none
    double precision, intent(in) :: rho_s0, tau
    double precision             :: tau_

    ! This is equation 2.3 of Yang et al. 2024; JCAP; 2; 32.
    tau_     =min(max(tau,0.0d0),1.0d0)
    get_rho_s=+rho_s0              &
         &    *(                   &
         &      +2.033d0           &
         &      +0.7381d0*tau_     &
         &      +7.2640d0*tau_**5  &
         &      -12.730d0*tau_**7  &
         &      +9.9150d0*tau_**9  &
         &      +(1.0d0-2.033d0)   &
         &      *log(tau_+0.001d0) &
         &      /log(    +0.001d0) &
         &     )
    return
  end function get_rho_s

  double precision function get_r_s(r_s0, tau)
    !!{
    Return the SIDM parametric-profile scale radius as a function of the dimensionless gravothermal time $\tau$, normalized to the
    initial NFW scale radius \mono{r\_s0}, using the fit of \cite{yang_parametric_2024}. \mono{tau} is clamped to the range $[0,1]$.
    !!}
    implicit none
    double precision, intent(in) :: r_s0, tau
    double precision             :: tau_

    ! This is equation 2.3 of Yang et al. 2024; JCAP; 2; 32.
    tau_   =min(max(tau,0.0d0),1.0d0)
    get_r_s=+r_s0               &
         & *(                   &
         &   +0.7178d0          &
         &   -0.1026d0*tau_     &
         &   +0.2474d0*tau_**2  &
         &   -0.4079d0*tau_**3  &
         &   +(1.0d0-0.7178d0)  &
         &   *log(tau_+0.001d0) &
         &   /log(    +0.001d0) &
         &  )
    return
  end function get_r_s

  double precision function get_r_c(r_s0, tau)
    !!{
    Return the SIDM parametric-profile core radius as a function of the dimensionless gravothermal time $\tau$, normalized to the
    initial NFW scale radius \mono{r\_s0}, using the fit of \cite{yang_parametric_2024}. \mono{tau} is clamped to the range $[0,1]$.
    !!}
    implicit none
    double precision, intent(in) :: r_s0, tau
    double precision             :: tau_

    ! This is equation 2.3 of Yang et al. 2024; JCAP; 2; 32.
    tau_   =min(max(tau,0.0d0),1.0d0)
    get_r_c=+r_s0                     &
         &  *(                        &
         &    +2.5550d0*sqrt(tau_)    &
         &    -3.6320d0*    tau_      &
         &    +2.1310d0*    tau_  **2 &
         &    -1.4150d0*    tau_  **3 &
         &    +0.4683d0*    tau_  **4 &
         &   )
    return
  end function get_r_c

end module SIDM_Parametric_Model
