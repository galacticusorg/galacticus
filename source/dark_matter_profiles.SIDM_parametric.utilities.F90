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
  public :: timescaleCollapse      , velocityMaximumRateTau, radiusMaximumRateTau, radiusMaximumNFW, velocityMaximumNFW, &
       &    radiusScaleNFW         , densityScaleNFW       , densityScale        , radiusScale     , radiusCore        , &
       &    radiusMaximumToScaleNFW                        , densityScaleFactorNFW

  ! Constants related to NFW profiles.
  !! The ratio of the radius of the peak of the rotation curve to the scale radius.
  double precision, parameter :: radiusMaximumToScaleNFW=2.1625815870646120d0
  !! The factor, f, appearing in the expression:
  !    ρₛ = (Vₘₐₓ/f rₛ)^2 / G
  !  and which is equal to:
  !    f = √{4π(log(2)-½)} Vₛ/Vₘₐₓ
  double precision, parameter :: densityScaleFactorNFW  =1.6483500453640068d0

contains

  double precision function timescaleCollapse(darkMatterParticle_,C,velocityMaximum,radiusVelocityMaximum,velocityMaximumSIDM)
    !!{
    Evaluate the gravothermal evolution timescale $t_\mathrm{c}$ following eqn.~(2.2) of \cite{yang_parametric_2024}, given the
    self-interacting dark matter particle \mono{darkMatterParticle\_}, the dimensionless model parameter \mono{C}, the maximum
    circular velocity \mono{velocityMaximum} and its radius \mono{radiusVelocityMaximum}, and the SIDM maximum circular velocity
    \mono{velocityMaximumSIDM}.
    !!}
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal, megaParsec                                 , massSolar, MpcPerKmPerSToGyr
    use :: Numerical_Constants_Prefixes    , only : centi                         , milli
    use :: Error                           , only : Error_Report
    use :: Dark_Matter_Particles           , only : darkMatterParticleClass       , darkMatterParticleSelfInteractingDarkMatter
    implicit none
    class           (darkMatterParticleClass), intent(inout) :: darkMatterParticle_
    double precision                         , intent(in   ) :: C                             , velocityMaximum    , &
         &                                                      radiusVelocityMaximum         , velocityMaximumSIDM
    ! Numerical coefficient in the gravothermal timescale (Yang et al. 2024; JCAP; 2; 32; eqn. 2.2).
    double precision                         , parameter     :: timescaleNormalization=150.0d0
    double precision                                         :: crossSectionEffective_        , radiusEffective    , &
         &                                                      densityEffective

    select type (darkMatterParticle_)
    class is (darkMatterParticleSelfInteractingDarkMatter)
       ! Get the effective cross section and convert to units of Mpc²/M☉.
       crossSectionEffective_=+darkMatterParticle_%crossSectionEffective(velocityMaximumSIDM) &
            &                 *(centi     **2/milli    )                                      &
            &                 /(megaParsec**2/massSolar)
    class default
       crossSectionEffective_=+0.0d0
       call Error_Report('unexpected class'//{introspection:location})
    end select
    radiusEffective  =+radiusVelocityMaximum   &
         &            /radiusMaximumToScaleNFW
    densityEffective = densityScaleNFW(radiusEffective,velocityMaximum)
    timescaleCollapse=+timescaleNormalization                &
         &            /C                                     &
         &            /crossSectionEffective_                &
         &            /densityEffective                      &
         &            /radiusEffective                       &
         &            /sqrt(                                 &
         &                   +4.0d0                          &
         &                   *Pi                             &
         &                   *gravitationalConstant_internal &
         &                   *densityEffective               &
         &                  )                                &
         &            *MpcPerKmPerSToGyr
    return
  end function timescaleCollapse

  double precision function velocityMaximumRateTau(tau,velocityMaximum)
    !!{
    Return the derivative $\mathrm{d}V_\mathrm{max}/\mathrm{d}\tau$ of the SIDM maximum circular
    velocity with respect to the dimensionless gravothermal time $\tau$, scaled by the CDM
    maximum circular velocity \mono{velocityMaximum}, using the polynomial fit of
    \cite[][eqn.~2.4]{yang_parametric_2024}. The result is zero for $\tau>1$ (the core-collapsed
    regime).
    !!}
    implicit none
    double precision, intent(in   ) :: tau, velocityMaximum

    if (tau > 1.0d0) then
       velocityMaximumRateTau=+   0.000000000000000d0
    else if (tau <= 1.0d0) then
       ! This is the derivative with respect to τ of equation 2.4 of Yang et al. 2024; JCAP; 2; 32. The exact values (more precise
       ! than those in the published paper) were taken from Daneng Yang's code.
       velocityMaximumRateTau=+velocityMaximum               &
            &                 *(+ 0.177749020000000d0        &
            &                   -13.195824689999998d0*tau**2 &
            &                   +66.620926760000000d0*tau**3 &
            &                   -94.337060499999990d0*tau**4 &
            &                   +63.537661110000000d0*tau**6 &
            &                   -21.925108889999997d0*tau**8 &
            &                 )
    end if
    return
  end function velocityMaximumRateTau

  double precision function radiusMaximumRateTau(tau,radiusMaximum)
    !!{
    Return the derivative $\mathrm{d}R_\mathrm{max}/\mathrm{d}\tau$ of the SIDM maximum-circular-velocity radius with respect to
    the dimensionless gravothermal time $\tau$, scaled by the CDM radius \mono{radiusMaximum}, using the polynomial fit of
    \cite{yang_parametric_2024}. The result is zero for $\tau>1$ (the core-collapsed regime).
    !!}
    implicit none
    double precision, intent(in   ) :: tau, radiusMaximum

    if (tau > 1.0d0) then
      radiusMaximumRateTau=+0.0d0
    else
       ! This is the derivative with respect to τ of equation 2.4 of Yang et al. 2024; JCAP; 2; 32. The exact values (more precise
       ! than those in the published paper) were taken from Daneng Yang's code.
       radiusMaximumRateTau=+radiusMaximum         &
            &               *(                     &
            &                 +0.00762288d0        &
            &                 -1.43996392d0*tau    &
            &                 +1.01282643d0*tau**2 &
            &                 -0.55015288d0*tau**3 &
            &                )
    end if
    return
  end function radiusMaximumRateTau

  double precision function radiusMaximumNFW(radiusMaximumSIDM,tau)
    !!{
    Map an SIDM maximum-circular-velocity radius \mono{radiusMaximumSIDM} back to the equivalent CDM (NFW) value by
    inverting the $R_\mathrm{max}(\tau)$ evolution fit of \cite{yang_parametric_2024}. \mono{tau} is clamped to
    the range $[0,1]$.
    !!}
    implicit none
    double precision, intent(in) :: radiusMaximumSIDM, tau
    double precision             :: tau_

    ! This is equation 2.4 of Yang et al. 2024; JCAP; 2; 32.
    tau_            =min(max(tau,0.0d0),1.0d0)
    radiusMaximumNFW=+radiusMaximumSIDM    &
         &           /(                    &
         &             +1.000000d0         &
         &             +0.007623d0*tau_    &
         &             -0.720000d0*tau_**2 &
         &             +0.337600d0*tau_**3 &
         &             -0.137500d0*tau_**4 &
         &            )
    return
  end function radiusMaximumNFW

  double precision function velocityMaximumNFW(velocityMaximumSIDM, tau)
    !!{
    Map an SIDM maximum circular velocity \mono{velocityMaximumSIDM} back to the equivalent CDM (NFW) value by inverting the
    $V_\mathrm{max}(\tau)$ evolution fit of \cite{yang_parametric_2024}. \mono{tau} is clamped to the range
    $[0,1]$.
    !!}
    implicit none
    double precision, intent(in) :: velocityMaximumSIDM, tau
    double precision             :: tau_

    ! This is equation 2.4 of Yang et al. 2024; JCAP; 2; 32.
    tau_              =min(max(tau,0.0d0),1.0d0)
    velocityMaximumNFW=+velocityMaximumSIDM &
         &             /(                   &
         &               +1.0000d0          &
         &               +0.1777d0*tau_     &
         &               -4.3990d0*tau_**3  &
         &               +16.660d0*tau_**4  &
         &               -18.870d0*tau_**5  &
         &               +9.0770d0*tau_**7  &
         &               -2.4360d0*tau_**9  &
         &              )
    return
  end function velocityMaximumNFW

  double precision function radiusScaleNFW(radiusMaximum)
    !!{
    Return the NFW scale radius corresponding to a maximum-circular-velocity radius \mono{radiusMaximum}, using the
    standard NFW relation $r_\mathrm{s} \approx R_\mathrm{max}/2.163$.
    !!}
    implicit none
    double precision, intent(in) :: radiusMaximum

    radiusScaleNFW=+radiusMaximum           &
         &         /radiusMaximumToScaleNFW
    return
  end function radiusScaleNFW

  double precision function densityScaleNFW(radiusScale_,velocityMaximum)
    !!{
    Return the NFW characteristic density corresponding to a given scale radius \mono{radiusScale\_} and maximum circular velocity
    \mono{velocityMaximum}.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    double precision, intent(in) :: radiusScale_, velocityMaximum

    densityScaleNFW=+(                              &
         &            +velocityMaximum              &
         &            /densityScaleFactorNFW        &
         &            /radiusScale_                 &
         &           )**2                           &
         &          /gravitationalConstant_internal
    return
  end function densityScaleNFW

  double precision function densityScale(densityScaleInitial,tau)
    !!{
    Return the SIDM parametric-profile characteristic density as a function of the dimensionless gravothermal time $\tau$,
    normalized to the initial NFW characteristic density \mono{densityScaleInitial}, using the fit of
    \cite{yang_parametric_2024}. \mono{tau} is clamped to the range $[0,1]$.
    !!}
    implicit none
    double precision, intent(in) :: densityScaleInitial, tau
    double precision             :: tau_

    ! This is equation 2.3 of Yang et al. 2024; JCAP; 2; 32.
    tau_        =min(max(tau,0.0d0),1.0d0)
    densityScale=+densityScaleInitial &
         &       *(                   &
         &         +2.033d0           &
         &         +0.7381d0*tau_     &
         &         +7.2640d0*tau_**5  &
         &         -12.730d0*tau_**7  &
         &         +9.9150d0*tau_**9  &
         &         +(1.0d0-2.033d0)   &
         &         *log(tau_+0.001d0) &
         &         /log(    +0.001d0) &
         &        )
    return
  end function densityScale

  double precision function radiusScale(radiusScaleInitial,tau)
    !!{
    Return the SIDM parametric-profile scale radius as a function of the dimensionless gravothermal time $\tau$, normalized to the
    initial NFW scale radius \mono{radiusScaleInitial}, using the fit of \cite{yang_parametric_2024}. \mono{tau} is clamped to the
    range $[0,1]$.
    !!}
    implicit none
    double precision, intent(in) :: radiusScaleInitial, tau
    double precision             :: tau_

    ! This is equation 2.3 of Yang et al. 2024; JCAP; 2; 32.
    tau_       =min(max(tau,0.0d0),1.0d0)
    radiusScale=+radiusScaleInitial  &
         &      *(                   &
         &        +0.7178d0          &
         &        -0.1026d0*tau_     &
         &        +0.2474d0*tau_**2  &
         &        -0.4079d0*tau_**3  &
         &        +(1.0d0-0.7178d0)  &
         &        *log(tau_+0.001d0) &
         &        /log(    +0.001d0) &
         &       )
    return
  end function radiusScale

  double precision function radiusCore(radiusScaleInitial,tau)
    !!{
    Return the SIDM parametric-profile core radius as a function of the dimensionless gravothermal time $\tau$, normalized to the
    initial NFW scale radius \mono{radiusScaleInitial}, using the fit of \cite{yang_parametric_2024}. \mono{tau} is clamped to the
    range $[0,1]$.
    !!}
    implicit none
    double precision, intent(in) :: radiusScaleInitial, tau
    double precision             :: tau_

    ! This is equation 2.3 of Yang et al. 2024; JCAP; 2; 32.
    tau_      =min(max(tau,0.0d0),1.0d0)
    radiusCore=+radiusScaleInitial       &
         &     *(                        &
         &       +2.5550d0*sqrt(tau_)    &
         &       -3.6320d0*    tau_      &
         &       +2.1310d0*    tau_  **2 &
         &       -1.4150d0*    tau_  **3 &
         &       +0.4683d0*    tau_  **4 &
         &      )
    return
  end function radiusCore

end module SIDM_Parametric_Model
