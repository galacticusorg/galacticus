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

!!{RST
Contains a module which implements the revised spherical collapse model for decaying dark matter (DDM)
cosmologies of \cite{montandon_decaying_2026}. This provides the mass-dependent critical overdensity
for collapse, $\delta_\mathrm{c}(M_0)$, and the mapping between the initial Lagrangian mass, $M_0$,
and the observed collapsed mass, $M_\mathrm{coll}$, both of which are used to build the DDM halo mass
function.

\emph{Note on unit conventions.} All physical computations in this module are performed in SI units
(the inputs---mass in $\mathrm{M}_\odot$, times in Gyr, and velocity kick in km/s---are converted to
SI on entry). Only dimensionless combinations ($\epsilon$, $\xi$, $\beta$, $\tilde{\Gamma}$), ratios
($M_\mathrm{coll}/M_0$, $\delta_\mathrm{c}$ fractions), and the transition mass scale $M_1$ (returned
in $\mathrm{M}_\odot$) are exposed.

\emph{Note on the Einstein--de Sitter (EdS) internal relations.} Following
\cite{montandon_decaying_2026}, the DDM corrections are evaluated using EdS relations \emph{internally}
even when the host cosmology is $\Lambda$CDM: the turnaround time is $t_\mathrm{ta}=t_\mathrm{coll}/2$,
the turnaround radius follows from Kepler's relation $G M_0 = \pi^2 R_\mathrm{ta}^3/(8 t_\mathrm{ta}^2)$,
the linear extrapolation uses the EdS growth factor, and $\delta_\mathrm{c}$ is normalized to
$\delta_\mathrm{c}^\mathrm{EdS}=(3/5)(3\pi/2)^{2/3}\approx1.686$. This is a deliberate choice, not an
approximation of convenience: the fitting constants ($A$, $\beta$, $\gamma$, $\nu$, $M_2/M_1$, and the
$M_1$ normalization) were calibrated against numerical solutions that fold \emph{all} DDM physics
through the EdS cycloid. Substituting the true $\Lambda$CDM turnaround/growth relations here would be
inconsistent with that calibration. The $\Lambda$CDM baseline instead re-enters only where these
functions are consumed---the DDM critical overdensity is applied as a multiplicative correction,
$\delta_\mathrm{c}^\mathrm{EdS}$-normalized, to a base $\Lambda$CDM critical overdensity, so that the
$\Lambda$CDM limit ($v_k\rightarrow0$ or $\Gamma\rightarrow0$) is recovered exactly.
!!}

module Decaying_Dark_Matter_Spherical_Collapse
  implicit none
  private
  public :: decayingDarkMatterEpsilon              , decayingDarkMatterGammaTilde     , &
       &    decayingDarkMatterCriticalOverdensityEdS, decayingDarkMatterJIntegral      , &
       &    decayingDarkMatterDeltaCLarge          , decayingDarkMatterDeltaCSmall    , &
       &    decayingDarkMatterMassScale1           , decayingDarkMatterDeltaCFit      , &
       &    decayingDarkMatterMassCollapsed

  ! Fitting constants from Montandon et al. (2026).
  ! Small-mass plateau [their eq. 42]: delta_c^small = delta_c^EdS + A gammaTilde^beta [ln(1+gammaTilde)]^(1-gamma).
  double precision, parameter :: fitSmallA               =2.3824d0
  double precision, parameter :: fitSmallBeta            =0.5818d0
  double precision, parameter :: fitSmallGamma           =0.5642d0
  ! Transition function [their eq. 43]: exponent nu and the (universal) mass-scale ratio M2/M1.
  double precision, parameter :: fitTransitionNu         =0.1484d0
  double precision, parameter :: fitTransitionMassRatio  =23.96031d0 ! = M2/M1 = 10^1.3795 ~ 24 (literal, as real exponentiation is not a constant expression).
  ! Transition mass scale M1 [their eqs. 44,45]. We adopt the dimensionally-explicit physical form of
  ! their eq. 45, M1 = coefficient * v_k^3 t_ta / G * gammaTilde^(-1/2), with the coefficient 2*sqrt(2)/pi
  ! derived from the condition that the kick velocity equals the orbital velocity at turnaround, times
  ! the empirical gammaTilde^(-1/2) scaling. (Their eq. 44 gives an equivalent best-fit constant B, but
  ! in unstated units; the physical form here is unit-transparent and validated against their Fig. 5.)
  double precision, parameter :: fitMassScale1Coefficient=0.9003163161571062d0 ! = 2*sqrt(2)/pi (literal, as sqrt is not a constant expression).

  ! Numbers of abscissae for the fixed midpoint quadratures. The integrands are smooth and bounded on
  ! their (fixed) intervals, so a midpoint rule---which also avoids the interval endpoints where some
  ! factors are individually singular but the integrand has a finite limit---is accurate and robust.
  integer         , parameter :: quadraturePointsJ       =2000 ! theta-integral for J(gammaTilde)     [once per epoch].
  integer         , parameter :: quadraturePointsTheta   = 400 ! theta-integral for fBoundBar         [once per mass ].
  integer         , parameter :: quadraturePointsU       = 500 ! u-integral     for fBound(beta,xi).

contains

  double precision function decayingDarkMatterEpsilon(velocityKick) result(epsilon)
    !!{RST
    Return the dimensionless mass-loss parameter $\epsilon=(v_k/c)/(1+v_k/c)$
    (\citealt{montandon_decaying_2026}, their eq. 11), where $v_k$ is the velocity kick imparted to the
    daughter particle. ``velocityKick`` is in km/s.
    !!}
    use :: Numerical_Constants_Physical , only : speedLight
    use :: Numerical_Constants_Prefixes , only : kilo
    implicit none
    double precision, intent(in   ) :: velocityKick
    double precision                :: velocityKickFractional

    velocityKickFractional=+velocityKick     &
         &                 *kilo             &
         &                 /speedLight
    epsilon               =+velocityKickFractional  &
         &                 /(                        &
         &                   +1.0d0                  &
         &                   +velocityKickFractional &
         &                  )
    return
  end function decayingDarkMatterEpsilon

  double precision function decayingDarkMatterGammaTilde(timeCollapse,lifetime) result(gammaTilde)
    !!{RST
    Return the dimensionless decay rate $\tilde{\Gamma}=\Gamma t_\mathrm{coll}/2$
    (\citealt{montandon_decaying_2026}), where $\Gamma$ is the decay rate (the reciprocal of the
    ``lifetime``). Both ``timeCollapse`` and ``lifetime`` are in Gyr.
    !!}
    implicit none
    double precision, intent(in   ) :: timeCollapse, lifetime

    gammaTilde=+0.5d0        &
         &     *timeCollapse &
         &     /lifetime
    return
  end function decayingDarkMatterGammaTilde

  double precision function decayingDarkMatterCriticalOverdensityEdS() result(deltaCEdS)
    !!{RST
    Return the Einstein--de Sitter critical overdensity for collapse,
    $\delta_\mathrm{c}^\mathrm{EdS}=(3/5)(3\pi/2)^{2/3}\approx1.686$
    (\citealt{montandon_decaying_2026}, their eq. 32).
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none

    deltaCEdS=+0.6d0                &
         &    *(                    &
         &      +1.5d0             &
         &      *Pi                &
         &     )**(2.0d0/3.0d0)
    return
  end function decayingDarkMatterCriticalOverdensityEdS

  double precision function decayingDarkMatterJIntegral(gammaTilde) result(jIntegral)
    !!{RST
    Return the integral $J(\tilde{\Gamma})$ (\citealt{montandon_decaying_2026}, their eq. 38) that
    controls the large-mass limit of the critical overdensity,
    \begin{equation}
     J(\tilde{\Gamma}) = -\int_0^{2\pi} \mathrm{d}\theta \, \frac{\sin\theta \left(1-\mathrm{e}^{-\tilde{\Gamma}\tilde{t}(\theta)}\right)}{(1-\cos\theta)^2} \left[6\pi+I(\theta)\right] ,
    \end{equation}
    with $\tilde{t}(\theta)=(\theta-\sin\theta)/\pi$ the EdS cycloid time and
    $I(\theta)=\sin\theta-3\theta+4\tan(\theta/2)$ (their eq. 36). The integrand is bounded on
    $(0,2\pi)$---the apparent singularities at $\theta=0$, $\pi$, and $2\pi$ are removable---and $J<0$,
    reflecting the delay of collapse relative to $\Lambda$CDM.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: gammaTilde
    double precision                :: theta      , cycloidTime, sinTheta, oneMinusCosTheta, &
         &                             iFunction  , integrand  , summation
    integer                         :: i

    summation=0.0d0
    do i=1,quadraturePointsJ
       ! Midpoint abscissa in (0,2*pi), avoiding the endpoints and (for even count) theta=pi.
       theta           =+(dble(i)-0.5d0) &
            &            /dble(quadraturePointsJ) &
            &            *2.0d0                   &
            &            *Pi
       sinTheta        =+sin(theta)
       oneMinusCosTheta=+1.0d0-cos(theta)
       cycloidTime     =+(theta-sinTheta)/Pi
       iFunction       =+sinTheta           &
            &            -3.0d0*theta        &
            &            +4.0d0*tan(0.5d0*theta)
       integrand       =+sinTheta                                  &
            &            *(1.0d0-exp(-gammaTilde*cycloidTime))     &
            &            /oneMinusCosTheta**2                      &
            &            *(6.0d0*Pi+iFunction)
       summation       =+summation+integrand
    end do
    ! Multiply by the midpoint interval width (2*pi/N) and apply the overall minus sign.
    jIntegral=-summation             &
         &    *2.0d0*Pi              &
         &    /dble(quadraturePointsJ)
    return
  end function decayingDarkMatterJIntegral

  double precision function decayingDarkMatterDeltaCLarge(gammaTilde,epsilon,jIntegral) result(deltaCLarge)
    !!{RST
    Return the large-mass limit of the critical overdensity for collapse,
    $\delta_\mathrm{c}^\mathrm{large}=\delta_\mathrm{c}^\mathrm{EdS}\left(1-\epsilon J(\tilde{\Gamma})/3\pi\right)$
    (\citealt{montandon_decaying_2026}, their eq. 40). ``jIntegral`` is $J(\tilde{\Gamma})$ as returned
    by {\normalfont \ttfamily decayingDarkMatterJIntegral}. In this limit all daughter particles are
    retained and $\delta_\mathrm{c}$ approaches (but slightly exceeds) the EdS value.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: gammaTilde, epsilon, jIntegral
    !$GLC attributes unused :: gammaTilde

    deltaCLarge=+decayingDarkMatterCriticalOverdensityEdS() &
         &      *(                                          &
         &        +1.0d0                                    &
         &        -epsilon                                  &
         &        *jIntegral                                &
         &        /(3.0d0*Pi)                               &
         &       )
    return
  end function decayingDarkMatterDeltaCLarge

  double precision function decayingDarkMatterDeltaCSmall(gammaTilde) result(deltaCSmall)
    !!{RST
    Return the small-mass limit (plateau) of the critical overdensity for collapse,
    $\delta_\mathrm{c}^\mathrm{small}(\tilde{\Gamma})=\delta_\mathrm{c}^\mathrm{EdS}+A\,\tilde{\Gamma}^\beta\left[\ln(1+\tilde{\Gamma})\right]^{1-\gamma}$
    (\citealt{montandon_decaying_2026}, their eq. 42). In this limit all massive daughters escape, the
    collapse is driven purely by the decaying parent, and $\delta_\mathrm{c}$ depends only on
    $\tilde{\Gamma}$ (not on the velocity kick).
    !!}
    implicit none
    double precision, intent(in   ) :: gammaTilde

    deltaCSmall=+decayingDarkMatterCriticalOverdensityEdS()          &
         &      +fitSmallA                                           &
         &      *gammaTilde**fitSmallBeta                            &
         &      *log(1.0d0+gammaTilde)**(1.0d0-fitSmallGamma)
    return
  end function decayingDarkMatterDeltaCSmall

  double precision function decayingDarkMatterMassScale1(velocityKick,gammaTilde,timeCollapse) result(massScale1)
    !!{RST
    Return the transition mass scale $M_1$ (in $\mathrm{M}_\odot$) between the small- and large-mass
    plateaux of the critical overdensity (\citealt{montandon_decaying_2026}, their eqs. 44,45),
    \begin{equation}
     M_1 = \frac{2\sqrt{2}}{\pi G}\, v_k^3\, t_\mathrm{ta}\, \tilde{\Gamma}^{-1/2} ,
    \end{equation}
    with $t_\mathrm{ta}=t_\mathrm{coll}/2$. Physically this is the mass scale at which the kick velocity
    equals the halo orbital velocity at turnaround. ``velocityKick`` is in km/s and ``timeCollapse`` in
    Gyr.
    !!}
    use :: Numerical_Constants_Physical    , only : gravitationalConstant
    use :: Numerical_Constants_Astronomical, only : massSolar            , gigaYear
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    double precision, intent(in   ) :: velocityKick, gammaTilde, timeCollapse
    double precision                :: velocityKickSI, timeTurnaroundSI

    velocityKickSI  =+velocityKick      &
         &           *kilo             ! km/s -> m/s.
    timeTurnaroundSI=+0.5d0            &
         &           *timeCollapse     &
         &           *gigaYear         ! Gyr  -> s.
    massScale1      =+fitMassScale1Coefficient    &
         &           *velocityKickSI**3           &
         &           *timeTurnaroundSI            &
         &           /gravitationalConstant       &
         &           /sqrt(gammaTilde)            &
         &           /massSolar                 ! kg -> Msun.
    return
  end function decayingDarkMatterMassScale1

  double precision function decayingDarkMatterDeltaCFit(mass0,deltaCLarge,deltaCSmall,massScale1) result(deltaCFit)
    !!{RST
    Return the mass-dependent critical overdensity for collapse, $\delta_\mathrm{c}^\mathrm{fit}(M_0)$,
    interpolating between the small- and large-mass plateaux (\citealt{montandon_decaying_2026}, their
    eq. 43),
    \begin{equation}
     \delta_\mathrm{c}^\mathrm{fit}(M_0) = \delta_\mathrm{c}^\mathrm{large} + \frac{\delta_\mathrm{c}^\mathrm{small}-\delta_\mathrm{c}^\mathrm{large}}{\left[\left(1+M_0/M_1\right)\left(1+\left(M_0/M_2\right)^4\right)\right]^\nu} ,
    \end{equation}
    with $M_2=(M_2/M_1)M_1$. As $M_0\rightarrow0$ the denominator tends to unity and
    $\delta_\mathrm{c}^\mathrm{fit}\rightarrow\delta_\mathrm{c}^\mathrm{small}$; as $M_0\rightarrow\infty$
    the denominator diverges and $\delta_\mathrm{c}^\mathrm{fit}\rightarrow\delta_\mathrm{c}^\mathrm{large}$.
    ``mass0`` and ``massScale1`` are in $\mathrm{M}_\odot$.
    !!}
    implicit none
    double precision, intent(in   ) :: mass0, deltaCLarge, deltaCSmall, massScale1
    double precision                :: massScale2, denominator

    massScale2  =+fitTransitionMassRatio &
         &       *massScale1
    denominator =+(                                     &
         &         +1.0d0                               &
         &         +mass0/massScale1                    &
         &        )                                     &
         &       *(                                     &
         &         +1.0d0                               &
         &         +(mass0/massScale2)**4               &
         &        )
    deltaCFit   =+deltaCLarge                           &
         &       +(                                     &
         &         +deltaCSmall                         &
         &         -deltaCLarge                         &
         &        )                                     &
         &       /denominator**fitTransitionNu
    return
  end function decayingDarkMatterDeltaCFit

  double precision function decayingDarkMatterFBound(beta,xi) result(fBound)
    !!{RST
    Return the volume-averaged fraction of daughter particles that remain gravitationally bound to the
    halo, $f_\mathrm{bound}(\beta,\xi)$ (\citealt{montandon_decaying_2026}, their eq. 22), for the
    dimensionless bulk-flow parameter $\beta$ and kick parameter $\xi$. This is evaluated as
    \begin{equation}
     f_\mathrm{bound} = \int_0^1 3u^2\, P_\mathrm{bound}(u)\, \mathrm{d}u ,
    \end{equation}
    where $P_\mathrm{bound}(u)=\max\left[0,\min\left(1,(1+C_\mathrm{bound}(u))/2\right)\right]$ (their
    eqs. 18--19, written here as a single clamp) and
    $C_\mathrm{bound}(u)=\left(3-\xi^2-u^2(1+\beta^2)\right)/(2\beta\xi u)$ (their eq. 18). The clamped
    form reproduces all three population branches (fully bound, partially bound, unbound) without
    explicit case analysis.
    !!}
    implicit none
    double precision, intent(in   ) :: beta      , xi
    double precision                :: u         , cBound, pBound, summation
    integer                         :: i

    summation=0.0d0
    do i=1,quadraturePointsU
       ! Midpoint abscissa in (0,1), avoiding u=0 where C_bound is (removably) singular.
       u        =+(dble(i)-0.5d0)          &
            &     /dble(quadraturePointsU)
       cBound   =+(                        &
            &      +3.0d0                  &
            &      -xi **2                 &
            &      -u  **2                 &
            &      *(1.0d0+beta**2)        &
            &     )                        &
            &     /(                       &
            &       +2.0d0                 &
            &       *beta                  &
            &       *xi                    &
            &       *u                     &
            &      )
       ! P_bound as a clamp of (1+C_bound)/2 to [0,1]; IEEE +/-Inf (from beta->0 at turnaround) clamps
       ! correctly to 1 or 0.
       pBound   =+max(0.0d0,min(1.0d0,0.5d0*(1.0d0+cBound)))
       summation=+summation           &
            &     +3.0d0*u**2*pBound
    end do
    fBound=+summation                &
         & /dble(quadraturePointsU)
    return
  end function decayingDarkMatterFBound

  double precision function decayingDarkMatterFBoundBar(gammaTilde,velocityKick,timeCollapse,mass0) result(fBoundBar)
    !!{RST
    Return the decay-time-averaged bound fraction of daughter mass, $\bar{f}_\mathrm{bound}$
    (\citealt{montandon_decaying_2026}, their eq. 48),
    \begin{equation}
     \bar{f}_\mathrm{bound} = \frac{\tilde{\Gamma}/\pi}{1-\mathrm{e}^{-2\tilde{\Gamma}}} \int_0^{2\pi} \mathrm{d}\theta\, (1-\cos\theta)\, f_\mathrm{bound}^\mathrm{EdS}(\theta)\, \exp\left[-\tilde{\Gamma}(\theta-\sin\theta)/\pi\right] ,
    \end{equation}
    where $f_\mathrm{bound}^\mathrm{EdS}(\theta)=f_\mathrm{bound}(\beta_\mathrm{EdS}(\theta),\xi_\mathrm{EdS}(\theta))$
    is evaluated along the EdS cycloid (their eq. 49) with
    $\beta_\mathrm{EdS}(\theta)=\sqrt{2}\left|\cos(\theta/2)\right|$ and
    $\xi_\mathrm{EdS}(\theta)=(2 v_k t_\mathrm{ta}/\pi R_\mathrm{ta})\sqrt{1-\cos\theta}$. The turnaround
    radius $R_\mathrm{ta}$ follows from Kepler's relation for the initial Lagrangian mass ``mass0``
    (in $\mathrm{M}_\odot$). This is the only mass-dependent ingredient of the mass mapping (through
    $R_\mathrm{ta}\propto M_0^{1/3}$). ``velocityKick`` is in km/s and ``timeCollapse`` in Gyr.
    !!}
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Physical    , only : gravitationalConstant
    use :: Numerical_Constants_Astronomical, only : massSolar            , gigaYear
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    double precision, intent(in   ) :: gammaTilde     , velocityKick    , timeCollapse, mass0
    double precision                :: velocityKickSI , timeTurnaroundSI, mass0SI     , radiusTurnaround, &
         &                             theta          , oneMinusCosTheta, betaEdS     , xiEdS           , &
         &                             weight         , summation       , normalization
    integer                         :: i

    velocityKickSI  =+velocityKick        &
         &           *kilo               ! km/s -> m/s.
    timeTurnaroundSI=+0.5d0              &
         &           *timeCollapse       &
         &           *gigaYear           ! Gyr  -> s.
    mass0SI         =+mass0              &
         &           *massSolar          ! Msun -> kg.
    ! Turnaround radius from Kepler's relation, G M_0 = pi^2 R_ta^3/(8 t_ta^2) [their eq. A2].
    radiusTurnaround=+(                                    &
         &             +8.0d0                              &
         &             *gravitationalConstant              &
         &             *mass0SI                            &
         &             *timeTurnaroundSI**2                &
         &             /Pi**2                              &
         &            )**(1.0d0/3.0d0)
    summation=0.0d0
    do i=1,quadraturePointsTheta
       theta           =+(dble(i)-0.5d0)              &
            &            /dble(quadraturePointsTheta) &
            &            *2.0d0                        &
            &            *Pi
       oneMinusCosTheta=+1.0d0-cos(theta)
       betaEdS         =+sqrt(2.0d0)                  &
            &            *abs(cos(0.5d0*theta))
       xiEdS           =+2.0d0                        &
            &            *velocityKickSI              &
            &            *timeTurnaroundSI            &
            &            /Pi                          &
            &            /radiusTurnaround            &
            &            *sqrt(oneMinusCosTheta)
       ! Decay-time weight along the cycloid: (1-cos theta) exp[-gammaTilde (theta-sin theta)/pi].
       weight          =+oneMinusCosTheta             &
            &            *exp(                         &
            &                 -gammaTilde              &
            &                 *(theta-sin(theta))      &
            &                 /Pi                      &
            &                )
       summation       =+summation                    &
            &            +weight                       &
            &            *decayingDarkMatterFBound(betaEdS,xiEdS)
    end do
    ! Midpoint interval width (2*pi/N) times the normalized prefactor.
    normalization=+gammaTilde                          &
         &        /Pi                                  &
         &        /(1.0d0-exp(-2.0d0*gammaTilde))
    fBoundBar    =+normalization                       &
         &        *summation                           &
         &        *2.0d0*Pi                            &
         &        /dble(quadraturePointsTheta)
    return
  end function decayingDarkMatterFBoundBar

  double precision function decayingDarkMatterMassCollapsed(mass0,timeCollapse,lifetime,velocityKick) result(massCollapsed)
    !!{RST
    Return the collapsed (observed) halo mass, $M_\mathrm{coll}$ (in $\mathrm{M}_\odot$), corresponding
    to an initial Lagrangian mass ``mass0`` (\citealt{montandon_decaying_2026}, their eq. 46),
    \begin{equation}
     \frac{M_\mathrm{coll}}{M_0} = \mathrm{e}^{-\Gamma t_\mathrm{coll}} + \sqrt{1-2\epsilon}\left(1-\mathrm{e}^{-\Gamma t_\mathrm{coll}}\right)\bar{f}_\mathrm{bound} ,
    \end{equation}
    i.e. the surviving parent mass plus the retained rest mass of gravitationally bound daughters. The
    mapping is monotonic in $M_0$ (tending to $\mathrm{e}^{-\Gamma t_\mathrm{coll}} M_0$ at low mass,
    where all daughters escape, and to $M_0$ at high mass, where all are retained). ``timeCollapse`` and
    ``lifetime`` are in Gyr and ``velocityKick`` in km/s.
    !!}
    implicit none
    double precision, intent(in   ) :: mass0        , timeCollapse, lifetime, velocityKick
    double precision                :: epsilon      , gammaTilde  , survivingParentFraction, fBoundBar

    epsilon                =decayingDarkMatterEpsilon   (velocityKick             )
    gammaTilde             =decayingDarkMatterGammaTilde(timeCollapse,lifetime    )
    ! Surviving (undecayed) parent fraction, e^{-Gamma t_coll} = e^{-2 gammaTilde}.
    survivingParentFraction=+exp(-2.0d0*gammaTilde)
    fBoundBar              =decayingDarkMatterFBoundBar(gammaTilde,velocityKick,timeCollapse,mass0)
    massCollapsed          =+mass0                                      &
         &                  *(                                          &
         &                    +survivingParentFraction                 &
         &                    +sqrt(1.0d0-2.0d0*epsilon)               &
         &                    *(1.0d0-survivingParentFraction)         &
         &                    *fBoundBar                               &
         &                   )
    return
  end function decayingDarkMatterMassCollapsed

end module Decaying_Dark_Matter_Spherical_Collapse
