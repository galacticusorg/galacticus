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
cosmologies of :cite:t:`montandon_decaying_2026`. This provides the mass-dependent critical overdensity
for collapse, :math:`\delta_\mathrm{c}(M_0)`, and the mapping between the initial Lagrangian mass, :math:`M_0`,
and the observed collapsed mass, :math:`M_\mathrm{coll}`, both of which are used to build the DDM halo mass
function.

*Note on unit conventions.* All physical computations in this module are performed in SI units
(the inputs---mass in :math:`\mathrm{M}_\odot`, times in Gyr, and velocity kick in km/s---are converted to
SI on entry). Only dimensionless combinations (:math:`\epsilon`, :math:`\xi`, :math:`\beta`, :math:`\tilde{\Gamma}`), ratios
(:math:`M_\mathrm{coll}/M_0`, :math:`\delta_\mathrm{c}` fractions), and the transition mass scale :math:`M_1` (returned
in :math:`\mathrm{M}_\odot`) are exposed.

*Note on the Einstein--de Sitter (EdS) internal relations.* Following
:cite:t:`montandon_decaying_2026`, the DDM corrections are evaluated using EdS relations *internally*
even when the host cosmology is :math:`\Lambda`CDM: the turnaround time is :math:`t_\mathrm{ta}=t_\mathrm{coll}/2`,
the turnaround radius follows from Kepler's relation :math:`G M_0 = \pi^2 R_\mathrm{ta}^3/(8 t_\mathrm{ta}^2)`,
the linear extrapolation uses the EdS growth factor, and :math:`\delta_\mathrm{c}` is normalized to
:math:`\delta_\mathrm{c}^\mathrm{EdS}=(3/5)(3\pi/2)^{2/3}\approx1.686`. This is a deliberate choice, not an
approximation of convenience: the fitting constants (:math:`A`, :math:`\beta`, :math:`\gamma`, :math:`\nu`, :math:`M_2/M_1`, and the
:math:`M_1` normalization) were calibrated against numerical solutions that fold *all* DDM physics
through the EdS cycloid. Substituting the true :math:`\Lambda`CDM turnaround/growth relations here would be
inconsistent with that calibration. The :math:`\Lambda`CDM baseline instead re-enters only where these
functions are consumed---the DDM critical overdensity is applied as a multiplicative correction,
:math:`\delta_\mathrm{c}^\mathrm{EdS}`-normalized, to a base :math:`\Lambda`CDM critical overdensity, so that the
:math:`\Lambda`CDM limit (:math:`v_k\rightarrow0` or :math:`\Gamma\rightarrow0`) is recovered exactly.

*Note on the calibration of* :math:`M_1`. The transition mass scale uses the calibrated fit of
:cite:t:`montandon_decaying_2026`, their eq. 44, and *not* the analytic estimate of their eq. 45---the
latter is an order-of-magnitude rationale for the :math:`v_k^3` scaling only, and over-estimates
:math:`M_1` by a factor of :math:`\sim200`. See the comment on ``fitMassScale1Coefficient`` below.
!!}

module Decaying_Dark_Matter_Spherical_Collapse
  use :: Numerical_Integration, only : integrator
  implicit none
  private
  public :: decayingDarkMatterEpsilon               , decayingDarkMatterGammaTilde , &
       &    decayingDarkMatterCriticalOverdensityEdS, decayingDarkMatterJIntegral  , &
       &    decayingDarkMatterDeltaCLarge           , decayingDarkMatterDeltaCSmall, &
       &    decayingDarkMatterMassScale1            , decayingDarkMatterDeltaCFit  , &
       &    decayingDarkMatterMassCollapsed

  ! Fitting constants from Montandon et al. (2026).
  ! Small-mass plateau [their eq. 42]: δ_c^small = δ_c^EdS + A Γ̃^β [ln(1+Γ̃)]^(1-γ).
  double precision, parameter :: fitSmallA             = 2.38240d0
  double precision, parameter :: fitSmallBeta          = 0.58180d0
  double precision, parameter :: fitSmallGamma         = 0.56420d0
  ! Transition function [their eq. 43]: exponent ν and the (universal) mass-scale ratio M₂/M₁.
  double precision, parameter :: fitTransitionNu       = 0.14840d0
  double precision, parameter :: fitTransitionMassRatio=10.0d0**1.3795d0 ! = M₂/M₁ = 10^1.3795
  ! Transition mass scale M₁ [their eq. 44]: M₁ = B vₖ³ Γ̃^(-1/2) t_ta, with log₁₀(B)=3.017 for
  ! M_2 in M☉, vₖ in km/s, and t_ta in Gyr. This is their *calibrated* fit (accurate to 10% across
  ! their grid of models and both redshifts), and is the expression used here.
  !
  ! Note that their eq. 45, M₁ ~ 2 √2 vₖ³ t_ta / π G, is only an order-of-magnitude analytic
  ! rationale for the vₖ³ scaling (it identifies the gravitating mass at turnaround with M₁ itself, via
  ! their "M_grav ~ M₁"), and is *not* the calibrated normalization: evaluated in these units its
  ! coefficient is 2.14e5, larger than B by a factor of ~206. Using it in place of B shifts the
  ! transition to far too high a mass and grossly over-suppresses the halo mass function. The value of B
  ! adopted here was verified by digitizing all fifteen δ_c(M_0) curves of their Fig. 5 (three
  ! lifetimes x five velocity kicks, at z=1.083) and fitting eq. 43 to each: B reproduces the resulting
  ! transition masses to ~13%, whereas the eq. 45 coefficient is high by a factor of ~230.
  double precision, parameter :: fitMassScale1Coefficient=10.0d0**3.017d0

  ! Tolerances for the adaptive quadratures. The integrands are smooth and bounded on their intervals, so
  ! modest relative tolerances suffice for the accuracy required of the derived quantities.
  double precision, parameter :: toleranceRelative       =1.0d-6, toleranceAbsolute=1.0d-10

  ! Parameters of the integrands, passed to the (argument-only) integrand functions via module-scope
  ! variables. These are threadprivate as the integrations may be in flight concurrently on multiple
  ! threads; the fBoundBar θ-integrand additionally nests the fBound u-integral, but that uses a disjoint
  ! set of variables so there is no interference.
  double precision :: jIntegrandGammaTilde
  double precision :: fBoundBeta               , fBoundXi
  double precision :: fBoundBarGammaTilde       , fBoundBarVelocityKickSI, fBoundBarTimeTurnaroundSI, &
       &              fBoundBarRadiusTurnaround
  !$omp threadprivate(jIntegrandGammaTilde,fBoundBeta,fBoundXi,fBoundBarGammaTilde,fBoundBarVelocityKickSI,fBoundBarTimeTurnaroundSI,fBoundBarRadiusTurnaround)

contains

  double precision function decayingDarkMatterEpsilon(velocityKick) result(epsilon)
    !!{RST
    Return the dimensionless mass-loss parameter :math:`\epsilon=(v_k/c)/(1+v_k/c)`
    (:cite:t:`montandon_decaying_2026`, their eq. 11), where :math:`v_k` is the velocity kick imparted to the
    daughter particle. ``velocityKick`` is in km/s.
    !!}
    use :: Numerical_Constants_Physical, only : speedLight
    use :: Numerical_Constants_Prefixes, only : kilo
    implicit none
    double precision, intent(in   ) :: velocityKick
    double precision                :: velocityKickFractional

    velocityKickFractional=+velocityKick             &
         &                 *kilo                     &
         &                 /speedLight
    epsilon               =+velocityKickFractional   &
         &                 /(                        &
         &                   +1.0d0                  &
         &                   +velocityKickFractional &
         &                  )
    return
  end function decayingDarkMatterEpsilon

  double precision function decayingDarkMatterGammaTilde(timeCollapse,lifetime) result(gammaTilde)
    !!{RST
    Return the dimensionless decay rate :math:`\tilde{\Gamma}=\Gamma t_\mathrm{coll}/2`
    (:cite:t:`montandon_decaying_2026`), where :math:`\Gamma` is the decay rate (the reciprocal of the
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
    :math:`\delta_\mathrm{c}^\mathrm{EdS}=(3/5)(3\pi/2)^{2/3}\approx1.686`
    (:cite:t:`montandon_decaying_2026`, their eq. 32).
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none

    deltaCEdS=+0.6d0            &
         &    *(                &
         &      +1.5d0          &
         &      *Pi             &
         &     )**(2.0d0/3.0d0)
    return
  end function decayingDarkMatterCriticalOverdensityEdS

  double precision function decayingDarkMatterJIntegral(gammaTilde) result(jIntegral)
    !!{RST
    Return the integral :math:`J(\tilde{\Gamma})` (:cite:t:`montandon_decaying_2026`, their eq. 38) that
    controls the large-mass limit of the critical overdensity,

    .. math::

       J(\tilde{\Gamma}) = -\int_0^{2\pi} \mathrm{d}\theta \, \frac{\sin\theta \left(1-\mathrm{e}^{-\tilde{\Gamma}\tilde{t}(\theta)}\right)}{(1-\cos\theta)^2} \left[6\pi+I(\theta)\right] ,

    with :math:`\tilde{t}(\theta)=(\theta-\sin\theta)/\pi` the EdS cycloid time and
    :math:`I(\theta)=\sin\theta-3\theta+4\tan(\theta/2)` (their eq. 36). The integrand is bounded on
    :math:`(0,2\pi)`---the apparent singularities at :math:`\theta=0`, :math:`\pi`, and :math:`2\pi` are removable---and :math:`J<0`,
    reflecting the delay of collapse relative to :math:`\Lambda`CDM.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: gammaTilde
    type            (integrator   ) :: integrator_

    ! Pass the decay rate to the integrand, and integrate. The integrand has removable singularities at
    ! θ=0, π, and 2π; the integral is split at θ=π so that all three fall on interval endpoints (which the
    ! Gauss-Kronrod rules never sample---unlike the interval midpoint θ=π on the full range, where the
    ! integrand would evaluate to a 0×∞ indeterminate form).
    jIntegrandGammaTilde=gammaTilde
    integrator_=integrator(decayingDarkMatterJIntegrand,toleranceRelative=toleranceRelative,toleranceAbsolute=toleranceAbsolute,hasSingularities=.true.)
    jIntegral  =-(                                          &
         &        +integrator_%integrate(0.0d0,       Pi )  &
         &        +integrator_%integrate(      Pi,2.0d0*Pi) &
         &       )
    return
  end function decayingDarkMatterJIntegral

  double precision function decayingDarkMatterJIntegrand(theta) result(integrand)
    !!{RST
    Integrand of :math:`J(\tilde{\Gamma})` (:cite:t:`montandon_decaying_2026`, their eq. 38); see
    ``decayingDarkMatterJIntegral``. The decay rate :math:`\tilde{\Gamma}` is passed via the module-scope
    ``jIntegrandGammaTilde``.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: theta
    double precision                :: sinTheta, oneMinusCosTheta, cycloidTime, iFunction

    sinTheta        =+sin(theta)
    oneMinusCosTheta=+1.0d0-cos(theta)
    cycloidTime     =+(theta-sinTheta)/Pi
    iFunction       =+sinTheta               &
         &           -3.0d0*theta            &
         &           +4.0d0*tan(0.5d0*theta)
    integrand       =+sinTheta                                       &
         &           *(1.0d0-exp(-jIntegrandGammaTilde*cycloidTime)) &
         &           /oneMinusCosTheta**2                            &
         &           *(6.0d0*Pi+iFunction)
    return
  end function decayingDarkMatterJIntegrand

  double precision function decayingDarkMatterDeltaCLarge(gammaTilde,epsilon,jIntegral) result(deltaCLarge)
    !!{RST
    Return the large-mass limit of the critical overdensity for collapse,
    :math:`\delta_\mathrm{c}^\mathrm{large}=\delta_\mathrm{c}^\mathrm{EdS}\left(1-\epsilon J(\tilde{\Gamma})/3\pi\right)`
    (:cite:t:`montandon_decaying_2026`, their eq. 40). ``jIntegral`` is :math:`J(\tilde{\Gamma})` as returned
    by ``decayingDarkMatterJIntegral``. In this limit all daughter particles are
    retained and :math:`\delta_\mathrm{c}` approaches (but slightly exceeds) the EdS value.
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
    :math:`\delta_\mathrm{c}^\mathrm{small}(\tilde{\Gamma})=\delta_\mathrm{c}^\mathrm{EdS}+A\,\tilde{\Gamma}^\beta\left[\ln(1+\tilde{\Gamma})\right]^{1-\gamma}`
    (:cite:t:`montandon_decaying_2026`, their eq. 42). In this limit all massive daughters escape, the
    collapse is driven purely by the decaying parent, and :math:`\delta_\mathrm{c}` depends only on
    :math:`\tilde{\Gamma}` (not on the velocity kick).
    !!}
    implicit none
    double precision, intent(in   ) :: gammaTilde

    deltaCSmall=+decayingDarkMatterCriticalOverdensityEdS()   &
         &      +fitSmallA                                    &
         &      *gammaTilde**fitSmallBeta                     &
         &      *log(1.0d0+gammaTilde)**(1.0d0-fitSmallGamma)
    return
  end function decayingDarkMatterDeltaCSmall

  double precision function decayingDarkMatterMassScale1(velocityKick,gammaTilde,timeCollapse) result(massScale1)
    !!{RST
    Return the transition mass scale :math:`M_1` (in :math:`\mathrm{M}_\odot`) between the small- and large-mass
    plateaux of the critical overdensity (:cite:t:`montandon_decaying_2026`, their eq. 44),

    .. math::

       M_1 = B\, v_k^3\, \tilde{\Gamma}^{-1/2}\, t_\mathrm{ta} ,

    with :math:`t_\mathrm{ta}=t_\mathrm{coll}/2` and :math:`\log_{10} B = 3.017` for :math:`M_1` in
    :math:`\mathrm{M}_\odot`, :math:`v_k` in km/s, and :math:`t_\mathrm{ta}` in Gyr. Physically this is
    (up to the calibrated normalization) the mass scale at which the kick velocity equals the halo
    orbital velocity at turnaround, with the residual :math:`\tilde{\Gamma}^{-1/2}` factor accounting for
    the mass lost by turnaround. As this is a fit in fixed units it requires no physical constants.
    ``velocityKick`` is in km/s and ``timeCollapse`` in Gyr.
    !!}
    implicit none
    double precision, intent(in   ) :: velocityKick, gammaTilde, timeCollapse
    double precision                :: timeTurnaround

    timeTurnaround=+0.5d0                    &
         &         *timeCollapse               ! Gyr.
    massScale1    =+fitMassScale1Coefficient &
         &         *velocityKick**3          &
         &         /sqrt(gammaTilde)         &
         &         *timeTurnaround
    return
  end function decayingDarkMatterMassScale1

  double precision function decayingDarkMatterDeltaCFit(mass0,deltaCLarge,deltaCSmall,massScale1) result(deltaCFit)
    !!{RST
    Return the mass-dependent critical overdensity for collapse, :math:`\delta_\mathrm{c}^\mathrm{fit}(M_0)`,
    interpolating between the small- and large-mass plateaux (:cite:t:`montandon_decaying_2026`, their
    eq. 43),

    .. math::

       \delta_\mathrm{c}^\mathrm{fit}(M_0) = \delta_\mathrm{c}^\mathrm{large} + \frac{\delta_\mathrm{c}^\mathrm{small}-\delta_\mathrm{c}^\mathrm{large}}{\left[\left(1+M_0/M_1\right)\left(1+\left(M_0/M_2\right)^4\right)\right]^\nu} ,

    with :math:`M_2=(M_2/M_1)M_1`. As :math:`M_0\rightarrow0` the denominator tends to unity and
    :math:`\delta_\mathrm{c}^\mathrm{fit}\rightarrow\delta_\mathrm{c}^\mathrm{small}`; as :math:`M_0\rightarrow\infty`
    the denominator diverges and :math:`\delta_\mathrm{c}^\mathrm{fit}\rightarrow\delta_\mathrm{c}^\mathrm{large}`.
    ``mass0`` and ``massScale1`` are in :math:`\mathrm{M}_\odot`.
    !!}
    implicit none
    double precision, intent(in   ) :: mass0      , deltaCLarge, &
         &                             deltaCSmall, massScale1
    double precision                :: massScale2 , denominator

    massScale2 =+fitTransitionMassRatio &
         &      *massScale1
    denominator=+(                            &
         &        +1.0d0                      &
         &        +mass0/massScale1           &
         &       )                            &
         &      *(                            &
         &        +1.0d0                      &
         &        +(mass0/massScale2)**4      &
         &       )
    deltaCFit  =+deltaCLarge                  &
         &      +(                            &
         &        +deltaCSmall                &
         &        -deltaCLarge                &
         &       )                            &
         &      /denominator**fitTransitionNu
    return
  end function decayingDarkMatterDeltaCFit

  double precision function decayingDarkMatterFBound(beta,xi) result(fBound)
    !!{RST
    Return the volume-averaged fraction of daughter particles that remain gravitationally bound to the
    halo, :math:`f_\mathrm{bound}(\beta,\xi)` (:cite:t:`montandon_decaying_2026`, their eq. 22), for the
    dimensionless bulk-flow parameter :math:`\beta` and kick parameter :math:`\xi`. This is evaluated as

    .. math::

       f_\mathrm{bound} = \int_0^1 3u^2\, P_\mathrm{bound}(u)\, \mathrm{d}u ,

    where :math:`P_\mathrm{bound}(u)=\max\left[0,\min\left(1,(1+C_\mathrm{bound}(u))/2\right)\right]` (their
    eqs. 18--19, written here as a single clamp) and
    :math:`C_\mathrm{bound}(u)=\left(3-\xi^2-u^2(1+\beta^2)\right)/(2\beta\xi u)` (their eq. 18). The clamped
    form reproduces all three population branches (fully bound, partially bound, unbound) without
    explicit case analysis.
    !!}
    implicit none
    double precision, intent(in   ) :: beta  , xi
    type            (integrator   ) :: integrator_

    ! Pass the bulk-flow and kick parameters to the integrand, and integrate over u∈(0,1). The integrand
    ! is bounded (and →0 at u=0), so plain adaptive Gauss-Kronrod quadrature suffices.
    fBoundBeta =beta
    fBoundXi   =xi
    integrator_=integrator(decayingDarkMatterFBoundIntegrand,toleranceRelative=toleranceRelative,toleranceAbsolute=toleranceAbsolute)
    fBound     =integrator_%integrate(0.0d0,1.0d0)
    return
  end function decayingDarkMatterFBound

  double precision function decayingDarkMatterFBoundIntegrand(u) result(integrand)
    !!{RST
    Integrand :math:`3u^2 P_\mathrm{bound}(u)` of :math:`f_\mathrm{bound}` (:cite:t:`montandon_decaying_2026`,
    their eq. 22); see ``decayingDarkMatterFBound``. The parameters :math:`\beta` and :math:`\xi` are
    passed via the module-scope ``fBoundBeta`` and ``fBoundXi``.
    !!}
    implicit none
    double precision, intent(in   ) :: u
    double precision                :: cBound, pBound

    cBound   =+(                     &
         &      +3.0d0               &
         &      -fBoundXi  **2       &
         &      -u         **2       &
         &      *(1.0d0+fBoundBeta**2) &
         &     )                     &
         &     /(                    &
         &       +2.0d0              &
         &       *fBoundBeta         &
         &       *fBoundXi           &
         &       *u                  &
         &      )
    ! P_bound as a clamp of (1+C_bound)/2 to [0,1]; IEEE +/-Inf (from β->0 at turnaround) clamps
    ! correctly to 1 or 0.
    pBound   =+max(0.0d0,min(1.0d0,0.5d0*(1.0d0+cBound)))
    integrand=+3.0d0*u**2*pBound
    return
  end function decayingDarkMatterFBoundIntegrand

  double precision function decayingDarkMatterFBoundBar(gammaTilde,velocityKick,timeCollapse,mass0) result(fBoundBar)
    !!{RST
    Return the decay-time-averaged bound fraction of daughter mass, :math:`\bar{f}_\mathrm{bound}`
    (:cite:t:`montandon_decaying_2026`, their eq. 48),

    .. math::

       \bar{f}_\mathrm{bound} = \frac{\tilde{\Gamma}/\pi}{1-\mathrm{e}^{-2\tilde{\Gamma}}} \int_0^{2\pi} \mathrm{d}\theta\, (1-\cos\theta)\, f_\mathrm{bound}^\mathrm{EdS}(\theta)\, \exp\left[-\tilde{\Gamma}(\theta-\sin\theta)/\pi\right] ,

    where :math:`f_\mathrm{bound}^\mathrm{EdS}(\theta)=f_\mathrm{bound}(\beta_\mathrm{EdS}(\theta),\xi_\mathrm{EdS}(\theta))`
    is evaluated along the EdS cycloid (their eq. 49) with
    :math:`\beta_\mathrm{EdS}(\theta)=\sqrt{2}\left|\cos(\theta/2)\right|` and
    :math:`\xi_\mathrm{EdS}(\theta)=(2 v_k t_\mathrm{ta}/\pi R_\mathrm{ta})\sqrt{1-\cos\theta}`. The turnaround
    radius :math:`R_\mathrm{ta}` follows from Kepler's relation for the initial Lagrangian mass ``mass0``
    (in :math:`\mathrm{M}_\odot`). This is the only mass-dependent ingredient of the mass mapping (through
    :math:`R_\mathrm{ta}\propto M_0^{1/3}`). ``velocityKick`` is in km/s and ``timeCollapse`` in Gyr.
    !!}
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Physical    , only : gravitationalConstant
    use :: Numerical_Constants_Astronomical, only : massSolar            , gigaYear
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    double precision, intent(in   ) :: gammaTilde    , velocityKick    , timeCollapse , mass0
    double precision                :: mass0SI       , normalization
    type            (integrator   ) :: integrator_

    ! Convert the (epoch- and mass-dependent) parameters to SI and store them for the integrand.
    fBoundBarGammaTilde      =+gammaTilde
    fBoundBarVelocityKickSI  =+velocityKick & ! km/s → m/s.
         &                    *kilo
    fBoundBarTimeTurnaroundSI=+0.5d0        & ! Gyr  → s.
         &                    *timeCollapse &
         &                    *gigaYear
    mass0SI                  =+mass0        & ! M☉   → kg.
         &                    *massSolar
    ! Turnaround radius from Kepler's relation, G M₀ = π² R_ta³/(8 t_ta²) [their eq. A2].
    fBoundBarRadiusTurnaround=+(                             &
         &                      +8.0d0                       &
         &                      *gravitationalConstant       &
         &                      *mass0SI                     &
         &                      *fBoundBarTimeTurnaroundSI**2 &
         &                      /Pi**2                       &
         &                     )**(1.0d0/3.0d0)
    ! Integrate over θ∈(0,2π), splitting at θ=π so that the endpoints θ=0, π, 2π---where β_EdS or ξ_EdS
    ! vanishes and the nested f_bound integrand is (removably) singular---are never sampled.
    integrator_=integrator(decayingDarkMatterFBoundBarIntegrand,toleranceRelative=toleranceRelative,toleranceAbsolute=toleranceAbsolute,hasSingularities=.true.)
    normalization=+gammaTilde                              &
         &        /Pi                                      &
         &        /(1.0d0-exp(-2.0d0*gammaTilde))
    fBoundBar    =+normalization                           &
         &        *(                                       &
         &          +integrator_%integrate(0.0d0,       Pi )  &
         &          +integrator_%integrate(      Pi,2.0d0*Pi) &
         &         )
    return
  end function decayingDarkMatterFBoundBar

  double precision function decayingDarkMatterFBoundBarIntegrand(theta) result(integrand)
    !!{RST
    Integrand :math:`(1-\cos\theta) f_\mathrm{bound}^\mathrm{EdS}(\theta) \exp[-\tilde{\Gamma}(\theta-\sin\theta)/\pi]`
    of :math:`\bar{f}_\mathrm{bound}` (:cite:t:`montandon_decaying_2026`, their eq. 48); see
    ``decayingDarkMatterFBoundBar``. The epoch- and mass-dependent parameters are passed via the
    module-scope ``fBoundBar*`` variables. This nests the ``fBound`` integral, which uses a disjoint set
    of module-scope variables.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: theta
    double precision                :: oneMinusCosTheta, betaEdS, xiEdS, weight

    oneMinusCosTheta=+1.0d0-cos(theta)
    betaEdS         =+sqrt(2.0d0)              &
         &           *abs(cos(0.5d0*theta))
    xiEdS           =+2.0d0                    &
         &           *fBoundBarVelocityKickSI  &
         &           *fBoundBarTimeTurnaroundSI &
         &           /Pi                       &
         &           /fBoundBarRadiusTurnaround &
         &           *sqrt(oneMinusCosTheta)
    ! Decay-time weight along the cycloid: (1-cos θ) exp[-Γ̃ (θ-sin θ)/π].
    weight          =+oneMinusCosTheta         &
         &           *exp(                     &
         &                -fBoundBarGammaTilde  &
         &                *(theta-sin(theta))   &
         &                /Pi                   &
         &               )
    integrand       =+weight                   &
         &           *decayingDarkMatterFBound(betaEdS,xiEdS)
    return
  end function decayingDarkMatterFBoundBarIntegrand

  double precision function decayingDarkMatterMassCollapsed(mass0,timeCollapse,lifetime,velocityKick) result(massCollapsed)
    !!{RST
    Return the collapsed (observed) halo mass, :math:`M_\mathrm{coll}` (in :math:`\mathrm{M}_\odot`), corresponding
    to an initial Lagrangian mass ``mass0`` (:cite:t:`montandon_decaying_2026`, their eq. 46),

    .. math::

       \frac{M_\mathrm{coll}}{M_0} = \mathrm{e}^{-\Gamma t_\mathrm{coll}} + \sqrt{1-2\epsilon}\left(1-\mathrm{e}^{-\Gamma t_\mathrm{coll}}\right)\bar{f}_\mathrm{bound} ,

    i.e. the surviving parent mass plus the retained rest mass of gravitationally bound daughters. The
    mapping is monotonic in :math:`M_0` (tending to :math:`\mathrm{e}^{-\Gamma t_\mathrm{coll}} M_0` at low mass,
    where all daughters escape, and to :math:`M_0` at high mass, where all are retained). ``timeCollapse`` and
    ``lifetime`` are in Gyr and ``velocityKick`` in km/s.
    !!}
    implicit none
    double precision, intent(in   ) :: mass0                  , timeCollapse, &
         &                             lifetime               , velocityKick
    double precision                :: epsilon                , gammaTilde  , &
         &                             survivingParentFraction, fBoundBar

    epsilon                = decayingDarkMatterEpsilon   (velocityKick         )
    gammaTilde             = decayingDarkMatterGammaTilde(timeCollapse,lifetime)
    ! Surviving (undecayed) parent fraction, e^{-Gamma t_coll} = e^{-2 gammaTilde}.
    survivingParentFraction=+exp(-2.0d0*gammaTilde)
    fBoundBar              = decayingDarkMatterFBoundBar(gammaTilde,velocityKick,timeCollapse,mass0)
    massCollapsed          =+mass0                             &
         &                  *(                                 &
         &                    +survivingParentFraction         &
         &                    +sqrt(1.0d0-2.0d0*epsilon)       &
         &                    *(1.0d0-survivingParentFraction) &
         &                    *fBoundBar                       &
         &                   )
    return
  end function decayingDarkMatterMassCollapsed

end module Decaying_Dark_Matter_Spherical_Collapse
