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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implementation of a gaseous dynamical friction class which uses the fitting formula of :cite:t:`kim_dynamical_2007` for a
  perturber on a circular orbit.
  !!}

  !![
  <dynamicalFrictionGaseous name="dynamicalFrictionGaseousKimKim2007" docformat="rst">
   <description>
   A gaseous dynamical friction class which uses the fitting formula of :cite:t:`kim_dynamical_2007` (their eqn. 15) for the
   azimuthal drag force on a perturber moving on a circular orbit of radius :math:`R_\mathrm{p}` through a uniform gaseous
   medium:

   .. math::

      \mathcal{I} = \left\{ \begin{array}{ll}
       0.7706 \ln\left(\frac{1+\mathcal{M}}{1.0004-0.9185\mathcal{M}}\right) - 1.4703\mathcal{M} &amp; \hbox{for } \mathcal{M} &lt; 1, \\
       \ln\left[ 330 (R_\mathrm{p}/r_\mathrm{min}) (\mathcal{M}-0.71)^{5.72} \mathcal{M}^{-9.58} \right] &amp; \hbox{for } 1 \le \mathcal{M} &lt; 4.4, \\
       \ln\left[ (R_\mathrm{p}/r_\mathrm{min}) / (0.11\mathcal{M}+1.65) \right] &amp; \hbox{for } 4.4 \le \mathcal{M},
      \end{array} \right.

   where :math:`\ln(R_\mathrm{p}/r_\mathrm{min})` is identified with the Coulomb logarithm, :math:`\ln\Lambda`. The drag force
   is applied antiparallel to the velocity of the perturber. The radial component of the force found by
   :cite:t:`kim_dynamical_2007`, which contributes little to orbital decay, is not included.

   The fits are discontinuous at :math:`\mathcal{M}=1` (except for :math:`R_\mathrm{p}/r_\mathrm{min} \approx 10`, the value
   to which they were calibrated) and at :math:`\mathcal{M}=4.4`. To avoid presenting a discontinuous force to the orbit
   integrator, :math:`\mathcal{I}` is interpolated linearly across the intervals :math:`[1-\delta,1+\delta]` and
   :math:`[4.4-\delta,4.4+\delta]`, where :math:`\delta=`\ ``[widthTransition]``. Setting ``[widthTransition]``\ :math:`=0`
   recovers the fits exactly. Since the fits can become negative for :math:`\mathcal{M} \ll 1` (where they are
   :math:`\lesssim 3\times 10^{-4}` in magnitude) and for very large :math:`\mathcal{M}`, :math:`\mathcal{I}` is limited to be
   non-negative.
   </description>
  </dynamicalFrictionGaseous>
  !!]
  type, extends(dynamicalFrictionGaseousClass) :: dynamicalFrictionGaseousKimKim2007
     !!{RST
     Implementation of a gaseous dynamical friction class which uses the fitting formula of :cite:t:`kim_dynamical_2007`.
     !!}
     private
     double precision :: widthTransition
   contains
     !![
     <methods docformat="rst">
       <method method="forceFactorSubsonic"     description="Return the fitting formula for the dimensionless drag force for subsonic perturbers."                                 />
       <method method="forceFactorTransonic"    description="Return the fitting formula for the dimensionless drag force for perturbers with :math:`1 \le \mathcal{M} &lt; 4.4`."/>
       <method method="forceFactorSupersonic"   description="Return the fitting formula for the dimensionless drag force for perturbers with :math:`\mathcal{M} \ge 4.4`."      />
     </methods>
     !!]
     procedure :: forceFactor           => kimKim2007ForceFactor
     procedure :: forceFactorSubsonic   => kimKim2007ForceFactorSubsonic
     procedure :: forceFactorTransonic  => kimKim2007ForceFactorTransonic
     procedure :: forceFactorSupersonic => kimKim2007ForceFactorSupersonic
  end type dynamicalFrictionGaseousKimKim2007

  interface dynamicalFrictionGaseousKimKim2007
     !!{RST
     Constructors for the :galacticus-class:`dynamicalFrictionGaseousKimKim2007` gaseous dynamical friction class.
     !!}
     module procedure kimKim2007ConstructorParameters
     module procedure kimKim2007ConstructorInternal
  end interface dynamicalFrictionGaseousKimKim2007

  ! Mach numbers at which the fitting formulae switch.
  double precision, parameter :: kimKim2007MachNumberTransonic=1.0d0, kimKim2007MachNumberSupersonic=4.4d0

contains

  function kimKim2007ConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`dynamicalFrictionGaseousKimKim2007` gaseous dynamical friction class which builds
    the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (dynamicalFrictionGaseousKimKim2007)                :: self
    type            (inputParameters                   ), intent(inout) :: parameters
    double precision                                                    :: widthTransition

    !![
    <inputParameter docformat="rst">
      <name>widthTransition</name>
      <defaultValue>0.05d0</defaultValue>
      <description>
      The half-width, :math:`\delta`, in Mach number of the intervals around :math:`\mathcal{M}=1` and :math:`\mathcal{M}=4.4`
      across which the :cite:t:`kim_dynamical_2007` fitting formulae are interpolated to make the drag force continuous.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=dynamicalFrictionGaseousKimKim2007(widthTransition)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function kimKim2007ConstructorParameters

  function kimKim2007ConstructorInternal(widthTransition) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`dynamicalFrictionGaseousKimKim2007` gaseous dynamical friction class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (dynamicalFrictionGaseousKimKim2007)                :: self
    double precision                                    , intent(in   ) :: widthTransition
    !![
    <constructorAssign variables="widthTransition"/>
    !!]

    ! The subsonic fit is valid only for M < 1.0892 (where its denominator vanishes), and the transonic fit only for M > 0.71,
    ! so the interpolation intervals must remain narrow.
    if (widthTransition < 0.0d0 .or. widthTransition >= 0.25d0) &
         & call Error_Report('[widthTransition] must be in the range [0,0.25)'//{introspection:location})
    return
  end function kimKim2007ConstructorInternal

  double precision function kimKim2007ForceFactor(self,machNumber,logarithmCoulomb) result(forceFactor)
    !!{RST
    Return the dimensionless drag force for the :cite:t:`kim_dynamical_2007` fitting formula.
    !!}
    implicit none
    class           (dynamicalFrictionGaseousKimKim2007), intent(inout) :: self
    double precision                                    , intent(in   ) :: machNumber           , logarithmCoulomb
    double precision                                                    :: machNumberLower      , machNumberUpper , &
         &                                                                 forceFactorLower     , forceFactorUpper

    if      (machNumber < kimKim2007MachNumberTransonic -self%widthTransition) then
       forceFactor=self%forceFactorSubsonic  (machNumber                 )
    else if (machNumber < kimKim2007MachNumberTransonic +self%widthTransition) then
       ! Interpolate across the transition from subsonic to transonic.
       machNumberLower =kimKim2007MachNumberTransonic-self%widthTransition
       machNumberUpper =kimKim2007MachNumberTransonic+self%widthTransition
       forceFactorLower=self%forceFactorSubsonic (machNumberLower                 )
       forceFactorUpper=self%forceFactorTransonic(machNumberUpper,logarithmCoulomb)
       forceFactor     =+forceFactorLower                         &
            &           +(forceFactorUpper-forceFactorLower)      &
            &           *(machNumber      -machNumberLower )      &
            &           /(machNumberUpper -machNumberLower )
    else if (machNumber < kimKim2007MachNumberSupersonic-self%widthTransition) then
       forceFactor=self%forceFactorTransonic (machNumber,logarithmCoulomb)
    else if (machNumber < kimKim2007MachNumberSupersonic+self%widthTransition) then
       ! Interpolate across the transition from transonic to supersonic.
       machNumberLower =kimKim2007MachNumberSupersonic-self%widthTransition
       machNumberUpper =kimKim2007MachNumberSupersonic+self%widthTransition
       forceFactorLower=self%forceFactorTransonic (machNumberLower,logarithmCoulomb)
       forceFactorUpper=self%forceFactorSupersonic(machNumberUpper,logarithmCoulomb)
       forceFactor     =+forceFactorLower                         &
            &           +(forceFactorUpper-forceFactorLower)      &
            &           *(machNumber      -machNumberLower )      &
            &           /(machNumberUpper -machNumberLower )
    else
       forceFactor=self%forceFactorSupersonic(machNumber,logarithmCoulomb)
    end if
    ! The fits can become (slightly) negative at very low and very high Mach numbers - limit to non-negative values.
    forceFactor=max(forceFactor,0.0d0)
    return
  end function kimKim2007ForceFactor

  double precision function kimKim2007ForceFactorSubsonic(self,machNumber) result(forceFactor)
    !!{RST
    Return the :cite:t:`kim_dynamical_2007` fitting formula for the dimensionless drag force for subsonic perturbers.
    !!}
    implicit none
    class           (dynamicalFrictionGaseousKimKim2007), intent(inout) :: self
    double precision                                    , intent(in   ) :: machNumber
    !$GLC attributes unused :: self

    forceFactor=+0.7706d0                             &
         &      *log(                                 &
         &           +(+1.0000d0+         machNumber) &
         &           /(+1.0004d0-0.9185d0*machNumber) &
         &          )                                 &
         &      -1.4703d0                             &
         &      *machNumber
    return
  end function kimKim2007ForceFactorSubsonic

  double precision function kimKim2007ForceFactorTransonic(self,machNumber,logarithmCoulomb) result(forceFactor)
    !!{RST
    Return the :cite:t:`kim_dynamical_2007` fitting formula for the dimensionless drag force for perturbers with
    :math:`1 \le \mathcal{M} < 4.4`.
    !!}
    implicit none
    class           (dynamicalFrictionGaseousKimKim2007), intent(inout) :: self
    double precision                                    , intent(in   ) :: machNumber, logarithmCoulomb
    !$GLC attributes unused :: self

    forceFactor=+logarithmCoulomb                   &
         &      +log(330.0d0)                       &
         &      +5.72d0*log(machNumber-0.71d0)      &
         &      -9.58d0*log(machNumber       )
    return
  end function kimKim2007ForceFactorTransonic

  double precision function kimKim2007ForceFactorSupersonic(self,machNumber,logarithmCoulomb) result(forceFactor)
    !!{RST
    Return the :cite:t:`kim_dynamical_2007` fitting formula for the dimensionless drag force for perturbers with
    :math:`\mathcal{M} \ge 4.4`.
    !!}
    implicit none
    class           (dynamicalFrictionGaseousKimKim2007), intent(inout) :: self
    double precision                                    , intent(in   ) :: machNumber, logarithmCoulomb
    !$GLC attributes unused :: self

    forceFactor=+logarithmCoulomb                   &
         &      -log(0.11d0*machNumber+1.65d0)
    return
  end function kimKim2007ForceFactorSupersonic
