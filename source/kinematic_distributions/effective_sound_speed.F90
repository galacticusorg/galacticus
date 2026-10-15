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
  Implementation of a kinematic distribution class for collisional gas with a uniform effective sound speed.
  !!}

  !![
  <kinematicsDistribution name="kinematicsDistributionEffectiveSoundSpeed" docformat="rst">
   <description>
   A kinematic distribution class for collisional gas, such as the interstellar medium, supported by a combination of thermal
   pressure and turbulence. The gas is described by a single, uniform, effective isothermal sound speed,
   :math:`c_\mathrm{s}=`\ ``[soundSpeedEffective]``, which includes both thermal and turbulent contributions. The 1D velocity
   dispersion is equal to :math:`c_\mathrm{s}`, and the temperature is that of an isothermal ideal gas of primordial composition
   with this sound speed, :math:`T = \mu m_\mathrm{u} c_\mathrm{s}^2 / k_\mathrm{B}`.
   </description>
  </kinematicsDistribution>
  !!]
  type, public, extends(kinematicsDistributionClass) :: kinematicsDistributionEffectiveSoundSpeed
     !!{RST
     A kinematic distribution for collisional gas with a uniform effective sound speed.
     !!}
     double precision :: soundSpeedEffective
   contains
     procedure :: isCollisional        => effectiveSoundSpeedIsCollisional
     procedure :: soundSpeed           => effectiveSoundSpeedSoundSpeed
     procedure :: temperature          => effectiveSoundSpeedTemperature
     procedure :: velocityDispersion1D => effectiveSoundSpeedVelocityDispersion1D
  end type kinematicsDistributionEffectiveSoundSpeed

  interface kinematicsDistributionEffectiveSoundSpeed
     !!{RST
     Constructors for the :galacticus-class:`kinematicsDistributionEffectiveSoundSpeed` kinematic distribution class.
     !!}
     module procedure effectiveSoundSpeedConstructorParameters
     module procedure effectiveSoundSpeedConstructorInternal
  end interface kinematicsDistributionEffectiveSoundSpeed

contains

  function effectiveSoundSpeedConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`kinematicsDistributionEffectiveSoundSpeed` kinematic distribution class which builds
    the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (kinematicsDistributionEffectiveSoundSpeed)                :: self
    type            (inputParameters                          ), intent(inout) :: parameters
    double precision                                                           :: soundSpeedEffective

    !![
    <inputParameter docformat="rst">
      <name>soundSpeedEffective</name>
      <description>
      The effective isothermal sound speed (in units of km/s), including both thermal and turbulent contributions, of the gas.
      </description>
      <defaultValue>10.0d0</defaultValue>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=kinematicsDistributionEffectiveSoundSpeed(soundSpeedEffective)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function effectiveSoundSpeedConstructorParameters

  function effectiveSoundSpeedConstructorInternal(soundSpeedEffective) result(self)
    !!{RST
    Constructor for the :galacticus-class:`kinematicsDistributionEffectiveSoundSpeed` kinematic distribution class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (kinematicsDistributionEffectiveSoundSpeed)                :: self
    double precision                                           , intent(in   ) :: soundSpeedEffective
    !![
    <constructorAssign variables="soundSpeedEffective"/>
    !!]

    if (soundSpeedEffective <= 0.0d0) call Error_Report('[soundSpeedEffective] must be positive'//{introspection:location})
    return
  end function effectiveSoundSpeedConstructorInternal

  logical function effectiveSoundSpeedIsCollisional(self)
    !!{RST
    Return true indicating that the effective sound speed kinematic distribution represents collisional particles.
    !!}
    implicit none
    class(kinematicsDistributionEffectiveSoundSpeed), intent(inout) :: self
    !$GLC attributes unused :: self

    effectiveSoundSpeedIsCollisional=.true.
    return
  end function effectiveSoundSpeedIsCollisional

  double precision function effectiveSoundSpeedSoundSpeed(self,coordinates)
    !!{RST
    Return the sound speed at the specified ``coordinates`` in an effective sound speed kinematic distribution.
    !!}
    implicit none
    class(kinematicsDistributionEffectiveSoundSpeed), intent(inout) :: self
    class(coordinate                               ), intent(in   ) :: coordinates
    !$GLC attributes unused :: coordinates

    effectiveSoundSpeedSoundSpeed=self%soundSpeedEffective
    return
  end function effectiveSoundSpeedSoundSpeed

  double precision function effectiveSoundSpeedTemperature(self,coordinates)
    !!{RST
    Return the temperature at the specified ``coordinates`` in an effective sound speed kinematic distribution. This is the
    temperature of an isothermal ideal gas of primordial composition with the effective sound speed.
    !!}
    use :: Numerical_Constants_Astronomical, only : meanAtomicMassPrimordial
    use :: Numerical_Constants_Atomic      , only : atomicMassUnit
    use :: Numerical_Constants_Physical    , only : boltzmannsConstant
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    class(kinematicsDistributionEffectiveSoundSpeed), intent(inout) :: self
    class(coordinate                               ), intent(in   ) :: coordinates
    !$GLC attributes unused :: coordinates

    effectiveSoundSpeedTemperature=+meanAtomicMassPrimordial      &
         &                         *atomicMassUnit                &
         &                         *(                             &
         &                           +self%soundSpeedEffective    &
         &                           *kilo                        &
         &                          )**2                          &
         &                         /boltzmannsConstant
    return
  end function effectiveSoundSpeedTemperature

  double precision function effectiveSoundSpeedVelocityDispersion1D(self,coordinates,massDistribution_,massDistributionEmbedding)
    !!{RST
    Return the 1D velocity dispersion at the specified ``coordinates`` in an effective sound speed kinematic distribution.
    !!}
    implicit none
    class(kinematicsDistributionEffectiveSoundSpeed), intent(inout)         :: self
    class(coordinate                               ), intent(in   )         :: coordinates
    class(massDistributionClass                    ), intent(inout), target :: massDistribution_, massDistributionEmbedding
    !$GLC attributes unused :: coordinates, massDistribution_, massDistributionEmbedding

    effectiveSoundSpeedVelocityDispersion1D=self%soundSpeedEffective
    return
  end function effectiveSoundSpeedVelocityDispersion1D
