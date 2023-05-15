!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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

  !!{
  Implementation of an isothermal kinematic distribution class.
  !!}

  !![
  <kinematicsDistribution name="kinematicsDistributionIsothermal">
   <description>An isothermal kinematic distribution class masses.</description>
  </kinematicsDistribution>
  !!]
  type, public, extends(kinematicsDistributionClass) :: kinematicsDistributionIsothermal
     !!{
     An isothermal kinematic distribution.
     !!}
     double precision :: temperature_, massAtomicMean
   contains
     procedure :: isCollisional                  => isothermalIsCollisional
     procedure :: temperature                    => isothermalTemperature
     procedure :: temperatureGradientLogarithmic => isothermalTemperatureGradientLogarithmic
     procedure :: velocityDispersion1D           => isothermalVelocityDispersion1D
  end type kinematicsDistributionIsothermal

  interface kinematicsDistributionIsothermal
     !!{
     Constructors for the {\normalfont \ttfamily isothermal} kinematic distribution class.
     !!}
     module procedure isothermalConstructorParameters
     module procedure isothermalConstructorInternal
  end interface kinematicsDistributionIsothermal

contains

  function isothermalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily iosothermal} kinematic distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (kinematicsDistributionIsothermal)                :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    double precision                                                  :: temperature, massAtomicMean

    !![
    <inputParameter>
      <name>temperature</name>
      <description>The temperature of the distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massAtomicMean</name>
      <description>The mean atomic mass (in atomic mass units) of the distribution.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=kinematicsDistributionIsothermal(temperature,massAtomicMean)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function isothermalConstructorParameters
  
  function isothermalConstructorInternal(temperature_,massAtomicMean) result(self)
    !!{
    Constructor for {\normalfont \ttfamily isothermal} kinematics distribution class.
    !!}
    implicit none
    type            (kinematicsDistributionIsothermal)                :: self
    double precision                                  , intent(in   ) :: temperature_, massAtomicMean
    !![
    <constructorAssign variables="temperature_, massAtomicMean"/>
    !!]

    return
  end function isothermalConstructorInternal

  logical function isothermalIsCollisional(self)
    !!{
    Return true indicating that the isothermal kinematic distribution represents collisional particles.
    !!}
    implicit none
    class(kinematicsDistributionIsothermal), intent(inout) :: self
    
    isothermalIsCollisional=.true.
    return
  end function isothermalIsCollisional

  double precision function isothermalTemperature(self,coordinates)
    !!{
    Return the temperature at the specified {\normalfont \ttfamily coordinates} in an isothermal kinematic distribution.
    !!}
    implicit none
    class(kinematicsDistributionIsothermal), intent(inout) :: self
    class(coordinate                      ), intent(in   ) :: coordinates
    !$GLC attributes unused :: coordinates

    isothermalTemperature=self%temperature_
    return
  end function isothermalTemperature

  double precision function isothermalTemperatureGradientLogarithmic(self,coordinates)
    !!{
    Return the logarithmic gradient of temperature at the specified {\normalfont \ttfamily coordinates} in an isothermal kinematic distribution.
    !!}
    implicit none
    class(kinematicsDistributionIsothermal), intent(inout) :: self
    class(coordinate                      ), intent(in   ) :: coordinates
    !$GLC attributes unused :: coordinates, self

    isothermalTemperatureGradientLogarithmic=0.0d0
    return
  end function isothermalTemperatureGradientLogarithmic

  double precision function isothermalVelocityDispersion1D(self,coordinates,massDistributionEmbedding)
    !!{
    Return the 1D velocity dispersion at the specified {\normalfont \ttfamily coordinates} in an isothermal kinematic distribution.
    !!}
    use :: Numerical_Constants_Atomic  , only : atomicMassUnit
    use :: Numerical_Constants_Physical, only : boltzmannsConstant
    use :: Numerical_Constants_Prefixes, only : kilo
    implicit none
    class(kinematicsDistributionIsothermal), intent(inout) :: self
    class(coordinate                      ), intent(in   ) :: coordinates
    class(massDistributionClass           ), intent(inout) :: massDistributionEmbedding
    !$GLC attributes unused :: coordinates, massDistributionEmbedding
    
    isothermalVelocityDispersion1D=+sqrt(                         &
         &                               +     boltzmannsConstant &
         &                               *self%temperature_       &
         &                               /self%massAtomicMean     &
         &                               /     atomicMassUnit     &
         &                              )                         &
         &                          /kilo
    return
  end function isothermalVelocityDispersion1D
