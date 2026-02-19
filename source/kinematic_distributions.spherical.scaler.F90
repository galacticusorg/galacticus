!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  Implementation of a spherical scaler kinematic distribution class.
  !!}

  !![
  <kinematicsDistribution name="kinematicsDistributionSphericalScaler">
   <description>A spherical scaler kinematic distribution class masses.</description>
  </kinematicsDistribution>
  !!]
  type, public, extends(kinematicsDistributionClass) :: kinematicsDistributionSphericalScaler
     !!{
     An spherical scaler kinematic distribution.
     !!}
     class           (kinematicsDistributionClass), pointer :: kinematicsDistribution_ => null()
     double precision                                       :: factorScalingLength              , factorScalingMass
   contains
     final     ::                         kinematicsSphericalScalerDestructor
     procedure :: isCollisional        => kinematicsSphericalScalerIsCollisional
     procedure :: velocityDispersion1D => kinematicsSphericalScalerVelocityDispersion1D
  end type kinematicsDistributionSphericalScaler

  interface kinematicsDistributionSphericalScaler
     !!{
     Constructors for the \refClass{kinematicsDistributionSphericalScaler} kinematic distribution class.
     !!}
     module procedure kinematicsSphericalScalerConstructorParameters
     module procedure kinematicsSphericalScalerConstructorInternal
  end interface kinematicsDistributionSphericalScaler

contains

  function kinematicsSphericalScalerConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{kinematicsDistributionSphericalScaler} kinematic distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (kinematicsDistributionSphericalScaler)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    class           (kinematicsDistributionClass          ), pointer       :: kinematicsDistribution_
    double precision                                                       :: factorScalingLength    , factorScalingMass

    !![
    <inputParameter>
      <name>factorScalingLength</name>
      <description>The factor by which to scale lengths.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>factorScalingMass</name>
      <description>The factor by which to scale the mass.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="kinematicsDistribution" name="kinematicsDistribution_" source="parameters"/>
    !!]
    self=kinematicsDistributionSphericalScaler(factorScalingLength,factorScalingMass,kinematicsDistribution_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="kinematicsDistribution_"/>
    !!]
    return
  end function kinematicsSphericalScalerConstructorParameters
    
  function kinematicsSphericalScalerConstructorInternal(factorScalingLength,factorScalingMass,kinematicsDistribution_) result(self)
    !!{
    Constructor for the \refClass{kinematicsDistributionSphericalScaler} convergence class.
    !!}
    implicit none
    type            (kinematicsDistributionSphericalScaler)                        :: self
    class           (kinematicsDistributionClass          ), intent(in   ), target :: kinematicsDistribution_
    double precision                                       , intent(in   )         :: factorScalingLength    , factorScalingMass
    !![
    <constructorAssign variables="factorScalingLength, factorScalingMass, *kinematicsDistribution_"/>
    !!]
 
    return
  end function kinematicsSphericalScalerConstructorInternal

  subroutine kinematicsSphericalScalerDestructor(self)
    !!{
    Destructor for the \refClass{kinematicsDistributionSphericalScaler} mass distribution class.
    !!}
    implicit none
    type(kinematicsDistributionSphericalScaler), intent(inout) :: self

    !![
    <objectDestructor name="self%kinematicsDistribution_"/>
    !!]
    return
  end subroutine kinematicsSphericalScalerDestructor

  logical function kinematicsSphericalScalerIsCollisional(self)
    !!{
    Return true indicating that the spherical scaler kinematic distribution represents collisional particles.
    !!}
    implicit none
    class(kinematicsDistributionSphericalScaler), intent(inout) :: self
    
    kinematicsSphericalScalerIsCollisional=.false.
    return
  end function kinematicsSphericalScalerIsCollisional

  double precision function kinematicsSphericalScalerVelocityDispersion1D(self,coordinates,massDistribution_,massDistributionEmbedding) result(velocityDispersion)
    !!{
    Return the 1D velocity dispersion at the specified {\normalfont \ttfamily coordinates} in a spherical scaler kinematic distribution.
    !!}
    implicit none
    class(kinematicsDistributionSphericalScaler), intent(inout)          :: self
    class(coordinate                           ), intent(in   )          :: coordinates
    class(massDistributionClass                ), intent(inout), target  :: massDistribution_  , massDistributionEmbedding
    class(massDistributionClass                )               , pointer :: massDistribution__

    massDistribution__ => massDistribution_
    if (associated(massDistribution__,massDistributionEmbedding)) then
       ! For the case of a self-gravitating scaled distribution we have an analytic solution for the velocity dispersion.
       select type (massDistributionEmbedding)
       class is (massDistributionSphericalScaler)
          velocityDispersion=+sqrt(                                                                                                                                                  &
               &                   +self%factorScalingMass                                                                                                                           &
               &                   /self%factorScalingLength                                                                                                                         &
               &                  )                                                                                                                                                  &
               &             *self%kinematicsDistribution_%velocityDispersion1D(coordinates,massDistributionEmbedding%massDistribution_,massDistributionEmbedding%massDistribution_)
       class default
          velocityDispersion=0.0d0
          call Error_Report('expecting a spherical scaler mass distribution, but received '//char(massDistributionEmbedding%objectType())//{introspection:location})
       end select
    else
       ! Our scaled distribution is embedded in another distribution. We must compute the velocity dispersion numerically.
       velocityDispersion=self%velocityDispersion1DNumerical(coordinates,massDistribution_,massDistributionEmbedding)
    end if
    return
  end function kinematicsSphericalScalerVelocityDispersion1D
