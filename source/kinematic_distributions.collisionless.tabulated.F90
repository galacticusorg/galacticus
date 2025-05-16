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
  Implementation of a kinematic distribution class for collisionless mass distributions using tabulated solutions.
  !!}

  !![
  <kinematicsDistribution name="kinematicsDistributionCollisionlessTabulated">
   <description>A kinematic distribution class for collisionless mass distributions using tabulated solutions.</description>
  </kinematicsDistribution>
  !!]
  type, public, extends(kinematicsDistributionClass) :: kinematicsDistributionCollisionlessTabulated
     !!{
     A kinematics distribution for collisionless distributions using tabulated solutions.
     !!}
   contains
     procedure :: isCollisional        => collisionlessTabulatedIsCollisional
     procedure :: velocityDispersion1D => collisionlessTabulatedVelocityDispersion1D
  end type kinematicsDistributionCollisionlessTabulated

  interface kinematicsDistributionCollisionlessTabulated
     !!{
     Constructors for the \refClass{kinematicsDistributionCollisionlessTabulated} kinematic distribution class.
     !!}
     module procedure collisionlessTabulatedConstructorParameters
     module procedure collisionlessTabulatedConstructorInternal
     module procedure collisionlessTabulatedConstructorDecorated
  end interface kinematicsDistributionCollisionlessTabulated

contains

  function collisionlessTabulatedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{kinematicsDistributionCollisionlessTabulated} kinematic distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (kinematicsDistributionCollisionlessTabulated)                :: self
    type            (inputParameters                             ), intent(inout) :: parameters
    double precision                                                              :: toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum

    !![
    <inputParameter>
      <name>toleranceRelativeVelocityDispersion</name>
      <defaultValue>1.0d-6</defaultValue>
      <source>parameters</source>
      <description>The relative tolerance to use in numerical solutions for the velocity dispersion in dark-matter-only density profiles.</description>
    </inputParameter>
    <inputParameter>
      <name>toleranceRelativeVelocityDispersionMaximum</name>
      <defaultValue>1.0d-3</defaultValue>
      <source>parameters</source>
      <description>The maximum relative tolerance to use in numerical solutions for the velocity dispersion in dark-matter-only density profiles.</description>
    </inputParameter>
    !!]
    self=kinematicsDistributionCollisionlessTabulated(toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function collisionlessTabulatedConstructorParameters

  function collisionlessTabulatedConstructorInternal(toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum) result(self)
    !!{
    Internal constructor for the \refClass{kinematicsDistributionCollisionlessTabulated} kinematic distribution class.
    !!}
    implicit none
    type            (kinematicsDistributionCollisionlessTabulated)                          :: self
    double precision                                              , intent(in   ), optional :: toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum
    !![
    <constructorAssign variables="toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum"/>
    !!]

    return
  end function collisionlessTabulatedConstructorInternal
  
  function collisionlessTabulatedConstructorDecorated(kinematicsDistribution_) result(self)
    !!{
    Internal constructor for the \refClass{kinematicsDistributionCollisionlessTabulated} kinematic distribution class.
    !!}
    implicit none
    type (kinematicsDistributionCollisionlessTabulated)                :: self
    class(kinematicsDistributionClass                 ), intent(in   ) :: kinematicsDistribution_

    self%toleranceRelativeVelocityDispersion       =kinematicsDistribution_%toleranceRelativeVelocityDispersion
    self%toleranceRelativeVelocityDispersionMaximum=kinematicsDistribution_%toleranceRelativeVelocityDispersionMaximum
    return
  end function collisionlessTabulatedConstructorDecorated
  
  logical function collisionlessTabulatedIsCollisional(self)
    !!{
    Return false indicating that the tabulated collisionless kinematic distribution represents collisionless particles.
    !!}
    implicit none
    class(kinematicsDistributionCollisionlessTabulated), intent(inout) :: self
    
    collisionlessTabulatedIsCollisional=.false.
    return
  end function collisionlessTabulatedIsCollisional

  double precision function collisionlessTabulatedVelocityDispersion1D(self,coordinates,massDistribution_,massDistributionEmbedding) result(velocityDispersion)
    !!{
    Return the 1D velocity dispersion at the specified {\normalfont \ttfamily coordinates} in a tabulated collisionless kinematic distribution.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(kinematicsDistributionCollisionlessTabulated), intent(inout)          :: self
    class(coordinate                                  ), intent(in   )          :: coordinates
    class(massDistributionClass                       ), intent(inout), target  :: massDistribution_ , massDistributionEmbedding
    class(massDistributionClass                       )               , pointer :: massDistribution__

    massDistribution__ => massDistribution_
    if (associated(massDistribution__,massDistributionEmbedding)) then
       select type (massDistributionEmbedding)
       class is (massDistributionSphericalTabulated)
          if (.not.massDistributionEmbedding%isTabulating()) then
             velocityDispersion=massDistributionEmbedding%velocityDispersion1D         (coordinates                          )
          else
             velocityDispersion=self                     %velocityDispersion1DNumerical(coordinates,massDistribution_,massDistributionEmbedding)
          end if
       class default
          velocityDispersion   =0.0d0
          call Error_Report('expecting a tabulated mass distribution, but received '//char(massDistributionEmbedding%objectType())//{introspection:location})
       end select
    else
       ! Our tabulated distribution is embedded in another distribution. We must compute the velocity dispersion numerically.
       velocityDispersion      =self                     %velocityDispersion1DNumerical(coordinates,massDistribution_,massDistributionEmbedding)
    end if
    return
  end function collisionlessTabulatedVelocityDispersion1D
