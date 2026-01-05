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

  !!{
  Implementation of a kinematic distribution class for collisionless mass distributions.
  !!}

  !![
  <kinematicsDistribution name="kinematicsDistributionCollisionless">
   <description>A kinematic distribution class for collisionless mass distributions.</description>
  </kinematicsDistribution>
  !!]
  type, public, extends(kinematicsDistributionClass) :: kinematicsDistributionCollisionless
     !!{
     A kinematics distribution for collisionless distributions.
     !!}
   contains
     procedure :: isCollisional        => collisionlessIsCollisional
     procedure :: velocityDispersion1D => collisionlessVelocityDispersion1D
  end type kinematicsDistributionCollisionless

  interface kinematicsDistributionCollisionless
     !!{
     Constructors for the \refClass{kinematicsDistributionCollisionless} kinematic distribution class.
     !!}
     module procedure collisionlessConstructorParameters
     module procedure collisionlessConstructorInternal
     module procedure collisionlessConstructorDecorated
  end interface kinematicsDistributionCollisionless

contains

  function collisionlessConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{kinematicsDistributionCollisionless} kinematic distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (kinematicsDistributionCollisionless)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    double precision                                                     :: toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum

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
    self=kinematicsDistributionCollisionless(toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function collisionlessConstructorParameters

  function collisionlessConstructorInternal(toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum) result(self)
    !!{
    Internal constructor for the \refClass{kinematicsDistributionCollisionless} kinematic distribution class.
    !!}
    implicit none
    type            (kinematicsDistributionCollisionless)                          :: self
    double precision                                     , intent(in   ), optional :: toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum
    !![
    <constructorAssign variables="toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum"/>
    !!]

    return
  end function collisionlessConstructorInternal
  
  function collisionlessConstructorDecorated(kinematicsDistribution_) result(self)
    !!{
    Internal constructor for the \refClass{kinematicsDistributionCollisionless} kinematic distribution class.
    !!}
    implicit none
    type (kinematicsDistributionCollisionless)                :: self
    class(kinematicsDistributionClass        ), intent(in   ) :: kinematicsDistribution_

    self%toleranceRelativeVelocityDispersion       =kinematicsDistribution_%toleranceRelativeVelocityDispersion
    self%toleranceRelativeVelocityDispersionMaximum=kinematicsDistribution_%toleranceRelativeVelocityDispersionMaximum
    return
  end function collisionlessConstructorDecorated
  
  logical function collisionlessIsCollisional(self)
    !!{
    Return false indicating that the collisionless kinematic distribution represents collisionless particles.
    !!}
    implicit none
    class(kinematicsDistributionCollisionless), intent(inout) :: self
    
    collisionlessIsCollisional=.false.
    return
  end function collisionlessIsCollisional

  double precision function collisionlessVelocityDispersion1D(self,coordinates,massDistribution_,massDistributionEmbedding) result(velocityDispersion)
    !!{
    Return the 1D velocity dispersion at the specified {\normalfont \ttfamily coordinates} in an collisionless kinematic distribution.
    !!}
    implicit none
    class(kinematicsDistributionCollisionless), intent(inout)          :: self
    class(coordinate                         ), intent(in   )          :: coordinates
    class(massDistributionClass              ), intent(inout), target  :: massDistribution_, massDistributionEmbedding

    velocityDispersion=self%velocityDispersion1DNumerical(coordinates,massDistribution_,massDistributionEmbedding)
    return
  end function collisionlessVelocityDispersion1D
