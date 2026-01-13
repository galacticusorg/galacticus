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

  !+    Contributions to this file made by: Yu Zhao

  !!{
  Implementation of a kinematic distribution class for the soliton mass distribution.
  !!}

  !![
  <kinematicsDistribution name="kinematicsDistributionSoliton">
    <description>
      A kinematic distribution class for the soliton mass distribution.
    </description>
  </kinematicsDistribution>
  !!]
  type, public, extends(kinematicsDistributionCollisionlessTabulated) :: kinematicsDistributionSoliton
     !!{
     A kinematics distribution for the Soliton mass distribution.
     !!}
   contains
     procedure :: velocityDispersion1D => solitonKinematicsVelocityDispersion1D
  end type kinematicsDistributionSoliton

  interface kinematicsDistributionSoliton
     !!{
     Constructors for the \refClass{kinematicsDistributionSoliton} kinematic distribution class.
     !!}
     module procedure solitonKinematicsConstructorParameters
     module procedure solitonKinematicsConstructorInternal
     module procedure solitonKinematicsConstructorDecorated
  end interface kinematicsDistributionSoliton

  ! Coefficient of the dimensionless radius in the soliton profile.
   double precision, parameter :: coefficientCore=0.091d0 ! Schive et al. (2014; https://ui.adsabs.harvard.edu/abs/2014PhRvL.113z1302S; equation 3).

contains

  function solitonKinematicsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{kinematicsDistributionSoliton} kinematic distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (kinematicsDistributionSoliton)                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    double precision                                               :: toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum

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
    self=kinematicsDistributionSoliton(toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function solitonKinematicsConstructorParameters

  function solitonKinematicsConstructorInternal(toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum) result(self)
    !!{
    Internal constructor for the \refClass{kinematicsDistributionSoliton} kinematic distribution class.
    !!}
    implicit none
    type            (kinematicsDistributionSoliton)                          :: self
    double precision                               , intent(in   ), optional :: toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum
    !![
    <constructorAssign variables="toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum"/>
    !!]

    return
  end function solitonKinematicsConstructorInternal
  
  function solitonKinematicsConstructorDecorated(kinematicsDistribution_) result(self)
    !!{
    Internal constructor for the \refClass{kinematicsDistributionSoliton} kinematic distribution class.
    !!}
    implicit none
    type (kinematicsDistributionSoliton)                :: self
    class(kinematicsDistributionClass  ), intent(in   ) :: kinematicsDistribution_

    self%toleranceRelativeVelocityDispersion       =kinematicsDistribution_%toleranceRelativeVelocityDispersion
    self%toleranceRelativeVelocityDispersionMaximum=kinematicsDistribution_%toleranceRelativeVelocityDispersionMaximum
    return
  end function solitonKinematicsConstructorDecorated
  
  logical function solitonKinematicsIsCollisional(self)
    !!{
    Return false indicating that the soliton distribution represents collisionless particles.
    !!}
    implicit none
    class(kinematicsDistributionSoliton), intent(inout) :: self
    !$GLC attributes unused :: self
    
    solitonKinematicsIsCollisional=.false.
    return
  end function solitonKinematicsIsCollisional

  double precision function solitonKinematicsVelocityDispersion1D(self,coordinates,massDistribution_,massDistributionEmbedding) result(velocityDispersion)
    !!{
    Return the 1D velocity dispersion at the specified {\normalfont \ttfamily coordinates} in a soliton kinematic distribution.
    !!}
    use :: Error      , only : Error_Report
    use :: Coordinates, only : coordinateSpherical, assignment(=)
    implicit none
    class(kinematicsDistributionSoliton), intent(inout)         :: self
    class(coordinate                   ), intent(in   )         :: coordinates
    class(massDistributionClass        ), intent(inout), target :: massDistribution_, massDistributionEmbedding

    select type (massDistributionEmbedding)
    class is (massDistributionSoliton)
        velocityDispersion=+self%velocityDispersion1DNumerical(coordinates,massDistribution_,massDistributionEmbedding)
    class default
        velocityDispersion=+0.0d0
        call Error_Report('expecting a soliton mass distribution, but received '//char(massDistributionEmbedding%objectType())//{introspection:location})
    end select
    return
  end function solitonKinematicsVelocityDispersion1D
