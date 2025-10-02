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
  Implementation of a kinematic distribution class for the soliton-NFW Heated mass distribution.
  !!}

  !![
  <kinematicsDistribution name="kinematicsDistributionSolitonNFWHeated">
   <description>A kinematic distribution class for the Soliton-NFW Heated mass distribution.</description>
  </kinematicsDistribution>
  !!]
  type, public, extends(kinematicsDistributionCollisionlessTabulated) :: kinematicsDistributionSolitonNFWHeated
     !!{
     A kinematics distribution for the Soliton-NFW Heated mass distribution.
     !!}
   contains
     procedure :: velocityDispersion1D => solitonNFWHeatedKinematicsVelocityDispersion1D
  end type kinematicsDistributionSolitonNFWHeated

  interface kinematicsDistributionSolitonNFWHeated
     !!{
     Constructors for the \refClass{kinematicsDistributionSolitonNFWHeated} kinematic distribution class.
     !!}
     module procedure solitonNFWHeatedKinematicsConstructorParameters
     module procedure solitonNFWHeatedKinematicsConstructorInternal
     module procedure solitonNFWHeatedKinematicsConstructorDecorated
  end interface kinematicsDistributionSolitonNFWHeated

contains

  function solitonNFWHeatedKinematicsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{kinematicsDistributionSolitonNFW} kinematic distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (kinematicsDistributionSolitonNFWHeated)                :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    double precision                                                  :: toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum

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
    self=kinematicsDistributionSolitonNFWHeated(toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function solitonNFWHeatedKinematicsConstructorParameters

  function solitonNFWHeatedKinematicsConstructorInternal(toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum) result(self)
    !!{
    Internal constructor for the \refClass{kinematicsDistributionSolitonNFWHeated} kinematic distribution class.
    !!}
    implicit none
    type            (kinematicsDistributionSolitonNFWHeated)                          :: self
    double precision                                  , intent(in   ), optional :: toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum
    !![
    <constructorAssign variables="toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum"/>
    !!]

    return
  end function solitonNFWHeatedKinematicsConstructorInternal
  
  function solitonNFWHeatedKinematicsConstructorDecorated(kinematicsDistribution_) result(self)
    !!{
    Internal constructor for the \refClass{kinematicsDistributionSolitonNFWHeated} kinematic distribution class.
    !!}
    implicit none
    type (kinematicsDistributionSolitonNFWHeated)                :: self
    class(kinematicsDistributionClass     ), intent(in   ) :: kinematicsDistribution_

    self%toleranceRelativeVelocityDispersion       =kinematicsDistribution_%toleranceRelativeVelocityDispersion
    self%toleranceRelativeVelocityDispersionMaximum=kinematicsDistribution_%toleranceRelativeVelocityDispersionMaximum
    return
  end function solitonNFWHeatedKinematicsConstructorDecorated
  
  logical function solitonNFWHeatedKinematicsIsCollisional(self)
    !!{
    Return false indicating that the soliton-NFW Heated distribution represents collisionless particles.
    !!}
    implicit none
    class(kinematicsDistributionSolitonNFWHeated), intent(inout) :: self
    !$GLC attributes unused :: self
    
    solitonNFWHeatedKinematicsIsCollisional=.false.
    return
  end function solitonNFWHeatedKinematicsIsCollisional

  double precision function solitonNFWHeatedKinematicsVelocityDispersion1D(self,coordinates,massDistribution_,massDistributionEmbedding) result(velocityDispersion)
    !!{
    Return the 1D velocity dispersion at the specified {\normalfont \ttfamily coordinates} in a soliton-NFW Heated kinematic distribution.
    !!}
    use :: Error      , only : Error_Report
    use :: Coordinates, only : coordinateSpherical, assignment(=)
    implicit none
    class           (kinematicsDistributionSolitonNFWHeated), intent(inout)          :: self
    class           (coordinate                      ), intent(in   )          :: coordinates
    class           (massDistributionClass           ), intent(inout), target  :: massDistribution_          , massDistributionEmbedding
    class           (massDistributionClass           )               , pointer :: massDistribution__
    double precision                                  , parameter              :: fractionSmall       =1.0d-2
    double precision                                                           :: radiusScaleFree            , radiusScaleFreeSmall
    type            (coordinateSpherical             )                         :: coordinatesReference

    massDistribution__ => massDistribution_
    if (associated(massDistribution__,massDistributionEmbedding)) then
       select type (massDistributionEmbedding)
       class is (massDistributionSolitonNFW)
          if (.not.massDistributionEmbedding%isTabulating()) then
             radiusScaleFree     =+coordinates              %rSpherical () &
                  &               /massDistributionEmbedding%radiusScale
             radiusScaleFreeSmall=+fractionSmall                           &
                  &               *massDistributionEmbedding%radiusCore    &
                  &               /massDistributionEmbedding%radiusScale
             if (radiusScaleFree < radiusScaleFreeSmall) then
                ! Small radius. Use the tabulated solution at the small radius boundary, extrapolated to the actual
                ! radius using the result for a power-law ρ(r) ∝ r⁰ profile.
                coordinatesReference=[radiusScaleFreeSmall*massDistributionEmbedding%radiusScale,0.0d0,0.0d0]
                velocityDispersion  =+massDistributionEmbedding%velocityDispersion1D         (coordinatesReference                                            ) &
                     &               *radiusScaleFree                                                                                                           &
                     &               /radiusScaleFreeSmall
             else
                velocityDispersion  =+massDistributionEmbedding%velocityDispersion1D         (coordinates                                                     )
             end if
          else
             velocityDispersion     =+self                     %velocityDispersion1DNumerical(coordinates         ,massDistribution_,massDistributionEmbedding)
          end if
       class default
          velocityDispersion        =+0.0d0
          call Error_Report('expecting a soliton-NFW Heated mass distribution, but received '//char(massDistributionEmbedding%objectType())//{introspection:location})
       end select
    else
       ! Our tabulated distribution is embedded in another distribution. We must compute the velocity dispersion numerically.
       velocityDispersion           =+self                     %velocityDispersion1DNumerical(coordinates         ,massDistribution_,massDistributionEmbedding)
    end if
    return
  end function solitonNFWHeatedKinematicsVelocityDispersion1D
