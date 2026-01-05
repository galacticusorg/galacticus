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
  Implementation of a kinematic distribution class for the cusp-NFW \citep{delos_cusp-halo_2025} mass distribution.
  !!}

  !![
  <kinematicsDistribution name="kinematicsDistributionCuspNFW">
   <description>A kinematic distribution class for the cusp-NFW \citep{delos_cusp-halo_2025} mass distribution.</description>
  </kinematicsDistribution>
  !!]
  type, public, extends(kinematicsDistributionCollisionlessTabulated) :: kinematicsDistributionCuspNFW
     !!{
     A kinematics distribution for the cusp-NFW \citep{delos_cusp-halo_2025} mass distribution.
     !!}
   contains
     procedure :: velocityDispersion1D => cuspNFWKinematicsVelocityDispersion1D
  end type kinematicsDistributionCuspNFW

  interface kinematicsDistributionCuspNFW
     !!{
     Constructors for the \refClass{kinematicsDistributionCuspNFW} kinematic distribution class.
     !!}
     module procedure cuspNFWKinematicsConstructorParameters
     module procedure cuspNFWKinematicsConstructorInternal
     module procedure cuspNFWKinematicsConstructorDecorated
  end interface kinematicsDistributionCuspNFW

contains

  function cuspNFWKinematicsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{kinematicsDistributionCuspNFW} kinematic distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (kinematicsDistributionCuspNFW)                :: self
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
    self=kinematicsDistributionCuspNFW(toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function cuspNFWKinematicsConstructorParameters

  function cuspNFWKinematicsConstructorInternal(toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum) result(self)
    !!{
    Internal constructor for the \refClass{kinematicsDistributionCuspNFW} kinematic distribution class.
    !!}
    implicit none
    type            (kinematicsDistributionCuspNFW)                          :: self
    double precision                               , intent(in   ), optional :: toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum
    !![
    <constructorAssign variables="toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum"/>
    !!]

    return
  end function cuspNFWKinematicsConstructorInternal
  
  function cuspNFWKinematicsConstructorDecorated(kinematicsDistribution_) result(self)
    !!{
    Internal constructor for the \refClass{kinematicsDistributionCuspNFW} kinematic distribution class.
    !!}
    implicit none
    type (kinematicsDistributionCuspNFW)                :: self
    class(kinematicsDistributionClass  ), intent(in   ) :: kinematicsDistribution_

    self%toleranceRelativeVelocityDispersion       =kinematicsDistribution_%toleranceRelativeVelocityDispersion
    self%toleranceRelativeVelocityDispersionMaximum=kinematicsDistribution_%toleranceRelativeVelocityDispersionMaximum
    return
  end function cuspNFWKinematicsConstructorDecorated
  
  logical function cuspNFWKinematicsIsCollisional(self)
    !!{
    Return false indicating that the cusp-NFW distribution represents collisionless particles.
    !!}
    implicit none
    class(kinematicsDistributionCuspNFW), intent(inout) :: self
    !$GLC attributes unused :: self
    
    cuspNFWKinematicsIsCollisional=.false.
    return
  end function cuspNFWKinematicsIsCollisional

  double precision function cuspNFWKinematicsVelocityDispersion1D(self,coordinates,massDistribution_,massDistributionEmbedding) result(velocityDispersion)
    !!{
    Return the 1D velocity dispersion at the specified {\normalfont \ttfamily coordinates} in a cusp-NFW kinematic distribution.
    !!}
    use :: Error      , only : Error_Report
    use :: Coordinates, only : coordinateSpherical, assignment(=)
    implicit none
    class           (kinematicsDistributionCuspNFW), intent(inout)          :: self
    class           (coordinate                   ), intent(in   )          :: coordinates
    class           (massDistributionClass        ), intent(inout), target  :: massDistribution_          , massDistributionEmbedding
    class           (massDistributionClass        )               , pointer :: massDistribution__
    double precision                                , parameter             :: fractionSmall       =1.0d-2
    double precision                                                        :: radiusScaleFree            , radiusScaleFreeSmall
    type            (coordinateSpherical          )                         :: coordinatesReference

    massDistribution__ => massDistribution_
    if (associated(massDistribution__,massDistributionEmbedding)) then
       select type (massDistributionEmbedding)
       class is (massDistributionCuspNFW)
          if (.not.massDistributionEmbedding%isTabulating()) then
             radiusScaleFree     =+coordinates              %rSpherical ()    &
                  &               /massDistributionEmbedding%radiusScale
             radiusScaleFreeSmall=+fractionSmall                              &
                  &               *massDistributionEmbedding%y            **2
             if (radiusScaleFree < radiusScaleFreeSmall) then
                ! Small radius. Use the tabulated solution at the small radius boundary, extrapolated to the actual
                ! radius using the result for a power-law ρ(r) ∝ r^{-3/2} profile.
                coordinatesReference=[radiusScaleFreeSmall*massDistributionEmbedding%radiusScale,0.0d0,0.0d0]
                velocityDispersion  =+massDistributionEmbedding%velocityDispersion1D         (coordinatesReference                                            ) &
                     &               *(                                                                                                                         &
                     &                 +radiusScaleFree                                                                                                         &
                     &                 /radiusScaleFreeSmall                                                                                                    &
                     &                )**0.25d0
             else
                velocityDispersion  =+massDistributionEmbedding%velocityDispersion1D         (coordinates                                                     )
             end if
          else
             velocityDispersion     =+self                     %velocityDispersion1DNumerical(coordinates         ,massDistribution_,massDistributionEmbedding)
          end if
       class default
          velocityDispersion        =+0.0d0
          call Error_Report('expecting a cusp-NFW mass distribution, but received '//char(massDistributionEmbedding%objectType())//{introspection:location})
       end select
    else
       ! Our tabulated distribution is embedded in another distribution. We must compute the velocity dispersion numerically.
       velocityDispersion           =+self                     %velocityDispersion1DNumerical(coordinates         ,massDistribution_,massDistributionEmbedding)
    end if
    return
  end function cuspNFWKinematicsVelocityDispersion1D
