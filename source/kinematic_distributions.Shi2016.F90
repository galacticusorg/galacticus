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
  Implementation of a kinematic distribution class for the \cite{shi_outer_2016} model of halo accretion flows.
  !!}

  !![
  <kinematicsDistribution name="kinematicsDistributionShi2016">
   <description>A kinematic distribution class for the \cite{shi_outer_2016} model of halo accretion flows.</description>
  </kinematicsDistribution>
  !!]
  type, public, extends(kinematicsDistributionCollisionless) :: kinematicsDistributionShi2016
     !!{
     A kinematics distribution for the \cite{shi_outer_2016} model of halo accretion flows.
     !!}
   contains
     procedure :: velocityRadial => shi2016KinematicsVelocityRadial
  end type kinematicsDistributionShi2016

  interface kinematicsDistributionShi2016
     !!{
     Constructors for the \refClass{kinematicsDistributionShi2016} kinematic distribution class.
     !!}
     module procedure shi2016KinematicsConstructorParameters
     module procedure shi2016KinematicsConstructorInternal
  end interface kinematicsDistributionShi2016

contains

  function shi2016KinematicsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{kinematicsDistributionShi2016} kinematic distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(kinematicsDistributionShi2016)                :: self
    type(inputParameters              ), intent(inout) :: parameters

    self=kinematicsDistributionShi2016()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function shi2016KinematicsConstructorParameters

  function shi2016KinematicsConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{kinematicsDistributionShi2016} kinematic distribution class.
    !!}
    implicit none
    type(kinematicsDistributionShi2016) :: self
    
    return
  end function shi2016KinematicsConstructorInternal
  
  double precision function shi2016KinematicsVelocityRadial(self,coordinates,massDistributionEmbedding) result(velocityRadial)
    !!{
    Return the radial velocity at the specified {\normalfont \ttfamily coordinates} in the \cite{shi_outer_2016} model for the accretion flow around a halo.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(kinematicsDistributionShi2016), intent(inout) :: self
    class(coordinate                   ), intent(in   ) :: coordinates
    class(massDistributionClass        ), intent(inout) :: massDistributionEmbedding

    select type (massDistributionEmbedding)
    class is (massDistributionShi2016)
       if      (coordinates%rSpherical() > massDistributionEmbedding%radiusMaximumPhysical) then
          ! Beyond the maximum radius for the flow just return the mean matter velocity.
          velocityRadial=+massDistributionEmbedding%cosmologyFunctions_         %hubbleParameterEpochal(massDistributionEmbedding%time        ) &
               &         *coordinates                                           %rSpherical            (                                      )
       else if (coordinates%rSpherical() < massDistributionEmbedding%radiusMinimumPhysical) then
          velocityRadial=+0.0d0
          call Error_Report('radius is less than minimum tabulated for accretion flow'//{introspection:location})
       else
          velocityRadial=+massDistributionEmbedding%interpolatorVelocityPhysical%interpolate           (coordinates              %rSpherical())
       end if
       velocityRadial=+velocityRadial                                 &
            &          *massDistributionEmbedding%scaleFactorVelocity
    class default
       velocityRadial=0.0d0
       call Error_Report('expecting a Shi2016 mass distribution, but received '//char(massDistributionEmbedding%objectType())//{introspection:location})
    end select
    return
  end function shi2016KinematicsVelocityRadial
