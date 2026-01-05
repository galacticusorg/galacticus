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
  Implementation of a kinematic distribution class for decorated mass distributions, that uses the undecorated mass distribution.
  !!}

  !![
  <kinematicsDistribution name="kinematicsDistributionUndecorator">
   <description>A kinematic distribution class for decorated mass distributions, that uses the undecorated mass distribution.</description>
  </kinematicsDistribution>
  !!]
  type, public, extends(kinematicsDistributionClass) :: kinematicsDistributionUndecorator
     !!{
     A kinematics distribution for decorated mass distributions, that uses the undecorated mass distribution.
     !!}
     class(kinematicsDistributionClass), pointer :: kinematicsDistribution_ => null()
   contains
     final     ::                         undecoratorDestructor
     procedure :: isCollisional        => undecoratorIsCollisional
     procedure :: velocityDispersion1D => undecoratorVelocityDispersion1D
  end type kinematicsDistributionUndecorator

  interface kinematicsDistributionUndecorator
     !!{
     Constructors for the \refClass{kinematicsDistributionUndecorator} kinematic distribution class.
     !!}
     module procedure undecoratorConstructorParameters
     module procedure undecoratorConstructorInternal
  end interface kinematicsDistributionUndecorator

contains

  function undecoratorConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{kinematicsDistributionUndecorator} kinematic distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (kinematicsDistributionUndecorator)                :: self
    type (inputParameters                  ), intent(inout) :: parameters
    class(kinematicsDistributionClass      ), pointer       :: kinematicsDistribution_

    !![
    <objectBuilder class="kinematicsDistribution" name="kinematicsDistribution_" source="parameters"/>
    !!]
    self=kinematicsDistributionUndecorator(kinematicsDistribution_)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function undecoratorConstructorParameters

  function undecoratorConstructorInternal(kinematicsDistribution_) result(self)
    !!{
    Internal constructor for the \refClass{kinematicsDistributionUndecorator} kinematic distribution class.
    !!}
    implicit none
    type (kinematicsDistributionUndecorator)                        :: self
    class(kinematicsDistributionClass      ), intent(in   ), target :: kinematicsDistribution_
    !![
    <constructorAssign variables="*kinematicsDistribution_"/>
    !!]

    return
  end function undecoratorConstructorInternal
  
  subroutine undecoratorDestructor(self)
    !!{
    Destructor for the \refClass{kinematicsDistributionUndecorator} kinematic distribution class.
    !!}
    implicit none
    type(kinematicsDistributionUndecorator), intent(inout) :: self

    !![
    <objectDestructor name="self%kinematicsDistribution_"/>
    !!]
    return
  end subroutine undecoratorDestructor

  logical function undecoratorIsCollisional(self)
    !!{
    Return whether undecorator kinematic distribution represents collisional particles.
    !!}
    implicit none
    class(kinematicsDistributionUndecorator), intent(inout) :: self
    
    undecoratorIsCollisional=self%kinematicsDistribution_%isCollisional()
    return
  end function undecoratorIsCollisional

  double precision function undecoratorVelocityDispersion1D(self,coordinates,massDistribution_,massDistributionEmbedding) result(velocityDispersion)
    !!{
    Return the 1D velocity dispersion at the specified {\normalfont \ttfamily coordinates} in an undecorator kinematic distribution.
    !!}
    implicit none
    class(kinematicsDistributionUndecorator), intent(inout)          :: self
    class(coordinate                       ), intent(in   )          :: coordinates
    class(massDistributionClass            ), intent(inout), target  :: massDistribution_ , massDistributionEmbedding
    class(massDistributionClass            )               , pointer :: massDistribution__

    massDistribution__ => massDistribution_
    if (associated(massDistribution__,massDistributionEmbedding)) then
        ! For the case of a self-gravitating distribution we can use the undecorated kinematic distribution in its own mass distribution.
       select type (massDistributionEmbedding)
       class is (massDistributionSphericalDecorator)
          velocityDispersion=self%kinematicsDistribution_%velocityDispersion1D(coordinates,massDistributionEmbedding%massDistribution_,massDistributionEmbedding%massDistribution_)
       class default
          velocityDispersion=+0.0d0
          call Error_Report('mass distribution must be of the `massDistributionSphericalDecorator` class but found `'//char(massDistributionEmbedding%objectType())//'`'//{introspection:location})
       end select
    else
       ! Our distribution is embedded in another distribution. We must compute the velocity dispersion numerically
       velocityDispersion=self%velocityDispersion1DNumerical(coordinates,massDistribution_,massDistributionEmbedding)
    end if
    return
  end function undecoratorVelocityDispersion1D
