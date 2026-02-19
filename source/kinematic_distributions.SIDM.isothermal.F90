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
  Provides a kinematic distribution class implementing the ``isothermal'' approximation to the effects of SIDM based on the model
  of \cite{jiang_semi-analytic_2023}.
  !!}

  !![
  <kinematicsDistribution name="kinematicsDistributionSIDMIsothermal">
    <description>
      A kinematic distribution class implementing the ``isothermal'' approximation to the effects of SIDM based on the model of
      \cite{jiang_semi-analytic_2023}.
    </description>
  </kinematicsDistribution>
  !!]
  type, public, extends(kinematicsDistributionClass) :: kinematicsDistributionSIDMIsothermal
     !!{
     A kinematic distribution class implementing the ``isothermal'' approximation to the effects of SIDM based on the model
     of \cite{jiang_semi-analytic_2023}.
     !!}
   contains
     procedure :: isCollisional        => sidmIsothermalIsCollisional
     procedure :: velocityDispersion1D => sidmIsothermalVelocityDispersion1D
  end type kinematicsDistributionSIDMIsothermal

  interface kinematicsDistributionSIDMIsothermal
     !!{
     Constructors for the \refClass{kinematicsDistributionSIDMIsothermal} kinematic distribution class.
     !!}
     module procedure sidmIsothermalConstructorParameters
     module procedure sidmIsothermalConstructorInternal
  end interface kinematicsDistributionSIDMIsothermal

contains

  function sidmIsothermalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{kinematicsDistributionSIDMIsothermal} kinematic distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(kinematicsDistributionSIDMIsothermal)                :: self
    type(inputParameters                     ), intent(inout) :: parameters
    
    self=kinematicsDistributionSIDMIsothermal()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function sidmIsothermalConstructorParameters

  function sidmIsothermalConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{kinematicsDistributionSIDMIsothermal} kinematic distribution class.
    !!}
    implicit none
    type(kinematicsDistributionSIDMIsothermal) :: self
    
    return
  end function sidmIsothermalConstructorInternal
  
  logical function sidmIsothermalIsCollisional(self)
    !!{
    Return true indicating that the sidmIsothermal kinematic distribution represents collisional particles.
    !!}
    implicit none
    class(kinematicsDistributionSIDMIsothermal), intent(inout) :: self
    
    sidmIsothermalIsCollisional=.false.
    return
  end function sidmIsothermalIsCollisional

  double precision function sidmIsothermalVelocityDispersion1D(self,coordinates,massDistribution_,massDistributionEmbedding) result(velocityDispersion)
    !!{
    Return the 1D velocity dispersion at the specified {\normalfont \ttfamily coordinates} in an SIDMIsothermal kinematic distribution.
    !!}
    use :: ISO_Varying_String, only : char
    implicit none
    class(kinematicsDistributionSIDMIsothermal), intent(inout)          :: self
    class(coordinate                          ), intent(in   )          :: coordinates
    class(massDistributionClass               ), intent(inout), target  :: massDistribution_ , massDistributionEmbedding
    class(massDistributionClass               )               , pointer :: massDistribution__
  
    massDistribution__ => massDistribution_
    if (associated(massDistribution__,massDistributionEmbedding)) then
       ! For the case of a self-gravitating SIDM isothermal distribution we have a constant velocity dispersion in the core region.
       select type (massDistributionEmbedding)
       class is (massDistributionSphericalSIDMIsothermal       )
          if (coordinates%rSpherical() > massDistributionEmbedding%radiusInteraction()) then
             velocityDispersion=self                     %velocityDispersion1DNumerical(coordinates,massDistribution_,massDistributionEmbedding)
          else
             velocityDispersion=massDistributionEmbedding%velocityDispersionCentral
          end if
       class is (massDistributionSphericalSIDMIsothermalBaryons)
          if (coordinates%rSpherical() > massDistributionEmbedding%radiusInteraction()) then
             velocityDispersion=self                     %velocityDispersion1DNumerical(coordinates,massDistribution_,massDistributionEmbedding)
          else
             velocityDispersion=massDistributionEmbedding%velocityDispersionCentral
          end if
       class default
          velocityDispersion=0.0d0
          call Error_Report('expecting an SIDMIsothermal mass distribution, but received '//char(massDistributionEmbedding%objectType())//{introspection:location})
       end select
    else
       ! Our SIDM isothermal distribution is embedded in another distribution. We must compute the velocity dispersion numerically.
       velocityDispersion=self%velocityDispersion1DNumerical(coordinates,massDistribution_,massDistributionEmbedding)
    end if
    return
  end function sidmIsothermalVelocityDispersion1D
