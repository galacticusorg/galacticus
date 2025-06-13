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
Implements an fixed critical overdensity class.
!!}

  !![
  <criticalOverdensity name="criticalOverdensityFixed">
   <description>
    A critical overdensity class in which the critical overdensity is set to a fixed number given by {\normalfont \ttfamily
    [criticalOverdensity]}.
   </description>
  </criticalOverdensity>
  !!]
  type, extends(criticalOverdensityClass) :: criticalOverdensityFixed
     !!{
     A fixed critical overdensity class.
     !!}
     private
     double precision :: criticalOverdensity_
    contains
     final     ::                    fixedDestructor
     procedure :: value           => fixedValue
     procedure :: gradientTime    => fixedGradientTime
     procedure :: gradientMass    => fixedGradientMass
     procedure :: isMassDependent => fixedIsMassDependent
     procedure :: isNodeDependent => fixedIsNodeDependent
     procedure :: isTreeDependent => fixedIsTreeDependent
  end type criticalOverdensityFixed

  interface criticalOverdensityFixed
     !!{
     Constructors for the fixed critical overdensity class.
     !!}
     module procedure fixedConstructorParameters
     module procedure fixedConstructorInternal
  end interface criticalOverdensityFixed

contains

  function fixedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the fixed critical overdensity class which takes a parameter set as input.
    !!}
    use :: Input_Parameters        , only : inputParameter, inputParameters
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (criticalOverdensityFixed     )                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    double precision                                               :: criticalOverdensity_
    class           (cosmologyFunctionsClass      ), pointer       :: cosmologyFunctions_
    class           (linearGrowthClass            ), pointer       :: linearGrowth_
    class           (cosmologicalMassVarianceClass), pointer       :: cosmologicalMassVariance_

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>criticalOverdensity</name>
      <variable>criticalOverdensity_</variable>
      <source>parameters</source>
      <defaultValue>(3.0d0/20.0d0)*(12.0d0*Pi)**(2.0d0/3.0d0)</defaultValue>
      <description>The value to use for the critical overdensity for collapse of dark matter halos when using a fixed value.</description>
    </inputParameter>
    <objectBuilder class="linearGrowth"             name="linearGrowth_"             source="parameters"/>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !!]
    self=criticalOverdensityFixed(criticalOverdensity_,linearGrowth_,cosmologyFunctions_,cosmologicalMassVariance_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="linearGrowth_"            />
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="cosmologicalMassVariance_"/>
    !!]
    return
  end function fixedConstructorParameters

  function fixedConstructorInternal(criticalOverdensity_,linearGrowth_,cosmologyFunctions_,cosmologicalMassVariance_) result(self)
    !!{
    Internal constructor for the fixed critical overdensity class.
    !!}
    implicit none
    type            (criticalOverdensityFixed     )                        :: self
    double precision                                       , intent(in   ) :: criticalOverdensity_
    class           (cosmologyFunctionsClass      ), target, intent(in   ) :: cosmologyFunctions_
    class           (linearGrowthClass            ), target, intent(in   ) :: linearGrowth_
    class           (cosmologicalMassVarianceClass), target, intent(in   ) :: cosmologicalMassVariance_
    !![
    <constructorAssign variables="criticalOverdensity_, *linearGrowth_, *cosmologyFunctions_, *cosmologicalMassVariance_"/>
    !!]

    return
  end function fixedConstructorInternal

  subroutine fixedDestructor(self)
    !!{
    Destructor for the fixed critical overdensity class.
    !!}
    implicit none
    type(criticalOverdensityFixed), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%linearGrowth_"            />
    !!]
    return
  end subroutine fixedDestructor

  double precision function fixedValue(self,time,expansionFactor,collapsing,mass,node)
    !!{
    Return the critical overdensity at the given time and mass.
    !!}
    implicit none
    class           (criticalOverdensityFixed), intent(inout)           :: self
    double precision                          , intent(in   ), optional :: time      , expansionFactor
    logical                                   , intent(in   ), optional :: collapsing
    double precision                          , intent(in   ), optional :: mass
    type            (treeNode                ), intent(inout), optional :: node
    !$GLC attributes unused :: mass, node, time, expansionFactor, collapsing

    fixedValue=+self%criticalOverdensity_
    return
  end function fixedValue

  double precision function fixedGradientTime(self,time,expansionFactor,collapsing,mass,node)
    !!{
    Return the gradient with respect to time of critical overdensity at the given time and mass.
    !!}
    implicit none
    class           (criticalOverdensityFixed), intent(inout)           :: self
    double precision                          , intent(in   ), optional :: time      , expansionFactor
    logical                                   , intent(in   ), optional :: collapsing
    double precision                          , intent(in   ), optional :: mass
    type            (treeNode                ), intent(inout), optional :: node
    !$GLC attributes unused :: self, mass, node, time, expansionFactor, collapsing

    fixedGradientTime=0.0d0
    return
  end function fixedGradientTime

  double precision function fixedGradientMass(self,time,expansionFactor,collapsing,mass,node)
    !!{
    Return the gradient with respect to mass of critical overdensity at the given time and mass.
    !!}
    implicit none
    class           (criticalOverdensityFixed), intent(inout)           :: self
    double precision                          , intent(in   ), optional :: time      , expansionFactor
    logical                                   , intent(in   ), optional :: collapsing
    double precision                          , intent(in   ), optional :: mass
    type            (treeNode                ), intent(inout), optional :: node
    !$GLC attributes unused :: self, time, expansionFactor, collapsing, mass, node

    fixedGradientMass=0.0d0
    return
  end function fixedGradientMass

  logical function fixedIsMassDependent(self)
    !!{
    Return whether the critical overdensity is mass dependent.
    !!}
    implicit none
    class(criticalOverdensityFixed), intent(inout) :: self
    !$GLC attributes unused :: self

    fixedIsMassDependent=.false.
    return
  end function fixedIsMassDependent

  logical function fixedIsNodeDependent(self)
    !!{
    Return whether the critical overdensity is node dependent.
    !!}
    implicit none
    class(criticalOverdensityFixed), intent(inout) :: self
    !$GLC attributes unused :: self

    fixedIsNodeDependent=.false.
    return
  end function fixedIsNodeDependent

  logical function fixedIsTreeDependent(self)
    !!{
    Return whether the critical overdensity is tree dependent.
    !!}
    implicit none
    class(criticalOverdensityFixed), intent(inout) :: self
    !$GLC attributes unused :: self

    fixedIsTreeDependent=.false.
    return
  end function fixedIsTreeDependent
