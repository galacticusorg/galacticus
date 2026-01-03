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
Implements an peak-background split critical overdensity class.
!!}

  !![
  <criticalOverdensity name="criticalOverdensityPeakBackgroundSplit">
   <description>The critical overdensity is given by some other critical overdensity class offset by the halo environmental overdensity.</description>
  </criticalOverdensity>
  !!]
  type, extends(criticalOverdensityClass) :: criticalOverdensityPeakBackgroundSplit
     !!{
     A peak-background split critical overdensity class.
     !!}
     private
     class(criticalOverdensityClass), pointer :: criticalOverdensity_ => null()
     class(haloEnvironmentClass    ), pointer :: haloEnvironment_     => null()
    contains
     final     ::                    peakBackgroundSplitDestructor
     procedure :: value           => peakBackgroundSplitValue
     procedure :: gradientTime    => peakBackgroundSplitGradientTime
     procedure :: gradientMass    => peakBackgroundSplitGradientMass
     procedure :: isMassDependent => peakBackgroundSplitIsMassDependent
     procedure :: isNodeDependent => peakBackgroundSplitIsNodeDependent
     procedure :: isTreeDependent => peakBackgroundSplitIsTreeDependent
  end type criticalOverdensityPeakBackgroundSplit

  interface criticalOverdensityPeakBackgroundSplit
     !!{
     Constructors for the \refClass{criticalOverdensityPeakBackgroundSplit} critical overdensity class.
     !!}
     module procedure peakBackgroundSplitConstructorParameters
     module procedure peakBackgroundSplitConstructorInternal
  end interface criticalOverdensityPeakBackgroundSplit

contains

  function peakBackgroundSplitConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{criticalOverdensityPeakBackgroundSplit} critical overdensity class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (criticalOverdensityPeakBackgroundSplit)                :: self
    type (inputParameters                       ), intent(inout) :: parameters
    class(criticalOverdensityClass              ), pointer       :: criticalOverdensity_
    class(haloEnvironmentClass                  ), pointer       :: haloEnvironment_
    class(cosmologyFunctionsClass               ), pointer       :: cosmologyFunctions_
    class(cosmologicalMassVarianceClass         ), pointer       :: cosmologicalMassVariance_
    class(linearGrowthClass                     ), pointer       :: linearGrowth_

    ! Check and read parameters.
    !![
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="haloEnvironment"          name="haloEnvironment_"          source="parameters"/>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="linearGrowth"             name="linearGrowth_"             source="parameters"/>
    !!]
    self=criticalOverdensityPeakBackgroundSplit(criticalOverdensity_,haloEnvironment_,cosmologyFunctions_,cosmologicalMassVariance_,linearGrowth_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="haloEnvironment_"         />
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="linearGrowth_"            />
    !!]
    return
  end function peakBackgroundSplitConstructorParameters

  function peakBackgroundSplitConstructorInternal(criticalOverdensity_,haloEnvironment_,cosmologyFunctions_,cosmologicalMassVariance_,linearGrowth_) result(self)
    !!{
    Internal constructor for the \refClass{criticalOverdensityPeakBackgroundSplit} critical overdensity class.
    !!}
    implicit none
    type (criticalOverdensityPeakBackgroundSplit)                        :: self
    class(criticalOverdensityClass              ), target, intent(in   ) :: criticalOverdensity_
    class(haloEnvironmentClass                  ), target, intent(in   ) :: haloEnvironment_
    class(cosmologyFunctionsClass               ), target, intent(in   ) :: cosmologyFunctions_
    class(cosmologicalMassVarianceClass         ), target, intent(in   ) :: cosmologicalMassVariance_
    class(linearGrowthClass                     ), target, intent(in   ) :: linearGrowth_
    !![
    <constructorAssign variables="*criticalOverdensity_, *haloEnvironment_, *cosmologyFunctions_, *cosmologicalMassVariance_, *linearGrowth_"/>
    !!]

    return
  end function peakBackgroundSplitConstructorInternal

  subroutine peakBackgroundSplitDestructor(self)
    !!{
    Destructor for the \refClass{criticalOverdensityPeakBackgroundSplit} critical overdensity class.
    !!}
    implicit none
    type(criticalOverdensityPeakBackgroundSplit), intent(inout) :: self

    !![
    <objectDestructor name="self%criticalOverdensity_"     />
    <objectDestructor name="self%haloEnvironment_"         />
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    <objectDestructor name="self%linearGrowth_"            />
    !!]
    return
  end subroutine peakBackgroundSplitDestructor

  double precision function peakBackgroundSplitValue(self,time,expansionFactor,collapsing,mass,node)
    !!{
    Return the critical overdensity for collapse at the given time and mass.
    !!}
    implicit none
    class           (criticalOverdensityPeakBackgroundSplit), intent(inout)           :: self
    double precision                                        , intent(in   ), optional :: time      , expansionFactor
    logical                                                 , intent(in   ), optional :: collapsing
    double precision                                        , intent(in   ), optional :: mass
    type            (treeNode                              ), intent(inout), optional :: node

    ! Get the critical overdensity at zero environmental overdensity.
    peakBackgroundSplitValue=+self%criticalOverdensity_%value(time,expansionFactor,collapsing,mass,node)
    ! Offset the critical overdensity by the environmental overdensity. Note that the convention used with Galacticus is that the
    ! critical overdensity always applies to perturbations at the current epoch, and σ(M,t) grows with time. The environmental
    ! overdensity function provides the overdensity at the current epoch.
    if (present(node))                                                                &
         & peakBackgroundSplitValue=+peakBackgroundSplitValue                         &
         &                          -self%haloEnvironment_   %overdensityLinear(node)
    return
  end function peakBackgroundSplitValue

  double precision function peakBackgroundSplitGradientTime(self,time,expansionFactor,collapsing,mass,node)
    !!{
    Return the gradient with respect to time of critical overdensity at the given time and mass.
    !!}
    implicit none
    class           (criticalOverdensityPeakBackgroundSplit), intent(inout)           :: self
    double precision                                        , intent(in   ), optional :: time      , expansionFactor
    logical                                                 , intent(in   ), optional :: collapsing
    double precision                                        , intent(in   ), optional :: mass
    type            (treeNode                              ), intent(inout), optional :: node

    peakBackgroundSplitGradientTime=+self%criticalOverdensity_%gradientTime                 (time,expansionFactor,collapsing,mass,node) &
         &                          -self%haloEnvironment_    %overdensityLinearGradientTime(                                     node)
    return
  end function peakBackgroundSplitGradientTime

  double precision function peakBackgroundSplitGradientMass(self,time,expansionFactor,collapsing,mass,node)
    !!{
    Return the gradient with respect to mass of critical overdensity at the given time and mass.
    !!}
    implicit none
    class           (criticalOverdensityPeakBackgroundSplit), intent(inout)           :: self
    double precision                                        , intent(in   ), optional :: time      , expansionFactor
    logical                                                 , intent(in   ), optional :: collapsing
    double precision                                        , intent(in   ), optional :: mass
    type            (treeNode                              ), intent(inout), optional :: node

    peakBackgroundSplitGradientMass=self%criticalOverdensity_%gradientMass(time,expansionFactor,collapsing,mass,node)
    return
  end function peakBackgroundSplitGradientMass

  logical function peakBackgroundSplitIsMassDependent(self)
    !!{
    Return whether the critical overdensity is mass dependent.
    !!}
    implicit none
    class(criticalOverdensityPeakBackgroundSplit), intent(inout) :: self

    peakBackgroundSplitIsMassDependent=self%criticalOverdensity_%isMassDependent()
    return
  end function peakBackgroundSplitIsMassDependent

  logical function peakBackgroundSplitIsNodeDependent(self)
    !!{
    Return whether the critical overdensity is node dependent.
    !!}
    implicit none
    class(criticalOverdensityPeakBackgroundSplit), intent(inout) :: self

    peakBackgroundSplitIsNodeDependent=self%haloEnvironment_%isNodeDependent()
    return
  end function peakBackgroundSplitIsNodeDependent

  logical function peakBackgroundSplitIsTreeDependent(self)
    !!{
    Return whether the critical overdensity is tree dependent.
    !!}
    implicit none
    class(criticalOverdensityPeakBackgroundSplit), intent(inout) :: self
    !$GLC attributes unused :: self

    peakBackgroundSplitIsTreeDependent=self%haloEnvironment_%isTreeDependent()
    return
  end function peakBackgroundSplitIsTreeDependent
