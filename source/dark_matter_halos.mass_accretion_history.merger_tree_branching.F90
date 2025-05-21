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
  An implementation of dark matter halo mass accretion histories computed from merger tree branching rates.
  !!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass      , criticalOverdensityClass
  use :: Merger_Tree_Branching     , only : mergerTreeBranchingProbabilityClass

  !![
  <darkMatterHaloMassAccretionHistory name="darkMatterHaloMassAccretionHistoryMergerTreeBranching">
   <description>Dark matter halo mass accretion histories computed from merger tree branching rates.</description>
  </darkMatterHaloMassAccretionHistory>
  !!]
  type, extends(darkMatterHaloMassAccretionHistoryClass) :: darkMatterHaloMassAccretionHistoryMergerTreeBranching
     !!{
     A dark matter halo mass accretion history class computed from merger tree branching rates.
     !!}
     private
     class(cosmologicalMassVarianceClass      ), pointer :: cosmologicalMassVariance_       => null()
     class(criticalOverdensityClass           ), pointer :: criticalOverdensity_            => null()
     class(mergerTreeBranchingProbabilityClass), pointer :: mergerTreeBranchingProbability_ => null()
   contains
     final     ::                      mergerTreeBranchingDestructor
     procedure :: time              => mergerTreeBranchingTime
     procedure :: massAccretionRate => mergerTreeBranchingMassAccretionRate
  end type darkMatterHaloMassAccretionHistoryMergerTreeBranching

  interface darkMatterHaloMassAccretionHistoryMergerTreeBranching
     !!{
     Constructors for the \refClass{darkMatterHaloMassAccretionHistoryMergerTreeBranching} dark matter halo mass accretion history class.
     !!}
     module procedure mergerTreeBranchingConstructorParameters
     module procedure mergerTreeBranchingConstructorInternal
  end interface darkMatterHaloMassAccretionHistoryMergerTreeBranching

contains

  function mergerTreeBranchingConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterHaloMassAccretionHistoryMergerTreeBranching} dark matter halo mass accretion history class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (darkMatterHaloMassAccretionHistoryMergerTreeBranching)                :: self
    type (inputParameters                                      ), intent(inout) :: parameters
    class(cosmologicalMassVarianceClass                        ), pointer       :: cosmologicalMassVariance_
    class(criticalOverdensityClass                             ), pointer       :: criticalOverdensity_
    class(mergerTreeBranchingProbabilityClass                  ), pointer       :: mergerTreeBranchingProbability_

    !![
    <objectBuilder class="criticalOverdensity"            name="criticalOverdensity_"            source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance"       name="cosmologicalMassVariance_"       source="parameters"/>
    <objectBuilder class="mergerTreeBranchingProbability" name="mergerTreeBranchingProbability_" source="parameters"/>
    !!]
    self=darkMatterHaloMassAccretionHistoryMergerTreeBranching(criticalOverdensity_,cosmologicalMassVariance_,mergerTreeBranchingProbability_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="criticalOverdensity_"           />
    <objectDestructor name="cosmologicalMassVariance_"      />
    <objectDestructor name="mergerTreeBranchingProbability_"/>
    !!]
    return
  end function mergerTreeBranchingConstructorParameters

  function mergerTreeBranchingConstructorInternal(criticalOverdensity_,cosmologicalMassVariance_,mergerTreeBranchingProbability_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterHaloMassAccretionHistoryMergerTreeBranching} dark matter halo mass accretion history class.
    !!}
    implicit none
    type (darkMatterHaloMassAccretionHistoryMergerTreeBranching)                        :: self
    class(mergerTreeBranchingProbabilityClass                  ), intent(in   ), target :: mergerTreeBranchingProbability_
    class(cosmologicalMassVarianceClass                        ), intent(in   ), target :: cosmologicalMassVariance_
    class(criticalOverdensityClass                             ), intent(in   ), target :: criticalOverdensity_
    !![
    <constructorAssign variables="*criticalOverdensity_, *cosmologicalMassVariance_, *mergerTreeBranchingProbability_"/>
    !!]

    return
  end function mergerTreeBranchingConstructorInternal

  subroutine mergerTreeBranchingDestructor(self)
    !!{
    Destructor for the \refClass{darkMatterHaloMassAccretionHistoryMergerTreeBranching} dark matter halo mass accretion history class.
    !!}
    implicit none
    type(darkMatterHaloMassAccretionHistoryMergerTreeBranching), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologicalMassVariance_"      />
    <objectDestructor name="self%criticalOverdensity_"           />
    <objectDestructor name="self%mergerTreeBranchingProbability_"/>
    !!]
    return
  end subroutine mergerTreeBranchingDestructor

  double precision function mergerTreeBranchingTime(self,node,mass)
    !!{
    Compute the time corresponding to {\normalfont \ttfamily mass} in the mass accretion history.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (darkMatterHaloMassAccretionHistoryMergerTreeBranching), intent(inout), target :: self
    type            (treeNode                                             ), intent(inout), target :: node
    double precision                                                       , intent(in   )         :: mass
    !$GLC attributes unused :: self, node, mass

    mergerTreeBranchingTime=0.0d0
    call Error_Report('"time" method is not supported'//{introspection:location})
    return
  end function mergerTreeBranchingTime

  double precision function mergerTreeBranchingMass(self,node,time)
    !!{
    Compute the mass corresponding to {\normalfont \ttfamily time} in the mass accretion history.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (darkMatterHaloMassAccretionHistoryMergerTreeBranching), intent(inout), target :: self
    type            (treeNode                                             ), intent(inout), target :: node
    double precision                                                       , intent(in   )         :: time
    !$GLC attributes unused :: self, node, time

    mergerTreeBranchingMass=0.0d0
    call Error_Report('"mass" method is not supported'//{introspection:location})
    return
  end function mergerTreeBranchingMass

  double precision function mergerTreeBranchingMassAccretionRate(self,node,time)
    !!{
    Compute the mass accretion rate at the given {\normalfont \ttfamily time} in the mass accretion history of
    {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (darkMatterHaloMassAccretionHistoryMergerTreeBranching), intent(inout) :: self
    type            (treeNode                                             ), intent(inout) :: node
    double precision                                                       , intent(in   ) :: time
    class           (nodeComponentBasic                                   ), pointer       :: basic
    double precision                                                                       :: overdensityCritical              , overdensityCriticalGrowthRate, &
         &                                                                                    rootVarianceLogarithmicGrowthRate, barrierEffectiveGrowthRate

    basic => node%basic()
    ! Compute critical overdensity and its growth rate for the host halo.
    overdensityCritical                 =+self%criticalOverdensity_     %value                              (time=basic%time(),mass=basic%mass(),node=node)
    overdensityCriticalGrowthRate       =+self%criticalOverdensity_     %gradientTime                       (time=basic%time(),mass=basic%mass(),node=node)
    ! Compute the growth rate of the mass variance. 
    rootVarianceLogarithmicGrowthRate   =+self%cosmologicalMassVariance_%rootVarianceLogarithmicGradientTime(time=basic%time(),mass=basic%mass()          )
    ! Compute absolute value of the rate of change of the effective barrier.
    barrierEffectiveGrowthRate          =+abs(                                           &
         &                                    +      overdensityCriticalGrowthRate       &
         &                                    -      overdensityCritical                 &
         &                                    *      rootVarianceLogarithmicGrowthRate   &
         &                                    /basic%time                             () &
         &                                   )
    ! Compute the mass growth rate.
    mergerTreeBranchingMassAccretionRate=+self %mergerTreeBranchingProbability_%fractionSubresolution(basic%mass(),overdensityCritical,basic%time(),0.5d0*basic%mass(),node) &
         &                               *basic%mass                                                 (                                                                     ) &
         &                               *      barrierEffectiveGrowthRate
    return
  end function mergerTreeBranchingMassAccretionRate

