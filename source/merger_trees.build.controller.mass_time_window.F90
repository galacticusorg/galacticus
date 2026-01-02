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
Implements a merger tree build controller class which follows branches only if they lie within a window of time and mass.
!!}

  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Cosmological_Density_Field, only : criticalOverdensityClass

  !![
  <mergerTreeBuildController name="mergerTreeBuildControllerMassTimeWindow">
   <description>A merger tree build controller class which follows branches only if they lie within a window of time and mass.</description>
  </mergerTreeBuildController>
  !!]
  type, extends(mergerTreeBuildControllerClass) :: mergerTreeBuildControllerMassTimeWindow
     !!{     
     A merger tree build controller class which follows branches only if they lie within a window of time and mass.
     !!}
     private
     class(cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
     class(criticalOverdensityClass), pointer :: criticalOverdensity_ => null()
     class(mergerTreeBranchingProbabilityClass), pointer :: mergerTreeBranchingProbability_ => null()
     double precision :: timeMinimum_, redshiftMaximum, &
          & massMinimum
   contains
     !![
     <methods>
       <method method="passes" description="Returns true if the given node lies within the allowed window."/>
     </methods>
     !!]
     final     ::                               massTimeWindowDestructor
     procedure :: control                    => massTimeWindowControl
     procedure :: branchingProbabilityObject => massTimeWindowBranchingProbabilityObject
     procedure :: nodesInserted              => massTimeWindowNodesInserted
     procedure :: passes                     => massTimeWindowPasses
  end type mergerTreeBuildControllerMassTimeWindow

  interface mergerTreeBuildControllerMassTimeWindow
     !!{
     Constructors for the \refClass{mergerTreeBuildControllerMassTimeWindow} merger tree build controller class.
     !!}
     module procedure massTimeWindowConstructorParameters
     module procedure massTimeWindowConstructorInternal
  end interface mergerTreeBuildControllerMassTimeWindow

contains

  function massTimeWindowConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeBuildControllerMassTimeWindow} merger tree build controller class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    use :: Error, only : Error_Report
    implicit none
    type (mergerTreeBuildControllerMassTimeWindow  )                :: self
    type (inputParameters                    ), intent(inout) :: parameters
     class(cosmologyFunctionsClass), pointer :: cosmologyFunctions_
     class(criticalOverdensityClass), pointer :: criticalOverdensity_
    class(mergerTreeBranchingProbabilityClass), pointer       :: mergerTreeBranchingProbability_
    double precision :: timeMinimum, redshiftMaximum, massMinimum
    
    !![
    <inputParameter>
      <name>massMinimum</name>
      <source>parameters</source>
      <description>The minimum mass to which branches should be followed.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    <objectBuilder class="criticalOverdensity" name="criticalOverdensity_" source="parameters"/>
    <objectBuilder class="mergerTreeBranchingProbability" name="mergerTreeBranchingProbability_" source="parameters"/>
    !!]
    if (parameters%isPresent('timeMinimum')) then
       if (parameters%isPresent('redshiftMaximum')) call Error_Report('specify only one of [timeMinimum] and [redshiftMaximum]'//{introspection:location})
       !![
       <inputParameter>
	 <name>timeMinimum</name>
	 <source>parameters</source>
	 <description>The minimum time to which branches should be followed.</description>
       </inputParameter>
       !!]
    else if (parameters%isPresent('redshiftMaximum')) then
       !![
       <inputParameter>
	 <name>redshiftMaximum</name>
	 <source>parameters</source>
	 <description>The maximum redshift to which branches should be followed.</description>
       </inputParameter>
       !!]
       timeMinimum=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftMaximum))
     else
       call Error_Report('specify one of [timeMinimum] and [redshiftMaximum]'//{introspection:location})
    end if
    self=mergerTreeBuildControllerMassTimeWindow(timeMinimum,massMinimum,cosmologyFunctions_,mergerTreeBranchingProbability_,criticalOverdensity_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="mergerTreeBranchingProbability_"/>
    <objectDestructor name="cosmologyFunctions_"/>
    <objectDestructor name="criticalOverdensity_"/>
    !!]
    return
  end function massTimeWindowConstructorParameters

  function massTimeWindowConstructorInternal(timeMinimum_,massMinimum,cosmologyFunctions_,mergerTreeBranchingProbability_,criticalOverdensity_) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeBuildControllerMassTimeWindow} merger tree build controller class .
    !!}
    implicit none
    type (mergerTreeBuildControllerMassTimeWindow  )                        :: self
    class(mergerTreeBranchingProbabilityClass), intent(in   ), target :: mergerTreeBranchingProbability_
    class(cosmologyFunctionsClass                ), intent(in   ), target :: cosmologyFunctions_
    class(criticalOverdensityClass                ), intent(in   ), target :: criticalOverdensity_
    double precision, intent(in   ) :: timeMinimum_, massMinimum
    !![
    <constructorAssign variables="timeMinimum_,massMinimum,*mergerTreeBranchingProbability_, *cosmologyFunctions_, *criticalOverdensity_"/>
    !!]

    self%redshiftMaximum=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(timeMinimum_))
    return
  end function massTimeWindowConstructorInternal

  subroutine massTimeWindowDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeBuildControllerMassTimeWindow} merger tree build controller class.
    !!}
    implicit none
    type(mergerTreeBuildControllerMassTimeWindow), intent(inout) :: self

    !![
    <objectDestructor name="self%mergerTreeBranchingProbability_"/>
    <objectDestructor name="self%cosmologyFunctions_"                />
    <objectDestructor name="self%criticalOverdensity_"                />
    !!]
    return
  end subroutine massTimeWindowDestructor

  logical function massTimeWindowControl(self,node,treeWalker_)
    !!{
    Skip side branches of a tree under construction.
    !!}
    implicit none
    class(mergerTreeBuildControllerMassTimeWindow), intent(inout)           :: self    
    type (treeNode                         ), intent(inout), pointer  :: node
    class(mergerTreeWalkerClass            ), intent(inout), optional :: treeWalker_

    massTimeWindowControl=.true.
    ! Move to the next node in the tree while such exists, and the current node does not pass the filter.
    do while (massTimeWindowControl.and..not.self%passes(node))
       if (present(treeWalker_)) then
          massTimeWindowControl=treeWalker_%next(node)
       else
          massTimeWindowControl=.false.
       end if
    end do
    return
  end function massTimeWindowControl

  logical function massTimeWindowPasses(self,node) result(passes)
    !!{
    Return true if the given node lies within the allowed window.
    !!}
    use :: Galacticus_Nodes, only : treeNode, nodeComponentBasic
    implicit none
    class           (mergerTreeBuildControllerMassTimeWindow), intent(inout)          :: self    
    type            (treeNode                               ), intent(inout), pointer :: node
    class           (nodeComponentBasic                     )               , pointer :: basic
    double precision                                                                  :: time

    basic  => node                     %basic         (                                                            )
    time   =  self%criticalOverdensity_%timeOfCollapse(criticalOverdensity=basic%time(),mass=basic%mass(),node=node)
    passes =  time         >= self%timeMinimum_ &
         &   .and.                              &
         &    basic%mass() >= self%massMinimum
    return
  end function massTimeWindowPasses

  function massTimeWindowBranchingProbabilityObject(self,node) result(mergerTreeBranchingProbability_)
    !!{
    Return a pointer the the merger tree branching probability object to use.
    !!}
    implicit none
    class(mergerTreeBranchingProbabilityClass), pointer       :: mergerTreeBranchingProbability_
    class(mergerTreeBuildControllerMassTimeWindow  ), intent(inout) :: self
    type (treeNode                           ), intent(inout) :: node
    !$GLC attributes unused :: node

    mergerTreeBranchingProbability_ => self%mergerTreeBranchingProbability_
    return
  end function massTimeWindowBranchingProbabilityObject

  subroutine massTimeWindowNodesInserted(self,nodeCurrent,nodeProgenitor1,nodeProgenitor2,didBranch)
    !!{
    Act on the insertion of nodes into the merger tree.
    !!}
    implicit none
    class  (mergerTreeBuildControllerMassTimeWindow), intent(inout)           :: self
    type   (treeNode                         ), intent(inout)           :: nodeCurrent    , nodeProgenitor1
    type   (treeNode                         ), intent(inout), optional :: nodeProgenitor2
    logical                                   , intent(in   ), optional :: didBranch
    !$GLC attributes unused :: self, nodeCurrent, nodeProgenitor1, nodeProgenitor2, didBranch

    ! Nothing to do.
    return
  end subroutine massTimeWindowNodesInserted
