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
Implements a merger tree build controller class which follows branches only if they pass a filter.
!!}
  
  use :: Galactic_Filters, only : galacticFilterClass

  !![
  <mergerTreeBuildController name="mergerTreeBuildControllerFiltered">
   <description>A merger tree build controller class which builds filtered trees.</description>
  </mergerTreeBuildController>
  !!]
  type, extends(mergerTreeBuildControllerClass) :: mergerTreeBuildControllerFiltered
     !!{     
     A merger tree build controller class which builds filtered trees. Specifically, only nodes which pass a provided filtered are
     followed.
     !!}
     private
     class(mergerTreeBranchingProbabilityClass), pointer :: mergerTreeBranchingProbability_ => null()
     class(galacticFilterClass                ), pointer :: galacticFilter_                 => null()
   contains
     final     ::                               filteredDestructor
     procedure :: control                    => filteredControl
     procedure :: branchingProbabilityObject => filteredBranchingProbabilityObject
     procedure :: nodesInserted              => filteredNodesInserted
  end type mergerTreeBuildControllerFiltered

  interface mergerTreeBuildControllerFiltered
     !!{
     Constructors for the {\normalfont \ttfamily filtered} merger tree build controller class.
     !!}
     module procedure filteredConstructorParameters
     module procedure filteredConstructorInternal
  end interface mergerTreeBuildControllerFiltered

contains

  function filteredConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily filtered} merger tree build controller class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (mergerTreeBuildControllerFiltered  )                :: self
    type (inputParameters                    ), intent(inout) :: parameters
    class(mergerTreeBranchingProbabilityClass), pointer       :: mergerTreeBranchingProbability_
    class(galacticFilterClass                ), pointer       :: galacticFilter_

    !![
    <objectBuilder class="mergerTreeBranchingProbability" name="mergerTreeBranchingProbability_" source="parameters"/>
    <objectBuilder class="galacticFilter"                 name="galacticFilter_"                 source="parameters"/>
    !!]
    self=mergerTreeBuildControllerFiltered(mergerTreeBranchingProbability_,galacticFilter_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="mergerTreeBranchingProbability_"/>
    !!]
    return
  end function filteredConstructorParameters

  function filteredConstructorInternal(mergerTreeBranchingProbability_,galacticFilter_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily filtered} merger tree build controller class .
    !!}
    implicit none
    type (mergerTreeBuildControllerFiltered  )                        :: self
    class(mergerTreeBranchingProbabilityClass), intent(in   ), target :: mergerTreeBranchingProbability_
    class(galacticFilterClass                ), intent(in   ), target :: galacticFilter_
    !![
    <constructorAssign variables="*mergerTreeBranchingProbability_, *galacticFilter_"/>
    !!]
    
    return
  end function filteredConstructorInternal

  subroutine filteredDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily filtered} merger tree build controller class.
    !!}
    implicit none
    type(mergerTreeBuildControllerFiltered), intent(inout) :: self

    !![
    <objectDestructor name="self%mergerTreeBranchingProbability_"/>
    <objectDestructor name="self%galacticFilter_"                />
    !!]
    return
  end subroutine filteredDestructor

  logical function filteredControl(self,node,treeWalker_)
    !!{
    Skip side branches of a tree under construction.
    !!}
    implicit none
    class(mergerTreeBuildControllerFiltered), intent(inout)           :: self    
    type (treeNode                         ), intent(inout), pointer  :: node
    class(mergerTreeWalkerClass            ), intent(inout), optional :: treeWalker_
    !$GLC attributes unused :: self

    filteredControl=.true.
    ! Move to the next node in the tree while such exists, and the current node does not pass the filter.
    do while (filteredControl.and..not.self%galacticFilter_%passes(node))
       if (present(treeWalker_)) then
          filteredControl=treeWalker_%next(node)
       else
          filteredControl=.false.
       end if
    end do
    return
  end function filteredControl

  function filteredBranchingProbabilityObject(self,node) result(mergerTreeBranchingProbability_)
    !!{
    Return a pointer the the merger tree branching probability object to use.
    !!}
    implicit none
    class(mergerTreeBranchingProbabilityClass), pointer       :: mergerTreeBranchingProbability_
    class(mergerTreeBuildControllerFiltered  ), intent(inout) :: self
    type (treeNode                           ), intent(inout) :: node
    !$GLC attributes unused :: node

    mergerTreeBranchingProbability_ => self%mergerTreeBranchingProbability_
    return
  end function filteredBranchingProbabilityObject

  subroutine filteredNodesInserted(self,nodeCurrent,nodeProgenitor1,nodeProgenitor2,didBranch)
    !!{
    Act on the insertion of nodes into the merger tree.
    !!}
    implicit none
    class  (mergerTreeBuildControllerFiltered), intent(inout)           :: self
    type   (treeNode                         ), intent(inout)           :: nodeCurrent    , nodeProgenitor1
    type   (treeNode                         ), intent(inout), optional :: nodeProgenitor2
    logical                                   , intent(in   ), optional :: didBranch
    !$GLC attributes unused :: self, nodeCurrent, nodeProgenitor1, nodeProgenitor2, didBranch

    ! Nothing to do.
    return
  end subroutine filteredNodesInserted
