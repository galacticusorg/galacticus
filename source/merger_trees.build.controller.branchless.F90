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
Implements a merger tree build controller class which builds branchless trees.
!!}

  !![
  <mergerTreeBuildController name="mergerTreeBuildControllerBranchless">
   <description>A merger tree build controller class which builds branchless trees.</description>
  </mergerTreeBuildController>
  !!]
  type, extends(mergerTreeBuildControllerClass) :: mergerTreeBuildControllerBranchless
     !!{     
     A merger tree build controller class which builds branchless trees. Specifically, whenever a branching occurs, only the
     primary progenitor is followed. This results in a tree which consists of the main branch, along with stubs of any branches
     off of the main branch, but those stubs do not grow full branches of their own.
     !!}
     private
     class(mergerTreeBranchingProbabilityClass), pointer :: mergerTreeBranchingProbability_ => null()
   contains
     final     ::                               branchlessDestructor
     procedure :: control                    => branchlessControl
     procedure :: branchingProbabilityObject => branchlessBranchingProbabilityObject
     procedure :: nodesInserted              => branchlessNodesInserted
  end type mergerTreeBuildControllerBranchless

  interface mergerTreeBuildControllerBranchless
     !!{
     Constructors for the \refClass{mergerTreeBuildControllerBranchless} merger tree build controller class.
     !!}
     module procedure branchlessConstructorParameters
     module procedure branchlessConstructorInternal
  end interface mergerTreeBuildControllerBranchless

contains

  function branchlessConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeBuildControllerBranchless} merger tree build controller class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (mergerTreeBuildControllerBranchless)                :: self
    type (inputParameters                    ), intent(inout) :: parameters
    class(mergerTreeBranchingProbabilityClass), pointer       :: mergerTreeBranchingProbability_

    !![
    <objectBuilder class="mergerTreeBranchingProbability" name="mergerTreeBranchingProbability_" source="parameters"/>
    !!]
    self=mergerTreeBuildControllerBranchless(mergerTreeBranchingProbability_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="mergerTreeBranchingProbability_"/>
    !!]
    return
  end function branchlessConstructorParameters

  function branchlessConstructorInternal(mergerTreeBranchingProbability_) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeBuildControllerBranchless} merger tree build controller class .
    !!}
    implicit none
    type (mergerTreeBuildControllerBranchless)                        :: self
    class(mergerTreeBranchingProbabilityClass), intent(in   ), target :: mergerTreeBranchingProbability_
    !![
    <constructorAssign variables="*mergerTreeBranchingProbability_"/>
    !!]
    
    return
  end function branchlessConstructorInternal

  subroutine branchlessDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeBuildControllerBranchless} merger tree build controller class.
    !!}
    implicit none
    type(mergerTreeBuildControllerBranchless), intent(inout) :: self

    !![
    <objectDestructor name="self%mergerTreeBranchingProbability_"/>
    !!]
    return
  end subroutine branchlessDestructor

  logical function branchlessControl(self,node,treeWalker_)
    !!{
    Skip side branches of a tree under construction.
    !!}
    implicit none
    class(mergerTreeBuildControllerBranchless), intent(inout)           :: self    
    type (treeNode                           ), intent(inout), pointer  :: node
    class(mergerTreeWalkerClass              ), intent(inout), optional :: treeWalker_
    !$GLC attributes unused :: self

    branchlessControl=.true.
    ! Move to the next node in the tree while such exists, and the current node is on a side branch.
    do while (branchlessControl.and.associated(node%parent).and..not.node%isPrimaryProgenitor())
       if (present(treeWalker_)) then
          branchlessControl=treeWalker_%next(node)
       else
          branchlessControl=.false.
       end if
    end do
    return
  end function branchlessControl

  function branchlessBranchingProbabilityObject(self,node) result(mergerTreeBranchingProbability_)
    !!{
    Return a pointer the the merger tree branching probability object to use.
    !!}
    implicit none
    class(mergerTreeBranchingProbabilityClass), pointer       :: mergerTreeBranchingProbability_
    class(mergerTreeBuildControllerBranchless), intent(inout) :: self
    type (treeNode                           ), intent(inout) :: node
    !$GLC attributes unused :: node

    mergerTreeBranchingProbability_ => self%mergerTreeBranchingProbability_
    return
  end function branchlessBranchingProbabilityObject

  subroutine branchlessNodesInserted(self,nodeCurrent,nodeProgenitor1,nodeProgenitor2,didBranch)
    !!{
    Act on the insertion of nodes into the merger tree.
    !!}
    implicit none
    class  (mergerTreeBuildControllerBranchless), intent(inout)           :: self
    type   (treeNode                           ), intent(inout)           :: nodeCurrent    , nodeProgenitor1
    type   (treeNode                           ), intent(inout), optional :: nodeProgenitor2
    logical                                     , intent(in   ), optional :: didBranch
    !$GLC attributes unused :: self, nodeCurrent, nodeProgenitor1, nodeProgenitor2, didBranch

    ! Nothing to do.
    return
  end subroutine branchlessNodesInserted
