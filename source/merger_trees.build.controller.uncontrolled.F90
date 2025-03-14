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
Implements a merger tree build controller class which provides no control.
!!}

  !![
  <mergerTreeBuildController name="mergerTreeBuildControllerUncontrolled">
   <description>A merger tree build controller class which provides no control.</description>
  </mergerTreeBuildController>
  !!]
  type, extends(mergerTreeBuildControllerClass) :: mergerTreeBuildControllerUncontrolled
     !!{
     A merger tree build controller class which provides no control.
     !!}
     private
     class(mergerTreeBranchingProbabilityClass), pointer :: mergerTreeBranchingProbability_ => null()
   contains
     final     ::                               uncontrolledDestructor
     procedure :: control                    => uncontrolledControl
     procedure :: branchingProbabilityObject => uncontrolledBranchingProbabilityObject
     procedure :: nodesInserted              => uncontrolledNodesInserted
  end type mergerTreeBuildControllerUncontrolled

  interface mergerTreeBuildControllerUncontrolled
     !!{
     Constructors for the ``uncontrolled'' merger tree build controller class.
     !!}
     module procedure uncontrolledConstructorParameters
     module procedure uncontrolledConstructorInternal
  end interface mergerTreeBuildControllerUncontrolled

contains

  function uncontrolledConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``uncontrolled'' merger tree build controller class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (mergerTreeBuildControllerUncontrolled)                :: self
    type (inputParameters                      ), intent(inout) :: parameters
    class(mergerTreeBranchingProbabilityClass  ), pointer       :: mergerTreeBranchingProbability_

    !![
    <objectBuilder class="mergerTreeBranchingProbability" name="mergerTreeBranchingProbability_" source="parameters"/>
    !!]
    self=mergerTreeBuildControllerUncontrolled(mergerTreeBranchingProbability_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="mergerTreeBranchingProbability_"/>
    !!]
    return
  end function uncontrolledConstructorParameters

  function uncontrolledConstructorInternal(mergerTreeBranchingProbability_) result(self)
    !!{
    Internal constructor for the ``uncontrolled'' merger tree build controller class .
    !!}
    implicit none
    type (mergerTreeBuildControllerUncontrolled)                     :: self
    class(mergerTreeBranchingProbabilityClass  ), intent(in), target :: mergerTreeBranchingProbability_
    !![
    <constructorAssign variables="*mergerTreeBranchingProbability_"/>
    !!]
    
    return
  end function uncontrolledConstructorInternal

  subroutine uncontrolledDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily uncontrolled} merger tree build controller class.
    !!}
    implicit none
    type(mergerTreeBuildControllerUncontrolled), intent(inout) :: self

    !![
    <objectDestructor name="self%mergerTreeBranchingProbability_"/>
    !!]
    return
  end subroutine uncontrolledDestructor

  logical function uncontrolledControl(self,node,treeWalker_)
    !!{
    Apply control to merger tree building.
    !!}
    implicit none
    class(mergerTreeBuildControllerUncontrolled), intent(inout)           :: self
    type (treeNode                             ), intent(inout), pointer  :: node
    class(mergerTreeWalkerClass                ), intent(inout), optional :: treeWalker_
    !$GLC attributes unused :: self, node, treeWalker_

    ! Always return true as we never want to halt tree building.
    uncontrolledControl=.true.
    return
  end function uncontrolledControl

  function uncontrolledBranchingProbabilityObject(self,node) result(mergerTreeBranchingProbability_)
    !!{
    Return a pointer the the merger tree branching probability object to use.
    !!}
    implicit none
    class(mergerTreeBranchingProbabilityClass  ), pointer       :: mergerTreeBranchingProbability_
    class(mergerTreeBuildControllerUncontrolled), intent(inout) :: self
    type (treeNode                             ), intent(inout) :: node
    !$GLC attributes unused :: node

    mergerTreeBranchingProbability_ => self%mergerTreeBranchingProbability_
    return
  end function uncontrolledBranchingProbabilityObject

  subroutine uncontrolledNodesInserted(self,nodeCurrent,nodeProgenitor1,nodeProgenitor2,didBranch)
    !!{
    Act on the insertion of nodes into the merger tree.
    !!}
    implicit none
    class  (mergerTreeBuildControllerUncontrolled), intent(inout)           :: self
    type   (treeNode                             ), intent(inout)           :: nodeCurrent     , nodeProgenitor1
    type   (treeNode                             ), intent(inout), optional :: nodeProgenitor2
    logical                                       , intent(in   ), optional :: didBranch
    !$GLC attributes unused :: self, nodeCurrent, nodeProgenitor1, nodeProgenitor2, didBranch

    ! Nothing to do.
    return
  end subroutine uncontrolledNodesInserted
