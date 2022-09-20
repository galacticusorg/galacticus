!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
Contains a module which implements a merger tree build controller class which builds constrained trees.
!!}

  !![
  <mergerTreeBuildController name="mergerTreeBuildControllerConstrained">
   <description>A merger tree build controller class which builds constrained trees.</description>
  </mergerTreeBuildController>
  !!]
  type, extends(mergerTreeBuildControllerClass) :: mergerTreeBuildControllerConstrained
     !!{     
     A merger tree build controller class which builds constrained trees.
     !!}
     private
     class  (mergerTreeBranchingProbabilityClass), pointer :: mergerTreeBranchingProbabilityUnconstrained_ => null(), mergerTreeBranchingProbabilityConstrained_ => null()
     integer                                               :: isConstrainedID
   contains
     final     ::                               constrainedDestructor
     procedure :: control                    => constrainedControl
     procedure :: branchingProbabilityObject => constrainedBranchingProbabilityObject
     procedure :: nodesInserted              => constrainedNodesInserted
  end type mergerTreeBuildControllerConstrained
  
  interface mergerTreeBuildControllerConstrained
     !!{
     Constructors for the ``constrained'' merger tree build controller class.
     !!}
     module procedure constrainedConstructorParameters
     module procedure constrainedConstructorInternal
  end interface mergerTreeBuildControllerConstrained

contains

  function constrainedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``constrained'' merger tree build controller class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (mergerTreeBuildControllerConstrained)               :: self
    type (inputParameters                    ), intent(inout) :: parameters
    class(mergerTreeBranchingProbabilityClass), pointer       :: mergerTreeBranchingProbabilityUnconstrained_, mergerTreeBranchingProbabilityConstrained_

    !![
    <objectBuilder class="mergerTreeBranchingProbability" name="mergerTreeBranchingProbabilityUnconstrained_" parameterName="mergerTreeBranchingProbabilityUnconstrained" source="parameters"/>
    <objectBuilder class="mergerTreeBranchingProbability" name="mergerTreeBranchingProbabilityConstrained_"   parameterName="mergerTreeBranchingProbabilityConstrained"   source="parameters"/>
    !!]
    self=mergerTreeBuildControllerConstrained(mergerTreeBranchingProbabilityUnconstrained_,mergerTreeBranchingProbabilityConstrained_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="mergerTreeBranchingProbabilityUnconstrained_"/>
    <objectDestructor name="mergerTreeBranchingProbabilityConstrained_"  />
    !!]
    return
  end function constrainedConstructorParameters

  function constrainedConstructorInternal(mergerTreeBranchingProbabilityUnconstrained_,mergerTreeBranchingProbabilityConstrained_) result(self)
    !!{
    Internal constructor for the ``constrained'' merger tree build controller class.
    !!}
    implicit none
    type (mergerTreeBuildControllerConstrained)                        :: self
    class(mergerTreeBranchingProbabilityClass ), intent(in   ), target :: mergerTreeBranchingProbabilityUnconstrained_, mergerTreeBranchingProbabilityConstrained_
    !![
    <constructorAssign variables="*mergerTreeBranchingProbabilityUnconstrained_, *mergerTreeBranchingProbabilityConstrained_"/>
    !!]

    !! NOTE: Here we add a "meta-property" to the "basic" component of each node. (The basic component is the thing that stores
    !! the total mass and current time of the node, so it's guaraneteed to already exist at this point. A meta-property is just
    !! some arbitrary data that we want to attach to each node in the tree. The directive below will obtain and store an ID value
    !! associated with this meta-property that we can then use to get and set it.) We'll use this meta-property to store the
    !! status of whether a node is on the constrained branch or not.
    !![
    <addMetaProperty component="basic" name="isConstrained" type="integer" id="self%isConstrainedID" isCreator="yes"/>
    !!]
    return
  end function constrainedConstructorInternal

  subroutine constrainedDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily constrained} merger tree build controller class.
    !!}
    implicit none
    type(mergerTreeBuildControllerConstrained), intent(inout) :: self

    !![
    <objectDestructor name="self%mergerTreeBranchingProbabilityUnconstrained_"/>
    <objectDestructor name="self%mergerTreeBranchingProbabilityConstrained_"  />
    !!]
    return
  end subroutine constrainedDestructor

  logical function constrainedControl(self,node,treeWalker_)
    !!{
    Apply control to merger tree building.
    !!}
    implicit none
    class(mergerTreeBuildControllerConstrained), intent(inout)          :: self
    type (treeNode                            ), intent(inout), pointer :: node
    class(mergerTreeWalkerClass               ), intent(inout)          :: treeWalker_
    !$GLC attributes unused :: self, node, treeWalker_

    ! Always return true as we never want to halt tree building.
    constrainedControl=.true.
    return
  end function constrainedControl

  function constrainedBranchingProbabilityObject(self,node) result(mergerTreeBranchingProbability_)
    !!{
    Return a pointer the the merger tree branching probability object to use.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class  (mergerTreeBranchingProbabilityClass ), pointer       :: mergerTreeBranchingProbability_
    class  (mergerTreeBuildControllerConstrained), intent(inout) :: self
    type   (treeNode                            ), intent(inout) :: node
    class  (nodeComponentBasic                  ), pointer       :: basic
    logical                                                      :: isConstrained

    basic         => node %basic                      (                   )
    isConstrained =  basic%integerRank0MetaPropertyGet(self%isConstrainedID) == 1
    !! NOTE: "isConstrained" will now be true if this node is on the constrained branch. So, we need an if/else here to return a
    !! pointer to the appropriate mergerTreeBranching object depending on the state of "isConstrained".  I think one further
    !! detail here is that we will need to also check if this node is the root node of the tree (which we can do with
    !! ".not.associated(node%parent)"). If it is, we also need to return the constrained branching probabilty, *and* mark it as
    !! being on the constrained branch.
    return
  end function constrainedBranchingProbabilityObject

  subroutine constrainedNodesInserted(self,nodeCurrent,nodeProgenitor1,nodeProgenitor2)
    !!{
    Act on the insertion of nodes into the merger tree.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class (mergerTreeBuildControllerConstrained), intent(inout)           :: self
    type  (treeNode                            ), intent(inout)           :: nodeCurrent     , nodeProgenitor1
    type  (treeNode                            ), intent(inout), optional :: nodeProgenitor2
    class (nodeComponentBasic                  ), pointer                 :: basicCurrent    , basicProgenitor1, &
         &                                                                   basicProgenitor2
    logical                                                               :: isConstrained

    !! NOTE: Here is where we should mark nodes as on the constrained branch or not.
    basicCurrent     => nodeCurrent    %basic                      (                    )
    basicProgenitor1 => nodeProgenitor1%basic                      (                    )
    isConstrained    =  basicCurrent   %integerRank0MetaPropertyGet(self%isConstrainedID) == 1
    if (isConstrained) then
       !! NOTE: Parent is on the constrained branch, so this progenitor also is - mark it as such using the "integerRank0MetaPropertySet" function.
    else
       !! NOTE: Parent is not on the constrained branch, so this progenitor also is not - mark it as such using the "integerRank0MetaPropertySet" function.
    end if
    ! If the second progenitor is present, mark it as not on the constrained branch.
    if (present(nodeProgenitor2)) then
       basicProgenitor2 => nodeProgenitor2%basic()
       !! NOTE: Need to mark this secondary progenitor as not on the main branch, e.g.:
       call basicProgenitor2%integerRank0MetaPropertySet(self%isConstrainedID,0)
    end if
    return
  end subroutine constrainedNodesInserted
