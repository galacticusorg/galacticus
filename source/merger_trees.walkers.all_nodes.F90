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
  Implements a depth-first merger tree walker over all all nodes.
  !!}
  use :: Galacticus_Nodes, only : mergerTree, treeNode

  !![
  <mergerTreeWalker name="mergerTreeWalkerAllNodes">
   <description>Provides a merger tree walker which iterates depth-first over all all nodes.</description>
  </mergerTreeWalker>
  !!]
  type, extends(mergerTreeWalkerClass) :: mergerTreeWalkerAllNodes
     !!{
     A merger tree walker which iterates depth-first over all all nodes.
     !!}
     private
     type   (mergerTree), pointer :: tree       => null(), treePrevious => null()
     type   (treeNode  ), pointer :: node       => null(), nodePrevious => null()
     logical                      :: spanForest          , nodesRemain_
   contains
     !![
     <methods>
       <method description="Step back to the previously visited node (if possible)."                     method="previous"/>
       <method description="Set the walker to the given node."                                           method="setNode" />
       <method description="Descend through the hierarchy to the deepest node along the current branch." method="descend" />
     </methods>
     !!]
     procedure :: next        => allNodesNext
     procedure :: previous    => allNodesPrevious
     procedure :: setNode     => allNodesSetNode
     procedure :: nodesRemain => allNodesNodesRemain
     procedure :: descend     => allNodesDescend
  end type mergerTreeWalkerAllNodes

  interface mergerTreeWalkerAllNodes
     !!{
     Constructors for the {\normalfont \ttfamily allNodes} merger tree walker class.
     !!}
     module procedure allNodesParameters
     module procedure allNodesInternal
  end interface mergerTreeWalkerAllNodes

contains

  function allNodesParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily allNodes} merger tree walker class which takes a parameter set as input.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(mergerTreeWalkerAllNodes)                :: self
    type(inputParameters         ), intent(inout) :: parameters
    !$GLC attributes unused :: self, parameters

    call Error_Report('this class can not be built from parameters'//{introspection:location})
    return
  end function allNodesParameters

  function allNodesInternal(tree,spanForest) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily allNodes} merger tree walker class.
    !!}
    implicit none
    type   (mergerTreeWalkerAllNodes)                          :: self
    type   (mergerTree              ), intent(in   ), target   :: tree
    logical                          , intent(in   ), optional :: spanForest
    !![
    <optionalArgument name="spanForest" defaultsTo=".false."/>
    !!]

    self%tree         => tree
    self%node         => null()
    self%treePrevious => null()
    self%nodePrevious => null()
    self%spanForest   =  spanForest_
    self%nodesRemain_ = .true.
    return
  end function allNodesInternal

  logical function allNodesNext(self,node)
    !!{
    This function will update to given {\normalfont \ttfamily node} to the next node which should be visited in a tree to
    perform a depth-first walk, including satellite nodes. Once the entire tree has been walked, a {\normalfont \ttfamily
    null()} pointer will be set, and a value of {\normalfont \ttfamily false} returned indicating that there are no more nodes
    to walk. Each node will be visited once and once only if the tree is walked in this way. Note that it is important that the
    walk descends to satellites before descending to children: the routines that destroy merger tree branches rely on this
    since child nodes are used in testing whether a node is a satellite---if they are destroyed prior to the test being made
    then problems with dangling pointers will occur.
    !!}
    implicit none
    class(mergerTreeWalkerAllNodes), intent(inout)          :: self
    type (treeNode                ), intent(inout), pointer :: node

    ! If we have already processed all nodes, simply return false and a null pointer.
    if (.not.self%nodesRemain_) then
       node         => null()
       allNodesNext =  .false.
       return
    end if
    ! Store the current node as the previous node.
    self%treePrevious => self%tree
    self%nodePrevious => self%node
    ! If the node is currently pointing to the base node of the tree, then attempt to move to the next tree (if we are spanning
    ! forests) - if this fails the tree walk is complete.
    if (associated(self%node,self%tree%nodeBase)) then
       if (self%spanForest) then
          do while (associated(self%tree))
             self%tree => self%tree%nextTree
             if (associated(self%tree)) then
                if (associated(self%tree%nodeBase)) exit
             end if
          end do
       else
          self%tree => null()
       end if
       if (associated(self%tree)) then
          ! Nullify the node which will cause the walk to begin on the new tree.
          self%node    => null()
       else
          ! No more trees exist to walk, so return a null node and a result of false.
          self%nodesRemain_ =  .false.
          node              => null()
          allNodesNext      =  .false.
          return
       end if
    end if
    ! If the node is currently null, set to the base node of the tree.
    if (.not.associated(self%node)) then
       if (associated(self%tree%nodeBase)) then
          self%node => self%tree%nodeBase
       else
          self%nodesRemain_ =  .false.
          node              => null()
          allNodesNext      =  .false.
          return
       end if
    end if
    ! Walk to the next node in the tree.
    if (associated(self%node,self%tree%nodeBase)) then
       ! This is the base of the merger tree.
       ! Descend through satellites and children.
       call self%descend()
    else
       if (associated(self%node%sibling)) then
          self%node => self%node%sibling
          call self%descend()
       else
          ! About to move back up the tree. Check if the node we're moving up from is a satellite.
          if (self%node%isSatellite()) then
             ! It is a satellite. Therefore, the parent may have children that have yet to be
             ! visited. Check if the parent has children.
             if (associated(self%node%parent%firstChild)) then
                ! Parent does have children, so move to the first one.
                self%node => self%node%parent%firstChild
                ! Descend through satellites and children.
                call self%descend()
             else
                ! Parent has no children, so move to the parent.
                self%node => self%node%parent
             end if
          else
             ! It is not a satellite, so all satellites and children have been processed.
             self%node => self%node%parent
          end if
       end if
    end if
    node         => self%node
    allNodesNext =  .true.
    return
  end function allNodesNext

  subroutine allNodesPrevious(self,node)
    !!{
    Step back to the previously visited node.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(mergerTreeWalkerAllNodes), intent(inout)          :: self
    type (treeNode                ), intent(inout), pointer :: node

    if (associated(self%treePrevious)) then
       self%tree         => self%treePrevious
       self%node         => self%nodePrevious
       node              => self%node
       self%treePrevious => null()
       self%nodePrevious => null()
    else
       call Error_Report('no previous node is available'//{introspection:location})
    end if
    return
  end subroutine allNodesPrevious

  subroutine allNodesDescend(self)
    !!{
    Descend to the deepest progenitor (satellites and children) of the current branch.
    !!}
    implicit none
    class(mergerTreeWalkerAllNodes), intent(inout) :: self

    ! Descend through satellites and children.
    do while (associated(self%node%firstSatellite).or.associated(self%node%firstChild))
       if (associated(self%node%firstSatellite)) then
          self%node => self%node%firstSatellite
       else
          self%node => self%node%firstChild
       end if
    end do
    return
  end subroutine allNodesDescend

  logical function allNodesNodesRemain(self)
    !!{
    Returns true if more nodes remain to be walked to.
    !!}
    implicit none
    class(mergerTreeWalkerAllNodes), intent(inout) :: self

    allNodesNodesRemain=self%nodesRemain_
    return
  end function allNodesNodesRemain

  subroutine allNodesSetNode(self,node)
    !!{
    Set the current node for the walker.
    !!}
    implicit none
    class(mergerTreeWalkerAllNodes), intent(inout)         :: self
    type (treeNode                ), intent(in   ), target :: node

    self%nodesRemain_ =  .true.
    self%node         => node
    self%tree         => node%hostTree
    self%treePrevious => null()
    self%nodePrevious => null()
    return
  end subroutine allNodesSetNode

