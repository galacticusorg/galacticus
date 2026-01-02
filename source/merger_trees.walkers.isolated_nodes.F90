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
Implements a depth-first merger tree walker over all isolated nodes.
!!}

  !![
  <mergerTreeWalker name="mergerTreeWalkerIsolatedNodes">
   <description>Provides a merger tree walker which iterates depth-first over all isolated nodes.</description>
  </mergerTreeWalker>
  !!]
  type, extends(mergerTreeWalkerAllNodes) :: mergerTreeWalkerIsolatedNodes
     !!{
     A merger tree walker which iterates depth-first over all isolated nodes.
     !!}
     private
   contains
     procedure :: next => isolatedNodesNext
  end type mergerTreeWalkerIsolatedNodes

  interface mergerTreeWalkerIsolatedNodes
     !!{
     Constructors for the \refClass{mergerTreeWalkerIsolatedNodes} merger tree walker class.
     !!}
     module procedure isolatedNodesParameters
     module procedure isolatedNodesInternal
  end interface mergerTreeWalkerIsolatedNodes

contains

  function isolatedNodesParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeWalkerIsolatedNodes} merger tree walker class which takes a parameter set as input.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(mergerTreeWalkerIsolatedNodes)                :: self
    type(inputParameters              ), intent(inout) :: parameters
    !$GLC attributes unused :: self, parameters

    call Error_Report('this class can not be built from parameters'//{introspection:location})
    return
  end function isolatedNodesParameters

  function isolatedNodesInternal(tree,spanForest) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeWalkerIsolatedNodes} merger tree walker class.
    !!}
    implicit none
    type(mergerTreeWalkerIsolatedNodes)                          :: self
    type(mergerTree                   ), intent(in   ), target   :: tree
    logical                            , intent(in   ), optional :: spanForest

    self%mergerTreeWalkerAllNodes=mergerTreeWalkerAllNodes(tree,spanForest)
    return
  end function isolatedNodesInternal

  logical function isolatedNodesNext(self,node)
    !!{
    This function will update to given {\normalfont \ttfamily node} to the next node which should be visited in a tree to
    perform a depth-first walk. Once the entire tree has been walked, a {\normalfont \ttfamily null()} pointer will be set, and
    a value of {\normalfont \ttfamily false} returned indicating that there are no more nodes to walk. Each node will be
    visited once and once only if the tree is walked in this way.
    !!}
    implicit none
    class(mergerTreeWalkerIsolatedNodes), intent(inout)          :: self
    type (treeNode                     ), intent(inout), pointer :: node

    ! If we have already processed all nodes, simply return false and a null pointer.
    if (.not.self%nodesRemain_) then
       node              => null()
       isolatedNodesNext =  .false.
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
             if (associated(self%tree).and.associated(self%tree%nodeBase)) exit
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
          isolatedNodesNext =  .false.
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
          isolatedNodesNext =  .false.
          return
       end if
    end if
    ! Walk to the next node in the tree.
    if (.not.associated(self%node%parent)) then
       ! This is the base of the merger tree.
       do while (associated(self%node%firstChild))
          self%node => self%node%firstChild
       end do
    else
       if (associated(self%node%sibling)) then
          self%node => self%node%sibling
          do while (associated(self%node%firstChild))
             self%node => self%node%firstChild
          end do
       else
          self%node => self%node%parent
       end if
    end if
    node              => self%node
    isolatedNodesNext =  .true.
    return
  end function isolatedNodesNext
