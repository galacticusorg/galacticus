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
Implements a depth-first merger tree walker over all nodes including formation nodes.
!!}

  !![
  <mergerTreeWalker name="mergerTreeWalkerAllAndFormationNodes">
   <description>Provides a merger tree walker which iterates depth-first over all nodes including formation nodes.</description>
  </mergerTreeWalker>
  !!]
  type, extends(mergerTreeWalkerAllNodes) :: mergerTreeWalkerAllAndFormationNodes
     !!{
     A merger tree walker which iterates depth-first over all nodes including formation nodes.
     !!}
     private
     type(treeNode), pointer :: nodeNonFormation => null()
   contains
     procedure :: next     => allAndFormationNodesNext
     procedure :: previous => allAndFormationNodesPrevious
  end type mergerTreeWalkerAllAndFormationNodes

  interface mergerTreeWalkerAllAndFormationNodes
     !!{
     Constructors for the \refClass{mergerTreeWalkerAllAndFormationNodes} merger tree walker class.
     !!}
     module procedure allAndFormationNodesParameters
     module procedure allAndFormationNodesInternal
  end interface mergerTreeWalkerAllAndFormationNodes

contains

  function allAndFormationNodesParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeWalkerAllAndFormationNodes} merger tree walker class which takes a parameter set as input.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(mergerTreeWalkerAllAndFormationNodes)                :: self
    type(inputParameters                     ), intent(inout) :: parameters
    !$GLC attributes unused :: self, parameters

    call Error_Report('this class can not be built from parameters'//{introspection:location})
    return
  end function allAndFormationNodesParameters

  function allAndFormationNodesInternal(tree,spanForest) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeWalkerAllAndFormationNodes} merger tree walker class.
    !!}
    implicit none
    type   (mergerTreeWalkerAllAndFormationNodes)                          :: self
    type   (mergerTree                          ), intent(in   ), target   :: tree
    logical                                      , intent(in   ), optional :: spanForest

    self%mergerTreeWalkerAllNodes=mergerTreeWalkerAllNodes(tree,spanForest)
    return
  end function allAndFormationNodesInternal

  logical function allAndFormationNodesNext(self,node)
    !!{
    Walk nodes of a tree including formation nodes.
    !!}
    implicit none
    class(mergerTreeWalkerAllAndFormationNodes), intent(inout)          :: self
    type (treeNode                            ), intent(inout), pointer :: node

    ! If we have already processed all nodes, simply return false and a null pointer.
    if (.not.self%nodesRemain_) then
       node                     => null()
       allAndFormationNodesNext =  .false.
       return
    end if
    ! Check for a formation node.
    if (associated(self%node).and.associated(self%node%formationNode)) then
       if (.not.associated(self%nodeNonFormation)) self%nodeNonFormation => self%node
       self%node                => self%node%formationNode
       node                     => self%node
       allAndFormationNodesNext =  .true.
    else
       if (associated(self%nodeNonFormation)) then
          self%node             => self%nodeNonFormation
          self%nodeNonFormation => null()
       end if
       allAndFormationNodesNext=self%mergerTreeWalkerAllNodes%next(node)
    end if
    return
  end function allAndFormationNodesNext

  subroutine allAndFormationNodesPrevious(self,node)
    !!{
    Step back to the previously visited node.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(mergerTreeWalkerAllAndFormationNodes), intent(inout)          :: self
    type (treeNode                            ), intent(inout), pointer :: node
    !$GLC attributes unused :: self, node

    call Error_Report('returning to previous node is not supported'//{introspection:location})
    return
  end subroutine allAndFormationNodesPrevious
