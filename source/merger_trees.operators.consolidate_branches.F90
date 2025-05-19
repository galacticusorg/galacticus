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
  Implements a merger tree operator that consolidates branches spanning a given amount of time or mass growth.
  !!}

  !![
  <mergerTreeOperator name="mergerTreeOperatorConsolidateBranches">
   <description>
    A merger tree operator class that consolidates branches spanning a given amount of time or mass growth. Starting from the tip
    of each branch, the branch is broken into segments for which the mass growth is less than $1+${\normalfont \ttfamily
    [fractionGrowthMass]} and the time growth is less than $1+${\normalfont \ttfamily [fractionGrowthTime]}. Any intermediate
    nodes in each segment are removed, with their siblings (if any) being made siblings of the node at the end of the
    segment. This reduces the time resolution along branches which can make evolution more efficient (at the cost of some loss of
    precision.
   </description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorConsolidateBranches
     !!{
     A merger tree operator class that consolidates branches spanning a given amount of time or mass growth.
     !!}
     private
     double precision :: fractionGrowthMass, fractionGrowthTime
   contains
     procedure :: operatePreEvolution => consolidateBranchesOperatePreEvolution
  end type mergerTreeOperatorConsolidateBranches

  interface mergerTreeOperatorConsolidateBranches
     !!{
     Constructors for the \refClass{mergerTreeOperatorConsolidateBranches} merger tree operator class.
     !!}
     module procedure consolidateBranchesConstructorParameters
     module procedure consolidateBranchesConstructorInternal
  end interface mergerTreeOperatorConsolidateBranches

contains

  function consolidateBranchesConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeOperatorConsolidateBranches} merger tree operator class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (mergerTreeOperatorConsolidateBranches)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    double precision                                                       :: fractionGrowthMass, fractionGrowthTime

    !![
    <inputParameter>
      <name>fractionGrowthMass</name>
      <description>The fraction of growth in mass over which branches may be consolidated.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>fractionGrowthTime</name>
      <description>The fraction of growth in time over which branches may be consolidated.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=mergerTreeOperatorConsolidateBranches(fractionGrowthMass,fractionGrowthTime)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function consolidateBranchesConstructorParameters

  function consolidateBranchesConstructorInternal(fractionGrowthMass,fractionGrowthTime) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeOperatorConsolidateBranches} merger tree operator class.
    !!}
    implicit none
    type            (mergerTreeOperatorConsolidateBranches)                :: self
    double precision                                       , intent(in   ) :: fractionGrowthMass, fractionGrowthTime
    !![
    <constructorAssign variables="fractionGrowthMass, fractionGrowthTime"/>
    !!]

    return
  end function consolidateBranchesConstructorInternal

  subroutine consolidateBranchesOperatePreEvolution(self,tree)
    !!{
    Perform a mass growth monotonizing operation on a merger tree.
    !!}
    use :: Galacticus_Nodes   , only : mergerTree                   , nodeComponentBasic, treeNode
    use :: Merger_Tree_Walkers, only : mergerTreeWalkerIsolatedNodes
    implicit none
    class           (mergerTreeOperatorConsolidateBranches), intent(inout), target :: self
    type            (mergerTree                           ), intent(inout), target :: tree
    type            (treeNode                             ), pointer               :: nodeDescendant  , nodeBase   , &
         &                                                                            node            , nodeNext   , &
         &                                                                            nodePriorSibling, nodeSibling, &
         &                                                                            nodeBranch
    class           (nodeComponentBasic                   ), pointer               :: basicDescendant , basic
    type            (mergerTree                           ), pointer               :: treeCurrent
    type            (mergerTreeWalkerIsolatedNodes        )                        :: treeWalker
    integer                                                                        :: countNodes

    ! Iterate over trees.
    treeCurrent => tree
    do while (associated(treeCurrent))       
       ! Walk the tree.
       treeWalker=mergerTreeWalkerIsolatedNodes(treeCurrent)
       do while (treeWalker%next(node))
          ! Ignore nodes that have children.
          if (associated(node%firstChild)) cycle
          ! Walk up the branch until sufficient growth in mass or time has occurred.
          nodeBase       => node
          nodeDescendant => node
          basic          => nodeBase%basic()
          countNodes     =  0
          do while (associated(nodeDescendant))
             basicDescendant => nodeDescendant%basic()
             countNodes      =  countNodes+1
             if     (                                                                       &
                  &   basicDescendant%mass() > (1.0d0+self%fractionGrowthMass)*basic%mass() &
                  &  .or.                                                                   &
                  &   basicDescendant%time() > (1.0d0+self%fractionGrowthTime)*basic%time() &
                  & ) then
                ! Back up to the previous node as we don't want to exceed the specified growth fractions.
                nodeDescendant => nodeDescendant%firstChild
                countNodes     =  countNodes-1
                ! If there are more than two nodes in the segment, there must be intermediate nodes that we can remove.
                if (countNodes > 2) then
                   basicDescendant => nodeDescendant%basic()
                   ! Find the final sibling of the final node in the segment.
                   nodePriorSibling => nodeDescendant
                   do while (associated(nodePriorSibling%sibling))
                      nodePriorSibling => nodePriorSibling%sibling
                   end do
                   ! Update parent pointers for any siblings of the first node in the segment.
                   nodeSibling => nodeBase%sibling
                   do while (associated(nodeSibling))
                      nodeSibling%parent => nodeDescendant
                      nodeSibling        => nodeSibling   %sibling
                   end do
                   ! Walk up the branch moving any siblings to the descendant node, and destroy the intermediate nodes.
                   nodeBranch => nodeBase%parent
                   do while (.not.associated(nodeBranch,nodeDescendant))
                      ! Find the next node in the branch for future reference.
                      nodeNext    => nodeBranch%parent
                      ! Update parent pointers for any siblings.
                      nodeSibling => nodeBranch%sibling
                      do while (associated(nodeSibling))
                         nodeSibling%parent => nodeDescendant%parent
                         nodeSibling        => nodeSibling   %sibling
                      end do
                      ! Move siblings to the node at the end of the segment.
                      nodePriorSibling%sibling => nodeBranch%sibling
                      do while (associated(nodePriorSibling%sibling))
                         nodePriorSibling => nodePriorSibling%sibling
                      end do
                      ! Destroy the intermediate node.
                      call nodeBranch%destroy()
                      deallocate(nodeBranch)
                      ! Move to the next node.
                      nodeBranch => nodeNext
                   end do
                   ! Link the nodes at the start and end of the segment.
                   nodeDescendant%firstChild => nodeBase
                   nodeBase      %parent     => nodeDescendant
                end if
                ! Go back to the parent node and reset the base point for further segments.
                nodeDescendant => nodeDescendant%parent
                nodeBase       => nodeDescendant
                basic          => nodeBase      %basic ()
                countNodes     =  0
             end if
             ! Move up the branch, stopping when we are no longer the primary progenitor.
             if (nodeDescendant%isPrimaryProgenitor()) then
                nodeDescendant => nodeDescendant%parent
             else
                nodeDescendant => null()
             end if
          end do
       end do
       ! Move to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do
    return
  end subroutine consolidateBranchesOperatePreEvolution

