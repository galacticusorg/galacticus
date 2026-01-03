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
  Implements a merger tree operator which prunes branches below a
  given level in the substructure hierarchy.
  !!}

  !![
  <mergerTreeOperator name="mergerTreeOperatorPruneHierarchy">
   <description>
    A merger tree operator class module which prunes branches below a given level in the substructure hierarchy. In any tree,
    the primary progenitor of the base node has substructure hierarchy depth 0. A branch which connects directly to this
    primary progenitor branch has substructure hierarchy depth 1, while a branch which connects directly to that branch has
    substructure hierarchy depth 2, and so on. The tree is pruned of all branches of hierarchy depth equal to or greater than
    the value provided by the {\normalfont \ttfamily [hierarchyDepth]} parameter.
   </description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorPruneHierarchy
     !!{
     A merger tree operator class which prunes branches below a given level in the
     substructure hierarchy.
     !!}
     private
     integer :: hierarchyDepth
   contains
     procedure :: operatePreEvolution => pruneHierarchyOperatePreEvolution
  end type mergerTreeOperatorPruneHierarchy

  interface mergerTreeOperatorPruneHierarchy
     !!{
     Constructors for the prune-hierarchy merger tree operator class.
     !!}
     module procedure pruneHierarchyConstructorParameters
     module procedure pruneHierarchyConstructorInternal
  end interface mergerTreeOperatorPruneHierarchy

contains

  function pruneHierarchyConstructorParameters(parameters) result(self)
    !!{
    Constructor for the prune-hierarchy merger tree operator class which takes a parameter set as input.
    !!}
    implicit none
    type   (mergerTreeOperatorPruneHierarchy)                :: self
    type   (inputParameters                 ), intent(inout) :: parameters
    integer                                                  :: hierarchyDepth

    !![
    <inputParameter>
      <name>hierarchyDepth</name>
      <source>parameters</source>
      <defaultValue>1</defaultValue>
      <description>The depth in the substructure hierarchy at which to prune a tree.</description>
    </inputParameter>
    !!]
    self=mergerTreeOperatorPruneHierarchy(hierarchyDepth)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function pruneHierarchyConstructorParameters

  function pruneHierarchyConstructorInternal(hierarchyDepth) result(self)
    !!{
    Internal constructor for the prune-hierarchy merger tree operator class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type   (mergerTreeOperatorPruneHierarchy)                :: self
    integer                                  , intent(in   ) :: hierarchyDepth

    if (hierarchyDepth < 1) call Error_Report('[hierarchyDepth] > 0 is required'//{introspection:location})
    self%hierarchyDepth=hierarchyDepth
    return
  end function pruneHierarchyConstructorInternal

  subroutine pruneHierarchyOperatePreEvolution(self,tree)
    !!{
    Perform a prune-hierarchy operation on a merger tree.
    !!}
    use :: Galacticus_Nodes              , only : treeNodeLinkedList
    use :: Merger_Tree_Walkers           , only : mergerTreeWalkerIsolatedNodes
    use :: Merger_Trees_Pruning_Utilities, only : Merger_Tree_Prune_Clean_Branch, Merger_Tree_Prune_Uniqueify_IDs, Merger_Tree_Prune_Unlink_Parent
    implicit none
    class  (mergerTreeOperatorPruneHierarchy), intent(inout), target  :: self
    type   (mergerTree                      ), intent(inout), target  :: tree
    type   (mergerTree                      )               , pointer :: treeCurrent
    type   (treeNode                        )               , pointer :: node          , nodeWork       , &
         &                                                               nodePrune     , nodePrevious
    type   (treeNodeLinkedList              )               , pointer :: nodeListHead  , nodeListCurrent, &
         &                                                               nodeListNext
    type   (mergerTreeWalkerIsolatedNodes   )                         :: treeWalker
    integer                                                           :: hierarchyDepth, i

    ! Iterate over nodes.
    treeCurrent => tree
    do while (associated(treeCurrent))
       ! Keep a record of the previous unpruned node, so that we can step back to this node after pruning a branch (and then
       ! continue the tree walk).
       nodePrevious => null()
       ! Walk the tree, pruning hierarchy.
       treeWalker=mergerTreeWalkerIsolatedNodes(treeCurrent)          
       do while (treeWalker%next(node))
          ! Find the depth in the hierarchy. Keep a list of branch nodes as we walk up the tree determining hierarchy level. We
          ! can then use this to determine which branch node corresponds to the hierarchy level at which we want to prune.
          hierarchyDepth  =  0
          nodeListHead    => null()
          nodelistCurrent => null()
          nodeWork        => node
          nodePrune       => null()
          do while (associated(nodeWork))
             ! Increment hierarchy depth if this node is not the main progenitor.
             if (.not.nodeWork%isPrimaryProgenitor().and.associated(nodeWork%parent)) then
                hierarchyDepth=hierarchyDepth+1
                ! Add this node to our list of branch nodes.
                if (associated(nodeListCurrent)) then
                   allocate(nodeListCurrent%next)
                   nodeListCurrent => nodeListCurrent%next
                else
                   allocate(nodeListHead)
                   nodeListCurrent => nodeListHead
                end if
                nodeListCurrent%node => nodeWork
                nodeListCurrent%next => null()
             end if
             nodeWork => nodeWork%parent
          end do
          ! If the node is sufficiently deep, find the node at the required hierarchy level and prune at that point.
          if (hierarchyDepth >= self%hierarchyDepth) then
             ! Find which node we should prune.
             nodeListCurrent => nodeListHead
             if (hierarchyDepth > self%hierarchyDepth) then
                do i=1,hierarchyDepth-self%hierarchyDepth
                   nodeListCurrent => nodeListCurrent%next
                end do
             end if
             nodePrune => nodeListCurrent%node
             ! Set the tree walker back to the previous unpruned node so that it will resume the walk from that point.
             call treeWalker%setNode(nodePrevious)
             ! Unlink the node.
             call Merger_Tree_Prune_Unlink_Parent(nodePrune,nodePrune%parent,.false.,.true.)
             ! Clean the branch.
             call Merger_Tree_Prune_Clean_Branch(nodePrune)
             ! Destroy the branch.
             call nodePrune%destroyBranch()
             deallocate(nodePrune)
          else
             ! Node is not being pruned, keep a record of this node so that we can step back to it after the next pruning.
             nodePrevious => node
          end if
          ! Clean up the list of branch nodes.
          if (associated(nodeListHead)) then
             nodeListCurrent => nodeListHead
             do while (associated(nodeListCurrent))
                nodeListNext => nodeListCurrent%next
                deallocate(nodeListCurrent)
                nodeListCurrent => nodeListNext
             end do
             nullify(nodeListHead)
          end if
       end do
       ! Move to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do
    ! Uniqueify nodes.
    call Merger_Tree_Prune_Uniqueify_IDs(tree)
    return
  end subroutine pruneHierarchyOperatePreEvolution

