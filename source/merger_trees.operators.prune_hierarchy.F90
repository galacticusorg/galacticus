!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !% Contains a module which implements a merger tree operator which prunes branches below a
  !% given level in the substructure hierarchy.

  !# <mergerTreeOperator name="mergerTreeOperatorPruneHierarchy">
  !#  <description>
  !#   Contains a module which implements a merger tree operator which prunes branches below a
  !#   given level in the substructure hierarchy. In any tree, the primary progenitor of the
  !#   base node has substructure hierarchy depth 0. A branch which connects directly to this
  !#   primary progenitor branch has substructure hierarchy depth 1, while a branch which
  !#   connects directly to that branch has substructure hierarchy depth 2, and so on. The tree
  !#   is pruned of all branches of hierarchy depth equal to or greater than the value provided
  !#   by the {\normalfont \ttfamily [hierarchyDepth]} parameter.
  !#  </description>
  !# </mergerTreeOperator>
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorPruneHierarchy
     !% A merger tree operator class which prunes branches below a given level in the
     !% substructure hierarchy.
     private
     integer :: hierarchyDepth
   contains
     procedure :: operatePreEvolution => pruneHierarchyOperatePreEvolution
  end type mergerTreeOperatorPruneHierarchy

  interface mergerTreeOperatorPruneHierarchy
     !% Constructors for the prune-hierarchy merger tree operator class.
     module procedure pruneHierarchyConstructorParameters
     module procedure pruneHierarchyConstructorInternal
  end interface mergerTreeOperatorPruneHierarchy

contains

  function pruneHierarchyConstructorParameters(parameters)
    !% Constructor for the prune-hierarchy merger tree operator class which takes a parameter set as input.
    implicit none
    type   (mergerTreeOperatorPruneHierarchy)                :: pruneHierarchyConstructorParameters
    type   (inputParameters                 ), intent(inout) :: parameters
    integer                                                  :: hierarchyDepth

    !# <inputParameter>
    !#   <name>hierarchyDepth</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1</defaultValue>
    !#   <description>The depth in the substructure hierarchy at which to prune a tree.</description>
    !#   <type>integer</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    pruneHierarchyConstructorParameters=pruneHierarchyConstructorInternal(hierarchyDepth)
    !# <inputParametersValidate source="parameters"/>
    return
  end function pruneHierarchyConstructorParameters

  function pruneHierarchyConstructorInternal(hierarchyDepth)
    !% Internal constructor for the prune-hierarchy merger tree operator class.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    type   (mergerTreeOperatorPruneHierarchy)                :: pruneHierarchyConstructorInternal
    integer                                  , intent(in   ) :: hierarchyDepth

    if (hierarchyDepth < 1) call Galacticus_Error_Report('[hierarchyDepth] > 0 is required'//{introspection:location})
    pruneHierarchyConstructorInternal%hierarchyDepth=hierarchyDepth
    return
  end function pruneHierarchyConstructorInternal

  subroutine pruneHierarchyOperatePreEvolution(self,tree)
    !% Perform a prune-hierarchy operation on a merger tree.
    use :: Merger_Tree_Walkers           , only : mergerTreeWalkerIsolatedNodes
    use :: Merger_Trees_Pruning_Utilities, only : Merger_Tree_Prune_Clean_Branch, Merger_Tree_Prune_Uniqueify_IDs, Merger_Tree_Prune_Unlink_Parent
    implicit none
    class  (mergerTreeOperatorPruneHierarchy), intent(inout), target  :: self
    type   (mergerTree                      ), intent(inout), target  :: tree
    type   (mergerTree                      )               , pointer :: treeCurrent
    type   (treeNode                        )               , pointer :: node          , nodeWork
    type   (mergerTreeWalkerIsolatedNodes   )                         :: treeWalker
    logical                                                           :: didPruning
    integer                                                           :: hierarchyDepth

    ! Iterate over nodes.
    treeCurrent => tree
    do while (associated(treeCurrent))
       didPruning=.true.
       do while (didPruning)
          didPruning=.false.
          ! Walk the tree, pruning hierarchy.
          treeWalker=mergerTreeWalkerIsolatedNodes(treeCurrent)
          do while (treeWalker%next(node))
             ! Find the depth in the hierarchy.
             hierarchyDepth =  0
             nodeWork       => node
             do while (associated(nodeWork))
                ! Increment hierarchy depth if this node is not the main progenitor.
                if (.not.nodeWork%isPrimaryProgenitor().and.associated(nodeWork%parent)) hierarchyDepth=hierarchyDepth+1
                nodeWork => nodeWork%parent
             end do
             ! Prune if this node is sufficiently deep in the hierarchy.
             if (hierarchyDepth >= self%hierarchyDepth) then
                didPruning=.true.
                ! Decouple from other nodes.
                hierarchyDepth =  0
                nodeWork       => node%parent
                do while (associated(nodeWork))
                   ! Increment hierarchy depth if this node is not the main progenitor.
                   if (.not.nodeWork%isPrimaryProgenitor().and.associated(nodeWork%parent)) hierarchyDepth=hierarchyDepth+1
                   nodeWork => nodeWork%parent
                end do
                ! Set the tree walker back to the base node so that it exits - we've changed the tree structure so need to
                ! begin the walk again.
                call treeWalker%setNode(treeCurrent%baseNode)
                ! Unlink the node.
                call Merger_Tree_Prune_Unlink_Parent(node,node%parent,hierarchyDepth >= self%hierarchyDepth,.true.)
                ! Clean the branch.
                call Merger_Tree_Prune_Clean_Branch(node)
                ! Destroy the branch.
                call node%destroyBranch()
                deallocate(node)
             end if
          end do
       end do
       ! Move to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do
    ! Uniqueify nodes.
    call Merger_Tree_Prune_Uniqueify_IDs(tree)
    return
  end subroutine pruneHierarchyOperatePreEvolution

