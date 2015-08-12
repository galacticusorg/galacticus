!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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
     final     ::            pruneHierarchyDestructor
     procedure :: operate => pruneHierarchyOperate
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
    type   (inputParameters                 ), intent(in   ) :: parameters
    integer                                                  :: hierarchyDepth
    !# <inputParameterList label="allowedParameterNames" />
        
    call parameters%checkParameters(allowedParameterNames)
    !# <inputParameter>
    !#   <name>hierarchyDepth</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1</defaultValue>
    !#   <description>The depth in the substructure hierarchy at which to prune a tree.</description>
    !#   <type>integer</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    pruneHierarchyConstructorParameters=pruneHierarchyConstructorInternal(hierarchyDepth)
    return
  end function pruneHierarchyConstructorParameters

  function pruneHierarchyConstructorInternal(hierarchyDepth)
    !% Internal constructor for the prune-hierarchy merger tree operator class.
    use Galacticus_Error
    implicit none
    type   (mergerTreeOperatorPruneHierarchy)                :: pruneHierarchyConstructorInternal
    integer                                  , intent(in   ) :: hierarchyDepth

    if (hierarchyDepth < 1) call Galacticus_Error_Report('pruneHierarchyConstructorInternal','[hierarchyDepth] > 0 is required')
    pruneHierarchyConstructorInternal%hierarchyDepth=hierarchyDepth
    return
  end function pruneHierarchyConstructorInternal

  elemental subroutine pruneHierarchyDestructor(self)
    !% Destructor for the merger tree operator function class.
    implicit none
    type(mergerTreeOperatorPruneHierarchy), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine pruneHierarchyDestructor

  subroutine pruneHierarchyOperate(self,tree)
    !% Perform a prune-hierarchy operation on a merger tree.
    use Merger_Trees_Pruning_Utilities
    implicit none
    class  (mergerTreeOperatorPruneHierarchy), intent(inout)          :: self
    type   (mergerTree                      ), intent(inout), target  :: tree
    type   (mergerTree                      )               , pointer :: treeCurrent
    type   (treeNode                        )               , pointer :: node          , nodePrevious, &
         &                                                               nodeWork      , nodeNext
    logical                                                           :: didPruning
    integer                                                           :: hierarchyDepth
    
    ! Iterate over trees.
    treeCurrent => tree
    do while (associated(treeCurrent))
       didPruning=.true.
       do while (didPruning)
          didPruning=.false.
          ! Get root node of the tree.
          node => treeCurrent%baseNode
          ! Walk the tree, pruning hierarchy.
          do while (associated(node))
             ! Record the parent node to which we will return.
             nodePrevious => node%parent
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
                if (node%isPrimaryProgenitorOf(nodePrevious)) then
                   nodePrevious%firstChild => node%sibling
                else
                   nodeNext => nodePrevious%firstChild
                   do while (.not.associated(nodeNext%sibling,node))
                      nodeNext => nodeNext%sibling
                   end do
                   nodeNext%sibling => node%sibling
                end if
                ! Decouple from other nodes.
                hierarchyDepth =  0
                nodeWork       => nodePrevious
                do while (associated(nodeWork))
                   ! Increment hierarchy depth if this node is not the main progenitor.
                   if (.not.nodeWork%isPrimaryProgenitor().and.associated(nodeWork%parent)) hierarchyDepth=hierarchyDepth+1
                   nodeWork => nodeWork%parent
                end do
                call Merger_Tree_Prune_Unlink_Parent(node,nodePrevious,hierarchyDepth >= self%hierarchyDepth,.true.)
                ! Clean the branch.
                call Merger_Tree_Prune_Clean_Branch(node)
                ! Destroy the branch.
                call treeCurrent%destroyBranch(node)
                ! Return to parent node.
                node => nodePrevious
             end if
             node => node%walkTree()
          end do
       end do
       ! Move to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do
    return
  end subroutine pruneHierarchyOperate
  
