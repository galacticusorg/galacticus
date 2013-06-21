!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which prunes hierarchy below a given depth in merger trees.

module Merger_Trees_Prune_Hierarchy
  !% Prunes hierarchy below a given depth in merger trees.
  implicit none
  private
  public :: Merger_Tree_Prune_Hierarchy

  ! Flag indicating if module is initialized.
  logical :: pruneHierarchyModuleInitialized=.false.

  ! Depth below which hierarchy should be pruned.
  integer :: mergerTreePruneHierarchyAtDepth

contains

  !# <mergerTreePreEvolveTask>
  !#   <unitName>Merger_Tree_Prune_Hierarchy</unitName>
  !# </mergerTreePreEvolveTask>
  subroutine Merger_Tree_Prune_Hierarchy(thisTree)
    !% Prune hierarchy from {\tt thisTree}.
    use Merger_Trees
    use Galacticus_Nodes
    use Input_Parameters
    implicit none
    type   (mergerTree), intent(in   ), target :: thisTree
    type   (treeNode  ), pointer               :: nextNode      , previousNode, thisNode, workNode
    type   (mergerTree), pointer               :: currentTree
    logical                                    :: didPruning
    integer                                    :: hierarchyDepth

    ! Check if module is initialized.
    if (.not.pruneHierarchyModuleInitialized) then
       !$omp critical (Merger_Tree_Prune_Hierarchy_Initialize)
       if (.not.pruneHierarchyModuleInitialized) then
          ! Get parameter specifying if pruning is required.
          !@ <inputParameter>
          !@   <name>mergerTreePruneHierarchyAtDepth</name>
          !@   <defaultValue>0</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The depth in the hierarchy at which to prune merger trees. (Zero indicates to not prune.)
          !@   </description>
          !@   <type>integer</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreePruneHierarchyAtDepth',mergerTreePruneHierarchyAtDepth,defaultValue=0)
          ! Flag that module is initialized.
          pruneHierarchyModuleInitialized=.true.
       end if
       !$omp end critical (Merger_Tree_Prune_Hierarchy_Initialize)
    end if

    ! Prune tree if necessary.
    if (mergerTreePruneHierarchyAtDepth > 0) then
       ! Iterate over trees.
       currentTree => thisTree
       do while (associated(currentTree))
          didPruning=.true.
          do while (didPruning)
             didPruning=.false.
             ! Get root node of the tree.
             thisNode => currentTree%baseNode
             ! Walk the tree, pruning hierarchy.
             do while (associated(thisNode))
                ! Record the parent node to which we will return.
                previousNode => thisNode%parent

                ! Find the depth in the hierarchy.
                hierarchyDepth=0
                workNode => thisNode
                do while (associated(workNode))
                   ! Increment hierarchy depth if this node is not the main progenitor.
                   if (.not.workNode%isPrimaryProgenitor().and.associated(workNode%parent)) hierarchyDepth=hierarchyDepth+1
                   workNode => workNode%parent
                end do

                ! Prune if this node is sufficiently deep in the hierarchy.
                if (hierarchyDepth >= mergerTreePruneHierarchyAtDepth) then
                   didPruning=.true.
                   ! Decouple from other nodes.
                   if (thisNode%isPrimaryProgenitorOf(previousNode)) then
                      previousNode%firstChild => thisNode%sibling
                   else
                      nextNode => previousNode%firstChild
                      do while (.not.associated(nextNode%sibling,thisNode))
                         nextNode => nextNode%sibling
                      end do
                      nextNode%sibling => thisNode%sibling
                   end if
                   call currentTree%destroyBranch(thisNode)
                   ! Return to parent node.
                   thisNode => previousNode
                end if
                call thisNode%walkTree(thisNode)
             end do
          end do
          ! Move to the next tree.
          currentTree => currentTree%nextTree
       end do
    end if

    return
  end subroutine Merger_Tree_Prune_Hierarchy

end module Merger_Trees_Prune_Hierarchy
