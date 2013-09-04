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

!% Contains a module which prunes branches in a tree that do not contain an ``essential'' node.

module Merger_Trees_Prune_Non_Essential
  !% Prunes branches in a tree that do not contain an ``essential'' node.
  use Kind_Numbers
  implicit none
  private
  public :: Merger_Tree_Non_Essential_Branches

  ! Flag indicating if module is initialized.
  logical                          :: pruneNonEssentialModuleInitialized=.false.

  ! Flag indicating if pruning is required.
  logical                          :: mergerTreePruneNonEssential

  ! Index of the essential node.
  integer         (kind=kind_int8) :: mergerTreePruningNonEssentialID

  ! Time at which the essential node is required.
  double precision                 :: mergerTreePruningNonEssentialTime

contains

  !# <mergerTreePreEvolveTask>
  !#   <unitName>Merger_Tree_Non_Essential_Branches</unitName>
  !# </mergerTreePreEvolveTask>
  subroutine Merger_Tree_Non_Essential_Branches(thisTree)
    !% Prune branches from {\tt thisTree}.
    use Galacticus_Nodes
    use Input_Parameters
    implicit none
    type   (mergerTree        ), intent(in   ), target :: thisTree
    type   (treeNode          ), pointer               :: nextNode   , previousNode, thisNode, essentialNode
    class  (nodeComponentBasic), pointer               :: thisBasic
    type   (mergerTree        ), pointer               :: currentTree

    ! Check if module is initialized.
    if (.not.pruneNonEssentialModuleInitialized) then
       !$omp critical (Merger_Tree_Non_Essential_Branches_Initialize)
       if (.not.pruneNonEssentialModuleInitialized) then
          ! Get parameter specifying if pruning is required.
          !@ <inputParameter>
          !@   <name>mergerTreePruneNonEssential</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not to prune non-essential merger trees branches prior to evolution.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreePruneNonEssential',mergerTreePruneNonEssential,defaultValue=.false.)
          if (mergerTreePruneNonEssential) then
             !@ <inputParameter>
             !@   <name>mergerTreePruningNonEssentialID</name>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     ID of the essential node to avoid pruning.
             !@   </description>
             !@   <type>integer</type>
             !@   <cardinality>1</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter('mergerTreePruningNonEssentialID',mergerTreePruningNonEssentialID)
             !@ <inputParameter>
             !@   <name>mergerTreePruningNonEssentialTime</name>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     Time of the essential node to avoid pruning.
             !@   </description>
             !@   <type>real</type>
             !@   <cardinality>1</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter('mergerTreePruningNonEssentialTime',mergerTreePruningNonEssentialTime)
          end if
          ! Flag that module is initialized.
          pruneNonEssentialModuleInitialized=.true.
       end if
       !$omp end critical (Merger_Tree_Non_Essential_Branches_Initialize)
    end if

    ! Prune tree if necessary.
    if (mergerTreePruneNonEssential) then
       ! Iterate over trees.
       currentTree => thisTree
       do while (associated(currentTree))
          ! Find the essential node.
          essentialNode => currentTree%getNode(mergerTreePruningNonEssentialID)
          ! Trace the essential node to the required time.
          thisBasic => essentialNode%basic()
          do while (thisBasic%time() < mergerTreePruningNonEssentialTime .and. associated(essentialNode%parent))
             essentialNode => essentialNode%parent
             thisBasic => essentialNode%basic()
          end do          
          
          ! Get root node of the tree.
          thisNode => currentTree%baseNode
          ! Walk the tree, pruning branches.
          do while (associated(thisNode))
             ! Record the parent node to which we will return.
             previousNode => thisNode%parent
             if     (                                               &
                  &  .not.                                          &
                  &   (                                             &
                  &     thisNode     %isProgenitorOf(essentialNode) &
                  &    .or.                                         &
                  &     essentialNode%isProgenitorOf(     thisNode) &
                  &   )                                             &
                  & ) then               
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
          ! Move to the next tree.
          currentTree => currentTree%nextTree
       end do
    end if
    return
  end subroutine Merger_Tree_Non_Essential_Branches

end module Merger_Trees_Prune_Non_Essential
