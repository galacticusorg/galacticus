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

!% Contains a module which prunes branches below a given mass threshold from merger trees.

module Merger_Trees_Prune_Branches
  !% Prunes branches below a given mass threshold from merger trees.
  implicit none
  private
  public :: Merger_Tree_Prune_Branches

  ! Flag indicating if module is initialized.
  logical          :: pruneBranchesModuleInitialized=.false.

  ! Flag indicating if pruning is required.
  logical          :: mergerTreePruneBranches

  ! Mass below which branches should be pruned.
  double precision :: mergerTreePruningMassThreshold

contains

  !# <mergerTreePreEvolveTask>
  !#   <unitName>Merger_Tree_Prune_Branches</unitName>
  !# </mergerTreePreEvolveTask>
  subroutine Merger_Tree_Prune_Branches(thisTree)
    !% Prune branches from {\tt thisTree}.
    use Merger_Trees
    use Galacticus_Nodes
    use Input_Parameters
    implicit none
    type (mergerTree        ), intent(in), target :: thisTree
    type (treeNode          ), pointer            :: thisNode,nextNode,destroyNode,previousNode
    class(nodeComponentBasic), pointer            :: thisBasicComponent
    type (mergerTree        ), pointer            :: currentTree
    logical                                       :: didPruning

    ! Check if module is initialized.
    if (.not.pruneBranchesModuleInitialized) then
       !$omp critical (Merger_Tree_Prune_Branches_Initialize)
       if (.not.pruneBranchesModuleInitialized) then
          ! Get parameter specifying if pruning is required.
          !@ <inputParameter>
          !@   <name>mergerTreePruneBranches</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not to prune merger trees prior to evolution.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreePruneBranches',mergerTreePruneBranches,defaultValue=.false.)
          !@ <inputParameter>
          !@   <name>mergerTreePruningMassThreshold</name>
          !@   <defaultValue>0</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Threshold mass below which merger tree branches should be pruned.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreePruningMassThreshold',mergerTreePruningMassThreshold,defaultValue=0.0d0)
          ! Flag that module is initialized.
          pruneBranchesModuleInitialized=.true.
       end if
       !$omp end critical (Merger_Tree_Prune_Branches_Initialize)
    end if

    ! Prune tree if necessary.
    if (mergerTreePruneBranches) then
       ! Iterate over trees.
       currentTree => thisTree
       do while (associated(currentTree))
          didPruning=.true.
          do while (didPruning)
             didPruning=.false.
             ! Get root node of the tree.       
             thisNode           => currentTree%baseNode
             thisBasicComponent => thisNode%basic()
             if (thisBasicComponent%mass() < mergerTreePruningMassThreshold) then
                ! Entire tree is below threshold. Destroy all but this base node. (Leaving just
                ! the base node makes the tree inert - i.e. it can not do anything.)
                thisNode => thisNode%firstChild
                do while (associated(thisNode))
                   nextNode => thisNode%sibling
                   call currentTree%destroyBranch(thisNode)
                   thisNode => nextNode
                end do
             else
                ! Walk the tree, pruning branches.
                do while (associated(thisNode))
                   thisBasicComponent => thisNode%basic()
                   ! Record the parent node to which we will return.
                   previousNode => thisNode%parent
                   if (thisBasicComponent%mass() < mergerTreePruningMassThreshold) then
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
                      ! Should call the "destroyBranch" method on "currentTree" here, but seems to cause a compiler crash under gFortran
                      ! v4.4. So, instead, do the branch destroy manually.
                      nextNode => thisNode
                      do while (associated(nextNode))
                         destroyNode => nextNode
                         call destroyNode%walkBranch(thisNode,nextNode)
                         call destroyNode%destroy()
                      end do
                      ! Return to parent node.
                      thisNode => previousNode
                   end if
                   call thisNode%walkTree(thisNode)
                end do
             end if
          end do
          ! Move to the next tree.
          currentTree => currentTree%nextTree
       end do
    end if

    return
  end subroutine Merger_Tree_Prune_Branches

end module Merger_Trees_Prune_Branches
