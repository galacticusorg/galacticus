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

!% Contains a module which prunes branches which can contain no baryons from merger trees.

module Merger_Trees_Prune_Baryons
  !% Prunes branches which can contain no baryons from merger trees.
  implicit none
  private
  public :: Merger_Tree_Prune_Baryons

  ! Flag indicating if module is initialized.
  logical :: pruneBaryonsModuleInitialized=.false.

  ! Flag indicating if pruning is required.
  logical :: mergerTreePruneBaryons

contains

  !# <mergerTreePreEvolveTask>
  !#   <unitName>Merger_Tree_Prune_Baryons</unitName>
  !# </mergerTreePreEvolveTask>
  subroutine Merger_Tree_Prune_Baryons(thisTree)
    !% Prune branches from {\normalfont \ttfamily thisTree}.
    use Galacticus_Nodes
    use Input_Parameters
    use Accretion_Halos
    implicit none
    type   (mergerTree        ), intent(in   ), target :: thisTree
    type   (treeNode          ), pointer               :: nextNode          , previousNode, thisNode
    class  (nodeComponentBasic), pointer               :: thisBasicComponent
    type   (mergerTree        ), pointer               :: currentTree
    class  (accretionHaloClass), pointer               :: accretionHalo_
    logical                                            :: didPruning

    ! Check if module is initialized.
    if (.not.pruneBaryonsModuleInitialized) then
       !$omp critical (Merger_Tree_Prune_Baryons_Initialize)
       if (.not.pruneBaryonsModuleInitialized) then
          ! Get parameter specifying if pruning is required.
          !@ <inputParameter>
          !@   <name>mergerTreePruneBaryons</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not to prune merger trees prior to evolution.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreePruneBaryons',mergerTreePruneBaryons,defaultValue=.false.)
          ! Flag that module is initialized.
          pruneBaryonsModuleInitialized=.true.
       end if
       !$omp end critical (Merger_Tree_Prune_Baryons_Initialize)
    end if

    ! Prune tree if necessary.
    if (mergerTreePruneBaryons) then
       ! Get required objects.
       accretionHalo_ => accretionHalo()
       ! Iterate over trees.
       currentTree => thisTree
       do while (associated(currentTree))
          didPruning=.true.
          do while (didPruning)
             didPruning=.false.
             ! Get root node of the tree.
             thisNode           => currentTree%baseNode
             thisBasicComponent => thisNode%basic()
             ! Check if the tree has any baryons.
             if (.not.accretionHalo_%branchHasBaryons(thisNode)) then
                ! Entire tree can be pruned. Destroy all but this base node. (Leaving just
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
                   if (.not.accretionHalo_%branchHasBaryons(thisNode)) then
                      didPruning=.true.
                      ! Decouple from other nodes.
                      previousBasicComponent => previousNode%basic()
                      call Merger_Tree_Prune_Unlink_Parent(thisNode,previousNode,.not.accretionHalo_%branchHasBaryons(previousNode))
                      ! Clean the branch.
                      call Merger_Tree_Prune_Clean_Branch(thisNode)
                      ! Destroy the branch.
                      call currentTree%destroyBranch(thisNode)
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
  end subroutine Merger_Tree_Prune_Baryons

end module Merger_Trees_Prune_Baryons
