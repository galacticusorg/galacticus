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

!% Defines functions bound to the merger tree object class.

subroutine Merger_Tree_Destroy(thisTree)
  !% Destroys the entire merger tree.
  implicit none
  class(mergerTree), intent(inout) :: thisTree

  select type (thisTree)
  type is (mergerTree)
     ! Destroy all nodes.
     if (associated(thisTree%baseNode)) call thisTree%destroyBranch(thisTree%baseNode)
     ! Destroy the HDF5 group associated with this tree.
     call thisTree%hdf5Group%destroy()
  end select
  return
end subroutine Merger_Tree_Destroy

subroutine Merger_Tree_Destroy_Branch(thisTree,thisNode)
  !% Destroy a branch of a tree which begins at {\tt thisNode}.
  implicit none
  class(mergerTree), intent(inout)          :: thisTree
  type (treeNode  ), intent(inout), pointer :: thisNode
  type (treeNode  )               , pointer :: destroyNode, nextNode

  ! Descend to the tip of the branch.
  call thisNode%walkBranchWithSatellites(thisNode,nextNode)
  ! Loop over all tree nodes.
  do while (associated(nextNode))
     ! Keep of a record of the current node, so that we can destroy it.
     destroyNode => nextNode

     ! Walk to the next node in the tree.
     call destroyNode%walkBranchWithSatellites(thisNode,nextNode)

     ! If the node about to be destroyed is the primary progenitor of its parent we must move the child pointer of the parent to
     ! point to the node's sibling. This is necessary as parent-child pointers are used to establish satellite status and so
     ! will be utilized when walking the tree. Failure to do this can result in attempts to use dangling pointers.
     if (associated(destroyNode%parent)) then
        if (associated(destroyNode%parent%firstChild,destroyNode)) destroyNode%parent%firstChild => destroyNode%sibling
     end if

     ! Destroy the current node.
     call destroyNode%destroy()
     deallocate(destroyNode)
  end do
  ! Destroy the base node of the branch.
  if (associated(thisNode%parent)) then
     if (associated(thisNode%parent%firstChild,thisNode)) thisNode%parent%firstChild => thisNode%sibling
  end if
  call thisNode%destroy()
  deallocate(thisNode)
  return
end subroutine Merger_Tree_Destroy_Branch

subroutine Tree_Node_Create_Initialize
  !% Initializes tree node create by calling all relevant initialization routines.
  use Input_Parameters
  implicit none

  if (.not.treeNodeCreateInitialized) then

     ! Initialize tree node methods.
     call Galacticus_Nodes_Initialize()

     ! Flag that tree node methods are now initialized.
     treeNodeCreateInitialized=.true.

  end if
  return
end subroutine Tree_Node_Create_Initialize

function Merger_Tree_Node_Get(thisTree,nodeIndex)
  !% Return a pointer to a node in {\tt thisTree} given the index of the node.
  implicit none
  class  (mergerTree    ), intent(in   ) :: thisTree
  integer(kind=kind_int8), intent(in   ) :: nodeIndex
  type   (treeNode      ), pointer       :: Merger_Tree_Node_Get, thisNode

  Merger_Tree_Node_Get => null()
  thisNode => thisTree%baseNode
  do while (associated(thisNode))
     if (thisNode%index() == nodeIndex) then
        Merger_Tree_Node_Get => thisNode
        return
     end if
     ! <gfortan 4.6> explicitly specify the target as thisNode since we can't use the "_Same_Node" tree walking procedures.
     call thisNode%walkTreeWithSatellites(thisNode)
  end do
  return
end function Merger_Tree_Node_Get
