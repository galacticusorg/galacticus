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

program Tests_Tree_Branch_Destroy
  use Unit_Tests
  use Memory_Management
  use Merger_Trees
  use Galacticus_Nodes
  use Kind_Numbers
  implicit none
  type   (mergerTree    ), pointer :: thisTree
  type   (treeNodeList  )          :: nodes   (5)
  integer(kind=kind_int8)          :: iNode

  ! Read in basic code memory usage.
  call Code_Memory_Usage('tests.tree_branch_destroy.size')

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Tree branch destruction: Avoid dangling pointers")

  ! Create tree.
  allocate(thisTree)

  ! Create nodes.
  do iNode=1,5
     call thisTree%createNode(nodes(iNode)%node)
  end do
  thisTree%baseNode => nodes(1)%node

  ! Set indices of nodes.
  call nodes(1)%node%indexSet(1_kind_int8)
  call nodes(2)%node%indexSet(2_kind_int8)
  call nodes(3)%node%indexSet(3_kind_int8)
  call nodes(4)%node%indexSet(4_kind_int8)
  call nodes(5)%node%indexSet(5_kind_int8)

  ! Set child nodes.
  nodes(1)%node%firstChild => nodes(2)%node
  nodes(3)%node%firstChild => nodes(4)%node

  ! Set parent nodes.
  nodes(2)%node%parent => nodes(1)%node
  nodes(3)%node%parent => nodes(1)%node
  nodes(4)%node%parent => nodes(3)%node
  nodes(5)%node%parent => nodes(3)%node

  ! Set sibling nodes.
  nodes(2)%node%sibling => nodes(3)%node
  nodes(4)%node%sibling => nodes(5)%node

  ! Destroy the branch rooted at node 2.
  call thisTree%destroyBranch(nodes(2)%node)

  ! The child node of node 1 should have been shifted to point to node 3.
  call Assert('Child node of node 1 updated',nodes(1)%node%firstChild%index(),3_kind_int8)

  ! Destroy the tree. If danling pointers are not fixed during tree destruction, this will trigger errors when running inside
  ! Valgrind.
  call thisTree%destroy()

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Tests_Tree_Branch_Destroy
