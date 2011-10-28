!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


program Tests_Tree_Branch_Destroy
  use Unit_Tests
  use Memory_Management
  use Merger_Trees
  use Tree_Nodes
  use Kind_Numbers
  implicit none
  type(mergerTree),       pointer :: thisTree
  type(treeNodeList)              :: nodes(5)
  integer(kind=kind_int8)         :: iNode

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
  nodes(1)%node%childNode => nodes(2)%node
  nodes(3)%node%childNode => nodes(4)%node
  
  ! Set parent nodes.
  nodes(2)%node%parentNode => nodes(1)%node
  nodes(3)%node%parentNode => nodes(1)%node
  nodes(4)%node%parentNode => nodes(3)%node
  nodes(5)%node%parentNode => nodes(3)%node

  ! Set sibling nodes.
  nodes(2)%node%siblingNode => nodes(3)%node
  nodes(4)%node%siblingNode => nodes(5)%node

  ! Destroy the branch rooted at node 2.
  call thisTree%destroyBranch(nodes(2)%node)
  
  ! The child node of node 1 should have been shifted to point to node 3.
  call Assert('Child node of node 1 updated',nodes(1)%node%childNode%index(),3_kind_int8)

  ! Destroy the tree. If danling pointers are not fixed during tree destruction, this will trigger errors when running inside
  ! Valgrind.
  call thisTree%destroy()

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Tests_Tree_Branch_Destroy
