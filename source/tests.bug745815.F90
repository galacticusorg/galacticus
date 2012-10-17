!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

program Tests_Bug745815
  !% Tests for regression of Bug \#745815 (http://bugs.launchpad.net/galacticus/+bug/745815): Skipping of a node during a tree
  !% walk.
  use Unit_Tests
  use Input_Parameters
  use ISO_Varying_String
  use Memory_Management
  use Merger_Trees
  use Tree_Nodes
  use Kind_Numbers
  implicit none
  type(varying_string)            :: parameterFile
  type(mergerTree)                :: thisTree
  type(treeNodeList)              :: nodes(5)
  logical                         :: nodeFound(5)
  type(treeNode),         pointer :: thisNode
  integer(kind=kind_int8)         :: iNode

  ! Read in basic code memory usage.
  call Code_Memory_Usage('tests.bug745815.size')

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Bug #745815: Node skip during tree-walk")

  ! Open the parameter file.
  parameterFile='testSuite/parameters/bug745815.xml'
  call Input_Parameters_File_Open(parameterFile)
  
  ! Create nodes.
  do iNode=1,5
     call thisTree%createNode(nodes(iNode)%node)
  end do

  ! Set indices of nodes.
  call nodes(1)%node%indexSet(100017990003559_kind_int8)
  call nodes(2)%node%indexSet(100017990003560_kind_int8)
  call nodes(3)%node%indexSet(100017990003561_kind_int8)
  call nodes(4)%node%indexSet(100017990003562_kind_int8)
  call nodes(5)%node%indexSet(100017990003571_kind_int8)
  
  ! Set child nodes.
  nodes(1)%node%childNode => nodes(2)%node
  nodes(2)%node%childNode => nodes(3)%node
  nodes(3)%node%childNode => nodes(4)%node
  
  ! Set parent nodes.
  nodes(2)%node%parentNode => nodes(1)%node
  nodes(3)%node%parentNode => nodes(2)%node
  nodes(4)%node%parentNode => nodes(3)%node
  nodes(5)%node%parentNode => nodes(3)%node
  
  ! Set satellite nodes.
  nodes(3)%node%satelliteNode => nodes(5)%node
  
  ! Walk the tree, with satellites.
  nodeFound=.false.
  thisNode => nodes(1)%node
  do while (associated(thisNode))
     do iNode=1,5
        if (nodes(iNode)%node%index() == thisNode%index()) nodeFound(iNode)=.true.
     end do
     call thisNode%walkTreeWIthSatellites(thisNode)
  end do
  call Assert('All nodes walked to',all(nodeFound),.true.)

  ! Destroy nodes.
  do iNode=1,5
     call nodes(iNode)%node%destroy()
  end do
  
  ! Close the parameter file.
  call Input_Parameters_File_Close  

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Tests_Bug745815
