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

!% Contains a module which implements dumping of a node's data to screen.

module Tree_Nodes_Dump
  !% Implements dumping of a node's data to screen.
  implicit none
  private
  public :: Node_Dump
  
contains

  subroutine Node_Dump(thisNode)
    !% Dump the properties of {\tt thisNode} to screen.
    use Tree_Nodes
    !# <include directive="nodeDumpTask" type="moduleUse">
    include 'objects.tree_node.dump.moduleUse.inc'
    !# </include>
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    write (0,'(a,i10,a)'   ) 'Dumping data for node [',thisNode%index(),']'
    write (0,'(1x,a20,i10)') 'parent node: ',thisNode%parentNode%index()
    write (0,'(1x,a20,i10)') 'sibling node: ',thisNode%siblingNode%index()
    write (0,'(1x,a20,i10)') 'child node: ',thisNode%childNode%index()
    !# <include directive="nodeDumpTask" type="code" action="subroutine">
    !#  <subroutineArgs>thisNode</subroutineArgs>
    include 'objects.tree_node.dump.inc'
    !# </include>
    
    return
  end subroutine Node_Dump

end module Tree_Nodes_Dump
