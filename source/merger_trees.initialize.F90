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

!% Contains a module which implements initialization of merger tree structures.

module Merger_Trees_Initialize
  !% Implements initialization of merger tree structures.
  implicit none
  private
  public :: Merger_Tree_Initialize

contains

  subroutine Merger_Tree_Initialize(thisTree)
    !% Walk through all nodes of a tree and call any routines that requested to perform initialization tasks.
    use Merger_Trees
    use Galacticus_Nodes
    !# <include directive="mergerTreeInitializeTask" type="moduleUse">
    include 'merger_trees.initialize.tasks.modules.inc'
    !# </include>
    implicit none
    type(mergerTree), intent(inout) :: thisTree  
    type(treeNode  ), pointer       :: thisNode  
                                              
    if (.not.thisTree%initialized) then
       thisNode => thisTree%baseNode
       do while (associated(thisNode))
          ! Call subroutines to perform any necessary initialization of this node.
          !# <include directive="mergerTreeInitializeTask" type="functionCall" functionType="void">
          !#  <functionArgs>thisNode</functionArgs>
          include 'merger_trees.initialize.tasks.inc'
          !# </include>

          call thisNode%walkTreeWithSatellites(thisNode)
       end do
       thisTree%initialized=.true.
    end if

    return
  end subroutine Merger_Tree_Initialize

end module Merger_Trees_Initialize
