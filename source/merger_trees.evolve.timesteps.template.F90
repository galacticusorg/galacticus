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

!% Contains a module which defines the template for tasks performed at the end of timesteps.

module Merger_Trees_Evolve_Timesteps_Template
  !% Defines the template for tasks performed at the end of timesteps.
  use Merger_Trees
  use Galacticus_Nodes
  implicit none
  public

  abstract interface
     subroutine End_Of_Timestep_Task_Template(thisTree,thisNode,deadlockStatus)
       import mergerTree, treeNode
       type(mergerTree), intent(in)             :: thisTree
       type(treeNode),   intent(inout), pointer :: thisNode
       integer,          intent(inout)          :: deadlockStatus
     end subroutine End_Of_Timestep_Task_Template
  end interface

end module Merger_Trees_Evolve_Timesteps_Template
