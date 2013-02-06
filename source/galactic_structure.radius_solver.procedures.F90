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

!% Contains a module which holds procedure pointers used by the galactic structure radii solver subsystem.

module Galactic_Structure_Radius_Solver_Procedures
  !% Holds procedure pointers used by the galactic structure radii solver subsystem.
  use Galacticus_Nodes
  implicit none
  public

  abstract interface
     double precision function Structure_Get_Template(thisNode)
       import treeNode
       type(treeNode), pointer, intent(inout) :: thisNode
     end function Structure_Get_Template
  end interface
  abstract interface
     subroutine Structure_Set_Template(thisNode,value)
       import treeNode
       type(treeNode),   pointer, intent(inout) :: thisNode
       double precision,          intent(in)    :: value
     end subroutine Structure_Set_Template
  end interface

end module Galactic_Structure_Radius_Solver_Procedures
