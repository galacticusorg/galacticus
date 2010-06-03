!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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
  use Tree_Nodes
  private
  public :: Radius_Get, Radius_Set, Velocity_Get, Velocity_Set

  ! Pointers to get and set procedures.
  procedure(Get_Template), pointer :: Radius_Get => null(), Velocity_Get => null()
  procedure(Set_Template), pointer :: Radius_Set => null(), Velocity_Set => null()
  !$omp threadprivate(Radius_Get,Radius_Set,Velocity_Get,Velocity_Set)

  abstract interface
     double precision function Get_Template(thisNode)
       import treeNode
       type(treeNode), pointer, intent(inout) :: thisNode
     end function Get_Template
  end interface
  abstract interface
     subroutine Set_Template(thisNode,value)
       import treeNode
       type(treeNode),   pointer, intent(inout) :: thisNode
       double precision,          intent(in)    :: value
     end subroutine Set_Template
  end interface

end module Galactic_Structure_Radius_Solver_Procedures
