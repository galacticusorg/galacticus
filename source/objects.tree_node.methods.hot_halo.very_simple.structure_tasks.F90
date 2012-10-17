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

!% Contains a module which implements structure tasks for the standard hot halo component.

module Tree_Node_Methods_Hot_Halo_Structure_Tasks_Very_Simple
  !% Implements structure tasks for the standard hot halo component.
  use Tree_Node_Methods_Hot_Halo_Data_Very_Simple
  implicit none
  private
  public :: Hot_Halo_Very_Simple_Density
  
contains

  !# <densityTask>
  !#  <unitName>Hot_Halo_Very_Simple_Density</unitName>
  !# </densityTask>
  subroutine Hot_Halo_Very_Simple_Density(thisNode,positionSpherical,massType,componentType,componentDensity)
    !% Computes the density at a given position in the hot halo.
    use Tree_Nodes
    use Hot_Halo_Density_Profile
    use Galactic_Structure_Options
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    integer,          intent(in)             :: massType,componentType
    double precision, intent(in)             :: positionSpherical(3)
    double precision, intent(out)            :: componentDensity
    
    componentDensity=0.0d0
    if (.not.methodSelected                           ) return
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeHotHalo)) return
    if (.not.thisNode%componentExists(componentIndex)) return
    select case (massType)
    case (massTypeAll,massTypeBaryonic,massTypeGaseous)
       componentDensity=Hot_Halo_Density(thisNode,positionSpherical(1))
    end select
    return
  end subroutine Hot_Halo_Very_Simple_Density

end module Tree_Node_Methods_Hot_Halo_Structure_Tasks_Very_Simple
