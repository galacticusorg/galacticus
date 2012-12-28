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

!% Contains a module which implements calculations of the density at a specific position.

module Galactic_Structure_Densities
  !% Implements calculations of the density at a specific position.
  implicit none
  private
  public :: Galactic_Structure_Density

contains

  double precision function Galactic_Structure_Density(thisNode,position,coordinateSystem,massType,componentType)
    !% Compute the density (of given {\tt massType}) at the specified {\tt position}. Assumes that galactic structure has already
    !% been computed.
    use Galacticus_Nodes
    use Galactic_Structure_Options
    use ISO_Varying_String
    use Galacticus_Error
    use Input_Parameters
    use Coordinate_Systems
    !# <include directive="densityTask" type="moduleUse">
    include 'galactic_structure.density.tasks.modules.inc'
    !# </include>
    implicit none
    type(treeNode),   intent(inout), pointer  :: thisNode
    integer,          intent(in),    optional :: massType,componentType,coordinateSystem
    double precision, intent(in)              :: position(3)
    integer                                   :: massTypeActual,coordinateSystemActual,componentTypeActual
    double precision                          :: positionSpherical(3),componentDensity

    ! Determine position in spherical coordinate system to use.
    if (present(coordinateSystem)) then
       coordinateSystemActual=coordinateSystem
    else
       coordinateSystemActual=coordinateSystemSpherical
    end if
    select case (coordinateSystemActual)
    case (coordinateSystemSpherical)
       positionSpherical=position
    case (coordinateSystemCylindrical)
       positionSpherical=Coordinates_Cylindrical_To_Spherical(position)
    case (coordinateSystemCartesian)
       positionSpherical=Coordinates_Cartesian_To_Spherical(position)
    case default
       call Galacticus_Error_Report('Galactic_Structure_Density','unknown coordinate system type')
    end select

    ! Determine which mass type to use.
    if (present(massType)) then
       massTypeActual=massType
    else
       massTypeActual=massTypeAll
    end if

    ! Determine which component type to use.
    if (present(componentType)) then
       componentTypeActual=componentType
    else
       componentTypeActual=componentTypeAll
    end if

    ! Initialize to zero density.
    Galactic_Structure_Density=0.0d0

    ! Call routines to supply the densities for all components.
    !# <include directive="densityTask" type="functionCall" functionType="void">
    !#  <functionArgs>thisNode,positionSpherical,massTypeActual,componentTypeActual,componentDensity</functionArgs>
    !#  <onReturn>Galactic_Structure_Density=Galactic_Structure_Density+componentDensity</onReturn>
    include 'galactic_structure.density.tasks.inc'
    !# </include>

    return
  end function Galactic_Structure_Density
  
end module Galactic_Structure_Densities
