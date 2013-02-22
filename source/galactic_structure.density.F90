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

  ! Module scope variables used in mapped function.
  double precision, dimension(3) :: positionSphericalShared
  integer                        :: componentTypeShared,massTypeShared,weightByShared,weightIndexShared
  logical                        :: haloLoadedShared
  !$omp threadprivate(positionSphericalShared,massTypeShared,componentTypeShared,weightByShared,weightIndexShared,haloLoadedShared)
 
contains

  double precision function Galactic_Structure_Density(thisNode,position,coordinateSystem,componentType,massType,haloLoaded)
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
    type            (treeNode         ), intent(inout), pointer  :: thisNode
    integer                            , intent(in),    optional :: componentType,massType,coordinateSystem
    logical                            , intent(in),    optional :: haloLoaded
    double precision                   , intent(in)              :: position(3)
    procedure       (Component_Density),                pointer  :: componentDensityFunction
    integer                                                      :: coordinateSystemActual
    double precision                                             :: componentDensity

    ! Determine position in spherical coordinate system to use.
    if (present(coordinateSystem)) then
       coordinateSystemActual=coordinateSystem
    else
       coordinateSystemActual=coordinateSystemSpherical
    end if
    select case (coordinateSystemActual)
    case (coordinateSystemSpherical)
       positionSphericalShared=position
    case (coordinateSystemCylindrical)
       positionSphericalShared=Coordinates_Cylindrical_To_Spherical(position)
    case (coordinateSystemCartesian)
       positionSphericalShared=Coordinates_Cartesian_To_Spherical  (position)
    case default
       call Galacticus_Error_Report('Galactic_Structure_Density','unknown coordinate system type')
    end select
    ! Determine which mass type to use.
    if (present(massType)) then
       massTypeShared=massType
    else
       massTypeShared=massTypeAll
    end if
    ! Determine which component type to use.
    if (present(componentType)) then
       componentTypeShared=componentType
    else
       componentTypeShared=componentTypeAll
    end if
    ! Determine which component type to use.
    if (present(haloLoaded)) then
       haloLoadedShared=haloLoaded
    else
       haloLoadedShared=.true.
    end if
    ! Call routines to supply the densities for all components.
    componentDensityFunction => Component_Density
    Galactic_Structure_Density=thisNode%mapDouble0(componentDensityFunction,reductionSummation)
    !# <include directive="densityTask" type="functionCall" functionType="function" returnParameter="componentDensity">
    !#  <functionArgs>thisNode,positionSphericalShared,massTypeShared,componentTypeShared,haloLoadedShared</functionArgs>
    !#  <onReturn>Galactic_Structure_Density=Galactic_Structure_Density+componentDensity</onReturn>
    include 'galactic_structure.density.tasks.inc'
    !# </include>
    return
  end function Galactic_Structure_Density
  
  double precision function Component_Density(component)
    !% Unary function returning the density in a component. Suitable for mapping over components.
    use Galacticus_Nodes
    implicit none
    class(nodeComponent), intent(inout) :: component
 
    Component_Density=component%density(positionSphericalShared,componentTypeShared&
         &,massTypeShared,haloLoadedShared)
    return
  end function Component_Density

end module Galactic_Structure_Densities
