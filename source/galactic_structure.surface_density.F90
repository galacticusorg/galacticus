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

!% Contains a module which implements calculations of the surface density at a specific position.

module Galactic_Structure_Surface_Densities
  !% Implements calculations of the surface density at a specific position.
  use Galacticus_Nodes
  use Galactic_Structure_Options
  implicit none
  private
  public :: Galactic_Structure_Surface_Density

  ! Module scope variables used in mapping over components.
  integer                        :: componentTypeShared      , massTypeShared
  logical                        :: haloLoadedShared
  double precision, dimension(3) :: positionCylindricalShared
  !$omp threadprivate(massTypeShared,componentTypeShared,haloLoadedShared,positionCylindricalShared)
contains

  double precision function Galactic_Structure_Surface_Density(thisNode,position,coordinateSystem,componentType,massType,haloLoaded)
    !% Compute the density (of given {\tt massType}) at the specified {\tt position}. Assumes that galactic structure has already
    !% been computed.
    use Galacticus_Error
    use Coordinate_Systems
    !# <include directive="surfaceDensityTask" type="moduleUse">
    include 'galactic_structure.surface_density.tasks.modules.inc'
    !# </include>
    implicit none
    type            (treeNode                 ), intent(inout)          , pointer :: thisNode
    integer                                    , intent(in   ), optional          :: componentType                     , coordinateSystem, &
         &                                                                           massType
    logical                                    , intent(in   ), optional          :: haloLoaded
    double precision                           , intent(in   )                    :: position                       (3)
    procedure       (Component_Surface_Density)                         , pointer :: componentSurfaceDensityFunction
    integer                                                                       :: coordinateSystemActual
    double precision                                                              :: componentDensity

    ! Determine position in cylindrical coordinate system to use.
    if (present(coordinateSystem)) then
       coordinateSystemActual=coordinateSystem
    else
       coordinateSystemActual=coordinateSystemCylindrical
    end if
    select case (coordinateSystemActual)
    case (coordinateSystemSpherical)
       positionCylindricalShared=Coordinates_Spherical_To_Cylindrical (position)
    case (coordinateSystemCylindrical)
       positionCylindricalShared=position
    case (coordinateSystemCartesian)
       positionCylindricalShared=Coordinates_Cartesian_To_Cylindrical(position)
    case default
       call Galacticus_Error_Report('Galactic_Structure_Surface_Density','unknown coordinate system type')
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
    ! Determine whether halo loading is to be used.
    if (present(haloLoaded)) then
       haloLoadedShared=haloLoaded
    else
       haloLoadedShared=.true.
    end if
    ! Call routines to supply the densities for all components.
    componentSurfaceDensityFunction => Component_Surface_Density
    Galactic_Structure_Surface_Density=thisNode%mapDouble0(componentSurfaceDensityFunction,reductionSummation)
    !# <include directive="surfaceDensityTask" type="functionCall" functionType="function" returnParameter="componentDensity">
    !#  <functionArgs>thisNode,positionCylindricalShared,massTypeShared,componentTypeShared,haloLoadedShared</functionArgs>
    !#  <onReturn>Galactic_Structure_Surface_Density=Galactic_Structure_Surface_Density+componentDensity</onReturn>
    include 'galactic_structure.surface_density.tasks.inc'
    !# </include>
    return
  end function Galactic_Structure_Surface_Density

  double precision function Component_Surface_Density(component)
    !% Unary function returning the surface density in a component. Suitable for mapping over components.
    implicit none
    class(nodeComponent), intent(inout) :: component

    Component_Surface_Density=component%surfaceDensity(positionCylindricalShared,componentTypeShared&
         &,massTypeShared,haloLoadedShared)
    return
  end function Component_Surface_Density

end module Galactic_Structure_Surface_Densities
