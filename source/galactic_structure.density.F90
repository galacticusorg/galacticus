!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
!!    Andrew Benson <abenson@carnegiescience.edu>
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
  integer                        :: componentTypeShared    , massTypeShared   , &
       &                            weightByShared         , weightIndexShared
  !$omp threadprivate(positionSphericalShared,massTypeShared,componentTypeShared,weightByShared,weightIndexShared)
contains

  double precision function Galactic_Structure_Density(node,position,coordinateSystem,componentType,massType,weightBy,weightIndex)
    !% Compute the density (of given {\normalfont \ttfamily massType}) at the specified {\normalfont \ttfamily position}. Assumes that galactic structure has already
    !% been computed.
    use :: Coordinate_Systems        , only : Coordinates_Cartesian_To_Spherical, Coordinates_Cylindrical_To_Spherical
    use :: Galactic_Structure_Options, only : componentTypeAll                  , coordinateSystemCartesian           , coordinateSystemCylindrical, coordinateSystemSpherical, &
          &                                   massTypeAll                       , weightByLuminosity                  , weightByMass
    use :: Galacticus_Error          , only : Galacticus_Error_Report
    use :: Galacticus_Nodes          , only : optimizeForDensitySummation       , reductionSummation                  , treeNode
    !# <include directive="densityTask" type="moduleUse">
    include 'galactic_structure.density.tasks.modules.inc'
    !# </include>
    implicit none
    type            (treeNode         ), intent(inout)           :: node
    integer                            , intent(in   ), optional :: componentType              , coordinateSystem, &
         &                                                          massType                   , weightBy        , &
         &                                                          weightIndex
    double precision                   , intent(in   )           :: position                (3)
    procedure       (Component_Density), pointer                 :: componentDensityFunction
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
       call Galacticus_Error_Report('unknown coordinate system type'//{introspection:location})
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
    ! Determine which weighting to use.
    if (present(weightBy)) then
       weightByShared=weightBy
       select case (weightByShared)
       case (weightByLuminosity)
          if (.not.present(weightIndex)) call Galacticus_Error_Report('weightIndex should be specified for luminosity weighting'//{introspection:location})
          weightIndexShared=weightIndex
       end select
    else
       weightByShared=weightByMass
    end if

    ! Call routines to supply the densities for all components.
    componentDensityFunction => Component_Density
    Galactic_Structure_Density=node%mapDouble0(componentDensityFunction,reductionSummation,optimizeFor=optimizeForDensitySummation)
    !# <include directive="densityTask" type="functionCall" functionType="function" returnParameter="componentDensity">
    !#  <functionArgs>node,positionSphericalShared,componentTypeShared,massTypeShared,weightByShared,weightIndexShared</functionArgs>
    !#  <onReturn>Galactic_Structure_Density=Galactic_Structure_Density+componentDensity</onReturn>
    include 'galactic_structure.density.tasks.inc'
    !# </include>
    return
  end function Galactic_Structure_Density

  double precision function Component_Density(component)
    !% Unary function returning the density in a component. Suitable for mapping over components.
    use :: Galacticus_Nodes, only : nodeComponent
    implicit none
    class(nodeComponent), intent(inout) :: component

    Component_Density=component%density(positionSphericalShared,componentTypeShared&
         &,massTypeShared,weightByShared,weightIndexShared)
    return
  end function Component_Density

end module Galactic_Structure_Densities
