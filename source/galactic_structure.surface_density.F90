!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

!% Contains a module which implements calculations of the surface density at a specific position.

module Galactic_Structure_Surface_Densities
  !% Implements calculations of the surface density at a specific position.
  use :: Galacticus_Nodes, only : treeNode
  implicit none
  private
  public :: Galactic_Structure_Surface_Density, Galactic_Structure_Radius_Enclosing_Surface_Density

  ! Module scope variables used in mapping over components.
  integer                                  :: componentTypeShared      , massTypeShared, weightByShared, weightIndexShared
  double precision                         :: surfaceDensityRoot
  double precision          , dimension(3) :: positionCylindricalShared
  type            (treeNode), pointer      :: activeNode
  !$omp threadprivate(massTypeShared,componentTypeShared,positionCylindricalShared,weightByShared,weightIndexShared,surfaceDensityRoot,activeNode)

contains

  double precision function Galactic_Structure_Surface_Density(node,position,coordinateSystem,componentType,massType,weightBy,weightIndex)
    !% Compute the density (of given {\normalfont \ttfamily massType}) at the specified {\normalfont \ttfamily position}. Assumes that galactic structure has already
    !% been computed.
    use :: Coordinate_Systems        , only : Coordinates_Cartesian_To_Cylindrical, Coordinates_Spherical_To_Cylindrical
    use :: Galactic_Structure_Options, only : componentTypeAll                    , coordinateSystemCartesian           , coordinateSystemCylindrical, coordinateSystemSpherical, &
          &                                   massTypeAll                         , weightByMass                        , weightIndexNull
    use :: Galacticus_Error          , only : Galacticus_Error_Report
    use :: Galacticus_Nodes          , only : optimizeForSurfaceDensitySummation  , optimizeforsurfacedensitysummation  , reductionSummation         , reductionsummation       , &
          &                                   treeNode
    implicit none
    type            (treeNode                 ), intent(inout)           :: node
    integer                                    , intent(in   ), optional :: componentType                     , coordinateSystem, &
         &                                                                  massType                          , weightBy        , &
         &                                                                  weightIndex
    double precision                           , intent(in   )           :: position                       (3)
    procedure       (Component_Surface_Density), pointer                 :: componentSurfaceDensityFunction
    integer                                                              :: coordinateSystemActual

    call Galactic_Structure_Surface_Density_Defaults(componentType,massType,weightBy,weightIndex)
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
       call Galacticus_Error_Report('unknown coordinate system type'//{introspection:location})
    end select
    ! Call routines to supply the densities for all components.
    componentSurfaceDensityFunction => Component_Surface_Density
    Galactic_Structure_Surface_Density=node%mapDouble0(componentSurfaceDensityFunction,reductionSummation,optimizeFor=optimizeForSurfaceDensitySummation)
    return
  end function Galactic_Structure_Surface_Density

  double precision function Component_Surface_Density(component)
    !% Unary function returning the surface density in a component. Suitable for mapping over components.
    use :: Galacticus_Nodes, only : nodeComponent
    implicit none
    class(nodeComponent), intent(inout) :: component

    Component_Surface_Density=component%surfaceDensity(positionCylindricalShared,componentTypeShared,massTypeShared,weightByShared,weightIndexShared)
    return
  end function Component_Surface_Density

  double precision function Galactic_Structure_Radius_Enclosing_Surface_Density(node,surfaceDensity,componentType,massType,weightBy,weightIndex)
    !% Return the radius enclosing a given mass (or fractional mass) in {\normalfont \ttfamily node}.
    use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScale      , darkMatterHaloScaleClass
    use :: Kind_Numbers           , only : kind_int8
    use :: Root_Finder            , only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    type            (treeNode                ), intent(inout), target   :: node
    integer                                   , intent(in   ), optional :: componentType                    , massType   , &
         &                                                                 weightBy                         , weightIndex
    double precision                          , intent(in   ), optional :: surfaceDensity
    class           (darkMatterHaloScaleClass), pointer                 :: darkMatterHaloScale_
    type            (rootFinder              ), save                    :: finder
    logical                                   , save                    :: finderConstructed   =.false.
    !$omp threadprivate(finder,finderConstructed)
    double precision                          , save                    :: radiusPrevious      =-huge(0.0d0)
    integer         (kind_int8               ), save                    :: uniqueIDPrevious    =-1_kind_int8
    !$omp threadprivate(radiusPrevious,uniqueIDPrevious)
    double precision                                                    :: radiusGuess

    ! Set default options.
    call Galactic_Structure_Surface_Density_Defaults(componentType,massType,weightBy,weightIndex)
    activeNode => node
    ! Initialize our root finder.
    if (.not.finderConstructed) then
       finder           =rootFinder(                                                             &
            &                       rootFunction                 =Surface_Density_Root         , &
            &                       rangeExpandDownward          =0.5d0                        , &
            &                       rangeExpandUpward            =2.0d0                        , &
            &                       rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &
            &                       rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
            &                       rangeExpandType              =rangeExpandMultiplicative    , &
            &                       toleranceAbsolute            =0.0d+0                       , &
            &                       toleranceRelative            =1.0d-6                         &
            &                      )
       finderConstructed=.true.
    end if
    ! Solve for the radius.
    activeNode         => node
    surfaceDensityRoot =  surfaceDensity
    if (node%uniqueID() == uniqueIDPrevious) then
       radiusGuess          =  radiusPrevious
    else
       darkMatterHaloScale_ => darkMatterHaloScale              (        )
       radiusGuess          =  darkMatterHaloScale_%virialRadius(node)
    end if
    Galactic_Structure_Radius_Enclosing_Surface_Density=finder%find(rootGuess=radiusGuess)
    uniqueIDPrevious                                   =node%uniqueID()
    radiusPrevious                                     =Galactic_Structure_Radius_Enclosing_Surface_Density
    return
  end function Galactic_Structure_Radius_Enclosing_Surface_Density

  double precision function Surface_Density_Root(radius)
    !% Root function used in solving for the radius that encloses a given surface density.
    use :: Galactic_Structure_Options, only : coordinateSystemCylindrical
    implicit none
    double precision, intent(in   ) :: radius

    ! Evaluate the root function.
    Surface_Density_Root=Galactic_Structure_Surface_Density(activeNode,[radius,0.0d0,0.0d0],coordinateSystemCylindrical,componentTypeShared,massTypeShared,weightByShared,weightIndexShared)-surfaceDensityRoot
  end function Surface_Density_Root

  subroutine Galactic_Structure_Surface_Density_Defaults(componentType,massType,weightBy,weightIndex)
    !% Set the default values for options in the surface density functions.
    use :: Galactic_Structure_Options, only : componentTypeAll       , massTypeAll, weightByLuminosity, weightByMass
    use :: Galacticus_Error          , only : Galacticus_Error_Report
    implicit none
    integer, intent(in   ), optional :: componentType, massType, weightBy, weightIndex

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
    return
  end subroutine Galactic_Structure_Surface_Density_Defaults

end module Galactic_Structure_Surface_Densities
