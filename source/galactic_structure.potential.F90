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

!+    Contributions to this file made by:  Stéphane Mangeon, Andrew Benson.

!% Contains a module which implements calculations of gravitationl potential.

module Galactic_Structure_Potentials
  !% Implements calculations of the gravitational potential.
  use Galacticus_Nodes
  use Galactic_Structure_Options
  implicit none
  private
  public :: Galactic_Structure_Potential

  ! Module scope variables used in mapping over components.
  integer          :: componentTypeShared, massTypeShared
  logical          :: haloLoadedShared
  double precision :: radiusShared
  !$omp threadprivate(massTypeShared,componentTypeShared,haloLoadedShared,radiusShared)
contains

  double precision function Galactic_Structure_Potential(thisNode,radius,componentType,massType,haloLoaded)
    !% Solve for the gravitational potential at a given radius. Assumes the galactic structure has already been computed.
    !# <include directive="potentialTask" type="moduleUse">
    include 'galactic_structure.potential.tasks.modules.inc'
    !# </include>
    implicit none
    type            (treeNode           ), intent(inout)          , pointer :: thisNode
    integer                              , intent(in   ), optional          :: componentType             , massType
    logical                              , intent(in   ), optional          :: haloLoaded
    double precision                     , intent(in   )                    :: radius
    procedure       (Component_Potential)                         , pointer :: componentPotentialFunction
    double precision                                                        :: componentPotential

    ! Determine which component type to use.
    if (present(componentType)) then
       componentTypeShared=componentType
    else
       componentTypeShared=componentTypeAll
    end if
    ! Determine which mass type to use.
    if (present(massType)) then
       massTypeShared     =massType
    else
       massTypeShared     =massTypeAll
    end if
    ! Determine whether halo loading is to be used.
    if (present(haloLoaded)) then
       haloLoadedShared=haloLoaded
    else
       ! Note that the default option is currently false, as adiabatic contraction effects on the dark matter potential are not
       ! accounted for.
       haloLoadedShared=.false.
    end if
    ! Store the radius.
    radiusShared=radius
    ! Call routines to supply the potential for all components.
    componentPotentialFunction => Component_Potential
    Galactic_Structure_Potential=thisNode%mapDouble0(componentPotentialFunction,reductionSummation)
    !# <include directive="potentialTask" type="functionCall" functionType="function" returnParameter="componentPotential">
    !#  <functionArgs>thisNode,radiusShared,componentTypeShared,massTypeShared,haloLoadedShared</functionArgs>
    !#  <onReturn>Galactic_Structure_Potential=Galactic_Structure_Potential+componentPotential</onReturn>
    include 'galactic_structure.potential.tasks.inc'
    !# </include>
    return
  end function Galactic_Structure_Potential

  double precision function Component_Potential(component)
    !% Unary function returning the potential in a component. Suitable for mapping over components.
    implicit none
    class(nodeComponent), intent(inout) :: component

    Component_Potential=component%potential(radiusShared,componentTypeShared,massTypeShared,haloLoadedShared)
    return
  end function Component_Potential

end module Galactic_Structure_Potentials












