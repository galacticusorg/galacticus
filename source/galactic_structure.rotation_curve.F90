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

!% Contains a module which implements calculations of the rotation curve as a specified radius.

module Galactic_Structure_Rotation_Curves
  !% Implements calculations of the rotation curve as a specified radius.
  use ISO_Varying_String
  use Tree_Nodes
  use Galactic_Structure_Options
  implicit none
  private
  public :: Galactic_Structure_Rotation_Curve

  ! Variables used in root finding.
  integer                  :: massTypeRoot,componentTypeRoot
  double precision         :: massRoot
  type(treeNode),  pointer :: activeNode
  !$omp threadprivate(massRoot,massTypeRoot,componentTypeRoot,activeNode)

contains

  double precision function Galactic_Structure_Rotation_Curve(thisNode,radius,massType,componentType)
    !% Solve for the rotation curve a given radius. Assumes that galacticus structure has already been solved for.,
    !# <include directive="rotationCurveTask" type="moduleUse">
    include 'galactic_structure.rotation_curve.tasks.modules.inc'
    !# </include>
    implicit none
    type(treeNode),   intent(inout), pointer  :: thisNode
    integer,          intent(in),    optional :: massType,componentType
    double precision, intent(in)              :: radius
    integer                                   :: massTypeActual,componentTypeActual
    double precision                          :: componentVelocity

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

    ! Initialize to zero mass.
    Galactic_Structure_Rotation_Curve=0.0d0

    ! Call routines to supply the velocities for all components.
    !# <include directive="rotationCurveTask" type="code" action="subroutine">
    !#  <subroutineArgs>thisNode,radius,massTypeActual,componentTypeActual,componentVelocity</subroutineArgs>
    !#  <subroutineAction>Galactic_Structure_Rotation_Curve=Galactic_Structure_Rotation_Curve+componentVelocity**2</subroutineAction>
    include 'galactic_structure.rotation_curve.tasks.inc'
    !# </include>

    ! We've added velocities in quadrature, so now take the square root.
    Galactic_Structure_Rotation_Curve=dsqrt(Galactic_Structure_Rotation_Curve)

    return
  end function Galactic_Structure_Rotation_Curve
  
end module Galactic_Structure_Rotation_Curves
