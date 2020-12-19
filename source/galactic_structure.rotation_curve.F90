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

!% Contains a module which implements calculations of the rotation curve as a specified radius.

module Galactic_Structure_Rotation_Curves
  !% Implements calculations of the rotation curve as a specified radius.
  implicit none
  private
  public :: Galactic_Structure_Rotation_Curve

  ! Module scope variables used in mapping over components.
  integer          :: componentTypeShared, massTypeShared
  double precision :: radiusShared
  !$omp threadprivate(massTypeShared,componentTypeShared,radiusShared)
contains

  double precision function Galactic_Structure_Rotation_Curve(node,radius,componentType,massType)
    !% Solve for the rotation curve a given radius. Assumes that galactic structure has already been solved for.
    use :: Galactic_Structure_Options, only : componentTypeAll                 , massTypeAll
    use :: Galacticus_Nodes          , only : optimizeForRotationCurveSummation, reductionSummation, treeNode
    !# <include directive="rotationCurveTask" type="moduleUse">
    include 'galactic_structure.rotation_curve.tasks.modules.inc'
    !# </include>
    implicit none
    type            (treeNode                ), intent(inout)                    :: node
    integer                                   , intent(in   ), optional          :: componentType                 , massType
    double precision                          , intent(in   )                    :: radius
    procedure       (Component_Rotation_Curve)                         , pointer :: componentRotationCurveFunction
    double precision                                                             :: componentVelocity             , rotationCurveSquared

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
    ! Store the radius.
    radiusShared=radius
    ! Call routines to supply the velocities for all components.
    componentRotationCurveFunction => Component_Rotation_Curve
    rotationCurveSquared=node%mapDouble0(componentRotationCurveFunction,reductionSummation,optimizeFor=optimizeForRotationCurveSummation)
    !# <include directive="rotationCurveTask" type="functionCall" functionType="function" returnParameter="componentVelocity">
    !#  <functionArgs>node,radiusShared,componentTypeShared,massTypeShared</functionArgs>
    !#  <onReturn>rotationCurveSquared=rotationCurveSquared+componentVelocity**2</onReturn>
    include 'galactic_structure.rotation_curve.tasks.inc'
    !# </include>
    ! We've added velocities in quadrature, so now take the square root.
    Galactic_Structure_Rotation_Curve=sqrt(rotationCurveSquared)
    return
  end function Galactic_Structure_Rotation_Curve

  double precision function Component_Rotation_Curve(component)
    !% Unary function returning the squared rotation curve in a component. Suitable for mapping over components.
    use :: Galacticus_Nodes, only : nodeComponent
    implicit none
    class(nodeComponent), intent(inout) :: component

    Component_Rotation_Curve=component%rotationCurve(radiusShared,componentTypeShared,massTypeShared)**2
    return
  end function Component_Rotation_Curve

end module Galactic_Structure_Rotation_Curves
