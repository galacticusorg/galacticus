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

!+    Contributions to this file made by:  Stéphane Mangeon, Andrew Benson.

!% Contains a module which implements calculations of the gradient of the rotation curve.

module Galactic_Structure_Rotation_Curve_Gradients
  !% Implements calculations of the rotation curve gradient
  private
  public :: Galactic_Structure_Rotation_Curve_Gradient

  ! Module scope variables used in mapping over components.
  integer          :: componentTypeShared, massTypeShared
  double precision :: radiusShared
  !$omp threadprivate(massTypeShared,componentTypeShared,radiusShared)

contains

  double precision function Galactic_Structure_Rotation_Curve_Gradient(node,radius,componentType,massType)
    !% Solve for the rotation curve gradient at a given radius. Assumes the galactic structure has already been computed.
    use :: Galactic_Structure_Options        , only : componentTypeAll                         , massTypeAll
    use :: Galactic_Structure_Rotation_Curves, only : Galactic_Structure_Rotation_Curve
    use :: Galacticus_Nodes                  , only : optimizeForRotationCurveGradientSummation, reductionSummation, treeNode
    !# <include directive="rotationCurveGradientTask" type="moduleUse">
    include 'galactic_structure.rotation_curve.gradient.tasks.modules.inc'
    !# </include>
    implicit none
    type            (treeNode                         ), intent(inout)           :: node
    integer                                            , intent(in   ), optional :: componentType                         , massType
    double precision                                   , intent(in   )           :: radius
    procedure       (Component_Rotation_Curve_Gradient)               , pointer  :: componentRotationCurveGradientFunction
    double precision                                                             :: componentRotationCurveGradient

    ! Determine which component type to use.
    if (present(componentType)) then
       componentTypeShared=componentType
    else
       componentTypeShared=componentTypeAll
    end if
    ! Determine which mass type to use
    if (present(massType)) then
       massTypeShared=massType
    else
       massTypeShared=massTypeAll
    end if
    ! Store the radius.
    radiusShared=radius
    ! Call routines to supply the gradient for all components' rotation curves. Specifically, the returned quantities are
    ! d(V^2)/dr so that they can be summed directly.
    componentRotationCurveGradientFunction => Component_Rotation_Curve_Gradient
    Galactic_Structure_Rotation_Curve_Gradient=node%mapDouble0(componentRotationCurveGradientFunction,reductionSummation,optimizeFor=optimizeForRotationCurveGradientSummation)
    !# <include directive="rotationCurveGradientTask" type="functionCall" functionType="function" returnParameter="componentRotationCurveGradient">
    !#  <functionArgs>node,radiusShared,massTypeShared,componentTypeShared</functionArgs>
    !#  <onReturn>Galactic_Structure_Rotation_Curve_Gradient=Galactic_Structure_Rotation_Curve_Gradient+componentRotationCurveGradient</onReturn>
    include 'galactic_structure.rotation_curve.gradient.tasks.inc'
    !# </include>

    ! Convert the summed d(V^2)/dr to dV/dr.
    Galactic_Structure_Rotation_Curve_Gradient=+0.5d0                                                                          &
         &                                     *Galactic_Structure_Rotation_Curve_Gradient                                     &
         &                                     /Galactic_Structure_Rotation_Curve         (node,radius,componentType,massType)
    return
  end function Galactic_Structure_Rotation_Curve_Gradient

  double precision function Component_Rotation_Curve_Gradient(component)
    !% Unary function returning the gradient of the squared rotation curve in a component. Suitable for mapping over components.
    use :: Galacticus_Nodes, only : nodeComponent
    implicit none
    class(nodeComponent), intent(inout) :: component

    Component_Rotation_Curve_Gradient=component%rotationCurveGradient(radiusShared,componentTypeShared,massTypeShared)
    return
  end function Component_Rotation_Curve_Gradient

end module Galactic_Structure_Rotation_Curve_Gradients












