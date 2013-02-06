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

!% Contains a module which implements calculations of the gradient of the rotation curve.

module Galactic_Structure_Rotation_Curve_Gradients
  !% Implements calculations of the rotation curve gradient
  use, intrinsic :: ISO_C_Binding
  use ISO_Varying_String
  use Galacticus_Nodes
  use Galactic_Structure_Options
  private
  public :: Galactic_Structure_Rotation_Curve_Gradient

  ! Variables used in root finding.
  integer                  :: massTypeRoot,componentTypeRoot,weightByRoot,weightIndexRoot
  double precision         :: massRoot
  type(treeNode),  pointer :: activeNode
  !$omp threadprivate(massRoot,massTypeRoot,componentTypeRoot,weightByRoot,weightIndexRoot,activeNode)

contains

  double precision function Galactic_Structure_Rotation_Curve_Gradient(thisNode,radius,massType,componentType,haloLoaded)
    !% Solve for the rotation curve gradient at a given radius. Assumes the galactic structure has already been computed.
    use Galacticus_Error
    use Galactic_Structure_Rotation_Curves
    use Input_Parameters
    use Dark_Matter_Profiles
    use Numerical_Constants_Physical
    use Numerical_Constants_Prefixes
    !# <include directive="rotationCurveGradientTask" type="moduleUse">
    include 'galactic_structure.rotation_curve.gradient.tasks.modules.inc'
    !# </include>
    implicit none
    type(treeNode),   intent(inout), pointer  :: thisNode
    integer,          intent(in),    optional :: componentType,massType  
    logical,          intent(in),    optional :: haloLoaded
    double precision, intent(in)              :: radius
    integer                                   :: massTypeActual,componentTypeActual  
    double precision                          :: componentRotationCurveGradient

    ! Determine which component type to use.
    if (present(componentType)) then
       componentTypeActual=componentType
    else
       componentTypeActual=componentTypeAll
    end if
    ! Determine which mass type to use
    if (present(massType)) then
       massTypeActual=massType
    else
       massTypeActual=massTypeAll
    end if
    ! Initialize to zero gradient.
    Galactic_Structure_Rotation_Curve_Gradient=0.0d0
    
    ! Call routines to supply the gradient for all components' rotation curves. Specifically, the returned quantities are
    ! d(V^2)/dr so that they can be summed directly.
    !# <include directive="rotationCurveGradientTask" type="functionCall" functionType="function" returnParameter="componentRotationCurveGradient">
    !#  <functionArgs>thisNode,radius,massTypeActual,componentTypeActual,haloLoaded</functionArgs>
    !#  <onReturn>Galactic_Structure_Rotation_Curve_Gradient=Galactic_Structure_Rotation_Curve_Gradient+componentRotationCurveGradient</onReturn>
    include 'galactic_structure.rotation_curve.gradient.tasks.inc'
    !# </include>
  
    ! Convert the summed d(V^2)/dr to dV/dr.
    Galactic_Structure_Rotation_Curve_Gradient= 0.5d0                                                                              &
         &                                     *Galactic_Structure_Rotation_Curve_Gradient                                         &
         &                                     /Galactic_Structure_Rotation_Curve         (thisNode,radius,massType,componentType)

    return
  end function Galactic_Structure_Rotation_Curve_Gradient

end module Galactic_Structure_Rotation_Curve_Gradients
    
    
    
    








