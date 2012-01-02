!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!+    Contributions to this file made by:  Stéphane Mangeon, Andrew Benson.

!% Contains a module which implements calculations of the gradient of the rotation curve.

module Galactic_Structure_Rotation_Curve_Gradients
  !% Implements calculations of the rotation curve gradient
  use, intrinsic :: ISO_C_Binding
  use ISO_Varying_String
  use Tree_Nodes
  use Galactic_Structure_Options
  private
  public :: Galactic_Structure_Rotation_Curve_Gradient

  ! Variables used in root finding.
  integer                  :: massTypeRoot,componentTypeRoot,weightByRoot,weightIndexRoot
  double precision         :: massRoot
  type(treeNode),  pointer :: activeNode
  !$omp threadprivate(massRoot,massTypeRoot,componentTypeRoot,weightByRoot,weightIndexRoot,activeNode)

contains

  double precision function Galactic_Structure_Rotation_Curve_Gradient(thisNode,radius,massType,componentType)
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
    !# <include directive="rotationCurveGradientTask" type="code" action="subroutine">
    !#  <subroutineArgs>thisNode,radius,massTypeActual,componentTypeActual,componentRotationCurveGradient</subroutineArgs>
    !#  <subroutineAction>Galactic_Structure_Rotation_Curve_Gradient=Galactic_Structure_Rotation_Curve_Gradient+componentRotationCurveGradient</subroutineAction>
    include 'galactic_structure.rotation_curve.gradient.tasks.inc'
    !# </include>
  
    ! Convert the summed d(V^2)/dr to dV/dr.
    Galactic_Structure_Rotation_Curve_Gradient= 0.5d0                                                                              &
         &                                     *Galactic_Structure_Rotation_Curve_Gradient                                         &
         &                                     /Galactic_Structure_Rotation_Curve         (thisNode,radius,massType,componentType)

    return
  end function Galactic_Structure_Rotation_Curve_Gradient

end module Galactic_Structure_Rotation_Curve_Gradients
    
    
    
    








