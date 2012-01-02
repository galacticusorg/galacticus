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

!% Contains a module which implements calculations of the velocity dispersions in isotropic spherical systems by solving the Jeans
!% equation.

module Galactic_Structure_Velocity_Dispersions
  !% Implements calculations of the velocity dispersions in isotropic spherical systems by solving the Jeans equation.
  use Tree_Nodes
  implicit none
  private
  public :: Galactic_Structure_Velocity_Dispersion
  
  ! Module scoped variables used in integrations.
  integer                 :: componentTypeGlobal
  type(treeNode), pointer :: activeNode
  !$omp threadprivate(componentTypeGlobal,activeNode)

contains

  double precision function Galactic_Structure_Velocity_Dispersion(thisNode,radius,radiusOuter,componentType)
    !% Returns the velocity dispersion of the specified {\tt componentType} in {\tt thisNode} at the given {\tt radius}.
    use FGSL
    use, intrinsic :: ISO_C_Binding
    use Numerical_Integration
    use Galactic_Structure_Options
    use Galactic_Structure_Densities
    type(treeNode),   intent(inout), pointer :: thisNode  
    double precision, intent(in)             :: radius,radiusOuter
    integer,          intent(in)             :: componentType
    double precision                         :: densityVelocityVariance,componentDensity
    type(c_ptr)                              :: parameterPointer
    type(fgsl_function)                      :: integrandFunction
    type(fgsl_integration_workspace)         :: integrationWorkspace

    activeNode         => thisNode
    componentTypeGlobal=  componentType
    densityVelocityVariance=Integrate(radius,radiusOuter,Velocity_Dispersion_Integrand,parameterPointer&
         &,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
    call Integrate_Done(integrandFunction,integrationWorkspace)

    ! Get the density at this radius.
    componentDensity=Galactic_Structure_Density(thisNode,[radius,0.0d0,0.0d0],componentType=componentType)

    ! Check for zero density.
    if (componentDensity <= 0.0d0) then
       Galactic_Structure_Velocity_Dispersion=0.0d0
    else
       Galactic_Structure_Velocity_Dispersion=sqrt(densityVelocityVariance/componentDensity)
    end if
    return
  end function Galactic_Structure_Velocity_Dispersion

  function Velocity_Dispersion_Integrand(radius,parameterPointer) bind(c)
    !% Integrand function used for finding velocity dispersions using Jeans equation.
    use, intrinsic :: ISO_C_Binding
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    use Galactic_Structure_Densities
    use Galactic_Structure_Enclosed_Masses
    implicit none
    real(c_double)        :: Velocity_Dispersion_Integrand
    real(c_double), value :: radius
    type(c_ptr),    value :: parameterPointer

    Velocity_Dispersion_Integrand= gravitationalConstantGalacticus                                                                     &
         &                        *Galactic_Structure_Enclosed_Mass(activeNode, radius                                               ) &
         &                        /radius**2                                                                                           &
         &                        *Galactic_Structure_Density      (activeNode,[radius,0.0d0,0.0d0],componentType=componentTypeGlobal)
    return
  end function Velocity_Dispersion_Integrand

end module Galactic_Structure_Velocity_Dispersions
