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

!+    Contributions to this file made by:  Stéphane Mangeon, Andrew Benson.

!% Contains a module which implements calculations of the velocity dispersions in isotropic spherical systems by solving the Jeans
!% equation.

module Galactic_Structure_Velocity_Dispersions
  !% Implements calculations of the velocity dispersions in isotropic spherical systems by solving the Jeans equation.
  use Galacticus_Nodes
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
