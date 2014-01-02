!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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
  integer                    :: componentTypeGlobal, massTypeGlobal
  type   (treeNode), pointer :: activeNode
  logical                    :: haloLoadedActual
  !$omp threadprivate(massTypeGlobal,componentTypeGlobal,activeNode,haloLoadedActual)
contains

  double precision function Galactic_Structure_Velocity_Dispersion(thisNode,radius,radiusOuter,componentType,massType,haloLoaded)
    !% Returns the velocity dispersion of the specified {\tt componentType} in {\tt thisNode} at the given {\tt radius}.
    use FGSL
    use, intrinsic :: ISO_C_Binding
    use Numerical_Integration
    use Galactic_Structure_Options
    use Galactic_Structure_Densities
    use Galactic_Structure_Enclosed_Masses
    type            (treeNode                  ), intent(inout), pointer  :: thisNode
    double precision                            , intent(in   )           :: radius              , radiusOuter
    integer                                     , intent(in   )           :: componentType       , massType
    logical                                     , intent(in   ), optional :: haloLoaded
    double precision                                                      :: componentDensity    , densityVelocityVariance, &
         &                                                                   massTotal
    type            (c_ptr                     )                          :: parameterPointer
    type            (fgsl_function             )                          :: integrandFunction
    type            (fgsl_integration_workspace)                          :: integrationWorkspace

    activeNode         => thisNode
    componentTypeGlobal=  componentType
    massTypeGlobal     =  massType
    haloLoadedActual   =  .true.
    if (present(haloLoaded)) haloLoadedActual=haloLoaded
    ! Find the total mass.
    massTotal=Galactic_Structure_Enclosed_Mass(thisNode,radiusLarge,componentType=componentType,massType=massType)
    ! Return with zero dispersion if the component is massless.
    if (massTotal <= 0.0d0) then
       Galactic_Structure_Velocity_Dispersion=0.0d0
       return
    end if
    ! Integrate the Jeans equation.
    densityVelocityVariance=Integrate(radius,radiusOuter,Velocity_Dispersion_Integrand,parameterPointer&
         &,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
    call Integrate_Done(integrandFunction,integrationWorkspace)
    ! Get the density at this radius.
    componentDensity=Galactic_Structure_Density(thisNode,[radius,0.0d0,0.0d0],componentType=componentType,massType=massType)
    ! Check for zero density.
    if (componentDensity <= 0.0d0) then
       Galactic_Structure_Velocity_Dispersion=0.0d0
    else
       Galactic_Structure_Velocity_Dispersion=sqrt(max(densityVelocityVariance,0.0d0)/componentDensity)
    end if
    return
  end function Galactic_Structure_Velocity_Dispersion

  function Velocity_Dispersion_Integrand(radius,parameterPointer) bind(c)
    !% Integrand function used for finding velocity dispersions using Jeans equation.
    use, intrinsic :: ISO_C_Binding
    use Numerical_Constants_Physical
    use Galactic_Structure_Densities
    use Galactic_Structure_Enclosed_Masses
    implicit none
    real(kind=c_double)        :: Velocity_Dispersion_Integrand
    real(kind=c_double), value :: radius
    type(c_ptr        ), value :: parameterPointer

    if (radius == 0.0d0) then
       Velocity_Dispersion_Integrand=0.0d0
    else
       Velocity_Dispersion_Integrand= gravitationalConstantGalacticus                                     &
            &                        *Galactic_Structure_Enclosed_Mass(                                   &
            &                                                          activeNode                       , &
            &                                                          radius                             &
            &                                                         )                                   &
            &                        /radius**2                                                           &
            &                        *Galactic_Structure_Density      (                                   &
            &                                                          activeNode                       , &
            &                                                          [radius,0.0d0,0.0d0]             , &
            &                                                          componentType=componentTypeGlobal, &
            &                                                          massType=massTypeGlobal          , &
            &                                                          haloLoaded=haloLoadedActual        &
            &                                                         )
    end if
    return
  end function Velocity_Dispersion_Integrand

end module Galactic_Structure_Velocity_Dispersions
