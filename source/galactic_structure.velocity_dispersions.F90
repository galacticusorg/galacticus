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

!% Contains a module which implements calculations of the velocity dispersions in isotropic spherical systems by solving the Jeans
!% equation.

module Galactic_Structure_Velocity_Dispersions
  !% Implements calculations of the velocity dispersions in isotropic spherical systems by solving the Jeans equation.
  use :: Galacticus_Nodes, only : treeNode
  implicit none
  private
  public :: Galactic_Structure_Velocity_Dispersion

  ! Module scoped variables used in integrations.
  integer                    :: componentTypeGlobal, massTypeGlobal
  type   (treeNode), pointer :: activeNode
  !$omp threadprivate(massTypeGlobal,componentTypeGlobal,activeNode)
contains

  double precision function Galactic_Structure_Velocity_Dispersion(node,radius,radiusOuter,componentType,massType)
    !% Returns the velocity dispersion of the specified {\normalfont \ttfamily componentType} in {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius}.
    use :: Galactic_Structure_Densities      , only : Galactic_Structure_Density
    use :: Galactic_Structure_Enclosed_Masses, only : Galactic_Structure_Enclosed_Mass
    use :: Galactic_Structure_Options        , only : radiusLarge
    use :: Numerical_Integration             , only : integrator
    type            (treeNode  ), intent(inout), target   :: node
    double precision            , intent(in   )           :: radius          , radiusOuter
    integer                     , intent(in   )           :: componentType   , massType
    double precision                                      :: componentDensity, densityVelocityVariance, &
         &                                                   massTotal
    type            (integrator)                          :: integrator_

    activeNode         => node
    componentTypeGlobal=  componentType
    massTypeGlobal     =  massType
    ! Find the total mass.
    massTotal=Galactic_Structure_Enclosed_Mass(node,radiusLarge,componentType=componentType,massType=massType)
    ! Return with zero dispersion if the component is massless.
    if (massTotal <= 0.0d0) then
       Galactic_Structure_Velocity_Dispersion=0.0d0
       return
    end if
    ! Integrate the Jeans equation.
    integrator_            =integrator           (Velocity_Dispersion_Integrand,toleranceRelative=1.0d-3)
    densityVelocityVariance=integrator_%integrate(radius                       ,radiusOuter             )
    ! Get the density at this radius.
    componentDensity=Galactic_Structure_Density(node,[radius,0.0d0,0.0d0],componentType=componentType,massType=massType)
    ! Check for zero density.
    if (componentDensity <= 0.0d0) then
       Galactic_Structure_Velocity_Dispersion=0.0d0
    else
       Galactic_Structure_Velocity_Dispersion=sqrt(max(densityVelocityVariance,0.0d0)/componentDensity)
    end if
    return
  end function Galactic_Structure_Velocity_Dispersion

  double precision function Velocity_Dispersion_Integrand(radius)
    !% Integrand function used for finding velocity dispersions using Jeans equation.
    use :: Galactic_Structure_Densities      , only : Galactic_Structure_Density
    use :: Galactic_Structure_Enclosed_Masses, only : Galactic_Structure_Enclosed_Mass
    use :: Numerical_Constants_Astronomical      , only : gravitationalConstantGalacticus
    implicit none
    double precision, intent(in   ) :: radius

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
            &                                                          massType=massTypeGlobal            &
            &                                                         )
    end if
    return
  end function Velocity_Dispersion_Integrand

end module Galactic_Structure_Velocity_Dispersions
