!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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






!% Contains a module which implements a \cite{white_galaxy_1991} cooling rate calculation.

module Cooling_Rates_White_Frenk
  !% Implements a \cite{white_galaxy_1991} cooling rate calculation.
  use, intrinsic :: ISO_C_Binding
  use Tree_Nodes
  private
  public :: Cooling_Rate_White_Frenk_Initialize

contains

  !# <coolingRateMethod>
  !#  <unitName>Cooling_Rate_White_Frenk_Initialize</unitName>
  !# </coolingRateMethod>
  subroutine Cooling_Rate_White_Frenk_Initialize(coolingRateMethod,Cooling_Rate_Get)
    !% Initializes the ``White + Frenk'' cooling rate module.
    use ISO_Varying_String
    implicit none
    type(varying_string),          intent(in)    :: coolingRateMethod
    procedure(),          pointer, intent(inout) :: Cooling_Rate_Get
    
    if (coolingRateMethod == 'White + Frenk') Cooling_Rate_Get => Cooling_Rate_White_Frenk
    return
  end subroutine Cooling_Rate_White_Frenk_Initialize

  double precision function Cooling_Rate_White_Frenk(thisNode)
    !% Computes the mass cooling rate in a hot gas halo utilizing the \cite{white_galaxy_1991} method.
    use Tree_Nodes
    use Tree_Node_Methods
    use Dark_Matter_Halo_Scales
    use Cooling_Times_Available
    use Cooling_Radii
    use Numerical_Constants_Math
    use Hot_Halo_Density_Profile
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision                         :: coolingRadius,coolingDensity,virialRadius,coolingRadiusGrowthRate

    ! Get the virial radius.
    virialRadius=Dark_Matter_Halo_Virial_Radius(thisNode)

    ! Get the cooling radius.
    coolingRadius=Cooling_Radius(thisNode)

    if (coolingRadius >= virialRadius) then
       ! Cooling radius exceeds the virial radius. Limit infall to the dynamical timescale.
       Cooling_Rate_White_Frenk=Tree_Node_Hot_Halo_Mass(thisNode)/Dark_Matter_Halo_Dynamical_Timescale(thisNode)
    else
       ! Find the density at the cooling radius.
       coolingDensity=Hot_Halo_Density(thisNode,coolingRadius)
       ! Find cooling radius growth rate.
       coolingRadiusGrowthRate=Cooling_Radius_Growth_Rate(thisNode)
       ! Compute the cooling rate.
       Cooling_Rate_White_Frenk=4.0d0*Pi*(coolingRadius**2)*coolingDensity*coolingRadiusGrowthRate
    end if
    return
  end function Cooling_Rate_White_Frenk

end module Cooling_Rates_White_Frenk
