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

!% Contains a module which implements an infall radius calculation in which the infall radius is the smaller of the cooling and
!% freefall radii.

module Infall_Radii_Cooling_Freefall
  !% Implements an infall radius calculation in which the infall radius is the smaller of the cooling and freefall radii.
  implicit none
  private
  public :: Infall_Radius_Cooling_Freefall_Initialize

contains

  !# <infallRadiusMethod>
  !#  <unitName>Infall_Radius_Cooling_Freefall_Initialize</unitName>
  !# </infallRadiusMethod>
  subroutine Infall_Radius_Cooling_Freefall_Initialize(infallRadiusMethod,Infall_Radius_Get,Infall_Radius_Growth_Rate_Get)
    !% Initializes the ``cooling and freefall'' infall radius module.
    use ISO_Varying_String
    use Abundances_Structure
    use Chemical_Abundances_Structure
    implicit none
    type(varying_string),                 intent(in)    :: infallRadiusMethod
    procedure(double precision), pointer, intent(inout) :: Infall_Radius_Get,Infall_Radius_Growth_Rate_Get
    
    if (infallRadiusMethod == 'coolingAndFreefall') then
       Infall_Radius_Get             => Infall_Radius_Cooling_Freefall
       Infall_Radius_Growth_Rate_Get => Infall_Radius_Growth_Rate_Cooling_Freefall
    end if
    return
  end subroutine Infall_Radius_Cooling_Freefall_Initialize

  double precision function Infall_Radius_Cooling_Freefall(thisNode)
    !% Return the infall radius in the ``cooling and freefall'' model in Mpc/Gyr.
    use Galacticus_Nodes
    use Cooling_Radii
    use Freefall_Radii
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision                         :: virialRadius,coolingRadius,freefallRadius
    logical                                  :: infallIsCoolingLimited    

    ! Get the virial radius.
    virialRadius  =Dark_Matter_Halo_Virial_Radius(thisNode)

    ! Get the cooling radius.
    coolingRadius =Cooling_Radius                (thisNode)

    ! Get the freefall radius.
    freefallRadius=Freefall_Radius               (thisNode)

    ! Compute the infall radius as the smaller of the cooling and freefall radii.
    infallIsCoolingLimited=(coolingRadius < freefallRadius)
    if (infallIsCoolingLimited) then
       Infall_Radius_Cooling_Freefall=coolingRadius
    else
       Infall_Radius_Cooling_Freefall=freefallRadius
    end if
    return
  end function Infall_Radius_Cooling_Freefall
  
  double precision function Infall_Radius_Growth_Rate_Cooling_Freefall(thisNode)
    !% Return the growth rate of the infall radius in the ``cooling and freefall'' model in Mpc/Gyr.
    use Galacticus_Nodes
    use Cooling_Radii
    use Freefall_Radii
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision                         :: virialRadius,coolingRadius,freefallRadius
    logical                                  :: infallIsCoolingLimited    

    ! Get the virial radius.
    virialRadius  =Dark_Matter_Halo_Virial_Radius(thisNode)

    ! Get the cooling radius.
    coolingRadius =Cooling_Radius                (thisNode)

    ! Get the freefall radius.
    freefallRadius=Freefall_Radius               (thisNode)

    ! Compute the infall radius as the smaller of the cooling and freefall radii.
    infallIsCoolingLimited=(coolingRadius < freefallRadius)
    if (infallIsCoolingLimited) then
       Infall_Radius_Growth_Rate_Cooling_Freefall=Cooling_Radius_Growth_Rate (thisNode)
    else
       Infall_Radius_Growth_Rate_Cooling_Freefall=Freefall_Radius_Growth_Rate(thisNode)
    end if
    return
  end function Infall_Radius_Growth_Rate_Cooling_Freefall
  
end module Infall_Radii_Cooling_Freefall
