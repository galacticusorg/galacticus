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

!% Contains a module which implements a simple infall radius calculation, simply assuming that the infall radius equals the cooling radius.

module Infall_Radii_Cooling_Radius
  !% Implements a simple infall radius calculation, simply assuming that the infall radius equals the cooling radius.
  implicit none
  private
  public :: Infall_Radius_Cooling_Radius_Initialize

contains

  !# <infallRadiusMethod>
  !#  <unitName>Infall_Radius_Cooling_Radius_Initialize</unitName>
  !# </infallRadiusMethod>
  subroutine Infall_Radius_Cooling_Radius_Initialize(infallRadiusMethod,Infall_Radius_Get,Infall_Radius_Growth_Rate_Get)
    !% Initializes the ``cooling radius'' infall radius module.
    use ISO_Varying_String
    use Abundances_Structure
    use Chemical_Abundances_Structure
    implicit none
    type     (varying_string                          ), intent(in   )          :: infallRadiusMethod             
    procedure(Infall_Radius_Cooling_Radius            ), intent(inout), pointer :: Infall_Radius_Get              
    procedure(Infall_Radius_Growth_Rate_Cooling_Radius), intent(inout), pointer :: Infall_Radius_Growth_Rate_Get  
                                                                                                               
    if (infallRadiusMethod == 'coolingRadius') then
       Infall_Radius_Get             => Infall_Radius_Cooling_Radius
       Infall_Radius_Growth_Rate_Get => Infall_Radius_Growth_Rate_Cooling_Radius
    end if
    return
  end subroutine Infall_Radius_Cooling_Radius_Initialize

  double precision function Infall_Radius_Cooling_Radius(thisNode)
    !% Return the growth rate of the infall radius in the ``cooling radius'' model in Mpc/Gyr.
    use Galacticus_Nodes
    use Cooling_Radii
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode  
                                                     
    Infall_Radius_Cooling_Radius=Cooling_Radius(thisNode)
    return
  end function Infall_Radius_Cooling_Radius
  
  double precision function Infall_Radius_Growth_Rate_Cooling_Radius(thisNode)
    !% Return the growth rate of the infall radius in the ``cooling radius'' model in Mpc/Gyr.
    use Galacticus_Nodes
    use Cooling_Radii
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode  
                                                     
    Infall_Radius_Growth_Rate_Cooling_Radius=Cooling_Radius_Growth_Rate(thisNode)
    return
  end function Infall_Radius_Growth_Rate_Cooling_Radius
  
end module Infall_Radii_Cooling_Radius
