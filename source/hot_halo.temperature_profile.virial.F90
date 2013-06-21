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

!% Contains a module which implements an isothermal (virial temperature) profile for hot gas halos.

module Hot_Halo_Temperature_Profile_Virial
  !% Implements an isothermal (virial temperature) profile for hot gas halos.
  implicit none
  private
  public :: Hot_Halo_Temperature_Virial

contains

  !# <hotHaloTemperatureMethod>
  !#  <unitName>Hot_Halo_Temperature_Virial</unitName>
  !# </hotHaloTemperatureMethod>
  subroutine Hot_Halo_Temperature_Virial(hotHaloTemperatureMethod,Hot_Halo_Temperature_Get,Hot_Halo_Temperature_Logarithmic_Slope_Get)
    !% Initialize the cored isothermal hot halo temperature profile module.
    use ISO_Varying_String
    implicit none
    type     (varying_string                                   ), intent(in   )          :: hotHaloTemperatureMethod
    procedure(Hot_Halo_Temperature_Virial_Get                  ), intent(inout), pointer :: Hot_Halo_Temperature_Get
    procedure(Hot_Halo_Temperature_Logarithmic_Slope_Virial_Get), intent(inout), pointer :: Hot_Halo_Temperature_Logarithmic_Slope_Get

    if (hotHaloTemperatureMethod == 'virial') then
       Hot_Halo_Temperature_Get                   => Hot_Halo_Temperature_Virial_Get
       Hot_Halo_Temperature_Logarithmic_Slope_Get => Hot_Halo_Temperature_Logarithmic_Slope_Virial_Get
    end if
    return
  end subroutine Hot_Halo_Temperature_Virial

  double precision function Hot_Halo_Temperature_Virial_Get(thisNode,radius)
    !% Compute the temperature at radius {\tt radius} in an isothermal (virial) temperature profile for {\tt thisNode}.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode
    double precision          , intent(in   )          :: radius

    Hot_Halo_Temperature_Virial_Get=Dark_Matter_Halo_Virial_Temperature(thisNode)
    return
  end function Hot_Halo_Temperature_Virial_Get

  double precision function Hot_Halo_Temperature_Logarithmic_Slope_Virial_Get(thisNode,radius)
    !% Compute the logarithmic slope of the temperature at radius {\tt radius} in an isothermal temperature profile
    !% for {\tt thisNode}.
    use Galacticus_Nodes
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode
    double precision          , intent(in   )          :: radius

    Hot_Halo_Temperature_Logarithmic_Slope_Virial_Get=0.0d0
    return
  end function Hot_Halo_Temperature_Logarithmic_Slope_Virial_Get

end module Hot_Halo_Temperature_Profile_Virial
