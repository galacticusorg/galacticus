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


!% Contains a module which implements an isothermal (virial temperature) profile for hot gas halos.

module Hot_Halo_Temperature_Profile_Virial
  !% Implements an isothermal (virial temperature) profile for hot gas halos.
  private
  public :: Hot_Halo_Temperature_Virial

contains

  !# <hotHaloTemperatureMethod>
  !#  <unitName>Hot_Halo_Temperature_Virial</unitName>
  !# </hotHaloTemperatureMethod>
  subroutine Hot_Halo_Temperature_Virial(hotHaloTemperatureMethod,Hot_Halo_Temperature_Get,Hot_Halo_Temperature_Logarithmic_Slope_Get)
    !% Initialize the cored isothermal hot halo temperature profile module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),          intent(in)    :: hotHaloTemperatureMethod
    procedure(),          pointer, intent(inout) :: Hot_Halo_Temperature_Get,Hot_Halo_Temperature_Logarithmic_Slope_Get
    
    if (hotHaloTemperatureMethod == 'virial') then
       Hot_Halo_Temperature_Get => Hot_Halo_Temperature_Virial_Get
       Hot_Halo_Temperature_Logarithmic_Slope_Get => Hot_Halo_Temperature_Logarithmic_Slope_Virial_Get
    end if
    return
  end subroutine Hot_Halo_Temperature_Virial

  double precision function Hot_Halo_Temperature_Virial_Get(thisNode,radius)
    !% Compute the temperature at radius {\tt radius} in an isothermal (virial) temperature profile for {\tt thisNode}.
    use Tree_Nodes
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: radius

    Hot_Halo_Temperature_Virial_Get=Dark_Matter_Halo_Virial_Temperature(thisNode)
    return
  end function Hot_Halo_Temperature_Virial_Get
  
  double precision function Hot_Halo_Temperature_Logarithmic_Slope_Virial_Get(thisNode,radius)
    !% Compute the logarithmic slope of the temperature at radius {\tt radius} in an isothermal temperature profile
    !% for {\tt thisNode}.
    use Tree_Nodes
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: radius

    Hot_Halo_Temperature_Logarithmic_Slope_Virial_Get=0.0d0
    return
  end function Hot_Halo_Temperature_Logarithmic_Slope_Virial_Get
  
end module Hot_Halo_Temperature_Profile_Virial
