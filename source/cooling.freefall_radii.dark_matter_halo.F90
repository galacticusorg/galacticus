!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements a simple cooling radius calculation (finds the radius at which the time available for
!% cooling equals the cooling time).

module Freefall_Radii_Dark_Matter_Halo
  !% Implements a simple cooling radius calculation (finds the radius at which the time available for cooling equals the cooling
  !% time).
  implicit none
  private
  public :: Freefall_Radius_Dark_Matter_Halo_Initialize

contains

  !# <freefallRadiusMethod>
  !#  <unitName>Freefall_Radius_Dark_Matter_Halo_Initialize</unitName>
  !# </freefallRadiusMethod>
  subroutine Freefall_Radius_Dark_Matter_Halo_Initialize(freefallRadiusMethod,Freefall_Radius_Get,Freefall_Radius_Growth_Rate_Get)
    !% Initializes the ``darkMatterHalo'' freefall radius module.
    use ISO_Varying_String
    implicit none
    type(varying_string),                 intent(in)    :: freefallRadiusMethod
    procedure(double precision), pointer, intent(inout) :: Freefall_Radius_Get,Freefall_Radius_Growth_Rate_Get
    
    if (freefallRadiusMethod == 'darkMatterHalo') then
       Freefall_Radius_Get             => Freefall_Radius_Dark_Matter_Halo
       Freefall_Radius_Growth_Rate_Get => Freefall_Radius_Growth_Rate_Dark_Matter_Halo
    end if
    return
  end subroutine Freefall_Radius_Dark_Matter_Halo_Initialize

  double precision function Freefall_Radius_Growth_Rate_Dark_Matter_Halo(thisNode)
    !% Return the growth rate of the freefall radius in the ``dark matter halo'' model in Mpc/Gyr.
    use Tree_Nodes
    use Cooling_Freefall_Times_Available
    use Dark_Matter_Profiles
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision                         :: timeAvailable,timeAvailableIncreaseRate

    ! Get the time available for freefall.
    timeAvailable=Cooling_Freefall_Time_Available(thisNode)
    
    ! Get the rate of increase of the time available for freefall.
    timeAvailableIncreaseRate=Cooling_Freefall_Time_Available_Increase_Rate(thisNode)

    ! Get freefall radius increase rate from dark matter profile.
    Freefall_Radius_Growth_Rate_Dark_Matter_Halo=Dark_Matter_Profile_Freefall_Radius_Increase_Rate(thisNode,timeAvailable)&
         &*timeAvailableIncreaseRate
    return
  end function Freefall_Radius_Growth_Rate_Dark_Matter_Halo

  double precision function Freefall_Radius_Dark_Matter_Halo(thisNode)
    !% Return the freefall radius in the ``dark matter halo'' model.
    use Tree_Nodes
    use Cooling_Freefall_Times_Available
    use Dark_Matter_Profiles
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision                         :: timeAvailable

    ! Get the time available for freefall.
    timeAvailable=Cooling_Freefall_Time_Available(thisNode)
    
    ! Get freefall radius from dark matter profile.
    Freefall_Radius_Dark_Matter_Halo=Dark_Matter_Profile_Freefall_Radius(thisNode,timeAvailable)
    return
  end function Freefall_Radius_Dark_Matter_Halo

end module Freefall_Radii_Dark_Matter_Halo
