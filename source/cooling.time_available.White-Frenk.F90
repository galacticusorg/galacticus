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


!% Contains a module which implements the \cite{white_galaxy_1991} method for computing the time available for cooling in hot
!% halos.

module Cooling_Time_Available_White_Frenk
  !% Implements the \cite{white_galaxy_1991} method for computing the time available for cooling in hot halos.
  private
  public :: Cooling_Time_Available_WF_Initialize

  ! Parameter which interpolates between age of Universe and dynamical time for the time available for cooling.
  double precision :: coolingTimeAvailableAgeFactor

contains

  !# <coolingTimeAvailableMethod>
  !#  <unitName>Cooling_Time_Available_WF_Initialize</unitName>
  !# </coolingTimeAvailableMethod>
  subroutine Cooling_Time_Available_WF_Initialize(coolingTimeAvailableMethod,Cooling_Time_Available_Get&
       &,Cooling_Time_Available_Increase_Rate_Get)
    !% Initialize the \cite{white_galaxy_1991} cooling time available module.
    use ISO_Varying_String
    use Input_Parameters
    use Galacticus_Error
    implicit none
    type(varying_string),          intent(in)    :: coolingTimeAvailableMethod
    procedure(double precision), pointer, intent(inout) :: Cooling_Time_Available_Get,Cooling_Time_Available_Increase_Rate_Get
    
    if (coolingTimeAvailableMethod.eq.'White-Frenk') then
       Cooling_Time_Available_Get => Cooling_Time_Available_WF
       Cooling_Time_Available_Increase_Rate_Get => Cooling_Time_Available_Increase_Rate_WF
       !@ <inputParameter>
       !@   <name>coolingTimeAvailableAgeFactor</name>
       !@   <defaultValue>0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Interpolates (geometrically) between the age of the Universe and the halo dynamical time for the time available for cooling in the {\tt White-Frenk} method.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('coolingTimeAvailableAgeFactor',coolingTimeAvailableAgeFactor,defaultValue=0.0d0)
       if (coolingTimeAvailableAgeFactor < 0.0d0 .or. coolingTimeAvailableAgeFactor > 1.0d0) call Galacticus_Error_Report('Cooling_Time_Available_WF_Initialize','0 <= coolingTimeAvailableAgeFactor <= 1 is required')

    end if
    return
  end subroutine Cooling_Time_Available_WF_Initialize

  double precision function Cooling_Time_Available_WF(thisNode)
    !% Compute the time available for cooling using the \cite{white_galaxy_1991} method. This is assumed to be equal to the
    !% dynamical timescale of the halo.
    use Tree_Nodes
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    ! Return the appropriate time.
    if (coolingTimeAvailableAgeFactor == 1.0d0) then
       ! Time available equals the age of the Universe, which is just the time for this node.
       Cooling_Time_Available_WF=Tree_Node_Time(thisNode)
    else if (coolingTimeAvailableAgeFactor == 0.0d0) then
       ! Time available equals the halo dynamical time.
       Cooling_Time_Available_WF=Dark_Matter_Halo_Dynamical_Timescale(thisNode)
    else
       ! Time is interpolated between age of Universe and dynamical time. Do the interpolation.
       Cooling_Time_Available_WF=dexp(dlog(Tree_Node_Time(thisNode))*coolingTimeAvailableAgeFactor&
            &+dlog(Dark_Matter_Halo_Dynamical_Timescale(thisNode))*(1.0d0-coolingTimeAvailableAgeFactor))
    end if
    return
  end function Cooling_Time_Available_WF
  
  double precision function Cooling_Time_Available_Increase_Rate_WF(thisNode)
    !% Compute the rate of increase of the time available for cooling using the \cite{white_galaxy_1991} method. We return a rate
    !% of 1, even though technically it can depend on halo properties.
    use Tree_Nodes
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    ! Simply return unit rate.
    Cooling_Time_Available_Increase_Rate_WF=1.0d0
    return
  end function Cooling_Time_Available_Increase_Rate_WF
  
end module Cooling_Time_Available_White_Frenk
