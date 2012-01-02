!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements the \cite{cole_hierarchical_2000} method for computing the time available for freefall in
!% cooling calculations in hot halos.

module Freefall_Times_Available_Halo_Formation
  !% Implements the \cite{cole_hierarchical_2000} method for computing the time available for freefall in
  !% cooling calculations in hot halos.
  implicit none
  private
  public :: Freefall_Time_Available_Halo_Formation_Initialize

contains

  !# <freefallTimeAvailableMethod>
  !#  <unitName>Freefall_Time_Available_Halo_Formation_Initialize</unitName>
  !# </freefallTimeAvailableMethod>
  subroutine Freefall_Time_Available_Halo_Formation_Initialize(freefallTimeAvailableMethod,Freefall_Time_Available_Get&
       &,Freefall_Time_Available_Increase_Rate_Get)
    !% Initialize the \cite{cole_hierarchical_2000} freefall time available module.
    use ISO_Varying_String
    use Input_Parameters
    use Galacticus_Error
    use Tree_Nodes
    implicit none
    type(varying_string),                 intent(in)    :: freefallTimeAvailableMethod
    procedure(double precision), pointer, intent(inout) :: Freefall_Time_Available_Get,Freefall_Time_Available_Increase_Rate_Get
    
    if (freefallTimeAvailableMethod == 'haloFormation') then
       ! Set pointers to our implementation.
       Freefall_Time_Available_Get               => Freefall_Time_Available_Halo_Formation
       Freefall_Time_Available_Increase_Rate_Get => Freefall_Time_Available_Increase_Rate_Halo_Formation
    end if
    return
  end subroutine Freefall_Time_Available_Halo_Formation_Initialize

  double precision function Freefall_Time_Available_Halo_Formation(thisNode)
    !% Compute the time available for freefall using the \cite{cole_hierarchical_2000} method. Specifically, the time available is
    !% assumed to be the time since the halo formation event. This function expects that the node passed in will be a formation
    !% node, such that its parent node gives the active node to which it is attached.
    use Tree_Nodes
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    Freefall_Time_Available_Halo_Formation=Tree_Node_Time(thisNode%parentNode)-Tree_Node_Time(thisNode)
    return
  end function Freefall_Time_Available_Halo_Formation
  
  double precision function Freefall_Time_Available_Increase_Rate_Halo_Formation(thisNode)
    !% Compute the rate of increase of the time available for freefall using the \cite{cole_hierarchical_2000} method. We return a rate
    !% of 1.
    use Tree_Nodes
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    ! Simply return unit rate.
    Freefall_Time_Available_Increase_Rate_Halo_Formation=1.0d0
    return
  end function Freefall_Time_Available_Increase_Rate_Halo_Formation
  
end module Freefall_Times_Available_Halo_Formation
