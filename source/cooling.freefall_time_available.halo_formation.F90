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
       ! Check that there is a gettable formation time property.
       if (.not.associated(Tree_Node_Formation_Time)) call Galacticus_Error_Report('Freefall_Time_Available_Halo_Formation_Initialize',"'haloFormation' method for time available for freefall requires a component that supports getting of the Formation_Time property")
    end if
    return
  end subroutine Freefall_Time_Available_Halo_Formation_Initialize

  double precision function Freefall_Time_Available_Halo_Formation(thisNode)
    !% Compute the time available for freefall using the \cite{cole_hierarchical_2000} method. Specifically, the time available is
    !% assumed to be the time since the halo formation event.
    use Tree_Nodes
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    Freefall_Time_Available_Halo_Formation=Tree_Node_Time(thisNode)-Tree_Node_Formation_Time(thisNode)
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
