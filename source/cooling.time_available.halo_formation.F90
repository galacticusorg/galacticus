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

!% Contains a module which implements the \cite{cole_hierarchical_2000} method for computing the time available for cooling in hot
!% halos.

module Cooling_Times_Available_Halo_Formation
  !% Implements the \cite{cole_hierarchical_2000} method for computing the time available for cooling in hot halos.
  use Galacticus_Nodes
  implicit none
  private
  public :: Cooling_Time_Available_Halo_Formation_Initialize

contains

  !# <coolingTimeAvailableMethod>
  !#  <unitName>Cooling_Time_Available_Halo_Formation_Initialize</unitName>
  !# </coolingTimeAvailableMethod>
  subroutine Cooling_Time_Available_Halo_Formation_Initialize(coolingTimeAvailableMethod,Cooling_Time_Available_Get&
       &,Cooling_Time_Available_Increase_Rate_Get)
    !% Initialize the \cite{cole_hierarchical_2000} cooling time available module.
    use ISO_Varying_String
    use Input_Parameters
    use Galacticus_Error
    implicit none
    type     (varying_string                                     ), intent(in   )          :: coolingTimeAvailableMethod
    procedure(Cooling_Time_Available_Halo_Formation              ), intent(inout), pointer :: Cooling_Time_Available_Get
    procedure(Cooling_Time_Available_Increase_Rate_Halo_Formation), intent(inout), pointer :: Cooling_Time_Available_Increase_Rate_Get

    if (coolingTimeAvailableMethod == 'haloFormation') then
       ! Set pointers to our implementation.
       Cooling_Time_Available_Get               => Cooling_Time_Available_Halo_Formation
       Cooling_Time_Available_Increase_Rate_Get => Cooling_Time_Available_Increase_Rate_Halo_Formation
       ! Check that there is a gettable formation time property.
       if (.not.defaultFormationTimeComponent%formationTimeIsGettable()) call Galacticus_Error_Report('Cooling_Time_Available_Halo_Formation_Initialize',"'haloFormation' method for time available for cooling requires a formationTime component that supports getting of the formationTime property")
    end if
    return
  end subroutine Cooling_Time_Available_Halo_Formation_Initialize

  double precision function Cooling_Time_Available_Halo_Formation(thisNode)
    !% Compute the time available for cooling using the \cite{cole_hierarchical_2000} method. Specifically, the time available is
    !% assumed to be the time since the halo formation event.
    implicit none
    type (treeNode                  ), intent(inout), pointer :: thisNode
    class(nodeComponentBasic        )               , pointer :: thisBasicComponent
    class(nodeComponentFormationTime)               , pointer :: thisFormationTimeComponent

    thisBasicComponent         => thisNode%basic        ()
    thisFormationTimeComponent => thisNode%formationTime()

    Cooling_Time_Available_Halo_Formation=thisBasicComponent%time()-thisFormationTimeComponent%formationTime()
    return
  end function Cooling_Time_Available_Halo_Formation

  double precision function Cooling_Time_Available_Increase_Rate_Halo_Formation(thisNode)
    !% Compute the rate of increase of the time available for cooling using the \cite{cole_hierarchical_2000} method. We return a rate
    !% of 1.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    ! Simply return unit rate.
    Cooling_Time_Available_Increase_Rate_Halo_Formation=1.0d0
    return
  end function Cooling_Time_Available_Increase_Rate_Halo_Formation

end module Cooling_Times_Available_Halo_Formation
