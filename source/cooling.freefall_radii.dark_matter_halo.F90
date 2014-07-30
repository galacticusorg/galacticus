!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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
    type     (varying_string                              ), intent(in   )          :: freefallRadiusMethod
    procedure(Freefall_Radius_Dark_Matter_Halo            ), intent(inout), pointer :: Freefall_Radius_Get
    procedure(Freefall_Radius_Growth_Rate_Dark_Matter_Halo), intent(inout), pointer :: Freefall_Radius_Growth_Rate_Get

    if (freefallRadiusMethod == 'darkMatterHalo') then
       Freefall_Radius_Get             => Freefall_Radius_Dark_Matter_Halo
       Freefall_Radius_Growth_Rate_Get => Freefall_Radius_Growth_Rate_Dark_Matter_Halo
    end if
    return
  end subroutine Freefall_Radius_Dark_Matter_Halo_Initialize

  double precision function Freefall_Radius_Growth_Rate_Dark_Matter_Halo(thisNode)
    !% Return the growth rate of the freefall radius in the ``dark matter halo'' model in Mpc/Gyr.
    use Galacticus_Nodes
    use Cooling_Freefall_Times_Available
    use Dark_Matter_Profiles
    implicit none
    type            (treeNode              ), intent(inout), pointer :: thisNode
    class           (darkMatterProfileClass)               , pointer :: darkMatterProfile_
    double precision                                                 :: timeAvailable     , timeAvailableIncreaseRate

    ! Get required objects.
    darkMatterProfile_ => darkMatterProfile()

    ! Get the time available for freefall.
    timeAvailable=Cooling_Freefall_Time_Available(thisNode)

    ! Get the rate of increase of the time available for freefall.
    timeAvailableIncreaseRate=Cooling_Freefall_Time_Available_Increase_Rate(thisNode)

    ! Get freefall radius increase rate from dark matter profile.
    Freefall_Radius_Growth_Rate_Dark_Matter_Halo=darkMatterProfile_%freefallRadiusIncreaseRate(thisNode,timeAvailable)&
         &*timeAvailableIncreaseRate
    return
  end function Freefall_Radius_Growth_Rate_Dark_Matter_Halo

  double precision function Freefall_Radius_Dark_Matter_Halo(thisNode)
    !% Return the freefall radius in the ``dark matter halo'' model.
    use Galacticus_Nodes
    use Cooling_Freefall_Times_Available
    use Dark_Matter_Profiles
    implicit none
    type            (treeNode              ), intent(inout), pointer :: thisNode
    class           (darkMatterProfileClass)               , pointer :: darkMatterProfile_
    double precision                                                 :: timeAvailable

    ! Get required objects.
    darkMatterProfile_ => darkMatterProfile()

    ! Get the time available for freefall.
    timeAvailable=Cooling_Freefall_Time_Available(thisNode)

    ! Get freefall radius from dark matter profile.
    Freefall_Radius_Dark_Matter_Halo=darkMatterProfile_%freefallRadius(thisNode,timeAvailable)
    return
  end function Freefall_Radius_Dark_Matter_Halo

end module Freefall_Radii_Dark_Matter_Halo
