!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
!!    Andrew Benson <abenson@carnegiescience.edu>
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

  !% Implementation of a simple freefall radius class.

  use Dark_Matter_Profiles
  
  !# <freefallRadius name="freefallRadiusDarkMatterHalo" defaultThreadPrivate="yes">
  !#  <description>
  !#   A freefall radius class which computes the freefall radius based on the freefall time in the dark matter halo.
  !#  </description>
  !# </freefallRadius>
  type, extends(freefallRadiusClass) :: freefallRadiusDarkMatterHalo
     !% Implementation of freefall radius class in which the freefall radius is based on the freefall time in the dark matter halo.
     private
     class(darkMatterProfileClass), pointer :: darkMatterProfile_
   contains
     final     ::                     darkMatterHaloDestructor
     procedure :: radius           => darkMatterHaloRadius
     procedure :: radiusGrowthRate => darkMatterHaloRadiusGrowthRate
  end type freefallRadiusDarkMatterHalo

  interface freefallRadiusDarkMatterHalo
     !% Constructors for the darkMatterHalo freefall radius class.
     module procedure darkMatterHaloConstructorParameters
     module procedure darkMatterHaloConstructorInternal
  end interface freefallRadiusDarkMatterHalo

contains

  function darkMatterHaloConstructorParameters(parameters) result(self)
    !% Constructor for the darkMatterHalo freefall radius class which builds the object from a parameter set.
    use Input_Parameters
    implicit none
    type (freefallRadiusDarkMatterHalo)                :: self
    type (inputParameters             ), intent(inout) :: parameters
    class(darkMatterProfileClass      ), pointer       :: darkMatterProfile_
    
    !# <objectBuilder class="darkMatterProfile" name="darkMatterProfile_" source="parameters"/>
    self=freefallRadiusDarkMatterHalo(darkMatterProfile_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function darkMatterHaloConstructorParameters

  function darkMatterHaloConstructorInternal(darkMatterProfile_) result(self)
    !% Internal constructor for the darkMatterHalo freefall radius class.
    implicit none
    type (freefallRadiusDarkMatterHalo)                        :: self
    class(darkMatterProfileClass      ), intent(in   ), target :: darkMatterProfile_
    !# <constructorAssign variables="*darkMatterProfile_"/>
    
    return
  end function darkMatterHaloConstructorInternal
  
  subroutine darkMatterHaloDestructor(self)
    !% Destructor for the darkMatterHalo freefall radius class.
    implicit none
    type(freefallRadiusDarkMatterHalo), intent(inout) :: self

    !# <objectDestructor name="self%darkMatterProfile_"/>
    return
  end subroutine darkMatterHaloDestructor

  double precision function darkMatterHaloRadiusGrowthRate(self,node)
    !% Returns the freefall radius growth rate (in Mpc/Gyr) in the hot atmosphere.
    use Cooling_Freefall_Times_Available
    implicit none
    class           (freefallRadiusDarkMatterHalo ), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    double precision                                               :: timeAvailable, timeAvailableIncreaseRate

    ! Get the time available for freefall.
    timeAvailable            =Cooling_Freefall_Time_Available              (node)
    ! Get the rate of increase of the time available for freefall.
    timeAvailableIncreaseRate=Cooling_Freefall_Time_Available_Increase_Rate(node)
    ! Get freefall radius increase rate from dark matter profile.
    darkMatterHaloRadiusGrowthRate=+self%darkMatterProfile_%freefallRadiusIncreaseRate(node,timeAvailable            ) &
         &                         *                                                        timeAvailableIncreaseRate
    return
  end function darkMatterHaloRadiusGrowthRate
  
  double precision function darkMatterHaloRadius(self,node)
    !% Return the freefall radius in the darkMatterHalo model.
    use Cooling_Freefall_Times_Available
    implicit none
    class           (freefallRadiusDarkMatterHalo ), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    double precision                                               :: timeAvailable

    ! Get the time available for freefall.
    timeAvailable       =Cooling_Freefall_Time_Available(node)
    ! Get freefall radius from dark matter profile.
    darkMatterHaloRadius=self%darkMatterProfile_%freefallRadius(node,timeAvailable)
    return
  end function darkMatterHaloRadius
