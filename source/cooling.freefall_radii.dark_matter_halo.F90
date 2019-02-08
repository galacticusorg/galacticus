!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  use Dark_Matter_Profiles            , only : darkMatterProfileClass    , darkMatterProfile
  use Cooling_Freefall_Times_Available, only : freefallTimeAvailableClass, freefallTimeAvailable

  !# <freefallRadius name="freefallRadiusDarkMatterHalo">
  !#  <description>
  !#   A freefall radius class which computes the freefall radius based on the freefall time in the dark matter halo.
  !#  </description>
  !# </freefallRadius>
  type, extends(freefallRadiusClass) :: freefallRadiusDarkMatterHalo
     !% Implementation of freefall radius class in which the freefall radius is based on the freefall time in the dark matter halo.
     private
     class(darkMatterProfileClass    ), pointer :: darkMatterProfile_ => null()
     class(freefallTimeAvailableClass), pointer :: freefallTimeAvailable_ => null()
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
    class(freefallTimeAvailableClass  ), pointer       :: freefallTimeAvailable_
    
    !# <objectBuilder class="darkMatterProfile"     name="darkMatterProfile_"     source="parameters"/>
    !# <objectBuilder class="freefallTimeAvailable" name="freefallTimeAvailable_" source="parameters"/>
    self=freefallRadiusDarkMatterHalo(darkMatterProfile_,freefallTimeAvailable_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="darkMatterProfile_"    />
    !# <objectDestructor name="freefallTimeAvailable_"/>
    return
  end function darkMatterHaloConstructorParameters

  function darkMatterHaloConstructorInternal(darkMatterProfile_,freefallTimeAvailable_) result(self)
    !% Internal constructor for the darkMatterHalo freefall radius class.
    implicit none
    type (freefallRadiusDarkMatterHalo)                        :: self
    class(darkMatterProfileClass      ), intent(in   ), target :: darkMatterProfile_
    class(freefallTimeAvailableClass  ), intent(in   ), target :: freefallTimeAvailable_
    !# <constructorAssign variables="*darkMatterProfile_, *freefallTimeAvailable_"/>
    
    return
  end function darkMatterHaloConstructorInternal
  
  subroutine darkMatterHaloDestructor(self)
    !% Destructor for the darkMatterHalo freefall radius class.
    implicit none
    type(freefallRadiusDarkMatterHalo), intent(inout) :: self

    !# <objectDestructor name="self%darkMatterProfile_"    />
    !# <objectDestructor name="self%freefallTimeAvailable_"/>
    return
  end subroutine darkMatterHaloDestructor

  double precision function darkMatterHaloRadiusGrowthRate(self,node)
    !% Returns the freefall radius growth rate (in Mpc/Gyr) in the hot atmosphere.
    implicit none
    class           (freefallRadiusDarkMatterHalo ), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    double precision                                               :: timeAvailable, timeAvailableIncreaseRate

    ! Get the time available for freefall.
    timeAvailable                 =+self%freefallTimeAvailable_%timeAvailable             (node                          )
    ! Get the rate of increase of the time available for freefall.
    timeAvailableIncreaseRate     =+self%freefallTimeAvailable_%timeAvailableIncreaseRate (node                          )
    ! Get freefall radius increase rate from dark matter profile.
    darkMatterHaloRadiusGrowthRate=+self%darkMatterProfile_    %freefallRadiusIncreaseRate(node,timeAvailable            ) &
         &                         *                                                            timeAvailableIncreaseRate
    return
  end function darkMatterHaloRadiusGrowthRate
  
  double precision function darkMatterHaloRadius(self,node)
    !% Return the freefall radius in the darkMatterHalo model.
    implicit none
    class           (freefallRadiusDarkMatterHalo ), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    double precision                                               :: timeAvailable

    ! Get the time available for freefall.
    timeAvailable       =self%freefallTimeAvailable_%timeAvailable (node              )
    ! Get freefall radius from dark matter profile.
    darkMatterHaloRadius=self%darkMatterProfile_    %freefallRadius(node,timeAvailable)
    return
  end function darkMatterHaloRadius
