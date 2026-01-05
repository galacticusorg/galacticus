!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

  !!{
  Implementation of a simple freefall radius class.
  !!}

  use :: Cooling_Freefall_Times_Available, only : freefallTimeAvailable, freefallTimeAvailableClass
  use :: Dark_Matter_Profiles_DMO        , only : darkMatterProfileDMO , darkMatterProfileDMOClass

  !![
  <freefallRadius name="freefallRadiusDarkMatterHalo">
   <description>
    A freefall radius class that assumes that the freefall radius corresponds to the radius at which the freefall time in the
    dark matter halo equals the time available for freefall (see \refPhysics{freefallTimeAvailable}).
   </description>
  </freefallRadius>
  !!]
  type, extends(freefallRadiusClass) :: freefallRadiusDarkMatterHalo
     !!{
     Implementation of freefall radius class in which the freefall radius is based on the freefall time in the dark matter halo.
     !!}
     private
     class(darkMatterProfileDMOClass ), pointer :: darkMatterProfileDMO_  => null()
     class(freefallTimeAvailableClass), pointer :: freefallTimeAvailable_ => null()
   contains
     final     ::                     darkMatterHaloDestructor
     procedure :: radius           => darkMatterHaloRadius
     procedure :: radiusGrowthRate => darkMatterHaloRadiusGrowthRate
  end type freefallRadiusDarkMatterHalo

  interface freefallRadiusDarkMatterHalo
     !!{
     Constructors for the darkMatterHalo freefall radius class.
     !!}
     module procedure darkMatterHaloConstructorParameters
     module procedure darkMatterHaloConstructorInternal
  end interface freefallRadiusDarkMatterHalo

contains

  function darkMatterHaloConstructorParameters(parameters) result(self)
    !!{
    Constructor for the darkMatterHalo freefall radius class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (freefallRadiusDarkMatterHalo)                :: self
    type (inputParameters             ), intent(inout) :: parameters
    class(darkMatterProfileDMOClass   ), pointer       :: darkMatterProfileDMO_
    class(freefallTimeAvailableClass  ), pointer       :: freefallTimeAvailable_

    !![
    <objectBuilder class="darkMatterProfileDMO"  name="darkMatterProfileDMO_"  source="parameters"/>
    <objectBuilder class="freefallTimeAvailable" name="freefallTimeAvailable_" source="parameters"/>
    !!]
    self=freefallRadiusDarkMatterHalo(darkMatterProfileDMO_,freefallTimeAvailable_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileDMO_" />
    <objectDestructor name="freefallTimeAvailable_"/>
    !!]
    return
  end function darkMatterHaloConstructorParameters

  function darkMatterHaloConstructorInternal(darkMatterProfileDMO_,freefallTimeAvailable_) result(self)
    !!{
    Internal constructor for the darkMatterHalo freefall radius class.
    !!}
    implicit none
    type (freefallRadiusDarkMatterHalo)                        :: self
    class(darkMatterProfileDMOClass   ), intent(in   ), target :: darkMatterProfileDMO_
    class(freefallTimeAvailableClass  ), intent(in   ), target :: freefallTimeAvailable_
    !![
    <constructorAssign variables="*darkMatterProfileDMO_, *freefallTimeAvailable_"/>
    !!]

    return
  end function darkMatterHaloConstructorInternal

  subroutine darkMatterHaloDestructor(self)
    !!{
    Destructor for the darkMatterHalo freefall radius class.
    !!}
    implicit none
    type(freefallRadiusDarkMatterHalo), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_" />
    <objectDestructor name="self%freefallTimeAvailable_"/>
    !!]
    return
  end subroutine darkMatterHaloDestructor

  double precision function darkMatterHaloRadiusGrowthRate(self,node) result(radiusGrowthRate)
    !!{
    Returns the freefall radius growth rate (in Mpc/Gyr) in the hot atmosphere.
    !!}
    use :: Mass_Distributions, only : massDistributionClass
    implicit none
    class           (freefallRadiusDarkMatterHalo ), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    class           (massDistributionClass        ), pointer       :: massDistribution_
    double precision                                               :: timeAvailable    , timeAvailableIncreaseRate

    ! Get the time available for freefall.
    timeAvailable                 =  +self             %freefallTimeAvailable_    %timeAvailable            (node                     )
    ! Get the rate of increase of the time available for freefall.
    timeAvailableIncreaseRate     =  +self             %freefallTimeAvailable_    %timeAvailableIncreaseRate(node                     )
    ! Get freefall radius increase rate from dark matter profile.
    massDistribution_             =>  self             %darkMatterProfileDMO_     %get                      (node                     )
    radiusGrowthRate              =  +massDistribution_%radiusFreefallIncreaseRate                          (timeAvailable            ) &
         &                           *                                                                       timeAvailableIncreaseRate
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end function darkMatterHaloRadiusGrowthRate

  double precision function darkMatterHaloRadius(self,node) result(radius)
    !!{
    Return the freefall radius in the darkMatterHalo model.
    !!}
    use :: Mass_Distributions, only : massDistributionClass
    implicit none
    class           (freefallRadiusDarkMatterHalo ), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    class           (massDistributionClass        ), pointer       :: massDistribution_
    double precision                                               :: timeAvailable

    ! Get the time available for freefall.
    timeAvailable     =  self             %freefallTimeAvailable_%timeAvailable(node         )
    ! Get freefall radius from dark matter profile.
    massDistribution_ => self             %darkMatterProfileDMO_ %get          (node         )
    radius            =  massDistribution_%radiusFreefall                      (timeAvailable)
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end function darkMatterHaloRadius
