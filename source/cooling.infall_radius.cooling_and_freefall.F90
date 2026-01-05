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
  Implementation of an infall radius calculation in which the infall radius is the smaller of the cooling and freefall radii.
  !!}

  use :: Cooling_Radii , only : coolingRadius , coolingRadiusClass
  use :: Freefall_Radii, only : freefallRadius, freefallRadiusClass

  !![
  <coolingInfallRadius name="coolingInfallRadiusCoolingFreefall">
   <description>
    A cooling infall radius class calculation which assumes that the infall radius is equal to the smaller of the cooling and
    freefall radii (see \refPhysics{coolingRadius} and \refPhysics{freefallRadius}).
   </description>
  </coolingInfallRadius>
  !!]
  type, extends(coolingInfallRadiusClass) :: coolingInfallRadiusCoolingFreefall
     !!{
     Implementation of an infall radius calculation in which the infall radius is the smaller of the cooling and freefall radii.
     !!}
     private
     class(coolingRadiusClass ), pointer :: coolingRadius_ => null()
     class(freefallRadiusClass), pointer :: freefallRadius_ => null()
   contains
     final     ::                       coolingFreefallDestructor
     procedure :: radius             => coolingFreefallRadius
     procedure :: radiusIncreaseRate => coolingFreefallRadiusIncreaseRate
  end type coolingInfallRadiusCoolingFreefall

  interface coolingInfallRadiusCoolingFreefall
     !!{
     Constructors for the cooling radius infall radii class.
     !!}
     module procedure coolingFreefallConstructorParameters
     module procedure coolingFreefallConstructorInternal
  end interface coolingInfallRadiusCoolingFreefall

contains

  function coolingFreefallConstructorParameters(parameters) result(self)
    !!{
    Constructor for the cooling radius infall radii class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (coolingInfallRadiusCoolingFreefall)                :: self
    type (inputParameters                   ), intent(inout) :: parameters
    class(coolingRadiusClass                ), pointer       :: coolingRadius_
    class(freefallRadiusClass               ), pointer       :: freefallRadius_

    !![
    <objectBuilder class="coolingRadius"  name="coolingRadius_"  source="parameters"/>
    <objectBuilder class="freefallRadius" name="freefallRadius_" source="parameters"/>
    !!]
    self=coolingInfallRadiusCoolingFreefall(coolingRadius_,freefallRadius_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="coolingRadius_" />
    <objectDestructor name="freefallRadius_"/>
    !!]
    return
  end function coolingFreefallConstructorParameters

  function coolingFreefallConstructorInternal(coolingRadius_,freefallRadius_) result(self)
    !!{
    Internal constructor for the cooling radius infall radii class.
    !!}
    implicit none
    type (coolingInfallRadiusCoolingFreefall)                        :: self
    class(coolingRadiusClass                ), intent(in   ), target :: coolingRadius_
    class(freefallRadiusClass               ), intent(in   ), target :: freefallRadius_
    !![
    <constructorAssign variables="*coolingRadius_, *freefallRadius_"/>
    !!]

    return
  end function coolingFreefallConstructorInternal

  subroutine coolingFreefallDestructor(self)
    !!{
    Destructor for the cooling radius infall radii class.
    !!}
    implicit none
    type(coolingInfallRadiusCoolingFreefall), intent(inout) :: self

    !![
    <objectDestructor name="self%coolingRadius_"  />
    <objectDestructor name="self%freefallRadius_" />
    !!]
    return
  end subroutine coolingFreefallDestructor

  double precision function coolingFreefallRadius(self,node)
    !!{
    Return the infall radius in the ``cooling radius'' model in Mpc/Gyr.
    !!}
    implicit none
    class           (coolingInfallRadiusCoolingFreefall), intent(inout) :: self
    type            (treeNode                          ), intent(inout) :: node
    double precision                                                    :: radiusCooling, radiusFreefall

    radiusCooling =self%coolingRadius_ %radius(node)
    radiusFreefall=self%freefallRadius_%radius(node)
    coolingFreefallRadius=min(radiusCooling,radiusFreefall)
    return
  end function coolingFreefallRadius

  double precision function coolingFreefallRadiusIncreaseRate(self,node)
    !!{
    Return the growth rate of the infall radius in the ``cooling radius'' model in Mpc/Gyr.
    !!}
    implicit none
    class           (coolingInfallRadiusCoolingFreefall), intent(inout) :: self
    type            (treeNode                          ), intent(inout) :: node
    double precision                                                    :: radiusCooling, radiusFreefall

    radiusCooling =self%coolingRadius_ %radius(node)
    radiusFreefall=self%freefallRadius_%radius(node)
    if (radiusCooling < radiusFreefall) then
       coolingFreefallRadiusIncreaseRate=self%coolingRadius_ %radiusGrowthRate(node)
    else
       coolingFreefallRadiusIncreaseRate=self%freefallRadius_%radiusGrowthRate(node)
    end if
    return
  end function coolingFreefallRadiusIncreaseRate
