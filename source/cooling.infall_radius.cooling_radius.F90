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
  Implementation of a simple infall radius calculation, simply assuming that the infall radius equals the cooling radius.
  !!}

  use :: Cooling_Radii, only : coolingRadius, coolingRadiusClass

  !![
  <coolingInfallRadius name="coolingInfallRadiusCoolingRadius">
   <description>
    A cooling infall radius class that assumes that the infall radius equals the cooling radius (see
    \refPhysics{coolingRadius}).
   </description>
  </coolingInfallRadius>
  !!]
  type, extends(coolingInfallRadiusClass) :: coolingInfallRadiusCoolingRadius
     !!{
     Implementation of a simple infall radius calculation, simply assuming that the infall radius equals the cooling radius.
     !!}
     private
     class(coolingRadiusClass), pointer :: coolingRadius_ => null()
   contains
     final     ::                       coolingRadiusDestructor
     procedure :: radius             => coolingRadiusRadius
     procedure :: radiusIncreaseRate => coolingRadiusRadiusIncreaseRate
  end type coolingInfallRadiusCoolingRadius

  interface coolingInfallRadiusCoolingRadius
     !!{
     Constructors for the cooling radius infall radii class.
     !!}
     module procedure coolingRadiusConstructorParameters
     module procedure coolingRadiusConstructorInternal
  end interface coolingInfallRadiusCoolingRadius

contains

  function coolingRadiusConstructorParameters(parameters) result(self)
    !!{
    Constructor for the cooling radius infall radii class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (coolingInfallRadiusCoolingRadius)                :: self
    type (inputParameters                 ), intent(inout) :: parameters
    class(coolingRadiusClass              ), pointer       :: coolingRadius_

    !![
    <objectBuilder class="coolingRadius" name="coolingRadius_" source="parameters"/>
    !!]
    self=coolingInfallRadiusCoolingRadius(coolingRadius_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="coolingRadius_"/>
    !!]
    return
  end function coolingRadiusConstructorParameters

  function coolingRadiusConstructorInternal(coolingRadius_) result(self)
    !!{
    Internal constructor for the cooling radius infall radii class.
    !!}
    implicit none
    type (coolingInfallRadiusCoolingRadius)                        :: self
    class(coolingRadiusClass              ), intent(in   ), target :: coolingRadius_
    !![
    <constructorAssign variables="*coolingRadius_"/>
    !!]

    return
  end function coolingRadiusConstructorInternal

  subroutine coolingRadiusDestructor(self)
    !!{
    Destructor for the cooling radius infall radii class.
    !!}
    implicit none
    type(coolingInfallRadiusCoolingRadius), intent(inout) :: self

    !![
    <objectDestructor name="self%coolingRadius_" />
    !!]
    return
  end subroutine coolingRadiusDestructor

  double precision function coolingRadiusRadius(self,node)
    !!{
    Return the infall radius in the ``cooling radius'' model in Mpc/Gyr.
    !!}
    implicit none
    class(coolingInfallRadiusCoolingRadius), intent(inout) :: self
    type (treeNode                        ), intent(inout) :: node

    coolingRadiusRadius=self%coolingRadius_%radius(node)
    return
  end function coolingRadiusRadius

  double precision function coolingRadiusRadiusIncreaseRate(self,node)
    !!{
    Return the growth rate of the infall radius in the ``cooling radius'' model in Mpc/Gyr.
    !!}
    implicit none
    class(coolingInfallRadiusCoolingRadius), intent(inout) :: self
    type (treeNode                        ), intent(inout) :: node

    coolingRadiusRadiusIncreaseRate=self%coolingRadius_%radiusGrowthRate(node)
    return
  end function coolingRadiusRadiusIncreaseRate
