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
  Implements a stellar population selector class which returns a fixed population.
  !!}

  !![
  <stellarPopulationSelector name="stellarPopulationSelectorFixed">
   <description>
    A stellar population selector class which selects a fixed population irrespective of physical conditions.
   </description>
  </stellarPopulationSelector>
  !!]
  type, extends(stellarPopulationSelectorClass) :: stellarPopulationSelectorFixed
     !!{
     A fixed stellar population selector class.
     !!}
     private
     class(stellarPopulationClass), pointer :: stellarPopulation_ => null()
   contains
     final     ::                                 fixedDestructor
     procedure :: select                       => fixedSelect
     procedure :: isStarFormationRateDependent => fixedIsStarFormationRateDependent
  end type stellarPopulationSelectorFixed

  interface stellarPopulationSelectorFixed
     !!{
     Constructors for the \refClass{stellarPopulationSelectorFixed} stellar population selector class.
     !!}
     module procedure fixedConstructorParameters
     module procedure fixedConstructorInternal
  end interface stellarPopulationSelectorFixed

contains

  function fixedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{stellarPopulationSelectorFixed} stellar population class which takes a parameter list as input.
    !!}
    use :: Input_Parameters   , only : inputParameter   , inputParameters
    use :: Stellar_Populations, only : stellarPopulation, stellarPopulationClass
    implicit none
    type (stellarPopulationSelectorFixed)                :: self
    type (inputParameters               ), intent(inout) :: parameters
    class(stellarPopulationClass        ), pointer       :: stellarPopulation_

    !![
    <objectBuilder class="stellarPopulation" name="stellarPopulation_" source="parameters"/>
    !!]
    self=stellarPopulationSelectorFixed(stellarPopulation_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="stellarPopulation_"/>
    !!]
    return
  end function fixedConstructorParameters

  function fixedConstructorInternal(stellarPopulation_) result(self)
    !!{
    Internal constructor for the \refClass{stellarPopulationSelectorFixed} stellar population selector class.
    !!}
    implicit none
    type (stellarPopulationSelectorFixed)                        :: self
    class(stellarPopulationClass        ), intent(in   ), target :: stellarPopulation_
    !![
    <constructorAssign variables="*stellarPopulation_"/>
    !!]

    return
  end function fixedConstructorInternal

  subroutine fixedDestructor(self)
    !!{
    Destructor for the \refClass{stellarPopulationSelectorFixed} stellar population selector class.
    !!}
    implicit none
    type(stellarPopulationSelectorFixed), intent(inout) :: self

    !![
    <objectDestructor name="self%stellarPopulation_"/>
    !!]
    return
  end subroutine fixedDestructor

  function fixedSelect(self,rateStarFormation,abundances_,component)
    !!{
    Return a fixed stellar population.
    !!}
    implicit none
    class           (stellarPopulationClass        ), pointer       :: fixedSelect
    class           (stellarPopulationSelectorFixed), intent(inout) :: self
    double precision                                , intent(in   ) :: rateStarFormation
    type            (abundances                    ), intent(in   ) :: abundances_
    class           (nodeComponent                 ), intent(in   ) :: component
    !$GLC attributes unused :: rateStarFormation, abundances_, component

    fixedSelect => self%stellarPopulation_
    return
  end function fixedSelect

  logical function fixedIsStarFormationRateDependent(self)
    !!{
    Return false indicating that stellar population selection is not dependent on star formation rate.
    !!}
    implicit none
    class(stellarPopulationSelectorFixed), intent(inout) :: self
    !$GLC attributes unused :: self

    fixedIsStarFormationRateDependent=.false.
    return
  end function fixedIsStarFormationRateDependent

