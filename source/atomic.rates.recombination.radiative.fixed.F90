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
  An implementation of atomic radiative recombination rates which uses a fixed rate coefficient.
  !!}

  !![
  <atomicRecombinationRateRadiative name="atomicRecombinationRateRadiativeFixed">
   <description>Atomic radiative recombination rates with a fixed rate coefficient.</description>
  </atomicRecombinationRateRadiative>
  !!]
  type, extends(atomicRecombinationRateRadiativeClass) :: atomicRecombinationRateRadiativeFixed
     !!{
     A radiative recombination rate class which uses a fixed rate coefficient.
     !!}
     private
     double precision :: rateCoefficient
   contains
     procedure :: rate => fixedRate
  end type atomicRecombinationRateRadiativeFixed

  interface atomicRecombinationRateRadiativeFixed
     !!{
     Constructors for the \refClass{atomicRecombinationRateRadiativeFixed} atomic radiative recombination class.
     !!}
     module procedure fixedConstructorParameters
     module procedure fixedConstructorInternal
  end interface atomicRecombinationRateRadiativeFixed

contains

  function fixedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{atomicRecombinationRateRadiativeFixed} atomic radiative recombination class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (atomicRecombinationRateRadiativeFixed)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    double precision                                                       :: rateCoefficient

    !![
    <inputParameter>
      <name>rateCoefficient</name>
      <description>The rate coefficient (in units of cm$^3$ s$^{-1}$) for radiative recombination.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=atomicRecombinationRateRadiativeFixed(rateCoefficient)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function fixedConstructorParameters

  function fixedConstructorInternal(rateCoefficient) result(self)
    !!{
    Internal constructor for the \refClass{atomicRecombinationRateRadiativeFixed} atomic radiative recombination class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (atomicRecombinationRateRadiativeFixed)                :: self
    double precision                                       , intent(in   ) :: rateCoefficient
    !![
    <constructorAssign variables="rateCoefficient"/>
    !!]
    
    return
  end function fixedConstructorInternal

  double precision function fixedRate(self,atomicNumber,ionizationState,temperature,level)
    !!{
    Returns a fixed rate coefficient.
    !!}
    implicit none
    class           (atomicRecombinationRateRadiativeFixed), intent(inout)           :: self
    integer                                                , intent(in   )           :: atomicNumber, ionizationState
    double precision                                       , intent(in   )           :: temperature
    type            (enumerationRecombinationCaseType     ), intent(in   ), optional :: level
    !$GLC attributes unused :: atomicNumber, ionizationState, temperature, level

    fixedRate=self%rateCoefficient
    return
  end function fixedRate
