!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  An implementation of calculations of chemical reaction rates which assumes zero rates.
  !!}

  !![
  <chemicalReactionRate name="chemicalReactionRateZero">
   <description>
    A chemical reaction rate class in which all rates are zero.
   </description>
  </chemicalReactionRate>
  !!]
  type, extends(chemicalReactionRateClass) :: chemicalReactionRateZero
     !!{
     A chemical reaction rate class in which all rates are zero.
     !!}
     private
   contains
     procedure :: rates => zeroRates
  end type chemicalReactionRateZero

  interface chemicalReactionRateZero
     !!{
     Constructors for the \refClass{chemicalReactionRateZero} chemical reaction rates class.
     !!}
     module procedure zeroConstructorParameters
  end interface chemicalReactionRateZero

contains

  function zeroConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{chemicalReactionRateZero} chemical reaction rates class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(chemicalReactionRateZero)                :: self
    type(inputParameters         ), intent(inout) :: parameters

    self=chemicalReactionRateZero()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function zeroConstructorParameters

  subroutine zeroRates(self,lengthColumn,temperature,chemicalDensity,factorClumping,radiation,chemicalRates,node)
    !!{
    Return zero rates of chemical reactions.
    !!}
    implicit none
    class           (chemicalReactionRateZero), intent(inout), target :: self
    type            (chemicalAbundances      ), intent(in   )         :: chemicalDensity
    double precision                          , intent(in   )         :: lengthColumn   , temperature, &
         &                                                               factorClumping
    class           (radiationFieldClass     ), intent(inout)         :: radiation
    type            (chemicalAbundances      ), intent(inout)         :: chemicalRates
    type            (treeNode                ), intent(inout)         :: node
    !$GLC attributes unused :: self, chemicalDensity, lengthColumn, temperature, factorClumping, radiation, node

    call chemicalRates%reset()
    return
  end subroutine zeroRates
