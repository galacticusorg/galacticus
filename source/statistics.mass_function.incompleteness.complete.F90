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
  Implements calculations of incompleteness assuming a complete sample.
  !!}

  !![
  <massFunctionIncompleteness name="massFunctionIncompletenessComplete">
   <description>
  A mass function incompleteness class which assumes a fully complete mass function.
   </description>
  </massFunctionIncompleteness>
  !!]
  type, extends(massFunctionIncompletenessClass) :: massFunctionIncompletenessComplete
     !!{
     A class implementing incompleteness calculations for a complete survey.
     !!}
     private
   contains
     procedure :: completeness => completeCompleteness
  end type massFunctionIncompletenessComplete

  interface massFunctionIncompletenessComplete
     !!{
     Constructors for the \refClass{massFunctionIncompletenessComplete} incompleteness class.
     !!}
     module procedure completeConstructorParameters
  end interface massFunctionIncompletenessComplete

contains

  function completeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{massFunctionIncompletenessComplete} incompleteness class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(massFunctionIncompletenessComplete)                :: self
    type(inputParameters                   ), intent(inout) :: parameters

    self=massFunctionIncompletenessComplete()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function completeConstructorParameters

  double precision function completeCompleteness(self,mass)
    !!{
    Return the completeness.
    !!}
    implicit none
    class           (massFunctionIncompletenessComplete), intent(inout) :: self
    double precision                                    , intent(in   ) :: mass
    !$GLC attributes unused :: self, mass

    completeCompleteness=1.0d0
    return
  end function completeCompleteness
