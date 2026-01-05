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
  Implements a \gls{cgm} heating class with zero heating.
  !!}

  !![
  <circumgalacticMediumHeating name="circumgalacticMediumHeatingZero">
   <description>
    A \gls{cgm} heating class with zero heating.
   </description>
  </circumgalacticMediumHeating>
  !!]
  type, extends(circumgalacticMediumHeatingClass) :: circumgalacticMediumHeatingZero
     !!{
     A \gls{cgm} heating class with zero heating.
     !!}
     private
   contains
     procedure :: heatingRate => zeroHeatingRate
  end type circumgalacticMediumHeatingZero
  
  interface circumgalacticMediumHeatingZero
     !!{
     Constructors for the \refClass{circumgalacticMediumHeatingZero} class.
     !!}
     module procedure zeroConstructorParameters
  end interface circumgalacticMediumHeatingZero

contains

  function zeroConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{circumgalacticMediumHeatingZero} class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(circumgalacticMediumHeatingZero)                :: self
    type(inputParameters                ), intent(inout) :: parameters
    
    self=circumgalacticMediumHeatingZero()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function zeroConstructorParameters

  double precision function zeroHeatingRate(self,node) result(rateHeating)
    !!{
    Compute the heating rate of the \gls{cgm}, assumed to be always zero.
    !!}
    implicit none
    class(circumgalacticMediumHeatingZero), intent(inout) :: self
    type (treeNode                       ), intent(inout) :: node

    rateHeating=+0.0d0
    return
  end function zeroHeatingRate
