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
  Implementation of a zero cooling rate class.
  !!}

  !![
  <coolingRate name="coolingRateZero">
   <description>A cooling rate class in which the cooling rate is always zero.</description>
  </coolingRate>
  !!]
  type, extends(coolingRateClass) :: coolingRateZero
     !!{
     Implementation of cooling rate class in which the cooling rate is always zero.
     !!}
     private
   contains
     procedure :: rate => zeroRate
  end type coolingRateZero

  interface coolingRateZero
     !!{
     Constructors for the zero cooling rate class.
     !!}
     module procedure zeroConstructorParameters
  end interface coolingRateZero

contains

  function zeroConstructorParameters(parameters) result(self)
    !!{
    Constructor for the zero cooling rate class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(coolingRateZero)                :: self
    type(inputParameters), intent(inout) :: parameters

    self=coolingRateZero()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function zeroConstructorParameters

  double precision function zeroRate(self,node)
    !!{
    Returns the cooling rate (in $M_\odot$ Gyr$^{-1}$) in the hot atmosphere for a model in which this rate is always zero.
    !!}
    implicit none
    class           (coolingRateZero), intent(inout) :: self
    type            (treeNode       ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    zeroRate=0.0d0
    return
  end function zeroRate
