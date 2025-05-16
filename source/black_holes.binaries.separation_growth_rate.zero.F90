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

!+    Contributions to this file made by:  Stéphane Mangeon, Andrew Benson.

  !!{
  Implements a black hole binary separation growth class in which the separation does not grow.
  !!}

  !![
  <blackHoleBinarySeparationGrowthRate name="blackHoleBinarySeparationGrowthRateZero">
   <description>
    A black hole binary separation growth class in which the separation does not grow.
   </description>
  </blackHoleBinarySeparationGrowthRate>
  !!]
  type, extends(blackHoleBinarySeparationGrowthRateClass) :: blackHoleBinarySeparationGrowthRateZero
     !!{
     A black hole binary separation growth class in which the separation does not grow.
     !!}
     private
   contains
     procedure :: growthRate => zeroGrowthRate
  end type blackHoleBinarySeparationGrowthRateZero

  interface blackHoleBinarySeparationGrowthRateZero
     !!{
     Constructors for the \refClass{blackHoleBinarySeparationGrowthRateZero} black hole binary separation growth rate class.
     !!}
     module procedure zeroConstructorParameters
  end interface blackHoleBinarySeparationGrowthRateZero

contains

  function zeroConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{blackHoleBinarySeparationGrowthRateZero} black hole binary separation growth rate class which takes a parameter
    set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(blackHoleBinarySeparationGrowthRateZero)                :: self
    type(inputParameters                        ), intent(inout) :: parameters

    self=blackHoleBinarySeparationGrowthRateZero()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function zeroConstructorParameters

  double precision function zeroGrowthRate(self,blackHole)
    !!{
    Returns a separation growth rate for a binary black hole that is always zero.
    !!}
    implicit none
    class(blackHoleBinarySeparationGrowthRateZero), intent(inout) :: self
    class(nodeComponentBlackHole                 ), intent(inout) :: blackHole
    !$GLC attributes unused :: self, blackHole

    zeroGrowthRate=0.0d0
    return
  end function zeroGrowthRate
