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

!+    Contributions to this file made by:  Anthony Pullen, Andrew Benson.

  !!{
  Contains a class which implements a tidal heating rate model in which the heating rate is always zero.
  !!}

  !![
  <satelliteTidalHeatingRate name="satelliteTidalHeatingRateZero">
   <description>A satellite tidal heating rate class which implements a tidal heating rate model in which the heating rate is always zero.</description>
  </satelliteTidalHeatingRate>
  !!]
  type, extends(satelliteTidalHeatingRateClass) :: satelliteTidalHeatingRateZero
     !!{
     A satellite tidal heating rate class which implements a tidal heating rate model in which the heating rate is always zero.
     !!}
     private
   contains
     procedure :: heatingRate => zeroHeatingRate
  end type satelliteTidalHeatingRateZero

  interface satelliteTidalHeatingRateZero
     !!{
     Constructors for the zero satellite tidal heating rate class.
     !!}
     module procedure zeroConstructorParameters
  end interface satelliteTidalHeatingRateZero

contains

  function zeroConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{satelliteTidalHeatingRateZero} satellite tidal heating rate class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(satelliteTidalHeatingRateZero)                :: self
    type(inputParameters              ), intent(inout) :: parameters

    self=satelliteTidalHeatingRateZero()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function zeroConstructorParameters

  double precision function zeroHeatingRate(self,node)
    !!{
    Return the the tidal heating rate for the given node.
    !!}
    implicit none
    class(satelliteTidalHeatingRateZero), intent(inout) :: self
    type (treeNode                     ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    zeroHeatingRate=0.0d0
    return
  end function zeroHeatingRate
