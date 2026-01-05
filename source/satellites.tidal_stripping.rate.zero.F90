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

!+    Contributions to this file made by:  Anthony Pullen, Andrew Benson.

  !!{
  Implementation of a satellite tidal stripping class in which the stripping rate is always zero.
  !!}


  !![
  <satelliteTidalStripping name="satelliteTidalStrippingZero">
   <description>A satellite tidal stripping class in which the stripping rate is always zero.</description>
  </satelliteTidalStripping>
  !!]
  type, extends(satelliteTidalStrippingClass) :: satelliteTidalStrippingZero
     !!{
     Implementation of a satellite tidal stripping class in which the stripping rate is always zero.
     !!}
     private
   contains
     procedure :: massLossRate => zeroMassLossRate
  end type satelliteTidalStrippingZero

  interface satelliteTidalStrippingZero
     !!{
     Constructors for the \refClass{satelliteTidalStrippingZero} satellite tidal stripping class.
     !!}
     module procedure zeroConstructorParameters
  end interface satelliteTidalStrippingZero

contains

  function zeroConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{satelliteTidalStrippingZero} satellite tidal stripping class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(satelliteTidalStrippingZero)                :: self
    type(inputParameters            ), intent(inout) :: parameters

    self=satelliteTidalStrippingZero()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function zeroConstructorParameters

  double precision function zeroMassLossRate(self,node)
    !!{
    Return a mass loss rate for satellites due to tidal stripping which is always zero.
    !!}
    implicit none
    class(satelliteTidalStrippingZero), intent(inout) :: self
    type (treeNode                   ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    zeroMassLossRate=0.0d0
    return
  end function zeroMassLossRate

