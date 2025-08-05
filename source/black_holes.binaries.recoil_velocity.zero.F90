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
  Implements a black hole binary recoil velocity class in which the recoil velocity is zero.
  !!}

  !![
  <blackHoleBinaryRecoil name="blackHoleBinaryRecoilZero">
   <description>
    A black hole binary recoil class in which the recoil velocity is always zero.
   </description>
  </blackHoleBinaryRecoil>
  !!]
  type, extends(blackHoleBinaryRecoilClass) :: blackHoleBinaryRecoilZero
     !!{
     A black hole binary recoil class in which the recoil velocity is always zero.
     !!}
     private
   contains
     procedure :: velocity => zeroVelocity
  end type blackHoleBinaryRecoilZero

  interface blackHoleBinaryRecoilZero
     !!{
     Constructors for the \refClass{blackHoleBinaryRecoilZero} black hole binary recoil class.
     !!}
     module procedure zeroConstructorParameters
  end interface blackHoleBinaryRecoilZero

contains

  function zeroConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{blackHoleBinaryRecoilZero} black hole binary recoil class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(blackHoleBinaryRecoilZero)                :: self
    type(inputParameters          ), intent(inout) :: parameters

    self=blackHoleBinaryRecoilZero()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function zeroConstructorParameters

  double precision function zeroVelocity(self,blackHole1,blackHole2)
    !!{
    Compute the recoil velocity for a pair of merging black holes.
    !!}
    implicit none
    class(blackHoleBinaryRecoilZero), intent(inout) :: self
    class(nodeComponentBlackHole   ), intent(inout) :: blackHole1, blackHole2
    !$GLC attributes unused :: self, blackHole1, blackHole2

    zeroVelocity=0.0d0
    return
  end function zeroVelocity
