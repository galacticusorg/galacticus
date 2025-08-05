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
  Implements a model of ram pressure stripping of hot halos which always returns zero force.
  !!}

  !![
  <hotHaloRamPressureForce name="hotHaloRamPressureForceZero">
   <description>
    A hot halo ram pressure force class which assumes a zero ram pressure force due to the hot halo.
   </description>
  </hotHaloRamPressureForce>
  !!]
  type, extends(hotHaloRamPressureForceClass) :: hotHaloRamPressureForceZero
     !!{
     Implementation of a hot halo ram pressure force class which always returns zero force.
     !!}
     private
   contains
     procedure :: force => zeroForce
  end type hotHaloRamPressureForceZero

  interface hotHaloRamPressureForceZero
     !!{
     Constructors for the \refClass{hotHaloRamPressureForceZero} hot halo ram pressure force class.
     !!}
     module procedure zeroConstructorParameters
  end interface hotHaloRamPressureForceZero

contains

  function zeroConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{hotHaloRamPressureForceZero} hot halo ram pressure force class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(hotHaloRamPressureForceZero)                :: self
    type(inputParameters            ), intent(inout) :: parameters

    self=hotHaloRamPressureForceZero()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function zeroConstructorParameters

  double precision function zeroForce(self,node)
    !!{
    Return a zero ram pressure force due to the hot halo.
    !!}
    implicit none
    class(hotHaloRamPressureForceZero), intent(inout) :: self
    type (treeNode                   ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    zeroForce=0.0d0
    return
  end function zeroForce
