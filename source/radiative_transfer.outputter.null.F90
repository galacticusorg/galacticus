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

  !![
  <radiativeTransferOutputter name="radiativeTransferOutputterNull">
   <description>A null radiative transfer outputter class.</description>
  </radiativeTransferOutputter>
  !!]
  type, extends(radiativeTransferOutputterClass) :: radiativeTransferOutputterNull
     !!{
     Implementation of a null radiative transfer outputter class.
     !!}
     private
   contains
  end type radiativeTransferOutputterNull

  interface radiativeTransferOutputterNull
     !!{
     Constructors for the \refClass{radiativeTransferOutputterNull} radiative transfer outputter packet class.
     !!}
     module procedure nullConstructorParameters
  end interface radiativeTransferOutputterNull
  
contains

  function nullConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{radiativeTransferOutputterNull} radiative transfer outputter class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(radiativeTransferOutputterNull)                :: self
    type(inputParameters               ), intent(inout) :: parameters
    
    self=radiativeTransferOutputterNull()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function nullConstructorParameters
