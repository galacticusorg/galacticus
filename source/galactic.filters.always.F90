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
Implements a galactic filter which always passes.
!!}

  !![
  <galacticFilter name="galacticFilterAlways">
   <description>A galactic filter which always passes. (Used mostly for testing purposes.)</description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterAlways
     !!{
     A galactic filter class which always passes.
     !!}
     private
   contains
     procedure :: passes => alwaysPasses
  end type galacticFilterAlways

  interface galacticFilterAlways
     !!{
     Constructors for the \refClass{galacticFilterAlways} galactic filter class.
     !!}
     module procedure alwaysConstructorParameters
  end interface galacticFilterAlways

contains

  function alwaysConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{galacticFilterAlways} galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(galacticFilterAlways)                :: self
    type(inputParameters     ), intent(inout) :: parameters

    self=galacticFilterAlways()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function alwaysConstructorParameters

  logical function alwaysPasses(self,node)
    !!{
    Implement an always-pass galactic filter.
    !!}
    implicit none
    class(galacticFilterAlways), intent(inout)         :: self
    type (treeNode            ), intent(inout), target :: node
    !$GLC attributes unused :: self, node

    alwaysPasses=.true.
    return
  end function alwaysPasses
