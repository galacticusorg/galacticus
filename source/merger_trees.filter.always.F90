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
Implements a merger tree filter which always passes.
!!}

  !![
  <mergerTreeFilter name="mergerTreeFilterAlways">
   <description>A merger tree filter which always passes. (Used mostly for testing purposes.)</description>
  </mergerTreeFilter>
  !!]
  type, extends(mergerTreeFilterClass) :: mergerTreeFilterAlways
     !!{
     A merger tree filter class which always passes.
     !!}
     private
   contains
     procedure :: passes => alwaysPasses
  end type mergerTreeFilterAlways

  interface mergerTreeFilterAlways
     !!{
     Constructors for the ``always'' merger tree filter class.
     !!}
     module procedure alwaysConstructorParameters
  end interface mergerTreeFilterAlways

contains

  function alwaysConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``always'' merger tree filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(mergerTreeFilterAlways)                :: self
    type(inputParameters       ), intent(inout) :: parameters

    self=mergerTreeFilterAlways()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function alwaysConstructorParameters

  logical function alwaysPasses(self,tree)
    !!{
    Implement an always-pass merger tree filter.
    !!}
    implicit none
    class(mergerTreeFilterAlways), intent(inout) :: self
    type (mergerTree            ), intent(in   ) :: tree
    !$GLC attributes unused :: self, tree

    alwaysPasses=.true.
    return
  end function alwaysPasses
