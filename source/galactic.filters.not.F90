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
Implements an inverting filter.
!!}

  !![
  <galacticFilter name="galacticFilterNot">
   <description>A filter which simply inverts the result of another filter.</description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterNot
     !!{
     A galactic filter which simply inverts the result of another filter.
     !!}
     private
     class(galacticFilterClass), pointer :: galacticFilter_ => null()
   contains
     final     ::           notDestructor
     procedure :: passes => notPasses
  end type galacticFilterNot

  interface galacticFilterNot
     !!{
     Constructors for the \refClass{galacticFilterNot} galactic filter class.
     !!}
     module procedure notConstructorParameters
     module procedure notConstructorInternal
  end interface galacticFilterNot

contains

  function notConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{galacticFilterNot} galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (galacticFilterNot  )                :: self
    type (inputParameters    ), intent(inout) :: parameters
    class(galacticFilterClass), pointer       :: galacticFilter_

    !![
    <objectBuilder class="galacticFilter" name="galacticFilter_" source="parameters"/>
    !!]
    self=galacticFilterNot(galacticFilter_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticFilter_"/>
    !!]
    return
  end function notConstructorParameters

  function notConstructorInternal(galacticFilter_) result(self)
    !!{
    Internal constructor for the \refClass{galacticFilterNot} galactic filter class.
    !!}
    implicit none
    type (galacticFilterNot  )                        :: self
    class(galacticFilterClass), intent(in   ), target :: galacticFilter_
    !![
    <constructorAssign variables="*galacticFilter_"/>
    !!]

    return
  end function notConstructorInternal

  subroutine notDestructor(self)
    !!{
    Destructor for the \refClass{galacticFilterNot} galactic filter class.
    !!}
    implicit none
    type(galacticFilterNot), intent(inout) :: self

    !![
    <objectDestructor name="self%galacticFilter_"/>
    !!]
    return
  end subroutine notDestructor

  logical function notPasses(self,node)
    !!{
    Implement a not galactic filter.
    !!}
    implicit none
    class(galacticFilterNot), intent(inout)         :: self
    type (treeNode         ), intent(inout), target :: node

    notPasses=.not.self%galacticFilter_%passes(node)
    return
  end function notPasses
