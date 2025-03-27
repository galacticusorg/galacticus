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
Implements a null filter.
!!}

  !![
  <galacticFilter name="galacticFilterNull">
   <description>A filter which simply returns the result of another filter. This is intended for use in filter pipelines where it may be useful to optionally switch in this filter or a \refClass{galacticFilterNull} filter (for example).</description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterNull
     !!{
     A galactic filter which simply returns the result of another filter.
     !!}
     private
     class(galacticFilterClass), pointer :: galacticFilter_ => null()
   contains
     final     ::           nullDestructor
     procedure :: passes => nullPasses
  end type galacticFilterNull

  interface galacticFilterNull
     !!{
     Constructors for the {\normalfont \ttfamily null} galactic filter class.
     !!}
     module procedure nullConstructorParameters
     module procedure nullConstructorInternal
  end interface galacticFilterNull

contains

  function nullConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily null} galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (galacticFilterNull )                :: self
    type (inputParameters    ), intent(inout) :: parameters
    class(galacticFilterClass), pointer       :: galacticFilter_

    !![
    <objectBuilder class="galacticFilter" name="galacticFilter_" source="parameters"/>
    !!]
    self=galacticFilterNull(galacticFilter_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticFilter_"/>
    !!]
    return
  end function nullConstructorParameters

  function nullConstructorInternal(galacticFilter_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily null} galactic filter class.
    !!}
    implicit none
    type (galacticFilterNull )                        :: self
    class(galacticFilterClass), intent(in   ), target :: galacticFilter_
    !![
    <constructorAssign variables="*galacticFilter_"/>
    !!]

    return
  end function nullConstructorInternal

  subroutine nullDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily null} galactic filter class.
    !!}
    implicit none
    type(galacticFilterNull), intent(inout) :: self

    !![
    <objectDestructor name="self%galacticFilter_"/>
    !!]
    return
  end subroutine nullDestructor

  logical function nullPasses(self,node)
    !!{
    Implement a null galactic filter.
    !!}
    implicit none
    class(galacticFilterNull), intent(inout)         :: self
    type (treeNode          ), intent(inout), target :: node

    nullPasses=self%galacticFilter_%passes(node)
    return
  end function nullPasses
