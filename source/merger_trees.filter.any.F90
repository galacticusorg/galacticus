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
Implements a merger tree filter class which is the ``any'' combination of a set of other filters.
!!}

  !![
  <mergerTreeFilter name="mergerTreeFilterAny">
   <description>A merger tree filter class which is the {\normalfont \ttfamily any} combination of a set of other filters.</description>
   <linkedList type="filterList" variable="filters" next="next" object="filter_" objectType="mergerTreeFilterClass"/>
  </mergerTreeFilter>
  !!]
  type, extends(mergerTreeFilterClass) :: mergerTreeFilterAny
     !!{
     A merger tree filter class which is the {\normalfont \ttfamily any} combination of a set of other filters.
     !!}
     private
     type(filterList), pointer :: filters => null()
  contains
     final     ::           anyDestructor
     procedure :: passes => anyPasses
  end type mergerTreeFilterAny

  interface mergerTreeFilterAny
     !!{
     Constructors for the any merger tree filter class.
     !!}
     module procedure anyConstructorParameters
     module procedure anyConstructorInternal
  end interface mergerTreeFilterAny

contains

  function anyConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeFilterAny} merger tree filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (mergerTreeFilterAny)                :: self
    type   (inputParameters    ), intent(inout) :: parameters
    type   (filterList         ), pointer       :: filter_
    integer                                     :: i

    self   %filters => null()
    filter_         => null()
    do i=1,parameters%copiesCount('mergerTreeFilter',zeroIfNotPresent=.true.)
       if (associated(filter_)) then
          allocate(filter_%next)
          filter_ => filter_%next
       else
          allocate(self%filters)
          filter_ => self%filters
       end if
       !![
       <objectBuilder class="mergerTreeFilter" name="filter_%filter_" source="parameters" copy="i" />
       !!]
    end do
    !![
    <inputParametersValidate source="parameters" multiParameters="mergerTreeFilter"/>
    !!]
    return
  end function anyConstructorParameters

  function anyConstructorInternal(filters) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeFilterAny} filter class.
    !!}
    implicit none
    type(mergerTreeFilterAny)                        :: self
    type(filterList         ), target, intent(in   ) :: filters
    type(filterList         ), pointer               :: filter_

    self   %filters => filters
    filter_         => filters
    do while (associated(filter_))
       !![
       <referenceCountIncrement owner="filter_" object="filter_"/>
       !!]
       filter_ => filter_%next
    end do
    return
  end function anyConstructorInternal

  subroutine anyDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeFilterAny} merger tree filter class.
    !!}
    implicit none
    type(mergerTreeFilterAny), intent(inout) :: self
    type(filterList         ), pointer       :: filter_, filterNext

    if (associated(self%filters)) then
       filter_ => self%filters
       do while (associated(filter_))
          filterNext => filter_%next
          !![
          <objectDestructor name="filter_%filter_"/>
          !!]
          deallocate(filter_)
          filter_ => filterNext
       end do
    end if
    return
  end subroutine anyDestructor

  logical function anyPasses(self,tree)
    !!{
    Apply a set of filters to a {\normalfont \ttfamily tree} combined with {\normalfont \ttfamily any} operations.
    !!}
    implicit none
    class(mergerTreeFilterAny), intent(inout) :: self
    type (mergerTree         ), intent(in   ) :: tree
    type (filterList         ), pointer       :: filter_

    ! Assume the tree fails to pass initially. Iterate through filters and evaluate each one. If any one evaluates to true, exit
    ! the iteration.
    anyPasses =  .false.
    filter_  => self%filters
    do while (associated(filter_))
       anyPasses=filter_%filter_%passes(tree)
       if (anyPasses) then
          filter_ => null()
       else
          filter_ => filter_%next
       end if
    end do
    return
  end function anyPasses
