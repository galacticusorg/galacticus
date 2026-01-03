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
Implements a merger tree filter class which is the ``all'' combination of a set of other filters.
!!}

  !![
  <mergerTreeFilter name="mergerTreeFilterAll">
   <description>A merger tree filter class which is the {\normalfont \ttfamily all} combination of a set of other filters.</description>
   <linkedList type="filterList" variable="filters" next="next" object="filter_" objectType="mergerTreeFilterClass"/>
  </mergerTreeFilter>
  !!]
  type, extends(mergerTreeFilterClass) :: mergerTreeFilterAll
     !!{
     A merger tree filter class which is the {\normalfont \ttfamily all} combination of a set of other filters.
     !!}
     private
     type(filterList), pointer :: filters => null()
  contains
     final     ::           allDestructor
     procedure :: passes => allPasses
  end type mergerTreeFilterAll

  interface mergerTreeFilterAll
     !!{
     Constructors for the all merger tree filter class.
     !!}
     module procedure allConstructorParameters
     module procedure allConstructorInternal
  end interface mergerTreeFilterAll

contains

  function allConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeFilterAll} merger tree filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (mergerTreeFilterAll)                :: self
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
  end function allConstructorParameters

  function allConstructorInternal(filters) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeFilterAll} filter class.
    !!}
    implicit none
    type(mergerTreeFilterAll)                        :: self
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
  end function allConstructorInternal

  subroutine allDestructor(self)
    !!{
    Destructor for the all merger tree filter class.
    !!}
    implicit none
    type(mergerTreeFilterAll), intent(inout) :: self
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
  end subroutine allDestructor

  logical function allPasses(self,tree)
    !!{
    Apply a set of filters to a {\normalfont \ttfamily tree} combined with ``all'' operations.
    !!}
    implicit none
    class(mergerTreeFilterAll), intent(inout) :: self
    type (mergerTree         ), intent(in   ) :: tree
    type (filterList         ), pointer       :: filter_

    ! Assume the tree passes initially. Iterate through filters and evaluate each one. If any one evaluates to false, exit the
    ! iteration.
    allPasses =  .true.
    filter_   => self%filters
    do while (associated(filter_))
       allPasses=filter_%filter_%passes(tree)
       if (allPasses) then
          filter_ => filter_%next
       else
          filter_ => null()
       end if
    end do
    return
  end function allPasses
