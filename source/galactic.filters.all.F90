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
Implements a galactic filter class which is the ``all'' combination of a set of other filters.
!!}

  !![
  <galacticFilter name="galacticFilterAll">
   <description>A galactic filter class which is the ``all'' combination of a set of other filters.</description>
   <linkedList type="filterList" variable="filters" next="next" object="filter_" objectType="galacticFilterClass"/>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterAll
     !!{
     A galactic filter class which is the ``all'' combination of a set of other filters.
     !!}
     private
     type(filterList), pointer :: filters => null()
  contains
     final     ::           allDestructor
     procedure :: passes => allPasses
  end type galacticFilterAll

  interface galacticFilterAll
     !!{
     Constructors for the all galactic filter class.
     !!}
     module procedure allConstructorParameters
     module procedure allConstructorInternal
  end interface galacticFilterAll

contains

  function allConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``all'' galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (galacticFilterAll)                :: self
    type   (inputParameters  ), intent(inout) :: parameters
    type   (filterList       ), pointer       :: filter_
    integer                                   :: i

    self   %filters => null()
    filter_         => null()
    do i=1,parameters%copiesCount('galacticFilter',zeroIfNotPresent=.true.)
       if (associated(filter_)) then
          allocate(filter_%next)
          filter_ => filter_%next
       else
          allocate(self%filters)
          filter_ => self%filters
       end if
       !![
       <objectBuilder class="galacticFilter" name="filter_%filter_" source="parameters" copy="i" />
       !!]
    end do
    !![
    <inputParametersValidate source="parameters" multiParameters="galacticFilter"/>
    !!]
    return
  end function allConstructorParameters

  function allConstructorInternal(filters) result(self)
    !!{
    Internal constructor for the ``all'' filter class.
    !!}
    implicit none
    type(galacticFilterAll)                        :: self
    type(filterList       ), target, intent(in   ) :: filters
    type(filterList       ), pointer               :: filter_

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
    Destructor for the all galactic filter class.
    !!}
    implicit none
    type(galacticFilterAll), intent(inout) :: self
    type(filterList       ), pointer       :: filter_, filterNext

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

  logical function allPasses(self,node)
    !!{
    Apply a set of filters to a {\normalfont \ttfamily node} combined with ``all'' operations.
    !!}
    implicit none
    class(galacticFilterAll), intent(inout)         :: self
    type (treeNode         ), intent(inout), target :: node
    type (filterList       ), pointer               :: filter_

    ! Assume the node passes initially. Iterate through filters and evaluate each one. If any one evaluates to false, exit the
    ! iteration.
    allPasses =  .true.
    filter_   => self%filters
    do while (associated(filter_))
       allPasses=filter_%filter_%passes(node)
       if (allPasses) then
          filter_ => filter_%next
       else
          filter_ => null()
       end if
    end do
    return
  end function allPasses
