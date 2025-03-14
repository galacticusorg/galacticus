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
Implements a galactic filter class which is the ``any'' combination of a set of other filters.
!!}

  !![
  <galacticFilter name="galacticFilterAny">
   <description>A galactic filter class which is the ``any'' combination of a set of other filters.</description>
   <linkedList type="filterList" variable="filters" next="next" object="filter_" objectType="galacticFilterClass"/>
  </galacticFilter>
  !!]

  type, extends(galacticFilterClass) :: galacticFilterAny
     !!{
     A galactic filter class which is the ``any'' combination of a set of other filters.
     !!}
     private
     type(filterList), pointer :: filters => null()
  contains
     final     ::           anyDestructor
     procedure :: passes => anyPasses
  end type galacticFilterAny

  interface galacticFilterAny
     !!{
     Constructors for the any galactic filter class.
     !!}
     module procedure anyConstructorParameters
     module procedure anyConstructorInternal
  end interface galacticFilterAny

contains

  function anyConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``any'' galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (galacticFilterAny)                :: self
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
  end function anyConstructorParameters

  function anyConstructorInternal(filters) result(self)
    !!{
    Internal constructor for the ``any'' filter class.
    !!}
    implicit none
    type(galacticFilterAny)                        :: self
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
  end function anyConstructorInternal

  subroutine anyDestructor(self)
    !!{
    Destructor for the ``any'' galactic filter class.
    !!}
    implicit none
    type(galacticFilterAny), intent(inout) :: self
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
  end subroutine anyDestructor

  logical function anyPasses(self,node)
    !!{
    Apply a set of filters to a {\normalfont \ttfamily node} combined with ``any'' operations.
    !!}
    implicit none
    class(galacticFilterAny), intent(inout)         :: self
    type (treeNode         ), intent(inout), target :: node
    type (filterList       ), pointer               :: filter_

    ! Assume the node fails to pass initially. Iterate through filters and evaluate each one. If any one evaluates to true, exit
    ! the iteration.
    anyPasses =  .false.
    filter_  => self%filters
    do while (associated(filter_))
       anyPasses=filter_%filter_%passes(node)
       if (anyPasses) then
          filter_ => null()
       else
          filter_ => filter_%next
       end if
    end do
    return
  end function anyPasses
