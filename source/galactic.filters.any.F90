!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

!% Contains a module which implements a galactic filter class which is the ``any'' combination of a set of other filters.

  !# <galacticFilter name="galacticFilterAny" defaultThreadPrivate="yes">
  !#  <description>A galactic filter class which is the ``any'' combination of a set of other filters.</description>
  !# </galacticFilter>

  type, extends(galacticFilterClass) :: galacticFilterAny
     !% A galactic filter class which is the ``any'' combination of a set of other filters.
     private
     type(filterList), pointer :: filters
  contains
     final     ::            anyDestructor
     procedure :: passes  => anyPasses
  end type galacticFilterAny

  interface galacticFilterAny
     !% Constructors for the any galactic filter class.
     module procedure anyConstructorParameters
     module procedure anyConstructorInternal
  end interface galacticFilterAny

contains

  function anyConstructorParameters(parameters)
    !% Constructor for the ``any'' galactic filter class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type   (galacticFilterAny)                :: anyConstructorParameters
    type   (inputParameters  ), intent(inout) :: parameters
    type   (filterList       ), pointer       :: filter_
    integer                                   :: i

    !$omp critical(galacticFilterAnyInitialize)
    anyConstructorParameters%filters => null()
    filter_                          => null()
    do i=1,parameters%copiesCount('galacticFilterMethod',zeroIfNotPresent=.true.)
       if (associated(filter_)) then
          allocate(filter_%next)
          filter_ => filter_%next
       else
          allocate(anyConstructorParameters%filters)
          filter_ => anyConstructorParameters%filters
       end if
       filter_%filter_ => galacticFilter(parameters,i)
    end do
    !$omp end critical(galacticFilterAnyInitialize)
    return
  end function anyConstructorParameters

  function anyConstructorInternal(filters)
    !% Internal constructor for the ``any'' filter class.
    implicit none
    type(galacticFilterAny)                        :: anyConstructorInternal
    type(filterList       ), target, intent(in   ) :: filters

    anyConstructorInternal%filters => filters
    return
  end function anyConstructorInternal

  elemental subroutine anyDestructor(self)
    !% Destructor for the ``any'' galactic filter class.
    implicit none
    type(galacticFilterAny), intent(inout) :: self
    type(filterList       ), pointer       :: filter_, filterNext

    if (associated(self%filters)) then
       filter_ => self%filters
       do while (associated(filter_))
          filterNext => filter_%next
          deallocate(filter_%filter_)
          deallocate(filter_          )
          filter_ => filterNext
       end do
    end if
    return
  end subroutine anyDestructor

  logical function anyPasses(self,node)
    !% Apply a set of filters to a {\normalfont \ttfamily node} combined with ``any'' operations.
    use Galacticus_Nodes
    implicit none
    class(galacticFilterAny), intent(inout) :: self
    type (treeNode        ), intent(inout) :: node
    type (filterList      ), pointer       :: filter_

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
