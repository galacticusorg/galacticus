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

!% Contains a module which implements a galactic filter class which is the ``all'' combination of a set of other filters.

  !# <galacticFilter name="galacticFilterAll" defaultThreadPrivate="yes">
  !#  <description>A galactic filter class which is the ``all'' combination of a set of other filters.</description>
  !# </galacticFilter>

  type, extends(galacticFilterClass) :: galacticFilterAll
     !% A galactic filter class which is the ``all'' combination of a set of other filters.
     private
     type(filterList), pointer :: filters
  contains
     final     ::            allDestructor
     procedure :: passes  => allPasses
  end type galacticFilterAll

  interface galacticFilterAll
     !% Constructors for the all galactic filter class.
     module procedure allConstructorParameters
     module procedure allConstructorInternal
  end interface galacticFilterAll

contains

  function allConstructorParameters(parameters)
    !% Constructor for the ``all'' galactic filter class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type   (galacticFilterAll)                :: allConstructorParameters
    type   (inputParameters  ), intent(inout) :: parameters
    type   (filterList       ), pointer       :: filter_
    integer                                   :: i

    !$omp critical(galacticFilterAllInitialize)
    allConstructorParameters%filters => null()
    filter_                          => null()
    do i=1,parameters%copiesCount('galacticFilterMethod',zeroIfNotPresent=.true.)
       if (associated(filter_)) then
          allocate(filter_%next)
          filter_ => filter_%next
       else
          allocate(allConstructorParameters%filters)
          filter_ => allConstructorParameters%filters
       end if
       filter_%filter_ => galacticFilter(parameters,i)
    end do
    !$omp end critical(galacticFilterAllInitialize)
    return
  end function allConstructorParameters

  function allConstructorInternal(filters)
    !% Internal constructor for the ``all'' filter class.
    implicit none
    type(galacticFilterAll)                        :: allConstructorInternal
    type(filterList       ), target, intent(in   ) :: filters

    allConstructorInternal%filters => filters
    return
  end function allConstructorInternal

  elemental subroutine allDestructor(self)
    !% Destructor for the all galactic filter class.
    implicit none
    type(galacticFilterAll), intent(inout) :: self
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
  end subroutine allDestructor

  logical function allPasses(self,node)
    !% Apply a set of filters to a {\normalfont \ttfamily node} combined with ``all'' operations.
    use Galacticus_Nodes
    implicit none
    class(galacticFilterAll), intent(inout) :: self
    type (treeNode         ), intent(inout) :: node
    type (filterList       ), pointer       :: filter_

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
