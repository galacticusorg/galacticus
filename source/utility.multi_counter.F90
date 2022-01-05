!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
  Contains a module which implements multi-counters - objects which iterate over all combinations of an arbitary number of
  counters, each with an arbitrary range.
  !!}

module Multi_Counters
  !!{
  Implements multi-counters - objects which iterate over all combinations of an arbitary number of counters, each with an
  arbitrary range.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_size_t
  implicit none
  private
  public :: multiCounter

  type :: multiCounter
     !!{
     Class providing multi-counters - objects which iterate over all combinations of an arbitary number of counters, each with
     an arbitrary range.
     !!}
     integer(c_size_t), allocatable, dimension(:) :: ranges, values
   contains
     !![
     <methods>
       <method description="Reset the multi-counter back to its initial count state, such that the next increment will return the first count." method="reset" />
       <method description="Return the number of counters configured in this multi-counter." method="count" />
       <method description="Append a new counter to the multi-counter, with the specified {\normalfont \ttfamily range}." method="append" />
       <method description="Increment the state of the multi-counter. Return {\normalfont \ttfamily .false.} if incrementing was not possible (i.e. counter was in the final state), {\normalfont \ttfamily .true.} otherwise." method="increment" />
       <method description="Return {\normalfont \ttfamily .true.} if the counter is in its final state, {\normalfont \ttfamily .false.} otherwise." method="isFinal" />
       <method description="Return the state of the {\normalfont \ttfamily i}$^\mathrm{th}$ counter." method="state" />
     </methods>
     !!]
     final     ::              multiCounterDestructor
     procedure :: append    => multiCounterAppend
     procedure :: reset     => multiCounterReset
     procedure :: count     => multiCounterCount
     procedure :: state     => multiCounterState
     procedure :: increment => multiCounterIncrement
     procedure :: isFinal   => multiCounterIsFinal
  end type multiCounter

  interface multiCounter
     !!{
     Constructors for multi-counters.
     !!}
     module procedure multiCounterConstructor
  end interface multiCounter

contains

  function multiCounterConstructor(ranges) result (self)
    !!{
    Constructor for multi-counters where the ranges are provided.
    !!}
    use :: Galacticus_Error , only : Galacticus_Error_Report
    use :: Memory_Management, only : allocateArray
    implicit none
    type   (multiCounter)                              :: self
    integer(c_size_t    ), intent(in   ), dimension(:) :: ranges

    ! Validate ranges.
    if (any(ranges < 1_c_size_t)) call Galacticus_Error_Report('ranges must be positive'//{introspection:location})
    ! Build the object.
    call allocateArray(self%ranges,shape(ranges))
    call allocateArray(self%values,shape(ranges))
    self%ranges=ranges
    call self%reset()
    return
  end function multiCounterConstructor

  subroutine multiCounterDestructor(self)
    !!{
    Destroy a multi-counter object.
    !!}
    use :: Memory_Management, only : deallocateArray
    type(multiCounter), intent(inout) :: self

    if (allocated(self%ranges)) call deallocateArray(self%ranges)
    if (allocated(self%values)) call deallocateArray(self%values)
    return
  end subroutine multiCounterDestructor

  subroutine multiCounterReset(self)
    !!{
    Reset the state of the multi-counter.
    !!}
    implicit none
    class(multiCounter), intent(inout) :: self

    self%values=0_c_size_t
    return
  end subroutine multiCounterReset

  function multiCounterCount(self)
    !!{
    Return the number of counters in the multi-counter.
    !!}
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    integer(c_size_t    )                :: multiCounterCount
    class  (multiCounter), intent(inout) :: self

    ! Validate input.
    if (.not.allocated(self%ranges)) call Galacticus_Error_Report('no counters defined'//{introspection:location})
    ! Return the count.
    multiCounterCount=size(self%ranges)
    return
  end function multiCounterCount

  function multiCounterState(self,i)
    !!{
    Return the state of the multi-counter.
    !!}
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    integer(c_size_t    )                :: multiCounterState
    class  (multiCounter), intent(inout) :: self
    integer(c_size_t    ), intent(in   ) :: i

    ! Validate input.
    if (               .not.allocated(self%ranges)) call Galacticus_Error_Report('no counters defined'//{introspection:location})
    if (i < 1 .or. i >           size(self%ranges)) call Galacticus_Error_Report('out of range'//{introspection:location})
    ! Return the state.
    multiCounterState=self%values(i)
    return
  end function multiCounterState

  subroutine multiCounterAppend(self,range)
    !!{
    Append a new counter with the given {\normalfont \ttfamily range}.
    !!}
    use :: Galacticus_Error , only : Galacticus_Error_Report
    use :: Memory_Management, only : allocateArray          , deallocateArray
    implicit none
    class  (multiCounter), intent(inout)               :: self
    integer(c_size_t    ), intent(in   )               :: range
    integer(c_size_t    ), allocatable  , dimension(:) :: rangesTmp

    ! Validate range.
    if (range < 1_c_size_t) call Galacticus_Error_Report('range must be positive'//{introspection:location})
    ! Expand the range.
    if (allocated(self%ranges)) then
       call move_alloc(self%ranges,rangesTmp)
       call allocateArray(self%ranges,shape(rangesTmp)+1)
       self%ranges(1:size(rangesTmp))=rangesTmp
       call deallocateArray(rangesTmp                     )
       call deallocateArray(self%values                   )
       call   allocateArray(self%values,shape(self%ranges))
    else
       call allocateArray(self%ranges,[1])
       call allocateArray(self%values,[1])
    end if
    self%ranges(size(self%ranges))=range
    call self%reset()
    return
  end subroutine multiCounterAppend

  logical function multiCounterIncrement(self)
    !!{
    Increment a multi-counter. Return true if increment was possible, false otherwise.
    !!}
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class  (multiCounter), intent(inout) :: self
    integer(c_size_t    )                :: i

    ! Validate.
    if (.not.allocated(self%ranges)) call Galacticus_Error_Report('no counters defined'//{introspection:location})
    ! Assume incrementing was possible.
    multiCounterIncrement=.true.
    ! Increment.
    if (all(self%values == 0_c_size_t)) then
       ! Counter is in initial state, put into first state.
       self%values=1_c_size_t
    else
       do i=1,size(self%ranges)
          self%values(i)=self%values(i)+1_c_size_t
          if (self%values(i) > self%ranges(i)) then
             self%values(i)=1_c_size_t
          else
             return
          end if
       end do
       ! All values exceeded their range - the counter has been fully iterated over.
       self%values          =0_c_size_t
       multiCounterIncrement=.false.
    end if
    return
  end function multiCounterIncrement

  logical function multiCounterIsFinal(self)
    !!{
    Return true if a multi-counter is in its final state, false otherwise.
    !!}
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class(multiCounter), intent(in   ) :: self

    ! Validate.
    if (.not.allocated(self%ranges)) call Galacticus_Error_Report('no counters defined'//{introspection:location})
    ! Determine if in final state.
    multiCounterIsFinal=all(self%values == self%ranges)
    return
  end function multiCounterIsFinal

end module Multi_Counters
