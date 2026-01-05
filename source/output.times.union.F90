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

  type, public :: outputTimesUnionList
     class(outputTimesClass    ), pointer :: outputTimes_ => null()
     type (outputTimesUnionList), pointer :: next         => null()
  end type outputTimesUnionList

  !![
  <outputTimes name="outputTimesUnion">
    <description>      
      An output times class which simply constructs the union of output times from a set of other \refClass{outputTimesClass} objects.
    </description>
    <linkedList type="outputTimesUnionList" variable="outputTimesUnion_" next="next" object="outputTimes_" objectType="outputTimesClass"/>
  </outputTimes>
  !!]
  type, extends(outputTimesList) :: outputTimesUnion
     !!{
     Implementation of an output times class which simply constructs the union of output times from a set of other \refClass{outputTimesClass} objects.
     !!}
     private
     type(outputTimesUnionList), pointer :: outputTimesUnion_=> null()
   contains
     !![
     <methods>
       <method method="initialize" description="Initialize the set of output times."/>
     </methods>
     !!]
     final     ::               unionDestructor
     procedure :: initialize => unionInitialize
  end type outputTimesUnion

  interface outputTimesUnion
     !!{
     Constructors for the \refClass{outputTimesUnion} output times class.
     !!}
     module procedure unionConstructorParameters
     module procedure unionConstructorInternal
  end interface outputTimesUnion

contains

  function unionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputTimesUnion} output times class which takes a parameter set as input.
    !!}
    use :: Input_Parameters , only : inputParameters
    implicit none
    type   (outputTimesUnion    )                :: self
    type   (inputParameters     ), intent(inout) :: parameters
    type   (outputTimesUnionList), pointer       :: outputTimesUnion_
    integer                                      :: i

    self%outputTimesUnion_ => null()
    outputTimesUnion_      => null()
    do i=1,parameters%copiesCount('outputTimesUnion',zeroIfNotPresent=.true.)
       if (associated(outputTimesUnion_)) then
          allocate(outputTimesUnion_%next)
          outputTimesUnion_ => outputTimesUnion_%next
       else
          allocate(self%outputTimesUnion_)
          outputTimesUnion_ => self             %outputTimesUnion_
       end if
       outputTimesUnion_%outputTimes_ => outputTimes(parameters,i)
    end do
    !![
    <inputParametersValidate source="parameters" multiParameters="outputTimesUnion"/>
    !!]
    call self%initialize()
    return
  end function unionConstructorParameters

  function unionConstructorInternal(outputTimesUnion_) result(self)
    !!{
    Constructor for the \refClass{outputTimesUnion} output times class which takes a parameter set as input.
    !!}
    implicit none
    type(outputTimesUnion    )                        :: self
    type(outputTimesUnionList), target, intent(in   ) :: outputTimesUnion_
    !![
    <constructorAssign variables="*outputTimesUnion_"/>
    !!]

    call self%initialize()
    return
  end function unionConstructorInternal

  subroutine unionDestructor(self)
    !!{
    Destructor for the \refClass{outputTimesUnion} output times class.
    !!}
    implicit none
    type(outputTimesUnion    ), intent(inout) :: self
    type(outputTimesUnionList), pointer       :: outputTimesUnion_, outputTimesUnionNext

    if (associated(self%outputTimesUnion_)) then
       outputTimesUnion_ => self%outputTimesUnion_
       do while (associated(outputTimesUnion_))
          outputTimesUnionNext => outputTimesUnion_   %next
          deallocate(outputTimesUnion_%outputTimes_)
          deallocate(outputTimesUnion_             )
          outputTimesUnion_    => outputTimesUnionNext
       end do
    end if
    return
  end subroutine unionDestructor

  subroutine unionInitialize(self)
    !!{
    Initialize the list of times and redshifts.
    !!}
    use :: Sorting, only : sortIndex
    implicit none
    class           (outputTimesUnion    ), intent(inout)             :: self
    type            (outputTimesUnionList), pointer                   :: outputTimesUnion_
    double precision                      , allocatable, dimension(:) :: times            , redshifts
    integer         (c_size_t            ), allocatable, dimension(:) :: order
    integer         (c_size_t            )                            :: countOutputs     , i
    
    if (associated(self%outputTimesUnion_)) then
       ! Count the number of output times.
       countOutputs=0_c_size_t
       outputTimesUnion_ => self%outputTimesUnion_
       do while (associated(outputTimesUnion_))
          countOutputs      =  +                               countOutputs   &
               &               +outputTimesUnion_%outputTimes_%count       ()
          outputTimesUnion_ => outputTimesUnion_%next
       end do
       ! Allocate arrays.
       allocate(     times    (countOutputs))
       allocate(self%times    (countOutputs))
       allocate(     redshifts(countOutputs))
       allocate(self%redshifts(countOutputs))
       ! Extract all times and redshifts.
       countOutputs      =  0_c_size_t
       outputTimesUnion_ => self%outputTimesUnion_
       do while (associated(outputTimesUnion_))
          do i=1_c_size_t,outputTimesUnion_%outputTimes_%count()
             countOutputs              =+countOutputs &
                  &                     +1_c_size_t
             times       (countOutputs)= outputTimesUnion_%outputTimes_%time    (i)
             redshifts   (countOutputs)= outputTimesUnion_%outputTimes_%redshift(i)
          end do
          outputTimesUnion_ => outputTimesUnion_%next
       end do
       ! Sort the times in order.
       order=sortIndex(times)
       ! Assign to self in order.
       do i=1_c_size_t,size(order)
          self%times    (i)=times    (order(i))
          self%redshifts(i)=redshifts(order(i))
       end do
    end if
    return
  end subroutine unionInitialize
