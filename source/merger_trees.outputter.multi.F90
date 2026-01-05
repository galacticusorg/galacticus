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
  Implements a merger trees outputter class which combines multiple other outputters.
  !!}

  type, public :: multiOutputterList
     class(mergerTreeOutputterClass), pointer :: outputter_ => null()
     type (multiOutputterList      ), pointer :: next       => null()
  end type multiOutputterList

  !![
  <mergerTreeOutputter name="mergerTreeOutputterMulti">
   <description>A merger tree outputter which combines multiple other outputters.</description>
   <linkedList type="multiOutputterList" variable="outputters" next="next" object="outputter_" objectType="mergerTreeOutputterClass"/>
  </mergerTreeOutputter>
  !!]
  type, extends(mergerTreeOutputterClass) :: mergerTreeOutputterMulti
     !!{
     Implementation of a merger tree outputter which combines multiple other outputters.
     !!}
     private
     type(multiOutputterList), pointer :: outputters => null()
   contains
     final     ::               multiDestructor
     procedure :: outputTree => multiOutputTree
     procedure :: outputNode => multiOutputNode
     procedure :: finalize   => multiFinalize
     procedure :: reduce     => multiReduce
  end type mergerTreeOutputterMulti

  interface mergerTreeOutputterMulti
     !!{
     Constructors for the \refClass{mergerTreeOutputterMulti} merger tree outputter.
     !!}
     module procedure multiConstructorParameters
     module procedure multiConstructorInternal
  end interface mergerTreeOutputterMulti

contains

  function multiConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeOutputterMulti} merger tree outputter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (mergerTreeOutputterMulti)                :: self
    type   (inputParameters         ), intent(inout) :: parameters
    type   (multiOutputterList      ), pointer       :: outputter_
    integer                                          :: i

    self      %outputters => null()
    outputter_           => null()
    do i=1,parameters%copiesCount('mergerTreeOutputter',zeroIfNotPresent=.true.)
       if (associated(outputter_)) then
          allocate(outputter_%next)
          outputter_ => outputter_%next
       else
          allocate(self%outputters)
          outputter_ => self%outputters
       end if
       !![
       <objectBuilder class="mergerTreeOutputter" name="outputter_%outputter_" source="parameters" copy="i" />
       !!]
    end do
    !![
    <inputParametersValidate source="parameters" multiParameters="mergerTreeOutputter"/>
    !!]
    return
  end function multiConstructorParameters

  function multiConstructorInternal(outputters) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeOutputterMulti} outputter class.
    !!}
    implicit none
    type(mergerTreeOutputterMulti)                        :: self
    type(multiOutputterList      ), target, intent(in   ) :: outputters
    type(multiOutputterList      ), pointer               :: outputter_

    self      %outputters => outputters
    outputter_            => outputters
    do while (associated(outputter_))
       !![
       <referenceCountIncrement owner="outputter_" object="outputter_"/>
       !!]
       outputter_ => outputter_%next
    end do
    return
  end function multiConstructorInternal

  subroutine multiDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeOutputterMulti} outputter class.
    !!}
    implicit none
    type(mergerTreeOutputterMulti), intent(inout) :: self
    type(multiOutputterList      ), pointer       :: outputter_, outputterNext

    if (associated(self%outputters)) then
       outputter_ => self%outputters
       do while (associated(outputter_))
          outputterNext => outputter_%next
          !![
          <objectDestructor name="outputter_%outputter_"/>
          !!]
          deallocate(outputter_)
          outputter_ => outputterNext
       end do
    end if
    return
  end subroutine multiDestructor

  subroutine multiOutputTree(self,tree,indexOutput,time)
    !!{
    Output from all outputters.
    !!}
    implicit none
    class           (mergerTreeOutputterMulti), intent(inout)         :: self
    type            (mergerTree              ), intent(inout), target :: tree
    integer         (c_size_t                ), intent(in   )         :: indexOutput
    double precision                          , intent(in   )         :: time
    type            (multiOutputterList      ), pointer               :: outputter_

    outputter_ => self%outputters
    do while (associated(outputter_))
       call outputter_%outputter_%outputTree(tree,indexOutput,time)
       outputter_ => outputter_%next
    end do
    return
  end subroutine multiOutputTree

  subroutine multiOutputNode(self,node,indexOutput)
    !!{
    Output from all outputters.
    !!}
    implicit none
    class  (mergerTreeOutputterMulti), intent(inout) :: self
    type   (treeNode                ), intent(inout) :: node
    integer(c_size_t                ), intent(in   ) :: indexOutput
    type   (multiOutputterList      ), pointer       :: outputter_

    outputter_ => self%outputters
    do while (associated(outputter_))
       call outputter_%outputter_%outputNode(node,indexOutput)
       outputter_ => outputter_%next
    end do
    return
  end subroutine multiOutputNode

  subroutine multiFinalize(self)
    !!{
    Finalize all outputters.
    !!}
    implicit none
    class(mergerTreeOutputterMulti), intent(inout) :: self
    type (multiOutputterList      ), pointer       :: outputter_

    outputter_ => self%outputters
    do while (associated(outputter_))
       call outputter_%outputter_%finalize()
       outputter_ => outputter_%next
    end do
    return
  end subroutine multiFinalize

  subroutine multiReduce(self,reduced)
    !!{
    Reduce over the outputter.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(mergerTreeOutputterMulti), intent(inout) :: self
    class(mergerTreeOutputterClass), intent(inout) :: reduced
    type (multiOutputterList      ), pointer       :: outputter_, outputterReduced_

    select type (reduced)
    type is (mergerTreeOutputterMulti)
       outputter_        => self   %outputters
       outputterReduced_ => reduced%outputters
       do while (associated(outputter_))
          call outputter_%outputter_%reduce(outputterReduced_%outputter_)
          outputter_        => outputter_%next
          outputterReduced_ => outputterReduced_%next
       end do
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine multiReduce
