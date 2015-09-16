!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a sequence of operators on merger trees.

  !# <mergerTreeOperator name="mergerTreeOperatorSequence" defaultThreadPrivate="yes">
  !#  <description>Provides a sequence of operators on merger trees.</description>
  !# </mergerTreeOperator>

  type, public :: operatorList
     class(mergerTreeOperatorClass), pointer :: operator_
     type (operatorList           ), pointer :: next      => null()
  end type operatorList

  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorSequence
     !% A sequence merger tree operator class.
     private
     type(operatorList), pointer :: operators
  contains
     final     ::             sequenceDestructor
     procedure :: operate  => sequenceOperate
     procedure :: finalize => sequenceFinalize
  end type mergerTreeOperatorSequence

  interface mergerTreeOperatorSequence
     !% Constructors for the sequence merger tree operator class.
     module procedure sequenceConstructorParameters
     module procedure sequenceConstructorInternal
  end interface mergerTreeOperatorSequence

contains

  function sequenceConstructorParameters(parameters)
    !% Constructor for the sequence merger tree operator class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type   (mergerTreeOperatorSequence)                :: sequenceConstructorParameters
    type   (inputParameters           ), intent(inout) :: parameters
    type   (operatorList              ), pointer       :: operator_
    integer                                            :: i

    !$omp critical(mergerTreeOperatorSequenceInitialize)
    operator_ => null()
    do i=1,parameters%copiesCount('mergerTreeOperatorMethod',zeroIfNotPresent=.true.)
       if (associated(operator_)) then
          allocate(operator_%next)
          operator_ => operator_%next
       else
          allocate(sequenceConstructorParameters%operators)
          operator_ => sequenceConstructorParameters%operators
       end if
       operator_%operator_ => mergerTreeOperator(parameters,i)
    end do
    !$omp end critical(mergerTreeOperatorSequenceInitialize)
    return
  end function sequenceConstructorParameters

  function sequenceConstructorInternal(operators)
    !% Internal constructor for the sequence merger tree operator class.
    implicit none
    type(mergerTreeOperatorSequence)                        :: sequenceConstructorInternal
    type(operatorList              ), target, intent(in   ) :: operators

    sequenceConstructorInternal%operators => operators
    return
  end function sequenceConstructorInternal

  elemental subroutine sequenceDestructor(self)
    !% Destructor for the merger tree operator function class.
    implicit none
    type(mergerTreeOperatorSequence), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine sequenceDestructor

  subroutine sequenceOperate(self,tree)
    !% Perform a sequence operation on a merger tree.
    implicit none
    class(mergerTreeOperatorSequence), intent(inout)         :: self
    type (mergerTree                ), intent(inout), target :: tree
    type (operatorList              ), pointer               :: operator_

    operator_ => self%operators
    do while (associated(operator_))
       call operator_%operator_%operate(tree)
       operator_ => operator_%next
    end do
    return
  end subroutine sequenceOperate

  subroutine sequenceFinalize(self)
    !% Perform a finalization on a sequence of operators on a merger tree.
    implicit none
    class(mergerTreeOperatorSequence), intent(inout)         :: self
    type (operatorList              ), pointer               :: operator_

    operator_ => self%operators
    do while (associated(operator_))
       call operator_%operator_%finalize()
       operator_ => operator_%next
    end do
    return
  end subroutine sequenceFinalize
