!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

!% Contains a module which implements an N-body data operator which applies a sequence of other operators.

  type, public :: nbodyOperatorList
     class(nbodyOperatorClass), pointer :: operator_
     type (nbodyOperatorList ), pointer :: next     => null()
  end type nbodyOperatorList

  !# <nbodyOperator name="nbodyOperatorSequence">
  !#  <description>An N-body data operator which applies a sequence of other operators.</description>
  !# </nbodyOperator>
  type, extends(nbodyOperatorClass) :: nbodyOperatorSequence
     !% An N-body data operator which applies a sequence of other operators.
     private
     type(nbodyOperatorList), pointer :: operators => null()
   contains
     final     ::             sequenceDestructor
     procedure :: operate  => sequenceOperate
     procedure :: deepCopy => sequenceDeepCopy
  end type nbodyOperatorSequence

  interface nbodyOperatorSequence
     !% Constructors for the ``sequence'' N-body operator class.
     module procedure sequenceConstructorParameters
     module procedure sequenceConstructorInternal
  end interface nbodyOperatorSequence

contains

  function sequenceConstructorParameters(parameters) result (self)
    !% Constructor for the ``sequence'' N-body operator class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nbodyOperatorSequence)                :: self
    type   (inputParameters      ), intent(inout) :: parameters
    type   (nbodyOperatorList    ), pointer       :: operator_
    integer                                       :: i

    self     %operators => null()
    operator_           => null()
    do i=1,parameters%copiesCount('nbodyOperatorMethod',zeroIfNotPresent=.true.)
       if (associated(operator_)) then
          allocate(operator_%next)
          operator_ => operator_%next
       else
          allocate(self%operators)
          operator_ => self%operators
       end if
       !# <objectBuilder class="nbodyOperator" name="operator_%operator_" source="parameters" copy="i" />
    end do
    return
  end function sequenceConstructorParameters

  function sequenceConstructorInternal(operators) result (self)
    !% Internal constructor for the ``sequence'' N-body operator class.
    implicit none
    type(nbodyOperatorSequence)                        :: self
    type(nbodyOperatorList    ), target, intent(in   ) :: operators
    type(nbodyOperatorList    ), pointer               :: operator_

    self     %operators => operators
    operator_           => operators
    do while (associated(operator_))
       !# <referenceCountIncrement owner="operator_" object="operator_"/>
       operator_ => operator_%next
    end do
    return
  end function sequenceConstructorInternal

  subroutine sequenceDestructor(self)
    !% Destructor for the sequence N-body operator class.
    implicit none
    type(nbodyOperatorSequence), intent(inout) :: self
    type(nbodyOperatorList    ), pointer       :: operator_, operatorNext

    if (associated(self%operators)) then
       operator_ => self%operators
       do while (associated(operator_))
          operatorNext => operator_%next
          !# <objectDestructor name="operator_%operator_"/>
          deallocate(operator_)
          operator_ => operatorNext
       end do
    end if
    return
  end subroutine sequenceDestructor

  subroutine sequenceOperate(self,simulation)
    !% Apply a sequence of N-body simulation operators.
    implicit none
    class(nbodyOperatorSequence), intent(inout) :: self
    type (nBodyData            ), intent(inout) :: simulation
    type (nbodyOperatorList    ), pointer       :: operator_

    operator_       => self%operators
    do while (associated(operator_))
       call operator_%operator_%operate(simulation)
       operator_ => operator_%next
    end do
    return
  end subroutine sequenceOperate

  subroutine sequenceDeepCopy(self,destination)
    !% Perform a deep copy for the {\normalfont \ttfamily sequence} N-body operator class.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class(nbodyOperatorSequence), intent(inout) :: self
    class(nbodyOperatorClass   ), intent(inout) :: destination
    type (nbodyOperatorList    ), pointer       :: operator_   , operatorDestination_, &
         &                                         operatorNew_

    call self%nbodyOperatorClass%deepCopy(destination)
    select type (destination)
    type is (nbodyOperatorSequence)
       destination%operators => null          ()
       operatorDestination_  => null          ()
       operator_             => self%operators
       do while (associated(operator_))
          allocate(operatorNew_)
          if (associated(operatorDestination_)) then
             operatorDestination_%next       => operatorNew_
             operatorDestination_            => operatorNew_
          else
             destination          %operators => operatorNew_
             operatorDestination_            => operatorNew_
          end if
          allocate(operatorNew_%operator_,mold=operator_%operator_)
          !# <deepCopy source="operator_%operator_" destination="operatorNew_%operator_"/>
          operator_ => operator_%next
       end do
    class default
       call Galacticus_Error_Report('destination and source types do not match'//{introspection:location})
    end select
    return
  end subroutine sequenceDeepCopy
