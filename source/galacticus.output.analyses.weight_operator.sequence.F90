!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains a module which implements a sequence output analysis weight operator class.

  type, public :: weightOperatorList
     class(outputAnalysisWeightOperatorClass), pointer :: operator_ => null()
     type (weightOperatorList               ), pointer :: next      => null()
  end type weightOperatorList

  !# <outputAnalysisWeightOperator name="outputAnalysisWeightOperatorSequence" defaultThreadPrivate="yes">
  !#  <description>A sequence output analysis weight operator class.</description>
  !# </outputAnalysisWeightOperator>
  type, extends(outputAnalysisWeightOperatorClass) :: outputAnalysisWeightOperatorSequence
     !% A sequence output weight operator class.
     private
     type(weightOperatorList), pointer :: operators
   contains
     !@ <objectMethods>
     !@   <object>outputAnalysisWeightOperatorSequence</object>
     !@   <objectMethod>
     !@     <method>prepend</method>
     !@     <arguments>\textcolor{red}{\textless class(outputAnalysisWeightOperatorClass)\textgreater} operator\_\argin</arguments>
     !@     <type>\void</type>
     !@     <description>Prepend an operator to a sequence of weight operators.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::             sequenceDestructor
     procedure :: operate  => sequenceOperate
     procedure :: prepend  => sequencePrepend
     procedure :: deepCopy => sequenceDeepCopy
  end type outputAnalysisWeightOperatorSequence

  interface outputAnalysisWeightOperatorSequence
     !% Constructors for the ``sequence'' output analysis class.
     module procedure sequenceConstructorParameters
     module procedure sequenceConstructorInternal
  end interface outputAnalysisWeightOperatorSequence

contains

  function sequenceConstructorParameters(parameters) result (self)
    !% Constructor for the ``sequence'' output analysis weight operator class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type   (outputAnalysisWeightOperatorSequence)                :: self
    type   (inputParameters                     ), intent(inout) :: parameters
    type   (weightOperatorList                  ), pointer       :: operator_
    integer                                                      :: i

    self     %operators => null()
    operator_           => null()
    do i=1,parameters%copiesCount('outputAnalysisWeightOperatorMethod',zeroIfNotPresent=.true.)
       if (associated(operator_)) then
          allocate(operator_%next)
          operator_ => operator_%next
       else
          allocate(self%operators)
          operator_ => self%operators
       end if
       operator_%operator_ => outputAnalysisWeightOperator(parameters,i)
    end do
    return
  end function sequenceConstructorParameters

  function sequenceConstructorInternal(operators) result (self)
    !% Internal constructor for the sequence merger tree normalizer class.
    implicit none
    type(outputAnalysisWeightOperatorSequence)                        :: self
    type(weightOperatorList                  ), target, intent(in   ) :: operators

    self%operators => operators
    return
  end function sequenceConstructorInternal

  elemental subroutine sequenceDestructor(self)
    !% Destructor for the sequence weight operator class.
    implicit none
    type(outputAnalysisWeightOperatorSequence), intent(inout) :: self
    type(weightOperatorList                  ), pointer       :: operator_, operatorNext

    if (associated(self%operators)) then
       operator_ => self%operators
       do while (associated(operator_))
          operatorNext => operator_%next
          deallocate(operator_%operator_)
          deallocate(operator_          )
          operator_ => operatorNext
       end do
    end if
    return
  end subroutine sequenceDestructor

  double precision function sequenceOperate(self,weightValue,node,propertyValue,propertyValueIntrinsic,propertyType,propertyQuantity,outputIndex)
    !% Implement an sequence output analysis weight operator.
    implicit none
    class           (outputAnalysisWeightOperatorSequence), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                      , intent(in   ) :: weightValue  , propertyValueIntrinsic, &
         &                                                                   propertyValue
    integer                                               , intent(in   ) :: propertyType , propertyQuantity
    integer         (c_size_t                            ), intent(in   ) :: outputIndex
    type            (weightOperatorList                  ), pointer       :: operator_

    sequenceOperate =  weightValue
    operator_       => self%operators
    do while (associated(operator_))
       sequenceOperate =  operator_%operator_%operate(sequenceOperate,node,propertyValue,propertyValueIntrinsic,propertyType,propertyQuantity,outputIndex)
       operator_       => operator_%next
    end do
    return
  end function sequenceOperate
  
  subroutine sequencePrepend(self,operator_)
    !% Prepend an operator to the sequence.
    implicit none
    class(outputAnalysisWeightOperatorSequence), intent(inout)          :: self
    class(outputAnalysisWeightOperatorClass   ), intent(in   ), target  :: operator_
    type (weightOperatorList                  )               , pointer :: operatorNew
    
    allocate(operatorNew)
    operatorNew%operator_ => operator_
    operatorNew%next      => self       %operators
    self       %operators => operatorNew
    return
  end subroutine sequencePrepend

  subroutine sequenceDeepCopy(self,destination)
    !% Perform a deep copy for the {\normalfont \ttfamily sequence} output analysis weight operator class.
    use Galacticus_Error
    implicit none
    class(outputAnalysisWeightOperatorSequence), intent(inout) :: self
    class(outputAnalysisWeightOperatorClass   ), intent(  out) :: destination
    type (weightOperatorList                  ), pointer       :: operator_   , operatorDestination_, &
         &                                                        operatorNew_

    call self%outputAnalysisWeightOperatorClass%deepCopy(destination)
    select type (destination)
    type is (outputAnalysisWeightOperatorSequence)
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
          call operator_%operator_%deepCopy(operatorNew_%operator_)
          operator_ => operator_%next
       end do       
    class default
       call Galacticus_Error_Report('destination and source types do not match'//{introspection:location})
    end select
    return
  end subroutine sequenceDeepCopy
