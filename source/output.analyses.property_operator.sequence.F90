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
Implements a sequence output analysis property operator class.
!!}

  type, public :: propertyOperatorList
     class(outputAnalysisPropertyOperatorClass), pointer :: operator_ => null()
     type (propertyOperatorList               ), pointer :: next      => null()
  end type propertyOperatorList

  !![
  <outputAnalysisPropertyOperator name="outputAnalysisPropertyOperatorSequence">
   <description>A sequence output analysis property operator class.</description>
   <linkedList type="propertyOperatorList" variable="operators" next="next" object="operator_" objectType="outputAnalysisPropertyOperatorClass"/>
  </outputAnalysisPropertyOperator>
  !!]
  type, extends(outputAnalysisPropertyOperatorClass) :: outputAnalysisPropertyOperatorSequence
     !!{
     A sequence output property operator class.
     !!}
     private
     type(propertyOperatorList), pointer :: operators => null()
   contains
     !![
     <methods>
       <method description="Prepend an operator to a sequence of property operators." method="prepend" />
       <method description="Apppend an operator to a sequence of property operators." method="apppend" />
     </methods>
     !!]
     final     ::            sequenceDestructor
     procedure :: operate => sequenceOperate
     procedure :: prepend => sequencePrepend
     procedure :: apppend => sequenceApppend
  end type outputAnalysisPropertyOperatorSequence

  interface outputAnalysisPropertyOperatorSequence
     !!{
     Constructors for the {\normalfont \ttfamily sequence} output analysis class.
     !!}
     module procedure sequenceConstructorParameters
     module procedure sequenceConstructorInternal
  end interface outputAnalysisPropertyOperatorSequence

contains

  function sequenceConstructorParameters(parameters) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily sequence} output analysis property operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (outputAnalysisPropertyOperatorSequence)                :: self
    type   (inputParameters                       ), intent(inout) :: parameters
    type   (propertyOperatorList                  ), pointer       :: operator_
    integer                                                        :: i

    self     %operators => null()
    operator_           => null()
    do i=1,parameters%copiesCount('outputAnalysisPropertyOperator',zeroIfNotPresent=.true.)
       if (associated(operator_)) then
          allocate(operator_%next)
          operator_ => operator_%next
       else
          allocate(self%operators)
          operator_ => self%operators
       end if
       !![
       <objectBuilder class="outputAnalysisPropertyOperator" name="operator_%operator_" source="parameters" copy="i" />
       !!]
    end do
    !![
    <inputParametersValidate source="parameters" multiParameters="outputAnalysisPropertyOperator"/>
    !!]
    return
  end function sequenceConstructorParameters

  function sequenceConstructorInternal(operators) result (self)
    !!{
    Internal constructor for the sequence merger tree normalizer class.
    !!}
    implicit none
    type(outputAnalysisPropertyOperatorSequence)                        :: self
    type(propertyOperatorList                  ), target, intent(in   ) :: operators
    type(propertyOperatorList                  ), pointer               :: operator_

    self     %operators => operators
    operator_           => operators
    do while (associated(operator_))
       !![
       <referenceCountIncrement owner="operator_" object="operator_"/>
       !!]
       operator_ => operator_%next
    end do
    return
  end function sequenceConstructorInternal

  subroutine sequenceDestructor(self)
    !!{
    Destructor for the sequence property operator class.
    !!}
    implicit none
    type(outputAnalysisPropertyOperatorSequence), intent(inout) :: self
    type(propertyOperatorList                  ), pointer       :: operator_, operatorNext

    if (associated(self%operators)) then
       operator_ => self%operators
       do while (associated(operator_))
          operatorNext => operator_%next
          !![
          <objectDestructor name="operator_%operator_"/>
          !!]
          deallocate(operator_)
          operator_ => operatorNext
       end do
    end if
    return
  end subroutine sequenceDestructor

  double precision function sequenceOperate(self,propertyValue,node,propertyType,outputIndex)
    !!{
    Implement an sequence output analysis property operator.
    !!}
    implicit none
    class           (outputAnalysisPropertyOperatorSequence   ), intent(inout)           :: self
    double precision                                           , intent(in   )           :: propertyValue
    type            (treeNode                                 ), intent(inout), optional :: node
    type            (enumerationOutputAnalysisPropertyTypeType), intent(inout), optional :: propertyType
    integer         (c_size_t                                 ), intent(in   ), optional :: outputIndex
    type            (propertyOperatorList                     ), pointer                 :: operator_

    sequenceOperate =  propertyValue
    operator_       => self%operators
    do while (associated(operator_))
       sequenceOperate =  operator_%operator_%operate(sequenceOperate,node,propertyType,outputIndex)
       operator_       => operator_%next
    end do
    return
  end function sequenceOperate

  subroutine sequencePrepend(self,operator_)
    !!{
    Prepend an operator to the sequence.
    !!}
    implicit none
    class(outputAnalysisPropertyOperatorSequence), intent(inout)          :: self
    class(outputAnalysisPropertyOperatorClass   ), intent(in   ), target  :: operator_
    type (propertyOperatorList                  )               , pointer :: operatorNew

    allocate(operatorNew)
    operatorNew%operator_ => operator_
    operatorNew%next      => self       %operators
    self       %operators => operatorNew
    return
  end subroutine sequencePrepend

  subroutine sequenceApppend(self,operator_)
    !!{
    Append an operator to the sequence.
    !!}
    implicit none
    class(outputAnalysisPropertyOperatorSequence), intent(inout)          :: self
    class(outputAnalysisPropertyOperatorClass   ), intent(in   ), target  :: operator_
    type (propertyOperatorList                  )               , pointer :: operatorNew, operatorCurrent

    allocate(operatorNew)
    operatorNew    %operator_ => operator_
    operatorNew    %next      => null()
    operatorCurrent           => self     %operators
    do while (associated(operatorCurrent%next))
       operatorCurrent => operatorCurrent%next
    end do
    operatorCurrent%next => operatorNew
    return
  end subroutine sequenceApppend
