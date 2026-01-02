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
Implements a sequence output analysis weight operator class.
!!}

  type, public :: weightOperatorList
     class(outputAnalysisWeightOperatorClass), pointer :: operator_ => null()
     type (weightOperatorList               ), pointer :: next      => null()
  end type weightOperatorList

  !![
  <outputAnalysisWeightOperator name="outputAnalysisWeightOperatorSequence">
   <description>A sequence output analysis weight operator class.</description>
   <linkedList type="weightOperatorList" variable="operators" next="next" object="operator_" objectType="outputAnalysisWeightOperatorClass"/>
  </outputAnalysisWeightOperator>
  !!]
  type, extends(outputAnalysisWeightOperatorClass) :: outputAnalysisWeightOperatorSequence
     !!{
     A sequence output weight operator class.
     !!}
     private
     type(weightOperatorList), pointer :: operators => null()
   contains
     !![
     <methods>
       <method description="Prepend an operator to a sequence of weight operators." method="prepend" />
     </methods>
     !!]
     final     ::            sequenceDestructor
     procedure :: operate => sequenceOperate
     procedure :: prepend => sequencePrepend
  end type outputAnalysisWeightOperatorSequence

  interface outputAnalysisWeightOperatorSequence
     !!{
     Constructors for the \refClass{outputAnalysisWeightOperatorSequence} output analysis class.
     !!}
     module procedure sequenceConstructorParameters
     module procedure sequenceConstructorInternal
  end interface outputAnalysisWeightOperatorSequence

contains

  function sequenceConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisWeightOperatorSequence} output analysis weight operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (outputAnalysisWeightOperatorSequence)                :: self
    type   (inputParameters                     ), intent(inout) :: parameters
    type   (weightOperatorList                  ), pointer       :: operator_
    integer                                                      :: i

    self     %operators => null()
    operator_           => null()
    do i=1,parameters%copiesCount('outputAnalysisWeightOperator',zeroIfNotPresent=.true.)
       if (associated(operator_)) then
          allocate(operator_%next)
          operator_ => operator_%next
       else
          allocate(self%operators)
          operator_ => self%operators
       end if
       !![
       <objectBuilder class="outputAnalysisWeightOperator" name="operator_%operator_" source="parameters" copy="i" />
       !!]
    end do
    !![
    <inputParametersValidate source="parameters" multiParameters="outputAnalysisWeightOperator"/>
    !!]
    return
  end function sequenceConstructorParameters

  function sequenceConstructorInternal(operators) result (self)
    !!{
    Internal constructor for the sequence output analysis weight operator class.
    !!}
    implicit none
    type(outputAnalysisWeightOperatorSequence)                        :: self
    type(weightOperatorList                  ), target, intent(in   ) :: operators
    type(weightOperatorList                  ), pointer               :: operator_

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
    Destructor for the sequence weight operator class.
    !!}
    implicit none
    type(outputAnalysisWeightOperatorSequence), intent(inout) :: self
    type(weightOperatorList                  ), pointer       :: operator_, operatorNext

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

  double precision function sequenceOperate(self,weightValue,node,propertyValue,propertyValueIntrinsic,propertyType,propertyQuantity,outputIndex)
    !!{
    Implement an sequence output analysis weight operator.
    !!}
    implicit none
    class           (outputAnalysisWeightOperatorSequence         ), intent(inout) :: self
    type            (treeNode                                     ), intent(inout) :: node
    double precision                                               , intent(in   ) :: weightValue     , propertyValueIntrinsic, &
         &                                                                            propertyValue
    type            (enumerationOutputAnalysisPropertyTypeType    ), intent(in   ) :: propertyType
    type            (enumerationOutputAnalysisPropertyQuantityType), intent(in   ) :: propertyQuantity
    integer         (c_size_t                                     ), intent(in   ) :: outputIndex
    type            (weightOperatorList                           ), pointer       :: operator_

    sequenceOperate =  weightValue
    operator_       => self%operators
    do while (associated(operator_))
       sequenceOperate =  operator_%operator_%operate(sequenceOperate,node,propertyValue,propertyValueIntrinsic,propertyType,propertyQuantity,outputIndex)
       operator_       => operator_%next
    end do
    return
  end function sequenceOperate

  subroutine sequencePrepend(self,operator_)
    !!{
    Prepend an operator to the sequence.
    !!}
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
