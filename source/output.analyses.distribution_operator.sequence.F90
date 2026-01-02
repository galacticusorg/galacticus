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
Implements a sequence output analysis distribution operator class.
!!}

  type, public :: distributionOperatorList
     class(outputAnalysisDistributionOperatorClass), pointer :: operator_ => null()
     type (distributionOperatorList               ), pointer :: next      => null()
  end type distributionOperatorList

  !![
  <outputAnalysisDistributionOperator name="outputAnalysisDistributionOperatorSequence">
   <description>A sequence output analysis distribution operator class.</description>
   <linkedList type="distributionOperatorList" variable="operators" next="next" object="operator_" objectType="outputAnalysisDistributionOperatorClass"/>
  </outputAnalysisDistributionOperator>
  !!]
  type, extends(outputAnalysisDistributionOperatorClass) :: outputAnalysisDistributionOperatorSequence
     !!{
     A sequence output distribution operator class.
     !!}
     private
     type(distributionOperatorList), pointer :: operators => null()
   contains
     !![
     <methods>
       <method description="Prepend an operator to a sequence of distribution operators." method="prepend" />
     </methods>
     !!]
     final     ::                        sequenceDestructor
     procedure :: operateScalar       => sequenceOperateScalar
     procedure :: operateDistribution => sequenceOperateDistribution
     procedure :: prepend             => sequencePrepend
  end type outputAnalysisDistributionOperatorSequence

  interface outputAnalysisDistributionOperatorSequence
     !!{
     Constructors for the \refClass{outputAnalysisDistributionOperatorSequence} output analysis class.
     !!}
     module procedure sequenceConstructorParameters
     module procedure sequenceConstructorInternal
  end interface outputAnalysisDistributionOperatorSequence

contains

  function sequenceConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisDistributionOperatorSequence} output analysis distribution operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (outputAnalysisDistributionOperatorSequence)                :: self
    type   (inputParameters                           ), intent(inout) :: parameters
    type   (distributionOperatorList                  ), pointer       :: operator_
    integer                                                            :: i

    self     %operators => null()
    operator_           => null()
    do i=1,parameters%copiesCount('outputAnalysisDistributionOperator',zeroIfNotPresent=.true.)
       if (associated(operator_)) then
          allocate(operator_%next)
          operator_ => operator_%next
       else
          allocate(self%operators)
          operator_ => self%operators
       end if
       !![
       <objectBuilder class="outputAnalysisDistributionOperator" name="operator_%operator_" source="parameters" copy="i" />
       !!]
    end do
    !![
    <inputParametersValidate source="parameters" multiParameters="outputAnalysisDistributionOperator"/>
    !!]
    return
  end function sequenceConstructorParameters

  function sequenceConstructorInternal(operators) result (self)
    !!{
    Internal constructor for the sequence merger tree normalizer class.
    !!}
    implicit none
    type(outputAnalysisDistributionOperatorSequence)                        :: self
    type(distributionOperatorList                  ), target, intent(in   ) :: operators
    type(distributionOperatorList                  ), pointer               :: operator_

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
    Destructor for the sequence distribution operator class.
    !!}
    implicit none
    type(outputAnalysisDistributionOperatorSequence), intent(inout) :: self
    type(distributionOperatorList                  ), pointer       :: operator_, operatorNext

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

  function sequenceOperateScalar(self,propertyValue,propertyType,propertyValueMinimum,propertyValueMaximum,outputIndex,node)
    !!{
    Implement an sequence output analysis distribution operator.
    !!}
    implicit none
    class           (outputAnalysisDistributionOperatorSequence), intent(inout)                                        :: self
    double precision                                            , intent(in   )                                        :: propertyValue
    type            (enumerationOutputAnalysisPropertyTypeType ), intent(in   )                                        :: propertyType
    double precision                                            , intent(in   ), dimension(:)                          :: propertyValueMinimum , propertyValueMaximum
    integer         (c_size_t                                  ), intent(in   )                                        :: outputIndex
    type            (treeNode                                  ), intent(inout)                                        :: node
    double precision                                                           , dimension(size(propertyValueMinimum)) :: sequenceOperateScalar
    type            (distributionOperatorList                  ), pointer                                              :: operator_

    operator_  => self%operators
    do while (associated(operator_))
       ! For first operator, apply to a scalar. Subsequent operators are applied to the distribution.
       if (associated(operator_,self%operators)) then
          sequenceOperateScalar=operator_%operator_%operateScalar      (propertyValue        ,propertyType,propertyValueMinimum,propertyValueMaximum,outputIndex,node)
       else
          sequenceOperateScalar=operator_%operator_%operateDistribution(sequenceOperateScalar,propertyType,propertyValueMinimum,propertyValueMaximum,outputIndex,node)
       end if
       operator_ => operator_%next
    end do
    return
  end function sequenceOperateScalar

  function sequenceOperateDistribution(self,distribution,propertyType,propertyValueMinimum,propertyValueMaximum,outputIndex,node)
    !!{
    Implement a random error output analysis distribution operator.
    !!}
    implicit none
    class           (outputAnalysisDistributionOperatorSequence), intent(inout)                                        :: self
    double precision                                            , intent(in   ), dimension(:)                          :: distribution
    type            (enumerationOutputAnalysisPropertyTypeType ), intent(in   )                                        :: propertyType
    double precision                                            , intent(in   ), dimension(:)                          :: propertyValueMinimum       , propertyValueMaximum
    integer         (c_size_t                                  ), intent(in   )                                        :: outputIndex
    type            (treeNode                                  ), intent(inout)                                        :: node
    double precision                                                           , dimension(size(propertyValueMinimum)) :: sequenceOperateDistribution
    type            (distributionOperatorList                  ), pointer                                              :: operator_

    operator_                   => self        %operators
    sequenceOperateDistribution =  distribution
    do while (associated(operator_))
       sequenceOperateDistribution=operator_%operator_%operateDistribution(sequenceOperateDistribution,propertyType,propertyValueMinimum,propertyValueMaximum,outputIndex,node)
       operator_ => operator_%next
    end do
   return
  end function sequenceOperateDistribution

  subroutine sequencePrepend(self,operator_)
    !!{
    Prepend an operator to the sequence.
    !!}
    implicit none
    class(outputAnalysisDistributionOperatorSequence), intent(inout)          :: self
    class(outputAnalysisDistributionOperatorClass   ), intent(in   ), target  :: operator_
    type (distributionOperatorList                  )               , pointer :: operatorNew

    allocate(operatorNew)
    operatorNew%operator_ => operator_
    operatorNew%next      => self       %operators
    self       %operators => operatorNew
    return
  end subroutine sequencePrepend
