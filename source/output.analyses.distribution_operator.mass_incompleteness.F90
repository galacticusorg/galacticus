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
  Implements the effects of incompleteness as a function of mass on the distribution.
  !!}

  use :: Mass_Function_Incompletenesses, only : massFunctionIncompletenessClass

  !![
  <outputAnalysisDistributionOperator name="outputAnalysisDistributionOperatorMassIncompleteness">
   <description>An output analysis distribution operator class which implements the effects of incompleteness as a function of mass on the distribution.</description>
  </outputAnalysisDistributionOperator>
  !!]
  type, extends(outputAnalysisDistributionOperatorClass) :: outputAnalysisDistributionOperatorMassIncompleteness
     !!{
     An output distribution operator class which implements the effects of incompleteness as a function of mass on the distribution.
     !!}
     private
     class (massFunctionIncompletenessClass), pointer :: massFunctionIncompleteness_ => null()
   contains
     final     ::                        massIncompletenessDestructor
     procedure :: operateScalar       => massIncompletenessOperateScalar
     procedure :: operateDistribution => massIncompletenessOperateDistribution
  end type outputAnalysisDistributionOperatorMassIncompleteness

  interface outputAnalysisDistributionOperatorMassIncompleteness
     !!{
     Constructors for the \refClass{outputAnalysisDistributionOperatorMassIncompleteness} output distribution operator class.
     !!}
     module procedure massIncompletenessConstructorParameters
     module procedure massIncompletenessConstructorInternal
  end interface outputAnalysisDistributionOperatorMassIncompleteness

contains

  function massIncompletenessConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisDistributionOperatorMassIncompleteness} output analysis distribution operator operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (outputAnalysisDistributionOperatorMassIncompleteness)                :: self
    type (inputParameters                                     ), intent(inout) :: parameters
    class(massFunctionIncompletenessClass                     ), pointer       :: massFunctionIncompleteness_

    !![
    <objectBuilder class="massFunctionIncompleteness" name="massFunctionIncompleteness_" source="parameters"/>
    !!]
    self=outputAnalysisDistributionOperatorMassIncompleteness(massFunctionIncompleteness_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="massFunctionIncompleteness_"/>
    !!]
    return
  end function massIncompletenessConstructorParameters

  function massIncompletenessConstructorInternal(massFunctionIncompleteness_) result(self)
    !!{
    Internal constructor for the \refClass{outputAnalysisDistributionOperatorMassIncompleteness} output analysis distribution operator class.
    !!}
    implicit none
    type (outputAnalysisDistributionOperatorMassIncompleteness)                        :: self
    class(massFunctionIncompletenessClass                     ), intent(in   ), target :: massFunctionIncompleteness_
    !![
    <constructorAssign variables="*massFunctionIncompleteness_"/>
    !!]

    return
  end function massIncompletenessConstructorInternal

  subroutine massIncompletenessDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisDistributionOperatorMassIncompleteness} output analysis distribution operator operator class.
    !!}
    implicit none
    type(outputAnalysisDistributionOperatorMassIncompleteness), intent(inout) :: self

    !![
    <objectDestructor name="self%massFunctionIncompleteness_"/>
    !!]
    return
  end subroutine massIncompletenessDestructor

  function massIncompletenessOperateScalar(self,propertyValue,propertyType,propertyValueMinimum,propertyValueMaximum,outputIndex,node)
    !!{
    Implement a mass incompleteness output analysis distribution operator.
    !!}
    use :: Error                  , only : Error_Report
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear, outputAnalysisPropertyTypeLog10
    implicit none
    class           (outputAnalysisDistributionOperatorMassIncompleteness), intent(inout)                                        :: self
    double precision                                                      , intent(in   )                                        :: propertyValue
    type            (enumerationOutputAnalysisPropertyTypeType           ), intent(in   )                                        :: propertyType
    double precision                                                      , intent(in   ), dimension(:)                          :: propertyValueMinimum            , propertyValueMaximum
    integer         (c_size_t                                            ), intent(in   )                                        :: outputIndex
    type            (treeNode                                            ), intent(inout)                                        :: node
    double precision                                                                     , dimension(size(propertyValueMinimum)) :: massIncompletenessOperateScalar
    integer                                                                                                                      :: i
    double precision                                                                                                             :: mass
    !$GLC attributes unused :: outputIndex, node

    do i=1,size(propertyValueMinimum)
       if (propertyValue >= propertyValueMinimum(i) .and. propertyValue < propertyValueMaximum(i)) then
          select case (propertyType%ID)
          case (outputAnalysisPropertyTypeLinear%ID)
             mass=        propertyValue
          case (outputAnalysisPropertyTypeLog10 %ID)
             mass=10.0d0**propertyValue
          case default
             mass= 0.0d0
             call Error_Report('unsupported property type'//{introspection:location})
          end select
          massIncompletenessOperateScalar(i)=self%massFunctionIncompleteness_%completeness(mass)
       end if
    end do
    return
  end function massIncompletenessOperateScalar

  function massIncompletenessOperateDistribution(self,distribution,propertyType,propertyValueMinimum,propertyValueMaximum,outputIndex,node)
    !!{
    Implement a mass incompleteness output analysis distribution operator.
    !!}
    use :: Error                  , only : Error_Report
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear, outputAnalysisPropertyTypeLog10
    implicit none
    class           (outputAnalysisDistributionOperatorMassIncompleteness), intent(inout)                                        :: self
    double precision                                                      , intent(in   ), dimension(:)                          :: distribution
    type            (enumerationOutputAnalysisPropertyTypeType           ), intent(in   )                                        :: propertyType
    double precision                                                      , intent(in   ), dimension(:)                          :: propertyValueMinimum                  , propertyValueMaximum
    integer         (c_size_t                                            ), intent(in   )                                        :: outputIndex
    type            (treeNode                                            ), intent(inout)                                        :: node
    double precision                                                                     , dimension(size(propertyValueMinimum)) :: massIncompletenessOperateDistribution
    integer                                                                                                                      :: i
    double precision                                                                                                             :: mass
    !$GLC attributes unused :: outputIndex, node

    do i=1,size(propertyValueMinimum)
       select case (propertyType%ID)
       case (outputAnalysisPropertyTypeLinear%ID)
          mass=         0.5d0*(propertyValueMinimum(i)+propertyValueMaximum(i))
       case (outputAnalysisPropertyTypeLog10%ID)
          mass=10.0d0**(0.5d0*(propertyValueMinimum(i)+propertyValueMaximum(i)))
       case default
          mass= 0.0d0
          call Error_Report('unsupported property type'//{introspection:location})
       end select
       massIncompletenessOperateDistribution(i)=+                                 distribution(i   ) &
            &                                   *self%massFunctionIncompleteness_%completeness(mass)
    end do
    return
  end function massIncompletenessOperateDistribution
