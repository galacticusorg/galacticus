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
Contains a module which implements a high-pass filter analysis weight operator class.
!!}

  !![
  <outputAnalysisWeightOperator name="outputAnalysisWeightOperatorFilterHighPass">
   <description>A high-pass filter analysis weight operator class.</description>
  </outputAnalysisWeightOperator>
  !!]
  type, extends(outputAnalysisWeightOperatorClass) :: outputAnalysisWeightOperatorFilterHighPass
     !!{
     A high-pass filter weight operator class.
     !!}
     private
     double precision :: filterThreshold
   contains
     procedure :: operate => filterHighPassOperate
  end type outputAnalysisWeightOperatorFilterHighPass

  interface outputAnalysisWeightOperatorFilterHighPass
     !!{
     Constructors for the ``filterHighPass'' output analysis class.
     !!}
     module procedure filterHighPassConstructorParameters
     module procedure filterHighPassConstructorInternal
  end interface outputAnalysisWeightOperatorFilterHighPass

contains

  function filterHighPassConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``filterHighPass'' output analysis weight operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(outputAnalysisWeightOperatorFilterHighPass)                :: self
    type(inputParameters                           ), intent(inout) :: parameters
    double precision                                                :: filterThreshold

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>filterThreshold</name>
      <source>parameters</source>
      <variable>filterThreshold</variable>
      <description>Threshold for the high-pass filter distribution operator.</description>
    </inputParameter>
    !!]
    self=outputAnalysisWeightOperatorFilterHighPass(filterThreshold)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function filterHighPassConstructorParameters

  function filterHighPassConstructorInternal(filterThreshold) result (self)
    !!{
    Internal constructor for the ``filterHighPass'' output analysis distribution operator class.
    !!}
    implicit none
    type            (outputAnalysisWeightOperatorFilterHighPass)                :: self
    double precision                                            , intent(in   ) :: filterThreshold
    !![
    <constructorAssign variables="filterThreshold"/>
    !!]

    return
  end function filterHighPassConstructorInternal

  double precision function filterHighPassOperate(self,weightValue,node,propertyValue,propertyValueIntrinsic,propertyType,propertyQuantity,outputIndex)
    !!{
    Implement an filterHighPass output analysis weight operator.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (outputAnalysisWeightOperatorFilterHighPass   ), intent(inout) :: self
    type            (treeNode                                     ), intent(inout) :: node
    double precision                                               , intent(in   ) :: propertyValue   , propertyValueIntrinsic, &
         &                                                                            weightValue
    type            (enumerationOutputAnalysisPropertyTypeType    ), intent(in   ) :: propertyType
    type            (enumerationOutputAnalysisPropertyQuantityType), intent(in   ) :: propertyQuantity
    integer         (c_size_t                                     ), intent(in   ) :: outputIndex
    !$GLC attributes unused :: node, propertyValueIntrinsic, propertyType, outputIndex, propertyQuantity

    if (propertyValue > self%filterThreshold) then
       filterHighPassOperate=weightValue
    else
       filterHighPassOperate=0.0d0
    end if
    return
  end function filterHighPassOperate
