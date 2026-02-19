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
Implements a subsampling weight operator class.
!!}

  !![
  <outputAnalysisWeightOperator name="outputAnalysisWeightOperatorSubsampling">
   <description>A subsampling weight operator class.</description>
  </outputAnalysisWeightOperator>
  !!]
  type, extends(outputAnalysisWeightOperatorClass) :: outputAnalysisWeightOperatorSubsampling
     !!{
     A subsampling weight operator class.
     !!}
     private
   contains
     procedure :: operate => subsamplingOperate
  end type outputAnalysisWeightOperatorSubsampling

  interface outputAnalysisWeightOperatorSubsampling
     !!{
     Constructors for the \refClass{outputAnalysisWeightOperatorSubsampling} output analysis class.
     !!}
     module procedure subsamplingConstructorParameters
  end interface outputAnalysisWeightOperatorSubsampling

contains

  function subsamplingConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisWeightOperatorSubsampling} output analysis weight operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(outputAnalysisWeightOperatorSubsampling)                :: self
    type(inputParameters                        ), intent(inout) :: parameters

    self=outputAnalysisWeightOperatorSubsampling()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function subsamplingConstructorParameters

  double precision function subsamplingOperate(self,weightValue,node,propertyValue,propertyValueIntrinsic,propertyType,propertyQuantity,outputIndex)
    !!{
    Implement an subsampling output analysis weight operator.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (outputAnalysisWeightOperatorSubsampling      ), intent(inout) :: self
    type            (treeNode                                     ), intent(inout) :: node
    double precision                                               , intent(in   ) :: propertyValue   , propertyValueIntrinsic, &
         &                                                                            weightValue
    type            (enumerationOutputAnalysisPropertyTypeType    ), intent(in   ) :: propertyType
    type            (enumerationOutputAnalysisPropertyQuantityType), intent(in   ) :: propertyQuantity
    integer         (c_size_t                                     ), intent(in   ) :: outputIndex
    !$GLC attributes unused :: self, propertyValue, propertyValueIntrinsic, propertyType, propertyQuantity, outputIndex

    subsamplingOperate=weightValue*node%subsamplingWeight()
    return
  end function subsamplingOperate
