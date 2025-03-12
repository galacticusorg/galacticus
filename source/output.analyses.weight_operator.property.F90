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
Implements an analysis weight operator class which weights by a property value.
!!}
  use :: Node_Property_Extractors          , only : nodePropertyExtractorClass
  use :: Output_Analysis_Property_Operators, only : outputAnalysisPropertyOperatorClass

  !![
  <outputAnalysisWeightOperator name="outputAnalysisWeightOperatorProperty">
   <description>An analysis weight operator class which weights by a property value.</description>
  </outputAnalysisWeightOperator>
  !!]
  type, extends(outputAnalysisWeightOperatorClass) :: outputAnalysisWeightOperatorProperty
     !!{
     An weight operator class which weights by a property value.
     !!}
     private
     class(nodePropertyExtractorClass         ), pointer :: extractor_ => null()
     class(outputAnalysisPropertyOperatorClass), pointer :: operator_  => null()
   contains
     final     ::            propertyDestructor
     procedure :: operate => propertyOperate
  end type outputAnalysisWeightOperatorProperty

  interface outputAnalysisWeightOperatorProperty
     !!{
     Constructors for the ``property'' output analysis class.
     !!}
     module procedure propertyConstructorParameters
     module procedure propertyConstructorInternal
  end interface outputAnalysisWeightOperatorProperty

contains

  function propertyConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``property'' output analysis weight operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (outputAnalysisWeightOperatorProperty)                :: self
    type (inputParameters                     ), intent(inout) :: parameters
    class(nodePropertyExtractorClass          ), pointer       :: extractor_
    class(outputAnalysisPropertyOperatorClass ), pointer       :: operator_

    ! Check and read parameters.
    !![
    <objectBuilder class="nodePropertyExtractor"          name="extractor_" source="parameters"/>
    <objectBuilder class="outputAnalysisPropertyOperator" name="operator_"  source="parameters"/>
    !!]
    ! Construct the object.
    self=outputAnalysisWeightOperatorProperty(extractor_,operator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="extractor_"/>
    <objectDestructor name="operator_" />
    !!]
    return
  end function propertyConstructorParameters

  function propertyConstructorInternal(extractor_,operator_) result(self)
    !!{
    Internal constructor for the ``property'' output analysis weight operator class.
    !!}
    use :: Error                   , only : Error_Report
    use :: Node_Property_Extractors, only : nodePropertyExtractorClass, nodePropertyExtractorScalar
    implicit none
    type (outputAnalysisWeightOperatorProperty)                        :: self
    class(nodePropertyExtractorClass          ), intent(in   ), target :: extractor_
    class(outputAnalysisPropertyOperatorClass ), intent(in   ), target :: operator_
    !![
    <constructorAssign variables="*extractor_, *operator_"/>
    !!]

    select type (extractor_)
    class is (nodePropertyExtractorScalar)
       ! This is acceptable.
    class default
       call Error_Report('property extrator must be of scalar class'//{introspection:location})
    end select
    return
  end function propertyConstructorInternal

  subroutine propertyDestructor(self)
    !!{
    Destructor for the ``property'' output analysis weight operator class.
    !!}
    implicit none
    type(outputAnalysisWeightOperatorProperty), intent(inout) :: self

    !![
    <objectDestructor name="self%extractor_"/>
    <objectDestructor name="self%operator_" />
    !!]
    return
  end subroutine propertyDestructor

  double precision function propertyOperate(self,weightValue,node,propertyValue,propertyValueIntrinsic,propertyType,propertyQuantity,outputIndex)
    !!{
    Implement an property output analysis weight operator.
    !!}
    use, intrinsic :: ISO_C_Binding           , only : c_size_t
    use            :: Node_Property_Extractors, only : nodePropertyExtractorScalar
    implicit none
    class           (outputAnalysisWeightOperatorProperty         ), intent(inout) :: self
    type            (treeNode                                     ), intent(inout) :: node
    double precision                                               , intent(in   ) :: propertyValue      , propertyValueIntrinsic, &
         &                                                                            weightValue
    type            (enumerationOutputAnalysisPropertyTypeType    ), intent(in   ) :: propertyType
    type            (enumerationOutputAnalysisPropertyQuantityType), intent(in   ) :: propertyQuantity
    integer         (c_size_t                                     ), intent(in   ) :: outputIndex
    double precision                                                               :: weightPropertyValue
    type            (enumerationOutputAnalysisPropertyTypeType    )                :: weightPropertyType
    !$GLC attributes unused :: propertyType, propertyValueIntrinsic, propertyValue, propertyQuantity

    select type (extractor_ => self%extractor_)
    class is (nodePropertyExtractorScalar)
       weightPropertyValue=+            extractor_%extract(node                                                   )
    class default
       weightPropertyValue=+0.0d0
    end select
    weightPropertyType    = self       %extractor_%type   (                                                       )
    propertyOperate       =+self       %operator_ %operate(weightPropertyValue,node,weightPropertyType,outputIndex) &
         &                 *weightValue
    return
  end function propertyOperate
