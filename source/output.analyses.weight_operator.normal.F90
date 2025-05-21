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
Implements a weight operator class in which the weight is multiplied by an integral over a normal distribution.
!!}

  use :: Node_Property_Extractors          , only : nodePropertyExtractorClass
  use :: Output_Analysis_Property_Operators, only : outputAnalysisPropertyOperatorClass

  !![
  <outputAnalysisWeightOperator name="outputAnalysisWeightOperatorNormal">
   <description>A weight operator class in which the weight is multiplied by an integral over a normal distribution.</description>
  </outputAnalysisWeightOperator>
  !!]
  type, extends(outputAnalysisWeightOperatorClass) :: outputAnalysisWeightOperatorNormal
     !!{
     A high-pass filter weight operator class.
     !!}
     private
     class           (nodePropertyExtractorClass          ), pointer :: nodePropertyExtractor_          => null()
     class           (outputAnalysisPropertyOperatorClass ), pointer :: outputAnalysisPropertyOperator_ => null()
     double precision                                                :: rangeLower                               , rangeUpper, &
          &                                                             rootVariance_
   contains
     !![
     <methods>
       <method description="Return the root-variance to use in the weight operator." method="rootVariance" />
     </methods>
     !!]
     final     ::                 normalDestructor
     procedure :: operate      => normalOperate
     procedure :: rootVariance => normalRootVariance
  end type outputAnalysisWeightOperatorNormal

  interface outputAnalysisWeightOperatorNormal
     !!{
     Constructors for the \refClass{outputAnalysisWeightOperatorNormal} output analysis class.
     !!}
     module procedure normalConstructorParameters
     module procedure normalConstructorInternal
  end interface outputAnalysisWeightOperatorNormal

contains

  function normalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisWeightOperatorNormal} output analysis weight operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (outputAnalysisWeightOperatorNormal  )                :: self
    type            (inputParameters                     ), intent(inout) :: parameters
    class           (nodePropertyExtractorClass          ), pointer       :: nodePropertyExtractor_
    class           (outputAnalysisPropertyOperatorClass ), pointer       :: outputAnalysisPropertyOperator_
    double precision                                                      :: rangeLower                      , rangeUpper, &
         &                                                                   rootVariance_

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>rangeLower</name>
      <source>parameters</source>
      <description>Lower integration limit for the normal distribution weight operator.</description>
    </inputParameter>
    <inputParameter>
      <name>rangeUpper</name>
      <source>parameters</source>
      <description>Upper integration limit for the normal distribution weight operator.</description>
    </inputParameter>
    <inputParameter>
      <name>rootVariance</name>
      <variable>rootVariance_</variable>
      <source>parameters</source>
      <description>Root variance for the normal distribution weight operator.</description>
    </inputParameter>
    <objectBuilder class="nodePropertyExtractor"          name="nodePropertyExtractor_"          source="parameters"/>
    <objectBuilder class="outputAnalysisPropertyOperator" name="outputAnalysisPropertyOperator_" source="parameters"/>
    !!]
    self=outputAnalysisWeightOperatorNormal(rangeLower,rangeUpper,rootVariance_,nodePropertyExtractor_,outputAnalysisPropertyOperator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="nodePropertyExtractor_"         />
    <objectDestructor name="outputAnalysisPropertyOperator_"/>
    !!]
    return
  end function normalConstructorParameters

  function normalConstructorInternal(rangeLower,rangeUpper,rootVariance_,nodePropertyExtractor_,outputAnalysisPropertyOperator_) result (self)
    !!{
    Internal constructor for the \refClass{outputAnalysisWeightOperatorNormal} output analysis distribution operator class.
    !!}
    use :: Error                   , only : Error_Report
    use :: Node_Property_Extractors, only : nodePropertyExtractorClass, nodePropertyExtractorScalar
    implicit none
    type            (outputAnalysisWeightOperatorNormal )                        :: self
    double precision                                     , intent(in   )         :: rangeLower                      , rangeUpper, &
         &                                                                          rootVariance_
    class           (nodePropertyExtractorClass         ), intent(in   ), target :: nodePropertyExtractor_
    class           (outputAnalysisPropertyOperatorClass), intent(in   ), target :: outputAnalysisPropertyOperator_
    !![
    <constructorAssign variables="rangeLower, rangeUpper, rootVariance_, *nodePropertyExtractor_, *outputAnalysisPropertyOperator_"/>
    !!]

    select type (nodePropertyExtractor_)
    class is (nodePropertyExtractorScalar)
       ! This is acceptable.
    class default
       call Error_Report('property extrator must be of scalar class'//{introspection:location})
    end select
    return
  end function normalConstructorInternal

  subroutine normalDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisWeightOperatorNormal} output analysis weight operator class.
    !!}
    type(outputAnalysisWeightOperatorNormal), intent(inout) :: self

    !![
    <objectDestructor name="self%nodePropertyExtractor_"         />
    <objectDestructor name="self%outputAnalysisPropertyOperator_"/>
    !!]
    return
  end subroutine normalDestructor

  double precision function normalRootVariance(self,node,propertyValue,propertyValueIntrinsic,propertyType,propertyQuantity,outputIndex)
    !!{
    Return the root variance for use in the {\normalfont \ttfamily normal} output analysis weight operator class.
    !!}
    implicit none
    class           (outputAnalysisWeightOperatorNormal           ), intent(inout) :: self
    type            (treeNode                                     ), intent(inout) :: node
    double precision                                               , intent(in   ) :: propertyValue, propertyValueIntrinsic
    type            (enumerationOutputAnalysisPropertyTypeType    ), intent(in   ) :: propertyType
    type            (enumerationOutputAnalysisPropertyQuantityType), intent(in   ) :: propertyQuantity
    integer         (c_size_t                                     ), intent(in   ) :: outputIndex
    !$GLC attributes unused :: node, propertyValue, propertyValueIntrinsic, propertyType, propertyQuantity, outputIndex

    normalRootVariance=self%rootVariance_
    return
  end function normalRootVariance

  double precision function normalOperate(self,weightValue,node,propertyValue,propertyValueIntrinsic,propertyType,propertyQuantity,outputIndex)
    !!{
    Implement an normal output analysis weight operator.
    !!}
    use            :: Error_Functions         , only : Error_Function
    use, intrinsic :: ISO_C_Binding           , only : c_size_t
    use            :: Node_Property_Extractors, only : nodePropertyExtractorScalar
    implicit none
    class           (outputAnalysisWeightOperatorNormal           ), intent(inout) :: self
    type            (treeNode                                     ), intent(inout) :: node
    double precision                                               , intent(in   ) :: propertyValue      , propertyValueIntrinsic, &
         &                                                                            weightValue
    type            (enumerationOutputAnalysisPropertyTypeType    ), intent(in   ) :: propertyType
    type            (enumerationOutputAnalysisPropertyQuantityType), intent(in   ) :: propertyQuantity
    integer         (c_size_t                                     ), intent(in   ) :: outputIndex
    double precision                                                               :: normalPropertyValue, rootVariance
    type            (enumerationOutputAnalysisPropertyTypeType    )                :: normalPropertyType
    !$GLC attributes unused :: propertyValue, propertyValueIntrinsic, propertyType, propertyQuantity

    ! Extract property and operate on it.
    normalPropertyType    =self%nodePropertyExtractor_          %type   (                                                       )
    select type (extractor_ => self%nodePropertyExtractor_)
    class is (nodePropertyExtractorScalar)
       normalPropertyValue=                           extractor_%extract(                    node                               )
       class default
       normalPropertyValue=+0.0d0
    end select
    normalPropertyValue   =self%outputAnalysisPropertyOperator_ %operate(normalPropertyValue,node,normalPropertyType,outputIndex)
    ! Multiply weight by integral over a normal distribution between given limits, using the property value as the mean.
    if     (                                     &
         &   normalPropertyValue == +huge(0.0d0) &
         &  .or.                                 &
         &   normalPropertyValue == -huge(0.0d0) &
         & ) then
       normalOperate=+0.0d0
    else
       rootVariance =+self%rootVariance(node,propertyValue,propertyValueIntrinsic,propertyType,propertyQuantity,outputIndex)
       if (rootVariance <= 0.0d0) then
          ! If the root variance is zero we assume that the normal distribution has the limiting case of a δ-function, such
          ! that the weight is multiplied by 0 or 1 depending on whether the property values lies inside or outside of our
          ! integration range.
          normalOperate=0.0d0
          if     (                                        &
               &   normalPropertyValue >= self%rangeLower &
               &  .and.                                   &
               &   normalPropertyValue <  self%rangeUpper &
               & ) normalOperate=weightValue
       else
          normalOperate=+weightValue                                                                      &
               &        *(                                                                                &
               &          +Error_Function((self%rangeUpper-normalPropertyValue)/sqrt(2.0d0)/rootVariance) &
               &          -Error_Function((self%rangeLower-normalPropertyValue)/sqrt(2.0d0)/rootVariance) &
               &         )                                                                                &
               &        /2.0d0
       end if
    end if
    return
  end function normalOperate
