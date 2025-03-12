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
Implements a weight operator class in which the weight is multiplied by an integral over the N-body
halo mass distribution.
!!}

  use :: Statistics_NBody_Halo_Mass_Errors, only : nbodyHaloMassErrorClass

  !![
  <outputAnalysisWeightOperator name="outputAnalysisWeightOperatorNbodyMass">
   <description>A weight operator class in which the weight is multiplied by an integral over the N-body halo mass distribution.</description>
  </outputAnalysisWeightOperator>
  !!]
  type, extends(outputAnalysisWeightOperatorNormal) :: outputAnalysisWeightOperatorNbodyMass
     !!{
     A weight operator class in which the weight is multiplied by an integral over the N-body halo mass distribution.
     !!}
     private
     class(nbodyHaloMassErrorClass), pointer :: nbodyHaloMassError_ => null()
   contains
     final     ::                 nbodyMassDestructor
     procedure :: rootVariance => nbodyMassRootVariance
  end type outputAnalysisWeightOperatorNbodyMass

  interface outputAnalysisWeightOperatorNbodyMass
     !!{
     Constructors for the ``nbodyMass'' output analysis class.
     !!}
     module procedure nbodyMassConstructorParameters
     module procedure nbodyMassConstructorInternal
  end interface outputAnalysisWeightOperatorNbodyMass

contains

  function nbodyMassConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``nbodyMass'' output analysis weight operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (outputAnalysisWeightOperatorNbodyMass)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    class           (nodePropertyExtractorClass           ), pointer       :: nodePropertyExtractor_
    class           (outputAnalysisPropertyOperatorClass  ), pointer       :: outputAnalysisPropertyOperator_
    class           (nbodyHaloMassErrorClass              ), pointer       :: nbodyHaloMassError_
    double precision                                                       :: rangeLower                      , rangeUpper

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>rangeLower</name>
      <source>parameters</source>
      <description>Lower integration limit for the nbodyMass distribution weight operator.</description>
    </inputParameter>
    <inputParameter>
      <name>rangeUpper</name>
      <source>parameters</source>
      <description>Upper integration limit for the nbodyMass distribution weight operator.</description>
    </inputParameter>
    <objectBuilder class="nodePropertyExtractor"           name="nodePropertyExtractor_"           source="parameters"/>
    <objectBuilder class="outputAnalysisPropertyOperator"  name="outputAnalysisPropertyOperator_"  source="parameters"/>
    <objectBuilder class="nbodyHaloMassError"              name="nbodyHaloMassError_"              source="parameters"/>
    !!]
    self=outputAnalysisWeightOperatorNbodyMass(rangeLower,rangeUpper,nodePropertyExtractor_,outputAnalysisPropertyOperator_,nbodyHaloMassError_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="nodePropertyExtractor_"/>
    <objectDestructor name="outputAnalysisPropertyOperator_" />
    <objectDestructor name="nbodyHaloMassError_"             />
    !!]
    return
  end function nbodyMassConstructorParameters

  function nbodyMassConstructorInternal(rangeLower,rangeUpper,nodePropertyExtractor_,outputAnalysisPropertyOperator_,nbodyHaloMassError_) result (self)
    !!{
    Internal constructor for the ``nbodyMass'' output analysis distribution operator class.
    !!}
    use :: Error                   , only : Error_Report
    use :: Node_Property_Extractors, only : nodePropertyExtractorClass, nodePropertyExtractorScalar
    implicit none
    type            (outputAnalysisWeightOperatorNbodyMass)                        :: self
    double precision                                       , intent(in   )         :: rangeLower                      , rangeUpper
    class           (nodePropertyExtractorClass           ), intent(in   ), target :: nodePropertyExtractor_
    class           (outputAnalysisPropertyOperatorClass  ), intent(in   ), target :: outputAnalysisPropertyOperator_
    class(nbodyHaloMassErrorClass                         ), intent(in   ), target :: nbodyHaloMassError_
    !![
    <constructorAssign variables="rangeLower, rangeUpper, *nodePropertyExtractor_, *outputAnalysisPropertyOperator_, *nbodyHaloMassError_"/>
    !!]

    select type (nodePropertyExtractor_)
    class is (nodePropertyExtractorScalar)
       ! This is acceptable.
    class default
       call Error_Report('property extrator must be of scalar class'//{introspection:location})
    end select
    return
  end function nbodyMassConstructorInternal

  subroutine nbodyMassDestructor(self)
    !!{
    Destructor for  the ``nbodyMass'' output analysis weight operator class.
    !!}
    type(outputAnalysisWeightOperatorNbodyMass), intent(inout) :: self

    !![
    <objectDestructor name="self%nbodyHaloMassError_"             />
    <objectDestructor name="self%nodePropertyExtractor_"/>
    <objectDestructor name="self%outputAnalysisPropertyOperator_" />
    !!]
    return
  end subroutine nbodyMassDestructor

  double precision function nbodyMassRootVariance(self,node,propertyValue,propertyValueIntrinsic,propertyType,propertyQuantity,outputIndex)
    !!{
    Return the root variance for use in the ``nbodyMass'' output analysis weight operator class.
    !!}
    use :: Node_Property_Extractors, only : nodePropertyExtractorScalar
    implicit none
    class           (outputAnalysisWeightOperatorNbodyMass        ), intent(inout) :: self
    type            (treeNode                                     ), intent(inout) :: node
    double precision                                               , intent(in   ) :: propertyValue         , propertyValueIntrinsic
    type            (enumerationOutputAnalysisPropertyTypeType    ), intent(in   ) :: propertyType
    type            (enumerationOutputAnalysisPropertyQuantityType), intent(in   ) :: propertyQuantity
    integer         (c_size_t                                     ), intent(in   ) :: outputIndex
    double precision                                                               :: nbodyMassPropertyValue
    type            (enumerationOutputAnalysisPropertyTypeType)                    :: nbodyMassPropertyType
    !$GLC attributes unused :: propertyValue, propertyValueIntrinsic, propertyType, propertyQuantity, outputIndex

    ! Extract property and operate on it.
    nbodyMassPropertyType   = self%nodePropertyExtractor_%type           (    )
    select type (extractor_ => self%nodePropertyExtractor_)
    class is (nodePropertyExtractorScalar)
       nbodyMassPropertyValue=+                extractor_%extract        (node)
    class default
       nbodyMassPropertyValue=+0.0d0
    end select
    nbodyMassRootVariance   =+self%nbodyHaloMassError_   %errorFractional(node) &
         &                   *nbodyMassPropertyValue
    return
  end function nbodyMassRootVariance
