!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

!% Contains a module which implements a weight operator class in which the weight is multiplied by an integral over a normal distribution.

  use Output_Analysis_Property_Extractions
  use Output_Analysis_Property_Operators

  !# <outputAnalysisWeightOperator name="outputAnalysisWeightOperatorNormal" defaultThreadPrivate="yes">
  !#  <description>A weight operator class in which the weight is multiplied by an integral over a normal distribution.</description>
  !# </outputAnalysisWeightOperator>
  type, extends(outputAnalysisWeightOperatorClass) :: outputAnalysisWeightOperatorNormal
     !% A high-pass filter weight operator class.
     private
     class           (outputAnalysisPropertyExtractorClass), pointer :: outputAnalysisPropertyExtractor_
     class           (outputAnalysisPropertyOperatorClass ), pointer :: outputAnalysisPropertyOperator_
     double precision                                                :: rangeLower                      , rangeUpper, &
          &                                                             rootVariance
   contains
     final     ::            normalDestructor
     procedure :: operate => normalOperate
  end type outputAnalysisWeightOperatorNormal

  interface outputAnalysisWeightOperatorNormal
     !% Constructors for the ``normal'' output analysis class.
     module procedure normalConstructorParameters
     module procedure normalConstructorInternal
  end interface outputAnalysisWeightOperatorNormal

contains

  function normalConstructorParameters(parameters) result(self)
    !% Constructor for the ``normal'' output analysis weight operator class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type            (outputAnalysisWeightOperatorNormal  )                :: self
    type            (inputParameters                     ), intent(inout) :: parameters
    class           (outputAnalysisPropertyExtractorClass), pointer       :: outputAnalysisPropertyExtractor_
    class           (outputAnalysisPropertyOperatorClass ), pointer       :: outputAnalysisPropertyOperator_
    double precision                                                      :: rangeLower                      , rangeUpper, &
         &                                                                   rootVariance

    ! Check and read parameters.
    !# <inputParameter>
    !#   <name>rangeLower</name>
    !#   <source>parameters</source>
    !#   <description>Lower integration limit for the normal distribution weight operator.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>rangeUpper</name>
    !#   <source>parameters</source>
    !#   <description>Upper integration limit for the normal distribution weight operator.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>rootVariance</name>
    !#   <source>parameters</source>
    !#   <description>Root variance for the normal distribution weight operator.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="outputAnalysisPropertyExtractor"      name="outputAnalysisPropertyExtractor_"      source="parameters"          />
    !# <objectBuilder class="outputAnalysisPropertyOperator"       name="outputAnalysisPropertyOperator_"       source="parameters"          />
    self=outputAnalysisWeightOperatorNormal(rangeLower,rangeUpper,rootVariance,outputAnalysisPropertyExtractor_,outputAnalysisPropertyOperator_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function normalConstructorParameters

  function normalConstructorInternal(rangeLower,rangeUpper,rootVariance,outputAnalysisPropertyExtractor_,outputAnalysisPropertyOperator_) result (self)
    !% Internal constructor for the ``normal'' output analysis distribution operator class.
    use Input_Parameters
    implicit none
    type            (outputAnalysisWeightOperatorNormal  )                        :: self
    double precision                                      , intent(in   )         :: rangeLower                      , rangeUpper, &
         &                                                                           rootVariance
    class           (outputAnalysisPropertyExtractorClass), intent(in   ), target :: outputAnalysisPropertyExtractor_
    class           (outputAnalysisPropertyOperatorClass ), intent(in   ), target :: outputAnalysisPropertyOperator_
    !# <constructorAssign variables="rangeLower, rangeUpper, rootVariance, *outputAnalysisPropertyExtractor_, *outputAnalysisPropertyOperator_"/>

    return
  end function normalConstructorInternal

  subroutine normalDestructor(self)
    !% Destructor for  the ``normal'' output analysis weight operator class.
    type(outputAnalysisWeightOperatorNormal), intent(inout) :: self
    
    !# <objectDestructor name="self%outputAnalysisPropertyExtractor_" />
    !# <objectDestructor name="self%outputAnalysisPropertyOperator_"  />
    return
  end subroutine normalDestructor
  
  double precision function normalOperate(self,weightValue,node,propertyValue,propertyValueIntrinsic,propertyType,propertyQuantity,outputIndex)
    !% Implement an normal output analysis weight operator.
    use, intrinsic :: ISO_C_Binding
    use               Error_Functions
    implicit none
    class           (outputAnalysisWeightOperatorNormal), intent(inout) :: self
    type            (treeNode                          ), intent(inout) :: node
    double precision                                    , intent(in   ) :: propertyValue      , propertyValueIntrinsic, &
         &                                                                 weightValue
    integer                                             , intent(in   ) :: propertyType       , propertyQuantity
    integer         (c_size_t                          ), intent(in   ) :: outputIndex
    double precision                                                    :: normalPropertyValue
    integer                                                             :: normalPropertyType
    !GCC$ attributes unused :: propertyValue,propertyValueIntrinsic, propertyType, propertyQuantity

    ! Extract property and operate on it.
    normalPropertyType =self%outputAnalysisPropertyExtractor_%type   (                                                       )
    normalPropertyValue=self%outputAnalysisPropertyExtractor_%extract(                    node                               )
    normalPropertyValue=self%outputAnalysisPropertyOperator_ %operate(normalPropertyValue,node,normalPropertyType,outputIndex)
    ! Multiply weight by integral over a normal distribution between given limits, using the property value as the mean.
    if     (                                     &
         &   normalPropertyValue == +huge(0.0d0) &
         &  .or.                                 &
         &   normalPropertyValue == -huge(0.0d0) &
         & ) then
       normalOperate=+0.0d0
    else
       normalOperate=+weightValue                                                                           &
            &        *(                                                                                     &
            &          +Error_Function((self%rangeUpper-normalPropertyValue)/sqrt(2.0d0)/self%rootVariance) &
            &          -Error_Function((self%rangeLower-normalPropertyValue)/sqrt(2.0d0)/self%rootVariance) &
            &         )                                                                                     &
            &        /2.0d0
    end if
    return
  end function normalOperate
