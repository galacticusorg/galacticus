!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains a module which implements a weight operator class in which the weight is multiplied by an integral over the N-body
!% halo mass distribution.

  use Statistics_NBody_Halo_Mass_Errors

  !# <outputAnalysisWeightOperator name="outputAnalysisWeightOperatorNbodyMass">
  !#  <description>A weight operator class in which the weight is multiplied by an integral over the N-body halo mass distribution.</description>
  !# </outputAnalysisWeightOperator>
  type, extends(outputAnalysisWeightOperatorNormal) :: outputAnalysisWeightOperatorNbodyMass
     !% A weight operator class in which the weight is multiplied by an integral over the N-body halo mass distribution.
     private
     class(nbodyHaloMassErrorClass), pointer :: nbodyHaloMassError_ => null()
   contains
     final     ::                 nbodyMassDestructor
     procedure :: rootVariance => nbodyMassRootVariance
  end type outputAnalysisWeightOperatorNbodyMass

  interface outputAnalysisWeightOperatorNbodyMass
     !% Constructors for the ``nbodyMass'' output analysis class.
     module procedure nbodyMassConstructorParameters
     module procedure nbodyMassConstructorInternal
  end interface outputAnalysisWeightOperatorNbodyMass

contains

  function nbodyMassConstructorParameters(parameters) result(self)
    !% Constructor for the ``nbodyMass'' output analysis weight operator class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type            (outputAnalysisWeightOperatorNbodyMass)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    class           (outputAnalysisPropertyExtractorClass ), pointer       :: outputAnalysisPropertyExtractor_
    class           (outputAnalysisPropertyOperatorClass  ), pointer       :: outputAnalysisPropertyOperator_
    class           (nbodyHaloMassErrorClass              ), pointer       :: nbodyHaloMassError_
    double precision                                                       :: rangeLower                      , rangeUpper

    ! Check and read parameters.
    !# <inputParameter>
    !#   <name>rangeLower</name>
    !#   <source>parameters</source>
    !#   <description>Lower integration limit for the nbodyMass distribution weight operator.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>rangeUpper</name>
    !#   <source>parameters</source>
    !#   <description>Upper integration limit for the nbodyMass distribution weight operator.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="outputAnalysisPropertyExtractor" name="outputAnalysisPropertyExtractor_" source="parameters"/>
    !# <objectBuilder class="outputAnalysisPropertyOperator"  name="outputAnalysisPropertyOperator_"  source="parameters"/>
    !# <objectBuilder class="nbodyHaloMassError"              name="nbodyHaloMassError_"              source="parameters"/>
    self=outputAnalysisWeightOperatorNbodyMass(rangeLower,rangeUpper,outputAnalysisPropertyExtractor_,outputAnalysisPropertyOperator_,nbodyHaloMassError_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="outputAnalysisPropertyExtractor_"/>
    !# <objectDestructor name="outputAnalysisPropertyOperator_" />
    !# <objectDestructor name="nbodyHaloMassError_"             />
    return
  end function nbodyMassConstructorParameters

  function nbodyMassConstructorInternal(rangeLower,rangeUpper,outputAnalysisPropertyExtractor_,outputAnalysisPropertyOperator_,nbodyHaloMassError_) result (self)
    !% Internal constructor for the ``nbodyMass'' output analysis distribution operator class.
    use Input_Parameters
    implicit none
    type            (outputAnalysisWeightOperatorNbodyMass)                        :: self
    double precision                                       , intent(in   )         :: rangeLower                      , rangeUpper
    class           (outputAnalysisPropertyExtractorClass ), intent(in   ), target :: outputAnalysisPropertyExtractor_
    class           (outputAnalysisPropertyOperatorClass  ), intent(in   ), target :: outputAnalysisPropertyOperator_
    class(nbodyHaloMassErrorClass                         ), intent(in   ), target :: nbodyHaloMassError_
    !# <constructorAssign variables="rangeLower, rangeUpper, *outputAnalysisPropertyExtractor_, *outputAnalysisPropertyOperator_, *nbodyHaloMassError_"/>

    return
  end function nbodyMassConstructorInternal

  subroutine nbodyMassDestructor(self)
    !% Destructor for  the ``nbodyMass'' output analysis weight operator class.
    type(outputAnalysisWeightOperatorNbodyMass), intent(inout) :: self
    
    !# <objectDestructor name="self%nbodyHaloMassError_"             />
    !# <objectDestructor name="self%outputAnalysisPropertyExtractor_"/>
    !# <objectDestructor name="self%outputAnalysisPropertyOperator_" />
    return
  end subroutine nbodyMassDestructor

  double precision function nbodyMassRootVariance(self,node,propertyValue,propertyValueIntrinsic,propertyType,propertyQuantity,outputIndex)
    !% Return the root variance for use in the ``nbodyMass'' output analysis weight operator class.
    implicit none
    class           (outputAnalysisWeightOperatorNbodyMass), intent(inout) :: self
    type            (treeNode                             ), intent(inout) :: node
    double precision                                       , intent(in   ) :: propertyValue         , propertyValueIntrinsic
    integer                                                , intent(in   ) :: propertyType          , propertyQuantity
    integer         (c_size_t                             ), intent(in   ) :: outputIndex
    double precision                                                       :: nbodyMassPropertyValue
    integer                                                                :: nbodyMassPropertyType
    !GCC$ attributes unused :: propertyValue, propertyValueIntrinsic, propertyType, propertyQuantity, outputIndex

    ! Extract property and operate on it.
    nbodyMassPropertyType =+self%outputAnalysisPropertyExtractor_%type           (    )
    nbodyMassPropertyValue=+self%outputAnalysisPropertyExtractor_%extract        (node)
    nbodyMassRootVariance =+self%nbodyHaloMassError_             %errorFractional(node) &
         &                 *nbodyMassPropertyValue
    return
  end function nbodyMassRootVariance
  
