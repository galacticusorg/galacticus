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

!% Contains a module which implements a ratio output analysis property extractor class.

  !# <outputAnalysisPropertyExtractor name="outputAnalysisPropertyExtractorRatio" defaultThreadPrivate="yes">
  !#  <description>A ratio output analysis property extractor class.</description>
  !# </outputAnalysisPropertyExtractor>
  type, extends(outputAnalysisPropertyExtractorClass) :: outputAnalysisPropertyExtractorRatio
     !% A ratio extractor output analysis class. This extractor extracts two other properties and takes their ratio.
     private
     class(outputAnalysisPropertyExtractorClass), pointer :: propertyNumerator_, propertyDenominator_
   contains
     final     ::             ratioDestructor
     procedure :: extract  => ratioExtract
     procedure :: type     => ratioType
  end type outputAnalysisPropertyExtractorRatio

  interface outputAnalysisPropertyExtractorRatio
     !% Constructors for the ``ratio'' output analysis class.
     module procedure ratioConstructorParameters
     module procedure ratioConstructorInternal
  end interface outputAnalysisPropertyExtractorRatio

contains

  function ratioConstructorParameters(parameters) result(self)
    !% Constructor for the ``ratio'' output analysis property extractor class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type (outputAnalysisPropertyExtractorRatio)                :: self
    type (inputParameters                     ), intent(inout) :: parameters
    class(outputAnalysisPropertyExtractorClass), pointer       :: propertyNumerator_, propertyDenominator_
    !GCC$ attributes unused :: parameters

    !# <objectBuilder class="outputAnalysisPropertyExtractor" name="propertyNumerator_"   parameterName="outputAnalysisPropertyExtractorNumerator"   source="parameters"/>
    !# <objectBuilder class="outputAnalysisPropertyExtractor" name="propertyDenominator_" parameterName="outputAnalysisPropertyExtractorDenominator" source="parameters"/>
    self=outputAnalysisPropertyExtractorRatio(propertyNumerator_,propertyDenominator_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function ratioConstructorParameters

  function ratioConstructorInternal(propertyNumerator_,propertyDenominator_) result(self)
    !% Internal constructor for the ``ratio'' output analysis property extractor class.
    use Input_Parameters
    implicit none
    type (outputAnalysisPropertyExtractorRatio)         :: self
    class(outputAnalysisPropertyExtractorClass), target :: propertyNumerator_, propertyDenominator_
    !# <constructorAssign variables="*propertyNumerator_, *propertyDenominator_"/>
    
    return
  end function ratioConstructorInternal
  
  subroutine ratioDestructor(self)
    !% Destructor for the ``ratio'' output analysis property extractor class.
    implicit none
    type(outputAnalysisPropertyExtractorRatio), intent(inout) :: self

    !# <objectDestructor name="self%propertyNumerator_"  />
    !# <objectDestructor name="self%propertyDenominator_"/>
    return
  end subroutine ratioDestructor

  double precision function ratioExtract(self,node)
    !% Implement a ratio output analysis.
    use Galacticus_Error
    implicit none
    class           (outputAnalysisPropertyExtractorRatio), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                                      :: numerator, denominator

    numerator   = self%propertyNumerator_  %extract(node)
    denominator = self%propertyDenominator_%extract(node)
    if (denominator == 0.0d0) then
       ratioExtract=+0.0d0
       if (numerator /= 0.0d0) call Galacticus_Error_Report('denominator is zero'//{introspection:location})
    else
       ratioExtract=+numerator   &
            &       /denominator
    end if
    return
  end function ratioExtract

  integer function ratioType(self)
    !% Return the type of the ratio property.
    use Output_Analyses_Options
    implicit none
    class(outputAnalysisPropertyExtractorRatio), intent(inout) :: self
    !GCC$ attributes unused :: self

    ratioType=outputAnalysisPropertyTypeLinear
    return
  end function ratioType
