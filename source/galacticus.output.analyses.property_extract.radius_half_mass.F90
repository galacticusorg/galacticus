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

!% Contains a module which implements a half-stellar mass radius output analysis property extractor class.


  !# <outputAnalysisPropertyExtractor name="outputAnalysisPropertyExtractorHalfMassRadius">
  !#  <description>A half-(stellar) mass output analysis property extractor class.</description>
  !# </outputAnalysisPropertyExtractor>
  type, extends(outputAnalysisPropertyExtractorClass) :: outputAnalysisPropertyExtractorHalfMassRadius
     !% A half-(stellar) mass property extractor output analysis class.
     private
   contains
     procedure :: extract  => halfMassRadiusExtract
     procedure :: type     => halfMassRadiusType
  end type outputAnalysisPropertyExtractorHalfMassRadius

  interface outputAnalysisPropertyExtractorHalfMassRadius
     !% Constructors for the ``halfMassRadius'' output analysis class.
     module procedure halfMassRadiusConstructorParameters
  end interface outputAnalysisPropertyExtractorHalfMassRadius

contains

  function halfMassRadiusConstructorParameters(parameters) result(self)
    !% Constructor for the ``halfMassRadius'' output analysis property extractor class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type(outputAnalysisPropertyExtractorHalfMassRadius)                :: self
    type(inputParameters                              ), intent(inout) :: parameters
    !GCC$ attributes unused :: parameters

    ! Build the object.
    self=outputAnalysisPropertyExtractorHalfMassRadius()
    return
  end function halfMassRadiusConstructorParameters

  double precision function halfMassRadiusExtract(self,node)
    !% Implement a half-mass output analysis.
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    implicit none
    class(outputAnalysisPropertyExtractorHalfMassRadius), intent(inout) :: self
    type (treeNode                                     ), intent(inout) :: node
    !GCC$ attributes unused :: self

    halfMassRadiusExtract=Galactic_Structure_Radius_Enclosing_Mass(node,fractionalMass=0.5d0,massType=massTypeStellar)
    return
  end function halfMassRadiusExtract

  integer function halfMassRadiusType(self)
    !% Return the type of the half-mass radius property.
    use Output_Analyses_Options
    implicit none
    class(outputAnalysisPropertyExtractorHalfMassRadius), intent(inout) :: self
    !GCC$ attributes unused :: self

    halfMassRadiusType=outputAnalysisPropertyTypeLinear
    return
  end function halfMassRadiusType
