!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Contains a module which implements a stellar mass-weighted morphology output analysis property extractor class.

  !# <outputAnalysisPropertyExtractor name="outputAnalysisPropertyExtractorMassStellarMorphology">
  !#  <description>A stellar mass-weighted morphology output analysis property extractor class.</description>
  !# </outputAnalysisPropertyExtractor>
  type, extends(outputAnalysisPropertyExtractorClass) :: outputAnalysisPropertyExtractorMassStellarMorphology
     !% A stelalr mass output analysis class.
     private
   contains
     procedure :: extract  => massStellarMorphologyExtract
     procedure :: type     => massStellarMorphologyType
  end type outputAnalysisPropertyExtractorMassStellarMorphology

  interface outputAnalysisPropertyExtractorMassStellarMorphology
     !% Constructors for the ``massStellarMorphology'' output analysis class.
     module procedure massStellarMorphologyConstructorParameters
  end interface outputAnalysisPropertyExtractorMassStellarMorphology

contains

  function massStellarMorphologyConstructorParameters(parameters) result(self)
    !% Constructor for the ``massStellarMorphology'' output analysis property extractor class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type(outputAnalysisPropertyExtractorMassStellarMorphology)                :: self
    type(inputParameters                                     ), intent(inout) :: parameters
    !GCC$ attributes unused :: parameters

    self=outputAnalysisPropertyExtractorMassStellarMorphology()
    return
  end function massStellarMorphologyConstructorParameters

  double precision function massStellarMorphologyExtract(self,node)
    !% Implement a stellar mass-weighted morphology output analysis.
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    implicit none
    class           (outputAnalysisPropertyExtractorMassStellarMorphology), intent(inout) :: self
    type            (treeNode                                            ), intent(inout) :: node
    double precision                                                                      :: massStellarDisk, massStellarSpheroid
    !GCC$ attributes unused :: self

    massStellarDisk    =Galactic_Structure_Enclosed_Mass(node,radiusLarge,massType=massTypeStellar,componentType=componentTypeDisk    )
    massStellarSpheroid=Galactic_Structure_Enclosed_Mass(node,radiusLarge,massType=massTypeStellar,componentType=componentTypeSpheroid)
    if (massStellarDisk+massStellarSpheroid > 0.0d0) then
       massStellarMorphologyExtract=+  massStellarSpheroid &
            &                       /(                     &
            &                         +massStellarDisk     &
            &                         +massStellarSpheroid &
            &                       )
    else
       massStellarMorphologyExtract=+0.0d0
    end if
    return
  end function massStellarMorphologyExtract

  integer function massStellarMorphologyType(self)
    !% Return the type of the stellar mass-weighted morphology property.
    use Output_Analyses_Options
    implicit none
    class(outputAnalysisPropertyExtractorMassStellarMorphology), intent(inout) :: self
    !GCC$ attributes unused :: self

    massStellarMorphologyType=outputAnalysisPropertyTypeLinear
    return
  end function massStellarMorphologyType
