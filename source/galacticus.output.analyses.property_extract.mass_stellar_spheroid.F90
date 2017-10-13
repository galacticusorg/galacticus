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

!% Contains a module which implements a spheroid stellar mass output analysis property extractor class.

  !# <outputAnalysisPropertyExtractor name="outputAnalysisPropertyExtractorMassStellarSpheroid">
  !#  <description>A spheroid stellar mass output analysis property extractor class.</description>
  !# </outputAnalysisPropertyExtractor>
  type, extends(outputAnalysisPropertyExtractorClass) :: outputAnalysisPropertyExtractorMassStellarSpheroid
     !% A stelalr mass output analysis class.
     private
   contains
     procedure :: extract  => massStellarSpheroidExtract
     procedure :: type     => massStellarSpheroidType
     procedure :: quantity => massStellarSpheroidQuantity
  end type outputAnalysisPropertyExtractorMassStellarSpheroid

  interface outputAnalysisPropertyExtractorMassStellarSpheroid
     !% Constructors for the ``massStellarSpheroid'' output analysis class.
     module procedure massStellarSpheroidConstructorParameters
  end interface outputAnalysisPropertyExtractorMassStellarSpheroid

contains

  function massStellarSpheroidConstructorParameters(parameters) result(self)
    !% Constructor for the ``massStellarSpheroid'' output analysis property extractor class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type(outputAnalysisPropertyExtractorMassStellarSpheroid)                :: self
    type(inputParameters                                   ), intent(inout) :: parameters
    !GCC$ attributes unused :: parameters

    self=outputAnalysisPropertyExtractorMassStellarSpheroid()
    return
  end function massStellarSpheroidConstructorParameters

  double precision function massStellarSpheroidExtract(self,node)
    !% Implement a stellar mass-weighted morphology output analysis.
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    implicit none
    class           (outputAnalysisPropertyExtractorMassStellarSpheroid), intent(inout) :: self
    type            (treeNode                                          ), intent(inout) :: node
    !GCC$ attributes unused :: self

    massStellarSpheroidExtract=Galactic_Structure_Enclosed_Mass(node,radiusLarge,massType=massTypeStellar,componentType=componentTypeSpheroid)
    return
  end function massStellarSpheroidExtract

  integer function massStellarSpheroidType(self)
    !% Return the type of the stellar mass-weighted morphology property.
    use Output_Analyses_Options
    implicit none
    class(outputAnalysisPropertyExtractorMassStellarSpheroid), intent(inout) :: self
    !GCC$ attributes unused :: self

    massStellarSpheroidType=outputAnalysisPropertyTypeLinear
    return
  end function massStellarSpheroidType

  integer function massStellarSpheroidQuantity(self)
    !% Return the class of the stellar luminosity property.
    use Output_Analyses_Options
    implicit none
    class(outputAnalysisPropertyExtractorMassStellarSpheroid), intent(inout) :: self
    !GCC$ attributes unused :: self

    massStellarSpheroidQuantity=outputAnalysisPropertyQuantityMass
    return
  end function massStellarSpheroidQuantity
