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

!% Contains a module which implements a stellar mass output analysis property extractor class.

  !# <outputAnalysisPropertyExtractor name="outputAnalysisPropertyExtractorMassStellar">
  !#  <description>A stellar mass output analysis property extractor class.</description>
  !# </outputAnalysisPropertyExtractor>
  type, extends(outputAnalysisPropertyExtractorClass) :: outputAnalysisPropertyExtractorMassStellar
     !% A stelalr mass output analysis class.
     private
   contains
     procedure :: extract  => massStellarExtract
  end type outputAnalysisPropertyExtractorMassStellar

  interface outputAnalysisPropertyExtractorMassStellar
     !% Constructors for the ``massStellar'' output analysis class.
     module procedure massStellarConstructorParameters
  end interface outputAnalysisPropertyExtractorMassStellar

contains

  function massStellarConstructorParameters(parameters)
    !% Constructor for the ``massStellar'' output analysis property extractor class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(outputAnalysisPropertyExtractorMassStellar)                :: massStellarConstructorParameters
    type(inputParameters                           ), intent(inout) :: parameters
    !GCC$ attributes unused :: parameters

    massStellarConstructorParameters=outputAnalysisPropertyExtractorMassStellar()
    return
  end function massStellarConstructorParameters

  double precision function massStellarExtract(self,node)
    !% Implement a massStellar output analysis.
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    implicit none
    class(outputAnalysisPropertyExtractorMassStellar), intent(inout) :: self
    type (treeNode                                  ), intent(inout) :: node
    !GCC$ attributes unused :: self

    massStellarExtract=Galactic_Structure_Enclosed_Mass(node,radiusLarge,massType=massTypeStellar)
    return
  end function massStellarExtract
