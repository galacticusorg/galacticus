!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a null output analysis class.

  !# <outputAnalysisPropertyExtractor name="outputAnalysisPropertyExtractorNull">
  !#  <description>A null output analysis property extractor class.</description>
  !# </outputAnalysisPropertyExtractor>
  type, extends(outputAnalysisPropertyExtractorClass) :: outputAnalysisPropertyExtractorNull
     !% A null output analysis class.
     private
   contains
     procedure :: extract  => nullExtract
  end type outputAnalysisPropertyExtractorNull

  interface outputAnalysisPropertyExtractorNull
     !% Constructors for the ``null'' output analysis class.
     module procedure nullConstructorParameters
  end interface outputAnalysisPropertyExtractorNull

contains

  function nullConstructorParameters(parameters)
    !% Constructor for the ``null'' output analysis property extractor class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(outputAnalysisPropertyExtractorNull)                :: nullConstructorParameters
    type(inputParameters                    ), intent(inout) :: parameters
    !GCC$ attributes unused :: parameters

    nullConstructorParameters=outputAnalysisPropertyExtractorNull()
    return
  end function nullConstructorParameters

  double precision function nullExtract(self,node)
    !% Implement a null output analysis.
    implicit none
    class(outputAnalysisPropertyExtractorNull), intent(inout) :: self
    type (treeNode                           ), intent(inout) :: node
    !GCC$ attributes unused :: self, node

    nullExtract=0.0d0
    return
  end function nullExtract
