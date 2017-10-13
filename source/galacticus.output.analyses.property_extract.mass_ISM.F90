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

!% Contains a module which implements an ISM mass output analysis property extractor class.

  !# <outputAnalysisPropertyExtractor name="outputAnalysisPropertyExtractorMassISM">
  !#  <description>An ISM mass output analysis property extractor class.</description>
  !# </outputAnalysisPropertyExtractor>
  type, extends(outputAnalysisPropertyExtractorClass) :: outputAnalysisPropertyExtractorMassISM
     !% A stelalr mass output analysis class.
     private
   contains
     procedure :: extract  => massISMExtract
     procedure :: type     => massISMType
     procedure :: quantity => massISMQuantity
  end type outputAnalysisPropertyExtractorMassISM

  interface outputAnalysisPropertyExtractorMassISM
     !% Constructors for the ``massISM'' output analysis class.
     module procedure massISMConstructorParameters
  end interface outputAnalysisPropertyExtractorMassISM

contains

  function massISMConstructorParameters(parameters)
    !% Constructor for the ``massISM'' output analysis property extractor class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type(outputAnalysisPropertyExtractorMassISM)                :: massISMConstructorParameters
    type(inputParameters                      ), intent(inout) :: parameters
    !GCC$ attributes unused :: parameters

    massISMConstructorParameters=outputAnalysisPropertyExtractorMassISM()
    return
  end function massISMConstructorParameters

  double precision function massISMExtract(self,node)
    !% Implement a massISM output analysis.
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    implicit none
    class(outputAnalysisPropertyExtractorMassISM), intent(inout) :: self
    type (treeNode                              ), intent(inout) :: node
    !GCC$ attributes unused :: self

    massISMExtract=+Galactic_Structure_Enclosed_Mass(node,radiusLarge,massType=massTypeGaseous,componentType=componentTypeDisk    ) &
         &         +Galactic_Structure_Enclosed_Mass(node,radiusLarge,massType=massTypeGaseous,componentType=componentTypeSpheroid)
    return
  end function massISMExtract

  integer function massISMType(self)
    !% Return the type of the stellar mass property.
    use Output_Analyses_Options
    implicit none
    class(outputAnalysisPropertyExtractorMassISM), intent(inout) :: self
    !GCC$ attributes unused :: self

    massISMType=outputAnalysisPropertyTypeLinear
    return
  end function massISMType

  integer function massISMQuantity(self)
    !% Return the class of the stellar luminosity property.
    use Output_Analyses_Options
    implicit none
    class(outputAnalysisPropertyExtractorMassISM), intent(inout) :: self
    !GCC$ attributes unused :: self

    massISMQuantity=outputAnalysisPropertyQuantityMass
    return
  end function massISMQuantity
