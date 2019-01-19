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

!% Contains a module which implements a black hole mass output analysis property extractor class.

  !# <outputAnalysisPropertyExtractor name="outputAnalysisPropertyExtractorMassBlackHole">
  !#  <description>An ISM mass output analysis property extractor class.</description>
  !# </outputAnalysisPropertyExtractor>
  type, extends(outputAnalysisPropertyExtractorClass) :: outputAnalysisPropertyExtractorMassBlackHole
     !% A black hole mass output analysis class.
     private
   contains
     procedure :: extract  => massBlackHoleExtract
     procedure :: type     => massBlackHoleType
     procedure :: quantity => massBlackHoleQuantity
  end type outputAnalysisPropertyExtractorMassBlackHole

  interface outputAnalysisPropertyExtractorMassBlackHole
     !% Constructors for the ``massBlackHole'' output analysis class.
     module procedure massBlackHoleConstructorParameters
  end interface outputAnalysisPropertyExtractorMassBlackHole

contains

  function massBlackHoleConstructorParameters(parameters)
    !% Constructor for the ``massBlackHole'' output analysis property extractor class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type(outputAnalysisPropertyExtractorMassBlackHole)                :: massBlackHoleConstructorParameters
    type(inputParameters                             ), intent(inout) :: parameters
    !GCC$ attributes unused :: parameters

    massBlackHoleConstructorParameters=outputAnalysisPropertyExtractorMassBlackHole()
    return
  end function massBlackHoleConstructorParameters

  double precision function massBlackHoleExtract(self,node)
    !% Implement a massBlackHole output analysis.
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    use Galacticus_Nodes                  , only : nodeComponentBlackHole
    implicit none
    class(outputAnalysisPropertyExtractorMassBlackHole), intent(inout) :: self
    type (treeNode                                    ), intent(inout) :: node
    class(nodeComponentBlackHole                      ), pointer       :: blackHole
    !GCC$ attributes unused :: self

    blackHole            => node     %blackHole()
    massBlackHoleExtract =  blackHole%mass     ()
    return
  end function massBlackHoleExtract

  integer function massBlackHoleType(self)
    !% Return the type of the stellar mass property.
    use Output_Analyses_Options
    implicit none
    class(outputAnalysisPropertyExtractorMassBlackHole), intent(inout) :: self
    !GCC$ attributes unused :: self

    massBlackHoleType=outputAnalysisPropertyTypeLinear
    return
  end function massBlackHoleType

  integer function massBlackHoleQuantity(self)
    !% Return the class of the stellar luminosity property.
    use Output_Analyses_Options
    implicit none
    class(outputAnalysisPropertyExtractorMassBlackHole), intent(inout) :: self
    !GCC$ attributes unused :: self

    massBlackHoleQuantity=outputAnalysisPropertyQuantityMass
    return
  end function massBlackHoleQuantity
