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

!% Contains a module which implements a spin parameter output analysis property extractor class.

  !# <outputAnalysisPropertyExtractor name="outputAnalysisPropertyExtractorSpin">
  !#  <description>A spin parameter output analysis property extractor class.</description>
  !# </outputAnalysisPropertyExtractor>
  type, extends(outputAnalysisPropertyExtractorClass) :: outputAnalysisPropertyExtractorSpin
     !% A spin parameter property extractor output analysis class.
     private
   contains
     procedure :: extract  => spinExtract
     procedure :: type     => spinType
  end type outputAnalysisPropertyExtractorSpin

  interface outputAnalysisPropertyExtractorSpin
     !% Constructors for the ``spin'' output property extractor class.
     module procedure spinConstructorParameters
  end interface outputAnalysisPropertyExtractorSpin

contains

  function spinConstructorParameters(parameters) result(self)
    !% Constructor for the ``spin'' output analysis property extractor class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type (outputAnalysisPropertyExtractorSpin)                :: self
    type (inputParameters                    ), intent(inout) :: parameters
    !GCC$ attributes unused :: parameters

    ! Build the object.
    self=outputAnalysisPropertyExtractorSpin()
    return
  end function spinConstructorParameters

  double precision function spinExtract(self,node)
    !% Implement a spin output property extractor.
    implicit none
    class(outputAnalysisPropertyExtractorSpin), intent(inout) :: self
    type (treeNode                           ), intent(inout) :: node
    class(nodeComponentSpin                  ), pointer       :: spin
    !GCC$ attributes unused :: self

    spin        => node%spin()
    spinExtract =  spin%spin()
    return
  end function spinExtract

  integer function spinType(self)
    !% Return the type of the spin parameter property.
    use Output_Analyses_Options
    implicit none
    class(outputAnalysisPropertyExtractorSpin), intent(inout) :: self
    !GCC$ attributes unused :: self

    spinType=outputAnalysisPropertyTypeLinear
    return
  end function spinType
