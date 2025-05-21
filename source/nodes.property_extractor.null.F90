!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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

!!{
Implements a null output analysis class.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorNull">
   <description>A null output analysis property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorClass) :: nodePropertyExtractorNull
     !!{
     A null output analysis class.
     !!}
     private
   contains
     procedure :: type => nullType
  end type nodePropertyExtractorNull

  interface nodePropertyExtractorNull
     !!{
     Constructors for the \refClass{nodePropertyExtractorNull} output analysis class.
     !!}
     module procedure nullConstructorParameters
  end interface nodePropertyExtractorNull

contains

  function nullConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorNull} output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorNull)                :: self
    type(inputParameters          ), intent(inout) :: parameters

    self=nodePropertyExtractorNull()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function nullConstructorParameters

  function nullType(self)
    !!{
    Return the type of the null property.
    !!}
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeUnknown
    implicit none
    type (enumerationOutputAnalysisPropertyTypeType)                :: nullType
    class(nodePropertyExtractorNull                ), intent(inout) :: self
    !$GLC attributes unused :: self

    nullType=outputAnalysisPropertyTypeUnknown
    return
  end function nullType
