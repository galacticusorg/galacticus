!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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
  Implements a node property extractor class that combines physical and index properties of galaxy mergers.
  !!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorGalaxyMergers">
   <description>An output extractor property extractor class that combines physical and index properties of galaxy mergers.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorMulti) :: nodePropertyExtractorGalaxyMergers
     !!{
     An output extractor property extractor class that combines physical and index properties of galaxy mergers.
     !!}
     private
   contains
  end type nodePropertyExtractorGalaxyMergers

  interface nodePropertyExtractorGalaxyMergers
     !!{
     Constructors for the \refClass{nodePropertyExtractorGalaxyMergers} output extractor class.
     !!}
     module procedure galaxyMergersConstructorParameters
     module procedure galaxyMergersConstructorInternal
  end interface nodePropertyExtractorGalaxyMergers

contains

  function galaxyMergersConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorGalaxyMergers} output extractor property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorGalaxyMergers)                :: self
    type(inputParameters                   ), intent(inout) :: parameters

    self=nodePropertyExtractorGalaxyMergers()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function galaxyMergersConstructorParameters

  function galaxyMergersConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorGalaxyMergers} output extractor property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorGalaxyMergers) :: self
    
    allocate(                                              self%extractors                )
    allocate(                                              self%extractors%next           )    
    allocate(nodePropertyExtractorGalaxyMergersIndices  :: self%extractors     %extractor_)
    allocate(nodePropertyExtractorGalaxyMergersPhysical :: self%extractors%next%extractor_)
    select type (extractor_ => self%extractors     %extractor_)
    type is (nodePropertyExtractorGalaxyMergersIndices )
       !![
       <referenceConstruct isResult="yes" owner="self" nameAssociated="extractor_" object="extractor_" constructor="nodePropertyExtractorGalaxyMergersIndices ()"/>
       !!]
    end select
    select type (extractor_ => self%extractors%next%extractor_)
    type is (nodePropertyExtractorGalaxyMergersPhysical)
       !![
       <referenceConstruct isResult="yes" owner="self" nameAssociated="extractor_" object="extractor_" constructor="nodePropertyExtractorGalaxyMergersPhysical()"/>
       !!]
    end select
    return
  end function galaxyMergersConstructorInternal
