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
  Implements a node property extractor class that combines all extractors needed for building galaxy merger trees.
  !!}

  public :: nodePropertyExtractorGalaxyMergerTreeSet
  
  !![
  <nodePropertyExtractor name="nodePropertyExtractorGalaxyMergerTree">
   <description>An output extractor property extractor class that combines all extractors needed for building galaxy merger trees.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorMulti) :: nodePropertyExtractorGalaxyMergerTree
     !!{
     An output extractor property extractor class that combines all extractors needed for building galaxy merger trees.
     !!}
     private
  end type nodePropertyExtractorGalaxyMergerTree

  interface nodePropertyExtractorGalaxyMergerTree
     !!{
     Constructors for the \refClass{nodePropertyExtractorGalaxyMergerTree} output extractor class.
     !!}
     module procedure galaxyMergerTreeConstructorParameters
     module procedure galaxyMergerTreeConstructorInternal
  end interface nodePropertyExtractorGalaxyMergerTree

contains

  function galaxyMergerTreeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorGalaxyMergerTree} output extractor property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorGalaxyMergerTree)                :: self
    type(inputParameters                      ), intent(inout) :: parameters

    self=nodePropertyExtractorGalaxyMergerTree()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function galaxyMergerTreeConstructorParameters

  function galaxyMergerTreeConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorGalaxyMergerTree} output extractor property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorGalaxyMergerTree) :: self
    
    allocate(self%extractors               )
    allocate(self%extractors%next          )
    allocate(self%extractors%next%next     )
    allocate(self%extractors%next%next%next)
    allocate(nodePropertyExtractorGalaxyMergerTreeIndices        :: self%extractors               %extractor_)
    allocate(nodePropertyExtractorGalaxyMergerTreePhysical       :: self%extractors%next          %extractor_)
    allocate(nodePropertyExtractorGalaxyMergerTreeMergerIndices  :: self%extractors%next%next     %extractor_)
    allocate(nodePropertyExtractorGalaxyMergerTreeMergerPhysical :: self%extractors%next%next%next%extractor_)
    select type(extractor_ => self%extractors               %extractor_)
    type is (nodePropertyExtractorGalaxyMergerTreeIndices       )
       !![
       <referenceConstruct isResult="yes" owner="self" nameAssociated="extractor_" object="extractor_" constructor="nodePropertyExtractorGalaxyMergerTreeIndices       ()"/>
       !!]
    end select
    select type(extractor_ => self%extractors%next          %extractor_)
    type is (nodePropertyExtractorGalaxyMergerTreePhysical      )
       !![
       <referenceConstruct isResult="yes" owner="self" nameAssociated="extractor_" object="extractor_" constructor="nodePropertyExtractorGalaxyMergerTreePhysical      ()"/>
       !!]
    end select
    select type(extractor_ => self%extractors%next%next     %extractor_)
    type is (nodePropertyExtractorGalaxyMergerTreeMergerIndices )
       !![
       <referenceConstruct isResult="yes" owner="self" nameAssociated="extractor_" object="extractor_" constructor="nodePropertyExtractorGalaxyMergerTreeMergerIndices ()"/>
       !!]
    end select
    select type(extractor_ => self%extractors%next%next%next%extractor_)
    type is (nodePropertyExtractorGalaxyMergerTreeMergerPhysical)
       !![
       <referenceConstruct isResult="yes" owner="self" nameAssociated="extractor_" object="extractor_" constructor="nodePropertyExtractorGalaxyMergerTreeMergerPhysical()"/>
       !!]
    end select
    return
  end function galaxyMergerTreeConstructorInternal
