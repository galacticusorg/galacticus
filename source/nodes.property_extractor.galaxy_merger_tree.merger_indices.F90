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

  !![
  <nodePropertyExtractor name="nodePropertyExtractorGalaxyMergerTreeMergerIndices">
   <description>
     A node property extractor which extracts the merger indices properties of galaxy merger trees.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorIntegerList) :: nodePropertyExtractorGalaxyMergerTreeMergerIndices
     !!{
     A property extractor which extracts the merger indices properties of galaxy merger trees.
     !!}
     private
     integer :: mergeIndexID, mergeTargetID
   contains
     procedure :: elementCount => galaxyMergerTreeMergerIndicesElementCount
     procedure :: extract      => galaxyMergerTreeMergerIndicesExtract
     procedure :: names        => galaxyMergerTreeMergerIndicesNames
     procedure :: descriptions => galaxyMergerTreeMergerIndicesDescriptions
     procedure :: unitsInSI    => galaxyMergerTreeMergerIndicesUnitsInSI
  end type nodePropertyExtractorGalaxyMergerTreeMergerIndices

  interface nodePropertyExtractorGalaxyMergerTreeMergerIndices
     !!{
     Constructors for the \refClass{nodePropertyExtractorGalaxyMergerTreeMergerIndices} output extractor class.
     !!}
     module procedure galaxyMergerTreeMergerIndicesConstructorParameters
     module procedure galaxyMergerTreeMergerIndicesConstructorInternal
  end interface nodePropertyExtractorGalaxyMergerTreeMergerIndices

contains

  function galaxyMergerTreeMergerIndicesConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorGalaxyMergerTreeMergerIndices} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(nodePropertyExtractorGalaxyMergerTreeMergerIndices)                :: self
    type(inputParameters                                   ), intent(inout) :: parameters

    self=nodePropertyExtractorGalaxyMergerTreeMergerIndices()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function galaxyMergerTreeMergerIndicesConstructorParameters

  function galaxyMergerTreeMergerIndicesConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorGalaxyMergerTreeMergerIndices} output extractor property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorGalaxyMergerTreeMergerIndices) :: self
    
    !![
    <addMetaProperty component="basic" type="longInteger" name="galaxyMergerTreeMergeIndex"  id="self%mergeIndexID"  rank="1" isCreator="no"/>
    <addMetaProperty component="basic" type="longInteger" name="galaxyMergerTreeMergeTarget" id="self%mergeTargetID" rank="1" isCreator="no"/>
    !!]
    return
  end function galaxyMergerTreeMergerIndicesConstructorInternal

  integer function galaxyMergerTreeMergerIndicesElementCount(self)
    !!{
    Return a count of the number of properties extracted.
    !!}
    implicit none
    class(nodePropertyExtractorGalaxyMergerTreeMergerIndices), intent(inout) :: self

    galaxyMergerTreeMergerIndicesElementCount=2
    return
  end function galaxyMergerTreeMergerIndicesElementCount

  function galaxyMergerTreeMergerIndicesExtract(self,node,instance) result(galaxyMergerTree)
    !!{
    Implement a galaxyMergerTreeMergerIndices output extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    integer(kind_int8                                         ), dimension(:,:), allocatable :: galaxyMergerTree
    class  (nodePropertyExtractorGalaxyMergerTreeMergerIndices), intent(inout)               :: self
    type   (treeNode                                          ), intent(inout)               :: node
    type   (multiCounter                                      ), intent(inout) , optional    :: instance
    class  (nodeComponentBasic                                )                , pointer     :: basic
    integer(kind_int8                                         ), dimension(:  ), allocatable :: mergeIndices    , mergeTargets
    !$GLC attributes unused :: instance
    !$GLC attributes initialized :: mergeIndices, mergeTargets

    basic        => node %basic                          (                  )
    mergeIndices =  basic%longIntegerRank1MetaPropertyGet(self% mergeIndexID)
    mergeTargets =  basic%longIntegerRank1MetaPropertyGet(self%mergeTargetID)
    allocate(galaxyMergerTree(size(mergeIndices),2))
    galaxyMergerTree(:,1)=mergeIndices
    galaxyMergerTree(:,2)=mergeTargets
    return
  end function galaxyMergerTreeMergerIndicesExtract
  
  subroutine galaxyMergerTreeMergerIndicesNames(self,names)
    !!{
    Return the names of the {\normalfont \ttfamily galaxyMergerTreeMergerIndices} properties.
    !!}
    implicit none
    class(nodePropertyExtractorGalaxyMergerTreeMergerIndices), intent(inout)                             :: self
    type (varying_string                                    ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self

    allocate(names(2))
    names(1)=var_str('galaxyMergerTreeMergeSatelliteIndex')
    names(2)=var_str('galaxyMergerTreeMergeCentralIndex'  )
    return
  end subroutine galaxyMergerTreeMergerIndicesNames

  subroutine galaxyMergerTreeMergerIndicesDescriptions(self,descriptions)
    !!{
    Return the descriptions of the {\normalfont \ttfamily galaxyMergerTreeMergerIndices} properties.
    !!}
    implicit none
    class(nodePropertyExtractorGalaxyMergerTreeMergerIndices), intent(inout)                             :: self
    type (varying_string                                    ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self

    allocate(descriptions(2))
    descriptions(1)=var_str('Galaxy merger tree merging satellite indices.')
    descriptions(2)=var_str('Galaxy merger tree merging central indices.'  )
    return
  end subroutine galaxyMergerTreeMergerIndicesDescriptions

  function galaxyMergerTreeMergerIndicesUnitsInSI(self) result(unitsInSI)
    !!{
    Return the units of the {\normalfont \ttfamily galaxyMergerTreeMergerIndices} properties in the SI system.
    !!}
    implicit none
    double precision                                                    , dimension(:) , allocatable :: unitsInSI
    class           (nodePropertyExtractorGalaxyMergerTreeMergerIndices), intent(inout)              :: self
    !$GLC attributes unused :: self

    allocate(unitsInSI(2))
    unitsInSI=0.0d0
    return
  end function galaxyMergerTreeMergerIndicesUnitsInSI
