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
  <nodePropertyExtractor name="nodePropertyExtractorGalaxyMergerTreeIndices">
   <description>
     A node property extractor which extracts the indices properties of galaxy merger trees.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorIntegerList) :: nodePropertyExtractorGalaxyMergerTreeIndices
     !!{
     A property extractor which extracts the indices properties of galaxy merger trees.
     !!}
     private
     integer :: nodeIndexID, countID
   contains
     procedure :: elementCount => galaxyMergerTreeIndicesElementCount
     procedure :: extract      => galaxyMergerTreeIndicesExtract
     procedure :: names        => galaxyMergerTreeIndicesNames
     procedure :: descriptions => galaxyMergerTreeIndicesDescriptions
     procedure :: unitsInSI    => galaxyMergerTreeIndicesUnitsInSI
  end type nodePropertyExtractorGalaxyMergerTreeIndices

  interface nodePropertyExtractorGalaxyMergerTreeIndices
     !!{
     Constructors for the {\normalfont \ttfamily galaxyMergerTreeIndices} output extractor class.
     !!}
     module procedure galaxyMergerTreeIndicesConstructorParameters
     module procedure galaxyMergerTreeIndicesConstructorInternal
  end interface nodePropertyExtractorGalaxyMergerTreeIndices

contains

  function galaxyMergerTreeIndicesConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily galaxyMergerTreeIndices} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(nodePropertyExtractorGalaxyMergerTreeIndices)                :: self
    type(inputParameters                             ), intent(inout) :: parameters

    self=nodePropertyExtractorGalaxyMergerTreeIndices()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function galaxyMergerTreeIndicesConstructorParameters

  function galaxyMergerTreeIndicesConstructorInternal() result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily galaxyMergerTreeIndices} output extractor property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorGalaxyMergerTreeIndices) :: self
    
    !![
    <addMetaProperty component="basic" type="longInteger" name="galaxyMergerTreeNodeIndex" id="self%nodeIndexID" rank="1" isCreator="no"/>
    <addMetaProperty component="basic" type="longInteger" name="galaxyMergerTreeCount"     id="self%countID"     rank="1" isCreator="no"/>
    !!]
    return
  end function galaxyMergerTreeIndicesConstructorInternal

  integer function galaxyMergerTreeIndicesElementCount(self)
    !!{
    Return a count of the number of properties extracted.
    !!}
    implicit none
    class(nodePropertyExtractorGalaxyMergerTreeIndices), intent(inout) :: self

    galaxyMergerTreeIndicesElementCount=2
    return
  end function galaxyMergerTreeIndicesElementCount

  function galaxyMergerTreeIndicesExtract(self,node,instance) result(galaxyMergerTree)
    !!{
    Implement a galaxyMergerTreeIndices output extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    integer(kind_int8                                   ), dimension(:,:), allocatable :: galaxyMergerTree
    class  (nodePropertyExtractorGalaxyMergerTreeIndices), intent(inout)               :: self
    type   (treeNode                                    ), intent(inout)               :: node
    type   (multiCounter                                ), intent(inout) , optional    :: instance
    class  (nodeComponentBasic                          )                , pointer     :: basic
    integer(kind_int8                                   ), dimension(:  ), allocatable :: nodeIndices     , counts
    !$GLC attributes unused :: instance
    !$GLC attributes initialized :: nodeIndices, counts

    basic       => node %basic                          (                )
    nodeIndices =  basic%longIntegerRank1MetaPropertyGet(self%nodeIndexID)
    counts      =  basic%longIntegerRank1MetaPropertyGet(self%    countID)
    allocate(galaxyMergerTree(size(nodeIndices),2))
    galaxyMergerTree(:,1)=nodeIndices
    galaxyMergerTree(:,2)=counts
    return
  end function galaxyMergerTreeIndicesExtract
  
  subroutine galaxyMergerTreeIndicesNames(self,names)
    !!{
    Return the names of the {\normalfont \ttfamily galaxyMergerTreeIndices} properties.
    !!}
    implicit none
    class(nodePropertyExtractorGalaxyMergerTreeIndices), intent(inout)                             :: self
    type (varying_string                              ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self

    allocate(names(2))
    names(1)=var_str('galaxyMergerTreeNodeIndex')
    names(2)=var_str('galaxyMergerTreeCount'    )
    return
  end subroutine galaxyMergerTreeIndicesNames

  subroutine galaxyMergerTreeIndicesDescriptions(self,descriptions)
    !!{
    Return the descriptions of the {\normalfont \ttfamily galaxyMergerTreeIndices} properties.
    !!}
    implicit none
    class(nodePropertyExtractorGalaxyMergerTreeIndices), intent(inout)                             :: self
    type (varying_string                              ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self

    allocate(descriptions(2))
    descriptions(1)=var_str('Galaxy merger tree node indices.'      )
    descriptions(2)=var_str('Galaxy merger tree sample time counts.')
    return
  end subroutine galaxyMergerTreeIndicesDescriptions

  function galaxyMergerTreeIndicesUnitsInSI(self) result(unitsInSI)
    !!{
    Return the units of the {\normalfont \ttfamily galaxyMergerTreeIndices} properties in the SI system.
    !!}
    implicit none
    double precision                                              , dimension(:) , allocatable :: unitsInSI
    class           (nodePropertyExtractorGalaxyMergerTreeIndices), intent(inout)              :: self
    !$GLC attributes unused :: self

    allocate(unitsInSI(2))
    unitsInSI=0.0d0
    return
  end function galaxyMergerTreeIndicesUnitsInSI
