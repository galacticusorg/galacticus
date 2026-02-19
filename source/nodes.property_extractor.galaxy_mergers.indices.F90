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
  <nodePropertyExtractor name="nodePropertyExtractorGalaxyMergersIndices">
   <description>
     A node property extractor which extracts the indices properties of galaxy-galaxy mergers.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorIntegerList) :: nodePropertyExtractorGalaxyMergersIndices
     !!{
     A property extractor which extracts the indices properties of galaxy-galaxy mergers.
     !!}
     private
     integer :: galaxyMergerSatelliteIndexID
   contains
     procedure :: elementCount => galaxyMergersIndicesElementCount
     procedure :: extract      => galaxyMergersIndicesExtract
     procedure :: names        => galaxyMergersIndicesNames
     procedure :: descriptions => galaxyMergersIndicesDescriptions
     procedure :: unitsInSI    => galaxyMergersIndicesUnitsInSI
  end type nodePropertyExtractorGalaxyMergersIndices

  interface nodePropertyExtractorGalaxyMergersIndices
     !!{
     Constructors for the \refClass{nodePropertyExtractorGalaxyMergersIndices} output extractor class.
     !!}
     module procedure galaxyMergersIndicesConstructorParameters
     module procedure galaxyMergersIndicesConstructorInternal
  end interface nodePropertyExtractorGalaxyMergersIndices

contains

  function galaxyMergersIndicesConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorGalaxyMergersIndices} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(nodePropertyExtractorGalaxyMergersIndices)                :: self
    type(inputParameters                          ), intent(inout) :: parameters

    self=nodePropertyExtractorGalaxyMergersIndices()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function galaxyMergersIndicesConstructorParameters

  function galaxyMergersIndicesConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorGalaxyMergersIndices} output extractor property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorGalaxyMergersIndices) :: self
    
    !![
    <addMetaProperty component="basic" name="galaxyMergerSatelliteIndex" id="self%galaxyMergerSatelliteIndexID" type="longInteger" rank="1" isCreator="no"/>
    !!]
    return
  end function galaxyMergersIndicesConstructorInternal

  integer function galaxyMergersIndicesElementCount(self)
    !!{
    Return a count of the number of properties extracted.
    !!}
    implicit none
    class(nodePropertyExtractorGalaxyMergersIndices), intent(inout) :: self

    galaxyMergersIndicesElementCount=1
    return
  end function galaxyMergersIndicesElementCount

  function galaxyMergersIndicesExtract(self,node,instance) result(galaxyMergers)
    !!{
    Implement a galaxyMergersIndices output extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    integer(kind_int8                                ), dimension(:,:), allocatable :: galaxyMergers
    class  (nodePropertyExtractorGalaxyMergersIndices), intent(inout)               :: self
    type   (treeNode                                 ), intent(inout)               :: node
    type   (multiCounter                             ), intent(inout) , optional    :: instance
    class  (nodeComponentBasic                       )                , pointer     :: basic
    integer(kind_int8                                ), dimension(:  ), allocatable :: indicesSatellite
    !$GLC attributes unused :: instance
    !$GLC attributes initialized :: indicesSatellite

    basic            => node %basic                          (                                 )
    indicesSatellite =  basic%longIntegerRank1MetaPropertyGet(self%galaxyMergerSatelliteIndexID)
    allocate(galaxyMergers(size(indicesSatellite),1))
    galaxyMergers(:,1)=indicesSatellite
    return
  end function galaxyMergersIndicesExtract
  
  subroutine galaxyMergersIndicesNames(self,names)
    !!{
    Return the names of the {\normalfont \ttfamily galaxyMergersIndices} properties.
    !!}
    implicit none
    class(nodePropertyExtractorGalaxyMergersIndices), intent(inout)                             :: self
    type (varying_string                           ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self

    allocate(names(1))
    names(1)=var_str('galaxyMergersMassSatelliteIndex')
    return
  end subroutine galaxyMergersIndicesNames

  subroutine galaxyMergersIndicesDescriptions(self,descriptions)
    !!{
    Return the descriptions of the {\normalfont \ttfamily galaxyMergersIndices} properties.
    !!}
    implicit none
    class(nodePropertyExtractorGalaxyMergersIndices), intent(inout)                             :: self
    type (varying_string                           ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self

    allocate(descriptions(1))
    descriptions(1)=var_str('Satellite galaxy index in galaxy-galaxy mergers.')
    return
  end subroutine galaxyMergersIndicesDescriptions

  function galaxyMergersIndicesUnitsInSI(self) result(unitsInSI)
    !!{
    Return the units of the {\normalfont \ttfamily galaxyMergersIndices} properties in the SI system.
    !!}
    implicit none
    double precision                                           , dimension(:) , allocatable :: unitsInSI
    class           (nodePropertyExtractorGalaxyMergersIndices), intent(inout)              :: self
    !$GLC attributes unused :: self

    allocate(unitsInSI(1))
    unitsInSI(1)=0.0d0
    return
  end function galaxyMergersIndicesUnitsInSI
