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
  <nodePropertyExtractor name="nodePropertyExtractorHierarchy">
   <description>
     A node property extractor which extracts meta-properties related to the position of a node within the (sub-)halo hierarchy.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorIntegerTuple) :: nodePropertyExtractorHierarchy
     !!{
     A property extractor which extracts meta-properties related to the position of a node within the (sub-)halo hierarchy.
     !!}
     private
     integer :: nodeHierarchyLevelID       , nodeHierarchyLevelDepthID, &
          &     nodeHierarchyLevelMaximumID
   contains
     procedure :: elementCount => hierarchyElementCount
     procedure :: extract      => hierarchyExtract
     procedure :: names        => hierarchyNames
     procedure :: descriptions => hierarchyDescriptions
     procedure :: unitsInSI    => hierarchyUnitsInSI
  end type nodePropertyExtractorHierarchy

  interface nodePropertyExtractorHierarchy
     !!{
     Constructors for the \refClass{nodePropertyExtractorHierarchy} output extractor class.
     !!}
     module procedure hierarchyConstructorParameters
     module procedure hierarchyConstructorInternal
  end interface nodePropertyExtractorHierarchy

contains

  function hierarchyConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorHierarchy} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorHierarchy)                :: self
    type (inputParameters               ), intent(inout) :: parameters

    self=nodePropertyExtractorHierarchy()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function hierarchyConstructorParameters

  function hierarchyConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorHierarchy} output extractor property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorHierarchy) :: self

    !![
    <addMetaProperty component="basic" name="nodeHierarchyLevel"        type="integer" id="self%nodeHierarchyLevelID"        isCreator="no"/>
    <addMetaProperty component="basic" name="nodeHierarchyLevelDepth"   type="integer" id="self%nodeHierarchyLevelDepthID"   isCreator="no"/>
    <addMetaProperty component="basic" name="nodeHierarchyLevelMaximum" type="integer" id="self%nodeHierarchyLevelMaximumID" isCreator="no"/>
    !!]
    return
  end function hierarchyConstructorInternal

  integer function hierarchyElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily hierarchy} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorHierarchy), intent(inout) :: self
    double precision                                , intent(in   ) :: time
    !$GLC attributes unused :: self, time

    hierarchyElementCount=3
    return
  end function hierarchyElementCount

  function hierarchyExtract(self,node,time,instance)
    !!{
    Implement a hierarchy output extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    integer         (kind_int8                     ), dimension(:) , allocatable :: hierarchyExtract
    class           (nodePropertyExtractorHierarchy), intent(inout)              :: self
    type            (treeNode                      ), intent(inout)              :: node
    double precision                                , intent(in   )              :: time
    type            (multiCounter                  ), intent(inout), optional    :: instance
    class           (nodeComponentBasic            )               , pointer     :: basic
    !$GLC attributes unused :: time, instance

    basic => node%basic()
    allocate(hierarchyExtract(3))
    hierarchyExtract=[                                                                     &
         &            basic%integerRank0MetaPropertyGet(self%nodeHierarchyLevelID       ), &
         &            basic%integerRank0MetaPropertyGet(self%nodeHierarchyLevelDepthID  ), &
         &            basic%integerRank0MetaPropertyGet(self%nodeHierarchyLevelMaximumID)  &
         &           ]
    return
  end function hierarchyExtract

  subroutine hierarchyNames(self,time,names)
    !!{
    Return the names of the {\normalfont \ttfamily hierarchy} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorHierarchy), intent(inout)                             :: self
    double precision                                , intent(in   )                             :: time
    type            (varying_string                ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self, time

    allocate(names(3))
    names(1)=var_str('nodeHierarchyLevel'       )
    names(2)=var_str('nodeHierarchyLevelDepth'  )
    names(3)=var_str('nodeHierarchyLevelMaximum')
    return
  end subroutine hierarchyNames

  subroutine hierarchyDescriptions(self,time,descriptions)
    !!{
    Return the descriptions of the {\normalfont \ttfamily hierarchy} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorHierarchy), intent(inout)                             :: self
    double precision                                , intent(in   )                             :: time
    type            (varying_string                ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self, time

    allocate(descriptions(3))
    descriptions(1)=var_str('The current level of the node in the tree hierarchy (0 for a (non-sub-)halo; 1 for a sub-halo; 2 for a sub-sub-halo; etc.).'                               )
    descriptions(2)=var_str('Maximum level of the node in the tree hierarchy that could possibly ever be reached (0 for a (non-sub-)halo; 1 for a sub-halo; 2 for a sub-sub-halo; etc.)')
    descriptions(3)=var_str('Maximum level that the node ever reached in the tree hierarchy (0 for a (non-sub-)halo; 1 for a sub-halo; 2 for a sub-sub-halo; etc.).'                    )
    return
  end subroutine hierarchyDescriptions

  function hierarchyUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily hierarchy} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear
    implicit none
    double precision                                , dimension(:) , allocatable :: hierarchyUnitsInSI
    class           (nodePropertyExtractorHierarchy), intent(inout)              :: self
    double precision                                , intent(in   )              :: time
   !$GLC attributes unused :: self, time

    allocate(hierarchyUnitsInSI(3))
    hierarchyUnitsInSI=0.0d0
    return
  end function hierarchyUnitsInSI

