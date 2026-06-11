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

  !![
  <nodePropertyExtractor name="nodePropertyExtractorHierarchy" docformat="rst">
   <description>
   Extracts meta-properties describing the structural position of a node within the halo hierarchy, such as its depth level in the subhalo nesting, enabling analysis of multi-level substructure in merger trees.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorIntegerTuple) :: nodePropertyExtractorHierarchy
     !!{RST
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
     procedure :: units       => hierarchyUnits
  end type nodePropertyExtractorHierarchy

  interface nodePropertyExtractorHierarchy
     !!{RST
     Constructors for the ``nodePropertyExtractorHierarchy`` property extractor class.
     !!}
     module procedure hierarchyConstructorParameters
     module procedure hierarchyConstructorInternal
  end interface nodePropertyExtractorHierarchy

contains

  function hierarchyConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the ``nodePropertyExtractorHierarchy`` property extractor class which takes a parameter set as input.
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
    !!{RST
    Internal constructor for the ``nodePropertyExtractorHierarchy`` property extractor class.
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
    !!{RST
    Return the number of elements in the ``hierarchy`` property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorHierarchy), intent(inout) :: self
    double precision                                , intent(in   ) :: time
    !$GLC attributes unused :: self, time

    hierarchyElementCount=3
    return
  end function hierarchyElementCount

  function hierarchyExtract(self,node,time,instance)
    !!{RST
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
    !!{RST
    Return the names of the ``hierarchy`` properties.
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
    !!{RST
    Return the descriptions of the ``hierarchy`` properties.
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
    !!{RST
    Return the units of the ``hierarchy`` properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear
    implicit none
    double precision                                , dimension(:) , allocatable :: hierarchyUnitsInSI
    class           (nodePropertyExtractorHierarchy), intent(inout)              :: self
    double precision                                , intent(in   )              :: time
   !$GLC attributes unused :: self, time

    allocate(hierarchyUnitsInSI(3))
    hierarchyUnitsInSI=1.0d0
    return
  end function hierarchyUnitsInSI

  function hierarchyUnits(self,time) result(units)
    !!{RST
    Return the units of the hierarchy properties.
    !!}
    use :: Units_MetaData, only : unitType
    implicit none
    type            (unitType                      ), dimension(:), allocatable :: units
    class           (nodePropertyExtractorHierarchy), intent(inout)             :: self
    double precision                                , intent(in   )             :: time
    integer                                                                     :: i
    !$GLC attributes unused :: self, time

    allocate(units(3))
    do i=1,3
       units(i)=unitType(1.0d0)
    end do
    return
  end function hierarchyUnits
