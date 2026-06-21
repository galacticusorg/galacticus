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
  <nodePropertyExtractor name="nodePropertyExtractorGalaxyMergerTreeMergerPhysical" docformat="rst">
   <description>
   Extracts physical (floating-point) properties associated with merger events in galaxy merger trees, such as masses, mass ratios, and times of mergers between progenitor galaxies.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorList) :: nodePropertyExtractorGalaxyMergerTreeMergerPhysical
     !!{RST
     A property extractor which extracts the physical properties of galaxy merger trees.
     !!}
     private
     integer :: timeMergeID
   contains
     procedure :: elementCount => galaxyMergerTreeMergerPhysicalElementCount
     procedure :: extract      => galaxyMergerTreeMergerPhysicalExtract
     procedure :: names        => galaxyMergerTreeMergerPhysicalNames
     procedure :: descriptions => galaxyMergerTreeMergerPhysicalDescriptions
     procedure :: unitsInSI    => galaxyMergerTreeMergerPhysicalUnitsInSI
     procedure :: units        => galaxyMergerTreeMergerPhysicalUnits
  end type nodePropertyExtractorGalaxyMergerTreeMergerPhysical

  interface nodePropertyExtractorGalaxyMergerTreeMergerPhysical
     !!{RST
     Constructors for the :galacticus-class:`nodePropertyExtractorGalaxyMergerTreeMergerPhysical` property extractor class.
     !!}
     module procedure galaxyMergerTreeMergerPhysicalConstructorParameters
     module procedure galaxyMergerTreeMergerPhysicalConstructorInternal
  end interface nodePropertyExtractorGalaxyMergerTreeMergerPhysical

contains

  function galaxyMergerTreeMergerPhysicalConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodePropertyExtractorGalaxyMergerTreeMergerPhysical` property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(nodePropertyExtractorGalaxyMergerTreeMergerPhysical)                :: self
    type(inputParameters                                    ), intent(inout) :: parameters

    self=nodePropertyExtractorGalaxyMergerTreeMergerPhysical()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function galaxyMergerTreeMergerPhysicalConstructorParameters

  function galaxyMergerTreeMergerPhysicalConstructorInternal() result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodePropertyExtractorGalaxyMergerTreeMergerPhysical` property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorGalaxyMergerTreeMergerPhysical) :: self
    
    !![
    <addMetaProperty component="basic" name="galaxyMergerTreeTimeMerge" id="self%timeMergeID" rank="1" isEvolvable="no" isCreator="no"/>
      !!]
    return
  end function galaxyMergerTreeMergerPhysicalConstructorInternal

  integer function galaxyMergerTreeMergerPhysicalElementCount(self)
    !!{RST
    Return a count of the number of properties extracted.
    !!}
    implicit none
    class(nodePropertyExtractorGalaxyMergerTreeMergerPhysical), intent(inout) :: self

    galaxyMergerTreeMergerPhysicalElementCount=1
    return
  end function galaxyMergerTreeMergerPhysicalElementCount

  function galaxyMergerTreeMergerPhysicalExtract(self,node,instance) result(galaxyMergerTree)
    !!{RST
    Implement a galaxyMergerTreeMergerPhysical output extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    double precision                                                     , dimension(:,:), allocatable :: galaxyMergerTree
    class           (nodePropertyExtractorGalaxyMergerTreeMergerPhysical), intent(inout)               :: self
    type            (treeNode                                           ), intent(inout)               :: node
    type            (multiCounter                                       ), intent(inout) , optional    :: instance
    class           (nodeComponentBasic                                 )                , pointer     :: basic
    double precision                                                     , dimension(:  ), allocatable :: timesMerge
    !$GLC attributes unused :: instance
    !$GLC attributes initialized :: timesMerge
    
    basic      => node %basic                    (                )
    timesMerge =  basic%floatRank1MetaPropertyGet(self%timeMergeID)
    allocate(galaxyMergerTree(size(timesMerge),1))
    galaxyMergerTree(:,1)=timesMerge
    return
  end function galaxyMergerTreeMergerPhysicalExtract
  
  subroutine galaxyMergerTreeMergerPhysicalNames(self,names)
    !!{RST
    Return the names of the ``galaxyMergerTreeMergerPhysical`` properties.
    !!}
    implicit none
    class(nodePropertyExtractorGalaxyMergerTreeMergerPhysical), intent(inout)                             :: self
    type (varying_string                                     ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self

    allocate(names(1))
    names(1)=var_str('galaxyMergerTreeMergeTime'                )
    return
  end subroutine galaxyMergerTreeMergerPhysicalNames

  subroutine galaxyMergerTreeMergerPhysicalDescriptions(self,descriptions)
    !!{RST
    Return the descriptions of the ``galaxyMergerTreeMergerPhysical`` properties.
    !!}
    implicit none
    class(nodePropertyExtractorGalaxyMergerTreeMergerPhysical), intent(inout)                             :: self
    type (varying_string                                     ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self

    allocate(descriptions(1))
    descriptions(1)=var_str('Merger times in the galaxy merger tree.')
    return
  end subroutine galaxyMergerTreeMergerPhysicalDescriptions

  function galaxyMergerTreeMergerPhysicalUnitsInSI(self) result(unitsInSI)
    !!{RST
    Return the units of the ``galaxyMergerTreeMergerPhysical`` properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear
    implicit none
    double precision                                                     , dimension(:) , allocatable :: unitsInSI
    class           (nodePropertyExtractorGalaxyMergerTreeMergerPhysical), intent(inout)              :: self
    !$GLC attributes unused :: self

    allocate(unitsInSI(1))
    unitsInSI(1)=gigaYear
    return
  end function galaxyMergerTreeMergerPhysicalUnitsInSI

  function galaxyMergerTreeMergerPhysicalUnits(self) result(units)
    !!{RST
    Return the units of the galaxyMergerTreeMergerPhysical properties.
    !!}
    use :: Units_MetaData                  , only : unitType
    use :: Numerical_Constants_Astronomical, only : gigaYear
    implicit none
    type (unitType                                           ), dimension(:), allocatable :: units
    class(nodePropertyExtractorGalaxyMergerTreeMergerPhysical), intent(inout)             :: self
    !$GLC attributes unused :: self

    allocate(units(1))
    units(1)=unitType(gigaYear,description='Gyr',quantity='Gyr')
    return
  end function galaxyMergerTreeMergerPhysicalUnits
