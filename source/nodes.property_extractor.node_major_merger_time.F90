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
  <nodePropertyExtractor name="nodePropertyExtractorNodeMajorMergerTime" docformat="rst">
   <description>
   Extracts the cosmic time of the most recent major halo merger event for each node, where major mergers are defined by a configurable mass ratio threshold applied to the merging halo pair.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorNodeMajorMergerTime
     !!{RST
     A property extractor which extracts the time of the last major merger for each node.
     !!}
     private
     integer :: nodeMajorMergerTimeID
   contains
     procedure :: extract     => nodeMajorMergerTimeExtract
     procedure :: name        => nodeMajorMergerTimeName
     procedure :: description => nodeMajorMergerTimeDescription
     procedure :: unitsInSI   => nodeMajorMergerTimeUnitsInSI
     procedure :: units       => nodeMajorMergerTimeUnits
  end type nodePropertyExtractorNodeMajorMergerTime

  interface nodePropertyExtractorNodeMajorMergerTime
     !!{RST
     Constructors for the ``nodePropertyExtractorNodeMajorMergerTime`` property extractor class.
     !!}
     module procedure nodeMajorMergerTimeConstructorParameters
     module procedure nodeMajorMergerTimeConstructorInternal
  end interface nodePropertyExtractorNodeMajorMergerTime

contains

  function nodeMajorMergerTimeConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the ``nodePropertyExtractorNodeMajorMergerTime`` property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(nodePropertyExtractorNodeMajorMergerTime)                :: self
    type(inputParameters                         ), intent(inout) :: parameters

    self=nodePropertyExtractorNodeMajorMergerTime()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function nodeMajorMergerTimeConstructorParameters

  function nodeMajorMergerTimeConstructorInternal() result(self)
    !!{RST
    Internal constructor for the ``nodePropertyExtractorNodeMajorMergerTime`` property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorNodeMajorMergerTime) :: self
    
    !![
    <addMetaProperty component="basic" name="nodeMajorMergerTime" id="self%nodeMajorMergerTimeID"/>
    !!]
    return
  end function nodeMajorMergerTimeConstructorInternal

  double precision function nodeMajorMergerTimeExtract(self,node,instance)
    !!{RST
    Implement a nodeMajorMergerTime output extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(nodePropertyExtractorNodeMajorMergerTime), intent(inout), target   :: self
    type (treeNode                                ), intent(inout), target   :: node
    type (multiCounter                            ), intent(inout), optional :: instance
    class(nodeComponentBasic                      )               , pointer  :: basic
    !$GLC attributes unused :: instance

    basic                      => node %basic                    (                          )
    nodeMajorMergerTimeExtract =  basic%floatRank0MetaPropertyGet(self%nodeMajorMergerTimeID)
    return
  end function nodeMajorMergerTimeExtract
  
  function nodeMajorMergerTimeName(self)
    !!{RST
    Return the names of the ``nodeMajorMergerTime`` properties.
    !!}
    implicit none
    type (varying_string                          )                :: nodeMajorMergerTimeName
    class(nodePropertyExtractorNodeMajorMergerTime), intent(inout) :: self
    !$GLC attributes unused :: self

    nodeMajorMergerTimeName=var_str('nodeMajorMergerTime')
    return
  end function nodeMajorMergerTimeName

  function nodeMajorMergerTimeDescription(self)
    !!{RST
    Return the descriptions of the ``nodeMajorMergerTime`` properties.
    !!}
    implicit none
    type (varying_string                          )                :: nodeMajorMergerTimeDescription
    class(nodePropertyExtractorNodeMajorMergerTime), intent(inout) :: self
    !$GLC attributes unused :: self

    nodeMajorMergerTimeDescription=var_str('Time of the last node major merger.')
    return
  end function nodeMajorMergerTimeDescription

  double precision function nodeMajorMergerTimeUnitsInSI(self)
    !!{RST
    Return the units of the ``nodeMajorMergerTime`` properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear
    implicit none
    class(nodePropertyExtractorNodeMajorMergerTime), intent(inout) :: self
    !$GLC attributes unused :: self

    nodeMajorMergerTimeUnitsInSI=gigaYear
    return
  end function nodeMajorMergerTimeUnitsInSI

  function nodeMajorMergerTimeUnits(self) result(units)
    !!{RST
    Return the units of the node major merger time property.
    !!}
    use :: Units_MetaData, only : unitType
    implicit none
    type (unitType                                )                :: units
    class(nodePropertyExtractorNodeMajorMergerTime), intent(inout) :: self
    !$GLC attributes unused :: self

    units=unitType(self%unitsInSI(),description='Gyr',quantity='Gyr')
    return
  end function nodeMajorMergerTimeUnits
