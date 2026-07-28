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

!!{RST
Implements a node unique ID property extractor.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorUniqueID" docformat="rst">
   <description>
   Extracts the internal unique global identifier (the {\normalfont \ttfamily uniqueID}) assigned to each node when it is created. This is the value drawn from the process-global unique-ID counter; under MPI it is guaranteed to be distinct across processes.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorIntegerScalar) :: nodePropertyExtractorUniqueID
     !!{RST
     A node unique ID property extractor.
     !!}
     private
   contains
     procedure :: extract     => uniqueIDExtract
     procedure :: name        => uniqueIDName
     procedure :: description => uniqueIDDescription
     procedure :: units       => uniqueIDUnits
  end type nodePropertyExtractorUniqueID

  interface nodePropertyExtractorUniqueID
     !!{RST
     Constructors for the :galacticus-class:`nodePropertyExtractorUniqueID` property extractor class.
     !!}
     module procedure uniqueIDConstructorParameters
  end interface nodePropertyExtractorUniqueID

contains

  function uniqueIDConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodePropertyExtractorUniqueID` property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorUniqueID)                :: self
    type(inputParameters              ), intent(inout) :: parameters

    self=nodePropertyExtractorUniqueID()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function uniqueIDConstructorParameters

  function uniqueIDExtract(self,node,time,instance)
    !!{RST
    Implement a ``uniqueID`` node property extractor.
    !!}
    implicit none
    integer         (kind_int8                    )                          :: uniqueIDExtract
    class           (nodePropertyExtractorUniqueID), intent(inout)           :: self
    type            (treeNode                     ), intent(inout), target   :: node
    double precision                               , intent(in   )           :: time
    type            (multiCounter                 ), intent(inout), optional :: instance
    !$GLC attributes unused :: self, time, instance

    uniqueIDExtract=node%uniqueID()
    return
  end function uniqueIDExtract

  function uniqueIDName(self)
    !!{RST
    Return the name of the unique ID property.
    !!}
    implicit none
    type (varying_string               )                :: uniqueIDName
    class(nodePropertyExtractorUniqueID), intent(inout) :: self
    !$GLC attributes unused :: self

    uniqueIDName=var_str('nodeUniqueID')
    return
  end function uniqueIDName

  function uniqueIDDescription(self)
    !!{RST
    Return a description of the unique ID property.
    !!}
    implicit none
    type (varying_string               )                :: uniqueIDDescription
    class(nodePropertyExtractorUniqueID), intent(inout) :: self
    !$GLC attributes unused :: self

    uniqueIDDescription=var_str('The internal unique global identifier of the node.')
    return
  end function uniqueIDDescription

  function uniqueIDUnits(self) result(units)
    !!{RST
    Return the units of the unique ID property.
    !!}
    use :: Units_MetaData, only : unitType
    implicit none
    type (unitType                     )                :: units
    class(nodePropertyExtractorUniqueID), intent(inout) :: self
    !$GLC attributes unused :: self

    units=unitType(1.0d0)
    return
  end function uniqueIDUnits
