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
Implements a node branch tip index property extractor.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorUniqueIDBranchTip" docformat="rst">
   <description>
   Extracts the unique global identifier of the tip (earliest progenitor) node on the current merger tree branch, providing a persistent cross-snapshot identifier that enables tracking of branch origins across different output times and tree realizations.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorIntegerScalar) :: nodePropertyExtractorUniqueIDBranchTip
     !!{RST
     A node branch tip index property extractor.
     !!}
     private
     integer :: uniqueIDBranchTipID
   contains
     procedure :: extract     => uniqueIDBranchTipExtract
     procedure :: name        => uniqueIDBranchTipName
     procedure :: description => uniqueIDBranchTipDescription
     procedure :: units       => uniqueIDBranchTipUnits
  end type nodePropertyExtractorUniqueIDBranchTip

  interface nodePropertyExtractorUniqueIDBranchTip
     !!{RST
     Constructors for the :galacticus-class:`nodePropertyExtractorUniqueIDBranchTip` property extractor class.
     !!}
     module procedure uniqueIDBranchTipConstructorParameters
     module procedure uniqueIDBranchTipConstructorInternal
  end interface nodePropertyExtractorUniqueIDBranchTip

contains

  function uniqueIDBranchTipConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodePropertyExtractorUniqueIDBranchTip` property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorUniqueIDBranchTip)                :: self
    type(inputParameters                       ), intent(inout) :: parameters

    self=nodePropertyExtractorUniqueIDBranchTip()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function uniqueIDBranchTipConstructorParameters

  function uniqueIDBranchTipConstructorInternal() result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodePropertyExtractorUniqueIDBranchTip` property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorUniqueIDBranchTip) :: self

    !![
    <addMetaProperty component="basic" name="nodeUniqueIDBranchTip" type="longInteger" id="self%uniqueIDBranchTipID" isCreator="no"/>
    !!]
    return
  end function uniqueIDBranchTipConstructorInternal

  function uniqueIDBranchTipExtract(self,node,time,instance)
    !!{RST
    Implement a ``uniqueIDBranchTip`` node property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    integer         (kind_int8                             )                          :: uniqueIDBranchTipExtract
    class           (nodePropertyExtractorUniqueIDBranchTip), intent(inout)           :: self
    type            (treeNode                              ), intent(inout), target   :: node
    double precision                                        , intent(in   )           :: time
    type            (multiCounter                          ), intent(inout), optional :: instance
    class           (nodeComponentBasic                    )               , pointer  :: basic
    !$GLC attributes unused :: instance, time

    basic                 => node %basic                          (                     )
    uniqueIDBranchTipExtract =  basic%longIntegerRank0MetaPropertyGet(self%uniqueIDBranchTipID)
    return
  end function uniqueIDBranchTipExtract


  function uniqueIDBranchTipName(self)
    !!{RST
    Return the name of the branch tip index property.
    !!}
    implicit none
    type (varying_string                        )                :: uniqueIDBranchTipName
    class(nodePropertyExtractorUniqueIDBranchTip), intent(inout) :: self
    !$GLC attributes unused :: self

    uniqueIDBranchTipName=var_str('nodeUniqueIDBranchTip')
    return
  end function uniqueIDBranchTipName

  function uniqueIDBranchTipDescription(self)
    !!{RST
    Return a description of the branch tip index property.
    !!}
    implicit none
    type (varying_string                        )                :: uniqueIDBranchTipDescription
    class(nodePropertyExtractorUniqueIDBranchTip), intent(inout) :: self
    !$GLC attributes unused :: self

    uniqueIDBranchTipDescription=var_str('Unique ID of the node at the tip of this branch.')
    return
  end function uniqueIDBranchTipDescription

  function uniqueIDBranchTipUnits(self) result(units)
    !!{RST
    Return the units of the uniqueIDBranchTip property.
    !!}
    use :: Units_MetaData, only : unitType
    implicit none
    type (unitType                              )                :: units
    class(nodePropertyExtractorUniqueIDBranchTip), intent(inout) :: self
    !$GLC attributes unused :: self

    units=unitType(1.0d0)
    return
  end function uniqueIDBranchTipUnits
