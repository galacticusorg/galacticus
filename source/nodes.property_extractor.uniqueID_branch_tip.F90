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
Implements a node branch tip index property extractor.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorUniqueIDBranchTip">
   <description>A node branch tip index property extractor.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorIntegerScalar) :: nodePropertyExtractorUniqueIDBranchTip
     !!{
     A node branch tip index property extractor.
     !!}
     private
     integer :: uniqueIDBranchTipID
   contains
     procedure :: extract     => uniqueIDBranchTipExtract
     procedure :: name        => uniqueIDBranchTipName
     procedure :: description => uniqueIDBranchTipDescription
  end type nodePropertyExtractorUniqueIDBranchTip

  interface nodePropertyExtractorUniqueIDBranchTip
     !!{
     Constructors for the \refClass{nodePropertyExtractorUniqueIDBranchTip} output analysis class.
     !!}
     module procedure uniqueIDBranchTipConstructorParameters
     module procedure uniqueIDBranchTipConstructorInternal
  end interface nodePropertyExtractorUniqueIDBranchTip

contains

  function uniqueIDBranchTipConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorUniqueIDBranchTip} node property extractor class which takes a parameter set as input.
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
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorUniqueIDBranchTip} node property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorUniqueIDBranchTip) :: self

    !![
    <addMetaProperty component="basic" name="nodeUniqueIDBranchTip" type="longInteger" id="self%uniqueIDBranchTipID" isCreator="no"/>
    !!]
    return
  end function uniqueIDBranchTipConstructorInternal

  function uniqueIDBranchTipExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily uniqueIDBranchTip} node property extractor.
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
    !!{
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
    !!{
    Return a description of the branch tip index property.
    !!}
    implicit none
    type (varying_string                        )                :: uniqueIDBranchTipDescription
    class(nodePropertyExtractorUniqueIDBranchTip), intent(inout) :: self
    !$GLC attributes unused :: self

    uniqueIDBranchTipDescription=var_str('Unique ID of the node at the tip of this branch.')
    return
  end function uniqueIDBranchTipDescription

