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

!!{
Implements a node branch tip index property extractor.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorIndexBranchTip">
   <description>A node branch tip index property extractor.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorIntegerScalar) :: nodePropertyExtractorIndexBranchTip
     !!{
     A node branch tip index property extractor.
     !!}
     private
     integer :: indexBranchTipID
   contains
     procedure :: extract     => indexBranchTipExtract
     procedure :: name        => indexBranchTipName
     procedure :: description => indexBranchTipDescription
  end type nodePropertyExtractorIndexBranchTip

  interface nodePropertyExtractorIndexBranchTip
     !!{
     Constructors for the \refClass{nodePropertyExtractorIndexBranchTip} output analysis class.
     !!}
     module procedure indexBranchTipConstructorParameters
     module procedure indexBranchTipConstructorInternal
  end interface nodePropertyExtractorIndexBranchTip

contains

  function indexBranchTipConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorIndexBranchTip} node property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorIndexBranchTip)                :: self
    type(inputParameters                    ), intent(inout) :: parameters

    self=nodePropertyExtractorIndexBranchTip()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function indexBranchTipConstructorParameters

  function indexBranchTipConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorIndexBranchTip} node property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorIndexBranchTip) :: self

    !![
    <addMetaProperty component="basic" name="nodeIndexBranchTip" type="longInteger" id="self%indexBranchTipID" isCreator="no"/>
    !!]
    return
  end function indexBranchTipConstructorInternal

  function indexBranchTipExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily indexBranchTip} node property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    integer         (kind_int8                          )                          :: indexBranchTipExtract
    class           (nodePropertyExtractorIndexBranchTip), intent(inout)           :: self
    type            (treeNode                           ), intent(inout), target   :: node
    double precision                                     , intent(in   )           :: time
    type            (multiCounter                       ), intent(inout), optional :: instance
    class           (nodeComponentBasic                 )               , pointer  :: basic
    !$GLC attributes unused :: instance, time

    basic                 => node %basic                          (                     )
    indexBranchTipExtract =  basic%longIntegerRank0MetaPropertyGet(self%indexBranchTipID)
    return
  end function indexBranchTipExtract


  function indexBranchTipName(self)
    !!{
    Return the name of the branch tip index property.
    !!}
    implicit none
    type (varying_string                     )                :: indexBranchTipName
    class(nodePropertyExtractorIndexBranchTip), intent(inout) :: self
    !$GLC attributes unused :: self

    indexBranchTipName=var_str('nodeIndexBranchTip')
    return
  end function indexBranchTipName

  function indexBranchTipDescription(self)
    !!{
    Return a description of the branch tip index property.
    !!}
    implicit none
    type (varying_string                     )                :: indexBranchTipDescription
    class(nodePropertyExtractorIndexBranchTip), intent(inout) :: self
    !$GLC attributes unused :: self

    indexBranchTipDescription=var_str('Index of the node at the tip of this branch.')
    return
  end function indexBranchTipDescription

