!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
Contains a module which implements an ISM mass output analysis property extractor class.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorIndicesTree">
   <description>An ISM mass output analysis property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorIntegerScalar) :: nodePropertyExtractorIndicesTree
     !!{
     A stelalr mass output analysis class.
     !!}
     private
   contains
     procedure :: extract     => indicesTreeExtract
     procedure :: type        => indicesTreeType
     procedure :: name        => indicesTreeName
     procedure :: description => indicesTreeDescription
  end type nodePropertyExtractorIndicesTree

  interface nodePropertyExtractorIndicesTree
     !!{
     Constructors for the ``indicesTree'' output analysis class.
     !!}
     module procedure indicesTreeConstructorParameters
  end interface nodePropertyExtractorIndicesTree

contains

  function indicesTreeConstructorParameters(parameters)
    !!{
    Constructor for the {\normalfont \ttfamily indicesTree} node property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorIndicesTree)                :: indicesTreeConstructorParameters
    type(inputParameters                 ), intent(inout) :: parameters

    indicesTreeConstructorParameters=nodePropertyExtractorIndicesTree()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function indicesTreeConstructorParameters

  function indicesTreeExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily indicesTree} node property extractor.
    !!}
    implicit none
    integer         (kind_int8                       )                          :: indicesTreeExtract
    class           (nodePropertyExtractorIndicesTree), intent(inout)           :: self
    type            (treeNode                        ), intent(inout), target   :: node
    double precision                                  , intent(in   )           :: time
    type            (multiCounter                    ), intent(inout), optional :: instance
    !$GLC attributes unused :: self, instance, time

    indicesTreeExtract=node%hostTree%index
    return
  end function indicesTreeExtract

  integer function indicesTreeType(self)
    !!{
    Return the type of the stellar mass property.
    !!}
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear
    implicit none
    class(nodePropertyExtractorIndicesTree), intent(inout) :: self
    !$GLC attributes unused :: self

    indicesTreeType=outputAnalysisPropertyTypeLinear
    return
  end function indicesTreeType

  function indicesTreeName(self)
    !!{
    Return the name of the indicesTree property.
    !!}
    implicit none
    type (varying_string                  )                :: indicesTreeName
    class(nodePropertyExtractorIndicesTree), intent(inout) :: self
    !$GLC attributes unused :: self

    indicesTreeName=var_str('mergerTreeIndex')
    return
  end function indicesTreeName

  function indicesTreeDescription(self)
    !!{
    Return a description of the indicesTree property.
    !!}
    implicit none
    type (varying_string                  )                :: indicesTreeDescription
    class(nodePropertyExtractorIndicesTree), intent(inout) :: self
    !$GLC attributes unused :: self

    indicesTreeDescription=var_str('Tree index for this node.')
    return
  end function indicesTreeDescription

