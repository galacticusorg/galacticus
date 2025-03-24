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
  <nodePropertyExtractor name="nodePropertyExtractorNodeFormationTime">
   <description>
     A node property extractor which extracts the formation time of each node.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorNodeFormationTime
     !!{
     A property extractor which extracts the formation time of each node.
     !!}
     private
     integer :: nodeFormationTimeID
   contains
     procedure :: extract     => nodeFormationTimeExtract
     procedure :: name        => nodeFormationTimeName
     procedure :: description => nodeFormationTimeDescription
     procedure :: unitsInSI   => nodeFormationTimeUnitsInSI
  end type nodePropertyExtractorNodeFormationTime

  interface nodePropertyExtractorNodeFormationTime
     !!{
     Constructors for the {\normalfont \ttfamily nodeFormationTime} output extractor class.
     !!}
     module procedure nodeFormationTimeConstructorParameters
     module procedure nodeFormationTimeConstructorInternal
  end interface nodePropertyExtractorNodeFormationTime

contains

  function nodeFormationTimeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily nodeFormationTime} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(nodePropertyExtractorNodeFormationTime)                :: self
    type(inputParameters                       ), intent(inout) :: parameters

    self=nodePropertyExtractorNodeFormationTime()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function nodeFormationTimeConstructorParameters

  function nodeFormationTimeConstructorInternal() result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily nodeFormationTime} output extractor property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorNodeFormationTime) :: self
    
    !![
    <addMetaProperty component="basic" name="nodeFormationTime" id="self%nodeFormationTimeID" isEvolvable="no" isCreator="no"/>
    !!]
    return
  end function nodeFormationTimeConstructorInternal

  double precision function nodeFormationTimeExtract(self,node,instance)
    !!{
    Implement a nodeFormationTime output extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(nodePropertyExtractorNodeFormationTime), intent(inout), target   :: self
    type (treeNode                              ), intent(inout), target   :: node
    type (multiCounter                          ), intent(inout), optional :: instance
    class(nodeComponentBasic                    )               , pointer  :: basic
    !$GLC attributes unused :: instance

    basic                    => node %basic                    (                        )
    nodeFormationTimeExtract =  basic%floatRank0MetaPropertyGet(self%nodeFormationTimeID)
    return
  end function nodeFormationTimeExtract
  
  function nodeFormationTimeName(self)
    !!{
    Return the names of the {\normalfont \ttfamily nodeFormationTime} properties.
    !!}
    implicit none
    type (varying_string                        )                :: nodeFormationTimeName
    class(nodePropertyExtractorNodeFormationTime), intent(inout) :: self
    !$GLC attributes unused :: self

    nodeFormationTimeName=var_str('nodeFormationTime')
    return
  end function nodeFormationTimeName

  function nodeFormationTimeDescription(self)
    !!{
    Return the descriptions of the {\normalfont \ttfamily nodeFormationTime} properties.
    !!}
    implicit none
    type (varying_string                        )                :: nodeFormationTimeDescription
    class(nodePropertyExtractorNodeFormationTime), intent(inout) :: self
    !$GLC attributes unused :: self

    nodeFormationTimeDescription=var_str('Time of node formation.')
    return
  end function nodeFormationTimeDescription

  double precision function nodeFormationTimeUnitsInSI(self)
    !!{
    Return the units of the {\normalfont \ttfamily nodeFormationTime} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear
    implicit none
    class(nodePropertyExtractorNodeFormationTime), intent(inout) :: self
    !$GLC attributes unused :: self

    nodeFormationTimeUnitsInSI=gigaYear
    return
  end function nodeFormationTimeUnitsInSI
