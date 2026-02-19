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
  <nodePropertyExtractor name="nodePropertyExtractorNodeMajorMergerTime">
   <description>
     A node property extractor which extracts the time of the last major merger for each node.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorNodeMajorMergerTime
     !!{
     A property extractor which extracts the time of the last major merger for each node.
     !!}
     private
     integer :: nodeMajorMergerTimeID
   contains
     procedure :: extract     => nodeMajorMergerTimeExtract
     procedure :: name        => nodeMajorMergerTimeName
     procedure :: description => nodeMajorMergerTimeDescription
     procedure :: unitsInSI   => nodeMajorMergerTimeUnitsInSI
  end type nodePropertyExtractorNodeMajorMergerTime

  interface nodePropertyExtractorNodeMajorMergerTime
     !!{
     Constructors for the \refClass{nodePropertyExtractorNodeMajorMergerTime} output extractor class.
     !!}
     module procedure nodeMajorMergerTimeConstructorParameters
     module procedure nodeMajorMergerTimeConstructorInternal
  end interface nodePropertyExtractorNodeMajorMergerTime

contains

  function nodeMajorMergerTimeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorNodeMajorMergerTime} property extractor class which takes a parameter set as input.
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
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorNodeMajorMergerTime} output extractor property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorNodeMajorMergerTime) :: self
    
    !![
    <addMetaProperty component="basic" name="nodeMajorMergerTime" id="self%nodeMajorMergerTimeID"/>
    !!]
    return
  end function nodeMajorMergerTimeConstructorInternal

  double precision function nodeMajorMergerTimeExtract(self,node,instance)
    !!{
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
    !!{
    Return the names of the {\normalfont \ttfamily nodeMajorMergerTime} properties.
    !!}
    implicit none
    type (varying_string                          )                :: nodeMajorMergerTimeName
    class(nodePropertyExtractorNodeMajorMergerTime), intent(inout) :: self
    !$GLC attributes unused :: self

    nodeMajorMergerTimeName=var_str('nodeMajorMergerTime')
    return
  end function nodeMajorMergerTimeName

  function nodeMajorMergerTimeDescription(self)
    !!{
    Return the descriptions of the {\normalfont \ttfamily nodeMajorMergerTime} properties.
    !!}
    implicit none
    type (varying_string                          )                :: nodeMajorMergerTimeDescription
    class(nodePropertyExtractorNodeMajorMergerTime), intent(inout) :: self
    !$GLC attributes unused :: self

    nodeMajorMergerTimeDescription=var_str('Time of the last node major merger.')
    return
  end function nodeMajorMergerTimeDescription

  double precision function nodeMajorMergerTimeUnitsInSI(self)
    !!{
    Return the units of the {\normalfont \ttfamily nodeMajorMergerTime} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear
    implicit none
    class(nodePropertyExtractorNodeMajorMergerTime), intent(inout) :: self
    !$GLC attributes unused :: self

    nodeMajorMergerTimeUnitsInSI=gigaYear
    return
  end function nodeMajorMergerTimeUnitsInSI
