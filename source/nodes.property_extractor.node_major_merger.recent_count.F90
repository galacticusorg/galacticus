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
  Implements a node property extractor which extracts the number of recent node major mergers.
  !!}

  use :: Output_Times, only : outputTimesClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorNodeMajorMergerRecentCount">
   <description>
   Implements a node property extractor which extracts the number of recent node major mergers.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorIntegerScalar) :: nodePropertyExtractorNodeMajorMergerRecentCount
     !!{
     A  node property extractor which extracts the number of recent node major mergers.
     !!}
     private
     class  (outputTimesClass), pointer :: outputTimes_                 => null()
     integer                            :: nodeMajorMergerRecentCountID
   contains
     final     ::                nodeMajorMergerRecentCountDestructor
     procedure :: extract     => nodeMajorMergerRecentCountExtract
     procedure :: name        => nodeMajorMergerRecentCountName
     procedure :: description => nodeMajorMergerRecentCountDescription
  end type nodePropertyExtractorNodeMajorMergerRecentCount

  interface nodePropertyExtractorNodeMajorMergerRecentCount
     !!{
     Constructors for the {\normalfont \ttfamily nodeMajorMergerRecentCount} output analysis class.
     !!}
     module procedure nodeMajorMergerRecentCountConstructorParameters
     module procedure nodeMajorMergerRecentCountConstructorInternal
  end interface nodePropertyExtractorNodeMajorMergerRecentCount

contains

  function nodeMajorMergerRecentCountConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily nodeMajorMergerRecentCount} node property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorNodeMajorMergerRecentCount)                :: self
    type (inputParameters                                ), intent(inout) :: parameters
    class(outputTimesClass                               ), pointer       :: outputTimes_

    !![
    <objectBuilder class="outputTimes" name="outputTimes_" source="parameters"/>
    !!]
    self=nodePropertyExtractorNodeMajorMergerRecentCount(outputTimes_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputTimes_"/>
    !!]
    return
  end function nodeMajorMergerRecentCountConstructorParameters

  function nodeMajorMergerRecentCountConstructorInternal(outputTimes_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily nodeMajorMergerRecentCount} node property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorNodeMajorMergerRecentCount)                        :: self
    class(outputTimesClass                               ), intent(in   ), target :: outputTimes_
    !![
    <constructorAssign variables="*outputTimes_"/>
    !!]
 
    !![
    <addMetaProperty component="basic" name="nodeMajorMergerRecentCount" id="self%nodeMajorMergerRecentCountID" type="integer" rank="1"/>
    !!]
    return
  end function nodeMajorMergerRecentCountConstructorInternal

  subroutine nodeMajorMergerRecentCountDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily nodeMajorMergerRecentCount} node operator class.
    !!}
    implicit none
    type(nodePropertyExtractorNodeMajorMergerRecentCount), intent(inout) :: self
     
    !![
    <objectDestructor name="self%outputTimes_"/>
    !!]
    return
  end subroutine nodeMajorMergerRecentCountDestructor

  function nodeMajorMergerRecentCountExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily nodeMajorMergerRecentCount} node property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    integer         (kind_int8                                      )                             :: nodeMajorMergerRecentCountExtract
    class           (nodePropertyExtractorNodeMajorMergerRecentCount), intent(inout)              :: self
    type            (treeNode                                       ), intent(inout), target      :: node
    double precision                                                 , intent(in   )              :: time
    type            (multiCounter                                   ), intent(inout), optional    :: instance
    class           (nodeComponentBasic                             )               , pointer     :: basic
    integer                                                          , dimension(:) , allocatable :: countMergersRecent
    !$GLC attributes unused :: instance

    allocate(countMergersRecent(self%outputTimes_%count()))
    basic                             => node %basic                      (                                 )
    countMergersRecent                =  basic%integerRank1MetaPropertyGet(self%nodeMajorMergerRecentCountID)
    nodeMajorMergerRecentCountExtract =  countMergersRecent(self%outputTimes_%index(time,findClosest=.true.))
   return
  end function nodeMajorMergerRecentCountExtract


  function nodeMajorMergerRecentCountName(self)
    !!{
    Return the name of the recent node major merger count property.
    !!}
    implicit none
    type (varying_string                                 )                :: nodeMajorMergerRecentCountName
    class(nodePropertyExtractorNodeMajorMergerRecentCount), intent(inout) :: self
    !$GLC attributes unused :: self

    nodeMajorMergerRecentCountName=var_str('nodeMajorMergerRecentCount')
    return
  end function nodeMajorMergerRecentCountName

  function nodeMajorMergerRecentCountDescription(self)
    !!{
    Return a description of the recent node major merger count property.
    !!}
    implicit none
    type (varying_string                                 )                :: nodeMajorMergerRecentCountDescription
    class(nodePropertyExtractorNodeMajorMergerRecentCount), intent(inout) :: self
    !$GLC attributes unused :: self

    nodeMajorMergerRecentCountDescription=var_str('The number of recent major mergers experienced by this node.')
    return
  end function nodeMajorMergerRecentCountDescription
