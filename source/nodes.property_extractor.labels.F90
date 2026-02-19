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
  <nodePropertyExtractor name="nodePropertyExtractorLabels">
   <description>A node property extractor which extracts labels associated with nodes.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorIntegerTuple) :: nodePropertyExtractorLabels
     !!{
     A property extractor which extracts labels associated with nodes.
     !!}
     private
   contains
     procedure :: elementCount => labelsElementCount
     procedure :: extract      => labelsExtract
     procedure :: names        => labelsNames
     procedure :: descriptions => labelsDescriptions
     procedure :: unitsInSI    => labelsUnitsInSI
  end type nodePropertyExtractorLabels

  interface nodePropertyExtractorLabels
     !!{
     Constructors for the \refClass{nodePropertyExtractorLabels} output extractor class.
     !!}
     module procedure labelsConstructorParameters
  end interface nodePropertyExtractorLabels

contains

  function labelsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorLabels} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorLabels)                :: self
    type (inputParameters            ), intent(inout) :: parameters

    self=nodePropertyExtractorLabels()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function labelsConstructorParameters

  integer function labelsElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily labels} property extractors.
    !!}
    use :: Nodes_Labels, only : nodeLabelCount
    implicit none
    class           (nodePropertyExtractorLabels), intent(inout) :: self
    double precision                             , intent(in   ) :: time
    !$GLC attributes unused :: self, time

    labelsElementCount=nodeLabelCount()
    return
  end function labelsElementCount

  function labelsExtract(self,node,time,instance)
    !!{
    Implement a labels output extractor.
    !!}
    use :: Nodes_Labels, only : nodeLabelList
    implicit none
    integer         (kind_int8                  ), dimension(:) , allocatable :: labelsExtract
    class           (nodePropertyExtractorLabels), intent(inout)              :: self
    type            (treeNode                   ), intent(inout)              :: node
    double precision                             , intent(in   )              :: time
    type            (multiCounter               ), intent(inout), optional    :: instance
    !$GLC attributes unused :: self, time, instance

    call nodeLabelList(node,labelsExtract)
    return
  end function labelsExtract

  subroutine labelsNames(self,time,names)
    !!{
    Return the names of the {\normalfont \ttfamily labels} properties.
    !!}
    use :: Nodes_Labels   , only : nodeLabelNames
    use :: String_Handling, only : String_Upper_Case_First
    implicit none
    class           (nodePropertyExtractorLabels), intent(inout)                             :: self
    double precision                             , intent(in   )                             :: time
    type            (varying_string             ), intent(inout), dimension(:) , allocatable :: names
    integer                                                                                  :: i
    !$GLC attributes unused :: self, time

    call nodeLabelNames(names)
    do i=1,size(names)
       names(i)="nodeLabel"//String_Upper_Case_First(char(names(i)))
    end do
    return
  end subroutine labelsNames

  subroutine labelsDescriptions(self,time,descriptions)
    !!{
    Return the descriptions of the {\normalfont \ttfamily labels} properties.
    !!}
    use :: Nodes_Labels, only : nodeLabelDescriptions
    implicit none
    class           (nodePropertyExtractorLabels), intent(inout)                             :: self
    double precision                             , intent(in   )                             :: time
    type            (varying_string             ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self, time

    call nodeLabelDescriptions(descriptions)
    return
  end subroutine labelsDescriptions

  function labelsUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily labels} properties in the SI system.
    !!}
    implicit none
    double precision                             , dimension(:) , allocatable :: labelsUnitsInSI
    class           (nodePropertyExtractorLabels), intent(inout)              :: self
    double precision                             , intent(in   )              :: time

    allocate(labelsUnitsInSI(self%elementCount(time)))
    labelsUnitsInSI=0.0d0
    return
  end function labelsUnitsInSI

