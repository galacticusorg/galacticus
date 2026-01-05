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

  use :: Hashes, only : doubleHash, rank1DoubleHash

  !![
  <nodePropertyExtractor name="nodePropertyExtractorList" abstract="yes">
   <description>An abstract output analysis property extractor class which provides a list of floating point properties.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorClass), abstract :: nodePropertyExtractorList
     !!{
     A list property extractor.
     !!}
     private
   contains
     !![
     <methods>
       <method method="elementCount" description="Return a count of the number of properties extracted."              />
       <method method="extract"      description="Extract the properties from the given {\normalfont \ttfamily node}."/>
       <method method="names"        description="Return the name of the properties extracted."                       />
       <method method="descriptions" description="Return a description of the properties extracted."                  />
       <method method="unitsInSI"    description="Return the units of the properties extracted in the SI system."     />
       <method method="metaData"     description="Populate a hash with meta-data for the property."                   />
     </methods>
     !!]
     procedure(listElementCount), deferred :: elementCount
     procedure(listExtract     ), deferred :: extract
     procedure(listNames       ), deferred :: names
     procedure(listDescriptions), deferred :: descriptions
     procedure(listUnitsInSI   ), deferred :: unitsInSI
     procedure                             :: metaData     => listMetaData
  end type nodePropertyExtractorList

  abstract interface
     function listElementCount(self)
       !!{
       Interface for list property count.
       !!}
       import nodePropertyExtractorList
       integer                                           :: listElementCount
       class  (nodePropertyExtractorList), intent(inout) :: self
     end function listElementCount
  end interface

  abstract interface
     function listExtract(self,node,instance)
       !!{
       Interface for list property extraction.
       !!}
       import nodePropertyExtractorList, treeNode, multiCounter
       double precision                           , dimension(:,:), allocatable :: listExtract
       class           (nodePropertyExtractorList), intent(inout)               :: self
       type            (treeNode                 ), intent(inout)               :: node
       type            (multiCounter             ), intent(inout) , optional    :: instance
     end function listExtract
  end interface

  abstract interface
     subroutine listNames(self,names)
       !!{
       Interface for list names.
       !!}
       import varying_string, nodePropertyExtractorList
       class(nodePropertyExtractorList), intent(inout)                             :: self
       type (varying_string           ), intent(inout), dimension(:) , allocatable :: names
     end subroutine listNames
  end interface

  abstract interface
     subroutine listDescriptions(self,descriptions)
       !!{
       Interface for list descriptions.
       !!}
       import varying_string, nodePropertyExtractorList
       class(nodePropertyExtractorList), intent(inout)                             :: self
       type (varying_string           ), intent(inout), dimension(:) , allocatable :: descriptions
     end subroutine listDescriptions
  end interface

  abstract interface
     function listUnitsInSI(self)
       !!{
       Interface for list property units.
       !!}
       import nodePropertyExtractorList
       double precision                           , dimension(:) , allocatable :: listUnitsInSI
       class           (nodePropertyExtractorList), intent(inout)              :: self
     end function listUnitsInSI
  end interface

contains
  
  subroutine listMetaData(self,node,metaDataRank0,metaDataRank1)
    !!{
    Interface for list property meta-data.
    !!}
    implicit none
    class(nodePropertyExtractorList), intent(inout) :: self
    type (treeNode                 ), intent(inout) :: node
    type (doubleHash               ), intent(inout) :: metaDataRank0
    type (rank1DoubleHash          ), intent(inout) :: metaDataRank1
    !$GLC attributes unused :: self, node, metaDataRank0, metaDataRank1
    
    return
  end subroutine listMetaData
