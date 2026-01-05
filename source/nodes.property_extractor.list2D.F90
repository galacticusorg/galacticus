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
  <nodePropertyExtractor name="nodePropertyExtractorList2D" abstract="yes">
   <description>An abstract output analysis property extractor class which provides a 2D list of floating point properties.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorClass), abstract :: nodePropertyExtractorList2D
     !!{
     A 2D list property extractor.
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
     procedure(list2DElementCount), deferred :: elementCount
     procedure(list2DExtract     ), deferred :: extract
     procedure(list2DNames       ), deferred :: names
     procedure(list2DDescriptions), deferred :: descriptions
     procedure(list2DUnitsInSI   ), deferred :: unitsInSI
     procedure                               :: metaData     => list2DMetaData
  end type nodePropertyExtractorList2D

  abstract interface
     function list2DElementCount(self)
       !!{
       Interface for list2D property count.
       !!}
       import nodePropertyExtractorList2D
       integer                                             :: list2DElementCount
       class  (nodePropertyExtractorList2D), intent(inout) :: self
     end function list2DElementCount
  end interface

  abstract interface
     function list2DExtract(self,node,instance)
       !!{
       Interface for list2D property extraction.
       !!}
       import nodePropertyExtractorList2D, treeNode, multiCounter
       double precision                             , dimension(:,:,:), allocatable :: list2DExtract
       class           (nodePropertyExtractorList2D), intent(inout)                 :: self
       type            (treeNode                   ), intent(inout)                 :: node
       type            (multiCounter               ), intent(inout)   , optional    :: instance
     end function list2DExtract
  end interface

  abstract interface
     subroutine list2DNames(self,names)
       !!{
       Interface for list2D names.
       !!}
       import varying_string, nodePropertyExtractorList2D
       class(nodePropertyExtractorList2D), intent(inout)                             :: self
       type (varying_string             ), intent(inout), dimension(:) , allocatable :: names
     end subroutine list2DNames
  end interface

  abstract interface
     subroutine list2DDescriptions(self,descriptions)
       !!{
       Interface for list2D descriptions.
       !!}
       import varying_string, nodePropertyExtractorList2D
       class(nodePropertyExtractorList2D), intent(inout)                             :: self
       type (varying_string             ), intent(inout), dimension(:) , allocatable :: descriptions
     end subroutine list2DDescriptions
  end interface

  abstract interface
     function list2DUnitsInSI(self)
       !!{
       Interface for list2D property units.
       !!}
       import nodePropertyExtractorList2D
       double precision                             , dimension(:) , allocatable :: list2DUnitsInSI
       class           (nodePropertyExtractorList2D), intent(inout)              :: self
     end function list2DUnitsInSI
  end interface

contains
  
  subroutine list2DMetaData(self,node,time,metaDataRank0,metaDataRank1)
    !!{
    Interface for list2D property meta-data.
    !!}
    implicit none
    class           (nodePropertyExtractorList2D), intent(inout) :: self
    type            (treeNode                   ), intent(inout) :: node
    double precision                             , intent(in   ) :: time
    type            (doubleHash                 ), intent(inout) :: metaDataRank0
    type            (rank1DoubleHash            ), intent(inout) :: metaDataRank1
    !$GLC attributes unused :: self, node, time, metaDataRank0, metaDataRank1
    
    return
  end subroutine list2DMetaData
