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

  use :: Kind_Numbers, only : kind_int8
  use :: Hashes      , only : doubleHash, rank1DoubleHash

  !![
  <nodePropertyExtractor name="nodePropertyExtractorIntegerList" abstract="yes">
   <description>An abstract output analysis property extractor class which provides a list of integer properties.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorClass), abstract :: nodePropertyExtractorIntegerList
     !!{
     An integer list property extractor.
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
     procedure(integerListElementCount), deferred :: elementCount
     procedure(integerListExtract     ), deferred :: extract
     procedure(integerListNames       ), deferred :: names
     procedure(integerListDescriptions), deferred :: descriptions
     procedure(integerListUnitsInSI   ), deferred :: unitsInSI
     procedure                                    :: metaData     => integerListMetaData
  end type nodePropertyExtractorIntegerList

  abstract interface
     function integerListElementCount(self)
       !!{
       Interface for integer list property count.
       !!}
       import nodePropertyExtractorIntegerList
       integer                                                  :: integerListElementCount
       class  (nodePropertyExtractorIntegerList), intent(inout) :: self
     end function integerListElementCount
  end interface

  abstract interface
     function integerListExtract(self,node,instance)
       !!{
       Interface for integer list property extraction.
       !!}
       import nodePropertyExtractorIntegerList, treeNode, multiCounter, kind_int8
       integer(kind_int8)                       , dimension(:,:), allocatable :: integerListExtract
       class  (nodePropertyExtractorIntegerList), intent(inout)               :: self
       type   (treeNode                        ), intent(inout)               :: node
       type   (multiCounter                    ), intent(inout) , optional    :: instance
     end function integerListExtract
  end interface

  abstract interface
     subroutine integerListNames(self,names)
       !!{
       Interface for list names.
       !!}
       import varying_string, nodePropertyExtractorIntegerList
       class(nodePropertyExtractorIntegerList), intent(inout)                             :: self
       type (varying_string                  ), intent(inout), dimension(:) , allocatable :: names
     end subroutine integerListNames
  end interface

  abstract interface
     subroutine integerListDescriptions(self,descriptions)
       !!{
       Interface for list descriptions.
       !!}
       import varying_string, nodePropertyExtractorIntegerList
       class(nodePropertyExtractorIntegerList), intent(inout)                             :: self
       type (varying_string                  ), intent(inout), dimension(:) , allocatable :: descriptions
     end subroutine integerListDescriptions
  end interface

  abstract interface
     function integerListUnitsInSI(self)
       !!{
       Interface for list property units.
       !!}
       import nodePropertyExtractorIntegerList
       double precision                                  , dimension(:) , allocatable :: integerListUnitsInSI
       class           (nodePropertyExtractorIntegerList), intent(inout)              :: self
     end function integerListUnitsInSI
  end interface

contains
  
  subroutine integerListMetaData(self,node,metaDataRank0,metaDataRank1)
    !!{
    Interface for list property meta-data.
    !!}
    implicit none
    class(nodePropertyExtractorIntegerList), intent(inout) :: self
    type (treeNode                        ), intent(inout) :: node
    type (doubleHash                      ), intent(inout) :: metaDataRank0
    type (rank1DoubleHash                 ), intent(inout) :: metaDataRank1
    !$GLC attributes unused :: self, node, metaDataRank0, metaDataRank1
    
    return
  end subroutine integerListMetaData
