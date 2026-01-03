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
  <nodePropertyExtractor name="nodePropertyExtractorTuple" abstract="yes">
   <description>An abstract output analysis property extractor class which provides a tuple of floating point properties.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorClass), abstract :: nodePropertyExtractorTuple
     !!{
     A tuple property extractor.
     !!}
     private
   contains
     !![
     <methods>
       <method method="elementCount" description="Return the number of properties in the tuple."                      />
       <method method="extract"      description="Extract the properties from the given {\normalfont \ttfamily node}."/>
       <method method="names"        description="Return the names of the properties extracted."                      />
       <method method="descriptions" description="Return descriptions of the properties extracted."                   />
       <method method="unitsInSI"    description="Return the units of the properties extracted in the SI system."     />
       <method method="metaData"     description="Populate a hash with meta-data for the property."                   />
     </methods>
     !!]
     procedure(tupleElementCount), deferred :: elementCount
     procedure(tupleExtract     ), deferred :: extract
     procedure(tupleNames       ), deferred :: names
     procedure(tupleDescriptions), deferred :: descriptions
     procedure(tupleUnitsInSI   ), deferred :: unitsInSI
     procedure                              :: metaData     => tupleMetaData
  end type nodePropertyExtractorTuple

  abstract interface
     function tupleExtract(self,node,time,instance)
       !!{
       Interface for tuple property extraction.
       !!}
       import nodePropertyExtractorTuple, treeNode, multiCounter
       double precision                            , dimension(:) , allocatable :: tupleExtract
       class           (nodePropertyExtractorTuple), intent(inout), target      :: self
       type            (treeNode                  ), intent(inout), target      :: node
       double precision                            , intent(in   )              :: time
       type            (multiCounter              ), intent(inout), optional    :: instance
     end function tupleExtract
  end interface

  abstract interface
     subroutine tupleNames(self,time,names)
       !!{
       Interface for tuple property names.
       !!}
       import varying_string, nodePropertyExtractorTuple
       class           (nodePropertyExtractorTuple), intent(inout)                             :: self
       double precision                            , intent(in   )                             :: time
       type            (varying_string            ), intent(inout), dimension(:) , allocatable :: names
     end subroutine tupleNames
  end interface

  abstract interface
     subroutine tupleDescriptions(self,time,descriptions)
       !!{
       Interface for tuple property names.
       !!}
       import varying_string, nodePropertyExtractorTuple
       class           (nodePropertyExtractorTuple), intent(inout)                             :: self
       double precision                            , intent(in   )                             :: time
       type            (varying_string            ), intent(inout), dimension(:) , allocatable :: descriptions
     end subroutine tupleDescriptions
  end interface

  abstract interface
     function tupleUnitsInSI(self,time)
       !!{
       Interface for tuple property units.
       !!}
       import nodePropertyExtractorTuple
       double precision                            , dimension(:) , allocatable :: tupleUnitsInSI
       class           (nodePropertyExtractorTuple), intent(inout)              :: self
       double precision                            , intent(in   )              :: time
     end function tupleUnitsInSI
  end interface

  abstract interface
     integer function tupleElementCount(self,time)
       !!{
       Interface for tuple element count.
       !!}
       import nodePropertyExtractorTuple
       class           (nodePropertyExtractorTuple), intent(inout) :: self
       double precision                            , intent(in   ) :: time
     end function tupleElementCount
  end interface

contains
  
  subroutine tupleMetaData(self,node,indexProperty,metaDataRank0,metaDataRank1)
    !!{
    Interface for tuple property meta-data.
    !!}
    implicit none
    class  (nodePropertyExtractorTuple), intent(inout) :: self
    type   (treeNode                  ), intent(inout) :: node
    integer                            , intent(in   ) :: indexProperty
    type   (doubleHash                ), intent(inout) :: metaDataRank0
    type   (rank1DoubleHash           ), intent(inout) :: metaDataRank1
    !$GLC attributes unused :: self, node, indexProperty, metaDataRank0, metaDataRank1
    
    return
  end subroutine tupleMetaData
