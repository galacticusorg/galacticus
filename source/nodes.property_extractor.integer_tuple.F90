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
  <nodePropertyExtractor name="nodePropertyExtractorIntegerTuple" abstract="yes">
   <description>An abstract output analysis property extractor class which provides a tuple of integer properties.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorClass), abstract :: nodePropertyExtractorIntegerTuple
     !!{
     A integerTuple property extractor.
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
     procedure(integerTupleElementCount), deferred :: elementCount
     procedure(integerTupleExtract     ), deferred :: extract
     procedure(integerTupleNames       ), deferred :: names
     procedure(integerTupleDescriptions), deferred :: descriptions
     procedure(integerTupleUnitsInSI   ), deferred :: unitsInSI
     procedure                                     :: metaData    => integerTupleMetaData
  end type nodePropertyExtractorIntegerTuple

  abstract interface
     function integerTupleExtract(self,node,time,instance)
       !!{
       Interface for {\normalfont \ttfamily integerTuple} property extraction.
       !!}
       import nodePropertyExtractorIntegerTuple, treeNode, multiCounter, kind_int8
       integer         (kind_int8                        ), dimension(:) , allocatable :: integerTupleExtract
       class           (nodePropertyExtractorIntegerTuple), intent(inout)              :: self
       type            (treeNode                         ), intent(inout)              :: node
       double precision                                   , intent(in   )              :: time
       type            (multiCounter                     ), intent(inout), optional    :: instance
     end function integerTupleExtract
  end interface

  abstract interface
     subroutine integerTupleNames(self,time,names)
       !!{
       Interface for {\normalfont \ttfamily integerTuple} property names.
       !!}
       import varying_string, nodePropertyExtractorIntegerTuple
       class           (nodePropertyExtractorIntegerTuple), intent(inout)                             :: self
       double precision                                   , intent(in   )                             :: time
       type            (varying_string                   ), intent(inout), dimension(:) , allocatable :: names
     end subroutine integerTupleNames
  end interface

  abstract interface
     subroutine integerTupleDescriptions(self,time,descriptions)
       !!{
       Interface for {\normalfont \ttfamily integerTuple} property descriptions.
       !!}
       import varying_string, nodePropertyExtractorIntegerTuple
       class           (nodePropertyExtractorIntegerTuple), intent(inout)                             :: self
       double precision                                   , intent(in   )                             :: time
       type            (varying_string                   ), intent(inout), dimension(:) , allocatable :: descriptions
     end subroutine integerTupleDescriptions
  end interface

  abstract interface
     function integerTupleUnitsInSI(self,time)
       !!{
       Interface for {\normalfont \ttfamily integerTuple property units.
       !!}
       import nodePropertyExtractorIntegerTuple
       double precision                                   , dimension(:) , allocatable :: integerTupleUnitsInSI
       class           (nodePropertyExtractorIntegerTuple), intent(inout)              :: self
       double precision                                   , intent(in   )              :: time
     end function integerTupleUnitsInSI
  end interface

  abstract interface
     integer function integerTupleElementCount(self,time)
       !!{
       Interface for {\normalfont \ttfamily integerTuple} element count.
       !!}
       import nodePropertyExtractorIntegerTuple
       class           (nodePropertyExtractorIntegerTuple), intent(inout) :: self
       double precision                                   , intent(in   ) :: time
     end function integerTupleElementCount
  end interface

contains

  subroutine integerTupleMetaData(self,node,indexProperty,metaDataRank0,metaDataRank1)
    !!{
    Interface for integerTuple property meta-data.
    !!}
    implicit none
    class  (nodePropertyExtractorIntegerTuple), intent(inout) :: self
    type   (treeNode                         ), intent(inout) :: node
    integer                                   , intent(in   ) :: indexProperty
    type   (doubleHash                       ), intent(inout) :: metaDataRank0
    type   (rank1DoubleHash                  ), intent(inout) :: metaDataRank1
    !$GLC attributes unused :: self, node, indexProperty, metaDataRank0, metaDataRank1
    
    return
  end subroutine integerTupleMetaData
