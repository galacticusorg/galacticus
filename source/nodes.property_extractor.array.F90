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
  <nodePropertyExtractor name="nodePropertyExtractorArray" abstract="yes">
   <description>An abstract output analysis property extractor class which provides a array of floating point properties.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorClass), abstract :: nodePropertyExtractorArray
     !!{
     A array property extractor.
     !!}
     private
   contains
     !![
     <methods>
       <method method="columnDescriptions" description="Return a description of the columns."                               />
       <method method="size"               description="Return the number of elements in the array."                        />
       <method method="elementCount"       description="Return the number of properties in the array."                      />
       <method method="extract"            description="Extract the properties from the given {\normalfont \ttfamily node}."/>
       <method method="names"              description="Return the name of the properties extracted."                       />
       <method method="descriptions"       description="Return a description of the properties extracted."                  />
       <method method="unitsInSI"          description="Return the units of the properties extracted in the SI system."     />
       <method method="metaData"           description="Populate a hash with meta-data for the property."                   />
     </methods>
     !!]
     procedure(arrayColumns     ), deferred :: columnDescriptions
     procedure(arraySize        ), deferred :: size
     procedure(arrayElementCount), deferred :: elementCount
     procedure(arrayExtract     ), deferred :: extract
     procedure(arrayNames       ), deferred :: names
     procedure(arrayDescriptions), deferred :: descriptions
     procedure(arrayUnitsInSI   ), deferred :: unitsInSI
     procedure                              :: metaData           => arrayMetaData
  end type nodePropertyExtractorArray

  abstract interface
     function arrayExtract(self,node,time,instance)
       !!{
       Interface for array property extraction.
       !!}
       import nodePropertyExtractorArray, treeNode, multiCounter
       double precision                            , dimension(:,:), allocatable :: arrayExtract
       class           (nodePropertyExtractorArray), intent(inout) , target      :: self
       type            (treeNode                  ), intent(inout) , target      :: node
       double precision                            , intent(in   )               :: time
       type            (multiCounter              ), intent(inout) , optional    :: instance
     end function arrayExtract
  end interface

  abstract interface
     subroutine arrayNames(self,names,time)
       !!{
       Interface for array names.
       !!}
       import varying_string, nodePropertyExtractorArray
       class           (nodePropertyExtractorArray), intent(inout)                            :: self
       double precision                            , intent(in   ), optional                  :: time
       type            (varying_string            ), intent(inout), allocatable, dimension(:) :: names
    end subroutine arrayNames
  end interface

  abstract interface
     subroutine arrayDescriptions(self,descriptions,time)
       !!{
       Interface for array descriptions.
       !!}
       import varying_string, nodePropertyExtractorArray
       class           (nodePropertyExtractorArray), intent(inout)                            :: self
       double precision                            , intent(in   ), optional                  :: time
       type            (varying_string            ), intent(inout), allocatable, dimension(:) :: descriptions
    end subroutine arrayDescriptions
  end interface

  abstract interface
     subroutine arrayColumns(self,descriptions,values,valuesDescription,valuesUnitsInSI,time)
       !!{
       Interface for array column descriptions.
       !!}
       import varying_string, nodePropertyExtractorArray
       class           (nodePropertyExtractorArray), intent(inout)                            :: self
       double precision                            , intent(in   ), optional                  :: time
       type            (varying_string            ), intent(inout), allocatable, dimension(:) :: descriptions
       double precision                            , intent(inout), allocatable, dimension(:) :: values
       type            (varying_string            ), intent(  out)                            :: valuesDescription
       double precision                            , intent(  out)                            :: valuesUnitsInSI
     end subroutine arrayColumns
  end interface

  abstract interface
     function arrayUnitsInSI(self,time)
       !!{
       Interface for array property units.
       !!}
       import nodePropertyExtractorArray
       double precision                            , allocatable  , dimension(:) :: arrayUnitsInSI
       class           (nodePropertyExtractorArray), intent(inout)               :: self
       double precision                            , intent(in   ), optional     :: time
     end function arrayUnitsInSI
  end interface

  abstract interface
     integer function arrayElementCount(self,time)
       !!{
       Interface for array element count.
       !!}
       import nodePropertyExtractorArray
       class           (nodePropertyExtractorArray), intent(inout) :: self
       double precision                            , intent(in   ) :: time
     end function arrayElementCount
  end interface

  abstract interface
     function arraySize(self,time)
       !!{
       Interface for array element count.
       !!}
       import nodePropertyExtractorArray, c_size_t
       integer         (c_size_t                  )                :: arraySize
       class           (nodePropertyExtractorArray), intent(inout) :: self
       double precision                            , intent(in   ) :: time
     end function arraySize
  end interface

contains
  
  subroutine arrayMetaData(self,node,indexProperty,metaDataRank0,metaDataRank1)
    !!{
    Interface for array property meta-data.
    !!}
    implicit none
    class  (nodePropertyExtractorArray), intent(inout) :: self
    type   (treeNode                  ), intent(inout) :: node
    integer                            , intent(in   ) :: indexProperty
    type   (doubleHash                ), intent(inout) :: metaDataRank0
    type   (rank1DoubleHash           ), intent(inout) :: metaDataRank1
    !$GLC attributes unused :: self, node, indexProperty, metaDataRank0, metaDataRank1
    
    return
  end subroutine arrayMetaData
