!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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

  use :: Hashes, only : doubleHash

  !![
  <nodePropertyExtractor name="nodePropertyExtractorList" abstract="yes">
   <description>An abstract output analysis property extractor class which provieds a list of floating point properties.</description>
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
       <method method="extract"     description="Extract the properties from the given {\normalfont \ttfamily node}."/>
       <method method="name"        description="Return the name of the properties extracted."                       />
       <method method="description" description="Return a description of the properties extracted."                  />
       <method method="unitsInSI"   description="Return the units of the properties extracted in the SI system."     />
       <method method="metaData"    description="Populate a hash with meta-data for the property."                   />
     </methods>
     !!]
     procedure(listExtract    ), deferred :: extract
     procedure(listName       ), deferred :: name
     procedure(listDescription), deferred :: description
     procedure(listUnitsInSI  ), deferred :: unitsInSI
     procedure                            :: metaData    => listMetaData
  end type nodePropertyExtractorList

  abstract interface
     function listExtract(self,node,instance)
       !!{
       Interface for list property extraction.
       !!}
       import nodePropertyExtractorList, treeNode, multiCounter
       double precision                           , dimension(:) , allocatable :: listExtract
       class           (nodePropertyExtractorList), intent(inout)              :: self
       type            (treeNode                 ), intent(inout)              :: node
       type            (multiCounter             ), intent(inout), optional    :: instance
     end function listExtract
  end interface

  abstract interface
     function listName(self)
       !!{
       Interface for list names.
       !!}
       import varying_string, nodePropertyExtractorList
       type (varying_string           )                :: listName
       class(nodePropertyExtractorList), intent(inout) :: self
     end function listName
  end interface

  abstract interface
     function listDescription(self)
       !!{
       Interface for list descriptions.
       !!}
       import varying_string, nodePropertyExtractorList
       type (varying_string           )                :: listDescription
       class(nodePropertyExtractorList), intent(inout) :: self
     end function listDescription
  end interface

  abstract interface
     double precision function listUnitsInSI(self)
       !!{
       Interface for list property units.
       !!}
       import nodePropertyExtractorList
       class (nodePropertyExtractorList), intent(inout) :: self
     end function listUnitsInSI
  end interface

contains
  
  subroutine listMetaData(self,metaData)
    !!{
    Interface for list property meta-data.
    !!}
    implicit none
    class(nodePropertyExtractorList), intent(inout) :: self
    type (doubleHash               ), intent(inout) :: metaData
    !$GLC attributes unused :: self, metaData
    
    return
  end subroutine listMetaData
