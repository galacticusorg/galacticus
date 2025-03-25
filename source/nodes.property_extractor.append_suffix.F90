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
  Implements a node property extractor that appends a suffix to property names.
  !!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorAppendSuffix">
   <description>A node property extractor that appends a suffix to property names.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorMulti) :: nodePropertyExtractorAppendSuffix
     !!{
     A node property extractor that appends a suffix to property names.
     !!}
     private
     type(varying_string) :: suffix
   contains
     procedure :: names => appendSuffixNames
  end type nodePropertyExtractorAppendSuffix

  interface nodePropertyExtractorAppendSuffix
     !!{
     Constructors for the {\normalfont \ttfamily appendSuffix} output extractor class.
     !!}
     module procedure appendSuffixConstructorParameters
     module procedure appendSuffixConstructorInternal
  end interface nodePropertyExtractorAppendSuffix

contains

  function appendSuffixConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily appendSuffix} output extractor property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(nodePropertyExtractorAppendSuffix)                :: self
    type(inputParameters                  ), intent(inout) :: parameters

    self%nodePropertyExtractorMulti=nodePropertyExtractorMulti(parameters)
    !![
    <inputParameter>
      <name>suffix</name>
      <variable>self%suffix</variable>
      <description>The suffix to append to parameter names.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParametersValidate source="parameters" multiParameters="nodePropertyExtractor"/>
    !!]
    return
  end function appendSuffixConstructorParameters

  function appendSuffixConstructorInternal(suffix,extractors) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily appendSuffix} output extractor property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorAppendSuffix)                :: self
    type(multiExtractorList               ), intent(in   ) :: extractors
    type(varying_string                   ), intent(in   ) :: suffix
    !![
    <constructorAssign variables="suffix"/>
    !!]

    self%nodePropertyExtractorMulti=nodePropertyExtractorMulti(extractors)
    return
  end function appendSuffixConstructorInternal

  subroutine appendSuffixNames(self,elementType,time,names)
    !!{
    Return the names of the suffixed properties.
    !!}
    implicit none
    class           (nodePropertyExtractorAppendSuffix), intent(inout)                             :: self
    type            (enumerationElementTypeType       ), intent(in   )                             :: elementType
    double precision                                   , intent(in   )                             :: time
    type            (varying_string                   ), intent(inout), dimension(:) , allocatable :: names
    integer                                                                                        :: i
    
    call self%nodePropertyExtractorMulti%names(elementType,time,names)
    do i=1,size(names)
       names(i)=names(i)//self%suffix
    end do
    return
  end subroutine appendSuffixNames
