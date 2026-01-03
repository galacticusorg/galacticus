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

  !![
  <nodePropertyExtractor name="nodePropertyExtractorPresetNamedIntegers">
   <description>
    A node property extractor which extracts ``preset'' named integer quantities. These are typically used to provide additional
    quantities read from N-body merger trees.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorIntegerTuple) :: nodePropertyExtractorPresetNamedIntegers
     !!{
     A property extractor which extracts presetNamedIntegers properties.
     !!}
     private
     type   (varying_string), allocatable, dimension(:) :: presetNames
     integer                , allocatable, dimension(:) :: indexPresetNames
   contains
     procedure :: elementCount => presetNamedIntegersElementCount
     procedure :: extract      => presetNamedIntegersExtract
     procedure :: names        => presetNamedIntegersNames
     procedure :: descriptions => presetNamedIntegersDescriptions
     procedure :: unitsInSI    => presetNamedIntegersUnitsInSI
  end type nodePropertyExtractorPresetNamedIntegers

  interface nodePropertyExtractorPresetNamedIntegers
     !!{
     Constructors for the \refClass{nodePropertyExtractorPresetNamedIntegers} output extractor class.
     !!}
     module procedure presetNamedIntegersConstructorParameters
     module procedure presetNamedIntegersConstructorInternal
  end interface nodePropertyExtractorPresetNamedIntegers

contains

  function presetNamedIntegersConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorPresetNamedIntegers} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorPresetNamedIntegers)                              :: self
    type(inputParameters                         ), intent(inout)               :: parameters
    type(varying_string                          ), allocatable  , dimension(:) :: presetNames

    allocate(presetNames(parameters%count('presetNames')))
    !![
    <inputParameter>
      <name>presetNames</name>
      <description>The names of preset properties to extract.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=nodePropertyExtractorPresetNamedIntegers(presetNames)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function presetNamedIntegersConstructorParameters

  function presetNamedIntegersConstructorInternal(presetnames) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorPresetNamedIntegers} output extractor property extractor class.
    !!}
    implicit none
    type   (nodePropertyExtractorPresetNamedIntegers)                              :: self
    type   (varying_string                          ), intent(in   ), dimension(:) :: presetNames
    integer                                                                        :: i
    !![
    <constructorAssign variables="presetNames"/>
    !!]

    allocate(self%indexPresetNames(size(self%presetNames)))
    do i=1,size(self%presetNames)
       !![
       <addMetaProperty component="basic" name="'preset:'//char(self%presetNames(i))" type="longInteger" id="self%indexPresetNames(i)" isCreator="no"/>
       !!]
    end do
    return
  end function presetNamedIntegersConstructorInternal

  integer function presetNamedIntegersElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily presetNamedIntegers} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorPresetNamedIntegers), intent(inout) :: self
    double precision                                          , intent(in   ) :: time
    !$GLC attributes unused :: time

    presetNamedIntegersElementCount=size(self%presetNames)
    return
  end function presetNamedIntegersElementCount

  function presetNamedIntegersExtract(self,node,time,instance)
    !!{
    Implement a presetNamedIntegers output extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    integer         (kind_int8                               ), dimension(:) , allocatable :: presetNamedIntegersExtract
    class           (nodePropertyExtractorPresetNamedIntegers), intent(inout)              :: self
    type            (treeNode                                ), intent(inout)              :: node
    double precision                                          , intent(in   )              :: time
    type            (multiCounter                            ), intent(inout), optional    :: instance
    class           (nodeComponentBasic                      )               , pointer     :: basic
    integer                                                                                :: i
    !$GLC attributes unused :: time, instance

    basic => node%basic()
    allocate(presetNamedIntegersExtract(size(self%presetNames)))
    do i=1,size(self%presetNames)
       presetNamedIntegersExtract(i)=basic%longIntegerRank0MetaPropertyGet(self%indexPresetNames(i))
    end do
    return
  end function presetNamedIntegersExtract

  subroutine presetNamedIntegersNames(self,time,names)
    !!{
    Return the names of the {\normalfont \ttfamily presetNamedIntegers} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorPresetNamedIntegers), intent(inout)                             :: self
    double precision                                          , intent(in   )                             :: time
    type            (varying_string                          ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: time

    allocate(names(size(self%presetNames)))
    names='preset:'//self%presetNames
    return
  end subroutine presetNamedIntegersNames

  subroutine presetNamedIntegersDescriptions(self,time,descriptions)
    !!{
    Return the descriptions of the {\normalfont \ttfamily presetNamedIntegers} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorPresetNamedIntegers), intent(inout)                             :: self
    double precision                                          , intent(in   )                             :: time
    type            (varying_string                          ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: time

    allocate(descriptions(size(self%presetNames)))
    descriptions=var_str('Preset named property.')
    return
  end subroutine presetNamedIntegersDescriptions

  function presetNamedIntegersUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily presetNamedIntegers} properties in the SI system.
    !!}
    implicit none
    double precision                                          , dimension(:) , allocatable :: presetNamedIntegersUnitsInSI
    class           (nodePropertyExtractorPresetNamedIntegers), intent(inout)              :: self
    double precision                                          , intent(in   )              :: time
   !$GLC attributes unused :: time

    allocate(presetNamedIntegersUnitsInSI(size(self%presetNames)))
    presetNamedIntegersUnitsInSI=-1.0d0
    return
  end function presetNamedIntegersUnitsInSI

