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
  <nodePropertyExtractor name="nodePropertyExtractorPresetNamedReals">
   <description>
    A node property extractor which extracts ``preset'' named real quantities. These are typically used to provide additional
    quantities read from N-body merger trees.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorPresetNamedReals
     !!{
     A property extractor which extracts presetNamedReals properties.
     !!}
     private
     type   (varying_string), allocatable, dimension(:) :: presetNames
     integer                , allocatable, dimension(:) :: indexPresetNames
   contains
     procedure :: elementCount => presetNamedRealsElementCount
     procedure :: extract      => presetNamedRealsExtract
     procedure :: names        => presetNamedRealsNames
     procedure :: descriptions => presetNamedRealsDescriptions
     procedure :: unitsInSI    => presetNamedRealsUnitsInSI
  end type nodePropertyExtractorPresetNamedReals

  interface nodePropertyExtractorPresetNamedReals
     !!{
     Constructors for the \refClass{nodePropertyExtractorPresetNamedReals} output extractor class.
     !!}
     module procedure presetNamedRealsConstructorParameters
     module procedure presetNamedRealsConstructorInternal
  end interface nodePropertyExtractorPresetNamedReals

contains

  function presetNamedRealsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorPresetNamedReals} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorPresetNamedReals)                              :: self
    type(inputParameters                      ), intent(inout)               :: parameters
    type(varying_string                       ), allocatable  , dimension(:) :: presetNames

    allocate(presetNames(parameters%count('presetNames')))
    !![
    <inputParameter>
      <name>presetNames</name>
      <description>The names of preset properties to extract.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=nodePropertyExtractorPresetNamedReals(presetNames)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function presetNamedRealsConstructorParameters

  function presetNamedRealsConstructorInternal(presetnames) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorPresetNamedReals} output extractor property extractor class.
    !!}
    implicit none
    type   (nodePropertyExtractorPresetNamedReals)                              :: self
    type   (varying_string                       ), intent(in   ), dimension(:) :: presetNames
    integer                                                                     :: i
    !![
    <constructorAssign variables="presetNames"/>
    !!]

    allocate(self%indexPresetNames(size(self%presetNames)))
    do i=1,size(self%presetNames)
       !![
       <addMetaProperty component="basic" name="'preset:'//char(self%presetNames(i))" type="float" id="self%indexPresetNames(i)" isCreator="no"/>
       !!]
    end do
    return
  end function presetNamedRealsConstructorInternal

  integer function presetNamedRealsElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily presetNamedReals} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorPresetNamedReals), intent(inout) :: self
    double precision                                       , intent(in   ) :: time
    !$GLC attributes unused :: time

    presetNamedRealsElementCount=size(self%presetNames)
    return
  end function presetNamedRealsElementCount

  function presetNamedRealsExtract(self,node,time,instance)
    !!{
    Implement a presetNamedReals output extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    double precision                                       , dimension(:) , allocatable :: presetNamedRealsExtract
    class           (nodePropertyExtractorPresetNamedReals), intent(inout), target      :: self
    type            (treeNode                             ), intent(inout), target      :: node
    double precision                                       , intent(in   )              :: time
    type            (multiCounter                         ), intent(inout), optional    :: instance
    class           (nodeComponentBasic                   )               , pointer     :: basic
    integer                                                                             :: i
    !$GLC attributes unused :: time, instance

    basic => node%basic()
    allocate(presetNamedRealsExtract(size(self%presetNames)))
    do i=1,size(self%presetNames)
       presetNamedRealsExtract(i)=basic%floatRank0MetaPropertyGet(self%indexPresetNames(i))
    end do
    return
  end function presetNamedRealsExtract

  subroutine presetNamedRealsNames(self,time,names)
    !!{
    Return the names of the {\normalfont \ttfamily presetNamedReals} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorPresetNamedReals), intent(inout)                             :: self
    double precision                                       , intent(in   )                             :: time
    type            (varying_string                       ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: time

    allocate(names(size(self%presetNames)))
    names='preset:'//self%presetNames
    return
  end subroutine presetNamedRealsNames

  subroutine presetNamedRealsDescriptions(self,time,descriptions)
    !!{
    Return the descriptions of the {\normalfont \ttfamily presetNamedReals} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorPresetNamedReals), intent(inout)                             :: self
    double precision                                       , intent(in   )                             :: time
    type            (varying_string                       ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: time

    allocate(descriptions(size(self%presetNames)))
    descriptions=var_str('Preset named property.')
    return
  end subroutine presetNamedRealsDescriptions

  function presetNamedRealsUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily presetNamedReals} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    double precision                                       , dimension(:) , allocatable :: presetNamedRealsUnitsInSI
    class           (nodePropertyExtractorPresetNamedReals), intent(inout)              :: self
    double precision                                       , intent(in   )              :: time
   !$GLC attributes unused :: time

    allocate(presetNamedRealsUnitsInSI(size(self%presetNames)))
    presetNamedRealsUnitsInSI=-1.0d0
    return
  end function presetNamedRealsUnitsInSI

