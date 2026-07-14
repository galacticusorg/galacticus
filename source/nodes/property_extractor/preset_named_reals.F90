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
  <nodePropertyExtractor name="nodePropertyExtractorPresetNamedReals" docformat="rst">
   <description>
   A node property extractor which extracts "preset" named real quantities. These are typically used to provide additional quantities read from N-body merger trees.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorPresetNamedReals
     !!{RST
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
     procedure :: units       => presetNamedRealsUnits
  end type nodePropertyExtractorPresetNamedReals

  interface nodePropertyExtractorPresetNamedReals
     !!{RST
     Constructors for the :galacticus-class:`nodePropertyExtractorPresetNamedReals` property extractor class.
     !!}
     module procedure presetNamedRealsConstructorParameters
     module procedure presetNamedRealsConstructorInternal
  end interface nodePropertyExtractorPresetNamedReals

contains

  function presetNamedRealsConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodePropertyExtractorPresetNamedReals` property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorPresetNamedReals)                              :: self
    type(inputParameters                      ), intent(inout)               :: parameters
    type(varying_string                       ), allocatable  , dimension(:) :: presetNames

    allocate(presetNames(parameters%count('presetNames')))
    !![
    <inputParameter docformat="rst">
      <name>presetNames</name>
      <description>
      The names of preset properties to extract.
      </description>
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
    !!{RST
    Internal constructor for the :galacticus-class:`nodePropertyExtractorPresetNamedReals` property extractor class.
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
    !!{RST
    Return the number of elements in the ``presetNamedReals`` property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorPresetNamedReals), intent(inout) :: self
    double precision                                       , intent(in   ) :: time
    !$GLC attributes unused :: time

    presetNamedRealsElementCount=size(self%presetNames)
    return
  end function presetNamedRealsElementCount

  function presetNamedRealsExtract(self,node,time,instance)
    !!{RST
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
    !!{RST
    Return the names of the ``presetNamedReals`` properties.
    !!}
    implicit none
    class           (nodePropertyExtractorPresetNamedReals), intent(inout)                             :: self
    double precision                                       , intent(in   )                             :: time
    type            (varying_string                       ), intent(inout), dimension(:) , allocatable :: names
    integer                                                                                            :: i
    !$GLC attributes unused :: time

    allocate(names(size(self%presetNames)))
    ! Note: an explicit loop is used here (rather than a whole-array elemental assignment) as the elemental concatenation is
    ! miscompiled under link-time optimization (observed with gfortran 16), resulting in a double-free of the temporary - see
    ! https://github.com/galacticusorg/galacticus/issues/1216.
    do i=1,size(self%presetNames)
       names(i)='preset:'//self%presetNames(i)
    end do
    return
  end subroutine presetNamedRealsNames

  subroutine presetNamedRealsDescriptions(self,time,descriptions)
    !!{RST
    Return the descriptions of the ``presetNamedReals`` properties.
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
    !!{RST
    Return the units of the ``presetNamedReals`` properties in the SI system.
    !!}
    implicit none
    double precision                                       , dimension(:) , allocatable :: presetNamedRealsUnitsInSI
    class           (nodePropertyExtractorPresetNamedReals), intent(inout)              :: self
    double precision                                       , intent(in   )              :: time
    !$GLC attributes unused :: time

    allocate(presetNamedRealsUnitsInSI(size(self%presetNames)))
    presetNamedRealsUnitsInSI=-1.0d0
    return
  end function presetNamedRealsUnitsInSI

  function presetNamedRealsUnits(self,time) result(units)
    !!{RST
    Return the units of the presetNamedReals properties.
    !!}
    use :: Units_MetaData, only : unitType
    implicit none
    type            (unitType                             ), dimension(:) , allocatable :: units
    class           (nodePropertyExtractorPresetNamedReals), intent(inout)              :: self
    double precision                                       , intent(in   )              :: time
    double precision                                       , dimension(:) , allocatable :: siValues
    integer                                                                             :: i

    siValues=self%unitsInSI(time)
    allocate(units(size(siValues)))
    do i=1,size(siValues)
       units(i)=unitType(-1.0d0)
    end do
    return
  end function presetNamedRealsUnits
