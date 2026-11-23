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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Implements a node property extractor for the index, in the catalog from which it was read, of a postprocessed galaxy.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorIndexPostprocessSource" docformat="rst">
   <description>
   Extracts the index (counting from zero) which a galaxy had in the ``nodeData`` datasets of the output file from which it was
   read by the :galacticus-class:`taskPostprocessOutputs` task. Postprocessed galaxies are not guaranteed to be output in the
   order in which they were read, and may be a subset of those in the original file, so this index is what allows each to be
   matched to its entry in the original file.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorIntegerScalar) :: nodePropertyExtractorIndexPostprocessSource
     !!{RST
     A property extractor for the index of a postprocessed galaxy in the catalog from which it was read.
     !!}
     private
     integer :: indexPostprocessSourceID
   contains
     procedure :: extract     => indexPostprocessSourceExtract
     procedure :: name        => indexPostprocessSourceName
     procedure :: description => indexPostprocessSourceDescription
  end type nodePropertyExtractorIndexPostprocessSource

  interface nodePropertyExtractorIndexPostprocessSource
     !!{RST
     Constructors for the :galacticus-class:`nodePropertyExtractorIndexPostprocessSource` property extractor class.
     !!}
     module procedure indexPostprocessSourceConstructorParameters
     module procedure indexPostprocessSourceConstructorInternal
  end interface nodePropertyExtractorIndexPostprocessSource

contains

  function indexPostprocessSourceConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodePropertyExtractorIndexPostprocessSource` property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorIndexPostprocessSource)                :: self
    type(inputParameters                            ), intent(inout) :: parameters

    self=nodePropertyExtractorIndexPostprocessSource()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function indexPostprocessSourceConstructorParameters

  function indexPostprocessSourceConstructorInternal() result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodePropertyExtractorIndexPostprocessSource` property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorIndexPostprocessSource) :: self

    !![
    <addMetaProperty component="basic" name="nodeIndexPostprocessSource" type="longInteger" id="self%indexPostprocessSourceID" isCreator="no"/>
    !!]
    return
  end function indexPostprocessSourceConstructorInternal

  function indexPostprocessSourceExtract(self,node,time,instance)
    !!{RST
    Implement an ``indexPostprocessSource`` node property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    integer         (kind_int8                                  )                          :: indexPostprocessSourceExtract
    class           (nodePropertyExtractorIndexPostprocessSource), intent(inout)           :: self
    type            (treeNode                                   ), intent(inout), target   :: node
    double precision                                             , intent(in   )           :: time
    type            (multiCounter                               ), intent(inout), optional :: instance
    class           (nodeComponentBasic                         )               , pointer  :: basic
    !$GLC attributes unused :: instance, time

    basic                         => node %basic                          (                             )
    indexPostprocessSourceExtract =  basic%longIntegerRank0MetaPropertyGet(self%indexPostprocessSourceID)
    return
  end function indexPostprocessSourceExtract

  function indexPostprocessSourceName(self)
    !!{RST
    Return the name of the source index property.
    !!}
    implicit none
    type (varying_string                             )                :: indexPostprocessSourceName
    class(nodePropertyExtractorIndexPostprocessSource), intent(inout) :: self
    !$GLC attributes unused :: self

    indexPostprocessSourceName=var_str('nodeIndexPostprocessSource')
    return
  end function indexPostprocessSourceName

  function indexPostprocessSourceDescription(self)
    !!{RST
    Return a description of the source index property.
    !!}
    implicit none
    type (varying_string                             )                :: indexPostprocessSourceDescription
    class(nodePropertyExtractorIndexPostprocessSource), intent(inout) :: self
    !$GLC attributes unused :: self

    indexPostprocessSourceDescription=var_str('Index (counting from zero) of this galaxy in the `nodeData` datasets of the file from which it was postprocessed.')
    return
  end function indexPostprocessSourceDescription
