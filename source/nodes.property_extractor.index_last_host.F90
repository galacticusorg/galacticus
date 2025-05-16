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
Implements a node property extractor for the index of the last host node.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorIndexLastHost">
   <description>A last host node index property extractor.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorIntegerScalar) :: nodePropertyExtractorIndexLastHost
     !!{
     A last host node index property extractor.
     !!}
     private
     integer :: indexLastHostID
   contains
     procedure :: extract     => indexLastHostExtract
     procedure :: name        => indexLastHostName
     procedure :: description => indexLastHostDescription
  end type nodePropertyExtractorIndexLastHost

  interface nodePropertyExtractorIndexLastHost
     !!{
     Constructors for the \refClass{nodePropertyExtractorIndexLastHost} output analysis class.
     !!}
     module procedure indexLastHostConstructorParameters
     module procedure indexLastHostConstructorInternal
  end interface nodePropertyExtractorIndexLastHost

contains

  function indexLastHostConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorIndexLastHost} node property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorIndexLastHost)                :: self
    type(inputParameters                   ), intent(inout) :: parameters

    self=nodePropertyExtractorIndexLastHost()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function indexLastHostConstructorParameters

  function indexLastHostConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorIndexLastHost} node property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorIndexLastHost) :: self

    !![
    <addMetaProperty component="basic" name="nodeIndexLastHost" type="longInteger" id="self%indexLastHostID" isCreator="no"/>
    !!]
    return
  end function indexLastHostConstructorInternal

  function indexLastHostExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily indexLastHost} node property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    integer         (kind_int8                         )                          :: indexLastHostExtract
    class           (nodePropertyExtractorIndexLastHost), intent(inout)           :: self
    type            (treeNode                          ), intent(inout), target   :: node
    double precision                                    , intent(in   )           :: time
    type            (multiCounter                      ), intent(inout), optional :: instance
    class           (nodeComponentBasic                )               , pointer  :: basic
    !$GLC attributes unused :: instance, time

    basic                => node %basic                          (                     )
    indexLastHostExtract =  basic%longIntegerRank0MetaPropertyGet(self%indexLastHostID)
    return
  end function indexLastHostExtract

  function indexLastHostName(self)
    !!{
    Return the name of the branch tip index property.
    !!}
    implicit none
    type (varying_string                    )                :: indexLastHostName
    class(nodePropertyExtractorIndexLastHost), intent(inout) :: self
    !$GLC attributes unused :: self

    indexLastHostName=var_str('nodeIndexLastHost')
    return
  end function indexLastHostName

  function indexLastHostDescription(self)
    !!{
    Return a description of the branch tip index property.
    !!}
    implicit none
    type (varying_string                    )                :: indexLastHostDescription
    class(nodePropertyExtractorIndexLastHost), intent(inout) :: self
    !$GLC attributes unused :: self

    indexLastHostDescription=var_str('Index of the node in whcih this node was last a satellite.')
    return
  end function indexLastHostDescription

