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
Implements a host index output analysis property extractor class.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorIndicesHost">
   <description>
    A node property extractor which extracts the index of the node which hosts a given node. For unhosted nodes (i.e. nodes
    which are not subhalos), a value of $-1$ is extracted instead.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorIntegerScalar) :: nodePropertyExtractorIndicesHost
     !!{
     A host index output analysis class.
     !!}
     private
     logical :: topLevel
   contains
     procedure :: extract     => indicesHostExtract
     procedure :: name        => indicesHostName
     procedure :: description => indicesHostDescription
  end type nodePropertyExtractorIndicesHost

  interface nodePropertyExtractorIndicesHost
     !!{
     Constructors for the \refClass{nodePropertyExtractorIndicesHost} output analysis class.
     !!}
     module procedure indicesHostConstructorParameters
     module procedure indicesHostConstructorInternal
  end interface nodePropertyExtractorIndicesHost

contains

  function indicesHostConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorIndicesHost} node property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nodePropertyExtractorIndicesHost)                :: self
    type   (inputParameters                 ), intent(inout) :: parameters
    logical                                                  :: topLevel

    !![
    <inputParameter>
      <name>topLevel</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, output the index of the host at the top level of the hierarchy, otherwise output the index of the direct host.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
     self=nodePropertyExtractorIndicesHost(topLevel)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function indicesHostConstructorParameters

  function indicesHostConstructorInternal(topLevel) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorIndicesHost} node property extractor class.
    !!}
    implicit none
    type   (nodePropertyExtractorIndicesHost)                :: self
    logical                                  , intent(in   ) :: topLevel
    !![
    <constructorAssign variables="topLevel"/>
    !!]

    return
  end function indicesHostConstructorInternal

  function indicesHostExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily indicesHost} node property extractor.
    !!}
    implicit none
    integer         (kind_int8                       )                          :: indicesHostExtract
    class           (nodePropertyExtractorIndicesHost), intent(inout)           :: self
    type            (treeNode                        ), intent(inout), target   :: node
    double precision                                  , intent(in   )           :: time
    type            (multiCounter                    ), intent(inout), optional :: instance
    type            (treeNode                        ), pointer                 :: nodeHost
    !$GLC attributes unused :: self, instance, time

    if (node%isSatellite()) then
       nodeHost => node%parent
       if (self%topLevel) then
          do while (nodeHost%isSatellite())
             nodeHost => nodeHost%parent
          end do
       end if
    else
       nodeHost => node
    end if
    indicesHostExtract=nodeHost%index()
    return
  end function indicesHostExtract


  function indicesHostName(self)
    !!{
    Return the name of the indicesHost property.
    !!}
    implicit none
    type (varying_string                  )                :: indicesHostName
    class(nodePropertyExtractorIndicesHost), intent(inout) :: self
    !$GLC attributes unused :: self

    indicesHostName=var_str('hostIndex')
    return
  end function indicesHostName

  function indicesHostDescription(self)
    !!{
    Return a description of the indicesHost property.
    !!}
    implicit none
    type (varying_string                  )                :: indicesHostDescription
    class(nodePropertyExtractorIndicesHost), intent(inout) :: self
    !$GLC attributes unused :: self

    indicesHostDescription=var_str('ID of the node which hosts this node (or -1 is there is no host).')
    return
  end function indicesHostDescription

