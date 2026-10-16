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
Implements a node property extractor for the index of the host node.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorIndicesHost" docformat="rst">
   <description>
   A node property extractor which extracts the index of the node which hosts a given node. For a node which is not a subhalo the node is considered to be its own host, and so its own index is extracted---this matches the convention used for the ``hostIndex`` dataset in merger tree files, and means that this extractor never returns :math:`-1`.

   By default the index of the *direct* host is extracted, so that a sub-subhalo reports the subhalo in which it resides. If ``topLevel`` is set to ``true`` the hierarchy is instead followed upward until an isolated halo is reached, so that every node reports the isolated halo which ultimately contains it.

   Note that this differs from the ``parentIndex`` property of the :galacticus-class:`nodePropertyExtractorNodeIndices` extractor, which is the host only for subhalos (for an isolated halo it is instead the halo into which that halo grows at the next timestep of the merger tree, or :math:`-1` at the base of the tree). These conventions are described in `Merger Tree Structure in the Output &lt;https://galacticus.readthedocs.io/en/latest/manuals/user-guide/output-tree-structure.html&gt;`_.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorIntegerScalar) :: nodePropertyExtractorIndicesHost
     !!{RST
     A host node index property extractor class.
     !!}
     private
     logical :: topLevel
   contains
     procedure :: extract     => indicesHostExtract
     procedure :: name        => indicesHostName
     procedure :: description => indicesHostDescription
  end type nodePropertyExtractorIndicesHost

  interface nodePropertyExtractorIndicesHost
     !!{RST
     Constructors for the :galacticus-class:`nodePropertyExtractorIndicesHost` property extractor class.
     !!}
     module procedure indicesHostConstructorParameters
     module procedure indicesHostConstructorInternal
  end interface nodePropertyExtractorIndicesHost

contains

  function indicesHostConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodePropertyExtractorIndicesHost` property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nodePropertyExtractorIndicesHost)                :: self
    type   (inputParameters                 ), intent(inout) :: parameters
    logical                                                  :: topLevel

    !![
    <inputParameter docformat="rst">
      <name>topLevel</name>
      <defaultValue>.false.</defaultValue>
      <description>
      If true, output the index of the host at the top level of the hierarchy, otherwise output the index of the direct host.
      </description>
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
    !!{RST
    Internal constructor for the :galacticus-class:`nodePropertyExtractorIndicesHost` property extractor class.
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
    !!{RST
    Implement a ``indicesHost`` node property extractor.
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
    !!{RST
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
    !!{RST
    Return a description of the indicesHost property.
    !!}
    implicit none
    type (varying_string                  )                :: indicesHostDescription
    class(nodePropertyExtractorIndicesHost), intent(inout) :: self
    !$GLC attributes unused :: self

    indicesHostDescription=var_str('Index of the node which hosts this node - equal to the index of this node itself if it is not a subhalo.')
    return
  end function indicesHostDescription
