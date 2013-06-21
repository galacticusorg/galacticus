!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which handles outputting of tree link data to the \glc\ output file.

module Galacticus_Output_Trees_Links
  !% Handles outputting of tree link data to the \glc\ output file.
  use Galacticus_Nodes
  implicit none
  private
  public :: Galacticus_Output_Tree_Links, Galacticus_Output_Tree_Links_Property_Count, Galacticus_Output_Tree_Links_Names

  ! Number of link properties.
  integer, parameter :: linkPropertyCount=6

contains

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Links_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Links</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Links_Names(thisNode,integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty&
       &,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set the names of link properties to be written to the \glc\ output file.
    implicit none
    type            (treeNode)              , intent(inout), pointer :: thisNode
    double precision                        , intent(in   )          :: time
    integer                                 , intent(inout)          :: doubleProperty         , integerProperty
    character       (len=*   ), dimension(:), intent(inout)          :: doublePropertyComments , doublePropertyNames   , &
         &                                                              integerPropertyComments, integerPropertyNames
    double precision          , dimension(:), intent(inout)          :: doublePropertyUnitsSI  , integerPropertyUnitsSI

    integerProperty=integerProperty+1
    !@ <outputProperty>
    !@   <name>nodeIndex</name>
    !@   <datatype>integer</datatype>
    !@   <cardinality>0..1</cardinality>
    !@   <description>Tree-unique ID for this node.</description>
    !@   <label>???</label>
    !@   <outputType>nodeData</outputType>
    !@ </outputProperty>
    integerPropertyNames   (integerProperty)='nodeIndex'
    integerPropertyComments(integerProperty)='Tree-unique ID for this node.'
    integerPropertyUnitsSI (integerProperty)=0.0d0
    integerProperty=integerProperty+1
    !@ <outputProperty>
    !@   <name>parentNode</name>
    !@   <datatype>integer</datatype>
    !@   <cardinality>0..1</cardinality>
    !@   <description>ID of parent node.</description>
    !@   <label>???</label>
    !@   <outputType>nodeData</outputType>
    !@ </outputProperty>
    integerPropertyNames   (integerProperty)='parentNode'
    integerPropertyComments(integerProperty)='ID of parent node.'
    integerPropertyUnitsSI (integerProperty)=0.0d0
    integerProperty=integerProperty+1
    !@ <outputProperty>
    !@   <name>childNode</name>
    !@   <datatype>integer</datatype>
    !@   <cardinality>0..1</cardinality>
    !@   <description>ID of primary child node.</description>
    !@   <label>???</label>
    !@   <outputType>nodeData</outputType>
    !@ </outputProperty>
    integerPropertyNames   (integerProperty)='childNode'
    integerPropertyComments(integerProperty)='ID of primary child node.'
    integerPropertyUnitsSI (integerProperty)=0.0d0
    integerProperty=integerProperty+1
    !@ <outputProperty>
    !@   <name>siblingNode</name>
    !@   <datatype>integer</datatype>
    !@   <cardinality>0..1</cardinality>
    !@   <description>ID of sibling node.</description>
    !@   <label>???</label>
    !@   <outputType>nodeData</outputType>
    !@ </outputProperty>
    integerPropertyNames   (integerProperty)='siblingNode'
    integerPropertyComments(integerProperty)='ID of sibling node.'
    integerPropertyUnitsSI (integerProperty)=0.0d0
    integerProperty=integerProperty+1
    !@ <outputProperty>
    !@   <name>satelliteNode</name>
    !@   <datatype>integer</datatype>
    !@   <cardinality>0..1</cardinality>
    !@   <description>ID of first satellite node.</description>
    !@   <label>???</label>
    !@   <outputType>nodeData</outputType>
    !@ </outputProperty>
    integerPropertyNames   (integerProperty)='satelliteNode'
    integerPropertyComments(integerProperty)='ID of first satellite node.'
    integerPropertyUnitsSI (integerProperty)=0.0d0
    integerProperty=integerProperty+1
    !@ <outputProperty>
    !@   <name>nodeIsIsolated</name>
    !@   <datatype>integer</datatype>
    !@   <cardinality>0..1</cardinality>
    !@   <description>Is the node isolated (0|1)?</description>
    !@   <label>???</label>
    !@   <outputType>nodeData</outputType>
    !@ </outputProperty>
    integerPropertyNames   (integerProperty)='nodeIsIsolated'
    integerPropertyComments(integerProperty)='Is the node isolated (0|1)?'
    integerPropertyUnitsSI (integerProperty)=0.0d0
    return
  end subroutine Galacticus_Output_Tree_Links_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Links_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Links</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Links_Property_Count(thisNode,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of link properties to be written to the \glc\ output file.
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode
    double precision          , intent(in   )          :: time
    integer                   , intent(inout)          :: doublePropertyCount, integerPropertyCount

    integerPropertyCount=integerPropertyCount+linkPropertyCount
    return
  end subroutine Galacticus_Output_Tree_Links_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Links</unitName>
  !#  <sortName>Galacticus_Output_Tree_Links</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Links(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store link properties in the \glc\ output file buffers.
    use Kind_Numbers
    implicit none
    double precision                , intent(in   )          :: time
    type            (treeNode      ), intent(inout), pointer :: thisNode
    integer                         , intent(inout)          :: doubleBufferCount     , doubleProperty, integerBufferCount, &
         &                                                      integerProperty
    integer         (kind=kind_int8), intent(inout)          :: integerBuffer    (:,:)
    double precision                , intent(inout)          :: doubleBuffer     (:,:)

    integerProperty=integerProperty+1
    integerBuffer(integerBufferCount,integerProperty)=thisNode               %index()
    integerProperty=integerProperty+1
    integerBuffer(integerBufferCount,integerProperty)=thisNode%parent        %index()
    integerProperty=integerProperty+1
    integerBuffer(integerBufferCount,integerProperty)=thisNode%firstChild    %index()
    integerProperty=integerProperty+1
    integerBuffer(integerBufferCount,integerProperty)=thisNode%sibling       %index()
    integerProperty=integerProperty+1
    integerBuffer(integerBufferCount,integerProperty)=thisNode%firstSatellite%index()
    integerProperty=integerProperty+1
    select case (thisNode%isSatellite())
    case (.true.)
       integerBuffer(integerBufferCount,integerProperty)=0
    case (.false.)
       integerBuffer(integerBufferCount,integerProperty)=1
    end select
    return
  end subroutine Galacticus_Output_Tree_Links

end module Galacticus_Output_Trees_Links
