!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which handles outputting of tree index data to the \glc\ output file.

module Galacticus_Output_Trees_Tree_Indices
  !% Handles outputting of tree index data to the \glc\ output file.
  use Galacticus_Nodes
  implicit none
  private
  public :: Galacticus_Output_Tree_Tree_Indices, Galacticus_Output_Tree_Tree_Indices_Property_Count, Galacticus_Output_Tree_Tree_Indices_Names

  ! Number of properties.
  integer, parameter :: indexPropertyCount   =1

  ! Flag indicating whether or not tree index is to be output.
  logical            :: outputTreeIndices

  ! Flag indicating whether or not this module has been initialized.
  logical            :: outputTreeInitialized=.false.

contains

  subroutine Galacticus_Output_Tree_Indices_Initialize()
    !% Initializes the module by determining whether or not tree index data should be output.
    use Input_Parameters
    use Galacticus_Error
    implicit none

    if (.not.outputTreeInitialized) then
       !$omp critical(Galacticus_Output_Tree_Indices_Initialize)
       if (.not.outputTreeInitialized) then
          !@ <inputParameter>
          !@   <name>outputTreeIndices</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not descendent indices (i.e. index of the node at the next output) should be included in the output.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('outputTreeIndices',outputTreeIndices,defaultValue=.false.)
          ! Flag that module is now initialized.
          outputTreeInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Tree_Indices_Initialize)
    end if
    return
  end subroutine Galacticus_Output_Tree_Indices_Initialize

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Tree_Indices_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Tree_Indices</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Tree_Indices_Names(thisNode,integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty&
       &,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set the names of tree index properties to be written to the \glc\ output file.
    implicit none
    type            (treeNode)              , intent(inout), pointer :: thisNode
    double precision                        , intent(in   )          :: time
    integer                                 , intent(inout)          :: doubleProperty         , integerProperty
    character       (len=*   ), dimension(:), intent(inout)          :: doublePropertyComments , doublePropertyNames   , &
         &                                                              integerPropertyComments, integerPropertyNames
    double precision          , dimension(:), intent(inout)          :: doublePropertyUnitsSI  , integerPropertyUnitsSI

    ! Initialize the module.
    call Galacticus_Output_Tree_Indices_Initialize()

    if (outputTreeIndices) then
       integerProperty=integerProperty+1
       !@ <outputProperty>
       !@   <name>mergerTreeIndex</name>
       !@   <datatype>integer</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Tree index for this node.</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       integerPropertyNames   (integerProperty)='mergerTreeIndex'
       integerPropertyComments(integerProperty)='Tree index for this node.'
       integerPropertyUnitsSI (integerProperty)=0.0d0
    end if
    return
  end subroutine Galacticus_Output_Tree_Tree_Indices_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Tree_Indices_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Tree_Indices</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Tree_Indices_Property_Count(thisNode,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of link properties to be written to the \glc\ output file.
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode
    double precision          , intent(in   )          :: time
    integer                   , intent(inout)          :: doublePropertyCount, integerPropertyCount

    ! Initialize the module.
    call Galacticus_Output_Tree_Indices_Initialize()

    if (outputTreeIndices) integerPropertyCount=integerPropertyCount+indexPropertyCount
    return
  end subroutine Galacticus_Output_Tree_Tree_Indices_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Tree_Indices</unitName>
  !#  <sortName>Galacticus_Output_Tree_Tree_Indices</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Tree_Indices(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
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

    ! Initialize the module.
    call Galacticus_Output_Tree_Indices_Initialize()

    if (outputTreeIndices) then
       integerProperty=integerProperty+1
       integerBuffer(integerBufferCount,integerProperty)=thisNode%hostTree%index
    end if
    return
  end subroutine Galacticus_Output_Tree_Tree_Indices

end module Galacticus_Output_Trees_Tree_Indices
