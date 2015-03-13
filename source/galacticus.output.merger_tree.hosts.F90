!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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
!!    along with Galacticus.  I not, see <http://www.gnu.org/licenses/>.

!% Contains a module which handles outputting of tree host data to the \glc\ output file.

module Galacticus_Output_Trees_Hosts
  !% Handles outputting of tree host data to the \glc\ output file.
  use Galacticus_Nodes
  implicit none
  private
  public :: Galacticus_Output_Tree_Hosts, Galacticus_Output_Tree_Hosts_Property_Count, Galacticus_Output_Tree_Hosts_Names

  ! Number of host properties.
  integer, parameter :: hostPropertyCount     =1

  ! Flag indicating whether or not host index information is to be output.
  logical            :: outputHostIndices

  ! Flag indicating whether or not this module has been initialized.
  logical            :: outputHostsInitialized=.false.

contains

  subroutine Galacticus_Output_Tree_Hosts_Initialize()
    !% Initializes the module by determining whether or not host data should be output.
    use Input_Parameters
    use ISO_Varying_String
    use Galacticus_Error
    implicit none

    if (.not.outputHostsInitialized) then
       !$omp critical(Galacticus_Output_Tree_Hosts_Initialize)
       if (.not.outputHostsInitialized) then
          !@ <inputParameter>
          !@   <name>outputHostIndices</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not host indices should be included in the output.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('outputHostIndices',outputHostIndices,defaultValue=.false.)
          ! Flag that module is now initialized.
          outputHostsInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Tree_Hosts_Initialize)
    end if
    return
  end subroutine Galacticus_Output_Tree_Hosts_Initialize

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Hosts_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Hosts</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Hosts_Names(thisNode,integerProperty,integerPropertyNames,integerPropertyComments&
       &,integerPropertyUnitsSI,doubleProperty ,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set the names of host properties to be written to the \glc\ output file.
    implicit none
    type            (treeNode)              , intent(inout), pointer :: thisNode
    double precision                        , intent(in   )          :: time
    integer                                 , intent(inout)          :: doubleProperty         , integerProperty
    character       (len=*   ), dimension(:), intent(inout)          :: doublePropertyComments , doublePropertyNames   , &
         &                                                              integerPropertyComments, integerPropertyNames
    double precision          , dimension(:), intent(inout)          :: doublePropertyUnitsSI  , integerPropertyUnitsSI

    ! Initialize the module.
    call Galacticus_Output_Tree_Hosts_Initialize()

    ! Return property names if we are outputting virial data.
    if (outputHostIndices) then
       integerProperty=integerProperty+1
       !@ <outputProperty>
       !@   <name>hostIndex</name>
       !@   <datatype>integer</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>ID of the node which hosts this node (or $-1$ if there is no host).</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       integerPropertyNames   (integerProperty)='hostIndex'
       integerPropertyComments(integerProperty)='ID of the node which hosts this node (or -1 is there is no host).'
       integerPropertyUnitsSI (integerProperty)=0.0d0
    end if
    return
  end subroutine Galacticus_Output_Tree_Hosts_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Hosts_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Hosts</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Hosts_Property_Count(thisNode,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of host properties to be written to the \glc\ output file.
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode
    double precision          , intent(in   )          :: time
    integer                   , intent(inout)          :: doublePropertyCount, integerPropertyCount

    ! Initialize the module.
    call Galacticus_Output_Tree_Hosts_Initialize()
    ! Return property names if we are outputting virial data.
    if (outputHostIndices) integerPropertyCount=integerPropertyCount+hostPropertyCount
    return
  end subroutine Galacticus_Output_Tree_Hosts_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Hosts</unitName>
  !#  <sortName>Galacticus_Output_Tree_Hosts</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Hosts(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store host properties in the \glc\ output file buffers.
    use Kind_Numbers
    use Galacticus_Output_Times
    implicit none
    double precision                , intent(in   )          :: time
    type            (treeNode      ), intent(inout), pointer :: thisNode
    integer                         , intent(inout)          :: doubleBufferCount     , doubleProperty, &
         &                                                      integerBufferCount    , integerProperty
    integer         (kind=kind_int8), intent(inout)          :: integerBuffer    (:,:)
    double precision                , intent(inout)          :: doubleBuffer     (:,:)
    type            (treeNode      )               , pointer :: hostNode

    ! Initialize the module.
    call Galacticus_Output_Tree_Hosts_Initialize()
    ! Return property names if we are outputting virial data.
    if (outputHostIndices) then
       if (thisNode%isSatellite()) then
          ! Node is a satellite.
          hostNode => thisNode%parent
       else
          ! Node is not a satellite - treat as self-hosting.
          hostNode => thisNode
       end if
       ! Store the host index.
       integerProperty                                    =integerProperty  +1
       integerBuffer  (integerBufferCount,integerProperty)=hostNode%index ()
    end if
    return
  end subroutine Galacticus_Output_Tree_Hosts

end module Galacticus_Output_Trees_Hosts
