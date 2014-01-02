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

!% Contains a module which handles outputting of tree final descendent data to the \glc\ output file.

module Galacticus_Output_Trees_Final_Descendents
  !% Handles outputting of tree final descendent data to the \glc\ output file.
  use Galacticus_Nodes
  implicit none
  private
  public :: Galacticus_Output_Tree_Final_Descendents, Galacticus_Output_Tree_Final_Descendents_Property_Count,&
       & Galacticus_Output_Tree_Final_Descendents_Names

  ! Number of descendent properties.
  integer                     , parameter :: descendentPropertyCount=1

  ! Flag indicating whether or not descendent index information is to be output.
  logical                                 :: outputFinalDescendentIndices

  ! Flag indicating whether or not this module has been initialized.
  logical                                 :: outputFinalDescendentsInitialized=.false.

contains

  subroutine Galacticus_Output_Tree_Final_Descendents_Initialize
    !% Initializes the module by determining whether or not descendent data should be output.
    use Input_Parameters
    use Galacticus_Error
    use ISO_Varying_String
    implicit none

    if (.not.outputFInalDescendentsInitialized) then
       !$omp critical(Galacticus_Output_Tree_Final_Descendents_Initialize)
       if (.not.outputFinalDescendentsInitialized) then
          !@ <inputParameter>
          !@   <name>outputFinalDescendentIndices</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not final descendent indices (i.e. index of the node at the base of the tree) should be included in the output.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('outputFinalDescendentIndices',outputFinalDescendentIndices,defaultValue=.false.)
          ! Ensure that the satellite timeOfMerging property is gettable.
          if     (                                                                                                                &
               &   outputFinalDescendentIndices                                                                                   &
               &  .and.                                                                                                           &
               &   .not.defaultSatelliteComponent%timeOfMergingIsGettable()                                                       &
               & ) call Galacticus_Error_Report                                                                                   &
               &        (                                                                                                         &
               &         'Galacticus_Output_Tree_Final_Descendents_Initialize'                                                  , &
               &         'the satellite timeOfMerging property must be gettable to output descendent indices.'//                  &
               &         Galacticus_Component_List(                                                                               &
               &                                   'satellite'                                                                  , &
               &                                   defaultSatelliteComponent%timeOfMergingAttributeMatch(requireGettable=.true.)  &
               &                                  )                                                                               &
               &        )
          ! Flag that module is now initialized.
          outputFinalDescendentsInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Tree_Final_Descendents_Initialize)
    end if
    return
  end subroutine Galacticus_Output_Tree_Final_Descendents_Initialize

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Final_Descendents_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Final_Descendents</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Final_Descendents_Names(thisNode,integerProperty,integerPropertyNames,integerPropertyComments&
       &,integerPropertyUnitsSI,doubleProperty ,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set the names of descendent properties to be written to the \glc\ output file.
    implicit none
    type     (treeNode), intent(inout), pointer      :: thisNode
    double precision   , intent(in   )               :: time
    integer            , intent(inout)               :: integerProperty,doubleProperty
    character(len=*   ), intent(inout), dimension(:) :: integerPropertyNames,integerPropertyComments,doublePropertyNames &
         &,doublePropertyComments
    double precision   , intent(inout), dimension(:) :: integerPropertyUnitsSI,doublePropertyUnitsSI

    ! Initialize the module.
    call Galacticus_Output_Tree_Final_Descendents_Initialize()

    ! Return property names if we are outputting final descendant data.
    if (outputFinalDescendentIndices) then
       integerProperty=integerProperty+1
       !@ <outputProperty>
       !@   <name>finalDescendentIndex</name>
       !@   <datatype>integer</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>ID of the node which this node will have descended into at the base of the tree.</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       integerPropertyNames   (integerProperty)='finalFescendentIndex'
       integerPropertyComments(integerProperty)='ID of the node which this node will have descended into at the base of the tree.'
       integerPropertyUnitsSI (integerProperty)=0.0d0
    end if
    return
  end subroutine Galacticus_Output_Tree_Final_Descendents_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Final_Descendents_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Final_Descendents</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Final_Descendents_Property_Count(thisNode,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of descendent properties to be written to the \glc\ output file.
    use Galacticus_Nodes
    implicit none
    type   (treeNode              ), intent(inout), pointer :: thisNode
    double precision               , intent(in   )          :: time
    integer                        , intent(inout)          :: integerPropertyCount,doublePropertyCount

    ! Initialize the module.
    call Galacticus_Output_Tree_Final_Descendents_Initialize()

    ! Return property names if we are outputting virial data.
    if (outputFinalDescendentIndices) integerPropertyCount=integerPropertyCount+descendentPropertyCount
    return
  end subroutine Galacticus_Output_Tree_Final_Descendents_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Final_Descendents</unitName>
  !#  <sortName>Galacticus_Output_Tree_Final_Descendents</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Final_Descendents(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store descendent properties in the \glc\ output file buffers.
    use Galacticus_Nodes
    use Kind_Numbers
    use Galacticus_Output_Times
    implicit none
    double precision                        , intent(in   )          :: time
    type            (treeNode              ), intent(inout), pointer :: thisNode
    integer                                 , intent(inout)          :: integerProperty,integerBufferCount,doubleProperty&
         &,doubleBufferCount
    integer         (kind=kind_int8        ), intent(inout)          :: integerBuffer(:,:)
    double precision                        , intent(inout)          :: doubleBuffer (:,:)
    type            (treeNode              ),                pointer :: descendentNode
    class           (nodeComponentSatellite),                pointer :: thisSatelliteComponent

    ! Initialize the module.
    call Galacticus_Output_Tree_Final_Descendents_Initialize()

    ! If outputting indices, find them,
    if (outputFinalDescendentIndices) then
       integerProperty=integerProperty+1
       descendentNode => thisNode
       do while (associated(descendentNode%parent))
          ! Get the satellite component.
          thisSatelliteComponent => descendentNode%satellite()
          if (descendentNode%isSatellite()) then
             ! Move to the node with which the satellite merges.
             descendentNode => descendentNode%mergesWith()
          else
             descendentNode => descendentNode%parent
          end if
       end do
       ! Store the final descendent index.
       integerBuffer(integerBufferCount,integerProperty)=descendentNode%index()
    end if
    return
  end subroutine Galacticus_Output_Tree_Final_Descendents

end module Galacticus_Output_Trees_Final_Descendents
