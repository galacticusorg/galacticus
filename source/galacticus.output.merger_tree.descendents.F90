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

!% Contains a module which handles outputting of tree descendent data to the \glc\ output file.

module Galacticus_Output_Trees_Descendents
  !% Handles outputting of tree descendent data to the \glc\ output file.
  use Galacticus_Nodes
  implicit none
  private
  public :: Galacticus_Output_Tree_Descendents, Galacticus_Output_Tree_Descendents_Property_Count, Galacticus_Output_Tree_Descendents_Names

  ! Number of descendent properties.
  integer, parameter :: descendentPropertyCount     =1

  ! Flag indicating whether or not descendent index information is to be output.
  logical            :: outputDescendentIndices

  ! Flag indicating whether or not this module has been initialized.
  logical            :: outputDescendentsInitialized=.false.

contains

  subroutine Galacticus_Output_Tree_Descendents_Initialize
    !% Initializes the module by determining whether or not descendent data should be output.
    use Input_Parameters
    use ISO_Varying_String
    use Galacticus_Error
    implicit none

    if (.not.outputDescendentsInitialized) then
       !$omp critical(Galacticus_Output_Tree_Descendents_Initialize)
       if (.not.outputDescendentsInitialized) then
          !@ <inputParameter>
          !@   <name>outputDescendentIndices</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not descendent indices (i.e. index of the node at the next output) should be included in the output.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('outputDescendentIndices',outputDescendentIndices,defaultValue=.false.)
          ! Ensure that the satellite timeOfMerging property is gettable.
          if     (                                                                                                                &
               &   outputDescendentIndices                                                                                        &
               &  .and.                                                                                                           &
               &   .not.defaultSatelliteComponent%timeOfMergingIsGettable()                                                       &
               & ) call Galacticus_Error_Report                                                                                   &
               &        (                                                                                                         &
               &         'Galacticus_Output_Tree_Descendents_Initialize'                                                        , &
               &         'the satellite timeOfMerging property must be gettable to output descendent indices.'//                  &
               &         Galacticus_Component_List(                                                                               &
               &                                   'satellite'                                                                  , &
               &                                   defaultSatelliteComponent%timeOfMergingAttributeMatch(requireGettable=.true.)  &
               &                                  )                                                                               &
               &        )
          ! Flag that module is now initialized.
          outputDescendentsInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Tree_Descendents_Initialize)
    end if
    return
  end subroutine Galacticus_Output_Tree_Descendents_Initialize

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Descendents_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Descendents</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Descendents_Names(thisNode,integerProperty,integerPropertyNames,integerPropertyComments&
       &,integerPropertyUnitsSI,doubleProperty ,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set the names of descendent properties to be written to the \glc\ output file.
    implicit none
    type            (treeNode)              , intent(inout), pointer :: thisNode
    double precision                        , intent(in   )          :: time
    integer                                 , intent(inout)          :: doubleProperty         , integerProperty
    character       (len=*   ), dimension(:), intent(inout)          :: doublePropertyComments , doublePropertyNames   , &
         &                                                              integerPropertyComments, integerPropertyNames
    double precision          , dimension(:), intent(inout)          :: doublePropertyUnitsSI  , integerPropertyUnitsSI

    ! Initialize the module.
    call Galacticus_Output_Tree_Descendents_Initialize()

    ! Return property names if we are outputting virial data.
    if (outputDescendentIndices) then
       integerProperty=integerProperty+1
       !@ <outputProperty>
       !@   <name>descendentIndex</name>
       !@   <datatype>integer</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>ID of the node which this node will have descended into by the next timestep.</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       integerPropertyNames   (integerProperty)='descendentIndex'
       integerPropertyComments(integerProperty)='ID of the node which this node will have descended into by the next timestep.'
       integerPropertyUnitsSI (integerProperty)=0.0d0
    end if
    return
  end subroutine Galacticus_Output_Tree_Descendents_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Descendents_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Descendents</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Descendents_Property_Count(thisNode,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of descendent properties to be written to the \glc\ output file.
    implicit none
    type            (treeNode              ), intent(inout), pointer :: thisNode
    double precision                        , intent(in   )          :: time
    integer                                 , intent(inout)          :: doublePropertyCount, integerPropertyCount

    ! Initialize the module.
    call Galacticus_Output_Tree_Descendents_Initialize()

    ! Return property names if we are outputting virial data.
    if (outputDescendentIndices) integerPropertyCount=integerPropertyCount+descendentPropertyCount
    return
  end subroutine Galacticus_Output_Tree_Descendents_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Descendents</unitName>
  !#  <sortName>Galacticus_Output_Tree_Descendents</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Descendents(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store descendent properties in the \glc\ output file buffers.
    use Kind_Numbers
    use Galacticus_Output_Times
    implicit none
    double precision                        , intent(in   )          :: time
    type            (treeNode              ), intent(inout), pointer :: thisNode
    integer                                 , intent(inout)          :: doubleBufferCount          , doubleProperty, integerBufferCount, &
         &                                                              integerProperty
    integer         (kind=kind_int8        ), intent(inout)          :: integerBuffer         (:,:)
    double precision                        , intent(inout)          :: doubleBuffer          (:,:)
    type            (treeNode              )               , pointer :: descendentNode
    class           (nodeComponentBasic    )               , pointer :: thisBasicComponent
    class           (nodeComponentSatellite)               , pointer :: thisSatelliteComponent
    double precision                                                 :: outputTimeNext
    logical                                                          :: foundDescendent

    ! Initialize the module.
    call Galacticus_Output_Tree_Descendents_Initialize()

    ! Return property names if we are outputting virial data.
    if (outputDescendentIndices) then

       ! Get the satellite component.
       thisSatelliteComponent => thisNode%satellite()

       ! Get the time of the next output.
       outputTimeNext=Galacticus_Next_Output_Time(time)

       foundDescendent=.false.
       integerProperty=integerProperty+1
       if (outputTimeNext < 0.0d0) then
          ! There is no next output time.
          integerBuffer(integerBufferCount,integerProperty)=thisNode%index()
          foundDescendent=.true.
       else if (thisNode%isSatellite()) then
          ! Node is a satellite, so it's node index will remain unchanged.
          if (thisSatelliteComponent%timeOfMerging() > outputTimeNext) then
             ! Satellite will not have merged prior to the next output time, so retains its own index.
             integerBuffer(integerBufferCount,integerProperty)=thisNode%index()
             foundDescendent=.true.
          else
             ! Satellite will merge prior to the next output time - find the node it merges with.
             descendentNode => thisNode%mergesWith()
          end if
       else
          ! Node is not a satellite, so set the initial descendent to itself.
          descendentNode => thisNode
       end if
       ! Check if we still need to find the descendent.
       if (.not.foundDescendent) then
          ! No descendent has yet been found, so trace forward in time until we find one.
          ! Continue until the tree base is reached, or the next output time is reached.
          do while (.not.foundDescendent)
             ! Get the satellite component.
             thisSatelliteComponent => descendentNode%satellite()
             thisBasicComponent     => descendentNode%basic    ()
             ! If the next output time has been surpassed, then we are finished.
             if (thisBasicComponent%time() >= outputTimeNext) then
                foundDescendent=.true.
             else
                ! Test whether this node is the primary progenitor.
                if (descendentNode%isPrimaryProgenitor()) then
                   ! It is, so simply move to the parent node.
                   descendentNode => descendentNode%parent
                else
                   ! It is not, so it becomes a satellite. Test whether it has a merge target associated with it.
                   if (associated(descendentNode%mergeTarget)) then
                      ! It does. If merging occurs before the next output time, jump to that node. Otherwise we are finished.
                      if (thisSatelliteComponent%timeOfMerging() <= outputTimeNext) then
                         descendentNode => descendentNode%mergeTarget
                      else
                         foundDescendent=.true.
                      end if
                   else
                      ! We no longer can tell if this node will exist as a separate entity at the next output time. Assume that it
                      ! will, and therefore we are finished.
                      foundDescendent=.true.
                   end if
                end if
             end if
          end do
          ! If the descendent exists after the next output time, then we've gone one step too far - back up a step.
          if (thisBasicComponent%time() > outputTimeNext .and. associated(descendentNode%firstChild)) descendentNode => descendentNode%firstChild
          ! Store the descendent index.
          integerBuffer(integerBufferCount,integerProperty)=descendentNode%index()
       end if
    end if
    return
  end subroutine Galacticus_Output_Tree_Descendents

end module Galacticus_Output_Trees_Descendents
