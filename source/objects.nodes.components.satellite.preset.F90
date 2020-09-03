!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

!+    Contributions to this file made by:  Andrew Benson, Jianling Gan.

!% Contains a module which implements a preset satellite orbit component.

module Node_Component_Satellite_Preset
  !% Implements a preset satellite orbit component.
  implicit none
  private
  public :: Node_Component_Satellite_Preset_Satellite_Host_Change , Node_Component_Satellite_Preset_Inter_Tree_Attach, &
       &    Node_Component_Satellite_Preset_Inter_Tree_Insert     , Node_Component_Satellite_Preset_Rate_Compute     , &
       &    Node_Component_Satellite_Preset_Inter_Tree_Postprocess, Node_Component_Satellite_Preset_Thread_Initialize, &
       &    Node_Component_Satellite_Preset_Thread_Uninitialize

  !# <component>
  !#  <class>satellite</class>
  !#  <name>preset</name>
  !#  <isDefault>false</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>isOrphan</name>
  !#     <type>logical</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </property>
  !#   <property>
  !#     <name>mergeTime</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" isVirtual="true" />
  !#     <classDefault>-1.0d0</classDefault>
  !#     <getFunction>Node_Component_Satellite_Preset_Merge_Time</getFunction>
  !#     <setFunction>Node_Component_Satellite_Preset_Merge_Time_Set</setFunction>
  !#     <output unitsInSI="gigaYear" comment="Time until satellite merges."/>
  !#   </property>
  !#   <property>
  !#     <name>timeOfMerging</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <classDefault>-1.0d0</classDefault>
  !#   </property>
  !#   <property>
  !#     <name>boundMass</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="false" isGettable="true" isEvolvable="false" isVirtual="true" />
  !#     <getFunction>SatellitePresetMergeBoundMass</getFunction>
  !#     <classDefault>selfBasicComponent%mass()</classDefault>
  !#     <output unitsInSI="massSolar" comment="Bound mass of the node."/>
  !#   </property>
  !#   <property>
  !#     <name>boundMassHistory</name>
  !#     <type>history</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </property>
  !#   <property>
  !#     <name>virialOrbit</name>
  !#     <type>keplerOrbit</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </property>
  !#   <property>
  !#     <name>nodeIndexHistory</name>
  !#     <type>longIntegerHistory</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </property>
  !#   <property>
  !#     <name>nodeIndex</name>
  !#     <type>longInteger</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="false" isGettable="true" isEvolvable="false" />
  !#     <getFunction>SatellitePresetNodeIndex</getFunction>
  !#     <output unitsInSI="0.0d0" comment="Index of the satellite node used in preset satellite evolution."/>
  !#     <classDefault>self%hostNode%index()</classDefault>
  !#   </property>
  !#  </properties>
  !#  <functions>objects.nodes.components.satellite.preset.bound_functions.inc</functions>
  !# </component>

contains

  !# <nodeComponentThreadInitializationTask>
  !#  <unitName>Node_Component_Satellite_Preset_Thread_Initialize</unitName>
  !# </nodeComponentThreadInitializationTask>
  subroutine Node_Component_Satellite_Preset_Thread_Initialize(parameters_)
    !% Initializes the tree node scale dark matter profile module.
    use :: Events_Hooks    , only : nodePromotionEvent       , openMPThreadBindingAtLevel
    use :: Galacticus_Nodes, only : defaultSatelliteComponent
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_
    !$GLC attributes unused :: parameters_

    if (defaultSatelliteComponent%presetIsActive()) &
         call nodePromotionEvent%attach(defaultSatelliteComponent,nodePromotion,openMPThreadBindingAtLevel,label='nodeComponentSatellitePreset')
    return
  end subroutine Node_Component_Satellite_Preset_Thread_Initialize

  !# <nodeComponentThreadUninitializationTask>
  !#  <unitName>Node_Component_Satellite_Preset_Thread_Uninitialize</unitName>
  !# </nodeComponentThreadUninitializationTask>
  subroutine Node_Component_Satellite_Preset_Thread_Uninitialize()
    !% Uninitializes the tree node scale dark matter profile module.
    use :: Events_Hooks    , only : nodePromotionEvent
    use :: Galacticus_Nodes, only : defaultSatelliteComponent
    implicit none

    if (defaultSatelliteComponent%presetIsActive()) &
         & call nodePromotionEvent%detach(defaultSatelliteComponent,nodePromotion)
    return
  end subroutine Node_Component_Satellite_Preset_Thread_Uninitialize

  subroutine nodePromotion(self,node)
    !% Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply copy any preset satellite orbit
    !% from the parent.
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    class(*       ), intent(inout)         :: self
    type (treeNode), intent(inout), target :: node
    !$GLC attributes unused :: self
    
    ! Move the satellite orbit from the parent node.
    call node%parent%satelliteMove(node,overwrite=.true.)
    return
  end subroutine nodePromotion

  !# <interTreeSatelliteInsert>
  !#  <unitName>Node_Component_Satellite_Preset_Inter_Tree_Insert</unitName>
  !# </interTreeSatelliteInsert>
  subroutine Node_Component_Satellite_Preset_Inter_Tree_Insert(node,replaceNode)
    !% A satellite node is being moved between trees, and being added as a new satellite. Its (future-)histories will have been
    !% assigned to the {\normalfont \ttfamily replaceNode} so must be transferred.
    use :: Galacticus_Nodes, only : defaultSatelliteComponent, nodeComponentBasic, nodeComponentSatellite, treeNode
    use :: Histories       , only : history                  , longIntegerHistory
    implicit none
    type (treeNode              ), intent(inout), pointer :: node            , replaceNode
    class(nodeComponentSatellite)               , pointer :: satellite       , replaceSatellite
    class(nodeComponentBasic    )               , pointer :: basic
    type (longIntegerHistory    )                         :: historyIndex    , replaceHistoryIndex, &
         &                                                   moveHistoryIndex
    type (history               )                         :: historyMass     , replaceHistoryMass , &
         &                                                   moveHistoryMass

    ! Return immediately if the preset satellite implementation is not active.
    if (.not.defaultSatelliteComponent%presetIsActive()) return
    ! Get the basic component of the pulled node.
    basic               =>            node%basic            ()
    ! Get the node index histories attached to both the pulled node and its parent.
    satellite           =>            node%satellite        ()
    replaceSatellite    =>     replaceNode%satellite        ()
    ! Transfer node index history.
    historyIndex        =         satellite%nodeIndexHistory()
    replaceHistoryIndex =  replaceSatellite%nodeIndexHistory()
    ! Cut off history in node being replaced subsequent to current time.
    call replaceHistoryIndex%trimForward(basic%time(),moveHistoryIndex)
    ! Append removed history to pulled node.
    call historyIndex       %append     (             moveHistoryIndex)
    ! Set the histories.
    call        satellite%nodeIndexHistorySet(       historyIndex)
    call replaceSatellite%nodeIndexHistorySet(replaceHistoryIndex)
    ! Transfer subhalo mass history.
    historyMass         =         satellite%boundMassHistory()
    replaceHistoryMass  =  replaceSatellite%boundMassHistory()
    ! Cut off history in node being replaced subsequent to current time.
    call replaceHistoryMass %trimForward(basic%time(),moveHistoryMass )
    ! Append removed history to pulled node.
    call historyMass        %append     (             moveHistoryMass )
    ! Set the histories.
    call        satellite%boundMassHistorySet(       historyMass )
    call replaceSatellite%boundMassHistorySet(replaceHistoryMass )
    return
  end subroutine Node_Component_Satellite_Preset_Inter_Tree_Insert

  !# <interTreeSatelliteAttach>
  !#  <unitName>Node_Component_Satellite_Preset_Inter_Tree_Attach</unitName>
  !# </interTreeSatelliteAttach>
  subroutine Node_Component_Satellite_Preset_Inter_Tree_Attach(node)
    !% A satellite node is being moved between trees and attached as the primary progenitor of an existing satellite node. Ensure
    !% that preset satellite properties are correctly handled.
    use :: Galacticus_Nodes, only : defaultSatelliteComponent, nodeComponentBasic, nodeComponentSatellite, treeNode
    use :: Histories       , only : history                  , longIntegerHistory
    implicit none
    type (treeNode              ), intent(inout), pointer :: node
    class(nodeComponentSatellite)               , pointer :: pullSatellite   , attachSatellite
    class(nodeComponentBasic    )               , pointer :: attachBasic
    type (longIntegerHistory    )                         :: pullHistoryIndex, attachHistoryIndex
    type (history               )                         :: pullHistory     , attachHistory

    ! Return immediately if the preset satellite implementation is not active.
    if (.not.defaultSatelliteComponent%presetIsActive()) return
    ! Get the node index histories attached to both the pulled node and its parent.
    pullSatellite   => node                  %satellite       ()
    attachSatellite => node           %parent%satellite       ()
    ! Combined node index histories.
    pullHistoryIndex     =  pullSatellite         %nodeIndexHistory()
    attachHistoryIndex   =  attachSatellite       %nodeIndexHistory()
    ! Combine the histories.
    if (attachHistoryIndex%exists()) then
       ! The node has a history - combine it with our own.
       call pullHistoryIndex%append(attachHistoryIndex)
    else
       ! The node attached to has no node index history. But we must still add itself as an entry at the end of the pulled node's
       ! history.
       attachBasic => node%parent%basic()
       call pullHistoryIndex%append(attachBasic%time(),[attachSatellite%nodeIndex()])
    end if
    ! Set this history in both the pulled node and the attachment node. This ensures that the history will persist after the
    ! pulled node is promoted to the attachment node.
    call pullSatellite  %nodeIndexHistorySet(pullHistoryIndex)
    call attachSatellite%nodeIndexHistorySet(pullHistoryIndex)
    ! Combined bound mass histories.
    pullHistory     =  pullSatellite         %boundMassHistory()
    attachHistory   =  attachSatellite       %boundMassHistory()
    ! Combine the histories.
    if (attachHistory%exists()) then
       ! The node has a history - combined it with our own.
       call pullHistory%append(attachHistory)
    else
       ! The node attached to has no history. But we must still add itself as an entry at the end of the pulled node's history.
       attachBasic => node%parent%basic()
       call pullHistory%append(attachBasic%time(),[attachSatellite%boundMass()])
    end if
    ! Set this history in both the pulled node and the attachment node. This ensures that the history will persist after the
    ! pulled node is promoted to the attachment node.
    call pullSatellite  %boundMassHistorySet(pullHistory)
    call attachSatellite%boundMassHistorySet(pullHistory)
    ! The merge time for the pulled node is reset to that of the attachment node.
    call pullSatellite%timeOfMergingSet(attachSatellite%timeOfMerging())
    call pullSatellite%  virialOrbitSet(attachSatellite%virialOrbit  ())
    return
  end subroutine Node_Component_Satellite_Preset_Inter_Tree_Attach

  !# <interTreePostProcess>
  !#  <unitName>Node_Component_Satellite_Preset_Inter_Tree_Postprocess</unitName>
  !# </interTreePostProcess>
  !# <subhaloPromotionPostProcess>
  !#  <unitName>Node_Component_Satellite_Preset_Inter_Tree_Postprocess</unitName>
  !# </subhaloPromotionPostProcess>
  !# <branchJumpPostProcess>
  !#  <unitName>Node_Component_Satellite_Preset_Inter_Tree_Postprocess</unitName>
  !# </branchJumpPostProcess>
  subroutine Node_Component_Satellite_Preset_Inter_Tree_Postprocess(node)
    !% For inter-tree node transfers, ensure that any orphaned mergees of the transferred node are transferred over to the new
    !% branch.
    use :: Galacticus_Display, only : Galacticus_Display_Message, Galacticus_Verbosity_Level, verbosityInfo
    use :: Galacticus_Nodes  , only : nodeComponentSatellite    , treeNode                  , treeNodeLinkedList
    use :: ISO_Varying_String, only : var_str                   , varying_string            , operator(//)
    use :: String_Handling   , only : operator(//)
    implicit none
    type (treeNode              ), intent(inout), pointer :: node
    type (treeNode              )               , pointer :: mergee         , nodeWork
    type (treeNodeLinkedList    )               , pointer :: nodeStack      , nodeNext, &
         &                                                   nodeNew
    class(nodeComponentSatellite)               , pointer :: satelliteMergee
    type (varying_string        )                         :: message

    nodeStack => null()
    mergee => node%firstMergee
    do while (associated(mergee))
       satelliteMergee => mergee%satellite()
       if (satelliteMergee%isOrphan()) then
          if (Galacticus_Verbosity_Level() >= verbosityInfo) then
             message=var_str('Satellite node [')//mergee%index()//'] will be orphanized due to event'
             call Galacticus_Display_Message(message)
          end if
          allocate(nodeNew)
          nodeNew  %node => mergee
          nodeNew  %next => nodeStack
          nodeStack      => nodeNew
       end if
       mergee => mergee%siblingMergee
    end do
    ! Process the stack.
    do while (associated(nodeStack))
       ! Pop a node from the stack.
       nodeWork => nodeStack%node
       nodeNext => nodeStack%next
       deallocate(nodeStack)
       nodeStack => nodeNext
       ! Push any mergees onto the stack.
       mergee => nodeWork%firstMergee
       do while (associated(mergee))
          ! Only push orphaned nodes onto the stack
          satelliteMergee => mergee%satellite()
          if (satelliteMergee%isOrphan()) then
             allocate(nodeNew)
             nodeNew  %node => mergee
             nodeNew  %next => nodeStack
             nodeStack      => nodeNew
          end if
          mergee => mergee%siblingMergee
       end do
       ! Process the node.
       call Node_Component_Satellite_Preset_Orphanize(nodeWork)
    end do
    return
  end subroutine Node_Component_Satellite_Preset_Inter_Tree_Postprocess

  !# <satelliteHostChangeTask>
  !#  <unitName>Node_Component_Satellite_Preset_Satellite_Host_Change</unitName>
  !# </satelliteHostChangeTask>
  subroutine Node_Component_Satellite_Preset_Satellite_Host_Change(node)
    !% For satellite host changes, if the satellite is an orphan with a merge target ensure it remains in the branch of its merge
    !% target.
    use :: Galacticus_Display, only : Galacticus_Display_Message, Galacticus_Verbosity_Level, verbosityInfo
    use :: Galacticus_Nodes  , only : nodeComponentSatellite    , treeNode                  , treeNodeLinkedList
    use :: ISO_Varying_String, only : var_str                   , varying_string            , operator(//)
    use :: String_Handling   , only : operator(//)
    implicit none
    type (treeNode              ), intent(inout), target  :: node
    class(nodeComponentSatellite)               , pointer :: satellite
    type (treeNodeLinkedList    )               , pointer :: nodeStack, nodeNext, &
         &                                                   nodeNew
    type (treeNode              )               , pointer :: nodeWork , mergee
    type (varying_string        )                         :: message

    satellite => node%satellite()
    if (satellite%isOrphan().and.associated(node%mergeTarget)) then
       if (Galacticus_Verbosity_Level() >= verbosityInfo) then
          message=var_str('Satellite node [')//node%index()//'] will be orphanized due to host change'
          call Galacticus_Display_Message(message)
       end if
       ! Initialize a stack of nodes to allow us to process all mergees.
       allocate(nodeStack)
       nodeStack%node => node
       nodeStack%next => null()
       ! Process the stack.
       do while (associated(nodeStack))
          ! Pop a node from the stack.
          nodeWork => nodeStack%node
          nodeNext => nodeStack%next
          deallocate(nodeStack)
          nodeStack => nodeNext
          ! Push any mergees onto the stack.
          mergee => nodeWork%firstMergee
          do while (associated(mergee))
             ! Only push orphaned nodes onto the stack
             satellite => mergee%satellite()
             if (satellite%isOrphan()) then
                allocate(nodeNew)
                nodeNew  %node => mergee
                nodeNew  %next => nodeStack
                nodeStack      => nodeNew
             end if
             mergee => mergee%siblingMergee
          end do
          ! Process the node.
          call Node_Component_Satellite_Preset_Orphanize(nodeWork)
       end do
    end if
    return
  end subroutine Node_Component_Satellite_Preset_Satellite_Host_Change

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Satellite_Preset_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Satellite_Preset_Rate_Compute(node,interrupt,interruptProcedure,propertyType)
    !% Interrupt differential evolution when a preset satellite becomes an orphan.
    use :: Galacticus_Nodes, only : interruptTask, nodeComponentBasic       , nodeComponentSatellite, propertyTypeInactive, &
          &                         treeNode     , defaultSatelliteComponent
    use :: Histories       , only : history
    implicit none
    type     (treeNode              ), intent(inout), pointer :: node
    logical                          , intent(inout)          :: interrupt
    procedure(interruptTask         ), intent(inout), pointer :: interruptProcedure
    integer                          , intent(in   )          :: propertyType
    class    (nodeComponentBasic    )               , pointer :: basic
    class    (nodeComponentSatellite)               , pointer :: satellite
    type     (history               )                         :: historyBoundMass
    logical                                                   :: exceedsHistoryTime
    
    ! Return immediately if the preset satellite implementation is not active.
    if (.not.defaultSatelliteComponent%presetIsActive()) return
    ! Return immediately if inactive variables are requested.
    if (propertyType == propertyTypeInactive) return
    basic            => node     %basic           ()
    satellite        => node     %satellite       ()
    historyBoundMass =  satellite%boundMassHistory()

    if (historyBoundMass%exists()) then
       exceedsHistoryTime=basic%time() >= historyBoundMass%time(size(historyBoundMass%time))
    else
       exceedsHistoryTime=.true.
    end if

    if (.not.satellite%isOrphan() .and. node%isSatellite() .and. associated(node%mergeTarget) .and. exceedsHistoryTime) then
       interrupt          =  .true.
       interruptProcedure => Node_Component_Satellite_Preset_Orphanize
    end if
    return
  end subroutine Node_Component_Satellite_Preset_Rate_Compute

  subroutine Node_Component_Satellite_Preset_Orphanize(node)
    !% Handle orphanization of a preset satellite component. The satellite should be moved to the branch of its target node.
    use :: Galacticus_Display, only : Galacticus_Display_Message, Galacticus_Verbosity_Level, verbosityInfo
    use :: Galacticus_Nodes  , only : nodeComponentBasic        , nodeComponentSatellite    , treeNode
    use :: ISO_Varying_String, only : var_str                   , varying_string            , operator(//)
    use :: String_Handling   , only : operator(//)
    implicit none
    type (treeNode              ), intent(inout), target  :: node
    type (treeNode              )               , pointer :: nodeHost
    class(nodeComponentBasic    )               , pointer :: basic    , basicHost
    class(nodeComponentSatellite)               , pointer :: satellite
    type (varying_string        )                         :: message

    satellite => node    %satellite  ()
    call satellite%isOrphanSet(.true.)
    nodeHost  => node    %mergeTarget
    ! For satellite merge targets, step up through parents until an isolated host is found.
    do while (nodeHost%isSatellite())
       nodeHost => nodeHost%parent
    end do
    basic     => node    %basic      ()
    basicHost => nodeHost%basic      ()
    ! Trace the merge target progenitors back until one is found which exists at the time of the orphaned satellite.
    do while (basicHost%time() > basic%time())
       ! For satellite merge targets, step up through parents until an isolated host is found.
       do while (nodeHost%isSatellite())
          nodeHost => nodeHost%parent
       end do
       ! If a progenitor exists, move to it.
       if (associated(nodeHost%firstChild)) then
          nodeHost  => nodeHost%firstChild
          basicHost => nodeHost%basic     ()
       else
          ! No further progenitors exist, so stop here.
          exit
       end if
    end do
    ! Report.
    if (Galacticus_Verbosity_Level() >= verbosityInfo) then
       message=var_str('Satellite node [')//node%index()//'] is being orphanized'
       if (associated(node%parent,nodeHost)) then
          message=message//' - remains in same host ['//nodeHost%index()//']'
       else
          message=message//' - moves from host ['//node%parent%index()//'] to host ['//nodeHost%index()//']'
       end if
       call Galacticus_Display_Message(message)
    end if
    ! Move to the new host. (If the new host is the same as the current host, do nothing.)
    if (.not.associated(node%parent,nodeHost)) then
       if (associated(node%parent)) call node%removeFromHost()
       node    %sibling        => nodeHost%firstSatellite
       node    %parent         => nodeHost
       nodeHost%firstSatellite => node
    end if
    return
  end subroutine Node_Component_Satellite_Preset_Orphanize

end module Node_Component_Satellite_Preset
