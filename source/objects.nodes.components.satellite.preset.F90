!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Satellite_Preset_Promote, Node_Component_Satellite_Preset_Inter_Tree_Attach, &
       &    Node_Component_Satellite_Preset_Inter_Tree_Insert

  !# <component>
  !#  <class>satellite</class>
  !#  <name>preset</name>
  !#  <isDefault>false</isDefault>
  !#  <properties>
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

  !# <nodePromotionTask>
  !#  <unitName>Node_Component_Satellite_Preset_Promote</unitName>
  !# </nodePromotionTask>
  subroutine Node_Component_Satellite_Preset_Promote(thisNode)
    !% Ensure that {\normalfont \ttfamily thisNode} is ready for promotion to its parent. In this case, we simply copy any preset satellite orbit
    !% from the parent.
    implicit none
    type (treeNode), intent(inout), pointer :: thisNode

    ! Return immediately if the preset satellite implementation is not active.
    if (.not.defaultSatelliteComponent%presetIsActive()) return
    ! Move the satellite orbit from the parent node.
    call thisnode%parent%satelliteMove(thisNode,overwrite=.true.)
    return
  end subroutine Node_Component_Satellite_Preset_Promote

  !# <interTreeSatelliteInsert>
  !#  <unitName>Node_Component_Satellite_Preset_Inter_Tree_Insert</unitName>
  !# </interTreeSatelliteInsert>
  subroutine Node_Component_Satellite_Preset_Inter_Tree_Insert(node,replaceNode)
    !% A satellite node is being moved between trees, and being added as a new satellite. Its (future-)histories will have been
    !% assigned to the {\normalfont \ttfamily replaceNode} so must be transferred.
    use Histories
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
    use Histories
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
       ! The node has a history - combined it with our own.
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
    call attachSatellite%timeOfMergingSet(pullSatellite%timeOfMerging())
    return
  end subroutine Node_Component_Satellite_Preset_Inter_Tree_Attach
  
end module Node_Component_Satellite_Preset
