!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

!% Contains a module which implements a preset position component.

module Node_Component_Position_Preset
  !% Implements a preset position component.
  implicit none
  private
  public :: Node_Component_Position_Preset_Initialize         , Node_Component_Position_Preset_Inter_Tree_Insert, &
       &    Node_Component_Position_Preset_Move               , Node_Component_Position_Preset_Thread_Initialize, &
       &    Node_Component_Position_Preset_Thread_Uninitialize

  !# <component>
  !#  <class>position</class>
  !#  <name>preset</name>
  !#  <isDefault>false</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>position</name>
  !#     <type>double</type>
  !#     <rank>1</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <getFunction bindsTo="component">PositionPresetPosition</getFunction>
  !#     <output labels="[X,Y,Z]" unitsInSI="megaParsec" comment="Position of the node (in physical coordinates)."/>
  !#     <classDefault>[0.0d0,0.0d0,0.0d0]</classDefault>
  !#   </property>
  !#   <property>
  !#     <name>velocity</name>
  !#     <type>double</type>
  !#     <rank>1</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <getFunction bindsTo="component">PositionPresetVelocity</getFunction>
  !#     <output labels="[X,Y,Z]" unitsInSI="kilo" comment="Velocity of the node (in physical coordinates)."/>
  !#     <classDefault>[0.0d0,0.0d0,0.0d0]</classDefault>
  !#   </property>
  !#   <property>
  !#     <name>positionHistory</name>
  !#     <type>history</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </property>
  !#  </properties>
  !#  <functions>objects.nodes.components.position.preset.bound_functions.inc</functions>
  !# </component>

  ! Options.
  logical :: positionsPresetSatelliteToHost

contains

  !# <nodeComponentInitializationTask>
  !#  <unitName>Node_Component_Position_Preset_Initialize</unitName>
  !# </nodeComponentInitializationTask>
  subroutine Node_Component_Position_Preset_Initialize(parameters_)
    use :: Galacticus_Nodes, only : defaultPositionComponent
    use :: Input_Parameters, only : inputParameter          , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_

    ! Initialize the module if necessary.
    if (defaultPositionComponent%presetIsActive()) then
      ! Read parameters controlling the physical implementation.
       !# <inputParameter>
       !#   <name>positionsPresetSatelliteToHost</name>
       !#   <defaultValue>.false.</defaultValue>
       !#   <description>If true, the position of satellite halos will be adjusted to match that of their host halo.</description>
       !#   <source>parameters_</source>
       !# </inputParameter>
    end if
    return
  end subroutine Node_Component_Position_Preset_Initialize

  !# <nodeComponentThreadInitializationTask>
  !#  <unitName>Node_Component_Position_Preset_Thread_Initialize</unitName>
  !# </nodeComponentThreadInitializationTask>
  subroutine Node_Component_Position_Preset_Thread_Initialize(parameters_)
    !% Initializes the tree node scale dark matter profile module.
    use :: Events_Hooks    , only : nodePromotionEvent   , openMPThreadBindingAtLevel
    use :: Galacticus_Nodes, only : defaultPositionComponent
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_
    !$GLC attributes unused :: parameters_

    if (defaultPositionComponent%presetIsActive()) &
         call nodePromotionEvent%attach(defaultPositionComponent,nodePromotion,openMPThreadBindingAtLevel,label='nodeComponentPositionPreset')
    return
  end subroutine Node_Component_Position_Preset_Thread_Initialize

  !# <nodeComponentThreadUninitializationTask>
  !#  <unitName>Node_Component_Position_Preset_Thread_Uninitialize</unitName>
  !# </nodeComponentThreadUninitializationTask>
  subroutine Node_Component_Position_Preset_Thread_Uninitialize()
    !% Uninitializes the tree node scale dark matter profile module.
    use :: Events_Hooks    , only : nodePromotionEvent
    use :: Galacticus_Nodes, only : defaultPositionComponent
    implicit none

    if (defaultPositionComponent%presetIsActive()) &
         & call nodePromotionEvent%detach(defaultPositionComponent,nodePromotion)
    return
  end subroutine Node_Component_Position_Preset_Thread_Uninitialize

  subroutine nodePromotion(self,node)
    !% Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, update the position of {\normalfont \ttfamily
    !% node} to that of the parent.
    use :: Galacticus_Nodes, only : nodeComponentPosition, nodeComponentPositionPreset, treeNode
    implicit none
    class(*                    ), intent(inout)          :: self
    type (treeNode             ), intent(inout), target  :: node
    class(nodeComponentPosition)               , pointer :: positionParent   , position, &
         &                                                  positionSatellite
    type (treeNode             )               , pointer :: nodeSatellite
    !$GLC attributes unused :: self
    
    position       => node       %position()
    positionParent => node%parent%position()
    select type (positionParent)
    class is (nodeComponentPositionPreset)
       call position%       positionSet(positionParent%position       ())
       call position%       velocitySet(positionParent%velocity       ())
       call position%positionHistorySet(positionParent%positionHistory())
    end select
    if (positionsPresetSatelliteToHost) then
       nodeSatellite => node%firstSatellite
       do while (associated(nodeSatellite))
          positionSatellite => nodeSatellite%position()
          call positionSatellite%positionSet(positionParent%position())
          call positionSatellite%velocitySet(positionParent%velocity())
          nodeSatellite => nodeSatellite%sibling
       end do
    end if
    return
  end subroutine nodePromotion

  !# <nodeMergerTask>
  !#  <unitName>Node_Component_Position_Preset_Move</unitName>
  !# </nodeMergerTask>
  !# <satelliteHostChangeTask>
  !#  <unitName>Node_Component_Position_Preset_Move</unitName>
  !# </satelliteHostChangeTask>
  subroutine Node_Component_Position_Preset_Move(node)
    !% Optionally move a satellite to coincide with the postion of its host.
    use :: Galacticus_Nodes, only : defaultPositionComponent, nodeComponentPosition, treeNode
    implicit none
    type (treeNode             ), intent(inout) :: node
    type (treeNode             ), pointer       :: nodeHost
    class(nodeComponentPosition), pointer       :: position, positionHost

    ! Return immediately if this method is not active or if positions are not to be reset.
    if (.not.defaultPositionComponent%presetIsActive() .or. .not.positionsPresetSatelliteToHost) return

    ! Move the satellite to the position of its host.
    nodeHost     => node    %parent
    position     => node    %position()
    positionHost => nodeHost%position()
    call position%positionSet(positionHost%position())
    call position%velocitySet(positionHost%velocity())
    return
  end subroutine Node_Component_Position_Preset_Move

  !# <interTreePositionInsert>
  !#  <unitName>Node_Component_Position_Preset_Inter_Tree_Insert</unitName>
  !# </interTreePositionInsert>
  subroutine Node_Component_Position_Preset_Inter_Tree_Insert(node,replaceNode)
    !% A satellite node is being moved between trees, and being added as a new satellite. Its (future-)histories will have been
    !% assigned to the {\normalfont \ttfamily replaceNode} so must be transferred.
    use :: Galacticus_Nodes, only : defaultPositionComponent, nodeComponentBasic, nodeComponentPosition, treeNode
    use :: Histories       , only : history
    implicit none
    type (treeNode             ), intent(inout), pointer :: node               , replaceNode
    class(nodeComponentPosition)               , pointer :: position           , replacePosition
    class(nodeComponentBasic   )               , pointer :: basic
    type (history              )                         :: historyPosition    , replaceHistoryPosition, &
         &                                                  moveHistoryPosition

    ! Return immediately if the preset position implementation is not active.
    if (.not.defaultPositionComponent%presetIsActive()) return
    ! Get the basic component of the pulled node.
    basic                  =>           node%basic          ()
    ! Get the position components to both nodes.
    position               =>           node%position       ()
    replacePosition        =>    replaceNode%position       ()
     ! Transfer subhalo mass history.
    historyPosition        =        position%positionHistory()
    replaceHistoryPosition = replacePosition%positionHistory()
    ! Cut off history in node being replaced subsequent to current time.
    call replaceHistoryPosition%trimForward       (basic%time(),   moveHistoryPosition)
    ! Append removed history to pulled node.
    call historyPosition       %append            (                moveHistoryPosition)
    ! Set the histories.
    call        position       %positionHistorySet(                    historyPosition)
    call replacePosition       %positionHistorySet(             replaceHistoryPosition)
    return
  end subroutine Node_Component_Position_Preset_Inter_Tree_Insert

end module Node_Component_Position_Preset
