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

!% Contains a module which implements a preset position component.

module Node_Component_Position_Preset
  !% Implements a preset position component.
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Position_Preset_Node_Promotion, Node_Component_Position_Preset_Initialize       , &
  &         Node_Component_Position_Preset_Move          , Node_Component_Position_Preset_Inter_Tree_Insert

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

  ! Record of whether this module has been initialized.
  logical :: moduleInitialized             =.false.

  ! Options.
  logical :: positionsPresetSatelliteToHost
  
contains
  
  !# <nodeComponentInitializationTask>
  !#  <unitName>Node_Component_Position_Preset_Initialize</unitName>
  !# </nodeComponentInitializationTask>
  subroutine Node_Component_Position_Preset_Initialize()
    use Input_Parameters
    implicit none
    
    ! Initialize the module if necessary.
    !$omp critical (Node_Component_Position_Preset_Initialize)
    if (defaultPositionComponent%presetIsActive().and..not.moduleInitialized) then
      ! Read parameters controlling the physical implementation.
       !@ <inputParameter>
       !@   <name>positionsPresetSatelliteToHost</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    If true, the position of satellite halos will be adjusted to match that of their host halo.
       !@   </description>
       !@   <type>bool</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('positionsPresetSatelliteToHost',positionsPresetSatelliteToHost,defaultValue=.false.)
       moduleInitialized=.true.
    end if
    !$omp end critical (Node_Component_Position_Preset_Initialize)
    return
  end subroutine Node_Component_Position_Preset_Initialize
  
  !# <nodePromotionTask>
  !#  <unitName>Node_Component_Position_Preset_Node_Promotion</unitName>
  !# </nodePromotionTask>
  subroutine Node_Component_Position_Preset_Node_Promotion(thisNode)
    !% Ensure that {\normalfont \ttfamily thisNode} is ready for promotion to its parent. In this case, update the position of {\tt
    !% thisNode} to that of the parent.
    implicit none
    type (treeNode             ), intent(inout), pointer :: thisNode
    class(nodeComponentPosition)               , pointer :: parentPositionComponent, thisPositionComponent, &
         &                                                  satellitePosition
    type (treeNode             )               , pointer :: satelliteNode

    thisPositionComponent => thisNode%position()
    select type (thisPositionComponent)
    class is (nodeComponentPositionPreset)
       parentPositionComponent => thisNode%parent%position()
       select type (parentPositionComponent)
       class is (nodeComponentPositionPreset)
          call thisPositionComponent%       positionSet(parentPositionComponent%position       ())
          call thisPositionComponent%       velocitySet(parentPositionComponent%velocity       ())
          call thisPositionComponent%positionHistorySet(parentPositionComponent%positionHistory())
       end select
       if (positionsPresetSatelliteToHost) then
          satelliteNode => thisNode%firstSatellite
          do while (associated(satelliteNode))
             satellitePosition => satelliteNode%position()
             call satellitePosition%positionSet(parentPositionComponent%position())
             call satellitePosition%velocitySet(parentPositionComponent%velocity())
             satelliteNode => satelliteNode%sibling
          end do
       end if
    end select
    return
  end subroutine Node_Component_Position_Preset_Node_Promotion
  
  !# <nodeMergerTask>
  !#  <unitName>Node_Component_Position_Preset_Move</unitName>
  !# </nodeMergerTask>
  !# <satelliteHostChangeTask>
  !#  <unitName>Node_Component_Position_Preset_Move</unitName>
  !# </satelliteHostChangeTask>
  subroutine Node_Component_Position_Preset_Move(node)
    !% Optinally move a satellite to coincide with the postion of its host.
    implicit none
    type (treeNode             ), intent(inout), pointer :: node
    type (treeNode             )               , pointer :: nodeHost
    class(nodeComponentPosition)               , pointer :: position, positionHost

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
    use Histories
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
