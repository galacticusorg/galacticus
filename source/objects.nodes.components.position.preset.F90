!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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
  public :: Node_Component_Position_Preset_Node_Promotion
  
  !# <component>
  !#  <class>position</class>
  !#  <name>preset</name>
  !#  <isDefault>no</isDefault>
  !#  <methods>
  !#   <method>
  !#     <name>position</name>
  !#     <type>real</type>
  !#     <rank>1</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <getFunction bindsTo="component">PositionPresetPosition</getFunction>
  !#     <output labels="[X,Y,Z]" unitsInSI="megaParsec" comment="Position of the node."/>
  !#     <classDefault>[0.0d0,0.0d0,0.0d0]</classDefault>
  !#   </method>
  !#   <method>
  !#     <name>velocity</name>
  !#     <type>real</type>
  !#     <rank>1</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <getFunction bindsTo="component">PositionPresetVelocity</getFunction>
  !#     <output labels="[X,Y,Z]" unitsInSI="kilo" comment="Velocity of the node."/>
  !#     <classDefault>[0.0d0,0.0d0,0.0d0]</classDefault>
  !#   </method>
  !#   <method>
  !#     <name>positionHistory</name>
  !#     <type>history</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </method>
  !#  </methods>
  !#  <functions>objects.nodes.components.position.preset.custom_methods.inc</functions>
  !# </component>

contains

  !# <nodePromotionTask>
  !#  <unitName>Node_Component_Position_Preset_Node_Promotion</unitName>
  !# </nodePromotionTask>
  subroutine Node_Component_Position_Preset_Node_Promotion(thisNode)
    !% Ensure that {\tt thisNode} is ready for promotion to its parent. In this case, update the position of {\tt
    !% thisNode} to that of the parent.
    implicit none
    type (treeNode             ), pointer, intent(inout) :: thisNode
    class(nodeComponentPosition), pointer                :: thisPositionComponent,parentPositionComponent

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
    end select
    return
  end subroutine Node_Component_Position_Preset_Node_Promotion
  
end module Node_Component_Position_Preset
