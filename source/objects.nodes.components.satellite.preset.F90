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

!+    Contributions to this file made by:  Andrew Benson, Jianling Gan.

!% Contains a module which implements a preset satellite orbit component.

module Node_Component_Satellite_Preset
  !% Implements a preset satellite orbit component.
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Satellite_Preset_Promote

  !# <component>
  !#  <class>satellite</class>
  !#  <name>preset</name>
  !#  <isDefault>no</isDefault>
  !#  <methods>
  !#   <method>
  !#     <name>mergeTime</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <isVirtual>yes</isVirtual>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <classDefault>-1.0d0</classDefault>
  !#     <getFunction>Node_Component_Satellite_Preset_Merge_Time</getFunction>
  !#     <setFunction>Node_Component_Satellite_Preset_Merge_Time_Set</setFunction>
  !#     <output unitsInSI="gigaYear" comment="Time until satellite merges."/>
  !#   </method>
  !#   <method>
  !#     <name>timeOfMerging</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <classDefault>-1.0d0</classDefault>
  !#   </method>
  !#   <method>
  !#     <name>boundMass</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <isVirtual>yes</isVirtual>
  !#     <attributes isSettable="false" isGettable="true" isEvolvable="false" />
  !#     <getFunction>SatellitePresetMergeBoundMass</getFunction>
  !#     <output unitsInSI="massSolar" comment="Bound mass of the node."/>
  !#   </method>
  !#   <method>
  !#     <name>boundMassHistory</name>
  !#     <type>history</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </method>
  !#   <method>
  !#     <name>virialOrbit</name>
  !#     <type>keplerOrbit</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </method>
  !#  </methods>
  !#  <functions>objects.nodes.components.satellite.preset.bound_functions.inc</functions>
  !# </component>

contains

  !# <nodePromotionTask>
  !#  <unitName>Node_Component_Satellite_Preset_Promote</unitName>
  !# </nodePromotionTask>
  subroutine Node_Component_Satellite_Preset_Promote(thisNode)
    !% Ensure that {\tt thisNode} is ready for promotion to its parent. In this case, we simply copy any preset satellite orbit
    !% from the parent.
    implicit none
    type (treeNode              ), pointer, intent(inout) :: thisNode
    type (treeNode              ), pointer                :: parentNode
    class(nodeComponentSatellite), pointer                :: thisSatelliteComponent,parentSatelliteComponent

    ! Return immediately if the preset satellite implementation is not active.
    if (.not.defaultSatelliteComponent%presetIsActive()) return
    ! Get the satellite component and check if it is of preset class.
    thisSatelliteComponent   => thisNode  %satellite(autoCreate=.true.)
    ! Get the parent node of this node.
    parentNode               => thisNode  %parent
    parentSatelliteComponent => parentNode%satellite(                 )
    ! Copy the satellite orbit from the parent node.
    select type (parentSatelliteComponent)
    class is (nodeComponentSatellitePreset)
       thisSatelliteComponent=parentSatelliteComponent
    end select
    return
  end subroutine Node_Component_Satellite_Preset_Promote

end module Node_Component_Satellite_Preset
