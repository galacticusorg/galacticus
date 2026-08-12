!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

!!{RST
Contains a module which implements a preset satellite orbit component.
!!}

module Node_Component_Satellite_Preset
  !!{RST
  Implements a preset satellite orbit component.
  !!}
  implicit none
  private
  public :: Node_Component_Satellite_Preset_Thread_Initialize, Node_Component_Satellite_Preset_Thread_Uninitialize

  !![
  <component docformat="rst">
   <class>satellite</class>
   <name>preset</name>
   <description>
   A satellite whose properties are read from the merger tree rather than computed. Bound mass and host node index are stored as time-series histories, from which the instantaneous values are interpolated, and the time of merging is likewise preset. ``isOrphan`` records whether the satellite has outlived its tracked history. Intended for use with merger trees read from N-body simulations, where the subhalo trajectory is already known.
   </description>
   <isDefault>false</isDefault>
   <properties>
    <property>
      <name>isOrphan</name>
      <type>logical</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
    </property>
    <property>
      <name>timeOfMerging</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
      <classDefault>huge(0.0d0)</classDefault>
    </property>
    <property>
      <name>boundMass</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="false" isGettable="true" isEvolvable="false" isVirtual="true" />
      <getFunction>SatellitePresetMergeBoundMass</getFunction>
      <classDefault>selfBasic%mass()</classDefault>
      <output unitsInSI="massSolar" unitsDescription="Solar masses" unitsQuantity="solMass" comment="Bound mass of the node."/>
    </property>
    <property>
      <name>boundMassHistory</name>
      <type>history</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
    </property>
    <property>
      <name>virialOrbit</name>
      <type>keplerOrbit</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
    </property>
    <property>
      <name>nodeIndexHistory</name>
      <type>longIntegerHistory</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
    </property>
    <property>
      <name>nodeIndex</name>
      <type>longInteger</type>
      <rank>0</rank>
      <attributes isSettable="false" isGettable="true" isEvolvable="false" />
      <getFunction>SatellitePresetNodeIndex</getFunction>
      <output unitsInSI="0.0d0" comment="Index of the satellite node used in preset satellite evolution."/>
      <classDefault>self%hostNode%index()</classDefault>
    </property>
   </properties>
   <functions>objects/nodes/components/satellite/preset/bound_functions.inc</functions>
  </component>
  !!]

  ! A threadprivate object used to track to which thread events are attached.
  integer :: thread
  !$omp threadprivate(thread)
  
contains

  !![
  <nodeComponentThreadInitializationTask function="Node_Component_Satellite_Preset_Thread_Initialize"/>
  !!]
  subroutine Node_Component_Satellite_Preset_Thread_Initialize(parameters_)
    !!{RST
    Initializes the tree node scale dark matter profile module.
    !!}
    use :: Events_Hooks    , only : nodePromotionEvent       , openMPThreadBindingAtLevel
    use :: Galacticus_Nodes, only : defaultSatelliteComponent
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_
    !$GLC attributes unused :: parameters_

    if (defaultSatelliteComponent%presetIsActive()) &
         & call nodePromotionEvent%attach(thread,nodePromotion,openMPThreadBindingAtLevel,label='nodeComponentSatellitePreset')
    return
  end subroutine Node_Component_Satellite_Preset_Thread_Initialize

  !![
  <nodeComponentThreadUninitializationTask function="Node_Component_Satellite_Preset_Thread_Uninitialize"/>
  !!]
  subroutine Node_Component_Satellite_Preset_Thread_Uninitialize()
    !!{RST
    Uninitializes the tree node scale dark matter profile module.
    !!}
    use :: Events_Hooks    , only : nodePromotionEvent
    use :: Galacticus_Nodes, only : defaultSatelliteComponent
    implicit none

    if (defaultSatelliteComponent%presetIsActive()) then
       if (nodePromotionEvent%isAttached(thread,nodePromotion)) call nodePromotionEvent%detach(thread,nodePromotion)
    end if
    return
  end subroutine Node_Component_Satellite_Preset_Thread_Uninitialize

  subroutine nodePromotion(self,node)
    !!{RST
    Ensure that ``node`` is ready for promotion to its parent. In this case, we simply copy any preset satellite orbit from the parent.
    !!}
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    class(*       ), intent(inout)         :: self
    type (treeNode), intent(inout), target :: node
    !$GLC attributes unused :: self
    
    ! Move the satellite orbit from the parent node.
    call node%parent%satelliteMove(node,overwrite=.true.)
    return
  end subroutine nodePromotion

end module Node_Component_Satellite_Preset
