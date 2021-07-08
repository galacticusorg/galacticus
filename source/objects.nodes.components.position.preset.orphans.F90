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

!!{
Contains a module which implements a preset position component.
!!}

module Node_Component_Position_Preset_Orphans
  !!{
  Implements a preset position component with placement of orphan galaxies.
  !!}
  use :: Satellite_Oprhan_Distributions, only : satelliteOrphanDistributionClass
  implicit none
  private
  public :: Node_Component_Position_Preset_Orphans_Initialize         , Node_Component_Position_Preset_Orphans_Thread_Initialize, &
       &    Node_Component_Position_Preset_Orphans_Thread_Uninitialize

  !![
  <component>
   <class>position</class>
   <name>presetOrphans</name>
   <extends>
    <class>position</class>
    <name>preset</name>
   </extends>
   <isDefault>false</isDefault>
   <properties>
    <property>
      <name>position</name>
      <type>double</type>
      <rank>1</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
      <getFunction bindsTo="component">PositionPresetOrphansPosition</getFunction>
      <output labels="[X,Y,Z]" unitsInSI="megaParsec" comment="Position of the node (in physical coordinates)."/>
      <classDefault>[0.0d0,0.0d0,0.0d0]</classDefault>
    </property>
    <property>
      <name>positionOrphan</name>
      <type>double</type>
      <rank>1</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" isDeferred="get" />
    </property>
    <property>
      <name>velocity</name>
      <type>double</type>
      <rank>1</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
      <getFunction bindsTo="component">PositionPresetOrphansVelocity</getFunction>
      <output labels="[X,Y,Z]" unitsInSI="megaParsec" comment="Velocity of the node (peculiar, in physical coordinates)."/>
      <classDefault>[0.0d0,0.0d0,0.0d0]</classDefault>
    </property>
    <property>
      <name>velocityOrphan</name>
      <type>double</type>
      <rank>1</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" isDeferred="get" />
    </property>
    <property>
      <name>timeAssign</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
      <classDefault>-1.0d0</classDefault>
    </property>
   </properties>
   <functions>objects.nodes.components.position.preset.orphans.bound_functions.inc</functions>
  </component>
  !!]

  ! Objects used by this component.
  class(satelliteOrphanDistributionClass), pointer :: satelliteOrphanDistribution_
  !$omp threadprivate(satelliteOrphanDistribution_)

contains

  !![
  <nodeComponentInitializationTask>
   <unitName>Node_Component_Position_Preset_Orphans_Initialize</unitName>
  </nodeComponentInitializationTask>
  !!]
  subroutine Node_Component_Position_Preset_Orphans_Initialize(parameters_)
    use :: Galacticus_Nodes, only : defaultPositionComponent, nodeComponentPositionPresetOrphans
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(inputParameters                   ), intent(inout) :: parameters_
    type(nodeComponentPositionPresetOrphans)                :: position
    !$GLC attributes unused :: parameters_

    ! Initialize the module if necessary.
    if (defaultPositionComponent%presetIsActive()) then
       ! Bind the position get function.
       call position%positionOrphanFunction(Node_Component_Position_Preset_Orphans_Position_Orphan)
       call position%velocityOrphanFunction(Node_Component_Position_Preset_Orphans_Velocity_Orphan)
    end if
    return
  end subroutine Node_Component_Position_Preset_Orphans_Initialize

  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Node_Component_Position_Preset_Orphans_Thread_Initialize</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Node_Component_Position_Preset_Orphans_Thread_Initialize(parameters_)
    !!{
    Initializes the tree node preset orphans position module.
    !!}
    use :: Galacticus_Nodes, only : defaultPositionComponent
    use :: Input_Parameters, only : inputParameter          , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_

    if (defaultPositionComponent%presetOrphansIsActive()) then
       !![
       <objectBuilder class="satelliteOrphanDistribution" name="satelliteOrphanDistribution_" source="parameters_"/>
       !!]
    end if
    return
  end subroutine Node_Component_Position_Preset_Orphans_Thread_Initialize

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Node_Component_Position_Preset_Orphans_Thread_Uninitialize</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Node_Component_Position_Preset_Orphans_Thread_Uninitialize()
    !!{
    Uninitializes the tree node preset orphans position module.
    !!}
    use :: Galacticus_Nodes, only : defaultPositionComponent
    implicit none

    if (defaultPositionComponent%presetOrphansIsActive()) then
       !![
       <objectDestructor name="satelliteOrphanDistribution_"/>
       !!]
    end if
    return
  end subroutine Node_Component_Position_Preset_Orphans_Thread_Uninitialize

  function Node_Component_Position_Preset_Orphans_Position_Orphan(self)
    !!{
    Return the position of the orphan node.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentPositionPresetOrphans, treeNode
    implicit none
    double precision                                    , allocatable  , dimension(:) :: Node_Component_Position_Preset_Orphans_Position_Orphan
    class           (nodeComponentPositionPresetOrphans), intent(inout)               :: self
    type            (treeNode                          ), pointer                     :: node
    class           (nodeComponentBasic                ), pointer                     :: basic

    node  => self%host ()
    basic => node%basic()
    if (basic%time() /= self%timeAssign()) then
       call self%    timeAssignSet(basic                       %time    (    ))
       call self%positionOrphanSet(satelliteOrphanDistribution_%position(node))
       call self%velocityOrphanSet(satelliteOrphanDistribution_%velocity(node))
    end if
    Node_Component_Position_Preset_Orphans_Position_Orphan=self%positionOrphanValue()
    return
  end function Node_Component_Position_Preset_Orphans_Position_Orphan

  function Node_Component_Position_Preset_Orphans_Velocity_Orphan(self)
    !!{
    Return the velocity of the orphan node.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentPositionPresetOrphans, treeNode
    implicit none
    double precision                                    , allocatable  , dimension(:) :: Node_Component_Position_Preset_Orphans_Velocity_Orphan
    class           (nodeComponentPositionPresetOrphans), intent(inout)               :: self
    type            (treeNode                          ), pointer                     :: node
    class           (nodeComponentBasic                ), pointer                     :: basic

    node  => self%host ()
    basic => node%basic()
    if (basic%time() /= self%timeAssign()) then
       call self%    timeAssignSet(basic                       %time    (    ))
       call self%positionOrphanSet(satelliteOrphanDistribution_%position(node))
       call self%velocityOrphanSet(satelliteOrphanDistribution_%velocity(node))
    end if
    Node_Component_Position_Preset_Orphans_Velocity_Orphan=self%velocityOrphanValue()
    return
  end function Node_Component_Position_Preset_Orphans_Velocity_Orphan

end module Node_Component_Position_Preset_Orphans
