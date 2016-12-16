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

module Node_Component_Position_Preset_Orphans
  !% Implements a preset position component with placement of orphan galaxies.
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Position_Preset_Orphans_Initialize

  !# <component>
  !#  <class>position</class>
  !#  <name>presetOrphans</name>
  !#  <extends>
  !#   <class>position</class>
  !#   <name>preset</name>
  !#  </extends>
  !#  <isDefault>false</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>position</name>
  !#     <type>double</type>
  !#     <rank>1</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <getFunction bindsTo="component">PositionPresetOrphansPosition</getFunction>
  !#     <output labels="[X,Y,Z]" unitsInSI="megaParsec" comment="Position of the node (in physical coordinates)."/>
  !#     <classDefault>[0.0d0,0.0d0,0.0d0]</classDefault>
  !#   </property>
  !#   <property>
  !#     <name>positionOrphan</name>
  !#     <type>double</type>
  !#     <rank>1</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" isDeferred="get" />
  !#   </property>
  !#   <property>
  !#     <name>timeAssign</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <classDefault>-1.0d0</classDefault>
  !#   </property>
  !#  </properties>
  !#  <functions>objects.nodes.components.position.preset.orphans.bound_functions.inc</functions>
  !# </component>

  ! Record of whether this module has been initialized.
  logical :: moduleInitialized=.false.

contains
  
  !# <nodeComponentInitializationTask>
  !#  <unitName>Node_Component_Position_Preset_Orphans_Initialize</unitName>
  !# </nodeComponentInitializationTask>
  subroutine Node_Component_Position_Preset_Orphans_Initialize()
    implicit none
    type(nodeComponentPositionPresetOrphans) :: position

    ! Initialize the module if necessary.
    !$omp critical (Node_Component_Position_Preset_Orphans_Initialize)
    if (defaultPositionComponent%presetIsActive().and..not.moduleInitialized) then
       ! Bind the position get function.
       call position%positionOrphanFunction(Node_Component_Position_Preset_Orphans_Position_Orphan)
       moduleInitialized=.true.
    end if
    !$omp end critical (Node_Component_Position_Preset_Orphans_Initialize)
    return
  end subroutine Node_Component_Position_Preset_Orphans_Initialize
  
  function Node_Component_Position_Preset_Orphans_Position_Orphan(self)
    !% Return the position of the orphan node.
    use Satellite_Oprhan_Distributions
    implicit none
    double precision                                    , allocatable  , dimension(:) :: Node_Component_Position_Preset_Orphans_Position_Orphan
    class           (nodeComponentPositionPresetOrphans), intent(inout)               :: self
    type            (treeNode                          ), pointer                     :: node
    class           (nodeComponentBasic                ), pointer                     :: basic
    class           (satelliteOrphanDistributionClass  ), pointer                     :: satelliteOrphanDistribution_

    node  => self%host ()
    basic => node%basic()
    if (basic%time() /= self%timeAssign()) then
       satelliteOrphanDistribution_ => satelliteOrphanDistribution()
       call self%    timeAssignSet(basic                       %time    (    ))
       call self%positionOrphanSet(satelliteOrphanDistribution_%position(node))
    end if
    Node_Component_Position_Preset_Orphans_Position_Orphan=self%positionOrphanValue()
    return
  end function Node_Component_Position_Preset_Orphans_Position_Orphan

end module Node_Component_Position_Preset_Orphans
