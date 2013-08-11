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

!% Contains a module which implements a very simple satellite orbit component.

module Node_Component_Satellite_Very_Simple
  !% Implements a very simple satellite orbit component.
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Satellite_Very_Simple_Halo_Formation_Task, Node_Component_Satellite_Very_Simple_Create, &
       &    Node_Component_Satellite_Very_Simple_Tree_Initialize

  !# <component>
  !#  <class>satellite</class>
  !#  <name>verySimple</name>
  !#  <isDefault>no</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>mergeTime</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <classDefault>-1.0d0</classDefault>
  !#     <getFunction>Node_Component_Satellite_Very_Simple_Merge_Time</getFunction>
  !#     <output unitsInSI="gigaYear" comment="Time until satellite merges."/>
  !#   </property>
  !#   <property>
  !#     <name>timeOfMerging</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <isVirtual>true</isVirtual>
  !#     <attributes isSettable="false" isGettable="true" isEvolvable="false" />
  !#     <classDefault>-1.0d0</classDefault>
  !#     <getFunction>Node_Component_Satellite_Very_Simple_Time_Of_Merging</getFunction>
  !#   </property>
  !#  </properties>
  !#  <functions>objects.nodes.components.satellite.very_simple.bound_functions.inc</functions>
  !# </component>

  ! Flag indicating whether or not to reset satellite orbits on halo formation events.
  logical            :: satelliteOrbitResetOnHaloFormation

  ! Option controlling whether or not unbound virial orbits are acceptable.
  logical, parameter :: acceptUnboundOrbits               =.false.

  ! Record of whether this module has been initialized.
  logical            :: moduleInitialized                 =.false.

contains

  subroutine Node_Component_Satellite_Very_Simple_Initialize()
    !% Initializes the tree node satellite orbit methods module.
    use Input_Parameters
    implicit none

    ! Test whether module is already initialize.
    !$omp critical (Node_Component_Satellite_Very_Simple_Initialize)
    if (.not.moduleInitialized) then
       ! Determine if satellite orbits are to be reset on halo formation events.
       !@ <inputParameter>
       !@   <name>satelliteOrbitResetOnHaloFormation</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether satellite virial orbital parameters should be reset on halo formation events.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('satelliteOrbitResetOnHaloFormation',satelliteOrbitResetOnHaloFormation,defaultValue=.false.)
       ! Record that the module is now initialized.
       moduleInitialized=.true.
    end if
    !$omp end critical (Node_Component_Satellite_Very_Simple_Initialize)
    return
  end subroutine Node_Component_Satellite_Very_Simple_Initialize

  !# <haloFormationTask>
  !#  <unitName>Node_Component_Satellite_Very_Simple_Halo_Formation_Task</unitName>
  !# </haloFormationTask>
  subroutine Node_Component_Satellite_Very_Simple_Halo_Formation_Task(thisNode)
    !% Reset the orbits of satellite galaxies on halo formation events.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode
    type(treeNode)               , pointer :: satelliteNode

    ! Return immediately if this method is not active.
    if (.not.defaultSatelliteComponent%verySimpleIsActive()) return

    ! Return immediately if orbits are not to be reset.
    if (.not.satelliteOrbitResetOnHaloFormation) return

    ! Loop over all satellites.
    satelliteNode => thisNode%firstSatellite
    do while (associated(satelliteNode))
       ! Create a new orbit for this satellite.
       call Node_Component_Satellite_Very_Simple_Create(satelliteNode)
       satelliteNode => satelliteNode%sibling
    end do
    return
  end subroutine Node_Component_Satellite_Very_Simple_Halo_Formation_Task

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Satellite_Very_Simple_Tree_Initialize</unitName>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Satellite_Very_Simple_Tree_Initialize(thisNode)
    !% Initialize the very simple satellite component.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    if (thisNode%isSatellite()) call Node_Component_Satellite_Very_Simple_Create(thisNode)
    return
  end subroutine Node_Component_Satellite_Very_Simple_Tree_Initialize

  !# <nodeMergerTask>
  !#  <unitName>Node_Component_Satellite_Very_Simple_Create</unitName>
  !# </nodeMergerTask>
  !# <satelliteHostChangeTask>
  !#  <unitName>Node_Component_Satellite_Very_Simple_Create</unitName>
  !# </satelliteHostChangeTask>
  subroutine Node_Component_Satellite_Very_Simple_Create(thisNode)
    !% Create a satellite orbit component and assign a time until merging and a bound mass equal initially to the total halo mass.
    use Virial_Orbits
    use Satellite_Merging_Timescales
    use Kepler_Orbits
    implicit none
    type            (treeNode              ), intent(inout), pointer :: thisNode
    type            (treeNode              )               , pointer :: hostNode
    class           (nodeComponentSatellite)               , pointer :: satelliteComponent
    logical                                                          :: isNewSatellite
    double precision                                                 :: mergeTime
    type            (keplerOrbit           )                         :: thisOrbit

    ! Return immediately if this method is not active.
    if (.not.defaultSatelliteComponent%verySimpleIsActive()) return

    ! Get the satellite component.
    satelliteComponent => thisNode%satellite()
    ! Determine if the satellite component exists already.
    isNewSatellite=.false.
    select type (satelliteComponent)
    type is (nodeComponentSatellite)
       isNewSatellite=.true.
    end select

    ! If this is a new satellite, create the component.
    if (isNewSatellite) satelliteComponent => thisNode%satellite(autoCreate=.true.)

    select type (satelliteComponent)
    class is (nodeComponentSatelliteVerySimple)
       ! Ensure the module has been initialized.
       call Node_Component_Satellite_Very_Simple_Initialize()
       ! Get an orbit for this satellite.
       hostNode => thisNode%parent
       thisOrbit=Virial_Orbital_Parameters(thisNode,hostNode,acceptUnboundOrbits)
       ! Compute and store a time until merging.
       mergeTime=Satellite_Time_Until_Merging(thisNode,thisOrbit)
       if (mergeTime >= 0.0d0) call satelliteComponent%mergeTimeSet(mergeTime)
    end select
    return
  end subroutine Node_Component_Satellite_Very_Simple_Create

end module Node_Component_Satellite_Very_Simple
