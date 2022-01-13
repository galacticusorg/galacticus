!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
Contains a module which implements a very simple satellite orbit component.
!!}

module Node_Component_Satellite_Very_Simple
  !!{
  Implements a very simple satellite orbit component.
  !!}
  use :: Satellite_Merging_Timescales, only : satelliteMergingTimescalesClass
  use :: Virial_Orbits               , only : virialOrbitClass
  implicit none
  private
  public :: Node_Component_Satellite_Very_Simple_Halo_Formation_Task, Node_Component_Satellite_Very_Simple_Create             , &
       &    Node_Component_Satellite_Very_Simple_Tree_Initialize    , Node_Component_Satellite_Very_Simple_Rate_Compute       , &
       &    Node_Component_Satellite_Very_Simple_Thread_Initialize  , Node_Component_Satellite_Very_Simple_Thread_Uninitialize, &
       &    Node_Component_Satellite_Very_Simple_State_Store        , Node_Component_Satellite_Very_Simple_State_Restore      , &
       &    Node_Component_Satellite_Very_Simple_Initialize

  !![
  <component>
   <class>satellite</class>
   <name>verySimple</name>
   <isDefault>false</isDefault>
   <properties>
    <property>
      <name>mergeTime</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <classDefault>-1.0d0</classDefault>
      <getFunction>Node_Component_Satellite_Very_Simple_Merge_Time</getFunction>
      <output unitsInSI="gigaYear" comment="Time until satellite merges."/>
    </property>
    <property>
      <name>timeOfMerging</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="false" isGettable="true" isEvolvable="false" isVirtual="true" />
      <classDefault>-1.0d0</classDefault>
      <getFunction>Node_Component_Satellite_Very_Simple_Time_Of_Merging</getFunction>
    </property>
   </properties>
   <functions>objects.nodes.components.satellite.very_simple.bound_functions.inc</functions>
  </component>
  !!]

  ! Objects used by this component.
  class(virialOrbitClass               ), pointer :: virialOrbit_
  class(satelliteMergingTimescalesClass), pointer :: satelliteMergingTimescales_
  !$omp threadprivate(virialOrbit_,satelliteMergingTimescales_)

  ! Flag indicating whether or not to reset satellite orbits on halo formation events.
  logical            :: satelliteOrbitResetOnHaloFormation

  ! Option controlling whether or not unbound virial orbits are acceptable.
  logical, parameter :: acceptUnboundOrbits               =.false.

contains

  !![
  <nodeComponentInitializationTask>
   <unitName>Node_Component_Satellite_Very_Simple_Initialize</unitName>
  </nodeComponentInitializationTask>
  !!]
  subroutine Node_Component_Satellite_Very_Simple_Initialize(parameters_)
    !!{
    Initializes the tree node satellite orbit methods module.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_

    ! Determine if satellite orbits are to be reset on halo formation events.
    !![
    <inputParameter>
      <name>satelliteOrbitResetOnHaloFormation</name>
      <defaultValue>.false.</defaultValue>
      <description>Specifies whether satellite virial orbital parameters should be reset on halo formation events.</description>
      <source>parameters_</source>
    </inputParameter>
    !!]
    return
  end subroutine Node_Component_Satellite_Very_Simple_Initialize

  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Node_Component_Satellite_Very_Simple_Thread_Initialize</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Node_Component_Satellite_Very_Simple_Thread_Initialize(parameters_)
    !!{
    Initializes the tree node very simple satellite module.
    !!}
    use :: Galacticus_Nodes, only : defaultSatelliteComponent
    use :: Input_Parameters, only : inputParameter           , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_

    if (defaultSatelliteComponent%verySimpleIsActive()) then
       !![
       <objectBuilder class="virialOrbit"                name="virialOrbit_"                source="parameters_"/>
       <objectBuilder class="satelliteMergingTimescales" name="satelliteMergingTimescales_" source="parameters_"/>
       !!]
    end if
    return
  end subroutine Node_Component_Satellite_Very_Simple_Thread_Initialize

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Node_Component_Satellite_Very_Simple_Thread_Uninitialize</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Node_Component_Satellite_Very_Simple_Thread_Uninitialize()
    !!{
    Uninitializes the tree node very simple satellite module.
    !!}
    use :: Galacticus_Nodes, only : defaultSatelliteComponent
    implicit none

    if (defaultSatelliteComponent%verySimpleIsActive()) then
       !![
       <objectDestructor name="virialOrbit_"               />
       <objectDestructor name="satelliteMergingTimescales_"/>
       !!]
    end if
    return
  end subroutine Node_Component_Satellite_Very_Simple_Thread_Uninitialize

  !![
  <haloFormationTask>
   <unitName>Node_Component_Satellite_Very_Simple_Halo_Formation_Task</unitName>
  </haloFormationTask>
  !!]
  subroutine Node_Component_Satellite_Very_Simple_Halo_Formation_Task(node)
    !!{
    Reset the orbits of satellite galaxies on halo formation events.
    !!}
    use :: Galacticus_Nodes, only : defaultSatelliteComponent, treeNode
    implicit none
    type(treeNode), intent(inout) :: node
    type(treeNode), pointer       :: satelliteNode

    ! Return immediately if this method is not active.
    if (.not.defaultSatelliteComponent%verySimpleIsActive()) return

    ! Return immediately if orbits are not to be reset.
    if (.not.satelliteOrbitResetOnHaloFormation) return

    ! Loop over all satellites.
    satelliteNode => node%firstSatellite
    do while (associated(satelliteNode))
       ! Create a new orbit for this satellite.
       call Node_Component_Satellite_Very_Simple_Create(satelliteNode)
       satelliteNode => satelliteNode%sibling
    end do
    return
  end subroutine Node_Component_Satellite_Very_Simple_Halo_Formation_Task

  !![
  <rateComputeTask>
   <unitName>Node_Component_Satellite_Very_Simple_Rate_Compute</unitName>
  </rateComputeTask>
  !!]
  subroutine Node_Component_Satellite_Very_Simple_Rate_Compute(node,interrupt,interruptProcedure,propertyType)
    !!{
    Compute the time until satellite merging rate of change.
    !!}
    use :: Galacticus_Nodes, only : defaultSatelliteComponent, nodeComponentSatellite, nodeComponentSatelliteVerySimple, propertyTypeInactive, &
          &                         treeNode
    implicit none
    type            (treeNode              ), intent(inout)          :: node
    logical                                 , intent(inout)          :: interrupt
    procedure       (                      ), intent(inout), pointer :: interruptProcedure
    integer                                 , intent(in   )          :: propertyType
    class           (nodeComponentSatellite)               , pointer :: satellite
    !$GLC attributes unused :: interrupt, interruptProcedure

    ! Return immediately if inactive variables are requested.
    if (propertyType == propertyTypeInactive) return
    ! Return immediately if this class is not in use.
    if (.not.defaultSatelliteComponent%verySimpleIsActive()) return
    ! Get the satellite component.
    satellite => node%satellite()
    ! Ensure that it is of the standard class.
    select type (satellite)
    class is (nodeComponentSatelliteVerySimple)
       if (node%isSatellite()) call satellite%mergeTimeRate(-1.0d0)
    end select
    return
  end subroutine Node_Component_Satellite_Very_Simple_Rate_Compute

  !![
  <mergerTreeInitializeTask>
   <unitName>Node_Component_Satellite_Very_Simple_Tree_Initialize</unitName>
   <after>darkMatterProfile</after>
  </mergerTreeInitializeTask>
  !!]
  subroutine Node_Component_Satellite_Very_Simple_Tree_Initialize(node)
    !!{
    Initialize the very simple satellite component.
    !!}
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    type(treeNode), intent(inout), pointer :: node

    if (node%isSatellite()) call Node_Component_Satellite_Very_Simple_Create(node)
    return
  end subroutine Node_Component_Satellite_Very_Simple_Tree_Initialize

  !![
  <nodeMergerTask>
   <unitName>Node_Component_Satellite_Very_Simple_Create</unitName>
  </nodeMergerTask>
  <satelliteHostChangeTask>
   <unitName>Node_Component_Satellite_Very_Simple_Create</unitName>
  </satelliteHostChangeTask>
  !!]
  subroutine Node_Component_Satellite_Very_Simple_Create(node)
    !!{
    Create a satellite orbit component and assign a time until merging and a bound mass equal initially to the total halo mass.
    !!}
    use :: Galacticus_Nodes, only : defaultSatelliteComponent, nodeComponentSatellite, nodeComponentSatelliteVerySimple, treeNode
    use :: Kepler_Orbits   , only : keplerOrbit
    implicit none
    type            (treeNode              ), intent(inout) :: node
    type            (treeNode              ), pointer       :: nodeHost
    class           (nodeComponentSatellite), pointer       :: satellite
    logical                                                 :: isNewSatellite
    double precision                                        :: mergeTime
    type            (keplerOrbit           )                :: orbit

    ! Return immediately if this method is not active.
    if (.not.defaultSatelliteComponent%verySimpleIsActive()) return

    ! Get the satellite component.
    satellite => node%satellite()
    ! Determine if the satellite component exists already.
    isNewSatellite=.false.
    select type (satellite)
    type is (nodeComponentSatellite)
       isNewSatellite=.true.
    end select

    ! If this is a new satellite, create the component.
    if (isNewSatellite) satellite => node%satellite(autoCreate=.true.)

    select type (satellite)
    class is (nodeComponentSatelliteVerySimple)
       ! Get an orbit for this satellite.
       if (node%isSatellite()) then
          nodeHost => node%parent
       else
          nodeHost => node%parent%firstChild
       end if
       orbit=virialOrbit_%orbit(node,nodeHost,acceptUnboundOrbits)
       ! Compute and store a time until merging.
       mergeTime=satelliteMergingTimescales_%timeUntilMerging(node,orbit)
       if (mergeTime >= 0.0d0) call satellite%mergeTimeSet(mergeTime)
    end select
    return
  end subroutine Node_Component_Satellite_Very_Simple_Create

  !![
  <galacticusStateStoreTask>
   <unitName>Node_Component_Satellite_Very_Simple_State_Store</unitName>
  </galacticusStateStoreTask>
  !!]
  subroutine Node_Component_Satellite_Very_Simple_State_Store(stateFile,gslStateFile,stateOperationID)
    !!{
    Store object state,
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Storing state for: componentSatellite -> verySimple',verbosity=verbosityLevelInfo)
    !![
    <stateStore variables="virialOrbit_ satelliteMergingTimescales_"/>
    !!]
    return
  end subroutine Node_Component_Satellite_Very_Simple_State_Store

  !![
  <galacticusStateRetrieveTask>
   <unitName>Node_Component_Satellite_Very_Simple_State_Restore</unitName>
  </galacticusStateRetrieveTask>
  !!]
  subroutine Node_Component_Satellite_Very_Simple_State_Restore(stateFile,gslStateFile,stateOperationID)
    !!{
    Retrieve object state.
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Retrieving state for: componentSatellite -> verySimple',verbosity=verbosityLevelInfo)
    !![
    <stateRestore variables="virialOrbit_ satelliteMergingTimescales_"/>
    !!]
    return
  end subroutine Node_Component_Satellite_Very_Simple_State_Restore

end module Node_Component_Satellite_Very_Simple
