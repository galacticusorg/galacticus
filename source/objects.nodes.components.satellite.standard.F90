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
Contains a module of satellite orbit tree node methods.
!!}
module Node_Component_Satellite_Standard
  !!{
  Implements the standard satellite component.
  !!}
  use :: Dark_Matter_Halos_Mass_Loss_Rates, only : darkMatterHaloMassLossRateClass
  use :: Kepler_Orbits                    , only : keplerOrbit
  use :: Satellite_Merging_Timescales     , only : satelliteMergingTimescalesClass
  use :: Virial_Orbits                    , only : virialOrbit                    , virialOrbitClass
  implicit none
  private
  public :: Node_Component_Satellite_Standard_Scale_Set          , Node_Component_Satellite_Standard_Create           , &
       &    Node_Component_Satellite_Standard_Rate_Compute       , Node_Component_Satellite_Standard_Initialize       , &
       &    Node_Component_Satellite_Standard_Halo_Formation_Task, Node_Component_Satellite_Standard_Tree_Initialize  , &
       &    Node_Component_Satellite_Standard_Inactive           , Node_Component_Satellite_Standard_Thread_Initialize, &
       &    Node_Component_Satellite_Standard_Thread_Uninitialize, Node_Component_Satellite_Standard_State_Store      , &
       &    Node_Component_Satellite_Standard_State_Restore

  !![
  <component>
   <class>satellite</class>
   <name>standard</name>
   <isDefault>true</isDefault>
   <properties>
    <property>
      <name>mergeTime</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <classDefault>-1.0d0</classDefault>
      <getFunction>Node_Component_Satellite_Standard_Merge_Time</getFunction>
      <output unitsInSI="gigaYear" comment="Time until satellite merges."/>
    </property>
    <property>
      <name>timeOfMerging</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" isVirtual="true" />
      <classDefault>-1.0d0</classDefault>
      <getFunction>Node_Component_Satellite_Standard_Time_Of_Merging</getFunction>
      <setFunction>Node_Component_Satellite_Standard_Time_Of_Merging_Set</setFunction>
    </property>
    <property>
      <name>boundMass</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <classDefault>selfBasicComponent%mass()</classDefault>
      <output unitsInSI="massSolar" comment="Bound mass of the node."/>
    </property>
    <property>
      <name>virialOrbit</name>
      <type>keplerOrbit</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" isDeferred="set:get" />
    </property>
   </properties>
   <functions>objects.nodes.components.satellite.standard.bound_functions.inc</functions>
  </component>
  !!]

  ! Objects used by this component.
  class(darkMatterHaloMassLossRateClass), pointer :: darkMatterHaloMassLossRate_
  class(virialOrbitClass               ), pointer :: virialOrbit_
  class(satelliteMergingTimescalesClass), pointer :: satelliteMergingTimescales_
  !$omp threadprivate(darkMatterHaloMassLossRate_,virialOrbit_,satelliteMergingTimescales_)

  ! Option indicating whether or not satellite virial orbital parameters will be stored.
  logical            :: satelliteOrbitStoreOrbitalParameters

  ! Option indicating whether or not to reset satellite orbits on halo formation events.
  logical            :: satelliteOrbitResetOnHaloFormation

  ! Record of whether satellite bound mass is an inactive variable.
  logical            :: satelliteBoundMassIsInactive

  ! Option controlling whether or not unbound virial orbits are acceptable.
  logical, parameter :: acceptUnboundOrbits                 =.false.

contains

  !![
  <nodeComponentInitializationTask>
   <unitName>Node_Component_Satellite_Standard_Initialize</unitName>
  </nodeComponentInitializationTask>
  !!]
   subroutine Node_Component_Satellite_Standard_Initialize(parameters_)
     !!{
     Initializes the standard satellite orbit component module.
     !!}
     use :: Galacticus_Nodes, only : nodeComponentSatelliteStandard
     use :: Input_Parameters, only : inputParameter                , inputParameters
     implicit none
     type(inputParameters               ), intent(inout) :: parameters_
     type(nodeComponentSatelliteStandard)                :: satellite

     if (satellite%standardIsActive()) then
        ! Determine if satellite orbits are to be stored.
        !![
        <inputParameter>
          <name>satelliteOrbitStoreOrbitalParameters</name>
          <defaultValue>.true.</defaultValue>
          <description>Specifies whether satellite virial orbital parameters should be stored (otherwise they are computed
             again---possibly at random---each time they are requested).</description>
          <source>parameters_</source>
        </inputParameter>
        !!]
        ! Determine if satellite orbits are to be reset on halo formation events.
        !![
        <inputParameter>
          <name>satelliteOrbitResetOnHaloFormation</name>
          <defaultValue>.false.</defaultValue>
          <description>Specifies whether satellite virial orbital parameters should be reset on halo formation events.</description>
          <source>parameters_</source>
        </inputParameter>
        !!]
        ! Determine if bound mass is an inactive variable.
        !![
        <inputParameter>
          <name>satelliteBoundMassIsInactive</name>
          <defaultValue>.false.</defaultValue>
          <description>Specifies whether or not the bound mass variable of the standard satellite component is inactive (i.e. does not appear in any ODE being solved).</description>
          <source>parameters_</source>
        </inputParameter>
        !!]
         ! Specify the function to use for setting virial orbits.
        call satellite%virialOrbitSetFunction(Node_Component_Satellite_Standard_Virial_Orbit_Set)
        call satellite%virialOrbitFunction   (Node_Component_Satellite_Standard_Virial_Orbit    )
     end if
     return
   end subroutine Node_Component_Satellite_Standard_Initialize

   !![
   <nodeComponentThreadInitializationTask>
    <unitName>Node_Component_Satellite_Standard_Thread_Initialize</unitName>
   </nodeComponentThreadInitializationTask>
   !!]
   subroutine Node_Component_Satellite_Standard_Thread_Initialize(parameters_)
     !!{
     Initializes the tree node standard satellite module.
     !!}
     use :: Galacticus_Nodes, only : defaultSatelliteComponent
     use :: Input_Parameters, only : inputParameter           , inputParameters
     implicit none
     type(inputParameters), intent(inout) :: parameters_

     if (defaultSatelliteComponent%standardIsActive()) then
        !![
        <objectBuilder class="darkMatterHaloMassLossRate" name="darkMatterHaloMassLossRate_" source="parameters_"/>
        <objectBuilder class="virialOrbit"                name="virialOrbit_"                source="parameters_"/>
        <objectBuilder class="satelliteMergingTimescales" name="satelliteMergingTimescales_" source="parameters_"/>
        !!]
     end if
     return
   end subroutine Node_Component_Satellite_Standard_Thread_Initialize

   !![
   <nodeComponentThreadUninitializationTask>
    <unitName>Node_Component_Satellite_Standard_Thread_Uninitialize</unitName>
   </nodeComponentThreadUninitializationTask>
   !!]
   subroutine Node_Component_Satellite_Standard_Thread_Uninitialize()
     !!{
     Uninitializes the tree node standard satellite module.
     !!}
     use :: Galacticus_Nodes, only : defaultSatelliteComponent
     implicit none

     if (defaultSatelliteComponent%standardIsActive()) then
        !![
        <objectDestructor name="darkMatterHaloMassLossRate_"/>
        <objectDestructor name="virialOrbit_"               />
        <objectDestructor name="satelliteMergingTimescales_"/>
        !!]
     end if
     return
   end subroutine Node_Component_Satellite_Standard_Thread_Uninitialize

  !![
  <inactiveSetTask>
   <unitName>Node_Component_Satellite_Standard_Inactive</unitName>
  </inactiveSetTask>
  !!]
  subroutine Node_Component_Satellite_Standard_Inactive(node)
    !!{
    Set Jacobian zero status for properties of {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite, nodeComponentSatelliteStandard, treeNode
    implicit none
    type (treeNode              ), intent(inout), pointer :: node
    class(nodeComponentSatellite)               , pointer :: satellite

    ! Get the satellite component.
    satellite => node%satellite()
    ! Check if an standard satellite component exists.
    select type (satellite)
    class is (nodeComponentSatelliteStandard)
       if (satelliteBoundMassIsInactive) call satellite%boundMassInactive()
    end select
    return
  end subroutine Node_Component_Satellite_Standard_Inactive

  !![
  <mergerTreeInitializeTask>
   <unitName>Node_Component_Satellite_Standard_Tree_Initialize</unitName>
  </mergerTreeInitializeTask>
  !!]
  subroutine Node_Component_Satellite_Standard_Tree_Initialize(node)
    !!{
    Initialize the standard satellite component.
    !!}
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    type(treeNode), intent(inout), pointer :: node

    if     (                                                                           &
         &                     node                            %isSatellite        ()  &
         &  .or.                                                                       &
         &   (                                                                         &
         &     .not.           node                            %isPrimaryProgenitor()  &
         &    .and.                                                                    &
         &          associated(node                            %parent               ) &
         &    .and.                                                                    &
         &                     satelliteOrbitStoreOrbitalParameters                    &
         &   )                                                                         &
         & ) call Node_Component_Satellite_Standard_Create(node)
    return
  end subroutine Node_Component_Satellite_Standard_Tree_Initialize

  !![
  <rateComputeTask>
   <unitName>Node_Component_Satellite_Standard_Rate_Compute</unitName>
  </rateComputeTask>
  !!]
  subroutine Node_Component_Satellite_Standard_Rate_Compute(node,interrupt,interruptProcedure,propertyType)
    !!{
    Compute the time until satellite merging rate of change.
    !!}
    use :: Galacticus_Nodes, only : defaultSatelliteComponent, nodeComponentSatellite, nodeComponentSatelliteStandard, propertyEvaluate, &
          &                         treeNode
    implicit none
    type            (treeNode              ), intent(inout)          :: node
    logical                                 , intent(inout)          :: interrupt
    procedure       (                      ), intent(inout), pointer :: interruptProcedure
    integer                                 , intent(in   )          :: propertyType
    class           (nodeComponentSatellite)               , pointer :: satellite
    double precision                                                 :: massLossRate
    !$GLC attributes unused :: interrupt, interruptProcedure

    ! Return immediately if this class is not in use.
    if (.not.defaultSatelliteComponent%standardIsActive()) return
    ! Get the satellite component.
    satellite => node%satellite()
    ! Ensure that it is of the standard class.
    select type (satellite)
    class is (nodeComponentSatelliteStandard)
       if (node%isSatellite()) then
          call satellite%mergeTimeRate(-1.0d0)
          ! Compute mass loss rate if necessary.
          if (propertyEvaluate(propertyType,satelliteBoundMassIsInactive)) then
             massLossRate=darkMatterHaloMassLossRate_%rate(node)
             call satellite%boundMassRate(massLossRate)
          end if
       end if
    end select
    return
  end subroutine Node_Component_Satellite_Standard_Rate_Compute

  function Node_Component_Satellite_Standard_Virial_Orbit(self) result(orbit)
    !!{
    Return the orbit of the satellite at the virial radius.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatelliteStandard, treeNode
    implicit none
    type (keplerOrbit                   )                :: orbit
    class(nodeComponentSatelliteStandard), intent(inout) :: self
    type (treeNode                      ), pointer       :: hostNode    , selfNode

    selfNode => self%host()
    if (selfNode%isSatellite().or.(.not.selfNode%isPrimaryProgenitor().and.associated(selfNode%parent))) then
       if (satelliteOrbitStoreOrbitalParameters) then
          orbit=self%virialOrbitValue()
          if (.not.orbit%isDefined()) then
             ! Orbit has not been defined - define it now.
             call Node_Component_Satellite_Standard_Create(selfNode)
             orbit=self%virialOrbitValue()
          end if
       else
          if (selfNode%isSatellite()) then
             hostNode => selfNode%parent
          else
             hostNode => selfNode%parent%firstChild
          end if
          orbit=virialOrbit_%orbit(selfNode,hostNode,acceptUnboundOrbits)
       end if
    else
       call orbit%reset()
    end if
    return
  end function Node_Component_Satellite_Standard_Virial_Orbit

  subroutine Node_Component_Satellite_Standard_Virial_Orbit_Set(self,orbit)
    !!{
    Set the orbit of the satellite at the virial radius.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite, nodeComponentSatelliteStandard, treeNode
    implicit none
    class           (nodeComponentSatellite         ), intent(inout) :: self
    type            (keplerOrbit                    ), intent(in   ) :: orbit
    type            (treeNode                       ), pointer       :: selfNode
    double precision                                                 :: mergeTime
    type            (keplerOrbit                    )                :: virialOrbit

    select type (self)
    class is (nodeComponentSatelliteStandard)
       ! Ensure the orbit is defined.
       call orbit%assertIsDefined()
       ! Get the node.
       selfNode => self%host()
       ! Update the stored time until merging to reflect the new orbit.
       virialOrbit=orbit
       mergeTime  =satelliteMergingTimescales_%timeUntilMerging(selfNode,virialOrbit)
       if (mergeTime >= 0.0d0) call self%mergeTimeSet(mergeTime)
       ! Store the orbit.
       call self%virialOrbitSetValue(orbit)
    end select
    return
  end subroutine Node_Component_Satellite_Standard_Virial_Orbit_Set

  !![
  <scaleSetTask>
   <unitName>Node_Component_Satellite_Standard_Scale_Set</unitName>
  </scaleSetTask>
  !!]
  subroutine Node_Component_Satellite_Standard_Scale_Set(node)
    !!{
    Set scales for properties of {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentSatellite, nodeComponentSatelliteStandard, treeNode
    implicit none
    type            (treeNode              ), intent(inout), pointer :: node
    class           (nodeComponentSatellite)               , pointer :: satellite
    class           (nodeComponentBasic    )               , pointer :: basic
    double precision                        , parameter              :: massScaleFractional=1.0d-6, timeScale=1.0d-3

    ! Get the satellite component.
    satellite => node%satellite()
    ! Ensure that it is of the standard class.
    select type (satellite)
    class is (nodeComponentSatelliteStandard)
       ! Get the basic component.
       basic => node%basic()
       ! Set scale for time.
       call satellite%mergeTimeScale(timeScale                       )
       ! Set scale for bound mass.
       call satellite%boundMassScale(massScaleFractional*basic%mass())
    end select
    return
  end subroutine Node_Component_Satellite_Standard_Scale_Set

  !![
  <haloFormationTask>
   <unitName>Node_Component_Satellite_Standard_Halo_Formation_Task</unitName>
  </haloFormationTask>
  !!]
  subroutine Node_Component_Satellite_Standard_Halo_Formation_Task(node)
    !!{
    Reset the orbits of satellite galaxies on halo formation events.
    !!}
    use :: Galacticus_Nodes, only : defaultSatelliteComponent, treeNode
    implicit none
    type(treeNode), intent(inout) :: node
    type(treeNode), pointer       :: satelliteNode

    ! Return immediately if this method is not active.
    if (.not.defaultSatelliteComponent%standardIsActive()) return

    ! Return immediately if orbits are not to be reset.
    if (.not.satelliteOrbitResetOnHaloFormation) return

    ! Loop over all satellites.
    satelliteNode => node%firstSatellite
    do while (associated(satelliteNode))
       ! Create a new orbit for this satellite.
       call Node_Component_Satellite_Standard_Create(satelliteNode)
       satelliteNode => satelliteNode%sibling
    end do

    return
  end subroutine Node_Component_Satellite_Standard_Halo_Formation_Task

  !![
  <nodeMergerTask>
   <unitName>Node_Component_Satellite_Standard_Create</unitName>
  </nodeMergerTask>
  <satelliteHostChangeTask>
   <unitName>Node_Component_Satellite_Standard_Create</unitName>
  </satelliteHostChangeTask>
  !!]
  subroutine Node_Component_Satellite_Standard_Create(node)
    !!{
    Create a satellite orbit component and assign a time until merging and a bound mass equal initially to the total halo mass.
    !!}
    use :: Galacticus_Nodes, only : defaultSatelliteComponent, nodeComponentBasic           , nodeComponentSatellite, nodeComponentSatelliteStandard, &
         &                          treeNode
    use :: Kepler_Orbits   , only : keplerOrbitMasses        , keplerOrbitRadius            , keplerOrbitTheta      , keplerOrbitPhi                , &
         &                          keplerOrbitVelocityRadial, keplerOrbitVelocityTangential
    implicit none
    type            (treeNode              ), intent(inout) :: node
    type            (treeNode              ), pointer       :: hostNode
    class           (nodeComponentSatellite), pointer       :: satellite
    class           (nodeComponentBasic    ), pointer       :: basic
    logical                                                 :: isNewSatellite
    double precision                                        :: mergeTime
    type            (keplerOrbit           )                :: orbit

    ! Return immediately if this method is not active.
    if (.not.defaultSatelliteComponent%standardIsActive()) return

    ! Get the satellite component.
    satellite => node%satellite()
    ! Determine if the satellite component exists already.
    isNewSatellite=.false.
    select type (satellite)
    type is (nodeComponentSatellite)
       isNewSatellite=.true.
    end select

    ! If this is a new satellite, create the component and set the bound mass.
    if (isNewSatellite) then
       satellite => node%satellite(autoCreate=.true.)
       select type (satellite)
       class is (nodeComponentSatelliteStandard)
          ! Set the bound mass of the satellite.
          basic => node%basic()
          call satellite%boundMassSet(basic%mass())
       end select
    end if

    select type (satellite)
    class is (nodeComponentSatelliteStandard)
       ! Get an orbit for this satellite.
       if (node%isSatellite()) then
          hostNode => node%parent
       else
          hostNode => node%parent%firstChild
       end if
       ! Test if the orbit is already defined. (This can happen if some other component requested the orbit during its initialization.)
       orbit=satellite%virialOrbitValue()
       if (orbit%isDefined()) then
          ! The orbit has been previous defined, reset derived properties.
          call orbit%reset(keep=[keplerOrbitMasses,keplerOrbitRadius,keplerOrbitTheta,keplerOrbitPhi,keplerOrbitVelocityRadial,keplerOrbitVelocityTangential])
       else
          ! The orbit is not previously defined, so choose an orbit now.
          orbit=virialOrbit_%orbit(node,hostNode,acceptUnboundOrbits)
       end if
       ! Store the orbit if necessary.
       if (satelliteOrbitStoreOrbitalParameters) call satellite%virialOrbitSet(orbit)
       ! Compute and store a time until merging.
       mergeTime=satelliteMergingTimescales_%timeUntilMerging(node,orbit)
       if (mergeTime >= 0.0d0) call satellite%mergeTimeSet(mergeTime)
    end select
    return
  end subroutine Node_Component_Satellite_Standard_Create

  !![
  <stateStoreTask>
   <unitName>Node_Component_Satellite_Standard_State_Store</unitName>
  </stateStoreTask>
  !!]
  subroutine Node_Component_Satellite_Standard_State_Store(stateFile,gslStateFile,stateOperationID)
    !!{
    Store object state,
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Storing state for: componentSatellite -> standard',verbosity=verbosityLevelInfo)
    !![
    <stateStore variables="darkMatterHaloMassLossRate_ virialOrbit_ satelliteMergingTimescales_"/>
    !!]
    return
  end subroutine Node_Component_Satellite_Standard_State_Store

  !![
  <stateRetrieveTask>
   <unitName>Node_Component_Satellite_Standard_State_Restore</unitName>
  </stateRetrieveTask>
  !!]
  subroutine Node_Component_Satellite_Standard_State_Restore(stateFile,gslStateFile,stateOperationID)
    !!{
    Retrieve object state.
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Retrieving state for: componentSatellite -> standard',verbosity=verbosityLevelInfo)
    !![
    <stateRestore variables="darkMatterHaloMassLossRate_ virialOrbit_ satelliteMergingTimescales_"/>
    !!]
    return
  end subroutine Node_Component_Satellite_Standard_State_Restore

end module Node_Component_Satellite_Standard
