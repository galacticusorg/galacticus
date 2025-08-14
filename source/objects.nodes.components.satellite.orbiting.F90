!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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

!+    Contributions to this file made by:  Anthony Pullen, Andrew Benson.

!!{
Contains a module of satellite orbit tree node methods.
!!}
module Node_Component_Satellite_Orbiting
  !!{
  Implements the orbiting satellite component.
  !!}
  use :: Cosmology_Parameters   , only : cosmologyParametersClass
  use :: Cosmology_Functions    , only : cosmologyFunctionsClass
  use :: Virial_Density_Contrast, only : virialDensityContrastClass
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  use :: Kepler_Orbits          , only : keplerOrbit
  use :: Tensors                , only : tensorRank2Dimension3Symmetric
  use :: Virial_Orbits          , only : virialOrbit                    , virialOrbitClass
  implicit none
  private
  public :: Node_Component_Satellite_Orbiting_Scale_Set        , Node_Component_Satellite_Orbiting_Create             , &
       &    Node_Component_Satellite_Orbiting_Tree_Initialize  , Node_Component_Satellite_Orbiting_Initialize         , &
       &    Node_Component_Satellite_Orbiting_Thread_Initialize, Node_Component_Satellite_Orbiting_Thread_Uninitialize, &
       &    Node_Component_Satellite_Orbiting_State_Store      , Node_Component_Satellite_Orbiting_State_Restore
  
  !![
  <component>
   <class>satellite</class>
   <name>orbiting</name>
   <isDefault>false</isDefault>
   <createFunction isDeferred="true" />
   <properties>
    <property>
      <name>position</name>
      <type>double</type>
      <rank>1</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output labels="[X,Y,Z]" unitsInSI="megaParsec" comment="Orbital position of the node relative to its immediate host (sub-)halo."/>
      <classDefault>[0.0d0,0.0d0,0.0d0]</classDefault>
    </property>
    <property>
      <name>velocity</name>
      <type>double</type>
      <rank>1</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output labels="[X,Y,Z]" unitsInSI="kilo" comment="Orbital velocity of the node relative to its immediate host (sub-)halo."/>
      <classDefault>[0.0d0,0.0d0,0.0d0]</classDefault>
    </property>
    <property>
      <name>timeOfMerging</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
      <classDefault>huge(0.0d0)</classDefault>
    </property>
    <property>
      <name>destructionTime</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
      <classDefault>-1.0d0</classDefault>
    </property>
    <property>
      <name>boundMass</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <classDefault>selfBasic%mass()</classDefault>
      <output unitsInSI="massSolar" comment="Bound mass of the node."/>
    </property>
    <property>
      <name>virialOrbit</name>
      <type>keplerOrbit</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" isDeferred="set:get" />
    </property>
    <property>
      <name>tidalTensorPathIntegrated</name>
      <type>tensorRank2Dimension3Symmetric</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
    </property>
    <property>
      <name>tidalHeatingNormalized</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output unitsInSI="kilo**2/megaParsec**2" comment="Energy/radius^2 of satellite."/>
    </property>
   </properties>
  </component>
  !!]

  !![
  <enumeration>
   <name>initializationTypeMassBound</name>
   <description>Specify how to initialize the bound mass of a satellite halo.</description>
   <encodeFunction>yes</encodeFunction>
   <entry label="basicMass"      />
   <entry label="maximumRadius"  />
   <entry label="densityContrast"/>
  </enumeration>
  !!]

  ! Objects used by this module.
  class(virialDensityContrastClass), pointer :: virialDensityContrast_
  class(cosmologyParametersClass  ), pointer :: cosmologyParameters_
  class(cosmologyFunctionsClass   ), pointer :: cosmologyFunctions_
  class(darkMatterHaloScaleClass  ), pointer :: darkMatterHaloScale_
  class(virialOrbitClass          ), pointer :: virialOrbit_
  !$omp threadprivate(darkMatterHaloScale_,virialOrbit_,virialDensityContrast_,cosmologyParameters_,cosmologyFunctions_)

  ! Option controlling whether or not unbound virial orbits are acceptable.
  logical                                                      :: acceptUnboundOrbits

  ! Option controlling how to initialize the bound mass of satellite halos.
  type            (enumerationInitializationTypeMassBoundType) :: initializationTypeMassBound
  double precision                                             :: radiusMaximumOverRadiusVirial, densityContrastMassBound

  ! A threadprivate object used to track to which thread events are attached.
  integer :: thread
  !$omp threadprivate(thread)

contains

  !![
  <nodeComponentInitializationTask>
   <unitName>Node_Component_Satellite_Orbiting_Initialize</unitName>
  </nodeComponentInitializationTask>
  !!]
  subroutine Node_Component_Satellite_Orbiting_Initialize(parameters)
    !!{
    Initializes the orbiting satellite methods module.
    !!}
    use :: Error             , only : Error_Report
    use :: Galacticus_Nodes  , only : defaultSatelliteComponent, nodeComponentSatelliteOrbiting
    use :: ISO_Varying_String, only : char                     , var_str                       , varying_string
    use :: Input_Parameters  , only : inputParameter           , inputParameters
    implicit none
    type(inputParameters               ), intent(inout) :: parameters
    type(nodeComponentSatelliteOrbiting)                :: satellite
    type(varying_string                )                :: initializationTypeMassBoundText
    type(inputParameters               )                :: subParameters

    ! Initialize the module if necessary.
    if (defaultSatelliteComponent%orbitingIsActive()) then
       ! Find our parameters.
       subParameters=parameters%subParameters('componentSatellite')
       !![
       <inputParameter>
         <name>initializationTypeMassBound</name>
         <defaultValue>var_str('basicMass')</defaultValue>
         <description>Specify how to initialize the bound mass of a satellite halo. By default, the initial bound mass of a satellite halo is set to the node mass.</description>
         <source>subParameters</source>
         <variable>initializationTypeMassBoundText</variable>
       </inputParameter>
       !!]
       initializationTypeMassBound=enumerationInitializationTypeMassBoundEncode(char(initializationTypeMassBoundText),includesPrefix=.false.)
       !![
       <inputParameter>
         <name>radiusMaximumOverRadiusVirial</name>
         <defaultValue>1.0d0</defaultValue>
         <description>The maximum radius of the satellite halo in units of its virial radius. If {\normalfont \ttfamily [initializationTypeMassBound]} is set to 'maximumRadius', this value will be used to compute the initial bound mass of the satellite halo assuming that its density profile is 0 beyond this maximum radius.</description>
         <source>subParameters</source>
       </inputParameter>
       <inputParameter>
         <name>densityContrastMassBound</name>
         <defaultValue>200.0d0</defaultValue>
         <description>The density contrast of the satellite halo. If {\normalfont \ttfamily [initializationTypeMassBound]} is set to 'densityContrast', this value will be used to compute the initial bound mass of the satellite halo.</description>
         <source>subParameters</source>
       </inputParameter>
       <inputParameter>
         <name>acceptUnboundOrbits</name>
         <defaultValue>.false.</defaultValue>
         <description>If true, accept unbound virial orbits for satellites, otherwise reject them.</description>
         <source>subParameters</source>
       </inputParameter>
       !!]
       ! Validate the parameters.
       select case (initializationTypeMassBound%ID)
       case (initializationTypeMassBoundMaximumRadius  %ID)
          if (radiusMaximumOverRadiusVirial <= 0.0d0) then
             call Error_Report('specify a positive maximum radius for the satellite'//{introspection:location})
          end if
       case (initializationTypeMassBoundDensityContrast%ID)
          if (densityContrastMassBound               <= 0.0d0) then
             call Error_Report('specify a positive density contrast for the satellite'//{introspection:location})
          end if
       case default
          ! Do nothing.
       end select
       ! Specify the function to use for setting virial orbits.
       call satellite%virialOrbitSetFunction(Node_Component_Satellite_Orbiting_Virial_Orbit_Set)
       call satellite%virialOrbitFunction   (Node_Component_Satellite_Orbiting_Virial_Orbit    )
       ! Bind a creation function.
       call satellite%createFunctionSet     (Node_Component_Satellite_Orbiting_Initializor     )
    end if
    return
  end subroutine Node_Component_Satellite_Orbiting_Initialize

  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Node_Component_Satellite_Orbiting_Thread_Initialize</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Node_Component_Satellite_Orbiting_Thread_Initialize(parameters)
    !!{
    Initializes the tree node orbiting satellite module.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : defaultSatelliteComponent
    use :: Input_Parameters, only : inputParameter             , inputParameters
    use :: Events_Hooks    , only : satellitePreHostChangeEvent, nodePromotionEvent, openMPThreadBindingAtLevel, subhaloPromotionEvent, &
         &                          dependencyDirectionBefore  , dependencyExact

    implicit none
    type(inputParameters), intent(inout) :: parameters
    type(inputParameters)                :: subParameters
    type(dependencyExact), dimension(1)  :: dependenciesSubhaloPromotion

    if (defaultSatelliteComponent%orbitingIsActive()) then
       ! Find our parameters.
       subParameters=parameters%subParameters('componentSatellite')
       !![
       <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="subParameters"/>
       <objectBuilder class="cosmologyParameters"   name="cosmologyParameters_"   source="subParameters"/>
       <objectBuilder class="virialDensityContrast" name="virialDensityContrast_" source="subParameters"/>
       <objectBuilder class="darkMatterHaloScale"   name="darkMatterHaloScale_"   source="subParameters"/>
       <objectBuilder class="virialOrbit"           name="virialOrbit_"           source="subParameters"/>
       !!]
       dependenciesSubhaloPromotion(1)=dependencyExact(dependencyDirectionBefore,'mergerTreeNodeEvolver')
       call       subhaloPromotionEvent%attach(thread,subhaloPromotion      ,openMPThreadBindingAtLevel,label='nodeComponentSatelliteOrbiting',dependencies=dependenciesSubhaloPromotion)
       call          nodePromotionEvent%attach(thread,nodePromotion         ,openMPThreadBindingAtLevel,label='nodeComponentSatelliteOrbiting'                                          )
       call satellitePreHostChangeEvent%attach(thread,satellitePreHostChange,openMPThreadBindingAtLevel,label='nodeComponentSatelliteOrbiting'                                          )
       ! Check that the virial orbit class supports setting of angular coordinates.
       if (.not.virialOrbit_%isAngularlyResolved()) call Error_Report('"orbiting" satellite component requires a virialOrbit class which provides angularly-resolved orbits'//{introspection:location})
    end if
    return
  end subroutine Node_Component_Satellite_Orbiting_Thread_Initialize

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Node_Component_Satellite_Orbiting_Thread_Uninitialize</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Node_Component_Satellite_Orbiting_Thread_Uninitialize()
    !!{
    Uninitializes the tree node orbiting satellite module.
    !!}
    use :: Galacticus_Nodes, only : defaultSatelliteComponent
    use :: Events_Hooks    , only : satellitePreHostChangeEvent, nodePromotionEvent, subhaloPromotionEvent
    implicit none

    if (defaultSatelliteComponent%orbitingIsActive()) then
       !![
       <objectDestructor name="cosmologyFunctions_"   />
       <objectDestructor name="cosmologyParameters_"  />
       <objectDestructor name="virialDensityContrast_"/>
       <objectDestructor name="darkMatterHaloScale_"  />
       <objectDestructor name="virialOrbit_"          />
       !!]
       if (      subhaloPromotionEvent%isAttached(thread,subhaloPromotion      )) call       subhaloPromotionEvent%detach(thread,subhaloPromotion      )
       if (satellitePreHostChangeEvent%isAttached(thread,satellitePreHostChange)) call satellitePreHostChangeEvent%detach(thread,satellitePreHostChange)
       if (         nodePromotionEvent%isAttached(thread,nodePromotion         )) call          nodePromotionEvent%detach(thread,nodePromotion         )
   end if
    return
  end subroutine Node_Component_Satellite_Orbiting_Thread_Uninitialize

  !![
  <mergerTreeInitializeTask>
   <unitName>Node_Component_Satellite_Orbiting_Tree_Initialize</unitName>
  </mergerTreeInitializeTask>
  !!]
  subroutine Node_Component_Satellite_Orbiting_Tree_Initialize(node)
    !!{
    Initialize the orbiting satellite component.
    !!}
    use :: Galacticus_Nodes, only : treeNode, nodeComponentSatellite
    implicit none
    type (treeNode              ), pointer, intent(inout) :: node
    class(nodeComponentSatellite), pointer                :: satellite

    if (node%isSatellite()) then
       satellite => node%satellite()
       select type (satellite)
       type is (nodeComponentSatellite)
          call Node_Component_Satellite_Orbiting_Create(node)
       end select
    end if
    return
  end subroutine Node_Component_Satellite_Orbiting_Tree_Initialize

  !![
  <scaleSetTask>
   <unitName>Node_Component_Satellite_Orbiting_Scale_Set</unitName>
  </scaleSetTask>
  !!]
  subroutine Node_Component_Satellite_Orbiting_Scale_Set(node)
    !!{
    Set scales for properties of {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentSatellite, nodeComponentSatelliteOrbiting, nodeComponentBasic, treeNode
    use :: Numerical_Constants_Astronomical, only : gigaYear              , megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Tensors                         , only : tensorUnitR2D3Sym
    implicit none
    type            (treeNode              ), pointer  , intent(inout) :: node
    class           (nodeComponentBasic    ), pointer                  :: basic
    class           (nodeComponentSatellite), pointer                  :: satellite
    double precision                        , parameter                :: positionScaleFractional                 =1.0d-2                              , &
         &                                                                velocityScaleFractional                 =1.0d-2                              , &
         &                                                                boundMassScaleFractional                =1.0d-6                              , &
         &                                                                tidalTensorPathIntegratedScaleFractional=1.0d-2                              , &
         &                                                                tidalHeatingNormalizedScaleFractional   =1.0d-2
    double precision                                                   :: virialRadius                                   , virialVelocity              , &
         &                                                                virialIntegratedTidalTensor                    , virialTidalHeatingNormalized, &
         &                                                                massSatellite

    ! Get the satellite component.
    satellite => node%satellite()
    ! Ensure that it is of the orbiting class.
    select type (satellite)
    class is (nodeComponentSatelliteOrbiting)
       basic                       => node                %basic         (    )
       massSatellite               =  basic               %mass          (    )
       virialRadius                =  darkMatterHaloScale_%radiusVirial  (node)
       virialVelocity              =  darkMatterHaloScale_%velocityVirial(node)
       virialIntegratedTidalTensor =   virialVelocity/virialRadius*megaParsec/kilo/gigaYear
       virialTidalHeatingNormalized=  (virialVelocity/virialRadius)**2
       call satellite%positionScale                 (                                          &
            &                                        +[1.0d0,1.0d0,1.0d0]                      &
            &                                        *virialRadius                             &
            &                                        *                 positionScaleFractional &
            &                                       )
       call satellite%velocityScale                 (                                          &
            &                                        +[1.0d0,1.0d0,1.0d0]                      &
            &                                        *virialVelocity                           &
            &                                        *                 velocityScaleFractional &
            &                                       )
       call satellite%boundMassScale                (                                          &
            &                                        +massSatellite                            &
            &                                        *                boundMassScaleFractional &
            &                                       )
       call satellite%tidalTensorPathIntegratedScale(                                          &
            &                                        +tensorUnitR2D3Sym                        &
            &                                        *virialIntegratedTidalTensor              &
            &                                        *tidalTensorPathIntegratedScaleFractional &
            &                                       )
       call satellite%tidalHeatingNormalizedScale   (                                          &
            &                                        +virialTidalHeatingNormalized             &
            &                                        *   tidalHeatingNormalizedScaleFractional &
            &                                       )
    end select
    return
  end subroutine Node_Component_Satellite_Orbiting_Scale_Set

  subroutine Node_Component_Satellite_Orbiting_Initializor(self,timeEnd)
    !!{
    Initializes an orbiting satellite component.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatelliteOrbiting
    implicit none
    type            (nodeComponentSatelliteOrbiting), intent(inout)           :: self
    double precision                                , intent(in   ), optional :: timeEnd
    !$GLC attributes unused :: timeEnd

    call Node_Component_Satellite_Orbiting_Create(self%hostNode)
    return
  end subroutine Node_Component_Satellite_Orbiting_Initializor

  !![
  <nodeMergerTask>
   <unitName>Node_Component_Satellite_Orbiting_Create</unitName>
  </nodeMergerTask>
  !!]
  subroutine Node_Component_Satellite_Orbiting_Create(node)
    !!{
    Create a satellite orbit component and assign initial position, velocity, orbit, and tidal heating quantities. (The initial
    bound mass is automatically set to the original halo mass by virtue of that being the class default).
    !!}
    use :: Galacticus_Nodes, only : defaultSatelliteComponent, nodeComponentSatellite, nodeComponentSatelliteOrbiting, treeNode
    implicit none
    type            (treeNode              ), intent(inout) :: node
    type            (treeNode              ), pointer       :: nodeHost
    class           (nodeComponentSatellite), pointer       :: satellite
    logical                                                 :: isNewSatellite
    type            (keplerOrbit           )                :: orbit

    ! Return immediately if this method is not active.
    if (.not.defaultSatelliteComponent%orbitingIsActive()) return
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
       ! Set the initial bound mass of this satellite.
       call Node_Component_Satellite_Orbiting_Bound_Mass_Initialize(satellite,node)
    end if
    ! Create an orbit for the satellite if needed.
    select type (satellite)
    class is (nodeComponentSatelliteOrbiting)
       orbit=satellite%virialOrbitValue()
       if (.not.orbit%isDefined()) then
          ! Get an orbit for this satellite.
          if (node%isSatellite()) then
             nodeHost => node%parent
          else
             nodeHost => node%parent%firstChild
          end if
          orbit=virialOrbit_%orbit(node,nodeHost,acceptUnboundOrbits)          
          ! Store the orbit.
          call satellite%virialOrbitSet(orbit)          
       end if
    end select
    return
  end subroutine Node_Component_Satellite_Orbiting_Create
  
  subroutine nodePromotion(self,node)
    !!{
    Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply copy any preexisting satellite orbit
    from the parent.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : treeNode    , nodeComponentSatellite, nodeComponentSatelliteOrbiting
    implicit none
    class(*                     ), intent(inout)          :: self
    type (treeNode              ), intent(inout), target  :: node
    class(nodeComponentSatellite)               , pointer :: satellite, satelliteParent
    type  (keplerOrbit          )                         :: orbit    , orbitParent
    !$GLC attributes unused :: self

    satelliteParent => node%parent%satellite()
    select type (satelliteParent)
    type is (nodeComponentSatelliteOrbiting)
       satellite => node%satellite()
       select type (satellite)
       type is (nodeComponentSatellite)
          ! This is as expected - nothing to do. 
       type is (nodeComponentSatelliteOrbiting)
          orbitParent=satelliteParent%virialOrbit()
          orbit      =satellite      %virialOrbit()
          if (.not. orbit == orbitParent) call Error_Report('multiple satellite components defined on branch have differnt orbits'//{introspection:location})
       class default
          call Error_Report('multiple satellite components defined on branch'//{introspection:location})
       end select
       call node%parent%satelliteMove(node,overwrite=.true.)
    end select
    return
  end subroutine nodePromotion

  subroutine subhaloPromotion(self,node,nodePromotion)
    !!{
    Remove the satellite component from the subhalo about to be promoted to an isolated halo (which should have no satellite component).    
    !!}
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    class(*                     ), intent(inout)          :: self
    type (treeNode              ), intent(inout), pointer :: node, nodePromotion
     !$GLC attributes unused :: self, nodePromotion
    
    call node%satelliteRemove(1)
    return
  end subroutine subhaloPromotion

  subroutine satellitePreHostChange(self,node,nodeHostNew)
    !!{
    A satellite is about to move to a new host, adjust its position and velocity appropriately
    !!}
    use :: Galacticus_Nodes, only : defaultSatelliteComponent, nodeComponentSatellite, treeNode
    implicit none
    class           (*                     ), intent(inout)         :: self
    type            (treeNode              ), intent(inout), target :: node             , nodeHostNew
    type            (treeNode              ), pointer               :: nodeHost         , nodeHostNew_
    class           (nodeComponentSatellite), pointer               :: satellite        , satelliteHost
    double precision                        , dimension(3)          :: positionSatellite, velocitySatellite, &
         &                                                             positionHost     , velocityHost
    !$GLC attributes unused :: self
    
    ! Return immediately if this method is not active.
    if (.not.defaultSatelliteComponent%orbitingIsActive()) return
    ! Extract current position and velocity.
    satellite        =>  node     %satellite()
    positionSatellite =  satellite%position ()
    velocitySatellite =  satellite%velocity ()
    ! Walk up through hosts until the new host is found.
    nodeHost     => node       %parent
    nodeHostNew_ => nodeHostNew
    do while (associated(nodeHost))
       ! Extract position and velocity of this host.
       satelliteHost =>      nodeHost%satellite()
       positionHost  =  satelliteHost%position ()
       velocityHost  =  satelliteHost%velocity ()
       ! Shift the current node position and velocity by those of the host.
       positionSatellite=positionSatellite+positionHost
       velocitySatellite=velocitySatellite+velocityHost
       ! Move to the next host - if we arrive at the new host we are finished.
       nodeHost => nodeHost%parent
       if (associated(nodeHost,nodeHostNew_)) exit
    end do
    ! Update the position and velocity of the node.
    call satellite%positionSet(positionSatellite)
    call satellite%velocitySet(velocitySatellite)
    return
  end subroutine satellitePreHostChange
  
  function Node_Component_Satellite_Orbiting_Virial_Orbit(self) result(orbit)
    !!{
    Return the orbit of the satellite at the virial radius.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatelliteOrbiting, treeNode
    implicit none
    type (keplerOrbit                   )                :: orbit
    class(nodeComponentSatelliteOrbiting), intent(inout) :: self
    type (treeNode                      ), pointer       :: selfNode

    selfNode => self%host            ()
    orbit    =  self%virialOrbitValue()
    if (orbit%isDefined().or.selfNode%isSatellite().or.(.not.selfNode%isPrimaryProgenitor().and.associated(selfNode%parent))) then
       if (.not.orbit%isDefined()) then
          ! Orbit has not been defined - define it now.
          call Node_Component_Satellite_Orbiting_Create(selfNode)
          orbit=self%virialOrbitValue()
       end if
    else
       call orbit%reset()
    end if
    return
  end function Node_Component_Satellite_Orbiting_Virial_Orbit
  
  subroutine Node_Component_Satellite_Orbiting_Virial_Orbit_Set(self,orbit)
    !!{
    Set the orbit of the satellite at the virial radius.
    !!}
    use :: Coordinates     , only : assignment(=)
    use :: Galacticus_Nodes, only : nodeComponentSatellite, nodeComponentSatelliteOrbiting
    use :: Tensors         , only : tensorNullR2D3Sym
    implicit none
    class           (nodeComponentSatellite), intent(inout) :: self
    type            (keplerOrbit           ), intent(in   ) :: orbit
    double precision                        , dimension(3)  :: position   , velocity
    type            (keplerOrbit           )                :: virialOrbit

    select type (self)
    class is (nodeComponentSatelliteOrbiting)
       ! Ensure the orbit is defined.
       call orbit%assertIsDefined()
       ! Store the orbit.
       call self%virialOrbitSetValue(orbit)
       ! Store orbital position and velocity.
       virialOrbit=orbit
       position   =virialOrbit%position()
       velocity   =virialOrbit%velocity()
       call self%positionSet(position)
       call self%velocitySet(velocity)
       ! Set the merging/destruction time to -1 to indicate that we don't know when merging/destruction will occur.
       call self%destructionTimeSet          (           -1.0d0)
       call self%tidalTensorPathIntegratedSet(tensorNullR2D3Sym)
       call self%tidalHeatingNormalizedSet   (            0.0d0)       
    end select
    return
  end subroutine Node_Component_Satellite_Orbiting_Virial_Orbit_Set

  subroutine Node_Component_Satellite_Orbiting_Bound_Mass_Initialize(satellite,node)
    !!{
    Set the initial bound mass of the satellite.
    !!}
    use :: Mass_Distributions                  , only : massDistributionClass
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Error                               , only : Error_Report
    use :: Galacticus_Nodes                    , only : nodeComponentSatellite             , nodeComponentSatelliteOrbiting, treeNode
    implicit none
    class           (nodeComponentSatellite), intent(inout) :: satellite
    type            (treeNode              ), intent(inout) :: node
    class           (massDistributionClass ), pointer       :: massDistribution_
    double precision                                        :: massSatellite    , maximumRadius

    select type (satellite)
    class is (nodeComponentSatelliteOrbiting)
       select case (initializationTypeMassBound%ID)
       case (initializationTypeMassBoundBasicMass      %ID)
          ! Do nothing. The bound mass of this satellite is set to the node mass by default.
       case (initializationTypeMassBoundMaximumRadius  %ID)
          ! Set the initial bound mass of this satellite by integrating the density profile up to a maximum radius.
          maximumRadius     =  +                     radiusMaximumOverRadiusVirial                &
               &               *darkMatterHaloScale_%radiusVirial                 (node         )
          massDistribution_ =>  node                %massDistribution             (             )
          massSatellite     =   massDistribution_   %massEnclosedBySphere         (maximumRadius)
          call satellite%boundMassSet(massSatellite)
       case (initializationTypeMassBoundDensityContrast%ID)
          ! Set the initial bound mass of this satellite by assuming a specified density contrast.
          massSatellite=Dark_Matter_Profile_Mass_Definition(                                                 &
               &                                                                   node                    , &
               &                                                                   densityContrastMassBound, &
               &                                            cosmologyParameters_  =cosmologyParameters_    , &
               &                                            cosmologyFunctions_   =cosmologyFunctions_     , &
               &                                            virialDensityContrast_=virialDensityContrast_    &
               &                                           )
          call satellite%boundMassSet(massSatellite)
       case default
          call Error_Report('type of method to initialize the bound mass of satellites can not be recognized. Available options are "basicMass", "maximumRadius", "densityContrast"'//{introspection:location})
       end select
    end select
    return
  end subroutine Node_Component_Satellite_Orbiting_Bound_Mass_Initialize

  !![
  <stateStoreTask>
   <unitName>Node_Component_Satellite_Orbiting_State_Store</unitName>
  </stateStoreTask>
  !!]
  subroutine Node_Component_Satellite_Orbiting_State_Store(stateFile,gslStateFile,stateOperationID)
    !!{
    Store object state,
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Storing state for: componentSatellite -> orbiting',verbosity=verbosityLevelInfo)
    !![
    <stateStore variables="darkMatterHaloScale_ virialOrbit_ virialDensityContrast_ cosmologyParameters_ cosmologyFunctions_"/>
    !!]
    return
  end subroutine Node_Component_Satellite_Orbiting_State_Store

  !![
  <stateRetrieveTask>
   <unitName>Node_Component_Satellite_Orbiting_State_Restore</unitName>
  </stateRetrieveTask>
  !!]
  subroutine Node_Component_Satellite_Orbiting_State_Restore(stateFile,gslStateFile,stateOperationID)
    !!{
    Retrieve object state.
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Retrieving state for: componentSatellite -> orbiting',verbosity=verbosityLevelInfo)
    !![
    <stateRestore variables="darkMatterHaloScale_ virialOrbit_ virialDensityContrast_ cosmologyParameters_ cosmologyFunctions_"/>
    !!]
    return
  end subroutine Node_Component_Satellite_Orbiting_State_Restore

end module Node_Component_Satellite_Orbiting
