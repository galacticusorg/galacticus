!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

!% Contains a module of satellite orbit tree node methods.
module Node_Component_Satellite_Orbiting
  !% Implements the orbiting satellite component.
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  use :: Kepler_Orbits          , only : keplerOrbit
  use :: Tensors                , only : tensorRank2Dimension3Symmetric
  use :: Virial_Orbits          , only : virialOrbit                    , virialOrbitClass
  implicit none
  private
  public :: Node_Component_Satellite_Orbiting_Scale_Set        , Node_Component_Satellite_Orbiting_Create             , &
       &    Node_Component_Satellite_Orbiting_Pre_Host_Change  , Node_Component_Satellite_Orbiting_Initialize         , &
       &    Node_Component_Satellite_Orbiting_Thread_Initialize, Node_Component_Satellite_Orbiting_Thread_Uninitialize, &
       &    Node_Component_Satellite_Orbiting_Tree_Initialize
  
  !# <component>
  !#  <class>satellite</class>
  !#  <name>orbiting</name>
  !#  <isDefault>false</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>position</name>
  !#     <type>double</type>
  !#     <rank>1</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <output labels="[X,Y,Z]" unitsInSI="megaParsec" comment="Orbital position of the node."/>
  !#     <classDefault>[0.0d0,0.0d0,0.0d0]</classDefault>
  !#   </property>
  !#   <property>
  !#     <name>velocity</name>
  !#     <type>double</type>
  !#     <rank>1</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <output labels="[X,Y,Z]" unitsInSI="kilo" comment="Orbital velocity of the node."/>
  !#     <classDefault>[0.0d0,0.0d0,0.0d0]</classDefault>
  !#   </property>
  !#   <property>
  !#     <name>mergeTime</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <classDefault>-1.0d0</classDefault>
  !#   </property>
  !#   <property>
  !#     <name>timeOfMerging</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="false" isGettable="true" isEvolvable="false" isVirtual="true" />
  !#     <classDefault>-1.0d0</classDefault>
  !#     <getFunction>Node_Component_Satellite_Orbiting_Time_Of_Merging</getFunction>
  !#   </property>
  !#   <property>
  !#     <name>destructionTime</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <classDefault>-1.0d0</classDefault>
  !#   </property>
  !#   <property>
  !#     <name>boundMass</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <classDefault>selfBasicComponent%mass()</classDefault>
  !#     <output unitsInSI="massSolar" comment="Bound mass of the node."/>
  !#   </property>
  !#   <property>
  !#     <name>virialOrbit</name>
  !#     <type>keplerOrbit</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" isDeferred="set:get" />
  !#   </property>
  !#   <property>
  !#     <name>tidalTensorPathIntegrated</name>
  !#     <type>tensorRank2Dimension3Symmetric</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#   </property>
  !#   <property>
  !#     <name>tidalHeatingNormalized</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <output unitsInSI="kilo**2/megaParsec**2" comment="Energy/radius^2 of satellite."/>
  !#   </property>
  !#  </properties>
  !#  <functions>objects.nodes.components.satellite.orbiting.bound_functions.inc</functions>
  !# </component>

  !# <enumeration>
  !#  <name>satelliteBoundMassInitializeType</name>
  !#  <description>Specify how to initialize the bound mass of a satellite halo.</description>
  !#  <encodeFunction>yes</encodeFunction>
  !#  <entry label="basicMass"      />
  !#  <entry label="maximumRadius"  />
  !#  <entry label="densityContrast"/>
  !# </enumeration>

  ! Objects used by this module.
  class(darkMatterHaloScaleClass      ), pointer :: darkMatterHaloScale_
  class(virialOrbitClass              ), pointer :: virialOrbit_
  !$omp threadprivate(darkMatterHaloScale_,virialOrbit_)

  ! Option controlling whether or not unbound virial orbits are acceptable.
  logical         , parameter :: acceptUnboundOrbits=.false.

  ! Option controlling how to initialize the bound mass of satellite halos.
  integer                     :: satelliteBoundMassInitializeType
  double precision            :: satelliteMaximumRadiusOverVirialRadius      , satelliteDensityContrast

contains

  !# <nodeComponentInitializationTask>
  !#  <unitName>Node_Component_Satellite_Orbiting_Initialize</unitName>
  !# </nodeComponentInitializationTask>
  subroutine Node_Component_Satellite_Orbiting_Initialize(parameters_)
    !% Initializes the orbiting satellite methods module.
    use :: Galacticus_Error  , only : Galacticus_Error_Report
    use :: Galacticus_Nodes  , only : defaultSatelliteComponent, nodeComponentSatelliteOrbiting
    use :: ISO_Varying_String, only : char                     , var_str                       , varying_string
    use :: Input_Parameters  , only : inputParameter           , inputParameters
    implicit none
    type(inputParameters               ), intent(inout) :: parameters_
    type(nodeComponentSatelliteOrbiting)                :: satellite
    type(varying_string                )                :: satelliteBoundMassInitializeTypeText

    ! Initialize the module if necessary.
    if (defaultSatelliteComponent%orbitingIsActive()) then
       ! Create the spheroid mass distribution.
       !# <inputParameter>
       !#   <name>satelliteBoundMassInitializeType</name>
       !#   <defaultValue>var_str('basicMass')</defaultValue>
       !#   <description>Specify how to initialize the bound mass of a satellite halo. By default, the initial bound mass of a satellite halo is set to the node mass.</description>
       !#   <source>parameters_</source>
       !#   <variable>satelliteBoundMassInitializeTypeText</variable>
       !# </inputParameter>
       satelliteBoundMassInitializeType=enumerationSatelliteBoundMassInitializeTypeEncode(char(satelliteBoundMassInitializeTypeText),includesPrefix=.false.)
       !# <inputParameter>
       !#   <name>satelliteMaximumRadiusOverVirialRadius</name>
       !#   <defaultValue>1.0d0</defaultValue>
       !#   <description>The maximum radius of the satellite halo in units of its virial radius. If {\normalfont \ttfamily [satelliteBoundMassInitializeType]} is set to 'maximumRadius', this value will be used to compute the initial bound mass of the satellite halo assuming that its density profile is 0 beyond this maximum radius.</description>
       !#   <source>parameters_</source>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>satelliteDensityContrast</name>
       !#   <defaultValue>200.0d0</defaultValue>
       !#   <description>The density contrast of the satellite halo. If {\normalfont \ttfamily [satelliteBoundMassInitializeType]} is set to 'densityContrast', this value will be used to compute the initial bound mass of the satellite halo.</description>
       !#   <source>parameters_</source>
       !# </inputParameter>
       ! Validate the parameters.
       select case (satelliteBoundMassInitializeType)
       case (satelliteBoundMassInitializeTypeMaximumRadius  )
          if (satelliteMaximumRadiusOverVirialRadius <= 0.0d0) then
             call Galacticus_Error_Report('specify a positive maximum radius for the satellite'//{introspection:location})
          end if
       case (satelliteBoundMassInitializeTypeDensityContrast)
          if (satelliteDensityContrast               <= 0.0d0) then
             call Galacticus_Error_Report('specify a positive density contrast for the satellite'//{introspection:location})
          end if
       case default
          ! Do nothing.
       end select
       ! Specify the function to use for setting virial orbits.
       call satellite%virialOrbitSetFunction(Node_Component_Satellite_Orbiting_Virial_Orbit_Set)
       call satellite%virialOrbitFunction   (Node_Component_Satellite_Orbiting_Virial_Orbit    )
    end if
    return
  end subroutine Node_Component_Satellite_Orbiting_Initialize

  !# <nodeComponentThreadInitializationTask>
  !#  <unitName>Node_Component_Satellite_Orbiting_Thread_Initialize</unitName>
  !# </nodeComponentThreadInitializationTask>
  subroutine Node_Component_Satellite_Orbiting_Thread_Initialize(parameters_)
    !% Initializes the tree node orbiting satellite module.
    use :: Galacticus_Nodes, only : defaultSatelliteComponent
    use :: Input_Parameters, only : inputParameter           , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_

    if (defaultSatelliteComponent%orbitingIsActive()) then
       !# <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters_"/>
       !# <objectBuilder class="virialOrbit"         name="virialOrbit_"         source="parameters_"/>
    end if
    return
  end subroutine Node_Component_Satellite_Orbiting_Thread_Initialize

  !# <nodeComponentThreadUninitializationTask>
  !#  <unitName>Node_Component_Satellite_Orbiting_Thread_Uninitialize</unitName>
  !# </nodeComponentThreadUninitializationTask>
  subroutine Node_Component_Satellite_Orbiting_Thread_Uninitialize()
    !% Uninitializes the tree node orbiting satellite module.
    use :: Galacticus_Nodes, only : defaultSatelliteComponent
    implicit none

    if (defaultSatelliteComponent%orbitingIsActive()) then
       !# <objectDestructor name="darkMatterHaloScale_"/>
       !# <objectDestructor name="virialOrbit_"        />
    end if
    return
  end subroutine Node_Component_Satellite_Orbiting_Thread_Uninitialize

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Satellite_Orbiting_Tree_Initialize</unitName>
  !#  <after>darkMatterProfile</after>
  !#  <after>spin</after>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Satellite_Orbiting_Tree_Initialize(thisNode)
    !% Initialize the orbiting satellite component.
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode

    if (thisNode%isSatellite()) call Node_Component_Satellite_Orbiting_Create(thisNode)
    return
  end subroutine Node_Component_Satellite_Orbiting_Tree_Initialize

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Satellite_Orbiting_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Satellite_Orbiting_Scale_Set(node)
    !% Set scales for properties of {\normalfont \ttfamily node}.
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
         &                                                                satelliteMass

    ! Get the satellite component.
    satellite => node%satellite()
    ! Ensure that it is of the orbiting class.
    select type (satellite)
    class is (nodeComponentSatelliteOrbiting)
       basic                       => node                %basic         (    )
       satelliteMass               =  basic               %mass          (    )
       virialRadius                =  darkMatterHaloScale_%virialRadius  (node)
       virialVelocity              =  darkMatterHaloScale_%virialVelocity(node)
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
            &                                        +satelliteMass                            &
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

  !# <nodeMergerTask>
  !#  <unitName>Node_Component_Satellite_Orbiting_Create</unitName>
  !# </nodeMergerTask>
  subroutine Node_Component_Satellite_Orbiting_Create(node)
    !% Create a satellite orbit component and assign initial position, velocity, orbit, and tidal heating quantities. (The initial
    !% bound mass is automatically set to the original halo mass by virtue of that being the class default).
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
    ! Create an orbit for the satellite.   
    select type (satellite)
    class is (nodeComponentSatelliteOrbiting)
       ! Get an orbit for this satellite.
       if (node%isSatellite()) then
          nodeHost => node%parent
       else
          nodeHost => node%parent%firstChild
       end if
       orbit=virialOrbit_%orbit(node,nodeHost,acceptUnboundOrbits)
       ! Store the orbit.
       call satellite%virialOrbitSet(orbit)
    end select
    return
  end subroutine Node_Component_Satellite_Orbiting_Create
  
  !# <satellitePreHostChangeTask>
  !#  <unitName>Node_Component_Satellite_Orbiting_Pre_Host_Change</unitName>
  !# </satellitePreHostChangeTask>
  subroutine Node_Component_Satellite_Orbiting_Pre_Host_Change(node,nodeHostNew)
    !% A satellite is about to move to a new host, adjust its position and velocity appopriately
    use :: Galacticus_Nodes, only : defaultSatelliteComponent, nodeComponentSatellite, treeNode
    implicit none
    type            (treeNode              ), intent(inout), target :: node             , nodeHostNew
    type            (treeNode              ), pointer               :: nodeHost         , nodeHostNew_
    class           (nodeComponentSatellite), pointer               :: satellite        , satelliteHost
    double precision                        , dimension(3)          :: positionSatellite, velocitySatellite, &
         &                                                             positionHost     , velocityHost

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
  end subroutine Node_Component_Satellite_Orbiting_Pre_Host_Change
  
  function Node_Component_Satellite_Orbiting_Virial_Orbit(self) result(orbit)
    !% Return the orbit of the satellite at the virial radius.
    use :: Galacticus_Nodes, only : nodeComponentSatelliteOrbiting, treeNode
    implicit none
    type (keplerOrbit                   )                :: orbit
    class(nodeComponentSatelliteOrbiting), intent(inout) :: self
    type (treeNode                      ), pointer       :: selfNode

    selfNode => self%host()
    if (selfNode%isSatellite().or.(.not.selfNode%isPrimaryProgenitor().and.associated(selfNode%parent))) then
       orbit=self%virialOrbitValue()
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
  
  subroutine Node_Component_Satellite_Orbiting_Virial_Orbit_Set(self,thisOrbit)
    !% Set the orbit of the satellite at the virial radius.
    use :: Coordinates     , only : assignment(=)
    use :: Galacticus_Nodes, only : nodeComponentSatellite, nodeComponentSatelliteOrbiting
    use :: Tensors         , only : tensorNullR2D3Sym
    implicit none
    class           (nodeComponentSatellite), intent(inout) :: self
    type            (keplerOrbit           ), intent(in   ) :: thisOrbit
    double precision                        , dimension(3)  :: position   , velocity
    type            (keplerOrbit           )                :: virialOrbit

    select type (self)
    class is (nodeComponentSatelliteOrbiting)
       ! Ensure the orbit is defined.
       call thisOrbit%assertIsDefined()
       ! Store the orbit.
       call self%virialOrbitSetValue(thisOrbit)
       ! Store orbitial position and velocity.
       virialOrbit=thisOrbit
       position   =virialOrbit%position()
       velocity   =virialOrbit%velocity()
       call self%positionSet(position)
       call self%velocitySet(velocity)
       ! Set the merging/destruction time to -1 to indicate that we don't know when merging/destruction will occur.
       call self%mergeTimeSet                (           -1.0d0)
       call self%destructionTimeSet          (           -1.0d0)
       call self%tidalTensorPathIntegratedSet(tensorNullR2D3Sym)
       call self%tidalHeatingNormalizedSet   (            0.0d0)
    end select
    return
  end subroutine Node_Component_Satellite_Orbiting_Virial_Orbit_Set

  subroutine Node_Component_Satellite_Orbiting_Bound_Mass_Initialize(satelliteComponent,thisNode)
    !% Set the initial bound mass of the satellite.
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galactic_Structure_Enclosed_Masses  , only : Galactic_Structure_Enclosed_Mass
    use :: Galacticus_Error                    , only : Galacticus_Error_Report
    use :: Galacticus_Nodes                    , only : nodeComponentSatellite             , nodeComponentSatelliteOrbiting, treeNode
    implicit none
    class           (nodeComponentSatellite), intent(inout) :: satelliteComponent
    type            (treeNode              ), intent(inout) :: thisNode
    double precision                                        :: virialRadius      , maximumRadius, &
         &                                                     satelliteMass

    select type (satelliteComponent)
    class is (nodeComponentSatelliteOrbiting)
       select case (satelliteBoundMassInitializeType)
       case (satelliteBoundMassInitializeTypeBasicMass      )
          ! Do nothing. The bound mass of this satellite is set to the node mass by default.
       case (satelliteBoundMassInitializeTypeMaximumRadius  )
          ! Set the initial bound mass of this satellite by integrating the density profile up to a maximum radius.
          virialRadius =darkMatterHaloScale_%virialRadius  (thisNode                         )
          maximumRadius=satelliteMaximumRadiusOverVirialRadius*virialRadius
          satelliteMass=Galactic_Structure_Enclosed_Mass   (thisNode,maximumRadius           )
          call satelliteComponent%boundMassSet(satelliteMass)
       case (satelliteBoundMassInitializeTypeDensityContrast)
          ! Set the initial bound mass of this satellite by assuming a specified density contrast.
          satelliteMass=Dark_Matter_Profile_Mass_Definition(thisNode,satelliteDensityContrast)
          call satelliteComponent%boundMassSet(satelliteMass)
       case default
          call Galacticus_Error_Report('type of method to initialize the bound mass of satellites can not be recognized. Available options are "basicMass", "maximumRadius", "densityContrast"'//{introspection:location})
       end select
    end select
    return
  end subroutine Node_Component_Satellite_Orbiting_Bound_Mass_Initialize

end module Node_Component_Satellite_Orbiting
