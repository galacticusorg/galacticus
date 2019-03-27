!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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
  use Kepler_Orbits
  use Tensors
  use Satellite_Dynamical_Friction
  use Satellite_Tidal_Stripping
  use Satellite_Tidal_Heating
  use Dark_Matter_Halo_Scales
  use Virial_Orbits
  use ISO_Varying_String          , only : varying_string, var_str, operator(==)
  implicit none
  private
  public :: Node_Component_Satellite_Orbiting_Scale_Set        , Node_Component_Satellite_Orbiting_Create             , &
       &    Node_Component_Satellite_Orbiting_Rate_Compute     , Node_Component_Satellite_Orbiting_Tree_Initialize    , &
       &    Node_Component_Satellite_Orbiting_Trigger_Merger   , Node_Component_Satellite_Orbiting_Initialize         , &
       &    Node_Component_Satellite_Orbiting_Thread_Initialize, Node_Component_Satellite_Orbiting_Thread_Uninitialize

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
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" isDeferred="set" />
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
  
  ! Objects used by this module.
  class(darkMatterHaloScaleClass       ), pointer :: darkMatterHaloScale_
  class(satelliteDynamicalFrictionClass), pointer :: satelliteDynamicalFriction_
  class(satelliteTidalHeatingRateClass ), pointer :: satelliteTidalHeatingRate_
  class(satelliteTidalStrippingClass   ), pointer :: satelliteTidalStripping_
  class(virialOrbitClass               ), pointer :: virialOrbit_
  !$omp threadprivate(darkMatterHaloScale_,satelliteDynamicalFriction_,satelliteTidalHeatingRate_,satelliteTidalStripping_,virialOrbit_)

  ! Option controlling whether or not unbound virial orbits are acceptable.
  logical         , parameter :: acceptUnboundOrbits=.false.

  ! Option controlling minimum mass of satellite halos before a merger is triggered.
  double precision            :: satelliteOrbitingDestructionMass
  logical                     :: satelliteOrbitingDestructionMassIsFractional
  ! Option controlling how to initialize the bound mass of satellite halos.
  type(varying_string)        :: satelliteBoundMassInitializeType
  double precision            :: satelliteMaximumRadiusOverVirialRadius      , satelliteDensityContrast

contains

  !# <mergerTreePreTreeConstructionTask>
  !#  <unitName>Node_Component_Satellite_Orbiting_Initialize</unitName>
  !# </mergerTreePreTreeConstructionTask>
  subroutine Node_Component_Satellite_Orbiting_Initialize()
    !% Initializes the orbiting satellite methods module.
    use Input_Parameters
    use Galacticus_Nodes, only : nodeComponentSatelliteOrbiting, defaultSatelliteComponent
    implicit none
    type(nodeComponentSatelliteOrbiting) :: satelliteComponent

    ! Initialize the module if necessary.
    if (defaultSatelliteComponent%orbitingIsActive()) then
       ! Create the spheroid mass distribution.
       !# <inputParameter>
       !#   <name>satelliteOrbitingDestructionMassIsFractional</name>
       !#   <defaultValue>.true.</defaultValue>
       !#   <description>If true, then {\normalfont \ttfamily [satelliteOrbitingDestructionMass]} specifies the fractional mass a halo must reach before it is tidally destroyed. Otherwise, {\normalfont \ttfamily [satelliteOrbitingDestructionMass]} specifies an absolute mass.</description>
       !#   <type>double</type>
       !#   <cardinality>1</cardinality>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>satelliteOrbitingDestructionMass</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultValue>0.01d0</defaultValue>
       !#   <description>The mass (possibly fractional---see {\normalfont \ttfamily [satelliteOrbitingDestructionMassIsFractional]}) below which the satellite is considered to be tidally destroyed and merged with the central halo.</description>
       !#   <source>globalParameters</source>
       !#   <type>double</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>satelliteBoundMassInitializeType</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultValue>var_str('basicMass')</defaultValue>
       !#   <description>Specify how to initialize the bound mass of a satellite halo. By default, the initial bound mass of a satellite halo is set to the node mass.</description>
       !#   <source>globalParameters</source>
       !#   <type>string</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>satelliteMaximumRadiusOverVirialRadius</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultValue>1.0d0</defaultValue>
       !#   <description>The maximum radius of the satellite halo in units of its virial radius. If {\normalfont \ttfamily [satelliteBoundMassInitializeType]} is set to 'maximumRadius', this value will be used to compute the initial bound mass of the satellite halo assuming that its density profile is 0 beyond this maximum radius.</description>
       !#   <source>globalParameters</source>
       !#   <type>double</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>satelliteDensityContrast</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultValue>200.0d0</defaultValue>
       !#   <description>The density contrast of the satellite halo. If {\normalfont \ttfamily [satelliteBoundMassInitializeType]} is set to 'densityContrast', this value will be used to compute the initial bound mass of the satellite halo.</description>
       !#   <source>globalParameters</source>
       !#   <type>double</type>
       !# </inputParameter>
       ! Specify the function to use for setting virial orbits.
       call satelliteComponent%virialOrbitSetFunction(Node_Component_Satellite_Orbiting_Virial_Orbit_Set)
    end if
    return
  end subroutine Node_Component_Satellite_Orbiting_Initialize

  !# <nodeComponentThreadInitializationTask>
  !#  <unitName>Node_Component_Satellite_Orbiting_Thread_Initialize</unitName>
  !# </nodeComponentThreadInitializationTask>
  subroutine Node_Component_Satellite_Orbiting_Thread_Initialize(globalParameters_)
    !% Initializes the tree node orbiting satellite module.
    use Galacticus_Nodes, only : defaultSatelliteComponent
    use Input_Parameters
    implicit none
    type(inputParameters), intent(inout) :: globalParameters_

    if (defaultSatelliteComponent%orbitingIsActive()) then
       !# <objectBuilder class="darkMatterHaloScale"        name="darkMatterHaloScale_"        source="globalParameters_"/>
       !# <objectBuilder class="satelliteDynamicalFriction" name="satelliteDynamicalFriction_" source="globalParameters_"/>
       !# <objectBuilder class="satelliteTidalHeatingRate"  name="satelliteTidalHeatingRate_"  source="globalParameters_"/>
       !# <objectBuilder class="satelliteTidalStripping"    name="satelliteTidalStripping_"    source="globalParameters_"/>
       !# <objectBuilder class="virialOrbit"                name="virialOrbit_"                source="globalParameters_"/>
    end if
    return
  end subroutine Node_Component_Satellite_Orbiting_Thread_Initialize

  !# <nodeComponentThreadUninitializationTask>
  !#  <unitName>Node_Component_Satellite_Orbiting_Thread_Uninitialize</unitName>
  !# </nodeComponentThreadUninitializationTask>
  subroutine Node_Component_Satellite_Orbiting_Thread_Uninitialize()
    !% Uninitializes the tree node orbiting satellite module.
    use Galacticus_Nodes, only : defaultSatelliteComponent
    implicit none

    if (defaultSatelliteComponent%orbitingIsActive()) then
       !# <objectDestructor name="darkMatterHaloScale_"       />
       !# <objectDestructor name="satelliteDynamicalFriction_"/>
       !# <objectDestructor name="satelliteTidalHeatingRate_" />
       !# <objectDestructor name="satelliteTidalStripping_"   />
       !# <objectDestructor name="virialOrbit_"               />
    end if
    return
  end subroutine Node_Component_Satellite_Orbiting_Thread_Uninitialize

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Satellite_Orbiting_Tree_Initialize</unitName>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Satellite_Orbiting_Tree_Initialize(thisNode)
    !% Initialize the orbiting satellite component.
    use Galacticus_Nodes, only : treeNode
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode

    if (thisNode%isSatellite()) call Node_Component_Satellite_Orbiting_Create(thisNode)
    return
  end subroutine Node_Component_Satellite_Orbiting_Tree_Initialize
  
  !# <rateComputeTask>
  !#  <unitName>Node_Component_Satellite_Orbiting_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Satellite_Orbiting_Rate_Compute(thisNode,odeConverged,interrupt,interruptProcedure,propertyType)
    !% Compute rate of change for satellite properties.
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Physical
    use Numerical_Constants_Math
    use Galactic_Structure_Enclosed_Masses
    use Vectors
    use Galactic_Structure_Densities
    use Galactic_Structure_Options
    use Galacticus_Nodes                  , only : treeNode            , nodeComponentBasic, nodeComponentSatellite   , nodeComponentSatelliteOrbiting, &
         &                                         propertyTypeInactive, interruptTask     , defaultSatelliteComponent
    implicit none
    type            (treeNode                       ), pointer     , intent(inout) :: thisNode
    logical                                                        , intent(in   ) :: odeConverged
    logical                                                        , intent(inout) :: interrupt
    procedure       (interruptTask                  ), pointer     , intent(inout) :: interruptProcedure
    integer                                                        , intent(in   ) :: propertyType
    class           (nodeComponentSatellite         ), pointer                     :: satelliteComponent
    class           (nodeComponentBasic             ), pointer                     :: basicComponent
    type            (treeNode                       ), pointer                     :: hostNode
    double precision                                 , dimension(3)                :: position,velocity
    double precision                                 , parameter                   :: radiusVirialFraction = 1.0d-2
    double precision                                                               :: radius,halfMassRadiusSatellite
    double precision                                                               :: halfMassRadiusCentral,orbitalRadiusTest
    double precision                                                               :: radiusVirial
    double precision                                                               :: orbitalPeriod
    double precision                                                               :: parentDensity
    double precision                                                               :: parentEnclosedMass,satelliteMass,basicMass,massDestruction
    double precision                                                               :: tidalHeatingNormalized,angularFrequency,radialFrequency
    type            (tensorRank2Dimension3Symmetric )                              :: tidalTensor,tidalTensorPathIntegrated,positionTensor

    ! Return immediately if inactive variables are requested.
    if (propertyType == propertyTypeInactive) return
    ! Return immediately if this class is not in use.
    if (.not.defaultSatelliteComponent%orbitingIsActive()) return
    ! Get the satellite component.
    satelliteComponent => thisNode%satellite()
    ! Ensure that it is of the orbiting class.
    select type (satelliteComponent)
    class is (nodeComponentSatelliteOrbiting)
       ! Proceed only for satellites.
       if (thisNode%isSatellite()) then
          ! Get all required quantities.
          hostNode                 => thisNode          %mergesWith               ()
          position                 =  satelliteComponent%position                 ()
          velocity                 =  satelliteComponent%velocity                 ()
          tidalTensorPathIntegrated=  satelliteComponent%tidalTensorPathIntegrated()
          tidalHeatingNormalized   =  satelliteComponent%tidalHeatingNormalized   ()
          positionTensor           =  Vector_Outer_Product            (         position                          ,symmetrize=.true.)
          radius                   =  Vector_Magnitude                (         position                                            )
          ! Set rate of change of position.
          call satelliteComponent%positionRate                 (                                                     &
               &                                                +(kilo*gigaYear/megaParsec)                          &
               &                                                *satelliteComponent%velocity()                       &
               &                                               )
          ! Other rates are only non-zero if the satellite is at non-zero radius.
          if (radius > 0.0d0) then
             ! Calcluate tidal tensor and rate of change of integrated tidal tensor.
             parentDensity         =  Galactic_Structure_Density      (hostNode,position,coordinateSystemCartesian                  )
             parentEnclosedMass    =  Galactic_Structure_Enclosed_Mass(hostNode,radius                                              )
             tidalTensor           = -(gravitationalConstantGalacticus*parentEnclosedMass         /radius**3)*tensorIdentityR2D3Sym &
                  &                  +(gravitationalConstantGalacticus*parentEnclosedMass*3.0d0   /radius**5)*positionTensor        &
                  &                  -(gravitationalConstantGalacticus*parentDensity     *4.0d0*Pi/radius**2)*positionTensor
             ! Compute the orbital period.
             angularFrequency  =  Vector_Magnitude(Vector_Product(position,velocity)) &
                  &              /radius**2                                           &
                  &              *kilo                                                &
                  &              *gigaYear                                            &
                  &              /megaParsec
             radialFrequency   =  abs             (   Dot_Product(position,velocity)) &
                  &              /radius**2                                           &
                  &              *kilo                                                &
                  &              *gigaYear                                            &
                  &              /megaParsec
             ! Find the orbital period. We use the larger of the angular and radial frequencies to avoid numerical problems for purely
             ! radial or purely circular orbits.
             orbitalPeriod     = 2.0d0*Pi/max(angularFrequency,radialFrequency)
             ! Calculate position, velocity, mass loss, integrated tidal tensor, and heating rates. In the direct (i.e. non-dynamical
             ! friction) acceleration we include a factor (1+m_{sat}/m_{host})=m_{sat}/µ (where µ is the reduced mass) to convert
             ! from the two-body problem of satellite and host orbitting their common center of mass to the equivalent one-body
             ! problem (since we're solving for the motion of the satellite relative to the center of the host which is held fixed).
             call satelliteComponent%velocityRate                 (                                                     &
                  &                                                -(kilo*gigaYear/megaParsec)                          &
                  &                                                *gravitationalConstantGalacticus                     &
                  &                                                *parentEnclosedMass                                  &
                  &                                                *position                                            &
                  &                                                /radius**3                                           &
                  &                                                *(                                                   &
                  &                                                  +1.0d0                                             &
                  &                                                  +satelliteComponent%boundMass()                    &
                  &                                                  /parentEnclosedMass                                &
                  &                                                )                                                    &
                  &                                                +satelliteDynamicalFriction_%acceleration (thisNode) &
                  &                                               )
             call satelliteComponent%boundMassRate                (                                                     &
                  &                                                +satelliteTidalStripping_%massLossRate    (thisNode) &
                  &                                               )
             call satelliteComponent%tidalTensorPathIntegratedRate(                                                     &
                  &                                                +tidalTensor                                         &
                  &                                                -tidalTensorPathIntegrated                           &
                  &                                                /orbitalPeriod                                       &
                  &                                               )
             call satelliteComponent%tidalHeatingNormalizedRate   (                                                     &
                  &                                                +satelliteTidalHeatingRate_%heatingRate   (thisNode) &
                  &                                               )
          end if
          ! Get half-mass radii of central and satellite galaxies.
          if (Galactic_Structure_Enclosed_Mass(hostNode,massType=massTypeGalactic,radius=0.0d0) > 0.0d0) then
             halfMassRadiusCentral  =Galactic_Structure_Radius_Enclosing_Mass(hostNode,fractionalMass=0.5d0,massType=massTypeGalactic)
          else
             halfMassRadiusCentral  =0.0d0
          end if
          if (Galactic_Structure_Enclosed_Mass(thisNode,massType=massTypeGalactic,radius=0.0d0) > 0.0d0) then
             halfMassRadiusSatellite=Galactic_Structure_Radius_Enclosing_Mass(thisNode,fractionalMass=0.5d0,massType=massTypeGalactic)
          else
             halfMassRadiusSatellite=0.0d0
          end if
          satelliteMass     =  satelliteComponent%boundMass()
          basicComponent    => thisNode          %basic    ()
          basicMass         =  basicComponent    %mass     ()
          radiusVirial      =  darkMatterHaloScale_%virialRadius(hostNode)
          orbitalRadiusTest =  max(halfMassRadiusSatellite+halfMassRadiusCentral,radiusVirialFraction*radiusVirial)
          ! Test merging criterion.
          if (satelliteOrbitingDestructionMassIsFractional) then
             massDestruction=satelliteOrbitingDestructionMass*basicMass
          else
             massDestruction=satelliteOrbitingDestructionMass
          end if
          if     (                                        &
               &   odeConverged                           &
               &  .and.                                   &
               &   (                                      &
               &     (                                    &
               &       radius        > 0.0d0              &
               &      .and.                               &
               &       radius        <  orbitalRadiusTest &
               &     )                                    &
               &    .or.                                  &
               &       satelliteMass <  massDestruction   &
               &   )                                      &
               & ) then
             ! Merging criterion met - trigger an interrupt.
             interrupt=.true.
             interruptProcedure => Node_Component_Satellite_Orbiting_Trigger_Merger
             return
          end if
       end if
    end select
    return
  end subroutine Node_Component_Satellite_Orbiting_Rate_Compute

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Satellite_Orbiting_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Satellite_Orbiting_Scale_Set(thisNode)
    !% Set scales for properties of {\normalfont \ttfamily thisNode}.
    use Galactic_Structure_Enclosed_Masses
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    use Galacticus_Nodes                  , only : treeNode, nodeComponentSatellite, nodeComponentSatelliteOrbiting
    implicit none
    type            (treeNode              ), pointer  , intent(inout) :: thisNode
    class           (nodeComponentSatellite), pointer                  :: satelliteComponent
    double precision                        , parameter                :: positionScaleFractional                 =1.0d-2                              , &
         &                                                                velocityScaleFractional                 =1.0d-2                              , &
         &                                                                boundMassScaleFractional                =1.0d-2                              , &
         &                                                                tidalTensorPathIntegratedScaleFractional=1.0d-2                              , &
         &                                                                tidalHeatingNormalizedScaleFractional   =1.0d-2
    double precision                                                   :: virialRadius                                   , virialVelocity              , &
         &                                                                virialIntegratedTidalTensor                    , virialTidalHeatingNormalized, &
         &                                                                satelliteMass

    ! Get the satellite component.
    satelliteComponent => thisNode%satellite()
    ! Ensure that it is of the orbiting class.
    select type (satelliteComponent)
    class is (nodeComponentSatelliteOrbiting)
       ! Get characteristic scales.
       satelliteMass               =Galactic_Structure_Enclosed_Mass   (thisNode)
       virialRadius                =darkMatterHaloScale_%virialRadius  (thisNode)
       virialVelocity              =darkMatterHaloScale_%virialVelocity(thisNode)
       virialIntegratedTidalTensor = virialVelocity/virialRadius*megaParsec/kilo/gigaYear
       virialTidalHeatingNormalized=(virialVelocity/virialRadius)**2
       call satelliteComponent%positionScale                 (                                          &
            &                                                  [1.0d0,1.0d0,1.0d0]                      &
            &                                                 *virialRadius                             &
            &                                                 *                 positionScaleFractional &
            &                                                )
       call satelliteComponent%velocityScale                 (                                          &
            &                                                  [1.0d0,1.0d0,1.0d0]                      &
            &                                                 *virialVelocity                           &
            &                                                 *                 velocityScaleFractional &
            &                                                )
       call satelliteComponent%boundMassScale                (                                          &
            &                                                  satelliteMass                            &
            &                                                 *                boundMassScaleFractional &
            &                                                )
       call satelliteComponent%tidalTensorPathIntegratedScale(                                          &
            &                                                  tensorUnitR2D3Sym                        &
            &                                                 *virialIntegratedTidalTensor              &
            &                                                 *tidalTensorPathIntegratedScaleFractional &
            &                                                )
       call satelliteComponent%tidalHeatingNormalizedScale   (                                          &
            &                                                  virialTidalHeatingNormalized             &
            &                                                 *   tidalHeatingNormalizedScaleFractional &
            &                                                )
    end select
    return
  end subroutine Node_Component_Satellite_Orbiting_Scale_Set

  !# <nodeMergerTask>
  !#  <unitName>Node_Component_Satellite_Orbiting_Create</unitName>
  !# </nodeMergerTask>
  subroutine Node_Component_Satellite_Orbiting_Create(thisNode)
    !% Create a satellite orbit component and assign initial position, velocity, orbit, and tidal heating quantities. (The initial
    !% bound mass is automatically set to the original halo mass by virtue of that being the class default).
    use Satellite_Merging_Timescales
    use Galacticus_Nodes            , only : treeNode, nodeComponentSatellite, nodeComponentSatelliteOrbiting, defaultSatelliteComponent
    implicit none
    type            (treeNode              ), pointer, intent(inout) :: thisNode
    type            (treeNode              ), pointer                :: hostNode
    class           (nodeComponentSatellite), pointer                :: satelliteComponent
    logical                                                          :: isNewSatellite
    type            (keplerOrbit           )                         :: thisOrbit

    ! Return immediately if this method is not active.
    if (.not.defaultSatelliteComponent%orbitingIsActive()) return
    ! Get the satellite component.
    satelliteComponent => thisNode%satellite()
    ! Determine if the satellite component exists already.
    isNewSatellite=.false.
    select type (satelliteComponent)
    type is (nodeComponentSatellite)
       isNewSatellite=.true.
    end select
    ! Return if this is not a new satellite.
    if (.not.isNewSatellite) return
    ! Create the satellite component.
    satelliteComponent => thisNode%satellite(autoCreate=.true.)
    select type (satelliteComponent)
    class is (nodeComponentSatelliteOrbiting)
       ! Get an orbit for this satellite.
       if (thisNode%isSatellite()) then
          hostNode => thisNode%parent
       else
          hostNode => thisNode%parent%firstChild
       end if
       thisOrbit=virialOrbit_%orbit(thisNode,hostNode,acceptUnboundOrbits)
       ! Store the orbit.
       call satelliteComponent%virialOrbitSet(thisOrbit)
       ! Set the initial bound mass of this satellite.
       call Node_Component_Satellite_Orbiting_Bound_Mass_Initialize(satelliteComponent,thisNode)
    end select
    return
  end subroutine Node_Component_Satellite_Orbiting_Create

  subroutine Node_Component_Satellite_Orbiting_Trigger_Merger(thisNode)
    !% Trigger a merger of the satellite by setting the time until merging to zero.
    use Galacticus_Nodes, only : treeNode, nodeComponentSatellite
    implicit none
    type (treeNode              ), intent(inout), pointer :: thisNode
    class(nodeComponentSatellite)               , pointer :: satelliteComponent

    satelliteComponent => thisNode%satellite()
    call satelliteComponent%mergeTimeSet(0.0d0)
    return
  end subroutine Node_Component_Satellite_Orbiting_Trigger_Merger

  subroutine Node_Component_Satellite_Orbiting_Virial_Orbit_Set(self,thisOrbit)
    !% Set the orbit of the satellite at the virial radius.
    use Numerical_Constants_Math
    use Vectors
    use Galacticus_Nodes        , only : nodeComponentSatellite, nodeComponentSatelliteOrbiting
    implicit none
    class           (nodeComponentSatellite), intent(inout) :: self
    type            (keplerOrbit           ), intent(in   ) :: thisOrbit
    double precision                        , dimension(3)  :: radialVector                    , velocityRadialVector     , &
         &                                                     velocityTangentialVector1       , velocityTangentialVector2
    double precision                                        :: orbitalRadius                   , orbitalVelocityRadial    , &
         &                                                     orbitalVelocityTangential       , velocityPhi              , &
         &                                                     orbitalPositionPhi              , orbitalPositionTheta
    type            (keplerOrbit           )                :: virialOrbit

    select type (self)
    class is (nodeComponentSatelliteOrbiting)
       ! Ensure the orbit is defined.
       call thisOrbit%assertIsDefined()
       ! Store the orbit.
       call self%virialOrbitSetValue(thisOrbit)
       ! Compute and store orbitial position and velocity.
       virialOrbit              =thisOrbit
       orbitalRadius            =virialOrbit%radius()
       orbitalVelocityRadial    =virialOrbit%velocityRadial()
       orbitalVelocityTangential=virialOrbit%velocityTangential()
       orbitalPositionPhi       =     2.0d0*Pi*self%hostNode%hostTree%randomNumberGenerator%uniformSample()
       orbitalPositionTheta     =acos(2.0d0   *self%hostNode%hostTree%randomNumberGenerator%uniformSample()-1.0d0)
       radialVector             =[                                                   &
            &                     sin(orbitalPositionTheta)*cos(orbitalPositionPhi), &
            &                     sin(orbitalPositionTheta)*sin(orbitalPositionPhi), &
            &                     cos(orbitalPositionTheta)                          &
            &                    ]
       call self%positionSet(orbitalRadius*radialVector)
       velocityRadialVector     =orbitalVelocityRadial*radialVector
       velocityTangentialVector1=Vector_Product(radialVector,[1.0d0,0.0d0,0.0d0]      )
       ! Normalization.
       velocityTangentialVector1=velocityTangentialVector1/sqrt(sum(velocityTangentialVector1**2))
       velocityTangentialVector2=Vector_Product(radialVector,velocityTangentialVector1)
       velocityPhi              =self%hostNode%hostTree%randomNumberGenerator%uniformSample()*2.0d0*Pi
       velocityTangentialVector1=velocityTangentialVector1*orbitalVelocityTangential*cos(velocityPhi)
       velocityTangentialVector2=velocityTangentialVector2*orbitalVelocityTangential*sin(velocityPhi)
       call self%velocitySet(velocityRadialVector+velocityTangentialVector1+velocityTangentialVector2)
       ! Set the merging time to -1 to indicate that we don't know when merging will occur.
       call self%mergeTimeSet                (           -1.0d0)
       call self%tidalTensorPathIntegratedSet(tensorNullR2D3Sym)
       call self%tidalHeatingNormalizedSet   (            0.0d0)       
    end select
    return
  end subroutine Node_Component_Satellite_Orbiting_Virial_Orbit_Set

  subroutine Node_Component_Satellite_Orbiting_Bound_Mass_Initialize(satelliteComponent,thisNode)
    !% Set the initial bound mass of the satellite.
    use Galacticus_Error
    use Galacticus_Nodes                    , only : nodeComponentSatellite          , treeNode, nodeComponentSatelliteOrbiting
    use Galactic_Structure_Enclosed_Masses  , only : Galactic_Structure_Enclosed_Mass
    use Dark_Matter_Profile_Mass_Definitions
    implicit none
    class           (nodeComponentSatellite)              , intent(inout) :: satelliteComponent
    type            (treeNode              ), pointer     , intent(in   ) :: thisNode
    double precision                                                      :: virialRadius      , maximumRadius, &
         &                                                                   satelliteMass

    select type (satelliteComponent)
    class is (nodeComponentSatelliteOrbiting)
       if      (satelliteBoundMassInitializeType == 'basicMass'      ) then
          ! Do nothing. The bound mass of this satellite is set to the node mass by default.
          satelliteMass=satelliteComponent%boundMass()
       else if (satelliteBoundMassInitializeType == 'maximumRadius'  ) then
          ! Set the initial bound mass of this satellite by integrating the density profile up to a maximum radius.
          if (satelliteMaximumRadiusOverVirialRadius > 0.0d0) then
             virialRadius  = darkMatterHaloScale_%virialRadius  (thisNode                         )
             maximumRadius = satelliteMaximumRadiusOverVirialRadius*virialRadius
             satelliteMass = Galactic_Structure_Enclosed_Mass   (thisNode,maximumRadius           )
          else
             call Galacticus_Error_Report('specify a positive maximum radius for the satellite'//{introspection:location})
          end if
       else if (satelliteBoundMassInitializeType == 'densityContrast') then
          ! Set the initial bound mass of this satellite by assuming a specified density contrast.
          if (satelliteDensityContrast > 0.0d0) then
             satelliteMass = Dark_Matter_Profile_Mass_Definition(thisNode,satelliteDensityContrast)
          else
             call Galacticus_Error_Report('specify a positive density contrast for the satellite'//{introspection:location})
          end if
       else
          call Galacticus_Error_Report('tpye of method to initialize the bound mass of satellites can not be recognized. Available options are "basicMass", "maximumRadius", "densityContrast"'//{introspection:location})
       end if
       call satelliteComponent%boundMassSet(satelliteMass)
    end select
    return
  end subroutine Node_Component_Satellite_Orbiting_Bound_Mass_Initialize

end module Node_Component_Satellite_Orbiting
