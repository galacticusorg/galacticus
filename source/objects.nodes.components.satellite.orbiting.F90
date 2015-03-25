!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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
  use Galacticus_Nodes
  use Kepler_Orbits
  use Tensors
  implicit none
  private
  public :: Node_Component_Satellite_Orbiting_Scale_Set     , Node_Component_Satellite_Orbiting_Create         , &
       &    Node_Component_Satellite_Orbiting_Rate_Compute  , Node_Component_Satellite_Orbiting_Tree_Initialize, &
       &    Node_Component_Satellite_Orbiting_Trigger_Merger, Node_Component_Satellite_Orbiting_Initialize

  !# <component>
  !#  <class>satellite</class>
  !#  <name>orbiting</name>
  !#  <isDefault>no</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>position</name>
  !#     <type>real</type>
  !#     <rank>1</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <output labels="[X,Y,Z]" unitsInSI="megaParsec" comment="Orbital position of the node."/> 
  !#     <classDefault>[0.0d0,0.0d0,0.0d0]</classDefault>
  !#   </property>
  !#   <property>
  !#     <name>velocity</name>
  !#     <type>real</type>
  !#     <rank>1</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <output labels="[X,Y,Z]" unitsInSI="kilo" comment="Orbital velocity of the node."/>
  !#     <classDefault>[0.0d0,0.0d0,0.0d0]</classDefault>
  !#   </property>
  !#   <property>
  !#     <name>mergeTime</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <classDefault>-1.0d0</classDefault>
  !#   </property>
  !#   <property>
  !#     <name>timeOfMerging</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <isVirtual>true</isVirtual>
  !#     <attributes isSettable="false" isGettable="true" isEvolvable="false" />
  !#     <classDefault>-1.0d0</classDefault>
  !#     <getFunction>Node_Component_Satellite_Orbiting_Time_Of_Merging</getFunction>
  !#   </property>
  !#   <property>
  !#     <name>boundMass</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <classDefault>selfBasicComponent%mass()</classDefault>
  !#     <output unitsInSI="massSolar" comment="Bound mass of the node."/>
  !#   </property>
  !#   <property>
  !#     <name>virialOrbit</name>
  !#     <type>keplerOrbit</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </property>
  !#   <property>
  !#     <name>tidalTensorPathIntegrated</name>
  !#     <type>tensorRank2Dimension3Symmetric</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#   </property>
  !#   <property>
  !#     <name>tidalHeatingNormalized</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <output unitsInSI="kilo**2/megaParsec**2" comment="Energy/radius^2 of satellite."/>
  !#   </property>
  !#  </properties>
  !#  <functions>objects.nodes.components.satellite.orbiting.bound_functions.inc</functions>
  !# </component>

  ! Record of whether the module has been initialized.
  logical                     :: moduleInitialized=.false.

  ! Option controlling whether or not unbound virial orbits are acceptable.
  logical         , parameter :: acceptUnboundOrbits=.false.

  ! Option controlling minimum separation of satellite and host halo centers before a merger is triggered.
  double precision            :: satelliteOrbitingDestructionMassFraction

contains

  !# <mergerTreePreTreeConstructionTask>
  !#  <unitName>Node_Component_Satellite_Orbiting_Initialize</unitName>
  !# </mergerTreePreTreeConstructionTask>
  subroutine Node_Component_Satellite_Orbiting_Initialize()
    !% Initializes the orbiting satellite methods module.
    use Input_Parameters
    implicit none

    ! Initialize the module if necessary.
    !$omp critical (Node_Component_Satellite_Orbiting_Initialize)
    if (defaultSatelliteComponent%orbitingIsActive().and..not.moduleInitialized) then
       ! Create the spheroid mass distribution.
       !@ <inputParameter>
       !@   <name>satelliteOrbitingDestructionMassFraction</name>
       !@   <defaultValue>0.01</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The fraction of the satellite's initial mass below which the satellite is considered to be tidally destroyed and merged with the central halo.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('satelliteOrbitingDestructionMassFraction',satelliteOrbitingDestructionMassFraction,defaultValue=0.01d0)
       ! Record that the module is now initialized.
       moduleInitialized=.true.
    end if
    !$omp end critical (Node_Component_Satellite_Orbiting_Initialize)
    return
  end subroutine Node_Component_Satellite_Orbiting_Initialize

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Satellite_Orbiting_Tree_Initialize</unitName>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Satellite_Orbiting_Tree_Initialize(thisNode)
    !% Initialize the orbiting satellite component.
    implicit none
    type (treeNode), pointer, intent(inout) :: thisNode

    if (thisNode%isSatellite()) call Node_Component_Satellite_Orbiting_Create(thisNode)
    return
  end subroutine Node_Component_Satellite_Orbiting_Tree_Initialize
  
  !# <rateComputeTask>
  !#  <unitName>Node_Component_Satellite_Orbiting_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Satellite_Orbiting_Rate_Compute(thisNode,interrupt,interruptProcedure)
    !% Compute rate of change for satellite properties.
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Physical
    use Numerical_Constants_Math
    use Galactic_Structure_Enclosed_Masses
    use Vectors
    use Galactic_Structure_Densities
    use Galactic_Structure_Options
    use Satellite_Dynamical_Friction
    use Satellite_Tidal_Stripping
    use Satellite_Tidal_Heating
    implicit none
    type            (treeNode                      ), pointer     , intent(inout) :: thisNode
    logical                                                       , intent(inout) :: interrupt
    procedure       (Interrupt_Procedure_Template  ), pointer     , intent(inout) :: interruptProcedure
    class           (nodeComponentSatellite        ), pointer                     :: satelliteComponent
    class           (nodeComponentBasic            ), pointer                     :: basicComponent
    type            (treeNode                      ), pointer                     :: hostNode
    class           (darkMatterHaloScaleClass      ), pointer                     :: darkMatterHaloScale_
    double precision                                , dimension(3)                :: position,velocity
    double precision                                , parameter                   :: radiusVirialFraction = 1.0d-2
    double precision                                                              :: radius,halfMassRadiusSatellite
    double precision                                                              :: halfMassRadiusCentral,orbitalRadiusTest
    double precision                                                              :: radiusVirial
    double precision                                                              :: orbitalPeriod
    double precision                                                              :: angularVelocity,parentDensity
    double precision                                                              :: parentEnclosedMass,satelliteMass,basicMass
    double precision                                                              :: tidalHeatingNormalized
    type            (tensorRank2Dimension3Symmetric)                              :: tidalTensor,tidalTensorPathIntegrated,positionTensor

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
          positionTensor           =  Vector_Self_Outer_Product       (         position                          )
          radius                   =  Vector_Magnitude                (         position                          )
          parentDensity            =  Galactic_Structure_Density      (hostNode,position,coordinateSystemCartesian)
          parentEnclosedMass       =  Galactic_Structure_Enclosed_Mass(hostNode,radius                            )
          ! Calcluate tidal tensor and rate of change of integrated tidal tensor.
          tidalTensor              =                                                                            &
               & -(gravitationalConstantGalacticus*parentEnclosedMass         /radius**3)*tensorIdentityR2D3Sym &
               & +(gravitationalConstantGalacticus*parentEnclosedMass*3.0d0   /radius**5)*positionTensor        &
               & -(gravitationalConstantGalacticus*parentDensity     *4.0d0*Pi/radius**2)*positionTensor
          angularVelocity=Vector_Magnitude(Vector_Product(position,velocity))/radius**2*kilo*gigaYear/megaParsec
          orbitalPeriod  =2.0d0*Pi/angularVelocity
          ! Calculate position, velocity, mass loss, integrated tidal tensor, and heating rates.
          call satelliteComponent%positionRate                 (                                                            &
               &                                                +(kilo*gigaYear/megaParsec)                                 &
               &                                                *satelliteComponent%velocity()                              &
               &                                               )
          call satelliteComponent%velocityRate                 (                                                            &
               &                                                -(kilo*gigaYear/megaParsec)                                 &
               &                                                *gravitationalConstantGalacticus                            &
               &                                                *Galactic_Structure_Enclosed_Mass         (hostNode,radius) &
               &                                                *position                                                   &
               &                                                /radius**3                                                  &
               &                                                +Satellite_Dynamical_Friction_Acceleration(thisNode       ) &
               &                                               )
          call satelliteComponent%boundMassRate                (                                                            &
               &                                                +Satellite_Tidal_Stripping_Rate           (thisNode       ) &
               &                                               )
          call satelliteComponent%tidalTensorPathIntegratedRate(                                                            &
               &                                                +tidalTensor                                                &
               &                                                -tidalTensorPathIntegrated                                  &
               &                                                /orbitalPeriod                                              &
               &                                               )
          call satelliteComponent%tidalHeatingNormalizedRate   (                                                            &
               &                                                +Satellite_Tidal_Heating_Rate(thisNode                    ) &
               &                                               )         
          ! Get half-mass radii of central and satellite galaxies.
          halfMassRadiusCentral  =Galactic_Structure_Radius_Enclosing_Mass(hostNode,fractionalMass=0.5d0,massType=massTypeGalactic)
          halfMassRadiusSatellite=Galactic_Structure_Radius_Enclosing_Mass(thisNode,fractionalMass=0.5d0,massType=massTypeGalactic)
          ! Convert from initial to final radius resulting from expansion due to tidal heating.
          halfMassRadiusSatellite =                                                &
               &                    halfMassRadiusSatellite                        &
               &                   /(                                              &
               &                      1.0d0                                        &
               &                     -2.0d0                                        &
               &                     *halfMassRadiusSatellite**3                   &
               &                     *tidalHeatingNormalized                       &
               &                     /(                                            &
               &                        gravitationalConstantGalacticus            &
               &                       *0.5d0                                      &
               &                       *Galactic_Structure_Enclosed_Mass(thisNode) &
               &                      )                                            &
               &                    )
          darkMatterHaloScale_ => darkMatterHaloScale         ()
          satelliteMass        =  satelliteComponent%boundMass()
          basicComponent       => thisNode          %basic    ()
          basicMass            =  basicComponent    %mass     ()
          radiusVirial         =  darkMatterHaloScale_%virialRadius(hostNode)
          orbitalRadiusTest    =  max(halfMassRadiusSatellite+halfMassRadiusCentral,radiusVirialFraction*radiusVirial)
          ! Test merging criterion.          
          if (radius < orbitalRadiusTest .or. satelliteMass < satelliteOrbitingDestructionMassFraction*basicMass) then
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
    use Dark_Matter_Halo_Scales
    use Galactic_Structure_Enclosed_Masses
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    implicit none
    type            (treeNode                      ), pointer  , intent(inout) :: thisNode
    class           (nodeComponentSatellite        ), pointer                  :: satelliteComponent
    class           (darkMatterHaloScaleClass      ), pointer                  :: darkMatterHaloScale_
    double precision                                , parameter                :: positionScaleFractional                 =1.0d-2                              , &
         &                                                                        velocityScaleFractional                 =1.0d-2                              , &
         &                                                                        boundMassScaleFractional                =1.0d-2                              , &
         &                                                                        tidalTensorPathIntegratedScaleFractional=1.0d-2                              , &
         &                                                                        tidalHeatingNormalizedScaleFractional   =1.0d-2
    double precision                                                           :: virialRadius                                   , virialVelocity              , &
         &                                                                        virialIntegratedTidalTensor                    , virialTidalHeatingNormalized, &
         &                                                                        satelliteMass

    ! Get the satellite component.
    satelliteComponent => thisNode%satellite()
    ! Ensure that it is of the orbiting class.
    select type (satelliteComponent)
    class is (nodeComponentSatelliteOrbiting)
       ! Get characteristic scales.
       darkMatterHaloScale_        => darkMatterHaloScale                (        )
       satelliteMass               =  Galactic_Structure_Enclosed_Mass   (thisNode)
       virialRadius                =  darkMatterHaloScale_%virialRadius  (thisNode)
       virialVelocity              =  darkMatterHaloScale_%virialVelocity(thisNode)
       virialIntegratedTidalTensor =   virialVelocity/virialRadius*megaParsec/kilo/gigaYear
       virialTidalHeatingNormalized=  (virialVelocity/virialRadius)**2
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
            &                                                  tensorIdentityR2D3Sym                    &
            &                                                 *virialIntegratedTidalTensor              &
            &                                                 *tidalTensorPathIntegratedScaleFractional &
            &                                                )
       call satelliteComponent%tidalHeatingNormalizedScale   (                                          &
            &                                                  virialtidalHeatingNormalized             &
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
    use Numerical_Constants_Math
    use Virial_Orbits
    use Satellite_Merging_Timescales
    use Pseudo_Random
    use Vectors
    implicit none
    type            (treeNode              ), pointer     , intent(inout) :: thisNode
    type            (treeNode              ), pointer                     :: hostNode
    class           (nodeComponentSatellite), pointer                     :: satelliteComponent
    logical                                                               :: isNewSatellite
    type            (keplerOrbit           )                              :: thisOrbit
    double precision                        , dimension(3)                :: radialVector                    , velocityRadialVector     , &
         &                                                                   velocityTangentialVector1       , velocityTangentialVector2
    type            (fgsl_rng              ), save                        :: pseudoSequenceObject
    logical                                 , save                        :: resetSequence            =.true.
    !$omp threadprivate(pseudoSequenceObject,resetSequence)
    double precision                                                      :: orbitalRadius                   , orbitalVelocityRadial    , &
         &                                                                   orbitalVelocityTangential       , orbitalVelocityPhi       , &
         &                                                                   orbitalVelocityTheta            , velocityPhi

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
       hostNode => thisNode%parent
       thisOrbit=Virial_Orbital_Parameters(thisNode,hostNode,acceptUnboundOrbits)
       ! Store the orbit.
       call satelliteComponent%virialOrbitSet(thisOrbit)
       orbitalRadius            =thisOrbit%radius()
       orbitalVelocityRadial    =thisOrbit%velocityRadial()
       orbitalVelocityTangential=thisOrbit%velocityTangential()
       orbitalVelocityPhi       =Pseudo_Random_Get(pseudoSequenceObject,resetSequence)*2.0d0*Pi
       orbitalVelocityTheta     =Pseudo_Random_Get(pseudoSequenceObject,resetSequence)*      Pi
       radialVector             =[                                                   &
            &                     sin(orbitalVelocityTheta)*cos(orbitalVelocityPhi), &
            &                     sin(orbitalVelocityTheta)*sin(orbitalVelocityPhi), &
            &                     cos(orbitalVelocityTheta)                          &
            &                    ]
       call satelliteComponent%positionSet(orbitalRadius*radialVector)
       velocityRadialVector     =-orbitalVelocityRadial*radialVector
       velocityTangentialVector1=Vector_Product(radialVector,[1.0d0,0.0d0,0.0d0]      )
       velocityTangentialVector2=Vector_Product(radialVector,velocityTangentialVector1)
       velocityPhi              =Pseudo_Random_Get(pseudoSequenceObject,resetSequence)*2.0d0*Pi
       velocityTangentialVector1=velocityTangentialVector1*orbitalVelocityTangential*cos(velocityPhi)
       velocityTangentialVector2=velocityTangentialVector2*orbitalVelocityTangential*sin(velocityPhi)
       call satelliteComponent%velocitySet(velocityRadialVector+velocityTangentialVector1+velocityTangentialVector2)
       ! Set the merging time to -1 to indicate that we don't know when merging will occur.
       call satelliteComponent%mergeTimeSet                (           -1.0d0)
       call satelliteComponent%tidalTensorPathIntegratedSet(tensorNullR2D3Sym)
       call satelliteComponent%tidalHeatingNormalizedSet   (            0.0d0)       
    end select
    return
  end subroutine Node_Component_Satellite_Orbiting_Create

  subroutine Node_Component_Satellite_Orbiting_Trigger_Merger(thisNode)
    !% Trigger a merger of the satellite by setting the time until merging to zero.
    implicit none
    type (treeNode              ), intent(inout), pointer :: thisNode
    class(nodeComponentSatellite)               , pointer :: satelliteComponent

    satelliteComponent => thisNode%satellite()
    call satelliteComponent%mergeTimeSet(0.0d0)
    return
  end subroutine Node_Component_Satellite_Orbiting_Trigger_Merger

end module Node_Component_Satellite_Orbiting
