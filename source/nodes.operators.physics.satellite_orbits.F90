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

  !+    Contributions to this file made by: Andrew Benson, Xiaolong Du.

  !!{RST
  Implements a node operator class that propagates satellite halos along their orbits.
  !!}

  use :: Satellites_Bound_Mass_Initialize, only : satelliteMassBoundInitializorClass
  use :: Virial_Orbits                   , only : virialOrbitClass

  !![
  <nodeOperator name="nodeOperatorSatelliteOrbit" docformat="rst">
   <description>
   A node operator class that integrates the orbital motion of satellite halos through the potential of their host halo, updating position and velocity at each ODE timestep via the equations of motion in the host potential. ``trackPreInfallOrbit`` enables approximate orbit integration before formal infall, which is useful for modeling pre-infall tidal effects and environmental processes.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorSatelliteOrbit
     !!{RST
     A node operator class that propagates satellite halos along their orbits.
     !!}
     private
     class  (virialOrbitClass                   ), pointer :: virialOrbit_                   => null()
     class  (satelliteMassBoundInitializorClass ), pointer :: satelliteMassBoundInitializor_ => null()
     logical                                               :: trackPreInfallOrbit                     , acceptUnboundOrbits, &
          &                                                   initializeOnly
     integer                                               :: rateGrowthMassBoundID
   contains
     final     ::                          satelliteOrbitDestructor
     procedure :: autoHook              => satelliteOrbitAutoHook
     procedure :: nodeTreeInitialize    => satelliteOrbitNodeTreeInitialize
     procedure :: nodeInitialize        => satelliteOrbitNodeInitialize
     procedure :: nodePromote           => satelliteOrbitNodePromote
     procedure :: nodesMerge            => satelliteOrbitNodeMerge
     procedure :: differentialEvolution => satelliteOrbitDifferentialEvolution
  end type nodeOperatorSatelliteOrbit
  
  interface nodeOperatorSatelliteOrbit
     !!{RST
     Constructors for the ``nodeOperatorSatelliteOrbit`` node operator class.
     !!}
     module procedure satelliteOrbitConstructorParameters
     module procedure satelliteOrbitConstructorInternal
  end interface nodeOperatorSatelliteOrbit

  ! Allowed range of satellite/host mass ratios.
  double precision                            , parameter :: massRatioMinimum=0.0d0, massRatioMaximum=1.0d3

  ! Sub-module scope objects used in ODE solving.
  class           (nodeOperatorSatelliteOrbit), pointer    :: self_
  type            (treeNode                  ), pointer    :: nodeProgenitor
  !$omp threadprivate(nodeProgenitor,self_)
  
contains

  function satelliteOrbitConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the ``nodeOperatorSatelliteOrbit`` node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nodeOperatorSatelliteOrbit        )                :: self
    type   (inputParameters                   ), intent(inout) :: parameters
    class  (virialOrbitClass                  ), pointer       :: virialOrbit_
    class  (satelliteMassBoundInitializorClass), pointer       :: satelliteMassBoundInitializor_
    logical                                                    :: trackPreInfallOrbit           , acceptUnboundOrbits, &
         &                                                        initializeOnly

    !![
    <inputParameter docformat="rst">
      <name>trackPreInfallOrbit</name>
      <defaultValue>.false.</defaultValue>
      <description>
      If true, (approximately) track the orbits of halos prior to infall.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>acceptUnboundOrbits</name>
      <defaultValue>.false.</defaultValue>
      <description>
      If true, accept unbound virial orbits for satellites, otherwise reject them.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>initializeOnly</name>
      <defaultValue>.false.</defaultValue>
      <description>
      If true, orbits are initialized, but not evolved.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="virialOrbit"                   name="virialOrbit_"                   source="parameters"/>
    <objectBuilder class="satelliteMassBoundInitializor" name="satelliteMassBoundInitializor_" source="parameters"/>
    !!]
    self=nodeOperatorSatelliteOrbit(trackPreInfallOrbit,acceptUnboundOrbits,initializeOnly,virialOrbit_,satelliteMassBoundInitializor_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="virialOrbit_"                  />
    <objectDestructor name="satelliteMassBoundInitializor_"/>
    !!]
    return
  end function satelliteOrbitConstructorParameters
  
  function satelliteOrbitConstructorInternal(trackPreInfallOrbit,acceptUnboundOrbits,initializeOnly,virialOrbit_,satelliteMassBoundInitializor_) result(self)
    !!{RST
    Internal constructor for the ``nodeOperatorSatelliteOrbit`` node operator class.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nodeOperatorSatelliteOrbit        )                        :: self
    logical                                    , intent(in   )         :: trackPreInfallOrbit           , acceptUnboundOrbits, &
         &                                                                initializeOnly
    class  (virialOrbitClass                  ), intent(in   ), target :: virialOrbit_
    class  (satelliteMassBoundInitializorClass), intent(in   ), target :: satelliteMassBoundInitializor_
    !![
    <constructorAssign variables="trackPreInfallOrbit, acceptUnboundOrbits, initializeOnly, *virialOrbit_, *satelliteMassBoundInitializor_"/>
    !!]

    ! Validate.
    !! Check that the virial orbit class supports setting of angular coordinates.
    if (.not.self%virialOrbit_%isAngularlyResolved()) call Error_Report('this `nodeOperator` requires a `virialOrbit` implementation which provides angularly-resolved orbits'//{introspection:location})
    ! Add pre-infall properties if needed.
    if (self%trackPreInfallOrbit) then
       !![
       <addMetaProperty component="satellite" name="rateGrowthMassBound" id="self%rateGrowthMassBoundID" isEvolvable="no" isCreator="yes"/>
       !!]
    end if
    return
  end function satelliteOrbitConstructorInternal

  subroutine satelliteOrbitAutoHook(self)
    !!{RST
    Attach to the subhalo orbit initialization event.
    !!}
    use :: Events_Hooks, only : subhaloOrbitInitializationEvent, openMPThreadBindingAtLevel
    implicit none
    class(nodeOperatorSatelliteOrbit), intent(inout) :: self

    call subhaloOrbitInitializationEvent%attach(self,subhaloOrbitInitialize,openMPThreadBindingAtLevel,label='nodeOperatorSatelliteOrbits')
    return
  end subroutine satelliteOrbitAutoHook

  subroutine satelliteOrbitDestructor(self)
    !!{RST
    Destructor for the ``nodeOperatorSatelliteOrbit`` node operator function class.
    !!}
    use :: Events_Hooks, only : subhaloOrbitInitializationEvent
    implicit none
    type(nodeOperatorSatelliteOrbit), intent(inout) :: self

    !![
    <objectDestructor name="self%virialOrbit_"                  />
    <objectDestructor name="self%satelliteMassBoundInitializor_"/>
    !!]
    if (subhaloOrbitInitializationEvent%isAttached(self,subhaloOrbitInitialize)) call subhaloOrbitInitializationEvent%detach(self,subhaloOrbitInitialize)
    return
  end subroutine satelliteOrbitDestructor

  subroutine subhaloOrbitInitialize(self,node,orbitIsDefined)
    !!{RST
    Initialize a new subhalo orbit.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : treeNode    , nodeComponentSatellite
    use :: Kepler_Orbits   , only : keplerOrbit
    implicit none
    class  (*                     ), intent(inout)           :: self
    type   (treeNode              ), intent(inout)           :: node
    logical                        , intent(in   ), optional :: orbitIsDefined
    class  (nodeComponentSatellite), pointer                 :: satellite
    type   (treeNode              ), pointer                 :: nodeHost
    type   (keplerOrbit           )                          :: orbit
    logical                                                  :: isNewSatellite, orbitIsDefined_
    
    select type (self)
    class is (nodeOperatorSatelliteOrbit)
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
          call satellite%boundMassSet(self%satelliteMassBoundInitializor_%massBound(node))
       end if
       ! Create an orbit for the satellite if needed.
       if (present(orbitIsDefined)) then
          orbitIsDefined_=orbitIsDefined
       else
          orbit          =satellite%virialOrbit()
          orbitIsDefined_=orbit    %isDefined  ()
       end if
       if (.not.orbitIsDefined_) then
          ! Get an orbit for this satellite.
          if (node%isSatellite()) then
             nodeHost => node%parent
          else
             nodeHost => node%parent%firstChild
          end if
          orbit=self%virialOrbit_%orbit(node,nodeHost,self%acceptUnboundOrbits)
          ! Store the orbit.
          call satellite%virialOrbitSet(orbit)
       end if
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine subhaloOrbitInitialize

  subroutine satelliteOrbitNodeTreeInitialize(self,node)
    !!{RST
    Initialize orbits for any initial subhalos.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite
    implicit none
    class(nodeOperatorSatelliteOrbit), intent(inout), target  :: self
    type (treeNode                  ), intent(inout), target  :: node
    class(nodeComponentSatellite    )               , pointer :: satellite

    ! If this node is an initial subhalo - create a satellite component in it.
    if (node%isSatellite()) call subhaloOrbitInitialize(self,node)
    return
  end subroutine satelliteOrbitNodeTreeInitialize
  
  subroutine satelliteOrbitNodeInitialize(self,node)
    !!{RST
    Estimate the position of nodes relative to their hosts prior to infall.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentBasic, nodeComponentSatellite
    use :: Kepler_Orbits                   , only : keplerOrbit
    use :: Numerical_Constants_Astronomical, only : MpcPerKmPerSToGyr  
    use :: Numerical_ODE_Solvers           , only : odeSolver
    implicit none
    class           (nodeOperatorSatelliteOrbit), intent(inout), target  :: self
    type            (treeNode                  ), intent(inout), target  :: node
    type            (treeNode                  )               , pointer :: nodeHost
    class           (nodeComponentBasic        )               , pointer :: basicHost            , basicProgenitor
    class           (nodeComponentSatellite    )               , pointer :: satellite            , satelliteProgenitor
    double precision                            , dimension(3)           :: positionProgenitor   , velocityProgenitor, &
         &                                                                  positionDescendent   , velocityEffective
    double precision                            , dimension(6)           :: phaseSpaceCoordinates
    type            (odeSolver                 )                         :: solver
    type            (keplerOrbit               )                         :: orbit
    double precision                                                     :: time                 , timeProgenitor     , &
         &                                                                  timeDescendent       , massBoundDescendent, &
         &                                                                  massBound            , rateGrowthMassBound
    !$GLC attributes unused :: self

    ! Compute the pre-infall orbit only for leaf nodes not on the main branch.
    if (.not.self%trackPreInfallOrbit .or. node%isOnMainBranch() .or. associated(node%firstChild)) return
    ! Extract the position, velocity, and time at infall.
    nodeProgenitor => node
    do while (nodeProgenitor%isPrimaryProgenitor())
       nodeProgenitor => nodeProgenitor%parent
    end do
    nodeHost           => nodeProgenitor%parent
    satellite          => nodeProgenitor%satellite  (autoCreate=.true.)
    basicHost          => nodeHost      %basic      (                 )
    orbit              =  satellite     %virialOrbit(                 )
    positionProgenitor =  satellite     %position   (                 )
    velocityProgenitor =  satellite     %velocity   (                 )
    massBound          =  satellite     %boundMass  (                 )
    time               =  basicHost     %time       (                 )
    ! Integrate the orbit backward in time to each progenitor halo.
    solver         =  odeSolver(6_c_size_t,orbitalODEs,toleranceRelative=1.0d-3,toleranceAbsolute=1.0d-6)
    self_          => self
    do while (associated(nodeProgenitor))
       basicProgenitor     => nodeProgenitor %basic             (                 )
       satelliteProgenitor => nodeProgenitor %satellite         (autoCreate=.true.)
       timeProgenitor      =  basicProgenitor%time              (                 )
       positionDescendent  =                  positionProgenitor
       timeDescendent      =                  time
       massBoundDescendent =                  massBound
       ! Solve the ODEs here.
       phaseSpaceCoordinates(1:3)=positionProgenitor
       phaseSpaceCoordinates(4:6)=velocityProgenitor
       call solver%solve(time,timeProgenitor,phaseSpaceCoordinates)
       positionProgenitor=phaseSpaceCoordinates(1:3)
       velocityProgenitor=phaseSpaceCoordinates(4:6)
       ! Compute the constant velocity required to reproduce the position evolution) for this progenitor.
       velocityEffective=+(+positionDescendent-positionProgenitor) &
            &            /(+    timeDescendent-    timeProgenitor) &
            &            *MpcPerKmPerSToGyr
       ! Compute the bound mass growth rate.
       massBound=satelliteProgenitor%boundMass()
       if (nodeProgenitor%isPrimaryProgenitor()) then
          rateGrowthMassBound=+(+massBoundDescendent-massBound     ) &
               &              /(+timeDescendent     -timeProgenitor)
       else
          rateGrowthMassBound=+0.0d0
       end if
       ! Set all properties.
       call satelliteProgenitor%           virialOrbitSet(                           satellite%virialOrbit        ())
       call satelliteProgenitor%              positionSet(                                     positionProgenitor   )
       call satelliteProgenitor%              velocitySet(                                     velocityEffective    )
       call satelliteProgenitor%floatRank0MetaPropertySet(self%rateGrowthMassBoundID,          rateGrowthMassBound  )
       ! Move to the next progenitor.
       time           =  timeProgenitor
       nodeProgenitor => nodeProgenitor%firstChild
    end do
    return
  end subroutine satelliteOrbitNodeInitialize
  
  integer function orbitalODEs(time,phaseSpaceCoordinates,phaseSpaceCoordinatesRateOfChange)
    !!{RST
    ODEs describing a halo orbit.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentBasic
    use :: Galactic_Structure_Options      , only : componentTypeDarkMatterOnly, massTypeDark
    use :: Interface_GSL                   , only : GSL_Success
    use :: Numerical_Constants_Astronomical, only : gigaYear                   , megaParsec, gravitationalConstant_internal, MpcPerKmPerSToGyr
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Mass_Distributions              , only : massDistributionClass
    use :: Vectors                         , only : Vector_Magnitude
    implicit none
    double precision                                     , intent(in   ) :: time
    double precision                       , dimension(:), intent(in   ) :: phaseSpaceCoordinates
    double precision                       , dimension(:), intent(  out) :: phaseSpaceCoordinatesRateOfChange
    double precision                       , dimension(3)                :: position                         , velocity            , &
         &                                                                  acceleration
    type            (treeNode             ), pointer                     :: nodeHost                         , nodeDescendent
    class           (nodeComponentBasic   ), pointer                     :: basicProgenitor                  , basicDescendent
    class           (massDistributionClass), pointer                     :: massDistributionDescendent       , massDistributionHost
    double precision                                                     :: massSatellite                    , massHost            , &
         &                                                                  factorInterpolate                , radius              , &
         &                                                                  massRatio

    ! Extract orbital position and velocity.
    orbitalODEs =GSL_Success
    position    =phaseSpaceCoordinates(1:3)
    velocity    =phaseSpaceCoordinates(4:6)
    ! Find the satellite at this time and estimate its mass.
    if (nodeProgenitor%isPrimaryProgenitor()) then
       nodeDescendent    =>  nodeProgenitor%parent
       basicProgenitor   =>  nodeProgenitor%basic ()
       basicDescendent   =>  nodeDescendent%basic ()
       factorInterpolate =  +(+                time  -basicProgenitor%time()) &
            &               /(+basicDescendent%time()-basicProgenitor%time())
       massSatellite     =  +basicDescendent%mass ()*       factorInterpolate  &
            &               +basicProgenitor%mass ()*(1.0d0-factorInterpolate)
    else
       basicProgenitor   =>  nodeProgenitor %basic()
       massSatellite     =   basicProgenitor%mass ()
    end if
    ! Find the host halo at this time and estimate its mass.
    nodeHost => nodeProgenitor
    do while (nodeHost%isPrimaryProgenitor())
       nodeHost => nodeHost%parent
    end do
    nodeHost        => nodeHost%parent%firstChild
    basicProgenitor => nodeHost%basic()
    do while (associated(nodeHost) .and. basicProgenitor%time() > time)
       nodeHost => nodeHost%firstChild
       if (associated(nodeHost)) basicProgenitor => nodeHost%basic()
    end do
    if (associated(nodeHost)) then
       nodeDescendent             => nodeHost      %parent
       basicDescendent            => nodeDescendent%basic           (                                        )
       massDistributionHost       => nodeHost      %massDistribution(componentTypeDarkMatterOnly,massTypeDark)
       massDistributionDescendent => nodeDescendent%massDistribution(componentTypeDarkMatterOnly,massTypeDark)
       radius                     =  Vector_Magnitude(position)
       factorInterpolate          =  +(+                time  -basicProgenitor%time()) &
            &                        /(+basicDescendent%time()-basicProgenitor%time())
       massHost                   =  +massDistributionDescendent%massEnclosedBySphere(radius)*       factorInterpolate  &
            &                        +massDistributionHost      %massEnclosedBySphere(radius)*(1.0d0-factorInterpolate)
       massRatio                  =  min(                       &
            &                                +massRatioMaximum, &
            &                            max(                   &
            &                                +massRatioMinimum, &
            &                                +massSatellite     &
            &                                /massHost          &
            &                               )                   &
            &                           )
       acceleration               =  -gravitationalConstant_internal    &
            &                        *massHost                          &
            &                        *position                          &
            &                        /radius                        **3 &
            &                        /MpcPerKmPerSToGyr
       !![
       <objectDestructor name="massDistributionHost"      />
       <objectDestructor name="massDistributionDescendent"/>
       !!]
    else
       ! No host exists at this time, assume zero acceleration.
       acceleration=0.0d0
       massRatio   =1.0d0
    end if
    ! Find the mass ratio.
    ! Set evolution rates.
    phaseSpaceCoordinatesRateOfChange(1:3)=+velocity                &
         &                                 /MpcPerKmPerSToGyr
    phaseSpaceCoordinatesRateOfChange(4:6)=+acceleration            &
         &                                 *(                       &
         &                                   +1.0d0                 &
         &                                   +massRatio             &
         &                                  )
    return
  end function orbitalODEs

  subroutine satelliteOrbitNodePromote(self,node)
    !!{RST
    Act on promotion of a node. Set the position and velocity to those of the parent.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite
    implicit none
    class(nodeOperatorSatelliteOrbit), intent(inout) :: self
    type (treeNode                  ), intent(inout) :: node
    class(nodeComponentSatellite    ), pointer       :: satellite, satelliteParent

    if (.not.self%trackPreInfallOrbit  ) return
    if (     node%isOnMainBranch     ()) return
    if (     self%initializeOnly       ) return
    satellite       => node       %satellite()
    satelliteParent => node%parent%satellite()   
    call        satellite%              positionSet(                           satelliteParent%position                 (                          ))
    call        satellite%              velocitySet(                           satelliteParent%velocity                 (                          ))
    if (self%trackPreInfallOrbit) &
         & call satellite%floatRank0MetaPropertySet(self%rateGrowthMassBoundID,satelliteParent%floatRank0MetaPropertyGet(self%rateGrowthMassBoundID))
    return
  end subroutine satelliteOrbitNodePromote

  subroutine satelliteOrbitNodeMerge(self,node)
    !!{RST
    Act on a merger between nodes.
    !!}
    implicit none
    class(nodeOperatorSatelliteOrbit), intent(inout) :: self
    type (treeNode                  ), intent(inout) :: node

    call subhaloOrbitInitialize(self,node)
    return
  end subroutine satelliteOrbitNodeMerge

  subroutine satelliteOrbitDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{RST
    Perform evolution of a satellite orbit due to its velocity and the acceleration of its host's potential.
    !!}
    use :: Error                           , only : Error_Report
    use :: Galacticus_Nodes                , only : nodeComponentSatellite
    use :: Numerical_Constants_Astronomical, only : gigaYear              , megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Vectors                         , only : Vector_Magnitude
    use :: Mass_Distributions              , only : massDistributionClass
    use :: Coordinates                     , only : coordinateCartesian   , assignment(=)
    implicit none
    class           (nodeOperatorSatelliteOrbit), intent(inout), target  :: self
    type            (treeNode                  ), intent(inout), target  :: node
    logical                                     , intent(inout)          :: interrupt
    procedure       (interruptTask             ), intent(inout), pointer :: functionInterrupt
    integer                                     , intent(in   )          :: propertyType
    type            (treeNode                  )               , pointer :: nodeHost
    class           (nodeComponentSatellite    )               , pointer :: satellite
    class           (massDistributionClass     )               , pointer :: massDistribution_, massDistributionHost_
    double precision                            , dimension(3)           :: position         , velocity             , &
         &                                                                  acceleration
    double precision                                                     :: massEnclosedHost , massEnclosedSatellite, &
         &                                                                  radius           , massRatio
    type            (coordinateCartesian       )                         :: coordinates
    !$GLC attributes unused :: interrupt, functionInterrupt, propertyType
    
    ! Ignore the main branch, and non-satellites unless we are tracking pre-infall orbits.
    if     (                                   &
         &          node%isOnMainBranch     () &
         &  .or.                               &
         &          self%initializeOnly        &
         &  .or.                               &
         &   (                                 &
         &     .not.self%trackPreInfallOrbit   &
         &    .and.                            &
         &     .not.node%isSatellite        () &
         &   )                                 &
         & ) return
    ! Extract satellite properties and set the rate of change of position to equal the velocity.
    satellite => node     %satellite (        )
    nodeHost  => node     %mergesWith(        )
    position  =  satellite%position  (        )
    velocity  =  satellite%velocity  (        )
    radius    =  Vector_Magnitude    (position)
    call satellite%positionRate(            &
         &                      +kilo       &
         &                      *gigaYear   &
         &                      /megaParsec &
         &                      *velocity   &
         &                     )
    if (self%trackPreInfallOrbit)                                                                        &
         & call satellite%boundMassRate(satellite%floatRank0MetaPropertyGet(self%rateGrowthMassBoundID))
    ! If the node is not a satellite, we assume no acceleration (if tracking pre-infall orbits we are using a kick-drift approach,
    ! so the velocity remains constant between kicks).
    if (.not.node%isSatellite()) return
    if (radius <= 0.0d0) return ! If radius is non-positive, assume no acceleration.
    coordinates           =  position
    massDistribution_     => node    %massDistribution()
    massDistributionHost_ => nodeHost%massDistribution()
    massEnclosedSatellite =  max(                                                             &
         &                                                 0.0d0                            , &
         &                       min(                                                         &
         &                           massDistribution_    %massEnclosedBySphere(radius     ), &
         &                           satellite            %boundMass           (           )  &
         &                          )                                                         &
         &                      )
    massEnclosedHost      =          massDistributionHost_%massEnclosedBySphere(radius     )
    acceleration          =          massDistributionHost_%acceleration        (coordinates)
    !![
    <objectDestructor name="massDistribution_"    />
    <objectDestructor name="massDistributionHost_"/>
    !!]
    ! Include a factor (1+m_{sat}/m_{host})=m_{sat}/µ (where µ is the reduced mass) to convert from the two-body problem of
    ! satellite and host orbiting their common center of mass to the equivalent one-body problem (since we're solving for the
    ! motion of the satellite relative to the center of the host which is held fixed).
    if (massEnclosedHost > 0.0d0) then
       massRatio=min(                            &
            &            +massRatioMaximum     , &
            &        max(                        &
            &            +massRatioMinimum     , &
            &            +massEnclosedSatellite  &
            &            /massEnclosedHost       &
            &           )                        &
            &       )
       call satellite%velocityRate(              &
            &                      +acceleration &
            &                      *(            &
            &                        +1.0d0      &
            &                        +massRatio  &
            &                      )             &
            &                     )
    else
       ! Enclosed mass is zero - this is acceptable only if the acceleration is zero.
       if (any(acceleration /= 0.0d0)) call Error_Report('zero host mass but non-zero acceleration'//{introspection:location})          
    end if
    return
  end subroutine satelliteOrbitDifferentialEvolution
  
