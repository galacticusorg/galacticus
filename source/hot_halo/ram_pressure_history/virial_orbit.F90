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

  !+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implements a hot halo ram pressure history class which integrates the satellite's virial orbit through its host.
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  use :: Galacticus_Nodes       , only : treeNode
  use :: Kind_Numbers           , only : kind_int8
  use :: Mass_Distributions     , only : massDistributionClass

  !![
  <hotHaloRamPressureHistory name="hotHaloRamPressureHistoryVirialOrbit" docformat="rst">
   <description>
   A hot halo ram pressure history class which reconstructs the history of the ram pressure experienced by a satellite by
   integrating its orbit through its host, starting from its virial orbit at the time of infall. It therefore requires only the
   virial orbit to be stored by the satellite component, and so works with the ``standard`` satellite component, which does not
   track the satellite's position and velocity.

   The satellite is placed on its virial orbit at the time at which it was last an isolated halo, and its motion is integrated
   up to the present time in the mass distribution of its host as it is at the present time, treated as static. The
   gravitational acceleration is that used by :galacticus-class:`nodeOperatorSatelliteOrbit`, including the factor of
   :math:`(1+M_\mathrm{sat}/M_\mathrm{host})` which converts the two-body problem into the equivalent one-body problem. No
   dynamical friction, or loss of mass by the satellite, is included. Along this orbit the ram pressure is
   :math:`p=\rho_\mathrm{hot, host}(r) v^2`, where :math:`\rho_\mathrm{hot, host}(r)` is the density of the host's hot halo at
   the orbital radius, :math:`r`, and :math:`v` is the orbital speed, and the time integral, :math:`I`, and running maximum,
   :math:`P`, of :math:`p` are accumulated exactly as by the :galacticus-class:`nodeOperatorRamPressureHistory` node operator,
   including the decay over the orbital period controlled by ``[efficiencyDecay]``. In a static host, and in the absence of
   dynamical friction and mass loss, the two therefore agree. (In both, :math:`P` is the largest ram pressure found at the end of
   any ODE step, so it falls slightly short of the true pericentric peak, by an amount which depends on how closely a step ends at
   pericenter.)

   The orbit is integrated afresh each time the history of a satellite is requested at a new time, so this class is
   substantially more expensive than :galacticus-class:`hotHaloRamPressureHistoryAccumulated`.
   </description>
  </hotHaloRamPressureHistory>
  !!]
  type, extends(hotHaloRamPressureHistoryClass) :: hotHaloRamPressureHistoryVirialOrbit
     !!{RST
     Implementation of a hot halo ram pressure history class which integrates the satellite's virial orbit through its host.
     !!}
     private
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_        => null()
     double precision                                    :: efficiencyDecay                      , toleranceRelative
     integer         (kind_int8               )          :: uniqueIDPrevious
     double precision                                    :: timePrevious                         , forceMaximumPrevious, &
          &                                                 forceTimeIntegratedPrevious
   contains
     !![
     <methods>
       <method method="compute" description="Compute the ram pressure history of the given node."/>
     </methods>
     !!]
     final     ::                        virialOrbitDestructor
     procedure :: forceMaximum        => virialOrbitForceMaximum
     procedure :: forceTimeIntegrated => virialOrbitForceTimeIntegrated
     procedure :: compute             => virialOrbitCompute
  end type hotHaloRamPressureHistoryVirialOrbit

  interface hotHaloRamPressureHistoryVirialOrbit
     !!{RST
     Constructors for the :galacticus-class:`hotHaloRamPressureHistoryVirialOrbit` hot halo ram pressure history class.
     !!}
     module procedure virialOrbitConstructorParameters
     module procedure virialOrbitConstructorInternal
  end interface hotHaloRamPressureHistoryVirialOrbit

  ! Objects used in the ODE functions.
  class           (hotHaloRamPressureHistoryVirialOrbit), pointer :: self_
  class           (massDistributionClass               ), pointer :: massDistributionHost_     , massDistributionHotHalo_, &
       &                                                             massDistributionSatellite_
  double precision                                                :: massBoundSatellite_       , timescaleDynamicalHost_
  !$omp threadprivate(self_,massDistributionHost_,massDistributionHotHalo_,massDistributionSatellite_,massBoundSatellite_,timescaleDynamicalHost_)

  ! Indices of the variables in the ODE system.
  integer, parameter :: indexPosition=1, indexVelocity=3, indexTimeIntegrated=5, indexMaximum=6, countVariables=6

contains

  function virialOrbitConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`hotHaloRamPressureHistoryVirialOrbit` hot halo ram pressure history class which builds
    the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (hotHaloRamPressureHistoryVirialOrbit)                :: self
    type            (inputParameters                     ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass            ), pointer       :: darkMatterHaloScale_
    double precision                                                      :: efficiencyDecay     , toleranceRelative

    !![
    <inputParameter docformat="rst">
      <name>efficiencyDecay</name>
      <defaultValue>1.0d0</defaultValue>
      <description>
      Efficiency of the decay of the time integral and running maximum of the ram pressure, in units of the inverse orbital period.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>toleranceRelative</name>
      <defaultValue>1.0d-6</defaultValue>
      <description>
      The relative tolerance used in integrating the satellite's orbit.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=hotHaloRamPressureHistoryVirialOrbit(efficiencyDecay,toleranceRelative,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function virialOrbitConstructorParameters

  function virialOrbitConstructorInternal(efficiencyDecay,toleranceRelative,darkMatterHaloScale_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`hotHaloRamPressureHistoryVirialOrbit` hot halo ram pressure history class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (hotHaloRamPressureHistoryVirialOrbit)                        :: self
    class           (darkMatterHaloScaleClass            ), intent(in   ), target :: darkMatterHaloScale_
    double precision                                      , intent(in   )         :: efficiencyDecay     , toleranceRelative
    !![
    <constructorAssign variables="efficiencyDecay, toleranceRelative, *darkMatterHaloScale_"/>
    !!]

    if (efficiencyDecay   <  0.0d0) call Error_Report('[efficiencyDecay] must be non-negative'//{introspection:location})
    if (toleranceRelative <= 0.0d0) call Error_Report('[toleranceRelative] must be positive'  //{introspection:location})
    self%uniqueIDPrevious           =-1_kind_int8
    self%timePrevious               =-huge(0.0d0)
    self%forceMaximumPrevious       =+0.0d0
    self%forceTimeIntegratedPrevious=+0.0d0
    !![
    <componentPropertyAssert class="satellite" properties="virialOrbit" require="gettable"/>
    !!]
    return
  end function virialOrbitConstructorInternal

  subroutine virialOrbitDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`hotHaloRamPressureHistoryVirialOrbit` hot halo ram pressure history class.
    !!}
    implicit none
    type(hotHaloRamPressureHistoryVirialOrbit), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine virialOrbitDestructor

  double precision function virialOrbitForceMaximum(self,node)
    !!{RST
    Return the running maximum of the ram pressure force along the satellite's virial orbit.
    !!}
    implicit none
    class(hotHaloRamPressureHistoryVirialOrbit), intent(inout) :: self
    type (treeNode                            ), intent(inout) :: node

    call self%compute(node)
    virialOrbitForceMaximum=self%forceMaximumPrevious
    return
  end function virialOrbitForceMaximum

  double precision function virialOrbitForceTimeIntegrated(self,node)
    !!{RST
    Return the time integral of the ram pressure force along the satellite's virial orbit.
    !!}
    implicit none
    class(hotHaloRamPressureHistoryVirialOrbit), intent(inout) :: self
    type (treeNode                            ), intent(inout) :: node

    call self%compute(node)
    virialOrbitForceTimeIntegrated=self%forceTimeIntegratedPrevious
    return
  end function virialOrbitForceTimeIntegrated

  subroutine virialOrbitCompute(self,node)
    !!{RST
    Compute the history of the ram pressure force along the satellite's virial orbit, by integrating the orbit from the time of
    infall to the present time.
    !!}
    use            :: Galactic_Structure_Options      , only : componentTypeHotHalo, massTypeGaseous
    use            :: Galacticus_Nodes                , only : nodeComponentBasic  , nodeComponentSatellite
    use, intrinsic :: ISO_C_Binding                   , only : c_size_t
    use            :: Kepler_Orbits                   , only : keplerOrbit
    use            :: Numerical_Constants_Astronomical, only : MpcPerKmPerSToGyr
    use            :: Numerical_Constants_Math        , only : Pi
    use            :: Numerical_ODE_Solvers           , only : odeSolver
    implicit none
    class           (hotHaloRamPressureHistoryVirialOrbit), intent(inout), target     :: self
    type            (treeNode                            ), intent(inout)             :: node
    class           (nodeComponentBasic                  ), pointer                   :: basic                     , basicHost
    class           (nodeComponentSatellite              ), pointer                   :: satellite
    type            (treeNode                            ), pointer                   :: nodeHost
    double precision                                      , dimension(countVariables) :: variables                 , scales
    double precision                                      , parameter                 :: fractionStepInitial=1.0d-3
    type            (odeSolver                           )                            :: solver
    type            (keplerOrbit                         )                            :: orbit
    double precision                                                                  :: time                      , timeInfall    , &
         &                                                                               radiusVirial              , velocityVirial, &
         &                                                                               ramPressureVirial

    ! Non-satellites have no ram pressure history.
    if (.not.node%isSatellite()) then
       self%forceMaximumPrevious       =0.0d0
       self%forceTimeIntegratedPrevious=0.0d0
       self%uniqueIDPrevious           =-1_kind_int8
       return
    end if
    ! Return previously computed values if possible.
    basic => node%basic()
    time  =  basic%time()
    if (node%uniqueID() == self%uniqueIDPrevious .and. time == self%timePrevious) return
    self%uniqueIDPrevious           =node%uniqueID()
    self%timePrevious               =time
    self%forceMaximumPrevious       =0.0d0
    self%forceTimeIntegratedPrevious=0.0d0
    ! Nothing has accumulated if no time has passed since infall.
    timeInfall=basic%timeLastIsolated()
    if (time <= timeInfall) return
    ! Find the host and the satellite's virial orbit.
    nodeHost  => node     %parent
    basicHost => nodeHost %basic      ()
    satellite => node     %satellite  ()
    orbit     =  satellite%virialOrbit()
    ! Get the mass distributions needed in the ODE functions.
    self_                      => self
    massDistributionHost_      => nodeHost                      %massDistribution  (                                    )
    massDistributionHotHalo_   => nodeHost                      %massDistribution  (componentTypeHotHalo,massTypeGaseous)
    massDistributionSatellite_ => node                          %massDistribution  (                                    )
    massBoundSatellite_        =  satellite                     %boundMass         (                                    )
    timescaleDynamicalHost_    =  self     %darkMatterHaloScale_%timescaleDynamical(nodeHost                            )
    ! Set scales for the ODE variables. The ram pressure scale is the mean density within the host's virial radius times its virial
    ! velocity squared, and the scale for its time integral is that times the virial crossing time.
    radiusVirial     =+self%darkMatterHaloScale_%radiusVirial  (nodeHost)
    velocityVirial   =+self%darkMatterHaloScale_%velocityVirial(nodeHost)
    ramPressureVirial=+3.0d0               &
         &            *basicHost%mass()    &
         &            /4.0d0               &
         &            /Pi                  &
         &            /radiusVirial  **3   &
         &            *velocityVirial**2
    scales(indexPosition      :indexPosition+1)=radiusVirial
    scales(indexVelocity      :indexVelocity+1)=velocityVirial
    scales(indexTimeIntegrated                )=ramPressureVirial*radiusVirial/velocityVirial*MpcPerKmPerSToGyr
    scales(indexMaximum                       )=ramPressureVirial
    ! Place the satellite on its virial orbit, working in the orbital plane, with the satellite initially on the x-axis.
    variables(indexPosition      :indexPosition+1)=[orbit%radius        (),0.0d0                    ]
    variables(indexVelocity      :indexVelocity+1)=[orbit%velocityRadial(),orbit%velocityTangential()]
    variables(indexTimeIntegrated                )=0.0d0
    variables(indexMaximum                       )=0.0d0
    ! Integrate the orbit up to the present time.
    solver=odeSolver(                                            &
         &           int(countVariables,c_size_t)              , &
         &           virialOrbitODEs                           , &
         &           toleranceAbsolute=self%toleranceRelative  , &
         &           toleranceRelative=self%toleranceRelative  , &
         &           scale            =scales                  , &
         &           hStart           =fractionStepInitial       &
         &                             *timescaleDynamicalHost_, &
         &           postStep         =virialOrbitPostStep       &
         &          )
    call solver%solve(timeInfall,time,variables)
    ! The running maximum must include the final state.
    variables(indexMaximum)=max(variables(indexMaximum),virialOrbitRamPressure(variables))
    self%forceMaximumPrevious       =variables(indexMaximum       )
    self%forceTimeIntegratedPrevious=variables(indexTimeIntegrated)
    !![
    <objectDestructor name="massDistributionHost_"     />
    <objectDestructor name="massDistributionHotHalo_"  />
    <objectDestructor name="massDistributionSatellite_"/>
    !!]
    return
  end subroutine virialOrbitCompute

  integer function virialOrbitODEs(time,variables,variablesRateOfChange)
    !!{RST
    The equations of motion of the satellite, and the rates of change of the time integral and running maximum of the ram
    pressure.
    !!}
    use :: Coordinates                     , only : coordinateCartesian, assignment(=)
    use :: Interface_GSL                   , only : GSL_Success
    use :: Numerical_Constants_Astronomical, only : MpcPerKmPerSToGyr
    implicit none
    double precision                     , intent(in   )               :: time
    double precision                     , intent(in   ), dimension(:) :: variables
    double precision                     , intent(  out), dimension(:) :: variablesRateOfChange
    double precision                     , parameter                   :: massRatioMinimum     =0.0d0, massRatioMaximum=1.0d3
    double precision                     , dimension(3)                :: acceleration
    type            (coordinateCartesian)                              :: coordinates
    double precision                                                   :: radius                     , massEnclosedHost      , &
         &                                                                massEnclosedSatellite      , massRatio             , &
         &                                                                rateDecay
    !$GLC attributes unused :: time

    radius                                              =+sqrt(sum(variables(indexPosition:indexPosition+1)**2))
    variablesRateOfChange(indexPosition:indexPosition+1)=+variables(indexVelocity:indexVelocity+1) &
         &                                               /MpcPerKmPerSToGyr
    variablesRateOfChange(indexVelocity:indexVelocity+1)=+0.0d0
    if (radius > 0.0d0) then
       ! Gravitational acceleration, matching that applied by the satelliteOrbit node operator.
       coordinates          =[variables(indexPosition),variables(indexPosition+1),0.0d0]
       acceleration         =massDistributionHost_%acceleration        (coordinates)
       massEnclosedHost     =massDistributionHost_%massEnclosedBySphere(radius     )
       massEnclosedSatellite=max(                                                             &
            &                                         0.0d0                                 , &
            &                    min(                                                         &
            &                        massDistributionSatellite_%massEnclosedBySphere(radius), &
            &                        massBoundSatellite_                                      &
            &                       )                                                         &
            &                   )
       if (massEnclosedHost > 0.0d0) then
          massRatio=min(                            &
               &            +massRatioMaximum     , &
               &        max(                        &
               &            +massRatioMinimum     , &
               &            +massEnclosedSatellite  &
               &            /massEnclosedHost       &
               &           )                        &
               &       )
          variablesRateOfChange(indexVelocity:indexVelocity+1)=+acceleration(1:2) &
               &                                               *(                 &
               &                                                 +1.0d0           &
               &                                                 +massRatio       &
               &                                                )
       end if
    end if
    ! Ram pressure history, with decay over the orbital period.
    rateDecay                                 =+self_%efficiencyDecay              &
         &                                     /virialOrbitPeriod(variables)
    variablesRateOfChange(indexTimeIntegrated)=+virialOrbitRamPressure(variables)  &
         &                                     -rateDecay                          &
         &                                     *variables(indexTimeIntegrated)
    variablesRateOfChange(indexMaximum       )=-rateDecay                          &
         &                                     *variables(indexMaximum       )
    virialOrbitODEs                           = GSL_Success
    return
  end function virialOrbitODEs

  subroutine virialOrbitPostStep(time,variables,status)
    !!{RST
    Update the running maximum of the ram pressure after each step of the orbit integration.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_double    , c_int
    use            :: Interface_GSL, only : GSL_Continue, GSL_Success
    implicit none
    real            (c_double), intent(in   ), value        :: time
    real            (c_double), intent(inout), dimension(*) :: variables
    integer         (c_int   ), intent(inout)               :: status
    double precision                                        :: ramPressure
    !$GLC attributes unused :: time

    ramPressure=virialOrbitRamPressure(variables(1:countVariables))
    if (ramPressure > variables(indexMaximum)) then
       variables(indexMaximum)=ramPressure
       ! Indicate that the state was changed, so that the ODE solver does not re-use derivatives from the end of the step.
       if (status == GSL_Success) status=GSL_Continue
    end if
    return
  end subroutine virialOrbitPostStep

  double precision function virialOrbitRamPressure(variables)
    !!{RST
    Return the ram pressure at the given point on the orbit.
    !!}
    use :: Coordinates, only : coordinateSpherical, assignment(=)
    implicit none
    double precision                     , intent(in   ), dimension(:) :: variables
    type            (coordinateSpherical)                              :: coordinates

    coordinates           =[sqrt(sum(variables(indexPosition:indexPosition+1)**2)),0.0d0,0.0d0]
    virialOrbitRamPressure=+massDistributionHotHalo_%density(coordinates)    &
         &                 *sum(variables(indexVelocity:indexVelocity+1)**2)
    return
  end function virialOrbitRamPressure

  double precision function virialOrbitPeriod(variables)
    !!{RST
    Return the orbital period at the given point on the orbit, defined as in ``Satellite_Orbit_Period()``.
    !!}
    use :: Numerical_Constants_Astronomical, only : MpcPerKmPerSToGyr
    use :: Numerical_Constants_Math        , only : Pi
    implicit none
    double precision, intent(in   ), dimension(:) :: variables
    double precision, parameter                   :: frequencyFractionalTiny=1.0d-6
    double precision                              :: radius                        , frequencyAngular, &
         &                                           frequencyRadial               , frequencyOrbital

    radius          =+sqrt(sum(variables(indexPosition:indexPosition+1)**2))
    if (radius <= 0.0d0) then
       virialOrbitPeriod=timescaleDynamicalHost_
       return
    end if
    frequencyAngular=+abs(                                                       &
         &                +variables(indexPosition  )*variables(indexVelocity+1) &
         &                -variables(indexPosition+1)*variables(indexVelocity  ) &
         &               )                                                       &
         &           /radius**2                                                  &
         &           /MpcPerKmPerSToGyr
    frequencyRadial =+abs(sum(variables(indexPosition:indexPosition+1)*variables(indexVelocity:indexVelocity+1))) &
         &           /radius**2                                                                                   &
         &           /MpcPerKmPerSToGyr
    frequencyOrbital=max(frequencyAngular,frequencyRadial)
    if (frequencyOrbital > frequencyFractionalTiny/timescaleDynamicalHost_) then
       virialOrbitPeriod=+2.0d0            &
            &            *Pi               &
            &            /frequencyOrbital
    else
       virialOrbitPeriod=+timescaleDynamicalHost_
    end if
    return
  end function virialOrbitPeriod
