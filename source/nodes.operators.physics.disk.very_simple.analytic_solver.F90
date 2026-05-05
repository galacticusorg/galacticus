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

  !!{
  Implements a node operator class which provides an analytic solution for the evolution of satellite
  galaxies represented by the \mono{verySimple} disk component.
  !!}

  use :: Abundances_Structure         , only : abundances
  use :: Cosmology_Functions          , only : cosmologyFunctionsClass
  use :: Star_Formation_Rates_Disks   , only : starFormationRateDisksClass
  use :: Stellar_Feedback_Outflows    , only : stellarFeedbackOutflowsClass
  use :: Stellar_Population_Properties, only : stellarPopulationPropertiesClass

  !![
  <nodeOperator name="nodeOperatorDiskVerySimpleAnalyticSolver">
   <description>
     A node operator class that analytically integrates the evolution of satellite disks for the
     \refClass{nodeComponentDiskVerySimple} disk component. At the start of each ODE step the
     timescales for fuel depletion, star formation and outflow are computed; the disk gas and
     stellar masses (and optionally their abundances), together with the host hot halo's outflowed
     mass and abundances, are then provided as closed-form analytic solutions during the
     integration. Satellites that are too small to ever be of interest, and that will not merge
     before the present day, are pruned from the tree at this point.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorDiskVerySimpleAnalyticSolver
     !!{
     A node operator class that provides an analytic ODE solution for the \mono{verySimple} disk component.
     !!}
     private
     class           (cosmologyFunctionsClass         ), pointer :: cosmologyFunctions_          => null()
     class           (stellarPopulationPropertiesClass), pointer :: stellarPopulationProperties_ => null()
     class           (stellarFeedbackOutflowsClass    ), pointer :: stellarFeedbackOutflows_     => null()
     class           (starFormationRateDisksClass     ), pointer :: starFormationRateDisks_      => null()
     double precision                                            :: pruneMassGas                          , pruneMassStars         , &
          &                                                         timePresentDay
     logical                                                     :: trackAbundances
     ! Per-step cached state, set by differentialEvolutionPre and read by differentialEvolutionSolveAnalytics.
     logical                                                     :: active
     double precision                                            :: timeStart                             , massGasInitial         , &
          &                                                         massStellarInitial                    , timescaleFuel          , &
          &                                                         timescaleStellar                      , timescaleOutflow
     type            (abundances                      )          :: abundancesGasInitial                  , abundancesStellarInitial, &
          &                                                         yieldMassEffective
   contains
     final     ::                                        diskVerySimpleAnalyticSolverDestructor
     procedure :: differentialEvolutionPre            => diskVerySimpleAnalyticSolverPre
     procedure :: differentialEvolutionAnalytics      => diskVerySimpleAnalyticSolverAnalytics
     procedure :: differentialEvolutionSolveAnalytics => diskVerySimpleAnalyticSolverSolveAnalytics
  end type nodeOperatorDiskVerySimpleAnalyticSolver

  interface nodeOperatorDiskVerySimpleAnalyticSolver
     !!{
     Constructors for the \refClass{nodeOperatorDiskVerySimpleAnalyticSolver} node operator class.
     !!}
     module procedure diskVerySimpleAnalyticSolverConstructorParameters
     module procedure diskVerySimpleAnalyticSolverConstructorInternal
  end interface nodeOperatorDiskVerySimpleAnalyticSolver

contains

  function diskVerySimpleAnalyticSolverConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorDiskVerySimpleAnalyticSolver} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nodeOperatorDiskVerySimpleAnalyticSolver)                :: self
    type            (inputParameters                         ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                 ), pointer       :: cosmologyFunctions_
    class           (stellarPopulationPropertiesClass        ), pointer       :: stellarPopulationProperties_
    class           (stellarFeedbackOutflowsClass            ), pointer       :: stellarFeedbackOutflows_
    class           (starFormationRateDisksClass             ), pointer       :: starFormationRateDisks_
    double precision                                                          :: pruneMassGas                , pruneMassStars
    logical                                                                   :: trackAbundances

    !![
    <inputParameter>
      <name>pruneMassGas</name>
      <defaultValue>0.0d0</defaultValue>
      <description>Gas mass below which the analytic solver will prune a satellite from the tree.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>pruneMassStars</name>
      <defaultValue>0.0d0</defaultValue>
      <description>Stellar mass below which the analytic solver will prune a satellite from the tree.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>trackAbundances</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, also analytically evolve disk and outflowed abundances.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"          name="cosmologyFunctions_"          source="parameters"/>
    <objectBuilder class="stellarPopulationProperties" name="stellarPopulationProperties_" source="parameters"/>
    <objectBuilder class="stellarFeedbackOutflows"     name="stellarFeedbackOutflows_"     source="parameters"/>
    <objectBuilder class="starFormationRateDisks"      name="starFormationRateDisks_"      source="parameters"/>
    !!]
    self=nodeOperatorDiskVerySimpleAnalyticSolver(pruneMassGas,pruneMassStars,trackAbundances,cosmologyFunctions_,stellarPopulationProperties_,stellarFeedbackOutflows_,starFormationRateDisks_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"         />
    <objectDestructor name="stellarPopulationProperties_"/>
    <objectDestructor name="stellarFeedbackOutflows_"    />
    <objectDestructor name="starFormationRateDisks_"     />
    !!]
    return
  end function diskVerySimpleAnalyticSolverConstructorParameters

  function diskVerySimpleAnalyticSolverConstructorInternal(pruneMassGas,pruneMassStars,trackAbundances,cosmologyFunctions_,stellarPopulationProperties_,stellarFeedbackOutflows_,starFormationRateDisks_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorDiskVerySimpleAnalyticSolver} node operator class.
    !!}
    implicit none
    type            (nodeOperatorDiskVerySimpleAnalyticSolver)                        :: self
    double precision                                          , intent(in   )         :: pruneMassGas                , pruneMassStars
    logical                                                   , intent(in   )         :: trackAbundances
    class           (cosmologyFunctionsClass                 ), intent(in   ), target :: cosmologyFunctions_
    class           (stellarPopulationPropertiesClass        ), intent(in   ), target :: stellarPopulationProperties_
    class           (stellarFeedbackOutflowsClass            ), intent(in   ), target :: stellarFeedbackOutflows_
    class           (starFormationRateDisksClass             ), intent(in   ), target :: starFormationRateDisks_
    !![
    <constructorAssign variables="pruneMassGas, pruneMassStars, trackAbundances, *cosmologyFunctions_, *stellarPopulationProperties_, *stellarFeedbackOutflows_, *starFormationRateDisks_"/>
    !!]

    self%active        =.false.
    self%timePresentDay=self%cosmologyFunctions_%cosmicTime(1.0d0)
    return
  end function diskVerySimpleAnalyticSolverConstructorInternal

  subroutine diskVerySimpleAnalyticSolverDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorDiskVerySimpleAnalyticSolver} node operator class.
    !!}
    implicit none
    type(nodeOperatorDiskVerySimpleAnalyticSolver), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"         />
    <objectDestructor name="self%stellarPopulationProperties_"/>
    <objectDestructor name="self%stellarFeedbackOutflows_"    />
    <objectDestructor name="self%starFormationRateDisks_"     />
    !!]
    return
  end subroutine diskVerySimpleAnalyticSolverDestructor

  subroutine diskVerySimpleAnalyticSolverPre(self,node)
    !!{
    Capture the initial state and rates for an analytic solution at the start of each ODE step. If the
    satellite is too small to ever be of interest and will not merge before the present day, prune it
    from the tree, distributing its outflowed mass over future host halos.
    !!}
    use :: Abundances_Structure, only : abundances        , max                , operator(*)
    use :: Galacticus_Nodes    , only : nodeComponentBasic, nodeComponentDisk  , nodeComponentDiskVerySimple, nodeComponentHotHalo, &
         &                              nodeComponentSatellite
    implicit none
    class           (nodeOperatorDiskVerySimpleAnalyticSolver), intent(inout) :: self
    type            (treeNode                                ), intent(inout) :: node
    type            (treeNode                                ), pointer       :: nodeWalk
    class           (nodeComponentBasic                      ), pointer       :: basic                   , basicHost                   , &
         &                                                                       basicParentHost
    class           (nodeComponentDisk                       ), pointer       :: disk
    class           (nodeComponentHotHalo                    ), pointer       :: hotHaloHost
    class           (nodeComponentSatellite                  ), pointer       :: satellite
    double precision                                          , parameter     :: massTolerance        =1.0d-6
    double precision                                                          :: rateFuel                , rateStars                   , &
         &                                                                       rateOutflow             , massStellarAsymptotic       , &
         &                                                                       massOutflowed
    type            (abundances                              ), save          :: rateAbundanceFuel       , rateAbundanceStars          , &
         &                                                                       abundancesOutflowed
    !$omp threadprivate(rateAbundanceFuel,rateAbundanceStars,abundancesOutflowed)

    self%active=.false.
    if (.not.node%isSatellite()) return
    disk => node%disk()
    select type (disk)
    class is (nodeComponentDiskVerySimple)
       if (.not.disk%isInitialized()) return
       self%massGasInitial          =disk%massGas          ()
       self%massStellarInitial      =disk%massStellar      ()
       self%abundancesGasInitial    =disk%abundancesGas    ()
       self%abundancesStellarInitial=disk%abundancesStellar()
       basic       => node%basic    ()
       satellite   => node%satellite()
       self%timeStart=basic%time()
       if (self%massGasInitial <= massTolerance) return
       call diskVerySimpleAnalyticSolverRates(self,node,rateFuel,rateAbundanceFuel,rateStars,rateAbundanceStars,rateOutflow)
       if (rateFuel == 0.0d0 .and. rateStars == 0.0d0 .and. rateOutflow == 0.0d0) return
       ! Compute timescales.
       self%timescaleFuel     =self%massGasInitial/(+rateOutflow-rateFuel                            )
       self%timescaleStellar  =self%massGasInitial/(                              +rateStars         )
       self%timescaleOutflow  =self%massGasInitial/(+rateOutflow                                     )
       self%yieldMassEffective=self%timescaleFuel *(            +rateAbundanceFuel+rateAbundanceStars)
       ! Estimate asymptotic final stellar mass.
       massStellarAsymptotic=self%massStellarInitial+self%massGasInitial*(self%timescaleFuel/self%timescaleStellar)
       ! If the satellite is too small to ever be of interest and will not merge before the present day, prune it.
       if     (                                                  &
            &   massStellarAsymptotic     < self%pruneMassStars  &
            &  .and.                                             &
            &   self%massGasInitial       < self%pruneMassGas    &
            &  .and.                                             &
            &   satellite%timeOfMerging() > self%timePresentDay  &
            & ) then
          ! Distribute analytic outflow over each future host halo.
          nodeWalk => node%parent
          do while (associated(nodeWalk%parent))
             basicHost       => nodeWalk       %basic  (                 )
             hotHaloHost     => nodeWalk       %hotHalo(autoCreate=.true.)
             basicParentHost => nodeWalk%parent%basic  (                 )
             massOutflowed   =  +self%massGasInitial                                                            &
                  &             *(                                                                              &
                  &               +self%timescaleFuel                                                           &
                  &               /self%timescaleOutflow                                                        &
                  &              )                                                                              &
                  &             *(                                                                              &
                  &               +exp(-max(0.0d0,basicHost      %time()-self%timeStart)/self%timescaleFuel)    &
                  &               -exp(-max(0.0d0,basicParentHost%time()-self%timeStart)/self%timescaleFuel)    &
                  &              )
             call hotHaloHost%outflowedMassSet(                             &
                  &                            +hotHaloHost%outflowedMass() &
                  &                            +            massOutflowed   &
                  &                           )
             if (self%trackAbundances) then
                abundancesOutflowed=+(                                       &
                     &                -(                                     &
                     &                  +self%abundancesGasInitial           &
                     &                  +self%yieldMassEffective             &
                     &                  *(                                   &
                     &                    +1.0d0                             &
                     &                    +(                                 &
                     &                      +basicHost        %time()        &
                     &                      -self%timeStart                  &
                     &                     )                                 &
                     &                    /self%timescaleFuel                &
                     &                   )                                   &
                     &                 )                                     &
                     &                *exp(                                  &
                     &                     -(                                &
                     &                       +basicHost      %time()         &
                     &                       -self%timeStart                 &
                     &                      )                                &
                     &                     /self%timescaleFuel               &
                     &                    )                                  &
                     &                +(                                     &
                     &                  +self%abundancesGasInitial           &
                     &                  +self%yieldMassEffective             &
                     &                  *(                                   &
                     &                    +1.0d0                             &
                     &                    +(                                 &
                     &                      +basicParentHost%time()          &
                     &                      -self%timeStart                  &
                     &                     )                                 &
                     &                    /self%timescaleFuel                &
                     &                   )                                   &
                     &                 )                                     &
                     &                *exp(                                  &
                     &                     -(                                &
                     &                       +basicParentHost%time()         &
                     &                       -self%timeStart                 &
                     &                      )                                &
                     &                     /self%timescaleFuel               &
                     &                    )                                  &
                     &               )                                       &
                     &              *self%timescaleFuel                      &
                     &              /self%timescaleOutflow
                call hotHaloHost%outflowedAbundancesSet(                                   &
                     &                                  +hotHaloHost%outflowedAbundances() &
                     &                                  +            abundancesOutflowed   &
                     &                                 )
             end if
             nodeWalk => nodeWalk%parent
          end do
          ! Trigger destruction of this satellite via the standard mechanism. This requires
          ! `mergerTreeEvolveTimestep value="satelliteDestruction"` to be included in the parameter
          ! file (the parametersMigrate helper ensures this), and a satellite component that
          ! supports a settable destructionTime (e.g., `orbiting`).
          if (satellite%destructionTimeIsSettable()) then
             call satellite%destructionTimeSet(0.0d0)
          else
             ! The active satellite component (e.g., `mergeTime`) does not expose a destructionTime
             ! property, so fall back to detaching the node from its host's satellite list. The node
             ! object remains allocated until the tree itself is destroyed.
             call node     %removeFromHost   (     )
          end if
          return
       end if
       self%active=.true.
    end select
    return
  end subroutine diskVerySimpleAnalyticSolverPre

  subroutine diskVerySimpleAnalyticSolverAnalytics(self,node)
    !!{
    Mark the analytically-solvable properties of the disk and host hot halo.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, nodeComponentDiskVerySimple, nodeComponentHotHalo
    implicit none
    class(nodeOperatorDiskVerySimpleAnalyticSolver), intent(inout) :: self
    type (treeNode                                ), intent(inout) :: node
    class(nodeComponentDisk                       ), pointer       :: disk
    class(nodeComponentHotHalo                    ), pointer       :: hotHalo

    if (.not.self%active) return
    disk => node%disk()
    select type (disk)
    class is (nodeComponentDiskVerySimple)
       call disk%massGasAnalytic    ()
       call disk%massStellarAnalytic()
       if (self%trackAbundances) then
          call disk%abundancesGasAnalytic    ()
          call disk%abundancesStellarAnalytic()
       end if
    end select
    hotHalo => node%hotHalo()
    select type (hotHalo)
    type is (nodeComponentHotHalo)
       ! Hot halo does not exist - nothing to mark.
    class default
       call hotHalo%outflowedMassAnalytic()
       if (self%trackAbundances) call hotHalo%outflowedAbundancesAnalytic()
    end select
    return
  end subroutine diskVerySimpleAnalyticSolverAnalytics

  subroutine diskVerySimpleAnalyticSolverSolveAnalytics(self,node,time)
    !!{
    Set the values of the analytically-solvable disk and outflow properties at the given time.
    !!}
    use :: Abundances_Structure, only : abundances
    use :: Galacticus_Nodes    , only : nodeComponentDisk, nodeComponentDiskVerySimple, nodeComponentHotHalo
    implicit none
    class           (nodeOperatorDiskVerySimpleAnalyticSolver), intent(inout) :: self
    type            (treeNode                                ), intent(inout) :: node
    double precision                                          , intent(in   ) :: time
    class           (nodeComponentDisk                       ), pointer       :: disk
    class           (nodeComponentHotHalo                    ), pointer       :: hotHalo
    double precision                                                          :: timeStep         , exponentialFactor      , &
         &                                                                       massGasFinal     , massStellarFinal       , &
         &                                                                       massOutflowed
    type            (abundances                              ), save          :: abundancesGasFinal, abundancesStellarFinal, &
         &                                                                       abundancesOutflowed
    !$omp threadprivate(abundancesGasFinal,abundancesStellarFinal,abundancesOutflowed)

    if (.not.self%active) return
    disk    => node%disk   ()
    hotHalo => node%hotHalo()
    timeStep         =time-self%timeStart
    exponentialFactor=exp(-timeStep/self%timescaleFuel)
    massGasFinal     =                        +self%massGasInitial                                           *       exponentialFactor
    massStellarFinal =+self%massStellarInitial+self%massGasInitial*(self%timescaleFuel/self%timescaleStellar)*(1.0d0-exponentialFactor)
    massOutflowed    =                         self%massGasInitial*(self%timescaleFuel/self%timescaleOutflow)*(1.0d0-exponentialFactor)
    select type (disk)
    class is (nodeComponentDiskVerySimple)
       call disk%massGasSet    (massGasFinal    )
       call disk%massStellarSet(massStellarFinal)
       if (self%trackAbundances) then
          abundancesGasFinal    =+(                              &
               &                   +self%abundancesGasInitial    &
               &                   +self%yieldMassEffective      &
               &                   *timeStep                     &
               &                   /self%timescaleFuel           &
               &                  )                              &
               &                 *exponentialFactor
          abundancesStellarFinal=+self%abundancesStellarInitial  &
               &                 +(                              &
               &                   +(                            &
               &                     +self%abundancesGasInitial  &
               &                     +self%yieldMassEffective    &
               &                    )                            &
               &                   -(                            &
               &                     +self%abundancesGasInitial  &
               &                     +self%yieldMassEffective    &
               &                     *(                          &
               &                       +1.0d0                    &
               &                       +timeStep                 &
               &                       /self%timescaleFuel       &
               &                      )                          &
               &                    )                            &
               &                   *exponentialFactor            &
               &                  )                              &
               &                 *self%timescaleFuel             &
               &                 /self%timescaleStellar
          call disk%abundancesGasSet    (abundancesGasFinal    )
          call disk%abundancesStellarSet(abundancesStellarFinal)
       end if
    end select
    select type (hotHalo)
    type is (nodeComponentHotHalo)
       ! Hot halo does not exist - nothing to set.
    class default
       call hotHalo%outflowedMassSet(massOutflowed)
       if (self%trackAbundances) then
          abundancesOutflowed=+(                              &
               &                +(                            &
               &                  +self%abundancesGasInitial  &
               &                  +self%yieldMassEffective    &
               &                 )                            &
               &                -(                            &
               &                  +self%abundancesGasInitial  &
               &                  +self%yieldMassEffective    &
               &                  *(                          &
               &                    +1.0d0                    &
               &                    +timeStep                 &
               &                    /self%timescaleFuel       &
               &                   )                          &
               &                 )                            &
               &                *exponentialFactor            &
               &               )                              &
               &              *self%timescaleFuel             &
               &              /self%timescaleOutflow
          call hotHalo%outflowedAbundancesSet(abundancesOutflowed)
       end if
    end select
    return
  end subroutine diskVerySimpleAnalyticSolverSolveAnalytics

  subroutine diskVerySimpleAnalyticSolverRates(self,node,fuelMassRate,fuelAbundancesRate,stellarMassRate,stellarAbundancesRate,massOutflowRate)
    !!{
    Compute the rates of change of disk gas and stellar mass, abundances, and outflowing mass that are needed to construct the
    analytic solution.
    !!}
    use :: Abundances_Structure          , only : abundances         , zeroAbundances
    use :: Galacticus_Nodes              , only : nodeComponentDisk  , nodeComponentDiskVerySimple
    use :: Histories                     , only : history
    use :: Stellar_Luminosities_Structure, only : stellarLuminosities
    implicit none
    class           (nodeOperatorDiskVerySimpleAnalyticSolver), intent(inout) :: self
    type            (treeNode                                ), intent(inout) :: node
    double precision                                          , intent(  out) :: fuelMassRate            , stellarMassRate     , &
         &                                                                       massOutflowRate
    type            (abundances                              ), intent(inout) :: fuelAbundancesRate      , stellarAbundancesRate
    class           (nodeComponentDisk                       ), pointer       :: disk
    double precision                                                          :: energyInputRate         , starFormationRate    , &
         &                                                                       rateOutflowEjective     , rateOutflowExpulsive
    type            (abundances                              ), save          :: fuelAbundances
    type            (history                                 )                :: stellarHistoryRate
    type            (stellarLuminosities                     )                :: luminositiesStellarRates
    !$omp threadprivate(fuelAbundances)

    disk => node%disk()
    fuelMassRate         =0.0d0
    stellarMassRate      =0.0d0
    massOutflowRate      =0.0d0
    fuelAbundancesRate   =zeroAbundances
    stellarAbundancesRate=zeroAbundances
    if (disk%massGas() <= 0.0d0) return
    fuelAbundances=disk%abundancesGas()
    call fuelAbundances%massToMassFraction(disk%massGas())
    select type (disk)
    class is (nodeComponentDiskVerySimple)
       starFormationRate=self%starFormationRateDisks_%rate(node)
    end select
    stellarHistoryRate=disk%stellarPropertiesHistory()
    call self%stellarPopulationProperties_%rates(starFormationRate,fuelAbundances,disk,node,stellarHistoryRate, &
         &  stellarMassRate,fuelMassRate,energyInputRate,fuelAbundancesRate,stellarAbundancesRate,              &
         &  luminositiesStellarRates,computeRateLuminosityStellar=.false.)
    call self%stellarFeedbackOutflows_%outflowRate(disk,starFormationRate,energyInputRate,rateOutflowEjective,rateOutflowExpulsive)
    massOutflowRate=rateOutflowEjective
    return
  end subroutine diskVerySimpleAnalyticSolverRates
