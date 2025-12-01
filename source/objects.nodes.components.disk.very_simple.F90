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

!!{
Contains a module that implements a very simple disk component.
!!}

module Node_Component_Disk_Very_Simple
  !!{
  Implements a very simple disk component.
  !!}
  use :: Cosmology_Functions             , only : cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales         , only : darkMatterHaloScaleClass
  use :: Galacticus_Nodes                , only : treeNode
  use :: Math_Exponentiation             , only : fastExponentiator
  use :: Satellite_Merging_Mass_Movements, only : mergerMassMovementsClass
  use :: Star_Formation_Rates_Disks      , only : starFormationRateDisksClass
  use :: Stellar_Feedback_Outflows       , only : stellarFeedbackOutflowsClass
  use :: Stellar_Population_Properties   , only : stellarPopulationPropertiesClass
  implicit none
  private
  public :: Node_Component_Disk_Very_Simple_Scale_Set   , Node_Component_Disk_Very_Simple_Thread_Uninitialize, &
       &    Node_Component_Disk_Very_Simple_Initialize  , Node_Component_Disk_Very_Simple_Pre_Evolve         , &
       &    Node_Component_Disk_Very_Simple_Rates       , Node_Component_Disk_Very_Simple_Analytic_Solver    , &
       &    Node_Component_Disk_Very_Simple_Post_Step   , Node_Component_Disk_Very_Simple_Thread_Initialize  , &
       &    Node_Component_Disk_Very_Simple_State_Store , Node_Component_Disk_Very_Simple_State_Restore

  !![
  <component>
   <class>disk</class>
   <name>verySimple</name>
   <isDefault>false</isDefault>
   <properties>
    <property>
      <name>isInitialized</name>
      <type>logical</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
    </property>
    <property>
      <name>massStellar</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output unitsInSI="massSolar" comment="Mass of stars in the very simple disk."/>
    </property>
    <property>
      <name>abundancesStellar</name>
      <type>abundances</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output unitsInSI="massSolar" comment="Mass of metals in the stellar phase of the standard disk."/>
    </property>
    <property>
      <name>massGas</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
      <output unitsInSI="massSolar" comment="Mass of gas in the very simple disk."/>
    </property>
    <property>
      <name>abundancesGas</name>
      <type>abundances</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
      <output unitsInSI="massSolar" comment="Mass of metals in the gas phase of the standard disk."/>
    </property>
    <property>
      <name>stellarPropertiesHistory</name>
      <type>history</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
    </property>
    <property>
      <name>luminositiesStellar</name>
      <type>stellarLuminosities</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output unitsInSI="luminosityZeroPointAB" comment="Luminosity of disk stars."/>
    </property>
   </properties>
   <bindings>
    <binding method="massBaryonic" function="Node_Component_Disk_Very_Simple_Mass_Baryonic"/>
   </bindings>
   <functions>objects.nodes.components.disk.very_simple.bound_functions.inc</functions>
  </component>
  !!]

  ! Objects used by this component.
  class(cosmologyFunctionsClass         ), pointer :: cosmologyFunctions_
  class(stellarPopulationPropertiesClass), pointer :: stellarPopulationProperties_
  class(darkMatterHaloScaleClass        ), pointer :: darkMatterHaloScale_
  class(stellarFeedbackOutflowsClass    ), pointer :: stellarFeedbackOutflows_
  class(starFormationRateDisksClass     ), pointer :: starFormationRateDisks_
  class(mergerMassMovementsClass        ), pointer :: mergerMassMovements_
  !$omp threadprivate(cosmologyFunctions_,stellarPopulationProperties_,darkMatterHaloScale_,stellarFeedbackOutflows_,starFormationRateDisks_,mergerMassMovements_)

  ! Record of whether to use the simple disk analytic solver.
  logical                             :: useAnalyticSolver
  double precision                    :: pruneMassStars          , pruneMassGas     , &
       &                                 timePresentDay   =-1.0d0

  ! Parameters controlling the physical implementation.
  double precision                    :: scaleAbsoluteMass
  logical                             :: trackAbundances         , trackLuminosities

  ! A threadprivate object used to track to which thread events are attached.
  integer :: thread
  !$omp threadprivate(thread)

contains

  !![
  <nodeComponentInitializationTask>
   <unitName>Node_Component_Disk_Very_Simple_Initialize</unitName>
  </nodeComponentInitializationTask>
  !!]
  subroutine Node_Component_Disk_Very_Simple_Initialize(parameters)
    !!{
    Initializes the tree node very simple disk component module.
    !!}
    use :: Galacticus_Nodes, only : defaultDiskComponent, nodeComponentDiskVerySimple
    use :: Input_Parameters, only : inputParameter      , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters
    type(inputParameters)                :: subParameters

    ! Initialize the module if necessary.
    if (defaultDiskComponent%verySimpleIsActive()) then
       ! Find our parameters.
       subParameters=parameters%subParameters('componentDisk')
       ! Read parameters controlling the physical implementation.
       !![
       <inputParameter>
         <name>scaleAbsoluteMass</name>
         <defaultValue>100.0d0</defaultValue>
         <description>The absolute mass scale below which calculations in the very simple disk component are allowed to become inaccurate.</description>
         <source>parameters</source>
       </inputParameter>
       <inputParameter>
         <name>trackAbundances</name>
         <defaultValue>.false.</defaultValue>
         <description>Specifies whether or not to track abundances in the very simple disk component.</description>
         <source>parameters</source>
       </inputParameter>
       <inputParameter>
         <name>trackLuminosities</name>
         <defaultValue>.false.</defaultValue>
         <description>Specifies whether or not to track stellar luminosities in the very simple disk component.</description>
         <source>parameters</source>
       </inputParameter>
       <inputParameter>
         <name>useAnalyticSolver</name>
         <defaultValue>.false.</defaultValue>
         <description>If true, employ an analytic ODE solver when evolving satellites.</description>
         <source>parameters</source>
       </inputParameter>
       <inputParameter>
         <name>pruneMassGas</name>
         <defaultValue>0.0d0</defaultValue>
         <description>Gas mass below which the analytic solver will prune a galaxy.</description>
         <source>parameters</source>
       </inputParameter>
       <inputParameter>
         <name>pruneMassStars</name>
         <defaultValue>0.0d0</defaultValue>
         <description>Stellar mass below which the analytic solver will prune a galaxy.</description>
         <source>parameters</source>
       </inputParameter>
       !!]
    end if
    return
  end subroutine Node_Component_Disk_Very_Simple_Initialize

  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Node_Component_Disk_Very_Simple_Thread_Initialize</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Node_Component_Disk_Very_Simple_Thread_Initialize(parameters)
    !!{
    Initializes the tree node very simple disk profile module.
    !!}
    use :: Events_Hooks    , only : dependencyDirectionAfter, dependencyRegEx, openMPThreadBindingAtLevel, postEvolveEvent, &
          &                         satelliteMergerEvent
    use :: Galacticus_Nodes, only : defaultDiskComponent
    use :: Input_Parameters, only : inputParameter          , inputParameters
     implicit none
    type(inputParameters), intent(inout) :: parameters
    type(dependencyRegEx), dimension(1)  :: dependencies
    type(inputParameters)                :: subParameters

    if (defaultDiskComponent%verySimpleIsActive()) then
       dependencies(1)=dependencyRegEx(dependencyDirectionAfter,'^remnantStructure:')
       call satelliteMergerEvent%attach(thread,satelliteMerger,openMPThreadBindingAtLevel,label='nodeComponentDiskVerySimple',dependencies=dependencies)
       call postEvolveEvent     %attach(thread,postEvolve     ,openMPThreadBindingAtLevel,label='nodeComponentDiskVerySimple'                          )
       ! Find our parameters.
       subParameters=parameters%subParameters('componentDisk')
       !![
       <objectBuilder class="cosmologyFunctions"          name="cosmologyFunctions_"          source="subParameters"/>
       <objectBuilder class="stellarPopulationProperties" name="stellarPopulationProperties_" source="subParameters"/>
       <objectBuilder class="darkMatterHaloScale"         name="darkMatterHaloScale_"         source="subParameters"/>
       <objectBuilder class="stellarFeedbackOutflows"     name="stellarFeedbackOutflows_"     source="subParameters"/>
       <objectBuilder class="starFormationRateDisks"      name="starFormationRateDisks_"      source="subParameters"/>
       <objectBuilder class="mergerMassMovements"         name="mergerMassMovements_"         source="subParameters"/>
       !!]
       ! If using the analytic solver, find the time at the present day.
       !$omp critical (Node_Component_Disk_Very_Simple_Thread_Initialize)
       if (useAnalyticSolver.and.timePresentDay < 0.0d0) timePresentDay=cosmologyFunctions_%cosmicTime(1.0d0)
       !$omp end critical (Node_Component_Disk_Very_Simple_Thread_Initialize)
    end if
    return
  end subroutine Node_Component_Disk_Very_Simple_Thread_Initialize

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Node_Component_Disk_Very_Simple_Thread_Uninitialize</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Node_Component_Disk_Very_Simple_Thread_Uninitialize()
    !!{
    Uninitializes the tree node very simple disk profile module.
    !!}
    use :: Events_Hooks    , only : postEvolveEvent     , satelliteMergerEvent
    use :: Galacticus_Nodes, only : defaultDiskComponent
    implicit none
    
    if (defaultDiskComponent%verySimpleIsActive()) then
       !![
       <objectDestructor name="cosmologyFunctions_"         />
       <objectDestructor name="stellarPopulationProperties_"/>
       <objectDestructor name="darkMatterHaloScale_"        />
       <objectDestructor name="stellarFeedbackOutflows_"    />
       <objectDestructor name="starFormationRateDisks_"     />
       <objectDestructor name="mergerMassMovements_"        />
       !!]
       if (satelliteMergerEvent%isAttached(thread,satelliteMerger)) call satelliteMergerEvent%detach(thread,satelliteMerger)
       if (postEvolveEvent     %isAttached(thread,postEvolve     )) call postEvolveEvent     %detach(thread,postEvolve     )
    end if
    return
  end subroutine Node_Component_Disk_Very_Simple_Thread_Uninitialize

  !![
  <preEvolveTask>
  <unitName>Node_Component_Disk_Very_Simple_Pre_Evolve</unitName>
  </preEvolveTask>
  !!]
  subroutine Node_Component_Disk_Very_Simple_Pre_Evolve(node)
    !!{
    Ensure the disk has been initialized.
    !!}
    use :: Galacticus_Nodes, only : defaultDiskComponent, nodeComponentDisk, nodeComponentDiskVerySimple, treeNode
    implicit none
    type (treeNode         ), intent(inout), pointer :: node
    class(nodeComponentDisk)               , pointer :: disk

    ! Check if we are the default method.
    if (.not.defaultDiskComponent%verySimpleIsActive()) return
    ! Get the disk component.
    disk => node%disk()
    ! Check if a very simple disk component exists.
    select type (disk)
       class is (nodeComponentDiskVerySimple)
       ! Initialize the disk
       call Node_Component_Disk_Very_Simple_Create(node)
    end select
    return
  end subroutine Node_Component_Disk_Very_Simple_Pre_Evolve

  subroutine postEvolve(self,node)
    !!{
    Catch rounding errors in the very simple disk gas evolution.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDisk, nodeComponentDiskVerySimple, treeNode
    use :: Histories       , only : history
    implicit none
    class(*                 ), intent(inout) :: self
    type (treeNode          ), intent(inout) :: node
    class(nodeComponentDisk ), pointer       :: disk
    class(nodeComponentBasic), pointer       :: basic
    type (history           )                :: stellarPropertiesHistory
    !$GLC attributes unused :: self

    ! Get the disk component.
    disk => node%disk()
    ! Check if a very simple disk component exists.
    select type (disk)
    class is (nodeComponentDiskVerySimple)
       ! Trim the stellar populations properties future history.
       basic => node%basic()
       stellarPropertiesHistory=disk%stellarPropertiesHistory()
       call stellarPropertiesHistory%trim(basic%time())
       call disk%stellarPropertiesHistorySet(stellarPropertiesHistory)
    end select
    return
  end subroutine postEvolve

  !![
  <postStepTask>
    <unitName>Node_Component_Disk_Very_Simple_Post_Step</unitName>
  </postStepTask>
  !!]
  subroutine Node_Component_Disk_Very_Simple_Post_Step(node,status)
    !!{
    Catch rounding errors in the very simple disk gas evolution.
    !!}
    use :: Abundances_Structure          , only : abs                 , zeroAbundances
    use :: Display                       , only : displayMessage      , verbosityLevelWarn
    use :: Galacticus_Nodes              , only : defaultDiskComponent, nodeComponentDisk      , nodeComponentDiskVerySimple, treeNode
    use :: Interface_GSL                 , only : GSL_Success         , GSL_Continue
    use :: ISO_Varying_String            , only : assignment(=)       , operator(//)           , varying_string
    use :: Stellar_Luminosities_Structure, only : abs                 , zeroStellarLuminosities
    use :: String_Handling               , only : operator(//)
    implicit none
    type            (treeNode          ), intent(inout), pointer :: node
    integer                             , intent(inout)          :: status
    class           (nodeComponentDisk )               , pointer :: disk
    double precision                    , save                   :: fractionalErrorMaximum=0.0d0
    double precision                                             :: massDisk                    , fractionalError
    character       (len=20            )                         :: valueString
    type            (varying_string    ), save                   :: message
    !$omp threadprivate(message)

    ! Return immediately if this class is not in use.
    if (.not.defaultDiskComponent%verySimpleIsActive()) return
    ! Get the disk component.
    disk => node%disk()
    ! Check if a very simple disk component exists.
    select type (disk)
    class is (nodeComponentDiskVerySimple)
       ! Note that "status" is not set to failure as these changes in state of the disk should not change any calculation of
       ! differential evolution rates as a negative gas mass was unphysical anyway.
       !
       ! Trap negative gas masses.
       if (disk%massGas() < 0.0d0) then
          ! Check if this exceeds the maximum previously recorded error.
          fractionalError=   abs(disk%massGas    ()) &
               &          /(                         &
               &                 disk%massStellar()  &
               &            +abs(disk%massGas    ()) &
               &           )
          !$omp critical (Very_Simple_Disk_Post_Evolve_Check)
          if (fractionalError > fractionalErrorMaximum) then
             ! Report a warning.
             message='Warning: disk has negative gas mass (fractional error exceeds any previously reported):'//char(10)
             message=message//'  Node index        = '//node%index() //char(10)
             write (valueString,'(e12.6)') disk%massGas()
             message=message//'  Disk gas mass     = '//trim(valueString)//char(10)
             write (valueString,'(e12.6)') disk%massStellar()
             message=message//'  Disk stellar mass = '//trim(valueString)//char(10)
             write (valueString,'(e12.6)') fractionalError
             message=message//'  Error measure     = '//trim(valueString)//char(10)
             if (fractionalErrorMaximum == 0.0d0) then
                ! This is the first time this warning has been issued, so give some extra information.
                message=message//'  Gas mass will be reset to zero (in future cases also).'//char(10)
                message=message//'  Future cases will be reported only when they exceed the previous maximum error measure.'//char(10)
                message=message//'  Negative masses are due to numerical inaccuracy in the ODE solutions.'//char(10)
                message=message//'  If significant, consider using a higher tolerance in the ODE solver.'
             end if
             call displayMessage(message,verbosityLevelWarn)
             ! Store the new maximum fractional error.
             fractionalErrorMaximum=fractionalError
          end if
          !$omp end critical (Very_Simple_Disk_Post_Evolve_Check)
          ! Get the total mass of the disk material
          massDisk= disk%massGas    () &
               &   +disk%massStellar()
          if (massDisk == 0.0d0) then
             call disk%        massStellarSet(                  0.0d0)
             call disk%  abundancesStellarSet(         zeroAbundances)
             call disk%luminositiesStellarSet(zeroStellarLuminosities)
          end if
          ! Reset the gas mass of the disk.
          call disk%      massGasSet(         0.0d0)
          call disk%abundancesGasSet(zeroAbundances)
          ! Indicate that ODE evolution should continue after this state change.
          if (status == GSL_Success) status=GSL_Continue
       end if
    end select
    return
  end subroutine Node_Component_Disk_Very_Simple_Post_Step

  subroutine Node_Component_Disk_Very_Simple_Create(node)
    !!{
    Create properties in a very simple disk component.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, treeNode
    use :: Histories       , only : history
    implicit none
    type   (treeNode         ), intent(inout), target  :: node
    class  (nodeComponentDisk)               , pointer :: disk
    type   (history          )                         :: stellarPropertiesHistory
    logical                                            :: createStellarPropertiesHistory

    ! Get the disk component.
    disk => node%disk()
    ! Exit if already initialized.
    if (disk%isInitialized()) return
    ! Determine which histories must be created.
    stellarPropertiesHistory      =     disk                    %stellarPropertiesHistory()
    createStellarPropertiesHistory=.not.stellarPropertiesHistory%exists                  ()
    call                                stellarPropertiesHistory%destroy                 ()
    ! Create the stellar properties history.
    if (createStellarPropertiesHistory) then
       ! Create the stellar properties history.
       call stellarPopulationProperties_%historyCreate(node,stellarPropertiesHistory)
       call disk%stellarPropertiesHistorySet          (     stellarPropertiesHistory)
    end if
    ! Record that the disk has been initialized.
    call disk%isInitializedSet(.true.)
    return
  end subroutine Node_Component_Disk_Very_Simple_Create

  !![
  <analyticSolverTask>
   <unitName>Node_Component_Disk_Very_Simple_Analytic_Solver</unitName>
  </analyticSolverTask>
  !!]
  subroutine Node_Component_Disk_Very_Simple_Analytic_Solver(node,timeStart,timeEnd,solved)
    use :: Abundances_Structure          , only : abundances        , max                , operator(*)
    use :: Error                         , only : Error_Report
    use :: Galacticus_Nodes              , only : nodeComponentBasic, nodeComponentDisk  , nodeComponentHotHalo, nodeComponentSatellite, &
          &                                       treeNode
    use :: Histories                     , only : history
    use :: Stellar_Luminosities_Structure, only : max               , stellarLuminosities
    implicit none
    type            (treeNode              ), intent(inout), pointer   :: node
    double precision                        , intent(in   )            :: timeStart              , timeEnd
    logical                                 , intent(inout)            :: solved
    type            (treeNode              ),                pointer   :: nodeHost
    class           (nodeComponentBasic    )               , pointer   :: basic                  , basicHost               , &
         &                                                                basicParentHost
    class           (nodeComponentDisk     )               , pointer   :: disk
    class           (nodeComponentHotHalo  )               , pointer   :: hotHalo                , hotHaloHost
    class           (nodeComponentSatellite)               , pointer   :: satellite
    double precision                                       , parameter :: massTolerance=1.0d-6
    double precision                                                   :: massGasInitial         , massStellarInitial      , &
         &                                                                timescaleFuel          , timescaleOutflow        , &
         &                                                                timescaleStellar       , massGasFinal            , &
         &                                                                massStellarFinal       , massOutflowed           , &
         &                                                                rateFuel               , rateStars               , &
         &                                                                rateOutflow            , timeStep                , &
         &                                                                exponentialFactor      , massStellarAsymptotic
    type            (history               )                           :: stellarHistoryRate
    type            (abundances            ), save                     :: rateAbundanceFuel      , rateAbundanceStars      , &
         &                                                                abundancesGasInitial   , abundancesStellarInitial, &
         &                                                                yieldMassEffective     , abundancesStellarFinal  , &
         &                                                                abundancesGasFinal     , abundancesOutflowed
    type            (stellarLuminosities   ), save                     :: luminositiesStellarRates
    !$omp threadprivate(rateAbundanceFuel,rateAbundanceStars,abundancesGasInitial,abundancesStellarInitial,yieldMassEffective,abundancesStellarFinal,abundancesGasFinal,abundancesOutflowed,luminositiesStellarRates)

    if (useAnalyticSolver) then
       disk => node%disk()
       if (node%isSatellite().and.disk%isInitialized()) then
          ! Luminosities can not be computed analytically.
          if (trackLuminosities) call Error_Report('analytic solver does not support stellar luminosity calculation'//{introspection:location})
          ! Calculate analytic solution.
          timeStep          =timeEnd-timeStart
          massGasInitial          =disk%massGas          ()
          massStellarInitial      =disk%massStellar      ()
          abundancesGasInitial    =disk%abundancesGas    ()
          abundancesStellarInitial=disk%abundancesStellar()
          if (massGasInitial > massTolerance) then
             hotHalo   => node%hotHalo  ()
             satellite => node%satellite()
             call Node_Component_Disk_Very_Simple_Rates(node,rateFuel,rateAbundanceFuel,rateStars,rateAbundanceStars,rateOutflow,stellarHistoryRate,luminositiesStellarRates)
             ! If any rates are zero, return without an analytic solution
             if (rateFuel == 0.0d0 .and. rateStars == 0.0d0 .and. rateOutflow == 0.0d0) return
             ! Compute timescales.
             timescaleFuel     =massGasInitial/(+rateOutflow-rateFuel                            )
             timescaleStellar  =massGasInitial/(                              +rateStars         )
             timescaleOutflow  =massGasInitial/(+rateOutflow                                     )
             yieldMassEffective=timescaleFuel *(            +rateAbundanceFuel+rateAbundanceStars)
             ! Estimate asymptotic final stellar mass
             massStellarAsymptotic=massStellarInitial+massGasInitial*(timescaleFuel/timescaleStellar)
             ! Check if this galaxy can be removed because it will never be of sufficient mass to be interesting.
             if     (                                            &
                  &   massStellarAsymptotic     < pruneMassStars &
                  &  .and.                                       &
                  &   massGasInitial            < pruneMassGas   &
                  &  .and.                                       &
                  &   satellite%timeOfMerging() > timePresentDay &
                  & ) then
                ! Galaxy is too small to care about. Add gas and abundances which will outflow from this satellite to its future host halos.
                nodeHost => node%parent
                do while (associated(nodeHost%parent))
                   basicHost       => nodeHost%basic       (                 )
                   hotHaloHost     => nodeHost%hotHalo     (autoCreate=.true.)
                   basicParentHost => nodeHost%parent%basic(                 )
                   massOutflowed   =  +massGasInitial                                                    &
                        &             *(                                                                 &
                        &               +timescaleFuel                                                   &
                        &               /timescaleOutflow                                                &
                        &              )                                                                 &
                        &             *(                                                                 &
                        &               +exp(-max(0.0d0,basicHost      %time()-timeStart)/timescaleFuel) &
                        &               -exp(-max(0.0d0,basicParentHost%time()-timeStart)/timescaleFuel) &
                        &              )
                   call hotHaloHost%outflowedMassSet(                             &
                        &                            +hotHaloHost%outflowedMass() &
                        &                            +            massOutflowed   &
                        &                           )
                   if (trackAbundances) then
                      abundancesOutflowed=+(                                &
                           &                -(                              &
                           &                  +abundancesGasInitial         &
                           &                  +yieldMassEffective           &
                           &                  *(                            &
                           &                    +1.0d0                      &
                           &                    +(                          &
                           &                      +basicHost        %time() &
                           &                      -timeStart                &
                           &                     )                          &
                           &                    /timescaleFuel              &
                           &                   )                            &
                           &                 )                              &
                           &                *exp(                           &
                           &                     -(                         &
                           &                       +basicHost      %time()  &
                           &                       -timeStart               &
                           &                      )                         &
                           &                     /timescaleFuel             &
                           &                    )                           &
                           &                +(                              &
                           &                  +abundancesGasInitial         &
                           &                  +yieldMassEffective           &
                           &                  *(                            &
                           &                    +1.0d0                      &
                           &                    +(                          &
                           &                      +basicParentHost%time()   &
                           &                      -timeStart                &
                           &                     )                          &
                           &                    /timescaleFuel              &
                           &                   )                            &
                           &                 )                              &
                           &                *exp(                           &
                           &                     -(                         &
                           &                       +basicParentHost%time()  &
                           &                       -timeStart               &
                           &                      )                         &
                           &                     /timescaleFuel             &
                           &                    )                           &
                           &               )                                &
                           &              *timescaleFuel                    &
                           &              /timescaleOutflow
                      call hotHaloHost%outflowedAbundancesSet(                                   &
                           &                                  +hotHaloHost%outflowedAbundances() &
                           &                                  +            abundancesOutflowed   &
                           &                                 )
                   end if
                   nodeHost => nodeHost%parent
                end do
                ! Remove the node from the host and destroy it.
                call node%removeFromHost()
                call node%destroy       ()
                deallocate(node)
                nullify   (node)
             else
                ! Galaxy is sufficiently large (or will merge), so simply process it to the end time.
                exponentialFactor=exp(-timeStep/timescaleFuel)
                massGasFinal     =                   +massGasInitial                                 *       exponentialFactor
                massStellarFinal =+massStellarInitial+massGasInitial*(timescaleFuel/timescaleStellar)*(1.0d0-exponentialFactor)
                massOutflowed    =                   +massGasInitial*(timescaleFuel/timescaleOutflow)*(1.0d0-exponentialFactor)
                call hotHalo%outflowedMassSet(massOutflowed   )
                call disk   %massGasSet      (massGasFinal    )
                call disk   %massStellarSet  (massStellarFinal)
                if (trackAbundances) then
                   abundancesGasFinal    =+(                        &
                        &                   +abundancesGasInitial   &
                        &                   +yieldMassEffective     &
                        &                   *timeStep               &
                        &                   /timescaleFuel          &
                        &                  )                        &
                        &                 *exponentialFactor
                   abundancesStellarFinal=+abundancesStellarFinal   &
                        &                 +(                        &
                        &                   +(                      &
                        &                     +abundancesGasInitial &
                        &                     +yieldMassEffective   &
                        &                    )                      &
                        &                   +(                      &
                        &                     +abundancesGasInitial &
                        &                     +yieldMassEffective   &
                        &                     *(                    &
                        &                       +1.0d0              &
                        &                       +timeStep           &
                        &                       /timescaleFuel      &
                        &                      )                    &
                        &                    )                      &
                        &                   *exponentialFactor      &
                        &                  )                        &
                        &                 *timescaleFuel            &
                        &                 /timescaleStellar
                   abundancesOutflowed   =+(                        &
                        &                   +(                      &
                        &                     +abundancesGasInitial &
                        &                     +yieldMassEffective   &
                        &                    )                      &
                        &                   +(                      &
                        &                     +abundancesGasInitial &
                        &                     +yieldMassEffective   &
                        &                     *(                    &
                        &                       +1.0d0              &
                        &                       +timeStep           &
                        &                       /timescaleFuel      &
                        &                      )                    &
                        &                    )                      &
                        &                   *exponentialFactor      &
                        &                  )                        &
                        &                 *timescaleFuel            &
                        &                 /timescaleOutflow
                   call hotHalo%outflowedAbundancesSet(abundancesOutflowed   )
                   call disk   %abundancesGasSet      (abundancesGasFinal    )
                   call disk   %abundancesStellarSet  (abundancesStellarFinal)
                end if
             end if
          end if
          ! Update time and merging times.
          if (associated(node)) then
             basic     => node%basic    ()
             satellite => node%satellite()
             call basic%timeSet(timeEnd)
          end if
          ! Record that we solved this system analytically.
          solved=.true.
       end if
    end if
    return
  end subroutine Node_Component_Disk_Very_Simple_Analytic_Solver

  subroutine Node_Component_Disk_Very_Simple_Rates(node,fuelMassRate,fuelAbundancesRate,stellarMassRate,stellarAbundancesRate,massOutflowRate,stellarHistoryRate,luminositiesStellarRates)
    !!{
    Compute rates.
    !!}
    use :: Abundances_Structure          , only : abundances         , zeroAbundances
    use :: Galacticus_Nodes              , only : nodeComponentDisk  , nodeComponentDiskVerySimple, treeNode
    use :: Histories                     , only : history
    use :: Stellar_Luminosities_Structure, only : stellarLuminosities
    implicit none
    type            (treeNode           ), intent(inout), pointer :: node
    type            (history            ), intent(inout)          :: stellarHistoryRate
    double precision                     , intent(  out)          :: fuelMassRate            , stellarMassRate       , &
         &                                                           massOutflowRate
    type            (abundances         ), intent(inout)          :: fuelAbundancesRate      , stellarAbundancesRate
    type            (stellarLuminosities), intent(inout)          :: luminositiesStellarRates
    class           (nodeComponentDisk  )               , pointer :: disk
    double precision                                              :: energyInputRate         , starFormationRate     , &
         &                                                           rateOutflowEjective     , rateOutflowExpulsive
    type            (abundances         )               , save    :: fuelAbundances
    !$omp threadprivate (fuelAbundances)

    ! Get the disk.
    disk => node%disk()
    ! Initialize to zero rates.
    fuelMassRate         =0.0d0
    stellarMassRate      =0.0d0
    massOutflowRate      =0.0d0
    fuelAbundancesRate   =zeroAbundances
    stellarAbundancesRate=zeroAbundances
    ! Check for a realistic disk, return immediately if disk is unphysical.
    if (disk%massGas() <= 0.0d0) return
    ! Compute fuel abundances.
    fuelAbundances=disk%abundancesGas()
    call fuelAbundances%massToMassFraction(disk%massGas())
    ! Compute the star formation rate.
    select type (disk)
    class is (nodeComponentDiskVerySimple)
       starFormationRate=starFormationRateDisks_%rate(node)
    end select
    ! Find rates of change of stellar mass, and gas mass.
    stellarHistoryRate=disk%stellarPropertiesHistory()
    call stellarPopulationProperties_%rates(starFormationRate,fuelAbundances,disk,node,stellarHistoryRate&
            &,stellarMassRate,fuelMassRate,energyInputRate,fuelAbundancesRate,stellarAbundancesRate,luminositiesStellarRates,computeRateLuminosityStellar=.true.)
    ! Find rate of outflow of material from the disk.
    call stellarFeedbackOutflows_%outflowRate(disk,starFormationRate,energyInputRate,rateOutflowEjective,rateOutflowExpulsive)
    massOutflowRate=rateOutflowEjective
    return
  end subroutine Node_Component_Disk_Very_Simple_Rates

  !![
  <scaleSetTask>
   <unitName>Node_Component_Disk_Very_Simple_Scale_Set</unitName>
  </scaleSetTask>
  !!]
  subroutine Node_Component_Disk_Very_Simple_Scale_Set(node)
    !!{
    Set scales for properties of {\normalfont \ttfamily node}.
    !!}
    use :: Abundances_Structure          , only : abs                 , abundances       , max                        , unitAbundances         , &
          &                                       zeroAbundances
    use :: Galacticus_Nodes              , only : defaultDiskComponent, nodeComponentDisk, nodeComponentDiskVerySimple, treeNode
    use :: Histories                     , only : history
    use :: Stellar_Luminosities_Structure, only : abs                 , max              , stellarLuminosities        , unitStellarLuminosities
    implicit none
    type            (treeNode           ), intent(inout), pointer :: node
    class           (nodeComponentDisk  )               , pointer :: disk
    double precision                     , parameter              :: luminosityMinimum             =1.0d0
    double precision                                              :: mass
    type            (history            )                         :: stellarPopulationHistoryScales
    type            (abundances         )                         :: abundancesTotal
    type            (stellarLuminosities)                         :: stellarLuminositiesScale

    ! Check if we are the default method.
    if (.not.defaultDiskComponent%verySimpleIsActive()) return
    ! Get the disk component.
    disk => node%disk()
    ! Check if a very simple disk component exists.
    select type (disk)
    class is (nodeComponentDiskVerySimple)
       ! Set scale for gas and stellar mass.
       mass=disk%massGas()+disk%massStellar()
       call disk%massGasScale    (max(mass,scaleAbsoluteMass))
       call disk%massStellarScale(max(mass,scaleAbsoluteMass))
       ! Set scale for gas and stellar abundances.
       abundancesTotal=disk%abundancesGas()+disk%abundancesStellar()
       call disk%abundancesGasScale    (max(abundancesTotal,unitAbundances*scaleAbsoluteMass))
       call disk%abundancesStellarScale(max(abundancesTotal,unitAbundances*scaleAbsoluteMass))
       ! Set scales for stellar population properties and star formation histories.
       stellarPopulationHistoryScales=disk%stellarPropertiesHistory()
       call stellarPopulationProperties_%scales(disk%massStellar(),zeroAbundances,stellarPopulationHistoryScales)
       call disk%stellarPropertiesHistoryScale (                                  stellarPopulationHistoryScales)
       call stellarPopulationHistoryScales%destroy()
       ! Set scale for stellar luminosities.
       stellarLuminositiesScale=max(                                  &
            &                       +abs(disk%luminositiesStellar()), &
            &                       +unitStellarLuminosities          &
            &                       *luminosityMinimum                &
            &                      )
       call stellarLuminositiesScale%truncate                (disk                    %luminositiesStellar())
       call disk                    %luminositiesStellarScale(stellarLuminositiesScale                      )
    end select
    return
  end subroutine Node_Component_Disk_Very_Simple_Scale_Set

  subroutine satelliteMerger(self,node)
    !!{
    Transfer any very simple disk associated with {\normalfont \ttfamily node} to its host halo.
    !!}
    use :: Abundances_Structure            , only : zeroAbundances
    use :: Error                           , only : Error_Report
    use :: Galacticus_Nodes                , only : nodeComponentDisk      , nodeComponentDiskVerySimple, nodeComponentSpheroid           , treeNode
    use :: Histories                       , only : history
    use :: Satellite_Merging_Mass_Movements, only : destinationMergerDisk  , destinationMergerSpheroid  , enumerationDestinationMergerType
    use :: Stellar_Luminosities_Structure  , only : zeroStellarLuminosities
    implicit none
    class  (*                               ), intent(inout) :: self
    type   (treeNode                        ), intent(inout) :: node
    type   (treeNode                        ), pointer       :: nodeHost
    class  (nodeComponentDisk               ), pointer       :: diskHost               , disk
    class  (nodeComponentSpheroid           ), pointer       :: spheroidHost
    type   (enumerationDestinationMergerType)                :: destinationGasSatellite, destinationGasHost       , &
         &                                                      destinationStarsHost   , destinationStarsSatellite
    type   (history                         )                :: historyHost            , historyNode
    logical                                                  :: mergerIsMajor
    !$GLC attributes unused :: self

    ! Check that the disk is of the verySimple class.
    disk => node%disk()
    select type (disk)
    class is (nodeComponentDiskVerySimple)
       ! Find the node to merge with and its disk component (and spheroid if necessary).
       nodeHost => node    %mergesWith(                 )
       diskHost => nodeHost%disk      (autoCreate=.true.)
       ! Get mass movement descriptors.
       call mergerMassMovements_%get(node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
       if     (                                                        &
            &   destinationGasSatellite   == destinationMergerSpheroid &
            &  .or.                                                    &
            &   destinationStarsSatellite == destinationMergerSpheroid &
            & )                                                        &
            & spheroidHost => nodeHost%spheroid(autoCreate=.true.)
       ! Move the gas component of the very simple disk to the host.
       select case (destinationGasSatellite%ID)
       case (destinationMergerDisk%ID)
          call diskHost    %massGasSet          (                              &
               &                                  diskHost    %      massGas() &
               &                                 +disk        %      massGas() &
               &                                )
          call diskHost    %abundancesGasSet    (                              &
               &                                  diskHost    %abundancesGas() &
               &                                 +disk        %abundancesGas() &
               &                                )
       case (destinationMergerSpheroid%ID)
          call spheroidHost%massGasSet          (                              &
               &                                  spheroidHost%massGas      () &
               &                                 +disk        %massGas      () &
               &                                )
          call spheroidHost%abundancesGasSet    (                              &
               &                                  spheroidHost%abundancesGas() &
               &                                 +disk        %abundancesGas() &
               &                                )
       case default
          call Error_Report(                                    &
               &            'unrecognized movesTo descriptor'// &
               &            {introspection:location}            &
               &           )
       end select
       call    disk%massGasSet                  (                                        &
            &                                                             0.0d0          &
            &                                   )
       call    disk%abundancesGasSet            (                                        &
            &                                                    zeroAbundances          &
            &                                   )
       ! Move the stellar component of the very simple disk to the host.
       select case (destinationStarsSatellite%ID)
       case (destinationMergerDisk%ID)
          call diskHost    %massStellarSet        (                                    &
               &                                    diskHost    %        massStellar() &
               &                                   +disk        %        massStellar() &
               &                                  )
          call diskHost    %abundancesStellarSet  (                                    &
               &                                    diskHost    %  abundancesStellar() &
               &                                   +disk        %  abundancesStellar() &
               &                                  )
          call diskHost    %luminositiesStellarSet(                                    &
               &                                    diskHost    %luminositiesStellar() &
               &                                   +disk        %luminositiesStellar() &
               &                                  )
          ! Also add stellar properties histories.
          historyNode=disk    %stellarPropertiesHistory()
          historyHost=diskHost%stellarPropertiesHistory()
          call historyHost%interpolatedIncrement      (historyNode)
          call historyNode%reset                      (           )
          call diskHost   %stellarPropertiesHistorySet(historyHost)
          call disk       %stellarPropertiesHistorySet(historyNode)
        case (destinationMergerSpheroid%ID)
          call spheroidHost%massStellarSet        (                                    &
               &                                    spheroidHost%  massStellar      () &
               &                                   +disk        %  massStellar      () &
               &                                  )
          call spheroidHost%abundancesStellarSet  (                                    &
               &                                    spheroidHost%  abundancesStellar() &
               &                                   +disk        %  abundancesStellar() &
               &                                  )
          call spheroidHost%luminositiesStellarSet(                                    &
               &                                    spheroidHost%luminositiesStellar() &
               &                                   +disk        %luminositiesStellar() &
               &                                  )
          ! Also add stellar properties histories.
          historyNode=disk    %stellarPropertiesHistory()
          historyHost=spheroidHost%stellarPropertiesHistory()
          call historyHost %interpolatedIncrement      (historyNode)
          call historyNode %reset                      (           )
          call spheroidHost%stellarPropertiesHistorySet(historyHost)
          call disk        %stellarPropertiesHistorySet(historyNode)
       case default
          call Error_Report(                                    &
               &            'unrecognized movesTo descriptor'// &
               &            {introspection:location}            &
               &           )
       end select
       call    disk%         massStellarSet(                       &
            &                                                0.0d0 &
            &                              )
       call    disk%  abundancesStellarSet(                        &
            &                                       zeroAbundances &
            &                             )
       call    disk%luminositiesStellarSet(                        &
            &                              zeroStellarLuminosities &
            &                             )
    end select
    return
  end subroutine satelliteMerger

  !![
  <stateStoreTask>
   <unitName>Node_Component_Disk_Very_Simple_State_Store</unitName>
  </stateStoreTask>
  !!]
  subroutine Node_Component_Disk_Very_Simple_State_Store(stateFile,gslStateFile,stateOperationID)
    !!{
    Store object state,
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Storing state for: componentDisk -> verySimple',verbosity=verbosityLevelInfo)
    !![
    <stateStore variables="cosmologyFunctions_ stellarPopulationProperties_ darkMatterHaloScale_ stellarFeedbackOutflows_ starFormationRateDisks_ mergerMassMovements_"/>
    !!]
    return
  end subroutine Node_Component_Disk_Very_Simple_State_Store

  !![
  <stateRetrieveTask>
   <unitName>Node_Component_Disk_Very_Simple_State_Restore</unitName>
  </stateRetrieveTask>
  !!]
  subroutine Node_Component_Disk_Very_Simple_State_Restore(stateFile,gslStateFile,stateOperationID)
    !!{
    Retrieve object state.
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Retrieving state for: componentDisk -> verySimple',verbosity=verbosityLevelInfo)
    !![
    <stateRestore variables="cosmologyFunctions_ stellarPopulationProperties_ darkMatterHaloScale_ stellarFeedbackOutflows_ starFormationRateDisks_ mergerMassMovements_"/>
    !!]
    return
  end subroutine Node_Component_Disk_Very_Simple_State_Restore

end module Node_Component_Disk_Very_Simple
