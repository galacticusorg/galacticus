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

!% Contains a module that implements a very simple disk component.

module Node_Component_Disk_Very_Simple
  !% Implements a very simple disk component.
  use :: Cosmology_Functions             , only : cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales         , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO        , only : darkMatterProfileDMOClass
  use :: Galacticus_Nodes                , only : treeNode
  use :: Math_Exponentiation             , only : fastExponentiator
  use :: Satellite_Merging_Mass_Movements, only : mergerMassMovementsClass
  use :: Star_Formation_Feedback_Disks   , only : starFormationFeedbackDisksClass
  use :: Star_Formation_Rates_Disks      , only : starFormationRateDisksClass
  use :: Stellar_Population_Properties   , only : stellarPopulationPropertiesClass
  implicit none
  private
  public :: Node_Component_Disk_Very_Simple_Scale_Set   , Node_Component_Disk_Very_Simple_Thread_Uninitialize, &
       &    Node_Component_Disk_Very_Simple_Initialize  , Node_Component_Disk_Very_Simple_Pre_Evolve         , &
       &    Node_Component_Disk_Very_Simple_Rates       , Node_Component_Disk_Very_Simple_Analytic_Solver    , &
       &    Node_Component_Disk_Very_Simple_Post_Step   , Node_Component_Disk_Very_Simple_Thread_Initialize  , &
       &    Node_Component_Disk_Very_Simple_Rate_Compute

  !# <component>
  !#  <class>disk</class>
  !#  <name>verySimple</name>
  !#  <isDefault>false</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>isInitialized</name>
  !#     <type>logical</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </property>
  !#   <property>
  !#     <name>massStellar</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of stars in the very simple disk."/>
  !#   </property>
  !#   <property>
  !#     <name>abundancesStellar</name>
  !#     <type>abundances</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of metals in the stellar phase of the standard disk."/>
  !#   </property>
  !#   <property>
  !#     <name>massGas</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" makeGeneric="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of gas in the very simple disk."/>
  !#   </property>
  !#   <property>
  !#     <name>abundancesGas</name>
  !#     <type>abundances</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" makeGeneric="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of metals in the gas phase of the standard disk."/>
  !#   </property>
  !#   <property>
  !#     <name>stellarPropertiesHistory</name>
  !#     <type>history</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#   </property>
  !#   <property>
  !#     <name>luminositiesStellar</name>
  !#     <type>stellarLuminosities</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <output unitsInSI="luminosityZeroPointAB" comment="Luminosity of disk stars."/>
  !#   </property>
  !#  </properties>
  !#  <bindings>
  !#   <binding method="attachPipe"     function="Node_Component_Disk_Very_Simple_Attach_Pipe" description="Attach pipes to the very simple disk component." bindsTo="component" returnType="\void" arguments="" />
  !#   <binding method="enclosedMass"   function="Node_Component_Disk_Very_Simple_Enclosed_Mass"   bindsTo="component" />
  !#  </bindings>
  !#  <functions>objects.nodes.components.disk.very_simple.bound_functions.inc</functions>
  !# </component>

  ! Objects used by this component.
  class(cosmologyFunctionsClass         ), pointer :: cosmologyFunctions_
  class(stellarPopulationPropertiesClass), pointer :: stellarPopulationProperties_
  class(darkMatterHaloScaleClass        ), pointer :: darkMatterHaloScale_
  class(starFormationFeedbackDisksClass ), pointer :: starFormationFeedbackDisks_
  class(starFormationRateDisksClass     ), pointer :: starFormationRateDisks_
  class(darkMatterProfileDMOClass       ), pointer :: darkMatterProfileDMO_
  class(mergerMassMovementsClass        ), pointer :: mergerMassMovements_
  !$omp threadprivate(cosmologyFunctions_,stellarPopulationProperties_,darkMatterHaloScale_,starFormationFeedbackDisks_,starFormationRateDisks_,darkMatterProfileDMO_,mergerMassMovements_)

  ! Record of whether to use the simple disk analytic solver.
  logical                             :: diskVerySimpleUseAnalyticSolver
  double precision                    :: diskVerySimpleAnalyticSolverPruneMassStars         , diskVerySimpleAnalyticSolverPruneMassGas, &
       &                                 timePresentDay                              =-1.0d0

  ! Parameters controlling the physical implementation.
  double precision                    :: diskOutflowTimescaleMinimum                        , diskVerySimpleMassScaleAbsolute
  logical                             :: diskVerySimpleTrackAbundances                      , diskVerySimpleTrackLuminosities

contains

  !# <nodeComponentInitializationTask>
  !#  <unitName>Node_Component_Disk_Very_Simple_Initialize</unitName>
  !# </nodeComponentInitializationTask>
  subroutine Node_Component_Disk_Very_Simple_Initialize(parameters_)
    !% Initializes the tree node very simple disk component module.
    use :: Galacticus_Nodes, only : defaultDiskComponent, nodeComponentDiskVerySimple
    use :: Input_Parameters, only : inputParameter      , inputParameters
    implicit none
    type(inputParameters            ), intent(inout) :: parameters_
    type(nodeComponentDiskVerySimple)                :: diskVerySimpleComponent

    ! Initialize the module if necessary.
    if (defaultDiskComponent%verySimpleIsActive()) then
       ! Read parameters controlling the physical implementation.
       !# <inputParameter>
       !#   <name>diskVerySimpleMassScaleAbsolute</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultValue>100.0d0</defaultValue>
       !#   <description>The absolute mass scale below which calculations in the very simple disk component are allowed to become inaccurate.</description>
       !#   <source>parameters_</source>
       !#   <type>double</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>diskOutflowTimescaleMinimum</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultValue>1.0d-3</defaultValue>
       !#   <description>The minimum timescale (in units of the halo dynamical time) on which outflows may deplete gas in the disk.</description>
       !#   <source>parameters_</source>
       !#   <type>double</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>diskVerySimpleTrackAbundances</name>
       !#   <cardinality>0..1</cardinality>
       !#   <defaultValue>.false.</defaultValue>
       !#   <description>Specifies whether or not to track abundances in the very simple disk component.</description>
       !#   <source>parameters_</source>
       !#   <type>boolean</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>diskVerySimpleTrackLuminosities</name>
       !#   <cardinality>0..1</cardinality>
       !#   <defaultValue>.false.</defaultValue>
       !#   <description>Specifies whether or not to track stellar luminosities in the very simple disk component.</description>
       !#   <source>parameters_</source>
       !#   <type>boolean</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>diskVerySimpleUseAnalyticSolver</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultValue>.false.</defaultValue>
       !#   <description>If true, employ an analytic ODE solver when evolving satellites.</description>
       !#   <source>parameters_</source>
       !#   <type>boolean</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>diskVerySimpleAnalyticSolverPruneMassGas</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultValue>0.0d0</defaultValue>
       !#   <description>Gas mass below which the analytic solver will prune a galaxy.</description>
       !#   <source>parameters_</source>
       !#   <type>double</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>diskVerySimpleAnalyticSolverPruneMassStars</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultValue>0.0d0</defaultValue>
       !#   <description>Stellar mass below which the analytic solver will prune a galaxy.</description>
       !#   <source>parameters_</source>
       !#   <type>double</type>
       !# </inputParameter>
       ! Attach the cooling mass pipe from the hot halo component.
       call diskVerySimpleComponent%attachPipe()
    end if
    return
  end subroutine Node_Component_Disk_Very_Simple_Initialize

  !# <nodeComponentThreadInitializationTask>
  !#  <unitName>Node_Component_Disk_Very_Simple_Thread_Initialize</unitName>
  !# </nodeComponentThreadInitializationTask>
  subroutine Node_Component_Disk_Very_Simple_Thread_Initialize(parameters_)
    !% Initializes the tree node very simple disk profile module.
    use :: Events_Hooks                        , only : satelliteMergerEvent       , postEvolveEvent, openMPThreadBindingAtLevel, dependencyRegEx, &
         &                                              dependencyDirectionAfter
    use :: Galacticus_Nodes                    , only : defaultDiskComponent
    use :: Input_Parameters                    , only : inputParameter             , inputParameters
     implicit none
    type(inputParameters), intent(inout) :: parameters_
    type(dependencyRegEx), dimension(1)  :: dependencies

    if (defaultDiskComponent%verySimpleIsActive()) then
       dependencies(1)=dependencyRegEx(dependencyDirectionAfter,'^remnantStructure:')
       call satelliteMergerEvent%attach(defaultDiskComponent,satelliteMerger,openMPThreadBindingAtLevel,label='nodeComponentDiskVerySimple',dependencies=dependencies)
       call postEvolveEvent     %attach(defaultDiskComponent,postEvolve     ,openMPThreadBindingAtLevel,label='nodeComponentDiskVerySimple'                          )
       !# <objectBuilder class="cosmologyFunctions"          name="cosmologyFunctions_"          source="parameters_"/>
       !# <objectBuilder class="stellarPopulationProperties" name="stellarPopulationProperties_" source="parameters_"/>
       !# <objectBuilder class="darkMatterHaloScale"         name="darkMatterHaloScale_"         source="parameters_"/>
       !# <objectBuilder class="darkMatterProfileDMO"        name="darkMatterProfileDMO_"        source="parameters_"/>
       !# <objectBuilder class="starFormationFeedbackDisks"  name="starFormationFeedbackDisks_"  source="parameters_"/>
       !# <objectBuilder class="starFormationRateDisks"      name="starFormationRateDisks_"      source="parameters_"/>
       !# <objectBuilder class="mergerMassMovements"         name="mergerMassMovements_"         source="parameters_"/>
       ! If using the analytic solver, find the time at the present day.
       !$omp critical (Node_Component_Disk_Very_Simple_Thread_Initialize)
       if (diskVerySimpleUseAnalyticSolver.and.timePresentDay < 0.0d0) timePresentDay=cosmologyFunctions_%cosmicTime(1.0d0)
       !$omp end critical (Node_Component_Disk_Very_Simple_Thread_Initialize)
    end if
    return
  end subroutine Node_Component_Disk_Very_Simple_Thread_Initialize

  !# <nodeComponentThreadUninitializationTask>
  !#  <unitName>Node_Component_Disk_Very_Simple_Thread_Uninitialize</unitName>
  !# </nodeComponentThreadUninitializationTask>
  subroutine Node_Component_Disk_Very_Simple_Thread_Uninitialize()
    !% Uninitializes the tree node very simple disk profile module.
    use :: Events_Hooks    , only : satelliteMergerEvent, postEvolveEvent
    use :: Galacticus_Nodes, only : defaultDiskComponent
    implicit none
    
    if (defaultDiskComponent%verySimpleIsActive()) then
       !# <objectDestructor name="cosmologyFunctions_"         />
       !# <objectDestructor name="stellarPopulationProperties_"/>
       !# <objectDestructor name="darkMatterHaloScale_"        />
       !# <objectDestructor name="darkMatterProfileDMO_"       />
       !# <objectDestructor name="starFormationFeedbackDisks_" />
       !# <objectDestructor name="starFormationRateDisks_"     />
       !# <objectDestructor name="mergerMassMovements_"        />
       call satelliteMergerEvent%detach(defaultDiskComponent,satelliteMerger)
       call postEvolveEvent     %detach(defaultDiskComponent,postEvolve     )
    end if
    return
  end subroutine Node_Component_Disk_Very_Simple_Thread_Uninitialize

  !# <preEvolveTask>
  !# <unitName>Node_Component_Disk_Very_Simple_Pre_Evolve</unitName>
  !# </preEvolveTask>
  subroutine Node_Component_Disk_Very_Simple_Pre_Evolve(node)
    !% Ensure the disk has been initialized.
    use :: Galacticus_Nodes, only : nodeComponentDisk, nodeComponentDiskVerySimple, treeNode, defaultDiskComponent
    implicit none
    type (treeNode         ), intent(inout), pointer :: node
    class(nodeComponentDisk)               , pointer :: disk

    ! Check if we are the default method.
    if (.not.defaultDiskComponent%verySimpleIsActive()) return
    ! Get the disk component.
    disk => node%disk()
    ! Check if an exponential disk component exists.
    select type (disk)
       class is (nodeComponentDiskVerySimple)
       ! Initialize the disk
       call Node_Component_Disk_Very_Simple_Create(node)
    end select
    return
  end subroutine Node_Component_Disk_Very_Simple_Pre_Evolve

  subroutine postEvolve(self,node)
    !% Catch rounding errors in the very simple disk gas evolution.
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

  !# <postStepTask>
  !# <unitName>Node_Component_Disk_Very_Simple_Post_Step</unitName>
  !# </postStepTask>
  subroutine Node_Component_Disk_Very_Simple_Post_Step(node,status)
    !% Catch rounding errors in the very simple disk gas evolution.
    use :: Abundances_Structure          , only : abs                       , zeroAbundances
    use :: Galacticus_Display            , only : Galacticus_Display_Message, verbosityWarn
    use :: Galacticus_Nodes              , only : nodeComponentDisk         , nodeComponentDiskVerySimple, treeNode    , defaultDiskComponent
    use :: Interface_GSL                 , only : GSL_Failure
    use :: ISO_Varying_String            , only : varying_string            , assignment(=)              , operator(//)
    use :: Stellar_Luminosities_Structure, only : abs                       , zeroStellarLuminosities
    use :: String_Handling               , only : operator(//)
    implicit none
    type            (treeNode          ), intent(inout), pointer :: node
    integer                             , intent(inout)          :: status
    class           (nodeComponentDisk )               , pointer :: disk
    double precision                    , save                   :: fractionalErrorMaximum=0.0d0
    double precision                                             :: diskMass                    , fractionalError
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
             call Galacticus_Display_Message(message,verbosityWarn)
             ! Store the new maximum fractional error.
             fractionalErrorMaximum=fractionalError
          end if
          !$omp end critical (Very_Simple_Disk_Post_Evolve_Check)
          ! Get the total mass of the disk material
          diskMass= disk%massGas    () &
               &   +disk%massStellar()
          if (diskMass == 0.0d0) then
             call disk%        massStellarSet(                  0.0d0)
             call disk%  abundancesStellarSet(         zeroAbundances)
             call disk%luminositiesStellarSet(zeroStellarLuminosities)
          end if
          ! Reset the gas mass of the disk.
          call disk%      massGasSet(         0.0d0)
          call disk%abundancesGasSet(zeroAbundances)
          ! Record that state was changed.
          status=GSL_Failure
       end if
    end select
    return
  end subroutine Node_Component_Disk_Very_Simple_Post_Step

  subroutine Node_Component_Disk_Very_Simple_Create(node)
    !% Create properties in a very simple disk component.
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

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Disk_Very_Simple_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Disk_Very_Simple_Rate_Compute(node,odeConverged,interrupt,interruptProcedureReturn,propertyType)
    !% Compute the very simple disk node mass rate of change.
    use :: Abundances_Structure          , only : abundances
    use :: Galacticus_Nodes              , only : interruptTask       , nodeComponentDisk, nodeComponentDiskVerySimple, nodeComponentHotHalo, &
          &                                       propertyTypeInactive, treeNode         , defaultDiskComponent
    use :: Histories                     , only : history
    use :: Stellar_Luminosities_Structure, only : stellarLuminosities
    implicit none
    type            (treeNode            ), intent(inout), pointer :: node
    logical                               , intent(in   )          :: odeConverged
    class           (nodeComponentDisk   )               , pointer :: disk
    class           (nodeComponentHotHalo)               , pointer :: hotHalo
    logical                               , intent(inout)          :: interrupt
    procedure       (interruptTask       ), intent(inout), pointer :: interruptProcedureReturn
    procedure       (interruptTask       )               , pointer :: interruptProcedure
    integer                               , intent(in   )          :: propertyType
    double precision                                               :: stellarMassRate           , fuelMassRate         , &
         &                                                            massOutflowRate
    type            (history             )                         :: stellarHistoryRate
    type            (abundances          ), save                   :: fuelAbundancesRate        , stellarAbundancesRate, &
         &                                                            abundancesOutflowRate
    type            (stellarLuminosities ), save                   :: luminositiesStellarRates
    !$omp threadprivate(fuelAbundancesRate,stellarAbundancesRate,abundancesOutflowRate,luminositiesStellarRates)
    !$GLC attributes unused :: odeConverged

    ! Return immediately if inactive variables are requested.
    if (propertyType == propertyTypeInactive) return
    ! Get a local copy of the interrupt procedure.
    interruptProcedure => interruptProcedureReturn
    ! Return immediately if this class is not in use.
    if (.not.defaultDiskComponent%verySimpleIsActive()) return
    ! Get the disk and check that it is of our class.
    disk => node%disk()
    select type (disk)
    class is (nodeComponentDiskVerySimple)
       ! Check for a realistic disk, return immediately if disk is unphysical.
       if (disk%massGas() < 0.0d0) return
       ! Interrupt if the disk is not initialized.
       if (.not.disk%isInitialized()) then
          interrupt=.true.
          interruptProcedureReturn => Node_Component_Disk_Very_Simple_Create
          return
       end if
       ! Get rates.
       call Node_Component_Disk_Very_Simple_Rates(node,fuelMassRate,fuelAbundancesRate,stellarMassRate,stellarAbundancesRate,massOutflowRate,stellarHistoryRate,luminositiesStellarRates)
     if (massOutflowRate > 0.0d0) then
          ! Push to the hot halo.
          hotHalo => node%hotHalo      ()
          call hotHalo%outflowingMassRate(+massOutflowRate)
          call disk   %       massGasRate(-massOutflowRate)
          if (diskVerySimpleTrackAbundances) then
             abundancesOutflowRate=disk%abundancesGas()
             call abundancesOutflowRate%massToMassFraction(disk%massGas())
             abundancesOutflowRate=abundancesOutflowRate*massOutflowRate
             call hotHalo%outflowingAbundancesRate(+abundancesOutflowRate)
             call disk   %       abundancesGasRate(-abundancesOutflowRate)
          end if
       end if
    end select
    ! Return the procedure pointer.
    interruptProcedureReturn => interruptProcedure
    return
  end subroutine Node_Component_Disk_Very_Simple_Rate_Compute

  !# <analyticSolverTask>
  !#  <unitName>Node_Component_Disk_Very_Simple_Analytic_Solver</unitName>
  !# </analyticSolverTask>
  subroutine Node_Component_Disk_Very_Simple_Analytic_Solver(node,timeStart,timeEnd,solved)
    use :: Abundances_Structure          , only : abundances             , max                , operator(*)
    use :: Galacticus_Error              , only : Galacticus_Error_Report
    use :: Galacticus_Nodes              , only : nodeComponentBasic     , nodeComponentDisk  , nodeComponentHotHalo, nodeComponentSatellite, &
          &                                       treeNode
    use :: Histories                     , only : history
    use :: Stellar_Luminosities_Structure, only : max                    , stellarLuminosities
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

    if (diskVerySimpleUseAnalyticSolver) then
       disk => node%disk()
       if (node%isSatellite().and.disk%isInitialized()) then
          ! Luminosities can not be computed analytically.
          if (diskVerySimpleTrackLuminosities) call Galacticus_Error_Report('analytic solver does not support stellar luminosity calculation'//{introspection:location})
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
             if     (                                                                        &
                  &   massStellarAsymptotic     < diskVerySimpleAnalyticSolverPruneMassStars &
                  &  .and.                                                                   &
                  &   massGasInitial            < diskVerySimpleAnalyticSolverPruneMassGas   &
                  &  .and.                                                                   &
                  &   satellite%timeOfMerging() > timePresentDay                             &
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
                   if (diskVerySimpleTrackAbundances) then
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
                if (diskVerySimpleTrackAbundances) then
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
             call basic    %     timeSet(                      timeEnd )
             call satellite%mergeTimeSet(satellite%mergeTime()-timeStep)
          end if
          ! Record that we solved this system analytically.
          solved=.true.
       end if
    end if
    return
  end subroutine Node_Component_Disk_Very_Simple_Analytic_Solver

  subroutine Node_Component_Disk_Very_Simple_Rates(node,fuelMassRate,fuelAbundancesRate,stellarMassRate,stellarAbundancesRate,massOutflowRate,stellarHistoryRate,luminositiesStellarRates)
    !% Compute rates.
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
    double precision                                              :: diskDynamicalTime       , fuelMass              , &
         &                                                           energyInputRate         , starFormationRate
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
    ! Find rate of outflow of material from the disk and pipe it to the outflowed reservoir.
    massOutflowRate=starFormationFeedbackDisks_%outflowRate(node,energyInputRate,starFormationRate)
    if (massOutflowRate > 0.0d0) then
       ! Limit the outflow rate timescale to a multiple of the dynamical time.
       fuelMass         =disk                %massGas           (    )
       diskDynamicalTime=darkMatterHaloScale_%dynamicalTimescale(node)
       massOutflowRate  =min(massOutflowRate,fuelMass/diskOutflowTimescaleMinimum/diskDynamicalTime)
    end if
    return
  end subroutine Node_Component_Disk_Very_Simple_Rates

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Disk_Very_Simple_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Disk_Very_Simple_Scale_Set(node)
    !% Set scales for properties of {\normalfont \ttfamily node}.
    use :: Abundances_Structure          , only : abs              , abundances                 , max                , unitAbundances         , &
          &                                       zeroAbundances
    use :: Galacticus_Nodes              , only : nodeComponentDisk, nodeComponentDiskVerySimple, treeNode           , defaultDiskComponent
    use :: Histories                     , only : history
    use :: Stellar_Luminosities_Structure, only : abs              , max                        , stellarLuminosities, unitStellarLuminosities
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
       call disk%massGasScale    (max(mass,diskVerySimpleMassScaleAbsolute))
       call disk%massStellarScale(max(mass,diskVerySimpleMassScaleAbsolute))
       ! Set scale for gas and stellar abundances.
       abundancesTotal=disk%abundancesGas()+disk%abundancesStellar()
       call disk%abundancesGasScale    (max(abundancesTotal,unitAbundances*diskVerySimpleMassScaleAbsolute))
       call disk%abundancesStellarScale(max(abundancesTotal,unitAbundances*diskVerySimpleMassScaleAbsolute))
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
    !% Transfer any very simple disk associated with {\normalfont \ttfamily node} to its host halo.
    use :: Abundances_Structure                , only : zeroAbundances
    use :: Galacticus_Error                    , only : Galacticus_Error_Report
    use :: Galacticus_Nodes                    , only : nodeComponentDisk      , nodeComponentDiskVerySimple, nodeComponentSpheroid, treeNode
    use :: Satellite_Merging_Mass_Movements    , only : destinationMergerDisk  , destinationMergerSpheroid
    use :: Stellar_Luminosities_Structure      , only : zeroStellarLuminosities
    implicit none
    class  (*                    ), intent(inout) :: self
    type   (treeNode             ), intent(inout) :: node
    type   (treeNode             ), pointer       :: nodeHost
    class  (nodeComponentDisk    ), pointer       :: diskHost               , disk
    class  (nodeComponentSpheroid), pointer       :: spheroidHost
    integer                                       :: destinationGasSatellite, destinationGasHost       , &
         &                                           destinationStarsHost   , destinationStarsSatellite
    logical                                       :: mergerIsMajor
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
       select case (destinationGasSatellite)
       case (destinationMergerDisk)
          call diskHost    %massGasSet          (                              &
               &                                  diskHost    %      massGas() &
               &                                 +disk        %      massGas() &
               &                                )
          call diskHost    %abundancesGasSet    (                              &
               &                                  diskHost    %abundancesGas() &
               &                                 +disk        %abundancesGas() &
               &                                )
       case (destinationMergerSpheroid)
          call spheroidHost%massGasSet          (                              &
               &                                  spheroidHost%massGas      () &
               &                                 +disk        %massGas      () &
               &                                )
          call spheroidHost%abundancesGasSet    (                              &
               &                                  spheroidHost%abundancesGas() &
               &                                 +disk        %abundancesGas() &
               &                                )
       case default
          call Galacticus_Error_Report(                                    &
               &                       'unrecognized movesTo descriptor'// &
               &                       {introspection:location}            &
               &                      )
       end select
       call    disk%massGasSet                  (                                        &
            &                                                             0.0d0          &
            &                                   )
       call    disk%abundancesGasSet            (                                        &
            &                                                    zeroAbundances          &
            &                                   )
       ! Move the stellar component of the very simple disk to the host.
       select case (destinationStarsSatellite)
       case (destinationMergerDisk)
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
       case (destinationMergerSpheroid)
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
       case default
          call Galacticus_Error_Report(                                    &
               &                       'unrecognized movesTo descriptor'// &
               &                       {introspection:location}            &
               &                      )
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

end module Node_Component_Disk_Very_Simple
