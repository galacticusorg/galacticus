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

!% Contains a module that implements a very simple disk component.

module Node_Component_Disk_Very_Simple
  !% Implements a very simple disk component.
  use ISO_Varying_String
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Disk_Very_Simple_Post_Evolve  , Node_Component_Disk_Very_Simple_Rate_Compute         , &
       &    Node_Component_Disk_Very_Simple_Scale_Set    , Node_Component_Disk_Very_Simple_Satellite_Merging    , &
       &    Node_Component_Disk_Very_Simple_Initialize   , Node_Component_Disk_Very_Simple_Pre_Evolve           , &
       &    Node_Component_Disk_Very_Simple_Rates        , Node_Component_Disk_Very_Simple_Analytic_Solver

  !# <component>
  !#  <class>disk</class>
  !#  <name>verySimple</name>
  !#  <isDefault>no</isDefault>
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
  !#     <name>starFormationRate</name>
  !#     <attributes isSettable="false" isGettable="true" isEvolvable="false" isDeferred="get" />
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <isVirtual>true</isVirtual>
  !#     <output condition="[[diskOutputStarFormationRate]]" unitsInSI="massSolar/gigaYear" comment="Disk star formation rate."/>
  !#   </property>
  !#   <property>
  !#     <name>stellarPropertiesHistory</name>
  !#     <type>history</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#   </property>
  !#  </properties>
  !#  <bindings>
  !#   <binding method="attachPipe"   function="Node_Component_Disk_Very_Simple_Attach_Pipe" description="Attach pipes to the very simple disk component." bindsTo="component" returnType="\void" arguments="" />
  !#   <binding method="enclosedMass" function="Node_Component_Disk_Very_Simple_Enclosed_Mass" bindsTo="component" />
  !#  </bindings>
  !#  <functions>objects.nodes.components.disk.very_simple.bound_functions.inc</functions>
  !# </component>

  ! Record of whether this module has been initialized.
  logical          :: moduleInitialized                         =.false.

  ! Record of whether to use the simple disk analytic solver.
  logical          :: diskVerySimpleUseAnalyticSolver
  double precision :: diskVerySimpleAnalyticSolverPruneMassStars        , diskVerySimpleAnalyticSolverPruneMassGas, &
       &              timePresentDay
 
  ! Parameters controlling the physical implementation.
  double precision :: diskOutflowTimescaleMinimum                       , diskStarFormationTimescaleMinimum       , &
       &              diskVerySimpleMassScaleAbsolute                   , diskVerySimpleSurfaceDensityThreshold   , &
       &              diskVerySimpleSurfaceDensityVelocityExponent
  logical          :: diskVerySimpleTrackAbundances
  
contains

  !# <nodeComponentInitializationTask>
  !#  <unitName>Node_Component_Disk_Very_Simple_Initialize</unitName>
  !# </nodeComponentInitializationTask>
  subroutine Node_Component_Disk_Very_Simple_Initialize()
    !% Initializes the tree node very simple disk component module.
    use Input_Parameters
    use Cosmology_Functions
    implicit none
    type (nodeComponentDiskVerySimple)          :: diskVerySimpleComponent
    class(cosmologyFunctionsClass    ), pointer :: cosmologyFunctions_

    ! Initialize the module if necessary.
    !$omp critical (Node_Component_Disk_Very_Simple_Initialize)
    if (defaultDiskComponent%verySimpleIsActive().and..not.moduleInitialized) then
       ! Read parameters controlling the physical implementation.
       !@ <inputParameter>
       !@   <name>diskVerySimpleMassScaleAbsolute</name>
       !@   <defaultValue>$100 M_\odot$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The absolute mass scale below which calculations in the very simple disk component are allowed to become inaccurate.
       !@   </description>
       !@   <type>double</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('diskVerySimpleMassScaleAbsolute',diskVerySimpleMassScaleAbsolute,defaultValue=100.0d0)
       !@ <inputParameter>
       !@   <name>diskOutflowTimescaleMinimum</name>
       !@   <defaultValue>$10^{-3}$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The minimum timescale (in units of the halo dynamical time) on which outflows may deplete gas in the disk.
       !@   </description>
       !@   <type>double</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('diskOutflowTimescaleMinimum',diskOutflowTimescaleMinimum,defaultValue=1.0d-3)
       !@ <inputParameter>
       !@   <name>diskStarFormationTimescaleMinimum</name>
       !@   <defaultValue>$10^{-3}$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The minimum timescale (in units of the halo dynamical time) on which star formation may occur in the disk.
       !@   </description>
       !@   <type>double</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('diskStarFormationTimescaleMinimum',diskStarFormationTimescaleMinimum,defaultValue=1.0d-3)
       !@ <inputParameter>
       !@   <name>diskVerySimpleTrackAbundances</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Specifies whether or not to track abundances in the very simple disk component.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>0..1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('diskVerySimpleTrackAbundances',diskVerySimpleTrackAbundances,defaultValue=.false.)
       !@ <inputParameter>
       !@   <name>diskVerySimpleSurfaceDensityThreshold</name>
       !@   <defaultValue>0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The threshold gas surface denisty above this star formation occurs [$M_\odot$/Mpc$^2$].
       !@   </description>
       !@   <type>double</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('diskVerySimpleSurfaceDensityThreshold',diskVerySimpleSurfaceDensityThreshold,defaultValue=0.0d0)
       !@ <inputParameter>
       !@   <name>diskVerySimpleSurfaceDensityVelocityExponent</name>
       !@   <defaultValue>0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The exponent of velocity in the threshold gas surface denisty above this star formation occurs.
       !@   </description>
       !@   <type>double</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('diskVerySimpleSurfaceDensityVelocityExponent',diskVerySimpleSurfaceDensityVelocityExponent,defaultValue=0.0d0)
       !@ <inputParameter>
       !@   <name>diskVerySimpleUseAnalyticSolver</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    If true, employ an analytic ODE solver when evolving satellites.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('diskVerySimpleUseAnalyticSolver',diskVerySimpleUseAnalyticSolver,defaultValue=.false.)
       !@ <inputParameter>
       !@   <name>diskVerySimpleAnalyticSolverPruneMassGas</name>
       !@   <defaultValue>0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Gas mass below which the analytic solver will prune a galaxy.
       !@   </description>
       !@   <type>double</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('diskVerySimpleAnalyticSolverPruneMassGas',diskVerySimpleAnalyticSolverPruneMassGas,defaultValue=0.0d0)
       !@ <inputParameter>
       !@   <name>diskVerySimpleAnalyticSolverPruneMassStars</name>
       !@   <defaultValue>0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Stellar mass below which the analytic solver will prune a galaxy.
       !@   </description>
       !@   <type>double</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('diskVerySimpleAnalyticSolverPruneMassStars',diskVerySimpleAnalyticSolverPruneMassStars,defaultValue=0.0d0)
       ! If using the analytic solver, find the time at the present day.
       if (diskVerySimpleUseAnalyticSolver) then
          cosmologyFunctions_ => cosmologyFunctions            (     )
          timePresentDay      =  cosmologyFunctions_%cosmicTime(1.0d0)
       end if
       ! Attach the cooling mass pipe from the hot halo component.
       call diskVerySimpleComponent%attachPipe()
       ! Bind the star formation rate function.
       call diskVerySimpleComponent%starFormationRateFunction(Node_Component_Disk_Very_Simple_SFR)
       ! Record that the module is now initialized.
       moduleInitialized=.true.
    end if
    !$omp end critical (Node_Component_Disk_Very_Simple_Initialize)
    return
  end subroutine Node_Component_Disk_Very_Simple_Initialize

  !# <preEvolveTask>
  !# <unitName>Node_Component_Disk_Very_Simple_Pre_Evolve</unitName>
  !# </preEvolveTask>
  subroutine Node_Component_Disk_Very_Simple_Pre_Evolve(thisNode)
    !% Ensure the disk has been initialized.
    implicit none
    type (treeNode         ), intent(inout), pointer :: thisNode
    class(nodeComponentDisk)               , pointer :: thisDiskComponent

    ! Get the disk component.
    thisDiskComponent => thisNode%disk()
    ! Check if an exponential disk component exists.
    select type (thisDiskComponent)
       class is (nodeComponentDiskVerySimple)
          ! Initialize the disk
       call Node_Component_Disk_Very_Simple_Create(thisNode)
    end select
    return
  end subroutine Node_Component_Disk_Very_Simple_Pre_Evolve

  !# <postEvolveTask>
  !# <unitName>Node_Component_Disk_Very_Simple_Post_Evolve</unitName>
  !# </postEvolveTask>
  subroutine Node_Component_Disk_Very_Simple_Post_Evolve(thisNode)
    !% Catch rounding errors in the very simple disk gas evolution.
    use Galacticus_Display
    use String_Handling
    use Histories
    use Abundances_Structure
    implicit none
    type            (treeNode          ), intent(inout), pointer :: thisNode
    class           (nodeComponentDisk )               , pointer :: thisDiskComponent
    class           (nodeComponentBasic)               , pointer :: thisBasicComponent
    double precision                    , save                   :: fractionalErrorMaximum=0.0d0
    double precision                                             :: diskMass                    , fractionalError
    character       (len=20            )                         :: valueString
    type            (varying_string    )                         :: message
    type            (history           )                         :: stellarPropertiesHistory

    ! Get the disk component.
    thisDiskComponent => thisNode%disk()
    ! Check if a very simple disk component exists.
    select type (thisDiskComponent)
    class is (nodeComponentDiskVerySimple)
       ! Trim the stellar populations properties future history.
       thisBasicComponent => thisNode%basic()
       stellarPropertiesHistory=thisDiskComponent%stellarPropertiesHistory()
       call stellarPropertiesHistory%trim(thisBasicComponent%time())
       call thisDiskComponent%stellarPropertiesHistorySet(stellarPropertiesHistory)
       ! Trap negative gas masses.
       if (thisDiskComponent%massGas() < 0.0d0) then
          ! Check if this exceeds the maximum previously recorded error.
          fractionalError=   abs(thisDiskComponent%massGas    ()) &
               &          /(                                      &
               &                 thisDiskComponent%massStellar()  &
               &            +abs(thisDiskComponent%massGas    ()) &
               &           )
          !$omp critical (Very_Simple_Disk_Post_Evolve_Check)
          if (fractionalError > fractionalErrorMaximum) then
             ! Report a warning.
             message='Warning: disk has negative gas mass (fractional error exceeds any previously reported):'//char(10)
             message=message//'  Node index        = '//thisNode%index() //char(10)
             write (valueString,'(e12.6)') thisDiskComponent%massGas()
             message=message//'  Disk gas mass     = '//trim(valueString)//char(10)
             write (valueString,'(e12.6)') thisDiskComponent%massStellar()
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
          diskMass= thisDiskComponent%massGas    () &
               &   +thisDiskComponent%massStellar()
          if (diskMass == 0.0d0) then
             call thisDiskComponent%      massStellarSet(         0.0d0)
             call thisDiskComponent%abundancesStellarSet(zeroAbundances)
          end if
          ! Reset the gas mass of the disk.
          call thisDiskComponent%      massGasSet(         0.0d0)
          call thisDiskComponent%abundancesGasSet(zeroAbundances)
       end if
    end select
    return
  end subroutine Node_Component_Disk_Very_Simple_Post_Evolve

  subroutine Node_Component_Disk_Very_Simple_Create(thisNode)
    !% Create properties in a very simple disk component.
    use Histories
    use Stellar_Population_Properties
    implicit none
    type   (treeNode         ), intent(inout), pointer :: thisNode
    class  (nodeComponentDisk)               , pointer :: thisDiskComponent
    type   (history          )                         :: stellarPropertiesHistory
    logical                                            :: createStellarPropertiesHistory

    ! Get the disk component.
    thisDiskComponent => thisNode%disk()
    ! Exit if already initialized.
    if (thisDiskComponent%isInitialized()) return
    ! Determine which histories must be created.
    stellarPropertiesHistory      =thisDiskComponent%stellarPropertiesHistory        ()
    createStellarPropertiesHistory=.not.             stellarPropertiesHistory%exists ()
    call                                             stellarPropertiesHistory%destroy()
    ! Create the stellar properties history.
    if (createStellarPropertiesHistory) then
       ! Create the stellar properties history.
       call Stellar_Population_Properties_History_Create (thisNode,stellarPropertiesHistory)
       call thisDiskComponent%stellarPropertiesHistorySet(         stellarPropertiesHistory)
    end if
    ! Record that the disk has been initialized.
    call thisDiskComponent%isInitializedSet(.true.)
    return
  end subroutine Node_Component_Disk_Very_Simple_Create

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Disk_Very_Simple_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Disk_Very_Simple_Rate_Compute(node,interrupt,interruptProcedureReturn)
    !% Compute the very simple disk node mass rate of change.
    use Star_Formation_Feedback_Disks
    use Stellar_Feedback
    use Stellar_Population_Properties
    use Dark_Matter_Halo_Scales
    use Abundances_Structure
    use Galactic_Structure_Options
    use Histories
    use Stellar_Luminosities_Structure
    implicit none
    type            (treeNode                    ), intent(inout), pointer :: node
    class           (nodeComponentDisk           )               , pointer :: disk
    class           (nodeComponentHotHalo        )               , pointer :: hotHalo
    logical                                       , intent(inout)          :: interrupt
    procedure       (Interrupt_Procedure_Template), intent(inout), pointer :: interruptProcedureReturn
    procedure       (Interrupt_Procedure_Template)               , pointer :: interruptProcedure
    double precision                                                       :: stellarMassRate         , fuelMassRate         , &
         &                                                                    massOutflowRate
    type            (history                     )                         :: stellarHistoryRate
    type            (abundances                  )                         :: fuelAbundancesRate      , stellarAbundancesRate, &
         &                                                                    abundancesOutflowRate
    
    ! Get a local copy of the interrupt procedure.
    interruptProcedure => interruptProcedureReturn
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
       call Node_Component_Disk_Very_Simple_Rates(node,fuelMassRate,fuelAbundancesRate,stellarMassRate,stellarAbundancesRate,massOutflowRate,stellarHistoryRate)
       ! If not tracking abundances, zero abundance rates here.
       if (.not.diskVerySimpleTrackAbundances) then
          fuelAbundancesRate   =zeroAbundances
          stellarAbundancesRate=zeroAbundances
       end if
       ! Adjust rates.
       call                                  disk%             massStellarRate(      stellarMassRate)
       call                                  disk%                 massGasRate(         fuelMassRate)
       call                                  disk%       abundancesStellarRate(stellarAbundancesRate)
       call                                  disk%           abundancesGasRate(   fuelAbundancesRate)
       if (stellarHistoryRate%exists()) call disk%stellarPropertiesHistoryRate(   stellarHistoryRate)
       if (massOutflowRate > 0.0d0) then
          ! Push to the hot halo.
          hotHalo               => node%hotHalo      ()
          abundancesOutflowRate =  disk%abundancesGas()
          call abundancesOutflowRate%massToMassFraction(disk%massGas())
          abundancesOutflowRate =abundancesOutflowRate*massOutflowRate
          call hotHalo%      outflowingMassRate(+      massOutflowRate)
          call hotHalo%outflowingAbundancesRate(+abundancesOutflowRate)
          call disk   %             massGasRate(-      massOutflowRate)
          call disk   %       abundancesGasRate(-abundancesOutflowRate)
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
    use Histories
    use Abundances_Structure
    implicit none
    type            (treeNode              ), intent(inout), pointer   :: node
    double precision                        , intent(in   )            :: timeStart           , timeEnd
    logical                                 , intent(inout)            :: solved
    type            (treeNode              ),                pointer   :: hostNode
    class           (nodeComponentBasic    )               , pointer   :: basic               , hostBasic               , &
         &                                                                hostParentBasic
    class           (nodeComponentDisk     )               , pointer   :: disk
    class           (nodeComponentHotHalo  )               , pointer   :: hotHalo             , hostHotHalo
    class           (nodeComponentSatellite)               , pointer   :: satellite
    double precision                                       , parameter :: massTolerance=1.0d-6
    double precision                                                   :: massGasInitial      , massStellarInitial      , &
         &                                                                timescaleFuel       , timescaleOutflow        , &
         &                                                                timescaleStellar    , massGasFinal            , &
         &                                                                massStellarFinal    , massOutflowed           , &
         &                                                                rateFuel            , rateStars               , &
         &                                                                rateOutflow         , timeStep                , &
         &                                                                exponentialFactor   , massStellarAsymptotic
    type            (history               )                           :: stellarHistoryRate
    type            (abundances            )                           :: rateAbundanceFuel   , rateAbundanceStars      , &
         &                                                                abundancesGasInitial, abundancesStellarInitial, &
         &                                                                yieldMassEffective  , abundancesStellarFinal  , &
         &                                                                abundancesGasFinal  , abundancesOutflowed
    
    if (diskVerySimpleUseAnalyticSolver) then
       disk => node%disk()
       if (node%isSatellite().and.disk%isInitialized()) then
          ! Calculate analytic solution.
          timeStep          =timeEnd-timeStart
          massGasInitial          =disk%massGas          ()
          massStellarInitial      =disk%massStellar      ()
          abundancesGasInitial    =disk%abundancesGas    ()
          abundancesStellarInitial=disk%abundancesStellar()
          if (massGasInitial > massTolerance) then
             hotHalo   => node%hotHalo  ()
             satellite => node%satellite()
             call Node_Component_Disk_Very_Simple_Rates(node,rateFuel,rateAbundanceFuel,rateStars,rateAbundanceStars,rateOutflow,stellarHistoryRate)
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
                hostNode => node%parent
                do while (associated(hostNode%parent))
                   hostBasic       => hostNode%basic       ()
                   hostHotHalo     => hostNode%hotHalo     ()
                   hostParentBasic => hostNode%parent%basic()
                   massOutflowed   =  +massGasInitial                                                         &
                        &                  *(                                                                 &
                        &                    +timescaleFuel                                                   &
                        &                    /timescaleOutflow                                                &
                        &                  )                                                                  &
                        &                  *(                                                                 &
                        &                    +exp(-max(0.0d0,hostBasic      %time()-timeStart)/timescaleFuel) &
                        &                    -exp(-max(0.0d0,hostParentBasic%time()-timeStart)/timescaleFuel) &
                        &                  )
                   abundancesOutflowed=+(                               &
                        &                -(                             &
                        &                  +abundancesGasInitial        &
                        &                  +yieldMassEffective          &
                        &                  *(                           &
                        &                    +1.0d0                     &
                        &                    +(                         &
                        &                      +hostBasic        %time()&
                        &                      -timeStart               &
                        &                     )                         &
                        &                    /timescaleFuel             &
                        &                   )                           &
                        &                 )                             &
                        &                *exp(                          &
                        &                     -(                        &
                        &                       +hostBasic      %time() &
                        &                       -timeStart              &
                        &                      )                        &
                        &                     /timescaleFuel            &
                        &                    )                          &
                        &                +(                             &
                        &                  +abundancesGasInitial        &
                        &                  +yieldMassEffective          &
                        &                  *(                           &
                        &                    +1.0d0                     &
                        &                    +(                         &
                        &                      +hostParentBasic%time()  &
                        &                      -timeStart               &
                        &                     )                         &
                        &                    /timescaleFuel             &
                        &                   )                           &
                        &                 )                             &
                        &                *exp(                          &
                        &                     -(                        &
                        &                       +hostParentBasic%time() &
                        &                       -timeStart              &
                        &                      )                        &
                        &                     /timescaleFuel            &
                        &                    )                          &
                        &               )                               &
                        &              *timescaleFuel                   &
                        &              /timescaleOutflow                   
                   call hostHotHalo%outflowedMassSet      (                                   &
                        &                                  +hostHotHalo%outflowedMass      () &
                        &                                  +            massOutflowed         &
                        &                                 )
                   call hostHotHalo%outflowedAbundancesSet(                                   &
                        &                                  +hostHotHalo%outflowedAbundances() &
                        &                                  +            abundancesOutflowed   &
                        &                                 )                   
                   hostNode => hostNode%parent
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
  
  subroutine Node_Component_Disk_Very_Simple_Rates(node,fuelMassRate,fuelAbundancesRate,stellarMassRate,stellarAbundancesRate,massOutflowRate,stellarHistoryRate)
    !% Compute rates.
    use Star_Formation_Feedback_Disks
    use Stellar_Feedback
    use Stellar_Population_Properties
    use Dark_Matter_Halo_Scales
    use Abundances_Structure
    use Galactic_Structure_Options
    use Histories
    use Stellar_Luminosities_Structure
    implicit none
    type            (treeNode                ), intent(inout), pointer :: node
    type            (history                 ), intent(inout)          :: stellarHistoryRate
    double precision                          , intent(  out)          :: fuelMassRate            , stellarMassRate       , &
         &                                                                massOutflowRate
    type            (abundances              ), intent(  out)          :: fuelAbundancesRate      , stellarAbundancesRate
    class           (nodeComponentDisk       )               , pointer :: disk
    class           (darkMatterHaloScaleClass)               , pointer :: darkMatterHaloScale_
    double precision                                                   :: diskDynamicalTime       , fuelMass              , &
         &                                                                energyInputRate         , starFormationRate
    type            (stellarLuminosities     )                         :: luminositiesStellarRates
    type            (abundances              )                         :: fuelAbundances
    
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
       starFormationRate=Node_Component_Disk_Very_Simple_SFR(disk)
    end select
    ! Find rates of change of stellar mass, and gas mass.
    stellarHistoryRate=disk%stellarPropertiesHistory()
    call Stellar_Population_Properties_Rates(starFormationRate,fuelAbundances,componentTypeDisk,node,stellarHistoryRate &
         &,stellarMassRate,stellarAbundancesRate,luminositiesStellarRates,fuelMassRate,fuelAbundancesRate,energyInputRate)
    ! Find rate of outflow of material from the disk and pipe it to the outflowed reservoir.
    massOutflowRate=Star_Formation_Feedback_Disk_Outflow_Rate(node,starFormationRate,energyInputRate)
    if (massOutflowRate > 0.0d0) then
       ! Limit the outflow rate timescale to a multiple of the dynamical time.
       darkMatterHaloScale_ => darkMatterHaloScale()
       fuelMass             =  disk                %massGas           (    )
       diskDynamicalTime    =  darkMatterHaloScale_%dynamicalTimescale(node)
       massOutflowRate      =  min(massOutflowRate,fuelMass/diskOutflowTimescaleMinimum/diskDynamicalTime)
    end if
    return
  end subroutine Node_Component_Disk_Very_Simple_Rates
  
  !# <scaleSetTask>
  !#  <unitName>Node_Component_Disk_Very_Simple_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Disk_Very_Simple_Scale_Set(thisNode)
    !% Set scales for properties of {\normalfont \ttfamily thisNode}.
    use Abundances_Structure
    use Histories
    use Stellar_Population_Properties
    implicit none
    type            (treeNode         ), intent(inout), pointer :: thisNode
    class           (nodeComponentDisk)               , pointer :: thisDiskComponent
    double precision                                            :: mass
    type            (history          )                         :: stellarPopulationHistoryScales
    type            (abundances       )                         :: abundancesTotal

    ! Get the disk component.
    thisDiskComponent => thisNode%disk()
    ! Check if a very simple disk component exists.
    select type (thisDiskComponent)
    class is (nodeComponentDiskVerySimple)
       ! Set scale for gas and stellar mass.
       mass=thisDiskComponent%massGas()+thisDiskComponent%massStellar()
       call thisDiskComponent%massGasScale    (max(mass,diskVerySimpleMassScaleAbsolute))
       call thisDiskComponent%massStellarScale(max(mass,diskVerySimpleMassScaleAbsolute))
       ! Set scale for gas and stellar abundances.
       abundancesTotal=thisDiskComponent%abundancesGas()+thisDiskComponent%abundancesStellar()
       call thisDiskComponent%abundancesGasScale    (max(abundancesTotal,unitAbundances*diskVerySimpleMassScaleAbsolute))
       call thisDiskComponent%abundancesStellarScale(max(abundancesTotal,unitAbundances*diskVerySimpleMassScaleAbsolute))
       ! Set scales for stellar population properties and star formation histories.
       stellarPopulationHistoryScales=thisDiskComponent%stellarPropertiesHistory()
       call Stellar_Population_Properties_Scales           (stellarPopulationHistoryScales,thisDiskComponent%massStellar(),zeroAbundances)
       call thisDiskComponent%stellarPropertiesHistoryScale(stellarPopulationHistoryScales                                               )
       call stellarPopulationHistoryScales%destroy()
    end select
    return
  end subroutine Node_Component_Disk_Very_Simple_Scale_Set

  !# <satelliteMergerTask>
  !#  <unitName>Node_Component_Disk_Very_Simple_Satellite_Merging</unitName>
  !#  <after>Satellite_Merging_Mass_Movement_Store</after>
  !#  <after>Satellite_Merging_Remnant_Size</after>
  !# </satelliteMergerTask>
  subroutine Node_Component_Disk_Very_Simple_Satellite_Merging(thisNode)
    !% Transfer any very simple disk associated with {\normalfont \ttfamily thisNode} to its host halo.
    use Satellite_Merging_Mass_Movements_Descriptors
    use Galacticus_Error
    use Abundances_Structure
    implicit none
    type (treeNode             ), intent(inout), pointer :: thisNode
    type (treeNode             )               , pointer :: hostNode
    class(nodeComponentDisk    )               , pointer :: hostDiskComponent    , thisDiskComponent
    class(nodeComponentSpheroid)               , pointer :: hostSpheroidComponent
    
    ! Check that the disk is of the verySimple class.
    thisDiskComponent => thisNode%disk()
    select type (thisDiskComponent)
    class is (nodeComponentDiskVerySimple)

       ! Find the node to merge with and its disk component (and spheroid if necessary).
       hostNode              => thisNode%mergesWith(                 )
       hostDiskComponent     => hostNode%disk      (autoCreate=.true.)
       if     (                                          &
            &   thisMergerGasMovesTo  == movesToSpheroid &
            &  .or.                                      &
            &   thisMergerStarsMoveTo == movesToSpheroid &
            & )                                          &
            & hostSpheroidComponent => hostNode%spheroid(autoCreate=.true.)
       ! Move the gas component of the very simple disk to the host.
       select case (thisMergerGasMovesTo)
       case (movesToDisk)
          call hostDiskComponent    %massGasSet          (                                       &
               &                                           hostDiskComponent    %      massGas() &
               &                                          +thisDiskComponent    %      massGas() &
               &                                         )
          call hostDiskComponent    %abundancesGasSet    (                                       &
               &                                           hostDiskComponent    %abundancesGas() &
               &                                          +thisDiskComponent    %abundancesGas() &
               &                                         )
       case (movesToSpheroid)
          call hostSpheroidComponent%massGasSet          (                                       &
               &                                           hostSpheroidComponent%massGas      () &
               &                                          +thisDiskComponent    %massGas      () &
               &                                         )
          call hostSpheroidComponent%abundancesGasSet    (                                       &
               &                                           hostSpheroidComponent%abundancesGas() &
               &                                          +thisDiskComponent    %abundancesGas() &
               &                                         )
       case default
          call Galacticus_Error_Report(                                                      &
               &                       'Node_Component_Disk_Very_Simple_Satellite_Merging',  &
               &                       'unrecognized movesTo descriptor'                     &
               &                      )
       end select
       call    thisDiskComponent%massGasSet          (                                       &
            &                                                                  0.0d0         &
            &                                        )
       call    thisDiskComponent%abundancesGasSet    (                                       &
            &                                                         zeroAbundances         &
            &                                        )

       ! Move the stellar component of the very simple disk to the host.
       select case (thisMergerStarsMoveTo)
       case (movesToDisk)
          call hostDiskComponent    %massStellarSet      (                                           &
               &                                           hostDiskComponent    %      massStellar() &
               &                                          +thisDiskComponent    %      massStellar() &
               &                                         )
          call hostDiskComponent    %abundancesStellarSet(                                           &
               &                                           hostDiskComponent    %abundancesStellar() &
               &                                          +thisDiskComponent    %abundancesStellar() &
               &                                         )
       case (movesToSpheroid)
          call hostSpheroidComponent%massStellarSet      (                                           &
               &                                           hostSpheroidComponent%massStellar      () &
               &                                          +thisDiskComponent    %massStellar      () &
               &                                         )
          call hostSpheroidComponent%abundancesStellarSet(                                           &
               &                                           hostSpheroidComponent%abundancesStellar() &
               &                                          +thisDiskComponent    %abundancesStellar() &
               &                                         )
       case default
          call Galacticus_Error_Report(                                                      &
               &                       'Node_Component_Disk_Very_Simple_Satellite_Merging',  &
               &                       'unrecognized movesTo descriptor'                     &
               &                      )
       end select
       call    thisDiskComponent%massStellarSet      (                                       &
            &                                                               0.0d0            &
            &                                        )
       call    thisDiskComponent%abundancesStellarSet(                                       &
            &                                                      zeroAbundances            &
            &                                        )
    end select
    return
  end subroutine Node_Component_Disk_Very_Simple_Satellite_Merging

  double precision function Node_Component_Disk_Very_Simple_SFR(self)
    !% Return the star formation rate of the very simple disk.
    use Star_Formation_Timescales_Disks
    use Dark_Matter_Halo_Scales
    use Dark_Matter_Profiles
    use Numerical_Constants_Math
    implicit none
    class           (nodeComponentDiskVerySimple), intent(inout) :: self
    class           (darkMatterHaloScaleClass   ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileClass     ), pointer       :: darkMatterProfile_
    double precision                             , parameter     :: velocityNormalization  =200.0d0
    double precision                                             :: diskDynamicalTime              , gasMass              , &
         &                                                          starFormationTimescale         , surfaceDensityCentral, &
         &                                                          radiusThreshold                , massStarForming      , &
         &                                                          surfaceDensityThreshold

    ! Get the gas mass.
    gasMass=self%massGas()    
    ! Get the star formation timescale.
    starFormationTimescale=Star_Formation_Timescale_Disk(self%hostNode)
    ! Limit the star formation timescale to a multiple of the dynamical time.
    darkMatterHaloScale_   => darkMatterHaloScale()
    darkMatterProfile_     => darkMatterProfile  ()
    diskDynamicalTime      =darkMatterHaloScale_%dynamicalTimescale(self%hostNode)
    starFormationTimescale =max(starFormationTimescale,diskStarFormationTimescaleMinimum*diskDynamicalTime)
    ! If timescale is finite and gas mass is positive, then compute star formation rate.
    if (starFormationTimescale > 0.0d0 .and. gasMass > 0.0d0 .and. self%radius() > 0.0d0) then
       ! Find mass of gas actively involved in star formation.
       if (diskVerySimpleSurfaceDensityThreshold > 0.0d0) then
          surfaceDensityCentral  =+gasMass          &
               &                  /self%radius()**2 &
               &                  /2.0d0            &
               &                  /Pi
          surfaceDensityThreshold=+  diskVerySimpleSurfaceDensityThreshold                      &
               &                  *(                                                            &
               &                    +darkMatterProfile_ %circularVelocityMaximum(self%hostNode) &
               &                    /velocityNormalization                                      &
               &                  )**diskVerySimpleSurfaceDensityVelocityExponent
          if (surfaceDensityCentral > surfaceDensityThreshold) then
             radiusThreshold=-log(surfaceDensityThreshold/surfaceDensityCentral)
             massStarForming=gasMass*(1.0d0-(1.0d0+radiusThreshold)*exp(-radiusThreshold))
          else
             massStarForming=0.0d0
          end if
       else
          massStarForming=gasMass
       end if
       Node_Component_Disk_Very_Simple_SFR=massStarForming/starFormationTimescale
    else
       Node_Component_Disk_Very_Simple_SFR=0.0d0
    end if
    return
  end function Node_Component_Disk_Very_Simple_SFR

end module Node_Component_Disk_Very_Simple
