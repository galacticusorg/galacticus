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

!% Contains a module that implements a very simple spheroid component.

module Node_Component_Spheroid_Very_Simple
  !% Implements a very simple spheroid component.
  use ISO_Varying_String
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Spheroid_Very_Simple_Post_Evolve               , Node_Component_Spheroid_Very_Simple_Rate_Compute         , &
       &    Node_Component_Spheroid_Very_Simple_Scale_Set                 , Node_Component_Spheroid_Very_Simple_Satellite_Merging    , &
       &    Node_Component_Spheroid_Very_Simple_Initialize                , Node_Component_Spheroid_Very_Simple_Pre_Evolve           , &
       &    Node_Component_Spheroid_Very_Simple_Rates                     , Node_Component_Spheroid_Very_Simple_Radius_Solver        , &
       &    Node_Component_Spheroid_Very_Simple_Radius_Solver_Plausibility
  
  !# <component>
  !#  <class>spheroid</class>
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
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of stars in the very simple spheroid."/>
  !#   </property>
  !#   <property>
  !#     <name>abundancesStellar</name>
  !#     <type>abundances</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of metals in the stellar phase of the very simple spheroid."/>
  !#   </property>
  !#   <property>
  !#     <name>massGas</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of gas in the very simple spheroid."/>
  !#   </property>
  !#   <property>
  !#     <name>abundancesGas</name>
  !#     <type>abundances</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true"  />
  !#     <output unitsInSI="massSolar" comment="Mass of metals in the gas phase of the very simple spheroid."/>
  !#   </property>
  !#   <property>
  !#     <name>starFormationRate</name>
  !#     <attributes isSettable="false" isGettable="true" isEvolvable="false" isDeferred="get" createIfNeeded="true" />
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <isVirtual>true</isVirtual>
  !#     <output condition="[[spheroidOutputStarFormationRate]]" unitsInSI="massSolar/gigaYear" comment="Spheroid star formation rate."/>
  !#   </property>
  !#   <property>
  !#     <name>stellarPropertiesHistory</name>
  !#     <type>history</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#   </property>
  !#   <property>
  !#     <name>radius</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <output unitsInSI="megaparsec" comment="Radial scale length in the spheroid."/>
  !#   </property>
  !#   <property>
  !#     <name>halfMassRadius</name>
  !#     <attributes isSettable="false" isGettable="true" isEvolvable="false" />
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <isVirtual>true</isVirtual>
  !#     <getFunction>Node_Component_Spheroid_Very_Simple_Half_Mass_Radius</getFunction>
  !#   </property>
  !#   <property>
  !#     <name>velocity</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <output unitsInSI="kilo" comment="Circular velocity of the spheroid."/>
  !#   </property>
  !#  </properties>
  !#  <bindings>
  !#   <binding method="enclosedMass" function="Node_Component_Spheroid_Very_Simple_Enclosed_Mass" bindsTo="component" />
  !#  </bindings>
  !#  <functions>objects.nodes.components.spheroid.very_simple.bound_functions.inc</functions>
  !# </component>

  ! Record of whether this module has been initialized.
  logical          :: moduleInitialized                         =.false.
 
  ! Parameters controlling the physical implementation.
  double precision :: spheroidOutflowTimescaleMinimum                       , spheroidStarFormationTimescaleMinimum       , &
       &              spheroidVerySimpleMassScaleAbsolute                   , spheroidMassToleranceAbsolute
  logical          :: spheroidVerySimpleTrackAbundances
  
contains

  !# <nodeComponentInitializationTask>
  !#  <unitName>Node_Component_Spheroid_Very_Simple_Initialize</unitName>
  !# </nodeComponentInitializationTask>
  subroutine Node_Component_Spheroid_Very_Simple_Initialize()
    !% Initializes the tree node very simple spheroid component module.
    use Input_Parameters
    use Cosmology_Functions
    implicit none
    type (nodeComponentSpheroidVerySimple)          :: spheroidVerySimpleComponent
    class(cosmologyFunctionsClass        ), pointer :: cosmologyFunctions_

    ! Initialize the module if necessary.
    !$omp critical (Node_Component_Spheroid_Very_Simple_Initialize)
    if (defaultSpheroidComponent%verySimpleIsActive().and..not.moduleInitialized) then
       ! Read parameters controlling the physical implementation.
       !@ <inputParameter>
       !@   <name>spheroidVerySimpleMassScaleAbsolute</name>
       !@   <defaultValue>$100 M_\odot$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The absolute mass scale below which calculations in the very simple spheroid component are allowed to become inaccurate.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('spheroidVerySimpleMassScaleAbsolute',spheroidVerySimpleMassScaleAbsolute,defaultValue=100.0d0)
       !@ <inputParameter>
       !@   <name>spheroidOutflowTimescaleMinimum</name>
       !@   <defaultValue>$10^{-3}$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The minimum timescale (in units of the halo dynamical time) on which outflows may deplete gas in the spheroid.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('spheroidOutflowTimescaleMinimum',spheroidOutflowTimescaleMinimum,defaultValue=1.0d-3)
       !@ <inputParameter>
       !@   <name>spheroidStarFormationTimescaleMinimum</name>
       !@   <defaultValue>$10^{-3}$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The minimum timescale (in units of the halo dynamical time) on which star formation may occur in the spheroid.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('spheroidStarFormationTimescaleMinimum',spheroidStarFormationTimescaleMinimum,defaultValue=1.0d-3)
       !@ <inputParameter>
       !@   <name>spheroidVerySimpleTrackAbundances</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Specifies whether or not to track abundances in the very simple spheroid component.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>0..1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('spheroidVerySimpleTrackAbundances',spheroidVerySimpleTrackAbundances,defaultValue=.false.)
       !@ <inputParameter>
       !@   <name>spheroidMassToleranceAbsolute</name>
       !@   <defaultValue>$10^{-6} M_\odot$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The mass tolerance used to judge whether the spheroid is physically plausible.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('spheroidMassToleranceAbsolute',spheroidMassToleranceAbsolute,defaultValue=1.0d-6)
       ! Bind the star formation rate function.
       call spheroidVerySimpleComponent%starFormationRateFunction(Node_Component_Spheroid_Very_Simple_SFR)
       ! Record that the module is now initialized.
       moduleInitialized=.true.
    end if
    !$omp end critical (Node_Component_Spheroid_Very_Simple_Initialize)
    return
  end subroutine Node_Component_Spheroid_Very_Simple_Initialize

  !# <preEvolveTask>
  !# <unitName>Node_Component_Spheroid_Very_Simple_Pre_Evolve</unitName>
  !# </preEvolveTask>
  subroutine Node_Component_Spheroid_Very_Simple_Pre_Evolve(thisNode)
    !% Ensure the spheroid has been initialized.
    implicit none
    type (treeNode             ), intent(inout), pointer :: thisNode
    class(nodeComponentSpheroid)               , pointer :: thisSpheroidComponent

    ! Get the spheroid component.
    thisSpheroidComponent => thisNode%spheroid()
    ! Check if an exponential spheroid component exists.
    select type (thisSpheroidComponent)
       class is (nodeComponentSpheroidVerySimple)
          ! Initialize the spheroid
          call Node_Component_Spheroid_Very_Simple_Create(thisNode)
    end select
    return
  end subroutine Node_Component_Spheroid_Very_Simple_Pre_Evolve

  !# <postEvolveTask>
  !# <unitName>Node_Component_Spheroid_Very_Simple_Post_Evolve</unitName>
  !# </postEvolveTask>
  subroutine Node_Component_Spheroid_Very_Simple_Post_Evolve(thisNode)
    !% Catch rounding errors in the very simple spheroid gas evolution.
    use Galacticus_Display
    use String_Handling
    use Histories
    use Abundances_Structure
    implicit none
    type            (treeNode              ), intent(inout), pointer :: thisNode
    class           (nodeComponentSpheroid )               , pointer :: thisSpheroidComponent
    class           (nodeComponentBasic    )               , pointer :: thisBasicComponent
    double precision                        , save                   :: fractionalErrorMaximum  =0.0d0
    double precision                                                 :: spheroidMass                  , fractionalError
    character       (len=20                )                         :: valueString
    type            (varying_string        )                         :: message
    type            (history               )                         :: stellarPropertiesHistory

    ! Get the spheroid component.
    thisSpheroidComponent => thisNode%spheroid()
    ! Check if a very simple spheroid component exists.
    select type (thisSpheroidComponent)
    class is (nodeComponentSpheroidVerySimple)
       ! Trim the stellar populations properties future history.
       thisBasicComponent => thisNode%basic()
       stellarPropertiesHistory=thisSpheroidComponent%stellarPropertiesHistory()
       call stellarPropertiesHistory%trim(thisBasicComponent%time())
       call thisSpheroidComponent%stellarPropertiesHistorySet(stellarPropertiesHistory)
       ! Trap negative gas masses.
       if (thisSpheroidComponent%massGas() < 0.0d0) then
          ! Check if this exceeds the maximum previously recorded error.
          fractionalError=   abs(thisSpheroidComponent%massGas    ()) &
               &          /(                                          &
               &                 thisSpheroidComponent%massStellar()  &
               &            +abs(thisSpheroidComponent%massGas    ()) &
               &           )
          !$omp critical (Very_Simple_Spheroid_Post_Evolve_Check)
          if (fractionalError > fractionalErrorMaximum) then
             ! Report a warning.
             message='Warning: spheroid has negative gas mass (fractional error exceeds any previously reported):'//char(10)
             message=message//'  Node index        = '//thisNode%index() //char(10)
             write (valueString,'(e12.6)') thisSpheroidComponent%massGas()
             message=message//'  Spheroid gas mass     = '//trim(valueString)//char(10)
             write (valueString,'(e12.6)') thisSpheroidComponent%massStellar()
             message=message//'  Spheroid stellar mass = '//trim(valueString)//char(10)
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
          !$omp end critical (Very_Simple_Spheroid_Post_Evolve_Check)
          ! Get the total mass of the spheroid material
          spheroidMass=+thisSpheroidComponent%massGas    () &
               &       +thisSpheroidComponent%massStellar()
          if (spheroidMass == 0.0d0) then
             call thisSpheroidComponent%      massStellarSet(         0.0d0)
             call thisSpheroidComponent%abundancesStellarSet(zeroAbundances)
          end if
          ! Reset the gas mass of the spheroid.
          call thisSpheroidComponent%      massGasSet(         0.0d0)
          call thisSpheroidComponent%abundancesGasSet(zeroAbundances)
       end if
    end select
    return
  end subroutine Node_Component_Spheroid_Very_Simple_Post_Evolve

  subroutine Node_Component_Spheroid_Very_Simple_Create(thisNode)
    !% Create properties in a very simple spheroid component.
    use Histories
    use Stellar_Population_Properties
    implicit none
    type   (treeNode             ), intent(inout), pointer :: thisNode
    class  (nodeComponentSpheroid)               , pointer :: thisSpheroidComponent
    type   (history              )                         :: stellarPropertiesHistory
    logical                                                :: createStellarPropertiesHistory

    ! Get the spheroid component.
    thisSpheroidComponent => thisNode%spheroid()
    ! Exit if already initialized.
    if (thisSpheroidComponent%isInitialized()) return
    ! Determine which histories must be created.
    stellarPropertiesHistory      =thisSpheroidComponent%stellarPropertiesHistory        ()
    createStellarPropertiesHistory=.not.                 stellarPropertiesHistory%exists ()
    call                                                 stellarPropertiesHistory%destroy()
    ! Create the stellar properties history.
    if (createStellarPropertiesHistory) then
       ! Create the stellar properties history.
       call Stellar_Population_Properties_History_Create     (thisNode,stellarPropertiesHistory)
       call thisSpheroidComponent%stellarPropertiesHistorySet(         stellarPropertiesHistory)
    end if
    ! Record that the spheroid has been initialized.
    call thisSpheroidComponent%isInitializedSet(.true.)
    return
  end subroutine Node_Component_Spheroid_Very_Simple_Create

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Spheroid_Very_Simple_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Spheroid_Very_Simple_Rate_Compute(node,interrupt,interruptProcedureReturn)
    !% Compute the very simple spheroid node mass rate of change.
    use Star_Formation_Feedback_Spheroids
    use Stellar_Feedback
    use Stellar_Population_Properties
    use Dark_Matter_Halo_Scales
    use Abundances_Structure
    use Galactic_Structure_Options
    use Histories
    use Stellar_Luminosities_Structure
    implicit none
    type            (treeNode                    ), intent(inout), pointer :: node
    class           (nodeComponentSpheroid       )               , pointer :: spheroid
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
    ! Get the spheroid and check that it is of our class.
    spheroid => node%spheroid()
    select type (spheroid)
    class is (nodeComponentSpheroidVerySimple)
       ! Check for a realistic spheroid, return immediately if spheroid is unphysical.
       if (spheroid%massGas() < 0.0d0) return
       ! Interrupt if the spheroid is not initialized.
       if (.not.spheroid%isInitialized()) then
          interrupt=.true.
          interruptProcedureReturn => Node_Component_Spheroid_Very_Simple_Create
          return
       end if
       ! Get rates.
       call Node_Component_Spheroid_Very_Simple_Rates(node,fuelMassRate,fuelAbundancesRate,stellarMassRate,stellarAbundancesRate,massOutflowRate,stellarHistoryRate)
       ! If not tracking abundances, zero abundance rates here.
       if (.not.spheroidVerySimpleTrackAbundances) then
          fuelAbundancesRate   =zeroAbundances
          stellarAbundancesRate=zeroAbundances
       end if
       ! Adjust rates.
       call                                  spheroid%             massStellarRate(      stellarMassRate)
       call                                  spheroid%                 massGasRate(         fuelMassRate)
       call                                  spheroid%       abundancesStellarRate(stellarAbundancesRate)
       call                                  spheroid%           abundancesGasRate(   fuelAbundancesRate)
       if (stellarHistoryRate%exists()) call spheroid%stellarPropertiesHistoryRate(   stellarHistoryRate)
       if (massOutflowRate > 0.0d0) then
          ! Push to the hot halo.
          hotHalo               => node%hotHalo      ()
          abundancesOutflowRate =  spheroid%abundancesGas()
          call abundancesOutflowRate%massToMassFraction(spheroid%massGas())
          abundancesOutflowRate =abundancesOutflowRate*massOutflowRate
          call hotHalo %      outflowingMassRate(+      massOutflowRate)
          call hotHalo %outflowingAbundancesRate(+abundancesOutflowRate)
          call spheroid%             massGasRate(-      massOutflowRate)
          call spheroid%       abundancesGasRate(-abundancesOutflowRate)
       end if
    end select
    ! Return the procedure pointer.
    interruptProcedureReturn => interruptProcedure
    return
  end subroutine Node_Component_Spheroid_Very_Simple_Rate_Compute

  subroutine Node_Component_Spheroid_Very_Simple_Rates(node,fuelMassRate,fuelAbundancesRate,stellarMassRate,stellarAbundancesRate,massOutflowRate,stellarHistoryRate)
    !% Compute rates.
    use Star_Formation_Feedback_Spheroids
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
    class           (nodeComponentSpheroid   )               , pointer :: spheroid
    class           (darkMatterHaloScaleClass)               , pointer :: darkMatterHaloScale_
    double precision                                                   :: spheroidDynamicalTime   , fuelMass              , &
         &                                                                energyInputRate         , starFormationRate
    type            (stellarLuminosities     )                         :: luminositiesStellarRates
    type            (abundances              )                         :: fuelAbundances
    
    ! Get the spheroid.
    spheroid => node%spheroid()
    ! Initialize to zero rates.
    fuelMassRate         =0.0d0
    stellarMassRate      =0.0d0
    massOutflowRate      =0.0d0
    fuelAbundancesRate   =zeroAbundances
    stellarAbundancesRate=zeroAbundances
    ! Check for a realistic spheroid, return immediately if spheroid is unphysical.
    if (spheroid%massGas() <= 0.0d0) return
    ! Compute fuel abundances.
    fuelAbundances=spheroid%abundancesGas()
    call fuelAbundances%massToMassFraction(spheroid%massGas())    
    ! Compute the star formation rate.
    select type (spheroid)
    class is (nodeComponentSpheroidVerySimple)
       starFormationRate=Node_Component_Spheroid_Very_Simple_SFR(spheroid)
    end select
    ! Find rates of change of stellar mass, and gas mass.
    stellarHistoryRate=spheroid%stellarPropertiesHistory()
    call Stellar_Population_Properties_Rates(starFormationRate,fuelAbundances,componentTypeSpheroid,node,stellarHistoryRate &
         &,stellarMassRate,stellarAbundancesRate,luminositiesStellarRates,fuelMassRate,fuelAbundancesRate,energyInputRate)
    ! Find rate of outflow of material from the spheroid and pipe it to the outflowed reservoir.
    massOutflowRate=Star_Formation_Feedback_Spheroid_Outflow_Rate(node,starFormationRate,energyInputRate)
    if (massOutflowRate > 0.0d0) then
       ! Limit the outflow rate timescale to a multiple of the dynamical time.
       darkMatterHaloScale_ => darkMatterHaloScale                    (    )
       fuelMass             =  spheroid            %massGas           (    )
       spheroidDynamicalTime=  darkMatterHaloScale_%dynamicalTimescale(node)
       massOutflowRate      =  min(massOutflowRate,fuelMass/spheroidOutflowTimescaleMinimum/spheroidDynamicalTime)
    end if
    return
  end subroutine Node_Component_Spheroid_Very_Simple_Rates
  
  !# <scaleSetTask>
  !#  <unitName>Node_Component_Spheroid_Very_Simple_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Spheroid_Very_Simple_Scale_Set(thisNode)
    !% Set scales for properties of {\normalfont \ttfamily thisNode}.
    use Abundances_Structure
    use Histories
    use Stellar_Population_Properties
    implicit none
    type            (treeNode             ), intent(inout), pointer :: thisNode
    class           (nodeComponentSpheroid)               , pointer :: thisSpheroidComponent
    double precision                                                :: mass
    type            (history              )                         :: stellarPopulationHistoryScales
    type            (abundances           )                         :: abundancesTotal

    ! Get the spheroid component.
    thisSpheroidComponent => thisNode%spheroid()
    ! Check if a very simple spheroid component exists.
    select type (thisSpheroidComponent)
    class is (nodeComponentSpheroidVerySimple)
       ! Set scale for gas and stellar mass.
       mass=thisSpheroidComponent%massGas()+thisSpheroidComponent%massStellar()
       call thisSpheroidComponent%massGasScale    (max(mass,spheroidVerySimpleMassScaleAbsolute))
       call thisSpheroidComponent%massStellarScale(max(mass,spheroidVerySimpleMassScaleAbsolute))
       ! Set scale for gas and stellar abundances.
       abundancesTotal=thisSpheroidComponent%abundancesGas()+thisSpheroidComponent%abundancesStellar()
       call thisSpheroidComponent%abundancesGasScale    (max(abundancesTotal,unitAbundances*spheroidVerySimpleMassScaleAbsolute))
       call thisSpheroidComponent%abundancesStellarScale(max(abundancesTotal,unitAbundances*spheroidVerySimpleMassScaleAbsolute))
       ! Set scales for stellar population properties and star formation histories.
       stellarPopulationHistoryScales=thisSpheroidComponent%stellarPropertiesHistory()
       call Stellar_Population_Properties_Scales               (stellarPopulationHistoryScales,thisSpheroidComponent%massStellar(),zeroAbundances)
       call thisSpheroidComponent%stellarPropertiesHistoryScale(stellarPopulationHistoryScales                                                   )
       call stellarPopulationHistoryScales%destroy()
    end select
    return
  end subroutine Node_Component_Spheroid_Very_Simple_Scale_Set

  !# <satelliteMergerTask>
  !#  <unitName>Node_Component_Spheroid_Very_Simple_Satellite_Merging</unitName>
  !#  <after>Satellite_Merging_Mass_Movement_Store</after>
  !#  <after>Satellite_Merging_Remnant_Size</after>
  !# </satelliteMergerTask>
  subroutine Node_Component_Spheroid_Very_Simple_Satellite_Merging(thisNode)
    !% Transfer any very simple spheroid associated with {\normalfont \ttfamily thisNode} to its host halo.
    use Satellite_Merging_Mass_Movements_Descriptors
    use Galacticus_Error
    use Abundances_Structure
    use Histories
    implicit none
    type (treeNode             ), intent(inout), pointer :: thisNode
    type (treeNode             )               , pointer :: hostNode
    class(nodeComponentSpheroid)               , pointer :: hostSpheroidComponent, thisSpheroidComponent
    class(nodeComponentDisk    )               , pointer :: hostDiskComponent
    type (history              )                         :: historyDisk          , historySpheroid      , &
         &                                                  thisHistory          , hostHistory

    ! Check that the very simple spheroid is active.
    if (defaultSpheroidComponent%verySimpleIsActive()) then
       
       ! Get the spheroid component, creating it if need be.
       thisSpheroidComponent => thisNode%spheroid(autoCreate=.true.)
       select type (thisSpheroidComponent)
       class is (nodeComponentSpheroidVerySimple)
          
          ! Find the node to merge with and its spheroid component.
          hostNode              => thisNode%mergesWith(                 )
          hostSpheroidComponent => hostNode%spheroid  (autoCreate=.true.)
          hostDiskComponent     => hostNode%disk      (autoCreate=.true.)

          ! Move gas material within the host if necessary.
          select case (thisHostGasMovesTo)
          case (movesToDisk)
             call hostDiskComponent    %          massGasSet(                                       &
                  &                                          +hostDiskComponent    %massGas      () &
                  &                                          +hostSpheroidComponent%massGas      () &
                  &                                         )
             call hostDiskComponent    %    abundancesGasSet(                                       &
                  &                                          +hostDiskComponent    %abundancesGas() &
                  &                                          +hostSpheroidComponent%abundancesGas() &
                  &                                         )
             call hostSpheroidComponent%          massGasSet(                                       &
                  &                                           0.0d0                                 &
                  &                                         )
             call hostSpheroidComponent%    abundancesGasSet(                                       &
                  &                                           zeroAbundances                        &
                  &                                         )
          case (movesToSpheroid)
             call hostSpheroidComponent%          massGasSet(                                       &
                  &                                          +hostSpheroidComponent%massGas      () &
                  &                                          +hostDiskComponent    %massGas      () &
                  &                                         )
             call hostSpheroidComponent%    abundancesGasSet(                                       &
                  &                                          +hostSpheroidComponent%abundancesGas() &
                  &                                          +hostDiskComponent    %abundancesGas() &
                  &                                         )
             call hostDiskComponent    %          massGasSet(                                       &
                  &                                           0.0d0                                 &
                  &                                         )
             call hostDiskComponent    %    abundancesGasSet(                                       &
                  &                                           zeroAbundances                        &
                  &                                         )
          case (doesNotMove)
             ! Do nothing.
          case default
             call Galacticus_Error_Report(                                                         &
                  &                       'Node_Component_Spheroid_Very_Simple_Satellite_Merging', &
                  &                       'unrecognized movesTo descriptor'                        &
                  &                      )
          end select
          
          ! Move stellar material within the host if necessary.
          select case (thisHostStarsMoveTo)
          case (movesToDisk)
             call hostDiskComponent    %      massStellarSet(                                           &
                  &                                          +hostDiskComponent    %      massStellar() &
                  &                                          +hostSpheroidComponent%      massStellar() &
                  &                                         )
             call hostDiskComponent    %abundancesStellarSet(                                           &
                  &                                           hostDiskComponent    %abundancesStellar() &
                  &                                          +hostSpheroidComponent%abundancesStellar() &
                  &                                         )
             call hostSpheroidComponent%      massStellarSet(                                           &
                  &                                           0.0d0                                     &
                  &                                         )
             call hostSpheroidComponent%abundancesStellarSet(                                           &
                  &                                           zeroAbundances                            &
                  &                                         )
             ! Also add stellar properties histories.
             historyDisk    =    hostDiskComponent%stellarPropertiesHistory()
             historySpheroid=hostSpheroidComponent%stellarPropertiesHistory()
             call historyDisk          %increment                  (historySpheroid)
             call historySpheroid      %reset                      (               )
             call hostDiskComponent    %stellarPropertiesHistorySet(historyDisk    )
             call hostSpheroidComponent%stellarPropertiesHistorySet(historySpheroid)
          case (movesToSpheroid)
             call hostSpheroidComponent%      massStellarSet  (                                         &
                  &                                           hostSpheroidComponent%      massStellar() &
                  &                                          +hostDiskComponent    %      massStellar() &
                  &                                         )
             call hostSpheroidComponent%abundancesStellarSet(                                           &
                  &                                           hostSpheroidComponent%abundancesStellar() &
                  &                                          +hostDiskComponent    %abundancesStellar() &
                  &                                         )
             call hostDiskComponent    %      massStellarSet(                                           &
                  &                                           0.0d0                                     &
                  &                                         )
             call hostDiskComponent    %abundancesStellarSet(                                           &
                  &                                           zeroAbundances                            &
                  &                                         )
             ! Also add stellar properties histories.
             historyDisk    =    hostDiskComponent%stellarPropertiesHistory()
             historySpheroid=hostSpheroidComponent%stellarPropertiesHistory()
             call historySpheroid      %increment                  (historyDisk    )
             call historyDisk          %reset                      (               )
             call hostSpheroidComponent%stellarPropertiesHistorySet(historySpheroid)
             call hostDiskComponent    %stellarPropertiesHistorySet( historyDisk   )
          case (doesNotMove)
             ! Do nothing.
          case default
             call Galacticus_Error_Report(                                                         &
                  &                       'Node_Component_Spheroid_Very_Simple_Satellite_Merging', &
                  &                       'unrecognized movesTo descriptor'                        &
                  &                      )
          end select
          
          ! Move the gas component of the very simple spheroid to the host.
          select case (thisMergerGasMovesTo)
          case (movesToDisk    )
             call hostDiskComponent%massGasSet              (                                           &
                  &                                          +hostDiskComponent    %      massGas()     &
                  &                                          +thisSpheroidComponent%      massGas()     &
                  &                                         )
             call hostDiskComponent%abundancesGasSet        (                                           &
                  &                                          +hostDiskComponent    %abundancesGas()     &
                  &                                          +thisSpheroidComponent%abundancesGas()     &
                  &                                         )
          case (movesToSpheroid)
             call hostSpheroidComponent%massGasSet          (                                           &
                  &                                          +hostSpheroidComponent%      massGas()     &
                  &                                          +thisSpheroidComponent%      massGas()     &
                  &                                         )
             call hostSpheroidComponent%abundancesGasSet    (                                           &
                  &                                          +hostSpheroidComponent%abundancesGas()     &
                  &                                          +thisSpheroidComponent%abundancesGas()     &
                  &                                         )
          case default
             call Galacticus_Error_Report(                                                              &
                  &                       'Node_Component_Spheroid_Very_Simple_Satellite_Merging',      &
                  &                       'unrecognized movesTo descriptor'                             &
                  &                      )
          end select
          call    thisSpheroidComponent%massGasSet          (                                           &
               &                                                                      0.0d0             &
               &                                            )
          call    thisSpheroidComponent%abundancesGasSet    (                                           &
               &                                                             zeroAbundances             &
               &                                            )
          
          ! Move the stellar component of the very simple spheroid to the host.
          select case (thisMergerStarsMoveTo)
          case (movesToDisk    )
             call hostDiskComponent%massStellarSet          (                                           &
                  &                                          +hostDiskComponent%          massStellar() &
                  &                                          +thisSpheroidComponent%      massStellar() &
                  &                                         )
             call hostDiskComponent%abundancesStellarSet    (                                           &
                  &                                          +hostDiskComponent    %abundancesStellar() &
                  &                                          +thisSpheroidComponent%abundancesStellar() &
                  &                                         )
             ! Also add stellar properties histories.
             thisHistory=thisSpheroidComponent%stellarPropertiesHistory()
             hostHistory=hostDiskComponent    %stellarPropertiesHistory()
             call hostHistory          %increment                  (thisHistory)
             call thisHistory          %reset                      (           )
             call hostDiskComponent    %stellarPropertiesHistorySet(hostHistory)
             call thisSpheroidComponent%stellarPropertiesHistorySet(thisHistory)
          case (movesToSpheroid)
             call hostSpheroidComponent%massStellarSet      (                                           &
                  &                                          +hostSpheroidComponent%      massStellar() &
                  &                                          +thisSpheroidComponent%      massStellar() &
                  &                                         )
             call hostSpheroidComponent%abundancesStellarSet(                                           &
                  &                                          +hostSpheroidComponent%abundancesStellar() &
                  &                                          +thisSpheroidComponent%abundancesStellar() &
                  &                                         )
             
             ! Also add stellar properties histories.
             thisHistory=thisSpheroidComponent%stellarPropertiesHistory()
             hostHistory=hostSpheroidComponent%stellarPropertiesHistory()
             call hostHistory          %increment                  (thisHistory)
             call thisHistory          %reset                      (           )
             call hostDiskComponent    %stellarPropertiesHistorySet(hostHistory)
             call thisSpheroidComponent%stellarPropertiesHistorySet(thisHistory)
          case default
             call Galacticus_Error_Report(                                                              &
                  &                       'Node_Component_Spheroid_Very_Simple_Satellite_Merging',      &
                  &                       'unrecognized movesTo descriptor'                             &
                  &                      )
          end select
          call    thisSpheroidComponent%massStellarSet      (                                           &
               &                                                                   0.0d0                &
               &                                            )
          call    thisSpheroidComponent%abundancesStellarSet(                                           &
               &                                                          zeroAbundances                &
               &                                            )
       end select
    end if
    return
  end subroutine Node_Component_Spheroid_Very_Simple_Satellite_Merging

  double precision function Node_Component_Spheroid_Very_Simple_SFR(self)
    !% Return the star formation rate of the very simple spheroid.
    use Star_Formation_Timescales_Spheroids
    use Dark_Matter_Halo_Scales
    implicit none
    class           (nodeComponentSpheroidVerySimple), intent(inout) :: self
    class           (darkMatterHaloScaleClass       ), pointer       :: darkMatterHaloScale_
    double precision                                                 :: spheroidDynamicalTime , gasMass, &
         &                                                              starFormationTimescale

    ! Get the gas mass.
    gasMass=self%massGas()    
    ! Get the star formation timescale.
    starFormationTimescale=Star_Formation_Timescale_Spheroid(self%hostNode)
    ! Limit the star formation timescale to a multiple of the dynamical time.
    darkMatterHaloScale_   => darkMatterHaloScale                    (             )
    spheroidDynamicalTime  =  darkMatterHaloScale_%dynamicalTimescale(self%hostNode)
    starFormationTimescale =max(starFormationTimescale,spheroidStarFormationTimescaleMinimum*spheroidDynamicalTime)
    ! If timescale is finite and gas mass is positive, then compute star formation rate.
    if (starFormationTimescale > 0.0d0 .and. gasMass > 0.0d0) then
       Node_Component_Spheroid_Very_Simple_SFR=gasMass/starFormationTimescale
    else
       Node_Component_Spheroid_Very_Simple_SFR=0.0d0
    end if
    return
  end function Node_Component_Spheroid_Very_Simple_SFR

  !# <radiusSolverPlausibility>
  !#  <unitName>Node_Component_Spheroid_Very_Simple_Radius_Solver_Plausibility</unitName>
  !# </radiusSolverPlausibility>
  subroutine Node_Component_Spheroid_Very_Simple_Radius_Solver_Plausibility(thisNode,galaxyIsPhysicallyPlausible)
    !% Determines whether the spheroid is physically plausible for radius solving tasks. Require that it have non-zero mass.
    use Dark_Matter_Halo_Scales
    implicit none
    type   (treeNode             ), intent(inout), pointer :: thisNode
    logical                       , intent(inout)          :: galaxyIsPhysicallyPlausible
    class  (nodeComponentSpheroid)               , pointer :: thisSpheroid

    ! Return immediately if our method is not selected.
    if (.not.defaultSpheroidComponent%verySimpleIsActive()) return

     ! Determine the plausibility of the current spheroid.
     thisSpheroid => thisNode%spheroid()
     select type (thisSpheroid)
     class is (nodeComponentSpheroidVerySimple)
        galaxyIsPhysicallyPlausible=(thisSpheroid%massStellar()+thisSpheroid%massGas() >= -spheroidMassToleranceAbsolute)
     end select
    return
  end subroutine Node_Component_Spheroid_Very_Simple_Radius_Solver_Plausibility

  !# <radiusSolverTask>
  !#  <unitName>Node_Component_Spheroid_Very_Simple_Radius_Solver</unitName>
  !# </radiusSolverTask>
  subroutine Node_Component_Spheroid_Very_Simple_Radius_Solver(thisNode,componentActive,specificAngularMomentumRequired,specificAngularMomentum,Radius_Get&
       &,Radius_Set,Velocity_Get,Velocity_Set)
    !% Interface for the size solver algorithm.
    use Dark_Matter_Halo_Spins
    implicit none
    type            (treeNode                                      ), intent(inout), pointer :: thisNode
    logical                                                         , intent(  out)          :: componentActive
    logical                                                         , intent(in   )          :: specificAngularMomentumRequired
    double precision                                                , intent(  out)          :: specificAngularMomentum
    procedure       (Node_Component_Spheroid_Very_Simple_Radius    ), intent(  out), pointer :: Radius_Get                     , Velocity_Get
    procedure       (Node_Component_Spheroid_Very_Simple_Radius_Set), intent(  out), pointer :: Radius_Set                     , Velocity_Set
    class           (nodeComponentSpheroid                         )               , pointer :: thisSpheroid
    class           (nodeComponentBasic                            )               , pointer :: thisBasic

    ! Determine if thisNode has an active spheroid component supported by this module.
    componentActive =  .false.
    thisSpheroid        => thisNode%spheroid()
    select type (thisSpheroid)
    class is (nodeComponentSpheroidVerySimple)
       componentActive        =  .true.
       if (specificAngularMomentumRequired) then
          thisBasic              => thisNode%basic()
          specificAngularMomentum=  Dark_Matter_Halo_Angular_Momentum(thisNode)/thisBasic%mass()
       end if
       ! Associate the pointers with the appropriate property routines.
       Radius_Get   => Node_Component_Spheroid_Very_Simple_Radius
       Radius_Set   => Node_Component_Spheroid_Very_Simple_Radius_Set
       Velocity_Get => Node_Component_Spheroid_Very_Simple_Velocity
       Velocity_Set => Node_Component_Spheroid_Very_Simple_Velocity_Set
    end select
    return
  end subroutine Node_Component_Spheroid_Very_Simple_Radius_Solver

  double precision function Node_Component_Spheroid_Very_Simple_Radius(thisNode)
    !% Return the radius of the spheroid used in structure solvers.
    implicit none
    type (treeNode             ), intent(inout), pointer :: thisNode
    class(nodeComponentSpheroid)               , pointer :: thisSpheroid

    thisSpheroid => thisNode%spheroid()
    Node_Component_Spheroid_Very_Simple_Radius=thisSpheroid%radius()
    return
  end function Node_Component_Spheroid_Very_Simple_Radius

  subroutine Node_Component_Spheroid_Very_Simple_Radius_Set(thisNode,radius)
    !% Set the radius of the spheroid used in structure solvers.
    implicit none
    type            (treeNode         ), intent(inout), pointer :: thisNode
    double precision                   , intent(in   )          :: radius
    class           (nodeComponentSpheroid)               , pointer :: thisSpheroid

    thisSpheroid => thisNode%spheroid()
    call thisSpheroid%radiusSet(radius)
    return
  end subroutine Node_Component_Spheroid_Very_Simple_Radius_Set

  double precision function Node_Component_Spheroid_Very_Simple_Velocity(thisNode)
    !% Return the circular velocity of the spheroid.
    implicit none
    type (treeNode             ), intent(inout), pointer :: thisNode
    class(nodeComponentSpheroid)               , pointer :: thisSpheroid

    thisSpheroid => thisNode%spheroid()
    Node_Component_Spheroid_Very_Simple_Velocity=thisSpheroid%velocity()
    return
  end function Node_Component_Spheroid_Very_Simple_Velocity

  subroutine Node_Component_Spheroid_Very_Simple_Velocity_Set(thisNode,velocity)
    !% Set the circular velocity of the spheroid.
    implicit none
    type            (treeNode             ), intent(inout), pointer :: thisNode
    double precision                       , intent(in   )          :: velocity
    class           (nodeComponentSpheroid)               , pointer :: thisSpheroid

    thisSpheroid => thisNode%spheroid()
    call thisSpheroid%velocitySet(velocity)
    return
  end subroutine Node_Component_Spheroid_Very_Simple_Velocity_Set

end module Node_Component_Spheroid_Very_Simple
