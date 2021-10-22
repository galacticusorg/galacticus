!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
Contains a module which implements the standard hot halo node component.
!!}

module Node_Component_Hot_Halo_Standard
  !!{
  Implements the standard hot halo node component.
  !!}
  use :: Accretion_Halos                           , only : accretionHaloClass
  use :: Chemical_Reaction_Rates                   , only : chemicalReactionRateClass
  use :: Chemical_States                           , only : chemicalStateClass
  use :: Cooling_Infall_Radii                      , only : coolingInfallRadiusClass
  use :: Cooling_Rates                             , only : coolingRateClass
  use :: Cooling_Specific_Angular_Momenta          , only : coolingSpecificAngularMomentumClass
  use :: Cosmology_Functions                       , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters                      , only : cosmologyParametersClass
  use :: Dark_Matter_Halo_Scales                   , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO                  , only : darkMatterProfileDMOClass
  use :: Hot_Halo_Mass_Distributions               , only : hotHaloMassDistributionClass
  use :: Hot_Halo_Outflows_Reincorporations        , only : hotHaloOutflowReincorporationClass
  use :: Hot_Halo_Ram_Pressure_Stripping           , only : hotHaloRamPressureStrippingClass
  use :: Hot_Halo_Ram_Pressure_Stripping_Timescales, only : hotHaloRamPressureTimescaleClass
  use :: Kind_Numbers                              , only : kind_int8
  use :: Radiation_Fields                          , only : radiationFieldClass                , radiationFieldCosmicMicrowaveBackground, radiationFieldList, radiationFieldSummation
  implicit none
  private
  public :: Node_Component_Hot_Halo_Standard_Initialize  , Node_Component_Hot_Halo_Standard_Thread_Initialize  , &
       &    Node_Component_Hot_Halo_Standard_Scale_Set   , Node_Component_Hot_Halo_Standard_Tree_Initialize    , &
       &    Node_Component_Hot_Halo_Standard_Node_Merger , Node_Component_Hot_Halo_Standard_Thread_Uninitialize, &
       &    Node_Component_Hot_Halo_Standard_Post_Step   , Node_Component_Hot_Halo_Standard_Formation          , &
       &    Node_Component_Hot_Halo_Standard_Rate_Compute, Node_Component_Hot_Halo_Standard_Pre_Evolve         , &
       &    Node_Component_Hot_Halo_Standard_Reset

  !![
  <component>
   <class>hotHalo</class>
   <name>standard</name>
   <isDefault>true</isDefault>
   <createFunction isDeferred="true" />
   <properties>
    <property>
      <name>isInitialized</name>
      <type>logical</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
    </property>
    <property>
      <name>mass</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
      <output unitsInSI="massSolar" comment="Mass of gas in the hot halo."/>
    </property>
    <property>
      <name>abundances</name>
      <type>abundances</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
      <output unitsInSI="massSolar" comment="Mass of metals in the hot phase of the hot halo."/>
    </property>
    <property>
      <name>chemicals</name>
      <type>chemicalAbundances</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
    </property>
    <property>
      <name>angularMomentum</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
      <output unitsInSI="massSolar*megaParsec*kilo" comment="Angular momentum of gas in the hot halo."/>
    </property>
    <property>
      <name>outflowedMass</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output unitsInSI="massSolar" comment="Mass of outflowed gas in the hot halo."/>
    </property>
    <property>
      <name>outflowedAngularMomentum</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output unitsInSI="massSolar*megaParsec*kilo" comment="Angular momentum of outflowed gas in the hot halo."/>
    </property>
    <property>
      <name>outflowedAbundances</name>
      <type>abundances</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output unitsInSI="massSolar" comment="Mass of metals in the outflowed phase of the hot halo."/>
    </property>
    <property>
      <name>outflowingMass</name>
      <attributes isSettable="false" isGettable="false" isEvolvable="true" isDeferred="rate" isVirtual="true" />
      <type>double</type>
      <rank>0</rank>
    </property>
    <property>
      <name>outflowingAngularMomentum</name>
      <attributes isSettable="false" isGettable="false" isEvolvable="true" isDeferred="rate" isVirtual="true" />
      <type>double</type>
      <rank>0</rank>
    </property>
    <property>
      <name>outflowingAbundances</name>
      <attributes isSettable="false" isGettable="false" isEvolvable="true" isDeferred="rate" isVirtual="true" />
      <type>abundances</type>
      <rank>0</rank>
    </property>
    <property>
      <name>unaccretedMass</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
      <output unitsInSI="massSolar" comment="Mass of gas that failed to accrete into the hot halo."/>
    </property>
    <property>
      <name>outerRadius</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" isDeferred="get" />
      <output unitsInSI="megaParsec" comment="Outer radius of the hot halo."/>
    </property>
    <property>
      <name>strippedMass</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
    </property>
    <property>
      <name>strippedAbundances</name>
      <type>abundances</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
    </property>
    <property>
      <name>hotHaloCoolingMass</name>
      <attributes isSettable="false" isGettable="false" isEvolvable="true" isDeferred="rate" bindsTo="top" isVirtual="true" />
      <type>double</type>
      <rank>0</rank>
    </property>
    <property>
      <name>hotHaloCoolingAngularMomentum</name>
      <attributes isSettable="false" isGettable="false" isEvolvable="true" isDeferred="rate" bindsTo="top" isVirtual="true" />
      <type>double</type>
      <rank>0</rank>
    </property>
    <property>
      <name>hotHaloCoolingAbundances</name>
      <attributes isSettable="false" isGettable="false" isEvolvable="true" isDeferred="rate" bindsTo="top" isVirtual="true" />
      <type>abundances</type>
      <rank>0</rank>
    </property>
    <property>
      <name>massSink</name>
      <attributes isSettable="false" isGettable="false" isEvolvable="true" isDeferred="rate" isVirtual="true" />
      <type>double</type>
      <rank>0</rank>
    </property>
    <property>
      <name>heatSource</name>
      <attributes isSettable="false" isGettable="false" isEvolvable="true" isDeferred="rate" isVirtual="true" />
      <type>double</type>
      <rank>0</rank>
    </property>
    <property>
      <name>massTotal</name>
      <attributes isSettable="false" isGettable="true" isEvolvable="false" isVirtual="true" />
      <type>double</type>
      <rank>0</rank>
      <getFunction>Node_Component_Hot_Halo_Standard_Mass_Total</getFunction>
    </property>
   </properties>
   <bindings>
     <binding method="massRemovalRate" function="Node_Component_Hot_Halo_Standard_Mass_Removal_Rate" description="Called whenever the standard hot halo component removes mass from the halo." returnType="\void" arguments="" bindsTo="component" />
     <binding method="outflowReturn" bindsTo="component" isDeferred="true" >
      <interface>
       <type>void</type>
       <self pass="true" intent="inout" />
       <argument>logical, intent(inout) :: interrupt</argument>
       <argument>procedure(interruptTask), intent(inout), pointer :: interruptProcedure</argument>
      </interface>
     </binding>
     <binding method="outerRadiusGrowthRate" bindsTo="component" isDeferred="true" >
      <interface>
       <type>double</type>
       <rank>0</rank>
       <self pass="true" intent="inout" />
      </interface>
     </binding>
   </bindings>
   <functions>objects.nodes.components.hot_halo.standard.bound_functions.inc</functions>
  </component>
  !!]

  ! Objects used by this component.
  class(cosmologyFunctionsClass            ), pointer :: cosmologyFunctions_
  class(darkMatterHaloScaleClass           ), pointer :: darkMatterHaloScale_
  class(darkMatterProfileDMOClass          ), pointer :: darkMatterProfileDMO_
  class(coolingSpecificAngularMomentumClass), pointer :: coolingSpecificAngularMomentum_
  class(coolingInfallRadiusClass           ), pointer :: coolingInfallRadius_
  class(hotHaloMassDistributionClass       ), pointer :: hotHaloMassDistribution_
  class(accretionHaloClass                 ), pointer :: accretionHalo_
  class(chemicalReactionRateClass          ), pointer :: chemicalReactionRate_
  class(hotHaloRamPressureStrippingClass   ), pointer :: hotHaloRamPressureStripping_
  class(hotHaloRamPressureTimescaleClass   ), pointer :: hotHaloRamPressureTimescale_
  class(hotHaloOutflowReincorporationClass ), pointer :: hotHaloOutflowReincorporation_
  class(chemicalStateClass                 ), pointer :: chemicalState_
  class(coolingRateClass                   ), pointer :: coolingRate_
  class(cosmologyParametersClass           ), pointer :: cosmologyParameters_
  !$omp threadprivate(cosmologyFunctions_,darkMatterHaloScale_,darkMatterProfileDMO_,coolingSpecificAngularMomentum_,coolingInfallRadius_,hotHaloMassDistribution_,accretionHalo_,chemicalReactionRate_,chemicalState_,hotHaloRamPressureStripping_,hotHaloRamPressureTimescale_,coolingRate_,cosmologyParameters_,hotHaloOutflowReincorporation_)

  ! Internal count of abundances and chemicals.
  integer                                                             :: abundancesCount                            , chemicalsCount

  ! Configuration variables.
  logical                                                             :: hotHaloExcessHeatDrivesOutflow             , hotHaloAngularMomentumAlwaysGrows
  double precision                                                    :: hotHaloExpulsionRateMaximum                , hotHaloOutflowStrippingEfficiency

  ! Quantities stored to avoid repeated computation.
  integer         (kind=kind_int8                        )            :: uniqueIDPrevious
  logical                                                             :: gotAngularMomentumCoolingRate      =.false., gotCoolingRate                   =.false., &
       &                                                                 gotOuterRadiusGrowthRate           =.false.
  double precision                                                    :: angularMomentumHeatingRateRemaining        , rateCooling                              , &
       &                                                                 massHeatingRateRemaining                   , outerRadiusGrowthRateStored
  !$omp threadprivate(gotCoolingRate,gotAngularMomentumCoolingRate,gotOuterRadiusGrowthRate,rateCooling,massHeatingRateRemaining,angularMomentumHeatingRateRemaining,outerRadiusGrowthRateStored,uniqueIDPrevious)
  ! Radiation structure.
  type           (radiationFieldSummation                ), pointer   :: radiation
  type           (radiationFieldCosmicMicrowaveBackground), pointer   :: radiationCosmicMicrowaveBackground
  class          (radiationFieldClass                    ), pointer   :: radiationIntergalacticBackground
  type           (radiationFieldList                     ), pointer   :: radiationFieldList_
  !$omp threadprivate(radiation,radiationCosmicMicrowaveBackground,radiationIntergalacticBackground,radiationFieldList_)

  ! Tracked properties control.
  logical                                                             :: hotHaloTrackStrippedGas

  ! Parameters controlling absolute tolerance scales.
  double precision                                        , parameter :: scaleMassRelative                 =1.0d-3
  double precision                                        , parameter :: scaleRadiusRelative               =1.0d-1

contains

  !![
  <nodeComponentInitializationTask>
   <unitName>Node_Component_Hot_Halo_Standard_Initialize</unitName>
  </nodeComponentInitializationTask>
  !!]
  subroutine Node_Component_Hot_Halo_Standard_Initialize(parameters_)
    !!{
    Initializes the standard hot halo component module.
    !!}
    use :: Abundances_Structure                 , only : Abundances_Property_Count           , abundances
    use :: Chemical_Abundances_Structure        , only : Chemicals_Property_Count
    use :: Galacticus_Error                     , only : Galacticus_Error_Report
    use :: Galacticus_Nodes                     , only : defaultHotHaloComponent             , nodeComponentHotHaloStandard
    use :: ISO_Varying_String                   , only : var_str                             , varying_string                 , char
    use :: Input_Parameters                     , only : inputParameter                      , inputParameters
    use :: Node_Component_Hot_Halo_Standard_Data, only : currentNode                         , formationNode                  , hotHaloAngularMomentumLossFraction, hotHaloCoolingFromNode   , &
          &                                              hotHaloNodeMergerLimitBaryonFraction, hotHaloOutflowReturnOnFormation, starveSatellites                  , starveSatellitesOutflowed
    implicit none
    type(inputParameters             ), intent(inout) :: parameters_
    type(varying_string              )                :: hotHaloCoolingFromText
    type(nodeComponentHotHaloStandard)                :: hotHalo

    ! Initialize the module if necessary.
    !$omp critical (Node_Component_Hot_Halo_Standard_Initialize)
    if (defaultHotHaloComponent%standardIsActive()) then

       ! Get numbers of abundance and chemicals properties.
       abundancesCount=Abundances_Property_Count()
       chemicalsCount =Chemicals_Property_Count ()

       ! Determine whether satellite nodes will be starved of gas.
       !![
       <inputParameter>
         <name>starveSatellites</name>
         <defaultValue>.false.</defaultValue>
         <description>Specifies whether or not the hot halo should be removed (``starved'') when a node becomes a satellite.</description>
         <source>parameters_</source>
       </inputParameter>
       !!]

       !![
       <inputParameter>
         <name>starveSatellitesOutflowed</name>
         <defaultValue>.false.</defaultValue>
         <description>Specifies whether or not the outflowed hot halo should be removed (``starved'') when a node becomes a satellite.</description>
         <source>parameters_</source>
       </inputParameter>
       !!]

       ! Determine whether stripped material should be tracked.
       !![
       <inputParameter>
         <name>hotHaloTrackStrippedGas</name>
         <defaultValue>.true.</defaultValue>
         <description>Specifies whether or not gas stripped from the hot halo should be tracked.</description>
         <source>parameters_</source>
       </inputParameter>
       !!]

       ! Determine whether outflowed gas should be restored to the hot reservoir on halo formation events.
       !![
       <inputParameter>
         <name>hotHaloOutflowReturnOnFormation</name>
         <defaultValue>.false.</defaultValue>
         <description>Specifies whether or not outflowed gas should be returned to the hot reservoir on halo formation events.</description>
         <source>parameters_</source>
       </inputParameter>
       !!]

       ! Determine whether negative angular momentum accretion rates onto the halo should be treated as positive for the purposes
       ! of computing the hot halo angular momentum.
       !![
       <inputParameter>
         <name>hotHaloAngularMomentumAlwaysGrows</name>
         <defaultValue>.false.</defaultValue>
         <description>Specifies whether or not negative rates of accretion of angular momentum into the hot halo will be treated as positive
            for the purposes of computing the hot halo angular momentum.</description>
         <source>parameters_</source>
       </inputParameter>
       !!]

       ! Determine whether the angular momentum of cooling gas should be computed from the "current node" or the "formation node".
       !![
       <inputParameter>
         <name>hotHaloCoolingFromNode</name>
         <defaultValue>var_str('currentNode')</defaultValue>
         <description>Specifies whether the angular momentum of cooling gas should be computed from the ``current node'' or the ``formation node''.</description>
         <source>parameters_</source>
         <variable>hotHaloCoolingFromText</variable>
       </inputParameter>
       !!]
       select case (char(hotHaloCoolingFromText))
       case ("currentNode"  )
          hotHaloCoolingFromNode=currentNode
       case ("formationNode")
          hotHaloCoolingFromNode=formationNode
       case default
          call Galacticus_Error_Report('hotHaloCoolingFromNode must be one of "currentNode" or "formationNode"'//{introspection:location})
       end select

       ! Determine whether excess heating of the halo will drive an outflow.
       !![
       <inputParameter>
         <name>hotHaloExcessHeatDrivesOutflow</name>
         <defaultValue>.true.</defaultValue>
         <description>Specifies whether heating of the halo in excess of its cooling rate will drive an outflow from the halo.</description>
         <source>parameters_</source>
       </inputParameter>
       !!]

       ! Get efficiency with which outflowing gas is stripped from the hot halo.
       !![
       <inputParameter>
         <name>hotHaloOutflowStrippingEfficiency</name>
         <defaultValue>0.1d0</defaultValue>
         <description>Specifies the efficiency with which outflowing gas is stripped from the hot halo, following the prescription of \citeauthor{font_colours_2008}~(\citeyear{font_colours_2008}; i.e. this is the parameter $\epsilon_\mathrm{strip}$ in their eqn.~6).</description>
         <source>parameters_</source>
       </inputParameter>
       !!]

       ! Get the maximum rate (in units of halo inverse dynamical time) at which gas can be expelled from the halo.
       !![
       <inputParameter>
         <name>hotHaloExpulsionRateMaximum</name>
         <defaultValue>1.0d0</defaultValue>
         <description>Specifies the maximum rate at which mass can be expelled from the hot halo in units of the inverse halo dynamical time.</description>
         <source>parameters_</source>
       </inputParameter>
       !!]

       ! Get fraction of angular momentum that is lost during cooling/infall.
       !![
       <inputParameter>
         <name>hotHaloAngularMomentumLossFraction</name>
         <defaultValue>0.3d0</defaultValue>
         <description>Specifies the fraction of angular momentum that is lost from cooling/infalling gas.</description>
         <source>parameters_</source>
       </inputParameter>
       !!]

       ! Get option controlling limiting of baryon fraction during node mergers.
       !![
       <inputParameter>
         <name>hotHaloNodeMergerLimitBaryonFraction</name>
         <defaultValue>.false.</defaultValue>
         <description>Controls whether the hot gas content of nodes should be limited to not exceed the universal baryon fraction at node
           merger events. If set to {\normalfont \ttfamily true}, hot gas (and angular momentum, abundances, and chemicals proportionally) will be
           removed from the merged halo to the unaccreted gas reservoir to limit the baryonic mass to the universal baryon
           fraction where possible.</description>
         <source>parameters_</source>
       </inputParameter>
       !!]
       ! Bind the outer radius get function.
       call hotHalo%                  outerRadiusFunction(Node_Component_Hot_Halo_Standard_Outer_Radius              )
       ! Bind the heat source pipe to the function that will handle heat input to the hot halo.
       call hotHalo%               heatSourceRateFunction(Node_Component_Hot_Halo_Standard_Heat_Source               )
       ! Bind outflowing material pipes to the functions that will handle input of outflowing material to the hot halo.
       call hotHalo%           outflowingMassRateFunction(Node_Component_Hot_Halo_Standard_Outflowing_Mass_Rate      )
       call hotHalo%outflowingAngularMomentumRateFunction(Node_Component_Hot_Halo_Standard_Outflowing_Ang_Mom_Rate   )
       call hotHalo%     outflowingAbundancesRateFunction(Node_Component_Hot_Halo_Standard_Outflowing_Abundances_Rate)
       ! Bind a creation function.
       call hotHalo%                    createFunctionSet(Node_Component_Hot_Halo_Standard_Initializor               )
       ! Bind the mass sink function.
       call hotHalo%                 massSinkRateFunction(Node_Component_Hot_Halo_Standard_Mass_Sink                 )
       ! Bind the outflow return function.
       call hotHalo%                outflowReturnFunction(Node_Component_Hot_Halo_Standard_Outflow_Return            )
       ! Bind the outer radius growth rate function.
       call hotHalo%        outerRadiusGrowthRateFunction(Node_Component_Hot_Halo_Standard_Outer_Radius_Growth_Rate  )
    end if
    !$omp end critical (Node_Component_Hot_Halo_Standard_Initialize)
    return
  end subroutine Node_Component_Hot_Halo_Standard_Initialize

  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Node_Component_Hot_Halo_Standard_Thread_Initialize</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Node_Component_Hot_Halo_Standard_Thread_Initialize(parameters_)
    !!{
    Initializes the tree node hot halo methods module.
    !!}
    use :: Events_Hooks    , only : nodePromotionEvent                   , satelliteMergerEvent    , postEvolveEvent, openMPThreadBindingAtLevel, &
         &                          dependencyRegEx                      , dependencyDirectionAfter
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: Galacticus_Nodes, only : defaultHotHaloComponent
    use :: Input_Parameters, only : inputParameter                       , inputParameters
    use :: Radiation_Fields, only : radiationFieldIntergalacticBackground
    implicit none
    type(inputParameters), intent(inout) :: parameters_
    type(dependencyRegEx), dimension(1)  :: dependencies

    ! Check if this implementation is selected. Define the radiation component to include both the CMB and the intergalactic background if it is.
    if (defaultHotHaloComponent%standardIsActive()) then
       !![
       <objectBuilder class="cosmologyParameters"            name="cosmologyParameters_"            source="parameters_"/>
       <objectBuilder class="cosmologyFunctions"             name="cosmologyFunctions_"             source="parameters_"/>
       <objectBuilder class="darkMatterHaloScale"            name="darkMatterHaloScale_"            source="parameters_"/>
       <objectBuilder class="darkMatterProfileDMO"           name="darkMatterProfileDMO_"           source="parameters_"/>
       <objectBuilder class="coolingSpecificAngularMomentum" name="coolingSpecificAngularMomentum_" source="parameters_"/>
       <objectBuilder class="coolingInfallRadius"            name="coolingInfallRadius_"            source="parameters_"/>
       <objectBuilder class="hotHaloMassDistribution"        name="hotHaloMassDistribution_"        source="parameters_"/>
       <objectBuilder class="accretionHalo"                  name="accretionHalo_"                  source="parameters_"/>
       <objectBuilder class="chemicalReactionRate"           name="chemicalReactionRate_"           source="parameters_"/>
       <objectBuilder class="chemicalState"                  name="chemicalState_"                  source="parameters_"/>
       <objectBuilder class="hotHaloRamPressureStripping"    name="hotHaloRamPressureStripping_"    source="parameters_"/>
       <objectBuilder class="hotHaloRamPressureTimescale"    name="hotHaloRamPressureTimescale_"    source="parameters_"/>
       <objectBuilder class="hotHaloOutflowReincorporation"  name="hotHaloOutflowReincorporation_"  source="parameters_"/>
       <objectBuilder class="coolingRate"                    name="coolingRate_"                    source="parameters_"/>
       !!]
       dependencies(1)=dependencyRegEx(dependencyDirectionAfter,'^remnantStructure:')
       call nodePromotionEvent  %attach(defaultHotHaloComponent,nodePromotion  ,openMPThreadBindingAtLevel,label='nodeComponentHotHaloStandard'                          )
       call satelliteMergerEvent%attach(defaultHotHaloComponent,satelliteMerger,openMPThreadBindingAtLevel,label='nodeComponentHotHaloStandard',dependencies=dependencies)
       call postEvolveEvent     %attach(defaultHotHaloComponent,postEvolve     ,openMPThreadBindingAtLevel,label='nodeComponentHotHaloStandard'                          )
       allocate(radiation                         )
       allocate(radiationFieldList_               )
       allocate(radiationCosmicMicrowaveBackground)
       !![
       <referenceConstruct object="radiationCosmicMicrowaveBackground" constructor="radiationFieldCosmicMicrowaveBackground(cosmologyFunctions_)"/>
       !!]
       radiationFieldList_%radiationField_ => radiationCosmicMicrowaveBackground
       if (parameters_%isPresent('radiationFieldIntergalacticBackground')) then
          !![
          <objectBuilder class="radiationField" name="radiationIntergalacticBackground" parameterName="radiationFieldIntergalacticBackground" source="parameters_"/>
          !!]
          select type (radiationIntergalacticBackground)
          class is (radiationFieldIntergalacticBackground)
             ! This is as expected.
          class default
             call Galacticus_Error_Report('radiation field is not of the intergalactic background class'//{introspection:location})
          end select
          allocate(radiationFieldList_%next)
          radiationFieldList_%next%radiationField_ => radiationIntergalacticBackground
       else
          radiationIntergalacticBackground => null()
       end if
       !![
       <referenceConstruct object="radiation" constructor="radiationFieldSummation(radiationFieldList_)"/>
       !!]
       nullify(radiationFieldList_)
    end if
    return
  end subroutine Node_Component_Hot_Halo_Standard_Thread_Initialize

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Node_Component_Hot_Halo_Standard_Thread_Uninitialize</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Node_Component_Hot_Halo_Standard_Thread_Uninitialize()
    !!{
    Uninitializes the tree node hot halo methods module.
    !!}
    use :: Events_Hooks    , only : nodePromotionEvent     , satelliteMergerEvent, postEvolveEvent
    use :: Galacticus_Nodes, only : defaultHotHaloComponent
    implicit none

    if (defaultHotHaloComponent%standardIsActive()) then
       !![
       <objectDestructor name="cosmologyParameters_"              />
       <objectDestructor name="cosmologyFunctions_"               />
       <objectDestructor name="darkMatterHaloScale_"              />
       <objectDestructor name="darkMatterProfileDMO_"             />
       <objectDestructor name="coolingSpecificAngularMomentum_"   />
       <objectDestructor name="coolingInfallRadius_"              />
       <objectDestructor name="hotHaloMassDistribution_"          />
       <objectDestructor name="accretionHalo_"                    />
       <objectDestructor name="chemicalReactionRate_"             />
       <objectDestructor name="chemicalState_"                    />
       <objectDestructor name="hotHaloRamPressureStripping_"      />
       <objectDestructor name="hotHaloRamPressureTimescale_"      />
       <objectDestructor name="hotHaloOutflowReincorporation_"    />
       <objectDestructor name="coolingRate_"                      />
       <objectDestructor name="radiationIntergalacticBackground"  />
       <objectDestructor name="radiationCosmicMicrowaveBackground"/>
       <objectDestructor name="radiation"                         />
       !!]
       if (nodePromotionEvent  %isAttached(defaultHotHaloComponent,nodePromotion  )) call nodePromotionEvent  %detach(defaultHotHaloComponent,nodePromotion  )
       if (satelliteMergerEvent%isAttached(defaultHotHaloComponent,satelliteMerger)) call satelliteMergerEvent%detach(defaultHotHaloComponent,satelliteMerger)
       if (postEvolveEvent     %isAttached(defaultHotHaloComponent,postEvolve     )) call postEvolveEvent     %detach(defaultHotHaloComponent,postEvolve     )
    end if
    return
  end subroutine Node_Component_Hot_Halo_Standard_Thread_Uninitialize

  !![
  <calculationResetTask>
  <unitName>Node_Component_Hot_Halo_Standard_Reset</unitName>
  </calculationResetTask>
  !!]
  subroutine Node_Component_Hot_Halo_Standard_Reset(node)
    !!{
    Remove memory of stored computed values as we're about to begin computing derivatives anew.
    !!}
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    type(treeNode), intent(inout) :: node

    uniqueIDPrevious             =node%uniqueID()
    gotCoolingRate               =.false.
    gotAngularMomentumCoolingRate=.false.
    gotOuterRadiusGrowthRate     =.false.
    return
  end subroutine Node_Component_Hot_Halo_Standard_Reset

  double precision function Node_Component_Hot_Halo_Standard_Outer_Radius(self)
    !!{
    Return the outer radius in the standard hot halo.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHaloStandard, treeNode
    implicit none
    class           (nodeComponentHotHaloStandard), intent(inout) :: self
    type            (treeNode                    ), pointer       :: selfHost
    double precision                                              :: radiusVirial

    selfHost     => self%host                        (        )
    radiusVirial =  darkMatterHaloScale_%virialRadius(selfHost)
    Node_Component_Hot_Halo_Standard_Outer_Radius=max(                                      &
         &                                            min(                                  &
         &                                                self%outerRadiusValue()         , &
         &                                                                    radiusVirial  &
         &                                               )                                , &
         &                                                scaleRadiusRelative*radiusVirial  &
         &                                            )
    return
  end function Node_Component_Hot_Halo_Standard_Outer_Radius

  !![
  <postStepTask>
   <unitName>Node_Component_Hot_Halo_Standard_Post_Step</unitName>
  <after>Node_Component_Basic_Standard_Post_Step</after>
  </postStepTask>
  !!]
  subroutine Node_Component_Hot_Halo_Standard_Post_Step(node,status)
    !!{
    Do processing of the node required after evolution.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHalo, nodeComponentHotHaloStandard, treeNode, defaultHotHaloComponent
    use :: Interface_GSL   , only : GSL_Failure
    implicit none
    type   (treeNode            ), intent(inout), pointer :: node
    integer                      , intent(inout)          :: status
    class  (nodeComponentHotHalo)               , pointer :: hotHalo

    ! Return immediately if this class is not in use.
    if (.not.defaultHotHaloComponent%standardIsActive()) return
    ! Limit hot gas mass, and outer radius to be non-negative.
    hotHalo => node%hotHalo()
    select type (hotHalo)
    class is (nodeComponentHotHaloStandard)
       if (hotHalo%       mass() < 0.0d0) then
          status=GSL_Failure
          call hotHalo%       massSet(0.0d0)
       end if
       if (hotHalo%outerRadius() < 0.0d0) then
          status=GSL_Failure
          call hotHalo%outerRadiusSet(0.0d0)
       end if
    end select
    return
  end subroutine Node_Component_Hot_Halo_Standard_Post_Step

  subroutine postEvolve(self,node)
    !!{
    Do processing of the node required after evolution.
    !!}
    use :: Abundances_Structure                 , only : zeroAbundances
    use :: Galacticus_Nodes                     , only : nodeComponentHotHalo, nodeComponentHotHaloStandard, nodeComponentSpin, treeNode
    use :: Node_Component_Hot_Halo_Standard_Data, only : starveSatellites    , starveSatellitesOutflowed
    implicit none
    class(*                   ), intent(inout) :: self
    type (treeNode            ), intent(inout) :: node
    type (treeNode            ), pointer       :: nodeParent
    class(nodeComponentHotHalo), pointer       :: hotHaloParent, hotHalo
    class(nodeComponentSpin   ), pointer       :: spinParent
    !$GLC attributes unused :: self

    ! Process hot gas for satellites.
    if (node%isSatellite()) then
       if (starveSatellites.or.starveSatellitesOutflowed) then
          hotHalo => node%hotHalo(autoCreate=.true.)
          select type (hotHalo)
          class is (nodeComponentHotHaloStandard)
             ! Transfer any outflowed gas to the hot halo of the parent node.
             nodeParent => node%parent
             do while (nodeParent%isSatellite())
                nodeParent => nodeParent%parent
             end do
             hotHaloParent => nodeParent%hotHalo (autoCreate=.true.)
             spinParent    => nodeParent%spin    (                 )
             call hotHaloParent%outflowedAngularMomentumSet(                                                           &
                  &                                         +hotHaloParent       %outflowedAngularMomentum(          ) &
                  &                                         +hotHalo             %outflowedMass           (          ) &
                  &                                         *spinParent          %spin                    (          ) &
                  &                                         *darkMatterHaloScale_%virialRadius            (nodeParent) &
                  &                                         *darkMatterHaloScale_%virialVelocity          (nodeParent) &
                  &                                        )
             call hotHalo      %outflowedAngularMomentumSet(                                                           &
                  &                                         +0.0d0                                                     &
                  &                                        )
             call hotHaloParent%outflowedMassSet           (                                                           &
                  &                                         +hotHaloParent       %outflowedMass           (          ) &
                  &                                         +hotHalo             %outflowedMass           (          ) &
                  &                                        )
             call hotHalo      %outflowedMassSet           (                                                           &
                  &                                         +0.0d0                                                     &
                  &                                        )
             call hotHaloParent%outflowedAbundancesSet     (                                                           &
                  &                                         +hotHaloParent       %outflowedAbundances     (          ) &
                  &                                         +hotHalo             %outflowedAbundances     (          ) &
                  &                                        )
             call hotHalo      %outflowedAbundancesSet     (                                                           &
                  &                                         +zeroAbundances                                            &
                  &                                        )
          end select
       end if
       ! Check if stripped mass is being tracked.
       if (hotHaloTrackStrippedGas) then
          hotHalo => node%hotHalo()
          select type (hotHalo)
          class is (nodeComponentHotHaloStandard)
             ! Transfer any stripped gas to the host halo.
             nodeParent => node%parent
             do while (nodeParent%isSatellite())
                nodeParent => nodeParent%parent
             end do
             call Node_Component_Hot_Halo_Standard_Create(nodeParent)
             hotHaloParent => nodeParent%hotHalo (autoCreate=.true.)
             spinParent    => nodeParent%spin    (                 )
             call hotHaloParent%outflowedAngularMomentumSet(                                                           &
                  &                                         +hotHaloParent       %outflowedAngularMomentum(          ) &
                  &                                         +hotHalo             %strippedMass            (          ) &
                  &                                         *spinParent          %spin                    (          ) &
                  &                                         *darkMatterHaloScale_%virialRadius            (nodeParent) &
                  &                                         *darkMatterHaloScale_%virialVelocity          (nodeParent) &
                  &                                        )
             call hotHaloParent%outflowedMassSet           (                                                           &
                  &                                         +hotHaloParent       %outflowedMass           (          ) &
                  &                                         +hotHalo             %strippedMass            (          ) &
                  &                                        )
             call   hotHalo%strippedMassSet                (                                                           &
                  &                                         +0.0d0                                                     &
                  &                                        )
             call hotHaloParent%outflowedAbundancesSet     (                                                           &
                  &                                         +hotHaloParent       %outflowedAbundances     (          ) &
                  &                                         +hotHalo             %strippedAbundances      (          ) &
                  &                                        )
             call hotHalo%strippedAbundancesSet            (                                                           &
                  &                                         +zeroAbundances                                            &
                  &                                        )
          end select
       end if
    end if
    return
  end subroutine postEvolve

  subroutine Node_Component_Hot_Halo_Standard_Strip_Gas_Rate(node,gasMassRate,interrupt,interruptProcedure)
    !!{
    Add gas stripped from the hot halo to the stripped gas reservoirs under the assumption of uniformly distributed properties
    (e.g. fully-mixed metals).
    !!}
    use :: Galacticus_Nodes, only : interruptTask, nodeComponentHotHalo, nodeComponentHotHaloStandard, treeNode
    implicit none
    type            (treeNode            ), intent(inout)          :: node
    double precision                      , intent(in   )          :: gasMassRate
    logical                               , intent(inout)          :: interrupt
    procedure       (interruptTask       ), intent(inout), pointer :: interruptProcedure
    class           (nodeComponentHotHalo)               , pointer :: hotHalo
    double precision                                               :: gasMass

    ! Exit immediately for zero rate.
    if (gasMassRate == 0.0d0) return

    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    select type (hotHalo)
    class is (nodeComponentHotHaloStandard)
       ! Get the gas mass present.
       gasMass=hotHalo%mass()
       ! If gas is present, adjust the rates.
       if (gasMass > 0.0d0) then
          ! Mass.
          call hotHalo%      strippedMassRate(                     gasMassRate        ,interrupt,interruptProcedure)
          ! Metal abundances.
          call hotHalo%strippedAbundancesRate(hotHalo%abundances()*gasMassRate/gasMass,interrupt,interruptProcedure)
       end if
    end select
    return
  end subroutine Node_Component_Hot_Halo_Standard_Strip_Gas_Rate

  subroutine Node_Component_Hot_Halo_Standard_Heat_Source(hotHalo,rate,interrupt,interruptProcedure)
    !!{
    An incoming pipe for sources of heating to the hot halo.
    !!}
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: Galacticus_Nodes, only : interruptTask          , nodeComponentHotHalo, nodeComponentHotHaloStandard, treeNode
    implicit none
    class           (nodeComponentHotHalo), intent(inout)                    :: hotHalo
    double precision                      , intent(in   )                    :: rate
    logical                               , intent(inout), optional          :: interrupt
    procedure       (interruptTask       ), intent(inout), optional, pointer :: interruptProcedure
    type            (treeNode            )                         , pointer :: node
    double precision                                                         :: excessMassHeatingRate, inputMassHeatingRate, massHeatingRate

     ! Trap cases where an attempt is made to remove energy via this input function.
     if (rate < 0.0d0) call Galacticus_Error_Report('attempt to remove energy via heat source pipe to hot halo'//{introspection:location})

     ! Get the node associated with this hot halo component.
     node => hotHalo%host()

     ! Ensure that the cooling rate has been computed.
     call Node_Component_Hot_Halo_Standard_Cooling_Rate(node)

     ! Compute the input mass heating rate from the input energy heating rate.
     inputMassHeatingRate=rate/darkMatterHaloScale_%virialVelocity(node)**2

     ! Limit the mass heating rate such that it never exceeds the remaining budget.
     massHeatingRate=min(inputMassHeatingRate,massHeatingRateRemaining)

     ! Update the remaining budget of allowed mass heating rate.
     if (massHeatingRateRemaining-massHeatingRate <= 0.0d0) massHeatingRate=massHeatingRateRemaining
     massHeatingRateRemaining=max(massHeatingRateRemaining-massHeatingRate,0.0d0)

     ! Call routine to apply this mass heating rate to all hot halo cooling pipes.
     call Node_Component_Hot_Halo_Standard_Push_To_Cooling_Pipes(node,-massHeatingRate,interrupt,interruptProcedure)

     ! If requested, compute the rate at which an outflow is driven from the halo by excess heating.
     if (hotHaloExcessHeatDrivesOutflow) then

        ! Compute the excess mass heating rate (i.e. that beyond which is being used to offset the cooling rate).
        excessMassHeatingRate=inputMassHeatingRate-massHeatingRate

        ! Remove any excess mass heating rate from the halo.
        select type (hotHalo)
        class is (nodeComponentHotHaloStandard)
           call Node_Component_Hot_Halo_Standard_Push_From_Halo(hotHalo,excessMassHeatingRate)
        end select

     end if

     return
   end subroutine Node_Component_Hot_Halo_Standard_Heat_Source

  subroutine Node_Component_Hot_Halo_Standard_Push_To_Cooling_Pipes(node,massRate,interrupt,interruptProcedure)
    !!{
    Push mass through the cooling pipes (along with appropriate amounts of metals and angular momentum) at the given rate.
    !!}
    use :: Abundances_Structure                 , only : abundances             , max                 , operator(*)
    use :: Galacticus_Error                     , only : Galacticus_Error_Report
    use :: Galacticus_Nodes                     , only : interruptTask          , nodeComponentHotHalo, nodeComponentHotHaloStandard      , treeNode
    use :: Node_Component_Hot_Halo_Standard_Data, only : currentNode            , formationNode       , hotHaloAngularMomentumLossFraction, hotHaloCoolingFromNode
    implicit none
    type            (treeNode                ), intent(inout)          , target  :: node
    double precision                          , intent(in   )                    :: massRate
    logical                                   , intent(inout), optional          :: interrupt
    procedure       (interruptTask           ), intent(inout), optional, pointer :: interruptProcedure
    type            (treeNode                )                         , pointer :: nodeCooling
    class           (nodeComponentHotHalo    )                         , pointer :: hotHaloCooling            , hotHalo
    type            (abundances              ), save                             :: abundancesCoolingRate
    !$omp threadprivate(abundancesCoolingRate)
    double precision                                                             :: angularMomentumCoolingRate, infallRadius

    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    select type (hotHalo)
    class is (nodeComponentHotHaloStandard)

       ! Ignore zero rates.
       if (massRate /= 0.0d0 .and. hotHalo%mass() > 0.0d0 .and. hotHalo%angularMomentum() > 0.0d0) then
          ! Remove mass from the hot component.
          call    hotHalo%massRate       (-massRate)
          call    hotHalo%massRemovalRate(+massRate)
          ! Pipe the mass rate to whichever component claimed it.
          if (hotHalo%hotHaloCoolingMassRateIsAttached()) &
               & call hotHalo%hotHaloCoolingMassRate(+massRate,interrupt,interruptProcedure)
          ! Find the node to use for cooling calculations.
          select case (hotHaloCoolingFromNode)
          case (currentNode  )
             nodeCooling => node
          case (formationNode)
             nodeCooling => node%formationNode
          case default
             nodeCooling => null()
             call Galacticus_Error_Report('unknown "hotHaloCoolingFromNode" - this should not happen'//{introspection:location})
          end select
          infallRadius              =coolingInfallRadius_%radius(node)
          angularMomentumCoolingRate=massRate*coolingSpecificAngularMomentum_%angularMomentumSpecific(nodeCooling,infallRadius)
          if (.not.gotAngularMomentumCoolingRate) then
             angularMomentumHeatingRateRemaining=rateCooling*coolingSpecificAngularMomentum_%angularMomentumSpecific(nodeCooling,infallRadius)
             gotAngularMomentumCoolingRate=.true.
          end if
          if (massRate < 0.0d0) then
             if (massHeatingRateRemaining == 0.0d0) then
                angularMomentumCoolingRate=-angularMomentumHeatingRateRemaining
                angularMomentumHeatingRateRemaining=0.0d0
             else
                angularMomentumCoolingRate=max(angularMomentumCoolingRate,-angularMomentumHeatingRateRemaining)
                angularMomentumHeatingRateRemaining=angularMomentumHeatingRateRemaining+angularMomentumCoolingRate
             end if
          end if
          call    hotHalo%angularMomentumRate       (     -angularMomentumCoolingRate                                                                                  )
          ! Pipe the cooling rate to which ever component claimed it.
          if (hotHalo%hotHaloCoolingAngularMomentumRateIsAttached()) &
               & call hotHalo%hotHaloCoolingAngularMomentumRate(sign(+angularMomentumCoolingRate*(1.0d0-hotHaloAngularMomentumLossFraction),massRate),interrupt,interruptProcedure)
          ! Get the rate of change of abundances.
          hotHaloCooling => nodeCooling%hotHalo()
          abundancesCoolingRate=hotHaloCooling%abundances()
          abundancesCoolingRate=massRate*abundancesCoolingRate/hotHaloCooling%mass()
          call    hotHalo%abundancesRate       (-abundancesCoolingRate )
          ! Pipe the cooling rate to which ever component claimed it.
          if (hotHalo%hotHaloCoolingAbundancesRateIsAttached()) &
               & call hotHalo%hotHaloCoolingAbundancesRate(+abundancesCoolingRate,interrupt,interruptProcedure)
       end if
    end select
    return
  end subroutine Node_Component_Hot_Halo_Standard_Push_To_Cooling_Pipes

  subroutine Node_Component_Hot_Halo_Standard_Push_From_Halo(hotHalo,massRate)
    !!{
    Push mass from the hot halo into an infinite sink (along with appropriate amounts of metals, chemicals and angular momentum) at the given rate.
    !!}
    use :: Abundances_Structure         , only : abundances
    use :: Chemical_Abundances_Structure, only : chemicalAbundances
    use :: Galacticus_Nodes             , only : nodeComponentHotHaloStandard, treeNode
    implicit none
    class           (nodeComponentHotHaloStandard), intent(inout) :: hotHalo
    double precision                              , intent(in   ) :: massRate
    type            (treeNode                    ), pointer       :: node
    type            (abundances                  ), save          :: abundancesRates
    type            (chemicalAbundances          ), save          :: chemicalsRates
    !$omp threadprivate(abundancesRates,chemicalsRates)
    double precision                                              :: angularMomentumRate, massRateLimited

    ! Ignore zero rates.
    if (massRate /= 0.0d0 .and. hotHalo%mass() > 0.0d0) then
       ! Limit the mass expulsion rate to a fraction of the halo dynamical timescale.
       node => hotHalo%hostNode
       massRateLimited=min(massRate,hotHaloExpulsionRateMaximum*hotHalo%mass()/darkMatterHaloScale_%dynamicalTimescale(node))
       ! Get the rate of change of abundances, chemicals, and angular momentum.
       abundancesRates    =hotHalo%abundances     ()*massRateLimited/hotHalo%mass()
       angularMomentumRate=hotHalo%angularMomentum()*massRateLimited/hotHalo%mass()
       chemicalsRates     =hotHalo%chemicals      ()*massRateLimited/hotHalo%mass()
       call hotHalo%    massRemovalRate(+    massRateLimited)
       call hotHalo%           massRate(-    massRateLimited)
       call hotHalo%     abundancesRate(-    abundancesRates)
       call hotHalo%angularMomentumRate(-angularMomentumRate)
       call hotHalo%      chemicalsRate(-     chemicalsRates)
       ! If this node is a satellite and stripped gas is being tracked, move mass and abundances to the stripped reservoir.
       if (node%isSatellite().and.hotHaloTrackStrippedGas) then
          call hotHalo%      strippedMassRate(+massRateLimited)
          call hotHalo%strippedAbundancesRate(+abundancesRates)
       end if
      ! Trigger an event to allow other processes to respond to this action.
      !![
      <eventHook name="hotHaloMassEjection">
	<import>
	  <module name="Galacticus_Nodes" symbols="nodeComponentHotHalo"/>
	</import>
	<interface>
	  class           (nodeComponentHotHalo), intent(inout) :: hotHalo
	  double precision                      , intent(in   ) :: massRateLimited
	</interface>
	<callWith>hotHalo,massRateLimited</callWith>
      </eventHook>
      !!]
    end if
    return
  end subroutine Node_Component_Hot_Halo_Standard_Push_From_Halo

  double precision function Node_Component_Hot_Halo_Standard_Outflow_Stripped_Fraction(node,hotHalo)
    !!{
    Compute the fraction of material outflowing into the hot halo of {\normalfont \ttfamily node} which is susceptible to being stripped
    away.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHaloStandard, treeNode
    implicit none
    type            (treeNode                    ), intent(inout), pointer :: node
    class           (nodeComponentHotHaloStandard)                         :: hotHalo
    double precision                                                       :: massOuter  , massVirial  , &
         &                                                                    radiusOuter, radiusVirial

    radiusOuter =hotHalo                 %outerRadius (                 )
    radiusVirial=darkMatterHaloScale_    %virialRadius(node             )
    massOuter   =hotHaloMassDistribution_%enclosedMass(node,radiusOuter )
    massVirial  =hotHaloMassDistribution_%enclosedMass(node,radiusVirial)
    if (massVirial > 0.0d0) then
       Node_Component_Hot_Halo_Standard_Outflow_Stripped_Fraction=hotHaloOutflowStrippingEfficiency*(1.0d0-massOuter/massVirial)
    else
       Node_Component_Hot_Halo_Standard_Outflow_Stripped_Fraction=hotHaloOutflowStrippingEfficiency
    end if
    return
  end function Node_Component_Hot_Halo_Standard_Outflow_Stripped_Fraction

  subroutine Node_Component_Hot_Halo_Standard_Outflowing_Mass_Rate(self,rate,interrupt,interruptProcedure)
    !!{
    Accept outflowing gas from a galaxy and deposit it into the outflowed and stripped reservoirs.
    !!}
    use :: Galacticus_Nodes, only : interruptTask, nodeComponentHotHalo, nodeComponentHotHaloStandard, treeNode
    implicit none
    class           (nodeComponentHotHalo), intent(inout)                    :: self
    double precision                      , intent(in   )                    :: rate
    logical                               , intent(inout), optional          :: interrupt
    procedure       (interruptTask       ), intent(inout), optional, pointer :: interruptProcedure
    type            (treeNode            )                         , pointer :: node
    double precision                                                         :: strippedOutflowFraction
    !$GLC attributes unused :: interrupt, interruptProcedure

    select type (self)
    class is (nodeComponentHotHaloStandard)
       ! Get the host node.
       node => self%host()
       if (node%isSatellite().and.hotHaloTrackStrippedGas) then
          strippedOutflowFraction=Node_Component_Hot_Halo_Standard_Outflow_Stripped_Fraction(node,self)
          call self% strippedMassRate(rate*       strippedOutflowFraction )
       else
          strippedOutflowFraction=0.0d0
       end if
       ! Funnel the outflow gas into the outflowed and stripped reservoirs in the computed proportions.
       call    self%outflowedMassRate(rate*(1.0d0-strippedOutflowFraction))
    end select
    return
  end subroutine Node_Component_Hot_Halo_Standard_Outflowing_Mass_Rate

  subroutine Node_Component_Hot_Halo_Standard_Outflowing_Ang_Mom_Rate(self,rate,interrupt,interruptProcedure)
    !!{
    Accept outflowing gas angular momentum from a galaxy and deposit it into the outflowed reservoir.
    !!}
    use :: Galacticus_Nodes                     , only : interruptTask                     , nodeComponentHotHalo, nodeComponentHotHaloStandard, treeNode
    use :: Node_Component_Hot_Halo_Standard_Data, only : hotHaloAngularMomentumLossFraction
    implicit none
    class           (nodeComponentHotHalo), intent(inout)                    :: self
    double precision                      , intent(in   )                    :: rate
    logical                               , intent(inout), optional          :: interrupt
    procedure       (interruptTask       ), intent(inout), optional, pointer :: interruptProcedure
    type            (treeNode            )                         , pointer :: node
    double precision                                                         :: strippedOutflowFraction
    !$GLC attributes unused :: interrupt, interruptProcedure

    select type (self)
    class is (nodeComponentHotHaloStandard)
       ! Get the host node.
       node => self%host()
       if (node%isSatellite().and.hotHaloTrackStrippedGas) then
          strippedOutflowFraction=Node_Component_Hot_Halo_Standard_Outflow_Stripped_Fraction(node,self)
       else
          strippedOutflowFraction=0.0d0
       end if
       call    self%outflowedAngularMomentumRate(rate*(1.0d0-strippedOutflowFraction)/(1.0d0-hotHaloAngularMomentumLossFraction))
    end select
    return
  end subroutine Node_Component_Hot_Halo_Standard_Outflowing_Ang_Mom_Rate

  subroutine Node_Component_Hot_Halo_Standard_Outflowing_Abundances_Rate(self,rate,interrupt,interruptProcedure)
    !!{
    Accept outflowing gas abundances from a galaxy and deposit it into the outflowed reservoir.
    !!}
    use :: Abundances_Structure, only : abundances
    use :: Galacticus_Nodes    , only : interruptTask, nodeComponentHotHalo, nodeComponentHotHaloStandard, treeNode
    implicit none
    class           (nodeComponentHotHalo), intent(inout)                    :: self
    type            (abundances          ), intent(in   )                    :: rate
    logical                               , intent(inout), optional          :: interrupt
    procedure       (interruptTask       ), intent(inout), optional, pointer :: interruptProcedure
    type            (treeNode            )                         , pointer :: node
    double precision                                                         :: strippedOutflowFraction
    !$GLC attributes unused :: interrupt, interruptProcedure

    select type (self)
    class is (nodeComponentHotHaloStandard)
       ! Get the host node.
       node => self%host()
       if (node%isSatellite().and.hotHaloTrackStrippedGas) then
          strippedOutflowFraction=Node_Component_Hot_Halo_Standard_Outflow_Stripped_Fraction(node,self)
          call    self%strippedAbundancesRate(rate*       strippedOutflowFraction )
       else
          strippedOutflowFraction=0.0d0
       end if
       call    self%outflowedAbundancesRate(rate*(1.0d0-strippedOutflowFraction))
    end select
    return
  end subroutine Node_Component_Hot_Halo_Standard_Outflowing_Abundances_Rate

  !![
  <preEvolveTask>
  <unitName>Node_Component_Hot_Halo_Standard_Pre_Evolve</unitName>
  </preEvolveTask>
  !!]
  subroutine Node_Component_Hot_Halo_Standard_Pre_Evolve(node)
    !!{
    Ensure the standard hot halo has been initialized.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHalo, nodeComponentHotHaloStandard, treeNode, defaultHotHaloComponent
    implicit none
    type (treeNode            ), intent(inout), pointer :: node
    class(nodeComponentHotHalo)               , pointer :: hotHalo

    ! Check if we are the default method.
    if (.not.defaultHotHaloComponent%standardIsActive()) return
    ! Get the spheroid component.
    hotHalo => node%hotHalo()
    ! Check if a standard hot halo component exists.
    select type (hotHalo)
    class is (nodeComponentHotHaloStandard)
       ! Initialize the hot halo.
       call Node_Component_Hot_Halo_Standard_Initializor(hotHalo)
    end select
    return
  end subroutine Node_Component_Hot_Halo_Standard_Pre_Evolve

  !![
  <rateComputeTask>
   <unitName>Node_Component_Hot_Halo_Standard_Rate_Compute</unitName>
  </rateComputeTask>
  !!]
  subroutine Node_Component_Hot_Halo_Standard_Rate_Compute(node,interrupt,interruptProcedure,propertyType)
    !!{
    Compute the hot halo node mass rate of change.
    !!}
    use :: Abundances_Structure                 , only : abs
    use :: Accretion_Halos                      , only : accretionModeHot                             , accretionModeTotal
    use :: Chemical_Abundances_Structure        , only : chemicalAbundances
    use :: Chemical_Reaction_Rates_Utilities    , only : Chemicals_Mass_To_Density_Conversion
    use :: Dark_Matter_Halo_Spins               , only : Dark_Matter_Halo_Angular_Momentum_Growth_Rate
    use :: Galacticus_Nodes                     , only : defaultHotHaloComponent                      , interruptTask       , nodeComponentBasic, nodeComponentHotHalo, &
          &                                              nodeComponentHotHaloStandard                 , propertyTypeInactive, treeNode
    use :: Node_Component_Hot_Halo_Standard_Data, only : outerRadiusOverVirialRadiusMinimum
    use :: Numerical_Constants_Astronomical     , only : gigaYear
    use :: Numerical_Constants_Math             , only : Pi
    use :: Radiation_Fields                     , only : radiationFieldIntergalacticBackground
    implicit none
    type            (treeNode              )           , intent(inout)          :: node
    logical                                            , intent(inout)          :: interrupt
    procedure       (interruptTask         )           , intent(inout), pointer :: interruptProcedure
    integer                                            , intent(in   )          :: propertyType
    class           (nodeComponentHotHalo  )                          , pointer :: hotHalo
    class           (nodeComponentBasic    )                          , pointer :: basic
    type            (chemicalAbundances    ), save                              :: chemicalDensitiesRates      , chemicalMasses         , &
         &                                                                         chemicalMassesRates         , chemicalsCoolingRate   , &
         &                                                                         chemicalDensities
    !$omp threadprivate(chemicalMasses,chemicalDensities,chemicalDensitiesRates,chemicalMassesRates,chemicalsCoolingRate)
    double precision                                                            :: angularMomentumAccretionRate, temperature            , &
         &                                                                         densityAtOuterRadius        , failedMassAccretionRate, &
         &                                                                         massLossRate                , massToDensityConversion, &
         &                                                                         outerRadius                 , outerRadiusGrowthRate  , &
         &                                                                         massAccretionRate

    ! Return immediately if inactive variables are requested.
    if (propertyType == propertyTypeInactive) return
    ! Return immediately if this class is not in use.
    if (.not.defaultHotHaloComponent%standardIsActive()) return
    ! Reset calculations if necessary.
    if (node%uniqueID() /= uniqueIDPrevious) call Node_Component_Hot_Halo_Standard_Reset(node)
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Ensure that the standard hot halo implementation is active.
    if (defaultHotHaloComponent%standardIsActive()) then
       ! Find the rate of gas mass accretion onto the halo. We take all of the unaccreted gas, including any which is nominally
       ! cold mode, since we do not care what mode it should be in (since it is not actually accreted). It will be assigned to the
       ! relevant mode later when reaccreted. Negative accretion is allowed here as (for example) we may need to record the amount
       ! of negative accretion during periods where a halo is declining in mass, such that this loss of mass must be accounted for
       ! when the halo starts growing again before we allow gas to actually accrete into the hot component.
       massAccretionRate      =accretionHalo_%accretionRate      (node,accretionModeHot  )
       failedMassAccretionRate=accretionHalo_%failedAccretionRate(node,accretionModeTotal)
       ! Get the basic component.
       basic => node%basic()
       ! Apply accretion rates.
       call hotHalo%          massRate(      massAccretionRate,interrupt,interruptProcedure)
       call hotHalo%unaccretedMassRate(failedMassAccretionRate,interrupt,interruptProcedure)
       ! Next compute the cooling rate in this halo.
       call Node_Component_Hot_Halo_Standard_Cooling_Rate(node)
       ! Pipe the cooling rate to which ever component claimed it.
       call Node_Component_Hot_Halo_Standard_Push_To_Cooling_Pipes(node,rateCooling,interrupt,interruptProcedure)
       ! Get the rate at which abundances are accreted onto this halo.
       call hotHalo%abundancesRate(accretionHalo_%accretionRateMetals(node,accretionModeHot),interrupt,interruptProcedure)
       ! Next block of tasks occur only if the accretion rate is non-zero.
       if (basic%accretionRate() /= 0.0d0) then
          ! Compute the rate of accretion of angular momentum.
          angularMomentumAccretionRate=Dark_Matter_Halo_Angular_Momentum_Growth_Rate(node,darkMatterProfileDMO_)*(massAccretionRate &
               &/basic%accretionRate())
             if (hotHaloAngularMomentumAlwaysGrows) angularMomentumAccretionRate=abs(angularMomentumAccretionRate)
          call hotHalo%angularMomentumRate(angularMomentumAccretionRate,interrupt,interruptProcedure)
       end if
       ! Next block of tasks occur only if chemicals are being tracked.
       if (chemicalsCount > 0) then
          ! Get the rate at which chemicals are accreted onto this halo.
          call hotHalo%chemicalsRate(accretionHalo_%accretionRateChemicals(node,accretionModeHot),interrupt,interruptProcedure)
          ! For non-zero cooling rate, adjust the rates of chemical masses.
          if (rateCooling > 0.0d0) then
             ! Compute the rate at which chemicals are lost via cooling.
             chemicalsCoolingRate=hotHalo%chemicals()*rateCooling/hotHalo%mass()
             ! Adjust the rates of chemical masses accordingly.
             call hotHalo%chemicalsRate(-chemicalsCoolingRate,interrupt,interruptProcedure)
          end if
          ! Compute the rates of change of chemical masses due to chemical reactions.
          ! Get the temperature of the hot reservoir.
          temperature=darkMatterHaloScale_%virialTemperature(node)
          ! Set the radiation background.
          call radiationCosmicMicrowaveBackground%timeSet(basic%time())
          if (associated(radiationIntergalacticBackground)) then
             select type (radiationIntergalacticBackground)
             class is (radiationFieldIntergalacticBackground)
                call radiationIntergalacticBackground%timeSet(basic%time())
             end select
          end if
          ! Get the masses of chemicals.
          chemicalMasses=hotHalo%chemicals()
          ! Truncate masses to zero to avoid unphysical behavior.
          call chemicalMasses%enforcePositive()
          ! Scale all chemical masses by their mass in atomic mass units to get a number density.
          call chemicalMasses%massToNumber(chemicalDensities)
          ! Compute factor converting mass of chemicals in (M_Solar/M_Atomic) to number density in cm^-3.
          massToDensityConversion=Chemicals_Mass_To_Density_Conversion(darkMatterHaloScale_%virialRadius(node))
          ! Convert to number density.
          chemicalDensities=chemicalDensities*massToDensityConversion
          ! Compute the chemical reaction rates.
          call chemicalReactionRate_%rates(temperature,chemicalDensities,radiation,chemicalDensitiesRates,node)
          ! Convert to mass change rates.
          call chemicalDensitiesRates%numberToMass(chemicalMassesRates)
          chemicalMassesRates=chemicalMassesRates*gigaYear/massToDensityConversion
          ! Adjust rates appropriately.
          call hotHalo%chemicalsRate(chemicalMassesRates,interrupt,interruptProcedure)
       end if
       ! Perform return of outflowed material.
       select type (hotHalo)
       class is (nodeComponentHotHaloStandard)
          call hotHalo%outflowReturn(interrupt,interruptProcedure)
          ! Test whether this halo is a satellite or not.
          if (node%isSatellite()) then
             ! For satellites, get the current ram pressure stripping radius for this hot halo.
             outerRadiusGrowthRate=hotHalo%outerRadiusGrowthRate()
             outerRadius          =hotHalo%outerRadius          ()
             if     (                                                                                                      &
                  &   outerRadiusGrowthRate   /= 0.0d0                                                                     &
                  &  .and.                                                                                                 &
                  &   hotHalo%mass         () >  0.0d0                                                                     &
                  &  .and.                                                                                                 &
                  &   outerRadius             <=                                   darkMatterHaloScale_%virialRadius(node) &
                  &  .and.                                                                                                 &
                  &   outerRadius             > outerRadiusOverVirialRadiusMinimum*darkMatterHaloScale_%virialRadius(node) &
                  & ) then
                densityAtOuterRadius=hotHaloMassDistribution_%density(node,outerRadius)
                massLossRate        =4.0d0*Pi*densityAtOuterRadius*outerRadius**2*outerRadiusGrowthRate
                call hotHalo%outerRadiusRate(+outerRadiusGrowthRate,interrupt,interruptProcedure)
                call hotHalo%   massSinkRate(+         massLossRate,interrupt,interruptProcedure)
                call Node_Component_Hot_Halo_Standard_Strip_Gas_Rate(node,-massLossRate,interrupt,interruptProcedure)
             end if
          else
             ! For isolated halos, the outer radius should grow with the virial radius.
             call hotHalo%outerRadiusRate(darkMatterHaloScale_%virialRadiusGrowthRate(node),interrupt,interruptProcedure)
          end if
       end select
    end if
    return
  end subroutine Node_Component_Hot_Halo_Standard_Rate_Compute

  double precision function Node_Component_Hot_Halo_Standard_Outer_Radius_Growth_Rate(self)
    !!{
    Compute the growth rate of the outer radius of the hot halo.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHaloStandard, treeNode
    implicit none
    class           (nodeComponentHotHaloStandard), intent(inout) :: self
    type            (treeNode                    ), pointer       :: node
    double precision                                              :: ramPressureRadius, outerRadius

    ! Compute the outer radius growth rate if necessary.
    if (.not.gotOuterRadiusGrowthRate) then
       node              => self                        %hostNode
       ramPressureRadius =  hotHaloRamPressureStripping_%radiusStripped(node)
       outerRadius       =  self                        %outerRadius   (    )
       ! Test whether the ram pressure radius is smaller than the current outer radius of the hot gas profile.
       if     (                                           &
            &  ramPressureRadius      < outerRadius .and. &
            &  self%angularMomentum() >       0.0d0       &
            & ) then
          ! The ram pressure stripping radius is within the outer radius. Cause the outer radius to shrink to the ram pressure
          ! stripping radius on the halo dynamical timescale.
          outerRadiusGrowthRateStored=                                &
               &  (ramPressureRadius-outerRadius)                     &
               & /hotHaloRamPressureTimescale_%timescale(node)
       else
          outerRadiusGrowthRateStored=0.0d0
       end if
       ! Record that outer radius growth rate is now computed.
       gotOuterRadiusGrowthRate=.true.
    end if
    ! Return the pre-computed value.
    Node_Component_Hot_Halo_Standard_Outer_Radius_Growth_Rate=outerRadiusGrowthRateStored
    return
  end function Node_Component_Hot_Halo_Standard_Outer_Radius_Growth_Rate

  subroutine Node_Component_Hot_Halo_Standard_Outflow_Return(self,interrupt,interruptProcedure)
    !!{
    Return outflowed gas to the hot halo.
    !!}
    use :: Abundances_Structure                 , only : abundances                           , max
    use :: Chemical_Abundances_Structure        , only : chemicalAbundances
    use :: Chemical_Reaction_Rates_Utilities    , only : Chemicals_Mass_To_Density_Conversion
    use :: Galacticus_Error                     , only : Galacticus_Error_Report
    use :: Galacticus_Nodes                     , only : interruptTask                        , nodeComponentBasic       , nodeComponentHotHaloStandard, treeNode
    use :: Node_Component_Hot_Halo_Standard_Data, only : starveSatellites                     , starveSatellitesOutflowed
    use :: Numerical_Constants_Atomic           , only : atomicMassHydrogen
    use :: Numerical_Constants_Math             , only : Pi
    use :: Radiation_Fields                     , only : radiationFieldIntergalacticBackground
    implicit none
    class           (nodeComponentHotHaloStandard), intent(inout)          :: self
    logical                                       , intent(inout)          :: interrupt
    procedure       (interruptTask               ), intent(inout), pointer :: interruptProcedure
    type            (treeNode                    )               , pointer :: node
    class           (nodeComponentBasic          )               , pointer :: basic
    double precision                                                       :: outflowedMass            , massReturnRate         , &
         &                                                                    angularMomentumReturnRate, massToDensityConversion, &
         &                                                                    hydrogenByMass           , temperature            , &
         &                                                                    numberDensityHydrogen    , densityAtOuterRadius   , &
         &                                                                    radiusVirial             , outerRadius            , &
         &                                                                    densityMinimum
    type            (abundances                  ), save                   :: abundancesReturnRate     , outflowedAbundances
    !$omp threadprivate(abundancesReturnRate,outflowedAbundances)
    type            (chemicalAbundances          ), save                   :: chemicalDensities        , chemicalMasses         , &
         &                                                                    chemicalMassesRates
    !$omp threadprivate(chemicalDensities,chemicalMassesRates,chemicalMasses)

    ! Get the hosting node.
    node => self%hostNode
    ! Next tasks occur only for systems in which outflowed gas is being recycled.
    if (.not.(starveSatellites.or.starveSatellitesOutflowed).or..not.node%isSatellite()) then
       outflowedMass =self                          %outflowedMass(    )
       massReturnRate=hotHaloOutflowReincorporation_%rate         (node)
       call self%              outflowedMassRate(-           massReturnRate,interrupt,interruptProcedure)
       call self%                       massRate(+           massReturnRate,interrupt,interruptProcedure)
       if (outflowedMass /= 0.0d0) then
          angularMomentumReturnRate=self%outflowedAngularMomentum()*massReturnRate/outflowedMass
          abundancesReturnRate     =self%outflowedAbundances     ()*massReturnRate/outflowedMass
          call self%outflowedAngularMomentumRate(-angularMomentumReturnRate,interrupt,interruptProcedure)
          call self%         angularMomentumRate(+angularMomentumReturnRate,interrupt,interruptProcedure)
          call self%     outflowedAbundancesRate(-     abundancesReturnRate,interrupt,interruptProcedure)
          call self%              abundancesRate(+     abundancesReturnRate,interrupt,interruptProcedure)
       end if
       ! The outer radius must be increased as the halo fills up with gas.
       outerRadius =self                %outerRadius (    )
       radiusVirial=darkMatterHaloScale_%virialRadius(node)
       if (outerRadius < radiusVirial) then
          basic                => node                    %basic  (                )
          densityAtOuterRadius =  hotHaloMassDistribution_%density(node,outerRadius)
          ! If the outer radius and density are non-zero we can expand the outer radius at a rate determined by the current
          ! density profile.
          if (outerRadius > 0.0d0 .and. densityAtOuterRadius > 0.0d0) then
             ! Limit the density at the outer radius to one third of the mean virial density (for baryons, assuming a
             ! universal baryon fraction) to prevent arbitrarily rapid growth of the outer radius in halos containing almost
             ! no gas.
             densityMinimum=(cosmologyParameters_%omegaBaryon()/cosmologyParameters_%omegaMatter())*basic%mass()/radiusVirial**3/4.0d0/Pi
             call self%outerRadiusRate(                           &
                  &                     massReturnRate            &
                  &                    /4.0d0                     &
                  &                    /Pi                        &
                  &                    /outerRadius**2            &
                  &                    /max(                      &
                  &                         densityAtOuterRadius, &
                  &                         densityMinimum        &
                  &                        )                      &
                  &                   )
             ! Otherwise, if we have a positive rate of mass return, simply grow the radius at the virial velocity.
          else if (massReturnRate > 0.0d0) then
             ! Force some growth here so the radius is not trapped at zero.
             call self%outerRadiusRate(                       &
                  &                    +massReturnRate        &
                  &                    /basic         %mass() &
                  &                    *radiusVirial          &
                  &                   )
          end if
       end if
       ! If we have a non-zero return rate, compute associated chemical rates.
       if (chemicalsCount > 0 .and. massReturnRate /= 0.0d0) then
          ! Get required objects.
          basic => node%basic()
          ! Compute coefficient in conversion of mass to density for this node.
          massToDensityConversion=Chemicals_Mass_To_Density_Conversion(darkMatterHaloScale_%virialRadius(node))/3.0d0
          ! Get the abundances of the outflowed material.
          outflowedAbundances    =self%outflowedAbundances()/outflowedMass
          ! Get the hydrogen mass fraction in outflowed gas.
          hydrogenByMass         =outflowedAbundances%hydrogenMassFraction()
          ! Compute the temperature and density of material in the hot halo.
          temperature            =darkMatterHaloScale_%virialTemperature(node)
          numberDensityHydrogen  =hydrogenByMass*outflowedMass*massToDensityConversion/atomicMassHydrogen
          ! Set the radiation field.
          call radiationCosmicMicrowaveBackground%timeSet(basic%time())
          if (associated(radiationIntergalacticBackground)) then
             select type (radiationIntergalacticBackground)
             class is (radiationFieldIntergalacticBackground)
                call radiationIntergalacticBackground%timeSet(basic%time())
             end select
          end if
          ! Get the chemical densities.
          call chemicalState_%chemicalDensities(chemicalDensities,numberDensityHydrogen,temperature,outflowedAbundances,radiation)
          ! Convert from densities to masses.
          call chemicalDensities%numberToMass(chemicalMasses)
          chemicalMassesRates=chemicalMasses*massReturnRate*hydrogenByMass/numberDensityHydrogen/atomicMassHydrogen
          ! Compute the rate at which chemicals are returned to the hot reservoir.
          call self%chemicalsRate(chemicalMassesRates,interrupt,interruptProcedure)
       end if
    end if
    return
  end subroutine Node_Component_Hot_Halo_Standard_Outflow_Return

  subroutine Node_Component_Hot_Halo_Standard_Mass_Sink(self,setValue,interrupt,interruptProcedure)
    !!{
    Account for a sink of gaseous material in the standard hot halo hot gas.
    !!}
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: Galacticus_Nodes, only : interruptTask          , nodeComponentHotHalo, nodeComponentHotHaloStandard
    implicit none
    class           (nodeComponentHotHalo), intent(inout)                    :: self
    double precision                      , intent(in   )                    :: setValue
    logical                               , intent(inout), optional          :: interrupt
    procedure       (interruptTask       ), intent(inout), optional, pointer :: interruptProcedure

    select type (self)
    class is (nodeComponentHotHaloStandard)
       ! Trap cases where an attempt is made to add gas via this sink function.
       if (setValue > 0.0d0) call Galacticus_Error_Report('attempt to add mass via sink in hot halo'//{introspection:location})
       ! Proportionally adjust the rates of all components of the hot gas reservoir.
       call Node_Component_Hot_Halo_Standard_Hot_Gas_All_Rate(self,setValue,interrupt,interruptProcedure)
    end select
    return
  end subroutine Node_Component_Hot_Halo_Standard_Mass_Sink

  subroutine Node_Component_Hot_Halo_Standard_Hot_Gas_All_Rate(self,gasMassRate,interrupt,interruptProcedure)
    !!{
    Adjusts the rates of all components of the hot gas reservoir under the assumption of uniformly distributed properties
    (e.g. fully-mixed metals).
    !!}
    use :: Galacticus_Nodes, only : interruptTask, nodeComponentHotHaloStandard
    implicit none
    class           (nodeComponentHotHaloStandard), intent(inout)                    :: self
    double precision                              , intent(in   )                    :: gasMassRate
    logical                                       , intent(inout), optional          :: interrupt
    procedure       (interruptTask               ), intent(inout), optional, pointer :: interruptProcedure
    double precision                                                                 :: gasMass

    ! Exit immediately for zero rate.
    if (gasMassRate == 0.0d0) return
    ! Get the gas mass present.
    gasMass=self%mass()
    ! If gas is present, adjust the rates.
    if (gasMass > 0.0d0) then
       ! Mass.
       call self%           massRate(                       gasMassRate        ,interrupt,interruptProcedure)
       ! Angular momentum.
       call self%angularMomentumRate(self%angularMomentum()*gasMassRate/gasMass,interrupt,interruptProcedure)
       ! Metal abundances.
       call self%     abundancesRate(self%abundances     ()*gasMassRate/gasMass,interrupt,interruptProcedure)
       ! Chemical abundances.
       call self%      chemicalsRate(self%chemicals      ()*gasMassRate/gasMass,interrupt,interruptProcedure)
    end if
    return
  end subroutine Node_Component_Hot_Halo_Standard_Hot_Gas_All_Rate

  !![
  <scaleSetTask>
   <unitName>Node_Component_Hot_Halo_Standard_Scale_Set</unitName>
  </scaleSetTask>
  !!]
  subroutine Node_Component_Hot_Halo_Standard_Scale_Set(node)
    !!{
    Set scales for properties of {\normalfont \ttfamily node}.
    !!}
    use :: Abundances_Structure         , only : unitAbundances
    use :: Chemical_Abundances_Structure, only : unitChemicalAbundances
    use :: Galacticus_Nodes             , only : nodeComponentBasic     , nodeComponentHotHalo, nodeComponentHotHaloStandard, treeNode, &
         &                                       defaultHotHaloComponent
    implicit none
    type            (treeNode            ), intent(inout), pointer :: node
    class           (nodeComponentHotHalo)               , pointer :: hotHalo
    class           (nodeComponentBasic  )               , pointer :: basic
    double precision                                               :: massVirial    , radiusVirial, &
         &                                                            velocityVirial

    ! Check if we are the default method.
    if (.not.defaultHotHaloComponent%standardIsActive()) return
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Ensure that it is of the standard class.
    select type (hotHalo)
    class is (nodeComponentHotHaloStandard)
       ! The the basic component.
       basic => node%basic()
       ! Get virial properties.
       massVirial    =basic%mass()
       radiusVirial  =darkMatterHaloScale_%virialRadius  (node)
       velocityVirial=darkMatterHaloScale_%virialVelocity(node)
       call    hotHalo%                    massScale(                       massVirial                               *scaleMassRelative  )
       call    hotHalo%           outflowedMassScale(                       massVirial                               *scaleMassRelative  )
       call    hotHalo%          unaccretedMassScale(                       massVirial                               *scaleMassRelative  )
       call    hotHalo%              abundancesScale(unitAbundances        *massVirial                               *scaleMassRelative  )
       call    hotHalo%     outflowedAbundancesScale(unitAbundances        *massVirial                               *scaleMassRelative  )
       call    hotHalo%               chemicalsScale(unitChemicalAbundances*massVirial                               *scaleMassRelative  )
       call    hotHalo%         angularMomentumScale(                       massVirial*radiusVirial*velocityVirial   *scaleMassRelative  )
       call    hotHalo%outflowedAngularMomentumScale(                       massVirial*radiusVirial*velocityVirial   *scaleMassRelative  )
       call    hotHalo%             outerRadiusScale(                                  radiusVirial                  *scaleRadiusRelative)
       if (hotHaloTrackStrippedGas) then
          call hotHalo%            strippedMassScale(                       massVirial                               *scaleMassRelative  )
          call hotHalo%      strippedAbundancesScale(unitAbundances        *massVirial                               *scaleMassRelative  )
       end if
    end select
    return
  end subroutine Node_Component_Hot_Halo_Standard_Scale_Set

  !![
  <mergerTreeInitializeTask>
   <unitName>Node_Component_Hot_Halo_Standard_Tree_Initialize</unitName>
   <after>darkMatterProfile</after>
  </mergerTreeInitializeTask>
  !!]
  subroutine Node_Component_Hot_Halo_Standard_Tree_Initialize(node)
    !!{
    Initialize the contents of the hot halo component for any sub-resolution accretion (i.e. the gas that would have been
    accreted if the merger tree had infinite resolution).
    !!}
    use :: Accretion_Halos       , only : accretionModeHot                 , accretionModeTotal
    use :: Dark_Matter_Halo_Spins, only : Dark_Matter_Halo_Angular_Momentum
    use :: Galacticus_Nodes      , only : defaultHotHaloComponent          , nodeComponentBasic, nodeComponentHotHalo, nodeEvent, &
          &                               nodeEventSubhaloPromotion        , treeNode
    implicit none
    type            (treeNode              ), intent(inout), pointer :: node
    class           (nodeComponentHotHalo  )               , pointer :: currentHotHaloComponent, hotHalo
    class           (nodeComponentBasic    )               , pointer :: basic
    class           (nodeEvent             )               , pointer :: event
    double precision                                                 :: angularMomentum        , failedMass, &
         &                                                              hotHaloMass

    ! If the node has a child or the standard hot halo is not active, then return immediately.
    if (associated(node%firstChild).or..not.defaultHotHaloComponent%standardIsActive()) return

    ! Search for a subhalo promotion events associated with this node.
    event => node%event
    do while (associated(event))
       ! Check if this event:
       !  a) is a subhalo promotion event;
       !  b) has no associated task (which means this is the node being promoted to, not the node being promoted itself).
       ! Do not assign any mass to such nodes, as they should receive gas from the node which is promoted to them.
       select type (event)
       type is (nodeEventSubhaloPromotion)
          if (.not.associated(event%task)) return
       end select
       event => event%next
    end do

    ! Get the hot halo component.
    currentHotHaloComponent => node%hotHalo()
    ! Ensure that it is of unspecified class.
    select type (currentHotHaloComponent)
    type is (nodeComponentHotHalo)
      ! Get the mass of hot gas accreted and the mass that failed to accrete.
       hotHaloMass=accretionHalo_%accretedMass      (node,accretionModeHot  )
       failedMass =accretionHalo_%failedAccretedMass(node,accretionModeTotal)
       ! If either is non-zero, then create a hot halo component and add these masses to it.
       if (hotHaloMass > 0.0d0 .or. failedMass > 0.0d0) then
          call Node_Component_Hot_Halo_Standard_Create(node)
          hotHalo => node%hotHalo()
          basic   => node%basic  ()
          call hotHalo%           massSet(hotHaloMass)
          call hotHalo% unaccretedMassSet( failedMass)
          ! Also add the appropriate angular momentum.
          angularMomentum=hotHaloMass*Dark_Matter_Halo_Angular_Momentum(node,darkMatterProfileDMO_)/basic%mass()
          call hotHalo%angularMomentumSet(angularMomentum   )
          ! Add the appropriate abundances.
          call hotHalo%abundancesSet(accretionHalo_%accretedMassMetals   (node,accretionModeHot))
          ! Also add the appropriate chemical masses.
          call hotHalo%chemicalsSet (accretionHalo_%accretedMassChemicals(node,accretionModeHot))
       end if
    end select
    return
  end subroutine Node_Component_Hot_Halo_Standard_Tree_Initialize

  !![
  <nodeMergerTask>
   <unitName>Node_Component_Hot_Halo_Standard_Node_Merger</unitName>
  </nodeMergerTask>
  !!]
  subroutine Node_Component_Hot_Halo_Standard_Node_Merger(node)
    !!{
    Starve {\normalfont \ttfamily node} by transferring its hot halo to its parent.
    !!}
    use :: Abundances_Structure                 , only : abundances                          , operator(*)            , zeroAbundances              , operator(>)
    use :: Accretion_Halos                      , only : accretionModeHot                    , accretionModeTotal
    use :: Chemical_Abundances_Structure        , only : chemicalAbundances                  , operator(*)            , zeroChemicalAbundances      , operator(>)
    use :: Galactic_Structure_Enclosed_Masses   , only : Galactic_Structure_Enclosed_Mass
    use :: Galactic_Structure_Options           , only : componentTypeAll                    , massTypeBaryonic       , radiusLarge
    use :: Galacticus_Nodes                     , only : nodeComponentBasic                  , nodeComponentHotHalo   , nodeComponentHotHaloStandard, nodeComponentSpin, &
          &                                              treeNode                            , defaultHotHaloComponent
    use :: Galacticus_Error                     , only : Galacticus_Error_Report 
    use :: Node_Component_Hot_Halo_Standard_Data, only : hotHaloNodeMergerLimitBaryonFraction, starveSatellites       , starveSatellitesOutflowed
    implicit none
    type            (treeNode            ), intent(inout) :: node
    type            (treeNode            ), pointer       :: nodeParent
    class           (nodeComponentHotHalo), pointer       :: hotHaloParent          , hotHalo
    class           (nodeComponentSpin   ), pointer       :: spinParent
    class           (nodeComponentBasic  ), pointer       :: parentBasic            , basic
    double precision                                      :: baryonicMassCurrent    , baryonicMassMaximum      , &
         &                                                   fractionRemove         , massAccreted             , &
         &                                                   massUnaccreted         , massReaccreted           , &
         &                                                   fractionAccreted       , angularMomentumAccreted  , &
         &                                                   massAccretedHot
    logical                                               :: massTotalNonZero
    type            (abundances          ), save          :: massMetalsAccreted     , fractionMetalsAccreted   , &
         &                                                   massMetalsReaccreted
    !$omp threadprivate(massMetalsAccreted,fractionMetalsAccreted,massMetalsReaccreted)
    type            (chemicalAbundances  ), save          :: massChemicalsAccreted  , fractionChemicalsAccreted, &
         &                                                   massChemicalsReaccreted
    !$omp threadprivate(massChemicalsAccreted,fractionChemicalsAccreted,massChemicalsReaccreted)

    ! Return immediately if this class is not in use.
    if (.not.defaultHotHaloComponent%standardIsActive()) return
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Ensure that it is of unspecified class.
    select type (hotHalo)
    class is (nodeComponentHotHaloStandard)
       ! Find the parent node and its hot halo and spin components.
       nodeParent             => node  %parent
       call Node_Component_Hot_Halo_Standard_Create(nodeParent)
       hotHaloParent => nodeParent%hotHalo(autoCreate=.true.)
       spinParent    => nodeParent%spin   (                 )
       parentBasic   => nodeParent%basic  (                 )
       basic         => node      %basic  (                 )
       ! Any gas that failed to be accreted by this halo is always transferred to the parent.
       call hotHaloParent%unaccretedMassSet(                                &
            &                                hotHaloParent%unaccretedMass() &
            &                               +hotHalo      %unaccretedMass() &
            &                              )
       call hotHalo      %unaccretedMassSet(                                &
            &                                0.0d0                          &
            &                              )
       ! Since the parent node is undergoing mass growth through this merger we potentially return some of the unaccreted gas to
       ! the hot phase.
       !! First, find the masses of hot and failed mass the node would have if it formed instantaneously.
       massAccretedHot =accretionHalo_%      accretedMass(nodeParent,accretionModeHot  )
       massAccreted    =accretionHalo_%      accretedMass(nodeParent,accretionModeTotal)
       massUnaccreted  =accretionHalo_%failedAccretedMass(nodeParent,accretionModeTotal)
       massTotalNonZero=+massAccreted+massUnaccreted > 0.0d0
       !! Find the fraction of mass that would be successfully accreted.
       if (massAccretedHot > 0.0d0) then
          if (.not.massTotalNonZero) call Galacticus_Error_Report('mass of hot-mode gas accreted is non-zero, but total mass is zero'//{introspection:location})
          fractionAccreted=+  massAccretedHot &
               &           /(                 &
               &             +massAccreted    & 
               &             +massUnaccreted  &
               &            )
          !! Find the change in the unaccreted mass.
          massReaccreted=+hotHaloParent   %unaccretedMass() &
               &         *fractionAccreted                  &
               &         *basic           %          mass() &
               &         /parentBasic     %          mass() 
          !! Reaccrete the gas.
          call hotHaloParent%unaccretedMassSet(hotHaloParent%unaccretedMass()-massReaccreted)
          call hotHaloParent%          massSet(hotHaloParent%          mass()+massReaccreted)
          ! Compute the reaccreted angular momentum.
          angularMomentumAccreted=+massreaccreted                                  &
               &                  *spinParent          %spin          (          ) &
               &                  *darkMatterHaloScale_%virialRadius  (nodeParent) &
               &                  *darkMatterHaloScale_%virialVelocity(nodeParent)
          call hotHaloParent%angularMomentumSet(hotHaloParent%angularMomentum()+angularMomentumAccreted)
       end if
       ! Compute the reaccreted metals.
       !! First, find the metal mass the node would have if it formed instantaneously.
       massMetalsAccreted=accretionHalo_%accretedMassMetals(nodeParent,accretionModeHot)
       !! Find the mass fraction of metals that would be successfully accreted.
       if (massMetalsAccreted > zeroAbundances) then
          if (.not.massTotalNonZero) call Galacticus_Error_Report('mass of hot-mode metals accreted is non-zero, but total mass is zero'//{introspection:location})
          fractionMetalsAccreted=+  massMetalsAccreted &
               &                 /(                    &
               &                   +massAccreted       &
               &                   +massUnaccreted     &
               &                  )
          !! Find the change in the unaccreted mass.
          massMetalsReaccreted=+hotHaloParent   %unaccretedMass() &
               &               *fractionMetalsAccreted            &
               &               *basic           %          mass() &
               &               /parentBasic     %          mass()
          !! Reaccrete the metals.
          call hotHaloParent%abundancesSet(hotHaloParent%abundances()+massMetalsReaccreted)
       end if
       ! Compute the reaccreted chemicals.
       !! First, find the chemicals mass that would be successfully accreted.
       massChemicalsAccreted=accretionHalo_%accretedMassChemicals(nodeParent,accretionModeHot)
       !! Find the mass fraction of chemicals that would be successfully accreted.
       if (massChemicalsAccreted > zeroChemicalAbundances) then
          if (.not.massTotalNonZero) call Galacticus_Error_Report('mass of hot-mode chemicals accreted is non-zero, but total mass is zero'//{introspection:location})
          fractionChemicalsAccreted=+  massChemicalsAccreted &
               &                    /(                       &
               &                      +massAccreted          &
               &                      +massUnaccreted        &
               &                     )
          !! Find the change in the unaccreted mass.
          massChemicalsReaccreted=+hotHaloParent   %unaccretedMass() &
               &                  *fractionChemicalsAccreted         &
               &                  *basic           %          mass() &
               &                  /parentBasic     %          mass()
          !! Reaccrete the metals.
          call hotHaloParent%chemicalsSet(hotHaloParent%chemicals()+massChemicalsReaccreted)
       end if
       ! Determine if starvation is to be applied.
       if (starveSatellites.or.starveSatellitesOutflowed) then
          ! Move the hot halo to the parent. We leave the hot halo in place even if it is starved, since outflows will accumulate to
          ! this hot halo (and will be moved to the parent at the end of the evolution timestep).
          if (starveSatellites) then
             call hotHaloParent%                    massSet(                                                           &
                  &                                          hotHaloParent       %mass                    (          ) &
                  &                                         +hotHalo             %mass                    (          ) &
                  &                                        )
             call hotHaloParent%         angularMomentumSet(                                                           &
                  &                                          hotHaloParent       %angularMomentum         (          ) &
                  &                                         +hotHalo             %mass                    (          ) &
                  &                                         *spinParent          %spin                    (          ) &
                  &                                         *darkMatterHaloScale_%virialRadius            (nodeParent) &
                  &                                         *darkMatterHaloScale_%virialVelocity          (nodeParent) &
                  &                                        )
          end if
          call    hotHaloParent%           outflowedMassSet(                                                           &
               &                                             hotHaloParent       %outflowedMass           (          ) &
               &                                            +hotHalo             %outflowedMass           (          ) &
               &                                           )
          call    hotHaloParent%outflowedAngularMomentumSet(                                                           &
               &                                             hotHaloParent       %outflowedAngularMomentum(          ) &
               &                                            +hotHalo             %outflowedMass           (          ) &
               &                                            *   spinParent       %spin                    (          ) &
               &                                            *darkMatterHaloScale_%virialRadius            (nodeParent) &
               &                                            *darkMatterHaloScale_%virialVelocity          (nodeParent) &
               &                                           )
          if (starveSatellites) then
             call hotHalo      %                    massSet(                                     &
                  &                                          0.0d0                               &
                  &                                        )
             call hotHalo      %         angularMomentumSet(                                     &
                  &                                          0.0d0                               &
                  &                                        )
          end if
          call    hotHalo      %           outflowedMassSet(                                     &
               &                                             0.0d0                               &
               &                                           )
          call    hotHalo      %outflowedAngularMomentumSet(                                     &
               &                                             0.0d0                               &
               &                                           )
          if (starveSatellites) then
             call hotHaloParent%              abundancesSet(                                     &
                  &                                          hotHaloParent%abundances         () &
                  &                                         +hotHalo      %abundances         () &
                  &                                        )
             call hotHalo      %              abundancesSet(                                     &
                  &                                          zeroAbundances                      &
                  &                                        )
          end if
          call    hotHaloParent%     outflowedAbundancesSet(                                     &
               &                                             hotHaloParent%outflowedAbundances() &
               &                                            +hotHalo      %outflowedAbundances() &
               &                                           )
          call    hotHalo      %     outflowedAbundancesSet(                                     &
               &                                             zeroAbundances                      &
               &                                           )
          if (starveSatellites) then
             call hotHaloParent%               chemicalsSet(                                     &
                  &                                          hotHaloParent%chemicals          () &
                  &                                         +hotHalo      %chemicals          () &
                  &                                        )
             call hotHalo      %               chemicalsSet(                                     &
                  &                                          zeroChemicalAbundances              &
                  &                                        )
          end if
          ! Check if the baryon fraction in the parent hot halo exceeds the universal value. If it does, mitigate this by moving
          ! some of the mass to the failed accretion reservoir.
          if (hotHaloNodeMergerLimitBaryonFraction) then
             ! Get the default cosmology.
             baryonicMassMaximum=parentBasic%mass()*cosmologyParameters_%OmegaBaryon()/cosmologyParameters_%OmegaMatter()
             baryonicMassCurrent=Galactic_Structure_Enclosed_Mass(nodeParent,radiusLarge,massType=massTypeBaryonic&
                  &,componentType =componentTypeAll)
             if (baryonicMassCurrent > baryonicMassMaximum .and. hotHaloParent%mass() > 0.0d0) then
                fractionRemove=min((baryonicMassCurrent-baryonicMassMaximum)/hotHaloParent%massTotal(),1.0d0)
                call hotHaloParent% unaccretedMassSet(                                                       &
                     &                                 hotHaloParent%unaccretedMass ()                       &
                     &                                +hotHaloParent%mass           ()*       fractionRemove &
                     &                               )
                call hotHaloParent%           massSet( hotHaloParent%mass           ()*(1.0d0-fractionRemove))
                call hotHaloParent%angularMomentumSet( hotHaloParent%angularMomentum()*(1.0d0-fractionRemove))
                call hotHaloParent%     abundancesSet( hotHaloParent%abundances     ()*(1.0d0-fractionRemove))
                call hotHaloParent%     chemicalsSet ( hotHaloParent%chemicals      ()*(1.0d0-fractionRemove))
             end if
          end if
       end if
    end select
    return
  end subroutine Node_Component_Hot_Halo_Standard_Node_Merger

  subroutine satelliteMerger(self,node)
    !!{
    Remove any hot halo associated with {\normalfont \ttfamily node} before it merges with its host halo.
    !!}
    use :: Abundances_Structure                 , only : abundances            , zeroAbundances
    use :: Chemical_Abundances_Structure        , only : zeroChemicalAbundances
    use :: Galacticus_Nodes                     , only : nodeComponentHotHalo  , nodeComponentHotHaloStandard, nodeComponentSpin, treeNode
    use :: Node_Component_Hot_Halo_Standard_Data, only : starveSatellites
    implicit none
    class(*                   ), intent(inout) :: self
    type (treeNode            ), intent(inout) :: node
    type (treeNode            ), pointer       :: nodeHost
    class(nodeComponentHotHalo), pointer       :: hotHaloHost, hotHalo
    class(nodeComponentSpin   ), pointer       :: spinHost
    !$GLC attributes unused :: self

    ! Return immediately if satellites are starved, as in that case there is no hot halo to transfer.
    if (starveSatellites) return
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Ensure that it is of unspecified class.
    select type (hotHalo)
    class is (nodeComponentHotHaloStandard)
       ! Find the node to merge with.
       nodeHost             => node    %mergesWith(                 )
       hotHaloHost          => nodeHost%hotHalo   (autoCreate=.true.)
       spinHost             => nodeHost%spin      (                 )
       ! Move the hot halo to the host.
       call hotHaloHost%                    massSet(                                                         &
            &                                        hotHaloHost         %mass                    (        ) &
            &                                       +hotHalo             %mass                    (        ) &
            &                                      )
       call hotHaloHost%          unaccretedMassSet(                                                         &
            &                                        hotHaloHost         %unaccretedMass          (        ) &
            &                                       +hotHalo             %unaccretedMass          (        ) &
            &                                      )
       call hotHaloHost%         angularMomentumSet(                                                         &
            &                                        hotHaloHost         %angularMomentum         (        ) &
            &                                       +hotHalo             %mass                    (        ) &
            &                                       *spinHost            %spin                    (        ) &
            &                                       *darkMatterHaloScale_%virialRadius            (nodeHost) &
            &                                       *darkMatterHaloScale_%virialVelocity          (nodeHost) &
            &                                      )
       call hotHaloHost%           outflowedMassSet(                                                         &
            &                                        hotHaloHost         %outflowedMass           (        ) &
            &                                       +hotHalo             %outflowedMass           (        ) &
            &                                      )
       call hotHaloHost%outflowedAngularMomentumSet(                                                         &
            &                                        hotHaloHost         %outflowedAngularMomentum(        ) &
            &                                       +hotHalo             %outflowedMass           (        ) &
            &                                       *spinHost            %spin                    (        ) &
            &                                       *darkMatterHaloScale_%virialRadius            (nodeHost) &
            &                                       *darkMatterHaloScale_%virialVelocity          (nodeHost) &
            &                                      )
       call hotHalo    %                    massSet(                                                         &
            &                                        0.0d0                                                   &
            &                                      )
       call hotHalo    %         angularMomentumSet(                                                         &
            &                                        0.0d0                                                   &
            &                                      )
       call hotHalo    %           outflowedMassSet(                                                         &
            &                                        0.0d0                                                   &
            &                                      )
       call hotHalo    %outflowedAngularMomentumSet(                                                         &
            &                                        0.0d0                                                   &
            &                                      )
       call hotHaloHost%              abundancesSet(                                                         &
            &                                        hotHaloHost         %abundances              (        ) &
            &                                       +hotHalo             %abundances              (        ) &
            &                                      )
       call hotHalo    %              abundancesSet(                                                         &
            &                                        zeroAbundances                                          &
            &                                      )
       call hotHaloHost%     outflowedAbundancesSet(                                                         &
            &                                        hotHaloHost         %outflowedAbundances     (        ) &
            &                                       +hotHalo             %outflowedAbundances     (        ) &
            &                                      )
       call hotHalo    %     outflowedAbundancesSet(                                                         &
            &                                        zeroAbundances                                          &
            &                                      )
       call hotHaloHost%               chemicalsSet(                                                         &
            &                                        hotHaloHost         %chemicals               (        ) &
            &                                       +hotHalo             %chemicals               (        ) &
            &                                      )
       call hotHalo    %               chemicalsSet(                                                         &
            &                                        zeroChemicalAbundances                                  &
            &                                      )
    end select
    return
  end subroutine satelliteMerger

  subroutine nodePromotion(self,node)
    !!{
    Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply update the hot halo mass of {\normalfont \ttfamily
    node} to account for any hot halo already in the parent.
    !!}
    use :: Abundances_Structure         , only : abundances            , zeroAbundances
    use :: Chemical_Abundances_Structure, only : zeroChemicalAbundances
    use :: Galacticus_Nodes             , only : nodeComponentHotHalo  , nodeComponentHotHaloStandard, treeNode
    implicit none
    class(*                   ), intent(inout) :: self
    type (treeNode            ), intent(inout) :: node
    type (treeNode            ), pointer       :: nodeParent
    class(nodeComponentHotHalo), pointer       :: hotHaloParent, hotHalo
    !$GLC attributes unused :: self
    
    hotHalo => node%hotHalo(autoCreate=.true.)
    ! Ensure that it is of specified class.
    select type (hotHalo)
    class is (nodeComponentHotHaloStandard)
       nodeParent    => node      %parent
       hotHaloParent => nodeParent%hotHalo(autoCreate=.true.)
       ! Update the outer radius to match the virial radius of the parent halo.
       call hotHalo%outerRadiusSet(                                              &
            &                      darkMatterHaloScale_%virialRadius(nodeParent) &
            &                     )
       ! If the parent node has a hot halo component, then add it to that of this node, and perform other changes needed prior to
       ! promotion.
       select type (hotHaloParent)
          class is (nodeComponentHotHaloStandard)
             ! If (outflowed) mass is non-positive, set mass and all related quantities to zero.
          if (hotHalo%         mass() <= 0.0d0) then
             call hotHalo%massSet           (0.0d0                 )
             call hotHalo%angularMomentumSet(0.0d0                 )
             call hotHalo%abundancesSet     (zeroAbundances        )
             call hotHalo%chemicalsSet      (zeroChemicalAbundances)
          end if
          if (hotHalo%outflowedMass() <= 0.0d0) then
             call hotHalo%outflowedMassSet           (         0.0d0)
             call hotHalo%outflowedAngularMomentumSet(         0.0d0)
             call hotHalo%outflowedAbundancesSet     (zeroAbundances)
          end if
          call    hotHalo%                    massSet(                                          &
               &                                       hotHalo      %mass                    () &
               &                                      +hotHaloParent%mass                    () &
               &                                     )
          call    hotHalo%         angularMomentumSet(                                          &
               &                                       hotHalo      %angularMomentum         () &
               &                                      +hotHaloParent%angularMomentum         () &
               &                                     )
          call    hotHalo%           outflowedMassSet(                                          &
               &                                       hotHalo      %outflowedMass           () &
               &                                      +hotHaloParent%outflowedMass           () &
               &                                     )
          call    hotHalo%outflowedAngularMomentumSet(                                          &
               &                                       hotHalo      %outflowedAngularMomentum() &
               &                                      +hotHaloParent%outflowedAngularMomentum() &
               &                                     )
          call    hotHalo%              abundancesSet(                                          &
               &                                       hotHalo      %abundances              () &
               &                                      +hotHaloParent%abundances              () &
               &                                     )
          call    hotHalo%     outflowedAbundancesSet(                                          &
               &                                       hotHalo      %outflowedAbundances     () &
               &                                      +hotHaloParent%outflowedAbundances     () &
               &                                     )
          call    hotHalo%               chemicalsSet(                                          &
               &                                       hotHalo      %chemicals               () &
               &                                      +hotHaloParent%chemicals               () &
               &                                     )
          call    hotHalo%          unaccretedMassSet(                                          &
               &                                       hotHalo      %unaccretedMass          () &
               &                                      +hotHaloParent%unaccretedMass          () &
               &                                     )
       end select
    end select
    return
  end subroutine nodePromotion

  subroutine Node_Component_Hot_Halo_Standard_Cooling_Rate(node)
    !!{
    Get and store the cooling rate for {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHalo, treeNode
    implicit none
    type (treeNode            ), intent(inout) :: node
    class(nodeComponentHotHalo), pointer       :: hotHalo

    if (.not.gotCoolingRate) then
       ! Get the hot halo component.
       hotHalo => node%hotHalo()
       if     (                                   &
            &   hotHalo%mass           () > 0.0d0 &
            &  .and.                              &
            &   hotHalo%angularMomentum() > 0.0d0 &
            &  .and.                              &
            &   hotHalo%outerRadius    () > 0.0d0 &
            & ) then
          ! Get the cooling rate.
          rateCooling=coolingRate_%rate(node)
       else
          rateCooling=0.0d0
       end if

       ! Store a copy of this cooling rate as the remaining mass heating rate budget. This is used to ensure that we never heat
       ! gas at a rate greater than it is cooling.
       massHeatingRateRemaining=rateCooling

       ! Flag that cooling rate has now been computed.
       gotCoolingRate=.true.
    end if
    return
  end subroutine Node_Component_Hot_Halo_Standard_Cooling_Rate

  subroutine Node_Component_Hot_Halo_Standard_Create(node)
    !!{
    Creates a hot halo component for {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHalo, nodeComponentHotHaloStandard, treeNode
    implicit none
    type (treeNode            ), intent(inout), pointer :: node
    class(nodeComponentHotHalo)               , pointer :: hotHalo

    ! Create the component.
    hotHalo => node%hotHalo(autoCreate=.true.)
    ! Initalize.
    select type (hotHalo)
    class is (nodeComponentHotHaloStandard)
       call Node_Component_Hot_Halo_Standard_Initializor(hotHalo)
    end select
    return
  end subroutine Node_Component_Hot_Halo_Standard_Create

  subroutine Node_Component_Hot_Halo_Standard_Initializor(self)
    !!{
    Initializes a standard hot halo component.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHaloStandard, treeNode
    implicit none
    type (nodeComponentHotHaloStandard)          :: self
    type (treeNode                    ), pointer :: node

    ! Return if already initialized.
    if (self%isInitialized()) return
    ! Get the hosting node.
    node => self%hostNode
    ! Initialize the outer boundary to the virial radius.
    call self%outerRadiusSet(darkMatterHaloScale_%virialRadius(node))
    ! Record that the spheroid has been initialized.
    call self%isInitializedSet(.true.)
    return
  end subroutine Node_Component_Hot_Halo_Standard_Initializor

  !![
  <haloFormationTask>
   <unitName>Node_Component_Hot_Halo_Standard_Formation</unitName>
  </haloFormationTask>
  !!]
  subroutine Node_Component_Hot_Halo_Standard_Formation(node)
    !!{
    Updates the hot halo gas distribution at a formation event, if requested.
    !!}
    use :: Abundances_Structure                 , only : abundances                           , zeroAbundances
    use :: Chemical_Abundances_Structure        , only : chemicalAbundances
    use :: Chemical_Reaction_Rates_Utilities    , only : Chemicals_Mass_To_Density_Conversion
    use :: Galacticus_Nodes                     , only : nodeComponentBasic                   , nodeComponentHotHalo, nodeComponentHotHaloStandard, treeNode
    use :: Node_Component_Hot_Halo_Standard_Data, only : hotHaloOutflowReturnOnFormation
    use :: Numerical_Constants_Atomic           , only : atomicMassHydrogen
    use :: Radiation_Fields                     , only : radiationFieldIntergalacticBackground
    implicit none
    type            (treeNode            ), intent(inout) :: node
    class           (nodeComponentBasic  ), pointer       :: basic
    class           (nodeComponentHotHalo), pointer       :: hotHalo
    type            (abundances          ), save          :: outflowedAbundances
    type            (chemicalAbundances  ), save          :: chemicalDensities    , chemicalMasses
    !$omp threadprivate(outflowedAbundances,chemicalDensities,chemicalMasses)
    double precision                                      :: hydrogenByMass       , massToDensityConversion, &
         &                                                   numberDensityHydrogen, temperature

    ! Return immediately if return of outflowed gas on formation events is not requested.
    if (.not.hotHaloOutflowReturnOnFormation) return

    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Ensure that it is of unspecified class.
    select type (hotHalo)
    class is (nodeComponentHotHaloStandard)

       ! Compute mass of chemicals transferred to the hot halo.
       if (chemicalsCount > 0 .and. hotHalo%outflowedMass() > 0.0d0) then
          ! Compute coefficient in conversion of mass to density for this node.
          massToDensityConversion=Chemicals_Mass_To_Density_Conversion(darkMatterHaloScale_%virialRadius(node))/3.0d0
          ! Get abundance mass fractions of the outflowed material.
          outflowedAbundances= hotHalo%outflowedAbundances() &
               &              /hotHalo%outflowedMass      ()
          ! Get the hydrogen mass fraction in outflowed gas.
          hydrogenByMass=outflowedAbundances%hydrogenMassFraction()
          ! Compute the temperature and density of material in the hot halo.
          temperature          =darkMatterHaloScale_%virialTemperature(node)
          numberDensityHydrogen=hydrogenByMass*hotHalo%outflowedMass()*massToDensityConversion/atomicMassHydrogen
          ! Set the radiation field.
          basic => node%basic()
          call radiationCosmicMicrowaveBackground%timeSet(basic%time())
          if (associated(radiationIntergalacticBackground)) then
             select type (radiationIntergalacticBackground)
             class is (radiationFieldIntergalacticBackground)
                call radiationIntergalacticBackground%timeSet(basic%time())
             end select
          end if
           ! Get the chemical densities.
          call chemicalState_%chemicalDensities(chemicalDensities,numberDensityHydrogen,temperature,outflowedAbundances,radiation)
          ! Convert from densities to masses.
          call chemicalDensities%numberToMass(chemicalMasses)
          chemicalMasses=chemicalMasses*hotHalo%outflowedMass()*hydrogenByMass/numberDensityHydrogen/atomicMassHydrogen
          ! Add chemicals to the hot component.
          call hotHalo%chemicalsSet(                      &
               &                      hotHalo%chemicals() &
               &                     +chemicalMasses      &
               &                    )
       end if
       ! Transfer mass, angular momentum and abundances.
       call hotHalo%                    massSet(&
            &                                    hotHalo%         mass           () &
            &                                   +hotHalo%outflowedMass           () &
            &                                  )
       call hotHalo%         angularMomentumSet(                                    &
            &                                    hotHalo%         angularMomentum() &
            &                                   +hotHalo%outflowedAngularMomentum() &
            &                                  )
       call hotHalo%              abundancesSet(                                    &
            &                                    hotHalo%         abundances     () &
            &                                   +hotHalo%outflowedAbundances     () &
            &                                  )
       call hotHalo%           outflowedMassSet(                                    &
            &                                    0.0d0                              &
            &                                  )
       call hotHalo%outflowedAngularMomentumSet(                                    &
            &                                    0.0d0                              &
            &                                  )
       call hotHalo%     outflowedAbundancesSet(                                    &
            &                                    zeroAbundances                     &
            &                                  )
    end select

    return
  end subroutine Node_Component_Hot_Halo_Standard_Formation

end module Node_Component_Hot_Halo_Standard
