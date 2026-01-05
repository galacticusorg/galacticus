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
Contains a module which implements the standard hot halo node component.
!!}

module Node_Component_Hot_Halo_Standard
  !!{
  Implements the standard hot halo node component.
  !!}
  use :: Accretion_Halos                           , only : accretionHaloClass
  use :: Chemical_Reaction_Rates                   , only : chemicalReactionRateClass
  use :: Chemical_States                           , only : chemicalStateClass
  use :: Cosmology_Functions                       , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters                      , only : cosmologyParametersClass
  use :: Cooling_Infall_Torques                    , only : coolingInfallTorqueClass
  use :: Dark_Matter_Halo_Scales                   , only : darkMatterHaloScaleClass
  use :: Hot_Halo_Mass_Distributions               , only : hotHaloMassDistributionClass
  use :: Hot_Halo_Temperature_Profiles             , only : hotHaloTemperatureProfileClass
  use :: Hot_Halo_Outflows_Reincorporations        , only : hotHaloOutflowReincorporationClass
  use :: Hot_Halo_Outflows_Stripping               , only : hotHaloOutflowStrippingClass
  use :: Hot_Halo_Ram_Pressure_Stripping           , only : hotHaloRamPressureStrippingClass
  use :: Hot_Halo_Ram_Pressure_Stripping_Timescales, only : hotHaloRamPressureTimescaleClass
  use :: Kind_Numbers                              , only : kind_int8
  use :: Radiation_Fields                          , only : radiationFieldClass
  implicit none
  private
  public :: Node_Component_Hot_Halo_Standard_Initialize         , Node_Component_Hot_Halo_Standard_Thread_Initialize, &
       &    Node_Component_Hot_Halo_Standard_Scale_Set          , Node_Component_Hot_Halo_Standard_Tree_Initialize  , &
       &    Node_Component_Hot_Halo_Standard_Post_Step          , Node_Component_Hot_Halo_Standard_Reset            , &
       &    Node_Component_Hot_Halo_Standard_Rate_Compute       , Node_Component_Hot_Halo_Standard_Pre_Evolve       , &
       &    Node_Component_Hot_Halo_Standard_State_Store        , Node_Component_Hot_Halo_Standard_State_Restore    , &
       &    Node_Component_Hot_Halo_Standard_Thread_Uninitialize

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
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" isNonNegative="true" />
      <output unitsInSI="massSolar" comment="Mass of gas in the hot halo."/>
    </property>
    <property>
      <name>abundances</name>
      <type>abundances</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" isNonNegative="true" />
      <output unitsInSI="massSolar" comment="Mass of metals in the hot phase of the hot halo."/>
    </property>
    <property>
      <name>chemicals</name>
      <type>chemicalAbundances</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" isNonNegative="true" />
      <output unitsInSI="massSolar" comment="Mass of chemicals in the hot phase of the hot halo."/>
    </property>
    <property>
      <name>angularMomentum</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" isNonNegative="true" />
      <output unitsInSI="massSolar*megaParsec*kilo" comment="Angular momentum of gas in the hot halo."/>
    </property>
    <property>
      <name>outflowedMass</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" isNonNegative="true" />
      <output unitsInSI="massSolar" comment="Mass of outflowed gas in the hot halo."/>
    </property>
    <property>
      <name>outflowedAngularMomentum</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" isNonNegative="true" />
      <output unitsInSI="massSolar*megaParsec*kilo" comment="Angular momentum of outflowed gas in the hot halo."/>
    </property>
    <property>
      <name>outflowedAbundances</name>
      <type>abundances</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" isNonNegative="true" />
      <output unitsInSI="massSolar" comment="Mass of metals in the outflowed phase of the hot halo."/>
    </property>
    <property>
      <name>outflowedChemicals</name>
      <type>chemicalAbundances</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" isNonNegative="true" />
      <output unitsInSI="massSolar" comment="Mass of chemicals in the outflowed phase of the hot halo."/>
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
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" isNonNegative="true" />
      <output unitsInSI="massSolar" comment="Mass of gas that failed to accrete into the hot halo."/>
    </property>
    <property>
      <name>unaccretedAbundances</name>
      <type>abundances</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" isNonNegative="true" />
      <output unitsInSI="massSolar" comment="Mass of metals that failed to accrete into the hot halo."/>
    </property>
    <property>
      <name>outerRadius</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" isDeferred="get" isNonNegative="true" />
      <output unitsInSI="megaParsec" comment="Outer radius of the hot halo."/>
    </property>
    <property>
      <name>strippedMass</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" isNonNegative="true" />
    </property>
    <property>
      <name>strippedAbundances</name>
      <type>abundances</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" isNonNegative="true" />
    </property>
    <property>
      <name>strippedChemicals</name>
      <type>chemicalAbundances</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" isNonNegative="true" />
    </property>
    <property>
      <name>massSink</name>
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
     <binding method="outflowReturn" isDeferred="true" >
      <interface>
       <type>void</type>
       <self pass="true" intent="inout" />
       <argument>logical                 , intent(inout)          :: interrupt</argument>
       <argument>procedure(interruptTask), intent(inout), pointer :: interruptProcedure</argument>
      </interface>
     </binding>
     <binding method="outerRadiusGrowthRate" isDeferred="true" >
      <interface>
       <type>double</type>
       <rank>0</rank>
       <self pass="true" intent="inout" />
      </interface>
     </binding>
     <binding method="massDistribution" isDeferred="true" >
      <interface>
       <type>class(massDistributionClass), pointer</type>
       <rank>0</rank>
       <module>Galactic_Structure_Options, only : enumerationWeightByType, enumerationComponentTypeType, enumerationMassTypeType</module>
       <module>Mass_Distributions        , only : massDistributionClass                                                         </module>
       <self pass="true" intent="inout" />
       <argument>type   (enumerationComponentTypeType), intent(in   ), optional :: componentType</argument>
       <argument>type   (enumerationMassTypeType     ), intent(in   ), optional :: massType     </argument>
       <argument>type   (enumerationWeightByType     ), intent(in   ), optional :: weightBy     </argument>
       <argument>integer                              , intent(in   ), optional :: weightIndex  </argument>
      </interface>
     </binding>
     <binding method="massBaryonic" function="Node_Component_Hot_Halo_Standard_Mass_Baryonic"/>
   </bindings>
   <functions>objects.nodes.components.hot_halo.standard.bound_functions.inc</functions>
  </component>
  !!]

  ! Objects used by this component.
  class(cosmologyFunctionsClass            ), pointer :: cosmologyFunctions_
  class(darkMatterHaloScaleClass           ), pointer :: darkMatterHaloScale_
  class(hotHaloMassDistributionClass       ), pointer :: hotHaloMassDistribution_
  class(hotHaloTemperatureProfileClass     ), pointer :: hotHaloTemperatureProfile_
  class(accretionHaloClass                 ), pointer :: accretionHalo_
  class(coolingInfallTorqueClass           ), pointer :: coolingInfallTorque_
  class(hotHaloRamPressureStrippingClass   ), pointer :: hotHaloRamPressureStripping_
  class(hotHaloRamPressureTimescaleClass   ), pointer :: hotHaloRamPressureTimescale_
  class(hotHaloOutflowReincorporationClass ), pointer :: hotHaloOutflowReincorporation_
  class(hotHaloOutflowStrippingClass       ), pointer :: hotHaloOutflowStripping_
  class(chemicalStateClass                 ), pointer :: chemicalState_
  class(cosmologyParametersClass           ), pointer :: cosmologyParameters_
  !$omp threadprivate(cosmologyFunctions_,darkMatterHaloScale_,hotHaloMassDistribution_,hotHaloTemperatureProfile_,accretionHalo_,chemicalState_,hotHaloRamPressureStripping_,hotHaloRamPressureTimescale_,cosmologyParameters_,hotHaloOutflowReincorporation_,hotHaloOutflowStripping_,coolingInfallTorque_)

  ! Internal count of abundances and chemicals.
  integer                                                                         :: abundancesCount                                             , chemicalsCount

  ! Quantities stored to avoid repeated computation.
  integer         (kind=kind_int8                                    )            :: uniqueIDPrevious
  logical                                                                         :: gotOuterRadiusGrowthRate                            =.false.
  double precision                                                                :: outerRadiusGrowthRateStored
  !$omp threadprivate(gotOuterRadiusGrowthRate,outerRadiusGrowthRateStored,uniqueIDPrevious)
  ! Radiation structure.
  class           (radiationFieldClass                               ), pointer   :: radiation
  !$omp threadprivate(radiation)

  ! Parameters controlling absolute tolerance scales.
  double precision                                                    , parameter :: scaleMassRelative                                  =1.0d-3
  double precision                                                    , parameter :: scaleRadiusRelative                                =1.0d-1

  ! Procedure pointer to mass distribution function.
  procedure       (Node_Component_Hot_Halo_Standard_Mass_Distribution), pointer   :: Node_Component_Hot_Halo_Standard_Mass_Distribution_
  
  ! A threadprivate object used to track to which thread events are attached.
  integer                                                                         :: thread
  !$omp threadprivate(thread)

contains

  !![
  <nodeComponentInitializationTask>
   <unitName>Node_Component_Hot_Halo_Standard_Initialize</unitName>
  </nodeComponentInitializationTask>
  !!]
  subroutine Node_Component_Hot_Halo_Standard_Initialize(parameters)
    !!{
    Initializes the standard hot halo component module.
    !!}
    use :: Abundances_Structure                 , only : Abundances_Property_Count      , abundances
    use :: Chemical_Abundances_Structure        , only : Chemicals_Property_Count
    use :: Error                                , only : Error_Report
    use :: Galacticus_Nodes                     , only : defaultHotHaloComponent        , nodeComponentHotHaloStandard
    use :: ISO_Varying_String                   , only : var_str                        , varying_string              , char
    use :: Input_Parameters                     , only : inputParameter                 , inputParameters
    use :: Node_Component_Hot_Halo_Standard_Data, only : fractionBaryonLimitInNodeMerger, outflowReturnOnFormation    , starveSatellites           , starveSatellitesOutflowed, &
         &                                               angularMomentumAlwaysGrows
    implicit none
    type(inputParameters             ), intent(inout) :: parameters
    type(nodeComponentHotHaloStandard)                :: hotHalo
    type(inputParameters             )                :: subParameters

    ! Initialize the module if necessary.
    !$omp critical (Node_Component_Hot_Halo_Standard_Initialize)
    if (defaultHotHaloComponent%standardIsActive()) then

       ! Get numbers of abundance and chemicals properties.
       abundancesCount=Abundances_Property_Count()
       chemicalsCount =Chemicals_Property_Count ()

       ! Find our parameters.
       subParameters=parameters%subParameters('componentHotHalo')
       ! Determine whether satellite nodes will be starved of gas.
       !![
       <inputParameter>
         <name>starveSatellites</name>
         <defaultValue>.false.</defaultValue>
         <description>Specifies whether or not the hot halo should be removed (``starved'') when a node becomes a satellite.</description>
         <source>subParameters</source>
       </inputParameter>
       !!]

       !![
       <inputParameter>
         <name>starveSatellitesOutflowed</name>
         <defaultValue>.false.</defaultValue>
         <description>Specifies whether or not the outflowed hot halo should be removed (``starved'') when a node becomes a satellite.</description>
         <source>subParameters</source>
       </inputParameter>
       !!]

       ! Determine whether outflowed gas should be restored to the hot reservoir on halo formation events.
       !![
       <inputParameter>
         <name>outflowReturnOnFormation</name>
         <defaultValue>.false.</defaultValue>
         <description>Specifies whether or not outflowed gas should be returned to the hot reservoir on halo formation events.</description>
         <source>subParameters</source>
       </inputParameter>
       !!]

       ! Determine whether negative angular momentum accretion rates onto the halo should be treated as positive for the purposes
       ! of computing the hot halo angular momentum.
       !![
       <inputParameter>
         <name>angularMomentumAlwaysGrows</name>
         <defaultValue>.false.</defaultValue>
         <description>Specifies whether or not negative rates of accretion of angular momentum into the hot halo will be treated as positive
            for the purposes of computing the hot halo angular momentum.</description>
         <source>subParameters</source>
       </inputParameter>
       !!]

       ! Get option controlling limiting of baryon fraction during node mergers.
       !![
       <inputParameter>
         <name>fractionBaryonLimitInNodeMerger</name>
         <defaultValue>.false.</defaultValue>
         <description>Controls whether the hot gas content of nodes should be limited to not exceed the universal baryon fraction at node
           merger events. If set to {\normalfont \ttfamily true}, hot gas (and angular momentum, abundances, and chemicals proportionally) will be
           removed from the merged halo to the unaccreted gas reservoir to limit the baryonic mass to the universal baryon
           fraction where possible.</description>
         <source>subParameters</source>
       </inputParameter>
       !!]
       ! Bind the outer radius get function.
       call hotHalo%                  outerRadiusFunction(Node_Component_Hot_Halo_Standard_Outer_Radius              )
       ! Bind outflowing material pipes to the functions that will handle input of outflowing material to the hot halo.
       call hotHalo%           outflowingMassRateFunction(Node_Component_Hot_Halo_Standard_Outflowing_Mass_Rate      )
       call hotHalo%outflowingAngularMomentumRateFunction(Node_Component_Hot_Halo_Standard_Outflowing_Ang_Mom_Rate   )
       call hotHalo%     outflowingAbundancesRateFunction(Node_Component_Hot_Halo_Standard_Outflowing_Abundances_Rate)
       ! Bind a creation function.
       call hotHalo%                    createFunctionSet(Node_Component_Hot_Halo_Standard_Initializor               )
       ! Bind the mass distribution function.
       Node_Component_Hot_Halo_Standard_Mass_Distribution_ => Node_Component_Hot_Halo_Standard_Mass_Distribution
       call hotHalo%             massDistributionFunction(Node_Component_Hot_Halo_Standard_Mass_Distribution_        )
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
  subroutine Node_Component_Hot_Halo_Standard_Thread_Initialize(parameters)
    !!{
    Initializes the tree node hot halo methods module.
    !!}
    use :: Events_Hooks    , only : haloFormationEvent      , postEvolveEvent, openMPThreadBindingAtLevel, dependencyRegEx, &
         &                          dependencyDirectionAfter
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : defaultHotHaloComponent
    use :: Input_Parameters, only : inputParameter          , inputParameters
    use :: Radiation_Fields, only : radiationFieldNull
    implicit none
    type(inputParameters), intent(inout) :: parameters
    type(dependencyRegEx), dimension(1)  :: dependencies
    type(inputParameters)                :: subParameters

    ! Check if this implementation is selected. Define the radiation component to include both the CMB and the intergalactic background if it is.
    if (defaultHotHaloComponent%standardIsActive()) then
       ! Find our parameters.
       subParameters=parameters%subParameters('componentHotHalo')
       !![
       <objectBuilder class="cosmologyParameters"           name="cosmologyParameters_"           source="subParameters"/>
       <objectBuilder class="cosmologyFunctions"            name="cosmologyFunctions_"            source="subParameters"/>
       <objectBuilder class="darkMatterHaloScale"           name="darkMatterHaloScale_"           source="subParameters"/>
       <objectBuilder class="hotHaloMassDistribution"       name="hotHaloMassDistribution_"       source="subParameters"/>
       <objectBuilder class="hotHaloTemperatureProfile"     name="hotHaloTemperatureProfile_"     source="subParameters"/>
       <objectBuilder class="accretionHalo"                 name="accretionHalo_"                 source="subParameters"/>
       <objectBuilder class="coolingInfallTorque"           name="coolingInfallTorque_"           source="subParameters"/>
       <objectBuilder class="chemicalState"                 name="chemicalState_"                 source="subParameters"/>
       <objectBuilder class="hotHaloRamPressureStripping"   name="hotHaloRamPressureStripping_"   source="subParameters"/>
       <objectBuilder class="hotHaloRamPressureTimescale"   name="hotHaloRamPressureTimescale_"   source="subParameters"/>
       <objectBuilder class="hotHaloOutflowReincorporation" name="hotHaloOutflowReincorporation_" source="subParameters"/>
       <objectBuilder class="hotHaloOutflowStripping"       name="hotHaloOutflowStripping_"       source="subParameters"/>
       !!]
       dependencies(1)=dependencyRegEx(dependencyDirectionAfter,'^remnantStructure:')
       call postEvolveEvent     %attach(thread,postEvolve     ,openMPThreadBindingAtLevel,label='nodeComponentHotHaloStandard'                          )
       call haloFormationEvent  %attach(thread,haloFormation  ,openMPThreadBindingAtLevel,label='nodeComponentHotHaloStandard'                          )
       if (parameters%isPresent('radiationFieldIntergalacticBackground')) then
          !![
          <objectBuilder class="radiationField" name="radiation" parameterName="radiationFieldIntergalacticBackground" source="parameters"/>
          !!]
       else
          allocate(radiationFieldNull :: radiation)
          select type (radiation)
          type is (radiationFieldNull)
             !![
	     <referenceConstruct object="radiation" constructor="radiationFieldNull()"/>
             !!]
          end select
       end if
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
    use :: Events_Hooks    , only : haloFormationEvent     , postEvolveEvent
    use :: Galacticus_Nodes, only : defaultHotHaloComponent
    implicit none

    if (defaultHotHaloComponent%standardIsActive()) then
       !![
       <objectDestructor name="cosmologyParameters_"          />
       <objectDestructor name="cosmologyFunctions_"           />
       <objectDestructor name="darkMatterHaloScale_"          />
       <objectDestructor name="hotHaloMassDistribution_"      />
       <objectDestructor name="hotHaloTemperatureProfile_"    />
       <objectDestructor name="accretionHalo_"                />
       <objectDestructor name="coolingInfallTorque_"          />
       <objectDestructor name="chemicalState_"                />
       <objectDestructor name="hotHaloRamPressureStripping_"  />
       <objectDestructor name="hotHaloRamPressureTimescale_"  />
       <objectDestructor name="hotHaloOutflowReincorporation_"/>
       <objectDestructor name="hotHaloOutflowStripping_"      />
       <objectDestructor name="radiation"                     />
       !!]
       if (postEvolveEvent   %isAttached(thread,postEvolve   )) call postEvolveEvent   %detach(thread,postEvolve   )
       if (haloFormationEvent%isAttached(thread,haloFormation)) call haloFormationEvent%detach(thread,haloFormation)
    end if
    return
  end subroutine Node_Component_Hot_Halo_Standard_Thread_Uninitialize

  !![
  <calculationResetTask>
  <unitName>Node_Component_Hot_Halo_Standard_Reset</unitName>
  </calculationResetTask>
  !!]
  subroutine Node_Component_Hot_Halo_Standard_Reset(node,uniqueID)
    !!{
    Remove memory of stored computed values as we're about to begin computing derivatives anew.
    !!}
    use :: Galacticus_Nodes, only : treeNode
    use :: Kind_Numbers    , only : kind_int8
    implicit none
    type   (treeNode ), intent(inout) :: node
    integer(kind_int8), intent(in   ) :: uniqueID
    !$GLC attributes unused :: node

    uniqueIDPrevious        =uniqueID
    gotOuterRadiusGrowthRate=.false.
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
    radiusVirial =  darkMatterHaloScale_%radiusVirial(selfHost)
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
  </postStepTask>
  !!]
  subroutine Node_Component_Hot_Halo_Standard_Post_Step(node,status)
    !!{
    Do processing of the node required after evolution.
    !!}
    use :: Chemical_Abundances_Structure, only : chemicalAbundances
    use :: Galacticus_Nodes             , only : nodeComponentHotHalo, nodeComponentHotHaloStandard, treeNode   , defaultHotHaloComponent
    use :: Interface_GSL                , only : GSL_Success         , GSL_Continue                , GSL_Failure
    implicit none
    type            (treeNode            ), intent(inout), pointer :: node
    integer                               , intent(inout)          :: status
    class           (nodeComponentHotHalo)               , pointer :: hotHalo
    type            (chemicalAbundances  ), save                   :: chemicalMasses
    !$omp threadprivate(chemicalMasses)
    integer                                                        :: i
    double precision                                               :: massChemicals , massChemicalsPositive

    ! Return immediately if this class is not in use.
    if (.not.defaultHotHaloComponent%standardIsActive()) return
    ! Limit hot gas mass, and outer radius to be non-negative.
    hotHalo => node%hotHalo()
    select type (hotHalo)
    class is (nodeComponentHotHaloStandard)
       ! Note that "status" is not set to failure as these changes in state of the hot halo should not change any calculation of
       ! differential evolution rates as a negative mass/outer radius was unphysical anyway.
       if (hotHalo%       mass() < 0.0d0) then
          call hotHalo%       massSet(0.0d0)
          ! Indicate that ODE evolution should continue after this state change.
          if (status == GSL_Success) status=GSL_Continue
       end if
       if (hotHalo%outerRadius() < 0.0d0) then
          call hotHalo%outerRadiusSet(0.0d0)
          ! Indicate that ODE evolution should continue after this state change.
          if (status == GSL_Success) status=GSL_Continue
       end if
       ! Truncate negative chemical species, keeping total mass of chemicals fixed.
       if (chemicalsCount > 0) then
          massChemicals        =0.0d0
          massChemicalsPositive=0.0d0
          chemicalMasses       =hotHalo%chemicals()
          do i=1,chemicalsCount
             massChemicals=+massChemicals               &
                  &        +chemicalMasses%abundance(i)
             if (chemicalMasses%abundance(i) < 0.0d0) then
                call chemicalMasses%abundanceSet(i,0.0d0)
             else
                massChemicalsPositive=+massChemicalsPositive       &
                     &                +chemicalMasses%abundance(i)
             end if
          end do
          if (massChemicalsPositive > massChemicals) then
             ! Mark ODE failure here to force derivatives to be recomputed.
             status=GSL_Failure
             if (massChemicalsPositive > 0.0d0) then
                call chemicalMasses%scale(                                  &
                     &                    +max(0.0d0,massChemicals        ) &
                     &                    /          massChemicalsPositive  &
                     &                   )
             else
                call chemicalMasses%reset()
             end if
             call hotHalo%chemicalsSet(chemicalMasses)
          end if
          if (chemicalMasses%sumOver() > hotHalo%mass()) then
             ! Ensure total mass of chemicals can not exceed the mass of the hot halo gas.
             call chemicalMasses%scale(                          &
                  &                    +hotHalo       %mass   () &
                  &                    /chemicalMasses%sumOver() &
                  &                   )
             call hotHalo%chemicalsSet(chemicalMasses)
             ! Mark ODE failure here to force derivatives to be recomputed.
             status=GSL_Failure
          end if
       end if
    end select
    return
  end subroutine Node_Component_Hot_Halo_Standard_Post_Step

  subroutine postEvolve(self,node)
    !!{
    Do processing of the node required after evolution.
    !!}
    use :: Abundances_Structure         , only : zeroAbundances
    use :: Chemical_Abundances_Structure, only : zeroChemicalAbundances 
    use :: Galacticus_Nodes             , only : nodeComponentHotHalo  , nodeComponentHotHaloStandard, nodeComponentSpin, nodeComponentBasic, &
         &                                       treeNode
    implicit none
    class(*                   ), intent(inout) :: self
    type (treeNode            ), intent(inout) :: node
    type (treeNode            ), pointer       :: nodeParent
    class(nodeComponentBasic  ), pointer       :: basicParent
    class(nodeComponentHotHalo), pointer       :: hotHaloParent, hotHalo
    class(nodeComponentSpin   ), pointer       :: spinParent
    !$GLC attributes unused :: self

    ! Process hot gas for satellites.
    if (node%isSatellite()) then
       ! Check if stripped mass is being tracked.
       if (.not.hotHaloOutflowStripping_%neverStripped(node)) then
          hotHalo => node%hotHalo()
          select type (hotHalo)
          class is (nodeComponentHotHaloStandard)
             ! Transfer any stripped gas to the host halo.
             nodeParent => node%parent
             do while (nodeParent%isSatellite())
                nodeParent => nodeParent%parent
             end do
             call Node_Component_Hot_Halo_Standard_Create(nodeParent)
             basicParent   => nodeParent%basic  (                 )
             hotHaloParent => nodeParent%hotHalo(autoCreate=.true.)
             spinParent    => nodeParent%spin   (                 )
             call hotHaloParent%outflowedAngularMomentumSet(                                           &
                  &                                          hotHaloParent %outflowedAngularMomentum() &
                  &                                         +hotHalo       %strippedMass            () &
                  &                                         *spinParent    %angularMomentum         () &
                  &                                         /basicParent   %mass                    () &
                  &                                        )
             call hotHaloParent%outflowedMassSet           (                                           &
                  &                                          hotHaloParent %outflowedMass           () &
                  &                                         +hotHalo       %strippedMass            () &
                  &                                        )
             call hotHalo      %strippedMassSet            (                                           &
                  &                                          0.0d0                                     &
                  &                                        )
             call hotHaloParent%outflowedAbundancesSet     (                                           &
                  &                                          hotHaloParent %outflowedAbundances     () &
                  &                                         +hotHalo       %strippedAbundances      () &
                  &                                        )
             call hotHalo      %strippedAbundancesSet      (                                           &
                  &                                          zeroAbundances                            &
                  &                                        )
             call hotHaloParent%outflowedChemicalsSet      (                                           &
                  &                                          hotHaloParent %outflowedChemicals      () &
                  &                                         +hotHalo       %strippedChemicals       () &
                  &                                        )
             call hotHalo      %strippedChemicalsSet       (                                           &
                  &                                          zeroChemicalAbundances                    &
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
          ! Chemical abundances.
          call hotHalo% strippedChemicalsRate(hotHalo%chemicals ()*gasMassRate/gasMass,interrupt,interruptProcedure)
       end if
    end select
    return
  end subroutine Node_Component_Hot_Halo_Standard_Strip_Gas_Rate

  subroutine Node_Component_Hot_Halo_Standard_Outflowing_Mass_Rate(self,rate,interrupt,interruptProcedure)
    !!{
    Accept outflowing gas from a galaxy and deposit it into the outflowed and stripped reservoirs.
    !!}
    use :: Abundances_Structure             , only : abundances
    use :: Chemical_Abundances_Structure    , only : chemicalAbundances
    use :: Chemical_Reaction_Rates_Utilities, only : Chemicals_Mass_To_Density_Conversion
    use :: Galacticus_Nodes                 , only : interruptTask                        , nodeComponentHotHalo, nodeComponentHotHaloStandard, nodeComponentBasic, &
         &                                           treeNode
    use :: Numerical_Constants_Atomic       , only : atomicMassHydrogen
    implicit none
    class           (nodeComponentHotHalo), intent(inout)                    :: self
    double precision                      , intent(in   )                    :: rate
    logical                               , intent(inout), optional          :: interrupt
    procedure       (interruptTask       ), intent(inout), optional, pointer :: interruptProcedure
    type            (treeNode            )                         , pointer :: node
    class           (nodeComponentBasic  )                         , pointer :: basic
    type            (abundances          ), save                             :: outflowedAbundances
    !$omp threadprivate(outflowedAbundances)
    type            (chemicalAbundances  ), save                             :: chemicalDensities      , chemicalMassesRates
    !$omp threadprivate(chemicalDensities,chemicalMassesRates)
    double precision                                                         :: fractionOutflowStripped, massToDensityConversion, &
         &                                                                      hydrogenByMass         , numberDensityHydrogen  , &
         &                                                                      temperature
    !$GLC attributes unused :: interrupt, interruptProcedure

    select type (self)
    class is (nodeComponentHotHaloStandard)
       ! Get the host node.
       node                    => self                    %host            (    )
       fractionOutflowStripped =  hotHaloOutflowStripping_%fractionStripped(node)
       ! Funnel the outflow gas into the outflowed and stripped reservoirs in the computed proportions.
       call self% strippedMassRate(rate*       fractionOutflowStripped )
       call self%outflowedMassRate(rate*(1.0d0-fractionOutflowStripped))
       ! If we have a non-zero return rate, compute associated chemical rates.
       if (chemicalsCount > 0 .and. rate /= 0.0d0 .and. self%outflowedMass() > 0.0d0) then
          ! Get required objects.
          basic => node%basic()
          ! Compute coefficient in conversion of mass to density for this node.
          massToDensityConversion=Chemicals_Mass_To_Density_Conversion(self%outerRadius())/3.0d0
          ! Get the abundances of the outflowed material.
          outflowedAbundances    =self%outflowedAbundances()/self%outflowedMass()
          ! Get the hydrogen mass fraction in outflowed gas.
          hydrogenByMass         =outflowedAbundances%hydrogenMassFraction()
          ! Compute the temperature and density of material in the hot halo.
          temperature            =darkMatterHaloScale_%temperatureVirial(node)
          numberDensityHydrogen  =hydrogenByMass*self%outflowedMass()*massToDensityConversion/atomicMassHydrogen
          ! Set the radiation field.
          call radiation%timeSet(basic%time())
          ! Get the chemical densities.
          call chemicalState_%chemicalDensities(chemicalDensities,numberDensityHydrogen,temperature,outflowedAbundances,radiation)
          ! Convert from densities to mass rates.
          call chemicalDensities  %numberToMass(chemicalMassesRates                                         )
          call chemicalMassesRates%scale       (rate*hydrogenByMass/numberDensityHydrogen/atomicMassHydrogen)
          ! Compute the rate at which chemicals are returned to the hot reservoir.
          if (.not.hotHaloOutflowStripping_%neverStripped(node)) &
               &  call self% strippedChemicalsRate(chemicalMassesRates*       fractionOutflowStripped ,interrupt,interruptProcedure)
          call         self%outflowedChemicalsRate(chemicalMassesRates*(1.0d0-fractionOutflowStripped),interrupt,interruptProcedure)
       end if
    end select
    return
  end subroutine Node_Component_Hot_Halo_Standard_Outflowing_Mass_Rate

  subroutine Node_Component_Hot_Halo_Standard_Outflowing_Ang_Mom_Rate(self,rate,interrupt,interruptProcedure)
    !!{
    Accept outflowing gas angular momentum from a galaxy and deposit it into the outflowed reservoir.
    !!}
    use :: Galacticus_Nodes, only : interruptTask, nodeComponentHotHalo, nodeComponentHotHaloStandard, treeNode
    implicit none
    class           (nodeComponentHotHalo), intent(inout)                    :: self
    double precision                      , intent(in   )                    :: rate
    logical                               , intent(inout), optional          :: interrupt
    procedure       (interruptTask       ), intent(inout), optional, pointer :: interruptProcedure
    type            (treeNode            )                         , pointer :: node
    double precision                                                         :: fractionOutflowStripped, fractionLossAngularMomentum
    !$GLC attributes unused :: interrupt, interruptProcedure

    select type (self)
    class is (nodeComponentHotHaloStandard)
       ! Get the host node.
       node                        => self                    %host                       (    )
       fractionOutflowStripped     =  hotHaloOutflowStripping_%fractionStripped           (node)
       fractionLossAngularMomentum =  coolingInfallTorque_    %fractionAngularMomentumLoss(node)
       call self%outflowedAngularMomentumRate(                                     &
            &                                 +rate                                &
            &                                 *(1.0d0-fractionOutflowStripped    ) &
            &                                 /(1.0d0-fractionLossAngularMomentum) &
            &                                )
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
    double precision                                                         :: fractionOutflowStripped
    !$GLC attributes unused :: interrupt, interruptProcedure

    select type (self)
    class is (nodeComponentHotHaloStandard)
       ! Get the host node.
       node                    => self                    %host            (    )
       fractionOutflowStripped =  hotHaloOutflowStripping_%fractionStripped(node)
       call self%outflowedAbundancesRate(rate*(1.0d0-fractionOutflowStripped))
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
    ! Get the hot halo component.
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
    use :: Abundances_Structure                 , only : abundances                        , abs
    use :: Accretion_Halos                      , only : accretionModeHot                  , accretionModeTotal
    use :: Galacticus_Nodes                     , only : defaultHotHaloComponent           , interruptTask             , nodeComponentHotHalo, &
          &                                              nodeComponentHotHaloStandard      , propertyInactive          , treeNode
    use :: Node_Component_Hot_Halo_Standard_Data, only : outerRadiusOverVirialRadiusMinimum
    use :: Numerical_Constants_Math             , only : Pi
    use :: Mass_Distributions                   , only : massDistributionClass
    use :: Coordinates                          , only : coordinateSpherical               , assignment(=)
    use :: Galactic_Structure_Options           , only : componentTypeHotHalo              , massTypeGaseous
    implicit none
    type            (treeNode             ), intent(inout)          :: node
    logical                                , intent(inout)          :: interrupt
    procedure       (interruptTask        ), intent(inout), pointer :: interruptProcedure
    integer                                , intent(in   )          :: propertyType
    class           (nodeComponentHotHalo )               , pointer :: hotHalo
    class           (massDistributionClass)               , pointer :: massDistribution_
    type            (coordinateSpherical  )                         :: coordinates
    double precision                                                :: outerRadius         , outerRadiusGrowthRate, &
         &                                                             densityAtOuterRadius, massLossRate

    ! Return immediately if inactive variables are requested.
    if (propertyInactive(propertyType)) return
    ! Return immediately if this class is not in use.
    if (.not.defaultHotHaloComponent%standardIsActive()) return
    ! Reset calculations if necessary.
    if (node%uniqueID() /= uniqueIDPrevious) call Node_Component_Hot_Halo_Standard_Reset(node,node%uniqueID())
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Ensure that the standard hot halo implementation is active.
    if (defaultHotHaloComponent%standardIsActive()) then
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
                  &   outerRadius             <=                                   darkMatterHaloScale_%radiusVirial(node) &
                  &  .and.                                                                                                 &
                  &   outerRadius             > outerRadiusOverVirialRadiusMinimum*darkMatterHaloScale_%radiusVirial(node) &
                  & ) then
                coordinates          =  [outerRadius,0.0d0,0.0d0]
                massDistribution_    => node             %massDistribution(componentTypeHotHalo,massTypeGaseous)
                densityAtOuterRadius =  massDistribution_%density         (coordinates                         )
                !![
		<objectDestructor name="massDistribution_"/>
                !!]
                massLossRate        =4.0d0*Pi*densityAtOuterRadius*outerRadius**2*outerRadiusGrowthRate
                call hotHalo%outerRadiusRate(+outerRadiusGrowthRate,interrupt,interruptProcedure)
                call hotHalo%   massSinkRate(+         massLossRate,interrupt,interruptProcedure)
                if (.not.hotHaloOutflowStripping_%neverStripped(node)) &
                     & call Node_Component_Hot_Halo_Standard_Strip_Gas_Rate(node,-massLossRate,interrupt,interruptProcedure)
             end if
          else
             ! For isolated halos, the outer radius should grow with the virial radius.
             call hotHalo%outerRadiusRate(darkMatterHaloScale_%radiusVirialGrowthRate(node),interrupt,interruptProcedure)
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
    use :: Galacticus_Nodes                     , only : interruptTask        , nodeComponentBasic       , nodeComponentHotHaloStandard, treeNode
    use :: Node_Component_Hot_Halo_Standard_Data, only : starveSatellites     , starveSatellitesOutflowed
    use :: Numerical_Constants_Math             , only : Pi
    use :: Mass_Distributions                   , only : massDistributionClass
    use :: Galactic_Structure_Options           , only : componentTypeHotHalo , massTypeGaseous
    use :: Coordinates                          , only : coordinateSpherical  , assignment(=)
    implicit none
    class           (nodeComponentHotHaloStandard), intent(inout)          :: self
    logical                                       , intent(inout)          :: interrupt
    procedure       (interruptTask               ), intent(inout), pointer :: interruptProcedure
    type            (treeNode                    )               , pointer :: node
    class           (nodeComponentBasic          )               , pointer :: basic
    class           (massDistributionClass       )               , pointer :: massDistribution_
    type            (coordinateSpherical         )                         :: coordinates
    double precision                                                       :: densityAtOuterRadius, massReturnRate, &
         &                                                                    radiusVirial        , outerRadius   , &
         &                                                                    densityMinimum
    !$GLC attributes unused :: interrupt, interruptProcedure

    ! Get the hosting node.
    node => self%hostNode
    ! Next tasks occur only for systems in which outflowed gas is being recycled.
    if (.not.(starveSatellites.or.starveSatellitesOutflowed).or..not.node%isSatellite()) then
       massReturnRate=hotHaloOutflowReincorporation_%rate         (node)
       ! The outer radius must be increased as the halo fills up with gas.
       outerRadius =self                %outerRadius (    )
       radiusVirial=darkMatterHaloScale_%radiusVirial(node)
       if (outerRadius < radiusVirial) then
          coordinates          =  [outerRadius,0.0d0,0.0d0]
          basic                => node             %basic           (                                    )
          massDistribution_    => node             %massDistribution(componentTypeHotHalo,massTypeGaseous)
          densityAtOuterRadius =  massDistribution_%density         (coordinates                         )
          !![
	  <objectDestructor name="massDistribution_"/>
          !!]
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
    end if
    return
  end subroutine Node_Component_Hot_Halo_Standard_Outflow_Return

  subroutine Node_Component_Hot_Halo_Standard_Mass_Sink(self,setValue,interrupt,interruptProcedure)
    !!{
    Account for a sink of gaseous material in the standard hot halo hot gas.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : interruptTask, nodeComponentHotHalo, nodeComponentHotHaloStandard
    implicit none
    class           (nodeComponentHotHalo), intent(inout)                    :: self
    double precision                      , intent(in   )                    :: setValue
    logical                               , intent(inout), optional          :: interrupt
    procedure       (interruptTask       ), intent(inout), optional, pointer :: interruptProcedure

    select type (self)
    class is (nodeComponentHotHaloStandard)
       ! Trap cases where an attempt is made to add gas via this sink function.
       if (setValue > 0.0d0) call Error_Report('attempt to add mass via sink in hot halo'//{introspection:location})
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
       call self%           massRate(                        gasMassRate         ,interrupt,interruptProcedure)
       ! Angular momentum.
       call self%angularMomentumRate(self%angularMomentum()*(gasMassRate/gasMass),interrupt,interruptProcedure)
       ! Metal abundances.
       call self%     abundancesRate(self%abundances     ()*(gasMassRate/gasMass),interrupt,interruptProcedure)
       ! Chemical abundances.
       call self%      chemicalsRate(self%chemicals      ()*(gasMassRate/gasMass),interrupt,interruptProcedure)
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
       radiusVirial  =darkMatterHaloScale_%radiusVirial  (node)
       velocityVirial=darkMatterHaloScale_%velocityVirial(node)
       call    hotHalo%                    massScale(                        massVirial                               *scaleMassRelative   )
       call    hotHalo%           outflowedMassScale(                        massVirial                               *scaleMassRelative   )
       call    hotHalo%          unaccretedMassScale(                        massVirial                               *scaleMassRelative   )
       call    hotHalo%              abundancesScale(unitAbundances        *(massVirial                               *scaleMassRelative  ))
       call    hotHalo%    unaccretedAbundancesScale(unitAbundances        *(massVirial                               *scaleMassRelative  ))
       call    hotHalo%     outflowedAbundancesScale(unitAbundances        *(massVirial                               *scaleMassRelative  ))
       call    hotHalo%               chemicalsScale(unitChemicalAbundances*(massVirial                               *scaleMassRelative  ))
       call    hotHalo%      outflowedChemicalsScale(unitChemicalAbundances*(massVirial                               *scaleMassRelative  ))
       call    hotHalo%         angularMomentumScale(                        massVirial*radiusVirial*velocityVirial   *scaleMassRelative   )
       call    hotHalo%outflowedAngularMomentumScale(                        massVirial*radiusVirial*velocityVirial   *scaleMassRelative   )
       call    hotHalo%             outerRadiusScale(                                   radiusVirial                  *scaleRadiusRelative )
       if (.not.hotHaloOutflowStripping_%neverStripped(node)) then
          call hotHalo%            strippedMassScale(                        massVirial                               *scaleMassRelative   )
          call hotHalo%      strippedAbundancesScale(unitAbundances        *(massVirial                               *scaleMassRelative  ))
          call hotHalo%       strippedChemicalsScale(unitChemicalAbundances*(massVirial                               *scaleMassRelative  ))
       end if
    end select
    return
  end subroutine Node_Component_Hot_Halo_Standard_Scale_Set

  !![
  <mergerTreeInitializeTask>
   <unitName>Node_Component_Hot_Halo_Standard_Tree_Initialize</unitName>
  </mergerTreeInitializeTask>
  !!]
  subroutine Node_Component_Hot_Halo_Standard_Tree_Initialize(node)
    !!{
    Initialize the contents of the hot halo component for any sub-resolution accretion (i.e. the gas that would have been
    accreted if the merger tree had infinite resolution).
    !!}
    use :: Accretion_Halos , only : accretionModeHot         , accretionModeTotal
    use :: Galacticus_Nodes, only : defaultHotHaloComponent  , nodeComponentBasic, nodeComponentHotHalo, nodeEvent, &
          &                         nodeEventSubhaloPromotion, treeNode          , nodeComponentSpin
    implicit none
    type            (treeNode            ), intent(inout), pointer :: node
    class           (nodeComponentHotHalo)               , pointer :: currentHotHaloComponent
    class           (nodeEvent           )               , pointer :: event
    double precision                                               :: failedMass             , massHotHalo

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
       massHotHalo=accretionHalo_%accretedMass      (node,accretionModeHot  )
       failedMass =accretionHalo_%failedAccretedMass(node,accretionModeTotal)
       ! If either is non-zero, then create a hot halo component and add these masses to it.
       if (massHotHalo > 0.0d0 .or. failedMass > 0.0d0) then
          ! NOTE: This call is required only to set up the initial outer radius. Once that is initialized by a nodeOperator, this
          !       call (and the remainder of this function) can be removed.
          call Node_Component_Hot_Halo_Standard_Create(node)
       end if
    end select
    return
  end subroutine Node_Component_Hot_Halo_Standard_Tree_Initialize

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
    ! Initialize.
    select type (hotHalo)
    class is (nodeComponentHotHaloStandard)
       call Node_Component_Hot_Halo_Standard_Initializor(hotHalo)
    end select
    return
  end subroutine Node_Component_Hot_Halo_Standard_Create

  subroutine Node_Component_Hot_Halo_Standard_Initializor(self,timeEnd)
    !!{
    Initializes a standard hot halo component.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHaloStandard, treeNode
    implicit none
    type            (nodeComponentHotHaloStandard), intent(inout)           :: self
    double precision                              , intent(in   ), optional :: timeEnd
    type            (treeNode                    ), pointer                 :: node
    !$GLC attributes unused :: timeEnd
    
    ! Return if already initialized.
    if (self%isInitialized()) return
    ! Get the hosting node.
    node => self%hostNode
    ! Initialize the outer boundary to the virial radius.
    call self%outerRadiusSet(darkMatterHaloScale_%radiusVirial(node))
    ! Record that the spheroid has been initialized.
    call self%isInitializedSet(.true.)
    return
  end subroutine Node_Component_Hot_Halo_Standard_Initializor

  subroutine haloFormation(self,node)
    !!{
    Updates the hot halo gas distribution at a formation event, if requested.
    !!}
    use :: Abundances_Structure                 , only : abundances                          , zeroAbundances
    use :: Chemical_Abundances_Structure        , only : chemicalAbundances                  , zeroChemicalAbundances
    use :: Chemical_Reaction_Rates_Utilities    , only : Chemicals_Mass_To_Density_Conversion
    use :: Galacticus_Nodes                     , only : nodeComponentBasic                  , nodeComponentHotHalo  , nodeComponentHotHaloStandard, treeNode
    use :: Node_Component_Hot_Halo_Standard_Data, only : outflowReturnOnFormation
    use :: Numerical_Constants_Atomic           , only : atomicMassHydrogen
    implicit none
    class           (*                   ), intent(inout) :: self
    type            (treeNode            ), intent(inout) :: node
    class           (nodeComponentBasic  ), pointer       :: basic
    class           (nodeComponentHotHalo), pointer       :: hotHalo
    type            (abundances          ), save          :: outflowedAbundances
    type            (chemicalAbundances  ), save          :: chemicalDensities    , chemicalMasses
    !$omp threadprivate(outflowedAbundances,chemicalDensities,chemicalMasses)
    double precision                                      :: hydrogenByMass       , massToDensityConversion, &
         &                                                   numberDensityHydrogen, temperature
    !$GLC attributes unused :: self

    ! Return immediately if return of outflowed gas on formation events is not requested.
    if (.not.outflowReturnOnFormation) return

    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Ensure that it is of unspecified class.
    select type (hotHalo)
    class is (nodeComponentHotHaloStandard)

       ! Compute mass of chemicals transferred to the hot halo.
       if (chemicalsCount > 0 .and. hotHalo%outflowedMass() > 0.0d0) then
          ! Compute coefficient in conversion of mass to density for this node.
          massToDensityConversion=Chemicals_Mass_To_Density_Conversion(hotHalo%outerRadius())/3.0d0
          ! Get abundance mass fractions of the outflowed material.
          outflowedAbundances= hotHalo%outflowedAbundances() &
               &              /hotHalo%outflowedMass      ()
          ! Get the hydrogen mass fraction in outflowed gas.
          hydrogenByMass=outflowedAbundances%hydrogenMassFraction()
          ! Compute the temperature and density of material in the hot halo.
          temperature          =darkMatterHaloScale_%temperatureVirial(node)
          numberDensityHydrogen=hydrogenByMass*hotHalo%outflowedMass()*massToDensityConversion/atomicMassHydrogen
          ! Set the radiation field.
          basic => node%basic()
          call radiation%timeSet(basic%time())
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
       call hotHalo%      outflowedChemicalsSet(                                    &
            &                                    zeroChemicalAbundances             &
            &                                  )
    end select

    return
  end subroutine haloFormation

  !![
  <stateStoreTask>
   <unitName>Node_Component_Hot_Halo_Standard_State_Store</unitName>
  </stateStoreTask>
  !!]
  subroutine Node_Component_Hot_Halo_Standard_State_Store(stateFile,gslStateFile,stateOperationID)
    !!{
    Store object state,
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Storing state for: componentHotHalo -> standard',verbosity=verbosityLevelInfo)
    !![
    <stateStore variables="cosmologyFunctions_ darkMatterHaloScale_ hotHaloMassDistribution_ hotHaloTemperatureProfile_ accretionHalo_ chemicalState_ hotHaloRamPressureStripping_ hotHaloRamPressureTimescale_ cosmologyParameters_ hotHaloOutflowReincorporation_ hotHaloOutflowStripping_ coolingInfallTorque_"/>
    !!]
    return
  end subroutine Node_Component_Hot_Halo_Standard_State_Store

  !![
  <stateRetrieveTask>
   <unitName>Node_Component_Hot_Halo_Standard_State_Restore</unitName>
  </stateRetrieveTask>
  !!]
  subroutine Node_Component_Hot_Halo_Standard_State_Restore(stateFile,gslStateFile,stateOperationID)
    !!{
    Retrieve object state.
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Retrieving state for: componentHotHalo -> standard',verbosity=verbosityLevelInfo)
    !![
    <stateRestore variables="cosmologyFunctions_ darkMatterHaloScale_ hotHaloMassDistribution_ hotHaloTemperatureProfile_ accretionHalo_ chemicalState_ hotHaloRamPressureStripping_ hotHaloRamPressureTimescale_ cosmologyParameters_ hotHaloOutflowReincorporation_ hotHaloOutflowStripping_ coolingInfallTorque_"/>
    !!]
    return
  end subroutine Node_Component_Hot_Halo_Standard_State_Restore

  function Node_Component_Hot_Halo_Standard_Mass_Distribution(self,componentType,massType,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the mass distribution associated with the hot halo.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentHotHaloStandard
    use :: Galactic_Structure_Options, only : enumerationWeightByType     , enumerationComponentTypeType, enumerationMassTypeType , massTypeGaseous, &
         &                                    componentTypeHotHalo
    use :: Mass_Distributions        , only : massDistributionClass       , kinematicsDistributionClass , massDistributionMatches_
    implicit none
    class  (massDistributionClass       ), pointer                 :: massDistribution_
    class  (kinematicsDistributionClass ), pointer                 :: kinematicsDistribution_
    class  (nodeComponentHotHaloStandard), intent(inout)           :: self
    type   (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type   (enumerationMassTypeType     ), intent(in   ), optional :: massType
    type   (enumerationWeightByType     ), intent(in   ), optional :: weightBy
    integer                              , intent(in   ), optional :: weightIndex
    !$GLC attributes unused :: weightIndex

    if (massDistributionMatches_(componentTypeHotHalo,massTypeGaseous,componentType,massType)) then
       massDistribution_ => hotHaloMassDistribution_%get(self%hostNode,weightBy,weightIndex)
       if (associated(massDistribution_)) then
          kinematicsDistribution_ => hotHaloTemperatureProfile_%get(self%hostNode)
          call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
          !![
          <objectDestructor name="kinematicsDistribution_"/>
          !!]
        end if
    else
       massDistribution_ => null()
    end if
    return
  end function Node_Component_Hot_Halo_Standard_Mass_Distribution

end module Node_Component_Hot_Halo_Standard
