!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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
Contains a module which implements the standard \glc{nsc} node component.
!!}

module Node_Component_NSC_Standard
  !!{
  Implements the standard \glc{nsc} node component.
  !!}
  use :: Dark_Matter_Halo_Scales         , only : darkMatterHaloScaleClass
  use :: Histories                       , only : history
  use :: Satellite_Merging_Mass_Movements, only : mergerMassMovementsClass
  use :: Satellite_Merging_Remnant_Sizes , only : mergerRemnantSizeClass
  use :: Star_Formation_Histories        , only : starFormationHistory            , starFormationHistoryClass
  use :: Stellar_Population_Properties   , only : stellarPopulationPropertiesClass
  implicit none
  private
  public :: Node_Component_NSC_Standard_Scale_Set        , Node_Component_NSC_Standard_Pre_Evolve         , &
       &    Node_Component_NSC_Standard_Post_Step        , Node_Component_NSC_Standard_Thread_Uninitialize, &
       &    Node_Component_NSC_Standard_Initialize       , Node_Component_NSC_Standard_Inactive           , &
       &    Node_Component_NSC_Standard_State_Store      , Node_Component_NSC_Standard_State_Retrieve     , &
       &    Node_Component_NSC_Standard_Thread_Initialize                     

  !![
  <component>
   <class>NSC</class>
   <name>standard</name>
   <isDefault>false</isDefault>
   <createFunction isDeferred="true" />
   <properties>
    <property>
      <name>isInitialized</name>
      <type>logical</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
    </property>
    <property>
      <name>isCollapsed</name>
      <type>logical</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false"/>
    </property>
    <property>
      <name>massStellar</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output unitsInSI="massSolar" comment="Mass of stars in the standard nuclear star cluster."/>
    </property>
    <property>
      <name>massStellarFormed</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
    </property>
    <property>
      <name>abundancesStellar</name>
      <type>abundances</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output unitsInSI="massSolar" comment="Mass of metals in the stellar phase of the standard nuclear star cluster."/>
    </property>
    <property>
      <name>massGas</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" makeGeneric="false" />
      <output unitsInSI="massSolar" comment="Mass of gas in the standard nuclear star cluster."/>
    </property>
    <property>
      <name>abundancesGas</name>
      <type>abundances</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" makeGeneric="false" />
      <output unitsInSI="massSolar" comment="Mass of metals in the gas phase of the standard nuclear star cluster."/>
    </property>
    <property>
      <name>angularMomentum</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" makeGeneric="false" />
      <output unitsInSI="massSolar*megaParsec*kilo" comment="Angular momentum of the standard nuclear star cluster."/>
      <getFunction>Node_Component_NSC_Standard_Angular_Momentum</getFunction>
    </property>
    <property>
      <name>radius</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
      <output unitsInSI="megaParsec" comment="Radial scale length in the standard nuclear star cluster."/>
      <getFunction>Node_Component_NSC_Standard_Radius</getFunction>
    </property>
    <property>
      <name>halfMassRadius</name>
      <attributes isSettable="false" isGettable="true" isEvolvable="false" isVirtual="true" />
      <type>double</type>
      <rank>0</rank>
      <output unitsInSI="megaParsec" comment="Radial scale length in the standard nuclear star cluster."/>
      <getFunction>Node_Component_NSC_Standard_Half_Mass_Radius</getFunction>
    </property>
    <property>
      <name>velocity</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
      <getFunction>Node_Component_NSC_Standard_Velocity</getFunction>
      <output unitsInSI="kilo" comment="Circular velocity of the standard nuclear star cluster at scale length."/>
    </property>
    <property>
      <name>luminositiesStellar</name>
      <type>stellarLuminosities</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output unitsInSI="luminosityZeroPointAB" comment="Luminosity of nuclear star cluster stars."/>
    </property>
    <property>
      <name>stellarPropertiesHistory</name>
      <type>history</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" isDeferred="rate" createIfNeeded="true" />
    </property>
    <property>
      <name>starFormationHistory</name>
      <type>history</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" isDeferred="rate" createIfNeeded="true" />
    </property>
    <property>
      <name>massGasSink</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="false" isGettable="false" isEvolvable="true" isDeferred="rate" isVirtual="true" />
    </property>
   </properties>
   <bindings>
    <binding method="massDistribution" function="Node_Component_NSC_Standard_Mass_Distribution" bindsTo="component"/>
    <binding method="massBaryonic"     function="Node_Component_NSC_Standard_Mass_Baryonic"     bindsTo="component"/>
   </bindings>
   <functions>objects.nodes.components.nuclear_star_cluster.standard.bound_functions.inc</functions>
  </component>
  !!]

  ! Objects used by this component.
  class(stellarPopulationPropertiesClass), pointer :: stellarPopulationProperties_
  class(darkMatterHaloScaleClass        ), pointer :: darkMatterHaloScale_
  class(starFormationHistoryClass       ), pointer :: starFormationHistory_
  class(mergerMassMovementsClass        ), pointer :: mergerMassMovements_
  class(mergerRemnantSizeClass          ), pointer :: mergerRemnantSize_
  !$omp threadprivate(stellarPopulationProperties_,darkMatterHaloScale_,starFormationHistory_,mergerMassMovements_,mergerRemnantSize_)

  ! Internal count of abundances.
  integer                                     :: abundancesCount
  
  ! Storage for the star formation and stellar properties histories time range, used when extending this range.
  double precision, allocatable, dimension(:) :: starFormationHistoryTemplate       , stellarPropertiesHistoryTemplate
  !$omp threadprivate(starFormationHistoryTemplate,stellarPropertiesHistoryTemplate)

  ! Parameters controlling the physical implementation.
  double precision                            :: toleranceAbsoluteMass              , toleranceRelativeMetallicity          
  logical                                     :: inactiveLuminositiesStellar

  ! The largest and smallest angular momentum, in units of that of a circular orbit at the virial radius, considered to be
  ! physically plausible for a nuclear star cluster.
  double precision, parameter                 :: angularMomentumMaximum      =1.0d+1
  double precision, parameter                 :: angularMomentumMinimum      =1.0d-6

  ! A threadprivate object used to track to which thread events are attached.
  integer                                     :: thread
  !$omp threadprivate(thread)

contains

  !![
  <nodeComponentInitializationTask>
   <unitName>Node_Component_NSC_Standard_Initialize</unitName>
  </nodeComponentInitializationTask>
  !!]
  subroutine Node_Component_NSC_Standard_Initialize(parameters)
    !!{
    Initializes the tree node standard nuclear star cluster methods module.
    !!}
    use :: Abundances_Structure            , only : Abundances_Property_Count
    use :: Error                           , only : Error_Report
    use :: Galacticus_Nodes                , only : defaultNSCComponent      , nodeComponentNSCStandard
    use :: Input_Parameters                , only : inputParameter           , inputParameters
    use :: Node_Component_NSC_Standard_Data, only : radiusNormalization
    implicit none
    type(inputParameters         ), intent(inout) :: parameters
    type(nodeComponentNSCStandard)                :: NSCStandardComponent
    type(inputParameters         )                :: subParameters

    if (defaultNSCComponent%standardIsActive()) then
       ! Get number of abundance properties.
       abundancesCount  =Abundances_Property_Count            ()
       ! Bind deferred functions.
       call NSCStandardComponent%             massGasSinkRateFunction(Node_Component_NSC_Standard_Mass_Gas_Sink_Rate         )
       call NSCStandardComponent%    starFormationHistoryRateFunction(Node_Component_NSC_Standard_Star_Formation_History_Rate)
       call NSCStandardComponent%stellarPropertiesHistoryRateFunction(Node_Component_NSC_Standard_Stellar_Prprts_History_Rate)
       ! Find our parameters.
       subParameters=parameters%subParameters('componentNSC')
       ! Read parameters controlling the physical implementation.
       !![
       <inputParameter>
          <name>radiusNormalization</name>
          <defaultValue>3.3d-6</defaultValue>
          <description>The initial value appearing in the radius-mass relation</description>
          <source>subParameters</source>
       </inputParameter>
       <inputParameter>
         <name>toleranceAbsoluteMass</name>
         <defaultValue>1.0d-6</defaultValue>
         <description>The mass tolerance used to judge whether the nuclear star cluster is physically plausible.</description>
         <source>subParameters</source>
       </inputParameter>
       <inputParameter>
         <name>toleranceRelativeMetallicity</name>
         <defaultValue>1.0d-4</defaultValue>
         <description>The metallicity tolerance for ODE solution.</description>
         <source>subParameters</source>
       </inputParameter>
       <inputParameter>
         <name>inactiveLuminositiesStellar</name>
         <defaultValue>.false.</defaultValue>
         <description>Specifies whether or not nuclear star cluster stellar luminosities are inactive properties (i.e. do not appear in any ODE being solved).</description>
         <source>subParameters</source>
       </inputParameter>

       !!]
    end if
    return
  end subroutine Node_Component_NSC_Standard_Initialize

  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Node_Component_NSC_Standard_Thread_Initialize</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Node_Component_NSC_Standard_Thread_Initialize(parameters)
    !!{
    Initializes the standard nuclear star cluster component module for each thread.
    !!}
    use :: Events_Hooks                    , only : dependencyDirectionAfter       , dependencyRegEx            , openMPThreadBindingAtLevel, &
          &                                         postEvolveEvent                , satelliteMergerEvent       , mergerTreeExtraOutputEvent
    use :: Error                           , only : Error_Report
    use :: Galacticus_Nodes                , only : defaultNSCComponent
    use :: Input_Parameters                , only : inputParameter                 , inputParameters
    use :: Mass_Distributions              , only : massDistributionSpherical      , kinematicsDistributionLocal
    use :: Node_Component_NSC_Standard_Data, only : massDistributionStellar_       , massDistributionGas_       , kinematicDistribution_                
    use :: Galactic_Structure_Options      , only : componentTypeNuclearStarCluster, massTypeStellar            , massTypeGaseous
    implicit none
    type(inputParameters), intent(inout) :: parameters
    type(dependencyRegEx), dimension(1)  :: dependencies
    type(inputParameters)                :: subParameters

    ! Check if this implementation is selected. If so, initialize the mass distribution.
    if (defaultNSCComponent%standardIsActive()) then
       dependencies(1)=dependencyRegEx(dependencyDirectionAfter,'^remnantStructure:')
       call satelliteMergerEvent      %attach(thread,satelliteMerger      ,openMPThreadBindingAtLevel,label='nodeComponentNSCStandard',dependencies=dependencies)
       call postEvolveEvent           %attach(thread,postEvolve           ,openMPThreadBindingAtLevel,label='nodeComponentNSCStandard'                          )
       call mergerTreeExtraOutputEvent%attach(thread,mergerTreeExtraOutput,openMPThreadBindingAtLevel,label='nodeComponentNSCStandard'                          )
       ! Find our parameters.
       subParameters=parameters%subParameters('componentNSC')
       !![
       <objectBuilder class="stellarPopulationProperties"                                     name="stellarPopulationProperties_" source="subParameters"                    />
       <objectBuilder class="darkMatterHaloScale"                                             name="darkMatterHaloScale_"         source="subParameters"                    />
       <objectBuilder class="starFormationHistory"                                            name="starFormationHistory_"        source="subParameters"                    />
       <objectBuilder class="mergerMassMovements"                                             name="mergerMassMovements_"         source="subParameters"                    />
       <objectBuilder class="mergerRemnantSize"                                               name="mergerRemnantSize_"           source="subParameters"                    />
       <objectBuilder class="massDistribution"            parameterName="massDistributionNSC" name="massDistributionStellar_"     source="subParameters" threadPrivate="yes" >
        <default>
         <massDistributionNSC value="hernquist">
          <dimensionless value="true"/>
         </massDistributionNSC>
        </default>
       </objectBuilder>
       !!]
       ! Validate the NSC mass distribution
       select type (massDistributionStellar_)
       class is (massDistributionSpherical)
          ! The NSC mass distribution must have spherical symmetry. So, this is acceptable.
       class default 
          call Error_Report('only spehrically symmetric mass distributions are allowed'//{introspection:location})
       end select
       if (.not.massDistributionStellar_%isDimensionless()) call Error_Report('Nuclear star cluster mass distribution must be dimensionless'//{introspection:location})
       ! Duplicate the dimensionless mass distribution to use for the gas component, and set component and mass type in both.
       !$omp critical(NSCStandardDeepCopy)
       allocate(massDistributionGas_,mold=massDistributionStellar_)
       !![
       <deepCopyReset variables="massDistributionStellar_"/>
       <deepCopy source="massDistributionStellar_" destination="massDistributionGas_"/>
       <deepCopyFinalize variables="massDistributionGas_"/>
       !!]
       !$omp end critical(NSCStandardDeepCopy)
       call massDistributionStellar_%setTypes(componentTypeNuclearStarCluster,massTypeStellar)
       call massDistributionGas_    %setTypes(componentTypeNuclearStarCluster,massTypeGaseous)
       ! Construct the kinematic distribution
       allocate(kinematicDistribution_)
       !![
       <referenceConstruct object="kinematicDistribution_" constructor="kinematicsDistributionLocal(alpha=1.0d0/sqrt(2.0d0))"/>
       !!]    
    end if
    return
  end subroutine Node_Component_NSC_Standard_Thread_Initialize

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Node_Component_NSC_Standard_Thread_Uninitialize</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Node_Component_NSC_Standard_Thread_Uninitialize()
    !!{
    Uninitializes the standard nuclear star cluster component module for each thread.
    !!}
    use :: Events_Hooks                    , only : postEvolveEvent         , satelliteMergerEvent, mergerTreeExtraOutputEvent
    use :: Galacticus_Nodes                , only : defaultNSCComponent
    use :: Node_Component_NSC_Standard_Data, only : massDistributionStellar_, massDistributionGas_, kinematicDistribution_
    implicit none

    if (defaultNSCComponent%standardIsActive()) then
       if (satelliteMergerEvent      %isAttached(thread,satelliteMerger      )) call satelliteMergerEvent      %detach(thread,satelliteMerger      )
       if (postEvolveEvent           %isAttached(thread,postEvolve           )) call postEvolveEvent           %detach(thread,postEvolve           )
       if (mergerTreeExtraOutputEvent%isAttached(thread,mergerTreeExtraOutput)) call mergerTreeExtraOutputEvent%detach(thread,mergerTreeExtraOutput)
       !![
       <objectDestructor name="stellarPopulationProperties_"/>
       <objectDestructor name="darkMatterHaloScale_"        />
       <objectDestructor name="starFormationHistory_"       />
       <objectDestructor name="mergerMassMovements_"        />
       <objectDestructor name="mergerRemnantSize_"          />
       <objectDestructor name="massDistributionStellar_"    />
       <objectDestructor name="massDistributionGas_"        />
       <objectDestructor name="kinematicDistribution_"      />
       !!]
    end if
    return
  end subroutine Node_Component_NSC_Standard_Thread_Uninitialize

  !![
  <preEvolveTask>
  <unitName>Node_Component_NSC_Standard_Pre_Evolve</unitName>
  </preEvolveTask>
  !!]
  subroutine Node_Component_NSC_Standard_Pre_Evolve(node)
    !!{
    Ensure the nuclear star cluster has been initialized.
    !!}
    use :: Galacticus_Nodes, only : defaultNSCComponent, nodeComponentNSC, nodeComponentNSCStandard, treeNode
    implicit none
    type (treeNode        ), intent(inout), pointer :: node
    class(nodeComponentNSC)               , pointer :: nuclearStarCluster

    ! Check if we are the default method.
    if (.not.defaultNSCComponent%standardIsActive()) return
    ! Get the nuclear star cluster component.
    nuclearStarCluster => node%NSC()
    ! Check if an standard nuclear star cluster component exists.
    select type (nuclearStarCluster)
    class is (nodeComponentNSCStandard)
       ! Initialize the NSC.
       call Node_Component_NSC_Standard_Create(node)
    end select
    return
  end subroutine Node_Component_NSC_Standard_Pre_Evolve

  subroutine postEvolve(self,node)
    !!{
    Trim histories attached to the nuclear star cluster.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentNSC, nodeComponentNSCStandard, treeNode
    use :: Histories       , only : history
    implicit none
    class(*                 ), intent(inout) :: self
    type (treeNode          ), intent(inout) :: node
    class(nodeComponentNSC  ), pointer       :: nuclearStarCluster
    class(nodeComponentBasic), pointer       :: basic
    type (history           )                :: stellarPropertiesHistory
    !$GLC attributes unused :: self

    ! Get the nuclear star cluster component.
    nuclearStarCluster => node%NSC()
    ! Check if an standard nuclear star cluster component exists.
    select type (nuclearStarCluster)
    class is (nodeComponentNSCStandard)
       ! Trim the stellar populations properties future history.
       basic                    => node              %basic                   ()
       stellarPropertiesHistory =  nuclearStarCluster%stellarPropertiesHistory()
       call stellarPropertiesHistory%trim(basic%time())
       call nuclearStarCluster%stellarPropertiesHistorySet(stellarPropertiesHistory)
    end select
    return
  end subroutine postEvolve

  !![
  <postStepTask>
    <unitName>Node_Component_NSC_Standard_Post_Step</unitName>
  </postStepTask>
  !!]
  subroutine Node_Component_NSC_Standard_Post_Step(node,status)
    !!{
    Trim histories attached to the nuclear star cluster.
    !!}
    use :: Abundances_Structure          , only : abs                , zeroAbundances
    use :: Display                       , only : displayMessage     , verbosityLevelWarn
    use :: Galacticus_Nodes              , only : defaultNSCComponent, nodeComponentNSC   , nodeComponentNSCStandard, treeNode
    use :: Histories                     , only : history
    use :: Interface_GSL                 , only : GSL_Success        , GSL_Continue
    use :: ISO_Varying_String            , only : assignment(=)      , operator(//)       , varying_string
    use :: Stellar_Luminosities_Structure, only : abs                , stellarLuminosities
    use :: String_Handling               , only : operator(//)
    implicit none
    type            (treeNode           ), intent(inout), pointer :: node
    integer                              , intent(inout)          :: status
    class           (nodeComponentNSC   )               , pointer :: nuclearStarCluster
    double precision                     , parameter              :: angularMomentumTolerance=1.0d-2
    double precision                     , save                   :: fractionalErrorMaximum  =0.0d+0
    double precision                                              :: massNuclearStarCluster         , fractionalError, &
         &                                                           specificAngularMomentum
    character       (len=20             )                         :: valueString
    type            (varying_string     ), save                   :: message
    !$omp threadprivate(message)
    type            (stellarLuminosities), save                   :: luminositiesStellar
    !$omp threadprivate(luminositiesStellar)
    type            (history            ), save                   :: historyStellar
    !$omp threadprivate(historyStellar)

    ! Return immediately if this class is not in use.
    if (.not.defaultNSCComponent%standardIsActive()) return
    ! Get the nuclear star cluster component.
    nuclearStarCluster => node%NSC()
    ! Check if an standard nuclear star cluster component exists.
    select type (nuclearStarCluster)
    class is (nodeComponentNSCStandard)
       ! Note that "status" is not set to failure as these changes in state of the nuclear star cluster should not change any calculation of
       ! differential evolution rates as a negative gas/stellar mass was unphysical anyway.
       !
       ! Trap negative gas masses.
       if (nuclearStarCluster%massGas() < 0.0d0) then
          ! Check if this exceeds the maximum previously recorded error.
          fractionalError=+  abs(nuclearStarCluster%massGas    ()) &
               &          /(                                       &
               &            +abs(nuclearStarCluster%massGas    ()) &
               &            +abs(nuclearStarCluster%massStellar()) &
               &           )
          !$omp critical (Standard_NSC_Post_Evolve_Check)
          if (fractionalError > fractionalErrorMaximum) then
             ! Report a warning.
             message='Warning: Nuclear star cluster has negative gas mass (fractional error exceeds any previously reported):'//char(10)
             message=message//'  Node index                        = '//node%index()     //char(10)
             write (valueString,'(e12.6)') nuclearStarCluster%massGas    ()
             message=message//'  Nuclear star cluster gas mass     = '//trim(valueString)//char(10)
             write (valueString,'(e12.6)') nuclearStarCluster%massStellar()
             message=message//'  Nuclear star Cluster stellar mass = '//trim(valueString)//char(10)
             write (valueString,'(e12.6)') fractionalError
             message=message//'  Error measure                     = '//trim(valueString)//char(10)
             if (fractionalErrorMaximum == 0.0d0) then
                ! This is the first time this warning has been issued, so give some extra information.
                message=message//'  Gas mass will be reset to zero (in future cases also).'                                 //char(10)
                message=message//'  Future cases will be reported only when they exceed the previous maximum error measure.'//char(10)
                message=message//'  Negative masses are due to numerically inaccuracy in the ODE solutions.'                //char(10)
                message=message//'  If significant, consider using a higher tolerance in the ODE solver.'
             end if
             call displayMessage(message,verbosityLevelWarn)
             ! Store the new maximum fractional error.
             fractionalErrorMaximum=fractionalError
          end if
          !$omp end critical (Standard_NSC_Post_Evolve_Check)
          ! Get the specific angular momentum of the nuclear star cluster material
          massNuclearStarCluster=+nuclearStarCluster%massGas    () &
               &                 +nuclearStarCluster%massStellar()
          if (massNuclearStarCluster == 0.0d0) then
             specificAngularMomentum=0.0d0
             call nuclearStarCluster%         massStellarSet    (                  0.0d0)
             call nuclearStarCluster%   abundancesStellarSet    (         zeroAbundances)
             historyStellar         =nuclearStarCluster%starFormationHistory    ()
             call historyStellar    %reset()
             call nuclearStarCluster%starFormationHistorySet    (         historyStellar)
             historyStellar         =nuclearStarCluster%stellarPropertiesHistory()
             call historyStellar    %reset()
             call nuclearStarCluster%stellarPropertiesHistorySet(         historyStellar)
             ! We need to reset the stellar luminosities to zero. We can't simply use the "zeroStellarLuminosities" instance since
             ! our luminosities may have been truncated. If we were to use "zeroStellarLuminosities" then the number of stellar
             ! luminosities associated with the nuclear star cluster would change - but we are in the middle of differential evolution here and we
             ! cannot change the number of evolvable properties as doing so will lead to invalid memory accesses during
             ! deserialization of properties from the ODE solver.
             call luminositiesStellar%destroy               (                   )
             luminositiesStellar=nuclearStarCluster%luminositiesStellar()
             call luminositiesStellar%reset                 (                   )
             call nuclearStarCluster %luminositiesStellarSet(luminositiesStellar)
          else
             specificAngularMomentum=nuclearStarCluster%angularMomentum()/massNuclearStarCluster
             if (specificAngularMomentum < 0.0d0)                          &
                  & specificAngularMomentum=+nuclearStarCluster%radius  () &
                  &                         *nuclearStarCluster%velocity()
          end if
          ! Reset the gas, abundances and angular momentum of the nuclear star cluster.
          call nuclearStarCluster%        massGasSet(                                                   0.0d0)
          call nuclearStarCluster%  abundancesGasSet(                                          zeroAbundances)
          call nuclearStarCluster%angularMomentumSet(specificAngularMomentum*nuclearStarCluster%massStellar())
          ! Indicate that ODE evolution should continue after this state change.
          if (status == GSL_Success) status=GSL_Continue
       end if
       ! Trap negative stellar masses.
       if (nuclearStarCluster%massStellar() < 0.0d0) then
          ! Check if this exceeds the maximum previously recorded error.
          fractionalError=+  abs(nuclearStarCluster%massStellar()) &
               &          /(                                       &
               &            +abs(nuclearStarCluster%massGas    ()) &
               &            +abs(nuclearStarCluster%massStellar()) &
               &           )
          !$omp critical (Standard_NSC_Post_Evolve_Check)
          if (fractionalError > fractionalErrorMaximum) then
             ! Report a warning.
             message='Warning: Nuclear star cluster has negative stellar mass (fractional error exceeds any previously reported):'//char(10)
             message=message//'  Node index                        = '//node%index()     //char(10)
             write (valueString,'(e12.6)') nuclearStarCluster%massGas    ()
             message=message//'  Nuclear star cluster gas mass     = '//trim(valueString)//char(10)
             write (valueString,'(e12.6)') nuclearStarCluster%massStellar()
             message=message//'  Nuclear star cluster stellar mass = '//trim(valueString)//char(10)
             write (valueString,'(e12.6)') fractionalError
             message=message//'  Error measure                     = '//trim(valueString)//char(10)
             if (fractionalErrorMaximum == 0.0d0) then
                ! This is the first time this warning has been issued, so give some extra information.
                message=message//'  Stellar mass will be reset to zero (in future cases also).'                             //char(10)
                message=message//'  Future cases will be reported only when they exceed the previous maximum error measure.'//char(10)
                message=message//'  Negative masses are due to numerically inaccuracy in the ODE solutions.'                //char(10)
                message=message//'  If significant, consider using a higher tolerance in the ODE solver.'
             end if
             call displayMessage(message,verbosityLevelWarn)
             ! Store the new maximum fractional error.
             fractionalErrorMaximum=fractionalError
          end if
          !$omp end critical (Standard_NSC_Post_Evolve_Check)
          ! Get the specific angular momentum of the nuclear star cluster material
          massNuclearStarCluster=+nuclearStarCluster%massGas    () &
               &                 +nuclearStarCluster%massStellar()
          if (massNuclearStarCluster == 0.0d0) then
             specificAngularMomentum=0.0d0
             call nuclearStarCluster%      massGasSet(         0.0d0)
             call nuclearStarCluster%abundancesGasSet(zeroAbundances)
          else
             specificAngularMomentum       =+nuclearStarCluster%angularMomentum       () &
                  &                         /                   massNuclearStarCluster
             if (specificAngularMomentum < 0.0d0) &
                  & specificAngularMomentum=+nuclearStarCluster%radius  () &
                  &                         *nuclearStarCluster%velocity()
          end if
          ! Reset the stellar, abundances and angular momentum of the nuclear star cluster.
          call nuclearStarCluster%      massStellarSet(                                               0.0d0)
          call nuclearStarCluster%abundancesStellarSet(                                      zeroAbundances)
          call nuclearStarCluster%  angularMomentumSet(specificAngularMomentum*nuclearStarCluster%massGas())
          historyStellar=nuclearStarCluster%starFormationHistory ()
          call historyStellar    %reset                      (              )
          call nuclearStarCluster%starFormationHistorySet    (historyStellar)
          historyStellar=nuclearStarCluster%stellarPropertiesHistory()
          call historyStellar    %reset                      (              )
          call nuclearStarCluster%stellarPropertiesHistorySet(historyStellar)
          ! Indicate that ODE evolution should continue after this state change.
          if (status == GSL_Success) status=GSL_Continue
       end if
       ! Trap negative angular momentum.
       if (nuclearStarCluster%angularMomentum() < 0.0d0) then
          ! Estimate a reasonable specific angular momentum.
          specificAngularMomentum=+nuclearStarCluster%radius     () &
               &                  *nuclearStarCluster%velocity   ()
          ! Get the mass of the NSC.
          massNuclearStarCluster =+nuclearStarCluster%massGas    () &
              &                   +nuclearStarCluster%massStellar()
          call nuclearStarCluster%angularMomentumSet(specificAngularMomentum*massNuclearStarCluster)
          ! Indicate that ODE evolution should continue after this state change.
          if (status == GSL_Success) status=GSL_Continue
       end if
    end select
    return
  end subroutine Node_Component_NSC_Standard_Post_Step

  subroutine Node_Component_NSC_Standard_Mass_Gas_Sink_Rate(self,rate,interrupt,interruptProcedure)
    !!{
    Account for a sink of gaseous material in the standard nuclear star cluster.
    !!}
    use :: Abundances_Structure, only : operator(*)
    use :: Error               , only : Error_Report
    use :: Galacticus_Nodes    , only : interruptTask, nodeComponentNSC
    implicit none
    class           (nodeComponentNSC), intent(inout)                    :: self
    logical                           , intent(inout), optional          :: interrupt
    procedure       (interruptTask   ), intent(inout), optional, pointer :: interruptProcedure
    double precision                  , intent(in   )                    :: rate
    double precision                                                     :: gasMass           , stellarMass
    !$GLC attributes unused :: interrupt, interruptProcedure

    ! Trap cases where an attempt is made to add gas via this sink function.
    if (rate > 0.0d0) call Error_Report(                                                                  &
         &                              'attempt to add mass via sink in standard nuclear star cluster'// &
         &                              {introspection:location}                                          &
         &                             )
    ! Return if no adjustment is being made.
    if (rate == 0.0d0) return
    ! Get the gas mass present.
    gasMass    =self%massGas    ()
    ! Get the stellar mass present.
    stellarMass=self%massStellar()
    ! If gas is present, adjust the rates.
    if (gasMass > 0.0d0 .and. gasMass+stellarMass > 0.0d0) then
       call self%        massGasRate( rate                                              )
       call self%angularMomentumRate((rate/(gasMass+stellarMass))*self%angularMomentum())
       call self%  abundancesGasRate((rate/ gasMass             )*self%abundancesGas  ())
    end if
    return
  end subroutine Node_Component_NSC_Standard_Mass_Gas_Sink_Rate

  subroutine Node_Component_NSC_Standard_Star_Formation_History_Rate(self,rate,interrupt,interruptProcedure)
    !!{
    Adjust the rates for the star formation history.
    !!}
    use :: Galacticus_Nodes , only : interruptTask, nodeComponentNSC, nodeComponentNSCStandard
    implicit none
    class    (nodeComponentNSC), intent(inout)                    :: self
    type     (history         ), intent(in   )                    :: rate
    logical                    , intent(inout), optional          :: interrupt
    procedure(interruptTask   ), intent(inout), optional, pointer :: interruptProcedure
    type     (history         )                                   :: starFormationHistory

    ! Get the star formation history in the spheroid.
    starFormationHistory=self%starFormationHistory()
    ! Check if the range of this history is sufficient.
    if (starFormationHistory_%rangeIsSufficient(starFormationHistory,rate)) then
       ! Range is sufficient, call the intrinsic rate function.
       select type (self)
       class is (nodeComponentNSCStandard)
          call self%starFormationHistoryRateIntrinsic(rate)
       end select
    else
       ! Range is insufficient.
       if (allocated(starFormationHistoryTemplate)) deallocate(starFormationHistoryTemplate)
       allocate(starFormationHistoryTemplate(size(rate%time)))
       starFormationHistoryTemplate =  rate%time
       interrupt                    =  .true.
       interruptProcedure           => Node_Component_NSC_Standard_Star_Formation_History_Extend
    end if
    return
  end subroutine Node_Component_NSC_Standard_Star_Formation_History_Rate

  subroutine Node_Component_NSC_Standard_Stellar_Prprts_History_Rate(self,rate,interrupt,interruptProcedure)
    !!{
    Adjust the rates for the stellar properties history.
    !!}
    use :: Error            , only : Error_Report
    use :: Galacticus_Nodes , only : interruptTask, nodeComponentNSC, nodeComponentNSCStandard
    implicit none
    class    (nodeComponentNSC), intent(inout)                    :: self
    type     (history         ), intent(in   )                    :: rate
    logical                    , intent(inout), optional          :: interrupt
    procedure(interruptTask   ), intent(inout), optional, pointer :: interruptProcedure
    type     (history         )                                   :: stellarPropertiesHistory

    ! Get the star formation history in the spheroid.
    stellarPropertiesHistory=self%stellarPropertiesHistory()
    ! Ensure that the history already exists.
    if (.not.stellarPropertiesHistory%exists())                                         &
         & call Error_Report(                                                           &
         &                   'no star formation history has been created in spheroid'// &
         &                   {introspection:location}                                   &
         &                  )
    ! Check if the star formation history in the spheroid spans a sufficient range to accept the input rates.
    if     (                                                                                                     &
         &       rate%time(              1) < stellarPropertiesHistory%time(                                  1) &
         &  .or. rate%time(size(rate%time)) > stellarPropertiesHistory%time(size(stellarPropertiesHistory%time)) &
         & ) then
       ! It does not, so interrupt evolution and extend the history.
       if (allocated(stellarPropertiesHistoryTemplate)) deallocate(stellarPropertiesHistoryTemplate)
       allocate(stellarPropertiesHistoryTemplate(size(rate%time)))
       stellarPropertiesHistoryTemplate=rate%time
       interrupt=.true.
       interruptProcedure => Node_Component_NSC_Standard_Stellar_Prprts_History_Extend
       return
    end if
    ! Adjust the rate.
    select type (self)
    class is (nodeComponentNSCStandard)
       call self%stellarPropertiesHistoryRateIntrinsic(rate)
    end select
    return
  end subroutine Node_Component_NSC_Standard_Stellar_Prprts_History_Rate

  !![
  <scaleSetTask>
   <unitName>Node_Component_NSC_Standard_Scale_Set</unitName>
  </scaleSetTask>
  !!]
  subroutine Node_Component_NSC_Standard_Scale_Set(node)
    !!{
    Set scales for properties of {\normalfont \ttfamily node}.
    !!}
    use :: Abundances_Structure          , only : abs                   , abundances      , max                     , operator(*)            , &
          &                                       unitAbundances
    use :: Galacticus_Nodes              , only : defaultNSCComponent   , nodeComponentNSC, nodeComponentNSCStandard, nodeComponentSpheroid  , &
          &                                       treeNode
    use :: Histories                     , only : history
    use :: Stellar_Luminosities_Structure, only : abs                   , max             , stellarLuminosities     , unitStellarLuminosities
    implicit none
    type            (treeNode             ), intent(inout), pointer :: node
    class           (nodeComponentNSC     )               , pointer :: nuclearStarCluster
    class           (nodeComponentSpheroid)               , pointer :: spheroid
    double precision                       , parameter              :: massMinimum                   =1.0d0
    double precision                       , parameter              :: scaleMassRelative             =1.0d-3
    double precision                       , parameter              :: angularMomentumMinimum        =1.0d-1
    double precision                       , parameter              :: fractionTolerance             =1.0d-4
    double precision                       , parameter              :: luminosityMinimum             =1.0d+0
    double precision                                                :: angularMomentum                       , mass
    type            (history              )                         :: stellarPopulationHistoryScales
    type            (stellarLuminosities  )                         :: stellarLuminositiesScale
    type            (abundances           )                         :: abundancesScale

    ! Check if we are the default method.
    if (.not.defaultNSCComponent%standardIsActive()) return
    ! Get the nuclear star cluster component.
    nuclearStarCluster => node%NSC()
    ! Check if an standard nuclear star cluster component exists.
    select type (nuclearStarCluster)
    class is (nodeComponentNSCStandard)
       ! Get spheroid component.
       spheroid  => node%spheroid()
       ! Set scale for angular momentum.
       angularMomentum=+abs(nuclearStarCluster    %angularMomentum()) &
            &          +abs(spheroid%angularMomentum())
       call nuclearStarCluster%angularMomentumScale(max(angularMomentum,angularMomentumMinimum))
       ! Set scale for masses.
       ! The scale here (and for other quantities below) combines the mass of nuclear star cluster and spheroid.
        mass     =    max(     scaleMassRelative*spheroid          %massStellar(), &
            &              max(                  nuclearStarCluster%massStellar(), &
            &                  massMinimum                                         &
            &                 )                                                    &
            &             ) 
       call nuclearStarCluster%massGasScale          (mass)
       call nuclearStarCluster%massStellarScale      (mass)
       call nuclearStarCluster%massStellarFormedScale(mass)
       ! Set scales for abundances if necessary.
       if (abundancesCount > 0) then
          ! Set scale for abundances.
          abundancesScale=+max(                                &
               &               +mass                           &
               &               *toleranceRelativeMetallicity , &
               &               +massMinimum                    &
               &              )                                &
               &          *unitAbundances
          ! Set scale for gas abundances.
          call nuclearStarCluster%abundancesGasScale    (abundancesScale)
          ! Set scale for stellar abundances.
          call nuclearStarCluster%abundancesStellarScale(abundancesScale)
       end if
       ! Set scale for stellar luminosities.
       stellarLuminositiesScale=max(                                                 &
            &                       +abs(nuclearStarCluster%luminositiesStellar())   &
            &                       +abs(spheroid          %luminositiesStellar()) , &
            &                           +                   unitStellarLuminosities  &
            &                           *                   luminosityMinimum        &
            &                      )
       call stellarLuminositiesScale%truncate                (nuclearStarCluster%luminositiesStellar     ())
       call nuclearStarCluster      %luminositiesStellarScale(                   stellarLuminositiesScale  )
       ! Set scales for stellar population properties and star formation histories.
       stellarPopulationHistoryScales=nuclearStarCluster%stellarPropertiesHistory()
       call stellarPopulationProperties_  %scales                       (nuclearStarCluster%massStellar(),nuclearStarCluster%abundancesStellar(),stellarPopulationHistoryScales)
       call nuclearStarCluster            %stellarPropertiesHistoryScale(                                                                        stellarPopulationHistoryScales)
       call stellarPopulationHistoryScales%destroy                      (                                                                                                      )
       stellarPopulationHistoryScales=nuclearStarCluster%starFormationHistory    ()
       call starFormationHistory_         %scales                       (stellarPopulationHistoryScales,node,nuclearStarCluster%massStellar(),nuclearStarCluster%massGas(),nuclearStarCluster%abundancesStellar())
       call nuclearStarCluster            %starFormationHistoryScale    (stellarPopulationHistoryScales                                                                                                          )
       call stellarPopulationHistoryScales%destroy                      (                                                                                                                                        )
    end select
    return
  end subroutine Node_Component_NSC_Standard_Scale_Set

  subroutine Node_Component_NSC_Standard_Create(node)
    !!{
    Create properties in an standard nuclear star cluster component.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentNSC, nodeComponentSpheroid ,treeNode
    use :: Histories       , only : history
    implicit none
    type   (treeNode             ), intent(inout) , target  :: node
    class  (nodeComponentNSC     )                , pointer :: nuclearStarCluster
    class  (nodeComponentSpheroid)                , pointer :: spheroid
    class  (nodeComponentBasic   )                , pointer :: basic
    type   (history              )                          :: historyStarFormation        , stellarPropertiesHistory      , &
         &                                                     spheroidStarFormationHistory   
    logical                                                 :: createStarFormationHistory  , createStellarPropertiesHistory
    double precision                                        :: timeBegin

    ! Get the nuclear star cluster component.
    nuclearStarCluster => node%NSC()
    ! Exit if already initialized.
    if (nuclearStarCluster%isInitialized()) return
    ! Determine which histories must be created.
    historyStarFormation          =nuclearStarCluster%starFormationHistory           ()
    createStarFormationHistory    =.not.             historyStarFormation    %exists ()
    call                                             historyStarFormation    %destroy()
    stellarPropertiesHistory      =nuclearStarCluster%stellarPropertiesHistory       ()
    createStellarPropertiesHistory=.not.             stellarPropertiesHistory%exists ()
    call                                             stellarPropertiesHistory%destroy()
    ! Create the stellar properties history.
    if (createStellarPropertiesHistory) then
       ! Create the stellar properties history.
       call stellarPopulationProperties_%historyCreate(node,stellarPropertiesHistory)
       call nuclearStarCluster%stellarPropertiesHistorySet(stellarPropertiesHistory)
    end if
    ! Create the star formation history.
    if (createStarFormationHistory) then
       spheroid => node%spheroid()
       spheroidStarFormationHistory=spheroid%starFormationHistory()
       if (spheroidStarFormationHistory%exists()) then
          timeBegin = spheroidStarFormationHistory%time(1)
       else
          basic    => node %basic()
          timeBegin=  basic%time ()
       end if
    call starFormationHistory_%create                 (node,historyStarFormation,timeBegin)
    call nuclearStarCluster   %starFormationHistorySet(     historyStarFormation          )
    end if
    ! Record that the nuclear star cluster has been initialized.
    call nuclearStarCluster%isInitializedSet(.true. )
    ! Allow the nuclear star cluster to form a new black hole seed.
    call nuclearStarCluster%isCollapsedSet  (.false.)
    return
  end subroutine Node_Component_NSC_Standard_Create

  !![
  <inactiveSetTask>
   <unitName>Node_Component_NSC_Standard_Inactive</unitName>
  </inactiveSetTask>
  !!]
  subroutine Node_Component_NSC_Standard_Inactive(node)
    !!{
    Set Jacobian zero status for properties of {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentNSC, nodeComponentNSCStandard, treeNode
    implicit none
    type (treeNode        ), intent(inout), pointer :: node
    class(nodeComponentNSC)               , pointer :: nuclearStarCluster

    ! Get the nuclear star cluster component.
    nuclearStarCluster => node%NSC()
    ! Check if an standard nuclear star cluster component exists.
    select type (nuclearStarCluster)
    class is (nodeComponentNSCStandard)
       if (inactiveLuminositiesStellar) call nuclearStarCluster%luminositiesStellarInactive()
    end select
    return
  end subroutine Node_Component_NSC_Standard_Inactive

  subroutine satelliteMerger(self,node)
    !!{
    Transfer any standard nuclear star cluster associated with {\normalfont \ttfamily node} to its host halo.
    !!}
    use :: Abundances_Structure            , only : zeroAbundances
    use :: Error                           , only : Error_Report
    use :: Galacticus_Nodes                , only : nodeComponentNSC         , nodeComponentNSCStandard, nodeComponentSpheroid   , nodeComponentDisk               , &
       &                                            treeNode
    use :: Satellite_Merging_Mass_Movements, only : destinationMergerSpheroid, destinationMergerDisk   , destinationMergerUnmoved, enumerationDestinationMergerType
    use :: Satellite_Merging_Remnant_Sizes , only : remnantNoChange
    use :: Stellar_Luminosities_Structure  , only : zeroStellarLuminosities
    implicit none
    class           (*                               ), intent(inout) :: self
    type            (treeNode                        ), intent(inout) :: node
    type            (treeNode                        ), pointer       :: nodeHost
    class           (nodeComponentNSC                ), pointer       :: nuclearStarCluster
    class           (nodeComponentSpheroid           ), pointer       :: spheroidHost           
    class           (nodeComponentDisk               ), pointer       :: diskHost               
    type            (history                         )                :: historyDisk                              , historySpheroid            , &
         &                                                               historyNuclearStarCluster                , history_
    double precision                                                  :: spheroidSpecificAngularMomentum          , diskSpecificAngularMomentum, &
         &                                                               angularMomentumSpecificNuclearStarCluster, massNuclearStarCluster
    type            (enumerationDestinationMergerType)                :: destinationGasSatellite                  , destinationGasHost         , &
         &                                                               destinationStarsHost                     , destinationStarsSatellite
    logical                                                           :: mergerIsMajor
    !$GLC attributes unused :: self

    ! Check that the nuclear star cluster is of the standard class.
    nuclearStarCluster => node%NSC()
    select type (nuclearStarCluster)
    class is (nodeComponentNSCStandard)
       ! Find the node to merge with.
       nodeHost    => node    %mergesWith(                 )
       spheroidHost=> nodeHost%spheroid  (autoCreate=.true.)
       diskHost    => nodeHost%disk      (autoCreate=.true.)
       ! Get specific angular momentum of the nuclear star cluster material.
       if (nuclearStarCluster%massGas()+nuclearStarCluster%massStellar() > 0.0d0) then
          angularMomentumSpecificNuclearStarCluster=+  nuclearStarCluster%angularMomentum() &
            &                                       /(                                      &
            &                                         +nuclearStarCluster%massGas        () &
            &                                         +nuclearStarCluster%massStellar    () &
            &                                        )     
       else
          angularMomentumSpecificNuclearStarCluster=+0.0d0
       end if
       if (spheroidHost      %massGas()+spheroidHost      %massStellar() > 0.0d0) then
          spheroidSpecificAngularMomentum          =+  spheroidHost      %angularMomentum() &
               &                                    /(                                      &
               &                                      +spheroidHost      %massGas        () &
               &                                      +spheroidHost      %massStellar    () &
               &                                     )
       else
          spheroidSpecificAngularMomentum          =+0.0d0
       end if
       if (diskHost          %massGas()+diskHost          %massStellar() > 0.0d0) then
          diskSpecificAngularMomentum              =+  diskHost          %angularMomentum() &
               &                                    /(                                      &
               &                                      +diskHost          %massGas        () &
               &                                      +diskHost          %massStellar    () &
               &                                     )
       else
          diskSpecificAngularMomentum              =+0.0d0
       end if
       ! Get mass movement descriptors.
       call mergerMassMovements_%get(node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
       ! Move gas material within the host if necessary
       select case (destinationGasHost%ID)
       case (destinationMergerDisk%ID)
          call diskHost              %        massGasSet(                                                                &
               &                                         +diskHost          %      massGas                            () &
               &                                         +nuclearStarCluster%      massGas                            () &
               &                                        )
          call diskHost              %  abundancesGasSet(                                                                &
               &                                         +diskHost          %abundancesGas                            () &
               &                                         +nuclearStarCluster%abundancesGas                            () &
               &                                        )
          call diskHost              %angularMomentumSet(                                                                &
               &                                         +diskHost          %angularMomentum                          () &
               &                                         +nuclearStarCluster%      massGas                            () &
               &                                         *                   angularMomentumSpecificNuclearStarCluster   &
               &                                        )
          call nuclearStarCluster    %angularMomentumSet(                                                                &
               &                                         +nuclearStarCluster%angularMomentum                          () &
               &                                         -nuclearStarCluster%      massGas                            () &
               &                                         *                   angularMomentumSpecificNuclearStarCluster   &        
               &                                        )
          call nuclearStarCluster    %        massGasSet(                                                                &
               &                                         +0.0d0                                                          &
               &                                        )
          call nuclearStarCluster    %  abundancesGasSet(                                                                &
               &                                          zeroAbundances                                                 &
               &                                        )     
       case (destinationMergerSpheroid%ID)
          call spheroidHost          %        massGasSet(                                                                &
               &                                         +spheroidHost      %      massGas                            () &
               &                                         +nuclearStarCluster%      massGas                            () &
               &                                        )
          call spheroidHost          %  abundancesGasSet(                                                                &
               &                                         +spheroidHost      %abundancesGas                            () &
               &                                         +nuclearStarCluster%abundancesGas                            () &
               &                                        )
          call spheroidHost          %angularMomentumSet(                                                                &
               &                                         +spheroidHost      %angularMomentum                          () &
               &                                         +nuclearStarCluster%      massGas                            () &
               &                                         *                   angularMomentumSpecificNuclearStarCluster   &
               &                                        )
          call nuclearStarCluster%            massGasSet(                                                                &
               &                                         +0.0d0                                                          &
               &                                        )
          call nuclearStarCluster%      abundancesGasSet(                                                                &
               &                                          zeroAbundances                                                 &
               &                                        )  
       case (destinationMergerUnmoved%ID)
          ! Do nothing.
       case default
          call Error_Report('unrecognized movesTo descriptor'//{introspection:location})
       end select
       ! Move stellar material within the host if necessary
       select case (destinationStarsHost%ID)
       case (destinationMergerDisk%ID)
          call diskHost          %        massStellarSet(                                                                &
               &                                         +diskHost          %        massStellar                      () &
               &                                         +nuclearStarCluster%        massStellar                      () &
               &                                        )
          call diskHost          %  abundancesStellarSet(                                                                &
               &                                         +diskHost          %  abundancesStellar                      () &
               &                                         +nuclearStarCluster%  abundancesStellar                      () &
               &                                        )
          call diskHost          %luminositiesStellarSet(                                                                &
               &                                         +diskHost          %luminositiesStellar                      () &
               &                                         +nuclearStarCluster%luminositiesStellar                      () &
               &                                        )
          call diskHost          %    angularMomentumSet(                                                                &
               &                                         +diskHost          %    angularMomentum                      () &
               &                                         +nuclearStarCluster%        massStellar                      () &
               &                                         *                   angularMomentumSpecificNuclearStarCluster   &
               &                                        )
          call nuclearStarCluster%    angularMomentumSet(                                                                &
               &                                         +nuclearStarCluster%    angularMomentum                      () &
               &                                         -nuclearStarCluster%        massStellar                      () &
               &                                         *                   angularMomentumSpecificNuclearStarCluster   &
               &                                        )
          call nuclearStarCluster%        massStellarSet(                                                                &
               &                                         +0.0d0                                                          &
               &                                        )
          call nuclearStarCluster%  abundancesStellarSet(                                                                &
               &                                          zeroAbundances                                                 &
               &                                        )
          call nuclearStarCluster%luminositiesStellarSet(                                                                &
               &                                          zeroStellarLuminosities                                        &
               &                                        )
          ! Also add stellar properties histories.
          historyDisk              =diskHost          %stellarPropertiesHistory()
          historyNuclearStarCluster=nuclearStarCluster%stellarPropertiesHistory()
          call historyDisk              %interpolatedIncrement      (historyNuclearStarCluster)
          call historyNuclearStarCluster%reset                      (                         )
          call diskHost                 %stellarPropertiesHistorySet(historyDisk              )
          call nuclearStarCluster       %stellarPropertiesHistorySet(historyNuclearStarCluster)
          ! Also add star formation histories.
          historyDisk              =diskHost          %starFormationHistory()
          historyNuclearStarCluster=nuclearStarCluster%starFormationHistory()
          call starFormationHistory_    %move                   (nodeHost,nodeHost,historyDisk,historyNuclearStarCluster)
          call diskHost                 %starFormationHistorySet(                  historyDisk                          )
          call nuclearStarCluster       %starFormationHistorySet(                              historyNuclearStarCluster)
          call historyDisk              %destroy                (                                                       )
          call historyNuclearStarCluster%destroy                (                                                       )              
       case (destinationMergerSpheroid%ID)
          call spheroidHost      %        massStellarSet(                                                                &
               &                                         +spheroidHost      %        massStellar                      () &
               &                                         +nuclearStarCluster%        massStellar                      () &
               &                                        )
          call spheroidHost      %    angularMomentumSet(                                                                &
               &                                         +spheroidHost      %    angularMomentum                      () &
               &                                         +nuclearStarCluster%        massStellar                      () &
               &                                         *                   angularMomentumSpecificNuclearStarCluster   &
               &                                        )
          call spheroidHost      %  abundancesStellarSet(                                                                &
               &                                         +spheroidHost      %  abundancesStellar                      () &
               &                                         +nuclearStarCluster%  abundancesStellar                      () &
               &                                        )
          call spheroidHost      %luminositiesStellarSet(                                                                &
               &                                         +spheroidHost      %luminositiesStellar                      () &
               &                                         +nuclearStarCluster%luminositiesStellar                      () &
               &                                        )
          call nuclearStarCluster%        massStellarSet(                                                                &
               &                                         +0.0d0                                                          &
               &                                        )
          call nuclearStarCluster%  abundancesStellarSet(                                                                &
               &                                          zeroAbundances                                                 &
               &                                        )
          call nuclearStarCluster%luminositiesStellarSet(                                                                &
               &                                          zeroStellarLuminosities                                        &
               &                                        )
          ! Also add stellar properties histories.
          historyNuclearStarCluster=nuclearStarCluster%stellarPropertiesHistory()
          historySpheroid          =spheroidHost      %stellarPropertiesHistory()
          call historySpheroid          %interpolatedIncrement      (historyNuclearStarCluster)
          call historyNuclearStarCluster%reset                      (                         )
          call spheroidHost             %stellarPropertiesHistorySet(historySpheroid          )
          call nuclearStarCluster       %stellarPropertiesHistorySet(historyNuclearStarCluster)
          ! Also add star formation histories.
          historyNuclearStarCluster=nuclearStarCluster%starFormationHistory   ()
          historySpheroid          =spheroidHost      %starFormationHistory   ()
          call starFormationHistory_    %move                   (nodeHost,nodeHost,historySpheroid,historyNuclearStarCluster)
          call spheroidHost             %starFormationHistorySet(                  historySpheroid                          )
          call nuclearStarCluster       %starFormationHistorySet(                                  historyNuclearStarCluster)
          call historyNuclearStarCluster%destroy                (                                                           )
          call historySpheroid          %destroy                (                                                           )
          historyNuclearStarCluster=nuclearStarCluster%starFormationHistory   ()
          historySpheroid          =spheroidHost      %starFormationHistory   ()
       case (destinationMergerUnmoved%ID)
          ! Do nothing
       case default
          call Error_Report('unrecognized movesTo descriptor'//{introspection:location})
       end select
       ! Get specific angular momentum of the nuclear star cluster material.
       massNuclearStarCluster=+nuclearStarCluster%massGas    () &
            &                 +nuclearStarCluster%massStellar()
       if (massNuclearStarCluster > 0.0d0) then
          angularMomentumSpecificNuclearStarCluster=nuclearStarCluster%angularMomentum()/massNuclearStarCluster
          ! Move the gas component of the standard nuclear star cluster to the host.
          select case (destinationGasSatellite%ID)
          case (destinationMergerDisk%ID)
             call diskHost    %       massGasSet     (                                                                &
                  &                                   +diskHost          %massGas                                  () &
                  &                                   +nuclearStarCluster%massGas                                  () &
                  &                                  )
             call diskHost    %  abundancesGasSet    (                                                                &
                  &                                   +diskHost          %abundancesGas                            () &
                  &                                   +nuclearStarCluster%abundancesGas                            () &
                  &                                  )
             call diskHost    %angularMomentumSet    (                                                                &
                  &                                   +diskHost          %angularMomentum                          () &
                  &                                   +nuclearStarCluster%massGas                                  () &
                  &                                   *                   angularMomentumSpecificNuclearStarCluster   &
                  &                                  )
          case (destinationMergerSpheroid%ID)
             call spheroidHost%massGasSet            (                                                                &
                  &                                   +spheroidHost      %massGas                                  () &
                  &                                   +nuclearStarCluster%massGas                                  () &
                  &                                  )
             call spheroidHost%abundancesGasSet      (                                                                &
                  &                                   +spheroidHost      %abundancesGas                            () &
                  &                                   +nuclearStarCluster%abundancesGas                            () &
                  &                                  )
             call spheroidHost%angularMomentumSet    (                                                                &
                  &                                   +spheroidHost      %angularMomentum                          () &
                  &                                   +nuclearStarCluster%massGas                                  () &
                  &                                   *                   angularMomentumSpecificNuclearStarCluster   &
                  &                                  )
          case default
             call Error_Report('unrecognized movesTo descriptor'//{introspection:location})
          end select
          call nuclearStarCluster%      massGasSet(         0.0d0)
          call nuclearStarCluster%abundancesGasSet(zeroAbundances)
          ! Move the stellar component of the standard nuclear star cluster to the host.
          select case (destinationStarsSatellite%ID)
          case (destinationMergerDisk%ID)
             call diskHost    %massStellarSet        (                                                                &
                  &                                   +diskHost          %        massStellar                      () &
                  &                                   +nuclearStarCluster%        massStellar                      () &
                  &                                  )
             call diskHost    %abundancesStellarSet  (                                                                &
                  &                                   +diskHost          %  abundancesStellar                      () &
                  &                                   +nuclearStarCluster%  abundancesStellar                      () &
                  &                                  )
             call diskHost    %luminositiesStellarSet(                                                                &
                  &                                   +diskHost          %luminositiesStellar                      () &
                  &                                   +nuclearStarCluster%luminositiesStellar                      () &
                  &                                  )
             call diskHost    %angularMomentumSet    (                                                                &
                  &                                   +diskHost          %angularMomentum                          () &
                  &                                   +nuclearStarCluster%        massStellar                      () &
                  &                                   *                   angularMomentumSpecificNuclearStarCluster   & 
                  &                                  )
             ! Also add stellar properties histories.
             historyNuclearStarCluster=nuclearStarCluster%stellarPropertiesHistory()
             history_                 =diskHost          %stellarPropertiesHistory()
             call history_                 %interpolatedIncrement      (historyNuclearStarCluster)
             call historyNuclearStarCluster%reset                      (                         )
             call diskHost                 %stellarPropertiesHistorySet(history_                 )
             call nuclearStarCluster       %stellarPropertiesHistorySet(historyNuclearStarCluster)
             ! Also add star formation histories.
             historyNuclearStarCluster=nuclearStarCluster%starFormationHistory    ()
             history_                 =diskHost          %starFormationHistory    ()
             call StarFormationHistory_    %move                   (nodeHost,node,history_,historyNuclearStarCluster)
             call diskHost                 %starFormationHistorySet(              history_                          )
             call nuclearStarCluster       %starFormationHistorySet(                       historyNuclearStarCluster)
             call history_                 %destroy                (                                                )
             call historyNuclearStarCluster%destroy                (                                                )
          case (destinationMergerSpheroid%ID)
             call spheroidHost%massStellarSet        (                                          &
                  &                                   +spheroidHost      %massStellar        () &
                  &                                   +nuclearStarCluster%massStellar        () &
                  &                                  )
             call spheroidHost%abundancesStellarSet  (                                          &
                  &                                   +spheroidHost      %abundancesStellar  () &
                  &                                   +nuclearStarCluster%abundancesStellar  () &
                  &                                  )
             call spheroidHost%luminositiesStellarSet(                                          &
                  &                                   +spheroidHost      %luminositiesStellar() &
                  &                                   +nuclearStarCluster%luminositiesStellar() &
                  &                                  )
             ! Also add stellar properties histories.
             historyNuclearStarCluster=nuclearStarCluster%stellarPropertiesHistory()
             history_                 =spheroidHost      %stellarPropertiesHistory()
             call history_                 %interpolatedIncrement      (historyNuclearStarCluster)
             call historyNuclearStarCluster%reset                      (                         )
             call spheroidHost             %stellarPropertiesHistorySet(history_                 )
             call nuclearStarCluster       %stellarPropertiesHistorySet(historyNuclearStarCluster)
             ! Also add star formation histories.
             historyNuclearStarCluster=nuclearStarCluster%starFormationHistory    ()
             history_                 =spheroidHost      %starFormationHistory    ()
             call starFormationHistory_    %move                   (nodeHost,node,history_,historyNuclearStarCluster)
             call spheroidHost             %starFormationHistorySet(              history_                          )
             call nuclearStarCluster       %starFormationHistorySet(                       historyNuclearStarCluster)
             call history_                 %destroy                (                                                )
             call historyNuclearStarCluster%destroy                (                                                )
          case default
             call Error_Report('unrecognized movesTo descriptor'//{introspection:location})
          end select
          call nuclearStarCluster%        massStellarSet(                  0.0d0)
          call nuclearStarCluster%  abundancesStellarSet(         zeroAbundances)
          call nuclearStarCluster%luminositiesStellarSet(zeroStellarLuminosities)
          call nuclearStarCluster%    angularMomentumSet(                  0.0d0)
          call nuclearStarCluster%        isCollapsedSet(                .false.)
       end if
    end select
    return
  end subroutine satelliteMerger

  subroutine Node_Component_NSC_Standard_Star_Formation_History_Extend(node,timeEnd)
    !!{
    Extend the range of a star formation history in a standard nuclear star cluster component for {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentNSC, treeNode
    implicit none
    type            (treeNode        ), intent(inout), target   :: node
    double precision                  , intent(in   ), optional :: timeEnd
    class           (nodeComponentNSC)               , pointer  :: nuclearStarCluster
    type            (history         )                          :: historyStarFormation
    !$GLC attributes unused :: timeEnd
    
    ! Get the nuclear star cluster component.
    nuclearStarCluster => node%NSC()
    ! Extend the range as necessary.
    historyStarFormation=nuclearStarCluster%starFormationHistory()
    call starFormationHistory_%extend                 (historyStarFormation,starFormationHistoryTemplate)
    call nuclearStarCluster   %starFormationHistorySet(historyStarFormation                             )
    return
  end subroutine Node_Component_NSC_Standard_Star_Formation_History_Extend

  subroutine Node_Component_NSC_Standard_Stellar_Prprts_History_Extend(node,timeEnd)
    !!{
    Extend the range of a stellar properties history in a standard nuclear star cluster component for {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentNSC, treeNode
    implicit none
    type            (treeNode        ), intent(inout), target   :: node
    double precision                  , intent(in   ), optional :: timeEnd
    class           (nodeComponentNSC)               , pointer  :: nuclearStarCluster
    type            (history         )                          :: stellarPropertiesHistory
    !$GLC attributes unused :: timeEnd
    
    ! Get the nuclear star cluster component.
    nuclearStarCluster => node%NSC()
    ! Extend the range as necessary.
    stellarPropertiesHistory=nuclearStarCluster%stellarPropertiesHistory()
    call stellarPropertiesHistory%extend                     (times=stellarPropertiesHistoryTemplate)
    call nuclearStarCluster      %stellarPropertiesHistorySet(      stellarPropertiesHistory        )
    return
  end subroutine Node_Component_NSC_Standard_Stellar_Prprts_History_Extend

  subroutine mergerTreeExtraOutput(self,node,iOutput,treeIndex,nodePassesFilter,treeLock)
   !!{
     Update the star formation history after an output time is reached.
     !!}
     use            :: Galacticus_Nodes, only : defaultNSCComponent, nodeComponentNSC, nodeComponentNSCStandard, treeNode
     use            :: Histories       , only : history
     use, intrinsic :: ISO_C_Binding   , only : c_size_t
     use            :: Kind_Numbers    , only : kind_int8
     use            :: Locks           , only : ompLock
     implicit none
     class  (*               ), intent(inout)          :: self
     type   (treeNode        ), intent(inout)          :: node
     integer(c_size_t        ), intent(in   )          :: iOutput
     integer(kind=kind_int8  ), intent(in   )          :: treeIndex
     logical                  , intent(in   )          :: nodePassesFilter
     type   (ompLock         ), intent(inout)          :: treeLock
     class  (nodeComponentNSC)               , pointer :: nuclearStarCluster
     type   (history         )                         :: historyStarFormation
     !$GLC attributes unused :: self, treeIndex, nodePassesFilter, treeLock

     ! Check if we are the default method.
     if (.not.defaultNSCComponent%standardIsActive()) return
     ! Output the star formation history if a nuclear star cluster exists for this component.
     nuclearStarCluster   => node%NSC                               ()
     historyStarFormation =  nuclearStarCluster%starFormationHistory()
     call starFormationHistory_%update(node,iOutput,historyStarFormation)
     ! Update the star formation history only if a nuclear star cluster exists.
     select type (nuclearStarCluster)
     class is (nodeComponentNSCStandard)
        call nuclearStarCluster%starFormationHistorySet(historyStarFormation)
     end select
     return
   end subroutine mergerTreeExtraOutput

  !![
  <stateStoreTask>
   <unitName>Node_Component_NSC_Standard_State_Store</unitName>
  </stateStoreTask>
  !!]
  subroutine Node_Component_NSC_Standard_State_Store(stateFile,gslStateFile,stateOperationID)
    !!{
    Write the state to file.
    !!}
    use            :: Display                         , only : displayMessage          , verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding                   , only : c_ptr                   , c_size_t
    use            :: Node_Component_NSC_Standard_Data, only : massDistributionStellar_, massDistributionGas_, kinematicDistribution_
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Storing state for: componentNSC -> standard',verbosity=verbosityLevelInfo)
    !![
    <stateStore variables="massDistributionStellar_ massDistributionGas_ kinematicDistribution_ stellarPopulationProperties_ darkMatterHaloScale_ starFormationHistory_ mergerMassMovements_ mergerRemnantSize_"/>
    !!]
    return
  end subroutine Node_Component_NSC_Standard_State_Store

  !![
  <stateRetrieveTask>
   <unitName>Node_Component_NSC_Standard_State_Retrieve</unitName>
  </stateRetrieveTask>
  !!]
  subroutine Node_Component_NSC_Standard_State_Retrieve(stateFile,gslStateFile,stateOperationID)
    !!{
    Retrieve the tabulation state from the file.
    !!}
    use            :: Display                         , only : displayMessage          , verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding                   , only : c_ptr                   , c_size_t
    use            :: Node_Component_NSC_Standard_Data, only : massDistributionStellar_, massDistributionGas_, kinematicDistribution_
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Retrieving state for: componentNSC -> standard',verbosity=verbosityLevelInfo)
    !![
    <stateRestore variables="massDistributionStellar_ massDistributionGas_ kinematicDistribution_ stellarPopulationProperties_ darkMatterHaloScale_ starFormationHistory_ mergerMassMovements_ mergerRemnantSize_"/>
    !!]
    return
  end subroutine Node_Component_NSC_Standard_State_Retrieve

end module Node_Component_NSC_Standard
