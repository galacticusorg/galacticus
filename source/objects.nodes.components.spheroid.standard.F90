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
Contains a module which implements the standard spheroid component.
!!}

module Node_Component_Spheroid_Standard
  !!{
  Implements the standard spheroid component.
  !!}
  use :: Dark_Matter_Halo_Scales         , only : darkMatterHaloScaleClass
  use :: Histories                       , only : history
  use :: Satellite_Merging_Mass_Movements, only : mergerMassMovementsClass
  use :: Satellite_Merging_Remnant_Sizes , only : mergerRemnantSizeClass
  use :: Star_Formation_Histories        , only : starFormationHistory            , starFormationHistoryClass
  use :: Stellar_Population_Properties   , only : stellarPopulationPropertiesClass
  implicit none
  private
  public :: Node_Component_Spheroid_Standard_Initialize         , Node_Component_Spheroid_Standard_Scale_Set                 , &
       &    Node_Component_Spheroid_Standard_Pre_Evolve         , Node_Component_Spheroid_Standard_Radius_Solver_Plausibility, &
       &    Node_Component_Spheroid_Standard_Thread_Uninitialize, Node_Component_Spheroid_Standard_Thread_Initialize         , &
       &    Node_Component_Spheroid_Standard_State_Store        , Node_Component_Spheroid_Standard_State_Retrieve            , &
       &    Node_Component_Spheroid_Standard_Inactive           , Node_Component_Spheroid_Standard_Post_Step                 , &
       &    Node_Component_Spheroid_Standard_Radius_Solver

  !![
  <component>
   <class>spheroid</class>
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
      <name>massStellar</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" isNonNegative="true"/>
      <output unitsInSI="massSolar" comment="Mass of stars in the standard spheroid."/>
    </property>
    <property>
      <name>massStellarFormed</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" isNonNegative="true" />
    </property>
    <property>
      <name>abundancesStellar</name>
      <type>abundances</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" isNonNegative="true" />
      <output unitsInSI="massSolar" comment="Mass of metals in the stellar phase of the standard spheroid."/>
    </property>
    <property>
      <name>massGas</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" isNonNegative="true"/>
      <output unitsInSI="massSolar" comment="Mass of gas in the standard spheroid."/>
    </property>
    <property>
      <name>abundancesGas</name>
      <type>abundances</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" isNonNegative="true" />
      <output unitsInSI="massSolar" comment="Mass of metals in the gas phase of the standard spheroid."/>
    </property>
    <property>
      <name>angularMomentum</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" isNonNegative="true"/>
      <output unitsInSI="massSolar*megaParsec*kilo" comment="Angular momentum of the standard spheroid."/>
    </property>
    <property>
      <name>radius</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
      <output unitsInSI="megaParsec" comment="Scale length of the standard spheroid."/>
    </property>
    <property>
      <name>halfMassRadius</name>
      <attributes isSettable="false" isGettable="true" isEvolvable="false" isVirtual="true" />
      <type>double</type>
      <rank>0</rank>
      <getFunction>Node_Component_Spheroid_Standard_Half_Mass_Radius</getFunction>
    </property>
    <property>
      <name>velocity</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
      <output unitsInSI="kilo" comment="Circular velocity at the scale length of the standard spheroid."/>
    </property>
    <property>
      <name>luminositiesStellar</name>
      <type>stellarLuminosities</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" isNonNegative="true" />
      <output unitsInSI="luminosityZeroPointAB" comment="Luminosity of spheroid stars."/>
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
      <attributes isSettable="true" isGettable="true" isEvolvable="true" isDeferred="rate" createIfNeeded="true" isNonNegative="true" />
    </property>
    <property>
      <name>massGasSink</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="false" isGettable="false" isEvolvable="true" isDeferred="rate" isVirtual="true" />
    </property>
    <property>
      <name>energyGasInput</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="false" isGettable="false" isEvolvable="true" isDeferred="rate" isVirtual="true" />
    </property>
   </properties>
   <bindings>
    <binding method="massDistribution" function="Node_Component_Spheroid_Standard_Mass_Distribution"/>
    <binding method="massBaryonic"     function="Node_Component_Spheroid_Standard_Mass_Baryonic"    />
   </bindings>
   <functions>objects.nodes.components.spheroid.standard.bound_functions.inc</functions>
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
  double precision, allocatable, dimension(:) :: starFormationHistoryTemplate   , stellarPropertiesHistoryTemplate
  !$omp threadprivate(starFormationHistoryTemplate,stellarPropertiesHistoryTemplate)

  ! Parameters controlling the physical implementation.
  double precision                            :: efficiencyEnergeticOutflow     , toleranceAbsoluteMass           , &
       &                                         toleranceRelativeMetallicity
  logical                                     :: inactiveLuminositiesStellar    , postStepZeroNegativeMasses

  ! Spheroid structural parameters.
  double precision                            :: ratioAngularMomentumScaleRadius

  ! Minimum absolute scales for physically plausible spheroids.
  double precision, parameter                 :: radiusMinimum                       =1.0d-12
  double precision, parameter                 :: massMinimum                         =1.0d-06
  double precision, parameter                 :: angularMomentumMinimum              =1.0d-20

  ! A threadprivate object used to track to which thread events are attached.
  integer :: thread
  !$omp threadprivate(thread)

contains

  !![
  <nodeComponentInitializationTask>
   <unitName>Node_Component_Spheroid_Standard_Initialize</unitName>
  </nodeComponentInitializationTask>
  !!]
  subroutine Node_Component_Spheroid_Standard_Initialize(parameters)
    !!{
    Initializes the tree node standard spheroid methods module.
    !!}
    use :: Abundances_Structure, only : Abundances_Property_Count
    use :: Galacticus_Nodes    , only : defaultSpheroidComponent , nodeComponentSpheroidStandard
    use :: Input_Parameters    , only : inputParameter           , inputParameters
    implicit none
    type(inputParameters              ), intent(inout) :: parameters
    type(nodeComponentSpheroidStandard)                :: spheroidStandardComponent
    type(inputParameters              )                :: subParameters

    ! Initialize the module if necessary.
    if (defaultSpheroidComponent%standardIsActive()) then

       ! Get number of abundance properties.
       abundancesCount  =Abundances_Property_Count            ()

       ! Bind deferred functions.
       call spheroidStandardComponent%          energyGasInputRateFunction(Node_Component_Spheroid_Standard_Energy_Gas_Input_Rate      )
       call spheroidStandardComponent%             massGasSinkRateFunction(Node_Component_Spheroid_Standard_Mass_Gas_Sink_Rate         )
       call spheroidStandardComponent%    starFormationHistoryRateFunction(Node_Component_Spheroid_Standard_Star_Formation_History_Rate)
       call spheroidStandardComponent%stellarPropertiesHistoryRateFunction(Node_Component_Spheroid_Standard_Stellar_Prprts_History_Rate)
       call spheroidStandardComponent%                   createFunctionSet(Node_Component_Spheroid_Standard_Initializor                )

       ! Find our parameters.
       subParameters=parameters%subParameters('componentSpheroid')
       ! Read parameters controlling the physical implementation.
       !![
       <inputParameter>
         <name>efficiencyEnergeticOutflow</name>
         <defaultValue>1.0d-2</defaultValue>
         <description>The proportionality factor relating mass outflow rate from the spheroid to the energy input rate divided by $V_\mathrm{spheroid}^2$.</description>
         <source>subParameters</source>
       </inputParameter>
       <inputParameter>
         <name>toleranceRelativeMetallicity</name>
         <defaultValue>1.0d-4</defaultValue>
         <description>The metallicity tolerance for ODE solution.</description>
         <source>subParameters</source>
       </inputParameter>
       <inputParameter>
         <name>toleranceAbsoluteMass</name>
         <defaultValue>1.0d-6</defaultValue>
         <description>The mass tolerance used to judge whether the spheroid is physically plausible.</description>
         <source>subParameters</source>
       </inputParameter>
       <inputParameter>
         <name>inactiveLuminositiesStellar</name>
         <defaultValue>.false.</defaultValue>
         <description>Specifies whether or not spheroid stellar luminosities are inactive properties (i.e. do not appear in any ODE being solved).</description>
         <source>subParameters</source>
       </inputParameter>
       <inputParameter>
         <name>postStepZeroNegativeMasses</name>
         <defaultValue>.true.</defaultValue>
         <description>If true, negative masses will be zeroed after each ODE step. Note that this can lead to non-conservation of mass.</description>
         <source>subParameters</source>
       </inputParameter>
       !!]
    end if
    return
  end subroutine Node_Component_Spheroid_Standard_Initialize

  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Node_Component_Spheroid_Standard_Thread_Initialize</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Node_Component_Spheroid_Standard_Thread_Initialize(parameters)
    !!{
    Initializes the standard spheroid module for each thread.
    !!}
    use :: Events_Hooks                         , only : dependencyDirectionAfter , dependencyRegEx           , openMPThreadBindingAtLevel, postEvolveEvent, &
          &                                              satelliteMergerEvent     , mergerTreeExtraOutputEvent
    use :: Error                                , only : Error_Report
    use :: Galacticus_Nodes                     , only : defaultSpheroidComponent
    use :: Input_Parameters                     , only : inputParameter           , inputParameters
    use :: Mass_Distributions                   , only : massDistributionSpherical, kinematicsDistributionLocal
    use :: Node_Component_Spheroid_Standard_Data, only : massDistributionStellar_ , massDistributionGas_       , kinematicDistribution_
    use :: Galactic_Structure_Options           , only : componentTypeSpheroid    , massTypeStellar            , massTypeGaseous
    implicit none
    type            (inputParameters), intent(inout) :: parameters
    type            (dependencyRegEx), dimension(3)  :: dependencies
    logical                                          :: densityMoment2IsInfinite                , densityMoment3IsInfinite
    double precision                                 :: massDistributionSpheroidDensityMomentum2, massDistributionSpheroidDensityMomentum3, &
         &                                              ratioAngularMomentumScaleRadiusDefault
    type            (inputParameters)                :: subParameters

    ! Check if this implementation is selected. If so, initialize the mass distribution.
    if (defaultSpheroidComponent%standardIsActive()) then
       dependencies(1)=dependencyRegEx(dependencyDirectionAfter,'^remnantStructure:')
       dependencies(2)=dependencyRegEx(dependencyDirectionAfter,'^preAnalysis:'     )
       dependencies(3)=dependencyRegEx(dependencyDirectionAfter,'^nodeComponentDisk')
       call satelliteMergerEvent      %attach(thread,satelliteMerger      ,openMPThreadBindingAtLevel,label='nodeComponentSpheroidStandard',dependencies=dependencies)
       call mergerTreeExtraOutputEvent%attach(thread,mergerTreeExtraOutput,openMPThreadBindingAtLevel,label='nodeComponentSpheroidStandard'                          )
       call postEvolveEvent           %attach(thread,postEvolve           ,openMPThreadBindingAtLevel,label='nodeComponentSpheroidStandard'                          ) 
       ! Find our parameters.
       subParameters=parameters%subParameters('componentSpheroid')
       !![
       <objectBuilder class="stellarPopulationProperties"                                          name="stellarPopulationProperties_" source="subParameters"                    />
       <objectBuilder class="darkMatterHaloScale"                                                  name="darkMatterHaloScale_"         source="subParameters"                    />
       <objectBuilder class="starFormationHistory"                                                 name="starFormationHistory_"        source="subParameters"                    />
       <objectBuilder class="mergerMassMovements"                                                  name="mergerMassMovements_"         source="subParameters"                    />
       <objectBuilder class="mergerRemnantSize"                                                    name="mergerRemnantSize_"           source="subParameters"                    />
       <objectBuilder class="massDistribution"            parameterName="massDistributionSpheroid" name="massDistributionStellar_"     source="subParameters" threadPrivate="yes" >
        <default>
         <massDistributionSpheroid value="hernquist">
          <dimensionless value="true"/>
         </massDistributionSpheroid>
        </default>
       </objectBuilder>
       !!]
       ! Validate the disk mass distribution.
       select type (massDistributionStellar_)
       class is (massDistributionSpherical)
          ! The spheroid mass distribution must have spherical symmetry. So, this is acceptable.        
       class default
          call Error_Report('only spehrically symmetric mass distributions are allowed'//{introspection:location})
       end select
       if (.not.massDistributionStellar_%isDimensionless()) call Error_Report('spheroid mass distribution must be dimensionless'//{introspection:location})
       ! Duplicate the dimensionless mass distribution to use for the gas component, and set component and mass types in both.
       !$omp critical(spheroidStandardDeepCopy)
       allocate(massDistributionGas_,mold=massDistributionStellar_)
       !![
       <deepCopyReset variables="massDistributionStellar_"/>
       <deepCopy source="massDistributionStellar_" destination="massDistributionGas_"/>
       <deepCopyFinalize variables="massDistributionGas_"/>  
       !!]
       !$omp end critical(spheroidStandardDeepCopy)
       call massDistributionStellar_%setTypes(componentTypeSpheroid,massTypeStellar)
       call massDistributionGas_    %setTypes(componentTypeSpheroid,massTypeGaseous)
       ! Construct the kinematic distribution.
       allocate(kinematicDistribution_)
       !![
       <referenceConstruct object="kinematicDistribution_" constructor="kinematicsDistributionLocal(alpha=1.0d0/sqrt(2.0d0))"/>
       !!]
       ! Determine the specific angular momentum at the scale radius in units of the mean specific angular
       ! momentum of the spheroid. This is equal to the ratio of the 2nd to 3rd radial moments of the density
       ! distribution (assuming a flat rotation curve).
       massDistributionSpheroidDensityMomentum2=massDistributionStellar_%densityRadialMoment(2.0d0,isInfinite=densityMoment2IsInfinite)
       massDistributionSpheroidDensityMomentum3=massDistributionStellar_%densityRadialMoment(3.0d0,isInfinite=densityMoment3IsInfinite)
       if (densityMoment2IsInfinite.or.densityMoment3IsInfinite) then
          ! One of the moments is infinite, so we can not compute the appropriate ratio. Simply assume a value
          ! of 0.5 as a default.
          ratioAngularMomentumScaleRadiusDefault=+0.5d0
       else
          ! Moments are well-defined, so compute their ratio.
          ratioAngularMomentumScaleRadiusDefault=+massDistributionSpheroidDensityMomentum2 &
               &                                 /massDistributionSpheroidDensityMomentum3
       end if
       !$omp critical (spheroidStandardInitializeAngularMomentum)
       !![
       <inputParameter>
         <name>ratioAngularMomentumScaleRadius</name>
         <defaultSource>($I_2/I_3$ where $I_n=\int_0^\infty \rho(r) r^n \mathrm{d}r$, where $\rho(r)$ is the spheroid density profile, unless either $I_2$ or $I_3$ is infinite, in which case a default of $1/2$ is used instead.)</defaultSource>
         <defaultValue>ratioAngularMomentumScaleRadiusDefault</defaultValue>
         <description>The assumed ratio of the specific angular momentum at the scale radius to the mean specific angular momentum of the standard spheroid component.</description>
         <source>subParameters</source>
       </inputParameter>
       !!]
       !$omp end critical (spheroidStandardInitializeAngularMomentum)
    end if
    return
  end subroutine Node_Component_Spheroid_Standard_Thread_Initialize

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Node_Component_Spheroid_Standard_Thread_Uninitialize</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Node_Component_Spheroid_Standard_Thread_Uninitialize()
    !!{
    Uninitializes the standard spheroid module for each thread.
    !!}
    use :: Events_Hooks                         , only : postEvolveEvent         , satelliteMergerEvent, mergerTreeExtraOutputEvent
    use :: Galacticus_Nodes                     , only : defaultSpheroidComponent
    use :: Node_Component_Spheroid_Standard_Data, only : massDistributionStellar_, massDistributionGas_, kinematicDistribution_
    implicit none

    if (defaultSpheroidComponent%standardIsActive()) then
       if (postEvolveEvent           %isAttached(thread,postEvolve           )) call postEvolveEvent           %detach(thread,postEvolve           )
       if (satelliteMergerEvent      %isAttached(thread,satelliteMerger      )) call satelliteMergerEvent      %detach(thread,satelliteMerger      )
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
  end subroutine Node_Component_Spheroid_Standard_Thread_Uninitialize

  !![
  <preEvolveTask>
  <unitName>Node_Component_Spheroid_Standard_Pre_Evolve</unitName>
  </preEvolveTask>
  !!]
  subroutine Node_Component_Spheroid_Standard_Pre_Evolve(node)
    !!{
    Ensure the spheroid has been initialized.
    !!}
    use :: Galacticus_Nodes, only : defaultSpheroidComponent, nodeComponentSpheroid, nodeComponentSpheroidStandard, treeNode
    implicit none
    type (treeNode             ), intent(inout), pointer :: node
    class(nodeComponentSpheroid)               , pointer :: spheroid

    ! Check if we are the default method.
    if (.not.defaultSpheroidComponent%standardIsActive()) return
    ! Get the spheroid component.
    spheroid => node%spheroid()
    ! Check if a standard spheroid component exists.
    select type (spheroid)
    class is (nodeComponentSpheroidStandard)
       ! Initialize the spheroid.
       call Node_Component_Spheroid_Standard_Initializor(spheroid)
    end select
    return
  end subroutine Node_Component_Spheroid_Standard_Pre_Evolve

  subroutine postEvolve(self,node)
    !!{
    Trim histories attached to the spheroid.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentSpheroid, nodeComponentSpheroidStandard, treeNode
    use :: Histories       , only : history
    implicit none
    class(*                    ), intent(inout) :: self
    type (treeNode             ), intent(inout) :: node
    class(nodeComponentSpheroid), pointer       :: spheroid
    class(nodeComponentBasic   ), pointer       :: basic
    type (history              )                :: stellarPropertiesHistory
    !$GLC attributes unused :: self

    ! Get the spheroid component.
    spheroid => node%spheroid()
    ! Check if an exponential spheroid component exists.
    select type (spheroid)
    class is (nodeComponentSpheroidStandard)
       ! Trim the stellar populations properties future history.
       basic => node%basic()
       stellarPropertiesHistory=spheroid%stellarPropertiesHistory()
       call stellarPropertiesHistory%trim(basic%time())
       call spheroid%stellarPropertiesHistorySet(stellarPropertiesHistory)
    end select
    return
  end subroutine postEvolve

  !![
  <postStepTask>
    <unitName>Node_Component_Spheroid_Standard_Post_Step</unitName>
  </postStepTask>
  !!]
  subroutine Node_Component_Spheroid_Standard_Post_Step(node,status)
    !!{
    Trim histories attached to the spheroid.
    !!}
    use :: Abundances_Structure          , only : abs                     , zeroAbundances
    use :: Display                       , only : displayMessage          , verbosityLevelWarn
    use :: Galacticus_Nodes              , only : defaultSpheroidComponent, nodeComponentSpheroid, nodeComponentSpheroidStandard, treeNode
    use :: Histories                     , only : history
    use :: Interface_GSL                 , only : GSL_Success             , GSL_Continue
    use :: ISO_Varying_String            , only : assignment(=)           , operator(//)         , varying_string
    use :: Stellar_Luminosities_Structure, only : abs                     , stellarLuminosities
    use :: String_Handling               , only : operator(//)
    implicit none
    type            (treeNode             ), intent(inout), pointer :: node
    integer                                , intent(inout)          :: status
    class           (nodeComponentSpheroid)               , pointer :: spheroid
    double precision                       , parameter              :: angularMomentumTolerance=1.0d-2
    double precision                       , parameter              :: massTolerance           =1.0d+0
    double precision                       , save                   :: fractionalErrorMaximum  =0.0d+0
    double precision                                                :: fractionalError                , specificAngularMomentum, &
         &                                                             massSpheroid
    character       (len=20               )                         :: valueString
    type            (varying_string       ), save                   :: message
    !$omp threadprivate(message)
    type            (stellarLuminosities  ), save                   :: luminositiesStellar
    !$omp threadprivate(luminositiesStellar)
    type            (history              ), save                   :: historyStellar
    !$omp threadprivate(historyStellar)

    ! Return immediately if this class is not in use or if masses are not to be zeroed.
    if (.not.defaultSpheroidComponent%standardIsActive().or..not.postStepZeroNegativeMasses) return
    ! Get the spheroid component.
    spheroid => node%spheroid()
    ! Check if an standard spheroid component exists.
    select type (spheroid)
    class is (nodeComponentSpheroidStandard)
       ! Note that "status" is not set to failure as these changes in state of the spheroid should not change any calculation of
       ! differential evolution rates as a negative gas/stellar mass was unphysical anyway.
       !
       ! Trap negative gas masses.
       if (spheroid%massGas() < 0.0d0) then
          ! Check if this exceeds the maximum previously recorded error.
          fractionalError=   abs(spheroid%massGas    ()) &
               &          /(                             &
               &             abs(spheroid%massGas    ()) &
               &            +abs(spheroid%massStellar()) &
               &           )
          !$omp critical (Standard_Spheroid_Post_Evolve_Check)
          if (fractionalError > fractionalErrorMaximum) then
             ! Report a warning.
             message='Warning: spheroid has negative gas mass (fractional error exceeds any previously reported):'//char(10)
             message=message//'  Node index            = '//node%index() //char(10)
             write (valueString,'(e12.6)') spheroid%massGas    ()
             message=message//'  Spheroid gas mass     = '//trim(valueString)//char(10)
             write (valueString,'(e12.6)') spheroid%massStellar()
             message=message//'  Spheroid stellar mass = '//trim(valueString)//char(10)
             write (valueString,'(e12.6)') fractionalError
             message=message//'  Error measure         = '//trim(valueString)//char(10)
             if (fractionalErrorMaximum == 0.0d0) then
                ! This is the first time this warning has been issued, so give some extra information.
                message=message//'  Gas mass will be reset to zero (in future cases also).'//char(10)
                message=message//'  Future cases will be reported only when they exceed the previous maximum error measure.'//char(10)
                message=message//'  Negative masses are due to numerically inaccuracy in the ODE solutions.'//char(10)
                message=message//'  If significant, consider using a higher tolerance in the ODE solver.'
             end if
             call displayMessage(message,verbosityLevelWarn)
             ! Store the new maximum fractional error.
             fractionalErrorMaximum=fractionalError
          end if
          !$omp end critical (Standard_Spheroid_Post_Evolve_Check)
          ! Get the specific angular momentum of the spheroid material
          massSpheroid= spheroid%massGas    () &
               &       +spheroid%massStellar()
          if (massSpheroid == 0.0d0) then
             specificAngularMomentum=0.0d0
             call spheroid%        massStellarSet(                  0.0d0)
             call spheroid%  abundancesStellarSet(         zeroAbundances)
             historyStellar=spheroid%starFormationHistory    ()
             call historyStellar%reset()
             call spheroid      %starFormationHistorySet   (historyStellar)
             historyStellar=spheroid%stellarPropertiesHistory()
             call historyStellar%reset()
             call spheroid      %stellarPropertiesHistorySet(historyStellar)
             ! We need to reset the stellar luminosities to zero. We can't simply use the "zeroStellarLuminosities" instance since
             ! our luminosities may have been truncated. If we were to use "zeroStellarLuminosities" then the number of stellar
             ! luminosities associated with the spheroid would change - but we are in the middle of differential evolution here and we
             ! cannot change the number of evolvable properties as doing so will lead to invalid memory accesses during
             ! deserialization of properties from the ODE solver.
             call luminositiesStellar%destroy()
             luminositiesStellar=spheroid%luminositiesStellar()
             call luminositiesStellar%reset()
             call spheroid%luminositiesStellarSet(luminositiesStellar)
          else
             specificAngularMomentum=spheroid%angularMomentum()/massSpheroid
          end if
          ! Reset the gas, abundances and angular momentum of the spheroid.
          call spheroid%        massGasSet(                                         0.0d0)
          call spheroid%  abundancesGasSet(                                zeroAbundances)
          call spheroid%angularMomentumSet(specificAngularMomentum*spheroid%massStellar())
          ! Indicate that ODE evolution should continue after this state change.
          if (status == GSL_Success) status=GSL_Continue
       end if
       ! Trap negative stellar masses.
       if (spheroid%massStellar() < 0.0d0) then
          ! Check if this exceeds the maximum previously recorded error.
          fractionalError=   abs(spheroid%massStellar()) &
               &          /(                             &
               &             abs(spheroid%massGas    ()) &
               &            +abs(spheroid%massStellar()) &
               &           )
          !$omp critical (Standard_Spheroid_Post_Evolve_Check)
          if (fractionalError > fractionalErrorMaximum) then
             ! Report a warning.
             message='Warning: spheroid has negative stellar mass (fractional error exceeds any previously reported):'//char(10)
             message=message//'  Node index            = '//node%index() //char(10)
             write (valueString,'(e12.6)') spheroid%massGas    ()
             message=message//'  Spheroid gas mass     = '//trim(valueString)//char(10)
             write (valueString,'(e12.6)') spheroid%massStellar()
             message=message//'  Spheroid stellar mass = '//trim(valueString)//char(10)
             write (valueString,'(e12.6)') fractionalError
             message=message//'  Error measure         = '//trim(valueString)//char(10)
             if (fractionalErrorMaximum == 0.0d0) then
                ! This is the first time this warning has been issued, so give some extra information.
                message=message//'  Stellar mass will be reset to zero (in future cases also).'//char(10)
                message=message//'  Future cases will be reported only when they exceed the previous maximum error measure.'//char(10)
                message=message//'  Negative masses are due to numerically inaccuracy in the ODE solutions.'//char(10)
                message=message//'  If significant, consider using a higher tolerance in the ODE solver.'
             end if
             call displayMessage(message,verbosityLevelWarn)
             ! Store the new maximum fractional error.
             fractionalErrorMaximum=fractionalError
          end if
          !$omp end critical (Standard_Spheroid_Post_Evolve_Check)
          ! Get the specific angular momentum of the spheroid material
          massSpheroid= spheroid%massGas    () &
               &       +spheroid%massStellar()
          if (massSpheroid == 0.0d0) then
             specificAngularMomentum=0.0d0
             call spheroid%      massGasSet(         0.0d0)
             call spheroid%abundancesGasSet(zeroAbundances)
          else
             specificAngularMomentum=spheroid%angularMomentum()/massSpheroid
          end if
          ! Reset the stellar mass, abundances and angular momentum of the spheroid.
          call spheroid%        massStellarSet(                                 0.0d0)
          call spheroid%  abundancesStellarSet(                        zeroAbundances)
          call spheroid%angularMomentumSet(specificAngularMomentum*spheroid%massGas())
          historyStellar=spheroid%starFormationHistory    ()
          call historyStellar%reset()
          call spheroid      %starFormationHistorySet   (historyStellar)
          historyStellar=spheroid%stellarPropertiesHistory()
          call historyStellar%reset()
          call spheroid      %stellarPropertiesHistorySet(historyStellar)
          ! Indicate that ODE evolution should continue after this state change.
          if (status == GSL_Success) status=GSL_Continue
       end if
       ! Trap negative angular momentum.
       if (spheroid%angularMomentum() < 0.0d0) then
          ! Estimate a reasonable specific angular momentum.
          specificAngularMomentum=+spheroid%radius                         () &
               &                  *spheroid%velocity                       () &
               &                  /         ratioAngularMomentumScaleRadius
          ! Get the mass of the spheroid.
          massSpheroid           =+spheroid%massGas                        () &
               &                  +spheroid%massStellar                    ()
          ! Reset the angular momentum of the spheroid.
          call spheroid%angularMomentumSet(specificAngularMomentum*massSpheroid)
          ! Indicate that ODE evolution should continue after this state change.
          if (status == GSL_Success) status=GSL_Continue
       end if
    end select
    return
  end subroutine Node_Component_Spheroid_Standard_Post_Step

  subroutine Node_Component_Spheroid_Standard_Mass_Gas_Sink_Rate(self,rate,interrupt,interruptProcedure)
    !!{
    Account for a sink of gaseous material in the standard spheroid.
    !!}
    use :: Abundances_Structure, only : operator(*)
    use :: Error               , only : Error_Report
    use :: Galacticus_Nodes    , only : interruptTask, nodeComponentSpheroid
    implicit none
    class           (nodeComponentSpheroid), intent(inout)                    :: self
    logical                                , intent(inout), optional          :: interrupt
    procedure       (interruptTask        ), intent(inout), optional, pointer :: interruptProcedure
    double precision                       , intent(in   )                    :: rate
    double precision                                                          :: gasMass           , stellarMass
    !$GLC attributes unused :: interrupt, interruptProcedure

    ! Trap cases where an attempt is made to add gas via this sink function.
    if (rate > 0.0d0) call Error_Report(                                                      &
         &                              'attempt to add mass via sink in standard spheroid'// &
         &                              {introspection:location}                              &
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
  end subroutine Node_Component_Spheroid_Standard_Mass_Gas_Sink_Rate

  subroutine Node_Component_Spheroid_Standard_Energy_Gas_Input_Rate(self,rate,interrupt,interruptProcedure)
    !!{
    Handles input of energy into the spheroid gas from other components (e.g. black holes). The energy input rate should be in
    units of $M_\odot$ km$^2$ s$^{-2}$ Gyr$^{-1}$.
    !!}
    use :: Abundances_Structure, only : abundances   , operator(*)
    use :: Error               , only : Error_Report
    use :: Galacticus_Nodes    , only : interruptTask, nodeComponentHotHalo, nodeComponentSpheroid, treeNode
    implicit none
    class           (nodeComponentSpheroid        ), intent(inout)                    :: self
    logical                                        , intent(inout), optional          :: interrupt
    procedure       (interruptTask                ), intent(inout), optional, pointer :: interruptProcedure
    double precision                               , intent(in   )                    :: rate
    class           (nodeComponentHotHalo         )                         , pointer :: hotHalo
    type            (treeNode                     )                         , pointer :: node
    type            (abundances                   )                                   :: abundancesOutflowRate
    double precision                                                                  :: angularMomentumOutflowRate, gasMass         , &
         &                                                                               massOutflowRate           , spheroidVelocity, &
         &                                                                               stellarMass
    !$GLC attributes unused :: interrupt, interruptProcedure

    ! Trap cases where an attempt is made to remove energy via this input function.
    if (rate < 0.0d0) call Error_Report(                                                                 &
         &                              'attempt to remove energy via input pipe in standard spheroid'// &
         &                              {introspection:location}                                         &
         &                             )

    ! Return if no adjustment is being made.
    if (rate == 0.0d0) return
    ! Get the gas mass present.
    gasMass    =self%massGas    ()
    ! Get the stellar mass present.
    stellarMass=self%massStellar()
    ! If gas is present, adjust the rates.
    if (gasMass > 0.0d0 .and. gasMass+stellarMass > 0.0d0) then
       ! Compute outflow rates of quantities and adjust rates in the spheroid appropriately.
       spheroidVelocity=self%velocity()
       if (spheroidVelocity > 0.0d0) then
          massOutflowRate           =efficiencyEnergeticOutflow*rate/spheroidVelocity**2
          angularMomentumOutflowRate=(massOutflowRate/(gasMass+stellarMass))*self%angularMomentum()
          abundancesOutflowRate     =(massOutflowRate/ gasMass             )*self%abundancesGas  ()
          call self%        massGasRate(-           massOutflowRate)
          call self%angularMomentumRate(-angularMomentumOutflowRate)
          call self%  abundancesGasRate(-     abundancesOutflowRate)
          ! Add outflowing rates to the hot halo component.
          node    => self%host   ()
          hotHalo => node%hotHalo()
          call hotHalo%           outflowingMassRate(           massOutflowRate)
          call hotHalo%outflowingAngularMomentumRate(angularMomentumOutflowRate)
          call hotHalo%     outflowingAbundancesRate(     abundancesOutflowRate)
       end if
    end if
    return
  end subroutine Node_Component_Spheroid_Standard_Energy_Gas_Input_Rate

  subroutine Node_Component_Spheroid_Standard_Star_Formation_History_Rate(self,rate,interrupt,interruptProcedure)
    !!{
    Adjust the rates for the star formation history.
    !!}
    use :: Galacticus_Nodes , only : interruptTask, nodeComponentSpheroid, nodeComponentSpheroidStandard
    implicit none
    class    (nodeComponentSpheroid), intent(inout)                    :: self
    type     (history              ), intent(in   )                    :: rate
    logical                         , intent(inout), optional          :: interrupt
    procedure(interruptTask        ), intent(inout), optional, pointer :: interruptProcedure
    type     (history              )                                   :: starFormationHistory

    ! Get the star formation history in the spheroid.
    starFormationHistory=self%starFormationHistory()
    ! Check if the range of this history is sufficient.
    if (starFormationHistory_%rangeIsSufficient(starFormationHistory,rate)) then
       ! Range is sufficient, call the intrinsic rate function.
       select type (self)
       class is (nodeComponentSpheroidStandard)
          call self%starFormationHistoryRateIntrinsic(rate)
       end select
    else
       ! Range is insufficient.
       if (allocated(starFormationHistoryTemplate)) deallocate(starFormationHistoryTemplate)
       allocate(starFormationHistoryTemplate(size(rate%time)))
       starFormationHistoryTemplate =  rate%time
       interrupt                    =  .true.
       interruptProcedure           => Node_Component_Spheroid_Standard_Star_Formation_History_Extend
    end if
    return
  end subroutine Node_Component_Spheroid_Standard_Star_Formation_History_Rate

  subroutine Node_Component_Spheroid_Standard_Stellar_Prprts_History_Rate(self,rate,interrupt,interruptProcedure)
    !!{
    Adjust the rates for the stellar properties history.
    !!}
    use :: Error            , only : Error_Report
    use :: Galacticus_Nodes , only : interruptTask, nodeComponentSpheroid, nodeComponentSpheroidStandard
    implicit none
    class    (nodeComponentSpheroid), intent(inout)                    :: self
    type     (history              ), intent(in   )                    :: rate
    logical                         , intent(inout), optional          :: interrupt
    procedure(interruptTask        ), intent(inout), optional, pointer :: interruptProcedure
    type     (history              )                                   :: stellarPropertiesHistory

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
       interruptProcedure => Node_Component_Spheroid_Standard_Stellar_Prprts_History_Extend
       return
    end if
    ! Adjust the rate.
    select type (self)
    class is (nodeComponentSpheroidStandard)
       call self%stellarPropertiesHistoryRateIntrinsic(rate)
    end select
    return
  end subroutine Node_Component_Spheroid_Standard_Stellar_Prprts_History_Rate

  !![
  <scaleSetTask>
   <unitName>Node_Component_Spheroid_Standard_Scale_Set</unitName>
  </scaleSetTask>
  !!]
  subroutine Node_Component_Spheroid_Standard_Scale_Set(node)
    !!{
    Set scales for properties of {\normalfont \ttfamily node}. Note that gas masses get an additional scaling down since they can approach
    zero and we'd like to prevent them from becoming negative.
    !!}
    use :: Abundances_Structure          , only : abs                     , max              , operator(*)          , unitAbundances               , &
         &                                        abundances
    use :: Galacticus_Nodes              , only : defaultSpheroidComponent, nodeComponentDisk, nodeComponentSpheroid, nodeComponentSpheroidStandard, &
         &                                        treeNode
    use :: Histories                     , only : history                 , operator(*)
    use :: Stellar_Luminosities_Structure, only : abs                     , max              , operator(*)          , stellarLuminosities          , &
         &                                        unitStellarLuminosities
    implicit none
    type            (treeNode             ), intent(inout), pointer :: node
    class           (nodeComponentSpheroid)               , pointer :: spheroid
    class           (nodeComponentDisk    )               , pointer :: disk
    double precision                       , parameter              :: massMinimum                   =1.0d+0
    double precision                       , parameter              :: angularMomentumMinimum        =0.1d+0
    double precision                       , parameter              :: luminosityMinimum             =1.0d+0
    double precision                                                :: angularMomentum                      , mass
    type            (history              )                         :: stellarPopulationHistoryScales
    type            (stellarLuminosities  )                         :: stellarLuminositiesScale
    type            (abundances           )                         :: abundancesScale

    ! Check if we are the default method.
    if (.not.defaultSpheroidComponent%standardIsActive()) return
    ! Get the spheroid component.
    spheroid => node%spheroid()
    ! Check if a standard spheroid component exists.
    select type (spheroid)
    class is (nodeComponentSpheroidStandard)
       ! Get the disk component.
       disk => node%disk()
       ! Set scale for angular momentum.
       !! The scale here (and for other quantities below) combines the mass of disk and spheroid. This avoids attempts to solve
       !! tiny spheroids to high precision in massive disk galaxies. This is particularly important here as dynamical processes
       !! (such as bar instabilities) can transfer mass/angular momentum/metals from such a disk to the spheroid. In cases where
       !! such a process creates the spheroid then, by definition, the spheroid initially has no mass/angular momentum.
       angularMomentum=+abs(disk    %angularMomentum()) &
            &          +abs(spheroid%angularMomentum())
       call spheroid%angularMomentumScale  (max(angularMomentum,angularMomentumMinimum))
       ! Set scale for masses.
       mass           =max(                                                      &
            &              +abs(disk%massGas    ())+abs(spheroid%massGas    ())  &
            &              +abs(disk%massStellar())+abs(spheroid%massStellar()), &
            &              +massMinimum                                          &
            &             )
       call spheroid%          massGasScale(mass)
       call spheroid%      massStellarScale(mass)
       call spheroid%massStellarFormedScale(mass)
       ! Set scales for abundances if necessary.
       if (abundancesCount > 0) then
          ! Set scale for abundances.
          abundancesScale=+max(                                     &
               &               +abs(+disk    %abundancesGas    ())  &
               &               +abs(+disk    %abundancesStellar())  &
               &               +abs(+spheroid%abundancesGas    ())  &
               &               +abs(+spheroid%abundancesStellar()), &
               &               +max(                                &
               &                    +mass                           &
               &                    *toleranceRelativeMetallicity , &
               &                    +massMinimum                    &
               &                   )                                &
               &                    *unitAbundances                 &
               &              )
          ! Set scale for gas abundances.
          call spheroid%abundancesGasScale    (                     &
               &                               +abundancesScale     &
               &                              )
          ! Set scale for stellar abundances.
          call spheroid%abundancesStellarScale(                     &
               &                               +abundancesScale     &
               &                              )
       end if
       ! Set scales for stellar luminosities.
       stellarLuminositiesScale=max(                                      &
            &                       +abs(disk    %luminositiesStellar())  &
            &                       +abs(spheroid%luminositiesStellar()), &
            &                           +unitStellarLuminosities          &
            &                           *luminosityMinimum                &
            &                      )
       call stellarLuminositiesScale%truncate                (spheroid                %luminositiesStellar())
       call spheroid                %luminositiesStellarScale(stellarLuminositiesScale                      )
       ! Set scales for stellar population properties history.
       stellarPopulationHistoryScales=spheroid%stellarPropertiesHistory()
       call stellarPopulationProperties_%scales   (spheroid%massStellar(),spheroid%abundancesStellar(),stellarPopulationHistoryScales)
       call spheroid%stellarPropertiesHistoryScale(                                                    stellarPopulationHistoryScales)
       call stellarPopulationHistoryScales%destroy()
       stellarPopulationHistoryScales=spheroid%starFormationHistory()
       call starFormationHistory_%scales          (stellarPopulationHistoryScales,node,spheroid%massStellar(),spheroid%massGas(),spheroid%abundancesStellar())
       call spheroid%starFormationHistoryScale    (stellarPopulationHistoryScales                                                         )
       call stellarPopulationHistoryScales%destroy()
    end select
    return
  end subroutine Node_Component_Spheroid_Standard_Scale_Set

  !![
  <inactiveSetTask>
   <unitName>Node_Component_Spheroid_Standard_Inactive</unitName>
  </inactiveSetTask>
  !!]
  subroutine Node_Component_Spheroid_Standard_Inactive(node)
    !!{
    Set Jacobian zero status for properties of {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSpheroid, nodeComponentSpheroidStandard, treeNode
    implicit none
    type (treeNode             ), intent(inout), pointer :: node
    class(nodeComponentSpheroid)               , pointer :: spheroid

    ! Get the spheroid component.
    spheroid => node%spheroid()
    ! Check if an standard spheroid component exists.
    select type (spheroid)
    class is (nodeComponentSpheroidStandard)
       if (inactiveLuminositiesStellar) call spheroid%luminositiesStellarInactive()
    end select
    return
  end subroutine Node_Component_Spheroid_Standard_Inactive

  subroutine satelliteMerger(self,node)
    !!{
    Transfer any standard spheroid associated with {\normalfont \ttfamily node} to its host halo.
    !!}
    use :: Abundances_Structure            , only : zeroAbundances
    use :: Error                           , only : Error_Report
    use :: Galacticus_Nodes                , only : nodeComponentDisk      , nodeComponentSpheroid    , nodeComponentSpheroidStandard, treeNode
    use :: Satellite_Merging_Mass_Movements, only : destinationMergerDisk  , destinationMergerSpheroid, destinationMergerUnmoved     , enumerationDestinationMergerType
    use :: Satellite_Merging_Remnant_Sizes , only : remnantNoChange
    use :: Stellar_Luminosities_Structure  , only : zeroStellarLuminosities
    implicit none
    class           (*                               ), intent(inout) :: self
    type            (treeNode                        ), intent(inout) :: node
    type            (treeNode                        ), pointer       :: nodeHost
    class           (nodeComponentDisk               ), pointer       :: diskHost                      , disk
    class           (nodeComponentSpheroid           ), pointer       :: spheroidHost                  , spheroid
    type            (history                         )                :: historyDisk                   , historySpheroid                , &
         &                                                               history_
    double precision                                                  :: angularMomentum               , diskSpecificAngularMomentum    , &
         &                                                               massSpheroid                  , spheroidSpecificAngularMomentum, &
         &                                                               radiusRemnant                 , velocityCircularRemnant        , &
         &                                                               angularMomentumSpecificRemnant
    type            (enumerationDestinationMergerType)                :: destinationGasSatellite       , destinationGasHost             , &
         &                                                               destinationStarsHost          , destinationStarsSatellite
    logical                                                           :: mergerIsMajor
    !$GLC attributes unused :: self

    ! Get the spheroid component, creating it if need be.
    spheroid => node%spheroid(autoCreate=.true.)
    select type (spheroid)
    class is (nodeComponentSpheroidStandard)
       disk => node%disk()
       ! Find the node to merge with.
       nodeHost     => node    %mergesWith(                 )
       diskHost     => nodeHost%disk      (autoCreate=.true.)
       spheroidHost => nodeHost%spheroid  (autoCreate=.true.)
       ! Get specific angular momentum of the host spheroid and disk material.
       if (spheroidHost%massGas()+spheroidHost%massStellar() > 0.0d0) then
          spheroidSpecificAngularMomentum=   spheroidHost%angularMomentum() &
               &                          /(                                &
               &                             spheroidHost%massGas        () &
               &                            +spheroidHost%massStellar    () &
               &                           )
       else
          spheroidSpecificAngularMomentum=0.0d0
       end if
       if (diskHost%massGas()+diskHost%massStellar() > 0.0d0) then
          diskSpecificAngularMomentum=   diskHost%angularMomentum() &
               &                      /(                            &
               &                         diskHost%massGas        () &
               &                        +diskHost%massStellar    () &
               &                       )
       else
          diskSpecificAngularMomentum=0.0d0
       end if
       ! Get mass movement descriptors.
       call mergerMassMovements_%get(node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
       ! Move gas material within the host if necessary.
       select case (destinationGasHost%ID)
       case (destinationMergerDisk   %ID)
          call diskHost    %        massGasSet(                                 &
               &                                diskHost    %massGas        ()  &
               &                               +spheroidHost%massGas        ()  &
               &                              )
          call diskHost    %  abundancesGasSet(                                 &
               &                                diskHost    %abundancesGas  ()  &
               &                               +spheroidHost%abundancesGas  ()  &
               &                              )
          call diskHost    %angularMomentumSet(                                 &
               &                                diskHost    %angularMomentum()  &
               &                               +spheroidHost%massGas        ()  &
               &                               *spheroidSpecificAngularMomentum &
               &                              )
          call spheroidHost%angularMomentumSet(                                 &
               &                                spheroidHost%angularMomentum()  &
               &                               -spheroidHost%massGas        ()  &
               &                               *spheroidSpecificAngularMomentum &
               &                              )
          call spheroidHost%        massGasSet(                                 &
               &                                0.0d0                           &
               &                              )
          call spheroidHost%  abundancesGasSet(                                 &
               &                                zeroAbundances                  &
               &                              )
       case (destinationMergerSpheroid%ID)
          call spheroidHost%        massGasSet(                                 &
               &                                spheroidHost%massGas        ()  &
               &                               +diskHost    %massGas        ()  &
               &                              )
          call diskHost    %angularMomentumSet(                                 &
               &                                diskHost%angularMomentum    ()  &
               &                               -diskHost%massGas            ()  &
               &                               *diskSpecificAngularMomentum     &
               &                              )
          call spheroidHost%  abundancesGasSet(                                 &
               &                                spheroidHost%abundancesGas  ()  &
               &                               +diskHost    %abundancesGas  ()  &
               &                              )
          call diskHost    %        massGasSet(                                 &
               &                                0.0d0                           &
               &                              )
          call diskHost    %  abundancesGasSet(                                 &
               &                                zeroAbundances                  &
               &                              )
       case (destinationMergerUnmoved%ID)
          ! Do nothing.
       case default
          call Error_Report('unrecognized movesTo descriptor'//{introspection:location})
       end select

       ! Move stellar material within the host if necessary.
       select case (destinationStarsHost%ID)
       case (destinationMergerDisk   %ID)
          call     diskHost%        massStellarSet(                                    &
               &                                    diskHost    %        massStellar() &
               &                                   +spheroidHost%        massStellar() &
               &                                  )
          call     diskHost%  abundancesStellarSet(                                    &
               &                                    diskHost    %  abundancesStellar() &
               &                                   +spheroidHost%  abundancesStellar() &
               &                                  )
          call     diskHost%luminositiesStellarSet(                                    &
               &                                    diskHost    %luminositiesStellar() &
               &                                   +spheroidHost%luminositiesStellar() &
               &                                  )
          call     diskHost%    angularMomentumSet(                                    &
               &                                    diskHost    %    angularMomentum() &
               &                                   +spheroidHost%        massStellar() &
               &                                   *spheroidSpecificAngularMomentum    &
               &                                  )
          call spheroidHost%    angularMomentumSet(                                    &
               &                                    spheroidHost%    angularMomentum() &
               &                                   -spheroidHost%        massStellar() &
               &                                   *spheroidSpecificAngularMomentum    &
               &                                  )
          call spheroidHost%        massStellarSet(                                    &
               &                                    0.0d0                              &
               &                                  )
          call spheroidHost%  abundancesStellarSet(                                    &
               &                                    zeroAbundances                     &
               &                                  )
          call spheroidHost%luminositiesStellarSet(                                    &
               &                                    zeroStellarLuminosities            &
               &                                  )
          ! Also add stellar properties histories.
          historyDisk    =    diskHost%stellarPropertiesHistory()
          historySpheroid=spheroidHost%stellarPropertiesHistory()
          call historyDisk    %interpolatedIncrement(historySpheroid    )
          call historySpheroid%reset    (                   )
          call diskHost    %stellarPropertiesHistorySet(historyDisk    )
          call spheroidHost%stellarPropertiesHistorySet(historySpheroid)
          ! Also add star formation histories.
          historyDisk    =    diskHost%starFormationHistory()
          historySpheroid=spheroidHost%starFormationHistory()
          call starFormationHistory_%move                   (nodeHost,nodeHost,historyDisk,historySpheroid)
          call diskHost             %starFormationHistorySet(                  historyDisk                )
          call spheroidHost         %starFormationHistorySet(                              historySpheroid)
          call historyDisk          %destroy                (                                             )
          call historySpheroid      %destroy                (                                             )
       case (destinationMergerSpheroid%ID)
          call spheroidHost%        massStellarSet(                                    &
               &                                    spheroidHost%        massStellar() &
               &                                   +diskHost    %        massStellar() &
               &                                  )
          call     diskHost%    angularMomentumSet( diskHost    %    angularMomentum() &
               &                                   -diskHost    %        massStellar() &
               &                                   *diskSpecificAngularMomentum        &
               &                                  )
          call spheroidHost%  abundancesStellarSet(                                    &
               &                                    spheroidHost%  abundancesStellar() &
               &                                   +diskHost    %  abundancesStellar() &
               &                                  )
          call spheroidHost%luminositiesStellarSet( spheroidHost%luminositiesStellar() &
               &                                   +diskHost    %luminositiesStellar() &
               &                                  )
          call     diskHost%        massStellarSet(                                    &
               &                                    0.0d0                              &
               &                                  )
          call     diskHost%  abundancesStellarSet(                                    &
               &                                    zeroAbundances                     &
               &                                  )
          call     diskHost%luminositiesStellarSet(                                    &
               &                                    zeroStellarLuminosities            &
               &                                  )
          ! Also add stellar properties histories.
          historyDisk    =    diskHost%stellarPropertiesHistory()
          historySpheroid=spheroidHost%stellarPropertiesHistory()
          call historySpheroid%interpolatedIncrement(historyDisk)
          call historyDisk    %reset    (           )
          call spheroidHost   %stellarPropertiesHistorySet(historySpheroid)
          call diskHost       %stellarPropertiesHistorySet(historyDisk    )
          ! Also add star formation histories.
          historyDisk    =diskHost    %starFormationHistory()
          historySpheroid=spheroidHost%starFormationHistory()
          call starFormationHistory_%move                   (nodeHost,nodeHost,historySpheroid,historyDisk)
          call spheroidHost         %starFormationHistorySet(                  historySpheroid            )
          call diskHost             %starFormationHistorySet(                                  historyDisk)
          call historyDisk          %destroy                (                                             )
          call historySpheroid      %destroy                (                                             )
          historyDisk    =diskHost    %starFormationHistory()
          historySpheroid=spheroidHost%starFormationHistory()
       case (destinationMergerUnmoved%ID)
          ! Do nothing.
       case default
          call Error_Report('unrecognized movesTo descriptor'//{introspection:location})
       end select
       ! If the entire host disk/spheroid (gas plus stars) was moved to the spheroid/disk, ensure that the
       ! corresponding angular momentum is precisely zero.
       if (diskHost    %massStellar()+diskHost    %massGas() == 0.0d0 .and. diskHost%angularMomentum() /= 0.0d0) &
            & call diskHost%angularMomentumSet    (0.0d0)
       if (spheroidHost%massStellar()+spheroidHost%massGas() == 0.0d0                                          ) &
            & call spheroidHost%angularMomentumSet(0.0d0)

       ! Get specific angular momentum of the spheroid material.
       massSpheroid=spheroid%massGas()+spheroid%massStellar()
       if (massSpheroid > 0.0d0) then
          spheroidSpecificAngularMomentum=spheroid%angularMomentum()/massSpheroid

          ! Move the gas component of the standard spheroid to the host.
          select case (destinationGasSatellite%ID)
          case (destinationMergerDisk   %ID)
             call     diskHost%        massGasSet(                                 &
                  &                                diskHost    %        massGas()  &
                  &                               +spheroid%            massGas()  &
                  &                              )
             call     diskHost%  abundancesGasSet(                                 &
                  &                                diskHost    %  abundancesGas()  &
                  &                               +spheroid    %  abundancesGas()  &
                  &                              )
             call     diskHost%angularMomentumSet( diskHost    %angularMomentum()  &
                  &                               +spheroid    %        massGas()  &
                  &                               *spheroidSpecificAngularMomentum &
                  &                              )
          case (destinationMergerSpheroid%ID)
             call spheroidHost%        massGasSet(                                 &
                  &                                spheroidHost%        massGas()  &
                  &                               +spheroid    %        massGas()  &
                  &                              )
             call spheroidHost%  abundancesGasSet(                                 &
                  &                                spheroidHost%  abundancesGas()  &
                  &                               +spheroid    %  abundancesGas()  &
                  &                              )
          case default
             call Error_Report('unrecognized movesTo descriptor'//{introspection:location})
          end select
          call spheroid%      massGasSet(0.0d0         )
          call spheroid%abundancesGasSet(zeroAbundances)
          ! Move the stellar component of the standard spheroid to the host.
          select case (destinationStarsSatellite%ID)
          case (destinationMergerDisk   %ID)
             call diskHost    %        massStellarSet(                                 &
                  &                                    diskHost%        massStellar()  &
                  &                                   +spheroid%        massStellar()  &
                  &                                  )
             call diskHost    %  abundancesStellarSet( diskHost%  abundancesStellar()  &
                  &                                   +spheroid%  abundancesStellar()  &
                  &                                  )
             call diskHost    %luminositiesStellarSet( diskHost%luminositiesStellar()  &
                  &                                   +spheroid%luminositiesStellar()  &
                  &                                  )
             call diskHost    %    angularMomentumSet( diskHost%    angularMomentum()  &
                  &                                   +spheroid%        massStellar()  &
                  &                                   *spheroidSpecificAngularMomentum &
                  &                                  )
             ! Also add stellar properties histories.
             historySpheroid=spheroid%stellarPropertiesHistory()
             history_       =diskHost%stellarPropertiesHistory()
             call history_       %interpolatedIncrement                  (historySpheroid)
             call historySpheroid%reset                      (               )
             call diskHost       %stellarPropertiesHistorySet(history_       )
             call spheroid       %stellarPropertiesHistorySet(historySpheroid)
             ! Also add star formation histories.
             historySpheroid=spheroid%starFormationHistory()
             history_       =diskHost%starFormationHistory()
             call starFormationHistory_%move                   (nodeHost,node,history_,historySpheroid)
             call diskHost             %starFormationHistorySet(              history_                )
             call spheroid             %starFormationHistorySet(                       historySpheroid)
             call history_             %destroy                (                                      )
             call historySpheroid      %destroy                (                                      )
          case (destinationMergerSpheroid%ID)
             call spheroidHost%        massStellarSet( spheroidHost%        massStellar() &
                  &                                   +spheroid    %        massStellar() &
                  &                                  )
             call spheroidHost%  abundancesStellarSet( spheroidHost%  abundancesStellar() &
                  &                                   +spheroid    %  abundancesStellar() &
                  &                                  )
             call spheroidHost%luminositiesStellarSet( spheroidHost%luminositiesStellar() &
                  &                                   +spheroid    %luminositiesStellar() &
                  &                                  )
             ! Also add stellar properties histories.
             historySpheroid=spheroid    %stellarPropertiesHistory()
             history_       =spheroidHost%stellarPropertiesHistory()
             call history_        %interpolatedIncrement(historySpheroid)
             call historySpheroid%reset                      (               )
             call spheroidHost   %stellarPropertiesHistorySet(history_       )
             call spheroid       %stellarPropertiesHistorySet(historySpheroid)
             ! Also add star formation histories.
             historySpheroid=spheroid    %starFormationHistory()
             history_       =spheroidHost%starFormationHistory()
             call starFormationHistory_%move                   (nodeHost,node,history_,historySpheroid)             
             call spheroidHost         %starFormationHistorySet(              history_                )
             call spheroid             %starFormationHistorySet(                       historySpheroid)
             call history_             %destroy                (                                      )
             call historySpheroid      %destroy                (                                      )
          case default
             call Error_Report('unrecognized movesTo descriptor'//{introspection:location})
          end select
          call spheroid%        massStellarSet(0.0d0                  )
          call spheroid%  abundancesStellarSet(zeroAbundances         )
          call spheroid%luminositiesStellarSet(zeroStellarLuminosities)
          call spheroid%    angularMomentumSet(0.0d0                  )
       end if
       ! Set the angular momentum of the spheroid.
       call mergerRemnantSize_%get(node,radiusRemnant,velocityCircularRemnant,angularMomentumSpecificRemnant)
       if (angularMomentumSpecificRemnant /= remnantNoChange) then
          ! Note that the remnant specific angular momentum computed by the merger remnant modules automatically gives the mean
          ! specific angular momentum of the component by virtue of the fact that it computes the ratio of the actual angular
          ! momentum to the contribution from the component's own rotation curve at its scale radius.
          angularMomentum=angularMomentumSpecificRemnant*(spheroidHost%massGas()+spheroidHost%massStellar())
          call spheroidHost%angularMomentumSet(angularMomentum)
       end if
    end select
    return
  end subroutine satelliteMerger

  !![
  <radiusSolverPlausibility>
   <unitName>Node_Component_Spheroid_Standard_Radius_Solver_Plausibility</unitName>
   <after>Node_Component_Basic_Standard_Plausibility</after>
  </radiusSolverPlausibility>
  !!]
  subroutine Node_Component_Spheroid_Standard_Radius_Solver_Plausibility(node)
    !!{
    Determines whether the spheroid is physically plausible for radius solving tasks. Require that it have non-zero mass and angular momentum.
    !!}
    use :: Galacticus_Nodes, only : defaultSpheroidComponent, nodeComponentSpheroid, nodeComponentSpheroidStandard, treeNode
    implicit none
    type            (treeNode                ), intent(inout) :: node
    class           (nodeComponentSpheroid   ), pointer       :: spheroid
    double precision                          , parameter     :: angularMomentumMaximum=1.0d+1
    double precision                          , parameter     :: angularMomentumMinimum=1.0d-6
    double precision                                          :: angularMomentumScale

    ! Return immediately if our method is not selected.
    if     (                                                &
         &  .not.                                           &
         &   (                                              &
         &     defaultSpheroidComponent%standardIsActive () &
         &    .and.                                         &
         &     node                %isPhysicallyPlausible   &
         &    .and.                                         &
         &     node                %isSolvable              &
         &   )                                              &
         & ) return

    ! Determine the plausibility of the current spheroid.
    spheroid => node%spheroid()
    select type (spheroid)
    class is (nodeComponentSpheroidStandard)
       if (spheroid%massStellar()+spheroid%massGas() >= 0.0d0 .and. spheroid%angularMomentum() > 0.0d0) then
          angularMomentumScale=(                                           &
               &                 spheroid%massStellar()                    &
               &                +spheroid%massGas    ()                    &
               &               )                                           &
               &               * darkMatterHaloScale_%radiusVirial  (node) &
               &               * darkMatterHaloScale_%velocityVirial(node)
          if     (                                                                          &
               &   spheroid%angularMomentum() > angularMomentumMaximum*angularMomentumScale &
               &  .or.                                                                      &
               &   spheroid%angularMomentum() < angularMomentumMinimum*angularMomentumScale &
               & ) then
             ! Ignore spheroids with angular momenta greatly exceeding that which would be expected if they had a radius
             ! comparable to the virial radius of their halo.
             node%isPhysicallyPlausible=.false.
          end if
       end if
    end select
    return
  end subroutine Node_Component_Spheroid_Standard_Radius_Solver_Plausibility

  double precision function Node_Component_Spheroid_Standard_Radius_Solve(node)
    !!{
    Return the circular radius of the standard spheroid.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSpheroid, treeNode
    implicit none
    type (treeNode             ), intent(inout) :: node
    class(nodeComponentSpheroid), pointer       :: spheroid

    spheroid => node%spheroid()
    Node_Component_Spheroid_Standard_Radius_Solve=spheroid%radius()
    return
  end function Node_Component_Spheroid_Standard_Radius_Solve

  double precision function Node_Component_Spheroid_Standard_Velocity_Solve(node)
    !!{
    Return the circular velocity of the standard spheroid.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSpheroid, treeNode
    implicit none
    type (treeNode             ), intent(inout) :: node
    class(nodeComponentSpheroid), pointer       :: spheroid

    spheroid => node%spheroid()
    Node_Component_Spheroid_Standard_Velocity_Solve=spheroid%velocity()
    return
  end function Node_Component_Spheroid_Standard_Velocity_Solve

  subroutine Node_Component_Spheroid_Standard_Radius_Solve_Set(node,radius)
    !!{
    Set the scale radius of the standard spheroid.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSpheroid, treeNode
    implicit none
    type            (treeNode             ), intent(inout) :: node
    double precision                       , intent(in   ) :: radius
    class           (nodeComponentSpheroid), pointer       :: spheroid

    spheroid => node%spheroid()
    call spheroid%radiusSet(max(radius,0.0d0))
    return
  end subroutine Node_Component_Spheroid_Standard_Radius_Solve_Set

  subroutine Node_Component_Spheroid_Standard_Velocity_Solve_Set(node,velocity)
    !!{
    Set the scale velocity of the standard spheroid.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSpheroid, treeNode
    implicit none
    type            (treeNode             ), intent(inout) :: node
    double precision                       , intent(in   ) :: velocity
    class           (nodeComponentSpheroid), pointer       :: spheroid

    spheroid => node%spheroid()
    call spheroid%velocitySet(max(velocity,0.0d0))
    return
  end subroutine Node_Component_Spheroid_Standard_Velocity_Solve_Set

  !![
  <radiusSolverTask>
   <unitName>Node_Component_Spheroid_Standard_Radius_Solver</unitName>
  </radiusSolverTask>
  !!]
  subroutine Node_Component_Spheroid_Standard_Radius_Solver(node,componentActive,specificAngularMomentumRequired,specificAngularMomentum,Radius_Get,Radius_Set,Velocity_Get&
       &,Velocity_Set)
    !!{
    Interface for the size solver algorithm.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSpheroid, nodeComponentSpheroidStandard, treeNode
    implicit none
    type            (treeNode                                         ), intent(inout)          :: node
    logical                                                            , intent(  out)          :: componentActive
    logical                                                            , intent(in   )          :: specificAngularMomentumRequired
    double precision                                                   , intent(  out)          :: specificAngularMomentum
    procedure       (Node_Component_Spheroid_Standard_Radius_Solve_Set), intent(  out), pointer :: Radius_Set                     , Velocity_Set
    procedure       (Node_Component_Spheroid_Standard_Radius_Solve    ), intent(  out), pointer :: Radius_Get                     , Velocity_Get
    class           (nodeComponentSpheroid                            )               , pointer :: spheroid
    double precision                                                                            :: angularMomentum                , specificAngularMomentumMean, &
         &                                                                                         massSpheroid

    ! Determine if node has an active disk component supported by this module.
    componentActive=.false.
    spheroid => node%spheroid()
    select type (spheroid)
       class is (nodeComponentSpheroidStandard)
       componentActive=.true.
       ! Get the angular momentum.
       if (specificAngularMomentumRequired) then
          angularMomentum=spheroid%angularMomentum()
          if (angularMomentum >= 0.0d0) then
             ! Compute the specific angular momentum at the scale radius, assuming a flat rotation curve.
             massSpheroid=spheroid%massGas()+spheroid%massStellar()
             if (massSpheroid > 0.0d0) then
                specificAngularMomentumMean=angularMomentum/massSpheroid
             else
                specificAngularMomentumMean=0.0d0
             end if
             specificAngularMomentum=ratioAngularMomentumScaleRadius*specificAngularMomentumMean
          end if
       end if
       ! Associate the pointers with the appropriate property routines.
       Radius_Get   => Node_Component_Spheroid_Standard_Radius_Solve
       Radius_Set   => Node_Component_Spheroid_Standard_Radius_Solve_Set
       Velocity_Get => Node_Component_Spheroid_Standard_Velocity_Solve
       Velocity_Set => Node_Component_Spheroid_Standard_Velocity_Solve_Set
    end select
    return
  end subroutine Node_Component_Spheroid_Standard_Radius_Solver

  subroutine Node_Component_Spheroid_Standard_Initializor(self,timeEnd)
    !!{
    Initializes a standard spheroid component.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDisk, nodeComponentSpheroidStandard, treeNode
    use :: Histories       , only : history
    implicit none
    type            (nodeComponentSpheroidStandard), intent(inout)           :: self
    double precision                               , intent(in   ), optional :: timeEnd
    type            (treeNode                     ), pointer                 :: node
    class           (nodeComponentDisk            ), pointer                 :: disk
    class           (nodeComponentBasic           ), pointer                 :: basic
    type            (history                      )                          :: historyStarFormation        , stellarPropertiesHistory      , &
         &                                                                      diskStarFormationHistory
    logical                                                                  :: createStarFormationHistory  , createStellarPropertiesHistory
    double precision                                                         :: timeBegin

    ! Return if already initialized.
    if (self%isInitialized()) return
    ! Get the associated node.
    node => self%host()
    ! Determine which histories must be created.
    historyStarFormation          =self%starFormationHistory            ()
    createStarFormationHistory    =.not.historyStarFormation    %exists ()
    call                                historyStarFormation    %destroy()
    stellarPropertiesHistory      =self%stellarPropertiesHistory        ()
    createStellarPropertiesHistory=.not.stellarPropertiesHistory%exists ()
    call                                stellarPropertiesHistory%destroy()
    ! Create the stellar properties history.
    if (createStellarPropertiesHistory) then
       call stellarPopulationProperties_%historyCreate(node,stellarPropertiesHistory)
       call self%stellarPropertiesHistorySet(               stellarPropertiesHistory)
    end if
    ! Create the star formation history.
    if (createStarFormationHistory) then
       disk                     => node%disk                ()
       diskStarFormationHistory =  disk%starFormationHistory()
       if (diskStarFormationHistory%exists()) then
          timeBegin = diskStarFormationHistory%time(1)
       else
          basic    => node %basic()
          timeBegin = basic%time ()
       end if
       call starFormationHistory_%create                 (node,historyStarFormation,timeBegin,timeEnd)
       call self                 %starFormationHistorySet(     historyStarFormation                  )
    end if
    ! Record that the spheroid has been initialized.
    call self%isInitializedSet(.true.)
    return
  end subroutine Node_Component_Spheroid_Standard_Initializor

  subroutine Node_Component_Spheroid_Standard_Star_Formation_History_Extend(node,timeEnd)
    !!{
    Extend the range of a star formation history in a standard spheroid component for {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSpheroid, treeNode
    implicit none
    type            (treeNode             ), intent(inout), target   :: node
    double precision                       , intent(in   ), optional :: timeEnd
    class           (nodeComponentSpheroid)               , pointer  :: spheroid
    type            (history              )                          :: historyStarFormation
    !$GLC attributes unused :: timeEnd
    
    ! Get the spheroid component.
    spheroid => node%spheroid()
    ! Extend the range as necessary.
    historyStarFormation=spheroid%starFormationHistory()
    call starFormationHistory_%extend(historyStarFormation,starFormationHistoryTemplate)
    call spheroid%starFormationHistorySet(historyStarFormation)
    return
  end subroutine Node_Component_Spheroid_Standard_Star_Formation_History_Extend

  subroutine Node_Component_Spheroid_Standard_Stellar_Prprts_History_Extend(node,timeEnd)
    !!{
    Extend the range of a stellar properties history in a standard spheroid component for {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSpheroid, treeNode
    implicit none
    type            (treeNode             ), intent(inout), target   :: node
    double precision                       , intent(in   ), optional :: timeEnd
    class           (nodeComponentSpheroid)               , pointer  :: spheroid
    type            (history              )                          :: stellarPropertiesHistory
    !$GLC attributes unused :: timeEnd
    
    ! Get the spheroid component.
    spheroid => node%spheroid()
    ! Extend the range as necessary.
    stellarPropertiesHistory=spheroid%stellarPropertiesHistory()
    call stellarPropertiesHistory%extend(times=stellarPropertiesHistoryTemplate)
    call spheroid%stellarPropertiesHistorySet(stellarPropertiesHistory)
    return
  end subroutine Node_Component_Spheroid_Standard_Stellar_Prprts_History_Extend

  subroutine mergerTreeExtraOutput(self,node,iOutput,treeIndex,nodePassesFilter,treeLock)
    !!{
    Update the star formation history after an output time is reached.
    !!}
    use            :: Galacticus_Nodes          , only : defaultSpheroidComponent, nodeComponentSpheroid, nodeComponentSpheroidStandard, treeNode
    use            :: Galactic_Structure_Options, only : componentTypeSpheroid
    use            :: Histories                 , only : history
    use, intrinsic :: ISO_C_Binding             , only : c_size_t
    use            :: Kind_Numbers              , only : kind_int8
    use            :: Locks                     , only : ompLock
    implicit none
    class  (*                    ), intent(inout)          :: self
    type   (treeNode             ), intent(inout)          :: node
    integer(c_size_t             ), intent(in   )          :: iOutput
    integer(kind=kind_int8       ), intent(in   )          :: treeIndex
    logical                       , intent(in   )          :: nodePassesFilter
    type   (ompLock              ), intent(inout)          :: treeLock
    class  (nodeComponentSpheroid)               , pointer :: spheroid
    type   (history              )                         :: historyStarFormation
    !$GLC attributes unused :: self, treeIndex, nodePassesFilter, treeLock

    ! Check if we are the default method.
    if (.not.defaultSpheroidComponent%standardIsActive()) return
    ! Output the star formation history.
    spheroid             => node    %spheroid            ()
    historyStarFormation =  spheroid%starFormationHistory()
    call starFormationHistory_%update(node,iOutput,historyStarFormation)
    ! Update the star formation history only if a spheroid exists.
    select type (spheroid)
    class is (nodeComponentSpheroidStandard)
       call spheroid%starFormationHistorySet(historyStarFormation)
    end select
    return
  end subroutine mergerTreeExtraOutput

  !![
  <stateStoreTask>
   <unitName>Node_Component_Spheroid_Standard_State_Store</unitName>
  </stateStoreTask>
  !!]
  subroutine Node_Component_Spheroid_Standard_State_Store(stateFile,gslStateFile,stateOperationID)
    !!{
    Write the tabulation state to file.
    !!}
    use            :: Display                              , only : displayMessage          , verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding                        , only : c_ptr                   , c_size_t
    use            :: Node_Component_Spheroid_Standard_Data, only : massDistributionStellar_, massDistributionGas_, kinematicDistribution_
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Storing state for: componentSpheroid -> standard',verbosity=verbosityLevelInfo)
    !![
    <stateStore variables="massDistributionStellar_ massDistributionGas_ kinematicDistribution_ stellarPopulationProperties_ darkMatterHaloScale_ starFormationHistory_ mergerMassMovements_ mergerRemnantSize_"/>
    !!]
    write (stateFile) ratioAngularMomentumScaleRadius
    return
  end subroutine Node_Component_Spheroid_Standard_State_Store

  !![
  <stateRetrieveTask>
   <unitName>Node_Component_Spheroid_Standard_State_Retrieve</unitName>
  </stateRetrieveTask>
  !!]
  subroutine Node_Component_Spheroid_Standard_State_Retrieve(stateFile,gslStateFile,stateOperationID)
    !!{
    Retrieve the tabulation state from the file.
    !!}
    use            :: Display                              , only : displayMessage          , verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding                        , only : c_ptr                   , c_size_t
    use            :: Node_Component_Spheroid_Standard_Data, only : massDistributionStellar_, massDistributionGas_, kinematicDistribution_
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Retrieving state for: componentSpheroid -> standard',verbosity=verbosityLevelInfo)
    !![
    <stateRestore variables="massDistributionStellar_ massDistributionGas_ kinematicDistribution_ stellarPopulationProperties_ darkMatterHaloScale_ starFormationHistory_ mergerMassMovements_ mergerRemnantSize_"/>
    !!]
    read (stateFile) ratioAngularMomentumScaleRadius
    return
  end subroutine Node_Component_Spheroid_Standard_State_Retrieve

end module Node_Component_Spheroid_Standard
