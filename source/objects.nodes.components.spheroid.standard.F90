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

!% Contains a module which implements the standard spheroid component.

module Node_Component_Spheroid_Standard
  !% Implements the standard spheroid component.
  use :: Dark_Matter_Halo_Scales                        , only : darkMatterHaloScaleClass
  use :: Histories                                      , only : history
  use :: Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroids, only : ramPressureStrippingSpheroidsClass
  use :: Satellite_Merging_Mass_Movements               , only : mergerMassMovementsClass
  use :: Satellite_Merging_Remnant_Sizes                , only : mergerRemnantSizeClass
  use :: Satellites_Tidal_Fields                        , only : satelliteTidalFieldClass
  use :: Star_Formation_Histories                       , only : starFormationHistory                        , starFormationHistoryClass
  use :: Stellar_Population_Properties                  , only : stellarPopulationPropertiesClass
  use :: Tidal_Stripping_Mass_Loss_Rate_Spheroids       , only : tidalStrippingSpheroidsClass
  implicit none
  private
  public :: Node_Component_Spheroid_Standard_Rate_Compute       , Node_Component_Spheroid_Standard_Scale_Set                    , &
       &    Node_Component_Spheroid_Standard_Radius_Solver      , Node_Component_Spheroid_Standard_Star_Formation_History_Output, &
       &    Node_Component_Spheroid_Standard_Pre_Evolve         , Node_Component_Spheroid_Standard_Radius_Solver_Plausibility   , &
       &    Node_Component_Spheroid_Standard_Thread_Uninitialize, Node_Component_Spheroid_Standard_Thread_Initialize            , &
       &    Node_Component_Spheroid_Standard_State_Store        , Node_Component_Spheroid_Standard_State_Retrieve               , &
       &    Node_Component_Spheroid_Standard_Inactive           , Node_Component_Spheroid_Standard_Post_Step                    , &
       &    Node_Component_Spheroid_Standard_Initialize

  !# <component>
  !#  <class>spheroid</class>
  !#  <name>standard</name>
  !#  <isDefault>true</isDefault>
  !#  <createFunction isDeferred="true" />
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
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of stars in the standard spheroid."/>
  !#   </property>
  !#   <property>
  !#     <name>massStellarFormed</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#   </property>
  !#   <property>
  !#     <name>abundancesStellar</name>
  !#     <type>abundances</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of metals in the stellar phase of the standard spheroid."/>
  !#   </property>
  !#   <property>
  !#     <name>massGas</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of gas in the standard spheroid."/>
  !#   </property>
  !#   <property>
  !#     <name>abundancesGas</name>
  !#     <type>abundances</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of metals in the gas phase of the standard spheroid."/>
  !#   </property>
  !#   <property>
  !#     <name>angularMomentum</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#     <output unitsInSI="massSolar*megaParsec*kilo" comment="Angular momentum of the standard spheroid."/>
  !#   </property>
  !#   <property>
  !#     <name>radius</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <output unitsInSI="megaParsec" comment="Scale length of the standard spheroid."/>
  !#   </property>
  !#   <property>
  !#     <name>halfMassRadius</name>
  !#     <attributes isSettable="false" isGettable="true" isEvolvable="false" isVirtual="true" />
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <getFunction>Node_Component_Spheroid_Standard_Half_Mass_Radius</getFunction>
  !#   </property>
  !#   <property>
  !#     <name>velocity</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <output unitsInSI="kilo" comment="Circular velocity at the scale length of the standard spheroid."/>
  !#   </property>
  !#   <property>
  !#     <name>luminositiesStellar</name>
  !#     <type>stellarLuminosities</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#     <output unitsInSI="luminosityZeroPointAB" comment="Luminosity of spheroid stars."/>
  !#   </property>
  !#   <property>
  !#     <name>stellarPropertiesHistory</name>
  !#     <type>history</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" isDeferred="rate" createIfNeeded="true" />
  !#   </property>
  !#   <property>
  !#     <name>starFormationHistory</name>
  !#     <type>history</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" isDeferred="rate" createIfNeeded="true" />
  !#   </property>
  !#   <property>
  !#     <name>massGasSink</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="false" isGettable="false" isEvolvable="true" isDeferred="rate" isVirtual="true" />
  !#   </property>
  !#   <property>
  !#     <name>energyGasInput</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="false" isGettable="false" isEvolvable="true" isDeferred="rate" isVirtual="true" />
  !#   </property>
  !#  </properties>
  !#  <bindings>
  !#   <binding method="enclosedMass"          function="Node_Component_Spheroid_Standard_Enclosed_Mass"           bindsTo="component" />
  !#   <binding method="acceleration"          function="Node_Component_Spheroid_Standard_Acceleration"            bindsTo="component" />
  !#   <binding method="tidalTensor"           function="Node_Component_Spheroid_Standard_Tidal_Tensor"            bindsTo="component" />
  !#   <binding method="chandrasekharIntegral" function="Node_Component_Spheroid_Standard_Chandrasekhar_Integral"  bindsTo="component" />
  !#   <binding method="density"               function="Node_Component_Spheroid_Standard_Density"                 bindsTo="component" />
  !#   <binding method="rotationCurve"         function="Node_Component_Spheroid_Standard_Rotation_Curve"          bindsTo="component" />
  !#   <binding method="rotationCurveGradient" function="Node_Component_Spheroid_Standard_Rotation_Curve_Gradient" bindsTo="component" />
  !#   <binding method="potential"             function="Node_Component_Spheroid_Standard_Potential"               bindsTo="component" />
  !#  </bindings>
  !#  <functions>objects.nodes.components.spheroid.standard.bound_functions.inc</functions>
  !# </component>

  ! Objects used by this component.
  class(stellarPopulationPropertiesClass            ), pointer :: stellarPopulationProperties_
  class(ramPressureStrippingSpheroidsClass          ), pointer :: ramPressureStrippingSpheroids_
  class(tidalStrippingSpheroidsClass                ), pointer :: tidalStrippingSpheroids_
  class(darkMatterHaloScaleClass                    ), pointer :: darkMatterHaloScale_
  class(satelliteTidalFieldClass                    ), pointer :: satelliteTidalField_
  class(starFormationHistoryClass                   ), pointer :: starFormationHistory_
  class(mergerMassMovementsClass                    ), pointer :: mergerMassMovements_
  class(mergerRemnantSizeClass                      ), pointer :: mergerRemnantSize_
  !$omp threadprivate(stellarPopulationProperties_,ramPressureStrippingSpheroids_,tidalStrippingSpheroids_,darkMatterHaloScale_,satelliteTidalField_,starFormationHistory_,mergerMassMovements_,mergerRemnantSize_)

  ! Internal count of abundances.
  integer                                     :: abundancesCount

  ! Storage for the star formation and stellar properties histories time range, used when extending this range.
  double precision, allocatable, dimension(:) :: starFormationHistoryTemplate, stellarPropertiesHistoryTemplate
  !$omp threadprivate(starFormationHistoryTemplate,stellarPropertiesHistoryTemplate)
  ! Parameters controlling the physical implementation.
  double precision                            :: spheroidEnergeticOutflowMassRate    , spheroidMassToleranceAbsolute
  logical                                     :: spheroidLuminositiesStellarInactive

  ! Spheroid structural parameters.
  double precision                            :: spheroidAngularMomentumAtScaleRadius

  ! Minimum absolute scales for physically plausible spheroids.
  double precision, parameter                 :: radiusMinimum                       =1.0d-12
  double precision, parameter                 :: massMinimum                         =1.0d-06
  double precision, parameter                 :: angularMomentumMinimum              =1.0d-20

contains

  !# <nodeComponentInitializationTask>
  !#  <unitName>Node_Component_Spheroid_Standard_Initialize</unitName>
  !# </nodeComponentInitializationTask>
  subroutine Node_Component_Spheroid_Standard_Initialize(parameters_)
    !% Initializes the tree node standard spheroid methods module.
    use :: Abundances_Structure, only : Abundances_Property_Count
    use :: Galacticus_Error    , only : Galacticus_Error_Report
    use :: Galacticus_Nodes    , only : defaultSpheroidComponent , nodeComponentSpheroidStandard
    use :: Input_Parameters    , only : inputParameter           , inputParameters
    implicit none
    type(inputParameters              ), intent(inout) :: parameters_
    type(nodeComponentSpheroidStandard)                :: spheroidStandardComponent

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

       ! Read parameters controlling the physical implementation.
       !# <inputParameter>
       !#   <name>spheroidEnergeticOutflowMassRate</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultValue>1.0d-2</defaultValue>
       !#   <description>The proportionallity factor relating mass outflow rate from the spheroid to the energy input rate divided by $V_\mathrm{spheroid}^2$.</description>
       !#   <source>parameters_</source>
       !#   <type>double</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>spheroidMassToleranceAbsolute</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultValue>1.0d-6</defaultValue>
       !#   <description>The mass tolerance used to judge whether the spheroid is physically plausible.</description>
       !#   <source>parameters_</source>
       !#   <type>double</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>spheroidLuminositiesStellarInactive</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultValue>.false.</defaultValue>
       !#   <description>Specifies whether or not spheroid stellar luminosities are inactive properties (i.e. do not appear in any ODE being solved).</description>
       !#   <source>parameters_</source>
       !#   <type>boolean</type>
       !# </inputParameter>
    end if
    return
  end subroutine Node_Component_Spheroid_Standard_Initialize

  !# <nodeComponentThreadInitializationTask>
  !#  <unitName>Node_Component_Spheroid_Standard_Thread_Initialize</unitName>
  !# </nodeComponentThreadInitializationTask>
  subroutine Node_Component_Spheroid_Standard_Thread_Initialize(parameters_)
    !% Initializes the standard spheroid module for each thread.
    use :: Events_Hooks                         , only : satelliteMergerEvent             , postEvolveEvent, openMPThreadBindingAtLevel, dependencyRegEx, &
         &                                               dependencyDirectionAfter
    use :: Galacticus_Error                     , only : Galacticus_Error_Report
    use :: Galacticus_Nodes                     , only : defaultSpheroidComponent
    use :: Input_Parameters                     , only : inputParameter                   , inputParameters
    use :: Mass_Distributions                   , only : massDistributionSymmetrySpherical
    use :: Node_Component_Spheroid_Standard_Data, only : spheroidMassDistribution
    implicit none
    type            (inputParameters), intent(inout) :: parameters_
    logical                                          :: densityMoment2IsInfinite                   , densityMoment3IsInfinite
    double precision                                 :: spheroidMassDistributionDensityMomentum2   , spheroidMassDistributionDensityMomentum3, &
         &                                              spheroidAngularMomentumAtScaleRadiusDefault
    type            (dependencyRegEx), dimension(2)  :: dependencies

    ! Check if this implementation is selected. If so, initialize the mass distribution.
    if (defaultSpheroidComponent%standardIsActive()) then
       call postEvolveEvent     %attach(defaultSpheroidComponent,postEvolve     ,openMPThreadBindingAtLevel,label='nodeComponentSpheroidStandard'                          )
       dependencies(1)=dependencyRegEx(dependencyDirectionAfter,'^remnantStructure:')
       dependencies(2)=dependencyRegEx(dependencyDirectionAfter,'^nodeComponentDisk')
       call satelliteMergerEvent%attach(defaultSpheroidComponent,satelliteMerger,openMPThreadBindingAtLevel,label='nodeComponentSpheroidStandard',dependencies=dependencies)
       !# <objectBuilder class="stellarPopulationProperties"                                            name="stellarPopulationProperties_"   source="parameters_"                    />
       !# <objectBuilder class="ramPressureStrippingSpheroids"                                          name="ramPressureStrippingSpheroids_" source="parameters_"                    />
       !# <objectBuilder class="tidalStrippingSpheroids"                                                name="tidalStrippingSpheroids_"       source="parameters_"                    />
       !# <objectBuilder class="darkMatterHaloScale"                                                    name="darkMatterHaloScale_"           source="parameters_"                    />
       !# <objectBuilder class="satelliteTidalField"                                                    name="satelliteTidalField_"           source="parameters_"                    />
       !# <objectBuilder class="starFormationHistory"                                                   name="starFormationHistory_"          source="parameters_"                    />
       !# <objectBuilder class="mergerMassMovements"                                                    name="mergerMassMovements_"           source="parameters_"                    />
       !# <objectBuilder class="mergerRemnantSize"                                                      name="mergerRemnantSize_"             source="parameters_"                    />
       !# <objectBuilder class="massDistribution"              parameterName="spheroidMassDistribution" name="spheroidMassDistribution"       source="parameters_" threadPrivate="yes" >
       !#  <default>
       !#   <spheroidMassDistribution value="hernquist">
       !#    <dimensionless value="true"/>
       !#   </spheroidMassDistribution>
       !#  </default>
       !# </objectBuilder>
       if (.not.spheroidMassDistribution%isDimensionless()                                     ) &
            & call Galacticus_Error_Report('spheroid mass distribution must be dimensionless'        //{introspection:location})
       if (.not.spheroidMassDistribution%symmetry       () == massDistributionSymmetrySpherical) &
            & call Galacticus_Error_Report('spheroid mass distribution must be spherically symmetric'//{introspection:location})
       ! Determine the specific angular momentum at the scale radius in units of the mean specific angular
       ! momentum of the spheroid. This is equal to the ratio of the 2nd to 3rd radial moments of the density
       ! distribution (assuming a flat rotation curve).
       spheroidMassDistributionDensityMomentum2=spheroidMassDistribution%densityRadialMoment(2.0d0,isInfinite=densityMoment2IsInfinite)
       spheroidMassDistributionDensityMomentum3=spheroidMassDistribution%densityRadialMoment(3.0d0,isInfinite=densityMoment3IsInfinite)
       if (densityMoment2IsInfinite.or.densityMoment3IsInfinite) then
          ! One of the moments is infinte, so we can not compute the appropriate ratio. Simply assume a value
          ! of 0.5 as a default.
          spheroidAngularMomentumAtScaleRadiusDefault=+0.5d0
       else
          ! Moments are well-defined, so compute their ratio.
          spheroidAngularMomentumAtScaleRadiusDefault=+spheroidMassDistributionDensityMomentum2 &
               &                                      /spheroidMassDistributionDensityMomentum3
       end if
       !$omp critical (spheroidStandardInitializeAngularMomentum)
       !# <inputParameter>
       !#   <name>spheroidAngularMomentumAtScaleRadius</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultSource>($I_2/I_3$ where $I_n=\int_0^\infty \rho(r) r^n \mathrm{d}r$, where $\rho(r)$ is the spheroid density profile, unless either $I_2$ or $I_3$ is infinite, in which case a default of $1/2$ is used instead.)</defaultSource>
       !#   <defaultValue>spheroidAngularMomentumAtScaleRadiusDefault</defaultValue>
       !#   <description>The assumed ratio of the specific angular momentum at the scale radius to the mean specific angular momentum of the standard spheroid component.</description>
       !#   <source>parameters_</source>
       !#   <type>double</type>
       !# </inputParameter>
       !$omp end critical (spheroidStandardInitializeAngularMomentum)
    end if
    return
  end subroutine Node_Component_Spheroid_Standard_Thread_Initialize

  !# <nodeComponentThreadUninitializationTask>
  !#  <unitName>Node_Component_Spheroid_Standard_Thread_Uninitialize</unitName>
  !# </nodeComponentThreadUninitializationTask>
  subroutine Node_Component_Spheroid_Standard_Thread_Uninitialize()
    !% Uninitializes the standard spheroid module for each thread.
    use :: Events_Hooks                         , only : satelliteMergerEvent    , postEvolveEvent
    use :: Galacticus_Nodes                     , only : defaultSpheroidComponent
    use :: Node_Component_Spheroid_Standard_Data, only : spheroidMassDistribution
    implicit none

    if (defaultSpheroidComponent%standardIsActive()) then
       call postEvolveEvent     %detach(defaultSpheroidComponent,postEvolve     )
       call satelliteMergerEvent%detach(defaultSpheroidComponent,satelliteMerger)
       !# <objectDestructor name="stellarPopulationProperties_"  />
       !# <objectDestructor name="ramPressureStrippingSpheroids_"/>
       !# <objectDestructor name="tidalStrippingSpheroids_"      />
       !# <objectDestructor name="darkMatterHaloScale_"          />
       !# <objectDestructor name="satelliteTidalField_"          />
       !# <objectDestructor name="starFormationHistory_"         />
       !# <objectDestructor name="mergerMassMovements_"          />
       !# <objectDestructor name="mergerRemnantSize_"            />
       !# <objectDestructor name="spheroidMassDistribution"      />
    end if
    return
  end subroutine Node_Component_Spheroid_Standard_Thread_Uninitialize

  !# <preEvolveTask>
  !# <unitName>Node_Component_Spheroid_Standard_Pre_Evolve</unitName>
  !# </preEvolveTask>
  subroutine Node_Component_Spheroid_Standard_Pre_Evolve(node)
    !% Ensure the spheroid has been initialized.
    use :: Galacticus_Nodes, only : nodeComponentSpheroid, nodeComponentSpheroidStandard, treeNode, defaultSpheroidComponent
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
    !% Trim histories attached to the spheroid.
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

  !# <postStepTask>
  !# <unitName>Node_Component_Spheroid_Standard_Post_Step</unitName>
  !# <after>Node_Component_Basic_Standard_Post_Step</after>
  !# </postStepTask>
  subroutine Node_Component_Spheroid_Standard_Post_Step(node,status)
    !% Trim histories attached to the spheroid.
    use :: Abundances_Structure          , only : abs                       , zeroAbundances
    use :: Galacticus_Display            , only : Galacticus_Display_Message, verbosityWarn
    use :: Galacticus_Error              , only : Galacticus_Error_Report
    use :: Galacticus_Nodes              , only : nodeComponentSpheroid     , nodeComponentSpheroidStandard, treeNode      , defaultSpheroidComponent
    use :: Interface_GSL                 , only : GSL_Failure
    use :: ISO_Varying_String            , only : assignment(=)             , operator(//)                 , varying_string
    use :: Stellar_Luminosities_Structure, only : stellarLuminosities       , abs
    use :: String_Handling               , only : operator(//)
    implicit none
    type            (treeNode             ), intent(inout), pointer :: node
    integer                                , intent(inout)          :: status
    class           (nodeComponentSpheroid)               , pointer :: spheroid
    double precision                       , parameter              :: angularMomentumTolerance=1.0d-2
    double precision                       , parameter              :: massTolerance           =1.0d+0
    double precision                       , save                   :: fractionalErrorMaximum  =0.0d+0
    double precision                                                :: fractionalError                , specificAngularMomentum, &
         &                                                             spheroidMass
    character       (len=20               )                         :: valueString
    type            (varying_string       ), save                   :: message
    !$omp threadprivate(message)
    type            (stellarLuminosities  ), save                   :: luminositiesStellar
    !$omp threadprivate(luminositiesStellar)

    ! Return immediately if this class is not in use.
    if (.not.defaultSpheroidComponent%standardIsActive()) return
    ! Get the spheroid component.
    spheroid => node%spheroid()
    ! Check if an standard spheroid component exists.
    select type (spheroid)
    class is (nodeComponentSpheroidStandard)
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
             call Galacticus_Display_Message(message,verbosityWarn)
             ! Store the new maximum fractional error.
             fractionalErrorMaximum=fractionalError
          end if
          !$omp end critical (Standard_Spheroid_Post_Evolve_Check)
          ! Get the specific angular momentum of the spheroid material
          spheroidMass= spheroid%massGas    () &
               &       +spheroid%massStellar()
          if (spheroidMass == 0.0d0) then
             specificAngularMomentum=0.0d0
             call spheroid%        massStellarSet(                  0.0d0)
             call spheroid%  abundancesStellarSet(         zeroAbundances)
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
             specificAngularMomentum=spheroid%angularMomentum()/spheroidMass
          end if
          ! Reset the gas, abundances and angular momentum of the spheroid.
          call spheroid%        massGasSet(                                                      0.0d0)
          call spheroid%  abundancesGasSet(                                             zeroAbundances)
          call spheroid%angularMomentumSet(specificAngularMomentum*spheroid%massStellar())
          status=GSL_Failure
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
             call Galacticus_Display_Message(message,verbosityWarn)
             ! Store the new maximum fractional error.
             fractionalErrorMaximum=fractionalError
          end if
          !$omp end critical (Standard_Spheroid_Post_Evolve_Check)
          ! Get the specific angular momentum of the spheroid material
          spheroidMass= spheroid%massGas    () &
               &       +spheroid%massStellar()
          if (spheroidMass == 0.0d0) then
             specificAngularMomentum=0.0d0
             call spheroid%      massGasSet(         0.0d0)
             call spheroid%abundancesGasSet(zeroAbundances)
          else
             specificAngularMomentum=spheroid%angularMomentum()/spheroidMass
          end if
          ! Reset the stellar, abundances and angular momentum of the spheroid.
          call spheroid%        massStellarSet(                                 0.0d0)
          call spheroid%  abundancesStellarSet(                        zeroAbundances)
          call spheroid%angularMomentumSet(specificAngularMomentum*spheroid%massGas())
          status=GSL_Failure
       end if
    end select
    return
  end subroutine Node_Component_Spheroid_Standard_Post_Step

  subroutine Node_Component_Spheroid_Standard_Mass_Gas_Sink_Rate(self,rate,interrupt,interruptProcedure)
    !% Account for a sink of gaseous material in the standard spheroid.
    use :: Abundances_Structure, only : operator(*)
    use :: Galacticus_Error    , only : Galacticus_Error_Report
    use :: Galacticus_Nodes    , only : interruptTask          , nodeComponentSpheroid
    implicit none
    class           (nodeComponentSpheroid), intent(inout)                    :: self
    logical                                , intent(inout), optional          :: interrupt
    procedure       (interruptTask        ), intent(inout), optional, pointer :: interruptProcedure
    double precision                       , intent(in   )                    :: rate
    double precision                                                          :: gasMass           , stellarMass
    !$GLC attributes unused :: interrupt, interruptProcedure

    ! Trap cases where an attempt is made to add gas via this sink function.
    if (rate > 0.0d0) call Galacticus_Error_Report(                                                      &
         &                                         'attempt to add mass via sink in standard spheroid'// &
         &                                         {introspection:location}                              &
         &                                        )
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
    !% Handles input of energy into the spheroid gas from other components (e.g. black holes). The energy input rate should be in
    !% units of $M_\odot$ km$^2$ s$^{-2}$ Gyr$^{-1}$.
    use :: Abundances_Structure, only : abundances             , operator(*)
    use :: Galacticus_Error    , only : Galacticus_Error_Report
    use :: Galacticus_Nodes    , only : interruptTask          , nodeComponentHotHalo, nodeComponentSpheroid, treeNode
    implicit none
    class           (nodeComponentSpheroid        ), intent(inout)                    :: self
    logical                                        , intent(inout), optional          :: interrupt
    procedure       (interruptTask                ), intent(inout), optional, pointer :: interruptProcedure
    double precision                               , intent(in   )                    :: rate
    class           (nodeComponentHotHalo         )                         , pointer :: selfHotHaloComponent
    type            (treeNode                     )                         , pointer :: selfNode
    type            (abundances                   )                                   :: abundancesOutflowRate
    double precision                                                                  :: angularMomentumOutflowRate, gasMass         , &
         &                                                                               massOutflowRate           , spheroidVelocity, &
         &                                                                               stellarMass
    !$GLC attributes unused :: interrupt, interruptProcedure

    ! Trap cases where an attempt is made to remove energy via this input function.
    if (rate < 0.0d0) call Galacticus_Error_Report(                                                                 &
         &                                         'attempt to remove energy via input pipe in standard spheroid'// &
         &                                         {introspection:location}                                         &
         &                                        )

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
          massOutflowRate           =spheroidEnergeticOutflowMassRate*rate/spheroidVelocity**2
          angularMomentumOutflowRate=(massOutflowRate/(gasMass+stellarMass))*self%angularMomentum()
          abundancesOutflowRate     =(massOutflowRate/ gasMass             )*self%abundancesGas  ()
          call self%        massGasRate(-           massOutflowRate)
          call self%angularMomentumRate(-angularMomentumOutflowRate)
          call self%  abundancesGasRate(-     abundancesOutflowRate)
          ! Add outflowing rates to the hot halo component.
          selfNode             => self    %host   ()
          selfHotHaloComponent => selfNode%hotHalo()
          call selfHotHaloComponent%           outflowingMassRate(           massOutflowRate)
          call selfHotHaloComponent%outflowingAngularMomentumRate(angularMomentumOutflowRate)
          call selfHotHaloComponent%     outflowingAbundancesRate(     abundancesOutflowRate)
       end if
    end if
    return
  end subroutine Node_Component_Spheroid_Standard_Energy_Gas_Input_Rate

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Spheroid_Standard_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Spheroid_Standard_Rate_Compute(node,interrupt,interruptProcedure,propertyType)
    !% Compute the standard spheroid node mass rate of change.
    use :: Abundances_Structure            , only : abs                     , abundances          , max                  , operator(*)
    use :: Galacticus_Error                , only : Galacticus_Error_Report
    use :: Galacticus_Nodes                , only : defaultSpheroidComponent, nodeComponentHotHalo, nodeComponentSpheroid, nodeComponentSpheroidStandard, &
         &                                          propertyTypeActive      , propertyTypeAll     , propertyTypeInactive , treeNode
    use :: Histories                       , only : history                 , operator(*)
    use :: Numerical_Constants_Astronomical, only : Mpc_per_km_per_s_To_Gyr
    use :: Stellar_Luminosities_Structure  , only : abs                     , max                 , operator(*)          , stellarLuminosities          , &
         &                                          zeroStellarLuminosities
    implicit none
    type            (treeNode             ), intent(inout), pointer :: node
    logical                                , intent(inout)          :: interrupt
    procedure       (                     ), intent(inout), pointer :: interruptProcedure
    integer                                , intent(in   )          :: propertyType
    class           (nodeComponentSpheroid)               , pointer :: spheroid
    class           (nodeComponentHotHalo )               , pointer :: hotHalo
    double precision                                                :: fractionGas             , fractionStellar, &
         &                                                             tidalField              , tidalTorque    , &
         &                                                             massLossRate
    type            (history              )                         :: historyTransferRate
    type            (stellarLuminosities  ), save                   :: luminositiesStellarRates
    !$omp threadprivate(luminositiesStellarRates)
    !$GLC attributes unused :: interrupt, interruptProcedure

    ! Return immediately if this class is not in use or only inactive properties are to be computed.
    if (.not.defaultSpheroidComponent%standardIsActive() .or. propertyType == propertyTypeInactive) return
    ! Get the disk and check that it is of our class.
    spheroid => node%spheroid()
    select type (spheroid)
    class is (nodeComponentSpheroidStandard)
       ! Check for a realistic spheroid, return immediately if spheroid is unphysical.
       if     (    spheroid%angularMomentum() < angularMomentumMinimum &
            & .or. spheroid%radius         () <          radiusMinimum &
            & .or. spheroid%massGas        () <            massMinimum &
            & ) return
       ! Apply mass loss rate due to ram pressure stripping.
       if (spheroid%massGas() > 0.0d0) then
          massLossRate=ramPressureStrippingSpheroids_%rateMassLoss(node)
          if (massLossRate > 0.0d0) then
             hotHalo => node%hotHalo()
             call spheroid%                  massGasRate(-massLossRate                                                                       )
             call spheroid%          angularMomentumRate(-massLossRate*spheroid%angularMomentum()/(spheroid%massGas()+spheroid%massStellar()))
             call spheroid%            abundancesGasRate(-massLossRate*spheroid%abundancesGas  ()/ spheroid%massGas()                        )
             call hotHalo %           outflowingMassRate(+massLossRate                                                                       )
             call hotHalo %outflowingAngularMomentumRate(+massLossRate*spheroid%angularMomentum()/(spheroid%massGas()+spheroid%massStellar()))
             call hotHalo %     outflowingAbundancesRate(+massLossRate*spheroid%abundancesGas  ()/ spheroid%massGas()                        )
          end if
       end if
       ! Apply mass loss rate due to tidal stripping.
       if (spheroid%massGas()+spheroid%massStellar() > 0.0d0) then
          massLossRate=tidalStrippingSpheroids_%rateMassLoss(node)
          if (massLossRate > 0.0d0) then
             hotHalo    => node%hotHalo()
             fractionGas    =  min(1.0d0,max(0.0d0,spheroid%massGas()/(spheroid%massGas()+spheroid%massStellar())))
             fractionStellar=  1.0d0-fractionGas
             if (fractionGas    > 0.0d0 .and. spheroid%massGas    () > 0.0d0) then
                call spheroid%                  massGasRate  (-fractionGas    *massLossRate                                                                         )
                call spheroid%            abundancesGasRate  (-fractionGas    *massLossRate*spheroid%abundancesGas    ()/ spheroid%massGas()                        )
                call hotHalo %           outflowingMassRate  (+fractionGas    *massLossRate                                                                         )
                call hotHalo %     outflowingAbundancesRate  (+fractionGas    *massLossRate*spheroid%abundancesGas    ()/ spheroid%massGas()                        )
                call hotHalo %outflowingAngularMomentumRate  (+fractionGas    *massLossRate*spheroid%angularMomentum  ()/(spheroid%massGas()+spheroid%massStellar()))
             end if
             if (fractionStellar > 0.0d0 .and. spheroid%massStellar() > 0.0d0) then
                ! If luminosities are being treated as inactive properties this is an error - they appear on the right-hand side
                ! of the following ODE terms so are not inactive. (An approach similar to what is used for transfer of
                ! luminosities to the spheroid by bar instabilities in the standard disk component could work here.)
                if (propertyType == propertyTypeActive .and. spheroidLuminositiesStellarInactive) call Galacticus_Error_Report('tidal mass loss not supported for inactive luminosity calculation'//{introspection:location})
                call spheroid%              massStellarRate  (-fractionStellar*massLossRate                                                                         )
                call spheroid%        abundancesStellarRate  (-fractionStellar*massLossRate*spheroid%abundancesStellar()/                    spheroid%massStellar() )
                luminositiesStellarRates=max(zeroStellarLuminosities,spheroid%luminositiesStellar())
                call spheroid%      luminositiesStellarRate(-fractionStellar*massLossRate*luminositiesStellarRates      /                    spheroid%massStellar() )
                ! Stellar properties history.
                historyTransferRate=spheroid%stellarPropertiesHistory()
                if (historyTransferRate%exists()) then
                   call spheroid%stellarPropertiesHistoryRate(-fractionStellar*massLossRate*historyTransferRate         /                    spheroid%massStellar() )
                end if
                call historyTransferRate%destroy()
                ! Star formation history.
                historyTransferRate=spheroid%starFormationHistory()
                if (historyTransferRate%exists()) then
                   call spheroid%starFormationHistoryRate     (-fractionStellar*massLossRate*historyTransferRate        /                    spheroid%massStellar() )
                end if
                call historyTransferRate%destroy()
             end if
             call       spheroid%          angularMomentumRate(-                massLossRate*spheroid%angularMomentum ()/(spheroid%massGas()+spheroid%massStellar()))
          end if
       end if

       ! Apply tidal heating.
       if (node%isSatellite() .and. spheroid%angularMomentum() < (spheroid%massGas()+spheroid%massStellar())*darkMatterHaloScale_%virialRadius(node)*darkMatterHaloScale_%virialVelocity(node) .and. spheroid%radius() < darkMatterHaloScale_%virialRadius(node)) then
          tidalField =satelliteTidalField_%tidalTensorRadial(node)
          tidalTorque=abs(tidalField)*(spheroid%massGas()+spheroid%massStellar())*spheroid%radius()**2
          call spheroid%angularMomentumRate(+tidalTorque)
       end if

    end select
    return
  end subroutine Node_Component_Spheroid_Standard_Rate_Compute

  subroutine Node_Component_Spheroid_Standard_Star_Formation_History_Rate(self,rate,interrupt,interruptProcedure)
    !% Adjust the rates for the star formation history.
    use :: Galacticus_Error , only : Galacticus_Error_Report
    use :: Galacticus_Nodes , only : interruptTask          , nodeComponentSpheroid, nodeComponentSpheroidStandard
    use :: Memory_Management, only : allocateArray          , deallocateArray
    implicit none
    class    (nodeComponentSpheroid), intent(inout)                    :: self
    type     (history              ), intent(in   )                    :: rate
    logical                         , intent(inout), optional          :: interrupt
    procedure(interruptTask        ), intent(inout), optional, pointer :: interruptProcedure
    type     (history              )                                   :: starFormationHistory

    ! Get the star formation history in the spheroid.
    starFormationHistory=self%starFormationHistory()
    ! Ensure that the history already exists.
    if (.not.starFormationHistory%exists())                                                        &
         & call Galacticus_Error_Report(                                                           &
         &                              'no star formation history has been created in spheroid'// &
         &                              {introspection:location}                                   &
         &                             )
    ! Check if the star formation history in the spheroid spans a sufficient range to accept the input rates.
    if     (                                                                                             &
         &       rate%time(              1) < starFormationHistory%time(                              1) &
         &  .or. rate%time(size(rate%time)) > starFormationHistory%time(size(starFormationHistory%time)) &
         & ) then
       ! It does not, so interrupt evolution and extend the history.
       if (allocated(starFormationHistoryTemplate)) call deallocateArray(starFormationHistoryTemplate)
       call allocateArray(starFormationHistoryTemplate,shape(rate%time))
       starFormationHistoryTemplate=rate%time
       interrupt=.true.
       interruptProcedure => Node_Component_Spheroid_Standard_Star_Formation_History_Extend
       return
    end if
    ! Adjust the rate.
    select type (self)
    class is (nodeComponentSpheroidStandard)
       call self%starFormationHistoryRateIntrinsic(rate)
    end select
    return
  end subroutine Node_Component_Spheroid_Standard_Star_Formation_History_Rate

  subroutine Node_Component_Spheroid_Standard_Stellar_Prprts_History_Rate(self,rate,interrupt,interruptProcedure)
    !% Adjust the rates for the stellar properties history.
    use :: Galacticus_Error , only : Galacticus_Error_Report
    use :: Galacticus_Nodes , only : interruptTask          , nodeComponentSpheroid, nodeComponentSpheroidStandard
    use :: Memory_Management, only : allocateArray          , deallocateArray
    implicit none
    class    (nodeComponentSpheroid), intent(inout)                    :: self
    type     (history              ), intent(in   )                    :: rate
    logical                         , intent(inout), optional          :: interrupt
    procedure(interruptTask        ), intent(inout), optional, pointer :: interruptProcedure
    type     (history              )                                   :: stellarPropertiesHistory

    ! Get the star formation history in the spheroid.
    stellarPropertiesHistory=self%stellarPropertiesHistory()
    ! Ensure that the history already exists.
    if (.not.stellarPropertiesHistory%exists())                                                    &
         & call Galacticus_Error_Report(                                                           &
         &                              'no star formation history has been created in spheroid'// &
         &                              {introspection:location}                                   &
         &                             )
    ! Check if the star formation history in the spheroid spans a sufficient range to accept the input rates.
    if     (                                                                                                     &
         &       rate%time(              1) < stellarPropertiesHistory%time(                                  1) &
         &  .or. rate%time(size(rate%time)) > stellarPropertiesHistory%time(size(stellarPropertiesHistory%time)) &
         & ) then
       ! It does not, so interrupt evolution and extend the history.
       if (allocated(stellarPropertiesHistoryTemplate)) call deallocateArray(stellarPropertiesHistoryTemplate)
       call allocateArray(stellarPropertiesHistoryTemplate,shape(rate%time))
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

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Spheroid_Standard_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Spheroid_Standard_Scale_Set(node)
    !% Set scales for properties of {\normalfont \ttfamily node}. Note that gas masses get an additional scaling down since they can approach
    !% zero and we'd like to prevent them from becoming negative.
    use :: Abundances_Structure          , only : abs                     , max                  , operator(*)                  , unitAbundances
    use :: Galacticus_Nodes              , only : nodeComponentDisk       , nodeComponentSpheroid, nodeComponentSpheroidStandard, treeNode           , &
         &                                        defaultSpheroidComponent
    use :: Histories                     , only : history                 , operator(*)
    use :: Stellar_Luminosities_Structure, only : abs                     , max                  , operator(*)                  , stellarLuminosities, &
          &                                       unitStellarLuminosities
    implicit none
    type            (treeNode             ), intent(inout), pointer :: node
    class           (nodeComponentSpheroid)               , pointer :: spheroid
    class           (nodeComponentDisk    )               , pointer :: disk
    double precision                       , parameter              :: massMinimum                   =1.0d0
    double precision                       , parameter              :: angularMomentumMinimum        =0.1d0
    double precision                       , parameter              :: gasMassScaling                =0.1d0
    double precision                       , parameter              :: luminosityMinimum             =1.0d0
    double precision                                                :: angularMomentum                     , mass
    type            (history              )                         :: stellarPopulationHistoryScales
    type            (stellarLuminosities  )                         :: stellarLuminositiesScale

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
       angularMomentum=abs(spheroid%angularMomentum())
       call spheroid%angularMomentumScale  (               max(angularMomentum,angularMomentumMinimum))

       ! Set scale for gas mass.
       mass           =abs(                        &
            &              +spheroid%massGas    () &
            &              +spheroid%massStellar() &
            &             )
       call spheroid%          massGasScale(gasMassScaling*max(           mass,           massMinimum))
       call spheroid%      massStellarScale(               max(           mass,           massMinimum))
       call spheroid%massStellarFormedScale(               max(           mass,           massMinimum))

       ! Set scales for abundances if necessary.
       if (abundancesCount > 0) then
          ! Set scale for gas abundances.
          call spheroid%abundancesGasScale    (                                                       &
               &                                            +gasMassScaling                           &
               &                                            *max(                                     &
               &                                                 +abs(+spheroid%abundancesGas    ())  &
               &                                                 +abs(+spheroid%abundancesStellar()), &
               &                                                      +massMinimum                    &
               &                                                      *unitAbundances                 &
               &                                                )                                     &
               &                                           )

          ! Set scale for stellar abundances.
          call spheroid%abundancesStellarScale(                                                      &
               &                                            max(                                     &
               &                                                 abs(+spheroid%abundancesStellar()), &
               &                                                     +massMinimum                    &
               &                                                     *unitAbundances                 &
               &                                               )                                     &
               &                                           )
       end if
       ! Set scales for stellar luminosities.
       stellarLuminositiesScale=max(                                       &
            &                       abs(spheroid  %luminositiesStellar()), &
            &                           +unitStellarLuminosities           &
            &                           *luminosityMinimum                 &
            &                      )
       call stellarLuminositiesScale%truncate                (spheroid   %luminositiesStellar())
       call spheroid   %luminositiesStellarScale(stellarLuminositiesScale                      )

       ! Set scales for stellar population properties history.
       stellarPopulationHistoryScales=spheroid%stellarPropertiesHistory()
       call stellarPopulationProperties_%scales   (spheroid%massStellar(),spheroid%abundancesStellar(),stellarPopulationHistoryScales)
       call spheroid%stellarPropertiesHistoryScale(                                                    stellarPopulationHistoryScales)
       call stellarPopulationHistoryScales%destroy()
       stellarPopulationHistoryScales=spheroid%starFormationHistory()
       call starFormationHistory_%scales          (stellarPopulationHistoryScales,spheroid%massStellar(),spheroid%abundancesStellar())
       call spheroid%starFormationHistoryScale    (stellarPopulationHistoryScales                                                    )
       call stellarPopulationHistoryScales%destroy()
    end select
    return
  end subroutine Node_Component_Spheroid_Standard_Scale_Set

  !# <inactiveSetTask>
  !#  <unitName>Node_Component_Spheroid_Standard_Inactive</unitName>
  !# </inactiveSetTask>
  subroutine Node_Component_Spheroid_Standard_Inactive(node)
    !% Set Jacobian zero status for properties of {\normalfont \ttfamily node}.
    use :: Galacticus_Nodes, only : nodeComponentSpheroid, nodeComponentSpheroidStandard, treeNode
    implicit none
    type (treeNode             ), intent(inout), pointer :: node
    class(nodeComponentSpheroid)               , pointer :: spheroid

    ! Get the spheroid component.
    spheroid => node%spheroid()
    ! Check if an standard spheroid component exists.
    select type (spheroid)
    class is (nodeComponentSpheroidStandard)
       if (spheroidLuminositiesStellarInactive) call spheroid%luminositiesStellarInactive()
    end select
    return
  end subroutine Node_Component_Spheroid_Standard_Inactive

  subroutine satelliteMerger(self,node)
    !% Transfer any standard spheroid associated with {\normalfont \ttfamily node} to its host halo.
    use :: Abundances_Structure            , only : zeroAbundances
    use :: Galacticus_Error                , only : Galacticus_Error_Report
    use :: Galacticus_Nodes                , only : treeNode                , nodeComponentDisk        , nodeComponentSpheroid   , nodeComponentSpheroidStandard
    use :: Satellite_Merging_Mass_Movements, only : destinationMergerDisk   , destinationMergerSpheroid, destinationMergerUnmoved
    use :: Satellite_Merging_Remnant_Sizes , only : remnantNoChange
    use :: Stellar_Luminosities_Structure  , only : zeroStellarLuminosities
    implicit none
    class           (*                    ), intent(inout) :: self
    type            (treeNode             ), intent(inout) :: node
    type            (treeNode             ), pointer       :: nodeHost
    class           (nodeComponentDisk    ), pointer       :: diskHost                      , disk
    class           (nodeComponentSpheroid), pointer       :: spheroidHost                  , spheroid
    type            (history              )                :: historyDisk                   , historySpheroid                , &
         &                                                    history_
    double precision                                       :: angularMomentum               , diskSpecificAngularMomentum    , &
         &                                                    spheroidMass                  , spheroidSpecificAngularMomentum, &
         &                                                    radiusRemnant                 ,velocityCircularRemnant         , &
         &                                                    angularMomentumSpecificRemnant
    integer                                                :: destinationGasSatellite       , destinationGasHost             , &
         &                                                    destinationStarsHost          , destinationStarsSatellite
    logical                                                :: mergerIsMajor
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
       select case (destinationGasHost)
       case (destinationMergerDisk)
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
       case (destinationMergerSpheroid)
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
       case (destinationMergerUnmoved)
          ! Do nothing.
       case default
          call Galacticus_Error_Report('unrecognized movesTo descriptor'//{introspection:location})
       end select

       ! Move stellar material within the host if necessary.
       select case (destinationStarsHost)
       case (destinationMergerDisk)
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
          call historyDisk    %increment(historySpheroid     ,autoExtend=.true.)
          call historySpheroid%reset  (                    )
          call diskHost    %    starFormationHistorySet(historyDisk    )
          call spheroidHost%    starFormationHistorySet(historySpheroid)
          call historyDisk    %destroy(recordMemory=.false.)
          call historySpheroid%destroy(recordMemory=.false.)
       case (destinationMergerSpheroid)
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
          call spheroidHost%stellarPropertiesHistorySet(historySpheroid)
          call diskHost    %stellarPropertiesHistorySet( historyDisk   )
          ! Also add star formation histories.
          historyDisk    =diskHost    %starFormationHistory()
          historySpheroid=spheroidHost%starFormationHistory()
          call historySpheroid%increment(historyDisk         ,autoExtend=.true.)
          call historyDisk    %reset  (                    )
          call spheroidHost%starFormationHistorySet    (historySpheroid)
          call diskHost    %starFormationHistorySet    (historyDisk    )
          call historyDisk    %destroy(recordMemory=.false.)
          call historySpheroid%destroy(recordMemory=.false.)
          historyDisk    =diskHost    %starFormationHistory()
          historySpheroid=spheroidHost%starFormationHistory()
       case (destinationMergerUnmoved)
          ! Do nothing.
       case default
          call Galacticus_Error_Report('unrecognized movesTo descriptor'//{introspection:location})
       end select
       ! If the entire host disk/spheroid (gas plus stars) was moved to the spheroid/disk, ensure that the
       ! corresponding angular momentum is precisely zero.
       if (diskHost    %massStellar()+diskHost    %massGas() == 0.0d0 .and. diskHost%angularMomentum() /= 0.0d0) &
            & call diskHost%angularMomentumSet    (0.0d0)
       if (spheroidHost%massStellar()+spheroidHost%massGas() == 0.0d0                                          ) &
            & call spheroidHost%angularMomentumSet(0.0d0)

       ! Get specific angular momentum of the spheroid material.
       spheroidMass=spheroid%massGas()+spheroid%massStellar()
       if (spheroidMass > 0.0d0) then
          spheroidSpecificAngularMomentum=spheroid%angularMomentum()/spheroidMass

          ! Move the gas component of the standard spheroid to the host.
          select case (destinationGasSatellite)
          case (destinationMergerDisk)
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
          case (destinationMergerSpheroid)
             call spheroidHost%        massGasSet(                                 &
                  &                                spheroidHost%        massGas()  &
                  &                               +spheroid    %        massGas()  &
                  &                              )
             call spheroidHost%  abundancesGasSet(                                 &
                  &                                spheroidHost%  abundancesGas()  &
                  &                               +spheroid    %  abundancesGas()  &
                  &                              )
          case default
             call Galacticus_Error_Report('unrecognized movesTo descriptor'//{introspection:location})
          end select
          call spheroid%      massGasSet(0.0d0         )
          call spheroid%abundancesGasSet(zeroAbundances)
          ! Move the stellar component of the standard spheroid to the host.
          select case (destinationStarsSatellite)
          case (destinationMergerDisk)
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
             call history_       %increment              (historySpheroid,autoExtend  =.true. )
             call historySpheroid%reset                  (                                    )
             call diskHost       %starFormationHistorySet(history_                            )
             call spheroid       %starFormationHistorySet(historySpheroid                     )
             call history_       %destroy                (                recordMemory=.false.)
             call historySpheroid%destroy                (                recordMemory=.false.)
          case (destinationMergerSpheroid)
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
             call history_       %increment              (historySpheroid,autoExtend  =.true. )
             call historySpheroid%reset                  (                                    )
             call spheroidHost   %starFormationHistorySet(history_                            )
             call spheroid       %starFormationHistorySet(historySpheroid                     )
             call history_       %destroy                (                recordMemory=.false.)
             call historySpheroid%destroy                (                recordMemory=.false.)
          case default
             call Galacticus_Error_Report('unrecognized movesTo descriptor'//{introspection:location})
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

  !# <radiusSolverPlausibility>
  !#  <unitName>Node_Component_Spheroid_Standard_Radius_Solver_Plausibility</unitName>
  !#  <after>Node_Component_Basic_Standard_Plausibility</after>
  !# </radiusSolverPlausibility>
  subroutine Node_Component_Spheroid_Standard_Radius_Solver_Plausibility(node)
    !% Determines whether the spheroid is physically plausible for radius solving tasks. Require that it have non-zero mass and angular momentum.
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
       ! Determine the plausibility of the current spheroid.
       if      (spheroid%angularMomentum()                    <                           0.0d0) &
            & node%isPhysicallyPlausible=.false.
       if      (spheroid%massStellar    ()+spheroid%massGas() <  -spheroidMassToleranceAbsolute) then
          node%isPhysicallyPlausible=.false.
       else if (spheroid%massStellar    ()+spheroid%massGas() >=                          0.0d0) then
          if      (                                                                              &
               &   spheroid%angularMomentum() < 0.0d0                                            &
               &  ) then
             node%isPhysicallyPlausible=.false.
          else
             angularMomentumScale=(                                           &
                  &                 spheroid%massStellar()                    &
                  &                +spheroid%massGas    ()                    &
                  &               )                                           &
                  &               * darkMatterHaloScale_%virialRadius  (node) &
                  &               * darkMatterHaloScale_%virialVelocity(node)
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
       end if

    end select
    return
  end subroutine Node_Component_Spheroid_Standard_Radius_Solver_Plausibility

  double precision function Node_Component_Spheroid_Standard_Radius_Solve(node)
    !% Return the circular radius of the standard spheroid.
    use :: Galacticus_Nodes, only : nodeComponentSpheroid, treeNode
    implicit none
    type (treeNode             ), intent(inout) :: node
    class(nodeComponentSpheroid), pointer       :: spheroid

    spheroid => node%spheroid()
    Node_Component_Spheroid_Standard_Radius_Solve=spheroid%radius()
    return
  end function Node_Component_Spheroid_Standard_Radius_Solve

  double precision function Node_Component_Spheroid_Standard_Velocity_Solve(node)
    !% Return the circular velocity of the standard spheroid.
    use :: Galacticus_Nodes, only : nodeComponentSpheroid, treeNode
    implicit none
    type (treeNode             ), intent(inout) :: node
    class(nodeComponentSpheroid), pointer       :: spheroid

    spheroid => node%spheroid()
    Node_Component_Spheroid_Standard_Velocity_Solve=spheroid%velocity()
    return
  end function Node_Component_Spheroid_Standard_Velocity_Solve

  subroutine Node_Component_Spheroid_Standard_Radius_Solve_Set(node,radius)
    !% Set the scale radius of the standard spheroid.
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
    !% Set the scale velocity of the standard spheroid.
    use :: Galacticus_Nodes, only : nodeComponentSpheroid, treeNode
    implicit none
    type            (treeNode             ), intent(inout) :: node
    double precision                       , intent(in   ) :: velocity
    class           (nodeComponentSpheroid), pointer       :: spheroid

    spheroid => node%spheroid()
    call spheroid%velocitySet(max(velocity,0.0d0))
    return
  end subroutine Node_Component_Spheroid_Standard_Velocity_Solve_Set

  !# <radiusSolverTask>
  !#  <unitName>Node_Component_Spheroid_Standard_Radius_Solver</unitName>
  !# </radiusSolverTask>
  subroutine Node_Component_Spheroid_Standard_Radius_Solver(node,componentActive,specificAngularMomentumRequired,specificAngularMomentum,Radius_Get,Radius_Set,Velocity_Get&
       &,Velocity_Set)
    !% Interface for the size solver algorithm.
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
         &                                                                                         spheroidMass

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
             spheroidMass=spheroid%massGas()+spheroid%massStellar()
             if (spheroidMass > 0.0d0) then
                specificAngularMomentumMean=angularMomentum/spheroidMass
             else
                specificAngularMomentumMean=0.0d0
             end if
             specificAngularMomentum=spheroidAngularMomentumAtScaleRadius*specificAngularMomentumMean
          end if
          ! Associate the pointers with the appropriate property routines.
          Radius_Get   => Node_Component_Spheroid_Standard_Radius_Solve
          Radius_Set   => Node_Component_Spheroid_Standard_Radius_Solve_Set
          Velocity_Get => Node_Component_Spheroid_Standard_Velocity_Solve
          Velocity_Set => Node_Component_Spheroid_Standard_Velocity_Solve_Set
       else
          call Node_Component_Spheroid_Standard_Radius_Solve_Set  (node,0.0d0)
          call Node_Component_Spheroid_Standard_Velocity_Solve_Set(node,0.0d0)
          componentActive=.false.
       end if
    end select
    return
  end subroutine Node_Component_Spheroid_Standard_Radius_Solver

  subroutine Node_Component_Spheroid_Standard_Initializor(self)
    !% Initializes a standard spheroid component.
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDisk, nodeComponentSpheroidStandard, treeNode
    use :: Histories       , only : history
    implicit none
    type            (nodeComponentSpheroidStandard)          :: self
    type            (treeNode                     ), pointer :: selfNode
    class           (nodeComponentDisk            ), pointer :: disk
    class           (nodeComponentBasic           ), pointer :: basic
    type            (history                      )          :: historyStarFormation        , stellarPropertiesHistory      , &
         &                                                      diskStarFormationHistory
    logical                                                  :: createStarFormationHistory  , createStellarPropertiesHistory
    double precision                                         :: timeBegin

    ! Return if already initialized.
    if (self%isInitialized()) return
    ! Get the associated node.
    selfNode => self%host()
    ! Determine which histories must be created.
    historyStarFormation          =self%starFormationHistory            ()
    createStarFormationHistory    =.not.historyStarFormation    %exists ()
    call                                historyStarFormation    %destroy()
    stellarPropertiesHistory      =self%stellarPropertiesHistory        ()
    createStellarPropertiesHistory=.not.stellarPropertiesHistory%exists ()
    call                                stellarPropertiesHistory%destroy()
    ! Create the stellar properties history.
    if (createStellarPropertiesHistory) then
       call stellarPopulationProperties_%historyCreate(selfNode,stellarPropertiesHistory)
       call self%stellarPropertiesHistorySet(                     stellarPropertiesHistory)
    end if
    ! Create the star formation history.
    if (createStarFormationHistory    ) then
       disk => selfNode%disk()
       diskStarFormationHistory=disk%starFormationHistory()
       if (diskStarFormationHistory%exists()) then
          timeBegin=  diskStarFormationHistory%time(1)
       else
          basic    => selfNode%basic()
          timeBegin=  basic   %time ()
       end if
       call starFormationHistory_%create (selfNode,historyStarFormation,timeBegin)
       call self% starFormationHistorySet(         historyStarFormation          )
    end if
    ! Record that the spheroid has been initialized.
    call self%isInitializedSet(.true.)
    return
  end subroutine Node_Component_Spheroid_Standard_Initializor

  subroutine Node_Component_Spheroid_Standard_Star_Formation_History_Extend(node)
    !% Extend the range of a star formation history in a standard spheroid component for {\normalfont \ttfamily node}.
    use :: Galacticus_Nodes, only : nodeComponentSpheroid, treeNode
    implicit none
    type (treeNode             ), intent(inout), target  :: node
    class(nodeComponentSpheroid)               , pointer :: spheroid
    type (history              )                         :: historyStarFormation

    ! Get the spheroid component.
    spheroid => node%spheroid()
    ! Extend the range as necessary.
    historyStarFormation=spheroid%starFormationHistory()
    call historyStarFormation%extend(times=starFormationHistoryTemplate)
    call spheroid%starFormationHistorySet(historyStarFormation)
    return
  end subroutine Node_Component_Spheroid_Standard_Star_Formation_History_Extend

  subroutine Node_Component_Spheroid_Standard_Stellar_Prprts_History_Extend(node)
    !% Extend the range of a stellar properties history in a standard spheroid component for {\normalfont \ttfamily node}.
    use :: Galacticus_Nodes, only : nodeComponentSpheroid, treeNode
    implicit none
    type (treeNode             ), intent(inout), target  :: node
    class(nodeComponentSpheroid)               , pointer :: spheroid
    type (history              )                         :: stellarPropertiesHistory

    ! Get the spheroid component.
    spheroid => node%spheroid()
    ! Extend the range as necessary.
    stellarPropertiesHistory=spheroid%stellarPropertiesHistory()
    call stellarPropertiesHistory%extend(times=stellarPropertiesHistoryTemplate)
    call spheroid%stellarPropertiesHistorySet(stellarPropertiesHistory)
    return
  end subroutine Node_Component_Spheroid_Standard_Stellar_Prprts_History_Extend

  !# <mergerTreeExtraOutputTask>
  !#  <unitName>Node_Component_Spheroid_Standard_Star_Formation_History_Output</unitName>
  !# </mergerTreeExtraOutputTask>
  subroutine Node_Component_Spheroid_Standard_Star_Formation_History_Output(node,iOutput,treeIndex,nodePassesFilter)
    !% Store the star formation history in the output file.
    use            :: Galacticus_Nodes, only : nodeComponentSpheroid, nodeComponentSpheroidStandard, treeNode, defaultSpheroidComponent
    use            :: Histories       , only : history
    use, intrinsic :: ISO_C_Binding   , only : c_size_t
    use            :: Kind_Numbers    , only : kind_int8
    implicit none
    type   (treeNode             ), intent(inout), pointer :: node
    integer(c_size_t             ), intent(in   )          :: iOutput
    integer(kind=kind_int8       ), intent(in   )          :: treeIndex
    logical                       , intent(in   )          :: nodePassesFilter
    class  (nodeComponentSpheroid)               , pointer :: spheroid
    type   (history              )                         :: historyStarFormation

    ! Check if we are the default method.
    if (.not.defaultSpheroidComponent%standardIsActive()) return
    ! Output the star formation history if a spheroid exists for this component.
    spheroid => node%spheroid()
    select type (spheroid)
    class is (nodeComponentSpheroidStandard)
       historyStarFormation=spheroid%starFormationHistory()
       call starFormationHistory_%output(node,nodePassesFilter,historyStarFormation,iOutput,treeIndex,'spheroid')
       call spheroid%starFormationHistorySet(historyStarFormation)
    end select
    return
  end subroutine Node_Component_Spheroid_Standard_Star_Formation_History_Output

  !# <galacticusStateStoreTask>
  !#  <unitName>Node_Component_Spheroid_Standard_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Node_Component_Spheroid_Standard_State_Store(stateFile,gslStateFile,stateOperatorID)
    !% Write the tablulation state to file.
    use            :: Galacticus_Display                   , only : Galacticus_Display_Message, verbosityInfo
    use, intrinsic :: ISO_C_Binding                        , only : c_size_t                  , c_ptr
    use            :: Node_Component_Spheroid_Standard_Data, only : spheroidMassDistribution
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperatorID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call Galacticus_Display_Message('Storing state for: treeNodeMethodSpheroid -> standard',verbosity=verbosityInfo)
    !# <workaround type="gfortran" PR="92836" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=92836">
    !#  <description>Internal file I/O in gfortran can be non-thread safe.</description>
    !# </workaround>
#ifdef THREADSAFEIO
    !$omp critical(gfortranInternalIO)
#endif
    write (stateFile) spheroidAngularMomentumAtScaleRadius
    write (stateFile) associated(spheroidMassDistribution)
#ifdef THREADSAFEIO
    !$omp end critical(gfortranInternalIO)
#endif
    if (associated(spheroidMassDistribution)) call spheroidMassDistribution%stateStore(stateFile,gslStateFile,stateOperatorID)
    return
  end subroutine Node_Component_Spheroid_Standard_State_Store

  !# <galacticusStateRetrieveTask>
  !#  <unitName>Node_Component_Spheroid_Standard_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Node_Component_Spheroid_Standard_State_Retrieve(stateFile,gslStateFile,stateOperationID)
    !% Retrieve the tabulation state from the file.
    use            :: Galacticus_Display                   , only : Galacticus_Display_Message, verbosityInfo
    use            :: Galacticus_Error                     , only : Galacticus_Error_Report
    use, intrinsic :: ISO_C_Binding                        , only : c_size_t                  , c_ptr
    use            :: Node_Component_Spheroid_Standard_Data, only : spheroidMassDistribution
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile
    logical                          :: wasAllocated

    call Galacticus_Display_Message('Retrieving state for: treeNodeMethodSpheroid -> standard',verbosity=verbosityInfo)
#ifdef THREADSAFEIO
    !$omp critical(gfortranInternalIO)
#endif
    read (stateFile) spheroidAngularMomentumAtScaleRadius
    read (stateFile) wasAllocated
#ifdef THREADSAFEIO
    !$omp end critical(gfortranInternalIO)
#endif
    if (wasAllocated) then
       if (.not.associated(spheroidMassDistribution)) call Galacticus_Error_Report('spheroidMassDistribution was stored, but is now not allocated'//{introspection:location})
       call spheroidMassDistribution%stateRestore(stateFile,gslStateFile,stateOperationID)
    end if
    return
  end subroutine Node_Component_Spheroid_Standard_State_Retrieve

end module Node_Component_Spheroid_Standard
