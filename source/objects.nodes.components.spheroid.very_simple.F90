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

!% Contains a module that implements a very simple spheroid component.

module Node_Component_Spheroid_Very_Simple
  !% Implements a very simple spheroid component.
  use :: Dark_Matter_Halo_Scales         , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO        , only : darkMatterProfileDMOClass
  use :: Satellite_Merging_Mass_Movements, only : mergerMassMovementsClass
  use :: Stellar_Population_Properties   , only : stellarPopulationPropertiesClass
  implicit none
  private
  public :: Node_Component_Spheroid_Very_Simple_Thread_Initialize         , Node_Component_Spheroid_Very_Simple_Post_Step          , &
       &    Node_Component_Spheroid_Very_Simple_Scale_Set                 , Node_Component_Spheroid_Very_Simple_Thread_Uninitialize, &
       &    Node_Component_Spheroid_Very_Simple_Initialize                , Node_Component_Spheroid_Very_Simple_Pre_Evolve         , &
       &    Node_Component_Spheroid_Very_Simple_Radius_Solver_Plausibility, Node_Component_Spheroid_Very_Simple_Radius_Solver

  !# <component>
  !#  <class>spheroid</class>
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
  !#     <type>double</type>
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
  !#     <name>stellarPropertiesHistory</name>
  !#     <type>history</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#   </property>
  !#   <property>
  !#     <name>luminositiesStellar</name>
  !#     <type>stellarLuminosities</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <output unitsInSI="luminosityZeroPointAB" comment="Luminosity of spheroid stars."/>
  !#   </property>
  !#   <property>
  !#     <name>radius</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <output unitsInSI="megaparsec" comment="Radial scale length in the spheroid."/>
  !#   </property>
  !#   <property>
  !#     <name>halfMassRadius</name>
  !#     <attributes isSettable="false" isGettable="true" isEvolvable="false" isVirtual="true" />
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <getFunction>Node_Component_Spheroid_Very_Simple_Half_Mass_Radius</getFunction>
  !#   </property>
  !#   <property>
  !#     <name>velocity</name>
  !#     <type>double</type>
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

  ! Objects used by this component.
  class(darkMatterHaloScaleClass        ), pointer :: darkMatterHaloScale_
  class(darkMatterProfileDMOClass       ), pointer :: darkMatterProfileDMO_
  class(mergerMassMovementsClass        ), pointer :: mergerMassMovements_
  class(stellarPopulationPropertiesClass), pointer :: stellarPopulationProperties_
  !$omp threadprivate(darkMatterHaloScale_,darkMatterProfileDMO_,mergerMassMovements_,stellarPopulationProperties_)

  ! Parameters controlling the physical implementation.
  double precision :: spheroidMassToleranceAbsolute      , spheroidVerySimpleMassScaleAbsolute
  logical          :: spheroidVerySimpleTrackAbundances  , spheroidVerySimpleTrackLuminosities

contains

  !# <nodeComponentInitializationTask>
  !#  <unitName>Node_Component_Spheroid_Very_Simple_Initialize</unitName>
  !# </nodeComponentInitializationTask>
  subroutine Node_Component_Spheroid_Very_Simple_Initialize(parameters_)
    !% Initializes the tree node very simple spheroid component module.
    use :: Galacticus_Nodes, only : defaultSpheroidComponent
    use :: Input_Parameters, only : inputParameter          , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_

    ! Initialize the module if necessary.
    if (defaultSpheroidComponent%verySimpleIsActive()) then
       ! Read parameters controlling the physical implementation.
       !# <inputParameter>
       !#   <name>spheroidVerySimpleMassScaleAbsolute</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultValue>100.0d0</defaultValue>
       !#   <description>The absolute mass scale below which calculations in the very simple spheroid component are allowed to become inaccurate.</description>
       !#   <source>parameters_</source>
       !#   <type>real</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>spheroidVerySimpleTrackAbundances</name>
       !#   <cardinality>0..1</cardinality>
       !#   <defaultValue>.false.</defaultValue>
       !#   <description>Specifies whether or not to track abundances in the very simple spheroid component.</description>
       !#   <source>parameters_</source>
       !#   <type>boolean</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>spheroidVerySimpleTrackLuminosities</name>
       !#   <cardinality>0..1</cardinality>
       !#   <defaultValue>.false.</defaultValue>
       !#   <description>Specifies whether or not to track stellar luminosities in the very simple disk component.</description>
       !#   <source>parameters_</source>
       !#   <type>boolean</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>spheroidMassToleranceAbsolute</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultValue>1.0d-6</defaultValue>
       !#   <description>The mass tolerance used to judge whether the spheroid is physically plausible.</description>
       !#   <source>parameters_</source>
       !#   <type>real</type>
       !# </inputParameter>
    end if
    return
  end subroutine Node_Component_Spheroid_Very_Simple_Initialize

  !# <nodeComponentThreadInitializationTask>
  !#  <unitName>Node_Component_Spheroid_Very_Simple_Thread_Initialize</unitName>
  !# </nodeComponentThreadInitializationTask>
  subroutine Node_Component_Spheroid_Very_Simple_Thread_Initialize(parameters_)
    !% Initializes the tree node very simple satellite module.
    use :: Events_Hooks    , only : satelliteMergerEvent    , postEvolveEvent, openMPThreadBindingAtLevel, dependencyRegEx, &
         &                          dependencyDirectionAfter
    use :: Galacticus_Nodes, only : defaultSpheroidComponent
    use :: Input_Parameters, only : inputParameter          , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_
    type(dependencyRegEx), dimension(1)  :: dependencies

    if (defaultSpheroidComponent%verySimpleIsActive()) then
       dependencies(1)=dependencyRegEx(dependencyDirectionAfter,'^remnantStructure:')
       call postEvolveEvent     %attach(defaultSpheroidComponent,postEvolve     ,openMPThreadBindingAtLevel,label='nodeComponentSpheroidVerySimple'                          )
       call satelliteMergerEvent%attach(defaultSpheroidComponent,satelliteMerger,openMPThreadBindingAtLevel,label='nodeComponentSpheroidVerySimple',dependencies=dependencies)
       !# <objectBuilder class="darkMatterHaloScale"         name="darkMatterHaloScale_"         source="parameters_"/>
       !# <objectBuilder class="darkMatterProfileDMO"        name="darkMatterProfileDMO_"        source="parameters_"/>
       !# <objectBuilder class="stellarPopulationProperties" name="stellarPopulationProperties_" source="parameters_"/>
       !# <objectBuilder class="mergerMassMovements"         name="mergerMassMovements_"         source="parameters_"/>
    end if
    return
  end subroutine Node_Component_Spheroid_Very_Simple_Thread_Initialize

  !# <nodeComponentThreadUninitializationTask>
  !#  <unitName>Node_Component_Spheroid_Very_Simple_Thread_Uninitialize</unitName>
  !# </nodeComponentThreadUninitializationTask>
  subroutine Node_Component_Spheroid_Very_Simple_Thread_Uninitialize()
    !% Uninitializes the tree node very simple satellite module.
    use :: Events_Hooks    , only : satelliteMergerEvent    , postEvolveEvent
    use :: Galacticus_Nodes, only : defaultSpheroidComponent
    implicit none

    if (defaultSpheroidComponent%verySimpleIsActive()) then
       call postEvolveEvent     %detach(defaultSpheroidComponent,postEvolve     )
       call satelliteMergerEvent%detach(defaultSpheroidComponent,satelliteMerger)
       !# <objectDestructor name="darkMatterHaloScale_"         />
       !# <objectDestructor name="darkMatterProfileDMO_"        />
       !# <objectDestructor name="stellarPopulationProperties_" />
       !# <objectDestructor name="mergerMassMovements_"         />
    end if
    return
  end subroutine Node_Component_Spheroid_Very_Simple_Thread_Uninitialize

  !# <preEvolveTask>
  !# <unitName>Node_Component_Spheroid_Very_Simple_Pre_Evolve</unitName>
  !# </preEvolveTask>
  subroutine Node_Component_Spheroid_Very_Simple_Pre_Evolve(node)
    !% Ensure the spheroid has been initialized.
    use :: Galacticus_Nodes, only : nodeComponentSpheroid, nodeComponentSpheroidVerySimple, treeNode, defaultSpheroidComponent
    implicit none
    type (treeNode             ), intent(inout), pointer :: node
    class(nodeComponentSpheroid)               , pointer :: spheroid

    ! Check if we are the default method.
    if (.not.defaultSpheroidComponent%verySimpleIsActive()) return
    ! Get the spheroid component.
    spheroid => node%spheroid()
    ! Check if an exponential spheroid component exists.
    select type (spheroid)
       class is (nodeComponentSpheroidVerySimple)
          ! Initialize the spheroid
          call Node_Component_Spheroid_Very_Simple_Create(node)
    end select
    return
  end subroutine Node_Component_Spheroid_Very_Simple_Pre_Evolve

  subroutine postEvolve(self,node)
    !% Catch rounding errors in the very simple spheroid gas evolution.
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentSpheroid, nodeComponentSpheroidVerySimple, treeNode
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
    ! Check if a very simple spheroid component exists.
    select type (spheroid)
    class is (nodeComponentSpheroidVerySimple)
       ! Trim the stellar populations properties future history.
       basic => node%basic()
       stellarPropertiesHistory=spheroid%stellarPropertiesHistory()
       call stellarPropertiesHistory%trim(basic%time())
       call spheroid%stellarPropertiesHistorySet(stellarPropertiesHistory)
    end select
    return
  end subroutine postEvolve

  !# <postStepTask>
  !# <unitName>Node_Component_Spheroid_Very_Simple_Post_Step</unitName>
  !# </postStepTask>
  subroutine Node_Component_Spheroid_Very_Simple_Post_Step(node,status)
    !% Catch rounding errors in the very simple spheroid gas evolution.
    use :: Abundances_Structure          , only : abs                       , zeroAbundances
    use :: Galacticus_Display            , only : Galacticus_Display_Message, verbosityWarn
    use :: Galacticus_Nodes              , only : nodeComponentSpheroid     , nodeComponentSpheroidVerySimple, treeNode    , defaultSpheroidComponent
    use :: Interface_GSL                 , only : GSL_Failure
    use :: ISO_Varying_String            , only : varying_string            , assignment(=)                  , operator(//)
    use :: Stellar_Luminosities_Structure, only : abs                       , zeroStellarLuminosities
    use :: String_Handling               , only : operator(//)
    implicit none
    type            (treeNode              ), intent(inout), pointer :: node
    integer                                 , intent(inout)          :: status
    class           (nodeComponentSpheroid )               , pointer :: spheroid
    double precision                        , save                   :: fractionalErrorMaximum  =0.0d0
    double precision                                                 :: spheroidMass                  , fractionalError
    character       (len=20                )                         :: valueString
    type            (varying_string        ), save                   :: message
    !$omp threadprivate(message)

    ! Return immediately if this class is not in use.
    if (.not.defaultSpheroidComponent%verySimpleIsActive()) return
    ! Get the spheroid component.
    spheroid => node%spheroid()
    ! Check if a very simple spheroid component exists.
    select type (spheroid)
    class is (nodeComponentSpheroidVerySimple)
       ! Trap negative gas masses.
       if (spheroid%massGas() < 0.0d0) then
          ! Check if this exceeds the maximum previously recorded error.
          fractionalError=   abs(spheroid%massGas    ()) &
               &          /(                             &
               &                 spheroid%massStellar()  &
               &            +abs(spheroid%massGas    ()) &
               &           )
          !$omp critical (Very_Simple_Spheroid_Post_Evolve_Check)
          if (fractionalError > fractionalErrorMaximum) then
             ! Report a warning.
             message='Warning: spheroid has negative gas mass (fractional error exceeds any previously reported):'//char(10)
             message=message//'  Node index        = '//node%index() //char(10)
             write (valueString,'(e12.6)') spheroid%massGas()
             message=message//'  Spheroid gas mass     = '//trim(valueString)//char(10)
             write (valueString,'(e12.6)') spheroid%massStellar()
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
          spheroidMass=+spheroid%massGas    () &
               &       +spheroid%massStellar()
          if (spheroidMass == 0.0d0) then
             call spheroid%        massStellarSet(                  0.0d0)
             call spheroid%  abundancesStellarSet(         zeroAbundances)
             call spheroid%luminositiesStellarSet(zeroStellarLuminosities)
          end if
          ! Reset the gas mass of the spheroid.
          call spheroid%      massGasSet(         0.0d0)
          call spheroid%abundancesGasSet(zeroAbundances)
          status=GSL_Failure
       end if
    end select
    return
  end subroutine Node_Component_Spheroid_Very_Simple_Post_Step

  subroutine Node_Component_Spheroid_Very_Simple_Create(node)
    !% Create properties in a very simple spheroid component.
    use :: Galacticus_Nodes, only : nodeComponentSpheroid, treeNode
    use :: Histories       , only : history
    implicit none
    type   (treeNode             ), intent(inout), target  :: node
    class  (nodeComponentSpheroid)               , pointer :: spheroid
    type   (history              )                         :: stellarPropertiesHistory
    logical                                                :: createStellarPropertiesHistory

    ! Get the spheroid component.
    spheroid => node%spheroid()
    ! Exit if already initialized.
    if (spheroid%isInitialized()) return
    ! Determine which histories must be created.
    stellarPropertiesHistory      =spheroid%stellarPropertiesHistory()
    createStellarPropertiesHistory=.not.stellarPropertiesHistory%exists ()
    call                                stellarPropertiesHistory%destroy()
    ! Create the stellar properties history.
    if (createStellarPropertiesHistory) then
       ! Create the stellar properties history.
       call stellarPopulationProperties_%historyCreate(node,stellarPropertiesHistory)
       call spheroid%stellarPropertiesHistorySet      (     stellarPropertiesHistory)
    end if
    ! Record that the spheroid has been initialized.
    call spheroid%isInitializedSet(.true.)
    return
  end subroutine Node_Component_Spheroid_Very_Simple_Create

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Spheroid_Very_Simple_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Spheroid_Very_Simple_Scale_Set(node)
    !% Set scales for properties of {\normalfont \ttfamily node}.
    use :: Abundances_Structure          , only : abs                  , abundances                     , max                , unitAbundances          , &
          &                                       zeroAbundances
    use :: Galacticus_Nodes              , only : nodeComponentSpheroid, nodeComponentSpheroidVerySimple, treeNode           , defaultSpheroidComponent
    use :: Histories                     , only : history
    use :: Stellar_Luminosities_Structure, only : abs                  , max                            , stellarLuminosities, unitStellarLuminosities
    implicit none
    type            (treeNode             ), intent(inout), pointer :: node
    class           (nodeComponentSpheroid)               , pointer :: spheroid
    double precision                       , parameter              :: luminosityMinimum             =1.0d0
    double precision                                                :: mass
    type            (history              )                         :: stellarPopulationHistoryScales
    type            (abundances           )                         :: abundancesTotal
    type            (stellarLuminosities  )                         :: stellarLuminositiesScale

    ! Check if we are the default method.
    if (.not.defaultSpheroidComponent%verySimpleIsActive()) return
    ! Get the spheroid component.
    spheroid => node%spheroid()
    ! Check if a very simple spheroid component exists.
    select type (spheroid)
    class is (nodeComponentSpheroidVerySimple)
       ! Set scale for gas and stellar mass.
       mass=spheroid%massGas()+spheroid%massStellar()
       call spheroid%massGasScale    (max(mass,spheroidVerySimpleMassScaleAbsolute))
       call spheroid%massStellarScale(max(mass,spheroidVerySimpleMassScaleAbsolute))
       ! Set scale for gas and stellar abundances.
       abundancesTotal=spheroid%abundancesGas()+spheroid%abundancesStellar()
       call spheroid%abundancesGasScale    (max(abundancesTotal,unitAbundances*spheroidVerySimpleMassScaleAbsolute))
       call spheroid%abundancesStellarScale(max(abundancesTotal,unitAbundances*spheroidVerySimpleMassScaleAbsolute))
       ! Set scales for stellar population properties and star formation histories.
       stellarPopulationHistoryScales=spheroid%stellarPropertiesHistory()
       call stellarPopulationProperties_%scales   (spheroid%massStellar(),zeroAbundances,stellarPopulationHistoryScales)
       call spheroid%stellarPropertiesHistoryScale(                                      stellarPopulationHistoryScales)
       call stellarPopulationHistoryScales%destroy()
       ! Set scale for stellar luminosities.
       stellarLuminositiesScale=max(                                      &
            &                       +abs(spheroid%luminositiesStellar()), &
            &                       +unitStellarLuminosities              &
            &                       *luminosityMinimum                    &
            &                      )
       call stellarLuminositiesScale%truncate                (spheroid                %luminositiesStellar())
       call spheroid                %luminositiesStellarScale(stellarLuminositiesScale                      )
    end select
    return
  end subroutine Node_Component_Spheroid_Very_Simple_Scale_Set

  subroutine satelliteMerger(self,node)
    !% Transfer any very simple spheroid associated with {\normalfont \ttfamily node} to its host halo.
    use :: Abundances_Structure                , only : zeroAbundances
    use :: Galacticus_Error                    , only : Galacticus_Error_Report
    use :: Galacticus_Nodes                    , only : nodeComponentDisk       , nodeComponentSpheroid    , nodeComponentSpheroidVerySimple, treeNode
    use :: Histories                           , only : history
    use :: Satellite_Merging_Mass_Movements    , only : destinationMergerDisk   , destinationMergerSpheroid, destinationMergerUnmoved
    use :: Stellar_Luminosities_Structure      , only : zeroStellarLuminosities
    implicit none
    class  (*                    ), intent(inout) :: self
    type   (treeNode             ), intent(inout) :: node
    type   (treeNode             ), pointer       :: nodeHost
    class  (nodeComponentSpheroid), pointer       :: spheroidHost           , spheroid
    class  (nodeComponentDisk    ), pointer       :: diskHost
    type   (history              )                :: historyDisk            , historySpheroid          , &
         &                                           historyHost
    integer                                       :: destinationGasSatellite, destinationGasHost       , &
         &                                           destinationStarsHost   , destinationStarsSatellite
    logical                                       :: mergerIsMajor
    !$GLC attributes unused :: self

    ! Get the spheroid component, creating it if need be.
    spheroid => node%spheroid(autoCreate=.true.)
    select type (spheroid)
    class is (nodeComponentSpheroidVerySimple)
       ! Find the node to merge with and its spheroid component.
       nodeHost     => node    %mergesWith(                 )
       spheroidHost => nodeHost%spheroid  (autoCreate=.true.)
       diskHost     => nodeHost%disk      (autoCreate=.true.)
       ! Get mass movement descriptors.
       call mergerMassMovements_%get(node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
       ! Move gas material within the host if necessary.
       select case (destinationGasHost)
       case (destinationMergerDisk)
          call diskHost    %          massGasSet(                                       &
               &                                          +diskHost    %massGas      () &
               &                                          +spheroidHost%massGas      () &
               &                                         )
          call diskHost    %    abundancesGasSet(                                       &
               &                                          +diskHost    %abundancesGas() &
               &                                          +spheroidHost%abundancesGas() &
               &                                         )
          call spheroidHost%          massGasSet(                                       &
               &                                           0.0d0                        &
               &                                         )
          call spheroidHost%    abundancesGasSet(                                       &
               &                                           zeroAbundances               &
               &                                         )
       case (destinationMergerSpheroid)
          call spheroidHost%          massGasSet(                                       &
               &                                          +spheroidHost%massGas      () &
               &                                          +diskHost    %massGas      () &
               &                                         )
          call spheroidHost%    abundancesGasSet(                                       &
               &                                          +spheroidHost%abundancesGas() &
               &                                          +diskHost    %abundancesGas() &
               &                                         )
          call diskHost    %          massGasSet(                                       &
               &                                           0.0d0                        &
               &                                         )
          call diskHost    %    abundancesGasSet(                                       &
               &                                           zeroAbundances               &
               &                                         )
       case (destinationMergerUnmoved)
          ! Do nothing.
       case default
          call Galacticus_Error_Report(                                    &
               &                       'unrecognized movesTo descriptor'// &
               &                       {introspection:location}            &
               &                      )
       end select
       ! Move stellar material within the host if necessary.
       select case (destinationStarsHost)
       case (destinationMergerDisk)
          call diskHost    %      massStellarSet  (                                    &
               &                                   +diskHost    %        massStellar() &
               &                                   +spheroidHost%        massStellar() &
               &                                  )
          call diskHost    %abundancesStellarSet  (                                    &
               &                                    diskHost    %  abundancesStellar() &
               &                                   +spheroidHost%  abundancesStellar() &
               &                                  )
          call diskHost    %luminositiesStellarSet(                                    &
               &                                    diskHost    %luminositiesStellar() &
               &                                   +spheroidHost%luminositiesStellar() &
               &                                  )
          call spheroidHost%      massStellarSet  (                                    &
               &                                    0.0d0                              &
               &                                  )
          call spheroidHost%abundancesStellarSet  (                                    &
               &                                    zeroAbundances                     &
               &                                  )
          call spheroidHost%luminositiesStellarSet(                                    &
               &                                    zeroStellarLuminosities            &
               &                                   )
          ! Also add stellar properties histories.
          historyDisk    =    diskHost%stellarPropertiesHistory()
          historySpheroid=spheroidHost%stellarPropertiesHistory()
          call historyDisk    %interpolatedIncrement                  (historySpheroid)
          call historySpheroid%reset                      (               )
          call diskHost       %stellarPropertiesHistorySet(historyDisk    )
          call spheroidHost   %stellarPropertiesHistorySet(historySpheroid)
       case (destinationMergerSpheroid)
          call spheroidHost%        massStellarSet(                                    &
               &                                    spheroidHost%        massStellar() &
               &                                   +diskHost    %        massStellar() &
               &                                  )
          call spheroidHost%  abundancesStellarSet(                                    &
               &                                    spheroidHost%  abundancesStellar() &
               &                                   +diskHost    %  abundancesStellar() &
               &                                  )
          call spheroidHost%luminositiesStellarSet(                                    &
               &                                    spheroidHost%luminositiesStellar() &
               &                                   +diskHost    %luminositiesStellar() &
               &                                  )
          call diskHost    %        massStellarSet(                                    &
               &                                    0.0d0                              &
               &                                  )
          call diskHost    %  abundancesStellarSet(                                    &
               &                                    zeroAbundances                     &
               &                                  )
          call diskHost    %luminositiesStellarSet(                                    &
               &                                    zeroStellarLuminosities            &
               &                                  )
          ! Also add stellar properties histories.
          historyDisk    =    diskHost%stellarPropertiesHistory()
          historySpheroid=spheroidHost%stellarPropertiesHistory()
          call historySpheroid%interpolatedIncrement                  (historyDisk    )
          call historyDisk    %reset                      (               )
          call spheroidHost   %stellarPropertiesHistorySet(historySpheroid)
          call diskHost       %stellarPropertiesHistorySet(historyDisk    )
       case (destinationMergerUnmoved)
          ! Do nothing.
       case default
          call Galacticus_Error_Report(                                    &
               &                       'unrecognized movesTo descriptor'// &
               &                       {introspection:location}            &
               &                      )
       end select
       ! Move the gas component of the very simple spheroid to the host.
       select case (destinationGasSatellite)
       case (destinationMergerDisk    )
          call diskHost%massGasSet              (                                       &
               &                                          +diskHost    %      massGas() &
               &                                          +spheroid    %      massGas() &
               &                                         )
          call diskHost%abundancesGasSet        (                                       &
               &                                          +diskHost    %abundancesGas() &
               &                                          +spheroid    %abundancesGas() &
               &                                         )
       case (destinationMergerSpheroid)
          call spheroidHost%massGasSet          (                                       &
               &                                          +spheroidHost%      massGas() &
               &                                          +spheroid    %      massGas() &
               &                                         )
          call spheroidHost%abundancesGasSet    (                                       &
               &                                          +spheroidHost%abundancesGas() &
               &                                          +spheroid    %abundancesGas() &
               &                                         )
       case default
          call Galacticus_Error_Report(                                    &
               &                       'unrecognized movesTo descriptor'// &
               &                       {introspection:location}            &
               &                      )
       end select
       call    spheroid%massGasSet          (                                           &
            &                                                                     0.0d0 &
            &                                            )
       call    spheroid%abundancesGasSet    (                                           &
            &                                                            zeroAbundances &
            &                                            )
       ! Move the stellar component of the very simple spheroid to the host.
       select case (destinationStarsSatellite)
       case (destinationMergerDisk    )
          call diskHost%        massStellarSet(                                &
               &                               +diskHost%        massStellar() &
               &                               +spheroid%        massStellar() &
               &                              )
          call diskHost%  abundancesStellarSet(                                &
               &                               +diskHost%  abundancesStellar() &
               &                               +spheroid%  abundancesStellar() &
               &                              )
          call diskHost%luminositiesStellarSet(                                &
               &                               +diskHost%luminositiesStellar() &
               &                               +spheroid%luminositiesStellar() &
               &                              )
          ! Also add stellar properties histories.
          historySpheroid=spheroid%stellarPropertiesHistory()
          historyHost    =diskHost%stellarPropertiesHistory()
          call historyHost    %interpolatedIncrement                  (historySpheroid)
          call historySpheroid%reset                      (               )
          call diskHost       %stellarPropertiesHistorySet(historyHost    )
          call spheroid       %stellarPropertiesHistorySet(historySpheroid)
       case (destinationMergerSpheroid)
          call spheroidHost%        massStellarSet(                                    &
               &                                   +spheroidHost%        massStellar() &
               &                                   +spheroid    %        massStellar() &
               &                                  )
          call spheroidHost%  abundancesStellarSet(                                    &
               &                                   +spheroidHost%  abundancesStellar() &
               &                                   +spheroid    %  abundancesStellar() &
               &                                  )
          call spheroidHost%luminositiesStellarSet(                                     &
               &                                   +spheroidHost%luminositiesStellar() &
               &                                   +spheroid    %luminositiesStellar() &
               &                                  )
          ! Also add stellar properties histories.
          historySpheroid=spheroid%stellarPropertiesHistory()
          historyHost=spheroidHost%stellarPropertiesHistory()
          call historyHost    %interpolatedIncrement                  (historySpheroid)
          call historySpheroid%reset                      (               )
          call diskHost       %stellarPropertiesHistorySet(historyHost    )
          call spheroid       %stellarPropertiesHistorySet(historySpheroid)
       case default
          call Galacticus_Error_Report(                                    &
               &                       'unrecognized movesTo descriptor'// &
               &                       {introspection:location}            &
               &                      )
       end select
       call    spheroid%         massStellarSet(                       &
            &                                   0.0d0                  &
            &                                  )
       call    spheroid%  abundancesStellarSet(                        &
            &                                  zeroAbundances          &
            &                                 )
       call    spheroid%luminositiesStellarSet(                        &
            &                                  zeroStellarLuminosities &
            &                                 )
    end select
    return
  end subroutine satelliteMerger

  !# <radiusSolverPlausibility>
  !#  <unitName>Node_Component_Spheroid_Very_Simple_Radius_Solver_Plausibility</unitName>
  !# </radiusSolverPlausibility>
  subroutine Node_Component_Spheroid_Very_Simple_Radius_Solver_Plausibility(node)
    !% Determines whether the spheroid is physically plausible for radius solving tasks. Require that it have non-zero mass.
    use :: Galacticus_Nodes, only : defaultSpheroidComponent, nodeComponentSpheroid, nodeComponentSpheroidVerySimple, treeNode
    implicit none
    type   (treeNode             ), intent(inout) :: node
    class  (nodeComponentSpheroid), pointer       :: spheroid

    ! Return immediately if our method is not selected.
    if (.not.defaultSpheroidComponent%verySimpleIsActive()) return
    ! Determine the plausibility of the current spheroid.
    spheroid => node%spheroid()
     select type (spheroid)
     class is (nodeComponentSpheroidVerySimple)
        if (spheroid%massStellar()+spheroid%massGas() < -spheroidMassToleranceAbsolute) node%isPhysicallyPlausible=.false.
     end select
    return
  end subroutine Node_Component_Spheroid_Very_Simple_Radius_Solver_Plausibility

  !# <radiusSolverTask>
  !#  <unitName>Node_Component_Spheroid_Very_Simple_Radius_Solver</unitName>
  !# </radiusSolverTask>
  subroutine Node_Component_Spheroid_Very_Simple_Radius_Solver(node,componentActive,specificAngularMomentumRequired,specificAngularMomentum,Radius_Get&
       &,Radius_Set,Velocity_Get,Velocity_Set)
    !% Interface for the size solver algorithm.
    use :: Dark_Matter_Halo_Spins, only : Dark_Matter_Halo_Angular_Momentum
    use :: Galacticus_Nodes      , only : nodeComponentBasic               , nodeComponentSpheroid, nodeComponentSpheroidVerySimple, treeNode
    implicit none
    type            (treeNode                                      ), intent(inout)          :: node
    logical                                                         , intent(  out)          :: componentActive
    logical                                                         , intent(in   )          :: specificAngularMomentumRequired
    double precision                                                , intent(  out)          :: specificAngularMomentum
    procedure       (Node_Component_Spheroid_Very_Simple_Radius    ), intent(  out), pointer :: Radius_Get                     , Velocity_Get
    procedure       (Node_Component_Spheroid_Very_Simple_Radius_Set), intent(  out), pointer :: Radius_Set                     , Velocity_Set
    class           (nodeComponentSpheroid                         )               , pointer :: spheroid
    class           (nodeComponentBasic                            )               , pointer :: basic

    ! Determine if node has an active spheroid component supported by this module.
    componentActive =  .false.
    spheroid        => node%spheroid()
    select type (spheroid)
    class is (nodeComponentSpheroidVerySimple)
       componentActive        =  .true.
       if (specificAngularMomentumRequired) then
          basic                   => node             %basic()
          specificAngularMomentum =  Dark_Matter_Halo_Angular_Momentum(node,darkMatterProfileDMO_)/basic%mass()
       end if
       ! Associate the pointers with the appropriate property routines.
       Radius_Get   => Node_Component_Spheroid_Very_Simple_Radius
       Radius_Set   => Node_Component_Spheroid_Very_Simple_Radius_Set
       Velocity_Get => Node_Component_Spheroid_Very_Simple_Velocity
       Velocity_Set => Node_Component_Spheroid_Very_Simple_Velocity_Set
    end select
    return
  end subroutine Node_Component_Spheroid_Very_Simple_Radius_Solver

  double precision function Node_Component_Spheroid_Very_Simple_Radius(node)
    !% Return the radius of the spheroid used in structure solvers.
    use :: Galacticus_Nodes, only : nodeComponentSpheroid, treeNode
    implicit none
    type (treeNode             ), intent(inout) :: node
    class(nodeComponentSpheroid), pointer       :: spheroid

    spheroid => node%spheroid()
    Node_Component_Spheroid_Very_Simple_Radius=spheroid%radius()
    return
  end function Node_Component_Spheroid_Very_Simple_Radius

  subroutine Node_Component_Spheroid_Very_Simple_Radius_Set(node,radius)
    !% Set the radius of the spheroid used in structure solvers.
    use :: Galacticus_Nodes, only : nodeComponentSpheroid, treeNode
    implicit none
    type            (treeNode             ), intent(inout) :: node
    double precision                       , intent(in   ) :: radius
    class           (nodeComponentSpheroid), pointer       :: spheroid

    spheroid => node%spheroid()
    call spheroid%radiusSet(radius)
    return
  end subroutine Node_Component_Spheroid_Very_Simple_Radius_Set

  double precision function Node_Component_Spheroid_Very_Simple_Velocity(node)
    !% Return the circular velocity of the spheroid.
    use :: Galacticus_Nodes, only : nodeComponentSpheroid, treeNode
    implicit none
    type (treeNode             ), intent(inout) :: node
    class(nodeComponentSpheroid), pointer       :: spheroid

    spheroid => node%spheroid()
    Node_Component_Spheroid_Very_Simple_Velocity=spheroid%velocity()
    return
  end function Node_Component_Spheroid_Very_Simple_Velocity

  subroutine Node_Component_Spheroid_Very_Simple_Velocity_Set(node,velocity)
    !% Set the circular velocity of the spheroid.
    use :: Galacticus_Nodes, only : nodeComponentSpheroid, treeNode
    implicit none
    type            (treeNode             ), intent(inout) :: node
    double precision                       , intent(in   ) :: velocity
    class           (nodeComponentSpheroid), pointer       :: spheroid

    spheroid => node%spheroid()
    call spheroid%velocitySet(velocity)
    return
  end subroutine Node_Component_Spheroid_Very_Simple_Velocity_Set

end module Node_Component_Spheroid_Very_Simple
