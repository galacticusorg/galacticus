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

!% Contains a module which implements the standard disk node component.

module Node_Component_Disk_Standard
  !% Implements the standard disk node component.
  use :: Dark_Matter_Halo_Scales                    , only : darkMatterHaloScaleClass
  use :: Galactic_Dynamics_Bar_Instabilities        , only : galacticDynamicsBarInstabilityClass
  use :: Ram_Pressure_Stripping_Mass_Loss_Rate_Disks, only : ramPressureStrippingDisksClass
  use :: Satellite_Merging_Mass_Movements           , only : mergerMassMovementsClass
  use :: Star_Formation_Histories                   , only : starFormationHistory                    , starFormationHistoryClass
  use :: Stellar_Population_Properties              , only : stellarPopulationPropertiesClass
  use :: Tidal_Stripping_Mass_Loss_Rate_Disks       , only : tidalStrippingDisksClass
  implicit none
  private
  public :: Node_Component_Disk_Standard_Scale_Set                    , Node_Component_Disk_Standard_Pre_Evolve         , &
       &    Node_Component_Disk_Standard_Radius_Solver_Plausibility   , Node_Component_Disk_Standard_Radius_Solver      , &
       &    Node_Component_Disk_Standard_Star_Formation_History_Output, Node_Component_Disk_Standard_Rate_Compute       , &
       &    Node_Component_Disk_Standard_Initialize                   , Node_Component_Disk_Standard_Calculation_Reset  , &
       &    Node_Component_Disk_Standard_State_Store                  , Node_Component_Disk_Standard_State_Retrieve     , &
       &    Node_Component_Disk_Standard_Thread_Initialize            , Node_Component_Disk_Standard_Inactive           , &
       &    Node_Component_Disk_Standard_Post_Step                    , Node_Component_Disk_Standard_Thread_Uninitialize

  !# <component>
  !#  <class>disk</class>
  !#  <name>standard</name>
  !#  <isDefault>true</isDefault>
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
  !#     <output unitsInSI="massSolar" comment="Mass of stars in the standard disk."/>
  !#   </property>
  !#   <property>
  !#     <name>massStellarFormed</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#   </property>
  !#   <property>
  !#     <name>fractionMassRetained</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
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
  !#     <output unitsInSI="massSolar" comment="Mass of gas in the standard disk."/>
  !#   </property>
  !#   <property>
  !#     <name>abundancesGas</name>
  !#     <type>abundances</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" makeGeneric="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of metals in the gas phase of the standard disk."/>
  !#   </property>
  !#   <property>
  !#     <name>angularMomentum</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" makeGeneric="true" />
  !#     <output unitsInSI="massSolar*megaParsec*kilo" comment="Angular momentum of the standard disk."/>
  !#   </property>
  !#   <property>
  !#     <name>radius</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <output unitsInSI="megaparsec" comment="Radial scale length in the standard disk."/>
  !#   </property>
  !#   <property>
  !#     <name>halfMassRadius</name>
  !#     <attributes isSettable="false" isGettable="true" isEvolvable="false" isVirtual="true" />
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <getFunction>Node_Component_Disk_Standard_Half_Mass_Radius</getFunction>
  !#   </property>
  !#   <property>
  !#     <name>velocity</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <output unitsInSI="kilo" comment="Circular velocity of the standard disk at scale length."/>
  !#   </property>
  !#   <property>
  !#     <name>luminositiesStellar</name>
  !#     <type>stellarLuminosities</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <output unitsInSI="luminosityZeroPointAB" comment="Luminosity of disk stars."/>
  !#   </property>
  !#   <property>
  !#     <name>stellarPropertiesHistory</name>
  !#     <type>history</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#   </property>
  !#   <property>
  !#     <name>starFormationHistory</name>
  !#     <type>history</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#   </property>
  !#  </properties>
  !#  <bindings>
  !#   <binding method="attachPipes" function="Node_Component_Disk_Standard_Attach_Pipes" description="Attach pipes to the standard disk component." returnType="\void" arguments="" bindsTo="component" />
  !#   <binding method="enclosedMass"          function="Node_Component_Disk_Standard_Enclosed_Mass"           bindsTo="component" />
  !#   <binding method="acceleration"          function="Node_Component_Disk_Standard_Acceleration"            bindsTo="component" />
  !#   <binding method="tidalTensor"           function="Node_Component_Disk_Standard_Tidal_Tensor"            bindsTo="component" />
  !#   <binding method="chandrasekharIntegral" function="Node_Component_Disk_Standard_Chandrasekhar_Integral"  bindsTo="component" />
  !#   <binding method="density"               function="Node_Component_Disk_Standard_Density"                 bindsTo="component" />
  !#   <binding method="potential"             function="Node_Component_Disk_Standard_Potential"               bindsTo="component" />
  !#   <binding method="rotationCurve"         function="Node_Component_Disk_Standard_Rotation_Curve"          bindsTo="component" />
  !#   <binding method="rotationCurveGradient" function="Node_Component_Disk_Standard_Rotation_Curve_Gradient" bindsTo="component" />
  !#   <binding method="surfaceDensity"        function="Node_Component_Disk_Standard_Surface_Density"         bindsTo="component" />
  !#  </bindings>
  !#  <functions>objects.nodes.components.disk.standard.bound_functions.inc</functions>
  !# </component>

  ! Objects used by this component.
  class(darkMatterHaloScaleClass           ), pointer :: darkMatterHaloScale_
  class(stellarPopulationPropertiesClass   ), pointer :: stellarPopulationProperties_
  class(galacticDynamicsBarInstabilityClass), pointer :: galacticDynamicsBarInstability_
  class(ramPressureStrippingDisksClass     ), pointer :: ramPressureStrippingDisks_
  class(tidalStrippingDisksClass           ), pointer :: tidalStrippingDisks_
  class(starFormationHistoryClass          ), pointer :: starFormationHistory_
  class(mergerMassMovementsClass           ), pointer :: mergerMassMovements_
  !$omp threadprivate(darkMatterHaloScale_,stellarPopulationProperties_,galacticDynamicsBarInstability_,ramPressureStrippingDisks_,tidalStrippingDisks_,starFormationHistory_,mergerMassMovements_)

  ! Internal count of abundances.
  integer                                     :: abundancesCount

  ! Parameters controlling the physical implementation.
  double precision                            :: diskMassToleranceAbsolute                         , diskStructureSolverRadius
  logical                                     :: diskNegativeAngularMomentumAllowed                , diskRadiusSolverCole2000Method       , &
       &                                         diskLuminositiesStellarInactive

  ! History of trial radii used to check for oscillations in the solution when solving for the structure of the disk.
  integer                                     :: radiusSolverIteration
  double precision, dimension(2)              :: radiusHistory
  !$omp threadprivate(radiusHistory,radiusSolverIteration)
  ! The largest and smallest angular momentum, in units of that of a circular orbit at the virial radius, considered to be physically plausible for a disk.
  double precision, parameter                 :: angularMomentumMaximum                    =1.0d+1
  double precision, parameter                 :: angularMomentumMinimum                    =1.0d-6

  ! Disk structural parameters.
  double precision                            :: diskStructureSolverSpecificAngularMomentum        , diskRadiusSolverFlatVsSphericalFactor
  !$omp threadprivate(diskStructureSolverSpecificAngularMomentum,diskRadiusSolverFlatVsSphericalFactor)

  ! Pipe attachment status.
  logical                                     :: pipesAttached                             =.false.

contains

  !# <nodeComponentInitializationTask>
  !#  <unitName>Node_Component_Disk_Standard_Initialize</unitName>
  !# </nodeComponentInitializationTask>
  subroutine Node_Component_Disk_Standard_Initialize(parameters_)
    !% Initializes the tree node standard disk methods module.
    use :: Abundances_Structure, only : Abundances_Property_Count
    use :: Galacticus_Error    , only : Galacticus_Error_Report
    use :: Galacticus_Nodes    , only : defaultDiskComponent     , nodeComponentDiskStandard
    use :: Input_Parameters    , only : inputParameter           , inputParameters
    implicit none
    type(inputParameters          ), intent(inout) :: parameters_
    type(nodeComponentDiskStandard)                :: diskStandardComponent

    if (defaultDiskComponent%standardIsActive()) then
       ! Get number of abundance properties.
       abundancesCount  =Abundances_Property_Count            ()
       ! Attach the cooling mass/angular momentum pipes from the hot halo component.
       if (.not.pipesAttached) then
          call diskStandardComponent%attachPipes()
          pipesAttached=.true.
       end if
       ! Read parameters controlling the physical implementation.
       !# <inputParameter>
       !#   <name>diskMassToleranceAbsolute</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultValue>1.0d-6</defaultValue>
       !#   <description>The mass tolerance used to judge whether the disk is physically plausible.</description>
       !#   <source>parameters_</source>
       !#   <type>double</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>diskStructureSolverRadius</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultValue>1.0d0</defaultValue>
       !#   <description>The radius (in units of the standard scale length) to use in solving for the size of the disk.</description>
       !#   <source>parameters_</source>
       !#   <type>double</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>diskRadiusSolverCole2000Method</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultValue>.false.</defaultValue>
       !#   <description></description>
       !#   <source>parameters_</source>
       !#   <type>boolean</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>diskNegativeAngularMomentumAllowed</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultValue>.true.</defaultValue>
       !#   <description>Specifies whether or not negative angular momentum is allowed for the disk.</description>
       !#   <source>parameters_</source>
       !#   <type>double</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>diskLuminositiesStellarInactive</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultValue>.false.</defaultValue>
       !#   <description>Specifies whether or not disk stellar luminosities are inactive properties (i.e. do not appear in any ODE being solved).</description>
       !#   <source>parameters_</source>
       !#   <type>boolean</type>
       !# </inputParameter>
    end if
    return
  end subroutine Node_Component_Disk_Standard_Initialize

  !# <nodeComponentThreadInitializationTask>
  !#  <unitName>Node_Component_Disk_Standard_Thread_Initialize</unitName>
  !# </nodeComponentThreadInitializationTask>
  subroutine Node_Component_Disk_Standard_Thread_Initialize(parameters_)
    !% Initializes the standard disk component module for each thread.
    use :: Events_Hooks                     , only : satelliteMergerEvent       , postEvolveEvent, openMPThreadBindingAtLevel, dependencyRegEx, &
         &                                           dependencyDirectionAfter
    use :: Galacticus_Error                 , only : Galacticus_Error_Report
    use :: Galacticus_Nodes                 , only : defaultDiskComponent
    use :: Input_Parameters                 , only : inputParameter             , inputParameters
    use :: Mass_Distributions               , only : massDistributionCylindrical
    use :: Node_Component_Disk_Standard_Data, only : diskMassDistribution
    implicit none
    type            (inputParameters), intent(inout) :: parameters_
    type            (dependencyRegEx), dimension(1)  :: dependencies
    double precision                                 :: diskMassDistributionDensityMoment1, diskMassDistributionDensityMoment2
    logical                                          :: surfaceDensityMoment1IsInfinite   , surfaceDensityMoment2IsInfinite

    ! Check if this implementation is selected. If so, initialize the mass distribution.
    if (defaultDiskComponent%standardIsActive()) then
       dependencies(1)=dependencyRegEx(dependencyDirectionAfter,'^remnantStructure:')
       call satelliteMergerEvent%attach(defaultDiskComponent,satelliteMerger,openMPThreadBindingAtLevel,label='nodeComponentDiskStandard',dependencies=dependencies)
       call postEvolveEvent     %attach(defaultDiskComponent,postEvolve     ,openMPThreadBindingAtLevel,label='nodeComponentDiskStandard'                          )
       !# <objectBuilder class="darkMatterHaloScale"                                                 name="darkMatterHaloScale_"            source="parameters_"                    />
       !# <objectBuilder class="stellarPopulationProperties"                                         name="stellarPopulationProperties_"    source="parameters_"                    />
       !# <objectBuilder class="galacticDynamicsBarInstability"                                      name="galacticDynamicsBarInstability_" source="parameters_"                    />
       !# <objectBuilder class="ramPressureStrippingDisks"                                           name="ramPressureStrippingDisks_"      source="parameters_"                    />
       !# <objectBuilder class="tidalStrippingDisks"                                                 name="tidalStrippingDisks_"            source="parameters_"                    />
       !# <objectBuilder class="starFormationHistory"                                                name="starFormationHistory_"           source="parameters_"                    />
       !# <objectBuilder class="mergerMassMovements"                                                 name="mergerMassMovements_"            source="parameters_"                    />
       !# <objectBuilder class="massDistribution"               parameterName="diskMassDistribution" name="diskMassDistribution"            source="parameters_" threadPrivate="yes" >
       !#  <default>
       !#   <diskMassDistribution value="exponentialDisk">
       !#    <dimensionless value="true"/>
       !#   </diskMassDistribution>
       !#  </default>
       !# </objectBuilder>
       if (.not.diskMassDistribution%isDimensionless()) call Galacticus_Error_Report('disk mass distribution must be dimensionless'//{introspection:location})
       ! Compute the specific angular momentum of the disk at this structure solver radius in units of the mean specific angular
       ! momentum of the disk assuming a flat rotation curve.
       select type (diskMassDistribution)
       class is (massDistributionCylindrical)
          ! Determine the specific angular momentum at the size solver radius in units of the mean specific angular
          ! momentum of the disk. This is equal to the ratio of the 1st to 2nd radial moments of the surface density
          ! distribution (assuming a flat rotation curve).
          diskMassDistributionDensityMoment1=diskMassDistribution%surfaceDensityRadialMoment(1.0d0,isInfinite=surfaceDensityMoment1IsInfinite)
          diskMassDistributionDensityMoment2=diskMassDistribution%surfaceDensityRadialMoment(2.0d0,isInfinite=surfaceDensityMoment2IsInfinite)
          if (surfaceDensityMoment1IsInfinite.or.surfaceDensityMoment2IsInfinite) then
             ! One or both of the moments are infinite. Simply assume a value of 0.5 as a default.
             diskStructureSolverSpecificAngularMomentum=0.5d0
          else
             diskStructureSolverSpecificAngularMomentum=  &
                  & +diskStructureSolverRadius            &
                  & /(                                    &
                  &   +diskMassDistributionDensityMoment2 &
                  &   /diskMassDistributionDensityMoment1 &
                  &  )
          end if
       class default
          call Galacticus_Error_Report('only cylindrically symmetric mass distributions are allowed'//{introspection:location})
       end select
       ! If necessary, compute the specific angular momentum correction factor to account for the difference between rotation
       ! curves for thin disk and a spherical mass distribution.
       if (diskRadiusSolverCole2000Method) then
          select type (diskMassDistribution)
             class is (massDistributionCylindrical)
             diskRadiusSolverFlatVsSphericalFactor=                                          &
                  & +diskMassDistribution%rotationCurve       (diskStructureSolverRadius)**2 &
                  & *                                          diskStructureSolverRadius     &
                  & -diskMassDistribution%massEnclosedBySphere(diskStructureSolverRadius)
          end select
       end if
    end if
    return
  end subroutine Node_Component_Disk_Standard_Thread_Initialize

  !# <nodeComponentThreadUninitializationTask>
  !#  <unitName>Node_Component_Disk_Standard_Thread_Uninitialize</unitName>
  !# </nodeComponentThreadUninitializationTask>
  subroutine Node_Component_Disk_Standard_Thread_Uninitialize()
    !% Uninitializes the standard disk component module for each thread.
    use :: Events_Hooks                     , only : satelliteMergerEvent, postEvolveEvent
    use :: Galacticus_Nodes                 , only : defaultDiskComponent
    use :: Node_Component_Disk_Standard_Data, only : diskMassDistribution
    implicit none

    if (defaultDiskComponent%standardIsActive()) then
       call satelliteMergerEvent%detach(defaultDiskComponent,satelliteMerger)
       call postEvolveEvent     %detach(defaultDiskComponent,postEvolve     )
       !# <objectDestructor name="darkMatterHaloScale_"           />
       !# <objectDestructor name="stellarPopulationProperties_"   />
       !# <objectDestructor name="galacticDynamicsBarInstability_"/>
       !# <objectDestructor name="ramPressureStrippingDisks_"     />
       !# <objectDestructor name="tidalStrippingDisks_"           />
       !# <objectDestructor name="starFormationHistory_"          />
       !# <objectDestructor name="mergerMassMovements_"           />
       !# <objectDestructor name="diskMassDistribution"           />
    end if
    return
  end subroutine Node_Component_Disk_Standard_Thread_Uninitialize

  !# <calculationResetTask>
  !#   <unitName>Node_Component_Disk_Standard_Calculation_Reset</unitName>
  !# </calculationResetTask>
  subroutine Node_Component_Disk_Standard_Calculation_Reset(node)
    !% Reset standard disk structure calculations.
    use :: Galacticus_Nodes                 , only : treeNode
    use :: Node_Component_Disk_Standard_Data, only : Node_Component_Disk_Standard_Reset
    implicit none
    type(treeNode), intent(inout) :: node

    call Node_Component_Disk_Standard_Reset(node%uniqueID())
    return
  end subroutine Node_Component_Disk_Standard_Calculation_Reset

  !# <preEvolveTask>
  !# <unitName>Node_Component_Disk_Standard_Pre_Evolve</unitName>
  !# </preEvolveTask>
  subroutine Node_Component_Disk_Standard_Pre_Evolve(node)
    !% Ensure the disk has been initialized.
    use :: Galacticus_Nodes, only : nodeComponentDisk, nodeComponentDiskStandard, treeNode, defaultDiskComponent
    implicit none
    type (treeNode         ), intent(inout), pointer :: node
    class(nodeComponentDisk)               , pointer :: disk

    ! Check if we are the default method.
    if (.not.defaultDiskComponent%standardIsActive()) return
    ! Get the disk component.
    disk => node%disk()
    ! Check if an standard disk component exists.
    select type (disk)
    class is (nodeComponentDiskStandard)
       ! Initialize the disk
       call Node_Component_Disk_Standard_Create(node)
    end select
    return
  end subroutine Node_Component_Disk_Standard_Pre_Evolve

  subroutine postEvolve(self,node)
    !% Trim histories attached to the disk.
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDisk, nodeComponentDiskStandard, treeNode
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
    ! Check if an standard disk component exists.
    select type (disk)
    class is (nodeComponentDiskStandard)
       ! Trim the stellar populations properties future history.
       basic => node%basic()
       stellarPropertiesHistory=disk%stellarPropertiesHistory()
       call stellarPropertiesHistory%trim(basic%time())
       call disk%stellarPropertiesHistorySet(stellarPropertiesHistory)
    end select
    return
  end subroutine postEvolve

  !# <postStepTask>
  !# <unitName>Node_Component_Disk_Standard_Post_Step</unitName>
  !# <after>Node_Component_Basic_Standard_Post_Step</after>
  !# </postStepTask>
  subroutine Node_Component_Disk_Standard_Post_Step(node,status)
    !% Trim histories attached to the disk.
    use :: Abundances_Structure          , only : abs                       , zeroAbundances
    use :: Galacticus_Display            , only : Galacticus_Display_Message, verbosityWarn
    use :: Galacticus_Error              , only : Galacticus_Error_Report
    use :: Galacticus_Nodes              , only : nodeComponentDisk         , nodeComponentDiskStandard, nodeComponentSpin, treeNode, &
         &                                        defaultDiskComponent
    use :: Interface_GSL                 , only : GSL_Failure
    use :: ISO_Varying_String            , only : assignment(=)             , operator(//)             , varying_string
    use :: Stellar_Luminosities_Structure, only : stellarLuminosities       , abs
    use :: String_Handling               , only : operator(//)
    implicit none
    type            (treeNode                ), intent(inout), pointer :: node
    integer                                   , intent(inout)          :: status
    class           (nodeComponentDisk       )               , pointer :: disk
    class           (nodeComponentSpin       )               , pointer :: spin
    double precision                          , parameter              :: angularMomentumTolerance=1.0d-2
    double precision                          , save                   :: fractionalErrorMaximum  =0.0d+0
    double precision                                                   :: diskMass                       , fractionalError, &
         &                                                                specificAngularMomentum
    character       (len=20                  )                         :: valueString
    type            (varying_string          ), save                   :: message
    !$omp threadprivate(message)
    type            (stellarLuminosities     ), save                   :: luminositiesStellar
    !$omp threadprivate(luminositiesStellar)

    ! Return immediately if this class is not in use.
    if (.not.defaultDiskComponent%standardIsActive()) return
    ! Get the disk component.
    disk => node%disk()
    ! Check if an standard disk component exists.
    select type (disk)
    class is (nodeComponentDiskStandard)
       ! Trap negative gas masses.
       if (disk%massGas() < 0.0d0) then
          ! Check if this exceeds the maximum previously recorded error.
          fractionalError=   abs(disk%massGas    ()) &
               &          /(                         &
               &             abs(disk%massGas    ()) &
               &            +abs(disk%massStellar()) &
               &           )
          !$omp critical (Standard_Disk_Post_Evolve_Check)
          if (fractionalError > fractionalErrorMaximum) then
             ! Report a warning.
             message='Warning: disk has negative gas mass (fractional error exceeds any previously reported):'//char(10)
             message=message//'  Node index        = '//node%index() //char(10)
             write (valueString,'(e12.6)') disk%massGas    ()
             message=message//'  Disk gas mass     = '//trim(valueString)//char(10)
             write (valueString,'(e12.6)') disk%massStellar()
             message=message//'  Disk stellar mass = '//trim(valueString)//char(10)
             write (valueString,'(e12.6)') fractionalError
             message=message//'  Error measure     = '//trim(valueString)//char(10)
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
          !$omp end critical (Standard_Disk_Post_Evolve_Check)
          ! Get the specific angular momentum of the disk material
          diskMass= disk%massGas    () &
               &   +disk%massStellar()
          if (diskMass == 0.0d0) then
             specificAngularMomentum=0.0d0
             call disk%        massStellarSet(                  0.0d0)
             call disk%  abundancesStellarSet(         zeroAbundances)
             ! We need to reset the stellar luminosities to zero. We can't simply use the "zeroStellarLuminosities" instance since
             ! our luminosities may have been truncated. If we were to use "zeroStellarLuminosities" then the number of stellar
             ! luminosities associated with the disk would change - but we are in the middle of differential evolution here and we
             ! cannot change the number of evolvable properties as doing so will lead to invalid memory accesses during
             ! deserialization of properties from the ODE solver.
             call luminositiesStellar%destroy()
             luminositiesStellar=disk%luminositiesStellar()
             call luminositiesStellar%reset()
             call disk%luminositiesStellarSet(luminositiesStellar)
          else
             specificAngularMomentum=disk%angularMomentum()/diskMass
             if (specificAngularMomentum < 0.0d0) specificAngularMomentum=disk%radius()*disk%velocity()
          end if
          ! Reset the gas, abundances and angular momentum of the disk.
          call disk%        massGasSet(                                     0.0d0)
          call disk%  abundancesGasSet(                            zeroAbundances)
          call disk%angularMomentumSet(specificAngularMomentum*disk%massStellar())
          status=GSL_Failure
       end if
       ! Trap negative stellar masses.
       if (disk%massStellar() < 0.0d0) then
          ! Check if this exceeds the maximum previously recorded error.
          fractionalError=   abs(disk%massStellar()) &
               &          /(                         &
               &             abs(disk%massGas    ()) &
               &            +abs(disk%massStellar()) &
               &           )
          !$omp critical (Standard_Disk_Post_Evolve_Check)
          if (fractionalError > fractionalErrorMaximum) then
             ! Report a warning.
             message='Warning: disk has negative stellar mass (fractional error exceeds any previously reported):'//char(10)
             message=message//'  Node index        = '//node%index() //char(10)
             write (valueString,'(e12.6)') disk%massGas    ()
             message=message//'  Disk gas mass     = '//trim(valueString)//char(10)
             write (valueString,'(e12.6)') disk%massStellar()
             message=message//'  Disk stellar mass = '//trim(valueString)//char(10)
             write (valueString,'(e12.6)') fractionalError
             message=message//'  Error measure     = '//trim(valueString)//char(10)
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
          !$omp end critical (Standard_Disk_Post_Evolve_Check)
          ! Get the specific angular momentum of the disk material
          diskMass= disk%massGas    () &
               &   +disk%massStellar()
          if (diskMass == 0.0d0) then
             specificAngularMomentum=0.0d0
             call disk%      massGasSet(         0.0d0)
             call disk%abundancesGasSet(zeroAbundances)
          else
             specificAngularMomentum=disk%angularMomentum()/diskMass
             if (specificAngularMomentum < 0.0d0) specificAngularMomentum=disk%radius()*disk%velocity()
          end if
          ! Reset the stellar, abundances and angular momentum of the disk.
          call disk%      massStellarSet(                                 0.0d0)
          call disk%abundancesStellarSet(                        zeroAbundances)
          call disk%  angularMomentumSet(specificAngularMomentum*disk%massGas())
          status=GSL_Failure
       end if
       ! Trap negative angular momentum.
       if (disk%angularMomentum() < 0.0d0) then
          spin => node%spin()
          if      (                           &
               &     disk  %massStellar    () &
               &    +disk  %massGas        () &
               &   <=                         &
               &     0.0d0                    &
               &  ) then
             call disk%angularMomentumSet(0.0d0)
          else if (.not.diskNegativeAngularMomentumAllowed) then
             if  (                                 &
                  &    abs(disk%angularMomentum()) &
                  &   /(                           &
                  &        disk%massStellar    ()  &
                  &     +  disk%massGas        ()  &
                  &    )                           &
                  &  <                             &
                  &    angularMomentumTolerance    &
                  &   *darkMatterHaloScale_%virialRadius  (node) &
                  &   *darkMatterHaloScale_%virialVelocity(node) &
                  &   *spin            %spin          (        ) &
                  & ) then
                call disk%angularMomentumSet(0.0d0)
             else
                message='negative angular momentum in disk with positive mass'
                write (valueString,'(e12.6)') disk%angularMomentum()
                message=message//char(10)//' -> angular momentum       = '//trim(valueString)
                write (valueString,'(e12.6)') disk%massStellar    ()
                message=message//char(10)//' -> stellar mass           = '//trim(valueString)
                write (valueString,'(e12.6)') disk%massGas        ()
                message=message//char(10)//' -> gas mass               = '//trim(valueString)
                write (valueString,'(e12.6)') +darkMatterHaloScale_%virialRadius  (node) &
                     &                        *darkMatterHaloScale_%virialVelocity(node) &
                     &                        *spin                %spin          (    )
                message=message//char(10)//' -> angular momentum scale = '//trim(valueString)
                call Galacticus_Error_Report(message//{introspection:location})
             end if
          end if
          status=GSL_Failure
       end if
    end select
    return
  end subroutine Node_Component_Disk_Standard_Post_Step

  subroutine Node_Component_Disk_Standard_Create(node)
    !% Create properties in an standard disk component.
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDisk, nodeComponentSpheroid, treeNode
    use :: Histories       , only : history
    implicit none
    type   (treeNode             ), intent(inout), target  :: node
    class  (nodeComponentDisk    )               , pointer :: disk
    class  (nodeComponentSpheroid)               , pointer :: spheroid
    class  (nodeComponentBasic   )               , pointer :: basic
    type   (history              )                         :: historyStarFormation        , stellarPropertiesHistory      , &
         &                                                    spheroidStarFormationHistory
    logical                                                :: createStarFormationHistory  , createStellarPropertiesHistory
    double precision                                       :: timeBegin

    ! Get the disk component.
    disk => node%disk()
    ! Exit if already initialized.
    if (disk%isInitialized()) return
    ! Determine which histories must be created.
    historyStarFormation          =disk%starFormationHistory            ()
    createStarFormationHistory    =.not.             historyStarFormation    %exists ()
    call                                             historyStarFormation    %destroy()
    stellarPropertiesHistory      =disk%stellarPropertiesHistory        ()
    createStellarPropertiesHistory=.not.             stellarPropertiesHistory%exists ()
    call                                             stellarPropertiesHistory%destroy()
    ! Set the fraction of mass retained.
    call disk%fractionMassRetainedSet(1.0d0)
    ! Create the stellar properties history.
    if (createStellarPropertiesHistory) then
       ! Create the stellar properties history.
       call stellarPopulationProperties_%historyCreate(node,stellarPropertiesHistory)
       call disk%stellarPropertiesHistorySet(stellarPropertiesHistory)
    end if
    ! Create the star formation history.
    if (createStarFormationHistory    ) then
       spheroid => node%spheroid()
       spheroidStarFormationHistory=spheroid%starFormationHistory()
       if (spheroidStarFormationHistory%exists()) then
          timeBegin=  spheroidStarFormationHistory%time(1)
       else
          basic    => node%basic()
          timeBegin=  basic   %time ()
       end if
       call starFormationHistory_%create(node,historyStarFormation,timeBegin)
       call disk%starFormationHistorySet(     historyStarFormation          )
    end if
    ! Record that the disk has been initialized.
    call disk%isInitializedSet(.true.)
    return
  end subroutine Node_Component_Disk_Standard_Create

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Disk_Standard_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Disk_Standard_Rate_Compute(node,odeConverged,interrupt,interruptProcedureReturn,propertyType)
    !% Compute the standard disk node mass rate of change.
    use :: Abundances_Structure            , only : abs                    , abundances           , max               , operator(*)              , &
          &                                         zeroAbundances
    use :: Galacticus_Error                , only : Galacticus_Error_Report
    use :: Galacticus_Nodes                , only : defaultDiskComponent   , interruptTask        , nodeComponentDisk , nodeComponentDiskStandard, &
          &                                         nodeComponentHotHalo   , nodeComponentSpheroid, propertyTypeActive, propertyTypeAll          , &
          &                                         propertyTypeInactive   , treeNode
    use :: Histories                       , only : history                , operator(*)
    use :: Numerical_Constants_Astronomical, only : Mpc_per_km_per_s_To_Gyr
    use :: Stellar_Luminosities_Structure  , only : abs                    , max                  , operator(*)       , stellarLuminosities      , &
         &                                          zeroStellarLuminosities
    implicit none
    type            (treeNode             ), intent(inout), pointer :: node
    logical                                , intent(in   )          :: odeConverged
    class           (nodeComponentDisk    )               , pointer :: disk
    class           (nodeComponentSpheroid)               , pointer :: spheroid
    class           (nodeComponentHotHalo )               , pointer :: hotHalo
    logical                                , intent(inout)          :: interrupt
    procedure       (interruptTask        ), intent(inout), pointer :: interruptProcedureReturn
    integer                                , intent(in   )          :: propertyType
    procedure       (interruptTask        )               , pointer :: interruptProcedure
    type            (abundances           ), save                   :: fuelAbundancesRates         , stellarAbundancesRates
    !$omp threadprivate(stellarAbundancesRates,fuelAbundancesRates)
    double precision                                                :: barInstabilitySpecificTorque, barInstabilityTimescale, &
         &                                                             fractionGas                 , fractionStellar        , &
         &                                                             massLossRate                , transferRate
    type            (history              )                         :: historyTransferRate
    type            (stellarLuminosities  ), save                   :: luminositiesTransferRate
    !$omp threadprivate(luminositiesTransferRate)
    !$GLC attributes unused :: odeConverged

    ! Return immediately if this class is not in use or only inactive properties are being computed.
    if (.not.defaultDiskComponent%standardIsActive() .or. propertyType == propertyTypeInactive) return
    ! Get a local copy of the interrupt procedure.
    interruptProcedure => interruptProcedureReturn
    ! Get the disk and check that it is of our class.
    disk => node%disk()
    select type (disk)
    class is (nodeComponentDiskStandard)
       ! Check for a realistic disk, return immediately if disk is unphysical.
       if     (     disk%angularMomentum() < 0.0d0 &
            &  .or. disk%radius         () < 0.0d0 &
            &  .or. disk%massGas        () < 0.0d0 &
            & ) return
       ! Interrupt if the disk is not initialized.
       if (.not.disk%isInitialized()) then
          interrupt=.true.
          interruptProcedureReturn => Node_Component_Disk_Standard_Create
          return
       end if
       ! Determine if the disk is bar unstable and, if so, the rate at which material is moved to the pseudo-bulge.
       if (node%isPhysicallyPlausible) then
          ! Disk has positive angular momentum, so compute an instability timescale.
          call galacticDynamicsBarInstability_%timescale(node,barInstabilityTimescale,barInstabilitySpecificTorque)
       else
          ! Disk has non-positive angular momentum, therefore it is unphysical. Do not compute an instability timescale in this
          ! case as the disk radius may be unphysical also.
          barInstabilityTimescale=-1.0d0
       end if
       ! Negative timescale indicates no bar instability.
       if (barInstabilityTimescale >= 0.0d0) then
          ! Disk is unstable, so compute rates at which material is transferred to the spheroid.
          spheroid => node%spheroid()
          ! Gas mass.
          transferRate               =max(         0.0d0         ,disk    %massGas             (                         ))/barInstabilityTimescale
          call                                      disk    %massGasRate             (-           transferRate                              )
          call                                      spheroid%massGasRate             (+           transferRate ,interrupt,interruptProcedure)
          ! Fraction of stellar mass transferred.
          transferRate               =max(         0.0d0         ,disk    %fractionMassRetained(                         ))/barInstabilityTimescale
          call                                      disk    %fractionMassRetainedRate(-           transferRate                              )
          ! Stellar mass.
          transferRate               =max(         0.0d0         ,disk    %massStellar         (                         ))/barInstabilityTimescale
          call                                      disk    %massStellarRate         (-           transferRate                              )
          call                                      spheroid%massStellarRate         (+           transferRate ,interrupt,interruptProcedure)
          ! Angular momentum.
          transferRate               =max(         0.0d0         ,disk    %angularMomentum     (                         ))/barInstabilityTimescale
          call                                      disk    %angularMomentumRate     (-           transferRate                              )
          call                                      spheroid%angularMomentumRate     (+           transferRate ,interrupt,interruptProcedure)
          ! Gas abundances.
          fuelAbundancesRates        =max(zeroAbundances         ,disk    %abundancesGas       (                         ))/barInstabilityTimescale
          call                                      disk    %abundancesGasRate       (-     fuelAbundancesRates                             )
          call                                      spheroid%abundancesGasRate       (+     fuelAbundancesRates,interrupt,interruptProcedure)
          ! Stellar abundances.
          stellarAbundancesRates     =max(zeroAbundances         ,disk    %abundancesStellar   (                         ))/barInstabilityTimescale
          call                                      disk    %abundancesStellarRate   (-  stellarAbundancesRates                             )
          call                                      spheroid%abundancesStellarRate   (+  stellarAbundancesRates,interrupt,interruptProcedure)
          ! Stellar luminosities.
          if (.not.diskLuminositiesStellarInactive .or. propertyType /= propertyTypeInactive) then
             luminositiesTransferRate=max(zeroStellarLuminosities,disk    %luminositiesStellar (                         ))/barInstabilityTimescale
             call                                   disk    %luminositiesStellarRate (-luminositiesTransferRate                             )
             call                                   spheroid%luminositiesStellarRate (+luminositiesTransferRate,interrupt,interruptProcedure)
          end if
          ! Stellar properties history.
          historyTransferRate=disk%stellarPropertiesHistory()
          if (historyTransferRate%exists()) then
             historyTransferRate=historyTransferRate/barInstabilityTimescale
             call disk    %stellarPropertiesHistoryRate(-historyTransferRate                             )
             call spheroid%stellarPropertiesHistoryRate(+historyTransferRate,interrupt,interruptProcedure)
          end if
          call historyTransferRate%destroy()
          ! Star formation history.
          historyTransferRate=disk%starFormationHistory()
          if (historyTransferRate%exists()) then
             historyTransferRate=historyTransferRate/barInstabilityTimescale
             call disk    %starFormationHistoryRate(-historyTransferRate                             )
             call spheroid%starFormationHistoryRate(+historyTransferRate,interrupt,interruptProcedure)
          end if
          call historyTransferRate%destroy()
          ! Additional external torque.
          if     (                                                                                                                                                            &
               &   spheroid%angularMomentum() < (spheroid%massGas()+spheroid%massStellar())*darkMatterHaloScale_%virialRadius(node)*darkMatterHaloScale_%virialVelocity(node) &
               &  .and.                                                                                                                                                       &
               &   spheroid%radius         () <                                             darkMatterHaloScale_%virialRadius(node)                                           &
               & ) then
             call spheroid%angularMomentumRate(+barInstabilitySpecificTorque*(spheroid%massGas()+spheroid%massStellar()),interrupt,interruptProcedure)
          end if
       end if

       ! Apply mass loss rate due to ram pressure stripping.
       if (disk%massGas() > 0.0d0) then
          massLossRate=ramPressureStrippingDisks_%rateMassLoss(node)
          if (massLossRate > 0.0d0) then
             hotHalo => node%hotHalo()
             call    disk%                  massGasRate(-massLossRate                                                           )
             call    disk%          angularMomentumRate(-massLossRate*disk%angularMomentum()/(disk%massGas()+disk%massStellar()))
             call    disk%            abundancesGasRate(-massLossRate*disk%abundancesGas  ()/ disk%massGas()                    )
             call hotHalo%           outflowingMassRate(+massLossRate                                                           )
             call hotHalo%outflowingAngularMomentumRate(+massLossRate*disk%angularMomentum()/(disk%massGas()+disk%massStellar()))
             call hotHalo%outflowingAbundancesRate     (+massLossRate*disk%abundancesGas  ()/ disk%massGas()                    )
          end if
       end if

       ! Apply mass loss rate due to tidal stripping.
       if (disk%massGas()+disk%massStellar() > 0.0d0) then
          massLossRate=tidalStrippingDisks_%rateMassLoss(node)
          if (massLossRate > 0.0d0) then
             hotHalo    => node%hotHalo()
             fractionGas    =  min(1.0d0,max(0.0d0,disk%massGas()/(disk%massGas()+disk%massStellar())))
             fractionStellar=  1.0d0-fractionGas
             if (fractionGas     > 0.0d0 .and. disk%massGas    () > 0.0d0) then
                call    disk%                  massGasRate(-fractionGas    *massLossRate                                                             )
                call    disk%            abundancesGasRate(-fractionGas    *massLossRate*disk%abundancesGas    ()/ disk%massGas()                    )
                call hotHalo%           outflowingMassRate(+fractionGas    *massLossRate                                                             )
                call hotHalo%outflowingAbundancesRate     (+fractionGas    *massLossRate*disk%abundancesGas    ()/ disk%massGas()                    )
                call hotHalo%outflowingAngularMomentumRate(+fractionGas    *massLossRate*disk%angularMomentum  ()/(disk%massGas()+disk%massStellar()))
             end if
             if (fractionStellar > 0.0d0 .and. disk%massStellar() > 0.0d0) then
                ! If luminosities are being treated as inactive properties this is an error - they appear on the right-hand side
                ! of the following ODE terms so are not inactive. (An approach similar to what is used for transfer of
                ! luminosities to the spheroid by bar instabilities could work here.)
                if (propertyType == propertyTypeActive .and. diskLuminositiesStellarInactive) call Galacticus_Error_Report('tidal mass loss not supported for inactive luminosity calculation'//{introspection:location})
                ! Stellar mass and metals.
                call    disk%              massStellarRate(-fractionStellar*massLossRate                                                             )
                call    disk%        abundancesStellarRate(-fractionStellar*massLossRate*disk%abundancesStellar()/                disk%massStellar() )
                ! Stellar luminosities.
                luminositiesTransferRate=max(zeroStellarLuminosities,disk%luminositiesStellar())
                call    disk%      luminositiesStellarRate(-fractionStellar*massLossRate*luminositiesTransferRate/                disk%massStellar() )
                ! Stellar properties history.
                historyTransferRate=disk%stellarPropertiesHistory()
                if (historyTransferRate%exists()) then
                   call disk%stellarPropertiesHistoryRate (-fractionStellar*massLossRate*historyTransferRate     /                disk%massStellar() )
                end if
                call historyTransferRate%destroy()
                ! Star formation history.
                historyTransferRate=disk%starFormationHistory()
                if (historyTransferRate%exists()) then
                   call disk    %starFormationHistoryRate (-fractionStellar*massLossRate*historyTransferRate     /                disk%massStellar() )
                end if
             end if
             call       disk%          angularMomentumRate(-                massLossRate*disk%angularMomentum  ()/(disk%massGas()+disk%massStellar()))
          end if
       end if
    end select

    ! Return the procedure pointer.
    interruptProcedureReturn => interruptProcedure

    return
  end subroutine Node_Component_Disk_Standard_Rate_Compute

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Disk_Standard_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Disk_Standard_Scale_Set(node)
    !% Set scales for properties of {\normalfont \ttfamily node}.
    use :: Abundances_Structure          , only : abs                 , abundances               , max                  , operator(*)            , &
          &                                       unitAbundances
    use :: Galacticus_Nodes              , only : nodeComponentDisk   , nodeComponentDiskStandard, nodeComponentSpheroid, treeNode               , &
         &                                        defaultDiskComponent
    use :: Histories                     , only : history
    use :: Stellar_Luminosities_Structure, only : abs                 , max                      , stellarLuminosities  , unitStellarLuminosities
    implicit none
    type            (treeNode                        ), intent(inout), pointer :: node
    class           (nodeComponentDisk               )               , pointer :: disk
    class           (nodeComponentSpheroid           )               , pointer :: spheroid
    double precision                                  , parameter              :: massMinimum                   =1.0d+00
    double precision                                  , parameter              :: angularMomentumMinimum        =1.0d-1
    double precision                                  , parameter              :: fractionTolerance             =1.0d-4
    double precision                                  , parameter              :: luminosityMinimum             =1.0d0
    double precision                                                           :: angularMomentum                      , mass
    type            (history                         )                         :: stellarPopulationHistoryScales
    type            (stellarLuminosities             )                         :: stellarLuminositiesScale
    type            (abundances                      )                         :: abundancesScale

    ! Check if we are the default method.
    if (.not.defaultDiskComponent%standardIsActive()) return
    ! Get the disk component.
    disk => node%disk()
    ! Check if an standard disk component exists.
    select type (disk)
    class is (nodeComponentDiskStandard)
       ! Get spheroid component.
       spheroid => node%spheroid()
       ! Set scale for angular momentum.
       angularMomentum=abs(disk%angularMomentum())
       call disk%angularMomentumScale(max(angularMomentum,angularMomentumMinimum))
       ! Set scale for masses.
       mass           =abs(                         &
            &              +abs(disk%massGas    ()) &
            &              +abs(disk%massStellar()) &
            &             )
       call disk%massGasScale                    (max(mass,massMinimum      ))
       call disk%massStellarScale                (max(mass,massMinimum      ))
       call disk%massStellarFormedScale          (max(mass,massMinimum      ))
       ! Set the scale for the retained stellar mass fraction.
       call disk%fractionMassRetainedScale(fractionTolerance*disk%fractionMassRetained())
       ! Set scales for abundances if necessary.
       if (abundancesCount > 0) then
          ! Set scale for abundances.
          abundancesScale=+max(                                     &
               &               +abs(+disk    %abundancesGas    ())  &
               &               +abs(+disk    %abundancesStellar()), &
               &                    +massMinimum                    &
               &                    *unitAbundances                 &
               &              )
          ! Set scale for gas abundances.
          call disk%abundancesGasScale    (abundancesScale)

          ! Set scale for stellar abundances.
          call disk%abundancesStellarScale(abundancesScale)
       end if
       ! Set scale for stellar luminosities.
       stellarLuminositiesScale=max(                                        &
            &                       abs(+disk      %luminositiesStellar()), &
            &                           +unitStellarLuminosities            &
            &                           *luminosityMinimum                  &
            &                      )
       call stellarLuminositiesScale%truncate                (disk       %luminositiesStellar())
       call disk       %luminositiesStellarScale(stellarLuminositiesScale                      )

       ! Set scales for stellar population properties and star formation histories.
       stellarPopulationHistoryScales=disk%stellarPropertiesHistory()
       call stellarPopulationProperties_%scales (disk%massStellar(),disk%abundancesStellar(),stellarPopulationHistoryScales)
       call disk%stellarPropertiesHistoryScale  (                                            stellarPopulationHistoryScales)
       call stellarPopulationHistoryScales%destroy()
       stellarPopulationHistoryScales=disk%starFormationHistory()
       call starFormationHistory_%scales        (stellarPopulationHistoryScales,disk%massStellar(),disk%abundancesStellar())
       call disk%starFormationHistoryScale      (stellarPopulationHistoryScales                                            )
       call stellarPopulationHistoryScales%destroy()

    end select
    return
  end subroutine Node_Component_Disk_Standard_Scale_Set

  !# <inactiveSetTask>
  !#  <unitName>Node_Component_Disk_Standard_Inactive</unitName>
  !# </inactiveSetTask>
  subroutine Node_Component_Disk_Standard_Inactive(node)
    !% Set Jacobian zero status for properties of {\normalfont \ttfamily node}.
    use :: Galacticus_Nodes, only : nodeComponentDisk, nodeComponentDiskStandard, treeNode
    implicit none
    type (treeNode         ), intent(inout), pointer :: node
    class(nodeComponentDisk)               , pointer :: disk

    ! Get the disk component.
    disk => node%disk()
    ! Check if an standard disk component exists.
    select type (disk)
    class is (nodeComponentDiskStandard)
       if (diskLuminositiesStellarInactive) call disk%luminositiesStellarInactive()
    end select
    return
  end subroutine Node_Component_Disk_Standard_Inactive

  subroutine satelliteMerger(self,node)
    !% Transfer any standard disk associated with {\normalfont \ttfamily node} to its host halo.
    use :: Abundances_Structure             , only : zeroAbundances
    use :: Galacticus_Error                 , only : Galacticus_Error_Report
    use :: Galacticus_Nodes                 , only : nodeComponentDisk      , nodeComponentDiskStandard, nodeComponentSpheroid, treeNode
    use :: Histories                        , only : history
    use :: Satellite_Merging_Mass_Movements , only : destinationMergerDisk  , destinationMergerSpheroid
    use :: Stellar_Luminosities_Structure   , only : zeroStellarLuminosities
    implicit none
    class           (*                    ), intent(inout) :: self
    type            (treeNode             ), intent(inout) :: node
    class           (nodeComponentDisk    ), pointer       :: diskHost               , disk
    class           (nodeComponentSpheroid), pointer       :: spheroidHost           , spheroid
    type            (treeNode             ), pointer       :: nodeHost
    type            (history              )                :: historyHost            , historyNode
    double precision                                       :: specificAngularMomentum
    integer                                                :: destinationGasSatellite, destinationGasHost       , &
         &                                                    destinationStarsHost   , destinationStarsSatellite
    logical                                                :: mergerIsMajor
    !$GLC attributes unused :: self

    ! Check that the disk is of the standard class.
    disk => node%disk()
    select type (disk)
    class is (nodeComponentDiskStandard)
       spheroid => node%spheroid()
       ! Find the node to merge with.
       nodeHost     => node%mergesWith  (                 )
       diskHost     => nodeHost%disk    (autoCreate=.true.)
       spheroidHost => nodeHost%spheroid(autoCreate=.true.)
       ! Get specific angular momentum of the disk material.
       if (                                               disk%massGas()+disk%massStellar() > 0.0d0) then
          specificAngularMomentum=disk%angularMomentum()/(disk%massGas()+disk%massStellar())
       else
          specificAngularMomentum=0.0d0
       end if
       ! Get mass movement descriptors.
       call mergerMassMovements_%get(node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
       ! Move the gas component of the standard disk to the host.
       select case (destinationGasSatellite)
       case (destinationMergerDisk)
          call diskHost    %massGasSet            (                                                                     &
               &                                             diskHost    %massGas            ()                         &
               &                                            +disk        %massGas            ()                         &
               &                                           )
          call diskHost    %abundancesGasSet      (                                                                     &
               &                                             diskHost    %abundancesGas      ()                         &
               &                                            +disk        %abundancesGas      ()                         &
               &                                           )
          call diskHost    %angularMomentumSet    (                                                                     &
               &                                             diskHost    %angularMomentum    ()                         &
               &                                            +disk        %massGas            ()*specificAngularMomentum &
               &                                           )
       case (destinationMergerSpheroid)
          call spheroidHost%massGasSet            (                                                                     &
               &                                             spheroidHost%massGas            ()                         &
               &                                            +disk        %massGas            ()                         &
               &                                           )
          call spheroidHost%abundancesGasSet      (                                                                     &
               &                                             spheroidHost%abundancesGas      ()                         &
               &                                            +disk        %abundancesGas      ()                         &
               &                                           )
       case default
          call Galacticus_Error_Report('unrecognized movesTo descriptor'//{introspection:location})
       end select
       call disk%      massGasSet(         0.0d0)
       call disk%abundancesGasSet(zeroAbundances)
       ! Move the stellar component of the standard disk to the host.
       select case (destinationStarsSatellite)
       case (destinationMergerDisk)
          call diskHost    %massStellarSet        (                                                                     &
               &                                             diskHost    %massStellar        ()                         &
               &                                            +disk        %massStellar        ()                         &
               &                                           )
          call diskHost    %abundancesStellarSet  (                                                                     &
               &                                             diskHost    %abundancesStellar  ()                         &
               &                                            +disk        %abundancesStellar  ()                         &
               &                                           )
          call diskHost    %luminositiesStellarSet(                                                                     &
               &                                             diskHost    %luminositiesStellar()                         &
               &                                            +disk        %luminositiesStellar()                         &
               &                                           )
          call diskHost    %angularMomentumSet    (                                                                     &
               &                                             diskHost    %angularMomentum    ()                         &
               &                                            +disk        %massStellar        ()*specificAngularMomentum &
               &                                           )
          ! Also add stellar properties histories.
          historyNode=disk    %stellarPropertiesHistory()
          historyHost=diskHost%stellarPropertiesHistory()
          call historyHost%interpolatedIncrement      (historyNode)
          call historyNode%reset                      (           )
          call diskHost   %stellarPropertiesHistorySet(historyHost)
          call disk       %stellarPropertiesHistorySet(historyNode)
          ! Also add star formation histories.
          historyNode=disk    %starFormationHistory    ()
          historyHost=diskHost%starFormationHistory    ()
          call historyHost%increment              (historyNode,autoExtend  =.true. )
          call historyNode%reset                  (                                )
          call diskHost   %starFormationHistorySet(historyHost                     )
          call disk       %starFormationHistorySet(historyNode                     )
          call historyNode%destroy                (            recordMemory=.false.)
          call historyHost%destroy                (            recordMemory=.false.)
       case (destinationMergerSpheroid)
          call spheroidHost%massStellarSet        (                                                                     &
               &                                             spheroidHost%massStellar        ()                         &
               &                                            +disk        %massStellar        ()                         &
               &                                           )
          call spheroidHost%abundancesStellarSet  (                                                                     &
               &                                             spheroidHost%abundancesStellar  ()                         &
               &                                            +disk        %abundancesStellar  ()                         &
               &                                            )
          call spheroidHost%luminositiesStellarSet(                                                                     &
               &                                             spheroidHost%luminositiesStellar()                         &
               &                                            +disk        %luminositiesStellar()                         &
               &                                           )
          ! Also add stellar properties histories.
          historyNode=disk    %stellarPropertiesHistory()
          historyHost=spheroidHost%stellarPropertiesHistory()
          call historyHost%interpolatedIncrement(historyNode)
          call historyNode%reset    (           )
          call spheroidHost%stellarPropertiesHistorySet(historyHost)
          call disk    %stellarPropertiesHistorySet(historyNode)
          ! Also add star formation histories.
          historyNode=disk    %starFormationHistory    ()
          historyHost=spheroidHost%starFormationHistory    ()
          call historyHost %increment              (historyNode,autoExtend  =.true. )
          call historyNode %reset                  (                                )
          call spheroidHost%starFormationHistorySet(historyHost                     )
          call disk        %starFormationHistorySet(historyNode                     )
          call historyNode %destroy                (            recordMemory=.false.)
          call historyHost %destroy                (            recordMemory=.false.)
       case default
          call Galacticus_Error_Report('unrecognized movesTo descriptor'//{introspection:location})
       end select
       call disk%        massStellarSet(                  0.0d0)
       call disk%  abundancesStellarSet(         zeroAbundances)
       call disk%luminositiesStellarSet(zeroStellarLuminosities)
       call disk%    angularMomentumSet(                  0.0d0)
    end select
    return
  end subroutine satelliteMerger

  !# <radiusSolverPlausibility>
  !#  <unitName>Node_Component_Disk_Standard_Radius_Solver_Plausibility</unitName>
  !#  <after>Node_Component_Basic_Standard_Plausibility</after>
  !# </radiusSolverPlausibility>
  subroutine Node_Component_Disk_Standard_Radius_Solver_Plausibility(node)
    !% Determines whether the disk is physically plausible for radius solving tasks. Require that it have non-zero mass and angular momentum.
    use :: Galacticus_Nodes, only : defaultDiskComponent, nodeComponentDisk, nodeComponentDiskStandard, treeNode
    implicit none
    type            (treeNode         ), intent(inout) :: node
    class           (nodeComponentDisk), pointer       :: disk
    double precision                                   :: angularMomentumScale

    ! Return immediately if our method is not selected.
    if     (                                                &
         &  .not.                                           &
         &   (                                              &
         &     defaultDiskComponent%standardIsActive     () &
         &    .and.                                         &
         &     node                %isPhysicallyPlausible   &
         &    .and.                                         &
         &     node                %isSolvable              &
         &   )                                              &
         & ) return

    ! Determine the plausibility of the current disk.
    disk => node%disk()
    select type (disk)
       class is (nodeComponentDiskStandard)
       if      (disk%angularMomentum()                             <                       0.0d0) &
            & node%isPhysicallyPlausible=.false.
       if      (disk%massStellar    ()+disk%massGas() <  -diskMassToleranceAbsolute) then
          node%isPhysicallyPlausible=.false.
       else if (disk%massStellar    ()+disk%massGas() >=                      0.0d0) then
          if      (                                                                               &
               &   disk%angularMomentum() < 0.0d0                                                 &
               &  ) then
             node%isPhysicallyPlausible=.false.
          else
             angularMomentumScale=(                                           &
                  &                 disk%massStellar()                        &
                  &                +disk%massGas    ()                        &
                  &               )                                           &
                  &               * darkMatterHaloScale_%virialRadius  (node) &
                  &               * darkMatterHaloScale_%virialVelocity(node)
             if     (                                                                      &
                  &   disk%angularMomentum() > angularMomentumMaximum*angularMomentumScale &
                  &  .or.                                                                  &
                  &   disk%angularMomentum() < angularMomentumMinimum*angularMomentumScale &
                  & ) then
                ! Ignore disks with angular momenta greatly exceeding that which would be expected if they had a radius comparable to the
                ! virial radius of their halo.
                node%isPhysicallyPlausible=.false.
             end if
          end if
       end if
    end select

    ! Reset the record of trial radii - negative values indicate that the entries have not yet been set to physically meaningful
    ! values.
    radiusHistory        =-1.0d0
    radiusSolverIteration= 0
    return
  end subroutine Node_Component_Disk_Standard_Radius_Solver_Plausibility

  double precision function Node_Component_Disk_Standard_Radius_Solve(node)
    !% Return the radius of the standard disk used in structure solvers.
    use :: Galacticus_Nodes, only : nodeComponentDisk, treeNode
    implicit none
    type (treeNode         ), intent(inout) :: node
    class(nodeComponentDisk), pointer       :: disk

    disk => node%disk()
    Node_Component_Disk_Standard_Radius_Solve=disk%radius()*diskStructureSolverRadius
    return
  end function Node_Component_Disk_Standard_Radius_Solve

  subroutine Node_Component_Disk_Standard_Radius_Solve_Set(node,radius)
    !% Set the radius of the standard disk used in structure solvers.
    use :: Galacticus_Nodes, only : nodeComponentDisk, treeNode
    implicit none
    type            (treeNode         ), intent(inout) :: node
    double precision                   , intent(in   ) :: radius
    class           (nodeComponentDisk), pointer       :: disk

    disk => node%disk()
    call disk%radiusSet(max(radius,0.0d0)/diskStructureSolverRadius)
    return
  end subroutine Node_Component_Disk_Standard_Radius_Solve_Set

  double precision function Node_Component_Disk_Standard_Velocity(node)
    !% Return the circular velocity of the standard disk.
    use :: Galacticus_Nodes, only : nodeComponentDisk, treeNode
    implicit none
    type (treeNode         ), intent(inout) :: node
    class(nodeComponentDisk), pointer       :: disk

    disk => node%disk()
    Node_Component_Disk_Standard_Velocity=disk%velocity()
    return
  end function Node_Component_Disk_Standard_Velocity

  subroutine Node_Component_Disk_Standard_Velocity_Set(node,velocity)
    !% Set the circular velocity of the standard disk.
    use :: Galacticus_Nodes, only : nodeComponentDisk, treeNode
    implicit none
    type            (treeNode         ), intent(inout) :: node
    double precision                   , intent(in   ) :: velocity
    class           (nodeComponentDisk), pointer       :: disk

    disk => node%disk()
    call disk%velocitySet(velocity)
    return
  end subroutine Node_Component_Disk_Standard_Velocity_Set

  !# <radiusSolverTask>
  !#  <unitName>Node_Component_Disk_Standard_Radius_Solver</unitName>
  !# </radiusSolverTask>
  subroutine Node_Component_Disk_Standard_Radius_Solver(node,componentActive,specificAngularMomentumRequired,specificAngularMomentum,Radius_Get,Radius_Set,Velocity_Get&
       &,Velocity_Set)
    !% Interface for the size solver algorithm.
    use :: Galacticus_Nodes            , only : nodeComponentDisk              , nodeComponentDiskStandard, treeNode
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    type            (treeNode                                     ), intent(inout)          :: node
    logical                                                        , intent(  out)          :: componentActive
    logical                                                        , intent(in   )          :: specificAngularMomentumRequired
    double precision                                               , intent(  out)          :: specificAngularMomentum
    procedure       (Node_Component_Disk_Standard_Radius_Solve    ), intent(  out), pointer :: Radius_Get                     , Velocity_Get
    procedure       (Node_Component_Disk_Standard_Radius_Solve_Set), intent(  out), pointer :: Radius_Set                     , Velocity_Set
    class           (nodeComponentDisk                            )               , pointer :: disk
    double precision                                                                        :: angularMomentum                , diskMass    , &
         &                                                                                     specificAngularMomentumMean

    ! Determine if node has an active disk component supported by this module.
    componentActive        =  .false.
    specificAngularMomentum=  0.0d0
    disk      => node%disk()
    select type (disk)
       class is (nodeComponentDiskStandard)
       componentActive=.true.
       ! Get the angular momentum.
       if (specificAngularMomentumRequired) then
          angularMomentum=disk%angularMomentum()
          if (angularMomentum >= 0.0d0) then
             ! Compute the specific angular momentum at the scale radius, assuming a flat rotation curve.
             diskMass= disk%massGas    () &
                  &   +disk%massStellar()
             if (diskMass > 0.0d0) then
                specificAngularMomentumMean=angularMomentum/diskMass
             else
                specificAngularMomentumMean=0.0d0
             end if
             specificAngularMomentum=specificAngularMomentumMean*diskStructureSolverSpecificAngularMomentum
             ! If using the Cole et al. (2000) method for disk radii, adjust the specific angular momentum to account for the
             ! difference between rotation curves for thin disk and a spherical mass distribution. Trap instances where this leads to
             ! imaginary specific angular momentum - this can happen as the radius solver explores the allowed range of radii when
             ! seeking a solution.
             if (diskRadiusSolverCole2000Method)                                                       &
                  & specificAngularMomentum=sqrt(                                                      &
                  &                               max(                                                 &
                  &                                    0.0d0,                                          &
                  &                                    specificAngularMomentum**2                      &
                  &                                   -diskRadiusSolverFlatVsSphericalFactor           &
                  &                                   *gravitationalConstantGalacticus                 &
                  &                                   *diskMass                                        &
                  &                                   *Node_Component_Disk_Standard_Radius_Solve(node) &
                  &                                  )                                                 &
                  )
          end if
          ! Associate the pointers with the appropriate property routines.
          Radius_Get   => Node_Component_Disk_Standard_Radius_Solve
          Radius_Set   => Node_Component_Disk_Standard_Radius_Solve_Set
          Velocity_Get => Node_Component_Disk_Standard_Velocity
          Velocity_Set => Node_Component_Disk_Standard_Velocity_Set
       else
          call Node_Component_Disk_Standard_Radius_Solve_Set(node,0.0d0)
          call Node_Component_Disk_Standard_Velocity_Set    (node,0.0d0)
          componentActive=.false.
       end if
    end select
    return
  end subroutine Node_Component_Disk_Standard_Radius_Solver

  !# <mergerTreeExtraOutputTask>
  !#  <unitName>Node_Component_Disk_Standard_Star_Formation_History_Output</unitName>
  !# </mergerTreeExtraOutputTask>
  subroutine Node_Component_Disk_Standard_Star_Formation_History_Output(node,iOutput,treeIndex,nodePassesFilter)
    !% Store the star formation history in the output file.
    use            :: Galacticus_Nodes, only : nodeComponentDisk, nodeComponentDiskStandard, treeNode, defaultDiskComponent
    use            :: Histories       , only : history
    use, intrinsic :: ISO_C_Binding   , only : c_size_t
    use            :: Kind_Numbers    , only : kind_int8
    implicit none
    type   (treeNode         ), intent(inout), pointer :: node
    integer(c_size_t         ), intent(in   )          :: iOutput
    integer(kind=kind_int8   ), intent(in   )          :: treeIndex
    logical                   , intent(in   )          :: nodePassesFilter
    class  (nodeComponentDisk)               , pointer :: disk
    type   (history          )                         :: historyStarFormation

    ! Check if we are the default method.
    if (.not.defaultDiskComponent%standardIsActive()) return
    ! Output the star formation history if a disk exists for this component.
    disk => node%disk()
    select type (disk)
       class is (nodeComponentDiskStandard)
       historyStarFormation=disk%starFormationHistory()
       call starFormationHistory_%output(node,nodePassesFilter,historyStarFormation,iOutput,treeIndex,'disk')
       call disk%starFormationHistorySet(historyStarFormation)
    end select
    return
  end subroutine Node_Component_Disk_Standard_Star_Formation_History_Output

  !# <galacticusStateStoreTask>
  !#  <unitName>Node_Component_Disk_Standard_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Node_Component_Disk_Standard_State_Store(stateFile,gslStateFile,stateOperationID)
    !% Write the tablulation state to file.
    use            :: Galacticus_Display               , only : Galacticus_Display_Message, verbosityInfo
    use, intrinsic :: ISO_C_Binding                    , only : c_size_t                  , c_ptr
    use            :: Node_Component_Disk_Standard_Data, only : diskMassDistribution
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call Galacticus_Display_Message('Storing state for: treeNodeMethodDisk -> standard',verbosity=verbosityInfo)
    !# <workaround type="gfortran" PR="92836" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=92836">
    !#  <description>Internal file I/O in gfortran can be non-thread safe.</description>
    !# </workaround>
#ifdef THREADSAFEIO
    !$omp critical(gfortranInternalIO)
#endif
    write (stateFile) diskStructureSolverSpecificAngularMomentum,diskRadiusSolverFlatVsSphericalFactor
    write (stateFile) associated(diskMassDistribution)
#ifdef THREADSAFEIO
    !$omp end critical(gfortranInternalIO)
#endif
    if (associated(diskMassDistribution)) call diskMassDistribution%stateStore(stateFile,gslStateFile,stateOperationID)
    return
  end subroutine Node_Component_Disk_Standard_State_Store

  !# <galacticusStateRetrieveTask>
  !#  <unitName>Node_Component_Disk_Standard_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Node_Component_Disk_Standard_State_Retrieve(stateFile,gslStateFile,stateOperationID)
    !% Retrieve the tabulation state from the file.
    use            :: Galacticus_Display               , only : Galacticus_Display_Message, verbosityInfo
    use            :: Galacticus_Error                 , only : Galacticus_Error_Report
    use, intrinsic :: ISO_C_Binding                    , only : c_size_t                  , c_ptr
    use            :: Node_Component_Disk_Standard_Data, only : diskMassDistribution
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile
    logical                          :: wasAllocated

    call Galacticus_Display_Message('Retrieving state for: treeNodeMethodDisk -> standard',verbosity=verbosityInfo)
    !# <workaround type="gfortran" PR="92836" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=92836">
    !#  <description>Internal file I/O in gfortran can be non-thread safe.</description>
    !# </workaround>
#ifdef THREADSAFEIO
    !$omp critical(gfortranInternalIO)
#endif
    read (stateFile) diskStructureSolverSpecificAngularMomentum,diskRadiusSolverFlatVsSphericalFactor
    read (stateFile) wasAllocated
#ifdef THREADSAFEIO
    !$omp end critical(gfortranInternalIO)
#endif
    if (wasAllocated) then
       if (.not.associated(diskMassDistribution)) call Galacticus_Error_Report('diskMassDistribution was stored, but is now not allocated'//{introspection:location})
       call diskMassDistribution%stateRestore(stateFile,gslStateFile,stateOperationID)
    end if
    return
  end subroutine Node_Component_Disk_Standard_State_Retrieve

end module Node_Component_Disk_Standard
