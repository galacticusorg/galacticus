!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
Contains a module which implements the standard disk node component.
!!}

module Node_Component_Disk_Standard
  !!{
  Implements the standard disk node component.
  !!}
  use :: Dark_Matter_Halo_Scales         , only : darkMatterHaloScaleClass
  use :: Satellite_Merging_Mass_Movements, only : mergerMassMovementsClass
  use :: Star_Formation_Histories        , only : starFormationHistory            , starFormationHistoryClass
  use :: Stellar_Population_Properties   , only : stellarPopulationPropertiesClass
  use :: Galactic_Structure              , only : galacticStructureClass
  implicit none
  private
  public :: Node_Component_Disk_Standard_Scale_Set                    , Node_Component_Disk_Standard_Pre_Evolve                  , &
       &    Node_Component_Disk_Standard_Radius_Solver_Plausibility   , Node_Component_Disk_Standard_Radius_Solver               , &
       &    Node_Component_Disk_Standard_Star_Formation_History_Output, Node_Component_Disk_Standard_Thread_Uninitialize         , &
       &    Node_Component_Disk_Standard_Initialize                   , Node_Component_Disk_Standard_Calculation_Reset           , &
       &    Node_Component_Disk_Standard_State_Store                  , Node_Component_Disk_Standard_State_Retrieve              , &
       &    Node_Component_Disk_Standard_Thread_Initialize            , Node_Component_Disk_Standard_Inactive                    , &
       &    Node_Component_Disk_Standard_Post_Step                    , Node_Component_Disk_Standard_Star_Formation_History_Flush

  !![
  <component>
   <class>disk</class>
   <name>standard</name>
   <isDefault>true</isDefault>
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
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output unitsInSI="massSolar" comment="Mass of stars in the standard disk."/>
    </property>
    <property>
      <name>massStellarFormed</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
    </property>
    <property>
      <name>fractionMassRetained</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
    </property>
    <property>
      <name>abundancesStellar</name>
      <type>abundances</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output unitsInSI="massSolar" comment="Mass of metals in the stellar phase of the standard disk."/>
    </property>
    <property>
      <name>massGas</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" makeGeneric="true" />
      <output unitsInSI="massSolar" comment="Mass of gas in the standard disk."/>
    </property>
    <property>
      <name>abundancesGas</name>
      <type>abundances</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" makeGeneric="true" />
      <output unitsInSI="massSolar" comment="Mass of metals in the gas phase of the standard disk."/>
    </property>
    <property>
      <name>angularMomentum</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" makeGeneric="true" />
      <output unitsInSI="massSolar*megaParsec*kilo" comment="Angular momentum of the standard disk."/>
    </property>
    <property>
      <name>radius</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
      <output unitsInSI="megaparsec" comment="Radial scale length in the standard disk."/>
    </property>
    <property>
      <name>halfMassRadius</name>
      <attributes isSettable="false" isGettable="true" isEvolvable="false" isVirtual="true" />
      <type>double</type>
      <rank>0</rank>
      <getFunction>Node_Component_Disk_Standard_Half_Mass_Radius</getFunction>
    </property>
    <property>
      <name>velocity</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
      <output unitsInSI="kilo" comment="Circular velocity of the standard disk at scale length."/>
    </property>
    <property>
      <name>luminositiesStellar</name>
      <type>stellarLuminosities</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output unitsInSI="luminosityZeroPointAB" comment="Luminosity of disk stars."/>
    </property>
    <property>
      <name>stellarPropertiesHistory</name>
      <type>history</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
    </property>
    <property>
      <name>starFormationHistory</name>
      <type>history</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
    </property>
   </properties>
   <bindings>
    <binding method="attachPipes"               function="Node_Component_Disk_Standard_Attach_Pipes"              bindsTo="component" description="Attach pipes to the standard disk component." returnType="\void" arguments=""/>
    <binding method="enclosedMass"              function="Node_Component_Disk_Standard_Enclosed_Mass"             bindsTo="component"                                                                                           />
    <binding method="acceleration"              function="Node_Component_Disk_Standard_Acceleration"              bindsTo="component"                                                                                           />
    <binding method="tidalTensor"               function="Node_Component_Disk_Standard_Tidal_Tensor"              bindsTo="component"                                                                                           />
    <binding method="density"                   function="Node_Component_Disk_Standard_Density"                   bindsTo="component"                                                                                           />
    <binding method="densitySphericalAverage"   function="Node_Component_Disk_Standard_Density_Spherical_Average" bindsTo="component"                                                                                           />
    <binding method="potential"                 function="Node_Component_Disk_Standard_Potential"                 bindsTo="component"                                                                                           />
    <binding method="rotationCurve"             function="Node_Component_Disk_Standard_Rotation_Curve"            bindsTo="component"                                                                                           />
    <binding method="rotationCurveGradient"     function="Node_Component_Disk_Standard_Rotation_Curve_Gradient"   bindsTo="component"                                                                                           />
    <binding method="surfaceDensity"            function="Node_Component_Disk_Standard_Surface_Density"           bindsTo="component"                                                                                           />
    <binding method="chandrasekharIntegral"   isDeferred="true"                                                   bindsTo="component"                                                                                            >
      <interface>
	<type>double</type>
	<shape>3</shape>
	<self pass="true" intent="inout" />
	<argument>type            (treeNode                    ), intent(inout)               :: nodeSatellite                       </argument>
	<argument>double precision                              , intent(in   ), dimension(3) :: positionCartesian, velocityCartesian</argument>
	<argument>type            (enumerationComponentTypeType), intent(in   )               :: componentType                       </argument>
	<argument>type            (enumerationMassTypeType     ), intent(in   )               :: massType                            </argument>
      </interface>
    </binding>
   </bindings>
   <functions>objects.nodes.components.disk.standard.bound_functions.inc</functions>
  </component>
  !!]

  ! Objects used by this component.
  class(darkMatterHaloScaleClass        ), pointer :: darkMatterHaloScale_
  class(stellarPopulationPropertiesClass), pointer :: stellarPopulationProperties_
  class(starFormationHistoryClass       ), pointer :: starFormationHistory_
  class(mergerMassMovementsClass        ), pointer :: mergerMassMovements_
  class(galacticStructureClass          ), pointer :: galacticStructure_
  !$omp threadprivate(darkMatterHaloScale_,stellarPopulationProperties_,starFormationHistory_,mergerMassMovements_,galacticStructure_)

  ! Internal count of abundances.
  integer                                     :: abundancesCount

  ! Parameters controlling the physical implementation.
  double precision                            :: toleranceAbsoluteMass                              , radiusStructureSolver               , &
       &                                         toleranceRelativeMetallicity
  logical                                     :: diskNegativeAngularMomentumAllowed                 , structureSolverUseCole2000Method    , &
       &                                         inactiveLuminositiesStellar

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

  ! A threadprivate object used to track to which thread events are attached.
  integer :: thread
  !$omp threadprivate(thread)

contains

  !![
  <nodeComponentInitializationTask>
   <unitName>Node_Component_Disk_Standard_Initialize</unitName>
  </nodeComponentInitializationTask>
  !!]
  subroutine Node_Component_Disk_Standard_Initialize(parameters)
    !!{
    Initializes the tree node standard disk methods module.
    !!}
    use :: Abundances_Structure, only : Abundances_Property_Count
    use :: Error               , only : Error_Report
    use :: Galacticus_Nodes    , only : defaultDiskComponent     , nodeComponentDiskStandard
    use :: Input_Parameters    , only : inputParameter           , inputParameters
    implicit none
    type(inputParameters          ), intent(inout) :: parameters
    type(nodeComponentDiskStandard)                :: diskStandardComponent
    type(inputParameters          )                :: subParameters

    if (defaultDiskComponent%standardIsActive()) then
       ! Get number of abundance properties.
       abundancesCount  =Abundances_Property_Count            ()
       ! Attach the cooling mass/angular momentum pipes from the hot halo component.
       if (.not.pipesAttached) then
          call diskStandardComponent%attachPipes()
          pipesAttached=.true.
       end if
       ! Bind the Chandrasekhar integral function.
       call diskStandardComponent%chandrasekharIntegralFunction(Node_Component_Disk_Standard_Chandrasekhar_Integral)
       ! Find our parameters.
       subParameters=parameters%subParameters('componentDisk')
       ! Read parameters controlling the physical implementation.
       !![
       <inputParameter>
         <name>toleranceAbsoluteMass</name>
         <defaultValue>1.0d-6</defaultValue>
         <description>The mass tolerance used to judge whether the disk is physically plausible.</description>
         <source>subParameters</source>
       </inputParameter>
       <inputParameter>
         <name>toleranceRelativeMetallicity</name>
         <defaultValue>1.0d-4</defaultValue>
         <description>The metallicity tolerance for ODE solution.</description>
         <source>subParameters</source>
       </inputParameter>
       <inputParameter>
         <name>radiusStructureSolver</name>
         <defaultValue>1.0d0</defaultValue>
         <description>The radius (in units of the standard scale length) to use in solving for the size of the disk.</description>
         <source>subParameters</source>
       </inputParameter>
       <inputParameter>
         <name>structureSolverUseCole2000Method</name>
         <defaultValue>.false.</defaultValue>
         <description></description>
         <source>subParameters</source>
       </inputParameter>
       <inputParameter>
         <name>diskNegativeAngularMomentumAllowed</name>
         <defaultValue>.true.</defaultValue>
         <description>Specifies whether or not negative angular momentum is allowed for the disk.</description>
         <source>subParameters</source>
       </inputParameter>
       <inputParameter>
         <name>inactiveLuminositiesStellar</name>
         <defaultValue>.false.</defaultValue>
         <description>Specifies whether or not disk stellar luminosities are inactive properties (i.e. do not appear in any ODE being solved).</description>
         <source>subParameters</source>
       </inputParameter>
       !!]
    end if
    return
  end subroutine Node_Component_Disk_Standard_Initialize

  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Node_Component_Disk_Standard_Thread_Initialize</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Node_Component_Disk_Standard_Thread_Initialize(parameters)
    !!{
    Initializes the standard disk component module for each thread.
    !!}
    use :: Events_Hooks                     , only : dependencyDirectionAfter   , dependencyRegEx      , openMPThreadBindingAtLevel, postEvolveEvent, &
          &                                          satelliteMergerEvent
    use :: Error                            , only : Error_Report
    use :: Galacticus_Nodes                 , only : defaultDiskComponent
    use :: Input_Parameters                 , only : inputParameter             , inputParameters
    use :: Mass_Distributions               , only : massDistributionCylindrical
    use :: Node_Component_Disk_Standard_Data, only : massDistributionDisk       , massDistributionDisk_        
    implicit none
    type            (inputParameters), intent(inout) :: parameters
    type            (dependencyRegEx), dimension(1)  :: dependencies
    double precision                                 :: massDistributionDiskDensityMoment1, massDistributionDiskDensityMoment2
    logical                                          :: surfaceDensityMoment1IsInfinite   , surfaceDensityMoment2IsInfinite
    type            (inputParameters)                :: subParameters

    ! Check if this implementation is selected. If so, initialize the mass distribution.
    if (defaultDiskComponent%standardIsActive()) then
       dependencies(1)=dependencyRegEx(dependencyDirectionAfter,'^remnantStructure:')
       call satelliteMergerEvent%attach(thread,satelliteMerger,openMPThreadBindingAtLevel,label='nodeComponentDiskStandard',dependencies=dependencies)
       call postEvolveEvent     %attach(thread,postEvolve     ,openMPThreadBindingAtLevel,label='nodeComponentDiskStandard'                          )
       ! Find our parameters.
       subParameters=parameters%subParameters('componentDisk')
       !![
       <objectBuilder class="darkMatterHaloScale"                                                 name="darkMatterHaloScale_"         source="subParameters"                    />
       <objectBuilder class="stellarPopulationProperties"                                         name="stellarPopulationProperties_" source="subParameters"                    />
       <objectBuilder class="starFormationHistory"                                                name="starFormationHistory_"        source="subParameters"                    />
       <objectBuilder class="mergerMassMovements"                                                 name="mergerMassMovements_"         source="subParameters"                    />
       <objectBuilder class="galacticStructure"                                                   name="galacticStructure_"           source="subParameters"                    />
       <objectBuilder class="massDistribution"               parameterName="massDistributionDisk" name="massDistributionDisk_"        source="subParameters" threadPrivate="yes" >
        <default>
         <massDistributionDisk value="exponentialDisk">
          <dimensionless value="true"/>
         </massDistributionDisk>
        </default>
       </objectBuilder>
       !!]
       ! Validate the disk mass distribution.
       select type (massDistributionDisk_)
       class is (massDistributionCylindrical)
          ! Since the disk must be cylindrical, deep-copy it to an object of that class. Then we do not need to perform
          ! type-guards elsewhere in the code.
          allocate(massDistributionDisk,mold=massDistributionDisk_)
          !$omp critical(diskStandardDeepCopy)
          !![
	  <deepCopyReset variables="massDistributionDisk_"/>
	  <deepCopy source="massDistributionDisk_" destination="massDistributionDisk"/>
	  <deepCopyFinalize variables="massDistributionDisk"/>  
          !!]
          !$omp end critical(diskStandardDeepCopy)
       class default
          call Error_Report('only cylindrically symmetric mass distributions are allowed'//{introspection:location})
       end select
       if (.not.massDistributionDisk%isDimensionless()) call Error_Report('disk mass distribution must be dimensionless'//{introspection:location})
       ! Compute the specific angular momentum of the disk at this structure solver radius in units of the mean specific angular
       ! momentum of the disk assuming a flat rotation curve.
       !! Determine the specific angular momentum at the size solver radius in units of the mean specific angular
       !! momentum of the disk. This is equal to the ratio of the 1st to 2nd radial moments of the surface density
       !! distribution (assuming a flat rotation curve).
       massDistributionDiskDensityMoment1=massDistributionDisk%surfaceDensityRadialMoment(1.0d0,isInfinite=surfaceDensityMoment1IsInfinite)
       massDistributionDiskDensityMoment2=massDistributionDisk%surfaceDensityRadialMoment(2.0d0,isInfinite=surfaceDensityMoment2IsInfinite)
       if (surfaceDensityMoment1IsInfinite.or.surfaceDensityMoment2IsInfinite) then
          ! One or both of the moments are infinite. Simply assume a value of 0.5 as a default.
          diskStructureSolverSpecificAngularMomentum=0.5d0
       else
          diskStructureSolverSpecificAngularMomentum=  &
               & +radiusStructureSolver            &
               & /(                                    &
               &   +massDistributionDiskDensityMoment2 &
               &   /massDistributionDiskDensityMoment1 &
               &  )
       end if
       ! If necessary, compute the specific angular momentum correction factor to account for the difference between rotation
       ! curves for thin disk and a spherical mass distribution.
       if (structureSolverUseCole2000Method) then
          diskRadiusSolverFlatVsSphericalFactor=                                          &
               & +massDistributionDisk%rotationCurve       (radiusStructureSolver)**2 &
               & *                                          radiusStructureSolver     &
               & -massDistributionDisk%massEnclosedBySphere(radiusStructureSolver)
       end if
    end if
    return
  end subroutine Node_Component_Disk_Standard_Thread_Initialize

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Node_Component_Disk_Standard_Thread_Uninitialize</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Node_Component_Disk_Standard_Thread_Uninitialize()
    !!{
    Uninitializes the standard disk component module for each thread.
    !!}
    use :: Events_Hooks                     , only : postEvolveEvent     , satelliteMergerEvent
    use :: Galacticus_Nodes                 , only : defaultDiskComponent
    use :: Node_Component_Disk_Standard_Data, only : massDistributionDisk, massDistributionDisk_
    implicit none

    if (defaultDiskComponent%standardIsActive()) then
       if (satelliteMergerEvent%isAttached(thread,satelliteMerger)) call satelliteMergerEvent%detach(thread,satelliteMerger)
       if (postEvolveEvent     %isAttached(thread,postEvolve     )) call postEvolveEvent     %detach(thread,postEvolve     )
       !![
       <objectDestructor name="darkMatterHaloScale_"        />
       <objectDestructor name="stellarPopulationProperties_"/>
       <objectDestructor name="starFormationHistory_"       />
       <objectDestructor name="mergerMassMovements_"        />
       <objectDestructor name="galacticStructure_"          />
       <objectDestructor name="massDistributionDisk"        />
       <objectDestructor name="massDistributionDisk_"       />
       !!]
    end if
    return
  end subroutine Node_Component_Disk_Standard_Thread_Uninitialize

  !![
  <calculationResetTask>
    <unitName>Node_Component_Disk_Standard_Calculation_Reset</unitName>
  </calculationResetTask>
  !!]
  subroutine Node_Component_Disk_Standard_Calculation_Reset(node)
    !!{
    Reset standard disk structure calculations.
    !!}
    use :: Galacticus_Nodes                 , only : treeNode
    use :: Node_Component_Disk_Standard_Data, only : Node_Component_Disk_Standard_Reset
    implicit none
    type(treeNode), intent(inout) :: node

    call Node_Component_Disk_Standard_Reset(node%uniqueID())
    return
  end subroutine Node_Component_Disk_Standard_Calculation_Reset

  !![
  <preEvolveTask>
  <unitName>Node_Component_Disk_Standard_Pre_Evolve</unitName>
  </preEvolveTask>
  !!]
  subroutine Node_Component_Disk_Standard_Pre_Evolve(node)
    !!{
    Ensure the disk has been initialized.
    !!}
    use :: Galacticus_Nodes, only : defaultDiskComponent, nodeComponentDisk, nodeComponentDiskStandard, treeNode
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
    !!{
    Trim histories attached to the disk.
    !!}
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

  !![
  <postStepTask>
    <unitName>Node_Component_Disk_Standard_Post_Step</unitName>
  </postStepTask>
  !!]
  subroutine Node_Component_Disk_Standard_Post_Step(node,status)
    !!{
    Trim histories attached to the disk.
    !!}
    use :: Abundances_Structure          , only : abs                 , zeroAbundances
    use :: Display                       , only : displayMessage      , verbosityLevelWarn
    use :: Error                         , only : Error_Report
    use :: Galacticus_Nodes              , only : defaultDiskComponent, nodeComponentDisk  , nodeComponentDiskStandard, nodeComponentSpin, &
          &                                       treeNode            , nodeComponentBasic
    use :: Interface_GSL                 , only : GSL_Success         , GSL_Continue
    use :: ISO_Varying_String            , only : assignment(=)       , operator(//)       , varying_string
    use :: Stellar_Luminosities_Structure, only : abs                 , stellarLuminosities
    use :: String_Handling               , only : operator(//)
    implicit none
    type            (treeNode           ), intent(inout), pointer :: node
    integer                              , intent(inout)          :: status
    class           (nodeComponentDisk  )               , pointer :: disk
    class           (nodeComponentBasic )               , pointer :: basic
    class           (nodeComponentSpin  )               , pointer :: spin
    double precision                     , parameter              :: angularMomentumTolerance=1.0d-2
    double precision                     , save                   :: fractionalErrorMaximum  =0.0d+0
    double precision                                              :: massDisk                       , fractionalError, &
         &                                                           specificAngularMomentum
    character       (len=20             )                         :: valueString
    type            (varying_string     ), save                   :: message
    !$omp threadprivate(message)
    type            (stellarLuminosities), save                   :: luminositiesStellar
    !$omp threadprivate(luminositiesStellar)

    ! Return immediately if this class is not in use.
    if (.not.defaultDiskComponent%standardIsActive()) return
    ! Get the disk component.
    disk => node%disk()
    ! Check if an standard disk component exists.
    select type (disk)
    class is (nodeComponentDiskStandard)
       ! Note that "status" is not set to failure as these changes in state of the disk should not change any calculation of
       ! differential evolution rates as a negative gas/stellar mass was unphysical anyway.
       !
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
             call displayMessage(message,verbosityLevelWarn)
             ! Store the new maximum fractional error.
             fractionalErrorMaximum=fractionalError
          end if
          !$omp end critical (Standard_Disk_Post_Evolve_Check)
          ! Get the specific angular momentum of the disk material
          massDisk= disk%massGas    () &
               &   +disk%massStellar()
          if (massDisk == 0.0d0) then
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
             specificAngularMomentum=disk%angularMomentum()/massDisk
             if (specificAngularMomentum < 0.0d0) specificAngularMomentum=disk%radius()*disk%velocity()
          end if
          ! Reset the gas, abundances and angular momentum of the disk.
          call disk%        massGasSet(                                     0.0d0)
          call disk%  abundancesGasSet(                            zeroAbundances)
          call disk%angularMomentumSet(specificAngularMomentum*disk%massStellar())
          ! Indicate that ODE evolution should continue after this state change.
          if (status == GSL_Success) status=GSL_Continue
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
             call displayMessage(message,verbosityLevelWarn)
             ! Store the new maximum fractional error.
             fractionalErrorMaximum=fractionalError
          end if
          !$omp end critical (Standard_Disk_Post_Evolve_Check)
          ! Get the specific angular momentum of the disk material
          massDisk= disk%massGas    () &
               &   +disk%massStellar()
          if (massDisk == 0.0d0) then
             specificAngularMomentum=0.0d0
             call disk%      massGasSet(         0.0d0)
             call disk%abundancesGasSet(zeroAbundances)
          else
             specificAngularMomentum=disk%angularMomentum()/massDisk
             if (specificAngularMomentum < 0.0d0) specificAngularMomentum=disk%radius()*disk%velocity()
          end if
          ! Reset the stellar, abundances and angular momentum of the disk.
          call disk%      massStellarSet(                                 0.0d0)
          call disk%abundancesStellarSet(                        zeroAbundances)
          call disk%  angularMomentumSet(specificAngularMomentum*disk%massGas())
          ! Indicate that ODE evolution should continue after this state change.
          if (status == GSL_Success) status=GSL_Continue
       end if
       ! Trap negative angular momentum.
       if (disk%angularMomentum() < 0.0d0) then
          spin  => node%spin ()
          basic => node%basic()
          if      (                           &
               &     disk  %massStellar    () &
               &    +disk  %massGas        () &
               &   <=                         &
               &     0.0d0                    &
               &  ) then
             call disk%angularMomentumSet(0.0d0)
          else if (.not.diskNegativeAngularMomentumAllowed) then
             if  (                                &
                  &    abs(disk%angularMomentum())&
                  &   /(                          &
                  &        disk%massStellar    () &
                  &     +  disk%massGas        () &
                  &    )                          &
                  &  <                            &
                  &    angularMomentumTolerance   &
                  &   *spin    %angularMomentum() &
                  &   /basic   %mass           () &
                  & ) then
                call disk%angularMomentumSet(0.0d0)
             else
                message='negative angular momentum in disk with positive mass'
                write (valueString,'(e12.6)') disk  %angularMomentum()
                message=message//char(10)//' -> angular momentum       = '//trim(valueString)
                write (valueString,'(e12.6)') disk  %massStellar    ()
                message=message//char(10)//' -> stellar mass           = '//trim(valueString)
                write (valueString,'(e12.6)') disk  %massGas        ()
                message=message//char(10)//' -> gas mass               = '//trim(valueString)
                write (valueString,'(e12.6)') +spin %angularMomentum() &
                     &                        /basic%mass           ()
                message=message//char(10)//' -> angular momentum scale = '//trim(valueString)
                call Error_Report(message//{introspection:location})
             end if
          end if
          ! Indicate that ODE evolution should continue after this state change.
          if (status == GSL_Success) status=GSL_Continue
       end if
    end select
    return
  end subroutine Node_Component_Disk_Standard_Post_Step

  subroutine Node_Component_Disk_Standard_Create(node)
    !!{
    Create properties in an standard disk component.
    !!}
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

  !![
  <scaleSetTask>
   <unitName>Node_Component_Disk_Standard_Scale_Set</unitName>
  </scaleSetTask>
  !!]
  subroutine Node_Component_Disk_Standard_Scale_Set(node)
    !!{
    Set scales for properties of {\normalfont \ttfamily node}.
    !!}
    use :: Abundances_Structure          , only : abs                 , abundances       , max                      , operator(*)            , &
          &                                       unitAbundances
    use :: Galacticus_Nodes              , only : defaultDiskComponent, nodeComponentDisk, nodeComponentDiskStandard, nodeComponentSpheroid  , &
          &                                       treeNode
    use :: Histories                     , only : history
    use :: Stellar_Luminosities_Structure, only : abs                 , max              , stellarLuminosities      , unitStellarLuminosities
    implicit none
    type            (treeNode                        ), intent(inout), pointer :: node
    class           (nodeComponentDisk               )               , pointer :: disk
    class           (nodeComponentSpheroid           )               , pointer :: spheroid
    double precision                                  , parameter              :: massMinimum                   =1.0d+0
    double precision                                  , parameter              :: angularMomentumMinimum        =1.0d-1
    double precision                                  , parameter              :: fractionTolerance             =1.0d-4
    double precision                                  , parameter              :: luminosityMinimum             =1.0d+0
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
       angularMomentum=+abs(disk    %angularMomentum()) &
            &          +abs(spheroid%angularMomentum())
       call disk%angularMomentumScale(max(angularMomentum,angularMomentumMinimum))
       ! Set scale for masses.
       !! The scale here (and for other quantities below) combines the mass of disk and spheroid. This avoids attempts to solve
       !! tiny disks to high precision in massive spheroidal galaxies.
        mass           =max(                                                      &
            &               +abs(disk%massGas    ())+abs(spheroid%massGas    ())  &
            &               +abs(disk%massStellar())+abs(spheroid%massStellar()), &
            &               +massMinimum                                          &
            &              )
       call disk%massGasScale          (mass)
       call disk%massStellarScale      (mass)
       call disk%massStellarFormedScale(mass)
       ! Set the scale for the retained stellar mass fraction.
       call disk%fractionMassRetainedScale(fractionTolerance*disk%fractionMassRetained())
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
          call disk%abundancesGasScale    (abundancesScale)
          ! Set scale for stellar abundances.
          call disk%abundancesStellarScale(abundancesScale)
       end if
       ! Set scale for stellar luminosities.
       stellarLuminositiesScale=max(                                      &
            &                       +abs(disk    %luminositiesStellar())  &
            &                       +abs(spheroid%luminositiesStellar()), &
            &                           +unitStellarLuminosities          &
            &                           *luminosityMinimum                &
            &                      )
       call stellarLuminositiesScale%truncate                (disk                    %luminositiesStellar())
       call disk                    %luminositiesStellarScale(stellarLuminositiesScale                      )
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

  !![
  <inactiveSetTask>
   <unitName>Node_Component_Disk_Standard_Inactive</unitName>
  </inactiveSetTask>
  !!]
  subroutine Node_Component_Disk_Standard_Inactive(node)
    !!{
    Set Jacobian zero status for properties of {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, nodeComponentDiskStandard, treeNode
    implicit none
    type (treeNode         ), intent(inout), pointer :: node
    class(nodeComponentDisk)               , pointer :: disk

    ! Get the disk component.
    disk => node%disk()
    ! Check if an standard disk component exists.
    select type (disk)
    class is (nodeComponentDiskStandard)
       if (inactiveLuminositiesStellar) call disk%luminositiesStellarInactive()
    end select
    return
  end subroutine Node_Component_Disk_Standard_Inactive

  subroutine satelliteMerger(self,node)
    !!{
    Transfer any standard disk associated with {\normalfont \ttfamily node} to its host halo.
    !!}
    use :: Abundances_Structure            , only : zeroAbundances
    use :: Error                           , only : Error_Report
    use :: Galacticus_Nodes                , only : nodeComponentDisk      , nodeComponentDiskStandard, nodeComponentSpheroid           , treeNode
    use :: Histories                       , only : history
    use :: Satellite_Merging_Mass_Movements, only : destinationMergerDisk  , destinationMergerSpheroid, enumerationDestinationMergerType
    use :: Stellar_Luminosities_Structure  , only : zeroStellarLuminosities
    implicit none
    class           (*                               ), intent(inout) :: self
    type            (treeNode                        ), intent(inout) :: node
    class           (nodeComponentDisk               ), pointer       :: diskHost               , disk
    class           (nodeComponentSpheroid           ), pointer       :: spheroidHost           , spheroid
    type            (treeNode                        ), pointer       :: nodeHost
    type            (history                         )                :: historyHost            , historyNode
    double precision                                                  :: specificAngularMomentum
    type            (enumerationDestinationMergerType)                :: destinationGasSatellite, destinationGasHost       , &
         &                                                               destinationStarsHost   , destinationStarsSatellite
    logical                                                           :: mergerIsMajor
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
       select case (destinationGasSatellite%ID)
       case (destinationMergerDisk%ID)
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
       case (destinationMergerSpheroid%ID)
          call spheroidHost%massGasSet            (                                                                     &
               &                                             spheroidHost%massGas            ()                         &
               &                                            +disk        %massGas            ()                         &
               &                                           )
          call spheroidHost%abundancesGasSet      (                                                                     &
               &                                             spheroidHost%abundancesGas      ()                         &
               &                                            +disk        %abundancesGas      ()                         &
               &                                           )
       case default
          call Error_Report('unrecognized movesTo descriptor'//{introspection:location})
       end select
       call disk%      massGasSet(         0.0d0)
       call disk%abundancesGasSet(zeroAbundances)
       ! Move the stellar component of the standard disk to the host.
       select case (destinationStarsSatellite%ID)
       case (destinationMergerDisk%ID)
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
          call historyNode%destroy                (                                )
          call historyHost%destroy                (                                )
       case (destinationMergerSpheroid%ID)
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
          call historyHost %interpolatedIncrement      (historyNode)
          call historyNode %reset                      (           )
          call spheroidHost%stellarPropertiesHistorySet(historyHost)
          call disk        %stellarPropertiesHistorySet(historyNode)
          ! Also add star formation histories.
          historyNode=disk        %starFormationHistory    ()
          historyHost=spheroidHost%starFormationHistory    ()
          call historyHost %increment              (historyNode,autoExtend  =.true. )
          call historyNode %reset                  (                                )
          call spheroidHost%starFormationHistorySet(historyHost                     )
          call disk        %starFormationHistorySet(historyNode                     )
          call historyNode %destroy                (                                )
          call historyHost %destroy                (                                )
       case default
          call Error_Report('unrecognized movesTo descriptor'//{introspection:location})
       end select
       call disk%        massStellarSet(                  0.0d0)
       call disk%  abundancesStellarSet(         zeroAbundances)
       call disk%luminositiesStellarSet(zeroStellarLuminosities)
       call disk%    angularMomentumSet(                  0.0d0)
    end select
    return
  end subroutine satelliteMerger

  !![
  <radiusSolverPlausibility>
   <unitName>Node_Component_Disk_Standard_Radius_Solver_Plausibility</unitName>
   <after>Node_Component_Basic_Standard_Plausibility</after>
  </radiusSolverPlausibility>
  !!]
  subroutine Node_Component_Disk_Standard_Radius_Solver_Plausibility(node)
    !!{
    Determines whether the disk is physically plausible for radius solving tasks. Require that it have non-zero mass and angular momentum.
    !!}
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
       if      (disk%angularMomentum()                <                   0.0d0) &
            & node%isPhysicallyPlausible=.false.
       if      (disk%massStellar    ()+disk%massGas() <  -toleranceAbsoluteMass) then
          node%isPhysicallyPlausible=.false.
       else if (disk%massStellar    ()+disk%massGas() >=                  0.0d0) then
          if      (                                                              &
               &   disk%angularMomentum() < 0.0d0                                &
               &  ) then
             node%isPhysicallyPlausible=.false.
          else
             angularMomentumScale=(                                           &
                  &                 disk%massStellar()                        &
                  &                +disk%massGas    ()                        &
                  &               )                                           &
                  &               * darkMatterHaloScale_%radiusVirial  (node) &
                  &               * darkMatterHaloScale_%velocityVirial(node)
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
    !!{
    Return the radius of the standard disk used in structure solvers.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, treeNode
    implicit none
    type (treeNode         ), intent(inout) :: node
    class(nodeComponentDisk), pointer       :: disk

    disk => node%disk()
    Node_Component_Disk_Standard_Radius_Solve=disk%radius()*radiusStructureSolver
    return
  end function Node_Component_Disk_Standard_Radius_Solve

  subroutine Node_Component_Disk_Standard_Radius_Solve_Set(node,radius)
    !!{
    Set the radius of the standard disk used in structure solvers.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, treeNode
    implicit none
    type            (treeNode         ), intent(inout) :: node
    double precision                   , intent(in   ) :: radius
    class           (nodeComponentDisk), pointer       :: disk

    disk => node%disk()
    call disk%radiusSet(max(radius,0.0d0)/radiusStructureSolver)
    return
  end subroutine Node_Component_Disk_Standard_Radius_Solve_Set

  double precision function Node_Component_Disk_Standard_Velocity(node)
    !!{
    Return the circular velocity of the standard disk.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, treeNode
    implicit none
    type (treeNode         ), intent(inout) :: node
    class(nodeComponentDisk), pointer       :: disk

    disk => node%disk()
    Node_Component_Disk_Standard_Velocity=disk%velocity()
    return
  end function Node_Component_Disk_Standard_Velocity

  subroutine Node_Component_Disk_Standard_Velocity_Set(node,velocity)
    !!{
    Set the circular velocity of the standard disk.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, treeNode
    implicit none
    type            (treeNode         ), intent(inout) :: node
    double precision                   , intent(in   ) :: velocity
    class           (nodeComponentDisk), pointer       :: disk

    disk => node%disk()
    call disk%velocitySet(velocity)
    return
  end subroutine Node_Component_Disk_Standard_Velocity_Set

  !![
  <radiusSolverTask>
   <unitName>Node_Component_Disk_Standard_Radius_Solver</unitName>
  </radiusSolverTask>
  !!]
  subroutine Node_Component_Disk_Standard_Radius_Solver(node,componentActive,specificAngularMomentumRequired,specificAngularMomentum,Radius_Get,Radius_Set,Velocity_Get&
       &,Velocity_Set)
    !!{
    Interface for the size solver algorithm.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentDisk              , nodeComponentDiskStandard, treeNode
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    type            (treeNode                                     ), intent(inout)          :: node
    logical                                                        , intent(  out)          :: componentActive
    logical                                                        , intent(in   )          :: specificAngularMomentumRequired
    double precision                                               , intent(  out)          :: specificAngularMomentum
    procedure       (Node_Component_Disk_Standard_Radius_Solve    ), intent(  out), pointer :: Radius_Get                     , Velocity_Get
    procedure       (Node_Component_Disk_Standard_Radius_Solve_Set), intent(  out), pointer :: Radius_Set                     , Velocity_Set
    class           (nodeComponentDisk                            )               , pointer :: disk
    double precision                                                                        :: angularMomentum                , massDisk    , &
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
             massDisk= disk%massGas    () &
                  &   +disk%massStellar()
             if (massDisk > 0.0d0) then
                specificAngularMomentumMean=angularMomentum/massDisk
             else
                specificAngularMomentumMean=0.0d0
             end if
             specificAngularMomentum=specificAngularMomentumMean*diskStructureSolverSpecificAngularMomentum
             ! If using the Cole et al. (2000) method for disk radii, adjust the specific angular momentum to account for the
             ! difference between rotation curves for thin disk and a spherical mass distribution. Trap instances where this leads to
             ! imaginary specific angular momentum - this can happen as the radius solver explores the allowed range of radii when
             ! seeking a solution.
             if (structureSolverUseCole2000Method)                                                       &
                  & specificAngularMomentum=sqrt(                                                      &
                  &                               max(                                                 &
                  &                                    0.0d0,                                          &
                  &                                    specificAngularMomentum**2                      &
                  &                                   -diskRadiusSolverFlatVsSphericalFactor           &
                  &                                   *gravitationalConstantGalacticus                 &
                  &                                   *massDisk                                        &
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

  !![
  <mergerTreeExtraOutputTask>
   <unitName>Node_Component_Disk_Standard_Star_Formation_History_Output</unitName>
  </mergerTreeExtraOutputTask>
  !!]
  subroutine Node_Component_Disk_Standard_Star_Formation_History_Output(node,iOutput,treeIndex,nodePassesFilter,treeLock)
    !!{
    Store the star formation history in the output file.
    !!}
    use            :: Galacticus_Nodes          , only : defaultDiskComponent, nodeComponentDisk, nodeComponentDiskStandard, treeNode
    use            :: Galactic_Structure_Options, only : componentTypeDisk
    use            :: Histories                 , only : history
    use, intrinsic :: ISO_C_Binding             , only : c_size_t
    use            :: Kind_Numbers              , only : kind_int8
    use            :: Locks                     , only : ompLock
    implicit none
    type   (treeNode         ), intent(inout), pointer :: node
    integer(c_size_t         ), intent(in   )          :: iOutput
    integer(kind=kind_int8   ), intent(in   )          :: treeIndex
    logical                   , intent(in   )          :: nodePassesFilter
    type   (ompLock          ), intent(inout)          :: treeLock
    class  (nodeComponentDisk)               , pointer :: disk
    type   (history          )                         :: historyStarFormation

    ! Check if we are the default method.
    if (.not.defaultDiskComponent%standardIsActive()) return
    ! Output the star formation history if a disk exists for this component.
    disk                 => node%disk                ()
    historyStarFormation =  disk%starFormationHistory()
    call starFormationHistory_%output(node,nodePassesFilter,historyStarFormation,iOutput,treeIndex,componentTypeDisk,treeLock)
    ! Update the star formation history only if a disk exists.
    select type (disk)
    class is (nodeComponentDiskStandard)
       call disk%starFormationHistorySet(historyStarFormation)
    end select
    return
  end subroutine Node_Component_Disk_Standard_Star_Formation_History_Output

  !![
  <mergerTreeExtraOutputFlush>
   <unitName>Node_Component_Disk_Standard_Star_Formation_History_Flush</unitName>
  </mergerTreeExtraOutputFlush>
  !!]
  subroutine Node_Component_Disk_Standard_Star_Formation_History_Flush(treeLock)
    !!{
    Flush star formation history data.
    !!}
    use :: Galacticus_Nodes          , only : defaultDiskComponent
    use :: Galactic_Structure_Options, only : componentTypeDisk
    use :: Locks                     , only : ompLock
    implicit none
    type(ompLock), intent(inout) :: treeLock

    ! Check if we are the default method.
    if (.not.defaultDiskComponent%standardIsActive()) return
    ! Flush the star formation history.
    call starFormationHistory_%outputFlush(componentTypeDisk,treeLock)
    return
  end subroutine Node_Component_Disk_Standard_Star_Formation_History_Flush

  !![
  <stateStoreTask>
   <unitName>Node_Component_Disk_Standard_State_Store</unitName>
  </stateStoreTask>
  !!]
  subroutine Node_Component_Disk_Standard_State_Store(stateFile,gslStateFile,stateOperationID)
    !!{
    Write the tabulation state to file.
    !!}
    use            :: Display                          , only : displayMessage      , verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding                    , only : c_ptr               , c_size_t
    use            :: Node_Component_Disk_Standard_Data, only : massDistributionDisk
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Storing state for: componentDisk -> standard',verbosity=verbosityLevelInfo)
    !![
    <stateStore variables="massDistributionDisk darkMatterHaloScale_ stellarPopulationProperties_ starFormationHistory_ mergerMassMovements_"/>
    !!]
    write (stateFile) diskStructureSolverSpecificAngularMomentum,diskRadiusSolverFlatVsSphericalFactor
    return
  end subroutine Node_Component_Disk_Standard_State_Store

  !![
  <stateRetrieveTask>
   <unitName>Node_Component_Disk_Standard_State_Retrieve</unitName>
  </stateRetrieveTask>
  !!]
  subroutine Node_Component_Disk_Standard_State_Retrieve(stateFile,gslStateFile,stateOperationID)
    !!{
    Retrieve the tabulation state from the file.
    !!}
    use            :: Display                          , only : displayMessage      , verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding                    , only : c_ptr               , c_size_t
    use            :: Node_Component_Disk_Standard_Data, only : massDistributionDisk
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Retrieving state for: componentDisk -> standard',verbosity=verbosityLevelInfo)
    !![
    <stateRestore variables="massDistributionDisk darkMatterHaloScale_ stellarPopulationProperties_ starFormationHistory_ mergerMassMovements_"/>
    !!]
    read (stateFile) diskStructureSolverSpecificAngularMomentum,diskRadiusSolverFlatVsSphericalFactor
    return
  end subroutine Node_Component_Disk_Standard_State_Retrieve
  
  function Node_Component_Disk_Standard_Chandrasekhar_Integral(self,nodeSatellite,positionCartesian,velocityCartesian,componentType,massType)
    !!{
    Computes the gravitational acceleration at a given position for a standard disk.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentDiskStandard        , treeNode
    use :: Galactic_Structure_Options      , only : componentTypeAll                 , componentTypeDisk           , massTypeAll            , weightByMass         , &
         &                                          weightIndexNull                  , enumerationComponentTypeType, enumerationMassTypeType
    use :: Numerical_Constants_Math        , only : Pi
    use :: Coordinates                     , only : assignment(=)                    , coordinateSpherical         , coordinateCartesian    , coordinateCylindrical
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Mass_Distributions              , only : massDistributionGaussianEllipsoid
    use :: Linear_Algebra                  , only : vector                           , matrix             , assignment(=)
    implicit none
    double precision                                                  , dimension(3) :: Node_Component_Disk_Standard_Chandrasekhar_Integral
    class           (nodeComponentDiskStandard        ), intent(inout)               :: self
    type            (treeNode                         ), intent(inout)               :: nodeSatellite
    double precision                                   , intent(in   ), dimension(3) :: positionCartesian                                          , velocityCartesian
    type            (enumerationComponentTypeType     ), intent(in   )               :: componentType
    type            (enumerationMassTypeType          ), intent(in   )               :: massType
    double precision                                   , parameter                   :: toomreQRadiusHalfMass                              =1.50d0  ! The Toomre Q-parameter at the disk half-mass radius (Benson et al.,
    ! 2004 , https://ui.adsabs.harvard.edu/abs/2004MNRAS.351.1215B, Appendix A).
    double precision                                   , parameter                   :: toomreQFactor                                      =3.36d0  ! The factor appearing in the definition of the Toomre Q-parameter for
    ! a stellar disk (Binney & Tremaine, eqn. 6.71).
    double precision                                                  , dimension(3) :: velocityDisk                                               , velocityRelative                , &
         &                                                                              positionSpherical                                          , positionSphericalMidplane       , &
         &                                                                              positionCartesianMidplane                                  , positionCylindricalMidplane     , &
         &                                                                              positionCylindricalHalfMass
    type            (massDistributionGaussianEllipsoid), save                        :: velocityDistribution
    logical                                            , save                        :: velocityDistributionInitialized                    =.false.
    !$omp threadprivate(velocityDistribution,velocityDistributionInitialized)
    type            (coordinateSpherical              )                              :: coordinatesSpherical
    type            (coordinateCartesian              )                              :: coordinatesCartesian
    type            (coordinateCylindrical            )                              :: coordinatesCylindrical
    double precision                                                                 :: velocityDispersionRadial                                   , velocityDispersionAzimuthal     , &
         &                                                                              velocityDispersionVertical                                 , velocityCircular                , &
         &                                                                              velocityCircularHalfMassRadius                             , velocityCircularSquaredGradient , &
         &                                                                              velocityCircularSquaredGradientHalfMassRadius              , density                         , &
         &                                                                              densityMidPlane                                            , densitySurface                  , &
         &                                                                              heightScale                                                , radiusMidplane                  , &
         &                                                                              frequencyCircular                                          , frequencyEpicyclic              , &
         &                                                                              frequencyCircularHalfMassRadius                            , frequencyEpicyclicHalfMassRadius, &
         &                                                                              densitySurfaceRadiusHalfMass                               , velocityDispersionRadialHalfMass, &
         &                                                                              velocityDispersionMaximum                                  , velocityRelativeMagnitude       , &
         &                                                                              factorSuppressionExtendedMass                              , radiusHalfMass
    type            (matrix                           )                              :: rotation

    ! Return if the disk component is not selected.
    Node_Component_Disk_Standard_Chandrasekhar_Integral=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDisk) .or. self%radius() <= 0.0d0) return
    ! Construct the velocity vector of the disk rotation.
    positionCartesianMidplane                    =[positionCartesian(1),positionCartesian(2),0.0d0]
    coordinatesCartesian                         = positionCartesian
    coordinatesSpherical                         = coordinatesCartesian
    positionSpherical                            = coordinatesSpherical
    coordinatesCartesian                         = positionCartesianMidplane
    coordinatesSpherical                         = coordinatesCartesian
    coordinatesCylindrical                       = coordinatesCartesian 
    positionSphericalMidplane                    = coordinatesSpherical
    positionCylindricalMidplane                  = coordinatesCylindrical
    positionCylindricalHalfMass                  =[self%halfMassRadius(),0.0d0,0.0d0]
    radiusMidplane                               = coordinatesCylindrical%r()
    velocityCircular                             =self%rotationCurve        (     radiusMidplane  ,componentType,massType)
    velocityCircularSquaredGradient              =self%rotationCurveGradient(     radiusMidplane  ,componentType,massType)
    velocityCircularHalfMassRadius               =self%rotationCurve        (self%halfMassRadius(),componentType,massType)
    velocityCircularSquaredGradientHalfMassRadius=self%rotationCurveGradient(self%halfMassRadius(),componentType,massType)
    velocityDisk                                =+[positionCartesianMidplane(2),-positionCartesianMidplane(1),0.0d0] &
         &                                        /radiusMidplane                                                     &
         &                                        *velocityCircular
    ! Compute epicyclic frequency.
    frequencyCircular               =velocityCircular              /     radiusMidplane
    frequencyCircularHalfMassRadius =velocityCircularHalfMassRadius/self%halfMassRadius()
    frequencyEpicyclic              =sqrt(velocityCircularSquaredGradient              /     radiusMidplane  +2.0d0*frequencyCircular              **2)
    frequencyEpicyclicHalfMassRadius=sqrt(velocityCircularSquaredGradientHalfMassRadius/self%halfMassRadius()+2.0d0*frequencyCircularHalfMassRadius**2)
    ! Get disk structural properties.
    density                     =+self%density       (positionSpherical          ,componentTypeDisk,massTypeAll,weightByMass,weightIndexNull)
    densityMidPlane             =+self%density       (positionSphericalMidplane  ,componentTypeDisk,massTypeAll,weightByMass,weightIndexNull)
    densitySurface              =+self%surfaceDensity(positionCylindricalMidplane,componentTypeDisk,massTypeAll,weightByMass,weightIndexNull)
    densitySurfaceRadiusHalfMass=+self%surfaceDensity(positionCylindricalHalfMass,componentTypeDisk,massTypeAll,weightByMass,weightIndexNull)
    if (density <= 0.0d0) return
    heightScale                 =+0.5d0           &
         &                       *densitySurface  &
         &                       /densityMidPlane
    ! Compute normalization of the radial velocity dispersion.
    velocityDispersionRadialHalfMass=+toomreQFactor                    &
         &                           *gravitationalConstantGalacticus  &
         &                           *densitySurfaceRadiusHalfMass     &
         &                           *toomreQRadiusHalfMass            &
         &                           /frequencyEpicyclicHalfMassRadius
    ! Find the velocity dispersion components of the disk.
    velocityDispersionRadial   =+velocityDispersionRadialHalfMass                &
         &                      *exp(-     radiusMidPlane  /self%radius()/2.0d0) &
         &                      /exp(-self%halfMassRadius()/self%radius()/2.0d0)
    velocityDispersionAzimuthal=+velocityDispersionRadial*frequencyEpicyclic/2.0d0/frequencyCircular
    velocityDispersionVertical =+sqrt(Pi*gravitationalConstantGalacticus*densitySurface*heightScale)
    velocityDispersionMaximum  =+maxval([velocityDispersionRadial,velocityDispersionAzimuthal,velocityDispersionVertical])
    velocityDispersionRadial   =+velocityDispersionRadial   /velocityDispersionMaximum
    velocityDispersionAzimuthal=+velocityDispersionAzimuthal/velocityDispersionMaximum
    velocityDispersionVertical =+velocityDispersionVertical /velocityDispersionMaximum
    if (any([velocityDispersionRadial,velocityDispersionAzimuthal,velocityDispersionVertical] <= 0.0d0)) return
    ! Find the relative velocity of the perturber and the disk.
    velocityRelative=(velocityCartesian-velocityDisk)/velocityDispersionMaximum
    ! Handle limiting case of large relative velocity.
    velocityRelativeMagnitude=sqrt(sum(velocityRelative**2))
    ! Initialize the velocity distribution.
    rotation=reshape(                                                                               &
         &            [                                                                             &
         &             +positionCartesianMidplane(1),-positionCartesianMidplane(2),+0.0d0         , &
         &             +positionCartesianMidplane(2),+positionCartesianMidplane(1),+0.0d0         , &
         &             +0.0d0                       ,+0.0d0                       ,+radiusMidplane  &
         &            ]                                                                             &
         &           /radiusMidplane                                                              , &
         &           [3,3]                                                                          &
         &          )
    coordinatesCartesian=velocityRelative
    if (.not.velocityDistributionInitialized) then
       velocityDistribution           =massDistributionGaussianEllipsoid(scaleLength=[1.0d0,1.0d0,1.0d0],rotation=rotation,mass=1.0d0,dimensionless=.true.)
       velocityDistributionInitialized=.true.
    end if
    call velocityDistribution%initialize(scaleLength=[velocityDispersionRadial,velocityDispersionAzimuthal,velocityDispersionVertical],rotation=rotation)
    ! Compute suppression factor due to satellite being an extended mass distribution. This is largely untested - it is meant to
    ! simply avoid extremely large accelerations for subhalo close to the disk plane when that subhalo is much more extended than
    ! the disk.
    radiusHalfMass=galacticStructure_%radiusEnclosingMass(                                 &
         &                                                               nodeSatellite   , &
         &                                                massFractional=0.5d0           , &
         &                                                componentType =componentTypeAll, &
         &                                                massType      =massTypeAll       &
         &                                               )
    if (radiusHalfMass > heightScale) then
       factorSuppressionExtendedMass=+heightScale    &
            &                        /radiusHalfMass
    else
       factorSuppressionExtendedMass=+1.0d0
    end if
    ! Evaluate the integral.
    Node_Component_Disk_Standard_Chandrasekhar_Integral=+density                                                             &
         &                                              *velocityDistribution         %acceleration(coordinatesCartesian)    &
         &                                              /velocityDispersionMaximum                                       **2 &
         &                                              *factorSuppressionExtendedMass
    return
  end function Node_Component_Disk_Standard_Chandrasekhar_Integral

end module Node_Component_Disk_Standard
