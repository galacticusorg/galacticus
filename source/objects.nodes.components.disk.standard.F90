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
  implicit none
  private
  public :: Node_Component_Disk_Standard_Scale_Set                 , Node_Component_Disk_Standard_Pre_Evolve                  , &
       &    Node_Component_Disk_Standard_Radius_Solver_Plausibility, Node_Component_Disk_Standard_Radius_Solver               , &
       &    Node_Component_Disk_Standard_Post_Step                 , Node_Component_Disk_Standard_Thread_Uninitialize         , &
       &    Node_Component_Disk_Standard_Initialize                , Node_Component_Disk_Standard_Thread_Initialize           , &
       &    Node_Component_Disk_Standard_State_Store               , Node_Component_Disk_Standard_State_Retrieve              , &
       &    Node_Component_Disk_Standard_Inactive

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
      <attributes isSettable="true" isGettable="true" isEvolvable="true" isNonNegative="true" />
      <output unitsInSI="massSolar" comment="Mass of stars in the standard disk."/>
    </property>
    <property>
      <name>massStellarFormed</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" isNonNegative="true" />
    </property>
    <property>
      <name>fractionMassRetained</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" isNonNegative="true" />
    </property>
    <property>
      <name>abundancesStellar</name>
      <type>abundances</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" isNonNegative="true" />
      <output unitsInSI="massSolar" comment="Mass of metals in the stellar phase of the standard disk."/>
    </property>
    <property>
      <name>massGas</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" isNonNegative="true" />
      <output unitsInSI="massSolar" comment="Mass of gas in the standard disk."/>
    </property>
    <property>
      <name>abundancesGas</name>
      <type>abundances</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" isNonNegative="true" />
      <output unitsInSI="massSolar" comment="Mass of metals in the gas phase of the standard disk."/>
    </property>
    <property>
      <name>angularMomentum</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" isNonNegative="true" />
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
      <attributes isSettable="true" isGettable="true" isEvolvable="true" isNonNegative="true" />
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
      <attributes isSettable="true" isGettable="true" isEvolvable="true" isNonNegative="true" />
    </property>
   </properties>
   <bindings>
    <binding method="massDistribution" function="Node_Component_Disk_Standard_Mass_Distribution"/>
    <binding method="massBaryonic"     function="Node_Component_Disk_Standard_Mass_Baryonic"    />
   </bindings>
   <functions>objects.nodes.components.disk.standard.bound_functions.inc</functions>
  </component>
  !!]

  ! Objects used by this component.
  class(darkMatterHaloScaleClass        ), pointer :: darkMatterHaloScale_
  class(stellarPopulationPropertiesClass), pointer :: stellarPopulationProperties_
  class(starFormationHistoryClass       ), pointer :: starFormationHistory_
  class(mergerMassMovementsClass        ), pointer :: mergerMassMovements_
  !$omp threadprivate(darkMatterHaloScale_,stellarPopulationProperties_,starFormationHistory_,mergerMassMovements_)

  ! Internal count of abundances.
  integer                                     :: abundancesCount

  ! Parameters controlling the physical implementation.
  double precision                            :: toleranceAbsoluteMass                              , radiusStructureSolver               , &
       &                                         toleranceRelativeMetallicity
  logical                                     :: diskNegativeAngularMomentumAllowed                 , structureSolverUseCole2000Method    , &
       &                                         inactiveLuminositiesStellar                        , postStepZeroNegativeMasses

  ! History of trial radii used to check for oscillations in the solution when solving for the structure of the disk.
  integer                                     :: radiusSolverIteration
  double precision, dimension(2)              :: radiusHistory
  !$omp threadprivate(radiusHistory,radiusSolverIteration)

  ! The largest and smallest angular momentum, in units of that of a circular orbit at the virial radius, considered to be
  ! physically plausible for a disk.
  double precision, parameter                 :: angularMomentumMaximum                    =1.0d+1
  double precision, parameter                 :: angularMomentumMinimum                    =1.0d-6

  ! Disk structural parameters.
  double precision                            :: ratioAngularMomentumSolverRadius        , diskRadiusSolverFlatVsSphericalFactor
  !$omp threadprivate(ratioAngularMomentumSolverRadius,diskRadiusSolverFlatVsSphericalFactor)

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
    type(inputParameters), intent(inout) :: parameters
    type(inputParameters)                :: subParameters

    if (defaultDiskComponent%standardIsActive()) then
       ! Get number of abundance properties.
       abundancesCount  =Abundances_Property_Count            ()
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
         <description>If true, use the method described in \cite{cole_hierarchical_2000} to correct for difference between thin disk and spherical mass distributions when solving for disk radii.</description>
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
       <inputParameter>
         <name>postStepZeroNegativeMasses</name>
         <defaultValue>.true.</defaultValue>
         <description>If true, negative masses will be zeroed after each ODE step. Note that this can lead to non-conservation of mass.</description>
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
    use :: Events_Hooks                     , only : dependencyDirectionAfter   , dependencyRegEx            , openMPThreadBindingAtLevel, postEvolveEvent, &
          &                                          satelliteMergerEvent       , mergerTreeExtraOutputEvent
    use :: Error                            , only : Error_Report
    use :: Galacticus_Nodes                 , only : defaultDiskComponent
    use :: Input_Parameters                 , only : inputParameter             , inputParameters
    use :: Node_Component_Disk_Standard_Data, only : massDistributionStellar_   , massDistributionGas_       , kinematicDistribution_
    use :: Mass_Distributions               , only : massDistributionCylindrical, kinematicsDistributionLocal
    use :: Galactic_Structure_Options       , only : componentTypeDisk          , massTypeStellar            , massTypeGaseous
    implicit none
    type            (inputParameters), intent(inout) :: parameters
    type            (dependencyRegEx), dimension(2)  :: dependencies
    double precision                                 :: massDistributionDiskDensityMoment1     , massDistributionDiskDensityMoment2, &
         &                                              ratioAngularMomentumSolverRadiusDefault
    logical                                          :: surfaceDensityMoment1IsInfinite        , surfaceDensityMoment2IsInfinite
    type            (inputParameters)                :: subParameters

    ! Check if this implementation is selected. If so, initialize the mass distribution.
    if (defaultDiskComponent%standardIsActive()) then
       dependencies(1)=dependencyRegEx(dependencyDirectionAfter,'^remnantStructure:')
       dependencies(2)=dependencyRegEx(dependencyDirectionAfter,'^preAnalysis:'     )
       call satelliteMergerEvent      %attach(thread,satelliteMerger      ,openMPThreadBindingAtLevel,label='nodeComponentDiskStandard',dependencies=dependencies)
       call postEvolveEvent           %attach(thread,postEvolve           ,openMPThreadBindingAtLevel,label='nodeComponentDiskStandard'                          )
       call mergerTreeExtraOutputEvent%attach(thread,mergerTreeExtraOutput,openMPThreadBindingAtLevel,label='nodeComponentDiskStandard'                          )
       ! Find our parameters.
       subParameters=parameters%subParameters('componentDisk')
       !![
       <objectBuilder class="darkMatterHaloScale"                                              name="darkMatterHaloScale_"         source="subParameters"                    />
       <objectBuilder class="stellarPopulationProperties"                                      name="stellarPopulationProperties_" source="subParameters"                    />
       <objectBuilder class="starFormationHistory"                                             name="starFormationHistory_"        source="subParameters"                    />
       <objectBuilder class="mergerMassMovements"                                              name="mergerMassMovements_"         source="subParameters"                    />
       <objectBuilder class="massDistribution"            parameterName="massDistributionDisk" name="massDistributionStellar_"     source="subParameters" threadPrivate="yes" >
        <default>
         <massDistributionDisk value="exponentialDisk">
          <dimensionless value="true"/>
         </massDistributionDisk>
        </default>
       </objectBuilder>
       !!]
       ! Validate the disk mass distribution.
       select type (massDistributionStellar_)
       class is (massDistributionCylindrical)
          ! The disk mass distribution must have cylindrical symmetry. So, this is acceptable.        
       class default
          call Error_Report('only cylindrically symmetric mass distributions are allowed'//{introspection:location})
       end select
       if (.not.massDistributionStellar_%isDimensionless()) call Error_Report('disk mass distribution must be dimensionless'//{introspection:location})
       ! Duplicate the dimensionless mass distribution to use for the gas component, and set component and mass types in both.
       !$omp critical(diskStandardDeepCopy)
       allocate(massDistributionGas_,mold=massDistributionStellar_)
       !![
       <deepCopyReset variables="massDistributionStellar_"/>
       <deepCopy source="massDistributionStellar_" destination="massDistributionGas_"/>
       <deepCopyFinalize variables="massDistributionGas_"/>  
       !!]
       !$omp end critical(diskStandardDeepCopy)
       call massDistributionStellar_%setTypes(componentTypeDisk,massTypeStellar)
       call massDistributionGas_    %setTypes(componentTypeDisk,massTypeGaseous)
       ! Construct the kinematic distribution.
       allocate(kinematicDistribution_)
       !![
       <referenceConstruct object="kinematicDistribution_" constructor="kinematicsDistributionLocal(alpha=1.0d0/sqrt(2.0d0))"/>
       !!]
       ! Compute the specific angular momentum of the disk at this structure solver radius in units of the mean specific angular
       ! momentum of the disk assuming a flat rotation curve.
       !! Determine the specific angular momentum at the size solver radius in units of the mean specific angular
       !! momentum of the disk. This is equal to the ratio of the 1st to 2nd radial moments of the surface density
       !! distribution (assuming a flat rotation curve).
       massDistributionDiskDensityMoment1=massDistributionStellar_%surfaceDensityRadialMoment(1.0d0,isInfinite=surfaceDensityMoment1IsInfinite)
       massDistributionDiskDensityMoment2=massDistributionStellar_%surfaceDensityRadialMoment(2.0d0,isInfinite=surfaceDensityMoment2IsInfinite)
       if (surfaceDensityMoment1IsInfinite.or.surfaceDensityMoment2IsInfinite) then
          ! One or both of the moments are infinite. Simply assume a value of 0.5 as a default.
          ratioAngularMomentumSolverRadiusDefault=+0.5d0
       else
          ratioAngularMomentumSolverRadiusDefault=+radiusStructureSolver                &
               &                                  /(                                    &
               &                                    +massDistributionDiskDensityMoment2 &
               &                                    /massDistributionDiskDensityMoment1 &
               &                                   )
       end if
       !$omp critical (diskStandardInitializeAngularMomentum)
       !![
       <inputParameter>
         <name>ratioAngularMomentumSolverRadius</name>
         <defaultSource>($I_1/I_2$ where $I_n=\int_0^\infty \Sigma(R) R^n \mathrm{d}R$, where $\Sigma(R)$ is the disk surface density profile, unless either $I_1$ or $I_2$ is infinite, in which case a default of $1/2$ is used instead.)</defaultSource>
         <defaultValue>ratioAngularMomentumSolverRadiusDefault</defaultValue>
         <description>The assumed ratio of the specific angular momentum at the structure solver radius to the mean specific angular momentum of the standard disk component.</description>
         <source>subParameters</source>
       </inputParameter>
       !!]
       !$omp end critical (diskStandardInitializeAngularMomentum)
       ! If necessary, compute the specific angular momentum correction factor to account for the difference between rotation
       ! curves for thin disk and a spherical mass distribution.
       if (structureSolverUseCole2000Method) then
          diskRadiusSolverFlatVsSphericalFactor=                                          &
               & +massDistributionStellar_%rotationCurve       (radiusStructureSolver)**2 &
               & *                                              radiusStructureSolver     &
               & -massDistributionStellar_%massEnclosedBySphere(radiusStructureSolver)
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
    use :: Events_Hooks                     , only : postEvolveEvent         , satelliteMergerEvent, mergerTreeExtraOutputEvent
    use :: Galacticus_Nodes                 , only : defaultDiskComponent
    use :: Node_Component_Disk_Standard_Data, only : massDistributionStellar_, massDistributionGas_, kinematicDistribution_
    implicit none

    if (defaultDiskComponent%standardIsActive()) then
       if (satelliteMergerEvent      %isAttached(thread,satelliteMerger      )) call satelliteMergerEvent      %detach(thread,satelliteMerger      )
       if (postEvolveEvent           %isAttached(thread,postEvolve           )) call postEvolveEvent           %detach(thread,postEvolve           )
       if (mergerTreeExtraOutputEvent%isAttached(thread,mergerTreeExtraOutput)) call mergerTreeExtraOutputEvent%detach(thread,mergerTreeExtraOutput)
       !![
       <objectDestructor name="darkMatterHaloScale_"        />
       <objectDestructor name="stellarPopulationProperties_"/>
       <objectDestructor name="starFormationHistory_"       />
       <objectDestructor name="mergerMassMovements_"        />
       <objectDestructor name="massDistributionStellar_"    />
       <objectDestructor name="massDistributionGas_"        />
       <objectDestructor name="kinematicDistribution_"      />
       !!]
    end if
    return
  end subroutine Node_Component_Disk_Standard_Thread_Uninitialize

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
    use :: Histories                     , only : history
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
    type            (history            ), save                   :: historyStellar
    !$omp threadprivate(historyStellar)
    
    ! Return immediately if this class is not in use or if masses are not to be zeroed.
    if (.not.defaultDiskComponent%standardIsActive().or..not.postStepZeroNegativeMasses) return
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
             historyStellar=disk%starFormationHistory    ()
             call historyStellar%reset()
             call disk         %starFormationHistorySet     (historyStellar)
             historyStellar=disk%stellarPropertiesHistory()
             call historyStellar%reset()
             call disk          %stellarPropertiesHistorySet(historyStellar)
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
          ! Reset the stellar mass, abundances and angular momentum of the disk.
          call disk%      massStellarSet(                                 0.0d0)
          call disk%abundancesStellarSet(                        zeroAbundances)
          call disk%  angularMomentumSet(specificAngularMomentum*disk%massGas())
          historyStellar=disk%starFormationHistory    ()
          call historyStellar%reset()
          call disk         %starFormationHistorySet     (historyStellar)
          historyStellar=disk%stellarPropertiesHistory()
          call historyStellar%reset()
          call disk          %stellarPropertiesHistorySet(historyStellar)
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
          basic    => node %basic()
          timeBegin=  basic%time ()
       end if
       call starFormationHistory_%create                 (node,historyStarFormation,timeBegin)
       call disk                 %starFormationHistorySet(     historyStarFormation          )
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
            &              +abs(disk%massGas    ())+abs(spheroid%massGas    ())  &
            &              +abs(disk%massStellar())+abs(spheroid%massStellar()), &
            &              +massMinimum                                          &
            &             )
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
       call starFormationHistory_%scales        (stellarPopulationHistoryScales,node,disk%massStellar(),disk%massGas(),disk%abundancesStellar())
       call disk%starFormationHistoryScale      (stellarPopulationHistoryScales                                                 )
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
    use :: Kind_NUmbers, only : kind_int8
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
          call starFormationHistory_%move                   (nodeHost,node,historyHost,historyNode)
          call diskHost             %starFormationHistorySet(              historyHost            )
          call disk                 %starFormationHistorySet(                          historyNode)
          call historyNode          %destroy                (                                     )
          call historyHost          %destroy                (                                     )
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
          historyNode=disk        %stellarPropertiesHistory()
          historyHost=spheroidHost%stellarPropertiesHistory()
          call historyHost %interpolatedIncrement      (historyNode)
          call historyNode %reset                      (           )
          call spheroidHost%stellarPropertiesHistorySet(historyHost)
          call disk        %stellarPropertiesHistorySet(historyNode)
          ! Also add star formation histories.
          historyNode=disk        %starFormationHistory    ()
          historyHost=spheroidHost%starFormationHistory    ()
          call starFormationHistory_%move                   (nodeHost,node,historyHost,historyNode)
          call spheroidHost         %starFormationHistorySet(              historyHost            )
          call disk                 %starFormationHistorySet(                          historyNode)
          call historyNode          %destroy                (                                     )
          call historyHost          %destroy                (                                     )
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
       if (disk%massStellar()+disk%massGas() >= 0.0d0 .and. disk%angularMomentum() > 0.0d0) then
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
  subroutine Node_Component_Disk_Standard_Radius_Solver(node,componentActive,component,specificAngularMomentumRequired,specificAngularMomentum,Radius_Get,Radius_Set,Velocity_Get&
       &,Velocity_Set)
    !!{
    Interface for the size solver algorithm.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentDisk             , nodeComponentDiskStandard, treeNode
    use :: Galactic_Structure_Options      , only : enumerationComponentTypeType  , componentTypeDisk
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    type            (treeNode                                     ), intent(inout)          :: node
    logical                                                        , intent(  out)          :: componentActive
    type            (enumerationComponentTypeType                 ), intent(  out)          :: component
    logical                                                        , intent(in   )          :: specificAngularMomentumRequired
    double precision                                               , intent(  out)          :: specificAngularMomentum
    procedure       (Node_Component_Disk_Standard_Radius_Solve    ), intent(  out), pointer :: Radius_Get                     , Velocity_Get
    procedure       (Node_Component_Disk_Standard_Radius_Solve_Set), intent(  out), pointer :: Radius_Set                     , Velocity_Set
    class           (nodeComponentDisk                            )               , pointer :: disk
    double precision                                                                        :: angularMomentum                , massDisk    , &
         &                                                                                     specificAngularMomentumMean

    ! Determine if node has an active disk component supported by this module.
    componentActive         =  .false.
    component               =  componentTypeDisk
    specificAngularMomentum =  0.0d0
    disk                    => node%disk()
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
             specificAngularMomentum=specificAngularMomentumMean*ratioAngularMomentumSolverRadius
             ! If using the Cole et al. (2000) method for disk radii, adjust the specific angular momentum to account for the
             ! difference between rotation curves for thin disk and a spherical mass distribution. Trap instances where this leads to
             ! imaginary specific angular momentum - this can happen as the radius solver explores the allowed range of radii when
             ! seeking a solution.
             if (structureSolverUseCole2000Method)                                                     &
                  & specificAngularMomentum=sqrt(                                                      &
                  &                               max(                                                 &
                  &                                    0.0d0,                                          &
                  &                                    specificAngularMomentum**2                      &
                  &                                   -diskRadiusSolverFlatVsSphericalFactor           &
                  &                                   *gravitationalConstant_internal                  &
                  &                                   *massDisk                                        &
                  &                                   *Node_Component_Disk_Standard_Radius_Solve(node) &
                  &                                  )                                                 &
                  )
          end if
       end if
       ! Associate the pointers with the appropriate property routines.
       Radius_Get   => Node_Component_Disk_Standard_Radius_Solve
       Radius_Set   => Node_Component_Disk_Standard_Radius_Solve_Set
       Velocity_Get => Node_Component_Disk_Standard_Velocity
       Velocity_Set => Node_Component_Disk_Standard_Velocity_Set
    end select
    return
  end subroutine Node_Component_Disk_Standard_Radius_Solver

  subroutine mergerTreeExtraOutput(self,node,iOutput,treeIndex,nodePassesFilter,treeLock)
    !!{
    Update the star formation history after an output time is reached.
    !!}
    use            :: Galacticus_Nodes          , only : defaultDiskComponent, nodeComponentDisk, nodeComponentDiskStandard, treeNode
    use            :: Galactic_Structure_Options, only : componentTypeDisk
    use            :: Histories                 , only : history
    use, intrinsic :: ISO_C_Binding             , only : c_size_t
    use            :: Kind_Numbers              , only : kind_int8
    use            :: Locks                     , only : ompLock
    implicit none
    class  (*                ), intent(inout)          :: self
    type   (treeNode         ), intent(inout)          :: node
    integer(c_size_t         ), intent(in   )          :: iOutput
    integer(kind=kind_int8   ), intent(in   )          :: treeIndex
    logical                   , intent(in   )          :: nodePassesFilter
    type   (ompLock          ), intent(inout)          :: treeLock
    class  (nodeComponentDisk)               , pointer :: disk
    type   (history          )                         :: historyStarFormation
    !$GLC attributes unused :: self, treeIndex, nodePassesFilter, treeLock
    
    ! Check if we are the default method.
    if (.not.defaultDiskComponent%standardIsActive()) return
    ! Output the star formation history if a disk exists for this component.
    disk                 => node%disk                ()
    historyStarFormation =  disk%starFormationHistory()
    call starFormationHistory_%update(node,iOutput,historyStarFormation)
    ! Update the star formation history only if a disk exists.
    select type (disk)
    class is (nodeComponentDiskStandard)
       call disk%starFormationHistorySet(historyStarFormation)
    end select
    return
  end subroutine mergerTreeExtraOutput

  !![
  <stateStoreTask>
   <unitName>Node_Component_Disk_Standard_State_Store</unitName>
  </stateStoreTask>
  !!]
  subroutine Node_Component_Disk_Standard_State_Store(stateFile,gslStateFile,stateOperationID)
    !!{
    Write the tabulation state to file.
    !!}
    use            :: Display                          , only : displayMessage          , verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding                    , only : c_ptr                   , c_size_t
    use            :: Node_Component_Disk_Standard_Data, only : massDistributionStellar_, massDistributionGas_, kinematicDistribution_
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Storing state for: componentDisk -> standard',verbosity=verbosityLevelInfo)
    !![
    <stateStore variables="massDistributionStellar_ massDistributionGas_ kinematicDistribution_ darkMatterHaloScale_ stellarPopulationProperties_ starFormationHistory_ mergerMassMovements_"/>
    !!]
    write (stateFile) ratioAngularMomentumSolverRadius,diskRadiusSolverFlatVsSphericalFactor
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
    use            :: Display                          , only : displayMessage          , verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding                    , only : c_ptr                   , c_size_t
    use            :: Node_Component_Disk_Standard_Data, only : massDistributionStellar_, massDistributionGas_, kinematicDistribution_
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Retrieving state for: componentDisk -> standard',verbosity=verbosityLevelInfo)
    !![
    <stateRestore variables="massDistributionStellar_ massDistributionGas_ kinematicDistribution_ darkMatterHaloScale_ stellarPopulationProperties_ starFormationHistory_ mergerMassMovements_"/>
    !!]
    read (stateFile) ratioAngularMomentumSolverRadius,diskRadiusSolverFlatVsSphericalFactor
    return
  end subroutine Node_Component_Disk_Standard_State_Retrieve
  
end module Node_Component_Disk_Standard
