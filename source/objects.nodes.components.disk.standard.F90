!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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
  use Galacticus_Nodes
  use ISO_Varying_String
  implicit none
  private
  public :: Node_Component_Disk_Standard_Scale_Set                    , Node_Component_Disk_Standard_Pre_Evolve       , &
       &    Node_Component_Disk_Standard_Radius_Solver_Plausibility   , Node_Component_Disk_Standard_Radius_Solver    , &
       &    Node_Component_Disk_Standard_Star_Formation_History_Output, Node_Component_Disk_Standard_Rate_Compute     , &
       &    Node_Component_Disk_Standard_Initialize                   , Node_Component_Disk_Standard_Post_Evolve      , &
       &    Node_Component_Disk_Standard_Satellite_Merging            , Node_Component_Disk_Standard_Calculation_Reset, &
       &    Node_Component_Disk_Standard_State_Store                  , Node_Component_Disk_Standard_State_Retrieve   , &
       &    Node_Component_Disk_Standard_Thread_Initialize

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
  !#     <name>starFormationRate</name>
  !#     <attributes isSettable="false" isGettable="true" isEvolvable="false" isDeferred="get" isVirtual="true" />
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <output condition="[[diskOutputStarFormationRate]]" unitsInSI="massSolar/gigaYear" comment="Disk star formation rate."/>
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
  !#   <binding method="density"               function="Node_Component_Disk_Standard_Density"                 bindsTo="component" />
  !#   <binding method="potential"             function="Node_Component_Disk_Standard_Potential"               bindsTo="component" />
  !#   <binding method="rotationCurve"         function="Node_Component_Disk_Standard_Rotation_Curve"          bindsTo="component" />
  !#   <binding method="rotationCurveGradient" function="Node_Component_Disk_Standard_Rotation_Curve_Gradient" bindsTo="component" />
  !#   <binding method="surfaceDensity"        function="Node_Component_Disk_Standard_Surface_Density"         bindsTo="component" />
  !#  </bindings>
  !#  <functions>objects.nodes.components.disk.standard.bound_functions.inc</functions>
  !# </component>

  ! Internal count of abundances.
  integer                                     :: abundancesCount

  ! Parameters controlling the physical implementation.
  double precision                            :: diskMassToleranceAbsolute                   , diskOutflowTimescaleMinimum          , &
       &                                         diskStructureSolverRadius
  logical                                     :: diskNegativeAngularMomentumAllowed          , diskRadiusSolverCole2000Method       , &
       &                                         diskStarFormationInSatellites

  ! History of trial radii used to check for oscillations in the solution when solving for the structure of the disk.
  integer                                     :: radiusSolverIteration
  double precision                            :: radiusHistory                     (2)
  !$omp threadprivate(radiusHistory,radiusSolverIteration)
  ! The largest and smallest angular momentum, in units of that of a circular orbit at the virial radius, considered to be physically plausible for a disk.
  double precision, parameter                 :: angularMomentumMaximum               =1.0d1
  double precision, parameter                 :: angularMomentumMinimum               =1.0d-6

  ! Disk structural parameters.
  type            (varying_string)            :: diskMassDistributionName
  double precision                            :: diskStructureSolverSpecificAngularMomentum  , diskRadiusSolverFlatVsSphericalFactor
  !$omp threadprivate(diskStructureSolverSpecificAngularMomentum,diskRadiusSolverFlatVsSphericalFactor)
  
  ! Record of whether this module has been initialized.
  logical                                     :: moduleInitialized                    =.false.

contains

  !# <nodeComponentInitializationTask>
  !#  <unitName>Node_Component_Disk_Standard_Initialize</unitName>
  !# </nodeComponentInitializationTask>
  subroutine Node_Component_Disk_Standard_Initialize()
    !% Initializes the tree node standard disk methods module.
    use Input_Parameters
    use Abundances_Structure
    use Galacticus_Error
    use Node_Component_Disk_Standard_Data
    implicit none
    type(nodeComponentDiskStandard) :: diskStandardComponent

    ! Initialize the module if necessary.
    !$omp critical (Node_Component_Disk_Standard_Initialize)
    if (defaultDiskComponent%standardIsActive().and..not.moduleInitialized) then

       ! Get number of abundance properties.
       abundancesCount  =Abundances_Property_Count            ()

       ! Attach the cooling mass/angular momentum pipes from the hot halo component.
       call diskStandardComponent%attachPipes()

       ! Bind the star formation rate function.
       call diskStandardComponent%starFormationRateFunction(Node_Component_Disk_Standard_Star_Formation_Rate)

       ! Read parameters controlling the physical implementation.
       !@ <inputParameter>
       !@   <name>diskMassToleranceAbsolute</name>
       !@   <defaultValue>$10^{-6} M_\odot$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The mass tolerance used to judge whether the disk is physically plausible.
       !@   </description>
       !@   <type>double</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('diskMassToleranceAbsolute',diskMassToleranceAbsolute,defaultValue=1.0d-6)
       !@ <inputParameter>
       !@   <name>diskOutflowTimescaleMinimum</name>
       !@   <defaultValue>$10^{-3}$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The minimum timescale (in units of the disk dynamical time) on which outflows may deplete gas in the disk.
       !@   </description>
       !@   <type>double</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('diskOutflowTimescaleMinimum',diskOutflowTimescaleMinimum,defaultValue=1.0d-3)
       !@ <inputParameter>
       !@   <name>diskStructureSolverRadius</name>
       !@   <defaultValue>1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The radius (in units of the standard scale length) to use in solving for the size of the disk.
       !@   </description>
       !@   <type>double</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('diskStructureSolverRadius',diskStructureSolverRadius,defaultValue=1.0d0)
       !@ <inputParameter>
       !@   <name>diskRadiusSolverCole2000Method</name>
       !@   <defaultValue>1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('diskRadiusSolverCole2000Method',diskRadiusSolverCole2000Method,defaultValue=.false.)
       !@ <inputParameter>
       !@   <name>heightToRadialScaleDisk</name>
       !@   <defaultValue>0.137 \citep{kregel_flattening_2002}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The ratio of scale height to scale radius for standard disks.
       !@   </description>
       !@   <type>double</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('heightToRadialScaleDisk',heightToRadialScaleDisk,defaultValue=0.137d0)
       !@ <inputParameter>
       !@   <name>diskNegativeAngularMomentumAllowed</name>
       !@   <defaultValue>true</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether or not negative angular momentum is allowed for the disk.
       !@   </description>
       !@   <type>double</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('diskNegativeAngularMomentumAllowed',diskNegativeAngularMomentumAllowed,defaultValue=.true.)
       !@ <inputParameter>
       !@   <name>diskStarFormationInSatellites</name>
       !@   <defaultValue>true</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether or not star formation occurs in disks in satellites.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('diskStarFormationInSatellites',diskStarFormationInSatellites,defaultValue=.true.)
       ! Create the disk mass distribution.
       !@ <inputParameter>
       !@   <name>diskMassDistribution</name>
       !@   <defaultValue>exponentialDisk</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The type of mass distribution to use for the standard disk component.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('diskMassDistribution',diskMassDistributionName,defaultValue="exponentialDisk")
       ! Record that the module is now initialized.
       moduleInitialized=.true.
    end if
    !$omp end critical (Node_Component_Disk_Standard_Initialize)
    return
  end subroutine Node_Component_Disk_Standard_Initialize

  !# <mergerTreeEvolveThreadInitialize>
  !#  <unitName>Node_Component_Disk_Standard_Thread_Initialize</unitName>
  !# </mergerTreeEvolveThreadInitialize>
  subroutine Node_Component_Disk_Standard_Thread_Initialize
    !% Initializes the tree node hot halo methods module.
    use Galacticus_Error
    use Node_Component_Disk_Standard_Data
    implicit none
    double precision :: diskMassDistributionDensityMoment1, diskMassDistributionDensityMoment2
    logical          :: surfaceDensityMoment1IsInfinite   , surfaceDensityMoment2IsInfinite

    ! Check if this implementation is selected. If so, initialize the mass distribution.
    if (defaultDiskComponent%standardIsActive()) then
       diskMassDistribution => Mass_Distribution_Create(char(diskMassDistributionName))
       select type (diskMassDistribution)
       type is (massDistributionExponentialDisk)
          call diskMassDistribution%initialize(scaleHeight=heightToRadialScaleDisk,isDimensionless=.true.)
       type is (massDistributionMiyamotoNagai  )
          call diskMassDistribution%initialize(b          =heightToRadialScaleDisk,isDimensionless=.true.)
       class default
          call Galacticus_Error_Report('Node_Component_Disk_Standard_Thread_Initialize','unsupported mass distribution')
       end select
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
          call Galacticus_Error_Report('Node_Component_Disk_Standard_Thread_Initialize','only cylcindrically symmetric mass distributions are allowed')
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

  !# <calculationResetTask>
  !#   <unitName>Node_Component_Disk_Standard_Calculation_Reset</unitName>
  !# </calculationResetTask>
  subroutine Node_Component_Disk_Standard_Calculation_Reset(node)
    !% Reset standard disk structure calculations.
    use Node_Component_Disk_Standard_Data
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
    implicit none
    type (treeNode         ), intent(inout), pointer :: node
    class(nodeComponentDisk)               , pointer :: disk

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

  !# <postEvolveTask>
  !# <unitName>Node_Component_Disk_Standard_Post_Evolve</unitName>
  !# </postEvolveTask>
  subroutine Node_Component_Disk_Standard_Post_Evolve(node)
    !% Trim histories attached to the disk.
    use Galacticus_Display
    use String_Handling
    use ISO_Varying_String
    use Abundances_Structure
    use Histories
    use Galacticus_Error
    use Dark_Matter_Halo_Scales
    use Stellar_Luminosities_Structure
    implicit none
    type            (treeNode                ), intent(inout), pointer :: node
    class           (nodeComponentDisk       )               , pointer :: disk
    class           (nodeComponentBasic      )               , pointer :: basic
    class           (nodeComponentSpin       )               , pointer :: spin
    class           (darkMatterHaloScaleClass)               , pointer :: darkMatterHaloScale_
    double precision                          , parameter              :: angularMomentumTolerance=1.0d-2
    double precision                          , save                   :: fractionalErrorMaximum  =0.0d0
    double precision                                                   :: diskMass                       , fractionalError, &
         &                                                                specificAngularMomentum
    character       (len=20                  )                         :: valueString
    type            (varying_string          )                         :: message
    type            (history                 )                         :: stellarPropertiesHistory

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
             call disk%luminositiesStellarSet(zeroStellarLuminosities)
          else
             specificAngularMomentum=disk%angularMomentum()/diskMass
             if (specificAngularMomentum < 0.0d0) specificAngularMomentum=disk%radius()*disk%velocity()
          end if

          ! Reset the gas, abundances and angular momentum of the disk.
          call disk%        massGasSet(                                     0.0d0)
          call disk%  abundancesGasSet(                            zeroAbundances)
          call disk%angularMomentumSet(specificAngularMomentum*disk%massStellar())

       end if
       ! Trap negative angular momentum.
       darkMatterHaloScale_ => darkMatterHaloScale()
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
                write (valueString,'(e12.6)')  darkMatterHaloScale_%virialRadius  (node) &
                     &                        *darkMatterHaloScale_%virialVelocity(node) &
                     &                        *spin                %spin          (    )
                message=message//char(10)//' -> angular momentum scale = '//trim(valueString)
                call Galacticus_Error_Report(                                                &
                     &                       'Node_Component_Disk_Standard_Post_Evolve',     &
                     &                       message                                         &
                     &                      )
             end if
          end if
       end if
    end select
    return
  end subroutine Node_Component_Disk_Standard_Post_Evolve

  subroutine Node_Component_Disk_Standard_Create(node)
    !% Create properties in an standard disk component.
    use Histories
    use Stellar_Population_Properties
    use Galacticus_Output_Star_Formation_Histories
    implicit none
    type   (treeNode             ), intent(inout), pointer :: node
    class  (nodeComponentDisk    )               , pointer :: disk
    class  (nodeComponentSpheroid)               , pointer :: spheroid
    class  (nodeComponentBasic   )               , pointer :: basic
    type   (history              )                         :: starFormationHistory        , stellarPropertiesHistory      , &
         &                                                    spheroidStarFormationHistory
    logical                                                :: createStarFormationHistory  , createStellarPropertiesHistory
    double precision                                       :: timeBegin
    
    ! Get the disk component.
    disk => node%disk()
    ! Exit if already initialized.
    if (disk%isInitialized()) return
    ! Determine which histories must be created.
    starFormationHistory          =disk%starFormationHistory            ()
    createStarFormationHistory    =.not.             starFormationHistory    %exists ()
    call                                             starformationhistory    %destroy()
    stellarPropertiesHistory      =disk%stellarPropertiesHistory        ()
    createStellarPropertiesHistory=.not.             stellarPropertiesHistory%exists ()
    call                                             stellarPropertiesHistory%destroy()
    ! Create the stellar properties history.
    if (createStellarPropertiesHistory) then
       ! Create the stellar properties history.
       call Stellar_Population_Properties_History_Create (node,stellarPropertiesHistory)
       call disk%stellarPropertiesHistorySet(         stellarPropertiesHistory)
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
       call Star_Formation_History_Create                (node,    starFormationHistory,timeBegin)
       call disk%    starFormationHistorySet(             starFormationHistory          )
    end if
    ! Record that the disk has been initialized.
    call disk%isInitializedSet(.true.)
    return
  end subroutine Node_Component_Disk_Standard_Create

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Disk_Standard_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Disk_Standard_Rate_Compute(node,odeConverged,interrupt,interruptProcedureReturn)
    !% Compute the standard disk node mass rate of change.
    use Abundances_Structure
    use Histories
    use Star_Formation_Feedback_Disks
    use Star_Formation_Feedback_Expulsion_Disks
    use Galactic_Structure_Options
    use Galactic_Dynamics_Bar_Instabilities
    use Galacticus_Output_Star_Formation_Histories
    use Stellar_Population_Properties
    use Numerical_Constants_Astronomical
    use Ram_Pressure_Stripping_Mass_Loss_Rate_Disks
    use Tidal_Stripping_Mass_Loss_Rate_Disks
    use Dark_Matter_Halo_Scales
    use Stellar_Luminosities_Structure
    implicit none
    type            (treeNode                    ), intent(inout), pointer :: node
    logical                                       , intent(in   )          :: odeConverged
    class           (nodeComponentDisk           )               , pointer :: disk
    class           (nodeComponentSpheroid       )               , pointer :: spheroid
    class           (nodeComponentHotHalo        )               , pointer :: hotHalo
    class           (darkMatterHaloScaleClass    )               , pointer :: darkMatterHaloScale_
    logical                                       , intent(inout)          :: interrupt
    procedure       (interruptTask               ), intent(inout), pointer :: interruptProcedureReturn
    procedure       (interruptTask               )               , pointer :: interruptProcedure
    type            (abundances                  ), save                   :: fuelAbundances              , fuelAbundancesRates       , &
         &                                                                    stellarAbundancesRates
    !$omp threadprivate(fuelAbundances,stellarAbundancesRates,fuelAbundancesRates)
    double precision                                                        :: angularMomentum             , angularMomentumOutflowRate, &
         &                                                                     barInstabilitySpecificTorque, barInstabilityTimescale   , &
         &                                                                     diskDynamicalTime           , diskMass                  , &
         &                                                                     energyInputRate             , fractionGas               , &
         &                                                                     fractionStellar             , fuelMass                  , &
         &                                                                     fuelMassRate                , gasMass                   , &
         &                                                                     massLossRate                , massOutflowRate           , &
         &                                                                     massOutflowRateFromHalo     , massOutflowRateToHotHalo  , &
         &                                                                     outflowToHotHaloFraction    , starFormationRate         , &
         &                                                                     stellarMassRate             , transferRate
    type            (history                      )                         :: historyTransferRate         , stellarHistoryRate
    type            (stellarLuminosities          ), save                   :: luminositiesStellarRates    , luminositiesTransferRate
    !$omp threadprivate(luminositiesStellarRates,luminositiesTransferRate)
    !GCC$ attributes unused :: odeConverged
    
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

       ! Compute the star formation rate.
       starFormationRate=disk%starFormationRate()

       ! Get the available fuel mass.
       fuelMass         =disk%massGas          ()

       ! Find the metallicity of the fuel supply.
       fuelAbundances   =disk%abundancesGas    ()
       call fuelAbundances%massToMassFraction(fuelMass)

       ! Find rates of change of stellar mass, gas mass, abundances and luminosities.
       stellarHistoryRate=disk%stellarPropertiesHistory()
       call Stellar_Population_Properties_Rates(starFormationRate,fuelAbundances,componentTypeDisk,node,stellarHistoryRate&
            &,stellarMassRate,stellarAbundancesRates,luminositiesStellarRates,fuelMassRate,fuelAbundancesRates,energyInputRate)
       if (stellarHistoryRate%exists()) call disk%stellarPropertiesHistoryRate(stellarHistoryRate)

       ! Adjust rates.
       call disk%        massStellarRate(         stellarMassRate)
       call disk%            massGasRate(            fuelMassRate)
       call disk%  abundancesStellarRate(  stellarAbundancesRates)
       call disk%      abundancesGasRate(     fuelAbundancesRates)
       call disk%luminositiesStellarRate(luminositiesStellarRates)

       ! Record the star formation history.
       stellarHistoryRate=disk%starFormationHistory()
       call stellarHistoryRate%reset()
       call Star_Formation_History_Record(node,stellarHistoryRate,fuelAbundances,starFormationRate)
       if (stellarHistoryRate%exists()) call disk%starFormationHistoryRate(stellarHistoryRate)

       ! Find rate of outflow of material from the disk and pipe it to the outflowed reservoir.
       massOutflowRateToHotHalo=Star_Formation_Feedback_Disk_Outflow_Rate          (node,starFormationRate,energyInputRate)
       massOutflowRateFromHalo =Star_Formation_Expulsive_Feedback_Disk_Outflow_Rate(node,starFormationRate,energyInputRate)
       massOutflowRate         =massOutflowRateToHotHalo+massOutflowRateFromHalo
       if (massOutflowRate > 0.0d0) then
          ! Find the fraction of material which outflows to the hot halo.
          outflowToHotHaloFraction=massOutflowRateToHotHalo/massOutflowRate

          ! Get the masses of the disk.
          gasMass =        disk%massGas    ()
          diskMass=gasMass+disk%massStellar()

          ! Limit the outflow rate timescale to a multiple of the dynamical time.
          diskDynamicalTime=Mpc_per_km_per_s_To_Gyr*disk%radius()/disk%velocity()
          massOutflowRate=min(massOutflowRate,gasMass/diskOutflowTimescaleMinimum/diskDynamicalTime)

          ! Compute the angular momentum outflow rate.
          if (diskMass > 0.0d0) then
             angularMomentum           =disk%angularMomentum()
             angularMomentumOutflowRate=angularMomentum*(massOutflowRate/diskMass)
             angularMomentumOutflowRate=sign(min(abs(angularMomentumOutflowRate),abs(angularMomentum/diskOutflowTimescaleMinimum &
                  &/diskDynamicalTime)),angularMomentumOutflowRate)
          else
             angularMomentumOutflowRate=0.0d0
          end if
          if (gasMass > 0.0d0) then
             fuelAbundancesRates=disk%abundancesGas()
             call fuelAbundancesRates%massToMassFraction(gasMass)
             fuelAbundancesRates=fuelAbundancesRates*massOutflowRate
          else
             fuelAbundancesRates=zeroAbundances
          end if
          hotHalo => node%hotHalo()
          call hotHalo%           outflowingMassRate( massOutflowRate           *outflowToHotHaloFraction)
          call disk   %                  massGasRate(-massOutflowRate                                    )
          call hotHalo%outflowingAngularMomentumRate( angularMomentumOutflowRate*outflowToHotHaloFraction)
          call disk   %          angularMomentumRate(-angularMomentumOutflowRate                         )
          call hotHalo%     outflowingAbundancesRate( fuelAbundancesRates       *outflowToHotHaloFraction)
          call disk   %            abundancesGasRate(-fuelAbundancesRates                                )
       end if

       ! Determine if the disk is bar unstable and, if so, the rate at which material is moved to the pseudo-bulge.
       if (node%isPhysicallyPlausible) then
          ! Disk has positive angular momentum, so compute an instability timescale.
          call Bar_Instability_Timescale(node,barInstabilityTimescale,barInstabilitySpecificTorque)
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
          transferRate            =max(         0.0d0         ,disk    %massGas                (                         ))/barInstabilityTimescale
          call                                      disk    %massGasRate            (-           transferRate                              )
          call                                      spheroid%massGasRate            (+           transferRate ,interrupt,interruptProcedure)
          ! Stellar mass.
          transferRate            =max(         0.0d0         ,disk    %massStellar            (                         ))/barInstabilityTimescale
          call                                      disk    %massStellarRate        (-           transferRate                              )
          call                                      spheroid%massStellarRate        (+           transferRate ,interrupt,interruptProcedure)
          ! Angular momentum.
          transferRate            =max(         0.0d0         ,disk    %angularMomentum        (                         ))/barInstabilityTimescale
          call                                      disk    %angularMomentumRate    (-           transferRate                              )
          call                                      spheroid%angularMomentumRate    (+           transferRate ,interrupt,interruptProcedure)
          ! Gas abundances.
          fuelAbundancesRates     =max(zeroAbundances         ,disk    %abundancesGas          (                         ))/barInstabilityTimescale
          call                                      disk    %abundancesGasRate      (-     fuelAbundancesRates                             )
          call                                      spheroid%abundancesGasRate      (+     fuelAbundancesRates,interrupt,interruptProcedure)
          ! Stellar abundances.
          stellarAbundancesRates  =max(zeroAbundances         ,disk    %abundancesStellar      (                         ))/barInstabilityTimescale
          call                                      disk    %abundancesStellarRate  (-  stellarAbundancesRates                             )
          call                                      spheroid%abundancesStellarRate  (+  stellarAbundancesRates,interrupt,interruptProcedure)
          ! Stellar luminosities.
          luminositiesTransferRate=max(zeroStellarLuminosities,disk    %luminositiesStellar    (                         ))/barInstabilityTimescale
          call                                      disk    %luminositiesStellarRate(-luminositiesTransferRate                             )
          call                                      spheroid%luminositiesStellarRate(+luminositiesTransferRate,interrupt,interruptProcedure)
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
          darkMatterHaloScale_ => darkMatterHaloScale()
          if (spheroid%angularMomentum() < (spheroid%massGas()+spheroid%massStellar())*darkMatterHaloScale_%virialRadius(node)*darkMatterHaloScale_%virialVelocity(node) .and. spheroid%radius() < darkMatterHaloScale_%virialRadius(node)) then
             call spheroid%angularMomentumRate(+barInstabilitySpecificTorque*(spheroid%massGas()+spheroid%massStellar()),interrupt,interruptProcedure)
          end if
       end if

       ! Apply mass loss rate due to ram pressure stripping.
       if (disk%massGas() > 0.0d0) then
          massLossRate=Ram_Pressure_Stripping_Mass_Loss_Rate_Disk(node)
          if (massLossRate > 0.0d0) then
             hotHalo => node%hotHalo()
             call    disk%                  massGasRate(-massLossRate                                                                       )
             call    disk%          angularMomentumRate(-massLossRate*disk%angularMomentum()/(disk%massGas()+disk%massStellar()))
             call    disk%            abundancesGasRate(-massLossRate*disk%abundancesGas  ()/ disk%massGas()                        )
             call hotHalo%           outflowingMassRate(+massLossRate                                                                       )
             call hotHalo%outflowingAngularMomentumRate(+massLossRate*disk%angularMomentum()/(disk%massGas()+disk%massStellar()))
             call hotHalo%outflowingAbundancesRate     (+massLossRate*disk%abundancesGas  ()/ disk%massGas()                        )
          end if
       end if

       ! Apply mass loss rate due to tidal stripping.
       if (disk%massGas()+disk%massStellar() > 0.0d0) then
          massLossRate=Tidal_Stripping_Mass_Loss_Rate_Disk(node)
          if (massLossRate > 0.0d0) then
             hotHalo    => node%hotHalo()
             fractionGas    =  min(1.0d0,max(0.0d0,disk%massGas()/(disk%massGas()+disk%massStellar())))
             fractionStellar=  1.0d0-fractionGas
             if (fractionGas     > 0.0d0 .and. disk%massGas    () > 0.0d0) then
                call    disk%                  massGasRate(-fractionGas    *massLossRate                                                                         )
                call    disk%            abundancesGasRate(-fractionGas    *massLossRate*disk%abundancesGas    ()/ disk%massGas()                        )
                call hotHalo%           outflowingMassRate(+fractionGas    *massLossRate                                                                         )
                call hotHalo%outflowingAbundancesRate     (+fractionGas    *massLossRate*disk%abundancesGas    ()/ disk%massGas()                        )
                call hotHalo%outflowingAngularMomentumRate(+fractionGas    *massLossRate*disk%angularMomentum  ()/(disk%massGas()+disk%massStellar()))
             end if
             if (fractionStellar > 0.0d0 .and. disk%massStellar() > 0.0d0) then
                call    disk%              massStellarRate(-fractionStellar*massLossRate                                                                         )
                call    disk%        abundancesStellarRate(-fractionStellar*massLossRate*disk%abundancesStellar()/                    disk%massStellar() )
                ! Stellar properties history.
                historyTransferRate=disk%stellarPropertiesHistory()         
                if (historyTransferRate%exists()) then
                   call disk%stellarPropertiesHistoryRate (-fractionStellar*massLossRate*historyTransferRate         /                    disk%massStellar() )
                end if
                call historyTransferRate%destroy()
                ! Star formation history.
                historyTransferRate=disk%starFormationHistory()
                if (historyTransferRate%exists()) then
                   call disk    %starFormationHistoryRate (-fractionStellar*massLossRate*historyTransferRate         /                    disk%massStellar() )
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
    use Histories
    use Stellar_Population_Properties
    use Galacticus_Output_Star_Formation_Histories
    use Abundances_Structure
    use Stellar_Luminosities_Structure
    implicit none
    type            (treeNode             ), intent(inout), pointer :: node
    class           (nodeComponentDisk    )               , pointer :: disk
    class           (nodeComponentSpheroid)               , pointer :: spheroid
    double precision                       , parameter              :: massMinimum                   =1.0d0
    double precision                       , parameter              :: angularMomentumMinimum        =1.0d-2
    double precision                       , parameter              :: luminosityMinimum             =1.0d0
    double precision                                                :: angularMomentum                      , mass
    type            (history              )                         :: stellarPopulationHistoryScales
    type            (stellarLuminosities  )                         :: stellarLuminositiesScale
    type            (abundances           )                         :: abundancesScale
    
    ! Get the disk component.
    disk => node%disk()
    ! Check if an standard disk component exists.
    select type (disk)
    class is (nodeComponentDiskStandard)
       ! Get spheroid component.
       spheroid => node%spheroid()
       ! Set scale for angular momentum.
       angularMomentum=disk%angularMomentum()+spheroid%angularMomentum()
       call disk%angularMomentumScale(max(angularMomentum,angularMomentumMinimum))
       ! Set scale for masses.
       mass           =abs(                                           &
            &              +disk%massGas    ()+spheroid%massGas    () &
            &              +disk%massStellar()+spheroid%massStellar() &
            &             )
       call disk%massGasScale        (max(           mass,           massMinimum))
       call disk%massStellarScale    (max(           mass,           massMinimum))
       ! Set scales for abundances if necessary.
       if (abundancesCount > 0) then
          ! Set scale for abundances.
          abundancesScale=+max(                                    &
               &               +abs(                               &
               &                    +disk    %abundancesGas    ()  &
               &                    +spheroid%abundancesGas    ()  &
               &                   )                               &
               &               +abs(                               &
               &                    +disk    %abundancesStellar()  &
               &                    +spheroid%abundancesStellar()  &
               &                   )                             , &
               &                    +massMinimum                   &
               &                    *unitAbundances                &
               &              )
          ! Set scale for gas abundances.
          call disk%abundancesGasScale    (abundancesScale)

          ! Set scale for stellar abundances.
          call disk%abundancesStellarScale(abundancesScale)
       end if
       ! Set scale for stellar luminosities.
       stellarLuminositiesScale=max(                                       &
            &                       abs(                                   &
            &                           +disk      %luminositiesStellar()  &
            &                           +spheroid  %luminositiesStellar()  &
            &                          )                                 , &
            &                           +unitStellarLuminosities           &
            &                           *luminosityMinimum                 &
            &                      )
       call stellarLuminositiesScale%truncate                (disk       %luminositiesStellar())
       call disk       %luminositiesStellarScale(stellarLuminositiesScale                      )

       ! Set scales for stellar population properties and star formation histories.
       stellarPopulationHistoryScales=disk%stellarPropertiesHistory()
       call Stellar_Population_Properties_Scales(stellarPopulationHistoryScales,disk%massStellar(),disk%abundancesStellar())
       call disk%stellarPropertiesHistoryScale  (stellarPopulationHistoryScales                                            )
       call stellarPopulationHistoryScales%destroy()
       stellarPopulationHistoryScales=disk%starFormationHistory()
       call Star_Formation_History_Scales       (stellarPopulationHistoryScales,disk%massStellar(),disk%abundancesStellar())
       call disk%starFormationHistoryScale      (stellarPopulationHistoryScales                                            )
       call stellarPopulationHistoryScales%destroy()

    end select
    return
  end subroutine Node_Component_Disk_Standard_Scale_Set

  !# <satelliteMergerTask>
  !#  <unitName>Node_Component_Disk_Standard_Satellite_Merging</unitName>
  !#  <after>Satellite_Merging_Mass_Movement_Store</after>
  !#  <after>Satellite_Merging_Remnant_Size</after>
  !# </satelliteMergerTask>
  subroutine Node_Component_Disk_Standard_Satellite_Merging(node)
    !% Transfer any standard disk associated with {\normalfont \ttfamily node} to its host halo.
    use Histories
    use Abundances_Structure
    use Satellite_Merging_Mass_Movements_Descriptors
    use Galacticus_Error
    use Stellar_Luminosities_Structure
    implicit none
    type            (treeNode             ), intent(inout), pointer :: node
    class           (nodeComponentDisk    )               , pointer :: diskHost               , disk
    class           (nodeComponentSpheroid)               , pointer :: spheroidHost           , spheroid
    type            (treeNode             )               , pointer :: nodeHost
    type            (history              )                         :: historyHost            , historyNode
    double precision                                                :: specificAngularMomentum

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

       ! Move the gas component of the standard disk to the host.
       select case (thisMergerGasMovesTo)
       case (movesToDisk)
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
       case (movesToSpheroid)
          call spheroidHost%massGasSet            (                                                                     &
               &                                             spheroidHost%massGas            ()                         &
               &                                            +disk        %massGas            ()                         &
               &                                           )
          call spheroidHost%abundancesGasSet      (                                                                     &
               &                                             spheroidHost%abundancesGas      ()                         &
               &                                            +disk        %abundancesGas      ()                         &
               &                                           )
       case default
          call Galacticus_Error_Report('Node_Component_Disk_Standard_Satellite_Merging','unrecognized movesTo descriptor')
       end select
       call disk%      massGasSet(         0.0d0)
       call disk%abundancesGasSet(zeroAbundances)

       ! Move the stellar component of the standard disk to the host.
       select case (thisMergerStarsMoveTo)
       case (movesToDisk)
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
          historyHost=diskHost    %stellarPropertiesHistory()
          call historyHost%increment(historyNode)
          call historyNode%reset    (           )
          call diskHost%stellarPropertiesHistorySet(historyHost)
          call disk%stellarPropertiesHistorySet(historyNode)
          ! Also add star formation histories.
          historyNode=disk    %starFormationHistory    ()
          historyHost=diskHost    %starFormationHistory    ()
          call historyHost%combine (historyNode)
          call historyNode%reset   (           )
          call diskHost    %starFormationHistorySet    (historyHost)
          call disk    %starFormationHistorySet    (historyNode)
          call historyNode%destroy(recordMemory=.false.)
          call historyHost%destroy(recordMemory=.false.)
       case (movesToSpheroid)
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
          call historyHost%increment(historyNode)
          call historyNode%reset    (           )
          call spheroidHost%stellarPropertiesHistorySet(historyHost)
          call disk    %stellarPropertiesHistorySet(historyNode)
          ! Also add star formation histories.
          historyNode=disk    %starFormationHistory    ()
          historyHost=spheroidHost%starFormationHistory    ()
          call historyHost%combine(historyNode)
          call historyNode%reset  (           )
          call spheroidHost%starFormationHistorySet(historyHost)
          call disk    %starFormationHistorySet(historyNode)
          call historyNode%destroy(recordMemory=.false.)
          call historyHost%destroy(recordMemory=.false.)
       case default
          call Galacticus_Error_Report('Node_Component_Disk_Standard_Satellite_Merging','unrecognized movesTo descriptor')
       end select
       call disk%        massStellarSet(                  0.0d0)
       call disk%  abundancesStellarSet(         zeroAbundances)
       call disk%luminositiesStellarSet(zeroStellarLuminosities)
       call disk%    angularMomentumSet(                  0.0d0)
    end select
    return
  end subroutine Node_Component_Disk_Standard_Satellite_Merging

  !# <radiusSolverPlausibility>
  !#  <unitName>Node_Component_Disk_Standard_Radius_Solver_Plausibility</unitName>
  !# </radiusSolverPlausibility>
  subroutine Node_Component_Disk_Standard_Radius_Solver_Plausibility(node,galaxyIsPhysicallyPlausible)
    !% Determines whether the disk is physically plausible for radius solving tasks. Require that it have non-zero mass and angular momentum.
    use Dark_Matter_Halo_Scales
    implicit none
    type            (treeNode                ), intent(inout) :: node
    logical                                   , intent(inout) :: galaxyIsPhysicallyPlausible
    class           (nodeComponentDisk       ), pointer       :: disk
    class           (darkMatterHaloScaleClass), pointer       :: darkMatterHaloScale_
    double precision                                          :: angularMomentumScale

    ! Return immediately if our method is not selected.
    if (.not.defaultDiskComponent%standardIsActive()) return

     ! Determine the plausibility of the current disk.
     disk => node%disk()
     select type (disk)
     class is (nodeComponentDiskStandard)
        if      (disk%angularMomentum()                             <                       0.0d0) &
             & galaxyIsPhysicallyPlausible=.false.
        if      (disk%massStellar    ()+disk%massGas() <  -diskMassToleranceAbsolute) then
           galaxyIsPhysicallyPlausible=.false.
        else if (disk%massStellar    ()+disk%massGas() >=                      0.0d0) then
           if      (                                                                               &
                &   disk%angularMomentum() < 0.0d0                                                 &
                &  ) then
              galaxyIsPhysicallyPlausible=.false.
           else
              darkMatterHaloScale_ => darkMatterHaloScale()
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
                 galaxyIsPhysicallyPlausible=.false.
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
    implicit none
    type (treeNode         ), intent(inout) :: node
    class(nodeComponentDisk), pointer       :: disk

    disk => node%disk()
    Node_Component_Disk_Standard_Radius_Solve=disk%radius()*diskStructureSolverRadius
    return
  end function Node_Component_Disk_Standard_Radius_Solve

  subroutine Node_Component_Disk_Standard_Radius_Solve_Set(node,radius)
    !% Set the radius of the standard disk used in structure solvers.
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
    implicit none
    type (treeNode         ), intent(inout) :: node
    class(nodeComponentDisk), pointer       :: disk

    disk => node%disk()
    Node_Component_Disk_Standard_Velocity=disk%velocity()
    return
  end function Node_Component_Disk_Standard_Velocity

  subroutine Node_Component_Disk_Standard_Velocity_Set(node,velocity)
    !% Set the circular velocity of the standard disk.
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
    use Tables
    use Node_Component_Disk_Standard_Data
    use Numerical_Constants_Physical
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

  double precision function Node_Component_Disk_Standard_Star_Formation_Rate(self)
    !% Return the star formation rate of the standard disk.
    use Star_Formation_Timescales_Disks
    implicit none
    class           (nodeComponentDiskStandard), intent(inout) :: self
    type            (treeNode                 ), pointer       :: node
    double precision                                           :: gasMass , starFormationTimescale

    ! Get the associated node.
    node => self%host()

    ! Get the star formation timescale.
    starFormationTimescale=Star_Formation_Timescale_Disk(node)

    ! Get the gas mass.
    gasMass=self%massGas()

    ! If timescale is finite and gas mass is positive, then compute star formation rate.
    if (starFormationTimescale > 0.0d0 .and. gasMass > 0.0d0 .and. (diskStarFormationInSatellites .or. .not.node%isSatellite())) then
       Node_Component_Disk_Standard_Star_Formation_Rate=gasMass/starFormationTimescale
    else
       Node_Component_Disk_Standard_Star_Formation_Rate=0.0d0
    end if
    return
  end function Node_Component_Disk_Standard_Star_Formation_Rate

  !# <mergerTreeExtraOutputTask>
  !#  <unitName>Node_Component_Disk_Standard_Star_Formation_History_Output</unitName>
  !# </mergerTreeExtraOutputTask>
  subroutine Node_Component_Disk_Standard_Star_Formation_History_Output(node,iOutput,treeIndex,nodePassesFilter)
    !% Store the star formation history in the output file.
    use, intrinsic :: ISO_C_Binding
    use Kind_Numbers
    use Histories
    use Galacticus_Output_Star_Formation_Histories
    implicit none
    type   (treeNode         ), intent(inout), pointer :: node
    integer(c_size_t         ), intent(in   )          :: iOutput
    integer(kind=kind_int8   ), intent(in   )          :: treeIndex
    logical                   , intent(in   )          :: nodePassesFilter
    class  (nodeComponentDisk)               , pointer :: disk
    type   (history          )                         :: starFormationHistory

    ! Output the star formation history if a disk exists for this component.
    disk => node%disk()
    select type (disk)
       class is (nodeComponentDiskStandard)
       starFormationHistory=disk%starFormationHistory()
       call Star_Formation_History_Output(node,nodePassesFilter,starFormationHistory,iOutput,treeIndex,'disk')
       call disk%starFormationHistorySet(starFormationHistory)
    end select
    return
  end subroutine Node_Component_Disk_Standard_Star_Formation_History_Output

  !# <galacticusStateStoreTask>
  !#  <unitName>Node_Component_Disk_Standard_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Node_Component_Disk_Standard_State_Store(stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    use Galacticus_Display
    use Node_Component_Disk_Standard_Data
    use FGSL
    implicit none
    integer           , intent(in   ) :: stateFile
    type   (fgsl_file), intent(in   ) :: fgslStateFile

    call Galacticus_Display_Message('Storing state for: treeNodeMethodDisk -> standard',verbosity=verbosityInfo)
    write (stateFile) diskStructureSolverSpecificAngularMomentum,diskRadiusSolverFlatVsSphericalFactor
    call diskMassDistribution%stateStore(stateFile,fgslStateFile)
    return
  end subroutine Node_Component_Disk_Standard_State_Store

  !# <galacticusStateRetrieveTask>
  !#  <unitName>Node_Component_Disk_Standard_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Node_Component_Disk_Standard_State_Retrieve(stateFile,fgslStateFile)
    !% Retrieve the tabulation state from the file.
    use Galacticus_Display
    use Node_Component_Disk_Standard_Data
    use FGSL
    implicit none
    integer           , intent(in   ) :: stateFile
    type   (fgsl_file), intent(in   ) :: fgslStateFile

    call Galacticus_Display_Message('Retrieving state for: treeNodeMethodDisk -> standard',verbosity=verbosityInfo)
    read (stateFile) diskStructureSolverSpecificAngularMomentum,diskRadiusSolverFlatVsSphericalFactor
    call diskMassDistribution%stateRestore(stateFile,fgslStateFile)
    return
  end subroutine Node_Component_Disk_Standard_State_Retrieve

end module Node_Component_Disk_Standard
