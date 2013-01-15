!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements the exponential disk node component.

module Node_Component_Disk_Exponential
  !% Implements the exponential disk node component.
  use Galacticus_Nodes
  use Tables
  implicit none
  private
  public :: Node_Component_Disk_Exponential_Scale_Set                 , Node_Component_Disk_Exponential_Density                      , &
       &    Node_Component_Disk_Exponential_Enclosed_Mass             , Node_Component_Disk_Exponential_Rotation_Curve               , &
       &    Node_Component_Disk_Exponential_Radius_Solver_Plausibility, Node_Component_Disk_Exponential_Radius_Solver                , &
       &    Node_Component_Disk_Exponential_State_Retrieve            , Node_Component_Disk_Exponential_Rate_Compute                 , &
       &    Node_Component_Disk_Exponential_Rotation_Curve_Gradient   , Node_Component_Disk_Exponential_Pre_Evolve                   , &
       &    Node_Component_Disk_Exponential_Satellite_Merging         , Node_Component_Disk_Exponential_Potential                    , &
       &    Node_Component_Disk_Exponential_State_Store               , Node_Component_Disk_Exponential_Star_Formation_History_Output, &
       &    Node_Component_Disk_Exponential_Initialize                , Node_Component_Disk_Exponential_Post_Evolve
  
  !# <component>
  !#  <class>disk</class>
  !#  <name>exponential</name>
  !#  <isDefault>yes</isDefault>
  !#  <methods>
  !#   <method>
  !#     <name>isInitialized</name>
  !#     <type>logical</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </method>
  !#   <method>
  !#     <name>massStellar</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of stars in the exponential disk."/>
  !#   </method>
  !#   <method>
  !#     <name>abundancesStellar</name>
  !#     <type>abundances</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of metals in the stellar phase of the exponential disk."/>
  !#   </method>
  !#   <method>
  !#     <name>massGas</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" makeGeneric="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of gas in the exponential disk."/>
  !#   </method>
  !#   <method>
  !#     <name>abundancesGas</name>
  !#     <type>abundances</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" makeGeneric="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of metals in the gas phase of the exponential disk."/>
  !#   </method>
  !#   <method>
  !#     <name>angularMomentum</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" makeGeneric="true" />
  !#     <output unitsInSI="massSolar*megaParsec*kilo" comment="Angular momentum of the exponential disk."/>
  !#   </method>
  !#   <method>
  !#     <name>radius</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <output unitsInSI="megaparsec" comment="Radial scale length in the exponential disk."/>
  !#   </method>
  !#   <method>
  !#     <name>halfMassRadius</name>
  !#     <attributes isSettable="false" isGettable="true" isEvolvable="false" />
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <isVirtual>yes</isVirtual>
  !#     <getFunction>Node_Component_Disk_Exponential_Half_Mass_Radius</getFunction>
  !#   </method>
  !#   <method>
  !#     <name>velocity</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <output unitsInSI="kilo" comment="Circular velocity of the exponential disk at scale length."/>
  !#   </method>
  !#   <method>
  !#     <name>starFormationRate</name>
  !#     <attributes isSettable="false" isGettable="true" isEvolvable="false" isDeferred="get" />
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <isVirtual>yes</isVirtual>
  !#     <output condition="[[diskOutputStarFormationRate]]" unitsInSI="massSolar/gigaYear" comment="Disk star formation rate."/>
  !#   </method>
  !#   <method>
  !#     <name>luminositiesStellar</name>
  !#     <type>real</type>
  !#     <rank>1</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <classDefault modules="Stellar_Population_Properties_Luminosities" count="Stellar_Population_Luminosities_Count()">0.0d0</classDefault>
  !#     <output labels="':'//Stellar_Population_Luminosities_Name({i})" count="Stellar_Population_Luminosities_Count()" condition="Stellar_Population_Luminosities_Output({i},time)" modules="Stellar_Population_Properties_Luminosities" unitsInSI="luminosityZeroPointAB" comment="Luminosity of disk stars."/>
  !#   </method>
  !#   <method>
  !#     <name>stellarPropertiesHistory</name>
  !#     <type>history</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#   </method>
  !#   <method>
  !#     <name>starFormationHistory</name>
  !#     <type>history</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#   </method>
  !#  </methods>
  !#  <bindings>
  !#   <binding method="attachPipes" function="Node_Component_Disk_Exponential_Attach_Pipes" type="void" bindsTo="component" />
  !#  </bindings>
  !#  <functions>objects.nodes.components.disk.exponential.custom_methods.inc</functions>
  !# </component>

  ! Internal count of abundances.
  integer                                     :: abundancesCount

  ! Internal count of luminosities and work arrays.
  integer                                     :: luminositiesCount
  double precision, allocatable, dimension(:) :: zeroLuminosities,luminositiesMinimum,luminositiesStellarRates,luminositiesTransferRate,luminositiesDisk
  !$omp threadprivate(zeroLuminosities,luminositiesMinimum,luminositiesStellarRates,luminositiesTransferRate,luminositiesDisk)

  ! Parameters controlling the physical implementation.
  double precision                            :: diskMassToleranceAbsolute,diskOutflowTimescaleMinimum,diskStructureSolverRadius&
       &,heightToRadialScaleDisk
  logical                                     :: diskRadiusSolverCole2000Method
  
  ! Tabulation of the exponential disk rotation curve.
  integer,          parameter                 :: rotationCurvePointsPerDecade=10
  integer                                     :: rotationCurvePointsCount
  logical                                     :: rotationCurveInitialized=.false.
  double precision, parameter                 :: rotationCurveHalfRadiusMinimumDefault=1.0d-6,rotationCurveHalfRadiusMaximumDefault=10.0d0
  double precision                            :: rotationCurveHalfRadiusMinimum=rotationCurveHalfRadiusMinimumDefault&
       &,rotationCurveHalfRadiusMaximum=rotationCurveHalfRadiusMaximumDefault
  double precision                            :: scaleLengthFactor,diskStructureSolverSpecificAngularMomentum,diskRadiusSolverFlatVsSphericalFactor
  logical                                     :: scaleLengthFactorSet=.false.
  type(table1DLogarithmicLinear)                    :: rotationCurveTable


  ! Tabulation of the exponential disk rotation curve gradient.
  integer,          parameter                 :: rotationCurveGradientPointsPerDecade=10
  integer                                     :: rotationCurveGradientPointsCount
  logical                                     :: rotationCurveGradientInitialized=.false.
  double precision, parameter                 :: rotationCurveGradientHalfRadiusMinimumDefault=1.0d-6,rotationCurveGradientHalfRadiusMaximumDefault=10.0d0
  double precision                            :: rotationCurveGradientHalfRadiusMinimum=rotationCurveGradientHalfRadiusMinimumDefault&
       &,rotationCurveGradientHalfRadiusMaximum=rotationCurveGradientHalfRadiusMaximumDefault
  type(table1DLogarithmicLinear)                    :: rotationCurveGradientTable

  ! History of trial radii used to check for oscillations in the solution when solving for the structure of the disk.
  integer                                     :: radiusSolverIteration
  double precision                            :: radiusHistory(2)
  !$omp threadprivate(radiusHistory,radiusSolverIteration)

  ! The largest and smallest angular momentum, in units of that of a circular orbit at the virial radius, considered to be physically plausible for a disk. 
  double precision, parameter                 :: angularMomentumMaximum=1.0d1
  double precision, parameter                 :: angularMomentumMinimum=1.0d-6

  ! The radius (in units of the disk scale length) beyond which the disk is treated as a point mass for the purposes of computing
  ! rotation curves.
  double precision, parameter                 :: fractionalRadiusMaximum=30.0d0
  
  ! Record of whether this module has been initialized.
  logical                                     :: moduleInitialized   =.false.
  logical                                     :: threadAllocationDone=.false.
  !$omp threadprivate(threadAllocationDone)

contains
  
  !# <mergerTreePreTreeConstructionTask>
  !#  <unitName>Node_Component_Disk_Exponential_Initialize</unitName>
  !# </mergerTreePreTreeConstructionTask>
  subroutine Node_Component_Disk_Exponential_Initialize()
    !% Initializes the tree node exponential disk methods module.
    use Input_Parameters
    use Abundances_Structure
    use Stellar_Population_Properties_Luminosities
    use Memory_Management
    implicit none
    type(nodeComponentDiskExponential) :: diskExponentialComponent
    
     ! Initialize the module if necessary.
     !$omp critical (Node_Component_Disk_Exponential_Initialize)
     if (defaultDiskComponent%exponentialIsActive().and..not.moduleInitialized) then
        
        ! Get number of abundance properties.
        abundancesCount  =Abundances_Property_Count            ()
        
        ! Get number of luminosity properties.
        luminositiesCount=Stellar_Population_Luminosities_Count()

        ! Attach the cooling mass/angular momentum pipes from the hot halo component.
        call diskExponentialComponent%attachPipes()

        ! Bind the star formation rate function.
        call diskExponentialComponent%starFormationRateFunction(Node_Component_Disk_Exponential_Star_Formation_Rate)

        ! Read parameters controlling the physical implementation.
        !@ <inputParameter>
        !@   <name>diskMassToleranceAbsolute</name>
        !@   <defaultValue>$10^{-6} M_\odot$</defaultValue>
        !@   <attachedTo>module</attachedTo>
        !@   <description>
        !@    The mass tolerance used to judge whether the disk is physically plausible.
        !@   </description>
        !@   <type>real</type>
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
        !@   <type>real</type>
        !@   <cardinality>1</cardinality>
        !@ </inputParameter>
        call Get_Input_Parameter('diskOutflowTimescaleMinimum',diskOutflowTimescaleMinimum,defaultValue=1.0d-3)
        !@ <inputParameter>
        !@   <name>diskStructureSolverRadius</name>
        !@   <defaultValue>1</defaultValue>
        !@   <attachedTo>module</attachedTo>
        !@   <description>
        !@    The radius (in units of the exponential scale length) to use in solving for the size of the disk.
        !@   </description>
        !@   <type>real</type>
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
        !@     The ratio of scale height to scale radius for exponential disks.
        !@   </description>
        !@   <type>real</type>
        !@   <cardinality>1</cardinality>
        !@ </inputParameter>
        call Get_Input_Parameter('heightToRadialScaleDisk',heightToRadialScaleDisk,defaultValue=0.137d0)

        ! Compute the specific angular momentum of the disk at this structure solver radius in units of the mean specific angular
        ! momentum of the disk assuming a flat rotation curve.
        diskStructureSolverSpecificAngularMomentum=diskStructureSolverRadius/2.0d0

        ! If necessary, compute the specific angular momentum correction factor to account for the difference between rotation
        ! curves for thin disk and a spherical mass distribution.
        if (diskRadiusSolverCole2000Method) diskRadiusSolverFlatVsSphericalFactor=2.0d0*diskStructureSolverRadius&
             &*Node_Component_Disk_Exponential_Rotation_Curve_Bessel_Factors(0.5d0*diskStructureSolverRadius)&
             &-Node_Component_Disk_Exponential_Enclosed_Mass_Dimensionless  (      diskStructureSolverRadius)

        ! Record that the module is now initialized.
        moduleInitialized=.true.
     end if
     !$omp end critical (Node_Component_Disk_Exponential_Initialize)

     ! Allocate work arrays for luminosities for this thread.
     if (.not.threadAllocationDone) then
        call Alloc_Array(luminositiesDisk        ,[luminositiesCount])
        call Alloc_Array(luminositiesTransferRate,[luminositiesCount])
        call Alloc_Array(luminositiesStellarRates,[luminositiesCount])
        call Alloc_Array(zeroLuminosities        ,[luminositiesCount])
        call Alloc_Array(luminositiesMinimum     ,[luminositiesCount])
        zeroLuminosities   =0.0d0
        luminositiesMinimum=1.0d0
        threadAllocationDone=.true.
     end if     
     return
   end subroutine Node_Component_Disk_Exponential_Initialize
        
  !# <preEvolveTask>
  !# <unitName>Node_Component_Disk_Exponential_Pre_Evolve</unitName>
  !# </preEvolveTask>
  subroutine Node_Component_Disk_Exponential_Pre_Evolve(thisNode)
    !% Ensure the disk has been initialized.
    implicit none
    type (treeNode         ), pointer, intent(inout) :: thisNode
    class(nodeComponentDisk), pointer                :: thisDiskComponent

    ! Get the disk component.
    thisDiskComponent => thisNode%disk()
    ! Check if an exponential disk component exists.
    select type (thisDiskComponent)
    class is (nodeComponentDiskExponential)
       ! Initialize the disk
       call Node_Component_Disk_Exponential_Create(thisNode)
    end select
    return
  end subroutine Node_Component_Disk_Exponential_Pre_Evolve

  !# <postEvolveTask>
  !# <unitName>Node_Component_Disk_Exponential_Post_Evolve</unitName>
  !# </postEvolveTask>
  subroutine Node_Component_Disk_Exponential_Post_Evolve(thisNode)
    !% Trim histories attached to the disk.
    use Galacticus_Display
    use String_Handling
    use ISO_Varying_String
    use Abundances_Structure
    use Histories
    implicit none
    type     (treeNode          ), pointer, intent(inout) :: thisNode
    class    (nodeComponentDisk ), pointer                :: thisDiskComponent
    class    (nodeComponentBasic), pointer                :: thisBasicComponent
    double precision             , save                   :: fractionalErrorMaximum=0.0d0
    double precision                                      :: specificAngularMomentum,fractionalError,diskMass
    character(len=20            )                         :: valueString
    type     (varying_string    )                         :: message
    type     (history           )                         :: stellarPropertiesHistory

    ! Get the disk component.
    thisDiskComponent => thisNode%disk()
    ! Check if an exponential disk component exists.
    select type (thisDiskComponent)
    class is (nodeComponentDiskExponential)
       ! Trim the stellar populations properties future history.
       thisBasicComponent => thisNode%basic()
       stellarPropertiesHistory=thisDiskComponent%stellarPropertiesHistory()
       call stellarPropertiesHistory%trim(thisBasicComponent%time())
       call thisDiskComponent%stellarPropertiesHistorySet(stellarPropertiesHistory)

       ! Trap negative gas masses.
       if (thisDiskComponent%massGas() < 0.0d0) then
          ! Check if this exceeds the maximum previously recorded error.
          fractionalError=   abs(thisDiskComponent%massGas    ()) &
               &          /(                                      &
               &             abs(thisDiskComponent%massGas    ()) &
               &            +abs(thisDiskComponent%massStellar()) &
               &           )
          !$omp critical (Exponential_Disk_Post_Evolve_Check)
          if (fractionalError > fractionalErrorMaximum) then
             ! Report a warning.          
             message='Warning: disk has negative gas mass (fractional error exceeds any previously reported):'//char(10)
             message=message//'  Node index        = '//thisNode%index() //char(10)
             write (valueString,'(e12.6)') thisDiskComponent%massGas    ()
             message=message//'  Disk gas mass     = '//trim(valueString)//char(10)
             write (valueString,'(e12.6)') thisDiskComponent%massStellar()
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
          !$omp end critical (Exponential_Disk_Post_Evolve_Check)
          
          ! Get the specific angular momentum of the disk material
          diskMass= thisDiskComponent%massGas    () &
               &   +thisDiskComponent%massStellar()
          if (diskMass == 0.0d0) then
             specificAngularMomentum=0.0d0
             call thisDiskComponent%        massStellarSet(           0.0d0)
             call thisDiskComponent%  abundancesStellarSet(  zeroAbundances)
             call thisDiskComponent%luminositiesStellarSet(zeroLuminosities)
          else
             specificAngularMomentum=thisDiskComponent%angularMomentum()/diskMass
          end if

          ! Reset the gas, abundances and angular momentum of the disk.
          call thisDiskComponent%        massGasSet(                                                  0.0d0)
          call thisDiskComponent%  abundancesGasSet(                                         zeroAbundances)
          call thisDiskComponent%angularMomentumSet(specificAngularMomentum*thisDiskComponent%massStellar())

       end if

    end select
    return
  end subroutine Node_Component_Disk_Exponential_Post_Evolve

  subroutine Node_Component_Disk_Exponential_Create(thisNode)
    !% Create properties in an exponential disk component.
    use Histories
    use Stellar_Population_Properties
    use Galacticus_Output_Star_Formation_Histories
    implicit none
    type   (treeNode         ), intent(inout), pointer :: thisNode
    class  (nodeComponentDisk),                pointer :: thisDiskComponent
    type   (history          )                         :: stellarPropertiesHistory,starFormationHistory
    logical                                            :: createStellarPropertiesHistory,createStarFormationHistory

    ! Get the disk component.
    thisDiskComponent => thisNode%disk()
    ! Exit if already initialized.
    if (thisDiskComponent%isInitialized()) return
    ! Determine which histories must be created.
    starFormationHistory          =thisDiskComponent%starFormationHistory            ()
    createStarFormationHistory    =.not.             starFormationHistory    %exists ()
    call                                             starformationhistory    %destroy()
    stellarPropertiesHistory      =thisDiskComponent%stellarPropertiesHistory        ()
    createStellarPropertiesHistory=.not.             stellarPropertiesHistory%exists ()
    call                                             stellarPropertiesHistory%destroy()
    ! Create the stellar properties history.
    if (createStellarPropertiesHistory) then
       ! Create the stellar properties history.
       call Stellar_Population_Properties_History_Create (thisNode,stellarPropertiesHistory)
       call thisDiskComponent%stellarPropertiesHistorySet(         stellarPropertiesHistory)
    end if
    ! Create the star formation history.
    if (createStarFormationHistory    ) then
       call Star_Formation_History_Create                (thisNode,    starFormationHistory)
       call thisDiskComponent%    starFormationHistorySet(             starFormationHistory)
    end if
    ! Record that the disk has been initialized.
    call thisDiskComponent%isInitializedSet(.true.)
    return
  end subroutine Node_Component_Disk_Exponential_Create
  
  !# <rateComputeTask>
  !#  <unitName>Node_Component_Disk_Exponential_Rate_Compute</unitName>
  !# </rateComputeTask>
   subroutine Node_Component_Disk_Exponential_Rate_Compute(thisNode,interrupt,interruptProcedureReturn)
     !% Compute the exponential disk node mass rate of change.
     use Abundances_Structure
     use Cosmological_Parameters
     use Histories
     use Cooling_Rates
     use Star_Formation_Feedback_Disks
     use Star_Formation_Feedback_Expulsion_Disks
     use Galactic_Structure_Options
     use Galactic_Dynamics_Bar_Instabilities
     use Galacticus_Output_Star_Formation_Histories 
     use Stellar_Population_Properties
     use Numerical_Constants_Astronomical
     implicit none
     type     (treeNode             ), pointer, intent(inout) :: thisNode
     class    (nodeComponentDisk    ), pointer                :: thisDiskComponent
     class    (nodeComponentSpheroid), pointer                :: thisSpheroidComponent
     class    (nodeComponentHotHalo ), pointer                :: thisHotHaloComponent
     logical                         ,          intent(inout) :: interrupt
     procedure(                     ), pointer, intent(inout) :: interruptProcedureReturn
     procedure(                     ), pointer                :: interruptProcedure
     type     (abundances           ), save                   :: fuelAbundances,stellarAbundancesRates,fuelAbundancesRates
     !$omp threadprivate(fuelAbundances,stellarAbundancesRates,fuelAbundancesRates)
     double precision                                         :: starFormationRate,stellarMassRate,fuelMassRate,fuelMass,massOutflowRate&
          &,diskMass,angularMomentumOutflowRate,transferRate,barInstabilityTimescale,gasMass,energyInputRate,diskDynamicalTime&
          &,massOutflowRateToHotHalo,massOutflowRateFromHalo,outflowToHotHaloFraction,angularMomentum
     type     (history              )                         :: historyTransferRate,stellarHistoryRate

    ! Get a local copy of the interrupt procedure.
    interruptProcedure => interruptProcedureReturn

    ! Get the disk and check that it is of our class.
    thisDiskComponent => thisNode%disk()
    select type (thisDiskComponent)
    class is (nodeComponentDiskExponential)

       ! Check for a realistic disk, return immediately if disk is unphysical.
       if     (     thisDiskComponent%angularMomentum() < 0.0d0 &
            &  .or. thisDiskComponent%radius         () < 0.0d0 &
            &  .or. thisDiskComponent%massGas        () < 0.0d0 &
            & ) return

       ! Interrupt if the disk is not initialized.
       if (.not.thisDiskComponent%isInitialized()) then
          interrupt=.true.
          interruptProcedureReturn => Node_Component_Disk_Exponential_Create
          return
       end if

       ! Compute the star formation rate.
       starFormationRate=thisDiskComponent%starFormationRate()
       
       ! Get the available fuel mass.
       fuelMass         =thisDiskComponent%massGas          ()
       
       ! Find the metallicity of the fuel supply.
       fuelAbundances   =thisDiskComponent%abundancesGas    ()
       call fuelAbundances%massToMassFraction(fuelMass)
       
       ! Find rates of change of stellar mass, gas mass, abundances and luminosities.
       stellarHistoryRate=thisDiskComponent%stellarPropertiesHistory()
       call Stellar_Population_Properties_Rates(starFormationRate,fuelAbundances,componentTypeDisk,thisNode,stellarHistoryRate&
            &,stellarMassRate,stellarAbundancesRates,luminositiesStellarRates,fuelMassRate,fuelAbundancesRates,energyInputRate)
       if (stellarHistoryRate%exists()) call thisDiskComponent%stellarPropertiesHistoryRate(stellarHistoryRate)

       ! Adjust rates.
       call thisDiskComponent%        massStellarRate(         stellarMassRate)
       call thisDiskComponent%            massGasRate(            fuelMassRate)
       call thisDiskComponent%  abundancesStellarRate(  stellarAbundancesRates)
       call thisDiskComponent%      abundancesGasRate(     fuelAbundancesRates)
       call thisDiskComponent%luminositiesStellarRate(luminositiesStellarRates)

       ! Record the star formation history.
       stellarHistoryRate=thisDiskComponent%starFormationHistory()
       call stellarHistoryRate%reset()
       call Star_Formation_History_Record(thisNode,stellarHistoryRate,fuelAbundances,starFormationRate)
       if (stellarHistoryRate%exists()) call thisDiskComponent%starFormationHistoryRate(stellarHistoryRate)

       ! Find rate of outflow of material from the disk and pipe it to the outflowed reservoir.
       massOutflowRateToHotHalo=Star_Formation_Feedback_Disk_Outflow_Rate          (thisNode,starFormationRate,energyInputRate)
       massOutflowRateFromHalo =Star_Formation_Expulsive_Feedback_Disk_Outflow_Rate(thisNode,starFormationRate,energyInputRate)
       massOutflowRate         =massOutflowRateToHotHalo+massOutflowRateFromHalo
       if (massOutflowRate > 0.0d0) then
          ! Find the fraction of material which outflows to the hot halo.
          outflowToHotHaloFraction=massOutflowRateToHotHalo/massOutflowRate

          ! Get the masses of the disk.
          gasMass =        thisDiskComponent%massGas    ()
          diskMass=gasMass+thisDiskComponent%massStellar()
          
          ! Limit the outflow rate timescale to a multiple of the dynamical time.
          diskDynamicalTime=Mpc_per_km_per_s_To_Gyr*thisDiskComponent%radius()/thisDiskComponent%velocity()   
          massOutflowRate=min(massOutflowRate,gasMass/diskOutflowTimescaleMinimum/diskDynamicalTime)
          
          ! Compute the angular momentum outflow rate.
          if (diskMass > 0.0d0) then
             angularMomentum           =thisDiskComponent%angularMomentum()
             angularMomentumOutflowRate=angularMomentum*(massOutflowRate/diskMass)
             angularMomentumOutflowRate=sign(min(abs(angularMomentumOutflowRate),abs(angularMomentum/diskOutflowTimescaleMinimum &
                  &/diskDynamicalTime)),angularMomentumOutflowRate)
          else
             angularMomentumOutflowRate=0.0d0
          end if
          if (gasMass > 0.0d0) then
             fuelAbundancesRates=thisDiskComponent%abundancesGas()
             call fuelAbundancesRates%massToMassFraction(gasMass)
             fuelAbundancesRates=fuelAbundancesRates*massOutflowRate
          else
             fuelAbundancesRates=zeroAbundances
          end if
          thisHotHaloComponent => thisNode%hotHalo()
          call thisHotHaloComponent%           outflowingMassRate( massOutflowRate           *outflowToHotHaloFraction)
          call thisDiskComponent   %                  massGasRate(-massOutflowRate                                    )
          call thisHotHaloComponent%outflowingAngularMomentumRate( angularMomentumOutflowRate*outflowToHotHaloFraction)
          call thisDiskComponent   %          angularMomentumRate(-angularMomentumOutflowRate                         )
          call thisHotHaloComponent%     outflowingAbundancesRate( fuelAbundancesRates       *outflowToHotHaloFraction)
          call thisDiskComponent   %            abundancesGasRate(-fuelAbundancesRates                                )
       end if

       ! Determine if the disk is bar unstable and, if so, the rate at which material is moved to the pseudo-bulge.
       if (thisNode%isPhysicallyPlausible) then
          ! Disk has positive angular momentum, so compute an instability timescale.
          barInstabilityTimescale=Bar_Instability_Timescale(thisNode)
       else
          ! Disk has non-positive angular momentum, therefore it is unphysical. Do not compute an instability timescale in this
          ! case as the disk radius may be unphysical also.
          barInstabilityTimescale=-1.0d0
       end if

       ! Negative timescale indicates no bar instability.
       if (barInstabilityTimescale >= 0.0d0) then
          ! Disk is unstable, so compute rates at which material is transferred to the spheroid.
          thisSpheroidComponent => thisNode%spheroid()
          ! Gas mass.
          transferRate          =max(         0.0d0,thisDiskComponent    %massGas                (                         ))/barInstabilityTimescale
          call                                      thisDiskComponent    %massGasRate            (-           transferRate                              )
          call                                      thisSpheroidComponent%massGasRate            (+           transferRate ,interrupt,interruptProcedure)
          ! Stellar mass.
          transferRate          =max(         0.0d0,thisDiskComponent    %massStellar            (                         ))/barInstabilityTimescale
          call                                      thisDiskComponent    %massStellarRate        (-           transferRate                              )
          call                                      thisSpheroidComponent%massStellarRate        (+           transferRate ,interrupt,interruptProcedure)
          ! Angular momentum.
          transferRate          =max(         0.0d0,thisDiskComponent    %angularMomentum        (                         ))/barInstabilityTimescale
          call                                      thisDiskComponent    %angularMomentumRate    (-           transferRate                              )
          call                                      thisSpheroidComponent%angularMomentumRate    (+           transferRate ,interrupt,interruptProcedure)
          ! Gas abundances.
          fuelAbundancesRates   =max(zeroAbundances,thisDiskComponent    %abundancesGas          (                         ))/barInstabilityTimescale
          call                                      thisDiskComponent    %abundancesGasRate      (-     fuelAbundancesRates                             )
          call                                      thisSpheroidComponent%abundancesGasRate      (+     fuelAbundancesRates,interrupt,interruptProcedure)
          ! Stellar abundances.
          stellarAbundancesRates=max(zeroAbundances,thisDiskComponent    %abundancesStellar      (                         ))/barInstabilityTimescale
          call                                      thisDiskComponent    %abundancesStellarRate  (-  stellarAbundancesRates                             )
          call                                      thisSpheroidComponent%abundancesStellarRate  (+  stellarAbundancesRates,interrupt,interruptProcedure)
          ! Stellar luminosities.
          luminositiesTransferRate=max(       0.0d0,thisDiskComponent    %luminositiesStellar    (                         ))/barInstabilityTimescale
          call                                      thisDiskComponent    %luminositiesStellarRate(-luminositiesTransferRate                             )
          call                                      thisSpheroidComponent%luminositiesStellarRate(+luminositiesTransferRate,interrupt,interruptProcedure)
          ! Stellar properties history.
          historyTransferRate=thisDiskComponent%stellarPropertiesHistory()
          if (historyTransferRate%exists()) then
             historyTransferRate=historyTransferRate/barInstabilityTimescale
             call thisDiskComponent    %stellarPropertiesHistoryRate(-historyTransferRate                             )
             call thisSpheroidComponent%stellarPropertiesHistoryRate(+historyTransferRate,interrupt,interruptProcedure)
          end if
          call historyTransferRate%destroy()
          ! Star formation history.
          historyTransferRate=thisDiskComponent%starFormationHistory()
          if (historyTransferRate%exists()) then
             historyTransferRate=historyTransferRate/barInstabilityTimescale
             call thisDiskComponent    %starFormationHistoryRate(-historyTransferRate                             )
             call thisSpheroidComponent%starFormationHistoryRate(+historyTransferRate,interrupt,interruptProcedure)
          end if

       end if

    end select

    ! Return the procedure pointer.
    interruptProcedureReturn => interruptProcedure
    
     return
   end subroutine Node_Component_Disk_Exponential_Rate_Compute

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Disk_Exponential_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Disk_Exponential_Scale_Set(thisNode)
    !% Set scales for properties of {\tt thisNode}.
    use Histories
    use Stellar_Population_Properties
    use Galacticus_Output_Star_Formation_Histories
    use Abundances_Structure
    implicit none
    type (treeNode             ), pointer,  intent(inout) :: thisNode
    class(nodeComponentDisk    ), pointer                 :: thisDiskComponent
    class(nodeComponentSpheroid), pointer                 :: thisSpheroidComponent
    double precision            , parameter               :: massMinimum           =1.0d0
    double precision            , parameter               :: angularMomentumMinimum=1.0d-2
    double precision                                      :: mass,angularMomentum
    type (history              )                          :: stellarPopulationHistoryScales

    ! Get the disk component.
    thisDiskComponent => thisNode%disk()
    ! Check if an exponential disk component exists.
    select type (thisDiskComponent)
    class is (nodeComponentDiskExponential)
       
       ! Get disk components.
       thisSpheroidComponent => thisNode%spheroid()

       ! Set scale for angular momentum.
       angularMomentum=thisDiskComponent%angularMomentum()+thisSpheroidComponent%angularMomentum()
       call thisDiskComponent%angularMomentumScale(max(angularMomentum,angularMomentumMinimum))

       ! Set scale for masses.
       mass           = thisDiskComponent%massGas       ()+thisSpheroidComponent%massGas        () &
            &          +thisDiskComponent%massStellar   ()+thisSpheroidComponent%massStellar    ()
       call thisDiskComponent%massGasScale              (max(           mass,           massMinimum))
       call thisDiskComponent%massStellarScale          (max(           mass,           massMinimum))

       ! Set scales for abundances if necessary.
       if (abundancesCount > 0) then
          ! Set scale for gas abundances.
          call thisDiskComponent%abundancesGasScale      (                                                  &
               &                                          max(                                              &
               &                                               thisDiskComponent    %abundancesGas      ()  &
               &                                              +thisSpheroidComponent%abundancesGas      (), &
               &                                               massMinimum*unitAbundances                   &
               &                                             )                                              &
               &                                         )
          
          ! Set scale for stellar abundances.
          call thisDiskComponent%abundancesStellarScale  (                                                  &
               &                                          max(                                              &
               &                                               thisDiskComponent    %abundancesStellar  ()  &
               &                                              +thisSpheroidComponent%abundancesStellar  (), &
               &                                               massMinimum*unitAbundances                   &
               &                                             )                                              &
               &                                         )
       end if

       ! Set scales for stellar luminosities if necessary.
       if (luminositiesCount > 0) then        
          ! Set scale for stellar luminosities.
          call thisDiskComponent%luminositiesStellarScale(                                                  &
               &                                          max(                                              &
               &                                               thisDiskComponent    %luminositiesStellar()  &
               &                                              +thisSpheroidComponent%luminositiesStellar(), &
               &                                               luminositiesMinimum                          &
               &                                             )                                              &
               &                                         )
       end if

       ! Set scales for stellar population properties and star formation histories.
       stellarPopulationHistoryScales=thisDiskComponent%stellarPropertiesHistory()
       call Stellar_Population_Properties_Scales           (stellarPopulationHistoryScales,thisDiskComponent%massStellar(),thisDiskComponent%abundancesStellar())
       call thisDiskComponent%stellarPropertiesHistoryScale(stellarPopulationHistoryScales                                                                      )
       call stellarPopulationHistoryScales%destroy()
       stellarPopulationHistoryScales=thisDiskComponent%starFormationHistory()
       call Star_Formation_History_Scales                  (stellarPopulationHistoryScales,thisDiskComponent%massStellar(),thisDiskComponent%abundancesStellar())
       call thisDiskComponent%starFormationHistoryScale    (stellarPopulationHistoryScales                                                                      )
       call stellarPopulationHistoryScales%destroy()

    end select
    return
  end subroutine Node_Component_Disk_Exponential_Scale_Set

   !# <satelliteMergerTask>
   !#  <unitName>Node_Component_Disk_Exponential_Satellite_Merging</unitName>
   !#  <after>Satellite_Merging_Mass_Movement_Store</after>
   !#  <after>Satellite_Merging_Remnant_Size</after>
   !# </satelliteMergerTask>
   subroutine Node_Component_Disk_Exponential_Satellite_Merging(thisNode)
     !% Transfer any exponential disk associated with {\tt thisNode} to its host halo.
     use Histories
     use Abundances_Structure
     use Satellite_Merging_Mass_Movements_Descriptors
     use Galacticus_Error
     implicit none
     type (treeNode             ), pointer, intent(inout)       :: thisNode
     class(nodeComponentDisk    ), pointer                      :: thisDiskComponent    ,hostDiskComponent
     class(nodeComponentSpheroid), pointer                      :: thisSpheroidComponent,hostSpheroidComponent
     type (treeNode             ), pointer                      :: hostNode
     double precision            , dimension(luminositiesCount) :: thisLuminosities,hostLuminosities
     type (history              )                               :: thisHistory,hostHistory
     double precision                                           :: specificAngularMomentum

     ! Check that the disk is of the exponential class.
     thisDiskComponent => thisNode%disk()
     select type (thisDiskComponent)
     class is (nodeComponentDiskExponential)
        thisSpheroidComponent => thisNode%spheroid()

        ! Find the node to merge with.
        hostNode              => thisNode%mergesWith(                 )
        hostDiskComponent     => hostNode%disk      (autoCreate=.true.)      
        hostSpheroidComponent => hostNode%spheroid  (autoCreate=.true.)      
 
        ! Get specific angular momentum of the disk material.
        if (                                                            thisDiskComponent%massGas()+thisDiskComponent%massStellar() > 0.0d0) then
           specificAngularMomentum=thisDiskComponent%angularMomentum()/(thisDiskComponent%massGas()+thisDiskComponent%massStellar())
        else
           specificAngularMomentum=0.0d0
        end if

        ! Move the gas component of the exponential disk to the host.
        select case (thisMergerGasMovesTo)
        case (movesToDisk)
           call hostDiskComponent    %massGasSet            (                                                                     &
                &                                             hostDiskComponent    %massGas            ()                         &
                &                                            +thisDiskComponent    %massGas            ()                         &
                &                                           )
           call hostDiskComponent    %abundancesGasSet      (                                                                     &
                &                                             hostDiskComponent    %abundancesGas      ()                         &
                &                                            +thisDiskComponent    %abundancesGas      ()                         &
                &                                           )
           call hostDiskComponent    %angularMomentumSet    (                                                                     &
                &                                             hostDiskComponent    %angularMomentum    ()                         &
                &                                            +thisDiskComponent    %massGas            ()*specificAngularMomentum &
                &                                           )
        case (movesToSpheroid)    
           call hostSpheroidComponent%massGasSet            (                                                                     &
                &                                             hostSpheroidComponent%massGas            ()                         &
                &                                            +thisDiskComponent    %massGas            ()                         &
                &                                           )
           call hostSpheroidComponent%abundancesGasSet      (                                                                     &
                &                                             hostSpheroidComponent%abundancesGas      ()                         &
                &                                            +thisDiskComponent    %abundancesGas      ()                         &
                &                                           )
        case default
           call Galacticus_Error_Report('Node_Component_Disk_Exponential_Satellite_Merging','unrecognized movesTo descriptor')
        end select
        call thisDiskComponent%      massGasSet(         0.0d0)
        call thisDiskComponent%abundancesGasSet(zeroAbundances)

        ! Move the stellar component of the exponential disk to the host.
        select case (thisMergerStarsMoveTo)
        case (movesToDisk)
           call hostDiskComponent    %massStellarSet        (                                                                     &
                &                                             hostDiskComponent    %massStellar        ()                         &
                &                                            +thisDiskComponent    %massStellar        ()                         &
                &                                           )
           call hostDiskComponent    %abundancesStellarSet  (                                                                     &
                &                                             hostDiskComponent    %abundancesStellar  ()                         &
                &                                            +thisDiskComponent    %abundancesStellar  ()                         &
                &                                           )
           call hostDiskComponent    %luminositiesStellarSet(                                                                     &
                &                                             hostDiskComponent    %luminositiesStellar()                         &
                &                                            +thisDiskComponent    %luminositiesStellar()                         &
                &                                           )
           call hostDiskComponent    %angularMomentumSet    (                                                                     &
                &                                             hostDiskComponent    %angularMomentum    ()                         &
                &                                            +thisDiskComponent    %massStellar        ()*specificAngularMomentum &
                &                                           )
           ! Also add stellar properties histories.
           thisHistory=thisDiskComponent    %stellarPropertiesHistory()
           hostHistory=hostDiskComponent    %stellarPropertiesHistory()
           call hostHistory%addRates(thisHistory)
           call thisHistory%reset   (           )
           call hostDiskComponent%stellarPropertiesHistorySet(hostHistory)
           call thisDiskComponent%stellarPropertiesHistorySet(thisHistory)
           ! Also add star formation histories.
           thisHistory=thisDiskComponent    %starFormationHistory    ()
           hostHistory=hostDiskComponent    %starFormationHistory    ()
           call hostHistory%combine (thisHistory)
           call thisHistory%reset   (           )
           call hostDiskComponent    %starFormationHistorySet    (hostHistory)
           call thisDiskComponent    %starFormationHistorySet    (thisHistory)
           call thisHistory%destroy(recordMemory=.false.)
           call hostHistory%destroy(recordMemory=.false.)
        case (movesToSpheroid)
           call hostSpheroidComponent%massStellarSet        (                                                                     &
                &                                             hostSpheroidComponent%massStellar        ()                         &
                &                                            +thisDiskComponent    %massStellar        ()                         &
                &                                           )
           call hostSpheroidComponent%abundancesStellarSet  (                                                                     &
                &                                             hostSpheroidComponent%abundancesStellar  ()                         &
                &                                            +thisDiskComponent    %abundancesStellar  ()                         &
                &                                            )
           call hostSpheroidComponent%luminositiesStellarSet(                                                                     &
                &                                             hostSpheroidComponent%luminositiesStellar()                         &
                &                                            +thisDiskComponent    %luminositiesStellar()                         &
                &                                           )
           ! Also add stellar properties histories.
           thisHistory=thisDiskComponent    %stellarPropertiesHistory()
           hostHistory=hostSpheroidComponent%stellarPropertiesHistory()
           call hostHistory%addRates(thisHistory)
           call thisHistory%reset   (           )
           call hostSpheroidComponent%stellarPropertiesHistorySet(hostHistory)
           call thisDiskComponent    %stellarPropertiesHistorySet(thisHistory)
           ! Also add star formation histories.
           thisHistory=thisDiskComponent    %starFormationHistory    ()
           hostHistory=hostSpheroidComponent%starFormationHistory    ()
           call hostHistory%combine(thisHistory)
           call thisHistory%reset  (           )
           call hostSpheroidComponent%starFormationHistorySet(hostHistory)
           call thisDiskComponent    %starFormationHistorySet(thisHistory)
           call thisHistory%destroy(recordMemory=.false.)
           call hostHistory%destroy(recordMemory=.false.)
        case default
           call Galacticus_Error_Report('Node_Component_Disk_Exponential_Satellite_Merging','unrecognized movesTo descriptor')
        end select
        call thisDiskComponent%        massStellarSet(           0.0d0)
        call thisDiskComponent%  abundancesStellarSet(  zeroAbundances)
        call thisDiskComponent%luminositiesStellarSet(zeroLuminosities)
        call thisDiskComponent%    angularMomentumSet(           0.0d0)
     end select
     return
   end subroutine Node_Component_Disk_Exponential_Satellite_Merging
  
  !# <enclosedMassTask>
  !#  <unitName>Node_Component_Disk_Exponential_Enclosed_Mass</unitName>
  !# </enclosedMassTask>
  subroutine Node_Component_Disk_Exponential_Enclosed_Mass(thisNode,radius,massType,componentType,weightBy,weightIndex,componentMass)
    !% Computes the mass within a given radius for an exponential disk.
    use Galactic_Structure_Options
    implicit none
    type (treeNode         ), intent(inout), pointer :: thisNode
    integer                 , intent(in   )          :: massType,componentType,weightBy,weightIndex
    double precision        , intent(in   )          :: radius
    double precision        , intent(  out)          :: componentMass
    class(nodeComponentDisk),                pointer :: thisDiskComponent
    double precision                                 :: fractionalRadius,diskRadius
    
    ! Return immediately if disk component is not requested.    
    componentMass=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDisk)) return

    ! Get the disk component and check that it is of the exponential class.
    thisDiskComponent => thisNode%disk()
    select type (thisDiskComponent)
    class is (nodeComponentDiskExponential)
    
       select case (weightBy)
       case (weightByMass      )
          select case (massType)
          case (massTypeAll,massTypeBaryonic,massTypeGalactic)
             componentMass=thisDiskComponent%massGas()+thisDiskComponent%massStellar()
          case (massTypeGaseous)
             componentMass=thisDiskComponent%massGas()
          case (massTypeStellar)
             componentMass=                            thisDiskComponent%massStellar()
          end select
       case (weightByLuminosity)
          select case (massType)
          case (massTypeAll,massTypeBaryonic,massTypeGalactic,massTypeStellar)
             luminositiesDisk=thisDiskComponent%luminositiesStellar()
             componentMass   =luminositiesDisk(weightIndex)
          end select
       end select
       ! Return if no mass.
       if (componentMass <=       0.0d0) return
       ! Return if the total mass was requested.
       if (radius        >= radiusLarge) return
       ! Compute the actual mass.
       diskRadius=thisDiskComponent%radius()
       if (diskRadius > 0.0d0) then
          fractionalRadius=radius/diskRadius
          componentMass=componentMass*Node_Component_Disk_Exponential_Enclosed_Mass_Dimensionless(fractionalRadius)
       end if
    end select
    return
  end subroutine Node_Component_Disk_Exponential_Enclosed_Mass

   double precision function Node_Component_Disk_Exponential_Enclosed_Mass_Dimensionless(radius)
     !% Returns the fractional mass enclosed within {\tt radius} in a dimensionless exponential disk.
     implicit none
     double precision, intent(in) :: radius
     
     Node_Component_Disk_Exponential_Enclosed_Mass_Dimensionless=1.0d0-(1.0d0+radius)*dexp(-radius)
     return
   end function Node_Component_Disk_Exponential_Enclosed_Mass_Dimensionless

  !# <rotationCurveTask>
  !#  <unitName>Node_Component_Disk_Exponential_Rotation_Curve</unitName>
  !# </rotationCurveTask>
  subroutine Node_Component_Disk_Exponential_Rotation_Curve(thisNode,radius,massType,componentType,componentVelocity)
    !% Computes the rotation curve at a given radius for an exponential disk.
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    implicit none
    type (treeNode         ), intent(inout), pointer :: thisNode
    integer                 , intent(in   )          :: massType,componentType
    double precision        , intent(in   )          :: radius
    double precision        , intent(  out)          :: componentVelocity
    class(nodeComponentDisk),                pointer :: thisDiskComponent
    double precision                                 :: fractionalRadius,fractionalRadiusFactor,diskRadius,componentMass,halfRadius

    ! Set to zero by default.
    componentVelocity=0.0d0

    ! Get the disk component and check that it is of the exponential class.
    thisDiskComponent => thisNode%disk()
    select type (thisDiskComponent)
    class is (nodeComponentDiskExponential)
           
       ! Get the mass of the disk.
       call Node_Component_Disk_Exponential_Enclosed_Mass(thisNode,radiusLarge,massType,componentType,weightByMass,weightIndexNull,componentMass)
       if (componentMass <= 0.0d0) return
       
       ! Compute the actual velocity.
       diskRadius=thisDiskComponent%radius()
       if (diskRadius > 0.0d0) then
          fractionalRadius=radius/diskRadius
          if (fractionalRadius > fractionalRadiusMaximum) then
             ! Beyond some maximum radius, approximate the disk as a spherical distribution to avoid evaluating Bessel functions for
             ! very large arguments.
             componentVelocity=dsqrt(gravitationalConstantGalacticus*componentMass/radius)
          else
             ! We are often called at precisely one scale length. Use pre-computed factors in that case.
             if (fractionalRadius == 1.0d0) then
                !$omp critical (ExponentialDiskFactorCompute)
                if (.not.scaleLengthFactorSet) then
                   halfRadius          =0.5d0
                   scaleLengthFactor   =Node_Component_Disk_Exponential_Rotation_Curve_Bessel_Factors(halfRadius)
                   scaleLengthFactorSet=.true.
                end if
                !$omp end critical (ExponentialDiskFactorCompute)
                fractionalRadiusFactor=scaleLengthFactor
             else
                halfRadius=0.5d0*fractionalRadius
                fractionalRadiusFactor=Node_Component_Disk_Exponential_Rotation_Curve_Bessel_Factors(halfRadius)
             end if
             componentVelocity=dsqrt(2.0d0*(gravitationalConstantGalacticus*componentMass/diskRadius)*fractionalRadiusFactor)
          end if
       end if
    end select
    return
  end subroutine Node_Component_Disk_Exponential_Rotation_Curve
   
   double precision function Node_Component_Disk_Exponential_Rotation_Curve_Bessel_Factors(halfRadius)
     !% Compute Bessel function factors appearing in the expression for an razor-thin exponential disk rotation curve.
     use Memory_Management
     use Numerical_Constants_Math
     use Numerical_Interpolation
     use Bessel_Functions
     implicit none
     double precision, intent(in) :: halfRadius
     double precision, parameter  :: halfRadiusSmall=1.0d-3
     integer                      :: iPoint
     double precision             :: x

     ! For small half-radii, use a series expansion for a more accurate result.
     if (halfRadius == 0.0d0) then
        Node_Component_Disk_Exponential_Rotation_Curve_Bessel_Factors=0.0d0
        return
     else if (halfRadius < halfRadiusSmall) then
        Node_Component_Disk_Exponential_Rotation_Curve_Bessel_Factors=(ln2-eulersConstant-0.5d0-log(halfRadius))*halfRadius**2
        return
     end if
     
     !$omp critical(Exponential_Disk_Rotation_Curve_Tabulate)
     if (.not.rotationCurveInitialized .or. halfRadius <  rotationCurveHalfRadiusMinimum .or. halfRadius > rotationCurveHalfRadiusMaximum) then
        ! Find the minimum and maximum half-radii to tabulate.
        rotationCurveHalfRadiusMinimum=min(rotationCurveHalfRadiusMinimum,0.5d0*halfRadius)
        rotationCurveHalfRadiusMaximum=max(rotationCurveHalfRadiusMaximum,2.0d0*halfRadius)
        
        ! Determine how many points to tabulate.
        rotationCurvePointsCount=int(dlog10(rotationCurveHalfRadiusMaximum/rotationCurveHalfRadiusMinimum)*dble(rotationCurvePointsPerDecade))+1
        
        ! Allocate table arrays.
        call rotationCurveTable%destroy()
        call rotationCurveTable%create(rotationCurveHalfRadiusMinimum,rotationCurveHalfRadiusMaximum,rotationCurvePointsCount)
        
        ! Compute Bessel factors.
        do iPoint=1,rotationCurvePointsCount
           x=rotationCurveTable%x(iPoint)
           call rotationCurveTable%populate(                                               &
                &                            x**2                                          &
                &                           *(                                             &
                &                              Bessel_Function_I0(x)*Bessel_Function_K0(x) &
                &                             -Bessel_Function_I1(x)*Bessel_Function_K1(x) &
                &                            ),                                            &
                &                           iPoint                                         &
                &                          )
        end do
        
        ! Flag that the rotation curve is now initialized.
        rotationCurveInitialized=.true.
     end if
     !$omp end critical(Exponential_Disk_Rotation_Curve_Tabulate)
     
     ! Interpolate in the tabulated function.
     !$omp critical(Exponential_Disk_Rotation_Curve_Interpolate)
Node_Component_Disk_Exponential_Rotation_Curve_Bessel_Factors=rotationCurveTable%interpolate(halfRadius)
     !$omp end critical(Exponential_Disk_Rotation_Curve_Interpolate)
     
     return
   end function Node_Component_Disk_Exponential_Rotation_Curve_Bessel_Factors

  !# <rotationCurveGradientTask>
  !#  <unitName>Node_Component_Disk_Exponential_Rotation_Curve_Gradient</unitName>
  !# </rotationCurveGradientTask>
  subroutine Node_Component_Disk_Exponential_Rotation_Curve_Gradient(thisNode,radius,massType,componentType,componentRotationCurveGradient)
    !% Computes the rotation curve gradient for an exponential disk.
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    use Numerical_Constants_Prefixes
    implicit none
    type (treeNode         ), intent(inout), pointer :: thisNode
    class(nodeComponentDisk),                pointer :: thisDiskComponent
    integer                 , intent(in   )          :: massType,componentType
    double precision        , intent(in   )          :: radius
    double precision        , intent(  out)          :: componentRotationCurveGradient
    double precision                                 :: componentMass,besselArgument,besselFactor

    ! Set to zero by default.
    componentRotationCurveGradient=0.0d0

    ! Return if radius is zero.
    if (radius <= 0.0d0) return
    
    ! Get the disk component and check that it is of the exponential class.
    thisDiskComponent => thisNode%disk()
    select type (thisDiskComponent)
    class is (nodeComponentDiskExponential)

       call Node_Component_Disk_Exponential_Enclosed_Mass(thisNode,radiusLarge,massType,componentType,weightByMass,weightIndexNull,componentMass)
       if (componentMass              <= 0.0d0) return
       if (thisDiskComponent%radius() <= 0.0d0) return
       besselArgument= radius                          &
            &         /2.0d0                           &
            &         /thisDiskComponent%radius()
    
       if (2.0d0*besselArgument > fractionalRadiusMaximum) then
          ! Beyond some maximum radius, approximate the disk as a point mass to avoid evaluating Bessel functions for
          ! very large arguments.
          componentRotationCurveGradient=-gravitationalConstantGalacticus*componentMass/radius**2
          return
       end if

       ! Checks for low radius and approximations.
       besselFactor=Node_Component_Disk_Exponential_Rttn_Crv_Grdnt_Bssl_Fctrs(besselArgument)
       componentRotationCurveGradient= gravitationalConstantGalacticus    &
            &                         *componentMass                      &
            &                         *besselFactor                       &
            &                         /thisDiskComponent%radius()**2 
    end select
    return
  end subroutine Node_Component_Disk_Exponential_Rotation_Curve_Gradient

  double precision function Node_Component_Disk_Exponential_Rttn_Crv_Grdnt_Bssl_Fctrs(halfRadius)
    !% Compute Bessel function factors appearing in the expression for a razor-thin exponential disk rotation curve gradient.
    use Memory_Management
    use Numerical_Constants_Math
    use Numerical_Interpolation
    use Bessel_Functions
    implicit none
    double precision, intent(in) :: halfRadius
    double precision, parameter  :: halfRadiusSmall=1.0d-3
    double precision, parameter  :: halfRadiusLarge=1.0d+2
    integer                      :: iPoint
    double precision             :: x

    ! For small and large half-radii, use a series expansion for a more accurate result.
    if (halfRadius == 0.0d0) then
       Node_Component_Disk_Exponential_Rttn_Crv_Grdnt_Bssl_Fctrs=0.0d0
       return
    else if (halfRadius < halfRadiusSmall) then
       Node_Component_Disk_Exponential_Rttn_Crv_Grdnt_Bssl_Fctrs=(ln2-log(halfRadius)-eulersConstant-1.0d0)*halfRadius &
            & +(1.5d0*ln2-1.5d0*log(halfRadius)+0.25d0-1.5d0*eulersConstant)*halfRadius**3
       return
    else if (halfRadius > halfRadiusLarge) then
       Node_Component_Disk_Exponential_Rttn_Crv_Grdnt_Bssl_Fctrs=-0.125d0-27.0d0/64.0d0/halfRadius**2
       return
    end if

    !$omp critical(Exponential_Disk_Rotation_Curve_Gradient_Tabulate)
    if (.not.rotationCurveGradientInitialized .or. halfRadius <  rotationCurveHalfRadiusMinimum .or. halfRadius > rotationCurveGradientHalfRadiusMaximum) then
       ! Find the minimum and maximum half-radii to tabulate.
       rotationCurveGradientHalfRadiusMinimum=min(rotationCurveGradientHalfRadiusMinimum,0.5d0*halfRadius)
       rotationCurveGradientHalfRadiusMaximum=max(rotationCurveGradientHalfRadiusMaximum,2.0d0*halfRadius)

       ! Determine how many points to tabulate.
       rotationCurveGradientPointsCount=int(dlog10(rotationCurveGradientHalfRadiusMaximum/rotationCurveGradientHalfRadiusMinimum)*dble(rotationCurveGradientPointsPerDecade))+1

       ! Allocate table arrays.
       call rotationCurveGradientTable%destroy()
       call rotationCurveGradientTable%create(rotationCurveGradientHalfRadiusMinimum,rotationCurveGradientHalfRadiusMaximum,rotationCurveGradientPointsCount)

       ! Compute Bessel factors.
       do iPoint=1,rotationCurveGradientPointsCount
          x=rotationCurveGradientTable%x(iPoint)
          call rotationCurveGradientTable%populate(                                                 &
               &                                    x**2                                            &
               &                                   *( x                                             &
               &                                     *  Bessel_Function_I0(x)*Bessel_Function_K0(x) &
               &                                     +x**2                                          &
               &                                     *( Bessel_Function_I1(x)*Bessel_Function_K0(x) &
               &                                       -Bessel_Function_I0(x)*Bessel_Function_K1(x) &
               &                                      )                                             &
               &                                    ),                                              &
               &                                   iPoint                                           &
               &                                  )
       end do
       
       ! Flag that the rotation curve is now initialized.
       rotationCurveGradientInitialized=.true.
    end if
    !$omp end critical(Exponential_Disk_Rotation_Curve_Gradient_Tabulate)
    
    ! Interpolate in the tabulated function.
    !$omp critical(Exponential_Disk_Rotation_Curve_Gradient_Interpolate)
    Node_Component_Disk_Exponential_Rttn_Crv_Grdnt_Bssl_Fctrs=rotationCurveGradientTable%interpolate(halfRadius)
    !$omp end critical(Exponential_Disk_Rotation_Curve_Gradient_Interpolate)
    return
  end function Node_Component_Disk_Exponential_Rttn_Crv_Grdnt_Bssl_Fctrs

  !# <potentialTask>
  !#  <unitName>Node_Component_Disk_Exponential_Potential</unitName>
  !# </potentialTask>
  subroutine Node_Component_Disk_Exponential_Potential(thisNode,radius,componentType,massType,componentPotential)
    !% Compute the gravitational potential due to an exponential disk.
    use Numerical_Constants_Physical 
    use Bessel_Functions
    use Galactic_Structure_Options
    implicit none
    type (treeNode         ), intent(inout), pointer :: thisNode
    class(nodeComponentDisk),                pointer :: thisDiskComponent
    integer                 , intent(in   )          :: componentType,massType
    double precision        , intent(in   )          :: radius
    double precision        , intent(  out)          :: componentPotential
    double precision                                 :: halfRadius,correctionSmallRadius,componentMass

    ! Return immediately if disk component is not requested.    
    componentPotential=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDisk)) return

    ! Avoid an arithmetic exception at radius zero.
    if (radius <= 0.0d0) return

    ! Get the disk component and check that it is of the exponential class.
    thisDiskComponent => thisNode%disk()
    select type (thisDiskComponent)
    class is (nodeComponentDiskExponential)
    
       ! Get the relevant mass of the disk.
       call Node_Component_Disk_Exponential_Enclosed_Mass(thisNode,radiusLarge,massType,componentType,weightByMass,weightIndexNull,componentMass)
       if (componentMass <= 0.0d0) return

       ! Compute the potential. If the radius is lower than the height then approximate the disk mass as being spherically distributed.
       if (radius > heightToRadialScaleDisk*thisDiskComponent%radius()) then
          halfRadius           =radius/2.0d0/thisDiskComponent%radius()
          correctionSmallRadius=0.0d0
       else
          halfRadius           =heightToRadialScaleDisk/2.0d0
          correctionSmallRadius=-Node_Component_Disk_Exponential_Enclosed_Mass_Dimensionless(heightToRadialScaleDisk)/heightToRadialScaleDisk/thisDiskComponent%radius()
          if (radius > 0.0d0) &
               & correctionSmallRadius=correctionSmallRadius+Node_Component_Disk_Exponential_Enclosed_Mass_Dimensionless(radius/thisDiskComponent%radius())/radius          
          correctionSmallRadius=-correctionSmallRadius*gravitationalConstantGalacticus*componentMass
       end if
       
       ! Compute the potential including the correction to small radii.
       componentPotential=-0.5d0                             &
            &             *gravitationalConstantGalacticus   &
            &             *componentMass                     &
            &             /thisDiskComponent%radius()        & 
            &             *(                                 &
            &                Bessel_Function_I0(halfRadius)  &
            &               *Bessel_Function_K1(halfRadius)  &
            &               -Bessel_Function_I1(halfRadius)  &
            &               *Bessel_Function_K0(halfRadius)  &
            &              )                                 &
            &             + correctionSmallRadius
    end select
    return
  end subroutine Node_Component_Disk_Exponential_Potential

  !# <densityTask>
  !#  <unitName>Node_Component_Disk_Exponential_Density</unitName>
  !# </densityTask>
  subroutine Node_Component_Disk_Exponential_Density(thisNode,positionSpherical,massType,componentType,componentDensity)
    !% Computes the density at a given position for an exponential disk.
    use Galactic_Structure_Options
    use Numerical_Constants_Math
    use Coordinate_Systems
    use Galacticus_Nodes
    implicit none
    type (treeNode         ), intent(inout), pointer :: thisNode
    integer                 , intent(in   )          :: massType,componentType
    double precision        , intent(in   )          :: positionSpherical(3)
    double precision        , intent(  out)          :: componentDensity
    class(nodeComponentDisk),                pointer :: thisDiskComponent
    double precision        , parameter              :: diskHeightToRadiusRatio=0.1d0
    double precision                                 :: fractionalRadius,fractionalHeight,positionCylindrical(3)

    ! Return immediately if disk component is not requested.    
    componentDensity=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDisk)) return

    ! Get the disk component and check that it is of the exponential class.
    thisDiskComponent => thisNode%disk()
    select type (thisDiskComponent)
    class is (nodeComponentDiskExponential)
    
       ! Determine mass type.
       select case (massType)
       case (massTypeAll,massTypeBaryonic,massTypeGalactic)
          componentDensity=thisDiskComponent%massGas()+thisDiskComponent%massStellar()
       case (massTypeGaseous)
          componentDensity=thisDiskComponent%massGas()
       case (massTypeStellar)
          componentDensity=                            thisDiskComponent%massStellar()
       end select
       ! Skip further calculation if mass or radius is zero.
       if (componentDensity > 0.0d0 .and. thisDiskComponent%radius() > 0.0d0) then
          ! Compute the actual density.
          positionCylindrical=Coordinates_Spherical_To_Cylindrical(positionSpherical)
          fractionalRadius=positionCylindrical(1)/                         thisDiskComponent%radius()
          fractionalHeight=positionCylindrical(2)/(diskHeightToRadiusRatio*thisDiskComponent%radius())
          componentDensity=componentDensity*exp(-fractionalRadius)/cosh(0.5d0*fractionalHeight)**2/4.0d0/Pi&
               &/thisDiskComponent%radius()**3/diskHeightToRadiusRatio
       end if
    end select
    return
  end subroutine Node_Component_Disk_Exponential_Density

  !# <radiusSolverPlausibility>
  !#  <unitName>Node_Component_Disk_Exponential_Radius_Solver_Plausibility</unitName>
  !# </radiusSolverPlausibility>
  subroutine Node_Component_Disk_Exponential_Radius_Solver_Plausibility(thisNode,galaxyIsPhysicallyPlausible)
     !% Determines whether the disk is physically plausible for radius solving tasks. Require that it have non-zero mass and angular momentum.
     use Dark_Matter_Halo_Scales
     implicit none
     type (treeNode         ), pointer, intent(inout) :: thisNode
     logical,                           intent(inout) :: galaxyIsPhysicallyPlausible
     class(nodeComponentDisk), pointer                :: thisDiskComponent
     double precision                                 :: angularMomentumScale

     ! Return immediately if our method is not selected.
     if (.not.defaultDiskComponent%exponentialIsActive()) return

     ! Determine the plausibility of the current disk.
     thisDiskComponent => thisNode%disk()
     select type (thisDiskComponent)
     class is (nodeComponentDiskExponential)
        if      (thisDiskComponent%massStellar()+thisDiskComponent%massGas() < -diskMassToleranceAbsolute) then
           galaxyIsPhysicallyPlausible=.false.
        else if (thisDiskComponent%massStellar()+thisDiskComponent%massGas() >=                     0.0d0) then
           if      (                                                                                                     &
                &   thisDiskComponent%angularMomentum() <  0.0d0                                                         &
                &  ) then
              galaxyIsPhysicallyPlausible=.false.
           else
              angularMomentumScale=(                                            &
                   &                 thisDiskComponent%massStellar()            &
                   &                +thisDiskComponent%massGas    ()            &
                   &               )                                            &
                   &               * Dark_Matter_Halo_Virial_Radius  (thisNode) &
                   &               * Dark_Matter_Halo_Virial_Velocity(thisNode)
              if     (                                                                                   &
                   &   thisDiskComponent%angularMomentum() > angularMomentumMaximum*angularMomentumScale &
                   &  .or.                                                                               &
                   &   thisDiskComponent%angularMomentum() < angularMomentumMinimum*angularMomentumScale &
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
   end subroutine Node_Component_Disk_Exponential_Radius_Solver_Plausibility
   
   double precision function Node_Component_Disk_Exponential_Radius_Solve(thisNode)
     !% Return the radius of the exponential disk used in structure solvers.
     implicit none
     type (treeNode         ), pointer, intent(inout) :: thisNode
     class(nodeComponentDisk), pointer                :: thisDiskComponent
     
     thisDiskComponent => thisNode%disk()
     Node_Component_Disk_Exponential_Radius_Solve=thisDiskComponent%radius()*diskStructureSolverRadius
     return
   end function Node_Component_Disk_Exponential_Radius_Solve
   
   subroutine Node_Component_Disk_Exponential_Radius_Solve_Set(thisNode,radius)
     !% Set the radius of the exponential disk used in structure solvers.
     implicit none
     type (treeNode         ), pointer, intent(inout) :: thisNode
     double precision        ,          intent(in   ) :: radius
     class(nodeComponentDisk), pointer                :: thisDiskComponent
     integer,                  parameter              :: iterationsForBisectionMinimum=10
     double precision                                 :: newRadius
     
     ! If using the Cole et al. (2000) method, check whether the solution is oscillating. This can happen as the effective
     ! angular momentum of the disk becomes radius dependent under this algorithm.
     newRadius=radius
     if (diskRadiusSolverCole2000Method) then
        if     (                                                                                                         &
             &             radiusSolverIteration                                        > iterationsForBisectionMinimum  &
             &  .and. all( radiusHistory                                                >= 0.0d0                       ) &
             &  .and.     (radiusHistory(2)-radiusHistory(1))*(radiusHistory(1)-radius) <  0.0d0                         &
             & ) then
           ! An oscillation has been detected - attempt to break out of it. The following heuristic has been found to work quite
           ! well - we bisect previous solutions in the oscillating sequence in a variety of different ways
           ! (arithmetic/geometric and using the current+previous or two previous solutions), alternating the bisection method
           ! sequentially. There's no guarantee that this will work in every situation however.
           select case (mod(radiusSolverIteration,4))
           case (0)
              newRadius=dsqrt (radius          *radiusHistory(1))
           case (1)
              newRadius=0.5d0*(radius          +radiusHistory(1))
           case (2)
              newRadius=dsqrt (radiusHistory(1)*radiusHistory(2))
           case (3)
              newRadius=0.5d0*(radiusHistory(1)+radiusHistory(2))
           end select
           radiusHistory=-1.0d0
        end if
        radiusSolverIteration=radiusSolverIteration+1
        radiusHistory(2)     =radiusHistory(1)
        radiusHistory(1)     =newRadius
     end if
     thisDiskComponent => thisNode%disk()
     call thisDiskComponent%radiusSet(max(newRadius,0.0d0)/diskStructureSolverRadius)
     return
   end subroutine Node_Component_Disk_Exponential_Radius_Solve_Set
   
  double precision function Node_Component_Disk_Exponential_Velocity(thisNode)
    !% Return the circular velocity of the exponential disk.
    implicit none
    type (treeNode         ), pointer, intent(inout) :: thisNode
    class(nodeComponentDisk), pointer                :: thisDiskComponent

    thisDiskComponent => thisNode%disk()
    Node_Component_Disk_Exponential_Velocity=thisDiskComponent%velocity()
    return
  end function Node_Component_Disk_Exponential_Velocity

  subroutine Node_Component_Disk_Exponential_Velocity_Set(thisNode,velocity)
    !% Set the circular velocity of the exponential disk.
    implicit none
    type (treeNode         ), pointer, intent(inout) :: thisNode
    double precision        ,          intent(in   ) :: velocity
    class(nodeComponentDisk), pointer                :: thisDiskComponent
    
    thisDiskComponent => thisNode%disk()
    call thisDiskComponent%velocitySet(velocity)
    return
  end subroutine Node_Component_Disk_Exponential_Velocity_Set
  
   !# <radiusSolverTask>
   !#  <unitName>Node_Component_Disk_Exponential_Radius_Solver</unitName>
   !# </radiusSolverTask>
   subroutine Node_Component_Disk_Exponential_Radius_Solver(thisNode,componentActive,specificAngularMomentum,Radius_Get,Radius_Set,Velocity_Get&
        &,Velocity_Set)
     !% Interface for the size solver algorithm.
     use Numerical_Constants_Physical
     implicit none
     type     (treeNode         ), pointer, intent(inout) :: thisNode
     logical,                               intent(  out) :: componentActive
     double precision,                      intent(  out) :: specificAngularMomentum
     procedure(double precision ), pointer, intent(  out) :: Radius_Get,Velocity_Get
     procedure(                 ), pointer, intent(  out) :: Radius_Set,Velocity_Set
     class    (nodeComponentDisk), pointer                :: thisDiskComponent
     double precision                                     :: specificAngularMomentumMean,angularMomentum,diskMass
     
     ! Determine if thisNode has an active disk component supported by this module.    
     componentActive=.false.
     thisDiskComponent => thisNode%disk()
     select type (thisDiskComponent)
     class is (nodeComponentDiskExponential)
        componentActive=.true.
        ! Get the angular momentum.
        angularMomentum=thisDiskComponent%angularMomentum()
        if (angularMomentum >= 0.0d0) then
           ! Compute the specific angular momentum at the scale radius, assuming a flat rotation curve.
           diskMass= thisDiskComponent%massGas    () &
                &   +thisDiskComponent%massStellar()
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
           if (diskRadiusSolverCole2000Method)                                                              &
                & specificAngularMomentum=dsqrt(                                                            &
                &                               max(                                                        &
                &                                    0.0d0,                                                 &
                &                                    specificAngularMomentum**2                             &
                &                                   -diskRadiusSolverFlatVsSphericalFactor                  &
                &                                   *gravitationalConstantGalacticus                        &
                &                                   *diskMass                                               &
                &                                   *Node_Component_Disk_Exponential_Radius_Solve(thisNode) &
                &                                  )                                                        &
                )

           ! Associate the pointers with the appropriate property routines.
           Radius_Get   => Node_Component_Disk_Exponential_Radius_Solve
           Radius_Set   => Node_Component_Disk_Exponential_Radius_Solve_Set
           Velocity_Get => Node_Component_Disk_Exponential_Velocity
           Velocity_Set => Node_Component_Disk_Exponential_Velocity_Set
        else
           call Node_Component_Disk_Exponential_Radius_Solve_Set(thisNode,0.0d0)
           call Node_Component_Disk_Exponential_Velocity_Set    (thisNode,0.0d0)
           componentActive=.false.
        end if
     end select
     return
   end subroutine Node_Component_Disk_Exponential_Radius_Solver
   
   double precision function Node_Component_Disk_Exponential_Star_Formation_Rate(self)
     !% Return the star formation rate of the exponential disk.
     use Star_Formation_Timescales_Disks
     implicit none
     class(nodeComponentDiskExponential), intent(inout) :: self
     type (treeNode                    ), pointer       :: thisNode
     double precision                                   :: starFormationTimescale,gasMass

     ! Get the associated node.
     thisNode => self%host()

     ! Get the star formation timescale.
     starFormationTimescale=Star_Formation_Timescale_Disk(thisNode)
      
     ! Get the gas mass.
     gasMass=self%massGas()
       
     ! If timescale is finite and gas mass is positive, then compute star formation rate.
     if (starFormationTimescale > 0.0d0 .and. gasMass > 0.0d0) then
        Node_Component_Disk_Exponential_Star_Formation_Rate=gasMass/starFormationTimescale
     else
        Node_Component_Disk_Exponential_Star_Formation_Rate=0.0d0
     end if
     return
   end function Node_Component_Disk_Exponential_Star_Formation_Rate

  !# <galacticusStateStoreTask>
  !#  <unitName>Node_Component_Disk_Exponential_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Node_Component_Disk_Exponential_State_Store(stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    use FGSL
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    write (stateFile) rotationCurveHalfRadiusMinimum,rotationCurveHalfRadiusMaximum,scaleLengthFactor,scaleLengthFactorSet,&
                    & rotationCurveGradientHalfRadiusMinimum,rotationCurveGradientHalfRadiusMaximum
    return
  end subroutine Node_Component_Disk_Exponential_State_Store
  
  !# <galacticusStateRetrieveTask>
  !#  <unitName>Node_Component_Disk_Exponential_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Node_Component_Disk_Exponential_State_Retrieve(stateFile,fgslStateFile)
    !% Retrieve the tabulation state from the file.
    use FGSL
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    ! Read the minimum and maximum tabulated times.
    read (stateFile) rotationCurveHalfRadiusMinimum,rotationCurveHalfRadiusMaximum,scaleLengthFactor,scaleLengthFactorSet,&
                   & rotationCurveGradientHalfRadiusMinimum,rotationCurveGradientHalfRadiusMaximum
    ! Flag that the table is now uninitialized.
    rotationCurveInitialized=.false.
    return
  end subroutine Node_Component_Disk_Exponential_State_Retrieve
  
   !# <mergerTreeExtraOutputTask>
   !#  <unitName>Node_Component_Disk_Exponential_Star_Formation_History_Output</unitName>
   !# </mergerTreeExtraOutputTask>
   subroutine Node_Component_Disk_Exponential_Star_Formation_History_Output(thisNode,iOutput,treeIndex,nodePassesFilter)
     !% Store the star formation history in the output file.
     use Kind_Numbers
     use Histories
     use Galacticus_Output_Star_Formation_Histories
     implicit none
     type   (treeNode         ), intent(inout), pointer :: thisNode
     integer                   , intent(in   )          :: iOutput
     integer(kind=kind_int8   ), intent(in   )          :: treeIndex
     logical                   , intent(in   )          :: nodePassesFilter
     class  (nodeComponentDisk),                pointer :: thisDiskComponent
     type   (history          )                         :: starFormationHistory

     ! Output the star formation history if a disk exists for this component.
     thisDiskComponent => thisNode%disk()
     select type (thisDiskComponent)
     class is (nodeComponentDiskExponential)
        starFormationHistory=thisDiskComponent%starFormationHistory()
        call Star_Formation_History_Output(thisNode,nodePassesFilter,starFormationHistory,iOutput,treeIndex,'disk')
        call thisDiskComponent%starFormationHistorySet(starFormationHistory)
     end select
     return
   end subroutine Node_Component_Disk_Exponential_Star_Formation_History_Output

end module Node_Component_Disk_Exponential
