!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module that implements a very simple disk component.

module Node_Component_Disk_Very_Simple
  !% Implements a very simple disk component.
  use ISO_Varying_String
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Disk_Very_Simple_Post_Evolve  , Node_Component_Disk_Very_Simple_Rate_Compute         , &
       &    Node_Component_Disk_Very_Simple_Scale_Set    , Node_Component_Disk_Very_Simple_Satellite_Merging    , &
       &    Node_Component_Disk_Very_Simple_Initialize   , Node_Component_Disk_Very_Simple_Pre_Evolve

  !# <component>
  !#  <class>disk</class>
  !#  <name>verySimple</name>
  !#  <isDefault>no</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>isInitialized</name>
  !#     <type>logical</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </property>
  !#   <property>
  !#     <name>massStellar</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of stars in the very simple disk."/>
  !#   </property>
  !#   <property>
  !#     <name>massGas</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" makeGeneric="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of gas in the very simple disk."/>
  !#   </property>
  !#   <property>
  !#     <name>stellarPropertiesHistory</name>
  !#     <type>history</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#   </property>
  !#  </properties>
  !#  <bindings>
  !#   <binding method="attachPipe"   function="Node_Component_Disk_Very_Simple_Attach_Pipe" description="Attach pipes to the very simple disk component." bindsTo="component" returnType="\void" arguments="" />
  !#   <binding method="enclosedMass" function="Node_Component_Disk_Very_Simple_Enclosed_Mass" bindsTo="component" />
  !#  </bindings>
  !#  <functions>objects.nodes.components.disk.very_simple.bound_functions.inc</functions>
  !# </component>

  ! Record of whether this module has been initialized.
  logical          :: moduleInitialized          =.false.

  ! Parameters controlling the physical implementation.
  double precision :: diskOutflowTimescaleMinimum        , diskStarFormationTimescaleMinimum, &
       &              diskVerySimpleMassScaleAbsolute

contains

  !# <nodeComponentInitializationTask>
  !#  <unitName>Node_Component_Disk_Very_Simple_Initialize</unitName>
  !# </nodeComponentInitializationTask>
  subroutine Node_Component_Disk_Very_Simple_Initialize()
    !% Initializes the tree node very simple disk component module.
    use Input_Parameters
    implicit none
    type(nodeComponentDiskVerySimple) :: diskVerySimpleComponent

    ! Initialize the module if necessary.
    !$omp critical (Node_Component_Disk_Very_Simple_Initialize)
    if (defaultDiskComponent%verySimpleIsActive().and..not.moduleInitialized) then
       ! Read parameters controlling the physical implementation.
       !@ <inputParameter>
       !@   <name>diskVerySimpleMassScaleAbsolute</name>
       !@   <defaultValue>$100 M_\odot$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The absolute mass scale below which calculations in the very simple disk component are allowed to become inaccurate.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('diskVerySimpleMassScaleAbsolute',diskVerySimpleMassScaleAbsolute,defaultValue=100.0d0)
       !@ <inputParameter>
       !@   <name>diskOutflowTimescaleMinimum</name>
       !@   <defaultValue>$10^{-3}$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The minimum timescale (in units of the halo dynamical time) on which outflows may deplete gas in the disk.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('diskOutflowTimescaleMinimum',diskOutflowTimescaleMinimum,defaultValue=1.0d-3)
       !@ <inputParameter>
       !@   <name>diskStarFormationTimescaleMinimum</name>
       !@   <defaultValue>$10^{-3}$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The minimum timescale (in units of the halo dynamical time) on which star formation may occur in the disk.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('diskStarFormationTimescaleMinimum',diskStarFormationTimescaleMinimum,defaultValue=1.0d-3)
       ! Attach the cooling mass pipe from the hot halo component.
       call diskVerySimpleComponent%attachPipe()
       ! Record that the module is now initialized.
       moduleInitialized=.true.
    end if
    !$omp end critical (Node_Component_Disk_Very_Simple_Initialize)
    return
  end subroutine Node_Component_Disk_Very_Simple_Initialize

  !# <preEvolveTask>
  !# <unitName>Node_Component_Disk_Very_Simple_Pre_Evolve</unitName>
  !# </preEvolveTask>
  subroutine Node_Component_Disk_Very_Simple_Pre_Evolve(thisNode)
    !% Ensure the disk has been initialized.
    implicit none
    type (treeNode         ), intent(inout), pointer :: thisNode
    class(nodeComponentDisk)               , pointer :: thisDiskComponent

    ! Get the disk component.
    thisDiskComponent => thisNode%disk()
    ! Check if an exponential disk component exists.
    select type (thisDiskComponent)
       class is (nodeComponentDiskVerySimple)
          ! Initialize the disk
       call Node_Component_Disk_Very_Simple_Create(thisNode)
    end select
    return
  end subroutine Node_Component_Disk_Very_Simple_Pre_Evolve

  !# <postEvolveTask>
  !# <unitName>Node_Component_Disk_Very_Simple_Post_Evolve</unitName>
  !# </postEvolveTask>
  subroutine Node_Component_Disk_Very_Simple_Post_Evolve(thisNode)
    !% Catch rounding errors in the very simple disk gas evolution.
    use Galacticus_Display
    use String_Handling
    use Histories
    implicit none
    type            (treeNode          ), intent(inout), pointer :: thisNode
    class           (nodeComponentDisk )               , pointer :: thisDiskComponent
    class           (nodeComponentBasic)               , pointer :: thisBasicComponent
    double precision                    , save                   :: fractionalErrorMaximum=0.0d0
    double precision                                             :: diskMass                    , fractionalError
    character       (len=20            )                         :: valueString
    type            (varying_string    )                         :: message
    type            (history           )                         :: stellarPropertiesHistory

    ! Get the disk component.
    thisDiskComponent => thisNode%disk()
    ! Check if a very simple disk component exists.
    select type (thisDiskComponent)
    class is (nodeComponentDiskVerySimple)
       ! Trim the stellar populations properties future history.
       thisBasicComponent => thisNode%basic()
       stellarPropertiesHistory=thisDiskComponent%stellarPropertiesHistory()
       call stellarPropertiesHistory%trim(thisBasicComponent%time())
       call thisDiskComponent%stellarPropertiesHistorySet(stellarPropertiesHistory)
       ! Trap negative gas masses.
       if (thisDiskComponent%massGas() < 0.0d0) then
          ! Check if this exceeds the maximum previously recorded error.
          fractionalError=   abs(thisDiskComponent%massGas    ()) &
               &          /(                                       &
               &                  thisDiskComponent%massStellar()  &
               &            +abs(thisDiskComponent%massGas    ()) &
               &           )
          !$omp critical (Very_Simple_Disk_Post_Evolve_Check)
          if (fractionalError > fractionalErrorMaximum) then
             ! Report a warning.
             message='Warning: disk has negative gas mass (fractional error exceeds any previously reported):'//char(10)
             message=message//'  Node index        = '//thisNode%index() //char(10)
             write (valueString,'(e12.6)') thisDiskComponent%massGas()
             message=message//'  Disk gas mass     = '//trim(valueString)//char(10)
             write (valueString,'(e12.6)') thisDiskComponent%massStellar()
             message=message//'  Disk stellar mass = '//trim(valueString)//char(10)
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
          !$omp end critical (Very_Simple_Disk_Post_Evolve_Check)
          ! Get the total mass of the disk material
          diskMass= thisDiskComponent%massGas    () &
               &   +thisDiskComponent%massStellar()
          if (diskMass == 0.0d0) call thisDiskComponent%massStellarSet(0.0d0)
          ! Reset the gas mass of the disk.
          call thisDiskComponent%massGasSet(0.0d0)
       end if
    end select
    return
  end subroutine Node_Component_Disk_Very_Simple_Post_Evolve

  subroutine Node_Component_Disk_Very_Simple_Create(thisNode)
    !% Create properties in a very simple disk component.
    use Histories
    use Stellar_Population_Properties
    implicit none
    type   (treeNode         ), intent(inout), pointer :: thisNode
    class  (nodeComponentDisk)               , pointer :: thisDiskComponent
    type   (history          )                         :: stellarPropertiesHistory
    logical                                            :: createStellarPropertiesHistory

    ! Get the disk component.
    thisDiskComponent => thisNode%disk()
    ! Exit if already initialized.
    if (thisDiskComponent%isInitialized()) return
    ! Determine which histories must be created.
    stellarPropertiesHistory      =thisDiskComponent%stellarPropertiesHistory        ()
    createStellarPropertiesHistory=.not.             stellarPropertiesHistory%exists ()
    call                                             stellarPropertiesHistory%destroy()
    ! Create the stellar properties history.
    if (createStellarPropertiesHistory) then
       ! Create the stellar properties history.
       call Stellar_Population_Properties_History_Create (thisNode,stellarPropertiesHistory)
       call thisDiskComponent%stellarPropertiesHistorySet(         stellarPropertiesHistory)
    end if
    ! Record that the disk has been initialized.
    call thisDiskComponent%isInitializedSet(.true.)
    return
  end subroutine Node_Component_Disk_Very_Simple_Create

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Disk_Very_Simple_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Disk_Very_Simple_Rate_Compute(node,interrupt,interruptProcedureReturn)
    !% Compute the very simple disk node mass rate of change.
    use Star_Formation_Feedback_Disks
    use Stellar_Feedback
    use Stellar_Population_Properties
    use Dark_Matter_Halo_Scales
    use Abundances_Structure
    use Galactic_Structure_Options
    use Histories
    use Stellar_Luminosities_Structure
    implicit none
    type            (treeNode                    ), intent(inout), pointer :: node
    class           (nodeComponentDisk           )               , pointer :: disk
    class           (nodeComponentHotHalo        )               , pointer :: hotHalo
    logical                                       , intent(inout)          :: interrupt
    procedure       (Interrupt_Procedure_Template), intent(inout), pointer :: interruptProcedureReturn
    procedure       (Interrupt_Procedure_Template)               , pointer :: interruptProcedure
    class           (darkMatterHaloScaleClass    )              , pointer :: darkMatterHaloScale_
    type            (abundances                  ), save                   :: fuelAbundancesRates     , stellarAbundancesRates
    !$omp threadprivate(stellarAbundancesRates,fuelAbundancesRates)
    double precision                                                       :: diskDynamicalTime       , diskMass              , &
         &                                                                    energyInputRate         , fuelMass              , &
         &                                                                    massOutflowRate         , starFormationRate     , &
         &                                                                    stellarMassRate         , fuelMassRate
    type            (stellarLuminosities         )                         :: luminositiesStellarRates
    type            (history                     )                         :: stellarHistoryRate
    
    ! Get a local copy of the interrupt procedure.
    interruptProcedure => interruptProcedureReturn
    ! Get the disk and check that it is of our class.
    disk => node%disk()
    select type (disk)
    class is (nodeComponentDiskVerySimple)
       ! Check for a realistic disk, return immediately if disk is unphysical.
       if (disk%massGas() < 0.0d0) return
       ! Interrupt if the disk is not initialized.
       if (.not.disk%isInitialized()) then
          interrupt=.true.
          interruptProcedureReturn => Node_Component_Disk_Very_Simple_Create
          return
       end if
       ! Compute the star formation rate.
       starFormationRate=Node_Component_Disk_Very_Simple_SFR(node)
       ! Get the available fuel mass.
       fuelMass         =disk%massGas()
       ! Find rates of change of stellar mass, and gas mass.
       stellarHistoryRate=disk%stellarPropertiesHistory()
       call Stellar_Population_Properties_Rates(starFormationRate,zeroAbundances,componentTypeDisk,node,stellarHistoryRate &
            &,stellarMassRate,stellarAbundancesRates,luminositiesStellarRates,fuelMassRate,fuelAbundancesRates,energyInputRate)
       if (stellarHistoryRate%exists()) call disk%stellarPropertiesHistoryRate(stellarHistoryRate)
       ! Adjust rates.
       call disk%massStellarRate(stellarMassRate)
       call disk%    massGasRate(   fuelMassRate)
       ! Find rate of outflow of material from the disk and pipe it to the outflowed reservoir.
       massOutflowRate=Star_Formation_Feedback_Disk_Outflow_Rate(node,starFormationRate,energyInputRate)
       if (massOutflowRate > 0.0d0) then
          ! Get the masses of the disk.
          diskMass=fuelMass+disk%massStellar()
          ! Limit the outflow rate timescale to a multiple of the dynamical time.
          darkMatterHaloScale_ => darkMatterHaloScale()
          diskDynamicalTime=darkMatterHaloScale_%dynamicalTimescale(node)
          massOutflowRate=min(massOutflowRate,fuelMass/diskOutflowTimescaleMinimum/diskDynamicalTime)
          ! Push to the hot halo.
          hotHalo => node%hotHalo()
          call hotHalo%outflowingMassRate(+massOutflowRate)
          call disk   %massGasRate       (-massOutflowRate)
       end if
    end select
    ! Return the procedure pointer.
    interruptProcedureReturn => interruptProcedure
    return
  end subroutine Node_Component_Disk_Very_Simple_Rate_Compute

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Disk_Very_Simple_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Disk_Very_Simple_Scale_Set(thisNode)
    !% Set scales for properties of {\tt thisNode}.
    use Abundances_Structure
    use Histories
    use Stellar_Population_Properties
    implicit none
    type            (treeNode         ), intent(inout), pointer :: thisNode
    class           (nodeComponentDisk)               , pointer :: thisDiskComponent
    double precision                                            :: mass
    type            (history          )                         :: stellarPopulationHistoryScales

    ! Get the disk component.
    thisDiskComponent => thisNode%disk()
    ! Check if a very simple disk component exists.
    select type (thisDiskComponent)
    class is (nodeComponentDiskVerySimple)
       ! Set scale for gas and stellar mass.
       mass=thisDiskComponent%massGas()+thisDiskComponent%massStellar()
       call thisDiskComponent%massGasScale    (max(mass,diskVerySimpleMassScaleAbsolute))
       call thisDiskComponent%massStellarScale(max(mass,diskVerySimpleMassScaleAbsolute))
       ! Set scales for stellar population properties and star formation histories.
       stellarPopulationHistoryScales=thisDiskComponent%stellarPropertiesHistory()
       call Stellar_Population_Properties_Scales           (stellarPopulationHistoryScales,thisDiskComponent%massStellar(),zeroAbundances)
       call thisDiskComponent%stellarPropertiesHistoryScale(stellarPopulationHistoryScales                                               )
       call stellarPopulationHistoryScales%destroy()
    end select
    return
  end subroutine Node_Component_Disk_Very_Simple_Scale_Set

  !# <satelliteMergerTask>
  !#  <unitName>Node_Component_Disk_Very_Simple_Satellite_Merging</unitName>
  !#  <after>Satellite_Merging_Mass_Movement_Store</after>
  !#  <after>Satellite_Merging_Remnant_Size</after>
  !# </satelliteMergerTask>
  subroutine Node_Component_Disk_Very_Simple_Satellite_Merging(thisNode)
    !% Transfer any very simple disk associated with {\tt thisNode} to its host halo.
    use Satellite_Merging_Mass_Movements_Descriptors
    use Galacticus_Error
    implicit none
    type (treeNode         ), intent(inout), pointer :: thisNode
    type (treeNode         )               , pointer :: hostNode
    class(nodeComponentDisk)               , pointer :: hostDiskComponent, thisDiskComponent

    ! Check that the disk is of the verySimple class.
    thisDiskComponent => thisNode%disk()
    select type (thisDiskComponent)
    class is (nodeComponentDiskVerySimple)

       ! Find the node to merge with and its disk component.
       hostNode          => thisNode%mergesWith()
       hostDiskComponent => hostNode%disk      ()

       ! Move the gas component of the very simple disk to the host.
       select case (thisMergerGasMovesTo)
       case (movesToDisk)
          call hostDiskComponent%massGasSet    (                                            &
               &                                 hostDiskComponent%    massGas()            &
               &                                +thisDiskComponent%    massGas()            &
               &                               )
       case (movesToSpheroid)
          call Galacticus_Error_Report(                                                     &
               &                       'Node_Component_Disk_Very_Simple_Satellite_Merging', &
               &                       'this component does not work with spheroids'        &
               &                      )
       case default
          call Galacticus_Error_Report(                                                     &
               &                       'Node_Component_Disk_Very_Simple_Satellite_Merging', &
               &                       'unrecognized movesTo descriptor'                    &
               &                      )
       end select
       call    thisDiskComponent%massGasSet    (                                            &
            &                                                              0.0d0            &
            &                                  )

       ! Move the stellar component of the very simple disk to the host.
       select case (thisMergerStarsMoveTo)
       case (movesToDisk)
          call hostDiskComponent%massStellarSet(                                            &
               &                                 hostDiskComponent%massStellar()            &
               &                                +thisDiskComponent%massStellar()            &
               &                               )
       case (movesToSpheroid)
          call Galacticus_Error_Report(                                                     &
               &                       'Node_Component_Disk_Very_Simple_Satellite_Merging', &
               &                       'this component does not work with spheroids'        &
               &                      )
       case default
          call Galacticus_Error_Report(                                                     &
               &                       'Node_Component_Disk_Very_Simple_Satellite_Merging', &
               &                       'unrecognized movesTo descriptor'                    &
               &                      )
       end select
       call    thisDiskComponent%massStellarSet(                                            &
            &                                                              0.0d0            &
            &                                  )
    end select
    return
  end subroutine Node_Component_Disk_Very_Simple_Satellite_Merging

  double precision function Node_Component_Disk_Very_Simple_SFR(thisNode)
    !% Return the star formation rate of the very simple disk.
    use Star_Formation_Timescales_Disks
    use Dark_Matter_Halo_Scales
    implicit none
    type            (treeNode          ), intent(inout), pointer :: thisNode
    class           (nodeComponentDisk )               , pointer :: thisDiskComponent
    class           (darkMatterHaloScaleClass)               , pointer :: darkMatterHaloScale_
    double precision                                             :: diskDynamicalTime, gasMass, starFormationTimescale

    ! Get the disk component.
    thisDiskComponent => thisNode%disk()

    ! Get the star formation timescale.
    starFormationTimescale=Star_Formation_Timescale_Disk(thisNode)

    ! Limit the star formation timescale to a multiple of the dynamical time.
    darkMatterHaloScale_   => darkMatterHaloScale()
    diskDynamicalTime      =darkMatterHaloScale_%dynamicalTimescale(thisNode)
    starFormationTimescale =max(starFormationTimescale,diskStarFormationTimescaleMinimum*diskDynamicalTime)

    ! Get the gas mass.
    gasMass=thisDiskComponent%massGas()

    ! If timescale is finite and gas mass is positive, then compute star formation rate.
    if (starFormationTimescale > 0.0d0 .and. gasMass > 0.0d0) then
       Node_Component_Disk_Very_Simple_SFR=gasMass/starFormationTimescale
    else
       Node_Component_Disk_Very_Simple_SFR=0.0d0
    end if
    return
  end function Node_Component_Disk_Very_Simple_SFR

end module Node_Component_Disk_Very_Simple
