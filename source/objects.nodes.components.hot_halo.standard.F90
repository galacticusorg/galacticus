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

!% Contains a module which implements the standard hot halo node component.

module Node_Component_Hot_Halo_Standard
  !% Implements the standard hot halo node component.
  use Galacticus_Nodes
  use Radiation_Structure
  implicit none
  private
  public :: Node_Component_Hot_Halo_Standard_Initialize  , Node_Component_Hot_Halo_Standard_Thread_Initialize, &
       &    Node_Component_Hot_Halo_Standard_Post_Evolve , Node_Component_Hot_Halo_Standard_Reset            , &
       &    Node_Component_Hot_Halo_Standard_Scale_Set   , Node_Component_Hot_Halo_Standard_Tree_Initialize  , &
       &    Node_Component_Hot_Halo_Standard_Node_Merger , Node_Component_Hot_Halo_Standard_Satellite_Merger , &
       &    Node_Component_Hot_Halo_Standard_Promote     , Node_Component_Hot_Halo_Standard_Formation        , &
       &    Node_Component_Hot_Halo_Standard_Rate_Compute, Node_Component_Hot_Halo_Standard_Pre_Evolve

  !# <component>
  !#  <class>hotHalo</class>
  !#  <name>standard</name>
  !#  <isDefault>yes</isDefault>
  !#  <createFunction isDeferred="true" />
  !#  <properties>
  !#   <property>
  !#     <name>isInitialized</name>
  !#     <type>logical</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </property>
  !#   <property>
  !#     <name>mass</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of gas in the hot halo."/>
  !#   </property>
  !#   <property>
  !#     <name>abundances</name>
  !#     <type>abundances</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of metals in the hot phase of the hot halo."/>
  !#   </property>
  !#   <property>
  !#     <name>chemicals</name>
  !#     <type>chemicals</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#   </property>
  !#   <property>
  !#     <name>angularMomentum</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#     <output unitsInSI="massSolar*megaParsec*kilo" comment="Angular momentum of gas in the hot halo."/>
  !#   </property>
  !#   <property>
  !#     <name>outflowedMass</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of outflowed gas in the hot halo."/>
  !#   </property>
  !#   <property>
  !#     <name>outflowedAngularMomentum</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <output unitsInSI="massSolar*megaParsec*kilo" comment="Angular momentum of outflowed gas in the hot halo."/>
  !#   </property>
  !#   <property>
  !#     <name>outflowedAbundances</name>
  !#     <type>abundances</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of metals in the outflowed phase of the hot halo."/>
  !#   </property>
  !#   <property>
  !#     <name>outflowingMass</name>
  !#     <attributes isSettable="false" isGettable="false" isEvolvable="true" isDeferred="rate" />
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <isVirtual>true</isVirtual>
  !#   </property>
  !#   <property>
  !#     <name>outflowingAngularMomentum</name>
  !#     <attributes isSettable="false" isGettable="false" isEvolvable="true" isDeferred="rate" />
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <isVirtual>true</isVirtual>
  !#   </property>
  !#   <property>
  !#     <name>outflowingAbundances</name>
  !#     <attributes isSettable="false" isGettable="false" isEvolvable="true" isDeferred="rate" />
  !#     <type>abundances</type>
  !#     <rank>0</rank>
  !#     <isVirtual>true</isVirtual>
  !#   </property>
  !#   <property>
  !#     <name>unaccretedMass</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of gas that failed to accrete into the hot halo."/>
  !#   </property>
  !#   <property>
  !#     <name>outerRadius</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" isDeferred="get" />
  !#     <output unitsInSI="megaParsec" comment="Outer radius of the hot halo."/>
  !#   </property>
  !#   <property>
  !#     <name>strippedMass</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#   </property>
  !#   <property>
  !#     <name>strippedAbundances</name>
  !#     <type>abundances</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#   </property>
  !#   <property>
  !#     <name>hotHaloCoolingMass</name>
  !#     <attributes isSettable="false" isGettable="false" isEvolvable="true" isDeferred="rate" bindsTo="top" />
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <isVirtual>true</isVirtual>
  !#   </property>
  !#   <property>
  !#     <name>hotHaloCoolingAngularMomentum</name>
  !#     <attributes isSettable="false" isGettable="false" isEvolvable="true" isDeferred="rate" bindsTo="top" />
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <isVirtual>true</isVirtual>
  !#   </property>
  !#   <property>
  !#     <name>hotHaloCoolingAbundances</name>
  !#     <attributes isSettable="false" isGettable="false" isEvolvable="true" isDeferred="rate" bindsTo="top" />
  !#     <type>abundances</type>
  !#     <rank>0</rank>
  !#     <isVirtual>true</isVirtual>
  !#   </property>
  !#   <property>
  !#     <name>massSink</name>
  !#     <attributes isSettable="false" isGettable="false" isEvolvable="true" isDeferred="rate" />
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <isVirtual>true</isVirtual>
  !#   </property>
  !#   <property>
  !#     <name>heatSource</name>
  !#     <attributes isSettable="false" isGettable="false" isEvolvable="true" isDeferred="rate" />
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <isVirtual>true</isVirtual>
  !#   </property>
  !#  </properties>
  !#  <bindings>
  !#    <binding method="outflowReturn" bindsTo="component" isDeferred="true" >
  !#     <interface>
  !#      <type>void</type>
  !#      <self pass="true" intent="inout" />
  !#      <argument>logical, intent(inout) :: interrupt</argument>
  !#      <argument>procedure(Interrupt_Procedure_Template), intent(inout), pointer :: interruptProcedure</argument>
  !#     </interface>
  !#    </binding>
  !#    <binding method="outerRadiusGrowthRate" bindsTo="component" isDeferred="true" >
  !#     <interface>
  !#      <type>real</type>
  !#      <rank>0</rank>
  !#      <self pass="true" intent="inout" />
  !#     </interface>
  !#    </binding>
  !#  </bindings>
  !# </component>

  ! Internal count of abundances and chemicals.
  integer                                                   :: abundancesCount                                 , chemicalsCount

  ! Configuration variables.
  logical                                                   :: hotHaloExcessHeatDrivesOutflow                  , hotHaloNodeMergerLimitBaryonFraction        , &
       &                                                       hotHaloAngularMomentumAlwaysGrows               , hotHaloOutflowReturnOnFormation
  double precision                                          :: hotHaloExpulsionRateMaximum                     , hotHaloOutflowStrippingEfficiency

  ! Quantities stored to avoid repeated computation.
  logical                                                   :: gotAngularMomentumCoolingRate           =.false., gotCoolingRate                      =.false.
  double precision                                          :: angularMomentumHeatingRateRemaining             , coolingRate                                 , &
       &                                                       massHeatingRateRemaining
  !$omp threadprivate(gotCoolingRate,gotAngularMomentumCoolingRate,coolingRate,massHeatingRateRemaining,angularMomentumHeatingRateRemaining)
  ! Radiation structure.
  type            (radiationStructure          )            :: radiation
  !$omp threadprivate(radiation)
  ! Record of whether this module has been initialized.
  logical                                                   :: moduleInitialized                       =.false.

  ! Tracked properties control.
  logical                                                   :: hotHaloTrackStrippedGas

contains

  !# <mergerTreePreTreeConstructionTask>
  !#  <unitName>Node_Component_Hot_Halo_Standard_Initialize</unitName>
  !# </mergerTreePreTreeConstructionTask>
  subroutine Node_Component_Hot_Halo_Standard_Initialize()
    !% Initializes the standard hot halo component module.
    use ISO_Varying_String
    use Input_Parameters
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Galacticus_Error
    use Node_Component_Hot_Halo_Standard_Data
    implicit none
    type(varying_string              ) :: hotHaloCoolingFromText
    type(nodeComponentHotHaloStandard) :: hotHaloComponent

    ! Initialize the module if necessary.
    !$omp critical (Node_Component_Hot_Halo_Standard_Initialize)
    if (defaultHotHaloComponent%standardIsActive().and..not.moduleInitialized) then

       ! Get numbers of abundance and chemicals properties.
       abundancesCount=Abundances_Property_Count()
       chemicalsCount =Chemicals_Property_Count ()

       ! Determine whether satellite nodes will be starved of gas.
       !@ <inputParameter>
       !@   <name>starveSatellites</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Specifies whether or not the hot halo should be removed (``starved'') when a node becomes a satellite.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('starveSatellites',starveSatellites,defaultValue=.false.)

       ! Determine whether stripped material should be tracked.
       !@ <inputParameter>
       !@   <name>hotHaloTrackStrippedGas</name>
       !@   <defaultValue>true</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Specifies whether or not gas stripped from the hot halo should be tracked.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('hotHaloTrackStrippedGas',hotHaloTrackStrippedGas,defaultValue=.true.)

       ! Determine whether outflowed gas should be restored to the hot reservoir on halo formation events.
       !@ <inputParameter>
       !@   <name>hotHaloOutflowReturnOnFormation</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Specifies whether or not outflowed gas should be returned to the hot reservoir on halo formation events.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('hotHaloOutflowReturnOnFormation',hotHaloOutflowReturnOnFormation,defaultValue=.false.)

       ! Determine whether negative angular momentum accretion rates onto the halo should be treated as positive for the purposes
       ! of computing the hot halo angular momentum.
       !@ <inputParameter>
       !@   <name>hotHaloAngularMomentumAlwaysGrows</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether or not negative rates of accretion of angular momentum into the hot halo will be treated as positive
       !@     for the purposes of computing the hot halo angular momentum.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('hotHaloAngularMomentumAlwaysGrows',hotHaloAngularMomentumAlwaysGrows,defaultValue=.false.)

       ! Determine whether the angular momentum of cooling gas should be computed from the "current node" or the "formation node".
       !@ <inputParameter>
       !@   <name>hotHaloCoolingFromNode</name>
       !@   <defaultValue>currentNode</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Specifies whether the angular momentum of cooling gas should be computed from the ``current node'' or the ``formation node''.
       !@   </description>
       !@   <type>integer</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('hotHaloCoolingFromNode',hotHaloCoolingFromText,defaultValue='currentNode')
       select case (char(hotHaloCoolingFromText))
       case ("currentNode"  )
          hotHaloCoolingFromNode=currentNode
       case ("formationNode")
          hotHaloCoolingFromNode=formationNode
       case default
          call Galacticus_Error_Report('Node_Component_Hot_Halo_Standard_Initialize','hotHaloCoolingFromNode must be one of "currentNode" or "formationNode"')
       end select

       ! Determine whether excess heating of the halo will drive an outflow.
       !@ <inputParameter>
       !@   <name>hotHaloExcessHeatDrivesOutflow</name>
       !@   <defaultValue>true</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Specifies whether heating of the halo in excess of its cooling rate will drive an outflow from the halo.
       !@   </description>
       !@   <type>integer</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('hotHaloExcessHeatDrivesOutflow',hotHaloExcessHeatDrivesOutflow,defaultValue=.true.)

       ! Get rate (in units of halo inverse dynamical time) at which outflowed gas returns to the hot gas reservoir.
       !@ <inputParameter>
       !@   <name>hotHaloOutflowReturnRate</name>
       !@   <defaultValue>5</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Specifies the rate at which reheated mass is returned to the hot phase in units of the inverse halo dynamical time.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('hotHaloOutflowReturnRate',hotHaloOutflowReturnRate,defaultValue=5.0d0)

       ! Get efficiency with which outflowing gas is stripped from the hot halo.
       !@ <inputParameter>
       !@   <name>hotHaloOutflowStrippingEfficiency</name>
       !@   <defaultValue>0.1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Specifies the efficiency with which outflowing gas is stripped from the hot halo, following the prescription of \citeauthor{font_colours_2008}~(\citeyear{font_colours_2008}; i.e. this is the parameter $\epsilon_{\rm strip}$ in their eqn.~6).
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('hotHaloOutflowStrippingEfficiency',hotHaloOutflowStrippingEfficiency,defaultValue=0.1d0)

       ! Get the maximum rate (in units of halo inverse dynamical time) at which gas can be expelled from the halo.
       !@ <inputParameter>
       !@   <name>hotHaloExpulsionRateMaximum</name>
       !@   <defaultValue>1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Specifies the maximum rate at which mass can be expelled from the hot halo in units of the inverse halo dynamical time.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('hotHaloExpulsionRateMaximum',hotHaloExpulsionRateMaximum,defaultValue=1.0d0)

       ! Get fraction of angular momentum that is lost during cooling/infall.
       !@ <inputParameter>
       !@   <name>hotHaloAngularMomentumLossFraction</name>
       !@   <defaultValue>0.3</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Specifies the fraction of angular momentum that is lost from cooling/infalling gas.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('hotHaloAngularMomentumLossFraction',hotHaloAngularMomentumLossFraction,defaultValue=0.3d0)

       ! Get option controlling limiting of baryon fraction during node mergers.
       !@ <inputParameter>
       !@   <name>hotHaloNodeMergerLimitBaryonFraction</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Controls whether the hot gas content of nodes should be limited to not exceed the universal baryon fraction at node
       !@    merger events. If set to {\tt true}, hot gas (and angular momentum, abundances, and chemicals proportionally) will be
       !@    removed from the merged halo to the unaccreted gas reservoir to limit the baryonic mass to the universal baryon
       !@    fraction where possible.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('hotHaloNodeMergerLimitBaryonFraction',hotHaloNodeMergerLimitBaryonFraction,defaultValue=.false.)

       ! Bind the outer radius get function.
       call hotHaloComponent%                  outerRadiusFunction(Node_Component_Hot_Halo_Standard_Outer_Radius              )

       ! Bind the heat source pipe to the function that will handle heat input to the hot halo.
       call hotHaloComponent%               heatSourceRateFunction(Node_Component_Hot_Halo_Standard_Heat_Source               )

       ! Bind outflowing material pipes to the functions that will handle input of outflowing material to the hot halo.
       call hotHaloComponent%           outflowingMassRateFunction(Node_Component_Hot_Halo_Standard_Outflowing_Mass_Rate      )
       call hotHaloComponent%outflowingAngularMomentumRateFunction(Node_Component_Hot_Halo_Standard_Outflowing_Ang_Mom_Rate   )
       call hotHaloComponent%     outflowingAbundancesRateFunction(Node_Component_Hot_Halo_Standard_Outflowing_Abundances_Rate)

       ! Bind a creation function.
       call hotHaloComponent%                    createFunctionSet(Node_Component_Hot_Halo_Standard_Initializor               )

       ! Bind the mass sink function.
       call hothaloComponent%                 massSinkRateFunction(Node_Component_Hot_Halo_Standard_Mass_Sink                 )

       ! Bind the outflow return function.
       call hotHaloComponent%                outflowReturnFunction(Node_Component_Hot_Halo_Standard_Outflow_Return            )

       ! Bind the outer radius growth rate function.
       call hotHaloComponent%        outerRadiusGrowthRateFunction(Node_Component_Hot_Halo_Standard_Outer_Radius_Growth_Rate  )

       ! Record that the module is now initialized.
       moduleInitialized=.true.

    end if
    !$omp end critical (Node_Component_Hot_Halo_Standard_Initialize)
    return
  end subroutine Node_Component_Hot_Halo_Standard_Initialize

  !# <mergerTreeEvolveThreadInitialize>
  !#  <unitName>Node_Component_Hot_Halo_Standard_Thread_Initialize</unitName>
  !# </mergerTreeEvolveThreadInitialize>
  subroutine Node_Component_Hot_Halo_Standard_Thread_Initialize
    !% Initializes the tree node hot halo methods module.
    implicit none

    ! Check if this implementation is selected. Define the radiation component to include both the CMB and the intergalactic background if it is.
    if (defaultHotHaloComponent%standardIsActive()) call radiation%define([radiationTypeCMB,radiationTypeIGB])
    return
  end subroutine Node_Component_Hot_Halo_Standard_Thread_Initialize

  !# <calculationResetTask>
  !# <unitName>Node_Component_Hot_Halo_Standard_Reset</unitName>
  !# </calculationResetTask>
  subroutine Node_Component_Hot_Halo_Standard_Reset(thisNode)
    !% Remove memory of stored computed values as we're about to begin computing derivatives anew.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    gotCoolingRate               =.false.
    gotAngularMomentumCoolingRate=.false.
    return
  end subroutine Node_Component_Hot_Halo_Standard_Reset

  double precision function Node_Component_Hot_Halo_Standard_Outer_Radius(self)
    !% Return the outer radius in the standard hot halo.
    use Dark_Matter_Halo_Scales
    implicit none
    class(nodeComponentHotHaloStandard), intent(inout) :: self
    type (treeNode                    ), pointer       :: selfHost

    selfHost => self%host()
    Node_Component_Hot_Halo_Standard_Outer_Radius=max(min(self%outerRadiusValue(),Dark_Matter_Halo_Virial_Radius(selfHost)),0.0d0)
    return
  end function Node_Component_Hot_Halo_Standard_Outer_Radius

  !# <postEvolveTask>
  !#  <unitName>Node_Component_Hot_Halo_Standard_Post_Evolve</unitName>
  !# </postEvolveTask>
  subroutine Node_Component_Hot_Halo_Standard_Post_Evolve(thisNode)
    !% Do processing of the node required after evolution.
    use Abundances_Structure
    use Dark_Matter_Halo_Scales
    use Node_Component_Hot_Halo_Standard_Data
    implicit none
    type (treeNode            ), intent(inout), pointer :: thisNode
    type (treeNode            )               , pointer :: parentNode
    class(nodeComponentHotHalo)               , pointer :: parentHotHaloComponent, thisHotHaloComponent
    class(nodeComponentSpin   )               , pointer :: parentSpinComponent
    
    ! Limit hot gas mass to be non-negative.
    thisHotHaloComponent => thisNode%hotHalo()
    select type (thisHotHaloComponent)
    class is (nodeComponentHotHaloStandard)
       if (thisHotHaloComponent%mass() < 0.0d0) call thisHotHaloComponent%massSet(0.0d0)
    end select
    ! Process hot gas for satellites.
    if (thisNode%isSatellite()) then
       if (starveSatellites) then
          thisHotHaloComponent => thisNode%hotHalo()
          select type (thisHotHaloComponent)
          class is (nodeComponentHotHaloStandard)
             ! Transfer any outflowed gas to the hot halo of the parent node.
             parentNode => thisNode%parent
             do while (parentNode%isSatellite())
                parentNode => parentNode%parent
             end do
             parentHotHaloComponent => parentNode%hotHalo()
             parentSpinComponent    => parentNode%spin   ()
             call parentHotHaloComponent%outflowedAngularMomentumSet(                                                   &
                  &                                                   parentHotHaloComponent%outflowedAngularMomentum() &
                  &                                                  +  thisHotHaloComponent%outflowedMass           () &
                  &                                                  *   parentSpinComponent%spin                    () &
                  &                                                  *Dark_Matter_Halo_Virial_Radius  (parentNode)      &
                  &                                                  *Dark_Matter_Halo_Virial_Velocity(parentNode)      &
                  &                                                 )
             call   thisHotHaloComponent%outflowedAngularMomentumSet(                                                   &
                  &                                                   0.0d0                                             &
                  &                                                 )
             call parentHotHaloComponent%outflowedMassSet           (                                                   &
                  &                                                   parentHotHaloComponent%outflowedMass           () &
                  &                                                  +  thisHotHaloComponent%outflowedMass           () &
                  &                                                 )
             call   thisHotHaloComponent%outflowedMassSet           (                                                   &
                  &                                                   0.0d0                                             &
                  &                                                 )
             call parentHotHaloComponent%outflowedAbundancesSet     (                                                   &
                  &                                                   parentHotHaloComponent%outflowedAbundances     () &
                  &                                                  +  thisHotHaloComponent%outflowedAbundances     () &
                  &                                                 )
             call thisHotHaloComponent%outflowedAbundancesSet       (                                                   &
                  &                                                   zeroAbundances                                    &
                  &                                                 )
          end select
       end if
       ! Check if stripped mass is being tracked.
       if (hotHaloTrackStrippedGas) then
          thisHotHaloComponent => thisNode%hotHalo()
          select type (thisHotHaloComponent)
          class is (nodeComponentHotHaloStandard)
             ! Transfer any stripped gas to the host halo.
             parentNode => thisNode%parent
             do while (parentNode%isSatellite())
                parentNode => parentNode%parent
             end do
             call Node_Component_Hot_Halo_Standard_Create(parentNode)
             parentHotHaloComponent => parentNode%hotHalo()
             parentSpinComponent    => parentNode%spin   ()
             call parentHotHaloComponent%outflowedAngularMomentumSet(                                                   &
                  &                                                   parentHotHaloComponent%outflowedAngularMomentum() &
                  &                                                  +  thisHotHaloComponent%strippedMass            () &
                  &                                                  *   parentSpinComponent%spin                    () &
                  &                                                  *Dark_Matter_Halo_Virial_Radius  (parentNode)      &
                  &                                                  *Dark_Matter_Halo_Virial_Velocity(parentNode)      &
                  &                                                 )
             call parentHotHaloComponent%outflowedMassSet           (                                                   &
                  &                                                   parentHotHaloComponent%outflowedMass           () &
                  &                                                  +  thisHotHaloComponent%strippedMass            () &
                  &                                                 )
             call   thisHotHaloComponent%strippedMassSet            (                                                   &
                  &                                                   0.0d0                                             &
                  &                                                 )
             call parentHotHaloComponent%outflowedAbundancesSet     (                                                   &
                  &                                                   parentHotHaloComponent%outflowedAbundances     () &
                  &                                                  +  thisHotHaloComponent%strippedAbundances      () &
                  &                                                 )
             call thisHotHaloComponent%strippedAbundancesSet        (                                                   &
                  &                                                   zeroAbundances                                    &
                  &                                                 )
          end select
       end if
    end if
    return
  end subroutine Node_Component_Hot_Halo_Standard_Post_Evolve

  subroutine Node_Component_Hot_Halo_Standard_Strip_Gas_Rate(thisNode,gasMassRate,interrupt,interruptProcedure)
    !% Add gas stripped from the hot halo to the stripped gas reservoirs under the assumption of uniformly distributed properties
    !% (e.g. fully-mixed metals).
    implicit none
    type            (treeNode                    ), intent(inout), pointer :: thisNode
    double precision                              , intent(in   )          :: gasMassRate
    logical                                       , intent(inout)          :: interrupt          
    procedure       (Interrupt_Procedure_Template), intent(inout), pointer :: interruptProcedure 
    class           (nodeComponentHotHalo        )               , pointer :: thisHotHaloComponent
    double precision                                                       :: gasMass

    ! Exit immediately for zero rate.
    if (gasMassRate == 0.0d0) return

    ! Get the hot halo component.
    thisHotHaloComponent => thisNode%hotHalo()
    select type (thisHotHaloComponent)
    class is (nodeComponentHotHaloStandard)
       ! Get the gas mass present.
       gasMass=thisHotHaloComponent%mass()
       ! If gas is present, adjust the rates.
       if (gasMass > 0.0d0) then
          ! Mass.
          call thisHotHaloComponent%      strippedMassRate(                                  gasMassRate        ,interrupt,interruptProcedure)
          ! Metal abundances.
          call thisHotHaloComponent%strippedAbundancesRate(thisHotHaloComponent%abundances()*gasMassRate/gasMass,interrupt,interruptProcedure)
       end if
    end select
    return
  end subroutine Node_Component_Hot_Halo_Standard_Strip_Gas_Rate

  subroutine Node_Component_Hot_Halo_Standard_Heat_Source(thisHotHaloComponent,rate,interrupt,interruptProcedure)
    !% An incoming pipe for sources of heating to the hot halo.
    use Galacticus_Error
    use Dark_Matter_Halo_Scales
    implicit none
    class           (nodeComponentHotHalo        ), intent(inout)                    :: thisHotHaloComponent
    double precision                              , intent(in   )                    :: rate
    logical                                       , intent(inout), optional          :: interrupt
    procedure       (Interrupt_Procedure_Template), intent(inout), optional, pointer :: interruptProcedure
    type            (treeNode                    )                         , pointer :: thisNode
    double precision                                                                 :: excessMassHeatingRate, inputMassHeatingRate, massHeatingRate

     ! Trap cases where an attempt is made to remove energy via this input function.
     if (rate < 0.0d0) call Galacticus_Error_Report('Node_Component_Hot_Halo_Standard_Heat_Source','attempt to remove energy via heat source pipe to hot halo')

     ! Get the node associated with this hot halo component.
     thisNode => thisHotHaloComponent%host()

     ! Ensure that the cooling rate has been computed.
     call Node_Component_Hot_Halo_Standard_Cooling_Rate(thisNode)

     ! Compute the input mass heating rate from the input energy heating rate.
     inputMassHeatingRate=rate/Dark_Matter_Halo_Virial_Velocity(thisNode)**2

     ! Limit the mass heating rate such that it never exceeds the remaining budget.
     massHeatingRate=min(inputMassHeatingRate,massHeatingRateRemaining)

     ! Update the remaining budget of allowed mass heating rate.
     if (massHeatingRateRemaining-massHeatingRate <= 0.0d0) massHeatingRate=massHeatingRateRemaining
     massHeatingRateRemaining=max(massHeatingRateRemaining-massHeatingRate,0.0d0)

     ! Call routine to apply this mass heating rate to all hot halo cooling pipes.
     call Node_Component_Hot_Halo_Standard_Push_To_Cooling_Pipes(thisNode,-massHeatingRate,interrupt,interruptProcedure)

     ! If requested, compute the rate at which an outflow is driven from the halo by excess heating.
     if (hotHaloExcessHeatDrivesOutflow) then

        ! Compute the excess mass heating rate (i.e. that beyond which is being used to offset the cooling rate).
        excessMassHeatingRate=inputMassHeatingRate-massHeatingRate

        ! Remove any excess mass heating rate from the halo.
        call Node_Component_Hot_Halo_Standard_Push_From_Halo(thisNode,excessMassHeatingRate)

     end if

     return
   end subroutine Node_Component_Hot_Halo_Standard_Heat_Source

  subroutine Node_Component_Hot_Halo_Standard_Push_To_Cooling_Pipes(thisNode,massRate,interrupt,interruptProcedure)
    !% Push mass through the cooling pipes (along with appropriate amounts of metals and angular momentum) at the given rate.
    use Cooling_Infall_Radii
    use Cooling_Specific_Angular_Momenta
    use Abundances_Structure
    use Node_Component_Hot_Halo_Standard_Data
    implicit none
    type            (treeNode                    ), intent(inout)          , pointer :: thisNode
    double precision                              , intent(in   )                    :: massRate
    logical                                       , intent(inout), optional          :: interrupt
    procedure       (Interrupt_Procedure_Template), intent(inout), optional, pointer :: interruptProcedure
    type            (treeNode                    )                         , pointer :: coolingFromNode
    class           (nodeComponentHotHalo        )                         , pointer :: coolingFromHotHaloComponent, thisHotHaloComponent
    type            (abundances                  ), save                             :: abundancesCoolingRate
    !$omp threadprivate(abundancesCoolingRate)
    double precision                                                                 :: angularMomentumCoolingRate , infallRadius

    ! Get the hot halo component.
    thisHotHaloComponent => thisNode%hotHalo()
    select type (thisHotHaloComponent)
    class is (nodeComponentHotHaloStandard)

       ! Ignore zero rates.
       if (massRate /= 0.0d0 .and. thisHotHaloComponent%mass() > 0.0d0 .and. thisHotHaloComponent%angularMomentum() > 0.0d0) then

          ! Remove mass from the hot component.
          call    thisHotHaloComponent%massRate       (-massRate                             )
          ! Pipe the mass rate to whichever component claimed it.
          if (thisHotHaloComponent%hotHaloCoolingMassRateIsAttached()) then
             call thisHotHaloComponent%hotHaloCoolingMassRate(+massRate,interrupt,interruptProcedure)
             if (interrupt) return
          end if

          ! Find the node to use for cooling calculations.
          select case (hotHaloCoolingFromNode)
          case (currentNode  )
             coolingFromNode => thisNode
          case (formationNode)
             coolingFromNode => thisNode%formationNode
          end select
          infallRadius=Infall_Radius(thisNode)
          angularMomentumCoolingRate=massRate*Cooling_Specific_Angular_Momentum(coolingFromNode,infallRadius)
          if (.not.gotAngularMomentumCoolingRate) then
             angularMomentumHeatingRateRemaining=coolingRate*Cooling_Specific_Angular_Momentum(coolingFromNode,infallRadius)
             gotAngularMomentumCoolingRate=.true.
          end if

          if (massRate < 0.0d0) then
             if (massHeatingRateRemaining == 0.0d0) then
                angularMomentumCoolingRate=-angularMomentumHeatingRateRemaining
                angularMomentumHeatingRateRemaining=0.0d0
             else
                angularMomentumCoolingRate=max(angularMomentumCoolingRate,-angularMomentumHeatingRateRemaining)
                angularMomentumHeatingRateRemaining=angularMomentumHeatingRateRemaining+angularMomentumCoolingRate
             end if
          end if
          call    thisHotHaloComponent%angularMomentumRate       (     -angularMomentumCoolingRate                                                                                  )
          ! Pipe the cooling rate to which ever component claimed it.
          if (thisHotHaloComponent%hotHaloCoolingAngularMomentumRateIsAttached()) then
             call thisHotHaloComponent%hotHaloCoolingAngularMomentumRate(sign(+angularMomentumCoolingRate*(1.0d0-hotHaloAngularMomentumLossFraction),massRate),interrupt,interruptProcedure)
             if (interrupt) return
          end if
          ! Get the rate of change of abundances.
          coolingFromHotHaloComponent => coolingFromNode%hotHalo()
          abundancesCoolingRate=coolingFromHotHaloComponent%abundances()
          abundancesCoolingRate=massRate*abundancesCoolingRate/coolingFromHotHaloComponent%mass()
          call    thisHotHaloComponent%abundancesRate       (-abundancesCoolingRate                             )
          ! Pipe the cooling rate to which ever component claimed it.
          if (thisHotHaloComponent%hotHaloCoolingAbundancesRateIsAttached()) then
             call thisHotHaloComponent%hotHaloCoolingAbundancesRate(+abundancesCoolingRate,interrupt,interruptProcedure)
             if (interrupt) return
          end if
       end if
    end select
    return
  end subroutine Node_Component_Hot_Halo_Standard_Push_To_Cooling_Pipes

  subroutine Node_Component_Hot_Halo_Standard_Push_From_Halo(thisNode,massRate)
    !% Push mass from the hot halo into an infinite sink (along with appropriate amounts of metals, chemicals and angular momentum) at the given rate.
    use Dark_Matter_Halo_Scales
    use Abundances_Structure
    use Chemical_Abundances_Structure
    implicit none
    type            (treeNode            )      , intent(inout), pointer :: thisNode
    double precision                            , intent(in   )          :: massRate
    class           (nodeComponentHotHalo)                     , pointer :: thisHotHalo
    type            (abundances          ), save                         :: abundancesRates
    type            (chemicalAbundances  ), save                         :: chemicalsRates
    !$omp threadprivate(abundancesRates,chemicalsRates)
    double precision                                                     :: angularMomentumRate , massRateLimited

    ! Get the hot halo component.
    thisHotHalo => thisNode%hotHalo()
    ! Ignore zero rates.
    if (massRate /= 0.0d0 .and. thisHotHalo%mass() > 0.0d0) then
       ! Limit the mass expulsion rate to a fraction of the halo dynamical timescale.
       massRateLimited=min(massRate,hotHaloExpulsionRateMaximum*thisHotHalo%mass()/Dark_Matter_Halo_Dynamical_Timescale(thisNode))
       ! Get the rate of change of abundances, chemicals, and angular momentum.
       abundancesRates    =thisHotHalo%abundances     ()*massRateLimited/thisHotHalo%mass()
       angularMomentumRate=thisHotHalo%angularMomentum()*massRateLimited/thisHotHalo%mass()
       chemicalsRates     =thisHotHalo%chemicals      ()*massRateLimited/thisHotHalo%mass()
       ! Remove mass, etc. from the hot component.
       call thisHotHalo%           massRate(-    massRateLimited)
       call thisHotHalo%     abundancesRate(-    abundancesRates)
       call thisHotHalo%angularMomentumRate(-angularMomentumRate)
       call thisHotHalo%      chemicalsRate(-     chemicalsRates)
       ! If this node is a satellite and stripped gas is being tracked, move mass and abundances to the stripped reservoir.
       if (thisNode%isSatellite().and.hotHaloTrackStrippedGas) then
          call thisHotHalo%      strippedMassRate(+massRateLimited)
          call thisHotHalo%strippedAbundancesRate(+abundancesRates)
       end if
    end if
    return
  end subroutine Node_Component_Hot_Halo_Standard_Push_From_Halo

  double precision function Node_Component_Hot_Halo_Standard_Outflow_Stripped_Fraction(thisNode,thisHotHaloComponent)
    !% Compute the fraction of material outflowing into the hot halo of {\tt thisNode} which is susceptible to being stripped
    !% away.
    use Hot_Halo_Mass_Distributions
    use Dark_Matter_Halo_Scales
    implicit none
    type            (treeNode                    ), intent(inout), pointer :: thisNode
    class           (nodeComponentHotHaloStandard)                         :: thisHotHaloComponent
    class           (hotHaloMassDistributionClass)               , pointer :: defaultHotHaloMassDistribution
    double precision                                                       :: massOuter                     , massVirial  , &
         &                                                                    radiusOuter                   , radiusVirial

    defaultHotHaloMassDistribution => hotHaloMassDistribution()
    radiusOuter =thisHotHaloComponent%outerRadius()
    radiusVirial=Dark_Matter_Halo_Virial_Radius             (thisNode             )
    massOuter   =defaultHotHaloMassDistribution%enclosedMass(thisNode,radiusOuter )
    massVirial  =defaultHotHaloMassDistribution%enclosedMass(thisNode,radiusVirial)
    if (massVirial > 0.0d0) then
       Node_Component_Hot_Halo_Standard_Outflow_Stripped_Fraction=hotHaloOutflowStrippingEfficiency*(1.0d0-massOuter/massVirial)
    else
       Node_Component_Hot_Halo_Standard_Outflow_Stripped_Fraction=hotHaloOutflowStrippingEfficiency
    end if
    return
  end function Node_Component_Hot_Halo_Standard_Outflow_Stripped_Fraction

  subroutine Node_Component_Hot_Halo_Standard_Outflowing_Mass_Rate(self,rate,interrupt,interruptProcedure)
    !% Accept outflowing gas from a galaxy and deposit it into the outflowed and stripped reservoirs.
    implicit none
    class           (nodeComponentHotHalo                                    ), intent(inout)                    :: self
    double precision                                                          , intent(in   )                    :: rate
    logical                                                                   , intent(inout), optional          :: interrupt
    procedure       (Interrupt_Procedure_Template                            ), intent(inout), optional, pointer :: interruptProcedure
    type            (treeNode                                                )                         , pointer :: selfNode
    double precision                                                                                             :: strippedOutflowFraction

    select type (self)
    class is (nodeComponentHotHaloStandard)
       ! Get the host node.
       selfNode => self%host()
       if (selfNode%isSatellite().and.hotHaloTrackStrippedGas) then
          strippedOutflowFraction=Node_Component_Hot_Halo_Standard_Outflow_Stripped_Fraction(selfNode,self)
          call self% strippedMassRate(rate*       strippedOutflowFraction )
       else
          strippedOutflowFraction=0.0d0
       end if
       ! Funnel the outflow gas into the outflowed and stripped reservoirs in the computed proportions.
       call    self%outflowedMassRate(rate*(1.0d0-strippedOutflowFraction))
    end select
    return
  end subroutine Node_Component_Hot_Halo_Standard_Outflowing_Mass_Rate

  subroutine Node_Component_Hot_Halo_Standard_Outflowing_Ang_Mom_Rate(self,rate,interrupt,interruptProcedure)
    !% Accept outflowing gas angular momentum from a galaxy and deposit it into the outflowed reservoir.
    use Node_Component_Hot_Halo_Standard_Data
    implicit none
    class           (nodeComponentHotHalo                                    ), intent(inout)                    :: self
    double precision                                                          , intent(in   )                    :: rate
    logical                                                                   , intent(inout), optional          :: interrupt
    procedure       (Interrupt_Procedure_Template                            ), intent(inout), optional, pointer :: interruptProcedure
    type            (treeNode                                                )                         , pointer :: selfNode
    double precision                                                                                             :: strippedOutflowFraction

    select type (self)
    class is (nodeComponentHotHaloStandard)
       ! Get the host node.
       selfNode => self%host()
       if (selfNode%isSatellite().and.hotHaloTrackStrippedGas) then
          strippedOutflowFraction=Node_Component_Hot_Halo_Standard_Outflow_Stripped_Fraction(selfNode,self)
       else
          strippedOutflowFraction=0.0d0
       end if
       call    self%outflowedAngularMomentumRate(rate*(1.0d0-strippedOutflowFraction)/(1.0d0-hotHaloAngularMomentumLossFraction))
    end select
    return
  end subroutine Node_Component_Hot_Halo_Standard_Outflowing_Ang_Mom_Rate

  subroutine Node_Component_Hot_Halo_Standard_Outflowing_Abundances_Rate(self,rate,interrupt,interruptProcedure)
    !% Accept outflowing gas abundances from a galaxy and deposit it into the outflowed reservoir.
    use Abundances_Structure
    implicit none
    class           (nodeComponentHotHalo                                    ), intent(inout)                    :: self
    type            (abundances                                              ), intent(in   )                    :: rate
    logical                                                                   , intent(inout), optional          :: interrupt
    procedure       (Interrupt_Procedure_Template                            ), intent(inout), optional, pointer :: interruptProcedure
    type            (treeNode                                                )                         , pointer :: selfNode
    double precision                                                                                             :: strippedOutflowFraction

    select type (self)
    class is (nodeComponentHotHaloStandard)
       ! Get the host node.
       selfNode => self%host()
       if (selfNode%isSatellite().and.hotHaloTrackStrippedGas) then
          strippedOutflowFraction=Node_Component_Hot_Halo_Standard_Outflow_Stripped_Fraction(selfNode,self)
          call    self%strippedAbundancesRate(rate*       strippedOutflowFraction )
       else
          strippedOutflowFraction=0.0d0
       end if
       call    self%outflowedAbundancesRate(rate*(1.0d0-strippedOutflowFraction))
    end select
    return
  end subroutine Node_Component_Hot_Halo_Standard_Outflowing_Abundances_Rate

  !# <preEvolveTask>
  !# <unitName>Node_Component_Hot_Halo_Standard_Pre_Evolve</unitName>
  !# </preEvolveTask>
  subroutine Node_Component_Hot_Halo_Standard_Pre_Evolve(thisNode)
    !% Ensure the standard hot halo has been initialized.
    implicit none
    type (treeNode            ), intent(inout), pointer :: thisNode
    class(nodeComponentHotHalo)               , pointer :: thisHotHaloComponent

    ! Get the spheroid component.
    thisHotHaloComponent => thisNode%hotHalo()
    ! Check if a standard hot halo component exists.
    select type (thisHotHaloComponent)
    class is (nodeComponentHotHaloStandard)
       ! Initialize the hot halo.
       call Node_Component_Hot_Halo_Standard_Initializor(thisHotHaloComponent)
    end select
    return
  end subroutine Node_Component_Hot_Halo_Standard_Pre_Evolve

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Hot_Halo_Standard_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Hot_Halo_Standard_Rate_Compute(thisNode,interrupt,interruptProcedure)
    !% Compute the hot halo node mass rate of change.
    use Abundances_Structure
    use Accretion_Halos
    use Accretion_Halos_Options
    use Dark_Matter_Halo_Spins
    use Dark_Matter_Halo_Scales
    use Chemical_States
    use Chemical_Abundances_Structure
    use Chemical_Reaction_Rates
    use Chemical_Reaction_Rates_Utilities
    use Numerical_Constants_Astronomical
    use Hot_Halo_Mass_Distributions
    use Hot_Halo_Ram_Pressure_Stripping_Timescales
    use Cosmology_Parameters
    use Node_Component_Hot_Halo_Standard_Data
    implicit none
    type            (treeNode                    )           , intent(inout), pointer :: thisNode
    logical                                                  , intent(inout)          :: interrupt
    procedure       (Interrupt_Procedure_Template)           , intent(inout), pointer :: interruptProcedure
    class           (nodeComponentHotHalo        )                          , pointer :: thisHotHaloComponent
    class           (nodeComponentBasic          )                          , pointer :: thisBasicComponent
    type            (abundances                  ), save                              :: accretionRateAbundances
    class           (cosmologyParametersClass    )                          , pointer :: thisCosmologyParameters
    class           (hotHaloMassDistributionClass)                          , pointer :: defaultHotHaloMassDistribution
    type            (chemicalAbundances          ), save                              :: accretionRateChemicals                   , chemicalDensities      , &
         &                                                                               chemicalDensitiesRates                   , chemicalMasses         , &
         &                                                                               chemicalMassesRates                      , chemicalsCoolingRate
    !$omp threadprivate(accretionRateAbundances,accretionRateChemicals,chemicalMasses)
    !$omp threadprivate(chemicalDensities,chemicalDensitiesRates,chemicalMassesRates,chemicalsCoolingRate)
    double precision                                                                  :: angularMomentumAccretionRate             , temperature            , &
         &                                                                               densityAtOuterRadius                     , failedMassAccretionRate, &
         &                                                                               massLossRate                             , massToDensityConversion, &
         &                                                                               outerRadius                              , outerRadiusGrowthRate  , &
         &                                                                               massAccretionRate                        , densityMinimum         , &
         &                                                                               radiusVirial

    ! Get the hot halo component.
    thisHotHaloComponent => thisNode%hotHalo()
    ! Ensure that the standard hot halo implementation is active.
    if (defaultHotHaloComponent%standardIsActive()) then
       ! Find the rate of gas mass accretion onto the halo.
       massAccretionRate      =Halo_Baryonic_Accretion_Rate       (thisNode,accretionModeHot)
       failedMassAccretionRate=Halo_Baryonic_Failed_Accretion_Rate(thisNode,accretionModeHot)
       ! Get the basic component.
       thisBasicComponent => thisNode%basic()
       ! Apply accretion rates.
       if (      massAccretionRate > 0.0d0 .or. thisHotHaloComponent%mass() > 0.0d0) &
            & call thisHotHaloComponent%          massRate(      massAccretionRate,interrupt,interruptProcedure)
       if (failedMassAccretionRate > 0.0d0 .or. thisHotHaloComponent%mass() > 0.0d0) &
            & call thisHotHaloComponent%unaccretedMassRate(failedMassAccretionRate,interrupt,interruptProcedure)
       ! Next compute the cooling rate in this halo.
       call Node_Component_Hot_Halo_Standard_Cooling_Rate(thisNode)
       ! Pipe the cooling rate to which ever component claimed it.
       call Node_Component_Hot_Halo_Standard_Push_To_Cooling_Pipes(thisNode,coolingRate,interrupt,interruptProcedure)
       ! Get the rate at which abundances are accreted onto this halo.
       call Halo_Baryonic_Accretion_Rate_Abundances(thisNode,accretionRateAbundances,accretionModeHot)
       call thisHotHaloComponent%abundancesRate(accretionRateAbundances,interrupt,interruptProcedure)
       ! Next block of tasks occur only if the accretion rate is non-zero.
          if (thisBasicComponent%accretionRate() /= 0.0d0) then
          ! Compute the rate of accretion of angular momentum.
          angularMomentumAccretionRate=Dark_Matter_Halo_Angular_Momentum_Growth_Rate(thisNode)*(massAccretionRate &
               &/thisBasicComponent%accretionRate())
             if (hotHaloAngularMomentumAlwaysGrows) angularMomentumAccretionRate=abs(angularMomentumAccretionRate)
          call thisHotHaloComponent%angularMomentumRate(angularMomentumAccretionRate,interrupt,interruptProcedure)
       end if
       ! Next block of tasks occur only if chemicals are being tracked.
       if (chemicalsCount > 0) then
          ! Get the rate at which chemicals are accreted onto this halo.
          call Halo_Baryonic_Accretion_Rate_Chemicals(thisNode,accretionRateChemicals,accretionModeHot)
          call thisHotHaloComponent%chemicalsRate(accretionRateChemicals,interrupt,interruptProcedure)
          ! For non-zero cooling rate, adjust the rates of chemical masses.
          if (coolingRate > 0.0d0) then
             ! Compute the rate at which chemicals are lost via cooling.
             chemicalsCoolingRate=thisHotHaloComponent%chemicals()*coolingRate/thisHotHaloComponent%mass()
             ! Adjust the rates of chemical masses accordingly.
             call thisHotHaloComponent%chemicalsRate(-chemicalsCoolingRate,interrupt,interruptProcedure)
          end if
          ! Compute the rates of change of chemical masses due to chemical reactions.
          ! Get the temperature of the hot reservoir.
          temperature=Dark_Matter_Halo_Virial_Temperature(thisNode)
          ! Set the radiation background.
          call radiation%set(thisNode)
          ! Get the masses of chemicals.
          chemicalMasses=thisHotHaloComponent%chemicals()
          ! Truncate masses to zero to avoid unphysical behavior.
          call chemicalMasses%enforcePositive()
          ! Scale all chemical masses by their mass in atomic mass units to get a number density.
          call chemicalMasses%massToNumber(chemicalDensities)
          ! Compute factor converting mass of chemicals in (M_Solar/M_Atomic) to number density in cm^-3.
          massToDensityConversion=Chemicals_Mass_To_Density_Conversion(Dark_Matter_Halo_Virial_Radius(thisNode))
          ! Convert to number density.
          chemicalDensities=chemicalDensities*massToDensityConversion
          ! Compute the chemical reaction rates.
          call Chemical_Reaction_Rate(chemicalDensitiesRates,temperature,chemicalDensities,radiation)
          ! Convert to mass change rates.
          call chemicalDensitiesRates%numberToMass(chemicalMassesRates)
          chemicalMassesRates=chemicalMassesRates*gigaYear/massToDensityConversion
          ! Adjust rates appropriately.
          call thisHotHaloComponent%chemicalsRate(chemicalMassesRates,interrupt,interruptProcedure)
       end if
       ! Perform return of outflowed material.
       select type (thisHotHaloComponent)
       class is (nodeComponentHotHaloStandard)
          call thisHotHaloComponent%outflowReturn(interrupt,interruptProcedure)
          ! Test whether this halo is a satellite or not.
          if (thisNode%isSatellite()) then
             ! For satellites, get the current ram pressure stripping radius for this hot halo.
             outerRadiusGrowthRate=thisHotHaloComponent%outerRadiusGrowthRate()
             outerRadius          =thisHotHaloComponent%outerRadius          ()
             if     (                                                                                                           &
                  &   outerRadiusGrowthRate                  /= 0.0d0                                                           &
                  &  .and.                                                                                                      &
                  &   thisHotHaloComponent%mass           () >  0.0d0                                                           &
                  &  .and.                                                                                                      &
                  &   outerRadius                 <=                                   Dark_Matter_Halo_Virial_Radius(thisNode) &
                  &  .and.                                                                                                      &
                  &   outerRadius                 > outerRadiusOverVirialRadiusMinimum*Dark_Matter_Halo_Virial_Radius(thisNode) &
                  & ) then
                defaultHotHaloMassDistribution => hotHaloMassDistribution()                
                densityAtOuterRadius =defaultHotHaloMassDistribution%density(thisNode,outerRadius)
                massLossRate=4.0d0*Pi*densityAtOuterRadius*outerRadius**2*outerRadiusGrowthRate
                call thisHotHaloComponent%outerRadiusRate(+outerRadiusGrowthRate,interrupt,interruptProcedure)
                call thisHotHaloComponent%   massSinkRate(+         massLossRate,interrupt,interruptProcedure)
                call Node_Component_Hot_Halo_Standard_Strip_Gas_Rate(thisNode,-massLossRate,interrupt,interruptProcedure)
             end if
          else
             ! For isolated halos, the outer radius should grow with the virial radius.
             call thisHotHaloComponent%outerRadiusRate(Dark_Matter_Halo_Virial_Radius_Growth_Rate(thisNode),interrupt,interruptProcedure)
          end if
       end select
    end if
    return
  end subroutine Node_Component_Hot_Halo_Standard_Rate_Compute
  
  double precision function Node_Component_Hot_Halo_Standard_Outer_Radius_Growth_Rate(self)
    !% Compute the growth rate of the outer radius of the hot halo.
    use Hot_Halo_Ram_Pressure_Stripping
    use Hot_Halo_Ram_Pressure_Stripping_Timescales
    implicit none
    class           (nodeComponentHotHaloStandard), intent(inout) :: self
    type            (treeNode                    ), pointer       :: selfNode
    double precision                                              :: ramPressureRadius, outerRadius

    selfNode          => self%hostNode
    ramPressureRadius =  Hot_Halo_Ram_Pressure_Stripping_Radius(selfNode)
    outerRadius       =  self%outerRadius()
    ! Test whether the ram pressure radius is smaller than the current outer radius of the hot gas profile.
    if     (                                           &
         &  ramPressureRadius      < outerRadius .and. &
         &  self%angularMomentum() >       0.0d0       &
         & ) then
       ! The ram pressure stripping radius is within the outer radius. Cause the outer radius to shrink to the ram pressure
       ! stripping radius on the halo dynamical timescale.
       Node_Component_Hot_Halo_Standard_Outer_Radius_Growth_Rate=  &
            &  (ramPressureRadius-outerRadius)                     &
            & /Hot_Halo_Ram_Pressure_Stripping_Timescale(selfNode)
    else
       Node_Component_Hot_Halo_Standard_Outer_Radius_Growth_Rate=0.0d0
    end if
    return
  end function Node_Component_Hot_Halo_Standard_Outer_Radius_Growth_Rate

  subroutine Node_Component_Hot_Halo_Standard_Outflow_Return(self,interrupt,interruptProcedure)
    !% Return outflowed gas to the hot halo.
    use Galacticus_Error
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Chemical_Reaction_Rates_Utilities
    use Chemical_States
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Atomic
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    use Node_Component_Hot_Halo_Standard_Data
    use Cosmology_Parameters
    use Hot_Halo_Mass_Distributions
    implicit none
    class           (nodeComponentHotHaloStandard), intent(inout)          :: self
    logical                                       , intent(inout)          :: interrupt
    procedure       (Interrupt_Procedure_Template), intent(inout), pointer :: interruptProcedure
    type            (treeNode                    ), pointer                :: selfNode
    class(nodeComponentBasic), pointer :: selfBasic
    class           (cosmologyParametersClass)               , pointer :: thisCosmologyParameters
    class           (hotHaloMassDistributionClass)               , pointer :: defaultHotHaloMassDistribution
    double precision                                                       :: outflowedMass            , massReturnRate         , &
         &                                                                    angularMomentumReturnRate, massToDensityConversion, &
         &                                                                    hydrogenByMass           , temperature            , &
         &                                                                    numberDensityHydrogen    , densityAtOuterRadius, radiusVirial, outerRadius, densityMinimum
    type            (abundances                  ), save                   :: abundancesReturnRate     , outflowedAbundances
    !$omp threadprivate(abundancesReturnRate,outflowedAbundances)
    type            (chemicalAbundances          ), save                   :: chemicalDensities        , chemicalMasses           , &
         &                                                                    chemicalMassesRates
    !$omp threadprivate(chemicalDensities,chemicalMassesRates,chemicalMasses)

    ! Get the hosting node.
    selfNode => self%hostNode
    ! Next tasks occur only for systems in which outflowed gas is being recycled.
    if (.not.starveSatellites.or..not.selfNode%isSatellite()) then
       outflowedMass            =self%outflowedMass()
       massReturnRate           =hotHaloOutflowReturnRate*outflowedMass                  /Dark_Matter_Halo_Dynamical_Timescale(selfNode)
       angularMomentumReturnRate=hotHaloOutflowReturnRate*self%outflowedAngularMomentum()/Dark_Matter_Halo_Dynamical_Timescale(selfNode)
       abundancesReturnRate     =hotHaloOutflowReturnRate*self%outflowedAbundances     ()/Dark_Matter_Halo_Dynamical_Timescale(selfNode)
       call self%           outflowedMassRate(-           massReturnRate,interrupt,interruptProcedure)
       call self%                    massRate(+           massReturnRate,interrupt,interruptProcedure)
       call self%outflowedAngularMomentumRate(-angularMomentumReturnRate,interrupt,interruptProcedure)
       call self%         angularMomentumRate(+angularMomentumReturnRate,interrupt,interruptProcedure)
       call self%     outflowedAbundancesRate(-     abundancesReturnRate,interrupt,interruptProcedure)
       call self%              abundancesRate(+     abundancesReturnRate,interrupt,interruptProcedure)
       ! The outer radius must be increased as the halo fills up with gas.
       outerRadius =self%outerRadius()
       radiusVirial=Dark_Matter_Halo_Virial_Radius(selfNode)
       if (outerRadius < radiusVirial) then 
          defaultHotHaloMassDistribution => hotHaloMassDistribution()
          densityAtOuterRadius=defaultHotHaloMassDistribution%density(selfNode,outerRadius)
          ! If the outer radius and density are non-zero we can expand the outer radius at a rate determined by the current
          ! density profile.
          if (outerRadius > 0.0d0 .and. densityAtOuterRadius > 0.0d0) then
             ! Limit the density at the outer radius to one third of the mean virial density (for baryons, assuming a
             ! universal baryon fraction) to prevent arbitrarily rapid growth of the outer radius in halos containing almost
             ! no gas.
             thisCosmologyParameters => cosmologyParameters()
             selfBasic => selfNode%basic()
             densityMinimum=(thisCosmologyParameters%omegaBaryon()/thisCosmologyParameters%omegaMatter())*selfBasic%mass()/radiusVirial**3/4.0d0/Pi
             call self%outerRadiusRate(                           &
                  &                     massReturnRate            &
                  &                    /4.0d0                     &
                  &                    /Pi                        &
                  &                    /outerRadius**2            &
                  &                    /max(                      &
                  &                         densityAtOuterRadius, &
                  &                         densityMinimum        &
                  &                        )                      &
                  &                   )
          ! Otherwise, if we have a positive rate of mass return, simply grow the radius at the virial velocity.
          else if (massReturnRate > 0.0d0) then
             ! Force some growth here so the radius is not trapped at zero.
             call self%outerRadiusRate(Dark_Matter_Halo_Virial_Velocity(selfNode)*kilo*gigaYear/megaParsec)
          end if
       end if
       ! If we have a non-zero return rate, compute associated chemical rates.
       if (chemicalsCount > 0 .and. massReturnRate /= 0.0d0) then
          ! Compute coefficient in conversion of mass to density for this node.
          massToDensityConversion=Chemicals_Mass_To_Density_Conversion(Dark_Matter_Halo_Virial_Radius(selfNode))/3.0d0
          ! Get the abundances of the outflowed material.
          outflowedAbundances    =self%outflowedAbundances()/outflowedMass
          ! Get the hydrogen mass fraction in outflowed gas.
          hydrogenByMass         =outflowedAbundances%hydrogenMassFraction()
          ! Compute the temperature and density of material in the hot halo.
          temperature            =Dark_Matter_Halo_Virial_Temperature(selfNode)
          numberDensityHydrogen  =hydrogenByMass*outflowedMass*massToDensityConversion/atomicMassHydrogen
          ! Set the radiation field.
          call radiation%set(selfNode)
          ! Get the chemical densities.
          call Chemical_Densities(chemicalDensities,temperature,numberDensityHydrogen,outflowedAbundances,radiation)
          ! Convert from densities to masses.
          call chemicalDensities%numberToMass(chemicalMasses)
          chemicalMassesRates=chemicalMasses*massReturnRate*hydrogenByMass/numberDensityHydrogen/atomicMassHydrogen
          ! Compute the rate at which chemicals are returned to the hot reservoir.
          call self%chemicalsRate(chemicalMassesRates,interrupt,interruptProcedure)
       end if
    end if
    return   
  end subroutine Node_Component_Hot_Halo_Standard_Outflow_Return

  subroutine Node_Component_Hot_Halo_Standard_Mass_Sink(self,setValue,interrupt,interruptProcedure)
    !% Account for a sink of gaseous material in the standard hot halo hot gas.
    use Galacticus_Error
    implicit none
    class           (nodeComponentHotHalo        ), intent(inout) :: self
    double precision                              , intent(in   ) :: setValue
    logical                                       , intent(inout), optional          :: interrupt          
    procedure       (Interrupt_Procedure_Template), intent(inout), optional, pointer :: interruptProcedure 

    select type (self)
    class is (nodeComponentHotHaloStandard)
       ! Trap cases where an attempt is made to add gas via this sink function.
       if (setValue > 0.0d0) call Galacticus_Error_Report('Node_Component_Hot_Halo_Standard_Mass_Sink','attempt to add mass via sink in hot halo')
       ! Proportionally adjust the rates of all components of the hot gas reservoir.
       call Node_Component_Hot_Halo_Standard_Hot_Gas_All_Rate(self,setValue,interrupt,interruptProcedure)
    end select
    return
  end subroutine Node_Component_Hot_Halo_Standard_Mass_Sink

  subroutine Node_Component_Hot_Halo_Standard_Hot_Gas_All_Rate(self,gasMassRate,interrupt,interruptProcedure)
    !% Adjusts the rates of all components of the hot gas reservoir under the assumption of uniformly distributed properties
    !% (e.g. fully-mixed metals).
    implicit none
    class           (nodeComponentHotHaloStandard), intent(inout)                    :: self
    double precision                              , intent(in   )                    :: gasMassRate
    logical                                       , intent(inout), optional          :: interrupt          
    procedure       (Interrupt_Procedure_Template), intent(inout), optional, pointer :: interruptProcedure 
    double precision                                                                 :: gasMass
    
    ! Exit immediately for zero rate.
    if (gasMassRate == 0.0d0) return
    ! Get the gas mass present.
    gasMass=self%mass()
    ! If gas is present, adjust the rates.
    if (gasMass > 0.0d0) then
       ! Mass.
       call self%           massRate(                       gasMassRate        ,interrupt,interruptProcedure)
       ! Angular momentum.
       call self%angularMomentumRate(self%angularMomentum()*gasMassRate/gasMass,interrupt,interruptProcedure)
       ! Metal abundances.
       call self%     abundancesRate(self%abundances     ()*gasMassRate/gasMass,interrupt,interruptProcedure)
       ! Chemical abundances.
       call self%      chemicalsRate(self%chemicals      ()*gasMassRate/gasMass,interrupt,interruptProcedure)
    end if
    return
  end subroutine Node_Component_Hot_Halo_Standard_Hot_Gas_All_Rate
  
  !# <scaleSetTask>
  !#  <unitName>Node_Component_Hot_Halo_Standard_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Hot_Halo_Standard_Scale_Set(thisNode)
    !% Set scales for properties of {\tt thisNode}.
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Dark_Matter_Halo_Scales
    implicit none
    type            (treeNode            ), intent(inout), pointer :: thisNode
    class           (nodeComponentHotHalo)               , pointer :: thisHotHaloComponent
    class           (nodeComponentBasic  )               , pointer :: thisBasicComponent
    double precision                      , parameter              :: scaleMassRelative   =1.0d-3
    double precision                      , parameter              :: scaleRadiusRelative =1.0d-1
    double precision                                               :: massVirial                 , radiusVirial, &
         &                                                            velocityVirial

    ! Get the hot halo component.
    thisHotHaloComponent => thisNode%hotHalo()
    ! Ensure that it is of the standard class.
    select type (thisHotHaloComponent)
    class is (nodeComponentHotHaloStandard)
       ! The the basic component.
       thisBasicComponent => thisNode%basic()
       ! Get virial properties.
       massVirial    =thisBasicComponent%mass()
       radiusVirial  =Dark_Matter_Halo_Virial_Radius  (thisNode)
       velocityVirial=Dark_Matter_Halo_Virial_Velocity(thisNode)
       call    thisHotHaloComponent%                    massScale(               massVirial                            *scaleMassRelative  )
       call    thisHotHaloComponent%           outflowedMassScale(               massVirial                            *scaleMassRelative  )
       call    thisHotHaloComponent%          unaccretedMassScale(               massVirial                            *scaleMassRelative  )
       call    thisHotHaloComponent%              abundancesScale(unitAbundances*massVirial                            *scaleMassRelative  )
       call    thisHotHaloComponent%     outflowedAbundancesScale(unitAbundances*massVirial                            *scaleMassRelative  )
       call    thisHotHaloComponent%               chemicalsScale(unitChemicals *massVirial                            *scaleMassRelative  )
       call    thisHotHaloComponent%         angularMomentumScale(               massVirial*radiusVirial*velocityVirial*scaleMassRelative  )
       call    thisHotHaloComponent%outflowedAngularMomentumScale(               massVirial*radiusVirial*velocityVirial*scaleMassRelative  )
       call    thisHotHaloComponent%             outerRadiusScale(                          radiusVirial               *scaleRadiusRelative)
       if (hotHaloTrackStrippedGas) then
          call thisHotHaloComponent%            strippedMassScale(               massVirial                            *scaleMassRelative  )
          call thisHotHaloComponent%      strippedAbundancesScale(unitAbundances*massVirial                            *scaleMassRelative  )
       end if
    end select
    return
  end subroutine Node_Component_Hot_Halo_Standard_Scale_Set

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Hot_Halo_Standard_Tree_Initialize</unitName>
  !#  <after>spin</after>
  !#  <after>darkMatterProfile</after>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Hot_Halo_Standard_Tree_Initialize(thisNode)
    !% Initialize the contents of the hot halo component for any sub-resolution accretion (i.e. the gas that would have been
    !% accreted if the merger tree had infinite resolution).
    use Accretion_Halos
    use Accretion_Halos_Options
    use Dark_Matter_Halo_Spins
    use Chemical_Abundances_Structure
    use Abundances_Structure
    implicit none
    type            (treeNode            ), intent(inout), pointer :: thisNode
    class           (nodeComponentHotHalo)               , pointer :: currentHotHaloComponent, thisHotHaloComponent
    class           (nodeComponentBasic  )               , pointer :: thisBasicComponent
    type            (abundances          ), save                   :: accretedAbundances
    type            (chemicalAbundances  ), save                   :: accretedChemicals
    !$omp threadprivate(accretedAbundances,accretedChemicals)
    double precision                                               :: angularMomentum        , failedHotHaloMass   , &
         &                                                            hotHaloMass

    ! If the node has a child or the standard hot halo is not active, then return immediately.
    if (associated(thisNode%firstChild).or..not.defaultHotHaloComponent%standardIsActive()) return

    ! Get the hot halo component.
    currentHotHaloComponent => thisNode%hotHalo()
    ! Ensure that it is of unspecified class.
    select type (currentHotHaloComponent)
    type is (nodeComponentHotHalo)
       ! Get the mass of hot gas accreted and the mass that failed to accrete.
       hotHaloMass      =Halo_Baryonic_Accreted_Mass       (thisNode,accretionModeHot)
       failedHotHaloMass=Halo_Baryonic_Failed_Accreted_Mass(thisNode,accretionModeHot)
       ! If either is non-zero, then create a hot halo component and add these masses to it.
       if (hotHaloMass > 0.0d0 .or. failedHotHaloMass > 0.0d0) then
          call Node_Component_Hot_Halo_Standard_Create(thisNode)
          thisHotHaloComponent => thisNode%hotHalo()
          thisBasicComponent   => thisNode%basic  ()
          call thisHotHaloComponent%           massSet(      hotHaloMass )
          call thisHotHaloComponent% unaccretedMassSet(failedHotHaloMass )
          ! Also add the appropriate angular momentum.
          angularMomentum=hotHaloMass*Dark_Matter_Halo_Angular_Momentum(thisNode)/thisBasicComponent%mass()
          call thisHotHaloComponent%angularMomentumSet(angularMomentum   )
          ! Add the appropriate abundances.
          call Halo_Baryonic_Accreted_Abundances(thisNode,accretedAbundances,accretionModeHot)
          call thisHotHaloComponent%     abundancesSet(accretedAbundances)
          ! Also add the appropriate chemical masses.
          call Halo_Baryonic_Accreted_Chemicals (thisNode,accretedChemicals ,accretionModeHot)
          call thisHotHaloComponent%      chemicalsSet(accretedChemicals )
       end if
    end select
    return
  end subroutine Node_Component_Hot_Halo_Standard_Tree_Initialize

  !# <nodeMergerTask>
  !#  <unitName>Node_Component_Hot_Halo_Standard_Node_Merger</unitName>
  !# </nodeMergerTask>
  subroutine Node_Component_Hot_Halo_Standard_Node_Merger(thisNode)
    !% Starve {\tt thisNode} by transferring its hot halo to its parent.
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Dark_Matter_Halo_Scales
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    use Cosmology_Parameters
    use Node_Component_Hot_Halo_Standard_Data
    implicit none
    type            (treeNode                ), intent(inout), pointer :: thisNode
    type            (treeNode                )               , pointer :: parentNode
    class           (nodeComponentHotHalo    )               , pointer :: parentHotHaloComponent , thisHotHaloComponent
    class           (nodeComponentSpin       )               , pointer :: parentSpinComponent
    class           (nodeComponentBasic      )               , pointer :: parentBasic
    class           (cosmologyParametersClass)               , pointer :: thisCosmologyParameters
    double precision                                                   :: baryonicMassCurrent    , baryonicMassMaximum , &
         &                                                                fractionRemove

    ! Get the hot halo component.
    thisHotHaloComponent => thisNode%hotHalo()
    ! Ensure that it is of unspecified class.
    select type (thisHotHaloComponent)
    class is (nodeComponentHotHaloStandard)

       ! Find the parent node and its hot halo and spin components.
       parentNode             => thisNode  %parent
       call Node_Component_Hot_Halo_Standard_Create(parentNode)
       parentHotHaloComponent => parentNode%hotHalo()
       parentSpinComponent    => parentNode%spin   ()

       ! Any gas that failed to be accreted by this halo is always transferred to the parent.
       call parentHotHaloComponent%unaccretedMassSet(                                         &
            &                                         parentHotHaloComponent%unaccretedMass() &
            &                                        +  thisHotHaloComponent%unaccretedMass() &
            &                                       )
       call   thisHotHaloComponent%unaccretedMassSet(                                         &
            &                                         0.0d0                                   &
            &                                       )

       ! Determine if starvation is to be applied.
       if (starveSatellites) then
          ! Move the hot halo to the parent. We leave the hot halo in place even if it is starved, since outflows will accumulate to
          ! this hot halo (and will be moved to the parent at the end of the evolution timestep).
          call parentHotHaloComponent%                    massSet(                                                   &
               &                                                   parentHotHaloComponent%mass                    () &
               &                                                  +  thisHotHaloComponent%mass                    () &
               &                                                 )
          call parentHotHaloComponent%         angularMomentumSet(                                                   &
               &                                                   parentHotHaloComponent%angularMomentum         () &
               &                                                  +  thisHotHaloComponent%mass                    () &
               &                                                  *   parentSpinComponent%spin                    () &
               &                                                  *Dark_Matter_Halo_Virial_Radius  (parentNode)      &
               &                                                  *Dark_Matter_Halo_Virial_Velocity(parentNode)      &
               &                                                 )
          call parentHotHaloComponent%           outflowedMassSet(                                                   &
               &                                                   parentHotHaloComponent%outflowedMass           () &
               &                                                  +  thisHotHaloComponent%outflowedMass           () &
               &                                                 )
          call parentHotHaloComponent%outflowedAngularMomentumSet(                                                   &
               &                                                   parentHotHaloComponent%outflowedAngularMomentum() &
               &                                                  +  thisHotHaloComponent%outflowedMass           () &
               &                                                  *   parentSpinComponent%spin                    () &
               &                                                  *Dark_Matter_Halo_Virial_Radius  (parentNode)      &
               &                                                  *Dark_Matter_Halo_Virial_Velocity(parentNode)      &
               &                                                 )
          call   thisHotHaloComponent%                    massSet(                                                   &
               &                                                   0.0d0                                             &
               &                                                 )
          call   thisHotHaloComponent%         angularMomentumSet(                                                   &
               &                                                   0.0d0                                             &
               &                                                 )
          call   thisHotHaloComponent%           outflowedMassSet(                                                   &
               &                                                   0.0d0                                             &
               &                                                 )
          call   thisHotHaloComponent%outflowedAngularMomentumSet(                                                   &
               &                                                   0.0d0                                             &
               &                                                 )
          call parentHotHaloComponent%              abundancesSet(                                                   &
               &                                                   parentHotHaloComponent%abundances              () &
               &                                                  +  thisHotHaloComponent%abundances              () &
               &                                                 )
          call   thisHotHaloComponent%              abundancesSet(                                                   &
               &                                                   zeroAbundances                                    &
               &                                                 )
          call parentHotHaloComponent%     outflowedAbundancesSet(                                                   &
               &                                                   parentHotHaloComponent%outflowedAbundances     () &
               &                                                  +  thisHotHaloComponent%outflowedAbundances     () &
               &                                                 )
          call   thisHotHaloComponent%     outflowedAbundancesSet(                                                   &
               &                                                   zeroAbundances                                    &
               &                                                 )
          call parentHotHaloComponent%               chemicalsSet(                                                   &
               &                                                   parentHotHaloComponent%chemicals               () &
               &                                                  +  thisHotHaloComponent%chemicals               () &
               &                                                 )
          call   thisHotHaloComponent%               chemicalsSet(                                                   &
               &                                                   zeroChemicals                                     &
               &                                                 )
          ! Check if the baryon fraction in the parent hot halo exceeds the universal value. If it does, mitigate this by moving
          ! some of the mass to the failed accretion reservoir.
          if (hotHaloNodeMergerLimitBaryonFraction) then
             ! Get the default cosmology.
             thisCosmologyParameters => cosmologyParameters()
             parentBasic        => parentNode%basic()
             baryonicMassMaximum=  parentBasic%mass()*thisCosmologyParameters%OmegaBaryon()/thisCosmologyParameters%OmegaMatter()
             baryonicMassCurrent=  Galactic_Structure_Enclosed_Mass(parentNode,radiusLarge,massType=massTypeBaryonic&
                  &,componentType =componentTypeAll)
             if (baryonicMassCurrent > baryonicMassMaximum .and. parentHotHaloComponent%mass() > 0.0d0) then
                fractionRemove=min((baryonicMassCurrent-baryonicMassMaximum)/parentHotHaloComponent%mass(),1.0d0)
                call parentHotHaloComponent% unaccretedMassSet(                                                                &
                     &                                          parentHotHaloComponent%unaccretedMass ()                       &
                     &                                         +parentHotHaloComponent%mass           ()*       fractionRemove &
                     &                                        )
                call parentHotHaloComponent%           massSet( parentHotHaloComponent%mass           ()*(1.0d0-fractionRemove))
                call parentHotHaloComponent%angularMomentumSet( parentHotHaloComponent%angularMomentum()*(1.0d0-fractionRemove))
                call parentHotHaloComponent%     abundancesSet( parentHotHaloComponent%abundances     ()*(1.0d0-fractionRemove))
                call parentHotHaloComponent%     chemicalsSet ( parentHotHaloComponent%chemicals      ()*(1.0d0-fractionRemove))
             end if
          end if
       end if
    end select
    return
  end subroutine Node_Component_Hot_Halo_Standard_Node_Merger

  !# <satelliteMergerTask>
  !#  <unitName>Node_Component_Hot_Halo_Standard_Satellite_Merger</unitName>
  !# </satelliteMergerTask>
  subroutine Node_Component_Hot_Halo_Standard_Satellite_Merger(thisNode)
    !% Remove any hot halo associated with {\tt thisNode} before it merges with its host halo.
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Dark_Matter_Halo_Scales
    use Node_Component_Hot_Halo_Standard_Data
    implicit none
    type (treeNode            ), intent(inout), pointer :: thisNode
    type (treeNode            )               , pointer :: hostNode
    class(nodeComponentHotHalo)               , pointer :: hostHotHaloComponent, thisHotHaloComponent
    class(nodeComponentSpin   )               , pointer :: hostSpinComponent

    ! Return immediately if satellites are starved, as in that case there is no hot halo to transfer.
    if (starveSatellites) return

    ! Get the hot halo component.
    thisHotHaloComponent => thisNode%hotHalo()
    ! Ensure that it is of unspecified class.
    select type (thisHotHaloComponent)
    class is (nodeComponentHotHaloStandard)

       ! Find the node to merge with.
       hostNode             => thisNode%mergesWith()
       hostHotHaloComponent => hostNode%hotHalo   ()
       hostSpinComponent    => hostNode%spin      ()

       ! Move the hot halo to the host.
       call hostHotHaloComponent%                    massSet(                                                 &
            &                                                 hostHotHaloComponent%mass                    () &
            &                                                +thisHotHaloComponent%mass                    () &
            &                                               )
       call hostHotHaloComponent%          unaccretedMassSet(                                                 &
            &                                                 hostHotHaloComponent%unaccretedMass          () &
            &                                                +thisHotHaloComponent%unaccretedMass          () &
            &                                               )
       call hostHotHaloComponent%         angularMomentumSet(                                                 &
            &                                                 hostHotHaloComponent%angularMomentum         () &
            &                                                +thisHotHaloComponent%mass                    () &
            &                                                *   hostSpinComponent%spin                    () &
            &                                                *Dark_Matter_Halo_Virial_Radius  (hostNode)      &
            &                                                *Dark_Matter_Halo_Virial_Velocity(hostNode)      &
            &                                               )
       call hostHotHaloComponent%           outflowedMassSet(                                                 &
            &                                                 hostHotHaloComponent%outflowedMass           () &
            &                                                +thisHotHaloComponent%outflowedMass           () &
            &                                               )
       call hostHotHaloComponent%outflowedAngularMomentumSet(                                                 &
            &                                                 hostHotHaloComponent%outflowedAngularMomentum() &
            &                                                +thisHotHaloComponent%outflowedMass           () &
            &                                                *   hostSpinComponent%spin                    () &
            &                                                *Dark_Matter_Halo_Virial_Radius  (hostNode)      &
            &                                                *Dark_Matter_Halo_Virial_Velocity(hostNode)      &
            &                                               )
       call thisHotHaloComponent%                    massSet(                                                 &
            &                                                 0.0d0                                           &
            &                                               )
       call thisHotHaloComponent%         angularMomentumSet(                                                 &
            &                                                 0.0d0                                           &
            &                                               )
       call thisHotHaloComponent%           outflowedMassSet(                                                 &
            &                                                 0.0d0                                           &
            &                                               )
       call thisHotHaloComponent%outflowedAngularMomentumSet(                                                 &
            &                                                 0.0d0                                           &
            &                                               )
       call hostHotHaloComponent%              abundancesSet(                                                 &
            &                                                 hostHotHaloComponent%abundances              () &
            &                                                +thisHotHaloComponent%abundances              () &
            &                                               )
       call thisHotHaloComponent%              abundancesSet(                                                 &
            &                                                 zeroAbundances                                  &
            &                                               )
       call hostHotHaloComponent%     outflowedAbundancesSet(                                                 &
            &                                                 hostHotHaloComponent%outflowedAbundances     () &
            &                                                +thisHotHaloComponent%outflowedAbundances     () &
            &                                               )
       call thisHotHaloComponent%     outflowedAbundancesSet(                                                 &
            &                                                 zeroAbundances                                  &
            &                                               )
       call hostHotHaloComponent%               chemicalsSet(                                                 &
            &                                                 hostHotHaloComponent%chemicals               () &
            &                                                +thisHotHaloComponent%chemicals               () &
            &                                               )
       call thisHotHaloComponent%               chemicalsSet(                                                 &
            &                                                 zeroChemicals                                   &
            &                                               )
    end select
    return
  end subroutine Node_Component_Hot_Halo_Standard_Satellite_Merger

  !# <nodePromotionTask>
  !#  <unitName>Node_Component_Hot_Halo_Standard_Promote</unitName>
  !# </nodePromotionTask>
  subroutine Node_Component_Hot_Halo_Standard_Promote(thisNode)
    !% Ensure that {\tt thisNode} is ready for promotion to its parent. In this case, we simply update the hot halo mass of {\tt
    !% thisNode} to account for any hot halo already in the parent.
    use Dark_Matter_Halo_Scales
    implicit none
    type (treeNode            ), intent(inout), pointer :: thisNode
    type (treeNode            )               , pointer :: parentNode
    class(nodeComponentHotHalo)               , pointer :: parentHotHaloComponent, thisHotHaloComponent

    ! Get the hot halo component.
    thisHotHaloComponent => thisNode%hotHalo()
    ! Ensure that it is of specified class.
    select type (thisHotHaloComponent)
    class is (nodeComponentHotHaloStandard)
       ! Get the parent node of this node and its hot halo component.
       parentNode             => thisNode  %parent
       parentHotHaloComponent => parentNode%hotHalo()
       ! If the parent node has a hot halo component, then add it to that of this node, and perform other changes needed prior to
       ! promotion.
       select type (parentHotHaloComponent)
       class is (nodeComponentHotHaloStandard)
          call thisHotHaloComponent%                    massSet(                                                   &
               &                                                   thisHotHaloComponent%mass                    () &
               &                                                +parentHotHaloComponent%mass                    () &
               &                                               )
          call thisHotHaloComponent%         angularMomentumSet(                                                   &
               &                                                   thisHotHaloComponent%angularMomentum         () &
               &                                                +parentHotHaloComponent%angularMomentum         () &
               &                                               )
          call thisHotHaloComponent%           outflowedMassSet(                                                   &
               &                                                   thisHotHaloComponent%outflowedMass           () &
               &                                                +parentHotHaloComponent%outflowedMass           () &
               &                                               )
          call thisHotHaloComponent%outflowedAngularMomentumSet(                                                   &
               &                                                   thisHotHaloComponent%outflowedAngularMomentum() &
               &                                                +parentHotHaloComponent%outflowedAngularMomentum() &
               &                                               )
          call thisHotHaloComponent%              abundancesSet(                                                   &
               &                                                   thisHotHaloComponent%abundances              () &
               &                                                +parentHotHaloComponent%abundances              () &
               &                                               )
          call thisHotHaloComponent%     outflowedAbundancesSet(                                                   &
               &                                                   thisHotHaloComponent%outflowedAbundances     () &
               &                                                +parentHotHaloComponent%outflowedAbundances     () &
               &                                               )
          call thisHotHaloComponent%               chemicalsSet(                                                   &
               &                                                   thisHotHaloComponent%chemicals               () &
               &                                                +parentHotHaloComponent%chemicals               () &
               &                                               )
          call thisHotHaloComponent%          unaccretedMassSet(                                                   &
               &                                                   thisHotHaloComponent%unaccretedMass          () &
               &                                                +parentHotHaloComponent%unaccretedMass          () &
               &                                               )
          ! Update the outer radius to match the virial radius of the parent halo.
          call thisHotHaloComponent%             outerRadiusSet(                                                   &
               &                                                Dark_Matter_Halo_Virial_Radius(parentNode)         &
               &                                               )
       end select
    end select
    return
  end subroutine Node_Component_Hot_Halo_Standard_Promote

  subroutine Node_Component_Hot_Halo_Standard_Cooling_Rate(thisNode)
    !% Get and store the cooling rate for {\tt thisNode}.
    use Cooling_Rates
    implicit none
    type (treeNode            ), intent(inout), pointer :: thisNode
    class(nodeComponentHotHalo)               , pointer :: thisHotHaloComponent

    if (.not.gotCoolingRate) then
       ! Get the hot halo component.
       thisHotHaloComponent => thisNode%hotHalo()
       if     (                                                &
            &   thisHotHaloComponent%mass           () > 0.0d0 &
            &  .and.                                           &
            &   thisHotHaloComponent%angularMomentum() > 0.0d0 &
            &  .and.                                           &
            &   thisHotHaloComponent%outerRadius    () > 0.0d0 &
            & ) then
          ! Get the cooling rate.
          coolingRate=Cooling_Rate(thisNode)
       else
          coolingRate=0.0d0
       end if

       ! Store a copy of this cooling rate as the remaining mass heating rate budget. This is used to ensure that we never heat
       ! gas at a rate greater than it is cooling.
       massHeatingRateRemaining=coolingRate

       ! Flag that cooling rate has now been computed.
       gotCoolingRate=.true.
    end if
    return
  end subroutine Node_Component_Hot_Halo_Standard_Cooling_Rate

  subroutine Node_Component_Hot_Halo_Standard_Create(thisNode)
    !% Creates a hot halo component for {\tt thisNode}.
    implicit none
    type (treeNode            ), intent(inout), pointer :: thisNode
    class(nodeComponentHotHalo)               , pointer :: thisHotHalo

    ! Create the component.
    thisHotHalo => thisNode%hotHalo(autoCreate=.true.)
    ! Initalize.
    select type (thisHotHalo)
    class is (nodeComponentHotHaloStandard)
       call Node_Component_Hot_Halo_Standard_Initializor(thisHotHalo)
    end select
    return
  end subroutine Node_Component_Hot_Halo_Standard_Create

  subroutine Node_Component_Hot_Halo_Standard_Initializor(self)
    !% Initializes a standard hot halo component.
    use Dark_Matter_Halo_Scales
    implicit none
    type(nodeComponentHotHaloStandard)          :: self
    type(treeNode                    ), pointer :: selfNode

    ! Return if already initialized.
    if (self%isInitialized()) return
    ! Get the hosting node.
    selfNode => self%hostNode
    ! Initialize the outer boundary to the virial radius.
    call self%outerRadiusSet(Dark_Matter_Halo_Virial_Radius(selfNode))
    ! Record that the spheroid has been initialized.
    call self%isInitializedSet(.true.)
    return
  end subroutine Node_Component_Hot_Halo_Standard_Initializor

  !# <haloFormationTask>
  !#  <unitName>Node_Component_Hot_Halo_Standard_Formation</unitName>
  !# </haloFormationTask>
  subroutine Node_Component_Hot_Halo_Standard_Formation(thisNode)
    !% Updates the hot halo gas distribution at a formation event, if requested.
    use Chemical_States
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Chemical_Reaction_Rates_Utilities
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Astronomical
    use Node_Component_Hot_Halo_Standard_Data
    implicit none
    type            (treeNode            ), intent(inout), pointer :: thisNode
    class           (nodeComponentHotHalo)               , pointer :: thisHotHaloComponent
    type            (abundances          ), save                   :: outflowedAbundances
    type            (chemicalAbundances  ), save                   :: chemicalDensities    , chemicalMasses
    !$omp threadprivate(outflowedAbundances,chemicalDensities,chemicalMasses)
    double precision                                               :: hydrogenByMass       , massToDensityConversion, &
         &                                                            numberDensityHydrogen, temperature

    ! Return immediately if return of outflowed gas on formation events is not requested.
    if (.not.hotHaloOutflowReturnOnFormation) return

    ! Get the hot halo component.
    thisHotHaloComponent => thisNode%hotHalo()
    ! Ensure that it is of unspecified class.
    select type (thisHotHaloComponent)
    class is (nodeComponentHotHaloStandard)

       ! Compute mass of chemicals transferred to the hot halo.
       if (chemicalsCount > 0 .and. thisHotHaloComponent%outflowedMass() > 0.0d0) then
          ! Compute coefficient in conversion of mass to density for this node.
          massToDensityConversion=Chemicals_Mass_To_Density_Conversion(Dark_Matter_Halo_Virial_Radius(thisNode))/3.0d0
          ! Get abundance mass fractions of the outflowed material.
          outflowedAbundances= thisHotHaloComponent%outflowedAbundances() &
               &              /thisHotHaloComponent%outflowedMass      ()
          ! Get the hydrogen mass fraction in outflowed gas.
          hydrogenByMass=outflowedAbundances%hydrogenMassFraction()
          ! Compute the temperature and density of material in the hot halo.
          temperature          =Dark_Matter_Halo_Virial_Temperature(thisNode)
          numberDensityHydrogen=hydrogenByMass*thisHotHaloComponent%outflowedMass()*massToDensityConversion/atomicMassHydrogen
          ! Set the radiation field.
          call radiation%set(thisNode)
          ! Get the chemical densities.
          call Chemical_Densities(chemicalDensities,temperature,numberDensityHydrogen,outflowedAbundances,radiation)
          ! Convert from densities to masses.
          call chemicalDensities%numberToMass(chemicalMasses)
          chemicalMasses=chemicalMasses*thisHotHaloComponent%outflowedMass()*hydrogenByMass/numberDensityHydrogen/atomicMassHydrogen
          ! Add chemicals to the hot component.
          call thisHotHaloComponent%chemicalsSet(                                  &
               &                                  thisHotHaloComponent%chemicals() &
               &                                 +chemicalMasses                   &
               &                                )
       end if
       ! Transfer mass, angular momentum and abundances.
       call thisHotHaloComponent%                    massSet(&
            &                                                 thisHotHaloComponent%         mass           () &
            &                                                +thisHotHaloComponent%outflowedMass           () &
            &)
       call thisHotHaloComponent%         angularMomentumSet(                                                 &
            &                                                 thisHotHaloComponent%         angularMomentum() &
            &                                                +thisHotHaloComponent%outflowedAngularMomentum() &
            &                                               )
       call thisHotHaloComponent%              abundancesSet(                                                 &
            &                                                 thisHotHaloComponent%         abundances     () &
            &                                                +thisHotHaloComponent%outflowedAbundances     () &
            &                                               )
       call thisHotHaloComponent%           outflowedMassSet(                                                 &
            &                                                 0.0d0                                           &
            &                                               )
       call thisHotHaloComponent%outflowedAngularMomentumSet(                                                 &
            &                                                 0.0d0                                           &
            &                                               )
       call thisHotHaloComponent%     outflowedAbundancesSet(                                                 &
            &                                                 zeroAbundances                                  &
            &                                               )
    end select

    return
  end subroutine Node_Component_Hot_Halo_Standard_Formation

end module Node_Component_Hot_Halo_Standard
