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

!% Contains a module which implements the simple black hole node component.

module Node_Component_Black_Hole_Simple
  !% Implements the simple black hole node component.
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Black_Hole_Simple_Initialize       , Node_Component_Black_Hole_Simple_Scale_Set              , &
       &    Node_Component_Black_Hole_Simple_Satellite_Merging, Node_Component_Black_Hole_Simple_Output_Names           , &
       &    Node_Component_Black_Hole_Simple_Output_Count     , Node_Component_Black_Hole_Simple_Output                 , &
       &    Node_Component_Black_Hole_Simple_Rate_Compute

  !# <component>
  !#  <class>blackHole</class>
  !#  <name>simple</name>
  !#  <isDefault>no</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>mass</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <classDefault>defaultBlackHoleComponent%massSeed()</classDefault>
  !#     <output unitsInSI="massSolar" comment="Mass of the black hole."/>
  !#   </property>
  !#   <property>
  !#     <name>massSeed</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="false" isGettable="true" isEvolvable="false" />
  !#     <getFunction>Node_Component_Black_Hole_Simple_Seed_Mass</getFunction>
  !#     <isVirtual>true</isVirtual>
  !#   </property>
  !#  </properties>
  !#  <bindings>
  !#     <binding method="enclosedMass" function="Node_Component_Black_Hole_Simple_Enclosed_Mass" bindsTo="component" />
  !#  </bindings>
  !#  <functions>objects.nodes.components.black_hole.simple.bound_functions.inc</functions>
  !# </component>

  ! Seed mass for black holes.
  double precision :: blackHoleSeedMass

  ! Feedback parameters.
  double precision :: blackHoleHeatingEfficiency                   , blackHoleJetEfficiency , &
       &              blackHoleToSpheroidStellarGrowthRatio        , blackHoleWindEfficiency
  logical          :: blackHoleAccretesFromHotHalo                 , blackHoleHeatsHotHalo

  ! Output options.
  logical          :: blackHoleOutputAccretion

  ! Record of whether this module has been initialized.
  logical          :: moduleInitialized                    =.false.

contains

  subroutine Node_Component_Black_Hole_Simple_Initialize()
    !% Initializes the simple black hole node component module.
    use Input_Parameters
    implicit none
    type(nodeComponentBlackHoleSimple) :: blackHoleSimple

    ! Initialize the module if necessary.
    !$omp critical (Node_Component_Black_Hole_Simple_Initialize)
    if (.not.moduleInitialized) then
       ! Get the black hole seed mass.
       blackHoleSeedMass=blackHoleSimple%massSeed()
       ! Get accretion rate enhancement factors.
       !@ <inputParameter>
       !@   <name>blackHoleToSpheroidStellarGrowthRatio</name>
       !@   <defaultValue>1.0d-3</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The ratio of the rates of black hole growth and spheroid stellar mass growth.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>blackHoles</group>
       !@ </inputParameter>
       call Get_Input_Parameter("blackHoleToSpheroidStellarGrowthRatio",blackHoleToSpheroidStellarGrowthRatio,defaultValue&
            &=1.0d-3)
       ! Options controlling AGN feedback.
       !@ <inputParameter>
       !@   <name>blackHoleHeatsHotHalo</name>
       !@   <defaultValue>true</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether or not the black hole should heat the hot halo.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@   <group>blackHoles</group>
       !@ </inputParameter>
       call Get_Input_Parameter("blackHoleHeatsHotHalo",blackHoleHeatsHotHalo,defaultValue=.true.)
       if (blackHoleHeatsHotHalo) then
          !@ <inputParameter>
          !@   <name>blackHoleAccretesFromHotHalo</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Controls whether the black hole additionally grows via accretion from the hot halo. If it does,
          !@     this accretion rate is used to determine AGN feedback power.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>blackHoles</group>
          !@ </inputParameter>
          call Get_Input_Parameter("blackHoleAccretesFromHotHalo",blackHoleAccretesFromHotHalo,defaultValue=.false.)
          !@ <inputParameter>
          !@   <name>blackHoleHeatingEfficiency</name>
          !@   <defaultValue>$10^{-3}$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The efficiency with which accretion onto a black hole heats the hot halo.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>blackHoles</group>
          !@ </inputParameter>
          call Get_Input_Parameter("blackHoleHeatingEfficiency",blackHoleHeatingEfficiency,defaultValue=1.0d-3)
          !@ <inputParameter>
          !@   <name>blackHoleJetEfficiency</name>
          !@   <defaultValue>$10^{-3}$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The efficiency with which accretion power onto a black hole is converted into jets.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>blackHoles</group>
          !@ </inputParameter>
          call Get_Input_Parameter("blackHoleJetEfficiency",blackHoleJetEfficiency,defaultValue=1.0d-3)
       else
          blackHoleHeatingEfficiency=0.0d0
          blackHoleJetEfficiency    =0.0d0
       end if
       ! Get options controlling winds.
       !@ <inputParameter>
       !@   <name>blackHoleWindEfficiency</name>
       !@   <defaultValue>$2.2157\times 10^{-3}$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The efficiency of the black hole accretion-driven wind.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>blackHoles</group>
       !@ </inputParameter>
       call Get_Input_Parameter("blackHoleWindEfficiency",blackHoleWindEfficiency,defaultValue=2.2157d-3)
       ! Get options controlling output.
       !@ <inputParameter>
       !@   <name>blackHoleOutputAccretion</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Determines whether or not accretion rates and jet powers will be output.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("blackHoleOutputAccretion",blackHoleOutputAccretion,defaultValue=.false.)
       ! Record that the module is now initialized.
       moduleInitialized=.true.
    end if
    !$omp end critical (Node_Component_Black_Hole_Simple_Initialize)
    return
  end subroutine Node_Component_Black_Hole_Simple_Initialize

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Black_Hole_Simple_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Black_Hole_Simple_Scale_Set(thisNode)
    !% Set scales for properties of {\tt thisNode}.
    implicit none
    type (treeNode              ), intent(inout), pointer :: thisNode
    class(nodeComponentBlackHole)               , pointer :: thisBlackHoleComponent
    class(nodeComponentSpheroid )               , pointer :: thisSpheroidComponent

    ! Get the black hole component.
    thisBlackHoleComponent => thisNode%blackHole()
    ! Ensure that it is of the standard class.
    select type (thisBlackHoleComponent)
    class is (nodeComponentBlackHoleSimple)
       ! Get the spheroid component.
       thisSpheroidComponent => thisNode%spheroid()
       ! Set scale for mass.
       call thisBlackHoleComponent%massScale(                                                                                &
            &                                max(                                                                            &
            &                                    thisSpheroidComponent %massStellar()*blackHoleToSpheroidStellarGrowthRatio, &
            &                                    thisBlackHoleComponent%mass       ()                                        &
            &                                   )                                                                            &
            &                               )
    end select
    return
  end subroutine Node_Component_Black_Hole_Simple_Scale_Set

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Black_Hole_Simple_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Black_Hole_Simple_Rate_Compute(thisNode,interrupt,interruptProcedure)
    !% Compute the black hole mass rate of change.
    use Cooling_Radii
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Physical
    implicit none
    type            (treeNode                                          )           , intent(inout), pointer :: thisNode
    logical                                                                        , intent(inout)          :: interrupt
    procedure       (Interrupt_Procedure_Template                      )           , intent(inout), pointer :: interruptProcedure
    class           (nodeComponentBlackHole                            )                          , pointer :: thisBlackHoleComponent
    class           (nodeComponentSpheroid                             )                          , pointer :: thisSpheroidComponent
    class           (nodeComponentHotHalo                              )                          , pointer :: thisHotHaloComponent
    class           (darkMatterHaloScaleClass)               , pointer :: darkMatterHaloScale_
    double precision                                                    , parameter                         :: coolingRadiusFractionalTransitionMinimum=0.9d0
    double precision                                                    , parameter                         :: coolingRadiusFractionalTransitionMaximum=1.0d0
    double precision                                                                                        :: coolingRadiusFractional                       , couplingEfficiency   , &
         &                                                                                                     energyInputRate                               , heatingRate          , &
         &                                                                                                     massAccretionRate                             , restMassAccretionRate, &
         &                                                                                                     x

    if (defaultBlackHoleComponent%simpleIsActive()) then

       ! Get the spheroid component.
       thisSpheroidComponent => thisNode%spheroid()

       ! Find the rate of rest mass accretion onto the black hole.
       restMassAccretionRate=blackHoleToSpheroidStellarGrowthRatio*thisSpheroidComponent%starFormationRate()

       ! Finish if there is no accretion.
       if (restMassAccretionRate <= 0.0d0) return

       ! Find the rate of increase in mass of the black hole.
       massAccretionRate=restMassAccretionRate*max((1.0d0-blackHoleHeatingEfficiency-blackHoleWindEfficiency),0.0d0)

       ! Get the black hole component.
       thisBlackHoleComponent => thisNode%blackHole()

       ! Detect black hole component type.
       select type (thisBlackHoleComponent)
       type is (nodeComponentBlackHole)
          ! Generic type - interrupt and create a simple black hole if accretion rate is non-zero.
          if (massAccretionRate /= 0.0d0) then
             interrupt=.true.
             interruptProcedure => Node_Component_Black_Hole_Simple_Create
          end if
          return
       class is (nodeComponentBlackHoleSimple)
          ! Simple type - continue processing.
          call thisBlackHoleComponent%massRate       (     massAccretionRate)
          ! Remove the accreted mass from the spheroid component.
          call thisSpheroidComponent %massGasSinkRate(-restMassAccretionRate)
          ! Add heating to the hot halo component.
          if (blackHoleHeatsHotHalo) then
             ! Compute jet coupling efficiency based on whether halo is cooling quasistatically.
             darkMatterHaloScale_ => darkMatterHaloScale()
             coolingRadiusFractional=Cooling_Radius(thisNode)/darkMatterHaloScale_%virialRadius(thisNode)
             if      (coolingRadiusFractional < coolingRadiusFractionalTransitionMinimum) then
                couplingEfficiency=1.0d0
             else if (coolingRadiusFractional > coolingRadiusFractionalTransitionMaximum) then
                couplingEfficiency=0.0d0
             else
                x=      (coolingRadiusFractional                 -coolingRadiusFractionalTransitionMinimum) &
                     & /(coolingRadiusFractionalTransitionMaximum-coolingRadiusFractionalTransitionMinimum)
                couplingEfficiency=x**2*(2.0d0*x-3.0d0)+1.0d0
             end if
             ! Compute the heating rate.
             heatingRate=couplingEfficiency*blackHoleHeatingEfficiency*restMassAccretionRate*(speedLight/kilo)**2
             ! Pipe this power to the hot halo.
             thisHotHaloComponent => thisNode%hotHalo()
             call thisHotHaloComponent%heatSourceRate(heatingRate,interrupt,interruptProcedure)
          end if
          ! Add energy to the spheroid component.
          if (blackHoleWindEfficiency > 0.0d0) then
             ! Compute the energy input and send it down the spheroid gas energy input pipe.
             energyInputRate=blackHoleWindEfficiency*restMassAccretionRate*(speedLight/kilo)**2
             call thisSpheroidComponent%energyGasInputRate(energyInputRate)
          end if
       end select
    end if
    return
  end subroutine Node_Component_Black_Hole_Simple_Rate_Compute

  !# <satelliteMergerTask>
  !#  <unitName>Node_Component_Black_Hole_Simple_Satellite_Merging</unitName>
  !# </satelliteMergerTask>
  subroutine Node_Component_Black_Hole_Simple_Satellite_Merging(thisNode)
    !% Merge (instantaneously) any simple black hole associated with {\tt thisNode} before it merges with its host halo.
    use Black_Hole_Binary_Mergers
    implicit none
    type            (treeNode              ), intent(inout), pointer :: thisNode
    type            (treeNode              )               , pointer :: hostNode
    class           (nodeComponentBlackHole)               , pointer :: hostBlackHoleComponent, thisBlackHoleComponent
    double precision                                                 :: blackHoleMassNew      , blackHoleSpinNew

    ! Check that the simple black hole is active.
    if (defaultBlackHoleComponent%simpleIsActive()) then
       ! Find the node to merge with.
       hostNode               => thisNode%mergesWith()
       ! Get the black holes.
       thisBlackHoleComponent => thisNode%blackHole ()
       hostBlackHoleComponent => hostNode%blackHole ()
       ! Compute the effects of the merger.
       call Black_Hole_Binary_Merger(thisBlackHoleComponent%mass(), &
            &                        hostBlackHoleComponent%mass(), &
            &                        0.0d0                        , &
            &                        0.0d0                        , &
            &                        blackHoleMassNew             , &
            &                        blackHoleSpinNew               &
            &                       )
       ! Move the black hole to the host.
       call hostBlackHoleComponent%massSet(blackHoleMassNew)
       call thisBlackHoleComponent%massSet(           0.0d0)
    end if
    return
  end subroutine Node_Component_Black_Hole_Simple_Satellite_Merging

  subroutine Node_Component_Black_Hole_Simple_Create(thisNode)
    !% Creates a simple black hole component for {\tt thisNode}.
    implicit none
    type (treeNode              ), intent(inout), pointer :: thisNode
    class(nodeComponentBlackHole)               , pointer :: thisBlackHoleComponent

    ! Create the component.
    thisBlackHoleComponent => thisNode%blackHole(autoCreate=.true.)
    ! Set the seed mass.
    call thisBlackHoleComponent%massSet(blackHoleSeedMass)
    return
  end subroutine Node_Component_Black_Hole_Simple_Create

  !# <mergerTreeOutputNames>
  !#  <unitName>Node_Component_Black_Hole_Simple_Output_Names</unitName>
  !#  <sortName>Node_Component_Black_Hole_Simple_Output</sortName>
  !# </mergerTreeOutputNames>
  subroutine Node_Component_Black_Hole_Simple_Output_Names(thisNode,integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI&
       &,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set names of black hole properties to be written to the \glc\ output file.
    use Numerical_Constants_Astronomical
    implicit none
    type            (treeNode)              , intent(inout), pointer :: thisNode
    double precision                        , intent(in   )          :: time
    integer                                 , intent(inout)          :: doubleProperty         , integerProperty
    character       (len=*   ), dimension(:), intent(inout)          :: doublePropertyComments , doublePropertyNames   , &
         &                                                              integerPropertyComments, integerPropertyNames
    double precision          , dimension(:), intent(inout)          :: doublePropertyUnitsSI  , integerPropertyUnitsSI

    ! Ensure that the black hole component is of the simple class.
    if (Node_Component_Black_Hole_Simple_Matches(thisNode)) then
       !@ <outputPropertyGroup>
       !@   <name>blackHole</name>
       !@   <description>Black hole properities</description>
       !@   <outputType>nodeData</outputType>
       !@ </outputPropertyGroup>
       if (blackHoleOutputAccretion) then
          doubleProperty=doubleProperty+1
          !@ <outputProperty>
          !@   <name>blackHoleAccretionRate</name>
          !@   <datatype>real</datatype>
          !@   <cardinality>0..1</cardinality>
          !@   <description>Rest-mass accretion rate onto the black hole.</description>
          !@   <label>???</label>
          !@   <outputType>nodeData</outputType>
          !@ </outputProperty>
          doublePropertyNames   (doubleProperty)='blackHoleAccretionRate'
          doublePropertyComments(doubleProperty)='Rest-mass accretion rate onto the black hole.'
          doublePropertyUnitsSI (doubleProperty)=massSolar/gigaYear
       end if
    end if
    return
  end subroutine Node_Component_Black_Hole_Simple_Output_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Node_Component_Black_Hole_Simple_Output_Count</unitName>
  !#  <sortName>Node_Component_Black_Hole_Simple_Output</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Node_Component_Black_Hole_Simple_Output_Count(thisNode,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of black hole properties to be written to the the \glc\ output file.
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode
    double precision          , intent(in   )          :: time
    integer                   , intent(inout)          :: doublePropertyCount  , integerPropertyCount
    integer                   , parameter              :: extraPropertyCount =1

    ! Ensure that the black hole component is of the simple class.
    if (Node_Component_Black_Hole_Simple_Matches(thisNode)) then
       if (blackHoleOutputAccretion) doublePropertyCount=doublePropertyCount+extraPropertyCount
    end if
    return
  end subroutine Node_Component_Black_Hole_Simple_Output_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Node_Component_Black_Hole_Simple_Output</unitName>
  !#  <sortName>Node_Component_Black_Hole_Simple_Output</sortName>
  !# </mergerTreeOutputTask>
  subroutine Node_Component_Black_Hole_Simple_Output(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store black hole properties in the \glc\ output file buffers.
    use Kind_Numbers
    implicit none
    double precision                        , intent(in   )          :: time
    type            (treeNode              ), intent(inout), pointer :: thisNode
    integer                                 , intent(inout)          :: doubleBufferCount          , doubleProperty, integerBufferCount, &
         &                                                              integerProperty
    integer         (kind=kind_int8        ), intent(inout)          :: integerBuffer         (:,:)
    double precision                        , intent(inout)          :: doubleBuffer          (:,:)
    class           (nodeComponentBlackHole)               , pointer :: thisBlackHoleComponent
    class           (nodeComponentSpheroid )               , pointer :: thisSpheroidComponent
    double precision                                                 :: restMassAccretionRate

    ! Ensure that the black hole component is of the simple class.
    if (Node_Component_Black_Hole_Simple_Matches(thisNode)) then
       ! Get the black hole component.
       thisBlackHoleComponent => thisNode%blackHole()
       ! Store the properties.
       if (blackHoleOutputAccretion) then
          ! Get the rest mass accretion rate.
          thisSpheroidComponent => thisNode%spheroid()
          restMassAccretionRate=blackHoleToSpheroidStellarGrowthRatio*thisSpheroidComponent%starFormationRate()
          doubleProperty=doubleProperty+1
          doubleBuffer(doubleBufferCount,doubleProperty)=restMassAccretionRate
       end if
    end if
    return
  end subroutine Node_Component_Black_Hole_Simple_Output

  logical function Node_Component_Black_Hole_Simple_Matches(thisNode)
    !% Return true if the black hole component of {\tt thisNode} is a match to the simple implementation.
    implicit none
    type (treeNode              ), intent(inout), pointer :: thisNode
    class(nodeComponentBlackHole)               , pointer :: thisBlackHoleComponent

    ! Get the black hole component.
    thisBlackHoleComponent => thisNode%blackHole()
    ! Ensure that it is of the simple class.
    Node_Component_Black_Hole_Simple_Matches=.false.
    select type (thisBlackHoleComponent)
    class is (nodeComponentBlackHoleSimple)
       Node_Component_Black_Hole_Simple_Matches=.true.
    type  is (nodeComponentBlackHole       )
       Node_Component_Black_Hole_Simple_Matches=defaultBlackHoleComponent%simpleIsActive()
    end select
    return
  end function Node_Component_Black_Hole_Simple_Matches

end module Node_Component_Black_Hole_Simple
