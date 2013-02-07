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

!% Contains a module which implements the simple black hole node component.

module Node_Component_Black_Hole_Simple
  !% Implements the simple black hole node component.
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Black_Hole_Simple_Initialize       , Node_Component_Black_Hole_Simple_Scale_Set              , &
       &    Node_Component_Black_Hole_Simple_Satellite_Merging, Node_Component_Black_Hole_Simple_Output_Names           , &
       &    Node_Component_Black_Hole_Simple_Output_Count     , Node_Component_Black_Hole_Simple_Output                 , &
       &    Node_Component_Black_Hole_Simple_Enclosed_Mass    , Node_Component_Black_Hole_Simple_Rotation_Curve         , &
       &    Node_Component_Black_Hole_Simple_Potential        , Node_Component_Black_Hole_Simple_Rotation_Curve_Gradient, &
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
  !#     <attributes isSettable="false" isGettable="true" isEvolvable="false" isVirtual="true" />
  !#     <getFunction>Node_Component_Black_Hole_Simple_Seed_Mass</getFunction>
  !#   </property>
  !#  </properties>
  !#  <functions>objects.nodes.components.black_hole.simple.bound_functions.inc</functions>
  !# </component>
  
  ! Seed mass for black holes.
  double precision :: blackHoleSeedMass

  ! Feedback parameters.
  double precision :: blackHoleToSpheroidStellarGrowthRatio,blackHoleWindEfficiency,blackHoleHeatingEfficiency,blackHoleJetEfficiency
  logical          :: blackHoleHeatsHotHalo,blackHoleAccretesFromHotHalo

  ! Output options.
  logical          :: blackHoleOutputAccretion

  ! Record of whether this module has been initialized.
  logical          :: moduleInitialized=.false.

contains

  subroutine Node_Component_Black_Hole_Simple_Initialize()
    !% Initializes the simple black hole node component module.
    use ISO_Varying_String
    use Input_Parameters
    use String_Handling
    use Galacticus_Display
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
    type (treeNode              ), pointer, intent(inout) :: thisNode
    class(nodeComponentBlackHole), pointer                :: thisBlackHoleComponent
    class(nodeComponentSpheroid ), pointer                :: thisSpheroidComponent

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
    use Accretion_Disks
    use Black_Hole_Fundamentals   
    use Cooling_Radii
    use Cooling_Rates
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Physical
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Math
    use Numerical_Constants_Astronomical
    implicit none
    type     (treeNode              ), pointer, intent(inout) :: thisNode
    logical                          ,          intent(inout) :: interrupt
    procedure(                      ), pointer, intent(inout) :: interruptProcedure
    class    (nodeComponentBlackHole), pointer                :: thisBlackHoleComponent
    class    (nodeComponentSpheroid ), pointer                :: thisSpheroidComponent
    class    (nodeComponentHotHalo  ), pointer                :: thisHotHaloComponent
    double precision                 , parameter              :: coolingRadiusFractionalTransitionMinimum=0.9d0
    double precision                 , parameter              :: coolingRadiusFractionalTransitionMaximum=1.0d0
    double precision                                          :: restMassAccretionRate,massAccretionRate,energyInputRate&
         &,heatingRate ,couplingEfficiency,coolingRadiusFractional,x

    if (defaultBlackHoleComponent%simpleIsActive()) then

       ! Get the spheroid component.
       thisSpheroidComponent => thisNode%spheroid()

       ! Find the rate of rest mass accretion onto the black hole.
       restMassAccretionRate=blackHoleToSpheroidStellarGrowthRatio*thisSpheroidComponent%starFormationRate()

       ! Finish if there is no accretion.
       if (restMassAccretionRate <= 0.0d0) return

       ! Find the rate of increase in mass of the black hole.
       massAccretionRate=restMassAccretionRate*(1.0d0-blackHoleHeatingEfficiency-blackHoleWindEfficiency)

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
             coolingRadiusFractional=Cooling_Radius(thisNode)/Dark_Matter_Halo_Virial_Radius(thisNode)
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
    type (treeNode              ), pointer, intent(inout) :: thisNode
    type (treeNode              ), pointer                :: hostNode
    class(nodeComponentBlackHole), pointer                :: thisBlackHoleComponent,hostBlackHoleComponent
    double precision                                      :: blackHoleMassNew,blackHoleSpinNew

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
    type (treeNode              ), pointer, intent(inout) :: thisNode
    class(nodeComponentBlackHole), pointer                :: thisBlackHoleComponent

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
    use ISO_Varying_String
    implicit none
    type            (treeNode), intent(inout), pointer      :: thisNode
    double precision          , intent(in   )               :: time
    integer                   , intent(inout)               :: integerProperty,doubleProperty
    character       (len=*   ), intent(inout), dimension(:) :: integerPropertyNames,integerPropertyComments,doublePropertyNames &
         &,doublePropertyComments
    double precision          , intent(inout), dimension(:) :: integerPropertyUnitsSI,doublePropertyUnitsSI
    
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
    integer                   , intent(inout)          :: integerPropertyCount,doublePropertyCount
    integer                   , parameter              :: extraPropertyCount=1

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
    use Galacticus_Nodes
    use Kind_Numbers
    use Accretion_Disks
    implicit none
    double precision                        , intent(in   )          :: time
    type            (treeNode              ), intent(inout), pointer :: thisNode
    integer                                 , intent(inout)          :: integerProperty,integerBufferCount,doubleProperty&
         &,doubleBufferCount
    integer         (kind=kind_int8        ), intent(inout)          :: integerBuffer(:,:)
    double precision                        , intent(inout)          :: doubleBuffer (:,:)
    class           (nodeComponentBlackHole),                pointer :: thisBlackHoleComponent
    class           (nodeComponentSpheroid ),                pointer :: thisSpheroidComponent
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
    class(nodeComponentBlackHole),                pointer :: thisBlackHoleComponent

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

  !# <enclosedMassTask>
  !#  <unitName>Node_Component_Black_Hole_Simple_Enclosed_Mass</unitName>
  !# </enclosedMassTask>
  double precision function Node_Component_Black_Hole_Simple_Enclosed_Mass(thisNode,radius,massType,componentType,weightBy,weightIndex,haloLoaded)
    !% Computes the mass within a given radius for a central black hole. Black hole is treated as a point mass.
    use Galactic_Structure_Options
    implicit none
    type            (treeNode              ), intent(inout), pointer  :: thisNode
    integer                                 , intent(in   )           :: massType,componentType,weightBy,weightIndex
    double precision                        , intent(in   )           :: radius
    logical                                 , intent(in   ), optional :: haloLoaded
    class           (nodeComponentBlackHole),                pointer  :: thisBlackHoleComponent

    ! Set zero enclosed mass by default.
    Node_Component_Black_Hole_Simple_Enclosed_Mass=0.0d0

    ! Return the black hole mass only if massType and componentType are of black hole type.
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeBlackHole)) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeBlackHole     )) return

    ! Get the spheroid component and check that it is of the Hernquist class.
    thisBlackHoleComponent => thisNode%blackHole()
    select type (thisBlackHoleComponent)
    class is (nodeComponentBlackHoleSimple)
       ! Compute the enclosed mass.
       if (radius >= 0.0d0) then
          select case (weightBy)
          case (weightByMass)
             Node_Component_Black_Hole_Simple_Enclosed_Mass=thisBlackHoleComponent%mass()
          end select
       end if
    end select
    return
  end function Node_Component_Black_Hole_Simple_Enclosed_Mass

  !# <rotationCurveTask>
  !#  <unitName>Node_Component_Black_Hole_Simple_Rotation_Curve</unitName>
  !# </rotationCurveTask>
  double precision function Node_Component_Black_Hole_Simple_Rotation_Curve(thisNode,radius,massType,componentType,haloLoaded)
    !% Computes the rotation curve for the central black hole. Assumes a point mass black hole with a Keplerian rotation curve,
    !% \emph{except} that the rotation speed is limited to never exceed the speed of light.
    use Galacticus_Nodes
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    use Numerical_Constants_Prefixes
    use Black_Hole_Fundamentals
    implicit none
    type            (treeNode              ), intent(inout), pointer  :: thisNode
    integer                                 , intent(in   )           :: massType,componentType
    double precision                        , intent(in   )           :: radius
    logical                                 , intent(in   ), optional :: haloLoaded
    class           (nodeComponentBlackHole),                pointer  :: thisBlackHoleComponent
    double precision                                                  :: componentMass

    ! Set to zero by default.
    Node_Component_Black_Hole_Simple_Rotation_Curve=0.0d0
    ! Get the black hole component and check that it is of the simple class.
    thisBlackHoleComponent => thisNode%blackHole()
    select type (thisBlackHoleComponent)
    class is (nodeComponentBlackHoleSimple)
       ! Check if the radius exceeds the gravitational radius.
       if (radius > Black_Hole_Gravitational_Radius(thisBlackHoleComponent)/(milli*speedLight)**2) then
          ! Radius is larger than the gravitational radius - compute the rotation speed.
          componentMass=Node_Component_Black_Hole_Simple_Enclosed_Mass(thisNode,radius,massType,componentType,weightByMass&
               &,weightIndexNull,haloLoaded)
          if (componentMass > 0.0d0) Node_Component_Black_Hole_Simple_Rotation_Curve=dsqrt(gravitationalConstantGalacticus&
               &*componentMass/radius)
       else
          ! Radius is less than the gravitational radius - return the speed of light.
          Node_Component_Black_Hole_Simple_Rotation_Curve=speedLight*milli
       end if
    end select
    return
  end function Node_Component_Black_Hole_Simple_Rotation_Curve

  !# <rotationCurveGradientTask>
  !#  <unitName>Node_Component_Black_Hole_Simple_Rotation_Curve_Gradient</unitName>
  !# </rotationCurveGradientTask>
  double precision function Node_Component_Black_Hole_Simple_Rotation_Curve_Gradient(thisNode,radius,massType,componentType&
       &,haloLoaded)
    !% Computes the rotation curve gradient for the central black hole. Assumes a point mass black hole with a Keplerian 
    !% rotation curve, \emph{except} that the rotation speed is limited to never exceed the speed of light.
    use Galacticus_Nodes
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    use Numerical_Constants_Prefixes
    use Black_Hole_Fundamentals
    implicit none
    type            (treeNode              ), intent(inout), pointer  :: thisNode
    integer                                 , intent(in   )           :: massType,componentType
    double precision                        , intent(in   )           :: radius
    logical                                 , intent(in   ), optional :: haloLoaded
    class           (nodeComponentBlackHole),                pointer  :: thisBlackHoleComponent
    double precision                                                  :: componentMass

    ! Set to zero by default.
    Node_Component_Black_Hole_Simple_Rotation_Curve_Gradient=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeBlackHole)) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeBlackHole     )) return  
    if (      radius        <= 0.0d0                                                         ) return
    ! Get the black hole component and check that it is of the simple class.
    thisBlackHoleComponent => thisNode%blackHole()
    select type (thisBlackHoleComponent)
    class is (nodeComponentBlackHoleSimple)
       componentMass=Node_Component_Black_Hole_Simple_Enclosed_Mass(thisNode,radius,massType,componentType,weightByMass,weightIndexNull&
            &,haloLoaded)
       if (componentMass ==0.0d0 ) return
       if (radius > Black_Hole_Gravitational_Radius(thisBlackHoleComponent)) then
          Node_Component_Black_Hole_Simple_Rotation_Curve_Gradient=       &
               &                         -gravitationalConstantGalacticus &
               &                         *componentMass                   &
               &                         /radius**2
       else
          Node_Component_Black_Hole_Simple_Rotation_Curve_Gradient=0.0d0
       end if
    end select
    return
  end function Node_Component_Black_Hole_Simple_Rotation_Curve_Gradient

  !# <potentialTask>
  !#  <unitName>Node_Component_Black_Hole_Simple_Potential</unitName>
  !# </potentialTask>
  double precision function Node_Component_Black_Hole_Simple_Potential(thisNode,radius,componentType,massType,haloLoaded)
    !% Compute the gravitational potential due to a black hole.
    use Galacticus_Nodes
    use Numerical_Constants_Physical 
    use Galactic_Structure_Options
    use Black_Hole_Fundamentals
    implicit none
    type            (treeNode              ), intent(inout), pointer  :: thisNode
    integer                                 , intent(in   )           :: componentType,massType
    double precision                        , intent(in   )           :: radius
    logical                                 , intent(in   ), optional :: haloLoaded
    class           (nodeComponentBlackHole),                pointer  :: thisBlackHoleComponent
    double precision                                                  :: componentMass

    ! Set to zero by default.
    Node_Component_Black_Hole_Simple_Potential=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeBlackHole)) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeBlackHole     )) return
    ! Get the black hole component and check that it is of the simple class.
    thisBlackHoleComponent => thisNode%blackHole()
    select type (thisBlackHoleComponent)
    class is (nodeComponentBlackHoleSimple)
       if (Black_Hole_Gravitational_Radius(thisBlackHoleComponent) <= 0.0d0) return
       ! Compute the potential - limit the radius to the gravitational radius to avoid divergent potentials.
       componentMass=Node_Component_Black_Hole_Simple_Enclosed_Mass(thisNode,radius,massType,componentType,weightByMass,weightIndexNull&
            &,haloLoaded)
       Node_Component_Black_Hole_Simple_Potential=-gravitationalConstantGalacticus*componentMass/max(radius &
            &,Black_Hole_Gravitational_Radius(thisBlackHoleComponent))
    end select
    return
  end function Node_Component_Black_Hole_Simple_Potential
  
end module Node_Component_Black_Hole_Simple
