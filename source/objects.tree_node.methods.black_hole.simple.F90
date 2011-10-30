!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module of black hole tree node methods.

module Tree_Node_Methods_Black_Hole_Simple
  !% Implement black hole tree node methods.
  use Tree_Nodes
  use Components
  implicit none
  private
  public :: Tree_Node_Methods_Black_Hole_Simple_Initialize, Galacticus_Output_Tree_Black_Hole_Simple,&
       & Galacticus_Output_Tree_Black_Hole_Simple_Property_Count, Galacticus_Output_Tree_Black_Hole_Simple_Names,&
       & Black_Hole_Simple_Satellite_Merging, Tree_Node_Methods_Black_Hole_Simple_Dump, Black_Hole_Simple_Scale_Set&
       &,Black_Hole_Enclosed_Mass_Simple,Black_Hole_Rotation_Curve_Simple, Black_Hole_Potential_Simple,&
       & Black_Hole_Rotation_Curve_Gradient_Simple
  
  ! The index used as a reference for this component.
  integer :: componentIndex=-1

  ! Property indices.
  integer, parameter :: propertyCount=1, dataCount=0, historyCount=0
  integer, parameter :: massIndex=1

  ! Define procedure pointers.
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Black_Hole_Mass</methodName>
  !# </treeNodeMethodsPointer>

  ! Flag to indicate if this method is selected.
  logical          :: methodSelected=.false.

  ! Seed mass for black holes.
  double precision :: blackHoleSeedMass

  ! Feedback parameters.
  double precision :: blackHoleToSpheroidStellarGrowthRatio,blackHoleWindEfficiency,blackHoleHeatingEfficiency
  logical          :: blackHoleHeatsHotHalo

  ! Output options.
  logical          :: blackHoleOutputAccretion

contains


  !# <treeNodeCreateInitialize>
  !#  <unitName>Tree_Node_Methods_Black_Hole_Simple_Initialize</unitName>
  !#  <optionName>treeNodeMethodBlackHole</optionName>
  !# </treeNodeCreateInitialize>
  subroutine Tree_Node_Methods_Black_Hole_Simple_Initialize(componentOption,componentTypeCount)
    !% Initializes the tree node simple black hole methods module.
    use ISO_Varying_String
    use Input_Parameters
    use String_Handling
    use Galacticus_Display
    implicit none
    type(varying_string), intent(in)    :: componentOption
    integer,              intent(inout) :: componentTypeCount
    type(varying_string)                :: message

    ! Check if this implementation is selected.
    if (componentOption == 'simple') then
       ! Record that method is selected.
       methodSelected=.true.

       ! Increment the component count and store the value for later reference.
       componentTypeCount=componentTypeCount+1
       componentIndex=componentTypeCount

       ! Display message.
       message='Simple black hole method selected [component index '
       message=message//componentIndex//']'
       call Galacticus_Display_Message(message,verbosityInfo)

       ! Set up procedure pointers.
       Tree_Node_Black_Hole_Mass              => Tree_Node_Black_Hole_Mass_Simple
       Tree_Node_Black_Hole_Mass_Set          => Tree_Node_Black_Hole_Mass_Set_Simple
       Tree_Node_Black_Hole_Mass_Rate_Adjust  => Tree_Node_Black_Hole_Mass_Rate_Adjust_Simple
       Tree_Node_Black_Hole_Mass_Rate_Compute => Tree_Node_Black_Hole_Mass_Rate_Compute_Simple
 
       ! Get the black hole seed mass.
       !@ <inputParameter>
       !@   <name>blackHoleSeedMass</name>
       !@   <defaultValue>100</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The mass of the seed black hole placed at the center of each newly formed galaxy.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>blackHoles</group>
       !@ </inputParameter>
       call Get_Input_Parameter("blackHoleSeedMass",blackHoleSeedMass,defaultValue=100.0d0)

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
       else
          blackHoleHeatingEfficiency=0.0d0
       end if

       ! Get options controlling output.
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

    end if
    return
  end subroutine Tree_Node_Methods_Black_Hole_Simple_Initialize

  double precision function Tree_Node_Black_Hole_Mass_Simple(thisNode,instance)
    !% Return the node black hole mass.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Black_Hole_Simple_Index(thisNode)
       Tree_Node_Black_Hole_Mass_Simple=thisNode%components(thisIndex)%instance(1)%properties(massIndex,propertyValue)
    else
       Tree_Node_Black_Hole_Mass_Simple=blackHoleSeedMass
    end if
    return
  end function Tree_Node_Black_Hole_Mass_Simple

  subroutine Tree_Node_Black_Hole_Mass_Set_Simple(thisNode,mass,instance)
    !% Set the node black hole mass.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: mass
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Black_Hole_Simple_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(massIndex,propertyValue)=mass
    return
  end subroutine Tree_Node_Black_Hole_Mass_Set_Simple

  subroutine Tree_Node_Black_Hole_Mass_Rate_Adjust_Simple(thisNode,interrupt,interruptProcedure,rateAdjustment,instance)
    !% Return the node black hole mass rate of change.
    use Cosmological_Parameters
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(),      pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment
    integer                                  :: thisIndex
    
    thisIndex=Tree_Node_Black_Hole_Simple_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(massIndex,propertyDerivative)=thisNode%components(thisIndex)%instance(1)%properties(massIndex&
         &,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Black_Hole_Mass_Rate_Adjust_Simple

  subroutine Tree_Node_Black_Hole_Mass_Rate_Compute_Simple(thisNode,interrupt,interruptProcedure)
    !% Compute the black hole node mass rate of change.
    use Accretion_Disks
    use Cooling_Radii
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Physical
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Math
    use Numerical_Constants_Astronomical
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(),      pointer, intent(inout) :: interruptProcedure
    procedure(),      pointer                :: interruptProcedurePassed
    double precision, parameter              :: coolingRadiusFractionalTransitionMinimum=0.9d0
    double precision, parameter              :: coolingRadiusFractionalTransitionMaximum=1.0d0
    double precision                         :: restMassAccretionRate,massAccretionRate,energyInputRate,heatingRate&
         &,couplingEfficiency,coolingRadiusFractional,x

    ! Get a local copy of the interrupt procedure.
    interruptProcedurePassed => interruptProcedure

    ! Find the rate of rest mass accretion onto the black hole.
    restMassAccretionRate=blackHoleToSpheroidStellarGrowthRatio*Tree_Node_Spheroid_SFR(thisNode)

    ! Finish if there is no accretion.
    if (restMassAccretionRate <= 0.0d0) return

    ! Find the rate of increase in mass of the black hole.
    massAccretionRate=restMassAccretionRate*(1.0d0-blackHoleHeatingEfficiency-blackHoleWindEfficiency)

    ! If no black hole component currently exists and we have some accretion then interrupt and create a black hole.
    if (.not.thisNode%componentExists(componentIndex)) then    
       if (massAccretionRate /= 0.0d0) then
          interrupt=.true.
          interruptProcedure => Black_Hole_Simple_Create
       end if
       return
    end if
    call Tree_Node_Black_Hole_Mass_Rate_Adjust_Simple(thisNode,interrupt,interruptProcedurePassed, massAccretionRate    )

    ! Remove the accreted mass from the spheroid component.
    call Tree_Node_Spheroid_Gas_Sink                 (thisNode,interrupt,interruptProcedurePassed,-restMassAccretionRate)

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
       call Tree_Node_Hot_Halo_Heat_Input(thisNode,interrupt,interruptProcedurePassed,heatingRate)

    end if

    ! Add energy to the spheroid component.
    if (blackHoleWindEfficiency > 0.0d0) then      
       ! Compute the energy input and send it down the spheroid gas energy input pipe.
       energyInputRate=blackHoleWindEfficiency*restMassAccretionRate*(speedLight/kilo)**2
       call Tree_Node_Spheroid_Gas_Energy_Input(thisNode,interrupt,interruptProcedurePassed,energyInputRate)
    end if

    ! Return our local copy of the interrupt procedure.
    interruptProcedure => interruptProcedurePassed

    return
  end subroutine Tree_Node_Black_Hole_Mass_Rate_Compute_Simple

  !# <scaleSetTask>
  !#  <unitName>Black_Hole_Simple_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Black_Hole_Simple_Scale_Set(thisNode)
    !% Set scales for properties of {\tt thisNode}.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex
 
    ! Determine if method is active and a black hole component exists.
    if (methodSelected.and.thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Black_Hole_Simple_Index(thisNode)
       ! Set scale for mass.
       thisNode%components(thisIndex)%instance(1)%properties(massIndex,propertyScale)=max(Tree_Node_Spheroid_Stellar_Mass(thisNode)&
            &*blackHoleToSpheroidStellarGrowthRatio,Tree_Node_Black_Hole_Mass(thisNode))
    end if
    return
  end subroutine Black_Hole_Simple_Scale_Set

  !# <satelliteMergerTask>
  !#  <unitName>Black_Hole_Simple_Satellite_Merging</unitName>
  !# </satelliteMergerTask>
  subroutine Black_Hole_Simple_Satellite_Merging(thisNode)
    !% Merge (instantaneously) any black hole associated with {\tt thisNode} before it merges with its host halo.
    use Black_Hole_Binary_Mergers
    implicit none
    type(treeNode),   pointer, intent(inout)     :: thisNode
    type(treeNode),   pointer                    :: hostNode
    double precision                             :: blackHoleMassNew,blackHoleSpinNew

    if (methodSelected.and.thisNode%componentExists(componentIndex)) then
       
       ! Find the node to merge with.
       call thisNode%mergesWith(hostNode)

       call Black_Hole_Binary_Merger(Tree_Node_Black_Hole_Mass_Simple(thisNode), &
            &                        Tree_Node_Black_Hole_Mass_Simple(hostNode), &
            &                        0.0d0                                     , &
            &                        0.0d0                                     , &
            &                        blackHoleMassNew                          , &
            &                        blackHoleSpinNew                            &
            &                       )

       ! Move the black hole to the host.
       call Tree_Node_Black_Hole_Mass_Set_Simple(hostNode,blackHoleMassNew)
       call Tree_Node_Black_Hole_Mass_Set_Simple(thisNode,0.0d0           )
    end if
    return
  end subroutine Black_Hole_Simple_Satellite_Merging
  
  integer function Tree_Node_Black_Hole_Simple_Index(thisNode)
    !% Ensure the black hole component exists and return its position in the components array.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    
    call thisNode%createComponent(componentIndex,propertyCount,dataCount,historyCount)
    Tree_Node_Black_Hole_Simple_Index=thisNode%componentIndex(componentIndex)
    return
  end function Tree_Node_Black_Hole_Simple_Index

  subroutine Black_Hole_Simple_Create(thisNode)
    !% Creates a black hole component for {\tt thisNode}.
    use ISO_Varying_String
    use Galacticus_Display
    use String_Handling
    implicit none
    type(treeNode),      pointer, intent(inout) :: thisNode
    type(varying_string)                        :: message

    ! Display a message.
    message='Creating black hole component for node '
    message=message//thisNode%index()
    call Galacticus_Display_Message(message,verbosityInfo)
    ! Create the component.
    call thisNode%createComponent(componentIndex,propertyCount,dataCount,historyCount)
    ! Set to the seed mass.
    call Tree_Node_Black_Hole_Mass_Set_Simple(thisNode,blackHoleSeedMass)
    return
  end subroutine Black_Hole_Simple_Create
  
  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Black_Hole_Simple_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Black_Hole_Simple</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Black_Hole_Simple_Names(integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI&
       &,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set names of black hole properties to be written to the \glc\ output file.
    use Numerical_Constants_Astronomical
    use ISO_Varying_String
    implicit none
    double precision, intent(in)                  :: time
    integer,          intent(inout)               :: integerProperty,doubleProperty
    character(len=*), intent(inout), dimension(:) :: integerPropertyNames,integerPropertyComments,doublePropertyNames &
         &,doublePropertyComments
    double precision, intent(inout), dimension(:) :: integerPropertyUnitsSI,doublePropertyUnitsSI
    
    if (methodSelected) then
       !@ <outputPropertyGroup>
       !@   <name>blackHole</name>
       !@   <description>Black hole properities</description>
       !@   <outputType>nodeData</outputType>
       !@ </outputPropertyGroup>
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>blackHoleMass</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Mass of the black hole.</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='blackHoleMass'
       doublePropertyComments(doubleProperty)='Mass of the black hole.'
       doublePropertyUnitsSI (doubleProperty)=massSolar
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
  end subroutine Galacticus_Output_Tree_Black_Hole_Simple_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Black_Hole_Simple_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Black_Hole_Simple</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Black_Hole_Simple_Property_Count(integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of black hole properties to be written to the the \glc\ output file.
    implicit none
    double precision, intent(in)    :: time
    integer,          intent(inout) :: integerPropertyCount,doublePropertyCount
    integer,          parameter     :: extraPropertyCount=1

    if (methodSelected) then
       doublePropertyCount=doublePropertyCount+propertyCount
       if (blackHoleOutputAccretion) doublePropertyCount=doublePropertyCount+extraPropertyCount
    end if
    return
  end subroutine Galacticus_Output_Tree_Black_Hole_Simple_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Black_Hole_Simple</unitName>
  !#  <sortName>Galacticus_Output_Tree_Black_Hole_Simple</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Black_Hole_Simple(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store black hole properties in the \glc\ output file buffers.
    use Tree_Nodes
    use Kind_Numbers
    use Accretion_Disks
    implicit none
    double precision,        intent(in)             :: time
    type(treeNode),          intent(inout), pointer :: thisNode
    integer,                 intent(inout)          :: integerProperty,integerBufferCount,doubleProperty,doubleBufferCount
    integer(kind=kind_int8), intent(inout)          :: integerBuffer(:,:)
    double precision,        intent(inout)          :: doubleBuffer(:,:)
    double precision                                :: restMassAccretionRate

    if (methodSelected) then
       ! Store the properties.
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Black_Hole_Mass_Simple(thisNode)
       if (blackHoleOutputAccretion) then
          ! Get the rest mass accretion rate.
          restMassAccretionRate=blackHoleToSpheroidStellarGrowthRatio*Tree_Node_Spheroid_SFR(thisNode)
          doubleProperty=doubleProperty+1
          doubleBuffer(doubleBufferCount,doubleProperty)=restMassAccretionRate
       end if
    end if
    return
  end subroutine Galacticus_Output_Tree_Black_Hole_Simple

  !# <nodeDumpTask>
  !#  <unitName>Tree_Node_Methods_Black_Hole_Simple_Dump</unitName>
  !# </nodeDumpTask>
  subroutine Tree_Node_Methods_Black_Hole_Simple_Dump(thisNode)
    !% Dump all properties of {\tt thisNode} to screen.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    if (methodSelected) then
       if (thisNode%componentExists(componentIndex)) then
          write (0,'(1x,a)'           ) 'black hole component -> properties:'
          write (0,'(2x,a50,1x,e12.6)') 'black hole mass:',Tree_Node_Black_Hole_Mass(thisNode)
       else
          write (0,'(1x,a)'           ) 'black hole component -> nonexistant'
       end if
    end if
    return
  end subroutine Tree_Node_Methods_Black_Hole_Simple_Dump

  !# <enclosedMassTask>
  !#  <unitName>Black_Hole_Enclosed_Mass_Simple</unitName>
  !# </enclosedMassTask>
  subroutine Black_Hole_Enclosed_Mass_Simple(thisNode,radius,massType,componentType,weightBy,weightIndex,componentMass)
    !% Computes the mass within a given radius for a central black hole. Black hole is treated as a point mass.
    use Galactic_Structure_Options
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    integer,          intent(in)             :: massType,componentType,weightBy,weightIndex
    double precision, intent(in)             :: radius
    double precision, intent(out)            :: componentMass

    ! Set zero enclosed mass by default.
    componentMass=0.0d0

    ! Return the black hole mass only if massType and componentType are of black hole type.
    if (.not.methodSelected                                                                  ) return
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeBlackHole)) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeBlackHole     )) return
    if (.not.thisNode%componentExists(componentIndex)                                        ) return

    ! Compute the enclosed mass.
    if (radius >= 0.0d0) then
       select case (weightBy)
       case (weightByMass)
          componentMass=Tree_Node_Black_Hole_Mass(thisNode,instance=1)
       end select
    end if
    return
  end subroutine Black_Hole_Enclosed_Mass_Simple

  !# <rotationCurveTask>
  !#  <unitName>Black_Hole_Rotation_Curve_Simple</unitName>
  !# </rotationCurveTask>
  subroutine Black_Hole_Rotation_Curve_Simple(thisNode,radius,massType,componentType,componentVelocity)
    !% Computes the rotation curve for the central black hole. Assumes a point mass black hole with a Keplerian rotation curve,
    !% \emph{except} that the rotation speed is limited to never exceed the speed of light.
    use Tree_Nodes
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    use Numerical_Constants_Prefixes
    use Black_Hole_Fundamentals
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    integer,          intent(in)             :: massType,componentType
    double precision, intent(in)             :: radius
    double precision, intent(out)            :: componentVelocity
    double precision                         :: componentMass

    ! Set to zero by default.
    componentVelocity=0.0d0

    ! Compute if a black hole is present.
    if (methodSelected .and. thisNode%componentExists(componentIndex)) then       
       ! Check if the radius exceeds the gravitational radius. Do the calculation here rather than calling the gravitational
       ! radius function since we want to enforce the calculation to always use the first black hole instance. <v0.9.1> This ugly
       ! solution should be solved by passing a black hole object directly to the gravitational radius function.
       if (radius > gravitationalConstantGalacticus*Tree_Node_Black_Hole_Mass(thisNode,instance=1)/(milli*speedLight)**2) then
          ! Radius is larger than the gravitational radius - compute the rotation speed.
          call Black_Hole_Enclosed_Mass_Simple(thisNode,radius,massType,componentType,weightByMass,weightIndexNull,componentMass)
          if (componentMass > 0.0d0) componentVelocity=dsqrt(gravitationalConstantGalacticus*componentMass/radius)
       else
          ! Radius is less than the gravitational radius - return the speed of light.
          componentVelocity=speedLight*milli
       end if
    end if
    return
  end subroutine Black_Hole_Rotation_Curve_Simple

  !# <rotationCurveGradientTask>
  !#  <unitName>Black_Hole_Rotation_Curve_Gradient_Simple</unitName>
  !# </rotationCurveGradientTask>
  subroutine Black_Hole_Rotation_Curve_Gradient_Simple(thisNode,radius,massType,componentType,componentRotationCurveGradient)
    !% Computes the rotation curve gradient for the central black hole. Assumes a point mass black hole with a Keplerian 
    !% rotation curve, \emph{except} that the rotation speed is limited to never exceed the speed of light.
    use Tree_Nodes
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    use Numerical_Constants_Prefixes
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    integer,          intent(in)             :: massType,componentType
    double precision, intent(in)             :: radius
    double precision, intent(out)            :: componentRotationCurveGradient
    double precision                         :: componentMass

    ! Set to zero by default.
    componentRotationCurveGradient=0.0d0
    if (.not.methodSelected                                                                  ) return
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeBlackHole)) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeBlackHole     )) return
    if (.not.thisNode%componentExists(componentIndex)                                        ) return
    if (radius <= 0.0d0) return
    call Black_Hole_Enclosed_Mass_Simple(thisNode,radius,massType,componentType,weightByMass,weightIndexNull,componentMass)
    if (componentMass ==0.0d0 ) return
    if (radius > gravitationalConstantGalacticus*componentMass/(milli*speedLight)**2) then
       componentRotationCurveGradient=-gravitationalConstantGalacticus &
            &                         *componentMass                   &
            &                         /radius**2
    else
       componentRotationCurveGradient=0.0d0
    end if
    return
  end subroutine Black_Hole_Rotation_Curve_Gradient_Simple

  !# <potentialTask>
  !#  <unitName>Black_Hole_Potential_Simple</unitName>
  !# </potentialTask>
  subroutine Black_Hole_Potential_Simple(thisNode,radius,componentType,massType,componentPotential)
    !% Compute the gravitational potential due to a black hole.
    use Tree_Nodes
    use Numerical_Constants_Physical 
    use Galactic_Structure_Options
    use Black_Hole_Fundamentals
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    integer,          intent(in)             :: componentType,massType
    double precision, intent(in)             :: radius
    double precision, intent(out)            :: componentPotential
    double precision                         :: componentMass

    componentPotential=0.0d0
    if (.not.methodSelected                                                                  ) return
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeBlackHole)) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeBlackHole     )) return
    if (.not.thisNode%componentExists(componentIndex)                                        ) return
    if (Black_Hole_Gravitational_Radius(thisNode) <=0.0d0) return
    ! Computes the potential - limit the radius to the gravitational radius to avoid divergent potentials.
    call Black_Hole_Enclosed_Mass_Simple(thisNode,radius,massType,componentType,weightByMass,weightIndexNull,componentMass)
    componentPotential=-gravitationalConstantGalacticus                       &
          &            *componentMass                                         &
          &            /max(radius,Black_Hole_Gravitational_Radius(thisNode))
    return
  end subroutine Black_Hole_Potential_Simple

end module Tree_Node_Methods_Black_Hole_Simple
