!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module of hot halo tree node methods.

module Tree_Node_Methods_Hot_Halo_Very_Simple
  !% Implement hot halo tree node methods.
  use Tree_Nodes
  use Components
  use Tree_Node_Methods_Hot_Halo_Data_Very_Simple
  implicit none
  private
  public :: Tree_Node_Methods_Hot_Halo_Initialize_Very_Simple, Tree_Node_Hot_Halo_Reset_Very_Simple_Very_Simple,&
       & Hot_Halo_Scale_Set_Very_Simple, Hot_Halo_Subresolution_Initialize_Very_Simple,&
       & Hot_Halo_Remove_Before_Satellite_Merging_Very_Simple, Tree_Node_Hot_Halo_Promote_Very_Simple,&
       & Galacticus_Output_Tree_Hot_Halo_Very_Simple_Names, Galacticus_Output_Tree_Hot_Halo_Very_Simple_Property_Count,&
       & Galacticus_Output_Tree_Hot_Halo_Very_Simple, Tree_Node_Methods_Hot_Halo_Very_Simple_Dump,&
       & Hot_Halo_Very_Simple_Property_Identifiers_Decode
  
  ! Property indices.
  integer, parameter :: propertyCount=1, dataCount=0, historyCount=0
  integer, parameter :: massIndex    =1

  ! Define procedure pointers.
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Hot_Halo_Mass</methodName>
  !# </treeNodeMethodsPointer>

  ! Pointer to procedure which will receive gas cooling out of the hot halo.
  !# <treeNodePipePointer>
  !#  <pipeName>Tree_Node_Hot_Halo_Cooling_Mass_To</pipeName>
  !# </treeNodePipePointer>
  !# <treeNodePipePointer>
  !#  <pipeName>Tree_Node_Hot_Halo_Outflow_Mass_To</pipeName>
  !# </treeNodePipePointer>

  ! Quantities stored to avoid repeated computation.
  logical                     :: gotCoolingRate=.false.
  double precision            :: coolingRate
  !$omp threadprivate(gotCoolingRate,coolingRate)

contains

  !# <treeNodeCreateInitialize>
  !#  <unitName>Tree_Node_Methods_Hot_Halo_Initialize_Very_Simple</unitName>
  !#  <optionName>treeNodeMethodHotHalo</optionName>
  !# </treeNodeCreateInitialize>
  subroutine Tree_Node_Methods_Hot_Halo_Initialize_Very_Simple(componentOption,componentTypeCount)
    !% Initializes the tree node hot halo methods module.
    use ISO_Varying_String
    use Input_Parameters
    use String_Handling
    use Galacticus_Display
    use Memory_Management
    use Galacticus_Error
    implicit none
    type(varying_string), intent(in)    :: componentOption
    integer,              intent(inout) :: componentTypeCount
    type(varying_string)                :: message

    ! Check if this implementation is selected.
    if (componentOption == 'verySimple') then
       ! Record that method is selected.
       methodSelected=.true.

       ! Increment the component count and store the value for later reference.
       componentTypeCount=componentTypeCount+1
       componentIndex=componentTypeCount

       ! Display message.
       message='Very simple hot halo method selected [component index '
       message=message//componentIndex//']'
       call Galacticus_Display_Message(message,verbosityInfo)

       ! Set up procedure pointers.
       ! Hot gas reservoir:
       Tree_Node_Hot_Halo_Mass                              => Tree_Node_Hot_Halo_Mass_Very_Simple
       Tree_Node_Hot_Halo_Mass_Set                          => Tree_Node_Hot_Halo_Mass_Set_Very_Simple
       Tree_Node_Hot_Halo_Mass_Rate_Adjust                  => Tree_Node_Hot_Halo_Mass_Rate_Adjust_Very_Simple
       Tree_Node_Hot_Halo_Mass_Rate_Compute                 => Tree_Node_Hot_Halo_Mass_Rate_Compute_Very_Simple

       ! Set internally connected pipes to our procedures.
       Tree_Node_Hot_Halo_Outflow_Mass_To                   => Tree_Node_Hot_Halo_Mass_Rate_Adjust_Very_Simple
       
    end if
    return
  end subroutine Tree_Node_Methods_Hot_Halo_Initialize_Very_Simple
    
  !# <calculationResetTask>
  !# <unitName>Tree_Node_Hot_Halo_Reset_Very_Simple_Very_Simple</unitName>
  !# </calculationResetTask>
  subroutine Tree_Node_Hot_Halo_Reset_Very_Simple_Very_Simple(thisNode)
    !% Remove memory of stored computed values as we're about to begin computing derivatives anew.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode

    gotCoolingRate=.false.
    return
  end subroutine Tree_Node_Hot_Halo_Reset_Very_Simple_Very_Simple

  double precision function Tree_Node_Hot_Halo_Mass_Very_Simple(thisNode,instance)
    !% Return the node hot halo mass.
    implicit none
    integer,        intent(in   ), optional :: instance
    type(treeNode), intent(inout), pointer  :: thisNode
    integer                                 :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Hot_Halo_Index(thisNode)
       Tree_Node_Hot_Halo_Mass_Very_Simple=thisNode%components(thisIndex)%instance(1)%properties(massIndex,propertyValue)
    else
       Tree_Node_Hot_Halo_Mass_Very_Simple=0.0d0
    end if
    return
  end function Tree_Node_Hot_Halo_Mass_Very_Simple

  subroutine Tree_Node_Hot_Halo_Mass_Set_Very_Simple(thisNode,mass,instance)
    !% Set the node hot halo mass.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: mass
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Hot_Halo_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(massIndex,propertyValue)=mass
    return
  end subroutine Tree_Node_Hot_Halo_Mass_Set_Very_Simple

  subroutine Hot_Halo_Very_Simple_Push_To_Cooling_Pipes(thisNode,interrupt,interruptProcedure,massRate)
    !% Push mass through the cooling pipes (along with appropriate amounts of metals and angular momentum) at the given rate.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(),      pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: massRate
    procedure(),      pointer                :: interruptProcedurePassed

    ! Ignore zero rates.
    if (massRate /= 0.0d0 .and. Tree_Node_Hot_Halo_Mass(thisNode) > 0.0d0) then
       
       ! Get a local copy of the interrupt procedure.
       interruptProcedurePassed => interruptProcedure

       ! Remove mass from the hot component.
       call Tree_Node_Hot_Halo_Mass_Rate_Adjust_Very_Simple(thisNode,interrupt,interruptProcedurePassed,-massRate)
       ! Pipe the mass rate to whatever component claimed it.

       if (associated(Tree_Node_Hot_Halo_Cooling_Mass_To)) call Tree_Node_Hot_Halo_Cooling_Mass_To(thisNode,interrupt &
            &,interruptProcedurePassed,massRate)

       ! Return our local copy of the interrupt procedure.
       interruptProcedure => interruptProcedurePassed

    end if
    return
  end subroutine Hot_Halo_Very_Simple_Push_To_Cooling_Pipes

  subroutine Tree_Node_Hot_Halo_Mass_Rate_Adjust_Very_Simple(thisNode,interrupt,interruptProcedure,rateAdjustment,instance)
    !% Return the node hot halo mass rate of change.
    use Cosmological_Parameters
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(), pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment
    integer                                  :: thisIndex
    
    thisIndex=Tree_Node_Hot_Halo_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(massIndex,propertyDerivative)=thisNode%components(thisIndex)%instance(1)%properties(massIndex&
         &,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Hot_Halo_Mass_Rate_Adjust_Very_Simple

  subroutine Tree_Node_Hot_Halo_Mass_Rate_Compute_Very_Simple(thisNode,interrupt,interruptProcedure)
    !% Compute the hot halo node mass rate of change.
    use Accretion_Halos
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    logical,                 intent(inout) :: interrupt
    procedure(),    pointer, intent(inout) :: interruptProcedure
    procedure(),    pointer                :: interruptProcedurePassed
    double precision                       :: massAccretionRate,failedMassAccretionRate

    ! Next compute the cooling rate in this halo.
    call Get_Cooling_Rate(thisNode)
    ! Pipe the cooling rate to which ever component claimed it.
    call Hot_Halo_Very_Simple_Push_To_Cooling_Pipes(thisNode,interrupt,interruptProcedurePassed,coolingRate)

    ! Return a copy of our local interrupt pointer.
    interruptProcedure => interruptProcedurePassed

    return
  end subroutine Tree_Node_Hot_Halo_Mass_Rate_Compute_Very_Simple

  !# <scaleSetTask>
  !#  <unitName>Hot_Halo_Scale_Set_Very_Simple</unitName>
  !# </scaleSetTask>
  subroutine Hot_Halo_Scale_Set_Very_Simple(thisNode)
    !% Set scales for properties of {\tt thisNode}.
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision, parameter              :: scaleMassRelative=1.0d-3
    integer                                  :: thisIndex
    double precision                         :: massVirial

    ! Determine if method is active and a hot halo component exists.
    if (methodSelected.and.thisNode%componentExists(componentIndex)) then
       thisIndex =Tree_Node_Hot_Halo_Index(thisNode)
       massVirial=Tree_Node_Mass          (thisNode)
       thisNode%components(thisIndex)%instance(1)%properties(massIndex,propertyScale)=massVirial*scaleMassRelative
    end if
    return
  end subroutine Hot_Halo_Scale_Set_Very_Simple

  !# <mergerTreeInitializeTask>
  !#  <unitName>Hot_Halo_Subresolution_Initialize_Very_Simple</unitName>
  !# </mergerTreeInitializeTask>
  subroutine Hot_Halo_Subresolution_Initialize_Very_Simple(thisNode)
    !% Initialize the contents of the hot halo component.
    use Cosmological_Parameters
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    type(treeNode),   pointer                :: childNode
    double precision                         :: hotHaloMass

    ! If this method is selected and the node has no child then initialize it.
    if (methodSelected) then
       ! Get the mass of hot gas.
       hotHaloMass=Tree_Node_Mass(thisNode)
       childNode => thisNode%childNode
       do while (associated(childNode))
          hotHaloMass=hotHaloMass-Tree_Node_Mass(childNode)
          childNode => childNode%siblingNode
       end do
       hotHaloMass=hotHaloMass*Omega_b()/Omega_Matter()
       ! If this is non-zero, then create a hot halo component and add this mass to it.
       if (hotHaloMass > 0.0d0) then
          call Hot_Halo_Create                        (thisNode            )
          call Tree_Node_Hot_Halo_Mass_Set_Very_Simple(thisNode,hotHaloMass)
       end if
    end if
    return
  end subroutine Hot_Halo_Subresolution_Initialize_Very_Simple

  !# <satelliteMergerTask>
  !#  <unitName>Hot_Halo_Remove_Before_Satellite_Merging_Very_Simple</unitName>
  !# </satelliteMergerTask>
  subroutine Hot_Halo_Remove_Before_Satellite_Merging_Very_Simple(thisNode)
    !% Remove any hot halo associated with {\tt thisNode} before it merges with its host halo.
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode),   pointer, intent(inout)     :: thisNode
    type(treeNode),   pointer                    :: hostNode

    ! Determine if starvation is to be applied.
    if (methodSelected.and.thisNode%componentExists(componentIndex)) then
       
       ! Find the node to merge with.
       call thisNode%mergesWith(hostNode)
       
       ! Move the hot halo to the host.
       call Tree_Node_Hot_Halo_Mass_Set_Very_Simple(hostNode,Tree_Node_Hot_Halo_Mass_Very_Simple(hostNode) &
            &+Tree_Node_Hot_Halo_Mass_Very_Simple(thisNode))     
    end if
    return
  end subroutine Hot_Halo_Remove_Before_Satellite_Merging_Very_Simple

  !# <nodePromotionTask>
  !#  <unitName>Tree_Node_Hot_Halo_Promote_Very_Simple</unitName>
  !# </nodePromotionTask>
  subroutine Tree_Node_Hot_Halo_Promote_Very_Simple(thisNode)
    !% Ensure that {\tt thisNode} is ready for promotion to its parent. In this case, we simply update the hot halo mass of {\tt
    !% thisNode} to account for any hot halo already in the parent.
    use Galacticus_Error
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode),   pointer, intent(inout)     :: thisNode
    type(treeNode),   pointer                    :: parentNode
    double precision                             :: hotHaloMass

    ! Check if this method is selected.
    if (methodSelected) then
       ! Get the parent node of this node.
       parentNode => thisNode%parentNode
       ! If the parent node has a hot halo component, then add it to that of this node, and perform other changes needed prior to
       ! promotion.
       if (parentNode%componentExists(componentIndex)) then
          hotHaloMass=Tree_Node_Hot_Halo_Mass_Very_Simple(thisNode)+Tree_Node_Hot_Halo_Mass_Very_Simple(parentNode)
          call Tree_Node_Hot_Halo_Mass_Set_Very_Simple(thisNode,hotHaloMass)
       end if
    end if
    return
  end subroutine Tree_Node_Hot_Halo_Promote_Very_Simple

  subroutine Get_Cooling_Rate(thisNode)
    !% Get and store the cooling rate for {\tt thisNode}.
    use Cooling_Rates
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode

    if (.not.gotCoolingRate) then
       if (Tree_Node_Hot_Halo_Mass_Very_Simple(thisNode) > 0.0d0) then
          ! Get the cooling time.
          coolingRate=Cooling_Rate(thisNode)
       else
          coolingRate=0.0d0
       end if

       ! Flag that cooling rate has now been computed.
       gotCoolingRate=.true.
    end if
    return
  end subroutine Get_Cooling_Rate

  integer function Tree_Node_Hot_Halo_Index(thisNode)
    !% Ensure the hot halo component exists and return its position in the components array.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    
    call thisNode%createComponent(componentIndex,propertyCount,dataCount,historyCount)
    Tree_Node_Hot_Halo_Index=thisNode%componentIndex(componentIndex)
    return
  end function Tree_Node_Hot_Halo_Index

  subroutine Hot_Halo_Create(thisNode)
    !% Creates a hot halo component for {\tt thisNode}.
    use ISO_Varying_String
    use Galacticus_Display
    use String_Handling
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode),      pointer, intent(inout) :: thisNode
    type(varying_string)                        :: message

    ! Display a message.
    if (Galacticus_Verbosity_Level() >= verbosityInfo) then
       message='Creating hot halo component for node '
       message=message//thisNode%index()
       call Galacticus_Display_Message(message,verbosityInfo)
    end if
    ! Create the component.
    call thisNode%createComponent(componentIndex,propertyCount,dataCount,historyCount)
    return
  end subroutine Hot_Halo_Create

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Hot_Halo_Very_Simple_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Hot_Halo_Very_Simple</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Hot_Halo_Very_Simple_Names(integerProperty,integerPropertyNames,integerPropertyComments&
       &,integerPropertyUnitsSI,doubleProperty ,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set names of hot halo properties to be written to the \glc\ output file.
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    use ISO_Varying_String
    implicit none
    double precision, intent(in)                  :: time
    integer,          intent(inout)               :: integerProperty,doubleProperty
    character(len=*), intent(inout), dimension(:) :: integerPropertyNames,integerPropertyComments,doublePropertyNames &
         &,doublePropertyComments
    integer                                       :: iAbundance
    double precision, intent(inout), dimension(:) :: integerPropertyUnitsSI,doublePropertyUnitsSI

    if (methodSelected) then
       !@ <outputPropertyGroup>
       !@   <name>hotHalo</name>
       !@   <description>Hot halo properities</description>
       !@   <outputType>nodeData</outputType>
       !@ </outputPropertyGroup>
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>hotHaloMass</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Mass of gas in the hot halo.</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@   <group>hotHalo</group>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='hotHaloMass'
       doublePropertyComments(doubleProperty)='Mass of gas in the hot halo.'
       doublePropertyUnitsSI (doubleProperty)=massSolar
    end if
    return
  end subroutine Galacticus_Output_Tree_Hot_Halo_Very_Simple_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Hot_Halo_Very_Simple_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Hot_Halo_Very_Simple</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Hot_Halo_Very_Simple_Property_Count(integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of hot halo properties to be written to the the \glc\ output file.
    implicit none
    double precision, intent(in)    :: time
    integer,          intent(inout) :: integerPropertyCount,doublePropertyCount

    if (methodSelected) doublePropertyCount=doublePropertyCount+propertyCount
    return
  end subroutine Galacticus_Output_Tree_Hot_Halo_Very_Simple_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Hot_Halo_Very_Simple</unitName>
  !#  <sortName>Galacticus_Output_Tree_Hot_Halo_Very_Simple</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Hot_Halo_Very_Simple(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store hot halo properties in the \glc\ output file buffers.
    use Tree_Nodes
    use Cooling_Radii
    use Kind_Numbers
    implicit none
    double precision,        intent(in)                 :: time
    type(treeNode),          intent(inout), pointer     :: thisNode
    integer,                 intent(inout)              :: integerProperty,integerBufferCount,doubleProperty,doubleBufferCount
    integer(kind=kind_int8), intent(inout)              :: integerBuffer(:,:)
    double precision,        intent(inout)              :: doubleBuffer(:,:)

    if (methodSelected) then
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Hot_Halo_Mass             (thisNode)
    end if
    return
  end subroutine Galacticus_Output_Tree_Hot_Halo_Very_Simple

  !# <nodeDumpTask>
  !#  <unitName>Tree_Node_Methods_Hot_Halo_Very_Simple_Dump</unitName>
  !# </nodeDumpTask>
  subroutine Tree_Node_Methods_Hot_Halo_Very_Simple_Dump(thisNode)
    !% Dump all properties of {\tt thisNode} to screen.
    use ISO_Varying_String
    implicit none
    type(treeNode),   intent(inout), pointer     :: thisNode

    if (methodSelected) then
       if (thisNode%componentExists(componentIndex)) then
          write (0,'(1x,a)'           ) 'hot halo component -> properties:'
          write (0,'(2x,a50,1x,e12.6)') 'hot halo mass:',Tree_Node_Hot_Halo_Mass_Very_Simple(thisNode)
       else
          write (0,'(1x,a)'           ) 'hot halo component -> nonexistant'
       end if
    end if
    return
  end subroutine Tree_Node_Methods_Hot_Halo_Very_Simple_Dump

  !# <decodePropertyIdentifiersTask>
  !#  <unitName>Hot_Halo_Very_Simple_Property_Identifiers_Decode</unitName>
  !# </decodePropertyIdentifiersTask>
  subroutine Hot_Halo_Very_Simple_Property_Identifiers_Decode(propertyComponent,propertyObject,propertyIndex,matchedProperty,propertyName)
    !% Decodes property identifiers to property names for the standard hot halo module.
    use ISO_Varying_String
    implicit none
    integer,              intent(in)    :: propertyComponent,propertyObject,propertyIndex
    logical,              intent(inout) :: matchedProperty
    type(varying_string), intent(inout) :: propertyName

    if (methodSelected.and..not.matchedProperty) then
       if (propertyComponent == componentIndex) then
          matchedProperty=.true.
          propertyName="hotHalo:"
          select case (propertyObject)
          case (objectTypeProperty)
             if (propertyIndex == massIndex) then
                propertyName=propertyName//":hotGasMass"
             end if
          end select
       end if
    end if

    return
  end subroutine Hot_Halo_Very_Simple_Property_Identifiers_Decode
  
end module Tree_Node_Methods_Hot_Halo_Very_Simple
