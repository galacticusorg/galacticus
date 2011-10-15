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


!% Contains a module with the standard implementation of basic tree node methods.

module Tree_Node_Methods_Basic_Standard
  !% The standard implementation of basic tree node methods.
  use Tree_Nodes
  use Components
  implicit none
  private
  public :: Tree_Node_Methods_Basic_Initialize_Standard, Tree_Node_Basic_Promote_Standard, Halo_Mass_Accretion_Rate_Standard, Tree_Node_Mass_Stop_Accretion_Standard,&
       & Galacticus_Output_Tree_Basic_Names_Standard, Galacticus_Output_Tree_Basic_Property_Count_Standard, Galacticus_Output_Tree_Basic_Standard,&
       & Tree_Node_Methods_Basic_Dump_Standard, Basic_Standard_Scale_Set
  
  ! The index used as a reference for this component.
  integer :: componentIndex=-1

  ! Property indices.
  integer, parameter :: propertyCount=2, dataCount=2, historyCount=0
  integer, parameter :: massIndex     =1
  integer, parameter :: timeIndex     =2
  integer, parameter :: rateIndex     =1
  integer, parameter :: isolationIndex=2

  ! Define procedure pointers.
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Mass</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Time</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Mass_Accretion_Rate</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Time_Last_Isolated</methodName>
  !# </treeNodeMethodsPointer>

  ! Flag to indicate if this method is selected.
  logical :: methodSelected=.false.

contains


  !# <treeNodeCreateInitialize>
  !#  <unitName>Tree_Node_Methods_Basic_Initialize_Standard</unitName>
  !#  <optionName default="standard">treeNodeMethodBasic</optionName>
  !# </treeNodeCreateInitialize>
  subroutine Tree_Node_Methods_Basic_Initialize_Standard(componentOption,componentTypeCount)
    !% Initializes the tree node basic methods module.
    use ISO_Varying_String
    use String_Handling
    use Galacticus_Display
    implicit none
    type(varying_string), intent(in)    :: componentOption
    integer,              intent(inout) :: componentTypeCount
    type(varying_string)                :: message

    ! Check if this implementation is selected.
    if (componentOption == 'standard') then
       ! Record that method is selected.
       methodSelected=.true.
       
       ! Increment the component count and store the value for later reference.
       componentTypeCount=componentTypeCount+1
       componentIndex=componentTypeCount

       ! Display message.
       message='Standard basic method selected [component index '
       message=message//componentIndex//']'
       call Galacticus_Display_Message(message,verbosityInfo)

       ! Set up procedure pointers.
       Tree_Node_Mass              => Tree_Node_Mass_Basic
       Tree_Node_Mass_Set          => Tree_Node_Mass_Set_Basic
       Tree_Node_Mass_Rate_Adjust  => Tree_Node_Mass_Rate_Adjust_Basic
       Tree_Node_Mass_Rate_Compute => Tree_Node_Mass_Rate_Compute_Basic
 
       Tree_Node_Mass_Accretion_Rate              => Tree_Node_Mass_Accretion_Rate_Basic
       Tree_Node_Mass_Accretion_Rate_Set          => null()
       Tree_Node_Mass_Accretion_Rate_Rate_Adjust  => null()
       Tree_Node_Mass_Accretion_Rate_Rate_Compute => Tree_Node_Rate_Rate_Compute_Dummy

       Tree_Node_Time              => Tree_Node_Time_Basic
       Tree_Node_Time_Set          => Tree_Node_Time_Set_Basic
       Tree_Node_Time_Rate_Adjust  => Tree_Node_Time_Rate_Adjust_Basic
       Tree_Node_Time_Rate_Compute => Tree_Node_Time_Rate_Compute_Basic

       Tree_Node_Time_Last_Isolated              => Tree_Node_Time_Last_Isolated_Basic
       Tree_Node_Time_Last_Isolated_Set          => null()
       Tree_Node_Time_Last_Isolated_Rate_Adjust  => null()
       Tree_Node_Time_Last_Isolated_Rate_Compute => Tree_Node_Rate_Rate_Compute_Dummy
      end if
    return
  end subroutine Tree_Node_Methods_Basic_Initialize_Standard

  double precision function Tree_Node_Mass_Basic(thisNode,instance)
    !% Return the node mass.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    thisIndex=Tree_Node_Basic_Index(thisNode)
    Tree_Node_Mass_Basic=thisNode%components(thisIndex)%instance(1)%properties(massIndex,propertyValue)
    return
  end function Tree_Node_Mass_Basic

  subroutine Tree_Node_Mass_Set_Basic(thisNode,mass,instance)
    !% Set the node mass.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: mass
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Basic_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(massIndex,propertyValue)=mass
    return
  end subroutine Tree_Node_Mass_Set_Basic

  subroutine Tree_Node_Mass_Rate_Adjust_Basic(thisNode,interrupt,interruptProcedure,rateAdjustment,instance)
    !% Return the node mass rate of change.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(), pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment
    integer                                  :: thisIndex

    ! Use the stored accretion rate value.
    thisIndex=Tree_Node_Basic_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(massIndex,propertyDerivative)=thisNode%components(thisIndex)%instance(1)%properties(massIndex&
         &,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Mass_Rate_Adjust_Basic

  subroutine Tree_Node_Mass_Rate_Compute_Basic(thisNode,interrupt,interruptProcedure)
    !% Compute the node mass rate of change.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    logical,                 intent(inout) :: interrupt
    procedure(), pointer, intent(inout) :: interruptProcedure
    integer                                :: thisIndex

    thisIndex=Tree_Node_Basic_Index(thisNode)
    call Tree_Node_Mass_Rate_Adjust_Basic(thisNode,interrupt,interruptProcedure,thisNode%components(thisIndex)%instance(1)%data(rateIndex))
    return
  end subroutine Tree_Node_Mass_Rate_Compute_Basic

  double precision function Tree_Node_Mass_Accretion_Rate_Basic(thisNode,instance)
    !% Returns the mass accretion rate for {\tt thisNode}.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    thisIndex=Tree_Node_Basic_Index(thisNode)
    Tree_Node_Mass_Accretion_Rate_Basic=thisNode%components(thisIndex)%instance(1)%data(rateIndex)
    return
  end function Tree_Node_Mass_Accretion_Rate_Basic

  double precision function Tree_Node_Time_Basic(thisNode,instance)
    !% Return the node time.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    thisIndex=Tree_Node_Basic_Index(thisNode)
    Tree_Node_Time_Basic=thisNode%components(thisIndex)%instance(1)%properties(timeIndex,propertyValue)
    return
  end function Tree_Node_Time_Basic

  subroutine Tree_Node_Time_Set_Basic(thisNode,time,instance)
    !% Set the node time.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: time
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Basic_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(timeIndex,propertyValue)=time
    return
  end subroutine Tree_Node_Time_Set_Basic

  subroutine Tree_Node_Time_Rate_Adjust_Basic(thisNode,interrupt,interruptProcedure,rateAdjustment,instance)
    !% Return the node time rate of change.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(), pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment
    integer                                  :: thisIndex
    
    thisIndex=Tree_Node_Basic_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(timeIndex,propertyDerivative)=thisNode%components(thisIndex)%instance(1)%properties(timeIndex&
         &,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Time_Rate_Adjust_Basic

  subroutine Tree_Node_Time_Rate_Compute_Basic(thisNode,interrupt,interruptProcedure)
    !% Compute the node time rate of change.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    logical,                 intent(inout) :: interrupt
    procedure(), pointer, intent(inout) :: interruptProcedure
    
    call Tree_Node_Time_Rate_Adjust_Basic(thisNode,interrupt,interruptProcedure,1.0d0) ! By definition.
    return
  end subroutine Tree_Node_Time_Rate_Compute_Basic

  double precision function Tree_Node_Time_Last_Isolated_Basic(thisNode,instance)
    !% Returns the time at which {\tt thisNode} was last an isolated halo.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    if (thisNode%isSatellite()) then
       thisIndex=Tree_Node_Basic_Index(thisNode)
       Tree_Node_Time_Last_Isolated_Basic=thisNode%components(thisIndex)%instance(1)%data(isolationIndex)
    else
       Tree_Node_Time_Last_Isolated_Basic=Tree_Node_Time_Basic(thisNode)
    end if
    return
  end function Tree_Node_Time_Last_Isolated_Basic

  !# <scaleSetTask>
  !#  <unitName>Basic_Standard_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Basic_Standard_Scale_Set(thisNode)
    !% Set scales for properties of {\tt thisNode}.
    implicit none
    type(treeNode),   pointer, intent(inout)       :: thisNode
    double precision, parameter                    :: timeScale        =1.0d-3
    double precision, parameter                    :: scaleMassRelative=1.0d-6
    integer                                        :: thisIndex
 
    ! Determine if method is active and a basic component exists.
    if (methodSelected.and.thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Basic_Index(thisNode)

       ! Set scale for time.
       thisNode%components(thisIndex)%instance(1)%properties(timeIndex,propertyScale)=timeScale

       ! Set scale for mass.
       thisNode%components(thisIndex)%instance(1)%properties(massIndex,propertyScale)=Tree_Node_Mass(thisNode)*scaleMassRelative

    end if
    return
  end subroutine Basic_Standard_Scale_Set

  !# <mergerTreeInitializeTask>
  !#  <unitName>Halo_Mass_Accretion_Rate_Standard</unitName>
  !# </mergerTreeInitializeTask>
  subroutine Halo_Mass_Accretion_Rate_Standard(thisNode)
    !% Set the mass accretion rate for {\tt thisNode}.
    implicit none
    type(treeNode),  pointer, intent(inout) :: thisNode
    integer                                 :: thisIndex
    double precision                        :: deltaTime,massUnresolved,progenitorMassTotal

    if (methodSelected) then
       ! Get component index for this node.
       thisIndex=Tree_Node_Basic_Index(thisNode)

       ! Determine if this node has a descendent.
       if (.not.associated(thisNode%parentNode)) then
          ! For parent-less nodes (i.e. the root node of the tree), the rate is set equal to that of the
          ! progenitor, if it has one.
          if (associated(thisNode%childNode)) then
             ! Ensure the child has a mass growth rate computed.
             call Halo_Mass_Accretion_Rate_Standard(thisNode%childNode)
             ! Get the growth rate of the child.
             thisNode%components(thisIndex)%instance(1)%data(rateIndex)=Tree_Node_Mass_Accretion_Rate_Basic(thisNode%childNode)
          else
             ! Parentless node has no child - set a zero growth rate.
             thisNode%components(thisIndex)%instance(1)%data(rateIndex)=0.0d0
          end if
       else
          ! Compute the unresolved mass.
          massUnresolved=Unresolved_Mass(thisNode%parentNode)
          if (massUnresolved > 0.0d0) then
             ! Positive mass growth - assume this occurs entirely in the main progenitor.
             if (thisNode%isPrimaryProgenitor()) then
                ! Main progenitor - compute required growth rate.
                deltaTime=Tree_Node_Time(thisNode%parentNode)-Tree_Node_Time(thisNode)
                if (deltaTime > 0.0d0) thisNode%components(thisIndex)%instance(1)%data(rateIndex)=massUnresolved/deltaTime
             else
                ! Non-main progenitor - assume zero growth rate.
                thisNode%components(thisIndex)%instance(1)%data(rateIndex)=0.0d0
             end if
          else
             ! Negative mass growth - assume all progenitors lose mass at proportionally equal rates.
             ! Compute the total mass in progenitors.
             progenitorMassTotal=Tree_Node_Mass_Basic(thisNode%parentNode)-massUnresolved
             ! Compute the time available for accretion.
             deltaTime=Tree_Node_Time(thisNode%parentNode)-Tree_Node_Time(thisNode)
             ! Compute mass growth rate.
             if (deltaTime > 0.0d0) thisNode%components(thisIndex)%instance(1)%data(rateIndex)=(massUnresolved/deltaTime)&
                  &*(Tree_Node_Mass_Basic(thisNode)/progenitorMassTotal)
          end if
      end if

    end if
    return
  end subroutine Halo_Mass_Accretion_Rate_Standard

  !# <nodeMergerTask>
  !#  <unitName>Tree_Node_Mass_Stop_Accretion_Standard</unitName>
  !# </nodeMergerTask>
  subroutine Tree_Node_Mass_Stop_Accretion_Standard(thisNode)
    !% Switch off accretion of new mass onto this node once it becomes a satellite.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    if (methodSelected) then
       thisIndex=Tree_Node_Basic_Index(thisNode)
       ! Shut down mass accretion onto the halo now that it is a satellite.
       thisNode%components(thisIndex)%instance(1)%data(rateIndex)=0.0d0
       ! Record the time at which the node became a satellite - used for computing halo scales etc.
       thisNode%components(thisIndex)%instance(1)%data(isolationIndex)=Tree_Node_Time(thisNode)
    end if
    return
  end subroutine Tree_Node_Mass_Stop_Accretion_Standard

  !# <nodePromotionTask>
  !#  <unitName>Tree_Node_Basic_Promote_Standard</unitName>
  !# </nodePromotionTask>
  subroutine Tree_Node_Basic_Promote_Standard(thisNode)
    !% Ensure that {\tt thisNode} is ready for promotion to its parent. In this case, we simply update the mass of {\tt thisNode}
    !% to be that of its parent.
    use Galacticus_Error
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    type(treeNode), pointer                :: parentNode
    integer                                :: thisIndex,parentIndex

    if (methodSelected) then
       parentNode => thisNode%parentNode
       if (Tree_Node_Time(thisNode) /= Tree_Node_Time(parentNode)) call Galacticus_Error_Report('Tree_Node_Basic_Promote','thisNode&
            & has not been evolved to its parent')
       ! Adjust the mass to that of the parent node.
       call Tree_Node_Mass_Set_Basic(thisNode,Tree_Node_Mass_Basic(parentNode))
       ! Adjust the accretion rate to that of the parent node.
       thisIndex  =Tree_Node_Basic_Index(thisNode)
       parentIndex=Tree_Node_Basic_Index(parentNode)
       thisNode%components(thisIndex)%instance(1)%data(rateIndex)=parentNode%components(parentIndex)%instance(1)%data(rateIndex)
    end if
    return
  end subroutine Tree_Node_Basic_Promote_Standard

  integer function Tree_Node_Basic_Index(thisNode)
    !% Ensure the basic component exists and return its position in the components array.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode

    call thisNode%createComponent(componentIndex,propertyCount,dataCount,historyCount)
    Tree_Node_Basic_Index=thisNode%componentIndex(componentIndex)
    return
  end function Tree_Node_Basic_Index

  double precision function Unresolved_Mass(thisNode)
    !% Return the unresolved mass for {\tt thisNode}.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    type(treeNode), pointer                :: childNode

    Unresolved_Mass=Tree_Node_Mass(thisNode)
    childNode => thisNode%childNode
    do while (associated(childNode))
       Unresolved_Mass=Unresolved_Mass-Tree_Node_Mass(childNode)
       childNode => childNode%siblingNode
    end do
    return
  end function Unresolved_Mass

  
  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Basic_Names_Standard</unitName>
  !#  <sortName>Galacticus_Output_Tree_Basic_Standard</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Basic_Names_Standard(integerProperty,integerPropertyNames,integerPropertyComments&
       &,integerPropertyUnitsSI,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set names of basic properties to be written to the \glc\ output file.
    use Numerical_Constants_Astronomical
    implicit none
    double precision, intent(in)                  :: time
    integer,          intent(inout)               :: integerProperty,doubleProperty
    character(len=*), intent(inout), dimension(:) :: integerPropertyNames,integerPropertyComments,doublePropertyNames &
         &,doublePropertyComments
    double precision, intent(inout), dimension(:) :: integerPropertyUnitsSI,doublePropertyUnitsSI

    if (methodSelected) then
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>nodeMass</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Total mass of the node, assuming univeral baryon fraction.</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='nodeMass'
       doublePropertyComments(doubleProperty)='Total mass of the node, assuming univeral baryon fraction.'
       doublePropertyUnitsSI (doubleProperty)=massSolar
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>nodeTimeLastIsolated</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Time at which node was last an isolated halo.</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='nodeTimeLastIsolated'
       doublePropertyComments(doubleProperty)='Time at which node was last an isolated halo.'
       doublePropertyUnitsSI (doubleProperty)=gigaYear
    end if
    return
  end subroutine Galacticus_Output_Tree_Basic_Names_Standard

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Basic_Property_Count_Standard</unitName>
  !#  <sortName>Galacticus_Output_Tree_Basic_Standard</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Basic_Property_Count_Standard(integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of basic properties to be written to the the \glc\ output file.
    implicit none
    double precision, intent(in)    :: time
    integer,          intent(inout) :: integerPropertyCount,doublePropertyCount

    if (methodSelected) doublePropertyCount=doublePropertyCount+2
    return
  end subroutine Galacticus_Output_Tree_Basic_Property_Count_Standard

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Basic_Standard</unitName>
  !#  <sortName>Galacticus_Output_Tree_Basic_Standard</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Basic_Standard(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store basic properties in the \glc\ output file buffers.
    use Tree_Nodes
    use Kind_Numbers
    implicit none
    double precision,        intent(in)             :: time
    type(treeNode),          intent(inout), pointer :: thisNode
    integer,                 intent(inout)          :: integerProperty,integerBufferCount,doubleProperty,doubleBufferCount
    integer(kind=kind_int8), intent(inout)          :: integerBuffer(:,:)
    double precision,        intent(inout)          :: doubleBuffer(:,:)

    if (methodSelected) then
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Mass(thisNode)
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Time_Last_Isolated_Basic(thisNode)
    end if
    return
  end subroutine Galacticus_Output_Tree_Basic_Standard

  !# <nodeDumpTask>
  !#  <unitName>Tree_Node_Methods_Basic_Dump_Standard</unitName>
  !# </nodeDumpTask>
  subroutine Tree_Node_Methods_Basic_Dump_Standard(thisNode)
    !% Dump all properties of {\tt thisNode} to screen.
    implicit none
    type(treeNode),   intent(inout), pointer     :: thisNode

    if (methodSelected) then
       if (thisNode%componentExists(componentIndex)) then
          write (0,'(1x,a)'           ) 'basic node component -> properties:'
          write (0,'(2x,a50,1x,e12.6)') 'node mass:'              ,Tree_Node_Mass(thisNode)
          write (0,'(2x,a50,1x,e12.6)') 'node time:'              ,Tree_Node_Time(thisNode)
          write (0,'(2x,a50,1x,e12.6)') 'node time last isolated:',Tree_Node_Time_Last_Isolated_Basic(thisNode)
       else
          write (0,'(1x,a)'           ) 'basic node component -> nonexistant'
       end if
    end if
    return
  end subroutine Tree_Node_Methods_Basic_Dump_Standard

end module Tree_Node_Methods_Basic_Standard
