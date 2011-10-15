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


!% Contains a module of preset spin tree node methods.

module Tree_Node_Methods_Spin_Preset
  !% Implement preset spin tree node method.
  use Tree_Nodes
  use Components
  implicit none
  private
  public :: Tree_Node_Methods_Spin_Preset_Initialize, Tree_Node_Methods_Spin_Preset_Initialize_Rates,&
       & Galacticus_Output_Tree_Spin_Preset, Galacticus_Output_Tree_Spin_Preset_Property_Count,&
       & Galacticus_Output_Tree_Spin_Preset_Names, Tree_Node_Methods_Spin_Preset_Dump, Tree_Node_Spin_Preset_Promote
  
  ! The index used as a reference for this component.
  integer :: componentIndex=-1

  ! Property indices.
  integer, parameter :: propertyCount=0, dataCount=2, historyCount=0
  integer, parameter :: spinIndex          =1
  integer, parameter :: spinGrowthRateIndex=2

  ! Define procedure pointers.
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Spin</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Spin_Growth_Rate</methodName>
  !# </treeNodeMethodsPointer>

  ! Flag to indicate if this method is selected.
  logical :: methodSelected=.false.

contains

  !# <treeNodeCreateInitialize>
  !#  <unitName>Tree_Node_Methods_Spin_Preset_Initialize</unitName>
  !#  <optionName>treeNodeMethodSpin</optionName>
  !# </treeNodeCreateInitialize>
  subroutine Tree_Node_Methods_Spin_Preset_Initialize(componentOption,componentTypeCount)
    !% Initializes the tree node preset spin method module.
    use ISO_Varying_String
    use String_Handling
    use Galacticus_Display
    implicit none
    type(varying_string), intent(in)    :: componentOption
    integer,              intent(inout) :: componentTypeCount
    type(varying_string)                :: message

    ! Check if this implementation is selected.
    if (componentOption == 'preset') then
       ! Record that method is selected.
       methodSelected=.true.

       ! Increment the component count and store the value for later reference.
       componentTypeCount=componentTypeCount+1
       componentIndex=componentTypeCount

       ! Display message.
       message='Preset spin method selected [component index '
       message=message//componentIndex//']'
       call Galacticus_Display_Message(message,verbosityInfo)

       ! Set up procedure pointers.
       Tree_Node_Spin                          => Tree_Node_Spin_Preset
       Tree_Node_Spin_Set                      => Tree_Node_Spin_Set_Preset
       Tree_Node_Spin_Rate_Adjust              => null()
       Tree_Node_Spin_Rate_Compute             => Tree_Node_Rate_Rate_Compute_Dummy
       Tree_Node_Spin_Growth_Rate              => Tree_Node_Spin_Growth_Rate_Preset
       Tree_Node_Spin_Growth_Rate_Set          => null()
       Tree_Node_Spin_Growth_Rate_Rate_Adjust  => null()
       Tree_Node_Spin_Growth_Rate_Rate_Compute => Tree_Node_Rate_Rate_Compute_Dummy
    end if
    return
  end subroutine Tree_Node_Methods_Spin_Preset_Initialize

  subroutine Tree_Node_Spin_Set_Preset(thisNode,spin,instance)
    !% Set the node spin mass.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: spin
    integer,          intent(in), optional   :: instance
    integer                                  :: thisIndex
    
    thisIndex=Tree_Node_Spin_Preset_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%data(spinIndex)=spin
    return
  end subroutine Tree_Node_Spin_Set_Preset

  double precision function Tree_Node_Spin_Preset(thisNode,instance)
    !% Return the node spin.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex
    
    thisIndex=Tree_Node_Spin_Preset_Index(thisNode)
    Tree_Node_Spin_Preset=thisNode%components(thisIndex)%instance(1)%data(spinIndex)
    return
  end function Tree_Node_Spin_Preset

  double precision function Tree_Node_Spin_Growth_Rate_Preset(thisNode,instance)
    !% Return the node spin growth rate.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    thisIndex=Tree_Node_Spin_Preset_Index(thisNode)
    Tree_Node_Spin_Growth_Rate_Preset=thisNode%components(thisIndex)%instance(1)%data(spinGrowthRateIndex)
    return
  end function Tree_Node_Spin_Growth_Rate_Preset

  !# <mergerTreeInitializeTask>
  !#  <unitName>Tree_Node_Methods_Spin_Preset_Initialize_Rates</unitName>
  !# </mergerTreeInitializeTask>
  subroutine Tree_Node_Methods_Spin_Preset_Initialize_Rates(thisNode)
    !% Initialize the spin of {\tt thisNode}.
    use Cosmological_Parameters
    use Halo_Spin_Distributions
    implicit none
    type(treeNode),  pointer, intent(inout) :: thisNode
    integer                                 :: thisIndex
    double precision                        :: deltaTime

    if (methodSelected) then
       thisIndex=Tree_Node_Spin_Preset_Index(thisNode)
       ! Check if this node is the primary progenitor.
       if (thisNode%isPrimaryProgenitor()) then
          ! It is, so compute the spin growth rate.
          deltaTime=Tree_Node_Time(thisNode%parentNode)-Tree_Node_Time(thisNode)
          if (deltaTime > 0.0d0) then
             thisNode%components(thisIndex)%instance(1)%data(spinGrowthRateIndex)=(Tree_Node_Spin_Preset(thisNode%parentNode) &
                  &-Tree_Node_Spin_Preset(thisNode))/deltaTime
          else
             thisNode%components(thisIndex)%instance(1)%data(spinGrowthRateIndex)=0.0d0
          end if
       else
          ! It is not, so set spin growth rate to zero.
          thisNode%components(thisIndex)%instance(1)%data(spinGrowthRateIndex)=0.0d0
       end if
    end if
    return
  end subroutine Tree_Node_Methods_Spin_Preset_Initialize_Rates

  integer function Tree_Node_Spin_Preset_Index(thisNode)
    !% Ensure the spin component exists and return its position in the components array.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    
    call thisNode%createComponent(componentIndex,propertyCount,dataCount,historyCount)
    Tree_Node_Spin_Preset_Index=thisNode%componentIndex(componentIndex)
    return
  end function Tree_Node_Spin_Preset_Index

  !# <nodePromotionTask>
  !#  <unitName>Tree_Node_Spin_Preset_Promote</unitName>
  !# </nodePromotionTask>
  subroutine Tree_Node_Spin_Preset_Promote(thisNode)
    !% Ensure that {\tt thisNode} is ready for promotion to its parent. In this case, we simply update the spin of {\tt thisNode}
    !% to be that of its parent.
    use Galacticus_Error
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    type(treeNode), pointer                :: parentNode
    integer                                :: thisIndex

    if (methodSelected) then
       parentNode => thisNode%parentNode
       if (Tree_Node_Time(thisNode) /= Tree_Node_Time(parentNode)) call Galacticus_Error_Report('Tree_Node_Spin_Preset_Promote','thisNode&
            & has not been evolved to its parent')
       ! Adjust the spin (and growth rate) to that of the parent node.
       thisIndex=Tree_Node_Spin_Preset_Index(thisNode)            
       thisNode%components(thisIndex)%instance(1)%data(spinIndex          )=Tree_Node_Spin_Preset            (parentNode)
       thisNode%components(thisIndex)%instance(1)%data(spinGrowthRateIndex)=Tree_Node_Spin_Growth_Rate_Preset(parentNode)
    end if
    return
  end subroutine Tree_Node_Spin_Preset_Promote

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Spin_Preset_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Spin_Preset</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Spin_Preset_Names(integerProperty,integerPropertyNames,integerPropertyComments&
       &,integerPropertyUnitsSI,doubleProperty ,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set names of spin properties to be written to the \glc\ output file.
    implicit none
    double precision, intent(in)                  :: time
    integer,          intent(inout)               :: integerProperty,doubleProperty
    character(len=*), intent(inout), dimension(:) :: integerPropertyNames,integerPropertyComments,doublePropertyNames &
         &,doublePropertyComments
    double precision, intent(inout), dimension(:) :: integerPropertyUnitsSI,doublePropertyUnitsSI

    if (methodSelected) then
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>nodeSpin</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Spin parameter of the node.</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='nodeSpin'
       doublePropertyComments(doubleProperty)='Spin parameter of the node.'
       doublePropertyUnitsSI (doubleProperty)=0.0d0
    end if
    return
  end subroutine Galacticus_Output_Tree_Spin_Preset_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Spin_Preset_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Spin_Preset</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Spin_Preset_Property_Count(integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of spin properties to be written to the the \glc\ output file.
    implicit none
    double precision, intent(in)    :: time
    integer,          intent(inout) :: integerPropertyCount,doublePropertyCount

    if (methodSelected) doublePropertyCount=doublePropertyCount+1
    return
  end subroutine Galacticus_Output_Tree_Spin_Preset_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Spin_Preset</unitName>
  !#  <sortName>Galacticus_Output_Tree_Spin_Preset</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Spin_Preset(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store spin properties in the \glc\ output file buffers.
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
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Spin_Preset(thisNode)
    end if
    return
  end subroutine Galacticus_Output_Tree_Spin_Preset

  !# <nodeDumpTask>
  !#  <unitName>Tree_Node_Methods_Spin_Preset_Dump</unitName>
  !# </nodeDumpTask>
  subroutine Tree_Node_Methods_Spin_Preset_Dump(thisNode)
    !% Dump all properties of {\tt thisNode} to screen.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode
 
    if (methodSelected) then
       if (thisNode%componentExists(componentIndex)) then
          write (0,'(1x,a)'           ) 'spin component -> properties:'
          write (0,'(2x,a50,1x,e12.6)') '            halo spin:',Tree_Node_Spin_Preset            (thisNode)
          write (0,'(2x,a50,1x,e12.6)') 'halo spin growth rate:',Tree_Node_Spin_Growth_Rate_Preset(thisNode)
       else
          write (0,'(1x,a)'           ) 'spin component -> nonexistant'
       end if
    end if
    return
  end subroutine Tree_Node_Methods_Spin_Preset_Dump

end module Tree_Node_Methods_Spin_Preset
