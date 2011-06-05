!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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


!% Contains a module of satellite orbit tree node methods.

module Tree_Node_Methods_Satellite_Orbit
  !% Implement satellite orbit tree node methods.
  use Tree_Nodes
  use Tree_Node_Methods
  use Components
  private
  public :: Tree_Node_Methods_Satellite_Orbit_Initialize, Satellite_Orbit_Create, Galacticus_Output_Tree_Satellite_Orbit_Simple,&
       & Galacticus_Output_Tree_Satellite_Orbit_Simple_Property_Count, Galacticus_Output_Tree_Satellite_Orbit_Simple_Names,&
       & Tree_Node_Methods_Satellite_Orbit_Simple_Dump
  
  ! The index used as a reference for this component.
  integer :: componentIndex=-1

  ! Property indices.
  integer, parameter :: propertyCount =1, dataCount=0, historyCount=0
  integer, parameter :: mergeTimeIndex=1

  ! Define procedure pointers.
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Satellite_Merge_Time</methodName>
  !# </treeNodeMethodsPointer>

  ! Procedure pointer for function that will be called to assign merging times to satellites.
  procedure(Satellite_Time_Until_Merging_Template), pointer :: Satellite_Time_Until_Merging => null()
  interface Satellite_Time_Until_Merging_Template
     double precision function Satellite_Time_Until_Merging_Template(thisNode)
       import treeNode
       type(treeNode), pointer, intent(in) :: thisNode
     end function Satellite_Time_Until_Merging_Template
  end interface

  ! Flag to indicate if this method is selected.
  logical :: methodSelected=.false.

contains


  !# <treeNodeCreateInitialize>
  !#  <unitName>Tree_Node_Methods_Satellite_Orbit_Initialize</unitName>
  !#  <optionName default="simple">treeNodeMethodSatelliteOrbit</optionName>
  !# </treeNodeCreateInitialize>
  subroutine Tree_Node_Methods_Satellite_Orbit_Initialize(componentOption,componentTypeCount)
    !% Initializes the tree node satellite orbit methods module.
    use ISO_Varying_String
    use Input_Parameters
    use Galacticus_Error
    use Galacticus_Display
    use String_Handling
    !# <include directive="satelliteMergingMethod" type="moduleUse">
    include 'objects.tree_node.methods.satellite_orbit.moduleUse.inc'
    !# </include>
    implicit none
    type(varying_string), intent(in)    :: componentOption
    integer,              intent(inout) :: componentTypeCount
    type(varying_string)                :: satelliteMergingMethod,message

    ! Check if this implementation is selected.
    if (componentOption.eq.'simple') then
       ! Record that method is selected.
       methodSelected=.true.

       ! Increment the component count and store the value for later reference.
       componentTypeCount=componentTypeCount+1
       componentIndex=componentTypeCount

       ! Display message.
       message='Simple satellite orbit method selected [component index '
       message=message//componentIndex//']'
       call Galacticus_Display_Message(message,verbosityInfo)

       ! Set up procedure pointers.
       Tree_Node_Satellite_Merge_Time              => Tree_Node_Satellite_Merge_Time_Simple
       Tree_Node_Satellite_Merge_Time_Set          => Tree_Node_Satellite_Merge_Time_Set_Simple
       Tree_Node_Satellite_Merge_Time_Rate_Adjust  => Tree_Node_Satellite_Merge_Time_Rate_Adjust_Simple
       Tree_Node_Satellite_Merge_Time_Rate_Compute => Tree_Node_Satellite_Merge_Time_Rate_Compute_Simple
    end if

    ! Get the satellite merging timescale method.
    !@ <inputParameter>
    !@   <name>satelliteMergingMethod</name>
    !@   <defaultValue>Jiang2008</defaultValue>
    !@   <attachedTo>module</attachedTo>
    !@   <description>
    !@     The name of the method to be used to compute satellite merging timescales.
    !@   </description>
    !@ </inputParameter>
    call Get_Input_Parameter('satelliteMergingMethod',satelliteMergingMethod,defaultValue='Jiang2008')
    ! Include file that makes calls to all available method initialization routines.
    !# <include directive="satelliteMergingMethod" type="code" action="subroutine">
    !#  <subroutineArgs>satelliteMergingMethod,Satellite_Time_Until_Merging</subroutineArgs>
    include 'objects.tree_node.methods.satellite_orbit.inc'
    !# </include>
    if (.not.associated(Satellite_Time_Until_Merging)) call&
         & Galacticus_Error_Report('Tree_Node_Methods_Satellite_Orbit_Initialize','method '//char(satelliteMergingMethod)//' is&
         & unrecognized')

    return
  end subroutine Tree_Node_Methods_Satellite_Orbit_Initialize
  

  double precision function Tree_Node_Satellite_Merge_Time_Simple(thisNode)
    !% Return the time until satellite merging.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    if (thisNode%componentExists(componentIndex)) then 
       thisIndex=Tree_Node_Satellite_Orbit_Index(thisNode)
       Tree_Node_Satellite_Merge_Time_Simple=thisNode%components(thisIndex)%properties(mergeTimeIndex,propertyValue)
    else
       Tree_Node_Satellite_Merge_Time_Simple=-1.0d0 ! Negative time indicates that this is not a satellite.
    end if
    return
  end function Tree_Node_Satellite_Merge_Time_Simple

  subroutine Tree_Node_Satellite_Merge_Time_Set_Simple(thisNode,mergeTime)
    !% Set the time until satellite merging.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: mergeTime
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Satellite_Orbit_Index(thisNode)
    thisNode%components(thisIndex)%properties(mergeTimeIndex,propertyValue)=mergeTime
    return
  end subroutine Tree_Node_Satellite_Merge_Time_Set_Simple

  subroutine Tree_Node_Satellite_Merge_Time_Rate_Adjust_Simple(thisNode,interrupt,interruptProcedure,rateAdjustment)
    !% Return the time until satellite merging rate of change.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(),      pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Satellite_Orbit_Index(thisNode)
    thisNode%components(thisIndex)%properties(mergeTimeIndex,propertyDerivative) &
         &=thisNode%components(thisIndex)%properties(mergeTimeIndex,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Satellite_Merge_Time_Rate_Adjust_Simple

  subroutine Tree_Node_Satellite_Merge_Time_Rate_Compute_Simple(thisNode,interrupt,interruptProcedure)
    !% Compute the time until satellite merging rate of change.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    logical,                 intent(inout) :: interrupt
    procedure(),    pointer, intent(inout) :: interruptProcedure
    
    if (thisNode%componentExists(componentIndex)) then 
       if (Tree_Node_Satellite_Merge_Time_Simple(thisNode) <= 0.0d0) then
          ! Merger has happened - interrupt the evolution.
          interrupt=.true.
          interruptProcedure => Satellite_Merger_Process
          return
       end if
       call Tree_Node_Satellite_Merge_Time_Rate_Adjust_Simple(thisNode,interrupt,interruptProcedure,-1.0d0)
    end if
    return
  end subroutine Tree_Node_Satellite_Merge_Time_Rate_Compute_Simple


  !# <nodeMergerTask>
  !#  <unitName>Satellite_Orbit_Create</unitName>
  !# </nodeMergerTask>
  !# <satelliteHostChangeTask>
  !#  <unitName>Satellite_Orbit_Create</unitName>
  !# </satelliteHostChangeTask>
  subroutine Satellite_Orbit_Create(thisNode)
    !% Create a satellite orbit component and assign a time until merging.
    implicit none
    type(treeNode),  pointer, intent(inout) :: thisNode
    double precision                        :: mergeTime

    if (methodSelected) then
       ! Create a satellite orbit component and assign a time until merging.
      call thisNode%createComponent(componentIndex,propertyCount,dataCount,historyCount)
       mergeTime=Satellite_Time_Until_Merging(thisNode)
       call Tree_Node_Satellite_Merge_Time_Set_Simple(thisNode,mergeTime)
    end if
    return
  end subroutine Satellite_Orbit_Create


  integer function Tree_Node_Satellite_Orbit_Index(thisNode)
    !% Ensure the satellite orbit component exists and return its position in the components array.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    
    ! Create the component.
    call thisNode%createComponent(componentIndex,propertyCount,dataCount,historyCount)
    Tree_Node_Satellite_Orbit_Index=thisNode%componentIndex(componentIndex)
    return
  end function Tree_Node_Satellite_Orbit_Index

  subroutine Satellite_Merger_Process(thisNode)
    !% Process a satellite node which has undergone a merger with its host node.
    !# <include directive="satelliteMergerTask" type="moduleUse">
    include 'objects.tree_node.methods.satellite_orbit.merger.process.moduleUse.inc'
    !# </include>
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode

    ! Allow arbitrary routines to process the merger.
    !# <include directive="satelliteMergerTask" type="code" action="subroutine">
    !#  <subroutineArgs>thisNode</subroutineArgs>
    include 'objects.tree_node.methods.satellite_orbit.merger.process.inc'
    !# </include>

    ! Finally remove the satellite node from the host and destroy it.
    call thisNode%removeFromHost()
    call thisNode%destroy
    return
  end subroutine Satellite_Merger_Process
  
  
  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Satellite_Orbit_Simple_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Satellite_Orbit_Simple</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Satellite_Orbit_Simple_Names(integerProperty,integerPropertyNames,integerPropertyComments&
       &,doubleProperty,doublePropertyNames,doublePropertyComments,time)
    !% Set names of satellite orbit properties to be written to the \glc\ output file.
    implicit none
    double precision, intent(in)                  :: time
    integer,          intent(inout)               :: integerProperty,doubleProperty
    character(len=*), intent(inout), dimension(:) :: integerPropertyNames,integerPropertyComments,doublePropertyNames &
         &,doublePropertyComments
    
    if (methodSelected) then
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='timeToMerge'
       doublePropertyComments(doubleProperty)='Time until satellite merges.'
    end if
    return
  end subroutine Galacticus_Output_Tree_Satellite_Orbit_Simple_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Satellite_Orbit_Simple_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Satellite_Orbit_Simple</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Satellite_Orbit_Simple_Property_Count(integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of satellite orbit properties to be written to the the \glc\ output file.
    implicit none
    double precision, intent(in)    :: time
    integer,          intent(inout) :: integerPropertyCount,doublePropertyCount

    if (methodSelected) doublePropertyCount=doublePropertyCount+1
    return
  end subroutine Galacticus_Output_Tree_Satellite_Orbit_Simple_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Satellite_Orbit_Simple</unitName>
  !#  <sortName>Galacticus_Output_Tree_Satellite_Orbit_Simple</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Satellite_Orbit_Simple(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store satellite orbit properties in the \glc\ output file buffers.
    use Tree_Nodes
    use Tree_Node_Methods
    implicit none
    double precision, intent(in)             :: time
    type(treeNode),   intent(inout), pointer :: thisNode
    integer,          intent(inout)          :: integerProperty,integerBufferCount,doubleProperty,doubleBufferCount
    integer,          intent(inout)          :: integerBuffer(:,:)
    double precision, intent(inout)          :: doubleBuffer(:,:)

    if (methodSelected) then
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Satellite_Merge_Time(thisNode)
    end if
    return
  end subroutine Galacticus_Output_Tree_Satellite_Orbit_Simple

  !# <nodeDumpTask>
  !#  <unitName>Tree_Node_Methods_Satellite_Orbit_Simple_Dump</unitName>
  !# </nodeDumpTask>
  subroutine Tree_Node_Methods_Satellite_Orbit_Simple_Dump(thisNode)
    !% Dump all properties of {\tt thisNode} to screen.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    if (methodSelected) then
       if (thisNode%componentExists(componentIndex)) then
          write (0,'(1x,a)'           ) 'satellite orbit component -> properties:'
          write (0,'(2x,a50,1x,e12.6)') 'time until merging:',Tree_Node_Satellite_Merge_Time(thisNode)
       else
          write (0,'(1x,a)'           ) 'satellite orbit component -> nonexistant'
       end if
    end if
    return
  end subroutine Tree_Node_Methods_Satellite_Orbit_Simple_Dump

end module Tree_Node_Methods_Satellite_Orbit
