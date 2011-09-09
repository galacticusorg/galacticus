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


!% Contains a module of merging statistics methods.

module Tree_Node_Methods_Merging_Stats_Standard
  !% Implement random spin tree node method.
  use Tree_Nodes
  use Components
  implicit none
  private
  public :: Tree_Node_Methods_Merging_Stats_Standard_Initialize, Galacticus_Output_Tree_Merging_Stats_Standard,&
       & Galacticus_Output_Tree_Merging_Stats_Standard_Property_Count, Galacticus_Output_Tree_Merging_Stats_Standard_Names,&
       & Tree_Node_Methods_Merging_Stats_Standard_Dump, Tree_Node_Merging_Stats_Standard_Record,&
       & Tree_Node_Merging_Stats_Standard_Node_Record, Tree_Node_Merging_Stats_Standard_Node_Promote,&
       & Tree_Node_Merging_Stats_Standard_Initialize
  
  ! The index used as a reference for this component.
  integer :: componentIndex=-1

  ! Property indices.
  integer, parameter :: propertyCount=0, dataCount=3, historyCount=0
  integer, parameter :: majorMergerTimeIndex    =1
  integer, parameter :: nodeMajorMergerTimeIndex=2
  integer, parameter :: nodeFormationTimeIndex  =3

  ! Define procedure pointers.
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Galaxy_Major_Merger_Time</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Node_Major_Merger_Time</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Node_Formation_Time</methodName>
  !# </treeNodeMethodsPointer>

  ! Flag to indicate if this method is selected.
  logical          :: methodSelected=.false.

  ! Flag to indicate whether or not to track merging statistics.
  logical          :: trackMergerStatistics=.false.

  ! Parameters controlling the statistics gathered.
  double precision :: nodeMajorMergerFraction,nodeFormationMassFraction

contains

  !# <treeNodeCreateInitialize>
  !#  <unitName>Tree_Node_Methods_Merging_Stats_Standard_Initialize</unitName>
  !#  <optionName default="standard">treeNodeMethodMergingStatistics</optionName>
  !# </treeNodeCreateInitialize>
  subroutine Tree_Node_Methods_Merging_Stats_Standard_Initialize(componentOption,componentTypeCount)
    !% Initializes the tree node random spin method module.
    use ISO_Varying_String
    use Input_Parameters
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
       message='Standard merging statistics method selected [component index '
       message=message//componentIndex//']'
       call Galacticus_Display_Message(message,verbosityInfo)

       ! Set up procedure pointers.
       Tree_Node_Galaxy_Major_Merger_Time              => Tree_Node_Galaxy_Major_Merger_Time_Merging_Stats_Standard
       Tree_Node_Galaxy_Major_Merger_Time_Set          => null()
       Tree_Node_Galaxy_Major_Merger_Time_Rate_Adjust  => null()
       Tree_Node_Galaxy_Major_Merger_Time_Rate_Compute => Tree_Node_Rate_Rate_Compute_Dummy
       Tree_Node_Node_Major_Merger_Time                => Tree_Node_Node_Major_Merger_Time_Merging_Stats_Standard
       Tree_Node_Node_Major_Merger_Time_Set            => null()
       Tree_Node_Node_Major_Merger_Time_Rate_Adjust    => null()
       Tree_Node_Node_Major_Merger_Time_Rate_Compute   => Tree_Node_Rate_Rate_Compute_Dummy
       Tree_Node_Node_Formation_Time                   => Tree_Node_Node_Formation_Time_Merging_Stats_Standard
       Tree_Node_Node_Formation_Time_Set               => null()
       Tree_Node_Node_Formation_Time_Rate_Adjust       => null()
       Tree_Node_Node_Formation_Time_Rate_Compute      => Tree_Node_Rate_Rate_Compute_Dummy

       !@ <inputParameter>
       !@   <name>trackMergerStatistics</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Indicates whether or not galaxy and node merger statistics should be tracked.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('trackMergerStatistics',trackMergerStatistics,defaultValue=.false.)
       !@ <inputParameter>
       !@   <name>nodeMajorMergerFraction</name>
       !@   <defaultValue>0.25</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The mass ratio ($M_2/M_1$ where $M_2 &lt; M_1$) of merging halos above which the merger should be considered to be ``major''.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('nodeMajorMergerFraction',nodeMajorMergerFraction,defaultValue=0.25d0)
       !@ <inputParameter>
       !@   <name>nodeFormationMassFraction</name>
       !@   <defaultValue>0.5</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The mass fraction in the main branch progenitor used to define the formation time of each halo.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('nodeFormationMassFraction',nodeFormationMassFraction,defaultValue=0.5d0)
    end if
    return
  end subroutine Tree_Node_Methods_Merging_Stats_Standard_Initialize

  double precision function Tree_Node_Galaxy_Major_Merger_Time_Merging_Stats_Standard(thisNode,instance)
    !% Return the time of the last major merger.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex
    
    if (trackMergerStatistics.and.thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Merging_Stats_Standard_Index(thisNode)
       Tree_Node_Galaxy_Major_Merger_Time_Merging_Stats_Standard=thisNode%components(thisIndex)%instance(1)%data(majorMergerTimeIndex)
    else
       Tree_Node_Galaxy_Major_Merger_Time_Merging_Stats_Standard=-1.0d0
    end if
    return
  end function Tree_Node_Galaxy_Major_Merger_Time_Merging_Stats_Standard

  double precision function Tree_Node_Node_Major_Merger_Time_Merging_Stats_Standard(thisNode,instance)
    !% Return the time of the last major merger of nodes.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex
    
    if (trackMergerStatistics.and.thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Merging_Stats_Standard_Index(thisNode)
       Tree_Node_Node_Major_Merger_Time_Merging_Stats_Standard=thisNode%components(thisIndex)%instance(1)%data(nodeMajorMergerTimeIndex)
    else
       Tree_Node_Node_Major_Merger_Time_Merging_Stats_Standard=-1.0d0
    end if
    return
  end function Tree_Node_Node_Major_Merger_Time_Merging_Stats_Standard

  double precision function Tree_Node_Node_Formation_Time_Merging_Stats_Standard(thisNode,instance)
    !% Return the time of the last major merger of nodes.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex
    
    if (trackMergerStatistics) then
       thisIndex=Tree_Node_Merging_Stats_Standard_Index(thisNode)
       Tree_Node_Node_Formation_Time_Merging_Stats_Standard=thisNode%components(thisIndex)%instance(1)%data(nodeFormationTimeIndex)
    else
       Tree_Node_Node_Formation_Time_Merging_Stats_Standard=-1.0d0
    end if
    return
  end function Tree_Node_Node_Formation_Time_Merging_Stats_Standard

  !# <mergerTreeInitializeTask>
  !#  <unitName>Tree_Node_Merging_Stats_Standard_Initialize</unitName>
  !# </mergerTreeInitializeTask>
  subroutine Tree_Node_Merging_Stats_Standard_Initialize(thisNode)
    !% Initialize the merging statistics component by creating components in nodes and computing formation times.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    ! Simply call the indexing function as this will create the component and initialize it.
    if (trackMergerStatistics) thisIndex=Tree_Node_Merging_Stats_Standard_Index(thisNode)
    return
  end subroutine Tree_Node_Merging_Stats_Standard_Initialize

  integer function Tree_Node_Merging_Stats_Standard_Index(thisNode)
    !% Ensure the merging statistics component exists and return its position in the components array.
    use Dark_Matter_Halo_Formation_Times
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode

    if (.not.thisNode%componentExists(componentIndex)) then
       call thisNode%createComponent(componentIndex,propertyCount,dataCount,historyCount)
       Tree_Node_Merging_Stats_Standard_Index=thisNode%componentIndex(componentIndex)
       thisNode%components(Tree_Node_Merging_Stats_Standard_Index)%instance(1)%data(majorMergerTimeIndex    )=-1.0d0
       thisNode%components(Tree_Node_Merging_Stats_Standard_Index)%instance(1)%data(nodeMajorMergerTimeIndex)=-1.0d0
       thisNode%components(Tree_Node_Merging_Stats_Standard_Index)%instance(1)%data(nodeFormationTimeIndex  )=Dark_Matter_Halo_Formation_Time(thisNode,nodeFormationMassFraction)
    else
       Tree_Node_Merging_Stats_Standard_Index=thisNode%componentIndex(componentIndex)
    end if
    return
  end function Tree_Node_Merging_Stats_Standard_Index
  
  !# <nodeMergerTask>
  !#  <unitName>Tree_Node_Merging_Stats_Standard_Node_Record</unitName>
  !# </nodeMergerTask>
  subroutine Tree_Node_Merging_Stats_Standard_Node_Record(thisNode)
    !% Record any major merger of {\tt thisNode}.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex
    
    if (methodSelected.and.trackMergerStatistics) then
       if (Tree_Node_Mass(thisNode) >= nodeMajorMergerFraction*Tree_Node_Mass(thisNode%parentNode)) then
          ! Record the merger time.
          thisIndex=Tree_Node_Merging_Stats_Standard_Index(thisNode%parentNode)
          thisNode%parentNode%components(thisIndex)%instance(1)%data(nodeMajorMergerTimeIndex)=Tree_Node_Time(thisNode)
       end if
    end if
    return
  end subroutine Tree_Node_Merging_Stats_Standard_Node_Record
  
  !# <nodePromotionTask>
  !#  <unitName>Tree_Node_Merging_Stats_Standard_Node_Promote</unitName>
  !# </nodePromotionTask>
  subroutine Tree_Node_Merging_Stats_Standard_Node_Promote(thisNode)
    !% Ensure that {\tt thisNode} is ready for promotion to its parent. In this case, we simply update the node merger time.
    use Galacticus_Error
    implicit none
    type(treeNode),   pointer, intent(inout)     :: thisNode
    type(treeNode),   pointer                    :: parentNode
    integer                                      :: thisIndex

    ! Check if this method is selected.
    if (methodSelected.and.trackMergerStatistics) then       
       ! Get the parent node of this node.
       parentNode => thisNode%parentNode
       ! If the parent node has a merging statistic component, then copy the time of the last node merger.
       if (parentNode%componentExists(componentIndex)) then
          thisIndex=Tree_Node_Merging_Stats_Standard_Index(thisNode)
          if (Tree_Node_Node_Major_Merger_Time(parentNode) > Tree_Node_Node_Major_Merger_Time(thisNode)) &
               & thisNode%components(thisIndex)%instance(1)%data(nodeMajorMergerTimeIndex)=Tree_Node_Node_Major_Merger_Time(parentNode)
          thisNode%components(thisIndex)%instance(1)%data(nodeFormationTimeIndex)=Tree_Node_Node_Formation_Time(parentNode)
       end if
    end if
    return
  end subroutine Tree_Node_Merging_Stats_Standard_Node_Promote

  !# <satelliteMergerTask>
  !#  <unitName>Tree_Node_Merging_Stats_Standard_Record</unitName>
  !#  <after>Satellite_Merging_Mass_Movement_Store</after>
  !# </satelliteMergerTask>
  subroutine Tree_Node_Merging_Stats_Standard_Record(thisNode)
    !% Record properties of a merging event for {\tt thisNode}.
    use Satellite_Merging_Mass_Movements_Descriptors
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    type(treeNode), pointer                :: hostNode
    integer                                :: thisIndex

    ! Determine if we should record.
    if (methodSelected.and.trackMergerStatistics) then       
       ! Record the time of this merger if it is a major merger.
       if (thisMergerIsMajor) then          
          ! Find the node to merge with.
          call thisNode%mergesWith(hostNode)
          ! Record the merger time.
          thisIndex=Tree_Node_Merging_Stats_Standard_Index(hostNode)
          hostNode%components(thisIndex)%instance(1)%data(majorMergerTimeIndex)=Tree_Node_Time(hostNode)
       end if
    end if
    return
  end subroutine Tree_Node_Merging_Stats_Standard_Record

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Merging_Stats_Standard_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Merging_Stats_Standard</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Merging_Stats_Standard_Names(integerProperty,integerPropertyNames,integerPropertyComments&
       &,integerPropertyUnitsSI,doubleProperty ,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set names of spin properties to be written to the \glc\ output file.
    use Numerical_Constants_Astronomical
    implicit none
    double precision, intent(in)                  :: time
    integer,          intent(inout)               :: integerProperty,doubleProperty
    character(len=*), intent(inout), dimension(:) :: integerPropertyNames,integerPropertyComments,doublePropertyNames &
         &,doublePropertyComments
    double precision, intent(inout), dimension(:) :: integerPropertyUnitsSI,doublePropertyUnitsSI

    if (methodSelected.and.trackMergerStatistics) then
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='majorMergerTimeLapse'
       doublePropertyComments(doubleProperty)='Time since the last major merger.'
       doublePropertyUnitsSI (doubleProperty)=gigaYear
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='nodeMajorMergerTimeLapse'
       doublePropertyComments(doubleProperty)='Time since the last node major merger.'
       doublePropertyUnitsSI (doubleProperty)=gigaYear
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='nodeFormationTime'
       doublePropertyComments(doubleProperty)='Formation time of the node.'
       doublePropertyUnitsSI (doubleProperty)=gigaYear
    end if
    return
  end subroutine Galacticus_Output_Tree_Merging_Stats_Standard_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Merging_Stats_Standard_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Merging_Stats_Standard</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Merging_Stats_Standard_Property_Count(integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of spin properties to be written to the the \glc\ output file.
    implicit none
    double precision, intent(in)    :: time
    integer,          intent(inout) :: integerPropertyCount,doublePropertyCount

    if (methodSelected.and.trackMergerStatistics) doublePropertyCount=doublePropertyCount+3
    return
  end subroutine Galacticus_Output_Tree_Merging_Stats_Standard_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Merging_Stats_Standard</unitName>
  !#  <sortName>Galacticus_Output_Tree_Merging_Stats_Standard</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Merging_Stats_Standard(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
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
    double precision                                :: mergerTime

    if (methodSelected.and.trackMergerStatistics) then
       doubleProperty=doubleProperty+1
       mergerTime=Tree_Node_Galaxy_Major_Merger_Time_Merging_Stats_Standard(thisNode)
       if (mergerTime >= 0.0d0) then
          doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Time(thisNode)-mergerTime
       else
          doubleBuffer(doubleBufferCount,doubleProperty)=                         mergerTime
       end if
       doubleProperty=doubleProperty+1
       mergerTime=Tree_Node_Node_Major_Merger_Time_Merging_Stats_Standard(thisNode)
       if (mergerTime >= 0.0d0) then
          doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Time(thisNode)-mergerTime
       else
          doubleBuffer(doubleBufferCount,doubleProperty)=                         mergerTime
       end if
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Node_Formation_Time(thisNode)
    end if
    return
  end subroutine Galacticus_Output_Tree_Merging_Stats_Standard

  !# <nodeDumpTask>
  !#  <unitName>Tree_Node_Methods_Merging_Stats_Standard_Dump</unitName>
  !# </nodeDumpTask>
  subroutine Tree_Node_Methods_Merging_Stats_Standard_Dump(thisNode)
    !% Dump all properties of {\tt thisNode} to screen.
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
 
    if (methodSelected.and.trackMergerStatistics) then
       if (thisNode%componentExists(componentIndex)) then
          write (0,'(1x,a)'           ) 'merging statistics component -> properties:'
          write (0,'(2x,a50,1x,e12.6)') '     time of last major merger:',Tree_Node_Galaxy_Major_Merger_Time(thisNode)
          write (0,'(2x,a50,1x,e12.6)') 'time of last node major merger:',Tree_Node_Node_Major_Merger_Time  (thisNode)
          write (0,'(2x,a50,1x,e12.6)') '        formation time of node:',Tree_Node_Node_Formation_Time     (thisNode)
       else
          write (0,'(1x,a)'           ) 'merging statistics component -> nonexistant'
       end if
    end if
    return
  end subroutine Tree_Node_Methods_Merging_Stats_Standard_Dump

end module Tree_Node_Methods_Merging_Stats_Standard
