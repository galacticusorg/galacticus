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


!% Contains a module of node position methods in which properties are preset.

module Tree_Node_Methods_Positions_Preset
  !% Implements node positions using preset data.
  use Tree_Nodes
  use Components
  private
  public :: Tree_Node_Methods_Position_Initialize_Preset, Galacticus_Output_Tree_Position_Preset,&
       & Galacticus_Output_Tree_Position_Preset_Property_Count, Galacticus_Output_Tree_Position_Preset_Names,&
       & Tree_Node_Methods_Position_Preset_Dump, Tree_Node_Position_Preset_Promote
  
  ! The index used as a reference for this component.
  integer :: componentIndex=-1

  ! Property indices.
  integer, parameter :: propertyCount=0, dataCount=6, historyCountInitial=0, historyCount=1
  integer, parameter :: xPositionIndex      =1
  integer, parameter :: yPositionIndex      =2
  integer, parameter :: zPositionIndex      =3
  integer, parameter :: xVelocityIndex      =4
  integer, parameter :: yVelocityIndex      =5
  integer, parameter :: zVelocityIndex      =6
  integer, parameter :: positionHistoryIndex=1

  ! Define procedure pointers.
  !# <treeNodeMethodsPointer type="array">
  !#  <methodName>Tree_Node_Position</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer type="array">
  !#  <methodName>Tree_Node_Velocity</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer type="history">
  !#  <methodName>Tree_Node_Position_6D_History</methodName>
  !# </treeNodeMethodsPointer>

  ! Flag to indicate if this method is selected.
  logical :: methodSelected=.false.

contains

  !# <treeNodeCreateInitialize>
  !#  <unitName>Tree_Node_Methods_Position_Initialize_Preset</unitName>
  !#  <optionName>treeNodeMethodPosition</optionName>
  !# </treeNodeCreateInitialize>
  subroutine Tree_Node_Methods_Position_Initialize_Preset(componentOption,componentTypeCount)
    !% Initializes the tree node preset position methods module.
    use ISO_Varying_String
    use Input_Parameters
    use Galacticus_Error
    use Galacticus_Display
    use String_Handling
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
       message='Preset position method selected [component index '
       message=message//componentIndex//']'
       call Galacticus_Display_Message(message,verbosityInfo)

       ! Set up procedure pointers.
       Tree_Node_Position                         => Tree_Node_Position_Preset
       Tree_Node_Position_Set                     => Tree_Node_Position_Set_Preset
       Tree_Node_Position_Rate_Adjust             => null()
       Tree_Node_Position_Rate_Compute            => Tree_Node_Rate_Rate_Compute_Dummy
       Tree_Node_Velocity                         => Tree_Node_Velocity_Preset
       Tree_Node_Velocity_Set                     => Tree_Node_Velocity_Set_Preset
       Tree_Node_Velocity_Rate_Adjust             => null()
       Tree_Node_Velocity_Rate_Compute            => Tree_Node_Rate_Rate_Compute_Dummy
       Tree_Node_Position_6D_History              => Tree_Node_Position_6D_History_Preset
       Tree_Node_Position_6D_History_Set          => Tree_Node_Position_6D_History_Set_Preset
       Tree_Node_Position_6D_History_Rate_Adjust  => null()
       Tree_Node_Position_6D_History_Rate_Compute => Tree_Node_Rate_Rate_Compute_Dummy
    end if

    return
  end subroutine Tree_Node_Methods_Position_Initialize_Preset
  
  subroutine Tree_Node_Position_Preset(thisNode,position)
    !% Return the position of the node.
    use Numerical_Interpolation
    use FGSL
    implicit none
    double precision, dimension(:), intent(out)   :: position
    type(treeNode),   pointer,      intent(inout) :: thisNode
    integer                                       :: thisIndex,iTime
    logical                                       :: usingHistory,interpolationReset
    type(fgsl_interp_accel)                       :: interpolationAccelerator

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Position_Index(thisNode)
       ! Check if this node has a position history attached to it.
       usingHistory=.false.
       if (allocated(thisNode%components(thisIndex)%histories)) then
          if (thisNode%components(thisIndex)%histories(positionHistoryIndex)%time(1) <= Tree_Node_Time(thisNode)) then
             interpolationReset=.true.
             iTime=Interpolate_Locate(size(thisNode%components(thisIndex)%histories(positionHistoryIndex)%time)&
                  &,thisNode%components(thisIndex)%histories(positionHistoryIndex)%time,interpolationAccelerator&
                  &,Tree_Node_Time(thisNode),reset=interpolationReset,closest=.true.)
             position=thisNode%components(thisIndex)%histories(positionHistoryIndex)%data(iTime,xPositionIndex:zPositionIndex)
             call Interpolate_Done(interpolationAccelerator=interpolationAccelerator,reset=interpolationReset)  
             usingHistory=.true.
          end if
       end if
       if (.not.usingHistory) position=thisNode%components(thisIndex)%data(xPositionIndex:zPositionIndex)
    else
       position=[0.0d0,0.0d0,0.0d0]
    end if
    return
  end subroutine Tree_Node_Position_Preset

  subroutine Tree_Node_Position_Set_Preset(thisNode,position)
    !% Set the node position.
    implicit none
    type(treeNode),   pointer,      intent(inout) :: thisNode
    double precision, dimension(:), intent(in)    :: position
    integer                                       :: thisIndex

    thisIndex=Tree_Node_Position_Index(thisNode)
    thisNode%components(thisIndex)%data(xPositionIndex:zPositionIndex)=position
    return
  end subroutine Tree_Node_Position_Set_Preset
  
  subroutine Tree_Node_Velocity_Preset(thisNode,velocity)
    !% Return the position of the node.
    use Numerical_Interpolation
    use FGSL
    implicit none
    type(treeNode),   pointer,      intent(inout) :: thisNode
    double precision, dimension(:), intent(out)   :: velocity
    integer                                       :: thisIndex,iTime
    logical                                       :: usingHistory,interpolationReset
    type(fgsl_interp_accel)                       :: interpolationAccelerator

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Position_Index(thisNode)
       ! Check if this node has a position history attached to it.
       usingHistory=.false.
       if (allocated(thisNode%components(thisIndex)%histories)) then
          if (thisNode%components(thisIndex)%histories(positionHistoryIndex)%time(1) <= Tree_Node_Time(thisNode)) then
             interpolationReset=.true.
             iTime=Interpolate_Locate(size(thisNode%components(thisIndex)%histories(positionHistoryIndex)%time)&
                  &,thisNode%components(thisIndex)%histories(positionHistoryIndex)%time,interpolationAccelerator&
                  &,Tree_Node_Time(thisNode),reset=interpolationReset,closest=.true.)
             velocity=thisNode%components(thisIndex)%histories(positionHistoryIndex)%data(iTime,xVelocityIndex:zVelocityIndex)
             call Interpolate_Done(interpolationAccelerator=interpolationAccelerator,reset=interpolationReset)  
             usingHistory=.true.
          end if
       end if
       if (.not.usingHistory) velocity=thisNode%components(thisIndex)%data(xVelocityIndex:zVelocityIndex)
    else
       velocity=[0.0d0,0.0d0,0.0d0]
    end if
    return
  end subroutine Tree_Node_Velocity_Preset

  subroutine Tree_Node_Velocity_Set_Preset(thisNode,position)
    !% Set the node position.
    implicit none
    type(treeNode),   pointer,      intent(inout) :: thisNode
    double precision, dimension(:), intent(in)    :: position
    integer                                       :: thisIndex

    thisIndex=Tree_Node_Position_Index(thisNode)
    thisNode%components(thisIndex)%data(xVelocityIndex:zVelocityIndex)=position
    return
  end subroutine Tree_Node_Velocity_Set_Preset

  function Tree_Node_Position_6D_History_Preset(thisNode)
    !% Return the satellite 6D position history.
    use Histories
    implicit none
    type(history)                          :: Tree_Node_Position_6D_History_Preset
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    ! Check if the component exists.
    if (thisNode%componentExists(componentIndex)) then
       ! It does, so get the index of the component.
       thisIndex=Tree_Node_Position_Index(thisNode)
       ! Check if the component has histories
       if (allocated(thisNode%components(thisIndex)%histories)) then
          ! It does, so return that history.
         Tree_Node_Position_6D_History_Preset =thisNode%components(thisIndex)%histories(positionHistoryIndex)
       else
          ! It does not, so return a null history.
          Tree_Node_Position_6D_History_Preset=nullHistory
       end if
    else
       ! Component does not exist, so return a null history.
      Tree_Node_Position_6D_History_Preset=nullHistory
    end if
    return
  end function Tree_Node_Position_6D_History_Preset

  subroutine Tree_Node_Position_6D_History_Set_Preset(thisNode,thisHistory)
    !% Set the node 6D position properties history.
    use Memory_Management
    use Histories
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    type(history),           intent(in)    :: thisHistory
    integer                                :: thisIndex

    ! Get the index of this component.
    thisIndex=Tree_Node_Position_Index(thisNode)
     ! If this component does not yet have a history associated with it, then create one.
     if (.not.allocated(thisNode%components(thisIndex)%histories)) then
        allocate(thisNode%components(thisIndex)%histories(historyCount))
        call Memory_Usage_Record(sizeof(thisNode%components(thisIndex)%histories),memoryType=memoryTypeNodes)
     end if
     ! Destroy the current history.
     call thisNode%components(thisIndex)%histories(positionHistoryIndex)%destroy()
     ! Assign the new history.
     thisNode%components(thisIndex)%histories(positionHistoryIndex)=thisHistory
    return
  end subroutine Tree_Node_Position_6D_History_Set_Preset

  integer function Tree_Node_Position_Index(thisNode)
    !% Ensure the position component exists and return its position in the components array.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    
    ! Create the component.
    call thisNode%createComponent(componentIndex,propertyCount,dataCount,historyCountInitial)
    Tree_Node_Position_Index=thisNode%componentIndex(componentIndex)
    return
  end function Tree_Node_Position_Index

  !# <nodePromotionTask>
  !#  <unitName>Tree_Node_Position_Preset_Promote</unitName>
  !# </nodePromotionTask>
  subroutine Tree_Node_Position_Preset_Promote(thisNode)
    !% Ensure that {\tt thisNode} is ready for promotion to its parent. In this case, update the position of {\tt
    !% thisNode} to that of the parent.
    use Histories
    implicit none
    type(treeNode),   pointer, intent(inout)     :: thisNode
    type(treeNode),   pointer                    :: parentNode
    double precision, dimension(3)               :: position,velocity
    integer                                      :: parentIndex
    type(history)                                :: positionHistory

    if (methodSelected) then
       parentNode => thisNode%parentNode
       if (parentNode%componentExists(componentIndex)) then
          call Tree_Node_Position    (parentNode,position)
          call Tree_Node_Position_Set(thisNode  ,position)
          call Tree_Node_Velocity    (parentNode,velocity)
          call Tree_Node_Velocity_Set(thisNode  ,velocity)
          ! Check if the parent node has a position history attached to it.
          parentIndex=Tree_Node_Position_Index(parentNode)
          if (allocated(parentNode%components(parentIndex)%histories)) then
             ! It does, so transfer it to the active node.
             positionHistory=Tree_Node_Position_6D_History(parentNode)
             call Tree_Node_Position_6D_History_Set(thisNode,positionHistory)
             call positionHistory%destroy()
          end if
       end if
    end if
    return
  end subroutine Tree_Node_Position_Preset_Promote
  
  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Position_Preset_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Position_Preset</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Position_Preset_Names(integerProperty,integerPropertyNames,integerPropertyComments&
       &,integerPropertyUnitsSI,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set names of position properties to be written to the \glc\ output file.
    use Numerical_Constants_Astronomical
    implicit none
    double precision, intent(in)                  :: time
    integer,          intent(inout)               :: integerProperty,doubleProperty
    character(len=*), intent(inout), dimension(:) :: integerPropertyNames,integerPropertyComments,doublePropertyNames &
         &,doublePropertyComments
    double precision, intent(inout), dimension(:) :: integerPropertyUnitsSI,doublePropertyUnitsSI
    
    if (methodSelected) then
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='positionX'
       doublePropertyComments(doubleProperty)='X position of the node.'
       doublePropertyUnitsSI (doubleProperty)=megaParsec
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='positionY'
       doublePropertyComments(doubleProperty)='Y position of the node.'
       doublePropertyUnitsSI (doubleProperty)=megaParsec
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='positionZ'
       doublePropertyComments(doubleProperty)='Z position of the node.'
       doublePropertyUnitsSI (doubleProperty)=megaParsec
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='velocityX'
       doublePropertyComments(doubleProperty)='X velocity of the node.'
       doublePropertyUnitsSI (doubleProperty)=kilo
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='velocityY'
       doublePropertyComments(doubleProperty)='Y velocity of the node.'
       doublePropertyUnitsSI (doubleProperty)=kilo
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='velocityZ'
       doublePropertyComments(doubleProperty)='Z velocity of the node.'
       doublePropertyUnitsSI (doubleProperty)=kilo
    end if
    return
  end subroutine Galacticus_Output_Tree_Position_Preset_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Position_Preset_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Position_Preset</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Position_Preset_Property_Count(integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of position properties to be written to the the \glc\ output file.
    implicit none
    double precision, intent(in)    :: time
    integer,          intent(inout) :: integerPropertyCount,doublePropertyCount

    if (methodSelected) doublePropertyCount=doublePropertyCount+6
    return
  end subroutine Galacticus_Output_Tree_Position_Preset_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Position_Preset</unitName>
  !#  <sortName>Galacticus_Output_Tree_Position_Preset</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Position_Preset(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store position properties in the \glc\ output file buffers.
    use Histories
    use Tree_Nodes
    implicit none
    double precision, intent(in)             :: time
    type(treeNode),   intent(inout), pointer :: thisNode
    integer,          intent(inout)          :: integerProperty,integerBufferCount,doubleProperty,doubleBufferCount
    integer,          intent(inout)          :: integerBuffer(:,:)
    double precision, intent(inout)          :: doubleBuffer(:,:)
    double precision, dimension(3)           :: position,velocity

    if (methodSelected) then
       call Tree_Node_Position(thisNode,position)
       call Tree_Node_Velocity(thisNode,velocity)
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=position(1)
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=position(2)
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=position(3)
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=velocity(1)
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=velocity(2)
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=velocity(3)
    end if
    return
  end subroutine Galacticus_Output_Tree_Position_Preset

  !# <nodeDumpTask>
  !#  <unitName>Tree_Node_Methods_Position_Preset_Dump</unitName>
  !# </nodeDumpTask>
  subroutine Tree_Node_Methods_Position_Preset_Dump(thisNode)
    !% Dump all properties of {\tt thisNode} to screen.
    use Histories
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, dimension(3)           :: position,velocity

    if (methodSelected) then
       if (thisNode%componentExists(componentIndex)) then
          write (0,'(1x,a)'           ) 'position component -> properties:'
          call Tree_Node_Position(thisNode,position)
          call Tree_Node_Velocity(thisNode,velocity)
          write (0,'(2x,a50,1x,e12.6)') 'position:',position
          write (0,'(2x,a50,1x,e12.6)') 'velocity:',velocity
       else
          write (0,'(1x,a)'           ) 'position component -> nonexistant'
       end if
    end if
    return
  end subroutine Tree_Node_Methods_Position_Preset_Dump

end module Tree_Node_Methods_Positions_Preset
