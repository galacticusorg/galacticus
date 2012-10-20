!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!+    Contributions to this file made by:  Andrew Benson, Jianling Gan.

!% Contains a module of satellite orbit tree node methods in which properties are preset.

module Tree_Node_Methods_Satellite_Orbit_Preset
  !% Implement satellite orbit tree node methods in which properties are preset.
  use Tree_Nodes
  implicit none
  private
  public :: Tree_Node_Methods_Satellite_Orbit_Initialize_Preset, Satellite_Orbit_Create_Preset, Galacticus_Output_Tree_Satellite_Orbit_Preset,&
       & Galacticus_Output_Tree_Satellite_Orbit_Preset_Property_Count, Galacticus_Output_Tree_Satellite_Orbit_Preset_Names,&
       & Tree_Node_Methods_Satellite_Orbit_Preset_Dump, Tree_Node_Satellite_Preset_Promote
  
  ! The index used as a reference for this component.
  integer :: componentIndex=-1

  ! Property indices.
  integer, parameter :: propertyCount =0, dataCount=5, historyCountInitial=0, historyCount=1
  ! Indices for data properties.
  integer, parameter :: mergeTimeIndex         =1
  integer, parameter :: velocityRadialIndex    =2
  integer, parameter :: velocityTangentialIndex=3
  integer, parameter :: virialRadiusIndex      =4
  integer, parameter :: hostMassIndex          =5
  ! Indices for histories.
  integer, parameter :: boundMassIndex         =1

  ! Define procedure pointers.
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Satellite_Merge_Time</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Satellite_Time_Of_Merging</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Bound_Mass</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer type="history">
  !#  <methodName>Tree_Node_Bound_Mass_History</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer type="keplerOrbit">
  !#  <methodName>Tree_Node_Satellite_Virial_Orbit</methodName>
  !# </treeNodeMethodsPointer>

  ! Flag to indicate if this method is selected.
  logical :: methodSelected=.false.

contains

  !# <treeNodeCreateInitialize>
  !#  <unitName>Tree_Node_Methods_Satellite_Orbit_Initialize_Preset</unitName>
  !#  <optionName>treeNodeMethodSatelliteOrbit</optionName>
  !# </treeNodeCreateInitialize>
  subroutine Tree_Node_Methods_Satellite_Orbit_Initialize_Preset(componentOption,componentTypeCount)
    !% Initializes the tree node satellite orbit methods module.
    use ISO_Varying_String
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
       message='Preset satellite orbit method selected [component index '
       message=message//componentIndex//']'
       call Galacticus_Display_Message(message,verbosityInfo)

       ! Set up procedure pointers.
       Tree_Node_Satellite_Merge_Time              => Tree_Node_Satellite_Merge_Time_Preset
       Tree_Node_Satellite_Merge_Time_Set          => Tree_Node_Satellite_Merge_Time_Set_Preset
       Tree_Node_Satellite_Merge_Time_Rate_Adjust  => null()
       Tree_Node_Satellite_Merge_Time_Rate_Compute => Tree_Node_Rate_Rate_Compute_Dummy

       Tree_Node_Satellite_Time_of_Merging              => Tree_Node_Satellite_Time_of_Merging_Preset
       Tree_Node_Satellite_Time_of_Merging_Set          => Tree_Node_Satellite_Time_of_Merging_Set_Preset
       Tree_Node_Satellite_Time_of_Merging_Rate_Adjust  => null()
       Tree_Node_Satellite_Time_of_Merging_Rate_Compute => Tree_Node_Rate_Rate_Compute_Dummy

       Tree_Node_Bound_Mass                        => Tree_Node_Bound_Mass_Preset
       Tree_Node_Bound_Mass_Set                    => null()
       Tree_Node_Bound_Mass_Rate_Adjust            => null()
       Tree_Node_Bound_Mass_Rate_Compute           => Tree_Node_Rate_Rate_Compute_Dummy
 
       Tree_Node_Bound_Mass_History                => Tree_Node_Bound_Mass_History_Preset
       Tree_Node_Bound_Mass_History_Set            => Tree_Node_Bound_Mass_History_Set_Preset
       Tree_Node_Bound_Mass_History_Rate_Adjust    => null()
       Tree_Node_Bound_Mass_History_Rate_Compute   => Tree_Node_Rate_Rate_Compute_Dummy

       Tree_Node_Satellite_Virial_Orbit            => Tree_Node_Satellite_Virial_Orbit_Preset
       Tree_Node_Satellite_Virial_Orbit_Set        => Tree_Node_Satellite_Virial_Orbit_Set_Preset
    end if

    return
  end subroutine Tree_Node_Methods_Satellite_Orbit_Initialize_Preset
  
  double precision function Tree_Node_Satellite_Merge_Time_Preset(thisNode,instance)
    !% Return the time until satellite merging.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    integer,        intent(in), optional   :: instance
    integer                                :: thisIndex

    if (thisNode%componentExists(componentIndex).and.thisNode%isSatellite()) then
       thisIndex=Tree_Node_Satellite_Orbit_Index(thisNode)
       Tree_Node_Satellite_Merge_Time_Preset=max(thisNode%components(thisIndex)%instance(1)%data(mergeTimeIndex)&
            &-Tree_Node_Time(thisNode),0.0d0)
    else
       Tree_Node_Satellite_Merge_Time_Preset=-1.0d0 ! Negative time indicates that this is not a satellite.
    end if
    return
  end function Tree_Node_Satellite_Merge_Time_Preset

  subroutine Tree_Node_Satellite_Merge_Time_Set_Preset(thisNode,mergeTime,instance)
    !% Set the time until satellite merging.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: mergeTime
    integer,        intent(in), optional   :: instance
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Satellite_Orbit_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%data(mergeTimeIndex)=mergeTime+Tree_Node_Time(thisNode)

    return
  end subroutine Tree_Node_Satellite_Merge_Time_Set_Preset

  double precision function Tree_Node_Satellite_Time_of_Merging_Preset(thisNode,instance)
    !% Return the time of satellite merging.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    integer,        intent(in), optional   :: instance
    integer                                :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Satellite_Orbit_Index(thisNode)
       Tree_Node_Satellite_Time_of_Merging_Preset=thisNode%components(thisIndex)%instance(1)%data(mergeTimeIndex)
    else
       Tree_Node_Satellite_Time_of_Merging_Preset=-1.0d0 ! Negative time indicates that this is not a satellite.
    end if
    return
  end function Tree_Node_Satellite_Time_of_Merging_Preset

  subroutine Tree_Node_Satellite_Time_of_Merging_Set_Preset(thisNode,mergeTime,instance)
    !% Set the time of satellite merging.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: mergeTime
    integer,        intent(in), optional   :: instance
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Satellite_Orbit_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%data(mergeTimeIndex)=mergeTime

    return
  end subroutine Tree_Node_Satellite_Time_of_Merging_Set_Preset

  double precision function Tree_Node_Bound_Mass_Preset(thisNode,instance)
    !% Return the satellite bound mass at the current time.
    use Histories
    use FGSL
    use Numerical_Interpolation
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    integer,        intent(in), optional   :: instance
    integer                                :: iTime
    logical                                :: interpolationReset
    type(history)                          :: thisHistory
    type(fgsl_interp_accel)                :: interpolationAccelerator

    ! Check if the component exists.
    if (thisNode%componentExists(componentIndex)) then
       ! It does, so get the preset bound mass history.
       thisHistory=Tree_Node_Bound_Mass_History_Preset(thisNode)
       ! Check if the history exists.
       if (thisHistory%exists()) then
          ! It does, so find the preset time closest to the present one and return the mass at that time.
          interpolationReset=.true.
          iTime=Interpolate_Locate(size(thisHistory%time),thisHistory%time,interpolationAccelerator,Tree_Node_Time(thisNode)&
               &,reset=interpolationReset,closest=.true.)
          Tree_Node_Bound_Mass_Preset=thisHistory%data(iTime,1)
          call Interpolate_Done(interpolationAccelerator=interpolationAccelerator,reset=interpolationReset)  
       else
          ! It does not, so return zero.
          Tree_Node_Bound_Mass_Preset=0.0d0
       end if
       call thisHistory%destroy()
    else
       ! Component does not exist, so return zero.
       Tree_Node_Bound_Mass_Preset=0.0d0
    end if
    return
  end function Tree_Node_Bound_Mass_Preset

  function Tree_Node_Bound_Mass_History_Preset(thisNode)
    !% Return the satellite bound mass history.
    use Histories
    implicit none
    type(history)                          :: Tree_Node_Bound_Mass_History_Preset
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    ! Check if the component exists.
    if (thisNode%componentExists(componentIndex)) then
       ! It does, so get the index of the component.
       thisIndex=Tree_Node_Satellite_Orbit_Index(thisNode)
       ! Check if the component has histories
       if (allocated(thisNode%components(thisIndex)%instance(1)%histories)) then
          ! It does, so return that history.
          Tree_Node_Bound_Mass_History_Preset=thisNode%components(thisIndex)%instance(1)%histories(boundMassIndex)
       else
          ! It does not, so return a null history.
          Tree_Node_Bound_Mass_History_Preset=nullHistory
       end if
    else
       ! Component does not exist, so return a null history.
       Tree_Node_Bound_Mass_History_Preset=nullHistory
    end if
    return
  end function Tree_Node_Bound_Mass_History_Preset

  subroutine Tree_Node_Bound_Mass_History_Set_Preset(thisNode,thisHistory)
    !% Set the node disk stellar properties history.
    use Memory_Management
    use Histories
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    type(history),           intent(in)    :: thisHistory
    integer                                :: thisIndex

    ! Get the index of this component.
    thisIndex=Tree_Node_Satellite_Orbit_Index(thisNode)
    ! If this component does not yet have a history associated with it, then create one.
    if (.not.allocated(thisNode%components(thisIndex)%instance(1)%histories)) then
       allocate(thisNode%components(thisIndex)%instance(1)%histories(historyCount))
       call Memory_Usage_Record(sizeof(thisNode%components(thisIndex)%instance(1)%histories),memoryType=memoryTypeNodes)
    end if
    ! Destroy the current history.
    call thisNode%components(thisIndex)%instance(1)%histories(boundMassIndex)%destroy()
    ! Assign the new history.
    thisNode%components(thisIndex)%instance(1)%histories(boundMassIndex)=thisHistory
    return
  end subroutine Tree_Node_Bound_Mass_History_Set_Preset

  subroutine Tree_Node_Satellite_Virial_Orbit_Set_Preset(thisNode,thisOrbit)
    !% Set the orbit of the satellite at the virial radius.
    use Kepler_Orbits_Structure
    implicit none
    type(keplerOrbit),         intent(inout) :: thisOrbit
    type(treeNode),   pointer, intent(inout) :: thisNode
    integer                                  :: thisIndex

    if (.not.thisNode%componentExists(componentIndex)) call Satellite_Orbit_Create_Preset(thisNode)
    thisIndex=Tree_Node_Satellite_Orbit_Index(thisNode)
    ! Ensure the orbit is defined.
    call thisOrbit%assertIsDefined()
    ! Extract components of the orbit.
    thisNode%components(thisIndex)%instance(1)%data(hostMassIndex          )=thisOrbit%hostMass          ()
    thisNode%components(thisIndex)%instance(1)%data(virialRadiusIndex      )=thisOrbit%radius            ()
    thisNode%components(thisIndex)%instance(1)%data(velocityRadialIndex    )=thisOrbit%velocityRadial    ()
    thisNode%components(thisIndex)%instance(1)%data(velocityTangentialIndex)=thisOrbit%velocityTangential()
    return
  end subroutine Tree_Node_Satellite_Virial_Orbit_Set_Preset

  function Tree_Node_Satellite_Virial_Orbit_Preset(thisNode) result (thisOrbit)
    !% Return the orbit of the satellite at the virial radius.
    use Kepler_Orbits_Structure
    implicit none
    type(keplerOrbit)                        :: thisOrbit
    type(treeNode),   pointer, intent(inout) :: thisNode
    integer                                  :: thisIndex

    if (thisNode%isSatellite().or..not.thisNode%isPrimaryProgenitor()) then
       if (.not.thisNode%componentExists(componentIndex)) call Satellite_Orbit_Create_Preset(thisNode)
       call thisOrbit%reset()
       thisIndex=Tree_Node_Satellite_Orbit_Index(thisNode)
       if (thisNode%components(thisIndex)%instance(1)%data(hostMassIndex) > 0.0d0) then
          call thisOrbit%massesSet            (                                                          &
               &                               Tree_Node_Mass(thisNode)                                , &
               &                               thisNode%components(thisIndex)%instance(1)%data(hostMassIndex          )  &
               &                              )
          call thisOrbit%radiusSet            (thisNode%components(thisIndex)%instance(1)%data(virialRadiusIndex      ))
          call thisOrbit%velocityRadialSet    (thisNode%components(thisIndex)%instance(1)%data(velocityRadialIndex    ))
          call thisOrbit%velocityTangentialSet(thisNode%components(thisIndex)%instance(1)%data(velocityTangentialIndex))
       end if
    else
       call thisOrbit%reset()
    end if
    return
  end function Tree_Node_Satellite_Virial_Orbit_Preset

  !# <nodeMergerTask>
  !#  <unitName>Satellite_Orbit_Create_Preset</unitName>
  !# </nodeMergerTask>
  !# <satelliteHostChangeTask>
  !#  <unitName>Satellite_Orbit_Create_Preset</unitName>
  !# </satelliteHostChangeTask>
  subroutine Satellite_Orbit_Create_Preset(thisNode)
    !% Create a satellite orbit component and assign a time until merging.
    implicit none
    type(treeNode),  pointer, intent(inout) :: thisNode

    if (methodSelected) then
       ! Create a satellite orbit component and assign a time until merging.
       call thisNode%createComponent(componentIndex,propertyCount,dataCount,historyCountInitial)
    end if
    return
  end subroutine Satellite_Orbit_Create_Preset

  !# <nodePromotionTask>
  !#  <unitName>Tree_Node_Satellite_Preset_Promote</unitName>
  !# </nodePromotionTask>
  subroutine Tree_Node_Satellite_Preset_Promote(thisNode)
    !% Ensure that {\tt thisNode} is ready for promotion to its parent. In this case, we simply copy any preset satellite orbit
    !% from the parent.
    use Kepler_Orbits_Structure
    use Histories
    implicit none
    type(treeNode   ),   pointer, intent(inout) :: thisNode
    type(treeNode   ),   pointer                :: parentNode
    type(keplerOrbit)                           :: thisOrbit
    type(history    )                           :: thisHistory

    ! Check if this method is selected.
    if (methodSelected) then
       ! Get the parent node of this node.
       parentNode => thisNode%parentNode
       ! If the parent node has a hot halo component, then add it to that of this node, and perform other changes needed prior to
       ! promotion.
       if (parentNode%componentExists(componentIndex)) then
          thisOrbit=Tree_Node_Satellite_Virial_Orbit(parentNode)
          thisHistory=Tree_Node_Bound_Mass_History  (parentNode)
          call                            Tree_Node_Satellite_Time_Of_Merging_Set(thisNode,Tree_Node_Satellite_Time_Of_Merging(parentNode))
          call                            Tree_Node_Bound_Mass_History_Set       (thisNode,thisHistory                                    )
          if (thisOrbit%isDefined()) call Tree_Node_Satellite_Virial_Orbit_Set   (thisNode,thisOrbit                                      )
          call thisHistory%destroy()
       end if
    end if
    return
  end subroutine Tree_Node_Satellite_Preset_Promote
  
  integer function Tree_Node_Satellite_Orbit_Index(thisNode)
    !% Ensure the satellite orbit component exists and return its position in the components array.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    
    ! Create the component.
    call thisNode%createComponent(componentIndex,propertyCount,dataCount,historyCount)
    Tree_Node_Satellite_Orbit_Index=thisNode%componentIndex(componentIndex)
    return
  end function Tree_Node_Satellite_Orbit_Index

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Satellite_Orbit_Preset_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Satellite_Orbit_Preset</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Satellite_Orbit_Preset_Names(integerProperty,integerPropertyNames,integerPropertyComments&
       &,integerPropertyUnitsSI,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set names of satellite orbit properties to be written to the \glc\ output file.
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
       !@   <name>timeToMerge</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Time until satellite merges.</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='timeToMerge'
       doublePropertyComments(doubleProperty)='Time until satellite merges.'
       doublePropertyUnitsSI (doubleProperty)=gigaYear
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>nodeBoundMass</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Bound mass of the node..</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='nodeBoundMass'
       doublePropertyComments(doubleProperty)='Bound mass of the node..'
       doublePropertyUnitsSI (doubleProperty)=massSolar
    end if
    return
  end subroutine Galacticus_Output_Tree_Satellite_Orbit_Preset_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Satellite_Orbit_Preset_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Satellite_Orbit_Preset</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Satellite_Orbit_Preset_Property_Count(integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of satellite orbit properties to be written to the the \glc\ output file.
    implicit none
    double precision, intent(in)    :: time
    integer,          intent(inout) :: integerPropertyCount,doublePropertyCount

    if (methodSelected) doublePropertyCount=doublePropertyCount+2
    return
  end subroutine Galacticus_Output_Tree_Satellite_Orbit_Preset_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Satellite_Orbit_Preset</unitName>
  !#  <sortName>Galacticus_Output_Tree_Satellite_Orbit_Preset</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Satellite_Orbit_Preset(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store satellite orbit properties in the \glc\ output file buffers.
    use Kind_Numbers
    implicit none
    double precision,        intent(in)             :: time
    type(treeNode),          intent(inout), pointer :: thisNode
    integer,                 intent(inout)          :: integerProperty,integerBufferCount,doubleProperty,doubleBufferCount
    integer(kind=kind_int8), intent(inout)          :: integerBuffer(:,:)
    double precision,        intent(inout)          :: doubleBuffer(:,:)

    if (methodSelected) then
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Satellite_Merge_Time(thisNode)
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Bound_Mass_Preset(thisNode)
    end if
    return
  end subroutine Galacticus_Output_Tree_Satellite_Orbit_Preset

  !# <nodeDumpTask>
  !#  <unitName>Tree_Node_Methods_Satellite_Orbit_Preset_Dump</unitName>
  !# </nodeDumpTask>
  subroutine Tree_Node_Methods_Satellite_Orbit_Preset_Dump(thisNode)
    !% Dump all properties of {\tt thisNode} to screen.
    use Histories
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode
    integer                                :: iTime
    type(history)                          :: thisHistory

    if (methodSelected) then
       if (thisNode%componentExists(componentIndex)) then
          write (0,'(1x,a)'           ) 'satellite orbit component -> properties:'
          write (0,'(2x,a50,1x,e12.6)') 'time until merging:',Tree_Node_Satellite_Merge_Time(thisNode)
          thisHistory=Tree_Node_Bound_Mass_History_Preset(thisNode)
          if (thisHistory%exists()) then
             write (0,'(2x,a50,1x,e12.6)') '        bound mass: [Time] [Mass]'
             do iTime=1,size(thisHistory%time)
                write (0,*) thisHistory%time(iTime),thisHistory%data(iTime,:)
             end do
             call thisHistory%destroy()
          end if
       else
          write (0,'(1x,a)'           ) 'satellite orbit component -> nonexistant'
       end if
    end if
    return
  end subroutine Tree_Node_Methods_Satellite_Orbit_Preset_Dump

end module Tree_Node_Methods_Satellite_Orbit_Preset
