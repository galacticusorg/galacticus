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

!% Contains a module of satellite orbit tree node methods.

module Tree_Node_Methods_Satellite_Orbit_Very_Simple
  !% Implement satellite orbit tree node methods.
  use Tree_Nodes
  use Kepler_Orbits_Structure
  use Components
  implicit none
  private
  public :: Tree_Node_Methods_Satellite_Orbit_Very_Simple_Initialize, Satellite_Orbit_Create_Very_Simple,&
       & Galacticus_Output_Tree_Satellite_Orbit_Very_Simple, Galacticus_Output_Tree_Satellite_Orbit_VSimple_Property_Count,&
       & Galacticus_Output_Tree_Satellite_Orbit_Very_Simple_Names, Tree_Node_Methods_Satellite_Orbit_Very_Simple_Dump,&
       & Satellite_Orbit_Very_Simple_Formation_Task, Satellite_Orbit_Very_Simple_Property_Identifiers_Decode
  
  ! The index used as a reference for this component.
  integer :: componentIndex=-1

  ! Property indices.
  integer, parameter :: propertyCount=0, dataCount=1, historyCount=0
  integer, parameter :: mergeTimeIndex=1

  ! Define procedure pointers.
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Satellite_Merge_Time</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Satellite_Time_Of_Merging</methodName>
  !# </treeNodeMethodsPointer>

  ! Flag to indicate if this method is selected.
  logical :: methodSelected=.false.

  ! Flag indicating whether or not to reset satellite orbits on halo formation events.
  logical :: satelliteOrbitResetOnHaloFormation

  ! Option controlling whether or not unbound virial orbits are acceptable.
  logical, parameter :: acceptUnboundOrbits=.false.

contains

  !# <treeNodeCreateInitialize>
  !#  <unitName>Tree_Node_Methods_Satellite_Orbit_Very_Simple_Initialize</unitName>
  !#  <optionName>treeNodeMethodSatelliteOrbit</optionName>
  !# </treeNodeCreateInitialize>
  subroutine Tree_Node_Methods_Satellite_Orbit_Very_Simple_Initialize(componentOption,componentTypeCount)
    !% Initializes the tree node satellite orbit methods module.
    use ISO_Varying_String
    use Input_Parameters
    use Galacticus_Error
    use Galacticus_Display
    use String_Handling
    implicit none
    type(varying_string), intent(in)    :: componentOption
    integer,              intent(inout) :: componentTypeCount
    type(varying_string)                :: satelliteMergingMethod,message

    ! Check if this implementation is selected.
    if (componentOption == 'verySimple') then
       ! Record that method is selected.
       methodSelected=.true.

       ! Increment the component count and store the value for later reference.
       componentTypeCount=componentTypeCount+1
       componentIndex=componentTypeCount

       ! Display message.
       message='Very simple satellite orbit method selected [component index '
       message=message//componentIndex//']'
       call Galacticus_Display_Message(message,verbosityInfo)

       ! Set up procedure pointers.
       Tree_Node_Satellite_Merge_Time                   => Tree_Node_Satellite_Merge_Time_Very_Simple
       Tree_Node_Satellite_Merge_Time_Set               => Tree_Node_Satellite_Merge_Time_Set_Very_Simple
       Tree_Node_Satellite_Merge_Time_Rate_Adjust       => null()
       Tree_Node_Satellite_Merge_Time_Rate_Compute      => Tree_Node_Rate_Rate_Compute_Dummy

       Tree_Node_Satellite_Time_of_Merging              => Tree_Node_Satellite_Time_of_Merging_Very_Simple
       Tree_Node_Satellite_Time_of_Merging_Set          => null()
       Tree_Node_Satellite_Time_of_Merging_Rate_Adjust  => null()
       Tree_Node_Satellite_Time_of_Merging_Rate_Compute => Tree_Node_Rate_Rate_Compute_Dummy

       ! Determine if satellite orbits are to be reset on halo formation events.
       !@ <inputParameter>
       !@   <name>satelliteOrbitResetOnHaloFormation</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether satellite virial orbital parameters should be reset on halo formation events.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('satelliteOrbitResetOnHaloFormation',satelliteOrbitResetOnHaloFormation,defaultValue=.false.)
    end if

    return
  end subroutine Tree_Node_Methods_Satellite_Orbit_Very_Simple_Initialize
  
  double precision function Tree_Node_Satellite_Merge_Time_Very_Simple(thisNode,instance)
    !% Return the time until satellite merging.
    implicit none
    type(treeNode), intent(inout), pointer  :: thisNode
    integer,        intent(in),    optional :: instance
    integer                                 :: thisIndex

    if (thisNode%componentExists(componentIndex).and.thisNode%isSatellite()) then
       thisIndex=Tree_Node_Satellite_Orbit_Index(thisNode)
       Tree_Node_Satellite_Merge_Time_Very_Simple=thisNode%components(thisIndex)%instance(1)%data(mergeTimeIndex)-Tree_Node_Time(thisNode)
    else
       Tree_Node_Satellite_Merge_Time_Very_Simple=-1.0d0 ! Negative time indicates that this is not a satellite.
    end if
    return
  end function Tree_Node_Satellite_Merge_Time_Very_Simple

  double precision function Tree_Node_Satellite_Time_Of_Merging_Very_Simple(thisNode,instance)
    !% Return the time until satellite merging.
    implicit none
    type(treeNode), intent(inout), pointer  :: thisNode
    integer,        intent(in),    optional :: instance
    integer                                 :: thisIndex

    if (thisNode%componentExists(componentIndex).and.thisNode%isSatellite()) then
       thisIndex=Tree_Node_Satellite_Orbit_Index(thisNode)
       Tree_Node_Satellite_Time_Of_Merging_Very_Simple=thisNode%components(thisIndex)%instance(1)%data(mergeTimeIndex)
    else
       Tree_Node_Satellite_Time_Of_Merging_Very_Simple=-1.0d0 ! Negative time indicates that this is not a satellite.
    end if
    return
  end function Tree_Node_Satellite_Time_Of_Merging_Very_Simple

  subroutine Tree_Node_Satellite_Merge_Time_Set_Very_Simple(thisNode,mergeTime,instance)
    !% Set the time until satellite merging.
    implicit none
    type(treeNode),   intent(inout), pointer  :: thisNode
    double precision, intent(in)              :: mergeTime
    integer,          intent(in),    optional :: instance
    integer                                   :: thisIndex

    thisIndex=Tree_Node_Satellite_Orbit_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%data(mergeTimeIndex)=mergeTime+Tree_Node_Time(thisNode)
    return
  end subroutine Tree_Node_Satellite_Merge_Time_Set_Very_Simple

  !# <haloFormationTask>
  !#  <unitName>Satellite_Orbit_Very_Simple_Formation_Task</unitName>
  !# </haloFormationTask>
  subroutine Satellite_Orbit_Very_Simple_Formation_Task(thisNode)
    !% Reset the orbits of satellite galaxies on halo formation events.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    type(treeNode), pointer                :: satelliteNode

    ! Return immediately if orbits are not to be reset.
    if (.not.satelliteOrbitResetOnHaloFormation) return

    ! Loop over all satellites.
    satelliteNode => thisNode%satelliteNode
    do while (associated(satelliteNode))
       ! Create a new orbit for this satellite.
       call Satellite_Orbit_Create_Very_Simple(satelliteNode)
       satelliteNode => satelliteNode%siblingNode
    end do

    return
  end subroutine Satellite_Orbit_Very_Simple_Formation_Task

  !# <nodeMergerTask>
  !#  <unitName>Satellite_Orbit_Create_Very_Simple</unitName>
  !# </nodeMergerTask>
  !# <satelliteHostChangeTask>
  !#  <unitName>Satellite_Orbit_Create_Very_Simple</unitName>
  !# </satelliteHostChangeTask>
  subroutine Satellite_Orbit_Create_Very_Simple(thisNode)
    !% Create a satellite orbit component and assign a time until merging and a bound mass equal initially to the total halo mass.
    use Numerical_Constants_Math
    use Dark_Matter_Halo_Scales
    use Virial_Orbits
    use Satellite_Merging_Timescales
    implicit none
    type(treeNode),     pointer, intent(inout) :: thisNode
    type(treeNode),     pointer                :: hostNode
    logical                                    :: isNewSatellite
    integer                                    :: thisIndex
    double precision                           :: mergeTime
    type(keplerOrbit)                          :: thisOrbit

    if (methodSelected) then
       ! Determine if this is a new satellite.
       isNewSatellite=.not.thisNode%componentExists(componentIndex)

       ! Create a satellite orbit component and assign a time until merging.
       call thisNode%createComponent(componentIndex,propertyCount,dataCount,historyCount)

       ! Get an orbit for this satellite.
       hostNode => thisNode%parentNode
       thisOrbit=Virial_Orbital_Parameters(thisNode,hostNode,acceptUnboundOrbits)

       ! Compute and store a time until merging.
       mergeTime=Satellite_Time_Until_Merging(thisNode,thisOrbit)
       if (mergeTime >= 0.0d0) call Tree_Node_Satellite_Merge_Time_Set_Very_Simple(thisNode,mergeTime)

    end if
    return
  end subroutine Satellite_Orbit_Create_Very_Simple

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
  !#  <unitName>Galacticus_Output_Tree_Satellite_Orbit_Very_Simple_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Satellite_Orbit_Very_Simple</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Satellite_Orbit_Very_Simple_Names(integerProperty,integerPropertyNames,integerPropertyComments&
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
    end if
    return
  end subroutine Galacticus_Output_Tree_Satellite_Orbit_Very_Simple_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Satellite_Orbit_VSimple_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Satellite_Orbit_VSimple</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Satellite_Orbit_VSimple_Property_Count(integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of satellite orbit properties to be written to the the \glc\ output file.
    implicit none
    double precision, intent(in)    :: time
    integer,          intent(inout) :: integerPropertyCount,doublePropertyCount

    if (methodSelected) doublePropertyCount=doublePropertyCount+1
    return
  end subroutine Galacticus_Output_Tree_Satellite_Orbit_VSimple_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Satellite_Orbit_Very_Simple</unitName>
  !#  <sortName>Galacticus_Output_Tree_Satellite_Orbit_Very_Simple</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Satellite_Orbit_Very_Simple(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store satellite orbit properties in the \glc\ output file buffers.
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
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Satellite_Merge_Time(thisNode)
    end if
    return
  end subroutine Galacticus_Output_Tree_Satellite_Orbit_Very_Simple

  !# <nodeDumpTask>
  !#  <unitName>Tree_Node_Methods_Satellite_Orbit_Very_Simple_Dump</unitName>
  !# </nodeDumpTask>
  subroutine Tree_Node_Methods_Satellite_Orbit_Very_Simple_Dump(thisNode)
    !% Dump all properties of {\tt thisNode} to screen.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    if (methodSelected) then
       if (thisNode%componentExists(componentIndex)) then
          write (0,'(1x,a)'           ) 'very simple satellite orbit component -> properties:'
          write (0,'(2x,a50,1x,e12.6)') 'time until merging:',Tree_Node_Satellite_Merge_Time(thisNode)
       else
          write (0,'(1x,a)'           ) 'satellite orbit component -> nonexistant'
       end if
    end if
    return
  end subroutine Tree_Node_Methods_Satellite_Orbit_Very_Simple_Dump

  !# <decodePropertyIdentifiersTask>
  !#  <unitName>Satellite_Orbit_Very_Simple_Property_Identifiers_Decode</unitName>
  !# </decodePropertyIdentifiersTask>
  subroutine Satellite_Orbit_Very_Simple_Property_Identifiers_Decode(propertyComponent,propertyObject,propertyIndex,matchedProperty,propertyName)
    !% Decodes property identifiers to property names for the standard satellite orbit module.
    use ISO_Varying_String
    implicit none
    integer,              intent(in)    :: propertyComponent,propertyObject,propertyIndex
    logical,              intent(inout) :: matchedProperty
    type(varying_string), intent(inout) :: propertyName

    if (methodSelected.and..not.matchedProperty) then
       if (propertyComponent == componentIndex) then
          matchedProperty=.true.
          propertyName="satelliteOrbit:"
          select case (propertyObject)
          case (objectTypeProperty)
             select case (propertyIndex)
             case (mergeTimeIndex)
                propertyName=propertyName//":mergeTime"
             end select
          end select
       end if
    end if

    return
  end subroutine Satellite_Orbit_Very_Simple_Property_Identifiers_Decode
  
end module Tree_Node_Methods_Satellite_Orbit_Very_Simple
