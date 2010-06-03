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






!% Contains a module of spin tree node methods.

module Tree_Node_Methods_Spin_Random
  !% Implement random spin tree node method.
  use Tree_Nodes
  use Tree_Node_Methods
  use Components
  private
  public :: Tree_Node_Methods_Spin_Random_Initialize, Tree_Node_Methods_Spin_Random_Initialize_Spin,&
       & Galacticus_Output_Tree_Spin_Random, Galacticus_Output_Tree_Spin_Random_Property_Count,&
       & Galacticus_Output_Tree_Spin_Random_Names, Tree_Node_Methods_Spin_Random_Dump
  
  ! The index used as a reference for this component.
  integer :: componentIndex=-1

  ! Property indices.
  integer, parameter :: propertyCount=0, dataCount=1, historyCount=0
  integer, parameter :: spinIndex    =1

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
  !#  <unitName>Tree_Node_Methods_Spin_Random_Initialize</unitName>
  !#  <optionName default="random">treeNodeMethodSpin</optionName>
  !# </treeNodeCreateInitialize>
  subroutine Tree_Node_Methods_Spin_Random_Initialize(componentOption,componentTypeCount)
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
    if (componentOption.eq.'random') then
       ! Record that method is selected.
       methodSelected=.true.

       ! Increment the component count and store the value for later reference.
       componentTypeCount=componentTypeCount+1
       componentIndex=componentTypeCount

       ! Display message.
       message='Random spin method selected [component index '
       message=message//componentIndex//']'
       call Galacticus_Display_Message(message,verbosityInfo)

       ! Set up procedure pointers.
       Tree_Node_Spin                         => Tree_Node_Spin_Random
       Tree_Node_Spin_Set                      => null()
       Tree_Node_Spin_Rate_Adjust              => null()
       Tree_Node_Spin_Rate_Compute             => Tree_Node_Rate_Rate_Compute_Dummy
       Tree_Node_Spin_Growth_Rate              => Tree_Node_Spin_Growth_Rate_Random
       Tree_Node_Spin_Growth_Rate_Set          => null()
       Tree_Node_Spin_Growth_Rate_Rate_Adjust  => null()
       Tree_Node_Spin_Growth_Rate_Rate_Compute => Tree_Node_Rate_Rate_Compute_Dummy
    end if
    return
  end subroutine Tree_Node_Methods_Spin_Random_Initialize
  
  
  double precision function Tree_Node_Spin_Random(thisNode)
    !% Return the node spin mass.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex
    
    ! Ensure the spin has been initialized.
    call Tree_Node_Methods_Spin_Random_Initialize_Spin(thisNode)
    thisIndex=Tree_Node_Spin_Random_Index(thisNode)
    Tree_Node_Spin_Random=thisNode%components(thisIndex)%data(spinIndex)
    return
  end function Tree_Node_Spin_Random

  double precision function Tree_Node_Spin_Growth_Rate_Random(thisNode)
    !% Return the node spin growth rate (always zero in this case).
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode

    Tree_Node_Spin_Growth_Rate_Random=0.0d0
    return
  end function Tree_Node_Spin_Growth_Rate_Random


  !# <mergerTreeInitializeTask>
  !#  <unitName>Tree_Node_Methods_Spin_Random_Initialize_Spin</unitName>
  !# </mergerTreeInitializeTask>
  subroutine Tree_Node_Methods_Spin_Random_Initialize_Spin(thisNode)
    !% Initialize the spin of {\tt thisNode}.
    use Cosmological_Parameters
    use Halo_Spin_Distributions
    implicit none
    type(treeNode),  pointer, intent(inout) :: thisNode
    type(treeNode),  pointer                :: relatedNode
    integer                                 :: thisIndex,relatedIndex

    if (methodSelected) then
       if (.not.thisNode%componentExists(componentIndex)) then
          thisIndex=Tree_Node_Spin_Random_Index(thisNode)
          ! Set the spin of the halo.
          thisNode%components(thisIndex)%data(spinIndex)=Halo_Spin_Distribution_Sample(thisNode)
          ! Propagate to any primary children of the halo.
          relatedNode => thisNode%childNode
          do while (associated(relatedNode))
             relatedIndex=Tree_Node_Spin_Random_Index(relatedNode)
             relatedNode%components(relatedIndex)%data(spinIndex)=thisNode%components(thisIndex)%data(spinIndex)
             relatedNode => relatedNode%childNode
          end do
          ! Propagate to any parents of the halo.
          relatedNode => thisNode
          do while (relatedNode%isPrimaryProgenitor())
             relatedNode => relatedNode%parentNode
             relatedIndex=Tree_Node_Spin_Random_Index(relatedNode)
             relatedNode%components(relatedIndex)%data(spinIndex)=thisNode%components(thisIndex)%data(spinIndex)
          end do
       end if
    end if
    return
  end subroutine Tree_Node_Methods_Spin_Random_Initialize_Spin
  

  integer function Tree_Node_Spin_Random_Index(thisNode)
    !% Ensure the spin component exists and return its position in the components array.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    
    call thisNode%createComponent(componentIndex,propertyCount,dataCount,historyCount)
    Tree_Node_Spin_Random_Index=thisNode%componentIndex(componentIndex)
    return
  end function Tree_Node_Spin_Random_Index

  
  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Spin_Random_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Spin_Random</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Spin_Random_Names(integerProperty,integerPropertyNames,integerPropertyComments,doubleProperty&
       &,doublePropertyNames,doublePropertyComments,time)
    !% Set names of spin properties to be written to the \glc\ output file.
    implicit none
    double precision, intent(in)                  :: time
    integer,          intent(inout)               :: integerProperty,doubleProperty
    character(len=*), intent(inout), dimension(:) :: integerPropertyNames,integerPropertyComments,doublePropertyNames &
         &,doublePropertyComments
    
    if (methodSelected) then
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='nodeSpin'
       doublePropertyComments(doubleProperty)='Spin parameter of the node.'
    end if
    return
  end subroutine Galacticus_Output_Tree_Spin_Random_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Spin_Random_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Spin_Random</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Spin_Random_Property_Count(integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of spin properties to be written to the the \glc\ output file.
    implicit none
    double precision, intent(in)    :: time
    integer,          intent(inout) :: integerPropertyCount,doublePropertyCount

    if (methodSelected) doublePropertyCount=doublePropertyCount+1
    return
  end subroutine Galacticus_Output_Tree_Spin_Random_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Spin_Random</unitName>
  !#  <sortName>Galacticus_Output_Tree_Spin_Random</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Spin_Random(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store spin properties in the \glc\ output file buffers.
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
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Spin_Random(thisNode)
    end if
    return
  end subroutine Galacticus_Output_Tree_Spin_Random

  !# <nodeDumpTask>
  !#  <unitName>Tree_Node_Methods_Spin_Random_Dump</unitName>
  !# </nodeDumpTask>
  subroutine Tree_Node_Methods_Spin_Random_Dump(thisNode)
    !% Dump all properties of {\tt thisNode} to screen.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode
 
    if (methodSelected) then
       if (thisNode%componentExists(componentIndex)) then
          write (0,'(1x,a)'           ) 'spin component -> properties:'
          write (0,'(2x,a50,1x,e12.6)') 'halo spin:',Tree_Node_Spin_Random(thisNode)
       else
          write (0,'(1x,a)'           ) 'spin component -> nonexistant'
       end if
    end if
    return
  end subroutine Tree_Node_Methods_Spin_Random_Dump

end module Tree_Node_Methods_Spin_Random
