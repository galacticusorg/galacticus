!% Contains a module which implements writing a merger tree to the \glc\ output file.

module Galacticus_Output_Merger_Tree
  !% Implements writing a merger tree to the \glc\ output file.
  use HDF5
  use IO_HDF5
  use ISO_Varying_String
  use Merger_Trees
  use Galacticus_Error
  use Galacticus_HDF5
  use Galacticus_HDF5_Groups
  private
  public :: Galacticus_Merger_Tree_Output

  ! Output groups.
  integer                                        :: outputGroupID
  integer                                        :: outputGroupsCount=0
  integer,             parameter                 :: outputGroupsIncrement=10
  integer(kind=HID_T), dimension(:), allocatable :: outputGroupsIDs

  ! Number of properties of each type.
  integer                                        :: integerPropertyCount=-1, doublePropertyCount=-1
  !$omp threadprivate(integerPropertyCount,doublePropertyCount)

  ! Buffers for properties.
  integer                                                      :: integerPropertiesWritten=0, doublePropertiesWritten=0
  integer                                                      :: integerBufferCount=0, doubleBufferCount=0
  integer,                         parameter                   :: bufferSize=1024, nameLengthMax=256, commentLengthMax=256
  integer,                         dimension(:,:), allocatable :: integerBuffer
  double precision,                dimension(:,:), allocatable :: doubleBuffer
  integer,                         dimension(:),   allocatable :: integerPropertyIDs,doublePropertyIDs
  character(len=nameLengthMax),    dimension(:),   allocatable :: integerPropertyNames,doublePropertyNames
  character(len=commentLengthMax), dimension(:),   allocatable :: integerPropertyComments,doublePropertyComments
  !$omp threadprivate(integerPropertiesWritten,doublePropertiesWritten,integerBufferCount,doubleBufferCount,integerBuffer)
  !$omp threadprivate(doubleBuffer,integerPropertyIDs,doublePropertyIDs,integerPropertyNames,doublePropertyNames)
  !$omp threadprivate(integerPropertyComments,doublePropertyComments)

contains

  subroutine Galacticus_Merger_Tree_Output(thisTree,iOutput,time,isLastOutput)
    !% Write properties of nodes in {\tt thisTree} to the \glc\ output file.
    use Tree_Nodes
    use Tree_Node_Methods
    use Galacticus_Output_Open
    !# <include directive="mergerTreeOutputTask" type="moduleUse">
    include 'galacticus.output.merger_tree.tasks.modules.inc'
    !# </include>
    implicit none
    type(mergerTree), intent(inout)          :: thisTree
    integer,          intent(in)             :: iOutput
    double precision, intent(in)             :: time
    logical,          intent(in),   optional :: isLastOutput
    type(treeNode),   pointer                :: thisNode
    integer                                  :: integerProperty,doubleProperty

    ! Ensure file is open.
    call Galacticus_Output_Open_File

    ! Create an output group.
    call Make_Output_Group(iOutput,time)

    ! Count up the number of properties to be output.
    call Count_Properties(time)

    ! Ensure buffers are allocated for temporary property storage.
    call Allocate_Buffers

    ! Get names for all properties to be output.
    call Establish_Property_Names(time)

    ! Ensure that a group has been made for this merger tree.
    call Galacticus_Merger_Tree_Output_Make_Group(thisTree,iOutput)

    ! Loop over all nodes in the tree.
    integerPropertiesWritten=0
    doublePropertiesWritten=0
    thisNode => thisTree%baseNode
    do while (associated(thisNode))
       if (Tree_Node_Time(thisNode) == time) then

          ! Establish node link properties.
          if (integerPropertyCount > 0) then
             integerProperty=0
             integerBufferCount=integerBufferCount+1     
          end if
          if (doublePropertyCount > 0) then
             doubleProperty=0
             doubleBufferCount=doubleBufferCount+1
          end if
          
          ! Establish all other properties.
          !# <include directive="mergerTreeOutputTask" type="code" action="subroutine">
          !#  <subroutineArgs>thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty,doubleBufferCount,doubleBuffer,time</subroutineArgs>
          include 'galacticus.output.merger_tree.tasks.inc'
          !# </include>
          
          ! If buffer is full, dump it to file.
          if (integerBufferCount == bufferSize) call Integer_Buffer_Dump(thisTree)
          if (doubleBufferCount  == bufferSize) call Double_Buffer_Dump (thisTree)
       end if
       call thisNode%walkTreeWithSatellites()
    end do
    if (integerPropertyCount > 0 .and. integerBufferCount > 0) call Integer_Buffer_Dump(thisTree)
    if (doublePropertyCount  > 0 .and. doubleBufferCount  > 0) call Double_Buffer_Dump (thisTree)
      
    ! Close down if this is the final output.
    if (present(isLastOutput)) then
       if (isLastOutput) call Galacticus_Output_Close_File
    end if
    return
  end subroutine Galacticus_Merger_Tree_Output

  subroutine Galacticus_Merger_Tree_Output_Make_Group(thisTree,iOutput)
    !% Make an group in the \glc\ file in which to store {\tt thisTree}.
    use String_Handling
    implicit none
    type(mergerTree),      intent(inout) :: thisTree
    integer,               intent(in)    :: iOutput
    type(varying_string)                 :: groupName,commentText
    integer(kind=HID_T)                  :: datasetID
    double precision                     :: datasetValue(1)

    ! Reset the group ID for this tree.
    thisTree%hdf5GroupID=0
    
    ! Create a name for the group.
    groupName='mergerTree'
    groupName=groupName//thisTree%index
    
    ! Create a comment for the group.
    commentText='Data for nodes within merger tree ID='
    commentText=commentText//thisTree%index
    
    ! Create a group for the tree.
    thisTree%hdf5GroupID=Galacticus_Output_Make_Group(groupName,commentText,outputGroupsIDs(iOutput))
    
    ! Add the merger tree weight to the group.
    datasetID=0
    datasetValue=[thisTree%volumeWeight]
    call Galacticus_Output_Dataset(thisTree%hdf5GroupID,datasetID,'volumeWeight','Number of trees per unit volume.'&
         &,datasetValue)
    return
  end subroutine Galacticus_Merger_Tree_Output_Make_Group

  subroutine Integer_Buffer_Dump(thisTree)
    !% Dump the contents of the integer properties buffer to the \glc\ output file.
    use Galacticus_HDF5_Groups
    implicit none
    type(mergerTree), intent(inout) :: thisTree
    integer                         :: iProperty

    ! Write integer data from the buffer.
    if (integerPropertyCount > 0) then
       do iProperty=1,integerPropertyCount
          call Galacticus_Output_Dataset(thisTree%hdf5GroupID,integerPropertyIDs(iProperty),integerPropertyNames(iProperty)&
               &,integerPropertyComments(iProperty),integerBuffer(1:integerBufferCount,iProperty),isExtendable=.true.)
       end do
       integerPropertiesWritten=integerPropertiesWritten+integerBufferCount
       integerBufferCount=0
    end if
    return
  end subroutine Integer_Buffer_Dump

  subroutine Double_Buffer_Dump(thisTree)
    !% Dump the contents of the double precision properties buffer to the \glc\ output file.
    use Galacticus_HDF5_Groups
    implicit none
    type(mergerTree),      intent(inout) :: thisTree
    integer                              :: iProperty

    ! Write double data from the buffer.
    if (doublePropertyCount > 0) then
       do iProperty=1,doublePropertyCount
          call Galacticus_Output_Dataset(thisTree%hdf5GroupID,doublePropertyIDs(iProperty),doublePropertyNames(iProperty)&
               &,doublePropertyComments(iProperty),doubleBuffer(1:doubleBufferCount,iProperty),isExtendable=.true.)
       end do
       doublePropertiesWritten=doublePropertiesWritten+doubleBufferCount
       doubleBufferCount=0
    end if
    return
  end subroutine Double_Buffer_Dump

  subroutine Count_Properties(time)
    !% Count up the number of properties that will be output.
    !# <include directive="mergerTreeOutputPropertyCount" type="moduleUse">
    include 'galacticus.output.merger_tree.property_count.modules.inc'
    !# </include>
    implicit none
    double precision, intent(in) :: time

    integerPropertyCount=0
    doublePropertyCount=0
    !# <include directive="mergerTreeOutputPropertyCount" type="code" action="subroutine">
    !#  <subroutineArgs>integerPropertyCount,doublePropertyCount,time</subroutineArgs>
    include 'galacticus.output.merger_tree.property_count.inc'
    !# </include>
    return
  end subroutine Count_Properties
  
  subroutine Allocate_Buffers
    !% Allocate buffers for storage of properties.
    use Memory_Management
    implicit none

    if (integerPropertyCount > 0 .and. (.not.allocated(integerBuffer) .or. integerPropertyCount > size(integerPropertyIDs)) ) then
       if (allocated(integerBuffer)) then
          call Dealloc_Array(integerBuffer          )
          call Dealloc_Array(integerPropertyIDs     )
          call Dealloc_Array(integerPropertyNames   )
          call Dealloc_Array(integerPropertyComments)
       end if
       call Alloc_Array(integerBuffer          ,bufferSize,integerPropertyCount,'integerBuffer'          )
       call Alloc_Array(integerPropertyIDs                ,integerPropertyCount,'integerPropertyIDs'     )
       call Alloc_Array(integerPropertyNames              ,integerPropertyCount,'integerPropertyNames'   )
       call Alloc_Array(integerPropertyComments           ,integerPropertyCount,'integerPropertyComments')
    end if
    integerPropertyIDs=0
    if (doublePropertyCount  > 0 .and. (.not.allocated(doubleBuffer ) .or. doublePropertyCount  > size(doublePropertyIDs ))) then
       if (allocated(doubleBuffer )) then
          call Dealloc_Array(doubleBuffer           )
          call Dealloc_Array(doublePropertyIDs      )
          call Dealloc_Array(doublePropertyNames    )
          call Dealloc_Array(doublePropertyComments )
       end if
       call Alloc_Array(doubleBuffer           ,bufferSize,doublePropertyCount ,'doubleBuffer'           )
       call Alloc_Array(doublePropertyIDs                 ,doublePropertyCount ,'doublePropertyIDs'      )
       call Alloc_Array(doublePropertyNames               ,doublePropertyCount ,'doublePropertyNames'    )
       call Alloc_Array(doublePropertyComments            ,doublePropertyCount ,'doublePropertyComments' )
    end if
    doublePropertyIDs=0
    return
  end subroutine Allocate_Buffers

  subroutine Establish_Property_Names(time)
    !% Set names for the properties.
    !# <include directive="mergerTreeOutputNames" type="moduleUse">
    include 'galacticus.output.merger_tree.names.modules.inc'
    !# </include>
    implicit none
    double precision, intent(in) :: time
    integer                      :: integerProperty,doubleProperty

    integerProperty=0
    doubleProperty =0
    !# <include directive="mergerTreeOutputNames" type="code" action="subroutine">
    !#  <subroutineArgs>integerProperty,integerPropertyNames,integerPropertyComments,doubleProperty,doublePropertyNames,doublePropertyComments,time</subroutineArgs>
    include 'galacticus.output.merger_tree.names.inc'
    !# </include>

    return
  end subroutine Establish_Property_Names

  subroutine Make_Output_Group(iOutput,time)
    !% Create a group in which to store this output.
    use String_Handling
    use Cosmology_Functions
    use Tree_Node_Methods
    implicit none
    integer,             intent(in)                :: iOutput      
    double precision,    intent(in)                :: time
    integer(kind=HID_T), dimension(:), allocatable :: outputGroupsIDsTemporary
    type(varying_string)                           :: groupName,commentText
    integer(kind=HID_T)                            :: datasetID
    double precision                               :: datasetValue(1)

    !$omp critical (Outputs_Group_Create)
    ! Ensure group ID space is large enough.
    if (iOutput > outputGroupsCount) then
       if (allocated(outputGroupsIDs)) then
          call Move_Alloc(outputGroupsIDs,outputGroupsIDsTemporary)
          outputGroupsCount=outputGroupsIncrement
          allocate(outputGroupsIDs(outputGroupsCount))
          outputGroupsIDs(1:size(outputGroupsIDsTemporary))=outputGroupsIDsTemporary
          outputGroupsIDs(size(outputGroupsIDsTemporary)+1:size(outputGroupsIDs))=0
          deallocate(outputGroupsIDsTemporary)
       else
          outputGroupsCount=outputGroupsIncrement
          allocate(outputGroupsIDs(outputGroupsCount))
          outputGroupsIDs=0
       end if
    end if

    ! Make the enclosing group if it has not been created.
    if (outputGroupID <= 0) then

       ! Create a name for the group.
       groupName='Outputs'

       ! Create a comment for the group.
       commentText='Contains all outputs from Galacticus.'

       ! Create a group for the tree.
       outputGroupID=Galacticus_Output_Make_Group(groupName,commentText,galacticusOutputID)
    end if

    ! Create the group if it has not been created.
    if (outputGroupsIDs(iOutput) <= 0) then
 
       ! Create a name for the group.
       groupName='Output'
       groupName=groupName//iOutput

       ! Create a comment for the group.
       commentText='Data for output number '
       commentText=commentText//iOutput

       ! Create a group for the tree.
       outputGroupsIDs(iOutput)=Galacticus_Output_Make_Group(groupName,commentText,outputGroupID)

       ! Add the time to this group.
       datasetID=0
       datasetValue=[time]
       call Galacticus_Output_Dataset(outputGroupsIDs(iOutput),datasetID,'outputTime','Time for this output.'&
            &,datasetValue)
       datasetID=0
       datasetValue=[Expansion_Factor(time)]
       call Galacticus_Output_Dataset(outputGroupsIDs(iOutput),datasetID,'outputExpansionFactor','Expansion factor for this&
            & output.',datasetValue)
    end if
    !$omp end critical (Outputs_Group_Create)
    return
  end subroutine Make_Output_Group

end module Galacticus_Output_Merger_Tree
