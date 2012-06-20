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


!% Contains a module which implements writing a merger tree to the \glc\ output file.

module Galacticus_Output_Merger_Tree
  !% Implements writing a merger tree to the \glc\ output file.
  use IO_HDF5
  use ISO_Varying_String
  use Merger_Trees
  use Galacticus_Error
  use Galacticus_HDF5
  use Kind_Numbers
  implicit none
  private
  public :: Galacticus_Merger_Tree_Output, Galacticus_Merger_Tree_Output_Finalize

  ! Output groups.
  type outputGroup
     !% Type used for output group information.
     logical                 :: opened,integerAttributesWritten,doubleAttributesWritten
     type(hdf5Object)        :: hdf5Group,nodeDataGroup
     integer(kind=kind_int8) :: length
  end type outputGroup
  type(hdf5Object)                                             :: outputsGroup
  logical                                                      :: outputsGroupOpened=.false.
  type(outputGroup),               dimension(:),   allocatable :: outputGroups
  integer                                                      :: outputGroupsCount=0
  integer,                         parameter                   :: outputGroupsIncrement=10

  ! Number of properties of each type.
  integer                                                      :: integerPropertyCount=-1, doublePropertyCount=-1
  !$omp threadprivate(integerPropertyCount,doublePropertyCount)

  ! Buffers for properties.
  integer                                                      :: integerPropertiesWritten=0, doublePropertiesWritten=0
  integer                                                      :: integerBufferCount      =0, doubleBufferCount      =0
  integer,                         parameter                   :: bufferSize=1024, nameLengthMax=256, commentLengthMax=256
  integer(kind=kind_int8),         dimension(:,:), allocatable :: integerBuffer
  double precision,                dimension(:,:), allocatable :: doubleBuffer
  character(len=nameLengthMax),    dimension(:),   allocatable :: integerPropertyNames   ,doublePropertyNames
  character(len=commentLengthMax), dimension(:),   allocatable :: integerPropertyComments,doublePropertyComments
  double precision,                dimension(:),   allocatable :: integerPropertyUnitsSI ,doublePropertyUnitsSI
  !$omp threadprivate(integerPropertiesWritten,doublePropertiesWritten,integerBufferCount,doubleBufferCount,integerBuffer)
  !$omp threadprivate(doubleBuffer,integerPropertyNames,doublePropertyNames,integerPropertyUnitsSI,doublePropertyUnitsSI)
  !$omp threadprivate(integerPropertyComments,doublePropertyComments)

  ! Flag indicating if module is initialized.
  logical                                                      :: mergerTreeOutputInitialized=.false.

  ! Flag indicating if merger tree references are to be output.
  logical                                                      :: mergerTreeOutputReferences

contains

  subroutine Galacticus_Merger_Tree_Output(thisTree,iOutput,time,isLastOutput)
    !% Write properties of nodes in {\tt thisTree} to the \glc\ output file.
    use Tree_Nodes
    use Galacticus_Output_Open
    use Galacticus_Merger_Tree_Output_Filters
    use Input_Parameters
    !# <include directive="mergerTreeOutputTask" type="moduleUse">
    include 'galacticus.output.merger_tree.tasks.modules.inc'
    !# </include>
    !# <include directive="mergerTreeExtraOutputTask" type="moduleUse">
    include 'galacticus.output.merger_tree.tasks.extra.modules.inc'
    !# </include>
    implicit none
    type(mergerTree),      intent(inout)          :: thisTree
    integer,               intent(in)             :: iOutput
    double precision,      intent(in)             :: time
    logical,               intent(in),   optional :: isLastOutput
    type(treeNode),        pointer                :: thisNode
    integer(kind=HSIZE_T), dimension(1)           :: referenceStart,referenceLength
    integer                                       :: integerProperty,doubleProperty,iProperty,iGroup
    logical                                       :: nodePassesFilter
    type(hdf5Object)                              :: toDataset

    ! Initialize if necessary.
    !$omp critical(Merger_Tree_Output)
    if (.not.mergerTreeOutputInitialized) then
       
       ! Ensure file is open.
       call Galacticus_Output_Open_File

       ! Get input parameters.
       !@ <inputParameter>
       !@   <name>mergerTreeOutputReferences</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Specifies whether or not references to individual merger tree datasets should be output.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@   <group>output</group>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeOutputReferences',mergerTreeOutputReferences,defaultValue=.false.)

       ! Flag that the module is now initialized.
       mergerTreeOutputInitialized=.true.

    end if
    
    ! Create an output group.
    call Make_Output_Group(iOutput,time)

    ! Ensure that the output filter subsystem is initialized.
    call Galacticus_Merger_Tree_Output_Filter_Initialize()

    ! Count up the number of properties to be output.
    call Count_Properties(time)

    ! Ensure buffers are allocated for temporary property storage.
    call Allocate_Buffers

    ! Get names for all properties to be output.
    call Establish_Property_Names(time)

    ! Loop over all nodes in the tree.
    integerPropertiesWritten=0
    doublePropertiesWritten =0
    thisNode => thisTree%baseNode
    do while (associated(thisNode))
       if (Tree_Node_Time(thisNode) == time) then
          nodePassesFilter=Galacticus_Merger_Tree_Output_Filter(thisNode)
          if (nodePassesFilter) then
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
             if (integerBufferCount == bufferSize) call Integer_Buffer_Dump(iOutput)
             if (doubleBufferCount  == bufferSize) call Double_Buffer_Dump (iOutput)
             
          end if

          ! Do any extra output tasks.
          !# <include directive="mergerTreeExtraOutputTask" type="code" action="subroutine">
          !#  <subroutineArgs>thisNode,iOutput,thisTree%index,nodePassesFilter</subroutineArgs>
          include 'galacticus.output.merger_tree.tasks.extra.inc'
          !# </include>

       end if
       call thisNode%walkTreeWithSatellites(thisNode)
    end do
    if (integerPropertyCount > 0 .and. integerBufferCount > 0) call Integer_Buffer_Dump(iOutput)
    if (doublePropertyCount  > 0 .and. doubleBufferCount  > 0) call Double_Buffer_Dump (iOutput)

    ! Compute the start and length of regions to reference.
    !$omp critical(HDF5_Access)
    referenceLength(1)=max(integerPropertiesWritten,doublePropertiesWritten)
    referenceStart (1)=outputGroups(iOutput)%length

    ! Create references to the datasets if requested.
    if (mergerTreeOutputReferences) then
       
       ! Ensure that a group has been made for this merger tree.
       call Galacticus_Merger_Tree_Output_Make_Group(thisTree,iOutput)
       
       ! Create references for this tree.
       if (integerPropertyCount > 0 .and. integerPropertiesWritten > 0) then
          do iProperty=1,integerPropertyCount
             toDataset=outputGroups(iOutput)%nodeDataGroup%openDataset(integerPropertyNames(iProperty))
             call thisTree%hdf5Group%createReference1D(toDataset,integerPropertyNames(iProperty),referenceStart+1,referenceLength)
             call toDataset%close()
          end do
       end if
       if (doublePropertyCount > 0  .and. doublePropertiesWritten  > 0) then
          do iProperty=1,doublePropertyCount
             toDataset=outputGroups(iOutput)%nodeDataGroup%openDataset(doublePropertyNames(iProperty))
             call thisTree%hdf5Group%createReference1D(toDataset,doublePropertyNames(iProperty),referenceStart+1,referenceLength)
             call toDataset%close()
          end do
       end if
       
       ! Close the tree group.
       call thisTree%hdf5Group%close()

    end if
    
    ! Store the start position and length of the node data for this tree, along with its volume weight.
    call outputGroups(iOutput)%hdf5Group%writeDataset([thisTree%index]       ,"mergerTreeIndex"     ,"Index of each merger tree."                                  ,appendTo=.true.)
    call outputGroups(iOutput)%hdf5Group%writeDataset(referenceStart         ,"mergerTreeStartIndex","Index in nodeData datasets at which each merger tree begins.",appendTo=.true.)
    call outputGroups(iOutput)%hdf5Group%writeDataset(referenceLength        ,"mergerTreeCount"     ,"Number of nodes in nodeData datasets for each merger tree."  ,appendTo=.true.)
    call outputGroups(iOutput)%hdf5Group%writeDataset([thisTree%volumeWeight],"mergerTreeWeight"    ,"Number density of each tree [Mpc⁻³]."                        ,appendTo=.true.)

    ! Increment the number of nodes written to this output group.
    outputGroups(iOutput)%length=outputGroups(iOutput)%length+referenceLength(1)
    !$omp end critical(HDF5_Access)

    ! Close down if this is the final output.
    if (present(isLastOutput)) then
       if (isLastOutput) then
          ! Close any open output groups.
          do iGroup=1,size(outputGroups)
             if (outputGroups(iGroup)%nodeDataGroup%isOpen()) call outputGroups(iGroup)%nodeDataGroup%close()
             if (outputGroups(iGroup)%hdf5Group    %isOpen()) call outputGroups(iGroup)%hdf5Group    %close()
          end do
          if (outputsGroup%isOpen()) call outputsGroup%close()
          ! Close the file.
          call Galacticus_Output_Close_File
       end if
    end if
    !$omp end critical(Merger_Tree_Output)
    return
  end subroutine Galacticus_Merger_Tree_Output

  subroutine Galacticus_Merger_Tree_Output_Finalize()
    !% Finalize merger tree output by closing any open groups.
    implicit none
    integer :: iGroup
    
    ! Close any open output groups.
    do iGroup=1,size(outputGroups)
       if (outputGroups(iGroup)%nodeDataGroup%isOpen()) call outputGroups(iGroup)%nodeDataGroup%close()
       if (outputGroups(iGroup)%hdf5Group    %isOpen()) call outputGroups(iGroup)%hdf5Group    %close()
    end do
    if (outputsGroup%isOpen()) call outputsGroup%close()
    return
  end subroutine Galacticus_Merger_Tree_Output_Finalize

  subroutine Galacticus_Merger_Tree_Output_Make_Group(thisTree,iOutput)
    !% Make an group in the \glc\ file in which to store {\tt thisTree}.
    use Numerical_Constants_Astronomical
    use String_Handling
    implicit none
    type(mergerTree),      intent(inout) :: thisTree
    integer,               intent(in)    :: iOutput
    type(varying_string)                 :: groupName,commentText

    ! Create a name for the group.
    groupName='mergerTree'
    groupName=groupName//thisTree%index
    
    ! Create a comment for the group.
    commentText='Data for nodes within merger tree ID='
    commentText=commentText//thisTree%index
    
    ! Create a group for the tree.
    thisTree%hdf5Group=outputGroups(iOutput)%hdf5Group%openGroup(char(groupName),char(commentText))
    
    ! Add the merger tree weight to the group.
    call thisTree%hdf5Group%writeAttribute(thisTree%volumeWeight,"volumeWeight"         )
    call thisTree%hdf5Group%writeAttribute(1.0d0/megaParsec**3  ,"volumeWeightUnitsInSI")
    return
  end subroutine Galacticus_Merger_Tree_Output_Make_Group

  subroutine Integer_Buffer_Dump(iOutput)
    !% Dump the contents of the integer properties buffer to the \glc\ output file.
    implicit none
    integer,          intent(in) :: iOutput
    integer                      :: iProperty
    type(hdf5Object)             :: thisDataset

    ! Write integer data from the buffer.
    if (integerPropertyCount > 0) then
       !$omp critical(HDF5_Access)
       do iProperty=1,integerPropertyCount
          call outputGroups(iOutput)%nodeDataGroup%writeDataset(integerBuffer(1:integerBufferCount,iProperty),integerPropertyNames(iProperty) &
               &,integerPropertyComments(iProperty),appendTo=.true.)
          if (.not.outputGroups(iOutput)%integerAttributesWritten.and.integerPropertyUnitsSI(iProperty) /= 0.0d0) then
             thisDataset=outputGroups(iOutput)%nodeDataGroup%openDataset(integerPropertyNames(iProperty))
             call thisDataset%writeAttribute(integerPropertyUnitsSI(iProperty),"unitsInSI")
             call thisDataset%close()
          end if
       end do
       integerPropertiesWritten=integerPropertiesWritten+integerBufferCount
       integerBufferCount=0
       outputGroups(iOutput)%integerAttributesWritten=.true.
       !$omp end critical(HDF5_Access)
    end if
    return
  end subroutine Integer_Buffer_Dump

  subroutine Double_Buffer_Dump(iOutput)
    !% Dump the contents of the double precision properties buffer to the \glc\ output file.
    implicit none
    integer,          intent(in) :: iOutput
    integer                      :: iProperty
    type(hdf5Object)             :: thisDataset

    ! Write double data from the buffer.
    if (doublePropertyCount > 0) then
       !$omp critical(HDF5_Access)
       do iProperty=1,doublePropertyCount
          call outputGroups(iOutput)%nodeDataGroup%writeDataset(doubleBuffer(1:doubleBufferCount,iProperty),doublePropertyNames(iProperty) &
               &,doublePropertyComments(iProperty),appendTo=.true.)
          if (.not.outputGroups(iOutput)%doubleAttributesWritten.and.doublePropertyUnitsSI(iProperty) /= 0.0d0) then
             thisDataset=outputGroups(iOutput)%nodeDataGroup%openDataset(doublePropertyNames(iProperty))
             call thisDataset%writeAttribute(doublePropertyUnitsSI(iProperty),"unitsInSI")
             call thisDataset%close()
          end if
       end do
       doublePropertiesWritten=doublePropertiesWritten+doubleBufferCount
       doubleBufferCount=0
       outputGroups(iOutput)%doubleAttributesWritten=.true.
       !$omp end critical(HDF5_Access)
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

    if (integerPropertyCount > 0 .and. (.not.allocated(integerBuffer) .or. integerPropertyCount > size(integerPropertyNames)) ) then
       if (allocated(integerBuffer)) then
          call Dealloc_Array(integerBuffer          )
          call Dealloc_Array(integerPropertyNames   )
          call Dealloc_Array(integerPropertyComments)
          call Dealloc_Array(integerPropertyUnitsSI )
       end if
       call Alloc_Array(integerBuffer          ,[bufferSize,integerPropertyCount])
       call Alloc_Array(integerPropertyNames              ,[integerPropertyCount])
       call Alloc_Array(integerPropertyComments           ,[integerPropertyCount])
       call Alloc_Array(integerPropertyUnitsSI            ,[integerPropertyCount])
    end if
    if (doublePropertyCount  > 0 .and. (.not.allocated(doubleBuffer ) .or. doublePropertyCount  > size(doublePropertyNames ))) then
       if (allocated(doubleBuffer )) then
          call Dealloc_Array(doubleBuffer           )
          call Dealloc_Array(doublePropertyNames    )
          call Dealloc_Array(doublePropertyComments )
          call Dealloc_Array(doublePropertyUnitsSI  )
       end if
       call Alloc_Array(doubleBuffer           ,[bufferSize,doublePropertyCount])
       call Alloc_Array(doublePropertyNames               ,[doublePropertyCount])
       call Alloc_Array(doublePropertyComments            ,[doublePropertyCount])
       call Alloc_Array(doublePropertyUnitsSI             ,[doublePropertyCount])
    end if
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
    !#  <subroutineArgs>integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time</subroutineArgs>
    include 'galacticus.output.merger_tree.names.inc'
    !# </include>

    return
  end subroutine Establish_Property_Names

  subroutine Make_Output_Group(iOutput,time)
    !% Create a group in which to store this output.
    use String_Handling
    use Cosmology_Functions
    use Memory_Management
    use Numerical_Constants_Astronomical
    !# <include directive="outputGroupOutputTask" type="moduleUse">
    include 'galacticus.output.merger_tree.outputGroup.tasks.modules.inc'
    !# </include>
    implicit none
    integer,              intent(in)                :: iOutput      
    double precision,     intent(in)                :: time
    type(outputGroup),    dimension(:), allocatable :: outputGroupsTemporary
    type(varying_string)                            :: groupName,commentText

    !$omp critical (HDF5_Access)
    ! Ensure group ID space is large enough.
    if (iOutput > outputGroupsCount) then
       if (allocated(outputGroups)) then
          call Move_Alloc(outputGroups,outputGroupsTemporary)
          outputGroupsCount=outputGroupsCount+outputGroupsIncrement
          allocate(outputGroups(outputGroupsCount))
          outputGroups(1:size(outputGroupsTemporary))=outputGroupsTemporary
          outputGroups(size(outputGroupsTemporary)+1:size(outputGroups))%opened                  =.false.
          outputGroups(size(outputGroupsTemporary)+1:size(outputGroups))%integerAttributesWritten=.false.
          outputGroups(size(outputGroupsTemporary)+1:size(outputGroups))%doubleAttributesWritten =.false.
          outputGroups(size(outputGroupsTemporary)+1:size(outputGroups))%length                  =0
          deallocate(outputGroupsTemporary)
          call Memory_Usage_Record(sizeof(outputGroups(1)),blockCount=0)
       else
          outputGroupsCount=outputGroupsIncrement
          allocate(outputGroups(outputGroupsCount))
          outputGroups%opened                  =.false.
          outputGroups%integerAttributesWritten=.false.
          outputGroups%doubleAttributesWritten =.false.
          outputGroups%length=0
          call Memory_Usage_Record(sizeof(outputGroups))
       end if
    end if

    ! Make the enclosing group if it has not been created.
    if (.not.outputsGroupOpened) then
       outputsGroup=galacticusOutputFile%openGroup('Outputs','Contains all outputs from Galacticus.')
       outputsGroupOpened=.true.
    end if

    ! Create the group if it has not been created.
    if (.not.outputGroups(iOutput)%opened) then
       ! Create a name for the group.
       groupName='Output'
       groupName=groupName//iOutput

       ! Create a comment for the group.
       commentText='Data for output number '
       commentText=commentText//iOutput

       ! Create a group for the tree.
       !@ <outputType>
       !@   <name>nodeData</name>
       !@   <description>A representation of the state of all nodes in the simulation at a given time. It consists of numerous datasets which gives the properties of nodes in all merger trees at that time.</description>
       !@ </outputType>
       outputGroups(iOutput)%hdf5Group    =outputsGroup                   %openGroup(char(groupName),char(commentText))
       outputGroups(iOutput)%nodeDataGroup=outputGroups(iOutput)%hdf5Group%openGroup("nodeData","Group containing data on all nodes at this output.")
       outputGroups(iOutput)%opened                  =.true.
       outputGroups(iOutput)%integerAttributesWritten=.false.
       outputGroups(iOutput)%doubleAttributesWritten =.false.

       ! Add the time to this group.
       call outputGroups(iOutput)%hdf5Group%writeAttribute(time                  ,'outputTime'           )
       call outputGroups(iOutput)%hdf5Group%writeAttribute(gigaYear              ,'timeUnitInSI'         )
       call outputGroups(iOutput)%hdf5Group%writeAttribute(Expansion_Factor(time),'outputExpansionFactor')

       ! Establish all other properties.
       !# <include directive="outputGroupOutputTask" type="code" action="subroutine">
       !#  <subroutineArgs>outputGroups(iOutput)%hdf5Group,time</subroutineArgs>
       include 'galacticus.output.merger_tree.outputGroup.tasks.inc'
       !# </include>

    end if
    !$omp end critical(HDF5_Access)
    return
  end subroutine Make_Output_Group

end module Galacticus_Output_Merger_Tree
