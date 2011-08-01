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


!% Contains a module which implements storing and restoring of the complete internal state of a merger tree.

module Merger_Trees_State_Store
  !% Implements storing and restoring of the complete internal state of a merger tree. Useful primarily for debugging purposes to
  !% begin running a tree from just prior to the point of failure.
  use Kind_Numbers
  use ISO_Varying_String
  implicit none
  private
  public :: Merger_Tree_State_Store, Merger_Tree_State_Store_Initialize

  ! Flag indicating if our tree has already been returned.
  logical :: doneTree=.false.

  ! File name from which stored trees should be read.
  type(varying_string) :: mergerTreeStateStoreFile

contains

  !# <mergerTreeConstructMethod>
  !#  <unitName>Merger_Tree_State_Store_Initialize</unitName>
  !# </mergerTreeConstructMethod>
  subroutine Merger_Tree_State_Store_Initialize(mergerTreeConstructMethod,Merger_Tree_Construct)
    !% Initialize the ``state restore'' method for constructing merger trees.
    use Input_Parameters
    implicit none
    type(varying_string),          intent(in)    :: mergerTreeConstructMethod
    procedure(),          pointer, intent(inout) :: Merger_Tree_Construct

    ! Check if our method is to be used.
    if (mergerTreeConstructMethod == 'state restore') then
       ! Assign pointer to our merger tree construction subroutine.
       Merger_Tree_Construct => Merger_Tree_State_Restore
       ! Get the name of the file from which to read the stored merger tree state.
       !@ <inputParameter>
       !@   <name>mergerTreeStateStoreFile</name>
       !@   <defaultValue>storedTree.dat</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of a file from which to restore a merger tree state.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeStateStoreFile',mergerTreeStateStoreFile,defaultValue='storedTree.dat')
    end if
    return
  end subroutine Merger_Tree_State_Store_Initialize

  subroutine Merger_Tree_State_Store(thisTree,storeFile)
    !% Store the complete internal state of a merger tree to file.
    use Merger_Trees
    use Tree_Nodes
    use Components
    use Galacticus_State
    use Galacticus_Error
    use File_Utilities
    implicit none
    type(mergerTree),        intent(in)                :: thisTree
    character(len=*),        intent(in)                :: storeFile
    type(treeNode),          pointer                   :: thisNode
    integer(kind=kind_int8), allocatable, dimension(:) :: nodeIndices
    integer                                            :: iComponent,iInstance,iHistory,nodeCount,fileUnit
    
    ! Take a snapshot of the internal state and store it.
    call Galacticus_State_Snapshot
    call Galacticus_State_Store

    ! Count nodes in the tree.
    nodeCount=0
    thisNode => thisTree%baseNode
    do while (associated(thisNode))
       nodeCount=nodeCount+1
       call thisNode%walkTreeWithSatellites(thisNode)
    end do

    ! Allocate and populate an array of node indices in the order in which the tree will be traversed.
    allocate(nodeIndices(nodeCount))
    thisNode => thisTree%baseNode
    nodeCount=0
    do while (associated(thisNode))
       nodeCount=nodeCount+1
       nodeIndices(nodeCount)=thisNode%index()
       call thisNode%walkTreeWithSatellites(thisNode)
    end do

    ! Open an output file.
    fileUnit=File_Units_Get()
    open(fileUnit,file=trim(storeFile),status='unknown',form='unformatted')
    ! Write basic tree information.
    write (fileUnit) thisTree%index,thisTree%volumeWeight,thisTree%initialized,nodeCount,Node_Array_Position(thisTree%baseNode%index(),nodeIndices)
    ! Start at the base of the tree.
    thisNode => thisTree%baseNode
    ! Loop over all nodes.
    do while (associated(thisNode))
       ! Write all node information.
       ! Indices.
       write (fileUnit) thisNode%index(),thisNode%uniqueID()
       ! Pointers to other nodes.
       write (fileUnit) &
            &  Node_Array_Position(thisNode%parentNode   %index(),nodeIndices) &
            & ,Node_Array_Position(thisNode%childNode    %index(),nodeIndices) &
            & ,Node_Array_Position(thisNode%siblingNode  %index(),nodeIndices) &
            & ,Node_Array_Position(thisNode%satelliteNode%index(),nodeIndices) &
            & ,Node_Array_Position(thisNode%mergeNode    %index(),nodeIndices) &
            & ,Node_Array_Position(thisNode%mergeeNode   %index(),nodeIndices) &
            & ,Node_Array_Position(thisNode%nextMergee   %index(),nodeIndices)
       ! Component indices.
       write (fileUnit) allocated(thisNode%componentIndex)
       if (allocated(thisNode%componentIndex)) then
          write (fileUnit) size(thisNode%componentIndex)
          write (fileUnit) thisNode%componentIndex
       end if
       ! Components.
       write (fileUnit) allocated(thisNode%components)
       if (allocated(thisNode%components)) then
          write (fileUnit) size(thisNode%components)
          do iComponent=1,size(thisNode%components)
             ! Instances.
             write (fileUnit) allocated(thisNode%components(iComponent)%instance)
             if (allocated(thisNode%components(iComponent)%instance)) then
                write (fileUnit) size(thisNode%components(iComponent)%instance)
                do iInstance=1,size(thisNode%components(iComponent)%instance)
                   ! Properties.
                   write (fileUnit) allocated(thisNode%components(iComponent)%instance(iInstance)%properties)
                   if (allocated(thisNode%components(iComponent)%instance(iInstance)%properties)) then
                      write (fileUnit) size(thisNode%components(iComponent)%instance(iInstance)%properties,dim=1)
                      write (fileUnit) thisNode%components(iComponent)%instance(iInstance)%properties(:,propertyValue)
                   end if
                   ! Data.
                   write (fileUnit) allocated(thisNode%components(iComponent)%instance(iInstance)%data)
                   if (allocated(thisNode%components(iComponent)%instance(iInstance)%data)) then
                      write (fileUnit) size(thisNode%components(iComponent)%instance(iInstance)%data)
                      write (fileUnit) thisNode%components(iComponent)%instance(iInstance)%data(:)
                   end if
                   ! Histories.
                   write (fileUnit) allocated(thisNode%components(iComponent)%instance(iInstance)%histories)
                   if (allocated(thisNode%components(iComponent)%instance(iInstance)%histories)) then
                      write (fileUnit) size(thisNode%components(iComponent)%instance(iInstance)%histories)
                      do iHistory=1,size(thisNode%components(iComponent)%instance(iInstance)%histories)
                         write (fileUnit) allocated(thisNode%components(iComponent)%instance(iInstance)%histories(iHistory)%time)
                         if (allocated(thisNode%components(iComponent)%instance(iInstance)%histories(iHistory)%time)) then
                            write (fileUnit) size(thisNode%components(iComponent)%instance(iInstance)%histories(iHistory)%time)&
                                 &,size(thisNode%components(iComponent)%instance(iInstance)%histories(iHistory)%data,dim=2)
                            write (fileUnit) thisNode%components(iComponent)%instance(iInstance)%histories(iHistory)%time
                            write (fileUnit) thisNode%components(iComponent)%instance(iInstance)%histories(iHistory)%data
                            write (fileUnit) thisNode%components(iComponent)%instance(iInstance)%histories(iHistory)%rangeType
                         end if
                      end do
                   end if
                end do
             end if
          end do
       end if
       ! Move to the next node in the tree.
       call thisNode%walkTreeWithSatellites(thisNode)
    end do
    close(fileUnit)
    ! Destroy the temporary array of indices.
    deallocate(nodeIndices)
    return
  end subroutine Merger_Tree_State_Store

  integer function Node_Array_Position(nodeIndex,nodeIndices)
    !% Returns the position of a node in the output list given its index.
    implicit none
    integer(kind=kind_int8), intent(in)               :: nodeIndex
    integer(kind=kind_int8), intent(in), dimension(:) :: nodeIndices
    
    if (nodeIndex == -1) then
       Node_Array_Position=-1
    else
       Node_Array_Position=1
       do while (nodeIndices(Node_Array_Position) /= nodeIndex)
          Node_Array_Position=Node_Array_Position+1
       end do
    end if
    return
  end function Node_Array_Position

  subroutine Merger_Tree_State_Restore(thisTree,skipTree)
    !% Restores the state of a merger tree from file.
    use Merger_Trees
    use Tree_Nodes
    use Components
    use Galacticus_State
    use File_Utilities
    use Galacticus_Error
    implicit none
    type(mergerTree),   intent(inout)             :: thisTree
    logical,            intent(in)                :: skipTree
    type(treeNodeList), allocatable, dimension(:) :: nodes
    integer                                       :: iComponent,iInstance,iHistory,nodeCount,nodeArrayIndex,iNode,parentNodeIndex&
         & ,childNodeIndex,siblingNodeIndex,satelliteNodeIndex,mergeNodeIndex,mergeeNodeIndex,nextMergeeIndex,arraySize&
         &,arraySize2 ,fileUnit
    integer(kind=kind_int8)                       :: nodeIndex,nodeUniqueID
    logical                                       :: isAllocated

    ! If the tree is to be skipped, or if we've already returned the tree, simply return.
    if (skipTree.or.doneTree) return

    ! Retrieve stored internal state if possible.
    call Galacticus_State_Retrieve

    ! Open an output file.
    fileUnit=File_Units_Get()
    open(unit=fileUnit,file=char(mergerTreeStateStoreFile),status='old',form='unformatted')
    ! Read basic tree information.
    read (fileUnit) thisTree%index,thisTree%volumeWeight,thisTree%initialized,nodeCount,nodeArrayIndex
    ! Allocate a list of nodes.
    allocate(nodes(nodeCount))
    ! Create nodes.
     do iNode=1,nodeCount
        call thisTree%createNode(nodes(iNode)%node)
    end do
    ! Assign the tree base node.
    thisTree%baseNode => nodes(nodeArrayIndex)%node
    ! Loop over all nodes.
    do iNode=1,nodeCount
       ! Read all node information.
       ! Indices/
       read (fileUnit) nodeIndex,nodeUniqueID
       call nodes(iNode)%node%indexSet(nodeIndex)
       call nodes(iNode)%node%uniqueIDSet(nodeUniqueID)
       ! Pointers to other nodes.
       read (fileUnit) parentNodeIndex,childNodeIndex,siblingNodeIndex,satelliteNodeIndex,mergeNodeIndex,mergeeNodeIndex,nextMergeeIndex
       nodes(iNode)%node%parentNode    => Pointed_At_Node(parentNodeIndex   ,nodes)
       nodes(iNode)%node%childNode     => Pointed_At_Node(childNodeIndex    ,nodes)
       nodes(iNode)%node%siblingNode   => Pointed_At_Node(siblingNodeIndex  ,nodes)
       nodes(iNode)%node%satelliteNode => Pointed_At_Node(satelliteNodeIndex,nodes)
       nodes(iNode)%node%mergeNode     => Pointed_At_Node(mergeNodeIndex    ,nodes)
       nodes(iNode)%node%mergeeNode    => Pointed_At_Node(mergeeNodeIndex   ,nodes)
       nodes(iNode)%node%nextMergee    => Pointed_At_Node(nextMergeeIndex   ,nodes)
       ! Component index array.
       read (fileUnit) isAllocated
       if (isAllocated) then
          read (fileUnit) arraySize
          if (arraySize /= size(nodes(iNode)%node%componentIndex)) call Galacticus_Error_Report('Merger_Tree_State_Restore','number of components has changed')
          read (fileUnit) nodes(iNode)%node%componentIndex
       end if
       ! Components.
       read (fileUnit) isAllocated
       if (isAllocated) then
          read (fileUnit) arraySize
          allocate(nodes(iNode)%node%components(arraySize))
          do iComponent=1,size(nodes(iNode)%node%components)        
             ! Instances.
             read (fileUnit) isAllocated
             if (isAllocated) then
                read (fileUnit) arraySize
                allocate(nodes(iNode)%node%components(iComponent)%instance(arraySize))
                do iInstance=1,size(nodes(iNode)%node%components(iComponent)%instance)
                   ! Properties.
                   read (fileUnit) isAllocated
                   if (isAllocated) then
                      read (fileUnit) arraySize
                      allocate(nodes(iNode)%node%components(iComponent)%instance(iInstance)%properties(arraySize,2))
                      read (fileUnit) nodes(iNode)%node%components(iComponent)%instance(iInstance)%properties(:,propertyValue)
                   end if
                   ! Data.
                   read (fileUnit) isAllocated
                   if (isAllocated) then
                      read (fileUnit) arraySize
                      allocate(nodes(iNode)%node%components(iComponent)%instance(iInstance)%data(arraySize))
                      read (fileUnit) nodes(iNode)%node%components(iComponent)%instance(iInstance)%data(:)
                   end if
                   ! Histories.
                   read (fileUnit) isAllocated
                   if (isAllocated) then
                      read (fileUnit) arraySize
                      allocate(nodes(iNode)%node%components(iComponent)%instance(iInstance)%histories(arraySize))
                      do iHistory=1,size(nodes(iNode)%node%components(iComponent)%instance(iInstance)%histories)
                         read (fileUnit) isAllocated
                         if (isAllocated) then
                            read (fileUnit) arraySize,arraySize2
                            allocate(nodes(iNode)%node%components(iComponent)%instance(iInstance)%histories(iHistory)%time(arraySize))
                            allocate(nodes(iNode)%node%components(iComponent)%instance(iInstance)%histories(iHistory)%data(arraySize,arraySize2))
                            read (fileUnit) nodes(iNode)%node%components(iComponent)%instance(iInstance)%histories(iHistory)%time
                            read (fileUnit) nodes(iNode)%node%components(iComponent)%instance(iInstance)%histories(iHistory)%data
                            read (fileUnit) nodes(iNode)%node%components(iComponent)%instance(iInstance)%histories(iHistory)%rangeType
                         end if
                      end do
                   end if
                end do
             end if
          end do
       end if
    end do
    close(fileUnit)
    ! Destroy the list of nodes.
    deallocate(nodes)
    ! Flag that the tree has been read and returned.
    doneTree=.true.
    return
  end subroutine Merger_Tree_State_Restore

  function Pointed_At_Node(nodeArrayIndex,nodes)
    !% Return a pointer to a node, given its position in the array of nodes. Return a null pointer if the array index is $-1$.
    use Tree_Nodes
    type(treeNode),     pointer                  :: Pointed_At_Node
    integer,            intent(in)               :: nodeArrayIndex
    type(treeNodeList), intent(in), dimension(:) :: nodes

    if (nodeArrayIndex == -1) then
       Pointed_At_Node => null()
    else
       Pointed_At_Node => nodes(nodeArrayIndex)%node
    end if
    return
  end function Pointed_At_Node
  
end module Merger_Trees_State_Store
