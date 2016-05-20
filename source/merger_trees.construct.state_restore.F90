!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@carnegiescience.edu>
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

!% Contains a module which implements storing and restoring of the complete internal state of a merger tree.

module Merger_Trees_State_Store
  !% Implements storing and restoring of the complete internal state of a merger tree. Useful primarily for debugging purposes to
  !% begin running a tree from just prior to the point of failure.
  use Kind_Numbers
  use ISO_Varying_String
  implicit none
  private
  public :: Merger_Tree_State_Store, Merger_Tree_State_From_File, Merger_Tree_State_Store_Initialize

  ! File name from which stored trees should be read.
  type   (varying_string) :: mergerTreeStateStoreFile

  ! File unit for the tree data file.
  integer                 :: treeDataUnit

contains

  !# <mergerTreeConstructMethod>
  !#  <unitName>Merger_Tree_State_Store_Initialize</unitName>
  !# </mergerTreeConstructMethod>
  subroutine Merger_Tree_State_Store_Initialize(mergerTreeConstructMethod,Merger_Tree_Construct)
    !% Initialize the ``state restore'' method for constructing merger trees.
    use Input_Parameters
    implicit none
    type     (varying_string           ), intent(in   )          :: mergerTreeConstructMethod
    procedure(Merger_Tree_State_Restore), intent(inout), pointer :: Merger_Tree_Construct

    ! Check if our method is to be used.
    if (mergerTreeConstructMethod == 'stateRestore') then
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
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeStateStoreFile',mergerTreeStateStoreFile,defaultValue='storedTree.dat')
       ! Open the tree file.
       open(newunit=treeDataUnit,file=char(mergerTreeStateStoreFile),status='old',form='unformatted')
    end if
    return
  end subroutine Merger_Tree_State_Store_Initialize

  subroutine Merger_Tree_State_Store(tree,storeFile,snapshot)
    !% Store the complete internal state of a merger tree to file.
    use Galacticus_Nodes
    use Galacticus_State
    use Galacticus_Error
    implicit none
    type     (mergerTree    ), intent(in   ), target       :: tree
    character(len=*         ), intent(in   )               :: storeFile
    logical                  , intent(in   ), optional     :: snapshot
    type     (mergerTree    ), pointer                     :: treeCurrent
    type     (treeNode      ), pointer                     :: currentNodeInTree, thisNode
    class    (nodeEvent     ), pointer                     :: event
    integer  (kind=kind_int8), allocatable  , dimension(:) :: nodeIndices
    integer                  , allocatable  , dimension(:) :: nodeCountTree
    type     (varying_string), save                        :: storeFilePrevious
    integer                                                :: fileUnit         , nodeCount, &
         &                                                    treeCount        , iTree    , &
         &                                                    eventCount
    !# <optionalArgument name="snapshot" defaultsTo=".true." />

    ! Take a snapshot of the internal state and store it.
    if (snapshot_) then
       call Galacticus_State_Snapshot()
       call Galacticus_State_Store   ()
    end if
    ! Open an output file. (Append to the old file if the file name has not changed.)
    if (trim(storeFile) == storeFilePrevious) then
       open(newunit=fileUnit,file=trim(storeFile),status='old'    ,form='unformatted',access='append')
    else
       storeFilePrevious=trim(storeFile)
       open(newunit=fileUnit,file=trim(storeFile),status='unknown',form='unformatted'                )
    end if    
    ! Count trees and check for events attached to trees.
    treeCount   =  0
    treeCurrent => tree
    do while (associated(treeCurrent))
       if (associated(treeCurrent%event)) call Galacticus_Error_Report('Merger_Tree_State_Store','tree events not current supported')
       treeCount   =  treeCount           +1
       treeCurrent => treeCurrent%nextTree
    end do
    write (fileUnit) treeCount
    ! Iterate over trees, counting nodes.
    allocate(nodeCountTree(treeCount))
    nodeCountTree =  0
    iTree         =  0
    treeCurrent   => tree
    do while (associated(treeCurrent))
       iTree             =  iTree               +1
       thisNode          => treeCurrent%baseNode
       currentNodeInTree => null()
       do while (associated(thisNode))
          nodeCountTree(iTree)=nodeCountTree(iTree)+1
          call Merger_Tree_State_Walk_Tree(thisNode,currentNodeInTree)
       end do
       treeCurrent => treeCurrent%nextTree
    end do
    nodeCount=sum(nodeCountTree)
    write (fileUnit) nodeCountTree
    ! Allocate and populate an array of node indices in the order in which the tree will be traversed.
    allocate(nodeIndices(nodeCount))
    nodeCount   =  0
    treeCurrent => tree
    do while (associated(treeCurrent))
       thisNode          => treeCurrent%baseNode
       currentNodeInTree => null()
       do while (associated(thisNode))
          nodeCount=nodeCount+1
          nodeIndices(nodeCount)=thisNode%uniqueID()
          call Merger_Tree_State_Walk_Tree(thisNode,currentNodeInTree)
       end do
       treeCurrent => treeCurrent%nextTree
    end do
    ! Write basic tree information.
    treeCurrent => tree
    do while (associated(treeCurrent))
       write (fileUnit) treeCurrent%index,treeCurrent%volumeWeight,treeCurrent%initializedUntil,Node_Array_Position(treeCurrent%baseNode%uniqueID(),nodeIndices)
       treeCurrent => treeCurrent%nextTree
    end do
    ! Output nodes.
    treeCurrent => tree
    do while (associated(treeCurrent))
       ! Start at the base of the tree.
       thisNode          => treeCurrent%baseNode
       currentNodeInTree => null()
       ! Loop over all nodes.
       do while (associated(thisNode))
          ! Write all node information.
          ! Indices.
          write (fileUnit) thisNode%index(),thisNode%uniqueID()
          ! Pointers to other nodes.
          write (fileUnit)                                                            &
               & Node_Array_Position(thisNode%parent        %uniqueID(),nodeIndices), &
               & Node_Array_Position(thisNode%firstChild    %uniqueID(),nodeIndices), &
               & Node_Array_Position(thisNode%sibling       %uniqueID(),nodeIndices), &
               & Node_Array_Position(thisNode%firstSatellite%uniqueID(),nodeIndices), &
               & Node_Array_Position(thisNode%mergeTarget   %uniqueID(),nodeIndices), &
               & Node_Array_Position(thisNode%firstMergee   %uniqueID(),nodeIndices), &
               & Node_Array_Position(thisNode%siblingMergee %uniqueID(),nodeIndices), &
               & Node_Array_Position(thisNode%formationNode %uniqueID(),nodeIndices)
          ! Store the node.
          call thisNode%dumpRaw(fileUnit)
          ! Store any events attached to the node.
          eventCount =  0
          event      => thisNode%event
          do while (associated(event))
             eventCount =  eventCount+1
             event      => event%next
          end do
          write (fileUnit) eventCount
          event => thisNode%event
          do while (associated(event))
             call event%dumpRaw(fileUnit)
             if (associated(event%node)) then
                write (fileUnit) Node_Array_Position(event%node%uniqueID(),nodeIndices)
             else
                write (fileUnit) -1
             end if
             event => event%next
          end do
          ! Move to the next node in the tree.
          call Merger_Tree_State_Walk_Tree(thisNode,currentNodeInTree)
       end do
       ! Move to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do
    close(fileUnit)
    ! Destroy the temporary array of indices.
    deallocate(nodeIndices)
    return
  end subroutine Merger_Tree_State_Store

  integer function Node_Array_Position(nodeIndex,nodeIndices)
    !% Returns the position of a node in the output list given its index.
    use Galacticus_Error
    use ISO_Varying_String
    use String_Handling
    implicit none
    integer(kind=kind_int8)              , intent(in   ) :: nodeIndex
    integer(kind=kind_int8), dimension(:), intent(in   ) :: nodeIndices
    type   (varying_string)                              :: message

    if (nodeIndex == -1) then
       Node_Array_Position=-1
    else
       Node_Array_Position=1
       do while (Node_Array_Position <= size(nodeIndices) .and. nodeIndices(Node_Array_Position) /= nodeIndex)
          Node_Array_Position=Node_Array_Position+1
       end do
       if (nodeIndices(Node_Array_Position) /= nodeIndex) then
          message="node ["
          message=message//nodeIndex//"] could not be found in merger tree"
          call Galacticus_Error_Report('Node_Array_Position',message)          
       end if
    end if
    return
  end function Node_Array_Position

  subroutine Merger_Tree_State_Restore(tree,skipTree)
    !% Restores the state of a merger tree from file.
    use Galacticus_Nodes
    use Galacticus_State
    use Galacticus_Error
    use String_Handling
    implicit none
    type   (mergerTree    ), intent(inout), target  :: tree
    logical                , intent(in   )          :: skipTree
    type   (mergerTree    )               , pointer :: treeCurrent, treeNext

    ! Retrieve stored internal state if possible.
    call Galacticus_State_Retrieve()
    ! Read the tree(s).
    call Merger_Tree_State_From_Unit(tree,skipTree,treeDataUnit)
    ! If the tree is to be skipped, destroy it.
    if (skipTree) then
       treeCurrent => tree
       do while (associated(treeCurrent))
          treeNext    => treeCurrent%nextTree
          call treeCurrent%destroy()
          treeCurrent%baseNode => null()
          treeCurrent => treeNext
       end do
    end if
    return
  end subroutine Merger_Tree_State_Restore

  function Pointed_At_Node(nodeArrayIndex,nodes)
    !% Return a pointer to a node, given its position in the array of nodes. Return a null pointer if the array index is $-1$.
    use Galacticus_Nodes
    type   (treeNode    ), pointer                     :: Pointed_At_Node
    integer                            , intent(in   ) :: nodeArrayIndex
    type   (treeNodeList), dimension(:), intent(in   ) :: nodes

    if (nodeArrayIndex == -1) then
       Pointed_At_Node => null()
    else
       Pointed_At_Node => nodes(nodeArrayIndex)%node
    end if
    return
  end function Pointed_At_Node

  subroutine Merger_Tree_State_Walk_Tree(thisNode,currentNodeInTree)
    !% Walk a merger tree for the purposes of storing the full state to file. Includes walking of formation nodes.
    use Galacticus_Nodes
    implicit none
    type(treeNode), intent(inout), pointer :: currentNodeInTree, thisNode

    if (associated(thisNode%formationNode)) then
       if (.not.associated(currentNodeInTree)) currentNodeInTree => thisNode
       thisNode => thisNode%formationNode
    else
       if (associated(currentNodeInTree)) then
          thisNode => currentNodeInTree
          currentNodeInTree => null()
       end if
       thisNode => thisNode%walkTreeWithSatellites()
    end if
    return
  end subroutine Merger_Tree_State_Walk_Tree

  subroutine Merger_Tree_State_From_File(tree,fileName,deleteAfterRead)
    !% Read the state of a merger tree from file.
    use Galacticus_Nodes
    use Galacticus_Error
    implicit none
    type     (mergerTree), intent(inout)           :: tree
    character(len=*     ), intent(in   )           :: fileName
    logical              , intent(in   ), optional :: deleteAfterRead
    integer                                        :: treeUnit      , ioStatus
    !# <optionalArgument name="deleteAfterRead" defaultsTo=".false." />

    ! Open the file.
    open(newUnit=treeUnit,file=fileName,status='old',form='unformatted',iostat=ioStatus)
    if (ioStatus /= 0) call Galacticus_Error_Report('Merger_Tree_State_From_File','unable to open file "'//trim(fileName)//'"')
    ! Read the tree(s).
    call Merger_Tree_State_From_Unit(tree,skipTree=.false.,unit=treeUnit)
    ! Close the file.
    if (deleteAfterRead_) then
       close(treeUnit,status='delete')
    else
       close(treeUnit                )
    end if
    return
  end subroutine Merger_Tree_State_From_File

  subroutine Merger_Tree_State_From_Unit(tree,skipTree,unit)
    !% Restores the state of a merger tree from file.
    use Galacticus_Nodes
    use Galacticus_Error
    use String_Handling
    implicit none
    type   (mergerTree    ), intent(inout), target       :: tree
    logical                , intent(in   )               :: skipTree
    integer                , intent(in   )               :: unit
    type   (mergerTree    )               , pointer      :: treeCurrent
    type   (treeNodeList  ), allocatable  , dimension(:) :: nodes
    integer                , allocatable  , dimension(:) :: nodeCountTree
    class  (nodeEvent     ), pointer                     :: event              , eventPrevious
    integer                                              :: fileStatus         , firstChildIndex   , firstMergeeIndex   , &
         &                                                  firstSatelliteIndex, formationNodeIndex, iNode              , &
         &                                                  mergeTargetIndex   , nodeArrayIndex    , nodeCount          , &
         &                                                  parentIndex        , siblingIndex      , siblingMergeeIndex , &
         &                                                  treeCount          , iTree             , iNodeTree          , &
         &                                                  eventCount         , iEvent            , eventNodeIndex
    integer(kind=kind_int8)                              :: nodeIndex          , nodeUniqueID      , nodeUniqueIDMaximum
    type   (varying_string)                              :: message

    ! Read number of trees.
    read (unit,iostat=fileStatus) treeCount
    if (fileStatus < 0) then
       close(unit)
       return
    end if
    ! Create trees
    if (treeCount > 1) then
       treeCurrent => tree
       do iTree=2,treeCount
          allocate(treeCurrent%nextTree)
          treeCurrent => treeCurrent%nextTree
       end do
    end if
    ! Read number of nodes.
    allocate(nodeCountTree(treeCount))
    read (unit) nodeCountTree
    nodeCount=sum(nodeCountTree)
    ! Allocate a list of nodes.
    allocate(nodes(nodeCount))
    ! Iterate over trees.
    treeCurrent => tree
    iNode       =  0
    do iTree=1,treeCount
       ! Read basic tree information.
       read (unit,iostat=fileStatus) treeCurrent%index,treeCurrent%volumeWeight,treeCurrent%initializedUntil,nodeArrayIndex
       ! Create nodes.
       do iNodeTree=1,nodeCountTree(iTree)
          iNode             =  iNode                         +1
          nodes(iNode)%node => treeNode(hostTree=treeCurrent)
       end do
       ! Assign the tree base node.
       if (.not.skipTree) treeCurrent%baseNode => nodes(nodeArrayIndex)%node
       ! Move to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do
    ! Loop over all nodes.
    nodeUniqueIDMaximum=-1
    do iNode=1,nodeCount
       ! Read all node information.
       ! Indices.
       read (unit) nodeIndex,nodeUniqueID
       call nodes(iNode)%node%indexSet   (nodeIndex   )
       call nodes(iNode)%node%uniqueIDSet(nodeUniqueID)
       ! Pointers to other nodes.
       read (unit) parentIndex,firstChildIndex,siblingIndex,firstSatelliteIndex,mergeTargetIndex,firstMergeeIndex&
            &,siblingMergeeIndex,formationNodeIndex
       nodes(iNode)%node%parent         => Pointed_At_Node(parentIndex        ,nodes)
       nodes(iNode)%node%firstChild     => Pointed_At_Node(firstChildIndex    ,nodes)
       nodes(iNode)%node%sibling        => Pointed_At_Node(siblingIndex       ,nodes)
       nodes(iNode)%node%firstSatellite => Pointed_At_Node(firstSatelliteIndex,nodes)
       nodes(iNode)%node%mergeTarget    => Pointed_At_Node(mergeTargetIndex   ,nodes)
       nodes(iNode)%node%firstMergee    => Pointed_At_Node(firstMergeeIndex   ,nodes)
       nodes(iNode)%node%siblingMergee  => Pointed_At_Node(siblingMergeeIndex ,nodes)
       nodes(iNode)%node%formationNode  => Pointed_At_Node(formationNodeIndex ,nodes)
       ! Read the node.
       call nodes(iNode)%node%readRaw(unit)
       ! Find the highest uniqueID.
       nodeUniqueIDMaximum=max(nodeUniqueIDMaximum,nodeUniqueID)
       ! Read any events attached to the node.
       eventPrevious => null()
       read (unit) eventCount
       do iEvent=1,eventCount
          ! Build an event of the correct type from the file.
          event => nodeEventBuildFromRaw(unit)
          ! Read and assign any node pointer that the event may have.
          read (unit) eventNodeIndex
          event%node => Pointed_At_Node(eventNodeIndex,nodes)
          ! Link the event to the node.
          if (iEvent == 1) then
             nodes(iNode)%node%event => event
          else
             eventPrevious%next => event
          end if
          eventPrevious => event
          event         => null()
       end do
    end do    
    ! Set the global maximum unique ID to the maximum found.
    call Galacticus_Nodes_Unique_ID_Set(nodeUniqueIDMaximum)
    ! Perform sanity checks.
    do iNode=1,nodeCount
       if (associated(nodes(iNode)%node%firstChild)) then
          if (.not.associated(nodes(iNode)%node,nodes(iNode)%node%firstChild%parent)) then
             message="child's parent is not self"
             message=message//char(10)//" -> self                : "//nodes(iNode)%node                  %uniqueID()
             message=message//char(10)//" -> self->child         : "//nodes(iNode)%node%firstChild       %uniqueID()
             message=message//char(10)//" -> self->child->parent : "//nodes(iNode)%node%firstChild%parent%uniqueID()
             call Galacticus_Error_Report('Merger_Tree_State_Restore',message)
          end if
       end if
    end do
    ! Destroy the list of nodes.
    deallocate(nodes)
    return
  end subroutine Merger_Tree_State_From_Unit
  
end module Merger_Trees_State_Store
