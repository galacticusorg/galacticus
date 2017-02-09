!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Contains a module which implements an object hierarchy for nodes in merger trees and all of their constituent physical
!% components.

module Galacticus_Nodes
  !% Implements an object hierarchy for nodes in merger trees and all of their constituent physical components.
  use Galacticus_Error
  use Memory_Management
  use ISO_Varying_String
  use Kepler_Orbits
  use Tensors
  use Abundances_Structure
  use Chemical_Abundances_Structure
  use Stellar_Luminosities_Structure
  use Histories
  use Numerical_Constants_Astronomical
  use IO_HDF5
  use Pseudo_Random
  private
  public :: nodeClassHierarchyInitialize, nodeClassHierarchyFinalize, Galacticus_Nodes_Unique_ID_Set, interruptTask

  type, public :: treeNodeList
     !% Type to give a list of treeNodes.
     type(treeNode), pointer :: node
  end type treeNodeList

  ! Include merger tree object.
  include "objects.merger_trees.type.inc"

  ! Include universe class.
  include "objects.universe.type.inc"

  ! Zero dimension arrays to be returned as defaults.
  integer                                            , dimension(0) :: nullInteger1d
  double precision                                   , dimension(0) :: nullDouble1d

  ! Labels for function mapping reduction types.
  integer                         , parameter, public               :: reductionSummation=1
  integer                         , parameter, public               :: reductionProduct  =2

  ! Unique ID counter.
  integer         (kind=kind_int8)                                  :: uniqueIdCount     =0

  ! Event ID counter.
  integer         (kind=kind_int8)                                  :: eventID           =0
  
  ! Define a constructor for treeNodes.
  interface treeNode
     module procedure Tree_Node_Constructor
  end interface treeNode

  ! Include node methods.
  !# <include directive="component" type="component">
  include 'objects.nodes.components.inc'
  !# </include>

  !
  ! Nodes functions.
  subroutine Galacticus_Nodes_Unique_ID_Set(uniqueID)
    !% Resets the global unique ID number.
    implicit none
    integer(kind=kind_int8), intent(in   ) :: uniqueID

    uniqueIdCount=uniqueID
    return
  end subroutine Galacticus_Nodes_Unique_ID_Set

  !
  ! Functions for treeNode class.
  function Tree_Node_Constructor(index,hostTree)
    !% Return a pointer to a newly created and initialized {\normalfont \ttfamily treeNode}.
    implicit none
    type   (treeNode      ), pointer                         :: Tree_Node_Constructor
    integer(kind=kind_int8), intent(in   ), optional         :: index
    type   (mergerTree    ), intent(in   ), optional, target :: hostTree
    integer                                                  :: allocErr

    ! Initialize tree node methods if necessary.
    call nodeClassHierarchyInitialize()

    ! Allocate the object.
    allocate(Tree_Node_Constructor,stat=allocErr)
    if (allocErr/=0) call Galacticus_Error_Report('Tree_Node_Constructor','unable to allocate node')
    call Memory_Usage_Record(sizeof(Tree_Node_Constructor),memoryType=memoryTypeNodes)

    ! Initialize the node.
    call Tree_Node_Constructor%initialize(index,hostTree)
    return
  end function Tree_Node_Constructor

  function Tree_Node_Type(self)
    !% Returns the name of a {\normalfont \ttfamily treeNode} object.
    implicit none
    class(treeNode      ), intent(in   ) :: self
    type (varying_string)                :: Tree_Node_Type
    !GCC$ attributes unused :: self
    
    Tree_Node_Type="treeNode"
    return
  end function Tree_Node_Type

  function Tree_Node_Index(self)
    !% Returns the index of a {\normalfont \ttfamily treeNode}.
    use Galacticus_Error
    implicit none
    class  (treeNode      ), intent(in   ), target :: self
    type   (treeNode      ), pointer               :: workNode
    integer(kind=kind_int8)                        :: Tree_Node_Index

    select type (self)
    type is (treeNode)
       workNode => self
    class default
       workNode => null()
       call Galacticus_Error_Report('Tree_Node_Index','treeNode of unknown class')
    end select
    if (associated(workNode)) then
       Tree_Node_Index=workNode%indexValue
    else
       Tree_Node_Index=-1
    end if
    return
  end function Tree_Node_Index

  subroutine Tree_Node_Index_Set(self,index)
    !% Sets the index of a {\normalfont \ttfamily treeNode}.
    implicit none
    class  (treeNode      ), intent(inout) :: self
    integer(kind=kind_int8), intent(in   ) :: index

    self%indexValue=index
    return
  end subroutine Tree_Node_Index_Set

  function Tree_Node_Unique_ID(self)
    !% Returns the unique ID of a {\normalfont \ttfamily treeNode}.
    use Galacticus_Error
    implicit none
    class  (treeNode      ), intent(in   ), target :: self
    type   (treeNode      ), pointer               :: workNode
    integer(kind=kind_int8)                        :: Tree_Node_Unique_ID

    select type (self)
    type is (treeNode)
       workNode => self
    class default
       workNode => null()
       call Galacticus_Error_Report('Tree_Node_Unique_ID','treeNode of unknown class')
    end select
    if (associated(workNode)) then
       Tree_Node_Unique_ID=workNode%uniqueIdValue
    else
       Tree_Node_Unique_ID=-1
    end if
    return
  end function Tree_Node_Unique_ID

  subroutine Tree_Node_Unique_ID_Set(self,uniqueID)
    !% Sets the index of a {\normalfont \ttfamily treeNode}.
    implicit none
    class  (treeNode      ), intent(inout)           :: self
    integer(kind=kind_int8), intent(in   ), optional :: uniqueID

    if (present(uniqueID)) then
       self%uniqueIdValue=uniqueID
    else
       !$omp critical(UniqueID_Assign)
       uniqueIDCount=uniqueIDCount+1
       if (uniqueIDCount <= 0) call Galacticus_Error_Report('Tree_Node_Unique_ID_Set','ran out of unique ID numbers')
       self%uniqueIdValue=uniqueIDCount
       !$omp end critical(UniqueID_Assign)
    end if
    return
  end subroutine Tree_Node_Unique_ID_Set

  double precision function Tree_Node_Time_Step(self)
    !% Returns the time-step last used by a {\normalfont \ttfamily treeNode}.
    implicit none
    class(treeNode), intent(in   ) :: self

    Tree_Node_Time_Step=self%timeStepValue
    return
  end function Tree_Node_Time_Step

  subroutine Tree_Node_Time_Step_Set(self,timeStep)
    !% Sets the time-step used by a {\normalfont \ttfamily treeNode}.
    implicit none
    class           (treeNode      ), intent(inout) :: self
    double precision                , intent(in   ) :: timeStep

    self%timeStepValue=timeStep
    return
  end subroutine Tree_Node_Time_Step_Set
  
 subroutine Tree_Node_Attach_Event(self,newEvent)
    !% Create a new event in a tree node.
    implicit none
    class(treeNode ), intent(inout)          :: self
    class(nodeEvent), intent(inout), pointer :: newEvent
    class(nodeEvent)               , pointer :: thisEvent

    nullify(newEvent%next)
    !$omp atomic
    eventID=eventID+1
    newEvent%ID=eventID
    if (associated(self%event)) then
       thisEvent => self%event
       do while (associated(thisEvent%next))
          thisEvent => thisEvent%next
       end do
       thisEvent%next => newEvent
    else
       self%event => newEvent
    end if
    return
  end subroutine Tree_Node_Attach_Event

  subroutine Tree_Node_Remove_Paired_Event(self,event)
    !% Removed a paired event from {\normalfont \ttfamily self}. Matching is done on the basis of event ID.
    implicit none
    class  (treeNode ), intent(inout) :: self
    class  (nodeEvent), intent(in   ) :: event
    class  (nodeEvent), pointer       :: lastEvent  , nextEvent, pairEvent
    logical                           :: pairMatched

    ! Locate the paired event in self and remove it.
    pairEvent => self%event
    lastEvent => self%event
    ! Iterate over all events.
    pairMatched=.false.
    do while (associated(pairEvent).and..not.pairMatched)
       ! Match the paired event ID with the current event ID.
       if (pairEvent%ID == event%ID) then
          pairMatched=.true.
          if (associated(pairEvent,self%event)) then
             self     %event => pairEvent%next
             lastEvent       => self     %event
          else
             lastEvent%next  => pairEvent%next
          end if
          nextEvent => pairEvent%next
          deallocate(pairEvent)
          pairEvent => nextEvent
       else
          lastEvent => pairEvent
          pairEvent => pairEvent%next
       end if
    end do
    if (.not.pairMatched) call Galacticus_Error_Report('Tree_Node_Remove_Paired_Event','unable to find paired event')
    return
  end subroutine Tree_Node_Remove_Paired_Event

  logical function Tree_Node_Is_Primary_Progenitor(self)
    !% Returns true if {\normalfont \ttfamily self} is the primary progenitor of its parent node.
    use Galacticus_Error
    implicit none
    class(treeNode), intent(inout) :: self

    select type (self)
    type is (treeNode)
       if (associated(self%parent)) then
          Tree_Node_Is_Primary_Progenitor=associated(self%parent%firstChild,self)
       else
          Tree_Node_Is_Primary_Progenitor=.false.
       end if
    class default
       Tree_Node_Is_Primary_Progenitor=.false.
       call Galacticus_Error_Report('Tree_Node_Is_Primary_Progenitor','treeNode is of unknown class')
    end select
    return
  end function Tree_Node_Is_Primary_Progenitor

  logical function Tree_Node_Is_Primary_Progenitor_Of_Index(self,targetNodeIndex)
    !% Return true if {\normalfont \ttfamily self} is a progenitor of the node with index {\normalfont \ttfamily targetNodeIndex}.
    use Galacticus_Error
    implicit none
    class  (treeNode      ), intent(in   ), target :: self
    integer(kind=kind_int8), intent(in   )         :: targetNodeIndex
    type   (treeNode      ), pointer               :: workNode

    Tree_Node_Is_Primary_Progenitor_Of_Index=.false.
    select type (self)
    type is (treeNode)
       workNode => self
    class default
       workNode => null()
       call Galacticus_Error_Report('Tree_Node_Is_Primary_Progenitor_Of_Index','treeNode of unknown class - this should not happen')
    end select
    do while (associated(workNode))
       if (workNode%index() == targetNodeIndex) then
          Tree_Node_Is_Primary_Progenitor_Of_Index=.true.
          return
       end if
       if (.not.workNode%isPrimaryProgenitor()) return
       workNode => workNode%parent
    end do
    return
  end function Tree_Node_Is_Primary_Progenitor_Of_Index

  logical function Tree_Node_Is_Primary_Progenitor_Of_Node(self,targetNode)
    !% Return true if {\normalfont \ttfamily self} is a progenitor of {\normalfont \ttfamily targetNode}.
    use Galacticus_Error
    implicit none
    class(treeNode), intent(in   ), target  :: self
    type (treeNode), intent(in   ), pointer :: targetNode
    type (treeNode)               , pointer :: workNode

    Tree_Node_Is_Primary_Progenitor_Of_Node=.false.
    select type (self)
    type is (treeNode)
       workNode => self
    class default
       workNode => null()
       call Galacticus_Error_Report('Tree_Node_Is_Primary_Progenitor_Of_Node','treeNode of unknown class - this should not happen')
    end select
    do while (associated(workNode))
       if (associated(workNode,targetNode)) then
          Tree_Node_Is_Primary_Progenitor_Of_Node=.true.
          return
       end if
       if (.not.workNode%isPrimaryProgenitor()) return
       workNode => workNode%parent
    end do
    return
  end function Tree_Node_Is_Primary_Progenitor_Of_Node

  logical function Tree_Node_Is_Progenitor_Of_Index(self,targetNodeIndex)
    !% Return true if {\normalfont \ttfamily self} is a progenitor of the node with index {\normalfont \ttfamily targetNodeIndex}.
    use Galacticus_Error
    implicit none
    class  (treeNode      ), intent(in   ), target :: self
    integer(kind=kind_int8), intent(in   )         :: targetNodeIndex
    type   (treeNode      ), pointer               :: workNode

    Tree_Node_Is_Progenitor_Of_Index=.false.
    select type (self)
    type is (treeNode)
       workNode => self
    class default
       workNode => null()
       call Galacticus_Error_Report('Tree_Node_Is_Progenitor_Of_Index','treeNode of unknown class')
    end select
    do while (associated(workNode))
       if (workNode%index() == targetNodeIndex) then
          Tree_Node_Is_Progenitor_Of_Index=.true.
          return
       end if
       workNode => workNode%parent
    end do
    return
  end function Tree_Node_Is_Progenitor_Of_Index

  logical function Tree_Node_Is_Progenitor_Of_Node(self,targetNode)
    !% Return true if {\normalfont \ttfamily self} is a progenitor of {\normalfont \ttfamily targetNode}.
    use Galacticus_Error
    implicit none
    class(treeNode), intent(in   ), target  :: self
    type (treeNode), intent(in   ), pointer :: targetNode
    type (treeNode)               , pointer :: workNode

    Tree_Node_Is_Progenitor_Of_Node=.false.
    select type (self)
    type is (treeNode)
       workNode => self
    class default
       workNode => null()
       call Galacticus_Error_Report('Tree_Node_Is_Progenitor_Of_Node','treeNode of unknown class - this should not happen')
    end select
    do while (associated(workNode))
       if (associated(workNode,targetNode)) then
          Tree_Node_Is_Progenitor_Of_Node=.true.
          return
       end if
       workNode => workNode%parent
    end do
    return
  end function Tree_Node_Is_Progenitor_Of_Node

  logical function Tree_Node_Is_On_Main_Branch(self)
    !% Returns true if {\normalfont \ttfamily self} is on the main branch.
    use Galacticus_Error
    implicit none
    class(treeNode), intent(inout), target :: self
    type (treeNode), pointer               :: workNode

    Tree_Node_Is_On_Main_Branch=.not.associated(self%parent)
    select type (self)
    type is (treeNode)
       workNode => self
    class default
       workNode => null()
       call Galacticus_Error_Report('Tree_Node_Is_On_Main_Branch','treeNode of unknown class - this should not happen')
    end select
    do while (associated(workNode%parent))
       if (.not.workNode%isPrimaryProgenitor()) return
       workNode => workNode%parent
    end do
    Tree_Node_Is_On_Main_Branch=.true.
    return
  end function Tree_Node_Is_On_Main_Branch

  logical function Tree_Node_Is_Satellite(self)
    !% Returns true if {\normalfont \ttfamily self} is a satellite.
    use Galacticus_Error
    implicit none
    class(treeNode), intent(in   ), target :: self
    type (treeNode), pointer               :: childNode, parentNode, selfActual

    select type (self)
    type is (treeNode)
       selfActual => self
    class default
       selfActual => null()
       call Galacticus_Error_Report('Tree_Node_Is_Satellite','treeNode of unknown class - this should not happen')
    end select
    parentNode => selfActual%parent
    select case (associated(parentNode))
    case (.false.)
       Tree_Node_Is_Satellite=.false.
       return
    case (.true.)
       childNode => parentNode%firstChild
       Tree_Node_Is_Satellite=.true.
       do while (associated(childNode))
          if (associated(childNode,selfActual)) then
             Tree_Node_Is_Satellite=.false.
             exit
          end if
          childNode => childNode%sibling
       end do
    end select
    return
  end function Tree_Node_Is_Satellite

  function Tree_Node_Get_Last_Satellite(self) result (satelliteNode)
    !% Returns a pointer to the final satellite node associated with {\normalfont \ttfamily self}.
    implicit none
    class(treeNode), intent(in   ) :: self
    type (treeNode), pointer       :: satelliteNode

    satelliteNode => self%firstSatellite
    do while (associated(satelliteNode%sibling))
       satelliteNode => satelliteNode%sibling
    end do
    return
  end function Tree_Node_Get_Last_Satellite

  function Tree_Node_Get_Earliest_Progenitor(self) result (progenitorNode)
    !% Returns a pointer to the earliest progenitor of {\normalfont \ttfamily self}.
    use Galacticus_Error
    implicit none
    type (treeNode), pointer       :: progenitorNode
    class(treeNode), intent(inout) :: self

    select type (self)
    type is (treeNode)
       progenitorNode => self
       do while (associated(progenitorNode%firstChild))
          progenitorNode => progenitorNode%firstChild
       end do
    class default
       progenitorNode => null()
       call Galacticus_Error_Report('Tree_Node_Get_Earliest_Progenitor','treeNode of unknown class - this should not happen')
    end select
    return
  end function Tree_Node_Get_Earliest_Progenitor

  function Tree_Node_Merges_With_Node(thisNode)
    !% Returns a pointer to the node with which {\normalfont \ttfamily thisNode} will merge.
    implicit none
    class(treeNode), intent(in   ) :: thisNode
    type (treeNode), pointer       :: Tree_Node_Merges_With_Node

    ! Check if a specific merge node has been set.
    if (associated(thisNode%mergeTarget)) then
       ! One has, so simply return it.
       Tree_Node_Merges_With_Node => thisNode%mergeTarget
    else
       ! No specific merge node has been set, assume merging with the parent node.
       Tree_Node_Merges_With_Node => thisNode%parent
    end if
    return
  end function Tree_Node_Merges_With_Node

  subroutine Tree_Node_Remove_From_Host(self)
    !% Remove {\normalfont \ttfamily self} from the linked list of its host node's satellites.
    use Galacticus_Display
    use Galacticus_Error
    use String_Handling
    implicit none
    class(treeNode      ), intent(in   ), target :: self
    type (treeNode      ), pointer               :: hostNode, previousNode, selfActual, thisNode
    type (varying_string)                        :: message

    select type (self)
    type is (treeNode)
       selfActual => self
    class default
       selfActual => null()
       call Galacticus_Error_Report('Tree_Node_Remove_From_Host','treeNode of unknown class')
    end select

    ! Remove from the parent node satellite list.
    hostNode => selfActual%parent
    message='Satellite node ['
    message=message//selfActual%index()//'] being removed from host node ['//hostNode%index()//']'
    call Galacticus_Display_Message(message,verbosityInfo)
    if (associated(hostNode%firstSatellite,selfActual)) then
       ! This is the first satellite, unlink it, and link to any sibling.
       hostNode%firstSatellite => selfActual%sibling
    else
       thisNode     => hostNode%firstSatellite
       previousNode => null()
       do while (associated(thisNode))
          if (associated(thisNode,selfActual)) then
             ! Found our node, link its older sibling to its younger sibling.
             previousNode%sibling => thisNode%sibling
             exit
          end if
          previousNode => thisNode
          thisNode     => thisNode%sibling
       end do
    end if
    return
  end subroutine Tree_Node_Remove_From_Host

  subroutine Tree_Node_Remove_from_Mergee(self)
    !% Remove {\normalfont \ttfamily self} from the linked list of its host node's satellites.
    use Galacticus_Display
    use Galacticus_Error
    use String_Handling
    implicit none
    class(treeNode      ), intent(in   ), target :: self
    type (treeNode      ), pointer               :: hostNode, previousNode, selfActual, thisNode
    type (varying_string)                        :: message

    select type (self)
    type is (treeNode)
       selfActual => self
    class default
       selfActual => null()
       call Galacticus_Error_Report('Tree_Node_Remove_from_Mergee','treeNode of unknown class')
    end select

    ! Remove from the mergee list of any merge target.
    if (associated(selfActual%mergeTarget)) then
       hostNode => selfActual%mergeTarget
       message='Mergee node ['
       message=message//selfActual%index()//'] being removed from merge target ['//hostNode%index()//']'
       call Galacticus_Display_Message(message,verbosityInfo)
       if (associated(hostNode%firstMergee,selfActual)) then
          ! This is the first mergee, unlink it, and link to any sibling.
          hostNode%firstMergee => selfActual%siblingMergee
       else
          thisNode     => hostNode%firstMergee
          previousNode => null()
          do while (associated(thisNode))
             if (associated(thisNode,selfActual)) then
                ! Found our node, link its older sibling to its younger sibling.
                previousNode%siblingMergee => thisNode%siblingMergee
                exit
             end if
             previousNode => thisNode
             thisNode     => thisNode%siblingMergee
          end do
       end if
    end if
    return
  end subroutine Tree_Node_Remove_from_Mergee

  function treeNodeWalkTree(self)
    !% This function provides a mechanism for walking through an entire merger tree. Given a
    !% pointer {\normalfont \ttfamily self} to a node of the tree, it will return the next node
    !% that should be visited in the tree. Thus, if {\normalfont \ttfamily self} is initially
    !% set to the base of the merger tree and {\normalfont \ttfamily Merger\_Tree\_Walk()} is
    !% called repeatedly it will walk through every node of the tree. Once the entire tree has
    !% been walked, a {\normalfont \ttfamily null()} pointer will be returned, indicating that
    !% there are no more nodes to walk. Each node will be visited once and once only if the tree
    !% is walked in this way.
    implicit none
    type (treeNode)               , pointer :: treeNodeWalkTree
    class(treeNode), intent(inout), target  :: self
    type (treeNode)               , pointer :: workNode

    workNode => self
    if (.not.associated(workNode%parent)) then
       ! This is the base of the merger tree.
       do while (associated(workNode%firstChild))
          workNode => workNode%firstChild
       end do
       if (associated(workNode,self)) nullify(workNode)
    else
       if (associated(workNode%sibling)) then
          workNode => workNode%sibling
          do while (associated(workNode%firstChild))
             workNode => workNode%firstChild
          end do
       else
          workNode => workNode%parent
          ! Terminate when back at tree base.
          if (.not.associated(workNode%parent)) workNode => null()
       end if
    end if
    treeNodeWalkTree => workNode
    return
  end function treeNodeWalkTree
  
  function treeNodeWalkTreeUnderConstruction(self)
    !% This function provides a mechanism for walking through a merger tree that is being built.
    implicit none
    type (treeNode)               , pointer :: treeNodeWalkTreeUnderConstruction
    class(treeNode), intent(inout), target  :: self
    type (treeNode)               , pointer :: workNode

    workNode => self
    if (associated(workNode%firstChild)) then
       ! Move to the primary child if one exists.
       do while (associated(workNode%firstChild))
          workNode => workNode%firstChild
       end do
    else
       if (associated(workNode%sibling)) then
          workNode => workNode%sibling
          do while (associated(workNode%firstChild))
             workNode => workNode%firstChild
          end do
       else
          do while (associated(workNode))
             if (associated(workNode%parent)) then
                workNode => workNode%parent
                if (associated(workNode%sibling)) then
                   workNode => workNode%sibling
                   do while (associated(workNode%firstChild))
                      workNode => workNode%firstChild
                   end do
                   exit
                end if
             else
                workNode => null()
             end if
          end do
       end if
    end if
    treeNodeWalkTreeUnderConstruction => workNode
    return
  end function treeNodeWalkTreeUnderConstruction

  function treeNodeWalkTreeWithSatellites(self)
    !% Merger tree walk function which also descends through satellite nodes. Note that it is
    !% important that the walk descends to satellites before descending to children: the
    !% routines that destroy merger tree branches rely on this since child nodes are used in
    !% testing whether a node is a satellite---if they are destroyed prior to the test being
    !% made then problems with dangling pointers will occur.
    implicit none
    type (treeNode)               , pointer :: treeNodeWalkTreeWithSatellites
    class(treeNode), intent(inout), target  :: self
    type (treeNode)               , pointer :: workNode

    workNode => self
    if (.not.associated(workNode%parent)) then
       ! This is the base of the merger tree.
       ! Descend through satellites and children.
       workNode => Merger_Tree_Walk_Descend_to_Progenitors(workNode)
       if (associated(workNode,self)) nullify(workNode)
    else
       if (associated(workNode%sibling)) then
          workNode => workNode%sibling
          ! Descend through satellites and children.
          workNode => Merger_Tree_Walk_Descend_to_Progenitors(workNode)
       else
          ! About to move back up the tree. Check if the node we're moving up from is a satellite.
          if (workNode%isSatellite()) then
             ! It is a satellite. Therefore, the parent may have children that have yet to be
             ! visited. Check if the parent has children.
             if (associated(workNode%parent%firstChild)) then
                ! Parent does have children, so move to the first one.
                workNode => workNode%parent%firstChild
                ! Descend through satellites and children.
                workNode => Merger_Tree_Walk_Descend_to_Progenitors(workNode)
             else
                ! Parent has no children, so move to the parent.
                workNode => workNode%parent
             end if
          else
             ! It is not a satellite, so all satellites and children have been processed.
             workNode => workNode%parent
          end if
          ! Terminate when back at tree base.
          if (.not.associated(workNode%parent)) workNode => null()
       end if
    end if
    treeNodeWalkTreeWithSatellites => workNode
    return
  end function treeNodeWalkTreeWithSatellites

  function treeNodeWalkBranch(self,startNode)
    !% This function provides a mechanism for walking through the branches of the merger
    !% tree. Given a pointer {\normalfont \ttfamily self} to a branch of the tree, it will
    !% return the next node that should be visited in the tree. Thus, if {\normalfont \ttfamily
    !% self} is initially set to the base of the merger tree and {\normalfont \ttfamily
    !% Merger\_Tree\_Walk\_Branch()} is called repeatedly it will walk through every node of the
    !% branch. Once the entire branch has been walked, a {\normalfont \ttfamily null()} pointer
    !% will be returned, indicating that there are no more nodes to walk. Each node will be
    !% visited once and once only if the branch is walked in this way.
    implicit none
    type (treeNode)               , pointer :: treeNodeWalkBranch
    class(treeNode), intent(inout), target  :: self
    type (treeNode), intent(inout), target  :: startNode
    type (treeNode)               , pointer :: workNode          , selfNode, &
         &                                     branchTip

    selfNode  => self
    workNode  => self
    branchTip => startNode
    if (associated(selfNode,branchTip)) then
       do while (associated(workNode%firstChild))
          workNode => workNode%firstChild
       end do
       if (associated(workNode,selfNode)) nullify(workNode)
    else
       if (associated(workNode%sibling)) then
          workNode => workNode%sibling
          do while (associated(workNode%firstChild))
             workNode => workNode%firstChild
          end do
       else
          workNode => workNode%parent
          ! Terminate when back at starting node.
          if (associated(workNode,branchTip)) workNode => null()
       end if
    end if
    treeNodeWalkBranch => workNode
    return
  end function treeNodeWalkBranch

  function treeNodeWalkBranchWithSatellites(self,startNode)
    !% This function provides a mechanism for walking through the branches of the merger
    !% tree. Given a pointer {\normalfont \ttfamily self} to a branch of the tree, it will
    !% return the next node that should be visited in the tree. Thus, if {\normalfont \ttfamily
    !% self} is initially set to the base of the merger tree and {\normalfont \ttfamily
    !% Merger\_Tree\_Walk\_Branch()} is called repeatedly it will walk through every node of the
    !% branch. Once the entire branch has been walked, a {\normalfont \ttfamily null()} pointer
    !% will be returned, indicating that there are no more nodes to walk. Each node will be
    !% visited once and once only if the branch is walked in this way. Note that it is important
    !% that the walk descends to satellites before descending to children: the routines that
    !% destroy merger tree branches rely on this since child nodes are used in testing whether a
    !% node is a satellite---if they are destroyed prior to the test being made then problems
    !% with dangling pointers will occur.
    implicit none
    type (treeNode)               , pointer :: treeNodeWalkBranchWithSatellites
    class(treeNode), intent(inout), target  :: self
    type (treeNode), intent(inout), pointer :: startNode
    type (treeNode)               , pointer :: workNode                        , selfNode

    selfNode => self
    workNode => self
    if (associated(selfNode,startNode)) then
       ! Descend through satellites and children.
       workNode => Merger_Tree_Walk_Descend_to_Progenitors(workNode)
       if (associated(workNode,selfNode)) nullify(workNode)
    else
       if (associated(workNode%sibling)) then
          workNode => workNode%sibling
          ! Descend through satellites and children.
          workNode => Merger_Tree_Walk_Descend_to_Progenitors(workNode)
       else
          ! About to move back up the tree. Check if the node we're moving up from is a satellite.
          if (workNode%isSatellite()) then
             ! It is a satellite. Therefore, the parent may have children that have yet to be
             ! visited. Check if the parent has children.
             if (associated(workNode%parent%firstChild)) then
                ! Parent does have children, so move to the first one.
                workNode => workNode%parent%firstChild
                ! Descend through satellites and children.
                workNode => Merger_Tree_Walk_Descend_to_Progenitors(workNode)
             else
                ! Parent has no satellites, so move to the parent.
                workNode => workNode%parent
             end if
          else
             ! It is not a satellite, so all satellites and children of the parent must have
             ! been processed. Therefore, move to the parent.
             workNode => workNode%parent
          end if
          ! Terminate when back at starting node.
          if (associated(workNode,startNode)) workNode => null()
       end if
    end if
    treeNodeWalkBranchWithSatellites => workNode
    return
  end function treeNodeWalkBranchWithSatellites
  
  function Merger_Tree_Walk_Descend_to_Progenitors(self) result (progenitorNode)
    !% Descend to the deepest progenitor (satellites and children) of {\normalfont \ttfamily self}.
    implicit none
    type(treeNode), intent(in   ), pointer :: self
    type(treeNode)               , pointer :: progenitorNode

    ! Begin at the input node.
    progenitorNode => self

    ! Descend through satellites and children.
    do while (associated(progenitorNode%firstSatellite).or.associated(progenitorNode%firstChild))
       if (associated(progenitorNode%firstSatellite)) then
          progenitorNode => progenitorNode%firstSatellite
       else
          progenitorNode => progenitorNode%firstChild
       end if
  end do
    return
  end function Merger_Tree_Walk_Descend_to_Progenitors

  subroutine treeNodeDestroyBranch(self)
    !% Destroy the tree branch rooted at this given node.
    implicit none
    class(treeNode), intent(inout), target  :: self
    type (treeNode)               , pointer :: nodeDestroy, nodeNext, &
         &                                     branchTip
    
    ! Descend to the tip of the branch.
    branchTip => self
    nodeNext  => branchTip%walkBranchWithSatellites(branchTip)
    ! Loop over all tree nodes.
    do while (associated(nodeNext))
       ! Keep of a record of the current node, so that we can destroy it.
       nodeDestroy => nodeNext
       ! Walk to the next node in the tree.
       nodeNext => nodeDestroy%walkBranchWithSatellites(branchTip)
       ! If the node about to be destroyed is the primary progenitor of its parent we must move the child pointer of the parent to
       ! point to the node's sibling. This is necessary as parent-child pointers are used to establish satellite status and so
       ! will be utilized when walking the tree. Failure to do this can result in attempts to use dangling pointers.
       if (associated(nodeDestroy%parent).and.associated(nodeDestroy%parent%firstChild,nodeDestroy)) &
            & nodeDestroy%parent%firstChild => nodeDestroy%sibling
       ! Destroy the current node.
       call nodeDestroy%destroy()
       deallocate(nodeDestroy)
    end do
    ! Destroy the base node of the branch.
    if (associated(self%parent).and.associated(self%parent%firstChild,self)) self%parent%firstChild => self%sibling
    call self%destroy()
    return
  end subroutine treeNodeDestroyBranch

  !
  ! Functions for nodeComponent class.
  function Node_Component_Generic_Type(self)
    !% Returns the name of a generic tree node component.
    implicit none
    class(nodeComponent ), intent(in   ) :: self
    type (varying_string)                :: Node_Component_Generic_Type
    !GCC$ attributes unused :: self

    Node_Component_Generic_Type="nodeComponent"
    return
  end function Node_Component_Generic_Type

  subroutine Node_Component_Generic_Destroy(self)
    !% Destroy a generic tree node component.
    implicit none
    class(nodeComponent), intent(inout) :: self
    !GCC$ attributes unused :: self
    
    ! Do nothing.
    return
  end subroutine Node_Component_Generic_Destroy

  subroutine Node_Component_ODE_Step_Initialize_Null(self)
    !% Initialize a generic tree node component for an ODE solver step.
    implicit none
    class(nodeComponent), intent(inout) :: self
    !GCC$ attributes unused :: self
    
    return
  end subroutine Node_Component_ODE_Step_Initialize_Null

  subroutine Node_Component_Dump_Null(self)
    !% Dump a generic tree node component.
    implicit none
    class(nodeComponent), intent(in   ) :: self
    !GCC$ attributes unused :: self

    return
  end subroutine Node_Component_Dump_Null

  subroutine Node_Component_Dump_XML_Null(self,fileHandle)
    !% Dump a generic tree node component to XML.
    implicit none
    class  (nodeComponent), intent(inout) :: self
    integer               , intent(in   ) :: fileHandle
    !GCC$ attributes unused :: self, fileHandle

    return
  end subroutine Node_Component_Dump_XML_Null

  subroutine Node_Component_Dump_Raw_Null(self,fileHandle)
    !% Dump a generic tree node component in binary.
    implicit none
    class  (nodeComponent), intent(in   ) :: self
    integer               , intent(in   ) :: fileHandle
    !GCC$ attributes unused :: self, fileHandle

    return
  end subroutine Node_Component_Dump_Raw_Null

  subroutine Node_Component_Read_Raw_Null(self,fileHandle)
    !% Read a generic tree node component in binary.
    implicit none
    class  (nodeComponent), intent(inout) :: self
    integer               , intent(in   ) :: fileHandle
    !GCC$ attributes unused :: self, fileHandle

    return
  end subroutine Node_Component_Read_Raw_Null

  subroutine Node_Component_Output_Count_Null(self,integerPropertyCount,doublePropertyCount,time,instance)
    !% Dump a generic tree node component.
    implicit none
    class           (nodeComponent), intent(inout) :: self
    integer                        , intent(inout) :: doublePropertyCount, integerPropertyCount
    double precision               , intent(in   ) :: time
    integer                        , intent(in   ) :: instance
    !GCC$ attributes unused :: self, integerPropertyCount, doublePropertyCount, time, instance

    return
  end subroutine Node_Component_Output_Count_Null

  subroutine Node_Component_Output_Names_Null(self,integerProperty,integerPropertyNames,integerPropertyComments &
       &,integerPropertyUnitsSI,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time,instance)
    !% Dump a generic tree node component.
    implicit none
    class           (nodeComponent )              , intent(inout) :: self
    double precision                              , intent(in   ) :: time
    integer                                       , intent(inout) :: doubleProperty         , integerProperty
    character       (len=*         ), dimension(:), intent(inout) :: doublePropertyComments , doublePropertyNames   , &
         &                                                           integerPropertyComments, integerPropertyNames
    double precision                , dimension(:), intent(inout) :: doublePropertyUnitsSI  , integerPropertyUnitsSI
    integer                                       , intent(in   ) :: instance
    !GCC$ attributes unused :: self, integerProperty, integerPropertyNames, integerPropertyComments, integerPropertyUnitsSI, doubleProperty, doublePropertyNames, doublePropertyComments, doublePropertyUnitsSI, time, instance
    
    return
  end subroutine Node_Component_Output_Names_Null

  subroutine Node_Component_Output_Null(self,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time,outputInstance,instance)
    !% Dump a generic tree node component.
    use Multi_Counters
    implicit none
    class           (nodeComponent    ), intent(inout) :: self
    double precision                   , intent(in   ) :: time
    integer                            , intent(inout) :: doubleBufferCount     , doubleProperty, integerBufferCount, &
         &                                                integerProperty
    integer         (kind=kind_int8   ), intent(inout) :: integerBuffer    (:,:)
    double precision                   , intent(inout) :: doubleBuffer     (:,:)
    type            (multiCounter     ), intent(in   ) :: outputInstance
    integer                            , intent(in   ) :: instance
    !GCC$ attributes unused :: self, integerProperty, integerBufferCount, integerBuffer, doubleProperty, doubleBufferCount, doubleBuffer, time, outputInstance, instance

    return
  end subroutine Node_Component_Output_Null

  integer function Node_Component_Serialize_Count_Zero(self)
    !% Return the serialization count of a generic tree node component.
    implicit none
    class(nodeComponent), intent(in   ) :: self
    !GCC$ attributes unused :: self
    
    Node_Component_Serialize_Count_Zero=0
    return
  end function Node_Component_Serialize_Count_Zero

  subroutine Node_Component_Serialization_Offsets(self,count)
    !% Return the serialization count of a generic tree node component.
    implicit none
    class  (nodeComponent), intent(in   ) :: self
    integer               , intent(inout) :: count
    !GCC$ attributes unused :: self, count
    
    return
  end subroutine Node_Component_Serialization_Offsets

  subroutine Node_Component_Serialize_Null(self,array)
    !% Serialize a generic tree node component.
    implicit none
    class           (nodeComponent)              , intent(in   ) :: self
    double precision               , dimension(:), intent(  out) :: array
    !GCC$ attributes unused :: self, array
    
    return
  end subroutine Node_Component_Serialize_Null

  subroutine Node_Component_Deserialize_Null(self,array)
    !% Deserialize a generic tree node component.
    implicit none
    class           (nodeComponent)              , intent(inout) :: self
    double precision               , dimension(:), intent(in   ) :: array
    !GCC$ attributes unused :: self, array
    
    return
  end subroutine Node_Component_Deserialize_Null

  function Node_Component_Host_Node(self)
    !% Return the host tree node of a tree node component.
    implicit none
    class(nodeComponent), intent(in   ) :: self
    type (treeNode     ), pointer       :: Node_Component_Host_Node

    Node_Component_Host_Node => self%hostNode
    return
  end function Node_Component_Host_Node

  subroutine Node_Component_Null_Void0_InOut(self)
    !% A null {\normalfont \ttfamily void} function for rank 0 {\normalfont \ttfamily nodeComponent} arrays.
    implicit none
    class(nodeComponent), intent(inout) :: self
    !GCC$ attributes unused :: self
    
    return
  end subroutine Node_Component_Null_Void0_InOut

  double precision function Node_Component_Null_Double0_InOut(self)
    !% A null {\normalfont \ttfamily double} function for rank 0 {\normalfont \ttfamily nodeComponent} arrays.
    implicit none
    class(nodeComponent), intent(inout) :: self
    !GCC$ attributes unused :: self

    Node_Component_Null_Double0_InOut=0.0d0
    return
  end function Node_Component_Null_Double0_InOut

  double precision function Node_Component_Enclosed_Mass_Null(self,radius,componentType,massType,weightBy,weightIndex,haloLoaded)
    !% A null implementation of the enclosed mass in a component. Always returns zero.
    implicit none
    class           (nodeComponent), intent(inout)           :: self
    integer                        , intent(in   )           :: componentType, massType, weightBy, weightIndex
    double precision               , intent(in   )           :: radius
    logical                        , intent(in   ), optional :: haloLoaded
    !GCC$ attributes unused :: self, radius, componentType, massType, weightBy, weightIndex, haloLoaded
    
    Node_Component_Enclosed_Mass_Null=0.0d0
    return
  end function Node_Component_Enclosed_Mass_Null

  double precision function Node_Component_Density_Null(self,positionSpherical,componentType,massType,weightBy,weightIndex,haloLoaded)
    !% A null implementation of the density in a component. Always returns zero.
    implicit none
    class           (nodeComponent)              , intent(inout)           :: self
    integer                                      , intent(in   )           :: componentType    , massType, weightBy, &
         &                                                                    weightIndex
    double precision               , dimension(3), intent(in   )           :: positionSpherical
    logical                                      , intent(in   ), optional :: haloLoaded
    !GCC$ attributes unused :: self, positionSpherical, componentType, massType, weightBy, weightIndex, haloLoaded
    
    Node_Component_Density_Null=0.0d0
    return
  end function Node_Component_Density_Null

  double precision function Node_Component_Surface_Density_Null(self,positionCylindrical,componentType,massType,weightBy,weightIndex,haloLoaded)
    !% A null implementation of the surface density in a component. Always returns zero.
    implicit none
    class           (nodeComponent)              , intent(inout)           :: self
    integer                                      , intent(in   )           :: componentType      , massType   , &
         &                                                                    weightBy           , weightIndex
    double precision               , dimension(3), intent(in   )           :: positionCylindrical
    logical                                      , intent(in   ), optional :: haloLoaded
    !GCC$ attributes unused :: self, positionCylindrical, componentType, massType, weightBy, weightIndex, haloLoaded
    
    Node_Component_Surface_Density_Null=0.0d0
    return
  end function Node_Component_Surface_Density_Null

  double precision function Node_Component_Potential_Null(self,radius,componentType,massType,haloLoaded,status)
    !% A null implementation of the gravitational potential in a component. Always returns zero.
    implicit none
    class           (nodeComponent), intent(inout)           :: self
    integer                        , intent(in   )           :: componentType, massType
    double precision               , intent(in   )           :: radius
    logical                        , intent(in   ), optional :: haloLoaded
    integer                        , intent(inout), optional :: status
    !GCC$ attributes unused :: self, radius, componentType, massType, haloLoaded, status
    
    Node_Component_Potential_Null=0.0d0
    return
  end function Node_Component_Potential_Null

  double precision function Node_Component_Rotation_Curve_Null(self,radius,componentType,massType,haloLoaded)
    !% A null implementation of the rotation curve due to a component. Always returns zero.
    implicit none
    class           (nodeComponent), intent(inout)           :: self
    integer                        , intent(in   )           :: componentType, massType
    double precision               , intent(in   )           :: radius
    logical                        , intent(in   ), optional :: haloLoaded
    !GCC$ attributes unused :: self, radius, componentType, massType, haloLoaded
    
    Node_Component_Rotation_Curve_Null=0.0d0
    return
  end function Node_Component_Rotation_Curve_Null

  double precision function Node_Component_Rotation_Curve_Gradient_Null(self,radius,componentType,massType,haloLoaded)
    !% A null implementation of the gradient of the rotation curve due to a component. Always returns zero.
    implicit none
    class           (nodeComponent), intent(inout)           :: self
    integer                        , intent(in   )           :: componentType, massType
    double precision               , intent(in   )           :: radius
    logical                        , intent(in   ), optional :: haloLoaded
    !GCC$ attributes unused :: self, radius, componentType, massType, haloLoaded
    
    Node_Component_Rotation_Curve_Gradient_Null=0.0d0
    return
  end function Node_Component_Rotation_Curve_Gradient_Null

  ! Simple Boolean functions.
  logical function Boolean_False()
    !% Returns Boolean false always.
    implicit none

    Boolean_False=.false.
    return
  end function Boolean_False

  logical function Boolean_True()
    !% Returns Boolean true always.
    implicit none

    Boolean_True=.true.
    return
  end function Boolean_True

  ! Include functions for the merger tree class.
  include "objects.merger_trees.functions.inc"

  ! Include functions for the universe class.
  include "objects.universe.functions.inc"

end module Galacticus_Nodes
