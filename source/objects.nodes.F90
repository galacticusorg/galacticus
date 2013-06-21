!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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
  use Kind_Numbers
  use Kepler_Orbits
  use Abundances_Structure
  use Chemical_Abundances_Structure
  use Histories
  use Numerical_Constants_Units
  use Numerical_Constants_Prefixes
  use Numerical_Constants_Astronomical
  use Numerical_Constants_Physical
  private
  public :: Galacticus_Nodes_Initialize, Galacticus_Nodes_Finalize, Interrupt_Procedure_Template
  !! <gfortran4.8> workaround
  public :: assignment(=)

  type, private :: nodeDataLogicalScalar
     !% Type describing a non-evolvable scalar logical property of a node component.
     logical :: value
  end type nodeDataLogicalScalar

  type, private :: nodeDataIntegerScalar
     !% Type describing a non-evolvable scalar integer property of a node component.
     integer :: value
  end type nodeDataIntegerScalar

  type, private :: nodeDataDoubleScalarEvolvable
     !% Type describing an evolvable scalar double property of a node component.
     double precision :: rate, scale, value
  end type nodeDataDoubleScalarEvolvable

  type, private :: nodeDataDoubleScalar
     !% Type describing a non-evolvable scalar double property of a node component.
     double precision :: value
  end type nodeDataDoubleScalar

  type, private :: nodeDataDouble1dEvolvable
     !% Type describing an evolvable 1-D double property of a node component.
     double precision, allocatable, dimension(:) :: rate, scale, value
  end type nodeDataDouble1dEvolvable

  type, private :: nodeDataDouble1d
     !% Type describing a non-evolvable 1-D double property of a node component.
     double precision, allocatable, dimension(:) :: value
  end type nodeDataDouble1d

  type, private :: nodeDataAbundancesScalarEvolvable
     !% Type describing an evolvable scalar {\\tt abundances} property of a node component.
     type(abundances        ) :: rate, scale, value
  end type nodeDataAbundancesScalarEvolvable

  type, private :: nodeDataChemicalAbundancesScalarEvolvable
     !% Type describing an evolvable scalar {\\tt chemicalAbundances} property of a node component.
     type(chemicalAbundances) :: rate, scale, value
  end type nodeDataChemicalAbundancesScalarEvolvable

  type, private :: nodeDataHistoryScalarEvolvable
     !% Type describing an evolvable scalar {\\tt history} property of a node component.
     type(history           ) :: rate, scale, value
  end type nodeDataHistoryScalarEvolvable

  type, private :: nodeDataHistoryScalar
     !% Type describing an non-evolvable scalar {\\tt history} property of a node component.
     type(history           ) :: value
  end type nodeDataHistoryScalar

  type, private :: nodeDataKeplerOrbitScalar
     !% Type describing a non-evolvable scalar {\\tt keplerOrbit} property of a node component.
     type(keplerOrbit       ) :: value
  end type nodeDataKeplerOrbitScalar

  type, public :: treeNodeList
     !% Type to give a list of treeNodes.
     type(treeNode), pointer :: node
  end type treeNodeList

  ! Zero dimension arrays to be returned as defaults.
  double precision                , dimension(0)         :: nullDouble1d

  ! Labels for function mapping reduction types.
  integer                         , parameter   , public :: reductionSummation=1
  integer                         , parameter   , public :: reductionProduct  =2

  ! Unique ID counter.
  integer         (kind=kind_int8)                       :: uniqueIdCount     =0

  ! Event ID counter.
  integer         (kind=kind_int8)                       :: eventID           =0

  ! Include node methods.
  !# <include directive="component" type="component">
  include 'objects.nodes.components.inc'
  !# </include>

  !
  ! Functions for treeNode class.
  function Tree_Node_Type(self)
    !% Returns the name of a {\tt treeNode} object.
    implicit none
    class(treeNode      ), intent(in   ) :: self
    type (varying_string)                :: Tree_Node_Type

    Tree_Node_Type="treeNode"
    return
  end function Tree_Node_Type

  function Tree_Node_Index(self)
    !% Returns the index of a {\tt treeNode}.
    implicit none
    class  (treeNode      ), intent(in   ), target :: self
    type   (treeNode      ), pointer               :: workNode
    integer(kind=kind_int8)                        :: Tree_Node_Index

    select type (self)
    type is (treeNode)
       workNode => self
    end select
    if (associated(workNode)) then
       Tree_Node_Index=workNode%indexValue
    else
       Tree_Node_Index=-1
    end if
    return
  end function Tree_Node_Index

  subroutine Tree_Node_Index_Set(self,index)
    !% Sets the index of a {\tt treeNode}.
    implicit none
    class  (treeNode      ), intent(inout) :: self
    integer(kind=kind_int8), intent(in   ) :: index

    self%indexValue=index
    return
  end subroutine Tree_Node_Index_Set

  function Tree_Node_Unique_ID(self)
    !% Returns the unique ID of a {\tt treeNode}.
    implicit none
    class  (treeNode      ), intent(in   ) :: self
    integer(kind=kind_int8)                :: Tree_Node_Unique_ID

    Tree_Node_Unique_ID=self%uniqueIdValue
    return
  end function Tree_Node_Unique_ID

  subroutine Tree_Node_Unique_ID_Set(self,uniqueID)
    !% Sets the index of a {\tt treeNode}.
    use Galacticus_Error
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

  function Tree_Node_Create_Event(self) result (newEvent)
    !% Create a new event in a tree node.
    implicit none
    type (nodeEvent), pointer       :: newEvent, thisEvent
    class(treeNode ), intent(inout) :: self

    allocate(newEvent)
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
  end function Tree_Node_Create_Event

  subroutine Tree_Node_Remove_Paired_Event(self,event)
    !% Removed a paired event from {\tt self}. Matching is done on the basis of event ID.
    implicit none
    class  (treeNode ), intent(inout) :: self
    type   (nodeEvent), intent(in   ) :: event
    type   (nodeEvent), pointer       :: lastEvent  , nextEvent, pairEvent
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
    !% Returns true if {\tt self} is the primary progenitor of its parent node.
    implicit none
    class(treeNode), intent(inout) :: self

    select type (self)
    type is (treeNode)
       if (associated(self%parent)) then
          Tree_Node_Is_Primary_Progenitor=associated(self%parent%firstChild,self)
       else
          Tree_Node_Is_Primary_Progenitor=.false.
       end if
    end select
    return
  end function Tree_Node_Is_Primary_Progenitor

  logical function Tree_Node_Is_Primary_Progenitor_Of_Index(self,targetNodeIndex)
    !% Return true if {\tt self} is a progenitor of the node with index {\tt targetNodeIndex}.
    implicit none
    class  (treeNode      ), intent(in   ), target :: self
    integer(kind=kind_int8), intent(in   )         :: targetNodeIndex
    type   (treeNode      ), pointer               :: workNode

    Tree_Node_Is_Primary_Progenitor_Of_Index=.false.
    select type (self)
    type is (treeNode)
       workNode => self
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
    !% Return true if {\tt self} is a progenitor of {\tt targetNode}.
    implicit none
    class(treeNode), intent(in   ), target  :: self
    type (treeNode), intent(in   ), pointer :: targetNode
    type (treeNode)               , pointer :: workNode

    Tree_Node_Is_Primary_Progenitor_Of_Node=.false.
    select type (self)
    type is (treeNode)
       workNode => self
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

  logical function Tree_Node_Is_On_Main_Branch(self)
    !% Returns true if {\tt self} is on the main branch.
    implicit none
    class(treeNode), intent(inout), target :: self
    type (treeNode), pointer               :: workNode

    Tree_Node_Is_On_Main_Branch=.not.associated(self%parent)
    select type (self)
    type is (treeNode)
       workNode => self
    end select
    do while (associated(workNode%parent))
       if (.not.workNode%isPrimaryProgenitor()) return
       workNode => workNode%parent
    end do
    Tree_Node_Is_On_Main_Branch=.true.
    return
  end function Tree_Node_Is_On_Main_Branch

  logical function Tree_Node_Is_Satellite(self)
    !% Returns true if {\tt self} is a satellite.
    implicit none
    class(treeNode), intent(in   ), target :: self
    type (treeNode), pointer               :: childNode, parentNode, selfActual

    select type (self)
    type is (treeNode)
       selfActual => self
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
    !% Returns a pointer to the final satellite node associated with {\tt self}.
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
    !% Returns a pointer to the earliest progenitor of {\tt self}.
    implicit none
    type (treeNode), pointer       :: progenitorNode
    class(treeNode), intent(inout) :: self

    select type (self)
    type is (treeNode)
       progenitorNode => self
       do while (associated(progenitorNode%firstChild))
          progenitorNode => progenitorNode%firstChild
       end do
    end select
    return
  end function Tree_Node_Get_Earliest_Progenitor

  function Tree_Node_Merges_With_Node(thisNode)
    !% Returns a pointer to the node with which {\tt thisNode} will merge.
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
    !% Remove {\tt self} from the linked list of its host node's satellites.
    use Galacticus_Display
    use ISO_Varying_String
    use String_Handling
    implicit none
    class(treeNode      ), intent(in   ), target :: self
    type (treeNode      ), pointer               :: hostNode, previousNode, selfActual, thisNode
    type (varying_string)                        :: message

    select type (self)
    type is (treeNode)
       selfActual => self
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
    !% Remove {\tt self} from the linked list of its host node's satellites.
    use Galacticus_Display
    use ISO_Varying_String
    use String_Handling
    implicit none
    class(treeNode      ), intent(in   ), target :: self
    type (treeNode      ), pointer               :: hostNode, previousNode, selfActual, thisNode
    type (varying_string)                        :: message

    select type (self)
    type is (treeNode)
       selfActual => self
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

  subroutine Tree_Node_Walk_Tree(self,nextNode)
    !% This function provides a mechanism for walking through an entire merger tree. Given a pointer {\tt self}
    !% to a node of the tree, it will return the next node that should be visited in the tree. Thus, if {\tt self} is
    !% initially set to the base of the merger tree and {\tt Merger\_Tree\_Walk()} is called repeatedly it will walk through every node
    !% of the tree. Once the entire tree has been walked, a {\tt null()} pointer will be returned, indicating that there
    !% are no more nodes to walk. Each node will be visited once and once only if the tree is walked in this way.
    implicit none
    class(treeNode), intent(inout), target  :: self
    type (treeNode), intent(inout), pointer :: nextNode
    type (treeNode)               , pointer :: selfActual, workNode

    select type (self)
    type is (treeNode)
       selfActual => self
    end select
    workNode => selfActual

    if (.not.associated(workNode%parent)) then
       ! This is the base of the merger tree.
       do while (associated(workNode%firstChild))
          workNode => workNode%firstChild
       end do
       if (associated(workNode,selfActual)) nullify(workNode)
    else
       if (associated(workNode%sibling)) then
          workNode => workNode%sibling
          do while (associated(workNode%firstChild))
             workNode => workNode%firstChild
          end do
       else
          workNode => workNode%parent
          if (.not.associated(workNode%parent)) workNode => null() ! Terminate when back at tree base.
       end if
    end if
    nextNode => workNode
    return
  end subroutine Tree_Node_Walk_Tree

  subroutine Tree_Node_Walk_Tree_Under_Construction(self,nextNode)
    !% This function provides a mechanism for walking through a merger tree that is being built.
    implicit none
    class(treeNode), intent(inout), target  :: self
    type (treeNode), intent(inout), pointer :: nextNode
    type (treeNode)               , pointer :: workNode

    select type (self)
    type is (treeNode)
       workNode => self
    end select

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
    nextNode => workNode
    return
  end subroutine Tree_Node_Walk_Tree_Under_Construction

  subroutine Tree_Node_Walk_Tree_With_Satellites(self,nextNode)
    !% Merger tree walk function which also descends through satellite nodes. Note that it is important that the walk descends to
    !% satellites before descending to children: the routines that destroy merger tree branches rely on this since child nodes are
    !% used in testing whether a node is a satellite---if they are destroyed prior to the test being made then problems with
    !% dangling pointers will occur.
    implicit none
    class(treeNode), intent(inout), target  :: self
    type (treeNode), intent(inout), pointer :: nextNode
    type (treeNode)               , pointer :: selfActual, workNode

    select type (self)
    type is (treeNode)
       selfActual => self
    end select
    workNode => selfActual
    if (.not.associated(workNode%parent)) then
       ! This is the base of the merger tree.
       ! Descend through satellites and children.
       workNode => Merger_Tree_Walk_Descend_to_Progenitors(workNode)
       if (associated(workNode,selfActual)) nullify(workNode)
    else
       if (associated(workNode%sibling)) then
          workNode => workNode%sibling
          ! Descend through satellites and children.
          workNode => Merger_Tree_Walk_Descend_to_Progenitors(workNode)
       else
          ! About to move back up the tree. Check if the node we're moving up from is a satellite.
          if (workNode%isSatellite()) then
             ! It is a satellite. Therefore, the parent may have children that have yet to be visited. Check if the parent has
             ! children.
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
          if (.not.associated(workNode%parent)) workNode => null() ! Terminate when back at tree base.
       end if
    end if
    nextNode => workNode
    return
  end subroutine Tree_Node_Walk_Tree_With_Satellites

  subroutine Tree_Node_Walk_Branch(self,startNode,nextNode)
    !% This function provides a mechanism for walking through the branches of the merger tree. Given a pointer {\tt self}
    !% to a branch of the tree, it will return the next node that should be visited in the tree. Thus, if {\tt self} is
    !% initially set to the base of the merger tree and {\tt Merger\_Tree\_Walk\_Branch()} is called repeatedly it will walk through every node
    !% of the branch. Once the entire branch has been walked, a {\tt null()} pointer will be returned, indicating that there
    !% are no more nodes to walk. Each node will be visited once and once only if the branch is walked in this way.
    implicit none
    class(treeNode), intent(inout), target  :: self
    type (treeNode), intent(inout), pointer :: nextNode  , startNode
    type (treeNode)               , pointer :: selfActual, workNode

    select type (self)
    type is (treeNode)
       selfActual => self
    end select
    workNode => selfActual

    if (associated(selfActual,startNode)) then
       do while (associated(workNode%firstChild))
          workNode => workNode%firstChild
       end do
       if (associated(workNode,selfActual)) nullify(workNode)
    else
       if (associated(workNode%sibling)) then
          workNode => workNode%sibling
          do while (associated(workNode%firstChild))
             workNode => workNode%firstChild
          end do
       else
          workNode => workNode%parent
          if (associated(workNode,startNode)) workNode => null() ! Terminate when back at starting node.
       end if
    end if
    nextNode => workNode
    return
  end subroutine Tree_Node_Walk_Branch

  subroutine Tree_Node_Walk_Branch_With_Satellites(self,startNode,nextNode)
    !% This function provides a mechanism for walking through the branches of the merger tree. Given a pointer {\tt self} to a
    !% branch of the tree, it will return the next node that should be visited in the tree. Thus, if {\tt self} is initially
    !% set to the base of the merger tree and {\tt Merger\_Tree\_Walk\_Branch()} is called repeatedly it will walk through every
    !% node of the branch. Once the entire branch has been walked, a {\tt null()} pointer will be returned, indicating that there
    !% are no more nodes to walk. Each node will be visited once and once only if the branch is walked in this way. Note that it
    !% is important that the walk descends to satellites before descending to children: the routines that destroy merger tree
    !% branches rely on this since child nodes are used in testing whether a node is a satellite---if they are destroyed prior to
    !% the test being made then problems with dangling pointers will occur.
    implicit none
    class(treeNode), intent(inout), target  :: self
    type (treeNode), intent(inout), pointer :: nextNode  , startNode
    type (treeNode)               , pointer :: selfActual, workNode

    select type (self)
    type is (treeNode)
       selfActual => self
    end select
    workNode => selfActual

    if (associated(selfActual,startNode)) then
       ! Descend through satellites and children.
       workNode => Merger_Tree_Walk_Descend_to_Progenitors(workNode)
       if (associated(workNode,selfActual)) nullify(workNode)
    else
       if (associated(workNode%sibling)) then
          workNode => workNode%sibling
          ! Descend through satellites and children.
          workNode => Merger_Tree_Walk_Descend_to_Progenitors(workNode)
       else
          ! About to move back up the tree. Check if the node we're moving up from is a satellite.
          if (workNode%isSatellite()) then
             ! It is a satellite. Therefore, the parent may have children that have yet to be visited. Check if the parent
             ! has children.
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
             ! It is not a satellite, so all satellites and children of the parent must have been processed. Therefore, move to
             ! the parent.
             workNode => workNode%parent
          end if
          if (associated(workNode,startNode)) workNode => null() ! Terminate when back at starting node.
       end if
    end if
    nextNode => workNode
    return
  end subroutine Tree_Node_Walk_Branch_With_Satellites

  function Merger_Tree_Walk_Descend_to_Progenitors(self) result (progenitorNode)
    !% Descend to the deepest progenitor (satellites and children) of {\tt self}.
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

  !
  ! Functions for nodeComponent class.
  function Node_Component_Generic_Type(self)
    !% Returns the name of a generic tree node component.
    implicit none
    class(nodeComponent ), intent(in   ) :: self
    type (varying_string)                :: Node_Component_Generic_Type

    Node_Component_Generic_Type="nodeComponent"
    return
  end function Node_Component_Generic_Type

  subroutine Node_Component_Generic_Destroy(self)
    !% Destroy a generic tree node component.
    implicit none
    class(nodeComponent), intent(inout) :: self

    ! Do nothing.
    return
  end subroutine Node_Component_Generic_Destroy

  subroutine Node_Component_ODE_Step_Initialize_Null(self)
    !% Initialize a generic tree node component for an ODE solver step.
    implicit none
    class(nodeComponent), intent(inout) :: self

    return
  end subroutine Node_Component_ODE_Step_Initialize_Null

  subroutine Node_Component_Dump_Null(self)
    !% Dump a generic tree node component.
    implicit none
    class(nodeComponent), intent(in   ) :: self

    return
  end subroutine Node_Component_Dump_Null

  subroutine Node_Component_Dump_Raw_Null(self,fileHandle)
    !% Dump a generic tree node component in binary.
    implicit none
    class  (nodeComponent), intent(in   ) :: self
    integer               , intent(in   ) :: fileHandle

    return
  end subroutine Node_Component_Dump_Raw_Null

  subroutine Node_Component_Output_Count_Null(self,integerPropertyCount,doublePropertyCount,time,instance)
    !% Dump a generic tree node component.
    implicit none
    class           (nodeComponent), intent(inout) :: self
    integer                        , intent(inout) :: doublePropertyCount, integerPropertyCount
    double precision               , intent(in   ) :: time
    integer                        , intent(in   ) :: instance

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

    return
  end subroutine Node_Component_Output_Names_Null

  subroutine Node_Component_Output_Null(self,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time,instance)
    !% Dump a generic tree node component.
    implicit none
    class           (nodeComponent    ), intent(inout) :: self
    double precision                   , intent(in   ) :: time
    integer                            , intent(inout) :: doubleBufferCount     , doubleProperty, integerBufferCount, &
         &                                                integerProperty
    integer         (kind=kind_int8   ), intent(inout) :: integerBuffer    (:,:)
    double precision                   , intent(inout) :: doubleBuffer     (:,:)
    integer                            , intent(in   ) :: instance

    return
  end subroutine Node_Component_Output_Null

  integer function Node_Component_Serialize_Count_Zero(self)
    !% Return the serialization count of a generic tree node component.
    implicit none
    class(nodeComponent), intent(in   ) :: self

    Node_Component_Serialize_Count_Zero=0
    return
  end function Node_Component_Serialize_Count_Zero

  subroutine Node_Component_Serialize_Null(self,array)
    !% Serialize a generic tree node component.
    implicit none
    class           (nodeComponent)              , intent(in   ) :: self
    double precision               , dimension(:), intent(  out) :: array

    return
  end subroutine Node_Component_Serialize_Null

  subroutine Node_Component_Deserialize_Null(self,array)
    !% Deserialize a generic tree node component.
    implicit none
    class           (nodeComponent)              , intent(inout) :: self
    double precision               , dimension(:), intent(in   ) :: array

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
    !% A null {\tt void} function for rank 0 {\tt nodeComponent} arrays.
    implicit none
    class(nodeComponent), intent(inout) :: self

    return
  end subroutine Node_Component_Null_Void0_InOut

  double precision function Node_Component_Null_Double0_InOut(self)
    !% A null {\tt double} function for rank 0 {\tt nodeComponent} arrays..
    implicit none
    class(nodeComponent), intent(inout) :: self

    return
  end function Node_Component_Null_Double0_InOut

  double precision function Node_Component_Enclosed_Mass_Null(self,radius,componentType,massType,weightBy,weightIndex,haloLoaded)
    !% A null implementation of the enclosed mass in a component. Always returns zero.
    implicit none
    class           (nodeComponent), intent(inout)           :: self
    integer                        , intent(in   )           :: componentType, massType, weightBy, weightIndex
    double precision               , intent(in   )           :: radius
    logical                        , intent(in   ), optional :: haloLoaded

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

    Node_Component_Density_Null=0.0d0
    return
  end function Node_Component_Density_Null

  double precision function Node_Component_Surface_Density_Null(self,positionCylindrical,componentType,massType,haloLoaded)
    !% A null implementation of the surface density in a component. Always returns zero.
    implicit none
    class           (nodeComponent)              , intent(inout)           :: self
    integer                                      , intent(in   )           :: componentType      , massType
    double precision               , dimension(3), intent(in   )           :: positionCylindrical
    logical                                      , intent(in   ), optional :: haloLoaded

    Node_Component_Surface_Density_Null=0.0d0
    return
  end function Node_Component_Surface_Density_Null

  double precision function Node_Component_Potential_Null(self,radius,componentType,massType,haloLoaded)
    !% A null implementation of the gravitational potential in a component. Always returns zero.
    implicit none
    class           (nodeComponent), intent(inout)           :: self
    integer                        , intent(in   )           :: componentType, massType
    double precision               , intent(in   )           :: radius
    logical                        , intent(in   ), optional :: haloLoaded

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

end module Galacticus_Nodes
