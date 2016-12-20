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

!% Contains a module which handles inter-tree nodes events.

module Node_Events_Inter_Tree
  !% Handles inter-tree node events.
  use, intrinsic :: ISO_C_Binding
  use               Kind_Numbers
  use               Galacticus_Nodes
  implicit none
  private
  public :: Node_Push_From_Tree, Node_Pull_From_Tree, Inter_Tree_Event_Post_Evolve

  type :: interTreeTransfer
     !% Type used for transfering nodes between trees.
     integer(c_size_t         )          :: splitForestUniqueID
     integer(kind_int8        )          :: pairedNodeID
     type   (treeNode         ), pointer :: node                => null()
     type   (interTreeTransfer), pointer :: next                => null()
  end type interTreeTransfer

  ! Head of linked list for inter-tree transfers.
  type   (interTreeTransfer), target :: interTreeWaitList
  integer(c_size_t         )         :: waitListSize                =0

  ! Record of warnings issued.
  logical                            :: warningNestedHierarchyIssued=.false.
  
contains

  logical function Node_Push_From_Tree(thisEvent,thisNode,deadlockStatus)
    !% Push a node from the tree.
    use ISO_Varying_String
    use Galacticus_Display
    use Galacticus_Error
    use String_Handling
    use Merger_Trees_Evolve_Deadlock_Status
    implicit none
    class    (nodeEvent        ), intent(in   )          :: thisEvent
    type     (treeNode         ), intent(inout), pointer :: thisNode
    integer                     , intent(inout)          :: deadlockStatus
    type     (interTreeTransfer)               , pointer :: waitListEntry
    integer  (c_size_t         )                         :: splitForestUniqueID
    integer  (kind_int8        )                         :: pairedNodeID
    type     (varying_string   )                         :: message
    character(len=12           )                         :: label

    ! If this node is not yet a subhalo, do not perform the event.
    if (.not.thisNode%isSatellite()) then
       Node_Push_From_Tree=.false.
       return
    end if    
    ! Report on activity and extract identifiers.
    select type (thisEvent)
    type is (nodeEventSubhaloPromotionInterTree)
       write (label,'(f12.6)') thisEvent%time
       message='Satellite node ['
       message=message//thisNode%index()//'] promoting to isolated node in another tree {event ID: '//thisEvent%ID//'} at time '//trim(label)//' Gyr'
       splitForestUniqueID=thisEvent%splitForestUniqueID
       pairedNodeID       =thisEvent%pairedNodeID
    type is (nodeEventBranchJumpInterTree      )
       write (label,'(f12.6)') thisEvent%time
       message='Satellite node ['
       message=message//thisNode%index()//'] jumping branch to another tree {event ID: '//thisEvent%ID//'} at time '//trim(label)//' Gyr'
       splitForestUniqueID=thisEvent%splitForestUniqueID
       pairedNodeID       =thisEvent%pairedNodeID
     class default
        splitForestUniqueID=-1
        pairedNodeID       =-1
       call Galacticus_Error_Report('Node_Push_From_Tree','unknown event type')
    end select
    call Galacticus_Display_Message(message,verbosityInfo)
    ! This is a subhalo jumping to another in another tree. Remove the node from its host.
    call thisNode%removeFromHost()
    ! Put the node onto the inter-tree wait list.
    !$omp critical(interTreeWaitList)
    if (.not.associated(interTreeWaitList%next)) allocate(interTreeWaitList%next)
    waitListEntry => interTreeWaitList%next
    do while (associated(waitListEntry%node))
       if (.not.associated(waitListEntry%next)) allocate(waitListEntry%next)
       waitListEntry => waitListEntry%next
    end do
    waitListEntry%node                => thisNode
    waitListEntry%next                => null()
    waitListEntry%splitForestUniqueID =  splitForestUniqueID
    waitListEntry%pairedNodeID        =  pairedNodeID
    waitListSize                      =  waitListSize       +1
    message='Inter-tree wait list now contains '
    message=message//waitListSize//' nodes'
    call Galacticus_Display_Message(message,verbosityInfo)
    !$omp end critical(interTreeWaitList)
    ! Since we changed the tree, record that the tree is not deadlocked.
    deadlockStatus=deadlockStatusIsNotDeadlocked
    ! Record that the task was performed.
    Node_Push_From_Tree=.true.
    return
  end function Node_Push_From_Tree

  logical function Node_Pull_From_Tree(thisEvent,thisNode,deadlockStatus)
    !% Pull a node from the tree.
    use ISO_Varying_String
    use Galacticus_Display
    use String_Handling
    use Galacticus_Error
    use Merger_Trees_Evolve_Deadlock_Status
    !# <include directive="interTreeSatelliteAttach" type="moduleUse">
    include 'events.inter_tree.satellite_attach.modules.inc'
    !# </include>
    !# <include directive="interTreeSatelliteInsert" type="moduleUse">
    include 'events.inter_tree.satellite_insert.modules.inc'
    !# </include>
    implicit none
    class           (nodeEvent         ), intent(in   )          :: thisEvent
    type            (treeNode          ), intent(inout), pointer :: thisNode
    integer                             , intent(inout)          :: deadlockStatus
    type            (interTreeTransfer )               , pointer :: waitListEntry              , waitListEntryPrevious
    type            (treeNode          )               , pointer :: pullNode                   , satelliteNode        , &
         &                                                          attachNode
    class           (nodeComponentBasic)               , pointer :: pullBasic                  , attachBasic
    class           (nodeEvent         )               , pointer :: attachedEvent              , lastEvent            , &
         &                                                          pairedEvent
    double precision                    , parameter              :: timeOffsetFractional=1.0d-6
    integer         (c_size_t          )                         :: splitForestUniqueID        , pairedNodeID
    type            (varying_string    )                         :: message
    character       (len=12            )                         :: label
    logical                                                      :: isPrimary                  , timeMatchRequired    , &
         &                                                          inSatellite

    ! Assume that the event cannot be performed by default.
    Node_Pull_From_Tree=.false.
    ! Assume that this event, if it can not be performed, at least makes the tree suspendable.
    if (deadlockStatus == deadlockStatusIsDeadlocked) deadlockStatus=deadlockStatusIsSuspendable
    ! Identify the event type.
    select type (thisEvent)
    type is (nodeEventSubhaloPromotionInterTree)
       write (label,'(f12.6)') thisEvent%time
       message='Searching for node to pull {promotion} to ['
       message=message//thisNode%index()//'] {event ID: '//thisEvent%ID//'} at time '//trim(label)//' Gyr'
       splitForestUniqueID=thisEvent%splitForestUniqueID
       pairedNodeID       =thisEvent%pairedNodeID
       isPrimary          =thisEvent%isPrimary
       timeMatchRequired  =.true.
    type is (nodeEventBranchJumpInterTree      )
       write (label,'(f12.6)') thisEvent%time
       message='Searching for node to pull {branch jump} to descendent ['
       message=message//thisNode%index()//'] in host ['//thisNode%parent%index()//'] {event ID: '//thisEvent%ID//'} at time '//trim(label)//' Gyr'
       splitForestUniqueID=thisEvent%splitForestUniqueID
       pairedNodeID       =thisEvent%pairedNodeID
       isPrimary          =thisEvent%isPrimary
       timeMatchRequired  =.false.
    class default
       splitForestUniqueID=-1
       pairedNodeID       =-1
       timeMatchRequired  =.false.
       isPrimary          =.false.
       call Galacticus_Error_Report('Node_Pull_From_Tree','unknown event type')
    end select
    call Galacticus_Display_Message(message,verbosityInfo)
    ! Search for the node to be pulled in the inter-tree wait list.
    !$omp critical(interTreeWaitList)
    if (associated(interTreeWaitList%next)) then
       waitListEntryPrevious => interTreeWaitList
       waitListEntry         => interTreeWaitList%next
       do while (associated(waitListEntry).and.associated(waitListEntry%node))
          if (waitListEntry%splitForestUniqueID == splitForestUniqueID .and. waitListEntry%pairedNodeID == pairedNodeID) then
             ! Grab the node from the waitlist.
             pullNode           => waitListEntry%node
             waitListEntry%node => null()
             message='Found node to pull to ['
             message=message//thisNode%index()//':'//pullNode%index()//'] {event ID: '//thisEvent%ID//'; primary? '
             if (isPrimary) then
                message=message//"yes"
             else
                message=message//"no"
             end if
             message=message//'} at time '//trim(label)//' Gyr'
             call Galacticus_Display_Message(message,verbosityInfo)
             ! Remove the node from the linked list.
             waitListEntryPrevious%next => waitListEntry        %next
             deallocate(waitListEntry)
             waitListEntry              => waitListEntryPrevious
             waitListSize               =  waitListSize              -1
             message='Inter-tree wait list now contains '
             message=message//waitListSize//' nodes'
             call Galacticus_Display_Message(message,verbosityInfo)
             ! Attach the pulled node.
             pullNode%sibling        => null()
             pullNode%firstChild     => null()
             pullNode%firstSatellite => null()
             pullNode%hostTree       => thisNode%hostTree
             if (isPrimary) then
                ! Handle primary progenitor cases.
                select type (thisEvent)
                type is (nodeEventSubhaloPromotionInterTree)                
                   ! Node being jumped to should not be a satellite in this case.
                   if (thisNode%isSatellite()) call Galacticus_Error_Report('Node_Pull_From_Tree','inter-tree primary subhalo promotion, but jumped-to node is a satellite - unexpected behavior')
                   ! Pulled node is the primary progenitor and a subhalo promotion. It is being pulled to a node that is a clone of its parent. Replace the
                   ! clone with the pulled node.                
                   pullNode     %parent                    => thisNode%parent
                   pullNode     %sibling                   => thisNode%sibling
                   pullNode     %firstSatellite            => thisNode%firstSatellite
                   pullNode     %event                     => thisNode%event
                   thisNode     %event                     => null()
                   thisNode     %parent        %firstChild => pullNode
                   satelliteNode                           => pullNode%firstSatellite
                   do while (associated(satelliteNode))
                      satelliteNode%parent => pullNode
                      satelliteNode        => satelliteNode%sibling
                   end do
                   ! Remove this event from those attached to the pulled node, and reattach to the clone node so that it will be
                   ! destroyed when we destroy the clone node.
                   if (pullNode%event%ID == thisEvent%ID) then
                      thisNode%event      => pullNode%event
                      pullNode%event      => pullNode%event%next
                      thisNode%event%next => null()
                   else
                      lastEvent     => null()
                      attachedEvent => pullNode%event
                      do while (associated(attachedEvent))
                         if (attachedEvent%ID == thisEvent%ID) then
                            thisNode %event      => attachedEvent
                            lastEvent      %next => attachedEvent%next
                            thisNode %event%next => null()
                            exit
                         end if
                         lastEvent     => attachedEvent
                         attachedEvent => attachedEvent%next
                      end do
                   end if
                   ! For any attached events, check if they have a paired event and redirect that paired event to this node.
                   attachedEvent => pullNode%event
                   do while (associated(attachedEvent))
                      if (associated(attachedEvent%node)) then
                         ! A paired event may exist. Search for it.
                         pairedEvent => attachedEvent%node%event
                         do while (associated(pairedEvent))
                            if (pairedEvent%ID == attachedEvent%ID) then
                               ! A paired event exists. Redirect the node associated with this paired event to our new node.
                               pairedEvent%node => pullNode
                            end if
                            pairedEvent => pairedEvent%next
                         end do
                      end if
                      attachedEvent => attachedEvent%next
                   end do
                   pullBasic   => pullNode       %basic()
                   attachBasic => pullNode%parent%basic()
                   if (timeMatchRequired .and. pullBasic%time() /= attachBasic%time()) then
                      message='pulled node does not match in time [primary]:'//char(10)
                      message=message//" event ID="//thisEvent%id//char(10)
                      write (label,'(f12.6)') attachBasic%time()
                      message=message//"  node ID="//thisNode%index()//"; time="//label//" Gyr"//char(10)
                      write (label,'(f12.6)')   pullBasic%time()
                      message=message//"  pull ID="//pullNode%index()//"; time="//label//" Gyr"
                      call Galacticus_Error_Report('Node_Pull_From_Tree',message)
                   end if
                   call pullBasic%timeSet(pullBasic%time()*(1.0d0-timeOffsetFractional))
                   ! Destroy the cloned node.
                   call thisNode %destroy()
                   deallocate(thisNode)
                type is (nodeEventBranchJumpInterTree      )
                   ! Pulled node is the primary progenitor, but is a branch jump event. If the target node has a timeLastIsolated
                   ! prior to the current time, that indicates that it had one or more progenitors (in this tree) which where
                   ! non-primary (since the pulled node is primary). One of those progenitors would be treated as primary. In that
                   ! case, we must retain that progenitor and add our pulled node as a satellite node in the host halo. Otherwise,
                   ! simply attach it as the primary progenitor of the node is jumps to.
                   attachBasic => thisNode%basic()
                   if (attachBasic%timeLastIsolated() < attachBasic%time()) then   
                      pullNode%parent                => thisNode%parent
                      pullNode%sibling               => thisNode%parent%firstSatellite
                      pullNode%parent%firstSatellite => pullNode
                      ! Allow any necessary manipulation of the nodes.
                      !# <include directive="interTreeSatelliteInsert" type="functionCall" functionType="void">
                      !#  <functionArgs>pullNode,thisNode</functionArgs>
                      include 'events.inter_tree.satellite_insert.inc'
                      !# </include>  
                   else
                      ! Attach pulled node as the primary progenitor of the target node as it is the primary (and only progenitor).
                      pullNode%parent     => thisNode
                      pullNode%sibling    => null()
                      thisNode%firstChild => pullNode
                      ! Reset the ID of the descendent to that of the progenitor to mimic what would
                      ! have occurred if trees were processed unsplit. We also reset the basic mass and time last isolated for the same
                      ! reason.
                      !# <include directive="interTreeSatelliteAttach" type="functionCall" functionType="void">
                      !#  <functionArgs>pullNode</functionArgs>
                      include 'events.inter_tree.satellite_attach.inc'
                      !# </include>
                      pullBasic   => pullNode%basic()
                      attachBasic => thisNode%basic()
                      call thisNode   %           indexSet(pullNode %           index())
                      call attachBasic%            massSet(pullBasic%            mass())
                      call attachBasic%timeLastIsolatedSet(pullBasic%timeLastIsolated())                   
                   end if
                end select
             else
                ! Pulled node is not the primary progenitor, so our node is actually a satellite in the node pulled to (or its
                ! host should it be a subhalo).
                inSatellite =  .false.
                attachNode  => thisNode
                do while(associated(attachNode))
                   if (attachNode%isSatellite()) then
                      inSatellite=.true.
                      exit
                   end if
                   attachNode => attachNode%parent
                end do
                if (inSatellite) then
                   ! Node being pulled is a satellite. Find its isolated host and attach as a satellite in there.
                   if (.not.warningNestedHierarchyIssued) then
                      message='nested hierarchy detected [node '
                      message=message//pullNode%index()//']'
                      message=message//char(10)//'ignoring as not currently supported'
                      message=message//char(10)//'warning will not be issued again'
                      call Galacticus_Display_Message(message,verbosityWarn)
                      warningNestedHierarchyIssued=.true.
                   end if
                   do while (attachNode%isSatellite())
                      attachNode => attachNode%parent
                   end do
                else
                   ! Node being pulled to is not a satellite, so attach as a satellite in it directly.
                   attachNode => thisNode
                end if
                pullNode   %sibling        => attachNode%firstSatellite
                attachNode %firstSatellite => pullNode
                pullNode   %parent         => attachNode
                pullBasic                  => pullNode                 %basic()
                if (associated(attachNode%parent)) then
                   attachBasic => attachNode%parent%basic()                
                   if (pullBasic%time() > attachBasic%time()) then
                      message='pulled node exceeds node in time [non-primary]:'//char(10)
                      message=message//" event ID="//thisEvent%id//char(10)
                      write (label,'(f12.6)') attachBasic%time()
                      message=message//"  node ID="//attachNode%index()//"; time="//label//" Gyr"//char(10)
                      write (label,'(f12.6)')   pullBasic%time()
                      message=message//"  pull ID="//pullNode  %index()//"; time="//label//" Gyr"
                      call Galacticus_Error_Report('Node_Pull_From_Tree',message)
                   end if
                end if
             end if
             ! Record that the event was performed, and set the deadlock status to not deadlocked since we changed the tree.
             deadlockStatus     =deadlockStatusIsNotDeadlocked
             Node_Pull_From_Tree=.true.
             exit
          end if
          waitListEntryPrevious => waitListEntry
          waitListEntry         => waitListEntry%next
       end do
    end if           
    !$omp end critical(interTreeWaitList)
    return
  end function Node_Pull_From_Tree

  !# <universePostEvolveTask>
  !#  <unitName>Inter_Tree_Event_Post_Evolve</unitName>
  !# </universePostEvolveTask>
  subroutine Inter_Tree_Event_Post_Evolve()
    !% Check that the inter-tree transfer list is empty after universe evolution.
    use ISO_Varying_String
    use Galacticus_Error
    use Galacticus_Display
    use String_Handling
    implicit none
    type(interTreeTransfer), pointer :: waitListEntry
    type(varying_string   )          :: message

    if (waitListSize > 0) then
       call Galacticus_Display_Indent('Nodes in inter-tree transfer wait list')
       waitListEntry => interTreeWaitList%next
       do while (associated(waitListEntry).and.associated(waitListEntry%node))
          message='Node ID: '
          message=message//waitListEntry%node%index()
          call Galacticus_Display_Message(message)
          waitListEntry => waitListEntry%next
       end do
       call Galacticus_Display_Unindent('done')
       call Galacticus_Error_Report('Inter_Tree_Event_Post_Evolve','nodes remain in the inter-tree transfer wait list - see preceeding report')
    end if
    return
  end subroutine Inter_Tree_Event_Post_Evolve
  
end module Node_Events_Inter_Tree
