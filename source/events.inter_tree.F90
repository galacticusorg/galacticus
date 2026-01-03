!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

!!{
Contains a module which handles inter-tree nodes events.
!!}

module Node_Events_Inter_Tree
  !!{
  Handles inter-tree node events.
  !!}
  use            :: Galacticus_Nodes, only : treeNode
  use, intrinsic :: ISO_C_Binding   , only : c_size_t
  use            :: Kind_Numbers    , only : kind_int8
  implicit none
  private
  public :: Node_Push_From_Tree, Node_Pull_From_Tree, Inter_Tree_Event_Post_Evolve

  type :: interTreeTransfer
     !!{
     Type used for transferring nodes between trees.
     !!}
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

  logical function Node_Push_From_Tree(event,node,deadlockStatus)
    !!{
    Push a node from the tree.
    !!}
    use :: Display                            , only : displayMessage               , verbosityLevelInfo
    use :: Error                              , only : Error_Report
    use :: Galacticus_Nodes                   , only : nodeComponentBasic           , nodeEvent                    , nodeEventBranchJumpInterTree, nodeEventSubhaloPromotionInterTree, &
          &                                            treeNode                     , treeNodeLinkedList
    use :: ISO_Varying_String                 , only : assignment(=)                , operator(//)                 , varying_string
    use :: Merger_Trees_Evolve_Deadlock_Status, only : deadlockStatusIsNotDeadlocked, enumerationDeadlockStatusType
    use :: String_Handling                    , only : operator(//)
    implicit none
    class    (nodeEvent                    ), intent(in   )          :: event
    type     (treeNode                     ), intent(inout), pointer :: node
    type     (enumerationDeadlockStatusType), intent(inout)          :: deadlockStatus
    type     (treeNode                     )               , pointer :: nodeMergee         , nodeWork
    class    (nodeComponentBasic           )               , pointer :: basic              , basicMergee
    type     (interTreeTransfer            )               , pointer :: waitListEntry
    type     (treeNodeLinkedList           )               , pointer :: nodeStack          , nodeNext, &
         &                                                              nodeNew
    integer  (c_size_t                     )                         :: splitForestUniqueID
    integer  (kind_int8                    )                         :: pairedNodeID
    type     (varying_string               )                         :: message
    character(len=12                       )                         :: label

    ! If this node is not yet a subhalo, do not perform the event.
    if (.not.node%isSatellite()) then
       Node_Push_From_Tree=.false.
       return
    end if
    ! Handle any mergees of the transferring node.
    if (associated(node%firstMergee)) then
       ! Any mergees must have reached the time of this halo before a transfer is allowed.
       basic      => node%basic      ()
       nodeMergee => node%firstMergee
       do while (associated(nodeMergee))
          basicMergee => nodeMergee%basic()
          if (basicMergee%time() < basic%time()) then
             Node_Push_From_Tree=.false.
             return
          end if
          nodeMergee => nodeMergee%siblingMergee
       end do
       ! Move mergees (and their mergees, etc.) to this node prior to inter-tree transfer.
       !! Build the initial stack.
       nodeStack  => null()
       nodeMergee => node%firstMergee
       do while (associated(nodeMergee))
          allocate(nodeNew)
          nodeNew   %node => nodeMergee
          nodeNew   %next => nodeStack
          nodeStack       => nodeNew
          nodeMergee      => nodeMergee%siblingMergee
       end do
       !! Process the stack.
       do while (associated(nodeStack))
          ! Pop a node from the stack.
          nodeWork => nodeStack%node
          nodeNext => nodeStack%next
          deallocate(nodeStack)
          nodeStack => nodeNext
          ! Push any mergees onto the stack.
          nodeMergee => nodeWork%firstMergee
          do while (associated(nodeMergee))
             allocate(nodeNew)
             nodeNew   %node => nodeMergee
             nodeNew   %next => nodeStack
             nodeStack       => nodeNew
             nodeMergee      => nodeMergee%siblingMergee
          end do
          ! Process the node.
          call nodeWork%removeFromHost()
          nodeWork%parent         => null()
          nodeWork%sibling        => node    %firstSatellite
          node    %firstSatellite => nodeWork
       end do
    end if
    ! Report on activity and extract identifiers.
    select type (event)
    type is (nodeEventSubhaloPromotionInterTree)
       write (label,'(f12.6)') event%time
       message='Satellite node ['
       message=message//node%index()//'] promoting to isolated node in another tree {event ID: '//event%ID//'} at time '//trim(label)//' Gyr'
       splitForestUniqueID=event%splitForestUniqueID
       pairedNodeID       =event%pairedNodeID
    type is (nodeEventBranchJumpInterTree      )
       write (label,'(f12.6)') event%time
       message='Satellite node ['
       message=message//node%index()//'] jumping branch to another tree {event ID: '//event%ID//'} at time '//trim(label)//' Gyr'
       splitForestUniqueID=event%splitForestUniqueID
       pairedNodeID       =event%pairedNodeID
     class default
        splitForestUniqueID=-1
        pairedNodeID       =-1
       call Error_Report('unknown event type'//{introspection:location})
    end select
    call displayMessage(message,verbosityLevelInfo)
    ! This is a subhalo jumping to another tree. Remove the node from its host, and explicitly nullify its parent pointer to
    ! prevent it being evolved any further until it is inserted into its new tree.
    call node%removeFromHost()
    node%parent => null()
    ! Put the node onto the inter-tree wait list.
    !$omp critical(interTreeWaitList)
    if (.not.associated(interTreeWaitList%next)) allocate(interTreeWaitList%next)
    waitListEntry => interTreeWaitList%next
    do while (associated(waitListEntry%node))
       if (.not.associated(waitListEntry%next)) allocate(waitListEntry%next)
       waitListEntry => waitListEntry%next
    end do
    waitListEntry%node                => node
    waitListEntry%next                => null()
    waitListEntry%splitForestUniqueID =  splitForestUniqueID
    waitListEntry%pairedNodeID        =  pairedNodeID
    waitListSize                      =  waitListSize       +1
    message='Inter-tree wait list now contains '
    message=message//waitListSize//' nodes'
    call displayMessage(message,verbosityLevelInfo)
    !$omp end critical(interTreeWaitList)
    ! Since we changed the tree, record that the tree is not deadlocked.
    deadlockStatus=deadlockStatusIsNotDeadlocked
    ! Record that the task was performed.
    Node_Push_From_Tree=.true.
    return
  end function Node_Push_From_Tree

  logical function Node_Pull_From_Tree(event,node,deadlockStatus)
    !!{
    Pull a node from the tree.
    !!}
    use :: Display                            , only : displayMessage            , verbosityLevelInfo           , verbosityLevelWarn
    use :: Error                              , only : Error_Report
    use :: Galacticus_Nodes                   , only : nodeComponentBasic        , nodeEvent                    , nodeEventBranchJumpInterTree, nodeEventSubhaloPromotionInterTree, &
          &                                            treeNode
    use :: ISO_Varying_String                 , only : assignment(=)             , operator(//)                 , varying_string
    use :: Merger_Trees_Evolve_Deadlock_Status, only : deadlockStatusIsDeadlocked, deadlockStatusIsNotDeadlocked, deadlockStatusIsSuspendable , enumerationDeadlockStatusType
    use :: String_Handling                    , only : operator(//)
    implicit none
    class           (nodeEvent                    ), intent(in   )          :: event
    type            (treeNode                     ), intent(inout), pointer :: node
    type            (enumerationDeadlockStatusType), intent(inout)          :: deadlockStatus
    type            (interTreeTransfer            )               , pointer :: waitListEntry              , waitListEntryPrevious
    type            (treeNode                     )               , pointer :: pullNode                   , satelliteNode        , &
         &                                                                     attachNode                 , hostNode             , &
         &                                                                     mergeeNode                 , mergeeNext
    class           (nodeComponentBasic           )               , pointer :: pullBasic                  , attachBasic
    class           (nodeEvent                    )               , pointer :: attachedEvent              , lastEvent            , &
         &                                                                     pairedEvent
    double precision                               , parameter              :: timeOffsetFractional=1.0d-6
    integer         (c_size_t                     )                         :: splitForestUniqueID        , pairedNodeID
    type            (varying_string               )                         :: message
    character       (len=12                       )                         :: label
    logical                                                                 :: isPrimary                  , timeMatchRequired    , &
         &                                                                     inSatellite

    ! Assume that the event cannot be performed by default.
    Node_Pull_From_Tree=.false.
    ! Assume that this event, if it can not be performed, at least makes the tree suspendable.
    if (deadlockStatus == deadlockStatusIsDeadlocked) deadlockStatus=deadlockStatusIsSuspendable
    ! Identify the event type.
    select type (event)
    type is (nodeEventSubhaloPromotionInterTree)
       write (label,'(f12.6)') event%time
       message='Searching for node to pull {promotion} to ['
       message=message//node%index()//'] {event ID: '//event%ID//'} at time '//trim(label)//' Gyr'
       splitForestUniqueID=event%splitForestUniqueID
       pairedNodeID       =event%pairedNodeID
       isPrimary          =event%isPrimary
       timeMatchRequired  =.true.
    type is (nodeEventBranchJumpInterTree      )
       write (label,'(f12.6)') event%time
       message='Searching for node to pull {branch jump} to descendeat ['
       message=message//node%index()//'] in host ['//node%parent%index()//'] {event ID: '//event%ID//'} at time '//trim(label)//' Gyr'
       splitForestUniqueID=event%splitForestUniqueID
       pairedNodeID       =event%pairedNodeID
       isPrimary          =event%isPrimary
       timeMatchRequired  =.false.
    class default
       splitForestUniqueID=-1
       pairedNodeID       =-1
       timeMatchRequired  =.false.
       isPrimary          =.false.
       call Error_Report('unknown event type'//{introspection:location})
    end select
    call displayMessage(message,verbosityLevelInfo)
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
             message=message//node%index()//':'//pullNode%index()//'] {event ID: '//event%ID//'; primary? '
             if (isPrimary) then
                message=message//"yes"
             else
                message=message//"no"
             end if
             message=message//'} at time '//trim(label)//' Gyr'
             call displayMessage(message,verbosityLevelInfo)
             ! Remove the node from the linked list.
             waitListEntryPrevious%next => waitListEntry        %next
             deallocate(waitListEntry)
             waitListEntry              => waitListEntryPrevious
             waitListSize               =  waitListSize              -1
             message='Inter-tree wait list now contains '
             message=message//waitListSize//' nodes'
             call displayMessage(message,verbosityLevelInfo)
             ! Attach the pulled node.
             pullNode     %sibling    => null()
             pullNode     %firstChild => null()
             pullNode     %hostTree   => node    %hostTree
             satelliteNode            => pullNode%firstSatellite
             do while (associated(satelliteNode))
                satelliteNode%hostTree => node         %hostTree
                satelliteNode          => satelliteNode%sibling
             end do
             if (isPrimary) then
                ! Handle primary progenitor cases.
                select type (event)
                type is (nodeEventSubhaloPromotionInterTree)
                   ! Node being jumped to should not be a satellite in this case.
                   if (node%isSatellite()) call Error_Report('inter-tree primary subhalo promotion, but jumped-to node is a satellite - unexpected behavior'//{introspection:location})
                   ! Pulled node is the primary progenitor and a subhalo promotion. It is being pulled to a node that is a clone of its parent. Replace the
                   ! clone with the pulled node.
                   if (associated(node%firstSatellite)) then
                      satelliteNode => node%firstSatellite
                      do while (associated(satelliteNode%sibling))
                         satelliteNode => satelliteNode%sibling
                      end do
                      satelliteNode%sibling => pullNode%firstSatellite
                   else
                      node%firstSatellite => pullNode%firstSatellite
                   end if
                   pullNode  %firstSatellite            => node    %firstSatellite
                   pullNode  %parent                    => node    %parent
                   pullNode  %sibling                   => node    %sibling
                   pullNode  %event                     => node    %event
                   node      %event                     => null()
                   node      %parent        %firstChild => pullNode
                   mergeeNode                           => node    %firstMergee
                   do while (associated(mergeeNode))
                      mergeeNode%mergeTarget   => pullNode
                      mergeeNext               => mergeeNode%siblingMergee
                      mergeeNode%siblingMergee => pullNode  %firstMergee
                      pullNode  %firstMergee   => mergeeNode
                      mergeeNode               => mergeeNext
                   end do
                   satelliteNode                        => pullNode%firstSatellite
                   do while (associated(satelliteNode))
                      satelliteNode%parent => pullNode
                      satelliteNode        => satelliteNode%sibling
                   end do
                   ! Remove this event from those attached to the pulled node, and reattach to the clone node so that it will be
                   ! destroyed when we destroy the clone node.
                   if (pullNode%event%ID == event%ID) then
                      node    %event      => pullNode%event
                      pullNode%event      => pullNode%event%next
                      node    %event%next => null()
                   else
                      lastEvent     => null()
                      attachedEvent => pullNode%event
                      do while (associated(attachedEvent))
                         if (attachedEvent%ID == event%ID) then
                            node %event      => attachedEvent
                            lastEvent      %next => attachedEvent%next
                            node %event%next => null()
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
                      message=message//" event ID="//event%id//char(10)
                      write (label,'(f12.6)') attachBasic%time()
                      message=message//"  node ID="//node%index()//"; time="//label//" Gyr"//char(10)
                      write (label,'(f12.6)')   pullBasic%time()
                      message=message//"  pull ID="//pullNode%index()//"; time="//label//" Gyr"
                      call Error_Report(message//{introspection:location})
                   end if
                   call pullBasic%timeSet(pullBasic%time()*(1.0d0-timeOffsetFractional))
                   ! Destroy the cloned node.
                   call node %destroy()
                   deallocate(node)
                type is (nodeEventBranchJumpInterTree      )
                   ! Pulled node is the primary progenitor, but is a branch jump event. If the target node has a timeLastIsolated
                   ! prior to the current time, that indicates that it had one or more progenitors (in this tree) which were
                   ! non-primary (since the pulled node is primary). One of those progenitors would be treated as primary. In that
                   ! case, we must retain that progenitor and add our pulled node as a satellite node in the host halo. Otherwise,
                   ! simply attach it as the primary progenitor of the node it jumps to.
                   attachBasic => node%basic()
                   if (attachBasic%timeLastIsolated() < attachBasic%time()) then
                      pullNode%parent                => node%parent
                      pullNode%sibling               => node%parent%firstSatellite
                      pullNode%parent%firstSatellite => pullNode
                      ! Allow any necessary manipulation of the nodes.
                      !![
		      <eventHook name="interTreeSatelliteInsert">
			<import>
			  <module name="Galacticus_Nodes" symbols="treeNode"/>
			</import>
			<interface>
			  type(treeNode), intent(inout), pointer :: pullNode, node
			</interface>
			<callWith>pullNode,node</callWith>
		      </eventHook>
                      !!]
                   else
                      ! Attach pulled node as the primary progenitor of the target node as it is the primary (and only progenitor).
                      pullNode%parent     => node
                      pullNode%sibling    => null()
                      node%firstChild => pullNode
                      ! Reset the ID of the descendant to that of the progenitor to mimic what would
                      ! have occurred if trees were processed unsplit. We also reset the basic mass and time last isolated for the same
                      ! reason.
                      !![
		      <eventHook name="interTreeSatelliteAttach">
			<import>
			  <module name="Galacticus_Nodes" symbols="treeNode"/>
			</import>
			<interface>
			  type(treeNode), intent(inout), pointer :: pullNode
			</interface>
			<callWith>pullNode</callWith>
		      </eventHook>
                      !!]
                      pullBasic   => pullNode%basic()
                      attachBasic => node%basic()
                      call node   %           indexSet(pullNode %           index())
                      call attachBasic%            massSet(pullBasic%            mass())
                      call attachBasic%timeLastIsolatedSet(pullBasic%timeLastIsolated())
                   end if
                end select
             else
                ! Pulled node is not the primary progenitor, so our node is actually a satellite in the node pulled to (or its
                ! host should it be a subhalo).
                inSatellite =  .false.
                attachNode  => node
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
                      call displayMessage(message,verbosityLevelWarn)
                      warningNestedHierarchyIssued=.true.
                   end if
                   do while (attachNode%isSatellite())
                      attachNode => attachNode%parent
                   end do
                else
                   ! Node being pulled to is not a satellite, so attach as a satellite in it directly.
                   attachNode => node
                end if
                pullNode   %sibling        => attachNode%firstSatellite
                attachNode %firstSatellite => pullNode
                pullNode   %parent         => attachNode
                pullBasic                  => pullNode                 %basic()
                if (associated(attachNode%parent)) then
                   attachBasic => attachNode%parent%basic()
                   if (pullBasic%time() > attachBasic%time()) then
                      message='pulled node exceeds node in time [non-primary]:'//char(10)
                      message=message//" event ID="//event%id//char(10)
                      write (label,'(f12.6)') attachBasic%time()
                      message=message//"  node ID="//attachNode%index()//"; time="//label//" Gyr"//char(10)
                      write (label,'(f12.6)')   pullBasic%time()
                      message=message//"  pull ID="//pullNode  %index()//"; time="//label//" Gyr"
                      call Error_Report(message//{introspection:location})
                   end if
                end if
                ! Assign a merging time to the new satellite if possible.
                select type (event)
                type is (nodeEventSubhaloPromotionInterTree)
                   attachBasic => attachNode%basic()
                   if (associated(event%mergeTimeSet)) call event%mergeTimeSet(event%creator,pullNode,attachNode)
                type is (nodeEventBranchJumpInterTree      )
                   attachBasic => attachNode%basic()
                   if (associated(event%mergeTimeSet)) call event%mergeTimeSet(event%creator,pullNode,attachNode)
                class default
                   call Error_Report('non-primary jump should be inter-tree branch jump'//{introspection:location})
                end select
             end if
             ! If the node or its parent are now satellites, and have their own satellites, transfer these satellites to the new
             ! host node.
             if ((pullNode%isSatellite().or.pullNode%parent%isSatellite()).and.associated(pullNode%firstSatellite)) then
                if (pullNode%isSatellite()) then
                   hostNode => pullNode%parent
                else if (pullNode%parent%isSatellite()) then
                   hostNode => pullNode%parent%parent
                else
                   hostNode => null()
                   call Error_Report('neither node nor parent are satellites - this should not happen'//{introspection:location})
                end if
                satelliteNode          => pullNode%firstSatellite
                satelliteNode%parent   => hostNode
                satelliteNode%hostTree => hostNode%hostTree
                do while (associated(satelliteNode%sibling))
                   satelliteNode          => satelliteNode%sibling
                   satelliteNode%parent   => hostNode
                   satelliteNode%hostTree => hostNode%hostTree
                end do
                satelliteNode%sibling   => hostNode%firstSatellite
                hostNode%firstSatellite => pullNode%firstSatellite
                pullNode%firstSatellite => null()
             end if
             ! Allow any postprocessing of the inter-tree transfer event that may be necessary.
             !![
	     <eventHook name="interTreePostProcess">
	       <import>
		 <module name="Galacticus_Nodes" symbols="treeNode"/>
	       </import>
	       <interface>
		 type(treeNode), intent(inout), pointer :: pullNode
	       </interface>
	       <callWith>pullNode</callWith>
	     </eventHook>
             !!]
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

  !![
  <universePostEvolveTask>
   <unitName>Inter_Tree_Event_Post_Evolve</unitName>
  </universePostEvolveTask>
  !!]
  subroutine Inter_Tree_Event_Post_Evolve()
    !!{
    Check that the inter-tree transfer list is empty after universe evolution.
    !!}
    use :: Display           , only : displayIndent, displayMessage, displayUnindent
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : assignment(=), varying_string
    use :: String_Handling   , only : operator(//)
    implicit none
    type(interTreeTransfer), pointer :: waitListEntry
    type(varying_string   )          :: message

    if (waitListSize > 0) then
       call displayIndent('Nodes in inter-tree transfer wait list')
       waitListEntry => interTreeWaitList%next
       do while (associated(waitListEntry).and.associated(waitListEntry%node))
          message='Node ID: '
          message=message//waitListEntry%node%index()
          call displayMessage(message)
          waitListEntry => waitListEntry%next
       end do
       call displayUnindent('done')
       call Error_Report('nodes remain in the inter-tree transfer wait list - see preceeding report'//{introspection:location})
    end if
    return
  end subroutine Inter_Tree_Event_Post_Evolve

end module Node_Events_Inter_Tree
