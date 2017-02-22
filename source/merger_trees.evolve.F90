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

!% Contains a module which implements evolution of merger trees.

module Merger_Trees_Evolve
  !% Implements evolution of merger trees.
  use Galacticus_Nodes
  use ISO_Varying_String
  use Kind_Numbers
  implicit none
  private
  public :: Merger_Tree_Evolve_To

  ! Flag indicating if evolver routine has been initialized.
  logical          :: mergerTreeEvolveToInitialized      =.false.
  logical          :: mergerTreeEvolveToThreadInitialized=.false.
  !$omp threadprivate(mergerTreeEvolveToThreadInitialized)

  ! Flag indicating whether or not to fail for trees which do not exist at the final output time.
  logical          :: allTreesExistAtFinalTime

  ! Flag indicating whether to dump merger tree structure after each evolutionary step.
  logical          :: mergerTreesDumpStructure

  ! Variables which limit extent to which satellites can evolve past their parent.
  logical          :: evolveToTimeInitialized      =.false.
  double precision :: timestepHostAbsolute                 , timestepHostRelative

  ! Structure used to store list of nodes for deadlock reporting.
  type :: deadlockList
     type   (deadlockList  ), pointer :: next     =>null()
     type   (treeNode      ), pointer :: nodeLock         , node
     integer(kind=kind_int8)          :: treeIndex
     type   (varying_string)          :: lockType
  end type deadlockList
  type(deadlockList), pointer :: deadlockHeadNode => null()
  !$omp threadprivate(deadlockHeadNode)
  
contains

  subroutine Merger_Tree_Evolve_To(thisTree,endTime,treeDidEvolve,suspendTree,deadlockReporting)
    !% Evolves all properties of a merger tree to the specified time.
    use Merger_Trees_Evolve_Node
    use Merger_Trees_Evolve_Timesteps_Template
    use Merger_Trees_Initialize
    use Merger_Trees_Dump
    use Galacticus_Error
    use Galacticus_Display
    use Input_Parameters
    use String_Handling
    use Merger_Trees_Evolve_Deadlock_Status
    !# <include directive="mergerTreeEvolveThreadInitialize" type="moduleUse">
    include 'merger_trees.evolve.threadInitialize.moduleUse.inc'
    !# </include>
    implicit none
    type            (mergerTree                   )           , target , intent(inout) :: tree
    double precision                                                   , intent(in   ) :: endTime
    logical                                                            , intent(  out) :: treeDidEvolve                     , suspendTree
    logical                                                            , intent(in   ) :: deadlockReporting
    class           (nodeEvent                    )           , pointer                :: event
    double precision                               , parameter                         :: timeTolerance              =1.0d-5
    double precision                               , parameter                         :: largeTime                  =1.0d10
    procedure       (interruptTask                )           , pointer                :: interruptProcedure
    procedure       (End_Of_Timestep_Task_Template)           , pointer                :: End_Of_Timestep_Task
    integer                                        , parameter                         :: verbosityLevel             =3
    class           (nodeComponentBasic           )           , pointer                :: basicBase                         , basicParent      , &
         &                                                                                basic
    type            (mergerTree                   )           , pointer                :: currentTree
    integer                                                                            :: deadlockStatus                    , nodesEvolvedCount, &
         &                                                                                nodesTotalCount                   , treeWalkCount    , &
         &                                                                                treeWalkCountPreviousOutput
    double precision                                                                   :: earliestTimeInTree                , endTimeThisNode  , &
         &                                                                                finalTimeInTree
    logical                                                                            :: didEvolve                         , interrupted
    character       (len=24                       )                                    :: label
    character       (len=35                       )                                    :: message
    type            (varying_string               )                                    :: lockType                          , vMessage
    logical                                                                            :: anyTreeExistsAtOutputTime         , hasIntertreeEvent   , &
         &                                                                                hasParent                         , treeLimited
    
    ! Check if this routine is initialized.
    if (.not.mergerTreeEvolveToInitialized) then
       !$omp critical (Merger_Tree_Evolve_To_Initialize)
       if (.not.mergerTreeEvolveToInitialized) then
          ! Read parameters.
          !@ <inputParameter>
          !@   <name>allTreesExistAtFinalTime</name>
          !@   <defaultValue>true</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not all merger trees are expected to exist at the final requested output time. If set to false,
          !@     then trees which finish before a given output time will be ignored.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('allTreesExistAtFinalTime',allTreesExistAtFinalTime,defaultValue=.true.)
          !@ <inputParameter>
          !@   <name>mergerTreesDumpStructure</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether merger tree structure should be dumped to a \href{http://www.graphviz.org/}{\normalfont \scshape dot} file.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreesDumpStructure',mergerTreesDumpStructure,defaultValue=.false.)

          ! Flag that this routine is now initialized.
          mergerTreeEvolveToInitialized=.true.
       end if
       !$omp end critical (Merger_Tree_Evolve_To_Initialize)
    end if
    
    ! Call routines to perform initializations which must occur for all threads if run in parallel.
    if (.not.mergerTreeEvolveToThreadInitialized) then
       !# <include directive="mergerTreeEvolveThreadInitialize" type="functionCall" functionType="void">
       include 'merger_trees.evolve.threadInitialize.inc'
       !# </include>
       mergerTreeEvolveToThreadInitialized=.true.
    end if

    ! Iterate through all trees.
    suspendTree               =  .false.
    anyTreeExistsAtOutputTime =  .false.
    treeDidEvolve             =  .false.
    currentTree               => thisTree
    do while (associated(currentTree))
       ! Skip empty trees.
       if (associated(currentTree%baseNode)) then
          ! Initialize the tree if necessary.
          call Merger_Tree_Initialize(currentTree,endTime)
          ! Check that the output time is not after the end time of this tree.
          basicBase => currentTree%baseNode%basic()
          if (endTime > basicBase%time()) then
             ! Final time is exceeded. Check if by a significant factor.
             if (endTime > basicBase%time()*(1.0d0+timeTolerance)) then
                ! Exceeded by a significant factor - report an error. Check if such behavior is expected.
                if (allTreesExistAtFinalTime) then
                   ! It is not, write an error and exit.
                   vMessage='requested time exceeds the final time in the tree'//char(10)
                   vMessage=vMessage//' HELP: If you expect that not all trees will exist at the latest requested'//char(10)
                   vMessage=vMessage//'       output time (this can happen when using trees extracted from N-body'//char(10)
                   vMessage=vMessage//'       simulations for example) set the following in your input parameter file:'//char(10)//char(10)
                   vMessage=vMessage//'         <allTreesExistAtFinalTime value="false" />'//char(10)
                   call Galacticus_Error_Report('Merger_Tree_Evolve_To',vMessage)
                end if
             else
                ! Not exceeded by a significant factor (can happen due to approximation errors). Unless there is an event
                ! associated with this node at the current time, simply reset to actual time requested.
                event => currentTree%baseNode%event
                do while (associated(event))
                   if (event%time == basicBase%time()) then
                      vMessage=          'requested time exceeds the final time in the tree by a small factor'  //char(10)
                      vMessage=vMessage//'refusing to adjust the final time in the tree due to associated event'//char(10)
                      write (label,'(e24.16)') endTime
                      vMessage=vMessage//'  requested time: '//trim(label)//' Gyr'//char(10)
                      write (label,'(e24.16)') basicBase%time()
                      vMessage=vMessage//'      final time: '//trim(label)//' Gyr'//char(10)
                      write (label,'(e24.16)') event%time
                      vMessage=vMessage//'      event time: '//trim(label)//' Gyr'//char(10)
                      vMessage=vMessage//'      event ID  : '//event%ID           //char(10)
                      vMessage=vMessage//' HELP: if you are reading merger trees from file and are attempting to'//char(10)
                      vMessage=vMessage//'       output at a "snapshot time" consider setting:'                  //char(10)
                      vMessage=vMessage//'           <mergerTreeReadOutputTimeSnapTolerance value="1.0e-3"/>'    //char(10)
                      vMessage=vMessage//'       or similar in your parameter file to ensure that nodes exist'   //char(10)
                      vMessage=vMessage//'       precisely at the output times you request'
                      call Galacticus_Error_Report('Merger_Tree_Evolve_To',vMessage)
                   end if
                   event => event%next
                end do
                call basicBase%timeSet(endTime)
                anyTreeExistsAtOutputTime=.true.
             end if
          else
             anyTreeExistsAtOutputTime=.true.
          end if
       end if
       ! Move to the next tree.
       currentTree => currentTree%nextTree
    end do

    ! If none of these trees exist at the output time, check if they contain any inter-tree events. If they do, we need to evolve
    ! the tree anyway, as it interacts with another tree that may exist at the output time. Otherwise, we can ignore this tree.
    if (.not.anyTreeExistsAtOutputTime) then
       hasInterTreeEvent=.false.
       ! Iterate over trees.
       currentTree => thisTree
       do while (associated(currentTree).and..not.hasIntertreeEvent)
          ! Iterate over nodes.
          thisNode => currentTree%baseNode
          do while (associated(thisNode).and..not.hasIntertreeEvent)                   
             ! Iterate over events.
             event => thisNode%event
             do while (associated(event).and..not.hasIntertreeEvent)
                select type (event)
                type is (nodeEventSubhaloPromotionInterTree)
                   hasIntertreeEvent=.true.
                type is (nodeEventBranchJumpInterTree      )
                   hasIntertreeEvent=.true.
                end select
                event => event%next
             end do
             thisNode => thisNode%walkTreeWithSatellites()
          end do
          currentTree => currentTree%nextTree
       end do
       if (.not.hasInterTreeEvent) then
          ! Mark the tree as evolved here, as the only reason that we did not evolve it was the given time target.
          treeDidEvolve=.true.
          return
       end if
    end if

    ! Outer loop: This causes the tree to be repeatedly walked and evolved until it has been evolved all the way to the specified
    ! end time. We stop when no nodes were evolved, which indicates that no further evolution is possible.
    didEvolve                  =.true.
    treeWalkCount              =0
    treeWalkCountPreviousOutput=0
    outerLoop: do while (didEvolve) ! Keep looping through the tree until we make a pass during which no nodes were evolved.
       
       ! Flag that no nodes have been evolved yet.
       didEvolve=.false.

       ! Increment tree walks counter.
       treeWalkCount=treeWalkCount+1

       ! Reset tree progress variables.
       nodesEvolvedCount =0
       nodesTotalCount   =0
       earliestTimeInTree=largeTime

       ! Set the deadlock status to deadlocked initially, unless we have been specifically asked for deadlock reporting via a
       ! function argument, in which case set to reporting status.
       deadlockStatus=deadlockStatusIsDeadlocked
       if (deadlockReporting) deadlockStatus=deadlockStatusIsReporting
       ! Enter loop for deadlock reporting.
       deadlock : do while (deadlockStatus /= deadlockStatusIsNotDeadlocked)

          ! Post a deadlocking message.
          if (deadlockStatus == deadlockStatusIsReporting) call Galacticus_Display_Indent("Deadlock report follows")
          
          ! Iterate through all trees.
          currentTree => tree
          treesLoop: do while (associated(currentTree))
             
             ! Skip empty trees.
             if (associated(currentTree%baseNode)) then

                ! Find the final time in this tree.
                thisNode           => currentTree%baseNode
                thisBasicComponent => thisNode%basic()
                finalTimeInTree    =  thisBasicComponent%time()

                ! Report on current tree if deadlocked.
                if (deadlockStatus == deadlockStatusIsReporting) then
                   vMessage="tree "
                   vMessage=vMessage//currentTree%index
                   call Galacticus_Display_Indent(vMessage)
                end if
                
                ! Point to the base of the tree.
                node => currentTree%baseNode
                
                ! Get the basic component of the node.
                basic => node%basic()
                
                ! Tree walk loop: Walk to each node in the tree and consider whether or not to evolve it.
                treeWalkLoop: do while (associated(node))                   
                   
                   ! Get the basic component of the node.
                   basic => node%basic()    
                   
                   ! Count nodes in the tree.
                   nodesTotalCount=nodesTotalCount+1
                   
                   ! Find the next node that we will process.
                   nodeNext => node%walkTreeWithSatellites()
                   
                   ! Evolve this node if it has a parent (or will transfer to another tree where it will have a parent), exists
                   ! before the output time, has no children (i.e. they've already all been processed), and either exists before
                   ! the final time in its tree, or exists precisely at that time and has some attached event yet to occur.
                   event       =>            thisNode%event
                   hasParent   =  associated(thisNode%parent)
                   treeLimited =  .true.
                   do while (associated(event).and.treeLimited)
                      ! Skip events which occur after the current evolution end time.
                      if (event%time <= endTime) then
                         ! Detect inter-tree events.
                         select type (event)
                         type is (nodeEventSubhaloPromotionInterTree)
                            hasParent  =.true.
                            treeLimited=.false.
                         type is (nodeEventBranchJumpInterTree      )
                            hasParent  =.true.
                            treeLimited=.false.
                         end select
                      end if
                      event => event%next
                   end do
                   evolveCondition: if (                                     &
                        &                     hasParent                                   &
                        &               .and.                                &
                        &                .not.associated(node%firstChild  )  &
                        &               .and.                                &
                        &                (                                   &
                        &                  thisBasicComponent%time() <  endTime           &
                        &                 .or.                               &
                        &                  (                                 &
                        &                    .not.treeLimited                             & ! For nodes that are not tree limited (i.e. have a node which
                        &                   .and.                                         & ! will jump to another tree), allow them to evolve if the node
                        &                    thisBasicComponent%time() == endTime         & ! is at the end time also, since the jump may occur at that time.
                        &                  )                                 &
                        &                )                                   &
                        &               .and.                                &
                        &                (                                   &
                        &                   .not.treeLimited                              &
                        &                .or.                                             &
                        &                .or.                                &
                        &                 (                                  &
                        &                  (                                 &
                        &                     associated(node%event      )   &
                        &                   .or.                             &
                        &                     associated(node%mergeTarget)   &
                        &                  )                                 &
                        &                  .and.                             &
                        &                    basic%time() <= finalTimeInTree &
                        &                 )                                  &
                        &                )                                   &
                        &              ) then

                      ! Flag that a node was evolved.
                      didEvolve=.true.

                      ! Update tree progress counter.
                      nodesEvolvedCount=nodesEvolvedCount+1
                      
                      ! Dump the merger tree structure for later plotting.
                      if (mergerTreesDumpStructure) call Merger_Tree_Dump(currentTree%index,currentTree%baseNode,[node%index()])
                      
                      ! Evolve the node, handling interrupt events. We keep on evolving it until no interrupt is returned (in which case
                      ! the node has reached the requested end time) or the node no longer exists (e.g. if it was destroyed).
                      interrupted=.true.
                      do while (interrupted.and.associated(node))
                         interrupted=.false.                         
                         ! Find maximum allowed end time for this particular node.
                         if (deadlockStatus == deadlockStatusIsReporting) then
                            vMessage="node "
                            write (label,'(e12.6)') basic%time()
                            vMessage=vMessage//node%index()//" (current:target times = "//label
                            write (label,'(e12.6)') endTime
                            vMessage=vMessage//":"//label//")"
                            call Galacticus_Display_Indent(vMessage)
                            endTimeThisNode=Evolve_To_Time(node,endTime,End_Of_Timestep_Task,report=.true.&
                                 &,nodeLock=nodeLock,lockType=lockType)
                            call Galacticus_Display_Unindent("end node")
                            call Deadlock_Add_Node(node,currentTree%index,nodeLock,lockType)
                         else
                            endTimeThisNode=Evolve_To_Time(node,endTime,End_Of_Timestep_Task,report=.false.)
                         end if
                         ! If this node is able to evolve by a finite amount, the tree is not deadlocked.
                         if (endTimeThisNode > basic%time()) deadlockStatus=deadlockStatusIsNotDeadlocked
                         ! Update record of earliest time in the tree.
                         earliestTimeInTree=min(earliestTimeInTree,endTimeThisNode)
                         ! Evolve the node to the next interrupt event, or the end time.
                         call Tree_Node_Evolve(currentTree,node,endTimeThisNode,interrupted,interruptProcedure)
                         ! Check for interrupt.
                         if (interrupted) then
                            ! If an interrupt occured call the specified procedure to handle it.
                            call interruptProcedure(node)
                            ! Something happened so the tree is not deadlocked.
                            deadlockStatus=deadlockStatusIsNotDeadlocked
                         else
                            ! Call routine to handle end of timestep processing.                            
                            if (associated(End_Of_Timestep_Task).and.associated(node)) call End_Of_Timestep_Task(currentTree,node,deadlockStatus)
                         end if
                      end do
                      ! If this halo has reached its parent halo, decide how to handle it.
                      if (associated(thisNode).and.associated(thisNode%parent)) then
                         nodeParent           => node%parent
                         basicParent => nodeParent%basic()
                         if (basic%time() >= basicParent%time()) then
                            ! Parent halo has been reached. Check if the node is the primary (major) progenitor of the parent node.
                            select case (node%isPrimaryProgenitor())
                            case (.false.)
                               ! It is not the major progenitor, so this could be a halo merger event unless the halo is already a
                               ! satellite. Check for satellite status and, if it's not a satellite, process this halo merging
                               ! event. Also record that the tree is not deadlocked, as we are changing the tree state.
                               if (.not.node%isSatellite()) then
                                  deadlockStatus=deadlockStatusIsNotDeadlocked
                                  call Events_Node_Merger(node)
                               end if
                            case (.true.)
                               ! This is the major progenitor, so promote the node to its parent providing that the node has no
                               ! siblings - this ensures that any siblings have already been evolved and become satellites of the
                               ! parent halo. Also record that the tree is not deadlocked, as we are changing the tree state.
                               if (.not.associated(thisNode%sibling).and..not.associated(thisNode%event)) then
                                  deadlockStatus=deadlockStatusIsNotDeadlocked
                                  call Tree_Node_Promote(node)
                               end if
                            end select
                         end if
                      end if
                   else
                      if (deadlockStatus == deadlockStatusIsReporting) then
                         vMessage="node "
                         write (label,'(e12.6)') basic%time()
                         vMessage=vMessage//node%index()//" (current:target times = "//label
                         write (label,'(e12.6)') endTime
                         vMessage=vMessage//":"//label//")"
                         call Galacticus_Display_Indent(vMessage)
                         call Galacticus_Display_Unindent("end node")
                         ! Determine why this node could not be evolved. We check the "has child" condition first as it's the only
                         ! one that provides additional connection between nodes, so leads to the most informative deadlock graph.
                         if      (associated(node%firstChild)) then
                            call Deadlock_Add_Node(node,currentTree%index,node%firstChild,var_str("has child"          ))
                         else if (.not.associated(thisNode%parent)) then
                            call Deadlock_Add_Node(node,currentTree%index,node           ,var_str("no parent"          ))
                         else if (basic%time() >= endTime) then
                            call Deadlock_Add_Node(node,currentTree%index,node           ,var_str("in future of output"))
                         else
                            call Deadlock_Add_Node(node,currentTree%index,node           ,var_str("in future of tree"  ))
                         end if
                      end if
                   end if evolveCondition

                   ! Step to the next node to consider.
                   node => nodeNext
                   
                end do treeWalkLoop
                
                ! Output tree progress information.
                if (treeWalkCount > int(treeWalkCountPreviousOutput*1.1d0)+1) then
                   if (Galacticus_Verbosity_Level() >= verbosityLevel) then
                      write (message,'(a,i9,a )') 'Evolving tree [',treeWalkCount,']'
                      call Galacticus_Display_Indent(message,verbosityLevel)
                      write (message,'(a,i9   )') 'Nodes in tree:         ',nodesTotalCount
                      call Galacticus_Display_Message(message,verbosityLevel)
                      write (message,'(a,i9   )') 'Nodes evolved:         ',nodesEvolvedCount
                      call Galacticus_Display_Message(message,verbosityLevel)
                      write (message,'(a,e10.4)') 'Earliest time in tree: ',earliestTimeInTree
                      call Galacticus_Display_Message(message,verbosityLevel)
                      call Galacticus_Display_Unindent('done',verbosityLevel)
                      treeWalkCountPreviousOutput=treeWalkCount
                   end if
                end if
                
                ! Report on current tree if deadlocked.
                if (deadlockStatus == deadlockStatusIsReporting) call Galacticus_Display_Unindent('end tree')
             end if
             ! Move to the next tree.
             currentTree => currentTree%nextTree
          end do treesLoop

          ! Perform any tree events.
          call Perform_Tree_Events(tree,deadlockStatus)
          
          ! Check deadlocking.
          if (didEvolve .and. deadlockStatus /= deadlockStatusIsNotDeadlocked) then
             if (deadlockStatus == deadlockStatusIsReporting) then
                call Galacticus_Display_Unindent("report done")
                call Deadlock_Tree_Output(endTime)
                if (.not.deadlockReporting) then
                   call Galacticus_Error_Report('Merger_Tree_Evolve_To','merger tree appears to be deadlocked (see preceding report) - check timestep criteria')
                else
                   return
                end if
             else
                ! Tree appears to be deadlocked. Check if it is suspendable.
                if (deadlockStatus == deadlockStatusIsSuspendable) then
                   ! Tree is suspendable, so do not attempt to process further, but simply return and flag it for suspension.
                   suspendTree  =.true.                   
                   return
                else
                   ! Tree is truly deadlocked. Switch to reporting mode and do one more pass through the tree.
                   deadlockStatus=deadlockStatusIsReporting
                end if
             end if
          else
             ! No evolution could occur, so the tree is not deadlocked.
             deadlockStatus=deadlockStatusIsNotDeadlocked
          end if
       end do deadlock
       ! Record tree evolution status.
       if (didEvolve) treeDidEvolve=.true.
    end do outerLoop
    return
  end subroutine Merger_Tree_Evolve_To

  recursive function Evolve_To_Time(thisNode,endTime,End_Of_Timestep_Task,report,lockNode,lockType) result(evolveToTime)
    !% Determine the time to which {\normalfont \ttfamily node} should be evolved.
    use Merger_Trees_Evolve_Timesteps_Template
    use Merger_Trees_Evolve_Node
    use Merger_Tree_Timesteps
    use Cosmology_Functions
    use Input_Parameters
    use Galacticus_Error
    use Galacticus_Display
    use String_Handling
    use Evolve_To_Time_Reports
    implicit none
    double precision                                                                  :: evolveToTime
    double precision                               , intent(in   )                    :: endTime
    type            (treeNode                     )                         , pointer :: nodeSatellite                , nodeSibling
    procedure       (End_Of_Timestep_Task_Template), intent(  out)          , pointer :: End_Of_Timestep_Task
    logical                                        , intent(in   )                    :: report
    type            (treeNode                     ), intent(  out), optional, pointer :: nodeLock
    type            (varying_string               ), intent(  out), optional          :: lockType
    procedure       (End_Of_Timestep_Task_Template)                         , pointer :: End_Of_Timestep_Task_Internal
    class           (nodeComponentBasic           )                         , pointer :: basicParent                  , basicSatellite    , &
         &                                                                               basicSibling                 , basic
    class           (nodeComponentSatellite       )                         , pointer :: satelliteSatellite
    class           (nodeEvent                    )                         , pointer :: thisEvent
    class           (treeEvent                    )                         , pointer :: thisTreeEvent
    class           (cosmologyFunctionsClass      )                         , pointer :: cosmologyFunctionsDefault
    double precision                                                                  :: expansionFactor              , expansionTimescale, &
         &                                                                               hostTimeLimit                , time              , &
         &                                                                               timeEarliest
    character       (len=9                        )                                   :: timeFormatted
    type            (varying_string               )                                   :: message
    
    ! Initialize if not yet done.
    !$omp critical (evolveToTimeInitialize)
    if (.not.evolveToTimeInitialized) then
       !@ <inputParameter>
       !@   <name>timestepHostRelative</name>
       !@   <defaultValue>0.1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The maximum allowed relative timestep for node evolution relative to the time of the host halo.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('timestepHostRelative',timestepHostRelative,defaultValue=0.1d0)
       !@ <inputParameter>
       !@   <name>timestepHostAbsolute</name>
       !@   <defaultValue>1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The maximum allowed absolute timestep (in Gyr) for node evolution relative to the time of the host halo.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('timestepHostAbsolute',timestepHostAbsolute,defaultValue=1.0d0)
       evolveToTimeInitialized=.true.
    end if
    !$omp end critical (evolveToTimeInitialize)

    ! Initially set to the global end time.
    evolveToTime=endTime
    if (report) call Evolve_To_Time_Report("start (target): ",evolveToTime)

    ! Initialize the lock node if present.
    if (present(nodeLock)) nodeLock => null()
    if (present(lockType)) lockType = 'null'

    ! Get the basic component of the node.
    basic => node%basic()

    ! Limit time based on satellite status.
    select case (node%isSatellite())
    case (.false.)
       ! Limit to the time of its parent node if this node is not a satellite.
       if (associated(node%parent)) then
          basicParent => node%parent%basic()
          if (parentBasicComponent%time() < evolveToTime) then
             if (present(nodeLock)) nodeLock => node%parent
             if (present(lockType)) lockType =  "promotion"
             evolveToTime=parentBasicComponent%time()
          end if
       end if
       if (report) call Evolve_To_Time_Report("promotion limit: ",evolveToTime)
    case (.true.)
       ! Do not let satellite evolve too far beyond parent.
       if (associated(node%parent%parent)) then
          ! The host halo has a parent, so use the host halo time to limit satellite evolution.
          basicParent => node%parent%basic()
          time=basicParent%time()
       else
          ! The host halo has no parent. The satellite must therefore be evolving to some event (e.g. a merger). We have to allow
          ! it to evolve ahead of the host halo in this case to avoid deadlocks.
          basicParent => node%parent%basic()
          time=max(basicParent%time(),basic%time())
       end if
       ! Check if the host has a child.
       select case (associated(node%parent%firstChild))
       case (.true. )
          ! Host still has a child - do not let the satellite evolve beyond the host.
          hostTimeLimit=max(time,basic%time())
          ! Check for any merge targets directed at this node.
          nodeSatellite => node%firstMergee
          timeEarliest=huge(1.0d0)
          do while (associated(nodeSatellite))
             satelliteSatellite => nodeSatellite%satellite()
             if (nodeSatellite%isSatellite().and.nodeSatellite%parent%isProgenitorOf(node%parent)) &
                  & timeEarliest=min(timeEarliest,satelliteSatellite%timeOfMerging())
             nodeSatellite => nodeSatellite%siblingMergee
          end do
          if (timeEarliest < huge(1.0d0)) hostTimeLimit=max(hostTimeLimit,timeEarliest)
       case (.false.)
          ! Find current expansion timescale. 
          if (timestepHostRelative > 0.0d0) then
             cosmologyFunctionsDefault => cosmologyFunctions()
             expansionFactor           =  cosmologyFunctionsDefault%expansionFactor(time)
             expansionTimescale        =  1.0d0/cosmologyFunctionsDefault%expansionRate(expansionFactor)
             hostTimeLimit             =  max(time+min(timestepHostRelative*expansionTimescale,timestepHostAbsolute),basic%time())
          else
             ! Avoid use of expansion timescale if host absolute timestep is non-positive. This allows static universe cases to be handled.
             hostTimeLimit=max(time+timestepHostAbsolute,basic%time())
          end if
       end select
       ! Limit to this time.
       if (hostTimeLimit < evolveToTime) then
          if (present(nodeLock)) nodeLock => node%parent
          if (present(lockType)) lockType =  "satellite in host"
          evolveToTime=hostTimeLimit
       end if
       if (report) call Evolve_To_Time_Report("satellite in host limit: ",evolveToTime,thisNode%parent%index())
    end select

    ! Also ensure that this node is not evolved beyond the time of any of its current satellites.
    nodeSatellite => node%firstSatellite
    time=basic%time()
    do while (associated(nodeSatellite))
       basicSatellite => nodeSatellite%basic()
       if (max(satelliteBasicComponent%time(),time) < evolveToTime) then
          if (present(nodeLock)) nodeLock => nodeSatellite
          if (present(lockType)) lockType =  "hosted satellite"
          evolveToTime=max(satelliteBasicComponent%time(),time)
       end if
       if (report) call Evolve_To_Time_Report("hosted satellite: ",evolveToTime,satelliteNode%index())
       nodeSatellite => nodeSatellite%sibling
    end do

    ! Also ensure that this node is not evolved beyond the time at which any of its mergees merge. In some cases, the node may
    ! already be in the future of a mergee. In such cases, simply freeze it at the current time.
    nodeSatellite => node%firstMergee
    time=basic%time()
    do while (associated(nodeSatellite))
       satelliteSatellite => nodeSatellite%satellite()
       if (max(satelliteSatelliteComponent%timeOfMerging(),time) < evolveToTime) then
          if (present(nodeLock)) nodeLock => nodeSatellite
          if (present(lockType)) then
             write (timeFormatted,'(f7.4)') max(satelliteSatellite%timeOfMerging(),time)
             lockType =  "mergee ("//trim(timeFormatted)//")"
          end if
          evolveToTime=max(satelliteSatelliteComponent%timeOfMerging(),time)
       end if
       if (report) call Evolve_To_Time_Report("mergee limit: ",evolveToTime,satelliteNode%index())
       nodeSatellite => nodeSatellite%siblingMergee
    end do

    ! Also ensure that a primary progenitor does not evolve in advance of siblings. This is important since we can not promote a
    ! primary progenitor into its parent until all siblings have become satellites in that parent.
    if (node%isPrimaryProgenitor()) then
       nodeSibling => node%sibling
       do while (associated(nodeSibling))
          basicSibling => nodeSibling%basic()
          if (max(thisBasicComponent%time(),siblingBasicComponent%time()) < evolveToTime) then
             if (present(nodeLock)) nodeLock => nodeSibling
             if (present(lockType)) lockType =  "sibling"
             evolveToTime=max(thisBasicComponent%time(),siblingBasicComponent%time())
          end if
          if (report) call Evolve_To_Time_Report("sibling: ",evolveToTime,siblingNode%index())
          nodeSibling => nodeSibling%sibling
       end do
    end if
    
    ! Also ensure that the timestep taken does not exceed the allowed timestep for this specific node.
    if (report) call Galacticus_Display_Indent("timestepping criteria")
    End_Of_Timestep_Task_Internal => null()
    evolveToTime=min(evolveToTime,thisBasicComponent%time()+Time_Step_Get(thisNode,evolveToTime,End_Of_Timestep_Task_Internal,report,lockNode,lockType))
    End_Of_Timestep_Task => End_Of_Timestep_Task_Internal    
    if (report) call Galacticus_Display_Unindent("done")

    ! Also ensure that the timestep doesn't exceed any event attached to the node.
    thisEvent => node%event
    do while (associated(thisEvent))      
       if (max(thisEvent%time,time) <= evolveToTime) then
          if (present(nodeLock)) nodeLock => thisEvent%node
          if (present(lockType)) then
             lockType =  "event ("
             select type (thisEvent)
             type is (nodeEventSubhaloPromotionInterTree)
                lockType=lockType//thisEvent%splitForestUniqueID//":"//thisEvent%pairedNodeID
             type is (nodeEventBranchJumpInterTree      )
                lockType=lockType//thisEvent%splitForestUniqueID//":"//thisEvent%pairedNodeID
             class default
                lockType=lockType//thisEvent%ID
             end select
             lockType=lockType//")"
          end if
          evolveToTime=max(thisEvent%time,time)
          End_Of_Timestep_Task => Perform_Node_Events
       end if
       if (report) then
          message="event ("
          message=message//thisEvent%ID//"): "
          call Evolve_To_Time_Report(char(message),evolveToTime,thisEvent%node%index())
       end if
       thisEvent => thisEvent%next
    end do

    ! Also ensure that the timestep doesn't exceed any event attached to the tree.
    treeEvent_ => node%hostTree%event
    do while (associated(treeEvent_))
       if (max(thisTreeEvent%time,time) <= evolveToTime) then
          if (present(nodeLock)) nodeLock => node%hostTree%baseNode
          if (present(lockType)) then
             lockType =  "tree event ("
             lockType=lockType//treeEvent_%ID//")"
          end if
          evolveToTime=max(thisTreeEvent%time,time)
       end if
       if (report) then
          message="tree event ("
          message=message//treeEvent_%ID//"): "
          call Evolve_To_Time_Report(char(message),evolveToTime,thisNode%hostTree%baseNode%index())
       end if
       treeEvent_ => treeEvent_%next
    end do

    ! Check that end time exceeds current time.
    if (evolveToTime < thisBasicComponent%time()) then
       message='end time ('
       write (timeFormatted,'(f7.4)') evolveToTime
       message=message//trim(timeFormatted)//' Gyr) is before current time ('
       write (timeFormatted,'(f7.4)') basic%time()
       message=message//trim(timeFormatted)//' Gyr) of node '
       message=message//node%index()
       message=message//' (time difference is '
       write (timeFormatted,'(e8.2)') thisBasicComponent%time()-evolveToTime
       message=message//trim(timeFormatted)
       if (.not.Tree_Node_Is_Accurate(thisBasicComponent%time(),evolveToTime)) then
          ! End time is well before current time. This is an error. Call ourself with reporting switched on to generate a report
          ! on the time limits.
          message=message//' Gyr)'
          if (.not.report) time=Evolve_To_Time(thisNode,endTime,End_Of_Timestep_Task,report=.true.)
          call Galacticus_Error_Report('Evolve_To_Time',message)
       else
          ! End time is before current time, but only by a small amount, simply reset the current time to the end time.
          message=message//' Gyr) - this should happen infrequently'
          call Galacticus_Display_Message(message,verbosityInfo)
          call thisBasicComponent%timeSet(evolveToTime)
       end if
    end if
    return
  end function Evolve_To_Time

  subroutine Deadlock_Add_Node(node,treeIndex,nodeLock,lockType)
    !% Add a node to the deadlocked nodes list.
    implicit none
    type   (treeNode      ), intent(in   ), pointer :: nodeLock        , node
    integer(kind=kind_int8), intent(in   )          :: treeIndex
    type   (varying_string), intent(in   )          :: lockType
    type   (deadlockList  )               , pointer :: deadlockThisNode

    ! Add a node to the deadlock linked list.
    if (associated(deadlockHeadNode)) then
       deadlockThisNode => deadlockHeadNode
       do while (associated(deadlockThisNode%next))
          deadlockThisNode => deadlockThisNode%next
       end do
       allocate(deadlockThisNode%next)
       deadlockThisNode => deadlockThisNode%next
    else
       allocate(deadlockHeadNode)
       deadlockThisNode => deadlockHeadNode
    end if
    ! Set properties.
    deadlockThisNode%node      => node
    deadlockThisNode%treeIndex =  treeIndex
    deadlockThisNode%nodeLock  => nodeLock
    deadlockThisNode%lockType  =  lockType
    return
  end subroutine Deadlock_Add_Node

  subroutine Deadlock_Tree_Output(endTime)
    !% Output the deadlocked nodes in {\normalfont \ttfamily dot} format.
    use String_Handling
    implicit none
    double precision                    , intent(in   ) :: endTime
    type            (deadlockList      ), pointer       :: nodeLock          , testNode, node
    class           (nodeComponentBasic), pointer       :: basic
    type            (treeNode          ), pointer       :: nodeParent
    logical                                             :: foundLockNode
    integer                                             :: treeUnit
    integer         (kind=kind_int8    )                :: uniqueID
    logical                                             :: inCycle           , nodesAdded
    character       (len=20            )                :: color             , style
    type            (varying_string    )                :: deadlockFileName

    ! If no deadlock list exists, simply return.
    if (.not.associated(deadlockHeadNode)) return
    ! Begin tree.
    deadlockFileName=var_str('galacticusDeadlockTree_')//deadlockHeadNode%node%hostTree%baseNode%uniqueID()//'.gv'
    open(newUnit=treeUnit,file=char(deadlockFileName),status='unknown',form='formatted')
    write (treeUnit,*) 'digraph Tree {'
    ! Find any nodes that cause a lock but which are not in our list.
    nodesAdded=.true.
    do while (nodesAdded)
       nodesAdded=.false.
       node => deadlockHeadNode
       do while (associated(node))
          if (associated(node%nodeLock)) then
             testNode => deadlockHeadNode
             foundLockNode=.false.
             do while (associated(testNode).and..not.foundLockNode)
                foundLockNode=(associated(node%nodeLock,testNode%node))
                testNode => testNode%next
             end do
             if (.not.foundLockNode) then
                nodesAdded =  .true.
                testNode   => deadlockHeadNode
                do while (associated(testNode%next))
                   testNode => testNode%next
                end do
                allocate(testNode%next)
                testNode => testNode%next
                ! Find root node.
                nodeParent => node%nodeLock
                do while (associated(nodeParent%parent))
                   nodeParent => nodeParent%parent
                end do
                ! Set properties.
                testNode%node      => node%nodeLock
                testNode%treeIndex =  nodeParent%index()
                testNode%nodeLock  => null()
                testNode%lockType  =  "unknown"
                basic => node%nodeLock%basic()
                if (associated(node%nodeLock%firstChild)) then
                   testNode%lockType = "child"
                   nodeLock        => deadlockHeadNode
                   do while (associated(nodeLock))
                      if (associated(node%nodeLock%firstChild,nodeLock%node)) then
                         testNode%nodeLock => node%nodeLock%firstChild
                         exit
                      end if
                      nodeLock => nodeLock%next
                   end do
                end if
                if (basic%time() >= endTime) testNode%lockType = "end time"
             end if
          end if
          node => node%next
       end do
    end do

    ! Iterate over all nodes visited.
    node => deadlockHeadNode
    do while (associated(node))
       ! Detect cycles.
       inCycle=.false.
       uniqueID=node%node%uniqueID()
       testNode => node%next
       do while (associated(testNode))
          if (testNode%node%uniqueID() == uniqueID) then
             inCycle=.true.
             exit
          end if
          testNode => testNode%next
       end do
       ! Output node.
       basic => node%node%basic()
       if (inCycle) then
          color="green"
          style="filled"
       else
          color="black"
          style="solid"
       end if
       write (treeUnit,'(a,i16.16,a,a,a,a,a,i16.16,a,i16.16,a,i16.16,a,f13.10,a,a,a)') '"',node%node%uniqueID(),'" [shape=circle, color=',trim(color),', style=',trim(style),' label="',node%node%uniqueid(),':',node%node%index(),'\ntree: ',node%treeIndex,'\ntime: ',basic%time(),'\n',char(node%lockType),'"];'
       if (associated(thisNode%lockNode)) write (treeUnit,'(a,i16.16,a,i16.16,a)') '"',thisNode%node%uniqueID(),'" -> "',thisNode%lockNode%uniqueID(),'"'
       node => node%next
    end do
    ! Close the tree.
    write (treeUnit,*) '}'
    close(treeUnit)
    ! Clean up the deadlock node list.
    thisNode => deadlockHeadNode
    do while (associated(thisNode))
       testNode => thisNode%next
       nullify(thisNode%node    )
       nullify(thisNode%lockNode)
       deallocate(thisNode)
       thisNode => testNode
    end do
    deadlockHeadNode => null()
    return
  end subroutine Deadlock_Tree_Output

  subroutine Perform_Node_Events(tree,node,deadlockStatus)
    !% Perform any events associated with {\normalfont \ttfamily node}.
    use Merger_Trees_Evolve_Deadlock_Status
    implicit none
    type            (mergerTree        ), intent(in   )          :: tree
    type            (treeNode          ), intent(inout), pointer :: node
    integer                             , intent(inout)          :: deadlockStatus
    class           (nodeEvent         )               , pointer :: lastEvent         , nextEvent, &
         &                                                          thisEvent
    double precision                                             :: nodeTime, timeEventEarliest
    logical                                                      :: taskDone
    !GCC$ attributes unused :: tree
    
    ! Get the current time.
    basic => node          %basic()
    nodeTime           =  basic%time ()

    ! Find the current earliest event.    
    thisEvent => thisNode%event
    timeEventEarliest=huge(1.0d0)
    do while (associated(thisEvent))
       timeEventEarliest=min(timeEventEarliest,thisEvent%time)
       thisEvent => thisEvent%next
    end do
    ! Get the first event.
    thisEvent => node%event
    lastEvent => node%event
    ! Iterate over all events.
    do while (associated(thisEvent))
       ! Process the event if it occurs at the present time.
       if (thisEvent%time <= nodeTime .and. thisEvent%time == timeEventEarliest .and. associated(thisEvent%task)) then
          taskDone=thisEvent%task(node,deadlockStatus)
          ! If the node is no longer associated, simply exit (as any events associated with it must have been processed already).
          if (.not.associated(node)) exit
          ! Move to the next event.
          if (taskDone) then
             ! The task was performed successfully, so remove it and move to the next event.
             if (associated(thisEvent,node%event)) then
                node%event => thisEvent%next
                lastEvent      => node %event
             else
                lastEvent%next => thisEvent%next
             end if
             nextEvent => thisEvent%next
             deallocate(thisEvent)
             thisEvent => nextEvent
          else
             ! The task was not performed, so simply move to the next event.
             lastEvent => thisEvent
             thisEvent => thisEvent%next
          end if
       else
          lastEvent => thisEvent
          thisEvent => thisEvent%next
       end if
    end do
    return
  end subroutine Perform_Node_Events

  subroutine Perform_Tree_Events(tree,deadlockStatus)
    !% Perform any events associated with {\normalfont \ttfamily node}.
    implicit none
    type            (mergerTree        ), intent(inout), target  :: tree
    integer                             , intent(inout)          :: deadlockStatus
    type            (treeEvent         )               , pointer :: lastEvent              , nextEvent, thisEvent
    double precision                                             :: treeTimeEarliest
    logical                                                      :: taskDone

    ! Find the earliest time in the tree.
    treeTimeEarliest=tree%earliestTime()
    ! Get the first event.
    thisEvent => tree%event
    lastEvent => tree%event
    ! Iterate over all events.
    do while (associated(thisEvent))
       ! Process the event if it occurs at the present time.
       if (thisEvent%time <= treeTimeEarliest .and. associated(thisEvent%task)) then
          taskDone=thisEvent%task(tree,deadlockStatus)
          ! Move to the next event.
          if (taskDone) then
             ! The task was performed successfully, so remove it and move to the next event.
             if (associated(thisEvent,tree%event)) then
                tree%event => thisEvent%next
                lastEvent      => tree %event
             else
                lastEvent%next => thisEvent%next
             end if
             nextEvent => thisEvent%next
             if (taskDone) deallocate(thisEvent)
             thisEvent => nextEvent
          else
             ! The task was not performed, so simply move to the next event.
             thisEvent => thisEvent%next
          end if
       else
          lastEvent => thisEvent
          thisEvent => thisEvent%next
       end if
    end do
    return
  end subroutine Perform_Tree_Events

end module Merger_Trees_Evolve
