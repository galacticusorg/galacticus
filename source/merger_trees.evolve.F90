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
  logical          :: mergerTreeEvolveToInitialized=.false.

  ! Flag indicating whether or not to fail for trees which do not exist at the final output time.
  logical          :: allTreesExistAtFinalTime

  ! Flag indicating whether to dump merger tree structure after each evolutionary step.
  logical          :: mergerTreesDumpStructure

  ! Variables which limit extent to which satellites can evolve past their parent.
  logical          :: evolveToTimeInitialized=.false.
  double precision :: timestepHostAbsolute,timestepHostRelative

  ! Structure used to store list of nodes for deadlock reporting.
  type :: deadlockList
     type   (deadlockList  ), pointer :: next => null()
     type   (treeNode      ), pointer :: node,lockNode
     integer(kind=kind_int8)          :: treeIndex
     type   (varying_string)          :: lockType
  end type deadlockList
  type(deadlockList), pointer :: deadlockHeadNode => null()
  !$omp threadprivate(deadlockHeadNode)

contains

  subroutine Merger_Tree_Evolve_To(thisTree,endTime)
    !% Evolves all properties of a merger tree to the specified time.
    use Merger_Trees_Evolve_Timesteps_Template
    use Merger_Trees
    use Merger_Trees_Initialize  
    use Merger_Trees_Dump
    use Events_Interrupts
    use Galacticus_Error
    use Galacticus_Display
    use Input_Parameters
    use ISO_Varying_String
    use String_Handling
    use Kind_Numbers
    use Merger_Trees_Evolve_Deadlock_Status
    !# <include directive="mergerTreeEvolveThreadInitialize" type="moduleUse">
    include 'merger_trees.evolve.threadInitialize.moduleUse.inc'
    !# </include>
    implicit none
    type            (mergerTree                   ), intent(inout), target :: thisTree
    double precision                               , intent(in   )         :: endTime
    type            (treeNode                     ), pointer               :: thisNode,nextNode,parentNode,lockNode
    double precision                               , parameter             :: timeTolerance=1.0d-5
    double precision                               , parameter             :: largeTime    =1.0d10
    procedure       (Interrupt_Procedure_Template ), pointer               :: interruptProcedure
    procedure       (End_Of_Timestep_Task_Template), pointer               :: End_Of_Timestep_Task
    integer                                        , parameter             :: verbosityLevel=3
    class           (nodeComponentBasic           ), pointer               :: thisBasicComponent,parentBasicComponent&
         &,baseNodeBasicComponent
    type            (mergerTree                   ), pointer               :: currentTree
    integer                                                                :: nodesEvolvedCount,nodesTotalCount,treeWalkCount&
         &,treeWalkCountPreviousOutput,deadlockStatus
    double precision                                                       :: endTimeThisNode,earliestTimeInTree,finalTimeInTree
    logical                                                                :: interrupted,didEvolve
    character       (len=12                       )                        :: label
    character       (len=35)                                               :: message
    type            (varying_string               )                        :: vMessage,lockType
    logical                                                                :: anyTreeExistsAtOutputTime

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
          !@     Specifies whether merger tree structure should be dumped to a \href{http://www.graphviz.org/}{\sc dot} file.
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
    !# <include directive="mergerTreeEvolveThreadInitialize" type="functionCall" functionType="void">
    include 'merger_trees.evolve.threadInitialize.inc'
    !# </include>

    ! Iterate through all trees.
    anyTreeExistsAtOutputTime=.false.
    currentTree => thisTree
    do while (associated(currentTree))
       ! Initialize the tree if necessary.
       call Merger_Tree_Initialize(currentTree)
       ! Check that the output time is not after the end time of this tree.
       baseNodeBasicComponent => currentTree%baseNode%basic()
       if (endTime > baseNodeBasicComponent%time()) then
          ! Final time is exceeded. Check if by a significant factor.
          if (endTime > baseNodeBasicComponent%time()*(1.0d0+timeTolerance)) then
             ! Exceeded by a significant factor - report an error. Check if such behavior is expected.
             if (allTreesExistAtFinalTime) then
                ! It is not, write an error and exit.
                vMessage='requested time exceeds the final time in the tree'//char(10)
                vMessage=vMessage//' HELP: If you expect that not all trees will exist at the latest requested'//char(10)
                vMessage=vMessage//'       output time (this can happen when using trees extracted from N-body'//char(10)
                vMessage=vMessage//'       simulations for example) set the following in your input parameter file:'//char(10)//char(10)
                vMessage=vMessage//'         <parameter>'//char(10)
                vMessage=vMessage//'          <name>allTreesExistAtFinalTime</name>'//char(10)
                vMessage=vMessage//'          <value>false</value>'//char(10)
                vMessage=vMessage//'         </parameter>'//char(10)
                call Galacticus_Error_Report('Merger_Tree_Evolve_To',vMessage)
             end if
          else
             ! Not exceeded by a significant factor (can happen due to approximation errors). Simply reset to actual time requested.
             call baseNodeBasicComponent%timeSet(endTime)
             anyTreeExistsAtOutputTime=.true.
          end if
       else
          anyTreeExistsAtOutputTime=.true.
       end if

       ! Move to the next tree.
       currentTree => currentTree%nextTree
    end do
    ! Return if none of these trees exist at the output time.
    if (.not.anyTreeExistsAtOutputTime) return

    ! Outer loop: This causes the tree to be repeatedly walked and evolved until it has been evolved all the way to the specified
    ! end time. We stop when no nodes were evolved, which indicates that no further evolution is possible.
    didEvolve=.true.
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

       ! Set the deadlock status to deadlocked initially.
       deadlockStatus=isDeadlocked
       ! Enter loop for deadlock reporting.
       deadlock : do while (deadlockStatus /= isNotDeadlocked)

          ! Post a deadlocking message.
          if (deadlockStatus == isReporting) call Galacticus_Display_Indent("Deadlock report follows")

          ! Iterate through all trees.
          currentTree => thisTree
          treesLoop: do while (associated(currentTree))

             ! Report on current tree if deadlocked.
             if (deadlockStatus == isReporting) then
                vMessage="tree "
                vMessage=vMessage//currentTree%index
                call Galacticus_Display_Indent(vMessage)
             end if

             ! Point to the base of the tree.
             thisNode => currentTree%baseNode

             ! Get the basic component of the node.
             thisBasicComponent => thisNode%basic()

             ! Record the final time in this tree.
             finalTimeInTree=thisBasicComponent%time()

             ! Tree walk loop: Walk to each node in the tree and consider whether or not to evolve it.
             treeWalkLoop: do while (associated(thisNode))

                ! Get the basic component of the node.
                thisBasicComponent => thisNode%basic()

                ! Count nodes in the tree.
                nodesTotalCount=nodesTotalCount+1

                ! Find the next node that we will process.
                call thisNode%walkTreeWithSatellites(nextNode)

                ! Evolve this node if it has a parent, exists before the output time, has no children
                ! (i.e. they've already all been processed), and either exists before the final time
                ! in its tree, or exists precisely at that time and has some attached event yet to occur.
                evolveCondition: if (                                                 &
                     &                     associated(thisNode%parent      )          &
                     &               .and.                                            &
                     &                .not.associated(thisNode%firstChild  )          &
                     &               .and.                                            & 
                     &                   thisBasicComponent%time() <  endTime         &
                     &               .and.                                            &
                     &                (                                               &
                     &                   thisBasicComponent%time() <  finalTimeInTree &
                     &                .or.                                            &
                     &                 (                                              &
                     &                  (                                             &
                     &                     associated(thisNode%event      )           &
                     &                   .or.                                         &
                     &                     associated(thisNode%mergeTarget)           &
                     &                  )                                             &
                     &                  .and.                                         &
                     &                   thisBasicComponent%time() <= finalTimeInTree &
                     &                 )                                              &
                     &                )                                               &
                     &              ) then

                   ! Flag that a node was evolved.
                   didEvolve=.true.

                   ! Update tree progress counter.
                   nodesEvolvedCount=nodesEvolvedCount+1

                   ! Dump the merger tree structure for later plotting.
                   if (mergerTreesDumpStructure) call Merger_Tree_Dump(currentTree%index,currentTree%baseNode,[thisNode%index()])

                   ! Evolve the node, handling interrupt events. We keep on evolving it until no interrupt is returned (in which case
                   ! the node has reached the requested end time) or the node no longer exists (e.g. if it was destroyed).
                   interrupted=.true.
                   do while (interrupted.and.associated(thisNode))
                      interrupted=.false.

                      ! Find maximum allowed end time for this particular node.
                      if (deadlockStatus == isReporting) then
                         vMessage="node "
                         write (label,'(e12.6)') thisBasicComponent%time()
                         vMessage=vMessage//thisNode%index()//" (current:target times = "//label
                         write (label,'(e12.6)') endTime
                         vMessage=vMessage//":"//label//")"
                         call Galacticus_Display_Indent(vMessage)
                         endTimeThisNode=Evolve_To_Time(thisNode,endTime,End_Of_Timestep_Task,report=.true.&
                              &,lockNode=lockNode,lockType=lockType)
                         call Galacticus_Display_Unindent("end node")
                         call Deadlock_Add_Node(thisNode,currentTree%index,lockNode,lockType)
                      else
                         endTimeThisNode=Evolve_To_Time(thisNode,endTime,End_Of_Timestep_Task,report=.false.)
                      end if
                      ! If this node is able to evolve by a finite amount, the tree is not deadlocked.
                      if (endTimeThisNode > thisBasicComponent%time()) deadlockStatus=isNotDeadlocked

                      ! Update record of earliest time in the tree.
                      earliestTimeInTree=min(earliestTimeInTree,endTimeThisNode)
                      ! Evolve the node to the next interrupt event, or the end time.
                      call currentTree%evolveNode(thisNode,endTimeThisNode,interrupted,interruptProcedure)

                      ! Check for interrupt.
                      if (interrupted) then
                         ! If an interrupt occured call the specified procedure to handle it.
                         call interruptProcedure(thisNode)
                         ! Something happened so the tree is not deadlocked.
                         deadlockStatus=isNotDeadlocked
                      else
                         ! Call routine to handle end of timestep processing.
                         if (associated(End_Of_Timestep_Task)) call End_Of_Timestep_Task(currentTree,thisNode,deadlockStatus)
                      end if
                   end do

                   ! If this halo has reached its parent halo, decide how to handle it.
                   if (associated(thisNode)) then
                      parentNode           => thisNode%parent
                      parentBasicComponent => parentNode%basic()
                      if (thisBasicComponent%time() >= parentBasicComponent%time()) then
                         ! Parent halo has been reached. Check if the node is the primary (major) progenitor of the parent node.
                         select case (thisNode%isPrimaryProgenitor())
                         case (.false.)
                            ! It is not the major progenitor, so this could be a halo merger event unless the halo is already a
                            ! satellite. Check for satellite status and, if it's not a satellite, process this halo merging event.
                            if (.not.thisNode%isSatellite()) call currentTree%mergeNode(thisNode)
                         case (.true.)
                            ! This is the major progenitor, so promote the node to its parent as it is the main progenitor providing
                            ! that the node has no siblings - this ensures that any siblings have already been evolved and become
                            ! satellites of the parent halo.
                            if (.not.associated(thisNode%sibling)) call currentTree%promoteNode(thisNode)
                         end select
                      end if
                   end if
                end if evolveCondition

                ! Step to the next node to consider.
                thisNode => nextNode

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
             if (deadlockStatus == isReporting) call Galacticus_Display_Unindent('end tree')

             ! Move to the next tree.
             currentTree => currentTree%nextTree
          end do treesLoop
          if (didEvolve .and. deadlockStatus /= isNotDeadlocked) then
             if (deadlockStatus == isReporting) then
                call Galacticus_Display_Unindent("report done")
                call Deadlock_Tree_Output(endTime)
                call Galacticus_Error_Report('Merger_Tree_Evolve_To','merger tree appears to be deadlocked (see preceding report) - check timestep criteria')
             else
                ! Tree is deadlocked. Switch to reporting mode and do one more pass through the tree.
                deadlockStatus=isReporting
             end if
          else
             ! No evolution could occur, so the tree is not deadlocked.
             deadlockStatus=isNotDeadlocked
          end if

       end do deadlock

    end do outerLoop

    return
  end subroutine Merger_Tree_Evolve_To

  double precision function Evolve_To_Time(thisNode,endTime,End_Of_Timestep_Task,report,lockNode,lockType)
    !% Determine the time to which {\tt thisNode} should be evolved.
    use Merger_Tree_Timesteps
    use Cosmology_Functions
    use Input_Parameters
    use Galacticus_Error
    use Galacticus_Display
    use ISO_Varying_String
    use String_Handling
    use Merger_Trees
    use Evolve_To_Time_Reports
    use Kind_Numbers
    implicit none
    type(treeNode),            intent(inout), pointer           :: thisNode
    double precision,          intent(in)                       :: endTime
    type(treeNode),                           pointer           :: satelliteNode
    procedure(),               intent(out),   pointer           :: End_Of_Timestep_Task
    logical,                   intent(in)                       :: report
    type(treeNode),            intent(out),   pointer, optional :: lockNode
    type(varying_string),      intent(out),            optional :: lockType  
    procedure(),                              pointer           :: End_Of_Timestep_Task_Internal
    class(nodeComponentBasic),                pointer           :: thisBasicComponent,parentBasicComponent,satelliteBasicComponent,siblingBasicComponent
    class(nodeComponentSatellite),            pointer           :: satelliteSatelliteComponent
    type(nodeEvent),                          pointer           :: thisEvent
    double precision                                            :: time,expansionFactor,expansionTimescale,hostTimeLimit
    character(len=9)                                            :: timeFormatted
    type(varying_string)                                        :: message

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
    Evolve_To_Time=endTime
    if (report) call Evolve_To_Time_Report("start (target): ",Evolve_To_Time)
    
    ! Initialize the lock node if present.
    if (present(lockNode)) lockNode => null()
    if (present(lockType)) lockType ='null'

    ! Get the basic component of the node.
    thisBasicComponent => thisNode%basic()

    ! Limit time based on satellite status.
    select case (thisNode%isSatellite())
    case (.false.)
       ! Limit to the time of its parent node if this node is not a satellite.
       if (associated(thisNode%parent)) then
          parentBasicComponent => thisNode%parent%basic()
          if (parentBasicComponent%time() < Evolve_To_Time) then
             if (present(lockNode)) lockNode => thisNode%parent
             if (present(lockType)) lockType =  "promotion"
             Evolve_To_Time=parentBasicComponent%time()
          end if
       end if
       if (report) call Evolve_To_Time_Report("promotion limit: ",Evolve_To_Time)
    case (.true.)
       ! Do not let satellite evolve too far beyond parent.
       ! Get the parent basic component.
       parentBasicComponent => thisNode%parent%basic()
       ! Get current cosmic time.
       time=parentBasicComponent%time()
       ! Find current expansion timescale.
       expansionFactor=Expansion_Factor(time)
       expansionTimescale=1.0d0/Expansion_Rate(expansionFactor)
       ! Determine suitable timestep.
       hostTimeLimit=time+min(timestepHostRelative*expansionTimescale,timestepHostAbsolute)
       ! Limit to this time.
       if (hostTimeLimit < Evolve_To_Time) then
          if (present(lockNode)) lockNode => thisNode%parent
          if (present(lockType)) lockType =  "satellite in host"
          Evolve_To_Time=hostTimeLimit
       end if
       if (report) call Evolve_To_Time_Report("satellite in host limit: ",Evolve_To_Time)
    end select

    ! Also ensure that this node is not evolved beyond the time of any of its current satellites.
    satelliteNode => thisNode%firstSatellite
    time=thisBasicComponent%time()
    do while (associated(satelliteNode))
       satelliteBasicComponent => satelliteNode%basic()
       if (max(satelliteBasicComponent%time(),time) < Evolve_To_Time) then
          if (present(lockNode)) lockNode => satelliteNode
          if (present(lockType)) lockType =  "hosted satellite"
          Evolve_To_Time=max(satelliteBasicComponent%time(),time)
       end if
       if (report) call Evolve_To_Time_Report("hosted satellite: ",Evolve_To_Time,satelliteNode%index())
       satelliteNode => satelliteNode%sibling
    end do

    ! Also ensure that this node is not evolved beyond the time at which any of its mergees merge. In some cases, the node may
    ! already be in the future of a mergee. In such cases, simply freeze it at the current time.
    satelliteNode => thisNode%firstMergee
    time=thisBasicComponent%time()
    do while (associated(satelliteNode))
       satelliteSatelliteComponent => satelliteNode%satellite()
       if (max(satelliteSatelliteComponent%timeOfMerging(),time) < Evolve_To_Time) then
          if (present(lockNode)) lockNode => satelliteNode
          if (present(lockType)) then
             write (timeFormatted,'(f7.4)') max(satelliteSatelliteComponent%timeOfMerging(),time)
             lockType =  "mergee ("//trim(timeFormatted)//")"
          end if
          Evolve_To_Time=max(satelliteSatelliteComponent%timeOfMerging(),time)
       end if
       if (report) call Evolve_To_Time_Report("mergee limit: ",Evolve_To_Time,satelliteNode%index())
       satelliteNode => satelliteNode%siblingMergee
    end do

    ! Also ensure that a primary progenitor does not evolve in advance of siblings. This is important since we can not promote a
    ! primary progenitor into its parent until all siblings have become satellites in that parent.
    if (thisNode%isPrimaryProgenitor().and.associated(thisNode%sibling)) then
       siblingBasicComponent => thisNode%sibling%basic()
       if (max(thisBasicComponent%time(),siblingBasicComponent%time()) < Evolve_To_Time) then
          if (present(lockNode)) lockNode => thisNode%sibling
          if (present(lockType)) lockType =  "sibling"
          Evolve_To_Time=max(thisBasicComponent%time(),siblingBasicComponent%time())
       end if
       if (report) call Evolve_To_Time_Report("sibling: ",Evolve_To_Time,thisNode%sibling%index())
    end if

    ! Also ensure that the timestep taken does not exceed the allowed timestep for this specific node.
    if (report) call Galacticus_Display_Indent("timestepping criteria")
    End_Of_Timestep_Task_Internal => null()
    Evolve_To_Time=min(Evolve_To_Time,thisBasicComponent%time()+Time_Step_Get(thisNode,Evolve_To_Time,End_Of_Timestep_Task_Internal,report,lockNode,lockType))
    End_Of_Timestep_Task => End_Of_Timestep_Task_Internal
    if (report) call Galacticus_Display_Unindent("done")

    ! Also ensure that the timestep doesn't exceed any event attached to the node
    thisEvent => thisNode%event
    do while (associated(thisEvent))
       if (max(thisEvent%time,time) <= Evolve_To_Time) then
          if (present(lockNode)) lockNode => thisEvent%node
          if (present(lockType)) then
             lockType =  "event ("
             lockType=lockType//thisEvent%ID//")"
          end if
          Evolve_To_Time=max(thisEvent%time,time)
          End_Of_Timestep_Task => Perform_Node_Events
       end if
       if (report) then
          message="event ("
          message=message//thisEvent%ID//"): "
          call Evolve_To_Time_Report(char(message),Evolve_To_Time,thisEvent%node%index())
       end if
       thisEvent => thisEvent%next
    end do

    ! Check that end time exceeds current time.
    if (Evolve_To_Time < thisBasicComponent%time()) then
       message='end time ('
       write (timeFormatted,'(f7.4)') Evolve_To_Time
       message=message//trim(timeFormatted)//' Gyr) is before current time ('
       write (timeFormatted,'(f7.4)') thisBasicComponent%time()
       message=message//trim(timeFormatted)//' Gyr) of node '
       message=message//thisNode%index()
       message=message//' (time difference is '
       write (timeFormatted,'(e8.2)') thisBasicComponent%time()-Evolve_To_Time
       message=message//trim(timeFormatted)
       if (.not.Tree_Node_Is_Accurate(thisBasicComponent%time(),Evolve_To_Time)) then
          ! End time is well before current time. This is an error.
          message=message//' Gyr)'
          call Galacticus_Error_Report('Evolve_To_Time',message)
       else
          ! End time is before current time, but only by a small amount, simply reset the current time to the end time.
          message=message//' Gyr) - this should happen infrequently'
          call Galacticus_Display_Message(message,verbosityInfo)
          call thisBasicComponent%timeSet(Evolve_To_Time)
       end if
    end if
    return
  end function Evolve_To_Time

  subroutine Deadlock_Add_Node(thisNode,treeIndex,lockNode,lockType)
    !% Add a node to the deadlocked nodes list.
    implicit none
    type   (treeNode      ), pointer, intent(in   ) :: thisNode,lockNode
    integer(kind=kind_int8),          intent(in   ) :: treeIndex
    type   (varying_string),          intent(in   ) :: lockType
    type   (deadlockList  ), pointer                :: deadlockThisNode

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
    deadlockThisNode%node      => thisNode
    deadlockThisNode%treeIndex =  treeIndex
    deadlockThisNode%lockNode  => lockNode
    deadlockThisNode%lockType  =  lockType
    return
  end subroutine Deadlock_Add_Node

  subroutine Deadlock_Tree_Output(endTime)
    !% Output the deadlocked nodes in {\tt dot} format.
    implicit none
    double precision                    , intent(in   ) :: endTime
    type            (deadlockList      ), pointer       :: thisNode,testNode,lockNode
    class           (nodeComponentBasic), pointer       :: thisBasicComponent
    type            (treeNode          ), pointer       :: parentNode
    logical                                             :: foundLockNode
    integer                                             :: treeUnit
    integer         (kind=kind_int8    )                :: uniqueID
    logical                                             :: inCycle
    character       (len=20            )                :: color,style

    ! Begin tree.
    open(newUnit=treeUnit,file='galacticusDeadlockTree.gv',status='unknown',form='formatted')
    write (treeUnit,*) 'digraph Tree {'
    
    ! Find any nodes that cause a lock but which are not in our list.
    thisNode => deadlockHeadNode
    do while (associated(thisNode))
       if (associated(thisNode%lockNode)) then
          testNode => deadlockHeadNode
          foundLockNode=.false.
          do while (associated(testNode).and..not.foundLockNode)
             foundLockNode=(associated(thisNode%lockNode,testNode%node))
             testNode => testNode%next
          end do
          if (.not.foundLockNode) then
             testNode => deadlockHeadNode
             do while (associated(testNode%next))
                testNode => testNode%next
             end do
             allocate(testNode%next)
             testNode => testNode%next
             ! Find root node.
             parentNode => thisNode%lockNode
             do while (associated(parentNode%parent))
                parentNode => parentNode%parent
             end do
             ! Set properties.
             testNode%node      => thisNode%lockNode
             testNode%treeIndex =  parentNode%index()
             testNode%lockNode  => null()
             testNode%lockType  =  "unknown"
             thisBasicComponent => thisNode%lockNode%basic()
             if (associated(thisNode%lockNode%firstChild)) then
                testNode%lockType = "child"
                lockNode        => deadlockHeadNode
                do while (associated(lockNode))
                   if (associated(thisNode%lockNode%firstChild,lockNode%node)) then
                      testNode%lockNode => thisNode%lockNode%firstChild
                      exit
                   end if
                   lockNode => lockNode%next
                end do
             end if
             if (thisBasicComponent%time() >= endTime) testNode%lockType = "end time"
          end if
       end if
       thisNode => thisNode%next
    end do
 
    ! Iterate over all nodes visited.
    thisNode => deadlockHeadNode
    do while (associated(thisNode))
       ! Detect cycles.
       inCycle=.false.
       uniqueID=thisNode%node%uniqueID()
       testNode => thisNode%next
       do while (associated(testNode))
          if (testNode%node%uniqueID() == uniqueID) then
             inCycle=.true.
             exit
          end if
          testNode => testNode%next
       end do
       ! Output node.
       thisBasicComponent => thisNode%node%basic()
       if (inCycle) then
          color="green"
          style="filled"
       else
          color="black"
          style="solid"
       end if
       write (treeUnit,'(a,i16.16,a,a,a,a,a,i16.16,a,i16.16,a,f7.4,a,a,a)') '"',thisNode%node%index(),'" [shape=circle, color=',trim(color),', style=',trim(style),' label="',thisNode%node%index(),'\ntree: ',thisNode%treeIndex,'\ntime: ',thisBasicComponent%time(),'\n',char(thisNode%lockType),'"];'
       if (associated(thisNode%lockNode)) write (treeUnit,'(a,i16.16,a,i16.16,a)') '"',thisNode%node%index(),'" -> "',thisNode%lockNode%index(),'"' ;
       thisNode => thisNode%next
    end do
    ! Close the tree.
    write (treeUnit,*) '}'
    close(treeUnit)
    return
  end subroutine Deadlock_Tree_Output

  subroutine Perform_Node_Events(thisTree,thisNode,deadlockStatus)
    !% Perform any events associated with {\tt thisNode}.
    use Merger_Trees
    implicit none
    type (mergerTree        ), intent(in   )          :: thisTree
    type (treeNode          ), intent(inout), pointer :: thisNode
    integer                  , intent(inout)          :: deadlockStatus
    type (nodeEvent         ),                pointer :: thisEvent,lastEvent,nextEvent
    class(nodeComponentBasic),                pointer :: thisBasicComponent
    double precision                                  :: nodeTime
    logical                                           :: taskDone

    ! Get the current time.
    thisBasicComponent => thisNode          %basic()
    nodeTime           =  thisBasicComponent%time ()
    ! Get the first event.
    thisEvent => thisNode%event
    lastEvent => thisNode%event
    ! Iterate over all events.
    do while (associated(thisEvent))
       ! Process the event if it occurs at the present time.
       if (thisEvent%time <= nodeTime .and. associated(thisEvent%task)) then
          taskDone=thisEvent%task(thisNode,deadlockStatus)
          ! If the node is no longer associated, simply exit (as any events associated with it must have been processed already).
          if (.not.associated(thisNode)) exit
          ! Move to the next event.
          if (taskDone) then
             ! The task was performed successfully, so remove it and move to the next event.
             if (associated(thisEvent,thisNode%event)) then
                thisNode%event => thisEvent%next
                lastEvent      => thisNode %event
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
  end subroutine Perform_Node_Events

end module Merger_Trees_Evolve
