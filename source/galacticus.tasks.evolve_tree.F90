!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements the task of evolving merger trees.

module Galacticus_Tasks_Evolve_Tree
  !% Implements the task of evolving merger trees.
  implicit none
  private
  public :: Galacticus_Task_Evolve_Tree

  ! Flag to indicate if output times have been initialized.
  logical          :: treeEvolveInitialized       =.false.

  ! Parameters controlling which trees will be processed.
  integer          :: treeEvolveWorkerCount               , treeEvolveWorkerNumber

  ! Parameters controlling load average.
  logical          :: treeEvolveLimitLoadAverage
  double precision :: treeEvolveLoadAverageMaximum

contains

  !# <galacticusTask>
  !#  <unitName>Galacticus_Task_Evolve_Tree</unitName>
  !#  <after>Galacticus_Task_Start</after>
  !#  <before>Galacticus_Task_End</before>
  !# </galacticusTask>
  logical function Galacticus_Task_Evolve_Tree()
    !% Evolves the complete set of merger trees as specified.
    use ISO_Varying_String
    use String_Handling
    use Merger_Trees_Evolve
    use Galacticus_Output_Merger_Tree
    use Galacticus_Display
    use Galacticus_Nodes
    use Input_Parameters
    use Galacticus_Output_Times
    use Galacticus_Error
    use Memory_Management
    use System_Load
    !$ use omp_lib
    ! Include modules needed for pre- and post-evolution and pre-construction tasks.
    !# <include directive="mergerTreePreEvolveTask" type="moduleUse">
    include 'galacticus.tasks.evolve_tree.preEvolveTask.moduleUse.inc'
    !# </include>
    !# <include directive="mergerTreePostEvolveTask" type="moduleUse">
    include 'galacticus.tasks.evolve_tree.postEvolveTask.moduleUse.inc'
    !# </include>
    !# <include directive="universePreEvolveTask" type="moduleUse" functionType="void">
    include 'galacticus.tasks.evolve_tree.universePreEvolveTask.moduleUse.inc'
    !# </include>
    implicit none
    type            (mergerTree    ), pointer     , save :: thisTree
    logical                                       , save :: finished                        , skipTree               , &
         &                                                  treeIsNew
    integer                                       , save :: iOutput
    double precision                              , save :: evolveToTime                    , treeTimeEarliest       , &
         &                                                  universalEvolveToTime
    type            (varying_string)              , save :: message
    character       (len=20        )              , save :: label
    !$omp threadprivate(thisTree,finished,skipTree,iOutput,evolveToTime,message,label,treeIsNew,treeTimeEarliest,universalEvolveToTime)
    integer                                              :: iTree
    integer                                       , save :: activeTasks                     , totalTasks
    double precision                , dimension(3), save :: loadAverage
    logical                                       , save :: overloaded                      , treeCanEvolve          , &
         &                                                  treeIsFinished                  , evolutionIsEventLimited, &
         &                                                  success
    !$omp threadprivate(activeTasks,totalTasks,loadAverage,overloaded,treeCanEvolve,treeIsFinished,evolutionIsEventLimited,success)
    character       (len=32        )                     :: treeEvolveLoadAverageMaximumText
    type            (universe      )                     :: universeWaiting                 , universeProcessed
    type            (universeEvent ), pointer     , save :: thisEvent
    !$omp threadprivate(thisEvent)

    ! Initialize the task if necessary.
    if (.not.treeEvolveInitialized) then
       !$omp critical (Tasks_Evolve_Tree_Initialize)
       if (.not.treeEvolveInitialized) then

          ! Get parameters controlling which trees will be processed.
          !@ <inputParameter>
          !@   <name>treeEvolveWorkerCount</name>
          !@   <defaultValue>1</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The number of workers that will work on this calculation.
          !@   </description>
          !@   <type>integer</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('treeEvolveWorkerCount',treeEvolveWorkerCount,defaultValue=1)
          !@ <inputParameter>
          !@   <name>treeEvolveWorkerNumber</name>
          !@   <defaultValue>1</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The number of this worker.
          !@   </description>
          !@   <type>integer</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('treeEvolveWorkerNumber',treeEvolveWorkerNumber,defaultValue=1)
          !@ <inputParameter>
          !@   <name>treeEvolveLimitLoadAverage</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not to limit the load average
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('treeEvolveLimitLoadAverage',treeEvolveLimitLoadAverage,defaultValue=.false.)
          !@ <inputParameter>
          !@   <name>treeEvolveLoadAverageMaximum</name>
          !@   <defaultValue>processorCount</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The maximum load average for which new trees will be processed.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('treeEvolveLoadAverageMaximum',treeEvolveLoadAverageMaximumText,defaultValue="processorCount")
          if (treeEvolveLoadAverageMaximumText == "processorCount" ) then
             treeEvolveLoadAverageMaximum=dble(System_Processor_Count())
          else
             read (treeEvolveLoadAverageMaximumText,*) treeEvolveLoadAverageMaximum
          end if
          ! Flag that this task is now initialized.
          treeEvolveInitialized=.true.
       end if
       !$omp end critical (Tasks_Evolve_Tree_Initialize)
    end if

    ! Ensure the nodes objects are initialized.
    call Galacticus_Nodes_Initialize()

    ! The following processes merger trees, one at a time, to each successive output time, then dumps their contents to file. It
    ! allows for the possibility of "universal events" - events which require all merger trees to reach the same cosmic time. If
    ! such an event exists, each tree is processed up to that time and then pushed onto a stack where it waits to be
    ! processed. Once all trees reach the event time the stack of trees is passed to the event task. Tree processing then
    ! continues by popping trees off of the stack and processing them further (possibly to the next universal event).

    ! Allow events to be attached to the universe.
    !# <include directive="universePreEvolveTask" type="functionCall" functionType="void">
    !#  <functionArgs>universeWaiting</functionArgs>
    include 'galacticus.tasks.evolve_tree.universePreEvolveTask.inc'
    !# </include>

    ! Initialize tree counter and record that we are not finished processing trees.
    finished=.false.
    iTree   =0

    ! Initialize universes which will act as tree stacks. We use two stacks: one for trees waiting to be processed, one for trees
    ! that have already been processed.
    universeWaiting%trees   => null()
    universeProcessed%trees => null()

    ! Begin parallel processing of trees until all work is done.
    !$omp parallel copyin(finished)
    do while (.not.finished)

       ! Attempt to get a new tree to process. We first tree to get a new tree. If no new trees exist, we will look for a tree on
       ! the stack waiting to be processed.
       if (treeEvolveWorkerCount == 1) then
          call Get_Tree(iTree,skipTree,thisTree,finished)
       else
          !$omp critical(Tree_Sharing)
          call Get_Tree(iTree,skipTree,thisTree,finished)
          !$omp end critical(Tree_Sharing)
       end if
       treeIsNew=.not.finished
       ! If no new tree was available, attempt to pop one off the universe stack.
       if (finished) then
          !$omp critical(universeTransform)
          thisTree  => universeWaiting%popTree()
          skipTree  =  .false.
          treeIsNew =  .false.
          finished  =  .not.associated(thisTree)
          !$omp end critical(universeTransform)
       end if

       ! If we got a tree (i.e. we are not "finished") process it.
       if (.not.finished) then
          treeIsFinished=.false.

          ! Skip this tree if necessary.
          if (.not.skipTree) then

             ! Spin while the system is overloaded.
             overloaded=treeEvolveLimitLoadAverage
             do while (overloaded)
                ! Get the load average.
                call System_Load_Get(loadAverage,activeTasks,totalTasks)
                ! If load is above allowed tolerances, sleep for a while.
                overloaded=(loadAverage(1) > treeEvolveLoadAverageMaximum)
                if (overloaded)                         &
                     & call Sleep( 5                    &
                     !$ &         +omp_get_thread_num() &
                     &           )
             end do

             ! If this is a new tree, perform any pre-evolution tasks on it.
             if (treeIsNew) then
                !# <include directive="mergerTreePreEvolveTask" type="functionCall" functionType="void">
                !#  <functionArgs>thisTree</functionArgs>
                include 'galacticus.tasks.evolve_tree.preEvolveTask.inc'
                !# </include>
                message="Evolving tree number "
             else
                message="Resuming tree number "
             end if
             ! Display a message.
             message=message//thisTree%index
             call Galacticus_Display_Indent(message)
             
             ! Iterate evolving the tree until we can evolve no more.
             treeCanEvolve =.true.
             treeIsFinished=.false.
             treeEvolveLoop : do while (treeCanEvolve)
                ! We want to find the maximum time to which we can evolve this tree. This will be the minimum of the next output
                ! time (at which we must stop and output the tree) and the next universal event time (at which we must stop and
                ! perform the event task).  Find the earliest time in the tree.
                treeTimeEarliest=thisTree%earliestTime()
                ! Find the next output time.
                evolveToTime=Galacticus_Next_Output_Time (treeTimeEarliest)
                ! If the tree is at or beyond the final output time, we are done.
                if (evolveToTime < 0.0d0) then
                   treeCanEvolve =.false.
                   treeIsFinished=.true.
                   cycle
                end if
                ! Find the index of this output.
                iOutput=Galacticus_Output_Time_Index(evolveToTime)
                ! Find the earliest universe event.
                !$omp critical(universeTransform)
                thisEvent               => universeWaiting%event
                evolutionIsEventLimited =  .false.
                do while (associated(thisEvent))
                   if (thisEvent%time < evolveToTime) then
                      evolveToTime           =thisEvent%time
                      evolutionIsEventLimited=.true.
                      universalEvolveToTime  =evolveToTime
                   end if
                   thisEvent => thisEvent%next
                end do
                !$omp end critical(universeTransform)
                ! Evolve the tree to the computed time.
                call Merger_Tree_Evolve_To(thisTree,evolveToTime)
                ! Determine what limited evolution.
                if (evolutionIsEventLimited) then
                   ! Tree evolution was limited by a universal event. Therefore it can evolve no further until that event's task
                   ! is performed.
                   treeCanEvolve=.false.
                else
                   ! Tree reached an output time, so output it. We can then continue evolving.
                   write (label,'(f7.2)') evolveToTime
                   message="Output tree data at t="//trim(label)//" Gyr"
                   call Galacticus_Display_Message(message)
                   call Galacticus_Merger_Tree_Output(thisTree,iOutput,evolveToTime,.false.)
                end if
             end do treeEvolveLoop
             ! If tree could not evolve further, but is not finished, push it to the universe stack.
             if (.not.treeIsFinished) then
                !$omp critical(universeTransform)
                call universeProcessed%pushTree(thisTree)
                !$omp end critical(universeTransform)
                thisTree => null()
                ! Unindent messages.
                call Galacticus_Display_Unindent('Suspending tree')
             else
                ! Unindent messages.
                call Galacticus_Display_Unindent('Finished tree'  )
             end if
             
          end if

          ! Destroy the tree.
          if (associated(thisTree)) then
             call thisTree%destroy()
             ! Deallocate the tree.
             call Memory_Usage_Record(sizeof(thisTree),addRemove=-1,memoryType=memoryTypeNodes)
             deallocate(thisTree)
          end if

          ! Perform any post-evolution tasks on the tree.
          if (treeIsFinished) then
             !# <include directive="mergerTreePostEvolveTask" type="functionCall" functionType="void">
             include 'galacticus.tasks.evolve_tree.postEvolveTask.inc'
             !# </include>
          end if

       end if

       ! If any trees were pushed onto the processed stack, then there must be an event to process.
       if (finished) then
          !$omp barrier
          !$omp single
          !$omp critical(universeTransform)
          if (associated(universeProcessed%trees)) then
             ! Transfer processed trees back to the waiting universe.
             universeWaiting  %trees => universeProcessed%trees
             universeProcessed%trees => null()
             ! Find the event to process.
             thisEvent => universeWaiting%event
             do while (associated(thisEvent))
                if (thisEvent%time < universalEvolveToTime) then
                   call Galacticus_Error_Report('Galacticus_Task_Evolve_Tree','a universal event exists in the past - this should not happen')
                else if (thisEvent%time == universalEvolveToTime) then
                   success=thisEvent%task(universeWaiting)
                   if (success) call universeWaiting%removeEvent(thisEvent)
                   exit
                end if
                thisEvent => thisEvent%next
             end do
             ! Mark that there is more work to do.
             finished=.false.
          end if
          !$omp end critical(universeTransform)
          !$omp end single copyprivate(finished)
       end if

    end do
    !$omp end parallel

    Galacticus_Task_Evolve_Tree=.false.
    return
  end function Galacticus_Task_Evolve_Tree

  subroutine Get_Tree(iTree,skipTree,thisTree,finished)
    !% Get a tree to process.
    use Galacticus_Nodes
    use Merger_Tree_Construction
    !# <include directive="mergerTreePreTreeConstructionTask" type="moduleUse">
    include 'galacticus.tasks.evolve_tree.preConstructionTask.moduleUse.inc'
    !# </include>
    implicit none
    integer            , intent(inout)          :: iTree
    logical            , intent(  out)          :: skipTree
    logical            , intent(inout)          :: finished
    type   (mergerTree), intent(  out), pointer :: thisTree

    ! Increment the tree counter.
    iTree=iTree+1
    ! Decide whether or not to skip this tree.
    skipTree=.not.(modulo(iTree-1+(iTree-1)/treeEvolveWorkerCount,treeEvolveWorkerCount) == treeEvolveWorkerNumber-1)
    ! Perform any pre-tree construction tasks.
    !# <include directive="mergerTreePreTreeConstructionTask" type="functionCall" functionType="void">
    include 'galacticus.tasks.evolve_tree.preConstructionTask.inc'
    !# </include>

    ! Get a tree.
    thisTree => Merger_Tree_Create(skipTree)
    finished=finished.or..not.associated(thisTree)
    return
  end subroutine Get_Tree

end module Galacticus_Tasks_Evolve_Tree
