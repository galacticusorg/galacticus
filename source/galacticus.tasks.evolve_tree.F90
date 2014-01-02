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

  ! Parameters controlling load averaging and thread locking.
  logical          :: treeEvolveLimitLoadAverage          , treeEvolveThreadLock
  double precision :: treeEvolveLoadAverageMaximum
  integer          :: treeEvolveThreadsMaximum

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
    use Memory_Management
    use System_Load
    use Semaphores
    !$ use omp_lib
    ! Include modules needed for pre- and post-evolution and pre-construction tasks.
    !# <include directive="mergerTreePreEvolveTask" type="moduleUse">
    include 'galacticus.tasks.evolve_tree.preEvolveTask.moduleUse.inc'
    !# </include>
    !# <include directive="mergerTreePostEvolveTask" type="moduleUse">
    include 'galacticus.tasks.evolve_tree.postEvolveTask.moduleUse.inc'
    !# </include>
    implicit none
    type            (mergerTree    ), pointer     , save :: thisTree
    logical                                       , save :: finished                        , skipTree
    integer                                       , save :: iOutput
    double precision                              , save :: outputTime
    type            (varying_string)              , save :: message
    character       (len=20        )              , save :: label
    !$omp threadprivate(thisTree,finished,skipTree,iOutput,outputTime,message,label)
    integer                                              :: iTree
    integer                                       , save :: activeTasks                     , totalTasks
    double precision                , dimension(3), save :: loadAverage
    logical                                       , save :: overloaded
    !$omp threadprivate(activeTasks,totalTasks,loadAverage,overloaded)
    type            (semaphore     ), pointer            :: galacticusMutex
    character       (len=32        )                     :: treeEvolveLoadAverageMaximumText,treeEvolveThreadsMaximumText

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
          !@ <inputParameter>
          !@   <name>treeEvolveThreadLock</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not to limit the number of threads across all \glc\ processes.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('treeEvolveThreadLock',treeEvolveThreadLock,defaultValue=.true.)
          !@ <inputParameter>
          !@   <name>treeEvolveThreadsMaximum</name>
          !@   <defaultValue>processorCount</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The maximum number of active threads across all \glc\ processes.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('treeEvolveThreadsMaximum',treeEvolveThreadsMaximumText,defaultValue="processorCount")
          if (treeEvolveThreadsMaximumText == "processorCount") then
             treeEvolveThreadsMaximum=System_Processor_Count()
          else
             read (treeEvolveThreadsMaximumText,*) treeEvolveThreadsMaximum
          end if
          ! Flag that this task is now initialized.
          treeEvolveInitialized=.true.
       end if
       !$omp end critical (Tasks_Evolve_Tree_Initialize)
    end if

    ! Ensure the nodes objects are initialized.
    call Galacticus_Nodes_Initialize()

    ! Begin looping through available trees.
    finished=.false.
    iTree=0

    ! Create a semaphore if threads are being locked.
    if (treeEvolveThreadLock) galacticusMutex => Semaphore_Open("/galacticus",treeEvolveThreadsMaximum)
    
    !$omp parallel copyin(finished)
    do while (.not.finished)

       ! If locking threads, claim one.
       if (treeEvolveThreadLock) call galacticusMutex%wait()

       ! Increment the tree number.
       if (treeEvolveWorkerCount == 1) then
          call Get_Tree(iTree,skipTree,thisTree,finished)
       else
          !$omp critical(Tree_Sharing)
          call Get_Tree(iTree,skipTree,thisTree,finished)
          !$omp end critical(Tree_Sharing)
       end if

       if (.not.finished) then

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

             ! Perform any pre-evolution tasks on the tree.
             !# <include directive="mergerTreePreEvolveTask" type="functionCall" functionType="void">
             !#  <functionArgs>thisTree</functionArgs>
             include 'galacticus.tasks.evolve_tree.preEvolveTask.inc'
             !# </include>

             ! Display a message.
             message="Evolving tree number "
             message=message//thisTree%index
             call Galacticus_Display_Indent(message)

             ! Loop over output times.
             outputTimeLoop : do iOutput=1,Galacticus_Output_Time_Count()

                ! Get the output time.
                outputTime=Galacticus_Output_Time(iOutput)

                ! Evolve the tree to the output time.
                call Merger_Tree_Evolve_To(thisTree,outputTime)

                ! Output the merger tree.
                write (label,'(f7.2)') outputTime
                message="Output tree data at t="//trim(label)//" Gyr"
                call Galacticus_Display_Message(message)
                call Galacticus_Merger_Tree_Output(thisTree,iOutput,outputTime,.false.)

             end do outputTimeLoop

             ! Unindent messages.
             call Galacticus_Display_Unindent('Finished tree')

          end if

          ! Destroy the tree.
          if (associated(thisTree)) then
             call thisTree%destroy()
             ! Deallocate the tree.
             call Memory_Usage_Record(sizeof(thisTree),addRemove=-1,memoryType=memoryTypeNodes)
             deallocate(thisTree)
          end if

          ! Perform any post-evolution tasks on the tree.
          !# <include directive="mergerTreePostEvolveTask" type="functionCall" functionType="void">
          include 'galacticus.tasks.evolve_tree.postEvolveTask.inc'
          !# </include>

       end if

       ! If locking threads, release ours.
       if (treeEvolveThreadLock) call galacticusMutex%post()

    end do
    !$omp end parallel

    ! Close the semaphore.
    if (treeEvolveThreadLock) call galacticusMutex%close()

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
