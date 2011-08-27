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


!% Contains a module which implements evolution of merger trees.

module Merger_Trees_Evolve
  !% Implements evolution of merger trees. 
  use Tree_Nodes
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
    !# <include directive="mergerTreeEvolveThreadInitialize" type="moduleUse">
    include 'merger_trees.evolve.threadInitialize.moduleUse.inc'
    !# </include>
    implicit none
    type(mergerTree),                         intent(inout) :: thisTree
    double precision,                         intent(in)    :: endTime
    type(treeNode),                           pointer       :: thisNode,nextNode
    double precision,                         parameter     :: timeTolerance=1.0d-5
    double precision,                         parameter     :: largeTime    =1.0d10
    procedure(Interrupt_Procedure_Template),  pointer       :: interruptProcedure
    procedure(End_Of_Timestep_Task_Template), pointer       :: End_Of_Timestep_Task
    integer,                                  parameter     :: verbosityLevel=3
    integer                                                 :: nodesEvolvedCount,treeWalkCount,treeWalkCountPreviousOutput
    double precision                                        :: endTimeThisNode,earliestTimeInTree
    logical                                                 :: interrupted,didEvolve,treeIsDeadlocked
    character(len=35)                                       :: message

    ! Check if this routine is initialized.
    !$omp critical(Merger_Tree_Evolve_To_Initialize)
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
       !@ </inputParameter>
       call Get_Input_Parameter('allTreesExistAtFinalTime',allTreesExistAtFinalTime,defaultValue=.true.)
       !@ <inputParameter>
       !@   <name>mergerTreesDumpStructure</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether merger tree structure should be dumped to a \href{http://www.graphviz.org/}{\sc dot} file.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreesDumpStructure',mergerTreesDumpStructure,defaultValue=.false.)
  
       ! Flag that this routine is now initialized.
       mergerTreeEvolveToInitialized=.true.
    end if
    !$omp end critical(Merger_Tree_Evolve_To_Initialize)

    ! Call routines to perform initializations which must occur for all threads if run in parallel.
    !# <include directive="mergerTreeEvolveThreadInitialize" type="code" action="subroutine">
    include 'merger_trees.evolve.threadInitialize.inc'
    !# </include>

    ! Initialize the tree if necessary.
    call Merger_Tree_Initialize(thisTree)

    ! Check that the output time is not after the end time of this tree.
    if (endTime > Tree_Node_Time(thisTree%baseNode)) then
       ! Final time is exceeded. Check if by a significant factor.
       if (endTime > Tree_Node_Time(thisTree%baseNode)*(1.0d0+timeTolerance)) then
          ! Exceeded by a significant factor - report an error. Check if such behavior is expected.
          if (allTreesExistAtFinalTime) then
             ! It is not, write an error and exit.
             call Galacticus_Error_Report('Merger_Tree_Evolve_To','requested time exceeds the final time in the tree')
          else
             ! It is, so simply ignore this tree.
             return
          end if
       else
          ! Not exceeded by a significant factor (can happen due to approximation errors). Simply reset to actual time requested.
          call Tree_Node_Time_Set(thisTree%baseNode,endTime)
       end if
    end if
    
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
       earliestTimeInTree=largeTime

       ! Point to the base of the tree.
       thisNode => thisTree%baseNode

       ! Set the deadlock flag to true initially.
       treeIsDeadlocked=.true.

       ! Check that tree extends beyond end time.
       if (Tree_Node_Time(thisNode) < endTime) call Galacticus_Error_Report('Merger_Tree_Evolve_To','merger tree stops before&
            & requested end time')

       ! Tree walk loop: Walk to each node in the tree and consider whether or not to evolve it.
       treeWalkLoop: do while (associated(thisNode))

          ! Find the next node that we will process.
          call thisNode%walkTreeWithSatellites(nextNode)

          ! Evolve this node if it exists before the output time and has no children (i.e. they've already all been processed).
          evolveCondition: if (.not.associated(thisNode%childNode) .and. Tree_Node_Time(thisNode) < endTime) then

             ! Flag that a node was evolved.
             didEvolve=.true.

             ! Update tree progress counter.
             nodesEvolvedCount=nodesEvolvedCount+1

             ! Dump the merger tree structure for later plotting.
             if (mergerTreesDumpStructure) call Merger_Tree_Dump(thisTree%index,thisTree%baseNode,[thisNode%index()])

             ! Evolve the node, handling interrupt events. We keep on evolving it until no interrupt is returned (in which case
             ! the node has reached the requested end time) or the node no longer exists (e.g. if it was destroyed).
             interrupted=.true.
             do while (interrupted.and.associated(thisNode))
                interrupted=.false.

                ! Find maximum allowed end time for this particular node.
                endTimeThisNode=Evolve_To_Time(thisNode,endTime,End_Of_Timestep_Task)
             
                ! If this node is able to  evolve by a finite amount, the tree is not deadlocked.
                if (endTimeThisNode > Tree_Node_Time(thisNode)) treeIsDeadlocked=.false.

                ! Update record of earliest time in the tree.
                earliestTimeInTree=min(earliestTimeInTree,endTimeThisNode)

                ! Evolve the node to the next interrupt event, or the end time.
                call thisTree%evolveNode(thisNode,endTimeThisNode,interrupted,interruptProcedure)

                ! Check for interrupt.
                if (interrupted) then
                   ! If an interrupt occured call the specified procedure to handle it.
                   call interruptProcedure(thisNode)
                else                   
                   ! Call routine to handle end of timestep processing.
                   if (associated(End_Of_Timestep_Task)) call End_Of_Timestep_Task(thisTree,thisNode)
                end if
             end do

             ! If this halo has reached its parent halo, decide how to handle it.
             if (associated(thisNode)) then
                if (Tree_Node_Time(thisNode) >= Tree_Node_Time(thisNode%parentNode)) then
                   ! Parent halo has been reached. Check if the node is the primary (major) progenitor of the parent node.
                   select case (thisNode%isPrimaryProgenitor())
                   case (.false.)
                      ! It is not the major progenitor, so this could be a halo merger event unless the halo is already a
                      ! satellite. Check for satellite status and, if it's not a satellite, process this halo merging event.
                      if (.not.thisNode%isSatellite()) call thisTree%mergeNode(thisNode)
                   case (.true.)
                      ! This is the major progenitor, so promote the node to its parent as it is the main progenitor providing
                      ! that the node has no siblings - this ensures that any siblings have already been evolved and become
                      ! satellites of the parent halo.
                      if (.not.associated(thisNode%siblingNode)) call thisTree%promoteNode(thisNode)
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
             write (message,'(a,i9   )') 'Nodes evolved:         ',nodesEvolvedCount
             call Galacticus_Display_Message(message,verbosityLevel)
             write (message,'(a,e10.4)') 'Earliest time in tree: ',earliestTimeInTree
             call Galacticus_Display_Message(message,verbosityLevel)
             call Galacticus_Display_Unindent('done',verbosityLevel)
             treeWalkCountPreviousOutput=treeWalkCount
          end if
       end if

       ! If some nodes were potentially evolvable, but none actually did evolve, then the tree must be deadlocked.
       if (didEvolve.and.treeIsDeadlocked) call Galacticus_Error_Report('Merger_Tree_Evolve_To','merger tree appears to be deadlocked - check timestep criteria')
       
    end do outerLoop

    return
  end subroutine Merger_Tree_Evolve_To

  double precision function Evolve_To_Time(thisNode,endTime,End_Of_Timestep_Task)
    !% Determine the time to which {\tt thisNode} should be evolved.
    use Merger_Tree_Timesteps
    use Cosmology_Functions
    use Input_Parameters
    use Galacticus_Error
    use Galacticus_Display
    use ISO_Varying_String
    use String_Handling
    use Merger_Trees
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: endTime
    type(treeNode),                  pointer :: satelliteNode
    procedure(),      intent(out),   pointer :: End_Of_Timestep_Task
    procedure(),                     pointer :: End_Of_Timestep_Task_Internal
    double precision                         :: time,expansionFactor,expansionTimescale,hostTimeLimit
    character(len=9)                         :: timeFormatted
    type(varying_string)                     :: message

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
       !@ </inputParameter>
       call Get_Input_Parameter('timestepHostRelative',timestepHostRelative,defaultValue=0.1d0)
       !@ <inputParameter>
       !@   <name>timestepHostAbsolute</name>
       !@   <defaultValue>1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The maximum allowed absolute timestep (in Gyr) for node evolution relative to the time of the host halo.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('timestepHostAbsolute',timestepHostAbsolute,defaultValue=1.0d0)       
       evolveToTimeInitialized=.true.
    end if
    !$omp end critical (evolveToTimeInitialize)

    ! Initially set to the global end time.
    Evolve_To_Time=endTime

    ! Limit time based on satellite status.
    select case (thisNode%isSatellite())
    case (.false.)
       ! Limit to the time of its parent node if this node is not a satellite.
       if (associated(thisNode%parentNode)) Evolve_To_Time=min(Evolve_To_Time,Tree_Node_Time(thisNode%parentNode))
    case (.true.)
       ! Do not let satellite evolve too far beyond parent.
       ! Get current cosmic time.
       time=Tree_Node_Time(thisNode%parentNode)
       ! Find current expansion timescale.
       expansionFactor=Expansion_Factor(time)
       expansionTimescale=1.0d0/Expansion_Rate(expansionFactor)
       ! Determine suitable timestep.
       hostTimeLimit=time+min(timestepHostRelative*expansionTimescale,timestepHostAbsolute)
       ! Limit to this time.
       Evolve_To_Time=min(Evolve_To_Time,hostTimeLimit)
    end select

    ! Also ensure that this node is not evolved beyond the time of any of its current satellites.
    satelliteNode => thisNode%satelliteNode
    do while (associated(satelliteNode))
       Evolve_To_Time=min(Evolve_To_Time,Tree_Node_Time(satelliteNode))
       satelliteNode => satelliteNode%siblingNode
    end do

    ! Also ensure that this node is not evolved beyond the time at which any of its mergees merge. In some cases, the node may
    ! already be in the future of a mergee. In such cases, simply freeze it at the current time.
    satelliteNode => thisNode%mergeeNode
    time=Tree_Node_Time(thisNode)
    do while (associated(satelliteNode))
        Evolve_To_Time=min(Evolve_To_Time,max(Tree_Node_Satellite_Time_of_Merging(satelliteNode),time))
        satelliteNode => satelliteNode%nextMergee
    end do

    ! Also ensure that a primary progenitor does not evolve in advance of siblings. This is important since we can not promote a
    ! primary progenitor into its parent until all siblings have become satellites in that parent.
    if (thisNode%isPrimaryProgenitor().and.associated(thisNode%siblingNode)) Evolve_To_Time=min(Evolve_To_Time&
         &,max(Tree_Node_Time(thisNode),Tree_Node_Time(thisNode%siblingNode)))

    ! Also ensure that the timestep taken does not exceed the allowed timestep for this specific node.
    End_Of_Timestep_Task_Internal => null()
    Evolve_To_Time=min(Evolve_To_Time,Tree_Node_Time(thisNode)+Time_Step_Get(thisNode,Evolve_To_Time,End_Of_Timestep_Task_Internal))
    End_Of_Timestep_Task => End_Of_Timestep_Task_Internal

    ! Check that end time exceeds current time.
    if (Evolve_To_Time < Tree_Node_Time(thisNode)) then
       message='end time ('
       write (timeFormatted,'(f7.4)') Evolve_To_Time
       message=message//trim(timeFormatted)//' Gyr) is before current time ('
       write (timeFormatted,'(f7.4)') Tree_Node_Time(thisNode)
       message=message//trim(timeFormatted)//' Gyr) of node '
       message=message//thisNode%index()
       message=message//' (time difference is '
       write (timeFormatted,'(e8.2)') Tree_Node_Time(thisNode)-Evolve_To_Time
       message=message//trim(timeFormatted)
       if (.not.Tree_Node_Is_Accurate(Tree_Node_Time(thisNode),Evolve_To_Time)) then
          ! End time is well before current time. This is an error.
          message=message//' Gyr)'
          call Galacticus_Error_Report('Evolve_To_Time',message)
       else
          ! End time is before current time, but only by a small amount, simply reset the current time to the end time.
          message=message//' Gyr) - this should happen infrequently'
          call Galacticus_Display_Message(message,verbosityInfo)
          call Tree_Node_Time_Set(thisNode,Evolve_To_Time)
       end if
    end if
    return
  end function Evolve_To_Time

end module Merger_Trees_Evolve
