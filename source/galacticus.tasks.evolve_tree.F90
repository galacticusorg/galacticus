!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements the task of evolving merger trees.

module Galacticus_Tasks_Evolve_Tree
  !% Implements the task of evolving merger trees.
  implicit none
  private
  public :: Galacticus_Task_Evolve_Tree

  ! Flag to indicate if output times have been initialized.
  logical :: treeEvolveInitialized=.false.

  ! Parameters controlling which trees will be processed.
  integer :: treeEvolveWorkerCount,treeEvolveWorkerNumber
  
contains

  !# <galacticusTask>
  !#  <unitName>Galacticus_Task_Evolve_Tree</unitName>
  !#  <after>Galacticus_Task_Start</after>
  !#  <after>Galacticus_Task_Test</after>
  !#  <before>Galacticus_Task_End</before>
  !# </galacticusTask>
  logical function Galacticus_Task_Evolve_Tree()
    !% Evolves the complete set of merger trees as specified.
    use ISO_Varying_String
    use String_Handling
    use Merger_Trees
    use Merger_Tree_Construction
    use Merger_Trees_Evolve
    use Tree_Nodes
    use Galacticus_Output_Merger_Tree
    use Galacticus_Display
    use Input_Parameters
    use Galacticus_Output_Times
    use Merger_Tree_Active
    use Memory_Management
    ! Include modules needed for pre- and post-evolution and pre-construction tasks.
    !# <include directive="mergerTreePreEvolveTask" type="moduleUse">
    include 'galacticus.tasks.evolve_tree.preEvolveTask.moduleUse.inc'
    !# </include>
    !# <include directive="mergerTreePostEvolveTask" type="moduleUse">
    include 'galacticus.tasks.evolve_tree.postEvolveTask.moduleUse.inc'
    !# </include>
    !# <include directive="mergerTreePreTreeConstructionTask" type="moduleUse">
    include 'galacticus.tasks.evolve_tree.preConstructionTask.moduleUse.inc'
    !# </include>
    implicit none
    type(mergerTree),     save, pointer :: thisTree
    logical,              save          :: finished,skipTree
    integer,              save          :: iOutput
    double precision,     save          :: outputTime
    type(varying_string), save          :: message
    character(len=20),    save          :: label
    !$omp threadprivate(thisTree,finished,skipTree,iOutput,outputTime,message,label)
    integer                             :: iTree

    ! Initialize the task if necessary.
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

       ! Flag that this task is now initialized.
       treeEvolveInitialized=.true.
    end if
    !$omp end critical (Tasks_Evolve_Tree_Initialize)

    ! Begin looping through available trees.
    finished=.false.
    iTree=0

    !$omp parallel copyin(finished)
    do while (.not.finished)
       ! Increment the tree number.
       !$omp atomic
       iTree=iTree+1
       ! Decide whether or not to skip this tree.
       skipTree=.not.(modulo(iTree-1+(iTree-1)/treeEvolveWorkerCount,treeEvolveWorkerCount) == treeEvolveWorkerNumber-1)
       ! Perform any pre-tree construction tasks.
       !# <include directive="mergerTreePreTreeConstructionTask" type="code" action="subroutine">
       include 'galacticus.tasks.evolve_tree.preConstructionTask.inc'
       !# </include>
       
       ! Get a tree.
       !$omp critical(Tree_Sharing)
       thisTree => Merger_Tree_Create(skipTree)
       finished=finished.or..not.associated(thisTree)
       !$omp end critical(Tree_Sharing)
       if (.not.finished) then
          
          ! Skip this tree if necessary.
          if (.not.skipTree) then

             ! Set this as the active tree.
             activeTreeWeight=thisTree%volumeWeight

             ! Perform any pre-evolution tasks on the tree.
             !# <include directive="mergerTreePreEvolveTask" type="code" action="subroutine">
             !#  <subroutineArgs>thisTree</subroutineArgs>
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
          !# <include directive="mergerTreePostEvolveTask" type="code" action="subroutine">
          include 'galacticus.tasks.evolve_tree.postEvolveTask.inc'
          !# </include>
             
       end if

    end do
    !$omp end parallel

    Galacticus_Task_Evolve_Tree=.false.
    return
  end function Galacticus_Task_Evolve_Tree

end module Galacticus_Tasks_Evolve_Tree
