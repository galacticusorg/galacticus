!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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
  private
  public :: Galacticus_Task_Evolve_Tree

  ! Flag to indicate if output times have been initialized.
  logical :: outputsInitialized=.false.

  ! Array of output times.
  integer                                     :: outputCount
  double precision, allocatable, dimension(:) :: outputTimes
  
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
    use Cosmology_Functions
    use Merger_Trees
    use Merger_Tree_Construction
    use Merger_Trees_Evolve
    use Tree_Nodes
    use Tree_Node_Methods
    use Galacticus_Output_Merger_Tree
    use Galacticus_Display
    use Input_Parameters
    use Sort
    use Memory_Management
    use Histories
    ! Include modules needed for pre-evolution tasks.
    !# <include directive="mergerTreePreEvolveTask" type="moduleUse">
    include 'galacticus.tasks.evolve_tree.preEvolveTask.moduleUse.inc'
    !# </include>
    implicit none
    type(mergerTree),     save, pointer :: thisTree
    logical,              save          :: finished
    integer,              save          :: iOutput
    double precision,     save          :: outputTime
    type(varying_string), save          :: message
    character(len=20),    save          :: label
    double precision                    :: aExpansion
    !$omp threadprivate(thisTree,finished,iOutput,outputTime,message,label)

    ! Initialize the task if necessary.
    !$omp critical (Tasks_Evolve_Tree_Initialize)
    if (.not.outputsInitialized) then

       ! Get a list of output redshifts - stored temporarily in the outputTimes array.
       outputCount=max(Get_Input_Parameter_Array_Size('outputRedshifts'),1)
       call Alloc_Array(outputTimes,outputCount,'outputTimes')
       !@ <inputParameter>
       !@   <name>outputRedshifts</name>
       !@   <defaultValue>0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     A list of redshifts at which \glc\ results should be output.
       !@   </description>
       !@ </inputParameter>
       if (outputCount == 1) then
          ! If only one (or zero) output redshifts present, make redshift zero the default.
          call Get_Input_Parameter('outputRedshifts',outputTimes,defaultValue=[0.0d0])
       else
          call Get_Input_Parameter('outputRedshifts',outputTimes)
       end if

       ! Convert redshifts to times.
       do iOutput=1,outputCount
          aExpansion=Expansion_Factor_from_Redshift(outputTimes(iOutput))
          outputTimes(iOutput)=Cosmology_Age(aExpansion)
       end do

       ! Sort the times.
       call Sort_Do(outputTimes)

       ! Set history ranges to include these times.
       call History_Set_Times(timeEarliest=outputTimes(1),timeLatest=outputTimes(outputCount))

       ! Flag that this task is now initialized.
       outputsInitialized=.true.
    end if
    !$omp end critical (Tasks_Evolve_Tree_Initialize)

    ! Begin looping through available trees.
    finished=.false.
    !$omp parallel copyin(finished)
    do while (.not.finished)
       thisTree => Merger_Tree_Create()
       finished=finished.or..not.associated(thisTree)
       if (.not.finished) then
          
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
          outputTimeLoop : do iOutput=1,outputCount

             ! Get the output time.
             outputTime=outputTimes(iOutput)

             ! Evolve the tree to the output time.
             call Merger_Tree_Evolve_To(thisTree,outputTime)
          
             ! Output the merger tree.
             write (label,'(f7.2)') outputTime
             message="Output tree data at t="//trim(label)//" Gyr"          
             call Galacticus_Display_Message(message)
             call Galacticus_Merger_Tree_Output(thisTree,iOutput,outputTime,.false.)
          
          end do outputTimeLoop
             
          ! Destroy the tree.
          call thisTree%destroy()

          ! Unindent messages.
          call Galacticus_Display_Unindent('Finished tree')

       end if
    end do
    !$omp end parallel

    Galacticus_Task_Evolve_Tree=.false.
    return
  end function Galacticus_Task_Evolve_Tree

end module Galacticus_Tasks_Evolve_Tree
