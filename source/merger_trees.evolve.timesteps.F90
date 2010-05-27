!% Contains a module which implements calculations of timesteps for merger tree evolution.

module Merger_Tree_Timesteps
  !% Implements calculations of timesteps for merger tree evolution.
  private
  public :: Time_Step_Get

contains

  double precision function Time_Step_Get(thisNode,evolveToTime,End_Of_Timestep_Task)
    !% Computes a suitable timestep over which to evolve a node in a tree.
    use Tree_Nodes
    use Merger_Trees
    use Input_Parameters
    use Galacticus_Error
    use Tree_Node_Methods
    use Merger_Trees_Evolve_Timesteps_Template
    !# <include directive="timeStepsTask" type="moduleUse">
    include 'merger_trees.evolve.timesteps.moduleUse.inc'
    !# </include>
    implicit none
    type(treeNode),                           intent(inout), pointer :: thisNode
    double precision,                         intent(in)             :: evolveToTime
    procedure(),                              intent(out),   pointer :: End_Of_Timestep_Task
    procedure(End_Of_Timestep_Task_Template),                pointer :: End_Of_Timestep_Task_Internal

    ! Call the function to get the timestep.
    Time_Step_Get=evolveToTime-Tree_Node_Time(thisNode)
    End_Of_Timestep_Task_Internal => null()
    !# <include directive="timeStepsTask" type="code" action="subroutine">
    !#  <subroutineArgs>thisNode,Time_Step_Get,End_Of_Timestep_Task_Internal</subroutineArgs>
    include 'merger_trees.evolve.timesteps.inc'
    !# </include>
    End_Of_Timestep_Task => End_Of_Timestep_Task_Internal
    return
  end function Time_Step_Get

end module Merger_Tree_Timesteps
