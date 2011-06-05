!% Contains a module which defines the template for tasks performed at the end of timesteps.

module Merger_Trees_Evolve_Timesteps_Template
  !% Defines the template for tasks performed at the end of timesteps.
  use Merger_Trees
  use Tree_Nodes
  public

  interface End_Of_Timestep_Task_Interface
     subroutine End_Of_Timestep_Task_Template(thisTree,thisNode)
       import mergerTree, treeNode
       type(mergerTree), intent(in)             :: thisTree
       type(treeNode),   intent(inout), pointer :: thisNode
     end subroutine End_Of_Timestep_Task_Template
  end interface

end module Merger_Trees_Evolve_Timesteps_Template
