!% Contains a module which holds procedure pointers used by the galactic structure radii solver subsystem.

module Galactic_Structure_Radius_Solver_Procedures
  !% Holds procedure pointers used by the galactic structure radii solver subsystem.
  use Tree_Nodes
  private
  public :: Radius_Get, Radius_Set, Velocity_Get, Velocity_Set

  ! Pointers to get and set procedures.
  procedure(Get_Template), pointer :: Radius_Get => null(), Velocity_Get => null()
  procedure(Set_Template), pointer :: Radius_Set => null(), Velocity_Set => null()
  !$omp threadprivate(Radius_Get,Radius_Set,Velocity_Get,Velocity_Set)

  abstract interface
     double precision function Get_Template(thisNode)
       import treeNode
       type(treeNode), pointer, intent(inout) :: thisNode
     end function Get_Template
  end interface
  abstract interface
     subroutine Set_Template(thisNode,value)
       import treeNode
       type(treeNode),   pointer, intent(inout) :: thisNode
       double precision,          intent(in)    :: value
     end subroutine Set_Template
  end interface

end module Galactic_Structure_Radius_Solver_Procedures
