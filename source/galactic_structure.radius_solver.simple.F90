!% Contains a module which implements a simple galactic radii solver (no adiabatic contraction and no self-gravity of baryons).

module Galactic_Structure_Radii_Simple
  !% Implements a simple galactic radii solver (no adiabatic contraction and no self-gravity of baryons).
  use Tree_Nodes
  use Galactic_Structure_Radius_Solver_Procedures
  private
  public :: Galactic_Structure_Radii_Simple_Initialize

contains

  !# <galacticStructureRadiusSolverMethod>
  !#  <unitName>Galactic_Structure_Radii_Simple_Initialize</unitName>
  !# </galacticStructureRadiusSolverMethod>
  subroutine Galactic_Structure_Radii_Simple_Initialize(galacticStructureRadiusSolverMethod,Galactic_Structure_Radii_Solve_Do)
    !% Initializes the ``simple'' galactic radii solver module.
    use ISO_Varying_String
    implicit none
    type(varying_string),          intent(in)    :: galacticStructureRadiusSolverMethod
    procedure(),          pointer, intent(inout) :: Galactic_Structure_Radii_Solve_Do
    
    if (galacticStructureRadiusSolverMethod == 'simple') Galactic_Structure_Radii_Solve_Do =>&
         & Galactic_Structure_Radii_Solve_Simple
    return
  end subroutine Galactic_Structure_Radii_Simple_Initialize

  subroutine Galactic_Structure_Radii_Solve_Simple(thisNode)
    !% Find the radii of galactic components in {\tt thisNode} using the ``simple'' method.
    use Tree_Nodes
    !# <include directive="radiusSolverTask" type="moduleUse">
    include 'galactic_structure.radius_solver.tasks.modules.inc'
    !# </include>
    !# <include directive="radiusSolverPlausibility" type="moduleUse">
    include 'galactic_structure.radius_solver.plausible.modules.inc'
    !# </include>
    implicit none
    type(treeNode),          intent(inout), pointer :: thisNode
    logical                                         :: componentActive,galaxyIsPhysicallyPlausible
    double precision                                :: specificAngularMomentum

    ! Check that the galaxy is physical plausible. In this simple solver, we don't act on this.
    galaxyIsPhysicallyPlausible=.true.
    !# <include directive="radiusSolverPlausibility" type="code" action="subroutine">
    !#  <subroutineArgs>thisNode,galaxyIsPhysicallyPlausible</subroutineArgs>
    include 'galactic_structure.radius_solver.plausible.inc'
    !# </include>

    !# <include directive="radiusSolverTask" type="code" action="subroutine">
    !#  <subroutineArgs>thisNode,componentActive,specificAngularMomentum,Radius_Get,Radius_Set,Velocity_Get,Velocity_Set</subroutineArgs>
    !#  <subroutineAction>if (componentActive) call Solve_For_Radius(thisNode,specificAngularMomentum)</subroutineAction>
    include 'galactic_structure.radius_solver.tasks.inc'
    !# </include>

    return
  end subroutine Galactic_Structure_Radii_Solve_Simple
  
  subroutine Solve_For_Radius(thisNode,specificAngularMomentum)
    !% Solve for the equilibrium radius of the given component.
    use Dark_Matter_Profiles
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: specificAngularMomentum
    double precision                         :: radius,velocity

    ! Find the radius in the dark matter profile with the required specific angular momentum
    radius=Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum(thisNode,specificAngularMomentum)

    ! Find the velocity at this radius.
    velocity=Dark_Matter_Profile_Circular_Velocity(thisNode,radius)

    ! Set the component size to new radius and velocity.
    call Radius_Set  (thisNode,radius  )
    call Velocity_Set(thisNode,velocity)
    return
  end subroutine Solve_For_Radius

end module Galactic_Structure_Radii_Simple
