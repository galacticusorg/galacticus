!% Contains a module which implements calculations of sizes of galactic components (or more general components).

module Galactic_Structure_Radii
  !% Implements calculations of sizes of galactic components (or more general components).
  use ISO_Varying_String
  use Tree_Nodes
  private
  public :: Galactic_Structure_Radii_Solve

  ! Flag to indicate if this module has been initialized.  
  logical              :: galacticStructureRadiusSolverInitialized=.false.

  ! Name of cooling rate available method used.
  type(varying_string) :: galacticStructureRadiusSolverMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Galactic_Structure_Radii_Solve_Template), pointer :: Galactic_Structure_Radii_Solve_Do => null()
  abstract interface
     subroutine Galactic_Structure_Radii_Solve_Template(thisNode)
       import treeNode
       type(treeNode), intent(inout), pointer :: thisNode
     end subroutine Galactic_Structure_Radii_Solve_Template
  end interface
  
contains

  !# <preDerivativeComputeTask>
  !# <unitName>Galactic_Structure_Radii_Solve</unitName>
  !# </preDerivativeComputeTask>
  subroutine Galactic_Structure_Radii_Solve(thisNode)
    !% Solve for the radii of galactic components in {\tt thisNode}.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="galacticStructureRadiusSolverMethod" type="moduleUse">
    include 'galactic_structure.radius_solver.modules.inc'
    !# </include>
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    !$omp critical(Galactic_Structure_Radii_Initialization) 
    ! Initialize if necessary.
    if (.not.galacticStructureRadiusSolverInitialized) then
       ! Get the galactic structure radii solver method parameter.
       !@ <inputParameter>
       !@   <name>galacticStructureRadiusSolverMethod</name>
       !@   <defaultValue>simple</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Selects the method to be used for solving for galactic structure.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('galacticStructureRadiusSolverMethod',galacticStructureRadiusSolverMethod,defaultValue='simple')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="galacticStructureRadiusSolverMethod" type="code" action="subroutine">
       !#  <subroutineArgs>galacticStructureRadiusSolverMethod,Galactic_Structure_Radii_Solve_Do</subroutineArgs>
       include 'galactic_structure.radius_solver.inc'
       !# </include>
       if (.not.associated(Galactic_Structure_Radii_Solve_Do)) &
            & call Galacticus_Error_Report('Galactic_Structure_Radii','method '//char(galacticStructureRadiusSolverMethod)//' is unrecognized')
       galacticStructureRadiusSolverInitialized=.true.
    end if
    !$omp end critical(Galactic_Structure_Radii_Initialization) 

    ! Call the routine to solve for the radii.
    call Galactic_Structure_Radii_Solve_Do(thisNode)

    return
  end subroutine Galactic_Structure_Radii_Solve

end module Galactic_Structure_Radii
