!% Contains a module which implements calculations of bar instability in galactic disks.

module Galactic_Dynamics_Bar_Instabilities
  !% Implements calculations of bar instability in galactic disks.
  use ISO_Varying_String
  use Tree_Nodes
  private
  public :: Bar_Instability_Timescale

  ! Flag to indicate if this module has been initialized.  
  logical              :: barInstabilitiesInitialized=.false.

  ! Name of cooling rate available method used.
  type(varying_string) :: barInstabilityMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Bar_Instability_Template), pointer :: Bar_Instability_Timescale_Get => null()
  abstract interface
     double precision function Bar_Instability_Template(thisNode)
       import treeNode
       type(treeNode), intent(inout), pointer :: thisNode
     end function Bar_Instability_Template
  end interface

contains

  subroutine Galactic_Dynamics_Bar_Instability_Initialize
    !% Initialize the bar instability module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="barInstabilityMethod" type="moduleUse">
    include 'galactic_dynamics.bar_instability.modules.inc'
    !# </include>
    implicit none

    !$omp critical(Galactic_Dynamics_Bar_Instability_Initialize) 
    ! Initialize if necessary.
    if (.not.barInstabilitiesInitialized) then
       ! Get the halo spin distribution method parameter.
       !@ <inputParameter>
       !@   <name>barInstabilityMethod</name>
       !@   <defaultValue>null</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for bar instability calculations.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('barInstabilityMethod',barInstabilityMethod,defaultValue='null')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="barInstabilityMethod" type="code" action="subroutine">
       !#  <subroutineArgs>barInstabilityMethod,Bar_Instability_Timescale_Get</subroutineArgs>
       include 'galactic_dynamics.bar_instability.inc'
       !# </include>
       if (.not.associated(Bar_Instability_Timescale_Get)) &
            & call Galacticus_Error_Report('Galactic_Dynamics_Bar_Instability_Initialize','method ' //char(barInstabilityMethod)//' is unrecognized')
       barInstabilitiesInitialized=.true.
    end if
    !$omp end critical(Galactic_Dynamics_Bar_Instability_Initialize) 

    return
  end subroutine Galactic_Dynamics_Bar_Instability_Initialize

  double precision function Bar_Instability_Timescale(thisNode)
    !% Returns a timescale on which the bar instability depletes material from a disk into a pseudo-bulge. A negative value
    !% indicates no instability.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode

    ! Initialize the module.
    call Galactic_Dynamics_Bar_Instability_Initialize

    ! Get the timescale using the selected method.
    Bar_Instability_Timescale=Bar_Instability_Timescale_Get(thisNode)

    return
  end function Bar_Instability_Timescale

end module Galactic_Dynamics_Bar_Instabilities
