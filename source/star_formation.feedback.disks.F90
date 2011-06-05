!% Contains a module which implements calculations of feedback from star formation in disks.

module Star_Formation_Feedback_Disks
  !% Implements calculations of feedback from star formation in disks.
  use ISO_Varying_String
  use Tree_Nodes
  private
  public :: Star_Formation_Feedback_Disk_Outflow_Rate
  
  ! Flag to indicate if this module has been initialized.  
  logical              :: starFormationFeedbackDisksInitialized=.false.

  ! Name of cooling rate available method used.
  type(varying_string) :: starFormationFeedbackDisksMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Star_Formation_Feedback_Disk_Outflow_Rate_Template), pointer :: Star_Formation_Feedback_Disk_Outflow_Rate_Get => null()
  interface Star_Formation_Feedback_Disk_Outflow_Rate_Template
     double precision function Star_Formation_Feedback_Disk_Outflow_Rate_Template(thisNode,starFormationRate,energyInputRate)
       import treeNode
       type(treeNode),   intent(inout), pointer :: thisNode
       double precision, intent(in)             :: starFormationRate,energyInputRate
     end function Star_Formation_Feedback_Disk_Outflow_Rate_Template
  end interface

contains

  subroutine Star_Formation_Feedback_Disks_Initialize
    !% Initialize the disk star formation feedback module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="starFormationFeedbackDisksMethod" type="moduleUse">
    include 'star_formation.feedbacks.disks.modules.inc'
    !# </include>
    implicit none

    !$omp critical(Star_Formation_Feedback_Disks_Initialization) 
    ! Initialize if necessary.
    if (.not.starFormationFeedbackDisksInitialized) then
       ! Get the disk star formation feedback method parameter.
       !@ <inputParameter>
       !@   <name>starFormationFeedbackDisksMethod</name>
       !@   <defaultValue>power law</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for calculations of \SNe\ feedback in disks.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationFeedbackDisksMethod',starFormationFeedbackDisksMethod,defaultValue='power law')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="starFormationFeedbackDisksMethod" type="code" action="subroutine">
       !#  <subroutineArgs>starFormationFeedbackDisksMethod,Star_Formation_Feedback_Disk_Outflow_Rate_Get</subroutineArgs>
       include 'star_formation.feedbacks.disks.inc'
       !# </include>
       if (.not.associated(Star_Formation_Feedback_Disk_Outflow_Rate_Get)) call Galacticus_Error_Report('Star_Formation_Feedback_Disks'&
            &,'method ' //char(starFormationFeedbackDisksMethod)//' is unrecognized')
       starFormationFeedbackDisksInitialized=.true.
    end if
    !$omp end critical(Star_Formation_Feedback_Disks_Initialization) 

    return
  end subroutine Star_Formation_Feedback_Disks_Initialize

  double precision function Star_Formation_Feedback_Disk_Outflow_Rate(thisNode,starFormationRate,energyInputRate)
    !% Returns the outflow rate due to star formation in the disk component of {\tt thisNode}.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: starFormationRate,energyInputRate

    ! Initialize the module.
    call Star_Formation_Feedback_Disks_Initialize

    ! Get the energy using the selected method.
    Star_Formation_Feedback_Disk_Outflow_Rate=Star_Formation_Feedback_Disk_Outflow_Rate_Get(thisNode,starFormationRate,energyInputRate)

    return
  end function Star_Formation_Feedback_Disk_Outflow_Rate
  
end module Star_Formation_Feedback_Disks
