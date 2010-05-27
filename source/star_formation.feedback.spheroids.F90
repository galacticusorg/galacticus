!% Contains a module which implements calculations of feedback from star formation in spheroids.

module Star_Formation_Feedback_Spheroids
  !% Implements calculations of feedback from star formation in spheroids.
  use ISO_Varying_String
  use Tree_Nodes
  private
  public :: Star_Formation_Feedback_Spheroid_Outflow_Rate
  
  ! Flag to indicate if this module has been initialized.  
  logical              :: starFormationFeedbackSpheroidsInitialized=.false.

  ! Name of cooling rate available method used.
  type(varying_string) :: starFormationFeedbackSpheroidsMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Star_Formation_Feedback_Spheroid_Outflow_Rate_Template), pointer :: Star_Formation_Feedback_Spheroid_Outflow_Rate_Get => null()
  interface Star_Formation_Feedback_Spheroid_Outflow_Rate_Template
     double precision function Star_Formation_Feedback_Spheroid_Outflow_Rate_Template(thisNode,starFormationRate,energyInputRate)
       import treeNode
       type(treeNode),   intent(inout), pointer :: thisNode
       double precision, intent(in)             :: starFormationRate,energyInputRate
     end function Star_Formation_Feedback_Spheroid_Outflow_Rate_Template
  end interface

contains

  subroutine Star_Formation_Feedback_Spheroids_Initialize
    !% Initialize the spheroid star formation feedback module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="starFormationFeedbackSpheroidsMethod" type="moduleUse">
    include 'star_formation.feedbacks.spheroids.modules.inc'
    !# </include>
    implicit none

    !$omp critical(Star_Formation_Feedback_Spheroids_Initialization) 
    ! Initialize if necessary.
    if (.not.starFormationFeedbackSpheroidsInitialized) then
       ! Get the spheroid star formation feedback method parameter.
       !@ <inputParameter>
       !@   <name>starFormationFeedbackSpheroidsMethod</name>
       !@   <defaultValue>power law</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for calculations of \SNe\ feedback in spheroids.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationFeedbackSpheroidsMethod',starFormationFeedbackSpheroidsMethod,defaultValue='power law')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="starFormationFeedbackSpheroidsMethod" type="code" action="subroutine">
       !#  <subroutineArgs>starFormationFeedbackSpheroidsMethod,Star_Formation_Feedback_Spheroid_Outflow_Rate_Get</subroutineArgs>
       include 'star_formation.feedbacks.spheroids.inc'
       !# </include>
       if (.not.associated(Star_Formation_Feedback_Spheroid_Outflow_Rate_Get)) call Galacticus_Error_Report('Star_Formation_Feedback_Spheroids'&
            &,'method ' //char(starFormationFeedbackSpheroidsMethod)//' is unrecognized')
       starFormationFeedbackSpheroidsInitialized=.true.
    end if
    !$omp end critical(Star_Formation_Feedback_Spheroids_Initialization) 

    return
  end subroutine Star_Formation_Feedback_Spheroids_Initialize

  double precision function Star_Formation_Feedback_Spheroid_Outflow_Rate(thisNode,starFormationRate,energyInputRate)
    !% Returns the outflow rate due to star formation in the spheroid component of {\tt thisNode}.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: starFormationRate,energyInputRate

    ! Initialize the module.
    call Star_Formation_Feedback_Spheroids_Initialize

    ! Get the energy using the selected method.
    Star_Formation_Feedback_Spheroid_Outflow_Rate=Star_Formation_Feedback_Spheroid_Outflow_Rate_Get(thisNode,starFormationRate,energyInputRate)

    return
  end function Star_Formation_Feedback_Spheroid_Outflow_Rate
  
end module Star_Formation_Feedback_Spheroids
