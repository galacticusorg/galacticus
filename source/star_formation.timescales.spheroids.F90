!% Contains a module which implements calculations of star formation timescales for galactic spheroids.

module Star_Formation_Timescales_Spheroids
  !% Implements calculations of star formation timescales for galactic spheroids.
  use ISO_Varying_String
  use Tree_Nodes
  private
  public :: Star_Formation_Timescale_Spheroid
  
  ! Flag to indicate if this module has been initialized.  
  logical              :: starFormationTimescaleSpheroidsInitialized=.false.

  ! Name of cooling rate available method used.
  type(varying_string) :: starFormationTimescaleSpheroidsMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Star_Formation_Timescale_Spheroid_Template), pointer :: Star_Formation_Timescale_Spheroid_Get => null()
  abstract interface
     double precision function Star_Formation_Timescale_Spheroid_Template(thisNode)
       import treeNode
       type(treeNode), intent(inout), pointer :: thisNode
     end function Star_Formation_Timescale_Spheroid_Template
  end interface

contains

  subroutine Star_Formation_Timescale_Spheroids_Initialize
    !% Initialize the spheroid star formation timecale module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="starFormationTimescaleSpheroidsMethod" type="moduleUse">
    include 'star_formation.timescales.spheroids.modules.inc'
    !# </include>
    implicit none

    !$omp critical(Star_Formation_Timescale_Spheroids_Initialization) 
    ! Initialize if necessary.
    if (.not.starFormationTimescaleSpheroidsInitialized) then
       ! Get the spheroid star formation timescale method parameter.
       !@ <inputParameter>
       !@   <name>starFormationTimescaleSpheroidsMethod</name>
       !@   <defaultValue>dynamical time</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for computing star formation timescales in spheroids.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationTimescaleSpheroidsMethod',starFormationTimescaleSpheroidsMethod,defaultValue='dynamical time')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="starFormationTimescaleSpheroidsMethod" type="code" action="subroutine">
       !#  <subroutineArgs>starFormationTimescaleSpheroidsMethod,Star_Formation_Timescale_Spheroid_Get</subroutineArgs>
       include 'star_formation.timescales.spheroids.inc'
       !# </include>
       if (.not.associated(Star_Formation_Timescale_Spheroid_Get)) call Galacticus_Error_Report('Star_Formation_Timescale_Spheroids'&
            &,'method ' //char(starFormationTimescaleSpheroidsMethod)//' is unrecognized')
       starFormationTimescaleSpheroidsInitialized=.true.
    end if
    !$omp end critical(Star_Formation_Timescale_Spheroids_Initialization) 

    return
  end subroutine Star_Formation_Timescale_Spheroids_Initialize

  double precision function Star_Formation_Timescale_Spheroid(thisNode)
    !% Returns the timescale (in Gyr) for star formation in the spheroid component of {\tt thisNode}.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode

    ! Initialize the module.
    call Star_Formation_Timescale_Spheroids_Initialize

    ! Get the energy using the selected method.
    Star_Formation_Timescale_Spheroid=Star_Formation_Timescale_Spheroid_Get(thisNode)

    return
  end function Star_Formation_Timescale_Spheroid
  
end module Star_Formation_Timescales_Spheroids
