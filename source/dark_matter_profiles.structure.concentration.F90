!% Contains a module which implements calculations of dark matter halo density profile concentrations.

module Dark_Matter_Profiles_Concentrations
  !% Implements calculations of dark matter halo density profile concentrations.
  use ISO_Varying_String
  use Tree_Nodes
  private
  public :: Dark_Matter_Profile_Concentration

  ! Flag to indicate if this module has been initialized.  
  logical              :: darkMatterConcentrationInitialized=.false.

  ! Name of cooling rate available method used.
  type(varying_string) :: darkMatterConcentrationMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Dark_Matter_Profile_Template), pointer :: Dark_Matter_Profile_Concentration_Get => null()
  abstract interface
     double precision function Dark_Matter_Profile_Template(thisNode)
       import treeNode
       type(treeNode), intent(inout), pointer :: thisNode
     end function Dark_Matter_Profile_Template
  end interface

contains

  subroutine Dark_Matter_Concentrations_Initialize
    !% Initialize the dark matter profile module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="darkMatterConcentrationMethod" type="moduleUse">
    include 'dark_matter_profiles.structure.concentration.modules.inc'
    !# </include>
    implicit none

    !$omp critical(Dark_Matter_Concentrations_Initialization) 
    ! Initialize if necessary.
    if (.not.darkMatterConcentrationInitialized) then
       ! Get the halo spin distribution method parameter.
       !@ <inputParameter>
       !@   <name>darkMatterConcentrationMethod</name>
       !@   <defaultValue>Gao 2008</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for calculations of dark matter halo density profile concentrations.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('darkMatterConcentrationMethod',darkMatterConcentrationMethod,defaultValue='Gao 2008')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="darkMatterConcentrationMethod" type="code" action="subroutine">
       !#  <subroutineArgs>darkMatterConcentrationMethod,Dark_Matter_Profile_Concentration_Get</subroutineArgs>
       include 'dark_matter_profiles.structure.concentration.inc'
       !# </include>
       if (.not.associated(Dark_Matter_Profile_Concentration_Get)) &
            & call Galacticus_Error_Report('Dark_Matter_Concentrations_Initialize','method ' //char(darkMatterConcentrationMethod)//' is unrecognized')
       darkMatterConcentrationInitialized=.true.
    end if
    !$omp end critical(Dark_Matter_Concentrations_Initialization) 

    return
  end subroutine Dark_Matter_Concentrations_Initialize

  double precision function Dark_Matter_Profile_Concentration(thisNode)
    !% Returns the concentration of the dark matter profile of {\tt thisNode}.
    implicit none
    type(treeNode),pointer, intent(inout) :: thisNode

    ! Initialize the module.
    call Dark_Matter_Concentrations_Initialize

    ! Get the concentration using the selected method.
    Dark_Matter_Profile_Concentration=Dark_Matter_Profile_Concentration_Get(thisNode)

    return
  end function Dark_Matter_Profile_Concentration

end module Dark_Matter_Profiles_Concentrations
