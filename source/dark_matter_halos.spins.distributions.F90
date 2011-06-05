!% Contains a module that implements calculations of dark matter halo spin distributions

module Halo_Spin_Distributions
  !% Implements calculations of dark matter halo spin distributions
  use ISO_Varying_String
  use Tree_Nodes
  private
  public :: Halo_Spin_Distribution_Sample

  ! Flag to indicate if this module has been initialized.  
  logical              :: haloSpinDistributionInitialized=.false.

  ! Name of cooling rate available method used.
  type(varying_string) :: haloSpinDistributionMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Halo_Spin_Sample_Get_Template), pointer :: Halo_Spin_Sample_Get => null()
  interface Halo_Spin_Sample_Get_Template
     double precision function Halo_Spin_Sample_Get_Template(thisNode)
       import treeNode
       type(treeNode), intent(inout), pointer :: thisNode
     end function Halo_Spin_Sample_Get_Template
  end interface
  
contains

  double precision function Halo_Spin_Distribution_Sample(thisNode)
    !% Return a halo spin selected randomly from a distribution.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="haloSpinDistributionMethod" type="moduleUse">
    include 'dark_matter_halos.spins.distributions.modules.inc'
    !# </include>
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    !$omp critical(Halo_Spin_Distribution_Initialization) 
    ! Initialize if necessary.
    if (.not.haloSpinDistributionInitialized) then
       ! Get the halo spin distribution method parameter.
       !@ <inputParameter>
       !@   <name>haloSpinDistributionMethod</name>
       !@   <defaultValue>lognormal</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be use for computing halo spin distributions.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('haloSpinDistributionMethod',haloSpinDistributionMethod,defaultValue='lognormal')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="haloSpinDistributionMethod" type="code" action="subroutine">
       !#  <subroutineArgs>haloSpinDistributionMethod,Halo_Spin_Sample_Get</subroutineArgs>
       include 'dark_matter_halos.spins.distributions.inc'
       !# </include>
       if (.not.associated(Halo_Spin_Sample_Get)) call Galacticus_Error_Report('Halo_Spin_Distribution','method ' &
            &//char(haloSpinDistributionMethod)//' is unrecognized')
       haloSpinDistributionInitialized=.true.
    end if
    !$omp end critical(Halo_Spin_Distribution_Initialization) 

    ! Get the cooling rate using the selected method.
    Halo_Spin_Distribution_Sample=Halo_Spin_Sample_Get(thisNode)
    return
  end function Halo_Spin_Distribution_Sample

end module Halo_Spin_Distributions
