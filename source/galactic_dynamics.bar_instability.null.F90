!% Contains a module which implements a null calculation of bar instability.

module Galactic_Dynamics_Bar_Instabilities_Null
  !% Implements a null calculation of bar instability.
  use Tree_Nodes
  private
  public :: Galactic_Dynamics_Bar_Instabilities_Null_Initialize

contains

  !# <barInstabilityMethod>
  !#  <unitName>Galactic_Dynamics_Bar_Instabilities_Null_Initialize</unitName>
  !# </barInstabilityMethod>
  subroutine Galactic_Dynamics_Bar_Instabilities_Null_Initialize(barInstabilityMethod,Bar_Instability_Timescale_Get)
    !% Initializes the ``Null'' bar instability module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),          intent(in)    :: barInstabilityMethod
    procedure(),          pointer, intent(inout) :: Bar_Instability_Timescale_Get
    
    if (barInstabilityMethod == 'null') Bar_Instability_Timescale_Get => Bar_Instability_Timescale_Null

    return
  end subroutine Galactic_Dynamics_Bar_Instabilities_Null_Initialize

  double precision function Bar_Instability_Timescale_Null(thisNode)
    !% Assumes that disks are never bar unstable and so returns an infinite timescale for bar instability.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    ! Assume infinite timescale (i.e. no instability).
    Bar_Instability_Timescale_Null=-1.0d0

    return
  end function Bar_Instability_Timescale_Null
  
end module Galactic_Dynamics_Bar_Instabilities_Null
