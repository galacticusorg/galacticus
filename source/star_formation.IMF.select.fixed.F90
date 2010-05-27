!% Contains a module which implements a fixed choice of stellar initial mass function.

module Star_Formation_IMF_Select_Fixed
  !% Implements a fixed choice of stellar initial mass function.
  private
  public :: IMF_Select_Fixed_Initialize

  ! Store the index of the selected IMF.
  integer :: imfSelectedIndex

contains

  !# <imfSelectionMethod>
  !#  <unitName>IMF_Select_Fixed_Initialize</unitName>
  !# </imfSelectionMethod>
  subroutine IMF_Select_Fixed_Initialize(imfSelectionMethod,IMF_Select,imfNames)
    !% Initializes the ``fixed'' IMF selection module.
    use ISO_Varying_String
    use Input_Parameters
    use Star_Formation_IMF_Utilities
    implicit none
    type(varying_string),          intent(in)    :: imfSelectionMethod,imfNames(:)
    procedure(),          pointer, intent(inout) :: IMF_Select
    type(varying_string)                         :: imfSelectionFixed
    
    if (imfSelectionMethod == 'fixed') then
       IMF_Select => IMF_Select_Fixed
       ! Get IMF choice.
       !@ <inputParameter>
       !@   <name>imfSelectionFixed</name>
       !@   <defaultValue>Salpeter</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the initial mass function to use in the ``fixed initial mass function'' module.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('imfSelectionFixed',imfSelectionFixed,defaultValue='Salpeter')
       imfSelectedIndex=IMF_Index_Lookup(imfSelectionFixed,imfNames)
    end if
    return
  end subroutine IMF_Select_Fixed_Initialize

  integer function IMF_Select_Fixed(starFormationRate,fuelAbundances)
    !% Return our selection of stellar initial mass function.
    use Abundances_Structure
    implicit none
    double precision,          intent(in) :: starFormationRate
    type(abundancesStructure), intent(in) :: fuelAbundances
    
    IMF_Select_Fixed=imfSelectedIndex
    return
  end function IMF_Select_Fixed
  
end module Star_Formation_IMF_Select_Fixed
