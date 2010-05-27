!% Contains a module of useful utilities required by the IMF subsystem.

module Star_Formation_IMF_Utilities
  !% Contains useful utilities required by the IMF subsystem.
  private
  public :: IMF_Index_Lookup

contains

  integer function IMF_Index_Lookup(imfSelection,imfNames)
    !% Returns the internal index of a stellar initial mass function specified by name via {\tt imfSelection}.
    use ISO_Varying_String
    use Galacticus_Error
    implicit none
    type(varying_string), intent(in) :: imfSelection,imfNames(:)
    
    IMF_Index_Lookup=1
    do while (IMF_Index_Lookup < size(imfNames))
       if (imfSelection == imfNames(IMF_Index_Lookup)) exit
       IMF_Index_Lookup=IMF_Index_Lookup+1
    end do
    if (IMF_Index_Lookup > size(imfNames)) call Galacticus_Error_Report('IMF_Index_Lookup','unmatched IMF name')
    return
  end function IMF_Index_Lookup
  
end module Star_Formation_IMF_Utilities
