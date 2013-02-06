!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    Galacticus is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

!% Contains a module of useful utilities required by the IMF subsystem.

module Star_Formation_IMF_Utilities
  !% Contains useful utilities required by the IMF subsystem.
  implicit none
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
