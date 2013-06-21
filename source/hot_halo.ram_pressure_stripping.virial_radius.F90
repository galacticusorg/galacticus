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

!% Contains a module which implements a null hot halo ram pressure stripping calculation, by simply returning the virial radius as
!% the ram pressure stripping radius.

module Hot_Halo_Ram_Pressure_Stripping_Virial_Radii
  !% Implements a null hot halo ram pressure stripping calculation, by simply returning the virial radius as the ram pressure
  !% stripping radius.
  use, intrinsic :: ISO_C_Binding
  implicit none
  private
  public :: Hot_Halo_Ram_Pressure_Stripping_Virial_Radii_Initialize

contains

  !# <hotHaloRamPressureStrippingMethod>
  !#  <unitName>Hot_Halo_Ram_Pressure_Stripping_Virial_Radii_Initialize</unitName>
  !# </hotHaloRamPressureStrippingMethod>
  subroutine Hot_Halo_Ram_Pressure_Stripping_Virial_Radii_Initialize(hotHaloRamPressureStrippingMethod,Hot_Halo_Ram_Pressure_Stripping_Get)
    !% Initializes the ``virial radius'' hot halo ram pressure stripping module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type     (varying_string                               ), intent(in   )          :: hotHaloRamPressureStrippingMethod    
    procedure(Hot_Halo_Ram_Pressure_Stripping_Virial_Radius), intent(inout), pointer :: Hot_Halo_Ram_Pressure_Stripping_Get  
                                                                                                                          
    if (hotHaloRamPressureStrippingMethod == 'virialRadius') Hot_Halo_Ram_Pressure_Stripping_Get => Hot_Halo_Ram_Pressure_Stripping_Virial_Radius
    return
  end subroutine Hot_Halo_Ram_Pressure_Stripping_Virial_Radii_Initialize

  double precision function Hot_Halo_Ram_Pressure_Stripping_Virial_Radius(thisNode)
    !% Computes the hot halo ram pressure stripping radius, assuming a null calculation in which that radius always equals the
    !% virial radius.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode  
                                                     
    Hot_Halo_Ram_Pressure_Stripping_Virial_Radius=Dark_Matter_Halo_Virial_Radius(thisNode)
    return
  end function Hot_Halo_Ram_Pressure_Stripping_Virial_Radius
  
end module Hot_Halo_Ram_Pressure_Stripping_Virial_Radii
