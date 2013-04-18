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

!% Contains a module which implements a calculation of the core radius in the hot halo density profile that is a fixed fraction of
!% the virial radius.

module Hot_Halo_Density_Cored_Isothermal_Core_Radii_Virial_Fraction
  !% Implements a calculation of the core radius in the hot halo density profile that is a fixed fraction of
  !% the virial radius.
  implicit none
  private
  public :: Hot_Halo_Density_Cored_Isothermal_Core_Radii_VF_Initialize

  ! Parameters of the model.
  double precision :: isothermalCoreRadiusOverVirialRadius
  
contains

  !# <hotHaloCoredIsothermalCoreRadiiMethod>
  !#  <unitName>Hot_Halo_Density_Cored_Isothermal_Core_Radii_VF_Initialize</unitName>
  !# </hotHaloCoredIsothermalCoreRadiiMethod>
  subroutine Hot_Halo_Density_Cored_Isothermal_Core_Radii_VF_Initialize(hotHaloCoredIsothermalCoreRadiiMethod&
       &,Hot_Halo_Density_Cored_Isothermal_Core_Radius_Get)
    !% Initializes the ``virial radius fraction'' cored isothermal hot halo profile core radius module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),                 intent(in)    :: hotHaloCoredIsothermalCoreRadiiMethod
    procedure(Hot_Halo_Density_Cored_Isothermal_Core_Radius_Virial_Fraction), pointer, intent(inout) :: Hot_Halo_Density_Cored_Isothermal_Core_Radius_Get
    
    if (hotHaloCoredIsothermalCoreRadiiMethod == 'virialRadiusFraction') then
       Hot_Halo_Density_Cored_Isothermal_Core_Radius_Get => Hot_Halo_Density_Cored_Isothermal_Core_Radius_Virial_Fraction
       !@ <inputParameter>
       !@   <name>isothermalCoreRadiusOverVirialRadius</name>
       !@   <defaultValue>0.4123</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The core radius in the ``cored isothermal'' hot halo density profile in units of the virial radius.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('isothermalCoreRadiusOverVirialRadius',isothermalCoreRadiusOverVirialRadius,defaultValue=0.4123d0)
    end if
    return
  end subroutine Hot_Halo_Density_Cored_Isothermal_Core_Radii_VF_Initialize

  double precision function Hot_Halo_Density_Cored_Isothermal_Core_Radius_Virial_Fraction(thisNode)
    !% Returns the radius (in Mpc) of the core radius in a cored isothermal hot halo density profile. Assumes that the radius is
    !% a fixed fraction of the halo virial radius.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode
    
    ! Compute the core radius.
    Hot_Halo_Density_Cored_Isothermal_Core_Radius_Virial_Fraction=isothermalCoreRadiusOverVirialRadius*Dark_Matter_Halo_Virial_Radius(thisNode)
   return
  end function Hot_Halo_Density_Cored_Isothermal_Core_Radius_Virial_Fraction
  
end module Hot_Halo_Density_Cored_Isothermal_Core_Radii_Virial_Fraction
