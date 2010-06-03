!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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






!% Contains a module which implements a isothermal halo spin distribution.

module Dark_Matter_Profiles_Isothermal
  !% Implements a isothermal halo spin distribution.
  use Tree_Nodes
  private
  public :: Dark_Matter_Profile_Isothermal_Initialize

contains

  !# <darkMatterProfileMethod>
  !#  <unitName>Dark_Matter_Profile_Isothermal_Initialize</unitName>
  !# </darkMatterProfileMethod>
  subroutine Dark_Matter_Profile_Isothermal_Initialize(darkMatterProfileMethod,Dark_Matter_Profile_Energy_Get&
       &,Dark_Matter_Profile_Energy_Growth_Rate_Get,Dark_Matter_Profile_Rotation_Normalization_Get &
       &,Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum_Get,Dark_Matter_Profile_Circular_Velocity_Get&
       &,Dark_Matter_Profile_Potential_Get,Dark_Matter_Profile_Enclosed_Mass_Get)
    !% Initializes the ``Isothermal'' halo spin distribution module.
    use ISO_Varying_String
    implicit none
    type(varying_string),          intent(in)    :: darkMatterProfileMethod
    procedure(),          pointer, intent(inout) :: Dark_Matter_Profile_Energy_Get,Dark_Matter_Profile_Energy_Growth_Rate_Get&
         &,Dark_Matter_Profile_Rotation_Normalization_Get,Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum_Get&
         &,Dark_Matter_Profile_Circular_Velocity_Get,Dark_Matter_Profile_Potential_Get,Dark_Matter_Profile_Enclosed_Mass_Get
    
    if (darkMatterProfileMethod == 'isothermal') then
       Dark_Matter_Profile_Energy_Get                                => Dark_Matter_Profile_Energy_Isothermal
       Dark_Matter_Profile_Energy_Growth_Rate_Get                    => Dark_Matter_Profile_Energy_Growth_Rate_Isothermal
       Dark_Matter_Profile_Rotation_Normalization_Get                => Dark_Matter_Profile_Rotation_Normalization_Isothermal
       Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum_Get => Radius_from_Specific_Angular_Momentum_Isothermal
       Dark_Matter_Profile_Circular_Velocity_Get                     => Dark_Matter_Profile_Circular_Velocity_Isothermal
       Dark_Matter_Profile_Potential_Get                             => Dark_Matter_Profile_Potential_Isothermal
       Dark_Matter_Profile_Enclosed_Mass_Get                         => Dark_Matter_Profile_Enclosed_Mass_Isothermal
    end if
    return
  end subroutine Dark_Matter_Profile_Isothermal_Initialize

  double precision function Dark_Matter_Profile_Enclosed_Mass_Isothermal(thisNode,radius)
    !% Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\tt thisNode} at the given {\tt radius} (given in
    !% units of Mpc).
    use Tree_Nodes
    use Tree_Node_Methods
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: radius

    Dark_Matter_Profile_Enclosed_Mass_Isothermal=Tree_Node_Mass(thisNode)*(radius/Dark_Matter_Halo_Virial_Radius(thisNode))
    return
  end function Dark_Matter_Profile_Enclosed_Mass_Isothermal
  
  double precision function Dark_Matter_Profile_Potential_Isothermal(thisNode,radius)
    !% Returns the potential (in (km/s)$^2$) in the dark matter profile of {\tt thisNode} at the given {\tt radius} (given in
    !% units of Mpc).
    use Tree_Nodes
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: radius

    Dark_Matter_Profile_Potential_Isothermal=(-1.0d0+dlog(radius/Dark_Matter_Halo_Virial_Radius(thisNode)))&
         &*Dark_Matter_Halo_Virial_Velocity(thisNode)**2
    return
  end function Dark_Matter_Profile_Potential_Isothermal
  
  double precision function Dark_Matter_Profile_Circular_Velocity_Isothermal(thisNode,radius)
    !% Returns the circular velocity (in km/s) in the dark matter profile of {\tt thisNode} at the given {\tt radius} (given in
    !% units of Mpc). For an isothermal halo this is independent of radius and therefore equal to the virial velocity.
    use Tree_Nodes
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: radius

    Dark_Matter_Profile_Circular_Velocity_Isothermal=Dark_Matter_Halo_Virial_Velocity(thisNode)
    return
  end function Dark_Matter_Profile_Circular_Velocity_Isothermal
  
  double precision function Radius_from_Specific_Angular_Momentum_Isothermal(thisNode,specificAngularMomentum)
    !% Returns the radius (in Mpc) in {\tt thisNode} at which a circular orbit has the given {\tt specificAngularMomentum} (given
    !% in units of km s$^{-1}$ Mpc). For an isothermal halo, the circular velocity is constant (and therefore equal to the virial
    !% velocity). Therefore, $r = j/V_{\rm virial}$ where $j$(={\tt specificAngularMomentum}) is the specific angular momentum and
    !% $r$ the required radius.
    use Tree_Nodes
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: specificAngularMomentum

    Radius_from_Specific_Angular_Momentum_Isothermal=specificAngularMomentum/Dark_Matter_Halo_Virial_Velocity(thisNode)
    return
  end function Radius_from_Specific_Angular_Momentum_Isothermal
  
  double precision function Dark_Matter_Profile_Rotation_Normalization_Isothermal(thisNode)
    !% Return the normalization of the rotation velocity vs. specific angular momentum relation.
    use Tree_Nodes
    use Tree_Node_Methods
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Math
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    Dark_Matter_Profile_Rotation_Normalization_Isothermal=2.0d0/Dark_Matter_Halo_Virial_Radius(thisNode)
    return
  end function Dark_Matter_Profile_Rotation_Normalization_Isothermal
  
  double precision function Dark_Matter_Profile_Energy_Isothermal(thisNode)
    !% Return the energy of an isothermal halo density profile.
    use Tree_Nodes
    use Tree_Node_Methods
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    Dark_Matter_Profile_Energy_Isothermal=-0.5d0*Tree_Node_Mass(thisNode)*Dark_Matter_Halo_Virial_Velocity(thisNode)**2
    return
  end function Dark_Matter_Profile_Energy_Isothermal
  
  double precision function Dark_Matter_Profile_Energy_Growth_Rate_Isothermal(thisNode)
    !% Return the rate of change of the energy of an isothermal halo density profile.
    use Tree_Nodes
    use Tree_Node_Methods
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    Dark_Matter_Profile_Energy_Growth_Rate_Isothermal=Dark_Matter_Profile_Energy_Isothermal(thisNode)&
         &*(Tree_Node_Mass_Accretion_Rate(thisNode)/Tree_Node_Mass(thisNode)+2.0d0&
         &*Dark_Matter_Halo_Virial_Velocity_Growth_Rate(thisNode)/Dark_Matter_Halo_Virial_Velocity(thisNode))
    return
  end function Dark_Matter_Profile_Energy_Growth_Rate_Isothermal
  
end module Dark_Matter_Profiles_Isothermal
