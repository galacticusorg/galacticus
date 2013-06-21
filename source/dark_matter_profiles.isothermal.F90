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

!% Contains a module which implements a isothermal halo spin distribution.

module Dark_Matter_Profiles_Isothermal
  !% Implements a isothermal halo spin distribution.
  implicit none
  private
  public :: Dark_Matter_Profile_Isothermal_Initialize

contains

  !# <darkMatterProfileMethod>
  !#  <unitName>Dark_Matter_Profile_Isothermal_Initialize</unitName>
  !# </darkMatterProfileMethod>
  subroutine Dark_Matter_Profile_Isothermal_Initialize(darkMatterProfileMethod,Dark_Matter_Profile_Density_Get &
       &,Dark_Matter_Profile_Energy_Get ,Dark_Matter_Profile_Energy_Growth_Rate_Get &
       &,Dark_Matter_Profile_Rotation_Normalization_Get ,Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum_Get &
       &,Dark_Matter_Profile_Circular_Velocity_Get ,Dark_Matter_Profile_Potential_Get,Dark_Matter_Profile_Enclosed_Mass_Get &
       &,Dark_Matter_Profile_kSpace_Get,Dark_Matter_Profile_Freefall_Radius_Get&
       &,Dark_Matter_Profile_Freefall_Radius_Increase_Rate_Get)
    !% Initializes the ``Isothermal'' halo spin distribution module.
    use ISO_Varying_String
    implicit none
    type     (varying_string                                              ), intent(in   )          :: darkMatterProfileMethod                                        
    procedure(Dark_Matter_Profile_Density_Isothermal                      ), intent(inout), pointer :: Dark_Matter_Profile_Density_Get                                
    procedure(Dark_Matter_Profile_Energy_Isothermal                       ), intent(inout), pointer :: Dark_Matter_Profile_Energy_Get                                 
    procedure(Dark_Matter_Profile_Energy_Growth_Rate_Isothermal           ), intent(inout), pointer :: Dark_Matter_Profile_Energy_Growth_Rate_Get                     
    procedure(Dark_Matter_Profile_Rotation_Normalization_Isothermal       ), intent(inout), pointer :: Dark_Matter_Profile_Rotation_Normalization_Get                 
    procedure(Radius_from_Specific_Angular_Momentum_Isothermal            ), intent(inout), pointer :: Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum_Get  
    procedure(Dark_Matter_Profile_Circular_Velocity_Isothermal            ), intent(inout), pointer :: Dark_Matter_Profile_Circular_Velocity_Get                      
    procedure(Dark_Matter_Profile_Potential_Isothermal                    ), intent(inout), pointer :: Dark_Matter_Profile_Potential_Get                              
    procedure(Dark_Matter_Profile_Enclosed_Mass_Isothermal                ), intent(inout), pointer :: Dark_Matter_Profile_Enclosed_Mass_Get                          
    procedure(Dark_Matter_Profile_kSpace_Isothermal                       ), intent(inout), pointer :: Dark_Matter_Profile_kSpace_Get                                 
    procedure(Dark_Matter_Profile_Freefall_Radius_Isothermal              ), intent(inout), pointer :: Dark_Matter_Profile_Freefall_Radius_Get                        
    procedure(Dark_Matter_Profile_Freefall_Radius_Increase_Rate_Isothermal), intent(inout), pointer :: Dark_Matter_Profile_Freefall_Radius_Increase_Rate_Get          
                                                                                                                                                                   
    if (darkMatterProfileMethod == 'isothermal') then
       Dark_Matter_Profile_Density_Get                               => Dark_Matter_Profile_Density_Isothermal
       Dark_Matter_Profile_Energy_Get                                => Dark_Matter_Profile_Energy_Isothermal
       Dark_Matter_Profile_Energy_Growth_Rate_Get                    => Dark_Matter_Profile_Energy_Growth_Rate_Isothermal
       Dark_Matter_Profile_Rotation_Normalization_Get                => Dark_Matter_Profile_Rotation_Normalization_Isothermal
       Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum_Get => Radius_from_Specific_Angular_Momentum_Isothermal
       Dark_Matter_Profile_Circular_Velocity_Get                     => Dark_Matter_Profile_Circular_Velocity_Isothermal
       Dark_Matter_Profile_Potential_Get                             => Dark_Matter_Profile_Potential_Isothermal
       Dark_Matter_Profile_Enclosed_Mass_Get                         => Dark_Matter_Profile_Enclosed_Mass_Isothermal
       Dark_Matter_Profile_kSpace_Get                                => Dark_Matter_Profile_kSpace_Isothermal
       Dark_Matter_Profile_Freefall_Radius_Get                       => Dark_Matter_Profile_Freefall_Radius_Isothermal
       Dark_Matter_Profile_Freefall_Radius_Increase_Rate_Get         => Dark_Matter_Profile_Freefall_Radius_Increase_Rate_Isothermal
    end if
    return
  end subroutine Dark_Matter_Profile_Isothermal_Initialize

  double precision function Dark_Matter_Profile_Density_Isothermal(thisNode,radius)
    !% Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\tt thisNode} at the given {\tt radius} (given
    !% in units of Mpc).
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Math
    implicit none
    type            (treeNode          ), intent(inout), pointer :: thisNode            
    double precision                    , intent(in   )          :: radius              
    class           (nodeComponentBasic)               , pointer :: thisBasicComponent  
                                                                                     
    thisBasicComponent => thisNode%basic()
    Dark_Matter_Profile_Density_Isothermal=thisBasicComponent%mass()/4.0d0/Pi/Dark_Matter_Halo_Virial_Radius(thisNode)/radius**2
    return
  end function Dark_Matter_Profile_Density_Isothermal
  
  double precision function Dark_Matter_Profile_Enclosed_Mass_Isothermal(thisNode,radius)
    !% Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\tt thisNode} at the given {\tt radius} (given in
    !% units of Mpc).
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    implicit none
    type            (treeNode          ), intent(inout), pointer :: thisNode            
    double precision                    , intent(in   )          :: radius              
    class           (nodeComponentBasic)               , pointer :: thisBasicComponent  
                                                                                     
    thisBasicComponent => thisNode%basic()
    Dark_Matter_Profile_Enclosed_Mass_Isothermal=thisBasicComponent%mass()*(radius/Dark_Matter_Halo_Virial_Radius(thisNode))
    return
  end function Dark_Matter_Profile_Enclosed_Mass_Isothermal
  
  double precision function Dark_Matter_Profile_Potential_Isothermal(thisNode,radius)
    !% Returns the potential (in (km/s)$^2$) in the dark matter profile of {\tt thisNode} at the given {\tt radius} (given in
    !% units of Mpc).
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    use Galacticus_Error
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode  
    double precision          , intent(in   )          :: radius    
                                                                 
    if (radius <= 0.0d0) call Galacticus_Error_Report('Dark_Matter_Profile_Potential_Isothermal','isothermal profile potential is infinite at zero radius')
    Dark_Matter_Profile_Potential_Isothermal=(-1.0d0+log(radius/Dark_Matter_Halo_Virial_Radius(thisNode)))&
         &*Dark_Matter_Halo_Virial_Velocity(thisNode)**2
    return
  end function Dark_Matter_Profile_Potential_Isothermal
  
  double precision function Dark_Matter_Profile_Circular_Velocity_Isothermal(thisNode,radius)
    !% Returns the circular velocity (in km/s) in the dark matter profile of {\tt thisNode} at the given {\tt radius} (given in
    !% units of Mpc). For an isothermal halo this is independent of radius and therefore equal to the virial velocity.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode  
    double precision          , intent(in   )          :: radius    
                                                                 
    Dark_Matter_Profile_Circular_Velocity_Isothermal=Dark_Matter_Halo_Virial_Velocity(thisNode)
    return
  end function Dark_Matter_Profile_Circular_Velocity_Isothermal
  
  double precision function Radius_from_Specific_Angular_Momentum_Isothermal(thisNode,specificAngularMomentum)
    !% Returns the radius (in Mpc) in {\tt thisNode} at which a circular orbit has the given {\tt specificAngularMomentum} (given
    !% in units of km s$^{-1}$ Mpc). For an isothermal halo, the circular velocity is constant (and therefore equal to the virial
    !% velocity). Therefore, $r = j/V_{\rm virial}$ where $j$(={\tt specificAngularMomentum}) is the specific angular momentum and
    !% $r$ the required radius.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode                 
    double precision          , intent(in   )          :: specificAngularMomentum  
                                                                                
    Radius_from_Specific_Angular_Momentum_Isothermal=specificAngularMomentum/Dark_Matter_Halo_Virial_Velocity(thisNode)
    return
  end function Radius_from_Specific_Angular_Momentum_Isothermal
  
  double precision function Dark_Matter_Profile_Rotation_Normalization_Isothermal(thisNode)
    !% Return the normalization of the rotation velocity vs. specific angular momentum relation.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode  
                                                     
    Dark_Matter_Profile_Rotation_Normalization_Isothermal=2.0d0/Dark_Matter_Halo_Virial_Radius(thisNode)
    return
  end function Dark_Matter_Profile_Rotation_Normalization_Isothermal
  
  double precision function Dark_Matter_Profile_Energy_Isothermal(thisNode)
    !% Return the energy of an isothermal halo density profile.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    implicit none
    type (treeNode          ), intent(inout), pointer :: thisNode            
    class(nodeComponentBasic)               , pointer :: thisBasicComponent  
                                                                          
    thisBasicComponent => thisNode%basic()
    Dark_Matter_Profile_Energy_Isothermal=-0.5d0*thisBasicComponent%mass()*Dark_Matter_Halo_Virial_Velocity(thisNode)**2
    return
  end function Dark_Matter_Profile_Energy_Isothermal
  
  double precision function Dark_Matter_Profile_Energy_Growth_Rate_Isothermal(thisNode)
    !% Return the rate of change of the energy of an isothermal halo density profile.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    implicit none
    type (treeNode          ), intent(inout), pointer :: thisNode            
    class(nodeComponentBasic)               , pointer :: thisBasicComponent  
                                                                          
    thisBasicComponent => thisNode%basic()
    Dark_Matter_Profile_Energy_Growth_Rate_Isothermal=Dark_Matter_Profile_Energy_Isothermal(thisNode)&
         &*(thisBasicComponent%accretionRate()/thisBasicComponent%mass()+2.0d0&
         &*Dark_Matter_Halo_Virial_Velocity_Growth_Rate(thisNode)/Dark_Matter_Halo_Virial_Velocity(thisNode))
    return
  end function Dark_Matter_Profile_Energy_Growth_Rate_Isothermal

  double precision function Dark_Matter_Profile_kSpace_Isothermal(thisNode,waveNumber)
    !% Returns the Fourier transform of the isothermal density profile at the specified {\tt waveNumber} (given in Mpc$^{-1}$), using the
    !% expression given in \citeauthor{cooray_halo_2002}~(\citeyear{cooray_halo_2002}; table~1).
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    use Exponential_Integrals
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode                          
    double precision          , intent(in   )          :: waveNumber                        
    double precision                                   :: radiusScale, waveNumberScaleFree  
    
    ! Get the scale radius (for which we use the virial radius).                                                                                     
    radiusScale=Dark_Matter_Halo_Virial_Radius(thisNode)

    ! Get the dimensionless wavenumber.
    waveNumberScaleFree=waveNumber*radiusScale

    ! Compute the Fourier transformed profile.
    Dark_Matter_Profile_kSpace_Isothermal=Sine_Integral(waveNumberScaleFree)/waveNumberScaleFree

    return
  end function Dark_Matter_Profile_kSpace_Isothermal
  
  double precision function Dark_Matter_Profile_Freefall_Radius_Isothermal(thisNode,time)
    !% Returns the freefall radius in the isothermal density profile at the specified {\tt time} (given in Gyr). For an isothermal
    !% potential, the freefall radius, $r_{\rm ff}(t)$, is:
    !% \begin{equation}
    !% r_{\rm ff}(t) = \sqrt{{2 \over \pi}} V_{\rm virial} t.
    !% \end{equation}
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Astronomical
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode  
    double precision          , intent(in   )          :: time      
                                                                 
    Dark_Matter_Profile_Freefall_Radius_Isothermal=sqrt(2.0d0/Pi)*Dark_Matter_Halo_Virial_Velocity(thisNode)*time&
         &/Mpc_per_km_per_s_To_Gyr
    return
  end function Dark_Matter_Profile_Freefall_Radius_Isothermal
  
  double precision function Dark_Matter_Profile_Freefall_Radius_Increase_Rate_Isothermal(thisNode,time)
    !% Returns the rate of increase of the freefall radius in the isothermal density profile at the specified {\tt time} (given in
    !% Gyr). For an isothermal potential, the rate of increase of the freefall radius, $\dot{r}_{\rm ff}(t)$, is:
    !% \begin{equation}
    !% \dot{r}_{\rm ff}(t) = \sqrt{{2 \over \pi}} V_{\rm virial}.
    !% \end{equation}
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Astronomical
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode  
    double precision          , intent(in   )          :: time      
                                                                 
    Dark_Matter_Profile_Freefall_Radius_Increase_Rate_Isothermal=sqrt(2.0d0/Pi)*Dark_Matter_Halo_Virial_Velocity(thisNode)&
         &/Mpc_per_km_per_s_To_Gyr
    return
  end function Dark_Matter_Profile_Freefall_Radius_Increase_Rate_Isothermal
  
end module Dark_Matter_Profiles_Isothermal
