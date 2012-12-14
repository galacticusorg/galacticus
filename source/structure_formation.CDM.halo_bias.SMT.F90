!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements calculations of halo bias using the fitting function of \cite{sheth_ellipsoidal_2001}.

module Dark_Matter_Halo_Biases_SMT
  !% Implements calculations of halo bias using the fitting function of \cite{sheth_ellipsoidal_2001}.
  implicit none
  private
  public :: Dark_Matter_Halo_Bias_SMT_Initialize

contains

  !# <darkMatterHaloBiasMethod>
  !#  <unitName>Dark_Matter_Halo_Bias_SMT_Initialize</unitName>
  !# </darkMatterHaloBiasMethod>
  subroutine Dark_Matter_Halo_Bias_SMT_Initialize(darkMatterHaloBiasMethod,Dark_Matter_Halo_Bias_Node_Get,Dark_Matter_Halo_Bias_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    implicit none
    type(varying_string),                 intent(in)    :: darkMatterHaloBiasMethod
    procedure(double precision), pointer, intent(inout) :: Dark_Matter_Halo_Bias_Get,Dark_Matter_Halo_Bias_Node_Get

    if (darkMatterHaloBiasMethod == 'SMT') then
       Dark_Matter_Halo_Bias_Node_Get => Dark_Matter_Halo_Bias_Node_SMT
       Dark_Matter_Halo_Bias_Get      => Dark_Matter_Halo_Bias_SMT
    end if
    return
  end subroutine Dark_Matter_Halo_Bias_SMT_Initialize

  double precision function Dark_Matter_Halo_Bias_Node_SMT(thisNode)
    !% Computes the bias for a dark matter halo using the method of \cite{mo_analytic_1996}.
    use Galacticus_Nodes
    implicit none
    type (treeNode          ), intent(inout), pointer :: thisNode
    class(nodeComponentBasic),                pointer :: thisBasicComponent

    ! Compute halo bias.
    thisBasicComponent => thisNode%basic()
    Dark_Matter_Halo_Bias_Node_SMT=Dark_Matter_Halo_Bias_SMT(thisBasicComponent%mass(),thisBasicComponent%time())
    return
  end function Dark_Matter_Halo_Bias_Node_SMT

  double precision function Dark_Matter_Halo_Bias_SMT(mass,time)
    !% Computes the bias for a dark matter halo using the method of \cite{sheth_ellipsoidal_2001}.
    use Critical_Overdensity
    use CDM_Power_Spectrum
    implicit none
    double precision, intent(in) :: mass,time
    double precision, parameter  :: a=0.707d0, b=0.5d0, c=0.6d0
    double precision             :: deltaCritical,sigma,nu

    ! Get critical overdensity for collapse and root-variance, then compute peak height parameter, nu.
    deltaCritical=Critical_Overdensity_for_Collapse(time=time,mass=mass)
    sigma        =sigma_CDM(mass)
    nu           =deltaCritical/sigma
    
    ! Compute halo bias.
    Dark_Matter_Halo_Bias_SMT=1.0d0+(dsqrt(a)*(a*nu**2)+dsqrt(a)*b*(a*nu**2)**(1.0d0-c)-(a*nu**2)**c/((a*nu**2)**c+b*(1.0d0-c)&
         &*(1.0d0-c/2.0d0)))/dsqrt(a)/deltaCritical
    return
  end function Dark_Matter_Halo_Bias_SMT

end module Dark_Matter_Halo_Biases_SMT
