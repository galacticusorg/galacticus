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

!% Contains a module which implements calculations of halo bias using the Press-Schechter mass function \citep{mo_analytic_1996}.

module Dark_Matter_Halo_Biases_Press_Schechter
  !% Implements calculations of halo bias using the Press-Schechter mass function \citep{mo_analytic_1996}.
  implicit none
  private
  public :: Dark_Matter_Halo_Bias_Press_Schechter_Initialize

contains

  !# <darkMatterHaloBiasMethod>
  !#  <unitName>Dark_Matter_Halo_Bias_Press_Schechter_Initialize</unitName>
  !# </darkMatterHaloBiasMethod>
  subroutine Dark_Matter_Halo_Bias_Press_Schechter_Initialize(darkMatterHaloBiasMethod,Dark_Matter_Halo_Bias_Node_Get,Dark_Matter_Halo_Bias_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    implicit none
    type(varying_string),                 intent(in)    :: darkMatterHaloBiasMethod
    procedure(double precision), pointer, intent(inout) :: Dark_Matter_Halo_Bias_Node_Get,Dark_Matter_Halo_Bias_Get

    if (darkMatterHaloBiasMethod == 'Press-Schechter') then
       Dark_Matter_Halo_Bias_Node_Get => Dark_Matter_Halo_Bias_Node_Press_Schechter
       Dark_Matter_Halo_Bias_Get      => Dark_Matter_Halo_Bias_Press_Schechter
    end if
    return
  end subroutine Dark_Matter_Halo_Bias_Press_Schechter_Initialize

  double precision function Dark_Matter_Halo_Bias_Node_Press_Schechter(thisNode)
    !% Computes the bias for a dark matter halo using the method of \cite{mo_analytic_1996}.
    use Critical_Overdensity
    use Power_Spectra
    use Galacticus_Nodes
    implicit none
    type (treeNode          ), intent(inout), pointer :: thisNode
    class(nodeComponentBasic),                pointer :: thisBasicComponent

    ! Compute halo bias.
    thisBasicComponent => thisNode%basic()
    Dark_Matter_Halo_Bias_Node_Press_Schechter=Dark_Matter_Halo_Bias_Press_Schechter(thisBasicComponent%mass()&
         &,thisBasicComponent%time())
    return
  end function Dark_Matter_Halo_Bias_Node_Press_Schechter

  double precision function Dark_Matter_Halo_Bias_Press_Schechter(mass,time)
    !% Computes the bias for a dark matter halo using the method of \cite{mo_analytic_1996}.
    use Critical_Overdensity
    use Power_Spectra
    implicit none
    double precision, intent(in) :: mass,time
    double precision             :: deltaCritical,sigma,nu

    ! Get critical overdensity for collapse and root-variance, then compute peak height parameter, nu.
    deltaCritical=Critical_Overdensity_for_Collapse(time=time,mass=mass)
    sigma        =Cosmological_Mass_Root_Variance(mass)
    nu           =deltaCritical/sigma
    
    ! Compute halo bias.
    Dark_Matter_Halo_Bias_Press_Schechter=1.0d0+(nu**2-1.0d0)/deltaCritical
    return
  end function Dark_Matter_Halo_Bias_Press_Schechter

end module Dark_Matter_Halo_Biases_Press_Schechter
