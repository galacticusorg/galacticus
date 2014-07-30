!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements calculations of halo bias using the fitting function of \cite{tinker_large_2010}.

module Dark_Matter_Halo_Biases_Tinker2010
  !% Implements calculations of halo bias using the fitting function of \cite{tinker_large_2010}.
  implicit none
  private
  public :: Dark_Matter_Halo_Bias_Tinker2010_Initialize

contains

  !# <darkMatterHaloBiasMethod>
  !#  <unitName>Dark_Matter_Halo_Bias_Tinker2010_Initialize</unitName>
  !# </darkMatterHaloBiasMethod>
  subroutine Dark_Matter_Halo_Bias_Tinker2010_Initialize(darkMatterHaloBiasMethod,Dark_Matter_Halo_Bias_Node_Get,Dark_Matter_Halo_Bias_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    implicit none
    type     (varying_string                       ), intent(in   )          :: darkMatterHaloBiasMethod
    procedure(Dark_Matter_Halo_Bias_Tinker2010     ), intent(inout), pointer :: Dark_Matter_Halo_Bias_Get
    procedure(Dark_Matter_Halo_Bias_Node_Tinker2010), intent(inout), pointer :: Dark_Matter_Halo_Bias_Node_Get

    if (darkMatterHaloBiasMethod == 'Tinker2010') then
       Dark_Matter_Halo_Bias_Node_Get => Dark_Matter_Halo_Bias_Node_Tinker2010
       Dark_Matter_Halo_Bias_Get      => Dark_Matter_Halo_Bias_Tinker2010
    end if
    return
  end subroutine Dark_Matter_Halo_Bias_Tinker2010_Initialize

  double precision function Dark_Matter_Halo_Bias_Node_Tinker2010(thisNode)
    !% Computes the bias for a dark matter halo using the method of \cite{mo_analytic_1996}.
    use Galacticus_Nodes
    implicit none
    type (treeNode          ), intent(inout), pointer :: thisNode
    class(nodeComponentBasic)               , pointer :: thisBasicComponent

    ! Compute halo bias.
    thisBasicComponent => thisNode%basic()
    Dark_Matter_Halo_Bias_Node_Tinker2010=Dark_Matter_Halo_Bias_Tinker2010(thisBasicComponent%mass(),thisBasicComponent%time())
    return
  end function Dark_Matter_Halo_Bias_Node_Tinker2010

  double precision function Dark_Matter_Halo_Bias_Tinker2010(mass,time)
    !% Computes the bias for a dark matter halo using the method of \cite{tinker_large_2010}.
    use Critical_Overdensity
    use Power_Spectra
    use Virial_Density_Contrast
    implicit none
    double precision                            , intent(in   ) :: mass                        , time
    double precision                            , parameter     :: lowerB                =1.5d0, lowerC             =2.4d0 , upperB=0.183d0
    class           (virialDensityContrastClass), pointer       :: virialDensityContrast_
    double precision                                            :: deltaCritical               , haloDensityContrast       , nu            , &
         &                                                         sigma                       , y
    double precision, save                                      :: lowerA                      , timePrevious       =-1.0d0, upperA        , &
         &                                                         upperC
    !$omp threadprivate(timePrevious,lowerA,upperA,upperC)

    ! Get critical overdensity for collapse and root-variance, then compute peak height parameter, nu.
    deltaCritical=Critical_Overdensity_for_Collapse(time=time,mass=mass)
    sigma        =Cosmological_Mass_Root_Variance(mass)
    nu           =deltaCritical/sigma

    ! Update fitting parameters if the time has changed.
    if (time /= timePrevious) then
       ! Store the new time.
       timePrevious=time

       ! Compute halo density contrast and logarithm.
       virialDensityContrast_    => virialDensityContrast()
       haloDensityContrast=virialDensityContrast_%densityContrast(time)
       y=log10(haloDensityContrast)

       ! Compute parameters as a function of halo overdensity (from Table 2 of Tinker et al. 2010)
       upperA=1.0d0+0.24d0*y*exp(-(4.0d0/y)**4)
       lowerA=0.44d0*y-0.88d0
       upperC=0.019d0+0.107d0*y+0.19d0*exp(-(4.0d0/y)**4)
    end if

    ! Compute halo bias (equation 6 of Tinker et al. 2010).
    Dark_Matter_Halo_Bias_Tinker2010=1.0d0-upperA*nu**lowerA/(nu**lowerA+deltaCritical**lowerA)+upperB*nu**lowerB+upperC*nu&
         &**lowerC
    return
  end function Dark_Matter_Halo_Bias_Tinker2010

end module Dark_Matter_Halo_Biases_Tinker2010
