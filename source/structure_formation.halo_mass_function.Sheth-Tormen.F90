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

!% Contains a module which generates a tabulated Sheth-Tormen halo mass function.

module Halo_Mass_Function_Sheth_Tormen
  !% Implements generation of a tabulated power-law primordial power spectrum.
  implicit none
  private
  public :: Halo_Mass_Function_Sheth_Tormen_Initialize

  ! Parameters controlling the gridding of the power spectrum and default wavenumber range.
  integer, parameter :: nPointsPerDecade=1000

contains

  !# <haloMassFunctionMethod>
  !#  <unitName>Halo_Mass_Function_Sheth_Tormen_Initialize</unitName>
  !# </haloMassFunctionMethod>
  subroutine Halo_Mass_Function_Sheth_Tormen_Initialize(haloMassFunctionMethod,Halo_Mass_Function_Differential_Get)
    !% Initializes the ``Sheth-Tormen mass function'' module.
    use ISO_Varying_String
    implicit none
    type     (varying_string  ), intent(in   )          :: haloMassFunctionMethod
    procedure(Halo_Mass_Function_Sheth_Tormen_Differential), intent(inout), pointer :: Halo_Mass_Function_Differential_Get

    if (haloMassFunctionMethod == 'Sheth-Tormen') Halo_Mass_Function_Differential_Get => Halo_Mass_Function_Sheth_Tormen_Differential
    return
  end subroutine Halo_Mass_Function_Sheth_Tormen_Initialize

  double precision function Halo_Mass_Function_Sheth_Tormen_Differential(time,mass)
    !% Compute the Sheth-Tormen halo mass function.
    use Numerical_Constants_Math
    use Power_Spectra
    use Critical_Overdensity
    use Cosmology_Parameters
    implicit none
    double precision                          , intent(in   ) :: mass                           , time
    double precision                          , parameter     :: a                      =0.707d0, normalization=0.3221836349d0, &
         &                                                       p                      =0.3d0
    class           (cosmologyParametersClass), pointer       :: thisCosmologyParameters
    double precision                                          :: alpha                          , nu                          , &
         &                                                       nuPrime

    ! Get the default cosmology.
    thisCosmologyParameters => cosmologyParameters()
    ! Compute the mass function.
    nu=(Critical_Overdensity_for_Collapse(time=time,mass=mass)/Cosmological_Mass_Root_Variance(mass))**2
    nuPrime=a*nu
    alpha=abs(Cosmological_Mass_Root_Variance_Logarithmic_Derivative(mass))
    Halo_Mass_Function_Sheth_Tormen_Differential=(thisCosmologyParameters%OmegaMatter()*thisCosmologyParameters%densityCritical()/mass**2)*alpha*sqrt(2.0d0*nuPrime/Pi)*normalization&
         &*(1.0d0+1.0d0/nuPrime**p) *exp(-0.5d0*nuPrime)
    return
  end function Halo_Mass_Function_Sheth_Tormen_Differential

end module Halo_Mass_Function_Sheth_Tormen
