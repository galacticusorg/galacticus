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

!% Contains a module which generates a tabulated Press-Schechter halo mass function.

module Halo_Mass_Function_Press_Schechter
  !% Implements generation of a tabulated power-law primordial power spectrum.
  implicit none
  private
  public :: Halo_Mass_Function_Press_Schechter_Initialize

  ! Parameters controlling the gridding of the power spectrum and default wavenumber range.
  integer, parameter :: nPointsPerDecade=1000

contains

  !# <haloMassFunctionMethod>
  !#  <unitName>Halo_Mass_Function_Press_Schechter_Initialize</unitName>
  !# </haloMassFunctionMethod>
  subroutine Halo_Mass_Function_Press_Schechter_Initialize(haloMassFunctionMethod,Halo_Mass_Function_Differential_Get)
    !% Initializes the ``Press-Schechter mass functon'' module.
    use ISO_Varying_String
    implicit none
    type     (varying_string  ), intent(in   )          :: haloMassFunctionMethod
    procedure(double precision), intent(inout), pointer :: Halo_Mass_Function_Differential_Get

    if (haloMassFunctionMethod == 'Press-Schechter') Halo_Mass_Function_Differential_Get => Halo_Mass_Function_Differential_Press_Schechter
    return
  end subroutine Halo_Mass_Function_Press_Schechter_Initialize

 double precision function Halo_Mass_Function_Differential_Press_Schechter(time,mass)
    !% Compute the Press-Schechter halo mass function.
    use Power_Spectra
    use Cosmology_Parameters
    use Excursion_Sets_First_Crossings
    implicit none
    double precision                          , intent(in   ) :: mass                   , time
    class           (cosmologyParametersClass), pointer       :: thisCosmologyParameters
    double precision                                          :: alpha                  , variance

    ! Get the default cosmology.
    thisCosmologyParameters => cosmologyParameters()
    alpha   =abs(Cosmological_Mass_Root_Variance_Logarithmic_Derivative(mass))
    variance=Cosmological_Mass_Root_Variance(mass)**2
    Halo_Mass_Function_Differential_Press_Schechter=2.0d0*(thisCosmologyParameters%OmegaMatter()*thisCosmologyParameters%densityCritical()/mass**2)*alpha*variance&
         &*Excursion_Sets_First_Crossing_Probability(variance,time)
    return
  end function Halo_Mass_Function_Differential_Press_Schechter

end module Halo_Mass_Function_Press_Schechter
