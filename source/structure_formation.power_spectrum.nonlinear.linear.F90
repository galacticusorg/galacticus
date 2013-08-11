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

!% Contains a module which implements the nonlinear power spectrum as the linear power spectrum (useful mostly for testing).

module Power_Spectra_Nonlinear_Linear
  !% Implements the nonlinear power spectrum as the linear power spectrum (useful mostly for testing).
  implicit none
  private
  public :: Power_Spectrum_Nonlinear_Linear_Initialize

contains

  !# <powerSpectrumNonlinearMethod>
  !#  <unitName>Power_Spectrum_Nonlinear_Linear_Initialize</unitName>
  !# </powerSpectrumNonlinearMethod>
  subroutine Power_Spectrum_Nonlinear_Linear_Initialize(powerSpectrumNonlinearMethod,Power_Spectrum_Nonlinear_Get)
    !% Initializes the ``lienar'' nonlinear power spectrum module.
    use ISO_Varying_String
    implicit none
    type     (varying_string  ), intent(in   )          :: powerSpectrumNonlinearMethod
    procedure(double precision), intent(inout), pointer :: Power_Spectrum_Nonlinear_Get

    if (powerSpectrumNonlinearMethod == 'linear') Power_Spectrum_Nonlinear_Get => Power_Spectrum_Nonlinear_Linear
    return
  end subroutine Power_Spectrum_Nonlinear_Linear_Initialize

  double precision function Power_Spectrum_Nonlinear_Linear(waveNumber,time)
    !% Return a nonlinear power spectrum equal to the linear power spectrum. (Useful mostly for testing.)
    use Power_Spectra
    use Linear_Growth
    implicit none
    double precision, intent(in   ) :: time, waveNumber

    Power_Spectrum_Nonlinear_Linear=Power_Spectrum(wavenumber)*Linear_Growth_Factor(time)**2
    return
  end function Power_Spectrum_Nonlinear_Linear

end module Power_Spectra_Nonlinear_Linear
