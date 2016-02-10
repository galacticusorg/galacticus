!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements the cosmological power spectrum.

module Power_Spectra
  !% Implements the cosmological power spectrum.
  use Tables
  use ISO_Varying_String
  implicit none
  private
  public :: Power_Spectrum, Power_Spectrum_Logarithmic_Derivative, Power_Spectrum_Dimensionless
  
contains
  
  double precision function Power_Spectrum_Dimensionless(wavenumber)
    !% Return the dimensionless power spectrum, $\Delta^2(k)$, for $k=${\normalfont \ttfamily wavenumber} [Mpc$^{-1}$].
    use Numerical_Constants_Math
    implicit none
    double precision, intent(in   ) :: wavenumber

    Power_Spectrum_Dimensionless=+4.0d0                         &
         &                       *Pi                            &
         &                       *               wavenumber **3 &
         &                       *Power_Spectrum(wavenumber)    &
         &                       /(                             &
         &                         +2.0d0                       &
         &                         *Pi                          &
         &                        )**3
    return
  end function Power_Spectrum_Dimensionless

  double precision function Power_Spectrum_Logarithmic_Derivative(wavenumber)
    !% Return the logarithmic derivative of the power spectrum, ${\mathrm d}\ln P(k)/{\mathrm d}\ln k$, for $k=${\normalfont
    !% \ttfamily wavenumber} [Mpc$^{-1}$].
    use Power_Spectra_Primordial_Transferred
    implicit none
    double precision                                         , intent(in   ) :: wavenumber
    class           (powerSpectrumPrimordialTransferredClass), pointer       :: powerSpectrumPrimordialTransferred_

    powerSpectrumPrimordialTransferred_ => powerSpectrumPrimordialTransferred()
    Power_Spectrum_Logarithmic_Derivative=powerSpectrumPrimordialTransferred_%logarithmicDerivative(wavenumber)
    return
  end function Power_Spectrum_Logarithmic_Derivative

  double precision function Power_Spectrum(wavenumber)
    !% Return the cosmological power spectrum for $k=${\normalfont \ttfamily wavenumber} [Mpc$^{-1}$].
    use Cosmological_Mass_Variance
    use Power_Spectra_Primordial_Transferred
    implicit none
    double precision                                         , intent(in   ) :: wavenumber
    class           (cosmologicalMassVarianceClass          ), pointer       :: cosmologicalMassVariance_
    class           (powerSpectrumPrimordialTransferredClass), pointer       :: powerSpectrumPrimordialTransferred_

    ! Get required objects.
    powerSpectrumPrimordialTransferred_ => powerSpectrumPrimordialTransferred()
    cosmologicalMassVariance_           => cosmologicalMassVariance          ()
    ! Compute the power spectrum.    
    Power_Spectrum=powerSpectrumPrimordialTransferred_%power(wavenumber)*cosmologicalMassVariance_%powerNormalization()
    return
  end function Power_Spectrum

end module Power_Spectra
