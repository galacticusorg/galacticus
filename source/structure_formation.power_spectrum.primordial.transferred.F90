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

!% Contains a module which implements the primordial power spectrum transferred to late times.

module Primordial_Power_Spectra_Transferred
  !% Implements the primordial power spectrum transferred to late times.
  implicit none
  private
  public :: Primordial_Power_Spectrum_Transferred

contains

  double precision function Primordial_Power_Spectrum_Transferred(wavenumber)
    !% Return the primordial power spectrum transferred to late times for $k=${\tt wavenumber} [Mpc$^{-1}$].
    use Transfer_Functions
    use Primordial_Power_Spectra
    implicit none
    double precision, intent(in   ) :: wavenumber

    Primordial_Power_Spectrum_Transferred=(Transfer_Function(wavenumber)**2)*Primordial_Power_Spectrum(wavenumber)
    return
  end function Primordial_Power_Spectrum_Transferred

end module Primordial_Power_Spectra_Transferred
