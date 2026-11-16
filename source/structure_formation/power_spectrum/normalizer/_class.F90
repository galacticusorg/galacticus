!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
!!    Andrew Benson <abenson@carnegiescience.edu>
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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a module which provides a class that implements normalization of power spectra.
!!}

module Power_Spectrum_Normalizers
  !!{RST
  Provides a class that implements normalization of power spectra.
  !!}
  use :: Power_Spectrum_Normalizer_Interfaces, only : powerSpectrumNormalizerPowerUnnormalized
  private
  public :: powerSpectrumNormalizerPowerUnnormalized

  !![
  <functionClass docformat="rst">
   <name>powerSpectrumNormalizer</name>
   <descriptiveName>Power Spectrum Normalizers</descriptiveName>
   <description>
   Class providing the normalization of a transferred primordial power spectrum (see
   :galacticus-class:`powerSpectrumPrimordialTransferred`). Given the unnormalized power spectrum, the normalizer returns the
   factor by which it must be multiplied to give the normalized power spectrum, for example such that the linear theory
   root-variance of the density field at the present day in spheres of radius :math:`8h^{-1}` Mpc equals :math:`\sigma_8`.

   The normalizer is a member of the transferred primordial power spectrum, which it normalizes, so every quantity derived from
   that power spectrum (the power spectrum itself, :math:`\sigma(M)`, halo mass functions, etc.) is consistently normalized.
   </description>
   <default>sigma8</default>
   <method name="normalization" >
    <description>
    Return the factor by which the power spectrum given by ``powerUnnormalized`` must be multiplied to normalize it. The function
    ``powerUnnormalized(wavenumber,time)`` returns the unnormalized transferred power spectrum at the given ``wavenumber``
    (in units of Mpc\ :math:`^{-1}`) and ``time`` (in units of Gyr).
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>procedure(powerSpectrumNormalizerPowerUnnormalized) :: powerUnnormalized</argument>
   </method>
  </functionClass>
  !!]

end module Power_Spectrum_Normalizers
