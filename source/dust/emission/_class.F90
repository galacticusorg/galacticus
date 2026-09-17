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
Contains a module which provides a class implementing spectra of thermal emission from interstellar dust.
!!}

module Dust_Emission_Spectra
  !!{RST
  Provides a class implementing spectra of thermal emission from interstellar dust.

  Dust re-radiates, in the infrared, the energy it absorbs from starlight and from nebular and accretion-disk emission.
  A ``dustEmissionSpectrum`` object is given the luminosity a population of dust absorbs, and its mass, and returns the
  spectrum of that re-emission. The absorbed luminosity is supplied by the ``absorbedFractions`` of a
  ``dustAttenuation`` object, and the mass by a ``dustProperties`` object, so that the dust which attenuates and the
  dust which emits are the same.

  Spectra are normalized analytically, not over the wavelengths at which they happen to be evaluated, so that the
  energy radiated is exactly what was absorbed---the property on which energy balance rests. The infrared emission is
  assumed not to be absorbed again.
  !!}
  implicit none
  private

  !![
  <functionClass docformat="rst">
   <name>dustEmissionSpectrum</name>
   <descriptiveName>Dust Emission Spectra</descriptiveName>
   <description>
   Class providing the spectrum of thermal emission from a population of interstellar dust.
   </description>
   <default>blackBodyModified</default>
   <method name="luminosity" >
    <description>
    Return the luminosity per unit frequency, :math:`L_\nu` in :math:`L_\odot\,\hbox{Hz}^{-1}`, emitted by dust at each
    of the given rest-frame ``wavelengths`` (in Å), given that the dust absorbs a luminosity ``luminosityAbsorbed`` (in
    :math:`L_\odot`), has a mass ``massDust`` (in :math:`M_\odot`), and is seen at cosmic ``time`` (in Gyr).

    The spectrum is normalized analytically, so that :math:`\int L_\nu\,\mathrm{d}\nu` over all frequencies equals the
    luminosity emitted, whatever wavelengths are requested. That is ``luminosityAbsorbed`` unless the dust is also
    heated by some other source---the cosmic microwave background, for example---in which case the implementation
    documents what it includes.
    </description>
    <type>double precision, dimension(size(wavelengths))</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ), dimension(:) :: wavelengths                       </argument>
    <argument>double precision, intent(in   )               :: luminosityAbsorbed, massDust, time</argument>
   </method>
  </functionClass>
  !!]

end module Dust_Emission_Spectra
