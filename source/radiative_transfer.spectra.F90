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

!!{
Contains a module which provides spectra for radiative transfer calculations.
!!}

module Radiative_Transfer_Spectra
  !!{
  Provides a class that implements spectra for radiative transfer calculations.
  !!}
  private

  !![
  <functionClass>
   <name>radiativeTransferSpectrum</name>
   <descriptiveName>Radiative Transfer Spectra</descriptiveName>
   <description>Class providing spectra for radiative transfer calculations.</description>
   <default>blackBody</default>
   <method name="luminosity" >
    <description>Return the luminosity in the given wavelength range.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: wavelengthMinimum, wavelengthMaximum</argument>
   </method>
   <method name="spectrum" >
    <description>Return the spectrum (in units of $L_\odot$ \AA$^{-1}$) of the source at the given {\normalfont \ttfamily wavelength}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: wavelength</argument>
   </method>
  </functionClass>
  !!]

end module Radiative_Transfer_Spectra
