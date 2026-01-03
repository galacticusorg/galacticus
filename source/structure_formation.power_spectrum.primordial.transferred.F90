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
Contains a module which provides a class that implements the transferred primordial power spectrum.
!!}

module Power_Spectra_Primordial_Transferred
  !!{
  Provides a class that implements the transferred primordial power spectrum.
  !!}
  private

  !![
  <functionClass>
   <name>powerSpectrumPrimordialTransferred</name>
   <descriptiveName>Transferred Primordial Power Spectrum</descriptiveName>
   <description>Class providing the transferred primordial power spectrum.</description>
   <default>simple</default>
   <method name="power" >
    <description>Return the (unnormalized) power in the transferred primordial power spectrum at the given {\normalfont \ttfamily wavenumber} (specified in units of Mpc$^{-1}$.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: wavenumber, time</argument>
   </method>
   <method name="logarithmicDerivative" >
    <description>Return the logarithmic derivative with respect to wavenumber of the transferred primordial power spectrum at the given {\normalfont \ttfamily wavenumber} (specified in units of Mpc$^{-1}$.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: wavenumber, time</argument>
   </method>
   <method name="growthIsWavenumberDependent" >
    <description>Return true if the growth of the power spectrum is wavenumber-dependent.</description>
    <type>logical</type>
    <pass>yes</pass>
   </method>
  </functionClass>
  !!]

end module Power_Spectra_Primordial_Transferred
