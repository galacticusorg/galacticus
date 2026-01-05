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
Contains a module that provides a class implementing dust attenuation of stellar spectra.
!!}

module Stellar_Spectra_Dust_Attenuations
  !!{
  Provides a class implementing dust attenuation of stellar spectra.
  !!}
  private

  !![
  <functionClass>
   <name>stellarSpectraDustAttenuation</name>
   <descriptiveName>Stellar Spectra Dust Attenuation</descriptiveName>
   <description>Class implementing dust attenuation of stellar spectra.</description>
   <default>zero</default>
   <method name="attenuation" >
    <description>Return the attenuation, in magnitudes, of stellar spectra due to dust at the given wavelength, age, and V-band extinction.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: wavelength, age, vBandAttenuation</argument>
   </method>
   <method name="isAgeDependent" >
    <description>Return true if the attenuation may depend on the age of the stellar population.</description>
    <type>logical</type>
    <pass>yes</pass>
    <code>
     !$GLC attributes unused :: self
     stellarSpectraDustAttenuationIsAgeDependent=.false.
    </code>
   </method>
   <method name="isSeparable" >
    <description>Return true if the attenuation is separable into a product of functions of wavelength, age, and V-band attenuation.</description>
    <type>logical</type>
    <pass>yes</pass>
    <code>
     !$GLC attributes unused :: self
     stellarSpectraDustAttenuationIsSeparable=.true.
    </code>
   </method>
  </functionClass>
  !!]

end module Stellar_Spectra_Dust_Attenuations
