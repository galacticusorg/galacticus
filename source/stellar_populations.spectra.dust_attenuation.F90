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

!% Contains a module that provides a class implementing dust attenuation of stellar spectra.

module Stellar_Spectra_Dust_Attenuations
  !% Provides a class implementing dust attenuation of stellar spectra.
  use ISO_Varying_String
  !# <include directive="stellarSpectraDustAttenuation" type="functionModules" >
  include 'stellarSpectraDustAttenuation.functionModules.inc'
  !# </include>
  private

  !# <include directive="stellarSpectraDustAttenuation" type="function" >
  !#  <descriptiveName>Stellar Spectra Dust Attenuation</descriptiveName>
  !#  <description>Class implementing dust attenuation of stellar spectra..</description>
  !#  <default>null</default>
  !#  <defaultThreadPrivate>yes</defaultThreadPrivate>
  !#  <method name="attenuation" >
  !#   <description>Return the attenuation, in magnitudes, of stellar spectra due to dust at the given wavelength, age, and V-band extinction.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: wavelength, age, vBandAttenuation</argument>
  !#  </method>
  include 'stellarSpectraDustAttenuation.type.inc'
  !# </include>

end module Stellar_Spectra_Dust_Attenuations
