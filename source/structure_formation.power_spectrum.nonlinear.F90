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
Contains a module which implements the nonlinear power spectrum.
!!}

module Power_Spectra_Nonlinear
  !!{
  Implements the nonlinear power spectrum.
  !!}
  private

  !![
  <functionClass>
   <name>powerSpectrumNonlinear</name>
   <descriptiveName>Nonlinear Power Spectrum</descriptiveName>
   <description>Class providing non-linear power spectra.</description>
   <default>cosmicEmu</default>
   <method name="value" >
    <description>Return the nonlinear power spectrum for $k=${\normalfont \ttfamily wavenumber} [Mpc$^{-1}$] at cosmic time $t=${\normalfont \ttfamily time} [Gyr].</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: wavenumber, time</argument>
   </method>
  </functionClass>
  !!]

end module Power_Spectra_Nonlinear
