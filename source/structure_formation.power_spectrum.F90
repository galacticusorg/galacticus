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
Contains a module which implements a cosmological power spectra class.
!!}

module Power_Spectra
  !!{
  Implements linear-theory power spectra.
  !!}
  private

  !![
  <functionClass>
   <name>powerSpectrum</name>
   <descriptiveName>Linear Theory Power Spectrum</descriptiveName>
   <description>Class providing the linear theory matter power spectrum $P(k,t)$ as a function of wavenumber $k$ and cosmic
    time $t$. The power spectrum quantifies the amplitude of density fluctuations on different scales and is the fundamental
    statistical descriptor of the large-scale structure of the Universe. Both the dimensional form $P(k)$ and the dimensionless
    form $\Delta^2(k) = k^3 P(k) / 2\pi^2$ are available, along with the logarithmic slope $\mathrm{d}\ln P/\mathrm{d}\ln k$.</description>
   <default>standard</default>
   <method name="power" >
    <description>Return the linear power spectrum for $k=$\mono{wavenumber} [Mpc$^{-1}$].</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: wavenumber, time</argument>
   </method>
   <method name="powerDimensionless" >
    <description>Return the dimensionless linear power spectrum for $k=$\mono{wavenumber} [Mpc$^{-1}$].</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: wavenumber, time</argument>
   </method>
   <method name="powerLogarithmicDerivative" >
    <description>Return the logarithmic derivative of the power spectrum, $\mathrm{d}\ln P(k)/\mathrm{d}\ln k$, for $k=$\mono{wavenumber} [Mpc$^{-1}$].</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: wavenumber, time</argument>
   </method>
  </functionClass>
  !!]

end module Power_Spectra
