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
Contains a module which provides a class that implements window functions for computing the
variance of the power spectrum.
!!}

module Power_Spectrum_Window_Functions
  !!{
  Provides a class which implements window functions for computing the variance of the power
  spectrum.
  !!}
  private

  !![
  <functionClass>
   <name>powerSpectrumWindowFunction</name>
   <descriptiveName>Power Spectrum Window Functions</descriptiveName>
   <description>Class providing window functions for filtering of power spectra.</description>
   <default>topHat</default>
   <defaultThreadprivate>no</defaultThreadprivate>
   <method name="value" >
    <description> Returns the window function for power spectrum variance computation at the specified {\normalfont \ttfamily wavenumber} (in Mpc$^{-1}$) for a given {\normalfont \ttfamily smoothingMass} (in $M_\odot$).</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: wavenumber, smoothingMass</argument>
   </method>
   <method name="wavenumberMaximum" >
    <description>Returns the maximum wavenumber for which the window function for power spectrum variance computation is non-zero for a given {\normalfont \ttfamily smoothingMass} (in $M_\odot$).</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: smoothingMass</argument>
   </method>
   <method name="amplitudeIsMassIndependent" >
    <description>Should return true if, and only if, the amplitude of the window function below the maximum wavenumber is independent of the smoothing mass scale.</description>
    <type>logical</type>
    <pass>yes</pass>
    <code>
      !$GLC attributes unused :: self
      powerSpectrumWindowFunctionAmplitudeIsMassIndependent=.false.
    </code>
   </method>
  </functionClass>
  !!]

end module Power_Spectrum_Window_Functions
