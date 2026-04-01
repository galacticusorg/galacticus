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
Contains a module which provides an object that implements cosmological parameters.
!!}

module Cosmology_Parameters
  !!{
  Provides an object that implements cosmological parameters.
  !!}
  private

  ! Enumeration for Hubble constant units.
  !![
  <enumeration>
   <name>hubbleUnits</name>
   <description>Enumeration of supported unit systems for the Hubble constant: standard cosmolgical units (km/s/Mpc), time units (inverse time), or the dimensionless little-$h$ convention where $H_0 = 100 h$ km/s/Mpc.</description>
   <visibility>public</visibility>
   <entry label="standard" />
   <entry label="time"     />
   <entry label="littleH"  />
  </enumeration>
  !!]

  !![
  <functionClass>
   <name>cosmologyParameters</name>
   <descriptiveName>Cosmological Parameters</descriptiveName>
   <description>Object providing the fundamental cosmological parameters that describe the background Universe: matter density
    $\Omega_\mathrm{m}$, dark energy density $\Omega_\Lambda$, baryon density $\Omega_\mathrm{b}$, radiation density
    $\Omega_\mathrm{r}$, curvature $\Omega_k$, Hubble constant $H_0$, and CMB temperature $T_\mathrm{CMB}$.
    These parameters are used throughout \glc\ to compute cosmological distances, times, and densities.</description>
   <default>simple</default>
   <method name="OmegaMatter" >
    <description>Return the cosmological matter density in units of the critical density at the present day.</description>
    <type>double precision</type>
    <pass>yes</pass>
   </method>
   <method name="OmegaDarkEnergy" >
    <description>Return the cosmological dark energy density in units of the critical density at the present day.</description>
    <type>double precision</type>
    <pass>yes</pass>
   </method>
   <method name="OmegaBaryon" >
    <description>Return the cosmological baryon density in units of the critical density at the present day.</description>
    <type>double precision</type>
    <pass>yes</pass>
   </method>
   <method name="OmegaRadiation" >
    <description>Return the cosmological radiation density in units of the critical density at the present day.</description>
    <type>double precision</type>
    <pass>yes</pass>
   </method>
   <method name="OmegaCurvature" >
    <description>Return the cosmological curvature density in units of the critical density at the present day.</description>
    <type>double precision</type>
    <pass>yes</pass>
   </method>
   <method name="HubbleConstant" >
    <description>Return the Hubble constant at the present day. The optional \mono{units} argument specifies if the return value should be in units of km/s/Mpc (hubbleUnitsStandard), Gyr$^{-1}$ (hubbleUnitsTime), or 100 km/s/Mpc (hubbleUnitsLittleH).</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(enumerationHubbleUnitsType), intent(in   ), optional :: units</argument>
   </method>
   <method name="temperatureCMB" >
    <description>Return the temperature of the cosmic microwave background radiation (in units of Kelvin) at the present day.</description>
    <type>double precision</type>
    <pass>yes</pass>
   </method>
   <method name="densityCritical" >
    <description>Return the critical density at the present day in units of $\mathrm{M}_\odot/$Mpc$^3$.</description>
    <type>double precision</type>
    <pass>yes</pass>
   </method>
  </functionClass>
  !!]

end module Cosmology_Parameters
