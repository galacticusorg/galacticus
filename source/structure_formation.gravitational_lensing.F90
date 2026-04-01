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
Contains a module which implements gravitational lensing from large scale structure.
!!}

module Gravitational_Lensing
  !!{
  Implements gravitational lensing from large scale structure.
  !!}
  private

  !![
  <functionClass>
  <name>gravitationalLensing</name>
   <descriptiveName>Gravitational Lensing</descriptiveName>
   <description>
    Class providing models of the gravitational lensing magnification distribution due to intervening
    large-scale structure along a line of sight---the probability density and cumulative distribution
    of the magnification factor $\mu$ as a function of source redshift and angular source size. Strong
    lensing by massive halos can boost observed fluxes significantly, affecting number counts and
    luminosity functions at the bright end. Implementations typically follow fitting functions calibrated
    to ray-tracing simulations (e.g.\ \citealt{takahashi_full-sky_2017}) and depend on the matter
    power spectrum and cosmological model.
   </description>
   <default>takahashi2011</default>
   <method name="magnificationPDF" >
    <description>Returns the differential probability $\mathrm{d}P/\mathrm{d}\mu$ for a source at the given \mono{redshift} and angular size \mono{scaleSource} to be magnified by factor \mono{magnification} due to gravitational lensing.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: magnification, redshift, scaleSource</argument>
   </method>
   <method name="magnificationCDF" >
    <description>Returns the cumulative probability $P(\mu' \le \mu)$ that a source at the given \mono{redshift} and angular size \mono{scaleSource} has a gravitational lensing magnification less than or equal to \mono{magnification}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: magnification, redshift, scaleSource</argument>
   </method>
  </functionClass>
  !!]

end module Gravitational_Lensing
