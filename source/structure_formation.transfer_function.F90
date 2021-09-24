!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
Contains a module which provides a class that implements transfer functions.
!!}

module Transfer_Functions
  !!{
  Provides an object that implements transfer functions.
  !!}
  private

  !![
  <functionClass>
   <name>transferFunction</name>
   <descriptiveName>Transfer Function</descriptiveName>
   <description>Class providing transfer functions for power spectra.</description>
   <default>eisensteinHu1999</default>
   <method name="value" >
    <description>Return the transfer function for $k=${\normalfont \ttfamily wavenumber} [Mpc$^{-1}$].</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: wavenumber</argument>
   </method>
   <method name="logarithmicDerivative" >
    <description>Return the logarithmic derivative of the transfer function for $k=${\normalfont \ttfamily wavenumber} [Mpc$^{-1}$].</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: wavenumber</argument>
   </method>
   <method name="epochTime" >
    <description>Return the cosmic time corresponding to the epoch for which this transfer function is defined.</description>
    <type>double precision</type>
    <pass>yes</pass>
   </method>
   <method name="halfModeMass" >
    <description>Return the mass (in $M_\odot$) corresponding to the wavenumber at which the transfer function is suppressed by a factor of two due to small-scale dark matter particle physics.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>integer, intent(  out), optional :: status</argument>
   </method>
   <method name="quarterModeMass" >
    <description>Return the mass (in $M_\odot$) corresponding to the wavenumber at which the transfer function is suppressed by a factor of four due to small-scale dark matter particle physics.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>integer, intent(  out), optional :: status</argument>
   </method>
   <method name="fractionModeMass" >
    <description>Return the mass (in $M_\odot$) corresponding to the wavenumber at which the transfer function is suppressed is reduced by {\normalfont \ttfamily fraction} due to small-scale dark matter particle physics.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>double precision, intent(in   )           :: fraction</argument>
    <argument>integer         , intent(  out), optional :: status</argument>
   </method>
  </functionClass>
  !!]

end module Transfer_Functions
