!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Contains a module which provides a class that implements cosmological structure mass variance.

module Cosmological_Mass_Variance
  !% Provides an object that implements cosmological structure mass variance.
  use FGSL
  
  !# <functionClass>
  !#  <name>cosmologicalMassVariance</name>
  !#  <descriptiveName>Mass Variance of Cosmological Density Field</descriptiveName>
  !#  <description>Object providing mass variance of the cosmological density field.</description>
  !#  <default>filteredPower</default>
  !#  <stateful>yes</stateful>
  !#  <defaultThreadPrivate>yes</defaultThreadPrivate>
  !#  <method name="powerNormalization" >
  !#   <description>Return the normalization of the power spectrum.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#  </method>
  !#  <method name="sigma8" >
  !#   <description>Return the value fo $\sigma_8$.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#  </method>
  !#  <method name="rootVariance" >
  !#   <description>Return the root-variance of the cosmological density field.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: mass</argument>
  !#  </method>
  !#  <method name="rootVarianceLogarithmicGradient" >
  !#   <description>Return the logarithmic gradient of the root-variance of the cosmological density field with respect to mass.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: mass</argument>
  !#  </method>
  !#  <method name="rootVarianceAndLogarithmicGradient" >
  !#   <description>Return the value and logarithmic gradient with respect to mass of the root-variance of the cosmological density field.</description>
  !#   <type>void</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: mass</argument>
  !#   <argument>double precision, intent(  out) :: rootVariance, rootVarianceLogarithmicGradient</argument>
  !#  </method>
  !#  <method name="mass" >
  !#   <description>Return the mass corresponding to the given {\normalfont \ttfamily rootVariance} of the cosmological density field.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: rootVariance</argument>
  !#  </method>
  !# </functionClass>

end module Cosmological_Mass_Variance
