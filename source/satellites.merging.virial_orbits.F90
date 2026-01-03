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
Contains a module which provides a class implementing satellite orbital parameters at virial radius crossing.
!!}

module Virial_Orbits
  !!{
  Provides a class implementing satellite orbital parameters at virial radius crossing.
  !!}
  use :: Galacticus_Nodes       , only : treeNode
  use :: Kepler_Orbits          , only : keplerOrbit
  use :: Virial_Density_Contrast, only : virialDensityContrast, virialDensityContrastClass
  private

  !![
  <functionClass>
   <name>virialOrbit</name>
   <descriptiveName>Virial Orbits</descriptiveName>
   <description>
    Class providing orbital parameters of satellite halos at the time when they first enter the virial radius of their host.
   </description>
   <default>benson2005</default>
   <method name="orbit" >
    <description>Returns an orbit object.</description>
    <type>type(keplerOrbit)</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>type(treeNode), intent(inout) :: node               , host</argument>
    <argument>logical       , intent(in   ) :: acceptUnboundOrbits</argument>
   </method>
   <method name="velocityDistributionFunction" >
    <description>Returns the distribution function of orbital velocity $\mathrm{d}^2p/\mathrm{d}v_r/\mathrm{d}v_\theta(v_r,v_\theta)$.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout) :: node          , host</argument>
    <argument>double precision          , intent(in   ) :: velocityRadial, velocityTangential</argument>
   </method>
   <method name="isAngularlyResolved" >
    <description>Returns true if this orbit class provides resolution of the orbital angular coordinates ($\theta$,$\phi$) when setting orbits, false otherwise.</description>
    <type>logical</type>
    <pass>yes</pass>
    <code>
     !$GLC attributes unused :: self
     virialOrbitIsAngularlyResolved=.false.
   </code>
   </method>
   <method name="densityContrastDefinition" >
    <description>Returns a {\normalfont \ttfamily virialDensityContrast} object describing the virial density contrast used to define this orbit class.</description>
    <type>class(virialDensityContrastClass)</type>
    <pass>yes</pass>
   </method>
   <method name="velocityTangentialMagnitudeMean" >
    <description>Returns the mean of the magnitude of tangential velocity averaged over all orbits.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node, host</argument>
   </method>
   <method name="velocityTangentialVectorMean" >
    <description>Returns the mean vector of the vector tangential velocity averaged over all orbits.</description>
    <type>double precision, dimension(3)</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node, host</argument>
   </method>
   <method name="angularMomentumMagnitudeMean" >
    <description>Returns the mean of the magnitude of the angular momentum averaged over all orbits.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node, host</argument>
   </method>
   <method name="angularMomentumVectorMean" >
    <description>Returns the mean vector of the angular momentum averaged over all orbits.</description>
    <type>double precision, dimension(3)</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node, host</argument>
   </method>
   <method name="velocityTotalRootMeanSquared" >
    <description>Returns the square root of the mean of the squared total velocity averaged over all orbits.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node, host</argument>
   </method>
   <method name="energyMean" >
    <description>Returns the square root of the mean of the energy averaged over all orbits.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node, host</argument>
   </method>
  </functionClass>
  !!]

end module Virial_Orbits
