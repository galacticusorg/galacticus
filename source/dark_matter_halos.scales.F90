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
Contains a module which provides a class implementing scales of dark matter halo scales.
!!}

module Dark_Matter_Halo_Scales
  !!{
  Provides a class implementing scales of dark matter halo scales.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>darkMatterHaloScale</name>
   <descriptiveName>Dark Matter Halo Scales</descriptiveName>
   <description>Class providing the characteristic physical scales of dark matter halos, including the virial radius,
    virial velocity, virial temperature, mean density, and dynamical timescale. These scales define the boundary and
    characteristic properties of a dark matter halo and are required by many other classes that compute cooling rates,
    star formation, and satellite orbital dynamics.</description>
   <default>virialDensityContrastDefinition</default>
   <method name="timescaleDynamical" >
    <description>The characteristic dynamical timescale of a dark matter halo.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
   <method name="velocityVirial" >
    <description>Returns the virial velocity (in km/s) of the dark matter halo associated with \mono{node}, defined as the circular velocity at the virial radius and providing a characteristic velocity scale for the halo.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
   <method name="velocityVirialGrowthRate" >
    <description>The growth rate of the virial velocity of a dark matter halo.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
   <method name="temperatureVirial" >
    <description>Returns the virial temperature (in K) of the dark matter halo associated with \mono{node}, i.e. the characteristic gas temperature corresponding to the virial velocity, below which gas can be thermally supported against gravitational collapse.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
   <method name="radiusVirial" >
    <description>Returns the virial radius (in Mpc) of the dark matter halo associated with \mono{node}, defined as the radius within which the mean interior density equals a specified overdensity threshold times the critical or mean density.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
   <method name="radiusVirialGradientLogarithmicMass" >
    <description>The logarithmic gradient of virial radius of a dark matter halo with halo mass at fixed epoch.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
   <method name="radiusVirialGrowthRate" >
    <description>Returns the time derivative of the virial radius (in Mpc/Gyr) of the dark matter halo associated with \mono{node}, indicating how rapidly the halo is growing or shrinking in physical size.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
   <method name="densityMean" >
    <description>Returns the mean interior density (in $\mathrm{M}_\odot$ Mpc$^{-3}$) of the dark matter halo associated with \mono{node}, computed as the halo mass divided by its virial volume, representing the average density within the virial radius.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
   <method name="densityMeanGrowthRate" >
    <description>Returns the time derivative of the mean interior density (in $\mathrm{M}_\odot$ Mpc$^{-3}$ Gyr$^{-1}$) of the dark matter halo associated with \mono{node}, reflecting how the balance between mass accretion and volume growth changes over time.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
  </functionClass>
  !!]

end module Dark_Matter_Halo_Scales
