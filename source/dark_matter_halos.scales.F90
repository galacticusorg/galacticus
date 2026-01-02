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
   <description>Class providing dark matter halo scales.</description>
   <default>virialDensityContrastDefinition</default>
   <method name="timescaleDynamical" >
    <description>The characteristic dynamical timescale of a dark matter halo.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
   <method name="velocityVirial" >
    <description>The virial velocity of a dark matter halo</description>
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
    <description>The virial temperature of a dark matter halo</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
   <method name="radiusVirial" >
    <description>The virial radius of a dark matter halo.</description>
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
    <description>The growth rate of the virial radius of a dark matter halo.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
   <method name="densityMean" >
    <description>The mean density of a dark matter halo.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
   <method name="densityMeanGrowthRate" >
    <description>The growth rate of the mean density of a dark matter halo.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
  </functionClass>
  !!]

end module Dark_Matter_Halo_Scales
