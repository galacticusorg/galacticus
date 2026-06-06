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

  !+    Contributions to this file made by: Matías Liempi

!!{
Contains a module which implements a class for black hole seeds.
!!}

module Black_Hole_Seeds
  !!{
  Implements a class for black hole seeds.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  implicit none
  private

  !![
  <enumeration>
    <name>blackHoleFormationChannel</name>
    <description>Enumeration of the physical channels through which black hole seeds can form, including undetermined origin and formation via stellar cluster collapse leading to a massive seed black hole.</description>
    <decodeFunction>yes</decodeFunction>
    <validator>yes</validator>
    <visibility>public</visibility>
    <entry label="undetermined"       />
    <entry label="starClusterCollapse"/>
  </enumeration>
  !!]
  
  !![
  <functionClass>
   <name>blackHoleSeeds</name>
   <descriptiveName>Black Hole Seeds</descriptiveName>
   <description>
    Class providing models of the initial seed masses for supermassive black holes and their formation channel.
    Black hole seeds are the initial conditions for black hole growth via accretion and mergers, and their masses
    and formation channels (e.g. stellar collapse, direct collapse, star cluster collapse) are set when a new
    black hole component is initialized in a halo.
   </description>
   <default>fixed</default>
   <method name="mass" >
    <description>Computes the mass of the black hole seed in the given \mono{node}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
   <method name="spin" >
    <description>Computes the spin of the black hole seed in the given \mono{node}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
   <method name="formationChannel">
    <description>Returns the formation channel of the seed in the given \mono{node}.</description>
    <type>type(enumerationBlackHoleFormationChannelType)</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
    <code>
      !$GLC attributes unused :: self, node
      blackHoleSeedsFormationChannel=blackHoleFormationChannelUndetermined
    </code>
   </method>
  </functionClass>
  !!]
  
end module Black_Hole_Seeds
