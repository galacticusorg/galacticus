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
Contains a module that provides a class for distributions of orphan satellites.
!!}

module Satellite_Oprhan_Distributions
  !!{
  Provides a class for dark matter halo spin distributions.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>satelliteOrphanDistribution</name>
   <descriptiveName>Satellite Orphan Distributions</descriptiveName>
   <description>Class providing distributions for orphan satellites.</description>
   <default>traceDarkMatter</default>
   <method name="extent" >
    <description>The maximum extent of the distribution, i.e. the radius of a sphere centered on the host halo which encompasses all orphan satellites.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
   <method name="position" >
    <description>Return the position of the given orphan in physical coordinates.</description>
    <type>double precision, dimension(3)</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
   <method name="velocity" >
    <description>Return the peculiar velocity of the given orphan in physical coordinates.</description>
    <type>double precision, dimension(3)</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
  </functionClass>
  !!]

end module Satellite_Oprhan_Distributions
