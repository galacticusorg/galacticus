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
Contains a module that implements calculations of the time available for cooling.
!!}

module Cooling_Times_Available
  !!{
  Provides a class that implements calculations of the time available for cooling.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>coolingTimeAvailable</name>
   <descriptiveName>Cooling Time Available</descriptiveName>
   <description>Class providing models of the time available for cooling (in Gyr)---the elapsed time since gas
    in the hot halo was first able to start cooling. This quantity, together with the cooling time, determines
    whether gas has had enough time to cool and fall in. Implementations typically anchor this time to halo
    formation or to the time since the last major merger, and different choices lead to significantly different
    predictions for the cold gas supply and star formation history of galaxies.</description>
   <default>whiteFrenk1991</default>
   <method name="timeAvailable" >
    <description>Return the time available for cooling in \mono{node} in units of Gyr.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
   <method name="timeAvailableIncreaseRate" >
    <description>Return the rate at which the time available for cooling increases in \mono{node} (dimensionless).</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
  </functionClass>
  !!]

end module Cooling_Times_Available
