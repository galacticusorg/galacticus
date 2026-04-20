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
Contains a module that implements a class for initializing the bound mass of satellite halos.
!!}

module Satellites_Bound_Mass_Initialize
  !!{
  Implements a class for initializing the bound mass of satellite halos.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>satelliteMassBoundInitializor</name>
   <descriptiveName>Satellite Bound Mass Initializor</descriptiveName>
   <description>Class providing the initial bound mass of satellite halos when they first become satellites.</description>
   <default>basicMass</default>
   <method name="massBound">
    <description>Returns the initial bound mass of the satellite halo associated with \mono{node} (in units of $\mathrm{M}_\odot$).</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
  </functionClass>
  !!]

end module Satellites_Bound_Mass_Initialize
