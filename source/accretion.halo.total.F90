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
Contains a module which provides a class for calculations of the intergalactic medium thermal and ionization state.
!!}

module Accretion_Halo_Totals
  !!{
  Provides a class for calculations of the total accretion rate onto halos for use by the halo accretion classes which compute
  the accretion rates of baryonic material.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>accretionHaloTotal</name>
   <descriptiveName>Halo Total Accretion Rates</descriptiveName>
   <description>
    Class providing total accretion rates onto halos, i.e. the mass which would be accreted in a dark matter-only universe.
   </description>
   <default>simple</default>
   <method name="accretionRate" >
    <description>Return the total accretion rate onto the given {\normalfont \ttfamily node}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
   <method name="accretedMass" >
    <description>Return the total accreted mass in the given {\normalfont \ttfamily node}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
  </functionClass>
  !!]

end module Accretion_Halo_Totals
