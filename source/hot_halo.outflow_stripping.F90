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
Contains a module which provides a class that implements stripping of outflowed mass in the hot halo.
!!}

module Hot_Halo_Outflows_Stripping
  !!{
  Provides a class that implements stripping of outflowed mass in the hot halo.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>hotHaloOutflowStripping</name>
   <descriptiveName>Hot Halo Outflow Reincorporation</descriptiveName>
   <description>Class providing models of the stripping of outflowed (ejected) gas from the hot halo of a
    satellite galaxy as it orbits within its host halo. When a galaxy becomes a satellite the ejected gas
    reservoir associated with it may be stripped by tidal or ram pressure forces from the host halo. The
    fraction of outflowed gas stripped (returned to the host halo's hot gas) depends on the orbital position
    and the relative pressures, affecting the satellite's subsequent star formation and stellar mass.</description>
   <default>zero</default>
   <method name="neverStripped" >
    <description>Return true if outflows are never stripped from this halo.</description>
    <type>logical</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
   <method name="fractionStripped" >
    <description>Return the fraction of outflowing material that should be stripped from the halo.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
  </functionClass>
  !!]

end module Hot_Halo_Outflows_Stripping
