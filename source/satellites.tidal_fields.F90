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
Contains a module that implements calculations of tidal fields acting on satellites.
!!}

module Satellites_Tidal_Fields
  !!{
  Implements calculations of tidal fields acting on satellites.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>satelliteTidalField</name>
   <descriptiveName>Satellite halo tidal field models.</descriptiveName>
   <description>
    Class providing models of tidal fields experienced by satellite halos.
   </description>
   <default>null</default>
   <method name="tidalTensorRadial" >
    <description>Returns the radial component, $\Phi_\mathrm{rr}$, of the tidal tensor, $\Phi_\mathrm{ab}$.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
  </functionClass>
  !!]

end module Satellites_Tidal_Fields
