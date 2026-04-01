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

!+    Contributions to this file made by:  Anthony Pullen, Andrew Benson.

!!{
Contains a module that provides a class to perform calculations of the tidal stripping radius for satellites.
!!}

module Satellite_Tidal_Stripping_Radii
  !!{
  Provides a class to perform calculations of the tidal stripping radius for satellites.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>satelliteTidalStrippingRadius</name>
   <descriptiveName>Tidal Stripping Radii</descriptiveName>
   <description>Class providing models of tidal stripping radii for satellites---the radius (in Mpc) within which
    material remains gravitationally bound to the satellite against the tidal field of the host halo. Beyond this
    radius the tidal force exceeds the satellite's self-gravity, so mass is stripped away. The tidal radius depends
    on the satellite mass, the host density at the satellite's orbital position, and the choice of tidal criterion
    (e.g.\ King 1962 or Jacobi radius). It sets the outer boundary used by tidal stripping rate models.</description>
   <default>king1962</default>
   <method name="radius" >
    <description>Returns the tidal stripping radius (in units of Mpc) for the given satellite node, i.e., the radius within which material remains gravitationally bound against the tidal field of the host halo.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>type(treeNode), intent(inout), target :: node</argument>
   </method>
  </functionClass>
  !!]

end module Satellite_Tidal_Stripping_Radii
