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
  use :: Tensors         , only : tensorRank2Dimension3Symmetric
  private

  !![
  <functionClass>
   <name>satelliteTidalField</name>
   <descriptiveName>Satellite Tidal Fields</descriptiveName>
   <description>Class providing models of tidal fields experienced by satellite halos---the second-order spatial
    derivatives of the gravitational potential that characterize the differential force stretching or compressing
    the satellite. The key quantity returned is the radial component $\Phi_\mathrm{rr}$ of the tidal tensor
    $\Phi_\mathrm{ab} = \partial^2\Phi/\partial x_a \partial x_b$, evaluated at the satellite's orbital position.
    This tidal field drives tidal heating and mass loss, and is used by tidal stripping rate calculations to
    determine how rapidly material is removed from the satellite.</description>
   <default>standard</default>
   <method name="tidalTensor" >
    <description>Returns the tidal tensor, $\Phi_\mathrm{ab}$.</description>
    <type>type(tensorRank2Dimension3Symmetric)</type>
    <pass>yes</pass>
    <argument>type   (treeNode), intent(inout)                   :: node                                        </argument>
    <argument>type   (treeNode), intent(inout), optional, target :: nodeHost                                    </argument>
    <argument>logical          , intent(in   ), optional         :: atPericenter, includeCentrifugalAcceleration</argument>
   </method>
   <method name="tidalTensorRadial" >
    <description>Returns the radial component, $\Phi_\mathrm{rr}$, of the tidal tensor, $\Phi_\mathrm{ab}$.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type   (treeNode), intent(inout)                   :: node                                        </argument>
    <argument>type   (treeNode), intent(inout), optional, target :: nodeHost                                    </argument>
    <argument>logical          , intent(in   ), optional         :: atPericenter, includeCentrifugalAcceleration</argument>
   </method>
   <method name="tidalTensorDominant" >
    <description>Returns the dominant eigenvalue of the tidal tensor, $\Phi_\mathrm{ab}$.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type   (treeNode), intent(inout)                   :: node                                        </argument>
    <argument>type   (treeNode), intent(inout), optional, target :: nodeHost                                    </argument>
    <argument>logical          , intent(in   ), optional         :: atPericenter, includeCentrifugalAcceleration</argument>
   </method>
  </functionClass>
  !!]

end module Satellites_Tidal_Fields
