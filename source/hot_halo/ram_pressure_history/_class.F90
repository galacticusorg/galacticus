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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a module that implements a class which provides the history of the ram pressure force experienced by satellites.
!!}

module Hot_Halo_Ram_Pressure_Histories
  !!{RST
  Implements a class which provides the history of the ram pressure force experienced by satellites.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass docformat="rst">
   <name>hotHaloRamPressureHistory</name>
   <descriptiveName>Hot Halo Ram Pressure History</descriptiveName>
   <description>
   Class providing the history of the ram pressure force (in :math:`(\mathrm{km\,s}^{-1})^2 \mathrm{M}_\odot \mathrm{Mpc}^{-3}`)
   exerted on a satellite galaxy by the hot gas atmosphere of its host halo, as needed by models of ram pressure stripping in which
   the outcome depends on the duration of the encounter as well as on its peak strength :cite:`koppen_ram_2018`. Two quantities
   are provided: the running maximum of the ram pressure, :math:`P`, and its time integral, :math:`I`, both accumulated since the
   satellite fell into its host. Implementations may apply a decay to both, so that they refer, effectively, to the most recent
   pericentric passage rather than to the whole orbital history. Both are zero for nodes which are not satellites.
   </description>
   <default>accumulated</default>
   <method name="forceMaximum" >
    <description>
    Return the (possibly decaying) running maximum of the ram pressure force experienced by ``node`` (in units of
    (km/s)\ :math:`^2\mathrm{M}_\odot`/Mpc\ :math:`^3`).
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
   <method name="forceTimeIntegrated" >
    <description>
    Return the (possibly decaying) time integral of the ram pressure force experienced by ``node`` (in units of
    (km/s)\ :math:`^2\mathrm{M}_\odot`/Mpc\ :math:`^3` Gyr).
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
  </functionClass>
  !!]

end module Hot_Halo_Ram_Pressure_Histories
