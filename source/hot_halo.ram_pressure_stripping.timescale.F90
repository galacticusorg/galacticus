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
Contains a module that implements a class for calculations of ram pressure stripping timescales for hot halos.
!!}

module Hot_Halo_Ram_Pressure_Stripping_Timescales
  !!{
  Implements a class for calculations of ram pressure stripping timescales for hot halos.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>hotHaloRamPressureTimescale</name>
   <descriptiveName>Hot Halo Ram Pressure Timescales</descriptiveName>
   <description>Class providing models of the ram pressure stripping timescale (in Gyr) for the hot gas
    atmosphere of a satellite galaxy orbiting in its host halo. Instead of computing an instantaneous
    stripping radius, the timescale approach allows for gradual stripping of hot gas over time, capturing
    the orbital history and the time required for ram pressure to overcome the satellite's self-gravity.
    The timescale may depend on the local ram pressure, the orbital velocity, or the halo dynamical time.</description>
   <default>ramPressureAcceleration</default>
   <method name="timescale" >
    <description>Return the ram pressure stripping timescale for \mono{node} (in units of Gyr).</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
  </functionClass>
  !!]

end module Hot_Halo_Ram_Pressure_Stripping_Timescales
