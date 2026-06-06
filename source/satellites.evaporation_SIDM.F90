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
Contains a module that implements calculations of the evaporation of orbiting satellites due to dark matter self-interactions.
!!}

module Satellite_Evaporation_SIDM
  !!{
  Implements calculation of the evaporation of orbiting satellites due to dark matter self-interactions.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>satelliteEvaporationSIDM</name>
   <descriptiveName>Dark Matter Self-interaction Satellite Evaporation</descriptiveName>
   <description>Class providing models of the evaporation of orbiting satellites due to dark matter self-interactions
    (SIDM)---the mass loss from satellite halos caused by energetic scattering events between satellite dark
    matter particles and those of the host halo. High-momentum-transfer scatters can eject satellite particles,
    gradually eroding the subhalo on top of tidal stripping. Implementations return the satellite mass loss rate
    (in $\mathrm{M}_\odot$~Gyr$^{-1}$) from SIDM-driven evaporation, which depends on the self-interaction cross-section,
    the satellite-host relative velocity, and the local number density of host particles.</description>
   <default>kummer2018</default>
   <method name="massLossRate" >
    <description>Returns the satellite mass loss rate due to dark matter self-interactions (in units of $\mathrm{M}_\odot$/Gyr).</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
  </functionClass>
  !!]

end module Satellite_Evaporation_SIDM
