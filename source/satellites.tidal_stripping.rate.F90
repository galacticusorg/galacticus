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
Contains a module that provides a class to perform calculations of the mass loss rate due to tidal stripping for satellites.
!!}

module Satellite_Tidal_Stripping
  !!{
  Provides a class to perform calculations of the mass loss rate due to tidal stripping for satellites.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>satelliteTidalStripping</name>
   <descriptiveName>Tidal Stripping of Satellites</descriptiveName>
   <description>Class providing models of tidal stripping for satellites---the gravitational removal of dark matter
    and stellar mass from satellite halos as they orbit through the tidal field of their host halo. The tidal force
    strips material outside the tidal radius, reducing the satellite mass over time at a rate (in $\mathrm{M}_\odot$~Gyr$^{-1}$)
    that depends on the satellite's orbit, concentration, and the host potential. This mass loss sets the subhalo
    abundance, the galaxy-to-halo mass ratio in satellites, and drives the evolution of satellite galaxies
    toward quiescence.</description>
   <default>zentner2005</default>
   <method name="massLossRate" >
    <description>Returns the rate of tidal mass loss for \mono{node} (in units of $\mathrm{M}_\odot$/Gyr).</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
  </functionClass>
  !!]

end module Satellite_Tidal_Stripping
