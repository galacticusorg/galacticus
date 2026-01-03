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
Contains a module that implements a class for computing tidal heating rates for satellites.
!!}

module Satellite_Tidal_Heating
  !!{
  Implements a class for calculations of tidal heating for satellites.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>satelliteTidalHeatingRate</name>
   <descriptiveName>Satellite halo tidal heating rate models.</descriptiveName>
   <description>
    Class providing models of tidal heating rates in satellite halos. Specifically, the integrated, normalized (i.e. the energy
    divided by radius squared) tidal heating energy, $Q_\mathrm{tidal}$.
   </description>
   <default>zero</default>
   <method name="heatingRate" >
    <description>Return the satellite tidal heating rate for {\normalfont \ttfamily node} (in units of (km/s/Mpc)$^2$/Gyr).</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
  </functionClass>
  !!]

end module Satellite_Tidal_Heating
