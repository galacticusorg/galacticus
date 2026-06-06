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
Contains a module which implements a class for calculations for progenitor properties for mergers.
!!}

module Satellite_Merging_Progenitor_Properties
  !!{
  Implements a class for calculations for progenitor properties for mergers.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  implicit none
  private

  !![
  <functionClass>
   <name>mergerProgenitorProperties</name>
   <descriptiveName>Merger Progenitor Properties</descriptiveName>
   <description>Class providing models of the effective properties of merger progenitors---the masses, radii,
    and angular momentum factors of the satellite and host galaxies immediately before a merger event,
    which are passed to the remnant size calculator. The progenitor properties determine the energy
    budget available to the merger remnant: the satellite mass, host spheroid mass, pre-merger host
    spheroid mass, effective radii, angular momentum factor, and the expected remnant spheroid and
    gas-spheroid masses used to compute the post-merger structure.</description>
   <default>standard</default>
   <method name="get" >
    <description>Calculates the effective masses, radii, and angular momentum factors of the satellite and host galaxy progenitors immediately before a merger event, providing the energy budget inputs required by remnant size calculators.</description>
    <type>void</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>type            (treeNode), intent(inout), target :: nodeSatellite            , nodeHost              </argument>
    <argument>double precision          , intent(  out)         :: massSatellite            , massHost              </argument>
    <argument>double precision          , intent(  out)         :: massSpheroidSatellite    , massSpheroidHost      </argument>
    <argument>double precision          , intent(  out)         :: massSpheroidHostPreMerger, radiusSatellite       </argument>
    <argument>double precision          , intent(  out)         :: radiusHost               , factorAngularMomentum </argument>
    <argument>double precision          , intent(  out)         :: massSpheroidRemnant      , massGasSpheroidRemnant</argument>
   </method>
  </functionClass>
  !!]

end module Satellite_Merging_Progenitor_Properties
