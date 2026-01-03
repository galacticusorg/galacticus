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
Contains a module that provides and object that implements satellite merging timescales.
!!}

module Satellite_Merging_Timescales
  !!{
  Provides and object that implements satellite merging timescales.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  use :: Kepler_Orbits   , only : keplerOrbit
  private

  ! Effective infinite time for merging. This is set to a fraction of the largest representable number. The fraction is included
  ! such that if perturbations are made around this value it does not cause floating point exceptions.
  double precision, public, parameter :: satelliteMergeTimeInfinite=1.0d-6*huge(1.0d0)

  !![
  <functionClass>
   <name>satelliteMergingTimescales</name>
   <descriptiveName>Satellite Merging Timescales</descriptiveName>
   <description>
    Object providing merging timescales for satellites.
   </description>
   <default>jiang2008</default>
   <method name="timeUntilMerging" >
    <description>Return the time (in Gyr) until the satellite will merge with its host given the current orbit.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode   ), intent(inout) :: node</argument>
    <argument>type(keplerOrbit), intent(inout) :: orbit</argument>
   </method>
  </functionClass>
  !!]

end module Satellite_Merging_Timescales
