!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
Contains a module which implements a class for determining how mass is moved around as a consequence of a satellite merging
event.
!!}

module Satellite_Merging_Mass_Movements
  !!{
  Implements a class for determining how mass is moved around as a consequence of a satellite merging event.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  implicit none
  private

  !![
  <functionClass>
   <name>mergerMassMovements</name>
   <descriptiveName>Merger Mass Movements</descriptiveName>
   <description>
    Class providing models of the movements of mass during mergers.
   </description>
   <default>simple</default>
   <method name="get" >
    <description>Determine movements of mass during mergers.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type   (treeNode                        ), intent(inout), target :: node                                                                                        </argument>
    <argument>type   (enumerationDestinationMergerType), intent(  out)         :: destinationGasSatellite, destinationStarsSatellite, destinationGasHost, destinationStarsHost</argument>
    <argument>logical                                  , intent(  out)         :: mergerIsMajor                                                                               </argument>
   </method>
  </functionClass>
  !!]

  !![
  <enumeration>
   <name>destinationMerger</name>
   <description>Enumeration of possible destinations for mass in mergers.</description>
   <encodeFunction>yes</encodeFunction>
   <entry label="unmoved" />
   <entry label="dominant"/>
   <entry label="disk"    />
   <entry label="spheroid"/>
  </enumeration>
  !!]

end module Satellite_Merging_Mass_Movements
