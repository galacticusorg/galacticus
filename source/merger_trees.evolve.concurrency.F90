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
Contains a module that implements a class for merger tree evolution concurrency.
!!}

module Merger_Trees_Evolve_Concurrency
  !!{
  Provides a class that implements merger tree evolution concurrency.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>mergerTreeEvolveConcurrency</name>
   <descriptiveName>Merger Tree Evolution Concurrency</descriptiveName>
   <description>Class providing logic for merger tree evolution concurrency.</description>
   <default>halosSubhalos</default>
   <method name="initializeTree">
    <description>Initialize concurrency for a new tree.</description>
    <type>void</type>
    <pass>yes</pass>
   </method>
   <method name="countPhases">
    <description>Return a count of the number of evolution phases.</description>
    <type>integer</type>
    <pass>yes</pass>
   </method>
   <method name="includeInEvolution">
    <description>Return true if the node should be included in the given phase of merger tree evolution.</description>
    <type>logical</type>
    <pass>yes</pass>
    <argument>integer          , intent(in   )         :: evolutionPhase</argument>
    <argument>type   (treeNode), intent(inout), target :: node          </argument>
   </method>
  </functionClass>
  !!]

end module Merger_Trees_Evolve_Concurrency
