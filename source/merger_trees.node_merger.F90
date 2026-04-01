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
Contains a module which provides a class that implements processing of mergers between nodes.
!!}

module Merger_Trees_Merge_Node
  !!{
  Provides a class that implements processing of mergers between nodes.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>mergerTreeNodeMerger</name>
   <descriptiveName>Merger Tree Node Merger Processing</descriptiveName>
   <description>Class providing models for processing merger tree nodes at the moment of merger---when
    a satellite halo's evolution time reaches that of its parent and the two are combined. The merger
    processor transfers mass, metals, angular momentum, and other properties between the merging
    components according to the chosen prescription, and decides how to restructure the tree (e.g.
    single-level vs.\ multi-level hierarchy). Different implementations correspond to different
    assumptions about satellite orbital decay and the timing of merger events.</description>
   <default>singleLevelHierarchy</default>
   <method name="process" >
    <description>Process the merger between \mono{node} and its parent node, then destroy it.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout), target :: node</argument>
   </method>
  </functionClass>
  !!]

end module Merger_Trees_Merge_Node
