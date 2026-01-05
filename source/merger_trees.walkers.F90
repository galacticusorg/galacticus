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

module Merger_Tree_Walkers
  !!{
  Provides a class of walker objects for merger trees.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>mergerTreeWalker</name>
   <descriptiveName>Merger Tree Walkers</descriptiveName>
   <description>Class providing walkers for merger trees. Walkers iterate over nodes in a tree.</description>
   <default>isolatedNodes</default>
   <method name="next" >
    <description>Update the pointer to the next node to visit. Returns true if such a node exists, returns false if no such node exists (i.e. if all nodes have been visited already).</description>
    <type>logical</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout), pointer :: node</argument>
   </method>
   <method name="nodesRemain" >
    <description>Returns true if more nodes remain to be walked to in the tree.</description>
    <type>logical</type>
    <pass>yes</pass>
   </method>
  </functionClass>
  !!]

end module Merger_Tree_Walkers
