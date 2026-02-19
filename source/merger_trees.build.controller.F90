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
Contains a module which provides controller objects for building merger trees.
!!}

module Merger_Tree_Build_Controllers
  !!{
  Provides controller objects for building merger trees.
  !!}
  use :: Galacticus_Nodes     , only : treeNode
  use :: Merger_Tree_Walkers  , only : mergerTreeWalkerClass
  use :: Merger_Tree_Branching, only : mergerTreeBranchingProbabilityClass
  private

  !![
  <functionClass>
   <name>mergerTreeBuildController</name>
   <descriptiveName>Merger Tree Build Controllers</descriptiveName>
   <description>Class providing merger tree build controllers.</description>
   <default>uncontrolled</default>
   <method name="control" >
    <description>Control the behavior of a tree build.</description>
    <type>logical</type>
    <pass>yes</pass>
    <argument>type (treeNode             ), intent(inout), pointer  :: node       </argument>
    <argument>class(mergerTreeWalkerClass), intent(inout), optional :: treeWalker_</argument>
   </method>
   <method name="timeMinimum" >
    <description>Return the minimum ``time'' (using the usual $w$ variable for merger tree building) allowed for this node.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout) :: node                                 </argument>
    <argument>double precision          , intent(in   ) :: massBranch, criticalOverdensityBranch</argument>
    <code>
      !$GLC attributes unused :: self, node, massBranch, criticalOverdensityBranch
      ! No limit to time by default.
      mergerTreeBuildControllerTimeMinimum=criticalOverdensityBranch
    </code>
   </method>   
   <method name="timeMaximum" >
    <description>Return the maximum ``time'' (using the usual $w$ variable for merger tree building) allowed for this node.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout) :: node                                                </argument>
    <argument>double precision          , intent(in   ) :: massBranch, criticalOverdensityBranch, timeReference</argument>
    <argument>logical                   , intent(  out) :: insertNode                                          </argument>
    <code>
      !$GLC attributes unused :: self, node, massBranch, criticalOverdensityBranch, timeReference
      ! No limit to time by default.
      mergerTreeBuildControllerTimeMaximum=huge(0.0d0)
      ! Do not force creation of a node by default.
      insertNode=.false.
    </code>
   </method>   
   <method name="controlTimeMaximum" >
    <description>Control the behavior of a tree build when the maximum time for a node is reached.</description>
    <type>logical</type>
    <pass>yes</pass>
    <argument>type            (treeNode ), intent(inout), target :: node                                 </argument>
    <argument>double precision           , intent(in   )         :: massBranch, criticalOverdensityBranch</argument>
   <argument> integer         (kind_int8), intent(inout)         :: nodeIndex                            </argument>
    <code>
      !$GLC attributes unused :: self, node, massBranch, criticalOverdensityBranch
      mergerTreeBuildControllerControlTimeMaximum=.true.
    </code>
   </method>
   <method name="branchingProbabilityObject" >
    <description>Return a branching probability object to use in tree building.</description>
    <type>class(mergerTreeBranchingProbabilityClass)</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
   <method name="nodesInserted" >
    <description>Alert the controller when new nodes are inserted into the tree.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type   (treeNode), intent(inout)           :: nodeCurrent    , nodeProgenitor1</argument>
    <argument>type   (treeNode), intent(inout), optional :: nodeProgenitor2                 </argument>
    <argument>logical          , intent(in   ), optional :: didBranch                       </argument>
   </method>
  </functionClass>
  !!]

end module Merger_Tree_Build_Controllers
