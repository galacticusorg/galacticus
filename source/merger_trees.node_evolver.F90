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
Contains a module which provides a class that implements evolution of nodes.
!!}

module Merger_Trees_Evolve_Node
  !!{
  Provides a class that implements evolution of nodes.
  !!}
  use :: Galactic_Structure_Solvers, only : galacticStructureSolverClass
  use :: Galacticus_Nodes          , only : interruptTask               , mergerTree, treeNode
  use :: Locks                     , only : ompLockClass
  private

  !![
  <functionClass>
   <name>mergerTreeNodeEvolver</name>
   <descriptiveName>Merger Tree Node Evolvers</descriptiveName>
   <description>Class providing evolvers for nodes in merger trees.</description>
   <default>standard</default>
   <method name="evolve" >
    <description>Evolve a node merger tree.</description>
    <type>void</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>type            (mergerTree                  ), intent(inout)           :: tree                     </argument>
    <argument>type            (treeNode                    ), intent(inout), pointer  :: node                     </argument>
    <argument>double precision                              , intent(in   )           :: timeEnd                  </argument>
    <argument>logical                                       , intent(  out)           :: interrupted              </argument>
    <argument>procedure       (interruptTask               ), intent(  out), pointer  :: functionInterrupt        </argument>
    <argument>class           (galacticStructureSolverClass), intent(in   ), target   :: galacticStructureSolver__</argument>
    <argument>class           (ompLockClass                ), intent(inout)           :: treeLock                 </argument>
    <argument>integer         (kind_int8                   ), intent(in   ), optional :: systemClockMaximum       </argument>
    <argument>integer                                       , intent(  out), optional :: status                   </argument>
   </method>
   <method name="promote" >
    <description>Promote {\normalfont \ttfamily node} to its parent node, then destroy it.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout), pointer :: node</argument>
   </method>
   <method name="merge" >
    <description>Handles instances where {\normalfont \ttfamily node} is about to merge with its parent node.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
   <method name="isAccurate" >
    <description>Return true if a tree node property is within expected accuracy of a given value.</description>
    <type>logical</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: valueNode, valueExpected</argument>
   </method>
  </functionClass>
  !!]

end module Merger_Trees_Evolve_Node
