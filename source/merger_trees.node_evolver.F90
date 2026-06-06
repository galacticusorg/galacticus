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
   <description>Class providing evolvers for individual nodes in a merger tree---the objects responsible
    for integrating the ODEs that govern the physical properties of a single halo or galaxy from one
    timestep to the next. The node evolver advances all state variables (mass, spin, metallicity, etc.)
    using the relevant physics modules and handles special events such as node promotion (when a node's
    time equals its parent's) and pre-merge processing. It also checks whether node properties have
    converged to the required accuracy.</description>
   <default>standard</default>
   <method name="evolve" >
    <description>Advance the physical properties of the given \mono{node} from its current time to \mono{timeEnd} by integrating the relevant ODE physics, setting \mono{interrupted} and \mono{functionInterrupt} if the evolution must be paused for an event such as node promotion or a merger.</description>
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
    <description>Promote \mono{node} to its parent node, then destroy it.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout), pointer :: node</argument>
   </method>
   <method name="merge" >
    <description>Handles instances where \mono{node} is about to merge with its parent node.</description>
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
