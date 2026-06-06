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
Contains a module which implements a class for computing radii of galactic components (or more general components).
!!}

module Galactic_Structure_Solvers
  !!{
  Implements a class for calculations of sizes of galactic components (or more general components).
  !!}
  use :: Galacticus_Nodes, only : treeNode
  implicit none
  private

  !![
  <functionClass>
   <name>galacticStructureSolver</name>
   <descriptiveName>Galactic Structure Solver</descriptiveName>
   <description>Class providing solvers for galactic structure---algorithms that determine the sizes
    (half-mass or scale radii) of galactic components (disc, bulge, nuclear star cluster), e.g. by finding
    the equilibrium configuration that satisfies the angular momentum conservation constraint, i.e.
    the radius at which the specific angular momentum of the component matches that of the infalling
    gas. Solving for the structure may require self-consistent accounting of the gravitational potential
    from all components including dark matter, gas, and stars.</description>
   <default>equilibrium</default>
   <method name="solve" >
    <description>Solves for the structure of components in the given \mono{node}.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type   (treeNode), intent(inout), target   :: node            </argument>
    <argument>logical          , intent(in   ), optional :: plausibilityOnly</argument>
   </method>
   <method name="revert" >
    <description>Revert the structure of components in the given \mono{node} (if necessary to ensure that the structure solver will give the same result when called consecutively).</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
  </functionClass>
  !!]

end module Galactic_Structure_Solvers
