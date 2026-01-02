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
Contains a module which provides a class that implements physical processes.
!!}

module Nodes_Operators
  !!{
  Provides a class that implements physical processes.
  !!}
  use :: Galacticus_Nodes, only : treeNode, interruptTask
  private

  !![
  <functionClass>
    <name>nodeOperator</name>
    <descriptiveName>Node Operators</descriptiveName>
    <description>Class providing operators acting on nodes.</description>
    <default>null</default>
    <method name="nodeInitialize" >
      <description>
	Perform initialization of node properties prior to tree evolution. This typically includes initializations related to the
	evolution of the tree (e.g. growth rates of scale radii, baryonic component initialization). Initializations related to
	the static structure of the tree (e.g. assigning scale radii, merging orbits, etc.) are typically handled by the
	{\normalfont \ttfamily nodeTreeInitialize} method.
      </description>
      <type>void</type>
      <pass>yes</pass>
      <selfTarget>yes</selfTarget>
      <argument>type(treeNode), intent(inout), target :: node</argument>
      <code>
	!$GLC attributes unused :: self, node
      </code>
    </method>
    <method name="nodeTreeInitialize" >
      <description>
	Perform initialization of node properties immediately after tree construction. This typically includes initializations
	related to the static structure of the tree (e.g. assign scale radii, merging orbits, etc.). Initializations related to
	evolution of the tree (e.g. growth rates of scale radii, baryonic component initialization) are typically handled by the
	{\normalfont \ttfamily nodeInitialize} method.
      </description>
      <type>void</type>
      <pass>yes</pass>
      <selfTarget>yes</selfTarget>
      <argument>type(treeNode), intent(inout), target :: node</argument>
      <code>
	!$GLC attributes unused :: self, node
      </code>
    </method>
    <method name="nodesMerge" >
      <description>Act on the merging of two nodes.</description>
      <type>void</type>
      <pass>yes</pass>
      <argument>type(treeNode), intent(inout) :: node</argument>
      <code>
	!$GLC attributes unused :: self, node
      </code>
    </method>
    <method name="nodePromote" >
      <description>Act on the promotion of a node to its parent.</description>
      <type>void</type>
      <pass>yes</pass>
      <argument>type(treeNode), intent(inout) :: node</argument>
      <code>
	!$GLC attributes unused :: self, node
      </code>
    </method>
    <method name="galaxiesMerge" >
      <description>Act on the merging of two galaxies.</description>
      <type>void</type>
      <pass>yes</pass>
      <argument>type(treeNode), intent(inout) :: node</argument>
      <code>
	!$GLC attributes unused :: self, node
      </code>
    </method>
    <method name="differentialEvolutionPre" >
      <description>Operate on a node before differential evolution.</description>
      <type>void</type>
      <pass>yes</pass>
      <argument>type(treeNode), intent(inout) :: node</argument>
      <code>
	!$GLC attributes unused :: self, node
      </code>
    </method>
    <method name="differentialEvolutionScales" >
      <description>Set scales for meta-properties evolved differentially.</description>
      <type>void</type>
      <pass>yes</pass>
      <argument>type(treeNode), intent(inout) :: node</argument>
      <code>
	!$GLC attributes unused :: self, node
      </code>
    </method>
    <method name="differentialEvolutionAnalytics" >
      <description>Mark (meta-)properties as analytically solved in differential evolution.</description>
      <type>void</type>
      <pass>yes</pass>
      <argument>type(treeNode), intent(inout) :: node</argument>
      <code>
	!$GLC attributes unused :: self, node
      </code>
    </method>
    <method name="differentialEvolutionInactives" >
      <description>Mark (meta-)properties as inactive in differential evolution.</description>
      <type>void</type>
      <pass>yes</pass>
      <argument>type(treeNode), intent(inout) :: node</argument>
      <code>
	!$GLC attributes unused :: self, node
      </code>
    </method>
    <method name="differentialEvolution" >
      <description>Operate on a node during differential evolution.</description>
      <type>void</type>
      <pass>yes</pass>
      <selfTarget>yes</selfTarget>
      <argument>type     (treeNode     ), intent(inout), target  :: node             </argument>
      <argument>logical                 , intent(inout)          :: interrupt        </argument>
      <argument>procedure(interruptTask), intent(inout), pointer :: functionInterrupt</argument>
      <argument>integer                 , intent(in   )          :: propertyType     </argument>
      <code>
	!$GLC attributes unused :: self, node, interrupt, functionInterrupt, propertyType
      </code>
    </method>
    <method name="differentialEvolutionSolveAnalytics" >
      <description>Set the values of analytically-solvable properties of a node during differential evolution.</description>
      <type>void</type>
      <pass>yes</pass>
      <argument>type            (treeNode), intent(inout) :: node</argument>
      <argument>double precision          , intent(in   ) :: time</argument>
      <code>
	!$GLC attributes unused :: self, node, time
      </code>
    </method>
    <method name="predeterminedSolveAnalytics" >
      <description>Set the values of analytically-solvable pre-determined properties of a node.</description>
      <type>void</type>
      <pass>yes</pass>
      <argument>type            (treeNode), intent(inout) :: node</argument>
      <argument>double precision          , intent(in   ) :: time</argument>
      <code>
	!$GLC attributes unused :: self, node, time
      </code>
    </method>
    <method name="differentialEvolutionStepFinalState" >
      <description>Operate on a node at the end of an ODE step.</description>
      <type>void</type>
      <pass>yes</pass>
      <argument>type(treeNode), intent(inout) :: node</argument>
      <code>
	!$GLC attributes unused :: self, node
      </code>
    </method>
    <method name="differentialEvolutionPost" >
      <description>Operate on a node after differential evolution.</description>
      <type>void</type>
      <pass>yes</pass>
      <argument>type(treeNode), intent(inout) :: node</argument>
      <code>
	!$GLC attributes unused :: self, node
      </code>
    </method>
    <method name="differentialEvolutionPostStep" >
      <description>Operate on a node after a differential evolution step.</description>
      <type>void</type>
      <pass>yes</pass>
      <argument>type   (treeNode), intent(inout) :: node</argument>
      <argument>integer          , intent(inout) :: status</argument>
      <code>
	!$GLC attributes unused :: self, node, status
      </code>
    </method>
  </functionClass>
  !!]

end module Nodes_Operators
