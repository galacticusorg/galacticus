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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a module which provides a class that implements physical processes.
!!}

module Nodes_Operators
  !!{RST
  Provides a class that implements physical processes.
  !!}
  use :: Galacticus_Nodes, only : treeNode, interruptTask
  private

  !![
  <enumeration docformat="rst">
   <name>nodeOperatorMethod</name>
   <description>
   Enumeration of the methods of the :galacticus-class:`nodeOperatorClass`, used to index the ``methodIsNull`` flags which record which methods fall through to the default (null) implementation for a given object.
   </description>
   <validator>yes</validator>
   <indexing>1</indexing>
   <entry label="nodeInitialize"                     />
   <entry label="nodeTreeInitialize"                 />
   <entry label="nodesMerge"                         />
   <entry label="nodePromote"                        />
   <entry label="galaxiesMerge"                      />
   <entry label="differentialEvolutionPre"           />
   <entry label="differentialEvolutionScales"        />
   <entry label="differentialEvolutionAnalytics"     />
   <entry label="differentialEvolutionInactives"     />
   <entry label="differentialEvolution"              />
   <entry label="differentialEvolutionSolveAnalytics"/>
   <entry label="predeterminedSolveAnalytics"        />
   <entry label="differentialEvolutionStepFinalState"/>
   <entry label="differentialEvolutionPost"          />
   <entry label="differentialEvolutionPostStep"      />
  </enumeration>
  !!]

  !![
  <functionClass docformat="rst">
    <name>nodeOperator</name>
    <descriptiveName>Node Operators</descriptiveName>
    <description>
    Class providing operators that act on merger tree nodes to implement physical processes. Node operators are called at key moments during tree evolution---at tree construction, node initialization, derivative evaluation, node merging, and promotion---to apply physics such as star formation, feedback, black hole growth, and structural evolution. Multiple operators can be stacked in a multi-operator list to apply several physical processes simultaneously.
    </description>
    <default>null</default>
    <method name="nodeInitialize" >
      <description>
      Perform initialization of node properties prior to tree evolution. This typically includes initializations related to the evolution of the tree (e.g. growth rates of scale radii, baryonic component initialization). Initializations related to the static structure of the tree (e.g. assigning scale radii, merging orbits, etc.) are typically handled by the ``nodeTreeInitialize`` method.
      </description>
      <type>void</type>
      <pass>yes</pass>
      <selfTarget>yes</selfTarget>
      <argument>type(treeNode), intent(inout), target :: node</argument>
      <code>
	!$GLC attributes unused :: node
	self%methodIsNull(nodeOperatorMethodNodeInitialize%ID)=.true.
      </code>
    </method>
    <method name="nodeTreeInitialize" >
      <description>
      Perform initialization of node properties immediately after tree construction. This typically includes initializations related to the static structure of the tree (e.g. assign scale radii, merging orbits, etc.). Initializations related to evolution of the tree (e.g. growth rates of scale radii, baryonic component initialization) are typically handled by the ``nodeInitialize`` method.
      </description>
      <type>void</type>
      <pass>yes</pass>
      <selfTarget>yes</selfTarget>
      <argument>type(treeNode), intent(inout), target :: node</argument>
      <code>
	!$GLC attributes unused :: node
	self%methodIsNull(nodeOperatorMethodNodeTreeInitialize%ID)=.true.
      </code>
    </method>
    <method name="nodesMerge" >
      <description>
      Apply node-operator-specific physics when two dark matter halo nodes merge, for example redistributing bound mass, transferring orbital energy, or recording merger event metadata.
      </description>
      <type>void</type>
      <pass>yes</pass>
      <argument>type(treeNode), intent(inout) :: node</argument>
      <code>
	!$GLC attributes unused :: node
	self%methodIsNull(nodeOperatorMethodNodesMerge%ID)=.true.
      </code>
    </method>
    <method name="nodePromote" >
      <description>
      Apply node-operator-specific physics when a node is promoted to replace its parent in the merger tree, for example inheriting parent properties or updating mass and structural parameters.
      </description>
      <type>void</type>
      <pass>yes</pass>
      <argument>type(treeNode), intent(inout) :: node</argument>
      <code>
	!$GLC attributes unused :: node
	self%methodIsNull(nodeOperatorMethodNodePromote%ID)=.true.
      </code>
    </method>
    <method name="galaxiesMerge" >
      <description>
      Apply node-operator-specific physics when two galaxies merge within a common dark matter halo, for example triggering starbursts, transferring cold gas, or updating spheroid properties.
      </description>
      <type>void</type>
      <pass>yes</pass>
      <argument>type(treeNode), intent(inout) :: node</argument>
      <code>
	!$GLC attributes unused :: node
	self%methodIsNull(nodeOperatorMethodGalaxiesMerge%ID)=.true.
      </code>
    </method>
    <method name="differentialEvolutionPre" >
      <description>
      Apply pre-evolution setup operations on a node immediately before the differential equation solver advances its properties, for example resetting accumulators or evaluating cached quantities.
      </description>
      <type>void</type>
      <pass>yes</pass>
      <argument>type(treeNode), intent(inout) :: node</argument>
      <code>
	!$GLC attributes unused :: node
	self%methodIsNull(nodeOperatorMethodDifferentialEvolutionPre%ID)=.true.
      </code>
    </method>
    <method name="differentialEvolutionScales" >
      <description>
      Set the characteristic scale values for any meta-properties that are evolved by the ODE solver, so that relative tolerances are computed correctly during integration.
      </description>
      <type>void</type>
      <pass>yes</pass>
      <argument>type(treeNode), intent(inout) :: node</argument>
      <code>
	!$GLC attributes unused :: node
	self%methodIsNull(nodeOperatorMethodDifferentialEvolutionScales%ID)=.true.
      </code>
    </method>
    <method name="differentialEvolutionAnalytics" >
      <description>
      Mark (meta-)properties as analytically solved in differential evolution.
      </description>
      <type>void</type>
      <pass>yes</pass>
      <argument>type(treeNode), intent(inout) :: node</argument>
      <code>
	!$GLC attributes unused :: node
	self%methodIsNull(nodeOperatorMethodDifferentialEvolutionAnalytics%ID)=.true.
      </code>
    </method>
    <method name="differentialEvolutionInactives" >
      <description>
      Mark (meta-)properties as inactive in differential evolution.
      </description>
      <type>void</type>
      <pass>yes</pass>
      <argument>type(treeNode), intent(inout) :: node</argument>
      <code>
	!$GLC attributes unused :: node
	self%methodIsNull(nodeOperatorMethodDifferentialEvolutionInactives%ID)=.true.
      </code>
    </method>
    <method name="differentialEvolution" >
      <description>
      Evaluate and add derivatives for any node meta-properties that this operator evolves during the ODE integration step, and signal if an interrupt is required.
      </description>
      <type>void</type>
      <pass>yes</pass>
      <selfTarget>yes</selfTarget>
      <argument>type     (treeNode     ), intent(inout), target  :: node             </argument>
      <argument>logical                 , intent(inout)          :: interrupt        </argument>
      <argument>procedure(interruptTask), intent(inout), pointer :: functionInterrupt</argument>
      <argument>integer                 , intent(in   )          :: propertyType     </argument>
      <code>
	!$GLC attributes unused :: node, interrupt, functionInterrupt, propertyType
	self%methodIsNull(nodeOperatorMethodDifferentialEvolution%ID)=.true.
      </code>
    </method>
    <method name="differentialEvolutionSolveAnalytics" >
      <description>
      Set the values of analytically-solvable properties of a node during differential evolution.
      </description>
      <type>void</type>
      <pass>yes</pass>
      <argument>type            (treeNode), intent(inout) :: node</argument>
      <argument>double precision          , intent(in   ) :: time</argument>
      <code>
	!$GLC attributes unused :: node, time
	self%methodIsNull(nodeOperatorMethodDifferentialEvolutionSolveAnalytics%ID)=.true.
      </code>
    </method>
    <method name="predeterminedSolveAnalytics" >
      <description>
      Set the values of analytically-solvable pre-determined properties of a node.
      </description>
      <type>void</type>
      <pass>yes</pass>
      <argument>type            (treeNode), intent(inout) :: node</argument>
      <argument>double precision          , intent(in   ) :: time</argument>
      <code>
	!$GLC attributes unused :: node, time
	self%methodIsNull(nodeOperatorMethodPredeterminedSolveAnalytics%ID)=.true.
      </code>
    </method>
    <method name="differentialEvolutionStepFinalState" >
      <description>
      Apply any required corrections or updates to node properties at the completion of each ODE integration step, before the solver accepts or rejects the step.
      </description>
      <type>void</type>
      <pass>yes</pass>
      <argument>type(treeNode), intent(inout) :: node</argument>
      <code>
	!$GLC attributes unused :: node
	self%methodIsNull(nodeOperatorMethodDifferentialEvolutionStepFinalState%ID)=.true.
      </code>
    </method>
    <method name="differentialEvolutionPost" >
      <description>
      Apply post-evolution updates to a node after the ODE solver has advanced all properties to a new time, for example recording derived quantities or checking thresholds.
      </description>
      <type>void</type>
      <pass>yes</pass>
      <argument>type(treeNode), intent(inout) :: node</argument>
      <code>
	!$GLC attributes unused :: node
	self%methodIsNull(nodeOperatorMethodDifferentialEvolutionPost%ID)=.true.
      </code>
    </method>
    <method name="differentialEvolutionPostStep" >
      <description>
      Apply post-step corrections or diagnostics to a node after each accepted ODE solver step, allowing fine-grained control of properties between accepted integration steps.
      </description>
      <type>void</type>
      <pass>yes</pass>
      <argument>type   (treeNode), intent(inout) :: node</argument>
      <argument>integer          , intent(inout) :: status</argument>
      <code>
	!$GLC attributes unused :: node, status
	self%methodIsNull(nodeOperatorMethodDifferentialEvolutionPostStep%ID)=.true.
      </code>
    </method>
    <!-- Flags which methods fall through to the default (null) implementation for this object. Each default implementation sets
         its flag when called, allowing nodeOperatorMulti to skip calling it thereafter. -->
    <data scope="self">logical, dimension(nodeOperatorMethodCount) :: methodIsNull=.false.</data>
  </functionClass>
  !!]

end module Nodes_Operators
