!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

!% Contains a module which provides a class that implements physical processes.

module Nodes_Operators
  !% Provides a class that implements physical processes.
  use :: Galacticus_Nodes, only : treeNode, interruptTask
  private

  !# <functionClass>
  !#  <name>nodeOperator</name>
  !#  <descriptiveName>Node Operators</descriptiveName>
  !#  <description>Class providing operators acting on nodes.</description>
  !#  <default>null</default>
  !#  <method name="nodePromote" >
  !#   <description>Act on the promotion of a node to its parent.</description>
  !#   <type>void</type>
  !#   <pass>yes</pass>
  !#   <argument>type(treeNode), intent(inout) :: node</argument>
  !#   <code>
  !#    !$GLC attributes unused :: self, node
  !#   </code>
  !#  </method>
  !#  <method name="galaxiesMerge" >
  !#   <description>Act on the merging of two galaxies.</description>
  !#   <type>void</type>
  !#   <pass>yes</pass>
  !#   <argument>type(treeNode), intent(inout) :: node</argument>
  !#   <code>
  !#    !$GLC attributes unused :: self, node
  !#   </code>
  !#  </method>
  !#  <method name="differentialEvolutionPre" >
  !#   <description>Operate on a node before differential evolution.</description>
  !#   <type>void</type>
  !#   <pass>yes</pass>
  !#   <argument>type(treeNode), intent(inout) :: node</argument>
  !#   <code>
  !#    !$GLC attributes unused :: self, node
  !#   </code>
  !#  </method>
  !#  <method name="differentialEvolution" >
  !#   <description>Operate on a node during differential evolution.</description>
  !#   <type>void</type>
  !#   <pass>yes</pass>
  !#   <argument>type     (treeNode     ), intent(inout)          :: node             </argument>
  !#   <argument>logical                 , intent(inout)          :: interrupt        </argument>
  !#   <argument>procedure(interruptTask), intent(inout), pointer :: functionInterrupt</argument>
  !#   <argument>integer                 , intent(in   )          :: propertyType     </argument>
  !#   <code>
  !#    !$GLC attributes unused :: self, node, interrupt, functionInterrupt, propertyType
  !#   </code>
  !#  </method>
  !#  <method name="differentialEvolutionStepFinalState" >
  !#   <description>Operate on a node at the end of an ODE step.</description>
  !#   <type>void</type>
  !#   <pass>yes</pass>
  !#   <argument>type(treeNode), intent(inout) :: node</argument>
  !#   <code>
  !#    !$GLC attributes unused :: self, node
  !#   </code>
  !#  </method>
  !#  <method name="differentialEvolutionPost" >
  !#   <description>Operate on a node after differential evolution.</description>
  !#   <type>void</type>
  !#   <pass>yes</pass>
  !#   <argument>type(treeNode), intent(inout) :: node</argument>
  !#   <code>
  !#    !$GLC attributes unused :: self, node
  !#   </code>
  !#  </method>
  !# </functionClass>
  
end module Nodes_Operators
