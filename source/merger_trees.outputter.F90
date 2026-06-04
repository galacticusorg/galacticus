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
Contains a module which provides a class that implements outputting of node data in a merger tree.
!!}

module Merger_Tree_Outputters
  !!{
  Provides a class that implements evolution of merger trees.
  !!}
  use            :: Galacticus_Nodes, only : mergerTree
  use, intrinsic :: ISO_C_Binding   , only : c_size_t
  private

  !![
  <enumeration>
   <name>outputGroupType</name>
   <description>Enumeration of the types of output groups.</description>
   <decodeFunction>yes</decodeFunction>
   <visibility>public</visibility>
   <entry label="tree"     />
   <entry label="node"     />
   <entry label="snapshot" />
   <entry label="lightcone"/>
  </enumeration>
  !!]

  !![
  <functionClass>
   <name>mergerTreeOutputter</name>
   <descriptiveName>Merger Tree Outputters</descriptiveName>
   <description>Class providing outputters for merger trees---objects that write merger tree and galaxy
    data to persistent storage at each requested output time. An outputter traverses the nodes of a
    merger tree and serializes the desired set of galaxy/halo properties to HDF5 or other formats.
    It also handles reduction across MPI processes (combining partial outputs from different CPU
    ranks) and finalization at the end of the simulation. The standard outputter writes node
    properties as defined by the active \refClass{nodePropertyExtractorClass} instances.</description>
   <default>standard</default>
   <method name="outputTree" >
    <description>Serialize all galaxy and halo properties from the given \mono{tree} to the output dataset corresponding to output index \mono{indexOutput} at cosmic time \mono{time}. Optionally provide an \mono{outputType} describing the type of output.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type            (mergerTree                    ), intent(inout), target   :: tree       </argument>
    <argument>integer         (c_size_t                      ), intent(in   )           :: indexOutput</argument>
    <argument>double precision                                , intent(in   )           :: time       </argument>
    <argument>type            (enumerationOutputGroupTypeType), intent(in   ), optional :: outputType </argument>
   </method>
   <method name="outputNode" >
    <description>Serialize galaxy and halo properties of the given individual \mono{node} to the output dataset corresponding to output index \mono{indexOutput}.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type   (treeNode                      ), intent(inout)           :: node       </argument>
    <argument>integer(c_size_t                      ), intent(in   )           :: indexOutput</argument>
    <argument>type   (enumerationOutputGroupTypeType), intent(in   ), optional :: outputType </argument>
   </method>
   <method name="finalize" >
    <description>Finalize the output of merger trees, flushing any buffered data to persistent storage and performing any post-processing required after all trees have been serialized. Optionally provide an \mono{outputType} describing the type of output.</description>
    <type>void</type>
    <pass>yes</pass>
    <code>
     !$GLC attributes unused :: self
    </code>
   </method>
   <method name="reduce" >
    <description>Reduce the object onto another object of the class.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>class(mergerTreeOutputterClass), intent(inout) :: reduced</argument>
    <code>
     !$GLC attributes unused :: self, reduced
    </code>
   </method>
  </functionClass>
  !!]

end module Merger_Tree_Outputters
