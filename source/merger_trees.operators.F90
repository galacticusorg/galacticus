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

module Merger_Tree_Operators
  !!{
  Provides an object that implements operators acting on merger trees.
  !!}
  use :: Galacticus_Nodes, only : mergerTree
  private

  !![
  <functionClass>
   <name>mergerTreeOperator</name>
   <descriptiveName>Merger Tree Operators</descriptiveName>
   <description>
    Class providing operators acting on merger trees.
   </description>
   <default>null</default>
   <method name="operatePreConstruction" >
    <description>Perform an operation on the merger tree prior to construction.</description>
    <type>void</type>
    <pass>yes</pass>
    <code>
     !$GLC attributes unused :: self
     ! Nothing to do.
     return
    </code>
   </method>
   <method name="operatePreInitialization" >
    <description>Perform an operation on the merger tree prior to initialization.</description>
    <type>void</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>type(mergerTree), intent(inout), target :: tree</argument>
    <code>
     !$GLC attributes unused :: self, tree
     ! Nothing to do.
     return
    </code>
   </method>
   <method name="operatePreEvolution" >
    <description>Perform an operation on the merger tree prior to evolution.</description>
    <type>void</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>type(mergerTree), intent(inout), target :: tree</argument>
    <code>
     !$GLC attributes unused :: self, tree
     ! Nothing to do.
     return
    </code>
   </method>
   <method name="operatePostEvolution" >
    <description>Perform an operation on the merger tree after evolution.</description>
    <type>void</type>
    <pass>yes</pass>
    <code>
     !$GLC attributes unused :: self
     ! Nothing to do.
     return
    </code>
   </method>
   <method name="finalize" >
    <description>Finalize a merger tree operator.</description>
    <type>void</type>
    <pass>yes</pass>
    <code>
     !$GLC attributes unused :: self
     ! Nothing to do.
     return
    </code>
   </method>
  </functionClass>
  !!]

end module Merger_Tree_Operators
