!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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
  !% Provides an object that implements operators acting on merger trees.
  use Galacticus_Nodes
  
  !# <functionClass>
  !#  <name>mergerTreeOperator</name>
  !#  <descriptiveName>Merger Tree Operators</descriptiveName>
  !#  <description>Class providing operators acting on merger trees.</description>
  !#  <default>null</default>
  !#  <method name="operate" >
  !#   <description>Perform an operation on the merger tree.</description>
  !#   <type>void</type>
  !#   <pass>yes</pass>
  !#   <argument>type(mergerTree), intent(inout), target :: tree</argument>
  !#  </method>
  !#  <method name="finalize" >
  !#   <description>Finalize a merger tree operator.</description>
  !#   <type>void</type>
  !#   <pass>yes</pass>
  !#   <code>
  !#    !GCC$ attributes unused :: self
  !#    ! Nothing to do.
  !#    return
  !#   </code>
  !#  </method>
  !# </functionClass>

end module Merger_Tree_Operators
