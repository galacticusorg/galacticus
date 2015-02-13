!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which provides a class of merger tree builders.

module Merger_Trees_Builders
  !% Provides a class of merger tree builders.
  use ISO_Varying_String
  use Galacticus_Nodes
  use FGSL
  !# <include directive="mergerTreeBuilder" type="functionModules" >
  include 'mergerTreeBuilder.functionModules.inc'
  !# </include>
  private

  !# <include directive="mergerTreeBuilder" type="function" >
  !#  <descriptiveName>Merger Tree Builders</descriptiveName>
  !#  <description>Class providing merger tree builders.</description>
  !#  <default>cole2000</default>
  !#  <defaultThreadPrivate>yes</defaultThreadPrivate>
  !#  <stateful>yes</stateful>
  !#  <method name="build" >
  !#   <description>Builds and returns a merger tree given the root {\normalfont \ttfamily node}.</description>
  !#   <type>void</type>
  !#   <pass>yes</pass>
  !#   <argument>type(mergerTree), intent(inout), target :: tree</argument>
  !#  </method>
  include 'mergerTreeBuilder.type.inc'
  !# </include>
  
end module Merger_Trees_Builders
