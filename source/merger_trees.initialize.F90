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
Contains a module which provides a class that implements initialization of merger trees.
!!}

module Merger_Tree_Initialization
  !!{
  Provides a class that implements initialization of merger trees.
  !!}
  use :: Galacticus_Nodes, only : mergerTree
  private

  !![
  <functionClass>
   <name>mergerTreeInitializor</name>
   <descriptiveName>Merger Tree Initializors</descriptiveName>
   <description>Class providing initialization of merger trees.</description>
   <default>standard</default>
   <method name="initialize" >
    <description>Initialize the given tree.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type            (mergerTree), intent(inout) :: tree   </argument>
    <argument>double precision            , intent(in   ) :: timeEnd</argument>
   </method>
  </functionClass>
  !!]

end module Merger_Tree_Initialization
