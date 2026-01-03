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
Contains a module which set random number seeds in merger trees.
!!}

module Merger_Tree_Seeds
  !!{
  Implements a class for setting random number seeds in merger trees.
  !!}
  use :: Galacticus_Nodes, only : mergerTree
  private

  !![
  <functionClass>
   <name>mergerTreeSeeds</name>
   <descriptiveName>Merger Tree Random Number Seeds</descriptiveName>
   <description>Class implementing setting of random number seeds in merger trees.</description>
   <default>incremental</default>
   <method name="set" >
    <description>Set the random number seed in the given {\normalfont \ttfamily tree}.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type(mergerTree), intent(inout) :: tree</argument>
   </method>
  </functionClass>
  !!]

end module Merger_Tree_Seeds
