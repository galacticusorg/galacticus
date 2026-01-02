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
Contains a module which provides a class that implements galactic filters.
!!}

module Galactic_Filters
  !!{
  Provides an object that implements galactic filters.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>galacticFilter</name>
   <descriptiveName>Galactic Filter</descriptiveName>
   <description>Object providing boolean filters acting on galaxies.</description>
   <default>always</default>
   <method name="passes" >
    <description>Return true if the given {\normalfont \ttfamily node} passes the filter.</description>
    <type>logical</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout), target :: node</argument>
   </method>
  </functionClass>
  !!]

  type, public :: filterList
     class(galacticFilterClass), pointer :: filter_ => null()
     type (filterList         ), pointer :: next    => null()
  end type filterList

end module Galactic_Filters
