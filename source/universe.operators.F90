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
Contains a module which provides a class that implements operators on universes.
!!}

module Universe_Operators
  !!{
  Provides a class that implements operators on universes.
  !!}
  use :: Galacticus_Nodes, only : universe
  private

  !![
  <functionClass>
   <name>universeOperator</name>
   <descriptiveName>Universe Operators</descriptiveName>
   <description>Class providing operators that act on \mono{universe} objects---top-level transformations applied
    to a fully-evolved \glc\ universe (containing the complete set of merger trees and their galaxies)
    prior to or after model evolution. Universe operators can perform global post-processing steps such as
    computing derived statistics, writing supplementary output, or modifying global properties. The interface
    provides a single \mono{operate} method that receives the universe object and performs whatever
    transformation the implementation requires.</description>
   <default>identity</default>
   <method name="operate" >
    <description>Operate on the universe.</description>
    <type>void</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>type(universe), intent(inout) :: universe_</argument>
   </method>
  </functionClass>
  !!]

end module Universe_Operators
