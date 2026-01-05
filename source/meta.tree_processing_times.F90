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
Contains a module providing a class for calculations of the time taken to process merger trees.
!!}

module Meta_Tree_Compute_Times
  !!{
  Provides a class for calculations of the time taken to process merger trees.
  !!}
  use :: Galacticus_Nodes, only : mergerTree
  implicit none
  private

  !![
  <functionClass>
   <name>metaTreeProcessingTime</name>
   <descriptiveName>Merger Tree Processing Times</descriptiveName>
   <description>Class providing estimates of processing times for merger trees.</description>
   <default>null</default>
   <method name="time" >
    <description>Return an estimate of the time needed to process a tree of the given mass.</description>
    <type>double precision</type>
    <argument>double precision, intent(in   ) :: massTree</argument>
    <pass>yes</pass>
    <code>
      ! Return a negative value indicating that no estimate is available.
      !$GLC attributes unused :: self, massTree
      metaTreeProcessingTimeTime=-1.0d0
    </code>
   </method>
   <method name="timeRemaining" >
    <description>Return an estimate of the remaining time needed to process a given tree.</description>
    <type>double precision</type>
    <argument>type(mergerTree), intent(inout) :: tree     </argument>
    <argument>double precision, intent(in   ) :: timeFinal</argument>
    <pass>yes</pass>
    <code>
      ! Return a negative value indicating that no estimate is available.
      !$GLC attributes unused :: self, tree, timeFinal
      metaTreeProcessingTimeTimeRemaining=-1.0d0
    </code>
   </method>
  </functionClass>
  !!]

end module Meta_Tree_Compute_Times
