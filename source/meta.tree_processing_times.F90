!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
  implicit none
  private

  !![
  <functionClass>
   <name>metaTreeProcessingTime</name>
   <descriptiveName>Merger Tree Processing Times</descriptiveName>
   <description>Class providing estimates of processing times for merger trees.</description>
   <default>file</default>
   <method name="time" >
    <description>Return an estimate of the time needed to process a tree of the given mass.</description>
    <type>double precision</type>
    <argument>double precision, intent(in   ) :: massTree</argument>
    <pass>yes</pass>
   </method>
  </functionClass>
  !!]

end module Meta_Tree_Compute_Times
