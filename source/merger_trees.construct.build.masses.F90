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
Contains a module which implements a class for creating sets of tree masses to use when building merger trees.
!!}

module Merger_Trees_Build_Masses
  !!{
  Implements a class for creating sets of tree masses to use when building merger trees.
  !!}
  private

  !![
  <functionClass>
   <name>mergerTreeBuildMasses</name>
   <descriptiveName>Merger Tree Build Masses</descriptiveName>
   <description>Class providing methods for creating sets of tree masses to use when building merger trees.</description>
   <default>sampledDistributionUniform</default>
   <method name="construct" >
    <description>Returns a set of merger tree masses (and either their weights, or corresponding mass interval) to be built.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   )                            :: time                                  </argument>
    <argument>double precision, intent(  out), allocatable, dimension(:) :: mass, massMinimum, massMaximum, weight</argument>
   </method>
  </functionClass>
  !!]

end module Merger_Trees_Build_Masses
