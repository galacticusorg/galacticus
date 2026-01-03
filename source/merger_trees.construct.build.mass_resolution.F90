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
Contains a module which provides a class of merger tree mass resolutions.
!!}

module Merger_Trees_Build_Mass_Resolution
  !!{
  Provides a class of merger tree mass resolutions.
  !!}
  use :: Galacticus_Nodes, only : mergerTree
  private

  !![
  <functionClass>
   <name>mergerTreeMassResolution</name>
   <descriptiveName>Merger Tree Building Mass Resolutions</descriptiveName>
   <description>
    Class providing mass resolutions to use when building merger trees.
   </description>
   <default>fixed</default>
   <method name="resolution" >
    <description>Gives the mass resolution to use for the given tree.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(mergerTree), intent(in   ) :: tree</argument>
   </method>
  </functionClass>
  !!]

end module Merger_Trees_Build_Mass_Resolution
