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

  !+    Contributions to this file made by: Matías Liempi

!!{
Contains a module which provides a class that implements rates of star formation in nuclear star clusters.
!!}

module Star_Formation_Rates_Nuclear_Star_Clusters
  !!{
  Provides a class that implements calculations of rates of formation in nuclear star clusters.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>starFormationRateNuclearStarClusters</name>
   <descriptiveName>Rates for star formation in nuclear star clusters</descriptiveName>
   <description>Class providing models of rates of star formation in nuclear star clusters.</description>
   <default>krumholz2009</default>
   <method name="rate" >
    <description>Returns the rate (in units of $\mathrm{M}_\odot$ Gyr$^{-1}$) for star formation in the nuclear star cluster component of {\normalfont \ttfamily node}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>type(treeNode)  , intent(inout) :: node</argument>
   </method>
  </functionClass>
  !!]

end module Star_Formation_Rates_Nuclear_Star_Clusters
