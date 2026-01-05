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
Contains a module that implements a class for calculations of ram pressure stripping timescales for hot halos.
!!}

module Hot_Halo_Ram_Pressure_Stripping_Timescales
  !!{
  Implements a class for calculations of ram pressure stripping timescales for hot halos.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>hotHaloRamPressureTimescale</name>
   <descriptiveName>Models of ram pressure stripping timescales due to the hot halo.</descriptiveName>
   <description>
    Class providing models of ram pressure stripping timescales due to the hot halo.
   </description>
   <default>ramPressureAcceleration</default>
   <method name="timescale" >
    <description>Return the ram pressure stripping timescale for {\normalfont \ttfamily node} (in units of Gyr).</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
  </functionClass>
  !!]

end module Hot_Halo_Ram_Pressure_Stripping_Timescales
