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
Contains a module that implements a class for calculations of ram pressure stripping of hot halos.
!!}

module Hot_Halo_Ram_Pressure_Stripping
  !!{
  Implements a class for calculations of ram pressure stripping of hot halos.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>hotHaloRamPressureStripping</name>
   <descriptiveName>Models of ram pressure stripping due to the hot halo.</descriptiveName>
   <description>
    Class providing models for the radius to which the hot halo is stripped by ram pressure forces.
   </description>
   <default>font2008</default>
   <method name="radiusStripped" >
    <description>Return the radius to which {\normalfont \ttfamily node} is stripped due to ram pressure from its host halo (in units of Mpc).</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>type(treeNode), intent(inout), target :: node</argument>
   </method>
  </functionClass>
  !!]

end module Hot_Halo_Ram_Pressure_Stripping
