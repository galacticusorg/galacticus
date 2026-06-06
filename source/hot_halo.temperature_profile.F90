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
Contains a module which provides a hot halo temperature profile class.
!!}

module Hot_Halo_Temperature_Profiles
  !!{
  Provides a hot halo temperature profile class.
  !!}
  use :: Galacticus_Nodes  , only : treeNode
  use :: Mass_Distributions, only : kinematicsDistributionClass
  private

  !![
  <functionClass>
   <name>hotHaloTemperatureProfile</name>
   <descriptiveName>Hot halo temperature profiles</descriptiveName>
   <description>Class providing the temperature profile of the hot gas atmosphere surrounding a galaxy,
    returned as a \refClass{kinematicsDistributionClass} object. The temperature profile enters the local
    cooling time calculation, determines the thermal pressure support against gravity, and sets the sound
    speed relevant for ram pressure estimates. Implementations range from an isothermal profile at the virial
    temperature to polytropic or observationally-motivated radially varying profiles.</description>
   <default>virial</default>
   <method name="get" >
    <description>Return the temperature distribution of the hot halo.</description>
    <type>class(kinematicsDistributionClass)</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
  </functionClass>
  !!]

end module Hot_Halo_Temperature_Profiles
