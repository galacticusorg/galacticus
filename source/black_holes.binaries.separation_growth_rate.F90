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

!+    Contributions to this file made by:  Stéphane Mangeon, Andrew Benson.

!!{
Contains a module which implements a class for calculations of black hole binary separation growth rate.
!!}

module Black_Hole_Binary_Separations
  !!{
  Implements a class for calculations of black hole binary separation growth rate.
  !!}
  use :: Galacticus_Nodes, only : nodeComponentBlackHole
  implicit none
  private

  !![
  <functionClass>
   <name>blackHoleBinarySeparationGrowthRate</name>
   <descriptiveName>Black Hole Binaries Separation Growth Rate</descriptiveName>
   <description>
    Class providing models of black hole binary separation growth rates.
   </description>
   <default>zero</default>
   <method name="growthRate" >
    <description>Computes the rate of growth of the separation of the given black hole and its binary companion in units of Mpc/Gyr.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>class(nodeComponentBlackHole), intent(inout) :: blackHole</argument>
   </method>
  </functionClass>
  !!]

end module Black_Hole_Binary_Separations
