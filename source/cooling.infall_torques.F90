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
Contains a module that implements calculations of the infall torques for cooling calculations.
!!}

module Cooling_Infall_Torques
  !!{
  Provides a class that implements calculations of the infall torques for cooling calculations.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  implicit none
  private

  !![
  <functionClass>
   <name>coolingInfallTorque</name>
   <descriptiveName>Cooling Infall Torque</descriptiveName>
   <description>Class providing models of the fraction of angular momentum lost by gas as it falls from the hot
    halo onto the galaxy. Torques from the dark matter halo, the existing galaxy, and gas dynamical effects can
    cause infalling gas to lose angular momentum before it joins the disk, thereby producing a more compact disc
    than would result from purely angular-momentum-conserving infall. The returned fraction of angular momentum
    loss directly modulates the specific angular momentum assigned to the cooling gas.</description>
   <default>fixed</default>
   <method name="fractionAngularMomentumLoss" >
    <description>Return the fraction of specific angular momentum lost by infalling gas as it travels from the infall radius to the galaxy disk, where torques from the halo, galaxy, or gas dynamics reduce the angular momentum of the accreting gas.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
  </functionClass>
  !!]

end module Cooling_Infall_Torques
