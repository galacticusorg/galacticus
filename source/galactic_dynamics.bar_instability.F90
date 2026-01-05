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
Contains a module which provides a class that implements calculations of bar instability in galactic disks.
!!}

module Galactic_Dynamics_Bar_Instabilities
  !!{
  Provides a class that implements calculations of bar instability in galactic disks.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>galacticDynamicsBarInstability</name>
   <descriptiveName>Bar instabilities in galactic disks</descriptiveName>
   <description>
    Class providing models of bar instabilities in galactic disks.
   </description>
   <default>efstathiou1982</default>
   <method name="timescale" >
    <description>
     Returns a timescale on which the bar instability depletes material from a disk into a pseudo-bulge. A negative value
     indicates no instability. Also returns the net torque due to any external force causing this instability, and the fraction of
     the angular momentum of the material depleted into the pseudo-bulge which is retained by the disk and the spheroid.
     </description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout) :: node</argument>
    <argument>double precision          , intent(  out) :: timescale, externalDrivingSpecificTorque, fractionAngularMomentumRetainedDisk, fractionAngularMomentumRetainedSpheroid</argument>
   </method>
  </functionClass>
  !!]

end module Galactic_Dynamics_Bar_Instabilities
