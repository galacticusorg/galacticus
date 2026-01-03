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
Contains a module which implements lightcone geometries.
!!}

module Geometry_Lightcones
  !!{
  Implements geometries of lightcones.
  !!}
  use            :: Galacticus_Nodes, only : treeNode
  use, intrinsic :: ISO_C_Binding   , only : c_size_t
  private

  !![
  <functionClass>
   <name>geometryLightcone</name>
   <descriptiveName>Lightcone Geometries</descriptiveName>
   <description>Class providing geometries of lightcones.</description>
   <default>null</default>
   <method name="timeMinimum" >
    <description>Returns the minimum time in the lightcone.</description>
    <type>double precision</type>
    <pass>yes</pass>
   </method>
   <method name="timeMaximum" >
    <description>Returns the maximum time in the lightcone.</description>
    <type>double precision</type>
    <pass>yes</pass>
   </method>
   <method name="isInLightcone" >
    <description>Returns true if the provided node lies within the lightcone.</description>
    <type>logical</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout)           :: node</argument>
    <argument>logical                   , intent(in   ), optional :: atPresentEpoch</argument>
    <argument>double precision          , intent(in   ), optional :: radiusBuffer</argument>
   </method>
   <method name="replicationCount" >
    <description>Returns the number of times the given nodes appears in the lightcone .</description>
    <type>integer(c_size_t)</type>
    <pass>yes</pass>
    <argument>type   (treeNode), intent(inout) :: node</argument>
   </method>
   <method name="solidAngle" >
    <description>Returns the solid angle subtended by the lightcone (in units of steradians).</description>
    <type>double precision</type>
    <pass>yes</pass>
   </method>
   <method name="position" >
    <description>Returns the position vector of a {\normalfont \ttfamily node} (in units of Mpc) in the lightcone coordinate system.</description>
    <type>double precision, dimension(3)</type>
    <pass>yes</pass>
    <argument>type   (treeNode), intent(inout), target :: node</argument>
    <argument>integer(c_size_t), intent(in   )         :: instance</argument>
   </method>
   <method name="velocity" >
    <description>Returns the velocity vector of a {\normalfont \ttfamily node} (in units of km/s) in the lightcone coordinate system.</description>
    <type>double precision, dimension(3)</type>
    <pass>yes</pass>
    <argument>type   (treeNode), intent(inout) :: node</argument>
    <argument>integer(c_size_t), intent(in   ) :: instance</argument>
   </method>
   <method name="timeLightconeCrossing" >
    <description>Returns the next time in the interval from the current node time to {\normalfont \ttfamily timeEnd} at which any replicant of this node will cross the lightcone. If no crossing occurs during this interval a very large value is returned instead.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout)                                      :: node                  </argument>
    <argument>double precision          , intent(in   )                                      :: timeStart    , timeEnd</argument>
    <argument>double precision          , intent(inout), dimension(:), allocatable, optional :: timesCrossing         </argument>
   </method>
   <method name="positionLightconeCrossing" >
    <description>Returns the position of the node at the time of lightcone crossing---which must have been previously identified via the {\normalfont \ttfamily timeLightconeCrossing} method.</description>
    <type>double precision, dimension(3)</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
   <method name="velocityLightconeCrossing" >
    <description>Returns the velocity of the node at the time of lightcone crossing---which must have been previously identified via the {\normalfont \ttfamily timeLightconeCrossing} method.</description>
    <type>double precision, dimension(3)</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
  </functionClass>
  !!]

end module Geometry_Lightcones
