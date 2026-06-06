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
Contains a module which implements a class for calculations of merger remnant sizes.
!!}

module Satellite_Merging_Remnant_Sizes
  !!{
  Implements a class for calculations of merger remnant sizes.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  implicit none
  private

  !![
  <functionClass>
   <name>mergerRemnantSize</name>
   <descriptiveName>Merger Remnant Sizes</descriptiveName>
   <description>Class providing models of merger remnant sizes---the half-mass radius, circular velocity, and
    specific angular momentum of the spheroid formed when two galaxies merge. Energy and angular momentum
    conservation arguments relate the remnant size to the progenitor radii, masses, and orbital parameters.
    These quantities set the structural properties of merger-built bulges and determine the subsequent
    evolution of the stellar velocity dispersion and black hole mass in the remnant galaxy.</description>
   <default>covington2008</default>
   <method name="get" >
    <description>Determine the half-mass radius, circular velocity, and specific angular momentum of the spheroidal remnant formed when two galaxies merge, using energy and angular momentum conservation arguments applied to the progenitor masses, radii, and orbital parameters.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type            (treeNode), intent(inout) :: node                                             </argument>
    <argument>double precision          , intent(  out) :: radius, velocityCircular, angularMomentumSpecific</argument>
   </method>
  </functionClass>
  !!]

  ! Value indicating that there was no change in the remnant spheroid size.
  double precision, parameter, public :: remnantNoChange=-1.0d0

end module Satellite_Merging_Remnant_Sizes
