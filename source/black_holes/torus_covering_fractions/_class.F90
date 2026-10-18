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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a module which provides a class implementing the covering fraction of the dusty torus around an AGN.
!!}

module Black_Hole_Torus_Covering_Fractions
  !!{RST
  Provides a class implementing the covering fraction of the dusty torus around an AGN.

  The covering fraction is the fraction of the sky, as seen from the accretion disk, which the torus obscures---equal,
  for randomly oriented AGN, to the fraction of them seen as Type 2. It sets the angular thickness of the torus in
  :galacticus-class:`dustAttenuationAGNTorus`, and so how much of the light of the accretion disk the torus absorbs and
  re-emits in the infrared.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  implicit none
  private

  !![
  <functionClass docformat="rst">
   <name>blackHoleTorusCoveringFraction</name>
   <descriptiveName>Black Hole Torus Covering Fractions</descriptiveName>
   <description>
   Class providing the covering fraction of the dusty torus around the central black hole of a galaxy.
   </description>
   <default>recedingTorus</default>
   <method name="coveringFraction" >
    <description>
    Return the fraction of the sky, as seen from the accretion disk of the central black hole of ``node``, which is
    covered by its dusty torus. The result lies between zero and one.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
  </functionClass>
  !!]

end module Black_Hole_Torus_Covering_Fractions
