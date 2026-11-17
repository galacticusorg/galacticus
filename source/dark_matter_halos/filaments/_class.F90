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
Contains a module which provides a class that implements the cosmic web filaments feeding dark matter halos.
!!}

module Dark_Matter_Halo_Filaments
  !!{RST
  Provides a class that implements the cosmic web filaments feeding dark matter halos.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass docformat="rst">
   <name>darkMatterHaloFilaments</name>
   <descriptiveName>Dark Matter Halo Filaments</descriptiveName>
   <description>
   Class providing the directions of the cosmic web filaments along which dark matter halos accrete. Each filament is a ray
   leaving the halo, so that a single filament passing through the halo is represented by two opposite directions. The halo has
   no orientation of its own, so the directions are expressed in the same arbitrary frame as the angles, :math:`(\theta,\phi)`,
   of the virial orbits of its satellites. They are used by :galacticus-class:`virialOrbitFilamentary`.

   Filaments persist: every node on a branch of the merger tree shares the filaments of the branch, so that a halo accretes along
   the same filaments throughout its history. Implementations store the directions of the filaments of a branch in a rank-1
   ``float`` meta-property, ``filamentAxes``, of the ``basic`` component of every node of the branch, as the three direction
   cosines of each filament in turn. Whichever object first requires the filaments of a branch draws them, and all others use
   those. So, when :galacticus-class:`nodeOperatorHaloFilaments` is used, a ``darkMatterHaloFilaments`` parameter should be
   given at the top level of the parameter file, where that operator and :galacticus-class:`virialOrbitFilamentary` share it.
   </description>
   <default>randomSequential</default>
   <method name="axes" >
    <description>
    Return the directions, as an array of shape :math:`(3,N)`, of the :math:`N` filaments feeding the halo of ``node``.
    </description>
    <type>double precision, allocatable, dimension(:,:)</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout), target :: node</argument>
   </method>
  </functionClass>
  !!]

end module Dark_Matter_Halo_Filaments
