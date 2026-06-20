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

!!{RST
Contains a module that provides a class for dark matter halo spin distributions.
!!}

module Halo_Spin_Distributions
  !!{RST
  Provides a class for dark matter halo spin distributions.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass docformat="rst">
   <name>haloSpinDistribution</name>
   <descriptiveName>Dark Matter Halo Spin Parameter Distributions</descriptiveName>
   <description>
   Class providing the distribution :math:`p(\lambda)` of the dark matter halo spin parameter :math:`\lambda \equiv J |E|^{1/2} / \mathrm{G} M^{5/2}`, where :math:`J` is the angular momentum, :math:`E` the total energy, and :math:`M` the halo mass. The spin distribution, sometimes approximated by a log-normal, determines the range of galaxy disk sizes that form from cooling gas conserving the specific angular momentum of the halo. Implementations return both the distribution function :math:`p(\lambda)` and random samples from it, enabling Monte Carlo assignment of spins to halos.
   </description>
   <default>bett2007</default>
   <method name="sample" >
    <description>
    Samples a spin parameter from the distribution for the given ``node``.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
   <method name="distribution" >
    <description>
    Return the spin distribution function, :math:`p(\lambda)`, for the given ``node``. It is assumed that ``node`` provides the value of the spin at which the distribution function should be evaluated.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
  </functionClass>
  !!]

end module Halo_Spin_Distributions
