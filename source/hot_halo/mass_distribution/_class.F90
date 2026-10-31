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
Contains a module which provides a hot halo mass distribution class.
!!}

module Hot_Halo_Mass_Distributions
  !!{RST
  Provides an object which provides a hot halo mass distribution class.
  !!}
  use :: Galacticus_Nodes          , only : treeNode
  use :: Mass_Distributions        , only : massDistributionClass
  use :: Galactic_Structure_Options, only : enumerationWeightByType
  private

  !![
  <functionClass docformat="rst">
   <name>hotHaloMassDistribution</name>
   <descriptiveName>Hot Halo Mass Distributions</descriptiveName>
   <description>
   Class providing the radial mass distribution of hot (virialized) gas in the halo, returned as a :galacticus-class:`massDistributionClass` object. The density profile of the hot atmosphere sets the local cooling rate and pressure support, and determines the ram pressure experienced by satellite galaxies. Common profiles include the :math:`\beta`-model and hydrostatic equilibrium solutions. The distribution can be weighted by mass or by other quantities for use in different physical calculations.
   </description>
   <default>betaProfile</default>
   <method name="get" >
    <description>
    Return the mass distribution of the hot halo. By default the distribution is truncated at the outer radius of the hot halo. If ``untruncated`` is true the same distribution (with the same shape and normalization) is returned without truncation, such that it extends beyond the outer radius. In that case, if the hot halo is empty, implementations should return the distribution normalized as for a hot halo of unit mass, as such requests are used only for the shape of the profile. Implementations which never truncate the distribution need apply only the empty-halo normalization.
    </description>
    <type>class(massDistributionClass)</type>
    <pass>yes</pass>
    <argument>type   (treeNode               ), intent(inout)           :: node       </argument>
    <argument>type   (enumerationWeightByType), intent(in   ), optional :: weightBy   </argument>
    <argument>integer                         , intent(in   ), optional :: weightIndex</argument>
    <argument>logical                         , intent(in   ), optional :: untruncated</argument>
   </method>
  </functionClass>
  !!]

end module Hot_Halo_Mass_Distributions
