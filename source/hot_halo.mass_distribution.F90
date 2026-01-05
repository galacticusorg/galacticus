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
Contains a module which provides a hot halo mass distribution class.
!!}

module Hot_Halo_Mass_Distributions
  !!{
  Provides an object which provides a hot halo mass distribution class.
  !!}
  use :: Galacticus_Nodes          , only : treeNode
  use :: Mass_Distributions        , only : massDistributionClass
  use :: Galactic_Structure_Options, only : enumerationWeightByType
  private

  !![
  <functionClass>
   <name>hotHaloMassDistribution</name>
   <descriptiveName>Hot Halo Mass Distributions</descriptiveName>
   <description>
    Object implementing hot halo mass distributions.
   </description>
   <default>betaProfile</default>
   <method name="get" >
    <description>Return the mass distribution of the hot halo.</description>
    <type>class(massDistributionClass)</type>
    <pass>yes</pass>
    <argument>type   (treeNode               ), intent(inout)           :: node       </argument>
    <argument>type   (enumerationWeightByType), intent(in   ), optional :: weightBy   </argument>
    <argument>integer                         , intent(in   ), optional :: weightIndex</argument>
   </method>
  </functionClass>
  !!]

end module Hot_Halo_Mass_Distributions
