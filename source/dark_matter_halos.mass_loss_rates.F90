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
Contains a module which implements calculations of mass loss rates from dark matter halos.
!!}

module Dark_Matter_Halos_Mass_Loss_Rates
  !!{
  Implements calculations of mass loss rates from dark matter halos.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  implicit none
  private

  !![
  <functionClass>
   <name>darkMatterHaloMassLossRate</name>
   <descriptiveName>Dark Matter Halo Mass Loss Rates</descriptiveName>
   <description>
    Class providing models of the rate of mass loss from dark matter (sub)halos due to tidal stripping and other
    processes. Returns the rate of change of halo mass (in $\mathrm{M}_\odot$ Gyr$^{-1}$) for a given node. This class is
    used in combination with tidal stripping models to track the evolution of subhalo masses as they orbit within a
    host halo.
   </description>
   <default>zero</default>
   <method name="rate" >
    <description>Returns the rate of mass loss (in $\mathrm{M}_\odot$/Gyr) from \mono{node}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
  </functionClass>
  !!]

end module Dark_Matter_Halos_Mass_Loss_Rates
