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
Contains a module that provides a class for dark matter halo spin distributions.
!!}

module Halo_Spin_Distributions
  !!{
  Provides a class for dark matter halo spin distributions.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>haloSpinDistribution</name>
   <descriptiveName>Dark Matter Halo Spin Parameter Distributions</descriptiveName>
   <description>
    Class providing dark matter halo spin parameter distributions.
   </description>
   <default>bett2007</default>
   <method name="sample" >
    <description>Samples a spin parameter from the distribution for the given {\normalfont \ttfamily node}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
   <method name="distribution" >
    <description>Return the spin distribution function, $p(\lambda)$, for the given {\normalfont \ttfamily node}. It is assumed that {\normalfont \ttfamily node} provides the value of the spin at which the distribution function should be evaluated.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
  </functionClass>
  !!]

end module Halo_Spin_Distributions
