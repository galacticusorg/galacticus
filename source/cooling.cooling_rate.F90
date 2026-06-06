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
Contains a module that implements calculations of the cooling rate.
!!}

module Cooling_Rates
  !!{
  Provides a class that implements calculations of the cooling rate.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>coolingRate</name>
   <descriptiveName>Cooling Rates</descriptiveName>
   <description>Class providing models of the mass cooling rate (in $\mathrm{M}_\odot$ Gyr$^{-1}$) at which gas cools out
    of the hot atmosphere of a dark matter halo and becomes available for accretion onto the central galaxy. The
    cooling rate drives the supply of cold gas for star formation and sets the growth rate of the galaxy disc.
    Implementations typically account for the interplay between the cooling time, the freefall time, the
    cooling radius, and the available mass of hot gas, following the two-regime picture of \glc{whiteFrenk1991}
    and subsequent refinements.</description>
   <default>whiteFrenk1991</default>
   <method name="rate" >
    <description>Returns the cooling rate of gas in the hot atmosphere surrounding the galaxy in \mono{node} in units of $\mathrm{M}_\odot/$Gyr.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
  </functionClass>
  !!]

end module Cooling_Rates
