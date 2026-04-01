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
Contains a module which provides an object that implements concentrations of dark matter halo profiles.
!!}

module Dark_Matter_Profile_Scales
  !!{
  Provides a class that implements scale radii dark matter halo profiles.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  implicit none
  private

  !![
  <functionClass>
   <name>darkMatterProfileScaleRadius</name>
   <descriptiveName>Dark Matter Profile Scale Radii</descriptiveName>
   <description>Class providing the scale radius $r_\mathrm{s}$ of dark matter halo density profiles. The scale radius
    sets the characteristic transition scale between the inner and outer slopes of the density profile (e.g. at $r_\mathrm{s}$
    the NFW profile transitions from $\rho \propto r^{-1}$ to $\rho \propto r^{-3}$). Implementations may derive the
    scale radius from a concentration parameter, from energy conservation arguments, or from other empirical relations.</description>
   <default>concentration</default>
   <method name="radius" >
    <description>Returns the scale radius $r_\mathrm{s}$ (in Mpc) of the dark matter halo density profile for the halo in \mono{node}, the characteristic radius at which the logarithmic slope of the density profile transitions between inner and outer power laws.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>type(treeNode), intent(inout), target :: node</argument>
   </method>
  </functionClass>
  !!]

end module Dark_Matter_Profile_Scales
