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
Contains a module which provides a class that implements concentrations of dark matter halo profiles.
!!}

module Dark_Matter_Profiles_Shape
  !!{
  Provides a class that implements shape parameters of dark matter halo profiles.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>darkMatterProfileShape</name>
   <descriptiveName>Dark Matter Profile Shapes</descriptiveName>
   <description>
    Class providing the shape parameter $\alpha$ of dark matter halo density profiles such as the Einasto profile,
    $\rho(r) \propto \exp\{-\frac{2}{\alpha}[(r/r_\mathrm{s})^\alpha - 1]\}$. The shape parameter controls the curvature
    of the profile near the center and is typically fit from $N$-body simulations as a function of halo mass and redshift.
   </description>
   <default>gao2008</default>
   <method name="shape" >
    <description>Returns the dimensionless shape parameter $\alpha$ of the Einasto dark matter density profile for the halo in \mono{node}, controlling the curvature of the inner profile with typical values in the range $0.1$--$0.3$ from N-body simulations.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
  </functionClass>
  !!]

end module Dark_Matter_Profiles_Shape
