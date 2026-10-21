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
Contains a module which provides a class that implements heating of the mass distributions of galactic components.
!!}

module Galactic_Component_Heatings
  !!{RST
  Provides a class that implements heating of the mass distributions of galactic components.
  !!}
  use :: Mass_Distributions, only : massDistributionSpherical
  private

  !![
  <functionClass docformat="rst">
   <name>galacticComponentHeating</name>
   <descriptiveName>Galactic Component Heatings</descriptiveName>
   <description>
   Class providing models of the heating of the mass distributions of spherical galactic components (spheroids and nuclear star
   clusters) by tidal shocks. Heating is applied after any truncation of the mass distribution (see
   :galacticus-class:`galacticComponentTruncationClass`), such that the truncated distribution is heated. Heating is applied
   only to satellites.
   </description>
   <default>none</default>
   <method name="isActive" >
    <description>
    Returns true if this class heats mass distributions.
    </description>
    <type>logical</type>
    <pass>yes</pass>
   </method>
   <method name="massDistribution" >
    <description>
    Returns the mass distribution ``massDistribution_`` heated consistently with the heated dark matter (only) mass distribution,
    ``massDistributionDarkMatter_``, of the host halo.
    </description>
    <type>class(massDistributionSpherical)</type>
    <pass>yes</pass>
    <argument>class(massDistributionSpherical), intent(inout), target :: massDistribution_, massDistributionDarkMatter_</argument>
   </method>
  </functionClass>
  !!]

end module Galactic_Component_Heatings
