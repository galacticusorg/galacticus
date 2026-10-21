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
Contains a module which provides a class that implements truncation of the mass distributions of galactic components.
!!}

module Galactic_Component_Truncations
  !!{RST
  Provides a class that implements truncation of the mass distributions of galactic components.
  !!}
  use :: Mass_Distributions, only : massDistributionSpherical
  private

  !![
  <functionClass docformat="rst">
   <name>galacticComponentTruncation</name>
   <descriptiveName>Galactic Component Truncations</descriptiveName>
   <description>
   Class providing models of the truncation of the mass distributions of spherical galactic components (spheroids and nuclear
   star clusters) by tidal stripping.

   The degree of truncation of a component is described by the fraction, :math:`f`, of the mass of its original, untruncated,
   mass distribution which it retains. This is tracked by the ``tidalTruncation`` meta-property of the component,
   :math:`x = -\log f`, which is created by the tidal mass loss node operators (e.g.
   :galacticus-class:`nodeOperatorTidalMassLossSpheroids`), and which evolves as :math:`\dot{x} = \dot{M}_\mathrm{strip}/M`,
   where :math:`\dot{M}_\mathrm{strip}` is the rate of tidal mass loss and :math:`M` is the mass of the component. Since only
   tidal stripping changes :math:`x`, all other processes (star formation, outflows, mergers, etc.) leave the shape of the
   mass distribution unchanged. The mass distribution of a truncated component is its dimensionless, untruncated mass
   distribution, truncated such that it retains a fraction :math:`f` of its mass, and then scaled to the current mass of the
   component.
   </description>
   <default>none</default>
   <method name="isActive" >
    <description>
    Returns true if this class truncates mass distributions.
    </description>
    <type>logical</type>
    <pass>yes</pass>
   </method>
   <method name="massDistribution" >
    <description>
    Returns the dimensionless mass distribution ``massDistribution_`` truncated such that it retains a fraction
    ``fractionMassRetained`` of its total mass.
    </description>
    <type>class(massDistributionSpherical)</type>
    <pass>yes</pass>
    <argument>class           (massDistributionSpherical), intent(inout), target :: massDistribution_   </argument>
    <argument>double precision                           , intent(in   )         :: fractionMassRetained</argument>
   </method>
  </functionClass>
  !!]

end module Galactic_Component_Truncations
