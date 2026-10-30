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
Contains a module which provides a class that implements tidal stripping.
!!}

module Tidal_Stripping_Mass_Loss_Rate
  !!{RST
  Provides a class that implements calculations of tidal stripping.
  !!}
  use :: Galacticus_Nodes, only : nodeComponent
  private

  !![
  <functionClass docformat="rst">
   <name>tidalStripping</name>
   <descriptiveName>Tidal Stripping Mass Loss Rates</descriptiveName>
   <description>
   Class providing models of tidal stripping of galactic components (the gas and stars of disks, spheroids, and nuclear star clusters)---the removal of material from a satellite galaxy as it orbits within its host halo's tidal field. Unlike dark matter tidal stripping (handled by :galacticus-class:`satelliteTidalStrippingClass`), this class operates on individual galaxy components and returns the mass loss rate (in :math:`\mathrm{M}_\odot` Gyr\ :math:`^{-1}`) for the given component. Optionally, it also returns the specific angular momentum of the stripped material in units of the mean specific angular momentum of the component.
   </description>
   <default>simple</default>
   <method name="rateMassLoss" >
    <description>
    Returns the rate of mass loss (in :math:`\mathrm{M}_\odot` Gyr\ :math:`^{-1}`) due to tidal stripping of ``component``. If present, ``ratioSpecificAngularMomentum`` is set to the specific angular momentum of the stripped material in units of the mean specific angular momentum of ``component``.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>class           (nodeComponent), intent(inout)           :: component                   </argument>
    <argument>double precision               , intent(  out), optional :: ratioSpecificAngularMomentum</argument>
   </method>
  </functionClass>
  !!]

end module Tidal_Stripping_Mass_Loss_Rate
