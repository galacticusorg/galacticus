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
Contains a module that implements the dimensionless drag force due to dynamical friction in a gaseous medium.
!!}

module Satellite_Dynamical_Friction_Gaseous
  !!{RST
  Implements the dimensionless drag force due to dynamical friction in a gaseous medium.
  !!}
  private

  !![
  <functionClass docformat="rst">
   <name>dynamicalFrictionGaseous</name>
   <descriptiveName>Gaseous Dynamical Friction</descriptiveName>
   <description>
   Class providing the dimensionless drag force, :math:`\mathcal{I}`, experienced by a perturber moving through a gaseous
   (collisional) medium. The drag force on a perturber of mass :math:`M` moving at speed :math:`V` through gas of density
   :math:`\rho` is

   .. math::

      \mathbf{F} = -\frac{4\pi \mathrm{G}^2 M^2 \rho}{V^2} \mathcal{I}(\mathcal{M},\ln\Lambda) \hat{\mathbf{v}},

   where :math:`\mathcal{M}=V/c_\mathrm{s}` is the Mach number of the perturber relative to the sound speed,
   :math:`c_\mathrm{s}`, of the gas, and :math:`\ln\Lambda` is the Coulomb logarithm. This takes the place of the factor
   :math:`\ln\Lambda[\mathrm{erf}(x)-2x\exp(-x^2)/\sqrt{\pi}]` in the :cite:t:`chandrasekhar_dynamical_1943` formula for a
   collisionless medium.

   This drag is applied to every mass distribution component whose kinematics are collisional. This includes the hot halo and,
   for the standard disk, spheroid, and nuclear star cluster components, the interstellar medium, whose kinematics are given by
   ``[kinematicsDistributionInterstellarMedium]`` (by default a uniform effective sound speed). For a disk, :math:`V` is the
   speed of the perturber relative to the rotating gas. Note that the drag is computed as for a uniform medium, which is a poor
   approximation for a perturber crossing a thin gaseous disk.
   </description>
   <default>kimKim2007</default>
   <method name="forceFactor" >
    <description>
    Return the dimensionless drag force, :math:`\mathcal{I}`, for a perturber moving at Mach number ``machNumber`` through a
    gaseous medium, given the Coulomb logarithm, ``logarithmCoulomb``.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: machNumber, logarithmCoulomb</argument>
   </method>
  </functionClass>
  !!]

end module Satellite_Dynamical_Friction_Gaseous
