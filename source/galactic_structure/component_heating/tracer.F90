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
  Implementation of a galactic component heating class which heats components as tracers of the dark matter.
  !!}

  !![
  <galacticComponentHeating name="galacticComponentHeatingTracer" docformat="rst">
    <description>
    A galactic component heating class which treats the material of galactic components as a tracer population heated along with
    the dark matter of their host halo (see :galacticus-class:`massDistributionSphericalHeatedTracer`): a shell of material
    initially at radius :math:`r` is moved to the same final radius as the shell of dark matter initially at :math:`r`. This
    requires that the dark matter of satellites is heated (e.g. by use of a ``heated`` or ``heatedMonotonic`` dark matter profile,
    with tidal heating of satellites), as the same heating then applies to both the dark matter and the galactic components.

    This is the tracer limit, in which the binding energy of each shell of the component is that given by the dark-matter-only
    mass distribution. That distribution treats the baryonic mass of the halo as being distributed like the dark matter, so
    compact, baryon-dominated components (such as compact spheroids or nuclear star clusters) are under-bound in their centers,
    and will expand too much in response to heating. This class should not be used together with the
    :galacticus-class:`nodeOperatorTidalHeatingSpheroids` node operator, which models the same heating differently.
    </description>
  </galacticComponentHeating>
  !!]
  type, extends(galacticComponentHeatingClass) :: galacticComponentHeatingTracer
     !!{RST
     A galactic component heating class which heats components as tracers of the dark matter.
     !!}
     private
   contains
     procedure :: isActive         => tracerIsActive
     procedure :: massDistribution => tracerMassDistribution
  end type galacticComponentHeatingTracer

  interface galacticComponentHeatingTracer
     !!{RST
     Constructors for the :galacticus-class:`galacticComponentHeatingTracer` galactic component heating class.
     !!}
     module procedure tracerConstructorParameters
  end interface galacticComponentHeatingTracer

contains

  function tracerConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`galacticComponentHeatingTracer` galactic component heating class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(galacticComponentHeatingTracer)                :: self
    type(inputParameters               ), intent(inout) :: parameters

    self=galacticComponentHeatingTracer()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function tracerConstructorParameters

  logical function tracerIsActive(self) result(isActive)
    !!{RST
    Return true, as this class heats mass distributions.
    !!}
    implicit none
    class(galacticComponentHeatingTracer), intent(inout) :: self
    !$GLC attributes unused :: self

    isActive=.true.
    return
  end function tracerIsActive

  function tracerMassDistribution(self,massDistribution_,massDistributionDarkMatter_) result(massDistributionHeated)
    !!{RST
    Return the mass distribution heated as a tracer of the dark matter.
    !!}
    use :: Mass_Distributions, only : massDistributionSphericalHeatedTracer
    implicit none
    class(massDistributionSpherical            ), pointer               :: massDistributionHeated
    class(galacticComponentHeatingTracer       ), intent(inout)         :: self
    class(massDistributionSpherical            ), intent(inout), target :: massDistribution_            , massDistributionDarkMatter_
    type (massDistributionSphericalHeatedTracer), pointer               :: massDistributionHeatedTracer_
    !$GLC attributes unused :: self

    allocate(massDistributionHeatedTracer_)
    !![
    <referenceConstruct object="massDistributionHeatedTracer_" constructor="massDistributionSphericalHeatedTracer(massDistribution_=massDistribution_,massDistributionReference_=massDistributionDarkMatter_)"/>
    !!]
    massDistributionHeated => massDistributionHeatedTracer_
    return
  end function tracerMassDistribution
