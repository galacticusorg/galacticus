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
  Implements the covering fraction of the dusty torus around an AGN in the receding torus model.
  !!}

  use :: Accretion_Disks           , only : accretionDisksClass
  use :: Black_Hole_Accretion_Rates, only : blackHoleAccretionRateClass

  !![
  <blackHoleTorusCoveringFraction name="blackHoleTorusCoveringFractionRecedingTorus" docformat="rst">
   <description>
   The covering fraction of the dusty torus around an AGN in the modified receding torus model of
   :cite:t:`simpson_luminosity_2005`, a generalization of that of :cite:t:`lawrence_relative_1991`:

   .. math::

      f_\mathrm{cov} = \left[ 1 + 3 \left( \frac{L}{L_0} \right)^{1-2\xi} \right]^{-1/2}.

   In the receding torus the inner edge of the torus lies at the dust sublimation radius, which grows as
   :math:`L^{1/2}`, while its height grows as :math:`L^\xi`; so the torus becomes thinner, and covers less of the sky, as
   the AGN becomes more luminous. :math:`\xi=0` is the original model of :cite:t:`lawrence_relative_1991`, and
   :math:`\xi=1/2` gives a covering fraction of :math:`1/2` independent of luminosity.

   :math:`L` is the bolometric luminosity of the accretion disk of the central black hole,
   :math:`\epsilon \dot{M} c^2`, where the accretion rate :math:`\dot{M}` (summed over all of its sources) is given by
   the ``[blackHoleAccretionRate]`` object and the radiative efficiency :math:`\epsilon` by the ``[accretionDisks]``
   object---the same luminosity from which the spectrum of the accretion disk is computed. A galaxy with no accreting
   black hole is fully covered.

   The defaults are the best fit of :cite:t:`simpson_luminosity_2005`, :math:`\xi=0.23` and
   :math:`L_0=10^{35.37}\,\hbox{W}` in [O III] luminosity, converted to a bolometric luminosity with the correction
   :math:`L/L_\mathrm{[O III]}=3500` of :cite:t:`heckman_present-day_2004`, which gives
   :math:`L_0=2.13\times 10^{12}\,L_\odot`.
   </description>
  </blackHoleTorusCoveringFraction>
  !!]
  type, extends(blackHoleTorusCoveringFractionClass) :: blackHoleTorusCoveringFractionRecedingTorus
     !!{RST
     The covering fraction of the dusty torus around an AGN in the receding torus model.
     !!}
     private
     class           (blackHoleAccretionRateClass), pointer :: blackHoleAccretionRate_  => null()
     class           (accretionDisksClass        ), pointer :: accretionDisks_          => null()
     double precision                                       :: luminosityCharacteristic          , exponentHeight
   contains
     final     ::                     recedingTorusDestructor
     procedure :: coveringFraction => recedingTorusCoveringFraction
  end type blackHoleTorusCoveringFractionRecedingTorus

  interface blackHoleTorusCoveringFractionRecedingTorus
     !!{RST
     Constructors for the :galacticus-class:`blackHoleTorusCoveringFractionRecedingTorus` torus covering fraction class.
     !!}
     module procedure recedingTorusConstructorParameters
     module procedure recedingTorusConstructorInternal
  end interface blackHoleTorusCoveringFractionRecedingTorus

contains

  function recedingTorusConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`blackHoleTorusCoveringFractionRecedingTorus` torus covering fraction class
    which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (blackHoleTorusCoveringFractionRecedingTorus)                :: self
    type            (inputParameters                            ), intent(inout) :: parameters
    class           (blackHoleAccretionRateClass                ), pointer       :: blackHoleAccretionRate_
    class           (accretionDisksClass                        ), pointer       :: accretionDisks_
    double precision                                                             :: luminosityCharacteristic, exponentHeight

    !![
    <inputParameter docformat="rst">
      <name>luminosityCharacteristic</name>
      <defaultValue>2.13d12</defaultValue>
      <defaultSource>:cite:t:`simpson_luminosity_2005` with the bolometric correction of :cite:t:`heckman_present-day_2004`</defaultSource>
      <description>
      The bolometric luminosity, :math:`L_0` in :math:`L_\odot`, at which the covering fraction of the torus is
      :math:`1/2`.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>exponentHeight</name>
      <defaultValue>0.23d0</defaultValue>
      <defaultSource>:cite:t:`simpson_luminosity_2005`</defaultSource>
      <description>
      The exponent, :math:`\xi`, of the dependence of the height of the torus on luminosity, :math:`h \propto L^\xi`.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="blackHoleAccretionRate" name="blackHoleAccretionRate_" source="parameters"/>
    <objectBuilder class="accretionDisks"         name="accretionDisks_"         source="parameters"/>
    !!]
    self=blackHoleTorusCoveringFractionRecedingTorus(luminosityCharacteristic,exponentHeight,blackHoleAccretionRate_,accretionDisks_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="blackHoleAccretionRate_"/>
    <objectDestructor name="accretionDisks_"        />
    !!]
    return
  end function recedingTorusConstructorParameters

  function recedingTorusConstructorInternal(luminosityCharacteristic,exponentHeight,blackHoleAccretionRate_,accretionDisks_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`blackHoleTorusCoveringFractionRecedingTorus` torus covering fraction
    class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (blackHoleTorusCoveringFractionRecedingTorus)                        :: self
    double precision                                             , intent(in   )         :: luminosityCharacteristic, exponentHeight
    class           (blackHoleAccretionRateClass                ), intent(in   ), target :: blackHoleAccretionRate_
    class           (accretionDisksClass                        ), intent(in   ), target :: accretionDisks_
    !![
    <constructorAssign variables="luminosityCharacteristic, exponentHeight, *blackHoleAccretionRate_, *accretionDisks_"/>
    !!]

    if (luminosityCharacteristic <= 0.0d0) call Error_Report('`luminosityCharacteristic` must be positive'//{introspection:location})
    return
  end function recedingTorusConstructorInternal

  subroutine recedingTorusDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`blackHoleTorusCoveringFractionRecedingTorus` torus covering fraction class.
    !!}
    implicit none
    type(blackHoleTorusCoveringFractionRecedingTorus), intent(inout) :: self

    !![
    <objectDestructor name="self%blackHoleAccretionRate_"/>
    <objectDestructor name="self%accretionDisks_"        />
    !!]
    return
  end subroutine recedingTorusDestructor

  double precision function recedingTorusCoveringFraction(self,node) result(coveringFraction)
    !!{RST
    Return the covering fraction of the dusty torus in the receding torus model.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentBlackHole
    use :: Numerical_Constants_Astronomical, only : gigaYear              , luminositySolar, massSolar
    use :: Numerical_Constants_Physical    , only : speedLight
    implicit none
    class           (blackHoleTorusCoveringFractionRecedingTorus), intent(inout) :: self
    type            (treeNode                                   ), intent(inout) :: node
    class           (nodeComponentBlackHole                     ), pointer       :: blackHole
    double precision                                                             :: rateAccretionSpheroid, rateAccretionHotHalo           , &
         &                                                                          rateAccretion        , rateAccretionNuclearStarCluster, &
         &                                                                          luminosityBolometric

    coveringFraction=1.0d0
    if (node%blackHoleCount() < 1) return
    blackHole => node%blackHole()
    call self%blackHoleAccretionRate_%rateAccretion(blackHole,rateAccretionSpheroid,rateAccretionHotHalo,rateAccretionNuclearStarCluster)
    rateAccretion=+rateAccretionSpheroid           &
         &        +rateAccretionHotHalo            &
         &        +rateAccretionNuclearStarCluster
    if (rateAccretion <= 0.0d0) return
    luminosityBolometric=+self%accretionDisks_%efficiencyRadiative(blackHole,rateAccretion) &
         &               *rateAccretion                                                     &
         &               *massSolar                                                         &
         &               *speedLight**2                                                     &
         &               /gigaYear                                                          &
         &               /luminositySolar
    if (luminosityBolometric <= 0.0d0) return
    coveringFraction=1.0d0/sqrt(1.0d0+3.0d0*(luminosityBolometric/self%luminosityCharacteristic)**(1.0d0-2.0d0*self%exponentHeight))
    return
  end function recedingTorusCoveringFraction
