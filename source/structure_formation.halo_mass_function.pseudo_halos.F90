!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
Contains a module which implements a dark matter halo mass function class which modifies another mass function by adding in a
population of pseudo-halos.
!!}

use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <haloMassFunction name="haloMassFunctionPseudoHalos">
   <description>
    The halo mass function is computed by adding a population of pseudo-halos to another halo mass function.
   </description>
  </haloMassFunction>
  !!]
  type, extends(haloMassFunctionClass) :: haloMassFunctionPseudoHalos
     !!{
     A halo mass function class that modifies another mass function by adding in a population of pseudo-halos.
     !!}
     private
     double precision                                   :: massZeroPointReference          , normalization       , &
          &                                                exponentMass                    , exponentRedshift    , &
         &                                                 massParticleReference           , exponentMassParticle, &
         &                                                 massParticle                    , massZeroPoint
     class           (haloMassFunctionClass  ), pointer :: massFunction_          => null()
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_    => null()
   contains
     final     ::                 pseudoHalosDestructor
     procedure :: differential => pseudoHalosDifferential
     procedure :: integrated   => pseudoHalosIntegrated
  end type haloMassFunctionPseudoHalos

  interface haloMassFunctionPseudoHalos
     !!{
     Constructors for the {\normalfont \ttfamily pseudoHalos} halo mass function class.
     !!}
     module procedure pseudoHalosConstructorParameters
     module procedure pseudoHalosConstructorInternal
  end interface haloMassFunctionPseudoHalos

contains

  function pseudoHalosConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily pseudoHalos} halo mass function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (haloMassFunctionPseudoHalos)                :: self
    type            (inputParameters            ), intent(inout) :: parameters
    class           (haloMassFunctionClass      ), pointer       :: massFunction_
    class           (cosmologyParametersClass   ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass    ), pointer       :: cosmologyFunctions_
    double precision                                             :: massZeroPointReference, normalization       , &
         &                                                          exponentMass          , exponentRedshift    , &
         &                                                          massParticleReference , exponentMassParticle, &
         &                                                          massParticle

    !![
    <inputParameter>
      <name>massZeroPointReference</name>
      <source>parameters</source>
      <description>The mass zero-point reference $M^\prime_0$ in the relation $M_0(m_\mathrm{p}) = M^\prime_0 (m_\mathrm{p}/m^\prime_\mathrm{p})^\gamma$.</description>
    </inputParameter>
    <inputParameter>
      <name>massParticleReference</name>
      <source>parameters</source>
      <description>The particle mass reference $m^\prime_\mathrm{p}$ in the relation $M_0(m_\mathrm{p}) = M^\prime_0 (m_\mathrm{p}/m^\prime_\mathrm{p})^\gamma$.</description>
    </inputParameter>
    <inputParameter>
      <name>massParticle</name>
      <source>parameters</source>
      <description>The particle mass $m_\mathrm{p}$ in the relation $M_0(m_\mathrm{p}) = M^\prime_0 (m_\mathrm{p}/m^\prime_\mathrm{p})^\gamma$.</description>
    </inputParameter>
    <inputParameter>
      <name>exponentMassParticle</name>
      <source>parameters</source>
      <description>The exponent $\gamma$ in the relation $M_0(m_\mathrm{p}) = M^\prime_0 (m_\mathrm{p}/m^\prime_\mathrm{p})^\gamma$.</description>
    </inputParameter>
    <inputParameter>
      <name>normalization</name>
      <source>parameters</source>
      <description>The normalization $n_0$ in the pseudo-halo mass function $n(M) = n_0 [M/M_0(m_\mathrm{p})])^\alpha (1+z)^\beta$.</description>
    </inputParameter>
    <inputParameter>
      <name>exponentMass</name>
      <source>parameters</source>
      <description>The exponent $\alpha$ in the pseudo-halo mass function $n(M) = n_0 [M/M_0(m_\mathrm{p})])^\alpha (1+z)^\beta$.</description>
    </inputParameter>
    <inputParameter>
      <name>exponentRedshift</name>
      <source>parameters</source>
      <description>The exponent $\beta$ in the pseudo-halo mass function $n(M) = n_0 [M/M_0(m_\mathrm{p})])^\alpha (1+z)^\beta$.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="haloMassFunction"    name="massFunction_"        source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    !!]
    self=haloMassFunctionPseudoHalos(normalization,exponentMass,exponentRedshift,massZeroPointReference,massParticleReference,massParticle,exponentMassParticle,massFunction_,cosmologyParameters_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="massFunction_"       />
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="cosmologyFunctions_" />
    !!]
    return
  end function pseudoHalosConstructorParameters

  function pseudoHalosConstructorInternal(normalization,exponentMass,exponentRedshift,massZeroPointReference,massParticleReference,massParticle,exponentMassParticle,massFunction_,cosmologyParameters_,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily pseudoHalos} halo mass function class.
    !!}
    implicit none
    type            (haloMassFunctionPseudoHalos)                        :: self
    class           (haloMassFunctionClass      ), target, intent(in   ) :: massFunction_
    class           (cosmologyParametersClass   ), target, intent(in   ) :: cosmologyParameters_
    class           (cosmologyFunctionsClass    ), target, intent(in   ) :: cosmologyFunctions_
    double precision                                     , intent(in   ) :: massZeroPointReference, normalization       , &
         &                                                                  exponentMass          , exponentRedshift    , &
         &                                                                  massParticleReference , exponentMassParticle, &
         &                                                                  massParticle
    !![
    <constructorAssign variables="normalization, exponentMass, exponentRedshift, massZeroPointReference, massParticleReference, massParticle, exponentMassParticle, *cosmologyParameters_, *cosmologyFunctions_, *massFunction_"/>
    !!]

    self%massZeroPoint=+  massZeroPointReference &
         &             *(                        &
         &               +massParticle           &
         &               /massParticleReference  &
         &              )**exponentMassParticle
    return
  end function pseudoHalosConstructorInternal

  subroutine pseudoHalosDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily pseudoHalos} halo mass function class.
    !!}
    implicit none
    type(haloMassFunctionPseudoHalos), intent(inout) :: self

    !![
    <objectDestructor name="self%massFunction_"       />
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%cosmologyFunctions_" />
    !!]
    return
  end subroutine pseudoHalosDestructor

  double precision function pseudoHalosDifferential(self,time,mass,node) result(massFunction)
    !!{
    Return the differential halo mass function at the given time and mass.
    !!}
    implicit none
    class           (haloMassFunctionPseudoHalos), intent(inout), target   :: self
    double precision                             , intent(in   )           :: time, mass
    type            (treeNode                   ), intent(inout), optional :: node

    massFunction=+   self%normalization                                                              &
         &       *(                                                                                  &
         &         +      mass                                                                       &
         &         / self%massZeroPoint                                                              &
         &        )**self%exponentMass                                                               &
         &       /   self%cosmologyFunctions_%expansionFactor(time          )**self%exponentRedshift &
         &       +   self%massFunction_      %differential   (time,mass,node)
    return
  end function pseudoHalosDifferential

  double precision function pseudoHalosIntegrated(self,time,massLow,massHigh,node,status) result(massFunction)
    !!{
    Return the integrated halo mass function at the given time and mass.
    !!}
    implicit none
    class           (haloMassFunctionPseudoHalos), intent(inout), target           :: self
    double precision                             , intent(in   )                   :: time       , massLow     , &
         &                                                                            massHigh
    type            (treeNode                   ), intent(inout), target, optional :: node 
    integer                                      , intent(  out)        , optional :: status
    double precision                                                               :: integralLow, integralHigh

    integralLow =+(massLow /self%massZeroPoint)**(self%exponentMass+1.0d0)
    integralHigh=+(massHigh/self%massZeroPoint)**(self%exponentMass+1.0d0)
    massFunction=+self%normalization                                                                                   &
         &       *self%massZeroPoint                                                                                   &
         &       *(                                                                                                    &
         &         +integralHigh                                                                                       &
         &         -integralLow                                                                                        &
         &        )                                                                                                    &
         &       /(                                                                                                    &
         &         +self%exponentMass                                                                                  &
         &         +1.0d0                                                                                              &
         &        )                                                                                                    &
         &       /  self%cosmologyFunctions_%expansionFactor(time                             )**self%exponentRedshift &
         &       +  self%massFunction_      %integrated     (time,massLow,massHigh,node,status)
    return
  end function pseudoHalosIntegrated
  
