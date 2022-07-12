!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
Contains a module which implements a decaying dark matter particle class.
!!}

  !![
  <darkMatterParticle name="darkMatterParticleDecayingDarkMatter">
   <description>Provides a decaying dark matter particle.</description>
  </darkMatterParticle>
  !!]
  type, extends(darkMatterParticleClass) :: darkMatterParticleDecayingDarkMatter
     !!{
     A decaying dark matter particle class.
     !!}
     private
     class           (darkMatterParticleClass), pointer :: darkMatterParticle_          => null()
     double precision                                   :: lifetime_, massSplitting_
   contains
     !![
     <methods>
       <method description="Return the self-interaction cross section, $\sigma$, of the dark matter particle in units of cm$^2$ g$^{-1}$." method="crossSectionSelfInteraction" />
       <method description="Return the differential self-interaction cross section, $\mathrm{d}\sigma/\mathrm{d}\Omega$, of the dark matter particle in units of cm$^2$ g$^{-1}$ ster$^{-1}$." method="crossSectionSelfInteractionDifferential" />
     </methods>
     !!]
     final     ::                         decayingDMDestructor
     procedure :: mass                 => decayingDMMass
     procedure :: lifetime             => decayingDMLifetime
     procedure :: massSplitting        => decayingDMMassSplitting
  end type darkMatterParticleDecayingDarkMatter

  interface darkMatterParticleDecayingDarkMatter
     !!{
     Constructors for the ``{\normalfont \ttfamily decayingDarkMatter}'' dark matter particle class.
     !!}
     module procedure decayingDMConstructorParameters
     module procedure decayingDMConstructorInternal
  end interface darkMatterParticleDecayingDarkMatter

contains

  function decayingDMConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``{\normalfont \ttfamily decayingDarkMatter}'' dark matter particle class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterParticleDecayingDarkMatter)                :: self
    type            (inputParameters                     ), intent(inout) :: parameters
    class           (darkMatterParticleClass             ), pointer       :: darkMatterParticle_
    double precision                                                      :: lifetime, massSplitting

    !![
    <inputParameter>
      <name>lifetime</name>
      <source>parameters</source>
      <variable>lifetime</variable>
      <description>Lifetime of the dark matter particle in Gyr.</description>
    </inputParameter>
    <inputParameter>
      <name>massSplitting</name>
      <source>parameters</source>
      <variable>massSplitting</variable>
      <description>Mass splitting of the dark matter decay.</description>
    </inputParameter>
    <objectBuilder class="darkMatterParticle"  name="darkMatterParticle_"  source="parameters"/>
    !!]
    self=darkMatterParticleSelfInteractingDarkMatter(lifetime,massSplitting,darkMatterParticle_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterParticle_" />
    !!]
    return
  end function decayingDMConstructorParameters

  function decayingDMConstructorInternal(lifetime,massSplitting,darkMatterParticle_) result(self)
    !!{
    Internal constructor for the ``{\normalfont \ttfamily decayingDarkMatter}'' dark matter particle class.
    !!}
    implicit none
    type            (darkMatterParticleDecayingDarkMatter)                        :: self
    class           (darkMatterParticleClass             ), intent(in   ), target :: darkMatterParticle_
    double precision                                      , intent(in   )         :: lifetime, massSplitting
    !![
    <constructorAssign variables="*darkMatterParticle_"/>
    !!]

    self%lifetime_=lifetime
    self%massSplitting_=massSplitting
    return
  end function decayingDMConstructorInternal

  subroutine decayingDMDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily decayingDarkMatter} dark matter particle class.
    !!}
    implicit none
    type(darkMatterParticleDecayingDarkMatter), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterParticle_" />
    !!]
    return
  end subroutine selfInteractingDMDestructor

  double precision function decayingDMMass(self)
    !!{
    Return the mass, in units of keV, of a decaying dark matter particle.
    !!}
    implicit none
    class(darkMatterParticleDecayingDarkMatter), intent(inout) :: self

    decayingDMMass=self%darkMatterParticle_%mass()
    return
  end function decayingDMMass

  double precision function decayingDMLifetime(self)
    !!{
    Return the lifetime, in units of Gyr, of a decaying dark matter particle.
    !!}
    implicit none
    class(darkMatterParticleDecayingDarkMatter), intent(inout) :: self

    decayingDMLifetime=self%lifetime_
    return
  end function decayingDMLifetime

  double precision function decayingDMMassSplitting(self)
    !!{
    Return the mass splitting of a decaying dark matter particle.
    !!}
    implicit none
    class(darkMatterParticleDecayingDarkMatter), intent(inout) :: self

    decayingDMMassSplitting=self%massSplitting_
    return
  end function decayingDMMassSplitting
