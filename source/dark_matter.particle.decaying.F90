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
     double precision                                   :: lifetime_, massSplitting_, gamma_
     logical                                            :: heating_, massLoss_
   contains
     !![
     <methods>
       <method description="Return the lifetime of the dark matter particle." method="lifetime" />
       <method description="Return the mass splitting of the decay." method="massSplitting" />
       <method description="Return boolean that indicates if heating from decays should be included." method="heating" />
       <method description="Return boolean that indicates if mass loss from decays should be included." method="massLoss" />
       <method description="Return free parameter in heating in response to mass loss." method="gamma" />
     </methods>
     !!]
     final     ::                         decayingDMDestructor
     procedure :: mass                 => decayingDMMass
     procedure :: lifetime             => decayingDMLifetime
     procedure :: massSplitting        => decayingDMMassSplitting
     procedure :: heating              => decayingDMHeating
     procedure :: massLoss             => decayingDMMassLoss
     procedure :: gamma                => decayingDMGamma
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
    double precision                                                      :: lifetime, massSplitting, gamma
    logical                                                               :: heating, massLoss

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
    <inputParameter>
      <name>heating</name>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
      <variable>heating</variable>
      <description>Boolean that indicates if heating from decays should be included.</description>
    </inputParameter>
    <inputParameter>
      <name>massLoss</name>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
      <variable>massLoss</variable>
      <description>Boolean that indicates if mass loss should be included.</description>
    </inputParameter>
    <inputParameter>
      <name>gamma</name>
      <defaultValue>1.0d0</defaultValue>
      <source>parameters</source>
      <variable>gamma</variable>
      <description>Free parameter in heating in response to mass loss.</description>
    </inputParameter>
    <objectBuilder class="darkMatterParticle"  name="darkMatterParticle_"  source="parameters"/>
    !!]
    self=darkMatterParticleDecayingDarkMatter(lifetime,massSplitting,heating,massLoss,gamma,darkMatterParticle_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterParticle_" />
    !!]
    return
  end function decayingDMConstructorParameters

  function decayingDMConstructorInternal(lifetime,massSplitting,heating,massLoss,gamma,darkMatterParticle_) result(self)
    !!{
    Internal constructor for the ``{\normalfont \ttfamily decayingDarkMatter}'' dark matter particle class.
    !!}
    implicit none
    type            (darkMatterParticleDecayingDarkMatter)                        :: self
    class           (darkMatterParticleClass             ), intent(in   ), target :: darkMatterParticle_
    double precision                                      , intent(in   )         :: lifetime, massSplitting, gamma
    logical                                               , intent(in   )         :: heating, massLoss
    !![
    <constructorAssign variables="*darkMatterParticle_"/>
    !!]
    self%lifetime_ = lifetime
    self%massSplitting_ = massSplitting
    self%heating_ = heating
    self%massLoss_ = massLoss
    self%gamma_ = gamma 
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
  end subroutine decayingDMDestructor

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

  logical function decayingDMHeating(self)
    !!{
    Return boolean that indicates if heating from decays is included.
    !!}
    implicit none
    class(darkMatterParticleDecayingDarkMatter), intent(inout) :: self

    decayingDMHeating=self%heating_
    return
  end function decayingDMHeating

  logical function decayingDMMassLoss(self)
    !!{
    Return boolean that indicates is mass loss is included.
    !!}
    implicit none
    class(darkMatterParticleDecayingDarkMatter), intent(inout) :: self

    decayingDMMassLoss=self%massLoss_
    return
  end function decayingDMMassLoss

  double precision function decayingDMGamma(self)
    !!{
    Return the free parameter in the heating in response to mass loss.
    !!}
    implicit none
    class(darkMatterParticleDecayingDarkMatter), intent(inout) :: self

    decayingDMGamma=self%gamma_
    return
  end function decayingDMGamma
