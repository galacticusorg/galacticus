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
Implements a decaying dark matter particle class.
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
     class           (darkMatterParticleClass), pointer :: darkMatterParticle_ => null()
     double precision                                   :: lifetime_                    , massSplitting_, &
          &                                                velocityKick_
   contains
     !![
     <methods>
       <method description="Return the lifetime of the dark matter particle." method="lifetime"     />
       <method description="Return the mass splitting of the decay."          method="massSplitting"/>
       <method description="Return the velocity kick imparted by the decay."  method="velocityKick" />
     </methods>
     !!]
     final     ::                  decayingDMDestructor
     procedure :: mass          => decayingDMMass
     procedure :: lifetime      => decayingDMLifetime
     procedure :: massSplitting => decayingDMMassSplitting
     procedure :: velocityKick  => decayingDMVelocityKick
  end type darkMatterParticleDecayingDarkMatter

  interface darkMatterParticleDecayingDarkMatter
     !!{
     Constructors for the {\normalfont \ttfamily decayingDarkMatter} dark matter particle class.
     !!}
     module procedure decayingDMConstructorParameters
     module procedure decayingDMConstructorInternal
  end interface darkMatterParticleDecayingDarkMatter

contains

  function decayingDMConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily decayingDarkMatter} dark matter particle class which takes a parameter set as input.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (darkMatterParticleDecayingDarkMatter)                :: self
    type            (inputParameters                     ), intent(inout) :: parameters
    class           (darkMatterParticleClass             ), pointer       :: darkMatterParticle_
    double precision                                                      :: lifetime           , massSplitting, &
         &                                                                   velocityKick

    !![
    <inputParameter>
      <name>lifetime</name>
      <source>parameters</source>
      <variable>lifetime</variable>
      <description>Lifetime of the dark matter particle in Gyr.</description>
    </inputParameter>
    !!]
    if      (parameters%isPresent('massSplitting')) then
       if (parameters%isPresent('velocityKick')) call Error_Report('specify [massSplitting] or [velocityKick], not both'//{introspection:location})
       !![
       <inputParameter>
	 <name>massSplitting</name>
	 <source>parameters</source>
	 <description>Mass splitting of the dark matter decay.</description>
       </inputParameter>
       !!]
    else if (parameters%isPresent('velocityKick' )) then
       !![
       <inputParameter>
	 <name>velocityKick</name>
	 <source>parameters</source>
	 <description>Velocity kick imparted by dark matter decay.</description>
       </inputParameter>
       !!]
    else
       call Error_Report('specify one of [massSplitting] and [velocityKick]'//{introspection:location})
    end if
    !![
    <objectBuilder class="darkMatterParticle"  name="darkMatterParticle_"  source="parameters"/>
    <conditionalCall>
      <call>self=darkMatterParticleDecayingDarkMatter(darkMatterParticle_,lifetime{conditions})</call>
      <argument name="massSplitting" value="massSplitting" parameterPresent="parameters"/>
      <argument name="velocityKick"  value="velocityKick"  parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterParticle_" />
    !!]
    return
  end function decayingDMConstructorParameters

  function decayingDMConstructorInternal(darkMatterParticle_,lifetime,massSplitting,velocityKick) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily decayingDarkMatter} dark matter particle class.
    !!}
    use :: Error                       , only : Error_Report
    use :: Numerical_Constants_Physical, only : speedLight
    use :: Numerical_Constants_Prefixes, only : kilo
    implicit none
    type            (darkMatterParticleDecayingDarkMatter)                          :: self
    class           (darkMatterParticleClass             ), intent(in   ), target   :: darkMatterParticle_
    double precision                                      , intent(in   )           :: lifetime
    double precision                                      , intent(in   ), optional :: massSplitting      , velocityKick
    !![
    <constructorAssign variables="*darkMatterParticle_"/>
    !!]

    self%lifetime_=lifetime
    if      (present(massSplitting)) then
       if (present(velocityKick)) call Error_Report('specify `massSplitting` or `velocityKick`, not both'//{introspection:location})
       self%massSplitting_=+massSplitting
       self%velocityKick_ =+massSplitting &
            &              *speedLight    &
            &              /kilo
    else if (present(velocityKick )) then
       self%velocityKick_ =+velocityKick
       self%massSplitting_=+velocityKick  &
            &              *kilo          &
            &              /speedLight
    else
       call Error_Report('specify one of `massSplitting` and `velocityKick`'//{introspection:location})
    end if
    return
  end function decayingDMConstructorInternal

  subroutine decayingDMDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily decayingDarkMatter} dark matter particle class.
    !!}
    implicit none
    type(darkMatterParticleDecayingDarkMatter), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterParticle_"/>
    !!]
    return
  end subroutine decayingDMDestructor

  double precision function decayingDMMass(self) result(mass)
    !!{
    Return the mass, in units of keV, of a decaying dark matter particle.
    !!}
    implicit none
    class(darkMatterParticleDecayingDarkMatter), intent(inout) :: self

    mass=self%darkMatterParticle_%mass()
    return
  end function decayingDMMass

  double precision function decayingDMLifetime(self) result(lifetime)
    !!{
    Return the lifetime, in units of Gyr, of a decaying dark matter particle.
    !!}
    implicit none
    class(darkMatterParticleDecayingDarkMatter), intent(inout) :: self

    lifetime=self%lifetime_
    return
  end function decayingDMLifetime

  double precision function decayingDMMassSplitting(self) result(massSplitting)
    !!{
    Return the mass splitting of a decaying dark matter particle.
    !!}
    implicit none
    class(darkMatterParticleDecayingDarkMatter), intent(inout) :: self

    massSplitting=self%massSplitting_
    return
  end function decayingDMMassSplitting

  double precision function decayingDMVelocityKick(self) result(velocityKick)
    !!{
    Return the velocity kick imparted by a decaying dark matter particle.
    !!}
    implicit none
    class(darkMatterParticleDecayingDarkMatter), intent(inout) :: self

    velocityKick=self%velocityKick_
    return
  end function decayingDMVelocityKick
