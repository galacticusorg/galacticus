!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
Implements a selfInteracting dark matter particle class.
!!}

  !![
  <darkMatterParticle name="darkMatterParticleSelfInteractingDarkMatter">
   <description>Provides a selfInteracting dark matter particle.</description>
  </darkMatterParticle>
  !!]
  type, extends(darkMatterParticleClass) :: darkMatterParticleSelfInteractingDarkMatter
     !!{
     A selfInteracting dark matter particle class.
     !!}
     private
     class           (darkMatterParticleClass), pointer :: darkMatterParticle_          => null()
     double precision                                   :: crossSectionSelfInteraction_
   contains
     !![
     <methods>
       <method description="Return the self-interaction cross section, $\sigma$, of the dark matter particle in units of cm$^2$ g$^{-1}$."                                                     method="crossSectionSelfInteraction"            />
       <method description="Return the differential self-interaction cross section, $\mathrm{d}\sigma/\mathrm{d}\Omega$, of the dark matter particle in units of cm$^2$ g$^{-1}$ ster$^{-1}$." method="crossSectionSelfInteractionDifferential"/>
     </methods>
     !!]
     final     ::                                            selfInteractingDMDestructor
     procedure :: mass                                    => selfInteractingDMMass
     procedure :: crossSectionSelfInteraction             => selfInteractingDMCrossSectionSelfInteraction
     procedure :: crossSectionSelfInteractionDifferential => selfInteractingDMCrossSectionSelfInteractionDifferential
  end type darkMatterParticleSelfInteractingDarkMatter

  interface darkMatterParticleSelfInteractingDarkMatter
     !!{
     Constructors for the \normalfont \ttfamily selfInteractingDarkMatter} dark matter particle class.
     !!}
     module procedure selfInteractingDMConstructorParameters
     module procedure selfInteractingDMConstructorInternal
  end interface darkMatterParticleSelfInteractingDarkMatter

contains

  function selfInteractingDMConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \normalfont \ttfamily selfInteractingDarkMatter} dark matter particle class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterParticleSelfInteractingDarkMatter)                :: self
    type            (inputParameters                            ), intent(inout) :: parameters
    class           (darkMatterParticleClass                    ), pointer       :: darkMatterParticle_
    double precision                                                             :: crossSectionSelfInteraction

    !![
    <inputParameter>
      <name>crossSectionSelfInteraction</name>
      <source>parameters</source>
      <description>The self-interaction cross section in units of cm$^2$ g$^{-1}$.</description>
    </inputParameter>
    <objectBuilder class="darkMatterParticle"  name="darkMatterParticle_"  source="parameters"/>
    !!]
    self=darkMatterParticleSelfInteractingDarkMatter(crossSectionSelfInteraction,darkMatterParticle_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterParticle_" />
    !!]
    return
  end function selfInteractingDMConstructorParameters

  function selfInteractingDMConstructorInternal(crossSectionSelfInteraction,darkMatterParticle_) result(self)
    !!{
    Internal constructor for the \Normalfont \ttfamily selfInteractingDarkMatter} dark matter particle class.
    !!}
    implicit none
    type            (darkMatterParticleSelfInteractingDarkMatter)                        :: self
    class           (darkMatterParticleClass                    ), intent(in   ), target :: darkMatterParticle_
    double precision                                             , intent(in   )         :: crossSectionSelfInteraction
    !![
    <constructorAssign variables="*darkMatterParticle_"/>
    !!]

    self%crossSectionSelfInteraction_=crossSectionSelfInteraction
    return
  end function selfInteractingDMConstructorInternal

  subroutine selfInteractingDMDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily selfInteractingDarkMatter} dark matter particle class.
    !!}
    implicit none
    type(darkMatterParticleSelfInteractingDarkMatter), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterParticle_" />
    !!]
    return
  end subroutine selfInteractingDMDestructor

  double precision function selfInteractingDMMass(self)
    !!{
    Return the mass, in units of keV, of a self-interacting dark matter particle.
    !!}
    implicit none
    class(darkMatterParticleSelfInteractingDarkMatter), intent(inout) :: self

    selfInteractingDMMass=self%darkMatterParticle_%mass()
    return
  end function selfInteractingDMMass

  double precision function selfInteractingDMCrossSectionSelfInteraction(self)
    !!{
    Return the self-interaction cross section, in units of cm$^2$ g$^{-1}$, of a self-interacting dark matter particle.
    !!}
    implicit none
    class(darkMatterParticleSelfInteractingDarkMatter), intent(inout) :: self

    selfInteractingDMCrossSectionSelfInteraction=self%crossSectionSelfInteraction_
    return
  end function selfInteractingDMCrossSectionSelfInteraction

  double precision function selfInteractingDMCrossSectionSelfInteractionDifferential(self,theta)
    !!{
    Return the differential self-interaction cross section, $\mathrm{d}\sigma/\mathrm{d}\theta$, in units of cm$^2$ g$^{-1}$
    ster$^{-1}$, of a self-interacting dark matter particle.
    !!}
    implicit none
    class           (darkMatterParticleSelfInteractingDarkMatter), intent(inout) :: self
    double precision                                             , intent(in   ) :: theta

    ! Currently isotropic scattering is assumed.
    selfInteractingDMCrossSectionSelfInteractionDifferential=+self%crossSectionSelfInteraction() &
         &                                                   *0.5d0                              &
         &                                                   *sin(theta)
    return
  end function selfInteractingDMCrossSectionSelfInteractionDifferential
