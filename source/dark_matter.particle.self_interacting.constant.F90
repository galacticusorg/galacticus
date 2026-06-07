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
Contains a module which implements a selfInteracting dark matter particle class.
!!}

  !![
  <darkMatterParticle name="darkMatterParticleSelfInteractingDarkMatterConstant">
   <description>Provides a selfInteracting dark matter particle.</description>
  </darkMatterParticle>
  !!]
  type, extends(darkMatterParticleSelfInteractingDarkMatter) :: darkMatterParticleSelfInteractingDarkMatterConstant
     !!{
     A selfInteracting dark matter particle class.
     !!}
     private
     class           (darkMatterParticleClass), pointer :: darkMatterParticle_          => null()
     double precision                                   :: crossSectionSelfInteraction_
   contains
     !![
     <methods>
       <method description="Return the self-interaction cross section, $\sigma$, of the dark matter particle in units of cm$^2$ g$^{-1}$." method="crossSectionSelfInteraction" />
       <method description="Return the differential self-interaction cross section, $\mathrm{d}\sigma/\mathrm{d}\Omega$, of the dark matter particle in units of cm$^2$ g$^{-1}$ ster$^{-1}$." method="crossSectionSelfInteractionDifferential"/>
       <method description="Return the momentum transfer self-interaction cross section, $\sigma$, of the dark matter particle in units of cm$^2$ g$^{-1}$." method="crossSectionSelfInteractionMomentumTransfer" />
       <method description="Return the viscosity self-interaction cross section, $\sigma$, of the dark matter particle in units of cm$^2$ g$^{-1}$." method="crossSectionSelfInteractionViscosity" />
     </methods>
     !!]
     final     ::                                                sidmConstantDestructor
     procedure :: mass                                        => sidmConstantMass
     procedure :: crossSectionSelfInteraction                 => sidmConstantCrossSectionSelfInteraction     
     procedure :: crossSectionSelfInteractionDifferential     => sidmConstantCrossSectionSelfInteractionDifferential
     procedure :: crossSectionSelfInteractionDifferentialCos  => sidmConstantCrossSectionSelfInteractionDifferentialCos
     procedure :: crossSectionSelfInteractionMomentumTransfer => sidmConstantCrossSectionMomentumTransfer
     procedure :: crossSectionSelfInteractionViscosity        => sidmConstantCrossSectionViscosity
  end type darkMatterParticleSelfInteractingDarkMatterConstant

  interface darkMatterParticleSelfInteractingDarkMatterConstant
     !!{
     Constructors for the ``{\normalfont \ttfamily selfInteractingDarkMatter}'' dark matter particle class.
     !!}
     module procedure sidmConstantConstructorParameters
     module procedure sidmConstantConstructorInternal
  end interface darkMatterParticleSelfInteractingDarkMatterConstant

contains

  function sidmConstantConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``{\normalfont \ttfamily selfInteractingDarkMatter}'' dark matter particle class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterParticleSelfInteractingDarkMatterConstant)        :: self
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
    self=darkMatterParticleSelfInteractingDarkMatterConstant(crossSectionSelfInteraction,darkMatterParticle_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterParticle_" />
    !!]
    return
  end function sidmConstantConstructorParameters

  function sidmConstantConstructorInternal(crossSectionSelfInteraction,darkMatterParticle_) result(self)
    !!{
    Internal constructor for the ``{\normalfont \ttfamily selfInteractingDarkMatter}'' dark matter particle class.
    !!}
    implicit none
    type            (darkMatterParticleSelfInteractingDarkMatterConstant)                :: self
    class           (darkMatterParticleClass                    ), intent(in   ), target :: darkMatterParticle_
    double precision                                             , intent(in   )         :: crossSectionSelfInteraction
    !![
    <constructorAssign variables="*darkMatterParticle_"/>
    !!]

    self%crossSectionSelfInteraction_=crossSectionSelfInteraction
    return
  end function sidmConstantConstructorInternal

  subroutine sidmConstantDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily selfInteractingDarkMatter} dark matter particle class.
    !!}
    implicit none
    type(darkMatterParticleSelfInteractingDarkMatterConstant), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterParticle_" />
    !!]
    return
  end subroutine sidmConstantDestructor

  double precision function sidmConstantMass(self)
    !!{
    Return the mass, in units of keV, of a self-interacting dark matter particle.
    !!}
    implicit none
    class(darkMatterParticleSelfInteractingDarkMatterConstant), intent(inout) :: self

    sidmConstantMass=self%darkMatterParticle_%mass()
    return
  end function sidmConstantMass

  double precision function sidmConstantCrossSectionSelfInteraction(self,velocityRelative)
    !!{
    Return the self-interaction cross section, in units of cm$^2$ g$^{-1}$, of a self-interacting dark matter particle.
    !!}
    implicit none
    class(darkMatterParticleSelfInteractingDarkMatterConstant), intent(inout) :: self
   
    double precision                                          , intent(in   ) :: velocityRelative

    sidmConstantCrossSectionSelfInteraction=self%crossSectionSelfInteraction_
    return
  end function sidmConstantCrossSectionSelfInteraction

  double precision function sidmConstantCrossSectionSelfInteractionDifferential(self,theta,velocityRelative)
    !!{
    Return the differential self-interaction cross section, $\mathrm{d}\sigma/\mathrm{d}\theta$, in units of cm$^2$ g$^{-1}$
    ster$^{-1}$, of a self-interacting dark matter particle.
    !!}
    implicit none
    class           (darkMatterParticleSelfInteractingDarkMatterConstant), intent(inout) :: self
    double precision                                                     , intent(in   ) :: velocityRelative
    double precision                                                     , intent(in   ) :: theta

    ! Currently isotropic scattering is assumed.
    sidmConstantCrossSectionSelfInteractionDifferential=+self%crossSectionSelfInteraction(velocityRelative) &
         &                                                   *0.5d0                              &
         &                                                   *sin(theta)
    return
  end function sidmConstantCrossSectionSelfInteractionDifferential

  double precision function sidmConstantCrossSectionSelfInteractionDifferentialCos(self,Costheta,velocityRelative)
    !!{
    Return the differential self-interaction cross section, $\mathrm{d}\sigma/\mathrm{d}\cos\theta$, in units of
    cm$^2$ g$^{-1}$ ster$^{-1}$, of a self-interacting dark matter particle.
    !!}
    implicit none
    class(darkMatterParticleSelfInteractingDarkMatterConstant), intent(inout) :: self

    double precision                                          , intent(in   ) :: velocityRelative, Costheta

    sidmConstantCrossSectionSelfInteractionDifferentialCos=self%crossSectionSelfInteraction_
    return
  end function sidmConstantCrossSectionSelfInteractionDifferentialCos

  double precision function sidmConstantCrossSectionMomentumTransfer(self,velocityRelative)
    !!{
    Return the momentum transfer self-interaction cross section, in units of cm$^2$ g$^{-1}$, of a self-interacting dark matter particle.
    !!}
    implicit none
    class(darkMatterParticleSelfInteractingDarkMatterConstant), intent(inout) :: self

    double precision                                          , intent(in   ) :: velocityRelative

    sidmConstantCrossSectionMomentumTransfer=self%crossSectionSelfInteraction_
    return
  end function sidmConstantCrossSectionMomentumTransfer

  double precision function sidmConstantCrossSectionViscosity(self,velocityRelative)
    !!{
    Return the viscosity self-interaction cross section, in units of cm$^2$ g$^{-1}$, of a self-interacting dark matter particle.
    !!}
    implicit none
    class(darkMatterParticleSelfInteractingDarkMatterConstant), intent(inout) :: self
    double precision                                          , intent(in   ) :: velocityRelative

    sidmConstantCrossSectionViscosity=self%crossSectionSelfInteraction_
    return
  end function sidmConstantCrossSectionViscosity
