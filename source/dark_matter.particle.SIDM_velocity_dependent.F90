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

!+    Contributions to this file made by: Niusha Ahvazi

!!{
Contains a module which implements a selfInteracting dark matter particle class.
!!}

  !![
  <darkMatterParticle name="darkMatterParticleSIDMVelocityDependent">
   <description>Provides a selfInteracting dark matter particle with velocity dependent cross section.</description>
  </darkMatterParticle>
  !!]
  type, extends(darkMatterParticleSelfInteractingDarkMatter) :: darkMatterParticleSIDMVelocityDependent
     !!{
     A selfInteracting dark matter particle class.
     !!}
     private
     class           (darkMatterParticleClass), pointer :: darkMatterParticle_ => null()
     double precision                                   :: sigma0                       , velocityCharacteristic_
   contains
     !![
     <methods>
       <method method="crossSectionSelfInteraction"                 description="Return the self-interaction cross section, $\sigma$, of the dark matter particle in units of cm$^2$ g$^{-1}$."                                                    />
       <method method="crossSectionSelfInteractionDifferential"     description="Return the differential self-interaction cross section, $\mathrm{d}\sigma/\mathrm{d}\Omega$, of the dark matter particle in units of cm$^2$ g$^{-1}$ ster$^{-1}$."/>
       <method method="crossSectionSelfInteractionMomentumTransfer" description="Return the momentum transfer self-interaction cross section, $\sigma$, of the dark matter particle in units of cm$^2$ g$^{-1}$."                                  />
       <method method="crossSectionSelfInteractionViscosity"        description="Return the viscosity self-interaction cross section, $\sigma$, of the dark matter particle in units of cm$^2$ g$^{-1}$."                                          />
       <method method="velocityCharacteristic"                      description="Return the characteristic velocity scale, $v_\mathrm{c}$, of the velocity-dependent cross section, in units of km s$^{-1}$."                                       />
     </methods>
     !!]
     final     ::                                                sidmVelocityDependentDestructor
     procedure :: mass                                        => sidmVelocityDependentMass
     procedure :: velocityCharacteristic                      => sidmVelocityDependentVelocityCharacteristic
     procedure :: crossSectionSelfInteraction                 => sidmVelocityDependentCrossSectionSelfInteraction
     procedure :: crossSectionSelfInteractionDifferential     => sidmVelocityDependentCrossSectionSelfInteractionDifferential
     procedure :: crossSectionSelfInteractionDifferentialCos  => sidmVelocityDependentCrossSectionSelfInteractionDifferentialCos
     procedure :: crossSectionSelfInteractionMomentumTransfer => sidmVelocityDependentCrossSectionMomentumTransfer
     procedure :: crossSectionSelfInteractionViscosity        => sidmVelocityDependentCrossSectionViscosity
  end type darkMatterParticleSIDMVelocityDependent

  interface darkMatterParticleSIDMVelocityDependent
     !!{
     Constructors for the ``\mono{selfInteractingDarkMatter}'' dark matter particle class.
     !!}
     module procedure sidmVelocityDependentConstructorParameters
     module procedure sidmVelocityDependentConstructorInternal
  end interface darkMatterParticleSIDMVelocityDependent

contains

  function sidmVelocityDependentConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``\mono{selfInteractingDarkMatter}'' dark matter particle class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterParticleSIDMVelocityDependent)                    :: self
    type            (inputParameters                            ), intent(inout) :: parameters
    class           (darkMatterParticleClass                    ), pointer       :: darkMatterParticle_
    double precision                                                             :: sigma0             , velocityCharacteristic

    !![
    <inputParameter>
      <name>velocityCharacteristic</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>The characteristic velocity scale, $v_\mathrm{c}$, above which the self-interaction becomes velocity-suppressed. The total cross section is $\sigma(v) = \sigma_0 \, v_\mathrm{c}^2/(v_\mathrm{c}^2+v^2)$, which is approximately constant for $v \ll v_\mathrm{c}$ and falls as $v^{-2}$ for $v \gg v_\mathrm{c}$.</description>
    </inputParameter>
    <inputParameter>
      <name>sigma0</name>
      <source>parameters</source>
      <defaultValue>2.4d4</defaultValue>
      <description>Below the characteristic velocity the scattering is roughly isotropic with $\sigma$ $\approx$ sigma0.</description>
    </inputParameter>
    <objectBuilder class="darkMatterParticle"  name="darkMatterParticle_"  source="parameters"/>
    !!]
    self=darkMatterParticleSIDMVelocityDependent(velocityCharacteristic,sigma0,darkMatterParticle_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterParticle_"/>
    !!]
    return
  end function sidmVelocityDependentConstructorParameters

  function sidmVelocityDependentConstructorInternal(velocityCharacteristic,sigma0,darkMatterParticle_) result(self)
    !!{
    Internal constructor for the ``\mono{selfInteractingDarkMatter}'' dark matter particle class.
    !!}
    implicit none
    type            (darkMatterParticleSIDMVelocityDependent)                        :: self
    class           (darkMatterParticleClass                ), intent(in   ), target :: darkMatterParticle_
    double precision                                         , intent(in   )         :: velocityCharacteristic, sigma0
    !![
    <constructorAssign variables="sigma0,*darkMatterParticle_"/>
    !!]
    self%velocityCharacteristic_=velocityCharacteristic
    return
  end function sidmVelocityDependentConstructorInternal

  subroutine sidmVelocityDependentDestructor(self)
    !!{
    Destructor for the \mono{selfInteractingDarkMatter} dark matter particle class.
    !!}
    implicit none
    type(darkMatterParticleSIDMVelocityDependent), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterParticle_" />
    !!]
    return
  end subroutine sidmVelocityDependentDestructor

  double precision function sidmVelocityDependentMass(self)
    !!{
    Return the mass, in units of keV, of a self-interacting dark matter particle.
    !!}
    implicit none
    class(darkMatterParticleSIDMVelocityDependent), intent(inout) :: self

    sidmVelocityDependentMass=self%darkMatterParticle_%mass()
    return
  end function sidmVelocityDependentMass

  double precision function sidmVelocityDependentVelocityCharacteristic(self)
    !!{
    Return the characteristic velocity scale, $v_\mathrm{c}$, in units of km s$^{-1}$, of the velocity-dependent
    self-interaction cross section.
    !!}
    implicit none
    class(darkMatterParticleSIDMVelocityDependent), intent(inout) :: self

    sidmVelocityDependentVelocityCharacteristic=self%velocityCharacteristic_
    return
  end function sidmVelocityDependentVelocityCharacteristic

  double precision function sidmVelocityDependentCrossSectionSelfInteraction(self,velocityRelative) result(crossSection)
    !!{
    Return the self-interaction cross section, in units of cm$^2$ g$^{-1}$, of a self-interacting dark matter particle.
    !!}
    implicit none
    class(darkMatterParticleSIDMVelocityDependent), intent(inout) :: self
    double precision                              , intent(in   ) :: velocityRelative

    crossSection=+  self%sigma0                    &
         &       *  self%velocityCharacteristic_**2 &
         &       /(                                &
         &         +self%velocityCharacteristic_**2 &
         &         +     velocityRelative      **2 &
         &        )
    return
  end function sidmVelocityDependentCrossSectionSelfInteraction

  double precision function sidmVelocityDependentCrossSectionSelfInteractionDifferential(self,theta,velocityRelative) result(crossSection)
    !!{
    Return the differential self-interaction cross section, $\mathrm{d}\sigma/\mathrm{d}\theta$, in units of cm$^2$ g$^{-1}$
    ster$^{-1}$, of a self-interacting dark matter particle.
    !!}
    implicit none
    class           (darkMatterParticleSIDMVelocityDependent), intent(inout) :: self
    double precision                                         , intent(in   ) :: velocityRelative
    double precision                                         , intent(in   ) :: theta

    ! Anisotropic (forward-peaked) differential cross section for the velocity-dependent model.
    crossSection=+  self%velocityCharacteristic_**4 &
         &       *  self%sigma0                    &
         &       *0.5d0                            &
         &       *         sin(theta)              &
         &       /(                                &
         &         +self%velocityCharacteristic_**2 &
         &         +0.5d0                          &
         &         *     velocityRelative      **2 &
         &         *(1.0d0-cos(theta))             &
         &        )**2
    return
  end function sidmVelocityDependentCrossSectionSelfInteractionDifferential

  double precision function sidmVelocityDependentCrossSectionSelfInteractionDifferentialCos(self,cosTheta,velocityRelative) result(crossSection)
    !!{
    Return the differential self-interaction cross section, $\mathrm{d}\sigma/\mathrm{d}\cos\theta$, as a function of $\cos\theta$, in units of cm$^2$ g$^{-1}$, of a self-interacting dark matter particle.
    !!}
    implicit none
    class           (darkMatterParticleSIDMVelocityDependent), intent(inout) :: self
    double precision                                         , intent(in   ) :: velocityRelative
    double precision                                         , intent(in   ) :: cosTheta

    ! Anisotropic (forward-peaked) differential cross section for the velocity-dependent model.
    crossSection=+  self%velocityCharacteristic_**4 &
         &       *  self%sigma0                    &
         &       *0.5d0                            &
         &       /(                                &
         &         +self%velocityCharacteristic_**2 &
         &         + 0.5d0                         &
         &         *     velocityRelative      **2 &
         &         *(1.0d0-cosTheta)               &
         &        )**2
    return
  end function sidmVelocityDependentCrossSectionSelfInteractionDifferentialCos

  double precision function sidmVelocityDependentCrossSectionMomentumTransfer(self,velocityRelative) result(crossSection)
    !!{
    Return the momentum transfer self-interaction cross section, in units of cm$^2$ g$^{-1}$, of a self-interacting dark matter particle.
    !!}
    implicit none
    class(darkMatterParticleSIDMVelocityDependent), intent(inout) :: self
    double precision                              , intent(in   ) :: velocityRelative

    crossSection=+2.0d0                                      &
         &       *self%sigma0                                &
         &       *(                                          &
         &         +       self%velocityCharacteristic_       &
         &         /            velocityRelative             &
         &        )**4                                       &
         &       *(                                          &
         &         +log(                                     &
         &              +1.0d0                               &
         &              +(                                   &
         &                +     velocityRelative             &
         &                /self%velocityCharacteristic_       &
         &               )**2                                &
         &             )                                     &
         &         -   (                                     &
         &              +          velocityRelative      **2 &
         &              /(                                   &
         &                +        velocityRelative      **2 &
         &                +   self%velocityCharacteristic_**2 &
         &               )                                   &
         &             )                                     &
         &        )
    return
  end function sidmVelocityDependentCrossSectionMomentumTransfer

  double precision function sidmVelocityDependentCrossSectionViscosity(self,velocityRelative) result(crossSection)
    !!{
    Return the viscosity self-interaction cross section, in units of cm$^2$ g$^{-1}$, of a self-interacting dark matter particle.
    !!}
    implicit none
    class(darkMatterParticleSIDMVelocityDependent), intent(inout) :: self
    double precision                              , intent(in   ) :: velocityRelative

    ! This is the analytic integral (3/2) * integral[ (1-cos²(θ)) dσ/dcos(θ) ], i.e. the viscosity cross section in the
    ! convention normalized such that isotropic scattering gives σ_V = σ. Note that this 3/2 normalization differs from the
    ! bare (1-cos²(theta)) weight used in the effective cross section integrand (see dark_matter.particle.self_interacting.F90).
    crossSection=+6.0d0                                                                &
         &       *self%sigma0                                                          &
         &       /                  (velocityRelative/self%velocityCharacteristic_)**4  &
         &       *(                                                                    &
         &         +   (1.0d0+2.0d0/(velocityRelative/self%velocityCharacteristic_)**2) &
         &         *log(1.0d0+      (velocityRelative/self%velocityCharacteristic_)**2) &
         &         -2.0d0                                                              &
         &        )
    return
  end function sidmVelocityDependentCrossSectionViscosity

