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

!!{
Implements a selfInteracting dark matter particle class.
!!}

  !![
  <darkMatterParticle name="darkMatterParticleSelfInteractingDarkMatter" abstract="yes">
    <description>
      Provides a self-interacting dark matter particle candidate in which dark matter undergoes elastic scattering, modifying halo
      density profiles on small scales. The elastic self-interaction cross section per unit mass in units of cm$^2$~g$^{-1}$ is
      set by \mono{[crossSectionSelfInteraction]}.
    </description>
  </darkMatterParticle>
  !!]
  type, abstract, extends(darkMatterParticleClass) :: darkMatterParticleSelfInteractingDarkMatter
     !!{
     A self-interacting dark matter particle class.
     !!}
     private
   contains
     !![
     <methods>
       <method method="crossSectionSelfInteraction"                 description="Return the self-interaction cross section, $\sigma$, of the dark matter particle in units of cm$^2$ g$^{-1}$."                                                                           />
       <method method="crossSectionSelfInteractionDifferential"     description="Return the differential self-interaction cross section, $\mathrm{d}\sigma/\mathrm{d}\Omega$, of the dark matter particle in units of cm$^2$ g$^{-1}$ ster$^{-1}$."                       />
       <method method="crossSectionSelfInteractionDifferentialCos"  description="Return the differential self-interaction cross section, $\mathrm{d}\sigma/\mathrm{d}\cos\theta$, of the dark matter particle as a function of $\cos\theta$, in units of cm$^2$ g$^{-1}$."/>
       <method method="crossSectionSelfInteractionMomentumTransfer" description="Return the momentum transfer self-interaction cross section, $\sigma$, of the dark matter particle in units of cm$^2$ g$^{-1}$."                                                         />
       <method method="crossSectionSelfInteractionViscosity"        description="Return the viscosity self-interaction cross section, $\sigma$, of the dark matter particle in units of cm$^2$ g$^{-1}$."                                                                 />
     </methods>
     !!]
     procedure(crossSectionSelfInteractionTemplate                ), deferred :: crossSectionSelfInteraction
     procedure(crossSectionSelfInteractionDifferentialTemplate    ), deferred :: crossSectionSelfInteractionDifferential
     procedure(crossSectionSelfInteractionDifferentialCosTemplate ), deferred :: crossSectionSelfInteractionDifferentialCos
     procedure(crossSectionSelfInteractionMomentumTransferTemplate), deferred :: crossSectionSelfInteractionMomentumTransfer
     procedure(crossSectionSelfInteractionViscosityTemplate       ), deferred :: crossSectionSelfInteractionViscosity
     procedure                                                                :: crossSectionEffective                       => crossSectionEffective_
  end type darkMatterParticleSelfInteractingDarkMatter

  abstract interface
     double precision function crossSectionSelfInteractionTemplate(self,velocityRelative)
       !!{
       Interface for self-interaction cross section, in units of cm$^2$ g$^{-1}$, of a self-interacting dark matter particle.
       !!}
       import darkMatterParticleSelfInteractingDarkMatter
       class(darkMatterParticleSelfInteractingDarkMatter), intent(inout) :: self
       double precision                                  , intent(in   ) :: velocityRelative
     end function crossSectionSelfInteractionTemplate
  end interface

  abstract interface
     double precision function crossSectionSelfInteractionDifferentialTemplate(self,theta,velocityRelative)
       !!{
       Interface for differential self-interaction cross section, $\mathrm{d}\sigma/\mathrm{d}\theta$, in units of cm$^2$ g$^{-1}$ ster$^{-1}$, of a self-interacting dark matter particle.
       !!}
       import darkMatterParticleSelfInteractingDarkMatter
       class (darkMatterParticleSelfInteractingDarkMatter), intent(inout) :: self
       double precision                                   , intent(in)    :: velocityRelative
       double precision                                   , intent(in   ) :: theta
     end function crossSectionSelfInteractionDifferentialTemplate

     double precision function crossSectionSelfInteractionDifferentialCosTemplate(self,Costheta,velocityRelative)
       !!{
       Interface for differential self-interaction cross section, $\mathrm{d}\sigma/\mathrm{d}\cos\theta$, as a function of $\cos\theta$, in units of cm$^2$ g$^{-1}$, of a self-interacting dark matter particle.
       !!}
       import darkMatterParticleSelfInteractingDarkMatter
       class (darkMatterParticleSelfInteractingDarkMatter), intent(inout) :: self
       double precision                                   , intent(in)    :: velocityRelative
       double precision                                   , intent(in   ) :: Costheta
     end function crossSectionSelfInteractionDifferentialCosTemplate

     double precision function crossSectionSelfInteractionMomentumTransferTemplate(self,velocityRelative)
       !!{
       Interface for momentum transfer self-interaction cross section, in units of cm$^2$ g$^{-1}$, of a self-interacting dark matter particle.
       !!}
       import darkMatterParticleSelfInteractingDarkMatter
       class(darkMatterParticleSelfInteractingDarkMatter), intent(inout) :: self
       double precision                                  , intent(in   ) :: velocityRelative
     end function crossSectionSelfInteractionMomentumTransferTemplate

     double precision function crossSectionSelfInteractionViscosityTemplate(self,velocityRelative)
       !!{
       Interface for viscosity self-interaction cross section, in units of cm$^2$ g$^{-1}$, of a self-interacting dark matter particle.
       !!}
       import darkMatterParticleSelfInteractingDarkMatter
       class(darkMatterParticleSelfInteractingDarkMatter), intent(inout) :: self
       double precision                                  , intent(in   ) :: velocityRelative
     end function crossSectionSelfInteractionViscosityTemplate
  end interface
   
  double precision :: velocityEffective_
  !$omp threadprivate (velocityEffective_)

  class (darkMatterParticleSelfInteractingDarkMatter) , pointer :: self_
  !$omp threadprivate (self_)

contains

  double precision function integrandNumerator(velocity,cosTheta)
    !!{
    Integrand of the numerator of the effective cross section, $\sigma_\mathrm{eff}$, of \cite{yang_parametric_2024} (their
    eqn.~1.1). The integrand is $(\mathrm{d}\sigma/\mathrm{d}\cos\theta)\,\sin^2\theta\,v^5\,f_\mathrm{MB}(v)$, where the kinetic-theory
    conductivity kernel contributes $\sin^2\theta\,v^5$ and the Maxwell--Boltzmann relative-velocity distribution contributes
    $f_\mathrm{MB}(v) \propto v^2 \exp[-v^2/(4 v_\mathrm{eff}^2)]$. The two factors of $v$ combine to give the $v^7$ below, and
    $\sin^2\theta = 1-\cos^2\theta$.
    !!}
    double precision, intent(in   ) :: velocity, cosTheta

    integrandNumerator=+self_%crossSectionSelfInteractionDifferentialCos(cosTheta,velocity) &
         &             *velocity**7                                                         &
         &             *(1.0d0-cosTheta**2)                                                 &
         &             *exp(                                                                &
         &                  +velocity          **2                                          &
         &                  /4.0d0                                                          &
         &                  /velocityEffective_**2                                          &
         &              )
    return
  end function integrandNumerator

  double precision function crossSectionEffective_(self,velocityMaximum)
    !!{
    Evaluate the effective cross section.
    !!}
    use :: Numerical_Integration_2D, only : integrator2D
    class           (darkMatterParticleSelfInteractingDarkMatter), intent(inout) , target  :: self
    double precision                                             , intent(in   )           :: velocityMaximum
    double precision                                             , parameter               :: factorVelocityEffective=0.64d0
    type            (integrator2D)                               , allocatable             :: integratorNumerator
    double precision                                             , dimension(2,2)          :: boundaries
    double precision                                                                       :: velocityEffective
    double precision                                                                       :: numeratorIntegral

    ! Effective velocity dispersion of the Maxwell-Boltzmann weighting: ν_eff = 0.64*Vmax for an NFW halo
    ! (Yang et al. 2024; JCAP; 2; 32).
    velocityEffective=factorVelocityEffective*velocityMaximum
    ! Set integration boundaries.
    boundaries(1,:)=[+0.0d0,+10.0d0*velocityEffective]
    boundaries(2,:)=[-1.0d0,+ 1.0d0                  ]
    ! Set sub-module scope copies.
    velocityEffective_ =  velocityEffective
    self_              => self
    ! Build an integrator.
    allocate(integratorNumerator)
    call integratorNumerator%setIntegrand(integrandNumerator)
    ! Compute the numerator of the expression.
    numeratorIntegral=integratorNumerator%integrate(boundaries)
    ! Normalize by the (cross-section-independent) denominator of eqn.~1.1 of Yang et al. (2024; JCAP; 2; 32), evaluated
    ! analytically: ½ * ∫ sin²(θ) v⁷ exp(-v²/4 Veff²) = ½ * (4/3) * 768 Veff⁸ = 512 Veff⁸, where the factor of ½ is the leading
    ! "2" in the numerator of that equation, (4/3) = integral over cos(θ) of sin²(θ), and 768 Veff⁸ = integral over v of v⁷
    ! exp(-v²/4 Veff²).
    crossSectionEffective_=numeratorIntegral/(512.0d0*velocityEffective_**8)
    return
  end function crossSectionEffective_

