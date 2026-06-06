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
   <description>Provides a self-interacting dark matter particle candidate in which dark matter undergoes elastic scattering, modifying halo density profiles on small scales. The elastic self-interaction cross section per unit mass in units of cm$^2$~g$^{-1}$ is set by \mono{[crossSectionSelfInteraction]}.</description>
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
       <method description="Return the self-interaction cross section, $\sigma$, of the dark matter particle in units of cm$^2$ g$^{-1}$." method="crossSectionSelfInteraction" />
       <method description="Return the differential self-interaction cross section, $\mathrm{d}\sigma/\mathrm{d}\Omega$, of the dark matter particle in units of cm$^2$ g$^{-1}$ ster$^{-1}$." method="crossSectionSelfInteractionDifferential"/>
       <method description="Return the differential self-interaction cross section, $\mathrm{d}\sigma/\mathrm{d}\Omega$, of the dark matter particle as a function of $\cos\theta$, in units of cm$^2$ g$^{-1}$ ster$^{-1}$." method="crossSectionSelfInteractionDifferentialCos"/>
       <method description="Return the momentum transfer self-interaction cross section, $\sigma$, of the dark matter particle in units of cm$^2$ g$^{-1}$." method="crossSectionSelfInteractionMomentumTransfer" />
       <method description="Return the viscosity self-interaction cross section, $\sigma$, of the dark matter particle in units of cm$^2$ g$^{-1}$." method="crossSectionSelfInteractionViscosity" />
     </methods>
     !!]
     procedure(crossSectionSelfInteractionTemplate                ), deferred :: crossSectionSelfInteraction
     procedure(crossSectionSelfInteractionDifferentialTemplate    ), deferred :: crossSectionSelfInteractionDifferential
     procedure(crossSectionSelfInteractionDifferentialCosTemplate ), deferred :: crossSectionSelfInteractionDifferentialCos
     procedure(crossSectionSelfInteractionMomentumTransferTemplate), deferred :: crossSectionSelfInteractionMomentumTransfer
     procedure(crossSectionSelfInteractionViscosityTemplate       ), deferred :: crossSectionSelfInteractionViscosity
     procedure                                                                :: effectiveCrossSection => effectiveCrossSection_
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
       Interface for differential self-interaction cross section, $\mathrm{d}\sigma/\mathrm{d}\theta$, in units of cm$^2$ g$^{-1}$ ster$^{-1}$, of a self-interacting dark matter particle.
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
   
  double precision :: Veff_
  !$omp threadprivate (Veff_)

  class (darkMatterParticleSelfInteractingDarkMatter) , pointer :: self_
  !$omp threadprivate (self_)

contains

  double precision function integrandNumerator(v,Costheta)
    double precision, intent(in   ) :: v, Costheta

    integrandNumerator = self_%crossSectionSelfInteractionDifferentialCos(Costheta,v)*v**7*(1.0d0-Costheta**2)*exp(-v**2/(4.0d0*Veff_**2))
    return
  end function integrandNumerator

  double precision function effectiveCrossSection_(self, Vmax)
    use :: SIDM_Effective_Cross_Section_Integrator, only : SIDMEffectiveCrossSectionIntegrator2D

    class           (darkMatterParticleSelfInteractingDarkMatter), intent(inout), target  :: self
    type            (SIDMEffectiveCrossSectionIntegrator2D)      , allocatable            :: integratorNumerator
    double precision                                             , intent(in   )          :: Vmax
    double precision                                             , dimension(2,2)         :: boundaries
    double precision                                                                      :: Veff
    double precision                                                                      :: numeratorIntegral

!    print *, 'begining of effective cross section function'
    
    Veff = 0.64d0*Vmax

    boundaries (1,:) = (/0.0d0,10.0d0*Veff/)
    boundaries (2,:) = (/-1.0d0,1.0d0/)
      
    Veff_=Veff
    self_ => self

    allocate(integratorNumerator)
    call integratorNumerator%setIntegrand(integrandNumerator)

!    print *, 'call intrgrator ...'

    numeratorIntegral = integratorNumerator%integrate(boundaries)

!    print *, 'numerator integral ...', numeratorIntegral

    effectiveCrossSection_ = numeratorIntegral/(512.0d0*Veff_**8)
    return
  end function effectiveCrossSection_

