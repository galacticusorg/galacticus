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
Tests of the self-interacting dark matter particle cross sections. For the velocity-dependent model the total, differential,
momentum-transfer, and viscosity cross sections are checked against their closed forms, and the effective cross section against
the analytic result for this model. For the constant model the differential and effective cross sections are checked. The
velocity-dependent model is $\sigma(v)=\sigma_0 w^2/(w^2+v^2)$ with characteristic velocity $w$; the analytic effective cross
section uses an effective velocity dispersion $v_\mathrm{eff}=0.64\,V_\mathrm{max}$.
!!}

program Tests_SIDM_Cross_Sections
  use :: Dark_Matter_Particles   , only : darkMatterParticleCDM, darkMatterParticleSIDMVelocityDependent, darkMatterParticleSelfInteractingDarkMatterConstant
  use :: Display                 , only : displayVerbositySet   , verbosityLevelStandard
  use :: Numerical_Constants_Math, only : Pi
  use :: Unit_Tests              , only : Assert                , Unit_Tests_Begin_Group                 , Unit_Tests_End_Group                              , Unit_Tests_Finish
  implicit none
  type            (darkMatterParticleCDM                             )            :: darkMatterParticleCDM_
  type            (darkMatterParticleSIDMVelocityDependent           )            :: particleVelocityDependent
  type            (darkMatterParticleSelfInteractingDarkMatterConstant)           :: particleConstant
  ! Velocity-dependent model parameters (sigma0 in cm² g⁻¹, characteristic velocity in km s⁻¹) and the constant-model cross
  ! section.
  double precision                                                    , parameter :: sigma0              =10.0d0, velocityCharacteristic=50.0d0, &
       &                                                                             crossSectionConstant=10.0d0

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("SIDM cross sections")
  ! Construct the dark matter particles.
  darkMatterParticleCDM_   =darkMatterParticleCDM                              (                                                                                                       )
  particleVelocityDependent=darkMatterParticleSIDMVelocityDependent           (velocityCharacteristic     =velocityCharacteristic,sigma0=sigma0,darkMatterParticle_=darkMatterParticleCDM_)
  particleConstant         =darkMatterParticleSelfInteractingDarkMatterConstant(crossSectionSelfInteraction=crossSectionConstant                ,darkMatterParticle_=darkMatterParticleCDM_)

  ! Velocity-dependent total cross section: sigma(v) = sigma0 w²/(w²+v²). Closed form, so the tolerance is tight.
  call Assert("velocity-dependent total cross section, v=10" ,particleVelocityDependent%crossSectionSelfInteraction( 10.0d0),9.6153846154d+00,relTol=1.0d-6)
  call Assert("velocity-dependent total cross section, v=50" ,particleVelocityDependent%crossSectionSelfInteraction( 50.0d0),5.0000000000d+00,relTol=1.0d-6)
  call Assert("velocity-dependent total cross section, v=200",particleVelocityDependent%crossSectionSelfInteraction(200.0d0),5.8823529412d-01,relTol=1.0d-6)
  ! Velocity-dependent differential cross sections. Note dsigma/dtheta = sin(theta) * dsigma/dcos(theta).
  call Assert("velocity-dependent dsigma/dtheta, theta=pi/2, v=50"      ,particleVelocityDependent%crossSectionSelfInteractionDifferential   (0.5d0*Pi, 50.0d0),2.2222222222d+00,relTol=1.0d-6)
  call Assert("velocity-dependent dsigma/dcostheta, costheta=0, v=50"   ,particleVelocityDependent%crossSectionSelfInteractionDifferentialCos( 0.0d0  , 50.0d0),2.2222222222d+00,relTol=1.0d-6)
  call Assert("velocity-dependent dsigma/dcostheta, costheta=-0.5, v=100",particleVelocityDependent%crossSectionSelfInteractionDifferentialCos(-0.5d0  ,100.0d0),3.1250000000d-01,relTol=1.0d-6)
  ! Velocity-dependent momentum-transfer cross section.
  call Assert("velocity-dependent momentum-transfer, v=10" ,particleVelocityDependent%crossSectionSelfInteractionMomentumTransfer( 10.0d0),9.4896836468d+00,relTol=1.0d-6)
  call Assert("velocity-dependent momentum-transfer, v=50" ,particleVelocityDependent%crossSectionSelfInteractionMomentumTransfer( 50.0d0),3.8629436112d+00,relTol=1.0d-6)
  call Assert("velocity-dependent momentum-transfer, v=200",particleVelocityDependent%crossSectionSelfInteractionMomentumTransfer(200.0d0),1.4781538074d-01,relTol=1.0d-6)
  ! Velocity-dependent viscosity cross section.
  call Assert("velocity-dependent viscosity, v=50" ,particleVelocityDependent%crossSectionSelfInteractionViscosity( 50.0d0),4.7664925008d+00,relTol=1.0d-6)
  call Assert("velocity-dependent viscosity, v=200",particleVelocityDependent%crossSectionSelfInteractionViscosity(200.0d0),2.7828867470d-01,relTol=1.0d-6)
  ! Velocity-dependent effective cross section, compared against the analytic result. The tolerance is looser here as the code
  ! evaluates this by 2-D integration (1% tolerance) followed by interpolation in a table.
  call Assert("velocity-dependent effective cross section, Vmax=20" ,particleVelocityDependent%crossSectionEffective( 20.0d0),4.9568318720d+00,relTol=3.0d-2)
  call Assert("velocity-dependent effective cross section, Vmax=80" ,particleVelocityDependent%crossSectionEffective( 80.0d0),3.9782568018d-01,relTol=3.0d-2)
  call Assert("velocity-dependent effective cross section, Vmax=320",particleVelocityDependent%crossSectionEffective(320.0d0),6.1145891038d-03,relTol=3.0d-2)
  ! Velocity-dependent effective cross section low-velocity limit: for Vmax << w the scattering is velocity-independent and
  ! sigma_eff -> sigma0.
  call Assert("velocity-dependent effective cross section, low-velocity limit",particleVelocityDependent%crossSectionEffective(1.0d0),sigma0,relTol=2.0d-2)
  ! Constant model: the differential cross section is sigma/2 per unit cos(theta), and the effective cross section equals the
  ! cross section itself.
  call Assert("constant total cross section"                ,particleConstant%crossSectionSelfInteraction          ( 50.0d0)        ,crossSectionConstant      ,relTol=1.0d-9)
  call Assert("constant dsigma/dcostheta = sigma/2"         ,particleConstant%crossSectionSelfInteractionDifferentialCos(0.3d0,50.0d0),crossSectionConstant*0.5d0,relTol=1.0d-9)
  call Assert("constant effective cross section = sigma"    ,particleConstant%crossSectionEffective                 (123.0d0)        ,crossSectionConstant      ,relTol=1.0d-9)

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_SIDM_Cross_Sections
