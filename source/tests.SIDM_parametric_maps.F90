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
Tests of the parametric self-interacting dark matter model fitting functions of \cite{yang_parametric_2024} (the
\mono{SIDM\_Parametric\_Model} module): the $V_\mathrm{max}$/$R_\mathrm{max}$ maps between SIDM and CDM (NFW) halos and their
$\tau$-derivatives, the characteristic-density/scale-radius/core-radius profile fits, the NFW scale relations, and the
gravothermal evolution timescale. Reference values are evaluated directly from the published polynomial coefficients; the
gravothermal timescale (whose absolute value depends on several physical-unit conversions) is checked via its scalings.
!!}

program Tests_SIDM_Parametric_Maps
  use :: Dark_Matter_Particles, only : darkMatterParticleCDM, darkMatterParticleSIDMVelocityDependent
  use :: Display              , only : displayVerbositySet  , verbosityLevelStandard
  use :: SIDM_Parametric_Model, only : velocityMaximumNFW, radiusMaximumNFW , velocityMaximumRateTau, radiusMaximumRateTau, &
       &                               densityScale      , radiusScale      , radiusCore            , radiusScaleNFW      , &
       &                               densityScaleNFW   , timescaleCollapse
  use :: Unit_Tests           , only : Assert               , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  type            (darkMatterParticleCDM                  ) :: darkMatterParticleCDM_
  type            (darkMatterParticleSIDMVelocityDependent) :: particle
  double precision                                         :: timescaleC1, timescaleC2, densityV1, densityV2, densityR2

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("SIDM parametric-model fitting functions")

  ! Inverse Vmax/Rmax maps (SIDM -> CDM), eqn. 2.4 of Yang et al. (2024).
  call Assert("velocityMaximumNFW(100,0.3)",velocityMaximumNFW(100.0d0,0.3d0), 9.7507123889d+01,relTol=1.0d-6)
  call Assert("velocityMaximumNFW(100,0.7)",velocityMaximumNFW(100.0d0,0.7d0), 9.1462296916d+01,relTol=1.0d-6)
  call Assert("radiusMaximumNFW( 10,0.3)",radiusMaximumNFW( 10.0d0,0.3d0), 1.0576544915d+01,relTol=1.0d-6)
  call Assert("radiusMaximumNFW( 10,0.7)",radiusMaximumNFW( 10.0d0,0.7d0), 1.3599537017d+01,relTol=1.0d-6)
  ! Derivatives dVmax/dtau and dRmax/dtau, eqn. 2.4.
  call Assert("velocityMaximumRateTau(0.3,100)",velocityMaximumRateTau(0.3d0,100.0d0), 6.9640078925d+00,relTol=1.0d-6)
  call Assert("velocityMaximumRateTau(0.7,100)",velocityMaximumRateTau(0.7d0,100.0d0), 1.2364796992d+01,relTol=1.0d-6)
  call Assert("radiusMaximumRateTau(0.3, 10)",radiusMaximumRateTau(0.3d0, 10.0d0),-3.4806604506d+00,relTol=1.0d-6)
  call Assert("radiusMaximumRateTau(0.7, 10)",radiusMaximumRateTau(0.7d0, 10.0d0),-6.9276935114d+00,relTol=1.0d-6)
  ! Profile-parameter fits rho_s, r_s, r_c, eqn. 2.3.
  call Assert("densityScale(1e7,0.3)",densityScale(1.0d7,0.3d0),2.0899456893d+07,relTol=1.0d-6)
  call Assert("densityScale(1e7,0.7)",densityScale(1.0d7,0.7d0),3.0691418253d+07,relTol=1.0d-6)
  call Assert("radiusScale(0.01,0.3)",radiusScale(0.01d0,0.3d0),7.4732221165d-03,relTol=1.0d-6)
  call Assert("radiusScale(0.01,0.7)",radiusScale(0.01d0,0.7d0),6.4180909178d-03,relTol=1.0d-6)
  call Assert("radiusCore(0.01,0.3)",radiusCore(0.01d0,0.3d0),4.6720936443d-03,relTol=1.0d-6)
  call Assert("radiusCore(0.01,0.7)",radiusCore(0.01d0,0.7d0),2.6655019779d-03,relTol=1.0d-6)
  ! NFW scale radius relation.
  call Assert("radiusScaleNFW(0.02)",radiusScaleNFW(0.02d0),9.2482059958d-03,relTol=1.0d-6)
  ! tau=0 (CDM) limit: the maps are the identity and the core radius vanishes.
  call Assert("velocityMaximumNFW tau=0 is the identity" ,velocityMaximumNFW (100.0d0,0.0d0),100.0d0,relTol=1.0d-9)
  call Assert("radiusMaximumNFW tau=0 is the identity" ,radiusMaximumNFW ( 10.0d0,0.0d0), 10.0d0,relTol=1.0d-9)
  call Assert("densityScale tau=0 is the identity",densityScale(1.0d7  ,0.0d0),1.0d7  ,relTol=1.0d-6)
  call Assert("radiusScale tau=0 is the identity"  ,radiusScale  (0.01d0 ,0.0d0),0.01d0 ,relTol=1.0d-6)
  call Assert("radiusCore tau=0 vanishes"         ,radiusCore  (0.01d0 ,0.0d0),0.0d0  ,absTol=1.0d-9)
  ! tau>1 (core-collapsed) limit: the velocity and radius derivatives vanish.
  call Assert("velocityMaximumRateTau tau>1 vanishes",velocityMaximumRateTau(1.5d0,100.0d0),0.0d0,absTol=1.0d-12)
  call Assert("radiusMaximumRateTau tau>1 vanishes",radiusMaximumRateTau(1.5d0, 10.0d0),0.0d0,absTol=1.0d-12)
  ! tau is clamped to [0,1]: evaluating beyond 1 returns the tau=1 result.
  call Assert("velocityMaximumNFW tau is clamped to [0,1]",velocityMaximumNFW(100.0d0,1.5d0),velocityMaximumNFW(100.0d0,1.0d0),relTol=1.0d-12)

  ! Characteristic density densityScaleNFW: proportional to Vmax^2 and to 1/Rs^2.
  densityV1=densityScaleNFW(0.01d0, 50.0d0)
  densityV2=densityScaleNFW(0.01d0,100.0d0)
  call Assert("densityScaleNFW proportional to Vmax^2",densityV2/densityV1,4.0d0,relTol=1.0d-9)
  densityR2=densityScaleNFW(0.02d0, 50.0d0)
  call Assert("densityScaleNFW proportional to 1/Rs^2",densityR2/densityV1,0.25d0,relTol=1.0d-9)

  ! Gravothermal timescale timescaleCollapse: build a velocity-dependent particle and check that t_c is inversely proportional to the
  ! calibration constant C (the absolute value depends on physical-unit conversions and is exercised end-to-end elsewhere).
  darkMatterParticleCDM_=darkMatterParticleCDM                  (                                                                                )
  particle              =darkMatterParticleSIDMVelocityDependent(velocityCharacteristic=50.0d0,sigma0=10.0d0,darkMatterParticle_=darkMatterParticleCDM_)
  timescaleC1=timescaleCollapse(particle,0.75d0,80.0d0,0.02d0,80.0d0)
  timescaleC2=timescaleCollapse(particle,1.50d0,80.0d0,0.02d0,80.0d0)
  call Assert("timescaleCollapse proportional to 1/C",timescaleC1/timescaleC2,2.0d0,relTol=1.0d-6)

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_SIDM_Parametric_Maps
