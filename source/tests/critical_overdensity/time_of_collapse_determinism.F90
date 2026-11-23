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

!+    Contributions to this file made by: Andrew Benson, Claude.

program Tests_Critical_Overdensity_Time_Of_Collapse_Determinism
  !!{RST
  Tests that the time of collapse returned by a critical overdensity class (from its tabulation of collapse time against
  critical overdensity) is independent of the history of requests made of it. Two otherwise identical critical overdensity
  objects are asked for the collapse times of the same set of critical overdensities, in opposite orders, each value being
  recorded as soon as it is requested. Since each request extends the tabulation (to larger or smaller critical overdensities),
  values recorded before an extension are compared with values recorded after it, and the results are required to agree
  *exactly*. Each collapse time is also checked against the defining condition, :math:`\delta_\mathrm{c}(t)/D(t)=\delta`.

  The same tests are then applied to a mass-dependent critical overdensity, for which collapse times are not tabulated, but are
  found by a root find for each request. There the defining condition is
  :math:`\delta_\mathrm{c}(t,M)\sigma(M,t_0)/\sigma(M,t)=\delta`.

  Finally, collapse times are found (by both the tabulated and the mass-dependent route) for critical overdensities in an
  overdense environment, using the peak-background split. There the effective critical overdensity at the present epoch is
  negative, which must be handled when constructing initial guesses for collapse times.
  !!}
  use :: Cosmological_Density_Field          , only : cosmologicalMassVarianceFilteredPower   , criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt, &
       &                                              criticalOverdensityBarkana2001WDM       , criticalOverdensityPeakBackgroundSplit                      , &
       &                                              haloEnvironmentFixed
  use :: Cosmology_Functions                 , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters                , only : cosmologyParametersSimple
  use :: Dark_Matter_Particles               , only : darkMatterParticleCDM                   , darkMatterParticleWDMThermal
  use :: Display                             , only : displayVerbositySet                     , verbosityLevelStandard
  use :: Events_Hooks                        , only : eventsHooksInitialize
  use :: Functions_Global_Utilities          , only : Functions_Global_Set
  use :: Linear_Growth                       , only : linearGrowthCollisionlessMatter
  use :: Power_Spectra_Primordial            , only : powerSpectrumPrimordialPowerLaw
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredSimple
  use :: Power_Spectrum_Normalizers          , only : powerSpectrumNormalizerSigma8
  use :: Power_Spectrum_Window_Functions     , only : powerSpectrumWindowFunctionTopHat
  use :: Transfer_Functions                  , only : transferFunctionIdentity
  use :: Unit_Tests                          , only : Assert                                  , Unit_Tests_Begin_Group                                      , &
       &                                              Unit_Tests_End_Group                    , Unit_Tests_Finish
  implicit none
  type            (cosmologyParametersSimple                                   )                          :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                              )                          :: cosmologyFunctions_
  type            (linearGrowthCollisionlessMatter                             )                          :: linearGrowth_
  type            (darkMatterParticleCDM                                       )                          :: darkMatterParticle_
  type            (darkMatterParticleWDMThermal                                )                          :: darkMatterParticleWDM_
  type            (powerSpectrumPrimordialPowerLaw                             )                          :: powerSpectrumPrimordial_
  type            (transferFunctionIdentity                                    )                          :: transferFunction_
  type            (powerSpectrumNormalizerSigma8                               )                          :: powerSpectrumNormalizer_
  type            (powerSpectrumPrimordialTransferredSimple                    )                          :: powerSpectrumPrimordialTransferred_
  type            (powerSpectrumWindowFunctionTopHat                           )                          :: powerSpectrumWindowFunction_
  type            (cosmologicalMassVarianceFilteredPower                       )                          :: cosmologicalMassVariance_
  type            (criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt)                          :: criticalOverdensityAscending_      , criticalOverdensityDescending_
  type            (criticalOverdensityBarkana2001WDM                           )                          :: criticalOverdensityMassAscending_  , criticalOverdensityMassDescending_
  type            (haloEnvironmentFixed                                        )                          :: haloEnvironment_
  type            (criticalOverdensityPeakBackgroundSplit                      )                          :: criticalOverdensityEnvironment_    , criticalOverdensityMassEnvironment_
  ! Linear overdensity (at the present epoch) of the environment, chosen to exceed the critical overdensity for collapse at the
  ! present epoch, so that the effective critical overdensity at that epoch is negative.
  double precision                                                              , parameter               :: overdensityEnvironment  =2.0d0
  ! Critical overdensity and mass for requests made in the overdense environment.
  double precision                                                              , parameter               :: overdensityInEnvironment=1.00003d0, massInEnvironment=4.71d11
  double precision                                                                                        :: time
  ! Critical overdensities at which collapse times are requested. These are deliberately not points of the lattice (of 10⁴
  ! points per unit, by default) on which collapse times are tabulated, so that interpolation in the tabulation is tested, and
  ! span collapse times from the present day to early times.
  double precision                                                              , dimension(6), parameter :: overdensities           =[1.75313d0,2.61777d0,6.10331d0,15.27193d0,41.68537d0,1.69991d0]
  ! Halo masses (in M☉) for the requests made of the mass-dependent critical overdensity. These span the Jeans mass of the warm
  ! dark matter particle, so that the critical overdensity is significantly mass-dependent.
  double precision                                                              , dimension(6), parameter :: masses                  =[3.17d10,1.23d8,4.71d11,7.29d9,2.58d12,8.66d8]
  double precision                                                              , dimension(6)            :: timesAscending                      , timesDescending                   , &
       &                                                                                                     timesRepeat
  logical                                                                                                 :: satisfiesCondition
  integer                                                                                                 :: i

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Initialize event hooks.
  call eventsHooksInitialize()
  call Functions_Global_Set ()
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Critical overdensity: determinism of collapse times")
  ! Construct the objects needed.
  cosmologyParameters_               =cosmologyParametersSimple                                   (                                                                              &
       &                                                                                           OmegaMatter                = 0.2750d0                                       , &
       &                                                                                           OmegaBaryon                = 0.0000d0                                       , &
       &                                                                                           OmegaDarkEnergy            = 0.7250d0                                       , &
       &                                                                                           temperatureCMB             = 2.7800d0                                       , &
       &                                                                                           HubbleConstant             =70.2000d0                                         &
       &                                                                                          )
  cosmologyFunctions_                =cosmologyFunctionsMatterLambda                              (                                                                              &
       &                                                                                           cosmologyParameters_       =cosmologyParameters_                              &
       &                                                                                          )
  linearGrowth_                      =linearGrowthCollisionlessMatter                             (                                                                              &
       &                                                                                           cosmologyParameters_       =cosmologyParameters_                            , &
       &                                                                                           cosmologyFunctions_        =cosmologyFunctions_                               &
       &                                                                                          )
  darkMatterParticle_                =darkMatterParticleCDM                                       (                                                                              &
       &                                                                                          )
  powerSpectrumPrimordial_           =powerSpectrumPrimordialPowerLaw                             (                                                                              &
       &                                                                                           index_                     =-1.0d0                                          , &
       &                                                                                           running                    =+0.0d0                                          , &
       &                                                                                           runningRunning             =+0.0d0                                          , &
       &                                                                                           wavenumberReference        =+1.0d0                                          , &
       &                                                                                           runningSmallScalesOnly     =.false.                                           &
       &                                                                                          )
  transferFunction_                  =transferFunctionIdentity                                    (                                                                              &
       &                                                                                           cosmologyParameters_       =cosmologyParameters_                            , &
       &                                                                                           time                       =13.8d0                                            &
       &                                                                                          )
  powerSpectrumNormalizer_           =powerSpectrumNormalizerSigma8                               (                                                                              &
       &                                                                                           sigma8                     =0.8d+0                                          , &
       &                                                                                           tolerance                  =1.0d-4                                          , &
       &                                                                                           integrationFailureIsFatal  =.true.                                          , &
       &                                                                                           cosmologyParameters_       =cosmologyParameters_                            , &
       &                                                                                           cosmologyFunctions_        =cosmologyFunctions_                             , &
       &                                                                                           transferFunction_          =transferFunction_                                 &
       &                                                                                          )
  powerSpectrumPrimordialTransferred_=powerSpectrumPrimordialTransferredSimple                    (                                                                              &
       &                                                                                           powerSpectrumPrimordial_   =powerSpectrumPrimordial_                        , &
       &                                                                                           transferFunction_          =transferFunction_                               , &
       &                                                                                           linearGrowth_              =linearGrowth_                                   , &
       &                                                                                           powerSpectrumNormalizer_   =powerSpectrumNormalizer_                          &
       &                                                                                          )
  powerSpectrumWindowFunction_       =powerSpectrumWindowFunctionTopHat                           (                                                                              &
       &                                                                                           cosmologyParameters_       =cosmologyParameters_                              &
       &                                                                                          )
  cosmologicalMassVariance_          =cosmologicalMassVarianceFilteredPower                       (                                                                              &
       &                                                                                           tolerance                               =1.0d-4                             , &
       &                                                                                           toleranceTopHat                         =1.0d-4                             , &
       &                                                                                           rootVarianceLogarithmicGradientTolerance=1.0d-9                             , &
       &                                                                                           integrationFailureIsFatal               =.true.                             , &
       &                                                                                           storeTabulations                        =.false.                            , &
       &                                                                                           nonMonotonicIsFatal                     =.true.                             , &
       &                                                                                           monotonicInterpolation                  =.false.                            , &
       &                                                                                           truncateAtParticleHorizon               =.false.                            , &
       &                                                                                           cosmologyParameters_                    =cosmologyParameters_               , &
       &                                                                                           cosmologyFunctions_                     =cosmologyFunctions_                , &
       &                                                                                           linearGrowth_                           =linearGrowth_                      , &
       &                                                                                           powerSpectrumPrimordialTransferred_     =powerSpectrumPrimordialTransferred_, &
       &                                                                                           powerSpectrumWindowFunction_            =powerSpectrumWindowFunction_         &
       &                                                                                          )
  ! Table storage is disabled so that the test neither depends upon, nor perturbs, any tabulations cached in the dynamic data
  ! path by other runs.
  criticalOverdensityAscending_      =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt(                                                                              &
       &                                                                                           linearGrowth_              =linearGrowth_                                   , &
       &                                                                                           cosmologyFunctions_        =cosmologyFunctions_                             , &
       &                                                                                           cosmologicalMassVariance_  =cosmologicalMassVariance_                       , &
       &                                                                                           darkMatterParticle_        =darkMatterParticle_                             , &
       &                                                                                           tableStore                 =.false.                                           &
       &                                                                                          )
  criticalOverdensityDescending_     =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt(                                                                              &
       &                                                                                           linearGrowth_              =linearGrowth_                                   , &
       &                                                                                           cosmologyFunctions_        =cosmologyFunctions_                             , &
       &                                                                                           cosmologicalMassVariance_  =cosmologicalMassVariance_                       , &
       &                                                                                           darkMatterParticle_        =darkMatterParticle_                             , &
       &                                                                                           tableStore                 =.false.                                           &
       &                                                                                          )
  ! Request the collapse times in opposite orders, recording each as soon as it is requested.
  do i=1,size(overdensities)
     timesAscending (                      i)=criticalOverdensityAscending_ %timeOfCollapse(criticalOverdensity=overdensities(                      i))
     timesDescending(size(overdensities)+1-i)=criticalOverdensityDescending_%timeOfCollapse(criticalOverdensity=overdensities(size(overdensities)+1-i))
  end do
  call Assert('collapse time is independent of the order in which critical overdensities are requested',timesAscending,timesDescending,absTol=0.0d0)
  ! Request the collapse times again, now that the tabulation spans all of them. They must be unchanged by the extensions of the
  ! tabulation made since they were first requested.
  do i=1,size(overdensities)
     timesRepeat(i)=criticalOverdensityAscending_%timeOfCollapse(criticalOverdensity=overdensities(i))
  end do
  call Assert('collapse time is unchanged by subsequent extension of the tabulation'                   ,timesRepeat   ,timesAscending ,absTol=0.0d0)
  ! Check that each collapse time satisfies the defining condition, δ_c(t)/D(t)=δ. The tolerance allows for the linear
  ! interpolation of the tabulated collapse times.
  satisfiesCondition=.true.
  do i=1,size(overdensities)
     satisfiesCondition=satisfiesCondition                                                                                            &
          &             .and.                                                                                                         &
          &              abs(                                                                                                         &
          &                  +criticalOverdensityAscending_%value(time=timesAscending(i))/linearGrowth_%value(time=timesAscending(i)) &
          &                  /overdensities(i)                                                                                        &
          &                  -1.0d0                                                                                                   &
          &                 )                                                                                                         &
          &             <                                                                                                             &
          &              1.0d-8
  end do
  call Assert('collapse time satisfies δ_c(t)/D(t)=δ'                                                ,satisfiesCondition,.true.)
  ! Construct mass-dependent critical overdensities, wrapping the (mass-independent) critical overdensity used above.
  darkMatterParticleWDM_             =darkMatterParticleWDMThermal     (                                                           &
       &                                                                mass                       =1.0d0                        , &
       &                                                                degreesOfFreedomEffective  =1.5d0                        , &
       &                                                                cosmologyParameters_       =cosmologyParameters_           &
       &                                                               )
  criticalOverdensityMassAscending_  =criticalOverdensityBarkana2001WDM(                                                           &
       &                                                                criticalOverdensityCDM     =criticalOverdensityAscending_, &
       &                                                                cosmologyParameters_       =cosmologyParameters_         , &
       &                                                                cosmologyFunctions_        =cosmologyFunctions_          , &
       &                                                                cosmologicalMassVariance_  =cosmologicalMassVariance_    , &
       &                                                                darkMatterParticle_        =darkMatterParticleWDM_       , &
       &                                                                linearGrowth_              =linearGrowth_                , &
       &                                                                useFittingFunction         =.true.                         &
       &                                                               )
  criticalOverdensityMassDescending_ =criticalOverdensityBarkana2001WDM(                                                           &
       &                                                                criticalOverdensityCDM     =criticalOverdensityAscending_, &
       &                                                                cosmologyParameters_       =cosmologyParameters_         , &
       &                                                                cosmologyFunctions_        =cosmologyFunctions_          , &
       &                                                                cosmologicalMassVariance_  =cosmologicalMassVariance_    , &
       &                                                                darkMatterParticle_        =darkMatterParticleWDM_       , &
       &                                                                linearGrowth_              =linearGrowth_                , &
       &                                                                useFittingFunction         =.true.                         &
       &                                                               )
  ! Request the collapse times in opposite orders, recording each as soon as it is requested.
  do i=1,size(overdensities)
     timesAscending (                      i)=criticalOverdensityMassAscending_ %timeOfCollapse(criticalOverdensity=overdensities(                      i),mass=masses(                      i))
     timesDescending(size(overdensities)+1-i)=criticalOverdensityMassDescending_%timeOfCollapse(criticalOverdensity=overdensities(size(overdensities)+1-i),mass=masses(size(overdensities)+1-i))
  end do
  call Assert('mass-dependent collapse time is independent of the order in which it is requested',timesAscending,timesDescending,absTol=0.0d0)
  ! Request the collapse times again, after all other requests have been made.
  do i=1,size(overdensities)
     timesRepeat(i)=criticalOverdensityMassAscending_%timeOfCollapse(criticalOverdensity=overdensities(i),mass=masses(i))
  end do
  call Assert('mass-dependent collapse time is unchanged by intervening requests'                ,timesRepeat   ,timesAscending ,absTol=0.0d0)
  ! Check that each collapse time satisfies the defining condition, δ_c(t,M)σ(M,t₀)/σ(M,t)=δ. No interpolation is involved, so
  ! the tolerance need only allow for that of the root finder.
  satisfiesCondition=.true.
  do i=1,size(overdensities)
     satisfiesCondition=satisfiesCondition                                                                                              &
          &             .and.                                                                                                           &
          &              abs(                                                                                                           &
          &                  +criticalOverdensityMassAscending_%value       (time=timesAscending(i)                    ,mass=masses(i)) &
          &                  *cosmologicalMassVariance_        %rootVariance(time=cosmologyFunctions_%cosmicTime(1.0d0),mass=masses(i)) &
          &                  /cosmologicalMassVariance_        %rootVariance(time=timesAscending(i)                    ,mass=masses(i)) &
          &                  /overdensities(i)                                                                                          &
          &                  -1.0d0                                                                                                     &
          &                 )                                                                                                           &
          &             <                                                                                                               &
          &              1.0d-8
  end do
  call Assert('mass-dependent collapse time satisfies δ_c(t,M)σ(M,t₀)/σ(M,t)=δ'                       ,satisfiesCondition,.true.)
  ! Construct critical overdensities in an overdense environment, using the peak-background split. The environment is
  ! independent of node and tree, so the mass-independent case uses the tabulation of collapse times.
  haloEnvironment_                   =haloEnvironmentFixed                  (                                                               &
       &                                                                     cosmologyFunctions_        =cosmologyFunctions_              , &
       &                                                                     linearGrowth_              =linearGrowth_                    , &
       &                                                                     overdensity                =overdensityEnvironment           , &
       &                                                                     massEnvironment            =1.0d15                             &
       &                                                                    )
  criticalOverdensityEnvironment_    =criticalOverdensityPeakBackgroundSplit(                                                               &
       &                                                                     criticalOverdensity_       =criticalOverdensityAscending_    , &
       &                                                                     haloEnvironment_           =haloEnvironment_                 , &
       &                                                                     cosmologyFunctions_        =cosmologyFunctions_              , &
       &                                                                     cosmologicalMassVariance_  =cosmologicalMassVariance_        , &
       &                                                                     linearGrowth_              =linearGrowth_                      &
       &                                                                    )
  criticalOverdensityMassEnvironment_=criticalOverdensityPeakBackgroundSplit(                                                               &
       &                                                                     criticalOverdensity_       =criticalOverdensityMassAscending_, &
       &                                                                     haloEnvironment_           =haloEnvironment_                 , &
       &                                                                     cosmologyFunctions_        =cosmologyFunctions_              , &
       &                                                                     cosmologicalMassVariance_  =cosmologicalMassVariance_        , &
       &                                                                     linearGrowth_              =linearGrowth_                      &
       &                                                                    )
  ! Check that the collapse times satisfy the defining conditions.
  time=criticalOverdensityEnvironment_    %timeOfCollapse(criticalOverdensity=overdensityInEnvironment)
  call Assert(                                                                                                                     &
       &      'collapse time in an overdense environment satisfies δ_c(t)/D(t)=δ'                                                , &
       &      +criticalOverdensityEnvironment_%value(time=time)/linearGrowth_%value(time=time)/overdensityInEnvironment          , &
       &      1.0d0                                                                                                              , &
       &      relTol=1.0d-8                                                                                                        &
       &     )
  time=criticalOverdensityMassEnvironment_%timeOfCollapse(criticalOverdensity=overdensityInEnvironment,mass=massInEnvironment)
  call Assert(                                                                                                                     &
       &      'mass-dependent collapse time in an overdense environment satisfies δ_c(t,M)σ(M,t₀)/σ(M,t)=δ'                      , &
       &      +criticalOverdensityMassEnvironment_%value       (time=time                                 ,mass=massInEnvironment) &
       &      *cosmologicalMassVariance_          %rootVariance(time=cosmologyFunctions_%cosmicTime(1.0d0),mass=massInEnvironment) &
       &      /cosmologicalMassVariance_          %rootVariance(time=time                                 ,mass=massInEnvironment) &
       &      /overdensityInEnvironment                                                                                          , &
       &      1.0d0                                                                                                              , &
       &      relTol=1.0d-8                                                                                                        &
       &     )
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_Critical_Overdensity_Time_Of_Collapse_Determinism
