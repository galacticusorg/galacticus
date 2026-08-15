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

program Tests_Merger_Tree_Branching
  !!{RST
  Tests of merger tree branching rates.
  !!}
  use :: Cosmological_Density_Field          , only : cosmologicalMassVarianceFilteredPower          , criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt
  use :: Cosmology_Functions                 , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters                , only : cosmologyParametersSimple
  use :: Dark_Matter_Particles               , only : darkMatterParticleCDM
  use :: Display                             , only : displayVerbositySet                            , verbosityLevelWorking
  use :: Events_Hooks                        , only : eventsHooksInitialize
  use :: Excursion_Sets_Barriers             , only : excursionSetBarrierCriticalOverdensity
  use :: Excursion_Sets_First_Crossings      , only : excursionSetFirstCrossingFarahiMidpoint        , excursionSetFirstCrossingFarahiMidpointBrownianBridge       , excursionSetFirstCrossingLinearBarrier, excursionSetFirstCrossingLinearBarrierBrownianBridge
  use :: Galacticus_Nodes                    , only : treeNode
  use :: ISO_Varying_String                  , only : var_str
  use :: Linear_Growth                       , only : linearGrowthCollisionlessMatter
  use, intrinsic :: ISO_C_Binding             , only : c_long
  use :: Error                               , only : Error_Report                                   , Error_Handler_Register
  use :: Merger_Tree_Branching               , only : mergerTreeBranchingProbabilityGnrlzdPrssSchchtr, mergerTreeBranchingProbabilityParkinsonColeHelly            , mergerTreeBranchingProbabilityPCHPlus, mergerTreeBranchingProbabilityClass                  , &
       &                                              mergerTreeBranchingBoundLower                  , mergerTreeBranchingBoundUpper
  use :: Merger_Tree_Branching_Modifiers     , only : mergerTreeBranchingProbabilityModifierIdentity , mergerTreeBranchingProbabilityModifierPCHPlus
  use :: Numerical_Constants_Math            , only : Pi
  use :: Numerical_Random_Numbers            , only : randomNumberGeneratorGSL
  use :: Power_Spectra_Primordial            , only : powerSpectrumPrimordialPowerLaw
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredSimple
  use :: Power_Spectrum_Window_Functions     , only : powerSpectrumWindowFunctionSharpKSpace
  use :: Transfer_Functions                  , only : transferFunctionIdentity
  use :: Unit_Tests                          , only : Assert                                         , Unit_Tests_Begin_Group                                      , Unit_Tests_End_Group                 , Unit_Tests_Finish                                    , &
       &                                              compareLessThan                                , compareLessThanOrEqual                                      , compareGreaterThan
  implicit none
  type            (cosmologyParametersSimple                                   )                           :: cosmologyParametersSimple_
  type            (cosmologyFunctionsMatterLambda                              )                           :: cosmologyFunctionsMatterLambda_
  type            (cosmologicalMassVarianceFilteredPower                       )                           :: cosmologicalMassVarianceFilteredPower_
  type            (powerSpectrumWindowFunctionSharpKSpace                      )                           :: powerSpectrumWindowFunctionSharpKSpace_
  type            (powerSpectrumPrimordialPowerLaw                             )                           :: powerSpectrumPrimordialPowerLaw_
  type            (transferFunctionIdentity                                    )                           :: transferFunctionIdentity_
  type            (powerSpectrumPrimordialTransferredSimple                    )                           :: powerSpectrumPrimordialTransferredSimple_
  type            (linearGrowthCollisionlessMatter                             )                           :: linearGrowthCollisionlessMatter_
  type            (excursionSetFirstCrossingLinearBarrier                      )                           :: excursionSetFirstCrossingLinearBarrier_
  type            (excursionSetFirstCrossingLinearBarrierBrownianBridge        )                           :: excursionSetFirstCrossingLinearBarrierBrownianBridge_
  type            (excursionSetFirstCrossingFarahiMidpoint                     )                           :: excursionSetFirstCrossingFarahiMidpoint_
  type            (excursionSetFirstCrossingFarahiMidpointBrownianBridge       )                           :: excursionSetFirstCrossingFarahiMidpointBrownianBridge_
  type            (excursionSetBarrierCriticalOverdensity                      )                           :: excursionSetBarrierCriticalOverdensity_
  type            (darkMatterParticleCDM                                       )                           :: darkMatterParticleCDM_
  type            (treeNode                                                    )                           :: node
  type            (mergerTreeBranchingProbabilityModifierIdentity              )                           :: mergerTreeBranchingProbabilityModifierIdentity_
  type            (criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt)                           :: criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_
  type            (mergerTreeBranchingProbabilityParkinsonColeHelly            ), dimension(            3) :: mergerTreeBranchingProbabilityParkinsonColeHelly_
  type            (mergerTreeBranchingProbabilityGnrlzdPrssSchchtr             ), dimension(            5) :: mergerTreeBranchingProbabilityGnrlzdPrssSchchtr_
  double precision                                                                                         :: time                                                                              , rootVarianceParent                             , &
       &                                                                                                      rootVarianceResolution                                                            , branchingProbabilityRate                       , &
       &                                                                                                      accretionRate                                                                     , criticalOverdensity_                           , &
       &                                                                                                      expansionFactor                                                                   , timeNow                                        , &
       &                                                                                                      branchingProbabilityRateTargetGeneral                                             , accretionRateTargetGeneral                     , &
       &                                                                                                      smoothAccretionRateTargetGeneral                                                  , smoothAccretionRate                            , &
       &                                                                                                      varianceParent                                                                    , varianceResolution                             , &
       &                                                                                                      errorMaximumUnconstrained                                                         , errorMaximumConstrained
  double precision                                                              , parameter                :: branchingProbabilityRateTarget                                =2.498324530964044d0, accretionRateTarget       =4.181139013841312d-1
  double precision                                                              , dimension(            2) :: redshift                                                      =[0.0d0,1.0d0]
  double precision                                                              , parameter                :: massParent                                                    =1.0d12             , massResolution            =1.0d11              , &
       &                                                                                                      fractionalTimeStep                                            =1.0d-3             , varianceConstrained       =2.0d02              , &
       &                                                                                                      criticalOverdensityConstrained                                =4.0d00             , toleranceVariance         =1.0d-1
  integer                                                                       , parameter                :: countVariance                                                 =1000
  double precision                                                              , dimension(countVariance) :: branchingRateUnconstrainedAnalytic                                                , branchingRateUnconstrainedNumerical            , &
       &                                                                                                      branchingRateConstrainedAnalytic                                                  , branchingRateConstrainedNumerical              , &
       &                                                                                                      variance_
  integer                                                                                                  :: i                                                                                 , j
  character       (len=16                                                      )                           :: label
  ! Objects and data used by the analytic-identity tests which follow the tests above. For an n=-1 power spectrum filtered with a
  ! sharp k-space window function σ(M) ∝ M^{-1/3}, so α ≡ dlnσ/dlnM is a constant. Every identity tested in that section rests on
  ! that constancy - see the comments at the head of each group.
  integer                                                                       , parameter                :: countIdentity             =6      , countRatio     = 6, &
       &                                                                                                      countClosedForm           =6      , countX         =16, &
       &                                                                                                      countAlpha                =5      , countGrid      =50, &
       &                                                                                                      countSample               =20000   , countInverse   = 3
  type            (mergerTreeBranchingProbabilityParkinsonColeHelly            ), dimension(countIdentity  ) :: branchingIdentity_
  type            (mergerTreeBranchingProbabilityParkinsonColeHelly            ), dimension(countClosedForm) :: branchingClosedForm_
  type            (mergerTreeBranchingProbabilityParkinsonColeHelly            )                           :: branchingOrdering_                , branchingPole_    , &
       &                                                                                                      branchingGeneric_                 , branchingSampler_ , &
       &                                                                                                      branchingFresh_                   , branchingReference_, &
       &                                                                                                      branchingSamplerGeneric_
  type            (mergerTreeBranchingProbabilityPCHPlus                       ), dimension(            3) :: branchingPCHPlus_
  type            (randomNumberGeneratorGSL                                    )                           :: randomNumberGenerator_
  ! Parameters of the configurations swept over. The values of γ₁ are chosen so that the effective exponent γ₁-1/α (=γ₁+3 here)
  ! avoids the odd integers at which the hypergeometric functions used to bound the branching probability are singular; γ₁=0,
  ! which does sit on such a point, is tested separately through `branchingPole_`.
  double precision                                                              , dimension(countIdentity  ), parameter :: G0Identity          =[ 0.57d0, 0.57d0, 1.00d0, 1.00d0, 0.80d0, 0.80d0]                    , &
       &                                                                                                                  gamma1Identity      =[ 0.38d0, 0.38d0,-0.50d0,-0.50d0, 0.80d0, 0.80d0]                    , &
       &                                                                                                                  gamma2Identity      =[-0.01d0,-0.01d0, 0.00d0, 0.00d0,-0.30d0,-0.30d0]
  logical                                                                       , dimension(countIdentity  ), parameter :: tabulateIdentity    =[.true. ,.false.,.true. ,.false.,.true. ,.false.]
  double precision                                                              , dimension(countClosedForm), parameter :: gamma1ClosedForm    =[ 0.00d0, 0.00d0,-1.00d0,-1.00d0,-3.00d0,-3.00d0]
  logical                                                                       , dimension(countClosedForm), parameter :: tabulateClosedForm  =[.true. ,.false.,.true. ,.false.,.true. ,.false.]
  double precision                                                              , dimension(countRatio     ), parameter :: ratioResolution     =[2.5d0,1.0d1,1.0d2,1.0d3,1.0d6,1.0d8]
  ! The three PCH+ configurations tested: a null one, one exercising only γ₄ and γ₅, and one exercising only γ₃.
  double precision                                                              , dimension(            3  ), parameter :: gamma3PCHPlus       =[0.0d0,0.0d0,0.3d0]                    , &
       &                                                                                                                  gamma4PCHPlus       =[0.0d0,0.5d0,0.0d0]                    , &
       &                                                                                                                  gamma5PCHPlus       =[0.0d0,0.3d0,0.0d0]
  double precision                                                              , dimension(countAlpha     ), parameter :: massAlpha           =[1.0d8,1.0d10,1.0d12,1.0d14,1.0d16]
  double precision                                                              , dimension(countInverse   ), parameter :: fractionInverse     =[0.25d0,0.50d0,0.75d0]
  double precision                                                              , dimension(countX         ), parameter :: xTarget             =[1.00d-7,1.00d-6,1.00d-5,1.00d-4,1.00d-3,1.00d-2,5.00d-2,9.90d-2,1.00d-1,1.01d-1,5.00d-1,1.00d0,5.00d0,1.25d1,1.00d2,1.00d3]
  ! Tolerances. These are set from the accuracy actually achieved - see the comment at each assertion for what limits it.
  double precision                                                              , parameter                :: toleranceIdentityTabulated=1.0d-3 , toleranceIdentityDirect   =3.0d-3, &
       &                                                                                                      toleranceRateIntegral     =3.0d-3 , toleranceMassBranchInverse=1.0d-3, &
       &                                                                                                      toleranceReduction        =1.0d-9
  ! Critical value of the Kolmogorov-Smirnov statistic. For a sample of this size the probability that a correctly drawn sample
  ! exceeds this is 2exp(-2Nd²) ≈ 2·10⁻⁷, so this cannot flap - and, the random sequence being seeded deterministically, the
  ! statistic is in any case identical from run to run of the same executable. The values actually attained are 4.9·10⁻³ and
  ! 6.3·10⁻³ for the two samplers, against the 8.6·10⁻³/√N ≈ 6.1·10⁻³ expected of a correctly drawn sample of this size, so both
  ! agree with the analytic distribution to within the sampling noise and this limit is left with a margin of some three times.
  double precision                                                              , parameter                :: toleranceKolmogorovSmirnov=2.0d-2
  ! Tolerances for the closed-form comparisons, one per configuration. The tabulated evaluation interpolates linearly between
  ! twenty points per decade, so its error grows as the square of the logarithmic gradient of the tabulated function: that
  ! function falls as x^{-1} for γ₁=0 but as x^{-4} for γ₁=-3, which is what makes the last of these so much the loosest. Values
  ! of γ₁ used in practice are of order a few tenths, where the tabulation behaves like the first of these.
  double precision                                                              , dimension(countClosedForm), parameter :: toleranceClosedForm =[1.0d-3,1.0d-4,2.0d-3,1.0d-4,2.0d-2,1.0d-4]
  double precision                                                                                         :: timeAnalytic                      , deltaAnalytic     , &
       &                                                                                                      sigmaParentAnalytic               , sigmaHalfAnalytic , &
       &                                                                                                      sigmaResolutionAnalytic           , massResolution_   , &
       &                                                                                                      probabilityFull                   , &
       &                                                                                                      massBranchTest                    , massTargetCurrent , &
       &                                                                                                      probabilityPartial                , stepExpected      , &
       &                                                                                                      modifierExpected
  double precision                                                              , dimension(countRatio     ) :: probabilityIdentity             , boundLowerIdentity, &
       &                                                                                                      boundUpperIdentity
  double precision                                                              , dimension(countX         ) :: fractionActual                  , fractionExpected  , &
       &                                                                                                      xActual
  double precision                                                              , dimension(countAlpha     ) :: alphaMeasured                   , sigmaScratch
  double precision                                                              , dimension(countInverse   ) :: massBranchReturned              , massBranchTarget
  double precision                                                              , dimension(countGrid      ) :: massGrid                        , cdfAnalytic       , &
       &                                                                                                      cdfSampled
  double precision                                                              , dimension(countSample    ) :: massSample
  integer                                                                                                  :: k

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelWorking)
  ! Register error handlers. This installs Galacticus' own GSL error handler, without which GSL aborts through its default handler
  ! irrespective of any `status` argument requesting that the failure be returned instead. Several branching probability code paths
  ! rely on being able to detect and recover from a failed evaluation of a hypergeometric function, and those paths cannot be
  ! reached at all unless this handler is in place.
  call Error_Handler_Register()
  ! Initialize event hooks.
  call eventsHooksInitialize()
  ! Build all objects needed for these tests.
  darkMatterParticleCDM_                                        =darkMatterParticleCDM                                        (                                                                                                        &
       &                                                                                                                      )
  cosmologyParametersSimple_                                    =cosmologyParametersSimple                                    (                                                                                                        &
       &                                                                                                                       OmegaMatter                             = 1.00d0                                                      , &
       &                                                                                                                       OmegaBaryon                             = 0.00d0                                                      , &
       &                                                                                                                       OmegaDarkEnergy                         = 0.00d0                                                      , &
       &                                                                                                                       temperatureCMB                          = 2.78d0                                                      , &
       &                                                                                                                       HubbleConstant                          =70.00d0                                                        &
       &                                                                                                                      )
  cosmologyFunctionsMatterLambda_                               =cosmologyFunctionsMatterLambda                               (                                                                                                        &
       &                                                                                                                       cosmologyParameters_                    =cosmologyParametersSimple_                                     &
       &                                                                                                                      )
  linearGrowthCollisionlessMatter_                              =linearGrowthCollisionlessMatter                              (                                                                                                        &
       &                                                                                                                       cosmologyParameters_                    =cosmologyParametersSimple_                                   , &
       &                                                                                                                       cosmologyFunctions_                     =cosmologyFunctionsMatterLambda_                                &
       &                                                                                                                      )
  powerSpectrumPrimordialPowerLaw_                              =powerSpectrumPrimordialPowerLaw                              (                                                                                                        &
       &                                                                                                                       index_                                  =-1.0d0                                                       , &
       &                                                                                                                       running                                 =+0.0d0                                                       , &
       &                                                                                                                       runningRunning                          =+0.0d0                                                       , &
       &                                                                                                                       wavenumberReference                     =+1.0d0                                                       , &
       &                                                                                                                       runningSmallScalesOnly                  =.false.                                                        &
       &                                                                                                                      )
  transferFunctionIdentity_                                     =transferFunctionIdentity                                     (                                                                                                        &
       &                                                                                                                       cosmologyParameters_                    =cosmologyParametersSimple_                                   , &
       &                                                                                                                       time                                    =13.8d0                                                         &
       &                                                                                                                      )
  powerSpectrumPrimordialTransferredSimple_                     =powerSpectrumPrimordialTransferredSimple                     (                                                                                                        &
       &                                                                                                                       powerSpectrumPrimordial_                =powerSpectrumPrimordialPowerLaw_                             , &
       &                                                                                                                       transferFunction_                       =transferFunctionIdentity_                                    , &
       &                                                                                                                       linearGrowth_                           =linearGrowthCollisionlessMatter_                               &
       &                                                                                                                      )
  powerSpectrumWindowFunctionSharpKSpace_                       =powerSpectrumWindowFunctionSharpKSpace                       (                                                                                                        &
       &                                                                                                                       cosmologyParameters_                    =cosmologyParametersSimple_                                   , &
       &                                                                                                                       normalization                           =0.0d0                                                          &
       &                                                                                                                      )

  cosmologicalMassVarianceFilteredPower_                        =cosmologicalMassVarianceFilteredPower                        (                                                                                                        &
       &                                                                                                                       sigma8                                  =1.0d+0                                                       , &
       &                                                                                                                       tolerance                               =1.0d-4                                                       , &
       &                                                                                                                       toleranceTopHat                         =1.0d-4                                                       , &
       &                                                                                                                       rootVarianceLogarithmicGradientTolerance=1.0d-9                                                       , &
       &                                                                                                                       integrationFailureIsFatal               =.true.                                                       , &
       &                                                                                                                       storeTabulations                        =.true.                                                       , &
       &                                                                                                                       nonMonotonicIsFatal                     =.true.                                                       , &
       &                                                                                                                       monotonicInterpolation                  =.false.                                                      , &
       &                                                                                                                       truncateAtParticleHorizon               =.false.                                                      , &
       &                                                                                                                       cosmologyParameters_                    =cosmologyParametersSimple_                                   , &
       &                                                                                                                       cosmologyFunctions_                     =cosmologyFunctionsMatterLambda_                              , &
       &                                                                                                                       linearGrowth_                           =linearGrowthCollisionlessMatter_                             , &
       &                                                                                                                       powerSpectrumPrimordialTransferred_     =powerSpectrumPrimordialTransferredSimple_                    , &
       &                                                                                                                       powerSpectrumWindowFunction_            =powerSpectrumWindowFunctionSharpKSpace_                        &
       &                                                                                                                      )
  criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_=criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt  (                                                                                                        &
       &                                                                                                                       linearGrowth_                           =linearGrowthCollisionlessMatter_                             , &
       &                                                                                                                       cosmologyFunctions_                     =cosmologyFunctionsMatterLambda_                              , &
       &                                                                                                                       cosmologicalMassVariance_               =cosmologicalMassVarianceFilteredPower_                       , &
       &                                                                                                                       darkMatterParticle_                     =darkMatterParticleCDM_                                       , &
       &                                                                                                                       tableStore                              =.true.                                                         &
       &                                                                                                                      )
  mergerTreeBranchingProbabilityModifierIdentity_               =mergerTreeBranchingProbabilityModifierIdentity               (                                                                                                        &
       &                                                                                                                      )
  mergerTreeBranchingProbabilityParkinsonColeHelly_(1)          =mergerTreeBranchingProbabilityParkinsonColeHelly             (                                                                                                        &
       &                                                                                                                       G0                                      =1.0d+0                                                       , &
       &                                                                                                                       gamma1                                  =0.0d+0                                                       , &
       &                                                                                                                       gamma2                                  =0.0d+0                                                       , &
       &                                                                                                                       accuracyFirstOrder                      =1.0d-1                                                       , &
       &                                                                                                                       precisionHypergeometric                 =1.0d-6                                                       , &
       &                                                                                                                       hypergeometricTabulate                  =.true.                                                       , &
       &                                                                                                                       cdmAssumptions                          =.true.                                                       , &
       &                                                                                                                       tolerateRoundOffErrors                  =.false.                                                      , &
       &                                                                                                                       cosmologicalMassVariance_               =cosmologicalMassVarianceFilteredPower_                       , &
       &                                                                                                                       criticalOverdensity_                    =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_  &
       &                                                                                                                      )
  mergerTreeBranchingProbabilityParkinsonColeHelly_(2)          =mergerTreeBranchingProbabilityParkinsonColeHelly             (                                                                                                        &
       &                                                                                                                       G0                                      =1.0d+0                                                       , &
       &                                                                                                                       gamma1                                  =0.0d+0                                                       , &
       &                                                                                                                       gamma2                                  =0.0d+0                                                       , &
       &                                                                                                                       accuracyFirstOrder                      =1.0d-1                                                       , &
       &                                                                                                                       precisionHypergeometric                 =1.0d-6                                                       , &
       &                                                                                                                       hypergeometricTabulate                  =.false.                                                      , &
       &                                                                                                                       cdmAssumptions                          =.true.                                                       , &
       &                                                                                                                       tolerateRoundOffErrors                  =.false.                                                      , &
       &                                                                                                                       cosmologicalMassVariance_               =cosmologicalMassVarianceFilteredPower_                       , &
       &                                                                                                                       criticalOverdensity_                    =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_  &
       &                                                                                                                      )
  mergerTreeBranchingProbabilityParkinsonColeHelly_(3)          =mergerTreeBranchingProbabilityParkinsonColeHelly             (                                                                                                        &
       &                                                                                                                       G0                                      =1.0d+0                                                       , &
       &                                                                                                                       gamma1                                  =0.0d+0                                                       , &
       &                                                                                                                       gamma2                                  =0.0d+0                                                       , &
       &                                                                                                                       accuracyFirstOrder                      =1.0d-1                                                       , &
       &                                                                                                                       precisionHypergeometric                 =1.0d-6                                                       , &
       &                                                                                                                       hypergeometricTabulate                  =.false.                                                      , &
       &                                                                                                                       cdmAssumptions                          =.false.                                                      , &
       &                                                                                                                       tolerateRoundOffErrors                  =.false.                                                      , &
       &                                                                                                                       cosmologicalMassVariance_               =cosmologicalMassVarianceFilteredPower_                       , &
       &                                                                                                                       criticalOverdensity_                    =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_  &
       &                                                                                                                      )
  excursionSetBarrierCriticalOverdensity_                       =excursionSetBarrierCriticalOverdensity                       (                                                                                                        &
       &                                                                                                                       criticalOverdensity_                    =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_, &
       &                                                                                                                       cosmologicalMassVariance_               =cosmologicalMassVarianceFilteredPower_                         &
       &                                                                                                                      )
  excursionSetFirstCrossingLinearBarrier_                       =excursionSetFirstCrossingLinearBarrier                       (                                                                                                        &
       &                                                                                                                       fractionalTimeStep                      =fractionalTimeStep                                           , &
       &                                                                                                                       excursionSetBarrier_                    =excursionSetBarrierCriticalOverdensity_                      , &
       &                                                                                                                       cosmologicalMassVariance_               =cosmologicalMassVarianceFilteredPower_                         &
       &                                                                                                                      )
  excursionSetFirstCrossingLinearBarrierBrownianBridge_         =excursionSetFirstCrossingLinearBarrierBrownianBridge         (                                                                                                        &
       &                                                                                                                       varianceConstrained                     =varianceConstrained                                          , &
       &                                                                                                                       criticalOverdensityConstrained          =criticalOverdensityConstrained                               , &
       &                                                                                                                       fractionalTimeStep                      =fractionalTimeStep                                           , &
       &                                                                                                                       cosmologyFunctions_                     =cosmologyFunctionsMatterLambda_                              , &
       &                                                                                                                       excursionSetFirstCrossing_              =excursionSetFirstCrossingLinearBarrier_                      , &
       &                                                                                                                       excursionSetBarrier_                    =excursionSetBarrierCriticalOverdensity_                      , &
       &                                                                                                                       cosmologicalMassVariance_               =cosmologicalMassVarianceFilteredPower_                       , &
       &                                                                                                                       criticalOverdensity_                    =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_, &
       &                                                                                                                       linearGrowth_                           =linearGrowthCollisionlessMatter_                               &
       &                                                                                                                      )
  excursionSetFirstCrossingFarahiMidpoint_                      =excursionSetFirstCrossingFarahiMidpoint                      (                                                                                                        &
       &                                                                                                                       fractionalTimeStep                      =fractionalTimeStep                                           , &
       &                                                                                                                       fileName                                =var_str("auto")                                              , &
       &                                                                                                                       varianceNumberPerUnitProbability        =100                                                          , &
       &                                                                                                                       varianceNumberPerUnit                   = 16                                                          , &
       &                                                                                                                       varianceNumberPerDecade                 = 32                                                          , &
       &                                                                                                                       varianceNumberPerDecadeNonCrossing      =  8                                                          , &
       &                                                                                                                       timeNumberPerDecade                     = 16                                                          , &
       &                                                                                                                       varianceIsUnlimited                     =.true.                                                       , &
       &                                                                                                                       cosmologyFunctions_                     =cosmologyFunctionsMatterLambda_                              , &
       &                                                                                                                       excursionSetBarrier_                    =excursionSetBarrierCriticalOverdensity_                      , &
       &                                                                                                                       cosmologicalMassVariance_               =cosmologicalMassVarianceFilteredPower_                         &
       &                                                                                                                      )
  excursionSetFirstCrossingFarahiMidpointBrownianBridge_        =excursionSetFirstCrossingFarahiMidpointBrownianBridge        (                                                                                                        &
       &                                                                                                                       varianceConstrained                     =varianceConstrained                                          , &
       &                                                                                                                       criticalOverdensityConstrained          =  criticalOverdensityConstrained                             , &
       &                                                                                                                       fractionalTimeStep                      =fractionalTimeStep                                           , &
       &                                                                                                                       fileName                                =var_str("auto")                                              , &
       &                                                                                                                       varianceNumberPerUnitProbability        =100                                                          , &
       &                                                                                                                       varianceNumberPerUnit                   = 32                                                          , &
       &                                                                                                                       varianceNumberPerDecade                 = 32                                                          , &
       &                                                                                                                       varianceNumberPerDecadeNonCrossing      =  8                                                          , &
       &                                                                                                                       timeNumberPerDecade                     = 64                                                          , &
       &                                                                                                                       varianceIsUnlimited                     =.false.                                                      , &
       &                                                                                                                       cosmologyFunctions_                     =cosmologyFunctionsMatterLambda_                              , &
       &                                                                                                                       excursionSetBarrier_                    =excursionSetBarrierCriticalOverdensity_                      , &
       &                                                                                                                       cosmologicalMassVariance_               =cosmologicalMassVarianceFilteredPower_                       , &
       &                                                                                                                       criticalOverdensity_                    =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_, &
       &                                                                                                                       linearGrowth_                           =linearGrowthCollisionlessMatter_                             , &
       &                                                                                                                       excursionSetFirstCrossing_              =excursionSetFirstCrossingLinearBarrier_                        &
       &                                                                                                                      )
  mergerTreeBranchingProbabilityGnrlzdPrssSchchtr_(1)           =mergerTreeBranchingProbabilityGnrlzdPrssSchchtr              (                                                                                                        &
       &                                                                                                                       deltaStepMaximum                        =1.0d-1                                                       , &
       &                                                                                                                       massMinimum                             =1.0d+0                                                       , &
       &                                                                                                                       smoothAccretion                         =.false.                                                      , &
       &                                                                                                                       distributionFunctionLowerHalfOnly       =.true.                                                       , &
       &                                                                                                                       distributionFunctionNormalize           =.true.                                                       , &
       &                                                                                                                       cosmologyFunctions_                     =cosmologyFunctionsMatterLambda_                              , &
       &                                                                                                                       criticalOverdensity_                    =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_, &
       &                                                                                                                       cosmologicalMassVariance_               =cosmologicalMassVarianceFilteredPower_                       , &
       &                                                                                                                       excursionSetFirstCrossing_              =excursionSetFirstCrossingLinearBarrier_                      , &
       &                                                                                                                       mergerTreeBranchingProbabilityModifier_ =mergerTreeBranchingProbabilityModifierIdentity_                &
       &                                                                                                                      )
  mergerTreeBranchingProbabilityGnrlzdPrssSchchtr_(2)           =mergerTreeBranchingProbabilityGnrlzdPrssSchchtr              (                                                                                                        &
       &                                                                                                                       deltaStepMaximum                        =1.0d-1                                                       , &
       &                                                                                                                       massMinimum                             =1.0d-1*massResolution                                        , &
       &                                                                                                                       smoothAccretion                         =.false.                                                      , &
       &                                                                                                                       distributionFunctionLowerHalfOnly       =.true.                                                       , &
       &                                                                                                                       distributionFunctionNormalize           =.true.                                                       , &
       &                                                                                                                       cosmologyFunctions_                     =cosmologyFunctionsMatterLambda_                              , &
       &                                                                                                                       criticalOverdensity_                    =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_, &
       &                                                                                                                       cosmologicalMassVariance_               =cosmologicalMassVarianceFilteredPower_                       , &
       &                                                                                                                       excursionSetFirstCrossing_              =excursionSetFirstCrossingLinearBarrier_                      , &
       &                                                                                                                       mergerTreeBranchingProbabilityModifier_ =mergerTreeBranchingProbabilityModifierIdentity_                &
       &                                                                                                                      )
  mergerTreeBranchingProbabilityGnrlzdPrssSchchtr_(3)           =mergerTreeBranchingProbabilityGnrlzdPrssSchchtr              (                                                                                                        &
       &                                                                                                                       deltaStepMaximum                        =1.0d-1                                                       , &
       &                                                                                                                       massMinimum                             =1.0d-1*massResolution                                        , &
       &                                                                                                                       smoothAccretion                         =.false.                                                      , &
       &                                                                                                                       distributionFunctionLowerHalfOnly       =.true.                                                       , &
       &                                                                                                                       distributionFunctionNormalize           =.true.                                                       , &
       &                                                                                                                       cosmologyFunctions_                     =cosmologyFunctionsMatterLambda_                              , &
       &                                                                                                                       criticalOverdensity_                    =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_, &
       &                                                                                                                       cosmologicalMassVariance_               =cosmologicalMassVarianceFilteredPower_                       , &
       &                                                                                                                       excursionSetFirstCrossing_              =excursionSetFirstCrossingFarahiMidpoint_                     , &
       &                                                                                                                       mergerTreeBranchingProbabilityModifier_ =mergerTreeBranchingProbabilityModifierIdentity_                &
       &                                                                                                                      )
  mergerTreeBranchingProbabilityGnrlzdPrssSchchtr_(4)           =mergerTreeBranchingProbabilityGnrlzdPrssSchchtr              (                                                                                                        &
       &                                                                                                                       deltaStepMaximum                        =1.0d-1                                                       , &
       &                                                                                                                       massMinimum                             =1.0d-1*massResolution                                        , &
       &                                                                                                                       smoothAccretion                         =.true.                                                       , &
       &                                                                                                                       distributionFunctionLowerHalfOnly       =.true.                                                       , &
       &                                                                                                                       distributionFunctionNormalize           =.true.                                                       , &
       &                                                                                                                       cosmologyFunctions_                     =cosmologyFunctionsMatterLambda_                              , &
       &                                                                                                                       criticalOverdensity_                    =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_, &
       &                                                                                                                       cosmologicalMassVariance_               =cosmologicalMassVarianceFilteredPower_                       , &
       &                                                                                                                       excursionSetFirstCrossing_              =excursionSetFirstCrossingLinearBarrier_                      , &
       &                                                                                                                       mergerTreeBranchingProbabilityModifier_ =mergerTreeBranchingProbabilityModifierIdentity_                &
       &                                                                                                                      )
  mergerTreeBranchingProbabilityGnrlzdPrssSchchtr_(5)           =mergerTreeBranchingProbabilityGnrlzdPrssSchchtr              (                                                                                                        &
       &                                                                                                                       deltaStepMaximum                        =1.0d-1                                                       , &
       &                                                                                                                       massMinimum                             =1.0d-1*massResolution                                        , &
       &                                                                                                                       smoothAccretion                         =.true.                                                       , &
       &                                                                                                                       distributionFunctionLowerHalfOnly       =.true.                                                       , &
       &                                                                                                                       distributionFunctionNormalize           =.true.                                                       , &
       &                                                                                                                       cosmologyFunctions_                     =cosmologyFunctionsMatterLambda_                              , &
       &                                                                                                                       criticalOverdensity_                    =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_, &
       &                                                                                                                       cosmologicalMassVariance_               =cosmologicalMassVarianceFilteredPower_                       , &
       &                                                                                                                       excursionSetFirstCrossing_              =excursionSetFirstCrossingFarahiMidpoint_                     , &
       &                                                                                                                       mergerTreeBranchingProbabilityModifier_ =mergerTreeBranchingProbabilityModifierIdentity_                &
       &                                                                                                                      )
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Merger tree branching")
  ! Set up physical system to be tested.
  timeNow               =cosmologyFunctionsMatterLambda_%cosmicTime(1.0d0)
  ! Iterate over redshifts.
  do i=1,size(redshift)
     expansionFactor       =+cosmologyFunctionsMatterLambda_                              %expansionFactorFromRedshift(               redshift       (i))
     time                  =+cosmologyFunctionsMatterLambda_                              %cosmicTime                 (               expansionFactor   )
     criticalOverdensity_  =+criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_%value                      (               time              ) &
          &                 *cosmologicalMassVarianceFilteredPower_                       %rootVariance               (massResolution,timeNow           ) &
          &                 /cosmologicalMassVarianceFilteredPower_                       %rootVariance               (massResolution,time              )
     rootVarianceParent    =+cosmologicalMassVarianceFilteredPower_                       %rootVariance               (massParent    ,time              )
     rootVarianceResolution=+cosmologicalMassVarianceFilteredPower_                       %rootVariance               (massResolution,time              )
     write (label,'(a,f3.1)') "Redshift z = ",redshift(i)
     call Unit_Tests_Begin_Group(label)
     ! For an n=-1 power spectrum, σ(M) ∝ M^{-⅓}.
     call Assert('σ(M) mass scaling',rootVarianceParent/rootVarianceResolution,(massParent/massResolution)**(-1.0d0/3.0d0),relTol=3.0d-5)
     ! Test branching rates and accretion rates.
     !  * For an n=-1 power spectrum the branching probability rate can be found to be:
     !      (1/σ₂) [√(2/π) ₂F₁[-1/2,2,1/2,1-x^⅔]] / sqrt{1-x₁₂^⅔};
     !  * For an n=-1 power spectrum the accretion rate can be found to be:
     !      (1/σ₂)  √(2/π) x₁₂^⅓ / sqrt{1-x₁₂^⅔};
     ! where x₁₂=M₁/M₂, M₁ is the mass resolution, and M₂ is the parent mass. Numerical result was evaluated using Mathematica.
     call Unit_Tests_Begin_Group("Parkinson-Cole-Helly branching rates"       )
     call Unit_Tests_Begin_Group("Tabulated ₂F₁; CDM assumptions"             )
     branchingProbabilityRate=mergerTreeBranchingProbabilityParkinsonColeHelly_(1)%probability          (massParent,criticalOverdensity_,time,massResolution,node)
     accretionRate           =mergerTreeBranchingProbabilityParkinsonColeHelly_(1)%fractionSubresolution(massParent,criticalOverdensity_,time,massResolution,node)
     call Assert('Branching probability rate',branchingProbabilityRate,branchingProbabilityRateTarget/rootVarianceParent,relTol=1.0d-4)
     call Assert('Accretion rate'            ,accretionRate           ,accretionRateTarget           /rootVarianceParent,relTol=2.0d-3)
     call Unit_Tests_End_Group  (                                             )
     call Unit_Tests_Begin_Group("Computed ₂F₁; CDM assumptions"              )
     branchingProbabilityRate=mergerTreeBranchingProbabilityParkinsonColeHelly_(2)%probability          (massParent,criticalOverdensity_,time,massResolution,node)
     accretionRate           =mergerTreeBranchingProbabilityParkinsonColeHelly_(2)%fractionSubresolution(massParent,criticalOverdensity_,time,massResolution,node)
     call Assert('Branching probability rate',branchingProbabilityRate,branchingProbabilityRateTarget/rootVarianceParent,relTol=1.0d-4)
     call Assert('Accretion rate'            ,accretionRate           ,accretionRateTarget           /rootVarianceParent,relTol=1.0d-4)
     call Unit_Tests_End_Group  (                                             )
     call Unit_Tests_Begin_Group("Computed ₂F₁; no CDM assumptions"           )
     branchingProbabilityRate=mergerTreeBranchingProbabilityParkinsonColeHelly_(3)%probability          (massParent,criticalOverdensity_,time,massResolution,node)
     accretionRate           =mergerTreeBranchingProbabilityParkinsonColeHelly_(3)%fractionSubresolution(massParent,criticalOverdensity_,time,massResolution,node)
     call Assert('Branching probability rate',branchingProbabilityRate,branchingProbabilityRateTarget/rootVarianceParent,relTol=1.0d-4)
     call Assert('Accretion rate'            ,accretionRate           ,accretionRateTarget           /rootVarianceParent,relTol=1.0d-4)
     call Unit_Tests_End_Group  (                                             )
     call Unit_Tests_End_Group  (                                             )
     call Unit_Tests_Begin_Group("Generalized Press-Schechter linear barrier branching rates")
     branchingProbabilityRate=mergerTreeBranchingProbabilityGnrlzdPrssSchchtr_ (1)%probability          (massParent,criticalOverdensity_,time,massResolution,node)
     accretionRate           =mergerTreeBranchingProbabilityGnrlzdPrssSchchtr_ (1)%fractionSubresolution(massParent,criticalOverdensity_,time,massResolution,node)
     call Assert('Branching probability rate',branchingProbabilityRate,branchingProbabilityRateTarget/rootVarianceParent,relTol=2.0d-3)
     call Assert('Accretion rate'            ,accretionRate           ,accretionRateTarget           /rootVarianceParent,relTol=2.0d-3)
     call Unit_Tests_End_Group  (                                             )
     call Unit_Tests_Begin_Group("Generalized Press-Schechter general barrier branching rates")
     smoothAccretionRateTargetGeneral     =mergerTreeBranchingProbabilityGnrlzdPrssSchchtr_(4)%fractionSubresolution(massParent,criticalOverdensity_,time,1.0d-1*massResolution,node)
     accretionRateTargetGeneral           =mergerTreeBranchingProbabilityGnrlzdPrssSchchtr_(2)%fractionSubresolution(massParent,criticalOverdensity_,time,       massResolution,node)
     branchingProbabilityRateTargetGeneral=mergerTreeBranchingProbabilityGnrlzdPrssSchchtr_(2)%probability          (massParent,criticalOverdensity_,time,       massResolution,node)
     smoothAccretionRate                  =mergerTreeBranchingProbabilityGnrlzdPrssSchchtr_(5)%fractionSubresolution(massParent,criticalOverdensity_,time,1.0d-1*massResolution,node)
     accretionRate                        =mergerTreeBranchingProbabilityGnrlzdPrssSchchtr_(3)%fractionSubresolution(massParent,criticalOverdensity_,time,       massResolution,node)
     branchingProbabilityRate             =mergerTreeBranchingProbabilityGnrlzdPrssSchchtr_(3)%probability          (massParent,criticalOverdensity_,time,       massResolution,node)
     call Assert('Branching probability rate',branchingProbabilityRate,branchingProbabilityRateTargetGeneral,relTol=2.5d-2)
     call Assert('Accretion rate'            ,accretionRate           ,accretionRateTargetGeneral           ,relTol=2.5d-2)
     call Assert('Smooth accretion rate'     ,smoothAccretionRate     ,smoothAccretionRateTargetGeneral     ,relTol=2.5d-2)
     call Unit_Tests_End_Group  (                                             )
     call Unit_Tests_End_Group  (                                             )
     call Unit_Tests_Begin_Group("First crossing distribution (numerical)"    )
     varianceParent    =cosmologicalMassVarianceFilteredPower_%rootVariance(massParent    ,time)**2
     varianceResolution=cosmologicalMassVarianceFilteredPower_%rootVariance(massResolution,time)**2
     do j=1,countVariance
        variance_                          (j)=+varianceResolution    &
             &                                 +(                     &
             &                                   +varianceParent      &
             &                                   -varianceResolution  &
             &                                  )                     &
             &                                 *dble(j            -1) &
             &                                 /dble(countVariance  )
        branchingRateUnconstrainedAnalytic (j)=excursionSetFirstCrossingLinearBarrier_               %rate(varianceParent,variance_(j),time,node)
        branchingRateUnconstrainedNumerical(j)=excursionSetFirstCrossingFarahiMidpoint_              %rate(varianceParent,variance_(j),time,node)
        branchingRateConstrainedNumerical  (j)=excursionSetFirstCrossingFarahiMidpointBrownianBridge_%rate(varianceParent,variance_(j),time,node)
        branchingRateConstrainedAnalytic   (j)=excursionSetFirstCrossingLinearBarrierBrownianBridge_ %rate(varianceParent,variance_(j),time,node)
     end do
     errorMaximumUnconstrained=maxval(                                                                                                                                         &
          &                                abs(branchingRateUnconstrainedNumerical-branchingRateUnconstrainedAnalytic)/branchingRateUnconstrainedAnalytic                      &
          &                          )
     errorMaximumConstrained  =maxval(                                                                                                                                         &
          &                                abs(branchingRateConstrainedNumerical  -branchingRateConstrainedAnalytic  )/branchingRateConstrainedAnalytic                      , &
          &                           mask= variance_                         < varianceConstrained*linearGrowthCollisionlessMatter_%value(time)**2*(1.0d0-toleranceVariance)  &
          &                                .and.                                                                                                                               &
          &                                 branchingRateConstrainedAnalytic > 0.0d0                                                                                           &
          &                          )     
     ! Tolerances here are set from the measured spread of these metrics, not from the accuracy actually achieved. The
     ! unconstrained metric is not reproducible: two runs of the same binary differ by of order one percent, and the spread
     ! between build environments is larger again, so a limit set just above the achieved value flaps. The constrained metric is
     ! instead controlled by `toleranceVariance` above, which excludes the region in which the analytic solution diverges.
     call Assert('Unconstrained case',errorMaximumUnconstrained,3.6d-2,compareLessThan)
     call Assert('Constrained case'  ,errorMaximumConstrained  ,3.3d-2,compareLessThan)
     call Unit_Tests_End_Group()
     call Unit_Tests_End_Group()
  end do

  ! ==================================================================================================================================
  ! Analytic identities of the Parkinson-Cole-Helly branching probability class.
  !
  ! For an n=-1 power spectrum filtered with a sharp k-space window function σ(M) ∝ M^{-1/3}, so that α ≡ dlnσ/dlnM = -1/3 exactly
  ! and, crucially, is independent of mass. That constancy makes several identities exact which hold only approximately for a
  ! general power spectrum, and so provides analytic targets for methods which otherwise have none.
  ! ==================================================================================================================================
  timeAnalytic           =cosmologyFunctionsMatterLambda_                              %cosmicTime  (1.0d0                                             )
  ! Use the critical overdensity evaluated at this epoch, rather than the rescaled value used above, so that the epoch recovered
  ! from it by `timeOfCollapse` - which is how the mass sampling algorithm obtains its epoch - is this same epoch. Every other
  ! method instead takes the epoch from its `time` argument.
  deltaAnalytic          =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_%value       (time=timeAnalytic,mass=massParent      ,node=node )
  sigmaParentAnalytic    =cosmologicalMassVarianceFilteredPower_                       %rootVariance(     massParent      ,timeAnalytic                )
  sigmaHalfAnalytic      =cosmologicalMassVarianceFilteredPower_                       %rootVariance(0.5d0*massParent     ,timeAnalytic                )
  call Unit_Tests_Begin_Group("Parkinson-Cole-Helly analytic identities (power-law σ(M))")
  ! Confirm the premise on which everything below rests.
  call Unit_Tests_Begin_Group("Power-law premise"                                        )
  do i=1,countAlpha
     call cosmologicalMassVarianceFilteredPower_%rootVarianceAndLogarithmicGradient(massAlpha(i),timeAnalytic,sigmaScratch(i),alphaMeasured(i))
  end do
  call Assert('α = dlnσ/dlnM = -⅓, independent of mass',alphaMeasured    ,spread(-1.0d0/3.0d0,1,countAlpha),relTol=1.0d-4)
  call Assert('σ(M/2)/σ(M) = 2^⅓'                      ,sigmaHalfAnalytic,sigmaParentAnalytic*2.0d0**(1.0d0/3.0d0),relTol=3.0d-5)
  call Unit_Tests_End_Group  (                                                           )

  ! ----------------------------------------------------------------------------------------------------------------------------------
  ! Bounds on the branching probability.
  !
  ! The bounds are obtained by holding the M′ appearing in the M/M′ factor of the branching integrand fixed - at the mass
  ! resolution for the upper bound, and at half the parent mass for the lower bound - which renders the integral analytic in terms
  ! of hypergeometric functions. When α is constant that substitution is not an approximation at all: M/M′ = (σ′/σ)^{1/α} exactly,
  ! and the factor (σ_res/σ_p)^{1/α} M/M_res carried by the CDM-assumptions form of the bound collapses to unity. Both bounds
  ! therefore reduce to the branching probability itself, and must equal the value obtained by integrating it numerically.
  !
  ! This is a strong test: a single comparison exercises the hypergeometric evaluations, both of their tabulations, the effective
  ! exponent γ₁-1/α, and the numerical integrator - against each other.
  ! ----------------------------------------------------------------------------------------------------------------------------------
  call Unit_Tests_Begin_Group("Probability bounds are exact for constant α")
  do i=1,countIdentity
     branchingIdentity_(i)=mergerTreeBranchingProbabilityParkinsonColeHelly(                                                                             &
          &                                                                 G0                       =G0Identity                                  (i) , &
          &                                                                 gamma1                   =gamma1Identity                              (i) , &
          &                                                                 gamma2                   =gamma2Identity                              (i) , &
          &                                                                 accuracyFirstOrder       =1.0d-1                                          , &
          &                                                                 precisionHypergeometric  =1.0d-6                                          , &
          &                                                                 hypergeometricTabulate   =tabulateIdentity                            (i) , &
          &                                                                 cdmAssumptions           =.true.                                          , &
          &                                                                 tolerateRoundOffErrors   =.false.                                         , &
          &                                                                 cosmologicalMassVariance_=cosmologicalMassVarianceFilteredPower_           , &
          &                                                                 criticalOverdensity_     =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_ &
          &                                                                )
     do j=1,countRatio
        massResolution_       =massParent/ratioResolution(j)
        probabilityIdentity(j)=branchingIdentity_(i)%probability     (massParent,deltaAnalytic,timeAnalytic,massResolution_                              ,node)
        boundLowerIdentity (j)=branchingIdentity_(i)%probabilityBound(massParent,deltaAnalytic,timeAnalytic,massResolution_,mergerTreeBranchingBoundLower,node)
        boundUpperIdentity (j)=branchingIdentity_(i)%probabilityBound(massParent,deltaAnalytic,timeAnalytic,massResolution_,mergerTreeBranchingBoundUpper,node)
     end do
     write (label,'(a,f5.2)') "γ₁ = ",gamma1Identity(i)
     if (tabulateIdentity(i)) then
        call Unit_Tests_Begin_Group(trim(label)//"; tabulated ₂F₁")
        ! The upper bound is interpolated in a table of thirty points per decade in mass, which limits the agreement here; the
        ! lower bound is computed directly even when tabulation is requested, and so is limited only by the integrator.
        call Assert('lower bound = probability',boundLowerIdentity,probabilityIdentity,relTol=toleranceIdentityDirect   )
        call Assert('upper bound = probability',boundUpperIdentity,probabilityIdentity,relTol=toleranceIdentityTabulated)
     else
        call Unit_Tests_Begin_Group(trim(label)//"; computed ₂F₁"  )
        ! Both bounds are computed directly here, so the agreement is limited by the tolerance of the integrator used to evaluate
        ! the branching probability against which they are compared.
        call Assert('lower bound = probability',boundLowerIdentity,probabilityIdentity,relTol=toleranceIdentityDirect   )
        call Assert('upper bound = probability',boundUpperIdentity,probabilityIdentity,relTol=toleranceIdentityDirect   )
     end if
     call Assert('upper bound = lower bound'   ,boundUpperIdentity,boundLowerIdentity ,relTol=toleranceIdentityTabulated)
     call Unit_Tests_End_Group  (                                 )
  end do
  call Unit_Tests_End_Group  (                                    )

  ! ----------------------------------------------------------------------------------------------------------------------------------
  ! Bounds bound the probability.
  !
  ! Where the identity above does not apply the bounds must at least do what their name promises. This is the property on which
  ! the tree building algorithm depends: branching events are drawn at the rate given by the upper bound and then rejected using
  ! the exact rate, so an upper bound which is not one silently biases every tree built. It is checked here both without the CDM
  ! assumptions, where the bounds are genuinely inequalities, and for γ₁=0, where the effective exponent γ₁-1/α lands exactly on
  ! the singular point of the hypergeometric function and the evaluation must fall back to the more robust form.
  ! ----------------------------------------------------------------------------------------------------------------------------------
  call Unit_Tests_Begin_Group("Probability bounds bracket the probability")
  branchingOrdering_=mergerTreeBranchingProbabilityParkinsonColeHelly(                                                                                                     &
       &                                                              G0                       =0.57d0                                                                   , &
       &                                                              gamma1                   =0.38d0                                                                   , &
       &                                                              gamma2                   =-0.01d0                                                                  , &
       &                                                              accuracyFirstOrder       =1.0d-1                                                                   , &
       &                                                              precisionHypergeometric  =1.0d-6                                                                   , &
       &                                                              hypergeometricTabulate   =.false.                                                                  , &
       &                                                              cdmAssumptions           =.false.                                                                  , &
       &                                                              tolerateRoundOffErrors   =.false.                                                                  , &
       &                                                              cosmologicalMassVariance_=cosmologicalMassVarianceFilteredPower_                                    , &
       &                                                              criticalOverdensity_     =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_              &
       &                                                             )
  ! γ₁=0 places the effective exponent γ₁-1/α at exactly 3, an odd integer, at which the third parameter of the hypergeometric
  ! function vanishes and its evaluation fails. Tabulation is disabled here because only the direct evaluation carries the
  ! fallback to the more robust bound; the tabulation of the upper bound has no such fallback.
  branchingPole_    =mergerTreeBranchingProbabilityParkinsonColeHelly(                                                                                                     &
       &                                                              G0                       =1.0d0                                                                    , &
       &                                                              gamma1                   =0.0d0                                                                    , &
       &                                                              gamma2                   =0.0d0                                                                    , &
       &                                                              accuracyFirstOrder       =1.0d-1                                                                   , &
       &                                                              precisionHypergeometric  =1.0d-6                                                                   , &
       &                                                              hypergeometricTabulate   =.false.                                                                  , &
       &                                                              cdmAssumptions           =.true.                                                                   , &
       &                                                              tolerateRoundOffErrors   =.false.                                                                  , &
       &                                                              cosmologicalMassVariance_=cosmologicalMassVarianceFilteredPower_                                    , &
       &                                                              criticalOverdensity_     =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_              &
       &                                                             )
  do j=1,countRatio
     massResolution_       =massParent/ratioResolution(j)
     probabilityIdentity(j)=branchingOrdering_%probability     (massParent,deltaAnalytic,timeAnalytic,massResolution_                              ,node)
     boundLowerIdentity (j)=branchingOrdering_%probabilityBound(massParent,deltaAnalytic,timeAnalytic,massResolution_,mergerTreeBranchingBoundLower,node)
     boundUpperIdentity (j)=branchingOrdering_%probabilityBound(massParent,deltaAnalytic,timeAnalytic,massResolution_,mergerTreeBranchingBoundUpper,node)
  end do
  call Unit_Tests_Begin_Group("No CDM assumptions"                     )
  call Assert('lower bound ≤ probability',boundLowerIdentity ,probabilityIdentity,compareLessThanOrEqual)
  call Assert('probability ≤ upper bound',probabilityIdentity,boundUpperIdentity ,compareLessThanOrEqual)
  call Unit_Tests_End_Group  (                                         )
  do j=1,countRatio
     massResolution_       =massParent/ratioResolution(j)
     probabilityIdentity(j)=branchingPole_%probability     (massParent,deltaAnalytic,timeAnalytic,massResolution_                              ,node)
     boundLowerIdentity (j)=branchingPole_%probabilityBound(massParent,deltaAnalytic,timeAnalytic,massResolution_,mergerTreeBranchingBoundLower,node)
     boundUpperIdentity (j)=branchingPole_%probabilityBound(massParent,deltaAnalytic,timeAnalytic,massResolution_,mergerTreeBranchingBoundUpper,node)
  end do
  ! γ₁=0 places the effective exponent exactly on the odd integer 3. That was once a singular point of the evaluation, at which
  ! the bounds had to fall back to a looser form; the integral is now evaluated between its limits directly, which has no such
  ! point, so the bounds here are as exact as anywhere else and are checked against the same identity.
  call Unit_Tests_Begin_Group("Formerly singular effective exponent (γ₁-1/α = 3)")
  call Assert('lower bound = probability',boundLowerIdentity,probabilityIdentity,relTol=toleranceIdentityDirect)
  call Assert('upper bound = probability',boundUpperIdentity,probabilityIdentity,relTol=toleranceIdentityDirect)
  call Unit_Tests_End_Group  (                                                   )
  call Unit_Tests_End_Group  (                                         )

  ! ----------------------------------------------------------------------------------------------------------------------------------
  ! Subresolution accretion rate against elementary closed forms.
  !
  ! The subresolution accretion rate is √(2/π) G₀ (δ_p/σ_p)^{γ₂}/σ_p multiplied by a factor which, in general, requires a
  ! hypergeometric function. At γ₁ = 0, -1 and -3 that function degenerates to an elementary one - see `factorAnalytic` below - so
  ! the rate can be checked against a closed form containing no special functions at all.
  !
  ! The sweep is made in x ≡ σ(M_res)/σ(M) - 1, which is the argument on which that factor depends and which approaches the
  ! singular point of the hypergeometric function as x → 0. It spans fourteen decades and straddles x = 1/10, at which the
  ! evaluation switches from a direct one to a series in 1-z, and the seed bounds of the tabulation at x = 10⁻⁹ and x = 12.5. Note
  ! that the expected value is built from σ evaluated by the same object used by the code, so that any inaccuracy in σ itself
  ! cancels and does not limit the comparison.
  ! ----------------------------------------------------------------------------------------------------------------------------------
  call Unit_Tests_Begin_Group("Subresolution accretion rate has an elementary closed form")
  do i=1,countClosedForm
     branchingClosedForm_(i)=mergerTreeBranchingProbabilityParkinsonColeHelly(                                                                             &
          &                                                                   G0                       =1.0d0                                           , &
          &                                                                   gamma1                   =gamma1ClosedForm                            (i) , &
          &                                                                   gamma2                   =0.0d0                                           , &
          &                                                                   accuracyFirstOrder       =1.0d-1                                          , &
          &                                                                   precisionHypergeometric  =1.0d-6                                          , &
          &                                                                   hypergeometricTabulate   =tabulateClosedForm                          (i) , &
          &                                                                   cdmAssumptions           =.false.                                         , &
          &                                                                   tolerateRoundOffErrors   =.false.                                         , &
          &                                                                   cosmologicalMassVariance_=cosmologicalMassVarianceFilteredPower_           , &
          &                                                                   criticalOverdensity_     =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_ &
          &                                                                  )
     do j=1,countX
        massResolution_        =massParent/(1.0d0+xTarget(j))**3
        sigmaResolutionAnalytic=cosmologicalMassVarianceFilteredPower_%rootVariance(massResolution_,timeAnalytic)
        xActual         (j)    =sigmaResolutionAnalytic/sigmaParentAnalytic-1.0d0
        fractionActual  (j)    =branchingClosedForm_(i)%fractionSubresolution(massParent,deltaAnalytic,timeAnalytic,massResolution_,node)
        fractionExpected(j)    =+sqrt(2.0d0/Pi)                                              &
             &                  /sigmaParentAnalytic                                         &
             &                  *factorAnalytic(xActual(j),gamma1ClosedForm(i))
     end do
     write (label,'(a,f5.2)') "γ₁ = ",gamma1ClosedForm(i)
     if (tabulateClosedForm(i)) then
        ! The tabulated evaluation is exact at its abscissae but is interpolated between them at twenty points per decade, which
        ! is what limits the agreement here.
        call Unit_Tests_Begin_Group(trim(label)//"; tabulated ₂F₁")
        call Assert('subresolution accretion rate',fractionActual,fractionExpected,relTol=toleranceClosedForm(i))
     else
        ! The direct evaluation loses precision as x → 0, where the argument of the hypergeometric function approaches its
        ! singular point and 1-z is formed by a cancelling subtraction from unity. Only the tabulated route uses the series in
        ! 1-z which avoids that cancellation, so this is the looser of the two at small x.
        call Unit_Tests_Begin_Group(trim(label)//"; computed ₂F₁"  )
        call Assert('subresolution accretion rate',fractionActual,fractionExpected,relTol=toleranceClosedForm(i))
     end if
     call Unit_Tests_End_Group  (                                   )
  end do
  call Unit_Tests_End_Group  (                                      )

  ! ----------------------------------------------------------------------------------------------------------------------------------
  ! Mutual consistency of the branching rate, probability, and branch mass.
  !
  ! These identities follow from the definitions alone and so hold for any power spectrum: the branching probability is the
  ! integral of the branching rate over all resolved progenitor masses; the rate is symmetric under exchange of the two
  ! progenitors; and the branch mass returned for a given probability is the mass at which the partial integral of the rate
  ! reaches it. They tie together three methods which have no analytic solution of their own.
  ! ----------------------------------------------------------------------------------------------------------------------------------
  call Unit_Tests_Begin_Group("Rate, probability, and branch mass are mutually consistent")
  branchingGeneric_=mergerTreeBranchingProbabilityParkinsonColeHelly(                                                                                                     &
       &                                                             G0                       =0.57d0                                                                   , &
       &                                                             gamma1                   =0.38d0                                                                   , &
       &                                                             gamma2                   =-0.01d0                                                                  , &
       &                                                             accuracyFirstOrder       =1.0d-1                                                                   , &
       &                                                             precisionHypergeometric  =1.0d-6                                                                   , &
       &                                                             hypergeometricTabulate   =.false.                                                                  , &
       &                                                             cdmAssumptions           =.false.                                                                  , &
       &                                                             tolerateRoundOffErrors   =.false.                                                                  , &
       &                                                             cosmologicalMassVariance_=cosmologicalMassVarianceFilteredPower_                                    , &
       &                                                             criticalOverdensity_     =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_              &
       &                                                            )
  massResolution_ =massParent/1.0d3
  probabilityFull =branchingGeneric_%probability(massParent,deltaAnalytic,timeAnalytic,massResolution_,node)
  call Assert('∫ rate dM′ over resolved progenitors = probability',rateIntegral(branchingGeneric_,massResolution_,0.5d0*massParent),probabilityFull,relTol=toleranceRateIntegral)
  massBranchTest  =massParent/1.0d2
  call Assert('rate(M-m) = rate(m)'                              ,branchingGeneric_%rate(massParent,deltaAnalytic,timeAnalytic,massParent-massBranchTest,node),branchingGeneric_%rate(massParent,deltaAnalytic,timeAnalytic,massBranchTest,node),relTol=1.0d-12)
  ! The branch mass is found by inverting the partial integral of the rate, so feeding it a partial integral computed
  ! independently must return the mass at which that integral was truncated.
  do j=1,countInverse
     massTargetCurrent    =massResolution_*(0.5d0*massParent/massResolution_)**fractionInverse(j)
     probabilityPartial   =rateIntegral(branchingGeneric_,massResolution_,massTargetCurrent)
     massBranchTarget  (j)=massTargetCurrent
     massBranchReturned(j)=branchingGeneric_%massBranch(massParent,deltaAnalytic,timeAnalytic,massResolution_,probabilityPartial,randomNumberGenerator_,node)
  end do
  call Assert('branch mass inverts ∫ rate dM′'                   ,massBranchReturned,massBranchTarget,relTol=toleranceMassBranchInverse)
  ! Supplying the full branching probability must return exactly half the parent mass, the largest branch mass permitted.
  call Assert('branch mass at full probability = M/2'            ,branchingGeneric_%massBranch(massParent,deltaAnalytic,timeAnalytic,massResolution_,probabilityFull,randomNumberGenerator_,node),0.5d0*massParent,relTol=toleranceMassBranchInverse)
  call Unit_Tests_End_Group  (                                   )

  ! ----------------------------------------------------------------------------------------------------------------------------------
  ! The branch mass sampler reproduces the branching rate.
  !
  ! Under CDM assumptions branch masses are drawn by rejection sampling from the function S(q) of Parkinson et al. (2008) rather
  ! than by inverting the integral of the rate, so the two are entirely independent implementations of the same distribution. The
  ! cumulative distribution of a sample drawn from the former is compared here with that implied by the latter.
  ! ----------------------------------------------------------------------------------------------------------------------------------
  call Unit_Tests_Begin_Group("Branch mass sampler reproduces the branching rate")
  branchingSampler_     =mergerTreeBranchingProbabilityParkinsonColeHelly(                                                                                                &
       &                                                             G0                       =0.57d0                                                                   , &
       &                                                             gamma1                   =0.38d0                                                                   , &
       &                                                             gamma2                   =-0.01d0                                                                  , &
       &                                                             accuracyFirstOrder       =1.0d-1                                                                   , &
       &                                                             precisionHypergeometric  =1.0d-6                                                                   , &
       &                                                             hypergeometricTabulate   =.false.                                                                  , &
       &                                                             cdmAssumptions           =.true.                                                                   , &
       &                                                             tolerateRoundOffErrors   =.false.                                                                  , &
       &                                                             cosmologicalMassVariance_=cosmologicalMassVarianceFilteredPower_                                    , &
       &                                                             criticalOverdensity_     =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_              &
       &                                                            )
  branchingSamplerGeneric_=mergerTreeBranchingProbabilityParkinsonColeHelly(                                                                                              &
       &                                                             G0                       =0.57d0                                                                   , &
       &                                                             gamma1                   =0.38d0                                                                   , &
       &                                                             gamma2                   =-0.01d0                                                                  , &
       &                                                             accuracyFirstOrder       =1.0d-1                                                                   , &
       &                                                             precisionHypergeometric  =1.0d-6                                                                   , &
       &                                                             hypergeometricTabulate   =.false.                                                                  , &
       &                                                             cdmAssumptions           =.false.                                                                  , &
       &                                                             tolerateRoundOffErrors   =.false.                                                                  , &
       &                                                             cosmologicalMassVariance_=cosmologicalMassVarianceFilteredPower_                                    , &
       &                                                             criticalOverdensity_     =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_              &
       &                                                            )
  ! The cumulative distribution of the branch mass requires no special functions of its own: the branching probability evaluated
  ! for a resolution m is precisely the integral of the branching rate from m up to half the parent mass, so the fraction of
  ! branches below m is one minus the ratio of that to the same quantity evaluated at the true resolution.
  do j=1,countGrid
     massGrid   (j)=massResolution_*(0.5d0*massParent/massResolution_)**(dble(j)/dble(countGrid+1))
     cdfAnalytic(j)=+1.0d0                                                                                                  &
          &         -branchingGeneric_%probability(massParent,deltaAnalytic,timeAnalytic,massGrid(j),node)/probabilityFull
  end do
  randomNumberGenerator_=randomNumberGeneratorGSL(8213_c_long)
  ! Under the CDM assumptions the branch mass is drawn by rejection sampling from the function S(q) of Parkinson et al. (2008),
  ! which does not reference the branching rate at all - so this compares two entirely independent constructions.
  do k=1,countSample
     massSample(k)=branchingSampler_       %massBranch(massParent,deltaAnalytic,timeAnalytic,massResolution_,0.0d0,randomNumberGenerator_,node)
  end do
  do j=1,countGrid
     cdfSampled(j)=dble(count(massSample < massGrid(j)))/dble(countSample)
  end do
  call Assert('Kolmogorov-Smirnov statistic, S(q) rejection sampler',maxval(abs(cdfSampled-cdfAnalytic)),toleranceKolmogorovSmirnov,compareLessThan)
  ! Without them the branch mass is instead found by inverting the integral of the rate, so a uniform deviate scaled by the total
  ! branching probability samples the same distribution by a wholly different route.
  do k=1,countSample
     massSample(k)=branchingSamplerGeneric_%massBranch(massParent,deltaAnalytic,timeAnalytic,massResolution_,randomNumberGenerator_%uniformSample()*probabilityFull,randomNumberGenerator_,node)
  end do
  do j=1,countGrid
     cdfSampled(j)=dble(count(massSample < massGrid(j)))/dble(countSample)
  end do
  call Assert('Kolmogorov-Smirnov statistic, rate inversion sampler',maxval(abs(cdfSampled-cdfAnalytic)),toleranceKolmogorovSmirnov,compareLessThan)
  call Unit_Tests_End_Group  (                                   )

  ! ----------------------------------------------------------------------------------------------------------------------------------
  ! Reductions of the PCH+ class.
  !
  ! PCH+ adds three parameters to the Parkinson-Cole-Helly algorithm. With all three set to zero it must reproduce that algorithm
  ! exactly; with γ₃ zero but γ₄ and γ₅ non-zero it must reproduce it multiplied by a factor which depends only on the parent halo
  ! and so can be written down; and with γ₃ non-zero the bounds identity above must continue to hold, since the substitution on
  ! which it rests is unaffected by γ₃. Together these check the fused form in which the class evaluates the progenitor mass
  ! function against the unfused form used by its parent.
  ! ----------------------------------------------------------------------------------------------------------------------------------
  call Unit_Tests_Begin_Group("PCH+ reductions")
  branchingReference_=mergerTreeBranchingProbabilityParkinsonColeHelly(                                                                                                   &
       &                                                             G0                       =0.57d0                                                                   , &
       &                                                             gamma1                   =0.38d0                                                                   , &
       &                                                             gamma2                   =-0.01d0                                                                  , &
       &                                                             accuracyFirstOrder       =1.0d-1                                                                   , &
       &                                                             precisionHypergeometric  =1.0d-6                                                                   , &
       &                                                             hypergeometricTabulate   =.false.                                                                  , &
       &                                                             cdmAssumptions           =.true.                                                                   , &
       &                                                             tolerateRoundOffErrors   =.false.                                                                  , &
       &                                                             cosmologicalMassVariance_=cosmologicalMassVarianceFilteredPower_                                    , &
       &                                                             criticalOverdensity_     =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_              &
       &                                                            )
  do i=1,3
     branchingPCHPlus_(i)=mergerTreeBranchingProbabilityPCHPlus(                                                                                                          &
          &                                                     G0                       =0.57d0                                                                       , &
          &                                                     gamma1                   =0.38d0                                                                       , &
          &                                                     gamma2                   =-0.01d0                                                                      , &
          &                                                     gamma3                   =gamma3PCHPlus                                                          (i) , &
          &                                                     gamma4                   =gamma4PCHPlus                                                          (i) , &
          &                                                     gamma5                   =gamma5PCHPlus                                                          (i) , &
          &                                                     accuracyFirstOrder       =1.0d-1                                                                       , &
          &                                                     precisionHypergeometric  =1.0d-6                                                                       , &
          &                                                     hypergeometricTabulate   =.false.                                                                      , &
          &                                                     cdmAssumptions           =.true.                                                                       , &
          &                                                     tolerateRoundOffErrors   =.false.                                                                      , &
          &                                                     cosmologicalMassVariance_=cosmologicalMassVarianceFilteredPower_                                        , &
          &                                                     criticalOverdensity_     =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_                 , &
          &                                                     linearGrowth_            =linearGrowthCollisionlessMatter_                                               &
          &                                                    )
  end do
  massResolution_=massParent/1.0d3
  call Unit_Tests_Begin_Group("γ₃ = γ₄ = γ₅ = 0 reduces to Parkinson-Cole-Helly")
  call Assert('branching probability'      ,branchingPCHPlus_(1)%probability          (massParent,deltaAnalytic,timeAnalytic,massResolution_                              ,node),branchingReference_%probability          (massParent,deltaAnalytic,timeAnalytic,massResolution_                              ,node),relTol=toleranceReduction)
  call Assert('subresolution accretion'    ,branchingPCHPlus_(1)%fractionSubresolution(massParent,deltaAnalytic,timeAnalytic,massResolution_                              ,node),branchingReference_%fractionSubresolution(massParent,deltaAnalytic,timeAnalytic,massResolution_                              ,node),relTol=toleranceReduction)
  call Assert('lower bound'                ,branchingPCHPlus_(1)%probabilityBound     (massParent,deltaAnalytic,timeAnalytic,massResolution_,mergerTreeBranchingBoundLower,node),branchingReference_%probabilityBound     (massParent,deltaAnalytic,timeAnalytic,massResolution_,mergerTreeBranchingBoundLower,node),relTol=toleranceReduction)
  call Assert('upper bound'                ,branchingPCHPlus_(1)%probabilityBound     (massParent,deltaAnalytic,timeAnalytic,massResolution_,mergerTreeBranchingBoundUpper,node),branchingReference_%probabilityBound     (massParent,deltaAnalytic,timeAnalytic,massResolution_,mergerTreeBranchingBoundUpper,node),relTol=toleranceReduction)
  call Assert('branching rate'             ,branchingPCHPlus_(1)%rate                 (massParent,deltaAnalytic,timeAnalytic,massBranchTest                              ,node),branchingReference_%rate                 (massParent,deltaAnalytic,timeAnalytic,massBranchTest                              ,node),relTol=toleranceReduction)
  call Unit_Tests_End_Group  (                                                        )
  ! With γ₃=0 the parameters γ₄ and γ₅ multiply the branching probability by a factor built from the logarithmic derivative of the
  ! growth factor and the logarithmic gradient of σ at the parent mass, both of which are constant across the integral.
  modifierExpected=+     linearGrowthCollisionlessMatter_      %logarithmicDerivativeExpansionFactor(             timeAnalytic) **0.5d0 &
       &           *abs(cosmologicalMassVarianceFilteredPower_%rootVarianceLogarithmicGradient     (massParent  ,timeAnalytic))**0.3d0
  call Unit_Tests_Begin_Group("γ₄ and γ₅ multiply by a constant factor"        )
  call Assert('branching probability'      ,branchingPCHPlus_(2)%probability          (massParent,deltaAnalytic,timeAnalytic,massResolution_                              ,node),modifierExpected*branchingReference_%probability          (massParent,deltaAnalytic,timeAnalytic,massResolution_                              ,node),relTol=toleranceReduction)
  call Assert('subresolution accretion'    ,branchingPCHPlus_(2)%fractionSubresolution(massParent,deltaAnalytic,timeAnalytic,massResolution_                              ,node),modifierExpected*branchingReference_%fractionSubresolution(massParent,deltaAnalytic,timeAnalytic,massResolution_                              ,node),relTol=toleranceReduction)
  call Unit_Tests_End_Group  (                                                        )
  ! The bounds identity is unaffected by γ₃, which changes the first parameter of the hypergeometric function but not the
  ! substitution which makes the bounds exact for constant α.
  do j=1,countRatio
     massResolution_       =massParent/ratioResolution(j)
     probabilityIdentity(j)=branchingPCHPlus_(3)%probability     (massParent,deltaAnalytic,timeAnalytic,massResolution_                              ,node)
     boundLowerIdentity (j)=branchingPCHPlus_(3)%probabilityBound(massParent,deltaAnalytic,timeAnalytic,massResolution_,mergerTreeBranchingBoundLower,node)
     boundUpperIdentity (j)=branchingPCHPlus_(3)%probabilityBound(massParent,deltaAnalytic,timeAnalytic,massResolution_,mergerTreeBranchingBoundUpper,node)
  end do
  call Unit_Tests_Begin_Group("γ₃ ≠ 0 preserves the exactness of the bounds"   )
  call Assert('lower bound = probability'  ,boundLowerIdentity,probabilityIdentity,relTol=toleranceIdentityDirect)
  call Assert('upper bound = probability'  ,boundUpperIdentity,probabilityIdentity,relTol=toleranceIdentityDirect)
  call Unit_Tests_End_Group  (                                                        )
  call Unit_Tests_End_Group  (                                                        )

  ! ----------------------------------------------------------------------------------------------------------------------------------
  ! Maximum step, limiting cases, and sentinel values.
  !
  ! The maximum step is checked on an object on which nothing else has been called, since it must not depend on state left behind
  ! by some earlier call. The remaining assertions pin the behavior at the edges of the domain: a halo at twice the mass
  ! resolution has no resolved progenitors, and a mass resolution above the halo mass admits no subresolution accretion rate at
  ! all, which is signalled by a negative return value on which the tree builder branches.
  ! ----------------------------------------------------------------------------------------------------------------------------------
  call Unit_Tests_Begin_Group("Maximum step, limiting cases, and sentinels")
  branchingFresh_=mergerTreeBranchingProbabilityParkinsonColeHelly(                                                                                                       &
       &                                                             G0                       =0.57d0                                                                   , &
       &                                                             gamma1                   =0.38d0                                                                   , &
       &                                                             gamma2                   =-0.01d0                                                                  , &
       &                                                             accuracyFirstOrder       =1.0d-1                                                                   , &
       &                                                             precisionHypergeometric  =1.0d-6                                                                   , &
       &                                                             hypergeometricTabulate   =.false.                                                                  , &
       &                                                             cdmAssumptions           =.false.                                                                  , &
       &                                                             tolerateRoundOffErrors   =.false.                                                                  , &
       &                                                             cosmologicalMassVariance_=cosmologicalMassVarianceFilteredPower_                                    , &
       &                                                             criticalOverdensity_     =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_              &
       &                                                            )
  massResolution_=massParent/1.0d3
  stepExpected   =+1.0d-1                       &
       &          *sqrt(                        &
       &                +2.0d0                  &
       &                *(                      &
       &                  +sigmaHalfAnalytic**2 &
       &                  -sigmaParentAnalytic**2 &
       &                 )                      &
       &               )
  call Assert('maximum step, object not previously used'      ,branchingFresh_%stepMaximum(massParent,deltaAnalytic,timeAnalytic,massResolution_    ),stepExpected,relTol=1.0d-12)
  call Assert('maximum step is unbounded below 2M_res'        ,branchingFresh_%stepMaximum(massParent,deltaAnalytic,timeAnalytic,0.5d0*massParent   ),1.0d10      ,relTol=1.0d-12)
  call Assert('probability vanishes at M_res = M/2'           ,branchingFresh_%probability(massParent,deltaAnalytic,timeAnalytic,0.5d0*massParent   ,node),0.0d0 ,absTol=0.0d0  )
  call Assert('lower bound vanishes at M_res = M/2'           ,branchingFresh_%probabilityBound(massParent,deltaAnalytic,timeAnalytic,0.5d0*massParent,mergerTreeBranchingBoundLower,node),0.0d0,absTol=0.0d0)
  call Assert('upper bound vanishes at M_res = M/2'           ,branchingFresh_%probabilityBound(massParent,deltaAnalytic,timeAnalytic,0.5d0*massParent,mergerTreeBranchingBoundUpper,node),0.0d0,absTol=0.0d0)
  call Assert('subresolution accretion signals M_res ≥ M'     ,branchingFresh_%fractionSubresolution(massParent,deltaAnalytic,timeAnalytic,2.0d0*massParent,node),-1.0d0,relTol=1.0d-12)
  ! An extreme dynamic range between the halo and resolution masses, well beyond that used in practice.
  call Assert('probability is positive at M/M_res = 10¹⁰'     ,branchingFresh_%probability(massParent,deltaAnalytic,timeAnalytic,massParent/1.0d10   ,node),0.0d0 ,compareGreaterThan)
  call Unit_Tests_End_Group  (                                                        )
  call Unit_Tests_End_Group  (                                                        )

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

contains

  double precision function factorAnalytic(x,gamma)
    ! Elementary closed forms of the factor
    !
    !   (1+x)^{γ-1} ₂F₁(3/2,(1-γ)/2;(3-γ)/2;1/(1+x)²) / (1-γ)
    !
    ! which appears in the subresolution accretion rate. Writing u=1+x, the hypergeometric function admits the representation
    ! ₂F₁(3/2,b;b+1;z) = b z^{-b} B_z(b,-1/2) with b=(1-γ)/2, and the incomplete beta function reduces to elementary functions
    ! whenever b is a positive integer or half-integer. Taking γ = 0, -1 and -3 gives b = 1/2, 1 and 2 respectively, and hence the
    ! three forms below. Note that γ=0 is the only one of the three for which the constant term of the series representation used
    ! by the code vanishes, so the other two are needed to exercise it.
    !
    ! Each is written free of cancellation. In terms of s ≡ √(u²-1)/u = √(1-ε) with ε ≡ 1/u², the γ=-1 form is 1/s-1 and the γ=-3
    ! form collapses to (1-s)²/s; written directly both lose precision as x grows - the first to first order and the second to
    ! second order, so that by x=10³ the γ=-3 form has lost some twelve digits and is no longer accurate enough to serve as a
    ! reference. Writing 1-s as ε/(1+s) removes the cancellation entirely.
    implicit none
    double precision, intent(in   ) :: x, gamma
    double precision                :: u, rootDifference, &
         &                             s, epsilon_       , &
         &                             oneMinusS

    u             =+1.0d0+x
    ! √(u²-1), written in terms of x so that it suffers no cancellation as x → 0.
    rootDifference=sqrt(x*(2.0d0+x))
    epsilon_      =+1.0d0/u**2
    s             =+rootDifference/u
    oneMinusS     =+epsilon_/(1.0d0+s)
    if      (gamma == +0.0d0) then
       factorAnalytic=+1.0d0/rootDifference
    else if (gamma == -1.0d0) then
       factorAnalytic=+oneMinusS/s
    else if (gamma == -3.0d0) then
       factorAnalytic=+oneMinusS**2/s
    else
       factorAnalytic=0.0d0
       call Error_Report('no elementary closed form is available for this γ'//{introspection:location})
    end if
    return
  end function factorAnalytic

  double precision function rateIntegral(branching_,massLower,massUpper)
    ! Integrate the branching rate over progenitor mass using the composite Simpson rule in ln M′. The integrand is smooth over
    ! the ranges used here, so a fixed rule suffices - and, being independent of the adaptive integrator used inside the branching
    ! probability itself, it provides a genuinely independent evaluation of that integral.
    implicit none
    class           (mergerTreeBranchingProbabilityClass), intent(inout), target :: branching_
    double precision                                     , intent(in   )         :: massLower   , massUpper
    integer                                              , parameter             :: countInterval=2000
    double precision                                                             :: logMassLower, logMassUpper, &
         &                                                                          deltaLogMass, weight      , &
         &                                                                          massCurrent
    integer                                                                      :: l

    logMassLower=log(massLower)
    logMassUpper=log(massUpper)
    deltaLogMass=(logMassUpper-logMassLower)/dble(countInterval)
    rateIntegral=0.0d0
    do l=0,countInterval
       if      (l == 0 .or. l == countInterval) then
          weight=1.0d0
       else if (mod(l,2) == 1                 ) then
          weight=4.0d0
       else
          weight=2.0d0
       end if
       massCurrent =exp(logMassLower+deltaLogMass*dble(l))
       ! The rate is per unit progenitor mass, so a factor of the mass converts it to a rate per unit logarithmic mass.
       rateIntegral=+rateIntegral                                                                                 &
            &       +weight                                                                                       &
            &       *branching_%rate(massParent,deltaAnalytic,timeAnalytic,massCurrent,node)                      &
            &       *massCurrent
    end do
    rateIntegral=+rateIntegral  &
         &       *deltaLogMass  &
         &       /3.0d0
    return
  end function rateIntegral

end program Tests_Merger_Tree_Branching
