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
  !!{
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
  use :: Merger_Tree_Branching               , only : mergerTreeBranchingProbabilityGnrlzdPrssSchchtr, mergerTreeBranchingProbabilityParkinsonColeHelly
  use :: Merger_Tree_Branching_Modifiers     , only : mergerTreeBranchingProbabilityModifierIdentity , mergerTreeBranchingProbabilityModifierPCHPlus
  use :: Numerical_Constants_Math            , only : Pi
  use :: Power_Spectra_Primordial            , only : powerSpectrumPrimordialPowerLaw
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredSimple
  use :: Power_Spectrum_Window_Functions     , only : powerSpectrumWindowFunctionSharpKSpace
  use :: Transfer_Functions                  , only : transferFunctionIdentity
  use :: Unit_Tests                          , only : Assert                                         , Unit_Tests_Begin_Group                                      , Unit_Tests_End_Group                 , Unit_Tests_Finish                                    , &
       &                                              compareLessThan
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
       &                                                                                                      criticalOverdensityConstrained                                =4.0d00             , toleranceVariance         =4.0d-2
  integer                                                                       , parameter                :: countVariance                                                 =1000
  double precision                                                              , dimension(countVariance) :: branchingRateUnconstrainedAnalytic                                                , branchingRateUnconstrainedNumerical            , &
       &                                                                                                      branchingRateConstrainedAnalytic                                                  , branchingRateConstrainedNumerical              , &
       &                                                                                                      variance_
  integer                                                                                                  :: i                                                                                 , j
  character       (len=16                                                      )                           :: label

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelWorking)
  ! Initialize event hooks.
  call eventsHooksInitialize()
  ! Build all objects needed for these tests.
  darkMatterParticleCDM_                                        =darkMatterParticleCDM                                        (                                                                                                       &
       &                                                                                                                      )
  cosmologyParametersSimple_                                    =cosmologyParametersSimple                                    (                                                                                                       &
       &                                                                                                                       OmegaMatter                            = 1.00d0                                                      , &
       &                                                                                                                       OmegaBaryon                            = 0.00d0                                                      , &
       &                                                                                                                       OmegaDarkEnergy                        = 0.00d0                                                      , &
       &                                                                                                                       temperatureCMB                         = 2.78d0                                                      , &
       &                                                                                                                       HubbleConstant                         =70.00d0                                                        &
       &                                                                                                                      )
  cosmologyFunctionsMatterLambda_                               =cosmologyFunctionsMatterLambda                               (                                                                                                       &
       &                                                                                                                       cosmologyParameters_                   =cosmologyParametersSimple_                                     &
       &                                                                                                                      )
  linearGrowthCollisionlessMatter_                              =linearGrowthCollisionlessMatter                              (                                                                                                       &
       &                                                                                                                       cosmologyParameters_                   =cosmologyParametersSimple_                                   , &
       &                                                                                                                       cosmologyFunctions_                    =cosmologyFunctionsMatterLambda_                                &
       &                                                                                                                      )
  powerSpectrumPrimordialPowerLaw_                              =powerSpectrumPrimordialPowerLaw                              (                                                                                                       &
       &                                                                                                                       index_                                 =-1.0d0                                                       , &
       &                                                                                                                       running                                =+0.0d0                                                       , &
       &                                                                                                                       runningRunning                         =+0.0d0                                                       , &
       &                                                                                                                       wavenumberReference                    =+1.0d0                                                       , &
       &                                                                                                                       runningSmallScalesOnly                 =.false.                                                        &
       &                                                                                                                      )
  transferFunctionIdentity_                                     =transferFunctionIdentity                                     (                                                                                                       &
       &                                                                                                                       cosmologyParameters_                   =cosmologyParametersSimple_                                   , &
       &                                                                                                                       time                                   =13.8d0                                                         &
       &                                                                                                                      )
  powerSpectrumPrimordialTransferredSimple_                     =powerSpectrumPrimordialTransferredSimple                     (                                                                                                       &
       &                                                                                                                       powerSpectrumPrimordial_               =powerSpectrumPrimordialPowerLaw_                             , &
       &                                                                                                                       transferFunction_                      =transferFunctionIdentity_                                    , &
       &                                                                                                                       linearGrowth_                          =linearGrowthCollisionlessMatter_                               &
       &                                                                                                                      )
  powerSpectrumWindowFunctionSharpKSpace_                       =powerSpectrumWindowFunctionSharpKSpace                       (                                                                                                       &
       &                                                                                                                       cosmologyParameters_                   =cosmologyParametersSimple_                                   , &
       &                                                                                                                       normalization                          =0.0d0                                                          &
       &                                                                                                                      )

  cosmologicalMassVarianceFilteredPower_                        =cosmologicalMassVarianceFilteredPower                        (                                                                                                       &
       &                                                                                                                       sigma8                                 =1.0d+0                                                       , &
       &                                                                                                                       tolerance                              =1.0d-4                                                       , &
       &                                                                                                                       toleranceTopHat                        =1.0d-4                                                       , &
       &                                                                                                                       nonMonotonicIsFatal                    =.true.                                                       , &
       &                                                                                                                       monotonicInterpolation                 =.false.                                                      , &
       &                                                                                                                       truncateAtParticleHorizon              =.false.                                                      , &
       &                                                                                                                       cosmologyParameters_                   =cosmologyParametersSimple_                                   , &
       &                                                                                                                       cosmologyFunctions_                    =cosmologyFunctionsMatterLambda_                              , &
       &                                                                                                                       linearGrowth_                          =linearGrowthCollisionlessMatter_                             , &
       &                                                                                                                       powerSpectrumPrimordialTransferred_    =powerSpectrumPrimordialTransferredSimple_                    , &
       &                                                                                                                       powerSpectrumWindowFunction_           =powerSpectrumWindowFunctionSharpKSpace_                        &
       &                                                                                                                      )
  criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_=criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt  (                                                                                                       &
       &                                                                                                                       linearGrowth_                          =linearGrowthCollisionlessMatter_                             , &
       &                                                                                                                       cosmologyFunctions_                    =cosmologyFunctionsMatterLambda_                              , &
       &                                                                                                                       cosmologicalMassVariance_              =cosmologicalMassVarianceFilteredPower_                       , &
       &                                                                                                                       darkMatterParticle_                    =darkMatterParticleCDM_                                       , &
       &                                                                                                                       tableStore                             =.true.                                                         &
       &                                                                                                                      )
  mergerTreeBranchingProbabilityModifierIdentity_               =mergerTreeBranchingProbabilityModifierIdentity               (                                                                                                       &
       &                                                                                                                      )
  mergerTreeBranchingProbabilityParkinsonColeHelly_(1)          =mergerTreeBranchingProbabilityParkinsonColeHelly             (                                                                                                       &
       &                                                                                                                       G0                                     =1.0d+0                                                       , &
       &                                                                                                                       gamma1                                 =0.0d+0                                                       , &
       &                                                                                                                       gamma2                                 =0.0d+0                                                       , &
       &                                                                                                                       accuracyFirstOrder                     =1.0d-1                                                       , &
       &                                                                                                                       precisionHypergeometric                =1.0d-6                                                       , &
       &                                                                                                                       hypergeometricTabulate                 =.true.                                                       , &
       &                                                                                                                       cdmAssumptions                         =.true.                                                       , &
       &                                                                                                                       tolerateRoundOffErrors                 =.false.                                                      , &
       &                                                                                                                       cosmologicalMassVariance_              =cosmologicalMassVarianceFilteredPower_                       , &
       &                                                                                                                       criticalOverdensity_                   =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_  &
       &                                                                                                                      )
  mergerTreeBranchingProbabilityParkinsonColeHelly_(2)          =mergerTreeBranchingProbabilityParkinsonColeHelly             (                                                                                                       &
       &                                                                                                                       G0                                     =1.0d+0                                                       , &
       &                                                                                                                       gamma1                                 =0.0d+0                                                       , &
       &                                                                                                                       gamma2                                 =0.0d+0                                                       , &
       &                                                                                                                       accuracyFirstOrder                     =1.0d-1                                                       , &
       &                                                                                                                       precisionHypergeometric                =1.0d-6                                                       , &
       &                                                                                                                       hypergeometricTabulate                 =.false.                                                      , &
       &                                                                                                                       cdmAssumptions                         =.true.                                                       , &
       &                                                                                                                       tolerateRoundOffErrors                 =.false.                                                      , &
       &                                                                                                                       cosmologicalMassVariance_              =cosmologicalMassVarianceFilteredPower_                       , &
       &                                                                                                                       criticalOverdensity_                   =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_  &
       &                                                                                                                      )
  mergerTreeBranchingProbabilityParkinsonColeHelly_(3)          =mergerTreeBranchingProbabilityParkinsonColeHelly             (                                                                                                       &
       &                                                                                                                       G0                                     =1.0d+0                                                       , &
       &                                                                                                                       gamma1                                 =0.0d+0                                                       , &
       &                                                                                                                       gamma2                                 =0.0d+0                                                       , &
       &                                                                                                                       accuracyFirstOrder                     =1.0d-1                                                       , &
       &                                                                                                                       precisionHypergeometric                =1.0d-6                                                       , &
       &                                                                                                                       hypergeometricTabulate                 =.false.                                                      , &
       &                                                                                                                       cdmAssumptions                         =.false.                                                      , &
       &                                                                                                                       tolerateRoundOffErrors                 =.false.                                                      , &
       &                                                                                                                       cosmologicalMassVariance_              =cosmologicalMassVarianceFilteredPower_                       , &
       &                                                                                                                       criticalOverdensity_                   =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_  &
       &                                                                                                                      )
  excursionSetBarrierCriticalOverdensity_                       =excursionSetBarrierCriticalOverdensity                       (                                                                                                       &
       &                                                                                                                       criticalOverdensity_                   =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_, &
       &                                                                                                                       cosmologicalMassVariance_              =cosmologicalMassVarianceFilteredPower_                         &
       &                                                                                                                      )
  excursionSetFirstCrossingLinearBarrier_                       =excursionSetFirstCrossingLinearBarrier                       (                                                                                                       &
       &                                                                                                                       fractionalTimeStep                     =fractionalTimeStep                                           , &
       &                                                                                                                       excursionSetBarrier_                   =excursionSetBarrierCriticalOverdensity_                      , &
       &                                                                                                                       cosmologicalMassVariance_              =cosmologicalMassVarianceFilteredPower_                         &
       &                                                                                                                      )
  excursionSetFirstCrossingLinearBarrierBrownianBridge_         =excursionSetFirstCrossingLinearBarrierBrownianBridge         (                                                                                                       &
       &                                                                                                                       varianceConstrained                    =varianceConstrained                                          , &
       &                                                                                                                       criticalOverdensityConstrained         =criticalOverdensityConstrained                               , &
       &                                                                                                                       fractionalTimeStep                     =fractionalTimeStep                                           , &
       &                                                                                                                       cosmologyFunctions_                    =cosmologyFunctionsMatterLambda_                              , &
       &                                                                                                                       excursionSetFirstCrossing_             =excursionSetFirstCrossingLinearBarrier_                      , &
       &                                                                                                                       excursionSetBarrier_                   =excursionSetBarrierCriticalOverdensity_                      , &
       &                                                                                                                       cosmologicalMassVariance_              =cosmologicalMassVarianceFilteredPower_                       , &
       &                                                                                                                       criticalOverdensity_                   =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_, &
       &                                                                                                                       linearGrowth_                          =linearGrowthCollisionlessMatter_                               &
       &                                                                                                                      )
  excursionSetFirstCrossingFarahiMidpoint_                      =excursionSetFirstCrossingFarahiMidpoint                      (                                                                                                       &
       &                                                                                                                       fractionalTimeStep                     =fractionalTimeStep                                           , &
       &                                                                                                                       fileName                               =var_str("auto")                                              , &
       &                                                                                                                       varianceNumberPerUnitProbability       =100                                                          , &
       &                                                                                                                       varianceNumberPerUnit                  = 16                                                          , &
       &                                                                                                                       varianceNumberPerDecade                = 32                                                          , &
       &                                                                                                                       varianceNumberPerDecadeNonCrossing     =  8                                                          , &
       &                                                                                                                       timeNumberPerDecade                    = 16                                                          , &
       &                                                                                                                       varianceIsUnlimited                    =.true.                                                       , &
       &                                                                                                                       cosmologyFunctions_                    =cosmologyFunctionsMatterLambda_                              , &
       &                                                                                                                       excursionSetBarrier_                   =excursionSetBarrierCriticalOverdensity_                      , &
       &                                                                                                                       cosmologicalMassVariance_              =cosmologicalMassVarianceFilteredPower_                         &
       &                                                                                                                      )
  excursionSetFirstCrossingFarahiMidpointBrownianBridge_        =excursionSetFirstCrossingFarahiMidpointBrownianBridge        (                                                                                                       &
       &                                                                                                                       varianceConstrained                    =varianceConstrained                                          , &
       &                                                                                                                       criticalOverdensityConstrained         =  criticalOverdensityConstrained                             , &
       &                                                                                                                       fractionalTimeStep                     =fractionalTimeStep                                           , &
       &                                                                                                                       fileName                               =var_str("auto")                                              , &
       &                                                                                                                       varianceNumberPerUnitProbability       =100                                                          , &
       &                                                                                                                       varianceNumberPerUnit                  = 32                                                          , &
       &                                                                                                                       varianceNumberPerDecade                = 32                                                          , &
       &                                                                                                                       varianceNumberPerDecadeNonCrossing     =  8                                                          , &
       &                                                                                                                       timeNumberPerDecade                    = 64                                                          , &
       &                                                                                                                       varianceIsUnlimited                    =.false.                                                      , &
       &                                                                                                                       cosmologyFunctions_                    =cosmologyFunctionsMatterLambda_                              , &
       &                                                                                                                       excursionSetBarrier_                   =excursionSetBarrierCriticalOverdensity_                      , &
       &                                                                                                                       cosmologicalMassVariance_              =cosmologicalMassVarianceFilteredPower_                       , &
       &                                                                                                                       criticalOverdensity_                   =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_, &
       &                                                                                                                       linearGrowth_                          =linearGrowthCollisionlessMatter_                             , &
       &                                                                                                                       excursionSetFirstCrossing_             =excursionSetFirstCrossingLinearBarrier_                        &
       &                                                                                                                      )
  mergerTreeBranchingProbabilityGnrlzdPrssSchchtr_(1)           =mergerTreeBranchingProbabilityGnrlzdPrssSchchtr              (                                                                                                       &
       &                                                                                                                       deltaStepMaximum                       =1.0d-1                                                       , &
       &                                                                                                                       massMinimum                            =1.0d+0                                                       , &
       &                                                                                                                       smoothAccretion                        =.false.                                                      , &
       &                                                                                                                       distributionFunctionLowerHalfOnly      =.true.                                                       , &
       &                                                                                                                       distributionFunctionNormalize          =.true.                                                       , &
       &                                                                                                                       cosmologyFunctions_                    =cosmologyFunctionsMatterLambda_                              , &
       &                                                                                                                       criticalOverdensity_                   =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_, &
       &                                                                                                                       cosmologicalMassVariance_              =cosmologicalMassVarianceFilteredPower_                       , &
       &                                                                                                                       excursionSetFirstCrossing_             =excursionSetFirstCrossingLinearBarrier_                      , &
       &                                                                                                                       mergerTreeBranchingProbabilityModifier_=mergerTreeBranchingProbabilityModifierIdentity_                &
       &                                                                                                                      )
  mergerTreeBranchingProbabilityGnrlzdPrssSchchtr_(2)           =mergerTreeBranchingProbabilityGnrlzdPrssSchchtr              (                                                                                                       &
       &                                                                                                                       deltaStepMaximum                       =1.0d-1                                                       , &
       &                                                                                                                       massMinimum                            =1.0d-1*massResolution                                        , &
       &                                                                                                                       smoothAccretion                        =.false.                                                      , &
       &                                                                                                                       distributionFunctionLowerHalfOnly      =.true.                                                       , &
       &                                                                                                                       distributionFunctionNormalize          =.true.                                                       , &
       &                                                                                                                       cosmologyFunctions_                    =cosmologyFunctionsMatterLambda_                              , &
       &                                                                                                                       criticalOverdensity_                   =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_, &
       &                                                                                                                       cosmologicalMassVariance_              =cosmologicalMassVarianceFilteredPower_                       , &
       &                                                                                                                       excursionSetFirstCrossing_             =excursionSetFirstCrossingLinearBarrier_                      , &
       &                                                                                                                       mergerTreeBranchingProbabilityModifier_=mergerTreeBranchingProbabilityModifierIdentity_                &
       &                                                                                                                      )
  mergerTreeBranchingProbabilityGnrlzdPrssSchchtr_(3)           =mergerTreeBranchingProbabilityGnrlzdPrssSchchtr              (                                                                                                       &
       &                                                                                                                       deltaStepMaximum                       =1.0d-1                                                       , &
       &                                                                                                                       massMinimum                            =1.0d-1*massResolution                                        , &
       &                                                                                                                       smoothAccretion                        =.false.                                                      , &
       &                                                                                                                       distributionFunctionLowerHalfOnly      =.true.                                                       , &
       &                                                                                                                       distributionFunctionNormalize          =.true.                                                       , &
       &                                                                                                                       cosmologyFunctions_                    =cosmologyFunctionsMatterLambda_                              , &
       &                                                                                                                       criticalOverdensity_                   =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_, &
       &                                                                                                                       cosmologicalMassVariance_              =cosmologicalMassVarianceFilteredPower_                       , &
       &                                                                                                                       excursionSetFirstCrossing_             =excursionSetFirstCrossingFarahiMidpoint_                     , &
       &                                                                                                                       mergerTreeBranchingProbabilityModifier_=mergerTreeBranchingProbabilityModifierIdentity_                &
       &                                                                                                                      )
  mergerTreeBranchingProbabilityGnrlzdPrssSchchtr_(4)           =mergerTreeBranchingProbabilityGnrlzdPrssSchchtr              (                                                                                                       &
       &                                                                                                                       deltaStepMaximum                       =1.0d-1                                                       , &
       &                                                                                                                       massMinimum                            =1.0d-1*massResolution                                        , &
       &                                                                                                                       smoothAccretion                        =.true.                                                       , &
       &                                                                                                                       distributionFunctionLowerHalfOnly      =.true.                                                       , &
       &                                                                                                                       distributionFunctionNormalize          =.true.                                                       , &
       &                                                                                                                       cosmologyFunctions_                    =cosmologyFunctionsMatterLambda_                              , &
       &                                                                                                                       criticalOverdensity_                   =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_, &
       &                                                                                                                       cosmologicalMassVariance_              =cosmologicalMassVarianceFilteredPower_                       , &
       &                                                                                                                       excursionSetFirstCrossing_             =excursionSetFirstCrossingLinearBarrier_                      , &
       &                                                                                                                       mergerTreeBranchingProbabilityModifier_=mergerTreeBranchingProbabilityModifierIdentity_                &
       &                                                                                                                      )
  mergerTreeBranchingProbabilityGnrlzdPrssSchchtr_(5)           =mergerTreeBranchingProbabilityGnrlzdPrssSchchtr              (                                                                                                       &
       &                                                                                                                       deltaStepMaximum                       =1.0d-1                                                       , &
       &                                                                                                                       massMinimum                            =1.0d-1*massResolution                                        , &
       &                                                                                                                       smoothAccretion                        =.true.                                                       , &
       &                                                                                                                       distributionFunctionLowerHalfOnly      =.true.                                                       , &
       &                                                                                                                       distributionFunctionNormalize          =.true.                                                       , &
       &                                                                                                                       cosmologyFunctions_                    =cosmologyFunctionsMatterLambda_                              , &
       &                                                                                                                       criticalOverdensity_                   =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_, &
       &                                                                                                                       cosmologicalMassVariance_              =cosmologicalMassVarianceFilteredPower_                       , &
       &                                                                                                                       excursionSetFirstCrossing_             =excursionSetFirstCrossingFarahiMidpoint_                     , &
       &                                                                                                                       mergerTreeBranchingProbabilityModifier_=mergerTreeBranchingProbabilityModifierIdentity_                &
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
     call Assert('Unconstrained case',errorMaximumUnconstrained,3.3d-2,compareLessThan)
     call Assert('Constrained case'  ,errorMaximumConstrained  ,3.3d-2,compareLessThan)
     call Unit_Tests_End_Group()
     call Unit_Tests_End_Group()
  end do
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_Merger_Tree_Branching
