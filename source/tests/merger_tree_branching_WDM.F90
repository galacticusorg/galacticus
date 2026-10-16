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

program Tests_Merger_Tree_Branching_WDM
  !!{RST
  Tests of merger tree branching probabilities for a truncated power spectrum.

  The companion tests use power-law and :term:`CDM` power spectra, in both of which
  :math:`\alpha \equiv -\mathrm{d}\log\sigma/\mathrm{d}\log M` is comfortably positive everywhere. A truncated spectrum---here that
  of a thermal warm dark matter particle---is the case the branching probability finds hardest: below the cut-off :math:`\sigma(M)`
  flattens, :math:`\alpha` falls towards zero, and the assumptions which the :cite:t:`parkinson_generating_2008` algorithm may make
  for :term:`CDM` cease to hold. That is the regime in which the ``cdmAssumptions`` and ``tolerateRoundOffErrors`` options exist to
  be used, and it is exercised here.
  !!}
  use :: Cosmological_Density_Field          , only : cosmologicalMassVarianceFilteredPower           , criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt
  use :: Cosmology_Functions                 , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters                , only : cosmologyParametersSimple
  use :: Dark_Matter_Particles               , only : darkMatterParticleCDM                          , darkMatterParticleWDMThermal
  use :: Display                             , only : displayVerbositySet                            , verbosityLevelWorking
  use :: Error                               , only : Error_Handler_Register
  use :: Events_Hooks                        , only : eventsHooksInitialize
  use :: Galacticus_Nodes                    , only : treeNode
  use :: Linear_Growth                       , only : linearGrowthCollisionlessMatter
  use :: Merger_Tree_Branching               , only : mergerTreeBranchingProbabilityParkinsonColeHelly, mergerTreeBranchingBoundLower                              , mergerTreeBranchingBoundUpper
  use :: Power_Spectra_Primordial            , only : powerSpectrumPrimordialPowerLaw
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredSimple
  use :: Power_Spectrum_Window_Functions     , only : powerSpectrumWindowFunctionTopHat
  use :: Transfer_Functions                  , only : transferFunctionEisensteinHu1999               , transferFunctionBode2001                                   , scaleCutOffModelBode2001
  use :: Unit_Tests                          , only : Assert                                         , Unit_Tests_Begin_Group                                     , Unit_Tests_End_Group         , Unit_Tests_Finish, &
       &                                              compareLessThanOrEqual                         , compareGreaterThan
  implicit none
  type            (cosmologyParametersSimple                                   )                :: cosmologyParametersSimple_
  type            (cosmologyFunctionsMatterLambda                              )                :: cosmologyFunctionsMatterLambda_
  ! A second mass variance, built on the untruncated CDM transfer function, against which the suppression of α by the truncation
  ! is measured. Comparing against CDM rather than against an absolute threshold keeps the check meaningful whatever the
  ! particle mass and filter happen to give.
  type            (cosmologicalMassVarianceFilteredPower                       )                :: cosmologicalMassVarianceFilteredPower_   , cosmologicalMassVarianceCDM_
  type            (powerSpectrumWindowFunctionTopHat                           )                :: powerSpectrumWindowFunctionTopHat_
  type            (powerSpectrumPrimordialPowerLaw                             )                :: powerSpectrumPrimordialPowerLaw_
  type            (transferFunctionEisensteinHu1999                            )                :: transferFunctionEisensteinHu1999_
  type            (transferFunctionBode2001                                    )                :: transferFunctionBode2001_
  type            (powerSpectrumPrimordialTransferredSimple                    )                :: powerSpectrumPrimordialTransferredSimple_, powerSpectrumPrimordialTransferredCDM_
  type            (linearGrowthCollisionlessMatter                             )                :: linearGrowthCollisionlessMatter_
  type            (criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt)                :: criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_
  type            (darkMatterParticleCDM                                       )                :: darkMatterParticleCDM_
  type            (darkMatterParticleWDMThermal                                )                :: darkMatterParticleWDMThermal_
  type            (treeNode                                                    )                :: node
  type            (mergerTreeBranchingProbabilityParkinsonColeHelly            )                :: branchingWDM_                        , branchingWDMTabulated_
  ! Mass of the thermal warm dark matter particle, in keV. A light particle is chosen so that the cut-off falls well within the
  ! range of masses swept below, placing the branching probability firmly in the truncated regime rather than merely at its edge.
  double precision                                                              , parameter     :: massParticleWDM      =1.0d0
  integer                                                                       , parameter     :: countMass            =15
  double precision                                                              , parameter     :: massMinimum          =1.0d7         , massMaximum           =1.0d14, &
       &                                                                                           massResolutionFraction=1.0d-2
  double precision                                                              , dimension(countMass) :: mass                          , probabilityMass       , &
       &                                                                                                 boundLowerMass                , boundUpperMass        , &
       &                                                                                                 alphaMass                     , sigmaMass             , &
       &                                                                                                 fractionSubresolutionMass     , alphaMassCDM          , &
       &                                                                                                 sigmaMassCDM
  double precision                                                                               :: time                                , deltaCritical_        , &
       &                                                                                            massResolution_
  integer                                                                                        :: i

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelWorking)
  ! Register error handlers.
  call Error_Handler_Register()
  ! Initialize event hooks.
  call eventsHooksInitialize()
  ! Build a ΛCDM cosmology whose power spectrum is truncated by the free streaming of a thermal warm dark matter particle.
  darkMatterParticleCDM_                                       =darkMatterParticleCDM                                        (                                                                                          &
       &                                                                                                                     )
  cosmologyParametersSimple_                                   =cosmologyParametersSimple                                    (                                                                                          &
       &                                                                                                                      OmegaMatter                             = 0.3000d0                                      , &
       &                                                                                                                      OmegaBaryon                             = 0.0450d0                                      , &
       &                                                                                                                      OmegaDarkEnergy                         = 0.7000d0                                      , &
       &                                                                                                                      temperatureCMB                          = 2.7250d0                                      , &
       &                                                                                                                      HubbleConstant                          =70.0000d0                                        &
       &                                                                                                                     )
  darkMatterParticleWDMThermal_                                =darkMatterParticleWDMThermal                                 (                                                                                          &
       &                                                                                                                      mass                                    =massParticleWDM                                , &
       &                                                                                                                      degreesOfFreedomEffective               =1.5d0                                          , &
       &                                                                                                                      cosmologyParameters_                    =cosmologyParametersSimple_                       &
       &                                                                                                                     )
  cosmologyFunctionsMatterLambda_                              =cosmologyFunctionsMatterLambda                               (                                                                                          &
       &                                                                                                                      cosmologyParameters_                    =cosmologyParametersSimple_                       &
       &                                                                                                                     )
  linearGrowthCollisionlessMatter_                             =linearGrowthCollisionlessMatter                              (                                                                                          &
       &                                                                                                                      cosmologyParameters_                    =cosmologyParametersSimple_                     , &
       &                                                                                                                      cosmologyFunctions_                     =cosmologyFunctionsMatterLambda_                  &
       &                                                                                                                     )
  powerSpectrumPrimordialPowerLaw_                             =powerSpectrumPrimordialPowerLaw                              (                                                                                          &
       &                                                                                                                      index_                                  =+0.96d0                                        , &
       &                                                                                                                      running                                 =+0.00d0                                        , &
       &                                                                                                                      runningRunning                          =+0.00d0                                        , &
       &                                                                                                                      wavenumberReference                     =+1.00d0                                        , &
       &                                                                                                                      runningSmallScalesOnly                  =.false.                                          &
       &                                                                                                                     )
  transferFunctionEisensteinHu1999_                            =transferFunctionEisensteinHu1999                             (                                                                                          &
       &                                                                                                                      neutrinoNumberEffective                 =3.046d0                                        , &
       &                                                                                                                      neutrinoMassSummed                      =0.000d0                                        , &
       &                                                                                                                      darkMatterParticle_                     =darkMatterParticleCDM_                         , &
       &                                                                                                                      cosmologyParameters_                    =cosmologyParametersSimple_                     , &
       &                                                                                                                      cosmologyFunctions_                     =cosmologyFunctionsMatterLambda_                  &
       &                                                                                                                     )
  transferFunctionBode2001_                                    =transferFunctionBode2001                                     (                                                                                          &
       &                                                                                                                      transferFunctionCDM                     =transferFunctionEisensteinHu1999_              , &
       &                                                                                                                      scaleCutOffModel                        =scaleCutOffModelBode2001                       , &
       &                                                                                                                      epsilon                                 =0.359d0                                        , &
       &                                                                                                                      eta                                     =3.810d0                                        , &
       &                                                                                                                      nu                                      =1.100d0                                        , &
       &                                                                                                                      time                                    =13.8d0                                         , &
       &                                                                                                                      cosmologyParameters_                    =cosmologyParametersSimple_                     , &
       &                                                                                                                      darkMatterParticle_                     =darkMatterParticleWDMThermal_                  , &
       &                                                                                                                      cosmologyFunctions_                     =cosmologyFunctionsMatterLambda_                  &
       &                                                                                                                     )
  powerSpectrumPrimordialTransferredSimple_                    =powerSpectrumPrimordialTransferredSimple                     (                                                                                          &
       &                                                                                                                      powerSpectrumPrimordial_                =powerSpectrumPrimordialPowerLaw_               , &
       &                                                                                                                      transferFunction_                       =transferFunctionBode2001_                      , &
       &                                                                                                                      linearGrowth_                           =linearGrowthCollisionlessMatter_                 &
       &                                                                                                                     )
  powerSpectrumWindowFunctionTopHat_                           =powerSpectrumWindowFunctionTopHat                            (                                                                                          &
       &                                                                                                                      cosmologyParameters_                    =cosmologyParametersSimple_                       &
       &                                                                                                                     )
  cosmologicalMassVarianceFilteredPower_                       =cosmologicalMassVarianceFilteredPower                        (                                                                                          &
       &                                                                                                                      sigma8                                  =0.8d+0                                         , &
       &                                                                                                                      tolerance                               =1.0d-4                                         , &
       &                                                                                                                      toleranceTopHat                         =1.0d-4                                         , &
       &                                                                                                                      rootVarianceLogarithmicGradientTolerance=1.0d-6                                         , &
       &                                                                                                                      integrationFailureIsFatal               =.true.                                         , &
       &                                                                                                                      storeTabulations                        =.true.                                         , &
       &                                                                                                                      nonMonotonicIsFatal                     =.true.                                         , &
       &                                                                                                                      monotonicInterpolation                  =.false.                                        , &
       &                                                                                                                      truncateAtParticleHorizon               =.false.                                        , &
       &                                                                                                                      cosmologyParameters_                    =cosmologyParametersSimple_                     , &
       &                                                                                                                      cosmologyFunctions_                     =cosmologyFunctionsMatterLambda_                , &
       &                                                                                                                      linearGrowth_                           =linearGrowthCollisionlessMatter_               , &
       &                                                                                                                      powerSpectrumPrimordialTransferred_     =powerSpectrumPrimordialTransferredSimple_      , &
       &                                                                                                                      powerSpectrumWindowFunction_            =powerSpectrumWindowFunctionTopHat_               &
       &                                                                                                                     )
  powerSpectrumPrimordialTransferredCDM_                       =powerSpectrumPrimordialTransferredSimple                     (                                                                                          &
       &                                                                                                                      powerSpectrumPrimordial_                =powerSpectrumPrimordialPowerLaw_               , &
       &                                                                                                                      transferFunction_                       =transferFunctionEisensteinHu1999_              , &
       &                                                                                                                      linearGrowth_                           =linearGrowthCollisionlessMatter_                 &
       &                                                                                                                     )
  cosmologicalMassVarianceCDM_                                 =cosmologicalMassVarianceFilteredPower                        (                                                                                          &
       &                                                                                                                      sigma8                                  =0.8d+0                                         , &
       &                                                                                                                      tolerance                               =1.0d-4                                         , &
       &                                                                                                                      toleranceTopHat                         =1.0d-4                                         , &
       &                                                                                                                      rootVarianceLogarithmicGradientTolerance=1.0d-6                                         , &
       &                                                                                                                      integrationFailureIsFatal               =.true.                                         , &
       &                                                                                                                      storeTabulations                        =.true.                                         , &
       &                                                                                                                      nonMonotonicIsFatal                     =.true.                                         , &
       &                                                                                                                      monotonicInterpolation                  =.false.                                        , &
       &                                                                                                                      truncateAtParticleHorizon               =.false.                                        , &
       &                                                                                                                      cosmologyParameters_                    =cosmologyParametersSimple_                     , &
       &                                                                                                                      cosmologyFunctions_                     =cosmologyFunctionsMatterLambda_                , &
       &                                                                                                                      linearGrowth_                           =linearGrowthCollisionlessMatter_               , &
       &                                                                                                                      powerSpectrumPrimordialTransferred_     =powerSpectrumPrimordialTransferredCDM_         , &
       &                                                                                                                      powerSpectrumWindowFunction_            =powerSpectrumWindowFunctionTopHat_               &
       &                                                                                                                     )
  criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_=criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt (                                                                                          &
       &                                                                                                                      linearGrowth_                           =linearGrowthCollisionlessMatter_               , &
       &                                                                                                                      cosmologyFunctions_                     =cosmologyFunctionsMatterLambda_                , &
       &                                                                                                                      cosmologicalMassVariance_               =cosmologicalMassVarianceFilteredPower_         , &
       &                                                                                                                      darkMatterParticle_                     =darkMatterParticleCDM_                         , &
       &                                                                                                                      tableStore                              =.true.                                           &
       &                                                                                                                     )
  ! The CDM assumptions must not be made for a truncated spectrum - that is precisely what the option exists to disable - and
  ! round-off errors in the branching integrals are tolerated, as the documentation of that option anticipates for spectra with a
  ! cut-off.
  branchingWDM_         =mergerTreeBranchingProbabilityParkinsonColeHelly(                                                                                                    &
       &                                                                  G0                       =0.57d0                                                                  , &
       &                                                                  gamma1                   =0.38d0                                                                  , &
       &                                                                  gamma2                   =-0.01d0                                                                 , &
       &                                                                  accuracyFirstOrder       =1.0d-1                                                                  , &
       &                                                                  precisionHypergeometric  =1.0d-6                                                                  , &
       &                                                                  hypergeometricTabulate   =.false.                                                                 , &
       &                                                                  cdmAssumptions           =.false.                                                                 , &
       &                                                                  tolerateRoundOffErrors   =.true.                                                                  , &
       &                                                                  cosmologicalMassVariance_=cosmologicalMassVarianceFilteredPower_                                   , &
       &                                                                  criticalOverdensity_     =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_             &
       &                                                                 )
  branchingWDMTabulated_=mergerTreeBranchingProbabilityParkinsonColeHelly(                                                                                                    &
       &                                                                  G0                       =0.57d0                                                                  , &
       &                                                                  gamma1                   =0.38d0                                                                  , &
       &                                                                  gamma2                   =-0.01d0                                                                 , &
       &                                                                  accuracyFirstOrder       =1.0d-1                                                                  , &
       &                                                                  precisionHypergeometric  =1.0d-6                                                                  , &
       &                                                                  hypergeometricTabulate   =.true.                                                                  , &
       &                                                                  cdmAssumptions           =.false.                                                                 , &
       &                                                                  tolerateRoundOffErrors   =.true.                                                                  , &
       &                                                                  cosmologicalMassVariance_=cosmologicalMassVarianceFilteredPower_                                   , &
       &                                                                  criticalOverdensity_     =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_             &
       &                                                                 )
  do i=1,countMass
     mass(i)=massMinimum*(massMaximum/massMinimum)**(dble(i-1)/dble(countMass-1))
  end do
  time=cosmologyFunctionsMatterLambda_%cosmicTime(1.0d0)

  call Unit_Tests_Begin_Group("Merger tree branching (truncated power spectrum)")

  ! Establish that the masses swept really do straddle the cut-off: below it σ(M) flattens and α falls away, which is the
  ! circumstance the tests which follow are meant to place the branching probability in. Without this the remaining tests could
  ! pass while exercising nothing more than an ordinary CDM spectrum.
  do i=1,countMass
     call cosmologicalMassVarianceFilteredPower_%rootVarianceAndLogarithmicGradient(mass(i),time,sigmaMass   (i),alphaMass   (i))
     call cosmologicalMassVarianceCDM_          %rootVarianceAndLogarithmicGradient(mass(i),time,sigmaMassCDM(i),alphaMassCDM(i))
  end do
  call Unit_Tests_Begin_Group("Power spectrum is truncated")
  call Assert('α is strongly suppressed relative to CDM at the lightest mass',abs(alphaMass(1)       ),0.3d0*abs(alphaMassCDM(1       )),compareLessThanOrEqual)
  call Assert('α is unsuppressed at the heaviest mass'                       ,abs(alphaMass(countMass)),0.9d0*abs(alphaMassCDM(countMass)),compareGreaterThan   )
  call Unit_Tests_End_Group  (                                                       )

  ! The bounds must still bound, and the branching probability and subresolution accretion rate must remain finite and positive,
  ! through and below the cut-off.
  do i=1,countMass
     massResolution_              =mass(i)*massResolutionFraction
     deltaCritical_               =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_%value(time=time,mass=mass(i),node=node)
     probabilityMass          (i) =branchingWDM_%probability          (mass(i),deltaCritical_,time,massResolution_                              ,node)
     boundLowerMass           (i) =branchingWDM_%probabilityBound     (mass(i),deltaCritical_,time,massResolution_,mergerTreeBranchingBoundLower,node)
     boundUpperMass           (i) =branchingWDM_%probabilityBound     (mass(i),deltaCritical_,time,massResolution_,mergerTreeBranchingBoundUpper,node)
     fractionSubresolutionMass(i) =branchingWDM_%fractionSubresolution(mass(i),deltaCritical_,time,massResolution_                              ,node)
  end do
  call Unit_Tests_Begin_Group("Bounds bracket the probability")
  call Assert('lower bound ≤ probability',boundLowerMass ,probabilityMass,compareLessThanOrEqual)
  call Assert('probability ≤ upper bound',probabilityMass,boundUpperMass ,compareLessThanOrEqual)
  call Assert('probability is positive'  ,probabilityMass,spread(0.0d0,1,countMass),compareGreaterThan)
  call Unit_Tests_End_Group  (                          )

  ! Tabulating the hypergeometric functions must not change the answer. Below the cut-off the ratio σ(M_res)/σ(M) approaches unity,
  ! which drives the argument of those functions towards the singular point at which their evaluation is most delicate, so this
  ! also checks that the tabulation extends properly into that limit.
  do i=1,countMass
     massResolution_ =mass(i)*massResolutionFraction
     deltaCritical_  =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_%value(time=time,mass=mass(i),node=node)
     boundUpperMass(i)=branchingWDMTabulated_%fractionSubresolution(mass(i),deltaCritical_,time,massResolution_,node)
  end do
  call Unit_Tests_Begin_Group("Tabulation is accurate through the cut-off")
  call Assert('subresolution accretion rate, tabulated against computed',boundUpperMass,fractionSubresolutionMass,relTol=5.0d-3)
  call Unit_Tests_End_Group  (                                                       )

  call Unit_Tests_End_Group  (                                                       )
  call Unit_Tests_Finish     (                                                       )
end program Tests_Merger_Tree_Branching_WDM
