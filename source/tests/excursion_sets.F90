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

program Tests_Excursion_Sets
  !!{RST
  Tests of excursion set calculations.
  !!}
  use :: Cosmological_Density_Field          , only : cosmologicalMassVarianceFilteredPower          , criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt
  use :: Cosmology_Functions                 , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters                , only : cosmologyParametersSimple
  use :: Dark_Matter_Particles               , only : darkMatterParticleCDM
  use :: Display                             , only : displayVerbositySet                            , verbosityLevelStandard
  use :: Events_Hooks                        , only : eventsHooksInitialize
  use :: Functions_Global_Utilities          , only : Functions_Global_Set
  use :: Excursion_Sets_Barriers             , only : excursionSetBarrierCriticalOverdensity
  use :: File_Utilities                      , only : File_Exists                                    , File_Remove
  use :: Excursion_Sets_First_Crossings      , only : excursionSetFirstCrossingClass                 , excursionSetFirstCrossingFarahi                             , excursionSetFirstCrossingFarahiMidpoint, excursionSetFirstCrossingLinearBarrier, &
          &                                           excursionSetFirstCrossingZhangHui              , excursionSetFirstCrossingZhangHuiHighOrder
  use :: Galacticus_Nodes                    , only : treeNode
  use :: Input_Paths                         , only : inputPath                                      , pathTypeDataDynamic
  use :: ISO_Varying_String                  , only : var_str                                        , varying_string                                              , operator(//)
  use :: Linear_Growth                       , only : linearGrowthCollisionlessMatter
  use :: Merger_Tree_Branching               , only : mergerTreeBranchingProbabilityGeneralizedPressSchechter, mergerTreeBranchingProbabilityParkinsonColeHelly
  use :: Merger_Tree_Branching_Modifiers     , only : mergerTreeBranchingProbabilityModifierIdentity
  use :: Power_Spectra_Primordial            , only : powerSpectrumPrimordialPowerLaw
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredSimple
  use :: Power_Spectrum_Normalizers          , only : powerSpectrumNormalizerSigma8
  use :: Power_Spectrum_Window_Functions     , only : powerSpectrumWindowFunctionSharpKSpace
  use :: Transfer_Functions                  , only : transferFunctionIdentity
  use :: Unit_Tests                          , only : Assert                                         , Unit_Tests_Begin_Group                                      , Unit_Tests_End_Group                   , Unit_Tests_Finish
  implicit none
  type            (cosmologyParametersSimple                                   )               :: cosmologyParametersSimple_
  type            (cosmologyFunctionsMatterLambda                              )               :: cosmologyFunctionsMatterLambda_
  type            (cosmologicalMassVarianceFilteredPower                       )               :: cosmologicalMassVarianceFilteredPower_
  type            (powerSpectrumWindowFunctionSharpKSpace                      )               :: powerSpectrumWindowFunctionSharpKSpace_
  type            (powerSpectrumPrimordialPowerLaw                             )               :: powerSpectrumPrimordialPowerLaw_
  type            (transferFunctionIdentity                                    )               :: transferFunctionIdentity_
  type            (powerSpectrumPrimordialTransferredSimple                    )               :: powerSpectrumPrimordialTransferredSimple_
  type            (powerSpectrumNormalizerSigma8                               )               :: powerSpectrumNormalizerSigma8_
  type            (linearGrowthCollisionlessMatter                             )               :: linearGrowthCollisionlessMatter_
  type            (excursionSetFirstCrossingLinearBarrier                      )               :: excursionSetFirstCrossingLinearBarrier_
  type            (excursionSetFirstCrossingFarahiMidpoint                     ), target       :: excursionSetFirstCrossingFarahiMidpoint_
  type            (excursionSetFirstCrossingFarahi                             ), target       :: excursionSetFirstCrossingFarahi_
  type            (excursionSetFirstCrossingFarahiMidpoint                     )               :: excursionSetFirstCrossingFarahiMidpointWrite_                                 , excursionSetFirstCrossingFarahiMidpointRead_
  type            (excursionSetFirstCrossingZhangHui                           ), target       :: excursionSetFirstCrossingZhangHui_
  type            (excursionSetFirstCrossingZhangHuiHighOrder                  ), target       :: excursionSetFirstCrossingZhangHuiHighOrder_
  class           (excursionSetFirstCrossingClass                              ), pointer      :: excursionSetFirstCrossing_
  type            (excursionSetBarrierCriticalOverdensity                      )               :: excursionSetBarrierCriticalOverdensity_
  type            (darkMatterParticleCDM                                       )               :: darkMatterParticleCDM_
  type            (treeNode                                                    )               :: node
  type            (criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt)               :: criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_
  double precision                                                                             :: time
  double precision                                                              , dimension(2) :: variance                                                      =[0.50d0,1.00d0]
  double precision                                                              , dimension(2) :: varianceProgenitor                                            =[1.00d0,1.01d0]
  double precision                                                              , dimension(2) :: probabilityTarget                                                             , probability, &
       &                                                                                          rateTarget                                                                    , rate       , &
       &                                                                                          probabilityRead                                                               , rateRead   , &
       &                                                                                          probabilityRates
  double precision                                                              , parameter    :: varianceLarge                                                 =10.0d0
  integer                                                                                      :: i                                                                             , j
  character       (len=32                                                      )               :: label
  type            (varying_string                                              )               :: fileNameRoundTrip
  logical                                                                                      :: testRates

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  call eventsHooksInitialize()
  call Functions_Global_Set ()
  ! Build all objects needed for these tests.
  darkMatterParticleCDM_                                        =darkMatterParticleCDM                                        (                                                                                                         &
       &                                                                                                                      )
  cosmologyParametersSimple_                                    =cosmologyParametersSimple                                    (                                                                                                         &
       &                                                                                                                       OmegaMatter                             = 1.00d0                                                       , &
       &                                                                                                                       OmegaBaryon                             = 0.00d0                                                       , &
       &                                                                                                                       OmegaDarkEnergy                         = 0.00d0                                                       , &
       &                                                                                                                       temperatureCMB                          = 2.78d0                                                       , &
       &                                                                                                                       HubbleConstant                          =70.00d0                                                         &
       &                                                                                                                      )
  cosmologyFunctionsMatterLambda_                               =cosmologyFunctionsMatterLambda                               (                                                                                                         &
       &                                                                                                                       cosmologyParameters_                    =cosmologyParametersSimple_                                      &
       &                                                                                                                      )
  linearGrowthCollisionlessMatter_                              =linearGrowthCollisionlessMatter                              (                                                                                                         &
       &                                                                                                                       cosmologyParameters_                    =cosmologyParametersSimple_                                    , &
       &                                                                                                                       cosmologyFunctions_                     =cosmologyFunctionsMatterLambda_                                 &
       &                                                                                                                      )
  powerSpectrumPrimordialPowerLaw_                              =powerSpectrumPrimordialPowerLaw                              (                                                                                                         &
       &                                                                                                                       index_                                  =-1.0d0                                                        , &
       &                                                                                                                       running                                 =+0.0d0                                                        , &
       &                                                                                                                       runningRunning                          =+0.0d0                                                        , &
       &                                                                                                                       wavenumberReference                     =+1.0d0                                                        , &
       &                                                                                                                       runningSmallScalesOnly                  =.false.                                                         &
       &                                                                                                                      )
  transferFunctionIdentity_                                     =transferFunctionIdentity                                     (                                                                                                         &
       &                                                                                                                       cosmologyParameters_                    =cosmologyParametersSimple_                                    , &
       &                                                                                                                       time                                    =13.8d0                                                          &
       &                                                                                                                      )
  powerSpectrumNormalizerSigma8_                                =powerSpectrumNormalizerSigma8(                                                                                                                                         &
       &                                                                                                                       sigma8                                  =1.0d+0                                                        , &
       &                                                                                                                       tolerance                               =1.0d-4                                                        , &
       &                                                                                                                       integrationFailureIsFatal               =.true.                                                        , &
       &                                                                                                                       cosmologyParameters_                    =cosmologyParametersSimple_                                    , &
       &                                                                                                                       cosmologyFunctions_                     =cosmologyFunctionsMatterLambda_                               , &
       &                                                                                                                       transferFunction_                       =transferFunctionIdentity_                                       &
       &                                                                                                                      )
  powerSpectrumPrimordialTransferredSimple_                     =powerSpectrumPrimordialTransferredSimple                     (                                                                                                         &
       &                                                                                                                       powerSpectrumPrimordial_                =powerSpectrumPrimordialPowerLaw_                              , &
       &                                                                                                                       transferFunction_                       =transferFunctionIdentity_                                     , &
       &                                                                                                                       linearGrowth_                           =linearGrowthCollisionlessMatter_                              , &
       &                                                                                                                       powerSpectrumNormalizer_                =powerSpectrumNormalizerSigma8_                                  &
       &                                                                                                                      )
  powerSpectrumWindowFunctionSharpKSpace_                       =powerSpectrumWindowFunctionSharpKSpace                       (                                                                                                         &
       &                                                                                                                       cosmologyParameters_                    =cosmologyParametersSimple_                                    , &
       &                                                                                                                       normalization                           =0.0d0                                                           &
       &                                                                                                                      )

  cosmologicalMassVarianceFilteredPower_                        =cosmologicalMassVarianceFilteredPower                        (                                                                                                         &
       &                                                                                                                       tolerance                               =1.0d-4                                                        , &
       &                                                                                                                       toleranceTopHat                         =1.0d-4                                                        , &
       &                                                                                                                       rootVarianceLogarithmicGradientTolerance=1.0d-9                                                        , &
       &                                                                                                                       integrationFailureIsFatal               =.true.                                                        , &
       &                                                                                                                       storeTabulations                        =.true.                                                        , &
       &                                                                                                                       nonMonotonicIsFatal                     =.true.                                                        , &
       &                                                                                                                       monotonicInterpolation                  =.false.                                                       , &
       &                                                                                                                       truncateAtParticleHorizon               =.false.                                                       , &
       &                                                                                                                       cosmologyParameters_                    =cosmologyParametersSimple_                                    , &
       &                                                                                                                       cosmologyFunctions_                     =cosmologyFunctionsMatterLambda_                               , &
       &                                                                                                                       linearGrowth_                           =linearGrowthCollisionlessMatter_                              , &
       &                                                                                                                       powerSpectrumPrimordialTransferred_     =powerSpectrumPrimordialTransferredSimple_                     , &
       &                                                                                                                       powerSpectrumWindowFunction_            =powerSpectrumWindowFunctionSharpKSpace_                         &
       &                                                                                                                      )
  criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_=criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt(                                                                                                           &
       &                                                                                                                       linearGrowth_                           =linearGrowthCollisionlessMatter_                              , &
       &                                                                                                                       cosmologyFunctions_                     =cosmologyFunctionsMatterLambda_                               , &
       &                                                                                                                       cosmologicalMassVariance_               =cosmologicalMassVarianceFilteredPower_                        , &
       &                                                                                                                       darkMatterParticle_                     =darkMatterParticleCDM_                                        , &
       &                                                                                                                       tableStore                              =.true.                                                          &
       &                                                                                                                      )
  excursionSetBarrierCriticalOverdensity_                       =excursionSetBarrierCriticalOverdensity                       (                                                                                                         &
       &                                                                                                                       criticalOverdensity_                    =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_,  &
       &                                                                                                                       cosmologicalMassVariance_               =cosmologicalMassVarianceFilteredPower_                          &
       &                                                                                                                      )
  excursionSetFirstCrossingLinearBarrier_                       =excursionSetFirstCrossingLinearBarrier                       (                                                                                                         &
       &                                                                                                                       fractionalTimeStep                      =0.01d0                                                        , &
       &                                                                                                                       excursionSetBarrier_                    =excursionSetBarrierCriticalOverdensity_                       , &
       &                                                                                                                       cosmologicalMassVariance_               =cosmologicalMassVarianceFilteredPower_                          &
       &                                                                                                                      )
  excursionSetFirstCrossingFarahiMidpoint_                      =excursionSetFirstCrossingFarahiMidpoint                      (                                                                                                         &
       &                                                                                                                       fractionalTimeStep                      =0.01d0                                                        , &
       &                                                                                                                       fileName                                =var_str("auto")                                               , &
       &                                                                                                                       varianceNumberPerUnitProbability        =1000                                                          , &
       &                                                                                                                       varianceNumberPerUnit                   =  32                                                          , &
       &                                                                                                                       varianceNumberPerDecade                 =  64                                                          , &
       &                                                                                                                       varianceNumberPerDecadeNonCrossing      =   8                                                         , &
       &                                                                                                                       varianceIsUnlimited                     =.true.                                                        , &
       &                                                                                                                       timeNumberPerDecade                     =  10                                                          , &
       &                                                                                                                       cosmologyFunctions_                     =cosmologyFunctionsMatterLambda_                               , &
       &                                                                                                                       excursionSetBarrier_                    =excursionSetBarrierCriticalOverdensity_                       , &
       &                                                                                                                       cosmologicalMassVariance_               =cosmologicalMassVarianceFilteredPower_                          &
       &                                                                                                                      )
  excursionSetFirstCrossingFarahi_                              =excursionSetFirstCrossingFarahi                              (                                                                                                         &
       &                                                                                                                       fractionalTimeStep                      =0.01d0                                                        , &
       &                                                                                                                       fileName                                =var_str("auto")                                               , &
       &                                                                                                                       varianceNumberPerUnitProbability        =1000                                                          , &
       &                                                                                                                       varianceNumberPerUnit                   =  32                                                          , &
       &                                                                                                                       varianceNumberPerDecade                 =  64                                                          , &
       &                                                                                                                       varianceNumberPerDecadeNonCrossing      =   8                                                          , &
       &                                                                                                                       timeNumberPerDecade                     =  10                                                          , &
       &                                                                                                                       varianceIsUnlimited                     =.true.                                                        , &
       &                                                                                                                       cosmologyFunctions_                     =cosmologyFunctionsMatterLambda_                               , &
       &                                                                                                                       excursionSetBarrier_                    =excursionSetBarrierCriticalOverdensity_                       , &
       &                                                                                                                       cosmologicalMassVariance_               =cosmologicalMassVarianceFilteredPower_                          &
       &                                                                                                                      )
  excursionSetFirstCrossingZhangHui_                            =excursionSetFirstCrossingZhangHui                            (                                                                                                         &
       &                                                                                                                       excursionSetBarrier_                    =excursionSetBarrierCriticalOverdensity_                         &
       &                                                                                                                      )
  excursionSetFirstCrossingZhangHuiHighOrder_                   =excursionSetFirstCrossingZhangHuiHighOrder                   (                                                                                                         &
       &                                                                                                                       excursionSetBarrier_                    =excursionSetBarrierCriticalOverdensity_                         &
       &                                                                                                                      )
  ! Evaluate at a large variance initially to ensure the full range is tabulated.
  time          =cosmologyFunctionsMatterLambda_         %cosmicTime (                               1.0d0     )
  probability(1)=excursionSetFirstCrossingFarahiMidpoint_%probability(varianceLarge                 ,time ,node)
  probability(1)=excursionSetFirstCrossingFarahi_        %probability(varianceLarge                 ,time ,node)
  rate       (1)=excursionSetFirstCrossingFarahiMidpoint_%rate       (variance     (1),varianceLarge,time ,node)
  rate       (1)=excursionSetFirstCrossingFarahi_        %rate       (variance     (1),varianceLarge,time ,node)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Excursion set problem")
  do j=1,4
     select case (j)
     case (1)
        excursionSetFirstCrossing_ => excursionSetFirstCrossingFarahi_
        testRates                  =  .true.
        label                      =  'Farahi'
     case (2)
        excursionSetFirstCrossing_ => excursionSetFirstCrossingFarahiMidpoint_
        testRates                  =  .true.
        label                      =  'Farahi midpoint'
     case (3)
        excursionSetFirstCrossing_ => excursionSetFirstCrossingZhangHui_
        testRates                  =  .false.
        label                      =  'Zhang-Hui'
     case (4)
        excursionSetFirstCrossing_ => excursionSetFirstCrossingZhangHuiHighOrder_
        testRates                  =  .false.
        label                      =  'Zhang-Hui high order'
     end select
     call Unit_Tests_Begin_Group(trim(label)//" solver")
     do i=1,size(variance)
        probabilityTarget(i)=excursionSetFirstCrossingLinearBarrier_%probability(variance(i)                      ,time,node)
        probability      (i)=excursionSetFirstCrossing_             %probability(variance(i)                      ,time,node)
        if (testRates) then
           rateTarget    (i)=excursionSetFirstCrossingLinearBarrier_%rate       (variance(i),varianceProgenitor(i),time,node)
           rate          (i)=excursionSetFirstCrossing_             %rate       (variance(i),varianceProgenitor(i),time,node)
        end if
     end do
     call                Assert('Linear barrier crossing probability',probability,probabilityTarget,relTol=1.0d-2)
     if (testRates) call Assert('Linear barrier crossing rates'      ,rate       ,rateTarget       ,relTol=2.0d-2)
     call Unit_Tests_End_Group()
  end do
  ! Check that a tabulation read back from file reproduces the tabulation from which it was written. Two otherwise identical
  ! objects share a file, which is first removed, so that the first object must build and store its tabulations, and the second
  ! must read them.
  call Unit_Tests_Begin_Group("Farahi midpoint tabulation file round trip")
  fileNameRoundTrip=inputPath(pathTypeDataDynamic)//'largeScaleStructure/excursionSets/testsExcursionSetsRoundTrip.hdf5'
  if (File_Exists(fileNameRoundTrip)) call File_Remove(fileNameRoundTrip)
  excursionSetFirstCrossingFarahiMidpointWrite_=excursionSetFirstCrossingFarahiMidpoint(                                                                            &
       &                                                                                fractionalTimeStep                =0.01d0                                 , &
       &                                                                                fileName                          =fileNameRoundTrip                      , &
       &                                                                                varianceNumberPerUnitProbability  = 100                                   , &
       &                                                                                varianceNumberPerUnit             =   8                                   , &
       &                                                                                varianceNumberPerDecade           =  16                                   , &
       &                                                                                varianceNumberPerDecadeNonCrossing=   8                                   , &
       &                                                                                varianceIsUnlimited               =.true.                                 , &
       &                                                                                timeNumberPerDecade               =  10                                   , &
       &                                                                                cosmologyFunctions_               =cosmologyFunctionsMatterLambda_        , &
       &                                                                                excursionSetBarrier_              =excursionSetBarrierCriticalOverdensity_, &
       &                                                                                cosmologicalMassVariance_         =cosmologicalMassVarianceFilteredPower_   &
       &                                                                               )
  excursionSetFirstCrossingFarahiMidpointRead_ =excursionSetFirstCrossingFarahiMidpoint(                                                                            &
       &                                                                                fractionalTimeStep                =0.01d0                                 , &
       &                                                                                fileName                          =fileNameRoundTrip                      , &
       &                                                                                varianceNumberPerUnitProbability  = 100                                   , &
       &                                                                                varianceNumberPerUnit             =   8                                   , &
       &                                                                                varianceNumberPerDecade           =  16                                   , &
       &                                                                                varianceNumberPerDecadeNonCrossing=   8                                   , &
       &                                                                                varianceIsUnlimited               =.true.                                 , &
       &                                                                                timeNumberPerDecade               =  10                                   , &
       &                                                                                cosmologyFunctions_               =cosmologyFunctionsMatterLambda_        , &
       &                                                                                excursionSetBarrier_              =excursionSetBarrierCriticalOverdensity_, &
       &                                                                                cosmologicalMassVariance_         =cosmologicalMassVarianceFilteredPower_   &
       &                                                                               )
  ! Evaluate the first object at a large variance so that its probability tabulation spans the full range, and is written to
  ! file, then record its probabilities as computed directly.
  probability(1)=excursionSetFirstCrossingFarahiMidpointWrite_%probability(varianceLarge,time,node)
  do i=1,size(variance)
     probability(i)=excursionSetFirstCrossingFarahiMidpointWrite_%probability(variance(i),time,node)
  end do
  ! Tabulate rates using the first object. This reads the file (which now exists) before building the rate tabulation, so the
  ! first object's probability tabulation must be unaffected by having been read back from file.
  rate(1)=excursionSetFirstCrossingFarahiMidpointWrite_%rate(variance(1),varianceLarge,time,node)
  do i=1,size(variance)
     probabilityRates(i)=excursionSetFirstCrossingFarahiMidpointWrite_%probability(variance(i)                      ,time,node)
     rate            (i)=excursionSetFirstCrossingFarahiMidpointWrite_%rate       (variance(i),varianceProgenitor(i),time,node)
  end do
  ! Evaluate the second object, which reads its tabulations from file.
  do i=1,size(variance)
     probabilityRead(i)=excursionSetFirstCrossingFarahiMidpointRead_%probability(variance(i)                      ,time,node)
     rateRead       (i)=excursionSetFirstCrossingFarahiMidpointRead_%rate       (variance(i),varianceProgenitor(i),time,node)
  end do
  call Assert('Crossing probability after reading rates from file',probabilityRates,probability,relTol=1.0d-12)
  call Assert('Crossing probability read from file'               ,probabilityRead ,probability,relTol=1.0d-12)
  call Assert('Crossing rates read from file'                     ,rateRead        ,rate       ,relTol=1.0d-12)
  call Unit_Tests_End_Group()
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_Excursion_Sets
