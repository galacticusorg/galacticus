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
Contains a program which tests the \cite{correa_accretion_2015} halo mass formation history.
!!}

program Test_Correa2015_MAH
  !!{
  Tests the \cite{correa_accretion_2015} halo mass formation history algorithm.
  !!}
  use :: Cosmological_Density_Field               , only : cosmologicalMassVarianceFilteredPower
  use :: Cosmology_Parameters                     , only : cosmologyParametersSimple
  use :: Cosmology_Functions                      , only : cosmologyFunctionsMatterLambda
  use :: Dark_Matter_Particles                    , only : darkMatterParticleCDM
  use :: Dark_Matter_Halo_Mass_Accretion_Histories, only : darkMatterHaloMassAccretionHistoryCorrea2015
  use :: Display                                  , only : displayVerbositySet                         , verbosityLevelStandard
  use :: Events_Hooks                             , only : eventsHooksInitialize
  use :: Functions_Global_Utilities               , only : Functions_Global_Set
  use :: Galacticus_Nodes                         , only : nodeClassHierarchyInitialize                , nodeComponentBasic               , treeNode
  use :: Input_Parameters                         , only : inputParameters
  use :: Node_Components                          , only : Node_Components_Initialize                  , Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Unit_Tests                               , only : Assert                                      , Unit_Tests_Begin_Group           , Unit_Tests_End_Group               , Unit_Tests_Finish
  use :: Linear_Growth                            , only : linearGrowthCollisionlessMatter
  use :: Power_Spectra_Primordial                 , only : powerSpectrumPrimordialPowerLaw
  use :: Power_Spectra_Primordial_Transferred     , only : powerSpectrumPrimordialTransferredSimple
  use :: Power_Spectrum_Window_Functions          , only : powerSpectrumWindowFunctionTopHat
  use :: Transfer_Functions                       , only : transferFunctionEisensteinHu1999
  implicit none  
  type            (treeNode                                    ), pointer      :: node
  class           (nodeComponentBasic                          ), pointer      :: basic
  type            (darkMatterHaloMassAccretionHistoryCorrea2015)               :: darkMatterHaloMassAccretionHistory_
  type            (cosmologyParametersSimple                   )               :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda              )               :: cosmologyFunctions_
  type            (linearGrowthCollisionlessMatter             )               :: linearGrowth_
  type            (cosmologicalMassVarianceFilteredPower       )               :: cosmologicalMassVariance_
  type            (powerSpectrumWindowFunctionTopHat           )               :: powerSpectrumWindowFunction_
  type            (powerSpectrumPrimordialPowerLaw             )               :: powerSpectrumPrimordial_
  type            (transferFunctionEisensteinHu1999            )               :: transferFunction_
  type            (powerSpectrumPrimordialTransferredSimple    )               :: powerSpectrumPrimordialTransferred_
  type            (darkMatterParticleCDM                       )               :: darkMatterParticle_
  double precision                                              , dimension(3) :: time                               , mass          , &
       &                                                                          redshift                           , redshiftTarget, &
       &                                                                          massRecovered
  integer                                                                      :: i
  type            (inputParameters                             )               :: parameters
  
  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Correa et al. 2015 mass accretion history algorithms")
  ! Test Correa et al. 2015 algorithm.
  parameters=inputParameters()
  call eventsHooksInitialize()
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  ! Create a node.
  node  => treeNode      (                 )
  ! Get the basic component.
  basic => node    %basic(autoCreate=.true.)
  ! Get required objects.
  !![
  <referenceConstruct object="darkMatterParticle_"                >
   <constructor>
    darkMatterParticleCDM                                        (                                                                        &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyParameters_"               >
   <constructor>
    cosmologyParametersSimple                                    (                                                                        &amp;
     &amp;                                                        OmegaMatter                        =0.2815d0                          , &amp;
     &amp;                                                        OmegaBaryon                        =0.0465d0                          , &amp;
     &amp;                                                        OmegaDarkEnergy                    =0.7185d0                          , &amp;
     &amp;                                                        temperatureCMB                     =2.72548d0                         , &amp;
     &amp;                                                        HubbleConstant                     =69.3d0                              &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctions_"                >
   <constructor>
    cosmologyFunctionsMatterLambda                               (                                                                         &amp;
     &amp;                                                        cosmologyParameters_               =cosmologyParameters_                 &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="linearGrowth_"                      >
   <constructor>
    linearGrowthCollisionlessMatter                              (                                                                         &amp;
     &amp;                                                        cosmologyParameters_               =cosmologyParameters_               , &amp;
     &amp;                                                        cosmologyFunctions_                =cosmologyFunctions_                  &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumPrimordial_"           >
   <constructor>
    powerSpectrumPrimordialPowerLaw                              (                                                                         &amp;
     &amp;                                                        index_                             =0.971d0                            , &amp;
     &amp;                                                        running                            =+0.0d0                             , &amp;
     &amp;                                                        runningRunning                     =+0.0d0                             , &amp;
     &amp;                                                        wavenumberReference                =+1.0d0                             , &amp;
     &amp;                                                        runningSmallScalesOnly             =.false.                              &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="transferFunction_"                  >
   <constructor>
    transferFunctionEisensteinHu1999                             (                                                                         &amp;
     &amp;                                                        neutrinoNumberEffective            =3.046d0                            , &amp;
     &amp;                                                        neutrinoMassSummed                 =0.0d0                              , &amp;
     &amp;                                                        darkMatterParticle_                =darkMatterParticle_                , &amp;
     &amp;                                                        cosmologyParameters_               =cosmologyParameters_               , &amp;
     &amp;                                                        cosmologyFunctions_                =cosmologyFunctions_                  &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumPrimordialTransferred_">
   <constructor>
    powerSpectrumPrimordialTransferredSimple                     (                                                                         &amp;
     &amp;                                                        powerSpectrumPrimordial_           =powerSpectrumPrimordial_           , &amp;
     &amp;                                                        transferFunction_                  =transferFunction_                  , &amp;
     &amp;                                                        linearGrowth_                      =linearGrowth_                        &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumWindowFunction_"       >
   <constructor>
    powerSpectrumWindowFunctionTopHat                            (                                                                         &amp;
     &amp;                                                        cosmologyParameters_               =cosmologyParameters_                 &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologicalMassVariance_"          >
   <constructor>
    cosmologicalMassVarianceFilteredPower                        (                                                                         &amp;
     &amp;                                                        sigma8                             =0.82d0                             , &amp;
     &amp;                                                        tolerance                          =4.0d-6                             , &amp;
     &amp;                                                        toleranceTopHat                    =1.0d-6                             , &amp;
     &amp;                                                        nonMonotonicIsFatal                =.true.                             , &amp;
     &amp;                                                        monotonicInterpolation             =.false.                            , &amp;
     &amp;                                                        truncateAtParticleHorizon          =.false.                            , &amp;
     &amp;                                                        cosmologyParameters_               =cosmologyParameters_               , &amp;
     &amp;                                                        cosmologyFunctions_                =cosmologyFunctions_                , &amp;
     &amp;                                                        linearGrowth_                      =linearGrowth_                      , &amp;
     &amp;                                                        powerSpectrumPrimordialTransferred_=powerSpectrumPrimordialTransferred_, &amp;
     &amp;                                                        powerSpectrumWindowFunction_       =powerSpectrumWindowFunction_         &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterHaloMassAccretionHistory_">
   <constructor>
    darkMatterHaloMassAccretionHistoryCorrea2015                 (                                                                         &amp;
     &amp;                                                        cosmologyFunctions_                =cosmologyFunctions_                , &amp;
     &amp;                                                        linearGrowth_                      =linearGrowth_                      , &amp;
     &amp;                                                        cosmologicalMassVariance_          =cosmologicalMassVariance_            &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  !!]
  ! Specify halo masses and redshifts.
  mass=          [                  &
       &          1.00000000000d12, &
       &          5.64259219617d11, &
       &          2.96360394164d11  &
       &         ]
  redshiftTarget=[                  &
       &          0.00000000000d00, &
       &          1.00000000000d00, &
       &          2.00000000000d00  &
       &         ]
  ! Set node properties.
  call basic%massSet(                                                                    &
       &                                                               mass          (1) &
       &            )
  call basic%timeSet(                                                                    &
       &             cosmologyFunctions_ %cosmicTime                 (                   &
       &              cosmologyFunctions_%expansionFactorFromRedshift (                  &
       &                                                               redshiftTarget(1) &
       &                                                              )                  &
       &                                                             )                   &
       &            )
  ! Compute the mass accretion history.
  do i=1,size(mass)
     time         (i)=darkMatterHaloMassAccretionHistory_%time(node,mass(i))
     massRecovered(i)=darkMatterHaloMassAccretionHistory_%mass(node,time(i))
     redshift     (i)=cosmologyFunctions_ %redshiftFromExpansionFactor(         &
          &            cosmologyFunctions_%expansionFactor             (        &
          &                                                             time(i) &
          &                                                            )        &
          &                                                           )
  end do
  call Assert('mass accretion history'          ,1.0d0+redshift,1.0d0+redshiftTarget,relTol=1.0d-3)
  call Assert('mass accretion history inversion',mass          ,massRecovered       ,relTol=1.0d-3)
  ! End unit tests.
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
end program Test_Correa2015_MAH
