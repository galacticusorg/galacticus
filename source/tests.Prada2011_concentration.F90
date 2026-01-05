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
Contains a program which tests the \cite{prada_halo_2011} halo concentration algorithm.
!!}

program Test_Prada2011_Concentration
  !!{
  Tests the \cite{prada_halo_2011} halo concentration algorithm. Values of concentration were read from their Figure~12.
  !!}
  use :: Cosmological_Density_Field          , only : cosmologicalMassVarianceFilteredPower
  use :: Cosmology_Functions                 , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters                , only : cosmologyParametersSimple               , hubbleUnitsLittleH
  use :: Dark_Matter_Particles               , only : darkMatterParticleCDM
  use :: Dark_Matter_Profiles_Concentration  , only : darkMatterProfileConcentrationPrada2011
  use :: Display                             , only : displayVerbositySet                     , verbosityLevelStandard
  use :: Events_Hooks                        , only : eventsHooksInitialize
  use :: Functions_Global_Utilities          , only : Functions_Global_Set
  use :: Galacticus_Nodes                    , only : nodeClassHierarchyInitialize            , nodeComponentBasic               , treeNode
  use :: Linear_Growth                       , only : linearGrowthCollisionlessMatter
  use :: ISO_Varying_String                  , only : assignment(=)                           , char                             , varying_string
  use :: Input_Parameters                    , only : inputParameters
  use :: Node_Components                     , only : Node_Components_Initialize              , Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Power_Spectra_Primordial            , only : powerSpectrumPrimordialPowerLaw
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredSimple
  use :: Power_Spectrum_Window_Functions     , only : powerSpectrumWindowFunctionTopHat
  use :: Transfer_Functions                  , only : transferFunctionEisensteinHu1999
  use :: Unit_Tests                          , only : Assert                                  , Unit_Tests_Begin_Group           , Unit_Tests_End_Group               , Unit_Tests_Finish
  implicit none
  integer                                                                         , parameter:: massCount                          =4
  double precision                                          , dimension(massCount), parameter:: logMass                            =[11.000d0,12.000d0,13.000d0,14.000d0]
  double precision                                          , dimension(massCount), parameter:: pradaLogConcentration              =[ 0.966d0, 0.887d0, 0.804d0, 0.728d0]
  double precision                                          , dimension(massCount)           :: ourLogConcentration
  type            (treeNode                                )                      , pointer  :: node
  class           (nodeComponentBasic                      )                      , pointer  :: basic
  type            (cosmologyParametersSimple               )                                 :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda          )                                 :: cosmologyFunctions_
  type            (darkMatterParticleCDM                   )                                 :: darkMatterParticle_
  type            (linearGrowthCollisionlessMatter         )                                 :: linearGrowth_
  type            (powerSpectrumWindowFunctionTopHat       )                                 :: powerSpectrumWindowFunction_
  type            (powerSpectrumPrimordialPowerLaw         )                                 :: powerSpectrumPrimordial_
  type            (transferFunctionEisensteinHu1999        )                                 :: transferFunction_
  type            (powerSpectrumPrimordialTransferredSimple)                                 :: powerSpectrumPrimordialTransferred_
  type            (cosmologicalMassVarianceFilteredPower   )                                 :: cosmologicalMassVariance_
  type            (darkMatterProfileConcentrationPrada2011 )                                 :: darkMatterProfileConcentration_
  integer                                                                                    :: iMass
  type            (inputParameters                         )                                 :: parameters
  type            (varying_string                          )                                 :: message                                                                  , parameterFile

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Prada2011 halo concentration algorithm")
  ! Test Prada2011 halo concentration algorithm.
  ! Read in controlling parameters.
  parameterFile='testSuite/parameters/Prada2011HaloConcentration/testParameters.xml'
  parameters=inputParameters(parameterFile)
  call eventsHooksInitialize()
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  ! Create a node.
  node  => treeNode      (                 )
  basic => node    %basic(autoCreate=.true.)
  ! Construct required objects.
  !![
  <referenceConstruct object="darkMatterParticle_"                >
   <constructor>
    darkMatterParticleCDM                   (                                                                         &amp;
     &amp;                                  )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyParameters_"               >
   <constructor>
    cosmologyParametersSimple               (                                                                         &amp;
     &amp;                                   OmegaMatter                        = 0.2700d0                          , &amp;
     &amp;                                   OmegaBaryon                        = 0.0469d0                          , &amp;
     &amp;                                   OmegaDarkEnergy                    = 0.7300d0                          , &amp;
     &amp;                                   temperatureCMB                     = 2.7000d0                          , &amp;
     &amp;                                   HubbleConstant                     =70.0000d0                            &amp;
     &amp;                                  )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctions_"                >
   <constructor>
    cosmologyFunctionsMatterLambda          (                                                                         &amp;
     &amp;                                   cosmologyParameters_               =cosmologyParameters_                 &amp;
     &amp;                                  )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="linearGrowth_"                      >
   <constructor>
    linearGrowthCollisionlessMatter         (                                                                         &amp;
     &amp;                                   cosmologyParameters_               =cosmologyParameters_               , &amp;
     &amp;                                   cosmologyFunctions_                =cosmologyFunctions_                  &amp;
     &amp;                                  )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumPrimordial_"           >
   <constructor>
    powerSpectrumPrimordialPowerLaw         (                                                                         &amp;
     &amp;                                   index_                             =+0.95d0                            , &amp;
     &amp;                                   running                            =+0.00d0                            , &amp;
     &amp;                                   runningRunning                     =+0.00d0                            , &amp;
     &amp;                                   wavenumberReference                =+1.00d0                            , &amp;
     &amp;                                   runningSmallScalesOnly             =.false.                              &amp;
     &amp;                                  )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="transferFunction_"                  >
   <constructor>
    transferFunctionEisensteinHu1999        (                                                                         &amp;
     &amp;                                   neutrinoNumberEffective            =3.046d0                            , &amp;
     &amp;                                   neutrinoMassSummed                 =0.000d0                            , &amp;
     &amp;                                   darkMatterParticle_                =darkMatterParticle_                , &amp;
     &amp;                                   cosmologyParameters_               =cosmologyParameters_               , &amp;
     &amp;                                   cosmologyFunctions_                =cosmologyFunctions_                  &amp;
     &amp;                                  )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumPrimordialTransferred_">
   <constructor>
    powerSpectrumPrimordialTransferredSimple(                                                                         &amp;
     &amp;                                   powerSpectrumPrimordial_           =powerSpectrumPrimordial_           , &amp;
     &amp;                                   transferFunction_                  =transferFunction_                  , &amp;
     &amp;                                   linearGrowth_                      =linearGrowth_                        &amp;
     &amp;                                  )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumWindowFunction_"       >
   <constructor>
    powerSpectrumWindowFunctionTopHat       (                                                                         &amp;
     &amp;                                   cosmologyParameters_               =cosmologyParameters_                 &amp;
     &amp;                                  )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologicalMassVariance_"          >
   <constructor>
    cosmologicalMassVarianceFilteredPower   (                                                                         &amp;
     &amp;                                   sigma8                             =0.82d+0                            , &amp;
     &amp;                                   tolerance                          =1.00d-4                            , &amp;
     &amp;                                   toleranceTopHat                    =1.00d-4                            , &amp;
     &amp;                                   nonMonotonicIsFatal                =.true.                             , &amp;
     &amp;                                   monotonicInterpolation             =.false.                            , &amp;
     &amp;                                   truncateAtParticleHorizon          =.false.                            , &amp;
     &amp;                                   cosmologyParameters_               =cosmologyParameters_               , &amp;
     &amp;                                   cosmologyFunctions_                =cosmologyFunctions_                , &amp;
     &amp;                                   linearGrowth_                      =linearGrowth_                      , &amp;
     &amp;                                   powerSpectrumPrimordialTransferred_=powerSpectrumPrimordialTransferred_, &amp;
     &amp;                                   powerSpectrumWindowFunction_       =powerSpectrumWindowFunction_         &amp;
     &amp;                                  )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfileConcentration_"                  >
   <constructor>
    darkMatterProfileConcentrationPrada2011 (                                                                         &amp;
     &amp;                                   A                                  =2.881d0                            , &amp;
     &amp;                                   B                                  =1.257d0                            , &amp;
     &amp;                                   C                                  =1.022d0                            , &amp;
     &amp;                                   D                                  =0.060d0                            , &amp;
     &amp;                                   C0                                 =3.681d0                            , &amp;
     &amp;                                   C1                                 =5.033d0                            , &amp;
     &amp;                                   X0                                 =0.424d0                            , &amp;
     &amp;                                   X1                                 =0.526d0                            , &amp;
     &amp;                                   inverseSigma0                      =1.047d0                            , &amp;
     &amp;                                   inverseSigma1                      =1.646d0                            , &amp;
     &amp;                                   alpha                              =6.948d0                            , &amp;
     &amp;                                   beta                               =7.386d0                            , &amp;
     &amp;                                   cosmologyParameters_               =cosmologyParameters_               , &amp;
     &amp;                                   cosmologyFunctions_                =cosmologyFunctions_                , &amp;
     &amp;                                   cosmologicalMassVariance_          =cosmologicalMassVariance_            &amp;
     &amp;                                  )
   </constructor>
  </referenceConstruct>
  !!]
  ! Set the time for the node.
  call basic%timeSet(cosmologyFunctions_%cosmicTime(1.00d0))
  ! Loop over halo masses
  do iMass=1,massCount
     ! Set the mass of the original node.
     call basic%massSet(10.0d0**logMass(iMass)/cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH))
     ! Compute and compare concentration at z=0.
     ourLogConcentration(iMass)=log10(darkMatterProfileConcentration_%concentration(node))
  end do
  ! Check that results are as expected.
  message="Halo concentration at z=0"
  call Assert(char(message),ourLogConcentration,pradaLogConcentration,absTol=0.01d0)
  ! End unit tests.
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
end program Test_Prada2011_Concentration
