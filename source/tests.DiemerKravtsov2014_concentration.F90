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
Contains a program which tests the \cite{diemer_universal_2014} halo concentration algorithm.
!!}

program Test_DiemerKravtsov2014_Concentration
  !!{
  Tests the \cite{diemer_universal_2014} halo concentration algorithm. Values of concentration were taken from their website\footnote{File no longer available---was downloaded from {\normalfont \ttfamily http://www.benediktdiemer.com/wp-content/uploads/2014/07/Concentration\_WMAP7\_median.txt}}.
  !!}
  use :: Cosmological_Density_Field          , only : cosmologicalMassVarianceFilteredPower           , criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt
  use :: Cosmology_Functions                 , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters                , only : cosmologyParametersSimple                       , hubbleUnitsLittleH
  use :: Dark_Matter_Particles               , only : darkMatterParticleCDM
  use :: Dark_Matter_Profiles_Concentration  , only : darkMatterProfileConcentrationDiemerKravtsov2014
  use :: Display                             , only : displayVerbositySet                             , verbosityLevelStandard
  use :: Events_Hooks                        , only : eventsHooksInitialize
  use :: File_Utilities                      , only : File_Exists
  use :: Functions_Global_Utilities          , only : Functions_Global_Set
  use :: Error                               , only : Error_Report
  use :: Galacticus_Nodes                    , only : nodeClassHierarchyInitialize                    , nodeComponentBasic                                          , treeNode
  use :: Input_Paths                         , only : inputPath                                       , pathTypeExec
  use :: ISO_Varying_String                  , only : assignment(=)                                   , char                                                        , operator(//)                       , varying_string
  use :: Input_Parameters                    , only : inputParameters
  use :: Linear_Growth                       , only : linearGrowthCollisionlessMatter
  use :: Node_Components                     , only : Node_Components_Initialize                      , Node_Components_Thread_Initialize                           , Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Power_Spectra                       , only : powerSpectrumStandard
  use :: Power_Spectra_Primordial            , only : powerSpectrumPrimordialPowerLaw
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredSimple
  use :: Power_Spectrum_Window_Functions     , only : powerSpectrumWindowFunctionTopHat
  use :: System_Download                     , only : download
  use :: Transfer_Functions                  , only : transferFunctionEisensteinHu1999
  use :: Unit_Tests                          , only : Assert                                          , Unit_Tests_Begin_Group                                      , Unit_Tests_End_Group               , Unit_Tests_Finish
  implicit none
  type            (treeNode                                                    ), pointer :: node
  class           (nodeComponentBasic                                          ), pointer :: basic
  type            (cosmologyParametersSimple                                   )          :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                              )          :: cosmologyFunctions_
  type            (linearGrowthCollisionlessMatter                             )          :: linearGrowth_
  type            (cosmologicalMassVarianceFilteredPower                       )          :: cosmologicalMassVariance_
  type            (powerSpectrumWindowFunctionTopHat                           )          :: powerSpectrumWindowFunction_
  type            (powerSpectrumPrimordialPowerLaw                             )          :: powerSpectrumPrimordial_
  type            (transferFunctionEisensteinHu1999                            )          :: transferFunction_
  type            (powerSpectrumPrimordialTransferredSimple                    )          :: powerSpectrumPrimordialTransferred_
  type            (powerSpectrumStandard                                       )          :: powerSpectrum_
  type            (darkMatterParticleCDM                                       )          :: darkMatterParticle_
  type            (criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt)          :: criticalOverdensity_
  type            (darkMatterProfileConcentrationDiemerKravtsov2014            )          :: darkMatterProfileConcentration_
  type            (varying_string                                              )          :: parameterFile
  type            (inputParameters                                             )          :: parameters
  double precision                                                                        :: ourConcentration                   , differenceFractional, &
       &                                                                                     concentration                      , mass                , &
       &                                                                                     redshift                           , nu                  , &
       &                                                                                     differenceFractionalMaximum
  integer                                                                                 :: referenceUnit                      , ioStatus            , &
       &                                                                         i

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("DiemerKravtsov2014 halo concentration algorithm")
  ! Test DiemerKravtsov2014 halo concentration algorithm.
  ! Read in controlling parameters.
  parameterFile='testSuite/parameters/DiemerKravtsov2014HaloConcentration/testParameters.xml'
  parameters=inputParameters(parameterFile)
  call eventsHooksInitialize()
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)

  ! Get the data file if we don't have it.
  if (.not.File_Exists(inputPath(pathTypeExec)//"testSuite/data/diemerKravtsov2014Concentration.txt")) then
     call download(                                                                  &
          &        "http://www.benediktdiemer.com/wp-content/uploads/cM_WMAP7.txt",  &
          &        char(inputPath(pathTypeExec))                                  // &
          &        "testSuite/data/diemerKravtsov2014Concentration.txt"           ,  &
          &        status=ioStatus                                                   &
          &       )
     if (ioStatus /= 0 .or. .not.File_Exists(inputPath(pathTypeExec)//"testSuite/data/diemerKravtsov2014Concentration.txt")) &
          & call Error_Report('unable to retrieve reference dataset'//{introspection:location})
  end if
  ! Create a node.
  node                            => treeNode                                        (                 )
  ! Get the basic component.
  basic                           => node%basic                                      (autoCreate=.true.)
  ! Get required objects.
  !![
  <referenceConstruct object="darkMatterParticle_"                >
   <constructor>
    darkMatterParticleCDM                                        (                                                                          &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyParameters_"               >
   <constructor>
    cosmologyParametersSimple                                    (                                                                          &amp;
     &amp;                                                        OmegaMatter                        = 0.2743d0                           , &amp;
     &amp;                                                        OmegaBaryon                        = 0.0458d0                           , &amp;
     &amp;                                                        OmegaDarkEnergy                    = 0.7257d0                           , &amp;
     &amp;                                                        temperatureCMB                     = 2.7000d0                           , &amp;
     &amp;                                                        HubbleConstant                     =70.2000d0                             &amp;
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
     &amp;                                                        index_                             =+0.968d0                           , &amp;
     &amp;                                                        running                            =+0.000d0                           , &amp;
     &amp;                                                        runningRunning                     =+0.000d0                           , &amp;
     &amp;                                                        wavenumberReference                =+1.000d0                           , &amp;
     &amp;                                                        runningSmallScalesOnly             =.false.                              &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="transferFunction_"                  >
   <constructor>
    transferFunctionEisensteinHu1999                             (                                                                         &amp;
     &amp;                                                        neutrinoNumberEffective            =3.046d0                            , &amp;
     &amp;                                                        neutrinoMassSummed                 =0.000d0                            , &amp;
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
     &amp;                                                        sigma8                             =0.816d+0                           , &amp;
     &amp;                                                        tolerance                          =1.000d-4                           , &amp;
     &amp;                                                        toleranceTopHat                    =1.000d-4                           , &amp;
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
  <referenceConstruct object="powerSpectrum_"                     >
   <constructor>
    powerSpectrumStandard                                        (                                                                         &amp;
     &amp;                                                        cosmologicalMassVariance_          =cosmologicalMassVariance_          , &amp;
     &amp;                                                        powerSpectrumPrimordialTransferred_=powerSpectrumPrimordialTransferred_  &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="criticalOverdensity_"               >
   <constructor>
    criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt(                                                                          &amp;
     &amp;                                                        linearGrowth_                      =linearGrowth_                      , &amp;
     &amp;                                                        cosmologyFunctions_                =cosmologyFunctions_                , &amp;
     &amp;                                                        cosmologicalMassVariance_          =cosmologicalMassVariance_          , &amp;
     &amp;                                                        darkMatterParticle_                =darkMatterParticle_                , &amp;
     &amp;                                                        tableStore                         =.true.                               &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  !!]
  darkMatterProfileConcentration_=darkMatterProfileConcentrationDiemerKravtsov2014(                          &
       &                                                                           0.69d0                   , &
       &                                                                           6.58d0                   , &
       &                                                                           1.37d0                   , &
       &                                                                           6.82d0                   , &
       &                                                                           1.42d0                   , &
       &                                                                           1.12d0                   , &
       &                                                                           1.69d0                   , &
       &                                                                           0.00d0                   , &
       &                                                                           cosmologyFunctions_      , &
       &                                                                           cosmologyParameters_     , &
       &                                                                           criticalOverdensity_     , &
       &                                                                           cosmologicalMassVariance_, &
       &                                                                           powerSpectrum_             &
       &                                                                          )
  ! Read the reference file.
  differenceFractionalMaximum=0.0d0
  open(newUnit=referenceUnit,file=char(inputPath(pathTypeExec)//"testSuite/data/diemerKravtsov2014Concentration.txt"),status='old',form='formatted',iostat=ioStatus)
  do i=1,7
     read (referenceUnit,*,ioStat=ioStatus) ! Skip header.
  end do
  do while (ioStatus == 0)
     read (referenceUnit,*,ioStat=ioStatus) redshift,nu,mass,concentration
     if (ioStatus /= 0) exit
     ! Set the time for the node.
     call basic%timeSet            (cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift)))
     call basic%timeLastIsolatedSet(cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift)))
     ! Set the mass of the original node (Diemer & Kravtsov masses are in units of M☉/h, so
     ! we convert from that system).
     call basic%massSet(mass/cosmologyParameters_%HubbleConstant(units=hubbleUnitsLittleH))
     ! Compute and compare concentration at z=0.
     ourConcentration           =darkMatterProfileConcentration_%concentration(node)
     differenceFractional       =abs(ourConcentration-concentration)/concentration
     differenceFractionalMaximum=max(differenceFractionalMaximum,differenceFractional)
  end do
  ! Assert that the maximum fractional difference is not too large. The ~3% differences are
  ! presumably because we don't use precisely the same transfer function as do Diemer &
  ! Kravtsov.
  call Assert("Halo concentration in WMAP7 reference model",differenceFractionalMaximum,0.0d0,absTol=2.9d-2)
  ! End unit tests.
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
end program Test_DiemerKravtsov2014_Concentration
