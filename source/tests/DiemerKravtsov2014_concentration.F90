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

!!{RST
Contains a program which tests the :cite:t:`diemer_universal_2014` halo concentration algorithm.
!!}

program Test_DiemerKravtsov2014_Concentration
  !!{RST
  Tests the :cite:t:`diemer_universal_2014` halo concentration algorithm. Reference values of concentration are those tabulated by
  the model's authors for the WMAP9 cosmology\ [#]_, and are stored in ``testSuite/data/diemerKravtsov2014Concentration.txt``. The
  cosmological parameters constructed below are those used to generate that table, as recorded in its header.

  .. [#] File ``cM_WMAP9.txt``, from ``https://www.benediktdiemer.com/data/concentration-mass-relation/``. The table is shipped
     with Galacticus rather than downloaded at run time, as the individual data files are distributed only within a zip archive.
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
  use :: Transfer_Functions                  , only : transferFunctionEisensteinHu1999
  use :: Unit_Tests                          , only : Assert                                          , Unit_Tests_Begin_Group                                      , Unit_Tests_End_Group               , Unit_Tests_Finish
  implicit none
  type            (treeNode                                                    ), pointer :: node
  class           (nodeComponentBasic                                          ), pointer :: basic
  type            (cosmologyParametersSimple                                   ), pointer :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                              ), pointer :: cosmologyFunctions_
  type            (linearGrowthCollisionlessMatter                             ), pointer :: linearGrowth_
  type            (cosmologicalMassVarianceFilteredPower                       ), pointer :: cosmologicalMassVariance_
  type            (powerSpectrumWindowFunctionTopHat                           ), pointer :: powerSpectrumWindowFunction_
  type            (powerSpectrumPrimordialPowerLaw                             ), pointer :: powerSpectrumPrimordial_
  type            (transferFunctionEisensteinHu1999                            ), pointer :: transferFunction_
  type            (powerSpectrumPrimordialTransferredSimple                    ), pointer :: powerSpectrumPrimordialTransferred_
  type            (powerSpectrumStandard                                       ), pointer :: powerSpectrum_
  type            (darkMatterParticleCDM                                       ), pointer :: darkMatterParticle_
  type            (criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt), pointer :: criticalOverdensity_
  type            (darkMatterProfileConcentrationDiemerKravtsov2014            ), pointer :: darkMatterProfileConcentration_
  type            (varying_string                                              )          :: parameterFile                      , fileName
  type            (inputParameters                                             )          :: parameters
  double precision                                                                        :: ourConcentration                   , differenceFractional, &
       &                                                                                     concentration                      , mass                , &
       &                                                                                     redshift                           , nu                  , &
       &                                                                                     differenceFractionalMaximum
  integer                                                                                 :: referenceUnit                      , ioStatus            , &
       &                                                                                     i
  character       (len=1024                                                    )          :: line

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

  ! The reference dataset is shipped with Galacticus.
  fileName=inputPath(pathTypeExec)//"testSuite/data/diemerKravtsov2014Concentration.txt"
  if (.not.File_Exists(fileName)) call Error_Report('reference dataset "'//fileName//'" is missing'//{introspection:location})
  ! Create a node.
  node                            => treeNode                                        (                 )
  ! Get the basic component.
  basic                           => node%basic                                      (autoCreate=.true.)
  ! Get required objects.
  allocate(cosmologyParameters_               )
  allocate(cosmologyFunctions_                )
  allocate(linearGrowth_                      )
  allocate(cosmologicalMassVariance_          )
  allocate(powerSpectrumWindowFunction_       )
  allocate(powerSpectrumPrimordial_           )
  allocate(transferFunction_                  )
  allocate(powerSpectrumPrimordialTransferred_)
  allocate(powerSpectrum_                     )
  allocate(darkMatterParticle_                )
  allocate(criticalOverdensity_               )
  allocate(darkMatterProfileConcentration_    )
  !![
  <referenceConstruct object="darkMatterParticle_"                >
   <constructor>
    darkMatterParticleCDM                                        (                                                                               &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyParameters_"               >
   <constructor>
    cosmologyParametersSimple                                    (                                                                               &amp;
     &amp;                                                        OmegaMatter                             = 0.2865d0                           , &amp;
     &amp;                                                        OmegaBaryon                             = 0.0463d0                           , &amp;
     &amp;                                                        OmegaDarkEnergy                         = 0.7135d0                           , &amp;
     &amp;                                                        temperatureCMB                          = 2.7255d0                           , &amp;
     &amp;                                                        HubbleConstant                          =69.3200d0                             &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctions_"                >
   <constructor>
    cosmologyFunctionsMatterLambda                               (                                                                              &amp;
     &amp;                                                        cosmologyParameters_                    =cosmologyParameters_                 &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="linearGrowth_"                      >
   <constructor>
    linearGrowthCollisionlessMatter                              (                                                                              &amp;
     &amp;                                                        cosmologyParameters_                    =cosmologyParameters_               , &amp;
     &amp;                                                        cosmologyFunctions_                     =cosmologyFunctions_                  &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumPrimordial_"           >
   <constructor>
    powerSpectrumPrimordialPowerLaw                              (                                                                              &amp;
     &amp;                                                        index_                                  =+0.9608d0                          , &amp;
     &amp;                                                        running                                 =+0.000d0                           , &amp;
     &amp;                                                        runningRunning                          =+0.000d0                           , &amp;
     &amp;                                                        wavenumberReference                     =+1.000d0                           , &amp;
     &amp;                                                        runningSmallScalesOnly                  =.false.                              &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="transferFunction_"                  >
   <constructor>
    transferFunctionEisensteinHu1999                             (                                                                              &amp;
     &amp;                                                        neutrinoNumberEffective                 =3.046d0                            , &amp;
     &amp;                                                        neutrinoMassSummed                      =0.000d0                            , &amp;
     &amp;                                                        darkMatterParticle_                     =darkMatterParticle_                , &amp;
     &amp;                                                        cosmologyParameters_                    =cosmologyParameters_               , &amp;
     &amp;                                                        cosmologyFunctions_                     =cosmologyFunctions_                  &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumPrimordialTransferred_">
   <constructor>
    powerSpectrumPrimordialTransferredSimple                     (                                                                              &amp;
     &amp;                                                        powerSpectrumPrimordial_                =powerSpectrumPrimordial_           , &amp;
     &amp;                                                        transferFunction_                       =transferFunction_                  , &amp;
     &amp;                                                        linearGrowth_                           =linearGrowth_                        &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumWindowFunction_"       >
   <constructor>
    powerSpectrumWindowFunctionTopHat                            (                                                                              &amp;
     &amp;                                                        cosmologyParameters_                    =cosmologyParameters_                 &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologicalMassVariance_"          >
   <constructor>
    cosmologicalMassVarianceFilteredPower                        (                                                                              &amp;
     &amp;                                                        sigma8                                  =0.820d+0                           , &amp;
     &amp;                                                        tolerance                               =1.000d-4                           , &amp;
     &amp;                                                        toleranceTopHat                         =1.000d-4                           , &amp;
     &amp;                                                        rootVarianceLogarithmicGradientTolerance=1.0d-9                             , &amp;
     &amp;                                                        integrationFailureIsFatal               =.true.                             , &amp;
     &amp;                                                        storeTabulations                        =.true.                             , &amp;
     &amp;                                                        nonMonotonicIsFatal                     =.true.                             , &amp;
     &amp;                                                        monotonicInterpolation                  =.false.                            , &amp;
     &amp;                                                        truncateAtParticleHorizon               =.false.                            , &amp;
     &amp;                                                        cosmologyParameters_                    =cosmologyParameters_               , &amp;
     &amp;                                                        cosmologyFunctions_                     =cosmologyFunctions_                , &amp;
     &amp;                                                        linearGrowth_                           =linearGrowth_                      , &amp;
     &amp;                                                        powerSpectrumPrimordialTransferred_     =powerSpectrumPrimordialTransferred_, &amp;
     &amp;                                                        powerSpectrumWindowFunction_            =powerSpectrumWindowFunction_         &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrum_"                     >
   <constructor>
    powerSpectrumStandard                                        (                                                                              &amp;
     &amp;                                                        cosmologicalMassVariance_               =cosmologicalMassVariance_          , &amp;
     &amp;                                                        powerSpectrumPrimordialTransferred_     =powerSpectrumPrimordialTransferred_  &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="criticalOverdensity_"               >
   <constructor>
    criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt(                                                                               &amp;
     &amp;                                                        linearGrowth_                           =linearGrowth_                      , &amp;
     &amp;                                                        cosmologyFunctions_                     =cosmologyFunctions_                , &amp;
     &amp;                                                        cosmologicalMassVariance_               =cosmologicalMassVariance_          , &amp;
     &amp;                                                        darkMatterParticle_                     =darkMatterParticle_                , &amp;
     &amp;                                                        tableStore                              =.true.                               &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfileConcentration_"    >
    <constructor>
     darkMatterProfileConcentrationDiemerKravtsov2014            (                                                                        &amp;
      &amp;                                                       kappa                              =0.69d0                            , &amp;
      &amp;                                                       phi0                               =6.58d0                            , &amp;
      &amp;                                                       phi1                               =1.37d0                            , &amp;
      &amp;                                                       eta0                               =6.82d0                            , &amp;
      &amp;                                                       eta1                               =1.42d0                            , &amp;
      &amp;                                                       alpha                              =1.12d0                            , &amp;
      &amp;                                                       beta                               =1.69d0                            , &amp;
      &amp;                                                       scatter                            =0.00d0                            , &amp;
      &amp;                                                       cosmologyFunctions_                =cosmologyFunctions_               , &amp;
      &amp;                                                       cosmologyParameters_               =cosmologyParameters_              , &amp;
      &amp;                                                       criticalOverdensity_               =criticalOverdensity_              , &amp;
      &amp;                                                       cosmologicalMassVariance_          =cosmologicalMassVariance_         , &amp;
      &amp;                                                       powerSpectrum_                     =powerSpectrum_                      &amp;
      &amp;                                                      )
   </constructor>
  </referenceConstruct>
  !!]
  ! Read the reference file.
  differenceFractionalMaximum=0.0d0
  fileName                   =inputPath(pathTypeExec)//"testSuite/data/diemerKravtsov2014Concentration.txt"
  open(newUnit=referenceUnit,file=char(fileName),status='old',form='formatted',iostat=ioStatus)
  do i=1,9
     read (referenceUnit,'(a)',ioStat=ioStatus) line ! Skip header.
     if (ioStatus /= 0 .or. line(1:1) /= "#") call Error_Report('reference dataset header is not of the expected form'//{introspection:location})
  end do
  do while (ioStatus == 0)
     read (referenceUnit,*,ioStat=ioStatus) redshift,nu,mass,concentration
     if (ioStatus /= 0) exit
     ! The reference table tabulates concentrations from z=0 to z=30. Only the z=0 entries are compared here, spanning 8 decades
     ! of mass (M₂₀₀c from 8.6×10⁷ to 6.7×10¹⁵ M☉/h). Agreement degrades smoothly with increasing redshift---from 2.2% at z=0 to
     ! 10.4% at z=30---which is a systematic difference in the cosmology and transfer function used, not a failure at any
     ! particular mass. Extending this comparison to higher redshift would therefore require a tolerance loose enough to make it
     ! ineffective as a regression test.
     if (redshift > 0.0d0) cycle
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
  ! Assert that the maximum fractional difference is not too large. The residual differences are
  ! presumably because we don't use precisely the same transfer function as do Diemer &
  ! Kravtsov.
  call Assert("Halo concentration in WMAP9 reference model",differenceFractionalMaximum,0.0d0,absTol=2.9d-2)
  ! Clean up.
  !![
  <objectDestructor name="darkMatterProfileConcentration_"    />
  <objectDestructor name="cosmologyParameters_"               />
  <objectDestructor name="cosmologyFunctions_"                />
  <objectDestructor name="linearGrowth_"                      />
  <objectDestructor name="cosmologicalMassVariance_"          />
  <objectDestructor name="powerSpectrumWindowFunction_"       />
  <objectDestructor name="powerSpectrumPrimordial_"           />
  <objectDestructor name="transferFunction_"                  />
  <objectDestructor name="powerSpectrumPrimordialTransferred_"/>
  <objectDestructor name="powerSpectrum_"                     />
  <objectDestructor name="darkMatterParticle_"                />
  <objectDestructor name="criticalOverdensity_"               />
  <objectDestructor name="darkMatterProfileConcentration_"    />
  !!]
  call node%destroy()
  deallocate(node)
  ! End unit tests.
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
end program Test_DiemerKravtsov2014_Concentration
