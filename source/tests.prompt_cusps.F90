!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
Contains a program which tests the \cite{delos_cusp-halo_2025} prompt cusps model.
!!}

program Test_Prompt_Cusps
  !!{
  Tests the \cite{delos_cusp-halo_2025} prompt cusps model. Values of prompt cusp properties were computed using the
  \href{https://github.com/delos/cusp-halo-relation}{cusp\_halo\_relation} Python module.  
  !!}
  use :: Coordinates                         , only : coordinateSpherical                                           , assignment(=)
  use :: Cosmological_Density_Field          , only : cosmologicalMassVarianceFilteredPower
  use :: Cosmology_Functions                 , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Functions_Parameters      , only : requestTypeTime
  use :: Cosmology_Parameters                , only : cosmologyParametersSimple
  use :: Dark_Matter_Particles               , only : darkMatterParticleCDM                                         , darkMatterParticleWDMThermal
  use :: Dark_Matter_Halo_Scales             , only : darkMatterHaloScaleVirialDensityContrastDefinition
  use :: Dark_Matter_Profiles_DMO            , only : darkMatterProfileDMOCuspNFW            
  use :: Display                             , only : displayVerbositySet                                           , verbosityLevelStandard           , displayMessage
  use :: Events_Hooks                        , only : eventsHooksInitialize
  use :: Functions_Global_Utilities          , only : Functions_Global_Set
  use :: Galacticus_Nodes                    , only : nodeClassHierarchyInitialize                                  , nodeComponentBasic               , nodeComponentDarkMatterProfile                                   , treeNode
  use :: ISO_Varying_String                  , only : assignment(=)                                                 , char                                                                  , operator(//)                , varying_string
  use :: Input_Parameters                    , only : inputParameters
  use :: Linear_Growth                       , only : linearGrowthCollisionlessMatter
  use :: Node_Components                     , only : Node_Components_Initialize                                    , Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Nodes_Operators                     , only : nodeOperatorDarkMatterProfilePromptCusps
  use :: Node_Property_Extractors            , only : nodePropertyExtractorPromptCusps
  use :: Mass_Distributions                  , only : massDistributionClass
  use :: Power_Spectra                       , only : powerSpectrumStandard
  use :: Power_Spectra_Primordial            , only : powerSpectrumPrimordialPowerLaw
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredSimple
  use :: Power_Spectrum_Window_Functions     , only : powerSpectrumWindowFunctionTopHat
  use :: Virial_Density_Contrast             , only : virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt
  use :: Transfer_Functions                  , only : transferFunctionCAMB                                          , transferFunctionBode2001         , scaleCutOffModelVogel23SpinHalf
  use :: Unit_Tests                          , only : Assert                                                        , Unit_Tests_Begin_Group                                                , Unit_Tests_End_Group        , Unit_Tests_Finish
  implicit none
  type            (treeNode                                                      ), pointer     :: node
  class           (nodeComponentBasic                                            ), pointer     :: basic
  class           (nodeComponentDarkMatterProfile                                ), pointer     :: darkMatterProfile
  class           (massDistributionClass                                         ), pointer     :: massDistribution_
  type            (cosmologyParametersSimple                                     )              :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                                )              :: cosmologyFunctions_
  type            (cosmologicalMassVarianceFilteredPower                         )              :: cosmologicalMassVariance_
  type            (linearGrowthCollisionlessMatter                               )              :: linearGrowth_
  type            (darkMatterHaloScaleVirialDensityContrastDefinition            )              :: darkMatterHaloScale_
  type            (powerSpectrumWindowFunctionTopHat                             )              :: powerSpectrumWindowFunction_
  type            (powerSpectrumPrimordialPowerLaw                               )              :: powerSpectrumPrimordial_
  type            (transferFunctionCAMB                                          )              :: transferFunctionCDM_
  type            (transferFunctionBode2001                                      )              :: transferFunctionWDM_
  type            (powerSpectrumPrimordialTransferredSimple                      )              :: powerSpectrumPrimordialTransferred_
  type            (powerSpectrumStandard                                         )              :: powerSpectrum_
  type            (darkMatterParticleCDM                                         )              :: darkMatterParticleCDM_
  type            (darkMatterParticleWDMThermal                                  )              :: darkMatterParticleWDM_
  type            (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt)              :: virialDensityContrast_
  type            (darkMatterProfileDMOCuspNFW                                   )              :: darkMatterProfileDMO_
  type            (nodeOperatorDarkMatterProfilePromptCusps                      )              :: nodeOperator_
  type            (nodePropertyExtractorPromptCusps                              )              :: nodePropertyExtractor_
  type            (varying_string                                                )              :: parameterFile
  type            (inputParameters                                               )              :: parameters
  type            (coordinateSpherical                                           )              :: coordinates
  double precision                                                               , parameter    :: massHalo                           =3.15057d+8, redshiftHalo            =1.50000d+00, &
       &                                                                                           radiusSlopeMinus2                  =1.08798d-3, densityVirial200Critical=1.38923d+14
  double precision                                                               , dimension(3) :: propertiesPromptCusp
  double precision                                                                              :: radiusScale                                   , densityScale                        , &
       &                                                                                           radiusVirial                                  , massVirial                          , &
       &                                                                                           densitySlopeLogarithmic                       , densityVirialGalacticus             , &
       &                                                                                           amplitudeCusp                                 , massCusp                            , &
       &                                                                                           yCusp
  character       (len=64                                                        )              :: message

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Prompt cusps")
  ! Read in controlling parameters.
  parameterFile='testSuite/parameters/promptCusps.xml'
  parameters=inputParameters(parameterFile)
  call eventsHooksInitialize            (          )
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  ! Get required objects.
  !![
  <referenceConstruct object="cosmologyParameters_"               >
   <constructor>
    cosmologyParametersSimple                                     (                                                                                                              &amp;
     &amp;                                                         OmegaMatter                               = 0.30890d0                                                       , &amp;
     &amp;                                                         OmegaBaryon                               = 0.04886d0                                                       , &amp;
     &amp;                                                         OmegaDarkEnergy                           = 0.69110d0                                                       , &amp;
     &amp;                                                         temperatureCMB                            = 2.70000d0                                                       , &amp;
     &amp;                                                         HubbleConstant                            =67.36000d0                                                         &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctions_"                >
   <constructor>
    cosmologyFunctionsMatterLambda                                (                                                                                                              &amp;
     &amp;                                                         cosmologyParameters_                      =cosmologyParameters_                                               &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterParticleCDM_"             >
   <constructor>
    darkMatterParticleCDM                                         (                                                                                                              &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterParticleWDM_"             >
   <constructor>
     darkMatterParticleWDMThermal                                 (                                                                                                              &amp;
     &amp;                                                         mass                                      =10.0d0                                                           , &amp;
     &amp;                                                         degreesOfFreedomEffective                 = 1.5d0                                                           , &amp;
     &amp;                                                         cosmologyParameters_                      =cosmologyParameters_                                               &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="linearGrowth_"                      >
   <constructor>
    linearGrowthCollisionlessMatter                               (                                                                                                              &amp;
     &amp;                                                         cosmologyParameters_                      =cosmologyParameters_                                             , &amp;
     &amp;                                                         cosmologyFunctions_                       =cosmologyFunctions_                                                &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumPrimordial_"           >
   <constructor>
    powerSpectrumPrimordialPowerLaw                               (                                                                                                              &amp;
     &amp;                                                         index_                                    =+0.9649d0                                                        , &amp;
     &amp;                                                         running                                   =+0.0000d0                                                        , &amp;
     &amp;                                                         runningRunning                            =+0.0000d0                                                        , &amp;
     &amp;                                                         wavenumberReference                       =+1.0000d0                                                        , &amp;
     &amp;                                                         runningSmallScalesOnly                    =.false.                                                            &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="transferFunctionCDM_"               >
   <constructor>
    transferFunctionCAMB                                          (                                                                                                              &amp;
     &amp;                                                         redshift                                  =100.0d0                                                          , &amp;
     &amp;                                                         cambCountPerDecade                        =  0                                                              , &amp;
     &amp;                                                         darkMatterParticle_                       =darkMatterParticleCDM_                                           , &amp;
     &amp;                                                         cosmologyParameters_                      =cosmologyParameters_                                             , &amp;
     &amp;                                                         cosmologyFunctions_                       =cosmologyFunctions_                                                &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="transferFunctionWDM_"               >
   <constructor>
    transferFunctionBode2001                                      (                                                                                                              &amp;
     &amp;                                                         scaleCutOffModel                          =scaleCutOffModelVogel23SpinHalf                                  , &amp;
     &amp;                                                         epsilon                                   =1.000d0                                                          , &amp;
     &amp;                                                         eta                                       =5.000d0                                                          , &amp;
     &amp;                                                         nu                                        =1.049d0                                                          , &amp;
     &amp;                                                         time                                      =cosmologyFunctions_%equalityEpochMatterRadiation(requestTypeTime), &amp;
     &amp;                                                         transferFunctionCDM                       =transferFunctionCDM_                                             , &amp;
     &amp;                                                         darkMatterParticle_                       =darkMatterParticleWDM_                                           , &amp;
     &amp;                                                         cosmologyParameters_                      =cosmologyParameters_                                             , &amp;
     &amp;                                                         cosmologyFunctions_                       =cosmologyFunctions_                                                &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumPrimordialTransferred_">
   <constructor>
    powerSpectrumPrimordialTransferredSimple                      (                                                                                                              &amp;
     &amp;                                                         powerSpectrumPrimordial_                  =powerSpectrumPrimordial_                                         , &amp;
     &amp;                                                         transferFunction_                         =transferFunctionWDM_                                             , &amp;
     &amp;                                                         linearGrowth_                             =linearGrowth_                                                      &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumWindowFunction_"       >
   <constructor>
    powerSpectrumWindowFunctionTopHat                             (                                                                                                              &amp;
     &amp;                                                         cosmologyParameters_                      =cosmologyParameters_                                               &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologicalMassVariance_"          >
   <constructor>
    cosmologicalMassVarianceFilteredPower                         (                                                                                                              &amp;
     &amp;                                                         sigma8                                    =0.8086539d+0                                                     , &amp;
     &amp;                                                         tolerance                                 =1.0000000d-4                                                     , &amp;
     &amp;                                                         toleranceTopHat                           =1.0000000d-4                                                     , &amp;
     &amp;                                                         nonMonotonicIsFatal                       =.true.                                                           , &amp;
     &amp;                                                         monotonicInterpolation                    =.false.                                                          , &amp;
     &amp;                                                         truncateAtParticleHorizon                 =.false.                                                          , &amp;
     &amp;                                                         cosmologyParameters_                      =cosmologyParameters_                                             , &amp;
     &amp;                                                         cosmologyFunctions_                       =cosmologyFunctions_                                              , &amp;
     &amp;                                                         linearGrowth_                             =linearGrowth_                                                    , &amp;
     &amp;                                                         powerSpectrumPrimordialTransferred_       =powerSpectrumPrimordialTransferred_                              , &amp;
     &amp;                                                         powerSpectrumWindowFunction_              =powerSpectrumWindowFunction_                                       &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrum_"                     >
   <constructor>
    powerSpectrumStandard                                         (                                                                                                              &amp;
     &amp;                                                         cosmologicalMassVariance_                 =cosmologicalMassVariance_                                        , &amp;
     &amp;                                                         powerSpectrumPrimordialTransferred_       =powerSpectrumPrimordialTransferred_                                &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="virialDensityContrast_"             >
   <constructor>
    virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt(                                                                                                              &amp;
     &amp;                                                         cosmologyFunctions_                       =cosmologyFunctions_                                              , &amp;
     &amp;                                                         tableStore                                =.true.                                                             &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterHaloScale_"               >
   <constructor>
    darkMatterHaloScaleVirialDensityContrastDefinition            (                                                                                                              &amp;
     &amp;                                                         cosmologyParameters_                      =cosmologyParameters_                                             , &amp;
     &amp;                                                         cosmologyFunctions_                       =cosmologyFunctions_                                              , &amp;
     &amp;                                                         virialDensityContrast_                    =virialDensityContrast_                                             &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfileDMO_"              >
   <constructor>
    darkMatterProfileDMOCuspNFW                                   (                                                                                                              &amp;
     &amp;                                                         velocityDispersionUseSeriesExpansion      =.true.                                                           , &amp;
     &amp;                                                         toleranceRelativeVelocityDispersion       =1.0d-6                                                           , &amp;
     &amp;                                                         toleranceRelativeVelocityDispersionMaximum=1.0d-3                                                           , &amp;
     &amp;                                                         darkMatterHaloScale_                      =darkMatterHaloScale_                                               &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="nodeOperator_"                      >
   <constructor>
    nodeOperatorDarkMatterProfilePromptCusps                      (                                                                                                              &amp;
     &amp;                                                         nonConvergenceIsFatal                     =.true.                                                           , &amp;
     &amp;                                                         alpha                                     =24.0d0                                                           , &amp;
     &amp;                                                         beta                                      = 7.3d0                                                           , &amp;
     &amp;                                                         kappa                                     = 4.5d0                                                           , &amp;
     &amp;                                                         C                                         = 0.8d0                                                           , &amp;
     &amp;                                                         p                                         = 1.9d0                                                           , &amp;
     &amp;                                                         coefficientScatter                        = 0.0d0                                                           , &amp;
     &amp;                                                         linearGrowth_                             =linearGrowth_                                                    , &amp;
     &amp;                                                         powerSpectrum_                            =powerSpectrum_                                                   , &amp;
     &amp;                                                         cosmologyParameters_                      =cosmologyParameters_                                             , &amp;
     &amp;                                                         cosmologyFunctions_                       =cosmologyFunctions_                                              , &amp;
     &amp;                                                         virialDensityContrast_                    =virialDensityContrast_                                           , &amp;
     &amp;                                                         darkMatterHaloScale_                      =darkMatterHaloScale_                                             , &amp;
     &amp;                                                         darkMatterProfileDMO_                     =darkMatterProfileDMO_                                              &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="nodePropertyExtractor_"             >
   <constructor>
    nodePropertyExtractorPromptCusps                              (                                                                                                              &amp;
     &amp;                                                         darkMatterHaloScale_                      =darkMatterHaloScale_                                               &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  !!]
  ! Test calculation of power spectrum integrals, σ₀ and σ₂.
  call Assert("σ₀(z=0)",nodeOperator_%sigma(0,cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(0.0d0))),7.58104d0,relTol=6.0d-3)
  call Assert("σ₀(z=8)",nodeOperator_%sigma(0,cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(8.0d0))),1.07711d0,relTol=6.0d-3)
  call Assert("σ₂(z=0)",nodeOperator_%sigma(2,cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(0.0d0))),2.57244d4,relTol=5.0d-2)
  call Assert("σ₂(z=8)",nodeOperator_%sigma(2,cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(8.0d0))),3.00944d5,relTol=5.0d-2)
  ! Create a node and set properties.
  node              => treeNode                  (                 )
  basic             => node    %basic            (autoCreate=.true.)
  darkMatterProfile => node    %darkMatterProfile(autoCreate=.true.)
  call basic            %massSet            (                                                                               massHalo         )
  call basic            %timeSet            (cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftHalo     )))
  call basic            %timeLastIsolatedSet(cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftHalo     )))
  call darkMatterProfile%scaleSet           (                                                                               radiusSlopeMinus2)
  ! Report on the virial density assumed by Galacticus.
  densityVirialGalacticus=darkMatterHaloScale_%densityMean(node)
  write (message,'(a21,e11.5,a12)') 'ρᵥ {Galacticus} = ',densityVirialGalacticus,' M☉/Mpc³'
  call displayMessage(trim(message))
  ! Compute, extract, and test prompt cusp properties of this node.
  call nodeOperator_%nodeTreeInitialize(node)
  massDistribution_       =>  darkMatterProfileDMO_ %get    (node             )
  propertiesPromptCusp    =   nodePropertyExtractor_%extract(node,basic%time())
  amplitudeCusp           =   propertiesPromptCusp(1)
  massCusp                =   propertiesPromptCusp(2)
  yCusp                   =   propertiesPromptCusp(3)
  radiusScale             =   darkMatterProfile     %scale  (                 )
  coordinates             =  [radiusScale      ,0.0d0,0.0d0]
  densityScale            =  +massDistribution_%density               (coordinates                                ) &
       &                     *4.0d0                                                                                 &
       &                     /sqrt(                                                                                 &
       &                           +1.0d0                                                                           &
       &                           +yCusp**2                                                                        &
       &                          )
  coordinates             =  [radiusSlopeMinus2,0.0d0,0.0d0]
  densitySlopeLogarithmic =  +massDistribution_%densityGradientRadial (coordinates             ,logarithmic=.true.)
  radiusVirial            =  +massDistribution_%radiusEnclosingDensity(densityVirial200Critical                   )
  massvirial              =  +massDistribution_%massEnclosedBySphere  (radiusVirial                               )
  !![
  <objectDestructor name="massDistribution_"/>
  !!]
  call Assert("A"  ,amplitudeCusp         ,+1.13266d+11,relTol=2.5d-3)
  call Assert("m"  ,massCusp              ,+1.56615d+06,relTol=3.0d-2)
  call Assert("y"  ,yCusp                 ,+2.25290d-01,relTol=4.0d-3)
  call Assert("rₛ" ,radiusScale            ,+1.14456d-03,relTol=1.0d-3)
  call Assert("ρₛ" ,densityScale           ,+1.29838d+16,relTol=1.0d-3)
  call Assert("α₋₂",densitySlopeLogarithmic,-2.00000d+00,relTol=1.0d-3)
  call Assert("rᵥ" ,radiusVirial           ,+8.01838d-03,relTol=1.0d-3)
  call Assert("mᵥ" ,massVirial             ,+3.00000d+08,relTol=1.0d-3)
  ! End unit tests.
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
end program Test_Prompt_Cusps
