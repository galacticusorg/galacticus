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
Contains a program to test calculations for fuzzy dark matter halo density profiles.
!!}

program Test_Dark_Matter_Profiles_Fuzzy_Dark_Matter
  !!{
  Test calculations for fuzzy dark matter halo density profiles.
  !!}
  use, intrinsic :: ISO_C_Binding             , only : c_long
  use            :: Calculations_Resets       , only : Calculations_Reset
  use            :: Cosmology_Functions       , only : cosmologyFunctionsMatterLambda
  use            :: Cosmology_Parameters      , only : cosmologyParametersSimple
  use            :: Coordinates               , only : coordinateSpherical                                           , assignment(=)
  use            :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScaleVirialDensityContrastDefinition
  use            :: Dark_Matter_Particles     , only : darkMatterParticleFuzzyDarkMatter
  use            :: Virial_Density_Contrast   , only : virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt
  use            :: Statistics_Distributions  , only : distributionFunction1DNormal
  use            :: Dark_Matter_Profiles_DMO  , only : darkMatterProfileDMOSolitonNFW
  use            :: Mass_Distributions        , only : kinematicsDistributionSolitonNFW                              , massDistributionSolitonNFW
  use            :: Display                   , only : displayMessage                                                , displayVerbositySet              , &
    &                                                  verbosityLevelStandard
  use            :: Events_Hooks              , only : eventsHooksInitialize, openMPThreadBindingAllLevels
  use            :: Functions_Global_Utilities, only : Functions_Global_Set
  use            :: Galacticus_Nodes          , only : nodeClassHierarchyFinalize                                    , nodeClassHierarchyInitialize     , &
    &                                                  nodeComponentBasic                                            , nodeComponentDarkMatterProfile   , &
    &                                                  treeNode                                                      , mergerTree
  use            :: Input_Parameters          , only : inputParameters
  use            :: Node_Components           , only : Node_Components_Initialize                                    , Node_Components_Thread_Initialize, &
    &                                                  Node_Components_Thread_Uninitialize                           , Node_Components_Uninitialize
  use            :: Nodes_Operators           , only : nodeOperatorDarkMatterProfileSoliton
  use            :: Numerical_Random_Numbers  , only : randomNumberGeneratorGSL
  use            :: Unit_Tests                , only : Assert                                                        , Unit_Tests_Begin_Group           , &
    &                                                  Unit_Tests_End_Group                                          , Unit_Tests_Finish
  use            :: Error                     , only : Error_Report
  use            :: File_Utilities            , only : Count_Lines_in_File
  use            :: Input_Paths               , only : inputPath                                                     , pathTypeExec
  use            :: ISO_Varying_String        , only : varying_string                                                , char                             , &
    &                                                  operator(//)
  implicit none
  type            (treeNode                                                      ), pointer                   :: node_
  type            (mergerTree                                                    ), target                    :: tree
  class           (nodeComponentBasic                                            ), pointer                   :: basic_
  class           (nodeComponentDarkMatterProfile                                ), pointer                   :: darkMatterProfile_
  double precision                                                                , parameter                 :: radiusScale                   =1.0d-2    , massVirial           =6.6d+09, &
       &                                                                                                         massParticle                  =0.8d+0
  double precision                                                                , dimension(:), allocatable :: density                                  , densityTarget                , &
       &                                                                                                         densitySlope                             , densitySlopeNumerical
  type            (nodeOperatorDarkMatterProfileSoliton                          )                            :: nodeOperator_
  type            (darkMatterHaloScaleVirialDensityContrastDefinition            )                            :: darkMatterHaloScale_
  type            (cosmologyParametersSimple                                     )                            :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                                )                            :: cosmologyFunctions_
  type            (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt)                            :: virialDensityContrast_
  type            (darkMatterProfileDMOSolitonNFW                                )                            :: darkMatterProfileChowdhury2021_
  type            (massDistributionSolitonNFW                                    )                            :: massDistribution_
  type            (kinematicsDistributionSolitonNFW                              )                            :: kinematicsDistribution_
  type            (darkMatterParticleFuzzyDarkMatter                             )                            :: darkMatterParticle_
  type            (inputParameters                                               )                            :: parameters
  double precision                                                                                            :: radiusCoreChowdhury            =0.720d-03, radiusScaleChowdhury =1.0d-02, &
       &                                                                                                         radiusSolitonChowdhury         =1.944d-03, densityNFWChowdhury  =5.5d+14, &
       &                                                                                                         densityCoreChowdhury           =1.040d+17
  integer                                                                                                     :: i                                        , fileUnit                     , &
       &                                                                                                         totalLinesInFile                         , dataLinesInFile
  double precision                                                                                            :: potentialNumerical                       , potential                    , &
       &                                                                                                         radiusVirial_                            , radiusScale_                 , &
       &                                                                                                         radiusCore_                              , radiusSoliton_               , &
       &                                                                                                         densityCore_                             , densityScale_                , &
       &                                                                                                         massCore_                                , radius
  type            (coordinateSpherical                                           )                            :: coordinates                              , coordinatesReference
  type            (varying_string                                                )                            :: fileName
  call displayVerbositySet(verbosityLevelStandard)
  call Unit_Tests_Begin_Group("Chowdhury2021 dark matter profiles")
  parameters=inputParameters('testSuite/parameters/darkMatterProfilesChowdhury2021.xml')
  call eventsHooksInitialize()
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  !![
  <referenceConstruct object="cosmologyParameters_"           >
    <constructor>
      cosmologyParametersSimple                                     (                                                                    &amp;
      &amp;                                                          OmegaMatter                                = 0.31530d0            , &amp;
      &amp;                                                          OmegaBaryon                                = 0.04930d0            , &amp;
      &amp;                                                          OmegaDarkEnergy                            = 0.68470d0            , &amp;
      &amp;                                                          temperatureCMB                             = 2.72548d0            , &amp;
      &amp;                                                          HubbleConstant                             =67.36d0                 &amp;
      &amp;                                                         )
    </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctions_"            >
    <constructor>
      cosmologyFunctionsMatterLambda                                (                                                                    &amp;
      &amp;                                                          cosmologyParameters_                       =cosmologyParameters_    &amp;
      &amp;                                                         )
    </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterParticle_"            >
    <constructor>
      darkMatterParticleFuzzyDarkMatter                             (                                                                    &amp;
      &amp;                                                          mass                                       =massParticle          , &amp;
      &amp;                                                          densityFraction                            =1.0d0                   &amp;
      &amp;                                                         )
    </constructor>
  </referenceConstruct>
  <referenceConstruct object="virialDensityContrast_"         >
    <constructor>
      virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt(                                                                    &amp;
      &amp;                                                          tableStore                                 =.true.                , &amp;
      &amp;                                                          cosmologyFunctions_                        =cosmologyFunctions_     &amp;
      &amp;                                                         )
    </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterHaloScale_"           >
    <constructor>
      darkMatterHaloScaleVirialDensityContrastDefinition            (                                                                    &amp;
      &amp;                                                          cosmologyParameters_                       =cosmologyParameters_  , &amp;
      &amp;                                                          cosmologyFunctions_                        =cosmologyFunctions_   , &amp;
      &amp;                                                          virialDensityContrast_                     =virialDensityContrast_  &amp;
      &amp;                                                         )
    </constructor>
  </referenceConstruct>
  <referenceConstruct object="kinematicsDistribution_"        >
    <constructor>
      kinematicsDistributionSolitonNFW                              (                                                                    &amp;
      &amp;                                                          toleranceRelativeVelocityDispersion        =1.0d-6                , &amp;
      &amp;                                                          toleranceRelativeVelocityDispersionMaximum =1.0d-3                  &amp;
      &amp;                                                         )
    </constructor>
  </referenceConstruct>
  <referenceConstruct object="massDistribution_"              >
    <constructor>
      massDistributionSolitonNFW                                    (                                                                    &amp;
      &amp;                                                          radiusScale                                =radiusScaleChowdhury  , &amp;
      &amp;                                                          radiusCore                                 =radiusCoreChowdhury   , &amp;
      &amp;                                                          radiusSoliton                              =radiusSolitonChowdhury, &amp;
      &amp;                                                          densitySolitonCentral                      =densityCoreChowdhury  , &amp;
      &amp;                                                          densityNormalizationNFW                    =densityNFWChowdhury     &amp;
      &amp;                                                         )
    </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfileChowdhury2021_">
    <constructor>
      darkMatterProfileDMOSolitonNFW                                (                                                                    &amp;
      &amp;                                                          darkMatterParticle_                        =darkMatterParticle_   , &amp;
      &amp;                                                          darkMatterHaloScale_                       =darkMatterHaloScale_  , &amp;
      &amp;                                                          virialDensityContrast_                     =virialDensityContrast_, &amp;
      &amp;                                                          cosmologyFunctions_                        =cosmologyFunctions_   , &amp;
      &amp;                                                          cosmologyParameters_                       =cosmologyParameters_  , &amp;
      &amp;                                                          toleranceRelativeVelocityDispersion        =1.0d-6                , &amp;
      &amp;                                                          toleranceRelativeVelocityDispersionMaximum =1.0d-3                  &amp;
      &amp;                                                         )
    </constructor>
  </referenceConstruct>
  !!]
  call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
  ! Construct a node with properties matched to that in Chowdhury et al. (2021; https://ui.adsabs.harvard.edu/abs/2021ApJ...916...27D).
  node_              => treeNode                  (                 )
  node_%hostTree     => tree
  basic_             => node_   %basic            (autoCreate=.true.)
  darkMatterProfile_ => node_   %darkMatterProfile(autoCreate=.true.)
  call basic_    %timeSet            (cosmologyFunctions_%cosmicTime(1.0d0))
  call basic_    %timeLastIsolatedSet(cosmologyFunctions_%cosmicTime(1.0d0))
  call basic_    %massSet            (massVirial                           )
  call darkMatterProfile_%scaleSet(radiusScale)
  call Calculations_Reset(node_)
  allocate(randomNumberGeneratorGSL :: tree%randomNumberGenerator_)
  select type (randomNumberGenerator_ => tree%randomNumberGenerator_)
  type is (randomNumberGeneratorGSL)
     randomNumberGenerator_=randomNumberGeneratorGSL(seed_=8322_c_long)
  end select
  call tree%properties%initialize()
  ! Build the node operator to compute massCore.
  nodeOperator_=nodeOperatorDarkMatterProfileSoliton(                                               &
                                                     darkMatterHaloScale_  =darkMatterHaloScale_  , &
                                                     darkMatterParticle_   =darkMatterParticle_   , &
                                                     cosmologyFunctions_   =cosmologyFunctions_   , &
                                                     cosmologyParameters_  =cosmologyParameters_  , &
                                                     virialDensityContrast_=virialDensityContrast_  &
                                                    )
  call nodeOperator_%nodeTreeInitialize(node_)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Chowdhury2021 profile")
  call darkMatterProfileChowdhury2021_%computeProperties(node_,radiusVirial_,radiusScale_,radiusCore_,radiusSoliton_,densityCore_,densityScale_,massCore_)
  ! Test that the potential agrees with a numerical evaluation.
  coordinates         =[2.0d+0*radiusScale,0.0d0,0.0d0]
  coordinatesReference=[1.0d+3*radiusScale,0.0d0,0.0d0]
  potential           =+massDistribution_%potential         (coordinates=coordinates         ) &
       &               -massDistribution_%potential         (coordinates=coordinatesReference)
  potentialNumerical  =+massDistribution_%potentialNumerical(coordinates=coordinates         ) &
       &               -massDistribution_%potentialNumerical(coordinates=coordinatesReference)
  call Assert("Potential",potential,potentialNumerical,relTol=2.0d-3)
  ! Evaluate the density profile.
  fileName=inputPath(pathTypeExec)//'testSuite/data/densityProfileFuzzyDarkMatterChowdhury2021.txt'
  totalLinesInFile=Count_Lines_in_File(fileName    )
  dataLinesInFile =Count_Lines_in_File(fileName,'#')-1
  allocate(density              (dataLinesInFile))
  allocate(densityTarget        (dataLinesInFile))
  allocate(densitySlope         (dataLinesInFile))
  allocate(densitySlopeNumerical(dataLinesInFile))
  open(newunit=fileUnit,file=char(fileName),status='old',form='formatted')
  do i=1,(totalLinesInFile-dataLinesInFile)
     read (fileUnit,*)
  end do
  do i=1,dataLinesInFile
     read (fileUnit,*) radius,densityTarget(i)
     coordinates = [radius,0.0d0,0.0d0]
     density              (i)=massDistribution_%density                       (coordinates                   )
     densitySlope         (i)=massDistribution_%densityGradientRadial         (coordinates,logarithmic=.true.)
     densitySlopeNumerical(i)=massDistribution_%densityGradientRadialNumerical(coordinates,logarithmic=.true.)
     ! Very close to the soliton radius the numerical estimate of density slope is unreliable (due to the discontinuity at that
     ! radius). Patch over this by replacing the numerical esimate with the analytic result close to the soliton radius.
     if (abs(log(radius/radiusSolitonChowdhury)) < 0.015d0) densitySlopeNumerical(i)=densitySlope(i)
  end do
  close(fileUnit)
  call Assert("Density"         ,density     ,densityTarget        ,relTol=5.0d-1)
  call Assert("Density gradient",densitySlope,densitySlopeNumerical,relTol=1.0d-3)
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
end program Test_Dark_Matter_Profiles_Fuzzy_Dark_Matter
