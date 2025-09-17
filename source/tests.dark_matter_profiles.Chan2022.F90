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
Contains a program to test calculations for fuzzy dark matter soliton core mass evolution.
!!}

program Test_Dark_Matter_Profiles_Chan2022
  !!{
  Test the calculation for fuzzy dark matter soliton core mass evolution.
  !!}
  use, intrinsic :: ISO_C_Binding               , only : c_long
  use            :: Calculations_Resets         , only : Calculations_Reset
  use            :: Cosmology_Functions         , only : cosmologyFunctionsMatterLambda
  use            :: Cosmology_Parameters        , only : cosmologyParametersSimple
  use            :: Dark_Matter_Halo_Scales     , only : darkMatterHaloScaleVirialDensityContrastDefinition
  use            :: Dark_Matter_Particles       , only : darkMatterParticleFuzzyDarkMatter
  use            :: Virial_Density_Contrast     , only : virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt
  use            :: Dark_Matter_Profiles_DMO    , only : darkMatterProfileDMOSolitonNFW
  use            :: Display                     , only : displayVerbositySet                                           , verbosityLevelStandard
  use            :: Events_Hooks                , only : eventsHooksInitialize
  use            :: Functions_Global_Utilities  , only : Functions_Global_Set
  use            :: Galacticus_Nodes            , only : nodeClassHierarchyFinalize                                    , nodeClassHierarchyInitialize     , &
    &                                                    nodeComponentBasic                                            , nodeComponentDarkMatterProfile   , &
    &                                                    treeNode                                                      , mergerTree
  use            :: Input_Parameters            , only : inputParameters
  use            :: Node_Components             , only : Node_Components_Initialize                                    , Node_Components_Thread_Initialize, &
    &                                                    Node_Components_Thread_Uninitialize                           , Node_Components_Uninitialize
  use            :: Nodes_Operators             , only : nodeOperatorDarkMatterProfileSoliton
  use            :: Numerical_Random_Numbers    , only : randomNumberGeneratorGSL
  use            :: Numerical_Constants_Prefixes, only : kilo
  use            :: Cosmology_Parameters        , only : hubbleUnitsStandard                                           , hubbleUnitsLittleH
  use            :: Unit_Tests                  , only : Assert                                                        , Unit_Tests_Begin_Group           , &
    &                                                    Unit_Tests_End_Group                                          , Unit_Tests_Finish
  implicit none
  type            (treeNode                                                      ), pointer                   :: node_
  type            (mergerTree                                                    ), target                    :: tree
  class           (nodeComponentBasic                                            ), pointer                   :: basic_
  class           (nodeComponentDarkMatterProfile                                ), pointer                   :: darkMatterProfile_
  type            (nodeOperatorDarkMatterProfileSoliton                          )                            :: nodeOperator_
  type            (darkMatterHaloScaleVirialDensityContrastDefinition            )                            :: darkMatterHaloScale_
  type            (cosmologyParametersSimple                                     )                            :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                                )                            :: cosmologyFunctions_
  type            (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt)                            :: virialDensityContrast_
  type            (darkMatterProfileDMOSolitonNFW                                )                            :: darkMatterProfileChan2022_
  type            (darkMatterParticleFuzzyDarkMatter                             )                            :: darkMatterParticle_
  type            (inputParameters                                               )                            :: parameters
  double precision                                                                , parameter                 :: massVirial        =1.0d+12          , massParticle       =5.0d+0
  double precision                                                                , parameter                 :: alpha             =0.515            , beta               =8.0d6  , &
         &                                                                                                       gamma             =10.0d0**(-5.73d0)                                 ! Best-fitting parameters from Chan et al. (2022; MNRAS; 551; 943; https://ui.adsabs.harvard.edu/abs/2022MNRAS.511..943C).
  double precision                                                                                            :: massHalo                            , expansionFactor            , &
         &                                                                                                       zeta_0                              , zeta_z                     , &
         &                                                                                                       radiusVirial_                       , radiusScale_               , &
         &                                                                                                       radiusCore_                         , radiusSoliton_             , &
         &                                                                                                       densityCore_                        , densityScale_              , &
         &                                                                                                       massCore_                           , massCoreAnalytic
  call displayVerbositySet(verbosityLevelStandard)
  call Unit_Tests_Begin_Group("Chan2022 solitonic core mass test")
  parameters=inputParameters('testSuite/parameters/darkMatterProfilesChan2022.xml')
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
  <referenceConstruct object="darkMatterProfileChan2022_"     >
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
  
  node_              => treeNode                  (                 )
  node_%hostTree     => tree
  basic_             => node_   %basic            (autoCreate=.true.)
  darkMatterProfile_ => node_   %darkMatterProfile(autoCreate=.true.)
  call basic_    %timeSet            (cosmologyFunctions_%cosmicTime(1.0d0))
  call basic_    %timeLastIsolatedSet(cosmologyFunctions_%cosmicTime(1.0d0))
  call basic_    %massSet            (massVirial                           )
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
  call Unit_Tests_Begin_Group("Chan2022 massCore")
  call darkMatterProfileChan2022_%computeProperties(node_,radiusVirial_,radiusScale_,radiusCore_,radiusSoliton_,densityCore_,densityScale_,massCore_)

  ! Use Equation (15) of Chan et al. (2022) to calculate the massCore.
  expansionFactor      =+cosmologyFunctions_   %expansionFactor            (basic_%time          ())
  massHalo             =+basic_                %mass                       (                       )
  zeta_0               =+virialDensityContrast_%densityContrast            (massHalo,expansionFactor=1.0d0          )
  zeta_z               =+virialDensityContrast_%densityContrast            (massHalo,expansionFactor=expansionFactor)
  massCoreAnalytic     =+(                          & ! Equation (15) of Chan et al. (2022; MNRAS; 551; 943; https://ui.adsabs.harvard.edu/abs/2022MNRAS.511..943C).
         &                +beta                     &
         &                *(                        &
         &                  +darkMatterParticle_%mass()*kilo      &
         &                  /8.0d-23                &
         &                 )**(-1.5d0)              &
         &                +(                        &
         &                  +sqrt(                  &
         &                        +zeta_z           &
         &                        /zeta_0           &
         &                       )                  &
         &                  *massHalo               &
         &                  /gamma                  &
         &                 )**alpha                 &
         &                *(                        &
         &                  +darkMatterParticle_%mass()*kilo      &
         &                  /8.0d-23                &
         &                 )**(1.5d0*(alpha-1.0d0)) &
         &               )&
         &              /sqrt(expansionFactor)

  call Assert("Core mass",massCore_,massCoreAnalytic,relTol=1.0d-4*massCoreAnalytic)
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
end program Test_Dark_Matter_Profiles_Chan2022
