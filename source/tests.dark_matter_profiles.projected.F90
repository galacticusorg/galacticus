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
Contains a program to test calculations of projected dark matter profiles.
!!}

program Test_Dark_Matter_Profiles_Projected
  !!{
  Test calculations of projected dark matter profiles.
  !!}
  use, intrinsic :: ISO_C_Binding             , only : c_size_t
  use            :: Calculations_Resets       , only : Calculations_Reset
  use            :: Cosmology_Functions       , only : cosmologyFunctionsMatterLambda
  use            :: Cosmology_Parameters      , only : cosmologyParametersSimple
  use            :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScaleVirialDensityContrastDefinition
  use            :: Virial_Density_Contrast   , only : virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt
  use            :: Dark_Matter_Profiles_DMO  , only : darkMatterProfileDMONFW
  use            :: Node_Property_Extractors  , only : nodePropertyExtractorProjectedMass                            , nodePropertyExtractorProjectedDensity
  use            :: Display                   , only : displayMessage                                                , displayVerbositySet                  , verbosityLevelStandard
  use            :: Events_Hooks              , only : eventsHooksInitialize
  use            :: Functions_Global_Utilities, only : Functions_Global_Set
  use            :: Galacticus_Nodes          , only : nodeClassHierarchyFinalize                                    , nodeClassHierarchyInitialize         , nodeComponentBasic                 , nodeComponentDarkMatterProfile, &
       &                                               treeNode
  use            :: Input_Parameters          , only : inputParameters
  use            :: Node_Components           , only : Node_Components_Initialize                                    , Node_Components_Thread_Initialize    , Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use            :: Unit_Tests                , only : Assert                                                        , Skip                                 , Unit_Tests_Begin_Group             , Unit_Tests_End_Group          , &
       &                                               Unit_Tests_Finish
  use            :: Multi_Counters            , only : multiCounter
  use            :: ISO_Varying_String        , only : varying_string                                                , assignment(=)                        , var_str
  implicit none
  type            (darkMatterHaloScaleVirialDensityContrastDefinition            )                          :: darkMatterHaloScale_
  type            (cosmologyParametersSimple                                     )                          :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                                )                          :: cosmologyFunctions_
  type            (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt)                          :: virialDensityContrast_
  type            (darkMatterProfileDMONFW                                       )                          :: darkMatterProfileDMO_
  type            (nodePropertyExtractorProjectedMass                            )                          :: nodePropertyExtractorProjectedMass_
  type            (nodePropertyExtractorProjectedDensity                         )                          :: nodePropertyExtractorProjectedDensity_
  type            (treeNode                                                      ), pointer                 :: node__
  class           (nodeComponentBasic                                            ), pointer                 :: basic__
  class           (nodeComponentDarkMatterProfile                                ), pointer                 :: darkMatterProfile__
  double precision                                                                , parameter               :: concentration                         =+8.0d0, massVirial   =+1.0d12
  integer                                                                         , parameter               :: countRadii=3
  type            (varying_string                                                ), dimension(countRadii  ) :: radiusSpecifiers
  double precision                                                                , dimension(countRadii,1) :: densityProjected                             , massProjected
  type            (inputParameters                                               )                          :: parameters
  double precision                                                                                          :: radiusScale                                  , radiusVirial
  type            (multiCounter                                                  )                          :: instance
  
  call displayVerbositySet(verbosityLevelStandard)
  call Unit_Tests_Begin_Group("Projected dark matter profiles")
  parameters=inputParameters(var_str('testSuite/parameters/darkMatterProfilesProjected.xml'))
  call eventsHooksInitialize()
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  !![
  <referenceConstruct object="cosmologyParameters_"  >
   <constructor>
    cosmologyParametersSimple                                     (                                               &amp;
     &amp;                                                         OmegaMatter           = 0.30d0               , &amp;
     &amp;                                                         OmegaBaryon           = 0.00d0               , &amp;
     &amp;                                                         OmegaDarkEnergy       = 0.70d0               , &amp;
     &amp;                                                         temperatureCMB        = 2.78d0               , &amp;
     &amp;                                                         HubbleConstant        =70.00d0                 &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctions_"   >
   <constructor>
    cosmologyFunctionsMatterLambda                                (                                               &amp;
     &amp;                                                         cosmologyParameters_  =cosmologyParameters_    &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="virialDensityContrast_">
   <constructor>
    virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt(                                               &amp;
     &amp;                                                         tableStore            =.true.,                 &amp;
     &amp;                                                         cosmologyFunctions_   =cosmologyFunctions_     &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterHaloScale_"  >
   <constructor>
    darkMatterHaloScaleVirialDensityContrastDefinition            (                                               &amp;
     &amp;                                                         cosmologyParameters_  =cosmologyParameters_  , &amp;
     &amp;                                                         cosmologyFunctions_   =cosmologyFunctions_   , &amp;
     &amp;                                                         virialDensityContrast_=virialDensityContrast_  &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  !!]
  node__                            =>  treeNode                              (                 )
  basic__                           =>  node__              %basic            (autoCreate=.true.)
  darkMatterProfile__               =>  node__              %darkMatterProfile(autoCreate=.true.)
  call basic__    %timeSet            (cosmologyFunctions_%cosmicTime(1.0d0))
  call basic__    %timeLastIsolatedSet(cosmologyFunctions_%cosmicTime(1.0d0))
  call basic__    %massSet            (massVirial                           )
  radiusVirial                      =  +darkMatterHaloScale_%radiusVirial     (           node__)
  radiusScale                       =  +radiusVirial &
       &                               /concentration
  call darkMatterProfile__%scaleSet(radiusScale)
  call Calculations_Reset(node__)
  radiusSpecifiers(1)='darkMatterScaleRadius:all:all:0.5'
  radiusSpecifiers(2)='darkMatterScaleRadius:all:all:1.0'
  radiusSpecifiers(3)='darkMatterScaleRadius:all:all:2.0'
  !![
  <referenceConstruct object="darkMatterProfileDMO_"               >
   <constructor>
    darkMatterProfileDMONFW              (                                                            &amp;
     &amp;                                velocityDispersionUseSeriesExpansion=.false.              , &amp;
     &amp;                                darkMatterHaloScale_                =darkMatterHaloScale_   &amp;
     &amp;                               )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="nodePropertyExtractorProjectedMass_"  >
   <constructor>
    nodePropertyExtractorProjectedMass   (                                                            &amp;
     &amp;                                radiusSpecifiers                    =radiusSpecifiers     , &amp;
     &amp;                                includeRadii                        =.false.              , &amp;
     &amp;                                darkMatterHaloScale_                =darkMatterHaloScale_   &amp;
     &amp;                               )
   </constructor>
  </referenceConstruct>
 <referenceConstruct object="nodePropertyExtractorProjectedDensity_">
   <constructor>
    nodePropertyExtractorProjectedDensity(                                                            &amp;
     &amp;                                radiusSpecifiers                    =radiusSpecifiers     , &amp;
     &amp;                                includeRadii                        =.false.              , &amp;
     &amp;                                tolerateIntegrationFailures         =.false.              , &amp;
     &amp;                                darkMatterHaloScale_                =darkMatterHaloScale_   &amp;
     &amp;                               )
   </constructor>
  </referenceConstruct>
  !!]
  ! Begin tests.
  instance        =multiCounter([1_c_size_t])
  densityProjected=nodePropertyExtractorProjectedDensity_%extract(node__,basic__%time(),instance)
  massProjected   =nodePropertyExtractorProjectedMass_   %extract(node__,basic__%time(),instance)  
  call Assert("Projected density",densityProjected(:,1),[8.09887d13,3.89131d13,1.53862d13],relTol=2.0d-3)
  call Assert("Projected mass"   ,   massProjected(:,1),[1.02724d11,2.34537d11,4.62114d11],relTol=2.0d-3)
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
end program Test_Dark_Matter_Profiles_Projected
