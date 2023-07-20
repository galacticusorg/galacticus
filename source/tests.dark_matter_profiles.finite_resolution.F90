!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
Contains a program to test calculations for finite resolution dark matter profiles.
!!}

program Test_Dark_Matter_Profiles_Finite_Resolution
  !!{
  Test calculations for finite resolution dark matter profiles.
  !!}
  use :: Calculations_Resets       , only : Calculations_Reset
  use :: Cosmology_Functions       , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters      , only : cosmologyParametersSimple
  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScaleVirialDensityContrastDefinition
  use :: Virial_Density_Contrast   , only : virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt
  use :: Dark_Matter_Profiles_DMO  , only : darkMatterProfileDMOFiniteResolution                          , darkMatterProfileDMOFiniteResolutionNFW, darkMatterProfileDMONFW
  use :: Mass_Distributions        , only : nonAnalyticSolversNumerical
  use :: Display                   , only : displayMessage                                                , displayVerbositySet                    , verbosityLevelStandard
  use :: Events_Hooks              , only : eventsHooksInitialize
  use :: Functions_Global_Utilities, only : Functions_Global_Set
  use :: Galacticus_Nodes          , only : nodeClassHierarchyFinalize                                    , nodeClassHierarchyInitialize           , nodeComponentBasic                 , nodeComponentDarkMatterProfile, &
          &                                 treeNode                                                      , nodeComponentSatellite
  use :: Input_Parameters          , only : inputParameters
  use :: Node_Components           , only : Node_Components_Initialize                                    , Node_Components_Thread_Initialize      , Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Numerical_Ranges          , only : Make_Range                                                    , rangeTypeLogarithmic
  use :: Numerical_Constants_Math  , only : Pi
  use :: Unit_Tests                , only : Assert                                                        , Unit_Tests_Begin_Group                 , Unit_Tests_End_Group               , Unit_Tests_Finish
  implicit none
  type            (darkMatterHaloScaleVirialDensityContrastDefinition            )                        :: darkMatterHaloScale_
  type            (cosmologyParametersSimple                                     )                        :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                                )                        :: cosmologyFunctions_
  type            (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt)                        :: virialDensityContrast_
  type            (darkMatterProfileDMONFW                                       )                        :: darkMatterProfileNFW_
  type            (darkMatterProfileDMOFiniteResolution                          )                        :: darkMatterProfileFiniteResolution_
  type            (darkMatterProfileDMOFiniteResolutionNFW                       )                        :: darkMatterProfileFiniteResolutionNFW_
  type            (treeNode                                                      ), pointer               :: node_
  class           (nodeComponentBasic                                            ), pointer               :: basic_
  class           (nodeComponentDarkMatterProfile                                ), pointer               :: darkMatterProfile_
  type            (inputParameters                                               )                        :: parameters
  double precision                                                                , parameter             :: concentration                        =+8.0d0, massVirial        =+1.0d12
  integer                                                                         , parameter             :: countRadii                           =10
  double precision                                                                , parameter             :: radiiMinimum                         =1.0d-6, radiiMaximum      =+1.0d01
  double precision                                                                , dimension(countRadii) :: massNumerical                               , mass                      , &
       &                                                                                                     densityNumerical                            , density                   , &
       &                                                                                                     velocityDispersionNumerical                 , velocityDispersion        , &
       &                                                                                                     radii                                       , radiusEnclosingDensity    , &
       &                                                                                                     radiusEnclosingMass
  integer                                                                                                 :: i
  double precision                                                                                        :: radiusScale                                 , radiusVirial              , &
       &                                                                                                     potentialNumerical                          , potential                 , &
       &                                                                                                     energyNumerical                             , energy

  call displayVerbositySet(verbosityLevelStandard)
  call Unit_Tests_Begin_Group("Finite resolution dark matter profiles")
  parameters=inputParameters('testSuite/parameters/darkMatterProfilesFiniteResolution.xml')
  call eventsHooksInitialize()
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  node_              =>  treeNode                  (                 )
  basic_             =>  node_   %basic            (autoCreate=.true.)
  darkMatterProfile_ =>  node_   %darkMatterProfile(autoCreate=.true.)
  !![
  <referenceConstruct object="cosmologyParameters_"                 >
   <constructor>
    cosmologyParametersSimple                                     (                                                                 &amp;
     &amp;                                                         OmegaMatter                         = 0.30d0                   , &amp;
     &amp;                                                         OmegaBaryon                         = 0.00d0                   , &amp;
     &amp;                                                         OmegaDarkEnergy                     = 0.70d0                   , &amp;
     &amp;                                                         temperatureCMB                      = 2.78d0                   , &amp;
     &amp;                                                         HubbleConstant                      =70.00d0                     &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctions_"                  >
   <constructor>
    cosmologyFunctionsMatterLambda                                (                                                                 &amp;
     &amp;                                                         cosmologyParameters_                =cosmologyParameters_        &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="virialDensityContrast_"               >
   <constructor>
    virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt(                                                                 &amp;
     &amp;                                                         tableStore                          =.true.,                     &amp;
     &amp;                                                         cosmologyFunctions_                 =cosmologyFunctions_         &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterHaloScale_"                 >
   <constructor>
    darkMatterHaloScaleVirialDensityContrastDefinition            (                                                                 &amp;
     &amp;                                                         cosmologyParameters_                =cosmologyParameters_      , &amp;
     &amp;                                                         cosmologyFunctions_                 =cosmologyFunctions_       , &amp;
     &amp;                                                         virialDensityContrast_              =virialDensityContrast_      &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfileNFW_"                >
   <constructor>
    darkMatterProfileDMONFW                                       (                                                                 &amp;
     &amp;                                                         velocityDispersionUseSeriesExpansion=.false.                   , &amp;
     &amp;                                                         darkMatterHaloScale_                =darkMatterHaloScale_        &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfileFiniteResolution_"   >
   <constructor>
    darkMatterProfileDMOFiniteResolution                         (                                                                   &amp;
     &amp;                                                         lengthResolution                    =1.0d-2                     , &amp;
     &amp;                                                         massResolution                      =1.0d+6                     , &amp;
     &amp;                                                         resolutionIsComoving                =.false.                    , &amp;
     &amp;                                                         nonAnalyticSolver                   =nonAnalyticSolversNumerical, &amp;
     &amp;                                                         darkMatterProfileDMO_               =darkMatterProfileNFW_      , &amp;
     &amp;                                                         darkMatterHaloScale_                =darkMatterHaloScale_       , &amp;
     &amp;                                                         cosmologyFunctions_                 =cosmologyFunctions_          &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfileFiniteResolutionNFW_">
   <constructor>
    darkMatterProfileDMOFiniteResolutionNFW                      (                                                                   &amp;
     &amp;                                                         lengthResolution                    =1.0d-2                     , &amp;
     &amp;                                                         massResolution                      =1.0d+6                     , &amp;
     &amp;                                                         resolutionIsComoving                =.false.                    , &amp;
     &amp;                                                         nonAnalyticSolver                   =nonAnalyticSolversNumerical, &amp;
     &amp;                                                         darkMatterHaloScale_                =darkMatterHaloScale_       , &amp;
     &amp;                                                         cosmologyFunctions_                 =cosmologyFunctions_          &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  !!]
  call basic_    %timeSet            (cosmologyFunctions_%cosmicTime(1.0d0))
  call basic_    %timeLastIsolatedSet(cosmologyFunctions_%cosmicTime(1.0d0))
  call basic_    %massSet            (massVirial                           )
  radiusVirial =+darkMatterHaloScale_%radiusVirial(node_)
  radiusScale  =+radiusVirial &
       &        /concentration
  call darkMatterProfile_%scaleSet(radiusScale)
  call Calculations_Reset(node_)
  ! Begin tests.
  call Unit_Tests_Begin_Group("Finite resolution NFW profile")
  radii=Make_Range(radiiMinimum,radiiMaximum,countRadii,rangeTypeLogarithmic)*radiusScale
  do i=1,countRadii
     mass                       (i)=+darkMatterProfileFiniteResolutionNFW_%enclosedMass            (node_,radius =                       radii(i)   )
     massNumerical              (i)=+darkMatterProfileFiniteResolution_   %enclosedMass            (node_,radius =                       radii(i)   )
     density                    (i)=+darkMatterProfileFiniteResolutionNFW_%density                 (node_,radius =                       radii(i)   )
     densityNumerical           (i)=+darkMatterProfileFiniteResolution_   %density                 (node_,radius =                       radii(i)   )
     velocityDispersion         (i)=+darkMatterProfileFiniteResolutionNFW_%radialVelocityDispersion(node_,radius =                       radii(i)   )
     velocityDispersionNumerical(i)=+darkMatterProfileFiniteResolution_   %radialVelocityDispersion(node_,radius =                       radii(i)   )
     radiusEnclosingDensity     (i)=+darkMatterProfileFiniteResolutionNFW_%radiusEnclosingDensity  (node_,density=3.0d0*mass(i)/4.0d0/Pi/radii(i)**3)
     radiusEnclosingMass        (i)=+darkMatterProfileFiniteResolutionNFW_%radiusEnclosingMass     (node_,mass   =      mass(i)                     )
  end do
  call Assert(                                    &
       &             "Density"                  , &
       &             density                    , &
       &             densityNumerical           , &
       &      relTol=+2.0d-2                      &
       &     )
  call Assert(                                    &
       &             "Enclosed mass"            , &
       &             mass                       , &
       &             massNumerical              , &
       &      relTol=+2.0d-2                      &
       &     )
  call Assert(                                    &
       &             "Velocity dispersion"      , &
       &             velocityDispersion         , &
       &             velocityDispersionNumerical, &
       &      relTol=+2.0d-2                      &
       &     )
  call Assert(                                    &
       &             "Radius enclosing density" , &
       &              radiusEnclosingDensity    , &
       &              radii                     , &
       &      absTol=+1.0d-4                    , &
       &      relTol=+2.0d-2                      &
       &     )  
  call Assert(                                    &
       &             "Radius enclosing mass"    , &
       &              radiusEnclosingMass       , &
       &              radii                     , &
       &      relTol=+2.0d-2                      &
       &     )  
  !! When comparing to the numerical calculation of potential we take the potential relative to 100 times the virial radius, as
  !! that is the radius to which the numerical solution is integrated.
  potential         =+darkMatterProfileFiniteResolutionNFW_%potential(node_,radius=2.0d+0*radiusScale ) &
       &             -darkMatterProfileFiniteResolutionNFW_%potential(node_,radius=1.0d+2*radiusVirial)
  potentialNumerical=+darkMatterProfileFiniteResolution_   %potential(node_,radius=2.0d+0*radiusScale )
  call Assert(                           &
       &             "Potential"       , &
       &             potential         , &
       &             potentialNumerical, &
       &      relTol=+2.0d-2             &
       &     )
  !! Total energy.
  energy         =+darkMatterProfileFiniteResolutionNFW_%energy(node_)
  energyNumerical=+darkMatterProfileFiniteResolution_   %energy(node_)
  call Assert(                        &
       &             "Energy"       , &
       &             energy         , &
       &             energyNumerical, &
       &      relTol=+3.0d-3          &
       &     )
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
end program Test_Dark_Matter_Profiles_Finite_Resolution
