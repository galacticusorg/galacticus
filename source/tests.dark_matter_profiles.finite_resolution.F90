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
Contains a program to test calculations for finite resolution dark matter profiles.
!!}

program Test_Dark_Matter_Profiles_Finite_Resolution
  !!{
  Test calculations for finite resolution dark matter profiles.
  !!}
  use :: Calculations_Resets       , only : Calculations_Reset
  use :: Coordinates               , only : coordinateSpherical                                           , assignment(=)
  use :: Cosmology_Functions       , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters      , only : cosmologyParametersSimple
  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScaleVirialDensityContrastDefinition
  use :: Error                     , only : Error_Handler_Register
  use :: Virial_Density_Contrast   , only : virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt
  use :: Dark_Matter_Profiles_DMO  , only : darkMatterProfileDMOFiniteResolution                          , darkMatterProfileDMOFiniteResolutionNFW, darkMatterProfileDMONFW
  use :: Mass_Distributions        , only : nonAnalyticSolversNumerical                                   , massDistributionClass                  , kinematicsDistributionClass
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
  class           (massDistributionClass                                         ), pointer               :: massDistributionFiniteResolutionNFW_             , massDistributionFiniteResolution_              , &
       &                                                                                                     massDistribution_
  class           (kinematicsDistributionClass                                   ), pointer               :: kinematicsDistributionFiniteResolutionNFW_       , kinematicsDistributionFiniteResolution_
  type            (inputParameters                                               )                        :: parameters
  double precision                                                                , parameter             :: concentration                             =+8.0d0, massVirial                             =+1.0d12
  integer                                                                         , parameter             :: countRadii                                =10
  double precision                                                                , parameter             :: radiiMinimum                              =1.0d-6, radiiMaximum                           =+1.0d01
  double precision                                                                , dimension(countRadii) :: massNumerical                                    , mass                                           , &
       &                                                                                                     densityNumerical                                 , density                                        , &
       &                                                                                                     velocityDispersionNumerical                      , velocityDispersion                             , &
       &                                                                                                     radii                                            , radiusEnclosingDensity                         , &
       &                                                                                                     radiusEnclosingMass
  integer                                                                                                 :: i
  double precision                                                                                        :: radiusScale                                      , radiusVirial                                   , &
       &                                                                                                     potentialNumerical                               , potential                                      , &
       &                                                                                                     energyNumerical                                  , energy
  type            (coordinateSpherical                                           )                        :: coordinates                                      , coordinatesScale                               , &
       &                                                                                                     coordinatesVirial

  call displayVerbositySet(verbosityLevelStandard)
  call Unit_Tests_Begin_Group("Finite resolution dark matter profiles")
  parameters=inputParameters('testSuite/parameters/darkMatterProfilesFiniteResolution.xml')
  call eventsHooksInitialize            (          )
  call Functions_Global_Set             (          )
  call Error_Handler_Register           (          )
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
  massDistribution_                          => darkMatterProfileNFW_                %get                   (node_)
  massDistributionFiniteResolutionNFW_       => darkMatterProfileFiniteResolutionNFW_%get                   (node_)
  massDistributionFiniteResolution_          => darkMatterProfileFiniteResolution_   %get                   (node_)
  kinematicsDistributionFiniteResolutionNFW_ => massDistributionFiniteResolutionNFW_ %kinematicsDistribution(     )
  kinematicsDistributionFiniteResolution_    => massDistributionFiniteResolution_    %kinematicsDistribution(     )
  ! Begin tests.
  call Unit_Tests_Begin_Group("Finite resolution NFW profile")
  radii=Make_Range(radiiMinimum,radiiMaximum,countRadii,rangeTypeLogarithmic)*radiusScale
  do i=1,countRadii
     coordinates                   =[radii(i),0.0d0,0.0d0]
     mass                       (i)=+massDistributionFiniteResolutionNFW_      %massEnclosedBySphere    (radius     =                       radii(i)                                                                                                                         )
     massNumerical              (i)=+massDistributionFiniteResolution_         %massEnclosedBySphere    (radius     =                       radii(i)                                                                                                                         )
     density                    (i)=+massDistributionFiniteResolutionNFW_      %density                 (coordinates=                       coordinates                                                                                                                      )
     densityNumerical           (i)=+massDistributionFiniteResolution_         %density                 (coordinates=                       coordinates                                                                                                                      )
     velocityDispersion         (i)=+kinematicsDistributionFiniteResolutionNFW_%velocityDispersion1D    (coordinates=                       coordinates,massDistribution_=massDistributionFiniteResolutionNFW_,massDistributionEmbedding=massDistributionFiniteResolutionNFW_)
     velocityDispersionNumerical(i)=+kinematicsDistributionFiniteResolution_   %velocityDispersion1D    (coordinates=                       coordinates,massDistribution_=massDistributionFiniteResolution_   ,massDistributionEmbedding=massDistributionFiniteResolution_   )
     radiusEnclosingDensity     (i)=+massDistributionFiniteResolutionNFW_      %radiusEnclosingDensity  (density    =3.0d0*mass(i)/4.0d0/Pi/radii(i)**3                                                                                                                      )
     radiusEnclosingMass        (i)=+massDistributionFiniteResolutionNFW_      %radiusEnclosingMass     (mass       =      mass(i)                                                                                                                                           )
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
  coordinatesScale  =[2.0d+0*radiusScale ,0.0d0,0.0d0]
  coordinatesVirial =[1.0d+2*radiusVirial,0.0d0,0.0d0]
  potential         =+massDistributionFiniteResolutionNFW_%potential(coordinates=coordinatesScale ) &
       &             -massDistributionFiniteResolutionNFW_%potential(coordinates=coordinatesVirial)
  potentialNumerical=+massDistributionFiniteResolution_   %potential(coordinates=coordinatesScale ) &
       &             -massDistributionFiniteResolution_   %potential(coordinates=coordinatesVirial)
  call Assert(                           &
       &             "Potential"       , &
       &             potential         , &
       &             potentialNumerical, &
       &      relTol=+4.0d-2             &
       &     )
  !! Total energy.
  energy         =+massDistributionFiniteResolutionNFW_%energy(radiusVirial,massDistributionFiniteResolutionNFW_)
  energyNumerical=+massDistributionFiniteResolution_   %energy(radiusVirial,massDistributionFiniteResolution_   )
  call Assert(                        &
       &             "Energy"       , &
       &             energy         , &
       &             energyNumerical, &
       &      relTol=+3.0d-3          &
       &     )
  !![
  <objectDestructor name="massDistribution_"                         />
  <objectDestructor name="massDistributionFiniteResolutionNFW_"      />
  <objectDestructor name="massDistributionFiniteResolution_"         />
  <objectDestructor name="kinematicsDistributionFiniteResolutionNFW_"/>
  <objectDestructor name="kinematicsDistributionFiniteResolution_"   />
  !!]
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
end program Test_Dark_Matter_Profiles_Finite_Resolution
