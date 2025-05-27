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
Contains a program to test calculations for Chowdhury2021 dark matter profiles.
!!}

program Test_Dark_Matter_Profiles_Chowdhury2021
  !!{
  Test calculations for Chowdhury dark matter profiles.
  !!}
  use :: Calculations_Resets       , only : Calculations_Reset
  use :: Cosmology_Functions       , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters      , only : cosmologyParametersSimple
  use :: Coordinates               , only : coordinateSpherical                                           , assignment(=)
  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScaleVirialDensityContrastDefinition
  use :: Dark_Matter_Particles     , only : darkMatterParticleFuzzyDarkMatter
  use :: Virial_Density_Contrast   , only : virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt
  use :: Dark_Matter_Profiles_DMO  , only : darkMatterProfileDMOClass, darkMatterProfileDMOSolitonNFW
  use :: Mass_Distributions        , only : massDistributionClass, massDistributionSolitonNFW
  use :: Display                   , only : displayMessage, displayVerbositySet, verbosityLevelStandard
  use :: Events_Hooks              , only : eventsHooksInitialize
  use :: Functions_Global_Utilities, only : Functions_Global_Set
  use :: Galacticus_Nodes          , only : nodeClassHierarchyFinalize, nodeClassHierarchyInitialize, nodeComponentBasic, nodeComponentDarkMatterProfile, treeNode, nodeComponentSatellite
  use :: Input_Parameters          , only : inputParameters
  use :: Node_Components           , only : Node_Components_Initialize, Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Numerical_Ranges          , only : Make_Range, rangeTypeLogarithmic
  use :: Unit_Tests                , only : Assert, Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  use :: Error                     , only : Error_Report

  implicit none
  type            (darkMatterHaloScaleVirialDensityContrastDefinition            )                        :: darkMatterHaloScale_
  type            (cosmologyParametersSimple                                     )                        :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                                )                        :: cosmologyFunctions_
  type            (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt)                        :: virialDensityContrast_
  type            (darkMatterProfileDMOSolitonNFW                                ), target                :: darkMatterProfileChowdhury2021_
  class           (darkMatterProfileDMOClass                                     ), pointer               :: darkMatterProfileDMO_
  class           (massDistributionClass                                         ), pointer               :: massDistribution_
  type            (darkMatterParticleFuzzyDarkMatter                             )                        :: darkMatterParticle_
  type            (treeNode                                                      ), pointer               :: node_
  class           (nodeComponentBasic                                            ), pointer               :: basic_
  class           (nodeComponentDarkMatterProfile                                ), pointer               :: darkMatterProfile_
  type            (inputParameters                                               )                        :: parameters
  double precision                                                                , parameter             :: concentration=+4.9d0, massVirial=+6.6d9
  integer                                                                         , parameter             :: countRadii   =200
  double precision                                                                , parameter             :: radiiMinimum =1.0d-6, radiiMaximum=+1.0d3
  double precision                                                                                        :: r_c=0.72d-3, r_s=1.0d-2, r_sol=1.944d-3, rho_s=5.5d14, rho_0=1.04d17
  double precision                                                                                        :: radiusVirial, radiusScale
  double precision                                                                , dimension(countRadii) :: radii, rho
  integer                                                                                                 :: i
  double precision                                                                                        :: potentialNumerical                          , potential
  type            (coordinateSpherical                                           )                        :: coordinates                                 , coordinatesReference
  
  call displayVerbositySet(verbosityLevelStandard)
  call Unit_Tests_Begin_Group("Chowdhury2021 dark matter profiles")
  parameters=inputParameters('testSuite/parameters/darkMatterProfilesChowdhury2021.xml')
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
     &amp;                                                         OmegaMatter                         = 0.31530d0                , &amp;
     &amp;                                                         OmegaBaryon                         = 0.04930d0                , &amp;
     &amp;                                                         OmegaDarkEnergy                     = 0.68470d0                , &amp;
     &amp;                                                         temperatureCMB                      = 2.72548d0                , &amp;
     &amp;                                                         HubbleConstant                      =67.36d0                     &amp;
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
  <referenceConstruct object="darkMatterParticle_"                >
   <constructor>
    darkMatterParticleFuzzyDarkMatter                             (                                                                  &amp;
     &amp;                                                        mass                               = 0.8d0                       , &amp;
     &amp;                                                        densityFraction                    = 1.0d0                         &amp;
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
  <referenceConstruct object="darkMatterProfileChowdhury2021_"    >
   <constructor>
    darkMatterProfileDMOSolitonNFW                                (                                                                  &amp;
    &amp;                                                          paramsUseInput                      =.true.                     , &amp;
    &amp;                                                          darkMatterParticle_                 =darkMatterParticle_        , &amp;
    &amp;                                                          darkMatterHaloScale_                =darkMatterHaloScale_       , &amp;
    &amp;                                                          cosmologyFunctions_                 =cosmologyFunctions_        , &amp;
    &amp;                                                          cosmologyParameters_                =cosmologyParameters_         &amp;
    &amp;                                                        )
   </constructor>
  </referenceConstruct>
  !!]
  call basic_    %timeSet            (cosmologyFunctions_%cosmicTime(1.0d0))
  call basic_    %timeLastIsolatedSet(cosmologyFunctions_%cosmicTime(1.0d0))
  call basic_    %massSet            (massVirial                           )
  radiusVirial =+darkMatterHaloScale_%radiusVirial(node_)
  radiusScale  =+radiusVirial/ concentration
  call darkMatterProfile_%scaleSet(radiusScale)
  call Calculations_Reset(node_)
  
  call Unit_Tests_Begin_Group("Chowdhury2021 profile")
  radii=Make_Range(radiiMinimum,radiiMaximum,countRadii,rangeTypeLogarithmic)*radiusScale

  darkMatterProfileDMO_      => darkMatterProfileChowdhury2021_

  select type (darkMatterProfileDMO_)
  type is (darkMatterProfileDMOSolitonNFW)
     call darkMatterProfileDMO_                                %paraminput            ( r_c, r_s, r_sol, rho_0, rho_s)
  end select
  massDistribution_          => darkMatterProfileDMO_          %get                   (node_)

  select type (massDistribution_)
     type is (massDistributionSolitonNFW)
        coordinates         =[2.0d+0*radiusScale ,0.0d0,0.0d0]
        coordinatesReference=[1.0d+2*radiusVirial,0.0d0,0.0d0]
        potential           =+massDistribution_%potential         (coordinates=coordinates         ) &
             &               -massDistribution_%potential         (coordinates=coordinatesReference)
        potentialNumerical  =+massDistribution_%potentialNumerical(coordinates=coordinates         ) &
             &               -massDistribution_%potentialNumerical(coordinates=coordinatesReference)
        call Assert(                           &
             &             "Potential"       , &
             &             potential         , &
             &             potentialNumerical, &
             &      relTol=+2.1d-2             &
             &     )
        do i=1, countRadii
           coordinates = [radii(i) ,0.0d0, 0.0d0]
           rho(i)      = massDistribution_%density(coordinates)
        end do
        open(unit=99, file="density_profile.txt", status="replace")
           do i = 1, countRadii
           write(99,"(ES15.8E2, 2X, ES25.16E2)") radii(i), rho(i)
           end do
        close(99)
  end select

  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
end program Test_Dark_Matter_Profiles_Chowdhury2021
