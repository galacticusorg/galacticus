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
Contains a program to test calculations for Zhao1996 dark matter profiles.
!!}

program Test_Dark_Matter_Profiles_Zhao1996
  !!{
  Test calculations for Zhao1996 dark matter profiles.
  !!}
  use :: Calculations_Resets       , only : Calculations_Reset
  use :: Cosmology_Functions       , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters      , only : cosmologyParametersSimple
  use :: Coordinates               , only : coordinateSpherical                                           , assignment(=)
  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScaleVirialDensityContrastDefinition
  use :: Virial_Density_Contrast   , only : virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt
  use :: Dark_Matter_Profiles_DMO  , only : darkMatterProfileDMOClass                                     , darkMatterProfileDMOZhao1996
  use :: Mass_Distributions        , only : nonAnalyticSolversNumerical                                   , massDistributionClass            , kinematicsDistributionClass        , massDistributionSpherical
  use :: Display                   , only : displayMessage                                                , displayVerbositySet              , verbosityLevelStandard
  use :: Events_Hooks              , only : eventsHooksInitialize
  use :: Functions_Global_Utilities, only : Functions_Global_Set
  use :: Galacticus_Nodes          , only : nodeClassHierarchyFinalize                                    , nodeClassHierarchyInitialize     , nodeComponentBasic                 , nodeComponentDarkMatterProfile, &
          &                                 treeNode                                                      , nodeComponentSatellite
  use :: Input_Parameters          , only : inputParameters
  use :: Node_Components           , only : Node_Components_Initialize                                    , Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Numerical_Ranges          , only : Make_Range                                                    , rangeTypeLogarithmic
  use :: Numerical_Constants_Math  , only : Pi
  use :: Unit_Tests                , only : Assert                                                        , Unit_Tests_Begin_Group           , Unit_Tests_End_Group               , Unit_Tests_Finish
  use :: Error                     , only : Error_Report
  implicit none
  type            (darkMatterHaloScaleVirialDensityContrastDefinition            ), pointer               :: darkMatterHaloScale_
  type            (cosmologyParametersSimple                                     ), pointer               :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                                ), pointer               :: cosmologyFunctions_
  type            (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt), pointer               :: virialDensityContrast_
  type            (darkMatterProfileDMOZhao1996                                  ), pointer               :: darkMatterProfileZhao1996General_           , darkMatterProfileZhao1996NFW_        , &
       &                                                                                                     darkMatterProfileZhao1996CoredNFW_          , darkMatterProfileZhao1996Gamma0_5NFW_, &
       &                                                                                                     darkMatterProfileZhao1996Gamma1_5NFW_
  class           (darkMatterProfileDMOClass                                     ), pointer               :: darkMatterProfileZhao1996_
  class           (massDistributionClass                                         ), pointer               :: massDistribution_
  class           (kinematicsDistributionClass                                   ), pointer               :: kinematicsDistribution_
  type            (treeNode                                                      ), pointer               :: node_
  class           (nodeComponentBasic                                            ), pointer               :: basic_
  class           (nodeComponentDarkMatterProfile                                ), pointer               :: darkMatterProfile_
  type            (inputParameters                                               )                        :: parameters
  double precision                                                                , parameter             :: concentration                        =+8.0d0, massVirial        =+1.0d12
  integer                                                                         , parameter             :: countRadii                           =20
  double precision                                                                , parameter             :: radiiMinimum                         =1.0d-6, radiiMaximum      =+1.0d03
  double precision                                                                , dimension(countRadii) :: massNumerical                               , mass                                 , &
       &                                                                                                     velocityDispersionNumerical                 , velocityDispersion                   , &
       &                                                                                                     radii
  integer                                                                                                 :: i                                           , j  
  double precision                                                                                        :: radiusScale                                 , radiusVirial                         , &
       &                                                                                                     potentialNumerical                          , potential
  type            (coordinateSpherical                                           )                        :: coordinates                                 , coordinatesReference
  
  call displayVerbositySet(verbosityLevelStandard)
  call Unit_Tests_Begin_Group("Zhao1996 dark matter profiles")
  parameters=inputParameters('testSuite/parameters/darkMatterProfilesZhao1996.xml')
  call eventsHooksInitialize()
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  node_              =>  treeNode                  (                 )
  basic_             =>  node_   %basic            (autoCreate=.true.)
  darkMatterProfile_ =>  node_   %darkMatterProfile(autoCreate=.true.)
  allocate(darkMatterHaloScale_                 )
  allocate(cosmologyParameters_                 )
  allocate(cosmologyFunctions_                  )
  allocate(virialDensityContrast_               )
  allocate(darkMatterProfileZhao1996General_    )
  allocate(darkMatterProfileZhao1996NFW_        )  
  allocate(darkMatterProfileZhao1996CoredNFW_   )
  allocate(darkMatterProfileZhao1996Gamma0_5NFW_)
  allocate(darkMatterProfileZhao1996Gamma1_5NFW_)
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
  <referenceConstruct object="darkMatterProfileZhao1996General_"    >
   <constructor>
    darkMatterProfileDMOZhao1996                                  (                                                                 &amp;
     &amp;                                                         alpha                               =1.0d0                     , &amp;
     &amp;                                                         beta                                =4.0d0                     , &amp;
     &amp;                                                         gamma                               =1.0d0                     , &amp;
     &amp;                                                         darkMatterHaloScale_                =darkMatterHaloScale_        &amp;
     &amp;                                                        )
     </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfileZhao1996NFW_"        >
   <constructor>
    darkMatterProfileDMOZhao1996                                  (                                                                 &amp;
     &amp;                                                         alpha                               =1.0d0                     , &amp;
     &amp;                                                         beta                                =3.0d0                     , &amp;
     &amp;                                                         gamma                               =1.0d0                     , &amp;
     &amp;                                                         darkMatterHaloScale_                =darkMatterHaloScale_        &amp;
     &amp;                                                        )
     </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfileZhao1996CoredNFW_"   >
   <constructor>
    darkMatterProfileDMOZhao1996                                  (                                                                 &amp;
     &amp;                                                         alpha                               =1.0d0                     , &amp;
     &amp;                                                         beta                                =3.0d0                     , &amp;
     &amp;                                                         gamma                               =0.0d0                     , &amp;
     &amp;                                                         darkMatterHaloScale_                =darkMatterHaloScale_        &amp;
     &amp;                                                        )
     </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfileZhao1996Gamma0_5NFW_">
   <constructor>
    darkMatterProfileDMOZhao1996                                  (                                                                 &amp;
     &amp;                                                         alpha                               =1.0d0                     , &amp;
     &amp;                                                         beta                                =3.0d0                     , &amp;
     &amp;                                                         gamma                               =0.5d0                     , &amp;
     &amp;                                                         darkMatterHaloScale_                =darkMatterHaloScale_        &amp;
     &amp;                                                        )
     </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfileZhao1996Gamma1_5NFW_">
   <constructor>
    darkMatterProfileDMOZhao1996                                  (                                                                 &amp;
     &amp;                                                         alpha                               =1.0d0                     , &amp;
     &amp;                                                         beta                                =3.0d0                     , &amp;
     &amp;                                                         gamma                               =1.5d0                     , &amp;
     &amp;                                                         darkMatterHaloScale_                =darkMatterHaloScale_        &amp;
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
  call Unit_Tests_Begin_Group("Zhao1996 profile")
  radii=Make_Range(radiiMinimum,radiiMaximum,countRadii,rangeTypeLogarithmic)*radiusScale
  do i=1,5
     select case(i)
     case (1)
        call Unit_Tests_Begin_Group("General"    )
        darkMatterProfileZhao1996_ => darkMatterProfileZhao1996General_
     case (2)
        call Unit_Tests_Begin_Group("NFW"        )
        darkMatterProfileZhao1996_ => darkMatterProfileZhao1996NFW_
     case (3)
        call Unit_Tests_Begin_Group("CoredNFW"   )
        darkMatterProfileZhao1996_ => darkMatterProfileZhao1996CoredNFW_
     case (4)
        call Unit_Tests_Begin_Group("Gamma0_5NFW")
        darkMatterProfileZhao1996_ => darkMatterProfileZhao1996Gamma0_5NFW_
     case (5)
        call Unit_Tests_Begin_Group("Gamma1_5NFW")
        darkMatterProfileZhao1996_ => darkMatterProfileZhao1996Gamma1_5NFW_
     case default
        call Error_Report('unknown profile'//{introspection:location})
     end select
     massDistribution_       => darkMatterProfileZhao1996_%get                   (node_)
     kinematicsDistribution_ => massDistribution_         %kinematicsDistribution(     )
     select type (massDistribution_)
     class is (massDistributionSpherical)
        do j=1,countRadii
           coordinates                   =[radii(j),0.0d0,0.0d0]
           mass                       (j)=+massDistribution_      %massEnclosedBySphere         (radius     =      radii(j)                                                                                )
           massNumerical              (j)=+massDistribution_      %massEnclosedBySphereNumerical(radius     =      radii(j)                                                                                )
           velocityDispersion         (j)=+kinematicsDistribution_%velocityDispersion1D         (coordinates=coordinates   ,massDistribution_=massDistribution_,massDistributionEmbedding=massDistribution_)
           velocityDispersionNumerical(j)=+kinematicsDistribution_%velocityDispersion1DNumerical(coordinates=coordinates   ,massDistribution_=massDistribution_,massDistributionEmbedding=massDistribution_)
        end do
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
     end select
     !![
     <objectDestructor name="      massDistribution_"/>
     <objectDestructor name="kinematicsDistribution_"/>
     !!]
     call Unit_Tests_End_Group            ()
  end do
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
  ! Clean up.
  call node_%destroy()
  deallocate(node_)
  !![
  <objectDestructor name="darkMatterHaloScale_"                 />
  <objectDestructor name="cosmologyParameters_"                 />
  <objectDestructor name="cosmologyFunctions_"                  />
  <objectDestructor name="virialDensityContrast_"               />
  <objectDestructor name="darkMatterProfileZhao1996General_"    />
  <objectDestructor name="darkMatterProfileZhao1996NFW_"        />  
  <objectDestructor name="darkMatterProfileZhao1996CoredNFW_"   />
  <objectDestructor name="darkMatterProfileZhao1996Gamma0_5NFW_"/>
  <objectDestructor name="darkMatterProfileZhao1996Gamma1_5NFW_"/>
  !!]
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
end program Test_Dark_Matter_Profiles_Zhao1996
