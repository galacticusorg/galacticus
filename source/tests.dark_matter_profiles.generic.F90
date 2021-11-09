!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
Contains a program to test calculations for generic dark matter profiles.
!!}

program Test_Dark_Matter_Profiles_Generic
  !!{
  Tests that numerical differentiation functions work.
  !!}
  use :: Cosmology_Functions         , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters        , only : cosmologyParametersSimple
  use :: Dark_Matter_Halo_Scales     , only : darkMatterHaloScaleVirialDensityContrastDefinition
  use :: Virial_Density_Contrast     , only : virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt
  use :: Dark_Matter_Profiles        , only : darkMatterProfile                                             , darkMatterProfileDarkMatterOnly
  use :: Dark_Matter_Profiles_DMO    , only : darkMatterProfileDMOBurkert                                   , darkMatterProfileDMOEinasto             , darkMatterProfileDMOIsothermal     , darkMatterProfileDMONFW       , &
          &                                   darkMatterProfileDMOTruncated                                 , darkMatterProfileDMOTruncatedExponential, darkMatterProfileDMOZhao1996
  use :: Dark_Matter_Profiles_Generic, only : nonAnalyticSolversNumerical
  use :: Display                     , only : displayMessage                                                , displayVerbositySet                     , verbosityLevelStandard
  use :: Events_Hooks                , only : eventsHooksInitialize
  use :: Functions_Global_Utilities  , only : Functions_Global_Set
  use :: Galacticus_Nodes            , only : nodeClassHierarchyFinalize                                    , nodeClassHierarchyInitialize            , nodeComponentBasic                 , nodeComponentDarkMatterProfile, &
          &                                   treeNode
  use :: Input_Parameters            , only : inputParameters
  use :: Node_Components             , only : Node_Components_Initialize                                    , Node_Components_Thread_Initialize       , Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Unit_Tests                  , only : Assert                                                        , Skip                                    , Unit_Tests_Begin_Group             , Unit_Tests_End_Group          , &
          &                                   Unit_Tests_Finish
  implicit none
  type            (darkMatterHaloScaleVirialDensityContrastDefinition            )               :: darkMatterHaloScale_
  type            (cosmologyParametersSimple                                     )               :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                                )               :: cosmologyFunctions_
  type            (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt)               :: virialDensityContrast_
  type            (darkMatterProfileDMOIsothermal                                ), pointer      :: darkMatterProfileIsothermal_
  type            (darkMatterProfileDMONFW                                       ), pointer      :: darkMatterProfileNFW_
  type            (darkMatterProfileDMOEinasto                                   ), pointer      :: darkMatterProfileEinasto_
  type            (darkMatterProfileDMOBurkert                                   ), pointer      :: darkMatterProfileBurkert_
  type            (darkMatterProfileDMOTruncated                                 ), pointer      :: darkMatterProfileTruncated_
  type            (darkMatterProfileDMOTruncatedExponential                      ), pointer      :: darkMatterProfileTruncatedExponential_
  type            (darkMatterProfileDMOZhao1996                                  ), pointer      :: darkMatterProfileZhao1996_
  type            (darkMatterProfileDarkMatterOnly                               ), pointer      :: darkMatterProfileIsothermal__                   , darkMatterProfileNFW__                        , &
       &                                                                                            darkMatterProfileEinasto__                      , darkMatterProfileBurkert__                    , &
       &                                                                                            darkMatterProfileTruncated__                    , darkMatterProfileTruncatedExponential__       , &
       &                                                                                            darkMatterProfileZhao1996__
  type            (treeNode                                                      ), pointer      :: node_
  class           (nodeComponentBasic                                            ), pointer      :: basic_
  class           (nodeComponentDarkMatterProfile                                ), pointer      :: darkMatterProfile_
  double precision                                                                , parameter    :: concentration                             =8.0d+0, massVirial                            =1.0d12, &
       &                                                                                            shapeProfile                              =1.8d-1
  type            (inputParameters                                               )               :: parameters
  integer                                                                                        :: i
  double precision                                                                               :: radiusScale                                      , radiusVirial                                 , &
       &                                                                                            timeDynamical
  double precision                                                                , dimension(9) :: enclosedMass                                     , enclosedMassNumerical                        , &
       &                                                                                            potentialNumerical                               , potential                                    , &
       &                                                                                            velocityCircularNumerical                        , velocityCircular                             , &
       &                                                                                            radialVelocityDispersionNumerical                , radialVelocityDispersion                     , &
       &                                                                                            kSpaceNumerical                                  , kSpace                                       , &
       &                                                                                            freefallRadiusNumerical                          , freefallRadius                               , &
       &                                                                                            freefallRadiusIncreaseRateNumerical              , freefallRadiusIncreaseRate                   , &
       &                                                                                            radiusEnclosingDensityNumerical                  , radiusEnclosingDensity                       , &
       &                                                                                            radiusEnclosingMassNumerical                     , radiusEnclosingMass                          , &
       &                                                                                            radiusFromSpecificAngularMomentumNumerical       , radiusFromSpecificAngularMomentum            , &
       &                                                                                            densityLogSlopeNumerical                         , densityLogSlope                              , &
       &                                                                                            density                                          , densityReference
  double precision                                                                , dimension(9) :: scaleFractional=[0.01d0,0.03d0,0.10d0,0.30d0,1.00d0,3.00d0,10.0d0,30.0d0,100.0d0]

  call displayVerbositySet(verbosityLevelStandard)
  call Unit_Tests_Begin_Group("Generic dark matter profiles")
  parameters=inputParameters('testSuite/parameters/darkMatterProfilesGeneric.xml')
  call eventsHooksInitialize()
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  allocate(darkMatterProfileNFW_                  )
  allocate(darkMatterProfileIsothermal_           )
  allocate(darkMatterProfileEinasto_              )
  allocate(darkMatterProfileBurkert_              )
  allocate(darkMatterProfileZhao1996_             )
  allocate(darkMatterProfileTruncated_            )
  allocate(darkMatterProfileTruncatedExponential_ )
  allocate(darkMatterProfileNFW__                 )
  allocate(darkMatterProfileIsothermal__          )
  allocate(darkMatterProfileEinasto__             )
  allocate(darkMatterProfileBurkert__             )
  allocate(darkMatterProfileZhao1996__            )
  allocate(darkMatterProfileTruncated__           )
  allocate(darkMatterProfileTruncatedExponential__)
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
  darkMatterProfileIsothermal_            =   darkMatterProfileDMOIsothermal          (                                                                  &
       &                                                                               darkMatterHaloScale_                =darkMatterHaloScale_         &
       &                                                                              )
  darkMatterProfileNFW_                   =   darkMatterProfileDMONFW                 (                                                                  &
       &                                                                               velocityDispersionUseSeriesExpansion=.false.                    , &
       &                                                                               darkMatterHaloScale_                =darkMatterHaloScale_         &
       &                                                                              )
  darkMatterProfileEinasto_               =   darkMatterProfileDMOEinasto             (                                                                  &
       &                                                                               darkMatterHaloScale_                =darkMatterHaloScale_         &
       &                                                                              )
  darkMatterProfileBurkert_               =   darkMatterProfileDMOBurkert             (                                                                  &
       &                                                                               darkMatterHaloScale_                =darkMatterHaloScale_         &
       &                                                                              )
  darkMatterProfileZhao1996_              =   darkMatterProfileDMOZhao1996            (                                                                  &
       &                                                                               alpha                               =1.0d0                      , &
       &                                                                               beta                                =3.0d0                      , &
       &                                                                               gamma                               =1.0d0                      , &
       &                                                                               darkMatterHaloScale_                =darkMatterHaloScale_         &
       &                                                                              )
  darkMatterProfileTruncated_             =   darkMatterProfileDMOTruncated           (                                                                  &
       &                                                                               radiusFractionalTruncateMinimum     = 1.0d0                     , &
       &                                                                               radiusFractionalTruncateMaximum     =20.0d0                     , &
       &                                                                               nonAnalyticSolver                   =nonAnalyticSolversNumerical, &
       &                                                                               darkMatterProfileDMO_               =darkMatterProfileNFW_      , &
       &                                                                               darkMatterHaloScale_                =darkMatterHaloScale_         &
       &                                                                              )
  darkMatterProfileTruncatedExponential_  =   darkMatterProfileDMOTruncatedExponential(                                                                  &
       &                                                                               radiusFractionalDecay               = 3.0d0                     , &
       &                                                                               alpha                               = 1.0d0                     , &
       &                                                                               beta                                = 3.0d0                     , &
       &                                                                               gamma                               = 1.0d0                     , &
       &                                                                               nonAnalyticSolver                   =nonAnalyticSolversNumerical, &
       &                                                                               darkMatterProfileDMO_               =darkMatterProfileNFW_      , &
       &                                                                               darkMatterHaloScale_                =darkMatterHaloScale_         &
       &                                                                              )
  darkMatterProfileIsothermal__           =   darkMatterProfileDarkMatterOnly(cosmologyParameters_,darkMatterHaloScale_,darkMatterProfileIsothermal_           )
  darkMatterProfileNFW__                  =   darkMatterProfileDarkMatterOnly(cosmologyParameters_,darkMatterHaloScale_,darkMatterProfileNFW_                  )
  darkMatterProfileEinasto__              =   darkMatterProfileDarkMatterOnly(cosmologyParameters_,darkMatterHaloScale_,darkMatterProfileEinasto_              )
  darkMatterProfileBurkert__              =   darkMatterProfileDarkMatterOnly(cosmologyParameters_,darkMatterHaloScale_,darkMatterProfileBurkert_              )
  darkMatterProfileZhao1996__             =   darkMatterProfileDarkMatterOnly(cosmologyParameters_,darkMatterHaloScale_,darkMatterProfileZhao1996_             )
  darkMatterProfileTruncated__            =   darkMatterProfileDarkMatterOnly(cosmologyParameters_,darkMatterHaloScale_,darkMatterProfileTruncated_            )
  darkMatterProfileTruncatedExponential__ =   darkMatterProfileDarkMatterOnly(cosmologyParameters_,darkMatterHaloScale_,darkMatterProfileTruncatedExponential_ )
  node_                                   =>  treeNode                              (                 )
  basic_                                  =>  node_               %basic            (autoCreate=.true.)
  darkMatterProfile_                      =>  node_               %darkMatterProfile(autoCreate=.true.)
  call basic_            %timeSet            (cosmologyFunctions_%cosmicTime(1.0d0))
  call basic_            %timeLastIsolatedSet(cosmologyFunctions_%cosmicTime(1.0d0))
  call basic_            %massSet            (massVirial                           )
  radiusVirial =+darkMatterHaloScale_%virialRadius      (node_)
  timeDynamical=+darkMatterHaloScale_%dynamicalTimescale(node_)
  radiusScale  =+radiusVirial &
       &        /concentration
  call darkMatterProfile_%scaleSet(radiusScale )
  call darkMatterProfile_%shapeSet(shapeProfile)
  call Unit_Tests_Begin_Group("Isothermal profile"   )
  do i=1,size(scaleFractional)
     enclosedMass                              (i)=darkMatterProfileIsothermal__%enclosedMass                              (node_,                                                                                     scaleFractional(i)*radiusScale  )
     enclosedMassNumerical                     (i)=darkMatterProfileIsothermal__%enclosedMassNumerical                     (node_,                                                                                     scaleFractional(i)*radiusScale  )
     potential                                 (i)=darkMatterProfileIsothermal__%potential                                 (node_,                                                                                     scaleFractional(i)*radiusScale  )
     potentialNumerical                        (i)=darkMatterProfileIsothermal__%potentialNumerical                        (node_,                                                                                     scaleFractional(i)*radiusScale  )
     velocityCircular                          (i)=darkMatterProfileIsothermal__%circularVelocity                          (node_,                                                                                     scaleFractional(i)*radiusScale  )
     velocityCircularNumerical                 (i)=darkMatterProfileIsothermal__%circularVelocityNumerical                 (node_,                                                                                     scaleFractional(i)*radiusScale  )
     radialVelocityDispersion                  (i)=darkMatterProfileIsothermal__%radialVelocityDispersion                  (node_,                                                                                     scaleFractional(i)*radiusScale  )
     radialVelocityDispersionNumerical         (i)=darkMatterProfileIsothermal__%radialVelocityDispersionNumerical         (node_,                                                                                     scaleFractional(i)*radiusScale  )
     kSpace                                    (i)=darkMatterProfileIsothermal__%kSpace                                    (node_,                                                                                     scaleFractional(i)/radiusScale  )
     kSpaceNumerical                           (i)=darkMatterProfileIsothermal__%kSpaceNumerical                           (node_,                                                                                     scaleFractional(i)/radiusScale  )
     freefallRadius                            (i)=darkMatterProfileIsothermal__%freefallRadius                            (node_,                                                                                     scaleFractional(i)*timeDynamical)
     freefallRadiusNumerical                   (i)=darkMatterProfileIsothermal__%freefallRadiusNumerical                   (node_,                                                                                     scaleFractional(i)*timeDynamical)
     freefallRadiusIncreaseRate                (i)=darkMatterProfileIsothermal__%freefallRadiusIncreaseRate                (node_,                                                                                     scaleFractional(i)*timeDynamical)
     freefallRadiusIncreaseRateNumerical       (i)=darkMatterProfileIsothermal__%freefallRadiusIncreaseRateNumerical       (node_,                                                                                     scaleFractional(i)*timeDynamical)
     radiusEnclosingDensity                    (i)=darkMatterProfileIsothermal__%radiusEnclosingDensity                    (node_,darkMatterProfileIsothermal__%density         (node_,scaleFractional(i)*radiusScale)                                 )
     radiusEnclosingDensityNumerical           (i)=darkMatterProfileIsothermal__%radiusEnclosingDensityNumerical           (node_,darkMatterProfileIsothermal__%density         (node_,scaleFractional(i)*radiusScale)                                 )
     radiusEnclosingMass                       (i)=darkMatterProfileIsothermal__%radiusEnclosingMass                       (node_,darkMatterProfileIsothermal__%enclosedMass    (node_,scaleFractional(i)*radiusScale)                                 )
     radiusEnclosingMassNumerical              (i)=darkMatterProfileIsothermal__%radiusEnclosingMassNumerical              (node_,darkMatterProfileIsothermal__%enclosedMass    (node_,scaleFractional(i)*radiusScale)                                 )
     radiusFromSpecificAngularMomentum         (i)=darkMatterProfileIsothermal__%radiusFromSpecificAngularMomentum         (node_,darkMatterProfileIsothermal__%circularVelocity(node_,scaleFractional(i)*radiusScale)*scaleFractional(i)*radiusScale  )
     radiusFromSpecificAngularMomentumNumerical(i)=darkMatterProfileIsothermal__%radiusFromSpecificAngularMomentumNumerical(node_,darkMatterProfileIsothermal__%circularVelocity(node_,scaleFractional(i)*radiusScale)*scaleFractional(i)*radiusScale  )
     densityLogSlope                           (i)=darkMatterProfileIsothermal__%densityLogSlope                           (node_,                                                                                     scaleFractional(i)*radiusScale  )
     densityLogSlopeNumerical                  (i)=darkMatterProfileIsothermal__%densityLogSlopeNumerical                  (node_,                                                                                     scaleFractional(i)*radiusScale  )
  end do
  potential         =potential         -potential         (1)
  potentialNumerical=potentialNumerical-potentialNumerical(1)
  call Skip  ("Energy                          , E   ","isothermal profile assumes virial equilibrium")
  call Skip  ("Energy growth rate              , ̇E   ","isothermal profile assumes virial equilibrium")
  call Skip  ("Radial moment                   , ℛ₁ ","1ˢᵗ moment diverges for isothermal profile"    )
  call Assert("Radial moment                   , ℛ₂ ",darkMatterProfileIsothermal__%radialMomentNumerical           (node_,2.0d0,0.0d0,radiusVirial),darkMatterProfileIsothermal__%radialMoment           (node_,2.0d0,0.0d0,radiusVirial),relTol=1.0d-6              )
  call Assert("Radial moment                   , ℛ₃ ",darkMatterProfileIsothermal__%radialMomentNumerical           (node_,3.0d0,0.0d0,radiusVirial),darkMatterProfileIsothermal__%radialMoment           (node_,3.0d0,0.0d0,radiusVirial),relTol=1.0d-6              )
  call Assert("Rotation normalization          , A   ",darkMatterProfileIsothermal__%rotationNormalizationNumerical  (node_                         ),darkMatterProfileIsothermal__%rotationNormalization  (node_                         ),relTol=1.0d-3              )
  call Assert("Peak circular velocity          , Vmax",darkMatterProfileIsothermal__%circularVelocityMaximumNumerical(node_                         ),darkMatterProfileIsothermal__%circularVelocityMaximum(node_                         ),relTol=1.0d-3              )
  call Assert("Enclosed mass                   , M(r)",                              enclosedMassNumerical                                           ,                              enclosedMass                                           ,relTol=2.0d-3              )
  call Assert("Potential                       , Φ(r)",                              potentialNumerical                                              ,                              potential                                              ,relTol=1.0d-3              )
  call Assert("Circular velocity               , V(r)",                              velocityCircularNumerical                                       ,                              velocityCircular                                       ,relTol=1.0d-3              )
  call Assert("Radial velocity dispersion      , σ(r)",                              radialVelocityDispersionNumerical                               ,                              radialVelocityDispersion                               ,relTol=3.0d-3              )
  call Assert("Fourier transform               , u(k)",                              kSpaceNumerical                                                 ,                              kSpace                                                 ,relTol=1.0d-3,absTol=1.0d-4)
  call Assert("Freefall radius                 , r(t)",                              freefallRadiusNumerical                                         ,                              freefallRadius                                         ,relTol=1.0d-3              )
  call Assert("Freefall radius increase rate   , ̇r(t)",                             freefallRadiusIncreaseRateNumerical                             ,                              freefallRadiusIncreaseRate                             ,relTol=1.0d-3              )
  call Assert("Radius enclosing density        , r(ρ)",                              radiusEnclosingDensityNumerical                                 ,                              radiusEnclosingDensity                                 ,relTol=1.0d-3              )
  call Assert("Radius enclosing mass           , r(M)",                              radiusEnclosingMassNumerical                                    ,                              radiusEnclosingMass                                    ,relTol=1.0d-3              )
  call Assert("Radius-specific angular momentum, r(j)",                              radiusFromSpecificAngularMomentumNumerical                      ,                              radiusFromSpecificAngularMomentum                      ,relTol=1.0d-3              )
  call Assert("Density log gradient            , α(r)",                              densityLogSlopeNumerical                                        ,                              densityLogSlope                                        ,relTol=1.0d-3              )
  call Unit_Tests_End_Group  (                       )
  call Unit_Tests_Begin_Group("NFW profile"          )
  do i=1,size(scaleFractional)
     enclosedMass                              (i)=darkMatterProfileNFW__       %enclosedMass                              (node_,                                                                                     scaleFractional(i)*radiusScale  )
     enclosedMassNumerical                     (i)=darkMatterProfileNFW__       %enclosedMassNumerical                     (node_,                                                                                     scaleFractional(i)*radiusScale  )
     potential                                 (i)=darkMatterProfileNFW__       %potential                                 (node_,                                                                                     scaleFractional(i)*radiusScale  )
     potentialNumerical                        (i)=darkMatterProfileNFW__       %potentialNumerical                        (node_,                                                                                     scaleFractional(i)*radiusScale  )
     velocityCircular                          (i)=darkMatterProfileNFW__       %circularVelocity                          (node_,                                                                                     scaleFractional(i)*radiusScale  )
     velocityCircularNumerical                 (i)=darkMatterProfileNFW__       %circularVelocityNumerical                 (node_,                                                                                     scaleFractional(i)*radiusScale  )
     radialVelocityDispersion                  (i)=darkMatterProfileNFW__       %radialVelocityDispersion                  (node_,                                                                                     scaleFractional(i)*radiusScale  )
     radialVelocityDispersionNumerical         (i)=darkMatterProfileNFW__       %radialVelocityDispersionNumerical         (node_,                                                                                     scaleFractional(i)*radiusScale  )
     kSpace                                    (i)=darkMatterProfileNFW__       %kSpace                                    (node_,                                                                                     scaleFractional(i)/radiusScale  )
     kSpaceNumerical                           (i)=darkMatterProfileNFW__       %kSpaceNumerical                           (node_,                                                                                     scaleFractional(i)/radiusScale  )
     freefallRadius                            (i)=darkMatterProfileNFW__       %freefallRadius                            (node_,                                                                                     scaleFractional(i)*timeDynamical)
     freefallRadiusNumerical                   (i)=darkMatterProfileNFW__       %freefallRadiusNumerical                   (node_,                                                                                     scaleFractional(i)*timeDynamical)
     freefallRadiusIncreaseRate                (i)=darkMatterProfileNFW__       %freefallRadiusIncreaseRate                (node_,                                                                                     scaleFractional(i)*timeDynamical)
     freefallRadiusIncreaseRateNumerical       (i)=darkMatterProfileNFW__       %freefallRadiusIncreaseRateNumerical       (node_,                                                                                     scaleFractional(i)*timeDynamical)
     radiusEnclosingDensity                    (i)=darkMatterProfileNFW__       %radiusEnclosingDensity                    (node_,darkMatterProfileNFW__       %density         (node_,scaleFractional(i)*radiusScale)                                 )
     radiusEnclosingDensityNumerical           (i)=darkMatterProfileNFW__       %radiusEnclosingDensityNumerical           (node_,darkMatterProfileNFW__       %density         (node_,scaleFractional(i)*radiusScale)                                 )
     radiusEnclosingMass                       (i)=darkMatterProfileNFW__       %radiusEnclosingMass                       (node_,darkMatterProfileNFW__       %enclosedMass    (node_,scaleFractional(i)*radiusScale)                                 )
     radiusEnclosingMassNumerical              (i)=darkMatterProfileNFW__       %radiusEnclosingMassNumerical              (node_,darkMatterProfileNFW__       %enclosedMass    (node_,scaleFractional(i)*radiusScale)                                 )
     radiusFromSpecificAngularMomentum         (i)=darkMatterProfileNFW__       %radiusFromSpecificAngularMomentum         (node_,darkMatterProfileNFW__       %circularVelocity(node_,scaleFractional(i)*radiusScale)*scaleFractional(i)*radiusScale  )
     radiusFromSpecificAngularMomentumNumerical(i)=darkMatterProfileNFW__       %radiusFromSpecificAngularMomentumNumerical(node_,darkMatterProfileNFW__       %circularVelocity(node_,scaleFractional(i)*radiusScale)*scaleFractional(i)*radiusScale  )
     densityLogSlope                           (i)=darkMatterProfileNFW__       %densityLogSlope                           (node_,                                                                                     scaleFractional(i)*radiusScale  )
     densityLogSlopeNumerical                  (i)=darkMatterProfileNFW__       %densityLogSlopeNumerical                  (node_,                                                                                     scaleFractional(i)*radiusScale  )
  end do
  potential         =potential         -potential         (1)
  potentialNumerical=potentialNumerical-potentialNumerical(1)
  call Assert("Energy                          , E   ",darkMatterProfileNFW__       %energyNumerical                 (node_                         ),darkMatterProfileNFW__       %energy                 (node_                         ),relTol=1.0d-3              )
  call Assert("Radial moment                   , ℛ₁ ",darkMatterProfileNFW__       %radialMomentNumerical           (node_,1.0d0,0.0d0,radiusVirial),darkMatterProfileNFW__       %radialMoment           (node_,1.0d0,0.0d0,radiusVirial),relTol=1.0d-6              )
  call Assert("Radial moment                   , ℛ₂ ",darkMatterProfileNFW__       %radialMomentNumerical           (node_,2.0d0,0.0d0,radiusVirial),darkMatterProfileNFW__       %radialMoment           (node_,2.0d0,0.0d0,radiusVirial),relTol=1.0d-6              )
  call Assert("Radial moment                   , ℛ₃ ",darkMatterProfileNFW__       %radialMomentNumerical           (node_,3.0d0,0.0d0,radiusVirial),darkMatterProfileNFW__       %radialMoment           (node_,3.0d0,0.0d0,radiusVirial),relTol=1.0d-6              )
  call Assert("Rotation normalization          , A   ",darkMatterProfileNFW__       %rotationNormalizationNumerical  (node_                         ),darkMatterProfileNFW__       %rotationNormalization  (node_                         ),relTol=1.0d-3              )
  call Assert("Peak circular velocity          , Vmax",darkMatterProfileNFW__       %circularVelocityMaximumNumerical(node_                         ),darkMatterProfileNFW__       %circularVelocityMaximum(node_                         ),relTol=1.0d-3              )
  call Assert("Enclosed mass                   , M(r)",                              enclosedMassNumerical                                           ,                              enclosedMass                                           ,relTol=1.0d-2              )
  call Assert("Potential                       , Φ(r)",                              potentialNumerical                                              ,                              potential                                              ,relTol=1.0d-3              )
  call Assert("Circular velocity               , V(r)",                              velocityCircularNumerical                                       ,                              velocityCircular                                       ,relTol=1.0d-3              )
  call Assert("Radial velocity dispersion      , σ(r)",                              radialVelocityDispersionNumerical                               ,                              radialVelocityDispersion                               ,relTol=3.0d-3              )
  call Assert("Fourier transform               , u(k)",                              kSpaceNumerical                                                 ,                              kSpace                                                 ,relTol=1.0d-3,absTol=1.0d-4)
  call Assert("Freefall radius                 , r(t)",                              freefallRadiusNumerical                                         ,                              freefallRadius                                         ,relTol=1.0d-3              )
  call Assert("Freefall radius increase rate   , ̇r(t)",                             freefallRadiusIncreaseRateNumerical                             ,                              freefallRadiusIncreaseRate                             ,relTol=1.0d-2              )
  call Assert("Radius enclosing density        , r(ρ)",                              radiusEnclosingDensityNumerical                                 ,                              radiusEnclosingDensity                                 ,relTol=1.0d-3              )
  call Assert("Radius enclosing mass           , r(M)",                              radiusEnclosingMassNumerical                                    ,                              radiusEnclosingMass                                    ,relTol=1.0d-3              )
  call Assert("Radius-specific angular momentum, r(j)",                              radiusFromSpecificAngularMomentumNumerical                      ,                              radiusFromSpecificAngularMomentum                      ,relTol=1.0d-3              )
  call Assert("Density log gradient            , α(r)",                              densityLogSlopeNumerical                                        ,                              densityLogSlope                                        ,relTol=1.0d-3              )
  call Unit_Tests_End_Group  (                       )
  call Unit_Tests_Begin_Group("Einasto profile"      )
  do i=1,size(scaleFractional)
     enclosedMass                              (i)=darkMatterProfileEinasto__   %enclosedMass                              (node_,                                                                                     scaleFractional(i)*radiusScale  )
     enclosedMassNumerical                     (i)=darkMatterProfileEinasto__   %enclosedMassNumerical                     (node_,                                                                                     scaleFractional(i)*radiusScale  )
     potential                                 (i)=darkMatterProfileEinasto__   %potential                                 (node_,                                                                                     scaleFractional(i)*radiusScale  )
     potentialNumerical                        (i)=darkMatterProfileEinasto__   %potentialNumerical                        (node_,                                                                                     scaleFractional(i)*radiusScale  )
     velocityCircular                          (i)=darkMatterProfileEinasto__   %circularVelocity                          (node_,                                                                                     scaleFractional(i)*radiusScale  )
     velocityCircularNumerical                 (i)=darkMatterProfileEinasto__   %circularVelocityNumerical                 (node_,                                                                                     scaleFractional(i)*radiusScale  )
     radialVelocityDispersion                  (i)=darkMatterProfileEinasto__   %radialVelocityDispersion                  (node_,                                                                                     scaleFractional(i)*radiusScale  )
     radialVelocityDispersionNumerical         (i)=darkMatterProfileEinasto__   %radialVelocityDispersionNumerical         (node_,                                                                                     scaleFractional(i)*radiusScale  )
     kSpace                                    (i)=darkMatterProfileEinasto__   %kSpace                                    (node_,                                                                                     scaleFractional(i)/radiusScale  )
     kSpaceNumerical                           (i)=darkMatterProfileEinasto__   %kSpaceNumerical                           (node_,                                                                                     scaleFractional(i)/radiusScale  )
     freefallRadius                            (i)=darkMatterProfileEinasto__   %freefallRadius                            (node_,                                                                                     scaleFractional(i)*timeDynamical)
     freefallRadiusNumerical                   (i)=darkMatterProfileEinasto__   %freefallRadiusNumerical                   (node_,                                                                                     scaleFractional(i)*timeDynamical)
     radiusEnclosingDensity                    (i)=darkMatterProfileEinasto__   %radiusEnclosingDensity                    (node_,darkMatterProfileEinasto__   %density         (node_,scaleFractional(i)*radiusScale)                                 )
     radiusEnclosingDensityNumerical           (i)=darkMatterProfileEinasto__   %radiusEnclosingDensityNumerical           (node_,darkMatterProfileEinasto__   %density         (node_,scaleFractional(i)*radiusScale)                                 )
     radiusEnclosingMass                       (i)=darkMatterProfileEinasto__   %radiusEnclosingMass                       (node_,darkMatterProfileEinasto__   %enclosedMass    (node_,scaleFractional(i)*radiusScale)                                 )
     radiusEnclosingMassNumerical              (i)=darkMatterProfileEinasto__   %radiusEnclosingMassNumerical              (node_,darkMatterProfileEinasto__   %enclosedMass    (node_,scaleFractional(i)*radiusScale)                                 )
     radiusFromSpecificAngularMomentum         (i)=darkMatterProfileEinasto__   %radiusFromSpecificAngularMomentum         (node_,darkMatterProfileEinasto__   %circularVelocity(node_,scaleFractional(i)*radiusScale)*scaleFractional(i)*radiusScale  )
     radiusFromSpecificAngularMomentumNumerical(i)=darkMatterProfileEinasto__   %radiusFromSpecificAngularMomentumNumerical(node_,darkMatterProfileEinasto__   %circularVelocity(node_,scaleFractional(i)*radiusScale)*scaleFractional(i)*radiusScale  )
     densityLogSlope                           (i)=darkMatterProfileEinasto__   %densityLogSlope                           (node_,                                                                                     scaleFractional(i)*radiusScale  )
     densityLogSlopeNumerical                  (i)=darkMatterProfileEinasto__   %densityLogSlopeNumerical                  (node_,                                                                                     scaleFractional(i)*radiusScale  )
  end do
  potential         =potential         -potential         (1)
  potentialNumerical=potentialNumerical-potentialNumerical(1)
  call Assert("Energy                          , E   ",darkMatterProfileEinasto__   %energyNumerical                 (node_                         ),darkMatterProfileEinasto__   %energy                 (node_                         ),relTol=1.0d-3              )
  call Assert("Radial moment                   , ℛ₁ ",darkMatterProfileEinasto__   %radialMomentNumerical           (node_,1.0d0,0.0d0,radiusVirial),darkMatterProfileEinasto__   %radialMoment           (node_,1.0d0,0.0d0,radiusVirial),relTol=1.0d-6              )
  call Assert("Radial moment                   , ℛ₂ ",darkMatterProfileEinasto__   %radialMomentNumerical           (node_,2.0d0,0.0d0,radiusVirial),darkMatterProfileEinasto__   %radialMoment           (node_,2.0d0,0.0d0,radiusVirial),relTol=1.0d-6              )
  call Assert("Radial moment                   , ℛ₃ ",darkMatterProfileEinasto__   %radialMomentNumerical           (node_,3.0d0,0.0d0,radiusVirial),darkMatterProfileEinasto__   %radialMoment           (node_,3.0d0,0.0d0,radiusVirial),relTol=1.0d-6              )
  call Assert("Rotation normalization          , A   ",darkMatterProfileEinasto__   %rotationNormalizationNumerical  (node_                         ),darkMatterProfileEinasto__   %rotationNormalization  (node_                         ),relTol=1.0d-3              )
  call Assert("Peak circular velocity          , Vmax",darkMatterProfileEinasto__   %circularVelocityMaximumNumerical(node_                         ),darkMatterProfileEinasto__   %circularVelocityMaximum(node_                         ),relTol=1.0d-3              )
  call Assert("Enclosed mass                   , M(r)",                              enclosedMassNumerical                                           ,                              enclosedMass                                           ,relTol=1.0d-2              )
  call Assert("Potential                       , Φ(r)",                              potentialNumerical                                              ,                              potential                                              ,relTol=1.0d-3              )
  call Assert("Circular velocity               , V(r)",                              velocityCircularNumerical                                       ,                              velocityCircular                                       ,relTol=1.0d-3              )
  call Assert("Radial velocity dispersion      , σ(r)",                              radialVelocityDispersionNumerical                               ,                              radialVelocityDispersion                               ,relTol=1.0d-2              )
  call Assert("Fourier transform               , u(k)",                              kSpaceNumerical                                                 ,                              kSpace                                                 ,relTol=1.0d-3,absTol=1.0d-4)
  call Assert("Freefall radius                 , r(t)",                              freefallRadiusNumerical                                         ,                              freefallRadius                                         ,relTol=1.0d-3              )
  call Assert("Radius enclosing density        , r(ρ)",                              radiusEnclosingDensityNumerical                                 ,                              radiusEnclosingDensity                                 ,relTol=1.0d-3              )
  call Assert("Radius enclosing mass           , r(M)",                              radiusEnclosingMassNumerical                                    ,                              radiusEnclosingMass                                    ,relTol=1.0d-3              )
  call Assert("Radius-specific angular momentum, r(j)",                              radiusFromSpecificAngularMomentumNumerical                      ,                              radiusFromSpecificAngularMomentum                      ,relTol=1.0d-3              )
  call Assert("Density log gradient            , α(r)",                              densityLogSlopeNumerical                                        ,                              densityLogSlope                                        ,relTol=1.0d-3              )
  call Unit_Tests_End_Group  (                       )
  call Unit_Tests_Begin_Group("Burkert profile"      )
  do i=1,size(scaleFractional)
     enclosedMass                              (i)=darkMatterProfileBurkert__   %enclosedMass                              (node_,                                                                                     scaleFractional(i)*radiusScale  )
     enclosedMassNumerical                     (i)=darkMatterProfileBurkert__   %enclosedMassNumerical                     (node_,                                                                                     scaleFractional(i)*radiusScale  )
     potential                                 (i)=darkMatterProfileBurkert__   %potential                                 (node_,                                                                                     scaleFractional(i)*radiusScale  )
     potentialNumerical                        (i)=darkMatterProfileBurkert__   %potentialNumerical                        (node_,                                                                                     scaleFractional(i)*radiusScale  )
     velocityCircular                          (i)=darkMatterProfileBurkert__   %circularVelocity                          (node_,                                                                                     scaleFractional(i)*radiusScale  )
     velocityCircularNumerical                 (i)=darkMatterProfileBurkert__   %circularVelocityNumerical                 (node_,                                                                                     scaleFractional(i)*radiusScale  )
     radialVelocityDispersion                  (i)=darkMatterProfileBurkert__   %radialVelocityDispersion                  (node_,                                                                                     scaleFractional(i)*radiusScale  )
     radialVelocityDispersionNumerical         (i)=darkMatterProfileBurkert__   %radialVelocityDispersionNumerical         (node_,                                                                                     scaleFractional(i)*radiusScale  )
     kSpace                                    (i)=darkMatterProfileBurkert__   %kSpace                                    (node_,                                                                                     scaleFractional(i)/radiusScale  )
     kSpaceNumerical                           (i)=darkMatterProfileBurkert__   %kSpaceNumerical                           (node_,                                                                                     scaleFractional(i)/radiusScale  )
     ! Freefall radius is only evaluated for sufficiently large radii, as this profile is cored, which means that the freefall
     ! radius-time curve becomes asympototic and difficult to evaluate numerically.
     if (i >= 4) then
        freefallRadius                         (i)=darkMatterProfileBurkert__   %freefallRadius                            (node_,                                                                                     scaleFractional(i)*timeDynamical)
        freefallRadiusNumerical                (i)=darkMatterProfileBurkert__   %freefallRadiusNumerical                   (node_,                                                                                     scaleFractional(i)*timeDynamical)
     else
        freefallRadius                         (i)=0.0d0
        freefallRadiusNumerical                (i)=0.0d0
     end if
     radiusEnclosingDensity                    (i)=darkMatterProfileBurkert__   %radiusEnclosingDensity                    (node_,darkMatterProfileBurkert__   %density         (node_,scaleFractional(i)*radiusScale)                                 )
     radiusEnclosingDensityNumerical           (i)=darkMatterProfileBurkert__   %radiusEnclosingDensityNumerical           (node_,darkMatterProfileBurkert__   %density         (node_,scaleFractional(i)*radiusScale)                                 )
     radiusEnclosingMass                       (i)=darkMatterProfileBurkert__   %radiusEnclosingMass                       (node_,darkMatterProfileBurkert__   %enclosedMass    (node_,scaleFractional(i)*radiusScale)                                 )
     radiusEnclosingMassNumerical              (i)=darkMatterProfileBurkert__   %radiusEnclosingMassNumerical              (node_,darkMatterProfileBurkert__   %enclosedMass    (node_,scaleFractional(i)*radiusScale)                                 )
     radiusFromSpecificAngularMomentum         (i)=darkMatterProfileBurkert__   %radiusFromSpecificAngularMomentum         (node_,darkMatterProfileBurkert__   %circularVelocity(node_,scaleFractional(i)*radiusScale)*scaleFractional(i)*radiusScale  )
     radiusFromSpecificAngularMomentumNumerical(i)=darkMatterProfileBurkert__   %radiusFromSpecificAngularMomentumNumerical(node_,darkMatterProfileBurkert__   %circularVelocity(node_,scaleFractional(i)*radiusScale)*scaleFractional(i)*radiusScale  )
     densityLogSlope                           (i)=darkMatterProfileBurkert__   %densityLogSlope                           (node_,                                                                                     scaleFractional(i)*radiusScale  )
     densityLogSlopeNumerical                  (i)=darkMatterProfileBurkert__   %densityLogSlopeNumerical                  (node_,                                                                                     scaleFractional(i)*radiusScale  )
  end do
  potential         =potential         -potential         (1)
  potentialNumerical=potentialNumerical-potentialNumerical(1)
  call Assert("Energy                          , E   ",darkMatterProfileBurkert__   %energyNumerical                 (node_                         ),darkMatterProfileBurkert__   %energy                 (node_                         ),relTol=1.0d-3              )
  call Assert("Radial moment                   , ℛ₁ ",darkMatterProfileBurkert__   %radialMomentNumerical           (node_,1.0d0,0.0d0,radiusVirial),darkMatterProfileBurkert__   %radialMoment           (node_,1.0d0,0.0d0,radiusVirial),relTol=1.0d-6              )
  call Assert("Radial moment                   , ℛ₂ ",darkMatterProfileBurkert__   %radialMomentNumerical           (node_,2.0d0,0.0d0,radiusVirial),darkMatterProfileBurkert__   %radialMoment           (node_,2.0d0,0.0d0,radiusVirial),relTol=1.0d-6              )
  call Assert("Radial moment                   , ℛ₃ ",darkMatterProfileBurkert__   %radialMomentNumerical           (node_,3.0d0,0.0d0,radiusVirial),darkMatterProfileBurkert__   %radialMoment           (node_,3.0d0,0.0d0,radiusVirial),relTol=1.0d-6              )
  call Assert("Rotation normalization          , A   ",darkMatterProfileBurkert__   %rotationNormalizationNumerical  (node_                         ),darkMatterProfileBurkert__   %rotationNormalization  (node_                         ),relTol=1.0d-3              )
  call Assert("Peak circular velocity          , Vmax",darkMatterProfileBurkert__   %circularVelocityMaximumNumerical(node_                         ),darkMatterProfileBurkert__   %circularVelocityMaximum(node_                         ),relTol=1.0d-3              )
  call Assert("Enclosed mass                   , M(r)",                              enclosedMassNumerical                                           ,                              enclosedMass                                           ,relTol=1.0d-2              )
  call Assert("Potential                       , Φ(r)",                              potentialNumerical                                              ,                              potential                                              ,relTol=1.0d-3              )
  call Assert("Circular velocity               , V(r)",                              velocityCircularNumerical                                       ,                              velocityCircular                                       ,relTol=1.0d-3              )
  call Assert("Radial velocity dispersion      , σ(r)",                              radialVelocityDispersionNumerical                               ,                              radialVelocityDispersion                               ,relTol=3.0d-3              )
  call Assert("Fourier transform               , u(k)",                              kSpaceNumerical                                                 ,                              kSpace                                                 ,relTol=1.0d-3,absTol=1.0d-4)
  call Assert("Freefall radius                 , r(t)",                              freefallRadiusNumerical                                         ,                              freefallRadius                                         ,relTol=1.0d-3              )
  call Assert("Radius enclosing density        , r(ρ)",                              radiusEnclosingDensityNumerical                                 ,                              radiusEnclosingDensity                                 ,relTol=1.0d-3              )
  call Assert("Radius enclosing mass           , r(M)",                              radiusEnclosingMassNumerical                                    ,                              radiusEnclosingMass                                    ,relTol=1.0d-3              )
  call Assert("Radius-specific angular momentum, r(j)",                              radiusFromSpecificAngularMomentumNumerical                      ,                              radiusFromSpecificAngularMomentum                      ,relTol=1.0d-3              )
  call Assert("Density log gradient            , α(r)",                              densityLogSlopeNumerical                                        ,                              densityLogSlope                                        ,relTol=1.0d-3              )
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Begin_Group("Truncated NFW profile")
  do i=1,size(scaleFractional)
     densityLogSlope                           (i)=darkMatterProfileTruncated__     %densityLogSlope                           (node_,                                                                                         scaleFractional(i)*radiusScale  )
     densityLogSlopeNumerical                  (i)=darkMatterProfileTruncated__     %densityLogSlopeNumerical                  (node_,                                                                                         scaleFractional(i)*radiusScale  )
     radialVelocityDispersion                  (i)=darkMatterProfileTruncated__     %radialVelocityDispersion                  (node_,                                                                                         scaleFractional(i)*radiusScale  )
     radialVelocityDispersionNumerical         (i)=darkMatterProfileTruncated__     %radialVelocityDispersionNumerical         (node_,                                                                                         scaleFractional(i)*radiusScale  )
  end do
  call Skip  ("Energy                          , E   ","implemented numerically"                                                                                                                                                                                       )
  call Skip  ("Energy growth rate              , ̇E   ","implemented numerically"                                                                                                                                                                                       )
  call Skip  ("Radial moment                   , ℛ₁ ","implemented numerically"                                                                                                                                                                                       )
  call Skip  ("Radial moment                   , ℛ₂ ","implemented numerically"                                                                                                                                                                                       )
  call Skip  ("Radial moment                   , ℛ₃ ","implemented numerically"                                                                                                                                                                                       )
  call Skip  ("Rotation normalization          , A   ","implemented numerically"                                                                                                                                                                                       )
  call Skip  ("Peak circular velocity          , Vmax","implemented numerically"                                                                                                                                                                                       )
  call Skip  ("Potential                       , Φ(r)","implemented numerically"                                                                                                                                                                                       )
  call Skip  ("Circular velocity               , V(r)","implemented numerically"                                                                                                                                                                                       )
  call Assert("Radial velocity dispersion      , σ(r)",                              radialVelocityDispersionNumerical                               ,                              radialVelocityDispersion                               ,relTol=1.0d-2              )
  call Skip  ("Fourier transform               , u(k)","implemented numerically"                                                                                                                                                                                       )
  call Skip  ("Freefall radius                 , r(t)","implemented numerically"                                                                                                                                                                                       )
  call Skip  ("Freefall radius increase rate   , ̇r(t)","implemented numerically"                                                                                                                                                                                       )
  call Skip  ("Radius enclosing density        , r(ρ)","implemented numerically"                                                                                                                                                                                       )
  call Skip  ("Radius enclosing mass           , r(M)","implemented numerically"                                                                                                                                                                                       )
  call Skip  ("Radius-specific angular momentum, r(j)","implemented numerically"                                                                                                                                                                                       )
  call Assert("Density log gradient            , α(r)",                              densityLogSlopeNumerical                                        ,                              densityLogSlope                                        ,relTol=1.0d-3              )
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Begin_Group("Exponentially-truncated NFW profile")
  do i=1,size(scaleFractional)
     densityLogSlope                           (i)=darkMatterProfileTruncatedExponential__     %densityLogSlope                           (node_,                                                                                         scaleFractional(i)*radiusScale  )
     densityLogSlopeNumerical                  (i)=darkMatterProfileTruncatedExponential__     %densityLogSlopeNumerical                  (node_,                                                                                         scaleFractional(i)*radiusScale  )
     radialVelocityDispersion                  (i)=darkMatterProfileTruncatedExponential__     %radialVelocityDispersion                  (node_,                                                                                         scaleFractional(i)*radiusScale  )
     radialVelocityDispersionNumerical         (i)=darkMatterProfileTruncatedExponential__     %radialVelocityDispersionNumerical         (node_,                                                                                         scaleFractional(i)*radiusScale  )
  end do
  call Skip  ("Energy                          , E   ","implemented numerically"                                                                                                                                                                                       )
  call Skip  ("Energy growth rate              , ̇E   ","implemented numerically"                                                                                                                                                                                       )
  call Skip  ("Radial moment                   , ℛ₁ ","implemented numerically"                                                                                                                                                                                       )
  call Skip  ("Radial moment                   , ℛ₂ ","implemented numerically"                                                                                                                                                                                       )
  call Skip  ("Radial moment                   , ℛ₃ ","implemented numerically"                                                                                                                                                                                       )
  call Skip  ("Rotation normalization          , A   ","implemented numerically"                                                                                                                                                                                       )
  call Skip  ("Peak circular velocity          , Vmax","implemented numerically"                                                                                                                                                                                       )
  call Skip  ("Potential                       , Φ(r)","implemented numerically"                                                                                                                                                                                       )
  call Skip  ("Circular velocity               , V(r)","implemented numerically"                                                                                                                                                                                       )
  call Assert("Radial velocity dispersion      , σ(r)",                              radialVelocityDispersionNumerical                               ,                              radialVelocityDispersion                               ,relTol=3.0d-3              )
  call Skip  ("Fourier transform               , u(k)","implemented numerically"                                                                                                                                                                                       )
  call Skip  ("Freefall radius                 , r(t)","implemented numerically"                                                                                                                                                                                       )
  call Skip  ("Freefall radius increase rate   , ̇r(t)","implemented numerically"                                                                                                                                                                                       )
  call Skip  ("Radius enclosing density        , r(ρ)","implemented numerically"                                                                                                                                                                                       )
  call Skip  ("Radius enclosing mass           , r(M)","implemented numerically"                                                                                                                                                                                       )
  call Skip  ("Radius-specific angular momentum, r(j)","implemented numerically"                                                                                                                                                                                       )
  call Assert("Density log gradient            , α(r)",                              densityLogSlopeNumerical                                        ,                              densityLogSlope                                        ,relTol=1.0d-3              )
  call Unit_Tests_End_Group               ()
  ! The Zhao1996 profile is a generalized NFW-type profile. So, compare to the NFW profile.
  call Unit_Tests_Begin_Group("Zhao1996 profile"          )
  do i=1,size(scaleFractional)
     density                                   (i)=darkMatterProfileZhao1996__  %density                          (node_,                                                                                     scaleFractional(i)*radiusScale  )
     densityReference                          (i)=darkMatterProfileNFW__       %density                          (node_,                                                                                     scaleFractional(i)*radiusScale  )
     enclosedMass                              (i)=darkMatterProfileZhao1996__  %enclosedMass                     (node_,                                                                                     scaleFractional(i)*radiusScale  )
     enclosedMassNumerical                     (i)=darkMatterProfileNFW__       %enclosedMass                     (node_,                                                                                     scaleFractional(i)*radiusScale  )
     potential                                 (i)=darkMatterProfileZhao1996__  %potential                        (node_,                                                                                     scaleFractional(i)*radiusScale  )
     potentialNumerical                        (i)=darkMatterProfileNFW__       %potential                        (node_,                                                                                     scaleFractional(i)*radiusScale  )
     velocityCircular                          (i)=darkMatterProfileZhao1996__  %circularVelocity                 (node_,                                                                                     scaleFractional(i)*radiusScale  )
     velocityCircularNumerical                 (i)=darkMatterProfileNFW__       %circularVelocity                 (node_,                                                                                     scaleFractional(i)*radiusScale  )
     radialVelocityDispersion                  (i)=darkMatterProfileZhao1996__  %radialVelocityDispersion         (node_,                                                                                     scaleFractional(i)*radiusScale  )
     radialVelocityDispersionNumerical         (i)=darkMatterProfileNFW__       %radialVelocityDispersion         (node_,                                                                                     scaleFractional(i)*radiusScale  )
     kSpace                                    (i)=darkMatterProfileZhao1996__  %kSpace                           (node_,                                                                                     scaleFractional(i)/radiusScale  )
     kSpaceNumerical                           (i)=darkMatterProfileNFW__       %kSpace                           (node_,                                                                                     scaleFractional(i)/radiusScale  )
     freefallRadius                            (i)=darkMatterProfileZhao1996__  %freefallRadius                   (node_,                                                                                     scaleFractional(i)*timeDynamical)
     freefallRadiusNumerical                   (i)=darkMatterProfileNFW__       %freefallRadius                   (node_,                                                                                     scaleFractional(i)*timeDynamical)
     freefallRadiusIncreaseRate                (i)=darkMatterProfileZhao1996__  %freefallRadiusIncreaseRate       (node_,                                                                                     scaleFractional(i)*timeDynamical)
     freefallRadiusIncreaseRateNumerical       (i)=darkMatterProfileNFW__       %freefallRadiusIncreaseRate       (node_,                                                                                     scaleFractional(i)*timeDynamical)
     radiusEnclosingDensity                    (i)=darkMatterProfileZhao1996__  %radiusEnclosingDensity           (node_,darkMatterProfileZhao1996__  %density         (node_,scaleFractional(i)*radiusScale)                                 )
     radiusEnclosingDensityNumerical           (i)=darkMatterProfileNFW__       %radiusEnclosingDensityNumerical  (node_,darkMatterProfileNFW__       %density         (node_,scaleFractional(i)*radiusScale)                                 )
     radiusEnclosingMass                       (i)=darkMatterProfileZhao1996__  %radiusEnclosingMass              (node_,darkMatterProfileZhao1996__  %enclosedMass    (node_,scaleFractional(i)*radiusScale)                                 )
     radiusEnclosingMassNumerical              (i)=darkMatterProfileNFW__       %radiusEnclosingMass              (node_,darkMatterProfileNFW__       %enclosedMass    (node_,scaleFractional(i)*radiusScale)                                 )
     radiusFromSpecificAngularMomentum         (i)=darkMatterProfileZhao1996__  %radiusFromSpecificAngularMomentum(node_,darkMatterProfileZhao1996__  %circularVelocity(node_,scaleFractional(i)*radiusScale)*scaleFractional(i)*radiusScale  )
     radiusFromSpecificAngularMomentumNumerical(i)=darkMatterProfileNFW__       %radiusFromSpecificAngularMomentum(node_,darkMatterProfileNFW__       %circularVelocity(node_,scaleFractional(i)*radiusScale)*scaleFractional(i)*radiusScale  )
     densityLogSlope                           (i)=darkMatterProfileZhao1996__  %densityLogSlope                  (node_,                                                                                     scaleFractional(i)*radiusScale  )
     densityLogSlopeNumerical                  (i)=darkMatterProfileNFW__       %densityLogSlope                  (node_,                                                                                     scaleFractional(i)*radiusScale  )
  end do
  potential         =potential         -potential         (1)
  potentialNumerical=potentialNumerical-potentialNumerical(1)
  call Assert("Energy                          , E   ",darkMatterProfileNFW__%energy                          (node_                         ),darkMatterProfileZhao1996__  %energy                 (node_                         ),relTol=1.0d-3              )
  call Assert("Radial moment                   , ℛ₁ ",darkMatterProfileNFW__%radialMoment                    (node_,1.0d0,0.0d0,radiusVirial),darkMatterProfileZhao1996__  %radialMoment           (node_,1.0d0,0.0d0,radiusVirial),relTol=1.0d-6              )
  call Assert("Radial moment                   , ℛ₂ ",darkMatterProfileNFW__%radialMoment                    (node_,2.0d0,0.0d0,radiusVirial),darkMatterProfileZhao1996__  %radialMoment           (node_,2.0d0,0.0d0,radiusVirial),relTol=1.0d-6              )
  call Assert("Radial moment                   , ℛ₃ ",darkMatterProfileNFW__%radialMoment                    (node_,3.0d0,0.0d0,radiusVirial),darkMatterProfileZhao1996__  %radialMoment           (node_,3.0d0,0.0d0,radiusVirial),relTol=1.0d-6              )
  call Assert("Rotation normalization          , A   ",darkMatterProfileNFW__%rotationNormalization           (node_                         ),darkMatterProfileZhao1996__  %rotationNormalization  (node_                         ),relTol=1.0d-3              )
  call Assert("Peak circular velocity          , Vmax",darkMatterProfileNFW__%circularVelocityMaximum         (node_                         ),darkMatterProfileZhao1996__  %circularVelocityMaximum(node_                         ),relTol=1.0d-3              )
  call Assert("Enclosed mass                   , M(r)",                       enclosedMassNumerical                                           ,                              enclosedMass                                           ,relTol=1.0d-2              )
  call Assert("Potential                       , Φ(r)",                       potentialNumerical                                              ,                              potential                                              ,relTol=1.0d-3              )
  call Assert("Circular velocity               , V(r)",                       velocityCircularNumerical                                       ,                              velocityCircular                                       ,relTol=1.0d-3              )
  call Assert("Radial velocity dispersion      , σ(r)",                       radialVelocityDispersionNumerical                               ,                              radialVelocityDispersion                               ,relTol=3.0d-3              )
  call Assert("Fourier transform               , u(k)",                       kSpaceNumerical                                                 ,                              kSpace                                                 ,relTol=1.0d-3,absTol=1.0d-4)
  call Assert("Freefall radius                 , r(t)",                       freefallRadiusNumerical                                         ,                              freefallRadius                                         ,relTol=1.0d-3              )
  call Assert("Freefall radius increase rate   , ̇r(t)",                      freefallRadiusIncreaseRateNumerical                             ,                              freefallRadiusIncreaseRate                             ,relTol=1.0d-2              )
  call Assert("Radius enclosing density        , r(ρ)",                       radiusEnclosingDensityNumerical                                 ,                              radiusEnclosingDensity                                 ,relTol=1.0d-3              )
  call Assert("Radius enclosing mass           , r(M)",                       radiusEnclosingMassNumerical                                    ,                              radiusEnclosingMass                                    ,relTol=1.0d-3              )
  call Assert("Radius-specific angular momentum, r(j)",                       radiusFromSpecificAngularMomentumNumerical                      ,                              radiusFromSpecificAngularMomentum                      ,relTol=1.0d-3              )
  call Assert("Density                         , ρ(r)",                       densityReference                                                ,                              density                                                ,relTol=1.0d-3              )
  call Assert("Density log gradient            , α(r)",                       densityLogSlopeNumerical                                        ,                              densityLogSlope                                        ,relTol=1.0d-3              )
  call Unit_Tests_End_Group               ()
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
end program Test_Dark_Matter_Profiles_Generic
