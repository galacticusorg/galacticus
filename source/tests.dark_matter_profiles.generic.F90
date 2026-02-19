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
Contains a program to test calculations for generic dark matter profiles.
!!}

program Test_Dark_Matter_Profiles_Generic
  !!{
  Tests calculations for generic dark matter profiles.
  !!}
  use :: Coordinates               , only : coordinateSpherical                                           , assignment(=)
  use :: Cosmology_Functions       , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters      , only : cosmologyParametersSimple
  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScaleVirialDensityContrastDefinition
  use :: Virial_Density_Contrast   , only : virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt
  use :: Dark_Matter_Profiles      , only : darkMatterProfile                                             , darkMatterProfileDarkMatterOnly
  use :: Dark_Matter_Profiles_DMO  , only : darkMatterProfileDMOBurkert                                   , darkMatterProfileDMOEinasto             , darkMatterProfileDMOIsothermal     , darkMatterProfileDMONFW       , &
          &                                 darkMatterProfileDMOTruncated                                 , darkMatterProfileDMOTruncatedExponential, darkMatterProfileDMOZhao1996
  use :: Mass_Distributions        , only : nonAnalyticSolversNumerical                                   , massDistributionClass                   , kinematicsDistributionClass        , massDistributionSpherical
  use :: Display                   , only : displayMessage                                                , displayVerbositySet                     , verbosityLevelStandard
  use :: Events_Hooks              , only : eventsHooksInitialize
  use :: Functions_Global_Utilities, only : Functions_Global_Set
  use :: Galacticus_Nodes          , only : nodeClassHierarchyFinalize                                    , nodeClassHierarchyInitialize            , nodeComponentBasic                 , nodeComponentDarkMatterProfile, &
          &                                 treeNode
  use :: Input_Parameters          , only : inputParameters
  use :: Node_Components           , only : Node_Components_Initialize                                    , Node_Components_Thread_Initialize       , Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Unit_Tests                , only : Assert                                                        , Skip                                    , Unit_Tests_Begin_Group             , Unit_Tests_End_Group          , &
          &                                 Unit_Tests_Finish
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
  class           (massDistributionClass                                         ), pointer      :: massDistribution_
  class           (kinematicsDistributionClass                                   ), pointer      :: kinematicsDistribution_
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
       &                                                                                            densityLogSlopeNumerical                         , densityLogSlope
  double precision                                                                , dimension(9) :: scaleFractional=[0.01d0,0.03d0,0.10d0,0.30d0,1.00d0,3.00d0,10.0d0,30.0d0,100.0d0]
  type            (coordinateSpherical                                           )               :: coordinates                                      , coordinatesOuter
  
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
       &                                                                               nonAnalyticSolver                   =nonAnalyticSolversNumerical, &
       &                                                                               darkMatterProfileDMO_               =darkMatterProfileNFW_      , &
       &                                                                               darkMatterHaloScale_                =darkMatterHaloScale_         &
       &                                                                              )
  darkMatterProfileIsothermal__           =   darkMatterProfileDarkMatterOnly(.true.,cosmologyParameters_,darkMatterProfileIsothermal_          )
  darkMatterProfileNFW__                  =   darkMatterProfileDarkMatterOnly(.true.,cosmologyParameters_,darkMatterProfileNFW_                 )
  darkMatterProfileEinasto__              =   darkMatterProfileDarkMatterOnly(.true.,cosmologyParameters_,darkMatterProfileEinasto_             )
  darkMatterProfileBurkert__              =   darkMatterProfileDarkMatterOnly(.true.,cosmologyParameters_,darkMatterProfileBurkert_             )
  darkMatterProfileZhao1996__             =   darkMatterProfileDarkMatterOnly(.true.,cosmologyParameters_,darkMatterProfileZhao1996_            )
  darkMatterProfileTruncated__            =   darkMatterProfileDarkMatterOnly(.true.,cosmologyParameters_,darkMatterProfileTruncated_           )
  darkMatterProfileTruncatedExponential__ =   darkMatterProfileDarkMatterOnly(.true.,cosmologyParameters_,darkMatterProfileTruncatedExponential_)
  node_                                   =>  treeNode                              (                 )
  basic_                                  =>  node_               %basic            (autoCreate=.true.)
  darkMatterProfile_                      =>  node_               %darkMatterProfile(autoCreate=.true.)
  call basic_            %timeSet            (cosmologyFunctions_%cosmicTime(1.0d0))
  call basic_            %timeLastIsolatedSet(cosmologyFunctions_%cosmicTime(1.0d0))
  call basic_            %massSet            (massVirial                           )
  radiusVirial    =+darkMatterHaloScale_%radiusVirial      (node_)
  timeDynamical   =+darkMatterHaloScale_%timescaleDynamical(node_)
  radiusScale     =+radiusVirial &
       &           /concentration
  coordinatesOuter=[radiusVirial,0.0d0,0.0d0]
  call darkMatterProfile_%scaleSet(radiusScale )
  call darkMatterProfile_%shapeSet(shapeProfile)
  call Unit_Tests_Begin_Group("Isothermal profile"   )
  massDistribution_       => darkMatterProfileIsothermal__%get                   (node_)
  kinematicsDistribution_ => massDistribution_            %kinematicsDistribution(     )
  select type (massDistribution_)
  class is (massDistributionSpherical)
     do i=1,size(scaleFractional)
        coordinates                                  =[scaleFractional(i)*radiusScale,0.0d0,0.0d0]
        enclosedMass                              (i)=massDistribution_      %massEnclosedBySphere                      (                                                                                    scaleFractional(i)*radiusScale                                       )
        enclosedMassNumerical                     (i)=massDistribution_      %massEnclosedBySphereNumerical             (                                                                                    scaleFractional(i)*radiusScale                                       )
        potential                                 (i)=massDistribution_      %potentialDifference                       (                                                                                    coordinates                     ,coordinatesOuter                    )
        potentialNumerical                        (i)=massDistribution_      %potentialDifferenceNumerical              (                                                                                    coordinates                     ,coordinatesOuter                    )
        velocityCircular                          (i)=massDistribution_      %rotationCurve                             (                                                                                    scaleFractional(i)*radiusScale                                       )
        velocityCircularNumerical                 (i)=massDistribution_      %rotationCurveNumerical                    (                                                                                    scaleFractional(i)*radiusScale                                       )
        radialVelocityDispersion                  (i)=kinematicsDistribution_%velocityDispersion1D                      (                                                                                    coordinates                     ,massDistribution_ ,massDistribution_)
        radialVelocityDispersionNumerical         (i)=kinematicsDistribution_%velocityDispersion1DNumerical             (                                                                                    coordinates                     ,massDistribution_ ,massDistribution_)
        kSpace                                    (i)=massDistribution_      %fourierTransform                          (radiusVirial,                                                                       scaleFractional(i)/radiusScale                                       )
        kSpaceNumerical                           (i)=massDistribution_      %fourierTransformNumerical                 (radiusVirial,                                                                       scaleFractional(i)/radiusScale                                       )
        freefallRadius                            (i)=massDistribution_      %radiusFreeFall                            (                                                                                    scaleFractional(i)*timeDynamical                                     )
        freefallRadiusNumerical                   (i)=massDistribution_      %radiusFreefallNumerical                   (                                                                                    scaleFractional(i)*timeDynamical                                     )
        freefallRadiusIncreaseRate                (i)=massDistribution_      %radiusFreefallIncreaseRate                (                                                                                    scaleFractional(i)*timeDynamical                                     )
        freefallRadiusIncreaseRateNumerical       (i)=massDistribution_      %radiusFreefallIncreaseRateNumerical       (                                                                                    scaleFractional(i)*timeDynamical                                     )
        radiusEnclosingDensity                    (i)=massDistribution_      %radiusEnclosingDensity                    (             massDistribution_%density             (coordinates                   )                                                                      )
        radiusEnclosingDensityNumerical           (i)=massDistribution_      %radiusEnclosingDensityNumerical           (             massDistribution_%density             (coordinates                   )                                                                      )
        radiusEnclosingMass                       (i)=massDistribution_      %radiusEnclosingMass                       (             massDistribution_%massEnclosedBySphere(scaleFractional(i)*radiusScale)                                                                      )
        radiusEnclosingMassNumerical              (i)=massDistribution_      %radiusEnclosingMassNumerical              (             massDistribution_%massEnclosedBySphere(scaleFractional(i)*radiusScale)                                                                      )
        radiusFromSpecificAngularMomentum         (i)=massDistribution_      %radiusFromSpecificAngularMomentum         (             massDistribution_%rotationCurve       (scaleFractional(i)*radiusScale)*scaleFractional(i)*radiusScale                                       )
        radiusFromSpecificAngularMomentumNumerical(i)=massDistribution_      %radiusFromSpecificAngularMomentumNumerical(             massDistribution_%rotationCurve       (scaleFractional(i)*radiusScale)*scaleFractional(i)*radiusScale                                       )
        densityLogSlope                           (i)=massDistribution_      %densityGradientRadial                     (                                                                                    coordinates                     ,logarithmic=.true.                  )
        densityLogSlopeNumerical                  (i)=massDistribution_      %densityGradientRadialNumerical            (                                                                                    coordinates                     ,logarithmic=.true.                  )
     end do
     call Skip  ("Energy                          , E   ","isothermal profile assumes virial equilibrium")
     call Skip  ("Radial moment                   , ℛ₁ " ,"1ˢᵗ moment diverges for isothermal profile"    )
     call Assert("Radial moment                   , ℛ₂ " ,massDistribution_%densityRadialMomentNumerical       (2.0d0,0.0d0,radiusVirial),massDistribution_%densityRadialMoment       (2.0d0,0.0d0,radiusVirial),relTol=1.0d-6              )
     call Assert("Radial moment                   , ℛ₃ " ,massDistribution_%densityRadialMomentNumerical       (3.0d0,0.0d0,radiusVirial),massDistribution_%densityRadialMoment       (3.0d0,0.0d0,radiusVirial),relTol=1.0d-6              )
     call Assert("Radius of peak circular velocity, Rmax",massDistribution_%radiusRotationCurveMaximumNumerical(                        ),massDistribution_%radiusRotationCurveMaximum(                        ),relTol=1.0d-3              )
     call Assert("Enclosed mass                   , M(r)",                  enclosedMassNumerical                                        ,                  enclosedMass                                        ,relTol=2.0d-3              )
     call Assert("Potential                       , Φ(r)",                  potentialNumerical                                           ,                  potential                                           ,relTol=1.0d-3              )
     call Assert("Circular velocity               , V(r)",                  velocityCircularNumerical                                    ,                  velocityCircular                                    ,relTol=1.0d-3              )
     call Assert("Radial velocity dispersion      , σ(r)",                  radialVelocityDispersionNumerical                            ,                  radialVelocityDispersion                            ,relTol=3.0d-3              )
     call Assert("Fourier transform               , u(k)",                  kSpaceNumerical                                              ,                  kSpace                                              ,relTol=1.0d-3,absTol=1.0d-4)
     call Assert("Freefall radius                 , r(t)",                  freefallRadiusNumerical                                      ,                  freefallRadius                                      ,relTol=1.0d-3              )
     call Assert("Freefall radius increase rate   , ̇r(t)",                 freefallRadiusIncreaseRateNumerical                          ,                  freefallRadiusIncreaseRate                          ,relTol=5.0d-3              )
     call Assert("Radius enclosing density        , r(ρ)",                  radiusEnclosingDensityNumerical                              ,                  radiusEnclosingDensity                              ,relTol=1.0d-3              )
     call Assert("Radius enclosing mass           , r(M)",                  radiusEnclosingMassNumerical                                 ,                  radiusEnclosingMass                                 ,relTol=1.0d-3              )
     call Assert("Radius-specific angular momentum, r(j)",                  radiusFromSpecificAngularMomentumNumerical                   ,                  radiusFromSpecificAngularMomentum                   ,relTol=1.0d-3              )
     call Assert("Density log gradient            , α(r)",                  densityLogSlopeNumerical                                     ,                  densityLogSlope                                     ,relTol=1.0d-3              )
  end select
  !![
  <objectDestructor name="massDistribution_"      />
  <objectDestructor name="kinematicsDistribution_"/>
  !!]
  call Unit_Tests_End_Group  (                       )
  call Unit_Tests_Begin_Group("NFW profile"          )
  massDistribution_       => darkMatterProfileNFW__%get                   (node_)
  kinematicsDistribution_ => massDistribution_     %kinematicsDistribution(     )
  select type (massDistribution_)
  class is (massDistributionSpherical)
     do i=1,size(scaleFractional)
        coordinates                                  =[scaleFractional(i)*radiusScale,0.0d0,0.0d0]
        enclosedMass                              (i)=massDistribution_      %massEnclosedBySphere                      (                                                                                    scaleFractional(i)*radiusScale                                       )
        enclosedMassNumerical                     (i)=massDistribution_      %massEnclosedBySphereNumerical             (                                                                                    scaleFractional(i)*radiusScale                                       )
        potential                                 (i)=massDistribution_      %potentialDifference                       (                                                                                    coordinates                     ,coordinatesOuter                    )
        potentialNumerical                        (i)=massDistribution_      %potentialDifferenceNumerical              (                                                                                    coordinates                     ,coordinatesOuter                    )
        velocityCircular                          (i)=massDistribution_      %rotationCurve                             (                                                                                    scaleFractional(i)*radiusScale                                       )
        velocityCircularNumerical                 (i)=massDistribution_      %rotationCurveNumerical                    (                                                                                    scaleFractional(i)*radiusScale                                       )
        radialVelocityDispersion                  (i)=kinematicsDistribution_%velocityDispersion1D                      (                                                                                    coordinates                     ,massDistribution_ ,massDistribution_)
        radialVelocityDispersionNumerical         (i)=kinematicsDistribution_%velocityDispersion1DNumerical             (                                                                                    coordinates                     ,massDistribution_ ,massDistribution_)
        kSpace                                    (i)=massDistribution_      %fourierTransform                          (radiusVirial,                                                                       scaleFractional(i)/radiusScale                                       )
        kSpaceNumerical                           (i)=massDistribution_      %fourierTransformNumerical                 (radiusVirial,                                                                       scaleFractional(i)/radiusScale                                       )
        freefallRadius                            (i)=massDistribution_      %radiusFreeFall                            (                                                                                    scaleFractional(i)*timeDynamical                                     )
        freefallRadiusNumerical                   (i)=massDistribution_      %radiusFreefallNumerical                   (                                                                                    scaleFractional(i)*timeDynamical                                     )
        freefallRadiusIncreaseRate                (i)=massDistribution_      %radiusFreefallIncreaseRate                (                                                                                    scaleFractional(i)*timeDynamical                                     )
        freefallRadiusIncreaseRateNumerical       (i)=massDistribution_      %radiusFreefallIncreaseRateNumerical       (                                                                                    scaleFractional(i)*timeDynamical                                     )
        radiusEnclosingDensity                    (i)=massDistribution_      %radiusEnclosingDensity                    (             massDistribution_%density             (coordinates                   )                                                                      )
        radiusEnclosingDensityNumerical           (i)=massDistribution_      %radiusEnclosingDensityNumerical           (             massDistribution_%density             (coordinates                   )                                                                      )
        radiusEnclosingMass                       (i)=massDistribution_      %radiusEnclosingMass                       (             massDistribution_%massEnclosedBySphere(scaleFractional(i)*radiusScale)                                                                      )
        radiusEnclosingMassNumerical              (i)=massDistribution_      %radiusEnclosingMassNumerical              (             massDistribution_%massEnclosedBySphere(scaleFractional(i)*radiusScale)                                                                      )
        radiusFromSpecificAngularMomentum         (i)=massDistribution_      %radiusFromSpecificAngularMomentum         (             massDistribution_%rotationCurve       (scaleFractional(i)*radiusScale)*scaleFractional(i)*radiusScale                                       )
        radiusFromSpecificAngularMomentumNumerical(i)=massDistribution_      %radiusFromSpecificAngularMomentumNumerical(             massDistribution_%rotationCurve       (scaleFractional(i)*radiusScale)*scaleFractional(i)*radiusScale                                       )
        densityLogSlope                           (i)=massDistribution_      %densityGradientRadial                     (                                                                                    coordinates                     ,logarithmic=.true.                  )
        densityLogSlopeNumerical                  (i)=massDistribution_      %densityGradientRadialNumerical            (                                                                                    coordinates                     ,logarithmic=.true.                  )
     end do
     call Assert("Energy                          , E   ",massDistribution_%energyNumerical                    (            radiusVirial,massDistribution_),massDistribution_%energy                    (            radiusVirial,massDistribution_),relTol=1.0d-3              )
     call Assert("Radial moment                   , ℛ₁  ",massDistribution_%densityRadialMomentNumerical       (1.0d0,0.0d0,radiusVirial                  ),massDistribution_%densityRadialMoment       (1.0d0,0.0d0,radiusVirial                  ),relTol=1.0d-6              )
     call Assert("Radial moment                   , ℛ₂  ",massDistribution_%densityRadialMomentNumerical       (2.0d0,0.0d0,radiusVirial                  ),massDistribution_%densityRadialMoment       (2.0d0,0.0d0,radiusVirial                  ),relTol=1.0d-6              )
     call Assert("Radial moment                   , ℛ₃  ",massDistribution_%densityRadialMomentNumerical       (3.0d0,0.0d0,radiusVirial                  ),massDistribution_%densityRadialMoment       (3.0d0,0.0d0,radiusVirial                  ),relTol=1.0d-6              )
     call Assert("Radius of peak circular velocity, Rmax",massDistribution_%radiusRotationCurveMaximumNumerical(                                          ),massDistribution_%radiusRotationCurveMaximum(                                          ),relTol=1.0d-3              )
     call Assert("Enclosed mass                   , M(r)",                  enclosedMassNumerical                                                          ,                  enclosedMass                                                          ,relTol=2.0d-3              )
     call Assert("Potential                       , Φ(r)",                  potentialNumerical                                                             ,                  potential                                                             ,relTol=1.0d-3              )
     call Assert("Circular velocity               , V(r)",                  velocityCircularNumerical                                                      ,                  velocityCircular                                                      ,relTol=1.0d-3              )
     call Assert("Radial velocity dispersion      , σ(r)",                  radialVelocityDispersionNumerical                                              ,                  radialVelocityDispersion                                              ,relTol=3.0d-3              )
     call Assert("Fourier transform               , u(k)",                  kSpaceNumerical                                                                ,                  kSpace                                                                ,relTol=1.0d-3,absTol=1.0d-4)
     call Assert("Freefall radius                 , r(t)",                  freefallRadiusNumerical                                                        ,                  freefallRadius                                                        ,relTol=1.0d-3              )
     call Assert("Freefall radius increase rate   , ̇r(t)",                 freefallRadiusIncreaseRateNumerical                                            ,                  freefallRadiusIncreaseRate                                            ,relTol=5.0d-3              )
     call Assert("Radius enclosing density        , r(ρ)",                  radiusEnclosingDensityNumerical                                                ,                  radiusEnclosingDensity                                                ,relTol=1.0d-3              )
     call Assert("Radius enclosing mass           , r(M)",                  radiusEnclosingMassNumerical                                                   ,                  radiusEnclosingMass                                                   ,relTol=1.0d-3              )
     call Assert("Radius-specific angular momentum, r(j)",                  radiusFromSpecificAngularMomentumNumerical                                     ,                  radiusFromSpecificAngularMomentum                                     ,relTol=1.0d-3              )
     call Assert("Density log gradient            , α(r)",                  densityLogSlopeNumerical                                                       ,                  densityLogSlope                                                       ,relTol=1.0d-3              )
  end select
  !![
  <objectDestructor name="massDistribution_"      />
  <objectDestructor name="kinematicsDistribution_"/>
  !!]
  call Unit_Tests_End_Group  (                       )
  call Unit_Tests_Begin_Group("Einasto profile"      )
  massDistribution_       => darkMatterProfileEinasto__%get                   (node_)
  kinematicsDistribution_ => massDistribution_         %kinematicsDistribution(     )
  select type (massDistribution_)
  class is (massDistributionSpherical)
     do i=1,size(scaleFractional)
        coordinates                                  =[scaleFractional(i)*radiusScale,0.0d0,0.0d0]
        enclosedMass                              (i)=massDistribution_      %massEnclosedBySphere                      (                                                                                    scaleFractional(i)*radiusScale                                       )
        enclosedMassNumerical                     (i)=massDistribution_      %massEnclosedBySphereNumerical             (                                                                                    scaleFractional(i)*radiusScale                                       )
        potential                                 (i)=massDistribution_      %potentialDifference                       (                                                                                    coordinates                     ,coordinatesOuter                    )
        potentialNumerical                        (i)=massDistribution_      %potentialDifferenceNumerical              (                                                                                    coordinates                     ,coordinatesOuter                    )
        velocityCircular                          (i)=massDistribution_      %rotationCurve                             (                                                                                    scaleFractional(i)*radiusScale                                       )
        velocityCircularNumerical                 (i)=massDistribution_      %rotationCurveNumerical                    (                                                                                    scaleFractional(i)*radiusScale                                       )
        radialVelocityDispersion                  (i)=kinematicsDistribution_%velocityDispersion1D                      (                                                                                    coordinates                     ,massDistribution_ ,massDistribution_)
        radialVelocityDispersionNumerical         (i)=kinematicsDistribution_%velocityDispersion1DNumerical             (                                                                                    coordinates                     ,massDistribution_ ,massDistribution_)
        kSpace                                    (i)=massDistribution_      %fourierTransform                          (radiusVirial,                                                                       scaleFractional(i)/radiusScale                                       )
        kSpaceNumerical                           (i)=massDistribution_      %fourierTransformNumerical                 (radiusVirial,                                                                       scaleFractional(i)/radiusScale                                       )
        freefallRadius                            (i)=massDistribution_      %radiusFreeFall                            (                                                                                    scaleFractional(i)*timeDynamical                                     )
        freefallRadiusNumerical                   (i)=massDistribution_      %radiusFreefallNumerical                   (                                                                                    scaleFractional(i)*timeDynamical                                     )
        freefallRadiusIncreaseRate                (i)=massDistribution_      %radiusFreefallIncreaseRate                (                                                                                    scaleFractional(i)*timeDynamical                                     )
        freefallRadiusIncreaseRateNumerical       (i)=massDistribution_      %radiusFreefallIncreaseRateNumerical       (                                                                                    scaleFractional(i)*timeDynamical                                     )
        radiusEnclosingDensity                    (i)=massDistribution_      %radiusEnclosingDensity                    (             massDistribution_%density             (coordinates                   )                                                                      )
        radiusEnclosingDensityNumerical           (i)=massDistribution_      %radiusEnclosingDensityNumerical           (             massDistribution_%density             (coordinates                   )                                                                      )
        radiusEnclosingMass                       (i)=massDistribution_      %radiusEnclosingMass                       (             massDistribution_%massEnclosedBySphere(scaleFractional(i)*radiusScale)                                                                      )
        radiusEnclosingMassNumerical              (i)=massDistribution_      %radiusEnclosingMassNumerical              (             massDistribution_%massEnclosedBySphere(scaleFractional(i)*radiusScale)                                                                      )
        radiusFromSpecificAngularMomentum         (i)=massDistribution_      %radiusFromSpecificAngularMomentum         (             massDistribution_%rotationCurve       (scaleFractional(i)*radiusScale)*scaleFractional(i)*radiusScale                                       )
        radiusFromSpecificAngularMomentumNumerical(i)=massDistribution_      %radiusFromSpecificAngularMomentumNumerical(             massDistribution_%rotationCurve       (scaleFractional(i)*radiusScale)*scaleFractional(i)*radiusScale                                       )
        densityLogSlope                           (i)=massDistribution_      %densityGradientRadial                     (                                                                                    coordinates                     ,logarithmic=.true.                  )
        densityLogSlopeNumerical                  (i)=massDistribution_      %densityGradientRadialNumerical            (                                                                                    coordinates                     ,logarithmic=.true.                  )
     end do
     call Assert("Energy                          , E   ",massDistribution_%energyNumerical                    (            radiusVirial,massDistribution_),massDistribution_%energy                    (            radiusVirial,massDistribution_),relTol=1.0d-3              )
     call Assert("Radial moment                   , ℛ₁  ",massDistribution_%densityRadialMomentNumerical       (1.0d0,0.0d0,radiusVirial                  ),massDistribution_%densityRadialMoment       (1.0d0,0.0d0,radiusVirial                  ),relTol=1.0d-6              )
     call Assert("Radial moment                   , ℛ₂  ",massDistribution_%densityRadialMomentNumerical       (2.0d0,0.0d0,radiusVirial                  ),massDistribution_%densityRadialMoment       (2.0d0,0.0d0,radiusVirial                  ),relTol=1.0d-6              )
     call Assert("Radial moment                   , ℛ₃  ",massDistribution_%densityRadialMomentNumerical       (3.0d0,0.0d0,radiusVirial                  ),massDistribution_%densityRadialMoment       (3.0d0,0.0d0,radiusVirial                  ),relTol=1.0d-6              )
     call Assert("Radius of peak circular velocity, Rmax",massDistribution_%radiusRotationCurveMaximumNumerical(                                          ),massDistribution_%radiusRotationCurveMaximum(                                          ),relTol=1.0d-3              )
     call Assert("Enclosed mass                   , M(r)",                  enclosedMassNumerical                                                          ,                  enclosedMass                                                          ,relTol=2.0d-3              )
     call Assert("Potential                       , Φ(r)",                  potentialNumerical                                                             ,                  potential                                                             ,relTol=1.0d-3              )
     call Assert("Circular velocity               , V(r)",                  velocityCircularNumerical                                                      ,                  velocityCircular                                                      ,relTol=1.0d-3              )
     call Assert("Radial velocity dispersion      , σ(r)",                  radialVelocityDispersionNumerical                                              ,                  radialVelocityDispersion                                              ,relTol=3.0d-3              )
     call Assert("Fourier transform               , u(k)",                  kSpaceNumerical                                                                ,                  kSpace                                                                ,relTol=1.0d-3,absTol=1.0d-4)
     call Assert("Freefall radius                 , r(t)",                  freefallRadiusNumerical                                                        ,                  freefallRadius                                                        ,relTol=1.0d-3              )
     call Assert("Freefall radius increase rate   , ̇r(t)",                 freefallRadiusIncreaseRateNumerical                                            ,                  freefallRadiusIncreaseRate                                            ,relTol=6.0d-3              )
     call Assert("Radius enclosing density        , r(ρ)",                  radiusEnclosingDensityNumerical                                                ,                  radiusEnclosingDensity                                                ,relTol=1.0d-3              )
     call Assert("Radius enclosing mass           , r(M)",                  radiusEnclosingMassNumerical                                                   ,                  radiusEnclosingMass                                                   ,relTol=1.0d-3              )
     call Assert("Radius-specific angular momentum, r(j)",                  radiusFromSpecificAngularMomentumNumerical                                     ,                  radiusFromSpecificAngularMomentum                                     ,relTol=1.0d-3              )
     call Assert("Density log gradient            , α(r)",                  densityLogSlopeNumerical                                                       ,                  densityLogSlope                                                       ,relTol=1.0d-3              )
  end select
  !![
  <objectDestructor name="massDistribution_"      />
  <objectDestructor name="kinematicsDistribution_"/>
  !!]
  call Unit_Tests_End_Group  (                       )
  call Unit_Tests_Begin_Group("Burkert profile"      )
  massDistribution_       => darkMatterProfileBurkert__%get                   (node_)
  kinematicsDistribution_ => massDistribution_         %kinematicsDistribution(     )
  select type (massDistribution_)
  class is (massDistributionSpherical)
     do i=size(scaleFractional),1,-1
        coordinates                                  =[scaleFractional(i)*radiusScale,0.0d0,0.0d0]
        enclosedMass                              (i)=massDistribution_      %massEnclosedBySphere                      (                                                                                    scaleFractional(i)*radiusScale                                       )
        enclosedMassNumerical                     (i)=massDistribution_      %massEnclosedBySphereNumerical             (                                                                                    scaleFractional(i)*radiusScale                                       )
        potential                                 (i)=massDistribution_      %potentialDifference                       (                                                                                    coordinates                     ,coordinatesOuter                    )
        potentialNumerical                        (i)=massDistribution_      %potentialDifferenceNumerical              (                                                                                    coordinates                     ,coordinatesOuter                    )
        velocityCircular                          (i)=massDistribution_      %rotationCurve                             (                                                                                    scaleFractional(i)*radiusScale                                       )
        velocityCircularNumerical                 (i)=massDistribution_      %rotationCurveNumerical                    (                                                                                    scaleFractional(i)*radiusScale                                       )
        radialVelocityDispersion                  (i)=kinematicsDistribution_%velocityDispersion1D                      (                                                                                    coordinates                     ,massDistribution_ ,massDistribution_)
        radialVelocityDispersionNumerical         (i)=kinematicsDistribution_%velocityDispersion1DNumerical             (                                                                                    coordinates                     ,massDistribution_ ,massDistribution_)
        kSpace                                    (i)=massDistribution_      %fourierTransform                          (radiusVirial,                                                                       scaleFractional(i)/radiusScale                                       )
        kSpaceNumerical                           (i)=massDistribution_      %fourierTransformNumerical                 (radiusVirial,                                                                       scaleFractional(i)/radiusScale                                       )
        freefallRadius                            (i)=massDistribution_      %radiusFreeFall                            (                                                                                    scaleFractional(i)*timeDynamical                                     )
        freefallRadiusNumerical                   (i)=massDistribution_      %radiusFreefallNumerical                   (                                                                                    scaleFractional(i)*timeDynamical                                     )
        freefallRadiusIncreaseRate                (i)=massDistribution_      %radiusFreefallIncreaseRate                (                                                                                    scaleFractional(i)*timeDynamical                                     )
        freefallRadiusIncreaseRateNumerical       (i)=massDistribution_      %radiusFreefallIncreaseRateNumerical       (                                                                                    scaleFractional(i)*timeDynamical                                     )
        radiusEnclosingDensity                    (i)=massDistribution_      %radiusEnclosingDensity                    (             massDistribution_%density             (coordinates                   )                                                                      )
        radiusEnclosingDensityNumerical           (i)=massDistribution_      %radiusEnclosingDensityNumerical           (             massDistribution_%density             (coordinates                   )                                                                      )
        radiusEnclosingMass                       (i)=massDistribution_      %radiusEnclosingMass                       (             massDistribution_%massEnclosedBySphere(scaleFractional(i)*radiusScale)                                                                      )
        radiusEnclosingMassNumerical              (i)=massDistribution_      %radiusEnclosingMassNumerical              (             massDistribution_%massEnclosedBySphere(scaleFractional(i)*radiusScale)                                                                      )
        radiusFromSpecificAngularMomentum         (i)=massDistribution_      %radiusFromSpecificAngularMomentum         (             massDistribution_%rotationCurve       (scaleFractional(i)*radiusScale)*scaleFractional(i)*radiusScale                                       )
        radiusFromSpecificAngularMomentumNumerical(i)=massDistribution_      %radiusFromSpecificAngularMomentumNumerical(             massDistribution_%rotationCurve       (scaleFractional(i)*radiusScale)*scaleFractional(i)*radiusScale                                       )
        densityLogSlope                           (i)=massDistribution_      %densityGradientRadial                     (                                                                                    coordinates                     ,logarithmic=.true.                  )
        densityLogSlopeNumerical                  (i)=massDistribution_      %densityGradientRadialNumerical            (                                                                                    coordinates                     ,logarithmic=.true.                  )
     end do
     call Assert("Energy                          , E   ",massDistribution_%energyNumerical                    (            radiusVirial,massDistribution_),massDistribution_%energy                    (            radiusVirial,massDistribution_),relTol=1.0d-3              )
     call Assert("Radial moment                   , ℛ₁  ",massDistribution_%densityRadialMomentNumerical       (1.0d0,0.0d0,radiusVirial                  ),massDistribution_%densityRadialMoment       (1.0d0,0.0d0,radiusVirial                  ),relTol=1.0d-6              )
     call Assert("Radial moment                   , ℛ₂  ",massDistribution_%densityRadialMomentNumerical       (2.0d0,0.0d0,radiusVirial                  ),massDistribution_%densityRadialMoment       (2.0d0,0.0d0,radiusVirial                  ),relTol=1.0d-6              )
     call Assert("Radial moment                   , ℛ₃  ",massDistribution_%densityRadialMomentNumerical       (3.0d0,0.0d0,radiusVirial                  ),massDistribution_%densityRadialMoment       (3.0d0,0.0d0,radiusVirial                  ),relTol=1.0d-6              )
     call Assert("Radius of peak circular velocity, Rmax",massDistribution_%radiusRotationCurveMaximumNumerical(                                          ),massDistribution_%radiusRotationCurveMaximum(                                          ),relTol=1.0d-3              )
     call Assert("Enclosed mass                   , M(r)",                  enclosedMassNumerical                                                          ,                  enclosedMass                                                          ,relTol=2.5d-3              )
     call Assert("Potential                       , Φ(r)",                  potentialNumerical                                                             ,                  potential                                                             ,relTol=1.0d-3              )
     call Assert("Circular velocity               , V(r)",                  velocityCircularNumerical                                                      ,                  velocityCircular                                                      ,relTol=1.0d-3              )
     call Assert("Radial velocity dispersion      , σ(r)",                  radialVelocityDispersionNumerical                                              ,                  radialVelocityDispersion                                              ,relTol=3.0d-3              )
     call Assert("Fourier transform               , u(k)",                  kSpaceNumerical                                                                ,                  kSpace                                                                ,relTol=1.0d-3,absTol=1.0d-4)
     call Assert("Freefall radius                 , r(t)",                  freefallRadiusNumerical                                                        ,                  freefallRadius                                                        ,relTol=1.0d-3              )
     call Assert("Freefall radius increase rate   , ̇r(t)",                  freefallRadiusIncreaseRateNumerical                                            ,                  freefallRadiusIncreaseRate                                            ,relTol=5.0d-3              )
     call Assert("Radius enclosing density        , r(ρ)",                  radiusEnclosingDensityNumerical                                                ,                  radiusEnclosingDensity                                                ,relTol=1.0d-3              )
     call Assert("Radius enclosing mass           , r(M)",                  radiusEnclosingMassNumerical                                                   ,                  radiusEnclosingMass                                                   ,relTol=1.0d-3              )
     call Assert("Radius-specific angular momentum, r(j)",                  radiusFromSpecificAngularMomentumNumerical                                     ,                  radiusFromSpecificAngularMomentum                                     ,relTol=1.0d-3              )
     call Assert("Density log gradient            , α(r)",                  densityLogSlopeNumerical                                                       ,                  densityLogSlope                                                       ,relTol=1.0d-3              )
  end select
  !![
  <objectDestructor name="massDistribution_"      />
  <objectDestructor name="kinematicsDistribution_"/>
  !!]
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Begin_Group("Truncated NFW profile")
  massDistribution_       => darkMatterProfileTruncated__%get                   (node_)
  kinematicsDistribution_ => massDistribution_           %kinematicsDistribution(     )
  select type (massDistribution_)
  class is (massDistributionSpherical)
     do i=1,size(scaleFractional)
        coordinates                                  =[scaleFractional(i)*radiusScale,0.0d0,0.0d0]
        radialVelocityDispersion                  (i)=kinematicsDistribution_%velocityDispersion1D                      (                                                                                    coordinates                     ,massDistribution_ ,massDistribution_)
        radialVelocityDispersionNumerical         (i)=kinematicsDistribution_%velocityDispersion1DNumerical             (                                                                                    coordinates                     ,massDistribution_ ,massDistribution_)
        densityLogSlope                           (i)=massDistribution_      %densityGradientRadial                     (                                                                                    coordinates                     ,logarithmic=.true.                  )
        densityLogSlopeNumerical                  (i)=massDistribution_      %densityGradientRadialNumerical            (                                                                                    coordinates                     ,logarithmic=.true.                  )
     end do
     call Skip  ("Energy                          , E   ","implemented numerically"                                                                                                                                                                                       )
     call Skip  ("Energy growth rate              , ̇E   ","implemented numerically"                                                                                                                                                                                       )
     call Skip  ("Radial moment                   , ℛ₁  ","implemented numerically"                                                                                                                                                                                       )
     call Skip  ("Radial moment                   , ℛ₂  ","implemented numerically"                                                                                                                                                                                       )
     call Skip  ("Radial moment                   , ℛ₃  ","implemented numerically"                                                                                                                                                                                       )
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
  end select
  !![
  <objectDestructor name="massDistribution_"      />
  <objectDestructor name="kinematicsDistribution_"/>
  !!]
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Begin_Group("Exponentially-truncated NFW profile")
  massDistribution_       => darkMatterProfileTruncatedExponential__%get                   (node_)
  kinematicsDistribution_ => massDistribution_                      %kinematicsDistribution(     )
  select type (massDistribution_)
  class is (massDistributionSpherical)
     do i=1,size(scaleFractional)
        coordinates                                  =[scaleFractional(i)*radiusScale,0.0d0,0.0d0]
        radialVelocityDispersion                  (i)=kinematicsDistribution_%velocityDispersion1D                      (                                                                                    coordinates                     ,massDistribution_ ,massDistribution_)
        radialVelocityDispersionNumerical         (i)=kinematicsDistribution_%velocityDispersion1DNumerical             (                                                                                    coordinates                     ,massDistribution_ ,massDistribution_)
        densityLogSlope                           (i)=massDistribution_      %densityGradientRadial                     (                                                                                    coordinates                     ,logarithmic=.true.                  )
        densityLogSlopeNumerical                  (i)=massDistribution_      %densityGradientRadialNumerical            (                                                                                    coordinates                     ,logarithmic=.true.                  )
     end do
     call Skip  ("Energy                          , E   ","implemented numerically"                                                                                                                                                                                       )
     call Skip  ("Energy growth rate              , ̇E   ","implemented numerically"                                                                                                                                                                                       )
     call Skip  ("Radial moment                   , ℛ₁  ","implemented numerically"                                                                                                                                                                                       )
     call Skip  ("Radial moment                   , ℛ₂  ","implemented numerically"                                                                                                                                                                                       )
     call Skip  ("Radial moment                   , ℛ₃  ","implemented numerically"                                                                                                                                                                                       )
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
  end select
  !![
  <objectDestructor name="massDistribution_"      />
  <objectDestructor name="kinematicsDistribution_"/>
  !!]
  call Unit_Tests_End_Group               ()
  ! The Zhao1996 profile is a generalized NFW-type profile. So, compare to the NFW profile.
  call Unit_Tests_Begin_Group("Zhao1996 profile"          )
  massDistribution_       => darkMatterProfileZhao1996__%get                   (node_)
  kinematicsDistribution_ => massDistribution_          %kinematicsDistribution(     )
  select type (massDistribution_)
  class is (massDistributionSpherical)
     do i=1,size(scaleFractional)
        coordinates                                  =[scaleFractional(i)*radiusScale,0.0d0,0.0d0]
        enclosedMass                              (i)=massDistribution_      %massEnclosedBySphere                      (                                                                                    scaleFractional(i)*radiusScale                                       )
        enclosedMassNumerical                     (i)=massDistribution_      %massEnclosedBySphereNumerical             (                                                                                    scaleFractional(i)*radiusScale                                       )
        potential                                 (i)=massDistribution_      %potentialDifference                       (                                                                                    coordinates                     ,coordinatesOuter                    )
        potentialNumerical                        (i)=massDistribution_      %potentialDifferenceNumerical              (                                                                                    coordinates                     ,coordinatesOuter                    )
        velocityCircular                          (i)=massDistribution_      %rotationCurve                             (                                                                                    scaleFractional(i)*radiusScale                                       )
        velocityCircularNumerical                 (i)=massDistribution_      %rotationCurveNumerical                    (                                                                                    scaleFractional(i)*radiusScale                                       )
        radialVelocityDispersion                  (i)=kinematicsDistribution_%velocityDispersion1D                      (                                                                                    coordinates                     ,massDistribution_ ,massDistribution_)
        radialVelocityDispersionNumerical         (i)=kinematicsDistribution_%velocityDispersion1DNumerical             (                                                                                    coordinates                     ,massDistribution_ ,massDistribution_)
        kSpace                                    (i)=massDistribution_      %fourierTransform                          (radiusVirial,                                                                       scaleFractional(i)/radiusScale                                       )
        kSpaceNumerical                           (i)=massDistribution_      %fourierTransformNumerical                 (radiusVirial,                                                                       scaleFractional(i)/radiusScale                                       )
        freefallRadius                            (i)=massDistribution_      %radiusFreeFall                            (                                                                                    scaleFractional(i)*timeDynamical                                     )
        freefallRadiusNumerical                   (i)=massDistribution_      %radiusFreefallNumerical                   (                                                                                    scaleFractional(i)*timeDynamical                                     )
        freefallRadiusIncreaseRate                (i)=massDistribution_      %radiusFreefallIncreaseRate                (                                                                                    scaleFractional(i)*timeDynamical                                     )
        freefallRadiusIncreaseRateNumerical       (i)=massDistribution_      %radiusFreefallIncreaseRateNumerical       (                                                                                    scaleFractional(i)*timeDynamical                                     )
        radiusEnclosingDensity                    (i)=massDistribution_      %radiusEnclosingDensity                    (             massDistribution_%density             (coordinates                   )                                                                      )
        radiusEnclosingDensityNumerical           (i)=massDistribution_      %radiusEnclosingDensityNumerical           (             massDistribution_%density             (coordinates                   )                                                                      )
        radiusEnclosingMass                       (i)=massDistribution_      %radiusEnclosingMass                       (             massDistribution_%massEnclosedBySphere(scaleFractional(i)*radiusScale)                                                                      )
        radiusEnclosingMassNumerical              (i)=massDistribution_      %radiusEnclosingMassNumerical              (             massDistribution_%massEnclosedBySphere(scaleFractional(i)*radiusScale)                                                                      )
        radiusFromSpecificAngularMomentum         (i)=massDistribution_      %radiusFromSpecificAngularMomentum         (             massDistribution_%rotationCurve       (scaleFractional(i)*radiusScale)*scaleFractional(i)*radiusScale                                       )
        radiusFromSpecificAngularMomentumNumerical(i)=massDistribution_      %radiusFromSpecificAngularMomentumNumerical(             massDistribution_%rotationCurve       (scaleFractional(i)*radiusScale)*scaleFractional(i)*radiusScale                                       )
        densityLogSlope                           (i)=massDistribution_      %densityGradientRadial                     (                                                                                    coordinates                     ,logarithmic=.true.                  )
        densityLogSlopeNumerical                  (i)=massDistribution_      %densityGradientRadialNumerical            (                                                                                    coordinates                     ,logarithmic=.true.                  )
     end do
     call Assert("Energy                          , E   ",massDistribution_%energyNumerical                    (            radiusVirial,massDistribution_),massDistribution_%energy                    (            radiusVirial,massDistribution_),relTol=1.0d-3              )
     call Assert("Radial moment                   , ℛ₁  ",massDistribution_%densityRadialMomentNumerical       (1.0d0,0.0d0,radiusVirial                  ),massDistribution_%densityRadialMoment       (1.0d0,0.0d0,radiusVirial                  ),relTol=1.0d-6              )
     call Assert("Radial moment                   , ℛ₂  ",massDistribution_%densityRadialMomentNumerical       (2.0d0,0.0d0,radiusVirial                  ),massDistribution_%densityRadialMoment       (2.0d0,0.0d0,radiusVirial                  ),relTol=1.0d-6              )
     call Assert("Radial moment                   , ℛ₃  ",massDistribution_%densityRadialMomentNumerical       (3.0d0,0.0d0,radiusVirial                  ),massDistribution_%densityRadialMoment       (3.0d0,0.0d0,radiusVirial                  ),relTol=1.0d-6              )
     call Assert("Radius of peak circular velocity, Rmax",massDistribution_%radiusRotationCurveMaximumNumerical(                                          ),massDistribution_%radiusRotationCurveMaximum(                                          ),relTol=1.0d-3              )
     call Assert("Enclosed mass                   , M(r)",                  enclosedMassNumerical                                                          ,                  enclosedMass                                                          ,relTol=2.5d-3              )
     call Assert("Potential                       , Φ(r)",                  potentialNumerical                                                             ,                  potential                                                             ,relTol=1.0d-3              )
     call Assert("Circular velocity               , V(r)",                  velocityCircularNumerical                                                      ,                  velocityCircular                                                      ,relTol=1.0d-3              )
     call Assert("Radial velocity dispersion      , σ(r)",                  radialVelocityDispersionNumerical                                              ,                  radialVelocityDispersion                                              ,relTol=3.0d-3              )
     call Assert("Fourier transform               , u(k)",                  kSpaceNumerical                                                                ,                  kSpace                                                                ,relTol=1.0d-3,absTol=1.0d-4)
     call Assert("Freefall radius                 , r(t)",                  freefallRadiusNumerical                                                        ,                  freefallRadius                                                        ,relTol=1.0d-3              )
     call Assert("Freefall radius increase rate   , ̇r(t)",                 freefallRadiusIncreaseRateNumerical                                            ,                  freefallRadiusIncreaseRate                                            ,relTol=5.0d-3              )
     call Assert("Radius enclosing density        , r(ρ)",                  radiusEnclosingDensityNumerical                                                ,                  radiusEnclosingDensity                                                ,relTol=1.0d-3              )
     call Assert("Radius enclosing mass           , r(M)",                  radiusEnclosingMassNumerical                                                   ,                  radiusEnclosingMass                                                   ,relTol=1.0d-3              )
     call Assert("Radius-specific angular momentum, r(j)",                  radiusFromSpecificAngularMomentumNumerical                                     ,                  radiusFromSpecificAngularMomentum                                     ,relTol=1.0d-3              )
     call Assert("Density log gradient            , α(r)",                  densityLogSlopeNumerical                                                       ,                  densityLogSlope                                                       ,relTol=1.0d-3              )
  end select
  !![
  <objectDestructor name="massDistribution_"      />
  <objectDestructor name="kinematicsDistribution_"/>
  !!]
  call Unit_Tests_End_Group               ()
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
end program Test_Dark_Matter_Profiles_Generic
