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
Contains a program which tests dark matter profiles.
!!}

program Test_Dark_Matter_Profiles
  !!{
  Tests dark matter profiles.
  !!}
  use :: Calculations_Resets             , only : Calculations_Reset
  use :: Coordinates                     , only : coordinateSpherical                                           , assignment(=)
  use :: Cosmology_Parameters            , only : cosmologyParametersSimple
  use :: Cosmology_Functions             , only : cosmologyFunctionsMatterLambda
  use :: Dark_Matter_Particles           , only : darkMatterParticleSelfInteractingDarkMatter                   , darkMatterParticleCDM
  use :: Dark_Matter_Halo_Scales         , only : darkMatterHaloScaleVirialDensityContrastDefinition
  use :: Dark_Matter_Profiles_DMO        , only : darkMatterProfileDMOBurkert                                   , darkMatterProfileDMONFW             , darkMatterProfileDMOFiniteResolution, darkMatterProfileDMOSIDMCoreNFW, &
       &                                          darkMatterProfileDMOSIDMIsothermal                            , darkMatterProfileDMOZhao1996
  use :: Dark_Matter_Profiles            , only : darkMatterProfileSIDMIsothermal                               , darkMatterProfileAdiabaticGnedin2004
  use :: Virial_Density_Contrast         , only : virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt, virialDensityContrastFixed          , fixedDensityTypeCritical
  use :: Events_Hooks                    , only : eventsHooksInitialize
  use :: Functions_Global_Utilities      , only : Functions_Global_Set
  use :: Display                         , only : displayVerbositySet                                           , verbosityLevelStandard
  use :: Mass_Distributions              , only : massDistributionClass                                         , massDistributionSpherical           , kinematicsDistributionClass          , nonAnalyticSolversNumerical   , &
       &                                          massDistributionSphericalSIDM
  use :: Galacticus_Nodes                , only : nodeClassHierarchyFinalize                                    , nodeClassHierarchyInitialize        , nodeComponentBasic                   , nodeComponentDarkMatterProfile, &
          &                                       treeNode                                                      , nodeComponentSpheroid
  use :: Numerical_Constants_Math        , only : Pi
  use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal                                , MpcPerKmPerSToGyr
  use :: Input_Parameters                , only : inputParameters
  use :: ISO_Varying_String              , only : varying_string                                                , assignment(=)                       , var_str
  use :: Node_Components                 , only : Node_Components_Initialize                                    , Node_Components_Thread_Initialize   , Node_Components_Thread_Uninitialize  , Node_Components_Uninitialize
  use :: Unit_Tests                      , only : Assert                                                        , Unit_Tests_Begin_Group              , Unit_Tests_End_Group                 , Unit_Tests_Finish
  implicit none
  type            (treeNode                                                      ), pointer      :: node                                                                                                , &
       &                                                                                            nodePippin                                                                                          , &
       &                                                                                            nodeJiang
  class           (nodeComponentBasic                                            ), pointer      :: basic
  class           (nodeComponentDarkMatterProfile                                ), pointer      :: dmProfile
  class           (nodeComponentSpheroid                                         ), pointer      :: spheroid
  class           (massDistributionClass                                         ), pointer      :: massDistribution_
  class           (kinematicsDistributionClass                                   ), pointer      :: kinematicsDistribution_
  double precision                                                                , parameter    :: concentration                       = 8.0d+0                                                        , &
       &                                                                                            massVirial                          = 1.0d+0                                                        , &
       &                                                                                            massSmall                           = 1.0d-7                                                        , &
       ! Mass and concentration of Pippin halos (Jiang et al. 2022).
       &                                                                                            concentrationPippin                 =15.8d0                                                         , &
       &                                                                                            massVirialPippin                    =10.0d0**9.89d0                                                 , &
       ! Mass and concentration of example halo from Fangzhou Jiang (private communication).
       &                                                                                            concentrationJiang                  =15.0d0                                                         , &
       &                                                                                            timeJiang                           =10.0d0                                                         , &
       &                                                                                            massVirialJiang                     =1.0d11                                                         , &
       &                                                                                            fractionMassBaryonicJiang           =1.0d-2                                                         , &
       &                                                                                            fractionRadiusHalfMassJiang         =2.0d-2                                                         , &
       &                                                                                            radiusHalfMassDimensionlessHernquist=1.0d0/(sqrt(2.0d0)-1.0d0)
  double precision                                                                , dimension(7) :: radius                              =[0.125000d0,0.25000d0,0.50000d0,1.00000d0,2.00000d0,4.00000d0,8.00000d0]
  double precision                                                                , dimension(7) :: timeFreefall
  double precision                                                                , dimension(7) :: mass                                                                                                , &
       &                                                                                            density                                                                                             , &
       &                                                                                            radiusRecoveredFromMass                                                                             , &
       &                                                                                            radiusRecoveredFromDensity                                                                          , &
       &                                                                                            radiusRecoveredFromTimeFreefall                                                                     , &
       &                                                                                            radiusRecoveredFromSpecificAngularMomentum                                                          , &
       &                                                                                            potential                                                                                           , &
       &                                                                                            fourier                                                                                             , &
       &                                                                                            radialVelocityDispersion                                                                            , &
       &                                                                                            radialVelocityDispersionSeriesExpansion
  type            (darkMatterParticleCDM                                         ), pointer      :: darkMatterParticleCDM_
  type            (darkMatterParticleSelfInteractingDarkMatter                   ), pointer      :: darkMatterParticleSelfInteractingDarkMatter_                                                        , &
       &                                                                                            darkMatterParticleSelfInteractingDarkMatterJiang_
  type            (darkMatterProfileDMOBurkert                                   ), pointer      :: darkMatterProfileDMOBurkert_
  type            (darkMatterProfileDMOZhao1996                                  ), pointer      :: darkMatterProfileDMOZhao1996_
  type            (darkMatterProfileDMONFW                                       ), pointer      :: darkMatterProfileDMONFW_                                                                            , &
       &                                                                                            darkMatterProfileDMONFWPippin_
  type            (darkMatterProfileDMONFW                                       ), pointer      :: darkMatterProfileDMONFWSeriesExpansion_
  type            (darkMatterProfileDMOFiniteResolution                          ), pointer      :: darkMatterProfileDMOFiniteResolution_
  type            (darkMatterProfileDMOSIDMCoreNFW                               ), pointer      :: darkMatterProfileDMOSIDMCoreNFW_
  type            (darkMatterProfileDMOSIDMIsothermal                            ), pointer      :: darkMatterProfileDMOSIDMIsothermal_
  type            (darkMatterProfileAdiabaticGnedin2004                          ), pointer      :: darkMatterProfileAdiabaticPippin_
  type            (darkMatterProfileSIDMIsothermal                               ), pointer      :: darkMatterProfileSIDMIsothermal_
  type            (cosmologyParametersSimple                                     ), pointer      :: cosmologyParameters_                                                                                , &
       &                                                                                            cosmologyParametersPippin_
  type            (cosmologyFunctionsMatterLambda                                ), pointer      :: cosmologyFunctions_                                                                                 , &
       &                                                                                            cosmologyFunctionsPippin_
  type            (darkMatterHaloScaleVirialDensityContrastDefinition            ), pointer      :: darkMatterHaloScale_                                                                                , &
       &                                                                                            darkMatterHaloScalePippin_
  type            (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt), pointer      :: virialDensityContrast_
  type            (virialDensityContrastFixed                                    ), pointer      :: virialDensityContrastPippin_
  type            (coordinateSpherical                                           )               :: coordinates
  type            (inputParameters                                               )               :: parameters
  integer                                                                                        :: i
  double precision                                                                               :: radiusScale                                                                                         , &
       &                                                                                            timeScale                                                                                           , &
       &                                                                                            radiusVelocityMaximum                                                                               , &
       &                                                                                            velocityMaximum                                                                                     , &
       &                                                                                            velocityMaximumIndirect                                                                             , &
       &                                                                                            radiusVirial                                                                                        , &
       &                                                                                            energyPotential                                                                                     , &
       &                                                                                            energyKinetic                                                                                       , &
       &                                                                                            energyPotentialNumerical                                                                            , &
       &                                                                                            energyKineticNumerical                                                                              , &
       &                                                                                            densityNormalization                                                                                , &
       &                                                                                            radiusSmall


  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Initialize event hooks and global functions.
  call eventsHooksInitialize()
  call Functions_Global_Set ()
  ! Initialize node components.
  parameters=inputParameters(var_str('testSuite/parameters/darkMatterProfiles.xml'))
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  ! Build required objects.
  allocate(cosmologyParameters_                             )
  allocate(cosmologyFunctions_                              )
  allocate(virialDensityContrast_                           )
  allocate(darkMatterHaloScale_                             )
  allocate(darkMatterProfileDMOBurkert_                     )
  allocate(darkMatterProfileDMOZhao1996_                    )
  allocate(darkMatterProfileDMONFW_                         )
  allocate(darkMatterProfileDMONFWPippin_                   )
  allocate(darkMatterProfileDMONFWSeriesExpansion_          )
  allocate(darkMatterProfileDMOFiniteResolution_            )
  allocate(darkMatterProfileDMOSIDMCoreNFW_                 )
  allocate(darkMatterProfileDMOSIDMIsothermal_              )
  allocate(darkMatterProfileSIDMIsothermal_                 )
  allocate(darkMatterProfileAdiabaticPippin_                )
  allocate(darkMatterParticleSelfInteractingDarkMatter_     )
  allocate(darkMatterParticleSelfInteractingDarkMatterJiang_)
  allocate(darkMatterParticleCDM_                           )
  allocate(cosmologyParametersPippin_                       )
  allocate(cosmologyFunctionsPippin_                        )
  allocate(virialDensityContrastPippin_                     )
  allocate(darkMatterHaloScalePippin_                       )
  !![
  <referenceConstruct object="cosmologyParameters_"        >
   <constructor>
    cosmologyParametersSimple                                     (                                               &amp;
     &amp;                                                         OmegaMatter           = 0.2815d0             , &amp;
     &amp;                                                         OmegaBaryon           = 0.0465d0             , &amp;
     &amp;                                                         OmegaDarkEnergy       = 0.7185d0             , &amp;
     &amp;                                                         temperatureCMB        = 2.7800d0             , &amp;
     &amp;                                                         HubbleConstant        =69.3000d0               &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctions_"         >
   <constructor>
    cosmologyFunctionsMatterLambda                                (                                               &amp;
     &amp;                                                         cosmologyParameters_  =cosmologyParameters_    &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="virialDensityContrast_"      >
   <constructor>
    virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt(                                               &amp;
     &amp;                                                         tableStore            =.true.                , &amp;
     &amp;                                                         cosmologyFunctions_   =cosmologyFunctions_     &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterHaloScale_"        >
   <constructor>
    darkMatterHaloScaleVirialDensityContrastDefinition            (                                               &amp;
     &amp;                                                         cosmologyParameters_  =cosmologyParameters_  , &amp;
     &amp;                                                         cosmologyFunctions_   =cosmologyFunctions_   , &amp;
     &amp;                                                         virialDensityContrast_=virialDensityContrast_  &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  !!]
  ! Create a node.
  node      => treeNode                  (                 )
  ! Create components.
  basic     => node    %basic            (autoCreate=.true.)
  dmProfile => node    %darkMatterProfile(autoCreate=.true.)
  ! Set properties.
  call basic%timeSet(cosmologyFunctions_%cosmicTime(1.0d0))
  call basic%massSet(massVirial                           )
  ! Compute scale radius.
  radiusScale         =+darkMatterHaloScale_%radiusVirial(node)     &
       &               /concentration
  densityNormalization=+massVirial                                  &
       &               /radiusScale**3
  timeScale            =+1.0d0/sqrt(                                &
         &                          +gravitationalConstant_internal &
         &                          *densityNormalization           &
         &                         )                                &
         &        *MpcPerKmPerSToGyr
  call dmProfile%scaleSet(radiusScale)
  ! Build dark matter profiles.
  !![
  <referenceConstruct object="darkMatterProfileDMOBurkert_"           >
   <constructor>
    darkMatterProfileDMOBurkert         (                                                                  &amp;
     &amp;                               darkMatterHaloScale_                =darkMatterHaloScale_         &amp;
     &amp;                              )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfileDMONFW_"               >
   <constructor>
    darkMatterProfileDMONFW             (                                                                  &amp;
     &amp;                               velocityDispersionUseSeriesExpansion=.false.                    , &amp;
     &amp;                               darkMatterHaloScale_                =darkMatterHaloScale_         &amp;
     &amp;                              )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfileDMONFWSeriesExpansion_">
   <constructor>
    darkMatterProfileDMONFW             (                                                                  &amp;
     &amp;                               velocityDispersionUseSeriesExpansion=.true.                     , &amp;
     &amp;                               darkMatterHaloScale_                =darkMatterHaloScale_         &amp;
     &amp;                              )
   </constructor>
  </referenceConstruct>  
  <referenceConstruct object="darkMatterProfileDMOFiniteResolution_"  >
   <constructor>
    darkMatterProfileDMOFiniteResolution(                                                                  &amp;
     &amp;                               lengthResolution                    =0.5d0*radiusScale          , &amp;
     &amp;                               massResolution                      =0.0d0                      , &amp;
     &amp;                               resolutionIsComoving                =.false.                    , &amp;
     &amp;                               nonAnalyticSolver                   =nonAnalyticSolversNumerical, &amp;
     &amp;                               darkMatterProfileDMO_               =darkMatterProfileDMONFW_   , &amp;
     &amp;                               cosmologyFunctions_                 =cosmologyFunctions_          &amp;
     &amp;                              )
   </constructor>
  </referenceConstruct>
  !!]
  ! Begin unit tests.
  call Unit_Tests_Begin_Group('Dark matter profiles')
  ! Test Zhao1996 profile.
  call Unit_Tests_Begin_Group('Zhao1996 profile')
  !! Special case: NFW: (α,β,γ) = (1,3,1).
  call Unit_Tests_Begin_Group('(α,β,γ) = (1,3,1)')
  !![
  <referenceConstruct object="darkMatterProfileDMOZhao1996_">
   <constructor>
    darkMatterProfileDMOZhao1996 (                                          &amp;
     &amp;                       alpha               =1.0d0               , &amp;
     &amp;                       beta                =3.0d0               , &amp;
     &amp;                       gamma               =1.0d0               , &amp;
     &amp;                       darkMatterHaloScale_=darkMatterHaloScale_  &amp;
     &amp;                      )
   </constructor>
  </referenceConstruct>
  !!]
  massDistribution_       => darkMatterProfileDMOZhao1996_%get                   (node)
  kinematicsDistribution_ => massDistribution_            %kinematicsDistribution(    )
  radiusVirial            =  darkMatterHaloScale_         %radiusVirial          (node)
  timeFreefall            =  [0.864113d0,1.29807d0,2.04324d0,3.44421d0,6.32148d0,1.26727d1,2.74642d1]
  select type (massDistribution_)
  class is (massDistributionSpherical)
     do i=1,7
        coordinates=[radiusScale*radius(i),0.0d0,0.0d0]
        mass                                      (i)=massDistribution_      %massEnclosedBySphere             (                                                         radiusScale   *radius      (i)                                        )
        density                                   (i)=massDistribution_      %density                          (                                                                        coordinates                                            )*radiusScale**3
        radiusRecoveredFromMass                   (i)=massDistribution_      %radiusEnclosingMass              (                                                 mass(i)                                                                       )/radiusScale
        radiusRecoveredFromDensity                (i)=massDistribution_      %radiusEnclosingDensity           (                                  3.0d0/4.0d0/Pi*mass(i)/radiusScale**3/radius      (i)**3                                     )/radiusScale
        radiusRecoveredFromSpecificAngularMomentum(i)=massDistribution_      %radiusFromSpecificAngularMomentum(             sqrt(gravitationalConstant_internal*mass(i)*radiusScale   *radius      (i)   )                                    )/radiusScale
        radiusRecoveredFromTimeFreefall           (i)=massDistribution_      %radiusFreefall                   (                                                         timeScale     *timeFreefall(i)                                        )/radiusScale
        potential                                 (i)=massDistribution_      %potential                        (                                                                        coordinates                                            )*radiusScale/gravitationalConstant_internal
        fourier                                   (i)=massDistribution_      %fourierTransform                 (radiusVirial,1.0d0                                      /radiusScale   /radius      (i)                                        )
        radialVelocityDispersion                  (i)=kinematicsDistribution_%velocityDispersion1D             (                                                                        coordinates        ,massDistribution_,massDistribution_)
     end do
     radiusSmall             =massDistribution_%radiusEnclosingMass         (massSmall                              )
     radiusVelocityMaximum   =massDistribution_%radiusRotationCurveMaximum  (                                       )
     velocityMaximum         =massDistribution_%velocityRotationCurveMaximum(                                       )
     velocityMaximumIndirect =massDistribution_%rotationCurve               (radiusVelocityMaximum                  )
     energyPotential         =massDistribution_%energyPotential             (radiusVirial                           )
     energyPotentialNumerical=massDistribution_%energyPotentialNumerical    (radiusVirial                           )
     energyKinetic           =massDistribution_%energyKinetic               (radiusVirial         ,massDistribution_)
     energyKineticNumerical  =massDistribution_%energyKineticNumerical      (radiusVirial         ,massDistribution_)
  end select
  !![
  <objectDestructor name="massDistribution_"            />
  <objectDestructor name="kinematicsDistribution_"      />
  <objectDestructor name="darkMatterProfileDMOZhao1996_"/>
  !!]
  ! Radial velocity dispersion in units of virial velocity.
  radialVelocityDispersion=radialVelocityDispersion/darkMatterHaloScale_%velocityVirial(node)
  call Assert(                         &
       &      'enclosed mass'        , &
       &      mass                   , &
       &      [                        &
       &       5.099550982355504d-3  , &
       &       1.768930674181593d-2  , &
       &       5.513246746363203d-2  , &
       &       1.476281525188409d-1  , &
       &       3.301489257042704d-1  , &
       &       6.186775455118112d-1  , &
       &       1.000000000000000d+0    &
       &      ]                      , &
       &      relTol=1.0d-6            &
       &     )  
  call Assert(                         &
       &      'radius enclosing mass', &
       &      radiusRecoveredFromMass, &
       &      radius                 , &
       &      relTol=1.0d-6            &
       &     )
  call Assert(                                          &
       &      'radius enclosing mass (at small radii)', &
       &      radiusSmall                             , &
       &      1.6764798688849444d-9                   , &
       &      relTol=1.0d-5                             &
       &     )
  call Assert(                         &
       &      'density'              , &
       &      density                , &
       &      [                        &
       &       3.844641857939078d-1  , &
       &       1.557079952465327d-1  , &
       &       5.406527612726829d-2  , &
       &       1.520585891079421d-2  , &
       &       3.379079757954268d-3  , &
       &       6.082343564317684d-4  , &
       &       9.386332660984080d-5    &
       &      ]                      , &
       &      relTol=1.0d-6            &
       &     )
  call Assert(                            &
       &      'radius enclosing density', &
       &      radiusRecoveredFromDensity, &
       &      radius                    , &
       &      relTol=2.0d-4               &
       &     )
  call Assert(                                            &
       &      'radius from specific angular momentum'   , &
       &      radiusRecoveredFromSpecificAngularMomentum, &
       &      radius                                    , &
       &      relTol=2.0d-4                               &
       &     )
  call Assert(                                            &
       &      'radius from freefall time'               , &
       &      radiusRecoveredFromTimeFreefall           , &
       &      radius                                    , &
       &      relTol=2.0d-4                               &
       &     )
  call Assert(                       &
       &      'potential'          , &
       &      +potential             &
       &      -potential(1)        , &
       &      [                      &
       &       4.412912928902085d-2, &
       &       8.210873989889300d-2, &
       &       1.445116765163305d-1, &
       &       2.345367646465509d-1, &
       &       3.444787600350539d-1, &
       &       4.567944810866741d-1, &
       &       5.544042971829187d-1  &
       &      ]                      &
       &      -4.412912928902085d-2, &
       &      relTol=1.0d-6          &
       &     )
  call Assert(                         &
       &      'fourier'              , &
       &      fourier                , &
       &      [                        &
       &       1.099052711906094d-2  , &
       &       3.746284433831412d-2  , &
       &       1.127728075593513d-1  , &
       &       2.619030036621630d-1  , &
       &       5.434560723827092d-1  , &
       &       8.448608065763800d-1  , &
       &       9.579597044271230d-1    &
       &      ]                      , &
       &      relTol=1.0d-6            &
       &     )
  call Assert(                         &
       &     'peak rotation velocity', &
       &     velocityMaximum         , &
       &     velocityMaximumIndirect , &
       &     relTol=1.0d-6             &
       &     )
  call Assert(                         &
       &     'potential energy'      , &
       &     energyPotential         , &
       &     energyPotentialNumerical, &
       &     relTol=1.0d-6             &
       &     )
  call Assert(                         &
       &     'kinetic energy'        , &
       &     energyKinetic           , &
       &     energyKineticNumerical  , &
       &     relTol=1.0d-6             &
       &     )
  call Assert(                              &
       &      'radial velocity dispersion', &
       &      radialVelocityDispersion    , &
       &      [                             &
       &       6.285734096346791d-1       , &
       &       7.048421272327700d-1       , &
       &       7.510776824829392d-1       , &
       &       7.558757679348769d-1       , &
       &       7.172187600402157d-1       , &
       &       6.441447075213688d-1       , &
       &       5.520559866965048d-1         &
       &      ]                           , &
       &      relTol=1.0d-6                 &
       &     )
  call Unit_Tests_End_Group()  
  !! Special case: cored NFW: (α,β,γ) = (1,3,0).
  call Unit_Tests_Begin_Group('(α,β,γ) = (1,3,0)')
  allocate(darkMatterProfileDMOZhao1996_)
  !![
   <referenceConstruct object="darkMatterProfileDMOZhao1996_">
   <constructor>
    darkMatterProfileDMOZhao1996 (                                          &amp;
     &amp;                       alpha               =1.0d0               , &amp;
     &amp;                       beta                =3.0d0               , &amp;
     &amp;                       gamma               =0.0d0               , &amp;
     &amp;                       darkMatterHaloScale_=darkMatterHaloScale_  &amp;
     &amp;                      )
   </constructor>
  </referenceConstruct>
  !!]
  massDistribution_       => darkMatterProfileDMOZhao1996_%get                   (node)
  kinematicsDistribution_ => massDistribution_            %kinematicsDistribution(    )
  radiusVirial            =  darkMatterHaloScale_         %radiusVirial          (node)
  timeFreefall            =  [2.912172466390932d0,3.227873348829691d0,3.870044363719088d0,5.196579302974257d0,8.00677623636939d0,1.417018715413697d1,2.823093711816523d1]
  select type (massDistribution_)
  class is (massDistributionSpherical)
     do i=1,7
        coordinates=[radiusScale*radius(i),0.0d0,0.0d0]
        mass                                      (i)=massDistribution_      %massEnclosedBySphere             (                                                         radiusScale   *radius      (i)                                        )
        density                                   (i)=massDistribution_      %density                          (                                                                        coordinates                                            )*radiusScale**3
        radiusRecoveredFromMass                   (i)=massDistribution_      %radiusEnclosingMass              (                                                 mass(i)                                                                       )/radiusScale
        radiusRecoveredFromDensity                (i)=massDistribution_      %radiusEnclosingDensity           (                                  3.0d0/4.0d0/Pi*mass(i)/radiusScale**3/radius      (i)**3                                     )/radiusScale
        radiusRecoveredFromSpecificAngularMomentum(i)=massDistribution_      %radiusFromSpecificAngularMomentum(             sqrt(gravitationalConstant_internal*mass(i)*radiusScale   *radius      (i)   )                                    )/radiusScale
        radiusRecoveredFromTimeFreefall           (i)=massDistribution_      %radiusFreefall                   (                                                         timeScale     *timeFreefall(i)                                        )/radiusScale
        potential                                 (i)=massDistribution_      %potential                        (                                                                        coordinates                                            )*radiusScale/gravitationalConstant_internal
        fourier                                   (i)=massDistribution_      %fourierTransform                 (radiusVirial,1.0d0                                      /radiusScale   /radius      (i)                                        )
        radialVelocityDispersion                  (i)=kinematicsDistribution_%velocityDispersion1D             (                                                                        coordinates        ,massDistribution_,massDistribution_)
     end do
     radiusVelocityMaximum   =massDistribution_%radiusRotationCurveMaximum  (                                       )
     velocityMaximum         =massDistribution_%velocityRotationCurveMaximum(                                       )
     velocityMaximumIndirect =massDistribution_%rotationCurve               (radiusVelocityMaximum                  )
     energyPotential         =massDistribution_%energyPotential             (radiusVirial                           )
     energyPotentialNumerical=massDistribution_%energyPotentialNumerical    (radiusVirial                           )
     energyKinetic           =massDistribution_%energyKinetic               (radiusVirial         ,massDistribution_)
     energyKineticNumerical  =massDistribution_%energyKineticNumerical      (radiusVirial         ,massDistribution_)
  end select
  !![
  <objectDestructor name="massDistribution_"            />
  <objectDestructor name="kinematicsDistribution_"      />
  <objectDestructor name="darkMatterProfileDMOZhao1996_"/>
  !!]
  ! Radial velocity dispersion in units of virial velocity.
  radialVelocityDispersion=radialVelocityDispersion/darkMatterHaloScale_%velocityVirial(node)
  call Assert(                         &
       &      'enclosed mass'        , &
       &      mass                   , &
       &      [                        &
       &       5.464789985591522d-4  , &
       &       3.442068263974007d-3  , &
       &       1.815032503316619d-2  , &
       &       7.461855208928223d-2  , &
       &       2.296390885460238d-1  , &
       &       5.359157644285496d-1  , &
       &       1.000000000000000d+0    &
       &      ]                      , &
       &      relTol=1.0d-6            &
       &     )
  call Assert(                         &
       &      'radius enclosing mass', &
       &      radiusRecoveredFromMass, &
       &      radius                 , &
       &      relTol=1.0d-6            &
       &     )
  call Assert(                         &
       &      'density'              , &
       &      density                , &
       &      [                        &
       &       6.119719178912785d-2  , &
       &       4.461275281427421d-2  , &
       &       2.581756528603831d-2  , &
       &       1.089178535504741d-2  , &
       &       3.227195660754789d-3  , &
       &       6.970742627230344d-4  , &
       &       1.195257652131403d-4    &
       &      ]                      , &
       &      relTol=1.0d-6            &
       &     )
  call Assert(                            &
       &      'radius enclosing density', &
       &      radiusRecoveredFromDensity, &
       &      radius                    , &
       &      relTol=2.0d-4               &
       &     )
  call Assert(                                            &
       &      'radius from specific angular momentum'   , &
       &      radiusRecoveredFromSpecificAngularMomentum, &
       &      radius                                    , &
       &      relTol=2.0d-4                               &
       &     )
  call Assert(                                            &
       &      'radius from freefall time'               , &
       &      radiusRecoveredFromTimeFreefall           , &
       &      radius                                    , &
       &      relTol=2.0d-4                               &
       &     )
  call Assert(                       &
       &      'potential'          , &
       &      +potential             &
       &      -potential(1)        , &
       &      [                      &
       &       2.387190797876185d-3, &
       &       8.130960771876030d-3, &
       &       2.453055501081222d-2, &
       &       6.225165933429317d-2, &
       &       1.285052760355665d-1, &
       &       2.164088001372156d-1, &
       &       3.075774583263617d-1  &
       &      ]                      &
       &      -2.387190797876185d-3, &
       &      relTol=1.0d-6          &
       &     )
  call Assert(                         &
       &      'fourier'              , &
       &      fourier                , &
       &      [                        &
       &       9.697536800190850d-4  , &
       &       8.183597728241470d-3  , &
       &       5.062326038988727d-2  , &
       &       1.681631870688479d-1  , &
       &       4.620457662739311d-1  , &
       &       8.150847326475010d-1  , &
       &       9.497559618451080d-1    &
       &      ]                      , &
       &      relTol=1.0d-6            &
       &     )
  call Assert(                         &
       &     'peak rotation velocity', &
       &     velocityMaximum         , &
       &     velocityMaximumIndirect , &
       &     relTol=1.0d-6             &
       &     )
  call Assert(                         &
       &     'potential energy'      , &
       &     energyPotential         , &
       &     energyPotentialNumerical, &
       &     relTol=1.0d-6             &
       &     )
  call Assert(                         &
       &     'kinetic energy'        , &
       &     energyKinetic           , &
       &     energyKineticNumerical  , &
       &     relTol=1.0d-6             &
       &     )
  call Assert(                              &
       &      'radial velocity dispersion', &
       &      radialVelocityDispersion    , &
       &      [                             &
       &       5.214290018440262d-1       , &
       &       5.654096128890455d-1       , &
       &       6.175226284000780d-1       , &
       &       6.597522036704852d-1       , &
       &       6.699313413583514d-1       , &
       &       6.379354590281848d-1       , &
       &       5.714635036404890d-1         &
       &      ]                           , &
       &      relTol=1.0d-6                 &
       &     )
  call Unit_Tests_End_Group()
  !! Special case: cored NFW: (α,β,γ) = (1,3,½).
  call Unit_Tests_Begin_Group('(α,β,γ) = (1,3,½)')
  allocate(darkMatterProfileDMOZhao1996_)
  !![
   <referenceConstruct object="darkMatterProfileDMOZhao1996_">
   <constructor>
    darkMatterProfileDMOZhao1996 (                                          &amp;
     &amp;                       alpha               =1.0d0               , &amp;
     &amp;                       beta                =3.0d0               , &amp;
     &amp;                       gamma               =0.5d0               , &amp;
     &amp;                       darkMatterHaloScale_=darkMatterHaloScale_  &amp;
     &amp;                      )
   </constructor>
  </referenceConstruct>
  !!]
  massDistribution_       => darkMatterProfileDMOZhao1996_%get                   (node)
  kinematicsDistribution_ => massDistribution_            %kinematicsDistribution(    )
  radiusVirial            =  darkMatterHaloScale_         %radiusVirial          (node)
  timeFreefall            =  [1.596420780360787d0,2.059037739059005d0,2.826539348071723d0,4.248077776277237d0,7.134952215026308d0,1.342236886351155d1,2.785866097530184d1]
  select type (massDistribution_)
  class is (massDistributionSpherical)
     do i=1,7
        coordinates=[radiusScale*radius(i),0.0d0,0.0d0]
        mass                                      (i)=massDistribution_      %massEnclosedBySphere             (                                                         radiusScale   *radius      (i)                                        )
        density                                   (i)=massDistribution_      %density                          (                                                                        coordinates                                            )*radiusScale**3
        radiusRecoveredFromMass                   (i)=massDistribution_      %radiusEnclosingMass              (                                                 mass(i)                                                                       )/radiusScale
        radiusRecoveredFromDensity                (i)=massDistribution_      %radiusEnclosingDensity           (                                  3.0d0/4.0d0/Pi*mass(i)/radiusScale**3/radius      (i)**3                                     )/radiusScale
        radiusRecoveredFromSpecificAngularMomentum(i)=massDistribution_      %radiusFromSpecificAngularMomentum(             sqrt(gravitationalConstant_internal*mass(i)*radiusScale   *radius      (i)   )                                    )/radiusScale
        radiusRecoveredFromTimeFreefall           (i)=massDistribution_      %radiusFreefall                   (                                                         timeScale     *timeFreefall(i)                                        )/radiusScale
        potential                                 (i)=massDistribution_      %potential                        (                                                                        coordinates                                            )*radiusScale/gravitationalConstant_internal
        fourier                                   (i)=massDistribution_      %fourierTransform                 (radiusVirial,1.0d0                                      /radiusScale   /radius      (i)                                        )
        radialVelocityDispersion                  (i)=kinematicsDistribution_%velocityDispersion1D             (                                                                        coordinates        ,massDistribution_,massDistribution_)
     end do
     radiusVelocityMaximum   =massDistribution_%radiusRotationCurveMaximum  (                                       )
     velocityMaximum         =massDistribution_%velocityRotationCurveMaximum(                                       )
     velocityMaximumIndirect =massDistribution_%rotationCurve               (radiusVelocityMaximum                  )
     energyPotential         =massDistribution_%energyPotential             (radiusVirial                           )
     energyPotentialNumerical=massDistribution_%energyPotentialNumerical    (radiusVirial                           )
     energyKinetic           =massDistribution_%energyKinetic               (radiusVirial         ,massDistribution_)
     energyKineticNumerical  =massDistribution_%energyKineticNumerical      (radiusVirial         ,massDistribution_)
  end select
  !![
  <objectDestructor name="massDistribution_"            />
  <objectDestructor name="kinematicsDistribution_"      />
  <objectDestructor name="darkMatterProfileDMOZhao1996_"/>
  !!]
  ! Radial velocity dispersion in units of virial velocity.
  radialVelocityDispersion=radialVelocityDispersion/darkMatterHaloScale_%velocityVirial(node)
  call Assert(                         &
       &      'enclosed mass'        , &
       &      mass                   , &
       &      [                        &
       &       1.654826011427427d-3  , &
       &       7.739711640396846d-3  , &
       &       3.140778408165459d-2  , &
       &       1.043599712384584d-1  , &
       &       2.742860732070792d-1  , &
       &       5.747348550741195d-1  , &
       &       1.000000000000000d+0    &
       &      ]                      , &
       &      relTol=1.0d-6            &
       &     )
  call Assert(                         &
       &      'radius enclosing mass', &
       &      radiusRecoveredFromMass, &
       &      radius                 , &
       &      relTol=1.0d-6            &
       &     )
  call Assert(                         &
       &      'density'              , &
       &      density                , &
       &      [                        &
       &1.550807828980169d-1,&
       &       8.42653949330162d-2,&
       &       3.777297120800491d-2,&
       &       1.301125858993861d-2,&
       &       3.338690510843061d-3,&
       &       6.583233979141888d-4,&
       &       1.070885480653384d-4&
       &      ]                      , &
       &      relTol=1.0d-6            &
       &     )
  call Assert(                            &
       &      'radius enclosing density', &
       &      radiusRecoveredFromDensity, &
       &      radius                    , &
       &      relTol=2.0d-4               &
       &     )
  call Assert(                                            &
       &      'radius from specific angular momentum'   , &
       &      radiusRecoveredFromSpecificAngularMomentum, &
       &      radius                                    , &
       &      relTol=2.0d-4                               &
       &     )
  call Assert(                                            &
       &      'radius from freefall time'               , &
       &      radiusRecoveredFromTimeFreefall           , &
       &      radius                                    , &
       &      relTol=2.0d-4                               &
       &     )
  call Assert(                       &
       &      'potential'          , &
       &      +potential             &
       &      -potential(1)        , &
       &      [                      &
       &9.59892229601949d-3,&
       &       2.419272545370556d-2,&
       &       5.585172068801785d-2,&
       &       1.136457588954521d-1,&
       &       1.98498742003645d-1,&
       &       2.975288623538137d-1,&
       &       3.917543232803806d-1&
       &      ]                      &
       &      -9.59892229601949d-3, &
       &      relTol=1.0d-6          &
       &     )
  call Assert(                         &
       &      'fourier'              , &
       &      fourier                , &
       &      [                        &
       &3.707401665292136d-3,&
       &       1.826316873546129d-2,&
       &       7.540534190396314d-2,&
       &       2.092215103402299d-1,&
       &       4.994927051818294d-1,&
       &       8.28909754178526d-1,&
       &       9.53572760138475d-1&
       &      ]                      , &
       &      relTol=1.0d-6            &
       &     )
  call Assert(                         &
       &     'peak rotation velocity', &
       &     velocityMaximum         , &
       &     velocityMaximumIndirect , &
       &     relTol=1.0d-6             &
       &     )
  call Assert(                         &
       &     'potential energy'      , &
       &     energyPotential         , &
       &     energyPotentialNumerical, &
       &     relTol=1.0d-6             &
       &     )
  call Assert(                         &
       &     'kinetic energy'        , &
       &     energyKinetic           , &
       &     energyKineticNumerical  , &
       &     relTol=1.0d-6             &
       &     )
  call Assert(                              &
       &      'radial velocity dispersion', &
       &      radialVelocityDispersion    , &
       &      [                             &
       &       5.347095696836751d-1       , &
       &       6.083588653620841d-1       , &
       &       6.691204678148943d-1       , &
       &       7.007128495278605d-1       , &
       &       6.909164992407667d-1       , &
       &       6.403334146398749d-1       , &
       &       5.617664344971543d-1         &
       &      ]                           , &
       &      relTol=1.0d-6                 &
       &     )
  call Unit_Tests_End_Group()
  !! Special case: cored NFW: (α,β,γ) = (1,3,1½).
  call Unit_Tests_Begin_Group('(α,β,γ) = (1,3,1½)')
  allocate(darkMatterProfileDMOZhao1996_)
  !![
   <referenceConstruct object="darkMatterProfileDMOZhao1996_">
   <constructor>
    darkMatterProfileDMOZhao1996 (                                          &amp;
     &amp;                       alpha               =1.0d0               , &amp;
     &amp;                       beta                =3.0d0               , &amp;
     &amp;                       gamma               =1.5d0               , &amp;
     &amp;                       darkMatterHaloScale_=darkMatterHaloScale_  &amp;
     &amp;                      )
   </constructor>
  </referenceConstruct>
  !!]
  massDistribution_       => darkMatterProfileDMOZhao1996_%get                   (node)
  kinematicsDistribution_ => massDistribution_            %kinematicsDistribution(    )
  radiusVirial            =  darkMatterHaloScale_         %radiusVirial          (node)
  timeFreefall            =  [4.596100900700435d-1,8.05108081560643d-1,1.455969425989372d0,2.760391236512836d0,5.555351793104961d0,1.191020177556312d1,2.703975097809943d1]
  select type (massDistribution_)
  class is (massDistributionSpherical)
     do i=1,7
        coordinates=[radiusScale*radius(i),0.0d0,0.0d0]
        mass                                      (i)=massDistribution_      %massEnclosedBySphere             (                                                         radiusScale   *radius      (i)                                        )
        density                                   (i)=massDistribution_      %density                          (                                                                        coordinates                                            )*radiusScale**3
        radiusRecoveredFromMass                   (i)=massDistribution_      %radiusEnclosingMass              (                                                 mass(i)                                                                       )/radiusScale
        radiusRecoveredFromDensity                (i)=massDistribution_      %radiusEnclosingDensity           (                                  3.0d0/4.0d0/Pi*mass(i)/radiusScale**3/radius      (i)**3                                     )/radiusScale
        radiusRecoveredFromSpecificAngularMomentum(i)=massDistribution_      %radiusFromSpecificAngularMomentum(             sqrt(gravitationalConstant_internal*mass(i)*radiusScale   *radius      (i)   )                                    )/radiusScale
        radiusRecoveredFromTimeFreefall           (i)=massDistribution_      %radiusFreefall                   (                                                         timeScale     *timeFreefall(i)                                        )/radiusScale
        potential                                 (i)=massDistribution_      %potential                        (                                                                        coordinates                                            )*radiusScale/gravitationalConstant_internal
        fourier                                   (i)=massDistribution_      %fourierTransform                 (radiusVirial,1.0d0                                      /radiusScale   /radius      (i)                                        )
        radialVelocityDispersion                  (i)=kinematicsDistribution_%velocityDispersion1D             (                                                                        coordinates        ,massDistribution_,massDistribution_)
     end do
     radiusVelocityMaximum   =massDistribution_%radiusRotationCurveMaximum  (                                       )
     velocityMaximum         =massDistribution_%velocityRotationCurveMaximum(                                       )
     velocityMaximumIndirect =massDistribution_%rotationCurve               (radiusVelocityMaximum                  )
     energyPotential         =massDistribution_%energyPotential             (radiusVirial                           )
     energyPotentialNumerical=massDistribution_%energyPotentialNumerical    (radiusVirial                           )
     energyKinetic           =massDistribution_%energyKinetic               (radiusVirial         ,massDistribution_)
     energyKineticNumerical  =massDistribution_%energyKineticNumerical      (radiusVirial         ,massDistribution_)
  end select
  !![
  <objectDestructor name="massDistribution_"            />
  <objectDestructor name="kinematicsDistribution_"      />
  <objectDestructor name="darkMatterProfileDMOZhao1996_"/>
  !!]
  ! Radial velocity dispersion in units of virial velocity.
  radialVelocityDispersion=radialVelocityDispersion/darkMatterHaloScale_%velocityVirial(node)
  call Assert(                         &
       &      'enclosed mass'        , &
       &      mass                   , &
       &      [                        &
       &       1.614787314131082d-2  , &
       &       4.146438397463787d-2  , &
       &       9.894487896260260d-2  , &
       &       2.125365304218604d-1  , &
       &       4.021269908070108d-1  , &
       &       6.698167367007793d-1  , &
       &       1.000000000000000d+0    &
       &      ]                      , &
       &      relTol=1.0d-6            &
       &     )
  call Assert(                         &
       &      'radius enclosing mass', &
       &      radiusRecoveredFromMass, &
       &      radius                 , &
       &      relTol=1.0d-6            &
       &     )
  call Assert(                         &
       &      'density'              , &
       &      density                , &
       &      [                        &
       &       9.202064069700650d-1  , &
       &       2.777819507076467d-1  , &
       &       7.471144923386984d-2  , &
       &       1.715671205314034d-2  , &
       &       3.301810774096491d-3  , &
       &       5.425428724758725d-4  , &
       &       7.942922246824230d-5    &
       &      ]                      , &
       &      relTol=1.0d-6            &
       &     )
  call Assert(                            &
       &      'radius enclosing density', &
       &      radiusRecoveredFromDensity, &
       &      radius                    , &
       &      relTol=2.0d-4               &
       &     )
  call Assert(                                            &
       &      'radius from specific angular momentum'   , &
       &      radiusRecoveredFromSpecificAngularMomentum, &
       &      radius                                    , &
       &      relTol=2.0d-4                               &
       &     )
  call Assert(                                            &
       &      'radius from freefall time'               , &
       &      radiusRecoveredFromTimeFreefall           , &
       &      radius                                    , &
       &      relTol=2.0d-4                               &
       &     )
  call Assert(                       &
       &      'potential'          , &
       &      +potential             &
       &      -potential(1)        , &
       &      [                      &
       &       2.773517522337437d-1, &
       &       3.795660488783141d-1, &
       &       5.062490622313113d-1, &
       &       6.498538783125577d-1, &
       &       7.947391738552302d-1, &
       &       9.233929853785370d-1, &
       &       1.024853878312558d+0  &
       &      ]                      &
       &      -2.773517522337437d-1, &
       &      relTol=1.0d-6          &
       &     )
  call Assert(                         &
       &      'fourier'              , &
       &      fourier                , &
       &      [                        &
       &       3.012898538414511d-2  , &
       &       7.490502708109203d-2  , &
       &       1.710991422211459d-1  , &
       &       3.320442821775421d-1  , &
       &       5.969264990169439d-1  , &
       &       8.639011197411600d-1  , &
       &       9.631744296899940d-1    &
       &      ]                      , &
       &      relTol=1.0d-5            &
       &     )
  call Assert(                         &
       &     'peak rotation velocity', &
       &     velocityMaximum         , &
       &     velocityMaximumIndirect , &
       &     relTol=1.0d-6             &
       &     )
  call Assert(                         &
       &     'potential energy'      , &
       &     energyPotential         , &
       &     energyPotentialNumerical, &
       &     relTol=1.0d-6             &
       &     )
  call Assert(                         &
       &     'kinetic energy'        , &
       &     energyKinetic           , &
       &     energyKineticNumerical  , &
       &     relTol=1.0d-6             &
       &     )
  call Assert(                              &
       &      'radial velocity dispersion', &
       &      radialVelocityDispersion    , &
       &      [                             &
       &       8.363292318000940d-1       , &
       &       8.746095866642690d-1       , &
       &       8.733596657690690d-1       , &
       &       8.296377525840130d-1       , &
       &       7.505291252896317d-1       , &
       &       6.497937689634748d-1       , &
       &       5.421632188128948d-1         &
       &      ]                           , &
       &      relTol=1.0d-6                 &
       &     )
  call Unit_Tests_End_Group()    
  !! General case: (α,β,γ) = (1,2,1).
  call Unit_Tests_Begin_Group('(α,β,γ) = (1,2,1)')
  allocate(darkMatterProfileDMOZhao1996_)
  !![
   <referenceConstruct object="darkMatterProfileDMOZhao1996_">
   <constructor>
    darkMatterProfileDMOZhao1996 (                                          &amp;
     &amp;                       alpha               =1.0d0               , &amp;
     &amp;                       beta                =2.0d0               , &amp;
     &amp;                       gamma               =1.0d0               , &amp;
     &amp;                       darkMatterHaloScale_=darkMatterHaloScale_  &amp;
     &amp;                      )
   </constructor>
  </referenceConstruct>
  !!]
  massDistribution_       => darkMatterProfileDMOZhao1996_%get                   (node)
  kinematicsDistribution_ => massDistribution_            %kinematicsDistribution(    )
  radiusVirial            =  darkMatterHaloScale_         %radiusVirial          (node)
  timeFreefall            =  [1.760884513374176d0,2.56756519338419d0,3.83552259760776d0,5.9420179316033d0,9.64631439102321d0,1.646977873143257d1,2.944555422801916d1]
  select type (massDistribution_)
  class is (massDistributionSpherical)
     do i=1,7
        coordinates=[radiusScale*radius(i),0.0d0,0.0d0]
        mass                                      (i)=massDistribution_      %massEnclosedBySphere             (                                                         radiusScale   *radius      (i)                                        )
        density                                   (i)=massDistribution_      %density                          (                                                                        coordinates                                            )*radiusScale**3
        radiusRecoveredFromMass                   (i)=massDistribution_      %radiusEnclosingMass              (                                                 mass(i)                                                                       )/radiusScale
        radiusRecoveredFromDensity                (i)=massDistribution_      %radiusEnclosingDensity           (                                  3.0d0/4.0d0/Pi*mass(i)/radiusScale**3/radius      (i)**3                                     )/radiusScale
        radiusRecoveredFromSpecificAngularMomentum(i)=massDistribution_      %radiusFromSpecificAngularMomentum(             sqrt(gravitationalConstant_internal*mass(i)*radiusScale   *radius      (i)   )                                    )/radiusScale
        radiusRecoveredFromTimeFreefall           (i)=massDistribution_      %radiusFreefall                   (                                                         timeScale     *timeFreefall(i)                                        )/radiusScale
        potential                                 (i)=massDistribution_      %potential                        (                                                                        coordinates                                            )*radiusScale/gravitationalConstant_internal
        fourier                                   (i)=massDistribution_      %fourierTransform                 (radiusVirial,1.0d0                                      /radiusScale   /radius      (i)                                        )
        radialVelocityDispersion                  (i)=kinematicsDistribution_%velocityDispersion1D             (                                                                        coordinates        ,massDistribution_,massDistribution_)
     end do
     radiusVelocityMaximum   =massDistribution_%radiusRotationCurveMaximum  (                                       )
     velocityMaximum         =massDistribution_%velocityRotationCurveMaximum(                                       )
     velocityMaximumIndirect =massDistribution_%rotationCurve               (radiusVelocityMaximum                  )
     energyPotential         =massDistribution_%energyPotential             (radiusVirial                           )
     energyPotentialNumerical=massDistribution_%energyPotentialNumerical    (radiusVirial                           )
     energyKinetic           =massDistribution_%energyKinetic               (radiusVirial         ,massDistribution_)
     energyKineticNumerical  =massDistribution_%energyKineticNumerical      (radiusVirial         ,massDistribution_)
   end select
  !![
  <objectDestructor name="massDistribution_"            />
  <objectDestructor name="kinematicsDistribution_"      />
  <objectDestructor name="darkMatterProfileDMOZhao1996_"/>
  !!]
  ! Radial velocity dispersion in units of virial velocity.
  radialVelocityDispersion=radialVelocityDispersion/darkMatterHaloScale_%velocityVirial(node)
  call Assert(                         &
       &      'enclosed mass'        , &
       &      mass                   , &
       &      [                        &
       &       1.243709056088815d-3  , &
       &       4.628207492038648d-3  , &
       &       1.629132354883366d-2  , &
       &       5.288035415632079d-2  , &
       &       1.553373421641236d-1  , &
       &       4.119687414110722d-1  , &
       &       1.000000000000000d+0    &
       &      ]                      , &
       &      relTol=1.0d-6            &
       &     )
  call Assert(                         &
       &      'radius enclosing mass', &
       &      radiusRecoveredFromMass, &
       &      radius                 , &
       &      relTol=1.0d-6            &
       &     )
  call Assert(                         &
       &      'density'              , &
       &      density                , &
       &      [                        &
       &       9.751958345559170d-2  , &
       &       4.388381255501626d-2  , &
       &       1.828492189792344d-2  , &
       &       6.856845711721292d-3  , &
       &       2.285615237240430d-3  , &
       &       6.856845711721292d-4  , &
       &       1.904679364367025d-4    &
       &      ]                      , &
       &      relTol=1.0d-6            &
       &     )
  call Assert(                            &
       &      'radius enclosing density', &
       &      radiusRecoveredFromDensity, &
       &      radius                    , &
       &      relTol=2.0d-4               &
       &     )
  call Assert(                                            &
       &      'radius from specific angular momentum'   , &
       &      radiusRecoveredFromSpecificAngularMomentum, &
       &      radius                                    , &
       &      relTol=2.0d-4                               &
       &     )
  call Assert(                                            &
       &      'radius from freefall time'               , &
       &      radiusRecoveredFromTimeFreefall           , &
       &      radius                                    , &
       &      relTol=2.0d-4                               &
       &     )
  call Assert(                       &
       &      'potential'          , &
       &      +potential             &
       &      -potential(1)        , &
       &      [                      &
       &       1.034803460994295d-2, &
       &       1.994179476929133d-2, &
       &       3.729169381246816d-2, &
       &       6.657062060529668d-2, &
       &       1.116566445896912d-1, &
       &       1.743643889079128d-1, &
       &       2.536506313435059d-1  &
       &      ]                      &
       &      -1.034803460994295d-2, &
       &      relTol=1.0d-6          &
       &     )
  call Assert(                         &
       &      'fourier'              , &
       &      fourier                , &
       &      [                        &
       &       2.498947569002664d-3  , &
       &       8.859191321495360d-3  , &
       &       3.901409916856754d-2  , &
       &       1.078443248271426d-1  , &
       &       3.537730117747184d-1  , &
       &       7.703987792232572d-1  , &
       &       9.371207041223680d-1    &
       &      ]                      , &
       &      relTol=1.0d-5            &
       &     )
  call Assert(                         &
       &     'peak rotation velocity', &
       &     velocityMaximum         , &
       &     velocityMaximumIndirect , &
       &     relTol=1.0d-6             &
       &     )
  call Assert(                         &
       &     'potential energy'      , &
       &     energyPotential         , &
       &     energyPotentialNumerical, &
       &     relTol=1.0d-6             &
       &     )
  call Assert(                         &
       &     'kinetic energy'        , &
       &     energyKinetic           , &
       &     energyKineticNumerical  , &
       &     relTol=1.0d-6             &
       &     )
  call Assert(                              &
       &      'radial velocity dispersion', &
       &      radialVelocityDispersion    , &
       &      [                             &
       &       3.790874753794872d-1       , &
       &       4.539034405088716d-1       , &
       &       5.289441877216294d-1       , &
       &       5.997499368748768d-1       , &
       &       6.624381831647858d-1       , &
       &       7.141526014953340d-1       , &
       &       7.536975187158942d-1         &
       &      ]                           , &
       &      relTol=1.0d-3                 &
       &     )
  call Unit_Tests_End_Group()  
  call Unit_Tests_End_Group()
  ! Test Burkert profile.
  call Unit_Tests_Begin_Group('Burkert profile')
  massDistribution_       => darkMatterProfileDMOBurkert_%get                   (node)
  kinematicsDistribution_ => massDistribution_           %kinematicsDistribution(    )
  radiusVirial            =  darkMatterHaloScale_        %radiusVirial          (node)
  timeFreefall            =  [3.378196272398096d0,3.533466912782305d0,3.901024113095533d0,4.840727096354573d0,7.260451446070998d0,1.323150559822296d1,2.774028478732573d1]
  select type (massDistribution_)
  class is (massDistributionSpherical)
     do i=1,7
        coordinates=[radiusScale*radius(i),0.0d0,0.0d0]
        mass                                      (i)=massDistribution_      %massEnclosedBySphere             (                                                         radiusScale   *radius      (i)                                        )
        density                                   (i)=massDistribution_      %density                          (                                                                        coordinates                                            )*radiusScale**3
        radiusRecoveredFromMass                   (i)=massDistribution_      %radiusEnclosingMass              (                                                 mass(i)                                                                       )/radiusScale
        radiusRecoveredFromDensity                (i)=massDistribution_      %radiusEnclosingDensity           (                                  3.0d0/4.0d0/Pi*mass(i)/radiusScale**3/radius      (i)**3                                     )/radiusScale
        radiusRecoveredFromSpecificAngularMomentum(i)=massDistribution_      %radiusFromSpecificAngularMomentum(             sqrt(gravitationalConstant_internal*mass(i)*radiusScale   *radius      (i)   )                                    )/radiusScale
        radiusRecoveredFromTimeFreefall           (i)=massDistribution_      %radiusFreefall                   (                                                         timeScale     *timeFreefall(i)                                        )/radiusScale
        potential                                 (i)=massDistribution_      %potential                        (                                                                        coordinates                                            )*radiusScale/gravitationalConstant_internal
        fourier                                   (i)=massDistribution_      %fourierTransform                 (radiusVirial,1.0d0                                      /radiusScale   /radius      (i)                                        )
        radialVelocityDispersion                  (i)=kinematicsDistribution_%velocityDispersion1D             (                                                                        coordinates        ,massDistribution_,massDistribution_)
     end do
     radiusVelocityMaximum   =massDistribution_%radiusRotationCurveMaximum  (                     )
     velocityMaximum         =massDistribution_%velocityRotationCurveMaximum(                     )
     velocityMaximumIndirect =massDistribution_%rotationCurve               (radiusVelocityMaximum)
     energyPotential         =massDistribution_%energyPotential             (radiusVirial         )
     energyPotentialNumerical=massDistribution_%energyPotentialNumerical    (radiusVirial         )
  end select
  !![
  <objectDestructor name="massDistribution_"      />
  <objectDestructor name="kinematicsDistribution_"/>
  !!]
  ! Radial velocity dispersion in units of virial velocity.
  radialVelocityDispersion=radialVelocityDispersion/darkMatterHaloScale_%velocityVirial(node)
  call Assert(                        &
       &      'enclosed mass'       , &
       &      mass                  , &
       &      [                       &
       &       4.1583650166653620d-4, &
       &       2.9870571374971090d-3, &
       &       1.8812441757378840d-2, &
       &       8.9614051908432800d-2, &
       &       0.2805458115927276d+0, &
       &       0.5990982283029260d+0, &
       &       1.0000000000000000d+0  &
       &      ]                     , &
       &      relTol=1.0d-6           &
       &     )
  call Assert(                         &
       &      'radius enclosing mass', &
       &      radiusRecoveredFromMass, &
       &      radius                 , &
       &      relTol=1.0d-6            &
       &     )
  call Assert(                        &
       &      'density'             , &
       &      density               , &
       &      [                       &
       &       4.9082352873191440d-2, &
       &       4.2225259457083800d-2, &
       &       2.9909558782101030d-2, &
       &       1.4020105679109860d-2, &
       &       3.7386948477626290d-3, &
       &       6.5976967901693440d-4, &
       &       9.5863970455452000d-5  &
       &      ]                     , &
       &      relTol=1.0d-6           &
       &     )
  call Assert(                            &
       &      'radius enclosing density', &
       &      radiusRecoveredFromDensity, &
       &      radius                    , &
       &      relTol=2.0d-4               &
       &     )
  call Assert(                                            &
       &      'radius from specific angular momentum'   , &
       &      radiusRecoveredFromSpecificAngularMomentum, &
       &      radius                                    , &
       &      relTol=2.0d-4                               &
       &     )
  call Assert(                                            &
       &      'radius from freefall time'               , &
       &      radiusRecoveredFromTimeFreefall           , &
       &      radius                                    , &
       &      relTol=2.0d-4                               &
       &     )
  call Assert(                         &
       &     'peak rotation velocity', &
       &     velocityMaximum         , &
       &     velocityMaximumIndirect , &
       &     relTol=1.0d-6             &
       &     )
  call Assert(                        &
       &      'potential'          ,  &
       &      +potential              &
       &      -potential(1)        ,  &
       &      [                       &
       &       -5.517710030243666d-1, &
       &       -5.470649572158022d-1, &
       &       -5.313012274275472d-1, &
       &       -4.884797937826585d-1, &
       &       -4.072028257618435d-1, &
       &       -3.040428693416573d-1, &
       &       -2.075890931613922d-1  &
       &      ]                       &
       &       +5.517710030243666d-1, &
       &      relTol=1.0d-6           &
       &     )
  call Assert(                        &
       &      'fourier'             , &
       &      fourier               , &
       &      [                       &
       &       3.2941046717529910d-4, &
       &       6.7507368425877680d-3, &
       &       6.0387631952687390d-2, &
       &       0.2118984282100852d+0, &
       &       0.5171391325515731d+0, &
       &       0.8364966656849830d+0, &
       &       0.9557322027757890d+0  &
       &      ]                     , &
       &      relTol=1.0d-6           &
       &     )
  call Assert(                         &
       &     'potential energy'      , &
       &     energyPotential         , &
       &     energyPotentialNumerical, &
       &     relTol=1.0d-6             &
       &     )
  call Assert(                                         &
       &      'radial velocity dispersion (analytic)', &
       &      radialVelocityDispersion               , &
       &      [                                        &
       &       6.395640553962476d-1                  , &
       &       6.595971570604928d-1                  , &
       &       6.818328396903716d-1                  , &
       &       6.957079957624708d-1                  , &
       &       6.836738566635967d-1                  , &
       &       6.307711855060676d-1                  , &
       &       5.466303194157120d-1                    &
       &      ]                                      , &
       &      relTol=1.0d-6                            &
       &     )
  call Unit_Tests_End_Group       ()
  ! Test NFW profile.
  call Unit_Tests_Begin_Group('NFW profile')
  massDistribution_       => darkMatterProfileDMONFW_%get                   (node)
  kinematicsDistribution_ => massDistribution_       %kinematicsDistribution(    )
  radiusVirial            =  darkMatterHaloScale_    %radiusVirial          (node)
  timeFreefall            =  [0.864113d0,1.29807d0,2.04324d0,3.44421d0,6.32148d0,1.26727d1,2.74642d1]
  select type (massDistribution_)
  class is (massDistributionSpherical)
     do i=1,7
        coordinates=[radiusScale*radius(i),0.0d0,0.0d0]
        mass                                      (i)=massDistribution_      %massEnclosedBySphere             (                                                         radiusScale   *radius      (i)                                        )
        density                                   (i)=massDistribution_      %density                          (                                                                        coordinates                                            )*radiusScale**3
        radiusRecoveredFromMass                   (i)=massDistribution_      %radiusEnclosingMass              (                                                 mass(i)                                                                       )/radiusScale
        radiusRecoveredFromDensity                (i)=massDistribution_      %radiusEnclosingDensity           (                                  3.0d0/4.0d0/Pi*mass(i)/radiusScale**3/radius      (i)**3                                     )/radiusScale
        radiusRecoveredFromSpecificAngularMomentum(i)=massDistribution_      %radiusFromSpecificAngularMomentum(             sqrt(gravitationalConstant_internal*mass(i)*radiusScale   *radius      (i)   )                                    )/radiusScale
        radiusRecoveredFromTimeFreefall           (i)=massDistribution_      %radiusFreefall                   (                                                         timeScale     *timeFreefall(i)                                        )/radiusScale
        potential                                 (i)=massDistribution_      %potential                        (                                                                        coordinates                                            )*radiusScale/gravitationalConstant_internal
        fourier                                   (i)=massDistribution_      %fourierTransform                 (radiusVirial,1.0d0                                      /radiusScale   /radius      (i)                                        )
        radialVelocityDispersion                  (i)=kinematicsDistribution_%velocityDispersion1D             (                                                                        coordinates        ,massDistribution_,massDistribution_)
     end do
     energyPotential         =massDistribution_%energyPotential             (radiusVirial                           )
     energyPotentialNumerical=massDistribution_%energyPotentialNumerical    (radiusVirial                           )
     energyKinetic           =massDistribution_%energyKinetic               (radiusVirial         ,massDistribution_)
     energyKineticNumerical  =massDistribution_%energyKineticNumerical      (radiusVirial         ,massDistribution_)
     radiusVelocityMaximum   =massDistribution_%radiusRotationCurveMaximum  (                                       )
     velocityMaximum         =massDistribution_%velocityRotationCurveMaximum(                                       )
     velocityMaximumIndirect =massDistribution_%rotationCurve               (radiusVelocityMaximum                  )
  end select
  !![
  <objectDestructor name="massDistribution_"      />
  <objectDestructor name="kinematicsDistribution_"/>
  !!]
  ! Radial velocity dispersion in units of virial velocity.
  radialVelocityDispersion=radialVelocityDispersion/darkMatterHaloScale_%velocityVirial(node)
  call Assert(                         &
       &      'enclosed mass'        , &
       &      mass                   , &
       &      [                        &
       &       5.099550982355504d-3  , &
       &       1.768930674181593d-2  , &
       &       5.513246746363203d-2  , &
       &       1.476281525188409d-1  , &
       &       3.301489257042704d-1  , &
       &       6.186775455118112d-1  , &
       &       1.000000000000000d+0    &
       &      ]                      , &
       &      relTol=1.0d-6            &
       &     )
  call Assert(                         &
       &      'radius enclosing mass', &
       &      radiusRecoveredFromMass, &
       &      radius                 , &
       &      relTol=1.0d-6            &
       &     )
  call Assert(                         &
       &      'density'              , &
       &      density                , &
       &      [                        &
       &       3.844641857939078d-1  , &
       &       1.557079952465327d-1  , &
       &       5.406527612726829d-2  , &
       &       1.520585891079421d-2  , &
       &       3.379079757954268d-3  , &
       &       6.082343564317684d-4  , &
       &       9.386332660984080d-5    &
       &      ]                      , &
       &      relTol=1.0d-6            &
       &     )
  call Assert(                            &
       &      'radius enclosing density', &
       &      radiusRecoveredFromDensity, &
       &      radius                    , &
       &      relTol=2.0d-4               &
       &     )
  call Assert(                                            &
       &      'radius from specific angular momentum'   , &
       &      radiusRecoveredFromSpecificAngularMomentum, &
       &      radius                                    , &
       &      relTol=2.0d-4                               &
       &     )
  call Assert(                                            &
       &      'radius from freefall time'               , &
       &      radiusRecoveredFromTimeFreefall           , &
       &      radius                                    , &
       &      relTol=2.0d-4                               &
       &     )
  call Assert(                       &
       &      'potential'          , &
       &      +potential             &
       &      -potential(1)        , &
       &      [                      &
       &       4.412912928902085d-2, &
       &       8.210873989889300d-2, &
       &       1.445116765163305d-1, &
       &       2.345367646465509d-1, &
       &       3.444787600350539d-1, &
       &       4.567944810866741d-1, &
       &       5.544042971829187d-1  &
       &      ]                      &
       &      -4.412912928902085d-2, &
       &      relTol=1.0d-6          &
       &     )
  call Assert(                         &
       &      'fourier'              , &
       &      fourier                , &
       &      [                        &
       &       1.099052711906094d-2  , &
       &       3.746284433831412d-2  , &
       &       1.127728075593513d-1  , &
       &       2.619030036621630d-1  , &
       &       5.434560723827092d-1  , &
       &       8.448608065763800d-1  , &
       &       9.579597044271230d-1    &
       &      ]                      , &
       &      relTol=1.0d-6            &
       &     )
  call Assert(                         &
       &     'peak rotation velocity', &
       &     velocityMaximum         , &
       &     velocityMaximumIndirect , &
       &     relTol=1.0d-6             &
       &     )
  call Assert(                         &
       &     'potential energy'      , &
       &     energyPotential         , &
       &     energyPotentialNumerical, &
       &     relTol=1.0d-6             &
       &     )
  call Assert(                         &
       &     'kinetic energy'        , &
       &     energyKinetic           , &
       &     energyKineticNumerical  , &
       &     relTol=1.0d-6             &
       &     )
  call Assert(                                                     &
       &      'radial velocity dispersion (analytic            )', &
       &      radialVelocityDispersion                           , &
       &      [                                                    &
       &       6.285734096346791d-1                              , &
       &       7.048421272327700d-1                              , &
       &       7.510776824829392d-1                              , &
       &       7.558757679348769d-1                              , &
       &       7.172187600402157d-1                              , &
       &       6.441447075213688d-1                              , &
       &       5.520559866965048d-1                                &
       &      ]                                                  , &
       &      relTol=1.0d-6                                        &
       &     )
  massDistribution_       => darkMatterProfileDMONFWSeriesExpansion_%get                   (node)
  kinematicsDistribution_ => massDistribution_                      %kinematicsDistribution(    )
  select type (massDistribution_)
  class is (massDistributionSpherical)
     do i=1,7
        coordinates                               =[radiusScale*radius(i),0.0d0,0.0d0]
        radialVelocityDispersionSeriesExpansion(i)=kinematicsDistribution_%velocityDispersion1D(coordinates,massDistribution_,massDistribution_)
     end do
  end select
  !![
  <objectDestructor name="massDistribution_"      />
  <objectDestructor name="kinematicsDistribution_"/>
  !!]
  ! Radial velocity dispersion in units of virial velocity.
  radialVelocityDispersionSeriesExpansion=radialVelocityDispersionSeriesExpansion/darkMatterHaloScale_%velocityVirial(node)
  call Assert(                                                     &
       &      'radial velocity dispersion (series approximation)', &
       &      radialVelocityDispersionSeriesExpansion            , &
       &      [                                                    &
       &       6.285734096346791d-1                              , &
       &       7.048421272327700d-1                              , &
       &       7.510776824829392d-1                              , &
       &       7.558757679348769d-1                              , &
       &       7.172187600402157d-1                              , &
       &       6.441447075213688d-1                              , &
       &       5.520559866965048d-1                                &
       &      ]                                                  , &
       &      relTol=1.0d-5                                        &
       &     )
  call Unit_Tests_End_Group       ()
  ! Test finite resolution NFW profile.
  call Unit_Tests_Begin_Group('Finite resolution NFW profile')
  massDistribution_ => darkMatterProfileDMOFiniteResolution_%get         (node)
  radiusVirial      =  darkMatterHaloScale_                 %radiusVirial(node)
  do i=1,7
     coordinates=[radiusScale*radius(i),0.0d0,0.0d0]
     mass   (i)=massDistribution_%massEnclosedBySphere(radius     =      radiusScale*radius(i)                         )
     density(i)=massDistribution_%density             (coordinates=      coordinates                                   )*radiusScale**3
     fourier(i)=massDistribution_%fourierTransform    (wavenumber =1.0d0/radiusScale/radius(i),radiusOuter=radiusVirial)
  end do
  !![
  <objectDestructor name="massDistribution_"/>
  !!]
  call Assert(                        &
       &      'enclosed mass'       , &
       &      mass                  , &
       &      [                       &
       &       8.18351042282589d-4,5.312260121162361d-3,2.774735862318237d-2,1.039240379463412d-1,2.7600761583923d-1,5.599835640058711d-1,9.39743487427036d-1   &
       &      ]                     , &
       &      relTol=1.0d-2           &
       &     )
  call Assert(                        &
       &      'density'             , &
       &      density               , &
       &      [                       &
       &       9.32462616056103d-2,6.963473240229225d-2,3.822992337631457d-2,1.36005336723227d-2,3.278188884572237d-3,6.035374957444192d-4,9.36805351479331d-5  &
       &      ]                     , &
       &      relTol=1.0d-6           &
       &     )
  call Assert(                        &
       &      'fourier'             , &
       &      fourier               , &
       &      [                       &
       &       1.677154217681274d-3,1.573965679970162d-2,7.938043609021705d-2,2.256540529519775d-1,5.18327690023956d-1,8.36115193021317d-1,9.55577051598601d-1 &
       &      ]                     , &
       &      relTol=1.0d-2           &
       &     )
  call Unit_Tests_End_Group       ()
  ! Construct profiles for SIDM models.
  !![
  <referenceConstruct object="cosmologyParametersPippin_"                        >
   <constructor>
    cosmologyParametersSimple                         (                                                                                       &amp;
     &amp;                                             OmegaMatter                         = 0.2660d0                                       , &amp;
     &amp;                                             OmegaBaryon                         = 0.0000d0                                       , &amp;
     &amp;                                             OmegaDarkEnergy                     = 0.7340d0                                       , &amp;
     &amp;                                             temperatureCMB                      = 2.7800d0                                       , &amp;
     &amp;                                             HubbleConstant                      =71.0000d0                                         &amp;
     &amp;                                            )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctionsPippin_"                         >
   <constructor>
    cosmologyFunctionsMatterLambda                    (                                                                                       &amp;
     &amp;                                             cosmologyParameters_                =cosmologyParametersPippin_                        &amp;
     &amp;                                            )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="virialDensityContrastPippin_"                      >
   <constructor>
    virialDensityContrastFixed                        (                                                                                       &amp;
     &amp;                                             densityContrastValue                =200.0d0                                         , &amp;
     &amp;                                             densityType                         =fixedDensityTypeCritical                        , &amp;
     &amp;                                             turnAroundOverVirialRadius          =  2.0d0                                         , &amp;
     &amp;                                             cosmologyParameters_                =cosmologyParametersPippin_                      , &amp;
     &amp;                                             cosmologyFunctions_                 =cosmologyFunctionsPippin_                         &amp;
     &amp;                                            )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterHaloScalePippin_"                        >
   <constructor>
    darkMatterHaloScaleVirialDensityContrastDefinition(                                                                                       &amp;
     &amp;                                             cosmologyParameters_                =cosmologyParametersPippin_                      , &amp;
     &amp;                                             cosmologyFunctions_                 =cosmologyFunctionsPippin_                       , &amp;
     &amp;                                             virialDensityContrast_              =virialDensityContrastPippin_                      &amp;
     &amp;                                            )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterParticleCDM_"                            >
   <constructor>
    darkMatterParticleCDM                             (                                                                                       &amp;
     &amp;                                            )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterParticleSelfInteractingDarkMatter_"      >
   <constructor>
    darkMatterParticleSelfInteractingDarkMatter       (                                                                                       &amp;
     &amp;                                             crossSectionSelfInteraction         =1.0d0                                           , &amp;
     &amp;                                             darkMatterParticle_                 =darkMatterParticleCDM_                            &amp;
     &amp;                                            )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterParticleSelfInteractingDarkMatterJiang_">
   <!-- The cross-section is reduced here since Jiang use a halo age of 10 Gyr instead of the current age of the Universe. By
        offsetting the cross-section by this factor we get the correct interaction radius.                                    -->
   <constructor>
    darkMatterParticleSelfInteractingDarkMatter       (                                                                                       &amp;
     &amp;                                             crossSectionSelfInteraction         =+1.0d0                                            &amp;
     &amp;                                                                                  *timeJiang                                        &amp;
     &amp;                                                                                  /cosmologyFunctionsPippin_%cosmicTime(1.0d0)    , &amp;
     &amp;                                             darkMatterParticle_                 =darkMatterParticleCDM_                            &amp;
     &amp;                                            )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfileDMOSIDMCoreNFW_"                  >
   <constructor>
    darkMatterProfileDMOSIDMCoreNFW                   (                                                                                       &amp;
     &amp;                                             factorRadiusCore                    =0.45d0                                          , &amp;
     &amp;                                             darkMatterHaloScale_                =darkMatterHaloScalePippin_                      , &amp;
     &amp;                                             darkMatterParticle_                 =darkMatterParticleSelfInteractingDarkMatter_      &amp;
     &amp;                                            )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfileDMONFWPippin_"                    >
   <constructor>
    darkMatterProfileDMONFW                           (                                                                                       &amp;
     &amp;                                             velocityDispersionUseSeriesExpansion=.false.                                         , &amp;
     &amp;                                             darkMatterHaloScale_                =darkMatterHaloScalePippin_                        &amp;
     &amp;                                            )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfileAdiabaticPippin_"                  >
   <constructor>
     darkMatterProfileAdiabaticGnedin2004             (                                                                                        &amp;
     &amp;                                             A                                   =0.85d0                                           , &amp;
     &amp;                                             omega                               =0.80d0                                           , &amp;
     &amp;                                             radiusFractionalPivot               =1.00d0                                           , &amp;
     &amp;                                             toleranceRelative                   =1.0d-2                                           , &amp;
     &amp;                                             nonAnalyticSolver                   =nonAnalyticSolversNumerical                      , &amp;
     &amp;                                             cosmologyParameters_                =cosmologyParametersPippin_                       , &amp;
     &amp;                                             darkMatterHaloScale_                =darkMatterHaloScalePippin_                       , &amp;
     &amp;                                             darkMatterProfileDMO_               =darkMatterProfileDMONFWPippin_                     &amp;
     &amp;                                            )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfileDMOSIDMIsothermal_"               >
   <constructor>
    darkMatterProfileDMOSIDMIsothermal                (                                                                                        &amp;
     &amp;                                             darkMatterProfileDMO_               =darkMatterProfileDMONFWPippin_                   , &amp;
     &amp;                                             darkMatterParticle_                 =darkMatterParticleSelfInteractingDarkMatterJiang_  &amp;
     &amp;                                            )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfileSIDMIsothermal_"                  >
   <constructor>
    darkMatterProfileSIDMIsothermal                   (                                                                                        &amp;
     &amp;                                             darkMatterProfile_                  =darkMatterProfileAdiabaticPippin_                , &amp;
     &amp;                                             darkMatterParticle_                 =darkMatterParticleSelfInteractingDarkMatterJiang_  &amp;
     &amp;                                            )
   </constructor>
  </referenceConstruct>
  !!]
  ! Test CoreNFW self-interacting dark matter profile.
  call Unit_Tests_Begin_Group('CoreNFW self-interacting dark matter profile')
  !! Set properties to match the Pippin halos of Elbert et al. (2015; https://ui.adsabs.harvard.edu/abs/2015MNRAS.453...29E).
  nodePippin => treeNode                    (                 )
  basic      => nodePippin%basic            (autoCreate=.true.)
  dmProfile  => nodePippin%darkMatterProfile(autoCreate=.true.)
  call basic%timeSet(cosmologyFunctionsPippin_%cosmicTime(1.0d0))
  call basic%massSet(massVirialPippin                           )
  radiusScale=+darkMatterHaloScalePippin_%radiusVirial(nodePippin) &
       &      /concentrationPippin                                          
  call dmProfile%scaleSet(radiusScale)
  massDistribution_ => darkMatterProfileDMOSIDMCoreNFW_%get(nodePippin)
  do i=1,7
     coordinates   =[radiusScale*radius(i),0.0d0,0.0d0]
     mass       (i)=massDistribution_%massEnclosedBySphere(radiusScale*radius(i))
     density    (i)=massDistribution_%density             (coordinates          )
  end do
  ! Interaction radius estimated from Figure 2 of Jiang et al. (2022).
  select type (massDistribution_)
  class is (massDistributionSphericalSIDM)
     call Assert(                                       &
          &      'interaction radius'                 , &
          &      massDistribution_%radiusInteraction(), &
          &      2.337664390096387d-3                 , &
          &      relTol=1.0d-2                          &
          &     )
  end select
  ! Mass and density computed using Mathematica.
  call Assert(                        &
       &      'enclosed mass'       , &
       &      mass                  , &
       &      [                       &
       &       8.129152247051000d+06, &
       &       5.187526955035893d+07, &
       &       2.497029836166847d+08, &
       &       7.849597988938884d+08, &
       &       1.782425795084028d+09, &
       &       3.340544103197336d+09, &
       &       5.399491464861386d+09  &
       &      ]                     , &
       &      relTol=1.0d-2           &
       &     )
  call Assert(                        &
       &      'density'             , &
       &      density               , &
       &      [                       &
       &       5.505559342899632d+16, &
       &       3.962252050514667d+16, &
       &       1.865412896231910d+16, &
       &       5.093277224421450d+15, &
       &       1.087639992726731d+15, &
       &       1.955797749926653d+14, &
       &       3.018205906859390d+13  &
       &      ]                     , &
       &      relTol=1.0d-2           &
       &     )
  call Unit_Tests_End_Group       ()
  !![
  <objectDestructor name="massDistribution_"/>
  !!]
  ! Test isothermal self-interacting dark matter profile.
  call Unit_Tests_Begin_Group('Isothermal self-interacting dark matter profile')
  !! Set properties to match the example halo generated by Fangzhou Jiang (private communication).
  nodeJiang  => treeNode                   (                 )
  basic      => nodeJiang%basic            (autoCreate=.true.)
  dmProfile  => nodeJiang%darkMatterProfile(autoCreate=.true.)
  call basic%timeSet(cosmologyFunctionsPippin_%cosmicTime(1.0d0))
  call basic%massSet(massVirialJiang                            )
  radiusScale=+darkMatterHaloScalePippin_%radiusVirial(nodeJiang) &
       &      /concentrationJiang                                          
  call dmProfile%scaleSet(radiusScale)
  massDistribution_       => darkMatterProfileDMOSIDMIsothermal_%get                   (nodeJiang)
  kinematicsDistribution_ => massDistribution_                  %kinematicsDistribution(         )
  !! Target values were provided by Fangzhou Jiang (private communication).
  select type (massDistribution_)
  class is (massDistributionSphericalSIDM)
     coordinates=[0.0d0,0.0d0,0.0d0]
     call Assert(                                                                                               &
          &      'interaction radius'                                                                         , &
          &      massDistribution_      %radiusInteraction   (                                               ), &
          &      6.9732d-3                                                                                    , &
          &      relTol=1.0d-2                                                                                  &
          &     )
     call Assert(                                                                                               &
          &      'central density'                                                                            , &
          &      massDistribution_      %density             (coordinates                                    ), &
          &      4.1168d16                                                                                    , &
          &      relTol=1.0d-1                                                                                  &
          &     )
     call Assert(                                                                                               &
          &      'central velocity dispersion'                                                                , &
          &      kinematicsDistribution_%velocityDispersion1D(coordinates,massDistribution_,massDistribution_), &
          &      54.9811d0                                                                                    , &
          &      relTol=1.0d-2                                                                                  &
          &     )
  end select
  call Unit_Tests_End_Group       ()
  !![
  <objectDestructor name="massDistribution_"      />
  <objectDestructor name="kinematicsDistribution_"/>
  !!]
  ! Test isothermal self-interacting dark matter profile with adiabatic contraction.
  call Unit_Tests_Begin_Group('Isothermal self-interacting dark matter profile (with adiabatic contraction)')
  !! Set properties to match the example halo generated by Fangzhou Jiang (private communication).
  call Unit_Tests_Begin_Group('Dark matter only case')
  basic      => nodeJiang%basic            (autoCreate=.true.)
  dmProfile  => nodeJiang%darkMatterProfile(autoCreate=.true.)
  call basic%timeSet(cosmologyFunctionsPippin_%cosmicTime(1.0d0))
  call basic%massSet(massVirialJiang                            )
  radiusScale=+darkMatterHaloScalePippin_%radiusVirial(nodeJiang) &
       &      /concentrationJiang                                          
  call dmProfile%scaleSet(radiusScale)
  massDistribution_       => darkMatterProfileSIDMIsothermal_%get                   (nodeJiang)
  kinematicsDistribution_ => massDistribution_               %kinematicsDistribution(         )
  !! Target values were provided by Fangzhou Jiang (private communication).
  select type (massDistribution_)
  class is (massDistributionSphericalSIDM)
     coordinates=[0.0d0,0.0d0,0.0d0]
     call Assert(                                                                                               &
          &      'interaction radius'                                                                         , &
          &      massDistribution_      %radiusInteraction   (                                               ), &
          &      6.9732d-3                                                                                    , &
          &      relTol=1.0d-2                                                                                  &
          &     )
     call Assert(                                                                                               &
          &      'central density'                                                                            , &
          &      massDistribution_      %density             (coordinates                                    ), &
          &      4.1168d16                                                                                    , &
          &      relTol=1.0d-1                                                                                  &
          &     )
     call Assert(                                                                                               &
          &      'central velocity dispersion'                                                                , &
          &      kinematicsDistribution_%velocityDispersion1D(coordinates,massDistribution_,massDistribution_), &
          &      54.9811d0                                                                                    , &
          &      relTol=1.0d-2                                                                                  &
          &     )
  end select
  call Unit_Tests_End_Group       ()
  !![
  <objectDestructor name="massDistribution_"      />
  <objectDestructor name="kinematicsDistribution_"/>
  !!]
  !! Insert a spheroid.
  call Unit_Tests_Begin_Group('With baryons case')
  spheroid      => nodeJiang%spheroid            (autoCreate=.true.)
  call spheroid%massStellarSet(fractionMassBaryonicJiang  *                           massVirialJiang                                                )
  call spheroid%     radiusSet(fractionRadiusHalfMassJiang*darkMatterHaloScalePippin_%radiusVirial   (nodeJiang)/radiusHalfMassDimensionlessHernquist)
  call Calculations_Reset(nodeJiang)
  massDistribution_       => darkMatterProfileSIDMIsothermal_%get                   (nodeJiang)
  kinematicsDistribution_ => massDistribution_               %kinematicsDistribution(         )
  !! Target values were measured from Figure A1 of Jiang et al. (2022).
  select type (massDistribution_)
  class is (massDistributionSphericalSIDM)
     coordinates=[0.0d0,0.0d0,0.0d0]
     call Assert(                                                                                               &
          &      'interaction radius'                                                                         , &
          &      massDistribution_      %radiusInteraction   (                                               ), &
          &      6.9732d-3                                                                                    , &
          &      relTol=1.0d-2                                                                                  &
          &     )
     call Assert(                                                                                               &
          &      'central density'                                                                            , &
          &      massDistribution_      %density             (coordinates                                    ), &
          &      2.534d17                                                                                     , &
          &      relTol=2.0d-1                                                                                  &
          &     )
     call Assert(                                                                                               &
          &      'central velocity dispersion'                                                                , &
          &      kinematicsDistribution_%velocityDispersion1D(coordinates,massDistribution_,massDistribution_), &
          &      59.139d0                                                                                     , &
          &      relTol=5.0d-2                                                                                  &
          &     )
  end select
  call Unit_Tests_End_Group()
  !![
  <objectDestructor name="massDistribution_"      />
  <objectDestructor name="kinematicsDistribution_"/>
  !!]
  call Unit_Tests_End_Group()
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
  ! Uninitialize node components.
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
  ! Clean up objects.
  !![
  <objectDestructor name="cosmologyParameters_"                             />
  <objectDestructor name="cosmologyFunctions_"                              />
  <objectDestructor name="virialDensityContrast_"                           />
  <objectDestructor name="darkMatterHaloScale_"                             />
  <objectDestructor name="darkMatterProfileDMOBurkert_"                     />
  <objectDestructor name="darkMatterProfileDMONFW_"                         />
  <objectDestructor name="darkMatterProfileDMONFWPippin_"                   />
  <objectDestructor name="darkMatterProfileDMONFWSeriesExpansion_"          />
  <objectDestructor name="darkMatterProfileDMOFiniteResolution_"            />
  <objectDestructor name="darkMatterProfileDMOSIDMCoreNFW_"                 />
  <objectDestructor name="darkMatterProfileDMOSIDMIsothermal_"              />
  <objectDestructor name="darkMatterProfileSIDMIsothermal_"                 />
  <objectDestructor name="darkMatterParticleSelfInteractingDarkMatter_"     />
  <objectDestructor name="darkMatterParticleSelfInteractingDarkMatterJiang_"/>
  <objectDestructor name="darkMatterParticleCDM_"                           />
  <objectDestructor name="cosmologyParametersPippin_"                       />
  <objectDestructor name="cosmologyFunctionsPippin_"                        />
  <objectDestructor name="virialDensityContrastPippin_"                     />
  <objectDestructor name="darkMatterHaloScalePippin_"                       />
  !!]
end program Test_Dark_Matter_Profiles
