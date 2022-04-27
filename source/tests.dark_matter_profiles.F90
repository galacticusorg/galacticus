!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
  use :: Calculations_Resets         , only : Calculations_Reset
  use :: Cosmology_Parameters        , only : cosmologyParametersSimple
  use :: Cosmology_Functions         , only : cosmologyFunctionsMatterLambda
  use :: Dark_Matter_Particles       , only : darkMatterParticleSelfInteractingDarkMatter                   , darkMatterParticleCDM
  use :: Dark_Matter_Halo_Scales     , only : darkMatterHaloScaleVirialDensityContrastDefinition
  use :: Dark_Matter_Profiles_DMO    , only : darkMatterProfileDMOBurkert                                   , darkMatterProfileDMONFW             , darkMatterProfileDMOFiniteResolution, darkMatterProfileDMOSIDMCoreNFW, &
       &                                      darkMatterProfileDMOSIDMIsothermal
  use :: Dark_Matter_Profiles        , only : darkMatterProfileSIDMIsothermal                               , darkMatterProfileAdiabaticGnedin2004
  use :: Dark_Matter_Profiles_Generic, only : nonAnalyticSolversNumerical
  use :: Galactic_Structure          , only : galacticStructureStandard
  use :: Virial_Density_Contrast     , only : virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt, virialDensityContrastFixed          , fixedDensityTypeCritical
  use :: Events_Hooks                , only : eventsHooksInitialize
  use :: Functions_Global_Utilities  , only : Functions_Global_Set
  use :: Display                     , only : displayVerbositySet                                           , verbosityLevelStandard
  use :: Galacticus_Nodes            , only : nodeClassHierarchyFinalize                                    , nodeClassHierarchyInitialize        , nodeComponentBasic                   , nodeComponentDarkMatterProfile, &
          &                                   treeNode                                                      , nodeComponentSpheroid
  use :: Input_Parameters            , only : inputParameters
  use :: ISO_Varying_String          , only : varying_string                                                , assignment(=)                       , var_str
  use :: Node_Components             , only : Node_Components_Initialize                                    , Node_Components_Thread_Initialize   , Node_Components_Thread_Uninitialize  , Node_Components_Uninitialize
  use :: Unit_Tests                  , only : Assert                                                        , Unit_Tests_Begin_Group              , Unit_Tests_End_Group                 , Unit_Tests_Finish
  implicit none
  type            (treeNode                                                      ), pointer      :: node                                                                                                , &
       &                                                                                            nodePippin                                                                                          , &
       &                                                                                            nodeJiang
  class           (nodeComponentBasic                                            ), pointer      :: basic
  class           (nodeComponentDarkMatterProfile                                ), pointer      :: dmProfile
  class           (nodeComponentSpheroid                                         ), pointer      :: spheroid
  double precision                                                                , parameter    :: concentration                       = 8.0d0                                                         , &
       &                                                                                            massVirial                          = 1.0d0                                                         , &
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
  double precision                                                                , dimension(7) :: radius                              =[0.125d0, 0.250d0, 0.500d0, 1.000d0, 2.000d0, 4.000d0, 8.000d0]
  double precision                                                                , dimension(7) :: mass                                                                                                , &
       &                                                                                            density                                                                                             , &
       &                                                                                            fourier                                                                                             , &
       &                                                                                            radialVelocityDispersion                                                                            , &
       &                                                                                            radialVelocityDispersionSeriesExpansion
  type            (darkMatterParticleCDM                                         ), pointer      :: darkMatterParticleCDM_
  type            (darkMatterParticleSelfInteractingDarkMatter                   ), pointer      :: darkMatterParticleSelfInteractingDarkMatter_                                                        , &
       &                                                                                            darkMatterParticleSelfInteractingDarkMatterJiang_
  type            (darkMatterProfileDMOBurkert                                   ), pointer      :: darkMatterProfileDMOBurkert_
  type            (darkMatterProfileDMONFW                                       ), pointer      :: darkMatterProfileDMONFW_                                                                            , &
       &                                                                                            darkMatterProfileDMONFWPippin_
  type            (darkMatterProfileDMONFW                                       ), pointer      :: darkMatterProfileDMONFWSeriesExpansion_
  type            (darkMatterProfileDMOFiniteResolution                          ), pointer      :: darkMatterProfileDMOFiniteResolution_
  type            (darkMatterProfileDMOSIDMCoreNFW                               ), pointer      :: darkMatterProfileDMOSIDMCoreNFW_
  type            (darkMatterProfileDMOSIDMIsothermal                            ), pointer      :: darkMatterProfileDMOSIDMIsothermal_
  type            (darkMatterProfileAdiabaticGnedin2004                          ), pointer      :: darkMatterProfileAdiabaticPippin_
  type            (darkMatterProfileSIDMIsothermal                               ), pointer      :: darkMatterProfileSIDMIsothermal_
  type            (galacticStructureStandard                                     ), pointer      :: galacticStructureStandard_
  type            (cosmologyParametersSimple                                     ), pointer      :: cosmologyParameters_                                                                                , &
       &                                                                                            cosmologyParametersPippin_
  type            (cosmologyFunctionsMatterLambda                                ), pointer      :: cosmologyFunctions_                                                                                 , &
       &                                                                                            cosmologyFunctionsPippin_
  type            (darkMatterHaloScaleVirialDensityContrastDefinition            ), pointer      :: darkMatterHaloScale_                                                                                , &
       &                                                                                            darkMatterHaloScalePippin_
  type            (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt), pointer      :: virialDensityContrast_
  type            (virialDensityContrastFixed                                    ), pointer      :: virialDensityContrastPippin_
  type            (inputParameters                                               )               :: parameters
  integer                                                                                        :: i
  double precision                                                                               :: radiusScale

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
  allocate(darkMatterProfileDMONFW_                         )
  allocate(darkMatterProfileDMONFWPippin_                   )
  allocate(darkMatterProfileDMONFWSeriesExpansion_          )
  allocate(darkMatterProfileDMOFiniteResolution_            )
  allocate(darkMatterProfileDMOSIDMCoreNFW_                 )
  allocate(darkMatterProfileDMOSIDMIsothermal_              )
  allocate(darkMatterProfileSIDMIsothermal_                 )
  allocate(darkMatterProfileAdiabaticPippin_                )
  allocate(galacticStructureStandard_                       )
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
  radiusScale=+darkMatterHaloScale_%radiusVirial(node) &
       &      /concentration
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
     &amp;                               darkMatterHaloScale_                =darkMatterHaloScale_       , &amp;
     &amp;                               cosmologyFunctions_                 =cosmologyFunctions_          &amp;
     &amp;                              )
   </constructor>
  </referenceConstruct>
  !!]
  ! Begin unit tests.
  call Unit_Tests_Begin_Group('Dark matter profiles')
  ! Test Burkert profile.
  call Unit_Tests_Begin_Group('Burkert profile')
  do i=1,7
     mass   (i)=darkMatterProfileDMOBurkert_%enclosedMass(node,      radiusScale*radius(i))
     density(i)=darkMatterProfileDMOBurkert_%density     (node,      radiusScale*radius(i))*radiusScale**3
     fourier(i)=darkMatterProfileDMOBurkert_%kSpace      (node,1.0d0/radiusScale/radius(i))
  end do
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
  call Unit_Tests_End_Group       ()
  ! Test NFW profile.
  call Unit_Tests_Begin_Group('NFW profile')
  do i=1,7
     mass                                   (i)=darkMatterProfileDMONFW_               %enclosedMass            (node,      radiusScale*radius(i))
     density                                (i)=darkMatterProfileDMONFW_               %density                 (node,      radiusScale*radius(i))*radiusScale**3
     fourier                                (i)=darkMatterProfileDMONFW_               %kSpace                  (node,1.0d0/radiusScale/radius(i))
     radialVelocityDispersion               (i)=darkMatterProfileDMONFW_               %radialVelocityDispersion(node,      radiusScale*radius(i))
     radialVelocityDispersionSeriesExpansion(i)=darkMatterProfileDMONFWSeriesExpansion_%radialVelocityDispersion(node,      radiusScale*radius(i))
  end do
  ! Radial velocity dispersion in units of virial velocity.
  radialVelocityDispersion               =radialVelocityDispersion               /darkMatterHaloScale_%velocityVirial(node)
  radialVelocityDispersionSeriesExpansion=radialVelocityDispersionSeriesExpansion/darkMatterHaloScale_%velocityVirial(node)
  call Assert(                        &
       &      'enclosed mass'       , &
       &      mass                  , &
       &      [                       &
       &       5.099550982355504d-3 , &
       &       1.768930674181593d-2 , &
       &       5.513246746363203d-2 , &
       &       1.476281525188409d-1 , &
       &       3.301489257042704d-1 , &
       &       6.186775455118112d-1 , &
       &       1.000000000000000d+0   &
       &      ]                     , &
       &      relTol=1.0d-6           &
       &     )
  call Assert(                        &
       &      'density'             , &
       &      density               , &
       &      [                       &
       &       3.844641857939078d-1 , &
       &       1.557079952465327d-1 , &
       &       5.406527612726829d-2 , &
       &       1.520585891079421d-2 , &
       &       3.379079757954268d-3 , &
       &       6.082343564317684d-4 , &
       &       9.386332660984080d-5   &
       &      ]                     , &
       &      relTol=1.0d-6           &
       &     )
  call Assert(                        &
       &      'fourier'             , &
       &      fourier               , &
       &      [                       &
       &       1.099052711906094d-2 , &
       &       3.746284433831412d-2 , &
       &       1.127728075593513d-1 , &
       &       2.619030036621630d-1 , &
       &       5.434560723827092d-1 , &
       &       8.448608065763800d-1 , &
       &       9.579597044271230d-1   &
       &      ]                     , &
       &      relTol=1.0d-6           &
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
  do i=1,7
     mass   (i)=darkMatterProfileDMOFiniteResolution_%enclosedMass(node,      radiusScale*radius(i))
     density(i)=darkMatterProfileDMOFiniteResolution_%density     (node,      radiusScale*radius(i))*radiusScale**3
     fourier(i)=darkMatterProfileDMOFiniteResolution_%kSpace      (node,1.0d0/radiusScale/radius(i))
  end do
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
  <referenceConstruct object="galacticStructureStandard_"                        >
   <constructor>
    galacticStructureStandard                         (                                                                                        &amp;
     &amp;                                             cosmologyFunctions_                 =cosmologyFunctionsPippin_                        , &amp;
     &amp;                                             darkMatterHaloScale_                =darkMatterHaloScalePippin_                       , &amp;
     &amp;                                             darkMatterProfile_                  =darkMatterProfileSIDMIsothermal_                   &amp;
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
     &amp;                                             nonAnalyticSolver                   =nonAnalyticSolversNumerical                      , &amp;
     &amp;                                             cosmologyParameters_                =cosmologyParametersPippin_                       , &amp;
     &amp;                                             darkMatterHaloScale_                =darkMatterHaloScalePippin_                       , &amp;
     &amp;                                             darkMatterProfileDMO_               =darkMatterProfileDMONFWPippin_                   , &amp;
     &amp;                                             galacticStructure_                  =galacticStructureStandard_                         &amp;
     &amp;                                            )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfileDMOSIDMIsothermal_"               >
   <constructor>
    darkMatterProfileDMOSIDMIsothermal                (                                                                                        &amp;
     &amp;                                             darkMatterProfileDMO_               =darkMatterProfileDMONFWPippin_                   , &amp;
     &amp;                                             darkMatterHaloScale_                =darkMatterHaloScalePippin_                       , &amp;
     &amp;                                             darkMatterParticle_                 =darkMatterParticleSelfInteractingDarkMatterJiang_  &amp;
     &amp;                                            )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfileSIDMIsothermal_"                  >
   <constructor>
    darkMatterProfileSIDMIsothermal                   (                                                                                        &amp;
     &amp;                                             darkMatterProfile_                  =darkMatterProfileAdiabaticPippin_                , &amp;
     &amp;                                             darkMatterHaloScale_                =darkMatterHaloScalePippin_                       , &amp;
     &amp;                                             darkMatterParticle_                 =darkMatterParticleSelfInteractingDarkMatterJiang_, &amp;
     &amp;                                             galacticStructure_                  =galacticStructureStandard_                         &amp;
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
  do i=1,7
     mass   (i)=darkMatterProfileDMOSIDMCoreNFW_%enclosedMass(nodePippin,radiusScale*radius(i))
     density(i)=darkMatterProfileDMOSIDMCoreNFW_%density     (nodePippin,radiusScale*radius(i))
  end do
  ! Interaction radius estimated from Figure 2 of Jiang et al. (2022).
  call Assert(                                                                &
       &      'interaction radius'                                          , &
       &      darkMatterProfileDMOSIDMCoreNFW_%radiusInteraction(nodePippin), &
       &      2.337664390096387d-3                                          , &
       &      relTol=1.0d-2                                                   &
       &     )
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
  !! Target values were provided by Fangzhou Jiang (private communication).
  call Assert(                                                                                      &
       &      'interaction radius'                                                                , &
       &      darkMatterProfileDMOSIDMIsothermal_%radiusInteraction       (nodeJiang             ), &
       &      6.9732d-3                                                                           , &
       &      relTol=1.0d-2                                                                         &
       &     )
  call Assert(                                                                                      &
       &      'central density'                                                                   , &
       &      darkMatterProfileDMOSIDMIsothermal_%density                 (nodeJiang,radius=0.0d0), &
       &      4.1168d16                                                                           , &
       &      relTol=1.0d-1                                                                         &
       &     )
  call Assert(                                                                                      &
       &      'central velocity dispersion'                                                       , &
       &      darkMatterProfileDMOSIDMIsothermal_%radialVelocityDispersion(nodeJiang,radius=0.0d0), &
       &      54.9811d0                                                                           , &
       &      relTol=1.0d-2                                                                         &
       &     )
  call Unit_Tests_End_Group       ()
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
  !! Target values were provided by Fangzhou Jiang (private communication).
  call Assert(                                                                                   &
       &      'interaction radius'                                                             , &
       &      darkMatterProfileSIDMIsothermal_%radiusInteraction       (nodeJiang             ), &
       &      6.9732d-3                                                                        , &
       &      relTol=1.0d-2                                                                      &
       &     )
  call Assert(                                                                                   &
       &      'central density'                                                                , &
       &      darkMatterProfileSIDMIsothermal_%density                 (nodeJiang,radius=0.0d0), &
       &      4.1168d16                                                                        , &
       &      relTol=1.0d-1                                                                      &
       &     )
  call Assert(                                                                                   &
       &      'central velocity dispersion'                                                    , &
       &      darkMatterProfileSIDMIsothermal_%radialVelocityDispersion(nodeJiang,radius=0.0d0), &
       &      54.9811d0                                                                        , &
       &      relTol=1.0d-2                                                                      &
       &     )
  call Unit_Tests_End_Group       ()
  !! Insert a spheroid.
  call Unit_Tests_Begin_Group('With baryons case')
  spheroid      => nodeJiang%spheroid            (autoCreate=.true.)
  call spheroid%massStellarSet(fractionMassBaryonicJiang  *                           massVirialJiang                                                )
  call spheroid%     radiusSet(fractionRadiusHalfMassJiang*darkMatterHaloScalePippin_%radiusVirial   (nodeJiang)/radiusHalfMassDimensionlessHernquist)
  call Calculations_Reset(nodeJiang)
  !! Target values were measured from Figure A1 of Jiang et al. (2022).
  call Assert(                                                                                   &
       &      'interaction radius'                                                             , &
       &      darkMatterProfileSIDMIsothermal_%radiusInteraction       (nodeJiang             ), &
       &      6.9732d-3                                                                        , &
       &      relTol=1.0d-2                                                                      &
       &     )
  call Assert(                                                                                   &
       &      'central density'                                                                , &
       &      darkMatterProfileSIDMIsothermal_%density                 (nodeJiang,radius=0.0d0), &
       &      2.534d17                                                                         , &
       &      relTol=2.0d-1                                                                      &
       &     )
  call Assert(                                                                                   &
       &      'central velocity dispersion'                                                    , &
       &      darkMatterProfileSIDMIsothermal_%radialVelocityDispersion(nodeJiang,radius=0.0d0), &
       &      59.139d0                                                                         , &
       &      relTol=5.0d-2                                                                      &
       &     )
  call Unit_Tests_End_Group()
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
  <objectDestructor name="galacticStructureStandard_"                       />
  <objectDestructor name="darkMatterParticleSelfInteractingDarkMatter_"     />
  <objectDestructor name="darkMatterParticleSelfInteractingDarkMatterJiang_"/>
  <objectDestructor name="darkMatterParticleCDM_"                           />
  <objectDestructor name="cosmologyParametersPippin_"                       />
  <objectDestructor name="cosmologyFunctionsPippin_"                        />
  <objectDestructor name="virialDensityContrastPippin_"                     />
  <objectDestructor name="darkMatterHaloScalePippin_"                       />
  !!]
end program Test_Dark_Matter_Profiles
