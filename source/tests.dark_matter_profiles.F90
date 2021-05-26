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

!% Contains a program which tests dark matter profiles.

program Test_Dark_Matter_Profiles
  !% Tests dark matter profiles.
  use :: Cosmology_Parameters        , only : cosmologyParametersSimple
  use :: Cosmology_Functions         , only : cosmologyFunctionsMatterLambda
  use :: Dark_Matter_Halo_Scales     , only : darkMatterHaloScaleVirialDensityContrastDefinition
  use :: Dark_Matter_Profiles_DMO    , only : darkMatterProfileDMOBurkert                                   , darkMatterProfileDMONFW, darkMatterProfileDMOFiniteResolution
  use :: Dark_Matter_Profiles_Generic, only : nonAnalyticSolversNumerical
  use :: Virial_Density_Contrast     , only : virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt
  use :: Events_Hooks                , only : eventsHooksInitialize
  use :: Functions_Global_Utilities  , only : Functions_Global_Set
  use :: Display                     , only : displayVerbositySet                                           , verbosityLevelStandard
  use :: Galacticus_Nodes            , only : nodeClassHierarchyFinalize                                    , nodeClassHierarchyInitialize     , nodeComponentBasic                 , nodeComponentDarkMatterProfile, &
          &                                   treeNode
  use :: Input_Parameters            , only : inputParameters
  use :: ISO_Varying_String          , only : varying_string                                                , assignment(=)                    , var_str
  use :: Node_Components             , only : Node_Components_Initialize                                    , Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Unit_Tests                  , only : Assert                                                        , Unit_Tests_Begin_Group           , Unit_Tests_End_Group               , Unit_Tests_Finish
  implicit none
  type            (treeNode                                                      ), pointer      :: node
  class           (nodeComponentBasic                                            ), pointer      :: basic
  class           (nodeComponentDarkMatterProfile                                ), pointer      :: dmProfile
  double precision                                                                , parameter    :: concentration                   =8.0d0                                                          , &
       &                                                                                            massVirial                      =1.0d0
  double precision                                                                , dimension(7) :: radius                          =[0.125d0, 0.250d0, 0.500d0, 1.000d0, 2.000d0, 4.000d0, 8.000d0]
  double precision                                                                , dimension(7) :: mass                                                                                            , &
       &                                                                                            density                                                                                         , &
       &                                                                                            fourier                                                                                         , &
       &                                                                                            radialVelocityDispersion                                                                        , &
       &                                                                                            radialVelocityDispersionSeriesExpansion
  type            (darkMatterProfileDMOBurkert                                   ), pointer      :: darkMatterProfileDMOBurkert_
  type            (darkMatterProfileDMONFW                                       ), pointer      :: darkMatterProfileDMONFW_
  type            (darkMatterProfileDMONFW                                       ), pointer      :: darkMatterProfileDMONFWSeriesExpansion_
  type            (darkMatterProfileDMOFiniteResolution                          ), pointer      :: darkMatterProfileDMOFiniteResolution_
  type            (cosmologyParametersSimple                                     ), pointer      :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                                ), pointer      :: cosmologyFunctions_
  type            (darkMatterHaloScaleVirialDensityContrastDefinition            ), pointer      :: darkMatterHaloScale_
  type            (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt), pointer      :: virialDensityContrast_
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
  allocate(cosmologyParameters_        )
  allocate(cosmologyFunctions_         )
  allocate(virialDensityContrast_      )
  allocate(darkMatterHaloScale_        )
  allocate(darkMatterProfileDMOBurkert_)
  allocate(darkMatterProfileDMONFW_    )
  allocate(darkMatterProfileDMONFWSeriesExpansion_  )
  allocate(darkMatterProfileDMOFiniteResolution_    )
  !# <referenceConstruct object="cosmologyParameters_"        >
  !#  <constructor>
  !#   cosmologyParametersSimple                                     (                                               &amp;
  !#    &amp;                                                         OmegaMatter           = 0.2815d0             , &amp;
  !#    &amp;                                                         OmegaBaryon           = 0.0465d0             , &amp;
  !#    &amp;                                                         OmegaDarkEnergy       = 0.7185d0             , &amp;
  !#    &amp;                                                         temperatureCMB        = 2.7800d0             , &amp;
  !#    &amp;                                                         HubbleConstant        =69.3000d0               &amp;
  !#    &amp;                                                        )
  !#  </constructor>
  !# </referenceConstruct>
  !# <referenceConstruct object="cosmologyFunctions_"         >
  !#  <constructor>
  !#   cosmologyFunctionsMatterLambda                                (                                               &amp;
  !#    &amp;                                                         cosmologyParameters_  =cosmologyParameters_    &amp;
  !#    &amp;                                                        )
  !#  </constructor>
  !# </referenceConstruct>
  !# <referenceConstruct object="virialDensityContrast_"      >
  !#  <constructor>
  !#   virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt(                                               &amp;
  !#    &amp;                                                         tableStore            =.true.                , &amp;
  !#    &amp;                                                         cosmologyFunctions_   =cosmologyFunctions_     &amp;
  !#    &amp;                                                        )
  !#  </constructor>
  !# </referenceConstruct>
  !# <referenceConstruct object="darkMatterHaloScale_"        >
  !#  <constructor>
  !#   darkMatterHaloScaleVirialDensityContrastDefinition            (                                               &amp;
  !#    &amp;                                                         cosmologyParameters_  =cosmologyParameters_  , &amp;
  !#    &amp;                                                         cosmologyFunctions_   =cosmologyFunctions_   , &amp;
  !#    &amp;                                                         virialDensityContrast_=virialDensityContrast_  &amp;
  !#    &amp;                                                        )
  !#  </constructor>
  !# </referenceConstruct>
  ! Create a node.
  node      => treeNode                  (                 )
  ! Create components.
  basic     => node    %basic            (autoCreate=.true.)
  dmProfile => node    %darkMatterProfile(autoCreate=.true.)
  ! Set properties.
  call basic%timeSet(cosmologyFunctions_%cosmicTime(1.0d0))
  call basic%massSet(massVirial                           )
  ! Compute scale radius.
  radiusScale=+darkMatterHaloScale_%virialRadius(node) &
       &      /concentration
  call dmProfile%scaleSet(radiusScale)
  ! Build dark matter profiles.
  !# <referenceConstruct object="darkMatterProfileDMOBurkert_"           >
  !#  <constructor>
  !#   darkMatterProfileDMOBurkert         (                                                                  &amp;
  !#    &amp;                               darkMatterHaloScale_                =darkMatterHaloScale_         &amp;
  !#    &amp;                              )
  !#  </constructor>
  !# </referenceConstruct>
  !# <referenceConstruct object="darkMatterProfileDMONFW_"               >
  !#  <constructor>
  !#   darkMatterProfileDMONFW             (                                                                  &amp;
  !#    &amp;                               velocityDispersionUseSeriesExpansion=.false.                    , &amp;
  !#    &amp;                               darkMatterHaloScale_                =darkMatterHaloScale_         &amp;
  !#    &amp;                              )
  !#  </constructor>
  !# </referenceConstruct>
  !# <referenceConstruct object="darkMatterProfileDMONFWSeriesExpansion_">
  !#  <constructor>
  !#   darkMatterProfileDMONFW             (                                                                  &amp;
  !#    &amp;                               velocityDispersionUseSeriesExpansion=.true.                     , &amp;
  !#    &amp;                               darkMatterHaloScale_                =darkMatterHaloScale_         &amp;
  !#    &amp;                              )
  !#  </constructor>
  !# </referenceConstruct>  
  !# <referenceConstruct object="darkMatterProfileDMOFiniteResolution_"  >
  !#  <constructor>
  !#   darkMatterProfileDMOFiniteResolution(                                                                  &amp;
  !#    &amp;                               lengthResolution                    =0.5d0*radiusScale          , &amp;
  !#    &amp;                               massResolution                      =0.0d0                      , &amp;
  !#    &amp;                               resolutionIsComoving                =.false.                    , &amp;
  !#    &amp;                               nonAnalyticSolver                   =nonAnalyticSolversNumerical, &amp;
  !#    &amp;                               darkMatterProfileDMO_               =darkMatterProfileDMONFW_   , &amp;
  !#    &amp;                               darkMatterHaloScale_                =darkMatterHaloScale_       , &amp;
  !#    &amp;                               cosmologyFunctions_                 =cosmologyFunctions_          &amp;
  !#    &amp;                              )
  !#  </constructor>
  !# </referenceConstruct>
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
  radialVelocityDispersion               =radialVelocityDispersion               /darkMatterHaloScale_%virialVelocity(node)
  radialVelocityDispersionSeriesExpansion=radialVelocityDispersionSeriesExpansion/darkMatterHaloScale_%virialVelocity(node)
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
  ! End unit tests.
  call Unit_Tests_End_Group       ()
  call Unit_Tests_Finish          ()
  ! Uninitialize node components.
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
  ! Clean up objects.
  !# <objectDestructor name="cosmologyParameters_"                   />
  !# <objectDestructor name="cosmologyFunctions_"                    />
  !# <objectDestructor name="virialDensityContrast_"                 />
  !# <objectDestructor name="darkMatterHaloScale_"                   />
  !# <objectDestructor name="darkMatterProfileDMOBurkert_"           />
  !# <objectDestructor name="darkMatterProfileDMONFW_"               />
  !# <objectDestructor name="darkMatterProfileDMONFWSeriesExpansion_"/>
  !# <objectDestructor name="darkMatterProfileDMOFiniteResolution_"  />
end program Test_Dark_Matter_Profiles
