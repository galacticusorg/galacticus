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

!+    Contributions to this file made by: Xiaolong Du, Andrew Benson.

!!{
Contains a program which tests the calculation of dark matter halo radius enclosing a given mass.
!!}

program Test_Dark_Matter_Halo_Radius_Enclosing_Mass
  !!{
  Tests the calculation of dark matter halo radius enclosing a given mass.
  !!}
  use :: Cosmology_Parameters      , only : cosmologyParametersSimple
  use :: Cosmology_Functions       , only : cosmologyFunctionsMatterLambda
  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScale                                           , darkMatterHaloScaleClass                , darkMatterHaloScaleVirialDensityContrastDefinition
  use :: Virial_Density_Contrast   , only : virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt
  use :: Dark_Matter_Profiles_DMO  , only : darkMatterProfileDMOBurkert                                   , darkMatterProfileDMOClass               , darkMatterProfileDMOHeated                             , darkMatterProfileDMONFW       , &
          &                                 darkMatterProfileDMOTruncated                                 , darkMatterProfileDMOTruncatedExponential, darkMatterProfileHeatingTidal
  use :: Mass_Distributions        , only : nonAnalyticSolversFallThrough                                 , massDistributionClass
  use :: Display                   , only : displayVerbositySet                                           , verbosityLevelStandard
  use :: Events_Hooks              , only : eventsHooksInitialize
  use :: Functions_Global_Utilities, only : Functions_Global_Set
  use :: Error                     , only : Error_Report
  use :: Galacticus_Nodes          , only : nodeClassHierarchyFinalize                                    , nodeClassHierarchyInitialize            , nodeComponentBasic                                     , nodeComponentDarkMatterProfile, &
          &                                 nodeComponentSatellite                                        , treeNode
  use :: ISO_Varying_String        , only : assignment(=)                                                 , varying_string
  use :: Input_Parameters          , only : inputParameters
  use :: Node_Components           , only : Node_Components_Initialize                                    , Node_Components_Thread_Initialize       , Node_Components_Thread_Uninitialize                    , Node_Components_Uninitialize
  use :: Unit_Tests                , only : Assert                                                        , Unit_Tests_Begin_Group                  , Unit_Tests_End_Group                                   , Unit_Tests_Finish
  implicit none
  type            (treeNode                                                      )               :: node
  class           (nodeComponentBasic                                            ), pointer      :: basic
  class           (nodeComponentSatellite                                        ), pointer      :: satellite
  class           (nodeComponentDarkMatterProfile                                ), pointer      :: dmProfile
  class           (massDistributionClass                                         ), pointer      :: massDistribution_
  class           (darkMatterProfileDMOClass                                     ), pointer      :: darkMatterProfileDMO_
  type            (darkMatterProfileDMONFW                                       ), target       :: darkMatterProfileDMONFW_
  type            (darkMatterProfileDMOBurkert                                   ), target       :: darkMatterProfileDMOBurkert_
  type            (darkMatterProfileDMOTruncated                                 ), target       :: darkMatterProfileDMOTruncated_
  type            (darkMatterProfileDMOTruncatedExponential                      ), target       :: darkMatterProfileDMOTruncatedExponential_
  type            (darkMatterProfileDMOHeated                                    ), target       :: darkMatterProfileDMOHeated_
  type            (darkMatterProfileHeatingTidal                                 )               :: darkMatterProfileHeatingTidal_
  type            (cosmologyParametersSimple                                     )               :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                                )               :: cosmologyFunctions_
  type            (darkMatterHaloScaleVirialDensityContrastDefinition            )               :: darkMatterHaloScale_
  type            (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt)               :: virialDensityContrast_
  double precision                                                                , dimension(7) :: radiusOverVirialRadius                    =[0.125d0 , 0.250d0, 0.500d0, 1.000d0, 2.000d0, 4.000d0, 8.000d0]
  double precision                                                                , dimension(7) :: radius                                              , radiusRoot
  double precision                                                                , dimension(7) :: mass
  double precision                                                                , parameter    :: radiusFractionalTruncateMinimum           = 2.00d+00, radiusFractionalTruncateMaximum=8.0d0
  double precision                                                                , parameter    :: time                                      =13.80d+00
  double precision                                                                , parameter    :: massVirial                                = 1.00d+10, concentration                  =8.0d0
  double precision                                                                               :: radiusFractionalDecay                     = 0.06d+00
  double precision                                                                , parameter    :: heatingSpecific                           = 1.00d+06
  double precision                                                                , parameter    :: coefficientSecondOrder                    = 0.00d+00
  double precision                                                                , parameter    :: correlationVelocityRadius                 =-1.00d+00
  double precision                                                                , parameter    :: toleranceRelativeVelocityDispersion       = 1.00d-06
  double precision                                                                , parameter    :: toleranceRelativeVelocityDispersionMaximum= 1.00d-03
  double precision                                                                , parameter    :: toleranceRelativePotential                = 1.00d-03
  double precision                                                                , parameter    :: fractionRadiusFinalSmall                  = 1.00d-03
  logical                                                                         , parameter    :: tolerateVelocityMaximumFailure            =.false.
  logical                                                                         , parameter    :: toleratePotentialIntegrationFailure       =.false.
  logical                                                                         , parameter    :: tolerateEnclosedMassIntegrationFailure    =.false.
  logical                                                                         , parameter    :: tolerateVelocityDispersionFailure         =.false.
  double precision                                                                               :: radiusVirial                                        , radiusScale                           , &
       &                                                                                            toleranceRelative
  type            (varying_string                                                )               :: parameterFile
  type            (inputParameters                                               )               :: parameters
  integer                                                                                        :: i                                                  , j
  logical                                                                         , parameter    :: velocityDispersionUseSeriesExpansion     =.true.   , velocityDispersionApproximate  =.true.
  logical                                                                                        :: limitToVirialRadius

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group('Dark matter halo radius enclosing a give mass')
  ! Read in controlling parameters.
  parameterFile='testSuite/parameters/darkMatterHaloRadiusEnclosingMass.xml'
  parameters=inputParameters(parameterFile)
  ! Initialize event hooks and global functions.
  call eventsHooksInitialize()
  call Functions_Global_Set ()
  ! Initialize node components.
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  ! Create the dark matter profiles.
  !![
  <referenceConstruct object="cosmologyParameters_"  >
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
     &amp;                                                         tableStore            =.true.                , &amp;
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
  darkMatterProfileDMONFW_                  =  darkMatterProfileDMONFW                 (velocityDispersionUseSeriesExpansion,                                                                                              darkMatterHaloScale_ )
  darkMatterProfileDMOBurkert_              =  darkMatterProfileDMOBurkert             (                                                                                                                                   darkMatterHaloScale_ )
  darkMatterProfileDMOTruncated_            =  darkMatterProfileDMOTruncated           (radiusFractionalTruncateMinimum     ,radiusFractionalTruncateMaximum,                                                                                     &
       &                                                                                nonAnalyticSolversFallThrough       ,                                                                    darkMatterProfileDMONFW_ ,darkMatterHaloScale_ )
  darkMatterProfileDMOTruncatedExponential_ =  darkMatterProfileDMOTruncatedExponential(radiusFractionalDecay                                                                                                             ,                       &
       &                                                                                nonAnalyticSolversFallThrough       ,                                                                    darkMatterProfileDMONFW_ ,darkMatterHaloScale_ )
  darkMatterProfileHeatingTidal_            =  darkMatterProfileHeatingTidal           (coefficientSecondOrder              ,coefficientSecondOrder         ,coefficientSecondOrder             ,correlationVelocityRadius                      )
  darkMatterProfileDMOHeated_               =  darkMatterProfileDMOHeated              (nonAnalyticSolversFallThrough       ,velocityDispersionApproximate  ,tolerateEnclosedMassIntegrationFailure,tolerateVelocityDispersionFailure,tolerateVelocityMaximumFailure,toleratePotentialIntegrationFailure,fractionRadiusFinalSmall,toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum,toleranceRelativePotential,darkMatterProfileDMONFW_,                  &
       &                                                                                darkMatterProfileHeatingTidal_                                                                                                                          )
  ! Set up the node.
  basic     => node%basic                 (autoCreate=.true.)
  satellite => node%satellite             (autoCreate=.true.)
  dmProfile => node%darkMatterProfile     (autoCreate=.true.)
  call basic    %massSet                  (massVirial       )
  call basic    %timeSet                  (time             )
  call basic    %timeLastIsolatedSet      (time             )
  call satellite%tidalHeatingNormalizedSet(heatingSpecific  )
  ! Get the virial radius.
  radiusVirial=darkMatterHaloScale_%radiusVirial(node)
  ! Compute scale radius.
  radiusScale =radiusVirial/concentration
  call dmProfile%scaleSet(radiusScale)
  ! Test different dark matter profiles.
  radius      =radiusOverVirialRadius*radiusVirial
  do i=1,5
     limitToVirialRadius=.false.
     toleranceRelative  =1.0d-6
     select case (i)
     case (1)
        call Unit_Tests_Begin_Group('NFW profile'                       )
        darkMatterProfileDMO_ => darkMatterProfileDMONFW_
     case (2)
        call Unit_Tests_Begin_Group('Burkert profile'                   )
        darkMatterProfileDMO_ => darkMatterProfileDMOBurkert_
     case (3)
        call Unit_Tests_Begin_Group('Truncated profile'                 )
        darkMatterProfileDMO_ => darkMatterProfileDMOTruncated_
     case (4)
        call Unit_Tests_Begin_Group('Exponentially truncated profile'   )
        darkMatterProfileDMO_ => darkMatterProfileDMOTruncatedExponential_
        limitToVirialRadius   =  .true.
     case (5)
        call Unit_Tests_Begin_Group('Heated profile'                    )
        darkMatterProfileDMO_ => darkMatterProfileDMOHeated_
        toleranceRelative     =  3.0d-5
     case default
        call Error_Report('unknown profile'//{introspection:location})
     end select
     massDistribution_ => darkMatterProfileDMO_%get(node)
     do j=1,7
        mass      (j)=massDistribution_%massEnclosedBySphere(radius(j))
        radiusRoot(j)=massDistribution_%radiusEnclosingMass (mass  (j))
        if (limitToVirialRadius .and. radiusOverVirialRadius(j) > 1.0d0) radiusRoot(j)=radius(j)
     end do
     !![
     <objectDestructor name="massDistribution_"/>
     !!]
     call Assert('radius enclosing a given mass',radius,radiusRoot,relTol=toleranceRelative)
     call Unit_Tests_End_Group()
  end do
  ! End unit tests.
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  ! Uninitialize node components.
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
end program Test_Dark_Matter_Halo_Radius_Enclosing_Mass
