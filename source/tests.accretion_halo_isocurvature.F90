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
Contains a program which tests the correction to halo accretion due to isocurvature perturbations from \cite{jessop_ripples_2026}.
!!}

program Tests_Accretion_Halo_Isocurvature
  !!{
  Tests the correction to halo accretion due to isocurvature perturbations from \cite{jessop_ripples_2026}. Compares results for
  halos of mass $10^{11}\mathrm{M}_\odot$ to their stated result \citep[][eqn.~9]{jessop_ripples_2026}.
  !!}
  use :: Accretion_Halos                     , only : accretionHaloIsocurvature
  use :: Cosmology_Parameters                , only : cosmologyParametersSimple
  use :: Cosmology_Functions                 , only : cosmologyFunctionsMatterLambda
  use :: Cosmological_Density_Field          , only : criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt, cosmologicalMassVarianceFilteredPower
  use :: Dark_Matter_Particles               , only : darkMatterParticleCDM
  use :: Display                             , only : displayVerbositySet                                         , verbosityLevelStandard
  use :: Events_Hooks                        , only : eventsHooksInitialize
  use :: Functions_Global_Utilities          , only : Functions_Global_Set
  use :: Galacticus_Nodes                    , only : nodeClassHierarchyInitialize                                , nodeComponentBasic                   , treeNode                               , nodeClassHierarchyFinalize  , &
       &                                              nodeComponentDarkMatterProfile
  use :: Node_Components                     , only : Node_Components_Initialize                                  , Node_Components_Thread_Initialize    , Node_Components_Thread_Uninitialize    , Node_Components_Uninitialize
  use :: Input_Paths                         , only : inputPath                                                   , pathTypeExec
  use :: Input_Parameters                    , only : inputParameters
  use :: Linear_Growth                       , only : linearGrowthCollisionlessMatter
  use :: Power_Spectra_Primordial            , only : powerSpectrumPrimordialPowerLaw
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredSimple
  use :: Power_Spectrum_Window_Functions     , only : powerSpectrumWindowFunctionTopHat
  use :: Transfer_Functions                  , only : transferFunctionEisensteinHu1999
  use :: Unit_Tests                          , only : Assert                                                      , Unit_Tests_Begin_Group               , Unit_Tests_End_Group, Unit_Tests_Finish
  use :: ISO_Varying_String                  , only : assignment(=)                                               , operator(//)                         , varying_string
  implicit none
  type            (treeNode                                                    ), pointer :: node
  class           (nodeComponentBasic                                          ), pointer :: basic
  type            (accretionHaloIsocurvature                                   ), pointer :: accretionHalo_                    , accretionHaloNull => null()
  type            (cosmologyParametersSimple                                   ), pointer :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                              ), pointer :: cosmologyFunctions_
  type            (linearGrowthCollisionlessMatter                             ), pointer :: linearGrowth_
  type            (cosmologicalMassVarianceFilteredPower                       ), pointer :: cosmologicalMassVariance_
  type            (powerSpectrumWindowFunctionTopHat                           ), pointer :: powerSpectrumWindowFunction_
  type            (powerSpectrumPrimordialPowerLaw                             ), pointer :: powerSpectrumPrimordial_
  type            (transferFunctionEisensteinHu1999                            ), pointer :: transferFunction_
  type            (powerSpectrumPrimordialTransferredSimple                    ), pointer :: powerSpectrumPrimordialTransferred_
  type            (darkMatterParticleCDM                                       ), pointer :: darkMatterParticle_
  type            (criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt), pointer :: criticalOverdensity_
  double precision                                                                        :: fractionRelative
  type            (inputParameters                                             )          :: parameters
  type            (varying_string                                              )          :: parameterFile

  ! Initialize.
  parameterFile=inputPath(pathTypeExec)//'testSuite/parameters/accretionHaloIsocurvature.xml'
  parameters   =inputParameters(parameterFile)
  call displayVerbositySet              (verbosityLevelStandard)
  call eventsHooksInitialize            (                      )
  call Functions_Global_Set             (                      )
  call nodeClassHierarchyInitialize     (parameters            )
  call Node_Components_Initialize       (parameters            )
  call Node_Components_Thread_Initialize(parameters            )
  ! Construct all required objects.
  allocate(accretionHalo_                     )
  allocate(cosmologyParameters_               )
  allocate(cosmologyFunctions_                )
  allocate(linearGrowth_                      )
  allocate(cosmologicalMassVariance_          )
  allocate(powerSpectrumWindowFunction_       )
  allocate(powerSpectrumPrimordial_           )
  allocate(transferFunction_                  )
  allocate(powerSpectrumPrimordialTransferred_)
  allocate(darkMatterParticle_                )
  allocate(criticalOverdensity_               )
  !![
  <referenceConstruct object="darkMatterParticle_"                >
   <constructor>
    darkMatterParticleCDM                                       (                                                                             &amp;
     &amp;                                                      )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyParameters_"               >
   <constructor>
    cosmologyParametersSimple                                   (                                                                             &amp;
     &amp;                                                       OmegaMatter                        = 0.3060d0                              , &amp;
     &amp;                                                       OmegaBaryon                        = 0.0486d0                              , &amp;
     &amp;                                                       OmegaDarkEnergy                    = 0.6940d0                              , &amp;
     &amp;                                                       temperatureCMB                     = 2.7255d0                              , &amp;
     &amp;                                                       HubbleConstant                     =68.1000d0                                &amp;
     &amp;                                                      )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctions_"                >
   <constructor>
    cosmologyFunctionsMatterLambda                              (                                                                             &amp;
     &amp;                                                       cosmologyParameters_               =cosmologyParameters_                     &amp;
     &amp;                                                      )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="linearGrowth_"                      >
   <constructor>
    linearGrowthCollisionlessMatter                             (                                                                             &amp;
     &amp;                                                       cosmologyParameters_               =cosmologyParameters_                   , &amp;
     &amp;                                                       cosmologyFunctions_                =cosmologyFunctions_                      &amp;
     &amp;                                                      )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumPrimordial_"           >
   <constructor>
    powerSpectrumPrimordialPowerLaw                             (                                                                             &amp;
     &amp;                                                       index_                             =+0.967d0                               , &amp;
     &amp;                                                       running                            =+0.000d0                               , &amp;
     &amp;                                                       runningRunning                     =+0.000d0                               , &amp;
     &amp;                                                       wavenumberReference                =+1.000d0                               , &amp;
     &amp;                                                       runningSmallScalesOnly             =.false.                                  &amp;
     &amp;                                                      )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="transferFunction_"                  >
   <constructor>
    transferFunctionEisensteinHu1999                            (                                                                             &amp;
     &amp;                                                       neutrinoNumberEffective            =3.046d0                                , &amp;
     &amp;                                                       neutrinoMassSummed                 =0.0d0                                  , &amp;
     &amp;                                                       darkMatterParticle_                =darkMatterParticle_                    , &amp;
     &amp;                                                       cosmologyParameters_               =cosmologyParameters_                   , &amp;
     &amp;                                                       cosmologyFunctions_                =cosmologyFunctions_                      &amp;
     &amp;                                                      )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumPrimordialTransferred_">
   <constructor>
    powerSpectrumPrimordialTransferredSimple                    (                                                                             &amp;
     &amp;                                                       powerSpectrumPrimordial_           =powerSpectrumPrimordial_               , &amp;
     &amp;                                                       transferFunction_                  =transferFunction_                      , &amp;
     &amp;                                                       linearGrowth_                      =linearGrowth_                            &amp;
     &amp;                                                      )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumWindowFunction_"       >
   <constructor>
    powerSpectrumWindowFunctionTopHat                           (                                                                             &amp;
     &amp;                                                       cosmologyParameters_               =cosmologyParameters_                     &amp;
     &amp;                                                      )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologicalMassVariance_"          >
   <constructor>
    cosmologicalMassVarianceFilteredPower                       (                                                                             &amp;
     &amp;                                                       sigma8                             =0.807d+0                               , &amp;
     &amp;                                                       tolerance                          =1.000d-4                               , &amp;
     &amp;                                                       toleranceTopHat                    =1.000d-4                               , &amp;
     &amp;                                                       nonMonotonicIsFatal                =.true.                                 , &amp;
     &amp;                                                       monotonicInterpolation             =.false.                                , &amp;
     &amp;                                                       truncateAtParticleHorizon          =.false.                                , &amp;
     &amp;                                                       cosmologyParameters_               =cosmologyParameters_                   , &amp;
     &amp;                                                       cosmologyFunctions_                =cosmologyFunctions_                    , &amp;
     &amp;                                                       linearGrowth_                      =linearGrowth_                          , &amp;
     &amp;                                                       powerSpectrumPrimordialTransferred_=powerSpectrumPrimordialTransferred_    , &amp;
     &amp;                                                       powerSpectrumWindowFunction_       =powerSpectrumWindowFunction_             &amp;
     &amp;                                                      )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="criticalOverdensity_"               >
   <constructor>
    criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt(                                                                             &amp;
     &amp;                                                       cosmologyFunctions_                =cosmologyFunctions_                    , &amp;
     &amp;                                                       linearGrowth_                      =linearGrowth_                          , &amp;
     &amp;                                                       cosmologicalMassVariance_          =cosmologicalMassVariance_              , &amp;
     &amp;                                                       darkMatterParticle_                =darkMatterParticle_                    , &amp;
     &amp;                                                       normalization                      =1.0d0,                                   &amp;
     &amp;                                                       tableStore                         =.true.                                   &amp;
     &amp;                                                      )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="accretionHalo_"                     >
   <constructor>
    accretionHaloIsocurvature                                   (                                                                             &amp;
     &amp;                                                       countPerDecade                     =0                                      , &amp;
     &amp;                                                       accretionHalo_                     =accretionHaloNull                      , &amp;
     &amp;                                                       cosmologyParameters_               =cosmologyParameters_                   , &amp;
     &amp;                                                       criticalOverdensity_               =criticalOverdensity_                   , &amp;
     &amp;                                                       linearGrowth_                      =linearGrowth_                            &amp;
     &amp;                                                      )
   </constructor>
  </referenceConstruct>
  !!]
  ! Build a node at the present day with a mass of 10¹¹M☉.
  node             => treeNode                   (                 )
  basic             => node    %basic            (autoCreate=.true.)
  call basic%timeSet            (cosmologyFunctions_%cosmicTime(1.0d00))
  call basic%timeLastIsolatedSet(cosmologyFunctions_%cosmicTime(1.0d00))
  call basic%massSet            (                               1.0d11 )
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Accretion onto halos: isocurvature perturbation correction")
  ! Get the accreted fraction (relative to the universal baryon fraction).
  fractionRelative=accretionHalo_%fraction(node)
  call Assert('Relative baryon fraction for 10¹¹M☉ halo at z=0',fractionRelative,1.0d0-0.0069d0,relTol=1.0d-4)
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
  ! Clean up.
  call parameters%reset  ()
  call parameters%destroy()
  !![
  <objectDestructor name="accretionHalo_"                     />
  <objectDestructor name="cosmologyParameters_"               />
  <objectDestructor name="cosmologyFunctions_"                />
  <objectDestructor name="linearGrowth_"                      />
  <objectDestructor name="cosmologicalMassVariance_"          />
  <objectDestructor name="powerSpectrumWindowFunction_"       />
  <objectDestructor name="powerSpectrumPrimordial_"           />
  <objectDestructor name="transferFunction_"                  />
  <objectDestructor name="powerSpectrumPrimordialTransferred_"/>
  <objectDestructor name="darkMatterParticle_"                />
  <objectDestructor name="criticalOverdensity_"               />
  !!]
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
end program Tests_Accretion_Halo_Isocurvature
