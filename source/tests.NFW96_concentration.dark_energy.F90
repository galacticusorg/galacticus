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
Contains a program which tests the \cite{navarro_structure_1996} halo concentration algorithm in a dark energy Universe.
Comparisons are made to the {\normalfont \ttfamily charden} code written by Julio Navarro.
!!}

program Test_NFW96_Concentration_Dark_Energy
  !!{
  Tests the \cite{navarro_structure_1996} halo concentration algorithm in a dark energy Universe. Comparisons are made to the
  {\normalfont \ttfamily charden} code written by Julio Navarro.
  !!}
  use :: Cosmological_Density_Field          , only : criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt, cosmologicalMassVarianceFilteredPower
  use :: Cosmology_Functions                 , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters                , only : cosmologyParametersSimple                                   , hubbleUnitsLittleH
  use :: Dark_Matter_Particles               , only : darkMatterParticleCDM
  use :: Dark_Matter_Profiles_Concentration  , only : darkMatterProfileConcentrationNFW1996
  use :: Display                             , only : displayVerbositySet                                         , verbosityLevelStandard
  use :: Events_Hooks                        , only : eventsHooksInitialize
  use :: Functions_Global_Utilities          , only : Functions_Global_Set
  use :: Galacticus_Nodes                    , only : nodeClassHierarchyInitialize                                , nodeComponentBasic                   , treeNode
  use :: ISO_Varying_String                  , only : assignment(=)                                               , char                                 , operator(//)                       , varying_string              , &
       &                                              var_str
  use :: Input_Parameters                    , only : inputParameters
  use :: Linear_Growth                       , only : linearGrowthCollisionlessMatter
  use :: Node_Components                     , only : Node_Components_Initialize                                  , Node_Components_Thread_Initialize    , Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Power_Spectra                       , only : powerSpectrumStandard
  use :: Power_Spectra_Primordial            , only : powerSpectrumPrimordialPowerLaw
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredSimple
  use :: Power_Spectrum_Window_Functions     , only : powerSpectrumWindowFunctionTopHat
  use :: String_Handling                     , only : operator(//)                                                , String_Superscript
  use :: Transfer_Functions                  , only : transferFunctionBBKS
  use :: Unit_Tests                          , only : Assert                                                      , Unit_Tests_Begin_Group               , Unit_Tests_End_Group               , Unit_Tests_Finish
  use :: Virial_Density_Contrast             , only : virialDensityContrastFixed                                  , fixedDensityTypeCritical
  implicit none
  type            (treeNode                                                    ), pointer                 :: node
  class           (nodeComponentBasic                                          ), pointer                 :: basic
  integer                                                                       , dimension(6), parameter :: chardenLogHaloMass                 =[10,11,12,13,14,15]
  double precision                                                              , dimension(6), parameter :: chardenConcentrationZ0             =[10.2700200d00,9.0204391d00,7.8041310d00,6.6154380d00,5.4956946d00,4.4538398d00], &
       &                                                                                                     chardenConcentrationZ3             =[ 5.8715897d00,5.4417138d00,5.0239682d00,4.6186433d00,4.2366042d00,3.8884208d00]
  type            (cosmologyParametersSimple                                   )                          :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                              )                          :: cosmologyFunctions_
  type            (cosmologicalMassVarianceFilteredPower                       )                          :: cosmologicalMassVariance_
  type            (linearGrowthCollisionlessMatter                             )                          :: linearGrowth_
  type            (powerSpectrumWindowFunctionTopHat                           )                          :: powerSpectrumWindowFunction_
  type            (powerSpectrumPrimordialPowerLaw                             )                          :: powerSpectrumPrimordial_
  type            (transferFunctionBBKS                                        )                          :: transferFunction_
  type            (powerSpectrumPrimordialTransferredSimple                    )                          :: powerSpectrumPrimordialTransferred_
  type            (darkMatterParticleCDM                                       )                          :: darkMatterParticle_
  type            (criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt)                          :: criticalOverdensity_
  type            (virialDensityContrastFixed                                  )                          :: virialDensityContrast_
  type            (darkMatterProfileConcentrationNFW1996                       )                          :: darkMatterProfileConcentration_
  type            (varying_string                                              )                          :: message                                                                                                             , &
       &                                                                                                     parameterFile
  type            (inputParameters                                             )                          :: parameters
  integer                                                                                                 :: iMass
  double precision                                                                                        :: ourConcentration
  character       (len=2                                                       )                          :: label

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("NFW96 halo concentration algorithm: dark energy cosmology")
  ! Test NFW96 halo concentration algorithm in a dark energy universe.
  ! Read in controlling parameters.
  parameterFile='testSuite/parameters/NFW96HaloConcentration/darkEnergy.xml'
  parameters=inputParameters(parameterFile)
  call eventsHooksInitialize()
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  ! Create a node.
  node  => treeNode      (                 )
  ! Get the basic component.
  basic => node    %basic(autoCreate=.true.)
  ! Build required objects.
  !![
  <referenceConstruct object="darkMatterParticle_"                >
   <constructor>
    darkMatterParticleCDM                                        (                                                                         &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyParameters_"               >
   <constructor>
    cosmologyParametersSimple                                    (                                                                         &amp;
     &amp;                                                        OmegaMatter                        = 0.3d0                             , &amp;
     &amp;                                                        OmegaBaryon                        = 0.0296d0                          , &amp;
     &amp;                                                        OmegaDarkEnergy                    = 0.7d0                             , &amp;
     &amp;                                                        temperatureCMB                     = 2.7000d0                          , &amp;
     &amp;                                                        HubbleConstant                     =65.0d0                               &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctions_"                >
   <constructor>
    cosmologyFunctionsMatterLambda                               (                                                                         &amp;
     &amp;                                                        cosmologyParameters_               =cosmologyParameters_                 &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="linearGrowth_"                      >
   <constructor>
    linearGrowthCollisionlessMatter                              (                                                                         &amp;
     &amp;                                                        cosmologyParameters_               =cosmologyParameters_               , &amp;
     &amp;                                                        cosmologyFunctions_                =cosmologyFunctions_                  &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumPrimordial_"           >
   <constructor>
    powerSpectrumPrimordialPowerLaw                              (                                                                         &amp;
     &amp;                                                        index_                             =+1.000d0                           , &amp;
     &amp;                                                        running                            =+0.000d0                           , &amp;
     &amp;                                                        runningRunning                     =+0.000d0                           , &amp;
     &amp;                                                        wavenumberReference                =+1.000d0                           , &amp;
     &amp;                                                        runningSmallScalesOnly             =.false.                              &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="transferFunction_"                  >
   <constructor>
    transferFunctionBBKS                                         (                                                                         &amp;
     &amp;                                                        darkMatterParticle_                =darkMatterParticle_                , &amp;
     &amp;                                                        cosmologyParameters_               =cosmologyParameters_               , &amp;
     &amp;                                                        cosmologyFunctions_                =cosmologyFunctions_                  &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumPrimordialTransferred_">
   <constructor>
    powerSpectrumPrimordialTransferredSimple                     (                                                                         &amp;
     &amp;                                                        powerSpectrumPrimordial_           =powerSpectrumPrimordial_           , &amp;
     &amp;                                                        transferFunction_                  =transferFunction_                  , &amp;
     &amp;                                                        linearGrowth_                      =linearGrowth_                        &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumWindowFunction_"       >
   <constructor>
    powerSpectrumWindowFunctionTopHat                            (                                                                         &amp;
     &amp;                                                        cosmologyParameters_               =cosmologyParameters_                 &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologicalMassVariance_"          >
   <constructor>
    cosmologicalMassVarianceFilteredPower                        (                                                                         &amp;
     &amp;                                                        sigma8                             =0.900d+0                           , &amp;
     &amp;                                                        tolerance                          =1.000d-4                           , &amp;
     &amp;                                                        toleranceTopHat                    =1.000d-4                           , &amp;
     &amp;                                                        nonMonotonicIsFatal                =.true.                             , &amp;
     &amp;                                                        monotonicInterpolation             =.false.                            , &amp;
     &amp;                                                        truncateAtParticleHorizon          =.false.                            , &amp;
     &amp;                                                        cosmologyParameters_               =cosmologyParameters_               , &amp;
     &amp;                                                        cosmologyFunctions_                =cosmologyFunctions_                , &amp;
     &amp;                                                        linearGrowth_                      =linearGrowth_                      , &amp;
     &amp;                                                        powerSpectrumPrimordialTransferred_=powerSpectrumPrimordialTransferred_, &amp;
     &amp;                                                        powerSpectrumWindowFunction_       =powerSpectrumWindowFunction_         &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="criticalOverdensity_"               >
   <constructor>
    criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt(                                                                          &amp;
     &amp;                                                        linearGrowth_                      =linearGrowth_                      , &amp;
     &amp;                                                        cosmologyFunctions_                =cosmologyFunctions_                , &amp;
     &amp;                                                        cosmologicalMassVariance_          =cosmologicalMassVariance_          , &amp;
     &amp;                                                        darkMatterParticle_                =darkMatterParticle_                , &amp;
     &amp;                                                        tableStore                         =.true.                               &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>  
  <referenceConstruct object="virialDensityContrast_"               >
   <constructor>
    virialDensityContrastFixed(                                                                                                            &amp;
     &amp;                                                        densityContrastValue               =200.0d0                            , &amp;
     &amp;                                                        densityType                        =fixedDensityTypeCritical           , &amp;
     &amp;                                                        turnAroundOverVirialRadius         =  2.0d0                            , &amp;
     &amp;                                                        cosmologyParameters_               =cosmologyParameters_               , &amp;
     &amp;                                                        cosmologyFunctions_                =cosmologyFunctions_                  &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfileConcentration_"      >
   <constructor>
    darkMatterProfileConcentrationNFW1996                        (                                                                         &amp;
     &amp;                                                        f                                  =0.01d0                             , &amp;
     &amp;                                                        C                                  =3.00d3                             , &amp;
     &amp;                                                        cosmologyParameters_               =cosmologyParameters_               , &amp;
     &amp;                                                        cosmologyFunctions_                =cosmologyFunctions_                , &amp;
     &amp;                                                        criticalOverdensity_               =criticalOverdensity_               , &amp;
     &amp;                                                        cosmologicalMassVariance_          =cosmologicalMassVariance_          , &amp;
     &amp;                                                        virialDensityContrast_             =virialDensityContrast_               &amp;
     &amp;                                                       )
   </constructor>
  </referenceConstruct>
  !!]
  ! Loop over halo masses.
  do iMass=1,size(chardenLogHaloMass)
     ! Set the mass of the original node.
     call basic%massSet((10.0d0**chardenLogHaloMass(iMass))/cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH))
     ! Compute and compare concentration at z=0.
     call basic%timeSet(cosmologyFunctions_%cosmicTime(1.00d0))
     ourConcentration=darkMatterProfileConcentration_%concentration(node)
     write (label,'(i2)') chardenLogHaloMass(iMass)
     message=var_str("10")//String_Superscript(trim(adjustl(label)))//" M⊙/h halo concentration at z=0"
     call Assert(char(message),ourConcentration,chardenConcentrationZ0(iMass),relTol=0.02d0)
     ! Compute and compare concentration at z=3.
     call basic%timeSet(cosmologyFunctions_%cosmicTime(0.25d0))
     ourConcentration=darkMatterProfileConcentration_%concentration(node)
     message=var_str("10")//String_Superscript(trim(adjustl(label)))//" M⊙/h halo concentration at z=3"
     call Assert(char(message),ourConcentration,chardenConcentrationZ3(iMass),relTol=0.01d0)
  end do
  ! End unit tests.
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
end program Test_NFW96_Concentration_Dark_Energy
