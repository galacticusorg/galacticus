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
Contains a program which tests calculations of critical overdensity for collapse in warm dark matter model.
!!}

program Tests_Warm_Dark_Matter
  !!{
  Tests calculations of critical overdensity for collapse in warm dark matter model.
  !!}
  use :: Cosmology_Functions                 , only : cosmologyFunctionsMatterDarkEnergy
  use :: Cosmology_Parameters                , only : cosmologyParametersSimple
  use :: Cosmological_Density_Field          , only : criticalOverdensitySphericalCollapseClsnlssMttrDrkEnrgy, criticalOverdensityBarkana2001WDM, cosmologicalMassVarianceFilteredPower
  use :: Dark_Matter_Particles               , only : darkMatterParticleCDM                                  , darkMatterParticleWDMThermal
  use :: Display                             , only : displayVerbositySet                                    , verbosityLevelStandard
  use :: Events_Hooks                        , only : eventsHooksInitialize
  use :: Linear_Growth                       , only : linearGrowthCollisionlessMatter
  use :: Power_Spectra_Primordial            , only : powerSpectrumPrimordialPowerLaw
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredSimple
  use :: Power_Spectrum_Window_Functions     , only : powerSpectrumWindowFunctionTopHat
  use :: Transfer_Functions                  , only : transferFunctionEisensteinHu1999                       , transferFunctionBode2001                                                                , scaleCutOffModelBode2001
  use :: Unit_Tests                          , only : Assert                                                 , Unit_Tests_Begin_Group           , Unit_Tests_End_Group                 , Unit_Tests_Finish
  implicit none
  double precision                                                           , dimension(7) :: redshift                           =[0.0d0,1.0d00,3.0d00,7.0d00,15.0d00,31.0d0,63.0d0]
  double precision                                                           , dimension(5) :: mass                               =[1.0d9,1.0d10,1.0d11,1.0d12, 1.0d13              ]
  character       (len=1024                                                 )               :: message
  type            (cosmologyParametersSimple                                )               :: cosmologyParameters_
  type            (cosmologyFunctionsMatterDarkEnergy                       )               :: cosmologyFunctions_
  type            (linearGrowthCollisionlessMatter                          )               :: linearGrowth_
  type            (cosmologicalMassVarianceFilteredPower                    )               :: cosmologicalMassVariance_
  type            (powerSpectrumWindowFunctionTopHat                        )               :: powerSpectrumWindowFunction_
  type            (powerSpectrumPrimordialPowerLaw                          )               :: powerSpectrumPrimordial_
  type            (transferFunctionEisensteinHu1999                         )               :: transferFunctionCDM_
  type            (transferFunctionBode2001                                 )               :: transferFunction_
  type            (powerSpectrumPrimordialTransferredSimple                 )               :: powerSpectrumPrimordialTransferred_
  type            (darkMatterParticleCDM                                    )               :: darkMatterParticleCDM_
  type            (darkMatterParticleWDMThermal                             )               :: darkMatterParticle_
  type            (criticalOverdensitySphericalCollapseClsnlssMttrDrkEnrgy  )               :: criticalOverdensityCDM_
  type            (criticalOverdensityBarkana2001WDM                        )               :: criticalOverdensity_
  integer                                                                                   :: iExpansion                                                                            , iMass
  double precision                                                           , dimension(5) :: criticalOverdensityGradientMassValue                                                  , criticalOverdensityGradientMassNumeric, &
       &                                                                                       criticalOverdensityGradientTimeValue                                                  , criticalOverdensityGradientTimeNumeric
  double precision                                                                             age                                                                                   , expansionFactor                       , &
       &                                                                                       deltaT                                                                                , deltaMass

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Initialize event hooks.
  call eventsHooksInitialize()
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Critical overdensity for collapse: warm dark matter")
  ! Test critical overdensity for collapse in warm dark matter model.
  !![
  <referenceConstruct object="cosmologyParameters_"               >
   <constructor>
    cosmologyParametersSimple                                (                                                                             &amp;
     &amp;                                                    OmegaMatter                        = 0.3d0                                 , &amp;
     &amp;                                                    OmegaBaryon                        = 0.0d0                                 , &amp;
     &amp;                                                    OmegaDarkEnergy                    = 0.7d0                                 , &amp;
     &amp;                                                    temperatureCMB                     = 2.7d0                                 , &amp;
     &amp;                                                    HubbleConstant                     =70.0d0                                   &amp;
     &amp;                                                   )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctions_"                >
   <constructor>
    cosmologyFunctionsMatterDarkEnergy                       (                                                                             &amp;
     &amp;                                                    cosmologyParameters_               =cosmologyParameters_                  ,  &amp;
     &amp;                                                    darkEnergyEquationOfStateW0        =-1.0d0                                ,  &amp;
     &amp;                                                    darkEnergyEquationOfStateW1        =+0.0d0                                   &amp;
     &amp;                                                   )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterParticleCDM_"             >
   <constructor>
    darkMatterParticleCDM                                    (                                                                             &amp;
     &amp;                                                   )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterParticle_"                >
   <constructor>
    darkMatterParticleWDMThermal                             (                                                                             &amp;
     &amp;                                                    mass                               = 1.0d0                                 , &amp;
     &amp;                                                    degreesOfFreedomEffective          = 1.5d0                                 , &amp;
     &amp;                                                    cosmologyParameters_               =cosmologyParameters_                     &amp;
     &amp;                                                   )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="linearGrowth_"                      >
   <constructor>
    linearGrowthCollisionlessMatter                          (                                                                             &amp;
     &amp;                                                    cosmologyParameters_               =cosmologyParameters_                   , &amp;
     &amp;                                                    cosmologyFunctions_                =cosmologyFunctions_                      &amp;
     &amp;                                                   )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumPrimordial_"           >
   <constructor>
    powerSpectrumPrimordialPowerLaw                          (                                                                             &amp;
     &amp;                                                    index_                             =+1.0d0                                 , &amp;
     &amp;                                                    running                            =+0.0d0                                 , &amp;
     &amp;                                                    runningRunning                     =+0.0d0                                 , &amp;
     &amp;                                                    wavenumberReference                =+1.0d0                                 , &amp;
     &amp;                                                    runningSmallScalesOnly             =.false.                                  &amp;
     &amp;                                                   )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="transferFunctionCDM_"               >
   <constructor>
    transferFunctionEisensteinHu1999                         (                                                                             &amp;
     &amp;                                                    neutrinoNumberEffective            =3.046d0                                , &amp;
     &amp;                                                    neutrinoMassSummed                 =0.0d0                                  , &amp;
     &amp;                                                    darkMatterParticle_                =darkMatterParticleCDM_                 , &amp;
     &amp;                                                    cosmologyParameters_               =cosmologyParameters_                   , &amp;
     &amp;                                                    cosmologyFunctions_                =cosmologyFunctions_                      &amp;
     &amp;                                                   )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="transferFunction_"                  >
   <constructor>
    transferFunctionBode2001                                 (                                                                             &amp;
     &amp;                                                    transferFunctionCDM                =transferFunctionCDM_                   , &amp;
     &amp;                                                    scaleCutOffModel                   = scaleCutOffModelBode2001              , &amp;
     &amp;                                                    epsilon                            = 0.361d0                               , &amp;
     &amp;                                                    eta                                = 5.000d0                               , &amp;
     &amp;                                                    nu                                 = 1.200d0                               , &amp;
     &amp;                                                    time                               =13.462d0                               , &amp;
     &amp;                                                    cosmologyParameters_               =cosmologyParameters_                   , &amp;
     &amp;                                                    darkMatterParticle_                =darkMatterParticle_                    , &amp;
     &amp;                                                    cosmologyFunctions_                =cosmologyFunctions_                      &amp;
     &amp;                                                   )
  </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumPrimordialTransferred_">
   <constructor>
    powerSpectrumPrimordialTransferredSimple                 (                                                                             &amp;
     &amp;                                                    powerSpectrumPrimordial_           =powerSpectrumPrimordial_               , &amp;
     &amp;                                                    transferFunction_                  =transferFunction_                      , &amp;
     &amp;                                                    linearGrowth_                      =linearGrowth_                            &amp;
     &amp;                                                   )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumWindowFunction_"       >
   <constructor>
    powerSpectrumWindowFunctionTopHat                        (                                                                             &amp;
     &amp;                                                    cosmologyParameters_               =cosmologyParameters_                     &amp;
     &amp;                                                   )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologicalMassVariance_"          >
   <constructor>
    cosmologicalMassVarianceFilteredPower                    (                                                                             &amp;
     &amp;                                                    sigma8                             =0.8d0                                  , &amp;
     &amp;                                                    tolerance                          =1.0d-4                                 , &amp;
     &amp;                                                    toleranceTopHat                    =1.0d-4                                 , &amp;
     &amp;                                                    nonMonotonicIsFatal                =.true.                                 , &amp;
     &amp;                                                    monotonicInterpolation             =.false.                                , &amp;
     &amp;                                                    truncateAtParticleHorizon          =.false.                                , &amp;
     &amp;                                                    cosmologyParameters_               =cosmologyParameters_                   , &amp;
     &amp;                                                    cosmologyFunctions_                =cosmologyFunctions_                    , &amp;
     &amp;                                                    linearGrowth_                      =linearGrowth_                          , &amp;
     &amp;                                                    powerSpectrumPrimordialTransferred_=powerSpectrumPrimordialTransferred_    , &amp;
     &amp;                                                    powerSpectrumWindowFunction_       =powerSpectrumWindowFunction_             &amp;
     &amp;                                                   )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="criticalOverdensityCDM_"            >
   <constructor>
    criticalOverdensitySphericalCollapseClsnlssMttrDrkEnrgy  (                                                                             &amp;
     &amp;                                                    cosmologyFunctions_                =cosmologyFunctions_                    , &amp;
     &amp;                                                    linearGrowth_                      =linearGrowth_                          , &amp;
     &amp;                                                    cosmologicalMassVariance_          =cosmologicalMassVariance_              , &amp;
     &amp;                                                    darkMatterParticle_                =darkMatterParticleCDM_                 , &amp;
     &amp;                                                    normalization                      =1.0d0,                                   &amp;
     &amp;                                                    tableStore                         =.true.                                   &amp;
     &amp;                                                   )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="criticalOverdensity_"               >
   <constructor>
    criticalOverdensityBarkana2001WDM                        (                                                                             &amp;
     &amp;                                                    criticalOverdensityCDM             =criticalOverdensityCDM_                , &amp;
     &amp;                                                    cosmologyParameters_               =cosmologyParameters_                   , &amp;
     &amp;                                                    cosmologyFunctions_                =cosmologyFunctions_                    , &amp;
     &amp;                                                    cosmologicalMassVariance_          =cosmologicalMassVariance_              , &amp;
     &amp;                                                    darkMatterParticle_                =darkMatterParticle_                    , &amp;
     &amp;                                                    linearGrowth_                      =linearGrowth_                          , &amp;
     &amp;                                                    useFittingFunction                 =.true.                                   &amp;
     &amp;                                                   )
   </constructor>
  </referenceConstruct>
  !!]
  ! Gradient with respect to mass.
  do iExpansion=size(redshift),1,-1
     expansionFactor                                 =   cosmologyFunctions_ %expansionFactorFromRedshift    (     redshift       (iExpansion)                           )
     age                                             =   cosmologyFunctions_ %cosmicTime                     (     expansionFactor                                       )
     do iMass=1,size(mass)
        deltaMass                                    = 1.0d-9*mass(iMass)
        criticalOverdensityGradientMassValue  (iMass)=   criticalOverdensity_%gradientMass                   (time=age                        ,mass=mass(iMass)          )
        criticalOverdensityGradientMassNumeric(iMass)=+(                                                                                                                   &
             &                                          +criticalOverdensity_%value                          (time=age                        ,mass=mass(iMass)+deltaMass) &
             &                                          -criticalOverdensity_%value                          (time=age                        ,mass=mass(iMass)-deltaMass) &
             &                                         )                                                                                                                   &
             &                                        /2.0d0                                                                                                               &
             &                                        /deltaMass
     end do
     write (message,'(a,f6.1,a)') "gradient of critical density for collapse with respect to mass [z=",redshift(iExpansion),"]"
     call Assert(trim(message),criticalOverdensityGradientMassValue,criticalOverdensityGradientMassNumeric,absTol=1.0d-12,relTol=1.0d-3)
  end do
  ! Gradient with respect to time.
  do iExpansion=size(redshift),1,-1
     expansionFactor                                 =   cosmologyFunctions_ %expansionFactorFromRedshift    (     redshift       (iExpansion)                           )
     age                                             =   cosmologyFunctions_ %cosmicTime                     (     expansionFactor                                       )
     deltaT                                          = 1.0d-6*age
     do iMass=1,size(mass)
        criticalOverdensityGradientTimeValue  (iMass)=   criticalOverdensity_%gradientTime                   (time=age                        ,mass=mass(iMass)          )
        criticalOverdensityGradientTimeNumeric(iMass)=+(                                                                                                                   &
             &                                          +criticalOverdensity_%value                          (time=age+deltaT                 ,mass=mass(iMass)          ) &
             &                                          -criticalOverdensity_%value                          (time=age-deltaT                 ,mass=mass(iMass)          ) &
             &                                         )                                                                                                                   &
             &                                        /2.0d0                                                                                                               &
             &                                        /deltaT
     end do
     write (message,'(a,f6.1,a)') "gradient of critical density for collapse with respect to time [z=",redshift(iExpansion),"]"
     call Assert(trim(message),criticalOverdensityGradientTimeValue,criticalOverdensityGradientTimeNumeric,               relTol=1.0d-3)
  end do
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_Warm_Dark_Matter
