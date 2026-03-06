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
Contains a program which tests spherical collapse calculations for a dark energy Universe, specifically using an Einstein-de
Sitter cosmology (e.g. \citealt{kitayama_semianalytic_1996}; eqn.~A2).
!!}

program Tests_Spherical_Collapse_Dark_Energy_EdS
  !!{
  Tests spherical collapse calculations for a dark energy Universe, specifically using an Einstein-de
  Sitter cosmology. Compares results to the analytic solution.
  !!}
  use :: Cosmology_Functions                 , only : cosmologyFunctionsMatterDarkEnergy
  use :: Cosmology_Parameters                , only : cosmologyParametersSimple
  use :: Cosmological_Density_Field          , only : criticalOverdensitySphericalCollapseClsnlssMttrDrkEnrgy  , cosmologicalMassVarianceFilteredPower
  use :: Dark_Matter_Particles               , only : darkMatterParticleCDM
  use :: Display                             , only : displayVerbositySet                                      , verbosityLevelStandard
  use :: Events_Hooks                        , only : eventsHooksInitialize
  use :: Linear_Growth                       , only : linearGrowthCollisionlessMatter
  use :: Numerical_Constants_Math            , only : Pi
  use :: Power_Spectra_Primordial            , only : powerSpectrumPrimordialPowerLaw
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredSimple
  use :: Power_Spectrum_Window_Functions     , only : powerSpectrumWindowFunctionTopHat
  use :: Spherical_Collapse_Solvers          , only : cllsnlssMttrDarkEnergyFixedAtTurnaround
  use :: Transfer_Functions                  , only : transferFunctionEisensteinHu1999
  use :: Unit_Tests                          , only : Assert                                                   , Unit_Tests_Begin_Group               , Unit_Tests_End_Group, Unit_Tests_Finish
  use :: Virial_Density_Contrast             , only : virialDensityContrastSphericalCollapseClsnlssMttrDrkEnrgy
  implicit none
  double precision                                                           , dimension(7) :: redshift                           =[0.0d0,1.0d0,3.0d0,7.0d0,15.0d0,31.0d0,63.0d0]
  double precision                                                           , parameter    :: mass                               =1.0d12
  character       (len=1024                                                 )               :: message
  type            (cosmologyParametersSimple                                )               :: cosmologyParameters_
  type            (cosmologyFunctionsMatterDarkEnergy                       )               :: cosmologyFunctions_
  type            (linearGrowthCollisionlessMatter                          )               :: linearGrowth_
  type            (cosmologicalMassVarianceFilteredPower                    )               :: cosmologicalMassVariance_
  type            (powerSpectrumWindowFunctionTopHat                        )               :: powerSpectrumWindowFunction_
  type            (powerSpectrumPrimordialPowerLaw                          )               :: powerSpectrumPrimordial_
  type            (transferFunctionEisensteinHu1999                         )               :: transferFunction_
  type            (powerSpectrumPrimordialTransferredSimple                 )               :: powerSpectrumPrimordialTransferred_
  type            (darkMatterParticleCDM                                    )               :: darkMatterParticle_
  type            (virialDensityContrastSphericalCollapseClsnlssMttrDrkEnrgy)               :: virialDensityContrast_
  type            (criticalOverdensitySphericalCollapseClsnlssMttrDrkEnrgy  )               :: criticalOverdensity_
  integer                                                                                   :: iExpansion
  double precision                                                                          :: age                                                                               , criticalOverdensityValue     , &
       &                                                                                       criticalOverdensityExpected                                                       , expansionFactor              , &
       &                                                                                       virialDensityContrastActual                                                       , virialDensityContrastExpected

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Initialize event hooks.
  call eventsHooksInitialize()
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Spherical collapse: dark energy solver (Einstein-de Sitter cosmology)")
  ! Test spherical collapse in a flat universe.
  !![
  <referenceConstruct object="darkMatterParticle_"                >
   <constructor>
    darkMatterParticleCDM                                    (                                                                                  &amp;
     &amp;                                                   )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyParameters_"               >
   <constructor>
    cosmologyParametersSimple                                (                                                                                  &amp;
     &amp;                                                    OmegaMatter                             = 1.0d0                                 , &amp;
     &amp;                                                    OmegaBaryon                             = 0.0d0                                 , &amp;
     &amp;                                                    OmegaDarkEnergy                         = 0.0d0                                 , &amp;
     &amp;                                                    temperatureCMB                          = 2.7d0                                 , &amp;
     &amp;                                                    HubbleConstant                          =70.0d0                                   &amp;
     &amp;                                                   )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctions_"                >
   <constructor>
    cosmologyFunctionsMatterDarkEnergy                       (                                                                                  &amp;
     &amp;                                                    cosmologyParameters_                    =cosmologyParameters_                  ,  &amp;
     &amp;                                                    darkEnergyEquationOfStateW0             =-1.0d0                                ,  &amp;
     &amp;                                                    darkEnergyEquationOfStateW1             =+0.0d0                                   &amp;
     &amp;                                                   )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="linearGrowth_"                      >
   <constructor>
    linearGrowthCollisionlessMatter                          (                                                                                  &amp;
     &amp;                                                    cosmologyParameters_                    =cosmologyParameters_                   , &amp;
     &amp;                                                    cosmologyFunctions_                     =cosmologyFunctions_                      &amp;
     &amp;                                                   )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumPrimordial_"           >
   <constructor>
    powerSpectrumPrimordialPowerLaw                          (                                                                                  &amp;
     &amp;                                                    index_                                  =+1.0d0                                 , &amp;
     &amp;                                                    running                                 =+0.0d0                                 , &amp;
     &amp;                                                    runningRunning                          =+0.0d0                                 , &amp;
     &amp;                                                    wavenumberReference                     =+1.0d0                                 , &amp;
     &amp;                                                    runningSmallScalesOnly                  =.false.                                  &amp;
     &amp;                                                   )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="transferFunction_"                  >
   <constructor>
    transferFunctionEisensteinHu1999                         (                                                                                  &amp;
     &amp;                                                    neutrinoNumberEffective                 =3.046d0                                , &amp;
     &amp;                                                    neutrinoMassSummed                      =0.0d0                                  , &amp;
     &amp;                                                    darkMatterParticle_                     =darkMatterParticle_                    , &amp;
     &amp;                                                    cosmologyParameters_                    =cosmologyParameters_                   , &amp;
     &amp;                                                    cosmologyFunctions_                     =cosmologyFunctions_                      &amp;
     &amp;                                                   )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumPrimordialTransferred_">
   <constructor>
    powerSpectrumPrimordialTransferredSimple                 (                                                                                  &amp;
     &amp;                                                    powerSpectrumPrimordial_                =powerSpectrumPrimordial_               , &amp;
     &amp;                                                    transferFunction_                       =transferFunction_                      , &amp;
     &amp;                                                    linearGrowth_                           =linearGrowth_                            &amp;
     &amp;                                                   )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumWindowFunction_"       >
   <constructor>
    powerSpectrumWindowFunctionTopHat                        (                                                                                  &amp;
     &amp;                                                    cosmologyParameters_                    =cosmologyParameters_                     &amp;
     &amp;                                                   )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologicalMassVariance_"          >
   <constructor>
    cosmologicalMassVarianceFilteredPower                    (                                                                                  &amp;
     &amp;                                                    sigma8                                  =0.8d0                                  , &amp;
     &amp;                                                    tolerance                               =1.0d-4                                 , &amp;
     &amp;                                                    toleranceTopHat                         =1.0d-4                                 , &amp;
     &amp;                                                    rootVarianceLogarithmicGradientTolerance=1.0d-9                                 , &amp;
     &amp;                                                    integrationFailureIsFatal               =.true.                                 , &amp;
     &amp;                                                    storeTabulations                        =.true.                                 , &amp;
     &amp;                                                    nonMonotonicIsFatal                     =.true.                                 , &amp;
     &amp;                                                    monotonicInterpolation                  =.false.                                , &amp;
     &amp;                                                    truncateAtParticleHorizon               =.false.                                , &amp;
     &amp;                                                    cosmologyParameters_                    =cosmologyParameters_                   , &amp;
     &amp;                                                    cosmologyFunctions_                     =cosmologyFunctions_                    , &amp;
     &amp;                                                    linearGrowth_                           =linearGrowth_                          , &amp;
     &amp;                                                    powerSpectrumPrimordialTransferred_     =powerSpectrumPrimordialTransferred_    , &amp;
     &amp;                                                    powerSpectrumWindowFunction_            =powerSpectrumWindowFunction_             &amp;
     &amp;                                                   )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="criticalOverdensity_"               >
   <constructor>
    criticalOverdensitySphericalCollapseClsnlssMttrDrkEnrgy  (                                                                                  &amp;
     &amp;                                                    cosmologyFunctions_                     =cosmologyFunctions_                    , &amp;
     &amp;                                                    linearGrowth_                           =linearGrowth_                          , &amp;
     &amp;                                                    cosmologicalMassVariance_               =cosmologicalMassVariance_              , &amp;
     &amp;                                                    darkMatterParticle_                     =darkMatterParticle_                    , &amp;
     &amp;                                                    normalization                           =1.0d0                                  , &amp;
     &amp;                                                    tableStore                              =.true.                                   &amp;
     &amp;                                                   )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="virialDensityContrast_"             >
   <constructor>
    virialDensityContrastSphericalCollapseClsnlssMttrDrkEnrgy(                                                                                  &amp;
     &amp;                                                    cosmologyFunctions_                     =cosmologyFunctions_                    , &amp;
     &amp;                                                    energyFixedAt                           =cllsnlssMttrDarkEnergyFixedAtTurnaround, &amp;
     &amp;                                                    tableStore                              =.true.                                   &amp;
     &amp;                                                   )
   </constructor>
  </referenceConstruct>
  !!]
  do iExpansion=1,size(redshift)
     expansionFactor            =cosmologyFunctions_ %expansionFactorFromRedshift(redshift       (iExpansion))
     age                        =cosmologyFunctions_ %cosmicTime                 (expansionFactor            )
     criticalOverdensityValue   =criticalOverdensity_%value                      (age                        )
     criticalOverdensityExpected=3.0d0*(12.0d0*Pi)**(2.0d0/3.0d0)/20.0d0
     write (message,'(a,f6.1,a,f6.4,a)') "critical density for collapse [z=",redshift(iExpansion),";Ωₘ=",cosmologyFunctions_%omegaMatterEpochal(age),"]"
     call Assert(trim(message),criticalOverdensityValue,criticalOverdensityExpected,relTol=3.0d-5)
     virialDensityContrastActual  =virialDensityContrast_%densityContrast(mass,age)
     virialDensityContrastExpected=18.0d0*Pi**2
     write (message,'(a,f6.1,a,f6.4,a)') "virial density contrast       [z=",redshift(iExpansion),";Ωₘ=",cosmologyFunctions_%omegaMatterEpochal(age),"]"
     call Assert(trim(message),virialDensityContrastActual,virialDensityContrastExpected,relTol=3.0d-5)
  end do
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_Spherical_Collapse_Dark_Energy_EdS
