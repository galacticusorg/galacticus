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

program Tests_Sigma
  !!{
  Tests calculations of the mass variance, $\sigma(M)$.
  !!}
  use :: Cosmological_Density_Field          , only : cosmologicalMassVarianceFilteredPower
  use :: Cosmology_Functions                 , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters                , only : cosmologyParametersSimple               , hubbleUnitsLittleH
  use :: Dark_Matter_Particles               , only : darkMatterParticleCDM
  use :: Display                             , only : displayVerbositySet                     , verbosityLevelStandard
  use :: Error                               , only : Error_Handler_Register
  use :: Linear_Growth                       , only : linearGrowthCollisionlessMatter
  use :: Numerical_Constants_Math            , only : Pi
  use :: Numerical_Ranges                    , only : Make_Range                              , rangeTypeLogarithmic
  use :: Power_Spectra_Primordial            , only : powerSpectrumPrimordialPowerLaw
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredSimple
  use :: Power_Spectrum_Window_Functions     , only : powerSpectrumWindowFunctionSharpKSpace  , powerSpectrumWindowFunctionTopHat
  use :: Transfer_Functions                  , only : transferFunctionIdentity                , transferFunctionEisensteinHu1999
  use :: Unit_Tests                          , only : Assert                                  , Unit_Tests_Begin_Group           , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  integer                                                   , parameter            :: massCount                             =10
  double precision                                          , parameter            :: massMaximum                           =1.0d15, massMinimum                              =1.0d6
  double precision                                          , dimension(massCount) :: mass                                         , massFromSigma                                  , &
       &                                                                              sigma
  type            (darkMatterParticleCDM                   )                       :: darkMatterParticleCDM_
  type            (cosmologyParametersSimple               )                       :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda          )                       :: cosmologyFunctions_
  type            (linearGrowthCollisionlessMatter         )                       :: linearGrowth_
  type            (cosmologicalMassVarianceFilteredPower   )                       :: cosmologicalMassVarianceLCDM_                , cosmologicalMassVarianceFilteredPower_
  type            (powerSpectrumWindowFunctionSharpKSpace  )                       :: powerSpectrumWindowFunctionSharpKSpace_
  type            (powerSpectrumWindowFunctionTopHat       )                       :: powerSpectrumWindowFunctionLCDM_
  type            (powerSpectrumPrimordialPowerLaw         )                       :: powerSpectrumPrimordialLCDM_                 , powerSpectrumPrimordialPowerLaw_
  type            (transferFunctionIdentity                )                       :: transferFunctionIdentity_
  type            (transferFunctionEisensteinHu1999        )                       :: transferFunctionLCDM_
  type            (powerSpectrumPrimordialTransferredSimple)                       :: powerSpectrumPrimordialTransferredLCDM_      , powerSpectrumPrimordialTransferredSimple_
  integer                                                                          :: iMass
  double precision                                                                 :: mass8                                        , radius8                                        , &
       &                                                                              sigma8

  ! Initialize error handling.
  call Error_Handler_Register()
  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Power spectrum: σ(M)")
  ! Create an array of masses.
  mass=Make_Range(massMinimum,massMaximum,massCount,rangeType=rangeTypeLogarithmic)
  ! Construct required objects.
  !![
  <referenceConstruct object="darkMatterParticleCDM_"                 >
   <constructor>
    darkMatterParticleCDM                    (                                                                            &amp;
     &amp;                                   )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyParameters_"                   >
   <constructor>
    cosmologyParametersSimple                (                                                                            &amp;
     &amp;                                    OmegaMatter                        = 0.31530d0                            , &amp;
     &amp;                                    OmegaBaryon                        = 0.04930d0                            , &amp;
     &amp;                                    OmegaDarkEnergy                    = 0.68470d0                            , &amp;
     &amp;                                    temperatureCMB                     = 2.72548d0                            , &amp;
     &amp;                                    HubbleConstant                     =67.36000d0                              &amp;
     &amp;                                   )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctions_"                    >
   <constructor>
    cosmologyFunctionsMatterLambda           (                                                                            &amp;
     &amp;                                    cosmologyParameters_               =cosmologyParameters_                    &amp;
     &amp;                                   )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="linearGrowth_"                          >
   <constructor>
    linearGrowthCollisionlessMatter          (                                                                            &amp;
     &amp;                                    cosmologyParameters_               =cosmologyParameters_                  , &amp;
     &amp;                                    cosmologyFunctions_                =cosmologyFunctions_                     &amp;
     &amp;                                   )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumPrimordialLCDM_"           >
   <constructor>
    powerSpectrumPrimordialPowerLaw          (                                                                            &amp;
     &amp;                                    index_                             =+0.9649d0                             , &amp;
     &amp;                                    running                            =+0.0000d0                             , &amp;
     &amp;                                    runningRunning                     =+0.0000d0                             , &amp;
     &amp;                                    wavenumberReference                =+1.0000d0                             , &amp;
     &amp;                                    runningSmallScalesOnly             =.false.                                 &amp;
     &amp;                                   )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="transferFunctionLCDM_"                  >
   <constructor>
    transferFunctionEisensteinHu1999         (                                                                            &amp;
     &amp;                                    neutrinoNumberEffective            =3.046d0                               , &amp;
     &amp;                                    neutrinoMassSummed                 =0.000d0                               , &amp;
     &amp;                                    darkMatterParticle_                =darkMatterParticleCDM_                , &amp;
     &amp;                                    cosmologyParameters_               =cosmologyParameters_                  , &amp;
     &amp;                                    cosmologyFunctions_                =cosmologyFunctions_                     &amp;
     &amp;                                   )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumPrimordialTransferredLCDM_">
   <constructor>
    powerSpectrumPrimordialTransferredSimple(                                                                             &amp;
     &amp;                                   powerSpectrumPrimordial_           =powerSpectrumPrimordialLCDM_           , &amp;
     &amp;                                   transferFunction_                  =transferFunctionLCDM_                  , &amp;
     &amp;                                   linearGrowth_                      =linearGrowth_                            &amp;
     &amp;                                  )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumWindowFunctionLCDM_"       >
   <constructor>
    powerSpectrumWindowFunctionTopHat       (                                                                             &amp;
     &amp;                                   cosmologyParameters_               =cosmologyParameters_                     &amp;
     &amp;                                  )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologicalMassVarianceLCDM_"          >
   <constructor>
    cosmologicalMassVarianceFilteredPower   (                                                                             &amp;
     &amp;                                   sigma8                             =0.8111d+0                              , &amp;
     &amp;                                   tolerance                          =1.0000d-4                              , &amp;
     &amp;                                   toleranceTopHat                    =1.0000d-4                              , &amp;
     &amp;                                   nonMonotonicIsFatal                =.true.                                 , &amp;
     &amp;                                   monotonicInterpolation             =.false.                                , &amp;
     &amp;                                   truncateAtParticleHorizon          =.false.                                , &amp;
     &amp;                                   cosmologyParameters_               =cosmologyParameters_                   , &amp;
     &amp;                                   cosmologyFunctions_                =cosmologyFunctions_                    , &amp;
     &amp;                                   linearGrowth_                      =linearGrowth_                          , &amp;
     &amp;                                   powerSpectrumPrimordialTransferred_=powerSpectrumPrimordialTransferredLCDM_, &amp;
     &amp;                                   powerSpectrumWindowFunction_       =powerSpectrumWindowFunctionLCDM_         &amp;
     &amp;                                  )
   </constructor>
  </referenceConstruct>
  !!]
  ! Check that converting from mass to σ(M) and back to mass gives consistent answers at z=0.
  do iMass=1,massCount
     sigma        (iMass)=cosmologicalMassVarianceLCDM_%rootVariance(mass (iMass),cosmologyFunctions_%cosmicTime(1.0d0               ))
     massFromSigma(iMass)=cosmologicalMassVarianceLCDM_%mass        (sigma(iMass),cosmologyFunctions_%cosmicTime(1.0d0               ))
  end do
  call Assert('M -> σ(M) -> M conversion loop at z=0 ',mass,massFromSigma,relTol=1.0d-2)
  ! Check that converting from mass to σ(M) and back to mass gives consistent answers at z=10.
  do iMass=1,massCount
     sigma        (iMass)=cosmologicalMassVarianceLCDM_%rootVariance(mass (iMass),cosmologyFunctions_%cosmicTime(1.0d0/(1.0d0+10.0d0)))
     massFromSigma(iMass)=cosmologicalMassVarianceLCDM_%mass        (sigma(iMass),cosmologyFunctions_%cosmicTime(1.0d0/(1.0d0+10.0d0)))
  end do
  call Assert('M -> σ(M) -> M conversion loop at z=10',mass,massFromSigma,relTol=1.0d-2)
  ! Compute the mass corresponding to 8Mpc/h.
  radius8=8.0d0   /cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)
  mass8  =4.0d0*Pi*cosmologyParameters_%densityCritical()*cosmologyParameters_%OmegaMatter()*radius8**3/3.0d0
  sigma8=cosmologicalMassVarianceLCDM_%rootVariance(mass8,cosmologyFunctions_%cosmicTime(1.0d0))
  call Assert('σ₈ equals specified value',sigma8,cosmologicalMassVarianceLCDM_%sigma8(),relTol=1.0d-3)
  ! Check normalization of σ(M) when using a non-top-hat filter. Here, we use a simple power-law power spectrum with index
  ! n=-1, and a sharp k-space filter, and normalize to σ₈=1. For these, the normalization of σ(M₈) can be computed analytically to
  ! be (π√2/3)^{1/3}.
  powerSpectrumPrimordialPowerLaw_         =powerSpectrumPrimordialPowerLaw         (                                                                               &
       &                                                                             index_                             =-1.0d0                                   , &
       &                                                                             running                            =+0.0d0                                   , &
       &                                                                             runningRunning                     =+0.0d0                                   , &
       &                                                                             wavenumberReference                =+1.0d0                                   , &
       &                                                                             runningSmallScalesOnly             =.false.                                    &
       &                                                                            )
  transferFunctionIdentity_                =transferFunctionIdentity                (                                                                               &
       &                                                                             cosmologyParameters_               =cosmologyParameters_                     , &
       &                                                                             time                               =13.8d0                                     &
       &                                                                            )
  powerSpectrumPrimordialTransferredSimple_=powerSpectrumPrimordialTransferredSimple(                                                                               &
       &                                                                             powerSpectrumPrimordial_           =powerSpectrumPrimordialPowerLaw_         , &
       &                                                                             transferFunction_                  =transferFunctionIdentity_                , &
       &                                                                             linearGrowth_                      =linearGrowth_                              &
       &                                                                            )
  powerSpectrumWindowFunctionSharpKSpace_  =powerSpectrumWindowFunctionSharpKSpace  (                                                                               &
       &                                                                             cosmologyParameters_               =cosmologyParameters_                     , &
       &                                                                             normalization                      =0.0d0                                      &
       &                                                                            )
  cosmologicalMassVarianceFilteredPower_   =cosmologicalMassVarianceFilteredPower   (                                                                               &
       &                                                                             sigma8                             =1.0d+0                                   , &
       &                                                                             tolerance                          =1.0d-4                                   , &
       &                                                                             toleranceTopHat                    =1.0d-4                                   , &
       &                                                                             nonMonotonicIsFatal                =.true.                                   , &
       &                                                                             monotonicInterpolation             =.false.                                  , &
       &                                                                             truncateAtParticleHorizon          =.false.                                  , &
       &                                                                             cosmologyParameters_               =cosmologyParameters_                     , &
       &                                                                             cosmologyFunctions_                =cosmologyFunctions_                      , &
       &                                                                             linearGrowth_                      =linearGrowth_                            , &
       &                                                                             powerSpectrumPrimordialTransferred_=powerSpectrumPrimordialTransferredSimple_, &
       &                                                                             powerSpectrumWindowFunction_       =powerSpectrumWindowFunctionSharpKSpace_    &
       &                                                                            )
  sigma8=cosmologicalMassVarianceFilteredPower_%rootVariance(mass8,cosmologyFunctions_%cosmicTime(1.0d0))
  call Assert('σ(M₈) amplitude with sharp k-space window function',sigma8,(sqrt(2.0d0)*Pi/3.0d0)**(1.0d0/3.0d0),relTol=1.0d-4)
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_Sigma
