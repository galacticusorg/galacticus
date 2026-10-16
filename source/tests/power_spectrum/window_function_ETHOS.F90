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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a program that tests the :cite:t:`bohr_halo_2021` ETHOS power spectrum window function and its extended form.
!!}

program Tests_Power_Spectrum_Window_Function_ETHOS
  !!{RST
  Tests the :galacticus-class:`powerSpectrumWindowFunctionETHOSExtended` window function against exact results.

  The parameters of the extended window function depend on the logarithmic slope, :math:`n`, of the Gaussian-smoothed power
  spectrum as :math:`x = x_0 x_1^{n-n_0}`, with :math:`n_0=-2.6`. For a power-law power spectrum, :math:`P(k) \propto k^n`, the
  smoothed power spectrum is :math:`\tilde{P}(k) \propto k^n \exp(n^2\sigma^2/2)`, so its logarithmic slope is exactly :math:`n`
  at all wavenumbers and :math:`c_\mathrm{W}` and :math:`\beta` are known exactly. This is tested for two power-law indices. With
  :math:`c_\mathrm{W,1}=\beta_1=1` and :math:`x_\mathrm{min}=0` the extended window function must reduce exactly to the original
  ETHOS window function. Finally, the window function must be unity for :math:`kR/c_\mathrm{W} \le x_\mathrm{min}`, must equal
  :math:`1/2` at :math:`kR/c_\mathrm{W} = x_\mathrm{min}+1`, and must tend to unity as :math:`k \rightarrow 0`.
  !!}
  use :: Cosmology_Functions                 , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters                , only : cosmologyParametersSimple
  use :: Display                             , only : displayVerbositySet                     , verbosityLevelStandard
  use :: Linear_Growth                       , only : linearGrowthCollisionlessMatter
  use :: Numerical_Constants_Math            , only : Pi
  use :: Power_Spectra_Primordial            , only : powerSpectrumPrimordialPowerLaw
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredSimple
  use :: Power_Spectrum_Window_Functions     , only : powerSpectrumWindowFunctionETHOS        , powerSpectrumWindowFunctionETHOSExtended
  use :: Transfer_Functions                  , only : transferFunctionIdentity
  use :: Unit_Tests                          , only : Assert                                  , Unit_Tests_Begin_Group                  , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  type            (cosmologyParametersSimple               )               :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda          )               :: cosmologyFunctions_
  type            (linearGrowthCollisionlessMatter         )               :: linearGrowth_
  type            (transferFunctionIdentity                )               :: transferFunction_
  type            (powerSpectrumPrimordialPowerLaw         )               :: powerSpectrumPrimordialShallow_       , powerSpectrumPrimordialSteep_
  type            (powerSpectrumPrimordialTransferredSimple)               :: powerSpectrumTransferredShallow_      , powerSpectrumTransferredSteep_
  type            (powerSpectrumWindowFunctionETHOS        )               :: windowFunctionETHOS_
  type            (powerSpectrumWindowFunctionETHOSExtended)               :: windowFunctionUnity_                  , windowFunctionShallow_              , &
       &                                                                      windowFunctionSteep_
  ! Parameters of the extended window function, similar to those used in the reference models.
  double precision                                          , parameter    :: cW0                            =3.50d0, cW1                          =0.80d0, &
       &                                                                      beta0                          =4.50d0, beta1                        =1.40d0, &
       &                                                                      wavenumberScaledMinimum        =0.20d0, powerSpectrumSmoothingWidth  =0.68d0
  ! Power-law indices of the two power spectra, and the zero-point in slope used by the extended window function.
  double precision                                          , parameter    :: indexShallow                   =-2.0d0, indexSteep                   =-2.9d0, &
       &                                                                      indexReference                 =-2.6d0
  ! Smoothing mass [M☉].
  double precision                                          , parameter    :: massSmoothing                  =1.0d10
  ! Wavenumbers [Mpc⁻¹] at which to test.
  double precision                                          , dimension(4) :: wavenumber                     =[1.0d-1,1.0d0,1.0d1,1.0d2]
  double precision                                          , dimension(4) :: cWShallow                             , betaShallow                         , &
       &                                                                      cWSteep                               , betaSteep                           , &
       &                                                                      valueETHOS                            , valueUnity
  double precision                                                         :: time                                  , radius                              , &
       &                                                                      cWExpected                            , wavenumberScaled
  integer                                                                  :: i

  call displayVerbositySet   (verbosityLevelStandard                                           )
  call Unit_Tests_Begin_Group("Power spectrum window functions: ETHOS and extended ETHOS")
  !![
  <referenceConstruct object="cosmologyParameters_"           >
   <constructor>
    cosmologyParametersSimple               (                                                                                 &amp;
     &amp;                                   OmegaMatter                 = 0.28120d0                                        , &amp;
     &amp;                                   OmegaBaryon                 = 0.04611d0                                        , &amp;
     &amp;                                   OmegaDarkEnergy             = 0.71880d0                                        , &amp;
     &amp;                                   temperatureCMB              = 2.72548d0                                        , &amp;
     &amp;                                   HubbleConstant              =69.70000d0                                          &amp;
     &amp;                                  )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctions_"            >
   <constructor>
    cosmologyFunctionsMatterLambda          (                                                                                 &amp;
     &amp;                                   cosmologyParameters_        =cosmologyParameters_                                &amp;
     &amp;                                  )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="linearGrowth_"                  >
   <constructor>
    linearGrowthCollisionlessMatter         (                                                                                 &amp;
     &amp;                                   cosmologyParameters_        =cosmologyParameters_                              , &amp;
     &amp;                                   cosmologyFunctions_         =cosmologyFunctions_                                 &amp;
     &amp;                                  )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="transferFunction_"              >
   <constructor>
    transferFunctionIdentity                (                                                                                 &amp;
     &amp;                                   cosmologyParameters_        =cosmologyParameters_                              , &amp;
     &amp;                                   time                        =cosmologyFunctions_%cosmicTime(expansionFactor=1.0d0) &amp;
     &amp;                                  )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumPrimordialShallow_">
   <constructor>
    powerSpectrumPrimordialPowerLaw         (                                                                                 &amp;
     &amp;                                   index_                      =indexShallow                                      , &amp;
     &amp;                                   running                     =0.0d0                                             , &amp;
     &amp;                                   runningRunning              =0.0d0                                             , &amp;
     &amp;                                   wavenumberReference         =1.0d0                                             , &amp;
     &amp;                                   runningSmallScalesOnly      =.false.                                             &amp;
     &amp;                                  )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumPrimordialSteep_"  >
   <constructor>
    powerSpectrumPrimordialPowerLaw         (                                                                                 &amp;
     &amp;                                   index_                      =indexSteep                                        , &amp;
     &amp;                                   running                     =0.0d0                                             , &amp;
     &amp;                                   runningRunning              =0.0d0                                             , &amp;
     &amp;                                   wavenumberReference         =1.0d0                                             , &amp;
     &amp;                                   runningSmallScalesOnly      =.false.                                             &amp;
     &amp;                                  )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumTransferredShallow_">
   <constructor>
    powerSpectrumPrimordialTransferredSimple(                                                                                 &amp;
     &amp;                                   powerSpectrumPrimordial_    =powerSpectrumPrimordialShallow_                   , &amp;
     &amp;                                   transferFunction_           =transferFunction_                                 , &amp;
     &amp;                                   linearGrowth_               =linearGrowth_                                       &amp;
     &amp;                                  )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="powerSpectrumTransferredSteep_" >
   <constructor>
    powerSpectrumPrimordialTransferredSimple(                                                                                 &amp;
     &amp;                                   powerSpectrumPrimordial_    =powerSpectrumPrimordialSteep_                     , &amp;
     &amp;                                   transferFunction_           =transferFunction_                                 , &amp;
     &amp;                                   linearGrowth_               =linearGrowth_                                       &amp;
     &amp;                                  )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="windowFunctionETHOS_"           >
   <constructor>
    powerSpectrumWindowFunctionETHOS        (                                                                                 &amp;
     &amp;                                   cW_                         =cW0                                               , &amp;
     &amp;                                   beta_                       =beta0                                             , &amp;
     &amp;                                   cosmologyParameters_        =cosmologyParameters_                                &amp;
     &amp;                                  )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="windowFunctionUnity_"           >
   <constructor>
    powerSpectrumWindowFunctionETHOSExtended(                                                                                 &amp;
     &amp;                                   cW0                         =cW0                                               , &amp;
     &amp;                                   cW1                         =1.0d0                                             , &amp;
     &amp;                                   beta0                       =beta0                                             , &amp;
     &amp;                                   beta1                       =1.0d0                                             , &amp;
     &amp;                                   wavenumberScaledMinimum_    =0.0d0                                             , &amp;
     &amp;                                   powerSpectrumSmoothingWidth =powerSpectrumSmoothingWidth                       , &amp;
     &amp;                                   cosmologyParameters_        =cosmologyParameters_                              , &amp;
     &amp;                                   powerSpectrumPrimordialTransferred_=powerSpectrumTransferredShallow_             &amp;
     &amp;                                  )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="windowFunctionShallow_"         >
   <constructor>
    powerSpectrumWindowFunctionETHOSExtended(                                                                                 &amp;
     &amp;                                   cW0                         =cW0                                               , &amp;
     &amp;                                   cW1                         =cW1                                               , &amp;
     &amp;                                   beta0                       =beta0                                             , &amp;
     &amp;                                   beta1                       =beta1                                             , &amp;
     &amp;                                   wavenumberScaledMinimum_    =wavenumberScaledMinimum                           , &amp;
     &amp;                                   powerSpectrumSmoothingWidth =powerSpectrumSmoothingWidth                       , &amp;
     &amp;                                   cosmologyParameters_        =cosmologyParameters_                              , &amp;
     &amp;                                   powerSpectrumPrimordialTransferred_=powerSpectrumTransferredShallow_             &amp;
     &amp;                                  )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="windowFunctionSteep_"           >
   <constructor>
    powerSpectrumWindowFunctionETHOSExtended(                                                                                 &amp;
     &amp;                                   cW0                         =cW0                                               , &amp;
     &amp;                                   cW1                         =cW1                                               , &amp;
     &amp;                                   beta0                       =beta0                                             , &amp;
     &amp;                                   beta1                       =beta1                                             , &amp;
     &amp;                                   wavenumberScaledMinimum_    =0.0d0                                             , &amp;
     &amp;                                   powerSpectrumSmoothingWidth =powerSpectrumSmoothingWidth                       , &amp;
     &amp;                                   cosmologyParameters_        =cosmologyParameters_                              , &amp;
     &amp;                                   powerSpectrumPrimordialTransferred_=powerSpectrumTransferredSteep_               &amp;
     &amp;                                  )
   </constructor>
  </referenceConstruct>
  !!]
  time  =cosmologyFunctions_%cosmicTime(1.0d0)
  radius=+(                                           &
       &   +3.0d0                                     &
       &   /4.0d0                                     &
       &   /Pi                                        &
       &   *massSmoothing                             &
       &   /cosmologyParameters_%OmegaMatter    ()    &
       &   /cosmologyParameters_%densityCritical()    &
       &  )**(1.0d0/3.0d0)
  ! For power-law power spectra the parameters are known exactly at all wavenumbers. The smoothed slope is found by numerical
  ! integration to a relative tolerance of 10⁻⁶, so we allow an error ten times larger than this in the parameters.
  do i=1,size(wavenumber)
     cWShallow  (i)=windowFunctionShallow_%cW  (wavenumber(i),time)
     betaShallow(i)=windowFunctionShallow_%beta(wavenumber(i),time)
     cWSteep    (i)=windowFunctionSteep_  %cW  (wavenumber(i),time)
     betaSteep  (i)=windowFunctionSteep_  %beta(wavenumber(i),time)
  end do
  call Assert('c_W for P(k) ∝ k^-2.0',cWShallow  ,spread(cW0  *cW1  **(indexShallow-indexReference),1,size(wavenumber)),relTol=1.0d-5)
  call Assert('β   for P(k) ∝ k^-2.0',betaShallow,spread(beta0*beta1**(indexShallow-indexReference),1,size(wavenumber)),relTol=1.0d-5)
  call Assert('c_W for P(k) ∝ k^-2.9',cWSteep    ,spread(cW0  *cW1  **(indexSteep  -indexReference),1,size(wavenumber)),relTol=1.0d-5)
  call Assert('β   for P(k) ∝ k^-2.9',betaSteep  ,spread(beta0*beta1**(indexSteep  -indexReference),1,size(wavenumber)),relTol=1.0d-5)
  ! With c_W1=β₁=1 and x_min=0 the extended window function must be identical to the original ETHOS window function.
  do i=1,size(wavenumber)
     valueETHOS(i)=windowFunctionETHOS_%value(wavenumber(i),massSmoothing,time)
     valueUnity(i)=windowFunctionUnity_%value(wavenumber(i),massSmoothing,time)
  end do
  call Assert('reduces to ETHOS for c_W1=β₁=1',valueUnity,valueETHOS,relTol=1.0d-12)
  ! The window function must be unity below x_min, and one half at x=x_min+1.
  cWExpected      =cW0*cW1**(indexShallow-indexReference)
  wavenumberScaled=0.5d0*wavenumberScaledMinimum
  call Assert('W=1 for kR/c_W < x_min'  ,windowFunctionShallow_%value(wavenumberScaled*cWExpected/radius,massSmoothing,time),1.0d0,relTol=1.0d-12)
  wavenumberScaled=1.0d0+wavenumberScaledMinimum
  call Assert('W=½ at kR/c_W = x_min+1',windowFunctionShallow_%value(wavenumberScaled*cWExpected/radius,massSmoothing,time),0.5d0,relTol=1.0d-5 )
  ! The window function must tend to unity as k→0.
  call Assert('W→1 as k→0'              ,windowFunctionSteep_  %value(1.0d-8                              ,massSmoothing,time),1.0d0,relTol=1.0d-12)
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_Power_Spectrum_Window_Function_ETHOS
