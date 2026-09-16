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

!!{RST
Contains a program that tests transfer function calculations.
!!}

program Tests_Transfer_Functions
  !!{RST
  Tests transfer function calculations.
  !!}
  use :: Cosmology_Functions                 , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters                , only : cosmologyParametersSimple
  use :: Dark_Matter_Particles               , only : darkMatterParticleCDM
  use :: Display                             , only : displayVerbositySet                     , verbosityLevelStandard
  use :: Transfer_Functions                  , only : transferFunctionEisensteinHu1999        , transferFunctionEisensteinHu1998, transferFunctionCAMB, transferFunctionTypeTotal
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredSimple
  use :: Power_Spectra_Primordial            , only : powerSpectrumPrimordialPowerLaw
  use :: Linear_Growth                       , only : linearGrowthCollisionlessMatter
  use :: Unit_Tests                          , only : Assert                                  , Unit_Tests_Begin_Group          , Unit_Tests_End_Group, Unit_Tests_Finish
  use :: Dark_Matter_Particles               , only : darkMatterParticleWDMThermal
  use :: Transfer_Functions                  , only : transferFunctionBode2001                , enumerationScaleCutOffModelEncode
  use :: Numerical_Constants_Math            , only : Pi
  use :: Cosmology_Parameters                , only : hubbleUnitsLittleH
  use :: Transfer_Functions                  , only : transferFunctionAccelerator             , transferFunctionEnvelope
  implicit none
  type            (cosmologyParametersSimple               )                                             :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda          )                                             :: cosmologyFunctions_
  type            (transferFunctionEisensteinHu1999        )                                             :: transferFunctionEisensteinHu1999_                , transferFunctionEisensteinHu1999Massless_
  type            (transferFunctionEisensteinHu1998        )                                             :: transferFunctionEisensteinHu1998_
  type            (transferFunctionCAMB                    )                                             :: transferFunctionCAMB_
  type            (darkMatterParticleCDM                   )                                             :: darkMatterParticle_
  type            (linearGrowthCollisionlessMatter         )                                             :: linearGrowthCollisionlessMatter_
  type            (powerSpectrumPrimordialTransferredSimple)                                             :: powerSpectrumPrimordialTransferredSimple_
  type            (powerSpectrumPrimordialPowerLaw         )                                             :: powerSpectrumPrimordialPowerLaw_
  double precision                                          , parameter                                  :: stepLogarithmic                          =1.0d-3
  integer                                                   , parameter                                  :: wavenumberCount                          =1000
  double precision                                          , parameter                                  :: wavenumberMinimum                        =1.0d-3
  double precision                                          , parameter                                  :: wavenumberMaximum                        =1.0d+2
  double precision                                                     , dimension(wavenumberCount     ) :: transferFunctionLogarithmicDerivativeEH99        , transferFunctionLogarithmicDerivativeFiniteDifferenceEH99, &
       &                                                                                                    transferFunctionLogarithmicDerivativeEH98        , transferFunctionLogarithmicDerivativeFiniteDifferenceEH98, &
       &                                                                                                    powerSpectrumLogarithmicDerivativeEH98           , powerSpectrumLogarithmicDerivativeFiniteDifferenceEH98   , &
       &                                                                                                    transferFunctionValueEisensteinHu1999            , transferFunctionValueCAMB                                , &
       &                                                                                                    wavenumbers
  double precision                                                     , dimension(                   2) :: wavenumber                                       , powerSpectrumValueEH98_                                  , &
       &                                                                                                    transferFunctionValueEH98_                       , transferFunctionValueEH99_
  double precision                                                                                       :: timeNow
  integer                                                                                                :: i                                                , j
  ! Objects and workspace used to check that the tabulations built by the accelerator and envelope transfer functions do not
  ! depend on the order in which wavenumbers are requested of them.
  type            (transferFunctionAccelerator             )                                             :: transferFunctionAcceleratorAscending_            , transferFunctionAcceleratorDescending_
  type            (transferFunctionEnvelope                )                                             :: transferFunctionEnvelopeAscending_               , transferFunctionEnvelopeDescending_
  integer                                                   , parameter                                  :: wavenumberOrderCount                     =5
  double precision                                          , parameter, dimension(wavenumberOrderCount) :: wavenumbersOrder                         =[1.0d-2,1.0d-1,1.0d0,1.0d1,1.0d2]
  double precision                                                     , dimension(wavenumberOrderCount) :: acceleratorAscending                             , acceleratorDescending                                   , &
       &                                                                                                    envelopeAscending                                , envelopeDescending
  ! Objects and workspace for the Bode et al. (2001) warm dark matter modifier.
  type            (darkMatterParticleWDMThermal            )                                             :: darkMatterParticleWDM_                           , darkMatterParticleWDMHeavy_
  type            (transferFunctionBode2001                )                                             :: transferFunctionBode2001_                        , transferFunctionBode2001Heavy_
  ! Parameters of the modifier. These are the defaults of the class.
  double precision                                          , parameter                                  :: epsilonBode                              =0.359d0, etaBode                                  =3.81d0, &
       &                                                                                                    nuBode                                   =1.10d0
  ! Warm dark matter particle masses [keV] and effective degrees of freedom. The reference values against which the cut-off
  ! scale is scaled are those of the class.
  double precision                                          , parameter                                  :: massWDM                                  =3.0d0  , massWDMHeavy                             =6.0d0 , &
       &                                                                                                    degreesOfFreedomWDM                      =1.5d0  , degreesOfFreedomReferenceBode            =1.5d0 , &
       &                                                                                                    massReferenceBode                        =1.0d0
  integer                                                   , parameter                                  :: wavenumberCountBode                      =6
  double precision                                          , parameter, dimension(wavenumberCountBode)  :: wavenumbersBode                          =[1.0d-2,1.0d-1,1.0d0,5.0d0,1.0d1,5.0d1]
  double precision                                                     , dimension(wavenumberCountBode)  :: suppressionClass                                , suppressionExpected
  double precision                                                                                       :: scaleCutOffExpected                             , wavenumberHalfModeExpected                              , &
       &                                                                                                    massHalfModeExpected                            , densityMatter                                           , &
       &                                                                                                    scaleCutOffExpectedHeavy                        , wavenumberHalfModeExpectedHeavy                         , &
       &                                                                                                    massHalfModeExpectedHeavy

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Transfer functions")
  ! Construct required objects.
  cosmologyParameters_                     =cosmologyParametersSimple               (                                                            &
       &                                                                             OmegaMatter            = 0.300d0                          , &
       &                                                                             OmegaBaryon            = 0.045d0                          , &
       &                                                                             OmegaDarkEnergy        = 0.700d0                          , &
       &                                                                             temperatureCMB         = 2.700d0                          , &
       &                                                                             HubbleConstant         =70.0d0                              &
       &                                                                            )
  cosmologyFunctions_                      =cosmologyFunctionsMatterLambda          (                                                            &
       &                                                                             cosmologyParameters_   =cosmologyParameters_                &
       &                                                                            )
  darkMatterParticle_                      =darkMatterParticleCDM                   (                                                            &
       &                                                                            )
  transferFunctionEisensteinHu1999_        =transferFunctionEisensteinHu1999        (                                                            &
       &                                                                             neutrinoNumberEffective=3.046d0                           , &
       &                                                                             neutrinoMassSummed     =0.060d0                           , &
       &                                                                             darkMatterParticle_    =darkMatterParticle_               , &
       &                                                                             cosmologyParameters_   =cosmologyParameters_              , &
       &                                                                             cosmologyFunctions_    =cosmologyFunctions_                 &
       &                                                                            )
  transferFunctionEisensteinHu1999Massless_=transferFunctionEisensteinHu1999        (                                                            &
       &                                                                             neutrinoNumberEffective=3.046d0                           , &
       &                                                                             neutrinoMassSummed     =0.000d0                           , &
       &                                                                             darkMatterParticle_    =darkMatterParticle_               , &
       &                                                                             cosmologyParameters_   =cosmologyParameters_              , &
       &                                                                             cosmologyFunctions_    =cosmologyFunctions_                 &
       &                                                                            )
  transferFunctionEisensteinHu1998_        =transferFunctionEisensteinHu1998        (                                                            &
       &                                                                             darkMatterParticle_    =darkMatterParticle_               , &
       &                                                                             cosmologyParameters_   =cosmologyParameters_              , &
       &                                                                             cosmologyFunctions_    =cosmologyFunctions_                 &
       &                                                                            )
  powerSpectrumPrimordialPowerLaw_         =powerSpectrumPrimordialPowerLaw         (                                                            &
       &                                                                             index_                  =0.9667d0                         , &
       &                                                                             running                 =0.0000d0                         , &
       &                                                                             runningRunning          =0.0000d0                         , &
       &                                                                             wavenumberReference     =1.0000d0                         , &
       &                                                                             runningSmallScalesOnly  =.false.                            &
       &                                                                            )
  linearGrowthCollisionlessMatter_         =linearGrowthCollisionlessMatter         (                                                            &
       &                                                                             cosmologyParameters_    =cosmologyParameters_             , &
       &                                                                             cosmologyFunctions_     =cosmologyFunctions_                &
       &                                                                            )
  powerSpectrumPrimordialTransferredSimple_=powerSpectrumPrimordialTransferredSimple(                                                            &
       &                                                                             powerSpectrumPrimordial_=powerSpectrumPrimordialPowerLaw_ , &
       &                                                                             transferFunction_       =transferFunctionEisensteinHu1998_, &
       &                                                                             linearGrowth_           =linearGrowthCollisionlessMatter_   &
       &                                                                            )  
  transferFunctionCAMB_                    =transferFunctionCAMB                    (                                                            &
       &                                                                             darkMatterParticle_    =darkMatterParticle_               , &
       &                                                                             cosmologyParameters_   =cosmologyParameters_              , &
       &                                                                             cosmologyFunctions_    =cosmologyFunctions_               , &
       &                                                                             transferFunctionType   =transferFunctionTypeTotal         , &
       &                                                                             redshift               =0.0d0                             , &
       &                                                                             cambCountPerDecade     =0                                   &
       &                                                                            )
  ! Find the present time.
  timeNow=cosmologyFunctions_%cosmicTime(1.0d0)
  ! Iterate over reference wavenumbers.
  do j=1,wavenumberCount
     ! Compute logarithmic derivative of transfer function via finite difference.
     wavenumbers(j)=exp(log(wavenumberMinimum)+log(wavenumberMaximum/wavenumberMinimum)*dble(j-1)/dble(wavenumberCount-1))
     wavenumber (1)=wavenumbers(j)
     wavenumber (2)=wavenumbers(j)*exp(stepLogarithmic)
     do i=1,2
        powerSpectrumValueEH98_   (i)=powerSpectrumPrimordialTransferredSimple_%power(wavenumber(i),timeNow)
        transferFunctionValueEH98_(i)=transferFunctionEisensteinHu1998_        %value(wavenumber(i)        )
        transferFunctionValueEH99_(i)=transferFunctionEisensteinHu1999_        %value(wavenumber(i)        )
     end do
     transferFunctionValueEisensteinHu1999                    (j)=+transferFunctionEisensteinHu1999Massless_%value                (wavenumber(1)        )
     transferFunctionValueCAMB                                (j)=+transferFunctionCAMB_                    %value                (wavenumber(1)        )
     powerSpectrumLogarithmicDerivativeEH98                   (j)=+powerSpectrumPrimordialTransferredSimple_%logarithmicDerivative(wavenumber(1),timeNow)
     transferFunctionLogarithmicDerivativeEH98                (j)=+transferFunctionEisensteinHu1998_        %logarithmicDerivative(wavenumber(1)        )
     transferFunctionLogarithmicDerivativeEH99                (j)=+transferFunctionEisensteinHu1999_        %logarithmicDerivative(wavenumber(1)        )
     transferFunctionLogarithmicDerivativeFiniteDifferenceEH98(j)=+log(transferFunctionValueEH98_(2)/transferFunctionValueEH98_(1)) &
          &                                                       /log(wavenumber                (2)/wavenumber                (1))
     powerSpectrumLogarithmicDerivativeFiniteDifferenceEH98   (j)=+log(powerSpectrumValueEH98_   (2)/powerSpectrumValueEH98_   (1)) &
          &                                                       /log(wavenumber                (2)/wavenumber                (1))
     transferFunctionLogarithmicDerivativeFiniteDifferenceEH99(j)=+log(transferFunctionValueEH99_(2)/transferFunctionValueEH99_(1)) &
          &                                                       /log(wavenumber                (2)/wavenumber                (1))
  end do
  ! Normalize transfer functions to their large-scale values.
  transferFunctionValueEisensteinHu1999=+transferFunctionValueEisensteinHu1999    &
       &                                /transferFunctionValueEisensteinHu1999(1)
  transferFunctionValueCAMB            =+transferFunctionValueCAMB                &
       &                                /transferFunctionValueCAMB            (1)
  ! We expect agreement between Eisenstein & Hu (1999) and CAMB over only a limited range of wavenumbers. Outside of that range force them to be equal to avoid failed assertions.
  where(wavenumbers < 1.0d-1 .or. wavenumbers > 1.0d+0)
     transferFunctionValueEisensteinHu1999=transferFunctionValueCAMB
  end where  
  ! Test assertions.
  call Assert('Eisenstein-Hu 1998 P(k)'                                           ,powerSpectrumLogarithmicDerivativeEH98   ,powerSpectrumLogarithmicDerivativeFiniteDifferenceEH98   ,relTol=2.0d+1*stepLogarithmic,absTol=1.0d-2)
  call Assert('Eisenstein-Hu 1998 T(k)'                                           ,transferFunctionLogarithmicDerivativeEH98,transferFunctionLogarithmicDerivativeFiniteDifferenceEH98,relTol=2.0d+0*stepLogarithmic              )
  call Assert('Eisenstein-Hu 1999 T(k) log-derivative with non-zero neutrino mass',transferFunctionLogarithmicDerivativeEH99,transferFunctionLogarithmicDerivativeFiniteDifferenceEH99,relTol=2.0d+0*stepLogarithmic              )
  call Assert('CAMB vs. Eisenstein-Hu T(k) match'                                 ,transferFunctionValueEisensteinHu1999    ,transferFunctionValueCAMB                                ,relTol=4.0d-2                              )
  ! Check that the tabulations built by the accelerator and envelope transfer functions do not depend on the order in which
  ! wavenumbers are requested of them. Both build their tabulation on an absolute lattice, so the wavenumbers they evaluate
  ! depend only on which lattice points are spanned - two objects asked for the same set of wavenumbers must therefore agree
  ! exactly, whichever order they were asked in. Without pinning the bounds follow the first request, and the two disagree at
  ! the level of the interpolation error.
  transferFunctionAcceleratorAscending_ =transferFunctionAccelerator(transferFunctionEisensteinHu1998_,cosmologyParameters_,10)
  transferFunctionAcceleratorDescending_=transferFunctionAccelerator(transferFunctionEisensteinHu1998_,cosmologyParameters_,10)
  do i=1,wavenumberOrderCount
     acceleratorAscending (i)=transferFunctionAcceleratorAscending_ %value(wavenumbersOrder(i))
  end do
  do i=wavenumberOrderCount,1,-1
     acceleratorDescending(i)=transferFunctionAcceleratorDescending_%value(wavenumbersOrder(i))
  end do
  call Assert('accelerator T(k) is independent of the order in which wavenumbers are requested',acceleratorAscending,acceleratorDescending,absTol=0.0d0)
  transferFunctionEnvelopeAscending_ =transferFunctionEnvelope(100,1.0d-2,1.0d2,.false.,.false.,cosmologyParameters_,transferFunctionEisensteinHu1998_)
  transferFunctionEnvelopeDescending_=transferFunctionEnvelope(100,1.0d-2,1.0d2,.false.,.false.,cosmologyParameters_,transferFunctionEisensteinHu1998_)
  do i=1,wavenumberOrderCount
     envelopeAscending    (i)=transferFunctionEnvelopeAscending_    %value(wavenumbersOrder(i))
  end do
  do i=wavenumberOrderCount,1,-1
     envelopeDescending   (i)=transferFunctionEnvelopeDescending_   %value(wavenumbersOrder(i))
  end do
  call Assert('envelope T(k) is independent of the order in which wavenumbers are requested',envelopeAscending,envelopeDescending,absTol=0.0d0)
  ! The warm dark matter modifier of Bode et al. (2001). The class applies
  !
  !   T(k) -> T(k) [1 + (ε k R_c)^(2 ν)]^(-η/ν),
  !
  ! to a cold dark matter transfer function, with the cut-off scale R_c set here by equation (4) of Barkana et al. (2001) - the
  ! default - carrying the factor 0.932 which moves it to the epoch of matter-radiation equality. Both are written out
  ! independently below, so this checks the modifier, the cut-off scale, and the exponents together.
  call Unit_Tests_Begin_Group("Bode et al. (2001) warm dark matter modifier")
  darkMatterParticleWDM_       =darkMatterParticleWDMThermal(mass=massWDM     ,degreesOfFreedomEffective=degreesOfFreedomWDM,cosmologyParameters_=cosmologyParameters_)
  darkMatterParticleWDMHeavy_  =darkMatterParticleWDMThermal(mass=massWDMHeavy,degreesOfFreedomEffective=degreesOfFreedomWDM,cosmologyParameters_=cosmologyParameters_)
  transferFunctionBode2001_    =transferFunctionBode2001(                                                                                              &
       &                                                 transferFunctionCDM =transferFunctionEisensteinHu1998_                                      , &
       &                                                 scaleCutOffModel    =enumerationScaleCutOffModelEncode('barkana2001',includesPrefix=.false.), &
       &                                                 epsilon             =epsilonBode                                                            , &
       &                                                 eta                 =etaBode                                                                , &
       &                                                 nu                  =nuBode                                                                 , &
       &                                                 time                =timeNow                                                                , &
       &                                                 cosmologyParameters_=cosmologyParameters_                                                   , &
       &                                                 darkMatterParticle_ =darkMatterParticleWDM_                                                 , &
       &                                                 cosmologyFunctions_ =cosmologyFunctions_                                                      &
       &                                                )
  transferFunctionBode2001Heavy_=transferFunctionBode2001(                                                                                             &
       &                                                 transferFunctionCDM =transferFunctionEisensteinHu1998_                                      , &
       &                                                 scaleCutOffModel    =enumerationScaleCutOffModelEncode('barkana2001',includesPrefix=.false.), &
       &                                                 epsilon             =epsilonBode                                                            , &
       &                                                 eta                 =etaBode                                                                , &
       &                                                 nu                  =nuBode                                                                 , &
       &                                                 time                =timeNow                                                                , &
       &                                                 cosmologyParameters_=cosmologyParameters_                                                   , &
       &                                                 darkMatterParticle_ =darkMatterParticleWDMHeavy_                                            , &
       &                                                 cosmologyFunctions_ =cosmologyFunctions_                                                      &
       &                                                )
  ! Equation (4) of Barkana et al. (2001), with the 0.932 prefactor.
  scaleCutOffExpected          =+0.932d0                                                             &
       &                        *0.201d0                                                             &
       &                        *(                                                                   &
       &                          +(                                                                 &
       &                            +cosmologyParameters_%OmegaMatter   (                  )         &
       &                            -cosmologyParameters_%OmegaBaryon   (                  )         &
       &                           )                                                                 &
       &                          *  cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2      &
       &                          /0.15d0                                                            &
       &                         )                                                          **0.15d0 &
       &                        /(degreesOfFreedomWDM/degreesOfFreedomReferenceBode)        **0.29d0 &
       &                        /(massWDM            /massReferenceBode            )        **1.15d0
  scaleCutOffExpectedHeavy     =+scaleCutOffExpected            &
       &                        *(massWDM/massWDMHeavy)**1.15d0
  ! The suppression of the transfer function relative to the cold dark matter one.
  do i=1,wavenumberCountBode
     suppressionClass   (i)=+transferFunctionBode2001_        %value(wavenumbersBode(i)) &
          &                 /transferFunctionEisensteinHu1998_%value(wavenumbersBode(i))
     suppressionExpected(i)=+1.0d0                                                                    &
          &                 /(                                                                        &
          &                   +1.0d0                                                                  &
          &                   +(epsilonBode*wavenumbersBode(i)*scaleCutOffExpected)**(2.0d0  *nuBode) &
          &                  )                                                     **(etaBode/nuBode)
  end do
  call Assert('suppression of T(k)',suppressionClass,suppressionExpected,relTol=1.0d-9)
  ! The half-mode mass. The wavenumber at which the suppression is a factor of two follows in closed form from the modifier, and
  ! the mass from the convention R = lambda/2 = pi/k.
  densityMatter                 =+cosmologyParameters_%OmegaMatter    () &
       &                         *cosmologyParameters_%densityCritical()
  wavenumberHalfModeExpected    =+(                                       &
       &                           +2.0d0**(+nuBode/etaBode)              &
       &                           -1.0d0                                 &
       &                          )      **(+0.5d0 /nuBode )              &
       &                         /epsilonBode                             &
       &                         /scaleCutOffExpected
  massHalfModeExpected          =+4.0d0                                   &
       &                         *Pi                                      &
       &                         /3.0d0                                   &
       &                         *densityMatter                           &
       &                         *(Pi/wavenumberHalfModeExpected)**3
  wavenumberHalfModeExpectedHeavy=+(                                      &
       &                           +2.0d0**(+nuBode/etaBode)              &
       &                           -1.0d0                                 &
       &                          )      **(+0.5d0 /nuBode )              &
       &                         /epsilonBode                             &
       &                         /scaleCutOffExpectedHeavy
  massHalfModeExpectedHeavy     =+4.0d0                                   &
       &                         *Pi                                      &
       &                         /3.0d0                                   &
       &                         *densityMatter                           &
       &                         *(Pi/wavenumberHalfModeExpectedHeavy)**3
  call Assert('half-mode mass'                  ,transferFunctionBode2001_     %halfModeMass(),massHalfModeExpected     ,relTol=1.0d-9)
  call Assert('half-mode mass, heavier particle',transferFunctionBode2001Heavy_%halfModeMass(),massHalfModeExpectedHeavy,relTol=1.0d-9)
  ! The transfer function must indeed be suppressed by precisely a factor of two at the half-mode wavenumber, which ties the mass
  ! above to the definition it is meant to express.
  call Assert('T(k) is suppressed by two at the half-mode wavenumber'                                            &
       &     ,+transferFunctionBode2001_        %value(transferFunctionBode2001_%wavenumberAtSuppression(2.0d0)) &
       &      /transferFunctionEisensteinHu1998_%value(transferFunctionBode2001_%wavenumberAtSuppression(2.0d0)) &
       &     ,0.5d0,relTol=1.0d-9)
  ! Limits: the modifier leaves large scales untouched, suppresses monotonically, and gives a lower half-mode mass for a heavier
  ! particle - the cut-off scale going as m^-1.15, so the mass as m^-3.45.
  call Assert('T(k) is unmodified on large scales',suppressionClass(1),1.0d0,relTol=1.0d-6)
  call Assert('suppression is monotonic'          ,all(suppressionClass(2:wavenumberCountBode) < suppressionClass(1:wavenumberCountBode-1)),.true.)
  call Assert('a heavier particle has a lower half-mode mass',transferFunctionBode2001Heavy_%halfModeMass() < transferFunctionBode2001_%halfModeMass(),.true.)
  call Assert('half-mode mass scales as m^-3.45'  ,transferFunctionBode2001_%halfModeMass()/transferFunctionBode2001Heavy_%halfModeMass(),(massWDMHeavy/massWDM)**(3.0d0*1.15d0),relTol=1.0d-9)
  call Unit_Tests_End_Group()

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_Transfer_Functions
