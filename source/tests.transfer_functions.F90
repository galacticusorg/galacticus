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
Contains a program that tests transfer function calculations.
!!}

program Tests_Transfer_Functions
  !!{
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
  implicit none
  type            (cosmologyParametersSimple               )                             :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda          )                             :: cosmologyFunctions_
  type            (transferFunctionEisensteinHu1999        )                             :: transferFunctionEisensteinHu1999_               , transferFunctionEisensteinHu1999Massless_
  type            (transferFunctionEisensteinHu1998        )                             :: transferFunctionEisensteinHu1998_
  type            (transferFunctionCAMB                    )                             :: transferFunctionCAMB_
  type            (darkMatterParticleCDM                   )                             :: darkMatterParticle_
  type            (linearGrowthCollisionlessMatter         )                             :: linearGrowthCollisionlessMatter_
  type            (powerSpectrumPrimordialTransferredSimple)                             :: powerSpectrumPrimordialTransferredSimple_
  type            (powerSpectrumPrimordialPowerLaw         )                             :: powerSpectrumPrimordialPowerLaw_
  double precision                                          , parameter                  :: stepLogarithmic                          =1.0d-3
  integer                                                   , parameter                  :: wavenumberCount                          =1000
  double precision                                          , parameter                  :: wavenumberMinimum                        =1.0d-3
  double precision                                          , parameter                  :: wavenumberMaximum                        =1.0d+2
  double precision                                          , dimension(wavenumberCount) :: transferFunctionLogarithmicDerivativeEH99       , transferFunctionLogarithmicDerivativeFiniteDifferenceEH99, &
       &                                                                                    transferFunctionLogarithmicDerivativeEH98       , transferFunctionLogarithmicDerivativeFiniteDifferenceEH98, &
       &                                                                                    powerSpectrumLogarithmicDerivativeEH98          , powerSpectrumLogarithmicDerivativeFiniteDifferenceEH98   , &
       &                                                                                    transferFunctionValueEisensteinHu1999           , transferFunctionValueCAMB                                , &
       &                                                                                    wavenumbers
  double precision                                          , dimension(              2) :: wavenumber                                      , powerSpectrumValueEH98_                                  , &
       &                                                                                    transferFunctionValueEH98_                      , transferFunctionValueEH99_
  double precision :: timeNow
  integer                                                                                :: i                                               , j

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
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_Transfer_Functions
