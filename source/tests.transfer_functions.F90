!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
  use :: Cosmology_Functions  , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters , only : cosmologyParametersSimple
  use :: Dark_Matter_Particles, only : darkMatterParticleCDM
  use :: Display              , only : displayVerbositySet             , verbosityLevelStandard
  use :: Transfer_Functions   , only : transferFunctionEisensteinHu1999, transferFunctionCAMB
  use :: Unit_Tests           , only : Assert                          , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  type            (cosmologyParametersSimple       )                             :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda  )                             :: cosmologyFunctions_
  type            (transferFunctionEisensteinHu1999)                             :: transferFunctionEisensteinHu1999_           , transferFunctionEisensteinHu1999Massless_
  type            (transferFunctionCAMB            )                             :: transferFunctionCAMB_
  type            (darkMatterParticleCDM           )                             :: darkMatterParticle_
  double precision                                  , parameter                  :: stepLogarithmic                      =1.0d-3
  integer                                           , parameter                  :: wavenumberCount                      =10
  double precision                                  , parameter                  :: wavenumberMinimum                    =1.0d-2
  double precision                                  , parameter                  :: wavenumberMaximum                    =1.0d+1
  double precision                                  , dimension(wavenumberCount) :: transferFunctionLogarithmicDerivative       , transferFunctionLogarithmicDerivativeFiniteDifference, &
       &                                                                            transferFunctionValueEisensteinHu1999       , transferFunctionValueCAMB
  double precision                                  , dimension(              2) :: wavenumber                                  , transferFunctionValue
  integer                                                                        :: i                                           , j

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Transfer functions")
  ! Construct required objects.
  cosmologyParameters_                     =cosmologyParametersSimple       (                                              &
       &                                                                     OmegaMatter            = 0.300d0            , &
       &                                                                     OmegaBaryon            = 0.045d0            , &
       &                                                                     OmegaDarkEnergy        = 0.700d0            , &
       &                                                                     temperatureCMB         = 2.700d0            , &
       &                                                                     HubbleConstant         =70.0d0                &
       &                                                                    )
  cosmologyFunctions_                      =cosmologyFunctionsMatterLambda  (                                              &
       &                                                                     cosmologyParameters_   =cosmologyParameters_  &
       &                                                                    )
  darkMatterParticle_                      =darkMatterParticleCDM           (                                              &
       &                                                                    )
  transferFunctionEisensteinHu1999_        =transferFunctionEisensteinHu1999(                                              &
       &                                                                     neutrinoNumberEffective=3.046d0             , &
       &                                                                     neutrinoMassSummed     =0.060d0             , &
       &                                                                     darkMatterParticle_    =darkMatterParticle_ , &
       &                                                                     cosmologyParameters_   =cosmologyParameters_, &
       &                                                                     cosmologyFunctions_    =cosmologyFunctions_   &
       &                                                                    )
  transferFunctionEisensteinHu1999Massless_=transferFunctionEisensteinHu1999(                                              &
       &                                                                     neutrinoNumberEffective=3.046d0             , &
       &                                                                     neutrinoMassSummed     =0.000d0             , &
       &                                                                     darkMatterParticle_    =darkMatterParticle_ , &
       &                                                                     cosmologyParameters_   =cosmologyParameters_, &
       &                                                                     cosmologyFunctions_    =cosmologyFunctions_   &
       &                                                                    )
  transferFunctionCAMB_                    =transferFunctionCAMB            (                                              &
       &                                                                     darkMatterParticle_    =darkMatterParticle_ , &
       &                                                                     cosmologyParameters_   =cosmologyParameters_, &
       &                                                                     cosmologyFunctions_    =cosmologyFunctions_ , &
       &                                                                     redshift               =0.0d0               , &
       &                                                                     cambCountPerDecade     =0                     &
       &                                                                    )
  ! Iterate over reference wavenumbers.
  do j=1,wavenumberCount
     ! Compute logarithmic derivative of transfer function via finite difference.
     wavenumber(1)=exp(log(wavenumberMinimum)+log(wavenumberMaximum/wavenumberMinimum)*dble(j-1)/dble(wavenumberCount-1))
     wavenumber(2)=wavenumber(1)*exp(stepLogarithmic)
     do i=1,2
        transferFunctionValue(i)=transferFunctionEisensteinHu1999_%value(wavenumber(i))
     end do
     transferFunctionValueEisensteinHu1999                (j)=+transferFunctionEisensteinHu1999Massless_%value                (wavenumber(1))
     transferFunctionValueCAMB                            (j)=+transferFunctionCAMB_                    %value                (wavenumber(1))
     transferFunctionLogarithmicDerivative                (j)=+transferFunctionEisensteinHu1999_        %logarithmicDerivative(wavenumber(1))
     transferFunctionLogarithmicDerivativeFiniteDifference(j)=+log(transferFunctionValue(2)/transferFunctionValue(1)) &
          &                                                   /log(wavenumber           (2)/wavenumber           (1))
  end do
  ! Normalize transfer functions to their large-scale values.
  transferFunctionValueEisensteinHu1999=+transferFunctionValueEisensteinHu1999    &
       &                                /transferFunctionValueEisensteinHu1999(1)
  transferFunctionValueCAMB            =+transferFunctionValueCAMB                &
       &                                /transferFunctionValueCAMB            (1)
  ! Test assertions.
  call Assert('Eisenstein-Hu T(k) log-derivative with non-zero neutrino mass',transferFunctionLogarithmicDerivative,transferFunctionLogarithmicDerivativeFiniteDifference,relTol=stepLogarithmic)
  call Assert('CAMB vs. Eisenstein-Hu T(k) match'                            ,transferFunctionValueEisensteinHu1999,transferFunctionValueCAMB                            ,relTol=2.0d-2         )
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_Transfer_Functions
