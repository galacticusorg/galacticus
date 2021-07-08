!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
  use :: Transfer_Functions   , only : transferFunctionEisensteinHu1999
  use :: Unit_Tests           , only : Assert                          , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  type            (cosmologyParametersSimple       )               :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda  )               :: cosmologyFunctions_
  type            (transferFunctionEisensteinHu1999)               :: transferFunctionEisensteinHu1999_
  type            (darkMatterParticleCDM           )               :: darkMatterParticle_
  double precision                                  , parameter    :: stepLogarithmic                     =1.0d-3
  double precision                                  , dimension(2) :: wavenumber                                  , transferFunctionValue                                , &
       &                                                              wavenumberReference
  double precision                                                 :: transferFunctionLogarithmicDerivative       , transferFunctionLogarithmicDerivativeFiniteDifference
  integer                                                          :: i                                           , j
  character       (len=9                          )                :: label

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Transfer functions")
  ! Construct required objects.
  cosmologyParameters_             =cosmologyParametersSimple       (                                              &
       &                                                             OmegaMatter            = 0.300d0            , &
       &                                                             OmegaBaryon            = 0.045d0            , &
       &                                                             OmegaDarkEnergy        = 0.700d0            , &
       &                                                             temperatureCMB         = 2.700d0            , &
       &                                                             HubbleConstant         =70.0d0                &
       &                                                            )
  cosmologyFunctions_              =cosmologyFunctionsMatterLambda  (                                              &
       &                                                             cosmologyParameters_   =cosmologyParameters_  &
       &                                                            )
  darkMatterParticle_              =darkMatterParticleCDM           (                                              &
       &                                                            )
  transferFunctionEisensteinHu1999_=transferFunctionEisensteinHu1999(                                              &
       &                                                             neutrinoNumberEffective=3.046d0             , &
       &                                                             neutrinoMassSummed     =0.060d0             , &
       &                                                             darkMatterParticle_    =darkMatterParticle_ , &
       &                                                             cosmologyParameters_   =cosmologyParameters_, &
       &                                                             cosmologyFunctions_    =cosmologyFunctions_   &
       &                                                            )
  ! Iterate over reference wavenumbers.
  wavenumberReference=[1.0d-2,1.0d+0]
  do j=1,size(wavenumberReference)
     ! Compute logarithmic derivative of transfer function via finite difference.
     wavenumber(1)=wavenumberReference(j)
     wavenumber(2)=wavenumberReference(j)*exp(stepLogarithmic)
     do i=1,2
        transferFunctionValue(i)=transferFunctionEisensteinHu1999_%value(wavenumber(i))
     end do
     transferFunctionLogarithmicDerivative                =transferFunctionEisensteinHu1999_%logarithmicDerivative(wavenumber(1))
     transferFunctionLogarithmicDerivativeFiniteDifference=+log(transferFunctionValue(2)/transferFunctionValue(1)) &
          &                                                /log(wavenumber           (2)/wavenumber           (1))
     write (label,'(e8.2)') wavenumberReference(j)
     call Assert('Eisenstein-Hu T(k) log-derivative with non-zero neutrino mass at k='//trim(label)//'Mpc',transferFunctionLogarithmicDerivative,transferFunctionLogarithmicDerivativeFiniteDifference,relTol=stepLogarithmic)
  end do
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_Transfer_Functions
