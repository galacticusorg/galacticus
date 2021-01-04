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

program Tests_Sigma
  !% Tests
  use :: Cosmological_Density_Field          , only : cosmologicalMassVariance                , cosmologicalMassVarianceClass, cosmologicalMassVarianceFilteredPower
  use :: Cosmology_Functions                 , only : cosmologyFunctions                      , cosmologyFunctionsClass
  use :: Cosmology_Parameters                , only : cosmologyParameters                     , cosmologyParametersClass     , hubbleUnitsLittleH
  use :: Galacticus_Display                  , only : Galacticus_Verbosity_Level_Set          , verbosityStandard
  use :: ISO_Varying_String                  , only : varying_string                          , assignment(=)
  use :: Input_Parameters                    , only : inputParameters
  use :: Linear_Growth                       , only : linearGrowth                            , linearGrowthClass
  use :: Numerical_Constants_Math            , only : Pi
  use :: Numerical_Ranges                    , only : Make_Range                              , rangeTypeLogarithmic
  use :: Power_Spectra_Primordial            , only : powerSpectrumPrimordialPowerLaw
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredSimple
  use :: Power_Spectrum_Window_Functions     , only : powerSpectrumWindowFunctionSharpKSpace
  use :: Transfer_Functions                  , only : transferFunctionIdentity
  use :: Unit_Tests                          , only : Assert                                  , Unit_Tests_Begin_Group       , Unit_Tests_End_Group                 , Unit_Tests_Finish
  implicit none
  type            (varying_string                          )                       :: parameterFile
  integer                                                   , parameter            :: massCount                             =10
  double precision                                          , parameter            :: massMaximum                           =1.0d15, massMinimum  =1.0d6
  double precision                                          , dimension(massCount) :: mass                                         , massFromSigma      , &
       &                                                                              sigma
  class           (cosmologyParametersClass                ), pointer              :: cosmologyParameters_
  class           (cosmologyFunctionsClass                 ), pointer              :: cosmologyFunctions_
  class           (linearGrowthClass                       ), pointer              :: linearGrowth_
  class           (cosmologicalMassVarianceClass           ), pointer              :: cosmologicalMassVariance_
  type            (cosmologicalMassVarianceFilteredPower   )                       :: cosmologicalMassVarianceFilteredPower_
  type            (powerSpectrumWindowFunctionSharpKSpace  )                       :: powerSpectrumWindowFunctionSharpKSpace_
  type            (powerSpectrumPrimordialPowerLaw         )                       :: powerSpectrumPrimordialPowerLaw_
  type            (transferFunctionIdentity                )                       :: transferFunctionIdentity_
  type            (powerSpectrumPrimordialTransferredSimple)                       :: powerSpectrumPrimordialTransferredSimple_
  integer                                                                          :: iMass
  double precision                                                                 :: mass8                                        , radius8            , &
       &                                                                              sigma8
  type            (inputParameters                         )                       :: parameters

  ! Set verbosity level.
  call Galacticus_Verbosity_Level_Set(verbosityStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Power spectrum: σ(M)")
  ! Open the parameter file.
  parameterFile='parameters.xml'
  parameters=inputParameters(parameterFile)
  call parameters%markGlobal()
  ! Create an array of masses.
  mass=Make_Range(massMinimum,massMaximum,massCount,rangeType=rangeTypeLogarithmic)
  ! Get required objects.
  cosmologyParameters_      => cosmologyParameters     ()
  cosmologyFunctions_       => cosmologyFunctions      ()
  linearGrowth_             => linearGrowth            ()
  cosmologicalMassVariance_ => cosmologicalMassVariance()
  ! Check that converting from mass to sigma and back to mass gives consistent answers.
  do iMass=1,massCount
     sigma        (iMass)=cosmologicalMassVariance_%rootVariance(mass (iMass),cosmologyFunctions_%cosmicTime(1.0d0))
     massFromSigma(iMass)=cosmologicalMassVariance_%mass        (sigma(iMass),cosmologyFunctions_%cosmicTime(1.0d0))
  end do
  call Assert('M -> σ(M) -> M conversion loop',mass,massFromSigma,relTol=1.0d-2)
  ! Compute the mass corresponding to 8Mpc/h.
  radius8=8.0d0   /cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)
  mass8  =4.0d0*Pi*cosmologyParameters_%densityCritical()*cosmologyParameters_%OmegaMatter()*radius8**3/3.0d0
  sigma8=cosmologicalMassVariance_%rootVariance(mass8,cosmologyFunctions_%cosmicTime(1.0d0))
  call Assert('σ₈ equals specified value',sigma8,cosmologicalMassVariance_%sigma8(),relTol=1.0d-3)
  ! Check normalization of sigma(M) when using a non-top-hat filter. Here, we use a simple power-law power spectrum with index
  ! n=-1, and a sharp k-space filter, and normalize to σ₈=1. For these, the normalization of σ(M₈) can be computed analytically to
  ! be (π√2/3)^{1/3}.
  powerSpectrumPrimordialPowerLaw_         =powerSpectrumPrimordialPowerLaw         (                                                                               &
       &                                                                             index                              =-1.0d0                                   , &
       &                                                                             running                            =+0.0d0                                   , &
       &                                                                             wavenumberReference                =+1.0d0                                     &
       &                                                                            )
  transferFunctionIdentity_                =transferFunctionIdentity                (                                                                               &
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
