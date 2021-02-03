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

!% Contains a program that tests power spectrum calculations.

program Tests_Power_Spectrum
  !% Tests power spectrum calculations.
  use :: Cosmological_Density_Field          , only : cosmologicalMassVariance                , cosmologicalMassVarianceClass
  use :: Cosmology_Functions                 , only : cosmologyFunctions                      , cosmologyFunctionsClass
  use :: Cosmology_Parameters                , only : cosmologyParameters                     , cosmologyParametersClass     , hubbleUnitsLittleH
  use :: Dark_Matter_Particles               , only : darkMatterParticleCDM
  use :: Display                             , only : displayVerbositySet                     , verbosityLevelStandard
  use :: ISO_Varying_String                  , only : assignment(=)                           , varying_string
  use :: Input_Parameters                    , only : inputParameters
  use :: Linear_Growth                       , only : linearGrowth                            , linearGrowthClass
  use :: Numerical_Constants_Math            , only : Pi
  use :: Power_Spectra                       , only : powerSpectrum                           , powerSpectrumClass
  use :: Power_Spectra_Primordial            , only : powerSpectrumPrimordialPowerLaw
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredSimple
  use :: Transfer_Functions                  , only : transferFunctionEisensteinHu1999
  use :: Unit_Tests                          , only : Assert                                  , Unit_Tests_Begin_Group       , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  type            (inputParameters                         ), target       :: parameters
  double precision                                          , parameter    :: radiusNormalization                      =8.0d0 ! Radius for σ(M) normalization in Mpc/h.
  type            (darkMatterParticleCDM                   ), target       :: darkMatterParticleCDM_
  class           (cosmologyParametersClass                ), pointer      :: cosmologyParameters_
  class           (cosmologyFunctionsClass                 ), pointer      :: cosmologyFunctions_
  type            (powerSpectrumPrimordialPowerLaw         ), target       :: powerSpectrumPrimordialPowerLaw_
  type            (powerSpectrumPrimordialTransferredSimple), target       :: powerSpectrumPrimordialTransferredSimple_
  type            (transferFunctionEisensteinHu1999        ), target       :: transferFunctionEisensteinHu1999_
  class           (cosmologicalMassVarianceClass           ), pointer      :: cosmologicalMassVariance_
  class           (powerSpectrumClass                      ), pointer      :: powerSpectrum_
  class           (linearGrowthClass                       ), pointer      :: linearGrowth_
  double precision                                          , dimension(5) :: powerComputed                                                   , &
       &                                                                      transferComputed                                                , &
       &                                                                      powerTransferredComputed                                        , &
       &                                                                      wavenumber                               =[                       &
       &                                                                                                                 0.1000000000000000d+0, &
       &                                                                                                                 0.3000000000000000d+0, &
       &                                                                                                                 1.0000000000000000d+0, &
       &                                                                                                                 3.0000000000000000d+0, &
       &                                                                                                                 1.0000000000000000d+1  &
       &                                                                                                                ]                     , &
       &                                                                      transferExpected                         =[                       &
       &                                                                                                                 7.7501624542585204d-2, &
       &                                                                                                                 1.5825809489746621d-2, &
       &                                                                                                                 2.2732135149658283d-3, &
       &                                                                                                                 3.4629417179195417d-4, &
       &                                                                                                                 4.0676349951660953d-5  &
       &                                                                                                                ]
  integer                                                                  :: i
  type            (varying_string                          )               :: parameterFile
  double precision                                                         :: mass                                                            , &
       &                                                                      ratio                                                           , &
       &                                                                      sigma

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Read parameters.
  parameterFile='testSuite/parameters/powerSpectrum.xml'
  parameters=inputParameters(parameterFile)
  call parameters%markGlobal()
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Power spectra")
  ! Get required objects.
  darkMatterParticleCDM_    =  darkMatterParticleCDM   ()
  cosmologyParameters_      => cosmologyParameters     ()
  cosmologyFunctions_       => cosmologyFunctions      ()
  cosmologicalMassVariance_ => cosmologicalMassVariance()
  powerSpectrum_            => powerSpectrum           ()
  linearGrowth_             => linearGrowth            ()
  ! Build primordial power spectrum, transfer function, and transferred power spectrum.
  powerSpectrumPrimordialPowerLaw_         =                                                                   &
       & powerSpectrumPrimordialPowerLaw         (                                                             &
       &                                          index                   =1.000d0                           , &
       &                                          running                 =0.000d0                           , &
       &                                          wavenumberReference     =1.000d0                             &
       &                                         )
  transferFunctionEisensteinHu1999_        =                                                                   &
       & transferFunctionEisensteinHu1999        (                                                             &
       &                                          neutrinoNumberEffective =3.046d0                           , &
       &                                          neutrinoMassSummed      =0.000d0                           , &
       &                                          darkMatterParticle_     =darkMatterParticleCDM_            , &
       &                                          cosmologyParameters_    =cosmologyParameters_              , &
       &                                          cosmologyFunctions_     =cosmologyFunctions_                 &
       &                                         )
  powerSpectrumPrimordialTransferredSimple_=                                                                  &
       & powerSpectrumPrimordialTransferredSimple(                                                            &
       &                                          powerSpectrumPrimordial_=powerSpectrumPrimordialPowerLaw_ , &
       &                                          transferFunction_       =transferFunctionEisensteinHu1999_, &
       &                                          linearGrowth_           =linearGrowth_                      &
       &                                         )
  ! Test that σ₈ is correctly recovered.
  mass   =+4.0d0                                                      &
       &  /3.0d0                                                      &
       &  *Pi                                                         &
       &  *  cosmologyParameters_%OmegaMatter    (                  ) &
       &  *  cosmologyParameters_%densityCritical(                  ) &
       &  *(                                                          &
       &    +radiusNormalization                                      &
       &    /cosmologyParameters_%HubbleConstant (hubbleUnitsLittleH) &
       &   )**3
  call Assert('σ₈   consistency',cosmologicalMassVariance_%rootVariance(mass,cosmologyFunctions_%cosmicTime(1.0d0)),cosmologicalMassVariance_%sigma8(),relTol=1.0d-3)
  ! Test that σ(M) scales as expected.
  ratio=cosmologicalMassVariance_%rootVariance(1.0d10,cosmologyFunctions_%cosmicTime(1.0d0))/cosmologicalMassVariance_%rootVariance(1.0d12,cosmologyFunctions_%cosmicTime(1.0d0))
  call Assert('σ(M) scaling',ratio,100.0d0**((-1.0d0+3.0d0)/6.0d0),relTol=1.0d-6)
  ! Test power spectrum normalization. For a power-law n=-1 power spectrum, the integral over k²
  ! P(k) W²(k)/2π² can be computed analytically and is equal to 9/8π². We can therefore express
  ! this integral as 2π P(k₈) 9/8π²R₈³ where R₈=(8/h)Mpc and k₈=2π/R₈.
  sigma  =+sqrt(                                                                                &
       &        +powerSpectrum_%power(                                                          &
       &                              +2.0d0                                                    &
       &                              *Pi                                                       &
       &                              *cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)  &
       &                              /radiusNormalization                                    , &
       &                               cosmologyFunctions_ %cosmicTime    (1.0d0             )  &
       &                             )                                                          &
       &        *9.0d0                                                                          &
       &        /4.0d0                                                                          &
       &        /Pi                                                                             &
       &        /                    (                                                          &
       &                              +radiusNormalization                                      &
       &                              /cosmologyParameters_%hubbleConstant(hubbleUnitsLittleH)  &
       &                             )**3                                                       &
       &       )
  call Assert('P(k) normalization',sigma,cosmologicalMassVariance_%sigma8(),relTol=1.0d-6)
  ! Test that power spectrum scales as expected.
  ratio=powerSpectrum_%power(1.0d0,cosmologyFunctions_%cosmicTime(1.0d0))/powerSpectrum_%power(0.1d0,cosmologyFunctions_%cosmicTime(1.0d0))
  call Assert('P(k) scaling',ratio,0.1d0,relTol=1.0d-6)
  ! Do reproducibility tests.
  call Unit_Tests_Begin_Group("Reproducibility")
  do i=1,size(wavenumber)
     powerComputed           (i)=powerSpectrumPrimordialPowerLaw_         %power(wavenumber=wavenumber(i)                                           )
     transferComputed        (i)=transferFunctionEisensteinHu1999_        %value(wavenumber=wavenumber(i)                                           )
     powerTransferredComputed(i)=powerSpectrumPrimordialTransferredSimple_%power(wavenumber=wavenumber(i),time=cosmologyFunctions_%cosmicTime(1.0d0))
  end do
  call Assert('Pₚ(k)'      ,powerComputed           ,wavenumber                    ,relTol=1.0d-6)
  call Assert('      T (k)',transferComputed        ,transferExpected              ,relTol=1.0d-6)
  call Assert('Pₚ(k)·T²(k)',powerTransferredComputed,wavenumber*transferExpected**2,relTol=1.0d-6)
  call Unit_Tests_End_Group()
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()
end program Tests_Power_Spectrum
