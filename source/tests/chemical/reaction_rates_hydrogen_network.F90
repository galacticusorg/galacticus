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
Contains a program to test the primordial hydrogen chemistry network against the published fits it implements.
!!}

program Test_Chemical_Reaction_Rates_Hydrogen_Network
  !!{RST
  Tests the rate coefficients and photo cross-sections of
  :galacticus-class:`chemicalReactionRateHydrogenNetwork` against an independent transcription of the fits of
  :cite:t:`abel_modeling_1997` (``hydrogenNetworkCheck.py`` in `galacticusDevTools
  <https://github.com/galacticusorg/galacticusDevTools>`_).

  Four rate coefficients are covered - their reactions 7, 8, 14 and 16 - at temperatures spanning the branches of each fit,
  including the two-branch structure of reaction 7 at :math:`6000` K and the constant low-temperature value of reaction 8 below
  :math:`0.1` eV.

  The photo cross-sections of the same file are **not** covered, for a reason which is not a matter of choice: exporting any one
  of them segfaults gfortran 16 while generating code for their module, as they hold ``save``, ``!$omp threadprivate``
  derived-type interpolation tables. Their fitting formulae were instead checked against :cite:t:`abel_modeling_1997` by hand -
  reactions 24, 26 and 28, the last in both the Lyman and Werner bands and for both the para and ortho states - and agree
  exactly. Two further cross-sections could not be checked against that paper at all:
  :math:`\hbox{H}_2^++\gamma \rightarrow \hbox{H}+\hbox{H}^+`, which Galacticus takes from :cite:t:`shapiro_hydrogen_1987`
  while :cite:t:`abel_modeling_1997` fit different data, and :math:`\hbox{H}+\gamma \rightarrow \hbox{H}^++\hbox{e}^-`,
  which is delegated to the photoionization cross-section class rather than fitted here.
  !!}
  use :: Chemical_Reaction_Rates     , only : hydrogenNetworkH_Electron_to_Hminus_Photon_RateCoefficient   , hydrogenNetworkH_Hminus_to_H2_Electron_RateCoefficient, &
       &                                      hydrogenNetworkHminus_Electron_to_H_2Electron_RateCoefficient, hydrogenNetworkHminus_Hplus_to_2H_RateCoefficient
  use :: Display                     , only : displayVerbositySet                                          , verbosityLevelStandard
  use :: Numerical_Constants_Physical, only : boltzmannsConstant
  use :: Numerical_Constants_Units   , only : electronVolt
  use :: Unit_Tests                  , only : Assert                                                       , Unit_Tests_Begin_Group                               , &
       &                                      Unit_Tests_End_Group                                         , Unit_Tests_Finish
  implicit none
  integer         , parameter                                 :: countTemperaturesKelvin=6, countTemperaturesElectronVolts=5
  double precision, dimension(countTemperaturesKelvin       ) :: temperaturesKelvin         =[1.00000000000000d+01,1.00000000000000d+02,1.00000000000000d+03,6.00000000000000d+03,1.00000000000000d+04,1.00000000000000d+05]
  double precision, dimension(countTemperaturesKelvin       ) :: rateCoefficientK7Reference =[1.08790814632211d-17,1.06245911999513d-16,8.46744196646279d-16,2.76872355994436d-15,3.53118078441984d-15,6.34047051840884d-15]
  double precision, dimension(countTemperaturesKelvin       ) :: rateCoefficientK16Reference=[2.21359436211787d-07,7.00000000000000d-08,2.21359436211787d-08,9.03696114115064d-09,7.00000000000000d-09,2.21359436211787d-09]
  double precision, dimension(countTemperaturesElectronVolts) :: temperaturesElectronVolts  =[1.00000000000000d-02,1.00000000000000d-01,1.00000000000000d+00,1.00000000000000d+01,1.00000000000000d+02]
  double precision, dimension(countTemperaturesElectronVolts) :: rateCoefficientK8Reference =[1.42800000000000d-09,1.42810817222593d-09,1.92346234542515d-09,3.71980645412309d-09,7.68952118929211d-09]
  double precision, dimension(countTemperaturesElectronVolts) :: rateCoefficientK14Reference=[7.21516825664962d-43,1.62423969648126d-12,1.49509149247854d-08,6.29202054563722d-07,8.98005780257346d-07]
  double precision, dimension(countTemperaturesKelvin       ) :: rateCoefficientK7          , rateCoefficientK16
  double precision, dimension(countTemperaturesElectronVolts) :: rateCoefficientK8          , rateCoefficientK14
  ! Both implementations evaluate the same expressions, so they are limited only by the order of the arithmetic.
  double precision, parameter                                 :: tolerance                  =1.0d-9
  integer                                                     :: i

  call displayVerbositySet(verbosityLevelStandard)
  call Unit_Tests_Begin_Group("Primordial hydrogen network")
  do i=1,countTemperaturesKelvin
     rateCoefficientK7 (i)=hydrogenNetworkH_Electron_to_Hminus_Photon_RateCoefficient(temperaturesKelvin(i))
     rateCoefficientK16(i)=hydrogenNetworkHminus_Hplus_to_2H_RateCoefficient         (temperaturesKelvin(i))
  end do
  do i=1,countTemperaturesElectronVolts
     ! These fits are expressed in electron volts, but the interfaces take Kelvin.
     rateCoefficientK8 (i)=hydrogenNetworkH_Hminus_to_H2_Electron_RateCoefficient       (temperaturesElectronVolts(i)*electronVolt/boltzmannsConstant)
     rateCoefficientK14(i)=hydrogenNetworkHminus_Electron_to_H_2Electron_RateCoefficient(temperaturesElectronVolts(i)*electronVolt/boltzmannsConstant)
  end do
  call Assert('k₇  : H + e⁻ → H⁻ + γ'  ,rateCoefficientK7 ,rateCoefficientK7Reference ,relTol=tolerance)
  call Assert('k₈  : H⁻ + H → H₂ + e⁻' ,rateCoefficientK8 ,rateCoefficientK8Reference ,relTol=tolerance)
  call Assert('k₁₄ : H⁻ + e⁻ → H + 2e⁻',rateCoefficientK14,rateCoefficientK14Reference,relTol=tolerance)
  call Assert('k₁₆ : H⁻ + H⁺ → 2H'     ,rateCoefficientK16,rateCoefficientK16Reference,relTol=tolerance)
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Test_Chemical_Reaction_Rates_Hydrogen_Network
