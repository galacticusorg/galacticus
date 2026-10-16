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
Contains a program to test the :cite:t:`meiksin_colour_2006` algorithm for intergalactic attenuation.
!!}

program Test_Meiksin2006
  !!{RST
  Tests :galacticus-class:`stellarPopulationSpectraPostprocessorMeiksin2006` against an independent implementation of
  :cite:t:`meiksin_colour_2006` written from that paper (``meiksin2006Check.py`` in `galacticusDevTools
  <https://github.com/galacticusorg/galacticusDevTools>`_).

  Wavelengths are specified as a fraction of the Lyman limit wavelength, :math:`\lambda/\lambda_\mathrm{L}`, in the rest frame.
  Each redshift then covers the same three regimes: the Lyman continuum, where the optically thick and optically thin terms
  contribute, for fractions below unity; the Lyman series alone between :math:`1` and :math:`4/3`; and, beyond :math:`4/3`,
  redward of Lyman-:math:`\alpha`, where nothing absorbs and the transmission must be exactly unity. Expressing the wavelengths
  this way also makes the comparison independent of the value adopted for :math:`\lambda_\mathrm{L}`, which carries a
  reduced-mass correction.

  The optically thick term is the delicate part of the model, and three details of it are easy to get wrong in ways which
  partly mask one another: :math:`\Gamma(2-\beta,1)` is the incomplete :math:`\Gamma` function rather than
  :math:`\Gamma(2-\beta)`; the series alternate as :math:`(-1)^n`, which in Fortran must be written ``(-1)**n``, since ``**``
  binds more tightly than unary minus; and the two series begin at :math:`n=0` and :math:`n=1` respectively. With all three
  right the bracketed factor reduces analytically to :math:`\Gamma(2-\beta)`, which the reference script checks.

  The transmission is also required to lie in :math:`(0,1]`. That bound is not a formality: it is violated over a large part of
  the plane if the alternating signs are corrected without the other two details, so it discriminates between a correct
  treatment of these terms and a plausible-looking one.
  !!}
  use :: Display                               , only : displayVerbositySet                            , verbosityLevelStandard
  use :: Numerical_Constants_Atomic            , only : lymanSeriesLimitWavelengthHydrogen_atomic
  use :: Stellar_Population_Spectra_Postprocess, only : stellarPopulationSpectraPostprocessorMeiksin2006
  use :: Unit_Tests                            , only : Assert                                         , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish, &
       &                                                compareLessThanOrEqual
  implicit none
  type            (stellarPopulationSpectraPostprocessorMeiksin2006)                      :: postprocessor
  integer                                                           , parameter           :: countRedshifts       =4        , countFractions=6
  double precision                                                  , dimension(countRedshifts) :: redshifts      =[1.00000000000000d+00,2.00000000000000d+00,3.00000000000000d+00,5.00000000000000d+00]
  ! Rest-frame wavelengths, as a fraction of the Lyman limit wavelength. 4/3 is Lyman-α.
  double precision                                                  , dimension(countFractions) :: wavelengthFractions=[3.00000000000000d-01,6.00000000000000d-01,9.00000000000000d-01,1.05000000000000d+00,1.20000000000000d+00,1.40000000000000d+00]
  double precision                                                  , dimension(countFractions,countRedshifts) :: transmissionReference=reshape([                                                                                                &
       &                                                                                                                                        6.49663288410986d-01,4.63373075703678d-01,6.74668971494205d-01,9.79492845348441d-01,9.81602249025122d-01,1.00000000000000d+00, &
       &                                                                                                                                        3.37420885794525d-01,1.50699251736674d-01,3.72332785962650d-01,9.05766475450778d-01,9.20130940007982d-01,1.00000000000000d+00, &
       &                                                                                                                                        1.21746784235458d-01,2.69163581199674d-02,1.42152427165103d-01,7.40077958921068d-01,7.85586534209925d-01,1.00000000000000d+00, &
       &                                                                                                                                        4.56329566301189d-03,1.16691066077264d-04,4.87971330232800d-03,2.33896399696313d-01,3.17893022211617d-01,1.00000000000000d+00  &
       &                                                                                                                                       ],[countFractions,countRedshifts])
  double precision                                                  , dimension(countFractions) :: transmission                   , transmissionZeroRedshift
  double precision                                                                        :: transmissionMaximum
  ! The two implementations evaluate the same expressions, so they are limited only by the order of the arithmetic.
  double precision                                                  , parameter           :: tolerance            =1.0d-9
  integer                                                                                 :: i                             , j

  call displayVerbositySet(verbosityLevelStandard)
  call Unit_Tests_Begin_Group("Meiksin (2006) IGM attenuation model")
  postprocessor=stellarPopulationSpectraPostprocessorMeiksin2006()
  ! At zero redshift there is no intervening medium, so nothing is absorbed.
  do i=1,countFractions
     transmissionZeroRedshift(i)=postprocessor%multiplier(wavelengthFractions(i)*lymanSeriesLimitWavelengthHydrogen_atomic,0.0d0,0.0d0)
  end do
  call Assert('no attenuation at zero redshift',transmissionZeroRedshift,spread(1.0d0,1,countFractions),relTol=tolerance)
  ! Compare with the independent implementation, and require the transmission to be physical.
  transmissionMaximum=0.0d0
  do j=1,countRedshifts
     do i=1,countFractions
        transmission(i)=postprocessor%multiplier(wavelengthFractions(i)*lymanSeriesLimitWavelengthHydrogen_atomic,0.0d0,redshifts(j))
        transmissionMaximum=max(transmissionMaximum,transmission(i))
     end do
     call Assert('transmission against an independent implementation',transmission,transmissionReference(:,j),relTol=tolerance)
  end do
  call Assert('transmission never exceeds unity',transmissionMaximum,1.0d0,compareLessThanOrEqual)
  ! Redward of Lyman-α no transition absorbs, so the transmission is exactly unity.
  do j=1,countRedshifts
     transmission(1)=postprocessor%multiplier(1.4d0*lymanSeriesLimitWavelengthHydrogen_atomic,0.0d0,redshifts(j))
     call Assert('unit transmission redward of Lyman-α',transmission(1),1.0d0,relTol=tolerance)
  end do
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Test_Meiksin2006
