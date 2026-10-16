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
Contains a program which pins the :math:`r_\mathrm{LSS}` fitting function of the
:cite:t:`gnedin_effect_2000` filtering mass calculation.
!!}

program Test_Intergalactic_Medium_Filtering_Mass_Gnedin2000
  !!{RST
  Pins the :math:`r_\mathrm{LSS}` fitting function used by
  :galacticus-class:`intergalacticMediumFilteringMassGnedin2000`.

  The implementation exchanges the coefficients of the two inverse powers of expansion factor relative to the equation printed
  by :cite:t:`naoz_formation_2007`. That is deliberate, and it is their equation which is in error: as implemented the fit
  reproduces their figure precisely, whereas following their equation does not match that figure and gives
  :math:`|r_\mathrm{LSS}| > 1` at high redshift. The latter is unphysical here, because the ODE system integrated for the
  filtering mass uses :math:`(1+r_\mathrm{LSS})` as a factor multiplying the pressure term, so a value below :math:`-1` would
  reverse the sign of that term.

  This test therefore exists to stop the implementation being "corrected" back to the published ordering. It pins the values
  themselves, and separately asserts :math:`|r_\mathrm{LSS}| < 1` across the redshift and matter density ranges over which the
  fit is stated to be valid - the property the published ordering violates. The values below were evaluated directly from the
  coefficients of :cite:t:`naoz_formation_2007`, which the test shares; the load-bearing content is the pairing of each
  coefficient with its power of expansion factor, which it does not.
  !!}
  use :: Display                              , only : displayVerbositySet, verbosityLevelStandard
  use :: Intergalactic_Medium_Filtering_Masses, only : gnedin2000rLSS
  use :: Unit_Tests                           , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish, &
       &                                               compareLessThan
  implicit none
  integer                                  , parameter                                :: countDensities =3, countRedshifts=5
  double precision                         , dimension(               countDensities) :: densitiesMatter=[2.50000000000000d-01,3.00000000000000d-01,4.00000000000000d-01]
  double precision                         , dimension(countRedshifts               ) :: redshifts      =[7.00000000000000d+00,2.00000000000000d+01,5.00000000000000d+01,1.00000000000000d+02,1.50000000000000d+02]
  double precision                         , dimension(countRedshifts,countDensities) :: rLSSReference=reshape([                                                                                               &
       &                                                                                                       -2.38617008678272d-02,-5.92185167098079d-02,-1.34059478269306d-01,-2.44743361315984d-01,-3.42376469350414d-01, &
       &                                                                                                       -2.44107190325231d-02,-6.15496119855029d-02,-1.39857693249380d-01,-2.54977027569293d-01,-3.55798117701402d-01, &
       &                                                                                                       -2.50863982016537d-02,-6.47670112943312d-02,-1.47942748566348d-01,-2.69098196827547d-01,-3.74026588811937d-01  &
       &                                                                                                      ],[countRedshifts,countDensities])
  double precision                         , dimension(countRedshifts               ) :: rLSS
  double precision                                                                    :: rLSSMagnitudeMaximum
  double precision                         , parameter                                :: tolerance           =1.0d-12
  integer                                                                             :: i                           , j

  call displayVerbositySet(verbosityLevelStandard)
  call Unit_Tests_Begin_Group("Gnedin (2000) filtering mass: r_LSS")
  rLSSMagnitudeMaximum=0.0d0
  do j=1,countDensities
     do i=1,countRedshifts
        rLSS(i)             =gnedin2000rLSS(densitiesMatter(j),1.0d0/(1.0d0+redshifts(i)))
        rLSSMagnitudeMaximum=max(rLSSMagnitudeMaximum,abs(rLSS(i)))
     end do
     call Assert('r_LSS against the coefficients of Naoz & Barkana (2007)',rLSS,rLSSReference(:,j),relTol=tolerance)
  end do
  ! The published ordering of the coefficients violates this; the implemented ordering does not.
  call Assert('|r_LSS| < 1 over the validated range',rLSSMagnitudeMaximum,1.0d0,compareLessThan)
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Test_Intergalactic_Medium_Filtering_Mass_Gnedin2000
