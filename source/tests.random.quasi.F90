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
Contains a program to test quasi-random number functions.
!!}

program Test_Quasi_Random
  !!{
  Tests that quasi-random number functions work.
  !!}
  use :: Display                         , only : displayVerbositySet, verbosityLevelStandard
  use :: Numerical_Quasi_Random_Sequences, only : gsl_qrng_sobol     , quasiRandomNumberGenerator
  use :: Unit_Tests                      , only : Assert             , Unit_Tests_Begin_Group    , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  type            (quasiRandomNumberGenerator)               :: quasiRandomSequence
  double precision                            , dimension(7) :: r                  , rSobol
  integer                                                    :: i

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  quasiRandomSequence=quasiRandomNumberGenerator(gsl_qrng_sobol)
  call Unit_Tests_Begin_Group("quasi-random number sequences")
  do i=1,7
     r(i)=quasiRandomSequence%get()
  end do
  rSobol=[1.0d0/2.0d0,3.0d0/4.0d0,1.0d0/4.0d0,3.0d0/8.0d0,7.0d0/8.0d0,5.0d0/8.0d0,1.0d0/8.0d0]
  call Assert('Sobol sequence (1D)',r,rSobol,absTol=1.0d-6)
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Quasi_Random
