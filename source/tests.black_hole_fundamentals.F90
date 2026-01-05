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
Contains a program to test the black hole fundamental functions.
!!}

program Test_Black_Hole_Fundamentals
  !!{
  Tests of black hole fundamental functions.
  !!}
  use :: Black_Hole_Fundamentals, only : Black_Hole_Horizon_Radius, Black_Hole_ISCO_Radius, orbitPrograde
  use :: Display                , only : displayVerbositySet      , verbosityLevelStandard
  use :: Unit_Tests             , only : Assert                   , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish, &
          &                              compareEquals
  implicit none

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Black hole functions")

  ! ISCO radius for a Schwarzchild black hole should be 6 for prograde orbits.
  call Assert("Schwarzchild metric ISCO radius"   ,Black_Hole_ISCO_Radius   (0.0d0,orbitPrograde),6.0d0,compareEquals,1.0d-6)

  ! ISCO radius for an extreme Kerr black hole should be 1 for prograde orbits.
  call Assert("Extreme Kerr metric ISCO radius"   ,Black_Hole_ISCO_Radius   (1.0d0,orbitPrograde),1.0d0,compareEquals,1.0d-6)

  ! Horizon radius for a Schwarzchild black hole should be 2.
  call Assert("Schwarzchild metric horizon radius",Black_Hole_Horizon_Radius(0.0d0              ),2.0d0,compareEquals,1.0d-6)

  ! Horizon radius for an extreme Kerr black hole should be 1.
  call Assert("Extreme Kerr metric horizon radius",Black_Hole_Horizon_Radius(1.0d0              ),1.0d0,compareEquals,1.0d-6)

  ! End the unit testing.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Black_Hole_Fundamentals
