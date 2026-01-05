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
Contains a program to test mathematical arithmetic functions.
!!}

program Test_Math_Arithmetic
  !!{
  Tests of mathematical arithmetic functions.
  !!}
  use :: Display        , only : displayVerbositySet, verbosityLevelStandard
  use :: Kind_Numbers   , only : kind_int8
  use :: Math_Arithmetic, only : divideSafe
  use :: Unit_Tests     , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Math: arithmetic functions")
  ! A simple division that will not overflow.
  call Assert("1.0/0.001",divideSafe(+     1.0d+0 ,+1.0d-3               ),+     1.0d+3 ,relTol=1.0d-12)
  ! Dividing the largest representable number, by the largest representable number less than 1 should lead to an overflow, but
  ! should be caught by safe division and just return the largest representable number.
  call Assert("+H/(+1+ε)",divideSafe(+huge(1.0d+0),+1.0d+0-epsilon(1.0d0)),+huge(1.0d+0)               )
  ! As above, but ensuring that sign handling is correct.
  call Assert("-H/(+1+ε)",divideSafe(-huge(1.0d+0),+1.0d+0-epsilon(1.0d0)),-huge(1.0d+0)               )
  call Assert("-H/(-1-ε)",divideSafe(-huge(1.0d+0),-1.0d+0+epsilon(1.0d0)),+huge(1.0d+0)               )
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Test_Math_Arithmetic
