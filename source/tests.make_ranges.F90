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
Contains a program to test the numerical range making code.
!!}

program Test_Make_Ranges
  !!{
  Tests that numerical range making code works correctly.
  !!}
  use :: Array_Utilities , only : Array_Reverse
  use :: Display         , only : displayVerbositySet, verbosityLevelStandard
  use :: Numerical_Ranges, only : Make_Range         , rangeTypeLinear       , rangeTypeLogarithmic
  use :: Unit_Tests      , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  double precision, dimension(1:11) :: range1
  double precision, dimension(0:10) :: range2

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Numerical ranges")

  ! Create a linear range.
  range1=Make_Range(1.0d0,3.0d0,size(range1),rangeType=rangeTypeLinear)
  call Assert("linear range creation",range1,[1.0d0,1.2d0,1.4d0,1.6d0,1.8d0,2.0d0,2.2d0,2.4d0,2.6d0,2.8d0,3.0d0])

  ! Create a linear range in an offset array.
  range2=Make_Range(1.0d0,3.0d0,size(range2),rangeType=rangeTypeLinear)
  call Assert("linear range creation (offset array)",range2,[1.0d0,1.2d0,1.4d0,1.6d0,1.8d0,2.0d0,2.2d0,2.4d0,2.6d0,2.8d0,3.0d0])

  ! Create a logarithmic range.
  range1=Make_Range(10.0d0,1000.0d0,size(range1),rangeType=rangeTypeLogarithmic)
  call Assert("logarithmic range creation",range1,10.0d0**([1.0d0,1.2d0,1.4d0,1.6d0,1.8d0,2.0d0,2.2d0,2.4d0,2.6d0,2.8d0,3.0d0]),relTol=1.0d-6)

  ! Create a reverse order range.
  range1=Make_Range(3.0d0,1.0d0,size(range1),rangeType=rangeTypeLinear)
  call Assert("reverse linear range creation",range1,Array_Reverse([1.0d0,1.2d0,1.4d0,1.6d0,1.8d0,2.0d0,2.2d0,2.4d0,2.6d0,2.8d0,3.0d0]))

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Make_Ranges
