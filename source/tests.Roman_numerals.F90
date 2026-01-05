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
Contains a program to test Roman numeral conversion functions.
!!}

program Test_Roman_Numerals
  !!{
  Tests that Roman numeral conversion functions work.
  !!}
  use :: Display                 , only : displayVerbositySet, verbosityLevelStandard
  use :: ISO_Varying_String      , only : var_str
  use :: Numerical_Roman_Numerals, only : Roman_Numerals
  use :: Unit_Tests              , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none

  call displayVerbositySet(verbosityLevelStandard                                                    )
  call Unit_Tests_Begin_Group        ("Roman numeral conversion"                                           )
  call Assert                        ('1987 = MCMLXXXVII'       ,Roman_Numerals(1987),var_str('MCMLXXXVII'))
  call Unit_Tests_End_Group          (                                                                     )
  call Unit_Tests_Finish             (                                                                     )
end program Test_Roman_Numerals
