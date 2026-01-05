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
Contains a program to test debugging tools.
!!}

program Test_Debugging
  !!{
  Contains a program to test debugging tools.
  !!}
  use :: Display                 , only : displayVerbositySet, verbosityLevelStandard
  use :: Test_Debugging_Functions, only : dummyCaller
  use :: ISO_Varying_String      , only : varying_string     , char                  , var_str
#ifdef DEBUGGING
  use :: Unit_Tests              , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
#else
  use :: Unit_Tests              , only : Skip               , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
#endif  
  implicit none
#ifdef DEBUGGING
  type(varying_string) :: caller
#endif

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Debugging")
#ifdef DEBUGGING
  call dummyCaller(caller)
  call Assert("Identify calling function",caller,var_str("main")      )
#else
  call Skip  ("Identify calling function","not compiled for debugging")
#endif
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Test_Debugging
