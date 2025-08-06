!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
Contains a program which tests regular expression functionality.
!!}

program Tests_Regular_Expressions
  !!{
  Tests regular expression functionality.
  !!}
  use :: Display            , only : displayVerbositySet, verbosityLevelStandard
  use :: Regular_Expressions, only : regEx
  use :: Unit_Tests         , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Regular expressions")

  ! Test reg-ex matching.
  block
    type(regEx) :: regEx_
    regEx_=regEx("stellarPopulationFileFor[a-zA-Z0-9]+IMF")
    call Assert("Matches",regEx_%matches('nothingToSeeHereMoveAlong'          ),.false.)
    call Assert("Matches",regEx_%matches('stellarPopulationFileForSalpeterIMF'),.true. )
  end block

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

end program Tests_Regular_Expressions
