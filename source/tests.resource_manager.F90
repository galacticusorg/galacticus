!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
Contains a program for testing the {\normalfont \ttfamily resourceManager} class.
!!}
  
program Test_Resource_Manager
  !!{
  Test the {\normalfont \ttfamily resourceManager} class.
  !!}
  use :: Test_Resource_Manager_Wrapper, only : resourceHolder
  use :: Display                      , only : displayVerbositySet, verbosityLevelStandard
  use :: Unit_Tests                   , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none

  ! The expectation is that the following sequence should result in 5 calls to the destructor:
  !  1. Destruct of "self" in 1st-level constructor before assignment.
  !  2. Destruct of function result from 2nd-level constructor after assignment.
  !  3. Destruct of "holder" before assignment.
  !  4. Destruct of the "resourceHolder()" function result after assignment.
  !  5. Destruct of "holder" at the end of the block.
  call displayVerbositySet(verbosityLevelStandard)
  call Unit_Tests_Begin_Group("Resource manager")
  block
    ! Use a block here to ensure that "holder" goes out of scope before the end of the block.
    type(resourceHolder) :: holder
    holder=resourceHolder()
    call Assert("count of pointers to shared resource",holder%manager%count(),1)
  end block
  call Unit_Tests_Finish()

end program Test_Resource_Manager
