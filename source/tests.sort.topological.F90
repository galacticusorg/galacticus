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
Contains a program to test sorting functions.
!!}

program Test_Sort_Topological
  !!{
  Tests of topological sorting functions.
  !!}
  use :: Display            , only : displayVerbositySet, verbosityLevelStandard
  use :: Error   , only : errorStatusFail    , errorStatusSuccess
  use :: Sorting_Topological, only : Sort_Topological
  use :: Unit_Tests         , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  integer, parameter                      :: countObjects=10, countDependencies=5
  integer, dimension(countObjects       ) :: order
  integer, dimension(countDependencies,2) :: dependencies
  integer                                 :: countOrdered   , status             , &
       &                                     i              , j
  logical                                 :: success

  call displayVerbositySet(verbosityLevelStandard)
  call Unit_Tests_Begin_Group("Topological sorting")
  ! A case with circular dependencies.
  dependencies(1,:)=[1,2]
  dependencies(2,:)=[6,3]
  dependencies(3,:)=[6,4]
  dependencies(4,:)=[1,9]
  dependencies(5,:)=[9,1]
  call Sort_Topological(countObjects,countDependencies,dependencies,order,countOrdered,status)
  call Assert('circular dependency detected',status == errorStatusFail   ,.true.)
  ! A case without circular dependencies.
  dependencies(1,:)=[1,2]
  dependencies(2,:)=[6,3]
  dependencies(3,:)=[6,4]
  dependencies(4,:)=[1,9]
  dependencies(5,:)=[9,6]
  call Sort_Topological(countObjects,countDependencies,dependencies,order,countOrdered,status)
  call Assert('sort succeeds'               ,status == errorStatusSuccess,.true.)
  ! Test that ordering is correct.
  success=.true.
  do i=1,countDependencies
     j=1
     do while (order(j) /= dependencies(i,1))
        j=j+1
     end do
     if (.not.any(order(1:j-1) == dependencies(i,2))) success=.false.
  end do
  call Assert('ordering correct'            ,                     success,.true.)
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Test_Sort_Topological
