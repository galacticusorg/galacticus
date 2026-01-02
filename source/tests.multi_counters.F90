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
Contains a program to test multi-counters.
!!}

program Test_Multi_Counters
  !!{
  Tests of multi-counters.
  !!}
  use            :: Display       , only : displayVerbositySet, verbosityLevelStandard
  use, intrinsic :: ISO_C_Binding , only : c_size_t
  use            :: Multi_Counters, only : multiCounter
  use            :: Unit_Tests    , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  type   (multiCounter)                   :: counter1, counter2, &
       &                                     counter3
  logical              , dimension(2,5,3) :: state
  integer(c_size_t    ), dimension(0    ) :: zeroSize
  integer(c_size_t    )                   :: i       , j       , &
       &                                     k

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Multi-counters")

  ! Construct the counters.
  counter1=multiCounter([2_c_size_t,5_c_size_t,3_c_size_t])
  counter2=multiCounter(                                  )
  counter3=multiCounter(zeroSize                          )
  call counter2%append(2_c_size_t)
  call counter2%append(5_c_size_t)
  call counter2%append(3_c_size_t)

  ! Iterate over first counter.
  state=.false.
  do while (counter1%increment())
     i=counter1%state(1_c_size_t)
     j=counter1%state(2_c_size_t)
     k=counter1%state(3_c_size_t)
     state(i,j,k)=.true.
  end do
  call Assert('coverage preset ranges'  ,all(state),.true.)

  ! Iterate over second counter.
  state=.false.
  do while (counter2%increment())
     i=counter2%state(1_c_size_t)
     j=counter2%state(2_c_size_t)
     k=counter2%state(3_c_size_t)
     state(i,j,k)=.true.
  end do
  call Assert('coverage appended ranges',all(state),.true.)

  ! Iterate over third counter.
  i=0
  do while (counter3%increment())
     i=i+1
  end do
  call Assert('zero dimension counter',i,1_c_size_t)

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

end program Test_Multi_Counters
