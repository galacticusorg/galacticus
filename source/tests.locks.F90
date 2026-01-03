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
Contains a program to test OpenMP lock functions.
!!}

program Test_Locks
  !!{
  Tests of OpenMP locking functions.
  !!}
  use            :: Array_Utilities   , only : Array_Is_Monotonic   , directionIncreasing
  use            :: Display           , only : displayMessage       , displayVerbositySet   , verbosityLevelStandard
  use, intrinsic :: ISO_C_Binding     , only : c_size_t
  use            :: ISO_Varying_String, only : operator(//)         , var_str               , varying_string
  use            :: Locks             , only : ompIncrementalLock
  use            :: String_Handling   , only : operator(//)
  use            :: Unit_Tests        , only : Assert               , Unit_Tests_Begin_Group, Unit_Tests_End_Group  , Unit_Tests_Finish
  use            :: Events_Hooks      , only : eventsHooksInitialize
  implicit none
  integer         (c_size_t          ), parameter               :: elementCount   =100_c_size_t
  integer         (c_size_t          ), dimension(elementCount) :: orderedCount
  type            (ompIncrementalLock)                          :: incrementalLock
  integer         (c_size_t          )                          :: i                           , unorderedCounter       , &
       &                                                           orderedCounter              , unorderedCounterPrivate
  double precision                                              :: uniformRandom
  integer                                                       :: sleepTime
  type            (varying_string    )                          :: message

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Initialize event hooks.
  call eventsHooksInitialize()
  call Unit_Tests_Begin_Group("OpenMP")
  ! Test incremental locks. We generate a counter which is not guaranteed to be ordered in terms of OpenMP threads. Then we use an
  ! incremental lock to force the access back to being ordered by counter value and thereby construct and ordered array.
  call random_seed()
  incrementalLock =ompIncrementalLock()
  orderedCounter  =0_c_size_t
  unorderedCounter=0_c_size_t
  !$omp parallel do private(unorderedCounterPrivate,message,uniformRandom,sleepTime)
  do i=1_c_size_t,elementCount
     !$omp critical(increment)
     unorderedCounter       =unorderedCounter+1_c_size_t
     unorderedCounterPrivate=unorderedCounter
     !$omp end critical(increment)
     ! Sleep for a random duration to force desynchronization of threads.
     call random_number(uniformRandom)
     sleepTime=int(uniformRandom*3.0d0)
     message=var_str("unordered counter retrieved with value ")//unorderedCounterPrivate//var_str("; pausing for ")//sleepTime//var_str(" seconds")
     call displayMessage(message)
     call sleep        (sleepTime)
     call incrementalLock%set  (unorderedCounterPrivate)
     orderedCounter                =orderedCounter         +1_c_size_t
     message=var_str("ordered counter retrieved with value ")//orderedCounter
     call displayMessage(message)
     orderedCount  (orderedCounter)=unorderedCounterPrivate
     call incrementalLock%unset( )
  end do
  !$omp end parallel do
  call Assert('OpenMP incremental lock',Array_Is_Monotonic(orderedCount,direction=directionIncreasing,allowEqual=.false.),.true.)
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Test_Locks
