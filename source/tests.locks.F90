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
  use            :: OMP_Lib           , only : OMP_Get_Thread_Num   , OMP_Set_Max_Active_Levels, OMP_Get_Supported_Active_Levels, omp_get_ancestor_thread_num, omp_get_level  , &
       &                                       OMP_Get_Num_Threads
  use            :: Array_Utilities   , only : Array_Is_Monotonic   , directionIncreasing
  use            :: Display           , only : displayMessage       , displayVerbositySet      , verbosityLevelStandard         , displayIndent              , displayUnindent
  use, intrinsic :: ISO_C_Binding     , only : c_size_t
  use            :: ISO_Varying_String, only : operator(//)         , var_str                  , varying_string
  use            :: Locks             , only : ompIncrementalLock   , ompReadWriteLock         , mutex
  use            :: String_Handling   , only : operator(//)
  use            :: Unit_Tests        , only : Assert               , Unit_Tests_Begin_Group   , Unit_Tests_End_Group           , Unit_Tests_Finish
  use            :: Events_Hooks      , only : eventsHooksInitialize
  implicit none
  integer         (c_size_t          ), parameter               :: elementCount   =100_c_size_t
  integer         (c_size_t          ), dimension(elementCount) :: orderedCount
  type            (ompIncrementalLock)                          :: incrementalLock
  type            (ompReadWriteLock  )                          :: readWriteLock
  type            (mutex             )                          :: lock
  integer         (c_size_t          )                          :: i                           , unorderedCounter       , &
       &                                                           orderedCounter              , unorderedCounterPrivate
  double precision                                              :: uniformRandom
  integer                                                       :: sleepTime
  type            (varying_string    )                          :: message
  logical                                                       :: lockFailed                  , lockHeld

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Allow nested parallelism.
  call OMP_Set_Max_Active_Levels(OMP_Get_Supported_Active_Levels())
  ! Initialize event hooks.
  call eventsHooksInitialize()
  call Unit_Tests_Begin_Group("OpenMP")
  ! Test read/write locks. First obtain a write lock and ensure that only one thread at a time holds it.
  call displayIndent('Read/write lock...')
  readWriteLock=ompReadWriteLock()
  i            =0
  lockFailed   =.false.
  !$omp parallel private(message)
  call readWriteLock%setWrite(haveReadLock=.false.)
  message=var_str("thread ")//OMP_get_Thread_Num()//' has write lock'
  call displayMessage(message)
  !$omp atomic
  i=i+1
  call sleep(1)
  if (i /= 1) lockFailed=.true.
  !$omp atomic
  i=i-1
  message=var_str("thread ")//OMP_Get_Thread_Num()//' releasing write lock'
  call displayMessage(message)
  call readWriteLock%unsetWrite(haveReadLock=.false.)
  !$omp end parallel
  call Assert('OpenMP read/write lock [write]',lockFailed,.false.)
  lockFailed=.false.
  i         =0
  !$omp parallel
  call readWriteLock%setRead()
  message=var_str("thread ")//OMP_Get_Thread_Num()//' has read lock'
  call displayMessage(message)
  !$omp atomic
  i=i+1
  !$omp barrier
  !$omp single
  if (i /= OMP_Get_Num_Threads()) lockFailed=.true.
  !$omp end single
  message=var_str("thread ")//OMP_get_Thread_Num()//' has releasing read lock'
  call readWriteLock%unsetRead()
  !$omp end parallel
  call Assert('OpenMP read/write lock [read]',lockFailed,.false.)
  call displayUnindent('...done')
  ! Test incremental locks. We generate a counter which is not guaranteed to be ordered in terms of OpenMP threads. Then we use an
  ! incremental lock to force the access back to being ordered by counter value and thereby construct and ordered array.
  call displayIndent('Incremental lock...')
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
  call displayUnindent('...done')
  call Assert('OpenMP incremental lock',Array_Is_Monotonic(orderedCount,direction=directionIncreasing,allowEqual=.false.),.true.)
  ! Nested parallelism test.
  call displayUnindent('Nested parallelism lock...')
  lock      =mutex(recursiveLock=.true.)
  lockFailed=.false.
  lockHeld  =.false.
  !$omp parallel
  ! Make thread 0 obtain the lock.
  if (OMP_Get_Thread_Num() == 0) then
     call displayMessage("initial thread is acquiring the lock in outer parallel region")
     call lock%set()
     call displayMessage("initial thread has acquired the lock in outer parallel region")
     lockHeld=.true.
     call lock%  set()
     call displayMessage("initial thread has acquired the lock recursively in outer parallel region")
     call lock%unset()
     call displayMessage("initial thread has released the lock recursively in outer parallel region")
  end if
  !$omp barrier
  ! Start a nested parallel region.
  !$omp parallel private(message)
  if (OMP_Get_Ancestor_Thread_Num(OMP_Get_Level()-1) == 0 .and. OMP_Get_Thread_Num() == 0) then
     call displayMessage("initial thread is sleeping for 10s in inner parallel region")
     call sleep(10)
     lockHeld=.false.
     call lock%unset()
     call displayMessage("initial thread has released the lock in inner parallel region")
  end if
  call lock%set  ()
  message=var_str("thread ")//OMP_Get_Ancestor_Thread_Num(OMP_Get_Level()-1)//":"//OMP_Get_Thread_Num()//" has acquired the lock in inner parallel region"
  call displayMessage(message)
  if (lockHeld) lockFailed=.true.
  call lock%unset()
  !$omp end parallel
  !$omp barrier
  !$omp end parallel
  call displayUnindent('...done')
  call Assert('OpenMP lock is nested parallel region',lockFailed,.false.)
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Test_Locks
