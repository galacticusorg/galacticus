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
Contains a program to test perfect hashing algorithms.
!!}

program Test_Perfect_Hashes
  !!{
  Tests perfect hashing algorithms.
  !!}
  use :: Display          , only : displayVerbositySet, verbosityLevelStandard
  use :: Hashes_Perfect   , only : hashPerfect
  use :: Kind_Numbers     , only : kind_int8
  use :: Unit_Tests       , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  integer                , parameter                        :: keyCount   =11
  integer(kind=kind_int8)                                   :: i
  integer(kind=kind_int8)             , dimension(keyCount) :: keys          , retrievedValues, values
  integer                , allocatable, dimension(:)        :: bucketCount
  logical                             , dimension(keyCount) :: keyPresent
  type   (hashPerfect   )                                   :: hash

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Perfect hashes")

  ! Create a list of keys.
  keys  =[ 6,28,65,14,88,184,38,44,523,12,98]

  ! Create a list of values.
  values=[99,88,77,66,55, 44,33,22, 11,86,34]

  ! Create the hash function.
  call hash%create(keys,values)

  ! Allocate arrays.
  allocate(bucketCount(0:hash%size()-1))

  ! Look up indices and presence.
  bucketCount=0
  do i=1,keyCount
     keyPresent(i)=hash%isPresent(keys(i))
     bucketCount(hash%index(keys(i)))=bucketCount(hash%index(keys(i)))+1
     retrievedValues(i)=hash%value(keys(i))
  end do

  ! Make assertions.
  call Assert('All keys present in hash',all(keyPresent                             ),.true.)
  call Assert('No hash collisions'      ,all(bucketCount >= 0 .and. bucketCount <= 1),.true.)
  call Assert('All values correct'      ,retrievedValues                             ,values)

  ! Destroy the hash.
  call hash%destroy()

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

  ! Clean up.
  deallocate(bucketCount)
  
end program Test_Perfect_Hashes
