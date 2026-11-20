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

!+    Contributions to this file made by: Claude.

!!{RST
Contains a program to test nearest neighbor searches, including searches made concurrently from OpenMP threads.
!!}

program Test_Nearest_Neighbors
  !!{RST
  Tests of nearest neighbor searches. Searches are made serially, then concurrently from OpenMP threads (both of a shared tree,
  and of a tree built by each thread, as the N-body operators do), and the results compared. The ANN library keeps its working
  state in global variables, so concurrent use gives wrong results or crashes unless calls to it are serialized.
  !!}
  use :: Display          , only : displayVerbositySet, verbosityLevelStandard
  use :: Events_Hooks     , only : eventsHooksInitialize
  use :: Kind_Numbers     , only : kind_int8
  use :: Nearest_Neighbors, only : nearestNeighbors
  use :: OMP_Lib          , only : OMP_Set_Num_Threads
  use :: Unit_Tests       , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  integer                  , parameter                               :: countPoints      =4000, countQueries=4000, &
       &                                                                countNeighbors   =   8, countThreads=   8
  double precision         , parameter                               :: radius           =0.1d0
  double precision         , dimension(3,countPoints   )             :: points
  double precision         , dimension(3,countQueries  )             :: queries
  integer                  , dimension(countNeighbors,countQueries)  :: indexSerial       , indexShared       , &
       &                                                                indexPrivate
  double precision         , dimension(countNeighbors,countQueries)  :: distanceSerial    , distanceShared    , &
       &                                                                distancePrivate
  integer                  , dimension(               countQueries)  :: countSerial       , countShared       , &
       &                                                                countPrivate
  integer                  , dimension(:)            , allocatable   :: neighborIndex
  double precision         , dimension(:)            , allocatable   :: neighborDistance
  type            (nearestNeighbors)                                 :: neighborFinder
  type            (nearestNeighbors)                 , allocatable   :: neighborFinderPrivate
  integer                                                            :: i                 , j
  integer         (kind_int8      )                                  :: seed              =12345_kind_int8

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Initialize event hooks, which are used in OpenMP parallel regions.
  call eventsHooksInitialize()
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Nearest neighbors")
  ! Generate points and queries uniformly in the unit cube using a simple linear congruential generator, so that the test is
  ! deterministic.
  do i=1,countPoints
     do j=1,3
        points (j,i)=uniformDeviate()
     end do
  end do
  do i=1,countQueries
     do j=1,3
        queries(j,i)=uniformDeviate()
     end do
  end do
  ! Serial searches.
  neighborFinder=nearestNeighbors(transpose(points))
  do i=1,countQueries
     call neighborFinder%search           (queries(:,i),countNeighbors,0.0d0,indexSerial(:,i),distanceSerial(:,i))
     call neighborFinder%searchFixedRadius(queries(:,i),radius        ,0.0d0,countSerial(  i)                    )
  end do
  ! Check the serial searches against a brute-force search for the first few queries.
  call Assert('serial search matches brute force',all([(bruteForceNearest(queries(:,i)) == indexSerial(1,i),i=1,100)]),.true.)
  call Assert('serial fixed radius count matches brute force',all([(bruteForceCount(queries(:,i)) == countSerial(i),i=1,100)]),.true.)
  ! Concurrent searches of a shared tree.
  call OMP_Set_Num_Threads(countThreads)
  !$omp parallel do private(neighborIndex,neighborDistance) schedule(dynamic,1)
  do i=1,countQueries
     call neighborFinder%search           (queries(:,i),countNeighbors,0.0d0,indexShared(:,i),distanceShared(:,i)                )
     call neighborFinder%searchFixedRadius(queries(:,i),radius        ,0.0d0,countShared(  i),neighborIndex      ,neighborDistance)
  end do
  !$omp end parallel do
  call Assert('concurrent search of a shared tree: indices',all(indexShared == indexSerial),.true.)
  call Assert('concurrent search of a shared tree: distances',all(distanceShared == distanceSerial),.true.)
  call Assert('concurrent search of a shared tree: fixed radius counts',all(countShared == countSerial),.true.)
  ! Concurrent searches of trees built by each thread.
  !$omp parallel private(neighborFinderPrivate,neighborIndex,neighborDistance)
  allocate(neighborFinderPrivate)
  neighborFinderPrivate=nearestNeighbors(transpose(points))
  !$omp do schedule(dynamic,1)
  do i=1,countQueries
     call neighborFinderPrivate%search           (queries(:,i),countNeighbors,0.0d0,indexPrivate(:,i),distancePrivate(:,i)                )
     call neighborFinderPrivate%searchFixedRadius(queries(:,i),radius        ,0.0d0,countPrivate(  i),neighborIndex       ,neighborDistance)
  end do
  !$omp end do
  deallocate(neighborFinderPrivate)
  !$omp end parallel
  call Assert('concurrent search of per-thread trees: indices',all(indexPrivate == indexSerial),.true.)
  call Assert('concurrent search of per-thread trees: distances',all(distancePrivate == distanceSerial),.true.)
  call Assert('concurrent search of per-thread trees: fixed radius counts',all(countPrivate == countSerial),.true.)
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

contains

  double precision function uniformDeviate()
    !!{RST
    Return a uniform deviate in :math:`[0,1)` from a linear congruential generator.
    !!}
    implicit none

    seed          =modulo(16807_kind_int8*seed,2147483647_kind_int8)
    uniformDeviate=dble(seed)/2147483647.0d0
    return
  end function uniformDeviate

  integer function bruteForceNearest(query)
    !!{RST
    Return the index of the point nearest to the given query by brute force.
    !!}
    implicit none
    double precision, dimension(3), intent(in   ) :: query

    bruteForceNearest=minloc(sum((points-spread(query,2,countPoints))**2,dim=1),dim=1)
    return
  end function bruteForceNearest

  integer function bruteForceCount(query)
    !!{RST
    Return the number of points within the search radius of the given query by brute force.
    !!}
    implicit none
    double precision, dimension(3), intent(in   ) :: query

    bruteForceCount=count(sum((points-spread(query,2,countPoints))**2,dim=1) <= radius**2)
    return
  end function bruteForceCount

end program Test_Nearest_Neighbors
