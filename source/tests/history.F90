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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a program to test interpolation in history objects.
!!}

program Test_History
  !!{RST
  Tests interpolation in ``history`` objects.

  Two operations formerly constructed a full ``interpolator`` object - copying the time array and allocating GSL objects - in
  order to obtain a bracketing index (in ``stellarPopulationPropertiesNoninstantaneous``) or a pair of linear interpolation
  weights (in ``history%interpolatedIncrement``). Both now compute what they need directly. Since these values feed the ODE
  right-hand side, the replacements must be *bit-identical* to the ``interpolator``-based results, not merely equal to within
  some tolerance. Each test below therefore compares against a reference computed via the ``interpolator`` object with no
  tolerance supplied to ``Assert``, which requires exact equality.
  !!}
  use, intrinsic :: ISO_C_Binding          , only : c_size_t
  use            :: Arrays_Search          , only : searchArray
  use            :: Display                , only : displayVerbositySet, verbosityLevelStandard
  use            :: Histories              , only : history
  use            :: Numerical_Interpolation, only : interpolator
  use            :: Unit_Tests             , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, &
  &                                                 Unit_Tests_Finish
  implicit none
  ! A non-uniform time grid, of the kind a stellar population property history carries.
  double precision              , dimension(8)  :: timeHistory            =[0.1d0,0.2d0,0.4d0,0.8d0,1.6d0,3.2d0,6.4d0,13.8d0]
  ! Times at which to query the grid: interior points, points coincident with grid points, and points below and above the
  ! tabulated range.
  double precision                , dimension(13) :: timeQuery    =[                                                            &
       &                                                            -1.0d0, 0.0d0, 0.1d0, 0.15d0, 0.2d0, 0.7d0, 1.6d0, 1.61d0,  &
       &                                                             5.0d0, 6.4d0,13.8d0,13.81d0,20.0d0                         &
       &                                                           ]
  integer         (c_size_t    ), dimension(13) :: indexSearch                                                               , &
       &                                           indexInterpolatorFresh                                                    , &
       &                                           indexInterpolatorCached
  type            (interpolator)                :: interpolator_
  type            (history     )                :: history_                                                                  , &
       &                                           historyReference                                                          , &
       &                                           historyAdd                                                                , &
       &                                           historyOriginal
  integer                                       :: i                                                                         , &
       &                                           j
  double precision                              :: timeBase                                                                  , &
       &                                           timeAdd

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("history interpolation")

  ! Locating the bracketing index in a history time array. `stellarPopulationPropertiesNoninstantaneous` previously built an
  ! `interpolator` from the history times purely to call its `locate()` method once. Verify that `searchArray` returns exactly
  ! the same index, both against a freshly-constructed interpolator (which is what the old code built on every call) and
  ! against a persistent one, whose cached bracket must not alter the result.
  do i=1,size(timeQuery)
     indexSearch            (i)=searchArray  (timeHistory,timeQuery(i))
     interpolator_             =interpolator (timeHistory             )
     indexInterpolatorFresh (i)=interpolator_%locate     (timeQuery(i))
  end do
  interpolator_=interpolator(timeHistory)
  do i=1,size(timeQuery)
     indexInterpolatorCached(i)=interpolator_%locate     (timeQuery(i))
  end do
  call Assert('bracketing index from searchArray matches a freshly-built interpolator',int(indexSearch),int(indexInterpolatorFresh ))
  call Assert('bracketing index from searchArray matches a cached interpolator'       ,int(indexSearch),int(indexInterpolatorCached))

  ! Interpolated increment of one history onto the times of another. Build a base history and a history to add, on different
  ! time grids, with the times of the base history deliberately spanning the full range of the added history plus points
  ! outside it (which must be left untouched). Increment a copy using the reference (`interpolator`-based) algorithm, and
  ! compare with the result of `interpolatedIncrement`.
  call history_%create(3,11)
  do i=1,size(history_%time)
     timeBase          =0.05d0+0.2d0*dble(i-1)
     history_%time  (i)=timeBase
     do j=1,size(history_%data,dim=2)
        history_%data(i,j)=dble(j)*sin(3.0d0*timeBase)+0.5d0*dble(i)
     end do
  end do
  call historyAdd%create(3,6)
  do i=1,size(historyAdd%time)
     timeAdd             =0.3d0+0.31d0*dble(i-1)
     historyAdd%time  (i)=timeAdd
     do j=1,size(historyAdd%data,dim=2)
        historyAdd%data(i,j)=dble(j)*cos(2.0d0*timeAdd)-0.25d0*dble(i)
     end do
  end do
  historyOriginal =history_
  historyReference=history_
  call Interpolated_Increment_Reference(historyReference,historyAdd)
  call history_%interpolatedIncrement  (historyAdd      )
  call Assert('interpolated increment is bit-identical to the interpolator-based reference',history_%data,historyReference%data)
  ! Confirm that the comparison above is not vacuous: every base-history time lying within the range of the added history
  ! (points 3-10, given the grids constructed above) must have been changed by a non-zero amount, while those outside that
  ! range (points 1, 2 and 11) must have been left exactly as they were.
  call Assert('interpolated increment changes every point within the range of the added history' ,all(history_%data( 3:10,:) /= historyOriginal%data( 3:10,:)),.true.)
  call Assert('interpolated increment leaves points outside the range of the added history alone',history_%data( 1: 2,:)     ,historyOriginal%data( 1: 2,:)                )
  call Assert('interpolated increment leaves points beyond the range of the added history alone' ,history_%data(11:11,:)     ,historyOriginal%data(11:11,:)                )

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

contains

  subroutine Interpolated_Increment_Reference(history_,addHistory)
    !!{RST
    Reference implementation of ``history%interpolatedIncrement``, retaining the ``interpolator`` object that the production
    implementation no longer builds. Used only to check that the production implementation gives bit-identical results.
    !!}
    implicit none
    type            (history     ), intent(inout) :: history_
    type            (history     ), intent(in   ) :: addHistory
    double precision              , dimension(2)  :: interpolationFactors
    integer                                       :: iPoint              , iHistory
    integer         (c_size_t    )                :: interpolationPoint  , addHistoryPointCount
    type            (interpolator)                :: interpolator_

    addHistoryPointCount=size(addHistory%time)
    interpolationPoint  =1
    interpolator_       =interpolator(addHistory%time)
    do iPoint=1,size(history_%time)
       if (history_%time(iPoint) >= addHistory%time(1) .and. history_%time(iPoint) <= addHistory%time(addHistoryPointCount)) then
          do while (history_%time(iPoint) > addHistory%time(interpolationPoint) .and. interpolationPoint < addHistoryPointCount-1)
             interpolationPoint=interpolationPoint+1
          end do
          call interpolator_%linearWeights(history_%time(iPoint),interpolationPoint,interpolationFactors)
          forall(iHistory=1:size(history_%data,dim=2))
             history_%data (iPoint,iHistory)=history_%data (iPoint,iHistory)+addHistory%data(interpolationPoint,iHistory)&
                  &*interpolationFactors(1)+addHistory%data(interpolationPoint+1,iHistory)*interpolationFactors(2)
          end forall
       end if
    end do
    return
  end subroutine Interpolated_Increment_Reference

end program Test_History
