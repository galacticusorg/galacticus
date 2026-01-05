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
Contains a program to test tables.
!!}

program Test_Tables
  !!{
  Tests that tables work correctly.
  !!}
  use :: Array_Utilities, only : directionDecreasing         , directionIncreasing
  use :: Display        , only : displayVerbositySet         , verbosityLevelStandard
  use :: Tables         , only : table                       , table1D                 , table1DLinearCSpline              , table1DLinearLinear, &
          &                      table1DLinearMonotoneCSpline, table1DLogarithmicLinear, table1DNonUniformLinearLogarithmic, table2DLogLogLin
  use :: Unit_Tests     , only : Assert                      , Unit_Tests_Begin_Group  , Unit_Tests_End_Group              , Unit_Tests_Finish
  implicit none
  class           (table           ), allocatable :: myTable
  class           (table1D         ), allocatable :: myReversedTable
  type            (table2DLogLogLin)              :: myTable2D
  integer                                         :: i              , j
  double precision                                :: x              , y, &
       &                                             yPrevious
  logical                                         :: isMonotonic

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Tables")

  ! Allocate a table object.
  allocate(table1DLinearLinear :: myTable)
  select type (myTable)
  type is (table1DLinearLinear)
     ! Create a table.
     call myTable%create  (1.0d0,6.0d0,11                                                          )
     call myTable%populate([2.0d0,3.0d0,-23.0d0,4.0d0,6.0d0,-1.0d0,-5.0d0,-0.1d0,5.0d0,9.0d0,3.0d0])
     ! Test interpolation in 1-D table.
     call Assert(                                           &
          &      'linear interpolation in 1D linear table', &
          &      [                                          &
          &       myTable%interpolate(1.5d0)              , &
          &       myTable%interpolate(2.4d0)              , &
          &       myTable%interpolate(4.1d0)              , &
          &       myTable%interpolate(5.7d0)                &
          &      ]                                        , &
          &      [                                          &
          &        3.00d0                                 , &
          &       -1.40d0                                 , &
          &       -4.02d0                                 , &
          &        6.60d0                                   &
          &      ]                                        , &
          &      absTol=1.0d-6                              &
          &     )
     ! Test interpolation gradient in 1-D table.
     call Assert(                                           &
          &      'linear gradient in 1D linear table'     , &
          &      [                                          &
          &       myTable%interpolateGradient(1.5d0)      , &
          &       myTable%interpolateGradient(2.4d0)      , &
          &       myTable%interpolateGradient(4.1d0)      , &
          &       myTable%interpolateGradient(5.7d0)        &
          &      ]                                        , &
          &      [                                          &
          &       -52.0d0                                 , &
          &        54.0d0                                 , &
          &         9.8d0                                 , &
          &       -12.0d0                                   &
          &      ]                                        , &
          &      absTol=1.0d-6                              &
          &     )
     ! Destroy the table.
     call myTable%destroy()
  end select
  deallocate(myTable)

  ! Allocate a table object.
  allocate(table1DLogarithmicLinear :: myTable)
  select type (myTable)
  type is (table1DLogarithmicLinear)
     ! Create a table.
     call myTable%create  (1.0d1,1.0d6,11                                                          )
     call myTable%populate([2.0d0,3.0d0,-23.0d0,4.0d0,6.0d0,-1.0d0,-5.0d0,-0.1d0,5.0d0,9.0d0,3.0d0])
     ! Test interpolation in 1-D table.
     call Assert(                                                &
          &      'linear interpolation in 1D logarithmic table', &
          &      [                                               &
          &       myTable%interpolate(10.0d0**1.5d0)           , &
          &       myTable%interpolate(10.0d0**2.4d0)           , &
          &       myTable%interpolate(10.0d0**4.1d0)           , &
          &       myTable%interpolate(10.0d0**5.7d0)             &
          &      ]                                             , &
          &      [                                               &
          &        3.00d0                                      , &
          &       -1.40d0                                      , &
          &       -4.02d0                                      , &
          &        6.60d0                                        &
          &      ]                                             , &
          &      absTol=1.0d-6                                   &
          &     )
     ! Test interpolation in 1-D table.
     call Assert(                                                &
          &      'linear gradient in 1D logarithmic table'     , &
          &      [                                               &
          &       myTable%interpolateGradient(10.0d0**1.5d0)   , &
          &       myTable%interpolateGradient(10.0d0**2.4d0)   , &
          &       myTable%interpolateGradient(10.0d0**4.1d0)   , &
          &       myTable%interpolateGradient(10.0d0**5.7d0)     &
          &      ]                                             , &
          &      [                                               &
          &         2.0d0/10.0d0**1.5d0/log(10.0d0)            , &
          &        54.0d0/10.0d0**2.4d0/log(10.0d0)            , &
          &         9.8d0/10.0d0**4.1d0/log(10.0d0)            , &
          &       -12.0d0/10.0d0**5.7d0/log(10.0d0)              &
          &      ]                                             , &
          &      absTol=1.0d-6                                   &
          &     )
     ! Destroy the table.
     call myTable%destroy()
  end select
  deallocate(myTable)

  ! Allocate a table object.
  allocate(table1DLogarithmicLinear :: myTable)
  select type (myTable)
  type is (table1DLogarithmicLinear)
     ! Create a table.
     call myTable%create  (1.0d1,1.0d6,11                                                        )
     call myTable%populate([2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0,8.0d0,9.0d0,10.0d0,11.0d0,12.0d0])
     ! Reverse the table.
     call myTable%reverse(myReversedTable,precise=.true.)
     call Assert(                                                            &
          &      'reversed log-linear table consistent with original table', &
          &      [                                                           &
          &       myReversedTable%interpolate(+ 3.00d0)                    , &
          &       myReversedTable%interpolate(+ 4.80d0)                    , &
          &       myReversedTable%interpolate(+ 8.20d0)                    , &
          &       myReversedTable%interpolate(+11.40d0)                      &
          &      ],                                                          &
          &      [                                                           &
          &       10.0d0**1.5d0                                            , &
          &       10.0d0**2.4d0                                            , &
          &       10.0d0**4.1d0                                            , &
          &       10.0d0**5.7d0                                              &
          &      ],                                                          &
          &      absTol=1.0d-6                                               &
          &     )
     ! Destroy the table.
     call myTable        %destroy()
     call myReversedTable%destroy()
  end select
  deallocate(myTable        )
  deallocate(myReversedTable)

  ! Allocate a table object.
  allocate(table1DLinearLinear :: myTable)
  select type (myTable)
  type is (table1DLinearLinear)
     ! Create a table.
     call myTable%create  (1.0d0,6.0d0,11                                                        )
     call myTable%populate([2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0,8.0d0,9.0d0,10.0d0,11.0d0,12.0d0])
     ! Assert that the table is monotonic.
     call Assert('table is monotonic'                   ,myTable%isMonotonic(                             ),.true. )
     ! Assert that the table is monotonically increasing.
     call Assert('table is monotonically increasing'    ,myTable%isMonotonic(direction=directionIncreasing),.true. )
     ! Assert that the table is not monotonically decreasing.
     call Assert('table is not monotonically decreasing',myTable%isMonotonic(direction=directionDecreasing),.false.)
     ! Reverse the table.
     call myTable%reverse(myReversedTable)
     call Assert(                                                &
          &      'reverse table consistent with original table', &
          &      [                                               &
          &       myReversedTable%interpolate( 2.0d0),           &
          &       myReversedTable%interpolate( 3.0d0),           &
          &       myReversedTable%interpolate( 4.0d0),           &
          &       myReversedTable%interpolate( 5.0d0),           &
          &       myReversedTable%interpolate( 6.0d0),           &
          &       myReversedTable%interpolate( 7.0d0),           &
          &       myReversedTable%interpolate( 8.0d0),           &
          &       myReversedTable%interpolate( 9.0d0),           &
          &       myReversedTable%interpolate(10.0d0),           &
          &       myReversedTable%interpolate(11.0d0),           &
          &       myReversedTable%interpolate(12.0d0)            &
          &      ],                                              &
          &      [                                               &
          &                                    1.0d0 ,           &
          &                                    1.5d0 ,           &
          &                                    2.0d0 ,           &
          &                                    2.5d0 ,           &
          &                                    3.0d0 ,           &
          &                                    3.5d0 ,           &
          &                                    4.0d0 ,           &
          &                                    4.5d0 ,           &
          &                                    5.0d0 ,           &
          &                                    5.5d0 ,           &
          &                                    6.0d0             &
          &      ],                                              &
          &      absTol=1.0d-6                                   &
          &     )
     ! Destroy the tables.
     call myTable        %destroy()
     call myReversedTable%destroy()
  end select
  deallocate(myReversedTable)
  deallocate(myTable        )

  ! Allocate a table object.
  allocate(table1DLinearCSpline :: myTable)
  select type (myTable)
  type is (table1DLinearCSpline)
     ! Create a table.
     call myTable%create  (1.0d0,6.0d0,11                                                          )
     call myTable%populate([2.0d0,3.0d0,-23.0d0,4.0d0,6.0d0,-1.0d0,-5.0d0,-0.1d0,5.0d0,9.0d0,3.0d0])
     ! Test interpolation in 1-D table.
     call Assert(                                                 &
          &      'cubic spline interpolation in 1D linear table', &
          &      [                                                &
          &       myTable%interpolate(1.5d0)                    , &
          &       myTable%interpolate(2.4d0)                    , &
          &       myTable%interpolate(4.1d0)                    , &
          &       myTable%interpolate(5.7d0)                      &
          &      ]                                              , &
          &      [                                                &
          &        3.0000000d0                                  , &
          &       -1.8298492d0                                  , &
          &       -4.5767772d0                                  , &
          &        7.6134254d0                                    &
          &      ]                                              , &
          &      absTol=1.0d-6                                    &
          &     )
     ! Test gradient interpolation in 1-D table.
     call Assert(                                                 &
          &      'cubic spline gradient in 1D linear table'     , &
          &      [                                                &
          &       myTable%interpolateGradient(1.5d0)            , &
          &       myTable%interpolateGradient(2.4d0)            , &
          &       myTable%interpolateGradient(4.1d0)            , &
          &       myTable%interpolateGradient(5.7d0)              &
          &      ]                                              , &
          &      [                                                &
          &       -43.8943139d0                                 , &
          &        66.8232266d0                                 , &
          &         6.6091756d0                                 , &
          &       -11.5777394d0                                   &
          &      ]                                              , &
          &      absTol=1.0d-6                                    &
          &     )
  end select
  deallocate(myTable)

  ! Allocate a monotonic cubic spline interpolator table.
  allocate(table1DLinearMonotoneCSpline :: myTable)
  select type (myTable)
  type is (table1DLinearMonotoneCSpline)
     ! Create a table.
     call myTable%create  (1.0d0,6.0d0,11                                                       )
     call myTable%populate([2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0,8.0d0,9.0d0,9.0d0,11.0d0,12.0d0])
     ! Test interpolation in 1-D table.
     isMonotonic=.true.
     yPrevious  =0.0d0
     do i=1,10000
        x=1.0d0+5.0d0*dble(i-1)/dble(9999)
        y=myTable%interpolate(x)
        if (i > 1 .and. y < yPrevious) isMonotonic=.false.
        yPrevious=y
     end do
     call Assert('monotone cubic spline interpolator is monotonic',isMonotonic,.true.)
  end select
  deallocate(myTable)

  ! Allocate a table object.
  allocate(table1DNonUniformLinearLogarithmic :: myTable)
  select type (myTable)
  type is (table1DNonUniformLinearLogarithmic)
     ! Create a table.
     call myTable%create  ([1.0d0,1.5d0, 2.0d0,2.5d0,3.0d0,3.5d0,4.0d0,4.5d0,5.0d0,5.5d0,6.0d0])
     call myTable%populate([2.0d0,3.0d0,23.0d0,4.0d0,6.0d0,1.0d0,5.0d0,0.1d0,5.0d0,9.0d0,3.0d0])
     ! Test interpolation in 1-D table.
     call Assert(                                                 &
          &      'logarithmic interpolation in 1D linear table' , &
          &      [                                                &
          &       myTable%interpolate(1.5d0)                    , &
          &       myTable%interpolate(2.4d0)                    , &
          &       myTable%interpolate(4.1d0)                    , &
          &       myTable%interpolate(5.7d0)                      &
          &      ]                                              , &
          &      [                                                &
          &       +3.000000000d0                                , &
          &       +5.675361899d0                                , &
          &       +2.286525260d0                                , &
          &       +5.799546135d0                                  &
          &      ]                                              , &
          &      absTol=1.0d-6                                    &
          &     )

     ! Test gradient interpolation in 1-D table.
     call Assert(                                                 &
          &      'logarithmic gradient in 1D linear table'      , &
          &      [                                                &
          &       myTable%interpolateGradient(1.5d0)            , &
          &       myTable%interpolateGradient(2.4d0)            , &
          &       myTable%interpolateGradient(4.1d0)            , &
          &       myTable%interpolateGradient(5.7d0)              &
          &      ]                                              , &
          &      [                                                &
          &       +12.22129156d0                                , &
          &       -19.85468441d0                                , &
          &       -17.88987883d0                                , &
          &       -12.74290530d0                                  &
          &      ]                                              , &
          &      absTol=1.0d-6                                    &
          &     )
  end select
  deallocate(myTable)

  ! Create a 2D log-log-linear table.
  call myTable2D%create  (1.0d0,1.0d5,11,1.0d0,1.0d5,11,1)
  do i=1,11
     do j=1,11
        call myTable2D%populate(dble(i+j),i,j)
     end do
  end do
  ! Test interpolation in 2D table.
  call Assert(                                                     &
       &      'linear interpolation in 2D log-log table'         , &
       &      [                                                    &
       &       myTable2D%interpolate(10.0d0**1.5d0,10.0d0**2.4d0), &
       &       myTable2D%interpolate(10.0d0**2.4d0,10.0d0**3.1d0), &
       &       myTable2D%interpolate(10.0d0**3.1d0,10.0d0**4.7d0), &
       &       myTable2D%interpolate(10.0d0**4.7d0,10.0d0**4.9d0)  &
       &      ]                                                  , &
       &      [                                                    &
       &        9.8d0                                            , &
       &       13.0d0                                            , &
       &       17.6d0                                            , &
       &       21.2d0                                              &
       &      ]                                                  , &
       &      absTol=1.0d-6                                        &
       &     )
  ! Test gradient interpolation in 2D table.
  call Assert(                                                               &
       &      'linear gradient interpolation in 2D log-log table'          , &
       &      [                                                              &
       &       myTable2D%interpolateGradient(10.0d0**1.5d0,10.0d0**2.4d0,1), &
       &       myTable2D%interpolateGradient(10.0d0**2.4d0,10.0d0**3.1d0,1), &
       &       myTable2D%interpolateGradient(10.0d0**3.1d0,10.0d0**4.7d0,2), &
       &       myTable2D%interpolateGradient(10.0d0**4.7d0,10.0d0**4.9d0,2)  &
       &      ]                                                            , &
       &      [                                                              &
       &       2.0d0*0.1d0**1.5d0/log(10.0d0)                              , &
       &       2.0d0*0.1d0**2.4d0/log(10.0d0)                              , &
       &       2.0d0*0.1d0**4.7d0/log(10.0d0)                              , &
       &       2.0d0*0.1d0**4.9d0/log(10.0d0)                                &
       &      ]                                                            , &
       &      absTol=1.0d-6                                                  &
       &     )
  ! Destroy the table.
  call myTable2D%destroy()

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

end program Test_Tables
