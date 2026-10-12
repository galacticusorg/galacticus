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
Contains a program to test tables.
!!}

program Test_Tables
  !!{RST
  Tests that tables work correctly.
  !!}
  use :: Array_Utilities , only : directionDecreasing             , directionIncreasing
  use :: Display         , only : displayVerbositySet             , verbosityLevelStandard
  use :: Numerical_Ranges, only : Range_Pinned                    , rangeLattice                     , gridSchemePerOctave               , gridSchemePerDecade, &
          &                       gridSchemePerUnit               , Make_Range                       , rangeTypeLinear
  use :: Tables          , only : table                           , table1D                          , table1DNonUniformLinearLogarithmic, table1DLinearLinear, &
          &                       table1DLinearMonotoneCSpline    , table1DLogarithmicLinear         , table1DLinearCSpline              , table2DLogLogLin   , &
          &                       table1DLogarithmicCSpline       , table1DLogarithmicMonotoneCSpline, table2D                           , table2DLinLinLin   , &
          &                       table1DGeneric                  , table1DMonotoneCSpline
  use :: Table_Labels    , only : enumerationExtrapolationTypeType, extrapolationTypeFix             , extrapolationTypeZero
  use :: Unit_Tests      , only : Assert                          , Unit_Tests_Begin_Group           , Unit_Tests_End_Group              , Unit_Tests_Finish
  implicit none
  class           (table                           ), allocatable                 :: myTable
  class           (table1D                         ), allocatable                 :: myReversedTable
  type            (table2DLogLogLin                )                              :: myTable2D
  type            (rangeLattice                    )                              :: latticeNarrow        , latticeWide
  integer                                                                         :: i                    , j                   , &
       &                                                                             offset
  double precision                                                                :: x                    , y                   , &
       &                                                                             yPrevious
  logical                                                                         :: isMonotonic
  logical                                           , allocatable, dimension(:  ) :: isComputed
  double precision                                  , allocatable, dimension(:  ) :: xValuesNarrow        , xValuesWide         , &
       &                                                                             xValuesDirect
  double precision                                  , allocatable, dimension(:,:) :: yValuesNarrow        , yValuesWide         , &
       &                                                                             yValuesDirect
  type            (table2DLogLogLin                )                              :: myTable2DExtend
  type            (rangeLattice                    )                              :: latticeX2D           , latticeY2D
  logical                                           , allocatable, dimension(:,:) :: isComputed2D
  double precision                                  , allocatable, dimension(:,:) :: zValuesNarrow2D      , zValuesWide2D
  double precision                                  , allocatable, dimension(:  ) :: xValuesSpline        , interpolatedExtended, &
       &                                                                             interpolatedDirect
  double precision                                  , allocatable, dimension(:,:) :: yValuesSpline
  class           (table2D                         ), allocatable                 :: myTable2DGeneric
  class           (table1D                         ), allocatable                 :: myRangeTable
  type            (enumerationExtrapolationTypeType)              :: extrapolationRange
  character       (len=32                          )                              :: labelExtrapolation
  logical                                                                         :: logarithmic
  type            (table2DLinLinLin                )                              :: myTable2DLinearExtend
  double precision                                  , allocatable, dimension(:,:) :: zValuesNarrowLinear2D, zValuesWideLinear2D

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

  call Unit_Tests_End_Group()

  ! Exercise both concrete 2D table types through the `table2D` base class. Every assertion made by `Test_Table_2D` is written
  ! purely in terms of that base class - which is itself the demonstration that both types now provide the same interface.
  call Unit_Tests_Begin_Group("2D table interface")
  allocate(table2DLinLinLin :: myTable2DGeneric)
  select type (myTable2DGeneric)
  type is (table2DLinLinLin)
     call myTable2DGeneric%create(                                                      &
          &                       Make_Range(1.0d0,6.0d0,11,rangeType=rangeTypeLinear), &
          &                       Make_Range(2.0d0,9.0d0, 8,rangeType=rangeTypeLinear), &
          &                       tableCount=2                                          &
          &                      )
  end select
  call Test_Table_2D('linLinLin',myTable2DGeneric)
  deallocate(myTable2DGeneric)
  allocate(table2DLogLogLin :: myTable2DGeneric)
  select type (myTable2DGeneric)
  type is (table2DLogLogLin)
     call myTable2DGeneric%create(1.0d0,1.0d5,11,2.0d0,1.0d3,8,tableCount=2)
  end select
  call Test_Table_2D('logLogLin',myTable2DGeneric)
  deallocate(myTable2DGeneric)
  call Unit_Tests_End_Group()

  ! Test values and gradients of 1D tables beyond their range. With "fix" extrapolation the value beyond the range is that at the
  ! boundary, and so the gradient there is zero. With "zero" extrapolation (supported only by the linear and logarithmic tables)
  ! both are zero.
  call Unit_Tests_Begin_Group("Table range")
  do i=1,10
     extrapolationRange=extrapolationTypeFix
     logarithmic       =.false.
     select case (i)
     case ( 1)
        allocate(table1DGeneric                    :: myRangeTable)
     case ( 2)
        allocate(table1DLinearLinear               :: myRangeTable)
     case ( 3)
        allocate(table1DLinearLinear               :: myRangeTable)
        extrapolationRange=extrapolationTypeZero
     case ( 4)
        allocate(table1DLinearCSpline              :: myRangeTable)
     case ( 5)
        allocate(table1DMonotoneCSpline            :: myRangeTable)
     case ( 6)
        allocate(table1DLinearMonotoneCSpline      :: myRangeTable)
     case ( 7)
        allocate(table1DLogarithmicLinear          :: myRangeTable)
        logarithmic       =.true.
     case ( 8)
        allocate(table1DLogarithmicLinear          :: myRangeTable)
        logarithmic       =.true.
        extrapolationRange=extrapolationTypeZero
     case ( 9)
        allocate(table1DLogarithmicCSpline         :: myRangeTable)
        logarithmic       =.true.
     case (10)
        allocate(table1DLogarithmicMonotoneCSpline :: myRangeTable)
        logarithmic       =.true.
     end select
     ! Linear tables span 1 ≤ x ≤ 7, logarithmic tables span 1 ≤ x ≤ 1000, each with seven points.
     select type (myRangeTable)
     type is (table1DGeneric                   )
        call myRangeTable%create([1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0],extrapolationType=[extrapolationRange,extrapolationRange])
        call myRangeTable%populate([1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0])
     type is (table1DLinearLinear              )
        call myRangeTable%create(1.0d0,7.0d0,7,extrapolationType=[extrapolationRange,extrapolationRange])
        call myRangeTable%populate([1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0])
     type is (table1DLinearCSpline             )
        call myRangeTable%create(1.0d0,7.0d0,7,extrapolationType=[extrapolationRange,extrapolationRange])
        call myRangeTable%populate([1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0])
     type is (table1DMonotoneCSpline           )
        call myRangeTable%create([1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0],extrapolationType=[extrapolationRange,extrapolationRange])
        call myRangeTable%populate([1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0])
     type is (table1DLinearMonotoneCSpline     )
        call myRangeTable%create(1.0d0,7.0d0,7,extrapolationType=[extrapolationRange,extrapolationRange])
        call myRangeTable%populate([1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0])
     type is (table1DLogarithmicLinear         )
        call myRangeTable%create(1.0d0,1.0d3,7,extrapolationType=[extrapolationRange,extrapolationRange])
        call myRangeTable%populate([1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0])
     type is (table1DLogarithmicCSpline        )
        call myRangeTable%create(1.0d0,1.0d3,7,extrapolationType=[extrapolationRange,extrapolationRange])
        call myRangeTable%populate([1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0])
     type is (table1DLogarithmicMonotoneCSpline)
        call myRangeTable%create(1.0d0,1.0d3,7,extrapolationType=[extrapolationRange,extrapolationRange])
        call myRangeTable%populate([1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0])
     end select
     if (extrapolationRange == extrapolationTypeFix) then
        labelExtrapolation='fix'
     else
        labelExtrapolation='zero'
     end if
     call Test_Table_Range(trim(Table_Type_Label(i))//' ('//trim(labelExtrapolation)//')',myRangeTable,logarithmic,extrapolationRange == extrapolationTypeZero)
     call myRangeTable%destroy()
     deallocate(myRangeTable)
  end do
  call Unit_Tests_End_Group()

  ! Test values and gradients of 2D tables beyond their range, as for the 1D tables above.
  call Unit_Tests_Begin_Group("2D table range")
  do i=1,4
     if (i == 1 .or. i == 3) then
        extrapolationRange=extrapolationTypeFix
        labelExtrapolation='fix'
     else
        extrapolationRange=extrapolationTypeZero
        labelExtrapolation='zero'
     end if
     logarithmic=i <= 2
     if (logarithmic) then
        allocate(table2DLogLogLin :: myTable2DGeneric)
     else
        allocate(table2DLinLinLin :: myTable2DGeneric)
     end if
     select type (myTable2DGeneric)
     type is (table2DLogLogLin)
        call myTable2DGeneric%create(1.0d0,1.0d3,7,1.0d0,1.0d3,7,extrapolationTypeX=extrapolationRange,extrapolationTypeY=extrapolationRange)
        do j=1,7
           do offset=1,7
              call myTable2DGeneric%populate(log(myTable2DGeneric%x(j))+2.0d0*log(myTable2DGeneric%y(offset)),j,offset)
           end do
        end do
     type is (table2DLinLinLin)
        call myTable2DGeneric%create(Make_Range(1.0d0,7.0d0,7,rangeType=rangeTypeLinear),Make_Range(1.0d0,7.0d0,7,rangeType=rangeTypeLinear),extrapolationTypeX=extrapolationRange,extrapolationTypeY=extrapolationRange)
        do j=1,7
           do offset=1,7
              call myTable2DGeneric%populate(myTable2DGeneric%x(j)+2.0d0*myTable2DGeneric%y(offset),j,offset)
           end do
        end do
     end select
     if (logarithmic) then
        call Test_Table_2D_Range('logLogLin ('//trim(labelExtrapolation)//')',myTable2DGeneric,logarithmic,extrapolationRange == extrapolationTypeZero)
     else
        call Test_Table_2D_Range('linLinLin ('//trim(labelExtrapolation)//')',myTable2DGeneric,logarithmic,extrapolationRange == extrapolationTypeZero)
     end if
     call myTable2DGeneric%destroy()
     deallocate(myTable2DGeneric)
  end do
  call Unit_Tests_End_Group()

  ! Test extension of tables onto an absolute lattice.
  call Unit_Tests_Begin_Group("Table extension")

  ! Build a logarithmic table on a per-octave lattice, tabulating y=x².
  allocate(table1DLogarithmicLinear :: myTable)
  select type (myTable)
  type is (table1DLogarithmicLinear)
     latticeNarrow=Range_Pinned(3.0d0,4,gridSchemePerOctave)
     call myTable%extend(latticeNarrow,isComputed)
     call Assert('extension of an empty table requires every point to be computed',count(isComputed),0                  )
     call Assert('extension of an empty table gives the lattice point count'      ,myTable%size()   ,latticeNarrow%count)
     do i=1,myTable%size()
        call myTable%populate(myTable%x(i)**2,i)
     end do
     xValuesNarrow=myTable%xs()
     yValuesNarrow=myTable%ys()
     ! Extend the table to a wider range on the same lattice. Only the newly added points should require computation.
     latticeWide=Range_Pinned(30.0d0,4,gridSchemePerOctave,latticeCurrent=myTable%lattice)
     call myTable%extend(latticeWide,isComputed)
     offset=latticeNarrow%indexMinimum-latticeWide%indexMinimum
     call Assert('extension marks precisely the previously computed points as computed',count(isComputed)                                   ,latticeNarrow%count)
     call Assert('extension marks the correct window as computed'                      ,all(isComputed(offset+1:offset+latticeNarrow%count)),.true.             )
     xValuesWide=myTable%xs()
     yValuesWide=myTable%ys()
     ! Both abscissae and previously computed values must be preserved bit-for-bit.
     call Assert('extension preserves abscissae bit-for-bit'                                  , &
          &      all(xValuesWide(offset+1:offset+latticeNarrow%count  ) == xValuesNarrow     ), &
          &      .true.                                                                         &
          &     )
     call Assert('extension preserves values bit-for-bit'                                     , &
          &      all(yValuesWide(offset+1:offset+latticeNarrow%count,1) == yValuesNarrow(:,1)), &
          &      .true.                                                                         &
          &     )
     ! Compute the newly added points.
     do i=1,myTable%size()
        if (.not.isComputed(i)) call myTable%populate(myTable%x(i)**2,i)
     end do
     xValuesWide=myTable%xs()
     yValuesWide=myTable%ys()
     call myTable%destroy()
  end select
  deallocate(myTable)

  ! Build a second table directly on the wider lattice - it must be bit-identical to the extended table.
  allocate(table1DLogarithmicLinear :: myTable)
  select type (myTable)
  type is (table1DLogarithmicLinear)
     call myTable%extend(latticeWide,isComputed)
     do i=1,myTable%size()
        call myTable%populate(myTable%x(i)**2,i)
     end do
     xValuesDirect=myTable%xs()
     yValuesDirect=myTable%ys()
     call Assert('a table built directly has abscissae identical to one built by extension',all(xValuesDirect      == xValuesWide     ),.true.)
     call Assert('a table built directly has values identical to one built by extension'   ,all(yValuesDirect(:,1) == yValuesWide(:,1)),.true.)
     call myTable%destroy()
  end select
  deallocate(myTable)

  ! Test extension of a two-dimensional table, in which each axis is pinned independently and the previously computed values
  ! occupy a rectangular block of the extended table.
  latticeX2D=Range_Pinned(15.0d0,4,gridSchemePerDecade,anchorEvery=2)
  latticeY2D=Range_Pinned( 3.0d0,4,gridSchemePerDecade,anchorEvery=2)
  call myTable2DExtend%extend(latticeX2D,latticeY2D,isComputed2D)
  call Assert('2D extension of an empty table requires every point to be computed',count(isComputed2D),0)
  do i=1,latticeX2D%count
     do j=1,latticeY2D%count
        call myTable2DExtend%populate(myTable2DExtend%x(i)*myTable2DExtend%y(j),i,j)
     end do
  end do
  zValuesNarrow2D=myTable2DExtend%zs()
  ! Extend both axes and check that the previously computed block is preserved exactly.
  latticeX2D=Range_Pinned(1.5d3,4,gridSchemePerDecade,anchorEvery=2,latticeCurrent=myTable2DExtend%latticeX)
  latticeY2D=Range_Pinned(3.0d2,4,gridSchemePerDecade,anchorEvery=2,latticeCurrent=myTable2DExtend%latticeY)
  call myTable2DExtend%extend(latticeX2D,latticeY2D,isComputed2D)
  call Assert('2D extension preserves precisely the previously computed block',count(isComputed2D),size(zValuesNarrow2D,dim=1)*size(zValuesNarrow2D,dim=2))
  zValuesWide2D=myTable2DExtend%zs()
  call Assert('2D extension preserves the previously computed values bit-for-bit'                               , &
       &      all(zValuesWide2D(1:size(zValuesNarrow2D,dim=1),1:size(zValuesNarrow2D,dim=2)) == zValuesNarrow2D), &
       &      .true.                                                                                              &
       &     )
  call myTable2DExtend%destroy()

  ! The same for a linearly-spaced 2D table, whose axes are pinned to `perUnit` lattices.
  latticeX2D=Range_Pinned(15.0d0,4,gridSchemePerUnit,anchorEvery=2)
  latticeY2D=Range_Pinned( 3.0d0,4,gridSchemePerUnit,anchorEvery=2)
  call myTable2DLinearExtend%extend(latticeX2D,latticeY2D,isComputed2D)
  call Assert('2D extension of an empty linear table requires every point to be computed',count(isComputed2D),0)
  do i=1,latticeX2D%count
     do j=1,latticeY2D%count
        call myTable2DLinearExtend%populate(myTable2DLinearExtend%x(i)*myTable2DLinearExtend%y(j),i,j)
     end do
  end do
  zValuesNarrowLinear2D=myTable2DLinearExtend%zs()
  latticeX2D=Range_Pinned(45.0d0,4,gridSchemePerUnit,anchorEvery=2,latticeCurrent=myTable2DLinearExtend%latticeX)
  latticeY2D=Range_Pinned(11.0d0,4,gridSchemePerUnit,anchorEvery=2,latticeCurrent=myTable2DLinearExtend%latticeY)
  call myTable2DLinearExtend%extend(latticeX2D,latticeY2D,isComputed2D)
  call Assert('2D extension of a linear table preserves precisely the previously computed block',count(isComputed2D),size(zValuesNarrowLinear2D,dim=1)*size(zValuesNarrowLinear2D,dim=2))
  zValuesWideLinear2D=myTable2DLinearExtend%zs()
  call Assert('2D extension of a linear table preserves the previously computed values bit-for-bit'                                     , &
       &      all(zValuesWideLinear2D(1:size(zValuesNarrowLinear2D,dim=1),1:size(zValuesNarrowLinear2D,dim=2)) == zValuesNarrowLinear2D), &
       &      .true.                                                                                                                      &
       &     )
  call Assert('a linearly-spaced 2D table takes its spacing from the lattice'                                                           , &
       &      myTable2DLinearExtend%x(2)-myTable2DLinearExtend%x(1) == latticeX2D%step()                                                , &
       &      .true.                                                                                                                      &
       &     )
  call myTable2DLinearExtend%destroy()

  ! Test extension of a cubic-spline table. A cubic spline is not local - every coefficient depends on every tabulated value -
  ! so extension preserves the tabulated values but not the interpolant between them. What it must guarantee is that an extended
  ! table is indistinguishable from one built directly on the wider lattice, including in what it interpolates.
  allocate(table1DLogarithmicMonotoneCSpline :: myTable)
  select type (myTable)
  type is (table1DLogarithmicMonotoneCSpline)
     latticeNarrow=Range_Pinned(3.0d0,4,gridSchemePerOctave)
     call myTable%extend(latticeNarrow,isComputed)
     call Assert('spline extension of an empty table requires every point to be computed',count(isComputed),0                  )
     call Assert('spline extension of an empty table gives the lattice point count'      ,myTable%size()   ,latticeNarrow%count)
     do i=1,myTable%size()
        call myTable%populate(myTable%x(i)**2,i)
     end do
     xValuesSpline=myTable%xs()
     yValuesSpline=myTable%ys()
     ! Extend to a wider range on the same lattice.
     latticeWide=Range_Pinned(30.0d0,4,gridSchemePerOctave,latticeCurrent=myTable%lattice)
     call myTable%extend(latticeWide,isComputed)
     offset=latticeNarrow%indexMinimum-latticeWide%indexMinimum
     call Assert('spline extension marks precisely the previously computed points as computed',count(isComputed)                                   ,latticeNarrow%count)
     call Assert('spline extension marks the correct window as computed'                      ,all(isComputed(offset+1:offset+latticeNarrow%count)),.true.             )
     xValuesWide=myTable%xs()
     yValuesWide=myTable%ys()
     call Assert('spline extension preserves abscissae bit-for-bit'                            , &
          &      all(xValuesWide(offset+1:offset+latticeNarrow%count  ) == xValuesSpline      ), &
          &      .true.                                                                          &
          &     )
     call Assert('spline extension preserves tabulated values bit-for-bit'                     , &
          &      all(yValuesWide(offset+1:offset+latticeNarrow%count,1) == yValuesSpline(:,1) ), &
          &      .true.                                                                          &
          &     )
     ! Compute the newly added points, then record what the extended table holds and interpolates.
     do i=1,myTable%size()
        if (.not.isComputed(i)) call myTable%populate(myTable%x(i)**2,i)
     end do
     xValuesWide=myTable%xs()
     yValuesWide=myTable%ys()
     allocate(interpolatedExtended(latticeWide%count-1))
     do i=1,latticeWide%count-1
        interpolatedExtended(i)=myTable%interpolate(sqrt(myTable%x(i)*myTable%x(i+1)))
     end do
     call myTable%destroy()
  end select
  deallocate(myTable)

  ! Build the same table directly on the wider lattice. It must agree with the extended table not only in its tabulated values
  ! but in what it interpolates - the spline coefficients having been rebuilt over the whole range by both routes.
  allocate(table1DLogarithmicMonotoneCSpline :: myTable)
  select type (myTable)
  type is (table1DLogarithmicMonotoneCSpline)
     call myTable%extend(latticeWide,isComputed)
     do i=1,myTable%size()
        call myTable%populate(myTable%x(i)**2,i)
     end do
     xValuesDirect=myTable%xs()
     yValuesDirect=myTable%ys()
     call Assert('a spline table built directly has abscissae identical to one built by extension'       ,all(xValuesDirect      == xValuesWide     ),.true.)
     call Assert('a spline table built directly has tabulated values identical to one built by extension',all(yValuesDirect(:,1) == yValuesWide(:,1)),.true.)
     allocate(interpolatedDirect(latticeWide%count-1))
     do i=1,latticeWide%count-1
        interpolatedDirect(i)=myTable%interpolate(sqrt(myTable%x(i)*myTable%x(i+1)))
     end do
     call Assert('a spline table built directly interpolates identically to one built by extension'      ,all(interpolatedDirect == interpolatedExtended),.true.)
     call myTable%destroy()
  end select
  deallocate(myTable)

  ! A linearly-spaced cubic-spline table extends onto a `perUnit` lattice in the same way.
  allocate(table1DLinearCSpline :: myTable)
  select type (myTable)
  type is (table1DLinearCSpline)
     latticeNarrow=Range_Pinned(3.0d0,4,gridSchemePerUnit)
     call myTable%extend(latticeNarrow,isComputed)
     do i=1,myTable%size()
        call myTable%populate(myTable%x(i)**2,i)
     end do
     yValuesSpline=myTable%ys()
     latticeWide=Range_Pinned(9.0d0,4,gridSchemePerUnit,latticeCurrent=myTable%lattice)
     call myTable%extend(latticeWide,isComputed)
     offset=latticeNarrow%indexMinimum-latticeWide%indexMinimum
     yValuesWide=myTable%ys()
     call Assert('a linearly-spaced spline table preserves tabulated values bit-for-bit on extension' , &
          &      all(yValuesWide(offset+1:offset+latticeNarrow%count,1) == yValuesSpline(:,1)        ), &
          &      .true.                                                                                 &
          &     )
     call Assert('a linearly-spaced spline table takes its spacing from the lattice'                  , &
          &      myTable%x(2)-myTable%x(1) == latticeWide%step()                                      , &
          &      .true.                                                                                 &
          &     )
     call myTable%destroy()
  end select
  deallocate(myTable)

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

contains

  function Table_Type_Label(i)
    !!{RST
    Return a label for the ``i``\ :math:`^\mathrm{th}` table type tested in the "Table range" group.
    !!}
    implicit none
    character(len=32)                :: Table_Type_Label
    integer          , intent(in   ) :: i

    select case (i)
    case ( 1)
       Table_Type_Label='generic'
    case ( 2,3)
       Table_Type_Label='linear'
    case ( 4)
       Table_Type_Label='linear cubic spline'
    case ( 5)
       Table_Type_Label='monotone spline'
    case ( 6)
       Table_Type_Label='linear monotone spline'
    case ( 7,8)
       Table_Type_Label='logarithmic'
    case ( 9)
       Table_Type_Label='logarithmic cubic spline'
    case (10)
       Table_Type_Label='logarithmic monotone spline'
    case default
       Table_Type_Label='unknown'
    end select
    return
  end function Table_Type_Label

  subroutine Test_Table_Range(label,table_,logarithmic,zero)
    !!{RST
    Run a set of assertions over a one-dimensional table populated with :math:`y=1,\ldots,7` at its seven points. A linear table
    spans :math:`1 \le x \le 7`, so that :math:`y=x`. A logarithmic table spans :math:`1 \le x \le 1000`, so that :math:`y = 1 +
    s \log x` with :math:`s=6/\log(1000)`. Every interpolant reproduces these exactly within the range. Beyond the range, with
    "fix" extrapolation the value is that at the boundary and the gradient is zero, while with "zero" extrapolation (``zero`` is
    true) both are zero.
    !!}
    implicit none
    character       (len=*  ), intent(in   ) :: label
    class           (table1D), intent(inout) :: table_
    logical                  , intent(in   ) :: logarithmic           , zero
    double precision                         :: xWithin               , yWithin       , &
         &                                      gradientWithin        , xAbove        , &
         &                                      yAbove                , yBelow
    double precision         , parameter     :: xBelow        =0.5d0
    double precision         , parameter     :: absTol        =1.0d-12, relTol =1.0d-6

    if (logarithmic) then
       xWithin       =1.0d2
       yWithin       =1.0d0+6.0d0/log(1.0d3)*log(xWithin)
       gradientWithin=      6.0d0/log(1.0d3)/    xWithin
       xAbove        =2.0d3
    else
       xWithin       =4.5d0
       yWithin       =xWithin
       gradientWithin=1.0d0
       xAbove        =1.0d1
    end if
    if (zero) then
       yAbove=0.0d0
       yBelow=0.0d0
    else
       yAbove=7.0d0
       yBelow=1.0d0
    end if
    call Assert(label//': value within range'   ,table_%interpolate        (xWithin),yWithin       ,absTol=absTol,relTol=relTol)
    call Assert(label//': value above range'    ,table_%interpolate        (xAbove ),yAbove        ,absTol=absTol,relTol=relTol)
    call Assert(label//': value below range'    ,table_%interpolate        (xBelow ),yBelow        ,absTol=absTol,relTol=relTol)
    call Assert(label//': gradient within range',table_%interpolateGradient(xWithin),gradientWithin,absTol=absTol,relTol=relTol)
    call Assert(label//': gradient above range' ,table_%interpolateGradient(xAbove ),0.0d0         ,absTol=absTol,relTol=relTol)
    call Assert(label//': gradient below range' ,table_%interpolateGradient(xBelow ),0.0d0         ,absTol=absTol,relTol=relTol)
    return
  end subroutine Test_Table_Range

  subroutine Test_Table_2D_Range(label,table_,logarithmic,zero)
    !!{RST
    Run a set of assertions over a two-dimensional table populated with :math:`z=f(x)+2f(y)` at its :math:`7\times 7` points. A
    linear table spans :math:`1 \le x,y \le 7` with :math:`f(x)=x`. A logarithmic table spans :math:`1 \le x,y \le 1000` with
    :math:`f(x)=\log x`. Every interpolant reproduces these exactly within the range. Beyond the range of an axis, with "fix"
    extrapolation the value is that at the boundary of that axis and the gradient with respect to that axis is zero, while the
    gradient with respect to the other axis is unchanged. With "zero" extrapolation (``zero`` is true) the value and both
    gradients are zero.
    !!}
    implicit none
    character       (len=*  ), intent(in   ) :: label
    class           (table2D), intent(inout) :: table_
    logical                  , intent(in   ) :: logarithmic        , zero
    double precision                         :: xWithin            , yWithin         , &
         &                                      xAbove             , fWithinX        , &
         &                                      fWithinY           , fMinimum        , &
         &                                      fMaximum           , gradientX       , &
         &                                      gradientY          , beyond
    double precision         , parameter     :: xBelow     =0.5d0
    double precision         , parameter     :: absTol     =1.0d-12, relTol   =1.0d-6

    if (logarithmic) then
       xWithin  =1.0d2
       yWithin  =3.0d1
       xAbove   =2.0d3
       fWithinX =log(xWithin)
       fWithinY =log(yWithin)
       fMinimum =0.0d0
       fMaximum =log(1.0d3)
       gradientX=1.0d0/xWithin
       gradientY=2.0d0/yWithin
    else
       xWithin  =4.5d0
       yWithin  =2.5d0
       xAbove   =1.0d1
       fWithinX =xWithin
       fWithinY =yWithin
       fMinimum =1.0d0
       fMaximum =7.0d0
       gradientX=1.0d0
       gradientY=2.0d0
    end if
    ! Beyond the range every expected result is multiplied by this factor - unity under "fix" extrapolation, zero under "zero".
    if (zero) then
       beyond=0.0d0
    else
       beyond=1.0d0
    end if
    call Assert(label//': value within range'               ,table_%interpolate        (xWithin,yWithin  ),       fWithinX+2.0d0*fWithinY  ,absTol=absTol,relTol=relTol)
    call Assert(label//': x-gradient within range'          ,table_%interpolateGradient(xWithin,yWithin,1),       gradientX                ,absTol=absTol,relTol=relTol)
    call Assert(label//': y-gradient within range'          ,table_%interpolateGradient(xWithin,yWithin,2),       gradientY                ,absTol=absTol,relTol=relTol)
    call Assert(label//': value with x above range'         ,table_%interpolate        (xAbove ,yWithin  ),beyond*(fMaximum+2.0d0*fWithinY),absTol=absTol,relTol=relTol)
    call Assert(label//': x-gradient with x above range'    ,table_%interpolateGradient(xAbove ,yWithin,1),0.0d0                           ,absTol=absTol,relTol=relTol)
    call Assert(label//': y-gradient with x above range'    ,table_%interpolateGradient(xAbove ,yWithin,2),beyond*gradientY                ,absTol=absTol,relTol=relTol)
    call Assert(label//': value with x below range'         ,table_%interpolate        (xBelow ,yWithin  ),beyond*(fMinimum+2.0d0*fWithinY),absTol=absTol,relTol=relTol)
    call Assert(label//': x-gradient with x below range'    ,table_%interpolateGradient(xBelow ,yWithin,1),0.0d0                           ,absTol=absTol,relTol=relTol)
    call Assert(label//': y-gradient with x below range'    ,table_%interpolateGradient(xBelow ,yWithin,2),beyond*gradientY                ,absTol=absTol,relTol=relTol)
    call Assert(label//': value with y above range'         ,table_%interpolate        (xWithin,xAbove   ),beyond*(fWithinX+2.0d0*fMaximum),absTol=absTol,relTol=relTol)
    call Assert(label//': x-gradient with y above range'    ,table_%interpolateGradient(xWithin,xAbove ,1),beyond*gradientX                ,absTol=absTol,relTol=relTol)
    call Assert(label//': y-gradient with y above range'    ,table_%interpolateGradient(xWithin,xAbove ,2),0.0d0                           ,absTol=absTol,relTol=relTol)
    call Assert(label//': value with y below range'         ,table_%interpolate        (xWithin,xBelow   ),beyond*(fWithinX+2.0d0*fMinimum),absTol=absTol,relTol=relTol)
    call Assert(label//': x-gradient with y below range'    ,table_%interpolateGradient(xWithin,xBelow ,1),beyond*gradientX                ,absTol=absTol,relTol=relTol)
    call Assert(label//': y-gradient with y below range'    ,table_%interpolateGradient(xWithin,xBelow ,2),0.0d0                           ,absTol=absTol,relTol=relTol)
    call Assert(label//': value with x above, y below range',table_%interpolate        (xAbove ,xBelow   ),beyond*(fMaximum+2.0d0*fMinimum),absTol=absTol,relTol=relTol)
    call Assert(label//': x-gradient with x above, y below' ,table_%interpolateGradient(xAbove ,xBelow ,1),0.0d0                           ,absTol=absTol,relTol=relTol)
    call Assert(label//': y-gradient with x above, y below' ,table_%interpolateGradient(xAbove ,xBelow ,2),0.0d0                           ,absTol=absTol,relTol=relTol)
    return
  end subroutine Test_Table_2D_Range

  subroutine Test_Table_2D(label,table_)
    !!{RST
    Run a set of assertions over a two-dimensional table through the ``table2D`` base class. ``table_`` must already have been
    created, with two tables; it is populated, exercised, and destroyed here.

    The assertions are made without knowing which coordinate a given type interpolates in: interpolation at a tabulated point
    returns the tabulated value whatever the data, and the gradient is checked against a finite difference of the interpolant
    itself. The values tabulated vary quadratically with the index along each axis, so that the gradient differs from one cell
    to the next---which is what lets a gradient memoized at one point be told apart from the gradient at another.
    !!}
    implicit none
    character       (len=*  ), intent(in   )                 :: label
    class           (table2D), intent(inout)                 :: table_
    integer                                                  :: i               , j
    double precision         , allocatable  , dimension(:  ) :: xValues         , yValues
    double precision         , allocatable  , dimension(:,:) :: zValues
    double precision                                         :: delta            , gradient         , &
         &                                                      xMid             , yMid             , &
         &                                                      gradientMemoized , gradientReference, &
         &                                                      gradientDiscarded, valueDiscarded

    call Assert(label//': a created table reports itself initialized',table_%isInitialized()                 ,.true.)
    call Assert(label//': the table reports the size of each axis'   ,[table_%size(dim=1),table_%size(dim=2)],[11,8])
    ! Populate both tables, the second offset from the first so that the `table` argument is exercised.
    do i=1,table_%size(dim=1)
       do j=1,table_%size(dim=2)
          call table_%populate(        dble(i)**2+2.0d0*dble(j)   ,i,j        )
          call table_%populate(100.0d0+dble(i)   +2.0d0*dble(j)**2,i,j,table=2)
       end do
    end do
    ! The array accessors must agree with the element accessors.
    xValues=table_%xs(       )
    yValues=table_%ys(       )
    zValues=table_%zs(table=2)
    ! The abscissae are compared to within a few units in the last place rather than bit-for-bit: where the accessors apply a
    ! transformation to the internal coordinate - as the logarithmic types apply `exp` - the array accessor transforms a whole
    ! array and the element accessor a scalar, and at `-O3` those take vectorized and scalar paths through the mathematical
    ! library which need not agree in their final bit. The values themselves are copied rather than transformed, so those must
    ! agree exactly.
    call Assert(label//': `xs` agrees with `x`',xValues,[(table_%x(i),i=1,table_%size(dim=1))],relTol=1.0d-14)
    call Assert(label//': `ys` agrees with `y`',yValues,[(table_%y(j),j=1,table_%size(dim=2))],relTol=1.0d-14)
    call Assert(label//': `zs` agrees with `z`',all([((zValues(i,j) == table_%z(i,j,table=2),i=1,table_%size(dim=1)),j=1,table_%size(dim=2))]),.true.)
    ! Interpolation at a tabulated point must return the tabulated value, in either table.
    call Assert(label//': interpolation at a tabulated point returns the tabulated value'                        , &
         &      [table_%interpolate(table_%x(4),table_%y(3)),table_%interpolate(table_%x(4),table_%y(3),table=2)], &
         &      [table_%z          (        4  ,        3  ),table_%z          (        4  ,        3  ,table=2)], &
         &      absTol=1.0d-9                                                                                      &
         &     )
    ! The gradient returned must be that of the interpolant, in each dimension and in either table. Both points are taken well
    ! inside a single cell, where the interpolant is smooth.
    xMid    =table_%x(4)+0.5d0*(table_%x(5)-table_%x(4))
    yMid    =table_%y(3)+0.5d0*(table_%y(4)-table_%y(3))
    delta   =1.0d-4     *(table_%x(5)-table_%x(4))
    gradient=+(                                             &
         &     +table_%interpolate(xMid+delta,yMid,table=2) &
         &     -table_%interpolate(xMid-delta,yMid,table=2) &
         &    )                                             &
         &   /(2.0d0*delta)
    call Assert(label//': the gradient in the first dimension is that of the interpolant' ,table_%interpolateGradient(xMid,yMid,1,table=2),gradient,relTol=1.0d-6)
    delta   =1.0d-4     *(table_%y(4)-table_%y(3))
    gradient=+(                                             &
         &     +table_%interpolate(xMid,yMid+delta,table=2) &
         &     -table_%interpolate(xMid,yMid-delta,table=2) &
         &    )                                             &
         &   /(2.0d0*delta)
    call Assert(label//': the gradient in the second dimension is that of the interpolant',table_%interpolateGradient(xMid,yMid,2,table=2),gradient,relTol=1.0d-6)
    ! Value and gradient are distinct quantities computed from the same interpolation factors, and each table is distinct from
    ! the others, so each must be memoized under its own key. Interleaving them at a single point - where a shared key would
    ! report a hit - must therefore still answer each request with its own quantity, from its own table. A tabulated point is
    ! used, since interpolation there must return the tabulated value exactly.
    gradientMemoized=table_%interpolateGradient(table_%x(4),table_%y(3),1,table=1)
    call Assert(label//': a value interpolated after a gradient at the same point is not the memoized gradient', &
         &      table_%interpolate        (table_%x(4),table_%y(3)  ,table=1)                                  , &
         &      table_%z                  (         4 ,         3   ,table=1)                                  , &
         &      absTol=1.0d-9                                                                                    &
         &     )
    call Assert(label//': switching table returns the value from that table'                                   , &
         &      table_%interpolate        (table_%x(4),table_%y(3)  ,table=2)                                  , &
         &      table_%z                  (         4 ,         3   ,table=2)                                  , &
         &      absTol=1.0d-9                                                                                    &
         &     )
    call Assert(label//': switching back returns the value from the original table'                            , &
         &      table_%interpolate        (table_%x(4),table_%y(3)  ,table=1)                                  , &
         &      table_%z                  (         4 ,         3   ,table=1)                                  , &
         &      absTol=1.0d-9                                                                                    &
         &     )
    call Assert(label//': a gradient interpolated after a value at the same point is not the memoized value'   , &
         &      table_%interpolateGradient(table_%x(4),table_%y(3),1,table=1)                                  , &
         &      gradientMemoized                                                                                 &
         &     )
    ! A memoized value must be discarded when the values it was computed from are replaced.
    call table_%populate(1.0d0+table_%z(4,3,table=1),4,3,table=1)
    call Assert(label//': a memoized value is discarded when the table is repopulated'                         , &
         &      table_%interpolate        (table_%x(4),table_%y(3)  ,table=1)                                  , &
         &      table_%z                  (         4 ,         3   ,table=1)                                  , &
         &      absTol=1.0d-9                                                                                    &
         &     )
    ! Value and gradient share the interpolation factors from which both are computed, but must be memoized under keys of
    ! their own. Evaluating either at a point advances those shared factors, and must not thereby make the other - memoized at
    ! some earlier and different point - answerable at the new one. Tabulated points are used, so that interpolation at each
    ! returns the tabulated value exactly.
    valueDiscarded   =table_%interpolate        (table_%x(4),table_%y(3)  ,table=1)
    gradientDiscarded=table_%interpolateGradient(table_%x(6),table_%y(5),1,table=1)
    call Assert(label//': a value is not memoized across a gradient taken at another point'                    , &
         &      table_%interpolate        (table_%x(6),table_%y(5)  ,table=1)                                  , &
         &      table_%z                  (         6 ,         5   ,table=1)                                  , &
         &      absTol=1.0d-9                                                                                    &
         &     )
    ! The same in the other direction. The reference gradient is taken with the memoized state freshly discarded, which
    ! repopulating a point with the value it already holds achieves without altering the table.
    call table_%populate(table_%z(1,1,table=1),1,1,table=1)
    gradientReference=table_%interpolateGradient(table_%x(6),table_%y(5),1,table=1)
    call table_%populate(table_%z(1,1,table=1),1,1,table=1)
    gradientDiscarded=table_%interpolateGradient(table_%x(4),table_%y(3),1,table=1)
    valueDiscarded   =table_%interpolate        (table_%x(6),table_%y(5)  ,table=1)
    call Assert(label//': a gradient is not memoized across a value taken at another point'                    , &
         &      table_%interpolateGradient(table_%x(6),table_%y(5),1,table=1)                                  , &
         &      gradientReference                                                                                &
         &     )
    ! A destroyed table must no longer report itself initialized. Destruction also discards the memoized interpolation state
    ! along with the values from which it was computed, so that a later interpolation cannot be answered from the memo.
    call table_%destroy()
    call Assert(label//': a destroyed table no longer reports itself initialized',table_%isInitialized(),.false.)
    return
  end subroutine Test_Table_2D

end program Test_Tables
