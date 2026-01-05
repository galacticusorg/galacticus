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
Contains a program to test mathematical special functions.
!!}

program Test_Math_Fast
  !!{
  Tests of mathematical fast functions.
  !!}
  use :: Display            , only : displayVerbositySet, verbosityLevelStandard
  use :: Kind_Numbers       , only : kind_int8
  use :: Math_Exponentiation, only : cubeRoot           , fastExponentiator
  use :: Unit_Tests         , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  integer                            , parameter             :: pointCount  =10000
  double precision                   , dimension(pointCount) :: exact                  , fast                , &
       &                                                        x
  double precision                   , parameter             :: rangeMinimum=    1.00d0, rangeMaximum =10.0d0, &
       &                                                        exponent        =2.59d0
  type            (fastExponentiator)                        :: fastPower
  integer                                                    :: i
  integer         (kind_int8        )                        :: timeBeforeExact        , timeAfterExact      , &
       &                                                        timeBeforeFast         , timeAfterFast
  character       (len=7            )                        :: label

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Math: fast functions")

  ! Create array of x values at which to evaluate functions.
  fastPower=fastExponentiator(rangeMinimum,rangeMaximum,exponent,100.0d0,.true.)
  forall(i=1:pointCount)
     x(i)=rangeMinimum+(rangeMaximum-rangeMinimum)*dble(i-1)/dble(pointCount-1)
  end forall

  ! Test fast exponentiator.
  ! Time the intrinsic exponentiator.
  call system_clock(timeBeforeExact)
  exact=x**exponent
  call system_clock(timeAfterExact )
  ! Time the fast exponentiator.
  call system_clock(timeBeforeFast )
  do i=1,pointCount
     fast(i)=fastPower%exponentiate(x(i))
  end do
  call system_clock(timeAfterFast  )
  write (label,'(f7.3)') 100.0d0*(dble(timeAfterExact-timeBeforeExact)/dble(timeAfterFast-timeBeforeFast)-1.0d0)
  call Assert("Fast xʸ [speed up "//trim(label)//"%]",exact,fast,relTol=1.0d-4)

  ! Test cube root function.
  ! Time the intrinsic cube root function.
  call system_clock(timeBeforeExact)
  exact=x**(1.0d0/3.0d0)
  call system_clock(timeAfterExact )
  ! Time the fast cube root function.
  call system_clock(timeBeforeFast )
  forall(i=1:pointCount)
     fast(i)=cubeRoot(x(i))
  end forall
  call system_clock(timeAfterFast  )
  write (label,'(f7.3)') 100.0d0*(dble(timeAfterExact-timeBeforeExact)/dble(timeAfterFast-timeBeforeFast)-1.0d0)
  call Assert("Fast ∛x [speed up "//trim(label)//"%]",exact,fast,relTol=1.0d-6)

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Math_Fast
