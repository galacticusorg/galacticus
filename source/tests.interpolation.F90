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
Contains a program to test the numerical interpolation code.
!!}

program Test_Interpolation
  !!{
  Tests that numerical interpolation code works correctly.
  !!}
  use            :: Display                , only : displayVerbositySet         , verbosityLevelStandard
  use, intrinsic :: ISO_C_Binding          , only : c_size_t
  use            :: Numerical_Interpolation, only : interpolator                , interpolator2D
  use            :: Table_Labels           , only : extrapolationTypeExtrapolate, extrapolationTypeFix
  use            :: Unit_Tests             , only : Assert                      , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish, &
          &                                         compareGreaterThan          , compareLessThan
  implicit none
  type            (interpolator  ), allocatable     :: interpolator_
  type            (interpolator2D), allocatable     :: interpolator2D_
  double precision                , dimension(  10) :: xArray         =[1.0d0,3.0d0,  3.3d0,4.3d0,6.7d0, 7.2d0, 8.9d0, 9.1d0,12.0d0,13.0d0]
  double precision                , dimension(  10) :: yArray         =[2.0d0,3.0d0,-23.0d0,4.0d0,6.0d0,-1.0d0,-5.0d0,-0.1d0, 5.0d0, 9.0d0]
  double precision                , dimension(   3) :: x2Array        =[1.0d0,3.0d0,4.3d0]
  double precision                , dimension(   3) :: y2Array        =[2.0d0,3.0d0,6.0d0]
  double precision                , dimension(3, 3) :: z2Array        =reshape([1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0,8.0d0,9.0d0],[3,3])
  double precision                                  :: x                                                                                     , y, &
       &                                               z
  integer         (c_size_t      )                  :: i
  
  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Numerical interpolation")

  ! Begin tests of 1D interpolator.
  call Unit_Tests_Begin_Group("1D interpolator")

  ! Test location.
  allocate(interpolator_)
  interpolator_=interpolator(xArray,yArray)
  x=5.5d0
  i=interpolator_%locate(x)
  deallocate(interpolator_)
  call Assert("location (lower bound)",xArray(i           ),x,compare=compareLessThan   )
  call Assert("location (upper bound)",xArray(i+1_c_size_t),x,compare=compareGreaterThan)

  ! Test interpolations.
  allocate(interpolator_)
  interpolator_=interpolator(xArray,yArray)
  x=5.5d0
  y=interpolator_%interpolate(x)
  deallocate(interpolator_)
  call Assert("linear interpolation",y,5.0d0)

  ! Test derivative interpolations.
  allocate(interpolator_)
  interpolator_=interpolator(xArray,yArray)
  x=5.5d0
  y=interpolator_%derivative(x)
  deallocate(interpolator_)
  call Assert("linear derivative interpolation",y,2.0d0/2.4d0,relTol=1.0d-6)

  ! Test linear extrapolation.
  allocate(interpolator_)
  interpolator_=interpolator(xArray,yArray,extrapolationType=[extrapolationTypeExtrapolate,extrapolationTypeExtrapolate])
  x=15.0d0
  y=interpolator_%interpolate(x)
  deallocate(interpolator_)
  call Assert("linear extrapolation",y,17.0d0)

  ! Test fixed extrapolation.
  allocate(interpolator_)
  interpolator_=interpolator(xArray,yArray,extrapolationType=extrapolationTypeFix)
  x=15.0d0
  y=interpolator_%interpolate(x)
  deallocate(interpolator_)
  call Assert("fixed extrapolation",y,9.0d0)

  ! Done 1D.
  call Unit_Tests_End_Group()

  ! Begin tests of 2D interpolator.
  call Unit_Tests_Begin_Group("2D interpolator")

  ! Test interpolations.
  allocate(interpolator2D_)
  interpolator2D_=interpolator2D(x2Array,y2Array,z2Array)
  x=3.4d0
  y=2.7d0
  z=interpolator2D_%interpolate(x,y)
  deallocate(interpolator2D_)
  call Assert("linear interpolation",z,4.407692307692d0,relTol=1.0d-6)

  ! Done 2D.
  call Unit_Tests_End_Group()

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Interpolation
