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
Contains a program to test vector functions.
!!}

program Test_Vectors
  !!{
  Tests of vector functions.
  !!}
  use :: Display           , only : displayIndent         , displayMessage        , displayUnindent     , displayVerbositySet, &
          &                         verbosityLevelStandard
  use :: ISO_Varying_String, only : assignment(=)         , varying_string
  use :: Kind_Numbers      , only : kind_int8
  use :: Unit_Tests        , only : Assert                , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  use :: Vectors           , only : Vector_Magnitude      , Vector_Outer_Product  , Vector_Product
  implicit none
  double precision                , allocatable, dimension(:  ) :: vector1   , vector2
  double precision                , allocatable, dimension(:,:) :: matrix12
  type            (varying_string)                              :: message
  character       (len= 3        )                              :: units
  character       (len=24        )                              :: label
  integer         (kind=kind_int8)                              :: countStart, countEnd , countRate

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Vectors")

  ! Begin benchmarks.
  call displayIndent("Benchmarks:")

  ! Test vector magnitude functions.
  call Assert('vector magnitude',                                                                                                   &
       &      [Vector_Magnitude([1.0d0,2.0d0,3.0d0]),Vector_Magnitude([1.0d0,2.0d0,-3.0d0]),Vector_Magnitude([0.0d0,0.0d0,0.0d0])], &
       &      [sqrt(14.0d0)                         ,sqrt(14.0d0)                          ,0.0d0                                ]  &
       &     )

  ! Test vector products.
  call Assert('vector product',                                        &
       &      Vector_Product([1.0d0,2.0d0,3.0d0],[3.0d0,2.0d0,1.0d0]), &
       &      [-4.0d0,8.0d0,-4.0d0]                                    &
       &     )
  call Assert('vector product with self',                              &
       &      Vector_Product([1.0d0,2.0d0,3.0d0],[1.0d0,2.0d0,3.0d0]), &
       &      [0.0d0,0.0d0,0.0d0]                                      &
       &     )

  ! Test and benchmark vector outer products.
  allocate(vector1(1000))
  allocate(vector2(1000))
  vector1=2.0d0
  vector2=0.5d0
  call System_Clock(countStart,countRate)
  matrix12=Vector_Outer_Product(vector1,vector2)
  call System_Clock(countEnd  ,countRate)
  select case (countRate)
  case (      1000)
     units="ms"
  case (   1000000)
     units="μs"
  case (1000000000)
     units="ns"
  end select
  write (label,'(i20)') countEnd-countStart
  message="Vector outer product         : "//trim(label)//" "//trim(units)
  call displayMessage(message)
  call Assert('vector outer product'           , &
       &      all(abs(matrix12-1.0d0) < 1.0d-6), &
       &      .true.                             &
       &     )

  ! Test self vector outer product.
  matrix12=Vector_Outer_Product(vector1,symmetrize=.true.)
  call Assert('vector outer product'           , &
       &      all(abs(matrix12-4.0d0) < 1.0d-6), &
       &      .true.                             &
       &     )
  ! Clean up.
  deallocate(vector1,vector2,matrix12)

  ! End benchmarks.
  call displayUnindent("done")

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Vectors
