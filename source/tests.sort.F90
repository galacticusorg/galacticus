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
Contains a program to test sorting functions.
!!}

program Test_Sort
  !!{
  Tests of sorting functions.
  !!}
  use            :: Display      , only : displayVerbositySet, verbosityLevelStandard
  use, intrinsic :: ISO_C_Binding, only : c_size_t
  use            :: Kind_Numbers , only : kind_int8
  use            :: Sorting      , only : sort               , sortByIndex           , sortIndex           , sortSmallest     , &
       &                                  sortLargest        , sortSmallestIndex     , sortLargestIndex
  use            :: Unit_Tests   , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  integer         (kind=c_size_t ), dimension(19) :: indexArray
  integer                         , dimension(19) :: integerArray
  integer         (kind=kind_int8), dimension(19) :: longIntegerArray
  double precision                , dimension(19) :: doubleArray
  double precision                , dimension( 5) :: doubleArraySort
  integer         (kind=c_size_t ), dimension( 5) :: indexArraySort
  logical                         , dimension(19) :: mask

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Sorting")

  ! Test integer sorting.
  integerArray=[-3,-9,-4,-6,-7,-2,-8,-5,-1,6,4,9,8,1,5,7,0,2,3]
  indexArray=sortIndex(integerArray)
  call Assert("integer index",int(indexArray),[2,7,5,4,8,3,1,6,9,17,14,18,19,11,15,10,16,13,12])
  call sort(integerArray)
  call Assert("integer sort",integerArray,[-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9])
 
  ! Test long integer sorting.
  longIntegerArray=[                                &
       &            +1223347737404150611_kind_int8, &
       &            +1371817143140704696_kind_int8, &
       &            -1700307937423120610_kind_int8, &
       &            +2108237805665842044_kind_int8, &
       &            +1973667863223201888_kind_int8, &
       &            +1990247277598863726_kind_int8, &
       &            -5756193962970214969_kind_int8, &
       &            +2074912360708357675_kind_int8, &
       &            +3041147994740233394_kind_int8, &
       &            +1465611346823512051_kind_int8, &
       &            +1192828201664220458_kind_int8, &
       &            +1317135667165572363_kind_int8, &
       &            +1889488732752974024_kind_int8, &
       &            +5302951902808101060_kind_int8, &
       &            -1795924725765858256_kind_int8, &
       &            +1992422947416728096_kind_int8, &
       &            +2356833083276067511_kind_int8, &
       &            +1931332895789691698_kind_int8, &
       &            -1767936648278096670_kind_int8  &
       &           ]
  call sort(longIntegerArray)
  call Assert(                                 &
       &      "long integer sort"            , &
       &      longIntegerArray               , &
       &      [                                &
       &       -5756193962970214969_kind_int8, &
       &       -1795924725765858256_kind_int8, &
       &       -1767936648278096670_kind_int8, &
       &       -1700307937423120610_kind_int8, &
       &       +1192828201664220458_kind_int8, &
       &       +1223347737404150611_kind_int8, &
       &       +1317135667165572363_kind_int8, &
       &       +1371817143140704696_kind_int8, &
       &       +1465611346823512051_kind_int8, &
       &       +1889488732752974024_kind_int8, &
       &       +1931332895789691698_kind_int8, &
       &       +1973667863223201888_kind_int8, &
       &       +1990247277598863726_kind_int8, &
       &       +1992422947416728096_kind_int8, &
       &       +2074912360708357675_kind_int8, &
       &       +2108237805665842044_kind_int8, &
       &       +2356833083276067511_kind_int8, &
       &       +3041147994740233394_kind_int8, &
       &       +5302951902808101060_kind_int8  &
       &      ]                                &
       &     )

  ! Test double sorting.
  doubleArray=[-3.0d0,-9.0d0,-4.0d0,-6.0d0,-7.0d0,-2.0d0,-8.0d0,-5.0d0,-1.0d0,6.0d0,4.0d0,9.0d0,8.0d0,1.0d0,5.0d0,7.0d0,0.0d0,2.0d0,3.0d0]
  indexArray=sortIndex(doubleArray)
  call Assert("double index",int(indexArray),[2,7,5,4,8,3,1,6,9,17,14,18,19,11,15,10,16,13,12])
  call sort       (doubleArray           )
  call Assert("double sort",doubleArray,[-9.0d0,-8.0d0,-7.0d0,-6.0d0,-5.0d0,-4.0d0,-3.0d0,-2.0d0,-1.0d0,0.0d0,1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0,8.0d0,9.0d0])
  doubleArray=[-3.0d0,-9.0d0,-4.0d0,-6.0d0,-7.0d0,-2.0d0,-8.0d0,-5.0d0,-1.0d0,6.0d0,4.0d0,9.0d0,8.0d0,1.0d0,5.0d0,7.0d0,0.0d0,2.0d0,3.0d0]
  call sortByIndex(doubleArray,indexArray)
  call Assert("double sort by index",doubleArray,[-9.0d0,-8.0d0,-7.0d0,-6.0d0,-5.0d0,-4.0d0,-3.0d0,-2.0d0,-1.0d0,0.0d0,1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0,8.0d0,9.0d0])

  ! Test sorting k smallest and largest elements.
  doubleArray=[-3.0d0,-9.0d0,-4.0d0,-6.0d0,-7.0d0,-2.0d0,-8.0d0,-5.0d0,-1.0d0,6.0d0,4.0d0,9.0d0,8.0d0,1.0d0,5.0d0,7.0d0,0.0d0,2.0d0,3.0d0]
  doubleArraySort=sortSmallest(doubleArray,5_c_size_t)
  call Assert("double sort smallest",doubleArraySort,[-9.0d0,-8.0d0,-7.0d0,-6.0d0,-5.0d0])
  doubleArraySort=sortLargest (doubleArray,5_c_size_t)
  call Assert("double sort largest" ,doubleArraySort,[ 9.0d0, 8.0d0, 7.0d0, 6.0d0, 5.0d0])
  indexArraySort=sortSmallestIndex(doubleArray,5_c_size_t)
  call Assert("double sort smallest index",int(indexArraySort),[2,7,5,4,8])
  indexArraySort=sortLargestIndex (doubleArray,5_c_size_t)
  call Assert("double sort largest index" ,int(indexArraySort),[12,13,16,10,15])

  ! Test sorting k smallest and largest elements with mask.
  doubleArray=[-3.0d0,-9.0d0,-4.0d0,-6.0d0,-7.0d0,-2.0d0,-8.0d0,-5.0d0,-1.0d0,6.0d0,4.0d0,9.0d0,8.0d0,1.0d0,5.0d0,7.0d0,0.0d0,2.0d0,3.0d0]
  mask=[.true.,.false.,.true.,.false.,.true.,.false.,.true.,.false.,.true.,.false.,.true.,.false.,.true.,.false.,.true.,.false.,.true.,.false.,.true.]
  doubleArraySort=sortSmallest(doubleArray,5_c_size_t,mask=mask)
  call Assert("double sort smallest with mask",doubleArraySort,[-8.0d0,-7.0d0,-4.0d0,-3.0d0,-1.0d0])
  doubleArraySort=sortLargest (doubleArray,5_c_size_t,mask=mask)
  call Assert("double sort largest with mask" ,doubleArraySort,[ 8.0d0, 5.0d0, 4.0d0, 3.0d0, 0.0d0])
  indexArraySort=sortSmallestIndex(doubleArray,5_c_size_t,mask=mask)
  call Assert("double sort smallest index with mask",int(indexArraySort),[7,5,3,1,9])
  indexArraySort=sortLargestIndex (doubleArray,5_c_size_t,mask=mask)
  call Assert("double sort largest index with mask" ,int(indexArraySort),[13,15,11,19,17])

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Sort
