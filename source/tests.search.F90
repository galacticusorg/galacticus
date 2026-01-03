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
Contains a program to test array search functions.
!!}

program Test_Search
  !!{
  Tests that array search functions work.
  !!}
  use :: Arrays_Search     , only : searchArray        , searchArrayClosest
  use :: Display           , only : displayVerbositySet, verbosityLevelStandard
  use :: ISO_Varying_String, only : assignment(=)      , var_str               , varying_string
  use :: Kind_Numbers      , only : kind_int8
  use :: Unit_Tests        , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  double precision                , dimension(10) :: myArray    =[0.0d0,1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0,8.0d0,9.0d0]
  double precision                , dimension(10) :: mySearch   =[3.4d0,9.0d0,4.2d0,-1.0d0,10.0d0,5.5d0,5.999999d0,6.000001d0,1.1d0,7.5d0]
  integer         (kind=kind_int8), dimension(10) :: myIntArray =[0,1,2,3,4,4,5,6,7,8]
  integer         (kind=kind_int8), dimension(10) :: myIntSearch=[0,1,2,3,4,5,6,7,8,8]
  integer                         , dimension(10) :: myIntExpect=[1,2,3,4,6,7,8,9,10,10]
  integer                         , dimension(10) :: myIndices
  type            (varying_string), dimension(26) :: stringArray
  integer                                         :: i

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Define an array of varying strings.
  stringArray=[            &
       &       'Alpha   ', &
       &       'Bravo   ', &
       &       'Charlie ', &
       &       'Delta   ', &
       &       'Echo    ', &
       &       'Foxtrot ', &
       &       'Golf    ', &
       &       'Hotel   ', &
       &       'India   ', &
       &       'Juliet  ', &
       &       'Kilo    ', &
       &       'Lima    ', &
       &       'Mike    ', &
       &       'November', &
       &       'Oscar   ', &
       &       'Papa    ', &
       &       'Quebec  ', &
       &       'Romeo   ', &
       &       'Sierra  ', &
       &       'Tango   ', &
       &       'Uniform ', &
       &       'Victor  ', &
       &       'Whiskey ', &
       &       'X-ray   ', &
       &       'Yankee  ', &
       &       'Zulu    '  &
       &      ]

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("array search")

  ! Test searching of double arrays.
  do i=1,size(mySearch)
     myIndices(i)=int(searchArray(myArray,mySearch(i)))
  end do
  call Assert('search double array',myIndices,max(min(int(mySearch)+1,9),1))
  do i=1,size(mySearch)
     myIndices(i)=int(searchArrayClosest(myArray,mySearch(i)))
  end do
  call Assert('search double array for closest match',myIndices,max(min(nint(mySearch)+1,10),1))

  ! Test searching of long integer arrays.
  do i=1,size(mySearch)
     myIndices(i)=int(searchArray(myIntArray,myIntSearch(i)))
  end do
  call Assert('search long integer array',myIndices,myIntExpect)

  ! Test searching of varying string arrays.
  call Assert('search string array',int(                                               &
       &                                [                                              &
       &                                 searchArray(stringArray,var_str('Monkey'  )), &
       &                                 searchArray(stringArray,var_str('Unicorn' )), &
       &                                 searchArray(stringArray,var_str('Beaver'  )), &
       &                                 searchArray(stringArray,var_str('Fox'     )), &
       &                                 searchArray(stringArray,var_str('Oscar'   )), &
       &                                 searchArray(stringArray,var_str('Shoggoth'))  &
       &                                ]                                              &
       &                               ),                                              &
       &                                [                                              &
       &                                 13,                                           &
       &                                 20,                                           &
       &                                  1,                                           &
       &                                  5,                                           &
       &                                 15,                                           &
       &                                 18                                            &
       &                                ]                                              &
       &     )

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Search
