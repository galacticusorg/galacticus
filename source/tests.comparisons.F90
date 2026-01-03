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
Contains a program to test numerical comparison functions.
!!}

program Test_Comparison
  !!{
  Tests that numerical comparison functions work.
  !!}
  use :: Display             , only : displayVerbositySet, verbosityLevelStandard
  use :: Numerical_Comparison, only : Values_Agree       , Values_Differ
  use :: Unit_Tests          , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("numerical comparison")

  ! Check results.
  call Assert('double values agree' ,[                                                           &
       &                              Values_Agree ( 1.0d0, 1.000d0                            ), &
       &                              Values_Agree ( 1.0d0, 2.000d0                            ), &
       &                              Values_Agree ( 1.0d0, 1.005d0,absTol=0.01d0              ), &
       &                              Values_Agree ( 1.0d0, 1.015d0,absTol=0.01d0              ), &
       &                              Values_Agree (10.0d0,10.050d0,              relTol=0.01d0), &
       &                              Values_Agree (10.0d0,10.150d0,              relTol=0.01d0), &
       &                              Values_Agree (10.0d0,10.050d0,absTol=0.10d0,relTol=0.01d0), &
       &                              Values_Agree (10.0d0,10.150d0,absTol=0.20d0,relTol=0.01d0), &
       &                              Values_Agree (10.0d0,10.150d0,absTol=0.10d0,relTol=0.02d0), &
       &                              Values_Agree (10.0d0,10.250d0,absTol=0.10d0,relTol=0.02d0)  &
       &                             ],                                                           &
       &                             [                                                            &
       &                              .true. ,                                                    &
       &                              .false.,                                                    &
       &                              .true. ,                                                    &
       &                              .false.,                                                    &
       &                              .true. ,                                                    &
       &                              .false.,                                                    &
       &                              .true. ,                                                    &
       &                              .true. ,                                                    &
       &                              .true. ,                                                    &
       &                              .false.                                                     &
       &                             ]                                                            &
       &     )
  call Assert('double values differ',[                                                            &
       &                              Values_Differ( 1.0d0, 1.000d0                            ), &
       &                              Values_Differ( 1.0d0, 2.000d0                            ), &
       &                              Values_Differ( 1.0d0, 1.005d0,absTol=0.01d0              ), &
       &                              Values_Differ( 1.0d0, 1.015d0,absTol=0.01d0              ), &
       &                              Values_Differ(10.0d0,10.050d0,              relTol=0.01d0), &
       &                              Values_Differ(10.0d0,10.150d0,              relTol=0.01d0), &
       &                              Values_Differ(10.0d0,10.050d0,absTol=0.10d0,relTol=0.01d0), &
       &                              Values_Differ(10.0d0,10.150d0,absTol=0.20d0,relTol=0.01d0), &
       &                              Values_Differ(10.0d0,10.150d0,absTol=0.10d0,relTol=0.02d0), &
       &                              Values_Differ(10.0d0,10.250d0,absTol=0.10d0,relTol=0.02d0)  &
       &                             ],                                                           &
       &                             [                                                            &
       &                              .false.,                                                    &
       &                              .true. ,                                                    &
       &                              .false.,                                                    &
       &                              .true. ,                                                    &
       &                              .false.,                                                    &
       &                              .true. ,                                                    &
       &                              .false.,                                                    &
       &                              .true. ,                                                    &
       &                              .true. ,                                                    &
       &                              .true.                                                      &
       &                              ]                                                           &
       &     )
  call Assert('real values agree'   ,[                                                           &
       &                              Values_Agree ( 1.0e0, 1.000e0                            ), &
       &                              Values_Agree ( 1.0e0, 2.000e0                            ), &
       &                              Values_Agree ( 1.0e0, 1.005e0,absTol=0.01e0              ), &
       &                              Values_Agree ( 1.0e0, 1.015e0,absTol=0.01e0              ), &
       &                              Values_Agree (10.0e0,10.050e0,              relTol=0.01e0), &
       &                              Values_Agree (10.0e0,10.150e0,              relTol=0.01e0), &
       &                              Values_Agree (10.0e0,10.050e0,absTol=0.10e0,relTol=0.01e0), &
       &                              Values_Agree (10.0e0,10.150e0,absTol=0.20e0,relTol=0.01e0), &
       &                              Values_Agree (10.0e0,10.150e0,absTol=0.10e0,relTol=0.02e0), &
       &                              Values_Agree (10.0e0,10.250e0,absTol=0.10e0,relTol=0.02e0)  &
       &                             ],                                                           &
       &                             [                                                            &
       &                              .true. ,                                                    &
       &                              .false.,                                                    &
       &                              .true. ,                                                    &
       &                              .false.,                                                    &
       &                              .true. ,                                                    &
       &                              .false.,                                                    &
       &                              .true. ,                                                    &
       &                              .true. ,                                                    &
       &                              .true. ,                                                    &
       &                              .false.                                                     &
       &                             ]                                                            &
       &     )
  call Assert('real values differ'  ,[                                                            &
       &                              Values_Differ( 1.0e0, 1.000e0                            ), &
       &                              Values_Differ( 1.0e0, 2.000e0                            ), &
       &                              Values_Differ( 1.0e0, 1.005e0,absTol=0.01e0              ), &
       &                              Values_Differ( 1.0e0, 1.015e0,absTol=0.01e0              ), &
       &                              Values_Differ(10.0e0,10.050e0,              relTol=0.01e0), &
       &                              Values_Differ(10.0e0,10.150e0,              relTol=0.01e0), &
       &                              Values_Differ(10.0e0,10.050e0,absTol=0.10e0,relTol=0.01e0), &
       &                              Values_Differ(10.0e0,10.150e0,absTol=0.20e0,relTol=0.01e0), &
       &                              Values_Differ(10.0e0,10.150e0,absTol=0.10e0,relTol=0.02e0), &
       &                              Values_Differ(10.0e0,10.250e0,absTol=0.10e0,relTol=0.02e0)  &
       &                             ],                                                           &
       &                             [                                                            &
       &                              .false.,                                                    &
       &                              .true. ,                                                    &
       &                              .false.,                                                    &
       &                              .true. ,                                                    &
       &                              .false.,                                                    &
       &                              .true. ,                                                    &
       &                              .false.,                                                    &
       &                              .true. ,                                                    &
       &                              .true. ,                                                    &
       &                              .true.                                                      &
       &                              ]                                                           &
       &     )

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Comparison
