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
Contains a program to test convex hull functions.
!!}

program Test_Convex_Hulls
  !!{
  Tests of convex hull functions.
  !!}
  use :: Display           , only : displayVerbositySet, verbosityLevelStandard
  use :: Points_Convex_Hull, only : convexHull
  use :: Unit_Tests        , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish, &
       &                            Skip
  implicit none
  type(convexHull), allocatable :: hull

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Convex hulls")
  ! Test convex hull volume for a unit octahedron.
  call Unit_Tests_Begin_Group("Regular octahedron")
  allocate(hull)
#ifdef QHULLAVAIL
  hull=convexHull(                               &
       &          reshape(                       &
       &                  [                      &
       &                   +1.0d0,+0.0d0,+0.0d0, &
       &                   -1.0d0,+0.0d0,+0.0d0, &
       &                   +0.0d0,+1.0d0,+0.0d0, &
       &                   +0.0d0,-1.0d0,+0.0d0, &
       &                   +0.0d0,+0.0d0,+1.0d0, &
       &                   +0.0d0,+0.0d0,-1.0d0  &
       &                  ]                    , &
       &                  [                      &
       &                   3,6                   &
       &                  ]                      &
       &                 )                       &
       &         )
  call Assert("volume"                         ,hull%volume       (                          ),sqrt(2.0d0)*sqrt(2.0d0)**3/3.0d0)
  call Assert("point [0,0     ,0] inside hull" ,hull%pointIsInHull([0.0d0,0.0d0       ,0.0d0]),.true.                          )
  call Assert("point [2,0     ,0] outside hull",hull%pointIsInHull([2.0d0,0.0d0       ,0.0d0]),.false.                         )
  call Assert("point [½,½-10⁻⁶,0] inside hull" ,hull%pointIsInHull([0.5d0,0.5d0-1.0d-6,0.0d0]),.true.                          )
  call Assert("point [½,½+10⁻⁶,0] outside hull",hull%pointIsInHull([0.5d0,0.5d0+1.0d-6,0.0d0]),.false.                         )
#else
  call   Skip("volume"                         ,"qhull library unavailable")
  call   Skip("point [0,0     ,0] inside hull" ,"qhull library unavailable")
  call   Skip("point [2,0     ,0] outside hull","qhull library unavailable")
  call   Skip("point [½,½-10⁻⁶,0] inside hull" ,"qhull library unavailable")
  call   Skip("point [½,½+10⁻⁶,0] outside hull","qhull library unavailable")
#endif
  deallocate(hull)
  call Unit_Tests_End_Group()
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()
end program Test_Convex_Hulls
