!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

!% Contains a program to test linear algebra functions.

program Test_Math_Linear_Algebra
  !% Tests of linear algebra functions.
  use :: Galacticus_Display      , only : Galacticus_Verbosity_Level_Set   , verbosityStandard
  use :: Linear_Algebra          , only : vector                           , matrix                , assignment(=), matrixRotation
  use :: Numerical_Constants_Math, only : Pi
  use :: Unit_Tests              , only : Assert                           , Unit_Tests_Begin_Group               , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  double precision        , dimension(3  ) :: vectorComponents
  type            (vector), dimension(3  ) :: vectors                 , vectorsRotated
  double precision        , dimension(3,3) :: matrixRotationComponents
  type            (vector)                 :: vector1                 , vector2
  double precision                         :: angle

  ! Set verbosity level.
  call Galacticus_Verbosity_Level_Set(verbosityStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Math: linear algebra")

  ! Vector cross product.
  vector1=[+3.0d0,-1.0d0,+4.5d0]
  vector2=[+6.0d0,+2.3d0,+4.0d0]
  vectorComponents=vector1.cross.vector2
  call Assert("vector cross product",vectorComponents,[-287.0d0/20.0d0,15.0d0,129.0d0/10.0d0],relTol=1.0d-9)

  ! Test rotation matrix construction.
  vectors       (1)=[+1.0d0,+0.0d0,+0.0d0]
  vectors       (2)=[+0.0d0,+1.0d0,+0.0d0]
  vectors       (3)=[+0.0d0,+0.0d0,+1.0d0]
  vectorsRotated(1)=[+0.0d0,+1.0d0,+0.0d0]
  vectorsRotated(2)=[-1.0d0,+0.0d0,+0.0d0]
  vectorsRotated(3)=[+0.0d0,+0.0d0,+1.0d0]
  angle            =-Pi/2.0d0
  matrixRotationComponents=matrixRotation(vectors,vectorsRotated)
  call Assert("rotation matrix construct (-π/2, around z-axis)", &
       &      matrixRotationComponents                         , &
       &      reshape(                                           &
       &              [                                          &
       &               +cos(angle),-sin(angle),+0.0d0,           &
       &               +sin(angle),+cos(angle),+0.0d0,           &
       &               +    0.0d0 ,+    0.0d0 ,+1.0d0            &
       &              ]                                        , &
       &              [3,3]                                      &
       &             )                                         , &
       &      absTol=1.0d-9                                    , &
       &      relTol=1.0d-9                                      &
       &     )

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

end program Test_Math_Linear_Algebra
