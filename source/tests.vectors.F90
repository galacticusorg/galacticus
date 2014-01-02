!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a program to test vector functions.

program Test_Vectors
  !% Tests of vector functions.
  use Unit_Tests
  use Vectors
  implicit none

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Vectors")

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

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Vectors
