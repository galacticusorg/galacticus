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

!% Contains a program to test sorting functions.

program Test_Sort
  !% Tests of sorting functions.
  use Unit_Tests
  use Sort
  use, intrinsic :: ISO_C_Binding
  implicit none
  integer         (kind=c_size_t), dimension(19) :: indexArray
  integer                        , dimension(19) :: integerArray
  double precision               , dimension(19) :: doubleArray

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Sorting")

  ! Test integer sorting.
  integerArray=[-3,-9,-4,-6,-7,-2,-8,-5,-1,6,4,9,8,1,5,7,0,2,3]
  indexArray=Sort_Index_Do(integerArray)
  call Assert("integer index",int(indexArray),[2,7,5,4,8,3,1,6,9,17,14,18,19,11,15,10,16,13,12])
  call Sort_Do(integerArray)
  call Assert("integer sort",integerArray,[-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9])

  ! Test double sorting.
  doubleArray=[-3.0d0,-9.0d0,-4.0d0,-6.0d0,-7.0d0,-2.0d0,-8.0d0,-5.0d0,-1.0d0,6.0d0,4.0d0,9.0d0,8.0d0,1.0d0,5.0d0,7.0d0,0.0d0,2.0d0,3.0d0]
  indexArray=Sort_Index_Do(doubleArray)
  call Assert("double index",int(indexArray),[2,7,5,4,8,3,1,6,9,17,14,18,19,11,15,10,16,13,12])
  call Sort_Do(doubleArray)
  call Assert("double sort",doubleArray,[-9.0d0,-8.0d0,-7.0d0,-6.0d0,-5.0d0,-4.0d0,-3.0d0,-2.0d0,-1.0d0,0.0d0,1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0,8.0d0,9.0d0])

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Sort
