!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a program to test the numerical interpolation code.

program Test_Interpolation
  !% Tests that numerical interpolation code works correctly.
  use Unit_Tests
  use Numerical_Interpolation
  use FGSL
  implicit none
  type(fgsl_interp)                      :: interpolationObject
  type(fgsl_interp_accel)                :: interpolationAccelerator
  logical                                :: interpolationReset=.true.
  double precision,       dimension(0:9) :: xArray=[1.0d0,3.0d0  ,3.3d0,4.3d0,6.7d0, 7.2d0, 8.9d0, 9.1d0,12.0d0,13.0d0]
  double precision,       dimension(0:9) :: yArray=[2.0d0,3.0d0,-23.0d0,4.0d0,6.0d0,-1.0d0,-5.0d0,-0.1d0, 5.0d0, 9.0d0]
  double precision                       :: x,y

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Numerical interpolation")
  
  ! Test interpolations.
  x=5.5d0
  y=Interpolate(10,xArray,yArray,interpolationObject,interpolationAccelerator,x,reset=interpolationReset)
  call Assert("linear interpolation",y,5.0d0)
   
  ! Test derivative interpolations.
  x=5.5d0
  y=Interpolate_Derivative(10,xArray,yArray,interpolationObject,interpolationAccelerator,x,reset=interpolationReset)
  call Assert("linear derivative interpolation",y,2.0d0/2.4d0,relTol=1.0d-6)
   
  ! Test linear extrapolation.
  x=15.0d0
  y=Interpolate(10,xArray,yArray,interpolationObject,interpolationAccelerator,x,reset=interpolationReset,extrapolationType=extrapolationTypeLinear)
  call Assert("linear extrapolation",y,17.0d0)
   
  ! Test fixed extrapolation.
  x=15.0d0
  y=Interpolate(10,xArray,yArray,interpolationObject,interpolationAccelerator,x,reset=interpolationReset,extrapolationType=extrapolationTypeFixed)
  call Assert("fixed extrapolation",y,9.0d0)
   
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Interpolation
