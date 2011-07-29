!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


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
