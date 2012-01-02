!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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


!% Contains a program to test the numerical range making code.

program Test_Make_Ranges
  !% Tests that numerical range making code works correctly.
  use Unit_Tests
  use Numerical_Ranges
  use Array_Utilities
  implicit none
  double precision, dimension(1:11) :: range1
  double precision, dimension(0:10) :: range2

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Numerical ranges")
  
  ! Create a linear range.
  range1=Make_Range(1.0d0,3.0d0,size(range1),rangeType=rangeTypeLinear)
  call Assert("linear range creation",range1,[1.0d0,1.2d0,1.4d0,1.6d0,1.8d0,2.0d0,2.2d0,2.4d0,2.6d0,2.8d0,3.0d0])
  
  ! Create a linear range in an offset array.
  range2=Make_Range(1.0d0,3.0d0,size(range2),rangeType=rangeTypeLinear)
  call Assert("linear range creation (offset array)",range2,[1.0d0,1.2d0,1.4d0,1.6d0,1.8d0,2.0d0,2.2d0,2.4d0,2.6d0,2.8d0,3.0d0])
  
  ! Create a logarithmic range.
  range1=Make_Range(10.0d0,1000.0d0,size(range1),rangeType=rangeTypeLogarithmic)
  call Assert("logarithmic range creation",range1,10.0d0**([1.0d0,1.2d0,1.4d0,1.6d0,1.8d0,2.0d0,2.2d0,2.4d0,2.6d0,2.8d0,3.0d0]),relTol=1.0d-6)

  ! Create a reverse order range.
  range1=Make_Range(3.0d0,1.0d0,size(range1),rangeType=rangeTypeLinear)
  call Assert("reverse linear range creation",range1,Array_Reverse([1.0d0,1.2d0,1.4d0,1.6d0,1.8d0,2.0d0,2.2d0,2.4d0,2.6d0,2.8d0,3.0d0]))

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Make_Ranges
