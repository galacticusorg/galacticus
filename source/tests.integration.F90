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


!% Contains a program to test integration routines.

program Test_Integration
  !% Tests that numerical integration routines work.
  use, intrinsic :: ISO_C_Binding
  use Unit_Tests
  use Numerical_Integration
  use FGSL
  use Test_Integration_Functions
  use Numerical_Constants_Math
  implicit none
  double precision                 :: integral
  type(fgsl_function)              :: integrandFunction
  type(fgsl_integration_workspace) :: integrationWorkspace
  type(c_ptr)                      :: parameterPointer
  logical                          :: integrationReset

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Numerical integration")

  ! Test simple integrations.
  integrationReset=.true.
  integral=Integrate(0.0d0,1.0d0,Integrand1,parameterPointer,integrandFunction&
       &,integrationWorkspace,toleranceRelative=1.0d-6,reset=integrationReset)
  call Assert("integrate f(x)=x from 0 to 1",integral,0.5d0,relTol=1.0d-6)

  integrationReset=.true.
  integral=Integrate(0.0d0,2.0d0*Pi,Integrand2,parameterPointer,integrandFunction&
       &,integrationWorkspace,toleranceRelative=1.0d-6,reset=integrationReset)
  call Assert("integrate f(x)=sin(x) from 0 to 2 Pi",integral,0.0d0,absTol=1.0d-6)

  integrationReset=.true.
  integral=Integrate(0.0d0,10.0d0,Integrand3,parameterPointer,integrandFunction&
       &,integrationWorkspace,toleranceRelative=1.0d-6,reset=integrationReset)
  call Assert("integrate f(x)=1/sqrt(x) from 0 to 10",integral,2.0d0*sqrt(10.0d0),relTol=1.0d-6)

  ! Test 2D integrations.
  integrationReset=.true.
  integral=Integrate(0.0d0,2.0d0*Pi,Integrand4,parameterPointer,integrandFunction&
       &,integrationWorkspace,toleranceRelative=1.0d-6,reset=integrationReset)
  call Assert("integrate f(x,y)=cos(x)*y from x=0 to 2 Pi and y=0 to x",integral,2.0d0*Pi,relTol=1.0d-6)

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Integration
