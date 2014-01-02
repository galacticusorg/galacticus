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

!% Contains a program to test integration routines.

program Test_Integration
  !% Tests that numerical integration routines work.
  use, intrinsic :: ISO_C_Binding
  use Unit_Tests
  use Numerical_Integration
  use Test_Integration_Functions
  use Numerical_Constants_Math
  implicit none
  double precision                             :: integral
  type            (fgsl_function             ) :: integrandFunction
  type            (fgsl_integration_workspace) :: integrationWorkspace
  type            (c_ptr                     ) :: parameterPointer
  logical                                      :: integrationReset

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
