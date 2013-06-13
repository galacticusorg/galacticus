!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module of integrands for unit tests.

module Test_Integration_Functions
  !% Contains integrands for unit tests.
  use, intrinsic :: ISO_C_Binding
  use FGSL
  implicit none
  private
  public :: Integrand1, Integrand2, Integrand3, Integrand4

  type(fgsl_function),              save :: integrandFunction
  type(fgsl_integration_workspace), save :: integrationWorkspace
  logical,                          save :: integrationReset=.true.
  
contains
  
  function Integrand1(x,parameterPointer) bind(c)
    !% Integral for unit testing.
    implicit none
    real(c_double)        :: Integrand1
    real(c_double), value :: x
    type(c_ptr),    value :: parameterPointer
    
    Integrand1=x
    return
  end function Integrand1
  
  function Integrand2(x,parameterPointer) bind(c)
    !% Integral for unit testing.
    implicit none
    real(c_double)        :: Integrand2
    real(c_double), value :: x
    type(c_ptr),    value :: parameterPointer
    
    Integrand2=sin(x)
    return
  end function Integrand2
  
  function Integrand3(x,parameterPointer) bind(c)
    !% Integral for unit testing.
    implicit none
    real(c_double)        :: Integrand3
    real(c_double), value :: x
    type(c_ptr),    value :: parameterPointer
    
    Integrand3=1.0d0/sqrt(x)
    return
  end function Integrand3
  
  function Integrand4(x,parameterPointer) bind(c)
    !% Integral for unit testing.
    use Numerical_Integration
    implicit none
    real(c_double)        :: Integrand4
    real(c_double), value :: x
    type(c_ptr),    value :: parameterPointer
    
    Integrand4=cos(x)*Integrate(0.0d0,x,Integrand1,parameterPointer,integrandFunction&
       &,integrationWorkspace,toleranceRelative=1.0d-6,reset=integrationReset)
    return
  end function Integrand4
  
end module Test_Integration_Functions
