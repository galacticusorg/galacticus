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

!% Contains a module of functions for root finding unit tests.

module Test_Root_Finding_Functions
  !% Contains functions for root finding unit tests.
  use, intrinsic :: ISO_C_Binding
  implicit none
  private
  public :: Root_Function_1, Root_Function_2, Root_Function_3

contains

  function Root_Function_1(x,parameterPointer) bind(c)
    !% Function for root finding unit tests.
    real(c_double)          :: Root_Function_1
    real(c_double), value   :: x
    type(c_ptr),    value   :: parameterPointer
   
    Root_Function_1=x
  end function Root_Function_1

  function Root_Function_2(x,parameterPointer) bind(c)
    !% Function for root finding unit tests.
    real(c_double)          :: Root_Function_2
    real(c_double), value   :: x
    type(c_ptr),    value   :: parameterPointer
   
    Root_Function_2=x**2-5.0d0*x+1.0d0
  end function Root_Function_2

  function Root_Function_3(x,parameterPointer) bind(c)
    !% Function for root finding unit tests.
    real(c_double)          :: Root_Function_3
    real(c_double), value   :: x
    type(c_ptr),    value   :: parameterPointer
   
    Root_Function_3=x*exp(-x)+1.0d0
  end function Root_Function_3

end module Test_Root_Finding_Functions
