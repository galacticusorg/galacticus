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

!% Contains a module of functions for root finding unit tests.

module Test_Root_Finding_Functions
  !% Contains functions for root finding unit tests.
  implicit none
  private
  public :: Root_Function_1, Root_Function_2, Root_Function_2_Derivative, Root_Function_2_Both, Root_Function_3

contains

  double precision function Root_Function_1(x)
    !% Function for root finding unit tests.
    double precision, intent(in   ) :: x

    Root_Function_1=x
  end function Root_Function_1

  double precision function Root_Function_2(x)
    !% Function for root finding unit tests.
    double precision, intent(in   ) :: x

    Root_Function_2=x**2-5.0d0*x+1.0d0
  end function Root_Function_2

  double precision function Root_Function_2_Derivative(x)
    !% Function for root finding unit tests.
    double precision, intent(in   ) :: x

    Root_Function_2_Derivative=2.0d0*x-5.0d0
  end function Root_Function_2_Derivative

  subroutine  Root_Function_2_Both(x,f,df)
    !% Function for root finding unit tests.
    double precision, intent(in   ) :: x
    double precision, intent(  out) :: f, df

    f=x**2-5.0d0*x+1.0d0
    df=2.0d0*x-5.0d0
  end subroutine Root_Function_2_Both

  double precision function Root_Function_3(x)
    !% Function for root finding unit tests.
    double precision, intent(in   ) :: x

    Root_Function_3=x*exp(-x)+1.0d0
  end function Root_Function_3

end module Test_Root_Finding_Functions
