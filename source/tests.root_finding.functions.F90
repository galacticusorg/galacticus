!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
!!    Andrew Benson <abenson@carnegiescience.edu>
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

!!{
Contains a module of functions for root finding unit tests.
!!}

module Test_Root_Finding_Functions
  !!{
  Contains functions for root finding unit tests.
  !!}
  implicit none
  private
  public :: Root_Function_1, Root_Function_2, Root_Function_2_Derivative, Root_Function_2_Both, Root_Function_3, &
       & Root_Function_4, Root_Function_4_Derivative, Root_Function_4_Both

contains

  double precision function Root_Function_1(x)
    !!{
    Function for root finding unit tests.
    !!}
    double precision, intent(in   ) :: x

    Root_Function_1=x
  end function Root_Function_1

  double precision function Root_Function_2(x)
    !!{
    Function for root finding unit tests.
    !!}
    double precision, intent(in   ) :: x

    Root_Function_2=x**2-5.0d0*x+1.0d0
  end function Root_Function_2

  double precision function Root_Function_2_Derivative(x)
    !!{
    Function for root finding unit tests.
    !!}
    double precision, intent(in   ) :: x

    Root_Function_2_Derivative=2.0d0*x-5.0d0
  end function Root_Function_2_Derivative

  subroutine  Root_Function_2_Both(x,f,df)
    !!{
    Function for root finding unit tests.
    !!}
    double precision, intent(in   ) :: x
    double precision, intent(  out) :: f, df

    f=x**2-5.0d0*x+1.0d0
    df=2.0d0*x-5.0d0
  end subroutine Root_Function_2_Both

  double precision function Root_Function_3(x)
    !!{
    Function for root finding unit tests.
    !!}
    double precision, intent(in   ) :: x

    Root_Function_3=x*exp(-x)+1.0d0
  end function Root_Function_3

  double precision function Root_Function_4(x)
    !!{
    Function for root finding unit tests.
    !!}
    double precision, intent(in   ) :: x

    if (abs(x) > 1.0d0) then
       Root_Function_4=sign(1.0d0,x)
    else
       Root_Function_4=x
    end if
    return
  end function Root_Function_4

  double precision function Root_Function_4_Derivative(x)
    !!{
    Function for root finding unit tests.
    !!}
    double precision, intent(in   ) :: x

    if (abs(x) > 1.0d0) then
       Root_Function_4_Derivative=0.0d0
    else
       Root_Function_4_Derivative=1.0d0
    end if
    return
  end function Root_Function_4_Derivative

  subroutine  Root_Function_4_Both(x,f,df)
    !!{
    Function for root finding unit tests.
    !!}
    double precision, intent(in   ) :: x
    double precision, intent(  out) :: f, df

    if (abs(x) > 1.0d0) then
       f=sign(1.0d0,x)
       df=0.0d0
    else
       f=x
       df=1.0d0
    end if
    return
  end subroutine Root_Function_4_Both

end module Test_Root_Finding_Functions
