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
Contains a module which implements calculations of factorials.
!!}

! Add dependency on GSL library.
!; gsl

module Factorials
  !!{
  Implements calculations of factorials
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_double, c_int
  implicit none
  private
  public :: Factorial, Logarithmic_Factorial, Logarithmic_Double_Factorial

  interface
     function gsl_sf_fact(n) bind(c,name='gsl_sf_fact')
       !!{
       Template for the GSL factorial function.
       !!}
       import
       real   (c_double)        :: gsl_sf_fact
       integer(c_int   ), value :: n
     end function gsl_sf_fact
  end interface
  
  interface
     function gsl_sf_lnfact(n) bind(c,name='gsl_sf_lnfact')
       !!{
       Template for the GSL logarithm of the factorial function.
       !!}
       import
       real   (c_double)        :: gsl_sf_lnfact
       integer(c_int   ), value :: n
     end function gsl_sf_lnfact
  end interface
  
contains

  double precision function Factorial(argument)
    !!{
    Computes the factorial of {\normalfont \ttfamily argument}.
    !!}
    implicit none
    integer, intent(in   ) :: argument

    Factorial=GSL_SF_Fact(argument)
    return
  end function Factorial

  double precision function Logarithmic_Factorial(argument)
    !!{
    Computes the logarithmic of the factorial of {\normalfont \ttfamily argument}.
    !!}
    implicit none
    integer, intent(in   ) :: argument

    Logarithmic_Factorial=GSL_SF_LnFact(argument)
    return
  end function Logarithmic_Factorial

  double precision function Logarithmic_Double_Factorial(argument)
    !!{
    Computes the natural logarithm of the double factorial, $k!!$.
    !!}
    implicit none
    integer, intent(in   ) :: argument
    integer                :: i

    Logarithmic_Double_Factorial=0.0d0
    i=argument
    do while (i >= 2)
       Logarithmic_Double_Factorial=Logarithmic_Double_Factorial+log(dble(i))
       i=i-2
    end do
    return
  end function Logarithmic_Double_Factorial

end module Factorials
