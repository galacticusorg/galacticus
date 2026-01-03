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
Contains a module which implements trigonometric functions.
!!}

module Trigonometric_Functions
  !!{
  Implements trigonometric functions.
  !!}
  implicit none
  private
  public :: cot, cosec, hypotenuse

  interface cot
     module procedure cotDouble
     module procedure cotDoubleComplex
  end interface cot

  interface cosec
     module procedure cosecDouble
     module procedure cosecDoubleComplex
  end interface cosec

contains

  double precision function hypotenuse(x)
    !!{
    Compute the $N$-dimensional hypotenuse, $(\sum_{i=1}^N x_i^2)^{1/2}$ avoiding undue floating point overflow.
    !!}
    implicit none
    double precision, intent(in   ), dimension(:) :: x
    double precision                              :: xMaximum

    xMaximum  =maxval(abs(x))
    hypotenuse=xMaximum*sqrt(sum((x/xMaximum)**2))
    return
  end function hypotenuse
  
  double precision function cotDouble(x)
    !!{
    Implements cotangent for double precision {\normalfont \ttfamily x}.
    !!}
    implicit none
    double precision, intent(in   ) :: x

    cotDouble=1.0d0/tan(x)
    return
  end function cotDouble

  double complex function cotDoubleComplex(x)
    !!{
    Implements cotangent for double precision complex {\normalfont \ttfamily x}.
    !!}
    implicit none
    double complex, intent(in   ) :: x

    cotDoubleComplex=1.0d0/tan(x)
    return
  end function cotDoubleComplex

  double precision function cosecDouble(x)
    !!{
    Implements cosecant for double precision {\normalfont \ttfamily x}.
    !!}
    implicit none
    double precision, intent(in   ) :: x

    cosecDouble=1.0d0/sin(x)
    return
  end function cosecDouble

  double complex function cosecDoubleComplex(x)
    !!{
    Implements cosecant for double precision complex {\normalfont \ttfamily x}.
    !!}
    implicit none
    double complex, intent(in   ) :: x

    cosecDoubleComplex=1.0d0/sin(x)
    return
  end function cosecDoubleComplex

end module Trigonometric_Functions
