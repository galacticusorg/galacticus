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
Contains a module which implements exponential integrals.
!!}

! Add dependency on GSL library.
!; gsl

module Exponential_Integrals
  !!{
  Implements exponential integrals.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_double
  implicit none
  private
  public :: Sine_Integral, Cosine_Integral, Exponential_Integral

  interface Exponential_Integral
     module procedure Exponential_Integral_Double
     module procedure Exponential_Integral_Double_Complex
  end interface Exponential_Integral

  interface
     function gsl_sf_Si(x) bind(c,name='gsl_sf_Si')
       !!{
       Template for the GSL Sine integral function.
       !!}
       import c_double
       real(c_double)        :: gsl_sf_Si
       real(c_double), value :: x
     end function gsl_sf_Si
     function gsl_sf_Ci(x) bind(c,name='gsl_sf_Ci')
       !!{
       Template for the GSL Cosine integral function.
       !!}
       import c_double
       real(c_double)        :: gsl_sf_Ci
       real(c_double), value :: x
     end function gsl_sf_Ci
  end interface
  
contains

  double precision function Sine_Integral(x)
    !!{
    Evaluate the $\hbox{Si}(x)\equiv\int_0^x \d t \sin(t)/t$ sine integral.
    !!}
    implicit none
    double precision, intent(in   ) :: x

    Sine_Integral=GSL_SF_Si(x)
    return
  end function Sine_Integral

  double precision function Cosine_Integral(x)
    !!{
    Evaluate the $\hbox{Ci}(x)\equiv\int_0^x \d t \cos(t)/t$ cosine integral.
    !!}
    implicit none
    double precision, intent(in   ) :: x

    Cosine_Integral=GSL_SF_Ci(x)
    return
  end function Cosine_Integral

  double complex function Exponential_Integral_Double_Complex(z)
    !!{
    Exponential integral, $E_\mathrm{i}(z)$, for complex argument {\normalfont \ttfamily z}.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double complex, intent(in   ) :: z

    call e1z(-z,Exponential_Integral_Double_Complex)
    Exponential_Integral_Double_Complex=-Exponential_Integral_Double_Complex+dcmplx(0.0d0,1.0d0)*Pi
    return
  end function Exponential_Integral_Double_Complex

  double precision function Exponential_Integral_Double(x)
    !!{
    Exponential integral for real argument {\normalfont \ttfamily x}.
    !!}
    implicit none
    double precision, intent(in   ) :: x

    Exponential_Integral_Double=dreal(Exponential_Integral_Double_Complex(dcmplx(x,0.0d0)))
    return
  end function Exponential_Integral_Double

  ! Include functions.
  include "math.exponential_integrals.functions.inc"

end module Exponential_Integrals
