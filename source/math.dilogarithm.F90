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
Contains a module which implements dilogarithms.
!!}

! Add dependency on GSL library.
!; gsl

module Dilogarithms
  !!{
  Implements dilogarithms.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_double, c_int
  use            :: Interface_GSL, only : gsl_sf_result
  implicit none
  private
  public :: Dilogarithm

  interface
     function gsl_sf_dilog(x) bind(c,name='gsl_sf_dilog')
       !!{
       Template for the GSL dilogarithm function.
       !!}
       import c_double
       real(c_double)        :: gsl_sf_dilog
       real(c_double), value :: x
     end function gsl_sf_dilog

     function gsl_sf_complex_dilog_e(x,y,a,b) bind(c,name='gsl_sf_complex_dilog_e')
       !!{
       Template for the GSL complex dilogarithm C function.
       !!}
       import
       integer(c_int        )        :: gsl_sf_complex_dilog_e
       real   (c_double     ), value :: x                     , y
       type   (gsl_sf_result)        :: a                     , b
     end function gsl_sf_complex_dilog_e
  end interface
  
  interface Dilogarithm
     module procedure Dilogarithm_Real
     module procedure Dilogarithm_Complex
  end interface Dilogarithm

contains

  double precision function Dilogarithm_Real(x)
    !!{
    Evaluate the $\hbox{Si}(x)\equiv\int_0^x \d t \sin(t)/t$ sine integral.
    !!}
    implicit none
    double precision, intent(in   ) :: x

    Dilogarithm_Real=GSL_SF_Dilog(x)
    return
  end function Dilogarithm_Real

  double complex function Dilogarithm_Complex(x)
    !!{
    Evaluate the dilogarithm for complex argument.
    !!}
    implicit none
    double complex               , intent(in   ) :: x
    type          (gsl_sf_result)                :: a     , b
    real          (c_double     )                :: r     , theta
    integer       (c_int        )                :: status

    r                  =sqrt (imag(x)**2+real(x)**2)
    theta              =atan2(imag(x)   ,real(x)   )
    status             =GSL_SF_Complex_Dilog_E(r,theta,a,b)
    Dilogarithm_Complex=dcmplx(a%val,b%val)
    return
  end function Dilogarithm_Complex

end module Dilogarithms
