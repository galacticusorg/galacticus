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
Contains a module which implements calculations of Bessel functions.
!!}

! Add dependency on GSL library.
!; gsl

module Bessel_Functions
  !!{
  Implements calculations of Bessel functions.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_double, c_int
  implicit none
  private
  public :: Bessel_Function_J0, Bessel_Function_J1     , Bessel_Function_K0     , Bessel_Function_K1     , &
       &    Bessel_Function_I0, Bessel_Function_I1     , Bessel_Function_J0_Zero, Bessel_Function_J1_Zero, &
       &    Bessel_Function_Jn, Bessel_Function_Jn_Zero, Bessel_Function_In

  interface
     function gsl_sf_bessel_J0(x) bind(c,name='gsl_sf_bessel_J0')
       !!{
       Template for the GSL $J_0$ Bessel function.
       !!}
       import
       real(c_double)        :: gsl_sf_bessel_J0
       real(c_double), value :: x
     end function gsl_sf_bessel_J0

     function gsl_sf_bessel_J1(x) bind(c,name='gsl_sf_bessel_J1')
       !!{
       Template for the GSL $J_1$ Bessel function.
       !!}
       import
       real(c_double)        :: gsl_sf_bessel_J1
       real(c_double), value :: x
     end function gsl_sf_bessel_J1
     
     function gsl_sf_bessel_Jn(n,x) bind(c,name='gsl_sf_bessel_Jn')
       !!{
       Template for the GSL $J_\mathrm{n}$ Bessel function.
       !!}
       import
       real   (c_double)        :: gsl_sf_bessel_Jn
       integer(c_int   ), value :: n
       real   (c_double), value :: x
     end function gsl_sf_bessel_Jn
     
     function gsl_sf_bessel_zero_J0(s) bind(c,name='gsl_sf_bessel_zero_J0')
       !!{
       Template for the GSL zeros-of-the-$J_0$ Bessel function.
       !!}
       import
       real   (c_double)        :: gsl_sf_bessel_zero_J0
       integer(c_int   ), value :: s
     end function gsl_sf_bessel_zero_J0

     function gsl_sf_bessel_zero_J1(s) bind(c,name='gsl_sf_bessel_zero_J1')
       !!{
       Template for the GSL zeros-of-the-$J_1$ Bessel function.
       !!}
       import
       real   (c_double)        :: gsl_sf_bessel_zero_J1
       integer(c_int   ), value :: s
     end function gsl_sf_bessel_zero_J1

     function gsl_sf_bessel_zero_Jnu(nu,s) bind(c,name='gsl_sf_bessel_zero_Jnu')
       !!{
       Template for the GSL zeros-of-the-$J_\mathrm{n}$ Bessel function.
       !!}
       import
       real   (c_double)        :: gsl_sf_bessel_zero_Jnu
       real   (c_double), value :: nu
       integer(c_int   ), value :: s
     end function gsl_sf_bessel_zero_Jnu

     function gsl_sf_bessel_I0(x) bind(c,name='gsl_sf_bessel_I0')
       !!{
       Template for the GSL $I_0$ Bessel function.
       !!}
       import
       real(c_double)        :: gsl_sf_bessel_I0
       real(c_double), value :: x
     end function gsl_sf_bessel_I0

     function gsl_sf_bessel_I1(x) bind(c,name='gsl_sf_bessel_I1')
       !!{
       Template for the GSL $I_1$ Bessel function.
       !!}
       import
       real(c_double)        :: gsl_sf_bessel_I1
       real(c_double), value :: x
     end function gsl_sf_bessel_I1

     function gsl_sf_bessel_In(n,x) bind(c,name='gsl_sf_bessel_In')
       !!{
       Template for the GSL $I_\mathrm{n}$ Bessel function.
       !!}
       import
       real   (c_double)        :: gsl_sf_bessel_In
       integer(c_int   ), value :: n
       real   (c_double), value :: x
     end function gsl_sf_bessel_In

     function gsl_sf_bessel_Inu(nu,x) bind(c,name='gsl_sf_bessel_Inu')
       !!{
       Template for the GSL $I_\nu$ Bessel function.
       !!}
       import
       real(c_double)        :: gsl_sf_bessel_Inu
       real(c_double), value :: nu              , x
     end function gsl_sf_bessel_Inu

     function gsl_sf_bessel_K0(x) bind(c,name='gsl_sf_bessel_K0')
       !!{
       Template for the GSL $K_0$ Bessel function.
       !!}
       import
       real(c_double)        :: gsl_sf_bessel_K0
       real(c_double), value :: x
     end function gsl_sf_bessel_K0

     function gsl_sf_bessel_K1(x) bind(c,name='gsl_sf_bessel_K1')
       !!{
       Template for the GSL $K_1$ Bessel function.
       !!}
       import
       real(c_double)        :: gsl_sf_bessel_K1
       real(c_double), value :: x
     end function gsl_sf_bessel_K1
  end interface

  ! Generic interfaces.
  interface Bessel_Function_In
     module procedure Bessel_Function_In_Integer_Order
     module procedure Bessel_Function_In_Fractional_Order
  end interface Bessel_Function_In
     
contains

  double precision function Bessel_Function_J0(argument)
    !!{
    Computes the $J_0$ Bessel function.
    !!}
    implicit none
    double precision, intent(in   ) :: argument

    Bessel_Function_J0=GSL_SF_Bessel_J0(argument)
    return
  end function Bessel_Function_J0

  double precision function Bessel_Function_J1(argument)
    !!{
    Computes the $J_1$ Bessel function.
    !!}
    implicit none
    double precision, intent(in   ) :: argument

    Bessel_Function_J1=GSL_SF_Bessel_J1(argument)
    return
  end function Bessel_Function_J1

  double precision function Bessel_Function_Jn(n,argument)
    !!{
    Computes the $J_n$ Bessel function.
    !!}
    implicit none
    integer         , intent(in   ) :: n
    double precision, intent(in   ) :: argument

    Bessel_Function_Jn=GSL_SF_Bessel_Jn(n,argument)
    return
  end function Bessel_Function_Jn

  double precision function Bessel_Function_J0_Zero(s)
    !!{
    Computes the $s^\mathrm{th}$ zero of the $J_0$ Bessel function.
    !!}
    implicit none
    integer, intent(in   ) :: s

    Bessel_Function_J0_Zero=GSL_SF_Bessel_Zero_J0(s)
    return
  end function Bessel_Function_J0_Zero

  double precision function Bessel_Function_J1_Zero(s)
    !!{
    Computes the $s^\mathrm{th}$ zero of the $J_1$ Bessel function.
    !!}
    implicit none
    integer, intent(in   ) :: s

    Bessel_Function_J1_Zero=GSL_SF_Bessel_Zero_J1(s)
    return
  end function Bessel_Function_J1_Zero

  double precision function Bessel_Function_Jn_Zero(n,s)
    !!{
    Computes the $s^\mathrm{th}$ zero of the $J_1$ Bessel function.
    !!}
    implicit none
    double precision, intent(in   ) :: n
    integer         , intent(in   ) :: s

    Bessel_Function_Jn_Zero=GSL_SF_Bessel_Zero_Jnu(n,s)
    return
  end function Bessel_Function_Jn_Zero

  double precision function Bessel_Function_K0(argument)
    !!{
    Computes the $K_0$ Bessel function.
    !!}
    implicit none
    double precision, intent(in   ) :: argument

    Bessel_Function_K0=GSL_SF_Bessel_K0(argument)
    return
  end function Bessel_Function_K0

  double precision function Bessel_Function_K1(argument)
    !!{
    Computes the $K_1$ Bessel function.
    !!}
    implicit none
    double precision, intent(in   ) :: argument

    Bessel_Function_K1=GSL_SF_Bessel_K1(argument)
    return
  end function Bessel_Function_K1

  double precision function Bessel_Function_I0(argument)
    !!{
    Computes the $I_0$ Bessel function.
    !!}
    implicit none
    double precision, intent(in   ) :: argument

    Bessel_Function_I0=GSL_SF_Bessel_I0(argument)
    return
  end function Bessel_Function_I0

  double precision function Bessel_Function_I1(argument)
    !!{
    Computes the $I_1$ Bessel function.
    !!}
    implicit none
    double precision, intent(in   ) :: argument

    Bessel_Function_I1=GSL_SF_Bessel_I1(argument)
    return
  end function Bessel_Function_I1

  double precision function Bessel_Function_In_Integer_Order(n,argument) result(In)
    !!{
    Computes the $I_n$ Bessel function for integer order.
    !!}
    implicit none
    integer         , intent(in   ) :: n
    double precision, intent(in   ) :: argument

    In=GSL_SF_Bessel_In(n,argument)
    return
  end function Bessel_Function_In_Integer_Order

  double precision function Bessel_Function_In_Fractional_Order(n,argument) result(In)
    !!{
    Computes the $I_n$ Bessel function for fractional order.
    !!}
    implicit none
    double precision, intent(in   ) :: n, argument

    In=GSL_SF_Bessel_Inu(n,argument)
    return
  end function Bessel_Function_In_Fractional_Order

end module Bessel_Functions
