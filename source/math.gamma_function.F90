!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
Contains a module which implements calculations of Gamma functions.
!!}

! Add dependency on GSL library.
!; gsl

module Gamma_Functions
  !!{
  Implements calculations of Gamma functions.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_double
  implicit none
  private
  public :: Gamma_Function_Incomplete        , Gamma_Function_Incomplete_Complementary        , Gamma_Function_Logarithmic            , Gamma_Function, &
       &    Inverse_Gamma_Function_Incomplete, Inverse_Gamma_Function_Incomplete_Complementary, Gamma_Function_Incomplete_Unnormalized

  interface
     function gsl_sf_gamma_inc(a,x) bind(c,name='gsl_sf_gamma_inc')
       !!{
       Template for the GSL unnormalized incomplete Gamma function.
       !!}
       import
       real(c_double)        :: gsl_sf_gamma_inc
       real(c_double), value :: a               , x
     end function gsl_sf_gamma_inc

     function gsl_sf_gamma_inc_Q(a,x) bind(c,name='gsl_sf_gamma_inc_Q')
       !!{
       Template for the GSL incomplete Gamma function.
       !!}
       import
       real(c_double)        :: gsl_sf_gamma_inc_Q
       real(c_double), value :: a                 , x
     end function gsl_sf_gamma_inc_Q

     function gsl_sf_gamma_inc_P(a,x) bind(c,name='gsl_sf_gamma_inc_P')
       !!{
       Template for the GSL complementary incomplete Gamma function.
       !!}
       import
       real(c_double)        :: gsl_sf_gamma_inc_P
       real(c_double), value :: a                 , x
     end function gsl_sf_gamma_inc_P

     function gsl_sf_gamma(x) bind(c,name='gsl_sf_gamma')
       !!{
       Template for the GSL Gamma function.
       !!}
       import
       real(c_double)        :: gsl_sf_gamma
       real(c_double), value :: x
     end function gsl_sf_gamma

     function gsl_sf_lngamma(x) bind(c,name='gsl_sf_lngamma')
       !!{
       Template for the GSL log-of-the-Gamma function.
       !!}
       import
       real(c_double)        :: gsl_sf_lngamma
       real(c_double), value :: x
     end function gsl_sf_lngamma
  end interface
  
contains

  double precision function Gamma_Function_Incomplete_Unnormalized(exponent,argument)
    !!{
    Computes the unnormalized incomplete Gamma function.
    !!}
    implicit none
    double precision, intent(in   ) :: argument, exponent

    Gamma_Function_Incomplete_Unnormalized=GSL_SF_Gamma_Inc(exponent,argument)
    return
  end function Gamma_Function_Incomplete_Unnormalized

  double precision function Gamma_Function_Incomplete(exponent,argument)
    !!{
    Computes the incomplete Gamma function.
    !!}
    implicit none
    double precision, intent(in   ) :: argument, exponent

    Gamma_Function_Incomplete=GSL_SF_Gamma_Inc_Q(exponent,argument)
    return
  end function Gamma_Function_Incomplete

  double precision function Gamma_Function_Incomplete_Complementary(exponent,argument)
    !!{
    Computes the complementary incomplete Gamma function.
    !!}
    implicit none
    double precision, intent(in   ) :: argument, exponent

    Gamma_Function_Incomplete_Complementary=GSL_SF_Gamma_Inc_P(exponent,argument)
    return
  end function Gamma_Function_Incomplete_Complementary

  double precision function Gamma_Function(exponent)
    !!{
    Computes the Gamma function.
    !!}
    implicit none
    double precision, intent(in   ) :: exponent

    Gamma_Function=GSL_SF_Gamma(exponent)
    return
  end function Gamma_Function

  double precision function Gamma_Function_Logarithmic(exponent)
    !!{
    Computes the logarithm of the Gamma function.
    !!}
    implicit none
    double precision, intent(in   ) :: exponent

    Gamma_Function_Logarithmic=GSL_SF_lnGamma(exponent)
    return
  end function Gamma_Function_Logarithmic

  double precision function Inverse_Gamma_Function_Incomplete_Complementary(a,P)
    !!{
    Returns the inverse of the incomplete function. That is, it returns $x$ given $P(a,x)$.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : varying_string         , assignment(=), operator(//)
    use :: Incomplete_Gamma  , only : GamInv
    implicit none
    double precision                , intent(in   ) :: P         , a
    integer                                         :: errorState
    double precision                                :: Q
    type            (varying_string)                :: message

    Q=1.0d0-P
    call GamInv(a,Inverse_Gamma_Function_Incomplete_Complementary,-1.0d0,P,Q,errorState)
    if (errorState < 0) then
       select case (errorState)
       case (-2)
          message='input error: a <= 0'
       case (-3)
          message='no solution obtained: A/a is too large'
       case (-4)
          message='input error: P or Q < 0 or P+Q != 1'
       case (-6)
          message='20 iterations were performed - this should not happen'
       case (-7)
          message='iteration failed'
       case (-8)
          message='accuracy lost'
       end select
       call Error_Report(message//{introspection:location})
    end if
    return
  end function Inverse_Gamma_Function_Incomplete_Complementary

  double precision function Inverse_Gamma_Function_Incomplete(a,Q)
    !!{
    Returns the inverse of the incomplete function. That is, it returns $x$ given $Q(a,x)$.
    !!}
    implicit none
    double precision, intent(in   ) :: Q, a
    double precision                :: P

    P=1.0d0-Q
    Inverse_Gamma_Function_Incomplete=Inverse_Gamma_Function_Incomplete_Complementary(a,P)
    return
  end function Inverse_Gamma_Function_Incomplete

end module Gamma_Functions
